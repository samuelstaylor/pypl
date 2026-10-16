#!/usr/bin/env python3
"""Regression and sanity checks for the 1D CCD path.

Run this after touching pypl/config_coord_1d_solver.py (or the parsers it depends
on) to confirm nothing moved:

    python testing/check_1d_ccd.py            # quiet unless something fails
    python testing/check_1d_ccd.py -v         # print every check

Exit code is 0 if all checks pass, 1 otherwise. Reference values come from
doc/tutorials/pypl_002.ipynb; update REFERENCE below only when a change to the
physics is intended.
"""

import argparse
import io
import contextlib
import sys

import numpy as np

from harness import ccd_data, fits
from pypl.config_coord_1d_solver import config_coord_1d_solver
from pypl_ext.config_coord_1d_vof_solver import config_coord_1d_vof_solver

# Values produced by tutorial 002 for the NV center in diamond.
REFERENCE = {
    "delta_q": 0.695732,
    "gs_phonon": 59.08890,
    "es_phonon": 64.31191,
    "gs_hrf": 3.4211,
    "es_hrf": 3.7235,
}

ENE_AXIS = np.linspace(-200.0, 1200.0, 1401)
SIGMA = [8.0, 25.0]
GAMMA = 0.3
TEMP = 5.0
EZPL = 1945.0


class Checker:
    """Minimal pass/fail collector so a failure does not abort the remaining checks."""

    def __init__(self, verbose=False):
        self.verbose = verbose
        self.failures = []
        self.n = 0

    def check(self, name, ok, detail=""):
        self.n += 1
        if ok:
            if self.verbose:
                print("  PASS  %-46s %s" % (name, detail))
        else:
            print("  FAIL  %-46s %s" % (name, detail))
            self.failures.append(name)
        return ok

    def close(self, name, got, expected, tol, unit=""):
        detail = "got % .6g  expected % .6g  (tol %g%s)" % (got, expected, tol, unit)
        return self.check(name, abs(got - expected) <= tol, detail)

    def report(self):
        print("\n%d/%d checks passed." % (self.n - len(self.failures), self.n))
        for name in self.failures:
            print("  failed: %s" % name)
        return 1 if self.failures else 0


def quiet(fn, *a, **kw):
    """Call fn while swallowing the library's progress prints."""

    with contextlib.redirect_stdout(io.StringIO()):
        return fn(*a, **kw)


def make_solver(omega_i, omega_f, delta_q, ni, nf, temp=TEMP, eneaxis=ENE_AXIS):
    solver = config_coord_1d_solver(omega_i, omega_f, delta_q)
    solver.compute_franck_condon_integrals(ni=ni, nf=nf)
    quiet(solver.bulid_fc_lsp, eneaxis=eneaxis, temp=temp, sigma=SIGMA, zpl_lorentzian=True, gamma=GAMMA)
    return solver


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("-v", "--verbose", action="store_true", help="print passing checks too")
    args = p.parse_args()
    c = Checker(args.verbose)
    # The reference numbers are tutorial 002's, whatever PYPL_TESTING_DATASET says.
    ccd_data.use_dataset("nv_diamond")

    # -------------------------------------------------- CCD construction
    print("CCD construction (parsers + parabolic fits)")
    ccd = quiet(ccd_data.build_ccd)
    c.close("delta_q", ccd["delta_q"], REFERENCE["delta_q"], 1e-5, " amu^0.5 A")
    c.close("gs_phonon", ccd["gs_phonon"], REFERENCE["gs_phonon"], 1e-3, " meV")
    c.close("es_phonon", ccd["es_phonon"], REFERENCE["es_phonon"], 1e-3, " meV")
    c.close("gs_hrf", ccd["gs_hrf"], REFERENCE["gs_hrf"], 1e-3)
    c.close("es_hrf", ccd["es_hrf"], REFERENCE["es_hrf"], 1e-3)
    c.check("ES frequency exceeds GS frequency", ccd["es_phonon"] > ccd["gs_phonon"])

    delta_q, w_gs, w_es = ccd["delta_q"], ccd["gs_phonon"], ccd["es_phonon"]

    # -------------------------------------------------- Franck-Condon integrals
    print("\nFranck-Condon integrals")
    solver = config_coord_1d_solver(w_es, w_gs, delta_q)
    solver.compute_franck_condon_integrals(ni=10, nf=300)
    # Completeness: each initial state must have unit total overlap with the final basis.
    row_sums = np.sum(solver.fc_ints**2, axis=1)
    c.close("closure sum_j |F_0j|^2", row_sums[0], 1.0, 1e-8)
    c.check(
        "closure holds for i < 10", np.allclose(row_sums, 1.0, atol=1e-6), "max dev %.2e" % np.max(np.abs(row_sums - 1))
    )
    c.check("all FC integrals finite", np.all(np.isfinite(solver.fc_ints)))
    # |F_00|^2 is the Debye-Waller factor; for a displaced oscillator it is well below 1.
    c.check("0 < |F_00|^2 < 1", 0.0 < solver.fc_ints[0, 0] ** 2 < 1.0, "|F_00|^2 = %.4e" % solver.fc_ints[0, 0] ** 2)

    # Zero displacement must give the identity matrix (orthonormal, undisplaced oscillators).
    undisplaced = config_coord_1d_solver(w_gs, w_gs, 0.0)
    undisplaced.compute_franck_condon_integrals(ni=8, nf=8)
    c.check(
        "delta_q = 0 gives identity overlap",
        np.allclose(undisplaced.fc_ints, np.eye(8), atol=1e-10),
        "max dev %.2e" % np.max(np.abs(undisplaced.fc_ints - np.eye(8))),
    )

    # -------------------------------------------------- lineshape
    print("\nFC lineshape")
    pl = make_solver(w_es, w_gs, delta_q, 50, 60)
    integral = np.sum(pl.fc_lineshape) * (ENE_AXIS[1] - ENE_AXIS[0])
    c.close("lineshape integrates to 1", integral, 1.0, 0.05)
    c.check("lineshape is non-negative", np.all(pl.fc_lineshape >= 0.0), "min %.3e" % pl.fc_lineshape.min())
    c.close("Boltzmann weights sum to 1", float(np.sum(pl.prefactors)), 1.0, 1e-12)
    c.check(
        "T = 5 K populates only the lowest level",
        pl.prefactors[0] > 1.0 - 1e-6,
        "p_0 = %.8f" % pl.prefactors[0],
    )
    hot = make_solver(w_es, w_gs, delta_q, 50, 60, temp=800.0)
    c.check("T = 800 K populates excited levels", hot.prefactors[1] > 1e-3, "p_1 = %.4e" % hot.prefactors[1])

    # Convergence: raising the vibrational basis must not move the spectrum.
    converged = make_solver(w_es, w_gs, delta_q, 80, 100)
    rel = np.max(np.abs(converged.fc_lineshape - pl.fc_lineshape)) / np.max(pl.fc_lineshape)
    c.check("order (50,60) is converged vs (80,100)", rel < 1e-4, "max rel dev %.2e" % rel)

    # Equal frequencies collapse the recurrences to a = c = 0, b = -d, so the FC matrix
    # must satisfy F_ij = (-1)^(i+j) F_ji. This catches index slips in the recursion.
    eq = config_coord_1d_solver(w_gs, w_gs, delta_q)
    eq.compute_franck_condon_integrals(ni=30, nf=30)
    signs = np.fromfunction(lambda i, j: (-1.0) ** (i + j), eq.fc_ints.shape)
    dev = np.max(np.abs(eq.fc_ints - signs * eq.fc_ints.T))
    c.check("equal frequencies give F_ij = (-1)^(i+j) F_ji", dev < 1e-9, "max dev %.2e" % dev)

    # With equal frequencies the two transition directions are the same problem, so the
    # emission and absorption lineshapes must agree point by point.
    eq_pl = make_solver(w_gs, w_gs, delta_q, 60, 60)
    eq_ab = make_solver(w_gs, w_gs, delta_q, 60, 60)
    c.check(
        "equal frequencies give identical PL and Abs lineshapes",
        np.allclose(eq_pl.fc_lineshape, eq_ab.fc_lineshape, rtol=0, atol=1e-14),
    )

    # -------------------------------------------------- spectra
    print("\nPL and absorption spectra")
    pl_axis, pl_spec = pl.compute_spectrum(tdm=1.0, zpl=EZPL, spectrum_type="PL")
    ab = make_solver(w_gs, w_es, delta_q, 60, 50)
    ab_axis, ab_spec = ab.compute_spectrum(tdm=1.0, zpl=EZPL, spectrum_type="Abs")

    c.check("PL axis runs downward from the ZPL", pl_axis[0] > pl_axis[-1] and pl_axis.max() > EZPL)
    c.check("Abs axis runs upward from the ZPL", ab_axis[0] < ab_axis[-1] and ab_axis.min() < EZPL)
    c.check("PL is non-negative", np.all(pl_spec >= 0.0))
    c.check("Abs is non-negative", np.all(ab_spec >= 0.0))
    c.close("PL peaks at the ZPL", pl_axis[np.argmax(pl_spec)], EZPL, 1.0, " meV")
    c.close("Abs peaks at the ZPL", ab_axis[np.argmax(ab_spec)], EZPL, 1.0, " meV")

    # The phonon sideband sits below the ZPL in emission and above it in absorption.
    pl_norm = ccd_data.normalize(pl_axis, pl_spec)
    ab_norm = ccd_data.normalize(ab_axis, ab_spec)
    c.check("PL weight lies below the ZPL", np.sum(pl_norm[pl_axis < EZPL]) > np.sum(pl_norm[pl_axis > EZPL]))
    c.check("Abs weight lies above the ZPL", np.sum(ab_norm[ab_axis > EZPL]) > np.sum(ab_norm[ab_axis < EZPL]))

    try:
        pl.compute_spectrum(spectrum_type="bogus")
        c.check("invalid spectrum_type raises ValueError", False, "no exception raised")
    except ValueError:
        c.check("invalid spectrum_type raises ValueError", True)

    check_extension(c, ccd)
    check_datasets(c)
    return c.report()


def check_datasets(c):
    """Large-n harmonic functions (needed for large Huang-Rhys factors) and the NV- in MgO loader."""

    import os

    print()
    print("Large Huang-Rhys factors and the NV- in MgO dataset")
    sub = config_coord_1d_vof_solver
    x = np.linspace(-6, 6, 2001)
    chi = sub.harmonic_wavefunctions(x, 20, 59.0)
    xi = np.sqrt(sub.harmonic_alpha(59.0)) * x
    explicit = np.array(
        [
            (sub.harmonic_alpha(59.0) / np.pi) ** 0.25
            / np.sqrt(2.0**n * float(np.prod(np.arange(1, n + 1))))
            * np.polynomial.hermite.hermval(xi, np.eye(21)[n])
            * np.exp(-(xi**2) / 2)
            for n in range(21)
        ]
    )
    dev = np.max(np.abs(chi - explicit)) / np.max(np.abs(explicit))
    c.check("recurrence = explicit H_n formula (n <= 20)", dev < 1e-12, "max rel dev %.1e" % dev)
    q = np.linspace(-30, 30, 200001)
    chi = sub.harmonic_wavefunctions(q, 300, 43.8)
    gram = (chi * (q[1] - q[0])) @ chi.T
    dev = np.max(np.abs(gram - np.eye(301)))
    c.check("harmonic functions orthonormal up to n = 300", np.isfinite(dev) and dev < 1e-9, "max dev %.1e" % dev)
    s = sub(74.07, 43.77, 3.0979)
    s.compute_franck_condon_integrals_hermite(2, 200, npts=100001)
    c.check("hermite F_ij finite with 200 levels (no 2^n n! overflow)", np.all(np.isfinite(s.fc_ints)))
    c.close("hermite closure of row 0, S ~ 50", np.sum(s.fc_ints[0] ** 2), 1.0, 1e-10)

    if not os.path.isdir(ccd_data.MGO_DATA_DIR):
        print("  (doc/nv_mgo_west not present: MgO loader checks skipped)")
        return
    mgo = quiet(ccd_data.build_ccd, dataset="nv_mgo")
    c.close("MgO delta_q (ccd_nv_mgo.py: 3.097879)", mgo["delta_q"], 3.097879, 1e-6, " amu^0.5 A")
    c.check(
        "MgO images: every one, 30 x 1A1 (no sf +2.00), 33 x 3E",
        (mgo["n_gs_images"], mgo["n_es_images"]) == (30, 33),
    )
    q_gs, e_gs = ccd_data.samples(mgo, "gs")
    q_es, e_es = ccd_data.samples(mgo, "es")
    c.check(
        "MgO 3E minimum at Q = 0, 1A1 at Q = delta_q",
        (ccd_data.center(mgo, "es"), ccd_data.center(mgo, "gs")) == (0.0, mgo["delta_q"]),
    )
    c.check("MgO 3E lowest image at Q = 0 (q = +0.00)", q_es[np.argmin(e_es)] == 0.0)
    c.close("MgO 1A1 lowest image at Q = delta_q (q = +1.00)", q_gs[np.argmin(e_gs)], mgo["delta_q"], 1e-12)
    c.close("MgO 3E minimum vs 3A2 zero (ccd_nv_mgo.py: 3.1376 eV)", e_es.min(), 3.1376, 5e-5, " eV")
    c.close("MgO ZPL from the images (ccd_nv_mgo.py: 554.4 meV)", 1e3 * (e_es.min() - e_gs.min()), 554.4, 0.05, " meV")
    c.check("dataset choice leaves STATE_NAMES on the active one", ccd_data.STATE_NAMES == {"gs": "GS", "es": "ES"})
    q_es_path, q_gs_path = mgo["es_path_q"], mgo["gs_path_q"]
    par_es, par_gs = q_es_path[mgo["es_potential"].index], q_gs_path[mgo["gs_potential"].index]
    c.check(
        "MgO parabola images: 3E q = -0.10..+0.45, 1A1 q = +0.60..+1.15",
        (par_es.min(), par_es.max(), len(par_es), par_gs.min(), par_gs.max(), len(par_gs))
        == (-0.10, 0.45, 12, 0.60, 1.15, 12),
    )
    cub = {s: mgo["%s_path_q" % s][fits.make_potential(mgo, s, "cubic").index] for s in ccd_data.STATES}
    c.check(
        "MgO cubic-fit images: q = -0.10..+1.15 for both (26 each)",
        all((v.min(), v.max(), len(v)) == (-0.10, 1.15, 26) for v in cub.values()),
    )
    c.check(
        "MgO spline / Morse / polyN use every image",
        all(
            len(fits.make_potential(mgo, s, spec).index) == mgo["n_%s_images" % s]
            for s in ccd_data.STATES
            for spec in ("spline", "morse", "poly4")
        ),
    )

    # Smooth (C2) spline tails: the spline's own value, slope and curvature continued past the end images.
    smooth = fits.make_potential(mgo, "gs", "spline3-notaknot-smooth")
    plain = fits.make_potential(mgo, "gs", "spline3-notaknot")
    q_in = np.linspace(q_gs[0], q_gs[-1], 4001)
    c.check("smooth-tail spline = not-a-knot spline inside the data", np.array_equal(smooth(q_in), plain(q_in)))
    f = lambda x: float(smooth(np.array([x]))[0])  # noqa: E731
    h = 1e-4
    jumps = [
        max(
            abs((f(e + h) - f(e)) / h - (f(e) - f(e - h)) / h),
            abs((f(e + 2 * h) - 2 * f(e + h) + f(e)) / h**2 - (f(e) - 2 * f(e - h) + f(e - 2 * h)) / h**2) * h,
        )
        for e in (q_gs[0], q_gs[-1])
    ]
    c.check(
        "smooth tails are C2 at both end images",
        smooth.params["junction"] == "C2" and max(jumps) < 1e-3,
        "max jump %.1e" % max(jumps),
    )
    c.check(
        "smooth tails bounded on MgO (V'' > 0 at every end)",
        all(fits.make_potential(mgo, st, "spline3-notaknot-smooth").bounded for st in ccd_data.STATES),
    )
    try:
        fits.parse_spec("spline3-clamped-smooth")
        c.check("clamped + smooth tails is rejected", False, "no exception")
    except ValueError:
        c.check("clamped + smooth tails is rejected", True)
    c.check(
        "spline_spec(3, not-a-knot, tails=smooth) is canonical",
        fits.spline_spec(3, "not-a-knot", None, "smooth")
        == fits.canonical("spline3-notaknot-smooth")
        == "spline3-notaknot-smooth",
    )

    # build_fc_lsp_from_levels with harmonic ladders is bulid_fc_lsp, even with thermal populations.
    lv = sub(64.31, 59.09, 0.6957)
    lv.compute_franck_condon_integrals(ni=20, nf=25)
    axis = np.linspace(-200, 1200, 1401)
    kw = dict(eneaxis=axis, temp=300, sigma=[8, 25], zpl_lorentzian=True, gamma=0.3)
    ref_ls = quiet(lv.bulid_fc_lsp, **kw)[1].copy()
    lev_ls = quiet(
        lv.build_fc_lsp_from_levels, energies_i=64.31 * np.arange(20), energies_f=59.09 * np.arange(25), **kw
    )[1]
    c.check(
        "build_fc_lsp_from_levels(harmonic ladders) = bulid_fc_lsp",
        np.max(np.abs(lev_ls - ref_ls)) <= 1e-14 * ref_ls.max(),
        "max dev %.1e" % np.max(np.abs(lev_ls - ref_ls)),
    )

    # With the initial (3E) well on the left of the final one, fd must still give the solver-frame F_ij.
    w_es, w_gs, dq = mgo["es_phonon"], mgo["gs_phonon"], mgo["delta_q"]
    ref = sub(w_es, w_gs, dq)
    ref.compute_franck_condon_integrals_hermite(12, 12)
    fd = sub(w_es, w_gs, dq)
    fd.compute_franck_condon_integrals_fd(
        12, 12, mgo["es_potential"], mgo["gs_potential"], q_i=ccd_data.center(mgo, "es"), q_f=ccd_data.center(mgo, "gs")
    )
    dev = np.max(np.abs(fd.fc_ints - ref.fc_ints))
    c.check("fd = hermite with the 3E well left of 1A1 (MgO orientation)", dev < 1e-7, "max dF %.1e" % dev)


def check_extension(c, ccd):
    """pypl_ext must only add to PyPL: nothing overridden, inherited results bit-identical."""

    print()
    print("Extension (pypl_ext) leaves PyPL unchanged")
    base, sub = config_coord_1d_solver, config_coord_1d_vof_solver
    overridden = [k for k in vars(sub) if k in vars(base) and not k.startswith("__")]
    c.check("vof solver overrides no config_coord_1d_solver attribute", not overridden, "overridden: %s" % overridden)
    c.check("vof solver is a config_coord_1d_solver", issubclass(sub, base))

    w_gs, w_es, delta_q = ccd["gs_phonon"], ccd["es_phonon"], ccd["delta_q"]
    old, new = base(w_es, w_gs, delta_q), sub(w_es, w_gs, delta_q)
    old.compute_franck_condon_integrals(ni=50, nf=60)
    new.compute_franck_condon_integrals(ni=50, nf=60)
    c.check("inherited compute_franck_condon_integrals is bit-identical", np.array_equal(old.fc_ints, new.fc_ints))
    for solver in (old, new):
        quiet(solver.bulid_fc_lsp, eneaxis=ENE_AXIS, temp=TEMP, sigma=SIGMA, zpl_lorentzian=True, gamma=GAMMA)
    c.check("inherited bulid_fc_lsp is bit-identical", np.array_equal(old.fc_lineshape, new.fc_lineshape))
    same = all(np.array_equal(a, b) for a, b in zip(old.compute_spectrum(zpl=EZPL), new.compute_spectrum(zpl=EZPL)))
    c.check("inherited compute_spectrum is bit-identical", same)

    print()
    print("New F_ij routes (20 x 20 block)")
    ref = sub(w_es, w_gs, delta_q)
    ref.compute_franck_condon_integrals(ni=20, nf=20)
    herm = sub(w_es, w_gs, delta_q)
    herm.compute_franck_condon_integrals_hermite(ni=20, nf=20)
    dev = np.max(np.abs(herm.fc_ints - ref.fc_ints))
    c.check("hermite F_ij matches the Ruhoff recurrences", dev < 1e-10, "max dev %.2e" % dev)
    fd = sub(w_es, w_gs, delta_q)
    fd.compute_franck_condon_integrals_fd(20, 20, ccd["es_potential"], ccd["gs_potential"], q_i=delta_q, q_f=0.0)
    dev = np.max(np.abs(fd.fc_ints - ref.fc_ints))
    c.check("fd F_ij on the parabolas matches the Ruhoff recurrences", dev < 1e-8, "max dev %.2e" % dev)
    c.check(
        "fd eigenstates are bound (|psi| at the walls < 1e-12)", fd.fd_info["wall"] < 1e-12, "%.1e" % fd.fd_info["wall"]
    )

    print()
    print("Potential fits")
    spline = fits.make_potential(ccd, "gs", "spline")
    c.check(
        "spline passes through every image",
        np.max(np.abs(spline.residuals())) < 1e-12,
        "max residual %.1e eV" % np.max(np.abs(spline.residuals())),
    )
    c.check("degree-4 polynomial is flagged unbounded", not fits.make_potential(ccd, "gs", "poly4").bounded)
    c.close("parabola phonon_energy is the solver's omega", ccd["gs_potential"].phonon_energy(), w_gs, 1e-12, " meV")
    check_spline_options(c, ccd)


def check_spline_options(c, ccd):
    """Every spline order and end condition behaves as documented."""

    from scipy.interpolate import CubicSpline

    print()
    print("Spline options")
    q, e = ccd_data.samples(ccd, "gs")
    x = np.linspace(q.min(), q.max(), 2001)
    dev = np.max(np.abs(fits.make_potential(ccd, "gs", "spline")(x) - CubicSpline(q, e, bc_type="natural")(x)))
    c.check("default spline is the natural cubic spline", dev < 1e-12, "max dev %.1e eV" % dev)
    for spec in ("spline3-notaknot", "spline3-clamped", "spline4-notaknot", "spline5-natural", "spline5-clamped"):
        res = np.max(np.abs(fits.make_potential(ccd, "gs", spec).residuals()))
        c.check("%s passes through every image" % spec, res < 1e-12, "max residual %.1e eV" % res)
    c.check(
        "cubic clamped spline joins its tails in C1",
        fits.make_potential(ccd, "gs", "spline3-clamped").params["junction"] == "C1",
    )
    c.check(
        "quintic clamped spline joins its tails in C2",
        fits.make_potential(ccd, "gs", "spline5-clamped").params["junction"] == "C2",
    )
    rms = fits.make_potential(ccd, "gs", "spline3-s0.05").rms() * 1e3
    # FITPACK meets its sum-of-squares target to a relative 1e-3, i.e. the rms to sqrt(1.001).
    c.check("smoothing spline meets its rms target", rms <= 0.05 * np.sqrt(1.001), "rms %.5f meV (target 0.05)" % rms)
    for bad in ("spline4-natural", "spline2-clamped"):
        try:
            fits.parse_spec(bad)
            c.check("%s is rejected (even order)" % bad, False, "accepted")
        except ValueError:
            c.check("%s is rejected (even order)" % bad, True)
    check_cubic_correction(c, ccd)


def check_cubic_correction(c, ccd):
    """The first-order cubic correction: units, b -> 0 limit, norms, order of accuracy and sign."""

    from scipy import constants

    print()
    print("Cubic anharmonic correction")
    sub = config_coord_1d_vof_solver
    w_gs, w_es, delta_q = ccd["gs_phonon"], ccd["es_phonon"], ccd["delta_q"]
    gs, es = fits.make_potential(ccd, "gs", "cubic"), fits.make_potential(ccd, "es", "cubic")
    c.check("cubic fit keeps the parabola's a (omega unchanged)", gs.params["a"] == ccd["gs_params"][0])
    c.close("cubic fit phonon_energy is the harmonic omega", gs.phonon_energy(), w_gs, 1e-12, " meV")

    # xi against an independent SI evaluation of b / (hbar Omega) (hbar / 2 Omega)^(3/2).
    amu = constants.physical_constants["atomic mass constant"][0]
    omega = w_gs * 1e-3 * constants.eV / constants.hbar
    b_si = gs.params["b"] * constants.eV / (amu**1.5 * 1e-30)
    xi_si = b_si / (constants.hbar * omega) * (constants.hbar / (2 * omega)) ** 1.5
    xi = sub.anharmonic_strength(gs.params["b"], w_gs)
    c.check("xi matches the SI formula", abs(xi - xi_si) <= 1e-12 * abs(xi_si), "xi = %.6e" % xi)

    def solve(scale, sign=1.0, n=12):
        s = sub(w_es, w_gs, delta_q)
        with contextlib.redirect_stdout(io.StringIO()):
            s.compute_franck_condon_integrals_cubic(
                n, n, sign * scale * es.params["b"], sign * scale * gs.params["b"], q_i=delta_q, q_f=0.0
            )
        return s

    def exact(scale, n=12):
        def well(p, center):
            return lambda q: p["a"] * (q - center) ** 2 + scale * p["b"] * (q - center) ** 3

        s = sub(w_es, w_gs, delta_q)
        s.compute_franck_condon_integrals_fd(n, n, well(es.params, delta_q), well(gs.params, 0.0), q_i=delta_q, q_f=0.0)
        return s.fc_ints[:8, :8]

    herm = sub(w_es, w_gs, delta_q)
    herm.compute_franck_condon_integrals_hermite(12, 12)
    dev = np.max(np.abs(solve(0.0).fc_ints - herm.fc_ints))
    c.check("b = 0 reproduces the Hermite route", dev < 1e-13, "max dev %.1e" % dev)

    full = solve(1.0, n=40)
    predicted = (sub.cubic_correction_coefficients(40, full.cubic_xi_f) ** 2).sum(axis=1)
    dev = np.max(np.abs(full.cubic_norms_f - predicted))
    c.check("norms equal 1 + sum of squared PT coefficients", dev < 1e-12, "max dev %.1e" % dev)
    c.check("norms are 1 + O(xi^2) at n = 0", abs(full.cubic_norms_f[0] - 1) < 10 * full.cubic_xi_f**2,
            "norm - 1 = %.2e, xi^2 = %.2e" % (full.cubic_norms_f[0] - 1, full.cubic_xi_f**2))  # fmt: skip

    err_1 = np.max(np.abs(solve(1.0).fc_ints[:8, :8] - exact(1.0)))
    err_half = np.max(np.abs(solve(0.5).fc_ints[:8, :8] - exact(0.5)))
    c.check("PT error vs the exact cubic falls as b^2", 3.5 < err_1 / err_half < 4.5, "ratio %.2f" % (err_1 / err_half))
    err_flip = np.max(np.abs(solve(1.0, sign=-1.0).fc_ints[:8, :8] - exact(1.0)))
    c.check("emission-frame sign of b is right (wrong sign is much worse)", err_flip > 10 * err_1,
            "%.1e vs %.1e" % (err_flip, err_1))  # fmt: skip

    captured = io.StringIO()
    with contextlib.redirect_stdout(captured):
        sub(w_es, w_gs, delta_q).compute_franck_condon_integrals_cubic(5, 5, 50.0, 50.0, q_i=delta_q, q_f=0.0)
    c.check("|xi| > 0.1 prints the warning", "|xi| =" in captured.getvalue() and "_fd" in captured.getvalue())


if __name__ == "__main__":
    sys.exit(main())
