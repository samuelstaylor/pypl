#!/usr/bin/env python3
"""The 1D CCD with its vibrational wavefunctions, and the Franck-Condon lineshape (VOF) on a side panel.

    python validation/plot_ccd_wavefunctions.py                                  # Hermite wavefunctions, PyPL VOF
    python validation/plot_ccd_wavefunctions.py --vof-method pypl hermite fd     # three VOFs overlaid, labelled
    python validation/plot_ccd_wavefunctions.py --wf-method fd --vof-method fd   # finite differences throughout
    python validation/plot_ccd_wavefunctions.py --wf-method fd --fd-fit spline --vof-method pypl fd
    python validation/plot_ccd_wavefunctions.py --n-levels 5 --well gs --density

Left panel: the DFT images, the fitted wells and the lowest ``--n-levels`` wavefunctions of each well,
each on its own energy line, drawn as a dimgray fill behind the curves (the yujin_isc_fig1.py style).
``--wf-method`` picks where the wavefunctions come from:

    hermite   the analytic harmonic eigenfunctions (config_coord_1d_vof_solver.harmonic_wavefunction),
              at E_min + (v + 1/2) hbar omega
    fd        finite-difference eigenstates of the fitted wells (solve_schrodinger_fd, --fd-fit), at their
              own eigenvalues

``config_coord_1d_solver`` never forms wavefunctions, so the drawn ones are *checked*: their overlaps
<chi_i | chi_j> are compared with the solver's own F_ij, and the script refuses to plot if they
disagree (1e-10 for hermite, 1e-8 for fd). An anharmonic --fd-fit has no harmonic counterpart, so
there the mismatch is reported and the plot is drawn anyway.

Right panel: ``bulid_fc_lsp`` with F_ij from each ``--vof-method`` (the first filled, with the
discrete p_i |F_ij|^2 as ticks; the others dashed), on the CCD's total-energy axis.

Output: ``validation/outputs/plot_ccd_wavefunctions/<run>/ccd_wavefunctions.png`` with log.txt and index.md,
where <run> names the choices, e.g. ``wf-hermite_vof-pypl`` or ``wf-fd_vof-pypl-fd_fit-spline``.
"""

import argparse
import contextlib
import dataclasses
import io
from functools import partial

import matplotlib
import numpy as np

from harness import ccd_data, fits, plots
from harness.ccd_data import STATE_LABELS, STATE_MARKERS, STATE_NAMES, STATES
from harness.methods import (
    DIRECTION_LABELS,
    DIRECTIONS,
    METHODS,
    add_vof_args,
    compute_vof,
    fd_q_range,
    params_from_args,
)
from harness.output import OutputDir
from pypl.config_coord_1d_vof_solver import config_coord_1d_vof_solver, grid_spacing

WF_COLOR = "dimgray"
# The reference VOF of each direction is a solid, filled line in the direction's colour. Every further method
# gets the direction's own (colour, linestyle, width), so emission and absorption never look alike.
VOF_TEST_STYLES = {
    "emission": (("#8B0000", "--", 1.6), ("#F4A261", ":", 1.8), ("#5D1A12", (0, (6, 2, 1, 2)), 1.4)),
    "absorption": (("#0B3D91", "-.", 1.6), ("#00ACC1", (0, (1, 1)), 1.8), ("#1A237E", (0, (8, 3)), 1.4)),
}


# ------------------------------------------------------------------------- wavefunctions


def harmonic_levels(ccd, n_levels):
    """``{state: [(energy eV, callable chi_v(q), (q_left, q_right)), ...]}`` from the Hermite functions."""

    out = {}
    for state in STATES:
        phonon, center = ccd["%s_phonon" % state], ccd_data.center(ccd, state)
        e_min = ccd["%s_params" % state][1]
        levels = []
        for v in range(n_levels):
            turn = config_coord_1d_vof_solver.classical_turning_point(v, phonon)
            chi = partial(config_coord_1d_vof_solver.harmonic_wavefunction, n=v, omega=phonon, center=center)
            levels.append((e_min + (v + 0.5) * phonon * 1e-3, chi, (center - turn, center + turn)))
        out[state] = levels
    return out


def fd_levels(ccd, n_levels, params):
    """The same, from finite-difference eigenstates of the fitted wells on one shared grid of Q.

    Returns ``(q, levels)``; the energies are total energies because the fits carry their offsets.
    """

    tp = {s: config_coord_1d_vof_solver.classical_turning_point(n_levels - 1, ccd["%s_phonon" % s]) for s in STATES}
    q = np.linspace(*ccd_data.q_span(ccd, tp, params.fd_margin), params.fd_npts)
    out = {}
    for state in STATES:
        pot = fits.make_potential(ccd, state, params.fd_fit)
        energies, psi = config_coord_1d_vof_solver.solve_schrodinger_fd(pot(q), q, n_levels, params.fd_stencil)
        levels = []
        for v in range(n_levels):
            allowed = q[pot(q) <= energies[v]]
            levels.append((energies[v], partial(np.interp, xp=q, fp=psi[v]), (allowed.min(), allowed.max())))
        out[state] = levels
    return q, out


def cubic_levels(ccd, n_levels, params):
    """The same, from the first-order corrected eigenfunctions of the cubic wells (``--wf-method cubic``).

    Each well's b comes from ``params.cubic_fit``; the frequency is the harmonic one and the level energies
    stay harmonic (the first-order shift of a cubic perturbation vanishes). Drawn on the CCD's own Q axis,
    so b keeps its CCD sign here.
    """

    out = {}
    q = np.linspace(-20.0, 20.0, 40001)
    for state in STATES:
        phonon, center = ccd["%s_phonon" % state], ccd_data.center(ccd, state)
        pot = fits.make_potential(ccd, state, params.cubic_fit)
        xi = config_coord_1d_vof_solver.anharmonic_strength(pot.params["b"], phonon)
        e_min = ccd["%s_params" % state][1]
        levels = []
        for v in range(n_levels):
            chi = partial(config_coord_1d_vof_solver.perturbed_wavefunction, n=v, omega=phonon, xi=xi, center=center)
            if params.cubic_normalize:
                norm = np.sqrt(np.sum(chi(q) ** 2) * grid_spacing(q))
                chi = (lambda f, c: lambda x: f(x) / c)(chi, norm)
            energy = e_min + (v + 0.5) * phonon * 1e-3
            allowed = q[(pot(q) <= energy) & (np.abs(q - center) < 5)]
            if allowed.size == 0:
                # The cubic lies above the harmonic level everywhere near the minimum (e.g. a fitted a < 0):
                # mark the harmonic turning points instead.
                turn = config_coord_1d_vof_solver.classical_turning_point(v, phonon)
                allowed = np.array([center - turn, center + turn])
            levels.append((energy, chi, (allowed.min(), allowed.max())))
        out[state] = levels
    return out


def check_against_solver(ccd, levels, q, n_levels, simpson=False, cubic=None):
    """Overlap the drawn ES and GS functions and compare with the solver's F_ij (emission, ES -> GS).

    The drawn functions sit at their physical minima (ES at +dQ) while the solver centres the initial
    state at -dQ; mirroring Q maps one onto the other and multiplies the overlap by (-1)^(i+j).
    ``simpson`` integrates with Simpson weights (for the analytic functions); FD eigenvectors are
    orthonormal under the plain sum, so they use that. With ``cubic`` (the VOFParams) the reference is the
    cubic method's own F_ij for the same b, so the check ties the drawn corrected functions to that
    method's basis. Returns ``(numerical, solver_block, max_dev)``.
    """

    solver = config_coord_1d_vof_solver(ccd["es_phonon"], ccd["gs_phonon"], ccd["delta_q"])
    if cubic is None:
        solver.compute_franck_condon_integrals(ni=n_levels + 2, nf=n_levels + 2)
    else:
        b = {s: fits.make_potential(ccd, s, cubic.cubic_fit).params["b"] for s in STATES}
        with contextlib.redirect_stdout(io.StringIO()):
            solver.compute_franck_condon_integrals_cubic(
                n_levels,
                n_levels,
                b["es"],
                b["gs"],
                q_i=ccd_data.center(ccd, "es"),
                q_f=ccd_data.center(ccd, "gs"),
                normalize=cubic.cubic_normalize,
            )
    block = solver.fc_ints[:n_levels, :n_levels]
    weights = np.full(len(q), grid_spacing(q))
    if simpson:
        weights[1:-1:2] *= 4.0 / 3.0
        weights[2:-1:2] *= 2.0 / 3.0
        weights[[0, -1]] /= 3.0
    es = np.array([chi(q) for _, chi, _ in levels["es"]])
    gs = np.array([chi(q) for _, chi, _ in levels["gs"]])
    # The solver frame has the initial (ES) well on the left. A CCD with the ES on the right (diamond) is
    # its mirror image, so overlaps of the drawn functions pick up (-1)^(i+j); with the ES on the left (MgO)
    # the two frames coincide.
    if ccd_data.center(ccd, "es") > ccd_data.center(ccd, "gs"):
        parity = (-1.0) ** np.add.outer(np.arange(n_levels), np.arange(n_levels))
    else:
        parity = np.ones((n_levels, n_levels))
    # np.sum's pairwise summation, not a BLAS matrix product, whose rounding on 2e5 terms is ~1e-12.
    numerical = parity * np.array([[np.sum(a * b * weights) for b in gs] for a in es])
    return numerical, block, float(np.max(np.abs(numerical - block)))


def report_check(out, numerical, block, max_dev, tol, source):
    out.log("\nWavefunction check: <chi_i | chi_j> of the %d drawn levels vs the solver's F_ij" % len(block))
    out.log("  from the drawn %s:" % source)
    for row in numerical:
        out.log("    " + "  ".join("% .9f" % v for v in row))
    out.log("  solver fc_ints:")
    for row in block:
        out.log("    " + "  ".join("% .9f" % v for v in row))
    out.log("  %s (max deviation %.2e, tolerance %.0e)" % ("MATCH" if max_dev <= tol else "MISMATCH", max_dev, tol))


# --------------------------------------------------------------------------------- panels


def draw_ccd(ax, ccd, args, params, levels):
    """Fitted wells, DFT images and wavefunctions. Returns the (y_lo, y_hi) the diagram wants."""

    delta_q = ccd["delta_q"]
    tp = {
        s: config_coord_1d_vof_solver.classical_turning_point(args.n_levels - 1, ccd["%s_phonon" % s]) for s in STATES
    }
    lo, hi = ccd_data.q_span(ccd, tp, 0.15)
    # Reach every DFT image too: every fit but the parabola is made from all of them.
    q_all = np.concatenate([ccd_data.samples(ccd, s)[0] for s in STATES])
    q = np.linspace(min(lo, q_all.min() - 0.15), max(hi, q_all.max() + 0.15), 2001)

    for state in STATES:
        curves = [chi(q) ** 2 if args.density else chi(q) for _, chi, _ in levels[state]]
        # One scale per well keeps the relative amplitudes physical.
        scale = args.wf_scale * ccd["%s_phonon" % state] * 1e-3 / max(np.max(np.abs(c)) for c in curves)
        for (energy, _, (left, right)), curve in zip(levels[state], curves):
            ax.fill_between(q, energy, energy + curve * scale, color=WF_COLOR, alpha=args.alpha, lw=0, zorder=1)
            ax.hlines(energy, left, right, color=WF_COLOR, lw=0.9, alpha=0.75, zorder=2)

    spec = {"fd": params.fd_fit, "cubic": params.cubic_fit}.get(args.wf_method, fits.DEFAULT_FIT)
    for state in ("es", "gs"):
        color = plots.STATE_COLORS[state]
        suffix = "" if spec == fits.DEFAULT_FIT else " (%s)" % fits.label(spec)
        ax.plot(
            q,
            fits.make_potential(ccd, state, spec)(q),
            color=color,
            lw=1.4,
            zorder=3,
            label="%s fit%s\n[%s]" % (STATE_LABELS[state], suffix, fits.images(ccd, state, spec)),
        )
        q_s, e_s = ccd_data.samples(ccd, state)
        ax.plot(
            q_s,
            e_s,
            color=color,
            ls="",
            marker=STATE_MARKERS[state],
            ms=5,
            zorder=4,
            label="%s DFT energy" % STATE_LABELS[state],
        )
    wf_label = r"$|\chi_v(Q)|^2$" if args.density else r"$\chi_v(Q)$"
    ax.fill_between(
        [],
        [],
        [],
        color=WF_COLOR,
        alpha=args.alpha,
        label="%s, v = 0..%d (%s)" % (wf_label, args.n_levels - 1, args.wf_method),
    )
    ax.axvline(0.0, color="gray", ls="--", lw=1, zorder=0)
    ax.axvline(delta_q, color="gray", ls="--", lw=1, zorder=0)

    if args.well in STATES:
        e_min, spacing = ccd["%s_params" % args.well][1], ccd["%s_phonon" % args.well] * 1e-3
        y = (e_min - 0.4 * spacing, e_min + (args.n_levels - 0.5 + 1.2 * args.wf_scale) * spacing)
    else:
        top = ccd["es_params"][1] + (args.n_levels - 0.5) * ccd["es_phonon"] * 1e-3
        y = (
            ccd["gs_params"][1] - 0.07,
            max(top + args.wf_scale * ccd["es_phonon"] * 1e-3, ccd["es_params"][1] + 0.3) + 0.07,
        )
    ax.set_xlim((q[0], q[-1]))
    ax.set_xlabel(r"Q (amu$^{1/2}$ $\AA$)")
    ax.set_ylabel("Total energy (eV)")
    plots.style(ax, grid=False)
    return y


def draw_vof(ax, results):
    """VOFs of every method, normalized to the tallest sideband peak, with the first method's ticks."""

    peak = max(r.sideband_max for group in results.values() for r in group)
    for direction, group in results.items():
        color = plots.DIRECTION_COLORS[direction]
        for k, r in enumerate(group):
            if k == 0:
                line_color, ls, lw = color, "-", 1.3
            else:
                styles = VOF_TEST_STYLES[direction]
                line_color, ls, lw = styles[min(k - 1, len(styles) - 1)]
            env = r.lineshape / peak
            if k == 0:
                ax.fill_betweenx(r.energy_ev, 0.0, env, color=color, alpha=0.25, lw=0, zorder=1)
            label = "%s\n%s" % (DIRECTION_LABELS[direction].replace("->", r"$\rightarrow$"), r.label)
            ax.plot(env, r.energy_ev, color=line_color, ls=ls, lw=lw, zorder=2 + k, label=label)
        first = group[0]
        keep = first.weights > 1e-4 * np.max(first.weights)
        sticks = first.weights[keep] / np.max(first.weights) * first.sideband_max / peak
        ax.hlines(first.levels_ev[keep], 0.0, sticks, color=color, lw=0.8, alpha=0.5, zorder=1.5)
    ax.set_xlim((0.0, 1.18))
    ax.set_xlabel("VOF (norm.)")
    ax.tick_params(labelleft=False)
    plots.style(ax)
    ax.grid(axis="y", visible=False)


# ----------------------------------------------------------------------------------- main


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    g = p.add_argument_group("wavefunctions")
    g.add_argument(
        "--wf-method",
        choices=("hermite", "fd", "cubic"),
        default="hermite",
        help="source of the drawn wavefunctions: hermite, fd (eigenstates of --fd-fit), cubic (first-order corrected "
        "for --cubic-points / --cubic-free-a)",
    )
    g.add_argument("--n-levels", type=int, default=3, help="levels drawn per well (default 3)")
    g.add_argument("--density", action="store_true", help="draw |chi_v|^2 instead of chi_v")
    g.add_argument("--wf-scale", type=float, default=0.375, help="peak amplitude / level spacing (default 0.375)")
    g.add_argument("--alpha", type=float, default=0.3, help="fill opacity (default 0.3)")
    g.add_argument("--well", choices=("both", "gs", "es"), default="both", help="zoom the energy axis onto one well")
    g = p.add_argument_group("VOF panel")
    g.add_argument("--no-vof", action="store_true", help="omit the VOF panel")
    g.add_argument(
        "--vof-method",
        nargs="+",
        choices=sorted(METHODS),
        default=["pypl"],
        metavar="METHOD",
        help="F_ij source(s): %s (default pypl)" % ", ".join(METHODS),
    )
    add_vof_args(p)
    p.add_argument("--tag", default="", help="suffix of the run directory")
    p.add_argument("-v", "--verbose", action="store_true", help="print every method's settings")
    p.add_argument("--show", action="store_true", help="display the figure as well as saving it")
    return p.parse_args()


def main():
    args = parse_args()
    if args.n_levels < 1:
        raise SystemExit("--n-levels must be at least 1")
    params = params_from_args(args)
    if not args.show:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    uses_fd = args.wf_method == "fd" or (not args.no_vof and "fd" in args.vof_method)
    anharmonic = uses_fd and params.fd_fit != fits.DEFAULT_FIT
    # Every run that uses finite differences names its fit, the default parabola included.
    run = "wf-%s_%s" % (args.wf_method, "no-vof" if args.no_vof else "vof-" + "-".join(args.vof_method))
    uses_cubic = args.wf_method == "cubic" or (not args.no_vof and "cubic" in args.vof_method)
    anharmonic = anharmonic or uses_cubic
    run += "_fit-%s" % params.fd_fit if uses_fd else ""
    run += (
        ("_cubicfit-%s%s" % (params.cubic_fit, "-normalized" if params.cubic_normalize else "")) if uses_cubic else ""
    )
    run += "_%s" % args.tag if args.tag else ""
    ccd = ccd_data.build_ccd(verbose=False)
    if uses_fd:
        try:
            fits.check_bounded(ccd, params.fd_fit, fd_q_range(ccd, params))
        except ValueError as err:
            raise SystemExit("Cannot solve on --fd-fit %s: %s" % (params.fd_fit, err))
    out = OutputDir(
        "plot_ccd_wavefunctions",
        run,
        title="1D CCD with wavefunctions: %s" % run,
        verbose=args.verbose,
        keep_figures=args.show,
    )
    out.save_json("settings.json", dict(vars(args), params=dataclasses.asdict(params)))
    if uses_fd:
        out.log("Finite differences on fit %s (%s):" % (params.fd_fit, fits.label(params.fd_fit)))
        for state in STATES:
            out.log("  %s well: %s" % (STATE_NAMES[state], fits.make_potential(ccd, state, params.fd_fit).describe()))

    if args.wf_method == "hermite":
        levels = harmonic_levels(ccd, args.n_levels)
        tp = max(
            config_coord_1d_vof_solver.classical_turning_point(args.n_levels, ccd["%s_phonon" % s]) for s in STATES
        )
        q_check = np.linspace(*ccd_data.q_span(ccd, {"gs": tp, "es": tp}, 6.0), 200001)
        tol, source = 1e-10, "Hermite functions"
    elif args.wf_method == "cubic":
        levels = cubic_levels(ccd, args.n_levels, params)
        solver = config_coord_1d_vof_solver(ccd["es_phonon"], ccd["gs_phonon"], ccd["delta_q"])
        q_check, _ = solver._simpson_grid(args.n_levels + 2, args.n_levels + 2, 100001)
        tol = 1e-10
        source = "first-order corrected functions (%s), vs the cubic method" % fits.label(params.cubic_fit)
        for state in STATES:
            pot = fits.make_potential(ccd, state, params.cubic_fit)
            xi = config_coord_1d_vof_solver.anharmonic_strength(pot.params["b"], ccd["%s_phonon" % state])
            out.log("  %s well: %s; xi = %+.3e" % (STATE_NAMES[state], pot.describe(), xi))
    else:
        q_check, levels = fd_levels(ccd, args.n_levels, params)
        tol, source = 1e-8, "finite-difference eigenstates (%s)" % fits.label(params.fd_fit)
    simpson = args.wf_method in ("hermite", "cubic")
    cubic = params if args.wf_method == "cubic" else None
    numerical, block, max_dev = check_against_solver(ccd, levels, q_check, args.n_levels, simpson, cubic)
    report_check(out, numerical, block, max_dev, tol, source)
    if max_dev > tol:
        if not anharmonic or args.wf_method == "cubic":
            raise SystemExit("The drawn wavefunctions disagree with the solver; refusing to plot.")
        out.log(
            "  fit %s is anharmonic: a mismatch with the harmonic solver is expected; plotting anyway" % params.fd_fit
        )

    out.log("\nVibrational levels (eV)")
    for v in range(args.n_levels):
        out.log("  v = %d   GS % .6f   ES % .6f" % (v, levels["gs"][v][0], levels["es"][v][0]))

    if args.no_vof:
        fig, ax = plt.subplots(figsize=(6.6, 8.6))
    else:
        fig, (ax, ax_vof) = plt.subplots(
            1, 2, figsize=(8.8, 8.6), sharey=True, gridspec_kw={"width_ratios": [3, 1], "wspace": 0.05}
        )
    y_lo, y_hi = draw_ccd(ax, ccd, args, params, levels)
    subtitle = ""
    if not args.no_vof:
        results = {d: [compute_vof(ccd, d, m, params, log=out.detail) for m in args.vof_method] for d in DIRECTIONS}
        out.log("\nFranck-Condon lineshapes via bulid_fc_lsp")
        for direction, group in results.items():
            for r in group:
                out.log(
                    "  %-9s %-36s closure %.10f  sideband peak %.4f eV  centroid %.1f meV"
                    % (direction, r.label, r.closure, r.sideband_peak_ev, r.centroid_mev)
                )
            for r in group[1:]:
                dev = np.max(np.abs(r.lineshape - group[0].lineshape)) / np.max(group[0].lineshape)
                out.log("    %s vs %s: max |difference| / peak = %.2e" % (r.label, group[0].label, dev))
        draw_vof(ax_vof, results)
        subtitle = "\nVOF: " + "  |  ".join(r.label for r in results["emission"])
        if args.well == "both":
            absorption = results["absorption"][0]
            tail = absorption.lineshape > 0.03 * absorption.sideband_max
            y_hi = max(y_hi, float(absorption.energy_ev[tail].max()) + 0.04)
    ax.set_ylim((y_lo, y_hi))
    if args.no_vof:
        ax.legend(loc="center left" if args.well == "both" else "upper right", edgecolor="black", fontsize=9)
    else:
        # Both legends go below the panels, where they cannot cover the curves: CCD entries left, VOF entries right.
        h_ccd, l_ccd = ax.get_legend_handles_labels()
        h_vof, l_vof = ax_vof.get_legend_handles_labels()
        l_vof = ["VOF, " + label.replace("\n", ": ") for label in l_vof]
        pad = max(len(l_ccd), len(l_vof))
        blank = plt.Line2D([], [], color="none")
        fig.legend(
            h_ccd + [blank] * (pad - len(l_ccd)) + h_vof + [blank] * (pad - len(l_vof)),
            l_ccd + [""] * (pad - len(l_ccd)) + l_vof + [""] * (pad - len(l_vof)),
            loc="upper center",
            bbox_to_anchor=(0.5, 0.0),
            ncol=2,
            fontsize=8,
            edgecolor="black",
        )
    wf_source = {
        "hermite": "Hermite",
        "fd": "finite difference on %s" % fits.label(params.fd_fit),
        "cubic": "Hermite + cubic PT, %s" % fits.label(params.cubic_fit),
    }[args.wf_method]
    fit_line = ""
    if uses_cubic:
        b = {s: fits.make_potential(ccd, s, params.cubic_fit).params["b"] for s in STATES}
        xi = {s: config_coord_1d_vof_solver.anharmonic_strength(b[s], ccd["%s_phonon" % s]) for s in STATES}
        fit_line += "\ncubic: b GS %+.4f, ES %+.4f eV/(amu$^{3/2}$ $\\AA^3$);  $\\xi$ GS %+.1e, ES %+.1e" % (
            b["gs"],
            b["es"],
            xi["gs"],
            xi["es"],
        )
    if uses_fd:
        hw = {s: fits.make_potential(ccd, s, params.fd_fit).phonon_energy() for s in STATES}
        fit_line = "\nfd wells: %s (hbar w at the minimum: %s %.2f, %s %.2f meV)" % (
            fits.label(params.fd_fit),
            STATE_NAMES["gs"],
            hw["gs"],
            STATE_NAMES["es"],
            hw["es"],
        )
    fig.suptitle(
        "%s, 1D CCD, wavefunctions: %s\n"
        r"harmonic %s %.2f meV,  %s %.2f meV,  $\Delta Q$ = %.3f amu$^{1/2}\,\AA$%s%s"
        % (
            ccd_data.active().title,
            wf_source,
            STATE_NAMES["gs"],
            ccd["gs_phonon"],
            STATE_NAMES["es"],
            ccd["es_phonon"],
            ccd["delta_q"],
            fit_line,
            subtitle,
        ),
        fontsize=12,
    )
    out.save_fig(
        fig,
        "ccd_wavefunctions.png",
        "The CCD with %d wavefunctions per well (%s)%s."
        % (args.n_levels, wf_source, "" if args.no_vof else " and the VOF panel"),
    )
    out.add_summary(
        "Wavefunctions: %s. VOF methods: %s." % (wf_source, "none" if args.no_vof else ", ".join(args.vof_method))
    )
    if uses_fd:
        out.add_summary("", "Fit **%s** (%s):" % (params.fd_fit, fits.label(params.fd_fit)))
        out.add_summary(
            *["- %s well: %s" % (STATE_NAMES[s], fits.make_potential(ccd, s, params.fd_fit).describe()) for s in STATES]
        )
    out.add_summary(
        "",
        "Wavefunction check: max |<chi_i|chi_j> - F_ij| = %.2e (tolerance %.0e)%s."
        % (max_dev, tol, "; mismatch expected, the fit is anharmonic" if anharmonic and max_dev > tol else ""),
    )
    out.close()
    if args.show:
        plt.show()


if __name__ == "__main__":
    main()
