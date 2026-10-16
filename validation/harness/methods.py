"""Franck-Condon lineshape (VOF) methods, shared parameters and command-line flags.

A VOF is always ``bulid_fc_lsp`` applied to a Franck-Condon matrix F_ij. A *method* is only the
source of F_ij; every method runs on a ``config_coord_1d_vof_solver``, which is PyPL's
``config_coord_1d_solver`` with extra methods and nothing changed:

    pypl      compute_franck_condon_integrals           (inherited: Ruhoff recurrences)
    hermite   compute_franck_condon_integrals_hermite   (quadrature of harmonic eigenfunctions)
    fd        compute_franck_condon_integrals_fd        (finite differences on fitted wells, see fits.py)

``compute_vof(ccd, direction, method, params)`` returns a ``VOFResult`` carrying the solver, the
lineshape, timing and ``notes``: every setting the method used, as printable lines.
"""

import contextlib
import io
import time
from dataclasses import dataclass, field
from typing import Callable

import numpy as np

from harness import ccd_data, fits
from pypl.config_coord_1d_vof_solver import FD_STENCIL_NAMES, FD_STENCILS, KINETIC, config_coord_1d_vof_solver

DIRECTIONS = ("emission", "absorption")
# Follows the active dataset: "ES -> GS" for diamond, "3E -> 1A1" for MgO.
DIRECTION_LABELS = ccd_data.DIRECTION_LABELS


def direction_states(direction):
    """``(initial, final)`` electronic states of a direction."""

    return ("es", "gs") if direction == "emission" else ("gs", "es")


# ---------------------------------------------------------------------------- parameters


@dataclass(frozen=True)
class VOFParams:
    """Every parameter of a VOF calculation. Defaults reproduce tutorial 002."""

    temp: float = 5.0  # K
    sigma: tuple = (8.0, 25.0)  # meV, [sigma_min, sigma_max]; widens linearly with |E|
    gamma: float = 0.3  # meV, Lorentzian ZPL width
    zpl_lorentzian: bool = True  # False: the ZPL keeps its Gaussian
    order_es: int = 50  # ES vibrational levels
    order_gs: int = 60  # GS vibrational levels
    ene_range: tuple = (-200.0, 1200.0)  # meV, relative to the ZPL
    resol: int = 1401  # axis points (1 meV steps resolve the 0.3 meV ZPL)
    npts: int = 100001  # quadrature points (hermite)
    fd_npts: int = 20001  # grid points (fd)
    fd_stencil: int = 5  # 3- or 5-point second derivative (fd)
    fd_margin: float = 2.0  # amu^0.5 A of grid beyond the outermost turning point (fd)
    fd_fit: str = fits.DEFAULT_FIT  # fitted wells (fd), a fits.py spec
    cubic_fit: str = "cubic"  # cubic fit giving b (cubic method): cubic[-near][-free]
    cubic_normalize: bool = False  # renormalize the first-order states (cubic method)

    def axis(self):
        return np.linspace(self.ene_range[0], self.ene_range[1], self.resol)

    def orders(self, direction):
        """``(ni, nf)``: the basis sizes swap with the direction, as in tutorial 002."""

        return (self.order_es, self.order_gs) if direction == "emission" else (self.order_gs, self.order_es)

    def describe(self):
        """The settings shared by every method, as printable lines."""

        zpl = "Lorentzian, gamma = %g meV" % self.gamma if self.zpl_lorentzian else "Gaussian (no Lorentzian ZPL)"
        axis = self.axis()
        return [
            "bulid_fc_lsp: T = %g K, sigma = %s meV, ZPL %s" % (self.temp, list(self.sigma), zpl),
            "lineshape axis: %g..%g meV from the ZPL, %d points, step %.4g meV"
            % (axis[0], axis[-1], len(axis), axis[1] - axis[0]),
            "basis: %d ES x %d GS levels (emission ni x nf = %d x %d, absorption %d x %d)"
            % (self.order_es, self.order_gs, self.order_es, self.order_gs, self.order_gs, self.order_es),
        ]


def add_vof_args(parser):
    """Add the VOFParams flags to ``parser``. Pair with ``params_from_args``."""

    d = VOFParams()
    ccd_data.add_dataset_arg(parser)
    g = parser.add_argument_group("VOF parameters (identical in every script)")
    g.add_argument("--temp", "--vof-temp", dest="temp", type=float, default=d.temp, help="temperature in K")
    g.add_argument(
        "--sigma",
        "--vof-sigma",
        dest="sigma",
        type=float,
        nargs="+",
        default=list(d.sigma),
        metavar="SIGMA",
        help="Gaussian broadening in meV: one value is constant, two are the solver's [sigma_min, sigma_max]",
    )
    g.add_argument(
        "--gamma", "--vof-zpl-gamma", dest="gamma", type=float, default=d.gamma, help="Lorentzian ZPL width in meV"
    )
    g.add_argument(
        "--no-zpl", "--no-vof-zpl", dest="no_zpl", action="store_true", help="broaden the ZPL with the Gaussian"
    )
    g.add_argument("--order-es", type=int, default=d.order_es, help="ES vibrational levels")
    g.add_argument("--order-gs", type=int, default=d.order_gs, help="GS vibrational levels")
    g.add_argument(
        "--ene-range", type=float, nargs=2, default=list(d.ene_range), metavar=("LO", "HI"), help="axis in meV"
    )
    g.add_argument("--resol", type=int, default=d.resol, help="points on that axis")
    g.add_argument("--npts", type=int, default=d.npts, help="quadrature points (hermite method)")
    g.add_argument("--fd-npts", type=int, default=d.fd_npts, help="grid points (fd method)")
    g.add_argument(
        "--fd-stencil", type=int, choices=sorted(FD_STENCILS), default=d.fd_stencil, help="d^2/dQ^2 stencil (fd)"
    )
    g.add_argument("--fd-margin", type=float, default=d.fd_margin, help="grid padding beyond the turning points (fd)")
    g.add_argument(
        "--fd-fit",
        default=d.fd_fit,
        metavar="FIT",
        help="fitted wells for the fd method: " + "; ".join("%s = %s" % kv for kv in fits.FIT_HELP.items()),
    )
    s = parser.add_argument_group("spline options (with --fd-fit spline; they build the spline spec)")
    s.add_argument(
        "--spline-order", type=int, choices=(2, 3, 4, 5), help="degree of the pieces: 3 cubic (default), 4, 5"
    )
    s.add_argument(
        "--spline-bc",
        choices=("natural", "not-a-knot", "clamped"),
        help="end conditions: natural (default for odd order), not-a-knot (default for even), clamped = slope "
        "(and for order 5 curvature) matched to the harmonic tails; natural and clamped need an odd order",
    )
    s.add_argument(
        "--spline-tails",
        choices=("harmonic", "smooth"),
        help="continuation past the end images: harmonic (default) or smooth = C2, the spline's own value, slope "
        "and curvature (not with clamped)",
    )
    s.add_argument(
        "--spline-smoothing",
        type=float,
        metavar="RMS_MEV",
        help="smoothing spline with rms residual <= RMS_MEV meV instead of exact interpolation (no end conditions)",
    )
    c = parser.add_argument_group("cubic perturbation options (method cubic)")
    c.add_argument(
        "--cubic-points",
        choices=("all", "parabola"),
        default="all",
        help="images the cubic b is fitted to: all (default) or only the parabola's near-minimum ones",
    )
    c.add_argument(
        "--cubic-free-a",
        action="store_true",
        help="fit a too; by default a is fixed to the parabola's, so the unperturbed Hamiltonian is the harmonic "
        "one whose frequency every method uses",
    )
    c.add_argument(
        "--cubic-normalize", action="store_true", help="renormalize the first-order states (default: norm 1 + O(xi^2))"
    )
    return g


def cubic_spec(args):
    """The cubic fit spec from ``--cubic-points`` / ``--cubic-free-a``."""

    return "cubic" + ("-near" if args.cubic_points == "parabola" else "") + ("-free" if args.cubic_free_a else "")


def fd_q_range(ccd, params):
    """CCD Q range of the fd grid of a VOF solve (the same for emission and absorption).

    The grid spans the initial well's outermost harmonic turning point to the final well's, plus the
    margin: on the CCD, from the left turning point of the well at the smaller Q to the right turning point
    of the other (diamond: GS left, ES right; MgO: 3E left, 1A1 right).
    """

    tp = config_coord_1d_vof_solver.classical_turning_point
    half = {"gs": tp(params.order_gs - 1, ccd["gs_phonon"]), "es": tp(params.order_es - 1, ccd["es_phonon"])}
    return ccd_data.q_span(ccd, half, params.fd_margin)


def resolve_fd_fit(args):
    """The canonical fit spec from ``--fd-fit`` and the ``--spline-*`` flags (which refine a spline spec)."""

    try:
        kind, options = fits.parse_spec(args.fd_fit)
        flags = (args.spline_order, args.spline_bc, args.spline_smoothing, args.spline_tails)
        if all(f is None for f in flags):
            return fits.canonical(args.fd_fit)
        if kind != "spline":
            raise ValueError("--spline-order / --spline-bc / --spline-smoothing / --spline-tails need --fd-fit spline")
        order = args.spline_order or options["order"]
        smoothing = args.spline_smoothing if args.spline_smoothing is not None else options["smoothing_mev"]
        bc = args.spline_bc
        if bc is None and smoothing is None and args.spline_order is None:
            bc = options["bc"]
        return fits.spline_spec(order, bc, smoothing, args.spline_tails or options["tails"])
    except ValueError as err:
        raise SystemExit("Invalid fit: %s" % err)


def params_from_args(args):
    """Build a VOFParams from the namespace produced by a parser that called ``add_vof_args``."""

    if len(args.sigma) not in (1, 2):
        raise SystemExit("--sigma takes one or two values")
    ccd_data.use_dataset(args.dataset)
    args.fd_fit = resolve_fd_fit(args)
    return VOFParams(
        temp=args.temp,
        sigma=(args.sigma[0], args.sigma[-1]),
        gamma=args.gamma,
        zpl_lorentzian=not args.no_zpl,
        order_es=args.order_es,
        order_gs=args.order_gs,
        ene_range=tuple(args.ene_range),
        resol=args.resol,
        npts=args.npts,
        fd_npts=args.fd_npts,
        fd_stencil=args.fd_stencil,
        fd_margin=args.fd_margin,
        fd_fit=args.fd_fit,
        cubic_fit=cubic_spec(args),
        cubic_normalize=args.cubic_normalize,
    )


# -------------------------------------------------------------------------------- methods


def harmonic_source_lines(ccd, direction):
    """Where a harmonic method's omega_i and omega_f come from: the two parabola fits."""

    return [
        "omega_%s = hbar sqrt(2a / amu) = %.5f meV, %s well [%s]: %s"
        % (
            role,
            ccd["%s_phonon" % state],
            ccd_data.STATE_NAMES[state],
            fits.images(ccd, state, "parabola"),
            ccd["%s_potential" % state].describe(),
        )
        for role, state in zip(("i", "f"), direction_states(direction))
    ]


def run_pypl(solver, ni, nf, params, ccd, direction):
    solver.compute_franck_condon_integrals(ni=ni, nf=nf)
    return [
        "method: inherited config_coord_1d_solver.compute_franck_condon_integrals, Ruhoff recurrences "
        "(Chem. Phys. 186, 355 (1994)), float64; no grid",
    ] + harmonic_source_lines(ccd, direction)


def run_hermite(solver, ni, nf, params, ccd, direction):
    solver.compute_franck_condon_integrals_hermite(ni, nf, npts=params.npts)
    q = solver.hermite_grid
    return [
        "method: compute_franck_condon_integrals_hermite, <chi_i(Q + dQ) | chi_j(Q)> of harmonic eigenfunctions",
        "grid: %d points, Q = %.4f..%.4f amu^0.5 A, h = %.3e amu^0.5 A; composite Simpson"
        % (len(q), q[0], q[-1], q[1] - q[0]),
    ] + harmonic_source_lines(ccd, direction)


def run_fd(solver, ni, nf, params, ccd, direction):
    fits.check_bounded(ccd, params.fd_fit, fd_q_range(ccd, params))
    initial, final = direction_states(direction)
    pot_i, pot_f = fits.make_potential(ccd, initial, params.fd_fit), fits.make_potential(ccd, final, params.fd_fit)
    solver.compute_franck_condon_integrals_fd(
        ni,
        nf,
        pot_i,
        pot_f,
        q_i=ccd_data.center(ccd, initial),
        q_f=ccd_data.center(ccd, final),
        npts=params.fd_npts,
        stencil=params.fd_stencil,
        margin=params.fd_margin,
    )
    info = solver.fd_info
    lines = [
        "method: compute_franck_condon_integrals_fd, H = -(hbar^2/2) d^2/dQ^2 + V(Q), hbar^2/2 = %.6e eV amu A^2"
        % KINETIC,
        "fit: %s (%s)" % (params.fd_fit, fits.label(params.fd_fit)),
        "potential (initial, %s) [%s]: %s"
        % (ccd_data.STATE_NAMES[initial], fits.images(ccd, initial, params.fd_fit), pot_i.describe()),
        "potential (final, %s) [%s]: %s"
        % (ccd_data.STATE_NAMES[final], fits.images(ccd, final, params.fd_fit), pot_f.describe()),
        "grid: %d points, frame %s, x = %.4f..%.4f (Q = %.4f..%.4f amu^0.5 A), h = %.4e amu^0.5 A"
        % (info["npts"], info["frame"], *info["x_range"], *info["q_range"], info["h"]),
        "padding: %.2f amu^0.5 A beyond the harmonic turning points of v = %d (%.4f) and v = %d (%.4f)"
        % (info["margin"], ni - 1, info["turning_points"][0], nf - 1, info["turning_points"][1]),
        "stencil: %s; hard walls" % FD_STENCIL_NAMES[info["stencil"]],
        "eigensolver: scipy.sparse.linalg.eigsh, shift-invert at min(V) - 1 eV, lowest %d and %d states" % (ni, nf),
        "checks: |psi| at walls / peak = %.2e, orthonormality = %.2e, ladder |E_v - E_0 - v hbar w| = %.2e (i) "
        "/ %.2e (f) meV" % (info["wall"], info["orthonormality"], info["ladder_i"], info["ladder_f"]),
    ]
    if info["wall"] > 1e-8:
        lines.append("WARNING: eigenstates reach the grid walls; the fit has too few bound levels for ni / nf")
    return lines


def run_cubic(solver, ni, nf, params, ccd, direction):
    initial, final = direction_states(direction)
    pot_i, pot_f = fits.make_potential(ccd, initial, params.cubic_fit), fits.make_potential(
        ccd, final, params.cubic_fit
    )
    # The solver prints its validity warnings; capture them so they become labelled notes.
    captured = io.StringIO()
    with contextlib.redirect_stdout(captured):
        solver.compute_franck_condon_integrals_cubic(
            ni,
            nf,
            pot_i.params["b"],
            pot_f.params["b"],
            q_i=ccd_data.center(ccd, initial),
            q_f=ccd_data.center(ccd, final),
            npts=params.npts,
            normalize=params.cubic_normalize,
        )
    info, q = solver.cubic_info, solver.cubic_grid
    lines = [
        "method: compute_franck_condon_integrals_cubic, harmonic eigenfunctions + first-order correction for "
        "H' = b Q^3 (omega and HR factors unchanged; level energies stay harmonic)",
        "cubic fit: %s (%s)" % (params.cubic_fit, fits.label(params.cubic_fit)),
        "potential (initial, %s) [%s]: %s"
        % (ccd_data.STATE_NAMES[initial], fits.images(ccd, initial, params.cubic_fit), pot_i.describe()),
        "potential (final, %s) [%s]: %s"
        % (ccd_data.STATE_NAMES[final], fits.images(ccd, final, params.cubic_fit), pot_f.describe()),
        "solver frame: %s, so b_i = %+.6f, b_f = %+.6f eV/(amu^1.5 A^3)"
        % ("as on the CCD" if info["frame_sign"] > 0 else "mirrored (both b change sign)", info["b_i"], info["b_f"]),
        "anharmonic strength: xi_i = %+.3e, xi_f = %+.3e; largest mixing 3|xi|(n+1)^1.5 = %.3g (i) / %.3g (f); "
        "below 0.1 up to n = %d (i) / %d (f)"
        % (
            info["xi_i"],
            info["xi_f"],
            info["max_mixing_i"],
            info["max_mixing_f"],
            info["n_valid_i"],
            info["n_valid_f"],
        ),
        "norms <psi~|psi~> - 1: initial %.2e..%.2e, final %.2e..%.2e%s"
        % (
            solver.cubic_norms_i.min() - 1,
            solver.cubic_norms_i.max() - 1,
            solver.cubic_norms_f.min() - 1,
            solver.cubic_norms_f.max() - 1,
            " (states renormalized)" if params.cubic_normalize else " (not renormalized)",
        ),
        "grid: %d points, Q = %.4f..%.4f amu^0.5 A, h = %.3e; composite Simpson; harmonic levels up to n + 3"
        % (len(q), q[0], q[-1], q[1] - q[0]),
    ] + harmonic_source_lines(ccd, direction)
    return lines + captured.getvalue().splitlines()


@dataclass(frozen=True)
class VOFMethod:
    """A source of F_ij: ``run(solver, ni, nf, params, ccd, direction)`` fills the solver and returns notes.

    ``anharmonic`` methods solve a different (anharmonic) problem, so they are compared with the harmonic
    methods only on request and their differences are physics, not error.
    """

    name: str
    label: str
    run: Callable
    anharmonic: bool = False


METHODS = {
    m.name: m
    for m in (
        VOFMethod("pypl", "PyPL (Ruhoff)", run_pypl),
        VOFMethod("hermite", "Hermite (quadrature)", run_hermite),
        VOFMethod("fd", "Finite difference", run_fd),
        VOFMethod("cubic", "Hermite + cubic PT", run_cubic, anharmonic=True),
    )
}


def is_anharmonic(method, params):
    """Whether ``method`` with ``params`` departs from the harmonic wells the reference uses."""

    return METHODS[method].anharmonic or (method == "fd" and params.fd_fit != fits.DEFAULT_FIT)


def method_label(method, params):
    """Method label, naming the fit for a non-default fd fit and the options of the cubic method."""

    label = METHODS[method].label
    if method == "fd" and params.fd_fit != fits.DEFAULT_FIT:
        label += " [%s]" % fits.label(params.fd_fit)
    if method == "cubic":
        options = [fits.label(params.cubic_fit)] + (["normalized"] if params.cubic_normalize else [])
        label += " [%s]" % ", ".join(options)
    return label


# --------------------------------------------------------------------------------- result


@dataclass
class VOFResult:
    """One VOF: the solver that produced it, its lineshape and how it was made."""

    method: str
    label: str
    direction: str
    solver: config_coord_1d_vof_solver
    e_zpl_mev: float  # total energy of the final state's v = 0 level
    axis: np.ndarray  # meV, relative to the ZPL
    lineshape: np.ndarray
    fit: str = None  # the fits.py spec of the wells (fd only)
    notes: list = field(default_factory=list)
    seconds: float = 0.0  # wall time building F_ij

    @property
    def omega_i(self):
        return self.solver.omega_i

    @property
    def omega_f(self):
        return self.solver.omega_f

    @property
    def energy_ev(self):
        return (self.axis + self.e_zpl_mev) * 1e-3

    @property
    def weights(self):
        """The discrete p_i |F_ij|^2, flattened."""

        return (self.solver.prefactors[:, None] * self.solver.fc_ints**2).ravel()

    @property
    def levels_ev(self):
        return (self.solver.energy_v.ravel() + self.e_zpl_mev) * 1e-3

    @property
    def closure(self):
        return float(np.sum(self.weights))

    @property
    def integral(self):
        return float(np.sum(self.lineshape) * (self.axis[1] - self.axis[0]))

    @property
    def _sideband(self):
        return self.axis > 0.5 * self.omega_f

    @property
    def sideband_max(self):
        return float(np.max(self.lineshape[self._sideband]))

    @property
    def centroid_mev(self):
        """First moment of the lineshape (meV from the ZPL): the mean phonon energy released or absorbed."""

        return float(np.sum(self.axis * self.lineshape) / np.sum(self.lineshape))

    @property
    def sideband_peak_ev(self):
        sb = self._sideband
        return float((self.axis[sb][np.argmax(self.lineshape[sb])] + self.e_zpl_mev) * 1e-3)


def compute_vof(ccd, direction, method, params, log=print):
    """Franck-Condon lineshape for one direction with F_ij from ``method`` (a METHODS key).

    Emission is ES -> GS, absorption GS -> ES. ``log`` receives every note, prefixed with the
    direction and method; pass ``None`` to stay silent.
    """

    if direction not in DIRECTIONS:
        raise ValueError("direction must be one of %s" % (DIRECTIONS,))
    if method not in METHODS:
        raise ValueError("unknown VOF method %r; choose from %s" % (method, sorted(METHODS)))
    initial, final = direction_states(direction)
    ni, nf = params.orders(direction)

    solver = config_coord_1d_vof_solver(ccd["%s_phonon" % initial], ccd["%s_phonon" % final], ccd["delta_q"])
    start = time.perf_counter()
    notes = METHODS[method].run(solver, ni, nf, params, ccd, direction)
    seconds = time.perf_counter() - start

    # bulid_fc_lsp prints its own "Integral check"; capture it so it can carry a label.
    captured = io.StringIO()
    with contextlib.redirect_stdout(captured):
        axis, lineshape = solver.bulid_fc_lsp(
            eneaxis=params.axis(),
            temp=params.temp,
            sigma=list(params.sigma),
            zpl_lorentzian=params.zpl_lorentzian,
            gamma=params.gamma,
        )
    notes = notes + ["F_ij built in %.2f s" % seconds] + captured.getvalue().splitlines()
    label = method_label(method, params)
    if log is not None:
        for line in notes:
            log("  [%s, %s] %s" % (DIRECTION_LABELS[direction], label, line))

    e_min_f = ccd["%s_params" % final][1]
    return VOFResult(
        method=method,
        label=label,
        direction=direction,
        solver=solver,
        e_zpl_mev=e_min_f * 1e3 + 0.5 * solver.omega_f,
        axis=axis,
        lineshape=lineshape,
        fit={"fd": params.fd_fit, "cubic": params.cubic_fit}.get(method),
        notes=notes,
        seconds=seconds,
    )


def compare_results(ref, tst, tol_fc=1e-7, tol_vof=1e-8):
    """Error statistics of ``tst`` against ``ref`` (two VOFResults of the same direction)."""

    fc_dev = np.abs(tst.solver.fc_ints - ref.solver.fc_ints)
    i_max, j_max = np.unravel_index(np.argmax(fc_dev), fc_dev.shape)
    diff = tst.lineshape - ref.lineshape
    peak = float(np.max(ref.lineshape))
    stats = {
        "fc_max": float(np.max(fc_dev)),
        "fc_at": (int(i_max), int(j_max)),
        "fc_rms": float(np.sqrt(np.mean(fc_dev**2))),
        "fc_low": float(np.max(fc_dev[:30, :30])),
        "vof_max": float(np.max(np.abs(diff))) / peak,
        "vof_rms": float(np.sqrt(np.mean(diff**2))) / peak,
        "zpl_weight_ref": float(ref.solver.fc_ints[0, 0] ** 2),
        "zpl_weight": float(tst.solver.fc_ints[0, 0] ** 2),
        "centroid_ref_mev": ref.centroid_mev,
        "centroid_shift_mev": tst.centroid_mev - ref.centroid_mev,
    }
    stats["direction"] = ref.direction
    stats["reference"] = ref.label
    stats["test"] = tst.label
    stats["fc_at_text"] = "(%d, %d)" % stats["fc_at"]
    stats["closure"] = tst.closure
    stats["integral"] = tst.integral
    stats["seconds"] = tst.seconds
    stats["fc_ok"] = stats["fc_max"] <= tol_fc
    stats["vof_ok"] = stats["vof_max"] <= tol_vof
    stats["ok"] = stats["fc_ok"] and stats["vof_ok"]
    stats["verdict"] = "OK" if stats["ok"] else "FAIL"
    return stats


# Columns of a method-comparison table: (key in compare_results' dict, header, format).
COMPARISON_COLUMNS = [
    ("direction", "direction", "%s"),
    ("test", "test", "%s"),
    ("reference", "reference", "%s"),
    ("fc_max", "max dF_ij", "%.2e"),
    ("fc_at_text", "at (i, j)", "%s"),
    ("fc_low", "max dF_ij (i,j<30)", "%.2e"),
    ("vof_max", "max dVOF / peak", "%.2e"),
    ("vof_rms", "rms dVOF / peak", "%.2e"),
    ("seconds", "F_ij time (s)", "%.2f"),
    ("verdict", "verdict", "%s"),
]
CSV_KEYS = (
    "direction",
    "reference",
    "test",
    "fc_max",
    "fc_at_text",
    "fc_rms",
    "fc_low",
    "vof_max",
    "vof_rms",
    "zpl_weight_ref",
    "zpl_weight",
    "centroid_shift_mev",
    "closure",
    "integral",
    "seconds",
    "verdict",
)
