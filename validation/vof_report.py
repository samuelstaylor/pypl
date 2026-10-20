#!/usr/bin/env python3
"""Thorough, sectioned comparison of every route to the 1D CCD vibrational overlap function (VOF).

    python validation/vof_report.py                          # everything, ~8 min
    python validation/vof_report.py --skip-convergence       # without the grid-convergence scan
    python validation/vof_report.py --fits parabola spline morse --tag three-fits
    python validation/vof_report.py --temp 300 --tag 300K
    python validation/vof_report.py --wf-levels 0 1 2 3 --convergence-stencils 3 5 7 --tag stencils

Writes ``validation/outputs/vof_report/<tag or default>/``. Open ``index.md`` first: it lists the
headline numbers and every file with a one-line description. The sections are:

    01_potential_fits   each fit of the computed images: residuals, implied hbar*omega, and how far it
                        departs from the harmonic well where the high levels live
    02_fd_grid          the finite-difference grid drawn point by point, the wavefunctions it gives,
                        their decay to the walls, and convergence with grid points and stencil
    03_wavefunctions    vibrational wavefunctions: analytic Hermite vs finite differences on each fit
    04_fc_matrix        the Franck-Condon matrix of each method and its error against the reference
    05_vof_methods      VOF of every F_ij method on the harmonic wells (the like-for-like test)
    06_vof_fits         VOF of the finite-difference method on every fit, against the harmonic reference
    07_cubic_perturbation  the first-order cubic correction against the exact solution of the same cubic, or,
                        where the cubic turns over inside the grid, against the harmonic reference and fd on
                        the bounded polynomial fits (b -> 0 scaling on the scaled cubics that stay bounded)

The harmonic methods (``pypl``, ``hermite``) only know the parabola near each minimum, so in
sections 03 and 06 every non-parabola difference is the anharmonicity of that fit, not an error.
All VOF parameters (``--temp``, ``--sigma``, ``--fd-npts``, ...) are the shared flags of
``harness.methods`` and mean the same as in ``compare_vof.py`` and ``plot_ccd_wavefunctions.py``.
The report's own settings (levels drawn, convergence grids and stencils, cubic fits and b scales)
are the flags of the "report settings" group; every value used is printed and saved in settings.json.
"""

import argparse
import contextlib
import dataclasses
import io
import time

import matplotlib
import numpy as np
from matplotlib.patches import Patch

from harness import ccd_data, fits, plots
from harness.ccd_data import STATE_LABELS, STATE_NAMES, STATES
from harness.methods import COMPARISON_COLUMNS, CSV_KEYS, DIRECTIONS, METHODS, add_vof_args, compare_results
from harness.methods import compute_vof, direction_states, fd_q_range, params_from_args
from harness.output import OutputDir, markdown_table, text_table
from pypl.config_coord_1d_vof_solver import FD_STENCILS, KINETIC, config_coord_1d_vof_solver, grid_spacing

SECTIONS = {
    "01_potential_fits": "Each fit of the computed images: residuals, implied hbar*omega, and how far it departs from "
    "the harmonic well out where the high vibrational levels live.",
    "02_fd_grid": "The finite-difference grid point by point, the wavefunctions it gives, their decay to the walls, "
    "and convergence with grid points and stencil.",
    "03_wavefunctions": "Vibrational wavefunctions: analytic Hermite functions vs finite differences on each fit.",
    "04_fc_matrix": "The Franck-Condon matrix of each method and its error against the reference.",
    "05_vof_methods": "VOF of every F_ij method on the harmonic wells: the like-for-like test.",
    "06_vof_fits": "VOF of the finite-difference method on each fit, against the harmonic reference.",
    "07_cubic_perturbation": "The first-order cubic anharmonic correction (method cubic): the cubic fits, its "
    "validation against the exact solution of the same cubic (where one exists on the grid), wavefunctions, F_ij "
    "and VOF.",
}
# Defaults of the report settings (each one a command-line option, see parse_args).
SHOWN_LEVELS = (0, 1, 2, 5, 10, 20, 30, 40)  # levels drawn individually (clipped to the basis; the last is added)
CONVERGENCE_NPTS = (2001, 4001, 8001, 16001, 32001)
CONVERGENCE_STENCILS = (3, 5)
CUBIC_SPECS = ("cubic", "cubic-near", "cubic-free", "cubic-near-free")
CUBIC_SCALES = (1.0, 0.5, 0.25, 0.125, 0.0625)
MIN_EXACT_SCALES = 3  # b scales with an exact solution needed for the slopes; smaller ones are added by halving
SMALLEST_SCALE = 2.0**-12

# Short legend labels of the methods in the report's figures.
REPORT_LABELS = {
    "pypl": "PyPL (Ruhoff recurrences)",
    "hermite": "Hermite quadrature",
    "fd": "FD",
    "cubic": "First-order cubic PT",
}
WIDTH = 7.2  # in, double column
RES_FLOOR = 1e-6  # meV: rms residuals below this are a machine-precision interpolation
MINUS = "\N{MINUS SIGN}"


@dataclasses.dataclass
class Report:
    """Everything the sections share."""

    ccd: dict
    params: object
    out: OutputDir
    plt: object
    fits: list  # every fit spec compared in 01
    bounded: list  # the specs with bound states in both wells, used in 02-06
    reference: str
    tol_fc: float
    tol_vof: float
    shown_levels: tuple = SHOWN_LEVELS
    convergence_npts: tuple = CONVERGENCE_NPTS
    convergence_stencils: tuple = CONVERGENCE_STENCILS
    cubic_specs: tuple = CUBIC_SPECS
    cubic_scales: tuple = CUBIC_SCALES
    wells: dict = dataclasses.field(default_factory=dict)  # (state, spec) -> (q, energies, psi), from 03
    headline: list = dataclasses.field(default_factory=list)

    def order(self, state):
        return self.params.order_gs if state == "gs" else self.params.order_es

    def turning_point(self, state, v):
        return config_coord_1d_vof_solver.classical_turning_point(v, self.ccd["%s_phonon" % state])

    def fd_q_range(self):
        """CCD Q range of the finite-difference grid of the VOF solves (see ``harness.methods.fd_q_range``)."""

        return fd_q_range(self.ccd, self.params)

    def levels(self, n):
        return sorted({v for v in self.shown_levels if v < n} | {n - 1})


# ---------------------------------------------------------------------------------- figure helpers


def new_figure(plt, nrows, ncols, height, width=WIDTH, **kw):
    """A constrained-layout figure of ``nrows x ncols`` panels (``axes`` always 2D)."""

    return plt.subplots(nrows, ncols, figsize=(width, height), squeeze=False, constrained_layout=True, **kw)


def unique_handles(axes):
    """Legend handles and labels of every axis in ``axes``, each label once, in drawing order."""

    handles, labels = [], []
    for ax in np.ravel(axes):
        for h, lab in zip(*ax.get_legend_handles_labels()):
            if lab not in labels and not lab.startswith("_"):
                handles.append(h)
                labels.append(lab)
    return handles, labels


def shared_legend(fig, axes=None, ncol=3, handles=None, labels=None):
    """One legend for the figure, below the panels (``plots.legend_below``)."""

    if handles is None:
        handles, labels = unique_handles(axes)
    if handles:
        plots.legend_below(fig, handles, labels, ncol=ncol, y=0.0, frameon=False, columnspacing=1.4)


def symlog_if_wide(ax, values, threshold, linthresh):
    """Use a symmetric-log y axis when ``max |values|`` exceeds ``threshold`` (linear otherwise)."""

    values = np.asarray([v for v in np.ravel(values) if np.isfinite(v)])
    if values.size and np.max(np.abs(values)) > threshold:
        ax.set_yscale("symlog", linthresh=linthresh)
        tight_symlog_limits(ax, [values], linthresh)
        return True
    return False


def tight_symlog_limits(ax, values, linthresh):
    """y limits of a symlog axis that cover ``values`` without an empty half."""

    values = np.asarray([v for v in np.concatenate([np.ravel(v) for v in values]) if np.isfinite(v)])
    lo, hi = float(values.min()), float(values.max())
    bottom = 2.0 * lo if lo < -linthresh else (-linthresh if lo < -0.1 * linthresh else -0.1 * linthresh)
    top = 2.0 * hi if hi > linthresh else (linthresh if hi > 0.1 * linthresh else 0.1 * linthresh)
    ax.set_ylim(bottom, top)
    symlog_ticks(ax, linthresh)


def symlog_ticks(ax, linthresh):
    """At most ~4 decades labelled per sign on a symlog y axis, as 0 and +-10^k (no crowded or duplicate labels)."""

    lo, hi = ax.get_ylim()
    k0 = int(np.ceil(np.log10(linthresh)))
    k_hi = int(np.floor(np.log10(max(abs(lo), abs(hi), linthresh))))
    step = max(1, int(np.ceil((k_hi - k0 + 1) / 4)))
    ticks = [0.0]
    # +-linthresh itself sits right next to 0, so label it only when there is little else to label.
    for k in range(k_hi, k0 if k_hi > k0 + 1 else k0 - 1, -step):
        ticks += [t for t in (10.0**k, -(10.0**k)) if lo <= t <= hi]
    ax.yaxis.set_major_locator(matplotlib.ticker.FixedLocator(sorted(ticks)))
    ax.yaxis.set_major_formatter(
        matplotlib.ticker.FuncFormatter(
            lambda t, _: "0" if t == 0 else r"$%s10^{%d}$" % (MINUS if t < 0 else "", round(np.log10(abs(t))))
        )
    )
    ax.yaxis.set_minor_locator(matplotlib.ticker.NullLocator())


def sci_tex(x, digits=2):
    """``x`` as math text, e.g. ``-1.15\\times10^{-3}``."""

    if x == 0:
        return "0"
    k = int(np.floor(np.log10(abs(x))))
    return r"%.*f\times10^{%d}" % (digits, x / 10.0**k, k)


def direction_title(direction):
    """``Emission, ES -> GS`` with the dataset's math-text state labels (e.g. ``$^3E$ -> $^1A_1$``)."""

    initial, final = direction_states(direction)
    return r"%s, %s $\rightarrow$ %s" % (direction.capitalize(), STATE_LABELS[initial], STATE_LABELS[final])


def log_ticks_125(axis):
    """Ticks at 1, 2 and 5 x 10^k on a log axis spanning about a decade, labelled as math text."""

    def label(x, _):
        k = int(np.floor(np.log10(x) + 1e-9))
        m = x / 10.0**k
        return r"$10^{%d}$" % k if abs(m - 1) < 1e-6 else r"$%g\times10^{%d}$" % (round(m, 3), k)

    axis.set_major_locator(matplotlib.ticker.LogLocator(subs=(1.0, 2.0, 5.0), numticks=12))
    axis.set_major_formatter(matplotlib.ticker.FuncFormatter(label))
    axis.set_minor_locator(matplotlib.ticker.NullLocator())


def integer_maps(fig):
    """Integer level ticks on every image panel of ``fig`` (the shared ``plots.fc_error_figure``); returns ``fig``."""

    for ax in fig.axes:
        if ax.images:
            integer_ticks(ax, y=True)
    return fig


def integer_ticks(*axes, x=True, y=False):
    """Integer tick positions (vibrational levels)."""

    for ax in axes:
        if x:
            ax.xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(integer=True))
        if y:
            ax.yaxis.set_major_locator(matplotlib.ticker.MaxNLocator(integer=True))


def fraction_label(s):
    """``1``, ``1/2``, ``1/4`` ... for powers of two, otherwise ``%g``."""

    if s > 0 and s < 1 and abs(np.log2(1 / s) - round(np.log2(1 / s))) < 1e-9:
        return "1/%d" % round(1 / s)
    return "%g" % s


def method_label(method):
    return REPORT_LABELS.get(method, METHODS[method].label)


def fd_kwargs(params, npts=None, stencil=None):
    """Grid keywords of ``compute_franck_condon_integrals_fd`` from ``params``.

    ``npts`` / ``stencil`` override the parameters (convergence scan). The grid spacing, Q range and eigensolver shift
    are passed only when the harness's VOFParams has them (``fd_h``, ``fd_range``, ``fd_shift``) and they are set.
    """

    kw = {"stencil": stencil or params.fd_stencil, "margin": params.fd_margin}
    h = getattr(params, "fd_h", None)
    if npts is not None or not h:
        kw["npts"] = npts or params.fd_npts
    else:
        kw["h"] = h
    q_range = getattr(params, "fd_range", None)
    if q_range is not None:
        kw["q_range"] = tuple(q_range)
    shift = getattr(params, "fd_shift", None)
    if shift is not None:
        kw["shift"] = shift
    return kw


def quadrature_kwargs(params):
    """``npts`` and, when the harness has them, ``h`` / ``margin`` of the Hermite (Simpson) quadrature grid."""

    kw = {"npts": params.npts}
    if getattr(params, "hermite_h", None):
        kw["h"] = params.hermite_h
    if getattr(params, "hermite_margin", None) is not None:
        kw["margin"] = params.hermite_margin
    return kw


# ---------------------------------------------------------------------------------- 01


def fit_legend_label(ccd, spec, pots):
    """``fits.label(spec)``, noting the wells where the fit is unbounded."""

    unbounded = [STATE_LABELS[s] for s in STATES if not pots[s].bounded]
    if not unbounded:
        return fits.label(spec)
    return "%s (unbounded%s)" % (fits.label(spec), "" if len(unbounded) == len(STATES) else ": " + unbounded[0])


def section_potential_fits(rep):
    """Fit residuals, implied frequencies and the departure from harmonic beyond the sampled range."""

    name = "01_potential_fits"
    rep.out.section(name, SECTIONS[name])
    rep.out.log("\n[%s] %d fits of the computed images" % (name, len(rep.fits)))
    plt, ccd = rep.plt, rep.ccd
    pots = {spec: {s: fits.make_potential(ccd, s, spec) for s in STATES} for spec in rep.fits}
    labels = {spec: fit_legend_label(ccd, spec, pots[spec]) for spec in rep.fits}
    fig, axes = new_figure(plt, 4, 2, 9.4)
    rows = []
    q_lo, q_hi = rep.fd_q_range()
    q_grid = np.linspace(q_lo, q_hi, 2001)
    letters = iter("abcdefgh")
    for col, state in enumerate(STATES):
        q, e = ccd_data.samples(ccd, state)
        center = ccd_data.center(ccd, state)
        harmonic = fits.make_potential(ccd, state, "parabola")
        used = harmonic.index
        tp = rep.turning_point(state, rep.order(state) - 1)
        q_near = np.linspace(q.min() - 0.05, q.max() + 0.05, 600)
        q_far = np.linspace(center - tp - 0.3, center + tp + 0.3, 1200)
        others = np.setdiff1d(np.arange(len(q)), used)
        ax_fit, ax_full, ax_res, ax_far = axes[:, col]
        for ax in (ax_fit, ax_full):
            ax.plot(q[used], e[used], "o", color="black", ms=4, zorder=5, label="images used by the parabola")
            ax.plot(q[others], e[others], "o", mfc="white", mec="black", ms=4, mew=0.8, zorder=5, label="other images")
        for ax in (ax_full, ax_far):
            ax.axvspan(q.min(), q.max(), color=plots.OKABE_ITO["yellow"], alpha=0.3, lw=0, label="sampled range")
        residuals, departures = [], []
        for spec in rep.fits:
            pot = pots[spec][state]
            st = plots.fit_style(spec, rep.fits)
            kw = dict(color=st["color"], ls=st["ls"], lw=1.1, alpha=1.0 if pot.bounded else 0.55, label=labels[spec])
            ax_fit.plot(q_near, pot(q_near), **kw)
            ax_full.plot(q_grid, pot(q_grid), **kw)
            ax_res.plot(q, pot.residuals() * 1e3, marker="o", ms=2.5, **dict(kw, lw=0.8))
            dv = (pot(q_far) - harmonic(q_far)) * 1e3
            ax_far.plot(q_far, dv, **kw)
            residuals.append(pot.residuals() * 1e3)
            departures.append(dv)
            q_min, e_min = pot.minimum()
            rows.append(
                {
                    "state": STATE_NAMES[state],
                    "fit": spec,
                    "label": fits.label(spec),
                    "images": fits.images(ccd, state, spec),
                    "bounded": pot.bounded,
                    "fitted_images": "%d-%d of %d" % (pot.index[0] + 1, pot.index[-1] + 1, len(q)),
                    "rms_fitted_meV": pot.rms(pot.index) * 1e3,
                    "rms_all_images_meV": pot.rms() * 1e3,
                    "max_abs_residual_meV": np.max(np.abs(pot.residuals())) * 1e3,
                    "hbar_omega_at_min_meV": pot.phonon_energy(),
                    "q_min": q_min,
                    "e_min_eV": e_min,
                    "V_left_wall_eV": float(pot(q_lo)),
                    "V_right_wall_eV": float(pot(q_hi)),
                    "dV_left_turning_point_meV": (pot(center - tp) - harmonic(center - tp)) * 1e3,
                    "dV_right_turning_point_meV": (pot(center + tp) - harmonic(center + tp)) * 1e3,
                    "description": pot.describe(),
                }
            )
            rep.out.detail("  %s %-10s %s" % (STATE_NAMES[state], spec, pot.describe()))
            rep.out.detail(
                "  %s %-10s V at the FD walls: %.4g eV (Q = %.2f), %.4g eV (Q = %.2f)"
                % (STATE_NAMES[state], spec, pot(q_lo), q_lo, pot(q_hi), q_hi)
            )

        span = e.max() - e.min()
        ax_fit.set_ylim(e.min() - 0.12 * span, e.max() + 0.15 * span)
        ax_fit.set_xlim(q_near[0], q_near[-1])
        for wall in (q_lo, q_hi):
            ax_full.axvline(wall, color="black", lw=1.4, label="FD grid walls")
        top = 1.15 * max(harmonic(q_lo), harmonic(q_hi))
        ax_full.set_ylim(e.min() - 0.08 * (top - e.min()), top)
        ax_full.set_xlim(q_lo - 0.25, q_hi + 0.25)
        ax_res.axhline(0.0, color="0.5", lw=0.6)
        symlog_if_wide(ax_res, residuals, 20.0, 1.0)
        for side in (-1, 1):
            ax_far.axvline(center + side * tp, color="0.4", ls=":", lw=0.9, label="turning points of the highest level")
        ax_far.set_yscale("symlog", linthresh=1.0)
        tight_symlog_limits(ax_far, departures, 1.0)
        ax_far.yaxis.set_minor_locator(matplotlib.ticker.NullLocator())
        ax_far.set_xlim(q_far[0], q_far[-1])

        ax_fit.set_title("%s: fits near the images" % STATE_LABELS[state])
        ax_full.set_title("%s: full FD grid" % STATE_LABELS[state])
        ax_res.set_title("%s: residuals" % STATE_LABELS[state])
        ax_far.set_title("%s: departure from the parabola" % STATE_LABELS[state])
        ax_fit.set_ylabel(plots.ENERGY_LABEL)
        ax_full.set_ylabel(plots.ENERGY_LABEL)
        ax_res.set_ylabel("Fit %s computed (meV)" % MINUS)
        ax_far.set_ylabel(r"$V_\mathrm{fit} - V_\mathrm{parabola}$ (meV)")
        for a in axes[:, col]:
            a.set_xlabel(plots.Q_LABEL)
            plots.style(a)
    for a in axes.T.ravel():
        plots.panel_label(a, next(letters))
    shared_legend(fig, axes, ncol=3)
    rep.out.save_fig(
        fig,
        name + "/fits_overview.png",
        "Per well (columns): (a, e) the computed images (filled = used by the parabola) with every fit near the "
        "samples; (b, f) every fit over the full Q range of the finite-difference grid (%.2f..%.2f), with the sampled "
        "range shaded and the grid walls in black (V at the walls is in fit_summary.csv and log.txt); (c, g) the "
        "residual at each image (symmetric-log axis when a fit misses by more than 20 meV); (d, h) V_fit - "
        "V_parabola out to the classical turning points of the highest level (dotted), symmetric-log axis. "
        "Faded curves are unbounded fits; spline variants share one colour and differ by line style." % (q_lo, q_hi),
    )

    # rms residual per fit: horizontal bars with the fit names on the axis.
    n_fits = len(rep.fits)
    fig, axes = new_figure(plt, 1, 2, 0.24 * n_fits + 1.3, sharey=True)
    y = np.arange(n_fits)[::-1]
    for ax, state in zip(axes[0], STATES):
        for k, spec in enumerate(rep.fits):
            r = next(r for r in rows if r["fit"] == spec and r["state"] == STATE_NAMES[state])
            value, color = r["rms_all_images_meV"], plots.fit_style(spec, rep.fits)["color"]
            if value < RES_FLOOR:
                ax.plot(RES_FLOOR, y[k], marker="|", ms=8, mew=1.5, color=color)
                ax.text(
                    RES_FLOOR * 1.6, y[k], "exact interpolation", va="center", fontsize=plots.FIGURE["font_size"] - 2
                )
                continue
            ax.barh(
                y[k],
                value - RES_FLOOR / 3,
                left=RES_FLOOR / 3,
                height=0.7,
                color=color,
                alpha=1.0 if r["bounded"] else 0.45,
                hatch=None if r["bounded"] else "////",
                edgecolor="white" if r["bounded"] else color,
                lw=0.4,
            )
        ax.set_xscale("log")
        ax.set_xlim(RES_FLOOR / 3, None)
        ax.set_title("%s well (%d images)" % (STATE_LABELS[state], ccd["n_%s_images" % state]))
        ax.set_xlabel("rms residual over all images (meV)")
        plots.style(ax)
        ax.grid(axis="y", visible=False)
        ax.tick_params(axis="y", which="both", length=0)
    axes[0, 0].set_yticks(y)
    axes[0, 0].set_yticklabels([fits.label(s) for s in rep.fits])
    if any(not r["bounded"] for r in rows):
        shared_legend(
            fig,
            handles=[Patch(facecolor="0.6", alpha=0.45, hatch="////", edgecolor="0.6")],
            labels=["unbounded in that well (no bound states)"],
            ncol=1,
        )
    rep.out.save_fig(
        fig,
        name + "/fit_residual_rms.png",
        "rms residual over all images, per fit and well (log scale). Hatched = unbounded, so unusable for the "
        "finite-difference solve; fits with an rms below %g meV interpolate the images to machine precision and "
        "are marked 'exact interpolation'." % RES_FLOOR,
    )

    for state in STATES:
        q, e = ccd_data.samples(ccd, state)
        columns = {"Q": q, "E_eV": e}
        columns.update({"residual_%s_meV" % s: pots[s][state].residuals() * 1e3 for s in rep.fits})
        rep.out.save_csv(
            name + "/residuals_%s.csv" % state,
            columns,
            "%s: fit minus computed energy at every image, per fit." % STATE_NAMES[state],
        )
    rep.out.save_table(
        name + "/fit_summary.csv",
        rows,
        "Per well and fit: images fitted, rms / max residual, hbar*omega at the minimum, the minimum, V at the two "
        "walls of the finite-difference grid, the departure from harmonic at the outer turning points, and the full "
        "description.",
    )

    columns = [
        ("state", "well", "%s"),
        ("label", "fit", "%s"),
        ("fitted_images", "images", "%s"),
        ("rms_all_images_meV", "rms (meV)", "%.3f"),
        ("max_abs_residual_meV", "max |res| (meV)", "%.3f"),
        ("hbar_omega_at_min_meV", "hbar w at min (meV)", "%.3f"),
        ("dV_right_turning_point_meV", "dV at right turning pt (meV)", "%.3g"),
        ("V_left_wall_eV", "V left wall (eV)", "%.4g"),
        ("V_right_wall_eV", "V right wall (eV)", "%.4g"),
        ("bounded", "bounded", "%s"),
    ]
    for line in text_table(rows, columns):
        rep.out.log("  " + line)
    rep.headline += ["## Potential fits", ""] + markdown_table(rows, columns) + [""]


# ---------------------------------------------------------------------------------- 02


def fd_solve(rep, direction, npts=None, stencil=None):
    """The default fd solve (parabola wells) of one direction, optionally with another grid."""

    initial, final = direction_states(direction)
    ni, nf = rep.params.orders(direction)
    solver = config_coord_1d_vof_solver(
        rep.ccd["%s_phonon" % initial], rep.ccd["%s_phonon" % final], rep.ccd["delta_q"]
    )
    start = time.perf_counter()
    solver.compute_franck_condon_integrals_fd(
        ni,
        nf,
        fits.make_potential(rep.ccd, initial),
        fits.make_potential(rep.ccd, final),
        q_i=ccd_data.center(rep.ccd, initial),
        q_f=ccd_data.center(rep.ccd, final),
        **fd_kwargs(rep.params, npts, stencil),
    )
    return solver, time.perf_counter() - start


def draw_grid_direction(rep, direction):
    """Grid, wavefunctions, a zoom onto individual grid points, and the decay to the walls."""

    plt, ccd = rep.plt, rep.ccd
    solver, seconds = fd_solve(rep, direction)
    info, q = solver.fd_info, solver.fd_q
    fig, axes = new_figure(plt, 2, 3, 5.6, gridspec_kw={"width_ratios": [1.7, 1, 1]})
    stats = []
    tail_colors = (plots.OKABE_ITO["blue"], plots.OKABE_ITO["orange"], plots.OKABE_ITO["vermillion"])
    for row, (role, state, energies, psi) in enumerate(
        (
            ("initial", direction_states(direction)[0], solver.fd_energies_i, solver.fd_wavefunctions_i),
            ("final", direction_states(direction)[1], solver.fd_energies_f, solver.fd_wavefunctions_f),
        )
    ):
        pot = fits.make_potential(ccd, state)
        v_grid = pot(q)
        n = len(energies)
        spacing = ccd["%s_phonon" % state] * 1e-3
        color = plots.STATE_COLORS[state]
        ax_all, ax_zoom, ax_tail = axes[row]

        # (a) the shown levels at their eigenvalues, each scaled to the gap to the next shown level
        ax_all.plot(q, v_grid, color=color, lw=1.3, label="parabola, %s" % STATE_LABELS[state])
        shown = rep.levels(n)
        y_lo, y_hi = v_grid.min() - 1.0 * spacing, energies[-1] + 1.5 * spacing
        last_label = -np.inf
        label_x = ccd_data.center(ccd, state) + rep.turning_point(state, n - 1) + 0.4
        for k, v in enumerate(shown):
            gap = energies[shown[k + 1]] - energies[v] if k + 1 < len(shown) else energies[v] - energies[shown[k - 1]]
            amp = 0.45 * gap * psi[v] / np.max(np.abs(psi[v]))
            # Drawn only where the level lives (its turning points plus its decaying tails), clear of the labels.
            near = np.abs(q - ccd_data.center(ccd, state)) < rep.turning_point(state, v) + 0.3
            ax_all.fill_between(
                q[near],
                energies[v],
                energies[v] + amp[near],
                color="0.45",
                alpha=0.5,
                lw=0,
                label=r"$\psi_v$ drawn at $E_v$",
            )
            # Label only levels far enough apart to stay legible (the lowest ones crowd at this scale).
            if energies[v] - last_label > 0.06 * (y_hi - y_lo) or v == shown[-1]:
                ax_all.text(
                    label_x,
                    energies[v],
                    "$v = %d$" % v,
                    ha="left",
                    va="center",
                    fontsize=plots.FIGURE["font_size"] - 2,
                    # The padding colour (yellow at alpha 0.3 on white), so grid lines do not cross the label.
                    bbox={"facecolor": "#FBF7C6", "edgecolor": "none", "pad": 0.6},
                )
                last_label = energies[v]
        tp = rep.turning_point(state, n - 1)
        center = ccd_data.center(ccd, state)
        for side in (-1, 1):
            edge = center + side * tp
            outer = q.min() if side < 0 else q.max()
            ax_all.axvspan(
                min(edge, outer),
                max(edge, outer),
                color=plots.OKABE_ITO["yellow"],
                alpha=0.3,
                lw=0,
                label="padding beyond the turning points",
            )
        for a in (ax_all, ax_tail):
            a.axvline(q.min(), color="black", lw=1.4, label="grid walls")
            a.axvline(q.max(), color="black", lw=1.4)
        ax_all.set_ylim(y_lo, y_hi)
        ax_all.set_xlim(q.min() - 0.02 * np.ptp(q), q.max() + 0.02 * np.ptp(q))
        ax_all.set_title("%s state, %s" % (role.capitalize(), STATE_LABELS[state]))
        ax_all.set_ylabel(plots.ENERGY_LABEL)

        # (b) the highest level at the well centre, where its wavelength is shortest, point by point
        k_max = np.sqrt((energies[-1] - v_grid.min()) / KINETIC)
        wavelength = 2 * np.pi / k_max
        ppw = wavelength / info["h"]
        nodes = np.nonzero(np.diff(np.sign(psi[-1])))[0]
        node = nodes[np.argmin(np.abs(q[nodes] - center))] if len(nodes) else int(np.argmin(np.abs(q - center)))
        window = np.abs(q - q[node]) < wavelength / 8
        ax_zoom.plot(q[window], psi[-1][window], "-", color=color, lw=0.6, alpha=0.6)
        ax_zoom.plot(q[window], psi[-1][window], "o", color=color, ms=1.6, label="grid points")
        ax_zoom.axhline(0.0, color="0.5", lw=0.6)
        ax_zoom.set_title(r"$\psi_{%d}$ at the well centre" % (n - 1))
        ax_zoom.set_ylabel(r"$\psi_{%d}$" % (n - 1))
        ax_zoom.xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(3))

        # (c) decay to the walls on a log scale
        for (v, tag), c in zip(((0, "0"), (n // 2, "n/2"), (n - 1, "n-1")), tail_colors):
            ax_tail.semilogy(
                q, np.maximum(np.abs(psi[v]) / np.max(np.abs(psi[v])), 1e-20), lw=0.8, color=c, label="$v = %s$" % tag
            )
        ax_tail.set_ylim(1e-20, 3)
        ax_tail.yaxis.set_major_locator(matplotlib.ticker.LogLocator(numticks=6))
        ax_tail.yaxis.set_minor_locator(matplotlib.ticker.NullLocator())
        ax_tail.set_title("Decay to the walls")
        ax_tail.set_ylabel(r"$|\psi_v|\,/\,\max|\psi_v|$")
        for a in axes[row]:
            a.set_xlabel(plots.Q_LABEL)
            plots.style(a)
        stats.append(
            {
                "direction": direction,
                "state": STATE_NAMES[state],
                "role": role,
                "levels": n,
                "npts": info["npts"],
                "h": info["h"],
                "stencil": info["stencil"],
                "frame": info["frame"],
                "q_lo": q.min(),
                "q_hi": q.max(),
                "points_per_shortest_wavelength": ppw,
                "zoom_points": int(window.sum()),
                "wall_amplitude": info["wall"],
                "orthonormality": info["orthonormality"],
                "ladder_error_meV": info["ladder_i" if role == "initial" else "ladder_f"],
                "seconds": seconds,
            }
        )
    for a, letter in zip(axes.ravel(), "abcdef"):
        plots.panel_label(a, letter)
    shared_legend(fig, axes, ncol=4)
    fig.suptitle("FD grid, %s" % direction_title(direction))
    s_i, s_f = stats
    rep.out.save_fig(
        fig,
        "02_fd_grid/fd_grid_%s.png" % direction,
        "%s: (a, d) both wells on the shared grid (%d points, h = %.2e, %d-point stencil, frame %s) with selected "
        "levels drawn at their eigenvalues, the padding beyond the turning points of the highest level shaded, the "
        "walls in black; (b, e) a quarter of the shortest wavelength of the highest level at the well centre, point by "
        "point (%.0f / %.0f points per wavelength); (c, f) |psi| / peak decaying to the walls (%.1e at the walls; "
        "the ~1e-16 floor is the float64 precision of the eigensolver). F_ij took %.1f s."
        % (
            direction.capitalize(),
            info["npts"],
            info["h"],
            info["stencil"],
            info["frame"],
            s_i["points_per_shortest_wavelength"],
            s_f["points_per_shortest_wavelength"],
            info["wall"],
            seconds,
        ),
    )
    return stats


def section_fd_grid(rep, convergence=True):
    name = "02_fd_grid"
    rep.out.section(name, SECTIONS[name])
    rep.out.log("\n[%s] grid visualization%s" % (name, " and convergence scan" if convergence else ""))
    stats = []
    for direction in DIRECTIONS:
        stats += draw_grid_direction(rep, direction)
    rep.out.save_table(
        name + "/grid_summary.csv",
        stats,
        "Per direction and well: grid points, spacing, range, "
        "resolution of the highest level, wall amplitude, orthonormality and level-ladder error.",
    )
    for line in text_table(
        stats,
        [
            ("direction", "direction", "%s"),
            ("state", "well", "%s"),
            ("npts", "points", "%d"),
            ("h", "h", "%.3e"),
            ("points_per_shortest_wavelength", "pts/wavelength", "%.0f"),
            ("wall_amplitude", "|psi| at walls", "%.1e"),
            ("ladder_error_meV", "ladder err (meV)", "%.1e"),
            ("seconds", "F_ij time (s)", "%.1f"),
        ],
    ):
        rep.out.log("  " + line)
    if not convergence:
        return

    # Convergence of F_ij against the Hermite quadrature (accurate to ~1e-12), emission only.
    ref = config_coord_1d_vof_solver(rep.ccd["es_phonon"], rep.ccd["gs_phonon"], rep.ccd["delta_q"])
    ref.compute_franck_condon_integrals_hermite(*rep.params.orders("emission"), **quadrature_kwargs(rep.params))
    rep.out.log(
        "  convergence scan (emission, parabola wells): stencils %s, grid points %s; reference Hermite quadrature, "
        "%d points" % (list(rep.convergence_stencils), list(rep.convergence_npts), len(ref.hermite_grid))
    )
    rows = []
    for stencil in rep.convergence_stencils:
        for npts in rep.convergence_npts:
            solver, seconds = fd_solve(rep, "emission", npts=npts, stencil=stencil)
            rows.append(
                {
                    "stencil": stencil,
                    "npts": npts,
                    "h": solver.fd_info["h"],
                    "q_lo": float(solver.fd_q.min()),
                    "q_hi": float(solver.fd_q.max()),
                    "max_dF_vs_hermite": float(np.max(np.abs(solver.fc_ints - ref.fc_ints))),
                    "ladder_error_meV": max(solver.fd_info["ladder_i"], solver.fd_info["ladder_f"]),
                    "seconds": seconds,
                }
            )
            rep.out.detail(
                "  convergence: stencil %d, %6d points, max dF %.2e, ladder %.2e meV, %.1f s"
                % (stencil, npts, rows[-1]["max_dF_vs_hermite"], rows[-1]["ladder_error_meV"], seconds)
            )
    keys = ("max_dF_vs_hermite", "ladder_error_meV", "seconds")
    for stencil in rep.convergence_stencils:
        sub = [r for r in rows if r["stencil"] == stencil]
        h = np.array([r["h"] for r in sub])
        for key in keys:
            values = np.array([r[key] for r in sub])
            # The measured order is the log-log slope; skip points already at the float64 floor.
            ok = values > (1e-12 if key != "seconds" else 0)
            slope = float(np.polyfit(np.log(h[ok]), np.log(values[ok]), 1)[0]) if ok.sum() >= 2 else float("nan")
            for r in sub:
                r["slope_" + key] = slope
        rep.out.log(
            "  %d-point stencil: F_ij error slope %.2f, ladder error slope %.2f (expected %d)"
            % (stencil, sub[0]["slope_max_dF_vs_hermite"], sub[0]["slope_ladder_error_meV"], stencil - 1)
        )

    fig, axes = new_figure(plt=rep.plt, nrows=1, ncols=2, height=3.2)
    axes = axes[0]
    markers = ("s", "o", "^", "D", "v")
    for n_st, stencil in enumerate(rep.convergence_stencils):
        sub = [r for r in rows if r["stencil"] == stencil]
        h = np.array([r["h"] for r in sub])
        color = plots.OKABE_ITO[("blue", "vermillion", "green", "orange", "purple")[n_st % 5]]
        for ax, key in zip(axes, keys[:2]):
            ax.loglog(
                h,
                [r[key] for r in sub],
                marker=markers[n_st % len(markers)],
                color=color,
                ms=4,
                label="%d-point stencil" % stencil,
            )
        # Order guides, drawn 5x above the data so they do not hide it.
        order = stencil - 1
        anchor = 5 * sub[0]["max_dF_vs_hermite"]
        axes[0].loglog(
            h, anchor * (h / h[0]) ** order, color=color, ls=":", lw=0.9, alpha=0.8, label=r"$\propto h^{%d}$" % order
        )
    default = next(s for s in stats if s["direction"] == "emission")
    default_h = default["h"]
    span = rows[0]["h"] * (rows[0]["npts"] - 1)
    for ax in axes:
        if default_h:
            ax.axvline(default_h, color="black", ls="--", lw=0.8, label="default grid (%d points)" % default["npts"])
        ax.invert_xaxis()
        log_ticks_125(ax.xaxis)
        top = ax.secondary_xaxis(
            "top", functions=(lambda x: span / np.maximum(x, 1e-300), lambda n: span / np.maximum(n, 1e-300))
        )
        log_ticks_125(top.xaxis)
        top.set_xlabel("Grid points")
        top.tick_params(direction="in", which="both")
        ax.set_xlabel(r"Grid spacing $h$ (%s)" % plots.Q_UNIT)
        plots.style(ax)
        ax.xaxis.set_ticks_position("bottom")
    axes[0].set_ylabel(r"$\max\,|F_{ij}^\mathrm{FD} - F_{ij}^\mathrm{Hermite}|$")
    axes[1].set_ylabel(r"$\max_v\,|E_v - E_0 - v\hbar\omega|$ (meV)")
    for a, letter in zip(axes, "ab"):
        plots.panel_label(a, letter, y=1.17)
    shared_legend(fig, axes, ncol=4)
    rep.out.save_fig(
        fig,
        name + "/fd_convergence.png",
        "Emission, parabola wells: (a) max |F_ij(FD) - F_ij(Hermite quadrature, %d points)| and (b) the level-ladder "
        "error vs grid spacing, for the %s-point stencils, with O(h^n) guides (dotted). Measured log-log slopes "
        "(F_ij): %s. The F_ij timings are in fd_convergence.csv."
        % (
            len(ref.hermite_grid),
            ", ".join(str(s) for s in rep.convergence_stencils),
            ", ".join(
                "%d-point %.2f" % (s, next(r["slope_max_dF_vs_hermite"] for r in rows if r["stencil"] == s))
                for s in rep.convergence_stencils
            ),
        ),
    )
    rep.out.save_table(
        name + "/fd_convergence.csv", rows, "The convergence scan behind fd_convergence.png, with timings and slopes."
    )


# ---------------------------------------------------------------------------------- 03


def solve_well(rep, state, spec):
    """FD eigenstates of one well on its own physical grid: ``(q, energies in eV, psi)``."""

    n = rep.order(state)
    q = well_grid(rep, state)
    pot = fits.make_potential(rep.ccd, state, spec)
    shift = getattr(rep.params, "fd_shift", None)
    extra = {} if shift is None else {"shift": shift}
    energies, psi = config_coord_1d_vof_solver.solve_schrodinger_fd(pot(q), q, n, rep.params.fd_stencil, **extra)
    return q, energies, psi


def wavefunction_rows_figure(rep, state, shown, window_of, left, right, titles, suptitle):
    """One row per level: curves on the left, differences on the right, each row scaled by its own power of ten.

    ``left(v)`` / ``right(v)`` return ``[(x, y, style dict), ...]`` for level ``v``.
    """

    plt = rep.plt
    fig, axes = new_figure(plt, len(shown), 2, 0.95 * len(shown) + 1.0)
    for k, v in enumerate(shown):
        window = window_of(v)
        for x, y, kw in left(v):
            axes[k, 0].plot(x[window], y[window], **kw)
        curves = right(v)
        _, exp = plots.scaled(np.concatenate([y[window] for _, y, _ in curves]))
        for x, y, kw in curves:
            axes[k, 1].plot(x[window], y[window] / 10.0**exp, **kw)
        axes[k, 1].axhline(0.0, color="0.6", lw=0.5)
        axes[k, 0].set_ylabel("$v = %d$" % v)
        axes[k, 1].set_ylabel(r"$\times 10^{%d}$" % exp if exp else "")
        for a in axes[k]:
            a.yaxis.set_major_locator(matplotlib.ticker.MaxNLocator(3))
            a.xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(5))
            plots.style(a)
    axes[0, 0].set_title(titles[0])
    axes[0, 1].set_title(titles[1])
    for a in axes[-1]:
        a.set_xlabel(plots.Q_LABEL)
    shared_legend(fig, axes, ncol=3)
    fig.suptitle(suptitle)
    return fig


def section_wavefunctions(rep):
    name = "03_wavefunctions"
    rep.out.section(name, SECTIONS[name])
    rep.out.log("\n[%s] Hermite vs finite differences on %d fits" % (name, len(rep.bounded)))
    plt, ccd = rep.plt, rep.ccd
    rows = []
    for state in STATES:
        phonon = ccd["%s_phonon" % state]
        center = ccd_data.center(ccd, state)
        for spec in rep.bounded:
            rep.wells[(state, spec)] = solve_well(rep, state, spec)
        q = rep.wells[(state, rep.bounded[0])][0]
        h = grid_spacing(q)
        n = rep.order(state)
        chi = config_coord_1d_vof_solver.harmonic_wavefunctions(q, n - 1, phonon, center)
        for spec in rep.bounded:
            _, energies, psi = rep.wells[(state, spec)]
            for v in range(n):
                rows.append(
                    {
                        "state": STATE_NAMES[state],
                        "fit": spec,
                        "v": v,
                        "E_v_minus_E_0_meV": (energies[v] - energies[0]) * 1e3,
                        "ladder_dev_meV": (energies[v] - energies[0]) * 1e3 - v * phonon,
                        "max_abs_dpsi": float(np.max(np.abs(psi[v] - chi[v]))),
                        "max_abs_dpsi_rel": float(np.max(np.abs(psi[v] - chi[v])) / np.max(np.abs(chi[v]))),
                        "infidelity": float(1.0 - abs(np.sum(psi[v] * chi[v]) * h)),
                    }
                )

        def styled(spec, lw):
            st = plots.fit_style(spec, rep.bounded)
            return dict(color=st["color"], ls=st["ls"], lw=lw, label="FD, " + fits.label(spec))

        fig = wavefunction_rows_figure(
            rep,
            state,
            rep.levels(n),
            lambda v: np.abs(q - center) < rep.turning_point(state, v) + 0.45,
            lambda v: [(q, chi[v], dict(color="0.7", lw=2.6, label=r"Hermite (analytic), $\chi_v$"))]
            + [(q, rep.wells[(state, s)][2][v], styled(s, 0.9)) for s in rep.bounded],
            lambda v: [(q, rep.wells[(state, s)][2][v] - chi[v], styled(s, 0.9)) for s in rep.bounded],
            (r"$\chi_v$ and $\psi_v$", r"$\psi_v - \chi_v$"),
            "%s well: vibrational wavefunctions" % STATE_LABELS[state],
        )
        rep.out.save_fig(
            fig,
            name + "/wavefunctions_%s.png" % state,
            "%s: selected levels v (rows) from the analytic Hermite function (grey) and from finite differences on "
            "every bounded fit (left), and their difference psi_v - chi_v (right, each row scaled by the power of "
            "ten on its axis). On the parabola the difference is numerical (~1e-9, invisible next to the anharmonic "
            "fits); every per-level maximum is in wavefunction_errors.csv." % STATE_NAMES[state],
        )

    fig, axes = new_figure(plt, 2, 2, 5.0, sharex="col")
    for col, state in enumerate(STATES):
        for spec in rep.bounded:
            sub = [r for r in rows if r["fit"] == spec and r["state"] == STATE_NAMES[state]]
            v = [r["v"] for r in sub]
            st = plots.fit_style(spec, rep.bounded)
            kw = dict(color=st["color"], ls=st["ls"], lw=1.0, marker=".", ms=3, label="FD, " + fits.label(spec))
            axes[0, col].semilogy(v, [max(r["infidelity"], 1e-17) for r in sub], **kw)
            axes[1, col].semilogy(v, [max(r["max_abs_dpsi_rel"], 1e-17) for r in sub], **kw)
        axes[0, col].set_title("%s well" % STATE_LABELS[state])
        axes[1, col].set_xlabel("Vibrational level $v$")
        for a in axes[:, col]:
            a.set_ylim(bottom=1e-17)
            integer_ticks(a)
            a.yaxis.set_major_locator(matplotlib.ticker.LogLocator(numticks=6))
            a.yaxis.set_minor_locator(matplotlib.ticker.NullLocator())
            plots.style(a)
    axes[0, 0].set_ylabel(r"$1 - |\langle\chi_v|\psi_v\rangle|$")
    axes[1, 0].set_ylabel(r"$\max|\psi_v - \chi_v|\,/\,\max|\chi_v|$")
    for a, letter in zip(axes.ravel(), "abcd"):
        plots.panel_label(a, letter)
    shared_legend(fig, axes, ncol=2)
    fig.suptitle("Wavefunction differences from the Hermite functions")
    rep.out.save_fig(
        fig,
        name + "/wavefunction_errors.png",
        "Per level, well and fit: (a, b) the infidelity 1 - |<chi_v|psi_v>| and (c, d) max |psi_v - chi_v| relative to "
        "the analytic Hermite function (floored at 1e-17). On the parabola this is the numerical error of the "
        "finite-difference solve; on the other fits it is their anharmonicity.",
    )
    rep.out.save_table(
        name + "/wavefunction_errors.csv",
        rows,
        "Per well, fit and level: E_v - E_0, its departure "
        "from v hbar w, max |psi - chi| (absolute and / max |chi|) and the infidelity.",
    )

    summary = []
    for state in STATES:
        for spec in rep.bounded:
            sub = [r for r in rows if r["fit"] == spec and r["state"] == STATE_NAMES[state]]
            summary.append(
                {
                    "state": STATE_NAMES[state],
                    "fit": fits.label(spec),
                    "infid": max(r["infidelity"] for r in sub),
                    "dpsi": max(r["max_abs_dpsi_rel"] for r in sub),
                    "ladder": max(abs(r["ladder_dev_meV"]) for r in sub),
                }
            )
    columns = [
        ("state", "well", "%s"),
        ("fit", "fd on fit", "%s"),
        ("infid", "max infidelity", "%.2e"),
        ("dpsi", "max dpsi / max chi", "%.2e"),
        ("ladder", "max ladder dev (meV)", "%.3g"),
    ]
    for line in text_table(summary, columns):
        rep.out.log("  " + line)
    rep.headline += ["## Wavefunctions vs analytic Hermite (all levels)", ""] + markdown_table(summary, columns) + [""]


# ------------------------------------------------------------------------------- 04-05


def section_methods(rep):
    """VOF and F_ij of every method on the harmonic wells (sections 04 and 05)."""

    rep.out.section("04_fc_matrix", SECTIONS["04_fc_matrix"])
    rep.out.section("05_vof_methods", SECTIONS["05_vof_methods"])
    rep.out.log("\n[04_fc_matrix, 05_vof_methods] every method on the parabola wells, reference %s" % rep.reference)
    params = dataclasses.replace(rep.params, fd_fit=fits.DEFAULT_FIT)
    tests = [m for m in METHODS if m != rep.reference and not METHODS[m].anharmonic]
    comparison = {}
    for direction in DIRECTIONS:
        ref = compute_vof(rep.ccd, direction, rep.reference, params, log=rep.out.detail)
        ref.label = method_label(ref.method)
        tested = []
        for m in tests:
            tst = compute_vof(rep.ccd, direction, m, params, log=rep.out.detail)
            tst.label = method_label(m)
            tested.append((tst, compare_results(ref, tst, rep.tol_fc, rep.tol_vof)))
        comparison[direction] = (ref, tested)
        columns = {"E_minus_EZPL_meV": ref.axis, "VOF %s" % ref.label: ref.lineshape}
        columns.update({"VOF %s" % t.label: t.lineshape for t, _ in tested})
        columns.update({"%s - %s" % (t.label, ref.label): t.lineshape - ref.lineshape for t, _ in tested})
        rep.out.save_csv(
            "05_vof_methods/vof_%s.csv" % direction,
            columns,
            "%s: every method's lineshape and its " "difference from the reference." % direction,
        )
        for label, result in [(ref.method, ref)] + [(t.method, t) for t, _ in tested]:
            np_cols = {"i": np.arange(result.solver.ni)}
            np_cols.update({"j=%d" % j: result.solver.fc_ints[:, j] for j in range(result.solver.nf)})
            rep.out.save_csv(
                "04_fc_matrix/fc_%s_%s.csv" % (label, direction),
                np_cols,
                "F_ij from %s, %s (rows i = " "initial level, columns j = final level)." % (result.label, direction),
            )

    plt = rep.plt
    fig, axes = new_figure(plt, 1, 2, 3.2)
    vmin, vmax = -12.0, 0.0
    top = max(float(np.log10(np.max(np.abs(r.solver.fc_ints)))) for r, _ in comparison.values())
    for ax, (direction, (ref, _)) in zip(axes[0], comparison.items()):
        im = ax.imshow(
            np.log10(np.maximum(np.abs(ref.solver.fc_ints), 1e-18)),
            origin="lower",
            aspect="auto",
            cmap="viridis",
            vmin=vmin,
            vmax=vmax,
        )
        ax.set_title(direction_title(direction))
        ax.set_xlabel(r"$j$ (final level)")
        ax.set_ylabel(r"$i$ (initial level)")
        plots.style(ax, grid=False)
        integer_ticks(ax, y=True)
    fig.colorbar(im, ax=axes[0].tolist(), label=r"$\log_{10}|F_{ij}|$", extend="both" if top > vmax else "min")
    fig.suptitle("Franck%sCondon matrix, %s" % ("\N{EN DASH}", ref.label))
    rep.out.save_fig(
        fig,
        "04_fc_matrix/fc_matrix.png",
        "The reference Franck-Condon matrix, log10 |F_ij|, per direction (clipped at 1e-12%s)."
        % (
            "; values above 1, impossible for an overlap, are the float64 instability of the reference"
            if top > 0
            else ""
        ),
    )
    rep.out.save_fig(
        integer_maps(plots.fc_error_figure(plt, comparison, r"$F_{ij}$ error vs %s, parabola wells" % ref.label)),
        "04_fc_matrix/fc_error_maps.png",
        "log10 |F_ij(method) - F_ij(reference)| per method and direction; the cross marks the maximum.",
    )

    fig, axes = new_figure(plt, 1, 2, 3.0, sharey=True)
    for ax, (direction, (ref, tested)) in zip(axes[0], comparison.items()):
        for tst, _ in tested:
            dev = np.max(np.abs(tst.solver.fc_ints - ref.solver.fc_ints), axis=1)
            ax.semilogy(
                np.arange(len(dev)),
                np.maximum(dev, 1e-18),
                marker=".",
                ms=3,
                lw=1.0,
                color=plots.result_color(tst),
                ls=plots.METHOD_STYLES.get(tst.method, "-"),
                label=tst.label,
            )
        ax.set_title(direction_title(direction))
        ax.set_xlabel(r"$i$ (initial level)")
        integer_ticks(ax)
        plots.style(ax)
    axes[0, 0].set_ylabel(r"$\max_j\,|\Delta F_{ij}|$")
    for a, letter in zip(axes[0], "ab"):
        plots.panel_label(a, letter)
    shared_legend(fig, axes, ncol=3)
    fig.suptitle(r"$F_{ij}$ error by initial level, vs %s" % ref.label)
    rep.out.save_fig(
        fig,
        "04_fc_matrix/fc_error_by_level.png",
        "Largest F_ij difference from the reference in each row i, per method. With the PyPL reference the "
        "difference grows with the level: the float64 drift of its recurrence.",
    )

    title = "VOF on the harmonic wells, reference: %s" % comparison["emission"][0].label
    rep.out.save_fig(
        plots.vof_comparison_figure(plt, comparison, params, rep.tol_vof, title),
        "05_vof_methods/vof_methods.png",
        "Per direction: every method's VOF overlaid on the reference, and the signed and absolute error relative to "
        "the reference peak, with the tolerance (dashed). Every method's settings are in log.txt.",
    )
    rows = [stats for _, tested in comparison.values() for _, stats in tested]
    rep.out.save_table(
        "04_fc_matrix/fc_errors.csv",
        [{k: r[k] for k in CSV_KEYS} for r in rows],
        "Every F_ij and VOF error statistic per method and direction.",
    )
    rep.out.save_table(
        "05_vof_methods/summary.csv",
        [{k: r[k] for k in CSV_KEYS} for r in rows],
        "Every error statistic per method and direction, with OK / FAIL.",
    )
    for line in text_table(rows, COMPARISON_COLUMNS):
        rep.out.log("  " + line)
    rep.headline += ["## F_ij methods on the harmonic wells (reference %s)" % ref.label, ""]
    rep.headline += markdown_table(rows, COMPARISON_COLUMNS) + [""]
    return all(r["ok"] for r in rows)


# ---------------------------------------------------------------------------------- 06


def fits_vof_figure(rep, comparison, title):
    """Like ``plots.vof_comparison_figure``, with each fit drawn in its ``plots.fit_style`` (colour and line style)."""

    fig = rep.plt.figure(figsize=(WIDTH, 6.6), constrained_layout=True)
    grid = fig.add_gridspec(3, 2, height_ratios=(2.2, 1, 1))
    axes = []
    for col, (ref, tested) in enumerate(comparison.values()):
        top = fig.add_subplot(grid[0, col])
        mid = fig.add_subplot(grid[1, col], sharex=top)
        bottom = fig.add_subplot(grid[2, col], sharex=top)
        axes += [top, mid, bottom]
        peak = np.max(ref.lineshape)
        color = plots.DIRECTION_COLORS[ref.direction]
        top.fill_between(ref.axis, 0.0, ref.lineshape, color=color, alpha=0.18, linewidth=0)
        top.plot(ref.axis, ref.lineshape, color=color, lw=2.4, alpha=0.55, label="reference: %s" % ref.label)
        errors = [(tst.lineshape - ref.lineshape) / peak for tst, _ in tested]
        _, k = plots.scaled(np.concatenate(errors))
        for (tst, _), err in zip(tested, errors):
            st = plots.fit_style(tst.fit, rep.bounded)
            top.plot(tst.axis, tst.lineshape, color=st["color"], ls=st["ls"], lw=1.1, label=tst.label)
            mid.plot(ref.axis, err / 10.0**k, color=st["color"], ls=st["ls"], lw=0.9)
            bottom.semilogy(ref.axis, np.maximum(np.abs(err), 1e-18), color=st["color"], ls=st["ls"], lw=0.9)
        mid.axhline(0.0, color="0.6", lw=0.6)
        sideband = max([ref.sideband_max, 1e-30] + [t.sideband_max for t, _ in tested])
        top.set_ylim(0.0, 1.3 * sideband)
        if peak > top.get_ylim()[1]:
            plots.off_scale_note(top, peak)
        top.set_title(direction_title(ref.direction))
        bottom.set_ylim(1e-12, None)
        bottom.yaxis.set_major_locator(matplotlib.ticker.LogLocator(numticks=5))
        bottom.yaxis.set_minor_locator(matplotlib.ticker.NullLocator())
        if col == 0:
            top.set_ylabel(plots.VOF_LABEL)
            mid.set_ylabel(r"$\Delta$VOF / peak" + plots.times_ten(k))
            bottom.set_ylabel(r"$|\Delta$VOF$|$ / peak")
        bottom.set_xlabel(plots.E_ZPL_LABEL)
        for a in (top, mid, bottom):
            plots.style(a)
        for a in (top, mid):
            a.tick_params(labelbottom=False)
    shared_legend(fig, axes, ncol=2)
    fig.suptitle(title)
    return fig


def section_vof_fits(rep):
    name = "06_vof_fits"
    rep.out.section(name, SECTIONS[name])
    rep.out.log("\n[%s] fd on every bounded fit vs %s (harmonic)" % (name, rep.reference))
    plt, ccd = rep.plt, rep.ccd
    comparison, rows = {}, []
    for direction in DIRECTIONS:
        ref = compute_vof(ccd, direction, rep.reference, rep.params, log=rep.out.detail)
        ref.label = method_label(ref.method)
        tested = []
        for spec in rep.bounded:
            tst = compute_vof(ccd, direction, "fd", dataclasses.replace(rep.params, fd_fit=spec), log=rep.out.detail)
            tst.label = "FD, " + fits.label(spec)
            tst.fit = spec
            stats = compare_results(ref, tst, rep.tol_fc, rep.tol_vof)
            stats["fit"] = fits.label(spec)
            tested.append((tst, stats))
            rows.append(stats)
        comparison[direction] = (ref, tested)
        columns = {"E_minus_EZPL_meV": ref.axis, "VOF %s" % ref.label: ref.lineshape}
        columns.update({"VOF fd %s" % s["fit"]: t.lineshape for t, s in tested})
        rep.out.save_csv(
            name + "/vof_%s.csv" % direction, columns, "%s: VOF of the reference and of fd on every fit." % direction
        )

    rep.out.save_fig(
        fits_vof_figure(rep, comparison, "VOF of FD on each fit, reference: %s" % ref.label),
        name + "/vof_fits.png",
        "Per direction: FD on every bounded fit overlaid on the harmonic reference (%s), with the signed and absolute "
        "differences relative to the reference peak. The differences are the anharmonicity of each fit, not errors. "
        "Spline variants share one colour and differ by line style." % ref.label,
    )

    fig, axes = new_figure(plt, 1, 2, 3.2)
    for ax, state in zip(axes[0], STATES):
        phonon = ccd["%s_phonon" % state]
        values = []
        for spec in rep.bounded:
            energies = rep.wells[(state, spec)][1]
            v = np.arange(len(energies))
            dev = (energies - energies[0]) * 1e3 - v * phonon
            values.append(dev)
            st = plots.fit_style(spec, rep.bounded)
            ax.plot(v, dev, marker=".", ms=3, lw=1.0, color=st["color"], ls=st["ls"], label="FD, " + fits.label(spec))
        ax.axhline(0.0, color="0.5", lw=0.6)
        ax.set_yscale("symlog", linthresh=1.0)
        tight_symlog_limits(ax, values, 1.0)
        ax.yaxis.set_minor_locator(matplotlib.ticker.NullLocator())
        ax.set_title(r"%s well, $\hbar\omega$ = %.2f meV" % (STATE_LABELS[state], phonon))
        ax.set_xlabel("Vibrational level $v$")
        integer_ticks(ax)
        plots.style(ax)
    axes[0, 0].set_ylabel(r"$E_v - E_0 - v\hbar\omega$ (meV)")
    for a, letter in zip(axes[0], "ab"):
        plots.panel_label(a, letter)
    shared_legend(fig, axes, ncol=2)
    rep.out.save_fig(
        fig,
        name + "/level_ladders.png",
        "Per well and fit: how far the FD level ladder departs from the harmonic v hbar omega (symmetric-log axis, "
        "linear within +-1 meV). With the default line positions the VOF still places its lines on the harmonic "
        "ladder; only F_ij follows the fit.",
    )

    keys = (
        "direction",
        "fit",
        "fc_max",
        "vof_max",
        "vof_rms",
        "zpl_weight_ref",
        "zpl_weight",
        "centroid_ref_mev",
        "centroid_shift_mev",
        "closure",
        "seconds",
    )
    rep.out.save_table(
        name + "/summary.csv",
        [{k: r[k] for k in keys} for r in rows],
        "Per fit and direction: F_ij "
        "and VOF difference from the harmonic reference, ZPL weight F_00^2, shift of the lineshape centroid.",
    )
    columns = [
        ("direction", "direction", "%s"),
        ("fit", "fd on fit", "%s"),
        ("fc_max", "max dF_ij", "%.2e"),
        ("vof_max", "max dVOF / peak", "%.2e"),
        ("zpl_weight", "F_00^2", "%.5f"),
        ("centroid_shift_mev", "centroid shift (meV)", "%.2f"),
    ]
    for line in text_table(rows, columns):
        rep.out.log("  " + line)
    rep.headline += (
        [
            "## Finite differences on each fit vs the harmonic reference",
            "",
            "Harmonic F_00^2 = %.5f (emission). Differences here are the anharmonicity of each fit, not errors."
            % rows[0]["zpl_weight_ref"],
            "",
        ]
        + markdown_table(rows, columns)
        + [""]
    )


# ---------------------------------------------------------------------------------- 07


def cubic_potential(ccd, state, spec, scale=1.0):
    """The cubic fit ``spec`` of ``state`` as a callable, with its b multiplied by ``scale``."""

    p = fits.make_potential(ccd, state, spec).params
    center = ccd_data.center(ccd, state)
    return lambda q: p["a"] * (q - center) ** 2 + scale * p["b"] * (q - center) ** 3 + p["c"]


def cubic_exact_reason(rep, scale=1.0):
    """None if both cubic wells with b x ``scale`` are bound on the fd grid of the VOF solves, else why not.

    The turnover of a(Q-Q0)^2 + s b(Q-Q0)^3 is at Q0 - 2a/(3 s b), so it moves out as 1/s.
    """

    lo, hi = rep.fd_q_range()
    for state in STATES:
        p = fits.make_potential(rep.ccd, state, rep.params.cubic_fit).params
        center = ccd_data.center(rep.ccd, state)
        turnover = center + (p["turnover"] - center) / scale if scale else np.inf
        if lo < turnover < hi:
            return "the %s cubic%s turns over at Q = %.2f, inside the grid %.2f..%.2f: no bound states in the box" % (
                STATE_NAMES[state],
                "" if scale == 1 else " with b x %s" % fraction_label(scale),
                turnover,
                lo,
                hi,
            )
    return None


def cubic_exact_scales(rep):
    """The b scales of the b -> 0 scan whose cubic has an exact solution, and the skipped ones with the reason.

    Scales whose turnover lies inside the grid are dropped; if fewer than ``MIN_EXACT_SCALES`` remain, smaller
    scales are added by halving the smallest one, so the slopes can always be fitted.
    """

    exact = [s for s in rep.cubic_scales if cubic_exact_reason(rep, s) is None]
    skipped = {s: cubic_exact_reason(rep, s) for s in rep.cubic_scales if s not in exact}
    s = min(rep.cubic_scales)
    while len(exact) < MIN_EXACT_SCALES and s > SMALLEST_SCALE:
        s /= 2
        if cubic_exact_reason(rep, s) is None:
            exact.append(s)
    return sorted(exact, reverse=True), skipped


def well_grid(rep, state):
    """The physical grid of one well's fd solve in ``solve_well``."""

    n = rep.order(state)
    center = ccd_data.center(rep.ccd, state)
    span = rep.turning_point(state, n - 1) + rep.params.fd_margin
    h = getattr(rep.params, "fd_h", None)
    npts = int(np.ceil(2 * span / h)) + 1 if h else rep.params.fd_npts
    return np.linspace(center - span, center + span, npts)


def cubic_solve(rep, direction, scale=1.0, orders=None, exact=True):
    """Perturbative, exact (fd; None unless ``exact``) and harmonic solvers of ``direction`` on the cubic wells,
    b scaled by ``scale``."""

    ccd, spec = rep.ccd, rep.params.cubic_fit
    initial, final = direction_states(direction)
    ni, nf = orders or rep.params.orders(direction)
    q_i, q_f = ccd_data.center(ccd, initial), ccd_data.center(ccd, final)
    b = {s: scale * fits.make_potential(ccd, s, spec).params["b"] for s in STATES}
    pt = config_coord_1d_vof_solver(ccd["%s_phonon" % initial], ccd["%s_phonon" % final], ccd["delta_q"])
    with contextlib.redirect_stdout(io.StringIO()):
        pt.compute_franck_condon_integrals_cubic(
            ni,
            nf,
            b[initial],
            b[final],
            q_i=q_i,
            q_f=q_f,
            normalize=rep.params.cubic_normalize,
            **quadrature_kwargs(rep.params),
        )
    fd = None
    if exact:
        fd = config_coord_1d_vof_solver(ccd["%s_phonon" % initial], ccd["%s_phonon" % final], ccd["delta_q"])
        fd.compute_franck_condon_integrals_fd(
            ni,
            nf,
            cubic_potential(ccd, initial, spec, scale),
            cubic_potential(ccd, final, spec, scale),
            q_i=q_i,
            q_f=q_f,
            **fd_kwargs(rep.params),
        )
    harmonic = config_coord_1d_vof_solver(ccd["%s_phonon" % initial], ccd["%s_phonon" % final], ccd["delta_q"])
    harmonic.compute_franck_condon_integrals_hermite(ni, nf, **quadrature_kwargs(rep.params))
    return pt, fd, harmonic


def lineshape(rep, solver):
    with contextlib.redirect_stdout(io.StringIO()):
        return solver.bulid_fc_lsp(
            eneaxis=rep.params.axis(),
            temp=rep.params.temp,
            sigma=list(rep.params.sigma),
            zpl_lorentzian=rep.params.zpl_lorentzian,
            gamma=rep.params.gamma,
        )[1]


def cubic_fit_figure(rep, name):
    """Every cubic-fit option per well: near the samples, over the full fd grid, and their residuals."""

    plt, ccd = rep.plt, rep.ccd
    q_lo, q_hi = rep.fd_q_range()
    specs = list(rep.cubic_specs)
    fig, axes = new_figure(plt, 2, 3, 5.4)
    rows = []
    for row, state in enumerate(STATES):
        q, e = ccd_data.samples(ccd, state)
        center, phonon = ccd_data.center(ccd, state), ccd["%s_phonon" % state]
        q_near, q_grid = np.linspace(q.min() - 0.05, q.max() + 0.05, 600), np.linspace(q_lo, q_hi, 2001)
        harmonic = fits.make_potential(ccd, state, "parabola")
        ax_fit, ax_full, ax_res = axes[row]
        for ax, qq in ((ax_fit, q_near), (ax_full, q_grid)):
            ax.plot(q, e, "o", color="black", ms=3.5, zorder=5, label="images (computed energies)")
            ax.plot(qq, harmonic(qq), color="0.6", lw=2.2, alpha=0.6, label=r"parabola (harmonic $H_0$)")
        ax_full.axvspan(q.min(), q.max(), color=plots.OKABE_ITO["yellow"], alpha=0.3, lw=0, label="sampled range")
        residuals = []
        for spec in specs:
            pot = fits.make_potential(ccd, state, spec)
            p = pot.params
            xi = config_coord_1d_vof_solver.anharmonic_strength(p["b"], phonon)
            mixing = 3 * abs(xi) * np.arange(1, rep.order(state) + 1) ** 1.5
            st = plots.fit_style(spec, specs)
            used = spec == rep.params.cubic_fit
            kw = dict(color=st["color"], ls=st["ls"], lw=1.8 if used else 1.0)
            label = fits.label(spec) + (" (used)" if used else "")
            ax_fit.plot(q_near, pot(q_near), label=label, **kw)
            ax_full.plot(q_grid, pot(q_grid), **kw)
            ax_res.plot(q, pot.residuals() * 1e3, marker="o", ms=2.5, **dict(kw, lw=0.8))
            residuals.append(pot.residuals() * 1e3)
            rows.append(
                {
                    "state": STATE_NAMES[state],
                    "fit": spec,
                    "label": fits.label(spec),
                    "images": fits.images(ccd, state, spec),
                    "used_by_method": used,
                    "a": p["a"],
                    "a_fixed": p["a_fixed"],
                    "b": p["b"],
                    "turnover_q": p["turnover"],
                    "xi": xi,
                    "max_mixing": float(mixing.max()),
                    "n_valid": int(np.sum(mixing <= 0.1)) - 1,
                    "rms_all_images_meV": pot.rms() * 1e3,
                    "rms_fitted_meV": pot.rms(pot.index) * 1e3,
                    "description": pot.describe(),
                }
            )
            rep.out.detail("  %s %-16s %s; xi = %+.3e" % (STATE_NAMES[state], spec, pot.describe(), xi))
        for wall in (q_lo, q_hi):
            ax_full.axvline(wall, color="black", lw=1.4, label="FD grid walls")
        span = e.max() - e.min()
        ax_fit.set_ylim(e.min() - 0.12 * span, e.max() + 0.15 * span)
        ax_full.set_ylim(e.min() - 0.5, 1.15 * max(harmonic(q_lo), harmonic(q_hi)))
        ax_full.set_xlim(q_lo - 0.25, q_hi + 0.25)
        ax_res.axhline(0.0, color="0.5", lw=0.6)
        symlog_if_wide(ax_res, residuals, 20.0, 1.0)
        ax_fit.set_title("%s: near the images" % STATE_LABELS[state])
        ax_full.set_title("%s: full FD grid" % STATE_LABELS[state])
        ax_res.set_title("%s: residuals" % STATE_LABELS[state])
        ax_fit.set_ylabel(plots.ENERGY_LABEL)
        ax_res.set_ylabel("Fit %s computed (meV)" % MINUS)
        for a in axes[row]:
            a.set_xlabel(plots.Q_LABEL)
            plots.style(a)
    for a, letter in zip(axes.ravel(), "abcdef"):
        plots.panel_label(a, letter)
    shared_legend(fig, axes, ncol=3)
    fig.suptitle(r"Cubic fits $V = a(Q-Q_0)^2 + b(Q-Q_0)^3 + c$")
    rep.out.save_fig(
        fig,
        name + "/cubic_fits.png",
        "Per well: the cubic-fit options (%s; a fixed to the parabola's or fitted, all images or the parabola's) "
        "against the images (a, d), over the full FD grid with its walls (b, e), and their residuals (c, f). The fit "
        "used by the method (%s) is drawn thicker; b, xi, turnover and mixing are in cubic_fits.csv."
        % (", ".join(specs), rep.params.cubic_fit),
    )
    rep.out.save_table(
        name + "/cubic_fits.csv",
        rows,
        "Per well and cubic fit: a, b, turnover, xi, the largest "
        "first-order mixing 3|xi|(n+1)^1.5, the highest level where it stays <= 0.1, residuals.",
    )
    return rows


def scan_note(scales, skipped, breakdown):
    """Sentence on the b scales left out of the scan (turnover inside the grid), or "" if none were."""

    if not skipped:
        return ""
    return (
        " The cubic with b x %s turns over inside the fd grid and has no exact solution, so the scan uses b x %s, and "
        "(c) shows the exact eigenstates at b x %s, with the mixing at the fitted b dotted."
        % (", ".join(fraction_label(s) for s in skipped), ", ".join(fraction_label(s) for s in scales),
           fraction_label(breakdown))
    )  # fmt: skip


def cubic_validation_figure(rep, name):
    """b -> 0 scaling against the exact solution, the norms, and where first order breaks down."""

    plt, ccd = rep.plt, rep.ccd
    block = 8
    scan, solved = [], {}
    scales, skipped = cubic_exact_scales(rep)
    for scale, reason in skipped.items():
        rep.out.log("  b x %s left out of the b -> 0 scan: %s" % (fraction_label(scale), reason))
    added = [s for s in scales if s not in rep.cubic_scales]
    if added:
        rep.out.log("  b scales added so the scan has %d exact points: %s"
                    % (MIN_EXACT_SCALES, ", ".join(fraction_label(s) for s in added)))  # fmt: skip
    for scale in scales:
        small_pt, small_fd, small_h = cubic_solve(rep, "emission", scale, orders=(12, 12))
        pt, fd, harmonic = cubic_solve(rep, "emission", scale)
        solved[scale] = (pt, fd)
        ls_pt, ls_fd, ls_h = lineshape(rep, pt), lineshape(rep, fd), lineshape(rep, harmonic)
        peak = np.max(ls_fd)
        scan.append(
            {
                "b_scale": scale,
                "xi_i": pt.cubic_xi_i,
                "xi_f": pt.cubic_xi_f,
                "dF_pt_vs_exact_block": float(np.max(np.abs(small_pt.fc_ints - small_fd.fc_ints)[:block, :block])),
                "dF_harmonic_vs_exact_block": float(np.max(np.abs(small_h.fc_ints - small_fd.fc_ints)[:block, :block])),
                "dVOF_pt_vs_exact": float(np.max(np.abs(ls_pt - ls_fd)) / peak),
                "dVOF_harmonic_vs_exact": float(np.max(np.abs(ls_h - ls_fd)) / peak),
            }
        )
        rep.out.detail(
            "  b x %-6g dF(PT-exact) %.2e  dF(harmonic-exact) %.2e  dVOF(PT-exact) %.2e  dVOF(harmonic-exact) %.2e"
            % (scale, *[scan[-1][k] for k in list(scan[-1])[3:]])
        )
    zero_pt, _, zero_h = cubic_solve(rep, "emission", 0.0, exact=False)
    b0_dev = float(np.max(np.abs(zero_pt.fc_ints - zero_h.fc_ints)))
    full_pt = solved[1.0][0] if 1.0 in solved else cubic_solve(rep, "emission", 1.0, exact=False)[0]
    # Panel (c) needs exact eigenstates: at the fitted b if they exist, else at the largest scale that has them.
    breakdown = max(scales)
    breakdown_pt, breakdown_fd = solved[breakdown]

    fig, axes = new_figure(plt, 1, 3, 3.9)
    axes = axes[0]
    s = np.array([r["b_scale"] for r in scan])
    ax = axes[0]
    colors = (plots.OKABE_ITO["purple"], "0.45")
    for key, label, marker, color, fill in (
        ("dF_pt_vs_exact_block", r"$F_{ij}$, PT", "o", colors[0], True),
        ("dF_harmonic_vs_exact_block", r"$F_{ij}$, harmonic", "s", colors[1], True),
        ("dVOF_pt_vs_exact", "VOF, PT", "^", colors[0], False),
        ("dVOF_harmonic_vs_exact", "VOF, harmonic", "v", colors[1], False),
    ):
        y = np.array([r[key] for r in scan])
        slope = np.polyfit(np.log(s), np.log(y), 1)[0] if len(s) > 1 else float("nan")
        ax.loglog(
            s,
            y,
            marker=marker,
            color=color,
            mfc=color if fill else "white",
            ls="-" if fill else "--",
            lw=0.9,
            ms=4,
            label="%s (slope %.2f)" % (label, slope),
        )
    ax.xaxis.set_major_locator(matplotlib.ticker.FixedLocator(s))
    ax.xaxis.set_major_formatter(matplotlib.ticker.FixedFormatter([fraction_label(v) for v in s]))
    ax.xaxis.set_minor_locator(matplotlib.ticker.NullLocator())
    ax.set_xlabel(r"$b\,/\,b_\mathrm{fit}$")
    ax.set_ylabel("max |difference from exact|")
    ax.set_title(r"$b \rightarrow 0$ scaling")

    ax = axes[1]
    for role, norms, xi, state in (
        ("initial", full_pt.cubic_norms_i, full_pt.cubic_xi_i, "es"),
        ("final", full_pt.cubic_norms_f, full_pt.cubic_xi_f, "gs"),
    ):
        color = plots.STATE_COLORS[state]
        n = np.arange(len(norms))
        predicted = (config_coord_1d_vof_solver.cubic_correction_coefficients(len(norms), xi) ** 2).sum(axis=1) - 1
        ax.semilogy(n, norms - 1, "o", ms=2.5, color=color, label="%s: computed" % STATE_LABELS[state])
        ax.semilogy(n, predicted, "-", lw=0.9, color=color, label=r"%s: $\sum c^2 - 1$" % STATE_LABELS[state])
    ax.axhline(0.1, color="0.4", ls="--", lw=0.8, label="0.1")
    ax.set_xlabel("Vibrational level $n$")
    integer_ticks(ax)
    ax.set_ylabel(r"$\langle\tilde\psi_n|\tilde\psi_n\rangle - 1$")
    ax.set_title("Norms")

    ax = axes[2]
    h = grid_spacing(breakdown_fd.fd_grid)
    at = "" if breakdown == 1 else r", $b \times %s$" % fraction_label(breakdown)
    for role, xi, xi_fit, psi_fd, state in (
        ("initial", breakdown_pt.cubic_xi_i, full_pt.cubic_xi_i, breakdown_fd.fd_wavefunctions_i, "es"),
        ("final", breakdown_pt.cubic_xi_f, full_pt.cubic_xi_f, breakdown_fd.fd_wavefunctions_f, "gs"),
    ):
        color = plots.STATE_COLORS[state]
        n_levels = len(psi_fd)
        mixing = 3 * abs(xi) * np.arange(1, n_levels + 1) ** 1.5
        if breakdown != 1:
            ax.semilogy(np.arange(n_levels), 3 * abs(xi_fit) * np.arange(1, n_levels + 1) ** 1.5, ":", color=color,
                        lw=1.0, label=r"%s: $3|\xi|(n+1)^{3/2}$, fitted $b$" % STATE_LABELS[state])  # fmt: skip
        x = breakdown_fd.fd_grid
        center = -ccd["delta_q"] if role == "initial" else 0.0
        coeffs = config_coord_1d_vof_solver.cubic_correction_coefficients(n_levels, xi)
        chi = config_coord_1d_vof_solver.harmonic_wavefunctions(x, n_levels + 2, ccd["%s_phonon" % state], center)
        psi_pt = coeffs @ chi
        psi_pt /= np.sqrt((psi_pt**2).sum(axis=1) * h)[:, None]
        infidelity = 1 - np.abs((psi_pt * psi_fd).sum(axis=1) * h)
        ax.semilogy(
            np.arange(n_levels),
            mixing,
            "-",
            color=color,
            lw=1.0,
            label=r"%s: $3|\xi|(n+1)^{3/2}$%s" % (STATE_LABELS[state], at),
        )
        ax.semilogy(np.arange(n_levels), np.maximum(infidelity, 1e-16), "o", ms=2.5, color=color,
                    label="%s: infidelity%s" % (STATE_LABELS[state], at))  # fmt: skip
    ax.axhline(0.1, color="0.4", ls="--", lw=0.8, label="0.1 (warning threshold)")
    ax.set_xlabel("Vibrational level $n$")
    integer_ticks(ax)
    ax.set_ylabel("Mixing, infidelity")
    ax.set_title("Breakdown of first order" if breakdown == 1 else "Breakdown" + at)
    for a, letter in zip(axes, "abc"):
        plots.style(a)
        a.yaxis.set_minor_locator(matplotlib.ticker.NullLocator())
        plots.panel_label(a, letter)
        # Each panel has its own quantities, so its legend sits directly below it, clear of the data.
        a.legend(loc="upper center", bbox_to_anchor=(0.5, -0.2), ncol=1, frameon=False)
    fig.suptitle("First-order cubic correction, %s" % fits.label(rep.params.cubic_fit))
    rep.out.save_fig(
        fig,
        name + "/validation.png",
        "Emission, against the exact solution of the same cubic (FD): (a) max |difference| of F_ij (i, j < %d, basis "
        "12 x 12) and of the VOF / peak as b is scaled towards 0; the PT error should fall as b^2 (slope 2) and the "
        "anharmonic effect (harmonic vs exact) as b (slope 1); at b = 0, max |F_PT - F_Hermite| = %.1e. (b) norm - 1 of "
        "the corrected states vs the predicted sum of squared PT coefficients - 1. (c) the first-order mixing "
        "3|xi|(n+1)^1.5 and the infidelity against the exact eigenstates per level: first order breaks down where "
        "the mixing passes 0.1 (dashed).%s" % (block, b0_dev, scan_note(scales, skipped, breakdown)),
    )
    rep.out.save_table(name + "/validation_b_scaling.csv", scan, "The b-scaling scan behind validation.png.")
    return scan, b0_dev, full_pt


def cubic_wavefunction_figures(rep, name):
    """Harmonic, first-order corrected and exact (fd on the cubic) wavefunctions, per well.

    A well whose cubic turns over inside its grid has no exact states; it shows harmonic and PT, and PT - harmonic.
    """

    ccd = rep.ccd
    purple = plots.METHOD_COLORS["cubic"]
    for state in STATES:
        n = rep.order(state)
        phonon, center = ccd["%s_phonon" % state], ccd_data.center(ccd, state)
        pot = fits.make_potential(ccd, state, rep.params.cubic_fit)
        q = well_grid(rep, state)
        exact = not q[0] < pot.params["turnover"] < q[-1]
        if exact:
            q, _, psi_fd = solve_well(rep, state, rep.params.cubic_fit)
        else:
            rep.out.log(
                "  %s wavefunctions: the cubic turns over at Q = %.2f, inside the well's grid %.2f..%.2f; no exact "
                "states, PT is shown against the harmonic states" % (STATE_NAMES[state], pot.params["turnover"], q[0],
                                                                       q[-1])
            )  # fmt: skip
        xi = config_coord_1d_vof_solver.anharmonic_strength(pot.params["b"], phonon)
        chi = config_coord_1d_vof_solver.harmonic_wavefunctions(q, n + 2, phonon, center)
        psi_pt = config_coord_1d_vof_solver.cubic_correction_coefficients(n, xi) @ chi
        shown = rep.levels(n)
        rows = []
        for v in range(n):
            row = {"state": STATE_NAMES[state], "v": v, "mixing": 3 * abs(xi) * (v + 1) ** 1.5}
            if exact:
                row["max_abs_pt_minus_exact"] = float(np.max(np.abs(psi_pt[v] - psi_fd[v])))
                row["max_abs_harmonic_minus_exact"] = float(np.max(np.abs(chi[v] - psi_fd[v])))
            else:
                row["max_abs_pt_minus_harmonic"] = float(np.max(np.abs(psi_pt[v] - chi[v])))
                row["max_abs_harmonic"] = float(np.max(np.abs(chi[v])))
            rows.append(row)
        for r in rows:
            if r["v"] in shown:
                rep.out.detail(
                    "  %s v = %2d: mixing %.3g, " % (STATE_NAMES[state], r["v"], r["mixing"])
                    + (
                        "max |PT - exact| %.1e, max |harmonic - exact| %.1e"
                        % (r["max_abs_pt_minus_exact"], r["max_abs_harmonic_minus_exact"])
                        if exact
                        else "max |PT - harmonic| %.1e (max |harmonic| %.2f)"
                        % (r["max_abs_pt_minus_harmonic"], r["max_abs_harmonic"])
                    )
                )
        if exact:
            left = lambda v: [
                (q, chi[v], dict(color="0.7", lw=2.6, label="harmonic")),
                (q, psi_fd[v], dict(color="black", lw=1.0, label="exact (FD on the cubic)")),
                (q, psi_pt[v], dict(color=purple, ls="--", lw=1.0, label="first-order PT")),
            ]
            right = lambda v: [
                (q, chi[v] - psi_fd[v], dict(color="0.55", lw=0.9, label="harmonic %s exact" % MINUS)),
                (q, psi_pt[v] - psi_fd[v], dict(color=purple, lw=0.9, label="PT %s exact" % MINUS)),
            ]
            titles = (r"$\psi_v$", "Difference from exact")
            description = (
                "harmonic, first-order corrected and exact (FD on the same cubic) wavefunctions for selected levels "
                "(left) and their differences from exact"
            )
        else:
            left = lambda v: [
                (q, chi[v], dict(color="0.7", lw=2.6, label="harmonic")),
                (q, psi_pt[v], dict(color=purple, ls="--", lw=1.0, label="first-order PT")),
            ]
            right = lambda v: [(q, psi_pt[v] - chi[v], dict(color=purple, lw=0.9, label="PT %s harmonic" % MINUS))]
            titles = (r"$\psi_v$", "Correction (PT %s harmonic)" % MINUS)
            description = (
                "harmonic and first-order corrected wavefunctions for selected levels (left) and the correction PT - "
                "harmonic; the cubic turns over at Q = %.2f, inside the grid, so there are no exact states to compare "
                "with" % pot.params["turnover"]
            )
        fig = wavefunction_rows_figure(
            rep,
            state,
            shown,
            lambda v: np.abs(q - center) < rep.turning_point(state, v) + 0.45,
            left,
            right,
            titles,
            r"%s well: cubic correction, $\xi = %s$" % (STATE_LABELS[state], sci_tex(xi)),
        )
        rep.out.save_fig(
            fig,
            name + "/wavefunctions_%s.png" % state,
            "%s, %s: %s (right, each row scaled by the power of ten on its axis). The mixing 3|xi|(n+1)^1.5 and the "
            "maxima per level are in cubic_wavefunctions_%s.csv."
            % (STATE_NAMES[state], fits.label(rep.params.cubic_fit), description, state),
        )
        rep.out.save_table(
            name + "/cubic_wavefunctions_%s.csv" % state,
            rows,
            "%s: per level, the first-order mixing and %s."
            % (
                STATE_NAMES[state],
                "max |PT - exact| and |harmonic - exact|" if exact else "max |PT - harmonic| (no exact states)",
            ),
        )


def section_cubic(rep):
    """The first-order cubic anharmonic correction (method ``cubic``): fits, validation, wavefunctions, F_ij, VOF.

    The VOF and F_ij reference is the exact (fd) solution of the same cubic. If the cubic turns over inside the fd
    grid there is none; the reference is then the harmonic VOF, and fd on every bounded polynomial fit is shown
    alongside the PT as the anharmonic polynomial wells that do have exact solutions.
    """

    name = "07_cubic_perturbation"
    rep.out.section(name, SECTIONS[name])
    rep.out.log("\n[%s] first-order cubic correction, %s" % (name, fits.label(rep.params.cubic_fit)))
    rep.out.log(
        "  cubic fits drawn: %s; b scales: %s"
        % (", ".join(rep.cubic_specs), ", ".join("%g" % s for s in rep.cubic_scales))
    )
    no_exact = cubic_exact_reason(rep)
    polys = [s for s in rep.bounded if fits.parse_spec(s)[0] == "polynomial"]
    if no_exact:
        rep.out.log("  no exact solution of the cubic: %s" % no_exact)
        rep.out.log("  VOF / F_ij reference: the harmonic %s; also compared: fd on %s"
                    % (method_label(rep.reference), ", ".join(polys) or "no bounded polynomial fit"))  # fmt: skip
    fit_rows = cubic_fit_figure(rep, name)
    scan, b0_dev, _ = cubic_validation_figure(rep, name)
    cubic_wavefunction_figures(rep, name)

    if no_exact:
        ref_name, ref_column = "harmonic reference", "VOF harmonic (%s)" % method_label(rep.reference)
        tests = [("cubic", rep.params, method_label("cubic"))]
        tests += [("fd", dataclasses.replace(rep.params, fd_fit=s), "FD on %s" % fits.label(s)) for s in polys]
    else:
        ref_name, ref_column = "exact cubic solution", "VOF exact (fd on the cubic)"
        tests = [("cubic", rep.params, method_label("cubic")),
                 (rep.reference, rep.params, method_label(rep.reference) + ", harmonic")]  # fmt: skip
    comparison, rows = {}, []
    for direction in DIRECTIONS:
        if no_exact:
            ref = compute_vof(rep.ccd, direction, rep.reference, rep.params, log=rep.out.detail)
            ref.label = method_label(rep.reference) + ", harmonic"
        else:
            ref = compute_vof(rep.ccd, direction, "fd", dataclasses.replace(rep.params, fd_fit=rep.params.cubic_fit),
                              log=rep.out.detail)  # fmt: skip
            ref.label = "Exact (FD on the cubic)"
        tested = []
        for method, params, label in tests:
            tst = compute_vof(rep.ccd, direction, method, params, log=rep.out.detail)
            tst.label = label
            stats = compare_results(ref, tst, rep.tol_fc, rep.tol_vof)
            tested.append((tst, stats))
            rows.append(stats)
        comparison[direction] = (ref, tested)
        columns = {"E_minus_EZPL_meV": ref.axis, ref_column: ref.lineshape}
        columns.update({"VOF %s" % t.label: t.lineshape for t, _ in tested})
        rep.out.save_csv(name + "/vof_%s.csv" % direction, columns, "%s: VOF of the %s and of %s." % (
            direction, ref_name, ", ".join(t.label for t, _ in tested)))  # fmt: skip
    if no_exact:
        rep.out.save_fig(plots.vof_comparison_figure(rep.plt, comparison, rep.params, None,
                                                     "VOF of the first-order cubic correction, reference: harmonic"),
                         name + "/vof_cubic.png", "Per direction: the first-order cubic correction and FD on the "
                         "bounded polynomial fits (%s) against the harmonic VOF, with signed and absolute differences "
                         "relative to the harmonic peak. The cubic has no exact solution (%s)."
                         % (", ".join(polys) or "none", no_exact))  # fmt: skip
        rep.out.save_fig(integer_maps(plots.fc_error_figure(rep.plt, comparison, r"$F_{ij}$ difference from harmonic")),
                         name + "/fc_error_maps.png", "log10 |F_ij - F_ij(harmonic)| of the first-order correction "
                         "and of FD on the bounded polynomial fits.")  # fmt: skip
    else:
        rep.out.save_fig(plots.vof_comparison_figure(rep.plt, comparison, rep.params, None,
                                                     "VOF on the cubic wells, reference: exact FD"),
                         name + "/vof_cubic.png", "Per direction: the exact VOF on the cubic wells (FD), the "
                         "first-order correction and the harmonic reference, with signed and absolute differences "
                         "relative to the exact peak.")  # fmt: skip
        rep.out.save_fig(integer_maps(plots.fc_error_figure(rep.plt, comparison,
                                                            r"$F_{ij}$ error vs the exact cubic solution")),
                         name + "/fc_error_maps.png", "log10 |F_ij - F_ij(exact cubic)| of the first-order correction "
                         "and of the harmonic reference.")  # fmt: skip
    rep.out.save_table(name + "/summary.csv", [{k: r[k] for k in CSV_KEYS} for r in rows], "F_ij and VOF differences "
                       "from the %s, per method and direction." % ref_name)  # fmt: skip
    closure = []
    for direction in DIRECTIONS:
        ref, tested = comparison[direction]
        for res in [ref] + [t for t, _ in tested]:
            closure.append({"direction": direction, "method": res.label,
                            "sum_F0j2": float(np.sum(res.solver.fc_ints[0] ** 2)), "centroid_meV": res.centroid_mev})  # fmt: skip
    rep.out.save_table(name + "/closure.csv", closure, "Per direction and method: sum_j F_0j^2 (1 for a complete "
                       "orthonormal basis; first-order PT exceeds it where it breaks down) and the VOF centroid.")  # fmt: skip

    used = [r for r in fit_rows if r["used_by_method"]]
    columns = [("state", "well", "%s"), ("label", "cubic fit", "%s"), ("b", "b (eV/amu^1.5 A^3)", "%+.5f"),
               ("xi", "xi", "%+.3e"), ("max_mixing", "max mixing", "%.3g"), ("n_valid", "n_valid", "%d"),
               ("turnover_q", "turnover Q", "%.1f")]  # fmt: skip
    for line in text_table(fit_rows, columns):
        rep.out.log("  " + line)
    vof_columns = [("direction", "direction", "%s"),
                   ("test", "vs harmonic" if no_exact else "vs exact (fd on the cubic)", "%s"),
                   ("fc_max", "max dF_ij", "%.2e"), ("fc_low", "max dF_ij (i,j<30)", "%.2e"),
                   ("vof_max", "max dVOF / peak", "%.2e")]  # fmt: skip
    for line in text_table(rows, vof_columns):
        rep.out.log("  " + line)
    closure_columns = [("direction", "direction", "%s"), ("method", "method", "%s"),
                       ("sum_F0j2", "sum_j F_0j^2", "%.4f"), ("centroid_meV", "VOF centroid (meV)", "%.1f")]  # fmt: skip
    for line in text_table(closure, closure_columns):
        rep.out.log("  " + line)
    slopes = {}
    s = np.log([r["b_scale"] for r in scan])
    for key in ("dF_pt_vs_exact_block", "dVOF_pt_vs_exact", "dF_harmonic_vs_exact_block"):
        slopes[key] = float(np.polyfit(s, np.log([r[key] for r in scan]), 1)[0]) if len(s) > 1 else float("nan")
    scales = ", ".join(fraction_label(r["b_scale"]) for r in scan)
    rep.out.log("  b -> 0 (b x %s): PT error slope %.2f (F_ij, i,j<8) / %.2f (VOF), effect slope %.2f; at b = 0 "
                "|F_PT - F_hermite| = %.1e" % (scales, slopes["dF_pt_vs_exact_block"], slopes["dVOF_pt_vs_exact"],
                                               slopes["dF_harmonic_vs_exact_block"], b0_dev))  # fmt: skip
    rep.headline += ["## First-order cubic correction (%s)" % fits.label(rep.params.cubic_fit), ""]
    if no_exact:
        rep.headline += [
            "**No exact solution of the fitted cubic:** %s. The VOF and F_ij below are compared with the harmonic "
            "VOF, alongside FD on the bounded polynomial fits (%s); the b -> 0 validation uses the scaled cubics "
            "that stay bound on the grid." % (no_exact, ", ".join(polys) or "none"),
            "",
        ]
    rep.headline += markdown_table(used, columns) + [""] + markdown_table(rows, vof_columns) + [""]
    rep.headline += markdown_table(closure, closure_columns) + [""]
    valid = " / ".join("%d (%s)" % (r["n_valid"], r["state"]) for r in used)
    rep.headline += [
        "b -> 0 validation (b x %s): the PT error falls as b^%.2f (F_ij, i, j < 8) and b^%.2f (VOF), the anharmonic "
        "effect as b^%.2f; at b = 0 the method equals the Hermite route to %.1e. The mixing 3|xi|(n+1)^1.5 exceeds "
        "0.1 above n = %s, so the F_ij of higher levels are not perturbative; sum_j F_0j^2 in the closure table "
        "departs from 1 when the levels carrying the intensity are among them."
        % (scales, slopes["dF_pt_vs_exact_block"], slopes["dVOF_pt_vs_exact"], slopes["dF_harmonic_vs_exact_block"],
           b0_dev, valid),
        "",
    ]  # fmt: skip


# -------------------------------------------------------------------------------- main


def cubic_spec_arg(text):
    """argparse type of ``--cubic-specs``: a ``cubic[-near][-free]`` fit spec."""

    try:
        kind, _ = fits.parse_spec(text)
    except ValueError as err:
        raise argparse.ArgumentTypeError(str(err))
    if kind != "cubic":
        raise argparse.ArgumentTypeError("%r is not a cubic fit (cubic[-near][-free])" % text)
    return text


def positive(kind):
    def convert(text):
        value = kind(text)
        if value <= 0:
            raise argparse.ArgumentTypeError("must be positive, got %s" % text)
        return value

    return convert


def non_negative_int(text):
    value = int(text)
    if value < 0:
        raise argparse.ArgumentTypeError("must be >= 0, got %s" % text)
    return value


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument(
        "--fits",
        nargs="+",
        default=list(fits.REPORT_FITS),
        metavar="FIT",
        help="fits to compare (default: %s)" % " ".join(fits.REPORT_FITS),
    )
    p.add_argument("--reference", choices=sorted(METHODS), default="pypl", help="reference F_ij method (default pypl)")
    p.add_argument("--skip-convergence", action="store_true", help="skip the grid-convergence scan in 02_fd_grid")
    p.add_argument("--skip-cubic", action="store_true", help="skip 07_cubic_perturbation")
    r = p.add_argument_group("report settings")
    r.add_argument(
        "--wf-levels",
        nargs="+",
        type=non_negative_int,
        default=list(SHOWN_LEVELS),
        metavar="V",
        help="levels drawn individually in 02, 03 and 07, clipped to the basis; the highest level is always added "
        "(default: %s)" % " ".join(map(str, SHOWN_LEVELS)),
    )
    r.add_argument(
        "--convergence-npts",
        nargs="+",
        type=positive(int),
        default=list(CONVERGENCE_NPTS),
        metavar="N",
        help="grid sizes of the FD convergence scan in 02 (default: %s)" % " ".join(map(str, CONVERGENCE_NPTS)),
    )
    r.add_argument(
        "--convergence-stencils",
        nargs="+",
        type=int,
        choices=sorted(FD_STENCILS),
        default=list(CONVERGENCE_STENCILS),
        metavar="S",
        help="FD stencils of the convergence scan, from %s (default: %s)"
        % (", ".join(map(str, sorted(FD_STENCILS))), " ".join(map(str, CONVERGENCE_STENCILS))),
    )
    r.add_argument(
        "--cubic-specs",
        nargs="+",
        type=cubic_spec_arg,
        default=list(CUBIC_SPECS),
        metavar="SPEC",
        help="cubic fits drawn in 07 (default: %s)" % " ".join(CUBIC_SPECS),
    )
    r.add_argument(
        "--cubic-scales",
        nargs="+",
        type=positive(float),
        default=list(CUBIC_SCALES),
        metavar="S",
        help="factors on b of the b -> 0 validation scan in 07, at least two (default: %s)"
        % " ".join("%g" % s for s in CUBIC_SCALES),
    )
    add_vof_args(p)
    if not any(a.dest == "dpi" for a in p._actions):
        plots.add_figure_args(p)
    p.add_argument("--tol-fc", type=float, default=1e-7, help="max |dF_ij| for OK in 05 (default 1e-7)")
    p.add_argument("--tol-vof", type=float, default=1e-8, help="max |dVOF| / peak for OK in 05 (default 1e-8)")
    p.add_argument("--tag", default="", help="name of the run directory (default: 'default')")
    p.add_argument("-v", "--verbose", action="store_true", help="print every method's settings, not only summaries")
    args = p.parse_args()
    if len(args.cubic_scales) < 2:
        p.error("--cubic-scales needs at least two values (the slopes are fitted)")
    if len(args.convergence_npts) < 2:
        p.error("--convergence-npts needs at least two values (the slopes are fitted)")
    return args


def main():
    args = parse_args()
    params = params_from_args(args)
    plots.configure(args)
    try:
        args.fits = list(dict.fromkeys(fits.canonical(spec) for spec in args.fits))
    except ValueError as err:
        raise SystemExit("Invalid fit: %s" % err)
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    report_settings = {
        "wf_levels": sorted(set(args.wf_levels)),
        "convergence_npts": sorted(set(args.convergence_npts)),
        "convergence_stencils": sorted(set(args.convergence_stencils)),
        "cubic_specs": list(dict.fromkeys(args.cubic_specs)),
        "cubic_scales": sorted(set(args.cubic_scales), reverse=True),
        "skip_convergence": args.skip_convergence,
        "skip_cubic": args.skip_cubic,
    }
    start = time.perf_counter()
    out = OutputDir("vof_report", args.tag or "default", title="1D CCD VOF report", verbose=args.verbose)
    for name, description in SECTIONS.items():
        out.section(name, description)
    ccd = ccd_data.build_ccd(verbose=False)
    bounded = [s for s in args.fits if all(fits.make_potential(ccd, st, s).bounded for st in STATES)]
    skipped = [s for s in args.fits if s not in bounded]
    out.save_json(
        "settings.json",
        {
            "fits": args.fits,
            "bounded_fits": bounded,
            "reference": args.reference,
            "tol_fc": args.tol_fc,
            "tol_vof": args.tol_vof,
            "report": report_settings,
            "figures": {k: plots.FIGURE[k] for k in plots.FIGURE},
            "params": dataclasses.asdict(params),
        },
    )
    rep = Report(
        ccd,
        params,
        out,
        plt,
        args.fits,
        bounded,
        args.reference,
        args.tol_fc,
        args.tol_vof,
        shown_levels=tuple(report_settings["wf_levels"]),
        convergence_npts=tuple(report_settings["convergence_npts"]),
        convergence_stencils=tuple(report_settings["convergence_stencils"]),
        cubic_specs=tuple(report_settings["cubic_specs"]),
        cubic_scales=tuple(report_settings["cubic_scales"]),
    )

    out.log(
        "1D CCD VOF report: delta_q = %.6f amu^0.5 A, GS %.5f meV, ES %.5f meV"
        % (ccd["delta_q"], ccd["gs_phonon"], ccd["es_phonon"])
    )
    for line in params.describe():
        out.log("  " + line)
    out.log(
        "  report: levels drawn %s; convergence scan %s points x %s-point stencils%s; cubic fits %s, b scales %s%s"
        % (
            report_settings["wf_levels"],
            report_settings["convergence_npts"],
            "/".join(map(str, report_settings["convergence_stencils"])),
            " (skipped)" if args.skip_convergence else "",
            ", ".join(report_settings["cubic_specs"]),
            ", ".join("%g" % s for s in report_settings["cubic_scales"]),
            " (section 07 skipped)" if args.skip_cubic else "",
        )
    )
    out.log(
        "  figures: %g pt base font, %d dpi, %s"
        % (plots.FIGURE["font_size"], plots.FIGURE["dpi"], ", ".join(plots.FIGURE["formats"]))
    )
    if skipped:
        out.log("  unbounded fits (01 only, no bound states for 02-06): %s" % ", ".join(skipped))

    section_potential_fits(rep)
    section_fd_grid(rep, convergence=not args.skip_convergence)
    section_wavefunctions(rep)
    methods_ok = section_methods(rep)
    section_vof_fits(rep)
    if not args.skip_cubic:
        section_cubic(rep)

    verdict = "PASS" if methods_ok else "FAIL"
    out.add_summary(
        "CCD: delta_q = %.6f amu^0.5 A, GS hbar w = %.5f meV, ES hbar w = %.5f meV. %s."
        % (ccd["delta_q"], ccd["gs_phonon"], ccd["es_phonon"], "; ".join(params.describe())),
        "",
        "Report settings: levels drawn %s; convergence scan %s points, %s-point stencils; cubic fits %s; b scales %s."
        % (
            report_settings["wf_levels"],
            report_settings["convergence_npts"],
            "/".join(map(str, report_settings["convergence_stencils"])),
            ", ".join(report_settings["cubic_specs"]),
            ", ".join("%g" % s for s in report_settings["cubic_scales"]),
        ),
        "",
        "**Like-for-like method test (section 05): %s.** Fits compared: %s%s."
        % (
            verdict,
            ", ".join(args.fits),
            "; unbounded, so only in section 01: %s" % ", ".join(skipped) if skipped else "",
        ),
        "",
        "Sections: " + ", ".join("`%s/`" % s for s in SECTIONS),
        "",
    )
    out.add_summary(*rep.headline)
    out.log("\nMethod test (05): %s. Report finished in %.0f s." % (verdict, time.perf_counter() - start))
    out.close()
    raise SystemExit(0 if methods_ok else 1)


if __name__ == "__main__":
    main()
