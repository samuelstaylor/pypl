#!/usr/bin/env python3
"""Thorough, sectioned comparison of every route to the 1D CCD vibrational overlap function (VOF).

    python testing/vof_report.py                          # everything, ~4 min
    python testing/vof_report.py --skip-convergence       # without the grid-convergence scan, ~2.5 min
    python testing/vof_report.py --fits parabola spline morse --tag three-fits
    python testing/vof_report.py --temp 300 --tag 300K

Writes ``testing/outputs/vof_report/<tag or default>/``. Open ``index.md`` first: it lists the
headline numbers and every file with a one-line description. The sections are:

    01_potential_fits   each fit of the DFT images: residuals, implied hbar*omega, and how far it
                        departs from the harmonic well where the high levels live
    02_fd_grid          the finite-difference grid drawn point by point, the wavefunctions it gives,
                        their decay to the walls, and convergence with grid points and stencil
    03_wavefunctions    vibrational wavefunctions: analytic Hermite vs finite differences on each fit
    04_fc_matrix        the Franck-Condon matrix of each method and its error against the reference
    05_vof_methods      VOF of every F_ij method on the harmonic wells (the like-for-like test)
    06_vof_fits         VOF of the finite-difference method on every fit, against the harmonic reference

The harmonic methods (``pypl``, ``hermite``) only know the parabola near each minimum, so in
sections 03 and 06 every non-parabola difference is the anharmonicity of that fit, not an error.
All VOF parameters (``--temp``, ``--sigma``, ``--fd-npts``, ...) are the shared flags of
``harness.methods`` and mean the same as in ``compare_vof.py`` and ``plot_ccd_wavefunctions.py``.
"""

import argparse
import contextlib
import dataclasses
import io
import time

import matplotlib
import numpy as np

from harness import ccd_data, fits, plots
from harness.ccd_data import STATE_LABELS, STATE_NAMES, STATES
from harness.methods import COMPARISON_COLUMNS, CSV_KEYS, DIRECTIONS, METHODS, add_vof_args, compare_results
from harness.methods import compute_vof, direction_states, fd_q_range, params_from_args
from harness.output import OutputDir, markdown_table, text_table
from pypl_ext.config_coord_1d_vof_solver import FD_STENCILS, KINETIC, config_coord_1d_vof_solver, grid_spacing

SECTIONS = {
    "01_potential_fits": "Each fit of the DFT images: residuals, implied hbar*omega, and how far it departs from "
    "the harmonic well out where the high vibrational levels live.",
    "02_fd_grid": "The finite-difference grid point by point, the wavefunctions it gives, their decay to the walls, "
    "and convergence with grid points and stencil.",
    "03_wavefunctions": "Vibrational wavefunctions: analytic Hermite functions vs finite differences on each fit.",
    "04_fc_matrix": "The Franck-Condon matrix of each method and its error against the reference.",
    "05_vof_methods": "VOF of every F_ij method on the harmonic wells: the like-for-like test.",
    "06_vof_fits": "VOF of the finite-difference method on each fit, against the harmonic reference.",
    "07_cubic_perturbation": "The first-order cubic anharmonic correction (method cubic): the cubic fits, its "
    "validation against the exact solution of the same cubic, wavefunctions, F_ij and VOF.",
}
# Levels drawn individually in the grid and wavefunction figures (clipped to the basis size).
SHOWN_LEVELS = (0, 1, 2, 5, 10, 20, 30, 40)
CONVERGENCE_NPTS = (2001, 4001, 8001, 16001, 32001)


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
        return sorted({v for v in SHOWN_LEVELS if v < n} | {n - 1})


# ---------------------------------------------------------------------------------- 01


def section_potential_fits(rep):
    """Fit residuals, implied frequencies and the departure from harmonic beyond the sampled range."""

    name = "01_potential_fits"
    rep.out.section(name, SECTIONS[name])
    rep.out.log("\n[%s] %d fits of the DFT images" % (name, len(rep.fits)))
    plt, ccd = rep.plt, rep.ccd
    fig, axes = plt.subplots(2, 4, figsize=(23, 9.5), gridspec_kw={"width_ratios": [1, 1.25, 1, 1]})
    rows = []
    q_lo, q_hi = rep.fd_q_range()
    q_grid = np.linspace(q_lo, q_hi, 2001)
    for row, state in enumerate(STATES):
        q, e = ccd_data.samples(ccd, state)
        center = ccd_data.center(ccd, state)
        harmonic = fits.make_potential(ccd, state, "parabola")
        used = harmonic.index
        tp = rep.turning_point(state, rep.order(state) - 1)
        q_near = np.linspace(q.min() - 0.05, q.max() + 0.05, 600)
        q_far = np.linspace(center - tp - 0.3, center + tp + 0.3, 1200)
        others = np.setdiff1d(np.arange(len(q)), used)

        ax_fit, ax_full, ax_res, ax_far = axes[row]
        for ax in (ax_fit, ax_full):
            ax.plot(
                q[used],
                e[used],
                "o",
                color="black",
                ms=6,
                zorder=5,
                label="DFT images used by the parabola [%s]" % fits.images(ccd, state, "parabola"),
            )
            ax.plot(q[others], e[others], "o", mfc="white", mec="black", ms=6, zorder=5, label="other DFT images")
        for ax in (ax_full, ax_far):
            ax.axvspan(q.min(), q.max(), color="gold", alpha=0.18, lw=0, label="sampled by DFT")
        wall_values = []
        for spec in rep.fits:
            pot = fits.make_potential(ccd, state, spec)
            c, ls = plots.fit_color(spec), "-" if pot.bounded else "--"
            tag = "" if pot.bounded else " (unbounded)"
            ax_fit.plot(q_near, pot(q_near), color=c, ls=ls, lw=1.2, label=fits.state_label(ccd, state, spec) + tag)
            ax_full.plot(q_grid, pot(q_grid), color=c, ls=ls, lw=1.2, label=fits.state_label(ccd, state, spec) + tag)
            ax_res.plot(q, pot.residuals() * 1e3, marker="o", ms=3.5, color=c, ls=ls, lw=0.9, label=fits.label(spec))
            ax_far.plot(q_far, (pot(q_far) - harmonic(q_far)) * 1e3, color=c, ls=ls, lw=1.2, label=fits.label(spec))
            wall_values.append("%-9s %10.4g %10.4g" % (spec, pot(q_lo), pot(q_hi)))
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
                    "dV_left_turning_point_meV": (pot(center - tp) - harmonic(center - tp)) * 1e3,
                    "dV_right_turning_point_meV": (pot(center + tp) - harmonic(center + tp)) * 1e3,
                    "description": pot.describe(),
                }
            )
            rep.out.detail("  %s %-10s %s" % (STATE_NAMES[state], spec, pot.describe()))

        # Full FD range: the parabola sets the scale; diverging fits run off it, so print V at both walls.
        for wall in (q_lo, q_hi):
            ax_full.axvline(wall, color="black", lw=2)
        top = 1.15 * max(harmonic(q_lo), harmonic(q_hi))
        ax_full.set_ylim(e.min() - 0.08 * (top - e.min()), top)
        ax_full.set_xlim(q_lo - 0.25, q_hi + 0.25)
        ax_full.text(
            0.5,
            0.97,
            "V at the FD walls (eV)\n%-9s %10s %10s\n" % ("fit", "Q=%.2f" % q_lo, "Q=%.2f" % q_hi)
            + "\n".join(wall_values),
            transform=ax_full.transAxes,
            ha="center",
            va="top",
            fontsize=6.5,
            family="monospace",
            bbox={"facecolor": "white", "edgecolor": "gray", "alpha": 0.9},
        )
        ax_full.set_title(
            "%s: every fit over the full FD grid, Q = %.2f..%.2f (walls in black)" % (STATE_LABELS[state], q_lo, q_hi),
            fontsize=10,
        )
        ax_full.set_ylabel("E (eV)")
        for a in (ax_far,):
            a.axvline(center - tp, color="gray", ls=":", lw=1)
            a.axvline(center + tp, color="gray", ls=":", lw=1)
        ax_res.axhline(0.0, color="gray", lw=0.8)
        ax_far.set_yscale("symlog", linthresh=1.0)
        ax_fit.set_title("%s: DFT images and fits" % STATE_LABELS[state])
        ax_res.set_title("%s: residual at every image" % STATE_LABELS[state])
        ax_far.set_title(
            "%s: V_fit - V_parabola out to the v = %d turning points (dotted)"
            % (STATE_LABELS[state], rep.order(state) - 1),
            fontsize=10,
        )
        ax_fit.set_ylabel("E (eV)")
        ax_res.set_ylabel("fit - DFT (meV)")
        ax_far.set_ylabel("V_fit - V_parabola (meV, symlog)")
        for a in axes[row]:
            a.set_xlabel(r"Q (amu$^{1/2}$ $\AA$)")
            plots.style(a)
        ax_fit.legend(fontsize=7, edgecolor="black")
        ax_far.legend(fontsize=7, edgecolor="black", loc="lower center")
    fig.suptitle("Potential fits of the 1D CCD wells", fontsize=13)
    fig.tight_layout()
    rep.out.save_fig(
        fig,
        name + "/fits_overview.png",
        "Per well: the DFT images with every fit near the samples (filled points are the parabola's); every fit "
        "over the full Q range of the finite-difference grid, with the DFT images, the sampled range (gold), the "
        "walls (black) and V at the walls; the residual at each image; and V_fit - V_parabola out to the turning "
        "points of the highest level (symlog). Dashed = unbounded.",
    )

    fig, ax = plt.subplots(figsize=(10, 4.5))
    width = 0.8 / len(rep.fits)
    for k, spec in enumerate(rep.fits):
        values = [r["rms_all_images_meV"] for r in rows if r["fit"] == spec]
        bounded = [r["bounded"] for r in rows if r["fit"] == spec]
        bars = ax.bar(
            np.arange(2) + (k - (len(rep.fits) - 1) / 2) * width,
            values,
            width,
            color=plots.fit_color(spec),
            label=fits.label(spec),
        )
        for bar, ok in zip(bars, bounded):
            if not ok:
                bar.set_hatch("//")
                bar.set_alpha(0.5)
    ax.set_xticks(np.arange(2))
    ax.set_xticklabels(["GS well (%d images)" % ccd["n_gs_images"], "ES well (%d images)" % ccd["n_es_images"]])
    ax.set_yscale("log")
    ax.set_ylabel("rms residual over all images (meV)")
    ax.set_title("Fit quality (hatched = unbounded, unusable for the Schroedinger solve)")
    ax.legend(fontsize=8, edgecolor="black", ncol=2)
    plots.style(ax)
    rep.out.save_fig(fig, name + "/fit_residual_rms.png", "rms residual over all images, per fit and well.")

    for state in STATES:
        q, e = ccd_data.samples(ccd, state)
        columns = {"Q": q, "E_DFT_eV": e}
        columns.update({"residual_%s_meV" % s: fits.make_potential(ccd, state, s).residuals() * 1e3 for s in rep.fits})
        rep.out.save_csv(
            name + "/residuals_%s.csv" % state, columns, "%s: fit - DFT at every image, per fit." % STATE_NAMES[state]
        )
    rep.out.save_table(
        name + "/fit_summary.csv",
        rows,
        "Per well and fit: images fitted, rms / max residual, hbar*omega "
        "at the minimum, the minimum, the departure from harmonic at the outer turning points, and the "
        "full description.",
    )

    columns = [
        ("state", "well", "%s"),
        ("label", "fit", "%s"),
        ("fitted_images", "images", "%s"),
        ("rms_all_images_meV", "rms (meV)", "%.3f"),
        ("max_abs_residual_meV", "max |res| (meV)", "%.3f"),
        ("hbar_omega_at_min_meV", "hbar w at min (meV)", "%.3f"),
        ("dV_right_turning_point_meV", "dV at right turning pt (meV)", "%.3g"),
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
        npts=npts or rep.params.fd_npts,
        stencil=stencil or rep.params.fd_stencil,
        margin=rep.params.fd_margin,
    )
    return solver, time.perf_counter() - start


def draw_grid_direction(rep, direction):
    """Grid, wavefunctions, a zoom onto individual grid points, and the decay to the walls."""

    plt, ccd = rep.plt, rep.ccd
    solver, seconds = fd_solve(rep, direction)
    info, q = solver.fd_info, solver.fd_q
    fig, axes = plt.subplots(2, 3, figsize=(18, 10), gridspec_kw={"width_ratios": [2.2, 1, 1.3]})
    stats = []
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
        ax_all.plot(q, v_grid, color=color, lw=1.4, label="%s parabola sampled on the grid" % STATE_LABELS[state])
        shown = rep.levels(n)
        for k, v in enumerate(shown):
            gap = energies[shown[k + 1]] - energies[v] if k + 1 < len(shown) else energies[v] - energies[shown[k - 1]]
            amp = 0.45 * gap * psi[v] / np.max(np.abs(psi[v]))
            ax_all.fill_between(q, energies[v], energies[v] + amp, color="dimgray", alpha=0.45, lw=0)
            ax_all.text(q.max(), energies[v], " v=%d" % v, fontsize=7, va="center")
        tp = rep.turning_point(state, n - 1)
        center = ccd_data.center(ccd, state)
        for side in (-1, 1):
            edge = center + side * tp
            outer = q.min() if side < 0 else q.max()
            ax_all.axvspan(min(edge, outer), max(edge, outer), color="gold", alpha=0.12, lw=0)
        ax_all.axvline(q.min(), color="black", lw=2)
        ax_all.axvline(q.max(), color="black", lw=2)
        ax_all.set_ylim(v_grid.min() - 0.3 * spacing, energies[-1] + 1.5 * spacing)
        ax_all.set_title(
            "%s state (%s): %d levels on %d points, walls in black, padding beyond the v = %d turning "
            "points in gold" % (role, STATE_LABELS[state], n, info["npts"], n - 1),
            fontsize=9,
        )
        ax_all.set_ylabel("E (eV)")

        # (b) the highest level at the well centre, where its wavelength is shortest, point by point
        k_max = np.sqrt((energies[-1] - v_grid.min()) / KINETIC)
        wavelength = 2 * np.pi / k_max
        ppw = wavelength / info["h"]
        nodes = np.nonzero(np.diff(np.sign(psi[-1])))[0]
        node = nodes[np.argmin(np.abs(q[nodes] - center))]
        window = np.abs(q - q[node]) < wavelength / 8
        ax_zoom.plot(q[window], psi[-1][window], "-", color=color, lw=0.6, alpha=0.6)
        ax_zoom.plot(q[window], psi[-1][window], "o", color=color, ms=3.2, label="grid points")
        ax_zoom.axhline(0.0, color="gray", lw=0.6)
        ax_zoom.set_title(
            "v = %d at the well centre, a quarter of its shortest wavelength\nh = %.2e: %.0f points "
            "per wavelength (%d shown)" % (n - 1, info["h"], ppw, window.sum()),
            fontsize=9,
        )
        ax_zoom.set_ylabel(r"$\psi_{%d}$ on the grid points" % (n - 1))
        ax_zoom.legend(fontsize=7, edgecolor="black")

        # (c) decay to the walls on a log scale
        for v in sorted({0, n // 2, n - 1}):
            ax_tail.semilogy(q, np.maximum(np.abs(psi[v]) / np.max(np.abs(psi[v])), 1e-20), lw=0.8, label="v = %d" % v)
        ax_tail.axvline(q.min(), color="black", lw=2)
        ax_tail.axvline(q.max(), color="black", lw=2)
        ax_tail.set_ylim(1e-20, 3)
        ax_tail.set_title(
            "|psi| / peak, decay to the walls: %.1e at the walls\n(the ~1e-16 floor is the "
            "eigensolver's float64 precision)" % info["wall"],
            fontsize=9,
        )
        ax_tail.legend(fontsize=7, edgecolor="black")
        for a in axes[row]:
            a.set_xlabel(r"Q (amu$^{1/2}$ $\AA$)")
            plots.style(a)
        stats.append(
            {
                "direction": direction,
                "state": STATE_NAMES[state],
                "role": role,
                "levels": n,
                "npts": info["npts"],
                "h": info["h"],
                "q_lo": q.min(),
                "q_hi": q.max(),
                "points_per_shortest_wavelength": ppw,
                "wall_amplitude": info["wall"],
                "orthonormality": info["orthonormality"],
                "ladder_error_meV": info["ladder_i" if role == "initial" else "ladder_f"],
                "seconds": seconds,
            }
        )
    fig.suptitle(
        "Finite-difference grid, %s (%s): %s, %s, F_ij in %.1f s"
        % (
            direction,
            "ES to GS" if direction == "emission" else "GS to ES",
            info["frame"],
            "%d-point stencil" % info["stencil"],
            seconds,
        ),
        fontsize=13,
    )
    fig.tight_layout()
    rep.out.save_fig(
        fig,
        "02_fd_grid/fd_grid_%s.png" % direction,
        "%s: both wells on the grid with levels at their "
        "eigenvalues, the highest level point by point near its turning point, and |psi| decaying to the "
        "walls." % direction,
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
        ],
    ):
        rep.out.log("  " + line)
    if not convergence:
        return

    # Convergence of F_ij against the Hermite quadrature (accurate to ~1e-12), emission only.
    ref = config_coord_1d_vof_solver(rep.ccd["es_phonon"], rep.ccd["gs_phonon"], rep.ccd["delta_q"])
    ref.compute_franck_condon_integrals_hermite(*rep.params.orders("emission"), npts=rep.params.npts)
    rows = []
    for stencil in sorted(FD_STENCILS):
        for npts in CONVERGENCE_NPTS:
            solver, seconds = fd_solve(rep, "emission", npts=npts, stencil=stencil)
            rows.append(
                {
                    "stencil": stencil,
                    "npts": npts,
                    "h": solver.fd_info["h"],
                    "max_dF_vs_hermite": float(np.max(np.abs(solver.fc_ints - ref.fc_ints))),
                    "ladder_error_meV": max(solver.fd_info["ladder_i"], solver.fd_info["ladder_f"]),
                    "seconds": seconds,
                }
            )
            rep.out.detail(
                "  convergence: stencil %d, %6d points, max dF %.2e, %.1f s"
                % (stencil, npts, rows[-1]["max_dF_vs_hermite"], seconds)
            )
    fig, axes = rep.plt.subplots(1, 3, figsize=(17, 5.2))
    keys = ("max_dF_vs_hermite", "ladder_error_meV", "seconds")
    for stencil, marker in zip(sorted(FD_STENCILS), ("s", "o")):
        sub = [r for r in rows if r["stencil"] == stencil]
        h = np.array([r["h"] for r in sub])
        for ax, key in zip(axes, keys):
            values = np.array([r[key] for r in sub])
            # The measured order is the log-log slope; skip points already at the float64 floor.
            ok = values > (1e-12 if key != "seconds" else 0)
            slope = np.polyfit(np.log(h[ok]), np.log(values[ok]), 1)[0]
            ax.loglog(h, values, marker=marker, label="%d-point stencil: slope %.2f" % (stencil, slope))
            for r in sub:
                r["slope_" + key] = slope
        order = 2 if stencil == 3 else 4
        # The data follow the guides so closely that they would hide them, so draw each 5x higher.
        anchor = 5 * sub[-1]["max_dF_vs_hermite"]
        axes[0].loglog(h, anchor * (h / h[-1]) ** order, color="gray", ls=":", lw=1.0, label="h^%d guide" % order)
    default_h = next(s["h"] for s in stats if s["direction"] == "emission")
    h_ticks = sorted({r["h"] for r in rows}, reverse=True)
    npts_ticks = [next(r["npts"] for r in rows if r["h"] == hv) for hv in h_ticks]
    labels = ("max |F_ij(fd) - F_ij(hermite)|", "max |E_v - E_0 - v hbar w| (meV)", "F_ij time (s)")
    for ax, label in zip(axes, labels):
        ax.axvline(default_h, color="black", ls="--", lw=0.8, label="default: %d points" % rep.params.fd_npts)
        ax.invert_xaxis()
        # A labelled tick at every measured h, and the matching grid size on a top axis.
        ax.set_xticks(h_ticks)
        ax.set_xticklabels(["%.1e" % hv for hv in h_ticks], rotation=30, fontsize=8)
        ax.xaxis.set_minor_locator(matplotlib.ticker.NullLocator())
        top = ax.secondary_xaxis("top")
        top.set_xticks(h_ticks)
        top.set_xticklabels([str(n) for n in npts_ticks], fontsize=8)
        top.xaxis.set_minor_locator(matplotlib.ticker.NullLocator())
        top.set_xlabel("grid points", fontsize=9)
        ax.set_xlabel(r"grid spacing h (amu$^{1/2}$ $\AA$)")
        ax.set_ylabel(label)
        ax.legend(fontsize=8, edgecolor="black")
        plots.style(ax)
        ax.xaxis.set_ticks_position("bottom")
    fig.suptitle(
        "Finite-difference convergence, emission (reference: Hermite quadrature, %d points)" % len(ref.hermite_grid)
    )
    fig.tight_layout()
    rep.out.save_fig(
        fig,
        name + "/fd_convergence.png",
        "Emission: F_ij error against the Hermite quadrature, level "
        "ladder error and time vs grid spacing, for the 3- and 5-point stencils, with O(h^2) / O(h^4) guides.",
    )
    rep.out.save_table(name + "/fd_convergence.csv", rows, "The convergence scan behind fd_convergence.png.")


# ---------------------------------------------------------------------------------- 03


def solve_well(rep, state, spec):
    """FD eigenstates of one well on its own physical grid: ``(q, energies in eV, psi)``."""

    n = rep.order(state)
    center = ccd_data.center(rep.ccd, state)
    span = rep.turning_point(state, n - 1) + rep.params.fd_margin
    q = np.linspace(center - span, center + span, rep.params.fd_npts)
    pot = fits.make_potential(rep.ccd, state, spec)
    energies, psi = config_coord_1d_vof_solver.solve_schrodinger_fd(pot(q), q, n, rep.params.fd_stencil)
    return q, energies, psi


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

        shown = rep.levels(n)
        fig, axes = plt.subplots(len(shown), 2, figsize=(14, 2.1 * len(shown)), squeeze=False)
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
                        "max_abs_dpsi_rel": float(np.max(np.abs(psi[v] - chi[v])) / np.max(np.abs(chi[v]))),
                        "infidelity": float(1.0 - abs(np.sum(psi[v] * chi[v]) * h)),
                    }
                )
            for k, v in enumerate(shown):
                window = np.abs(q - center) < rep.turning_point(state, v) + 0.45
                axes[k, 0].plot(
                    q[window], psi[v][window], color=plots.fit_color(spec), lw=1.0, label="fd, " + fits.label(spec)
                )
                axes[k, 1].plot(
                    q[window],
                    psi[v][window] - chi[v][window],
                    color=plots.fit_color(spec),
                    lw=0.9,
                    label="fd, " + fits.label(spec),
                )
        for k, v in enumerate(shown):
            window = np.abs(q - center) < rep.turning_point(state, v) + 0.45
            axes[k, 0].plot(q[window], chi[v][window], color="#DB4437", lw=2.4, alpha=0.4, label="Hermite (analytic)")
            axes[k, 0].set_ylabel(r"$\chi_{%d}$" % v)
            axes[k, 1].set_ylabel(r"$\psi_{%d} - \chi_{%d}$" % (v, v))
            axes[k, 1].axhline(0.0, color="gray", lw=0.6)
            # The parabola's difference is ~1e-9 and invisible at this scale, so print every fit's maximum.
            maxima = [
                "%s: %.1e" % (spec, np.max(np.abs(rep.wells[(state, spec)][2][v] - chi[v]))) for spec in rep.bounded
            ]
            axes[k, 1].text(
                0.01,
                0.95,
                "max |diff|  " + ",  ".join(maxima),
                transform=axes[k, 1].transAxes,
                fontsize=7,
                va="top",
                family="monospace",
                bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.8},
            )
            for a in axes[k]:
                plots.style(a)
        axes[0, 0].legend(fontsize=7, edgecolor="black", ncol=2)
        axes[0, 0].set_title("%s wavefunctions (pypl has none: its F_ij never forms them)" % STATE_LABELS[state])
        axes[0, 1].set_title("difference from the analytic Hermite function")
        for a in axes[-1]:
            a.set_xlabel(r"Q (amu$^{1/2}$ $\AA$)")
        fig.tight_layout()
        rep.out.save_fig(
            fig,
            name + "/wavefunctions_%s.png" % state,
            "%s: selected levels from the analytic Hermite "
            "function and from finite differences on every bounded fit, and their difference." % STATE_NAMES[state],
        )

    fig, axes = plt.subplots(2, 2, figsize=(14, 8), sharex="col")
    for col, state in enumerate(STATES):
        for spec in rep.bounded:
            sub = [r for r in rows if r["fit"] == spec and r["state"] == STATE_NAMES[state]]
            v = [r["v"] for r in sub]
            c = plots.fit_color(spec)
            axes[0, col].semilogy(
                v, [max(r["infidelity"], 1e-17) for r in sub], marker=".", color=c, label=fits.label(spec)
            )
            axes[1, col].semilogy(
                v, [max(r["max_abs_dpsi_rel"], 1e-17) for r in sub], marker=".", color=c, label=fits.label(spec)
            )
        axes[0, col].set_title("%s: 1 - |<chi_v | psi_v>|" % STATE_LABELS[state])
        axes[1, col].set_title("%s: max |psi_v - chi_v| / max |chi_v|" % STATE_LABELS[state])
        axes[1, col].set_xlabel("vibrational level v")
        for a in axes[:, col]:
            a.legend(fontsize=7, edgecolor="black")
            plots.style(a)
    fig.suptitle(
        "Wavefunction differences from the analytic Hermite functions (parabola = numerical error; "
        "other fits = anharmonicity)",
        fontsize=12,
    )
    fig.tight_layout()
    rep.out.save_fig(
        fig,
        name + "/wavefunction_errors.png",
        "Per level, well and fit: infidelity 1 - |<chi|psi>| and " "max |psi - chi| relative to the Hermite function.",
    )
    rep.out.save_table(
        name + "/wavefunction_errors.csv",
        rows,
        "Per well, fit and level: E_v - E_0, its departure "
        "from v hbar w, max |psi - chi| / max |chi| and the infidelity.",
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
        tested = []
        for m in tests:
            tst = compute_vof(rep.ccd, direction, m, params, log=rep.out.detail)
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
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.6))
    for ax, (direction, (ref, _)) in zip(axes, comparison.items()):
        im = ax.imshow(
            np.log10(np.maximum(np.abs(ref.solver.fc_ints), 1e-18)),
            origin="lower",
            aspect="auto",
            cmap="magma",
            vmin=-12,
            vmax=0,
        )
        ax.set_title("%s: log10 |F_ij| from %s" % (direction, ref.label), fontsize=10)
        ax.set_xlabel("j (final level)")
        ax.set_ylabel("i (initial level)")
        fig.colorbar(im, ax=ax)
    fig.tight_layout()
    rep.out.save_fig(
        fig, "04_fc_matrix/fc_matrix.png", "The reference Franck-Condon matrix, log10 |F_ij|, per direction."
    )
    rep.out.save_fig(
        plots.fc_error_figure(plt, comparison, "F_ij error vs %s (parabola wells)" % ref.label),
        "04_fc_matrix/fc_error_maps.png",
        "log10 |F_ij(method) - F_ij(reference)| per method and direction.",
    )

    fig, axes = plt.subplots(1, 2, figsize=(12, 4.2))
    for ax, (direction, (ref, tested)) in zip(axes, comparison.items()):
        for tst, _ in tested:
            dev = np.max(np.abs(tst.solver.fc_ints - ref.solver.fc_ints), axis=1)
            ax.semilogy(
                np.arange(len(dev)), np.maximum(dev, 1e-18), marker=".", color=plots.result_color(tst), label=tst.label
            )
        ax.set_title("%s: max_j |dF_ij| vs initial level i" % direction, fontsize=10)
        ax.set_xlabel("i (initial level)")
        ax.set_ylabel("max_j |dF_ij|")
        ax.legend(fontsize=8, edgecolor="black")
        plots.style(ax)
    fig.suptitle("Where the F_ij differences sit: they grow with the level (float64 drift of the reference recurrence)")
    fig.tight_layout()
    rep.out.save_fig(fig, "04_fc_matrix/fc_error_by_level.png", "Largest F_ij difference in each row i, per method.")

    title = "VOF on the harmonic wells: every method vs %s" % comparison["emission"][0].label
    rep.out.save_fig(
        plots.vof_comparison_figure(plt, comparison, params, rep.tol_vof, title),
        "05_vof_methods/vof_methods.png",
        "Per direction: every method's VOF overlaid, signed and log error "
        "against the reference, and every method's settings.",
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


def section_vof_fits(rep):
    name = "06_vof_fits"
    rep.out.section(name, SECTIONS[name])
    rep.out.log("\n[%s] fd on every bounded fit vs %s (harmonic)" % (name, rep.reference))
    plt, ccd = rep.plt, rep.ccd
    comparison, rows = {}, []
    for direction in DIRECTIONS:
        ref = compute_vof(ccd, direction, rep.reference, rep.params, log=rep.out.detail)
        tested = []
        for spec in rep.bounded:
            tst = compute_vof(ccd, direction, "fd", dataclasses.replace(rep.params, fd_fit=spec), log=rep.out.detail)
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

    title = "VOF of finite differences on each fit vs %s (differences = anharmonicity of the fit)" % ref.label
    rep.out.save_fig(
        plots.vof_comparison_figure(plt, comparison, rep.params, rep.tol_vof, title, with_settings=False),
        name + "/vof_fits.png",
        "Per direction: fd on every fit overlaid on the harmonic reference, with " "signed and log differences.",
    )

    fig, axes = plt.subplots(1, 2, figsize=(13, 4.6))
    for ax, state in zip(axes, STATES):
        phonon = ccd["%s_phonon" % state]
        for spec in rep.bounded:
            energies = rep.wells[(state, spec)][1]
            v = np.arange(len(energies))
            ax.plot(
                v,
                (energies - energies[0]) * 1e3 - v * phonon,
                marker=".",
                color=plots.fit_color(spec),
                label=fits.label(spec),
            )
        ax.axhline(0.0, color="gray", lw=0.8)
        ax.set_yscale("symlog", linthresh=1.0)
        ax.set_title(
            "%s: E_v - E_0 - v hbar w (hbar w = %.3f meV), the anharmonicity of each fit"
            % (STATE_LABELS[state], phonon),
            fontsize=10,
        )
        ax.set_xlabel("vibrational level v")
        ax.set_ylabel("meV (symlog)")
        ax.legend(fontsize=8, edgecolor="black")
        plots.style(ax)
    fig.tight_layout()
    rep.out.save_fig(
        fig,
        name + "/level_ladders.png",
        "Per well and fit: how far the FD level ladder departs from "
        "the harmonic v hbar w. bulid_fc_lsp still uses the harmonic positions; only F_ij follows the fit.",
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

CUBIC_SPECS = ("cubic", "cubic-near", "cubic-free", "cubic-near-free")
CUBIC_SCALES = (1.0, 0.5, 0.25, 0.125, 0.0625)


def cubic_potential(ccd, state, spec, scale=1.0):
    """The cubic fit ``spec`` of ``state`` as a callable, with its b multiplied by ``scale``."""

    p = fits.make_potential(ccd, state, spec).params
    center = ccd_data.center(ccd, state)
    return lambda q: p["a"] * (q - center) ** 2 + scale * p["b"] * (q - center) ** 3 + p["c"]


def cubic_solve(rep, direction, scale=1.0, orders=None):
    """Perturbative and exact (fd) solvers of ``direction`` on the cubic wells, b scaled by ``scale``."""

    ccd, spec = rep.ccd, rep.params.cubic_fit
    initial, final = direction_states(direction)
    ni, nf = orders or rep.params.orders(direction)
    q_i, q_f = ccd_data.center(ccd, initial), ccd_data.center(ccd, final)
    b = {s: scale * fits.make_potential(ccd, s, spec).params["b"] for s in STATES}
    pt = config_coord_1d_vof_solver(ccd["%s_phonon" % initial], ccd["%s_phonon" % final], ccd["delta_q"])
    with contextlib.redirect_stdout(io.StringIO()):
        pt.compute_franck_condon_integrals_cubic(
            ni, nf, b[initial], b[final], q_i=q_i, q_f=q_f, npts=rep.params.npts, normalize=rep.params.cubic_normalize
        )
    fd = config_coord_1d_vof_solver(ccd["%s_phonon" % initial], ccd["%s_phonon" % final], ccd["delta_q"])
    fd.compute_franck_condon_integrals_fd(
        ni,
        nf,
        cubic_potential(ccd, initial, spec, scale),
        cubic_potential(ccd, final, spec, scale),
        q_i=q_i,
        q_f=q_f,
        npts=rep.params.fd_npts,
        stencil=rep.params.fd_stencil,
        margin=rep.params.fd_margin,
    )
    harmonic = config_coord_1d_vof_solver(ccd["%s_phonon" % initial], ccd["%s_phonon" % final], ccd["delta_q"])
    harmonic.compute_franck_condon_integrals_hermite(ni, nf, npts=rep.params.npts)
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
    fig, axes = plt.subplots(2, 3, figsize=(19, 9.5))
    rows = []
    for row, state in enumerate(STATES):
        q, e = ccd_data.samples(ccd, state)
        center, phonon = ccd_data.center(ccd, state), ccd["%s_phonon" % state]
        q_near, q_grid = np.linspace(q.min() - 0.05, q.max() + 0.05, 600), np.linspace(q_lo, q_hi, 2001)
        harmonic = fits.make_potential(ccd, state, "parabola")
        ax_fit, ax_full, ax_res = axes[row]
        for ax in (ax_fit, ax_full):
            ax.plot(q, e, "o", color="black", ms=5, zorder=5, label="DFT images")
            ax.plot(q_near if ax is ax_fit else q_grid, harmonic(q_near if ax is ax_fit else q_grid), color="gray",
                    lw=2.2, alpha=0.5, label="parabola (harmonic H0) [%s]" % fits.images(ccd, state, "parabola"))  # fmt: skip
        ax_full.axvspan(q.min(), q.max(), color="gold", alpha=0.18, lw=0, label="sampled by DFT")
        for spec in CUBIC_SPECS:
            pot = fits.make_potential(ccd, state, spec)
            p = pot.params
            xi = config_coord_1d_vof_solver.anharmonic_strength(p["b"], phonon)
            mixing = 3 * abs(xi) * np.arange(1, rep.order(state) + 1) ** 1.5
            c = plots.fit_color(spec)
            ls = "-" if spec == rep.params.cubic_fit else "--"
            ax_fit.plot(q_near, pot(q_near), color=c, ls=ls, lw=1.2, label=fits.state_label(ccd, state, spec))
            ax_full.plot(q_grid, pot(q_grid), color=c, ls=ls, lw=1.2, label=fits.state_label(ccd, state, spec))
            ax_res.plot(q, pot.residuals() * 1e3, marker="o", ms=3.5, color=c, ls=ls, lw=0.9, label=fits.label(spec))
            rows.append(
                {
                    "state": STATE_NAMES[state],
                    "fit": spec,
                    "label": fits.label(spec),
                    "images": fits.images(ccd, state, spec),
                    "used_by_method": spec == rep.params.cubic_fit,
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
            ax_full.axvline(wall, color="black", lw=2)
        ax_full.set_ylim(e.min() - 0.5, 1.15 * max(harmonic(q_lo), harmonic(q_hi)))
        ax_res.axhline(0.0, color="gray", lw=0.8)
        ax_fit.set_title("%s: cubic fits near the DFT images (solid = used by the method)" % STATE_LABELS[state],
                         fontsize=10)  # fmt: skip
        ax_full.set_title("%s: over the full fd grid, Q = %.2f..%.2f (walls in black)" % (STATE_LABELS[state], q_lo, q_hi),
                          fontsize=10)  # fmt: skip
        ax_res.set_title("%s: residual at every image" % STATE_LABELS[state], fontsize=10)
        ax_fit.set_ylabel("E (eV)")
        ax_full.set_ylabel("E (eV)")
        ax_res.set_ylabel("fit - DFT (meV)")
        for a in axes[row]:
            a.set_xlabel(r"Q (amu$^{1/2}$ $\AA$)")
            plots.style(a)
        ax_fit.legend(fontsize=7, edgecolor="black")
    fig.suptitle("Cubic fits V = a Q^2 + b Q^3 + c for the perturbative correction (a fixed = the parabola's H0)",
                 fontsize=13)  # fmt: skip
    fig.tight_layout()
    rep.out.save_fig(fig, name + "/cubic_fits.png", "Per well: the four cubic-fit options (a fixed or fitted; all "
                     "images or the parabola's) against the DFT images, over the full fd grid, and their residuals.")  # fmt: skip
    rep.out.save_table(name + "/cubic_fits.csv", rows, "Per well and cubic fit: a, b, turnover, xi, the largest "
                       "first-order mixing 3|xi|(n+1)^1.5, the highest level where it stays <= 0.1, residuals.")  # fmt: skip
    return rows


def cubic_validation_figure(rep, name):
    """b -> 0 scaling against the exact solution, the norms, and where first order breaks down."""

    plt, ccd = rep.plt, rep.ccd
    block = 8
    scan = []
    for scale in CUBIC_SCALES:
        small_pt, small_fd, small_h = cubic_solve(rep, "emission", scale, orders=(12, 12))
        pt, fd, harmonic = cubic_solve(rep, "emission", scale)
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
    zero_pt, _, zero_h = cubic_solve(rep, "emission", 0.0)
    b0_dev = float(np.max(np.abs(zero_pt.fc_ints - zero_h.fc_ints)))
    full_pt, full_fd, _ = cubic_solve(rep, "emission", 1.0)

    fig, axes = plt.subplots(1, 3, figsize=(19, 5.4))
    s = np.array([r["b_scale"] for r in scan])
    ax = axes[0]
    for key, label, marker in (
        ("dF_pt_vs_exact_block", "F_ij (i, j < %d): first-order PT vs exact" % block, "o"),
        ("dF_harmonic_vs_exact_block", "F_ij (i, j < %d): harmonic vs exact" % block, "s"),
        ("dVOF_pt_vs_exact", "VOF / peak: first-order PT vs exact", "^"),
        ("dVOF_harmonic_vs_exact", "VOF / peak: harmonic vs exact", "v"),
    ):
        y = np.array([r[key] for r in scan])
        slope = np.polyfit(np.log(s), np.log(y), 1)[0]
        ax.loglog(s, y, marker=marker, label="%s: slope %.2f" % (label, slope))
    ax.set_xlabel("b scale factor (1 = the fitted b)")
    ax.set_ylabel("max |difference| (exact = fd on the same cubic)")
    ax.set_title("b -> 0 limit: PT error ~ b^2 (slope 2), the effect itself ~ b (slope 1)\nat b = 0: max |F_PT - "
                 "F_hermite| = %.1e" % b0_dev, fontsize=10)  # fmt: skip
    ax.legend(fontsize=7.5, edgecolor="black")
    plots.style(ax)

    ax = axes[1]
    for role, norms, xi, color in (
        ("initial (ES)", full_pt.cubic_norms_i, full_pt.cubic_xi_i, plots.STATE_COLORS["es"]),
        ("final (GS)", full_pt.cubic_norms_f, full_pt.cubic_xi_f, plots.STATE_COLORS["gs"]),
    ):
        n = np.arange(len(norms))
        predicted = (config_coord_1d_vof_solver.cubic_correction_coefficients(len(norms), xi) ** 2).sum(axis=1) - 1
        ax.semilogy(n, norms - 1, "o", ms=3.5, color=color, label="%s: computed <psi~|psi~> - 1" % role)
        ax.semilogy(n, predicted, "-", lw=0.9, color=color, label="%s: sum of squared PT coefficients" % role)
    ax.axhline(0.1, color="gray", ls="--", lw=0.9)
    ax.set_xlabel("vibrational level n")
    ax.set_ylabel("norm - 1")
    ax.set_title("Norms are 1 + O(xi^2), exactly as the coefficients predict\n(emission frame, fitted b)", fontsize=10)
    ax.legend(fontsize=7.5, edgecolor="black")
    plots.style(ax)

    ax = axes[2]
    h = grid_spacing(full_fd.fd_grid)
    for role, xi, psi_fd, color, state in (
        ("initial (ES)", full_pt.cubic_xi_i, full_fd.fd_wavefunctions_i, plots.STATE_COLORS["es"], "es"),
        ("final (GS)", full_pt.cubic_xi_f, full_fd.fd_wavefunctions_f, plots.STATE_COLORS["gs"], "gs"),
    ):
        n_levels = len(psi_fd)
        mixing = 3 * abs(xi) * np.arange(1, n_levels + 1) ** 1.5
        x = full_fd.fd_grid
        center = -ccd["delta_q"] if role.startswith("initial") else 0.0
        coeffs = config_coord_1d_vof_solver.cubic_correction_coefficients(n_levels, xi)
        chi = config_coord_1d_vof_solver.harmonic_wavefunctions(x, n_levels + 2, ccd["%s_phonon" % state], center)
        psi_pt = coeffs @ chi
        psi_pt /= np.sqrt((psi_pt**2).sum(axis=1) * h)[:, None]
        infidelity = 1 - np.abs((psi_pt * psi_fd).sum(axis=1) * h)
        ax.semilogy(np.arange(n_levels), mixing, "-", color=color, lw=1.0, label="%s: 3|xi|(n+1)^1.5" % role)
        ax.semilogy(np.arange(n_levels), np.maximum(infidelity, 1e-16), "o", ms=3, color=color,
                    label="%s: 1 - |<psi~|psi_exact>|" % role)  # fmt: skip
    ax.axhline(0.1, color="gray", ls="--", lw=0.9, label="0.1 (warning threshold)")
    ax.set_xlabel("vibrational level n")
    ax.set_title("Where first order breaks down: mixing grows as n^1.5\n(infidelity vs the exact eigenstates)",
                 fontsize=10)  # fmt: skip
    ax.legend(fontsize=7.5, edgecolor="black")
    plots.style(ax)
    fig.suptitle("Validation of the first-order cubic correction (%s)" % fits.label(rep.params.cubic_fit), fontsize=13)
    fig.tight_layout()
    rep.out.save_fig(fig, name + "/validation.png", "b -> 0 scaling of the PT error (b^2) and of the effect (b) "
                     "against fd on the same cubic; norms vs the predicted 1 + sum c^2; mixing and infidelity per "
                     "level, showing where first order breaks down.")  # fmt: skip
    rep.out.save_table(name + "/validation_b_scaling.csv", scan, "The b-scaling scan behind validation.png.")
    return scan, b0_dev, full_pt


def cubic_wavefunction_figures(rep, name):
    """Harmonic, first-order corrected and exact (fd on the cubic) wavefunctions, per well."""

    plt, ccd = rep.plt, rep.ccd
    for state in STATES:
        n = rep.order(state)
        phonon, center = ccd["%s_phonon" % state], ccd_data.center(ccd, state)
        q, _, psi_fd = solve_well(rep, state, rep.params.cubic_fit)
        xi = config_coord_1d_vof_solver.anharmonic_strength(fits.make_potential(ccd, state, rep.params.cubic_fit).params["b"],
                                                            phonon)  # fmt: skip
        chi = config_coord_1d_vof_solver.harmonic_wavefunctions(q, n + 2, phonon, center)
        psi_pt = config_coord_1d_vof_solver.cubic_correction_coefficients(n, xi) @ chi
        shown = rep.levels(n)
        fig, axes = plt.subplots(len(shown), 2, figsize=(14, 2.1 * len(shown)), squeeze=False)
        for k, v in enumerate(shown):
            window = np.abs(q - center) < rep.turning_point(state, v) + 0.45
            mix = 3 * abs(xi) * (v + 1) ** 1.5
            axes[k, 0].plot(q[window], chi[v][window], color="gray", lw=2.4, alpha=0.5, label="harmonic")
            axes[k, 0].plot(q[window], psi_fd[v][window], color="black", lw=1.0, label="exact (fd on the cubic)")
            axes[k, 0].plot(q[window], psi_pt[v][window], color=plots.METHOD_COLORS["cubic"], ls="--", lw=1.0,
                            label="first-order PT")  # fmt: skip
            axes[k, 1].plot(q[window], psi_pt[v][window] - psi_fd[v][window], color=plots.METHOD_COLORS["cubic"],
                            lw=0.9, label="PT - exact")  # fmt: skip
            axes[k, 1].plot(
                q[window], chi[v][window] - psi_fd[v][window], color="gray", lw=0.9, label="harmonic - exact"
            )
            axes[k, 1].axhline(0.0, color="gray", lw=0.6)
            axes[k, 1].text(0.01, 0.95, "mixing 3|xi|(n+1)^1.5 = %.3g;  max |PT - exact| %.1e, |harmonic - exact| %.1e"
                            % (mix, np.max(np.abs(psi_pt[v] - psi_fd[v])), np.max(np.abs(chi[v] - psi_fd[v]))),
                            transform=axes[k, 1].transAxes, fontsize=7, va="top", family="monospace",
                            bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.8})  # fmt: skip
            axes[k, 0].set_ylabel(r"$\psi_{%d}$" % v)
            for a in axes[k]:
                plots.style(a)
        axes[0, 0].legend(fontsize=7, edgecolor="black", ncol=3)
        axes[0, 1].legend(fontsize=7, edgecolor="black", ncol=2, loc="lower right")
        axes[0, 0].set_title("%s wavefunctions, xi = %+.2e (%s)" % (STATE_LABELS[state], xi,
                             fits.label(rep.params.cubic_fit)))  # fmt: skip
        axes[0, 1].set_title("difference from the exact eigenstates of the same cubic")
        for a in axes[-1]:
            a.set_xlabel(r"Q (amu$^{1/2}$ $\AA$)")
        fig.tight_layout()
        rep.out.save_fig(fig, name + "/wavefunctions_%s.png" % state, "%s: harmonic, first-order corrected and exact "
                         "(fd on the same cubic) wavefunctions for selected levels, and the differences from exact."
                         % STATE_NAMES[state])  # fmt: skip


def section_cubic(rep):
    """The first-order cubic anharmonic correction (method ``cubic``): fits, validation, wavefunctions, F_ij, VOF."""

    name = "07_cubic_perturbation"
    rep.out.section(name, SECTIONS[name])
    rep.out.log("\n[%s] first-order cubic correction, %s" % (name, fits.label(rep.params.cubic_fit)))
    fit_rows = cubic_fit_figure(rep, name)
    scan, b0_dev, _ = cubic_validation_figure(rep, name)
    cubic_wavefunction_figures(rep, name)

    exact_params = dataclasses.replace(rep.params, fd_fit=rep.params.cubic_fit)
    comparison, rows = {}, []
    for direction in DIRECTIONS:
        exact = compute_vof(rep.ccd, direction, "fd", exact_params, log=rep.out.detail)
        exact.label = "exact: fd on the cubic"
        tested = []
        for method, params in (("cubic", rep.params), (rep.reference, rep.params)):
            tst = compute_vof(rep.ccd, direction, method, params, log=rep.out.detail)
            stats = compare_results(exact, tst, rep.tol_fc, rep.tol_vof)
            tested.append((tst, stats))
            rows.append(stats)
        comparison[direction] = (exact, tested)
        columns = {"E_minus_EZPL_meV": exact.axis, "VOF exact (fd on the cubic)": exact.lineshape}
        columns.update({"VOF %s" % t.label: t.lineshape for t, _ in tested})
        rep.out.save_csv(name + "/vof_%s.csv" % direction, columns, "%s: VOF of the exact cubic solution, the "
                         "first-order correction and the harmonic reference." % direction)  # fmt: skip
    title = "VOF on the cubic wells: first-order PT and harmonic %s vs the exact (fd) solution" % rep.reference
    rep.out.save_fig(plots.vof_comparison_figure(rep.plt, comparison, rep.params, rep.tol_vof, title),
                     name + "/vof_cubic.png", "Per direction: the exact VOF on the cubic wells (fd), the first-order "
                     "correction and the harmonic reference, with signed / log differences and every setting.")  # fmt: skip
    rep.out.save_fig(plots.fc_error_figure(rep.plt, comparison, "F_ij error vs the exact cubic solution"),
                     name + "/fc_error_maps.png", "log10 |F_ij - F_ij(exact cubic)| of the first-order correction and "
                     "of the harmonic reference.")  # fmt: skip
    rep.out.save_table(name + "/summary.csv", [{k: r[k] for k in CSV_KEYS} for r in rows], "F_ij and VOF differences "
                       "from the exact cubic solution, per method and direction.")  # fmt: skip

    used = [r for r in fit_rows if r["used_by_method"]]
    columns = [("state", "well", "%s"), ("label", "cubic fit", "%s"), ("b", "b (eV/amu^1.5 A^3)", "%+.5f"),
               ("xi", "xi", "%+.3e"), ("max_mixing", "max mixing", "%.3g"), ("n_valid", "n_valid", "%d"),
               ("turnover_q", "turnover Q", "%.1f")]  # fmt: skip
    for line in text_table(fit_rows, columns):
        rep.out.log("  " + line)
    vof_columns = [("direction", "direction", "%s"), ("test", "vs exact (fd on the cubic)", "%s"),
                   ("fc_max", "max dF_ij", "%.2e"), ("fc_low", "max dF_ij (i,j<30)", "%.2e"),
                   ("vof_max", "max dVOF / peak", "%.2e")]  # fmt: skip
    for line in text_table(rows, vof_columns):
        rep.out.log("  " + line)
    slopes = {}
    s = np.log([r["b_scale"] for r in scan])
    for key in ("dF_pt_vs_exact_block", "dVOF_pt_vs_exact", "dF_harmonic_vs_exact_block"):
        slopes[key] = float(np.polyfit(s, np.log([r[key] for r in scan]), 1)[0])
    rep.out.log("  b -> 0: PT error slope %.2f (F_ij, i,j<8) / %.2f (VOF), effect slope %.2f; at b = 0 |F_PT - "
                "F_hermite| = %.1e" % (slopes["dF_pt_vs_exact_block"], slopes["dVOF_pt_vs_exact"],
                                       slopes["dF_harmonic_vs_exact_block"], b0_dev))  # fmt: skip
    rep.headline += ["## First-order cubic correction (%s)" % fits.label(rep.params.cubic_fit), ""]
    rep.headline += markdown_table(used, columns) + [""] + markdown_table(rows, vof_columns) + [""]
    rep.headline += [
        "b -> 0 validation: the PT error falls as b^%.2f (F_ij, i, j < 8) and b^%.2f (VOF), the anharmonic effect as "
        "b^%.2f; at b = 0 the method equals the Hermite route to %.1e. The mixing 3|xi|(n+1)^1.5 exceeds 0.1 above "
        "n = %s, so the high-level F_ij are not perturbative; the VOF is, because the intensity sits in low levels."
        % (slopes["dF_pt_vs_exact_block"], slopes["dVOF_pt_vs_exact"], slopes["dF_harmonic_vs_exact_block"], b0_dev,
           " / ".join("%d (%s)" % (r["n_valid"], r["state"]) for r in used)),
        "",
    ]  # fmt: skip


# -------------------------------------------------------------------------------- main


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
    add_vof_args(p)
    p.add_argument("--tol-fc", type=float, default=1e-7, help="max |dF_ij| for OK in 05 (default 1e-7)")
    p.add_argument("--tol-vof", type=float, default=1e-8, help="max |dVOF| / peak for OK in 05 (default 1e-8)")
    p.add_argument("--tag", default="", help="name of the run directory (default: 'default')")
    p.add_argument("-v", "--verbose", action="store_true", help="print every method's settings, not only summaries")
    return p.parse_args()


def main():
    args = parse_args()
    params = params_from_args(args)
    try:
        args.fits = list(dict.fromkeys(fits.canonical(spec) for spec in args.fits))
    except ValueError as err:
        raise SystemExit("Invalid fit: %s" % err)
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

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
            "params": dataclasses.asdict(params),
        },
    )
    rep = Report(ccd, params, out, plt, args.fits, bounded, args.reference, args.tol_fc, args.tol_vof)

    out.log(
        "1D CCD VOF report: delta_q = %.6f amu^0.5 A, GS %.5f meV, ES %.5f meV"
        % (ccd["delta_q"], ccd["gs_phonon"], ccd["es_phonon"])
    )
    for line in params.describe():
        out.log("  " + line)
    if skipped:
        out.log("  unbounded fits (01 only, no bound states for 02-06): %s" % ", ".join(skipped))

    section_potential_fits(rep)
    section_fd_grid(rep, convergence=not args.skip_convergence)
    section_wavefunctions(rep)
    methods_ok = section_methods(rep)
    section_vof_fits(rep)
    cubic_skipped = None
    if not args.skip_cubic:
        # Section 07 validates first-order PT against the exact (fd) solution of the same cubic, which exists only
        # if the cubic is a bound well on the grid (its turnover outside it).
        try:
            fits.check_bounded(ccd, params.cubic_fit, fd_q_range(ccd, params))
        except ValueError as err:
            cubic_skipped = "the %s fit has no exact solution to compare with: %s" % (params.cubic_fit, err)
            out.log("\n07_cubic_perturbation skipped: %s" % cubic_skipped)
            out.save_json(
                "07_cubic_perturbation/skipped.json",
                {"cubic_fit": params.cubic_fit, "reason": cubic_skipped},
                "Why this section was skipped (the cubic fit is unbounded on the fd grid).",
            )
        else:
            section_cubic(rep)

    verdict = "PASS" if methods_ok else "FAIL"
    if cubic_skipped:
        out.add_summary("**Section 07 skipped:** %s." % cubic_skipped, "")
    out.add_summary(
        "CCD: delta_q = %.6f amu^0.5 A, GS hbar w = %.5f meV, ES hbar w = %.5f meV. %s."
        % (ccd["delta_q"], ccd["gs_phonon"], ccd["es_phonon"], "; ".join(params.describe())),
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
