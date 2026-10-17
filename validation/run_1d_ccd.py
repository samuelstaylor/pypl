#!/usr/bin/env python3
"""End-to-end 1D CCD run: NV center in diamond, PL and absorption spectra.

Reproduces doc/tutorials/pypl_002.ipynb as a plain script with PyPL's own, unmodified
``config_coord_1d_solver``, and writes every result to ``validation/outputs/run_1d_ccd/<tag or default>/``:

    data/ccd_summary.csv      delta_q, frequencies, Huang-Rhys factors and the run parameters
    data/pl_spectrum.csv      PL spectrum (normalized)
    data/abs_spectrum.csv     absorption spectrum (normalized)
    figures/ccd_diagram.png   the computed images with their parabola fits
    figures/spectra.png       PL and absorption

    python validation/run_1d_ccd.py                       # defaults = tutorial 002
    python validation/run_1d_ccd.py --temp 300 --tag 300K
    python validation/run_1d_ccd.py --order-es 80 --order-gs 90 --tag hiorder
    python validation/run_1d_ccd.py --show                # also open the plots
"""

import argparse

import matplotlib
import numpy as np

from harness import ccd_data
from harness.output import OutputDir
from harness import plots
from harness.plots import ENERGY_LABEL, FIGURE, Q_LABEL, Q_UNIT, STATE_COLORS, legend_below, off_scale_note, style
from pypl.config_coord_1d_solver import config_coord_1d_solver

BLUE, RED = "#4285F4", "#DB4437"  # upper (ES) and lower (GS) state, as in tutorial 002


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--temp", type=float, default=5.0, help="temperature in K (default: 5)")
    p.add_argument("--ezpl", type=float, default=1945.0, help="ZPL energy in meV (default: 1945)")
    p.add_argument("--order-es", type=int, default=50, help="number of ES vibrational states (default: 50)")
    p.add_argument("--order-gs", type=int, default=60, help="number of GS vibrational states (default: 60)")
    p.add_argument(
        "--sigma",
        type=float,
        nargs=2,
        default=[8.0, 25.0],
        metavar=("SIGMA0", "SIGMA1"),
        help="Gaussian broadening range in meV (default: 8 25)",
    )
    p.add_argument("--gamma", type=float, default=0.3, help="ZPL Lorentzian broadening in meV (default: 0.3)")
    p.add_argument(
        "--ene-range",
        type=float,
        nargs=2,
        default=[-200.0, 1200.0],
        metavar=("LO", "HI"),
        help="lineshape energy range in meV (default: -200 1200)",
    )
    p.add_argument("--resol", type=int, default=1401, help="number of energy points (default: 1401)")
    p.add_argument("--tdm", type=float, default=1.0, help="transition dipole moment (default: 1.0)")
    p.add_argument("--tag", default="", help="name of the run directory (default: 'default')")
    p.add_argument("--show", action="store_true", help="display the figures in addition to saving them")
    ccd_data.add_dataset_arg(p)
    plots.add_figure_args(p)
    return p.parse_args()


def compute_spectrum(omega_i, omega_f, delta_q, ni, nf, args, eneaxis, spectrum_type):
    """Run the full FC -> lineshape -> spectrum chain for one direction of the transition."""

    solver = config_coord_1d_solver(omega_i, omega_f, delta_q)
    solver.compute_franck_condon_integrals(ni=ni, nf=nf)
    solver.bulid_fc_lsp(eneaxis=eneaxis, temp=args.temp, sigma=args.sigma, zpl_lorentzian=True, gamma=args.gamma)
    return solver.compute_spectrum(tdm=args.tdm, zpl=args.ezpl, spectrum_type=spectrum_type)


def draw_ccd(plt, ccd):
    """The two wells on a broken energy axis (upper state on top), with the images and the parabola fits."""

    fig, (top, bottom) = plt.subplots(
        2, 1, figsize=(4.8, 5.4), sharex=True, gridspec_kw={"hspace": 0.06, "height_ratios": (1, 1)}
    )
    q_all = np.concatenate([ccd_data.samples(ccd, s)[0] for s in ccd_data.STATES])
    x = np.linspace(q_all.min() - 0.1, q_all.max() + 0.1, 400)
    for ax, state, color, marker in ((top, "es", BLUE, "s"), (bottom, "gs", RED, "o")):
        q, e = ccd_data.samples(ccd, state)
        pot = ccd["%s_potential" % state]
        used = pot.index
        ax.plot(
            x,
            pot(x),
            color=color,
            lw=1.3,
            label="%s parabola (images %d-%d)" % (ccd_data.STATE_NAMES[state], used[0] + 1, used[-1] + 1),
        )
        ax.plot(q, e, color=color, ls="", marker=marker, ms=5, label=ccd_data.energy_label(state))
        ax.text(
            0.97 if state == "es" else 0.03,
            0.92,
            r"%s: $\hbar\omega$ = %.2f meV, $S$ = %.2f"
            % (ccd_data.STATE_NAMES[state], ccd["%s_phonon" % state], ccd["%s_hrf" % state]),
            color=color,
            transform=ax.transAxes,
            ha="right" if state == "es" else "left",
            va="top",
        )
        e_min = ccd["%s_params" % state][1]
        span = max(e.max(), pot(x).min()) - e_min
        ax.set_ylim(e_min - 0.08 * span, e_min + 1.25 * span)
        for xv in (0.0, ccd["delta_q"]):
            ax.axvline(xv, color="0.6", ls="--", lw=0.8, zorder=0)
        style(ax, grid=False)
    # Broken-axis marks.
    for ax, y in ((top, 0.0), (bottom, 1.0)):
        kw = dict(transform=ax.transAxes, color="k", clip_on=False, lw=0.8)
        ax.plot((-0.015, 0.015), (y - 0.015, y + 0.015), **kw)
        ax.plot((0.985, 1.015), (y - 0.015, y + 0.015), **kw)
    top.spines["bottom"].set_visible(False)
    bottom.spines["top"].set_visible(False)
    top.tick_params(bottom=False, labelbottom=False, labeltop=False)
    bottom.tick_params(top=False)
    top.tick_params(labeltop=False)
    bottom.set_xlabel(Q_LABEL)
    fig.supylabel(ENERGY_LABEL, fontsize=FIGURE["font_size"], x=-0.03)
    top.set_title(r"%s: 1D CCD, $\Delta Q$ = %.4f %s" % (ccd_data.active().title, ccd["delta_q"], Q_UNIT))
    handles, labels = [], []
    for ax in (top, bottom):
        h, l = ax.get_legend_handles_labels()
        handles += h
        labels += l
    legend_below(fig, handles, labels, ncol=1, y=0.0)
    return fig


def draw_spectra(plt, pl, pl_norm, ab, ab_norm, temp):
    fig, ax = plt.subplots(figsize=(7.2, 3.2))
    ax.plot(pl[0] * 1e-3, pl_norm * 1e3, color=RED, lw=1.2, label="Photoluminescence")
    ax.plot(ab[0] * 1e-3, ab_norm * 1e3, color=BLUE, lw=1.2, label="Absorption")
    ax.set_xlim((1.5, 2.4))
    window = (pl[0] * 1e-3 > 1.5) & (pl[0] * 1e-3 < 1.9)
    sideband = float(np.max(pl_norm[window] * 1e3))
    ax.set_ylim((0.0, 1.35 * sideband))
    zpl = max(float(np.max(pl_norm)), float(np.max(ab_norm))) * 1e3
    if zpl > 1.35 * sideband:
        off_scale_note(ax, zpl, x=pl[0][np.argmax(pl_norm)] * 1e-3)
    ax.legend(loc="upper right")
    ax.set_xlabel(r"$\hbar\omega$ (eV)")
    ax.set_ylabel(r"Normalized intensity (eV$^{-1}$)")
    ax.set_title(r"%s: 1D CCD spectra, $T$ = %g K" % (ccd_data.active().title, temp))
    style(ax)
    return fig


def main():
    args = parse_args()
    ccd_data.use_dataset(args.dataset)
    plots.configure(args)
    if not args.show:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out = OutputDir("run_1d_ccd", args.tag or "default", title="Tutorial 002 as a script", keep_figures=args.show)
    out.section("data", "Spectra and CCD numbers, as CSV.")
    out.section("figures", "The CCD diagram and the spectra.")
    out.save_json("settings.json", vars(args))

    out.log("[1/3] Constructing the 1D configurational coordinate diagram")
    ccd = ccd_data.build_ccd()
    delta_q = ccd["delta_q"]

    out.log("")
    out.log("[2/3] Computing spectra")
    out.log(
        "  temp=%g K  ezpl=%g meV  sigma=%s meV  gamma=%g meV  orders=(es %d, gs %d)"
        % (args.temp, args.ezpl, args.sigma, args.gamma, args.order_es, args.order_gs)
    )
    eneaxis = np.linspace(args.ene_range[0], args.ene_range[1], args.resol)
    # PL: the excited state is the initial state. Absorption swaps the roles.
    print("  PL  (ES -> GS):", end=" ")
    pl = compute_spectrum(
        ccd["es_phonon"], ccd["gs_phonon"], delta_q, args.order_es, args.order_gs, args, eneaxis, "PL"
    )
    print("  Abs (GS -> ES):", end=" ")
    ab = compute_spectrum(
        ccd["gs_phonon"], ccd["es_phonon"], delta_q, args.order_gs, args.order_es, args, eneaxis, "Abs"
    )
    pl_norm, ab_norm = ccd_data.normalize(*pl), ccd_data.normalize(*ab)
    pl_peak, ab_peak = pl[0][np.argmax(pl[1])], ab[0][np.argmax(ab[1])]
    out.log("  PL  peak at % .2f meV" % pl_peak)
    out.log("  Abs peak at % .2f meV" % ab_peak)

    out.log("")
    out.log("[3/3] Writing outputs")
    summary = [{"quantity": k, "value": ccd[k]} for k in ("delta_q", "gs_phonon", "es_phonon", "gs_hrf", "es_hrf")]
    summary += [
        {"quantity": k, "value": getattr(args, k)}
        for k in ("temp", "ezpl", "gamma", "tdm", "order_es", "order_gs", "resol")
    ]
    summary.append({"quantity": "sigma", "value": " ".join(str(v) for v in args.sigma)})
    out.save_table("data/ccd_summary.csv", summary, "delta_q (amu^0.5 A), frequencies (meV), HRFs, run parameters.")
    out.save_csv("data/pl_spectrum.csv", {"energy_meV": pl[0], "PL_normalized_per_meV": pl_norm}, "PL spectrum.")
    out.save_csv("data/abs_spectrum.csv", {"energy_meV": ab[0], "Abs_normalized_per_meV": ab_norm}, "Absorption.")
    out.save_fig(
        draw_ccd(plt, ccd),
        "figures/ccd_diagram.png",
        "Computed images (%s), parabola fits, frequencies, HRFs." % "; ".join(ccd_data.ENERGY_METHODS.values()),
    )
    out.save_fig(draw_spectra(plt, pl, pl_norm, ab, ab_norm, args.temp), "figures/spectra.png", "PL and absorption.")
    out.add_summary("PL peak at %.2f meV, absorption peak at %.2f meV." % (pl_peak, ab_peak))
    out.close()
    if args.show:
        plt.show()


if __name__ == "__main__":
    main()
