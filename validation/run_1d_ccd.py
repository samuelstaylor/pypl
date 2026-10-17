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
from harness.plots import STATE_COLORS, style
from pypl.config_coord_1d_solver import config_coord_1d_solver

BLUE, RED = STATE_COLORS["es"], STATE_COLORS["gs"]


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
    return p.parse_args()


def compute_spectrum(omega_i, omega_f, delta_q, ni, nf, args, eneaxis, spectrum_type):
    """Run the full FC -> lineshape -> spectrum chain for one direction of the transition."""

    solver = config_coord_1d_solver(omega_i, omega_f, delta_q)
    solver.compute_franck_condon_integrals(ni=ni, nf=nf)
    solver.bulid_fc_lsp(eneaxis=eneaxis, temp=args.temp, sigma=args.sigma, zpl_lorentzian=True, gamma=args.gamma)
    return solver.compute_spectrum(tdm=args.tdm, zpl=args.ezpl, spectrum_type=spectrum_type)


def draw_ccd(plt, ccd):
    fig, ax = plt.subplots(figsize=(5, 7))
    for state, color, marker in (("es", BLUE, "s"), ("gs", RED, "o")):
        name = ccd_data.STATE_NAMES[state]
        q, e = ccd_data.samples(ccd, state)
        ax.plot(q, e, color=color, ls="", marker=marker, ms=5, label=ccd_data.energy_label(state))
        ax.plot(q, ccd["%s_fit_energies" % state], color=color, ls="--", lw=1.0, label="%s fit" % name)
        ax.text(
            0.46,
            0.6 if state == "es" else 0.4,
            "%s phonon: %.2f meV\n%s HRF: %.2f" % (name, ccd["%s_phonon" % state], name, ccd["%s_hrf" % state]),
            color=color,
            transform=ax.transAxes,
        )
    ax.axvline(0.0, color="gray", ls="--", lw=1)
    ax.axvline(ccd["delta_q"], color="gray", ls="--", lw=1)
    ax.legend(loc="center left", edgecolor="black")
    ax.set_xlim((-0.2, 0.9))
    ax.set_xlabel(r"Q (amu$^{1/2}$ $\AA$)")
    ax.set_ylabel("Total energy (eV)")
    style(ax, grid=False)
    return fig


def draw_spectra(plt, pl, pl_norm, ab, ab_norm, temp):
    fig, ax = plt.subplots(figsize=(12, 4))
    ax.plot(pl[0] * 1e-3, pl_norm * 1e3, color=RED, lw=1, label="PL")
    ax.plot(ab[0] * 1e-3, ab_norm * 1e3, color=BLUE, lw=1, label="Abs")
    ax.set_xlim((1.5, 2.4))
    ax.set_ylim((0.0, 6))
    ax.legend(loc="upper right", edgecolor="black")
    ax.set_xlabel(r"$\hbar\omega$ (eV)")
    ax.set_ylabel("PL / Abs (arb. unit)")
    ax.set_title("NV$^-$ in diamond, 1D CCD, T = %g K" % temp)
    style(ax)
    return fig


def main():
    args = parse_args()
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
