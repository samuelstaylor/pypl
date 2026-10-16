#!/usr/bin/env python3
"""Quick PASS/FAIL comparison of the Franck-Condon lineshape (VOF) from different F_ij methods.

    python testing/compare_vof.py                         # PyPL (Ruhoff) vs every other method
    python testing/compare_vof.py --test fd               # PyPL vs finite differences only
    python testing/compare_vof.py --reference hermite --test fd pypl
    python testing/compare_vof.py --test fd --fd-fit spline   # anharmonic wells: FAIL by design
    python testing/compare_vof.py --test cubic                # cubic perturbation (anharmonic: FAIL by design)
    python testing/compare_vof.py --reference fd --fd-fit cubic --test cubic   # perturbation vs exact on the cubic
    python testing/compare_vof.py --temp 300 --tag 300K -v

Every lineshape comes from the same ``bulid_fc_lsp``; the methods (``harness.methods.METHODS``)
differ only in where F_ij comes from, so any disagreement is a disagreement in F_ij. For the
thorough, sectioned comparison (fits, grid, wavefunctions, F_ij, VOF), run ``vof_report.py``.

Output, in ``testing/outputs/compare_vof/<reference>_vs_<tests>[_fit-<fit>][_<tag>]/``:

    figures/vof_comparison.png   per direction: VOFs overlaid, signed error, log error, settings
    figures/fc_error_maps.png    log10 |F_ij(test) - F_ij(reference)|
    data/vof_<direction>.csv     axis, every lineshape and every difference
    data/summary.csv             every error statistic, OK / FAIL
    log.txt, settings.json, index.md

Exit status is 0 only if every test agrees with the reference within tolerance.
"""

import argparse
import dataclasses

import matplotlib
from harness import ccd_data, fits, plots
from harness.methods import (
    COMPARISON_COLUMNS,
    CSV_KEYS,
    DIRECTION_LABELS,
    DIRECTIONS,
    METHODS,
    add_vof_args,
    compare_results,
    compute_vof,
    is_anharmonic,
    method_label,
    params_from_args,
)
from harness.output import OutputDir, markdown_table, text_table


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    names = ", ".join("%s = %s" % (m.name, m.label) for m in METHODS.values())
    p.add_argument("--reference", choices=sorted(METHODS), default="pypl", help="reference method: " + names)
    p.add_argument(
        "--test",
        nargs="+",
        choices=sorted(METHODS),
        metavar="METHOD",
        help="methods compared with the reference (default: every other method)",
    )
    add_vof_args(p)
    g = p.add_argument_group("tolerances")
    g.add_argument(
        "--tol-fc", type=float, default=1e-7, help="max |dF_ij| (default 1e-7; the solver's own float64 drift is ~2e-8)"
    )
    g.add_argument("--tol-vof", type=float, default=1e-8, help="max |dVOF| / peak (default 1e-8)")
    p.add_argument("--tag", default="", help="suffix of the run directory")
    p.add_argument("-v", "--verbose", action="store_true", help="print every method's settings, not only the summary")
    p.add_argument("--show", action="store_true", help="display the figures as well as saving them")
    return p.parse_args()


def run_direction(ccd, direction, args, tests, params, out):
    """Reference and test VOFs for one direction: ``(ref, [(test, stats), ...])``; writes the CSV."""

    out.detail("")
    out.detail("=" * 100)
    ref = compute_vof(ccd, direction, args.reference, params, log=out.detail)
    tested = []
    for name in tests:
        tst = compute_vof(ccd, direction, name, params, log=out.detail)
        tested.append((tst, compare_results(ref, tst, args.tol_fc, args.tol_vof)))

    columns = {"E_minus_EZPL_meV": ref.axis, "VOF %s" % ref.label: ref.lineshape}
    columns.update({"VOF %s" % t.label: t.lineshape for t, _ in tested})
    columns.update({"%s - %s" % (t.label, ref.label): t.lineshape - ref.lineshape for t, _ in tested})
    description = "%s (%s): every lineshape and every difference from the reference."
    out.save_csv("data/vof_%s.csv" % direction, columns, description % (direction, DIRECTION_LABELS[direction]))
    return ref, tested


def main():
    args = parse_args()
    params = params_from_args(args)
    # By default only the methods that solve the same harmonic problem; anharmonic ones on request.
    default = [m for m in METHODS if not METHODS[m].anharmonic]
    tests = [t for t in dict.fromkeys(args.test or default) if t != args.reference]
    if not tests:
        raise SystemExit("--test names only the reference method; there is nothing to compare")
    if not args.show:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    anharmonic = [method_label(m, params) for m in [args.reference] + tests if is_anharmonic(m, params)]
    run = "%s_vs_%s" % (args.reference, "-".join(tests))
    if "fd" in [args.reference] + tests and params.fd_fit != fits.DEFAULT_FIT:
        run += "_fit-%s" % params.fd_fit
    if "cubic" in [args.reference] + tests:
        run += "_cubicfit-%s%s" % (params.cubic_fit, "-normalized" if params.cubic_normalize else "")
    run += "_%s" % args.tag if args.tag else ""
    out = OutputDir(
        "compare_vof", run, title="VOF method comparison: %s" % run, verbose=args.verbose, keep_figures=args.show
    )
    out.section("figures", "Comparison figures.")
    out.section("data", "Numbers behind the figures, as CSV.")
    settings = {"reference": args.reference, "tests": tests, "tol_fc": args.tol_fc, "tol_vof": args.tol_vof}
    out.save_json("settings.json", dict(settings, params=dataclasses.asdict(params)))

    ccd = ccd_data.build_ccd(verbose=False)
    ref_label = method_label(args.reference, params)
    out.log("Comparing VOFs: %s (reference) vs %s" % (ref_label, ", ".join(method_label(t, params) for t in tests)))
    for line in params.describe():
        out.log("  " + line)

    comparison = {d: run_direction(ccd, d, args, tests, params, out) for d in DIRECTIONS}
    rows = [stats for _, tested in comparison.values() for _, stats in tested]
    all_ok = all(r["ok"] for r in rows)
    verdict = "PASS" if all_ok else "FAIL"
    if {args.reference, *tests} == {"fd", "cubic"} and params.fd_fit == params.cubic_fit:
        verdict += "; same cubic potential: differences are the error of first-order perturbation theory"
    elif anharmonic:
        verdict += "; %s anharmonic, so differences are physics, not error" % " and ".join(anharmonic)

    title = "Franck-Condon lineshape via bulid_fc_lsp: %s (reference)   [%s]" % (ref_label, verdict)
    out.save_fig(
        plots.vof_comparison_figure(plt, comparison, params, args.tol_vof, title),
        "figures/vof_comparison.png",
        "Per direction: VOFs overlaid, (test - ref) / peak, log10 |test - ref| / peak with the tolerance, "
        "and every method's settings.",
    )
    out.save_fig(
        plots.fc_error_figure(plt, comparison, "Franck-Condon matrix error vs %s" % ref_label),
        "figures/fc_error_maps.png",
        "log10 |F_ij(test) - F_ij(reference)| per test and direction.",
    )
    out.save_table(
        "data/summary.csv",
        [{k: r[k] for k in CSV_KEYS} for r in rows],
        "Every error statistic per direction and test, with OK / FAIL.",
    )

    out.log("")
    for line in text_table(rows, COMPARISON_COLUMNS):
        out.log(line)
    out.log("")
    out.log("PASS: every lineshape matches the reference." if all_ok else "FAIL: the lineshapes differ.")
    methods = {args.reference, *tests}
    same_cubic = methods == {"fd", "cubic"} and params.fd_fit == params.cubic_fit
    if same_cubic:
        out.log(
            "NOTE: fd and cubic solve the same cubic potential (%s), fd exactly and cubic to first order in b, so "
            "these differences are the truncation error of the perturbation theory." % fits.label(params.cubic_fit)
        )
    elif anharmonic:
        out.log(
            "NOTE: %s solve(s) an anharmonic problem, not the harmonic wells pypl / hermite use; differences from "
            "them are the anharmonicity, not an error. For a like-for-like test use the default methods and "
            "--fd-fit parabola." % ", ".join(anharmonic)
        )
    tolerances = "max dF_ij <= %.0e, max dVOF / peak <= %.0e" % (args.tol_fc, args.tol_vof)
    out.add_summary("**%s** (tolerances: %s)" % (verdict, tolerances), "")
    out.add_summary(*markdown_table(rows, COMPARISON_COLUMNS))
    out.close()
    if args.show:
        plt.show()
    raise SystemExit(0 if all_ok else 1)


if __name__ == "__main__":
    main()
