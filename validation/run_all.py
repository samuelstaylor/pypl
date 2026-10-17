#!/usr/bin/env python3
"""Master test suite: run every VOF method and fit with one set of parameters, and compare them all.

    python validation/run_all.py                                   # validation/run_all.toml, every step (~15 min)
    python validation/run_all.py --config my_params.toml           # your own copy of run_all.toml
    python validation/run_all.py --set physics.temp=300 --name 300K
    python validation/run_all.py --set fd.npts=40001 --set fd.stencil=3 --name fd-3pt
    python validation/run_all.py --set 'fd.fits=["spline4", "morse"]' --only compare plots
    python validation/run_all.py --dry-run                         # list every command, run nothing
    python validation/run_all.py --jobs 2                          # run steps in parallel
    python validation/run_all.py --config validation/run_all_nv_mgo.toml   # NV- in MgO, 3E -> 1A1 (its own suite dir)

Every parameter (physics, which methods, the fd grid and stencil, the fits and spline options, the cubic
variants with a fixed or fitted, tolerances, which steps) lives in a TOML file; ``validation/run_all.toml``
is the commented default and documents every key. ``--set section.key=value`` overrides one key (the value
is TOML, so lists and strings work), and missing keys keep their defaults.

From that configuration the suite builds and runs, as separate processes:

    checks     check_1d_ccd.py: the 81 regression checks (they use fixed reference parameters)
    tutorial   run_1d_ccd.py: tutorial 002 with these physics parameters
    compare    compare_vof.py: the harmonic methods like for like (must PASS); fd on every fit and every cubic
               variant against the reference; every cubic variant against the exact fd solution of its cubic
    plots      plot_ccd_wavefunctions.py: the CCD with Hermite, fd (every fit) and cubic (every variant)
               wavefunctions and their VOFs
    report     vof_report.py: the full seven-section report

The dataset (``suite.dataset``: nv_diamond, or nv_mgo for NV- in MgO) is passed to every step. A suite
directory remembers its dataset, and the suite refuses to replace one made from a different dataset, so the
MgO suite (default name ``nv_mgo``) can never overwrite the diamond one (``default``).

Fits are validated first: an unknown spec stops the suite, and a fit without bound states (e.g. poly4 here)
is skipped with the reason. Re-running a suite never deletes the previous one: it is moved to
``validation/outputs/_archive/suite/<name>/<YYYYmmdd-HHMMSS>/`` first. Everything goes to
``validation/outputs/suite/<name>/``, one directory per step
exactly as each script would write it alone, plus

    index.md            status, time and output link of every step, and a method-comparison table
    suite_config.json   the resolved configuration (every value actually used)
    <config>.toml       a copy of the configuration file
    logs/               the full console output of every step

Exit status is 0 when nothing errored and every like-for-like step passed. Comparisons of anharmonic
methods or fits are *expected* to differ from the harmonic reference; they are reported as DIFFERS.
"""

import argparse
import copy
import csv
import datetime
import json
import os
import re
import shutil
import subprocess
import sys
import time
import tomllib
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass, field

from harness import REPO_ROOT, VALIDATION_DIR, ccd_data, fits
from harness.methods import METHODS, VOFParams, fd_q_range
from harness.output import archive_previous, markdown_table

DEFAULT_CONFIG = os.path.join(VALIDATION_DIR, "run_all.toml")
SUITE_ROOT = os.path.join(VALIDATION_DIR, "outputs", "suite")
STEP_GROUPS = ("checks", "tutorial", "compare", "plots", "report")
# Rows of the method-comparison table, in this order.
KIND_ORDER = ("like for like", "PT vs exact cubic", "anharmonic vs harmonic")

# Every key the configuration may set, with its default. A user file only overrides what it names.
DEFAULTS = {
    "suite": {
        "name": "default",
        "dataset": ccd_data.DEFAULT_DATASET,
        "jobs": 1,
        "steps": {g: True for g in STEP_GROUPS},
    },
    "physics": {
        "temp": 5.0,
        "sigma": [8.0, 25.0],
        "gamma": 0.3,
        "zpl_lorentzian": True,
        "order_es": 50,
        "order_gs": 60,
        "ene_range": [-200.0, 1200.0],
        "resol": 1401,
        "ezpl": 1945.0,
        "tdm": 1.0,
        "line_positions": "harmonic",
    },
    "methods": {"reference": "pypl", "harmonic": ["hermite", "fd"]},
    "hermite": {"npts": 100001, "h": 0.0, "margin": 6.0},
    "fd": {
        "npts": 20001,
        "h": 0.0,
        "range": [],
        "stencil": 5,
        "margin": 2.0,
        "shift": "auto",
        "fits": [],
        "splines": [],
    },
    "cubic": {"enabled": True, "warn": 0.1, "variants": [{"points": "all", "free_a": False, "normalize": False}]},
    "plots": {"n_levels": 3, "check_tol": [1e-10, 1e-8]},
    "report": {"convergence": True, "cubic": True, "fits": [], "options": []},
    "figures": {"dpi": 300, "formats": ["png"], "font_size": 9.0},
    "checks": {"tolerances": {}},
    "tolerances": {"fc": 1e-7, "vof": 1e-8},
}
LIST_OF_TABLES = {
    ("fd", "splines"): {"order", "bc", "smoothing", "tails"},
    ("cubic", "variants"): {"points", "free_a", "normalize"},
}


# ------------------------------------------------------------------------- configuration


def merge(base, override, path=""):
    """``base`` updated with ``override``, rejecting keys ``base`` does not have (typos)."""

    out = copy.deepcopy(base)
    for key, value in override.items():
        where = "%s.%s" % (path, key) if path else key
        if key not in base:
            raise SystemExit("Unknown configuration key %r; allowed here: %s" % (where, ", ".join(sorted(base))))
        if isinstance(base[key], dict) and isinstance(value, dict) and base[key]:
            out[key] = merge(base[key], value, where)
        else:
            out[key] = value
    return out


def set_override(config, assignment):
    """Apply one ``--set section.key=value``; the value is parsed as TOML (strings may be left unquoted)."""

    if "=" not in assignment:
        raise SystemExit("--set expects section.key=value, got %r" % assignment)
    path, raw = assignment.split("=", 1)
    try:
        value = tomllib.loads("v = %s" % raw)["v"]
    except tomllib.TOMLDecodeError:
        value = raw
    keys = path.strip().split(".")
    node = {}
    leaf = node
    for key in keys[:-1]:
        leaf[key] = {}
        leaf = leaf[key]
    leaf[keys[-1]] = value
    return merge(config, node)


def load_config(args):
    with open(args.config, "rb") as f:
        user = tomllib.load(f)
    config = merge(DEFAULTS, user)
    for assignment in args.set or []:
        config = set_override(config, assignment)
    if args.name:
        config["suite"]["name"] = args.name
    if args.jobs:
        config["suite"]["jobs"] = args.jobs
    for (section, key), allowed in LIST_OF_TABLES.items():
        for entry in config[section][key]:
            unknown = set(entry) - allowed
            if unknown:
                raise SystemExit(
                    "Unknown key(s) %s in [[%s.%s]]; allowed: %s" % (sorted(unknown), section, key, sorted(allowed))
                )
    if config["suite"]["dataset"] not in ccd_data.DATASETS:
        raise SystemExit(
            "suite.dataset = %r is not a dataset; choose from %s"
            % (config["suite"]["dataset"], ", ".join(ccd_data.DATASETS))
        )
    if config["suite"]["dataset"] != ccd_data.DEFAULT_DATASET and config["suite"]["steps"]["tutorial"]:
        raise SystemExit(
            "steps.tutorial is tutorial 002 (NV- in diamond PL / Abs) and has no %s version; set it to false"
            % config["suite"]["dataset"]
        )
    if config["methods"]["reference"] not in METHODS or METHODS[config["methods"]["reference"]].anharmonic:
        raise SystemExit("methods.reference must be a harmonic method: pypl, hermite or fd")
    for m in config["methods"]["harmonic"]:
        if m not in METHODS or METHODS[m].anharmonic:
            raise SystemExit("methods.harmonic may only list pypl, hermite and fd, got %r" % m)
    return config


def resolve_fits(config, ccd):
    """``(usable, skipped)``: canonical fd fit specs with bound states, and ``{spec: reason}`` for the rest."""

    specs = list(config["fd"]["fits"])
    for s in config["fd"]["splines"]:
        try:
            specs.append(
                fits.spline_spec(s.get("order", 3), s.get("bc"), s.get("smoothing"), s.get("tails", "harmonic"))
            )
        except ValueError as err:
            raise SystemExit("Invalid [[fd.splines]] entry %s: %s" % (s, err))
    params = grid_params(config)
    usable, skipped = [], {}
    for spec in specs:
        try:
            spec = fits.canonical(spec)
        except ValueError as err:
            raise SystemExit("Invalid fd fit in the configuration: %s" % err)
        if spec == fits.DEFAULT_FIT or spec in usable or spec in skipped:
            continue
        try:
            fits.check_bounded(ccd, spec, fd_q_range(ccd, params))
            usable.append(spec)
        except ValueError as err:
            skipped[spec] = str(err)
    return usable, skipped


# ---------------------------------------------------------------------------------- steps


@dataclass
class Step:
    group: str
    title: str
    script: str
    args: list
    expect: str = "pass"  # pass: must succeed; differ: an anharmonic comparison, a FAIL verdict is expected
    status: str = "not run"
    seconds: float = 0.0
    output: str = ""  # run directory, relative to the repository
    log: str = ""
    details: dict = field(default_factory=dict)

    def command(self):
        return " ".join(["python", "validation/" + self.script] + [str(a) for a in self.args])


def grid_params(config):
    """VOFParams with the basis and fd grid of ``config`` (enough to know the fd grid's Q range)."""

    fd = config["fd"]
    return VOFParams(
        order_es=config["physics"]["order_es"],
        order_gs=config["physics"]["order_gs"],
        fd_margin=fd["margin"],
        fd_range=tuple(fd["range"]) if fd["range"] else None,
    )


def figure_args(config):
    f = config["figures"]
    return ["--dpi", f["dpi"], "--fig-format", *f["formats"], "--font-size", f["font_size"]]


def vof_args(config, fd_fit="parabola", cubic=None):
    """The shared VOF flags of compare_vof.py / plot_ccd_wavefunctions.py / vof_report.py."""

    p, fd, herm = config["physics"], config["fd"], config["hermite"]
    args = ["--dataset", config["suite"]["dataset"]]
    args += ["--temp", p["temp"], "--sigma", *p["sigma"], "--gamma", p["gamma"]]
    args += [] if p["zpl_lorentzian"] else ["--no-zpl"]
    args += ["--order-es", p["order_es"], "--order-gs", p["order_gs"], "--ene-range", *p["ene_range"]]
    args += ["--resol", p["resol"], "--line-positions", p["line_positions"]]
    args += ["--npts", herm["npts"], "--hermite-margin", herm["margin"]]
    args += ["--hermite-h", herm["h"]] if herm["h"] else []
    args += ["--fd-stencil", fd["stencil"], "--fd-margin", fd["margin"], "--fd-fit", fd_fit]
    args += ["--fd-h", fd["h"]] if fd["h"] else ["--fd-npts", fd["npts"]]
    args += ["--fd-range", *fd["range"]] if fd["range"] else []
    args += [] if fd["shift"] == "auto" else ["--fd-shift", fd["shift"]]
    args += ["--cubic-warn", config["cubic"]["warn"]]
    args += figure_args(config)
    if cubic is not None:
        args += ["--cubic-points", cubic.get("points", "all")]
        args += ["--cubic-free-a"] if cubic.get("free_a", False) else []
        args += ["--cubic-normalize"] if cubic.get("normalize", False) else []
    return args


def cubic_spec(variant):
    return (
        "cubic"
        + ("-near" if variant.get("points", "all") == "parabola" else "")
        + ("-free" if variant.get("free_a", False) else "")
    )


def resolve_cubic_exact(config, ccd):
    """``{spec: reason}`` of the cubic variants with no exact fd solution (the cubic turns over inside the grid)."""

    params = grid_params(config)
    out = {}
    for variant in config["cubic"]["variants"] if config["cubic"]["enabled"] else []:
        spec = cubic_spec(variant)
        try:
            fits.check_bounded(ccd, spec, fd_q_range(ccd, params))
        except ValueError as err:
            out[spec] = str(err)
    return out


def build_steps(config, usable_fits, no_exact_cubic=()):
    ref = config["methods"]["reference"]
    harmonic = [m for m in config["methods"]["harmonic"] if m != ref]
    tol = ["--tol-fc", config["tolerances"]["fc"], "--tol-vof", config["tolerances"]["vof"]]
    variants = config["cubic"]["variants"] if config["cubic"]["enabled"] else []
    on = config["suite"]["steps"]
    steps = []

    if on["checks"]:
        tols = ["--tol=%s=%r" % kv for kv in sorted(config["checks"]["tolerances"].items())]
        steps.append(Step("checks", "81 regression and physics checks", "check_1d_ccd.py", ["-v"] + tols))
    if on["tutorial"]:
        p = config["physics"]
        args = ["--temp", p["temp"], "--ezpl", p["ezpl"], "--order-es", p["order_es"], "--order-gs", p["order_gs"]]
        args += ["--sigma", *p["sigma"][:2] * (2 if len(p["sigma"]) == 1 else 1), "--gamma", p["gamma"]]
        args += ["--ene-range", *p["ene_range"], "--resol", p["resol"], "--tdm", p["tdm"]]
        args += ["--dataset", config["suite"]["dataset"]] + figure_args(config)
        steps.append(Step("tutorial", "tutorial 002 spectra (PyPL's own solver)", "run_1d_ccd.py", args))

    if on["compare"]:
        if harmonic:
            steps.append(
                Step(
                    "compare",
                    "harmonic methods like for like: %s vs %s" % (ref, ", ".join(harmonic)),
                    "compare_vof.py",
                    ["--reference", ref, "--test", *harmonic] + vof_args(config) + tol,
                )
            )
        for spec in usable_fits:
            steps.append(
                Step(
                    "compare",
                    "fd on %s vs %s" % (fits.label(spec), ref),
                    "compare_vof.py",
                    ["--reference", ref, "--test", "fd"] + vof_args(config, spec) + tol,
                    "differ",
                )
            )
        for variant in variants:
            spec = cubic_spec(variant)
            norm = " (normalized)" if variant.get("normalize") else ""
            steps.append(
                Step(
                    "compare",
                    "cubic PT [%s]%s vs %s" % (fits.label(spec), norm, ref),
                    "compare_vof.py",
                    ["--reference", ref, "--test", "cubic"] + vof_args(config, cubic=variant) + tol,
                    "differ",
                )
            )
            if spec in no_exact_cubic:
                continue
            steps.append(
                Step(
                    "compare",
                    "cubic PT [%s]%s vs the exact fd solution of that cubic" % (fits.label(spec), norm),
                    "compare_vof.py",
                    ["--reference", "fd", "--test", "cubic"] + vof_args(config, spec, cubic=variant) + tol,
                    "differ",
                )
            )

    if on["plots"]:
        n = ["--n-levels", config["plots"]["n_levels"], "--check-tol", *config["plots"]["check_tol"]]
        steps.append(
            Step(
                "plots",
                "Hermite wavefunctions, VOF of %s" % ", ".join([ref] + harmonic),
                "plot_ccd_wavefunctions.py",
                ["--wf-method", "hermite", "--vof-method", ref, *harmonic] + vof_args(config) + n,
            )
        )
        for spec in [fits.DEFAULT_FIT] + usable_fits:
            steps.append(
                Step(
                    "plots",
                    "fd wavefunctions on %s" % fits.label(spec),
                    "plot_ccd_wavefunctions.py",
                    ["--wf-method", "fd", "--vof-method", ref, "fd"] + vof_args(config, spec) + n,
                    "pass" if spec == fits.DEFAULT_FIT else "differ",
                )
            )
        for variant in variants:
            steps.append(
                Step(
                    "plots",
                    "cubic PT wavefunctions [%s]" % fits.label(cubic_spec(variant)),
                    "plot_ccd_wavefunctions.py",
                    ["--wf-method", "cubic", "--vof-method", ref, "cubic"] + vof_args(config, cubic=variant) + n,
                    "differ",
                )
            )

    if on["report"]:
        rep = config["report"]
        args = ["--reference", ref] + vof_args(config, cubic=variants[0] if variants else None) + tol
        args += ["--fits", *rep["fits"]] if rep["fits"] else []
        args += [] if rep["convergence"] else ["--skip-convergence"]
        args += [] if (rep["cubic"] and variants) else ["--skip-cubic"]
        args += [str(a) for a in rep["options"]]
        steps.append(Step("report", "full seven-section report", "vof_report.py", args))
    return steps


# ----------------------------------------------------------------------------------- run


OUTPUT_RE = re.compile(r"Outputs in (.+?) \(see index\.md\)")


def run_step(step, number, total, suite_dir, dataset):
    env = dict(os.environ, PYPL_VALIDATION_OUTPUT_ROOT=suite_dir, PYPL_VALIDATION_DATASET=dataset, MPLBACKEND="Agg")
    start = time.perf_counter()
    proc = subprocess.run(
        [sys.executable, "-u", os.path.join(VALIDATION_DIR, step.script)] + [str(a) for a in step.args],
        cwd=REPO_ROOT,
        env=env,
        capture_output=True,
        text=True,
        encoding="utf-8",
        errors="replace",
    )
    step.seconds = time.perf_counter() - start
    text = proc.stdout + ("\n--- stderr ---\n" + proc.stderr if proc.stderr.strip() else "")
    log_name = "%02d_%s_%s.log" % (number, step.group, re.sub(r"[^A-Za-z0-9]+", "-", step.title).strip("-")[:60])
    step.log = os.path.join("logs", log_name)
    with open(os.path.join(suite_dir, step.log), "w", encoding="utf-8") as f:
        f.write("$ %s\n\n%s" % (step.command(), text))

    match = OUTPUT_RE.search(proc.stdout)
    step.output = match.group(1).replace("\\", "/") if match else ""
    crashed = "Traceback (most recent call last)" in text or (
        proc.returncode != 0 and not match and step.group != "checks"
    )
    if crashed:
        step.status = "ERROR"
    elif proc.returncode == 0:
        step.status = "PASS" if step.group in ("checks", "compare", "report") else "OK"
    else:
        step.status = "DIFFERS (expected)" if step.expect == "differ" else "FAIL"
    if step.group == "checks":
        summary = re.search(r"(\d+)/(\d+) checks passed", proc.stdout)
        step.details["checks"] = summary.group(0) if summary else "no summary"
    print(
        "[%2d/%d] %-8s %-66s %6.1f s  %s" % (number, total, step.group, step.title[:66], step.seconds, step.status),
        flush=True,
    )
    return step


def compare_rows(step):
    """The per-direction rows of a compare step's data/summary.csv."""

    path = os.path.join(REPO_ROOT, step.output, "data", "summary.csv")
    if not step.output or not os.path.exists(path):
        return []
    with open(path, encoding="utf-8") as f:
        return list(csv.DictReader(f))


# ----------------------------------------------------------------------------------- index


def link(suite_dir, step):
    if not step.output:
        return "`%s`" % step.log.replace("\\", "/")
    rel = os.path.relpath(os.path.join(REPO_ROOT, step.output), suite_dir).replace("\\", "/")
    return "[`%s/`](%s/index.md)" % (rel, rel)


def write_index(suite_dir, config, steps, skipped, started, config_path):
    lines = [
        "# PyPL 1D CCD test suite: %s (%s)"
        % (config["suite"]["name"], ccd_data.DATASETS[config["suite"]["dataset"]].title.replace("$^-$", "-")),
        "",
        "Started %s, %d steps in %.0f s. Configuration: `%s` (copied here), resolved values in `suite_config.json`."
        % (
            started.isoformat(timespec="seconds"),
            len(steps),
            sum(s.seconds for s in steps),
            os.path.basename(config_path),
        ),
        "",
    ]
    ds = ccd_data.DATASETS[config["suite"]["dataset"]]
    lines += ["**Dataset.**", ""] + ["- %s" % line for line in ds.describe()] + [""]
    p = config["physics"]
    lines += [
        "**Physics:** T = %g K, sigma = %s meV, ZPL %s, %d ES x %d GS levels, axis %g..%g meV (%d points). "
        "**Reference:** %s. **fd:** %d points, %d-point stencil, margin %g. **Hermite:** %d points."
        % (
            p["temp"],
            p["sigma"],
            "Lorentzian, gamma = %g meV" % p["gamma"] if p["zpl_lorentzian"] else "Gaussian",
            p["order_es"],
            p["order_gs"],
            *p["ene_range"],
            p["resol"],
            METHODS[config["methods"]["reference"]].label,
            config["fd"]["npts"],
            config["fd"]["stencil"],
            config["fd"]["margin"],
            config["hermite"]["npts"],
        ),
        "",
    ]
    bad = [s for s in steps if s.status in ("ERROR", "FAIL")]
    lines += [
        "**Result: %s.**"
        % (
            "all like-for-like steps passed, nothing errored"
            if not bad
            else "%d step(s) failed or errored: %s" % (len(bad), "; ".join(s.title for s in bad))
        ),
        "",
    ]
    if skipped:
        lines += ["Skipped (no bound states on the fd grid):", ""]
        lines += ["- `%s`: %s" % (spec, reason.split(" (")[0]) for spec, reason in skipped.items()] + [""]

    # Method comparison: one row per compare step and test.
    rows = []
    for s in steps:
        if s.group != "compare":
            continue
        for r in compare_rows(s):
            exact = s.expect == "differ" and r["reference"].startswith("Finite difference [cubic")
            kind = (
                "like for like" if s.expect == "pass" else ("PT vs exact cubic" if exact else "anharmonic vs harmonic")
            )
            rows.append(
                {
                    "test": r["test"],
                    "reference": r["reference"],
                    "direction": r["direction"],
                    "fc": float(r["fc_max"]),
                    "vof": float(r["vof_max"]),
                    "seconds": float(r["seconds"]),
                    "verdict": r["verdict"] if s.expect == "pass" else "differs (expected)",
                    "kind": kind,
                }
            )
    if rows:
        merged = {}
        for r in rows:
            key = (r["test"], r["reference"])
            m = merged.setdefault(key, dict(r, vof=0.0, fc=0.0, seconds=0.0, n=0))
            m["vof"], m["fc"] = max(m["vof"], r["vof"]), max(m["fc"], r["fc"])
            m["seconds"] += r["seconds"]
            m["n"] += 1
        table = sorted(merged.values(), key=lambda m: (KIND_ORDER.index(m["kind"]), m["reference"], m["vof"]))
        for m in table:
            m["seconds"] /= m["n"]
        lines += [
            "## Method comparison",
            "",
            "Largest difference over emission and absorption. **Like for like:** methods solving the same harmonic "
            "problem; they must agree within the tolerances (max dF_ij <= %g, max dVOF / peak <= %g). **PT vs exact "
            "cubic:** the error of first-order perturbation theory, against the exact fd solution of the same cubic. "
            "**Anharmonic vs harmonic:** how much each fit or correction changes the VOF from the harmonic reference "
            "(physics, not error), smallest first. See `validation/METHODS.md` for what each method does."
            % (config["tolerances"]["fc"], config["tolerances"]["vof"]),
            "",
        ]
        columns = [
            ("kind", "kind", "%s"),
            ("test", "method", "%s"),
            ("reference", "compared with", "%s"),
            ("vof", "max dVOF / peak", "%.2e"),
            ("fc", "max dF_ij", "%.2e"),
            ("seconds", "F_ij time (s)", "%.2f"),
            ("verdict", "verdict", "%s"),
        ]
        lines += markdown_table(table, columns) + [""]

    lines += ["## Steps", "", "| # | step | status | time (s) | output |", "| --- | --- | --- | --- | --- |"]
    for i, s in enumerate(steps, 1):
        extra = " (%s)" % s.details["checks"] if "checks" in s.details else ""
        lines.append(
            "| %d | %s: %s | %s%s | %.1f | %s |" % (i, s.group, s.title, s.status, extra, s.seconds, link(suite_dir, s))
        )
    lines += (
        [
            "",
            "## Commands",
            "",
            "Each step was run as (from the repository root, with `PYPL_VALIDATION_OUTPUT_ROOT` "
            "pointing at this directory and `PYPL_VALIDATION_DATASET=%s`):" % config["suite"]["dataset"],
            "",
            "```bash",
        ]
        + [s.command() for s in steps]
        + ["```", ""]
    )
    with open(os.path.join(suite_dir, "index.md"), "w", encoding="utf-8") as f:
        f.write("\n".join(lines))


# ------------------------------------------------------------------------------------ main


def existing_dataset(suite_dir):
    """The dataset of the suite already in ``suite_dir`` (suites from before datasets existed: nv_diamond)."""

    path = os.path.join(suite_dir, "suite_config.json")
    if not os.path.exists(path):
        return None
    with open(path, encoding="utf-8") as f:
        return json.load(f).get("suite", {}).get("dataset", ccd_data.DEFAULT_DATASET)


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--config", default=DEFAULT_CONFIG, help="TOML configuration (default: validation/run_all.toml)")
    p.add_argument("--set", action="append", metavar="SECTION.KEY=VALUE", help="override one key (repeatable)")
    p.add_argument("--name", help="suite name, i.e. the directory validation/outputs/suite/<name>/")
    p.add_argument("--only", nargs="+", choices=STEP_GROUPS, help="run only these step groups")
    p.add_argument("--jobs", type=int, help="steps run in parallel (default: suite.jobs)")
    p.add_argument("--dry-run", action="store_true", help="print every step's command and exit")
    return p.parse_args()


def main():
    args = parse_args()
    config = load_config(args)
    if args.only:
        config["suite"]["steps"] = {g: g in args.only for g in STEP_GROUPS}
    dataset = ccd_data.use_dataset(config["suite"]["dataset"])
    ccd = ccd_data.build_ccd(verbose=False)
    usable, skipped = resolve_fits(config, ccd)
    no_exact_cubic = resolve_cubic_exact(config, ccd)
    steps = build_steps(config, usable, no_exact_cubic)
    for spec, reason in no_exact_cubic.items():
        skipped["%s (exact fd check of cubic PT)" % spec] = reason

    name = "".join(c if c.isalnum() or c in "-_.+" else "_" for c in config["suite"]["name"]).strip(".") or "default"
    config["suite"]["name"] = name
    print(
        "PyPL 1D CCD test suite '%s' on %s (%s): %d steps; fd fits %s"
        % (name, dataset.name, dataset.transition, len(steps), ", ".join(usable) or "none")
    )
    for spec, reason in skipped.items():
        print("  skipping fit %s: %s" % (spec, reason.split(" (")[0]))
    if args.dry_run:
        for i, s in enumerate(steps, 1):
            print("[%2d] %-8s %s\n     %s" % (i, s.group, s.title, s.command()))
        return

    suite_dir = os.path.join(SUITE_ROOT, name)
    previous = existing_dataset(suite_dir)
    if previous is not None and previous != dataset.name:
        raise SystemExit(
            "%s holds a suite run on dataset %s; refusing to replace it with a %s run. Choose another --name."
            % (os.path.relpath(suite_dir, REPO_ROOT), previous, dataset.name)
        )
    # Never delete a previous suite: move it to validation/outputs/_archive/suite/<name>/<stamp>/.
    archived = archive_previous(suite_dir)
    if archived:
        print("Previous '%s' suite moved to %s" % (name, os.path.relpath(archived, REPO_ROOT)))
    os.makedirs(os.path.join(suite_dir, "logs"))
    shutil.copy(args.config, suite_dir)
    with open(os.path.join(suite_dir, "suite_config.json"), "w", encoding="utf-8") as f:
        json.dump(dict(config, resolved_fits=usable, skipped_fits=skipped), f, indent=2)
    print("Output: %s\n" % os.path.relpath(suite_dir, REPO_ROOT))

    started = datetime.datetime.now()
    jobs = max(1, int(config["suite"]["jobs"]))
    with ThreadPoolExecutor(max_workers=jobs) as pool:
        futures = [pool.submit(run_step, s, i, len(steps), suite_dir, dataset.name) for i, s in enumerate(steps, 1)]
        for future in futures:
            future.result()

    write_index(suite_dir, config, steps, skipped, started, args.config)
    bad = [s for s in steps if s.status in ("ERROR", "FAIL")]
    print(
        "\n%d steps in %.0f s: %d passed / ok, %d differ as expected, %d failed or errored."
        % (
            len(steps),
            sum(s.seconds for s in steps),
            sum(s.status in ("PASS", "OK") for s in steps),
            sum(s.status.startswith("DIFFERS") for s in steps),
            len(bad),
        )
    )
    for s in bad:
        print("  %s: %s  (log: %s)" % (s.status, s.title, os.path.join(os.path.relpath(suite_dir, REPO_ROOT), s.log)))
    print("Summary: %s" % os.path.relpath(os.path.join(suite_dir, "index.md"), REPO_ROOT))
    raise SystemExit(1 if bad else 0)


if __name__ == "__main__":
    main()
