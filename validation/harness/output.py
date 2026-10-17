"""Organized output for the test scripts.

Every run of a script writes into its own directory,

    validation/outputs/<script>/<run>/
        index.md        what was run, the headline results, and every file with a one-line description
        log.txt         the full console output, including every method's settings
        settings.json   every parameter of the run
        <section>/      figures (.png) and data (.csv) grouped by topic

so runs never clobber each other unless they share a run name. A re-run with the same name starts from
an empty directory (nothing stale is left behind), but the previous run is never deleted: it is moved to

    validation/outputs/_archive/<same path>/<YYYYmmdd-HHMMSS of the previous run>/

by ``archive_previous``, which run_all.py uses for whole suites too.

Runs on a dataset other than the default (``--dataset nv_mgo``) go to ``validation/outputs/<dataset>/<script>/<run>/``
instead, so they can never replace a diamond run; their index.md and log.txt open with the dataset.
"""

import csv
import datetime
import json
import os
import shutil
import sys

import numpy as np

from harness import REPO_ROOT, VALIDATION_DIR, ccd_data

# run_all.py points this at its suite directory so every script of a suite writes beneath it.
OUTPUT_ROOT = os.environ.get("PYPL_VALIDATION_OUTPUT_ROOT") or os.path.join(VALIDATION_DIR, "outputs")


ARCHIVE_ROOT = os.path.join(VALIDATION_DIR, "outputs", "_archive")


def archive_previous(path):
    """Move an existing run directory ``path`` out of the way instead of deleting it; return where it went.

    Under ``validation/outputs/`` it goes to ``validation/outputs/_archive/<path relative to outputs>/<stamp>/``;
    elsewhere (a custom ``PYPL_VALIDATION_OUTPUT_ROOT``) to ``<parent>/_archive/<name>/<stamp>/``. The stamp is
    the time the previous run last wrote a file, so archived runs sort chronologically. Returns None if
    ``path`` does not exist.
    """

    if not os.path.isdir(path):
        return None
    newest = max(
        [os.path.getmtime(path)]
        + [os.path.getmtime(os.path.join(d, f)) for d, _, files in os.walk(path) for f in files]
    )
    stamp = datetime.datetime.fromtimestamp(newest).strftime("%Y%m%d-%H%M%S")
    outputs = os.path.join(VALIDATION_DIR, "outputs")
    rel = os.path.relpath(os.path.abspath(path), outputs)
    if rel.startswith("..") or rel.split(os.sep)[0] == "_archive":
        base = os.path.join(os.path.dirname(os.path.abspath(path)), "_archive", os.path.basename(path))
    else:
        base = os.path.join(ARCHIVE_ROOT, rel)
    target, k = os.path.join(base, stamp), 1
    while os.path.exists(target):
        k += 1
        target = os.path.join(base, "%s-%d" % (stamp, k))
    os.makedirs(base, exist_ok=True)
    shutil.move(path, target)
    return target


def output_root():
    """Where run directories go: the suite directory if set, else outputs/ (default dataset) or outputs/<dataset>/."""

    if os.environ.get("PYPL_VALIDATION_OUTPUT_ROOT"):
        return OUTPUT_ROOT
    if ccd_data.ACTIVE == ccd_data.DEFAULT_DATASET:
        return OUTPUT_ROOT
    return os.path.join(OUTPUT_ROOT, ccd_data.ACTIVE)


class OutputDir:
    """One run directory: ``outputs/<script>/<run>/``, with a log, an index and topic subdirectories.

    Parameters
    ----------
    script : str
        Name of the producing script; the first directory level.
    run : str
        Name of this run (e.g. the method choice or a ``--tag``); the second level.
    title : str
        Heading of ``index.md``.
    verbose : bool
        Echo detail lines (``log(..., detail=True)``) to the console too; they always reach log.txt.
    keep_figures : bool
        Leave figures open after saving them (for ``--show``); by default they are closed.
    """

    def __init__(self, script, run="default", title=None, verbose=False, keep_figures=False):
        self.script = script
        # The run name becomes a directory name, so keep it a single safe path component.
        self.run = "".join(c if c.isalnum() or c in "-_.+" else "_" for c in (run or "default")).strip(".") or "default"
        self.path = os.path.join(output_root(), script, self.run)
        self.title = title or "%s: %s" % (script, self.run)
        self.dataset = ccd_data.active()
        if self.dataset.name != ccd_data.DEFAULT_DATASET:
            self.title += " [%s]" % self.dataset.title.replace("$^-$", "-")
        self.verbose = verbose
        self.keep_figures = keep_figures
        archived = archive_previous(self.path)
        os.makedirs(self.path)
        self._log = open(os.path.join(self.path, "log.txt"), "w", encoding="utf-8")
        self.files = []  # (section, relative path, description)
        self.sections = {}  # section -> description
        self.summary = []  # markdown lines for the top of index.md
        self.log("# %s, started %s" % (self.title, datetime.datetime.now().isoformat(timespec="seconds")), detail=True)
        self.log("# command: %s" % " ".join([os.path.basename(sys.executable)] + sys.argv), detail=True)
        if archived:
            self.log("# previous run with this name moved to %s" % os.path.relpath(archived, REPO_ROOT))
        # Every run records its dataset: data, state roles and how each state's energies were computed. On the
        # default dataset it goes to log.txt only, to keep the console short.
        for line in self.dataset.describe():
            self.log("# %s" % line, detail=self.dataset.name == ccd_data.DEFAULT_DATASET)

    # ------------------------------------------------------------------------ logging

    def log(self, line="", detail=False):
        """Write ``line`` to log.txt; print it as well unless it is a detail line and not verbose."""

        self._log.write(line + "\n")
        self._log.flush()
        if not detail or self.verbose:
            print(line)

    def detail(self, line=""):
        """A line for log.txt only (and the console when verbose)."""

        self.log(line, detail=True)

    # ------------------------------------------------------------------------- files

    def section(self, name, description=""):
        """Create (if needed) and return the absolute path of subdirectory ``name``."""

        path = os.path.join(self.path, name)
        os.makedirs(path, exist_ok=True)
        if name not in self.sections:
            self.sections[name] = description
        return path

    def _register(self, relpath, description):
        section = relpath.split("/")[0] if "/" in relpath else ""
        self.files.append((section, relpath, description))
        self.detail("  wrote %s" % os.path.relpath(os.path.join(self.path, relpath), REPO_ROOT))

    def _prepare(self, relpath):
        path = os.path.join(self.path, relpath)
        os.makedirs(os.path.dirname(path), exist_ok=True)
        return path

    def save_fig(self, fig, relpath, description=""):
        """Save a matplotlib figure at ``relpath`` (e.g. ``"05_vof_methods/vof_methods.png"``).

        The figure is written in every format of ``plots.FIGURE["formats"]`` (``--fig-format``) at
        ``plots.FIGURE["dpi"]`` (``--dpi``); ``relpath``'s extension is replaced by each format, and the first one is
        listed in index.md. On a dataset other than the default, a small caption at the top right names the dataset
        and the state roles.
        """

        from harness import plots

        if self.dataset.name != ccd_data.DEFAULT_DATASET:
            labels = self.dataset.state_labels or self.dataset.state_names
            fig.text(
                0.995,
                1.0,
                "%s: %s $\\rightarrow$ %s" % (self.dataset.title, labels["es"], labels["gs"]),
                ha="right",
                va="bottom",
                fontsize=plots.FIGURE["font_size"] - 2,
                color="0.4",
            )
        stem = os.path.splitext(relpath)[0]
        paths = ["%s.%s" % (stem, fmt) for fmt in plots.FIGURE["formats"]]
        for path in paths:
            fig.savefig(self._prepare(path), dpi=plots.FIGURE["dpi"], bbox_inches="tight", pad_inches=0.04)
        self._register(paths[0], description)
        for path in paths[1:]:
            self._register(
                path, "The same figure as %s, in %s format." % (os.path.basename(paths[0]), path.rsplit(".", 1)[1])
            )
        if not self.keep_figures:
            import matplotlib.pyplot as plt

            plt.close(fig)

    def save_csv(self, relpath, columns, description="", comments=()):
        """Write a CSV with one column per ``columns`` item (``{header: 1D array}``); ``comments`` lead with '#'."""

        headers = list(columns)
        data = [np.asarray(columns[h]) for h in headers]
        with open(self._prepare(relpath), "w", newline="", encoding="utf-8") as f:
            for c in comments:
                f.write("# %s\n" % c)
            writer = csv.writer(f)
            writer.writerow(headers)
            for row in zip(*data):
                writer.writerow([_fmt(v) for v in row])
        self._register(relpath, description)

    def save_table(self, relpath, rows, description="", comments=()):
        """Write a CSV from a list of dicts sharing the same keys (a summary table)."""

        with open(self._prepare(relpath), "w", newline="", encoding="utf-8") as f:
            for c in comments:
                f.write("# %s\n" % c)
            writer = csv.DictWriter(f, fieldnames=list(rows[0]))
            writer.writeheader()
            for row in rows:
                writer.writerow({k: _fmt(v) for k, v in row.items()})
        self._register(relpath, description)

    def save_json(self, relpath, obj, description=""):
        with open(self._prepare(relpath), "w", encoding="utf-8") as f:
            json.dump(obj, f, indent=2, default=_jsonable)
        self._register(relpath, description)

    def add_summary(self, *lines):
        """Markdown lines placed under the title of index.md."""

        self.summary.extend(lines)

    # --------------------------------------------------------------------------- close

    def close(self):
        """Write index.md and close log.txt. Returns the run directory."""

        lines = [
            "# %s" % self.title,
            "",
            "Generated %s by `%s`."
            % (datetime.datetime.now().isoformat(timespec="seconds"), " ".join(["python"] + sys.argv)),
            "",
        ]
        lines += ["**Dataset.**", ""] + ["- %s" % line for line in self.dataset.describe()] + [""]
        if self.summary:
            lines += self.summary + [""]
        lines += [
            "## Files",
            "",
            "| File | Content |",
            "| --- | --- |",
            "| `log.txt` | Full console output, " "including every method's settings. |",
            "| `settings.json` | Every parameter of the run. |",
        ]
        current = None
        for section, relpath, description in sorted(self.files, key=lambda f: (f[0], f[1])):
            if relpath == "settings.json":
                continue
            if section != current and section:
                current = section
                lines.append("| **`%s/`** | %s |" % (section, self.sections.get(section, "")))
            lines.append("| `%s` | %s |" % (relpath, description.replace("|", "\\|")))
        with open(os.path.join(self.path, "index.md"), "w", encoding="utf-8") as f:
            f.write("\n".join(lines) + "\n")
        self.log("\nOutputs in %s (see index.md)" % os.path.relpath(self.path, REPO_ROOT))
        self._log.close()
        return self.path


def text_table(rows, columns):
    """A fixed-width console table. ``columns`` is a list of ``(key, header, format)``."""

    cells = [[fmt % row[key] for key, _, fmt in columns] for row in rows]
    widths = [max(len(header), *(len(c[k]) for c in cells)) for k, (_, header, _) in enumerate(columns)]
    lines = ["  ".join(header.ljust(w) for (_, header, _), w in zip(columns, widths))]
    lines.append("  ".join("-" * w for w in widths))
    lines += ["  ".join(c.ljust(w) for c, w in zip(row, widths)) for row in cells]
    return lines


def markdown_table(rows, columns):
    """The same table as Markdown, for index.md."""

    lines = ["| %s |" % " | ".join(header for _, header, _ in columns), "|" + " --- |" * len(columns)]
    lines += ["| %s |" % " | ".join(fmt % row[key] for key, _, fmt in columns) for row in rows]
    return lines


def _fmt(value):
    if isinstance(value, (float, np.floating)):
        return "%.10g" % value
    return value


def _jsonable(obj):
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    if isinstance(obj, (np.floating, np.integer)):
        return obj.item()
    return str(obj)
