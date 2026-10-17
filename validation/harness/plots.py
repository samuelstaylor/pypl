"""Shared plot style, labels, colours and the figures more than one script draws.

Every figure of the validation suite goes through this module, so that the whole suite has one look:

- ``apply_style()`` sets the matplotlib defaults (fonts, math text, line widths, colour cycle); ``style(ax)`` adds
  inward ticks on all four sides and a light grid. ``FIGURE`` holds the user-configurable settings (resolution,
  file formats, base font size), set from the command line by ``add_figure_args`` / ``configure``.
- Labels are math text with upright units: ``Q_LABEL``, ``VOF_LABEL``, ``HW``, ``fij``.
- Colours: the two states follow the active dataset (``STATE_COLORS``); methods and fits have fixed colours and line
  styles (``METHOD_COLORS``, ``fit_style``) from the colour-blind-safe Okabe-Ito palette, so a method or fit looks
  the same in every figure. All spline variants share one colour and differ by line style.
- Long settings and per-curve statistics belong in ``log.txt``, the CSV files and ``index.md``, not on the figures.
"""

import textwrap

import matplotlib as mpl
import numpy as np

from harness import ccd_data

# --------------------------------------------------------------------------------------------- settings

# User-configurable figure settings (command line: --dpi, --fig-format, --font-size).
FIGURE = {"dpi": 300, "formats": ("png",), "font_size": 9.0}

# Okabe & Ito (2008) colour-blind-safe palette.
OKABE_ITO = {
    "black": "#000000",
    "orange": "#E69F00",
    "sky": "#56B4E9",
    "green": "#009E73",
    "yellow": "#F0E442",
    "blue": "#0072B2",
    "vermillion": "#D55E00",
    "purple": "#CC79A7",
    "gray": "#7F7F7F",
}

# Directions follow the state colours of the datasets (emission = upper state's red, absorption = lower's blue).
DIRECTION_COLORS = {"emission": "#DB4437", "absorption": "#4285F4"}
# The colour of each well follows the active dataset (diamond: GS red, ES blue; MgO: 3E red, 1A1 blue).
STATE_COLORS = ccd_data.STATE_COLORS
METHOD_COLORS = {
    "pypl": OKABE_ITO["black"],
    "hermite": OKABE_ITO["blue"],
    "fd": OKABE_ITO["green"],
    "cubic": OKABE_ITO["purple"],
}
METHOD_STYLES = {"pypl": "-", "hermite": "--", "fd": ":", "cubic": "-."}
# Families of fits: one colour each; variants within a family differ by line style.
FIT_FAMILY_COLORS = {
    "parabola": OKABE_ITO["black"],
    "poly": OKABE_ITO["orange"],
    "morse": OKABE_ITO["sky"],
    "cubic": OKABE_ITO["purple"],
    "spline": OKABE_ITO["green"],
}
FAMILY_STYLES = ("-", "--", ":", "-.", (0, (5, 1, 1, 1)), (0, (1, 1)), (0, (3, 1, 1, 1, 1, 1)))
TEST_STYLES = ("--", ":", "-.", (0, (1, 3)))
NOTE_WIDTH = 118

# Labels (math text, upright units).
ANGSTROM = r"\mathrm{\AA}"
Q_LABEL = r"$Q$ (amu$^{1/2}\,%s$)" % ANGSTROM
Q_UNIT = r"amu$^{1/2}\,%s$" % ANGSTROM
HW = r"$\hbar\omega$"
VOF_LABEL = r"VOF (meV$^{-1}$)"
E_ZPL_LABEL = r"$E - E_\mathrm{ZPL}$ (meV)"
ENERGY_LABEL = "Total energy (eV)"


def fij(sub="ij"):
    """``$F_{ij}$``."""

    return r"$F_{%s}$" % sub


def apply_style():
    """Set the matplotlib defaults of the suite (call once, before drawing)."""

    base = FIGURE["font_size"]
    mpl.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["DejaVu Sans", "Arial", "Helvetica"],
            "mathtext.fontset": "dejavusans",
            "font.size": base,
            "axes.titlesize": base + 0.5,
            "axes.labelsize": base,
            "xtick.labelsize": base - 1,
            "ytick.labelsize": base - 1,
            "legend.fontsize": base - 1.5,
            "legend.title_fontsize": base - 1,
            "figure.titlesize": base + 2,
            "axes.linewidth": 0.8,
            "lines.linewidth": 1.3,
            "lines.markersize": 4.5,
            "xtick.major.width": 0.8,
            "ytick.major.width": 0.8,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "legend.frameon": True,
            "legend.edgecolor": "0.3",
            "legend.framealpha": 0.92,
            "legend.fancybox": False,
            "axes.prop_cycle": mpl.cycler(
                color=[OKABE_ITO[k] for k in ("blue", "vermillion", "green", "orange", "purple", "sky", "black")]
            ),
            "image.cmap": "viridis",
            "savefig.dpi": FIGURE["dpi"],
            "savefig.bbox": "tight",
            "savefig.pad_inches": 0.04,
            "figure.dpi": 100,
        }
    )


def add_figure_args(parser):
    """``--dpi``, ``--fig-format``, ``--font-size``: figure settings for every script."""

    g = parser.add_argument_group("figures")
    g.add_argument("--dpi", type=int, default=FIGURE["dpi"], help="resolution of raster figures (default 300)")
    g.add_argument(
        "--fig-format",
        nargs="+",
        default=list(FIGURE["formats"]),
        choices=("png", "pdf", "svg"),
        help="file formats of every figure, the first is the one linked in index.md (default png)",
    )
    g.add_argument("--font-size", type=float, default=FIGURE["font_size"], help="base font size in pt (default 9)")
    return g


def configure(args):
    """Take the figure settings from parsed arguments and apply the style."""

    FIGURE["dpi"] = getattr(args, "dpi", FIGURE["dpi"])
    FIGURE["formats"] = tuple(getattr(args, "fig_format", FIGURE["formats"]))
    FIGURE["font_size"] = getattr(args, "font_size", FIGURE["font_size"])
    apply_style()


def style(ax, grid=True):
    """Inward ticks on all four sides and a light grid."""

    ax.tick_params(direction="in", which="both")
    ax.xaxis.set_ticks_position("both")
    ax.yaxis.set_ticks_position("both")
    if grid:
        ax.grid(color="0.85", ls="-", lw=0.5)
        ax.set_axisbelow(True)


def panel_label(ax, letter, x=-0.02, y=1.02):
    """A bold panel letter, (a), (b), ..., at the top left of ``ax``."""

    ax.text(x, y, "(%s)" % letter, transform=ax.transAxes, fontweight="bold", va="bottom", ha="right")


def legend_below(fig, handles, labels, ncol=2, y=0.0, **kw):
    """One legend for the whole figure, centred below the panels."""

    return fig.legend(handles, labels, loc="upper center", bbox_to_anchor=(0.5, y), ncol=ncol, **kw)


def scaled(values):
    """``(values / 10^k, k)`` with ``10^k`` the order of magnitude of ``max |values|``, to avoid offset labels."""

    peak = float(np.max(np.abs(values))) if np.size(values) else 0.0
    k = int(np.floor(np.log10(peak))) if peak > 0 else 0
    return np.asarray(values) / 10.0**k, k


def times_ten(k):
    return "" if k == 0 else r" ($\times 10^{%d}$)" % k


# --------------------------------------------------------------------------------------------- colours


def fit_family(spec):
    for family in ("parabola", "poly", "morse", "cubic", "spline"):
        if spec.startswith(family):
            return family
    return "other"


def fit_style(spec, specs=None):
    """``{"color", "ls"}`` of a fit: its family's colour, and a line style that tells variants of a family apart.

    ``specs`` (all fits drawn together) fixes the order of the line styles within a family.
    """

    family = fit_family(spec)
    color = FIT_FAMILY_COLORS.get(family, OKABE_ITO["gray"])
    members = [s for s in (specs or [spec]) if fit_family(s) == family] or [spec]
    index = members.index(spec) if spec in members else 0
    return {"color": color, "ls": FAMILY_STYLES[index % len(FAMILY_STYLES)]}


def fit_color(spec):
    return fit_style(spec)["color"]


def result_color(result):
    """A VOFResult's colour: its fit's for a non-default fd fit, otherwise its method's."""

    if result.fit is not None and result.fit != "parabola":
        return fit_color(result.fit)
    return METHOD_COLORS.get(result.method, "black")


def transition_label(direction):
    """``'$^3E$ $\\rightarrow$ $^1A_1$'`` style label of a direction, from the dataset's math-text state labels."""

    initial, final = ("es", "gs") if direction == "emission" else ("gs", "es")
    labels = ccd_data.STATE_LABELS
    return r"%s $\rightarrow$ %s" % (labels[initial], labels[final])


def ref_direction_label(result):
    return transition_label(result.direction)


def direction_title(direction):
    return "%s, %s" % (direction.capitalize(), transition_label(direction))


# --------------------------------------------------------------------------------------------- text helpers


def wrap_lines(lines, width=NOTE_WIDTH):
    out = []
    for line in lines:
        out.extend(textwrap.wrap(line, width, subsequent_indent="    ") or [""])
    return "\n".join(out)


def settings_text(params, ref, tested):
    """Common settings, then every method's own notes (for log.txt and index.md, not for figures)."""

    ni, nf = ref.solver.ni, ref.solver.nf
    lines = [
        "%s: omega_i = %.5f meV, omega_f = %.5f meV, dQ = %.6f amu^0.5 A, ni x nf = %d x %d"
        % (ref.direction, ref.omega_i, ref.omega_f, ref.solver.delta_q, ni, nf)
    ] + params.describe()[:2]
    for r in [ref] + [t for t, _ in tested]:
        lines += ["", "[%s]%s" % (r.label, " reference" if r is ref else "")] + r.notes
    return wrap_lines(lines)


def off_scale_note(ax, value, x=0.0, label="ZPL"):
    """Mark a peak clipped by the y limits with a small arrow and its height."""

    top = ax.get_ylim()[1]
    ax.annotate(
        "%s %.2g\n(off scale)" % (label, value),
        xy=(x, top),
        xytext=(10, -4),
        textcoords="offset points",
        ha="left",
        va="top",
        fontsize=FIGURE["font_size"] - 2,
        color="0.35",
        arrowprops=dict(arrowstyle="-|>", color="0.35", lw=0.6, shrinkA=0, shrinkB=0),
    )


# --------------------------------------------------------------------------------------------- shared figures


def draw_vof_comparison(fig, grid, col, ref, tested, tol_vof):
    """One column of a VOF comparison: overlay, signed error, log error. Returns the legend handles and labels."""

    top = fig.add_subplot(grid[0, col])
    mid = fig.add_subplot(grid[1, col], sharex=top)
    bottom = fig.add_subplot(grid[2, col], sharex=top)
    peak = np.max(ref.lineshape)
    color = DIRECTION_COLORS[ref.direction]

    top.fill_between(ref.axis, 0.0, ref.lineshape, color=color, alpha=0.18, linewidth=0)
    top.plot(ref.axis, ref.lineshape, color=color, lw=2.4, alpha=0.55, label="reference: %s" % ref.label)
    errors = [(tst.lineshape - ref.lineshape) / peak for tst, _ in tested]
    scaled_err, k = scaled(np.concatenate(errors)) if errors else (np.zeros(0), 0)
    worst = 1e-30
    for n, ((tst, stats), err) in enumerate(zip(tested, errors)):
        c = result_color(tst)
        ls = TEST_STYLES[n % len(TEST_STYLES)]
        top.plot(tst.axis, tst.lineshape, color=c, ls=ls, lw=1.2, label=tst.label)
        mid.plot(ref.axis, err / 10.0**k, color=c, lw=0.9, ls="-")
        bottom.semilogy(ref.axis, np.maximum(np.abs(err), 1e-18), color=c, lw=0.9)
        worst = max(worst, stats["vof_max"])
    if tol_vof is not None:
        bottom.axhline(tol_vof, color="0.3", ls="--", lw=0.9, label="tolerance (%.0e)" % tol_vof)
    mid.axhline(0.0, color="0.6", lw=0.6)
    if k:
        mid.text(
            0.01,
            0.97,
            r"$\times 10^{%d}$" % k,
            transform=mid.transAxes,
            ha="left",
            va="top",
            fontsize=FIGURE["font_size"] - 1,
        )

    sideband = max(ref.sideband_max, 1e-30)
    top.set_ylim(0.0, 1.3 * sideband)
    if peak > 1.3 * sideband:
        off_scale_note(top, peak)
    top.set_title(direction_title(ref.direction))
    if col == 0:
        top.set_ylabel(VOF_LABEL)
        mid.set_ylabel(r"$\Delta$VOF / peak")
        bottom.set_ylabel(r"$|\Delta$VOF$|$ / peak")
    bottom.set_ylim(1e-18, max(10 * (tol_vof or 0.0), 10 * worst, 1e-16))
    bottom.set_xlabel(E_ZPL_LABEL)
    for a in (top, mid, bottom):
        style(a)
    for a in (top, mid):
        a.tick_params(labelbottom=False)
    handles, labels = top.get_legend_handles_labels()
    h2, l2 = bottom.get_legend_handles_labels()
    return handles + h2, labels + l2


def vof_comparison_figure(plt, comparison, params, tol_vof, title, with_settings=False):
    """``comparison`` = ``{direction: (ref, [(test, stats), ...])}``; one column per direction.

    The settings and statistics of every curve are written to log.txt / index.md / the CSVs, not drawn
    (``with_settings`` is accepted for compatibility and ignored). ``tol_vof = None`` draws no tolerance line.
    """

    fig = plt.figure(figsize=(7.2, 6.4))
    grid = fig.add_gridspec(3, 2, height_ratios=(2.2, 1, 1), hspace=0.08, wspace=0.24, top=0.9)
    legend = None
    for col, (ref, tested) in enumerate(comparison.values()):
        handles, labels = draw_vof_comparison(fig, grid, col, ref, tested, tol_vof)
        legend = legend or (handles, labels)
    if legend:
        legend_below(fig, *legend, ncol=2, y=0.02)
    fig.suptitle(title, y=0.965)
    return fig


def fc_error_figure(plt, comparison, title):
    """log10 |F_ij(test) - F_ij(ref)|: one row per test, one column per direction, maximum marked; one colour bar."""

    n_tests = len(next(iter(comparison.values()))[1])
    fig, axes = plt.subplots(n_tests, 2, figsize=(7.0, 2.9 * n_tests + 0.6), squeeze=False, constrained_layout=True)
    images = []
    vmax = -6.0
    devs = {}
    for col, (direction, (ref, tested)) in enumerate(comparison.items()):
        for row, (tst, stats) in enumerate(tested):
            dev = np.log10(np.maximum(np.abs(tst.solver.fc_ints - ref.solver.fc_ints), 1e-18))
            devs[row, col] = dev
            vmax = max(vmax, float(dev.max()))
    for col, (direction, (ref, tested)) in enumerate(comparison.items()):
        for row, (tst, stats) in enumerate(tested):
            ax = axes[row, col]
            im = ax.imshow(devs[row, col], origin="lower", aspect="auto", cmap="viridis", vmin=-17, vmax=vmax)
            images.append(im)
            i, j = stats["fc_at"]
            ax.plot(j, i, marker="x", color=OKABE_ITO["vermillion"], ms=6, mew=1.5, clip_on=False)
            ax.set_title("%s: %s" % (direction.capitalize(), tst.label), fontsize=FIGURE["font_size"] - 0.5)
            ax.set_xlabel(r"$j$ (final level)")
            ax.set_ylabel(r"$i$ (initial level)")
            style(ax, grid=False)
    fig.colorbar(
        images[0], ax=axes, label=r"$\log_{10}|\Delta F_{ij}|$", extend="min", shrink=0.9 if n_tests > 1 else 1.0
    )
    fig.suptitle(title)
    return fig


apply_style()
