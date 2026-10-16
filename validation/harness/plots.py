"""Shared plot style and the figures more than one script draws."""

import textwrap

import numpy as np

from harness import ccd_data

DIRECTION_COLORS = {"emission": "#DB4437", "absorption": "#4285F4"}
# The colour of each well follows the active dataset (diamond: GS red, ES blue; MgO: 3E red, 1A1 blue).
STATE_COLORS = ccd_data.STATE_COLORS
# One colour per F_ij method and per fit, so a method or fit looks the same in every figure.
METHOD_COLORS = {"pypl": "#555555", "hermite": "black", "fd": "#0F9D58", "cubic": "#8E24AA"}
FIT_COLORS = {
    "parabola": "black",
    "poly2": "#9E9E9E",
    "poly4": "#AB47BC",
    "poly6": "#F4B400",
    "poly8": "#FF7043",
    "morse": "#4285F4",
    # The cubic fits of the perturbation method in purples.
    "cubic": "#8E24AA",
    "cubic-near": "#BA68C8",
    "cubic-free": "#4A148C",
    "cubic-near-free": "#CE93D8",
    # Splines in greens and teals, so the variants read as one family.
    "spline": "#0F9D58",
    "spline3-notaknot": "#00ACC1",
    "spline3-clamped": "#1B5E20",
    "spline4-notaknot": "#9CCC65",
    "spline5-clamped": "#00897B",
    "spline3-s0.05": "#E91E63",
}
# Fits outside the table cycle through these.
EXTRA_COLORS = ("#795548", "#8D6E63", "#5C6BC0", "#26A69A", "#D4E157", "#EC407A")
TEST_STYLES = ("--", ":", "-.", (0, (1, 3)))
NOTE_WIDTH = 118


def fit_color(spec):
    if spec in FIT_COLORS:
        return FIT_COLORS[spec]
    return EXTRA_COLORS[sum(map(ord, spec)) % len(EXTRA_COLORS)]


def style(ax, grid=True):
    """Inward ticks on all four sides and a light dashed grid."""

    ax.tick_params(direction="in")
    ax.xaxis.set_ticks_position("both")
    ax.yaxis.set_ticks_position("both")
    if grid:
        ax.grid(color="gray", ls="--", lw=0.4, alpha=0.5)


def wrap_lines(lines, width=NOTE_WIDTH):
    out = []
    for line in lines:
        out.extend(textwrap.wrap(line, width, subsequent_indent="    ") or [""])
    return "\n".join(out)


def settings_text(params, ref, tested):
    """Common settings, then every method's own notes, for a figure's settings panel."""

    ni, nf = ref.solver.ni, ref.solver.nf
    lines = [
        "%s: omega_i = %.5f meV, omega_f = %.5f meV, dQ = %.6f amu^0.5 A, ni x nf = %d x %d"
        % (ref.direction, ref.omega_i, ref.omega_f, ref.solver.delta_q, ni, nf)
    ] + params.describe()[:2]
    for r in [ref] + [t for t, _ in tested]:
        lines += ["", "[%s]%s" % (r.label, " reference" if r is ref else "")] + r.notes
    return wrap_lines(lines)


def draw_vof_comparison(fig, grid, col, ref, tested, tol_vof, text):
    """One column of a VOF comparison: overlay, signed error, log error, settings panel."""

    top = fig.add_subplot(grid[0, col])
    mid = fig.add_subplot(grid[1, col], sharex=top)
    bottom = fig.add_subplot(grid[2, col], sharex=top)
    peak = np.max(ref.lineshape)
    color = DIRECTION_COLORS[ref.direction]

    top.fill_between(ref.axis, 0.0, ref.lineshape, color=color, alpha=0.2, linewidth=0)
    top.plot(ref.axis, ref.lineshape, color=color, lw=2.6, alpha=0.6, label=ref.label + " (reference)")
    worst = 1e-30
    for k, (tst, stats) in enumerate(tested):
        c = result_color(tst)
        ls = TEST_STYLES[k % len(TEST_STYLES)]
        top.plot(tst.axis, tst.lineshape, color=c, ls=ls, lw=1.1, label=tst.label)
        err = (tst.lineshape - ref.lineshape) / peak
        mid.plot(ref.axis, err, color=c, lw=0.9, label="%s (max %.1e)" % (tst.label, stats["vof_max"]))
        bottom.semilogy(ref.axis, np.maximum(np.abs(err), 1e-18), color=c, lw=0.8, label=tst.label)
        worst = max(worst, stats["vof_max"])
    bottom.axhline(tol_vof, color="gray", ls="--", lw=0.9, label="tolerance %.0e" % tol_vof)

    sideband = max(ref.sideband_max, 1e-30)
    top.set_ylim(0.0, 1.35 * sideband)
    top.annotate(
        "ZPL peak %.3g (off scale)" % peak,
        xy=(0.0, 1.35 * sideband),
        xytext=(8, -12),
        textcoords="offset points",
        fontsize=7,
        color="dimgray",
    )
    top.set_title("%s (%s)" % (ref.direction, ref_direction_label(ref)))
    top.set_ylabel("VOF (1/meV)")
    mid.set_ylabel("(test - ref) / peak")
    bottom.set_ylabel("|test - ref| / peak")
    bottom.set_ylim(1e-18, max(10 * tol_vof, 10 * worst))
    bottom.set_xlabel("E - E$_{ZPL}$ (meV)")
    for a in (top, mid, bottom):
        a.legend(edgecolor="black", fontsize=7, loc="upper right")
        style(a)
    for a in (top, mid):
        a.tick_params(labelbottom=False)
    if text is not None:
        text_ax = fig.add_subplot(grid[3, col])
        text_ax.axis("off")
        text_ax.text(
            0.0, 0.95, text, family="monospace", fontsize=5.6, va="top", ha="left", transform=text_ax.transAxes
        )


def ref_direction_label(result):
    return ccd_data.DIRECTION_LABELS[result.direction].replace("->", "to")


def result_color(result):
    """A VOFResult's colour: its fit's for a non-default fd fit, otherwise its method's."""

    if result.fit is not None and result.fit != "parabola":
        return fit_color(result.fit)
    return METHOD_COLORS.get(result.method, "black")


def vof_comparison_figure(plt, comparison, params, tol_vof, title, with_settings=True):
    """``comparison`` = ``{direction: (ref, [(test, stats), ...])}``; one column per direction."""

    height = (3, 1.3, 1.3, 2.6) if with_settings else (3, 1.3, 1.3)
    fig = plt.figure(figsize=(15, 15.5 if with_settings else 10))
    grid = fig.add_gridspec(len(height), 2, height_ratios=height, hspace=0.08, wspace=0.18, top=0.93)
    for col, (ref, tested) in enumerate(comparison.values()):
        text = settings_text(params, ref, tested) if with_settings else None
        draw_vof_comparison(fig, grid, col, ref, tested, tol_vof, text)
    fig.suptitle(title, fontsize=13, y=0.965)
    return fig


def fc_error_figure(plt, comparison, title):
    """log10 |F_ij(test) - F_ij(ref)|: one row per test, one column per direction, maximum marked."""

    n_tests = len(next(iter(comparison.values()))[1])
    fig, axes = plt.subplots(n_tests, 2, figsize=(10, 3.9 * n_tests), squeeze=False)
    for col, (direction, (ref, tested)) in enumerate(comparison.items()):
        for row, (tst, stats) in enumerate(tested):
            ax = axes[row, col]
            dev = np.log10(np.maximum(np.abs(tst.solver.fc_ints - ref.solver.fc_ints), 1e-18))
            im = ax.imshow(dev, origin="lower", aspect="auto", cmap="viridis", vmin=-17, vmax=max(-6, dev.max()))
            ax.plot(*stats["fc_at"][::-1], marker="x", color="red", ms=7)
            ax.set_title(
                "%s: %s vs %s\nmax %.2e at (i, j) = %s"
                % (direction, tst.label, ref.label, stats["fc_max"], stats["fc_at"]),
                fontsize=9,
            )
            ax.set_xlabel("j (final level)")
            ax.set_ylabel("i (initial level)")
            fig.colorbar(im, ax=ax, label=r"log$_{10}$ |$\Delta F_{ij}$|")
    fig.suptitle(title, fontsize=12)
    fig.tight_layout()
    return fig
