"""Potential fits of the CCD wells, used by the finite-difference method.

A fit is named by a short spec string:

    parabola   a (Q - Q0)^2 + c fitted to the images near the minimum (GS 1-5, ES 5-11). This is the
               curve the solver's omega comes from, so it is the default and the one to test with.
    polyN      even polynomial of degree N through every image (e.g. poly2, poly6, poly8)
    morse      Morse well through every image
    spline     spline through every image, with the parabola's curvature as harmonic tails
    cubic      a (Q - Q0)^2 + b (Q - Q0)^3 + c with a fixed to the parabola's, b and c fitted to every
               image: the potential of the cubic perturbation method. ``cubic-near`` fits only the
               parabola's images, ``cubic-free`` fits a too; they combine as ``cubic-near-free``.
               A cubic is unbounded beyond its turnover, so fd accepts it only when the turnover lies
               outside the grid (a bound well in the box).

A spline spec can set the order, the end conditions or a smoothing target:

    spline[K][-natural | -notaknot | -clamped | -sRMS][-smooth]

    K           polynomial degree of the pieces: 2, 3 (default), 4 or 5
    -natural    vanishing higher derivatives at the end images (default for odd K; odd K only)
    -notaknot   no end conditions imposed (default for even K; any K)
    -clamped    slope (and for K = 5 curvature) matched to the harmonic tails, a smooth join (odd K only)
    -sRMS       smoothing spline with rms residual <= RMS meV instead of exact interpolation
    -smooth     smooth (C2) tails: the spline's own value, slope and curvature continued past the end images,
                instead of the harmonic tails; nothing inside the data changes (not with -clamped)

e.g. ``spline`` (cubic, natural), ``spline5-clamped``, ``spline4``, ``spline3-s0.05``,
``spline3-notaknot-smooth``. The command-line flags ``--spline-order``, ``--spline-bc``, ``--spline-smoothing``
and ``--spline-tails`` build the same spec from ``--fd-fit spline``. ``canonical(spec)`` gives the one spelling used in run names and labels.

``make_potential(ccd, state, spec)`` returns a fitted ``config_coord_1d_potential``. To add a fit,
add a ``config_coord_1d_potential.fit_*`` method and a branch here, plus an entry in ``FIT_HELP``.
"""

import re

from harness import ccd_data
from pypl.config_coord_1d_potential import SPLINE_ORDER_NAMES, config_coord_1d_potential

DEFAULT_FIT = "parabola"
# Fits compared by vof_report.py. poly4 and poly8 are included on purpose: they are unbounded here.
REPORT_FITS = (
    "parabola",
    "poly2",
    "poly4",
    "poly6",
    "poly8",
    "morse",
    "spline",
    "spline3-notaknot",
    "spline3-clamped",
    "spline4-notaknot",
    "spline5-clamped",
    "spline3-s0.05",
)
FIT_HELP = {
    "parabola": "parabola near each minimum, the solver's own wells (default; use for testing)",
    "polyN": "even polynomial of degree N through every image, e.g. poly6",
    "morse": "Morse well through every image",
    "cubic[-near][-free]": "cubic a Q^2 + b Q^3 + c (a fixed to the parabola unless -free; every image unless "
    "-near); fd accepts it only if its turnover lies outside the grid",
    "spline[K][-natural|-notaknot|-clamped|-sRMS][-smooth]": "spline of degree K (2-5, default 3) through every "
    "image with harmonic tails (-smooth: C2 tails continuing the spline itself); natural / not-a-knot / "
    "tail-clamped ends, or a smoothing spline with rms <= RMS meV",
}
SPLINE_BC = {"natural": "natural", "notaknot": "not-a-knot", "clamped": "clamped"}
_SPLINE_RE = re.compile(r"spline([2-5])?(?:-(natural|notaknot|clamped)|-s(\d+(?:\.\d*)?|\.\d+))?(-smooth)?")


def parse_spec(spec):
    """``(kind, options)`` for a spec string; raises ValueError on an unknown or impossible one.

    ``options`` is ``{}`` for parabola and morse, ``{"degree": N}`` for a polynomial, and
    ``{"order", "bc", "smoothing_mev", "tails"}`` for a spline, with the defaults filled in.
    """

    if spec in ("parabola", "morse"):
        return spec, {}
    match = re.fullmatch(r"cubic(-near)?(-free)?", spec)
    if match:
        return "cubic", {"points": "parabola" if match.group(1) else "all", "fixed_a": not match.group(2)}
    match = re.fullmatch(r"poly(\d+)", spec)
    if match:
        return "polynomial", {"degree": int(match.group(1))}
    match = _SPLINE_RE.fullmatch(spec)
    if match:
        order = int(match.group(1) or 3)
        smoothing = float(match.group(3)) if match.group(3) is not None else None
        bc = None if smoothing is not None else SPLINE_BC[match.group(2) or ("natural" if order % 2 else "notaknot")]
        if bc in ("natural", "clamped") and order % 2 == 0:
            raise ValueError("fit %s: '%s' end conditions need an odd spline order; use -notaknot" % (spec, bc))
        tails = "smooth" if match.group(4) else "harmonic"
        if tails == "smooth" and bc == "clamped":
            raise ValueError(
                "fit %s: -clamped matches the harmonic tails; -smooth tails follow the spline instead" % spec
            )
        return "spline", {"order": order, "bc": bc, "smoothing_mev": smoothing, "tails": tails}
    raise ValueError(
        "unknown fit %r; use parabola, polyN (e.g. poly6), morse, or spline[K][-natural|-notaknot|-clamped|-sRMS]"
        % spec
    )


def spline_spec(order=3, bc=None, smoothing_mev=None, tails="harmonic"):
    """The canonical spec of a spline: ``spline`` for the default cubic natural one with harmonic tails."""

    if tails not in ("harmonic", "smooth"):
        raise ValueError("spline tails must be 'harmonic' or 'smooth', got %r" % (tails,))
    suffix = "-smooth" if tails == "smooth" else ""
    if smoothing_mev is not None:
        if bc is not None:
            raise ValueError("end conditions do not apply to a smoothing spline; give a bc or a smoothing, not both")
        return "spline%d-s%g%s" % (order, smoothing_mev, suffix)
    bc = bc or ("natural" if order % 2 else "not-a-knot")
    if order == 3 and bc == "natural" and not suffix:
        return "spline"
    spec = "spline%d-%s%s" % (order, {v: k for k, v in SPLINE_BC.items()}[bc], suffix)
    parse_spec(spec)
    return spec


def canonical(spec):
    """One spelling per fit, e.g. ``spline3-natural`` -> ``spline`` and ``spline4`` -> ``spline4-notaknot``."""

    kind, options = parse_spec(spec)
    return spline_spec(**options) if kind == "spline" else spec


def label(spec):
    """Human-readable name of a fit spec, for legends and tables."""

    kind, options = parse_spec(spec)
    if kind == "spline":
        name = SPLINE_ORDER_NAMES[options["order"]]
        tails = ", smooth C2 tails" if options["tails"] == "smooth" else ""
        if options["smoothing_mev"] is not None:
            return "%s smoothing spline (rms %g meV%s)" % (name, options["smoothing_mev"], tails)
        return "%s spline (%s%s)" % (name, options["bc"], tails)
    if kind == "polynomial":
        return "polynomial deg. %d" % options["degree"]
    if kind == "cubic":
        return "cubic (%s, %s)" % (
            "a fixed" if options["fixed_a"] else "a fitted",
            "near-min images" if options["points"] == "parabola" else ccd_data.active().cubic_images or "all images",
        )
    return {"parabola": "parabola (near min.)", "morse": "Morse"}[kind]


def images(ccd, state, spec):
    """Which images of ``state`` the fit ``spec`` uses, as printable text (see ``ccd_data.image_note``)."""

    kind, options = parse_spec(spec)
    if kind == "parabola" or (kind == "cubic" and options["points"] == "parabola"):
        return ccd_data.image_note(ccd, state, ccd["%s_potential" % state].index)
    if kind == "cubic":
        return ccd_data.image_note(ccd, state, ccd.get("%s_cubic_index" % state))
    return ccd_data.image_note(ccd, state)


def state_label(ccd, state, spec):
    """``label(spec)`` plus the images it used for ``state``, for legends: ``parabola (near min.) [images ...]``."""

    return "%s [%s]" % (label(spec), images(ccd, state, spec))


def make_potential(ccd, state, spec=DEFAULT_FIT):
    """A ``config_coord_1d_potential`` of ``state`` fitted as ``spec`` (it may be unbounded; check ``.bounded``)."""

    kind, options = parse_spec(spec)
    if kind == "parabola":
        return ccd["%s_potential" % state]
    q, e = ccd_data.samples(ccd, state)
    potential = config_coord_1d_potential(q, e, ccd_data.center(ccd, state))
    if kind == "polynomial":
        potential.fit_polynomial(options["degree"])
    elif kind == "cubic":
        index = (
            ccd["%s_potential" % state].index if options["points"] == "parabola" else ccd.get("%s_cubic_index" % state)
        )
        potential.fit_cubic(index=index, curvature=ccd["%s_params" % state][0] if options["fixed_a"] else None)
    elif kind == "morse":
        potential.fit_morse()
    else:
        smoothing = options["smoothing_mev"] * 1e-3 if options["smoothing_mev"] is not None else None
        potential.fit_spline(
            tail_curvature=ccd["%s_params" % state][0],
            order=options["order"],
            bc=options["bc"],
            smoothing=smoothing,
            tails=options["tails"],
        )
    return potential


def check_bounded(ccd, spec, q_range=None):
    """Raise ValueError if ``spec`` gives either well no bound states.

    A cubic counts as bound inside a finite grid ``q_range = (lo, hi)`` if its turnover lies outside it.
    """

    for state in ccd_data.STATES:
        potential = make_potential(ccd, state, spec)
        if potential.kind == "cubic" and q_range is not None:
            turnover = potential.params["turnover"]
            if not q_range[0] < turnover < q_range[1]:
                continue
            raise ValueError(
                "the %s cubic turns over at Q = %.2f, inside the grid %.2f..%.2f: no bound states in the box"
                % (ccd_data.STATE_NAMES[state], turnover, *q_range)
            )
        if not potential.bounded:
            raise ValueError(
                "fit %s of the %s images is unbounded below and has no bound states (%s)"
                % (spec, ccd_data.STATE_NAMES[state], potential.describe())
            )
