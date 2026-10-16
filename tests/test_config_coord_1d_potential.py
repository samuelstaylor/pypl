"""Fits of a sampled potential energy curve: parabola, polynomial, Morse, cubic and splines."""

import numpy as np
import pytest
from scipy.interpolate import CubicSpline

from pypl.config_coord_1d_potential import SPLINE_ORDER_NAMES, config_coord_1d_potential


def fit(diamond, state, method, *args, **kw):
    d = diamond[state]
    pot = config_coord_1d_potential(d["q"], d["e"], q0=d["center"])
    getattr(pot, method)(*args, **kw)
    return pot


def test_parabola_reproduces_tutorial_002(diamond):
    gs, es = diamond["gs"]["parabola"], diamond["es"]["parabola"]
    assert gs.kind == "parabola" and gs.bounded
    assert gs.phonon_energy() == pytest.approx(59.08890, abs=1e-5)
    assert es.phonon_energy() == pytest.approx(64.31191, abs=1e-5)
    np.testing.assert_array_equal(gs.index, np.arange(5))
    assert gs(np.array([0.0]))[0] == pytest.approx(gs.params["c"])


def test_not_fitted_raises():
    pot = config_coord_1d_potential([0.0, 1.0], [0.0, 1.0])
    assert pot.describe() == "not fitted"
    with pytest.raises(RuntimeError):
        pot(0.0)


def test_polynomials(diamond):
    p2 = fit(diamond, "gs", "fit_polynomial", 2)
    assert p2.bounded and p2.rms() < diamond["gs"]["parabola"].rms()
    p4 = fit(diamond, "gs", "fit_polynomial", 4)
    assert not p4.bounded and "UNBOUNDED" in p4.describe()
    p3 = fit(diamond, "gs", "fit_polynomial", 3)
    assert not p3.bounded and "odd degree" in p3.describe()


def test_morse(diamond):
    m = fit(diamond, "es", "fit_morse")
    assert m.kind == "morse" and m.bounded and m.params["D"] > 0
    assert m.rms() < 1e-3
    assert m.phonon_energy() == pytest.approx(diamond["es"]["phonon"], rel=0.1)


def test_cubic(diamond):
    a = diamond["gs"]["parabola"].params["a"]
    c = fit(diamond, "gs", "fit_cubic", curvature=a)
    assert c.params["a"] == a and c.params["a_fixed"]
    assert c.phonon_energy() == pytest.approx(diamond["gs"]["phonon"], abs=1e-12)
    assert c.params["turnover"] == pytest.approx(-2 * a / (3 * c.params["b"]))
    assert not c.bounded and "turnover" in c.describe()
    free = fit(diamond, "gs", "fit_cubic")
    assert not free.params["a_fixed"]


def test_default_spline_is_natural_cubic(diamond):
    s = fit(diamond, "gs", "fit_spline")
    q, e = diamond["gs"]["q"], diamond["gs"]["e"]
    x = np.linspace(q.min(), q.max(), 2001)
    assert np.max(np.abs(s(x) - CubicSpline(q, e, bc_type="natural")(x))) < 1e-12
    assert s.params["bc"] == "natural" and s.params["tails"] == "harmonic"


@pytest.mark.parametrize(
    "order, bc",
    [(3, "not-a-knot"), (3, "clamped"), (4, "not-a-knot"), (5, "natural"), (5, "clamped"), (2, "not-a-knot")],
)
def test_interpolating_splines_pass_through_samples(diamond, order, bc):
    s = fit(diamond, "gs", "fit_spline", order=order, bc=bc)
    assert np.max(np.abs(s.residuals())) < 1e-12
    assert SPLINE_ORDER_NAMES[order] in s.describe()


def test_clamped_junctions(diamond):
    assert fit(diamond, "gs", "fit_spline", order=3, bc="clamped").params["junction"] == "C1"
    assert fit(diamond, "gs", "fit_spline", order=5, bc="clamped").params["junction"] == "C2"
    assert fit(diamond, "gs", "fit_spline", order=3, bc="natural").params["junction"] == "C0"


def test_smoothing_spline(diamond):
    s = fit(diamond, "gs", "fit_spline", smoothing=5e-5)
    # FITPACK meets its sum-of-squares target to a relative 1e-3.
    assert s.rms() <= 5e-5 * np.sqrt(1.001)


@pytest.mark.parametrize(
    "kw",
    [
        dict(order=4, bc="natural"),
        dict(order=2, bc="clamped"),
        dict(order=6),
        dict(bc="periodic"),
        dict(smoothing=1e-4, bc="natural"),
        dict(smoothing=-1.0),
        dict(tails="linear"),
        dict(bc="clamped", tails="smooth"),
    ],
)
def test_spline_rejects_invalid_options(diamond, kw):
    with pytest.raises(ValueError):
        fit(diamond, "gs", "fit_spline", **kw)


def test_smooth_tails_are_c2_and_leave_the_inside_unchanged(diamond):
    smooth = fit(diamond, "gs", "fit_spline", order=3, bc="not-a-knot", tails="smooth")
    plain = fit(diamond, "gs", "fit_spline", order=3, bc="not-a-knot")
    q = diamond["gs"]["q"]
    inside = np.linspace(q.min(), q.max(), 4001)
    assert np.array_equal(smooth(inside), plain(inside))
    assert smooth.params["junction"] == "C2" and smooth.bounded

    f = lambda x: float(smooth(np.array([x]))[0])  # noqa: E731
    h = 1e-4
    for end in (q.min(), q.max()):
        slope_jump = (f(end + h) - f(end)) / h - (f(end) - f(end - h)) / h
        curv_jump = (f(end + 2 * h) - 2 * f(end + h) + f(end)) / h**2 - (
            f(end) - 2 * f(end - h) + f(end - 2 * h)
        ) / h**2
        assert abs(slope_jump) < 1e-3 and abs(curv_jump) * h < 1e-3

    # Beyond the data the tail is the spline's own second-order Taylor expansion.
    hi, (v1, v2) = q.max(), (smooth.params["end_slopes"][1], smooth.params["end_curvatures"][1])
    x = hi + 0.5
    assert f(x) == pytest.approx(smooth.params["shift_hi"] + v1 * 0.5 + 0.5 * v2 * 0.25, rel=1e-12)


def test_smooth_tails_unbounded_when_end_curvature_negative():
    q = np.linspace(-1, 1, 9)
    pot = config_coord_1d_potential(q, q**2 - 0.8 * q**4)  # V'' < 0 at both ends
    pot.fit_spline(bc="not-a-knot", tails="smooth")
    assert not pot.bounded and "UNBOUNDED" in pot.describe()


def test_minimum_and_phonon_energy_of_spline(diamond):
    s = fit(diamond, "gs", "fit_spline", order=3, bc="not-a-knot", tails="smooth")
    q_min, e_min = s.minimum()
    assert abs(q_min) < 0.05 and e_min < 1e-3
    assert s.phonon_energy() == pytest.approx(diamond["gs"]["phonon"], rel=0.06)
