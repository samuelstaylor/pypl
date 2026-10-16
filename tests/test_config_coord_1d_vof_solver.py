"""The extended 1D CCD solver: extension contract, hermite / fd / cubic routes and line positions from levels."""

import numpy as np
import pytest
from scipy import constants

from pypl.config_coord_1d_potential import config_coord_1d_potential
from pypl.config_coord_1d_solver import config_coord_1d_solver
from pypl.config_coord_1d_vof_solver import KINETIC, config_coord_1d_vof_solver, grid_spacing

AXIS = np.linspace(-200.0, 1200.0, 1401)
LSP = dict(eneaxis=AXIS, temp=5.0, sigma=[8.0, 25.0], zpl_lorentzian=True, gamma=0.3)
SUB = config_coord_1d_vof_solver


@pytest.fixture(scope="module")
def pair(diamond):
    """Emission on tutorial 002: initial = ES (at delta_q), final = GS (at 0)."""

    return diamond["es"]["phonon"], diamond["gs"]["phonon"], diamond["delta_q"]


@pytest.fixture(scope="module")
def ruhoff20(pair):
    s = config_coord_1d_solver(*pair)
    s.compute_franck_condon_integrals(ni=20, nf=20)
    return s.fc_ints


# ----------------------------------------------------------------------------- extension contract


def test_subclass_overrides_nothing():
    base = config_coord_1d_solver
    assert issubclass(SUB, base)
    assert [k for k in vars(SUB) if k in vars(base) and not k.startswith("__")] == []


def test_inherited_methods_bit_identical(pair):
    old, new = config_coord_1d_solver(*pair), SUB(*pair)
    for s in (old, new):
        s.compute_franck_condon_integrals(ni=50, nf=60)
        s.bulid_fc_lsp(**LSP)
    assert np.array_equal(old.fc_ints, new.fc_ints)
    assert np.array_equal(old.fc_lineshape, new.fc_lineshape)
    for a, b in zip(old.compute_spectrum(zpl=1945.0), new.compute_spectrum(zpl=1945.0)):
        assert np.array_equal(a, b)


# ----------------------------------------------------------------------------- harmonic functions


def test_recurrence_matches_explicit_hermite_formula():
    x = np.linspace(-6, 6, 2001)
    alpha = SUB.harmonic_alpha(59.0)
    xi = np.sqrt(alpha) * x
    explicit = np.array(
        [
            (alpha / np.pi) ** 0.25
            / np.sqrt(2.0**n * float(np.prod(np.arange(1, n + 1))))
            * np.polynomial.hermite.hermval(xi, np.eye(21)[n])
            * np.exp(-(xi**2) / 2)
            for n in range(21)
        ]
    )
    chi = SUB.harmonic_wavefunctions(x, 20, 59.0)
    assert np.max(np.abs(chi - explicit)) / np.max(np.abs(explicit)) < 1e-12
    assert np.array_equal(SUB.harmonic_wavefunction(x, 7, 59.0), chi[7])


def test_harmonic_functions_orthonormal_at_large_n():
    q = np.linspace(-30, 30, 200001)
    chi = SUB.harmonic_wavefunctions(q, 300, 43.8)
    gram = (chi * grid_spacing(q)) @ chi.T
    assert np.max(np.abs(gram - np.eye(301))) < 1e-9


def test_huang_rhys_and_turning_point():
    s = SUB.huang_rhys(59.0, 0.7)
    omega = 59.0e-3 * constants.eV / constants.hbar
    amu = constants.physical_constants["atomic mass constant"][0]
    assert s == pytest.approx(omega * (0.7e-10) ** 2 * amu / (2 * constants.hbar), rel=1e-12)
    # E_n = (n + 1/2) hbar w = alpha x^2 hbar w / 2 at the turning point.
    alpha = SUB.harmonic_alpha(59.0)
    assert alpha * SUB.classical_turning_point(4, 59.0) ** 2 == pytest.approx(9.0, rel=1e-12)


def test_kinetic_prefactor():
    amu = constants.physical_constants["atomic mass constant"][0]
    assert KINETIC == pytest.approx(constants.hbar**2 / (2 * amu * 1e-20) / constants.eV, rel=1e-12)


# ----------------------------------------------------------------------------- F_ij routes


def test_hermite_matches_ruhoff(pair, ruhoff20):
    s = SUB(*pair)
    s.compute_franck_condon_integrals_hermite(ni=20, nf=20)
    assert s.fc_ints.shape == (20, 20) and (s.ni, s.nf) == (20, 20)
    assert np.max(np.abs(s.fc_ints - ruhoff20)) < 1e-10


def test_hermite_stable_for_large_huang_rhys():
    s = SUB(74.07, 43.77, 3.0979)  # S ~ 50
    s.compute_franck_condon_integrals_hermite(2, 200)
    assert np.all(np.isfinite(s.fc_ints))
    assert np.sum(s.fc_ints[0] ** 2) == pytest.approx(1.0, abs=1e-10)


def test_fd_on_parabolas_matches_ruhoff(pair, ruhoff20, diamond):
    s = SUB(*pair)
    s.compute_franck_condon_integrals_fd(20, 20, diamond["es"]["parabola"], diamond["gs"]["parabola"], q_i=pair[2])
    assert np.max(np.abs(s.fc_ints - ruhoff20)) < 1e-8
    info = s.fd_info
    assert info["wall"] < 1e-12 and info["orthonormality"] < 1e-10
    assert info["ladder_i"] < 1e-4 and info["ladder_f"] < 1e-4
    assert s.fd_energies_i.shape == (20,) and s.fd_wavefunctions_f.shape == (20, info["npts"])


def test_fd_three_point_stencil_is_less_accurate(pair, ruhoff20, diamond):
    s = SUB(*pair)
    pots = (diamond["es"]["parabola"], diamond["gs"]["parabola"])
    s.compute_franck_condon_integrals_fd(20, 20, *pots, q_i=pair[2], stencil=3)
    err3 = np.max(np.abs(s.fc_ints - ruhoff20))
    assert 1e-8 < err3 < 1e-3


def test_fd_is_orientation_independent(pair, diamond):
    """Mirroring the CCD (initial well on the left instead of the right) must not change F_ij."""

    w_i, w_f, dq = pair
    ref = SUB(*pair)
    ref.compute_franck_condon_integrals_hermite(12, 12)
    mirrored_i = lambda q: diamond["es"]["parabola"](dq - q)  # noqa: E731  initial minimum now at 0
    mirrored_f = lambda q: diamond["gs"]["parabola"](dq - q)  # noqa: E731  final minimum now at dq
    s = SUB(*pair)
    s.compute_franck_condon_integrals_fd(12, 12, mirrored_i, mirrored_f, q_i=0.0, q_f=dq)
    assert np.max(np.abs(s.fc_ints - ref.fc_ints)) < 1e-7


def test_solve_schrodinger_fd_harmonic_ladder():
    omega = 60.0
    a = (omega * 1e-3) ** 2 / (4 * KINETIC)  # V = a Q^2 has hbar w = 2 sqrt(a KINETIC)
    q = np.linspace(-4, 4, 8001)
    energies, psi = SUB.solve_schrodinger_fd(a * q**2, q, 6)
    np.testing.assert_allclose(energies * 1e3, omega * (np.arange(6) + 0.5), atol=1e-6)
    assert np.allclose((psi * grid_spacing(q)) @ psi.T, np.eye(6), atol=1e-10)
    # Signed like the Hermite functions, which are the eigenstates of this well.
    chi = SUB.harmonic_wavefunctions(q, 5, omega)
    overlaps = np.sum(psi * chi, axis=1) * grid_spacing(q)
    np.testing.assert_allclose(overlaps, 1.0, atol=1e-8)


# ----------------------------------------------------------------------------- cubic correction


@pytest.fixture(scope="module")
def cubics(diamond):
    out = {}
    for state in ("gs", "es"):
        d = diamond[state]
        pot = config_coord_1d_potential(d["q"], d["e"], q0=d["center"])
        pot.fit_cubic(curvature=d["parabola"].params["a"])
        out[state] = pot
    return out


def test_anharmonic_strength_matches_si(diamond, cubics):
    w, b = diamond["gs"]["phonon"], cubics["gs"].params["b"]
    amu = constants.physical_constants["atomic mass constant"][0]
    omega = w * 1e-3 * constants.eV / constants.hbar
    b_si = b * constants.eV / (amu**1.5 * 1e-30)
    xi_si = b_si / (constants.hbar * omega) * (constants.hbar / (2 * omega)) ** 1.5
    assert SUB.anharmonic_strength(b, w) == pytest.approx(xi_si, rel=1e-12)


def test_cubic_coefficients_structure():
    c = SUB.cubic_correction_coefficients(5, 0.01)
    assert c.shape == (5, 8)
    np.testing.assert_array_equal(np.diag(c), np.ones(5))
    assert c[0, 1] == pytest.approx(-0.03) and c[2, 1] == pytest.approx(3 * 0.01 * 2**1.5)
    assert c[3, 0] == pytest.approx(0.01 * np.sqrt(6) / 3)
    np.testing.assert_array_equal(SUB.cubic_correction_coefficients(5, 0.0), np.eye(5, 8))


def test_cubic_with_b_zero_is_hermite(pair):
    h, c = SUB(*pair), SUB(*pair)
    h.compute_franck_condon_integrals_hermite(12, 12)
    c.compute_franck_condon_integrals_cubic(12, 12, 0.0, 0.0, q_i=pair[2], q_f=0.0)
    assert np.max(np.abs(c.fc_ints - h.fc_ints)) < 1e-13


def test_cubic_norms_and_pt_order(pair, cubics):
    w_i, w_f, dq = pair
    b_es, b_gs = cubics["es"].params["b"], cubics["gs"].params["b"]

    def pt(scale, sign=1.0, n=12):
        s = SUB(*pair)
        s.compute_franck_condon_integrals_cubic(n, n, sign * scale * b_es, sign * scale * b_gs, q_i=dq, q_f=0.0)
        return s

    def exact(scale, n=12):
        def well(p, center):
            return lambda q: p["a"] * (q - center) ** 2 + scale * p["b"] * (q - center) ** 3

        s = SUB(*pair)
        s.compute_franck_condon_integrals_fd(
            n, n, well(cubics["es"].params, dq), well(cubics["gs"].params, 0.0), q_i=dq, q_f=0.0
        )
        return s.fc_ints[:8, :8]

    full = pt(1.0, n=40)
    predicted = (SUB.cubic_correction_coefficients(40, full.cubic_xi_f) ** 2).sum(axis=1)
    assert np.max(np.abs(full.cubic_norms_f - predicted)) < 1e-12
    err_1 = np.max(np.abs(pt(1.0).fc_ints[:8, :8] - exact(1.0)))
    err_half = np.max(np.abs(pt(0.5).fc_ints[:8, :8] - exact(0.5)))
    assert 3.5 < err_1 / err_half < 4.5  # second-order error
    err_flip = np.max(np.abs(pt(1.0, sign=-1.0).fc_ints[:8, :8] - exact(1.0)))
    assert err_flip > 10 * err_1  # the frame sign of b matters


def test_cubic_warning(pair, capsys):
    SUB(*pair).compute_franck_condon_integrals_cubic(5, 5, 50.0, 50.0, q_i=pair[2], q_f=0.0)
    out = capsys.readouterr().out
    assert "|xi| =" in out and "_fd" in out


# ----------------------------------------------------------------------------- line positions from levels


def test_levels_lineshape_equals_bulid_fc_lsp_for_harmonic_ladders():
    s = SUB(64.31, 59.09, 0.6957)
    s.compute_franck_condon_integrals(ni=20, nf=25)
    kw = dict(eneaxis=AXIS, temp=300, sigma=[8, 25], zpl_lorentzian=True, gamma=0.3)
    ref = s.bulid_fc_lsp(**kw)[1].copy()
    got = s.build_fc_lsp_from_levels(energies_i=64.31 * np.arange(20), energies_f=59.09 * np.arange(25), **kw)[1]
    assert np.max(np.abs(got - ref)) <= 1e-14 * ref.max()


def test_levels_lineshape_uses_fd_levels(pair, diamond):
    s = SUB(*pair)
    s.compute_franck_condon_integrals_fd(10, 12, diamond["es"]["parabola"], diamond["gs"]["parabola"], q_i=pair[2])
    s.build_fc_lsp_from_levels(eneaxis=AXIS, temp=5, sigma=[8, 25], zpl_lorentzian=False)
    expected = (s.fd_energies_f - s.fd_energies_f[0])[None, :] * 1e3 - (s.fd_energies_i - s.fd_energies_i[0])[
        :, None
    ] * 1e3
    np.testing.assert_allclose(s.energy_v, expected)


def test_levels_lineshape_errors(pair):
    s = SUB(*pair)
    s.compute_franck_condon_integrals(ni=5, nf=5)
    with pytest.raises(ValueError):
        s.build_fc_lsp_from_levels(eneaxis=AXIS)
    with pytest.raises(ValueError):
        s.build_fc_lsp_from_levels(eneaxis=AXIS, energies_i=np.arange(3.0), energies_f=np.arange(5.0))
