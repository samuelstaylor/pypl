"""Physics checks of the original 1D CCD solver (Ruhoff recurrences + FC lineshape) on tutorial 002."""

import numpy as np
import pytest

from pypl.config_coord_1d_solver import config_coord_1d_solver

AXIS = np.linspace(-200.0, 1200.0, 1401)
LSP = dict(eneaxis=AXIS, sigma=[8.0, 25.0], zpl_lorentzian=True, gamma=0.3)
EZPL = 1945.0


def solver(omega_i, omega_f, delta_q, ni, nf, temp=5.0):
    s = config_coord_1d_solver(omega_i, omega_f, delta_q)
    s.compute_franck_condon_integrals(ni=ni, nf=nf)
    s.bulid_fc_lsp(temp=temp, **LSP)
    return s


@pytest.fixture(scope="module")
def pl(diamond):
    return solver(diamond["es"]["phonon"], diamond["gs"]["phonon"], diamond["delta_q"], 50, 60)


def test_closure(diamond):
    s = config_coord_1d_solver(diamond["es"]["phonon"], diamond["gs"]["phonon"], diamond["delta_q"])
    s.compute_franck_condon_integrals(ni=10, nf=300)
    rows = np.sum(s.fc_ints**2, axis=1)
    assert rows[0] == pytest.approx(1.0, abs=1e-8)
    np.testing.assert_allclose(rows, 1.0, atol=1e-6)
    assert 0.0 < s.fc_ints[0, 0] ** 2 < 1.0


def test_identity_at_zero_displacement(diamond):
    s = config_coord_1d_solver(diamond["gs"]["phonon"], diamond["gs"]["phonon"], 0.0)
    s.compute_franck_condon_integrals(ni=8, nf=8)
    np.testing.assert_allclose(s.fc_ints, np.eye(8), atol=1e-10)


def test_equal_frequency_symmetry(diamond):
    w = diamond["gs"]["phonon"]
    s = config_coord_1d_solver(w, w, diamond["delta_q"])
    s.compute_franck_condon_integrals(ni=30, nf=30)
    signs = np.fromfunction(lambda i, j: (-1.0) ** (i + j), s.fc_ints.shape)
    np.testing.assert_allclose(s.fc_ints, signs * s.fc_ints.T, atol=1e-9)


def test_lineshape_normalization_and_populations(pl, diamond):
    # The 0.3 meV Lorentzian ZPL sampled on a 1 meV grid integrates to ~1.01, not 1 (not a normalization error).
    assert np.sum(pl.fc_lineshape) * (AXIS[1] - AXIS[0]) == pytest.approx(1.0, abs=0.05)
    assert np.all(pl.fc_lineshape >= 0.0)
    assert np.sum(pl.prefactors) == pytest.approx(1.0, abs=1e-12)
    assert pl.prefactors[0] > 1.0 - 1e-6
    hot = solver(diamond["es"]["phonon"], diamond["gs"]["phonon"], diamond["delta_q"], 50, 60, temp=800.0)
    assert hot.prefactors[1] > 1e-3


def test_gaussian_zpl_integrates_to_one(diamond):
    s = config_coord_1d_solver(diamond["es"]["phonon"], diamond["gs"]["phonon"], diamond["delta_q"])
    s.compute_franck_condon_integrals(ni=50, nf=60)
    s.bulid_fc_lsp(eneaxis=AXIS, temp=5.0, sigma=[8.0, 25.0], zpl_lorentzian=False)
    assert np.sum(s.fc_lineshape) * (AXIS[1] - AXIS[0]) == pytest.approx(1.0, abs=1e-4)


def test_basis_convergence(pl, diamond):
    big = solver(diamond["es"]["phonon"], diamond["gs"]["phonon"], diamond["delta_q"], 80, 100)
    assert np.max(np.abs(big.fc_lineshape - pl.fc_lineshape)) / pl.fc_lineshape.max() < 1e-4


def test_spectra(pl, diamond):
    pl_axis, pl_spec = pl.compute_spectrum(tdm=1.0, zpl=EZPL, spectrum_type="PL")
    ab = solver(diamond["gs"]["phonon"], diamond["es"]["phonon"], diamond["delta_q"], 60, 50)
    ab_axis, ab_spec = ab.compute_spectrum(tdm=1.0, zpl=EZPL, spectrum_type="Abs")
    assert pl_axis[0] > pl_axis[-1] and ab_axis[0] < ab_axis[-1]
    assert np.all(pl_spec >= 0.0) and np.all(ab_spec >= 0.0)
    assert pl_axis[np.argmax(pl_spec)] == pytest.approx(EZPL, abs=1.0)
    assert ab_axis[np.argmax(ab_spec)] == pytest.approx(EZPL, abs=1.0)
    assert np.sum(pl_spec[pl_axis < EZPL]) > np.sum(pl_spec[pl_axis > EZPL])
    assert np.sum(ab_spec[ab_axis > EZPL]) > np.sum(ab_spec[ab_axis < EZPL])
    with pytest.raises(ValueError):
        pl.compute_spectrum(spectrum_type="bogus")
