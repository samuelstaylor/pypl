"""The original PyPL API must reproduce tutorials 001 and 002 exactly as the original code did.

``tests/data/legacy_reference.npz`` was generated from the untouched original package (branch ``pypl-original``) by
``tests/tools/make_legacy_reference.py``. Every array the tutorials compute is compared against it. Do not regenerate
the reference to make this test pass: a failure means the behaviour of the original API changed.
"""

import os

import numpy as np
import pytest

import legacy_cases

REFERENCE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data", "legacy_reference.npz")

# Across platforms / numpy builds the FFT and BLAS may round differently; on one machine the match is bit for bit.
RTOL = 1e-10


@pytest.fixture(scope="module")
def reference():
    with np.load(REFERENCE) as data:
        return {k: data[k] for k in data.files}


@pytest.fixture(scope="module")
def current():
    return legacy_cases.all_cases()


def test_reference_comes_from_the_original_code(reference):
    assert str(reference["_source_commit"]) not in ("", "unknown")


def test_same_quantities(reference, current):
    assert set(current) == set(reference) - {"_source_commit"}


@pytest.mark.parametrize("key", sorted(k for k in np.load(REFERENCE).files if k != "_source_commit"))
def test_tutorial_array_unchanged(reference, current, key):
    got, expected = current[key], reference[key]
    assert got.shape == expected.shape
    scale = np.max(np.abs(expected)) if expected.size else 0.0
    np.testing.assert_allclose(got, expected, rtol=RTOL, atol=RTOL * scale, err_msg=key)


# Numbers printed in the stored outputs of doc/tutorials/pypl_001.ipynb and pypl_002.ipynb. The notebooks were run
# with other numpy / scipy versions, so they agree to their printed precision, not bit for bit.
NOTEBOOK_001 = {
    "total HRF, PL, displacements": ("t001_pl_dis_hrf", 2.985289200067),
    "total HRF, PL, forces": ("t001_pl_forces_hrf", 2.942036502131),
    "total HRF, Abs, displacements": ("t001_abs_dis_hrf", 3.102826963299),
}
NOTEBOOK_001_TOP_MODES = [  # (1-based index, freq in meV, HRF): the ten largest partial HR factors
    (64, 62.3239, 8.3595e-01),
    (6, 34.4922, 3.5594e-01),
    (36, 49.0021, 3.2950e-01),
    (252, 81.4371, 1.0780e-01),
    (30, 48.2382, 8.3555e-02),
    (144, 70.7789, 7.6982e-02),
    (5, 34.3174, 6.6586e-02),
    (17, 42.5984, 6.2820e-02),
    (122, 68.2145, 5.7055e-02),
    (321, 89.3907, 4.9803e-02),
]


@pytest.mark.parametrize("label", sorted(NOTEBOOK_001))
def test_notebook_001_total_huang_rhys(current, label):
    key, printed = NOTEBOOK_001[label]
    assert np.sum(current[key]) == pytest.approx(printed, rel=1e-8)


def test_notebook_001_largest_modes(current):
    from scipy import constants

    hrf, freqs = current["t001_pl_dis_hrf"], current["t001_pl_dis_freqs"]
    top = np.argsort(hrf)[-10:][::-1]
    for i, (index, freq, s) in zip(top, NOTEBOOK_001_TOP_MODES):
        assert i + 1 == index
        assert freqs[i] * constants.hbar / constants.eV * 1e3 == pytest.approx(freq, abs=1e-4)
        assert hrf[i] == pytest.approx(s, rel=1e-4)


def test_notebook_002_printed_values(current):
    from pypl.config_coord_1d_vof_solver import config_coord_1d_vof_solver

    delta_q = float(current["t002_delta_q"])
    gs_phonon, es_phonon = current["t002_phonons"]
    assert delta_q == pytest.approx(6.9573191451e-01, rel=1e-10)
    assert gs_phonon == pytest.approx(59.08890, abs=1e-5)
    assert es_phonon == pytest.approx(64.31191, abs=1e-5)
    assert config_coord_1d_vof_solver.huang_rhys(gs_phonon, delta_q) == pytest.approx(3.42111, abs=1e-5)
    assert config_coord_1d_vof_solver.huang_rhys(es_phonon, delta_q) == pytest.approx(3.72351, abs=1e-5)
