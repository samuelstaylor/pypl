"""Intersystem crossing rate from Fermi's golden rule."""

import numpy as np
import pytest
from scipy import constants

from pypl.isc import SOC_UNITS, compute_isc_lifetime, compute_isc_rate


def test_rate_formula():
    soc_ghz, vof = 204.9925, 2.70e-7
    lam = constants.h * soc_ghz * 1e9
    expected = 4 * np.pi / constants.hbar * lam**2 * vof / (constants.eV * 1e-3)
    assert compute_isc_rate(soc_ghz, vof) == pytest.approx(expected, rel=1e-14)
    assert compute_isc_rate(soc_ghz, vof) / 1e6 == pytest.approx(3.705, abs=1e-3)  # MHz


def test_units_are_consistent():
    lam_mev = constants.h * 100e9 / constants.eV * 1e3
    assert compute_isc_rate(lam_mev, 1e-4, soc_unit="meV") == pytest.approx(compute_isc_rate(100.0, 1e-4), rel=1e-12)
    assert compute_isc_rate(lam_mev * 1e3, 1e-4, soc_unit="ueV") == pytest.approx(
        compute_isc_rate(100.0, 1e-4), rel=1e-12
    )
    assert set(SOC_UNITS) == {"GHz", "meV", "ueV", "cm-1"}
    with pytest.raises(ValueError):
        compute_isc_rate(1.0, 1.0, soc_unit="Hz")


def test_array_and_lifetime():
    vof = np.array([0.0, 1e-5, 2e-5])
    rate = compute_isc_rate(200.0, vof)
    assert rate.shape == (3,) and rate[2] == pytest.approx(2 * rate[1])
    tau = compute_isc_lifetime(rate)
    assert np.isinf(tau[0]) and tau[1] == pytest.approx(1 / rate[1])
    assert compute_isc_lifetime(4.0) == 0.25
