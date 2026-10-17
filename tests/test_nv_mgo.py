"""NV- in MgO (tutorial 003): the WEST parser, the CCD and the 3E -> 1A1 VOFs at large Huang-Rhys factors."""

import os

import numpy as np
import pytest

from conftest import MGO_DIR
from pypl.config_coord_1d_solver import config_coord_1d_solver
from pypl.config_coord_1d_vof_solver import config_coord_1d_vof_solver
from pypl.utils import parse_west_summary, ry2ev

SUB = config_coord_1d_vof_solver


def test_parse_west_summary(tmp_path):
    f = tmp_path / "summary.json"
    f.write_text('{"plep": {"n_plep_eigen": 2, "eigenval": [0.1, 0.25]}}')
    np.testing.assert_allclose(parse_west_summary(str(f)), [0.1 * ry2ev, 0.25 * ry2ev])
    bad = tmp_path / "bad.json"
    bad.write_text('{"other": 1}')
    with pytest.raises(ValueError):
        parse_west_summary(str(bad))


def test_mgo_ccd(mgo):
    assert mgo["delta_q"] == pytest.approx(3.097879, abs=1e-6)
    assert (len(mgo["es"]["q"]), len(mgo["gs"]["q"])) == (33, 30)
    es, gs = mgo["es"], mgo["gs"]
    assert es["q_image"][np.argmin(es["e"])] == 0.0 and gs["q_image"][np.argmin(gs["e"])] == 1.0
    assert es["e"].min() == pytest.approx(3.1376, abs=5e-5)
    assert 1e3 * (es["e"].min() - gs["e"].min()) == pytest.approx(554.4, abs=0.05)
    # The harmonic parabolas of tutorial 003.
    assert es["phonon"] == pytest.approx(60.8332, abs=1e-4)
    assert gs["phonon"] == pytest.approx(31.5511, abs=1e-4)
    assert len(es["parabola"].index) == len(gs["parabola"].index) == 12
    assert es["spline"].bounded and gs["spline"].bounded


def test_ruhoff_row_zero_exact_at_large_huang_rhys(mgo):
    """Large S: row 0 of the recurrence is exact, higher rows lose accuracy (hermite is the reference there)."""

    w_i, w_f, dq = mgo["es"]["phonon"], mgo["gs"]["phonon"], mgo["delta_q"]
    ruhoff = config_coord_1d_solver(w_i, w_f, dq)
    ruhoff.compute_franck_condon_integrals(ni=4, nf=100)
    herm = SUB(w_i, w_f, dq)
    herm.compute_franck_condon_integrals_hermite(4, 100)
    assert np.max(np.abs(ruhoff.fc_ints[0] - herm.fc_ints[0])) < 1e-12
    assert np.sum(herm.fc_ints[0] ** 2) == pytest.approx(1.0, abs=1e-6)


def test_fd_matches_hermite_in_mgo_orientation(mgo):
    w_i, w_f, dq = mgo["es"]["phonon"], mgo["gs"]["phonon"], mgo["delta_q"]
    ref = SUB(w_i, w_f, dq)
    ref.compute_franck_condon_integrals_hermite(12, 12)
    fd = SUB(w_i, w_f, dq)
    fd.compute_franck_condon_integrals_fd(12, 12, mgo["es"]["parabola"], mgo["gs"]["parabola"], q_i=0.0, q_f=dq)
    assert np.max(np.abs(fd.fc_ints - ref.fc_ints)) < 1e-7


@pytest.mark.slow
def test_tutorial_003_isc_numbers(mgo):
    """The headline numbers of tutorial 003 (harmonic vs smooth not-a-knot spline VOF and ISC rate)."""

    from pypl.isc import compute_isc_rate

    axis = np.linspace(0, 2000, 2001)
    kw = dict(eneaxis=axis, temp=5, sigma=[6, 8], zpl_lorentzian=True, gamma=0.3)
    w_i, w_f, dq = mgo["es"]["phonon"], mgo["gs"]["phonon"], mgo["delta_q"]
    harm = config_coord_1d_solver(w_i, w_f, dq)
    harm.compute_franck_condon_integrals(ni=100, nf=100)
    harm.bulid_fc_lsp(**kw)
    spl = SUB(w_i, w_f, dq)
    spl.compute_franck_condon_integrals_fd(
        100, 100, mgo["es"]["spline"], mgo["gs"]["spline"], q_i=0.0, q_f=dq, npts=40001
    )
    spl.build_fc_lsp_from_levels(**kw)

    delta = 554.41449
    f_harm, f_spl = np.interp(delta, axis, harm.fc_lineshape), np.interp(delta, axis, spl.fc_lineshape)
    assert f_harm == pytest.approx(2.70e-7, rel=5e-3)
    assert f_spl == pytest.approx(3.97e-5, rel=5e-3)
    assert compute_isc_rate(204.9925, f_harm) / 1e6 == pytest.approx(3.711, rel=5e-3)
    assert compute_isc_rate(204.9925, f_spl) / 1e6 == pytest.approx(545.1, rel=5e-3)
    assert axis[np.argmax(harm.fc_lineshape)] == pytest.approx(1145, abs=2)
    assert axis[np.argmax(spl.fc_lineshape)] == pytest.approx(660, abs=2)


def test_data_present():
    assert os.path.isdir(MGO_DIR), "doc/tutorials/003_nv_mgo_isc is part of the repository"
