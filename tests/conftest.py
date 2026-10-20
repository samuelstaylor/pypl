"""Shared fixtures: the NV- in diamond CCD of tutorial 002 and the NV- in MgO CCD of tutorial 003."""

import os
import re
import sys

import numpy as np
import pytest

TESTS_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.dirname(TESTS_DIR)
# Test the working tree even when an older pypl is installed.
for _path in (REPO_ROOT, TESTS_DIR):
    if _path not in sys.path:
        sys.path.insert(0, _path)

from pypl.config_coord_1d_potential import config_coord_1d_potential  # noqa: E402
from pypl.utils import parse_atoms_qexml, parse_total_energy_qexml, parse_west_summary  # noqa: E402

DIAMOND_DIR = os.path.join(REPO_ROOT, "doc", "tutorials", "002_nv_diamond_1d_ccd")
MGO_DIR = os.path.join(REPO_ROOT, "doc", "tutorials", "003_nv_mgo_isc")


def pytest_configure(config):
    config.addinivalue_line("markers", "slow: takes more than a few seconds (deselect with -m 'not slow')")


@pytest.fixture(scope="session")
def diamond():
    """Tutorial 002: images, delta_q and the near-minimum parabolas of NV- in diamond (GS at 0, ES at delta_q)."""

    masses = {"C": 12.0107, "N": 14.0067}
    symbols, gs_coord, _ = parse_atoms_qexml(os.path.join(DIAMOND_DIR, "gs_dft", "pwscf.xml"))
    _, es_coord, _ = parse_atoms_qexml(os.path.join(DIAMOND_DIR, "es_cdft", "pwscf.xml"))
    m = np.array([masses[s] for s in symbols])
    delta_q = np.linalg.norm((es_coord - gs_coord) * m[:, None] ** 0.5)
    rel = np.linspace(-0.2, 1.2, 15)
    out = {"delta_q": delta_q}
    energies = {}
    for state in ("gs", "es"):
        q, e = [], []
        for i in range(15):
            name = os.path.join(DIAMOND_DIR, "%s_1d_ccd" % state, "Image-%d" % (i + 1), "pwscf.xml")
            if os.path.exists(name):
                q.append(delta_q * rel[i])
                e.append(parse_total_energy_qexml(name))
        energies[state] = (np.array(q), np.array(e))
    zero = energies["gs"][1].min()
    for state, center, index in (("gs", 0.0, slice(None, 5)), ("es", delta_q, slice(4, None))):
        q, e = energies[state][0], energies[state][1] - zero
        pot = config_coord_1d_potential(q, e, q0=center)
        pot.fit_parabola(index)
        out[state] = {"q": q, "e": e, "center": center, "parabola": pot, "phonon": pot.phonon_energy()}
    return out


def mgo_image_q(name):
    match = re.fullmatch(r"Image_([+-]?\d+\.\d+)", name)
    return float(match.group(1)) if match else None


@pytest.fixture(scope="session")
def mgo():
    """Tutorial 003: 3E (sc root 1, at Q = 0) and 1A1 (sf root 4, at Q = delta_q) of NV- in MgO."""

    if not os.path.isdir(MGO_DIR):
        pytest.skip("doc/tutorials/003_nv_mgo_isc is not present")
    xml = lambda image: os.path.join(MGO_DIR, "qe", image, "data-file-schema.xml")  # noqa: E731
    masses = {"N": 14.0067, "O": 15.999, "Mg": 24.305}
    symbols, c0, _ = parse_atoms_qexml(xml("Image_+0.00"))
    _, c1, _ = parse_atoms_qexml(xml("Image_+1.00"))
    m = np.array([masses[s] for s in symbols])
    delta_q = np.linalg.norm((c1 - c0) * m[:, None] ** 0.5)
    zero = min(
        parse_total_energy_qexml(xml(n))
        for n in os.listdir(os.path.join(MGO_DIR, "qe"))
        if mgo_image_q(n) is not None and n != "Image_-1.00"
    )

    def load(run, root, exclude=()):
        points = []
        for name in os.listdir(os.path.join(MGO_DIR, run)):
            q = mgo_image_q(name)
            summary = os.path.join(MGO_DIR, run, name, "summary.json")
            if q is None or name in exclude or not (os.path.exists(summary) and os.path.exists(xml(name))):
                continue
            points.append((q, parse_total_energy_qexml(xml(name)) + parse_west_summary(summary)[root] - zero))
        return np.array(sorted(points)).T

    out = {"delta_q": delta_q}
    for state, run, root, exclude, center, window in (
        ("es", "sc", 0, (), 0.0, (-0.10, 0.45)),
        ("gs", "sf", 3, (), delta_q, (0.60, 1.15)),
    ):
        q_img, e = load(run, root, exclude)
        q = q_img * delta_q
        parabola = config_coord_1d_potential(q, e, q0=center)
        parabola.fit_parabola(np.where((q_img >= window[0] - 1e-9) & (q_img <= window[1] + 1e-9))[0])
        spline = config_coord_1d_potential(q, e, q0=center)
        spline.fit_spline(order=3, bc="not-a-knot", tails="smooth")
        out[state] = {
            "q_image": q_img,
            "q": q,
            "e": e,
            "center": center,
            "parabola": parabola,
            "spline": spline,
            "phonon": parabola.phonon_energy(),
        }
    return out
