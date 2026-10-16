"""The computations of tutorials 001 and 002, as plain functions returning arrays.

The same functions produce the pinned reference (``tests/data/legacy_reference.npz``, generated from the original
PyPL code by ``tests/tools/make_legacy_reference.py``) and the values the test suite compares against it, so any
change in the behaviour of the original PyPL API shows up as a test failure.

Only the original API is used here: ``hr_solver``, ``config_coord_1d_solver`` and the parsers in ``utils``.
"""

import os

import numpy as np
from scipy import constants
from scipy.optimize import curve_fit

TUTORIALS = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "doc", "tutorials")
DATA_001 = os.path.join(TUTORIALS, "001_nv_diamond_abs_pl")
DATA_002 = os.path.join(TUTORIALS, "002_nv_diamond_1d_ccd")
MASS_LIST = {"C": 12.0107, "N": 14.0067}
EZPL = 1945  # meV


def tutorial_001():
    """Tutorial 001: NV- in diamond from full phonons, displacement and force routes, PL and absorption."""
    from pypl.hr_solver import hr_solver
    from pypl.utils import parse_atoms_qexml, parse_forces_qexml, parse_phonopy_h5

    out = {}
    gs_freqs, gs_modes = parse_phonopy_h5(os.path.join(DATA_001, "phonon", "gs_ph_mesh.hdf5"))
    es_freqs, es_modes = parse_phonopy_h5(os.path.join(DATA_001, "phonon", "es_ph_mesh.hdf5"))
    symbols, gs_coord, cell = parse_atoms_qexml(os.path.join(DATA_001, "gs_dft", "pwscf.xml"))
    _, es_coord, _ = parse_atoms_qexml(os.path.join(DATA_001, "es_cdft", "pwscf.xml"))
    _, forces = parse_forces_qexml(os.path.join(DATA_001, "gs_forces_at_es_coord", "pwscf.xml"))
    out["t001_gs_freqs"], out["t001_gs_coord"], out["t001_cell"], out["t001_forces"] = gs_freqs, gs_coord, cell, forces

    # PL from displacements.
    pl_dis = hr_solver()
    hrf = pl_dis.compute_hrf_dis(gs_freqs, gs_modes, symbols, gs_coord, es_coord, cell, mass_list=MASS_LIST)
    out["t001_pl_dis_freqs"], out["t001_pl_dis_hrf"] = hrf["freqs"], hrf["hr_factors"]
    for i, sigma in enumerate([[3, 2], [6, 2], [6, 6]]):
        axis, density = pl_dis.compute_spectral_density(hrf, energy_axis=np.linspace(0, 200, 201), sigma=sigma)
        out["t001_pl_dis_spectral_density_%d" % i] = density
    axis, lineshape = pl_dis.compute_lineshape_fft(hrf, temp=4, sigma=[6, 2], zpl_broadening=0.3)
    out["t001_pl_dis_lineshape_axis"], out["t001_pl_dis_lineshape"] = axis, lineshape
    energy, spectrum = pl_dis.compute_spectrum(EZPL, spectrum_type="PL", lineshape=(axis, lineshape))
    out["t001_pl_dis_spectrum_axis"], out["t001_pl_dis_spectrum"] = energy, spectrum

    # Time-domain integration (a coarser time grid than the default, to keep the test fast).
    axis, lineshape = pl_dis.compute_lineshape_numerical_integration(
        hrf, temp=4, sigma=[6, 2], zpl_broadening=0.3, time_resolution=20001, lineshape_energy_resolution=141
    )
    out["t001_pl_dis_lineshape_numint_axis"], out["t001_pl_dis_lineshape_numint"] = axis, lineshape

    # PL from forces, at several temperatures.
    pl_forces = hr_solver()
    hrf_f = pl_forces.compute_hrf_forces(gs_freqs, gs_modes, symbols, forces, MASS_LIST)
    out["t001_pl_forces_hrf"] = hrf_f["hr_factors"]
    for temp, gamma in zip([4, 150, 300], [0.3, 1.8, 4.8]):
        lineshape = pl_forces.compute_lineshape_fft(hrf_f, temp=temp, sigma=[6.0, 2.0], zpl_broadening=gamma)
        energy, spectrum = pl_forces.compute_spectrum(EZPL, spectrum_type="PL", lineshape=lineshape)
        out["t001_pl_forces_spectrum_%dK" % temp] = spectrum

    # Absorption from displacements with the excited-state phonons.
    abs_dis = hr_solver()
    hrf_a = abs_dis.compute_hrf_dis(es_freqs, es_modes, symbols, gs_coord, es_coord, cell, mass_list=MASS_LIST)
    out["t001_abs_dis_hrf"] = hrf_a["hr_factors"]
    lineshape = abs_dis.compute_lineshape_fft(hrf_a, temp=4, sigma=[6.0, 2.0], zpl_broadening=0.3)
    energy, spectrum = abs_dis.compute_spectrum(EZPL, spectrum_type="Abs", lineshape=lineshape)
    out["t001_abs_dis_spectrum_axis"], out["t001_abs_dis_spectrum"] = energy, spectrum
    return out


def tutorial_002():
    """Tutorial 002: NV- in diamond in the 1D configuration-coordinate model, PL and absorption."""
    from pypl.config_coord_1d_solver import config_coord_1d_solver
    from pypl.utils import parse_atoms_qexml, parse_total_energy_qexml

    out = {}
    symbols, gs_coord, _ = parse_atoms_qexml(os.path.join(DATA_002, "gs_dft", "pwscf.xml"))
    _, es_coord, _ = parse_atoms_qexml(os.path.join(DATA_002, "es_cdft", "pwscf.xml"))
    masses = np.array([MASS_LIST[s] for s in symbols])
    delta_q = np.linalg.norm((es_coord - gs_coord) * masses[:, None] ** 0.5)

    rel = np.linspace(-0.2, 1.2, 15)
    gs_idx, es_idx = [], []
    gs_e, es_e = np.zeros(rel.shape[0]), np.zeros(rel.shape[0])
    for i in range(rel.shape[0]):
        name = os.path.join(DATA_002, "gs_1d_ccd", "Image-%d" % (i + 1), "pwscf.xml")
        if os.path.exists(name):
            gs_idx.append(i)
            gs_e[i] = parse_total_energy_qexml(name)
        name = os.path.join(DATA_002, "es_1d_ccd", "Image-%d" % (i + 1), "pwscf.xml")
        if os.path.exists(name):
            es_idx.append(i)
            es_e[i] = parse_total_energy_qexml(name)
    ref = np.min(gs_e)
    gs_e, es_e = gs_e[gs_idx] - ref, es_e[es_idx] - ref
    gs_q, es_q = delta_q * rel[gs_idx], delta_q * rel[es_idx]

    gs_params = curve_fit(lambda x, a, c: a * x**2 + c, gs_q[:5], gs_e[:5])[0]
    es_params = curve_fit(lambda x, a, c: a * (x - delta_q) ** 2 + c, es_q[4:], es_e[4:])[0]

    def phonon(a):
        w2 = 2 * a * constants.eV / (1e-10**2 * constants.physical_constants["atomic mass constant"][0])
        return np.sqrt(w2 * constants.hbar**2 / constants.eV**2) * 1000

    gs_phonon, es_phonon = phonon(gs_params[0]), phonon(es_params[0])
    out["t002_delta_q"] = np.array(delta_q)
    out["t002_gs_energies"], out["t002_es_energies"] = gs_e, es_e
    out["t002_phonons"] = np.array([gs_phonon, es_phonon])

    axis = np.linspace(-200, 1200, 1401)
    pl = config_coord_1d_solver(es_phonon, gs_phonon, delta_q)
    pl.compute_franck_condon_integrals(ni=50, nf=60)
    pl.bulid_fc_lsp(eneaxis=axis, temp=5, sigma=[8, 25], zpl_lorentzian=True, gamma=0.3)
    out["t002_pl_fc_ints"], out["t002_pl_lineshape"] = pl.fc_ints, pl.fc_lineshape
    out["t002_pl_spectrum_axis"], out["t002_pl_spectrum"] = pl.compute_spectrum(tdm=1.0, zpl=EZPL, spectrum_type="PL")

    ab = config_coord_1d_solver(gs_phonon, es_phonon, delta_q)
    ab.compute_franck_condon_integrals(ni=60, nf=50)
    ab.bulid_fc_lsp(eneaxis=axis, temp=5, sigma=[8, 25], zpl_lorentzian=True, gamma=0.3)
    out["t002_abs_fc_ints"], out["t002_abs_lineshape"] = ab.fc_ints, ab.fc_lineshape
    out["t002_abs_spectrum_axis"], out["t002_abs_spectrum"] = ab.compute_spectrum(
        tdm=1.0, zpl=EZPL, spectrum_type="Abs"
    )

    # Thermal populations and a Gaussian ZPL exercise the other branches of bulid_fc_lsp.
    pl.bulid_fc_lsp(eneaxis=axis, temp=300, sigma=[8, 25], zpl_lorentzian=False)
    out["t002_pl_lineshape_300K_gaussian_zpl"] = pl.fc_lineshape
    return out


def all_cases():
    out = tutorial_001()
    out.update(tutorial_002())
    return {k: np.asarray(v) for k, v in out.items()}
