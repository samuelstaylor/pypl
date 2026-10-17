## **PyPL** code repository

PyPL is a Python package for computing vibrationally resolved optical spectra of point defects, using atomic geometries and phonons obtained from density functional theory (DFT) calculations.

## Features

- **Full-phonon lineshapes** (`hr_solver`): partial Huang-Rhys factors from displacements or forces, spectral densities, and temperature-dependent photoluminescence / absorption lineshapes by the generating-function method.
- **1D configuration-coordinate model** (`config_coord_1d_solver`): Franck-Condon integrals of two displaced harmonic oscillators by the Ruhoff recurrences and the broadened Franck-Condon lineshape.
- **Anharmonic 1D CCD** (`config_coord_1d_vof_solver`, a subclass of `config_coord_1d_solver` that only adds methods):
  - Franck-Condon integrals by quadrature of harmonic eigenfunctions (stable for hundreds of levels, i.e. large Huang-Rhys factors);
  - finite-difference eigenstates of any potential, e.g. a fit of the computed energies at the images;
  - a first-order cubic anharmonic correction;
  - a lineshape with lines at the actual vibrational level energies (`build_fc_lsp_from_levels`).
- **Potential fits** (`config_coord_1d_potential`): parabola, polynomial, Morse, cubic, and splines of degree 2-5 (natural, not-a-knot or clamped ends; harmonic or smooth C² tails).
- **Intersystem crossing** (`pypl.isc`): golden-rule ISC rate Γ(Δ) = 4π/ħ |λ⊥|² F(Δ) and lifetime from any vibrational overlap function.
- **Parsers** (`pypl.utils`): Quantum ESPRESSO XML (structure, forces, total energy), Phonopy (HDF5, YAML) and WEST TDDFT `summary.json`.

## Documentation

The tutorial and documentation are hosted on [GitHub Pages](https://miccompy.github.io/pypl/). The tutorials and their input data (`001_*`, `002_*`, `003_*`) are in `doc/tutorials/`:

| Notebook | Content |
| --- | --- |
| `pypl_001.ipynb` | PL and absorption of NV⁻ in diamond from full phonons (displacement and force routes) |
| `pypl_002.ipynb` | PL and absorption of NV⁻ in diamond in the 1D CCD model |
| `pypl_003.ipynb` | 3E → 1A1 intersystem crossing of NV⁻ in MgO: harmonic vs anharmonic (smooth spline) vibrational overlap and ISC rate, from WEST TDDFT data |
| `pypl_004.ipynb` | Cubic anharmonicity by first-order perturbation theory: corrected wavefunctions, Franck-Condon integrals and VOF, the regime of validity, and a supplementary derivation |
| `pypl_004a.ipynb` | The same analysis for the strongly anharmonic NV⁻ in MgO (3E → 1A1): where and by how much perturbation theory fails |
| `pypl_005.ipynb` | Running the test suite and the validation checks, and comparing potential fits and VOF methods |

## Installation

Run `pip install .` in the main folder. For development, `pip install -e ".[dev]"` adds pytest, black and nbsphinx.

## Testing

```bash
pytest                          # the test suite in tests/ (~1 min)
pytest -m "not slow"            # skip the slowest end-to-end test
python validation/check_1d_ccd.py   # the 1D CCD regression checks of the validation suite
```

`tests/test_legacy_tutorials.py` guarantees that the original API is unchanged: it reruns the computations of tutorials 001 and 002 and compares every array with a reference generated from the original code (branch `pypl-original`), and checks the numbers printed in the notebooks. `validation/` holds the method-comparison suite for the 1D CCD (`validation/README.md`, `validation/METHODS.md`). Everything that is tested, with references and tolerances, is documented in `doc/validation.rst`, and every formula and derivation in `doc/supporting_information.rst` (the *Validation and testing* and *Supporting information* pages of the documentation).

Formatting is checked with `black --check --line-length 120 pypl tests validation`.

## Authors

Yu Jin (jinyuchem@uchicago.edu)

Samuel S. Taylor (sstaylor@uchicago.edu)

Victor Yu (yuw@anl.gov)

Marco Govoni (mgovoni@unimore.it)

Giulia Galli (gagalli@uchicago.edu)

## Contact

Please use the GitHub [issue tracker](https://github.com/MICCoMPy/pypl/issues/) for bug reports. Contributions to new features are welcome.
