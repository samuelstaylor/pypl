# validation/

The validation suite of PyPL's **1D configuration-coordinate (1D CCD)** path. It checks that independent routes to
the Franck–Condon matrix F_ij and the vibrational overlap function (VOF) agree where they must, and measures how
much anharmonic fits of the potential energy surfaces (PES) change the result where they should.

**Full documentation:**
- [`doc/validation.rst`](../doc/validation.rst) (*Validation and testing* in the docs): every test and check, its
  reference and tolerance, the datasets, the scripts, the metrics, reference results and known findings.
- [`doc/supporting_information.rst`](../doc/supporting_information.rst): every formula and derivation.
- [`METHODS.md`](METHODS.md): a tutorial-style guide to each method and fit, with measured accuracy and when to use it.
- [`doc/tutorials/pypl_005.ipynb`](../doc/tutorials/pypl_005.ipynb): a hands-on walk-through.

## What is tested

| Question | Where | Pass criterion |
| --- | --- | --- |
| Does the original PyPL API still give the tutorial results? | `tests/test_legacy_tutorials.py` | 35 arrays equal to a reference from the original code |
| Are the solver's outputs physical (closure, symmetry, populations, spectra)? | `check_1d_ccd.py`, `tests/` | exact identities to 1e-6 … 1e-12 |
| Does the extension change nothing it inherits? | `check_1d_ccd.py`, `tests/` | bit-identical arrays |
| Do Ruhoff, Hermite quadrature and finite differences agree on harmonic wells? | `compare_vof.py`, `check_1d_ccd.py` | max \|ΔF_ij\| ≤ 1e-7, max \|ΔVOF\|/peak ≤ 1e-8 |
| Is the finite-difference solver converged and of the expected order? | `vof_report.py` section 02 | log–log slopes 2.00 (3-point), 4.00 (5-point) |
| Are the fits what their names say (interpolation, end conditions, tails, boundedness)? | `check_1d_ccd.py`, `tests/` | residuals ≤ 1e-12 eV, C1/C2 junctions, unbounded fits flagged |
| Is first-order cubic perturbation theory correct, and where is it valid? | `check_1d_ccd.py`, `vof_report.py` section 07, tutorials 004/004a | error ∝ b², effect ∝ b, mixing criterion |
| How much does anharmonicity change the VOF? | `compare_vof.py --fd-fit …`, `vof_report.py` section 06, `run_all.py` | reported as the effect (DIFFERS), not a failure |

## Datasets and naming of the data points

Two datasets (`--dataset`, or `PYPL_VALIDATION_DATASET`), defined in `harness/ccd_data.py`:

| Dataset | Transition | Lower state ("GS") | Upper state ("ES") |
| --- | --- | --- | --- |
| `nv_diamond` (default) | NV⁻ in diamond, 3E → 3A2 | 3A2: DFT (PBE) | 3E: ΔSCF, constrained occupations (PBE) |
| `nv_mgo` | NV⁻ in MgO, 3E → 1A1 | 1A1: spin-flip TDDFT@DDH (WEST, root 4) | 3E: spin-conserving TDDFT@DDH (WEST, root 1) |

The data points are **computed energies** at the images, not "DFT energies": the method differs between datasets and
between the two states (ground-state DFT, ΔSCF, spin-conserving or spin-flip TDDFT, and in other systems BSE or
embedding methods). Each `Dataset` records `energy_methods` (per state) and `energy_zero`; figures label the points
with them (`ccd_data.energy_label`) and every log and `index.md` prints them.

## Layout

```
validation/
├── harness/                      helpers shared by the scripts (not installed)
│   ├── ccd_data.py                 the datasets: images, energies and methods, parabola wells, roles and labels
│   ├── fits.py                     fit specifications (parabola, polyN, morse, cubic, spline...) -> fitted potentials
│   ├── methods.py                  VOFParams, shared flags, the F_ij methods, compute_vof, compare_results
│   ├── output.py                   OutputDir: one documented directory per run; archiving of reruns
│   └── plots.py                    shared style and comparison figures
├── check_1d_ccd.py               81 physics and regression checks (exit 0 on success, -v lists them)
├── compare_vof.py                PASS/FAIL: a reference method against test methods
├── plot_ccd_wavefunctions.py     CCD + wavefunctions + VOF panel
├── vof_report.py                 the seven-section report
├── run_1d_ccd.py                 tutorial 002 as a script
├── run_all.py                    the master suite, configured by run_all.toml / run_all_nv_mgo.toml
├── data/nv_mgo_west/             how the MgO data were collected on the HPC (the data: doc/tutorials/003_nv_mgo_isc/)
├── METHODS.md                    the method guide
└── outputs/                      results (reruns are archived to outputs/_archive/, never deleted)
```

## Running

No installation is needed for the scripts (`harness` puts the repository root on `sys.path`); the tests need
`pip install ".[test]"`.

```bash
pytest                                              # package tests, ~1 min
python validation/check_1d_ccd.py                   # 81 checks, ~15 s (-v: every check and its value)
python validation/compare_vof.py                    # pypl vs hermite and fd, ~30 s
python validation/compare_vof.py --test fd --fd-fit spline3-notaknot-smooth   # anharmonic effect of a fit
python validation/compare_vof.py --reference fd --fd-fit cubic --test cubic   # perturbation theory vs exact cubic
python validation/plot_ccd_wavefunctions.py --wf-method fd --fd-fit spline3-notaknot-smooth --vof-method pypl fd
python validation/vof_report.py                     # everything, ~8 min (--skip-convergence / --skip-cubic)
python validation/run_all.py                        # every method, fit and variant, NV- in diamond, ~20 min
python validation/run_all.py --config validation/run_all_nv_mgo.toml          # the same for NV- in MgO
python validation/run_all.py --set physics.temp=300 --name 300K               # override any key
python validation/run_all.py --dry-run              # list the commands without running them
```

Every setting is a flag or a `run_all.toml` key (table in `doc/validation.rst`, *Configuration*). Shared flags:
`--dataset --temp --sigma --gamma --no-zpl --order-es --order-gs --ene-range --resol --line-positions`, the
quadrature grid (`--npts` or `--hermite-h`, `--hermite-margin`), the finite-difference grid (`--fd-npts` or `--fd-h`,
`--fd-range QMIN QMAX` or `--fd-margin`, `--fd-stencil {3,5,7}`, `--fd-shift`), `--fd-fit`, figures (`--dpi`,
`--fig-format png pdf svg`, `--font-size`), pass criteria (`--tol-fc`, `--tol-vof`; `check_1d_ccd.py --tol NAME=VALUE`,
`--list-tolerances`; `plot_ccd_wavefunctions.py --check-tol`), the spline options (`--spline-order --spline-bc --spline-smoothing
--spline-tails`) and the cubic options (`--cubic-points --cubic-free-a --cubic-normalize`); `--tag` names the run and
`-v` prints every setting. Fit specs: `parabola`, `polyN`, `morse`, `cubic[-near][-free]`,
`spline[K][-natural|-notaknot|-clamped|-sRMS][-smooth]`; the recommended anharmonic fit is
`spline3-notaknot-smooth`.

`run_all.toml` documents every key (`[suite]`, `[physics]`, `[methods]`, `[hermite]`, `[fd]` with `[[fd.splines]]`,
`[[cubic.variants]]`, `[plots]`, `[report]`, `[tolerances]`). Misspelled keys and invalid fits are rejected before
anything runs; unbounded fits are skipped with the reason. Step statuses: PASS, OK, DIFFERS (expected), FAIL, ERROR.

## Outputs

Every run writes `outputs/<script>/<run>/` (a suite: `outputs/suite/<name>/<script>/<run>/`; another dataset outside a
suite: `outputs/<dataset>/<script>/<run>/`):

| File | Content |
| --- | --- |
| `index.md` | **Start here**: headline results and every file with a description |
| `log.txt` | full console output, including every method's settings |
| `settings.json` | every parameter |
| `data/`, `figures/`, `NN_section/` | CSV tables and the figures made from them |

A rerun under the same name moves the previous run to `outputs/_archive/<same path>/<timestamp>/`; nothing is deleted.
`PYPL_VALIDATION_OUTPUT_ROOT` redirects all outputs.

## Maintenance

- Whenever a method or fit is added or changed, update `METHODS.md`, `doc/validation.rst` and
  `doc/supporting_information.rst`, and add tests and checks.
- `check_1d_ccd.py` pins reference values in `REFERENCE`; change them only for an intended change of physics, never
  widen a tolerance to make a check pass.
- Black at line length 120, as CI checks it: `python -m black --check --line-length 120 pypl tests validation`.
