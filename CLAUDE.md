# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

PyPL computes vibrationally resolved optical spectra (photoluminescence and absorption) of point defects in solids, starting from geometries, forces, and phonons produced by DFT codes (Quantum ESPRESSO + Phonopy). It is a small pure-Python library — no compiled extensions, no runtime state beyond the solver objects.

## Commands

```bash
pip install -e ".[dev]"                          # editable install + pytest, black, nbsphinx (installed locally 2026-10-16)
pytest                                           # tests/, ~70 s; -m "not slow" skips the tutorial-003 end-to-end test
python validation/check_1d_ccd.py                # 81 1D CCD checks
black --check --diff --color --line-length 120 pypl tests validation   # lint exactly as CI does
cd doc && make html                              # build Sphinx docs into doc/_build (needs nbsphinx + pandoc)
```

CI (`.github/workflows/ci.yml`): job `build` = black check on `pypl tests validation` + `pip install .`; job `test` =
`pip install ".[test]"`, `pytest -W error::SyntaxWarning` and `validation/check_1d_ccd.py` on Python 3.9 / 3.11 / 3.12.
Black with `--line-length 120` is mandatory (`[tool.black]` in pyproject.toml also sets it and excludes
`validation/data/`). `black` 26.10.1 is installed locally (CI's conda black is newer than the 26.1 pin in `[dev]`).

**Branches (2026-10-16):** `main` and `pypl-original` = the original package (upstream 50bd02a), never modified.
`vof-extensions` = this work: commit 7ed4613 is the raw snapshot (testing/ + pypl_ext), later commits integrate it into
the package. A worktree of `pypl-original` was made at `../pypl_original_wt` to generate the legacy reference.

**`tests/` (pytest)** — the package's test suite:
- `test_legacy_tutorials.py`: `legacy_cases.py` reruns tutorials 001 / 002 (35 arrays: HR factors, spectral densities,
  FFT and time-integration lineshapes, spectra at 4/150/300 K, 1D CCD F_ij, lineshapes, spectra, a 300 K Gaussian-ZPL
  lineshape) and compares them with `tests/data/legacy_reference.npz`, generated FROM THE ORIGINAL CODE by
  `tests/tools/make_legacy_reference.py --pypl-root ../pypl_original_wt` (records the source commit). Locally all 35 are
  bit-identical; the test allows rtol 1e-10 for other platforms. Also checks the numbers printed in the notebooks
  (they were run with other numpy/scipy versions: total HRF 2.985289200067 vs 2.985289195896 now, both codes).
  **Never regenerate the reference from modified code.**
- `test_config_coord_1d_solver.py`, `test_config_coord_1d_vof_solver.py` (contract, harmonic functions, hermite / fd /
  cubic routes, orientation, levels lineshape), `test_config_coord_1d_potential.py`, `test_isc.py`, `test_nv_mgo.py`
  (WEST parser, MgO CCD, Ruhoff row-0 at large S, tutorial-003 numbers, marked slow). Fixtures `diamond` / `mgo` in
  `conftest.py` build both CCDs with the package API. 103 tests, all pass.

For the multi-mode path the legacy test pins every tutorial-001 array. `pypl/config_coord_1d_solver.py` also has an `if __name__ == "__main__"` block that reproduces an NV-center PL/Abs plot.

## Architecture

There is no `pypl/__init__.py`; `pypl` is an implicit namespace package. Import submodules directly: `from pypl.hr_solver import hr_solver`, `from pypl.utils import *`.

Two independent physical methods share only the parsers in `utils.py`:

**Multi-mode (full phonon) path** — the main workflow:
- `utils.py` — parsers producing "user units": `parse_atoms_qexml` / `parse_forces_qexml` / `parse_total_energy_qexml` (QE `pwscf.xml` → Å, eV/Å, eV), `parse_phonopy_h5` / `parse_phonopy_yaml` (→ THz, real eigenmodes).
- `hr_solver.py` (`hr_solver`) — the public façade. Every method takes user units, converts to SI, delegates, and returns user units. This is the only layer that should do unit conversion.
- `hr_factors.py` (`hr_factors`) — SI-only core. Computes partial Huang–Rhys factors `S_k = ω_k ΔQ_k² / 2ħ` from either forces (`compute_hrf_forces`) or GS↔ES displacements (`compute_hrf_dis`, with minimum-image wrapping through the cell). Owns the default IUPAC atomic-mass table.
- `lineshape.py` (`lineshape`) — SI-only core. Builds the time-domain generating function `G(t,T)` from `{ω_k, S_k}` and transforms it, either by Simpson integration (`compute_lineshape_numerical_integration`) or FFT (`compute_lineshape_fft`, faster, requires a symmetric energy range).

**1D configuration-coordinate path** — `config_coord_1d_solver.py` (`config_coord_1d_solver`), fully self-contained: Franck–Condon overlap integrals between two displaced harmonic oscillators via the Ruhoff recurrence relations, then a Gaussian/Lorentzian-broadened FC lineshape. Used when only two effective frequencies and a ΔQ are available rather than a full phonon spectrum.

Typical call chain: `compute_hrf_dis`/`compute_hrf_forces` → `hrf_dict` → `compute_spectral_density` and/or `compute_lineshape_fft` → `compute_spectrum(ezpl, spectrum_type="PL"|"Abs")`.

### Conventions that matter

- **Unit boundary.** `hr_factors` and `lineshape` assume strict SI internally (rad/s, m, J, kg); only masses are passed in amu and converted on entry. `hr_solver` and `config_coord_1d_solver` are the user-facing layer and speak THz / Å / eV / meV / K. Keep conversions at that boundary rather than scattering factors into the cores.
- **Data passed as dicts.** HR results travel as `{"freqs": ndarray (rad/s), "hr_factors": ndarray}`; lineshapes as the tuple `(energy_axis_meV, A)`. Several `hr_solver` methods also cache the last lineshape on `self.lineshape`.
- **Mode indexing.** `hr_factors` drops the first `nom_translational = 3` modes (acoustic/translational) plus any leading imaginary modes, zero-filling their HR factors, and reshapes `modes` to a square `(3N, 3N)` matrix — the input must be the complete Γ-point mode set for the same supercell as the coordinates/forces.
- **Broadening.** `sigma` is always a two-element `[σ_at_min_ω, σ_at_max_ω]` in meV, linearly interpolated across the phonon range (`f_sigma`, duplicated in both `hr_solver` and `lineshape`).
- Docstrings are NumPy-style with LaTeX math and are rendered by `doc/documentation.rst` via `sphinx.ext.autodoc` + `napoleon`; new public functions should follow the same format and be added to that file.

## Validation suite (`validation/`, formerly `testing/`)

> **Rule: whenever a VOF method or fit is added or changed, update `validation/METHODS.md`** (its section for that
> method and the comparison table in section 7), as well as this file and `validation/README.md`. `METHODS.md` is
> the tutorial-style explanation of every method: physics and derivation, how to run it, how to read its
> output, measured accuracy and when to use it.

A local harness for the **1D CCD** path on the NV-center data in `doc/tutorials/002_nv_diamond_1d_ccd/`.
It compares three independent routes to the Franck-Condon matrix `F_ij` and the VOF (`bulid_fc_lsp`
lineshape) built from it: `pypl` (PyPL's Ruhoff recurrences), `hermite` (quadrature of analytic
harmonic eigenfunctions) and `fd` (finite-difference Schrödinger solve on fitted CCD wells). Not
installed; CI runs `check_1d_ccd.py`. Renamed from `testing/` on 2026-10-16 (env var now `PYPL_VALIDATION_OUTPUT_ROOT`,
`harness.VALIDATION_DIR`); old outputs keep their old paths in their index.md / logs. `validation/README.md`
is the user-facing description; this section is the working memory.

```
validation/
├── harness/                          validation-only helpers, never moved into pypl/
│   ├── __init__.py                     puts the repo root and validation/ on sys.path
│   ├── ccd_data.py                     build_ccd(): the tutorial-002 CCD; samples(), center(), normalize()
│   ├── fits.py                         fit specs parabola | polyN | spline | morse -> make_potential()
│   ├── methods.py                      VOFParams, add_vof_args, METHODS, compute_vof, compare_results
│   ├── output.py                       OutputDir, text_table, markdown_table
│   └── plots.py                        style, colours, vof_comparison_figure, fc_error_figure
├── check_1d_ccd.py                   81 regression + physics checks (diamond + MgO), ~12 s
├── run_1d_ccd.py                     tutorial 002 as a script (unmodified config_coord_1d_solver)
├── compare_vof.py                    quick PASS/FAIL method comparison, ~30 s
├── plot_ccd_wavefunctions.py         CCD + wavefunctions (--wf-method) + VOF panel (--vof-method)
├── vof_report.py                     the thorough seven-section report, ~8 min
├── run_all.py                        master suite: every method/fit/variant from one TOML config, ~20 min
├── run_all.toml                      the commented default config (documents every key), NV- in diamond
├── run_all_nv_mgo.toml               the NV- in MgO suite config (suite name nv_mgo, dataset nv_mgo)
├── data/nv_mgo_west/             how the MgO data was collected on the HPC (README, collector); the data is doc/tutorials/003_nv_mgo_isc/
├── README.md
├── METHODS.md                        tutorial overview of every VOF method (update with every method change)
└── outputs/<script>/<run>/           every run's results (committed on vof-extensions); _archive/ = replaced runs; _legacy_flat/ = old flat files
```

### Docs pages and naming (2026-10-16)

- `doc/validation.rst` (*Validation and testing*): layers, terminology, datasets (with the method of each state's energies),
  every pytest test and validation check with assertion / tolerance / why, scripts, metrics, statuses, outputs, reference
  results, solver findings 1-7, how to extend. `doc/supporting_information.rst`: S1 units ... S12 metrics (every formula:
  Ruhoff recurrences with a..e and F00, Hermite recurrence, Simpson, FD stencils with truncation terms, fits incl. spline
  DOF / end conditions / tails, lineshape and sigma law, cubic PT, ISC, multi-mode path). Both in the index toctree; docs
  build without warnings. Update both whenever a method, fit or test changes.
- **Naming rule (user, 2026-10-16): never call the data points "DFT energies/images".** They are *computed energies* /
  *images*; each `ccd_data.Dataset` has `energy_methods` {"gs","es"} and `energy_zero` (diamond: DFT (PBE) / Delta-SCF
  constrained occupations (PBE); MgO: spin-flip TDDFT@DDH root 4 / spin-conserving TDDFT@DDH root 1, zero = 3A2 DDH
  minimum, PBE0 with 0.34 exact exchange). `ccd_data.ENERGY_METHODS`, `energy_label(state)` give legend text; every
  index.md / log.txt now has the Dataset block (diamond: log.txt only on the console side). CSV column `E_DFT_eV` is now
  `E_eV`. Fixed: `ccd_data` read the old `PYPL_TESTING_DATASET`; now `PYPL_VALIDATION_DATASET` (old name still accepted).

### The package extension (`pypl/config_coord_1d_vof_solver.py`, `config_coord_1d_potential.py`, `isc.py`)

Merged from `testing/pypl_ext/` on 2026-10-16. The user's requirement: **no existing PyPL function or
name may change behaviour; every tutorial must work exactly as before** (pinned by `tests/test_legacy_tutorials.py`). So:

- `config_coord_1d_vof_solver` *subclasses* `config_coord_1d_solver` and only **adds** methods. It has
  no `__init__` and overrides nothing; the new methods fill the same `fc_ints` / `ni` / `nf`, so the
  inherited `bulid_fc_lsp` / `compute_spectrum` work on them unchanged. `check_1d_ccd.py` enforces it:
  no overridden attribute, and inherited `compute_franck_condon_integrals`, `bulid_fc_lsp` and
  `compute_spectrum` are **bit-identical** (`np.array_equal`) on subclass and parent.
- PyPL style throughout: lowercase class names, NumPy docstrings with `$...$` LaTeX, user units (meV,
  Å·amu^½, eV), attributes set by `compute_*` methods. They import only `pypl`, numpy and scipy.
- Also added on the merge: `config_coord_1d_vof_solver.huang_rhys(omega, delta_q)` (tutorial-002 expression),
  `pypl/isc.py` (`compute_isc_rate(soc, vof, soc_unit="GHz"|"meV"|"ueV"|"cm-1")` = 4 pi / hbar |soc|^2 F in Hz,
  `compute_isc_lifetime`), `utils.parse_west_summary(fileName)` (plep/eigenval x ry2ev, eV), `pypl/__init__.py`
  (`__version__ = "1.1.0"`, read by pyproject and doc/conf.py). All in `doc/documentation.rst`; docs build without
  warnings. Only edits to old modules: raw strings for two printed `\Delta` / `\hbar` labels (same text; Python 3.12
  SyntaxWarning) and a stray `w` line removed from the `compute_spectral_density` docstring.

`config_coord_1d_vof_solver` adds: `harmonic_alpha(omega)`, `harmonic_wavefunctions(q, nmax, omega, center)`
(all levels by the stable normalized three-term recurrence; the old explicit 2^n n! form overflowed above n ~ 170),
`harmonic_wavefunction(q, n, omega,
center)`, `classical_turning_point(n, omega)`, `compute_franck_condon_integrals_hermite(ni, nf,
npts=100001)` (sets `hermite_grid`), `solve_schrodinger_fd(potential, q, n, stencil=5)` and
`compute_franck_condon_integrals_fd(ni, nf, potential_i, potential_f, q_i=None, q_f=0.0, npts=20001,
stencil=5, margin=2.0)`, plus the **optional cubic anharmonic correction**:
`anharmonic_strength(b, omega)` (xi = b/(hbar W) (hbar/2W)^1.5 = b / (omega eV) * (2 alpha)^-1.5),
`cubic_correction_coefficients(n, xi)` ((n, n+3) matrix C, corrected states = C @ psi_0..n+2),
`perturbed_wavefunction(q, n, omega, xi, center)` and
`compute_franck_condon_integrals_cubic(ni, nf, b_i, b_f, q_i=None, q_f=0.0, npts=100001, normalize=False,
warn_threshold=0.1)`. The last uses the Hermite Simpson grid (`_simpson_grid`, shared with the hermite method)
and flips both b when q_f < q_i (mirrored frame). It sets `cubic_xi_i/f`, `cubic_norms_i/f`, `cubic_grid`,
`cubic_wavefunctions_i/f` and `cubic_info` (b in the frame, xi, max mixing, `n_valid_i/f`, warnings), and it
prints a warning for |xi| > 0.1 and, separately, for mixing 3|xi|(n+1)^1.5 > 0.1 inside the basis, which is the
one that fires here (n_valid ~ 9 / 8). The FD method sets `fd_grid` (solver frame x), `fd_q` (physical Q),
`fd_energies_i/f` (eV), `fd_wavefunctions_i/f` and `fd_info` (npts, h, stencil, margin, frame, ranges,
turning points, wall, orthonormality, ladder_i/f). Module level: `grid_spacing(q)`, `KINETIC` =
ħ²/2 in eV·amu·Å², `FD_STENCILS`, `FD_STENCIL_NAMES`.

`config_coord_1d_potential(q, energies, q0)` has `fit_parabola(index)` (minimum fixed at q0, curve_fit
exactly as the notebook), `fit_cubic(index, curvature)` (a(Q-Q0)^2 + b(Q-Q0)^3 + c, a fixed when `curvature` is given;
`params` a, b, c, a_fixed, turnover; `phonon_energy` = the harmonic 2a value; unbounded beyond the turnover), `fit_polynomial(degree, index)`, `fit_morse(index)` (both signs of b tried),
`fit_spline(index, tail_curvature, order=3, bc=None, smoothing=None, tails="harmonic")` (spline of degree 2-5 with
continuous harmonic tails, or `tails="smooth"`: C2 Taylor tails continuing the spline's own value / slope / curvature,
spec suffix `-smooth`, flag `--spline-tails`, bound only if V'' > 0 at both ends, not with clamped; `bc` natural / not-a-knot / clamped, natural and clamped for odd order only; `smoothing` = target rms in
eV for a FITPACK smoothing spline; `params` records order, bc, n_knots, slope / curvature jumps and the C0/C1/C2
junction; `phonon_energy` uses the spline's analytic V''; `SPLINE_ORDER_NAMES` names the degrees), and `__call__`, `residuals`, `rms`, `minimum`, `phonon_energy`, `describe`, plus the attributes
`kind`, `index`, `params`, `bounded`.

**Numerical gotchas, each found the hard way:**
- **Grid spacing.** Always use `grid_spacing(q) = (q[-1]-q[0])/(len(q)-1)`, never `q[1]-q[0]`, which
  loses ~5 digits to cancellation at |Q| ~ 6. That slip gave the old Hermite method a hidden ~1e-11
  normalization error; fixing it took the Hermite-vs-PyPL VOF agreement from 4.2e-12 to 7e-16.
- **Eigensolver.** Use sparse shift-invert `eigsh` for the FD Hamiltonian. `scipy.linalg.eig_banded`
  forms an N x N transform and hangs at 5001 points.
- **Bit-identity of omega.** `phonon_energy()` deliberately uses the notebook's exact expression,
  `sqrt(2a eV/(1e-10**2 amu) hbar^2/eV^2) * 1000`, so the frequencies are bit-identical to the tutorial.
  An algebraically equal form moved PyPL's float64 drift from (49, 48) to (49, 56). The same applies to
  `huang_rhys` (`1e-10**2`, not `1e-20`).
- **Sign convention.** FD eigenvectors are signed so the right-hand tail is positive, which is the
  Hermite convention. The solver frame puts the final minimum at 0 and the initial one at -ΔQ, hence
  `x = sign(Q_f - Q_i)(Q - Q_f)`, mirrored for emission. Physical-frame overlaps need `(-1)^(i+j)`.
- **Summation in checks.** Sum large-N overlaps with `np.sum`, not a BLAS matrix product, whose rounding
  over 2e5 terms is ~1e-12.

### `harness/`

- `ccd_data.build_ccd()` builds the wells as `config_coord_1d_potential` parabolas fitted to
  `GS_FIT_SLICE = [:5]` (GS images 1-5 of 15) and `ES_FIT_SLICE = [4:]` (ES images 5-11 of 11). It returns
  `delta_q`, `gs/es_phonon`, `gs/es_hrf`, `gs/es_params` `(a, c)`, `gs/es_potential`, the samples
  `gs/es_1d_coord`, `gs/es_energies`, `gs/es_fit_energies` and `n_gs/es_images`. Reference values (pinned in
  `check_1d_ccd.REFERENCE`): ΔQ 0.695732 amu^½ Å, GS 59.08890 meV (HRF 3.4211), ES 64.31191 meV (HRF 3.7235).
- `fits`: a spec string (`parabola`, `polyN`, `spline`, `morse`) becomes a fitted potential through
  `make_potential(ccd, state, spec)`. `parabola` returns `ccd["<state>_potential"]` itself, and `spline`
  takes its tail curvature from that parabola. `check_bounded` raises on an unbounded fit.
  `REPORT_FITS` = parabola, poly2, poly4, poly6, poly8, morse, spline, spline3-notaknot, spline3-clamped,
  spline4-notaknot, spline5-clamped, spline3-s0.05. Spline specs: `spline[K][-natural|-notaknot|-clamped|-sRMS]`
  (RMS in meV); `canonical()` gives one spelling (`spline3-natural` -> `spline`, `spline4` -> `spline4-notaknot`)
  and `spline_spec(order, bc, smoothing_mev)` builds one. `methods.resolve_fd_fit` turns `--fd-fit spline`
  plus `--spline-order` / `--spline-bc` / `--spline-smoothing` into the canonical spec and rejects impossible
  combinations; every run name and label uses that canonical spec.
- Cubic: spec `cubic[-near][-free]` (default `cubic`: a fixed to the parabola, b and c fitted to all images).
  `VOFParams.cubic_fit` / `cubic_normalize` come from `--cubic-points {all,parabola}`, `--cubic-free-a` and
  `--cubic-normalize`. `METHODS['cubic']` ("Hermite + cubic PT") has `anharmonic=True`, so `compare_vof.py`'s
  default tests and report section 05 exclude it; `is_anharmonic(method, params)` drives the NOTEs. With
  `--fd-fit cubic`, fd solves the same cubic exactly; `check_bounded(ccd, spec, q_range)` allows a cubic only
  if its turnover lies outside the grid (`methods.fd_q_range`). compare_vof recognizes `--reference fd --fd-fit
  cubic --test cubic` as a like-for-like test of the PT truncation error. Run names carry `_cubicfit-<spec>`.
- `methods`: `VOFParams` holds every parameter, with tutorial defaults plus `npts`, `fd_npts`, `fd_stencil`,
  `fd_margin`, `fd_fit`; `params.orders(direction)` swaps the basis sizes for absorption.
  `add_vof_args` / `params_from_args` give every script the same flags. Each `METHODS` entry has a
  `run(solver, ni, nf, params, ccd, direction)` that fills a `config_coord_1d_vof_solver` and returns
  `notes`, i.e. every setting it used. pypl / hermite notes include the parabola fit (images, Q range, a,
  rms) that their omega comes from. `compute_vof(..., log=)` returns a `VOFResult`: solver, lineshape,
  notes, seconds, fit, closure, integral, centroid_mev, sideband peak. `compare_results(ref, tst)`
  returns the statistics, and `COMPARISON_COLUMNS` / `CSV_KEYS` are its table layouts.
- `output.OutputDir(script, run, title, verbose, keep_figures)` manages `outputs/<script>/<run>/`. The run
  name is sanitized. A re-run under the same name starts empty, but **nothing is ever deleted**:
  `archive_previous(path)` moves the old run to `outputs/_archive/<same path>/<YYYYmmdd-HHMMSS>/`, and
  `run_all.py` does the same for a whole suite (`_archive/suite/<name>/<stamp>/`). Never reintroduce
  `shutil.rmtree` on outputs (the user lost a suite to it on 2026-10-16). It provides
  `log` (console + log.txt), `detail` (log.txt only unless -v), `section`, `save_fig` / `save_csv` /
  `save_table` / `save_json` (each takes a description), `add_summary`, and `close()`, which writes
  index.md listing every file with its description.

### Scripts

| Script | Output directory | Contents |
| --- | --- | --- |
| `check_1d_ccd.py [-v]` | none | 75 checks (incl. 11 on spline options, 9 on the cubic correction, 16 on large-n Hermite functions and the MgO loader / windows / orientation; it always pins nv_diamond for the rest); exit 0 = pass. Tutorial numbers, FC closure, identity at ΔQ = 0, `F_ij = (-1)^(i+j) F_ji` at equal ω, lineshape and Boltzmann sanity, basis convergence, spectra, plus the extension contract above and quick hermite / fd / fit checks. **Update `REFERENCE` deliberately when a physics change is intended; do not widen tolerances.** |
| `run_1d_ccd.py` | `run_1d_ccd/<tag>/` | Tutorial 002 with the unmodified solver: data/ (summary, PL, Abs CSV) and figures/ (CCD, spectra). |
| `compare_vof.py` | `compare_vof/<ref>_vs_<tests>[_fit-X][_tag]/` | `--reference` (pypl) vs `--test` (default: all others). figures/vof_comparison.png (overlay, signed error, log error, settings panel), figures/fc_error_maps.png, data/vof_<direction>.csv, data/summary.csv. Exit 0 only if all pass. A non-parabola fit adds a NOTE that the FAIL is anharmonicity. |
| `plot_ccd_wavefunctions.py` | `plot_ccd_wavefunctions/wf-<m>_vof-<ms>[_fit-X][_tag]/`; `_fit-X` is added whenever fd is used, parabola included, and the title, legends and index.md name the fit and its hbar w | Left: DFT images, fits and wavefunctions (`--wf-method hermite` at E_min+(v+½)ħω, or `fd` at its own eigenvalues). The drawn functions are checked against the solver's F_ij (hermite 1e-10, measured 1.7e-16; fd 1e-8, measured 1.6e-12) and the script refuses to plot on a mismatch, except for an anharmonic fit, which is reported. Right: the VOF of each `--vof-method`. Flags: --n-levels, --density, --wf-scale, --alpha, --well, --no-vof. |
| `vof_report.py` | `vof_report/<tag or default>/` | The full comparison, in six sections (below). Flags: --fits, --reference, --skip-convergence, the VOF flags, --tag, -v. |

`vof_report.py` sections. Each holds PNGs and the CSVs behind them, and index.md carries the headline
tables:
1. `01_potential_fits`: `fits_overview.png`. Per well, four panels: the images (filled = used by the
   parabola) with every fit; every fit over the full FD-grid Q range (`Report.fd_q_range()`, -4.90..5.23,
   with the DFT points, walls and V at the walls printed); the residual at each image; and V_fit - V_parabola out to the turning points of the highest
   level (symlog, sampled range shaded, unbounded fits dashed). Also `fit_residual_rms.png`,
   `fit_summary.csv`, `residuals_<state>.csv`.
2. `02_fd_grid`: `fd_grid_<direction>.png`. Both wells on the grid with selected levels at their
   eigenvalues (padding beyond the turning points in gold, walls in black), a quarter-wavelength of the
   highest level point by point at the well centre, and log |ψ| decaying to the walls. Also
   `grid_summary.csv`, and `fd_convergence.png/.csv` (F_ij vs Hermite, ladder error and time vs h, 3- and
   5-point stencils, h^2 / h^4 guides; every measured h is a labelled tick with the grid size on a top axis, and
   the fitted log-log slope is in each legend and in the CSV: measured 2.00 / 4.00 for both the F_ij and
   ladder errors; time scales roughly linearly in N but is machine-noisy; skip with --skip-convergence).
3. `03_wavefunctions`: `wavefunctions_<state>.png`. Levels 0, 1, 2, 5, 10, 20, 30, 40 and the last, Hermite
   vs FD on every bounded fit, with the difference and each fit's max |diff| printed. Also
   `wavefunction_errors.png/.csv` (per level: infidelity 1-|<χ|ψ>|, max |ψ-χ|/max|χ|, ladder deviation).
4. `04_fc_matrix`: `fc_matrix.png` (log |F_ij| of the reference), `fc_error_maps.png`,
   `fc_error_by_level.png` (max_j |dF_ij| vs i), `fc_<method>_<direction>.csv` (every matrix),
   `fc_errors.csv`.
5. `05_vof_methods`: `vof_methods.png` (every method on the parabola wells vs the reference, with the
   settings panel), `vof_<direction>.csv`, `summary.csv`. Section 05 sets the report's exit code.
6. `06_vof_fits`: `vof_fits.png` (FD on every bounded fit vs the harmonic reference), `level_ladders.png`
   (E_v - E_0 - vħω per fit), `vof_<direction>.csv`, `summary.csv` (dF, dVOF, F_00², centroid shift).
7. `07_cubic_perturbation` (skip with --skip-cubic): `cubic_fits.png/.csv` (all four cubic specs: b, xi,
   turnover, max mixing, n_valid); `validation.png` (b scaled 1..1/16: PT-vs-exact error slope 2, effect slope 1,
   F_ij block i, j < 8 and VOF; norms vs 1 + sum c^2; mixing and infidelity per level vs exact eigenstates);
   `wavefunctions_<state>.png` (harmonic / PT / exact fd-on-cubic + differences); `fc_error_maps.png` and
   `vof_cubic.png` against the exact cubic VOF; `summary.csv`, `validation_b_scaling.csv`.

Shared VOF flags (`harness.methods.add_vof_args`): `--temp --sigma --gamma --no-zpl --order-es --order-gs
--ene-range --resol --npts --fd-npts --fd-stencil --fd-margin --fd-fit`. The old `--vof-temp`,
`--vof-sigma`, `--vof-zpl-gamma` and `--no-vof-zpl` spellings still work. Defaults are tutorial 002: 5 K,
sigma [8, 25] meV, Lorentzian ZPL with gamma 0.3 meV, 50 ES x 60 GS levels (swapped for absorption), and
an axis from -200 to 1200 meV with 1401 points.

### Master suite (`run_all.py` + `run_all.toml`)

`run_all.py` reads a TOML config (default `validation/run_all.toml`, with comments on every key). It deep-merges the
config over `DEFAULTS`, rejecting unknown keys and invalid fits; `--set section.key=value` overrides a single key
(the value is TOML), `--name`, `--only`, `--jobs` and `--dry-run` are flags. From the config it builds `Step`s,
each one run as a separate subprocess of an existing script:
- checks;
- tutorial;
- compare: harmonic like for like; fd per fit vs the reference; each cubic variant vs the reference and vs fd on
  its own cubic (`cubic_spec(variant)`);
- plots: hermite; fd on parabola plus every fit; cubic per variant;
- report.

fd fits are canonicalized and bound-checked first (`resolve_fits`, using `fd_q_range`); unbounded ones are
skipped with the reason. Subprocesses get `PYPL_VALIDATION_OUTPUT_ROOT=outputs/suite/<name>/`, which
`harness.output.OUTPUT_ROOT` honours, so every script writes `suite/<name>/<script>/<run>/`. The runner parses
each step's "Outputs in ... (see index.md)" line for its directory and `data/summary.csv` for the comparison
table. Status is PASS / OK / DIFFERS (expected) (an `expect="differ"` step with a FAIL verdict) / FAIL / ERROR
(a traceback, or a non-zero exit with no output). `index.md` gives the physics line, skipped fits, the method
comparison (kinds in `KIND_ORDER`: like for like, PT vs exact cubic, anharmonic vs harmonic), every step with a
link, and every command. `logs/` holds each step's console output. A smoke run with a small basis and
`--jobs 2` took 3 min; the default config is 28 steps.

### Current status (2026-10-15, after the restructure)

All runs below used `C:\Python312\python.exe` at the defaults. `check_1d_ccd.py` passes 59/59.
`compare_vof.py`, `run_1d_ccd.py`, `plot_ccd_wavefunctions.py` (all modes) and `vof_report.py` all
run cleanly, and black is clean on every testing file.

| Comparison (parabola wells unless noted) | max abs dF_ij | max dVOF / peak |
| --- | --- | --- |
| hermite vs pypl | 1.70e-8 (em.), 1.37e-8 (abs.) | 7.0e-16 |
| fd vs pypl | 1.80e-8 (em.), 1.51e-8 (abs.) | 1.3e-12 / 9.5e-13 |
| fd vs hermite | 3.5e-9 | ~1e-12 |
| fd, 3-point stencil, vs pypl | 5.5e-5 | 1.2e-6 |
| fd on spline vs pypl (anharmonic) | 1.8e-2 | 2.1e-2, mostly at the ZPL |
| fd on poly2 vs pypl (anharmonic) | 3.4e-3 | 2.4e-3 |
| fd on Morse vs pypl (anharmonic) | 4.0e-1 | 1.9e-2 |
| fd on poly6 vs pypl (anharmonic, extrapolated walls) | 7.1e-1 | 7.0e-2 |

- **The ~1.7e-8 F_ij residual is the PyPL recurrence's own float64 drift**, which sits only at the
  highest levels (`fc_error_by_level.png`). Within i, j < 30 the differences are 5e-11 (hermite) and
  6e-10 (fd).
- **FD grid at the defaults:** 20001 points, h = 5.07e-4, 302-318 points per shortest wavelength, wall
  amplitude ~4e-16, ladder within 2e-6 meV. The wavefunctions match Hermite to an infidelity of ~1e-15 and
  max |dψ| / max |χ| of 1.6e-8 at v = 59.
- **Fits (rms over all images).**
  - GS: parabola (1-5) 2.11 meV, poly2 0.15, poly6 0.024, spline 0 (exact), Morse 0.11.
  - ES: parabola (5-11) 0.84, poly2 0.077, poly6 0.006, spline 0, Morse 0.044.
  - Unbounded: poly4 and poly8, in both wells. Least squares leaves the leading coefficient negative
    (poly4 GS p_4 = -0.018, ES -0.0097; poly8 GS p_8 = -0.42, ES -3.25), so V -> -inf for large |Q|.
    The images span only ~1 (GS) / ~0.7 (ES) amu^0.5 A, so the highest coefficients are barely constrained
    and their sign is set by a few meV of softening plus noise. `describe()` now states the reason.
  - Out where the high levels live, poly6 is ~3e4 meV above the parabola. The spline's tails stay within
    5 meV and Morse within ~170 meV.
  - Wavefunction infidelity vs Hermite (max over levels): spline 8e-4 (GS) / 8e-5 (ES), poly2 ~2e-2,
    Morse 0.4-0.7, poly6 ~1 (unusable). So the spline is the only anharmonic fit that stays physical.
  - All six spline variants give the same wavefunctions, ladders and VOF to within a few percent of
    each other (VOF 2.1-2.2% of the peak from harmonic). **Correction (2026-10-16): most of that 2% is the harmonic tails' slope kink;
    with `-smooth` tails the not-a-knot spline is 1.9e-3 from harmonic (smooth vs harmonic tails: 2.0e-2, centroid
    -2.4 / -4.0 meV).** They differ only in hbar*omega at the minimum,
    GS 56.4-58.9 meV, where natural ends are the outlier because they force V'' = 0 at the first image,
    0.14 from the GS minimum.

- **Cubic correction (2026-10-15).**
  - b with a fixed, all images: GS +0.00916, ES -0.01265 eV/(amu^1.5 A^3); xi +1.03e-3 / -1.15e-3. The ES sign
    and the size depend strongly on the fit choice (near-min images: GS b = 0.042).
  - Validation: xi equals the SI formula, b = 0 gives the Hermite route to 4e-16, norms = 1 + sum c^2 to 4e-16.
    Against the exact cubic, the PT error slope is 1.99 (F_ij) / 2.00 (VOF) and the effect slope is 1.00; the
    wrong sign of b is 30x worse.
  - VOF: PT vs exact 3.3e-4 of the peak, harmonic vs exact 3.5e-2, so PT captures about 99% of the cubic
    effect on the VOF. But the mixing passes 0.1 above n ~ 9 / 8, norm - 1 reaches 4.1 at n = 59, and the
    high-level F_ij are wrong by up to 1.16: perturbative only where the intensity is (low levels at 5 K).

**Findings about the solver. None are fixed; nothing under `pypl/` was edited for this work.**

1. *Float64 drift in the Ruhoff recurrence.* F_ij is off by up to 1.7e-8, only at the highest indices
   of a 50 x 60 basis. The same recurrence in 60-digit arithmetic agrees with the Hermite quadrature to
   1e-12, so this is rounding, not physics, and it is invisible in any spectrum.
2. *ZPL replacement ignores the thermal weight.* With `zpl_lorentzian=True`, `bulid_fc_lsp` subtracts
   and adds `F_00²` without `p_0`. At 5 K `p_0 = 1`, so it is harmless. At 300 K `p_0 = 0.917` and the
   Lorentzian ZPL carries 0.02825 instead of `p_0 F_00²` = 0.02590, about 9% too much.
3. *`Integral check: 1.0101` is not a normalization error.* It is the 0.3 meV Lorentzian sampled on the
   1 meV grid. A Gaussian ZPL gives 1.000000, and the weights `p_i |F_ij|²` sum to 1.
4. The method is spelled `bulid_fc_lsp` in the package; scripts must use that spelling.
5. `bulid_fc_lsp` recomputes `energy_v` from `omega_i` / `omega_f`, so an anharmonic fit changes only the
   intensities, never the line positions.
6. *The Ruhoff recurrence is unstable at large Huang-Rhys factors.* On NV- in MgO (S ~ 36-85, 160 x 100 levels) row 0
   is exact (1e-15 vs hermite) but the error grows ~10x per initial level: > 1e-8 from i ~ 10-12, > 1 from i ~ 22-27,
   |F_ij| up to ~1e27-1e31 (an overlap cannot exceed 1). The 5 K VOF is still right (only i = 0 populated); any
   thermal or closure use of higher rows is wrong. Hence hermite is the MgO reference and pypl FAILs like for like.
7. *`bulid_fc_lsp` broadening grows linearly with |E| without a cap* (`sigma_factor = 160` meV); at the eV-scale
   sidebands of large-S systems it reaches hundreds of meV, so the MgO suite uses a constant sigma.

**The user's working preferences:**
- method-neutral names
- print and record every setting (grid, fit points, ranges, stencil, solver)
- organized, professional outputs with an index
- `parabola` stays the default and the testing fit; anharmonic fits are options, and their FAILs are
  labelled as anharmonicity
- the extension must never change existing PyPL behaviour

**NV- in MgO dataset (`--dataset nv_mgo`, data in `doc/tutorials/003_nv_mgo_isc/`, moved there
from `doc/nv_mgo_west/` on 2026-10-16; the transfer tarball was removed, it is in commit 7ed4613).** Collector layout: `{sc,sf}/Image_<q>/summary.json` (`plep.eigenval`, Ry), `qe/Image_<q>/data-file-schema.xml`
(34 images, from sc), `<state>/Image_<q>/inputs/` (logs), `reference/` (`ccd_nv_mgo.py` = the user's spec and the
reference for colours / labels; `sf_note.md`). `doc/tutorials/003_nv_mgo_isc/ccd_nv_mgo.py` is the user's plot script, edited to read
this layout (qe/ XMLs, local summary.json, script-relative paths, UTF-8 stdout, `SF_EXCLUDE = ["Image_+2.00"]`).
- E_state = qe etot + eigenval[root] x 13.605693 (3E: sc root 1; 1A1: sf root 4), zero = 3A2 (qe) minimum excluding
  Image_-1.00, exactly as ccd_nv_mgo.py. dQ = 3.097879 (Mg107 O108 N1, 216 atoms); every image is an exact linear
  interpolation of Image_+0.00 -> +1.00. ZPL 3E -> 1A1 = 554.4 meV (3E min 3.1376 eV, 1A1 min 2.5832 eV).
- **Layout the user asked for (2026-10-16):** Q = q dQ, **3E ('es', upper) at Q = 0 in red circles, 1A1 ('gs', lower)
  at Q = dQ in blue diamonds**, labels `$^3E$` / `$^1A_1$`. "Emission" = 3E -> 1A1, the transition of interest. The
  harness supports either orientation via `ccd["gs_center"/"es_center"]`, `ccd_data.center`, `ccd_data.q_span`; F_ij
  is orientation independent (checked: fd = hermite with 3E left).
- **Image windows the user specified (2026-10-16)**, q of Image_<q>, inclusive: parabola 3E -0.10..+0.45, 1A1
  +0.60..+1.15 (12 each); cubic fit (perturbative correction) -0.10..+1.15 for both (26 each; `ccd["<s>_cubic_index"]`);
  **every other fit (splines, Morse, polyN) uses ALL images** (3E 33, 1A1 30). The only image never used: sf
  Image_+2.00 (not a spin-flip run: wbse.out l_spin_flip = False, summary byte-identical to sc's). `fits.images()` /
  `fits.state_label()` / `ccd_data.image_note()` print the images in every legend, note, CSV and index.md.
- Harmonic (these parabolas): 3E 60.83 meV (S = 69.8), 1A1 31.55 meV (S = 36.2). sum F_0j^2 = 1 - 1e-6 at j = 55
  (emission) / 135 (absorption), hence `run_all_nv_mgo.toml`: 160 ES x 100 GS, axis -500..10000 meV at 5 meV,
  constant sigma 20 meV (PyPL's linear sigma would reach ~650 meV), Gaussian ZPL (F_00^2 ~ 1e-20), fd 40001 points,
  **reference = hermite** (see finding 6). Peak memory ~0.8 GB per step; machine has ~2.5 GB free, so jobs = 1.
- Data quirks: sc eig1/eig2 degenerate at +1.00 (E doublet), split 0.18 eV at Q = 0; 3E turns over (max ~4.81 eV near
  q = 0.95, second dip ~4.61 eV near q = 1.3); 1A1 (sf root 4) has a bump near q ~ 0 tracking 3E ~0.1 eV above it
  (likely a crossing / spin contamination; it sits where 3E v = 0 lives). sf ms0 is 0.15-0.17 eV above 3A2. XML_DIFFER
  images have identical SCF energies in sc and sf pw.out, so qe/ is exact for sf. `inputs/pwscf.xml` are nscf outputs
  with etot = 0: never read them.
- **Tutorial `doc/tutorials/pypl_003.ipynb`** (created 2026-10-16, executed, 34 cells; built and run by a scratch runner
  through jupyter_client because nbformat / nbclient are not installed; the `python3` kernelspec runs plain `python`,
  so put `C:\Python312` first on PATH). 3E -> 1A1 of NV- in MgO: harmonic VOF (PyPL Ruhoff, 60.8332 / 31.5511 meV,
  100 x 100, sigma [6, 8], gamma 0.3, 5 K, axis 0..2000 meV) vs smooth not-a-knot spline VOF (fd 40001 points +
  `build_fc_lsp_from_levels`), then Gamma = 4 pi / hbar |lambda_perp|^2 F(Delta) with lambda_perp(QDET) = 204.9925 GHz only
  (the user asked not to show the DFT lambda). Results at Delta_ad = 554.41 meV: harmonic F = 2.70e-7 /meV, Gamma 3.71
  MHz, tau 270 ns, FC shift 1145 meV; spline F = 3.97e-5 /meV, Gamma 545 MHz, tau 1.83 ns, FC shift 660 meV (x147).
  **The spline VOF is a single narrow peak: 92.5 % of |F_0j|^2 goes to 1A1 j = 18, 19, states with much of their weight in
  a 37 meV pocket of the sf root 4 near q ~ 0 (local min Q = -0.13, barrier Q = +0.25).** So the spline rate hinges on
  whether that shoulder is the true 1A1 PES; the notebook says so. 0-0 lines with ZPE: harmonic 579.3, spline 571.4 meV
  (F is evaluated from each VOF's own 0-0 line, as in the user's ISC script). Since the merge it uses only the package API
  (needs `pip install -e .`), is in `doc/tutorial.rst`, and was re-executed with identical results.
- **Recommended spline for MgO: `spline3-notaknot-smooth`** (user asked for smooth PESs, 2026-10-16). Clamped is
  wrong here: forcing the harmonic tail slope at the end images (3E Q = 6.2: +5.49 vs data +1.65) opens a 0.26-0.60 eV
  well with no data in the 3E gap Q = 4.65-6.20 and moves the 3E -> 1A1 VOF by 38-51 % of the peak. Smooth tails are
  C2, identical to not-a-knot inside the data, and move the VOF by 0.2 % (centroid 580.3 meV from the ZPL). In the
  dense images every spline cross-validates to ~2 meV. Its plot was added to the existing suite as
  `plot_ccd_wavefunctions/wf-fd_vof-hermite-fd_fit-spline3-notaknot-smooth/` (index.md has an "Added after this run"
  section); `run_all_nv_mgo.toml` now includes it. 80/80 checks.
- **Suite result (outputs/suite/nv_mgo/, 2026-10-16, 23 steps, 1053 s):** checks 75/75; fd vs hermite like for like
  dF 1.03e-8, dVOF 2e-9 (OK); pypl vs hermite dVOF 1e-14 but dF 2.7e27 (FAIL, finding 6); report FAIL only from pypl;
  section 07 skipped (no exact cubic). Skipped fits: poly6, poly8 (1A1 unbounded), cubic and the exact-fd checks of
  cubic PT (1A1 cubic turns over at Q = -1.70, cubic-free at -0.50; 3E at 3.24). Emission (3E -> 1A1) centroid shift
  vs harmonic (centroid ~1146 meV): splines -566..-591 meV, poly4 -186, Morse +42, poly2 -1043, cubic PT -34 (a fixed)
  / +20 (a fitted). Absorption shifts -2.1..-4.1 eV. The harmonic wells are far from the data over dQ (parabola rms over
  all images: 1A1 0.9 eV, 3E 5.6 eV with the old windows), so the anharmonic fits dominate.

**Next steps:**
1. DONE (2026-10-16): `config_coord_1d_vof_solver.build_fc_lsp_from_levels(eneaxis, temp, sigma, zpl_lorentzian, gamma,
   energies_i=None, energies_f=None)` places lines at actual level energies (default `fd_energies_*` x 1e3, meV) with
   Boltzmann weights from those levels; equals `bulid_fc_lsp` for harmonic ladders (checked at 300 K). The harness's
   compare / plot / report scripts still use `bulid_fc_lsp` (harmonic positions) for every method; switching them is open.
2. Image windows are per dataset constants (`MGO_PARABOLA_WINDOW`, `MGO_CUBIC_WINDOW` in `ccd_data.py`); harness
   flags to override them per run are still to do.
3. Fix solver finding 2 in `pypl/` (needs the user's approval because it changes PyPL behaviour), then pin
   it in `check_1d_ccd.py`.
4. DONE (2026-10-16): pypl_ext merged into `pypl/`.

### Environment notes

- Local interpreter is `C:\Python312\python.exe` (3.12.6). `numpy`, `scipy`, `matplotlib` and `h5py`
  were already present. **`pyyaml` was installed** because `pypl/utils.py` imports it, and **`black`
  26.10.1 was installed** on 2026-10-15; CI pins black 26.1.
- `pytest`, `sphinx`, `sphinx_rtd_theme` installed 2026-10-16, and pypl itself as an editable install
  (`pip install -e . --no-deps`), which the tutorials need.
- Everything (validation/, its outputs, doc/tutorials/003_nv_mgo_isc, CLAUDE.md) is tracked on `vof-extensions`; the user removed
  the testing ignore rules. `CLAUDE.md` is still listed in `.gitignore` but is tracked (added with -f).
- Heredocs in the Bash tool turn `\\n` into a real newline inside Python strings. For Python edits with
  backslashes, use the Edit / Write tools.

## Documentation and tutorials

`doc/tutorials/pypl_001.ipynb` (NV-center absorption/PL from full phonons, both the displacement and force routes) and `pypl_002.ipynb` (1D CCD) are the executable reference for the intended API usage; their input data lives in `doc/tutorials/001_*` and `002_*` and is real QE/Phonopy output.

Added 2026-10-16 (all in `doc/tutorial.rst`, executed and stored with outputs; nbsphinx does not re-run them):
- `pypl_003.ipynb`: 3E -> 1A1 ISC of NV- in MgO; data in `doc/tutorials/003_nv_mgo_isc/` (moved from `doc/nv_mgo_west/`, README.md added in the 001/002 style; the transfer tarball was dropped).
- `pypl_004.ipynb`: cubic perturbation theory on NV- diamond. Fits (a fixed), xi and n_valid, PT vs harmonic vs exact-fd
  wavefunctions (GS infidelity n=0: 5e-6 -> 2.5e-10), VOF (PT captures 99.1 % of the exact cubic effect), validity:
  b-scaling slopes 1.99 / 1.00, a universal infidelity-vs-mu_n plot on model wells (xi 1e-3..1e-2; harmonic ~mu^2, PT ~mu^3),
  NV- MgO as the failing case (|xi| ~ 0.01, n_valid = 1, turnovers inside the images). Supplementary S1-S4: ladder-operator
  derivation, coefficient table, E2 = -hbar w xi^2 (30n^2+30n+11) (checked against fd: ratio 1.002 at n=0, 1.039 at n=15),
  |Q_t| = l/(6|xi|), barrier = hbar w/(432 xi^2) (checked on the fits), n_valid formula.
- `pypl_004a.ipynb` (2026-10-16): the 004 analysis on NV- MgO. First figure = CCD + parabola/cubic/spline fits + turnovers
  (004's first figure is now also a CCD + fits plot; the per-well residual figure follows it). Cubic rms in window 108 / 40 meV,
  xi -0.0095 / +0.0089, n_valid 1, mixing at n = S 17 / 6; 1A1 PT vs cubic-in-a-box (wall Q = -1.3) and vs spline states;
  VOF: PT sum_j F_0j^2 = 4.52, exact cubic collapses (3E level -49 eV, grid past the turnover), spline peak 569 vs harmonic
  1145 meV; weight-averaged final-level mixing 6.1 (harmonic) / 2.3 (spline); b-scaling: PT captures 95/91/81/63/25 % at
  s = 1/64..1/4, s = 1/2 has a turnover inside the grid -> PT valid only for ~b/30; S4 table (E2(n=S) 13.5 / 3.2 hw,
  barrier 25 / 29 levels, turnovers inside the images).
- `pypl_005.ipynb`: runs pytest (-m "not slow") and check_1d_ccd.py via subprocess, method table (pypl/hermite/fd/cubic),
  fit table + VOF per bounded fit with build_fc_lsp_from_levels, compare_vof.py and plot_ccd_wavefunctions.py runs (shown
  figures), run_all --dry-run. Validation outputs go to a tempdir via PYPL_VALIDATION_OUTPUT_ROOT; paths are printed as <tmp>.
- Notebooks were generated from cell lists by a jupyter_client runner in the session scratchpad (nbformat/nbclient are not
  installed); edit the .ipynb directly from now on. The kernel needs pypl installed (`pip install -e .`).
- `validation/plot_ccd_wavefunctions.py`: the CCD and VOF legends now sit together below the panels (they used to cover
  the curves on MgO); existing figures in validation/outputs were not regenerated. Docs are built by `.github/workflows/doc.yml` on every push and deployed to GitHub Pages from `main` by `gh_pages.yml`; notebook rendering needs `nbsphinx` and `pandoc`, which are not project dependencies.
