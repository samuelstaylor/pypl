# validation/

Local test harness for PyPL's **1D configuration-coordinate (1D CCD)** path. It runs on two datasets
(`--dataset`): NV⁻ in diamond, 3E → 3A2, from `doc/tutorials/002_nv_diamond_1d_ccd/` (the default), and
NV⁻ in MgO, 3E → 1A1, from the WEST TDDFT-DDH data in `doc/tutorials/003_nv_mgo_isc/` (`nv_mgo`). It compares three independent routes to the
Franck–Condon matrix `F_ij` and the vibrational overlap function (VOF) built from it:

| Method | Where `F_ij` comes from |
| --- | --- |
| `pypl` | PyPL's own `config_coord_1d_solver.compute_franck_condon_integrals` (Ruhoff recurrences) |
| `hermite` | numerical quadrature of the analytic harmonic-oscillator eigenfunctions |
| `fd` | finite-difference solution of the Schrödinger equation on fitted CCD wells (`parabola`, `polyN`, `spline`, `morse`, `cubic`) |
| `cubic` | **optional, anharmonic:** harmonic eigenfunctions corrected to first order for a cubic term b Q³ of each well (same frequencies and HR factors) |

The classes it tests live in the package (`pypl/config_coord_1d_vof_solver.py`, `pypl/config_coord_1d_potential.py`);
the harness and scripts here are not installed. CI runs `check_1d_ccd.py`; the fast unit and regression tests of the
package are in `tests/` (`pytest`).

`doc/tutorials/pypl_005.ipynb` walks through the tests, the checks and these scripts, and compares fits and VOF
methods with the package API.

**[METHODS.md](METHODS.md)** is a tutorial-style overview of every method: the physics and derivations,
how to run it, how to read its output, measured accuracy, and when to use it.

## Layout

```
validation/
├── harness/                          validation-only helpers (not installed)
│   ├── ccd_data.py                     the datasets: nv_diamond (tutorial 002) and nv_mgo (doc/tutorials/003_nv_mgo_isc)
│   ├── fits.py                         fit specs: parabola, polyN, spline, morse
│   ├── methods.py                      VOFParams, shared CLI flags, METHODS, compute_vof, compare_results
│   ├── output.py                       OutputDir: one organized directory per run
│   └── plots.py                        shared style and comparison figures
├── check_1d_ccd.py                   regression + physics checks (81), exit 0 on success
├── run_1d_ccd.py                     tutorial 002 as a script (unmodified PyPL solver)
├── compare_vof.py                    quick PASS/FAIL comparison of F_ij methods
├── plot_ccd_wavefunctions.py         CCD + wavefunctions + VOF panel
├── vof_report.py                     the thorough, sectioned comparison report
├── run_all.py                        MASTER SUITE: runs every method / fit / variant with one parameter set
├── run_all.toml                      its commented configuration (every parameter and choice), NV⁻ in diamond
├── run_all_nv_mgo.toml               the same for NV⁻ in MgO, 3E -> 1A1 (suite name nv_mgo)
├── data/nv_mgo_west/                 how the MgO data was collected on the HPC (README, collector); data: doc/tutorials/003_nv_mgo_isc/
└── outputs/                          everything the scripts write; reruns archive to outputs/_archive/
```

## The extension never changes PyPL

`config_coord_1d_vof_solver` (in `pypl/`) is a **subclass** of PyPL's `config_coord_1d_solver` that only *adds*
methods (`compute_franck_condon_integrals_hermite`, `compute_franck_condon_integrals_fd`,
`solve_schrodinger_fd`, `harmonic_wavefunction`, `build_fc_lsp_from_levels`, ...). Every inherited method keeps its
name and behaviour, and the new methods fill the same attributes (`fc_ints`, `ni`, `nf`), so `bulid_fc_lsp`
and `compute_spectrum` work on their result unchanged. `check_1d_ccd.py` and `tests/` enforce this: no
inherited attribute is overridden, and `compute_franck_condon_integrals`, `bulid_fc_lsp` and
`compute_spectrum` give **bit-identical** results on the subclass and on the parent. `tests/test_legacy_tutorials.py`
also compares tutorials 001 and 002 with a reference generated from the original code.

Both classes follow PyPL's conventions: lowercase class names, NumPy docstrings with LaTeX, user units
(meV, Å·amu^½, eV).

```python
from pypl.config_coord_1d_potential import config_coord_1d_potential
from pypl.config_coord_1d_vof_solver import config_coord_1d_vof_solver

gs = config_coord_1d_potential(q_gs, e_gs, q0=0.0)        # DFT images: Q (amu^0.5 A), E (eV)
gs.fit_parabola(slice(0, 5))                              # or fit_polynomial(6), fit_spline(), fit_morse()
es = config_coord_1d_potential(q_es, e_es, q0=delta_q)
es.fit_parabola(slice(4, None))

solver = config_coord_1d_vof_solver(es.phonon_energy(), gs.phonon_energy(), delta_q)
solver.compute_franck_condon_integrals_fd(50, 60, es, gs, q_i=delta_q, q_f=0.0)   # emission ES -> GS
energy, vof = solver.bulid_fc_lsp(eneaxis, temp=5, sigma=[8, 25], gamma=0.3)       # inherited, unchanged
```

## Running

No `pip install .` needed: `harness` puts the repository root and `validation/` on `sys.path`.

### The master suite (start here)

`run_all.py` runs every method, every fit and every cubic variant with one set of parameters, as separate
steps:
- the regression checks;
- tutorial 002;
- every comparison;
- every CCD/wavefunction figure;
- the full report.

Everything goes to `validation/outputs/suite/<name>/`, one directory per step, exactly as each script would write
it alone. The suite's `index.md` gives the status, time and link of every step, plus a **method-comparison
table** for choosing a method.

```bash
python validation/run_all.py                                     # validation/run_all.toml, ~20 min
python validation/run_all.py --config my_params.toml             # a copy of run_all.toml with your choices
python validation/run_all.py --set physics.temp=300 --name 300K  # override any key from the command line
python validation/run_all.py --set fd.npts=40001 --set fd.stencil=3 --name fd-3pt
python validation/run_all.py --set 'fd.fits=["spline4", "morse"]' --only compare plots
python validation/run_all.py --dry-run                           # list every command without running
python validation/run_all.py --jobs 2                            # two steps in parallel (~1 GB each)
python validation/run_all.py --config validation/run_all_nv_mgo.toml   # NV- in MgO, 3E -> 1A1 -> outputs/suite/nv_mgo/
```

**NV⁻ in MgO.** `run_all_nv_mgo.toml` sets `suite.dataset = "nv_mgo"` and the suite name `nv_mgo`, so it writes
`outputs/suite/nv_mgo/` and never touches the diamond suite (`outputs/suite/default/`). A suite directory
records its dataset, and `run_all.py` refuses to replace one made from a different dataset. Its Huang-Rhys
factors (50 and 85) need 160 × 100 levels, a −0.5 … 10 eV axis and `hermite` as the reference, because PyPL's
recurrence is unstable there (see METHODS.md §2). The tutorial step does not apply (3E → 1A1 has no PL).
`doc/tutorials/pypl_003.ipynb` is the MgO tutorial: harmonic (PyPL, Ruhoff) vs smooth not-a-knot spline VOFs
of 3E → 1A1 and the ISC rate Γ(Δ) from each (λ⊥ from QDET). It uses only the package API
(`pypl.config_coord_1d_*`, `pypl.isc`, `pypl.utils.parse_west_summary`). The MgO CCD follows `doc/tutorials/003_nv_mgo_isc/reference/ccd_nv_mgo.py`: 3E (red) at Q = 0, 1A1 (blue) at Q = ΔQ. Image
windows (q of `Image_<q>`): parabola 3E −0.10 … +0.45 and 1A1 +0.60 … +1.15; cubic (perturbative correction)
−0.10 … +1.15; every other fit all images. Every legend, note and `index.md` names the images each fit used.

`run_all.toml` documents every key:

| Section | What it sets |
| --- | --- |
| `[suite]` | name, parallel jobs, which steps to run |
| `[physics]` | temperature, broadening, ZPL, basis sizes, energy axis |
| `[methods]` | the reference and the like-for-like methods |
| `[hermite]` | quadrature points |
| `[fd]` | grid points, stencil, margin, and the list of fits; `[[fd.splines]]` entries give extra splines by order / end conditions / smoothing |
| `[[cubic.variants]]` | cubic-correction variants (`points`, `free_a` = fixed or fitted a, `normalize`) |
| `[plots]` | wavefunctions drawn per well |
| `[report]` | convergence scan, section 07, fits |
| `[tolerances]` | what PASS means |

Misspelled keys and invalid fits are rejected before anything runs, and unbounded fits are skipped with the
reason. Steps comparing anharmonic methods are reported as **DIFFERS (expected)**, not as failures. The exit
status is 0 when nothing errored and every like-for-like step passed.

### Individual scripts

```bash
python validation/check_1d_ccd.py                  # 59 checks, ~5 s, quiet unless something fails (-v: all)
python validation/vof_report.py                    # the full report, ~8 min (--skip-convergence / --skip-cubic for less)
python validation/compare_vof.py                   # quick PASS/FAIL: pypl vs hermite and fd, ~30 s
python validation/compare_vof.py --test fd --fd-fit spline        # anharmonic wells (FAIL by design)
python validation/compare_vof.py --test fd --fd-fit spline --spline-order 5 --spline-bc clamped
python validation/compare_vof.py --test cubic                      # cubic perturbation vs harmonic (FAIL by design)
python validation/compare_vof.py --reference fd --fd-fit cubic --test cubic   # PT vs the exact cubic solution
python validation/plot_ccd_wavefunctions.py --wf-method cubic --vof-method pypl cubic
python validation/plot_ccd_wavefunctions.py --wf-method fd --vof-method pypl hermite fd
python validation/run_1d_ccd.py                    # tutorial 002 as a script
```

Every script takes `--tag NAME` (the run directory's name) and `-v` (print every method's settings,
which otherwise go only to `log.txt`). The VOF flags are shared by `vof_report.py`, `compare_vof.py`
and `plot_ccd_wavefunctions.py`: `--temp --sigma --gamma --no-zpl --order-es --order-gs --ene-range
--resol --npts --fd-npts --fd-stencil --fd-margin --fd-fit`.

## Spline options

`--fd-fit spline` is a natural cubic spline through every DFT image, continued outside them by the
harmonic parabola's curvature. Three flags configure it (or write the spec directly, e.g. `--fd-fit spline5-clamped`):

| Flag | Values | Meaning |
| --- | --- | --- |
| `--spline-order` | 2, 3 (default), 4, 5 | degree of the pieces: quadratic, cubic, quartic, quintic |
| `--spline-bc` | `natural`, `not-a-knot`, `clamped` | end conditions at the first and last image (below) |
| `--spline-smoothing RMS` | meV | a smoothing spline with rms residual ≤ RMS instead of exact interpolation; no end conditions apply |
| `--spline-tails` | `harmonic` (default), `smooth` | continuation past the end images: the harmonic parabola, or the spline's own C2 Taylor expansion (spec suffix `-smooth`) |

- `natural` (default for odd order): the higher derivatives vanish at the end images; V'' = 0 for a cubic, V''' = V'''' = 0 for a quintic. Odd order only.
- `not-a-knot` (default for even order): no end conditions; the knots next to the ends are dropped instead. Any order.
- `clamped`: slope (and for order 5, curvature) matched to the harmonic tails, so the spline joins its tails smoothly (C1 / C2) instead of with a kink. Odd order only.
- `-smooth` tails (any end condition but `clamped`): the tails continue the spline's own value, slope and curvature, so the PES is C2 everywhere while the spline inside the data is unchanged. Bound only if V'' > 0 at both end images. Recommended for NV⁻ in MgO: `--fd-fit spline3-notaknot-smooth` (see METHODS.md §5).

An even-order spline needs an odd number of end conditions, which cannot be split evenly between the two ends,
so `natural` and `clamped` are rejected for orders 2 and 4. Each spline's description in the logs states its
order, end conditions, number of knots, the continuity at the tails (C0/C1/C2) and the slope jump there. The
report compares six variants: cubic natural, not-a-knot and clamped, quartic not-a-knot, quintic clamped, and a
cubic smoothing spline with rms ≤ 0.05 meV.

Measured on this CCD, all six give the same physics. Their wavefunction infidelity against Hermite (7.6e-4 GS,
8e-5 ES), level ladders (2.1 / 1.0 meV) and VOF differences from the harmonic reference (2.1-2.2% of the peak)
agree to within a few percent of each other. Most of that 2% is the slope kink of the harmonic tails: with smooth
tails (`-smooth`) the not-a-knot spline is within 0.19% of the harmonic reference (METHODS.md §5.1). They differ mostly in hbar*omega at the minimum (GS 56.4-58.9,
ES 62.2-64.7 meV), because that is the curvature at a single point, the least well determined property of an
interpolant. `natural` is the outlier there: the GS minimum lies only 0.14 amu^0.5 A from the first image,
where natural ends force V'' = 0.

## Cubic anharmonic correction (method `cubic`)

An optional method: each well is fitted with V = a(Q−Q₀)² + b(Q−Q₀)³ + c, and the harmonic eigenfunctions are
corrected to first order in H′ = bQ³:

ψ̃ₙ = ψₙ + ξ[−3(n+1)^{3/2} ψₙ₊₁ − √((n+1)(n+2)(n+3))/3 ψₙ₊₃ + 3n^{3/2} ψₙ₋₁ + √(n(n−1)(n−2))/3 ψₙ₋₃],  ξ = (b/ħΩ)(ħ/2Ω)^{3/2}

The F_ij are integrated from these ψ̃ with the Hermite quadrature. Ω and the HR factors come from the parabola
as before, and the level energies stay harmonic, since the first-order shift of a cubic vanishes. When the method
isn't requested nothing changes.

| Flag | Default | Meaning |
| --- | --- | --- |
| `--cubic-points` | `all` | images b is fitted to: `all`, or `parabola` (the near-minimum ones) |
| `--cubic-free-a` | off | also fit a; by default a is fixed to the parabola's, so H₀ is exactly the harmonic Hamiltonian |
| `--cubic-normalize` | off | renormalize ψ̃ (by default the norm is 1 + O(ξ²)) |

Validated (`check_1d_ccd.py` and `vof_report.py` section 07):
- ξ matches an all-SI evaluation;
- b = 0 reproduces the Hermite route to 4e-16;
- the norms equal 1 + Σc² exactly;
- against the exact finite-difference solution of the same cubic (`--fd-fit cubic`), the error falls as b^2.00 while
  the effect falls as b^1.00;
- flipping the sign of b makes it 30× worse, so the mirror-frame sign is right.

**Validity.** ξ is only ~1e-3 here, far below the 0.1 warning threshold. But the mixing coefficients grow as
3|ξ|(n+1)^{3/2}: they pass 0.1 above n ≈ 9 and reach 1.4 at n = 59. The method therefore also warns when the
basis goes beyond that point, and reports the highest valid level. The VOF is still accurate, to 3e-4 of the
peak against the exact cubic solution, against an anharmonic effect of 3.5e-2. That is because at 5 K the
intensity sits in the low levels. The high-level F_ij themselves are not perturbative (they reach 1.16).

b also depends strongly on the fit choice. For the GS, b is 0.009 (all images, a fixed), 0.005 (a fitted) or 0.042
(near-minimum images). For the ES it even changes sign, because the near-harmonic data barely constrain it.

## Outputs

Each run writes one directory, `outputs/<script>/<run>/`. A re-run under the same name starts from an empty directory, but **nothing is ever deleted**: the previous run (or whole suite) is moved to `outputs/_archive/<same path>/<YYYYmmdd-HHMMSS>/`, stamped with the time that run last wrote a file.
A standalone run on another dataset (`--dataset nv_mgo`) writes `outputs/nv_mgo/<script>/<run>/` instead, and a
suite writes `outputs/suite/<name>/<script>/<run>/`, so diamond and MgO runs never overwrite each other. On MgO,
every `index.md` and `log.txt` opens with the dataset, its image selection and the role mapping, and every
figure carries a dataset caption:

| File | Content |
| --- | --- |
| `index.md` | **Start here.** Headline results and every file with a one-line description. |
| `log.txt` | The full console output, including every method's settings (grid, fit, stencil, ...). |
| `settings.json` | Every parameter of the run. |
| subdirectories | Figures (`.png`) and their data (`.csv`), grouped by topic. |

`vof_report.py` writes `outputs/vof_report/<tag or default>/` in six sections:

| Section | What it compares |
| --- | --- |
| `01_potential_fits/` | Every fit of the DFT images (parabola, polynomials of degree 2–8, spline, Morse): the fits near the DFT points, the fits over the **full Q range of the FD grid** (DFT points, sampled range and walls marked, V at the walls printed), the residual at each image, and V_fit − V_parabola out to the turning points of the highest level. Unbounded fits are dashed and flagged with the reason. |
| `02_fd_grid/` | The FD grid: both wells sampled, levels at their eigenvalues, the highest level point by point, its decay to the walls; plus F_ij / ladder error and time vs grid spacing for the 3- and 5-point stencils, with every grid size labelled and the measured convergence order (log–log slope) in the legend: 2.00 and 4.00. |
| `03_wavefunctions/` | Wavefunctions of every level: analytic Hermite vs FD on each fit, overlaid and differenced, with per-level infidelity and max error. |
| `04_fc_matrix/` | The reference F_ij, error maps of each method, error vs level, and every F_ij matrix as CSV. |
| `05_vof_methods/` | VOF of every method on the harmonic wells (the like-for-like test), with signed and log error. |
| `06_vof_fits/` | VOF of FD on every fit against the harmonic reference, and each fit's level ladder. |
| `07_cubic_perturbation/` | The cubic method: the four cubic-fit options (`cubic_fits.png`); validation (b-scaling vs the exact cubic, norms, mixing and infidelity per level); harmonic / first-order / exact wavefunctions per well; F_ij error maps and the VOF against the exact cubic solution. |

In sections 03 and 06 every non-parabola difference is the **anharmonicity of that fit**, not an error:
the harmonic methods only know the parabola near each minimum.

`outputs/_legacy_flat/` holds files written by the old flat layout; delete it when no longer needed.

## Formatting

Black at line length 120, as CI checks it (`validation/data/` is excluded in `pyproject.toml`):

```bash
python -m black --check --line-length 120 pypl tests validation
```
