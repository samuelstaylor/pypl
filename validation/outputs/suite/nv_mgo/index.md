# PyPL 1D CCD test suite: nv_mgo (NV- in MgO)

Started 2026-10-16T15:15:50, 23 steps in 1053 s. Configuration: `run_all_nv_mgo.toml` (copied here), resolved values in `suite_config.json`.

**Dataset.**

- dataset: nv_mgo (NV- in MgO), 3E -> 1A1 (WEST TDDFT-DDH)
- data: doc/nv_mgo_west
- roles: 'ES' = 3E (upper state, minimum at Q = 0), 'GS' = 1A1 (lower state, minimum at Q = delta_q); emission = ES -> GS = 3E -> 1A1
- Q = q delta_q, q from Image_<q>, as in reference/ccd_nv_mgo.py: 3E minimum at Q = 0, 1A1 minimum at Q = delta_q
- E = QE total energy (qe/) + WEST root (3E: sc root 1, 1A1: sf root 4) x 13.605693 eV/Ry, relative to the 3A2 (QE) minimum, as in ccd_nv_mgo.py
- image windows (q of Image_<q>, inclusive): parabola 3E q = -0.10..+0.45, 1A1 q = +0.60..+1.15; cubic a Q^2 + b Q^3 + c (perturbative correction) q = -0.10..+1.15 for both; every other fit (splines, Morse, polyN) ALL images
- all images: 3E all 33 sc images (q = -1.00..+2.00), 1A1 all 30 sf images except Image_+2.00 (not a spin-flip run, so not a 1A1 point)

**Physics:** T = 5 K, sigma = [20.0] meV, ZPL Gaussian, 160 ES x 100 GS levels, axis -500..10000 meV (2101 points). **Reference:** Hermite (quadrature). **fd:** 40001 points, 5-point stencil, margin 2. **Hermite:** 100001 points.

**Result: 2 step(s) failed or errored: harmonic methods like for like: hermite vs pypl, fd; full seven-section report.**

Skipped (no bound states on the fd grid):

- `poly6`: fit poly6 of the 1A1 images is unbounded below and has no bound states
- `poly8`: fit poly8 of the 1A1 images is unbounded below and has no bound states
- `cubic`: the 1A1 cubic turns over at Q = -1.70, inside the grid -6.68..10.23: no bound states in the box
- `cubic (exact fd check of cubic PT)`: the 1A1 cubic turns over at Q = -1.70, inside the grid -6.68..10.23: no bound states in the box
- `cubic-free (exact fd check of cubic PT)`: the 1A1 cubic turns over at Q = -0.50, inside the grid -6.68..10.23: no bound states in the box

## Method comparison

Largest difference over emission and absorption. **Like for like:** methods solving the same harmonic problem; they must agree within the tolerances (max dF_ij <= 1e-07, max dVOF / peak <= 1e-08). **PT vs exact cubic:** the error of first-order perturbation theory, against the exact fd solution of the same cubic. **Anharmonic vs harmonic:** how much each fit or correction changes the VOF from the harmonic reference (physics, not error), smallest first. See `testing/METHODS.md` for what each method does.

| kind | method | compared with | max dVOF / peak | max dF_ij | F_ij time (s) | verdict |
| --- | --- | --- | --- | --- | --- | --- |
| like for like | PyPL (Ruhoff) | Hermite (quadrature) | 9.62e-15 | 2.68e+27 | 0.03 | FAIL |
| like for like | Finite difference | Hermite (quadrature) | 1.95e-09 | 1.03e-08 | 11.46 | OK |
| anharmonic vs harmonic | Finite difference [Morse] | Hermite (quadrature) | 3.55e+00 | 6.39e-01 | 9.34 | differs (expected) |
| anharmonic vs harmonic | Hermite + cubic PT [cubic (a fixed, images q = -0.10..+1.15)] | Hermite (quadrature) | 4.79e+00 | 1.73e+02 | 1.55 | differs (expected) |
| anharmonic vs harmonic | Finite difference [polynomial deg. 4] | Hermite (quadrature) | 6.02e+00 | 6.45e-01 | 9.81 | differs (expected) |
| anharmonic vs harmonic | Finite difference [polynomial deg. 2] | Hermite (quadrature) | 7.62e+00 | 5.61e-01 | 12.01 | differs (expected) |
| anharmonic vs harmonic | Finite difference [quintic spline (clamped)] | Hermite (quadrature) | 7.66e+00 | 7.03e-01 | 10.49 | differs (expected) |
| anharmonic vs harmonic | Finite difference [quartic spline (not-a-knot)] | Hermite (quadrature) | 8.95e+00 | 7.73e-01 | 9.50 | differs (expected) |
| anharmonic vs harmonic | Finite difference [cubic smoothing spline (rms 0.05 meV)] | Hermite (quadrature) | 9.58e+00 | 7.79e-01 | 10.10 | differs (expected) |
| anharmonic vs harmonic | Hermite + cubic PT [cubic (a fitted, images q = -0.10..+1.15)] | Hermite (quadrature) | 9.77e+00 | 4.76e+02 | 1.66 | differs (expected) |
| anharmonic vs harmonic | Finite difference [cubic spline (natural)] | Hermite (quadrature) | 1.05e+01 | 8.64e-01 | 11.87 | differs (expected) |

## Steps

| # | step | status | time (s) | output |
| --- | --- | --- | --- | --- |
| 1 | checks: 75 regression and physics checks | PASS (75/75 checks passed) | 12.1 | `logs/01_checks_75-regression-and-physics-checks.log` |
| 2 | compare: harmonic methods like for like: hermite vs pypl, fd | FAIL | 32.4 | [`compare_vof/hermite_vs_pypl-fd/`](compare_vof/hermite_vs_pypl-fd/index.md) |
| 3 | compare: fd on cubic spline (natural) vs hermite | DIFFERS (expected) | 31.2 | [`compare_vof/hermite_vs_fd_fit-spline/`](compare_vof/hermite_vs_fd_fit-spline/index.md) |
| 4 | compare: fd on quintic spline (clamped) vs hermite | DIFFERS (expected) | 28.5 | [`compare_vof/hermite_vs_fd_fit-spline5-clamped/`](compare_vof/hermite_vs_fd_fit-spline5-clamped/index.md) |
| 5 | compare: fd on Morse vs hermite | DIFFERS (expected) | 26.2 | [`compare_vof/hermite_vs_fd_fit-morse/`](compare_vof/hermite_vs_fd_fit-morse/index.md) |
| 6 | compare: fd on polynomial deg. 2 vs hermite | DIFFERS (expected) | 31.6 | [`compare_vof/hermite_vs_fd_fit-poly2/`](compare_vof/hermite_vs_fd_fit-poly2/index.md) |
| 7 | compare: fd on polynomial deg. 4 vs hermite | DIFFERS (expected) | 27.1 | [`compare_vof/hermite_vs_fd_fit-poly4/`](compare_vof/hermite_vs_fd_fit-poly4/index.md) |
| 8 | compare: fd on quartic spline (not-a-knot) vs hermite | DIFFERS (expected) | 26.6 | [`compare_vof/hermite_vs_fd_fit-spline4-notaknot/`](compare_vof/hermite_vs_fd_fit-spline4-notaknot/index.md) |
| 9 | compare: fd on cubic smoothing spline (rms 0.05 meV) vs hermite | DIFFERS (expected) | 28.1 | [`compare_vof/hermite_vs_fd_fit-spline3-s0.05/`](compare_vof/hermite_vs_fd_fit-spline3-s0.05/index.md) |
| 10 | compare: cubic PT [cubic (a fixed, images q = -0.10..+1.15)] vs hermite | DIFFERS (expected) | 11.7 | [`compare_vof/hermite_vs_cubic_cubicfit-cubic/`](compare_vof/hermite_vs_cubic_cubicfit-cubic/index.md) |
| 11 | compare: cubic PT [cubic (a fitted, images q = -0.10..+1.15)] vs hermite | DIFFERS (expected) | 11.7 | [`compare_vof/hermite_vs_cubic_cubicfit-cubic-free/`](compare_vof/hermite_vs_cubic_cubicfit-cubic-free/index.md) |
| 12 | plots: Hermite wavefunctions, VOF of hermite, pypl, fd | OK | 33.9 | [`plot_ccd_wavefunctions/wf-hermite_vof-hermite-pypl-fd_fit-parabola/`](plot_ccd_wavefunctions/wf-hermite_vof-hermite-pypl-fd_fit-parabola/index.md) |
| 13 | plots: fd wavefunctions on parabola (near min.) | OK | 29.2 | [`plot_ccd_wavefunctions/wf-fd_vof-hermite-fd_fit-parabola/`](plot_ccd_wavefunctions/wf-fd_vof-hermite-fd_fit-parabola/index.md) |
| 14 | plots: fd wavefunctions on cubic spline (natural) | OK | 35.1 | [`plot_ccd_wavefunctions/wf-fd_vof-hermite-fd_fit-spline/`](plot_ccd_wavefunctions/wf-fd_vof-hermite-fd_fit-spline/index.md) |
| 15 | plots: fd wavefunctions on quintic spline (clamped) | OK | 32.8 | [`plot_ccd_wavefunctions/wf-fd_vof-hermite-fd_fit-spline5-clamped/`](plot_ccd_wavefunctions/wf-fd_vof-hermite-fd_fit-spline5-clamped/index.md) |
| 16 | plots: fd wavefunctions on Morse | OK | 30.3 | [`plot_ccd_wavefunctions/wf-fd_vof-hermite-fd_fit-morse/`](plot_ccd_wavefunctions/wf-fd_vof-hermite-fd_fit-morse/index.md) |
| 17 | plots: fd wavefunctions on polynomial deg. 2 | OK | 35.3 | [`plot_ccd_wavefunctions/wf-fd_vof-hermite-fd_fit-poly2/`](plot_ccd_wavefunctions/wf-fd_vof-hermite-fd_fit-poly2/index.md) |
| 18 | plots: fd wavefunctions on polynomial deg. 4 | OK | 27.7 | [`plot_ccd_wavefunctions/wf-fd_vof-hermite-fd_fit-poly4/`](plot_ccd_wavefunctions/wf-fd_vof-hermite-fd_fit-poly4/index.md) |
| 19 | plots: fd wavefunctions on quartic spline (not-a-knot) | OK | 33.9 | [`plot_ccd_wavefunctions/wf-fd_vof-hermite-fd_fit-spline4-notaknot/`](plot_ccd_wavefunctions/wf-fd_vof-hermite-fd_fit-spline4-notaknot/index.md) |
| 20 | plots: fd wavefunctions on cubic smoothing spline (rms 0.05 meV) | OK | 30.5 | [`plot_ccd_wavefunctions/wf-fd_vof-hermite-fd_fit-spline3-s0.05/`](plot_ccd_wavefunctions/wf-fd_vof-hermite-fd_fit-spline3-s0.05/index.md) |
| 21 | plots: cubic PT wavefunctions [cubic (a fixed, images q = -0.10..+1.15)] | OK | 10.9 | [`plot_ccd_wavefunctions/wf-cubic_vof-hermite-cubic_cubicfit-cubic/`](plot_ccd_wavefunctions/wf-cubic_vof-hermite-cubic_cubicfit-cubic/index.md) |
| 22 | plots: cubic PT wavefunctions [cubic (a fitted, images q = -0.10..+1.15)] | OK | 11.0 | [`plot_ccd_wavefunctions/wf-cubic_vof-hermite-cubic_cubicfit-cubic-free/`](plot_ccd_wavefunctions/wf-cubic_vof-hermite-cubic_cubicfit-cubic-free/index.md) |
| 23 | report: full seven-section report | FAIL | 474.8 | [`vof_report/default/`](vof_report/default/index.md) |

## Commands

Each step was run as (from the repository root, with `PYPL_TESTING_OUTPUT_ROOT` pointing at this directory and `PYPL_TESTING_DATASET=nv_mgo`):

```bash
python testing/check_1d_ccd.py -v
python testing/compare_vof.py --reference hermite --test pypl fd --dataset nv_mgo --temp 5.0 --sigma 20.0 --gamma 0.3 --no-zpl --order-es 160 --order-gs 100 --ene-range -500.0 10000.0 --resol 2101 --npts 100001 --fd-npts 40001 --fd-stencil 5 --fd-margin 2.0 --fd-fit parabola --tol-fc 1e-07 --tol-vof 1e-08
python testing/compare_vof.py --reference hermite --test fd --dataset nv_mgo --temp 5.0 --sigma 20.0 --gamma 0.3 --no-zpl --order-es 160 --order-gs 100 --ene-range -500.0 10000.0 --resol 2101 --npts 100001 --fd-npts 40001 --fd-stencil 5 --fd-margin 2.0 --fd-fit spline --tol-fc 1e-07 --tol-vof 1e-08
python testing/compare_vof.py --reference hermite --test fd --dataset nv_mgo --temp 5.0 --sigma 20.0 --gamma 0.3 --no-zpl --order-es 160 --order-gs 100 --ene-range -500.0 10000.0 --resol 2101 --npts 100001 --fd-npts 40001 --fd-stencil 5 --fd-margin 2.0 --fd-fit spline5-clamped --tol-fc 1e-07 --tol-vof 1e-08
python testing/compare_vof.py --reference hermite --test fd --dataset nv_mgo --temp 5.0 --sigma 20.0 --gamma 0.3 --no-zpl --order-es 160 --order-gs 100 --ene-range -500.0 10000.0 --resol 2101 --npts 100001 --fd-npts 40001 --fd-stencil 5 --fd-margin 2.0 --fd-fit morse --tol-fc 1e-07 --tol-vof 1e-08
python testing/compare_vof.py --reference hermite --test fd --dataset nv_mgo --temp 5.0 --sigma 20.0 --gamma 0.3 --no-zpl --order-es 160 --order-gs 100 --ene-range -500.0 10000.0 --resol 2101 --npts 100001 --fd-npts 40001 --fd-stencil 5 --fd-margin 2.0 --fd-fit poly2 --tol-fc 1e-07 --tol-vof 1e-08
python testing/compare_vof.py --reference hermite --test fd --dataset nv_mgo --temp 5.0 --sigma 20.0 --gamma 0.3 --no-zpl --order-es 160 --order-gs 100 --ene-range -500.0 10000.0 --resol 2101 --npts 100001 --fd-npts 40001 --fd-stencil 5 --fd-margin 2.0 --fd-fit poly4 --tol-fc 1e-07 --tol-vof 1e-08
python testing/compare_vof.py --reference hermite --test fd --dataset nv_mgo --temp 5.0 --sigma 20.0 --gamma 0.3 --no-zpl --order-es 160 --order-gs 100 --ene-range -500.0 10000.0 --resol 2101 --npts 100001 --fd-npts 40001 --fd-stencil 5 --fd-margin 2.0 --fd-fit spline4-notaknot --tol-fc 1e-07 --tol-vof 1e-08
python testing/compare_vof.py --reference hermite --test fd --dataset nv_mgo --temp 5.0 --sigma 20.0 --gamma 0.3 --no-zpl --order-es 160 --order-gs 100 --ene-range -500.0 10000.0 --resol 2101 --npts 100001 --fd-npts 40001 --fd-stencil 5 --fd-margin 2.0 --fd-fit spline3-s0.05 --tol-fc 1e-07 --tol-vof 1e-08
python testing/compare_vof.py --reference hermite --test cubic --dataset nv_mgo --temp 5.0 --sigma 20.0 --gamma 0.3 --no-zpl --order-es 160 --order-gs 100 --ene-range -500.0 10000.0 --resol 2101 --npts 100001 --fd-npts 40001 --fd-stencil 5 --fd-margin 2.0 --fd-fit parabola --cubic-points all --tol-fc 1e-07 --tol-vof 1e-08
python testing/compare_vof.py --reference hermite --test cubic --dataset nv_mgo --temp 5.0 --sigma 20.0 --gamma 0.3 --no-zpl --order-es 160 --order-gs 100 --ene-range -500.0 10000.0 --resol 2101 --npts 100001 --fd-npts 40001 --fd-stencil 5 --fd-margin 2.0 --fd-fit parabola --cubic-points all --cubic-free-a --tol-fc 1e-07 --tol-vof 1e-08
python testing/plot_ccd_wavefunctions.py --wf-method hermite --vof-method hermite pypl fd --dataset nv_mgo --temp 5.0 --sigma 20.0 --gamma 0.3 --no-zpl --order-es 160 --order-gs 100 --ene-range -500.0 10000.0 --resol 2101 --npts 100001 --fd-npts 40001 --fd-stencil 5 --fd-margin 2.0 --fd-fit parabola --n-levels 3
python testing/plot_ccd_wavefunctions.py --wf-method fd --vof-method hermite fd --dataset nv_mgo --temp 5.0 --sigma 20.0 --gamma 0.3 --no-zpl --order-es 160 --order-gs 100 --ene-range -500.0 10000.0 --resol 2101 --npts 100001 --fd-npts 40001 --fd-stencil 5 --fd-margin 2.0 --fd-fit parabola --n-levels 3
python testing/plot_ccd_wavefunctions.py --wf-method fd --vof-method hermite fd --dataset nv_mgo --temp 5.0 --sigma 20.0 --gamma 0.3 --no-zpl --order-es 160 --order-gs 100 --ene-range -500.0 10000.0 --resol 2101 --npts 100001 --fd-npts 40001 --fd-stencil 5 --fd-margin 2.0 --fd-fit spline --n-levels 3
python testing/plot_ccd_wavefunctions.py --wf-method fd --vof-method hermite fd --dataset nv_mgo --temp 5.0 --sigma 20.0 --gamma 0.3 --no-zpl --order-es 160 --order-gs 100 --ene-range -500.0 10000.0 --resol 2101 --npts 100001 --fd-npts 40001 --fd-stencil 5 --fd-margin 2.0 --fd-fit spline5-clamped --n-levels 3
python testing/plot_ccd_wavefunctions.py --wf-method fd --vof-method hermite fd --dataset nv_mgo --temp 5.0 --sigma 20.0 --gamma 0.3 --no-zpl --order-es 160 --order-gs 100 --ene-range -500.0 10000.0 --resol 2101 --npts 100001 --fd-npts 40001 --fd-stencil 5 --fd-margin 2.0 --fd-fit morse --n-levels 3
python testing/plot_ccd_wavefunctions.py --wf-method fd --vof-method hermite fd --dataset nv_mgo --temp 5.0 --sigma 20.0 --gamma 0.3 --no-zpl --order-es 160 --order-gs 100 --ene-range -500.0 10000.0 --resol 2101 --npts 100001 --fd-npts 40001 --fd-stencil 5 --fd-margin 2.0 --fd-fit poly2 --n-levels 3
python testing/plot_ccd_wavefunctions.py --wf-method fd --vof-method hermite fd --dataset nv_mgo --temp 5.0 --sigma 20.0 --gamma 0.3 --no-zpl --order-es 160 --order-gs 100 --ene-range -500.0 10000.0 --resol 2101 --npts 100001 --fd-npts 40001 --fd-stencil 5 --fd-margin 2.0 --fd-fit poly4 --n-levels 3
python testing/plot_ccd_wavefunctions.py --wf-method fd --vof-method hermite fd --dataset nv_mgo --temp 5.0 --sigma 20.0 --gamma 0.3 --no-zpl --order-es 160 --order-gs 100 --ene-range -500.0 10000.0 --resol 2101 --npts 100001 --fd-npts 40001 --fd-stencil 5 --fd-margin 2.0 --fd-fit spline4-notaknot --n-levels 3
python testing/plot_ccd_wavefunctions.py --wf-method fd --vof-method hermite fd --dataset nv_mgo --temp 5.0 --sigma 20.0 --gamma 0.3 --no-zpl --order-es 160 --order-gs 100 --ene-range -500.0 10000.0 --resol 2101 --npts 100001 --fd-npts 40001 --fd-stencil 5 --fd-margin 2.0 --fd-fit spline3-s0.05 --n-levels 3
python testing/plot_ccd_wavefunctions.py --wf-method cubic --vof-method hermite cubic --dataset nv_mgo --temp 5.0 --sigma 20.0 --gamma 0.3 --no-zpl --order-es 160 --order-gs 100 --ene-range -500.0 10000.0 --resol 2101 --npts 100001 --fd-npts 40001 --fd-stencil 5 --fd-margin 2.0 --fd-fit parabola --cubic-points all --n-levels 3
python testing/plot_ccd_wavefunctions.py --wf-method cubic --vof-method hermite cubic --dataset nv_mgo --temp 5.0 --sigma 20.0 --gamma 0.3 --no-zpl --order-es 160 --order-gs 100 --ene-range -500.0 10000.0 --resol 2101 --npts 100001 --fd-npts 40001 --fd-stencil 5 --fd-margin 2.0 --fd-fit parabola --cubic-points all --cubic-free-a --n-levels 3
python testing/vof_report.py --reference hermite --dataset nv_mgo --temp 5.0 --sigma 20.0 --gamma 0.3 --no-zpl --order-es 160 --order-gs 100 --ene-range -500.0 10000.0 --resol 2101 --npts 100001 --fd-npts 40001 --fd-stencil 5 --fd-margin 2.0 --fd-fit parabola --cubic-points all --tol-fc 1e-07 --tol-vof 1e-08
```

## Added after this run

Generated separately with the same settings, without rerunning the suite:

| added | content | output |
| --- | --- | --- |
| 2026-10-16 | fd wavefunctions and VOF on the cubic not-a-knot spline with smooth (C2) tails, `--fd-fit spline3-notaknot-smooth` (all images; tails continue the spline's own value, slope and curvature). The recommended spline for this dataset. | [`plot_ccd_wavefunctions/wf-fd_vof-hermite-fd_fit-spline3-notaknot-smooth/`](plot_ccd_wavefunctions/wf-fd_vof-hermite-fd_fit-spline3-notaknot-smooth/index.md) |
