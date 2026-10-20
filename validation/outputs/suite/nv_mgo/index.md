# PyPL 1D CCD test suite: nv_mgo (NV- in MgO)

Started 2026-10-17T23:49:29, 25 steps in 1807 s. Configuration: `run_all_nv_mgo.toml` (copied here), resolved values in `suite_config.json`.

**Dataset.**

- dataset: nv_mgo (NV- in MgO), 3E -> 1A1 (WEST TDDFT-DDH)
- data: doc/tutorials/003_nv_mgo_isc
- roles: 'ES' = 3E (upper state, minimum at Q = 0), 'GS' = 1A1 (lower state, minimum at Q = delta_q); emission = ES -> GS = 3E -> 1A1
- energies: 3E = spin-conserving TDDFT@DDH (WEST, root 1); 1A1 = spin-flip TDDFT@DDH (WEST, root 4); zero = the 3A2 minimum, DFT with the DDH functional (PBE0, exact exchange 0.34)
- Q = q delta_q, q from Image_<q>, as in reference/ccd_nv_mgo.py: 3E minimum at Q = 0, 1A1 minimum at Q = delta_q
- E = QE total energy (qe/) + WEST root (3E: sc root 1, 1A1: sf root 4) x 13.605693 eV/Ry, relative to the 3A2 (QE) minimum, as in ccd_nv_mgo.py
- image windows (q of Image_<q>, inclusive): parabola 3E q = -0.10..+0.45, 1A1 q = +0.60..+1.15; cubic a Q^2 + b Q^3 + c (perturbative correction) q = -0.10..+1.15 for both; every other fit (splines, Morse, polyN) ALL images
- all images: 3E all 33 sc images (q = -1.00..+2.00), 1A1 all 31 sf images (q = -1.00..+2.00)

**Physics:** T = 5 K, sigma = [20.0] meV, ZPL Gaussian, 160 ES x 100 GS levels, axis -500..10000 meV (2101 points). **Reference:** Hermite (quadrature). **fd:** 40001 points, 5-point stencil, margin 2. **Hermite:** 100001 points.

**Result: 2 step(s) failed or errored: harmonic methods like for like: hermite vs pypl, fd; full seven-section report.**

Skipped (no bound states on the fd grid):

- `poly6`: fit poly6 of the 1A1 images is unbounded below and has no bound states
- `poly8`: fit poly8 of the 1A1 images is unbounded below and has no bound states
- `cubic`: the 1A1 cubic turns over at Q = -1.70, inside the grid -6.68..10.23: no bound states in the box
- `cubic (exact fd check of cubic PT)`: the 1A1 cubic turns over at Q = -1.70, inside the grid -6.68..10.23: no bound states in the box
- `cubic-free (exact fd check of cubic PT)`: the 1A1 cubic turns over at Q = -0.50, inside the grid -6.68..10.23: no bound states in the box

## Method comparison

Largest difference over emission and absorption. **Like for like:** methods solving the same harmonic problem; they must agree within the tolerances (max dF_ij <= 1e-07, max dVOF / peak <= 1e-08). **PT vs exact cubic:** the error of first-order perturbation theory, against the exact fd solution of the same cubic. **Anharmonic vs harmonic:** how much each fit or correction changes the VOF from the harmonic reference (physics, not error), smallest first. See `validation/METHODS.md` for what each method does.

| kind | method | compared with | max dVOF / peak | max dF_ij | F_ij time (s) | verdict |
| --- | --- | --- | --- | --- | --- | --- |
| like for like | PyPL (Ruhoff) | Hermite (quadrature) | 9.62e-15 | 2.68e+27 | 0.11 | FAIL |
| like for like | Finite difference | Hermite (quadrature) | 1.95e-09 | 1.03e-08 | 23.55 | OK |
| anharmonic vs harmonic | Finite difference [Morse] | Hermite (quadrature) | 4.49e+00 | 5.79e-01 | 45.53 | differs (expected) |
| anharmonic vs harmonic | Hermite + cubic PT [cubic (a fixed, images q = -0.10..+1.15)] | Hermite (quadrature) | 4.79e+00 | 1.73e+02 | 4.30 | differs (expected) |
| anharmonic vs harmonic | Finite difference [polynomial deg. 4] | Hermite (quadrature) | 5.83e+00 | 6.55e-01 | 24.83 | differs (expected) |
| anharmonic vs harmonic | Finite difference [polynomial deg. 2] | Hermite (quadrature) | 7.10e+00 | 5.43e-01 | 92.35 | differs (expected) |
| anharmonic vs harmonic | Finite difference [quintic spline (clamped)] | Hermite (quadrature) | 7.60e+00 | 8.13e-01 | 16.90 | differs (expected) |
| anharmonic vs harmonic | Finite difference [quartic spline (not-a-knot)] | Hermite (quadrature) | 8.95e+00 | 7.85e-01 | 23.56 | differs (expected) |
| anharmonic vs harmonic | Finite difference [cubic spline (not-a-knot, smooth C2 tails)] | Hermite (quadrature) | 9.58e+00 | 7.95e-01 | 27.99 | differs (expected) |
| anharmonic vs harmonic | Finite difference [cubic smoothing spline (rms 0.05 meV)] | Hermite (quadrature) | 9.58e+00 | 7.76e-01 | 21.81 | differs (expected) |
| anharmonic vs harmonic | Hermite + cubic PT [cubic (a fitted, images q = -0.10..+1.15)] | Hermite (quadrature) | 9.77e+00 | 4.76e+02 | 5.02 | differs (expected) |
| anharmonic vs harmonic | Finite difference [cubic spline (natural)] | Hermite (quadrature) | 1.05e+01 | 7.63e-01 | 39.84 | differs (expected) |

## Steps

| # | step | status | time (s) | output |
| --- | --- | --- | --- | --- |
| 1 | checks: 81 regression and physics checks | PASS (81/81 checks passed) | 21.1 | `logs/01_checks_81-regression-and-physics-checks.log` |
| 2 | compare: harmonic methods like for like: hermite vs pypl, fd | FAIL | 61.9 | [`compare_vof/hermite_vs_pypl-fd/`](compare_vof/hermite_vs_pypl-fd/index.md) |
| 3 | compare: fd on cubic spline (natural) vs hermite | DIFFERS (expected) | 91.5 | [`compare_vof/hermite_vs_fd_fit-spline/`](compare_vof/hermite_vs_fd_fit-spline/index.md) |
| 4 | compare: fd on quintic spline (clamped) vs hermite | DIFFERS (expected) | 43.9 | [`compare_vof/hermite_vs_fd_fit-spline5-clamped/`](compare_vof/hermite_vs_fd_fit-spline5-clamped/index.md) |
| 5 | compare: fd on Morse vs hermite | DIFFERS (expected) | 108.2 | [`compare_vof/hermite_vs_fd_fit-morse/`](compare_vof/hermite_vs_fd_fit-morse/index.md) |
| 6 | compare: fd on polynomial deg. 2 vs hermite | DIFFERS (expected) | 220.2 | [`compare_vof/hermite_vs_fd_fit-poly2/`](compare_vof/hermite_vs_fd_fit-poly2/index.md) |
| 7 | compare: fd on polynomial deg. 4 vs hermite | DIFFERS (expected) | 59.6 | [`compare_vof/hermite_vs_fd_fit-poly4/`](compare_vof/hermite_vs_fd_fit-poly4/index.md) |
| 8 | compare: fd on cubic spline (not-a-knot, smooth C2 tails) vs hermite | DIFFERS (expected) | 68.2 | [`compare_vof/hermite_vs_fd_fit-spline3-notaknot-smooth/`](compare_vof/hermite_vs_fd_fit-spline3-notaknot-smooth/index.md) |
| 9 | compare: fd on quartic spline (not-a-knot) vs hermite | DIFFERS (expected) | 73.0 | [`compare_vof/hermite_vs_fd_fit-spline4-notaknot/`](compare_vof/hermite_vs_fd_fit-spline4-notaknot/index.md) |
| 10 | compare: fd on cubic smoothing spline (rms 0.05 meV) vs hermite | DIFFERS (expected) | 53.0 | [`compare_vof/hermite_vs_fd_fit-spline3-s0.05/`](compare_vof/hermite_vs_fd_fit-spline3-s0.05/index.md) |
| 11 | compare: cubic PT [cubic (a fixed, images q = -0.10..+1.15)] vs hermite | DIFFERS (expected) | 30.1 | [`compare_vof/hermite_vs_cubic_cubicfit-cubic/`](compare_vof/hermite_vs_cubic_cubicfit-cubic/index.md) |
| 12 | compare: cubic PT [cubic (a fitted, images q = -0.10..+1.15)] vs hermite | DIFFERS (expected) | 27.7 | [`compare_vof/hermite_vs_cubic_cubicfit-cubic-free/`](compare_vof/hermite_vs_cubic_cubicfit-cubic-free/index.md) |
| 13 | plots: Hermite wavefunctions, VOF of hermite, pypl, fd | OK | 33.8 | [`plot_ccd_wavefunctions/wf-hermite_vof-hermite-pypl-fd_fit-parabola/`](plot_ccd_wavefunctions/wf-hermite_vof-hermite-pypl-fd_fit-parabola/index.md) |
| 14 | plots: fd wavefunctions on parabola (near min.) | OK | 27.1 | [`plot_ccd_wavefunctions/wf-fd_vof-hermite-fd_fit-parabola/`](plot_ccd_wavefunctions/wf-fd_vof-hermite-fd_fit-parabola/index.md) |
| 15 | plots: fd wavefunctions on cubic spline (natural) | OK | 31.0 | [`plot_ccd_wavefunctions/wf-fd_vof-hermite-fd_fit-spline/`](plot_ccd_wavefunctions/wf-fd_vof-hermite-fd_fit-spline/index.md) |
| 16 | plots: fd wavefunctions on quintic spline (clamped) | OK | 30.3 | [`plot_ccd_wavefunctions/wf-fd_vof-hermite-fd_fit-spline5-clamped/`](plot_ccd_wavefunctions/wf-fd_vof-hermite-fd_fit-spline5-clamped/index.md) |
| 17 | plots: fd wavefunctions on Morse | OK | 30.7 | [`plot_ccd_wavefunctions/wf-fd_vof-hermite-fd_fit-morse/`](plot_ccd_wavefunctions/wf-fd_vof-hermite-fd_fit-morse/index.md) |
| 18 | plots: fd wavefunctions on polynomial deg. 2 | OK | 30.8 | [`plot_ccd_wavefunctions/wf-fd_vof-hermite-fd_fit-poly2/`](plot_ccd_wavefunctions/wf-fd_vof-hermite-fd_fit-poly2/index.md) |
| 19 | plots: fd wavefunctions on polynomial deg. 4 | OK | 29.4 | [`plot_ccd_wavefunctions/wf-fd_vof-hermite-fd_fit-poly4/`](plot_ccd_wavefunctions/wf-fd_vof-hermite-fd_fit-poly4/index.md) |
| 20 | plots: fd wavefunctions on cubic spline (not-a-knot, smooth C2 tails) | OK | 26.6 | [`plot_ccd_wavefunctions/wf-fd_vof-hermite-fd_fit-spline3-notaknot-smooth/`](plot_ccd_wavefunctions/wf-fd_vof-hermite-fd_fit-spline3-notaknot-smooth/index.md) |
| 21 | plots: fd wavefunctions on quartic spline (not-a-knot) | OK | 30.6 | [`plot_ccd_wavefunctions/wf-fd_vof-hermite-fd_fit-spline4-notaknot/`](plot_ccd_wavefunctions/wf-fd_vof-hermite-fd_fit-spline4-notaknot/index.md) |
| 22 | plots: fd wavefunctions on cubic smoothing spline (rms 0.05 meV) | OK | 31.3 | [`plot_ccd_wavefunctions/wf-fd_vof-hermite-fd_fit-spline3-s0.05/`](plot_ccd_wavefunctions/wf-fd_vof-hermite-fd_fit-spline3-s0.05/index.md) |
| 23 | plots: cubic PT wavefunctions [cubic (a fixed, images q = -0.10..+1.15)] | OK | 11.0 | [`plot_ccd_wavefunctions/wf-cubic_vof-hermite-cubic_cubicfit-cubic/`](plot_ccd_wavefunctions/wf-cubic_vof-hermite-cubic_cubicfit-cubic/index.md) |
| 24 | plots: cubic PT wavefunctions [cubic (a fitted, images q = -0.10..+1.15)] | OK | 11.3 | [`plot_ccd_wavefunctions/wf-cubic_vof-hermite-cubic_cubicfit-cubic-free/`](plot_ccd_wavefunctions/wf-cubic_vof-hermite-cubic_cubicfit-cubic-free/index.md) |
| 25 | report: full seven-section report | FAIL | 624.6 | [`vof_report/default/`](vof_report/default/index.md) |

## Commands

Each step was run as (from the repository root, with `PYPL_VALIDATION_OUTPUT_ROOT` pointing at this directory and `PYPL_VALIDATION_DATASET=nv_mgo`):

```bash
python validation/check_1d_ccd.py -v
python validation/compare_vof.py --reference hermite --test pypl fd --dataset nv_mgo --temp 5.0 --sigma 20.0 --gamma 0.3 --no-zpl --order-es 160 --order-gs 100 --ene-range -500.0 10000.0 --resol 2101 --line-positions harmonic --npts 100001 --hermite-margin 6.0 --fd-stencil 5 --fd-margin 2.0 --fd-fit parabola --fd-npts 40001 --cubic-warn 0.1 --dpi 300 --fig-format png --font-size 9.0 --tol-fc 1e-07 --tol-vof 1e-08
python validation/compare_vof.py --reference hermite --test fd --dataset nv_mgo --temp 5.0 --sigma 20.0 --gamma 0.3 --no-zpl --order-es 160 --order-gs 100 --ene-range -500.0 10000.0 --resol 2101 --line-positions harmonic --npts 100001 --hermite-margin 6.0 --fd-stencil 5 --fd-margin 2.0 --fd-fit spline --fd-npts 40001 --cubic-warn 0.1 --dpi 300 --fig-format png --font-size 9.0 --tol-fc 1e-07 --tol-vof 1e-08
python validation/compare_vof.py --reference hermite --test fd --dataset nv_mgo --temp 5.0 --sigma 20.0 --gamma 0.3 --no-zpl --order-es 160 --order-gs 100 --ene-range -500.0 10000.0 --resol 2101 --line-positions harmonic --npts 100001 --hermite-margin 6.0 --fd-stencil 5 --fd-margin 2.0 --fd-fit spline5-clamped --fd-npts 40001 --cubic-warn 0.1 --dpi 300 --fig-format png --font-size 9.0 --tol-fc 1e-07 --tol-vof 1e-08
python validation/compare_vof.py --reference hermite --test fd --dataset nv_mgo --temp 5.0 --sigma 20.0 --gamma 0.3 --no-zpl --order-es 160 --order-gs 100 --ene-range -500.0 10000.0 --resol 2101 --line-positions harmonic --npts 100001 --hermite-margin 6.0 --fd-stencil 5 --fd-margin 2.0 --fd-fit morse --fd-npts 40001 --cubic-warn 0.1 --dpi 300 --fig-format png --font-size 9.0 --tol-fc 1e-07 --tol-vof 1e-08
python validation/compare_vof.py --reference hermite --test fd --dataset nv_mgo --temp 5.0 --sigma 20.0 --gamma 0.3 --no-zpl --order-es 160 --order-gs 100 --ene-range -500.0 10000.0 --resol 2101 --line-positions harmonic --npts 100001 --hermite-margin 6.0 --fd-stencil 5 --fd-margin 2.0 --fd-fit poly2 --fd-npts 40001 --cubic-warn 0.1 --dpi 300 --fig-format png --font-size 9.0 --tol-fc 1e-07 --tol-vof 1e-08
python validation/compare_vof.py --reference hermite --test fd --dataset nv_mgo --temp 5.0 --sigma 20.0 --gamma 0.3 --no-zpl --order-es 160 --order-gs 100 --ene-range -500.0 10000.0 --resol 2101 --line-positions harmonic --npts 100001 --hermite-margin 6.0 --fd-stencil 5 --fd-margin 2.0 --fd-fit poly4 --fd-npts 40001 --cubic-warn 0.1 --dpi 300 --fig-format png --font-size 9.0 --tol-fc 1e-07 --tol-vof 1e-08
python validation/compare_vof.py --reference hermite --test fd --dataset nv_mgo --temp 5.0 --sigma 20.0 --gamma 0.3 --no-zpl --order-es 160 --order-gs 100 --ene-range -500.0 10000.0 --resol 2101 --line-positions harmonic --npts 100001 --hermite-margin 6.0 --fd-stencil 5 --fd-margin 2.0 --fd-fit spline3-notaknot-smooth --fd-npts 40001 --cubic-warn 0.1 --dpi 300 --fig-format png --font-size 9.0 --tol-fc 1e-07 --tol-vof 1e-08
python validation/compare_vof.py --reference hermite --test fd --dataset nv_mgo --temp 5.0 --sigma 20.0 --gamma 0.3 --no-zpl --order-es 160 --order-gs 100 --ene-range -500.0 10000.0 --resol 2101 --line-positions harmonic --npts 100001 --hermite-margin 6.0 --fd-stencil 5 --fd-margin 2.0 --fd-fit spline4-notaknot --fd-npts 40001 --cubic-warn 0.1 --dpi 300 --fig-format png --font-size 9.0 --tol-fc 1e-07 --tol-vof 1e-08
python validation/compare_vof.py --reference hermite --test fd --dataset nv_mgo --temp 5.0 --sigma 20.0 --gamma 0.3 --no-zpl --order-es 160 --order-gs 100 --ene-range -500.0 10000.0 --resol 2101 --line-positions harmonic --npts 100001 --hermite-margin 6.0 --fd-stencil 5 --fd-margin 2.0 --fd-fit spline3-s0.05 --fd-npts 40001 --cubic-warn 0.1 --dpi 300 --fig-format png --font-size 9.0 --tol-fc 1e-07 --tol-vof 1e-08
python validation/compare_vof.py --reference hermite --test cubic --dataset nv_mgo --temp 5.0 --sigma 20.0 --gamma 0.3 --no-zpl --order-es 160 --order-gs 100 --ene-range -500.0 10000.0 --resol 2101 --line-positions harmonic --npts 100001 --hermite-margin 6.0 --fd-stencil 5 --fd-margin 2.0 --fd-fit parabola --fd-npts 40001 --cubic-warn 0.1 --dpi 300 --fig-format png --font-size 9.0 --cubic-points all --tol-fc 1e-07 --tol-vof 1e-08
python validation/compare_vof.py --reference hermite --test cubic --dataset nv_mgo --temp 5.0 --sigma 20.0 --gamma 0.3 --no-zpl --order-es 160 --order-gs 100 --ene-range -500.0 10000.0 --resol 2101 --line-positions harmonic --npts 100001 --hermite-margin 6.0 --fd-stencil 5 --fd-margin 2.0 --fd-fit parabola --fd-npts 40001 --cubic-warn 0.1 --dpi 300 --fig-format png --font-size 9.0 --cubic-points all --cubic-free-a --tol-fc 1e-07 --tol-vof 1e-08
python validation/plot_ccd_wavefunctions.py --wf-method hermite --vof-method hermite pypl fd --dataset nv_mgo --temp 5.0 --sigma 20.0 --gamma 0.3 --no-zpl --order-es 160 --order-gs 100 --ene-range -500.0 10000.0 --resol 2101 --line-positions harmonic --npts 100001 --hermite-margin 6.0 --fd-stencil 5 --fd-margin 2.0 --fd-fit parabola --fd-npts 40001 --cubic-warn 0.1 --dpi 300 --fig-format png --font-size 9.0 --n-levels 3 --check-tol 1e-10 1e-08
python validation/plot_ccd_wavefunctions.py --wf-method fd --vof-method hermite fd --dataset nv_mgo --temp 5.0 --sigma 20.0 --gamma 0.3 --no-zpl --order-es 160 --order-gs 100 --ene-range -500.0 10000.0 --resol 2101 --line-positions harmonic --npts 100001 --hermite-margin 6.0 --fd-stencil 5 --fd-margin 2.0 --fd-fit parabola --fd-npts 40001 --cubic-warn 0.1 --dpi 300 --fig-format png --font-size 9.0 --n-levels 3 --check-tol 1e-10 1e-08
python validation/plot_ccd_wavefunctions.py --wf-method fd --vof-method hermite fd --dataset nv_mgo --temp 5.0 --sigma 20.0 --gamma 0.3 --no-zpl --order-es 160 --order-gs 100 --ene-range -500.0 10000.0 --resol 2101 --line-positions harmonic --npts 100001 --hermite-margin 6.0 --fd-stencil 5 --fd-margin 2.0 --fd-fit spline --fd-npts 40001 --cubic-warn 0.1 --dpi 300 --fig-format png --font-size 9.0 --n-levels 3 --check-tol 1e-10 1e-08
python validation/plot_ccd_wavefunctions.py --wf-method fd --vof-method hermite fd --dataset nv_mgo --temp 5.0 --sigma 20.0 --gamma 0.3 --no-zpl --order-es 160 --order-gs 100 --ene-range -500.0 10000.0 --resol 2101 --line-positions harmonic --npts 100001 --hermite-margin 6.0 --fd-stencil 5 --fd-margin 2.0 --fd-fit spline5-clamped --fd-npts 40001 --cubic-warn 0.1 --dpi 300 --fig-format png --font-size 9.0 --n-levels 3 --check-tol 1e-10 1e-08
python validation/plot_ccd_wavefunctions.py --wf-method fd --vof-method hermite fd --dataset nv_mgo --temp 5.0 --sigma 20.0 --gamma 0.3 --no-zpl --order-es 160 --order-gs 100 --ene-range -500.0 10000.0 --resol 2101 --line-positions harmonic --npts 100001 --hermite-margin 6.0 --fd-stencil 5 --fd-margin 2.0 --fd-fit morse --fd-npts 40001 --cubic-warn 0.1 --dpi 300 --fig-format png --font-size 9.0 --n-levels 3 --check-tol 1e-10 1e-08
python validation/plot_ccd_wavefunctions.py --wf-method fd --vof-method hermite fd --dataset nv_mgo --temp 5.0 --sigma 20.0 --gamma 0.3 --no-zpl --order-es 160 --order-gs 100 --ene-range -500.0 10000.0 --resol 2101 --line-positions harmonic --npts 100001 --hermite-margin 6.0 --fd-stencil 5 --fd-margin 2.0 --fd-fit poly2 --fd-npts 40001 --cubic-warn 0.1 --dpi 300 --fig-format png --font-size 9.0 --n-levels 3 --check-tol 1e-10 1e-08
python validation/plot_ccd_wavefunctions.py --wf-method fd --vof-method hermite fd --dataset nv_mgo --temp 5.0 --sigma 20.0 --gamma 0.3 --no-zpl --order-es 160 --order-gs 100 --ene-range -500.0 10000.0 --resol 2101 --line-positions harmonic --npts 100001 --hermite-margin 6.0 --fd-stencil 5 --fd-margin 2.0 --fd-fit poly4 --fd-npts 40001 --cubic-warn 0.1 --dpi 300 --fig-format png --font-size 9.0 --n-levels 3 --check-tol 1e-10 1e-08
python validation/plot_ccd_wavefunctions.py --wf-method fd --vof-method hermite fd --dataset nv_mgo --temp 5.0 --sigma 20.0 --gamma 0.3 --no-zpl --order-es 160 --order-gs 100 --ene-range -500.0 10000.0 --resol 2101 --line-positions harmonic --npts 100001 --hermite-margin 6.0 --fd-stencil 5 --fd-margin 2.0 --fd-fit spline3-notaknot-smooth --fd-npts 40001 --cubic-warn 0.1 --dpi 300 --fig-format png --font-size 9.0 --n-levels 3 --check-tol 1e-10 1e-08
python validation/plot_ccd_wavefunctions.py --wf-method fd --vof-method hermite fd --dataset nv_mgo --temp 5.0 --sigma 20.0 --gamma 0.3 --no-zpl --order-es 160 --order-gs 100 --ene-range -500.0 10000.0 --resol 2101 --line-positions harmonic --npts 100001 --hermite-margin 6.0 --fd-stencil 5 --fd-margin 2.0 --fd-fit spline4-notaknot --fd-npts 40001 --cubic-warn 0.1 --dpi 300 --fig-format png --font-size 9.0 --n-levels 3 --check-tol 1e-10 1e-08
python validation/plot_ccd_wavefunctions.py --wf-method fd --vof-method hermite fd --dataset nv_mgo --temp 5.0 --sigma 20.0 --gamma 0.3 --no-zpl --order-es 160 --order-gs 100 --ene-range -500.0 10000.0 --resol 2101 --line-positions harmonic --npts 100001 --hermite-margin 6.0 --fd-stencil 5 --fd-margin 2.0 --fd-fit spline3-s0.05 --fd-npts 40001 --cubic-warn 0.1 --dpi 300 --fig-format png --font-size 9.0 --n-levels 3 --check-tol 1e-10 1e-08
python validation/plot_ccd_wavefunctions.py --wf-method cubic --vof-method hermite cubic --dataset nv_mgo --temp 5.0 --sigma 20.0 --gamma 0.3 --no-zpl --order-es 160 --order-gs 100 --ene-range -500.0 10000.0 --resol 2101 --line-positions harmonic --npts 100001 --hermite-margin 6.0 --fd-stencil 5 --fd-margin 2.0 --fd-fit parabola --fd-npts 40001 --cubic-warn 0.1 --dpi 300 --fig-format png --font-size 9.0 --cubic-points all --n-levels 3 --check-tol 1e-10 1e-08
python validation/plot_ccd_wavefunctions.py --wf-method cubic --vof-method hermite cubic --dataset nv_mgo --temp 5.0 --sigma 20.0 --gamma 0.3 --no-zpl --order-es 160 --order-gs 100 --ene-range -500.0 10000.0 --resol 2101 --line-positions harmonic --npts 100001 --hermite-margin 6.0 --fd-stencil 5 --fd-margin 2.0 --fd-fit parabola --fd-npts 40001 --cubic-warn 0.1 --dpi 300 --fig-format png --font-size 9.0 --cubic-points all --cubic-free-a --n-levels 3 --check-tol 1e-10 1e-08
python validation/vof_report.py --reference hermite --dataset nv_mgo --temp 5.0 --sigma 20.0 --gamma 0.3 --no-zpl --order-es 160 --order-gs 100 --ene-range -500.0 10000.0 --resol 2101 --line-positions harmonic --npts 100001 --hermite-margin 6.0 --fd-stencil 5 --fd-margin 2.0 --fd-fit parabola --fd-npts 40001 --cubic-warn 0.1 --dpi 300 --fig-format png --font-size 9.0 --cubic-points all --tol-fc 1e-07 --tol-vof 1e-08
```
