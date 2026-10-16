# PyPL 1D CCD test suite: default

Started 2026-10-15T18:27:42, 28 steps in 1386 s. Configuration: `run_all.toml` (copied here), resolved values in `suite_config.json`.

**Physics:** T = 5 K, sigma = [8.0, 25.0] meV, ZPL Lorentzian, gamma = 0.3 meV, 50 ES x 60 GS levels, axis -200..1200 meV (1401 points). **Reference:** PyPL (Ruhoff). **fd:** 20001 points, 5-point stencil, margin 2. **Hermite:** 100001 points.

**Result: all like-for-like steps passed, nothing errored.**

## Method comparison

Largest difference over emission and absorption. **Like for like:** methods solving the same harmonic problem; they must agree within the tolerances (max dF_ij <= 1e-07, max dVOF / peak <= 1e-08). **PT vs exact cubic:** the error of first-order perturbation theory, against the exact fd solution of the same cubic. **Anharmonic vs harmonic:** how much each fit or correction changes the VOF from the harmonic reference (physics, not error), smallest first. See `testing/METHODS.md` for what each method does.

| kind | method | compared with | max dVOF / peak | max dF_ij | F_ij time (s) | verdict |
| --- | --- | --- | --- | --- | --- | --- |
| like for like | Hermite (quadrature) | PyPL (Ruhoff) | 6.95e-16 | 1.70e-08 | 1.26 | OK |
| like for like | Finite difference | PyPL (Ruhoff) | 1.29e-12 | 1.80e-08 | 2.51 | OK |
| PT vs exact cubic | Hermite + cubic PT [cubic (a fitted, all images)] | Finite difference [cubic (a fitted, all images)] | 2.71e-03 | 9.60e-01 | 4.98 | differs (expected) |
| PT vs exact cubic | Hermite + cubic PT [cubic (a fixed, all images)] | Finite difference [cubic (a fixed, all images)] | 3.25e-04 | 1.16e+00 | 4.56 | differs (expected) |
| anharmonic vs harmonic | Finite difference [polynomial deg. 2] | PyPL (Ruhoff) | 2.42e-03 | 3.36e-03 | 3.39 | differs (expected) |
| anharmonic vs harmonic | Finite difference [Morse] | PyPL (Ruhoff) | 1.95e-02 | 4.01e-01 | 3.14 | differs (expected) |
| anharmonic vs harmonic | Finite difference [cubic spline (natural)] | PyPL (Ruhoff) | 2.08e-02 | 1.79e-02 | 2.90 | differs (expected) |
| anharmonic vs harmonic | Finite difference [quintic spline (clamped)] | PyPL (Ruhoff) | 2.16e-02 | 1.79e-02 | 4.20 | differs (expected) |
| anharmonic vs harmonic | Finite difference [quartic spline (not-a-knot)] | PyPL (Ruhoff) | 2.18e-02 | 1.80e-02 | 28.74 | differs (expected) |
| anharmonic vs harmonic | Finite difference [cubic smoothing spline (rms 0.05 meV)] | PyPL (Ruhoff) | 2.21e-02 | 1.81e-02 | 40.09 | differs (expected) |
| anharmonic vs harmonic | Hermite + cubic PT [cubic (a fitted, all images)] | PyPL (Ruhoff) | 3.28e-02 | 9.81e-01 | 4.30 | differs (expected) |
| anharmonic vs harmonic | Finite difference [cubic (a fixed, all images)] | PyPL (Ruhoff) | 3.38e-02 | 4.89e-01 | 3.17 | differs (expected) |
| anharmonic vs harmonic | Hermite + cubic PT [cubic (a fixed, all images)] | PyPL (Ruhoff) | 3.41e-02 | 1.18e+00 | 4.65 | differs (expected) |
| anharmonic vs harmonic | Finite difference [polynomial deg. 6] | PyPL (Ruhoff) | 6.96e-02 | 7.06e-01 | 1.97 | differs (expected) |

## Steps

| # | step | status | time (s) | output |
| --- | --- | --- | --- | --- |
| 1 | checks: 59 regression and physics checks | PASS (59/59 checks passed) | 8.1 | `logs/01_checks_59-regression-and-physics-checks.log` |
| 2 | tutorial: tutorial 002 spectra (PyPL's own solver) | OK | 2.7 | [`run_1d_ccd/default/`](run_1d_ccd/default/index.md) |
| 3 | compare: harmonic methods like for like: pypl vs hermite, fd | PASS | 11.8 | [`compare_vof/pypl_vs_hermite-fd/`](compare_vof/pypl_vs_hermite-fd/index.md) |
| 4 | compare: fd on cubic spline (natural) vs pypl | DIFFERS (expected) | 9.8 | [`compare_vof/pypl_vs_fd_fit-spline/`](compare_vof/pypl_vs_fd_fit-spline/index.md) |
| 5 | compare: fd on quintic spline (clamped) vs pypl | DIFFERS (expected) | 12.6 | [`compare_vof/pypl_vs_fd_fit-spline5-clamped/`](compare_vof/pypl_vs_fd_fit-spline5-clamped/index.md) |
| 6 | compare: fd on Morse vs pypl | DIFFERS (expected) | 10.5 | [`compare_vof/pypl_vs_fd_fit-morse/`](compare_vof/pypl_vs_fd_fit-morse/index.md) |
| 7 | compare: fd on polynomial deg. 2 vs pypl | DIFFERS (expected) | 10.7 | [`compare_vof/pypl_vs_fd_fit-poly2/`](compare_vof/pypl_vs_fd_fit-poly2/index.md) |
| 8 | compare: fd on polynomial deg. 6 vs pypl | DIFFERS (expected) | 7.7 | [`compare_vof/pypl_vs_fd_fit-poly6/`](compare_vof/pypl_vs_fd_fit-poly6/index.md) |
| 9 | compare: fd on cubic (a fixed, all images) vs pypl | DIFFERS (expected) | 10.5 | [`compare_vof/pypl_vs_fd_fit-cubic/`](compare_vof/pypl_vs_fd_fit-cubic/index.md) |
| 10 | compare: fd on quartic spline (not-a-knot) vs pypl | DIFFERS (expected) | 67.7 | [`compare_vof/pypl_vs_fd_fit-spline4-notaknot/`](compare_vof/pypl_vs_fd_fit-spline4-notaknot/index.md) |
| 11 | compare: fd on cubic smoothing spline (rms 0.05 meV) vs pypl | DIFFERS (expected) | 90.4 | [`compare_vof/pypl_vs_fd_fit-spline3-s0.05/`](compare_vof/pypl_vs_fd_fit-spline3-s0.05/index.md) |
| 12 | compare: cubic PT [cubic (a fixed, all images)] vs pypl | DIFFERS (expected) | 24.2 | [`compare_vof/pypl_vs_cubic_cubicfit-cubic/`](compare_vof/pypl_vs_cubic_cubicfit-cubic/index.md) |
| 13 | compare: cubic PT [cubic (a fixed, all images)] vs the exact fd solution of that cubic | DIFFERS (expected) | 74.2 | [`compare_vof/fd_vs_cubic_fit-cubic_cubicfit-cubic/`](compare_vof/fd_vs_cubic_fit-cubic_cubicfit-cubic/index.md) |
| 14 | compare: cubic PT [cubic (a fitted, all images)] vs pypl | DIFFERS (expected) | 19.2 | [`compare_vof/pypl_vs_cubic_cubicfit-cubic-free/`](compare_vof/pypl_vs_cubic_cubicfit-cubic-free/index.md) |
| 15 | compare: cubic PT [cubic (a fitted, all images)] vs the exact fd solution of that cubic | DIFFERS (expected) | 90.3 | [`compare_vof/fd_vs_cubic_fit-cubic-free_cubicfit-cubic-free/`](compare_vof/fd_vs_cubic_fit-cubic-free_cubicfit-cubic-free/index.md) |
| 16 | plots: Hermite wavefunctions, VOF of pypl, hermite, fd | OK | 32.4 | [`plot_ccd_wavefunctions/wf-hermite_vof-pypl-hermite-fd_fit-parabola/`](plot_ccd_wavefunctions/wf-hermite_vof-pypl-hermite-fd_fit-parabola/index.md) |
| 17 | plots: fd wavefunctions on parabola (near min.) | OK | 8.3 | [`plot_ccd_wavefunctions/wf-fd_vof-pypl-fd_fit-parabola/`](plot_ccd_wavefunctions/wf-fd_vof-pypl-fd_fit-parabola/index.md) |
| 18 | plots: fd wavefunctions on cubic spline (natural) | OK | 8.8 | [`plot_ccd_wavefunctions/wf-fd_vof-pypl-fd_fit-spline/`](plot_ccd_wavefunctions/wf-fd_vof-pypl-fd_fit-spline/index.md) |
| 19 | plots: fd wavefunctions on quintic spline (clamped) | OK | 8.9 | [`plot_ccd_wavefunctions/wf-fd_vof-pypl-fd_fit-spline5-clamped/`](plot_ccd_wavefunctions/wf-fd_vof-pypl-fd_fit-spline5-clamped/index.md) |
| 20 | plots: fd wavefunctions on Morse | OK | 11.6 | [`plot_ccd_wavefunctions/wf-fd_vof-pypl-fd_fit-morse/`](plot_ccd_wavefunctions/wf-fd_vof-pypl-fd_fit-morse/index.md) |
| 21 | plots: fd wavefunctions on polynomial deg. 2 | OK | 68.7 | [`plot_ccd_wavefunctions/wf-fd_vof-pypl-fd_fit-poly2/`](plot_ccd_wavefunctions/wf-fd_vof-pypl-fd_fit-poly2/index.md) |
| 22 | plots: fd wavefunctions on polynomial deg. 6 | OK | 10.3 | [`plot_ccd_wavefunctions/wf-fd_vof-pypl-fd_fit-poly6/`](plot_ccd_wavefunctions/wf-fd_vof-pypl-fd_fit-poly6/index.md) |
| 23 | plots: fd wavefunctions on cubic (a fixed, all images) | OK | 60.8 | [`plot_ccd_wavefunctions/wf-fd_vof-pypl-fd_fit-cubic/`](plot_ccd_wavefunctions/wf-fd_vof-pypl-fd_fit-cubic/index.md) |
| 24 | plots: fd wavefunctions on quartic spline (not-a-knot) | OK | 53.0 | [`plot_ccd_wavefunctions/wf-fd_vof-pypl-fd_fit-spline4-notaknot/`](plot_ccd_wavefunctions/wf-fd_vof-pypl-fd_fit-spline4-notaknot/index.md) |
| 25 | plots: fd wavefunctions on cubic smoothing spline (rms 0.05 meV) | OK | 8.7 | [`plot_ccd_wavefunctions/wf-fd_vof-pypl-fd_fit-spline3-s0.05/`](plot_ccd_wavefunctions/wf-fd_vof-pypl-fd_fit-spline3-s0.05/index.md) |
| 26 | plots: cubic PT wavefunctions [cubic (a fixed, all images)] | OK | 17.8 | [`plot_ccd_wavefunctions/wf-cubic_vof-pypl-cubic_cubicfit-cubic/`](plot_ccd_wavefunctions/wf-cubic_vof-pypl-cubic_cubicfit-cubic/index.md) |
| 27 | plots: cubic PT wavefunctions [cubic (a fitted, all images)] | OK | 11.1 | [`plot_ccd_wavefunctions/wf-cubic_vof-pypl-cubic_cubicfit-cubic-free/`](plot_ccd_wavefunctions/wf-cubic_vof-pypl-cubic_cubicfit-cubic-free/index.md) |
| 28 | report: full seven-section report | PASS | 635.6 | [`vof_report/default/`](vof_report/default/index.md) |

## Commands

Each step was run as (from the repository root, with `PYPL_TESTING_OUTPUT_ROOT` pointing at this directory):

```bash
python testing/check_1d_ccd.py -v
python testing/run_1d_ccd.py --temp 5.0 --ezpl 1945.0 --order-es 50 --order-gs 60 --sigma 8.0 25.0 --gamma 0.3 --ene-range -200.0 1200.0 --resol 1401 --tdm 1.0
python testing/compare_vof.py --reference pypl --test hermite fd --temp 5.0 --sigma 8.0 25.0 --gamma 0.3 --order-es 50 --order-gs 60 --ene-range -200.0 1200.0 --resol 1401 --npts 100001 --fd-npts 20001 --fd-stencil 5 --fd-margin 2.0 --fd-fit parabola --tol-fc 1e-07 --tol-vof 1e-08
python testing/compare_vof.py --reference pypl --test fd --temp 5.0 --sigma 8.0 25.0 --gamma 0.3 --order-es 50 --order-gs 60 --ene-range -200.0 1200.0 --resol 1401 --npts 100001 --fd-npts 20001 --fd-stencil 5 --fd-margin 2.0 --fd-fit spline --tol-fc 1e-07 --tol-vof 1e-08
python testing/compare_vof.py --reference pypl --test fd --temp 5.0 --sigma 8.0 25.0 --gamma 0.3 --order-es 50 --order-gs 60 --ene-range -200.0 1200.0 --resol 1401 --npts 100001 --fd-npts 20001 --fd-stencil 5 --fd-margin 2.0 --fd-fit spline5-clamped --tol-fc 1e-07 --tol-vof 1e-08
python testing/compare_vof.py --reference pypl --test fd --temp 5.0 --sigma 8.0 25.0 --gamma 0.3 --order-es 50 --order-gs 60 --ene-range -200.0 1200.0 --resol 1401 --npts 100001 --fd-npts 20001 --fd-stencil 5 --fd-margin 2.0 --fd-fit morse --tol-fc 1e-07 --tol-vof 1e-08
python testing/compare_vof.py --reference pypl --test fd --temp 5.0 --sigma 8.0 25.0 --gamma 0.3 --order-es 50 --order-gs 60 --ene-range -200.0 1200.0 --resol 1401 --npts 100001 --fd-npts 20001 --fd-stencil 5 --fd-margin 2.0 --fd-fit poly2 --tol-fc 1e-07 --tol-vof 1e-08
python testing/compare_vof.py --reference pypl --test fd --temp 5.0 --sigma 8.0 25.0 --gamma 0.3 --order-es 50 --order-gs 60 --ene-range -200.0 1200.0 --resol 1401 --npts 100001 --fd-npts 20001 --fd-stencil 5 --fd-margin 2.0 --fd-fit poly6 --tol-fc 1e-07 --tol-vof 1e-08
python testing/compare_vof.py --reference pypl --test fd --temp 5.0 --sigma 8.0 25.0 --gamma 0.3 --order-es 50 --order-gs 60 --ene-range -200.0 1200.0 --resol 1401 --npts 100001 --fd-npts 20001 --fd-stencil 5 --fd-margin 2.0 --fd-fit cubic --tol-fc 1e-07 --tol-vof 1e-08
python testing/compare_vof.py --reference pypl --test fd --temp 5.0 --sigma 8.0 25.0 --gamma 0.3 --order-es 50 --order-gs 60 --ene-range -200.0 1200.0 --resol 1401 --npts 100001 --fd-npts 20001 --fd-stencil 5 --fd-margin 2.0 --fd-fit spline4-notaknot --tol-fc 1e-07 --tol-vof 1e-08
python testing/compare_vof.py --reference pypl --test fd --temp 5.0 --sigma 8.0 25.0 --gamma 0.3 --order-es 50 --order-gs 60 --ene-range -200.0 1200.0 --resol 1401 --npts 100001 --fd-npts 20001 --fd-stencil 5 --fd-margin 2.0 --fd-fit spline3-s0.05 --tol-fc 1e-07 --tol-vof 1e-08
python testing/compare_vof.py --reference pypl --test cubic --temp 5.0 --sigma 8.0 25.0 --gamma 0.3 --order-es 50 --order-gs 60 --ene-range -200.0 1200.0 --resol 1401 --npts 100001 --fd-npts 20001 --fd-stencil 5 --fd-margin 2.0 --fd-fit parabola --cubic-points all --tol-fc 1e-07 --tol-vof 1e-08
python testing/compare_vof.py --reference fd --test cubic --temp 5.0 --sigma 8.0 25.0 --gamma 0.3 --order-es 50 --order-gs 60 --ene-range -200.0 1200.0 --resol 1401 --npts 100001 --fd-npts 20001 --fd-stencil 5 --fd-margin 2.0 --fd-fit cubic --cubic-points all --tol-fc 1e-07 --tol-vof 1e-08
python testing/compare_vof.py --reference pypl --test cubic --temp 5.0 --sigma 8.0 25.0 --gamma 0.3 --order-es 50 --order-gs 60 --ene-range -200.0 1200.0 --resol 1401 --npts 100001 --fd-npts 20001 --fd-stencil 5 --fd-margin 2.0 --fd-fit parabola --cubic-points all --cubic-free-a --tol-fc 1e-07 --tol-vof 1e-08
python testing/compare_vof.py --reference fd --test cubic --temp 5.0 --sigma 8.0 25.0 --gamma 0.3 --order-es 50 --order-gs 60 --ene-range -200.0 1200.0 --resol 1401 --npts 100001 --fd-npts 20001 --fd-stencil 5 --fd-margin 2.0 --fd-fit cubic-free --cubic-points all --cubic-free-a --tol-fc 1e-07 --tol-vof 1e-08
python testing/plot_ccd_wavefunctions.py --wf-method hermite --vof-method pypl hermite fd --temp 5.0 --sigma 8.0 25.0 --gamma 0.3 --order-es 50 --order-gs 60 --ene-range -200.0 1200.0 --resol 1401 --npts 100001 --fd-npts 20001 --fd-stencil 5 --fd-margin 2.0 --fd-fit parabola --n-levels 3
python testing/plot_ccd_wavefunctions.py --wf-method fd --vof-method pypl fd --temp 5.0 --sigma 8.0 25.0 --gamma 0.3 --order-es 50 --order-gs 60 --ene-range -200.0 1200.0 --resol 1401 --npts 100001 --fd-npts 20001 --fd-stencil 5 --fd-margin 2.0 --fd-fit parabola --n-levels 3
python testing/plot_ccd_wavefunctions.py --wf-method fd --vof-method pypl fd --temp 5.0 --sigma 8.0 25.0 --gamma 0.3 --order-es 50 --order-gs 60 --ene-range -200.0 1200.0 --resol 1401 --npts 100001 --fd-npts 20001 --fd-stencil 5 --fd-margin 2.0 --fd-fit spline --n-levels 3
python testing/plot_ccd_wavefunctions.py --wf-method fd --vof-method pypl fd --temp 5.0 --sigma 8.0 25.0 --gamma 0.3 --order-es 50 --order-gs 60 --ene-range -200.0 1200.0 --resol 1401 --npts 100001 --fd-npts 20001 --fd-stencil 5 --fd-margin 2.0 --fd-fit spline5-clamped --n-levels 3
python testing/plot_ccd_wavefunctions.py --wf-method fd --vof-method pypl fd --temp 5.0 --sigma 8.0 25.0 --gamma 0.3 --order-es 50 --order-gs 60 --ene-range -200.0 1200.0 --resol 1401 --npts 100001 --fd-npts 20001 --fd-stencil 5 --fd-margin 2.0 --fd-fit morse --n-levels 3
python testing/plot_ccd_wavefunctions.py --wf-method fd --vof-method pypl fd --temp 5.0 --sigma 8.0 25.0 --gamma 0.3 --order-es 50 --order-gs 60 --ene-range -200.0 1200.0 --resol 1401 --npts 100001 --fd-npts 20001 --fd-stencil 5 --fd-margin 2.0 --fd-fit poly2 --n-levels 3
python testing/plot_ccd_wavefunctions.py --wf-method fd --vof-method pypl fd --temp 5.0 --sigma 8.0 25.0 --gamma 0.3 --order-es 50 --order-gs 60 --ene-range -200.0 1200.0 --resol 1401 --npts 100001 --fd-npts 20001 --fd-stencil 5 --fd-margin 2.0 --fd-fit poly6 --n-levels 3
python testing/plot_ccd_wavefunctions.py --wf-method fd --vof-method pypl fd --temp 5.0 --sigma 8.0 25.0 --gamma 0.3 --order-es 50 --order-gs 60 --ene-range -200.0 1200.0 --resol 1401 --npts 100001 --fd-npts 20001 --fd-stencil 5 --fd-margin 2.0 --fd-fit cubic --n-levels 3
python testing/plot_ccd_wavefunctions.py --wf-method fd --vof-method pypl fd --temp 5.0 --sigma 8.0 25.0 --gamma 0.3 --order-es 50 --order-gs 60 --ene-range -200.0 1200.0 --resol 1401 --npts 100001 --fd-npts 20001 --fd-stencil 5 --fd-margin 2.0 --fd-fit spline4-notaknot --n-levels 3
python testing/plot_ccd_wavefunctions.py --wf-method fd --vof-method pypl fd --temp 5.0 --sigma 8.0 25.0 --gamma 0.3 --order-es 50 --order-gs 60 --ene-range -200.0 1200.0 --resol 1401 --npts 100001 --fd-npts 20001 --fd-stencil 5 --fd-margin 2.0 --fd-fit spline3-s0.05 --n-levels 3
python testing/plot_ccd_wavefunctions.py --wf-method cubic --vof-method pypl cubic --temp 5.0 --sigma 8.0 25.0 --gamma 0.3 --order-es 50 --order-gs 60 --ene-range -200.0 1200.0 --resol 1401 --npts 100001 --fd-npts 20001 --fd-stencil 5 --fd-margin 2.0 --fd-fit parabola --cubic-points all --n-levels 3
python testing/plot_ccd_wavefunctions.py --wf-method cubic --vof-method pypl cubic --temp 5.0 --sigma 8.0 25.0 --gamma 0.3 --order-es 50 --order-gs 60 --ene-range -200.0 1200.0 --resol 1401 --npts 100001 --fd-npts 20001 --fd-stencil 5 --fd-margin 2.0 --fd-fit parabola --cubic-points all --cubic-free-a --n-levels 3
python testing/vof_report.py --reference pypl --temp 5.0 --sigma 8.0 25.0 --gamma 0.3 --order-es 50 --order-gs 60 --ene-range -200.0 1200.0 --resol 1401 --npts 100001 --fd-npts 20001 --fd-stencil 5 --fd-margin 2.0 --fd-fit parabola --cubic-points all --tol-fc 1e-07 --tol-vof 1e-08
```
