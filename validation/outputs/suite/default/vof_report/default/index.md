# 1D CCD VOF report

Generated 2026-10-15T18:50:48 by `python C:\Users\sammy\Research\Galli\programs\pypl\testing\vof_report.py --reference pypl --temp 5.0 --sigma 8.0 25.0 --gamma 0.3 --order-es 50 --order-gs 60 --ene-range -200.0 1200.0 --resol 1401 --npts 100001 --fd-npts 20001 --fd-stencil 5 --fd-margin 2.0 --fd-fit parabola --cubic-points all --tol-fc 1e-07 --tol-vof 1e-08`.

CCD: delta_q = 0.695732 amu^0.5 A, GS hbar w = 59.08890 meV, ES hbar w = 64.31191 meV. bulid_fc_lsp: T = 5 K, sigma = [8.0, 25.0] meV, ZPL Lorentzian, gamma = 0.3 meV; lineshape axis: -200..1200 meV from the ZPL, 1401 points, step 1 meV; basis: 50 ES x 60 GS levels (emission ni x nf = 50 x 60, absorption 60 x 50).

**Like-for-like method test (section 05): PASS.** Fits compared: parabola, poly2, poly4, poly6, poly8, morse, spline, spline3-notaknot, spline3-clamped, spline4-notaknot, spline5-clamped, spline3-s0.05; unbounded, so only in section 01: poly4, poly8.

Sections: `01_potential_fits/`, `02_fd_grid/`, `03_wavefunctions/`, `04_fc_matrix/`, `05_vof_methods/`, `06_vof_fits/`, `07_cubic_perturbation/`

## Potential fits

| well | fit | images | rms (meV) | max |res| (meV) | hbar w at min (meV) | dV at right turning pt (meV) | bounded |
| --- | --- | --- | --- | --- | --- | --- | --- |
| GS | parabola (near min.) | 1-5 of 15 | 2.114 | 5.079 | 59.089 | 0 | True |
| GS | polynomial deg. 2 | 1-15 of 15 | 0.152 | 0.288 | 59.703 | 70.1 | True |
| GS | polynomial deg. 4 | 1-15 of 15 | 0.052 | 0.080 | 58.691 | -589 | False |
| GS | polynomial deg. 6 | 1-15 of 15 | 0.024 | 0.054 | 58.753 | 3.3e+04 | True |
| GS | polynomial deg. 8 | 1-15 of 15 | 0.010 | 0.024 | 59.307 | -4.9e+05 | False |
| GS | Morse | 1-15 of 15 | 0.110 | 0.195 | 59.279 | 166 | True |
| GS | cubic spline (natural) | 1-15 of 15 | 0.000 | 0.000 | 56.425 | 5.08 | True |
| GS | cubic spline (not-a-knot) | 1-15 of 15 | 0.000 | 0.000 | 58.318 | 5.08 | True |
| GS | cubic spline (clamped) | 1-15 of 15 | 0.000 | 0.000 | 58.353 | 5.08 | True |
| GS | quartic spline (not-a-knot) | 1-15 of 15 | 0.000 | 0.000 | 58.556 | 5.08 | True |
| GS | quintic spline (clamped) | 1-15 of 15 | 0.000 | 0.000 | 58.579 | 5.08 | True |
| GS | cubic smoothing spline (rms 0.05 meV) | 1-15 of 15 | 0.050 | 0.093 | 58.873 | 5.1 | True |
| ES | parabola (near min.) | 5-11 of 11 | 0.842 | 2.141 | 64.312 | 0 | True |
| ES | polynomial deg. 2 | 1-11 of 11 | 0.077 | 0.142 | 65.029 | 76.8 | True |
| ES | polynomial deg. 4 | 1-11 of 11 | 0.042 | 0.082 | 64.630 | -626 | False |
| ES | polynomial deg. 6 | 1-11 of 11 | 0.006 | 0.013 | 64.357 | 1.44e+05 | True |
| ES | polynomial deg. 8 | 1-11 of 11 | 0.004 | 0.009 | 64.317 | -9.34e+06 | False |
| ES | Morse | 1-11 of 11 | 0.044 | 0.069 | 64.716 | -79.7 | True |
| ES | cubic spline (natural) | 1-11 of 11 | 0.000 | 0.000 | 62.157 | 0.326 | True |
| ES | cubic spline (not-a-knot) | 1-11 of 11 | 0.000 | 0.000 | 64.471 | 0.326 | True |
| ES | cubic spline (clamped) | 1-11 of 11 | 0.000 | 0.000 | 63.718 | 0.326 | True |
| ES | quartic spline (not-a-knot) | 1-11 of 11 | 0.000 | 0.000 | 64.526 | 0.326 | True |
| ES | quintic spline (clamped) | 1-11 of 11 | 0.000 | 0.000 | 63.471 | 0.326 | True |
| ES | cubic smoothing spline (rms 0.05 meV) | 1-11 of 11 | 0.044 | 0.069 | 64.715 | 0.305 | True |

## Wavefunctions vs analytic Hermite (all levels)

| well | fd on fit | max infidelity | max dpsi / max chi | max ladder dev (meV) |
| --- | --- | --- | --- | --- |
| GS | parabola (near min.) | 7.77e-16 | 1.59e-08 | 1.99e-06 |
| GS | polynomial deg. 2 | 2.42e-02 | 2.24e-01 | 36.2 |
| GS | polynomial deg. 6 | 9.98e-01 | 1.92e+00 | 3.1e+03 |
| GS | Morse | 6.63e-01 | 8.47e-01 | 9.85 |
| GS | cubic spline (natural) | 7.62e-04 | 4.24e-02 | 2.11 |
| GS | cubic spline (not-a-knot) | 7.62e-04 | 4.23e-02 | 2.12 |
| GS | cubic spline (clamped) | 7.62e-04 | 4.24e-02 | 2.12 |
| GS | quartic spline (not-a-knot) | 7.62e-04 | 4.23e-02 | 2.12 |
| GS | quintic spline (clamped) | 7.62e-04 | 4.25e-02 | 2.12 |
| GS | cubic smoothing spline (rms 0.05 meV) | 7.68e-04 | 4.25e-02 | 2.13 |
| ES | parabola (near min.) | 8.88e-16 | 8.20e-09 | 1.09e-06 |
| ES | polynomial deg. 2 | 2.05e-02 | 2.22e-01 | 35.1 |
| ES | polynomial deg. 6 | 9.98e-01 | 1.93e+00 | 2.87e+03 |
| ES | Morse | 4.39e-01 | 7.56e-01 | 18.6 |
| ES | cubic spline (natural) | 7.98e-05 | 1.74e-02 | 1.02 |
| ES | cubic spline (not-a-knot) | 7.98e-05 | 1.72e-02 | 1.03 |
| ES | cubic spline (clamped) | 7.98e-05 | 1.73e-02 | 1.03 |
| ES | quartic spline (not-a-knot) | 7.98e-05 | 1.72e-02 | 1.03 |
| ES | quintic spline (clamped) | 7.98e-05 | 1.73e-02 | 1.03 |
| ES | cubic smoothing spline (rms 0.05 meV) | 8.55e-05 | 1.76e-02 | 1.04 |

## F_ij methods on the harmonic wells (reference PyPL (Ruhoff))

| direction | test | reference | max dF_ij | at (i, j) | max dF_ij (i,j<30) | max dVOF / peak | rms dVOF / peak | F_ij time (s) | verdict |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| emission | Hermite (quadrature) | PyPL (Ruhoff) | 1.70e-08 | (49, 48) | 4.87e-11 | 6.95e-16 | 4.51e-17 | 3.73 | OK |
| emission | Finite difference | PyPL (Ruhoff) | 1.80e-08 | (49, 49) | 6.11e-10 | 1.28e-12 | 3.82e-13 | 43.13 | OK |
| absorption | Hermite (quadrature) | PyPL (Ruhoff) | 1.37e-08 | (50, 49) | 3.71e-11 | 6.95e-16 | 3.83e-17 | 1.67 | OK |
| absorption | Finite difference | PyPL (Ruhoff) | 1.51e-08 | (50, 49) | 6.11e-10 | 9.52e-13 | 3.37e-13 | 23.92 | OK |

## Finite differences on each fit vs the harmonic reference

Harmonic F_00^2 = 0.02825 (emission). Differences here are the anharmonicity of each fit, not errors.

| direction | fd on fit | max dF_ij | max dVOF / peak | F_00^2 | centroid shift (meV) |
| --- | --- | --- | --- | --- | --- |
| emission | parabola (near min.) | 1.80e-08 | 1.29e-12 | 0.02825 | 0.00 |
| emission | polynomial deg. 2 | 3.36e-03 | 2.42e-03 | 0.02818 | 0.07 |
| emission | polynomial deg. 6 | 7.06e-01 | 6.96e-02 | 0.03021 | -19.30 |
| emission | Morse | 4.01e-01 | 1.95e-02 | 0.02770 | 2.40 |
| emission | cubic spline (natural) | 1.79e-02 | 2.08e-02 | 0.02766 | 2.24 |
| emission | cubic spline (not-a-knot) | 1.80e-02 | 2.18e-02 | 0.02763 | 2.26 |
| emission | cubic spline (clamped) | 1.80e-02 | 2.17e-02 | 0.02763 | 2.26 |
| emission | quartic spline (not-a-knot) | 1.80e-02 | 2.18e-02 | 0.02763 | 2.26 |
| emission | quintic spline (clamped) | 1.79e-02 | 2.16e-02 | 0.02763 | 2.26 |
| emission | cubic smoothing spline (rms 0.05 meV) | 1.81e-02 | 2.21e-02 | 0.02762 | 2.27 |
| absorption | parabola (near min.) | 1.51e-08 | 9.57e-13 | 0.02825 | 0.00 |
| absorption | polynomial deg. 2 | 3.36e-03 | 2.42e-03 | 0.02818 | 0.26 |
| absorption | polynomial deg. 6 | 7.06e-01 | 6.96e-02 | 0.03021 | -24.27 |
| absorption | Morse | 4.01e-01 | 1.95e-02 | 0.02770 | 3.98 |
| absorption | cubic spline (natural) | 1.79e-02 | 2.08e-02 | 0.02766 | 1.87 |
| absorption | cubic spline (not-a-knot) | 1.80e-02 | 2.18e-02 | 0.02763 | 1.90 |
| absorption | cubic spline (clamped) | 1.80e-02 | 2.17e-02 | 0.02763 | 1.90 |
| absorption | quartic spline (not-a-knot) | 1.80e-02 | 2.18e-02 | 0.02763 | 1.90 |
| absorption | quintic spline (clamped) | 1.79e-02 | 2.16e-02 | 0.02763 | 1.91 |
| absorption | cubic smoothing spline (rms 0.05 meV) | 1.81e-02 | 2.21e-02 | 0.02762 | 1.91 |

## First-order cubic correction (cubic (a fixed, all images))

| well | cubic fit | b (eV/amu^1.5 A^3) | xi | max mixing | n_valid | turnover Q |
| --- | --- | --- | --- | --- | --- | --- |
| GS | cubic (a fixed, all images) | +0.00916 | +1.032e-03 | 1.44 | 9 | -30.4 |
| ES | cubic (a fixed, all images) | -0.01265 | -1.153e-03 | 1.22 | 8 | 26.8 |

| direction | vs exact (fd on the cubic) | max dF_ij | max dF_ij (i,j<30) | max dVOF / peak |
| --- | --- | --- | --- | --- |
| emission | Hermite + cubic PT [cubic (a fixed, all images)] | 1.16e+00 | 1.62e-01 | 3.25e-04 |
| emission | PyPL (Ruhoff) | 4.89e-01 | 2.73e-01 | 3.50e-02 |
| absorption | Hermite + cubic PT [cubic (a fixed, all images)] | 1.16e+00 | 1.62e-01 | 3.25e-04 |
| absorption | PyPL (Ruhoff) | 4.89e-01 | 2.73e-01 | 3.50e-02 |

b -> 0 validation: the PT error falls as b^1.99 (F_ij, i, j < 8) and b^2.00 (VOF), the anharmonic effect as b^1.00; at b = 0 the method equals the Hermite route to 6.7e-16. The mixing 3|xi|(n+1)^1.5 exceeds 0.1 above n = 9 (GS) / 8 (ES), so the high-level F_ij are not perturbative; the VOF is, because the intensity sits in low levels.


## Files

| File | Content |
| --- | --- |
| `log.txt` | Full console output, including every method's settings. |
| `settings.json` | Every parameter of the run. |
| **`01_potential_fits/`** | Each fit of the DFT images: residuals, implied hbar*omega, and how far it departs from the harmonic well out where the high vibrational levels live. |
| `01_potential_fits/fit_residual_rms.png` | rms residual over all images, per fit and well. |
| `01_potential_fits/fit_summary.csv` | Per well and fit: images fitted, rms / max residual, hbar*omega at the minimum, the minimum, the departure from harmonic at the outer turning points, and the full description. |
| `01_potential_fits/fits_overview.png` | Per well: the DFT images with every fit near the samples (filled points are the parabola's); every fit over the full Q range of the finite-difference grid, with the DFT images, the sampled range (gold), the walls (black) and V at the walls; the residual at each image; and V_fit - V_parabola out to the turning points of the highest level (symlog). Dashed = unbounded. |
| `01_potential_fits/residuals_es.csv` | ES: fit - DFT at every image, per fit. |
| `01_potential_fits/residuals_gs.csv` | GS: fit - DFT at every image, per fit. |
| **`02_fd_grid/`** | The finite-difference grid point by point, the wavefunctions it gives, their decay to the walls, and convergence with grid points and stencil. |
| `02_fd_grid/fd_convergence.csv` | The convergence scan behind fd_convergence.png. |
| `02_fd_grid/fd_convergence.png` | Emission: F_ij error against the Hermite quadrature, level ladder error and time vs grid spacing, for the 3- and 5-point stencils, with O(h^2) / O(h^4) guides. |
| `02_fd_grid/fd_grid_absorption.png` | absorption: both wells on the grid with levels at their eigenvalues, the highest level point by point near its turning point, and \|psi\| decaying to the walls. |
| `02_fd_grid/fd_grid_emission.png` | emission: both wells on the grid with levels at their eigenvalues, the highest level point by point near its turning point, and \|psi\| decaying to the walls. |
| `02_fd_grid/grid_summary.csv` | Per direction and well: grid points, spacing, range, resolution of the highest level, wall amplitude, orthonormality and level-ladder error. |
| **`03_wavefunctions/`** | Vibrational wavefunctions: analytic Hermite functions vs finite differences on each fit. |
| `03_wavefunctions/wavefunction_errors.csv` | Per well, fit and level: E_v - E_0, its departure from v hbar w, max \|psi - chi\| / max \|chi\| and the infidelity. |
| `03_wavefunctions/wavefunction_errors.png` | Per level, well and fit: infidelity 1 - \|<chi\|psi>\| and max \|psi - chi\| relative to the Hermite function. |
| `03_wavefunctions/wavefunctions_es.png` | ES: selected levels from the analytic Hermite function and from finite differences on every bounded fit, and their difference. |
| `03_wavefunctions/wavefunctions_gs.png` | GS: selected levels from the analytic Hermite function and from finite differences on every bounded fit, and their difference. |
| **`04_fc_matrix/`** | The Franck-Condon matrix of each method and its error against the reference. |
| `04_fc_matrix/fc_error_by_level.png` | Largest F_ij difference in each row i, per method. |
| `04_fc_matrix/fc_error_maps.png` | log10 \|F_ij(method) - F_ij(reference)\| per method and direction. |
| `04_fc_matrix/fc_errors.csv` | Every F_ij and VOF error statistic per method and direction. |
| `04_fc_matrix/fc_fd_absorption.csv` | F_ij from Finite difference, absorption (rows i = initial level, columns j = final level). |
| `04_fc_matrix/fc_fd_emission.csv` | F_ij from Finite difference, emission (rows i = initial level, columns j = final level). |
| `04_fc_matrix/fc_hermite_absorption.csv` | F_ij from Hermite (quadrature), absorption (rows i = initial level, columns j = final level). |
| `04_fc_matrix/fc_hermite_emission.csv` | F_ij from Hermite (quadrature), emission (rows i = initial level, columns j = final level). |
| `04_fc_matrix/fc_matrix.png` | The reference Franck-Condon matrix, log10 \|F_ij\|, per direction. |
| `04_fc_matrix/fc_pypl_absorption.csv` | F_ij from PyPL (Ruhoff), absorption (rows i = initial level, columns j = final level). |
| `04_fc_matrix/fc_pypl_emission.csv` | F_ij from PyPL (Ruhoff), emission (rows i = initial level, columns j = final level). |
| **`05_vof_methods/`** | VOF of every F_ij method on the harmonic wells: the like-for-like test. |
| `05_vof_methods/summary.csv` | Every error statistic per method and direction, with OK / FAIL. |
| `05_vof_methods/vof_absorption.csv` | absorption: every method's lineshape and its difference from the reference. |
| `05_vof_methods/vof_emission.csv` | emission: every method's lineshape and its difference from the reference. |
| `05_vof_methods/vof_methods.png` | Per direction: every method's VOF overlaid, signed and log error against the reference, and every method's settings. |
| **`06_vof_fits/`** | VOF of the finite-difference method on each fit, against the harmonic reference. |
| `06_vof_fits/level_ladders.png` | Per well and fit: how far the FD level ladder departs from the harmonic v hbar w. bulid_fc_lsp still uses the harmonic positions; only F_ij follows the fit. |
| `06_vof_fits/summary.csv` | Per fit and direction: F_ij and VOF difference from the harmonic reference, ZPL weight F_00^2, shift of the lineshape centroid. |
| `06_vof_fits/vof_absorption.csv` | absorption: VOF of the reference and of fd on every fit. |
| `06_vof_fits/vof_emission.csv` | emission: VOF of the reference and of fd on every fit. |
| `06_vof_fits/vof_fits.png` | Per direction: fd on every fit overlaid on the harmonic reference, with signed and log differences. |
| **`07_cubic_perturbation/`** | The first-order cubic anharmonic correction (method cubic): the cubic fits, its validation against the exact solution of the same cubic, wavefunctions, F_ij and VOF. |
| `07_cubic_perturbation/cubic_fits.csv` | Per well and cubic fit: a, b, turnover, xi, the largest first-order mixing 3\|xi\|(n+1)^1.5, the highest level where it stays <= 0.1, residuals. |
| `07_cubic_perturbation/cubic_fits.png` | Per well: the four cubic-fit options (a fixed or fitted; all images or the parabola's) against the DFT images, over the full fd grid, and their residuals. |
| `07_cubic_perturbation/fc_error_maps.png` | log10 \|F_ij - F_ij(exact cubic)\| of the first-order correction and of the harmonic reference. |
| `07_cubic_perturbation/summary.csv` | F_ij and VOF differences from the exact cubic solution, per method and direction. |
| `07_cubic_perturbation/validation.png` | b -> 0 scaling of the PT error (b^2) and of the effect (b) against fd on the same cubic; norms vs the predicted 1 + sum c^2; mixing and infidelity per level, showing where first order breaks down. |
| `07_cubic_perturbation/validation_b_scaling.csv` | The b-scaling scan behind validation.png. |
| `07_cubic_perturbation/vof_absorption.csv` | absorption: VOF of the exact cubic solution, the first-order correction and the harmonic reference. |
| `07_cubic_perturbation/vof_cubic.png` | Per direction: the exact VOF on the cubic wells (fd), the first-order correction and the harmonic reference, with signed / log differences and every setting. |
| `07_cubic_perturbation/vof_emission.csv` | emission: VOF of the exact cubic solution, the first-order correction and the harmonic reference. |
| `07_cubic_perturbation/wavefunctions_es.png` | ES: harmonic, first-order corrected and exact (fd on the same cubic) wavefunctions for selected levels, and the differences from exact. |
| `07_cubic_perturbation/wavefunctions_gs.png` | GS: harmonic, first-order corrected and exact (fd on the same cubic) wavefunctions for selected levels, and the differences from exact. |
