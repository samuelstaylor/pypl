# 1D CCD VOF report

Generated 2026-10-16T21:12:54 by `python C:\Users\sammy\Research\Galli\programs\pypl\validation\vof_report.py --reference pypl --dataset nv_diamond --temp 5.0 --sigma 8.0 25.0 --gamma 0.3 --order-es 50 --order-gs 60 --ene-range -200.0 1200.0 --resol 1401 --line-positions harmonic --npts 100001 --hermite-margin 6.0 --fd-stencil 5 --fd-margin 2.0 --fd-fit parabola --fd-npts 20001 --cubic-warn 0.1 --dpi 300 --fig-format png --font-size 9.0 --cubic-points all --tol-fc 1e-07 --tol-vof 1e-08`.

**Dataset.**

- dataset: nv_diamond (NV- in diamond), 3E -> 3A2
- data: doc/tutorials/002_nv_diamond_1d_ccd
- roles: 'ES' = ES (upper state, minimum at Q = delta_q), 'GS' = GS (lower state, minimum at Q = 0); emission = ES -> GS = ES -> GS
- energies: ES = ΔSCF, constrained occupations (PBE), 3E; GS = DFT (PBE), 3A2; zero = the GS (3A2) minimum image

CCD: delta_q = 0.695732 amu^0.5 A, GS hbar w = 59.08890 meV, ES hbar w = 64.31191 meV. lineshape: T = 5 K, sigma = [8.0, 25.0] meV, ZPL Lorentzian, gamma = 0.3 meV; lines at harmonic ladders (bulid_fc_lsp); lineshape axis: -200..1200 meV from the ZPL, 1401 points, step 1 meV; basis: 50 ES x 60 GS levels (emission ni x nf = 50 x 60, absorption 60 x 50).

Report settings: levels drawn [0, 1, 2, 5, 10, 20, 30, 40]; convergence scan [2001, 4001, 8001, 16001, 32001] points, 3/5-point stencils; cubic fits cubic, cubic-near, cubic-free, cubic-near-free; b scales 1, 0.5, 0.25, 0.125, 0.0625.

**Like-for-like method test (section 05): PASS.** Fits compared: parabola, poly2, poly4, poly6, poly8, morse, spline, spline3-notaknot, spline3-clamped, spline4-notaknot, spline5-clamped, spline3-s0.05; unbounded, so only in section 01: poly4, poly8.

Sections: `01_potential_fits/`, `02_fd_grid/`, `03_wavefunctions/`, `04_fc_matrix/`, `05_vof_methods/`, `06_vof_fits/`, `07_cubic_perturbation/`

## Potential fits

| well | fit | images | rms (meV) | max |res| (meV) | hbar w at min (meV) | dV at right turning pt (meV) | V left wall (eV) | V right wall (eV) | bounded |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| GS | parabola (near min.) | 1-5 of 15 | 2.114 | 5.079 | 59.089 | 0 | 10.03 | 11.43 | True |
| GS | polynomial deg. 2 | 1-15 of 15 | 0.152 | 0.288 | 59.703 | 70.1 | 10.25 | 11.67 | True |
| GS | polynomial deg. 4 | 1-15 of 15 | 0.052 | 0.080 | 58.691 | -589 | -4.101 | 2.066 | False |
| GS | polynomial deg. 6 | 1-15 of 15 | 0.024 | 0.054 | 58.753 | 3.3e+04 | 2813 | 1751 | True |
| GS | polynomial deg. 8 | 1-15 of 15 | 0.010 | 0.024 | 59.307 | -4.9e+05 | -2.671e+05 | -1.141e+05 | False |
| GS | Morse | 1-15 of 15 | 0.110 | 0.195 | 59.279 | 166 | 9.448 | 12.37 | True |
| GS | cubic spline (natural) | 1-15 of 15 | 0.000 | 0.000 | 56.425 | 5.08 | 10.03 | 11.44 | True |
| GS | cubic spline (not-a-knot) | 1-15 of 15 | 0.000 | 0.000 | 58.318 | 5.08 | 10.03 | 11.44 | True |
| GS | cubic spline (clamped) | 1-15 of 15 | 0.000 | 0.000 | 58.353 | 5.08 | 10.03 | 11.44 | True |
| GS | quartic spline (not-a-knot) | 1-15 of 15 | 0.000 | 0.000 | 58.556 | 5.08 | 10.03 | 11.44 | True |
| GS | quintic spline (clamped) | 1-15 of 15 | 0.000 | 0.000 | 58.579 | 5.08 | 10.03 | 11.44 | True |
| GS | cubic smoothing spline (rms 0.05 meV) | 1-15 of 15 | 0.050 | 0.093 | 58.873 | 5.1 | 10.03 | 11.44 | True |
| ES | parabola (near min.) | 5-11 of 11 | 0.842 | 2.141 | 64.312 | 0 | 17.2 | 11.89 | True |
| ES | polynomial deg. 2 | 1-11 of 11 | 0.077 | 0.142 | 65.029 | 76.8 | 17.54 | 12.13 | True |
| ES | polynomial deg. 4 | 1-11 of 11 | 0.042 | 0.082 | 64.630 | -626 | 10.61 | 6.394 | False |
| ES | polynomial deg. 6 | 1-11 of 11 | 0.006 | 0.013 | 64.357 | 1.44e+05 | 6230 | 3630 | True |
| ES | polynomial deg. 8 | 1-11 of 11 | 0.004 | 0.009 | 64.317 | -9.34e+06 | -2.368e+06 | -7.914e+05 | False |
| ES | Morse | 1-11 of 11 | 0.044 | 0.069 | 64.716 | -79.7 | 18.83 | 11.33 | True |
| ES | cubic spline (natural) | 1-11 of 11 | 0.000 | 0.000 | 62.157 | 0.326 | 17.21 | 11.89 | True |
| ES | cubic spline (not-a-knot) | 1-11 of 11 | 0.000 | 0.000 | 64.471 | 0.326 | 17.21 | 11.89 | True |
| ES | cubic spline (clamped) | 1-11 of 11 | 0.000 | 0.000 | 63.718 | 0.326 | 17.21 | 11.89 | True |
| ES | quartic spline (not-a-knot) | 1-11 of 11 | 0.000 | 0.000 | 64.526 | 0.326 | 17.21 | 11.89 | True |
| ES | quintic spline (clamped) | 1-11 of 11 | 0.000 | 0.000 | 63.471 | 0.326 | 17.21 | 11.89 | True |
| ES | cubic smoothing spline (rms 0.05 meV) | 1-11 of 11 | 0.044 | 0.069 | 64.715 | 0.305 | 17.21 | 11.89 | True |

## Wavefunctions vs analytic Hermite (all levels)

| well | fd on fit | max infidelity | max dpsi / max chi | max ladder dev (meV) |
| --- | --- | --- | --- | --- |
| GS | parabola (near min.) | 6.66e-16 | 1.59e-08 | 1.99e-06 |
| GS | polynomial deg. 2 | 2.42e-02 | 2.24e-01 | 36.2 |
| GS | polynomial deg. 6 | 9.98e-01 | 1.92e+00 | 3.1e+03 |
| GS | Morse | 6.63e-01 | 8.47e-01 | 9.85 |
| GS | cubic spline (natural) | 7.62e-04 | 4.24e-02 | 2.11 |
| GS | cubic spline (not-a-knot) | 7.62e-04 | 4.23e-02 | 2.12 |
| GS | cubic spline (clamped) | 7.62e-04 | 4.24e-02 | 2.12 |
| GS | quartic spline (not-a-knot) | 7.62e-04 | 4.23e-02 | 2.12 |
| GS | quintic spline (clamped) | 7.62e-04 | 4.25e-02 | 2.12 |
| GS | cubic smoothing spline (rms 0.05 meV) | 7.68e-04 | 4.25e-02 | 2.13 |
| ES | parabola (near min.) | 6.66e-16 | 8.20e-09 | 1.09e-06 |
| ES | polynomial deg. 2 | 2.05e-02 | 2.22e-01 | 35.1 |
| ES | polynomial deg. 6 | 9.98e-01 | 1.93e+00 | 2.87e+03 |
| ES | Morse | 4.39e-01 | 7.56e-01 | 18.6 |
| ES | cubic spline (natural) | 7.98e-05 | 1.74e-02 | 1.02 |
| ES | cubic spline (not-a-knot) | 7.98e-05 | 1.72e-02 | 1.03 |
| ES | cubic spline (clamped) | 7.98e-05 | 1.73e-02 | 1.03 |
| ES | quartic spline (not-a-knot) | 7.98e-05 | 1.72e-02 | 1.03 |
| ES | quintic spline (clamped) | 7.98e-05 | 1.73e-02 | 1.03 |
| ES | cubic smoothing spline (rms 0.05 meV) | 8.55e-05 | 1.76e-02 | 1.04 |

## F_ij methods on the harmonic wells (reference PyPL (Ruhoff recurrences))

| direction | test | reference | max dF_ij | at (i, j) | max dF_ij (i,j<30) | max dVOF / peak | rms dVOF / peak | F_ij time (s) | verdict |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| emission | Hermite quadrature | PyPL (Ruhoff recurrences) | 1.70e-08 | (49, 48) | 4.87e-11 | 6.95e-16 | 3.97e-17 | 0.68 | OK |
| emission | FD | PyPL (Ruhoff recurrences) | 1.80e-08 | (49, 49) | 6.11e-10 | 1.28e-12 | 3.82e-13 | 5.41 | OK |
| absorption | Hermite quadrature | PyPL (Ruhoff recurrences) | 1.37e-08 | (50, 49) | 3.71e-11 | 6.95e-16 | 4.46e-17 | 0.68 | OK |
| absorption | FD | PyPL (Ruhoff recurrences) | 1.51e-08 | (50, 49) | 6.11e-10 | 9.52e-13 | 3.37e-13 | 6.09 | OK |

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
| emission | First-order cubic PT | 1.16e+00 | 1.62e-01 | 3.25e-04 |
| emission | PyPL (Ruhoff recurrences), harmonic | 4.89e-01 | 2.73e-01 | 3.50e-02 |
| absorption | First-order cubic PT | 1.16e+00 | 1.62e-01 | 3.25e-04 |
| absorption | PyPL (Ruhoff recurrences), harmonic | 4.89e-01 | 2.73e-01 | 3.50e-02 |

b -> 0 validation: the PT error falls as b^1.99 (F_ij, i, j < 8) and b^2.00 (VOF), the anharmonic effect as b^1.00; at b = 0 the method equals the Hermite route to 8.3e-16. The mixing 3|xi|(n+1)^1.5 exceeds 0.1 above n = 9 (GS) / 8 (ES), so the high-level F_ij are not perturbative; the VOF is, because the intensity sits in low levels.


## Files

| File | Content |
| --- | --- |
| `log.txt` | Full console output, including every method's settings. |
| `settings.json` | Every parameter of the run. |
| **`01_potential_fits/`** | Each fit of the computed images: residuals, implied hbar*omega, and how far it departs from the harmonic well out where the high vibrational levels live. |
| `01_potential_fits/fit_residual_rms.png` | rms residual over all images, per fit and well (log scale). Hatched = unbounded, so unusable for the finite-difference solve; fits with an rms below 1e-06 meV interpolate the images to machine precision and are marked 'exact interpolation'. |
| `01_potential_fits/fit_summary.csv` | Per well and fit: images fitted, rms / max residual, hbar*omega at the minimum, the minimum, V at the two walls of the finite-difference grid, the departure from harmonic at the outer turning points, and the full description. |
| `01_potential_fits/fits_overview.png` | Per well (columns): (a, e) the computed images (filled = used by the parabola) with every fit near the samples; (b, f) every fit over the full Q range of the finite-difference grid (-4.90..5.23), with the sampled range shaded and the grid walls in black (V at the walls is in fit_summary.csv and log.txt); (c, g) the residual at each image (symmetric-log axis when a fit misses by more than 20 meV); (d, h) V_fit - V_parabola out to the classical turning points of the highest level (dotted), symmetric-log axis. Faded curves are unbounded fits; spline variants share one colour and differ by line style. |
| `01_potential_fits/residuals_es.csv` | ES: fit minus computed energy at every image, per fit. |
| `01_potential_fits/residuals_gs.csv` | GS: fit minus computed energy at every image, per fit. |
| **`02_fd_grid/`** | The finite-difference grid point by point, the wavefunctions it gives, their decay to the walls, and convergence with grid points and stencil. |
| `02_fd_grid/fd_convergence.csv` | The convergence scan behind fd_convergence.png, with timings and slopes. |
| `02_fd_grid/fd_convergence.png` | Emission, parabola wells: (a) max \|F_ij(FD) - F_ij(Hermite quadrature, 100001 points)\| and (b) the level-ladder error vs grid spacing, for the 3, 5-point stencils, with O(h^n) guides (dotted). Measured log-log slopes (F_ij): 3-point 2.00, 5-point 4.00. The F_ij timings are in fd_convergence.csv. |
| `02_fd_grid/fd_grid_absorption.png` | Absorption: (a, d) both wells on the shared grid (20001 points, h = 5.07e-04, 5-point stencil, frame x = Q - Q_f) with selected levels drawn at their eigenvalues, the padding beyond the turning points of the highest level shaded, the walls in black; (b, e) a quarter of the shortest wavelength of the highest level at the well centre, point by point (302 / 318 points per wavelength); (c, f) \|psi\| / peak decaying to the walls (5.8e-16 at the walls; the ~1e-16 floor is the float64 precision of the eigensolver). F_ij took 2.3 s. |
| `02_fd_grid/fd_grid_emission.png` | Emission: (a, d) both wells on the shared grid (20001 points, h = 5.07e-04, 5-point stencil, frame x = -(Q - Q_f)) with selected levels drawn at their eigenvalues, the padding beyond the turning points of the highest level shaded, the walls in black; (b, e) a quarter of the shortest wavelength of the highest level at the well centre, point by point (318 / 302 points per wavelength); (c, f) \|psi\| / peak decaying to the walls (3.8e-16 at the walls; the ~1e-16 floor is the float64 precision of the eigensolver). F_ij took 0.9 s. |
| `02_fd_grid/grid_summary.csv` | Per direction and well: grid points, spacing, range, resolution of the highest level, wall amplitude, orthonormality and level-ladder error. |
| **`03_wavefunctions/`** | Vibrational wavefunctions: analytic Hermite functions vs finite differences on each fit. |
| `03_wavefunctions/wavefunction_errors.csv` | Per well, fit and level: E_v - E_0, its departure from v hbar w, max \|psi - chi\| (absolute and / max \|chi\|) and the infidelity. |
| `03_wavefunctions/wavefunction_errors.png` | Per level, well and fit: (a, b) the infidelity 1 - \|<chi_v\|psi_v>\| and (c, d) max \|psi_v - chi_v\| relative to the analytic Hermite function (floored at 1e-17). On the parabola this is the numerical error of the finite-difference solve; on the other fits it is their anharmonicity. |
| `03_wavefunctions/wavefunctions_es.png` | ES: selected levels v (rows) from the analytic Hermite function (grey) and from finite differences on every bounded fit (left), and their difference psi_v - chi_v (right, each row scaled by the power of ten on its axis). On the parabola the difference is numerical (~1e-9, invisible next to the anharmonic fits); every per-level maximum is in wavefunction_errors.csv. |
| `03_wavefunctions/wavefunctions_gs.png` | GS: selected levels v (rows) from the analytic Hermite function (grey) and from finite differences on every bounded fit (left), and their difference psi_v - chi_v (right, each row scaled by the power of ten on its axis). On the parabola the difference is numerical (~1e-9, invisible next to the anharmonic fits); every per-level maximum is in wavefunction_errors.csv. |
| **`04_fc_matrix/`** | The Franck-Condon matrix of each method and its error against the reference. |
| `04_fc_matrix/fc_error_by_level.png` | Largest F_ij difference from the reference in each row i, per method. With the PyPL reference the difference grows with the level: the float64 drift of its recurrence. |
| `04_fc_matrix/fc_error_maps.png` | log10 \|F_ij(method) - F_ij(reference)\| per method and direction; the cross marks the maximum. |
| `04_fc_matrix/fc_errors.csv` | Every F_ij and VOF error statistic per method and direction. |
| `04_fc_matrix/fc_fd_absorption.csv` | F_ij from FD, absorption (rows i = initial level, columns j = final level). |
| `04_fc_matrix/fc_fd_emission.csv` | F_ij from FD, emission (rows i = initial level, columns j = final level). |
| `04_fc_matrix/fc_hermite_absorption.csv` | F_ij from Hermite quadrature, absorption (rows i = initial level, columns j = final level). |
| `04_fc_matrix/fc_hermite_emission.csv` | F_ij from Hermite quadrature, emission (rows i = initial level, columns j = final level). |
| `04_fc_matrix/fc_matrix.png` | The reference Franck-Condon matrix, log10 \|F_ij\|, per direction (clipped at 1e-12). |
| `04_fc_matrix/fc_pypl_absorption.csv` | F_ij from PyPL (Ruhoff recurrences), absorption (rows i = initial level, columns j = final level). |
| `04_fc_matrix/fc_pypl_emission.csv` | F_ij from PyPL (Ruhoff recurrences), emission (rows i = initial level, columns j = final level). |
| **`05_vof_methods/`** | VOF of every F_ij method on the harmonic wells: the like-for-like test. |
| `05_vof_methods/summary.csv` | Every error statistic per method and direction, with OK / FAIL. |
| `05_vof_methods/vof_absorption.csv` | absorption: every method's lineshape and its difference from the reference. |
| `05_vof_methods/vof_emission.csv` | emission: every method's lineshape and its difference from the reference. |
| `05_vof_methods/vof_methods.png` | Per direction: every method's VOF overlaid on the reference, and the signed and absolute error relative to the reference peak, with the tolerance (dashed). Every method's settings are in log.txt. |
| **`06_vof_fits/`** | VOF of the finite-difference method on each fit, against the harmonic reference. |
| `06_vof_fits/level_ladders.png` | Per well and fit: how far the FD level ladder departs from the harmonic v hbar omega (symmetric-log axis, linear within +-1 meV). With the default line positions the VOF still places its lines on the harmonic ladder; only F_ij follows the fit. |
| `06_vof_fits/summary.csv` | Per fit and direction: F_ij and VOF difference from the harmonic reference, ZPL weight F_00^2, shift of the lineshape centroid. |
| `06_vof_fits/vof_absorption.csv` | absorption: VOF of the reference and of fd on every fit. |
| `06_vof_fits/vof_emission.csv` | emission: VOF of the reference and of fd on every fit. |
| `06_vof_fits/vof_fits.png` | Per direction: FD on every bounded fit overlaid on the harmonic reference (PyPL (Ruhoff recurrences)), with the signed and absolute differences relative to the reference peak. The differences are the anharmonicity of each fit, not errors. Spline variants share one colour and differ by line style. |
| **`07_cubic_perturbation/`** | The first-order cubic anharmonic correction (method cubic): the cubic fits, its validation against the exact solution of the same cubic, wavefunctions, F_ij and VOF. |
| `07_cubic_perturbation/cubic_fits.csv` | Per well and cubic fit: a, b, turnover, xi, the largest first-order mixing 3\|xi\|(n+1)^1.5, the highest level where it stays <= 0.1, residuals. |
| `07_cubic_perturbation/cubic_fits.png` | Per well: the cubic-fit options (cubic, cubic-near, cubic-free, cubic-near-free; a fixed to the parabola's or fitted, all images or the parabola's) against the images (a, d), over the full FD grid with its walls (b, e), and their residuals (c, f). The fit used by the method (cubic) is drawn thicker; b, xi, turnover and mixing are in cubic_fits.csv. |
| `07_cubic_perturbation/cubic_wavefunctions_es.csv` | ES: per level, the first-order mixing and max \|PT - exact\| and \|harmonic - exact\|. |
| `07_cubic_perturbation/cubic_wavefunctions_gs.csv` | GS: per level, the first-order mixing and max \|PT - exact\| and \|harmonic - exact\|. |
| `07_cubic_perturbation/fc_error_maps.png` | log10 \|F_ij - F_ij(exact cubic)\| of the first-order correction and of the harmonic reference. |
| `07_cubic_perturbation/summary.csv` | F_ij and VOF differences from the exact cubic solution, per method and direction. |
| `07_cubic_perturbation/validation.png` | Emission, against the exact solution of the same cubic (FD): (a) max \|difference\| of F_ij (i, j < 8, basis 12 x 12) and of the VOF / peak as b is scaled towards 0; the PT error should fall as b^2 (slope 2) and the anharmonic effect (harmonic vs exact) as b (slope 1); at b = 0, max \|F_PT - F_Hermite\| = 8.3e-16. (b) norm - 1 of the corrected states vs the predicted sum of squared PT coefficients - 1. (c) the first-order mixing 3\|xi\|(n+1)^1.5 and the infidelity against the exact eigenstates per level: first order breaks down where the mixing passes 0.1 (dashed). |
| `07_cubic_perturbation/validation_b_scaling.csv` | The b-scaling scan behind validation.png. |
| `07_cubic_perturbation/vof_absorption.csv` | absorption: VOF of the exact cubic solution, the first-order correction and the harmonic reference. |
| `07_cubic_perturbation/vof_cubic.png` | Per direction: the exact VOF on the cubic wells (FD), the first-order correction and the harmonic reference, with signed and absolute differences relative to the exact peak. |
| `07_cubic_perturbation/vof_emission.csv` | emission: VOF of the exact cubic solution, the first-order correction and the harmonic reference. |
| `07_cubic_perturbation/wavefunctions_es.png` | ES, cubic (a fixed, all images): harmonic, first-order corrected and exact (FD on the same cubic) wavefunctions for selected levels (left) and their differences from exact (right, each row scaled by the power of ten on its axis). The mixing 3\|xi\|(n+1)^1.5 and the maxima per level are in cubic_wavefunctions_es.csv. |
| `07_cubic_perturbation/wavefunctions_gs.png` | GS, cubic (a fixed, all images): harmonic, first-order corrected and exact (FD on the same cubic) wavefunctions for selected levels (left) and their differences from exact (right, each row scaled by the power of ten on its axis). The mixing 3\|xi\|(n+1)^1.5 and the maxima per level are in cubic_wavefunctions_gs.csv. |
