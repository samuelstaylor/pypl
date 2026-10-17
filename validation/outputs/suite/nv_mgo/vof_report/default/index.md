# 1D CCD VOF report [NV- in MgO]

Generated 2026-10-16T21:44:07 by `python C:\Users\sammy\Research\Galli\programs\pypl\validation\vof_report.py --reference hermite --dataset nv_mgo --temp 5.0 --sigma 20.0 --gamma 0.3 --no-zpl --order-es 160 --order-gs 100 --ene-range -500.0 10000.0 --resol 2101 --line-positions harmonic --npts 100001 --hermite-margin 6.0 --fd-stencil 5 --fd-margin 2.0 --fd-fit parabola --fd-npts 40001 --cubic-warn 0.1 --dpi 300 --fig-format png --font-size 9.0 --cubic-points all --tol-fc 1e-07 --tol-vof 1e-08`.

**Dataset.**

- dataset: nv_mgo (NV- in MgO), 3E -> 1A1 (WEST TDDFT-DDH)
- data: doc/tutorials/003_nv_mgo_isc
- roles: 'ES' = 3E (upper state, minimum at Q = 0), 'GS' = 1A1 (lower state, minimum at Q = delta_q); emission = ES -> GS = 3E -> 1A1
- energies: 3E = spin-conserving TDDFT@DDH (WEST, root 1); 1A1 = spin-flip TDDFT@DDH (WEST, root 4); zero = the 3A2 minimum, DFT with the DDH functional (PBE0, exact exchange 0.34)
- Q = q delta_q, q from Image_<q>, as in reference/ccd_nv_mgo.py: 3E minimum at Q = 0, 1A1 minimum at Q = delta_q
- E = QE total energy (qe/) + WEST root (3E: sc root 1, 1A1: sf root 4) x 13.605693 eV/Ry, relative to the 3A2 (QE) minimum, as in ccd_nv_mgo.py
- image windows (q of Image_<q>, inclusive): parabola 3E q = -0.10..+0.45, 1A1 q = +0.60..+1.15; cubic a Q^2 + b Q^3 + c (perturbative correction) q = -0.10..+1.15 for both; every other fit (splines, Morse, polyN) ALL images
- all images: 3E all 33 sc images (q = -1.00..+2.00), 1A1 all 30 sf images except Image_+2.00 (not a spin-flip run, so not a 1A1 point)

**Section 07 skipped:** the cubic fit has no exact solution to compare with: the 1A1 cubic turns over at Q = -1.70, inside the grid -6.68..10.23: no bound states in the box.

CCD: delta_q = 3.097879 amu^0.5 A, GS hbar w = 31.55114 meV, ES hbar w = 60.83324 meV. lineshape: T = 5 K, sigma = [20.0, 20.0] meV, ZPL Gaussian (no Lorentzian ZPL); lines at harmonic ladders (bulid_fc_lsp); lineshape axis: -500..10000 meV from the ZPL, 2101 points, step 5 meV; basis: 160 ES x 100 GS levels (emission ni x nf = 160 x 100, absorption 100 x 160).

Report settings: levels drawn [0, 1, 2, 5, 10, 20, 30, 40]; convergence scan [2001, 4001, 8001, 16001, 32001] points, 3/5-point stencils; cubic fits cubic, cubic-near, cubic-free, cubic-near-free; b scales 1, 0.5, 0.25, 0.125, 0.0625.

**Like-for-like method test (section 05): FAIL.** Fits compared: parabola, poly2, poly4, poly6, poly8, morse, spline, spline3-notaknot, spline3-clamped, spline4-notaknot, spline5-clamped, spline3-s0.05; unbounded, so only in section 01: poly6, poly8.

Sections: `01_potential_fits/`, `02_fd_grid/`, `03_wavefunctions/`, `04_fc_matrix/`, `05_vof_methods/`, `06_vof_fits/`, `07_cubic_perturbation/`

## Potential fits

| well | fit | images | rms (meV) | max |res| (meV) | hbar w at min (meV) | dV at right turning pt (meV) | V left wall (eV) | V right wall (eV) | bounded |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1A1 | parabola (near min.) | 18-29 of 30 | 649.989 | 3130.185 | 31.551 | 0 | 13.98 | 8.657 | True |
| 1A1 | polynomial deg. 2 | 1-30 of 30 | 351.233 | 822.661 | 45.715 | 5.76e+03 | 22.2 | 18.52 | True |
| 1A1 | polynomial deg. 4 | 1-30 of 30 | 137.435 | 485.585 | 32.863 | 1.44e+04 | 56.52 | 52.81 | True |
| 1A1 | polynomial deg. 6 | 1-30 of 30 | 58.335 | 126.697 | 54.448 | -1.35e+05 | -211.1 | -747.9 | False |
| 1A1 | polynomial deg. 8 | 1-30 of 30 | 35.961 | 92.829 | 29.457 | -3.52e+03 | -1702 | -872.1 | False |
| 1A1 | Morse | 1-30 of 30 | 174.164 | 407.489 | 26.959 | -1.96e+03 | 53.41 | 4.104 | True |
| 1A1 | cubic spline (natural) | 1-30 of 30 | 0.000 | 0.000 | 40.696 | 537 | 17.11 | 9.194 | True |
| 1A1 | cubic spline (not-a-knot) | 1-30 of 30 | 0.000 | 0.000 | 40.852 | 537 | 17.11 | 9.194 | True |
| 1A1 | cubic spline (clamped) | 1-30 of 30 | 0.000 | 0.000 | 40.441 | 537 | 17.11 | 9.194 | True |
| 1A1 | quartic spline (not-a-knot) | 1-30 of 30 | 0.000 | 0.000 | 41.658 | 537 | 17.11 | 9.194 | True |
| 1A1 | quintic spline (clamped) | 1-30 of 30 | 0.000 | 0.000 | 42.045 | 537 | 17.11 | 9.194 | True |
| 1A1 | cubic smoothing spline (rms 0.05 meV) | 1-30 of 30 | 0.050 | 0.132 | 40.708 | 537 | 17.11 | 9.194 | True |
| 3E | parabola (near min.) | 4-15 of 33 | 3474.907 | 13840.022 | 60.833 | 0 | 22.92 | 49.51 | True |
| 3E | polynomial deg. 2 | 1-33 of 33 | 763.899 | 2281.609 | 36.829 | -7.39e+03 | 14.47 | 16.32 | True |
| 3E | polynomial deg. 4 | 1-33 of 33 | 189.518 | 635.517 | 61.253 | -7.92e+03 | 71.45 | 40.14 | True |
| 3E | polynomial deg. 6 | 1-33 of 33 | 55.269 | 131.229 | 70.119 | -8.25e+03 | -37.99 | 47.18 | False |
| 3E | polynomial deg. 8 | 1-33 of 33 | 5.973 | 15.139 | 71.747 | -8.12e+03 | 414.2 | -442.7 | True |
| 3E | Morse | 1-33 of 33 | 294.319 | 615.828 | 50.500 | -7.52e+03 | 125 | 6.365 | True |
| 3E | cubic spline (natural) | 1-33 of 33 | 0.000 | 0.000 | 87.696 | -8.13e+03 | 24.91 | 35.67 | True |
| 3E | cubic spline (not-a-knot) | 1-33 of 33 | 0.000 | 0.000 | 87.550 | -8.13e+03 | 24.91 | 35.67 | True |
| 3E | cubic spline (clamped) | 1-33 of 33 | 0.000 | 0.000 | 87.867 | -8.15e+03 | 24.91 | 35.67 | True |
| 3E | quartic spline (not-a-knot) | 1-33 of 33 | 0.000 | 0.000 | 84.527 | -8.13e+03 | 24.91 | 35.67 | True |
| 3E | quintic spline (clamped) | 1-33 of 33 | 0.000 | 0.000 | 84.571 | -8.15e+03 | 24.91 | 35.67 | True |
| 3E | cubic smoothing spline (rms 0.05 meV) | 1-33 of 33 | 0.050 | 0.127 | 87.328 | -8.13e+03 | 24.91 | 35.67 | True |

## Wavefunctions vs analytic Hermite (all levels)

| well | fd on fit | max infidelity | max dpsi / max chi | max ladder dev (meV) |
| --- | --- | --- | --- | --- |
| 1A1 | parabola (near min.) | 1.33e-15 | 5.46e-09 | 3.98e-07 |
| 1A1 | polynomial deg. 2 | 1.00e+00 | 1.90e+00 | 1.4e+03 |
| 1A1 | polynomial deg. 4 | 1.00e+00 | 1.86e+00 | 955 |
| 1A1 | Morse | 9.99e-01 | 1.95e+00 | 886 |
| 1A1 | cubic spline (natural) | 9.98e-01 | 1.88e+00 | 254 |
| 1A1 | cubic spline (not-a-knot) | 9.99e-01 | 1.87e+00 | 252 |
| 1A1 | cubic spline (clamped) | 1.00e+00 | 1.90e+00 | 248 |
| 1A1 | quartic spline (not-a-knot) | 9.99e-01 | 1.87e+00 | 274 |
| 1A1 | quintic spline (clamped) | 9.99e-01 | 1.90e+00 | 248 |
| 1A1 | cubic smoothing spline (rms 0.05 meV) | 9.99e-01 | 1.87e+00 | 252 |
| 3E | parabola (near min.) | 2.66e-15 | 5.95e-08 | 9.05e-06 |
| 3E | polynomial deg. 2 | 9.99e-01 | 1.73e+00 | 3.61e+03 |
| 3E | polynomial deg. 4 | 1.00e+00 | 1.86e+00 | 2.6e+03 |
| 3E | Morse | 1.00e+00 | 1.77e+00 | 2.75e+03 |
| 3E | cubic spline (natural) | 1.00e+00 | 1.81e+00 | 2.44e+03 |
| 3E | cubic spline (not-a-knot) | 1.00e+00 | 1.81e+00 | 2.46e+03 |
| 3E | cubic spline (clamped) | 1.00e+00 | 1.80e+00 | 2.5e+03 |
| 3E | quartic spline (not-a-knot) | 1.00e+00 | 1.81e+00 | 2.46e+03 |
| 3E | quintic spline (clamped) | 1.00e+00 | 1.88e+00 | 2.52e+03 |
| 3E | cubic smoothing spline (rms 0.05 meV) | 1.00e+00 | 1.81e+00 | 2.46e+03 |

## F_ij methods on the harmonic wells (reference Hermite quadrature)

| direction | test | reference | max dF_ij | at (i, j) | max dF_ij (i,j<30) | max dVOF / peak | rms dVOF / peak | F_ij time (s) | verdict |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| emission | PyPL (Ruhoff recurrences) | Hermite quadrature | 5.10e+26 | (159, 99) | 5.65e-05 | 1.98e-15 | 1.57e-16 | 0.05 | FAIL |
| emission | FD | Hermite quadrature | 1.03e-08 | (159, 95) | 1.97e-10 | 1.92e-10 | 3.36e-11 | 11.05 | OK |
| absorption | PyPL (Ruhoff recurrences) | Hermite quadrature | 2.68e+27 | (99, 159) | 1.26e-05 | 9.62e-15 | 1.89e-15 | 0.04 | FAIL |
| absorption | FD | Hermite quadrature | 1.03e-08 | (95, 159) | 1.97e-10 | 1.95e-09 | 6.05e-10 | 13.03 | OK |

## Finite differences on each fit vs the harmonic reference

Harmonic F_00^2 = 0.00000 (emission). Differences here are the anharmonicity of each fit, not errors.

| direction | fd on fit | max dF_ij | max dVOF / peak | F_00^2 | centroid shift (meV) |
| --- | --- | --- | --- | --- | --- |
| emission | parabola (near min.) | 1.03e-08 | 1.92e-10 | 0.00000 | 0.00 |
| emission | polynomial deg. 2 | 5.61e-01 | 2.18e+00 | 0.05473 | -1043.35 |
| emission | polynomial deg. 4 | 6.45e-01 | 7.46e-01 | 0.00000 | -186.04 |
| emission | Morse | 6.39e-01 | 4.95e-01 | 0.00000 | 42.37 |
| emission | cubic spline (natural) | 8.64e-01 | 5.70e+00 | 0.00000 | -572.47 |
| emission | cubic spline (not-a-knot) | 7.76e-01 | 5.09e+00 | 0.00000 | -565.85 |
| emission | cubic spline (clamped) | 8.14e-01 | 5.32e+00 | 0.00000 | -584.05 |
| emission | quartic spline (not-a-knot) | 7.73e-01 | 5.07e+00 | 0.00000 | -565.57 |
| emission | quintic spline (clamped) | 7.03e-01 | 4.82e+00 | 0.00000 | -590.54 |
| emission | cubic smoothing spline (rms 0.05 meV) | 7.79e-01 | 5.11e+00 | 0.00000 | -566.01 |
| absorption | parabola (near min.) | 1.03e-08 | 1.95e-09 | 0.00000 | 0.00 |
| absorption | polynomial deg. 2 | 5.61e-01 | 7.62e+00 | 0.05473 | -4094.87 |
| absorption | polynomial deg. 4 | 6.45e-01 | 6.02e+00 | 0.00000 | -2563.50 |
| absorption | Morse | 6.39e-01 | 3.55e+00 | 0.00000 | -3047.90 |
| absorption | cubic spline (natural) | 8.64e-01 | 1.05e+01 | 0.00000 | -2379.71 |
| absorption | cubic spline (not-a-knot) | 7.76e-01 | 9.58e+00 | 0.00000 | -2382.17 |
| absorption | cubic spline (clamped) | 8.14e-01 | 9.28e+00 | 0.00000 | -2173.41 |
| absorption | quartic spline (not-a-knot) | 7.73e-01 | 8.95e+00 | 0.00000 | -2385.15 |
| absorption | quintic spline (clamped) | 7.03e-01 | 7.66e+00 | 0.00000 | -2056.47 |
| absorption | cubic smoothing spline (rms 0.05 meV) | 7.79e-01 | 9.58e+00 | 0.00000 | -2382.17 |


## Files

| File | Content |
| --- | --- |
| `log.txt` | Full console output, including every method's settings. |
| `settings.json` | Every parameter of the run. |
| **`01_potential_fits/`** | Each fit of the computed images: residuals, implied hbar*omega, and how far it departs from the harmonic well out where the high vibrational levels live. |
| `01_potential_fits/fit_residual_rms.png` | rms residual over all images, per fit and well (log scale). Hatched = unbounded, so unusable for the finite-difference solve; fits with an rms below 1e-06 meV interpolate the images to machine precision and are marked 'exact interpolation'. |
| `01_potential_fits/fit_summary.csv` | Per well and fit: images fitted, rms / max residual, hbar*omega at the minimum, the minimum, V at the two walls of the finite-difference grid, the departure from harmonic at the outer turning points, and the full description. |
| `01_potential_fits/fits_overview.png` | Per well (columns): (a, e) the computed images (filled = used by the parabola) with every fit near the samples; (b, f) every fit over the full Q range of the finite-difference grid (-6.68..10.23), with the sampled range shaded and the grid walls in black (V at the walls is in fit_summary.csv and log.txt); (c, g) the residual at each image (symmetric-log axis when a fit misses by more than 20 meV); (d, h) V_fit - V_parabola out to the classical turning points of the highest level (dotted), symmetric-log axis. Faded curves are unbounded fits; spline variants share one colour and differ by line style. |
| `01_potential_fits/residuals_es.csv` | 3E: fit minus computed energy at every image, per fit. |
| `01_potential_fits/residuals_gs.csv` | 1A1: fit minus computed energy at every image, per fit. |
| **`02_fd_grid/`** | The finite-difference grid point by point, the wavefunctions it gives, their decay to the walls, and convergence with grid points and stencil. |
| `02_fd_grid/fd_convergence.csv` | The convergence scan behind fd_convergence.png, with timings and slopes. |
| `02_fd_grid/fd_convergence.png` | Emission, parabola wells: (a) max \|F_ij(FD) - F_ij(Hermite quadrature, 100001 points)\| and (b) the level-ladder error vs grid spacing, for the 3, 5-point stencils, with O(h^n) guides (dotted). Measured log-log slopes (F_ij): 3-point 2.00, 5-point 3.99. The F_ij timings are in fd_convergence.csv. |
| `02_fd_grid/fd_grid_absorption.png` | Absorption: (a, d) both wells on the shared grid (40001 points, h = 4.23e-04, 5-point stencil, frame x = -(Q - Q_f)) with selected levels drawn at their eigenvalues, the padding beyond the turning points of the highest level shaded, the walls in black; (b, e) a quarter of the shortest wavelength of the highest level at the well centre, point by point (383 / 218 points per wavelength); (c, f) \|psi\| / peak decaying to the walls (4.7e-16 at the walls; the ~1e-16 floor is the float64 precision of the eigensolver). F_ij took 33.0 s. |
| `02_fd_grid/fd_grid_emission.png` | Emission: (a, d) both wells on the shared grid (40001 points, h = 4.23e-04, 5-point stencil, frame x = Q - Q_f) with selected levels drawn at their eigenvalues, the padding beyond the turning points of the highest level shaded, the walls in black; (b, e) a quarter of the shortest wavelength of the highest level at the well centre, point by point (218 / 383 points per wavelength); (c, f) \|psi\| / peak decaying to the walls (6.4e-16 at the walls; the ~1e-16 floor is the float64 precision of the eigensolver). F_ij took 30.7 s. |
| `02_fd_grid/grid_summary.csv` | Per direction and well: grid points, spacing, range, resolution of the highest level, wall amplitude, orthonormality and level-ladder error. |
| **`03_wavefunctions/`** | Vibrational wavefunctions: analytic Hermite functions vs finite differences on each fit. |
| `03_wavefunctions/wavefunction_errors.csv` | Per well, fit and level: E_v - E_0, its departure from v hbar w, max \|psi - chi\| (absolute and / max \|chi\|) and the infidelity. |
| `03_wavefunctions/wavefunction_errors.png` | Per level, well and fit: (a, b) the infidelity 1 - \|<chi_v\|psi_v>\| and (c, d) max \|psi_v - chi_v\| relative to the analytic Hermite function (floored at 1e-17). On the parabola this is the numerical error of the finite-difference solve; on the other fits it is their anharmonicity. |
| `03_wavefunctions/wavefunctions_es.png` | 3E: selected levels v (rows) from the analytic Hermite function (grey) and from finite differences on every bounded fit (left), and their difference psi_v - chi_v (right, each row scaled by the power of ten on its axis). On the parabola the difference is numerical (~1e-9, invisible next to the anharmonic fits); every per-level maximum is in wavefunction_errors.csv. |
| `03_wavefunctions/wavefunctions_gs.png` | 1A1: selected levels v (rows) from the analytic Hermite function (grey) and from finite differences on every bounded fit (left), and their difference psi_v - chi_v (right, each row scaled by the power of ten on its axis). On the parabola the difference is numerical (~1e-9, invisible next to the anharmonic fits); every per-level maximum is in wavefunction_errors.csv. |
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
| `06_vof_fits/vof_fits.png` | Per direction: FD on every bounded fit overlaid on the harmonic reference (Hermite quadrature), with the signed and absolute differences relative to the reference peak. The differences are the anharmonicity of each fit, not errors. Spline variants share one colour and differ by line style. |
| **`07_cubic_perturbation/`** | The first-order cubic anharmonic correction (method cubic): the cubic fits, its validation against the exact solution of the same cubic, wavefunctions, F_ij and VOF. |
| `07_cubic_perturbation/skipped.json` | Why this section was skipped (the cubic fit is unbounded on the fd grid). |
