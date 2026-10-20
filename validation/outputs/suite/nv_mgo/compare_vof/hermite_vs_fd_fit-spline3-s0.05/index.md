# VOF method comparison: hermite_vs_fd_fit-spline3-s0.05 [NV- in MgO]

Generated 2026-10-18T00:02:49 by `python C:\Users\sammy\Research\Galli\programs\pypl\validation\compare_vof.py --reference hermite --test fd --dataset nv_mgo --temp 5.0 --sigma 20.0 --gamma 0.3 --no-zpl --order-es 160 --order-gs 100 --ene-range -500.0 10000.0 --resol 2101 --line-positions harmonic --npts 100001 --hermite-margin 6.0 --fd-stencil 5 --fd-margin 2.0 --fd-fit spline3-s0.05 --fd-npts 40001 --cubic-warn 0.1 --dpi 300 --fig-format png --font-size 9.0 --tol-fc 1e-07 --tol-vof 1e-08`.

**Dataset.**

- dataset: nv_mgo (NV- in MgO), 3E -> 1A1 (WEST TDDFT-DDH)
- data: doc/tutorials/003_nv_mgo_isc
- roles: 'ES' = 3E (upper state, minimum at Q = 0), 'GS' = 1A1 (lower state, minimum at Q = delta_q); emission = ES -> GS = 3E -> 1A1
- energies: 3E = spin-conserving TDDFT@DDH (WEST, root 1); 1A1 = spin-flip TDDFT@DDH (WEST, root 4); zero = the 3A2 minimum, DFT with the DDH functional (PBE0, exact exchange 0.34)
- Q = q delta_q, q from Image_<q>, as in reference/ccd_nv_mgo.py: 3E minimum at Q = 0, 1A1 minimum at Q = delta_q
- E = QE total energy (qe/) + WEST root (3E: sc root 1, 1A1: sf root 4) x 13.605693 eV/Ry, relative to the 3A2 (QE) minimum, as in ccd_nv_mgo.py
- image windows (q of Image_<q>, inclusive): parabola 3E q = -0.10..+0.45, 1A1 q = +0.60..+1.15; cubic a Q^2 + b Q^3 + c (perturbative correction) q = -0.10..+1.15 for both; every other fit (splines, Morse, polyN) ALL images
- all images: 3E all 33 sc images (q = -1.00..+2.00), 1A1 all 31 sf images (q = -1.00..+2.00)

**FAIL; Finite difference [cubic smoothing spline (rms 0.05 meV)] anharmonic, so differences are physics, not error** (tolerances: max dF_ij <= 1e-07, max dVOF / peak <= 1e-08)

| direction | test | reference | max dF_ij | at (i, j) | max dF_ij (i,j<30) | max dVOF / peak | rms dVOF / peak | F_ij time (s) | verdict |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| emission | Finite difference [cubic smoothing spline (rms 0.05 meV)] | Hermite (quadrature) | 7.76e-01 | (0, 18) | 7.76e-01 | 5.09e+00 | 3.76e-01 | 28.17 | FAIL |
| absorption | Finite difference [cubic smoothing spline (rms 0.05 meV)] | Hermite (quadrature) | 7.76e-01 | (18, 0) | 7.76e-01 | 9.58e+00 | 8.63e-01 | 15.45 | FAIL |

## Files

| File | Content |
| --- | --- |
| `log.txt` | Full console output, including every method's settings. |
| `settings.json` | Every parameter of the run. |
| **`data/`** | Numbers behind the figures, as CSV. |
| `data/summary.csv` | Every error statistic per direction and test, with OK / FAIL. |
| `data/vof_absorption.csv` | absorption (1A1 -> 3E): every lineshape and every difference from the reference. |
| `data/vof_emission.csv` | emission (3E -> 1A1): every lineshape and every difference from the reference. |
| **`figures/`** | Comparison figures. |
| `figures/fc_error_maps.png` | log10 \|F_ij(test) - F_ij(reference)\| per test (rows) and direction (columns); the cross marks the maximum. |
| `figures/vof_comparison.png` | Per direction (columns): the VOFs overlaid (reference thick and translucent, tests dashed / dotted), the signed difference (test - reference) / reference peak, and its absolute value on a log scale with the tolerance --tol-vof (dashed). Verdict: FAIL; Finite difference [cubic smoothing spline (rms 0.05 meV)] anharmonic, so differences are physics, not error. Every setting of every method is in log.txt. |
