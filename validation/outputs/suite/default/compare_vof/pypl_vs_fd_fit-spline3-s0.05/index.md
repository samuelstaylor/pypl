# VOF method comparison: pypl_vs_fd_fit-spline3-s0.05

Generated 2026-10-16T21:06:54 by `python C:\Users\sammy\Research\Galli\programs\pypl\validation\compare_vof.py --reference pypl --test fd --dataset nv_diamond --temp 5.0 --sigma 8.0 25.0 --gamma 0.3 --order-es 50 --order-gs 60 --ene-range -200.0 1200.0 --resol 1401 --line-positions harmonic --npts 100001 --hermite-margin 6.0 --fd-stencil 5 --fd-margin 2.0 --fd-fit spline3-s0.05 --fd-npts 20001 --cubic-warn 0.1 --dpi 300 --fig-format png --font-size 9.0 --tol-fc 1e-07 --tol-vof 1e-08`.

**Dataset.**

- dataset: nv_diamond (NV- in diamond), 3E -> 3A2
- data: doc/tutorials/002_nv_diamond_1d_ccd
- roles: 'ES' = ES (upper state, minimum at Q = delta_q), 'GS' = GS (lower state, minimum at Q = 0); emission = ES -> GS = ES -> GS
- energies: ES = ΔSCF, constrained occupations (PBE), 3E; GS = DFT (PBE), 3A2; zero = the GS (3A2) minimum image

**FAIL; Finite difference [cubic smoothing spline (rms 0.05 meV)] anharmonic, so differences are physics, not error** (tolerances: max dF_ij <= 1e-07, max dVOF / peak <= 1e-08)

| direction | test | reference | max dF_ij | at (i, j) | max dF_ij (i,j<30) | max dVOF / peak | rms dVOF / peak | F_ij time (s) | verdict |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| emission | Finite difference [cubic smoothing spline (rms 0.05 meV)] | PyPL (Ruhoff) | 1.81e-02 | (5, 8) | 1.81e-02 | 2.21e-02 | 9.39e-04 | 0.92 | FAIL |
| absorption | Finite difference [cubic smoothing spline (rms 0.05 meV)] | PyPL (Ruhoff) | 1.81e-02 | (8, 5) | 1.81e-02 | 2.21e-02 | 8.08e-04 | 0.98 | FAIL |

## Files

| File | Content |
| --- | --- |
| `log.txt` | Full console output, including every method's settings. |
| `settings.json` | Every parameter of the run. |
| **`data/`** | Numbers behind the figures, as CSV. |
| `data/summary.csv` | Every error statistic per direction and test, with OK / FAIL. |
| `data/vof_absorption.csv` | absorption (GS -> ES): every lineshape and every difference from the reference. |
| `data/vof_emission.csv` | emission (ES -> GS): every lineshape and every difference from the reference. |
| **`figures/`** | Comparison figures. |
| `figures/fc_error_maps.png` | log10 \|F_ij(test) - F_ij(reference)\| per test (rows) and direction (columns); the cross marks the maximum. |
| `figures/vof_comparison.png` | Per direction (columns): the VOFs overlaid (reference thick and translucent, tests dashed / dotted), the signed difference (test - reference) / reference peak, and its absolute value on a log scale with the tolerance --tol-vof (dashed). Verdict: FAIL; Finite difference [cubic smoothing spline (rms 0.05 meV)] anharmonic, so differences are physics, not error. Every setting of every method is in log.txt. |
