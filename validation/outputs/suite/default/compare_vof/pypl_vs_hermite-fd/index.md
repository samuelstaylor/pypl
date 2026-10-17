# VOF method comparison: pypl_vs_hermite-fd

Generated 2026-10-16T21:06:06 by `python C:\Users\sammy\Research\Galli\programs\pypl\validation\compare_vof.py --reference pypl --test hermite fd --dataset nv_diamond --temp 5.0 --sigma 8.0 25.0 --gamma 0.3 --order-es 50 --order-gs 60 --ene-range -200.0 1200.0 --resol 1401 --line-positions harmonic --npts 100001 --hermite-margin 6.0 --fd-stencil 5 --fd-margin 2.0 --fd-fit parabola --fd-npts 20001 --cubic-warn 0.1 --dpi 300 --fig-format png --font-size 9.0 --tol-fc 1e-07 --tol-vof 1e-08`.

**Dataset.**

- dataset: nv_diamond (NV- in diamond), 3E -> 3A2
- data: doc/tutorials/002_nv_diamond_1d_ccd
- roles: 'ES' = ES (upper state, minimum at Q = delta_q), 'GS' = GS (lower state, minimum at Q = 0); emission = ES -> GS = ES -> GS
- energies: ES = ΔSCF, constrained occupations (PBE), 3E; GS = DFT (PBE), 3A2; zero = the GS (3A2) minimum image

**PASS** (tolerances: max dF_ij <= 1e-07, max dVOF / peak <= 1e-08)

| direction | test | reference | max dF_ij | at (i, j) | max dF_ij (i,j<30) | max dVOF / peak | rms dVOF / peak | F_ij time (s) | verdict |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| emission | Hermite (quadrature) | PyPL (Ruhoff) | 1.70e-08 | (49, 48) | 4.87e-11 | 6.95e-16 | 3.97e-17 | 0.42 | OK |
| emission | Finite difference | PyPL (Ruhoff) | 1.80e-08 | (49, 49) | 6.11e-10 | 1.29e-12 | 3.81e-13 | 3.89 | OK |
| absorption | Hermite (quadrature) | PyPL (Ruhoff) | 1.37e-08 | (50, 49) | 3.71e-11 | 6.95e-16 | 4.46e-17 | 0.63 | OK |
| absorption | Finite difference | PyPL (Ruhoff) | 1.51e-08 | (50, 49) | 6.11e-10 | 9.55e-13 | 3.37e-13 | 0.87 | OK |

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
| `figures/vof_comparison.png` | Per direction (columns): the VOFs overlaid (reference thick and translucent, tests dashed / dotted), the signed difference (test - reference) / reference peak, and its absolute value on a log scale with the tolerance --tol-vof (dashed). Verdict: PASS. Every setting of every method is in log.txt. |
