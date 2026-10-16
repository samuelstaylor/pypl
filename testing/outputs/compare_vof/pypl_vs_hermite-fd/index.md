# VOF method comparison: pypl_vs_hermite-fd

Generated 2026-10-15T16:54:25 by `python compare_vof.py`.

**PASS** (tolerances: max dF_ij <= 1e-07, max dVOF / peak <= 1e-08)

| direction | test | reference | max dF_ij | at (i, j) | max dF_ij (i,j<30) | max dVOF / peak | rms dVOF / peak | F_ij time (s) | verdict |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| emission | Hermite (quadrature) | PyPL (Ruhoff) | 1.70e-08 | (49, 48) | 4.87e-11 | 6.95e-16 | 4.51e-17 | 1.33 | OK |
| emission | Finite difference | PyPL (Ruhoff) | 1.80e-08 | (49, 49) | 6.11e-10 | 1.29e-12 | 3.81e-13 | 2.72 | OK |
| absorption | Hermite (quadrature) | PyPL (Ruhoff) | 1.37e-08 | (50, 49) | 3.71e-11 | 6.95e-16 | 3.83e-17 | 1.61 | OK |
| absorption | Finite difference | PyPL (Ruhoff) | 1.51e-08 | (50, 49) | 6.11e-10 | 9.55e-13 | 3.37e-13 | 4.86 | OK |

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
| `figures/fc_error_maps.png` | log10 \|F_ij(test) - F_ij(reference)\| per test and direction. |
| `figures/vof_comparison.png` | Per direction: VOFs overlaid, (test - ref) / peak, log10 \|test - ref\| / peak with the tolerance, and every method's settings. |
