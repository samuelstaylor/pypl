# VOF method comparison: fd_vs_cubic_fit-cubic_cubicfit-cubic

Generated 2026-10-15T17:40:03 by `python compare_vof.py --reference fd --fd-fit cubic --test cubic`.

**FAIL; same cubic potential: differences are the error of first-order perturbation theory** (tolerances: max dF_ij <= 1e-07, max dVOF / peak <= 1e-08)

| direction | test | reference | max dF_ij | at (i, j) | max dF_ij (i,j<30) | max dVOF / peak | rms dVOF / peak | F_ij time (s) | verdict |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| emission | Hermite + cubic PT | Finite difference [cubic (a fixed, all images)] | 1.16e+00 | (48, 59) | 1.62e-01 | 3.25e-04 | 3.20e-05 | 1.50 | FAIL |
| absorption | Hermite + cubic PT | Finite difference [cubic (a fixed, all images)] | 1.16e+00 | (59, 48) | 1.62e-01 | 3.25e-04 | 3.59e-05 | 1.52 | FAIL |

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
