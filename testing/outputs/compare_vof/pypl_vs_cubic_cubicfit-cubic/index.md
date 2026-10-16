# VOF method comparison: pypl_vs_cubic_cubicfit-cubic

Generated 2026-10-15T17:39:47 by `python compare_vof.py --test cubic`.

**FAIL; Hermite + cubic PT anharmonic, so differences are physics, not error** (tolerances: max dF_ij <= 1e-07, max dVOF / peak <= 1e-08)

| direction | test | reference | max dF_ij | at (i, j) | max dF_ij (i,j<30) | max dVOF / peak | rms dVOF / peak | F_ij time (s) | verdict |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| emission | Hermite + cubic PT | PyPL (Ruhoff) | 1.18e+00 | (43, 59) | 3.48e-01 | 3.41e-02 | 1.54e-03 | 1.62 | FAIL |
| absorption | Hermite + cubic PT | PyPL (Ruhoff) | 1.18e+00 | (59, 43) | 3.48e-01 | 3.41e-02 | 1.54e-03 | 1.61 | FAIL |

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
