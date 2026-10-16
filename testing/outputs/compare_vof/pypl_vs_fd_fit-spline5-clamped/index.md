# VOF method comparison: pypl_vs_fd_fit-spline5-clamped

Generated 2026-10-15T16:25:21 by `python compare_vof.py --test fd --fd-fit spline --spline-order 5 --spline-bc clamped`.

**FAIL; fd uses fit spline5-clamped, so differences are anharmonicity, not error** (tolerances: max dF_ij <= 1e-07, max dVOF / peak <= 1e-08)

| direction | test | reference | max dF_ij | at (i, j) | max dF_ij (i,j<30) | max dVOF / peak | rms dVOF / peak | F_ij time (s) | verdict |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| emission | Finite difference [quintic spline (clamped)] | PyPL (Ruhoff) | 1.79e-02 | (5, 8) | 1.79e-02 | 2.16e-02 | 9.36e-04 | 1.60 | FAIL |
| absorption | Finite difference [quintic spline (clamped)] | PyPL (Ruhoff) | 1.79e-02 | (8, 5) | 1.79e-02 | 2.16e-02 | 8.02e-04 | 7.54 | FAIL |

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
