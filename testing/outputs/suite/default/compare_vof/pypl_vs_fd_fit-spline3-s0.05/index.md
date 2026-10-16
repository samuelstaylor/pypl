# VOF method comparison: pypl_vs_fd_fit-spline3-s0.05

Generated 2026-10-15T18:31:44 by `python C:\Users\sammy\Research\Galli\programs\pypl\testing\compare_vof.py --reference pypl --test fd --temp 5.0 --sigma 8.0 25.0 --gamma 0.3 --order-es 50 --order-gs 60 --ene-range -200.0 1200.0 --resol 1401 --npts 100001 --fd-npts 20001 --fd-stencil 5 --fd-margin 2.0 --fd-fit spline3-s0.05 --tol-fc 1e-07 --tol-vof 1e-08`.

**FAIL; Finite difference [cubic smoothing spline (rms 0.05 meV)] anharmonic, so differences are physics, not error** (tolerances: max dF_ij <= 1e-07, max dVOF / peak <= 1e-08)

| direction | test | reference | max dF_ij | at (i, j) | max dF_ij (i,j<30) | max dVOF / peak | rms dVOF / peak | F_ij time (s) | verdict |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| emission | Finite difference [cubic smoothing spline (rms 0.05 meV)] | PyPL (Ruhoff) | 1.81e-02 | (5, 8) | 1.81e-02 | 2.21e-02 | 9.39e-04 | 36.21 | FAIL |
| absorption | Finite difference [cubic smoothing spline (rms 0.05 meV)] | PyPL (Ruhoff) | 1.81e-02 | (8, 5) | 1.81e-02 | 2.21e-02 | 8.08e-04 | 43.98 | FAIL |

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
