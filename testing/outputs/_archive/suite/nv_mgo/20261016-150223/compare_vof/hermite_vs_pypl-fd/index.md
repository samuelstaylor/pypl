# VOF method comparison: hermite_vs_pypl-fd [NV- in MgO]

Generated 2026-10-16T15:00:14 by `python C:\Users\sammy\Research\Galli\programs\pypl\testing\compare_vof.py --reference hermite --test pypl fd --dataset nv_mgo --temp 5.0 --sigma 20.0 --gamma 0.3 --no-zpl --order-es 160 --order-gs 100 --ene-range -500.0 10000.0 --resol 2101 --npts 100001 --fd-npts 40001 --fd-stencil 5 --fd-margin 2.0 --fd-fit parabola --tol-fc 1e-07 --tol-vof 1e-08`.

**Dataset.**

- dataset: nv_mgo (NV- in MgO), 3E -> 1A1 (WEST TDDFT-DDH)
- data: doc/nv_mgo_west
- roles: 'ES' = 3E (upper state, minimum at Q = 0), 'GS' = 1A1 (lower state, minimum at Q = delta_q); emission = ES -> GS = 3E -> 1A1
- Q = q delta_q, q from Image_<q>, as in reference/ccd_nv_mgo.py: 3E minimum at Q = 0, 1A1 minimum at Q = delta_q
- E = QE total energy (qe/) + WEST root (3E: sc root 1, 1A1: sf root 4) x 13.605693 eV/Ry, relative to the 3A2 (QE) minimum, as in ccd_nv_mgo.py
- every fit except the parabola uses ALL images of each PES: 3E all 33 sc images (q = -1.00..+2.00), 1A1 all 30 sf images except Image_+2.00 (not a spin-flip run, so not 1A1)
- parabolas only: the five images nearest each minimum, 3E q = -0.10..+0.10, 1A1 q = +0.90..+1.10

**FAIL** (tolerances: max dF_ij <= 1e-07, max dVOF / peak <= 1e-08)

| direction | test | reference | max dF_ij | at (i, j) | max dF_ij (i,j<30) | max dVOF / peak | rms dVOF / peak | F_ij time (s) | verdict |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| emission | PyPL (Ruhoff) | Hermite (quadrature) | 1.68e+31 | (159, 99) | 4.21e-06 | 8.25e-15 | 1.64e-15 | 0.04 | FAIL |
| emission | Finite difference | Hermite (quadrature) | 1.29e-08 | (159, 96) | 1.96e-10 | 5.58e-10 | 1.26e-10 | 11.10 | OK |
| absorption | PyPL (Ruhoff) | Hermite (quadrature) | 2.63e+31 | (99, 159) | 1.50e-07 | 5.14e-15 | 1.19e-15 | 0.03 | FAIL |
| absorption | Finite difference | Hermite (quadrature) | 1.29e-08 | (96, 159) | 1.96e-10 | 3.64e-09 | 1.11e-09 | 14.61 | OK |

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
| `figures/fc_error_maps.png` | log10 \|F_ij(test) - F_ij(reference)\| per test and direction. |
| `figures/vof_comparison.png` | Per direction: VOFs overlaid, (test - ref) / peak, log10 \|test - ref\| / peak with the tolerance, and every method's settings. |
