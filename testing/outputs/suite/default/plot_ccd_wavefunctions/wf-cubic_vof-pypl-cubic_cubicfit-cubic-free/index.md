# 1D CCD with wavefunctions: wf-cubic_vof-pypl-cubic_cubicfit-cubic-free

Generated 2026-10-15T18:40:12 by `python C:\Users\sammy\Research\Galli\programs\pypl\testing\plot_ccd_wavefunctions.py --wf-method cubic --vof-method pypl cubic --temp 5.0 --sigma 8.0 25.0 --gamma 0.3 --order-es 50 --order-gs 60 --ene-range -200.0 1200.0 --resol 1401 --npts 100001 --fd-npts 20001 --fd-stencil 5 --fd-margin 2.0 --fd-fit parabola --cubic-points all --cubic-free-a --n-levels 3`.

Wavefunctions: Hermite + cubic PT, cubic (a fitted, all images). VOF methods: pypl, cubic.

Wavefunction check: max |<chi_i|chi_j> - F_ij| = 5.00e-16 (tolerance 1e-10).

## Files

| File | Content |
| --- | --- |
| `log.txt` | Full console output, including every method's settings. |
| `settings.json` | Every parameter of the run. |
| `ccd_wavefunctions.png` | The CCD with 3 wavefunctions per well (Hermite + cubic PT, cubic (a fitted, all images)) and the VOF panel. |
