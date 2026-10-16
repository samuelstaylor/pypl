# 1D CCD with wavefunctions: wf-fd_vof-pypl-fd_fit-poly2

Generated 2026-10-15T18:37:31 by `python C:\Users\sammy\Research\Galli\programs\pypl\testing\plot_ccd_wavefunctions.py --wf-method fd --vof-method pypl fd --temp 5.0 --sigma 8.0 25.0 --gamma 0.3 --order-es 50 --order-gs 60 --ene-range -200.0 1200.0 --resol 1401 --npts 100001 --fd-npts 20001 --fd-stencil 5 --fd-margin 2.0 --fd-fit poly2 --n-levels 3`.

Wavefunctions: finite difference on polynomial deg. 2. VOF methods: pypl, fd.

Fit **poly2** (polynomial deg. 2):
- GS well: degree-2 polynomial in (Q - 0.000000), coefficients (highest first) [0.4264, -0.00111, -0.0001572]; fitted to samples 1-15 of 15 (Q = -0.139..0.835); rms residual 1.52e-04 eV
- ES well: degree-2 polynomial in (Q - 0.695732), coefficients (highest first) [0.5058, 0.002201, 1.706]; fitted to samples 1-11 of 11 (Q = 0.139..0.835); rms residual 7.66e-05 eV

Wavefunction check: max |<chi_i|chi_j> - F_ij| = 6.01e-04 (tolerance 1e-08); mismatch expected, the fit is anharmonic.

## Files

| File | Content |
| --- | --- |
| `log.txt` | Full console output, including every method's settings. |
| `settings.json` | Every parameter of the run. |
| `ccd_wavefunctions.png` | The CCD with 3 wavefunctions per well (finite difference on polynomial deg. 2) and the VOF panel. |
