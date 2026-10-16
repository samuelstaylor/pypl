# 1D CCD with wavefunctions: wf-fd_vof-pypl-fd_fit-poly6

Generated 2026-10-15T18:37:41 by `python C:\Users\sammy\Research\Galli\programs\pypl\testing\plot_ccd_wavefunctions.py --wf-method fd --vof-method pypl fd --temp 5.0 --sigma 8.0 25.0 --gamma 0.3 --order-es 50 --order-gs 60 --ene-range -200.0 1200.0 --resol 1401 --npts 100001 --fd-npts 20001 --fd-stencil 5 --fd-margin 2.0 --fd-fit poly6 --n-levels 3`.

Wavefunctions: finite difference on polynomial deg. 6. VOF methods: pypl, fd.

Fit **poly6** (polynomial deg. 6):
- GS well: degree-6 polynomial in (Q - 0.000000), coefficients (highest first) [0.1337, -0.2934, 0.2023, -0.0296, 0.4127, 0.001382, 5.824e-06]; fitted to samples 1-15 of 15 (Q = -0.139..0.835); rms residual 2.37e-05 eV
- ES well: degree-6 polynomial in (Q - 0.695732), coefficients (highest first) [0.2858, 0.5198, 0.2939, 0.03818, 0.4955, 0.0007426, 1.706]; fitted to samples 1-11 of 11 (Q = 0.139..0.835); rms residual 6.13e-06 eV

Wavefunction check: max |<chi_i|chi_j> - F_ij| = 8.80e-02 (tolerance 1e-08); mismatch expected, the fit is anharmonic.

## Files

| File | Content |
| --- | --- |
| `log.txt` | Full console output, including every method's settings. |
| `settings.json` | Every parameter of the run. |
| `ccd_wavefunctions.png` | The CCD with 3 wavefunctions per well (finite difference on polynomial deg. 6) and the VOF panel. |
