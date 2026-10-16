# 1D CCD with wavefunctions: wf-fd_vof-pypl-fd_fit-spline

Generated 2026-10-15T16:07:28 by `python plot_ccd_wavefunctions.py --wf-method fd --vof-method pypl fd --fd-fit spline`.

Wavefunctions: finite difference on cubic spline. VOF methods: pypl, fd.

Fit **spline** (cubic spline):
- GS well: natural cubic spline, harmonic tails a = 0.417627 eV/(amu A^2) (offsets -0.000122 / 0.005065 eV); fitted to samples 1-15 of 15 (Q = -0.139..0.835); rms residual 0.00e+00 eV
- ES well: natural cubic spline, harmonic tails a = 0.494721 eV/(amu A^2) (offsets 1.708018 / 1.706203 eV); fitted to samples 1-11 of 11 (Q = 0.139..0.835); rms residual 6.69e-17 eV

Wavefunction check: max |<chi_i|chi_j> - F_ij| = 9.92e-03 (tolerance 1e-08); mismatch expected, the fit is anharmonic.

## Files

| File | Content |
| --- | --- |
| `log.txt` | Full console output, including every method's settings. |
| `settings.json` | Every parameter of the run. |
| `ccd_wavefunctions.png` | The CCD with 3 wavefunctions per well (finite difference on cubic spline) and the VOF panel. |
