# 1D CCD with wavefunctions: wf-fd_vof-pypl-hermite-fd_fit-parabola

Generated 2026-10-15T16:08:01 by `python plot_ccd_wavefunctions.py --wf-method fd --vof-method pypl hermite fd`.

Wavefunctions: finite difference on parabola (near min.). VOF methods: pypl, hermite, fd.

Fit **parabola** (parabola (near min.)):
- GS well: parabola a (Q - 0.000000)^2 + c, a = 0.417627 eV/(amu A^2), c = -0.000013 eV; fitted to samples 1-5 of 15 (Q = -0.139..0.139); rms residual 7.49e-05 eV
- ES well: parabola a (Q - 0.695732)^2 + c, a = 0.494721 eV/(amu A^2), c = 1.705877 eV; fitted to samples 5-11 of 11 (Q = 0.417..0.835); rms residual 1.49e-04 eV

Wavefunction check: max |<chi_i|chi_j> - F_ij| = 1.63e-12 (tolerance 1e-08).

## Files

| File | Content |
| --- | --- |
| `log.txt` | Full console output, including every method's settings. |
| `settings.json` | Every parameter of the run. |
| `ccd_wavefunctions.png` | The CCD with 3 wavefunctions per well (finite difference on parabola (near min.)) and the VOF panel. |
