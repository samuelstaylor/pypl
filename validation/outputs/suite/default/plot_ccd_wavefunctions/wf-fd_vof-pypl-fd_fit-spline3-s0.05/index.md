# 1D CCD with wavefunctions: wf-fd_vof-pypl-fd_fit-spline3-s0.05

Generated 2026-10-15T18:39:44 by `python C:\Users\sammy\Research\Galli\programs\pypl\testing\plot_ccd_wavefunctions.py --wf-method fd --vof-method pypl fd --temp 5.0 --sigma 8.0 25.0 --gamma 0.3 --order-es 50 --order-gs 60 --ene-range -200.0 1200.0 --resol 1401 --npts 100001 --fd-npts 20001 --fd-stencil 5 --fd-margin 2.0 --fd-fit spline3-s0.05 --n-levels 3`.

Wavefunctions: finite difference on cubic smoothing spline (rms 0.05 meV). VOF methods: pypl, fd.

Fit **spline3-s0.05** (cubic smoothing spline (rms 0.05 meV)):
- GS well: cubic spline (degree 3, smoothing to rms <= 5e-05 eV, 3 knots), harmonic tails a = 0.417627 eV/(amu A^2), C0 junction (slope jumps -2.02e-03 / -1.35e-02 eV/(amu^0.5 A)); fitted to samples 1-15 of 15 (Q = -0.139..0.835); rms residual 5.00e-05 eV
- ES well: cubic spline (degree 3, smoothing to rms <= 5e-05 eV, 2 knots), harmonic tails a = 0.494721 eV/(amu A^2), C0 junction (slope jumps 1.23e-02 / -3.12e-03 eV/(amu^0.5 A)); fitted to samples 1-11 of 11 (Q = 0.139..0.835); rms residual 4.41e-05 eV

Wavefunction check: max |<chi_i|chi_j> - F_ij| = 1.00e-02 (tolerance 1e-08); mismatch expected, the fit is anharmonic.

## Files

| File | Content |
| --- | --- |
| `log.txt` | Full console output, including every method's settings. |
| `settings.json` | Every parameter of the run. |
| `ccd_wavefunctions.png` | The CCD with 3 wavefunctions per well (finite difference on cubic smoothing spline (rms 0.05 meV)) and the VOF panel. |
