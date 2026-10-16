# 1D CCD with wavefunctions: wf-fd_vof-pypl-fd_fit-spline4-notaknot

Generated 2026-10-15T18:39:35 by `python C:\Users\sammy\Research\Galli\programs\pypl\testing\plot_ccd_wavefunctions.py --wf-method fd --vof-method pypl fd --temp 5.0 --sigma 8.0 25.0 --gamma 0.3 --order-es 50 --order-gs 60 --ene-range -200.0 1200.0 --resol 1401 --npts 100001 --fd-npts 20001 --fd-stencil 5 --fd-margin 2.0 --fd-fit spline4-notaknot --n-levels 3`.

Wavefunctions: finite difference on quartic spline (not-a-knot). VOF methods: pypl, fd.

Fit **spline4-notaknot** (quartic spline (not-a-knot)):
- GS well: quartic spline (degree 4, not-a-knot ends, 12 knots), harmonic tails a = 0.417627 eV/(amu A^2), C0 junction (slope jumps 6.94e-04 / -1.52e-02 eV/(amu^0.5 A)); fitted to samples 1-15 of 15 (Q = -0.139..0.835); rms residual 1.65e-17 eV
- ES well: quartic spline (degree 4, not-a-knot ends, 8 knots), harmonic tails a = 0.494721 eV/(amu A^2), C0 junction (slope jumps 8.74e-03 / -6.52e-03 eV/(amu^0.5 A)); fitted to samples 1-11 of 11 (Q = 0.139..0.835); rms residual 4.97e-16 eV

Wavefunction check: max |<chi_i|chi_j> - F_ij| = 9.87e-03 (tolerance 1e-08); mismatch expected, the fit is anharmonic.

## Files

| File | Content |
| --- | --- |
| `log.txt` | Full console output, including every method's settings. |
| `settings.json` | Every parameter of the run. |
| `ccd_wavefunctions.png` | The CCD with 3 wavefunctions per well (finite difference on quartic spline (not-a-knot)) and the VOF panel. |
