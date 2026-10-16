# 1D CCD with wavefunctions: wf-fd_vof-pypl-fd_fit-cubic

Generated 2026-10-15T18:38:42 by `python C:\Users\sammy\Research\Galli\programs\pypl\testing\plot_ccd_wavefunctions.py --wf-method fd --vof-method pypl fd --temp 5.0 --sigma 8.0 25.0 --gamma 0.3 --order-es 50 --order-gs 60 --ene-range -200.0 1200.0 --resol 1401 --npts 100001 --fd-npts 20001 --fd-stencil 5 --fd-margin 2.0 --fd-fit cubic --n-levels 3`.

Wavefunctions: finite difference on cubic (a fixed, all images). VOF methods: pypl, fd.

Fit **cubic** (cubic (a fixed, all images)):
- GS well: cubic a (Q - 0.000000)^2 + b (Q - 0.000000)^3 + c, a = 0.417627 eV/(amu A^2) (fixed), b = +0.009164 eV/(amu^1.5 A^3), c = 0.000052 eV, turnover at Q = -30.38; fitted to samples 1-15 of 15 (Q = -0.139..0.835); rms residual 1.60e-04 eV; unbounded beyond the turnover at Q = -30.38 (perturbative or in-box use only)
- ES well: cubic a (Q - 0.695732)^2 + b (Q - 0.695732)^3 + c, a = 0.494721 eV/(amu A^2) (fixed), b = -0.012655 eV/(amu^1.5 A^3), c = 1.705818 eV, turnover at Q = 26.76; fitted to samples 1-11 of 11 (Q = 0.139..0.835); rms residual 1.63e-04 eV; unbounded beyond the turnover at Q = 26.76 (perturbative or in-box use only)

Wavefunction check: max |<chi_i|chi_j> - F_ij| = 1.58e-02 (tolerance 1e-08); mismatch expected, the fit is anharmonic.

## Files

| File | Content |
| --- | --- |
| `log.txt` | Full console output, including every method's settings. |
| `settings.json` | Every parameter of the run. |
| `ccd_wavefunctions.png` | The CCD with 3 wavefunctions per well (finite difference on cubic (a fixed, all images)) and the VOF panel. |
