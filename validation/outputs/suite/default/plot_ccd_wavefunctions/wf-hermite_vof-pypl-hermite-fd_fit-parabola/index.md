# 1D CCD with wavefunctions: wf-hermite_vof-pypl-hermite-fd_fit-parabola

Generated 2026-10-15T18:35:44 by `python C:\Users\sammy\Research\Galli\programs\pypl\testing\plot_ccd_wavefunctions.py --wf-method hermite --vof-method pypl hermite fd --temp 5.0 --sigma 8.0 25.0 --gamma 0.3 --order-es 50 --order-gs 60 --ene-range -200.0 1200.0 --resol 1401 --npts 100001 --fd-npts 20001 --fd-stencil 5 --fd-margin 2.0 --fd-fit parabola --n-levels 3`.

Wavefunctions: Hermite. VOF methods: pypl, hermite, fd.

Fit **parabola** (parabola (near min.)):
- GS well: parabola a (Q - 0.000000)^2 + c, a = 0.417627 eV/(amu A^2), c = -0.000013 eV; fitted to samples 1-5 of 15 (Q = -0.139..0.139); rms residual 7.49e-05 eV
- ES well: parabola a (Q - 0.695732)^2 + c, a = 0.494721 eV/(amu A^2), c = 1.705877 eV; fitted to samples 5-11 of 11 (Q = 0.417..0.835); rms residual 1.49e-04 eV

Wavefunction check: max |<chi_i|chi_j> - F_ij| = 1.67e-16 (tolerance 1e-10).

## Files

| File | Content |
| --- | --- |
| `log.txt` | Full console output, including every method's settings. |
| `settings.json` | Every parameter of the run. |
| `ccd_wavefunctions.png` | The CCD with 3 wavefunctions per well (Hermite) and the VOF panel. |
