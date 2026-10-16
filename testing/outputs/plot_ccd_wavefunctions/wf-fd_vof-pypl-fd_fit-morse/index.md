# 1D CCD with wavefunctions: wf-fd_vof-pypl-fd_fit-morse

Generated 2026-10-15T16:07:47 by `python plot_ccd_wavefunctions.py --wf-method fd --vof-method pypl fd --fd-fit morse`.

Wavefunctions: finite difference on Morse. VOF methods: pypl, fd.

Fit **morse** (Morse):
- GS well: Morse D (1 - exp(-b (Q - Q_m)))^2 + c, D = 2255.18 eV, b = -0.0136521, Q_m = -0.000070, c = -0.000075 eV; fitted to samples 1-15 of 15 (Q = -0.139..0.835); rms residual 1.10e-04 eV
- ES well: Morse D (1 - exp(-b (Q - Q_m)))^2 + c, D = 2074.73 eV, b = 0.0155388, Q_m = 0.693884, c = 1.705824 eV; fitted to samples 1-11 of 11 (Q = 0.139..0.835); rms residual 4.41e-05 eV

Wavefunction check: max |<chi_i|chi_j> - F_ij| = 9.42e-03 (tolerance 1e-08); mismatch expected, the fit is anharmonic.

## Files

| File | Content |
| --- | --- |
| `log.txt` | Full console output, including every method's settings. |
| `settings.json` | Every parameter of the run. |
| `ccd_wavefunctions.png` | The CCD with 3 wavefunctions per well (finite difference on Morse) and the VOF panel. |
