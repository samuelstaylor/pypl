# 1D CCD with wavefunctions: wf-fd_vof-pypl-fd_fit-spline3-notaknot-smooth

Generated 2026-10-16T21:08:42 by `python C:\Users\sammy\Research\Galli\programs\pypl\validation\plot_ccd_wavefunctions.py --wf-method fd --vof-method pypl fd --dataset nv_diamond --temp 5.0 --sigma 8.0 25.0 --gamma 0.3 --order-es 50 --order-gs 60 --ene-range -200.0 1200.0 --resol 1401 --line-positions harmonic --npts 100001 --hermite-margin 6.0 --fd-stencil 5 --fd-margin 2.0 --fd-fit spline3-notaknot-smooth --fd-npts 20001 --cubic-warn 0.1 --dpi 300 --fig-format png --font-size 9.0 --n-levels 3 --check-tol 1e-10 1e-08`.

**Dataset.**

- dataset: nv_diamond (NV- in diamond), 3E -> 3A2
- data: doc/tutorials/002_nv_diamond_1d_ccd
- roles: 'ES' = ES (upper state, minimum at Q = delta_q), 'GS' = GS (lower state, minimum at Q = 0); emission = ES -> GS = ES -> GS
- energies: ES = ΔSCF, constrained occupations (PBE), 3E; GS = DFT (PBE), 3A2; zero = the GS (3A2) minimum image

Wavefunctions: FD wavefunctions on the cubic spline (not-a-knot, smooth C2 tails). VOF methods: pypl, fd.

Fit **spline3-notaknot-smooth** (cubic spline (not-a-knot, smooth C2 tails)):
- GS well: cubic spline (degree 3, not-a-knot ends, 13 knots), smooth C2 Taylor tails from the spline's own ends (V' = -0.116 / +0.713 eV/(amu^0.5 A), V'' = +0.859 / +0.898 eV/(amu A^2) at Q = -0.139 / 0.835); fitted to samples 1-15 of 15 (Q = -0.139..0.835); rms residual 8.78e-18 eV
- ES well: cubic spline (degree 3, not-a-knot ends, 9 knots), smooth C2 Taylor tails from the spline's own ends (V' = -0.560 / +0.144 eV/(amu^0.5 A), V'' = +0.968 / +1.071 eV/(amu A^2) at Q = 0.139 / 0.835); fitted to samples 1-11 of 11 (Q = 0.139..0.835); rms residual 1.16e-16 eV

Wavefunction check: max |<chi_i|chi_j> - F_ij| = 3.75e-03 (tolerance 1e-08); mismatch expected, the fit is anharmonic.

## Files

| File | Content |
| --- | --- |
| `log.txt` | Full console output, including every method's settings. |
| `settings.json` | Every parameter of the run. |
| `ccd_wavefunctions.png` | The CCD with 3 wavefunctions per well (FD wavefunctions on the cubic spline (not-a-knot, smooth C2 tails)) and the VOF panel. |
