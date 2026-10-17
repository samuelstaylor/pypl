# 1D CCD with wavefunctions: wf-fd_vof-pypl-fd_fit-poly2

Generated 2026-10-16T21:08:29 by `python C:\Users\sammy\Research\Galli\programs\pypl\validation\plot_ccd_wavefunctions.py --wf-method fd --vof-method pypl fd --dataset nv_diamond --temp 5.0 --sigma 8.0 25.0 --gamma 0.3 --order-es 50 --order-gs 60 --ene-range -200.0 1200.0 --resol 1401 --line-positions harmonic --npts 100001 --hermite-margin 6.0 --fd-stencil 5 --fd-margin 2.0 --fd-fit poly2 --fd-npts 20001 --cubic-warn 0.1 --dpi 300 --fig-format png --font-size 9.0 --n-levels 3 --check-tol 1e-10 1e-08`.

**Dataset.**

- dataset: nv_diamond (NV- in diamond), 3E -> 3A2
- data: doc/tutorials/002_nv_diamond_1d_ccd
- roles: 'ES' = ES (upper state, minimum at Q = delta_q), 'GS' = GS (lower state, minimum at Q = 0); emission = ES -> GS = ES -> GS
- energies: ES = ΔSCF, constrained occupations (PBE), 3E; GS = DFT (PBE), 3A2; zero = the GS (3A2) minimum image

Wavefunctions: FD wavefunctions on the polynomial deg. 2. VOF methods: pypl, fd.

Fit **poly2** (polynomial deg. 2):
- GS well: degree-2 polynomial in (Q - 0.000000), coefficients (highest first) [0.4264, -0.00111, -0.0001572]; fitted to samples 1-15 of 15 (Q = -0.139..0.835); rms residual 1.52e-04 eV
- ES well: degree-2 polynomial in (Q - 0.695732), coefficients (highest first) [0.5058, 0.002201, 1.706]; fitted to samples 1-11 of 11 (Q = 0.139..0.835); rms residual 7.66e-05 eV

Wavefunction check: max |<chi_i|chi_j> - F_ij| = 6.01e-04 (tolerance 1e-08); mismatch expected, the fit is anharmonic.

## Files

| File | Content |
| --- | --- |
| `log.txt` | Full console output, including every method's settings. |
| `settings.json` | Every parameter of the run. |
| `ccd_wavefunctions.png` | The CCD with 3 wavefunctions per well (FD wavefunctions on the polynomial deg. 2) and the VOF panel. |
