# 1D CCD with wavefunctions: wf-fd_vof-pypl-fd_fit-poly6

Generated 2026-10-16T21:08:34 by `python C:\Users\sammy\Research\Galli\programs\pypl\validation\plot_ccd_wavefunctions.py --wf-method fd --vof-method pypl fd --dataset nv_diamond --temp 5.0 --sigma 8.0 25.0 --gamma 0.3 --order-es 50 --order-gs 60 --ene-range -200.0 1200.0 --resol 1401 --line-positions harmonic --npts 100001 --hermite-margin 6.0 --fd-stencil 5 --fd-margin 2.0 --fd-fit poly6 --fd-npts 20001 --cubic-warn 0.1 --dpi 300 --fig-format png --font-size 9.0 --n-levels 3 --check-tol 1e-10 1e-08`.

**Dataset.**

- dataset: nv_diamond (NV- in diamond), 3E -> 3A2
- data: doc/tutorials/002_nv_diamond_1d_ccd
- roles: 'ES' = ES (upper state, minimum at Q = delta_q), 'GS' = GS (lower state, minimum at Q = 0); emission = ES -> GS = ES -> GS
- energies: ES = ΔSCF, constrained occupations (PBE), 3E; GS = DFT (PBE), 3A2; zero = the GS (3A2) minimum image

Wavefunctions: FD wavefunctions on the polynomial deg. 6. VOF methods: pypl, fd.

Fit **poly6** (polynomial deg. 6):
- GS well: degree-6 polynomial in (Q - 0.000000), coefficients (highest first) [0.1337, -0.2934, 0.2023, -0.0296, 0.4127, 0.001382, 5.824e-06]; fitted to samples 1-15 of 15 (Q = -0.139..0.835); rms residual 2.37e-05 eV
- ES well: degree-6 polynomial in (Q - 0.695732), coefficients (highest first) [0.2858, 0.5198, 0.2939, 0.03818, 0.4955, 0.0007426, 1.706]; fitted to samples 1-11 of 11 (Q = 0.139..0.835); rms residual 6.13e-06 eV

Wavefunction check: max |<chi_i|chi_j> - F_ij| = 8.80e-02 (tolerance 1e-08); mismatch expected, the fit is anharmonic.

## Files

| File | Content |
| --- | --- |
| `log.txt` | Full console output, including every method's settings. |
| `settings.json` | Every parameter of the run. |
| `ccd_wavefunctions.png` | The CCD with 3 wavefunctions per well (FD wavefunctions on the polynomial deg. 6) and the VOF panel. |
