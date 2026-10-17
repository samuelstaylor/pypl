# 1D CCD with wavefunctions: wf-fd_vof-pypl-fd_fit-parabola

Generated 2026-10-16T21:07:57 by `python C:\Users\sammy\Research\Galli\programs\pypl\validation\plot_ccd_wavefunctions.py --wf-method fd --vof-method pypl fd --dataset nv_diamond --temp 5.0 --sigma 8.0 25.0 --gamma 0.3 --order-es 50 --order-gs 60 --ene-range -200.0 1200.0 --resol 1401 --line-positions harmonic --npts 100001 --hermite-margin 6.0 --fd-stencil 5 --fd-margin 2.0 --fd-fit parabola --fd-npts 20001 --cubic-warn 0.1 --dpi 300 --fig-format png --font-size 9.0 --n-levels 3 --check-tol 1e-10 1e-08`.

**Dataset.**

- dataset: nv_diamond (NV- in diamond), 3E -> 3A2
- data: doc/tutorials/002_nv_diamond_1d_ccd
- roles: 'ES' = ES (upper state, minimum at Q = delta_q), 'GS' = GS (lower state, minimum at Q = 0); emission = ES -> GS = ES -> GS
- energies: ES = ΔSCF, constrained occupations (PBE), 3E; GS = DFT (PBE), 3A2; zero = the GS (3A2) minimum image

Wavefunctions: FD wavefunctions on the parabola (near min.). VOF methods: pypl, fd.

Fit **parabola** (parabola (near min.)):
- GS well: parabola a (Q - 0.000000)^2 + c, a = 0.417627 eV/(amu A^2), c = -0.000013 eV; fitted to samples 1-5 of 15 (Q = -0.139..0.139); rms residual 7.49e-05 eV
- ES well: parabola a (Q - 0.695732)^2 + c, a = 0.494721 eV/(amu A^2), c = 1.705877 eV; fitted to samples 5-11 of 11 (Q = 0.417..0.835); rms residual 1.49e-04 eV

Wavefunction check: max |<chi_i|chi_j> - F_ij| = 1.63e-12 (tolerance 1e-08).

## Files

| File | Content |
| --- | --- |
| `log.txt` | Full console output, including every method's settings. |
| `settings.json` | Every parameter of the run. |
| `ccd_wavefunctions.png` | The CCD with 3 wavefunctions per well (FD wavefunctions on the parabola (near min.)) and the VOF panel. |
