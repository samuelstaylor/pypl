# 1D CCD with wavefunctions: wf-cubic_vof-pypl-cubic_cubicfit-cubic

Generated 2026-10-16T21:09:09 by `python C:\Users\sammy\Research\Galli\programs\pypl\validation\plot_ccd_wavefunctions.py --wf-method cubic --vof-method pypl cubic --dataset nv_diamond --temp 5.0 --sigma 8.0 25.0 --gamma 0.3 --order-es 50 --order-gs 60 --ene-range -200.0 1200.0 --resol 1401 --line-positions harmonic --npts 100001 --hermite-margin 6.0 --fd-stencil 5 --fd-margin 2.0 --fd-fit parabola --fd-npts 20001 --cubic-warn 0.1 --dpi 300 --fig-format png --font-size 9.0 --cubic-points all --n-levels 3 --check-tol 1e-10 1e-08`.

**Dataset.**

- dataset: nv_diamond (NV- in diamond), 3E -> 3A2
- data: doc/tutorials/002_nv_diamond_1d_ccd
- roles: 'ES' = ES (upper state, minimum at Q = delta_q), 'GS' = GS (lower state, minimum at Q = 0); emission = ES -> GS = ES -> GS
- energies: ES = ΔSCF, constrained occupations (PBE), 3E; GS = DFT (PBE), 3A2; zero = the GS (3A2) minimum image

Wavefunctions: cubic-PT wavefunctions (cubic (a fixed, all images)). VOF methods: pypl, cubic.

Wavefunction check: max |<chi_i|chi_j> - F_ij| = 1.67e-16 (tolerance 1e-10).

## Files

| File | Content |
| --- | --- |
| `log.txt` | Full console output, including every method's settings. |
| `settings.json` | Every parameter of the run. |
| `ccd_wavefunctions.png` | The CCD with 3 wavefunctions per well (cubic-PT wavefunctions (cubic (a fixed, all images))) and the VOF panel. |
