# 1D CCD with wavefunctions: wf-hermite_vof-hermite-pypl-fd_fit-parabola [NV- in MgO]

Generated 2026-10-16T15:20:51 by `python C:\Users\sammy\Research\Galli\programs\pypl\testing\plot_ccd_wavefunctions.py --wf-method hermite --vof-method hermite pypl fd --dataset nv_mgo --temp 5.0 --sigma 20.0 --gamma 0.3 --no-zpl --order-es 160 --order-gs 100 --ene-range -500.0 10000.0 --resol 2101 --npts 100001 --fd-npts 40001 --fd-stencil 5 --fd-margin 2.0 --fd-fit parabola --n-levels 3`.

**Dataset.**

- dataset: nv_mgo (NV- in MgO), 3E -> 1A1 (WEST TDDFT-DDH)
- data: doc/nv_mgo_west
- roles: 'ES' = 3E (upper state, minimum at Q = 0), 'GS' = 1A1 (lower state, minimum at Q = delta_q); emission = ES -> GS = 3E -> 1A1
- Q = q delta_q, q from Image_<q>, as in reference/ccd_nv_mgo.py: 3E minimum at Q = 0, 1A1 minimum at Q = delta_q
- E = QE total energy (qe/) + WEST root (3E: sc root 1, 1A1: sf root 4) x 13.605693 eV/Ry, relative to the 3A2 (QE) minimum, as in ccd_nv_mgo.py
- image windows (q of Image_<q>, inclusive): parabola 3E q = -0.10..+0.45, 1A1 q = +0.60..+1.15; cubic a Q^2 + b Q^3 + c (perturbative correction) q = -0.10..+1.15 for both; every other fit (splines, Morse, polyN) ALL images
- all images: 3E all 33 sc images (q = -1.00..+2.00), 1A1 all 30 sf images except Image_+2.00 (not a spin-flip run, so not a 1A1 point)

Wavefunctions: Hermite. VOF methods: hermite, pypl, fd.

Fit **parabola** (parabola (near min.)):
- 1A1 well: parabola a (Q - 3.097879)^2 + c, a = 0.119071 eV/(amu A^2), c = 2.595544 eV; fitted to samples 18-29 of 30 (Q = 1.859..3.563); rms residual 1.05e-02 eV
- 3E well: parabola a (Q - 0.000000)^2 + c, a = 0.442649 eV/(amu A^2), c = 3.160211 eV; fitted to samples 4-15 of 33 (Q = -0.310..1.394); rms residual 1.89e-02 eV

Wavefunction check: max |<chi_i|chi_j> - F_ij| = 9.26e-23 (tolerance 1e-10).

## Files

| File | Content |
| --- | --- |
| `log.txt` | Full console output, including every method's settings. |
| `settings.json` | Every parameter of the run. |
| `ccd_wavefunctions.png` | The CCD with 3 wavefunctions per well (Hermite) and the VOF panel. |
