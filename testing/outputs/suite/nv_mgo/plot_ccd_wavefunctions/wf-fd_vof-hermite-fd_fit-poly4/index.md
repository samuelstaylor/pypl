# 1D CCD with wavefunctions: wf-fd_vof-hermite-fd_fit-poly4 [NV- in MgO]

Generated 2026-10-16T15:24:02 by `python C:\Users\sammy\Research\Galli\programs\pypl\testing\plot_ccd_wavefunctions.py --wf-method fd --vof-method hermite fd --dataset nv_mgo --temp 5.0 --sigma 20.0 --gamma 0.3 --no-zpl --order-es 160 --order-gs 100 --ene-range -500.0 10000.0 --resol 2101 --npts 100001 --fd-npts 40001 --fd-stencil 5 --fd-margin 2.0 --fd-fit poly4 --n-levels 3`.

**Dataset.**

- dataset: nv_mgo (NV- in MgO), 3E -> 1A1 (WEST TDDFT-DDH)
- data: doc/nv_mgo_west
- roles: 'ES' = 3E (upper state, minimum at Q = 0), 'GS' = 1A1 (lower state, minimum at Q = delta_q); emission = ES -> GS = 3E -> 1A1
- Q = q delta_q, q from Image_<q>, as in reference/ccd_nv_mgo.py: 3E minimum at Q = 0, 1A1 minimum at Q = delta_q
- E = QE total energy (qe/) + WEST root (3E: sc root 1, 1A1: sf root 4) x 13.605693 eV/Ry, relative to the 3A2 (QE) minimum, as in ccd_nv_mgo.py
- image windows (q of Image_<q>, inclusive): parabola 3E q = -0.10..+0.45, 1A1 q = +0.60..+1.15; cubic a Q^2 + b Q^3 + c (perturbative correction) q = -0.10..+1.15 for both; every other fit (splines, Morse, polyN) ALL images
- all images: 3E all 33 sc images (q = -1.00..+2.00), 1A1 all 30 sf images except Image_+2.00 (not a spin-flip run, so not a 1A1 point)

Wavefunctions: finite difference on polynomial deg. 4. VOF methods: hermite, fd.

Fit **poly4** (polynomial deg. 4):
- 1A1 well: degree-4 polynomial in (Q - 3.097879), coefficients (highest first) [0.009277, 0.04862, 0.1579, 0.06138, 2.633]; fitted to samples 1-30 of 30 (Q = -3.098..4.647); rms residual 1.37e-01 eV
- 3E well: degree-4 polynomial in (Q - 0.000000), coefficients (highest first) [0.01143, -0.1145, 0.297, 0.3031, 3.268]; fitted to samples 1-33 of 33 (Q = -3.098..6.196); rms residual 1.90e-01 eV

Wavefunction check: max |<chi_i|chi_j> - F_ij| = 1.75e-07 (tolerance 1e-08); mismatch expected, the fit is anharmonic.

## Files

| File | Content |
| --- | --- |
| `log.txt` | Full console output, including every method's settings. |
| `settings.json` | Every parameter of the run. |
| `ccd_wavefunctions.png` | The CCD with 3 wavefunctions per well (finite difference on polynomial deg. 4) and the VOF panel. |
