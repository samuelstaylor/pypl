# 1D CCD with wavefunctions: wf-fd_vof-hermite-fd_fit-spline [NV- in MgO]

Generated 2026-10-16T21:28:08 by `python C:\Users\sammy\Research\Galli\programs\pypl\validation\plot_ccd_wavefunctions.py --wf-method fd --vof-method hermite fd --dataset nv_mgo --temp 5.0 --sigma 20.0 --gamma 0.3 --no-zpl --order-es 160 --order-gs 100 --ene-range -500.0 10000.0 --resol 2101 --line-positions harmonic --npts 100001 --hermite-margin 6.0 --fd-stencil 5 --fd-margin 2.0 --fd-fit spline --fd-npts 40001 --cubic-warn 0.1 --dpi 300 --fig-format png --font-size 9.0 --n-levels 3 --check-tol 1e-10 1e-08`.

**Dataset.**

- dataset: nv_mgo (NV- in MgO), 3E -> 1A1 (WEST TDDFT-DDH)
- data: doc/tutorials/003_nv_mgo_isc
- roles: 'ES' = 3E (upper state, minimum at Q = 0), 'GS' = 1A1 (lower state, minimum at Q = delta_q); emission = ES -> GS = 3E -> 1A1
- energies: 3E = spin-conserving TDDFT@DDH (WEST, root 1); 1A1 = spin-flip TDDFT@DDH (WEST, root 4); zero = the 3A2 minimum, DFT with the DDH functional (PBE0, exact exchange 0.34)
- Q = q delta_q, q from Image_<q>, as in reference/ccd_nv_mgo.py: 3E minimum at Q = 0, 1A1 minimum at Q = delta_q
- E = QE total energy (qe/) + WEST root (3E: sc root 1, 1A1: sf root 4) x 13.605693 eV/Ry, relative to the 3A2 (QE) minimum, as in ccd_nv_mgo.py
- image windows (q of Image_<q>, inclusive): parabola 3E q = -0.10..+0.45, 1A1 q = +0.60..+1.15; cubic a Q^2 + b Q^3 + c (perturbative correction) q = -0.10..+1.15 for both; every other fit (splines, Morse, polyN) ALL images
- all images: 3E all 33 sc images (q = -1.00..+2.00), 1A1 all 30 sf images except Image_+2.00 (not a spin-flip run, so not a 1A1 point)

Wavefunctions: FD wavefunctions on the cubic spline (natural). VOF methods: hermite, fd.

Fit **spline** (cubic spline (natural)):
- 1A1 well: cubic spline (degree 3, natural ends, 30 knots), harmonic tails a = 0.119071 eV/(amu A^2), C0 junction (slope jumps 1.86e+00 / -5.45e-01 eV/(amu^0.5 A)); fitted to samples 1-30 of 30 (Q = -3.098..4.647); rms residual 2.69e-16 eV
- 3E well: cubic spline (degree 3, natural ends, 33 knots), harmonic tails a = 0.442649 eV/(amu A^2), C0 junction (slope jumps 1.04e+00 / 4.16e+00 eV/(amu^0.5 A)); fitted to samples 1-33 of 33 (Q = -3.098..6.196); rms residual 5.19e-16 eV

Wavefunction check: max |<chi_i|chi_j> - F_ij| = 2.11e-07 (tolerance 1e-08); mismatch expected, the fit is anharmonic.

## Files

| File | Content |
| --- | --- |
| `log.txt` | Full console output, including every method's settings. |
| `settings.json` | Every parameter of the run. |
| `ccd_wavefunctions.png` | The CCD with 3 wavefunctions per well (FD wavefunctions on the cubic spline (natural)) and the VOF panel. |
