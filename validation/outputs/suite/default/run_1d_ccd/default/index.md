# Tutorial 002 as a script

Generated 2026-10-16T21:05:57 by `python C:\Users\sammy\Research\Galli\programs\pypl\validation\run_1d_ccd.py --temp 5.0 --ezpl 1945.0 --order-es 50 --order-gs 60 --sigma 8.0 25.0 --gamma 0.3 --ene-range -200.0 1200.0 --resol 1401 --tdm 1.0 --dataset nv_diamond --dpi 300 --fig-format png --font-size 9.0`.

**Dataset.**

- dataset: nv_diamond (NV- in diamond), 3E -> 3A2
- data: doc/tutorials/002_nv_diamond_1d_ccd
- roles: 'ES' = ES (upper state, minimum at Q = delta_q), 'GS' = GS (lower state, minimum at Q = 0); emission = ES -> GS = ES -> GS
- energies: ES = ΔSCF, constrained occupations (PBE), 3E; GS = DFT (PBE), 3A2; zero = the GS (3A2) minimum image

PL peak at 1945.00 meV, absorption peak at 1945.00 meV.

## Files

| File | Content |
| --- | --- |
| `log.txt` | Full console output, including every method's settings. |
| `settings.json` | Every parameter of the run. |
| **`data/`** | Spectra and CCD numbers, as CSV. |
| `data/abs_spectrum.csv` | Absorption. |
| `data/ccd_summary.csv` | delta_q (amu^0.5 A), frequencies (meV), HRFs, run parameters. |
| `data/pl_spectrum.csv` | PL spectrum. |
| **`figures/`** | The CCD diagram and the spectra. |
| `figures/ccd_diagram.png` | Computed images (DFT (PBE), 3A2; ΔSCF, constrained occupations (PBE), 3E), parabola fits, frequencies, HRFs. |
| `figures/spectra.png` | PL and absorption. |
