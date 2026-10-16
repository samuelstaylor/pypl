# Tutorial 002 as a script

Generated 2026-10-15T18:27:52 by `python C:\Users\sammy\Research\Galli\programs\pypl\testing\run_1d_ccd.py --temp 5.0 --ezpl 1945.0 --order-es 50 --order-gs 60 --sigma 8.0 25.0 --gamma 0.3 --ene-range -200.0 1200.0 --resol 1401 --tdm 1.0`.

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
| `figures/ccd_diagram.png` | DFT images, parabola fits, frequencies, HRFs. |
| `figures/spectra.png` | PL and absorption. |
