Input and output files of Quantum ESPRESSO and WEST (TDDFT-DDH) calculations for the nitrogen-vacancy center (NV-) in MgO
along the 1D configuration coordinate between the 3E and 1A1 equilibrium structures (Mg107 O108 N1, 216 atoms)

- qe: ground-state (3A2) DFT calculation at every image, `qe/Image_<q>/data-file-schema.xml`
- sc: spin-conserving TDDFT at every image, `sc/Image_<q>/summary.json` (root 1 = 3E; excitation energies in Ry)
- sf: spin-flip TDDFT at every image, `sf/Image_<q>/summary.json` (roots 1-4 = ms0, 1B2, 1B1, 1A1)
- `<sc|sf>/Image_<q>/inputs`: input files, logs and small outputs of each calculation
- reference: the plotting script for the HPC directory layout, and notes
- ccd_nv_mgo.py: plots the CCD from this directory (ccd_nv_mgo.png)
- NOTES.md, MANIFEST.tsv: how and from where the data were collected

Image_<q> is the structure Q = q dQ on the linear path from the 3E (q = 0) to the 1A1 (q = 1) structure.
sf/Image_+2.00 was run without spin flip and is not used.
