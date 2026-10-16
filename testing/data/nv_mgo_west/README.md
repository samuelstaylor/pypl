# NV⁻ in MgO: 1D CCD from WEST TDDFT (DDH)

Test data for the PyPL harness: the negatively charged nitrogen-vacancy center in MgO, on a 1D
configuration coordinate with many images. Each image has a ground-state QE calculation and two WEST
(WBSE) TDDFT calculations with a DDH functional:

| Run | State index in `summary.json` (`plep.eigenval`) | State |
| --- | --- | --- |
| `sc/` spin-conserving TDDFT | 1 | ³E (minimum at Q = 0) |
| `sf/` spin-flip TDDFT | 1 | m_s = 0 triplet component (nearly flat) |
| | 2 | ¹B₂ |
| | 3 | ¹B₁ |
| | 4 | ¹A₁ (minimum at Q = ΔQ) |

The energy of each state on image Q is the ground-state DFT total energy plus that state's WBSE vertical
excitation energy: E_state(Q) = E_DFT(Q) + VEE_state(Q), with eigenvalues in Ry (× 13.605693 for eV). Images
are named `Image_<q>` (e.g. `Image_+0.05`, `Image_-0.10`), and Q = q·ΔQ. ΔQ is the mass-weighted distance
between `Image_+0.00` and `Image_+1.00`. `reference/ccd_nv_mgo.py` is the user's HPC script that builds this
CCD, kept as the specification.

## Layout (after the transfer)

```
nv_mgo_west/
├── README.md                    this file
├── collect_on_hpc.sh            the HPC-side collector (read-only on the source data)
├── reference/ccd_nv_mgo.py      how the CCD energies are assembled (the user's script)
├── MANIFEST.tsv                 per state and image: what was found, from where, checksums, flags
├── NOTES.md                     host, date, source path, counts, code versions, flags
├── qe/Image_<q>/data-file-schema.xml      ground-state QE XML (total energy + geometry)
├── sc/Image_<q>/summary.json              spin-conserving WBSE results (3E)
├── sc/Image_<q>/inputs/                   that image's small inputs, outputs and logs
├── sf/Image_<q>/summary.json              spin-flip WBSE results (ms0, 1B2, 1B1, 1A1)
└── sf/Image_<q>/inputs/
```

Only the files the CCD needs are transferred. Wavefunctions, charge densities, pseudopotentials and
anything over 50 MB stay on the HPC.

## Getting the data from the HPC

### Option A: let Claude Code on the HPC do it

Log in to the HPC and `cd` anywhere. Replace `/path/to/nv_mgo_ccd` with the directory that holds `sc/`, `sf/`
and `ccd_nv_mgo.py`, then paste:

```bash
claude "I need a small, read-only bundle of my NV- in MgO WEST 1D-CCD data so I can copy it to my workstation. The CCD directory is /path/to/nv_mgo_ccd. It contains sc/ (spin-conserving TDDFT, summary.json state 1 = 3E) and sf/ (spin-flip TDDFT, states 1-4 = ms0, 1B2, 1B1, 1A1), each with Image_<q> subdirectories such as Image_+0.05 and Image_-0.10, plus my plotting script ccd_nv_mgo.py, which reads the data. Do NOT modify, move or delete anything under /path/to/nv_mgo_ccd. First read ccd_nv_mgo.py and inspect a few image directories (including where the west.wbse.save symlinks point), and tell me what you find before copying. Then build ~/nv_mgo_west/ with exactly this layout: sc/<Image>/summary.json and sf/<Image>/summary.json, copied with symlinks dereferenced from <Image>/west.wbse.save/summary.json; qe/<Image>/data-file-schema.xml, the ground-state QE XML of that image, located exactly as resolve_xml() in ccd_nv_mgo.py does (the parent of the west.wbse.save symlink target, then pwscf.save/data-file-schema.xml; fallback <Image>/pwscf.save/data-file-schema.xml), searching sc/ first and sf/ as the fallback, as load_ground_state() does; <state>/<Image>/inputs/ with each image's small text inputs, outputs and logs (*.in, *.yaml, *.json, *.out, *.log, job scripts) under 50 MB; reference/ with ccd_nv_mgo.py and other small top-level scripts or notes. Never copy wavefunctions, charge densities, pseudopotentials or any file over 50 MB. Write MANIFEST.tsv with one row per state and image: q, whether summary.json exists, the number of eigenvalues and the first four in Ry, the real paths of summary.json and of the XML, the XML size and sha256, and notes. Write NOTES.md with host, date, source path, counts, the QE and WEST versions, and the DDH / TDDFT settings found in the inputs and outputs. Flag in both files: images missing summary.json or the XML; summary.json files in sc/ and sf/ that resolve to the same real file; and sc/ and sf/ XMLs for the same image that differ. Then create ~/nv_mgo_west.tar.gz, and finish by printing its size, its sha256, the counts and every flag."
```

Claude Code inspects the layout first, so it can adapt if the symlinks or file names differ from what
`ccd_nv_mgo.py` expects.

### Option B: the plain script

Copy `collect_on_hpc.sh` to the HPC (or paste its contents into a file there) and run:

```bash
bash collect_on_hpc.sh /path/to/nv_mgo_ccd ~      # writes ~/nv_mgo_west/ and ~/nv_mgo_west.tar.gz
```

It does the same as option A, deterministically, and prints the archive's size, its sha256 and any flags.

### Bringing it here

From Git Bash on the workstation (replace the user and host):

```bash
cd /c/Users/sammy/Research/Galli/programs/pypl/testing/data
scp USER@HPC_HOST:~/nv_mgo_west.tar.gz .
sha256sum nv_mgo_west.tar.gz                      # compare with the value printed on the HPC
tar -xzf nv_mgo_west.tar.gz                       # fills nv_mgo_west/ (this README is kept)
```

Run these from that directory, with relative paths: Git Bash's `tar` misreads `C:` in an absolute path as a
remote host.

## Status

Waiting for the data. The harness's loader (`testing/harness/ccd_data.py`) reads the tutorial-002 QE layout,
so this dataset will need its own loader. It must handle `Image_<q>` coordinates, E_DFT + VEE per state, the
Mg / O / N masses, and the energy zero at the ³A₂ minimum. It also needs a choice of which pair of states forms
the transition (³E at Q = 0 and ¹A₁ at Q = ΔQ in the reference script).
