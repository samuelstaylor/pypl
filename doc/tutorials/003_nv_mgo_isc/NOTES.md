# NV⁻ in MgO — WEST 1D-CCD bundle

- Host: `login36.chn` (`perlmutter`)
- Date: 2026-10-15T17:24:52-07:00
- Source: `/global/cfs/cdirs/m4382/sstaylor/projects/nv_mgo/tddft-ddh/phonons/ccd-extend` (read-only; nothing modified)
- Built by: `sstaylor`

## Layout

- `sc/<Image>/summary.json`, `sf/<Image>/summary.json` — dereferenced copies of `<Image>/west.wbse.save/summary.json`
- `qe/<Image>/data-file-schema.xml` — ground-state QE XML, located with `resolve_xml()`; sc/ searched first, sf/ fallback (as `load_ground_state()`)
- `<state>/<Image>/inputs/` — small inputs/outputs/logs (*.in, *.yaml, *.json, *.out, *.log, *.err, *.md, CRASH, pwscf.xml, *.xsf; <50 MB; symlinks dereferenced). Run scripts (*.sh) intentionally omitted.
- `reference/` — `ccd_nv_mgo.py`, `sf_note.md` (from `sf/note.md`)
- No wavefunctions, charge densities, pseudopotentials, `.save` directories or files >50 MB were copied.

State mapping (from ccd_nv_mgo.py): sc state 1 = ³E; sf states 1–4 = m_s=0, ¹B₂, ¹B₁, ¹A₁. Eigenvalues are in Ry.

## Counts

| | image dirs | with data (copied) | with summary.json | with XML | no data (skipped) | input files copied |
|---|---|---|---|---|---|---|
| sc | 46 | 34 | 33 | 34 | 12 | 332 |
| sf | 50 | 31 | 31 | 31 | 19 | 355 |

- qe/ XMLs: 34 (34 from sc, 0 from sf fallback)
- q range with summary.json: sc -1.00 … +2.00; sf -1.00 … +2.00

## Code versions

From `Program … v.` banners in the copied outputs (count = number of files):

- PWSCF v7.5 — 24 files
- WBSE v6.3.0 (v6.3.0) — 66 files
- WBSE_INIT v6.3.0 (v6.3.0) — 37 files
- WESTPP v6.3.0 (v6.3.0) — 158 files
- QE XML `<creator>` in qe/: `PWSCF 7.5` (×34)

## DFT (DDH) settings — from pw.in

DDH here means PBE0 form with the exact-exchange fraction set to α = 1/ε∞ (exx_fraction).

**sc/**

- `calculation`: `"scf"` (×34)
- `input_dft`: `'pbe0'` (×34)
- `exx_fraction`: `0.34` (×34)
- `ecutwfc`: `60` (×34)
- `ecutfock`: `60` (×34)
- `nspin`: `2` (×34)
- `tot_magnetization`: `2` (×34)
- `nbnd`: `900` (×34)
- `occupations`: `'fixed'` (×34)
- `nat`: `216` (×34)
- `ntyp`: `3` (×34)
- `K_POINTS`: `{gamma}` (×34)

**sf/**

- `calculation`: `"scf"` (×31)
- `input_dft`: `'pbe0'` (×31)
- `exx_fraction`: `0.34` (×31)
- `ecutwfc`: `60` (×31)
- `ecutfock`: `60` (×31)
- `nspin`: `2` (×31)
- `tot_magnetization`: `2` (×31)
- `nbnd`: `900` (×31)
- `occupations`: `'fixed'` (×31)
- `nat`: `216` (×31)
- `ntyp`: `3` (×31)
- `K_POINTS`: `{gamma}` (×31)

- `exx_fraction` stored in qe/ XMLs: `0.3400` (×34)
- sc WEST-reported EXX: use_ace `T` (×34), `F` (×7); alpha `0.340000` (×41); ecutfock `240.000000` (×41); exxdiv_treatment `gb` (×41)
- sf WEST-reported EXX: use_ace `T` (×31), `F` (×30); alpha `0.340000` (×61); ecutfock `240.000000` (×61); exxdiv_treatment `gb` (×61)

## TDDFT settings — echoed by WEST in wbse_init.out / wbse.out

| key | sc | sf |
|---|---|---|
| `wbse_init_calculation` | `S` (×41) | `S` (×61) |
| `solver` | `TDDFT` (×34) | `TDDFT` (×31) |
| `bse_method` | `PDEP` (×41) | `PDEP` (×61) |
| `n_pdep_eigen_to_use` | `1724` (×41) | `1724` (×61) |
| `localization` | `W` (×41) | `W` (×61) |
| `overlap_thr` | `0.001` (×41) | `0.001` (×61) |
| `chi_kernel` | `CHI` (×41) | `CHI` (×61) |
| `wbse_calculation` | `D` (×32), `d` (×2) | `D` (×19), `d` (×12) |
| `n_liouville_eigen` | `6` (×34) | `10` (×31) |
| `n_liouville_times` | `8` (×34) | `4` (×31) |
| `n_liouville_maxiter` | `1000` (×34) | `1000` (×31) |
| `trev_liouville` | `1e-06` (×34) | `1e-06` (×31) |
| `trev_liouville_rel` | `0.0001` (×34) | `0.0001` (×31) |
| `l_preconditioning` | `True` (×34) | `True` (×31) |
| `l_pre_shift` | `False` (×34) | `False` (×31) |
| `spin_excitation` | `S` (×34) | `S` (×31) |
| `l_spin_flip` | `False` (×34) | `True` (×31) |
| `l_spin_flip_kernel` | `False` (×34) | `True` (×31) |
| `l_spin_flip_alda0` | `False` (×34) | `False` (×31) |
| `spin_flip_cut` | `0.001` (×34) | `0.001` (×31) |
| `n_exx_lowrank` | `1294` (×34) | `1294` (×31) |
| `scissor_ope` | `0.0` (×34) | `0.0` (×31) |
| `wbse_epsinfty` | `1.0` (×34) | `1.0` (×31) |
| `l_forces` | `False` (×34) | `False` (×31) |

Input-file variants: restarted runs use `wbse_calculation: d` (and some `n_steps_write_restart: 10`) in sc/Image_-1.00, sf/Image_-0.20, sf/Image_-0.60, sf/Image_-1.00.

## FLAGS

- **SAME_SUMMARY_REALPATH**: none
- **MISSING_XML**: none
- **SF_NOT_SPIN_FLIP**: none
- **IDENTICAL_SUMMARY_CONTENT**: none
- **MISSING_SUMMARY**:
  - sc Image_+1.20: QE XML present but no west.wbse.save/summary.json
- **XML_DIFFER**:
  - sc+sf Image_-1.00: sc 1496eb2b81f3 (/global/cfs/cdirs/m4382/sstaylor/projects/nv_mgo/tddft-ddh/phonons/ccd-extend/sc/Image_-1.00/pwscf.save/data-file-schema.xml) vs sf e784fd6c876f (/global/cfs/cdirs/m4382/sstaylor/projects/nv_mgo/tddft-ddh/phonons/ccd-extend/sf/Image_-1.00/pwscf.save/data-file-schema.xml)
  - sc+sf Image_-0.60: sc e84cf8b5dff4 (/global/cfs/cdirs/m4382/sstaylor/projects/nv_mgo/tddft-ddh/phonons/ccd-extend/sc/Image_-0.60/pwscf.save/data-file-schema.xml) vs sf d69267a3ed7b (/global/cfs/cdirs/m4382/sstaylor/projects/nv_mgo/tddft-ddh/phonons/ccd-extend/sf/Image_-0.60/pwscf.save/data-file-schema.xml)
  - sc+sf Image_-0.20: sc 4346da065269 (/global/cfs/cdirs/m4382/sstaylor/projects/nv_mgo/tddft-ddh/phonons/ccd-extend/sc/Image_-0.20/pwscf.save/data-file-schema.xml) vs sf 7f63af2ae466 (/global/cfs/cdirs/m4382/sstaylor/projects/nv_mgo/tddft-ddh/phonons/ccd-extend/sf/Image_-0.20/pwscf.save/data-file-schema.xml)
  - sc+sf Image_+2.00: sc dbff669803cf (/global/cfs/cdirs/m4382/sstaylor/projects/nv_mgo/tddft-ddh/phonons/ccd-extend/sc/Image_+2.00/pwscf.save/data-file-schema.xml) vs sf 9494cf308e25 (/global/cfs/cdirs/m4382/sstaylor/projects/nv_mgo/tddft-ddh/phonons/ccd-extend/sf/Image_+2.00/pwscf.save/data-file-schema.xml)
- **XML_PATH_NOTE**:
  - sc Image_+1.10: XML resolves to /global/cfs/cdirs/m4382/sstaylor/projects/nv_mgo/tddft-ddh/phonons/ccd/sf/Image-25-og-attempt/pwscf.save/data-file-schema.xml
  - sf Image_+1.10: XML resolves to /global/cfs/cdirs/m4382/sstaylor/projects/nv_mgo/tddft-ddh/phonons/ccd/sf/Image-25-og-attempt/pwscf.save/data-file-schema.xml
- **GS_EXCLUDED_IN_SCRIPT**:
  - qe Image_-1.00: XML copied, but ccd_nv_mgo.py GROUND_STATE excludes this image
- **NO_DATA** (sc, 12 image dirs with only template inputs, skipped): Image_-2.00, Image_-1.50, Image_-1.40, Image_-1.30, Image_-1.20, Image_-1.10, Image_-0.90, Image_-0.80, Image_-0.70, Image_-0.50, Image_-0.40, Image_-0.30
- **NO_DATA** (sf, 19 image dirs with only template inputs, skipped): Image_-2.00, Image_-1.50, Image_-1.40, Image_-1.30, Image_-1.25, Image_-1.20, Image_-1.10, Image_-0.90, Image_-0.80, Image_-0.75, Image_-0.70, Image_-0.50, Image_-0.40, Image_-0.30, Image_-0.25, Image_+1.20, Image_+1.25, Image_+1.30, Image_+1.40

See MANIFEST.tsv for per-image eigenvalues, real paths and sha256 sums.
