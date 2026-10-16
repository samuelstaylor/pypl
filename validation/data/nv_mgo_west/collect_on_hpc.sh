#!/usr/bin/env bash
# collect_on_hpc.sh -- gather the NV-MgO WEST 1D-CCD data needed by the PyPL test harness.
#
# Run ON THE HPC, from anywhere:
#
#     bash collect_on_hpc.sh /path/to/nv_mgo_ccd            # the directory holding sc/ and sf/
#     bash collect_on_hpc.sh /path/to/nv_mgo_ccd ~/xfer      # choose where the bundle is written
#
# It never modifies or deletes anything under the source directory. It writes
#
#     <out>/nv_mgo_west/            the bundle (layout below)
#     <out>/nv_mgo_west.tar.gz      the same, compressed, to copy to the workstation
#
# Bundle layout (mirrors what ccd_nv_mgo.py reads):
#
#     sc/Image_+0.05/summary.json      west.wbse.save/summary.json of the spin-conserving run (3E = state 1)
#     sf/Image_+0.05/summary.json      west.wbse.save/summary.json of the spin-flip run (ms0, 1B2, 1B1, 1A1)
#     qe/Image_+0.05/data-file-schema.xml   ground-state QE XML (total energy + geometry), found exactly as
#                                      resolve_xml() does: via the west.wbse.save symlink's parent
#                                      pwscf.save/, else the image's own pwscf.save/. sc/ is searched
#                                      first, sf/ as the fallback, as in load_ground_state().
#     <state>/<image>/inputs/          small text files of the image (inputs, outputs, logs; < MAX_MB)
#     reference/                       small scripts / notes from the top of the source directory
#     MANIFEST.tsv                     one row per state and image: what was found, from where, checksums
#     NOTES.md                         host, date, source path, counts, code versions found in the outputs
#
# Wavefunctions, charge densities, pseudopotentials and anything larger than MAX_MB are NOT copied.

set -euo pipefail

SRC=${1:?usage: bash collect_on_hpc.sh <dir with sc/ and sf/> [out dir]}
OUT_PARENT=${2:-$PWD}
MAX_MB=${MAX_MB:-50}
SRC=$(cd "$SRC" && pwd -P)
OUT="$OUT_PARENT/nv_mgo_west"

for d in sc sf; do
    [ -d "$SRC/$d" ] || { echo "ERROR: $SRC/$d not found (expected sc/ and sf/ with Image_* directories)"; exit 1; }
done
case "$OUT" in "$SRC"/*) echo "ERROR: write the bundle outside the source tree"; exit 1;; esac
rm -rf "$OUT" "$OUT.tar.gz"
mkdir -p "$OUT/qe" "$OUT/reference"

sha() { sha256sum < "$1" | cut -d' ' -f1; }
size() { stat -L -c %s "$1"; }
small() { [ -f "$1" ] && [ "$(size "$1")" -le $((MAX_MB * 1024 * 1024)) ]; }

# resolve_xml(): print the path of the image's QE XML and how it was found, or nothing.
resolve_xml() {
    local img=$1 link="$1/west.wbse.save" raw cand
    if [ -L "$link" ]; then
        # Absolute or relative-to-the-link target ...
        cand="$(dirname "$(readlink -f "$link")")/pwscf.save/data-file-schema.xml"
        if [ -f "$cand" ]; then echo "$cand	symlink-parent"; return; fi
        # ... or, as os.readlink() in ccd_nv_mgo.py, a relative target taken from the source directory.
        raw=$(readlink "$link")
        cand="$SRC/$(dirname "$raw")/pwscf.save/data-file-schema.xml"
        if [ -f "$cand" ]; then echo "$cand	symlink-parent(relative to source dir)"; return; fi
    fi
    cand="$img/pwscf.save/data-file-schema.xml"
    if [ -f "$cand" ]; then echo "$cand	local-pwscf.save"; fi
}

# Number of eigenvalues and the first four (Ry) of a summary.json, if python3 is available.
PY_BIN=""
for p in python3 python; do
    if command -v "$p" >/dev/null && "$p" -c "import json" >/dev/null 2>&1; then PY_BIN=$p; break; fi
done
eig_info() {
    [ -n "$PY_BIN" ] || { echo "?	?"; return; }
    "$PY_BIN" - "$1" <<'PY' 2>/dev/null || echo "unreadable	?"
import json, sys
e = json.load(open(sys.argv[1]))["plep"]["eigenval"]
print("%d\t%s" % (len(e), ",".join("%.8f" % v for v in e[:4])))
PY
}

printf "state\timage\tq\thas_summary\tn_eigenval\tfirst_eigenval_Ry\tsummary_realpath\txml_found_via\txml_realpath\txml_bytes\txml_sha256\tnotes\n" > "$OUT/MANIFEST.tsv"
declare -A SUMMARY_OWNER=() XML_SHA=()
n_summary=0; n_xml=0; n_images=0

for state in sc sf; do
    for img in "$SRC/$state"/Image_*/; do
        [ -d "$img" ] || continue
        img=${img%/}; name=$(basename "$img"); n_images=$((n_images + 1))
        q=$(echo "$name" | sed -n 's/^Image_\([+-]\{0,1\}[0-9]*\.[0-9]*\)$/\1/p')
        dest="$OUT/$state/$name"; mkdir -p "$dest/inputs"
        notes=""

        summary="$img/west.wbse.save/summary.json"
        has_summary=no; neig="-"; eigs="-"; sreal="-"
        if [ -f "$summary" ]; then
            has_summary=yes; n_summary=$((n_summary + 1))
            cp -L "$summary" "$dest/summary.json"
            sreal=$(readlink -f "$summary")
            IFS=$'\t' read -r neig eigs <<< "$(eig_info "$summary")"
            if [ -n "${SUMMARY_OWNER[$sreal]:-}" ]; then
                notes+="SUMMARY_SHARED_WITH:${SUMMARY_OWNER[$sreal]};"
            fi
            SUMMARY_OWNER[$sreal]="$state/$name"
            for f in "$img"/west.wbse.save/*.json; do       # other small WBSE json files
                [ "$(basename "$f")" = summary.json ] && continue
                small "$f" && cp -L "$f" "$dest/inputs/" || true
            done
        fi

        xml_line=$(resolve_xml "$img" || true)
        xml=${xml_line%%	*}; via=${xml_line#*	}
        xreal="-"; xbytes="-"; xsha="-"
        if [ -n "$xml_line" ]; then
            xreal=$(readlink -f "$xml"); xbytes=$(size "$xml"); xsha=$(sha "$xml")
            if [ -z "${XML_SHA[$name]:-}" ]; then
                mkdir -p "$OUT/qe/$name"; cp -L "$xml" "$OUT/qe/$name/data-file-schema.xml"
                XML_SHA[$name]=$xsha; n_xml=$((n_xml + 1))
            elif [ "${XML_SHA[$name]}" != "$xsha" ]; then
                cp -L "$xml" "$dest/data-file-schema.xml"
                notes+="XML_DIFFERS_FROM_qe/$name;"
            fi
        else
            via="NOT_FOUND"
        fi

        # Small text files directly in the image directory (inputs, outputs, job logs).
        for f in "$img"/*; do
            [ -f "$f" ] || continue
            case "$(basename "$f")" in
                *.in|*.inp|*.yaml|*.yml|*.json|*.out|*.log|*.sh|*.slurm|*.pbs|slurm-*|*.err|*.txt|*.xml) ;;
                *) continue;;
            esac
            if small "$f"; then cp -L "$f" "$dest/inputs/"; else notes+="SKIPPED_LARGE:$(basename "$f");"; fi
        done
        rmdir "$dest/inputs" 2>/dev/null || true

        printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "$state" "$name" "${q:-?}" "$has_summary" \
            "$neig" "$eigs" "$sreal" "$via" "$xreal" "$xbytes" "$xsha" "${notes:--}" >> "$OUT/MANIFEST.tsv"
    done
done

# Scripts and notes at the top of the source directory (e.g. ccd_nv_mgo.py).
for f in "$SRC"/*.py "$SRC"/*.md "$SRC"/*.txt "$SRC"/*.sh "$SRC"/*.in; do
    [ -f "$f" ] && small "$f" && cp -L "$f" "$OUT/reference/" || true
done

{
    echo "# NV-MgO WEST 1D-CCD bundle"
    echo
    echo "- collected: $(date -Iseconds) on $(hostname) by $(whoami)"
    echo "- source: $SRC"
    echo "- image directories: $n_images (sc + sf), with summary.json: $n_summary, ground-state XMLs: $n_xml"
    echo "- sc images: $(ls -d "$SRC"/sc/Image_* 2>/dev/null | wc -l), sf images: $(ls -d "$SRC"/sf/Image_* 2>/dev/null | wc -l)"
    echo "- files larger than ${MAX_MB} MB were not copied (see SKIPPED_LARGE in MANIFEST.tsv)"
    echo
    echo "## Code versions found in the copied outputs"
    echo
    grep -rhoE "Program (PWSCF|WBSE|WSTAT|WEST)[^,]*v[0-9][0-9.a-zA-Z-]*|WEST version[^\"]*" "$OUT" 2>/dev/null | sort -u | sed 's/^/- /' || true
    echo
    echo "## Flags"
    echo
    grep -E "SUMMARY_SHARED|XML_DIFFERS|NOT_FOUND|SKIPPED_LARGE" "$OUT/MANIFEST.tsv" | cut -f1,2,8,12 | sed 's/^/- /' || echo "- none"
} > "$OUT/NOTES.md"

(cd "$OUT_PARENT" && tar -czf nv_mgo_west.tar.gz nv_mgo_west)
echo
echo "Bundle:   $OUT"
echo "Archive:  $OUT.tar.gz ($(du -h "$OUT.tar.gz" | cut -f1)), sha256 $(sha "$OUT.tar.gz")"
echo "Images: $n_images, with summary.json: $n_summary, ground-state XMLs: $n_xml"
echo "Flags:"
grep -cE "SUMMARY_SHARED|XML_DIFFERS|NOT_FOUND|SKIPPED_LARGE" "$OUT/MANIFEST.tsv" | sed 's/^/  rows flagged: /'
