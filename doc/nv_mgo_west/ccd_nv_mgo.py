"""
ccd_nv_mgo.py

1D configuration-coordinate diagram (CCD) for NV-MgO.

Spin-conserving (sc) TDDFT:
  state 1 → ^3E  (triplet excited state, Q=0 minimum)

Spin-flip (sf) TDDFT:
  state 1 → m_s=0 triplet component (low-lying, nearly flat)
  state 2 → ^1B_2
  state 3 → ^1B_1
  state 4 → ^1A_1 (singlet excited state, Q=ΔQ minimum)

Directory naming: Image_+0.00, Image_+0.05, Image_-0.10, …
The coordinate value is parsed directly from the directory name.

Only image directories that contain a completed WBSE summary.json are plotted.

Two directory layouts are supported:

  * Collected bundle (this directory):
      sc/Image_<q>/summary.json, sf/Image_<q>/summary.json
      qe/Image_<q>/data-file-schema.xml   (ground-state QE XML, shared by sc/sf)
  * Original HPC tree:
      <sc|sf>/Image_<q>/west.wbse.save/summary.json, with the QE XML located via
      the west.wbse.save symlink (fallback: local pwscf.save/data-file-schema.xml).

Paths are resolved relative to this script, so it can be run from anywhere.
"""

import os
import re
import sys
import json
import numpy as np
import matplotlib.pyplot as plt

# The messages use Δ, …, ^ etc.; Windows consoles default to cp1252, which cannot encode them.
if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

try:
    from pypl.utils import parse_atoms_qexml, parse_total_energy_qexml
except ImportError:
    # pypl not installed: fall back to the repository this script lives in (doc/nv_mgo_west → ../..)
    sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..")))
    from pypl.utils import parse_atoms_qexml, parse_total_energy_qexml


# ===========================================================================
# CONSTANTS  (do not edit)
# ===========================================================================

RY_TO_EV = 13.605693
EPS = 1e-12

MASS_LIST = {
    "H":  1.0079,
    "C":  12.0107,
    "N":  14.0067,
    "O":  15.999,
    "Mg": 24.305,
    "V":  50.9415,
}


# ===========================================================================
# USER INPUTS — COLORS
# ===========================================================================

COLOR_RED    = "#DB4437"
COLOR_BLUE   = "#4285F4"
COLOR_GREEN  = "#0F9D58"
COLOR_YELLOW = "#E6A817"   # darker gold — readable on white
COLOR_PURPLE = "#9C27B0"
COLOR_GRAY   = "#5F6368"
COLOR_TEAL   = "#00897B"
COLOR_ORANGE = "#F57C00"


# ===========================================================================
# USER INPUTS — PATHS
# ===========================================================================

BASE_DIR = os.path.dirname(os.path.abspath(__file__))   # directory holding sc/, sf/, qe/

SC_DIR = os.path.join(BASE_DIR, "sc")   # spin-conserving TDDFT image directories
SF_DIR = os.path.join(BASE_DIR, "sf")   # spin-flip     TDDFT image directories
QE_DIR = os.path.join(BASE_DIR, "qe")   # collected ground-state QE XMLs (qe/Image_*/data-file-schema.xml)

# Image directories used to compute ΔQ (the displacement between the two
# relaxed structures).  Must live inside SC_DIR.
INITIAL_IMAGE = "Image_+0.00"   # Q = 0    (^3E minimum)
FINAL_IMAGE   = "Image_+1.00"   # Q = ΔQ   (^1A_1 minimum; the ^3A_2 minimum is also near here)

# Global energy zero.
#   None  → auto-computed as the minimum DFT total energy across all images
#           (i.e. the ground-state minimum is set to 0 eV).
#   float → use this value in eV directly (e.g. a previously saved minimum).
GLOBAL_ENERGY_ZERO_EV = None


# ===========================================================================
# USER INPUTS — STATES TO PLOT
# ===========================================================================
#
# Each entry is a dict with:
#   key          : unique string identifier (used internally)
#   label        : LaTeX legend label
#   dir          : directory containing the Image_* folders for this state
#   color        : hex color string
#   marker       : matplotlib marker character
#   include_vee  : bool — add the WBSE vertical excitation energy on top of
#                  the DFT total energy for this state
#   state_index  : 1-based index into summary.json "eigenval" list
#                  (only used when include_vee=True)
#   exclude_dirs : list of "Image_*" names to skip
#   connect      : bool — draw a line through the raw points (None → global
#                  CONNECT_POINTS default)
# ===========================================================================

# Images to skip for every spin-flip state.  sf/Image_+2.00 was run without
# spin flip (wbse.out: l_spin_flip = False, 6 roots; its summary.json is
# byte-identical to sc/Image_+2.00), so its states 1-4 are NOT ms0/1B2/1B1/1A1.
SF_EXCLUDE = ["Image_+2.00"]

STATES = [
    # -----------------------------------------------------------------------
    # Ground state (^3A_2) — DFT total energy only, no VEE
    # Loaded separately via GROUND_STATE below; listed here for reference.
    # -----------------------------------------------------------------------

    # -----------------------------------------------------------------------
    # Spin-conserving TDDFT  (sc/)
    # -----------------------------------------------------------------------
    {
        "key":         "3E",
        "label":       r"$^3E$",
        "dir":         SC_DIR,
        "color":       COLOR_RED,
        "marker":      "o",
        "include_vee": True,
        "state_index": 1,          # SC state 1
        "exclude_dirs": [],
        "connect":     False,
    },

    # -----------------------------------------------------------------------
    # Spin-flip TDDFT  (sf/)
    # -----------------------------------------------------------------------
    {
        "key":         "ms0",
        "label":       r"$m_s{=}0$ triplet",
        "dir":         SF_DIR,
        "color":       COLOR_GRAY,
        "marker":      "x",
        "include_vee": True,
        "state_index": 1,          # SF state 1 (m_s=0 component, nearly flat)
        "exclude_dirs": SF_EXCLUDE,
        "connect":     False,
    },
    {
        "key":         "1B2",
        "label":       r"$^1B_2$",
        "dir":         SF_DIR,
        "color":       COLOR_GREEN,
        "marker":      "s",
        "include_vee": True,
        "state_index": 2,          # SF state 2
        "exclude_dirs": SF_EXCLUDE,
        "connect":     False,
    },
    {
        "key":         "1B1",
        "label":       r"$^1B_1$",
        "dir":         SF_DIR,
        "color":       COLOR_ORANGE,
        "marker":      "^",
        "include_vee": True,
        "state_index": 3,          # SF state 3
        "exclude_dirs": SF_EXCLUDE,
        "connect":     False,
    },
    {
        "key":         "1A1",
        "label":       r"$^1A_1$",
        "dir":         SF_DIR,
        "color":       COLOR_BLUE,
        "marker":      "D",
        "include_vee": True,
        "state_index": 4,          # SF state 4
        "exclude_dirs": SF_EXCLUDE,
        "connect":     False,
    },
]


# ===========================================================================
# USER INPUTS — GROUND STATE OVERLAY
# ===========================================================================
# DFT total energy only (no VEE).  The ground-state SCF is identical in sc/ and
# sf/, so each image is taken from the first directory in "dirs" that has a
# QE XML for it (sc/ first, sf/ as fallback for images not yet linked in sc/).

GROUND_STATE = {
    "enabled":     True,
    "label":       r"$^3A_2$",
    "dirs":        [SC_DIR, SF_DIR],
    "color":       COLOR_PURPLE,
    "marker":      "o",
    "exclude_dirs": ["Image_-1.00"],
    "connect":     False,
}


# ===========================================================================
# USER INPUTS — ZPL REFERENCE PAIR
# ===========================================================================
# ZPL energy = E(^1A_1 at Q=ΔQ) − E(^3E at Q=0).
# Set to None to skip ZPL / Franck-Condon analysis.

ZPL_INITIAL_KEY = "3E"    # key of the state whose Q=0 point defines the ZPL
ZPL_FINAL_KEY   = "1A1"   # key of the state whose Q=ΔQ point defines the ZPL


# ===========================================================================
# USER INPUTS — PLOT APPEARANCE
# ===========================================================================

FIGURE_SIZE    = (6, 7)
XLIM           = [None, None]       # x-axis limits in amu^0.5·Å  (None = auto from data)
YLIM           = [None, None]       # y-axis limits in eV           (None = auto from data)
AXIS_PAD       = 0.05               # fractional padding around data for auto limits
XLABEL         = r"Q ($\sqrt{\mathrm{amu}}$·Å)"
YLABEL         = r"Energy relative to $^3A_2$ minimum (eV)"
VLINE_COLOR    = "gray"
CONNECT_POINTS = False              # global default for per-state "connect"
OUTPUT_FILE    = os.path.join(BASE_DIR, "ccd_nv_mgo.png")
DPI            = 300


# ===========================================================================
# HELPER FUNCTIONS
# ===========================================================================

def parse_image_coord(img_name):
    """Return the float coordinate encoded in 'Image_+0.05', etc."""
    m = re.match(r"Image_([+-]?\d+\.\d+)$", img_name)
    return float(m.group(1)) if m else None


def resolve_xml(img_dir):
    """
    Locate data-file-schema.xml for an image directory.

    1. Collected bundle: QE_DIR/<Image_*>/data-file-schema.xml (same image name).
    2. Follow the west.wbse.save symlink to the original ccd Image-N directory
       and return its pwscf.save/data-file-schema.xml.
    3. Fall back to a local pwscf.save/data-file-schema.xml inside img_dir.

    Never use <img_dir>/inputs/pwscf.xml: those are nscf outputs with etot = 0.
    """
    xml = os.path.join(QE_DIR, os.path.basename(os.path.normpath(img_dir)), "data-file-schema.xml")
    if os.path.exists(xml):
        return xml
    wbse_link = os.path.join(img_dir, "west.wbse.save")
    if os.path.islink(wbse_link):
        parent = os.path.dirname(os.readlink(wbse_link))
        xml = os.path.join(parent, "pwscf.save", "data-file-schema.xml")
        if os.path.exists(xml):
            return xml
    xml = os.path.join(img_dir, "pwscf.save", "data-file-schema.xml")
    return xml if os.path.exists(xml) else None


def compute_delta_q(xml0, xml1):
    syms1, coords1, cell1 = parse_atoms_qexml(xml0)
    syms2, coords2, cell2 = parse_atoms_qexml(xml1)
    if syms1 != syms2:
        raise ValueError("Atomic species mismatch between initial and final structures.")
    if np.max(np.abs(cell1 - cell2)) > EPS:
        raise ValueError("Unit cell mismatch between initial and final structures.")
    masses = np.array([MASS_LIST[s] for s in syms1])
    return float(np.linalg.norm((coords2 - coords1) * np.sqrt(masses)[:, None]))


def resolve_summary(img_dir):
    """Return the WBSE summary.json of an image (bundle: img_dir/summary.json;
    HPC tree: img_dir/west.wbse.save/summary.json), or None if it is missing."""
    for path in (os.path.join(img_dir, "summary.json"),
                 os.path.join(img_dir, "west.wbse.save", "summary.json")):
        if os.path.exists(path):
            return path
    return None


def _image_dirs_sorted(data_dir):
    return sorted(
        [d for d in os.listdir(data_dir) if d.startswith("Image_")],
        key=lambda d: (parse_image_coord(d) or 0.0),
    )


def load_state(state_cfg, dq):
    """Load (Q, E_total+VEE) for one excited state. Skips dirs without summary.json."""
    data_dir    = state_cfg["dir"]
    incl_vee    = state_cfg.get("include_vee", False)
    si          = state_cfg.get("state_index", 1)
    excl        = set(state_cfg.get("exclude_dirs", []))
    label       = state_cfg.get("label", "")

    Q_list, E_list = [], []

    for img_name in _image_dirs_sorted(data_dir):
        if img_name in excl:
            continue
        q_rel = parse_image_coord(img_name)
        if q_rel is None:
            continue

        img_dir = os.path.join(data_dir, img_name)
        summary_json = resolve_summary(img_dir)
        if summary_json is None:
            continue   # skip silently — only plot dirs with completed WBSE

        xml_path = resolve_xml(img_dir)
        if xml_path is None:
            print(f"  [{label}] QE XML not found for {img_name}, skipping.")
            continue

        e_dft = parse_total_energy_qexml(xml_path)

        vee = 0.0
        if incl_vee:
            eigvals = json.load(open(summary_json))["plep"]["eigenval"]
            if si - 1 >= len(eigvals):
                print(f"  [{label}] state_index={si} out of range in {img_name}, skipping.")
                continue
            vee = eigvals[si - 1] * RY_TO_EV

        Q_list.append(q_rel * dq)
        E_list.append(e_dft + vee)

    Q, E = np.array(Q_list), np.array(E_list)
    print(f"  [{label}] {len(Q)} points")
    return Q, E


def load_ground_state(gs_cfg, dq):
    """Load (Q, E_DFT) for the ground state (no VEE). Includes all images with a QE XML.

    Images are searched in each directory of gs_cfg["dirs"] in order; the first
    directory providing a QE XML for a given image wins.
    """
    excl     = set(gs_cfg.get("exclude_dirs", []))
    label    = gs_cfg.get("label", "GS")

    found = {}   # img_name -> xml_path
    for data_dir in gs_cfg["dirs"]:
        for img_name in _image_dirs_sorted(data_dir):
            if img_name in excl or img_name in found:
                continue
            if parse_image_coord(img_name) is None:
                continue
            xml_path = resolve_xml(os.path.join(data_dir, img_name))
            if xml_path is not None:
                found[img_name] = xml_path

    Q_list, E_list = [], []
    for img_name in sorted(found, key=parse_image_coord):
        Q_list.append(parse_image_coord(img_name) * dq)
        E_list.append(parse_total_energy_qexml(found[img_name]))

    Q, E = np.array(Q_list), np.array(E_list)
    print(f"  [{label}] {len(Q)} points")
    return Q, E


# ===========================================================================
# COMPUTE ΔQ
# ===========================================================================

print("Computing ΔQ …")
_xml0 = resolve_xml(os.path.join(SC_DIR, INITIAL_IMAGE))
_xml1 = resolve_xml(os.path.join(SC_DIR, FINAL_IMAGE))

if _xml0 is None or _xml1 is None:
    raise FileNotFoundError(
        f"Cannot locate QE XML for {INITIAL_IMAGE} or {FINAL_IMAGE} in {SC_DIR}."
    )

delta_q = compute_delta_q(_xml0, _xml1)
print(f"  ΔQ = {delta_q:.6f}  amu^0.5·Å\n")


# ===========================================================================
# LOAD DATA
# ===========================================================================

print("Loading states …")
results = {}

for state in STATES:
    Q, E = load_state(state, delta_q)
    results[state["key"]] = {"cfg": state, "Q": Q, "E": E}

if GROUND_STATE.get("enabled", False):
    Q, E = load_ground_state(GROUND_STATE, delta_q)
    results["__gs__"] = {"cfg": GROUND_STATE, "Q": Q, "E": E}


# ===========================================================================
# GLOBAL ENERGY SHIFT
# ===========================================================================

if GLOBAL_ENERGY_ZERO_EV is None:
    # Use the minimum DFT total energy (ground state minimum → 0 eV).
    gs_E = results["__gs__"]["E"] if "__gs__" in results else None
    if gs_E is not None and len(gs_E):
        GLOBAL_ENERGY_ZERO_EV = float(np.min(gs_E))
    else:
        all_E = np.concatenate([r["E"] for r in results.values() if len(r["E"])])
        GLOBAL_ENERGY_ZERO_EV = float(np.min(all_E))
    print(f"\nGlobal energy zero (^3A_2 minimum): {GLOBAL_ENERGY_ZERO_EV:.6f} eV")
else:
    print(f"\nUser-supplied energy zero: {GLOBAL_ENERGY_ZERO_EV:.6f} eV")

for r in results.values():
    r["E"] = r["E"] - GLOBAL_ENERGY_ZERO_EV


# ===========================================================================
# ZPL (FROM RAW DATA)
# ===========================================================================

if ZPL_INITIAL_KEY and ZPL_FINAL_KEY:
    ri = results.get(ZPL_INITIAL_KEY)
    rf = results.get(ZPL_FINAL_KEY)
    if ri is not None and rf is not None:
        idx0 = np.where(np.isclose(ri["Q"], 0.0,      atol=1e-4))[0]
        idxf = np.where(np.isclose(rf["Q"], delta_q,  atol=1e-4))[0]
        if len(idx0) and len(idxf):
            zpl_raw = float(rf["E"][idxf[0]] - ri["E"][idx0[0]])
            print(f"ZPL (data):  {zpl_raw:.4f} eV  ({zpl_raw*1000:.1f} meV)")


# ===========================================================================
# PLOT
# ===========================================================================

fig, ax = plt.subplots(figsize=FIGURE_SIZE)

# Draw ground state and excited states.  Ground state goes first so excited
# states render on top.
plot_order = ["__gs__"] + [s["key"] for s in STATES]

for key in plot_order:
    r = results.get(key)
    if r is None or len(r["Q"]) == 0:
        continue
    cfg = r["cfg"]
    Q, E = r["Q"], r["E"]
    order = np.argsort(Q)

    ax.scatter(Q, E,
               marker=cfg.get("marker", "o"),
               color=cfg["color"],
               label=cfg["label"],
               zorder=3, s=28, linewidths=0.5)

    do_connect = cfg.get("connect")
    if do_connect is None:
        do_connect = CONNECT_POINTS
    if do_connect:
        ax.plot(Q[order], E[order], color=cfg["color"], linewidth=0.8, zorder=2)

# Axis limits: auto-scale to the plotted data (plus padding) unless overridden.
all_Q = np.concatenate([results[k]["Q"] for k in plot_order if k in results])
all_E = np.concatenate([results[k]["E"] for k in plot_order if k in results])
q_lo, q_hi = min(all_Q.min(), 0.0, delta_q), max(all_Q.max(), 0.0, delta_q)
e_lo, e_hi = all_E.min(), all_E.max()
q_pad, e_pad = AXIS_PAD * (q_hi - q_lo), AXIS_PAD * (e_hi - e_lo)
ax.set_xlim(XLIM[0] if XLIM[0] is not None else q_lo - q_pad,
            XLIM[1] if XLIM[1] is not None else q_hi + q_pad)
ax.set_ylim(YLIM[0] if YLIM[0] is not None else e_lo - e_pad,
            YLIM[1] if YLIM[1] is not None else e_hi + e_pad)

# Vertical dashed lines at Q=0 and Q=ΔQ, labelled just above the top spine
for xv, lbl in [(0, r"$Q_0$"), (delta_q, r"$Q_{\Delta}$")]:
    ax.axvline(xv, linestyle="--", color=VLINE_COLOR, linewidth=0.8, zorder=1)
    ax.text(xv, 1.01, lbl, transform=ax.get_xaxis_transform(),
            ha="center", va="bottom", fontsize=8, color=VLINE_COLOR)

ax.set_xlabel(XLABEL, fontsize=11)
ax.set_ylabel(YLABEL, fontsize=11)
ax.tick_params(direction="in", which="both", top=True, right=True)
plt.tight_layout()

# Legend: place it at the first location whose box contains no data points.
# Falls back to outside the axes (right side) if every inside spot is taken.
pts_disp = ax.transData.transform(np.column_stack([all_Q, all_E]))
renderer = fig.canvas.get_renderer()
leg_kw = dict(fontsize=8, framealpha=0.9, edgecolor="0.8")
leg = None
for loc in ["upper left", "upper right", "lower right", "lower left",
            "center right", "center left", "upper center", "lower center", "center"]:
    leg = ax.legend(loc=loc, **leg_kw)
    bbox = leg.get_window_extent(renderer).expanded(1.05, 1.05)
    inside = ((pts_disp[:, 0] >= bbox.x0) & (pts_disp[:, 0] <= bbox.x1) &
              (pts_disp[:, 1] >= bbox.y0) & (pts_disp[:, 1] <= bbox.y1))
    if not inside.any():
        break
    leg.remove()
    leg = None
if leg is None:
    ax.legend(loc="upper left", bbox_to_anchor=(1.02, 1.0), borderaxespad=0, **leg_kw)
    plt.tight_layout()

plt.savefig(OUTPUT_FILE, dpi=DPI, bbox_inches="tight")
print(f"\nSaved: {OUTPUT_FILE}")
