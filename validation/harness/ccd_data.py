"""The 1D configuration-coordinate diagrams the harness can run on.

Every script works on one *dataset*, chosen with ``--dataset`` (or the environment variable
``PYPL_TESTING_DATASET``; the default is ``nv_diamond``):

    nv_diamond   NV- in diamond, 3E -> 3A2, tutorial 002 (doc/tutorials/002_nv_diamond_1d_ccd)
    nv_mgo       NV- in MgO, 3E -> 1A1, WEST TDDFT-DDH (doc/nv_mgo_west)

``build_ccd()`` turns the active dataset into the three numbers ``config_coord_1d_solver`` needs
(delta_q, gs_phonon, es_phonon) plus the raw samples every other fit is made from. The harmonic wells
are ``config_coord_1d_potential`` parabolas fitted to the images near each minimum.

The two states are always called ``"gs"`` and ``"es"`` internally: ``"gs"`` is the lower state and ``"es"``
the upper one; emission is es -> gs. Where each minimum sits is the dataset's choice (``center``): diamond
has the GS at Q = 0 and the ES at delta_q, as in tutorial 002; MgO has 3E ("es") at Q = 0 and 1A1 ("gs") at
Q = delta_q, as in doc/nv_mgo_west/reference/ccd_nv_mgo.py. The solver methods take both minima, so every
F_ij is the same in either orientation.

``STATE_NAMES`` (plain text, for logs and CSV), ``STATE_LABELS`` (mathtext, for figures), ``STATE_COLORS`` and
``STATE_MARKERS`` describe the states of the active dataset; ``use_dataset`` updates them in place, so every
module that imported them follows.
"""

import json
import os
import re
from dataclasses import dataclass

import numpy as np
from scipy import constants

from harness import REPO_ROOT
from pypl.utils import parse_atoms_qexml, parse_total_energy_qexml
from pypl.config_coord_1d_potential import config_coord_1d_potential

STATES = ("gs", "es")
AMU = constants.physical_constants["atomic mass constant"][0]
RY_TO_EV = 13.605693

# ------------------------------------------------------------------------- NV- in diamond (tutorial 002)

DATA_DIR = os.path.join(REPO_ROOT, "doc", "tutorials", "002_nv_diamond_1d_ccd")

MASS_LIST = {"C": 12.0107, "N": 14.0067}
# Interpolation path of the images: 0.0 = GS geometry, 1.0 = ES geometry.
RELATIVE_COORDINATES = np.linspace(-0.2, 1.2, 15)
# Images the harmonic parabola is fitted to, as in the tutorial notebook.
GS_FIT_SLICE = slice(0, 5)
ES_FIT_SLICE = slice(4, None)
FIT_SLICES = {"gs": GS_FIT_SLICE, "es": ES_FIT_SLICE}

# ------------------------------------------------------------------------------ NV- in MgO (WEST)

MGO_DATA_DIR = os.path.join(REPO_ROOT, "doc", "nv_mgo_west")
MGO_MASS_LIST = {"N": 14.0067, "O": 15.999, "Mg": 24.305}
# Path coordinate q of the image directories (Image_<q>): 0 = 3E minimum, 1 = 1A1 minimum, and
# Q = q * delta_q, exactly as in doc/nv_mgo_west/reference/ccd_nv_mgo.py.
MGO_STATES = {
    # state: (WEST run directory, 0-based root in summary.json eigenval, physical name)
    "gs": ("sf", 3, "1A1"),
    "es": ("sc", 0, "3E"),
}
# Every image with a summary.json and a QE XML is used by every fit, however far from the minimum, except:
MGO_EXCLUDE = {
    # sf/Image_+2.00 was run without spin flip (wbse.out: l_spin_flip = False; its summary.json is
    # byte-identical to sc/Image_+2.00), so its root 4 is not 1A1: not a 1A1 data point at all.
    "gs": ("Image_+2.00",),
    "es": (),
}
# Image windows (path coordinate q of Image_<q>, inclusive) of the two restricted fits. Every other fit
# (splines, Morse, polyN, ...) uses every image.
MGO_PARABOLA_WINDOW = {"es": (-0.10, 0.45), "gs": (0.60, 1.15)}  # harmonic parabola: 3E / 1A1
MGO_CUBIC_WINDOW = {"es": (-0.10, 1.15), "gs": (-0.10, 1.15)}  # cubic a Q^2 + b Q^3 + c (perturbative correction)
# The QE (3A2) images whose lowest energy is the energy zero, as GROUND_STATE in ccd_nv_mgo.py.
MGO_ZERO_EXCLUDE = ("Image_-1.00",)


@dataclass(frozen=True)
class Dataset:
    """One 1D CCD the harness can run on."""

    name: str
    title: str  # plot titles, e.g. "NV$^-$ in diamond"
    transition: str  # the es -> gs transition, plain text
    state_names: dict  # {"gs": ..., "es": ...}, plain text
    data_dir: str
    builder: str  # name of the build function in this module
    notes: tuple = ()  # lines for logs and index.md
    es_at_zero: bool = False  # True: the upper state's minimum is at Q = 0 and the lower state's at delta_q
    state_labels: dict = None  # mathtext labels for figures (default: state_names)
    state_colors: dict = None
    state_markers: dict = None
    cubic_images: str = None  # the images of the cubic fit, if not all of them (for labels)

    def describe(self):
        """Printable lines: what the data is and which physical state plays which role."""

        return [
            "dataset: %s (%s), %s" % (self.name, self.title.replace("$^-$", "-"), self.transition),
            "data: %s" % os.path.relpath(self.data_dir, REPO_ROOT).replace("\\", "/"),
            "roles: 'ES' = %s (upper state, minimum at Q = %s), 'GS' = %s (lower state, minimum at Q = %s); "
            "emission = ES -> GS = %s -> %s"
            % (
                self.state_names["es"],
                "0" if self.es_at_zero else "delta_q",
                self.state_names["gs"],
                "delta_q" if self.es_at_zero else "0",
                self.state_names["es"],
                self.state_names["gs"],
            ),
        ] + list(self.notes)


DATASETS = {
    d.name: d
    for d in (
        Dataset(
            "nv_diamond",
            "NV$^-$ in diamond",
            "3E -> 3A2",
            {"gs": "GS", "es": "ES"},
            DATA_DIR,
            "_build_nv_diamond",
            state_colors={"gs": "#DB4437", "es": "#4285F4"},
            state_markers={"gs": "o", "es": "s"},
        ),
        Dataset(
            "nv_mgo",
            "NV$^-$ in MgO",
            "3E -> 1A1 (WEST TDDFT-DDH)",
            {"gs": "1A1", "es": "3E"},
            MGO_DATA_DIR,
            "_build_nv_mgo",
            (
                "Q = q delta_q, q from Image_<q>, as in reference/ccd_nv_mgo.py: 3E minimum at Q = 0, "
                "1A1 minimum at Q = delta_q",
                "E = QE total energy (qe/) + WEST root (3E: sc root 1, 1A1: sf root 4) x 13.605693 eV/Ry, "
                "relative to the 3A2 (QE) minimum, as in ccd_nv_mgo.py",
                "image windows (q of Image_<q>, inclusive): parabola 3E q = -0.10..+0.45, 1A1 q = +0.60..+1.15; "
                "cubic a Q^2 + b Q^3 + c (perturbative correction) q = -0.10..+1.15 for both; every other fit "
                "(splines, Morse, polyN) ALL images",
                "all images: 3E all 33 sc images (q = -1.00..+2.00), 1A1 all 30 sf images except Image_+2.00 "
                "(not a spin-flip run, so not a 1A1 point)",
            ),
            es_at_zero=True,
            state_labels={"gs": r"$^1A_1$", "es": r"$^3E$"},
            state_colors={"gs": "#4285F4", "es": "#DB4437"},
            state_markers={"gs": "D", "es": "o"},
            cubic_images="images q = -0.10..+1.15",
        ),
    )
}
DEFAULT_DATASET = "nv_diamond"
ACTIVE = os.environ.get("PYPL_TESTING_DATASET") or DEFAULT_DATASET
if ACTIVE not in DATASETS:
    raise SystemExit("PYPL_TESTING_DATASET=%s is not a dataset; choose from %s" % (ACTIVE, ", ".join(DATASETS)))
STATE_NAMES, STATE_LABELS, STATE_COLORS, STATE_MARKERS = {}, {}, {}, {}
# "ES -> GS" / "GS -> ES" for diamond, "3E -> 1A1" / "1A1 -> 3E" for MgO.
DIRECTION_LABELS = {}


def use_dataset(name):
    """Make ``name`` the active dataset (updating the ``STATE_*`` dicts in place) and return its ``Dataset``."""

    global ACTIVE
    if name not in DATASETS:
        raise SystemExit("Unknown dataset %r; choose from %s" % (name, ", ".join(DATASETS)))
    ACTIVE = name
    ds = DATASETS[name]
    for target, source in (
        (STATE_NAMES, ds.state_names),
        (STATE_LABELS, ds.state_labels or ds.state_names),
        (STATE_COLORS, ds.state_colors),
        (STATE_MARKERS, ds.state_markers),
    ):
        target.clear()
        target.update(source)
    DIRECTION_LABELS.clear()
    DIRECTION_LABELS.update(
        {
            "emission": "%s -> %s" % (ds.state_names["es"], ds.state_names["gs"]),
            "absorption": "%s -> %s" % (ds.state_names["gs"], ds.state_names["es"]),
        }
    )
    return ds


use_dataset(ACTIVE)


def active():
    """The active ``Dataset``."""

    return DATASETS[ACTIVE]


def add_dataset_arg(parser):
    """Add ``--dataset`` to ``parser`` (default: ``PYPL_TESTING_DATASET`` or nv_diamond)."""

    parser.add_argument(
        "--dataset",
        choices=sorted(DATASETS),
        default=ACTIVE,
        help="1D CCD to run on: %s (default %s)"
        % (
            "; ".join("%s = %s, %s" % (d.name, d.title.replace("$^-$", "-"), d.transition) for d in DATASETS.values()),
            ACTIVE,
        ),
    )


# ------------------------------------------------------------------------------------- helpers


def huang_rhys(delta_q, phonon_mev):
    """Single-mode Huang-Rhys factor S = omega delta_q^2 / (2 hbar)."""

    return delta_q**2 * 1e-10**2 * AMU * phonon_mev * 1e-3 * constants.eV / constants.hbar / (2 * constants.hbar)


def normalize(energy_axis, spectrum):
    """Normalize a spectrum to unit integral over its own (possibly descending) axis."""

    return spectrum / (np.sum(spectrum) * abs(energy_axis[1] - energy_axis[0]))


def samples(ccd, state):
    """``(Q, E)`` of the DFT images of ``"gs"`` or ``"es"``: Q in amu^0.5 A, E in eV relative to the GS minimum."""

    return ccd["%s_1d_coord" % state], ccd["%s_energies" % state]


def center(ccd, state):
    """Equilibrium geometry of ``state`` on the CCD: 0 or delta_q, as the dataset places it."""

    return ccd["%s_center" % state]


def q_span(ccd, half_widths, pad=0.0):
    """``(lo, hi)``: the Q range covering ``center +- half_widths[state]`` of both wells, plus ``pad``."""

    lo = min(center(ccd, s) - half_widths[s] for s in STATES) - pad
    hi = max(center(ccd, s) + half_widths[s] for s in STATES) + pad
    return lo, hi


def image_note(ccd, state, index=None):
    """Which images of ``state`` a fit used, e.g. ``images q = -0.10..+0.45 (12 of 33)`` or ``all 33 images``.

    q is the path coordinate of the images (diamond: the interpolation coordinate, 0 = GS, 1 = ES geometry;
    MgO: q of Image_<q>, 0 = 3E, 1 = 1A1 minimum).
    """

    q = ccd["%s_path_q" % state]
    if index is None or len(np.arange(len(q))[index]) == len(q):
        return "all %d images (q = %+.2f..%+.2f)" % (len(q), q.min(), q.max())
    sel = q[index]
    return "images q = %+.2f..%+.2f (%d of %d)" % (sel.min(), sel.max(), len(sel), len(q))


def _add_state(ccd, state, q, e, fit_index, path_q, cubic_index=None):
    """Fit the parabola of ``state`` to the images ``fit_index`` and store everything about it in ``ccd``.

    ``path_q`` is the path coordinate of every image (for ``image_note``); ``cubic_index`` the images of the
    cubic fit (None: all).
    """

    potential = config_coord_1d_potential(q, e, center(ccd, state))
    a, c = potential.fit_parabola(fit_index).values()
    phonon = potential.phonon_energy()
    ccd.update(
        {
            "%s_1d_coord" % state: q,
            "%s_energies" % state: e,
            "%s_params" % state: np.array([a, c]),
            "%s_potential" % state: potential,
            "%s_phonon" % state: phonon,
            "%s_hrf" % state: huang_rhys(ccd["delta_q"], phonon),
            "%s_fit_energies" % state: potential(q),
            "n_%s_images" % state: len(q),
            "%s_path_q" % state: np.asarray(path_q, dtype=float),
            "%s_cubic_index" % state: cubic_index,
        }
    )


def build_ccd(data_dir=None, verbose=True, dataset=None):
    """Construct the 1D CCD of ``dataset`` (default: the active one) from its DFT output.

    Returns a dict with the solver inputs ``delta_q``, ``gs_phonon``, ``es_phonon`` (meV), the
    Huang-Rhys factors ``gs_hrf`` / ``es_hrf``, the parabolas as ``gs_params`` / ``es_params``
    ``(a, c)`` and as fitted ``gs_potential`` / ``es_potential`` objects, the image counts, the
    samples ``gs_1d_coord`` / ``gs_energies`` / ``gs_fit_energies`` (and the ES analogues), and
    ``dataset``, the dataset's name. ``data_dir`` overrides the dataset's own directory.
    """

    ds = DATASETS[dataset or ACTIVE]
    ccd = globals()[ds.builder](data_dir or ds.data_dir)
    ccd["dataset"] = ds.name
    if verbose:
        print("1D CCD constructed from %s (%s)" % (os.path.relpath(data_dir or ds.data_dir, REPO_ROOT), ds.name))
        print(
            "  images (%s / %s)  : %d / %d"
            % (ds.state_names["gs"], ds.state_names["es"], ccd["n_gs_images"], ccd["n_es_images"])
        )
        print("  delta_q           : % .6f amu^0.5 Angstrom" % ccd["delta_q"])
        for state in STATES:
            print(
                "  %s phonon         : % .5f meV   (HRF % .4f)"
                % (ds.state_names[state], ccd["%s_phonon" % state], ccd["%s_hrf" % state])
            )
    return ccd


def _build_nv_diamond(data_dir):
    """NV- in diamond from the tutorial 002 Quantum ESPRESSO output, exactly as the tutorial notebook."""

    atomic_symbols, gs_coord, cell = parse_atoms_qexml(os.path.join(data_dir, "gs_dft", "pwscf.xml"))
    atomic_symbols_2, es_coord, cell_2 = parse_atoms_qexml(os.path.join(data_dir, "es_cdft", "pwscf.xml"))
    assert atomic_symbols == atomic_symbols_2, "Mismatch in atomic symbols between GS and ES structures."
    assert np.max(np.abs(cell - cell_2)) < 1e-12, "Mismatch in cell parameters."

    masses = np.array([MASS_LIST[sym] for sym in atomic_symbols])
    delta_q = np.linalg.norm((es_coord - gs_coord) * masses[:, None] ** 0.5)

    # Total energies along the path; some images have no pwscf.xml and are skipped.
    ccd = {"delta_q": delta_q, "gs_center": 0.0, "es_center": delta_q}
    raw = {}
    for state in STATES:
        idx, energies = [], []
        for i in range(len(RELATIVE_COORDINATES)):
            fname = os.path.join(data_dir, "%s_1d_ccd" % state, "Image-%d" % (i + 1), "pwscf.xml")
            if os.path.exists(fname):
                idx.append(i)
                energies.append(parse_total_energy_qexml(fname))
        raw[state] = (np.array(idx), np.array(energies))
    ref_energy = np.min(raw["gs"][1])

    for state in STATES:
        idx, energies = raw[state]
        _add_state(
            ccd,
            state,
            delta_q * RELATIVE_COORDINATES[idx],
            energies - ref_energy,
            FIT_SLICES[state],
            RELATIVE_COORDINATES[idx],
        )
    return ccd


def _mgo_image_q(name):
    match = re.fullmatch(r"Image_([+-]?\d+\.\d+)", name)
    return float(match.group(1)) if match else None


def mgo_images(data_dir, state):
    """``[(q, total energy in eV)]`` of one MgO state, sorted by q: every image with a summary.json and a QE
    XML, minus ``MGO_EXCLUDE``."""

    run, root, _ = MGO_STATES[state]
    out = []
    for name in os.listdir(os.path.join(data_dir, run)):
        q = _mgo_image_q(name)
        summary = os.path.join(data_dir, run, name, "summary.json")
        xml = os.path.join(data_dir, "qe", name, "data-file-schema.xml")
        if q is None or name in MGO_EXCLUDE[state] or not (os.path.exists(summary) and os.path.exists(xml)):
            continue
        with open(summary, encoding="utf-8") as f:
            eigenval = json.load(f)["plep"]["eigenval"]
        out.append((q, parse_total_energy_qexml(xml) + eigenval[root] * RY_TO_EV))
    return sorted(out)


def _build_nv_mgo(data_dir):
    """NV- in MgO, 3E (upper, "es") and 1A1 (lower, "gs"), from the collected WEST bundle in doc/nv_mgo_west.

    The same energies as ``doc/nv_mgo_west/ccd_nv_mgo.py``: the QE total energy of each image plus the
    WEST excitation energy of the state's root. delta_q is the mass-weighted distance between
    Image_+0.00 (3E minimum) and Image_+1.00 (1A1 minimum).
    """

    xml = lambda name: os.path.join(data_dir, "qe", name, "data-file-schema.xml")  # noqa: E731
    symbols, coord_3e, cell = parse_atoms_qexml(xml("Image_+0.00"))
    symbols_2, coord_1a1, cell_2 = parse_atoms_qexml(xml("Image_+1.00"))
    assert symbols == symbols_2, "Mismatch in atomic symbols between the 3E and 1A1 structures."
    assert np.max(np.abs(cell - cell_2)) < 1e-12, "Mismatch in cell parameters."
    masses = np.array([MGO_MASS_LIST[sym] for sym in symbols])
    delta_q = float(np.linalg.norm((coord_1a1 - coord_3e) * masses[:, None] ** 0.5))

    ccd = {"delta_q": delta_q, "es_center": 0.0, "gs_center": delta_q}
    ref_energy = min(
        parse_total_energy_qexml(xml(name))
        for name in os.listdir(os.path.join(data_dir, "qe"))
        if _mgo_image_q(name) is not None and name not in MGO_ZERO_EXCLUDE
    )
    for state in STATES:
        q_path, energies = np.array(mgo_images(data_dir, state)).T
        window = lambda lo, hi: np.where((q_path >= lo - 1e-9) & (q_path <= hi + 1e-9))[0]  # noqa: E731
        _add_state(
            ccd,
            state,
            q_path * delta_q,
            energies - ref_energy,
            window(*MGO_PARABOLA_WINDOW[state]),
            q_path,
            cubic_index=window(*MGO_CUBIC_WINDOW[state]),
        )
    return ccd
