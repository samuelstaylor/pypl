"""Test-harness helpers for the 1D CCD scripts in ``testing/``. Not part of PyPL and never moved into it.

    ccd_data   the NV-center 1D CCD from tutorial 002's Quantum ESPRESSO output
    fits       the registry of potential fits used by the finite-difference method
    methods    the registry of F_ij methods, VOFParams, the shared CLI flags, compute_vof()
    output     OutputDir: one organized directory per run, with a log, data/ and figures/
    plots      shared plot style and figure builders

Importing the package puts the repository root on ``sys.path`` so ``import pypl`` picks up the
working tree without ``pip install``, and ``testing/`` so ``import pypl_ext`` resolves.
"""

import os
import sys

TESTING_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
REPO_ROOT = os.path.dirname(TESTING_DIR)
for _path in (REPO_ROOT, TESTING_DIR):
    if _path not in sys.path:
        sys.path.insert(0, _path)
