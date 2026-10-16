"""Regenerate ``tests/data/legacy_reference.npz`` from the ORIGINAL PyPL code.

The reference pins the behaviour of the original API (tutorials 001 and 002). Generate it from the untouched
package, not from the working tree, e.g. with a worktree of the ``pypl-original`` branch::

    git worktree add ../pypl_original_wt pypl-original
    python tests/tools/make_legacy_reference.py --pypl-root ../pypl_original_wt

Regenerating it from modified code would defeat its purpose. Do it only when a behaviour change of the original
API has been deliberately approved, and say so in the commit message.
"""

import argparse
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
TESTS = os.path.dirname(HERE)
DEFAULT_OUT = os.path.join(TESTS, "data", "legacy_reference.npz")


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--pypl-root", required=True, help="checkout whose pypl/ package is imported (the original code)")
    p.add_argument("--out", default=DEFAULT_OUT)
    args = p.parse_args()

    root = os.path.abspath(args.pypl_root)
    sys.path.insert(0, root)
    sys.path.insert(1, TESTS)
    import numpy as np
    import pypl.config_coord_1d_solver

    loaded = os.path.dirname(os.path.dirname(os.path.abspath(pypl.config_coord_1d_solver.__file__)))
    if os.path.normcase(loaded) != os.path.normcase(root):
        sys.exit("pypl was imported from %s, not from --pypl-root %s" % (loaded, root))

    from legacy_cases import all_cases

    try:
        commit = subprocess.check_output(["git", "-C", root, "rev-parse", "HEAD"], text=True).strip()
    except (OSError, subprocess.CalledProcessError):
        commit = "unknown"
    cases = all_cases()
    cases["_source_commit"] = np.array(commit)
    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    np.savez_compressed(args.out, **cases)
    print("wrote %s: %d arrays from pypl at %s (commit %s)" % (args.out, len(cases) - 1, root, commit))


if __name__ == "__main__":
    main()
