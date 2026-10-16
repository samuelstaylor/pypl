Change Log
==========

v1.1.0 (unreleased)
-------------------

Every existing function behaves exactly as in v1.0; tutorials 001 and 002 give bit-identical results
(``tests/test_legacy_tutorials.py``).

- New ``config_coord_1d_vof_solver``, a subclass of ``config_coord_1d_solver`` that adds Franck-Condon integrals by
  Hermite quadrature, by finite differences on arbitrary potentials and with a first-order cubic anharmonic
  correction, ``build_fc_lsp_from_levels`` (lines at the actual vibrational levels), stable harmonic eigenfunctions
  for any quantum number, and ``huang_rhys``.
- New ``config_coord_1d_potential``: parabola, polynomial, Morse, cubic and spline fits of a sampled PES (spline
  degree 2-5; natural, not-a-knot or clamped ends; harmonic or smooth C2 tails).
- New ``pypl.isc``: golden-rule intersystem crossing rate and lifetime.
- New ``utils.parse_west_summary`` for WEST TDDFT excitation energies.
- New tutorial ``pypl_003``: 3E -> 1A1 ISC of NV- in MgO from WEST TDDFT-DDH data (``doc/nv_mgo_west``).
- pytest suite (``tests/``), validation suite for the 1D CCD (``validation/``), CI test job, ``pypl.__version__``.
- Docstring LaTeX markup fixed; invalid escape sequences in two printed strings made raw (same output).

v1.0 (2025/09/05)
-----------------

- Initial stable release of the PyPL code.
