.. _validation:

======================
Validation and testing
======================

This page describes everything that checks **PyPL**: what is tested, why, against what reference and with which
tolerance, how to run it and how to read the results. The formulas and derivations behind every method are collected
in :ref:`supporting_information`. A hands-on walk-through is Tutorial 005 (``doc/tutorials/pypl_005.ipynb``).

.. contents:: Contents
   :local:
   :depth: 2


Overview
========

There are three layers of checks. The first two are fast and run in continuous integration (CI) on every push; the
third produces documented comparisons of methods and fits.

.. list-table::
   :header-rows: 1
   :widths: 18 40 27 15

   * - Layer
     - What it is
     - How to run it
     - Time
   * - Test suite, ``tests/``
     - Unit and regression tests of the package (pytest), including the guarantee that the original API reproduces
       tutorials 001 and 002 bit for bit.
     - ``pytest``
     - ~1 min
   * - Validation checks
     - 81 physics and regression checks of the 1D configuration-coordinate (1D CCD) path on two datasets.
     - ``python validation/check_1d_ccd.py [-v]``
     - ~15 s
   * - Validation suite, ``validation/``
     - Scripts that compare every route to the Franck-Condon integrals and every fit of the potential energy
       surfaces, with figures, CSV tables and an ``index.md`` per run.
     - ``compare_vof.py``, ``plot_ccd_wavefunctions.py``, ``vof_report.py``, ``run_all.py``
     - 30 s to 20 min

CI (``.github/workflows/ci.yml``) runs the black formatting check on ``pypl``, ``tests`` and ``validation``, builds the
package, and runs ``pytest -W error::SyntaxWarning`` followed by ``validation/check_1d_ccd.py`` on Python 3.9, 3.11
and 3.12. Either failing fails the job.


Terminology
===========

**Image.** One configuration along the 1D configuration coordinate at which the total energy of an electronic state
was computed. The images of a 1D CCD are usually a linear interpolation between the equilibrium structures of the two
states, labelled by :math:`q` with :math:`Q = q\,\Delta Q`.

**Computed energies, not "DFT energies".** The data points of a potential energy surface (PES) are total energies of
one electronic state at each image, obtained by *some* electronic-structure method. That method differs between
datasets and even between the two states of one dataset:

.. list-table::
   :header-rows: 1
   :widths: 30 35 35

   * - Dataset
     - Lower state ("GS")
     - Upper state ("ES")
   * - NV\ :sup:`-` in diamond (tutorial 002)
     - :math:`^3A_2`: ground-state DFT (PBE)
     - :math:`^3E`: :math:`\Delta`\ SCF, constrained occupations (PBE)
   * - NV\ :sup:`-` in MgO (tutorial 003)
     - :math:`^1A_1`: spin-flip TDDFT on a DDH ground state (WEST, root 4)
     - :math:`^3E`: spin-conserving TDDFT on a DDH ground state (WEST, root 1)

Other systems may use the Bethe-Salpeter equation (BSE), quantum embedding (QDET), quantum chemistry or a machine
learned potential. "DFT" would therefore be wrong or misleading for most excited-state PESs: TDDFT and
:math:`\Delta`\ SCF are built on DFT, but the energy of the excited state comes from a different approximation with
different errors, and many-body methods use DFT only as a starting point or not at all. The package (which does not
know where the numbers come from) calls them *samples* or *images*; the validation suite labels them *computed
energies* and names the method of each state from the dataset's metadata (``Dataset.energy_methods`` in
``validation/harness/ccd_data.py``). The energy zero is recorded too (``Dataset.energy_zero``).

**GS / ES roles.** The validation code calls the lower state of a dataset "GS" and the upper one "ES", whatever their
physical names; *emission* is ES :math:`\rightarrow` GS and *absorption* GS :math:`\rightarrow` ES. Figures and logs
use the physical names (e.g. :math:`^3E \rightarrow {}^1A_1`).

**Franck-Condon matrix and VOF.** :math:`F_{ij} = \langle \chi_i | \chi_j \rangle` is the overlap of vibrational level
:math:`i` of the initial state with level :math:`j` of the final state. The vibrational overlap function (VOF, the
Franck-Condon lineshape) is :math:`I(E) = \sum_{ij} p_i |F_{ij}|^2\,G(E - E_{ij})`, see :ref:`si_lineshape`.

**Method.** A route to :math:`F_{ij}`: ``pypl`` (Ruhoff recurrences), ``hermite`` (quadrature of harmonic
eigenfunctions), ``fd`` (finite-difference eigenstates of any potential) or ``cubic`` (first-order cubic perturbation
theory). Every method fills the same ``fc_ints``, so the same lineshape code builds the VOF from any of them.

**Fit.** A model of one PES made from its images: parabola, polynomial, Morse, cubic or spline (see
:ref:`si_fits`).

**Like-for-like vs anharmonic comparisons.** On the same harmonic wells ``pypl``, ``hermite`` and ``fd`` solve the same
problem by independent routes, so their difference is numerical error and must be below tight tolerances (``PASS`` /
``OK``). An anharmonic fit or method compared with a harmonic reference differs because the physics differs; that
difference is reported as the anharmonic *effect* (``DIFFERS (expected)``), not as an error.


Datasets
========

.. list-table::
   :header-rows: 1
   :widths: 22 39 39

   * -
     - ``nv_diamond`` (default)
     - ``nv_mgo``
   * - System, transition
     - NV\ :sup:`-` in diamond, :math:`^3E \rightarrow {}^3A_2` (PL)
     - NV\ :sup:`-` in MgO, :math:`^3E \rightarrow {}^1A_1` (intersystem crossing)
   * - Data
     - ``doc/tutorials/002_nv_diamond_1d_ccd`` (Quantum ESPRESSO)
     - ``doc/tutorials/003_nv_mgo_isc`` (Quantum ESPRESSO + WEST)
   * - Energies
     - GS: DFT (PBE); ES: :math:`\Delta`\ SCF (PBE)
     - ground-state DDH (PBE0 with 34 % exact exchange) total energy + TDDFT excitation energy
       (:math:`^3E`: spin-conserving root 1; :math:`^1A_1`: spin-flip root 4)
   * - Energy zero
     - GS minimum image
     - :math:`^3A_2` minimum (DDH), excluding ``Image_-1.00``
   * - Minima
     - GS at :math:`Q = 0`, ES at :math:`\Delta Q`
     - :math:`^3E` at :math:`Q = 0` (red), :math:`^1A_1` at :math:`\Delta Q` (blue)
   * - :math:`\Delta Q` (amu\ :sup:`1/2` Å)
     - 0.695732
     - 3.097879
   * - Images
     - GS 15, ES 11 (:math:`q = -0.2 \dots 1.2`)
     - :math:`^3E` 33, :math:`^1A_1` 30 (:math:`q = -1.0 \dots 2.0`; ``sf/Image_+2.00`` is not a spin-flip run and is
       not used)
   * - Parabola images
     - GS 1-5, ES 5-11
     - :math:`^3E` :math:`q = -0.10 \dots 0.45`, :math:`^1A_1` :math:`q = 0.60 \dots 1.15` (12 each)
   * - Cubic-fit images
     - all
     - :math:`q = -0.10 \dots 1.15` (26 each)
   * - Other fits (splines, Morse, polynomials)
     - all images
     - all images
   * - Harmonic :math:`\hbar\omega` / Huang-Rhys :math:`S`
     - GS 59.089 meV / 3.42; ES 64.312 meV / 3.72
     - :math:`^1A_1` 31.551 meV / 36.2; :math:`^3E` 60.833 meV / 69.8
   * - Basis (ES :math:`\times` GS levels)
     - 50 :math:`\times` 60
     - 160 :math:`\times` 100
   * - Reference method
     - ``pypl``
     - ``hermite`` (the Ruhoff recurrence is unstable at large :math:`S`, finding 6 below)

A dataset is selected with ``--dataset`` or the environment variable ``PYPL_VALIDATION_DATASET``. Every log, figure
caption and ``index.md`` records the dataset, the state roles, the method behind each state's energies and the images
each fit used.


The test suite (``tests/``)
===========================

Run with ``pytest`` from the repository root (configured in ``pyproject.toml``); install the extras first with
``pip install ".[test]"``. ``pytest -m "not slow"`` skips the end-to-end tutorial-003 test. The fixtures in
``tests/conftest.py`` build both CCDs with the package API only.

``test_legacy_tutorials.py``: the original API is unchanged
-----------------------------------------------------------

``tests/legacy_cases.py`` reruns the computations of tutorials 001 and 002 as plain functions and returns 35 arrays.
``tests/data/legacy_reference.npz`` holds the same arrays computed by the **original** PyPL code (branch
``pypl-original``, commit recorded in the file), generated with
``python tests/tools/make_legacy_reference.py --pypl-root <checkout of the original code>``.

.. list-table::
   :header-rows: 1
   :widths: 35 65

   * - Test
     - What it checks
   * - ``test_tutorial_array_unchanged[key]`` (35)
     - Each array equals the reference to ``rtol = 1e-10`` (bit-identical on the machine that made the reference; the
       tolerance only absorbs FFT / BLAS rounding of other platforms). Tutorial 001: phonon frequencies, coordinates,
       cell, forces; partial Huang-Rhys factors from displacements and forces (PL) and for absorption; spectral
       densities at three broadenings; FFT lineshape and its axis; time-integration lineshape; PL spectra at 4, 150
       and 300 K; absorption spectrum. Tutorial 002: :math:`\Delta Q`, image energies, both phonon energies,
       :math:`F_{ij}`, lineshapes and spectra for PL and absorption, and a 300 K lineshape with a Gaussian ZPL.
   * - ``test_same_quantities``, ``test_reference_comes_from_the_original_code``
     - The reference holds exactly the computed quantities and records its source commit.
   * - ``test_notebook_001_total_huang_rhys``, ``test_notebook_001_largest_modes``,
       ``test_notebook_002_printed_values``
     - The numbers printed in the stored outputs of the two notebooks (total Huang-Rhys factors, the ten largest
       partial factors with their frequencies, :math:`\Delta Q`, phonon energies, Huang-Rhys factors), to their printed
       precision (the notebooks were run with other numpy / scipy versions).

**A failure here means the behaviour of the original API changed.** Never regenerate the reference from modified
code to make it pass.

``test_config_coord_1d_solver.py``: physics of the Ruhoff solver
----------------------------------------------------------------

.. list-table::
   :header-rows: 1
   :widths: 35 50 15

   * - Test
     - Assertion
     - Tolerance
   * - ``test_closure``
     - :math:`\sum_j |F_{ij}|^2 = 1` for :math:`i < 10` with 300 final levels; :math:`0 < |F_{00}|^2 < 1`
     - 1e-8 (row 0), 1e-6
   * - ``test_identity_at_zero_displacement``
     - :math:`\Delta Q = 0`, equal frequencies: :math:`F = I`
     - 1e-10
   * - ``test_equal_frequency_symmetry``
     - equal frequencies: :math:`F_{ij} = (-1)^{i+j} F_{ji}` (:ref:`si_fc`)
     - 1e-9
   * - ``test_lineshape_normalization_and_populations``
     - lineshape integral :math:`\approx 1`, non-negative; Boltzmann weights sum to 1, :math:`p_0 = 1` at 5 K,
       :math:`p_1 > 10^{-3}` at 800 K
     - 0.05, 1e-12
   * - ``test_gaussian_zpl_integrates_to_one``
     - with a Gaussian ZPL the lineshape integrates to 1 (the 1.01 of the Lorentzian ZPL is grid sampling)
     - 1e-4
   * - ``test_basis_convergence``
     - 50 :math:`\times` 60 levels converged against 80 :math:`\times` 100
     - 1e-4 of the peak
   * - ``test_spectra``
     - PL axis descends, Abs axis ascends; spectra non-negative; both peak at the ZPL; PL weight below and Abs weight
       above the ZPL; an invalid ``spectrum_type`` raises
     - 1 meV

``test_config_coord_1d_vof_solver.py``: the extended solver
-----------------------------------------------------------

.. list-table::
   :header-rows: 1
   :widths: 35 50 15

   * - Test
     - Assertion
     - Tolerance
   * - ``test_subclass_overrides_nothing``
     - ``config_coord_1d_vof_solver`` is a subclass and redefines no attribute of the parent
     - exact
   * - ``test_inherited_methods_bit_identical``
     - inherited ``compute_franck_condon_integrals``, ``bulid_fc_lsp``, ``compute_spectrum`` give identical arrays on
       the subclass and the parent
     - ``array_equal``
   * - ``test_recurrence_matches_explicit_hermite_formula``
     - normalized Hermite-function recurrence = explicit :math:`H_n` formula, :math:`n \le 20`
     - 1e-12 rel.
   * - ``test_harmonic_functions_orthonormal_at_large_n``
     - Gram matrix of :math:`\chi_0 \dots \chi_{300}`
     - 1e-9
   * - ``test_huang_rhys_and_turning_point``, ``test_kinetic_prefactor``
     - :math:`S = \omega \Delta Q^2/2\hbar` in SI; :math:`\alpha x_n^2 = 2n + 1`; :math:`\hbar^2/2` in eV amu Å\ :sup:`2`
     - 1e-12 rel.
   * - ``test_hermite_matches_ruhoff``
     - Hermite quadrature = Ruhoff, 20 :math:`\times` 20
     - 1e-10
   * - ``test_hermite_stable_for_large_huang_rhys``
     - 200 final levels at :math:`S \approx 50`: finite, row 0 closes
     - 1e-10
   * - ``test_fd_on_parabolas_matches_ruhoff``
     - finite differences on the parabolas = Ruhoff; walls, orthonormality and ladder diagnostics
     - 1e-8; 1e-12; 1e-10; 1e-4 meV
   * - ``test_fd_three_point_stencil_is_less_accurate``
     - the 3-point stencil error lies between 1e-8 and 1e-3
     - --
   * - ``test_fd_is_orientation_independent``
     - mirroring the CCD (initial well on the left) leaves :math:`F_{ij}` unchanged
     - 1e-7
   * - ``test_solve_schrodinger_fd_harmonic_ladder``
     - eigenvalues :math:`(n + \tfrac12)\hbar\omega`; orthonormal; signed like the Hermite functions
     - 1e-6 meV; 1e-10; 1e-8
   * - ``test_anharmonic_strength_matches_si``, ``test_cubic_coefficients_structure``
     - :math:`\xi` against an all-SI evaluation; the perturbation coefficients (:ref:`si_cubic`)
     - 1e-12 rel.
   * - ``test_cubic_with_b_zero_is_hermite``
     - :math:`b = 0` reproduces the Hermite route
     - 1e-13
   * - ``test_cubic_norms_and_pt_order``
     - norms = :math:`1 + \sum c^2`; PT error vs the exact cubic falls as :math:`b^2` (ratio 3.5-4.5 when :math:`b`
       halves); the wrong frame sign of :math:`b` is > 10 times worse
     - 1e-12
   * - ``test_cubic_warning``
     - :math:`|\xi| > 0.1` prints the warning
     - --
   * - ``test_levels_lineshape_equals_bulid_fc_lsp_for_harmonic_ladders``
     - ``build_fc_lsp_from_levels`` with harmonic ladders = ``bulid_fc_lsp`` at 300 K
     - 1e-14 of the peak
   * - ``test_levels_lineshape_uses_fd_levels``, ``test_levels_lineshape_errors``
     - line positions are the finite-difference level differences; missing / too few levels raise
     - 1e-7 rel.

``test_config_coord_1d_potential.py``: fits
-------------------------------------------

Parabola reproduces tutorial 002 (59.08890 / 64.31191 meV); an unfitted potential raises; polynomials (even degree
with positive leading coefficient bounded, odd degree and negative leading coefficient flagged); Morse (bounded,
rms < 1 meV, :math:`\hbar\omega` within 10 % of harmonic); cubic (:math:`a` kept, turnover
:math:`Q_0 - 2a/3b`, flagged unbounded); default spline = scipy's natural cubic spline (1e-12 eV); six order / end
condition combinations interpolate every image (1e-12 eV); clamped junctions C1 (cubic) and C2 (quintic); smoothing
spline meets its rms target; eight invalid option sets raise ``ValueError``; smooth tails are identical to the plain
spline inside the data, C2 at both ends, bounded, and equal to the second-order Taylor expansion beyond the data;
negative end curvature is flagged unbounded; the spline minimum and :math:`\hbar\omega` are sensible.

``test_isc.py`` and ``test_nv_mgo.py``
--------------------------------------

The ISC rate equals :math:`4\pi/\hbar\,|\lambda|^2 F` (1e-14), the units GHz / meV / µeV agree, an unknown unit
raises, arrays and lifetimes (:math:`\tau = 1/\Gamma`, infinite for :math:`\Gamma = 0`) behave. The WEST parser
converts Ry to eV and rejects a file without ``plep/eigenval``. The MgO CCD reproduces :math:`\Delta Q = 3.097879`,
the image counts, the minima, :math:`E(^3E) = 3.1376` eV, the adiabatic gap 554.4 meV and the parabola phonon
energies 60.8332 / 31.5511 meV; row 0 of the Ruhoff recurrence equals Hermite at :math:`S \approx 36`-70 (1e-12);
fd = Hermite in the MgO orientation (1e-7); and (``slow``) the headline numbers of tutorial 003: :math:`F(\Delta)`
2.70e-7 / 3.97e-5 meV\ :sup:`-1`, :math:`\Gamma` 3.711 / 545.1 MHz, peaks 1145 / 660 meV (0.5 %).


The validation checks (``validation/check_1d_ccd.py``)
======================================================

The script prints one header per group and a summary line; ``-v`` prints every check with its measured value. The exit
status is 0 only if all pass. Reference values are those of tutorial 002 and of the user's MgO analysis script;
``REFERENCE`` is changed only when a change of physics is intended, and tolerances are never widened to make a check
pass.

.. list-table::
   :header-rows: 1
   :widths: 24 56 20

   * - Group (checks)
     - What is asserted
     - Why
   * - CCD construction (6)
     - :math:`\Delta Q`, both :math:`\hbar\omega` and Huang-Rhys factors equal tutorial 002 (1e-5 / 1e-3); ES frequency
       above GS
     - the parsers and fits reproduce the tutorial
   * - Franck-Condon integrals (5)
     - closure of rows 0 (1e-8) and :math:`i < 10` (1e-6); all finite; :math:`0 < |F_{00}|^2 < 1`; identity at
       :math:`\Delta Q = 0`
     - exact properties of overlaps between complete sets
   * - FC lineshape (7)
     - integral, non-negativity, Boltzmann weights, populations at 5 and 800 K, convergence 50 :math:`\times` 60 vs
       80 :math:`\times` 100, :math:`F_{ij} = (-1)^{i+j}F_{ji}` and identical PL / Abs lineshapes at equal frequencies
     - catches index slips in the recursion and in the thermal sum
   * - PL and absorption spectra (9)
     - axis directions, non-negativity, peaks at the ZPL, sideband on the correct side, invalid type raises
     - the spectrum conventions
   * - Extension leaves PyPL unchanged (5)
     - no attribute overridden; inherited methods bit-identical
     - the extension contract
   * - New :math:`F_{ij}` routes (3)
     - Hermite and fd on the parabolas equal Ruhoff (1e-10, 1e-8); fd states vanish at the walls (1e-12)
     - three independent routes agree
   * - Potential fits (3), spline options (11)
     - splines interpolate; quartic flagged unbounded; parabola :math:`\hbar\omega` = solver's; natural cubic = scipy;
       every order / end condition interpolates; clamped junctions C1 / C2; smoothing target; even-order natural /
       clamped rejected
     - the fits are what their names say
   * - Cubic anharmonic correction (9)
     - :math:`a` kept; :math:`\hbar\omega` harmonic; :math:`\xi` vs SI; :math:`b = 0` = Hermite (1e-13); norms;
       :math:`1 + O(\xi^2)`; error :math:`\propto b^2`; frame sign of :math:`b`; warning
     - first-order perturbation theory is implemented correctly
   * - Large Huang-Rhys factors and NV\ :sup:`-` in MgO (23)
     - Hermite-function recurrence = explicit formula; orthonormal to :math:`n = 300`; 200 levels finite and closed at
       :math:`S \approx 50`; MgO :math:`\Delta Q`, images, minima, :math:`E(^3E)`, ZPL; image windows of every fit;
       smooth tails (inside = plain spline, C2, bounded, clamped + smooth rejected, canonical name);
       ``build_fc_lsp_from_levels`` = ``bulid_fc_lsp`` for harmonic ladders; fd = Hermite in the MgO orientation
     - large-:math:`S` stability and the MgO data path


The validation suite (``validation/``)
======================================

Scripts
-------

.. list-table::
   :header-rows: 1
   :widths: 25 45 30

   * - Script
     - What it does
     - Output directory
   * - ``run_1d_ccd.py``
     - Tutorial 002 as a script with the unmodified solver: CCD, PL and absorption.
     - ``run_1d_ccd/<tag>/``
   * - ``compare_vof.py``
     - ``--reference`` method against each ``--test`` method, both directions. PASS / FAIL with tolerances
       ``--tol-fc`` (max :math:`|\Delta F_{ij}|`, default 1e-7) and ``--tol-vof`` (max :math:`|\Delta \mathrm{VOF}|`
       / peak, default 1e-8). Exit 0 only if all pass.
     - ``compare_vof/<ref>_vs_<tests>[_fit-X]/``
   * - ``plot_ccd_wavefunctions.py``
     - The CCD with the computed images, the fitted wells and the lowest wavefunctions (``--wf-method``), and a VOF
       panel (``--vof-method``). The drawn wavefunctions are checked against the solver's :math:`F_{ij}` before
       plotting.
     - ``plot_ccd_wavefunctions/wf-<m>_vof-<ms>[_fit-X]/``
   * - ``vof_report.py``
     - Seven sections: (1) potential fits, (2) the finite-difference grid and its convergence, (3) wavefunctions,
       (4) :math:`F_{ij}` error maps, (5) VOF by method, (6) VOF by fit, (7) cubic perturbation theory.
     - ``vof_report/<tag>/``
   * - ``run_all.py``
     - Every method, fit and variant from one TOML configuration (``run_all.toml``, ``run_all_nv_mgo.toml``); each
       step is a subprocess of one of the scripts above.
     - ``suite/<name>/``

Configuration
-------------

Every setting can be changed without editing the source, either on the command line of a script or in the TOML
configuration of ``run_all.py`` (``validation/run_all.toml`` documents every key; ``--set section.key=value``
overrides one key). The defaults reproduce tutorial 002.

.. list-table::
   :header-rows: 1
   :widths: 26 34 40

   * - Setting
     - Script flag
     - ``run_all.toml`` key
   * - dataset
     - ``--dataset``
     - ``suite.dataset``
   * - temperature, broadening, ZPL
     - ``--temp``, ``--sigma``, ``--gamma``, ``--no-zpl``
     - ``physics.temp``, ``.sigma``, ``.gamma``, ``.zpl_lorentzian``
   * - basis, energy axis
     - ``--order-es``, ``--order-gs``, ``--ene-range``, ``--resol``
     - ``physics.order_es``, ``.order_gs``, ``.ene_range``, ``.resol``
   * - line positions of the VOF
     - ``--line-positions {harmonic,levels}``
     - ``physics.line_positions``
   * - quadrature grid (hermite, cubic)
     - ``--npts`` or ``--hermite-h``; ``--hermite-margin``
     - ``hermite.npts``, ``.h``, ``.margin``
   * - finite-difference grid
     - ``--fd-npts`` or ``--fd-h``; ``--fd-range QMIN QMAX`` or ``--fd-margin``; ``--fd-stencil {3,5,7}``;
       ``--fd-shift``
     - ``fd.npts``, ``.h``, ``.range``, ``.margin``, ``.stencil``, ``.shift``
   * - fits
     - ``--fd-fit SPEC``, ``--spline-*``
     - ``fd.fits``, ``[[fd.splines]]``
   * - cubic perturbation theory
     - ``--cubic-points``, ``--cubic-free-a``, ``--cubic-normalize``, ``--cubic-warn``
     - ``[[cubic.variants]]``, ``cubic.warn``
   * - pass criteria of comparisons
     - ``--tol-fc``, ``--tol-vof``
     - ``tolerances.fc``, ``tolerances.vof``
   * - pass criteria of the checks
     - ``check_1d_ccd.py --tol NAME=VALUE`` (``--list-tolerances``)
     - ``[checks.tolerances]``
   * - wavefunction check of the plots
     - ``plot_ccd_wavefunctions.py --check-tol ANALYTIC FD``
     - ``plots.check_tol``
   * - report contents
     - ``vof_report.py --fits``, ``--skip-*`` and the level / convergence / cubic-scan options (``--help``)
     - ``report.fits``, ``.convergence``, ``.cubic``, ``.options``
   * - figures
     - ``--dpi``, ``--fig-format png pdf svg``, ``--font-size``
     - ``figures.dpi``, ``.formats``, ``.font_size``
   * - output location
     - environment ``PYPL_VALIDATION_OUTPUT_ROOT``
     - ``suite.name``

Figures follow one style (``validation/harness/plots.py``): 9 pt sans-serif text with math-text symbols and upright
units, the colour-blind-safe Okabe-Ito palette for methods and fits (all spline variants share one colour and differ
by line style; the two states keep the dataset's colours), one shared legend per figure below the panels, viridis
for every colour map with a labelled colour bar, and 300 dpi by default. Settings and statistics are written to
``log.txt``, the CSV files and ``index.md`` rather than onto the figures.

**Fit specifications** for ``--fd-fit``: ``parabola``, ``polyN``, ``morse``, ``cubic[-near][-free]`` and
``spline[K][-natural|-notaknot|-clamped|-sRMS][-smooth]`` (``K`` = degree 2-5; ``-sRMS`` a smoothing spline with an
rms target in meV; ``-smooth`` C2 tails). Unbounded fits are refused, with the reason.

Comparison metrics
------------------

For a test result against a reference (same direction, same energy axis):

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Metric
     - Definition
   * - max :math:`|\Delta F_{ij}|`, at :math:`(i, j)`
     - largest absolute difference over the whole matrix, and where it occurs
   * - max :math:`|\Delta F_{ij}|`, :math:`i, j < 30`
     - the same over the low levels, which carry the intensity at low temperature
   * - max / rms :math:`|\Delta \mathrm{VOF}|` / peak
     - largest / rms difference of the lineshapes, divided by the reference peak
   * - centroid shift
     - :math:`\sum_E E\,I(E) / \sum_E I(E)` of the test minus that of the reference
   * - closure, integral
     - :math:`\sum_j |F_{0j}|^2` of the test; :math:`\int I(E)\,dE`
   * - verdict
     - ``OK`` if both maxima are within the tolerances, else ``FAIL``

``run_all.py`` gives each step a status: ``PASS`` (checks, comparisons, report within tolerance), ``OK`` (a plot or
script finished), ``DIFFERS (expected)`` (an anharmonic fit or method that differs from the harmonic reference, as it
must), ``FAIL`` (a like-for-like comparison out of tolerance) or ``ERROR`` (a crash). Its ``index.md`` lists the
dataset, the skipped fits with the reason, the comparison table, every step with a link and every command.

Outputs
-------

Every run writes ``validation/outputs/<script>/<run>/`` (or ``validation/outputs/suite/<name>/<script>/<run>/``):
``index.md`` (every file with a description, and the headline tables), ``log.txt`` (every setting and result),
``settings.json``, ``data/*.csv`` and ``figures/*.png``. A rerun under the same name **archives** the previous run to
``validation/outputs/_archive/<same path>/<timestamp>/``; nothing is ever deleted. ``PYPL_VALIDATION_OUTPUT_ROOT``
redirects all outputs (Tutorial 005 uses a temporary directory).


Reference results
=================

NV\ :sup:`-` in diamond, parabola wells, 5 K, 50 :math:`\times` 60 levels:

.. list-table::
   :header-rows: 1
   :widths: 46 27 27

   * - Comparison
     - max :math:`|\Delta F_{ij}|`
     - max :math:`|\Delta\mathrm{VOF}|` / peak
   * - hermite vs pypl
     - 1.7e-8 (5e-11 for :math:`i, j < 30`)
     - 7e-16
   * - fd (5-point) vs pypl
     - 1.8e-8 (6e-10)
     - 1.3e-12
   * - fd vs hermite
     - 3.5e-9
     - ~1e-12
   * - fd (3-point) vs pypl
     - 5.5e-5
     - 1.2e-6
   * - cubic PT vs exact cubic (fd)
     - 6e-3 (:math:`i, j < 10`)
     - 3.1e-4 (99 % of the 3.4e-2 cubic effect)
   * - fd on the smooth not-a-knot spline vs pypl (anharmonic effect)
     - --
     - 3.2e-3

The finite-difference errors fall as :math:`h^2` (3-point) and :math:`h^4` (5-point): measured log-log slopes 2.00 and
4.00. The 1.7e-8 residual of ``pypl`` sits at the highest indices and is the float64 rounding of the recurrence.

NV\ :sup:`-` in MgO, 160 :math:`\times` 100 levels, constant 20 meV broadening, ``hermite`` reference: fd vs hermite
1.0e-8 in :math:`F_{ij}` and 2e-9 in the VOF; pypl vs hermite 1e-14 in the 5 K VOF but up to 1e27 in :math:`F_{ij}`
(finding 6); the smooth-spline VOF peaks about 580 meV below the harmonic one (Tutorial 003).


Findings about the original solver
==================================

None of these is changed in the package (changing them would change PyPL's behaviour); they are documented and, where
relevant, avoided by the new methods.

1. **Float64 drift of the Ruhoff recurrence** (:math:`\sim 10^{-8}` at the highest indices of a 50 :math:`\times` 60
   basis); the same recurrence in 60-digit arithmetic agrees with Hermite quadrature to 1e-12. Invisible in spectra.
2. **The Lorentzian ZPL ignores the thermal weight**: ``bulid_fc_lsp`` replaces :math:`|F_{00}|^2` without
   :math:`p_0`; harmless at 5 K, about 9 % too much ZPL weight at 300 K.
3. **"Integral check: 1.0101"** is the 0.3 meV Lorentzian sampled on a 1 meV grid, not a normalization error.
4. The lineshape method is spelled ``bulid_fc_lsp``.
5. ``bulid_fc_lsp`` always places lines on the harmonic ladders; ``build_fc_lsp_from_levels`` uses actual levels.
6. **The Ruhoff recurrence is unstable at large Huang-Rhys factors**: on NV\ :sup:`-` in MgO row 0 is exact but the
   error grows about tenfold per initial level, exceeding 1 from :math:`i \approx 22`-27 (an overlap cannot exceed 1).
   The 5 K VOF is unaffected (only :math:`i = 0` is populated); use ``hermite`` for thermal or closure work.
7. **The Gaussian width grows linearly with** :math:`|E|` without a cap (:ref:`si_lineshape`); at eV-scale sidebands
   it reaches hundreds of meV, so the MgO suite uses a constant width.


Extending the validation
========================

- **A new method:** add a ``compute_franck_condon_integrals_*`` method to ``config_coord_1d_vof_solver`` (it must fill
  ``fc_ints``, ``ni``, ``nf`` and must not override anything), register it in ``validation/harness/methods.py``
  (``METHODS``), add tests in ``tests/`` and checks in ``check_1d_ccd.py``, and document it here, in
  :ref:`supporting_information` and in ``validation/METHODS.md``.
- **A new fit:** add a ``fit_*`` method to ``config_coord_1d_potential`` (set ``kind``, ``index``, ``params``,
  ``bounded``), add its spec to ``validation/harness/fits.py``, and test it.
- **A new dataset:** add a ``Dataset`` entry with its builder to ``validation/harness/ccd_data.py``, including
  ``energy_methods`` and ``energy_zero``, and a suite configuration ``validation/run_all_<name>.toml``.
