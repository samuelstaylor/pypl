.. _supporting_information:

======================
Supporting information
======================

Formulas and derivations behind every method of **PyPL** and of its validation framework (:ref:`validation`), in the
order a calculation uses them. Each section names the function that implements it. Symbols: :math:`Q` mass-weighted
configuration coordinate, :math:`\Omega` angular frequency, :math:`\hbar\omega` phonon energy, :math:`S` Huang-Rhys
factor, :math:`F_{ij}` Franck-Condon overlap, :math:`I(E)` lineshape (vibrational overlap function, VOF).

.. contents:: Contents
   :local:
   :depth: 1


.. _si_units:

S1. Units and conventions
=========================

The user-facing units are those of the 1D CCD literature: :math:`Q` in amu\ :sup:`1/2` Å, energies in eV (potentials)
or meV (phonons, lineshapes), temperatures in K. In the mass-weighted coordinate a harmonic well is

.. math::

   V(Q) = \tfrac12 \Omega^2 Q^2 = a\,Q^2, \qquad \Omega = \sqrt{2a}, \qquad \hbar\omega = \hbar\sqrt{2a},

with :math:`a` in eV/(amu Å\ :sup:`2`). Converted to SI, :math:`\Omega^2 = 2a\,e/(10^{-20}\,m_u)` (rad/s)\ :sup:`2`, with
:math:`e` the elementary charge and :math:`m_u` the atomic mass constant; ``config_coord_1d_potential.phonon_energy``
evaluates exactly the expression of tutorial 002, so the frequencies are bit-identical to it. The kinetic prefactor of
the Schrödinger equation in these units is

.. math::

   K = \frac{\hbar^2}{2} = \frac{\hbar^2}{2\, m_u\, 10^{-20}\,\mathrm{m^2}\; e} \ \text{eV amu Å}^2
   \approx 2.0904\times10^{-3}\ \text{eV amu Å}^2

(``KINETIC``). The inverse squared oscillator length is :math:`\alpha = \Omega/\hbar`
(``harmonic_alpha``), so :math:`\alpha = \hbar\omega / (2K)` and the zero-point length is
:math:`\ell = \sqrt{\hbar/2\Omega} = 1/\sqrt{2\alpha}`.

**Roles and frame.** The solver's frame puts the minimum of the final state at 0 and that of the initial state at
:math:`-\Delta Q`, i.e. :math:`F_{ij} = \langle \chi_i(Q + \Delta Q) | \chi_j(Q)\rangle`. A CCD with the minima the
other way round is mapped onto it by the mirror :math:`x = -(Q - Q_f)`; odd terms of a potential (e.g. :math:`bQ^3`)
change sign under it.


.. _si_ccd:

S2. The one-dimensional configuration coordinate
=================================================

For two equilibrium structures :math:`\{\mathbf R^{(1)}_\alpha\}`, :math:`\{\mathbf R^{(2)}_\alpha\}` with atomic masses
:math:`M_\alpha`,

.. math::

   \Delta Q = \Big[ \sum_\alpha M_\alpha\, |\mathbf R^{(2)}_\alpha - \mathbf R^{(1)}_\alpha|^2 \Big]^{1/2},

and the images are :math:`\mathbf R(q) = \mathbf R^{(1)} + q(\mathbf R^{(2)} - \mathbf R^{(1)})`, at
:math:`Q = q\,\Delta Q`. The total energy of each state is computed at the images by an electronic-structure method
(:ref:`validation`, Terminology). For TDDFT, :math:`E_\mathrm{state}(Q) = E_\mathrm{ref}(Q) + \Omega_\mathrm{root}(Q)`,
the ground-state total energy plus the excitation energy of the chosen root (``utils.parse_total_energy_qexml`` +
``utils.parse_west_summary``).

**Harmonic fit.** A parabola with its minimum fixed at the state's equilibrium geometry :math:`Q_0`,
:math:`V = a(Q - Q_0)^2 + c`, is fitted by least squares to a few images near the minimum (``fit_parabola``); it gives
:math:`\hbar\omega` (S1) and the single-mode Huang-Rhys factor

.. math::

   S = \frac{\Omega\,\Delta Q^2}{2\hbar} = \frac{\Delta Q^2}{2\ell^2},

the mean number of phonons emitted in a vertical transition (``config_coord_1d_vof_solver.huang_rhys``). The
Franck-Condon shift (peak of the emission sideband) is :math:`\approx S\hbar\omega_f` for nearly equal frequencies.


.. _si_ho:

S3. Harmonic oscillator eigenfunctions
======================================

.. math::

   \chi_n(Q) = \Big(\frac{\alpha}{\pi}\Big)^{1/4} \frac{1}{\sqrt{2^n n!}}\, H_n(\xi)\, e^{-\xi^2/2},
   \qquad \xi = \sqrt\alpha\,(Q - Q_c),
   \qquad E_n = \big(n + \tfrac12\big)\hbar\omega .

The explicit form overflows double precision for :math:`n \gtrsim 170` (:math:`2^n n!`). The normalized functions obey
the three-term recurrence, obtained from :math:`H_{n+1} = 2\xi H_n - 2nH_{n-1}` by dividing by
:math:`\sqrt{2^{n+1}(n+1)!}`,

.. math::

   \chi_0 = \Big(\frac{\alpha}{\pi}\Big)^{1/4} e^{-\xi^2/2}, \qquad
   \chi_1 = \sqrt2\,\xi\,\chi_0, \qquad
   \chi_{n+1} = \sqrt{\frac{2}{n+1}}\,\xi\,\chi_n - \sqrt{\frac{n}{n+1}}\,\chi_{n-1},

which never forms :math:`H_n` or :math:`2^n n!` and is stable for any :math:`n` (``harmonic_wavefunctions``;
validated against the explicit form to 1e-12 for :math:`n \le 20` and orthonormal to 1e-9 up to :math:`n = 300`). The
classical turning point of level :math:`n` is :math:`x_n = \sqrt{(2n+1)/\alpha}` (``classical_turning_point``), where
:math:`E_n = \tfrac12\Omega^2 x_n^2`.


.. _si_fc:

S4. Franck-Condon integrals: Ruhoff recurrences
===============================================

PyPL's original route (``config_coord_1d_solver.compute_franck_condon_integrals``; P. T. Ruhoff, Chem. Phys. 186, 355
(1994)). With :math:`\alpha_{i,f}` of the two wells (SI) and :math:`k = \Delta Q` (SI),

.. math::

   a = \frac{\alpha_i - \alpha_f}{\alpha_i + \alpha_f}, \quad
   b = \frac{2k\sqrt{\alpha_i}\,\alpha_f}{\alpha_i + \alpha_f}, \quad
   c = -a, \quad
   d = -\frac{2k\sqrt{\alpha_f}\,\alpha_i}{\alpha_i + \alpha_f}, \quad
   e = \frac{4\sqrt{\alpha_i\alpha_f}}{\alpha_i + \alpha_f},

.. math::

   F_{00} = \Big(\frac{e}{2}\Big)^{1/2} e^{bd/2e}
   = \Big(\frac{2\sqrt{\alpha_i\alpha_f}}{\alpha_i + \alpha_f}\Big)^{1/2}
     \exp\!\Big(-\frac{\alpha_i\alpha_f\,\Delta Q^2}{2(\alpha_i + \alpha_f)}\Big),

.. math::

   F_{0j} &= \tfrac{1}{\sqrt{2j}}\,d\,F_{0,j-1} + \sqrt{\tfrac{j-1}{j}}\,c\,F_{0,j-2}, \\
   F_{i0} &= \tfrac{1}{\sqrt{2i}}\,b\,F_{i-1,0} + \sqrt{\tfrac{i-1}{i}}\,a\,F_{i-2,0}, \\
   F_{ij} &= \tfrac{1}{\sqrt{2i}}\,b\,F_{i-1,j} + \sqrt{\tfrac{i-1}{i}}\,a\,F_{i-2,j}
            + \tfrac12\sqrt{\tfrac{j}{i}}\,e\,F_{i-1,j-1}.

For equal frequencies :math:`F_{00} = e^{-S/2}`, so :math:`|F_{00}|^2 = e^{-S}` is the Debye-Waller (zero-phonon)
weight.

**Exact properties used as checks.**

- *Closure.* :math:`\{\chi_j\}` is complete, so :math:`\sum_j |F_{ij}|^2 = 1` for every :math:`i` once enough final
  levels are kept; it fixes the basis size (e.g. :math:`1 - 10^{-6}` at :math:`j = 55` for MgO emission).
- *Identity.* At :math:`\Delta Q = 0` and :math:`\alpha_i = \alpha_f`: :math:`F_{ij} = \delta_{ij}`.
- *Symmetry.* At :math:`\alpha_i = \alpha_f`: :math:`a = c = 0` and :math:`b = -d`; exchanging the roles of the two
  wells reverses the displacement, :math:`\chi_n(-x) = (-1)^n\chi_n(x)`, hence :math:`F_{ij} = (-1)^{i+j}F_{ji}`, and
  PL and absorption lineshapes coincide.

**Numerical stability.** The recurrence accumulates rounding of :math:`\sim 10^{-8}` at the highest indices of a
50 :math:`\times` 60 basis for :math:`S \approx 3.5`. For :math:`S` of tens it becomes unstable in rows
:math:`i \gtrsim 10` (finding 6 of :ref:`validation`): row 0 stays exact, the higher rows grow without bound.


.. _si_hermite:

S5. Franck-Condon integrals by quadrature
=========================================

``compute_franck_condon_integrals_hermite`` evaluates

.. math::

   F_{ij} = \int \chi_i(Q + \Delta Q)\,\chi_j(Q)\,dQ \approx \sum_k w_k\,\chi_i(Q_k + \Delta Q)\,\chi_j(Q_k)

on a symmetric uniform grid :math:`|Q| \le 6 + \Delta Q + \max(x_{n_i-1}, x_{n_f-1})` (amu\ :sup:`1/2` Å) with an odd
number of points (default 100001) and composite Simpson weights :math:`w = \tfrac h3(1, 4, 2, 4, \dots, 2, 4, 1)`,
:math:`h = (Q_\mathrm{max} - Q_\mathrm{min})/(N - 1)`. (Computing :math:`h` as :math:`Q_1 - Q_0` loses about five digits
to cancellation.) The integrand is analytic and decays as a Gaussian, so the error is at machine precision; this route
is independent of S4 and stable at any :math:`S`.


.. _si_fd:

S6. Finite-difference eigenstates of any potential
==================================================

``solve_schrodinger_fd`` and ``compute_franck_condon_integrals_fd`` solve

.. math::

   \Big[-K\frac{d^2}{dx^2} + V(x)\Big]\psi_n = E_n\psi_n

on a uniform grid :math:`x_k`, with hard walls (:math:`\psi = 0` beyond the grid). Central differences,

.. math::

   \psi''_k &\approx \frac{\psi_{k-1} - 2\psi_k + \psi_{k+1}}{h^2} + O(h^2) \qquad \text{(3-point)}, \\
   \psi''_k &\approx \frac{-\psi_{k-2} + 16\psi_{k-1} - 30\psi_k + 16\psi_{k+1} - \psi_{k+2}}{12h^2} + O(h^4)
   \qquad \text{(5-point)},

follow from the Taylor expansion :math:`\psi_{k\pm m} = \sum_p (\pm mh)^p\psi^{(p)}_k/p!`: the 3-point combination
cancels all odd terms and leaves :math:`h^2\psi^{(4)}/12`; the 5-point weights also cancel the :math:`h^2` term,
leaving :math:`-h^4\psi^{(6)}/90`. The Hamiltonian is a symmetric banded matrix (bandwidth 1 or 2); its lowest
eigenpairs come from sparse shift-invert Lanczos (``eigsh`` with :math:`\sigma = \min V - 1` eV), so no dense matrix is
formed. Eigenvectors are normalized to :math:`\sum_k \psi_n(x_k)^2 h = 1` and signed positive on the right-hand tail,
the Hermite-function convention, so that :math:`F_{ij} = \sum_k \psi^{(i)}_i(x_k)\psi^{(f)}_j(x_k)\,h` has the same
signs as S4.

**Grid.** :math:`x \in [-\Delta Q - x^{(i)}_{n_i-1} - m,\; x^{(f)}_{n_f-1} + m]` with the harmonic turning points
:math:`x_n` (S3) and a margin :math:`m` (default 2 amu\ :sup:`1/2` Å), in the solver frame of S1;
:math:`Q = Q_f + s\,x`, :math:`s = \mathrm{sign}(Q_f - Q_i)`.

**Diagnostics** (``fd_info``): wall amplitude :math:`\max_n |\psi_n(x_\mathrm{edge})| / \max|\psi_n|` (bound states:
:math:`\sim 10^{-16}`); orthonormality :math:`\max|h\Psi\Psi^T - I|`; ladder deviation
:math:`\max_n |E_n - E_0 - n\hbar\omega|` (zero for a parabola). Convergence: the measured log-log slopes of the
:math:`F_{ij}` and ladder errors against :math:`h` are 2.00 (3-point) and 4.00 (5-point).


.. _si_fits:

S7. Fits of a sampled potential energy surface
==============================================

``config_coord_1d_potential(q, energies, q0)`` with ``fit_*`` methods. Every fit records the images it used, its
parameters and whether it is **bounded** (has bound states), which the finite-difference route requires.

**Parabola** :math:`a(Q - Q_0)^2 + c` (bounded if :math:`a > 0`).

**Polynomial** :math:`\sum_{k=0}^n p_k (Q - Q_0)^k`, least squares; bounded only for even :math:`n` with
:math:`p_n > 0`. A short sampled range leaves the highest coefficients poorly determined, and their sign can come out
negative (then :math:`V \to -\infty`).

**Morse** :math:`D(1 - e^{-b(Q - Q_m)})^2 + c`; :math:`V''(Q_m) = 2Db^2`, so the harmonic limit is
:math:`D \to \infty` at fixed :math:`Db^2`. Both signs of :math:`b` are tried.

**Cubic** :math:`a(Q - Q_0)^2 + b(Q - Q_0)^3 + c`, with :math:`a` optionally fixed to the parabola. It has a local
maximum (the *turnover*) where :math:`V' = 2a x + 3bx^2 = 0`, :math:`x_t = -2a/3b`, of height
:math:`V(x_t) - c = 4a^3/27b^2`; beyond it :math:`V \to -\infty`, so it is never bounded (S9).

**Splines** of degree :math:`k = 2 \dots 5` through the images (``fit_spline``). A spline of degree :math:`k` with
:math:`N` knots has :math:`N + k - 1` degrees of freedom; interpolating :math:`N` points leaves :math:`k - 1` conditions.
For odd :math:`k = 2m - 1` they are shared equally, :math:`m - 1` per end:

- *natural*: :math:`V^{(m)} = \dots = V^{(2m-2)} = 0` at both ends (cubic: :math:`V'' = 0`);
- *clamped*: :math:`V^{(1)} \dots V^{(m-1)}` equal to those of the tails (cubic: slope; quintic: slope and curvature);
- *not-a-knot*: no condition; the knots next to the ends are removed instead (any degree).

Even degrees allow only not-a-knot. With ``smoothing`` the spline is FITPACK's smoothing spline with the fewest knots
such that :math:`\sum_k (V(Q_k) - E_k)^2 \le N\,\mathrm{rms}^2`.

**Tails.** Beyond the first and last image :math:`Q_e` the potential must be continued:

- *harmonic* (default): :math:`a_\mathrm{tail}(Q - Q_0)^2 + \mathrm{const}`, continuous in value (C0); the slope jump
  :math:`2a_\mathrm{tail}(Q_e - Q_0) - V'(Q_e)` is recorded (zero only with *clamped* ends);
- *smooth*: the spline's own second-order Taylor expansion,
  :math:`V(Q_e) + V'(Q_e)(Q - Q_e) + \tfrac12 V''(Q_e)(Q - Q_e)^2`, continuous in value, slope and curvature (C2),
  leaving the spline inside the data unchanged; bounded iff :math:`V''(Q_e) > 0` at both ends.

The phonon energy at the minimum of any fit is :math:`\hbar\sqrt{V''(Q_\mathrm{min})}` (analytic for parabola, cubic
and splines; a central difference otherwise).


.. _si_lineshape:

S8. The Franck-Condon lineshape and spectra
===========================================

``bulid_fc_lsp`` (original) builds

.. math::

   I(E) = \sum_{i,j} p_i\,|F_{ij}|^2\,G(E; E_{ij}, \sigma_{ij}), \qquad
   E_{ij} = j\hbar\omega_f - i\hbar\omega_i, \qquad
   p_i = \frac{e^{-i\hbar\omega_i/k_BT}}{\sum_{i'} e^{-i'\hbar\omega_i/k_BT}},

with a normalized Gaussian whose width grows linearly with the line's distance from the zero-phonon line (ZPL),

.. math::

   G(E; \mu, \sigma) = \frac{1}{\sqrt{2\pi}\,\sigma}\,e^{-(E - \mu)^2/2\sigma^2}, \qquad
   \sigma = \sigma_0 + (\sigma_1 - \sigma_0)\,\frac{|\mu|}{160\ \mathrm{meV}},

where ``sigma = [σ0, σ1]`` (meV). With ``zpl_lorentzian`` the ZPL term :math:`|F_{00}|^2 G(E; 0, \sigma_0)` is replaced
by :math:`|F_{00}|^2 L(E; 0, \gamma)`, :math:`L = \frac{1}{\pi}\frac{\gamma}{(E - \mu)^2 + \gamma^2}`. Since
:math:`\sum_{ij} p_i|F_{ij}|^2 = 1` (closure), :math:`\int I\,dE = 1`; a Lorentzian narrower than the energy grid makes
the discrete sum slightly larger (1.0101 for :math:`\gamma = 0.3` meV on a 1 meV grid).

**Lines at actual levels.** ``build_fc_lsp_from_levels`` uses the same sum with
:math:`E_{ij} = (E^f_j - E^f_0) - (E^i_i - E^i_0)` and :math:`p_i \propto e^{-(E^i_i - E^i_0)/k_BT}` from the given
(e.g. finite-difference) levels; for harmonic ladders it reduces exactly to ``bulid_fc_lsp``.

**Spectra** (``compute_spectrum``), with :math:`E` measured from the ZPL:

.. math::

   I_\mathrm{PL}(\hbar\omega) \propto \mu^2\,(\hbar\omega)^3\,I(E_\mathrm{ZPL} - \hbar\omega), \qquad
   I_\mathrm{Abs}(\hbar\omega) \propto \mu^2\,\hbar\omega\,I(\hbar\omega - E_\mathrm{ZPL}).


.. _si_cubic:

S9. First-order cubic anharmonic correction
===========================================

Derived step by step and checked numerically in Tutorial 004 (supplementary S1-S4); summary:

.. math::

   \hat H = \hat H_0 + b\hat Q^3, \qquad \hat Q = \ell(\hat a + \hat a^\dagger), \qquad
   \hat H' = \hbar\Omega\,\xi\,(\hat a + \hat a^\dagger)^3, \qquad
   \xi = \frac{b\,\ell^3}{\hbar\Omega} = \frac{b}{\hbar\Omega}\Big(\frac{\hbar}{2\Omega}\Big)^{3/2}

(``anharmonic_strength``). The non-zero matrix elements of :math:`(\hat a + \hat a^\dagger)^3` from :math:`|n\rangle` are
:math:`3(n+1)^{3/2}` (:math:`n+1`), :math:`\sqrt{(n+1)(n+2)(n+3)}` (:math:`n+3`), :math:`3n^{3/2}` (:math:`n-1`) and
:math:`\sqrt{n(n-1)(n-2)}` (:math:`n-3`). First-order perturbation theory,
:math:`c_m = \langle m|\hat H'|n\rangle/(E_n - E_m)`, gives

.. math::

   \tilde\psi_n = \psi_n + \xi\Big[-3(n+1)^{3/2}\psi_{n+1} - \tfrac13\sqrt{(n+1)(n+2)(n+3)}\,\psi_{n+3}
   + 3n^{3/2}\psi_{n-1} + \tfrac13\sqrt{n(n-1)(n-2)}\,\psi_{n-3}\Big]

(``cubic_correction_coefficients``, ``perturbed_wavefunction``), with norm :math:`1 + \sum c_m^2`.
:math:`E^{(1)}_n = 0` (odd operator) and :math:`E^{(2)}_n = -\hbar\Omega\,\xi^2(30n^2 + 30n + 11)`.
``compute_franck_condon_integrals_cubic`` integrates :math:`F_{ij}` between the corrected states of both wells with the
quadrature of S5.

**Validity.** (1) Mixing :math:`\mu_n = 3|\xi|(n+1)^{3/2} \lesssim 0.1` for every level that carries weight,
:math:`n_\mathrm{valid} = (0.1/3|\xi|)^{2/3} - 1`; (2) :math:`|E^{(2)}_n| \ll \hbar\Omega`; (3) the level well inside the
turnover, :math:`|x_t| = \ell/6|\xi|`, and below the barrier :math:`\hbar\Omega/432\xi^2`. NV\ :sup:`-` in diamond
satisfies them (:math:`\xi \approx 10^{-3}`, Tutorial 004); NV\ :sup:`-` in MgO violates all three
(:math:`|\xi| \approx 10^{-2}`, :math:`S = 36`-70, turnovers inside the images; Tutorial 004a).


.. _si_isc:

S10. Intersystem crossing rate
==============================

Fermi's golden rule with a spin-orbit matrix element :math:`\lambda_\perp` and the VOF of the transition evaluated at
the electronic gap :math:`\Delta`,

.. math::

   \Gamma(\Delta) = \frac{4\pi}{\hbar}\,|\lambda_\perp|^2\,F(\Delta), \qquad \tau = 1/\Gamma

(``isc.compute_isc_rate``, ``compute_isc_lifetime``), with :math:`\lambda` converted to J (GHz: :math:`h\nu`; meV;
µeV; cm\ :sup:`-1`: :math:`hc\tilde\nu`) and :math:`F` from meV\ :sup:`-1` to J\ :sup:`-1`. The prefactor follows
Goldman et al., Phys. Rev. Lett. 114, 145502 (2015) and Jin et al., Phys. Rev. Lett. 135, 036401 (2025).


.. _si_multimode:

S11. The multi-mode (full phonon) path
======================================

For completeness, the formulas of ``hr_solver`` / ``hr_factors`` / ``lineshape``. With phonon frequencies
:math:`\omega_k` and normalized eigenvectors :math:`e_{k,\alpha i}` (atom :math:`\alpha`, direction :math:`i`),

.. math::

   \Delta Q_k = \sum_{\alpha i} \sqrt{M_\alpha}\,\Delta R_{\alpha i}\,e_{k,\alpha i}
   \quad \text{(displacements)}, \qquad
   \Delta Q_k = \frac{1}{\omega_k^2}\sum_{\alpha i} \frac{F_{\alpha i}}{\sqrt{M_\alpha}}\,e_{k,\alpha i}
   \quad \text{(forces)},

.. math::

   S_k = \frac{\omega_k\,\Delta Q_k^2}{2\hbar}, \qquad
   S(\hbar\omega) = \sum_k S_k\,G(\hbar\omega; \hbar\omega_k, \sigma_k),

with :math:`\sigma_k` interpolated linearly between ``sigma[0]`` at the lowest and ``sigma[1]`` at the highest
frequency. The lineshape follows from the generating function

.. math::

   G(t) = \exp\big[S(t) - S(0) + 2C(t) - 2C(0)\big], \quad
   S(t) = \sum_k S_k\,e^{-\sigma_k^2t^2/2\hbar^2}\,e^{i\omega_kt}, \quad
   C(t) = \sum_k n_k S_k\,e^{-\sigma_k^2t^2/2\hbar^2}\cos\omega_kt,

.. math::

   A(\hbar\omega) = \frac{1}{2\pi\hbar}\int dt\;e^{i\omega t - |t|\gamma_\mathrm{ZPL}/\hbar}\,G(t),

with Bose occupations :math:`n_k = (e^{\hbar\omega_k/k_BT} - 1)^{-1}`, by direct integration or FFT. Its
regression test is the tutorial-001 part of ``tests/test_legacy_tutorials.py``.


.. _si_metrics:

S12. Comparison metrics
=======================

For two lineshapes :math:`I_\mathrm{ref}`, :math:`I` on the same axis and two matrices :math:`F^\mathrm{ref}`, :math:`F`:

.. math::

   \delta_F = \max_{ij} |F_{ij} - F^\mathrm{ref}_{ij}|, \qquad
   \delta_\mathrm{VOF} = \frac{\max_E |I(E) - I_\mathrm{ref}(E)|}{\max_E I_\mathrm{ref}(E)}, \qquad
   \bar E = \frac{\sum_E E\,I(E)}{\sum_E I(E)} .

A wavefunction is compared by its infidelity :math:`1 - |\langle\psi|\phi\rangle|/(\|\psi\|\|\phi\|)`. An observed order
of accuracy is the slope of :math:`\log(\text{error})` against :math:`\log(h)` (finite differences) or
:math:`\log(b)` (perturbation theory: 2 for the error, 1 for the effect).
