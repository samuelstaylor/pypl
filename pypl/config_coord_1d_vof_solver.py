import numpy as np
import scipy.sparse as sparse
from scipy import constants
from scipy.sparse.linalg import eigsh

from pypl.config_coord_1d_solver import config_coord_1d_solver

AMU = constants.physical_constants["atomic mass constant"][0]

# hbar^2 / 2 in eV amu Angstrom^2: kinetic prefactor for a mass-weighted coordinate in Å·amu^{1/2}.
KINETIC = constants.hbar**2 / (2 * AMU * 1e-20 * constants.eV)

# Central-difference stencils for d^2/dQ^2: coefficients for offsets 0, 1, 2, ..., in units of 1 / h^2.
FD_STENCILS = {
    3: (-2.0, 1.0),
    5: (-30.0 / 12, 16.0 / 12, -1.0 / 12),
    7: (-490.0 / 180, 270.0 / 180, -27.0 / 180, 2.0 / 180),
}
FD_STENCIL_NAMES = {
    3: "3-point central, O(h^2): (1, -2, 1) / h^2",
    5: "5-point central, O(h^4): (-1, 16, -30, 16, -1) / (12 h^2)",
    7: "7-point central, O(h^6): (2, -27, 270, -490, 270, -27, 2) / (180 h^2)",
}
FD_DEFAULT_NPTS = 20001
QUADRATURE_DEFAULT_NPTS = 100001


def grid_spacing(q):
    """Spacing of a uniform grid. ``q[1] - q[0]`` loses ~5 digits to cancellation far from Q = 0."""

    return (q[-1] - q[0]) / (len(q) - 1)


def uniform_grid(lo, hi, npts=None, h=None, odd=False):
    r"""
    Uniform grid on :math:`[lo, hi]`, given either the number of points or the spacing.

    Parameters
    ----------
    lo, hi : float
        Grid limits.
    npts : int, optional
        Number of points. Exactly one of ``npts`` and ``h`` must be given.
    h : float, optional
        Grid spacing. The grid starts at ``lo`` with exactly this spacing and is extended past ``hi`` to the next
        whole step, so it always covers :math:`[lo, hi]`.
    odd : bool, optional
        Round the number of points up to an odd number (composite Simpson rule). Default is False.

    Returns
    -------
    ndarray
        The grid.
    """

    if (npts is None) == (h is None):
        raise ValueError("give exactly one of npts and h")
    if h is not None:
        if h <= 0:
            raise ValueError("the grid spacing h must be positive, got %r" % (h,))
        npts = int(np.ceil((hi - lo) / h - 1e-9)) + 1
        if odd and npts % 2 == 0:
            npts += 1
        return lo + h * np.arange(npts)
    npts = int(npts)
    if npts < 3:
        raise ValueError("a grid needs at least 3 points, got %d" % npts)
    if odd and npts % 2 == 0:
        npts += 1
    return np.linspace(lo, hi, npts)


class config_coord_1d_vof_solver(config_coord_1d_solver):
    """
    One-dimensional configuration coordinate diagram solver with independent routes to the
    Franck-Condon integrals.

    Extends :class:`config_coord_1d_solver` without changing it. Every inherited method
    (``compute_franck_condon_integrals``, ``bulid_fc_lsp``, ``compute_spectrum``, ``Gaussian``,
    ``Lorentzian``, ``BZ``) behaves exactly as in the parent. The methods added here fill the same
    attributes (``fc_ints``, ``ni``, ``nf``), so ``bulid_fc_lsp`` and ``compute_spectrum`` work on
    their result unchanged.

    Routes to :math:`F_{ij} = \\langle \\chi_i(Q + \\Delta Q) | \\chi_j(Q) \\rangle`:

    - ``compute_franck_condon_integrals`` (inherited): Ruhoff recurrence relations.
    - ``compute_franck_condon_integrals_hermite``: numerical quadrature of explicit harmonic
      oscillator eigenfunctions.
    - ``compute_franck_condon_integrals_fd``: eigenvectors of a finite-difference Hamiltonian
      on arbitrary potentials, e.g. fits of the 1D CCD (``config_coord_1d_potential``).
    - ``compute_franck_condon_integrals_cubic``: harmonic eigenfunctions corrected to first order
      for a cubic anharmonicity :math:`bQ^3` of each potential (optional; same frequencies).

    Parameters
    ----------
    omega_i : float
        Initial vibrational frequency (in meV).
    omega_f : float
        Final vibrational frequency (in meV).
    delta_q : float
        Displacement along the configuration coordinate (in Å·amu^{1/2}).
    prec : int, optional
        Precision for `decimal` operations, default is 100.
    """

    @staticmethod
    def harmonic_alpha(omega):
        r"""
        Inverse squared length of a harmonic oscillator, :math:`\alpha = \omega / \hbar`.

        Parameters
        ----------
        omega : float
            Phonon energy :math:`\hbar\omega` (in meV).

        Returns
        -------
        float
            :math:`\alpha` in (amu Å^2)^{-1}.
        """

        return omega * 1e-3 * constants.eV / constants.hbar**2 * 1e-20 * AMU

    @staticmethod
    def huang_rhys(omega, delta_q):
        r"""
        Single-mode Huang-Rhys factor of a harmonic potential.

        .. math::

            $S = \frac{\omega \, \Delta Q^2}{2\hbar}$

        Same expression as tutorial 002, so the values match it bit for bit.

        Parameters
        ----------
        omega : float
            Phonon energy :math:`\hbar\omega` (in meV).
        delta_q : float
            Displacement along the configuration coordinate (in Å·amu^{1/2}).

        Returns
        -------
        float
            :math:`S`, dimensionless.
        """

        return delta_q**2 * 1e-10**2 * AMU * omega * 1e-3 * constants.eV / constants.hbar / (2 * constants.hbar)

    @staticmethod
    def harmonic_wavefunction(q, n, omega, center=0.0):
        r"""
        Normalized harmonic oscillator eigenfunction.

        .. math::

            $\chi_n(Q) = \left(\frac{\alpha}{\pi}\right)^{1/4} \frac{1}{\sqrt{2^n n!}}
            H_n\!\left(\sqrt{\alpha}\,(Q - Q_c)\right) e^{-\alpha (Q - Q_c)^2 / 2}$

        Parameters
        ----------
        q : ndarray
            Configuration coordinate (in Å·amu^{1/2}).
        n : int
            Vibrational quantum number.
        omega : float
            Phonon energy :math:`\hbar\omega` (in meV).
        center : float, optional
            Center :math:`Q_c` of the oscillator (in Å·amu^{1/2}). Default is 0.0.

        Returns
        -------
        ndarray
            :math:`\chi_n(Q)` in (Å·amu^{1/2})^{-1/2}.

        Notes
        -----
        Evaluated with ``harmonic_wavefunctions``; to get several levels, call that directly.
        """

        return config_coord_1d_vof_solver.harmonic_wavefunctions(q, n, omega, center)[n]

    @staticmethod
    def harmonic_wavefunctions(q, nmax, omega, center=0.0):
        r"""
        Normalized harmonic oscillator eigenfunctions :math:`\chi_0, \dots, \chi_{n_{max}}`.

        Built with the three-term recurrence of the normalized Hermite functions,

        .. math::

            $\chi_{n+1}(Q) = \sqrt{\frac{2}{n+1}}\,\xi\,\chi_n(Q) - \sqrt{\frac{n}{n+1}}\,\chi_{n-1}(Q),
            \qquad \xi = \sqrt{\alpha}\,(Q - Q_c)$

        starting from :math:`\chi_0 = (\alpha/\pi)^{1/4} e^{-\xi^2/2}`. It never forms
        :math:`H_n` or :math:`2^n n!` separately, so it stays accurate for any :math:`n`;
        the explicit product overflows float64 above :math:`n \approx 170`, which large
        Huang-Rhys factors (hundreds of levels) need.

        Parameters
        ----------
        q : ndarray
            Configuration coordinate (in Å·amu^{1/2}).
        nmax : int
            Highest vibrational quantum number.
        omega : float
            Phonon energy :math:`\hbar\omega` (in meV).
        center : float, optional
            Center :math:`Q_c` of the oscillator (in Å·amu^{1/2}). Default is 0.0.

        Returns
        -------
        ndarray
            Shape ``(nmax + 1,) + q.shape``; row :math:`n` is :math:`\chi_n(Q)` in
            (Å·amu^{1/2})^{-1/2}.
        """

        alpha = config_coord_1d_vof_solver.harmonic_alpha(omega)
        xi = np.sqrt(alpha) * (np.asarray(q, dtype=float) - center)
        chi = np.empty((nmax + 1,) + xi.shape)
        chi[0] = (alpha / np.pi) ** 0.25 * np.exp(-(xi**2) / 2)
        if nmax >= 1:
            chi[1] = np.sqrt(2.0) * xi * chi[0]
        for n in range(1, nmax):
            chi[n + 1] = np.sqrt(2.0 / (n + 1)) * xi * chi[n] - np.sqrt(n / (n + 1)) * chi[n - 1]
        return chi

    @staticmethod
    def classical_turning_point(n, omega):
        r"""
        Classical turning point of harmonic level :math:`n`, :math:`\sqrt{(2n + 1) / \alpha}`.

        Parameters
        ----------
        n : int
            Vibrational quantum number.
        omega : float
            Phonon energy :math:`\hbar\omega` (in meV).

        Returns
        -------
        float
            Distance from the minimum (in Å·amu^{1/2}).
        """

        return np.sqrt((2 * n + 1) / config_coord_1d_vof_solver.harmonic_alpha(omega))

    def compute_franck_condon_integrals_hermite(self, ni, nf, npts=QUADRATURE_DEFAULT_NPTS, h=None, margin=6.0):
        r"""
        Compute Franck-Condon integrals by numerical quadrature of harmonic eigenfunctions.

        .. math::

            $F_{ij} = \int \chi_i(Q + \Delta Q)\, \chi_j(Q)\, dQ$

        integrated with the composite Simpson rule on a uniform grid wide enough for the
        highest level of either oscillator. Independent of the Ruhoff recurrences, so it checks
        ``compute_franck_condon_integrals``.

        Parameters
        ----------
        ni : int
            Number of vibrational states in the initial potential.
        nf : int
            Number of vibrational states in the final potential.
        npts : int, optional
            Quadrature points; rounded up to an odd number. Default is 100001. Ignored if ``h`` is given.
        h : float, optional
            Quadrature grid spacing (in Å·amu^{1/2}) instead of ``npts``.
        margin : float, optional
            Grid half-width beyond :math:`\Delta Q` plus the outermost classical turning point (in Å·amu^{1/2}).
            Default is 6.0.

        Notes
        -----
        Sets ``self.fc_ints`` (shape ``(ni, nf)``), ``self.ni``, ``self.nf`` and ``self.hermite_grid``.
        """

        q, weights = self._simpson_grid(ni - 1, nf - 1, npts, h, margin)
        chi_i = self.harmonic_wavefunctions(q, ni - 1, self.omega_i, -self.delta_q)
        chi_f = self.harmonic_wavefunctions(q, nf - 1, self.omega_f, 0.0)

        self.fc_ints = (chi_i * weights) @ chi_f.T
        self.ni = ni
        self.nf = nf
        self.hermite_grid = q

    def _simpson_grid(self, top_i, top_f, npts=QUADRATURE_DEFAULT_NPTS, h=None, margin=6.0):
        """Symmetric quadrature grid wide enough for levels ``top_i`` / ``top_f``, and its Simpson weights."""

        tp = max(self.classical_turning_point(top_i, self.omega_i), self.classical_turning_point(top_f, self.omega_f))
        span = margin + self.delta_q + tp
        if h is None:
            q = np.linspace(-span, span, npts + 1 - npts % 2)
        else:
            q = uniform_grid(-span, span, h=h, odd=True)
            q -= 0.5 * (q[0] + q[-1])  # keep the grid symmetric
        weights = np.ones(len(q))
        weights[1:-1:2] = 4.0
        weights[2:-1:2] = 2.0
        weights *= grid_spacing(q) / 3.0
        return q, weights

    @staticmethod
    def anharmonic_strength(b, omega):
        r"""
        Dimensionless strength of a cubic anharmonicity :math:`\hat{H}' = b\,Q^3`.

        .. math::

            $\xi = \frac{b}{\hbar\Omega} \left( \frac{\hbar}{2\Omega} \right)^{3/2}$

        is :math:`b\,Q^3` measured in units of :math:`\hbar\Omega`, with :math:`Q` in units of the
        oscillator's zero-point length :math:`\sqrt{\hbar / 2\Omega}` (in the mass-weighted
        coordinate, mass 1). The first-order correction to level :math:`n` mixes in
        :math:`|n \pm 1\rangle` with weight :math:`3\xi\,n^{3/2}`, so the treatment is reliable only
        while :math:`3|\xi|(n + 1)^{3/2} \ll 1` for every level that matters, a stricter condition
        than :math:`|\xi| \ll 1` alone.

        Parameters
        ----------
        b : float
            Cubic coefficient of the potential (in eV/(amu^{3/2} Å^3)).
        omega : float
            Harmonic phonon energy :math:`\hbar\Omega` (in meV).

        Returns
        -------
        float
            :math:`\xi`, dimensionless.
        """

        # hbar / (2 Omega) in amu Å^2 is 1 / (2 alpha).
        return b / (omega * 1e-3) * (2 * config_coord_1d_vof_solver.harmonic_alpha(omega)) ** -1.5

    @staticmethod
    def cubic_correction_coefficients(n_levels, xi):
        r"""
        First-order expansion of the corrected states in the harmonic basis.

        .. math::

            $\tilde{\psi}_n = \psi_n + \xi \left[ -3(n+1)^{3/2}\psi_{n+1} - \frac{\sqrt{(n+1)(n+2)(n+3)}}{3}\psi_{n+3}
            + 3n^{3/2}\psi_{n-1} + \frac{\sqrt{n(n-1)(n-2)}}{3}\psi_{n-3} \right]$

        from :math:`|\tilde{n}\rangle = |n\rangle + \sum_{m \ne n} \frac{\langle m|bQ^3|n\rangle}{(n - m)\hbar\Omega}|m\rangle`.
        Terms with a negative index vanish (their prefactor is zero).

        Parameters
        ----------
        n_levels : int
            Number of corrected states :math:`\tilde{\psi}_0 \dots \tilde{\psi}_{n_\text{levels}-1}`.
        xi : float
            Anharmonic strength from ``anharmonic_strength``.

        Returns
        -------
        ndarray of shape (n_levels, n_levels + 3)
            Row :math:`n` holds the coefficients of :math:`\psi_0 \dots \psi_{n_\text{levels}+2}`, so the
            corrected states are ``C @ psi`` for harmonic functions ``psi`` of levels 0 .. n_levels + 2.
        """

        coeffs = np.zeros((n_levels, n_levels + 3))
        for n in range(n_levels):
            coeffs[n, n] = 1.0
            coeffs[n, n + 1] = -3.0 * xi * (n + 1) ** 1.5
            coeffs[n, n + 3] = -xi * np.sqrt((n + 1) * (n + 2) * (n + 3)) / 3.0
            if n >= 1:
                coeffs[n, n - 1] = 3.0 * xi * n**1.5
            if n >= 3:
                coeffs[n, n - 3] = xi * np.sqrt(n * (n - 1) * (n - 2)) / 3.0
        return coeffs

    @staticmethod
    def perturbed_wavefunction(q, n, omega, xi, center=0.0):
        r"""
        First-order corrected eigenfunction :math:`\tilde{\psi}_n(Q)` of :math:`V = a(Q-Q_c)^2 + b(Q-Q_c)^3`.

        A linear combination of ``harmonic_wavefunctions`` with the coefficients of
        ``cubic_correction_coefficients``. Its norm is :math:`1 + O(\xi^2)`; it is not renormalized.

        Parameters
        ----------
        q : ndarray
            Configuration coordinate (in Å·amu^{1/2}).
        n : int
            Vibrational quantum number.
        omega : float
            Harmonic phonon energy :math:`\hbar\Omega` (in meV), from the quadratic coefficient.
        xi : float
            Anharmonic strength from ``anharmonic_strength``.
        center : float, optional
            Minimum :math:`Q_c` of the potential (in Å·amu^{1/2}). Default is 0.0.

        Returns
        -------
        ndarray
            :math:`\tilde{\psi}_n(Q)` in (Å·amu^{1/2})^{-1/2}.
        """

        row = config_coord_1d_vof_solver.cubic_correction_coefficients(n + 1, xi)[n]
        chi = config_coord_1d_vof_solver.harmonic_wavefunctions(q, len(row) - 1, omega, center)
        return sum(c * chi[m] for m, c in enumerate(row) if c != 0.0)

    def compute_franck_condon_integrals_cubic(
        self,
        ni,
        nf,
        b_i,
        b_f,
        q_i=None,
        q_f=0.0,
        npts=QUADRATURE_DEFAULT_NPTS,
        normalize=False,
        warn_threshold=0.1,
        h=None,
        margin=6.0,
    ):
        r"""
        Compute Franck-Condon integrals with a first-order cubic anharmonic correction.

        Each potential is :math:`V(Q) = a(Q - Q_0)^2 + b(Q - Q_0)^3 + c`. The harmonic frequencies
        ``omega_i`` / ``omega_f`` (from :math:`a`) and the Huang-Rhys factors are unchanged; the cubic
        term :math:`\hat{H}' = bQ^3` only corrects the vibrational states to first order
        (``perturbed_wavefunction``), and

        .. math::

            $F_{ij} = \int \tilde{\psi}^{(i)}_i(Q + \Delta Q)\, \tilde{\psi}^{(f)}_j(Q)\, dQ$

        is integrated with the same Simpson quadrature as ``compute_franck_condon_integrals_hermite``.
        The level energies stay harmonic: :math:`\langle n|Q^3|n\rangle = 0`, so the first-order
        energy shift vanishes, and ``bulid_fc_lsp`` places the lines at the harmonic positions.

        With ``b_i = b_f = 0`` this is exactly ``compute_franck_condon_integrals_hermite``.

        Validity: first-order perturbation theory needs the mixing coefficients to be small. The
        largest, :math:`3|\xi|(n + 1)^{3/2}`, grows with the level, so a small :math:`\xi` is not enough
        when high levels are kept. A warning is printed if :math:`|\xi| >` ``warn_threshold``, or if the
        coefficients exceed it within the basis; ``compute_franck_condon_integrals_fd`` solves the same
        potentials without approximation (inside a finite box, since a cubic is unbounded below).

        Parameters
        ----------
        ni : int
            Number of vibrational states in the initial potential.
        nf : int
            Number of vibrational states in the final potential.
        b_i : float
            Cubic coefficient of the initial-state potential (in eV/(amu^{3/2} Å^3)), on the same
            :math:`Q` axis as ``q_i`` / ``q_f``.
        b_f : float
            Cubic coefficient of the final-state potential.
        q_i : float, optional
            Minimum of the initial potential. Default is ``q_f - delta_q``.
        q_f : float, optional
            Minimum of the final potential. Default is 0.0. Only the orientation matters: the
            solver's frame has the initial minimum at :math:`-\Delta Q`, so when ``q_f < q_i`` the
            coordinate is mirrored and both cubic coefficients change sign.
        npts : int, optional
            Quadrature points; rounded up to an odd number. Default is 100001.
        normalize : bool, optional
            Renormalize each corrected state to 1. Default is False (pure first order, norm
            :math:`1 + O(\xi^2)`).
        warn_threshold : float, optional
            Limit on :math:`|\xi|` and on the mixing coefficients before a warning. Default is 0.1.
        h : float, optional
            Quadrature grid spacing (in Å·amu^{1/2}) instead of ``npts``.
        margin : float, optional
            Quadrature grid half-width beyond :math:`\Delta Q` and the outermost turning point. Default is 6.0.

        Notes
        -----
        Sets ``self.fc_ints`` (shape ``(ni, nf)``), ``self.ni``, ``self.nf`` and

        - ``self.cubic_xi_i``, ``self.cubic_xi_f`` : :math:`\xi` of each state in the solver frame.
        - ``self.cubic_norms_i``, ``self.cubic_norms_f`` : :math:`\langle\tilde{\psi}_n|\tilde{\psi}_n\rangle`
          before any renormalization.
        - ``self.cubic_grid``, ``self.cubic_wavefunctions_i``, ``self.cubic_wavefunctions_f`` : the
          corrected states on the quadrature grid.
        - ``self.cubic_info`` : dict of ``b`` in the solver frame, :math:`\xi`, the largest mixing
          coefficient, the highest level whose coefficients all stay below ``warn_threshold``
          (``n_valid``), and the warnings issued.
        """

        if q_i is None:
            q_i = q_f - self.delta_q
        sign = 1.0 if q_f >= q_i else -1.0
        b_i, b_f = sign * b_i, sign * b_f
        xi_i = self.anharmonic_strength(b_i, self.omega_i)
        xi_f = self.anharmonic_strength(b_f, self.omega_f)
        coeffs_i = self.cubic_correction_coefficients(ni, xi_i)
        coeffs_f = self.cubic_correction_coefficients(nf, xi_f)

        q, weights = self._simpson_grid(ni + 2, nf + 2, npts, h, margin)
        chi_i = self.harmonic_wavefunctions(q, ni + 2, self.omega_i, -self.delta_q)
        chi_f = self.harmonic_wavefunctions(q, nf + 2, self.omega_f, 0.0)
        psi_i, psi_f = coeffs_i @ chi_i, coeffs_f @ chi_f
        norms_i, norms_f = (psi_i**2 * weights).sum(axis=1), (psi_f**2 * weights).sum(axis=1)
        if normalize:
            psi_i /= np.sqrt(norms_i)[:, None]
            psi_f /= np.sqrt(norms_f)[:, None]

        warnings = []
        info = {"b_i": b_i, "b_f": b_f, "xi_i": xi_i, "xi_f": xi_f, "frame_sign": sign}
        for role, xi, n in (("initial", xi_i, ni), ("final", xi_f, nf)):
            mixing = 3 * abs(xi) * np.arange(1, n + 1) ** 1.5
            n_valid = int(np.sum(mixing <= warn_threshold)) - 1
            info["max_mixing_%s" % role[0]] = float(mixing.max())
            info["n_valid_%s" % role[0]] = n_valid
            if abs(xi) > warn_threshold:
                warnings.append(
                    "Warning: |xi| = %.3g > %g for the %s state: first-order perturbation theory is unreliable; "
                    "use compute_franck_condon_integrals_fd (the numerical Schrodinger solver) instead."
                    % (abs(xi), warn_threshold, role)
                )
            elif mixing.max() > warn_threshold:
                warnings.append(
                    "Warning: %s state: the first-order mixing 3|xi|(n+1)^1.5 reaches %.3g at n = %d (|xi| = %.2e); it "
                    "stays below %g only up to n = %d, so the correction is perturbative only for the lower levels; "
                    "compute_franck_condon_integrals_fd solves the full potential."
                    % (role, mixing.max(), n - 1, abs(xi), warn_threshold, n_valid)
                )
        for w in warnings:
            print(w)
        info["warnings"] = warnings

        self.fc_ints = (psi_i * weights) @ psi_f.T
        self.ni = ni
        self.nf = nf
        self.cubic_xi_i, self.cubic_xi_f = xi_i, xi_f
        self.cubic_norms_i, self.cubic_norms_f = norms_i, norms_f
        self.cubic_grid = q
        self.cubic_wavefunctions_i, self.cubic_wavefunctions_f = psi_i, psi_f
        self.cubic_info = info

    @staticmethod
    def solve_schrodinger_fd(potential, q, n, stencil=5, shift=None):
        r"""
        Lowest eigenstates of a 1D Hamiltonian by finite differences.

        .. math::

            $H = -\frac{\hbar^2}{2} \frac{d^2}{dQ^2} + V(Q)$

        in the mass-weighted coordinate :math:`Q`, discretized on the uniform grid ``q`` with a
        central-difference stencil and hard walls (:math:`\psi = 0` beyond the grid). The banded
        matrix is diagonalized by sparse shift-invert Lanczos below the potential minimum, so no
        dense ``len(q) x len(q)`` matrix is formed.

        Parameters
        ----------
        potential : ndarray
            :math:`V` on ``q`` (in eV).
        q : ndarray
            Uniform grid (in Å·amu^{1/2}).
        n : int
            Number of eigenstates.
        stencil : {3, 5, 7}, optional
            Points in the second-derivative stencil, O(h^2), O(h^4) or O(h^6). Default is 5.
        shift : float, optional
            Shift-invert target of the eigensolver (in eV). Default is :math:`\min V - 1` eV, below every state.

        Returns
        -------
        energies : ndarray of shape (n,)
            Eigenvalues (in eV), ascending.
        psi : ndarray of shape (n, len(q))
            Eigenfunctions normalized to :math:`\sum_k \psi(q_k)^2 h = 1` and signed like the
            Hermite functions: positive on the right-hand tail.
        """

        h = grid_spacing(q)
        coeffs = FD_STENCILS[stencil]
        bands = [np.full(len(q) - k, -KINETIC * c / h**2) for k, c in enumerate(coeffs)]
        bands[0] = bands[0] + potential
        offsets = list(range(1, len(coeffs)))
        hamiltonian = sparse.diags(bands[:0:-1] + bands, [-k for k in offsets[::-1]] + [0] + offsets, format="csc")

        if stencil not in FD_STENCILS:
            raise ValueError("stencil must be one of %s, got %r" % (sorted(FD_STENCILS), stencil))
        if n >= len(q) - 1:
            raise ValueError("%d eigenstates need more than %d grid points" % (n, len(q)))
        if shift is None:
            shift = float(np.min(potential)) - 1.0
        energies, vectors = eigsh(hamiltonian, k=n, sigma=shift, which="LM")
        order = np.argsort(energies)
        energies, psi = energies[order], vectors[:, order].T / np.sqrt(h)
        for v in range(n):
            # The last point above 1e-3 of the peak lies in the outer tail, past the last node.
            tail = np.nonzero(np.abs(psi[v]) > 1e-3 * np.max(np.abs(psi[v])))[0][-1]
            psi[v] *= np.sign(psi[v, tail])
        return energies, psi

    def fd_grid(self, ni, nf, q_i=None, q_f=0.0, npts=None, h=None, margin=2.0, q_range=None):
        r"""
        The finite-difference grid that ``compute_franck_condon_integrals_fd`` would use, without solving.

        By default the grid spans, in the solver frame (final minimum at 0, initial at :math:`-\Delta Q`), from the
        initial well's outermost classical turning point of level ``ni - 1`` to the final well's of level ``nf - 1``,
        plus ``margin`` on each side, with ``npts`` points (default 20001). ``h`` sets the spacing instead of the
        number of points, and ``q_range`` sets the limits on the :math:`Q` axis of the potentials instead of the
        turning points and the margin.

        Parameters
        ----------
        ni, nf : int
            Vibrational states of the initial and final potential.
        q_i, q_f : float, optional
            Minima of the initial and final potential. Default: ``q_f - delta_q`` and 0.0.
        npts : int, optional
            Grid points. Default is 20001 (if ``h`` is not given).
        h : float, optional
            Grid spacing (in Å·amu^{1/2}) instead of ``npts``.
        margin : float, optional
            Padding beyond the outermost harmonic turning points (in Å·amu^{1/2}). Default is 2.0.
        q_range : tuple of float, optional
            ``(q_min, q_max)``: explicit grid limits on the :math:`Q` axis of the potentials.

        Returns
        -------
        x : ndarray
            Grid in the solver frame.
        q : ndarray
            The same points on the :math:`Q` axis of the potentials.
        info : dict
            ``sign`` (+1 or -1, the frame), ``delta`` (:math:`|Q_f - Q_i|`), ``turning_points`` and ``grid_from``
            (``'turning points + margin'`` or ``'q_range'``).
        """

        if q_i is None:
            q_i = q_f - self.delta_q
        if npts is not None and h is not None:
            raise ValueError("give npts or h, not both")
        if npts is None and h is None:
            npts = FD_DEFAULT_NPTS
        sign = 1.0 if q_f >= q_i else -1.0
        delta = abs(q_f - q_i)
        turn_i = self.classical_turning_point(ni - 1, self.omega_i)
        turn_f = self.classical_turning_point(nf - 1, self.omega_f)
        if q_range is None:
            lo, hi = -delta - turn_i - margin, turn_f + margin
            grid_from = "turning points + margin"
        else:
            ends = sorted(sign * (np.asarray(q_range, dtype=float) - q_f))
            lo, hi = ends
            grid_from = "q_range"
            if not (lo < -delta and hi > 0.0):
                raise ValueError("q_range %s must contain both minima (%g and %g)" % (tuple(q_range), q_i, q_f))
        x = uniform_grid(lo, hi, npts=npts, h=h)
        return (
            x,
            q_f + sign * x,
            {"sign": sign, "delta": delta, "turning_points": (turn_i, turn_f), "grid_from": grid_from},
        )

    def compute_franck_condon_integrals_fd(
        self,
        ni,
        nf,
        potential_i,
        potential_f,
        q_i=None,
        q_f=0.0,
        npts=None,
        stencil=5,
        margin=2.0,
        h=None,
        q_range=None,
        shift=None,
    ):
        r"""
        Compute Franck-Condon integrals from finite-difference eigenstates of two potentials.

        Each potential is solved with ``solve_schrodinger_fd`` on one shared grid and

        .. math::

            $F_{ij} = \sum_k \psi^{(i)}_i(x_k)\, \psi^{(f)}_j(x_k)\, h$

        No harmonic approximation enters, so any potential can be used (e.g. an anharmonic fit
        of the 1D CCD). For harmonic potentials of frequencies ``omega_i``, ``omega_f`` this
        reproduces ``compute_franck_condon_integrals``.

        The grid is laid out in the solver's frame, in which the final minimum is at 0 and the
        initial one at :math:`-\Delta Q`: :math:`x = s\,(Q - Q_f)` with :math:`s = \text{sign}(Q_f - Q_i)`.
        So emission on a CCD with the excited state at larger :math:`Q` is solved in a mirrored
        frame, and the result follows the same sign convention as the parent class.

        Parameters
        ----------
        ni : int
            Number of vibrational states in the initial potential.
        nf : int
            Number of vibrational states in the final potential.
        potential_i : callable
            Initial-state potential :math:`V_i(Q)` (eV), :math:`Q` in Å·amu^{1/2}.
        potential_f : callable
            Final-state potential :math:`V_f(Q)` (eV).
        q_i : float, optional
            Minimum of the initial potential on the :math:`Q` axis of ``potential_i``.
            Default is ``q_f - delta_q``.
        q_f : float, optional
            Minimum of the final potential. Default is 0.0.
        npts : int, optional
            Grid points. Default is 20001 (if ``h`` is not given).
        stencil : {3, 5, 7}, optional
            Second-derivative stencil, O(h^2), O(h^4) or O(h^6). Default is 5.
        margin : float, optional
            Grid padding beyond the outermost harmonic turning point (in Å·amu^{1/2}). Default is 2.0.
        h : float, optional
            Grid spacing (in Å·amu^{1/2}) instead of ``npts``.
        q_range : tuple of float, optional
            ``(q_min, q_max)``: explicit grid limits on the :math:`Q` axis of the potentials, instead of the turning
            points and ``margin``. See ``fd_grid``.
        shift : float, optional
            Shift-invert target of the eigensolver (eV). Default is :math:`\min V - 1` eV.

        Notes
        -----
        Sets ``self.fc_ints`` (shape ``(ni, nf)``), ``self.ni``, ``self.nf`` and

        - ``self.fd_grid`` : grid in the solver frame :math:`x`.
        - ``self.fd_q`` : the same points on the :math:`Q` axis of the potentials.
        - ``self.fd_energies_i``, ``self.fd_energies_f`` : eigenvalues (in eV).
        - ``self.fd_wavefunctions_i``, ``self.fd_wavefunctions_f`` : eigenfunctions on ``fd_grid``.
        - ``self.fd_info`` : dict of grid and accuracy diagnostics, including ``points_per_wavelength``, the
          number of grid points per shortest local de Broglie wavelength :math:`2\pi\sqrt{K/(E_\mathrm{max} -
          V_\mathrm{min})}` of the highest level (a few tens or more for a converged result).
        """

        x, q, grid = self.fd_grid(ni, nf, q_i, q_f, npts=npts, h=h, margin=margin, q_range=q_range)
        sign, (turn_i, turn_f) = grid["sign"], grid["turning_points"]
        npts = len(x)
        h = grid_spacing(x)

        v_i, v_f = potential_i(q), potential_f(q)
        e_i, psi_i = self.solve_schrodinger_fd(v_i, x, ni, stencil, shift)
        e_f, psi_f = self.solve_schrodinger_fd(v_f, x, nf, stencil, shift)

        def per_wavelength(energies, potential):
            kinetic = max(float(energies[-1] - np.min(potential)), 1e-300)
            return float(2 * np.pi * np.sqrt(KINETIC / kinetic) / h)

        def wall(psi):
            return float(np.max(np.abs(psi[:, [0, -1]]) / np.max(np.abs(psi), axis=1)[:, None]))

        def ortho(psi):
            return float(np.max(np.abs((psi * h) @ psi.T - np.eye(len(psi)))))

        def ladder(energies, omega):
            return float(np.max(np.abs((energies - energies[0]) * 1e3 - omega * np.arange(len(energies)))))

        self.fc_ints = (psi_i * h) @ psi_f.T
        self.ni = ni
        self.nf = nf
        self.fd_grid = x
        self.fd_q = q
        self.fd_energies_i, self.fd_energies_f = e_i, e_f
        self.fd_wavefunctions_i, self.fd_wavefunctions_f = psi_i, psi_f
        self.fd_info = {
            "npts": npts,
            "h": h,
            "stencil": stencil,
            "margin": margin,
            "grid_from": grid["grid_from"],
            "points_per_wavelength": min(per_wavelength(e_i, v_i), per_wavelength(e_f, v_f)),
            "frame": "x = Q - Q_f" if sign > 0 else "x = -(Q - Q_f)",
            "x_range": (x[0], x[-1]),
            "q_range": (min(q[0], q[-1]), max(q[0], q[-1])),
            "turning_points": (turn_i, turn_f),
            "wall": max(wall(psi_i), wall(psi_f)),
            "orthonormality": max(ortho(psi_i), ortho(psi_f)),
            "ladder_i": ladder(e_i, self.omega_i),
            "ladder_f": ladder(e_f, self.omega_f),
        }

    def build_fc_lsp_from_levels(
        self, eneaxis=None, temp=4, sigma=10, zpl_lorentzian=True, gamma=1, energies_i=None, energies_f=None
    ):
        r"""
        Build the Franck-Condon lineshape with lines at given vibrational level energies.

        .. math::

            $I(E) = \sum_{i, j} p_i \, |F_{ij}|^2 \, G(E; E_{ij}, \sigma), \qquad
            E_{ij} = (E^f_j - E^f_0) - (E^i_i - E^i_0)$

        The same sum as ``bulid_fc_lsp`` (same Gaussian with linearly growing width, same
        Lorentzian ZPL replacement), except that line :math:`ij` sits at the difference of the
        actual level energies instead of :math:`j\hbar\omega_f - i\hbar\omega_i`, and the
        Boltzmann weights :math:`p_i` use the actual initial levels. With anharmonic wells
        (``compute_franck_condon_integrals_fd`` on a fitted potential) this keeps the line
        positions consistent with the potential that produced :math:`F_{ij}`; with harmonic
        ladders it equals ``bulid_fc_lsp``.

        Parameters
        ----------
        eneaxis : ndarray, optional
            Energy axis in meV, relative to the 0-0 line. Default is `np.linspace(-1000, 1000, 2001)`.
        temp : float, optional
            Temperature in Kelvin. Default is 4 K.
        sigma : float or list of float, optional
            Gaussian broadening in meV, as in ``bulid_fc_lsp``. Default is 10 meV.
        zpl_lorentzian : bool, optional
            Replace the zero-phonon line Gaussian with a Lorentzian. Default is True.
        gamma : float, optional
            Lorentzian width of the ZPL in meV. Default is 1 meV.
        energies_i, energies_f : ndarray, optional
            Level energies of the initial and final state in meV (only differences matter), at least
            ``ni`` / ``nf`` of them. Default: ``fd_energies_i`` / ``fd_energies_f`` from the last
            ``compute_franck_condon_integrals_fd`` (stored in eV).

        Returns
        -------
        eneaxis : ndarray
            Energy axis (meV).
        fc_lineshape : ndarray
            Lineshape intensity (1/meV).

        Notes
        -----
        Sets ``self.energy_v`` (meV), ``self.prefactors``, ``self.eneaxis``, ``self.fc_lineshape`` and the
        broadening attributes, as ``bulid_fc_lsp`` does.
        """

        if eneaxis is None:
            eneaxis = np.linspace(-1000, 1000, 2001)
        if energies_i is None or energies_f is None:
            if not hasattr(self, "fd_energies_i"):
                raise ValueError("give energies_i / energies_f, or run compute_franck_condon_integrals_fd first")
            energies_i = self.fd_energies_i * 1e3 if energies_i is None else energies_i
            energies_f = self.fd_energies_f * 1e3 if energies_f is None else energies_f
        e_i = np.asarray(energies_i, dtype=float)[: self.ni]
        e_f = np.asarray(energies_f, dtype=float)[: self.nf]
        if len(e_i) < self.ni or len(e_f) < self.nf:
            raise ValueError("need at least ni = %d initial and nf = %d final level energies" % (self.ni, self.nf))
        self.temp = temp
        self.sigma = sigma
        self.zpl_lorentzian = zpl_lorentzian
        self.gamma = gamma
        self.energy_v = -(e_i - e_i[0])[:, None] + (e_f - e_f[0])[None, :]
        prefactors = self.BZ(e_i - e_i[0], temp)
        prefactors /= np.sum(prefactors)
        self.prefactors = prefactors
        lineshape = self.fc_ints[:, :, None] ** 2 * self.Gaussian(
            eneaxis[None, None, :], self.energy_v[:, :, None], sigma
        )
        lineshape = np.sum(prefactors[:, None, None] * lineshape, axis=(0, 1))
        if zpl_lorentzian:
            lineshape -= self.fc_ints[0, 0] ** 2 * self.Gaussian(eneaxis, 0.0, sigma)
            lineshape += self.fc_ints[0, 0] ** 2 * self.Lorentzian(eneaxis, 0.0, gamma)
        print("Integral check:", np.sum(lineshape) * (eneaxis[1] - eneaxis[0]))
        self.eneaxis = eneaxis
        self.fc_lineshape = lineshape
        return self.eneaxis, self.fc_lineshape
