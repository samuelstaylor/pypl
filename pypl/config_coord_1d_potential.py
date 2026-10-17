import numpy as np
from scipy import constants
from scipy.optimize import curve_fit, minimize_scalar

AMU = constants.physical_constants["atomic mass constant"][0]
SPLINE_ORDER_NAMES = {2: "quadratic", 3: "cubic", 4: "quartic", 5: "quintic"}


class config_coord_1d_potential:
    """
    Fitted potential energy curve of one electronic state along a 1D configuration coordinate.

    The curve is sampled by total-energy calculations at a few configurations (the images of a
    1D CCD). A ``fit_*`` method turns those samples into a smooth potential :math:`V(Q)` that can
    be evaluated anywhere, e.g. on the grid of a finite-difference Schrödinger solve.

    Parameters
    ----------
    q : ndarray of shape (N,)
        Configuration coordinate of each sample (in Å·amu^{1/2}).
    energies : ndarray of shape (N,)
        Total energy of each sample (in eV).
    q0 : float, optional
        Equilibrium geometry of this state (in Å·amu^{1/2}). Parabola fits keep the minimum
        fixed here, and every fit is expanded around it. Default is 0.0.

    Notes
    -----
    After a fit, the following attributes are set:

    - ``kind`` : str, the fit used (``'parabola'``, ``'polynomial'``, ``'morse'`` or ``'spline'``).
    - ``index`` : ndarray of int, the samples the fit was made to.
    - ``params`` : dict, the fitted parameters.
    - ``bounded`` : bool, whether :math:`V(Q) \\to +\\infty` (or a finite plateau) on both sides,
      i.e. whether the potential has bound states.
    """

    def __init__(self, q, energies, q0=0.0):

        self.q = np.asarray(q, dtype=float)
        self.energies = np.asarray(energies, dtype=float)
        self.q0 = q0

        self.kind = None
        self.index = None
        self.params = {}
        self.bounded = None
        self._potential = None
        self._curvature = None

    def _select(self, index):
        """Indices of the samples to fit; ``None`` means all of them."""

        if index is None:
            return np.arange(len(self.q))
        return np.arange(len(self.q))[index]

    def fit_parabola(self, index=None):
        r"""
        Fit a parabola with its minimum fixed at ``q0``.

        .. math::

            $V(Q) = a\,(Q - Q_0)^2 + c$

        Parameters
        ----------
        index : slice or array of int, optional
            Samples to fit. Default is all of them; a few points around the minimum give the
            local harmonic curvature.

        Returns
        -------
        dict
            ``{'a': eV/(amu Å^2), 'c': eV}``.
        """

        self.index = self._select(index)
        x, y = self.q[self.index] - self.q0, self.energies[self.index]
        (a, c), _ = curve_fit(lambda x, a, c: a * x**2 + c, x, y)

        self.kind = "parabola"
        self.params = {"a": a, "c": c}
        self.bounded = a > 0
        self._potential = lambda q: a * (np.asarray(q) - self.q0) ** 2 + c
        return self.params

    def fit_polynomial(self, degree, index=None):
        r"""
        Fit a polynomial expanded around ``q0``.

        .. math::

            $V(Q) = \sum_{k=0}^{n} p_k \, (Q - Q_0)^k$

        Parameters
        ----------
        degree : int
            Polynomial degree :math:`n`. Only an even degree with :math:`p_n > 0` is bounded below.
        index : slice or array of int, optional
            Samples to fit. Default is all of them.

        Returns
        -------
        dict
            ``{'degree': int, 'coeffs': ndarray}``, highest order first (as ``numpy.polyval``).
        """

        self.index = self._select(index)
        x, y = self.q[self.index] - self.q0, self.energies[self.index]
        coeffs = np.polyfit(x, y, degree)

        self.kind = "polynomial"
        self.params = {"degree": degree, "coeffs": coeffs}
        self.bounded = degree % 2 == 0 and coeffs[0] > 0
        self._potential = lambda q: np.polyval(coeffs, np.asarray(q) - self.q0)
        return self.params

    def fit_morse(self, index=None):
        r"""
        Fit a Morse potential.

        .. math::

            $V(Q) = D \left( 1 - e^{-b (Q - Q_m)} \right)^2 + c$

        A negative :math:`b` softens the well toward smaller :math:`Q`. Both signs are tried,
        starting from the local harmonic curvature, and the better fit is kept. Nearly harmonic
        samples drive :math:`D` toward the harmonic limit.

        Parameters
        ----------
        index : slice or array of int, optional
            Samples to fit. Default is all of them.

        Returns
        -------
        dict
            ``{'D': eV, 'b': (amu^{1/2} Å)^{-1}, 'q_m': Å·amu^{1/2}, 'c': eV}``.
        """

        def morse(q, d, b, q_m, c):
            return d * (1.0 - np.exp(-b * (q - q_m))) ** 2 + c

        self.index = self._select(index)
        x, y = self.q[self.index], self.energies[self.index]
        a_guess = np.polyfit(x - self.q0, y, 2)[0]
        best = None
        for sign in (1.0, -1.0):
            try:
                p = curve_fit(
                    morse, x, y, p0=[5.0, sign * np.sqrt(abs(a_guess) / 5.0), self.q0, y.min()], maxfev=20000
                )[0]
            except RuntimeError:
                continue
            rms = np.sqrt(np.mean((morse(x, *p) - y) ** 2))
            if p[0] > 0 and (best is None or rms < best[1]):
                best = (p, rms)
        if best is None:
            raise RuntimeError("Morse fit did not converge to a bound well.")
        d, b, q_m, c = best[0]

        self.kind = "morse"
        self.params = {"D": d, "b": b, "q_m": q_m, "c": c}
        self.bounded = True
        self._potential = lambda q: morse(np.asarray(q), d, b, q_m, c)
        return self.params

    def fit_cubic(self, index=None, curvature=None):
        r"""
        Fit a cubic expanded around the fixed minimum ``q0``, for the perturbative anharmonic correction.

        .. math::

            $V(Q) = a\,(Q - Q_0)^2 + b\,(Q - Q_0)^3 + c$

        There is no linear term, so :math:`Q_0` stays a stationary point. :math:`a` sets the harmonic
        frequency and :math:`b` the anharmonicity :math:`\hat{H}' = bQ^3`. Passing ``curvature`` fixes
        :math:`a` (e.g. to the parabola fit that defines the solver's frequency), so that the
        unperturbed Hamiltonian is exactly the harmonic one and only :math:`b` and :math:`c` are
        fitted; otherwise all three are fitted. A cubic is unbounded below on one side, so it has
        no true bound states: use it perturbatively, or in a finite box well inside its turnover
        point :math:`Q_0 - 2a / 3b`.

        Parameters
        ----------
        index : slice or array of int, optional
            Samples to fit. Default is all of them; ``b`` needs samples away from the minimum.
        curvature : float, optional
            Fixed :math:`a` in eV/(amu Å^2). Default is None: fit :math:`a` too.

        Returns
        -------
        dict
            ``{'a': eV/(amu Å^2), 'b': eV/(amu^{3/2} Å^3), 'c': eV, 'a_fixed': bool, 'turnover': Å·amu^{1/2}}``,
            where ``turnover`` is the position of the local maximum (infinite for ``b = 0``).
        """

        self.index = self._select(index)
        x, y = self.q[self.index] - self.q0, self.energies[self.index]
        if curvature is None:
            a, b, c = np.linalg.lstsq(np.column_stack([x**2, x**3, np.ones_like(x)]), y, rcond=None)[0]
        else:
            a = curvature
            b, c = np.linalg.lstsq(np.column_stack([x**3, np.ones_like(x)]), y - a * x**2, rcond=None)[0]
        turnover = self.q0 - 2 * a / (3 * b) if b != 0 else np.inf

        self.kind = "cubic"
        self.params = {"a": a, "b": b, "c": c, "a_fixed": curvature is not None, "turnover": turnover}
        self.bounded = b == 0 and a > 0
        self._potential = lambda q: a * (np.asarray(q) - self.q0) ** 2 + b * (np.asarray(q) - self.q0) ** 3 + c
        return self.params

    def fit_spline(self, index=None, tail_curvature=None, order=3, bc=None, smoothing=None, tails="harmonic"):
        r"""
        Spline through (or near) the samples, continued harmonically outside them.

        Inside the sampled range :math:`V` is a piecewise polynomial of degree ``order``. Outside it,

        .. math::

            $V(Q) = a_\text{tail}\,(Q - Q_0)^2 + \text{const}$

        with the constant chosen on each side so that :math:`V` is continuous. ``bc`` sets the end
        conditions of an interpolating spline:

        - ``'natural'``: derivatives :math:`m \dots 2m-2` vanish at both ends, for ``order`` :math:`= 2m-1`
          (cubic: :math:`V'' = 0`; quintic: :math:`V''' = V'''' = 0`). Odd ``order`` only.
        - ``'not-a-knot'``: no end conditions; the knots next to the ends are dropped instead. Any ``order``
          (an even ``order`` uses FITPACK's interpolating spline, which places its knots the same way).
        - ``'clamped'``: derivatives :math:`1 \dots m-1` at both ends equal those of the harmonic tails, so the
          spline joins its tails smoothly (cubic: matching slope, C1; quintic: slope and curvature, C2).
          Odd ``order`` only.

        With ``smoothing`` the spline does not pass through every sample: it is FITPACK's smoothing
        spline with the fewest knots whose rms residual is at most ``smoothing``. End conditions do not
        apply to it.

        ``tails`` sets the continuation outside the sampled range:

        - ``'harmonic'`` (default): the parabola above, :math:`a_\text{tail}` from the near-minimum fit;
          continuous (C0) at the end samples unless ``bc='clamped'`` bends the spline to match its slope.
        - ``'smooth'``: the spline's own second-order Taylor expansion at each end sample,

          .. math::

              $V(Q) = V(Q_e) + V'(Q_e)\,(Q - Q_e) + \tfrac{1}{2} V''(Q_e)\,(Q - Q_e)^2$

          so value, slope and curvature are continuous (C2) and nothing inside the sampled range changes.
          The tails depend only on the data; ``tail_curvature`` is not used. The fit is bound only if
          :math:`V''(Q_e) > 0` at both ends.

        Parameters
        ----------
        index : slice or array of int, optional
            Samples to use. Default is all of them.
        tail_curvature : float, optional
            :math:`a_\text{tail}` in eV/(amu Å^2). Default is a parabola fitted to the same samples.
        order : int, optional
            Polynomial degree of the pieces, 2 to 5 (3 = cubic, 4 = quartic, 5 = quintic). Default is 3.
        bc : {'natural', 'not-a-knot', 'clamped'}, optional
            End conditions of an interpolating spline. Default is ``'natural'`` for an odd ``order`` and
            ``'not-a-knot'`` for an even one.
        smoothing : float, optional
            Target rms residual (in eV) of a smoothing spline. Default is None: interpolate exactly.
        tails : {'harmonic', 'smooth'}, optional
            Continuation outside the samples: the harmonic parabola (default) or the spline's own C2
            Taylor expansion at each end.

        Returns
        -------
        dict
            ``{'order', 'bc', 'smoothing', 'tails', 'tail_curvature', 'shift_lo', 'shift_hi', 'range', 'n_knots',
            'slope_jump', 'curvature_jump', 'junction', 'end_slopes', 'end_curvatures'}``: the jumps are
            V(tail) - V(spline) at the two ends (eV/(amu^{1/2} Å) and eV/(amu Å^2)), ``junction`` is the
            continuity there (C0, C1 or C2), and ``end_slopes`` / ``end_curvatures`` are the spline's V' and V''
            at the two end samples. With ``tails='smooth'``, ``tail_curvature`` is the pair of tail curvatures
            V''/2 (low end, high end).
        """

        from scipy.interpolate import UnivariateSpline, make_interp_spline

        if order not in (2, 3, 4, 5):
            raise ValueError("spline order must be 2, 3, 4 or 5, got %r" % (order,))
        if smoothing is not None:
            if bc is not None:
                raise ValueError("end conditions do not apply to a smoothing spline; leave bc unset")
            if smoothing < 0:
                raise ValueError("smoothing is a target rms residual and must be >= 0, got %r" % (smoothing,))
        elif bc is None:
            bc = "natural" if order % 2 else "not-a-knot"
        if bc not in (None, "natural", "not-a-knot", "clamped"):
            raise ValueError("bc must be 'natural', 'not-a-knot' or 'clamped', got %r" % (bc,))
        if tails not in ("harmonic", "smooth"):
            raise ValueError("tails must be 'harmonic' or 'smooth', got %r" % (tails,))
        if tails == "smooth" and bc == "clamped":
            raise ValueError("'clamped' ends match the harmonic tails; with tails='smooth' the tails match the spline")
        if bc in ("natural", "clamped") and order % 2 == 0:
            raise ValueError(
                "'%s' end conditions need an odd order: an even-order spline needs an odd number (%d) of end "
                "conditions, which cannot be shared equally between the two ends; use 'not-a-knot'" % (bc, order - 1)
            )

        self.index = self._select(index)
        x, y = self.q[self.index], self.energies[self.index]
        idx = np.argsort(x)
        x, y = x[idx], y[idx]
        if tail_curvature is None:
            tail_curvature = np.polyfit(x - self.q0, y, 2)[0]
        lo, hi = x[0], x[-1]

        def tail_derivative(q, n):
            """n-th derivative of the harmonic tail a (Q - Q0)^2 + const."""

            return {1: 2 * tail_curvature * (q - self.q0), 2: 2 * tail_curvature}.get(n, 0.0)

        m = (order + 1) // 2
        if smoothing is not None:
            spline = UnivariateSpline(x, y, k=order, s=len(x) * smoothing**2)
        elif order % 2 == 0:
            spline = UnivariateSpline(x, y, k=order, s=0)
        elif bc == "natural":
            spline = make_interp_spline(x, y, k=order, bc_type=([(n, 0.0) for n in range(m, order)],) * 2)
        elif bc == "clamped":
            ends = tuple([(n, tail_derivative(end, n)) for n in range(1, m)] for end in (lo, hi))
            spline = make_interp_spline(x, y, k=order, bc_type=ends)
        else:
            spline = make_interp_spline(x, y, k=order, bc_type=None)
        d1, d2 = spline.derivative(1), spline.derivative(2)
        end_values = (float(spline(lo)), float(spline(hi)))
        end_slopes = (float(d1(lo)), float(d1(hi)))
        end_curvatures = (float(d2(lo)), float(d2(hi)))
        if tails == "smooth":
            # C2 Taylor tails: value, slope and curvature of the spline itself at each end sample.
            def potential(q):
                q = np.asarray(q, dtype=float)
                out = np.asarray(spline(np.clip(q, lo, hi)), dtype=float)
                for k, (end, mask) in enumerate(((lo, q < lo), (hi, q > hi))):
                    dx = q[mask] - end
                    out[mask] = end_values[k] + end_slopes[k] * dx + 0.5 * end_curvatures[k] * dx**2
                return out

            knots = spline.get_knots() if hasattr(spline, "get_knots") else np.unique(spline.t)
            self.kind = "spline"
            self.params = {
                "order": order,
                "bc": bc,
                "smoothing": smoothing,
                "tails": tails,
                "tail_curvature": (0.5 * end_curvatures[0], 0.5 * end_curvatures[1]),
                "shift_lo": end_values[0],
                "shift_hi": end_values[1],
                "range": (lo, hi),
                "n_knots": len(knots),
                "slope_jump": (0.0, 0.0),
                "curvature_jump": (0.0, 0.0),
                "junction": "C2",
                "end_slopes": end_slopes,
                "end_curvatures": end_curvatures,
            }
            self.bounded = min(end_curvatures) > 0
            self._potential = potential
            self._curvature = lambda q: float(d2(q))
            return self.params
        shift_lo = float(spline(lo)) - tail_curvature * (lo - self.q0) ** 2
        shift_hi = float(spline(hi)) - tail_curvature * (hi - self.q0) ** 2
        slope_jump = tuple(tail_derivative(end, 1) - float(d1(end)) for end in (lo, hi))
        curvature_jump = tuple(tail_derivative(end, 2) - float(d2(end)) for end in (lo, hi))
        scale = max(abs(tail_derivative(lo, 1)), abs(tail_derivative(hi, 1)), 2 * abs(tail_curvature))
        junction = "C0"
        if max(abs(j) for j in slope_jump) <= 1e-8 * scale:
            junction = "C2" if max(abs(j) for j in curvature_jump) <= 1e-8 * scale else "C1"

        def potential(q):
            q = np.asarray(q, dtype=float)
            tails = tail_curvature * (q - self.q0) ** 2 + np.where(q < lo, shift_lo, shift_hi)
            return np.where((q >= lo) & (q <= hi), spline(np.clip(q, lo, hi)), tails)

        knots = spline.get_knots() if hasattr(spline, "get_knots") else np.unique(spline.t)
        self.kind = "spline"
        self.params = {
            "order": order,
            "bc": bc,
            "smoothing": smoothing,
            "tails": tails,
            "tail_curvature": tail_curvature,
            "shift_lo": shift_lo,
            "shift_hi": shift_hi,
            "range": (lo, hi),
            "n_knots": len(knots),
            "slope_jump": slope_jump,
            "curvature_jump": curvature_jump,
            "junction": junction,
            "end_slopes": end_slopes,
            "end_curvatures": end_curvatures,
        }
        self.bounded = tail_curvature > 0
        self._potential = potential
        self._curvature = lambda q: float(d2(q))
        return self.params

    def __call__(self, q):
        """Evaluate the fitted potential (in eV) at ``q`` (in Å·amu^{1/2})."""

        if self._potential is None:
            raise RuntimeError("No fit yet: call one of the fit_* methods first.")
        return self._potential(q)

    def residuals(self, index=None):
        """Fit minus sample energy (in eV) at the samples ``index`` (default: all samples)."""

        idx = self._select(index)
        return self(self.q[idx]) - self.energies[idx]

    def rms(self, index=None):
        """Root-mean-square residual (in eV) over the samples ``index`` (default: all samples)."""

        return float(np.sqrt(np.mean(self.residuals(index) ** 2)))

    def minimum(self):
        """
        Minimum of the fitted potential within the sampled range.

        Returns
        -------
        q_min : float
            Position of the minimum (in Å·amu^{1/2}).
        e_min : float
            Energy at the minimum (in eV).
        """

        res = minimize_scalar(lambda q: float(self(q)), bounds=(self.q.min(), self.q.max()), method="bounded")
        return float(res.x), float(res.fun)

    def phonon_energy(self, dq=1e-3):
        r"""
        Harmonic phonon energy at the minimum of the fit.

        .. math::

            $\hbar\omega = \hbar \sqrt{V''(Q_\text{min})}$

        with :math:`V''` analytic for a parabola and a spline, and from a central difference of step
        ``dq`` otherwise. For a parabola this is :math:`\hbar\sqrt{2a}`, exactly the frequency the
        harmonic solver uses.

        Returns
        -------
        float
            :math:`\hbar\omega` in meV.
        """

        if self.kind in ("parabola", "cubic"):
            # A cubic around its fixed minimum has V''(Q0) = 2a: the frequency is the quadratic one.
            curvature = 2 * self.params["a"]
        elif self.kind == "spline":
            curvature = self._curvature(self.minimum()[0])
        else:
            q_min, _ = self.minimum()
            curvature = (self(q_min + dq) - 2 * self(q_min) + self(q_min - dq)) / dq**2
        # Same expression as tutorial 002, so a parabola reproduces its frequencies bit for bit.
        omega_sq = curvature * constants.eV / (1e-10**2 * AMU)
        return float(np.sqrt(omega_sq * constants.hbar**2 / constants.eV**2) * 1000)

    def describe(self):
        """One line naming the fit, its parameters and the samples it used."""

        if self.kind is None:
            return "not fitted"
        idx = self.index
        points = "samples %d-%d of %d (Q = %.3f..%.3f)" % (
            idx[0] + 1,
            idx[-1] + 1,
            len(self.q),
            self.q[idx].min(),
            self.q[idx].max(),
        )
        if self.kind == "parabola":
            form = "parabola a (Q - %.6f)^2 + c, a = %.6f eV/(amu A^2), c = %.6f eV" % (
                self.q0,
                self.params["a"],
                self.params["c"],
            )
        elif self.kind == "polynomial":
            form = "degree-%d polynomial in (Q - %.6f), coefficients (highest first) [%s]" % (
                self.params["degree"],
                self.q0,
                ", ".join("%.4g" % v for v in self.params["coeffs"]),
            )
        elif self.kind == "cubic":
            p = self.params
            form = "cubic a (Q - %.6f)^2 + b (Q - %.6f)^3 + c, a = %.6f eV/(amu A^2) (%s), b = %+.6f eV/(amu^1.5 A^3), " % (
                self.q0,
                self.q0,
                p["a"],
                "fixed" if p["a_fixed"] else "fitted",
                p["b"],
            ) + "c = %.6f eV, turnover at Q = %.2f" % (
                p["c"],
                p["turnover"],
            )
        elif self.kind == "morse":
            p = self.params
            form = "Morse D (1 - exp(-b (Q - Q_m)))^2 + c, D = %.6g eV, b = %.6g, Q_m = %.6f, c = %.6f eV" % (
                p["D"],
                p["b"],
                p["q_m"],
                p["c"],
            )
        else:
            p = self.params
            ends = "smoothing to rms <= %.3g eV" % p["smoothing"] if p["smoothing"] is not None else p["bc"] + " ends"
            if p.get("tails") == "smooth":
                return "%s; fitted to %s; rms residual %.2e eV%s" % (
                    "%s spline (degree %d, %s, %d knots), smooth C2 Taylor tails from the spline's own ends "
                    "(V' = %+.3f / %+.3f eV/(amu^0.5 A), V'' = %+.3f / %+.3f eV/(amu A^2) at Q = %.3f / %.3f)"
                    % (
                        SPLINE_ORDER_NAMES[p["order"]],
                        p["order"],
                        ends,
                        p["n_knots"],
                        p["end_slopes"][0],
                        p["end_slopes"][1],
                        p["end_curvatures"][0],
                        p["end_curvatures"][1],
                        p["range"][0],
                        p["range"][1],
                    ),
                    points,
                    self.rms(idx),
                    self._unbounded_reason(),
                )
            form = (
                "%s spline (degree %d, %s, %d knots), harmonic tails a = %.6f eV/(amu A^2), %s junction "
                "(slope jumps %.2e / %.2e eV/(amu^0.5 A))"
                % (
                    SPLINE_ORDER_NAMES[p["order"]],
                    p["order"],
                    ends,
                    p["n_knots"],
                    p["tail_curvature"],
                    p["junction"],
                    p["slope_jump"][0],
                    p["slope_jump"][1],
                )
            )
        return "%s; fitted to %s; rms residual %.2e eV%s" % (form, points, self.rms(idx), self._unbounded_reason())

    def _unbounded_reason(self):
        """Why the fit has no bound states, or an empty string if it has them."""

        if self.bounded:
            return ""
        if self.kind == "polynomial":
            degree, lead = self.params["degree"], self.params["coeffs"][0]
            if degree % 2:
                return "; UNBOUNDED: odd degree %d, so V -> -inf on one side" % degree
            return "; UNBOUNDED: leading coefficient p_%d = %.3g < 0, so V -> -inf as |Q| -> inf" % (degree, lead)
        if self.kind == "cubic":
            return (
                "; unbounded beyond the turnover at Q = %.2f (perturbative or in-box use only)"
                % self.params["turnover"]
            )
        if self.kind == "spline" and self.params.get("tails") == "smooth":
            return "; UNBOUNDED: V'' <= 0 at an end sample (%.3g / %.3g), so a smooth tail turns down" % tuple(
                self.params["end_curvatures"]
            )
        return "; UNBOUNDED: curvature <= 0"
