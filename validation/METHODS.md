# VOF methods: a tutorial overview

This file explains every way the test harness can compute the **vibrational overlap function (VOF)** of the
1D configuration-coordinate (1D CCD) model: what each method does, the physics behind it, how to run it, how
to read its output, how accurate it is, and when to use it. Numbers are for the NV⁻ center in diamond
(tutorial 002) at the default settings, unless noted.

> **Keep this file current.** Whenever a method or fit is added or changed, update the matching section here
> (and the comparison table at the end).

Contents:

1. [The common framework](#1-the-common-framework)
2. [`pypl`: PyPL's Ruhoff recurrences (the reference)](#2-pypl-pypls-ruhoff-recurrences-the-reference)
3. [`hermite`: quadrature of harmonic wavefunctions](#3-hermite-quadrature-of-harmonic-wavefunctions)
4. [`fd`: finite-difference Schrödinger solver](#4-fd-finite-difference-schrödinger-solver)
5. [Potential fits for `fd`](#5-potential-fits-for-fd)
6. [`cubic`: first-order cubic anharmonic correction](#6-cubic-first-order-cubic-anharmonic-correction)
7. [Choosing a method](#7-choosing-a-method)

---

## 1. The common framework

### The quantity being computed

In the 1D CCD model each electronic state is a 1D potential well along the mass-weighted coordinate Q
(amu^½ Å). The ground state (GS) has its minimum at Q = 0, and the excited state (ES) has it at Q = ΔQ = 0.6957.
A transition from vibrational level i of the initial state to level j of the final state has intensity
proportional to the squared **Franck–Condon overlap**

F_ij = ⟨χ_i(Q + ΔQ) | χ_j(Q)⟩,

and the VOF is the broadened sum

I(E) = Σ_ij p_i |F_ij|² G(E; E_ij, σ),

where p_i is the Boltzmann population of initial level i, E_ij the transition energy, and G a Gaussian
(with a Lorentzian for the zero-phonon line). Emission is ES → GS; absorption is GS → ES.

### Datasets

Every script runs on one dataset, chosen with `--dataset` (default `nv_diamond`), defined in
`harness/ccd_data.py`:

| dataset | transition | 'ES' (upper) | 'GS' (lower) | ΔQ (amu^½ Å) | ħω (meV) / S, ES and GS | basis |
| --- | --- | --- | --- | --- | --- | --- |
| `nv_diamond` | 3E → 3A2, tutorial 002 | 3E, at Q = ΔQ | 3A2, at Q = 0 | 0.6957 | 64.31 / 3.72, 59.09 / 3.42 | 50 ES × 60 GS |
| `nv_mgo` | 3E → 1A1, WEST TDDFT-DDH (`doc/nv_mgo_west`) | 3E, at Q = 0 (red) | 1A1, at Q = ΔQ (blue) | 3.0979 | 60.83 / 69.8, 31.55 / 36.2 | 160 ES × 100 GS |

The roles are always called GS (lower) and ES (upper) in the code, and "emission" is ES → GS. Each dataset
places the two minima (`ccd_data.center`). MgO follows `doc/nv_mgo_west/reference/ccd_nv_mgo.py`: Q = q ΔQ from
`Image_<q>`, 3E at Q = 0, 1A1 at Q = ΔQ, energies from the 3A2 minimum, 3E red, 1A1 blue. So "emission" is
3E → 1A1, the transition of interest (higher → lower, ZPL ≈ 0.556 eV). The solver methods take both minima,
so F_ij is the same in either orientation (`check_1d_ccd.py` checks fd = hermite with the 3E well on the left).

**Images per fit** (q of `Image_<q>`, inclusive; every legend, note and `index.md` prints them):

| fit | 3E | 1A1 |
| --- | --- | --- |
| parabola (the harmonic wells) | q = −0.10 … +0.45 (12) | q = +0.60 … +1.15 (12) |
| cubic a Q² + b Q³ + c (the perturbative correction) | q = −0.10 … +1.15 (26) | q = −0.10 … +1.15 (26) |
| every other fit (splines, Morse, polyN) | all 33 (q = −1.00 … +2.00) | all 30 (q = −1.00 … +1.50) |

The one image never used is `sf/Image_+2.00`, which was not a spin-flip run, so its root 4 is not 1A1. With
these parabolas the harmonic wells are 3E 60.83 meV (S = 69.8) and 1A1 31.55 meV (S = 36.2). `validation/run_all_nv_mgo.toml` holds the settings the MgO
Huang-Rhys factors need: basis, a −0.5 … 10 eV axis, constant σ, and `hermite` as the reference.

### What a "method" is

**Every method computes only F_ij.** The VOF itself is always built by PyPL's own, unmodified
`config_coord_1d_solver.bulid_fc_lsp`, which also places every line at the harmonic energies
E_ij = j ħω_f − i ħω_i. So when two methods give different VOFs, the difference is entirely in F_ij.

The methods run on `config_coord_1d_vof_solver` (`pypl/config_coord_1d_vof_solver.py`). It is a subclass of PyPL's solver that
only *adds* methods, so all of PyPL's own functions behave exactly as before; `check_1d_ccd.py` verifies this
bit for bit.

### Basis size and directions

Each well keeps a finite number of levels: 50 ES and 60 GS. Emission uses (ni, nf) = (50, 60); absorption
swaps them to (60, 50). The thermal weights at 5 K put essentially all the population in initial level 0.

### Running and comparing methods

The simplest way is the master suite, which runs every method and fit with one parameter set and writes a
method-comparison table to `validation/outputs/suite/<name>/index.md`. Edit `validation/run_all.toml` (or pass
`--set section.key=value`) to choose the methods, the fd grid and stencil, the fits and spline options, and the
cubic variants:

```bash
python validation/run_all.py
```

Or run single scripts:

```bash
python validation/compare_vof.py                       # quick PASS/FAIL: pypl vs hermite and fd
python validation/compare_vof.py --test cubic          # any method(s) against the reference
python validation/plot_ccd_wavefunctions.py --wf-method fd --vof-method pypl hermite fd
python validation/vof_report.py                        # everything, in seven sections (~8 min)
python validation/check_1d_ccd.py                      # 59 regression checks (~5 s)
```

Each run writes `validation/outputs/<script>/<run>/`. Start with its `index.md`; every method's full settings
are in `log.txt`.

---

## 2. `pypl`: PyPL's Ruhoff recurrences (the reference)

**What it is.** PyPL's original method, `compute_franck_condon_integrals`. For two displaced harmonic
oscillators, F_ij satisfies three-term recurrence relations (Ruhoff, Chem. Phys. 186, 355 (1994)). Starting
from F_00, every other entry is built by recursion. No wavefunction is ever formed and no integral is taken.

**Inputs.** Only the two harmonic frequencies and ΔQ:

| | value | from |
| --- | --- | --- |
| ħω_GS | 59.089 meV | parabola fitted to GS images 1–5 of 15 (Q = −0.139 … 0.139) |
| ħω_ES | 64.312 meV | parabola fitted to ES images 5–11 of 11 (Q = 0.417 … 0.835) |
| ΔQ | 0.6957 amu^½ Å | GS and ES DFT geometries |

**Accuracy.** Exact for harmonic wells, apart from float64 rounding: the forward recursion accumulates an
error of up to 1.7e-8 at the highest indices, (49, 48) of a 50 × 60 basis. Below i, j < 30 the error is about
5e-11. It is invisible in any spectrum.

**Large Huang-Rhys factors: unstable.** The recursion is not stable when S is large. On NV⁻ in MgO
(`--dataset nv_mgo`: S = 50 and 85, 160 × 100 levels) row 0 is still exact (1e-15 against `hermite`), but
the error grows by about a decade per initial level: above 1e-8 from i ≈ 10–12, above 1e-3 from i ≈ 17–21,
above 1 from i ≈ 22–27, and |F_ij| reaches 1e31 although an overlap can never exceed 1. At 5 K only i = 0 is
populated, so the 5 K VOF is still right (1e-14), but any finite-temperature or closure use of the higher
rows is wrong. On that dataset `hermite` is the reference and `pypl` FAILs the like-for-like step on F_ij.

**When to use.** As the reference for harmonic wells with moderate S (diamond: S ≈ 3.5). It is also by far
the fastest (0.01 s). For S of tens or more, use `hermite` and check `pypl` against it.

---

## 3. `hermite`: quadrature of harmonic wavefunctions

**What it is.** It builds the harmonic eigenfunctions explicitly,

χ_n(Q) = (α/π)^{1/4} (2ⁿ n!)^{-1/2} Hₙ(√α (Q − Q_c)) e^{−α(Q−Q_c)²/2},  α = Ω/ħ,

evaluated (`harmonic_wavefunctions`) with the stable three-term recurrence of the normalized functions,
χ_{n+1} = √(2/(n+1)) ξ χ_n − √(n/(n+1)) χ_{n−1}, ξ = √α (Q − Q_c), from χ_0 = (α/π)^{1/4} e^{−ξ²/2}. It never
forms Hₙ or 2ⁿ n! separately. The explicit product overflows float64 above n ≈ 170, which the MgO dataset
(160 levels) reached. The recurrence equals the explicit formula to 1e-12 for n ≤ 20 and stays orthonormal to
1e-11 up to n = 400 (both pinned in `check_1d_ccd.py`). It then integrates F_ij = ∫ χ_i(Q + ΔQ) χ_j(Q) dQ with the composite Simpson rule on a 100001-point grid. It
answers the same question as `pypl` by a completely independent route, so it checks the recurrences.

**Run it.** `python validation/compare_vof.py --test hermite`, or
`python validation/plot_ccd_wavefunctions.py --wf-method hermite` to see the wavefunctions on the CCD.

**Accuracy.**
- VOF vs `pypl`: 7e-16 of the peak, i.e. machine precision.
- F_ij vs `pypl`: 1.7e-8, which is `pypl`'s own rounding, not Hermite's.

The grid spacing must be computed as (q[−1] − q[0]) / (N − 1). The earlier `q[1] − q[0]` lost about five
digits to cancellation and gave a hidden 1e-11 normalization error.

On NV⁻ in MgO (S = 50 / 85, 160 × 100 levels), `hermite` and `fd` agree to 1.3e-8 in F_ij and 4e-9 in
the VOF, which is what makes `hermite` the reference there.

**When to use.** To validate `pypl`, as the harmonic reference when S is large, or whenever you need the
harmonic wavefunctions themselves.

---

## 4. `fd`: finite-difference Schrödinger solver

**What it is.** It solves the 1D Schrödinger equation numerically on any potential V(Q):

H = −(ħ²/2) d²/dQ² + V(Q).

The steps are:

1. **Sample the potential.** V is one of the fits in section 5, sampled on a uniform grid (default 20001
   points).
2. **Discretize the second derivative.** A 5-point central stencil, (−1, 16, −30, 16, −1)/(12h²), accurate to
   O(h⁴). The 3-point stencil is O(h²). The walls are hard: ψ = 0 beyond the grid.
3. **Find the lowest states.** The lowest ni and nf eigenpairs come from a sparse shift-invert solver
   (`scipy.sparse.linalg.eigsh`). The dense alternative, `eig_banded`, hangs at 5001 points.
4. **Overlap them.** Both wells share one grid, so F_ij is a single product of eigenvector matrices.

No harmonic formula enters, so `fd` works for anharmonic potentials too.

**Grid.** It spans Q = −4.90 … 5.23 amu^½ Å: the outermost harmonic turning points (v = 59 GS, v = 49 ES) plus
2.0 amu^½ Å of padding. The spacing is h = 5.07e-4, which gives 302–318 points per shortest wavelength.
Emission is solved in a mirrored coordinate so that it uses the same sign convention as `pypl`.

**Run it.**

```bash
python validation/compare_vof.py --test fd                     # parabola wells: must match pypl
python validation/compare_vof.py --test fd --fd-fit spline     # anharmonic wells
python validation/plot_ccd_wavefunctions.py --wf-method fd --fd-fit spline --vof-method pypl fd
```

Flags: `--fd-npts`, `--fd-stencil {3,5}`, `--fd-margin`, `--fd-fit`.

**Accuracy (parabola wells).**
- vs `hermite`: 3.5e-9 in F_ij.
- vs `pypl`: 1.3e-12 of the VOF peak.
- Level ladder: within 2e-6 meV of v ħω for all 60 levels.
- At the walls: |ψ| ≈ 4e-16, the float64 floor.
- Convergence (report section 02): the measured log–log slopes are 2.00 for the 3-point stencil and 4.00 for
  the 5-point one, exactly the theoretical orders.

**When to use.** For any potential that isn't a parabola, or to check the harmonic methods. It is the
"exact" reference for anharmonic wells.

---

## 5. Potential fits for `fd`

The DFT images only sample Q = −0.14 … 0.83, but the FD grid, and the high vibrational levels, extend to
about ±5. **Every fit is therefore extrapolated over most of the grid**, and how a fit behaves outside the
data matters as much as how well it passes through it. Report section 01 shows each fit both near the
data and over the full grid.

| `--fd-fit` | potential | fitted to | rms GS / ES (meV) | VOF vs `pypl` (peak) | verdict |
| --- | --- | --- | --- | --- | --- |
| `parabola` (default) | a(Q−Q₀)² + c | GS 1–5, ES 5–11 | 2.11 / 0.84 | 1.3e-12 | the harmonic wells; use to test |
| `poly2` | quadratic | all images | 0.15 / 0.08 | 2.4e-3 | a different harmonic well |
| `poly4` | quartic | all images | 0.05 / 0.04 | — | **unbounded** (p₄ < 0) |
| `poly6` | sextic | all images | 0.024 / 0.006 | 7.0e-2 | bounded, but its walls are thousands of eV off; unusable |
| `poly8` | octic | all images | 0.010 / 0.004 | — | **unbounded** (p₈ < 0) |
| `spline` (+ variants) | piecewise polynomial | all images | 0 (exact) | 2.1–2.2e-2 | the only physical anharmonic fit |
| `morse` | D(1 − e^{−b(Q−Q_m)})² + c | all images | 0.11 / 0.04 | 1.9e-2 | collapses to the harmonic limit (D ≈ 2000 eV) |
| `cubic` | a Q² + b Q³ + c | all images, a fixed | 0.16 / 0.16 | 3.5e-2 | exact counterpart of method `cubic`; bound only in the box |

**Why polynomials of degree 4 and 8 are unbounded.** An even-degree polynomial has bound states only if its
leading coefficient is positive. Least squares does not enforce that, and here the leading coefficient comes
out negative (GS p₄ = −0.018, p₈ = −0.42), so V → −∞ at large |Q|. Over the short sampled range the top-order
term contributes only a few meV, so its sign is set by small softening and noise.

**Spline options.** `--fd-fit spline` with:

| flag | values |
| --- | --- |
| `--spline-order` | 2, 3, 4, 5 |
| `--spline-bc` | `natural`, `not-a-knot`, `clamped` |
| `--spline-smoothing RMS_meV` | smoothing target |
| `--spline-tails` | `harmonic` (default), `smooth` |

You can also name the variant directly in `--fd-fit`, e.g. `spline5-clamped`. The end conditions:
- **natural:** V'' = 0 at the end images. Odd orders only.
- **not-a-knot:** no conditions imposed. Any order.
- **clamped:** slope (and for order 5, curvature) matched to the harmonic tails, a smooth join. Odd orders only.

All six variants tested give the same physics to within a few percent of each other. They differ only in
ħω at the minimum, where the natural cubic is the outlier: it forces V'' = 0 at the first GS image, just
0.14 amu^½ Å from the minimum.

**Line positions of an anharmonic fit.** `bulid_fc_lsp` places line ij at j ħω_f − i ħω_i, the harmonic
ladder, whatever potential produced F_ij. `config_coord_1d_vof_solver.build_fc_lsp_from_levels(...)` builds the
same sum (same Gaussian, Lorentzian ZPL and Boltzmann weights) with each line at the actual level energies,
E_ij = (E^f_j − E^f_0) − (E^i_i − E^i_0), by default the `fd_energies_*` of the last fd solve. With harmonic
ladders it equals `bulid_fc_lsp` exactly (checked at 300 K); `bulid_fc_lsp` itself is unchanged. Use it whenever the
VOF's energy axis matters, e.g. an ISC rate F(Δ): on NV⁻ in MgO the smooth-spline ¹A₁ ladder is 373 meV at v = 10
against 315.5 meV harmonic. `doc/tutorials/pypl_003.ipynb` uses it.

### 5.1 Spline tails: `harmonic` or `smooth`

**The problem.** A spline is only defined between the first and last image. The finite-difference grid is
wider, because the highest vibrational levels reach past the data, so every spline needs a *tail* on each
side. What the tail is, and how it joins the spline, is a modelling choice with no data behind it.

**The two options.**

| tails | beyond the end image Q_e | junction at Q_e | depends on |
| --- | --- | --- | --- |
| `harmonic` (default) | a (Q − Q₀)² + const, a from the near-minimum parabola | C0: value only, a slope kink (C1 / C2 only with `-clamped`) | the parabola fit |
| `smooth` (`-smooth`) | V(Q_e) + V′(Q_e)(Q − Q_e) + ½ V″(Q_e)(Q − Q_e)² | C2: value, slope and curvature continuous | the spline's own ends, i.e. the data |

The `smooth` tail is the second-order Taylor expansion of the spline at its end image. It is the unique
quadratic that continues the curve with no jump in V, V′ or V″, so the potential is C² everywhere: inside the
data (a cubic spline is C² by construction) and across both junctions. That is the property an electronic-structure
PES has, smooth and twice differentiable. And because the tail adapts to the spline, not the reverse, **the spline
inside the data is untouched**: identical to the same spline with harmonic tails, to 0.0 eV.

**Why not `clamped`?** `clamped` also makes the junction smooth (C¹ for a cubic, C² for a quintic), but it does
so the other way round: it forces the spline's end slope to equal the harmonic tail's, 2a (Q_e − Q₀). That
slope comes from the parabola near the minimum, far from the end image, and it is usually wrong there:

| end image (NV⁻ in MgO) | slope of the data | slope `clamped` forces |
| --- | --- | --- |
| ³E, Q = 6.20 | +1.65 eV/(amu^½ Å) | +5.49 |
| ³E, Q = −3.10 | −4.50 | −2.74 |
| ¹A₁, Q = 4.65 | +1.16 | +0.37 |

To meet the steeper slope at Q = 6.20, the ³E spline swings down across its widest data gap (Q = 4.65–6.20)
and opens a well 0.26 eV (cubic) to 0.60 eV (quintic) deep that no image supports. A clamped spline is smooth,
but smooth *and wrong* inside the data. So `-smooth` cannot be combined with `-clamped`; the two are opposite
answers to the same question.

**Validity condition.** A smooth tail is a parabola with curvature ½ V″(Q_e). It rises, so the well stays bound,
only if V″ > 0 at both end images. The fit reports this: `bounded` is False otherwise, `describe()` gives the
end curvatures with an `UNBOUNDED` reason, and `fd` / `run_all.py` refuse the fit. Measured end curvatures, all
positive:

| dataset | state | V″ at the two end images (eV/(amu Å²)) |
| --- | --- | --- |
| NV⁻ in MgO | ³E | +2.024 (Q = −3.10), +0.805 (Q = 6.20) |
| NV⁻ in MgO | ¹A₁ | +0.179 (Q = −3.10), +0.875 (Q = 4.65) |
| NV⁻ in diamond | GS | +0.859, +0.898 |
| NV⁻ in diamond | ES | +0.968, +1.071 |

**Run it.** Any spline can take smooth tails, except `-clamped`:

```bash
python validation/plot_ccd_wavefunctions.py --dataset nv_mgo --wf-method fd --vof-method hermite fd \
    --fd-fit spline3-notaknot-smooth                              # the spec suffix -smooth
python validation/compare_vof.py --dataset nv_mgo --reference hermite --test fd \
    --fd-fit spline --spline-bc not-a-knot --spline-tails smooth  # the same fit from flags
```

In a suite config, add `tails = "smooth"` to a `[[fd.splines]]` entry (as `run_all_nv_mgo.toml` does).
In code: `config_coord_1d_potential.fit_spline(..., tails="smooth")`. Canonical names add `-smooth`:
`spline3-notaknot-smooth`, `spline3-natural-smooth`, `spline3-s0.05-smooth`, ...

**Reading the output.** The fit's description line names the tails and states the junction values, e.g.
`cubic spline (degree 3, not-a-knot ends, 28 knots), smooth C2 Taylor tails from the spline's own ends
(V' = -3.397 / +1.160, V'' = +0.179 / +0.875 at Q = -3.098 / 4.647)`. `params["junction"]` is `C2`, and the
slope and curvature jumps are zero. Legends call it `cubic spline (not-a-knot, smooth C2 tails)`.

**Measured effect.**
- *NV⁻ in MgO, 3E → 1A1:* `spline3-notaknot-smooth` changes the VOF by 0.2% of its peak relative to
  `spline3-notaknot` (centroid 580.3 meV from the ZPL, −566 meV from harmonic). Its ¹A₁ ladder differs only
  from v ≈ 30 (1027.5 vs 1004.5 meV), above the end image at 0.84 eV, where 5 K emission puts no weight.
  `clamped` changes the same VOF by 38–51%.
- *NV⁻ in diamond:* the tails matter more, because the images span only Q = −0.14 … 0.83 and the 50–60 levels
  reach far into them. Smooth tails change the VOF by 2.0% of the peak relative to harmonic tails (centroid
  −2.4 meV in emission, −4.0 meV in absorption). With them, the not-a-knot spline is within 0.19% of the
  harmonic reference, versus ~2% with harmonic tails. **Most of the ~2% "spline vs harmonic" difference
  reported for diamond (README, "Spline options") is therefore the slope kink of the harmonic tail, not anharmonicity of the
  data.**

**When to use it.** Whenever the PES should be smooth everywhere, which is always for an electronic-structure
surface, and in particular whenever levels reach beyond the data. The recommended spline is
**`spline3-notaknot-smooth`**: not-a-knot imposes nothing at the end images (the end shape comes from the data),
and the smooth tails add nothing artificial at the junctions. Report `spline4-notaknot-smooth` or
`spline3-s0.05-smooth` as robustness checks. Keep `harmonic` tails only to reproduce older results.

**Wavefunctions.** FD on the parabola reproduces the Hermite functions to an infidelity of about 1e-15. FD
on the spline changes them by 8e-4 (GS) and 8e-5 (ES); that change is the real anharmonicity of the DFT data.

---

## 6. `cubic`: first-order cubic anharmonic correction

**What it is.** An optional, analytic anharmonic correction. Each well is fitted with a cubic,

V(Q) = a(Q − Q₀)² + b(Q − Q₀)³ + c,

and the cubic term H′ = bQ³ is treated as a small perturbation of the harmonic Hamiltonian H₀. The harmonic
frequencies (from a), the HR factors and the level energies are unchanged; only the wavefunctions are
corrected. With b = 0 it is exactly the `hermite` method.

### 6.1 Fitting b: "a fixed" vs "a fitted"

- **a fixed** (default): a is set to the parabola's curvature, the one that defines ħω for every method.
  Only b and c are fitted, to all DFT images. Then H₀ is exactly the harmonic Hamiltonian, and nothing in the
  potential is dropped.
- **a fitted** (`--cubic-free-a`): a, b and c are all fitted together. a then differs slightly (GS 0.4208
  vs 0.4176), and the leftover (a_fit − a_parabola)Q² falls into neither H₀ nor the perturbation.
- `--cubic-points parabola` fits only the near-minimum images.

b is poorly determined by these nearly harmonic data. For the GS it is 0.0092 (a fixed, all images), 0.0053
(a fitted) or 0.042 (near-minimum images), and the ES b even changes sign.

### 6.2 Derivation of the correction

First-order perturbation theory writes the corrected state as the harmonic state plus admixtures of the
others:

|ñ⟩ = |n⟩ + Σ_{m≠n} cₘ |m⟩,  cₘ = ⟨m|bQ³|n⟩ / ((n − m)ħΩ).

In the mass-weighted coordinate, Q = √(ħ/2Ω)(a + a†). Expanding (a + a†)³, only orderings with a net change of
±1 or ±3 connect |n⟩ and |m⟩. For m = n + 1, three orderings contribute:

| ordering | amplitude |
| --- | --- |
| a†a†a\|n⟩ | n√(n+1) |
| a†aa†\|n⟩ | (n+1)√(n+1) |
| aa†a†\|n⟩ | (n+2)√(n+1) |
| **sum** | **3(n+1)^{3/2}** |

With the dimensionless strength ξ = (b/ħΩ)(ħ/2Ω)^{3/2}, the corrected wavefunction is

ψ̃ₙ = ψₙ + ξ[ −3(n+1)^{3/2} ψₙ₊₁ − √((n+1)(n+2)(n+3))/3 · ψₙ₊₃ + 3n^{3/2} ψₙ₋₁ + √(n(n−1)(n−2))/3 · ψₙ₋₃ ].

Terms with negative indices vanish. F_ij is then integrated between the corrected wavefunctions of *both*
states, ∫ ψ̃ᵢ(Q + ΔQ) ψ̃ⱼ(Q) dQ, with the Hermite quadrature. Building ψ̃ₙ needs harmonic functions up to n + 3.

Here ξ_GS = +1.03e-3 and ξ_ES = −1.15e-3. In emission, which is solved in the mirrored coordinate, both b
change sign.

### 6.3 Mixing: why the correction weakens at high n

The **mixing** is the coefficient of a neighbouring state in ψ̃ₙ. The largest is 3|ξ|(n+1)^{3/2}, the weight
of ψₙ₊₁. It grows with n:
- the matrix element ⟨n+1|Q³|n⟩ grows as n^{3/2}, because higher levels extend further out, where bQ³ is larger;
- the energy gap it is divided by stays ħΩ.

So the *same* b perturbs every higher level more. ξ alone only measures the effect on n = 0.

First-order perturbation theory is the first term of a series in these coefficients. The terms it drops are
of size about (mixing)². With mixing below 0.1 the dropped terms are under about 1%; with mixing near 1 they
are as large as the kept ones, and the result is no longer trustworthy. The norm shows this directly:
⟨ψ̃|ψ̃⟩ = 1 + Σc², which is 1 + 1e-5 at n = 0 but 1 + 4.1 at n = 59.

Setting 3|ξ|(n+1)^{3/2} = 0.1 gives the highest reliable level, `n_valid`: **n ≈ 9 (GS) and 8 (ES)**. The 0.1 is
a convention, the same number as the |ξ| > 0.1 rule.

### 6.4 Fidelity: how close the corrected state is to the truth

Fidelity is |⟨ψ̃|ψ_exact⟩|, the overlap of the normalized corrected state with the exact eigenstate of the
same cubic potential (computed by `fd`). 1 means identical; 1 − fidelity is the infidelity. For the GS:

| n | mixing | norm − 1 | 1 − fidelity, harmonic | 1 − fidelity, first-order PT |
| --- | --- | --- | --- | --- |
| 0 | 0.003 | 1e-5 | 5.2e-6 | 2.5e-10 |
| 9 | 0.098 | 0.017 | 8.4e-3 | 1.5e-5 |
| 20 | 0.30 | 0.17 | 8.2e-2 | 1.3e-3 |
| 40 | 0.81 | 1.3 | 0.56 | 0.079 |
| 59 | 1.44 | 4.1 | 0.71 | 0.54 |

The correction helps enormously at low n, and its accuracy degrades smoothly above `n_valid`. Even at high n it
stays closer to the exact states than the harmonic functions do, but by n ≈ 50–60 both are poor, and some
high-level F_ij exceed 1, which a true overlap cannot.

### 6.5 What the warnings mean

- `|xi| > 0.1`: the warning your spec asked for. It never fires here (ξ ≈ 1e-3).
- `mixing ... reaches X at n = N; it stays below 0.1 only up to n = n_valid`: fires here. It means first-order
  perturbation theory is an inaccurate description of the levels above `n_valid`.

**Neither warning changes the calculation.** Every level is still corrected, and the full ni × nf matrix goes
into the VOF. A warning is a flag to check whether those levels matter for your result. Here they barely do,
because the intensity sits in low levels:

| direction | T | lineshape weight from levels above `n_valid` |
| --- | --- | --- |
| emission | 5 K | 0.25 % |
| emission | 300 K | 0.7 % |
| absorption | 5 K | 2.3 % |
| absorption | 300 K | 4.1 % |

### 6.6 Validation

| test | result |
| --- | --- |
| ξ vs an all-SI evaluation | identical |
| b = 0 vs `hermite` | 4e-16 |
| norms vs 1 + Σc² | 4e-16 |
| error vs the exact cubic as b → 0 | falls as b^1.99 (F_ij) and b^2.00 (VOF): correct first order |
| size of the effect as b → 0 | falls as b^1.00 |
| wrong sign of b | 30× worse, so the mirror-frame sign is right |
| VOF vs the exact cubic solution | 3.3e-4 of the peak, against an effect of 3.5e-2 (≈ 99 % captured) |

### 6.7 Run it

```bash
python validation/compare_vof.py --test cubic                                  # vs harmonic pypl
python validation/compare_vof.py --reference fd --fd-fit cubic --test cubic    # vs the exact cubic
python validation/plot_ccd_wavefunctions.py --wf-method cubic --vof-method pypl cubic
```

Flags: `--cubic-points {all,parabola}`, `--cubic-free-a`, `--cubic-normalize`. Report section
`07_cubic_perturbation/` has the fits, the validation figure, the wavefunctions (harmonic, PT and exact),
the F_ij error maps and the VOF against the exact solution.

In `wavefunctions_<state>.png`:
- **PT − exact** is ψ̃ₙ minus the exact eigenstate of the same cubic: the error of the correction.
- **harmonic − exact** is the full anharmonic effect.

The PT curves there are not renormalized, so at high n part of PT − exact is just excess norm.

**When to use.** As a cheap, analytic estimate of cubic anharmonicity when the important levels stay below
`n_valid`: low temperature, small b. Otherwise use `fd --fd-fit cubic`, which solves the same potential
exactly, or `fd` on a spline, which follows the DFT data instead of a cubic model.

---

## 7. Choosing a method

| method | potential | approximation | speed | use it for |
| --- | --- | --- | --- | --- |
| `pypl` | harmonic | none (float64 drift 1.7e-8; unstable rows i >~ 10 at S ~ 50+) | 0.01 s | the reference at moderate S; production harmonic spectra |
| `hermite` | harmonic | quadrature (exact to 1e-16) | ~3 s | validating `pypl`; the reference at large S (MgO); harmonic wavefunctions |
| `fd` + `parabola` | harmonic | grid, O(h⁴) | ~7 s | validating both of the above |
| `fd` + `spline3-notaknot-smooth` | anharmonic, follows the DFT data, C² everywhere | grid; tails continue the spline | ~7–12 s | **the recommended realistic anharmonic VOF** (§5.1) |
| `fd` + `spline` | anharmonic, follows the DFT data | grid; harmonic tails with a slope kink | ~7 s | older results; the kink moves the diamond VOF by ~2% |
| `fd` + `cubic` | cubic model | grid; box well inside the turnover | ~7 s | exact solution of the cubic model |
| `cubic` | cubic model | first-order perturbation | ~2 s | quick anharmonic estimate for low levels |

Two limits apply to every anharmonic option:
- **Line positions stay harmonic.** `bulid_fc_lsp` puts the lines at harmonic energies, so anharmonicity
  changes only the intensities.
- **The DFT data cover only Q = −0.14 … 0.83** (diamond; MgO: Q = −3.10 … 6.20 for both states).
  Anything a fit does beyond that is a model choice, not data.
