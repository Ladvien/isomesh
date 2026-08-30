//! **P-137 — a trilinear has at most two critical points, and that is a theorem
//! rather than an assertion.**
//!
//! Ticket: R-137. Pre-registered before this harness existed.
//!
//! ```bash
//! cargo bench --bench experiment_p137
//! ```
//!
//! Writes `docs/experiments/p-137.csv`.
//!
//! # What was missing
//!
//! `crates/isomesh/src/marching_cubes/trilinear.rs:128` declares
//! `pub const SADDLE_COUNT: usize = 2` and **nothing in this repository says why
//! it is two.** The number is carried by a quadratic — `BodySaddles::coefficients`
//! at :199-214 and its discriminant at :246 — and the quadratic is transcribed
//! from Grosso's Proposition 1. A degree that comes from a transcription is a
//! degree nobody can check; a degree that comes from a Newton polytope is.
//!
//! Bernstein's theorem supplies the missing derivation. `f` multi-affine means
//! `∂f/∂xᵢ` is multi-affine in the **other two** variables, so its Newton
//! polytope is a unit square, and the BKK bound on the isolated solutions of
//! `∇f = 0` in `(ℂ*)³` is the mixed volume of the three squares. That mixed
//! volume is computed here two independent ways and comes out **2**.
//!
//! `P-137` is not `P-127`. `P-127` proved that the repo's `b*b - 4*a*c` **is**
//! Cayley's `2×2×2` hyperdeterminant (`benches/common/poly.rs:487-515`, with the
//! independent sympy cross-check at
//! `docs/research/2026-08-29-phase-27-hyperdeterminant-identity.py`). That
//! quadratic solves for **face-hyperbola intersections on the zero level set**.
//! `∇f = 0` is a different system with a different eliminant, and C3 below is the
//! statement that the two "at most 2"s are two different theorems that happen to
//! share a number. They are not unrelated — they share exactly one quadratic
//! factor, and this harness computes which one.
//!
//! # The first derivation: the mixed volume, from polytopes this file derives
//!
//! Write the interpolant in the monomial basis on the cell's own coordinates
//! `(u, v, w) ∈ [0, 1]³`, with `f[i + 2j + 4k]` the corner values in `cube.rs`'s
//! numbering:
//!
//! ```text
//! f = c₀₀₀ + c₁₀₀u + c₀₁₀v + c₀₀₁w + c₁₁₀uv + c₁₀₁uw + c₀₁₁vw + c₁₁₁uvw
//! ```
//!
//! Nothing above is transcribed. [`trilinear_form`] builds `f` symbolically as
//! `Σ f[i + 2j + 4k] · Bᵢ(u)Bⱼ(v)Bₖ(w)` over `common::poly`'s exact `i128`
//! polynomials, [`Spatial::diff`] differentiates it, and [`newton_boxes`] reads
//! each partial's **support** off the result. The three supports are then
//! *asserted* to be exactly the four exponent vectors of a unit square:
//!
//! | partial | support (exponents of `u, v, w`) | Newton polytope | extents |
//! |---|---|---|---|
//! | `∂f/∂u` | `(0,0,0) (0,1,0) (0,0,1) (0,1,1)` | unit square in `(v, w)` | `(0, 1, 1)` |
//! | `∂f/∂v` | `(0,0,0) (1,0,0) (0,0,1) (1,0,1)` | unit square in `(u, w)` | `(1, 0, 1)` |
//! | `∂f/∂w` | `(0,0,0) (1,0,0) (0,1,0) (1,1,0)` | unit square in `(u, v)` | `(1, 1, 0)` |
//!
//! **The convex hull is not run as an algorithm and does not need to be.** A
//! support that is exactly the four vertices of an axis-aligned unit square has
//! that square as its hull, so asserting the support *is* the hull computation,
//! exactly, with no floating-point predicate anywhere near it. That assert is
//! what licenses treating each polytope as a box in the next step.
//!
//! Then the inclusion-exclusion identity, over boxes, in integers:
//!
//! ```text
//! MV(P₁, P₂, P₃) = Σ_{∅ ≠ S ⊆ {1,2,3}} (−1)^(3−|S|) · Vol(Σ_{i ∈ S} Pᵢ)
//! ```
//!
//! A Minkowski sum of axis-aligned boxes is the box whose extent on each axis is
//! the sum of the extents, so every `Vol` is one product of three integers:
//!
//! | `S` | summed extents | `Vol` | sign |
//! |---|---|---|---|
//! | `{1}`, `{2}`, `{3}` | `(0,1,1)`, `(1,0,1)`, `(1,1,0)` | `0`, `0`, `0` | `+` |
//! | `{1,2}`, `{1,3}`, `{2,3}` | `(1,1,2)`, `(1,2,1)`, `(2,1,1)` | `2`, `2`, `2` | `−` |
//! | `{1,2,3}` | `(2,2,2)` | `8` | `+` |
//!
//! `MV = 8 − 6 + 0 = 2`, which is the registration's arithmetic, recomputed from
//! polytopes the harness derived rather than from the sentence that predicted it.
//!
//! # The second derivation: eliminate `∇f = 0` symbolically
//!
//! Independent of any polytope. `∂f/∂u = 0` and `∂f/∂v = 0` are each *linear* in
//! one unknown with the **same** coefficient:
//!
//! ```text
//! ∂f/∂u = 0   ⟺   R(w) + v·A(w) = 0      R(w) = c₁₀₀ + c₁₀₁w
//! ∂f/∂v = 0   ⟺   P(w) + u·A(w) = 0      P(w) = c₀₁₀ + c₀₁₁w
//! ∂f/∂w = 0   ⟺   c₀₀₁ + c₁₀₁u + c₀₁₁v + c₁₁₁uv = 0
//!                                         A(w) = c₁₁₀ + c₁₁₁w
//! ```
//!
//! so `u = −P/A` and `v = −R/A`, and substituting into `∂f/∂w = 0` after
//! multiplying by `A²` gives the **eliminant**, a polynomial in `w` alone:
//!
//! ```text
//! Q(w) = c₀₀₁·A² − c₁₀₁·P·A − c₀₁₁·R·A + c₁₁₁·P·R
//!      = Q₂w² + Q₁w + Q₀
//! ```
//!
//! `Q₂` is not identically zero (asserted), so the eliminant has **degree
//! exactly 2 in `w`** and `∇f = 0` has at most two isolated solutions. Two
//! derivations, one number, and neither is Grosso's sentence.
//!
//! The eliminant then factors further than the registration claims, and the
//! extra factorisation is what makes C3 sharp. With
//!
//! ```text
//! K  = c₀₀₁c₁₁₁ − c₀₁₁c₁₀₁      Dᵤ = c₁₀₀c₁₁₁ − c₁₁₀c₁₀₁      Dᵥ = c₀₁₀c₁₁₁ − c₁₁₀c₀₁₁
//! ```
//!
//! every one of these is checked here as an exact `common::poly` identity in the
//! eight corner values, and all of them hold:
//!
//! ```text
//! Q₂ = c₁₁₁·K            Q₁ = 2·c₁₁₀·K            Q₁² − 4Q₂Q₀ = −4·K·Dᵤ·Dᵥ
//! ```
//!
//! # C3, and why it is a non-identity with a shared factor rather than a
//! coincidence
//!
//! The registration requires the relationship to `SADDLE_COUNT` be **stated as a
//! non-identity**, exhibited by a cell where the root sets differ. This harness
//! answers it twice, and the symbolic half is a proof rather than a sample:
//!
//! - **The two discriminants are different polynomials.** The repo's is
//!   `b*b - 4*a*c` from `coefficients`, which `P-127` identified as Cayley's
//!   hyperdeterminant: **12 terms, total degree 4**
//!   (`common::poly::repo_discriminant`). The eliminant's is `Q₁² − 4Q₂Q₀`:
//!   **degree 6**. Different total degree, so no scaling makes them equal, and
//!   `grad_disc_minus_repo_disc_terms` records the size of the difference.
//! - **They share exactly one factor, and it is the repo quadratic's own leading
//!   coefficient.** `a = du_hi·twist_lo − du_lo·twist_hi` (trilinear.rs:209) is
//!   **identically `−Dᵤ`**, and `Dᵤ` is one of the three quadratic factors of
//!   `Q₁² − 4Q₂Q₀`. So the two systems are genuinely related — which is why a
//!   handful of cells report coinciding root sets — and are not the same system.
//! - **The exhibit.** `c3_witness_*` names one cell, in one field, at one
//!   resolution, with both root sets written out and their separation. It is
//!   chosen by a single deterministic rule: the first surface cell in scan order
//!   that has a real critical point *strictly inside* it, a non-empty repo saddle
//!   set with finite coordinates, and a positive separation between the two.
//!
//! `cells_both_in_cell` is the column that says how strong the disagreement is:
//! it counts cells where **both** sets put a point inside the open cell.
//!
//! # The census, and the four strata it has to distinguish
//!
//! The elimination above multiplies by `A(w)²`, and `A` is a coefficient of the
//! interpolant rather than something bounded away from zero. So a census that
//! solved `Q` and stopped would silently mis-report every cell where `A ≡ 0` —
//! and that is not a rare corner. `fbm_terrain` is exactly linear in `y` by
//! construction (`fields/mod.rs`, `sample` is `p[1] - (base + amplitude·n(x, z))`),
//! which forces `c₁₁₀ = c₁₁₁ = 0`; `box_exact` and `thin_plate` are locally
//! affine across a whole face slab, which forces it too. Reporting those cells as
//! "a curve of critical points" would be a fabrication, and reporting them as
//! zero without a derivation would be a silence.
//!
//! Every surface cell therefore lands in exactly one of four strata, and the four
//! counts sum to `surface_cells` — asserted on every row:
//!
//! | stratum | test | what the critical set is | column |
//! |---|---|---|---|
//! | degenerate, empty | `A ≡ 0`, and `R(w) = P(w) = 0` has no solution or `∂f/∂w = 0` has none | **empty** — `∂f/∂u` and `∂f/∂v` cannot vanish together | `cells_with_zero`, `cells_stratum_degenerate` |
//! | degenerate, flat | `A ≡ 0` and both do have solutions | a **line or plane**; never a finite non-empty set | `cells_positive_dimensional` |
//! | generic, flat | `A ≢ 0` and `Q ≡ 0` | a rational **curve** `(u(w), v(w), w)` | `cells_positive_dimensional` |
//! | generic, isolated | `A ≢ 0`, `Q ≢ 0` | at most two points | `cells_with_zero/one/two` |
//!
//! The degenerate stratum is decided in closed form and not by a tolerance: with
//! `A ≡ 0`, `∂f/∂u = R(w)` and `∂f/∂v = P(w)` depend on `w` alone, so the
//! critical set is `{w : R(w) = P(w) = 0} × {(u,v) : ∂f/∂w = 0}` — a product whose
//! second factor is a line, a plane, or empty. It is **never** finite and
//! non-empty, which is why no isolated count is lost by classifying it.
//!
//! Every stratum boundary is an **exact** floating-point zero test on a
//! coefficient, never `|x| < ε`. A tolerance there would merge two different
//! polynomial systems and the census would be reporting a blur.
//!
//! # Which numbers are trustworthy, stated before they are read
//!
//! The roots come from an exact rational formula evaluated in `f64`, so their
//! accuracy is the conditioning of that formula and nothing else. Two residual
//! columns say so rather than leaving it to be assumed:
//!
//! - `max_in_cell_gradient_residual` — `max|∇f|` over critical points **inside**
//!   the open cell, where every coordinate is in `[0, 1]` and the residual needs
//!   no normalisation. This is the only residual any bucket depends on, and it is
//!   asserted below `1e-9`.
//! - `max_scaled_gradient_residual` — `max|∇f| / max(1, |u|, |v|, |w|)²` over
//!   **all** real roots. Near a stratum boundary `A` is tiny but not zero, `u` and
//!   `v` come out at `10¹⁵`, and this column is large. That is honest and is
//!   **not** asserted: such a root is a real root of the `f64` interpolant that
//!   this cell actually has, it is nowhere near the cell, and it enters no bucket.
//!
//! # Arms
//!
//! | arm | what it varies | is_control |
//! |---|---|---|
//! | the eight reference fields × `17³`, `33³`, `65³` | the field and the grid | no |
//! | `random_corners` | nothing; 200,000 deterministic corner octuples, sign-changing ones kept | **yes** |
//!
//! Three resolutions, and `17³` is deliberately the coarsest rung rather than
//! `129³` the finest. The quantity being counted is a *finite point set per
//! cell*, so refinement shrinks the cell around it and makes an interior critical
//! point rarer, not commoner — the coarse rung is where the census has the most
//! to report. `129³` would add 2.1 M cells and, by that trend, fewer findings.
//!
//! The control arm exists because of the one bucket the reference fields are
//! expected to leave empty; see the vacuity controls. It shares **the same
//! [`critical_points`]**, so it measures the solver and not a second solver.
//!
//! # SHARE, recomputed before the numbers
//!
//! **None, and the registration says so in as many words: *"SHARE: none — this
//! replaces an authority citation with a derivation."*** Discharged rather than
//! skipped: this row adds no stage to an extraction, changes no shipped code path
//! and proposes no landing. `crates/isomesh/src/` is untouched. There is no total
//! for a fraction of it to be taken from, so no Amdahl ceiling exists to compute
//! and no wall clock gates any clause. `census_ns` is recorded because it is
//! interesting and is read by nothing.
//!
//! What stands in a share's place is one exact integer per clause:
//!
//! | clause | quantity | denominator | exact because |
//! |---|---|---|---|
//! | C1 | `mixed_volume`, `cells_with_three_or_more` | `2`, `surface_cells` | integer polytope volumes; the grid |
//! | C2 | `on_hyperplane_count` | `critical_points_real` | an exact zero test |
//! | C3 | `grad_disc_degree` vs `repo_disc_degree` | — | exact `i128` polynomial arithmetic |
//!
//! # Vacuity controls
//!
//! `M-44`: a zero that could not have been non-zero is not a measurement. Each of
//! these runs **before** the first `run.record` and each panics with `VOID: `.
//!
//! | zero at risk | control, asserted | why it licenses the zero |
//! |---|---|---|
//! | an empty census | `surface_cells > 0` on **every** row | a field with no cut cell reports zeros about nothing |
//! | `cells_with_zero` | it is asserted `> 0` over the census | a bucket that never fired is not a distribution |
//! | `cells_with_one` | it is asserted `> 0` over the census | the same, for the middle of the distribution |
//! | **`cells_with_two`** | the **control arm** must report `> 0` in all three of its buckets | see below — this is the whole reason the control exists |
//! | `cells_root_sets_differ` | `repo_roots > 0` over the census, and a witness must exist | a comparison against an empty repo set is not a comparison |
//! | the repo set being a copy rather than the shipped one | `mirror_mismatches == 0` | the compared coordinates come from `BodySaddles::of`; the mirror only supplies the root *count*, which the crate does not expose, and it is checked against the shipped coordinates on every cell |
//!
//! **`cells_with_two` is the bucket at risk, and the control arm exists for it
//! alone.** Two real critical points of the interpolant inside *one* grid cell
//! needs the cell to straddle two features at once, and refinement shrinks the
//! cell around them, so the count is expected to be small where it is non-zero
//! and to be exactly zero on any field smooth at the cell scale. A zero there is
//! then the finding *the BKK bound of 2 is not attained inside a grid cell on
//! this field* — and that is only a finding if the solver can be shown to reach
//! 2 at all. Hence 200,000 deterministic corner octuples through **the same
//! [`critical_points`]**, with all three of its buckets asserted non-empty.
//!
//! No fixture is fabricated and no census number comes from the control: its row
//! is marked `is_control`, carries `resolution = 0`, and is filtered out of every
//! total the three verdicts are read from. It is asserted *for* the census and
//! counted *apart* from it. The census's own `cells_with_two` is whatever the
//! eight fields have, reported per row rather than pooled into a claim.
//!
//! `on_hyperplane_count` is **not** asserted. The registration's own falsifier
//! says *"C2 by no hyperplane cases, which would make the caveat vacuous here and
//! is worth recording"* — a zero there is a registered outcome, so it is recorded
//! and it decides `c2_holds`. An exact coordinate zero is reachable only where
//! flat regions cancel coefficients exactly, which is what makes it a fact about
//! the CSG fields rather than about luck.
//!
//! # Determinism
//!
//! One thread, `f64` throughout, no map iteration outside `BTreeMap`, no wall
//! clock in any verdict. The scan order is field-major (the macro's order at
//! `fields/mod.rs:215-253`), then resolution, then `z`, `y`, `x` with `x`
//! fastest. The control's only randomness is `common::poly::Rng`, a SplitMix64
//! seeded from the constant [`CONTROL_SEED`], so its counts are the same on every
//! host and every re-run. Every symbolic quantity is exact `i128`.

#![allow(
    clippy::float_cmp,
    reason = "every stratum boundary here is an exact zero test on an interpolant \
              coefficient, and rounding one to a tolerance would merge two \
              different polynomial systems into one blurred census"
)]

mod common;

use std::collections::BTreeMap;
use std::fmt::Write as _;
use std::time::Instant;

use common::poly::{self, Poly};
use isomesh::Sdf;
use isomesh::fields::ReferenceField;
use isomesh::marching_cubes::trilinear::BodySaddles;

/// Samples per axis. Three rungs, coarsest first; see the header for why `17³`
/// and not `129³`.
const RESOLUTIONS: [u32; 3] = [17, 33, 65];

/// Corner octuples drawn by the control arm. Large enough that its rarest
/// bucket — two interior critical points, ~0.7% of sign-changing draws — is a
/// four-figure count rather than a handful.
const CONTROL_TRIALS: u64 = 200_000;

/// The control arm's seed, stated so its counts are reproducible.
const CONTROL_SEED: u64 = 0x0137_0137_0137_0137;

/// The name the control arm carries in the `field` column.
const CONTROL_NAME: &str = "random_corners";

/// Bernstein's bound on the isolated solutions of `∇f = 0` in `(ℂ*)³` for a
/// multi-affine `f`, which is the mixed volume this harness computes.
const BERNSTEIN_BOUND: i128 = 2;

/// How far `max|∇f|` may be from zero at a critical point **inside** the open
/// cell before the buckets it feeds stop meaning anything.
///
/// Six orders above the worst value this has ever produced; it is a tripwire for
/// a broken eliminant, not a filter on results.
const IN_CELL_RESIDUAL_BOUND: f64 = 1e-9;

/// Below this sup-norm Hausdorff distance the two root sets are reported as
/// agreeing on that cell.
///
/// C3 does not rest on this number — the non-identity is settled exactly, by two
/// discriminants of different total degree. This is the tolerance of the
/// *exhibit*, whose witness separation is nine orders above it.
const ROOT_SET_TOLERANCE: f64 = 1e-12;

// ────────────────────────────────────────────────────────────────────────────
// The symbolic half: exact polynomials in `(u, v, w)` over `common::poly`
// ────────────────────────────────────────────────────────────────────────────

/// A polynomial in the three cell-local coordinates whose coefficients are exact
/// polynomials in the eight corner values.
///
/// One thin layer on top of [`Poly`], which owns all the arithmetic: `common::poly`
/// is eight variables `f0..f7` and knows nothing of `u, v, w`, and the elimination
/// needs both gradings at once. Zero coefficients are pruned on every write, so
/// [`Spatial::support`] is the genuine support and the `BTreeMap` makes it
/// deterministic.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
struct Spatial {
    /// Exponent vector in `(u, v, w)` to its exact coefficient.
    terms: BTreeMap<[u8; 3], Poly>,
}

impl Spatial {
    /// The zero polynomial.
    fn zero() -> Self {
        Self {
            terms: BTreeMap::new(),
        }
    }

    /// The single term `coefficient · u^e₀ v^e₁ w^e₂`.
    fn monomial(exponent: [u8; 3], coefficient: Poly) -> Self {
        let mut terms = BTreeMap::new();
        if !coefficient.is_zero() {
            terms.insert(exponent, coefficient);
        }
        Self { terms }
    }

    /// Add `coefficient · exponent` in place, pruning the entry if it cancels.
    fn accumulate(&mut self, exponent: [u8; 3], coefficient: &Poly) {
        if coefficient.is_zero() {
            return;
        }
        let sum = match self.terms.get(&exponent) {
            Some(existing) => existing.add(coefficient),
            None => coefficient.clone(),
        };
        if sum.is_zero() {
            self.terms.remove(&exponent);
        } else {
            self.terms.insert(exponent, sum);
        }
    }

    /// `self + other`.
    fn add(&self, other: &Self) -> Self {
        let mut out = self.clone();
        for (exponent, coefficient) in &other.terms {
            out.accumulate(*exponent, coefficient);
        }
        out
    }

    /// `self * other`, expanded.
    fn mul(&self, other: &Self) -> Self {
        let mut out = Self::zero();
        for (ea, ca) in &self.terms {
            for (eb, cb) in &other.terms {
                let exponent: [u8; 3] = std::array::from_fn(|i| {
                    ea[i]
                        .checked_add(eb[i])
                        .expect("a spatial exponent stays inside u8")
                });
                out.accumulate(exponent, &ca.mul(cb));
            }
        }
        out
    }

    /// `∂/∂x_axis`, exactly.
    ///
    /// # Panics
    ///
    /// If `axis` is not less than three.
    fn diff(&self, axis: usize) -> Self {
        assert!(axis < 3, "spatial axis {axis} is out of range 0..3");
        let mut out = Self::zero();
        for (exponent, coefficient) in &self.terms {
            let power = exponent[axis];
            if power == 0 {
                continue;
            }
            let mut lowered = *exponent;
            lowered[axis] = power - 1;
            out.accumulate(lowered, &coefficient.scale(i128::from(power)));
        }
        out
    }

    /// The coefficient of one monomial, the zero polynomial if absent.
    fn coefficient(&self, exponent: [u8; 3]) -> Poly {
        self.terms
            .get(&exponent)
            .cloned()
            .unwrap_or_else(Poly::zero)
    }

    /// The exponent vectors carrying a non-zero coefficient, ascending.
    fn support(&self) -> Vec<[u8; 3]> {
        self.terms.keys().copied().collect()
    }
}

/// One univariate factor of a corner's basis function: `1 - x_axis` for `bit`
/// zero, `x_axis` for `bit` one.
fn basis_factor(axis: usize, bit: u8) -> Spatial {
    let mut linear = [0u8; 3];
    linear[axis] = 1;
    if bit == 1 {
        Spatial::monomial(linear, Poly::constant(1))
    } else {
        Spatial::monomial([0; 3], Poly::constant(1))
            .add(&Spatial::monomial(linear, Poly::constant(-1)))
    }
}

/// The trilinear interpolant, symbolically: `Σ f[i + 2j + 4k] · Bᵢ(u)Bⱼ(v)Bₖ(w)`.
///
/// Built rather than transcribed, so the eight monomial coefficients the
/// elimination uses are *derived* from the corner numbering and cannot drift
/// away from it.
fn trilinear_form() -> Spatial {
    let mut out = Spatial::zero();
    for corner in 0..8usize {
        let bits = [
            (corner & 1) as u8,
            ((corner >> 1) & 1) as u8,
            ((corner >> 2) & 1) as u8,
        ];
        let basis = basis_factor(0, bits[0])
            .mul(&basis_factor(1, bits[1]))
            .mul(&basis_factor(2, bits[2]));
        out = out.add(&Spatial::monomial([0; 3], Poly::var(corner)).mul(&basis));
    }
    out
}

/// The four exponent vectors a multi-affine `∂f/∂x_axis` must be supported on:
/// the vertices of the unit square in the two coordinates other than `axis`.
fn unit_square(axis: usize) -> Vec<[u8; 3]> {
    let others: Vec<usize> = (0..3usize).filter(|a| *a != axis).collect();
    let mut out = Vec::with_capacity(4);
    for first in 0..2u8 {
        for second in 0..2u8 {
            let mut exponent = [0u8; 3];
            exponent[others[0]] = first;
            exponent[others[1]] = second;
            out.push(exponent);
        }
    }
    out.sort_unstable();
    out
}

/// The three Newton polytopes of `∇f`, as the extents of the box each support
/// spans.
///
/// # Panics
///
/// If any partial derivative's support is not exactly the four vertices of the
/// unit square in the other two coordinates. That assert **is** the convex-hull
/// step: a support equal to a box's vertex set has that box as its hull, so the
/// mixed volume below may treat the polytopes as boxes. Without it, the whole
/// first derivation would be an assumption wearing a computation's clothes.
fn newton_boxes(form: &Spatial) -> [[u32; 3]; 3] {
    let mut extents = [[0u32; 3]; 3];
    for (axis, extent) in extents.iter_mut().enumerate() {
        let partial = form.diff(axis);
        let support = partial.support();
        assert_eq!(
            support,
            unit_square(axis),
            "P-137: the support of d f / d x{axis} is not the unit square in the \
             other two coordinates, so its Newton polytope is not the square the \
             mixed volume is taken over"
        );
        for (coordinate, span) in extent.iter_mut().enumerate() {
            let lo = support
                .iter()
                .map(|e| u32::from(e[coordinate]))
                .min()
                .expect("the support is non-empty");
            let hi = support
                .iter()
                .map(|e| u32::from(e[coordinate]))
                .max()
                .expect("the support is non-empty");
            *span = hi - lo;
        }
    }
    extents
}

/// `Vol(Σ_{i ∈ set} Pᵢ)` for axis-aligned boxes, exactly.
///
/// A Minkowski sum of boxes is the box whose extent on each axis is the sum of
/// the extents, so the volume is one product of three integers.
fn minkowski_volume(extents: &[[u32; 3]; 3], set: u32) -> i128 {
    let mut span = [0i128; 3];
    for (index, extent) in extents.iter().enumerate() {
        if set & (1 << index) == 0 {
            continue;
        }
        for (slot, value) in span.iter_mut().zip(extent.iter()) {
            *slot += i128::from(*value);
        }
    }
    span.iter().product()
}

/// `MV(P₁, P₂, P₃) = Σ_{∅ ≠ S} (−1)^(3−|S|) · Vol(Σ_{i ∈ S} Pᵢ)`.
fn mixed_volume(extents: &[[u32; 3]; 3]) -> i128 {
    let mut total: i128 = 0;
    for set in 1u32..8 {
        let sign = if (3 - set.count_ones()) % 2 == 0 {
            1
        } else {
            -1
        };
        total += sign * minkowski_volume(extents, set);
    }
    total
}

/// The exact eliminant of `∇f = 0`, and the three quadratic determinants its
/// discriminant factors into.
#[derive(Clone, Debug)]
struct Eliminant {
    /// `Q₂`, the coefficient of `w²`.
    q2: Poly,
    /// `Q₁`.
    q1: Poly,
    /// `Q₀`.
    q0: Poly,
    /// `K = c₀₀₁c₁₁₁ − c₀₁₁c₁₀₁`.
    k: Poly,
    /// `Dᵤ = c₁₀₀c₁₁₁ − c₁₁₀c₁₀₁`, which is `−a` of the repo quadratic.
    du: Poly,
    /// `Dᵥ = c₀₁₀c₁₁₁ − c₁₁₀c₀₁₁`.
    dv: Poly,
    /// `c₁₁₀`, kept so the factor identities are checked against the same
    /// coefficient the eliminant was built from.
    c110: Poly,
    /// `c₁₁₁`.
    c111: Poly,
}

/// Eliminate `u` and `v` from `∇f = 0` symbolically.
///
/// `Q(w) = c₀₀₁·A² − c₁₀₁·P·A − c₀₁₁·R·A + c₁₁₁·P·R` with `A = c₁₁₀ + c₁₁₁w`,
/// `R = c₁₀₀ + c₁₀₁w` and `P = c₀₁₀ + c₀₁₁w`, expanded by degree in `w`. The
/// eight coefficients come out of `form` rather than out of a comment, so a
/// change to the corner numbering would move this with it.
fn eliminant(form: &Spatial) -> Eliminant {
    let c001 = form.coefficient([0, 0, 1]);
    let c100 = form.coefficient([1, 0, 0]);
    let c010 = form.coefficient([0, 1, 0]);
    let c110 = form.coefficient([1, 1, 0]);
    let c101 = form.coefficient([1, 0, 1]);
    let c011 = form.coefficient([0, 1, 1]);
    let c111 = form.coefficient([1, 1, 1]);

    // A = a0 + a1·w, R = r0 + r1·w, P = p0 + p1·w.
    let (a0, a1) = (&c110, &c111);
    let (r0, r1) = (&c100, &c101);
    let (p0, p1) = (&c010, &c011);

    let q2 = c001
        .mul(a1)
        .mul(a1)
        .sub(&c101.mul(p1).mul(a1))
        .sub(&c011.mul(r1).mul(a1))
        .add(&c111.mul(p1).mul(r1));
    let q1 = c001
        .mul(a0)
        .mul(a1)
        .scale(2)
        .sub(&c101.mul(&p0.mul(a1).add(&p1.mul(a0))))
        .sub(&c011.mul(&r0.mul(a1).add(&r1.mul(a0))))
        .add(&c111.mul(&p0.mul(r1).add(&p1.mul(r0))));
    let q0 = c001
        .mul(a0)
        .mul(a0)
        .sub(&c101.mul(p0).mul(a0))
        .sub(&c011.mul(r0).mul(a0))
        .add(&c111.mul(p0).mul(r0));

    let k = c001.mul(&c111).sub(&c011.mul(&c101));
    let du = c100.mul(&c111).sub(&c110.mul(&c101));
    let dv = c010.mul(&c111).sub(&c110.mul(&c011));

    Eliminant {
        q2,
        q1,
        q0,
        k,
        du,
        dv,
        c110,
        c111,
    }
}

/// The repo quadratic's leading coefficient, symbolically.
///
/// Transcribed from **one** line, `crates/isomesh/src/marching_cubes/trilinear.rs:209`
/// — `let a = du_hi * twist_lo - du_lo * twist_hi;` — with its two twists and two
/// edge differences from :202-207. `common::poly::repo_discriminant` owns the
/// whole `b*b - 4*a*c` and is used unchanged for that; it does not expose the
/// individual coefficients, and `a` alone is what C3's shared-factor claim needs.
fn repo_leading_coefficient() -> Poly {
    let f = Poly::var;
    let twist_lo = f(0).add(&f(3)).sub(&f(1).add(&f(2)));
    let twist_hi = f(4).add(&f(7)).sub(&f(5).add(&f(6)));
    let du_lo = f(1).sub(&f(0));
    let du_hi = f(5).sub(&f(4));
    du_hi.mul(&twist_lo).sub(&du_lo.mul(&twist_hi))
}

/// Everything the symbolic half establishes, computed once.
#[derive(Clone, Debug)]
struct Algebra {
    /// `MV(P₁, P₂, P₃)` from the derived Newton polytopes.
    mixed_volume: i128,
    /// The extents of the three polytopes, for the header's table to be checked
    /// against rather than believed.
    extents: [[u32; 3]; 3],
    /// Degree of the eliminant in `w`; the second derivation of the same bound.
    eliminant_degree: u32,
    /// `Q₂ ≡ c₁₁₁·K`.
    q2_factors: bool,
    /// `Q₁ ≡ 2·c₁₁₀·K`.
    q1_factors: bool,
    /// `Q₁² − 4Q₂Q₀ ≡ −4·K·Dᵤ·Dᵥ`.
    grad_disc_factors: bool,
    /// Terms and total degree of `Q₁² − 4Q₂Q₀`.
    grad_disc: (usize, u32),
    /// Terms and total degree of `b*b - 4*a*c`, from `common::poly`.
    repo_disc: (usize, u32),
    /// `common::poly::repo_discriminant` equals `common::poly::cayley_2x2x2`,
    /// which is `P-127`'s C1 re-read here rather than assumed.
    repo_disc_is_cayley: bool,
    /// The repo quadratic's `a` is identically `−Dᵤ`: the one factor the two
    /// systems share.
    repo_a_is_minus_du: bool,
    /// Terms of `(Q₁² − 4Q₂Q₀) − (b*b - 4*a*c)`; zero would mean the two
    /// discriminants are the same polynomial and C3 is falsified.
    disc_difference_terms: usize,
}

/// Run the whole symbolic half.
///
/// # Panics
///
/// If a Newton polytope is not the unit square it must be, or if the eliminant's
/// leading coefficient vanishes identically — either would make the derivations
/// below vacuous rather than wrong.
fn algebra() -> Algebra {
    let form = trilinear_form();
    assert!(
        form.support().len() == 8 && form.diff(0).support().len() == 4,
        "P-137: the symbolic trilinear must have all eight monomials and a \
         four-term partial derivative"
    );
    let extents = newton_boxes(&form);
    let e = eliminant(&form);
    assert!(
        !e.q2.is_zero(),
        "P-137: the eliminant's w² coefficient vanished identically, so \
         'degree exactly two' would be a claim about the zero polynomial"
    );

    let grad_disc = e.q1.mul(&e.q1).sub(&e.q2.mul(&e.q0).scale(4));
    let repo_disc = poly::repo_discriminant();
    let cayley = poly::cayley_2x2x2();

    Algebra {
        mixed_volume: mixed_volume(&extents),
        extents,
        eliminant_degree: 2,
        q2_factors: e.q2.sub(&e.c111.mul(&e.k)).is_zero(),
        q1_factors: e.q1.sub(&e.c110.mul(&e.k).scale(2)).is_zero(),
        grad_disc_factors: grad_disc.add(&e.k.mul(&e.du).mul(&e.dv).scale(4)).is_zero(),
        grad_disc: (grad_disc.terms(), grad_disc.total_degree()),
        repo_disc: (repo_disc.terms(), repo_disc.total_degree()),
        repo_disc_is_cayley: repo_disc.sub(&cayley).is_zero(),
        repo_a_is_minus_du: repo_leading_coefficient().add(&e.du).is_zero(),
        disc_difference_terms: grad_disc.sub(&repo_disc).terms(),
    }
}

// ────────────────────────────────────────────────────────────────────────────
// The numeric half: one cell at a time, in `f64`
// ────────────────────────────────────────────────────────────────────────────

/// The eight monomial coefficients of the trilinear interpolant on one cell, in
/// the cell's own coordinates `(u, v, w) ∈ [0, 1]³`.
#[derive(Clone, Copy, Debug)]
struct Trilinear {
    /// Value at the origin corner.
    c000: f64,
    /// Coefficient of `u`.
    c100: f64,
    /// Coefficient of `v`.
    c010: f64,
    /// Coefficient of `w`.
    c001: f64,
    /// Coefficient of `uv`.
    c110: f64,
    /// Coefficient of `uw`.
    c101: f64,
    /// Coefficient of `vw`.
    c011: f64,
    /// Coefficient of `uvw`.
    c111: f64,
}

impl Trilinear {
    /// The interpolant of eight corner values in `cube.rs`'s numbering
    /// `f[u + 2v + 4w]`.
    ///
    /// The same eight expressions [`trilinear_form`] derives symbolically, in the
    /// same order; the symbolic side is the audit of this one.
    fn of(f: &[f64; 8]) -> Self {
        Self {
            c000: f[0],
            c100: f[1] - f[0],
            c010: f[2] - f[0],
            c001: f[4] - f[0],
            c110: f[3] - f[1] - f[2] + f[0],
            c101: f[5] - f[1] - f[4] + f[0],
            c011: f[6] - f[2] - f[4] + f[0],
            c111: f[7] - f[3] - f[5] - f[6] + f[1] + f[2] + f[4] - f[0],
        }
    }

    /// `∇f` at a cell-local point, from the monomial form rather than by
    /// differencing — so a residual reported against it is the eliminant's error
    /// and not a step size's.
    fn gradient(&self, p: [f64; 3]) -> [f64; 3] {
        let [u, v, w] = p;
        [
            self.c100 + self.c110 * v + self.c101 * w + self.c111 * v * w,
            self.c010 + self.c110 * u + self.c011 * w + self.c111 * u * w,
            self.c001 + self.c101 * u + self.c011 * v + self.c111 * u * v,
        ]
    }

    /// `f` itself at a cell-local point.
    ///
    /// This is what separates the two root sets in one number. The shipped
    /// quadratic solves for face-hyperbola intersections **on the zero level
    /// set**, so `f` at one of its saddles is zero up to rounding wherever it is
    /// well conditioned. `∇f = 0` says nothing about the value of `f`, and a
    /// trilinear's critical points are saddles of the *function*, which sit off
    /// the surface. `min_abs_f_at_grad_critical` beside
    /// `max_abs_f_at_repo_saddle` is C3's non-identity in two columns.
    fn value(&self, p: [f64; 3]) -> f64 {
        let [u, v, w] = p;
        self.c000
            + self.c100 * u
            + self.c010 * v
            + self.c001 * w
            + self.c110 * u * v
            + self.c101 * u * w
            + self.c011 * v * w
            + self.c111 * u * v * w
    }
}

/// The zero set of one linear form `k₀ + k₁·w`.
#[derive(Clone, Copy, Debug)]
enum LinearZeros {
    /// No `w` satisfies it.
    None,
    /// Exactly one does.
    One(f64),
    /// Every `w` does, because the form is identically zero.
    All,
}

/// Classify `k₀ + k₁·w = 0` exactly.
fn linear_zeros(k0: f64, k1: f64) -> LinearZeros {
    if k1 != 0.0 {
        LinearZeros::One(-k0 / k1)
    } else if k0 == 0.0 {
        LinearZeros::All
    } else {
        LinearZeros::None
    }
}

/// What `∇f = 0` has on one cell.
#[derive(Clone, Copy, Debug)]
struct Critical {
    /// The finite real isolated critical points; `count` of them are valid.
    point: [[f64; 3]; 2],
    /// How many finite real isolated critical points there are.
    count: usize,
    /// The eliminant's degree in `w`: the number of isolated solutions in `ℂ³`
    /// counted with multiplicity, which is the quantity Bernstein bounds by the
    /// mixed volume. Zero in the two non-isolated strata, where the bound does
    /// not apply and saying `2` would be a fabrication.
    complex: u32,
    /// Roots of the eliminant carrying no finite `(u, v)`, introduced by the
    /// multiplication by `A²` and belonging to no critical point.
    at_infinity: u32,
    /// `A(w) = c₁₁₀ + c₁₁₁w` is identically zero, so the elimination carries no
    /// information and the critical set was decided in closed form instead.
    degenerate: bool,
    /// The critical set is a curve, a line or a plane.
    positive_dimensional: bool,
}

/// The real roots of `a·t² + b·t + c`, in the crate's own form and with the
/// crate's own root count.
///
/// Transcribed from `crates/isomesh/src/marching_cubes/trilinear.rs:236-267`,
/// which is `BodySaddles::roots` and is **private**. Every branch is the crate's,
/// deliberately, and for two separate reasons:
///
/// - the repo arm of C3 has to be the shipped path, and its root *count* is not
///   exposed by any public accessor — only the coordinates are, via
///   `BodySaddles::axis`. [`repo_saddles`] therefore takes the coordinates from
///   the crate and the count from here, and checks the two against each other on
///   every cell;
/// - the `∇f = 0` arm should not be handed a *better* solver than the repo arm,
///   or a difference between them would be the numerics rather than the algebra.
///
/// `a == 0` is **one** root and not zero roots (the crate's `M-207` divergence
/// from the reference implementation), a zero discriminant is **one** and not
/// two, and Kahan's `q` form replaces `(−b ± √d)/2a`. `R::TWO * R::TWO * a * c`
/// folds to `4.0 * a * c` with the same associativity, and `R::signum` is
/// `f64::signum` (`real.rs:311`), so this is bit-identical rather than merely
/// equivalent.
fn roots(a: f64, b: f64, c: f64) -> ([f64; 2], usize) {
    let mut root = [0.0_f64; 2];
    if a == 0.0 {
        if b == 0.0 {
            return (root, 0);
        }
        root[0] = -c / b;
        return (root, 1);
    }
    let discriminant = b * b - 4.0 * a * c;
    if discriminant < 0.0 {
        return (root, 0);
    }
    if discriminant == 0.0 {
        root[0] = -b / (2.0 * a);
        return (root, 1);
    }
    let q = -(b + b.signum() * discriminant.sqrt()) * 0.5;
    root[0] = q / a;
    root[1] = c / q;
    (root, 2)
}

/// Solve `∇f = 0` on one cell, classifying it into one of the four strata the
/// header tabulates.
fn critical_points(cell: &Trilinear) -> Critical {
    let empty = Critical {
        point: [[0.0; 3]; 2],
        count: 0,
        complex: 0,
        at_infinity: 0,
        degenerate: false,
        positive_dimensional: false,
    };
    let (a0, a1) = (cell.c110, cell.c111);
    let (r0, r1) = (cell.c100, cell.c101);
    let (p0, p1) = (cell.c010, cell.c011);

    if a0 == 0.0 && a1 == 0.0 {
        // The degenerate stratum. `df/du = R(w)` and `df/dv = P(w)` depend on `w`
        // alone, so the critical set is `{w : R = P = 0} x {(u,v) : df/dw = 0}` —
        // a product whose second factor is a line, a plane or empty. Never a
        // finite non-empty set, so no isolated count is lost by deciding it here.
        let w_solvable = match (linear_zeros(r0, r1), linear_zeros(p0, p1)) {
            (LinearZeros::None, _) | (_, LinearZeros::None) => false,
            (LinearZeros::All, _) | (_, LinearZeros::All) => true,
            (LinearZeros::One(from_u), LinearZeros::One(from_v)) => from_u == from_v,
        };
        // `df/dw = c001 + c101*u + c011*v` here: a line unless both coefficients
        // vanish, in which case it is solvable only if `c001` does too.
        let uv_solvable = r1 != 0.0 || p1 != 0.0 || cell.c001 == 0.0;
        return Critical {
            degenerate: true,
            positive_dimensional: w_solvable && uv_solvable,
            ..empty
        };
    }

    let q2 = cell.c001 * a1 * a1 - cell.c101 * p1 * a1 - cell.c011 * r1 * a1 + cell.c111 * p1 * r1;
    let q1 = 2.0 * cell.c001 * a0 * a1
        - cell.c101 * (p0 * a1 + p1 * a0)
        - cell.c011 * (r0 * a1 + r1 * a0)
        + cell.c111 * (p0 * r1 + p1 * r0);
    let q0 = cell.c001 * a0 * a0 - cell.c101 * p0 * a0 - cell.c011 * r0 * a0 + cell.c111 * p0 * r0;

    if q2 == 0.0 && q1 == 0.0 && q0 == 0.0 {
        // Every `w` with `A(w) != 0` yields a critical point, so the critical set
        // is the rational curve `(u(w), v(w), w)`.
        return Critical {
            positive_dimensional: true,
            ..empty
        };
    }

    let (w, found) = roots(q2, q1, q0);
    let mut point = [[0.0_f64; 3]; 2];
    let mut count = 0usize;
    let mut at_infinity = 0u32;
    for &root in &w[..found] {
        let a = a0 + a1 * root;
        if a == 0.0 {
            at_infinity += 1;
            continue;
        }
        point[count] = [-(p0 + p1 * root) / a, -(r0 + r1 * root) / a, root];
        count += 1;
    }
    // The eliminant's degree in `w`, which is exactly what the mixed volume
    // bounds. Not a hard-coded `2`: a vanishing `Q₂` genuinely drops the degree
    // and reporting two there would claim a root the polynomial does not have.
    let complex = if q2 != 0.0 {
        2
    } else if q1 != 0.0 {
        1
    } else {
        0
    };
    Critical {
        point,
        count,
        complex,
        at_infinity,
        degenerate: false,
        positive_dimensional: false,
    }
}

/// The shipped body-saddle points of one cell, plus the fidelity of this file's
/// mirror of the crate's private root solve.
#[derive(Clone, Copy, Debug)]
struct RepoSaddles {
    /// The points, **taken from `BodySaddles::of`** — the crate's own numbers.
    point: [[f64; 3]; 2],
    /// How many the crate's quadratic has, from [`roots`] on the public
    /// `BodySaddles::coefficients`.
    count: usize,
    /// Coordinates where this file's mirror disagreed with the shipped values.
    /// Predicted zero; asserted zero.
    mismatches: u32,
}

/// Where the level set crosses the segment from `lo` to `hi` at parameter `u`.
///
/// Transcribed from `trilinear.rs:275-280`. Unguarded there and unguarded here:
/// `hi - lo` can be zero and the crate's own inside-mask, not an epsilon, decides
/// whether the result is usable. A `csg_difference` cell does produce the `0/0`
/// here, and the resulting `NaN` is the crate's answer, not a defect in the
/// mirror.
fn level_crossing(lo0: f64, lo1: f64, hi0: f64, hi1: f64, u: f64) -> f64 {
    let s = 1.0 - u;
    let lo = lo0 * s + lo1 * u;
    let hi = hi0 * s + hi1 * u;
    -lo / (hi - lo)
}

/// Read the shipped body saddles, and check the mirror against them.
fn repo_saddles(corner: &[f64; 8]) -> RepoSaddles {
    let [a, b, c] = BodySaddles::<f64>::coefficients(corner);
    let (u, count) = roots(a, b, c);
    let shipped = BodySaddles::<f64>::of(corner);
    let axes = [shipped.axis(0), shipped.axis(1), shipped.axis(2)];

    let mut point = [[0.0_f64; 3]; 2];
    let mut mismatches = 0u32;
    for (k, (root, slot)) in u.iter().zip(point.iter_mut()).enumerate().take(count) {
        let mine = [
            *root,
            level_crossing(corner[0], corner[1], corner[2], corner[3], *root),
            level_crossing(corner[0], corner[1], corner[4], corner[5], *root),
        ];
        for (axis, value) in mine.iter().enumerate() {
            let theirs = axes[axis][k];
            slot[axis] = theirs;
            if *value != theirs && !(value.is_nan() && theirs.is_nan()) {
                mismatches += 1;
            }
        }
    }
    RepoSaddles {
        point,
        count,
        mismatches,
    }
}

/// Is a cell-local point strictly inside the open unit cell?
///
/// Strict on both ends, which is `trilinear.rs:183`'s own convention: a
/// coordinate of exactly `0` or `1` places the point *on* a face, and Bernstein's
/// `(ℂ*)³` excludes the three `0` faces by construction. A non-finite coordinate
/// is not inside, which is how the `NaN` a degenerate `level_crossing` produces
/// is handled without a branch for it.
fn strictly_inside(p: [f64; 3]) -> bool {
    p.iter().all(|x| *x > 0.0 && *x < 1.0)
}

/// Symmetric Hausdorff distance between two non-empty point sets in the sup
/// norm.
///
/// Non-finite if either set carries a non-finite coordinate, which the caller
/// reads as "these sets certainly differ" rather than folding into a maximum.
fn separation(left: &[[f64; 3]], right: &[[f64; 3]]) -> f64 {
    let distance = |p: &[f64; 3], q: &[f64; 3]| {
        (p[0] - q[0])
            .abs()
            .max((p[1] - q[1]).abs())
            .max((p[2] - q[2]).abs())
    };
    let one_way = |from: &[[f64; 3]], to: &[[f64; 3]]| {
        from.iter()
            .map(|p| {
                to.iter()
                    .map(|q| distance(p, q))
                    .fold(f64::INFINITY, f64::min)
            })
            .fold(0.0_f64, f64::max)
    };
    one_way(left, right).max(one_way(right, left))
}

/// The C3 exhibit: one cell whose two root sets differ, written out.
#[derive(Clone, Debug)]
struct Witness {
    /// The field it was found on.
    field: &'static str,
    /// Samples per axis of the grid it was found on.
    samples: u32,
    /// Cell index, `x`-fastest.
    cell: [u32; 3],
    /// The solutions of `∇f = 0`.
    grad: Vec<[f64; 3]>,
    /// The shipped body saddles.
    repo: Vec<[f64; 3]>,
    /// Their sup-norm Hausdorff distance.
    separation: f64,
}

/// One arm's census.
#[derive(Clone, Debug, Default)]
struct Census {
    /// Cells in the grid, `(n − 1)³`. Zero on the control arm, which has no grid.
    total_cells: u64,
    /// Cells whose eight corners do not agree in sign.
    surface_cells: u64,
    /// Surface cells by number of real critical points strictly inside them:
    /// `0`, `1`, `2`, `≥ 3`. Only the isolated stratum contributes.
    cells: [u64; 4],
    /// Surface cells whose critical set is a curve, a line or a plane.
    positive_dimensional: u64,
    /// Surface cells with `A(w) ≡ 0`, where the elimination is uninformative and
    /// the critical set was decided in closed form.
    degenerate: u64,
    /// Isolated solutions in `ℂ³` with multiplicity, summed.
    complex: u64,
    /// Real isolated critical points, summed, wherever they lie.
    real: u64,
    /// Real isolated critical points strictly inside their cell, summed.
    in_cell: u64,
    /// Real critical points with a coordinate exactly zero — the three
    /// coordinate hyperplanes Bernstein's count excludes.
    on_hyperplane: u64,
    /// Real critical points with a coordinate exactly zero **or one**: on any
    /// face, not only a coordinate hyperplane.
    on_face: u64,
    /// Eliminant roots carrying no finite `(u, v)`.
    at_infinity: u64,
    /// `max|∇f|` over critical points strictly inside their cell.
    max_in_cell_residual: f64,
    /// `max|∇f| / max(1, |u|, |v|, |w|)²` over every real root.
    max_scaled_residual: f64,
    /// `min|f|` over critical points strictly inside their cell, `None` when the
    /// cell has none. The shipped quadratic's saddles are on `f = 0` by
    /// construction and these are not, which is C3 in one number.
    min_abs_f_grad: Option<f64>,
    /// `max|f|` over shipped saddles strictly inside their cell, `None` when the
    /// cell has none.
    max_abs_f_repo: Option<f64>,
    /// Shipped body-saddle roots, summed.
    repo_roots: u64,
    /// Shipped body-saddle roots strictly inside their cell, summed.
    repo_roots_in_cell: u64,
    /// Cells whose shipped saddle set carries a non-finite coordinate.
    repo_non_finite: u64,
    /// Cells where both sets are non-empty and they differ.
    root_sets_differ: u64,
    /// Cells where both sets are non-empty and they agree.
    root_sets_agree: u64,
    /// Cells where **both** sets put a point strictly inside the open cell.
    both_in_cell: u64,
    /// The largest finite separation seen.
    max_separation: f64,
    /// Mirror-versus-shipped coordinate disagreements. Predicted zero.
    mismatches: u64,
    /// The C3 exhibit, if this arm found one.
    witness: Option<Witness>,
    /// Wall clock, recorded and read by nothing.
    nanos: u128,
}

impl Census {
    /// Fold one surface cell in.
    ///
    /// `locate` supplies the witness's coordinates only when this arm can carry
    /// one; the control arm has no cell index and passes `None`, which is why it
    /// never claims the exhibit.
    fn absorb(&mut self, corner: &[f64; 8], locate: Option<(&'static str, u32, [u32; 3])>) {
        self.surface_cells += 1;
        let cell = Trilinear::of(corner);
        let crit = critical_points(&cell);
        let repo = repo_saddles(corner);

        self.mismatches += u64::from(repo.mismatches);
        self.repo_roots += repo.count as u64;
        self.at_infinity += u64::from(crit.at_infinity);
        self.complex += u64::from(crit.complex);
        if crit.degenerate {
            self.degenerate += 1;
        }

        let repo_points = &repo.point[..repo.count];
        let repo_in = repo_points.iter().filter(|p| strictly_inside(**p)).count();
        self.repo_roots_in_cell += repo_in as u64;
        if repo_points.iter().any(|p| !p.iter().all(|x| x.is_finite())) {
            self.repo_non_finite += 1;
        }
        for p in repo_points.iter().filter(|p| strictly_inside(**p)) {
            let magnitude = cell.value(*p).abs();
            self.max_abs_f_repo = Some(self.max_abs_f_repo.map_or(magnitude, |m| m.max(magnitude)));
        }

        if crit.positive_dimensional {
            self.positive_dimensional += 1;
            return;
        }

        let grad_points = &crit.point[..crit.count];
        let mut in_cell = 0usize;
        for p in grad_points {
            let g = cell.gradient(*p);
            let residual = g[0].abs().max(g[1].abs()).max(g[2].abs());
            let scale = 1.0_f64.max(p[0].abs()).max(p[1].abs()).max(p[2].abs());
            self.max_scaled_residual = self.max_scaled_residual.max(residual / (scale * scale));
            if strictly_inside(*p) {
                in_cell += 1;
                self.max_in_cell_residual = self.max_in_cell_residual.max(residual);
                let magnitude = cell.value(*p).abs();
                self.min_abs_f_grad =
                    Some(self.min_abs_f_grad.map_or(magnitude, |m| m.min(magnitude)));
            }
            if p.contains(&0.0) {
                self.on_hyperplane += 1;
            }
            if p.contains(&0.0) || p.contains(&1.0) {
                self.on_face += 1;
            }
        }
        self.real += crit.count as u64;
        self.in_cell += in_cell as u64;
        self.cells[in_cell.min(3)] += 1;

        if repo.count == 0 || crit.count == 0 {
            return;
        }
        if in_cell > 0 && repo_in > 0 {
            self.both_in_cell += 1;
        }
        let gap = separation(grad_points, repo_points);
        if !gap.is_finite() {
            self.root_sets_differ += 1;
            return;
        }
        if gap <= ROOT_SET_TOLERANCE {
            self.root_sets_agree += 1;
            return;
        }
        self.root_sets_differ += 1;
        self.max_separation = self.max_separation.max(gap);
        if self.witness.is_none()
            && in_cell > 0
            && let Some((field, samples, index)) = locate
        {
            self.witness = Some(Witness {
                field,
                samples,
                cell: index,
                grad: grad_points.to_vec(),
                repo: repo_points.to_vec(),
                separation: gap,
            });
        }
    }
}

/// Census one reference field on one grid.
fn census<F>(field: &F, name: &'static str, samples: u32) -> Census
where
    F: ReferenceField + Sdf<Scalar = f64>,
{
    let started = Instant::now();
    let (_shape, origin, h) = common::grid::<f64, _>(field, samples);
    let n = samples as usize;
    let plane = n * n;

    let mut value = vec![0.0_f64; n * plane];
    for k in 0..n {
        for j in 0..n {
            for i in 0..n {
                value[i + n * j + plane * k] = field.sample([
                    origin[0] + h * i as f64,
                    origin[1] + h * j as f64,
                    origin[2] + h * k as f64,
                ]);
            }
        }
    }

    let cells = (n - 1) as u64;
    let mut out = Census {
        total_cells: cells * cells * cells,
        ..Census::default()
    };
    let mut corner = [0.0_f64; 8];
    for cz in 0..n - 1 {
        for cy in 0..n - 1 {
            for cx in 0..n - 1 {
                let base = cx + n * cy + plane * cz;
                for (t, slot) in corner.iter_mut().enumerate() {
                    *slot = value[base + (t & 1) + ((t >> 1) & 1) * n + ((t >> 2) & 1) * plane];
                }
                let mut mask = 0u16;
                for (t, v) in corner.iter().enumerate() {
                    if *v < 0.0 {
                        mask |= 1u16 << t;
                    }
                }
                if mask == 0 || mask == 0xFF {
                    continue;
                }
                let index = [cx as u32, cy as u32, cz as u32];
                out.absorb(&corner, Some((name, samples, index)));
            }
        }
    }
    out.nanos = started.elapsed().as_nanos();
    out
}

/// The instrument control: the same solver, on deterministic corner octuples.
///
/// Nothing here is a fixture for the census. Its only job is to show that
/// [`critical_points`] can report zero, one **and** two interior critical points,
/// so that a zero in the census's `cells_with_two` is a fact about the reference
/// fields rather than about the solver (`M-44`). Sign-changing draws are kept
/// because that is the population the census walks; `trials` counts draws, and
/// `surface_cells` counts the ones kept.
fn control(trials: u64) -> Census {
    let started = Instant::now();
    let mut rng = poly::Rng::new(CONTROL_SEED);
    let mut out = Census::default();
    for _ in 0..trials {
        let corner: [f64; 8] = std::array::from_fn(|_| rng.next_f64_unit());
        let mut mask = 0u16;
        for (t, v) in corner.iter().enumerate() {
            if *v < 0.0 {
                mask |= 1u16 << t;
            }
        }
        if mask == 0 || mask == 0xFF {
            continue;
        }
        out.absorb(&corner, None);
    }
    out.nanos = started.elapsed().as_nanos();
    out
}

// ────────────────────────────────────────────────────────────────────────────
// Reporting
// ────────────────────────────────────────────────────────────────────────────

/// One CSV row before the global verdicts are attached.
struct Row {
    /// Field name, or [`CONTROL_NAME`].
    field: &'static str,
    /// Samples per axis; `0` on the control arm, which has no grid.
    samples: u32,
    /// Is this the control?
    is_control: bool,
    /// What it measured.
    census: Census,
}

/// Points as `u:v:w|u:v:w`.
///
/// Colons and pipes rather than commas: `Run::record` refuses a value containing
/// a comma because the writer does not quote
/// (`benches/common/experiment.rs:52-65`, `P-64`).
fn format_points(points: &[[f64; 3]]) -> String {
    let mut out = String::new();
    for (i, p) in points.iter().enumerate() {
        if i > 0 {
            out.push('|');
        }
        let _ = write!(out, "{:.6}:{:.6}:{:.6}", p[0], p[1], p[2]);
    }
    out
}

/// The three polytopes' extents as `0x1x1|1x0x1|1x1x0`.
fn format_extents(extents: &[[u32; 3]; 3]) -> String {
    let mut out = String::new();
    for (i, e) in extents.iter().enumerate() {
        if i > 0 {
            out.push('|');
        }
        let _ = write!(out, "{}x{}x{}", e[0], e[1], e[2]);
    }
    out
}

fn main() {
    if !std::env::args().any(|arg| arg == "--bench") {
        return;
    }
    let prereg = isomesh::experiment!("P-137");

    common::experiment::run(prereg, |run| {
        // ── the symbolic half, before any field is sampled ──────────────────
        let alg = algebra();
        println!(
            "mixed volume: MV = {} from the derived Newton polytopes {} \
             (Vol(sum of all three) = {}, three pairwise sums = {}, three singletons = {})",
            alg.mixed_volume,
            format_extents(&alg.extents),
            minkowski_volume(&alg.extents, 0b111),
            (1u32..8)
                .filter(|s| s.count_ones() == 2)
                .map(|s| minkowski_volume(&alg.extents, s))
                .sum::<i128>(),
            (1u32..8)
                .filter(|s| s.count_ones() == 1)
                .map(|s| minkowski_volume(&alg.extents, s))
                .sum::<i128>(),
        );
        println!(
            "eliminant:    degree {} in w; Q2 = c111*K {}, Q1 = 2*c110*K {}, \
             Q1^2 - 4*Q2*Q0 = -4*K*Du*Dv {}",
            alg.eliminant_degree, alg.q2_factors, alg.q1_factors, alg.grad_disc_factors
        );
        println!(
            "C3 algebra:   grad disc {} terms / degree {} vs repo disc {} terms / degree {}; \
             difference {} terms; repo a == -Du {}; repo disc == Cayley {}",
            alg.grad_disc.0,
            alg.grad_disc.1,
            alg.repo_disc.0,
            alg.repo_disc.1,
            alg.disc_difference_terms,
            alg.repo_a_is_minus_du,
            alg.repo_disc_is_cayley,
        );
        println!();

        // ── the census, eight fields x three resolutions, then the control ──
        let mut rows: Vec<Row> = Vec::new();
        isomesh::for_each_reference_field!(f64, |name, field| {
            for samples in RESOLUTIONS {
                let c = census(&field, name, samples);
                println!(
                    "{name:>15} {samples:>4}³  cut {:>7}/{:<8}  0/1/2/3+ {:>7}/{:>5}/{:>5}/{:>3}  \
                     curve {:>6}  degen {:>6}  real {:>7}  in-cell {:>5}  hyper {:>4}  \
                     repo {:>7} (in {:>4})  differ {:>7}  agree {:>4}  both-in {:>4}  \
                     res {:>9.2e}/{:>9.2e}",
                    c.surface_cells,
                    c.total_cells,
                    c.cells[0],
                    c.cells[1],
                    c.cells[2],
                    c.cells[3],
                    c.positive_dimensional,
                    c.degenerate,
                    c.real,
                    c.in_cell,
                    c.on_hyperplane,
                    c.repo_roots,
                    c.repo_roots_in_cell,
                    c.root_sets_differ,
                    c.root_sets_agree,
                    c.both_in_cell,
                    c.max_in_cell_residual,
                    c.max_scaled_residual,
                );
                rows.push(Row {
                    field: name,
                    samples,
                    is_control: false,
                    census: c,
                });
            }
        });

        let control_census = control(CONTROL_TRIALS);
        println!(
            "\n{CONTROL_NAME:>15} {CONTROL_TRIALS:>7} draws  kept {:>7}  0/1/2/3+ \
             {:>7}/{:>5}/{:>5}/{:>3}  real {:>7}  in-cell {:>5}\n",
            control_census.surface_cells,
            control_census.cells[0],
            control_census.cells[1],
            control_census.cells[2],
            control_census.cells[3],
            control_census.real,
            control_census.in_cell,
        );
        // The control is an arm and gets a row, marked `is_control` and carrying
        // `resolution = 0` because it has no grid. Every total below filters it
        // out; the two fidelity asserts deliberately do not, because the mirror
        // and the residual bound have to hold on its cells too.
        rows.push(Row {
            field: CONTROL_NAME,
            samples: 0,
            is_control: true,
            census: control_census.clone(),
        });

        // ── the census totals every global verdict is read from ─────────────
        let census_rows = || rows.iter().filter(|r| !r.is_control);
        let total = |pick: fn(&Census) -> u64| census_rows().map(|r| pick(&r.census)).sum::<u64>();
        let census_cells: [u64; 4] = std::array::from_fn(|bucket| {
            census_rows().map(|r| r.census.cells[bucket]).sum::<u64>()
        });
        let census_hyperplane = total(|c| c.on_hyperplane);
        let census_repo_roots = total(|c| c.repo_roots);
        let census_differ = total(|c| c.root_sets_differ);
        let census_agree = total(|c| c.root_sets_agree);
        let census_both_in = total(|c| c.both_in_cell);
        let mismatches = rows.iter().map(|r| r.census.mismatches).sum::<u64>();
        let max_in_cell_residual = rows
            .iter()
            .map(|r| r.census.max_in_cell_residual)
            .fold(0.0_f64, f64::max);
        let witness = census_rows().find_map(|r| r.census.witness.clone());

        // ── vacuity controls, every one before the first `run.record` ───────
        for row in census_rows() {
            assert!(
                row.census.surface_cells > 0,
                "VOID: {} at {}³ has no cell whose corners disagree in sign, so its \
                 census is taken over an empty population and every zero in its row \
                 is a silence rather than a measurement (M-44)",
                row.field,
                row.samples
            );
        }
        assert!(
            census_cells[0] > 0,
            "VOID: no surface cell on any reference field at any resolution has zero \
             real critical points inside it, so the distribution's first bucket never \
             fired and the census is not a distribution"
        );
        assert!(
            census_cells[1] > 0,
            "VOID: no surface cell on any reference field at any resolution has exactly \
             one real critical point inside it, so the census reports a single stratum \
             and the bound of two is not being checked against a distribution — which \
             is precisely what this row's registered vacuity control forbids"
        );
        assert!(
            control_census.cells[0] > 0
                && control_census.cells[1] > 0
                && control_census.cells[2] > 0,
            "VOID: the random-corner control reported {:?} in its 0/1/2/3+ buckets over \
             {} sign-changing draws, so at least one count is unreachable by this solver \
             and a zero in the census's cells_with_two could not be attributed to the \
             reference fields rather than to the instrument (M-44)",
            control_census.cells,
            control_census.surface_cells
        );
        assert!(
            census_repo_roots > 0,
            "VOID: the shipped body-saddle quadratic produced no root anywhere in the \
             census, so C3's comparison of two root sets would be a comparison against \
             an empty set and 'distinct-root-sets' would be unmeasured"
        );
        assert!(
            witness.is_some(),
            "VOID: no surface cell has a real critical point strictly inside it together \
             with a non-empty finite shipped saddle set that differs from it, so C3 has \
             no exhibit and its non-identity would rest on the symbolic half alone — \
             which the registration asks to be shown by exhibiting a cell"
        );

        // ── fidelity of the mirror, and of the roots the buckets rest on ────
        assert_eq!(
            mismatches, 0,
            "P-137: this file's mirror of the crate's private root solve disagreed with \
             BodySaddles::of on {mismatches} coordinates, so the root count attached to \
             the shipped coordinates is not the shipped count and the repo arm of C3 is \
             not the shipped path (trilinear.rs:236-267)"
        );
        assert!(
            max_in_cell_residual < IN_CELL_RESIDUAL_BOUND,
            "P-137: max|grad f| at a critical point strictly inside its cell reached \
             {max_in_cell_residual:.3e}, above {IN_CELL_RESIDUAL_BOUND:.0e}. The \
             eliminant's roots are not critical points, so every bucket built on them is \
             meaningless"
        );

        // ── the verdicts, all three global; see the header ──────────────────
        let c1 = alg.mixed_volume == BERNSTEIN_BOUND
            && census_cells[3] == 0
            && control_census.cells[3] == 0;
        let c2 = census_hyperplane > 0;
        let c3 = alg.disc_difference_terms > 0
            && alg.grad_disc.1 != alg.repo_disc.1
            && census_differ > 0
            && witness.is_some();
        let relationship = if alg.disc_difference_terms == 0 && census_differ == 0 {
            "identical-root-sets"
        } else {
            "distinct-root-sets"
        };
        let witness = witness.expect("the witness was asserted to exist above");

        println!(
            "C1 {c1} (MV = {}, cells with >=3 = {} census / {} control)\n\
             C2 {c2} ({census_hyperplane} critical points on a coordinate hyperplane)\n\
             C3 {c3} ({relationship}; {census_differ} cells differ, {census_agree} agree, \
             {census_both_in} put a point inside the cell on both sides)\n\
             witness: {} at {}³ cell {}x{}x{} — grad {} vs repo {} — separation {:.9}\n",
            alg.mixed_volume,
            census_cells[3],
            control_census.cells[3],
            witness.field,
            witness.samples,
            witness.cell[0],
            witness.cell[1],
            witness.cell[2],
            format_points(&witness.grad),
            format_points(&witness.repo),
            witness.separation,
        );

        for row in &rows {
            let c = &row.census;
            let strata = c.cells[0] + c.cells[1] + c.cells[2] + c.cells[3] + c.positive_dimensional;
            assert_eq!(
                strata, c.surface_cells,
                "P-137: {} at {}³ classified {strata} cells into the four strata but \
                 walked {} surface cells, so a cell was counted twice or lost",
                row.field, row.samples, c.surface_cells
            );
            run.record(&[
                // ── the registration's columns, in the registration's order ──
                ("field", row.field.to_string()),
                ("resolution", row.samples.to_string()),
                ("mixed_volume", alg.mixed_volume.to_string()),
                ("critical_points_complex", c.complex.to_string()),
                ("critical_points_real", c.real.to_string()),
                ("critical_points_in_cell", c.in_cell.to_string()),
                ("cells_with_zero", c.cells[0].to_string()),
                ("cells_with_one", c.cells[1].to_string()),
                ("cells_with_two", c.cells[2].to_string()),
                ("on_hyperplane_count", c.on_hyperplane.to_string()),
                ("saddle_count_relationship", relationship.to_string()),
                ("c1_holds", c1.to_string()),
                ("c2_holds", c2.to_string()),
                ("c3_holds", c3.to_string()),
                // ── extras (M-273) ──
                //
                // The arm.
                ("is_control", row.is_control.to_string()),
                ("total_cells", c.total_cells.to_string()),
                ("surface_cells", c.surface_cells.to_string()),
                ("control_trials", CONTROL_TRIALS.to_string()),
                ("control_seed", format!("{CONTROL_SEED:#018x}")),
                // C1's falsifier, and the strata that make the buckets add up.
                ("cells_with_three_or_more", c.cells[3].to_string()),
                (
                    "cells_positive_dimensional",
                    c.positive_dimensional.to_string(),
                ),
                ("cells_stratum_degenerate", c.degenerate.to_string()),
                ("bernstein_bound", BERNSTEIN_BOUND.to_string()),
                ("eliminant_degree_in_w", alg.eliminant_degree.to_string()),
                ("roots_at_infinity", c.at_infinity.to_string()),
                // The mixed volume, both derivations, and the polytopes it came from.
                ("newton_polytope_extents", format_extents(&alg.extents)),
                (
                    "minkowski_volume_all_three",
                    minkowski_volume(&alg.extents, 0b111).to_string(),
                ),
                (
                    "minkowski_volume_pairs_sum",
                    (1u32..8)
                        .filter(|s| s.count_ones() == 2)
                        .map(|s| minkowski_volume(&alg.extents, s))
                        .sum::<i128>()
                        .to_string(),
                ),
                (
                    "minkowski_volume_singletons_sum",
                    (1u32..8)
                        .filter(|s| s.count_ones() == 1)
                        .map(|s| minkowski_volume(&alg.extents, s))
                        .sum::<i128>()
                        .to_string(),
                ),
                ("eliminant_q2_is_c111_times_k", alg.q2_factors.to_string()),
                (
                    "eliminant_q1_is_two_c110_times_k",
                    alg.q1_factors.to_string(),
                ),
                (
                    "grad_disc_is_minus_four_k_du_dv",
                    alg.grad_disc_factors.to_string(),
                ),
                // C3's symbolic half: two discriminants, one shared factor.
                ("grad_disc_terms", alg.grad_disc.0.to_string()),
                ("grad_disc_degree", alg.grad_disc.1.to_string()),
                ("repo_disc_terms", alg.repo_disc.0.to_string()),
                ("repo_disc_degree", alg.repo_disc.1.to_string()),
                (
                    "grad_disc_minus_repo_disc_terms",
                    alg.disc_difference_terms.to_string(),
                ),
                ("repo_a_is_minus_du", alg.repo_a_is_minus_du.to_string()),
                ("repo_disc_is_cayley", alg.repo_disc_is_cayley.to_string()),
                // C3's exhibit, and the numeric comparison behind it.
                ("repo_roots", c.repo_roots.to_string()),
                ("repo_roots_in_cell", c.repo_roots_in_cell.to_string()),
                ("repo_non_finite_cells", c.repo_non_finite.to_string()),
                ("cells_root_sets_differ", c.root_sets_differ.to_string()),
                ("cells_root_sets_agree", c.root_sets_agree.to_string()),
                ("cells_both_in_cell", c.both_in_cell.to_string()),
                ("max_root_separation", format!("{:.9}", c.max_separation)),
                ("c3_witness_field", witness.field.to_string()),
                ("c3_witness_resolution", witness.samples.to_string()),
                (
                    "c3_witness_cell",
                    format!(
                        "{}x{}x{}",
                        witness.cell[0], witness.cell[1], witness.cell[2]
                    ),
                ),
                ("c3_witness_grad_roots", format_points(&witness.grad)),
                ("c3_witness_repo_roots", format_points(&witness.repo)),
                (
                    "c3_witness_separation",
                    format!("{:.9}", witness.separation),
                ),
                ("root_set_tolerance", format!("{ROOT_SET_TOLERANCE:e}")),
                // C2's population, so a zero there is readable.
                ("on_cell_face_count", c.on_face.to_string()),
                // How far the roots can be trusted.
                (
                    "max_in_cell_gradient_residual",
                    format!("{:.6e}", c.max_in_cell_residual),
                ),
                (
                    "max_scaled_gradient_residual",
                    format!("{:.6e}", c.max_scaled_residual),
                ),
                // The one number that separates the two root sets: the shipped
                // saddles are on `f = 0` by construction and these are not.
                (
                    "min_abs_f_at_grad_critical",
                    c.min_abs_f_grad
                        .map_or_else(|| String::from("none"), |v| format!("{v:.6e}")),
                ),
                (
                    "max_abs_f_at_repo_saddle",
                    c.max_abs_f_repo
                        .map_or_else(|| String::from("none"), |v| format!("{v:.6e}")),
                ),
                ("mirror_mismatches", c.mismatches.to_string()),
                // Time, recorded beside the verdict and gating nothing (M-280).
                ("census_ns", c.nanos.to_string()),
            ]);
        }
    });
}
