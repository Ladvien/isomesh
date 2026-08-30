//! **P-158 — a null registered on purpose: sparse grids need a regularity our
//! sharp fields have not got.**
//!
//! Ticket: R-158. Pre-registered before this harness existed.
//!
//! ```bash
//! cargo bench --bench experiment_p158
//! ```
//!
//! Writes `docs/experiments/p-158.csv`, one row per reference field.
//!
//! # What was missing
//!
//! Sparse grids are the one method in Phase 27 whose *hypothesis* is cheaper to
//! test than the method, and the repository had never tested either. The corpus
//! audit records the whole field as absent — *"Sparse grids / Smolyak | ABSENT
//! (0.597) | 0 | Registered expecting a null"* (`FINDINGS.md:25051`,
//! `docs/research/2026-08-29-phase-27-fifty-experiments-from-unmined-mathematics.md:85`)
//! — and no bench, test or document in the tree contains the words `Smolyak`,
//! `H²_mix` or `mixed derivative`. So there was no number anywhere saying whether
//! the eight reference fields satisfy the smoothness the `O(h⁻¹·|log h|^(d−1))`
//! point count is conditional on.
//!
//! There was no *instrument* either. `experiment_p43` and `experiment_p44` fit a
//! log-log slope to a Hausdorff distance over the resolution sweep
//! (`experiment_p43.rs:394-430`, `decay_exponent`) and that shape is reused here,
//! but neither differentiates the field: both read the *mesh's* error, not the
//! field's regularity. `validate::field_bound_report` asks whether `|∇f| ≈ 1`,
//! which is a first-derivative question. Nothing in the crate looks at a second
//! derivative, let alone a mixed sixth one, and `crates/isomesh/src/**` is
//! read-only for this phase — so the differentiator is bench-local and driven
//! through `Sdf::sample` alone.
//!
//! Three facts already in the tree are load-bearing for the reading of the
//! numbers below and are cited rather than rediscovered:
//!
//! - **`FbmTerrain::sample` is affine in `y`** (`fields/mod.rs:1352-1362`:
//!   `p[1] − (base_height + amplitude·n([p0, 0, p2]))`, the noise never sees
//!   `p[1]`). So `∂²u/∂y² ≡ 0` and the mixed sixth derivative is *structurally*
//!   zero — the roughest surface in the roster is the one field that is
//!   trivially in `H²_mix`, and for a reason that has nothing to do with its
//!   surface.
//! - **`Gyroid::sample` is a sum of three two-variable terms**
//!   (`fields/mod.rs:1022-1024`: `sin a·cos b + sin b·cos c + sin c·cos a`), so
//!   each term is annihilated by the third second derivative and the *nodal*
//!   gyroid also has an identically zero mixed sixth derivative. What survives
//!   in `capped_gyroid` is the `max` against the bounding sphere
//!   (`fields/mod.rs:744-747`), i.e. the seam, not the surface.
//! - **`ThinPlate` is a box** (`fields/mod.rs:651`, straight to `box_sample`)
//!   `0.025` thick, sub-voxel at every resolution this crate measures on, and
//!   `M-266`/`M-72` recorded that its survival under LOD is *alignment* rather
//!   than chance (`fields/mod.rs:606-616`). Its edges are box edges, so C1 reads
//!   on it exactly as it reads on `box_exact`.
//!
//! # Arms
//!
//! | arm | what it varies | is_control |
//! |---|---|---|
//! | `mixed_derivative` | the resolution `n ∈ {33, 49, 65, 97}`, per field | no |
//! | `stencil_smooth` | nothing — `x²y²z²`, whose exact `∂⁶` is `8` | **yes** |
//! | `stencil_kink` | nothing — `max(x+y+z, 0)`, whose exact stencil is `6/h⁵` | **yes** |
//! | `full` | nothing — the true field on `33³`, C2's baseline | no |
//! | `sparse` | the Smolyak level `L ∈ {5, 6, 7, 8}` | no |
//! | `extraction_floor` | nothing — the true field on the sparse arm's `65³` grid | **yes** |
//!
//! # C1: which norm, measured where, and why the band shrinks
//!
//! `∂⁶u/∂x²∂y²∂z²` is taken as the tensor cube of the three-point central second
//! difference: the `3×3×3` stencil `(1, −2, 1)⊗(1, −2, 1)⊗(1, −2, 1)` divided by
//! `hx²·hy²·hz²`. Its generating polynomial is `((z−1)²/z)³ = (z−1)⁶/z³`, so it
//! annihilates every polynomial of degree ≤ 5 in the stencil offset — which is
//! why the two calibration arms have *exact* closed forms and are asserted
//! against them rather than eyeballed.
//!
//! The norm is a maximum over the **surface band** `|u(p)| ≤ 2h`, not over the
//! whole domain, and both are recorded. The band is the right region for C1
//! because C1's own words scope the claim to the surface — *"unbounded on the CSG
//! fields … **where the surface has sharp edges**"* — and an SDF is singular on
//! its medial axis whether or not its surface is smooth. `sphere`'s medial axis
//! is the single point `r = 0`, two units from its surface; `box_exact`'s
//! singular set *is* its edge and corner set, which lies **on** the surface.
//! Measuring in a band that shrinks with `h` is therefore the discriminator: a
//! positive growth exponent means the singularity survives every shrinkage, i.e.
//! it is on the surface.
//!
//! The whole-domain norm is recorded beside it as `norm_domain_finest` /
//! `norm_domain_exponent` because *that* is the quantity Bungartz & Griebel's
//! theorem needs, over the box the interpolant is built on. It is the column that
//! explains C2, and the two columns are expected to disagree.
//!
//! Four resolutions, all odd, `33/49/65/97`. Odd matters: over the `[-2, 2]³`
//! domain `3/h = 3(n−1)/4` is `24, 36, 48, 72` and `2/h = (n−1)/2` is
//! `16, 24, 32, 48`, so the box corners at `±1` **and** the domain centre are
//! exact grid points at every rung. The alignment of the singular set is then
//! constant across the ladder and the fitted exponent is a property of the field
//! rather than of the grid phase — which is `M-266`'s lesson applied to a
//! derivative instead of a triangle count.
//!
//! ## The roundoff floor, and why a regime column exists
//!
//! A sixth difference divides by `h⁶`, so it amplifies the `O(ε·|u|)` rounding
//! already in the sampled values by `1/h⁶`. The stencil's `ℓ1` weight is
//! `(1+2+1)³ = 64`, giving a floor of `64·ε·max|u| / (hx²hy²hz²)` below which the
//! instrument measures its own arithmetic. For a field whose exact mixed
//! derivative is zero the measured value sits **at** that floor and its apparent
//! growth is the floor's own `h⁻⁶` — a divergence that says nothing about the
//! field. So `norm_regime` reports `resolved` or `at_roundoff_floor` from
//! `norm_over_floor_min`, and an unresolved field is recorded as bounded, which
//! is the truth: its sixth mixed derivative is numerically indistinguishable
//! from zero. `fbm_terrain` is the field this exists for, and the reason is
//! structural (affine in `y`), not numerical.
//!
//! # C2: the Smolyak combination, and what it is charged for
//!
//! The reconstruction is the **Smolyak combination of tensor-product multilinear
//! interpolants**, built bench-locally and stated here rather than assumed:
//!
//! ```text
//! A_L f = Σ_{|i| ∈ [L−2, L]} (−1)^(L−|i|) · C(2, L−|i|) · T_i f
//! ```
//!
//! with `i ∈ ℕ³`, `T_i` the trilinear interpolant on the nested one-dimensional
//! grids of `2^{i_k} + 1 = m_{i_k}` points per axis, and `d = 3` so the
//! combination coefficients are `(+1, −2, +1)`. Their numerical agreement with
//! the difference stencil above is a coincidence of `d = 3` and the two constants
//! are kept separate on purpose.
//!
//! Nestedness makes the sample set exact rather than estimated: `grid_i ⊂
//! grid_{i'}` whenever `i ≤ i'` componentwise, so `∪_{|i| ≤ L} grid_i =
//! ∪_{|i| = L} grid_i` and every point any term of the combination reads lies in
//! that union. `sparse_grid_points` is therefore the size of one shared cache
//! keyed by the finest-grid integer triple — the count of **distinct field
//! evaluations**, not a formula. It comes out `705, 1649, 3809, 8705` for
//! `L = 5, 6, 7, 8`, independent of the field.
//!
//! - **Full arm:** Marching Cubes on the true field at `33³`, the resolution a
//!   chunk in this repository is actually meshed at. `full_grid_points = 35 937`.
//! - **Sparse arm:** the headline level is `L = 8`, the largest level whose point
//!   count still fits C2's *"at least 2× fewer samples"*: `8705 ≤ 35 937/2`.
//!   `point_ratio = 4.128`, so C2 is offered twice the budget advantage it asks
//!   for and a negative verdict cannot be blamed on the budget. That arithmetic
//!   is a vacuity control, not a comment.
//! - The sparse arm is extracted on **`65³`**, twice as fine as the arm it is
//!   compared against, because its cost is counted in field evaluations and the
//!   reconstruction is defined everywhere. `extraction_floor` measures the *true*
//!   field on that same `65³` grid, which pins empirically how much of
//!   `hausdorff_sparse` the extraction grid could possibly be responsible for. A
//!   generous instrument makes a negative result stronger.
//!
//! Both arms are scored by the same instrument, `validate::accuracy` against the
//! **true** field, and `symmetric_hausdorff` is `max(mesh→field, field→mesh)`
//! (`validate/accuracy.rs:236`). A reconstruction that loses the surface
//! altogether produces zero samples in both directions, at which point
//! `symmetric_hausdorff` reads `0.0` and a lost surface would look like a perfect
//! one — the trap `experiment_p100.rs:1559-1570` caught by hand. Here it is
//! closed by definition: no coverage means no measurement, and an unmeasured
//! distance is recorded as `inf`, never as zero. `sparse_coverage` carries the
//! predicate.
//!
//! `ablation.rs:200` and `shootout.rs:183` gate a Hausdorff on
//! `field.bound().is_exact()`, and three fields in the roster are not exact
//! (`gyroid` is `Lipschitz`, `csg_difference` is `Underestimate`, `fbm_terrain`
//! and `noise_cavity` are `Unbounded`). This row does **not** skip them: `bound_exact`
//! is a recorded column, both arms of a row use the identical instrument on the
//! identical grid, and C2 reads a *comparison* between two arms rather than a
//! certified distance. Where `bound_exact` is false the number is a residual
//! rather than a distance, and the comparison is still the comparison C2 asks
//! for.
//!
//! ## When `hausdorff_sparse < hausdorff_full`, and why it is not evidence
//!
//! Three mechanisms can invert the two columns on a **non-qualifying** row, and
//! `hausdorff_sparse_floor` is the column that separates them from a real win.
//!
//! 1. The sparse arm extracts on `65³` and the arm it is compared with on `33³`.
//!    On a field whose error is dominated by extraction-grid resolution rather
//!    than by reconstruction error, the finer grid simply wins.
//! 2. `box_exact` and `ThinPlate` are piecewise linear on their faces, so the
//!    *multilinear* interpolant is exact there and the reconstruction throws
//!    nothing away over most of the surface.
//! 3. `accuracy`'s mesh→field direction is an **over**-estimate near a sharp or
//!    concave seam, because the Newton flow can land further away than the true
//!    nearest point (`validate/accuracy.rs:59-63`). Two meshes that place a box
//!    corner differently are not scored on the same footing there.
//!
//! The discriminator is `hausdorff_sparse < hausdorff_sparse_floor`: the floor
//! arm is the **true** field extracted on the sparse arm's own grid, so a
//! reconstruction beating it is beating the field it approximates. That is a
//! statement about the instrument at a sharp feature, not about Smolyak. C1
//! excludes those rows from C2's scope anyway, which is what SHARE means by
//! *"on qualifying fields only"* — but the columns are recorded on all eight so
//! that the exclusion can be checked rather than trusted.
//!
//! # SHARE, recomputed before the numbers
//!
//! SHARE: *"C2 moves the field-evaluation stage on qualifying fields only, and
//! the qualifying set is C1's output."* Discharged, and negatively, by arithmetic
//! available before any Hausdorff is looked at.
//!
//! The sparse arm saves `35 937 − 8705 = 27 232` field evaluations. To spend
//! them it evaluates a 109-term combination at every one of the `65³ = 274 625`
//! points of its extraction grid, and each term is one trilinear interpolation —
//! eight array reads. That is `8 · 109 · 274 625 = 239 473 000` reads, so a single
//! `Sdf::sample` would have to cost more than `239 473 000 / 27 232 ≈ 8794`
//! memory reads before the trade is even neutral. The most expensive field in the
//! roster is a four-octave lattice-hash fbm evaluating no transcendental at all
//! (`fields/mod.rs:1294-1295`). The mechanism cannot move the field-evaluation
//! stage on this crate's fields at any accuracy, and `break_even_reads_per_eval`
//! records the number. Sparse grids are for fields whose evaluation is a PDE
//! solve; an SDF is not one.
//!
//! # Verdicts: which are global and which are per row
//!
//! **`c1_holds` is global** and carries the same value on every row, because C1
//! is a statement about the roster: the norm is measured on all eight fields and
//! is unbounded on `box_exact` **and** `csg_difference`. **`c2_holds` is per
//! row** — the row's own field is in C1's qualifying set and clears both halves
//! of C2 — and the conjunction over the qualifying set is recorded once as the
//! extra `c2_global`. `c2_reachable` records, per row, whether
//! `hausdorff_full > 0`: a zero baseline makes *"matched Hausdorff error"*
//! arithmetically unreachable on that field, and an unreachable clause is
//! recorded with its arithmetic rather than dropped.
//!
//! No clause here is timed. `field_wall_seconds` is recorded because a reader
//! will ask, and is read by nothing — `M-280`'s rule, and `P-126`'s precedent.
//!
//! # Vacuity controls
//!
//! - **`stencil_smooth`** must reproduce the exact value `8` on `x²y²z²` at every
//!   resolution, to `1e-3` relative. Column: `control_smooth_value`. The stencil
//!   factorises on a separable input, so `8` is a closed form and not a fit; a
//!   miss means every norm below is measuring the instrument.
//! - **`stencil_kink`** must reproduce the exact constant `6` on
//!   `max(x+y+z, 0)`, whose stencil sum over the offsets is
//!   `Σ_s (−1)^{3−s} C(6, s+3)·max(s, 0) = 15 − 12 + 3 = 6`, giving `6/h⁵` at any
//!   grid point on the kink plane. Columns: `control_kink_constant` (asserted to
//!   `1e-6` relative) and `control_kink_exponent` (asserted to be `5` within
//!   `0.01`, which is also the check on the least-squares fit the field verdicts
//!   depend on). Together the two arms prove the instrument returns the right
//!   number for a smooth function **and** detects a Lipschitz kink at the right
//!   rate — the dynamic range C1 needs in order to say anything.
//! - **`sphere` must show a bounded norm** — the registered control. Column:
//!   `norm_finite` on the `sphere` row. An unbounded `sphere` means the
//!   measurement is broken rather than the fields being rough.
//! - **The band must be non-empty** on every field at every resolution. Column:
//!   `band_points_finest`. A maximum over an empty set is not a norm (`M-44`).
//! - **The full arm must extract something and must have accuracy coverage** on
//!   every field. Columns: `full_triangles`, `hausdorff_full`. Without a baseline
//!   there is nothing for C2 to match.
//! - **C1's qualifying set must be non-empty**, or C2 is scoped to the empty set
//!   and its verdict is vacuous whichever way it comes out. Column:
//!   `fields_qualifying`.
//! - **The headline sparse level must offer at least the `2×` C2 registers.**
//!   Column: `point_ratio`. If it did not, a negative C2 would be a statement
//!   about this harness's budget rather than about the method.

#![allow(
    clippy::cast_possible_truncation,
    clippy::cast_precision_loss,
    clippy::too_many_lines
)]

mod common;

use std::collections::BTreeMap;
use std::time::Instant;

use isomesh::fields::ReferenceField;
use isomesh::marching_cubes::MarchingCubes;
use isomesh::validate::{AccuracyConfig, accuracy};
use isomesh::{MeshBuffer, RuntimeShape3, Sdf};

// ─── constants ──────────────────────────────────────────────────────────────

/// The eight names `for_each_reference_field!` walks.
const FIELD_COUNT: usize = 8;

/// C1's resolution ladder, ascending. All odd, so the box corners at `±1` and
/// the domain centre are exact grid points at every rung.
const C1_RESOLUTIONS: [u32; 4] = [33, 49, 65, 97];

/// Band half-width for the surface norm, in cells: `|u(p)| ≤ BAND_CELLS·h`.
const BAND_CELLS: f64 = 2.0;

/// The three-point central second difference, before the `h²`.
const STENCIL: [f64; 3] = [1.0, -2.0, 1.0];

/// `ℓ1` norm of the tensor cube of [`STENCIL`]: `(1 + 2 + 1)³`.
const STENCIL_L1: f64 = 64.0;

/// A norm must clear this multiple of its own roundoff floor to be a
/// measurement of a derivative rather than of the arithmetic.
const FLOOR_MARGIN: f64 = 4.0;

/// Growth exponents below this are converging, at or above are diverging. The
/// weakest genuine divergence a sixth difference can show is `h⁻¹` (a bounded
/// jump in the fifth derivative), so a midpoint between `0` and `1` separates
/// the regimes without needing to be tuned.
const DIVERGENCE_EXPONENT: f64 = 0.5;

/// The exact `∂⁶/∂x²∂y²∂z²` of the smooth calibration field `x²y²z²`.
const CONTROL_SMOOTH_EXACT: f64 = 8.0;

/// Relative tolerance on [`CONTROL_SMOOTH_EXACT`]. The floor at the finest rung
/// is `2.2e-5` of `8`, so this is roughly fifty times the arithmetic.
const CONTROL_SMOOTH_TOLERANCE: f64 = 1e-3;

/// The exact stencil sum of `max(x+y+z, 0)` at a grid point on the kink plane:
/// `Σ_s (−1)^{3−s} C(6, s+3)·max(s, 0) = 15·1 − 6·2 + 1·3`.
const CONTROL_KINK_CONSTANT: f64 = 6.0;

/// Relative tolerance on [`CONTROL_KINK_CONSTANT`].
const CONTROL_KINK_TOLERANCE: f64 = 1e-6;

/// The exact growth exponent of a Lipschitz kink under a sixth difference: the
/// function is homogeneous of degree one about the kink, so the stencil sum is
/// `O(h)` and the quotient is `O(h⁻⁵)`.
const CONTROL_KINK_EXPONENT: f64 = 5.0;

/// Tolerance on [`CONTROL_KINK_EXPONENT`], which is also the check on the
/// least-squares fit the field verdicts are read from.
const CONTROL_EXPONENT_TOLERANCE: f64 = 0.01;

/// C2's baseline: samples per axis of the full grid, the resolution a chunk in
/// this repository is meshed at.
const FULL_SAMPLES: u32 = 33;

/// Samples per axis the sparse reconstruction is extracted on. Twice the full
/// arm's, deliberately; `extraction_floor` measures what that buys.
const SPARSE_EXTRACTION_SAMPLES: u32 = 65;

/// Smolyak levels, ascending. `L − 2 ≥ 0` is required by the combination.
const SPARSE_LEVELS: [u32; 4] = [5, 6, 7, 8];

/// The level C2's registered columns are read at: the largest in
/// [`SPARSE_LEVELS`] whose point count still clears [`POINT_RATIO_FLOOR`].
const HEADLINE_LEVEL: u32 = 8;

/// C2's registered sample advantage.
const POINT_RATIO_FLOOR: f64 = 2.0;

/// Smolyak combination coefficients `(−1)^k·C(d−1, k)` for `d = 3`, `k = L−|i|`.
///
/// Numerically equal to [`STENCIL`] and for an unrelated reason — one is
/// `(z−1)²`, the other is `(1−z)^{d−1}` at `d = 3`. Kept separate so that
/// changing `d` cannot silently change the difference operator.
const COMBINATION: [f64; 3] = [1.0, -2.0, 1.0];

/// The convergence rate in `L` that `H²_mix` regularity would buy: the sparse
/// interpolant's error falls like `2^{−2L}` up to logarithms.
const SMOLYAK_THEORY_RATE: f64 = -2.0;

// ─── the sixth mixed difference ─────────────────────────────────────────────

/// One resolution's reading of `‖∂⁶u/∂x²∂y²∂z²‖`.
#[derive(Clone, Copy, Debug)]
struct Rung {
    /// Samples per axis.
    samples: u32,
    /// Grid spacing, geometric mean over the three axes.
    h: f64,
    /// Maximum `|D|` over the surface band `|u| ≤ BAND_CELLS·h`.
    band_norm: f64,
    /// Maximum `|D|` over every interior grid point.
    domain_norm: f64,
    /// Band population, so the maximum is known not to be over nothing.
    band_points: u64,
    /// `STENCIL_L1·ε·max|u| / (hx²hy²hz²)`: the smallest value this instrument
    /// can distinguish from its own rounding at this spacing.
    floor: f64,
}

/// Sample `field` on an `n³` grid and take the sixth mixed difference on it.
///
/// The stencil reads the value grid rather than the field, so the `3×3×3`
/// neighbourhood is the grid's own neighbourhood and the spacing in the
/// denominator is exactly the spacing that was sampled. One field evaluation per
/// grid point and no more.
fn mixed_rung<S: Sdf<Scalar = f64>>(field: &S, lo: [f64; 3], hi: [f64; 3], samples: u32) -> Rung {
    let n = samples as usize;
    let steps = f64::from(samples - 1);
    let h = [
        (hi[0] - lo[0]) / steps,
        (hi[1] - lo[1]) / steps,
        (hi[2] - lo[2]) / steps,
    ];

    let mut values = vec![0.0f64; n * n * n];
    let mut u_scale = 0.0f64;
    for k in 0..n {
        for j in 0..n {
            for i in 0..n {
                let p = [
                    lo[0] + i as f64 * h[0],
                    lo[1] + j as f64 * h[1],
                    lo[2] + k as f64 * h[2],
                ];
                let v = field.sample(p);
                values[i + n * (j + n * k)] = v;
                u_scale = u_scale.max(v.abs());
            }
        }
    }
    assert!(
        u_scale > 0.0,
        "the field is identically zero on its own {samples}³ grid, so there is \
         no scale to state a roundoff floor against"
    );

    let denom = h[0] * h[0] * h[1] * h[1] * h[2] * h[2];
    let spacing = (h[0] * h[1] * h[2]).cbrt();
    let band = BAND_CELLS * spacing;

    let mut band_norm = 0.0f64;
    let mut domain_norm = 0.0f64;
    let mut band_points = 0u64;
    for k in 1..n - 1 {
        for j in 1..n - 1 {
            for i in 1..n - 1 {
                let mut acc = 0.0f64;
                for (dk, wk) in STENCIL.iter().enumerate() {
                    for (dj, wj) in STENCIL.iter().enumerate() {
                        let row = i + n * ((j + dj - 1) + n * (k + dk - 1));
                        let along_x = STENCIL[0] * values[row - 1]
                            + STENCIL[1] * values[row]
                            + STENCIL[2] * values[row + 1];
                        acc += wk * wj * along_x;
                    }
                }
                let d = (acc / denom).abs();
                domain_norm = domain_norm.max(d);
                if values[i + n * (j + n * k)].abs() <= band {
                    band_norm = band_norm.max(d);
                    band_points += 1;
                }
            }
        }
    }

    Rung {
        samples,
        h: spacing,
        band_norm,
        domain_norm,
        band_points,
        floor: STENCIL_L1 * f64::EPSILON * u_scale / denom,
    }
}

/// Least-squares slope of `log2(v)` against `log2(1/h)` over the ladder.
///
/// Positive means the quantity grows as the grid refines — a divergence. A
/// reading below its own [`Rung::floor`] is fitted **as** the floor, because the
/// instrument cannot resolve anything smaller and `log2(0)` is not a number.
fn growth_exponent(rungs: &[Rung], pick: impl Fn(&Rung) -> f64) -> f64 {
    let count = rungs.len() as f64;
    let xs: Vec<f64> = rungs.iter().map(|r| -r.h.log2()).collect();
    let ys: Vec<f64> = rungs.iter().map(|r| pick(r).max(r.floor).log2()).collect();
    let mx = xs.iter().sum::<f64>() / count;
    let my = ys.iter().sum::<f64>() / count;
    let mut sxy = 0.0f64;
    let mut sxx = 0.0f64;
    for (&x, &y) in xs.iter().zip(&ys) {
        sxy += (x - mx) * (y - my);
        sxx += (x - mx) * (x - mx);
    }
    if sxx <= 0.0 { 0.0 } else { sxy / sxx }
}

/// C1's verdict for one field, and the evidence behind it.
#[derive(Clone, Copy, Debug)]
struct NormVerdict {
    /// Band norm at the finest rung: the registered `mixed_derivative_norm`.
    norm: f64,
    /// Growth exponent of the band norm.
    exponent: f64,
    /// Whole-domain norm at the finest rung — the quantity the sparse-grid
    /// theorem needs, over the box the interpolant is built on.
    domain_norm: f64,
    /// Growth exponent of the whole-domain norm.
    domain_exponent: f64,
    /// Worst signal-to-floor ratio over the ladder.
    over_floor_min: f64,
    /// Every rung cleared [`FLOOR_MARGIN`] times its floor.
    resolved: bool,
    /// The registered `norm_finite`.
    finite: bool,
}

impl NormVerdict {
    /// `resolved` / `at_roundoff_floor`, the regime the verdict was read in.
    fn regime(self) -> &'static str {
        if self.resolved {
            "resolved"
        } else {
            "at_roundoff_floor"
        }
    }
}

/// Classify one field's ladder.
///
/// An unresolved ladder is bounded and not divergent: below the roundoff floor
/// the measured growth is the floor's own `h⁻⁶`, so the field's exact sixth
/// mixed derivative is numerically indistinguishable from zero.
fn norm_verdict(rungs: &[Rung]) -> NormVerdict {
    let finest = rungs.last().expect("the resolution ladder is not empty");
    let over_floor_min = rungs
        .iter()
        .map(|r| r.band_norm / r.floor)
        .fold(f64::INFINITY, f64::min);
    let resolved = over_floor_min > FLOOR_MARGIN;
    let exponent = growth_exponent(rungs, |r| r.band_norm);
    NormVerdict {
        norm: finest.band_norm,
        exponent,
        domain_norm: finest.domain_norm,
        domain_exponent: growth_exponent(rungs, |r| r.domain_norm),
        over_floor_min,
        resolved,
        finite: !resolved || exponent < DIVERGENCE_EXPONENT,
    }
}

// ─── the two stencil calibration fields ─────────────────────────────────────

/// `x²y²z²`, whose exact `∂⁶/∂x²∂y²∂z²` is `2·2·2 = 8` everywhere.
///
/// Separable, so the stencil factorises into three exact second differences of
/// a quadratic — `((x+h)² − 2x² + (x−h)²)/h² = 2` — and the closed form is not
/// an asymptotic statement about small `h` but an identity at every `h`.
#[derive(Clone, Copy, Debug)]
struct SmoothControl;

impl Sdf for SmoothControl {
    type Scalar = f64;
    fn sample(&self, p: [f64; 3]) -> f64 {
        p[0] * p[0] * p[1] * p[1] * p[2] * p[2]
    }
}

/// `max(x+y+z, 0)`: a Lipschitz kink across a plane no axis is aligned with.
///
/// Homogeneous of degree one about the kink, so the stencil sum at a grid point
/// on the plane is `h·Σ_s a_s·max(s, 0)` with `a_s = (−1)^{3−s}C(6, s+3)`, i.e.
/// `6h`, and the quotient is exactly `6/h⁵`.
#[derive(Clone, Copy, Debug)]
struct KinkControl;

impl Sdf for KinkControl {
    type Scalar = f64;
    fn sample(&self, p: [f64; 3]) -> f64 {
        (p[0] + p[1] + p[2]).max(0.0)
    }
}

// ─── the Smolyak sparse grid ────────────────────────────────────────────────

/// One tensor-product multilinear interpolant of the combination.
#[derive(Clone, Debug)]
struct Term {
    /// The combination coefficient `(−1)^{L−|i|}·C(2, L−|i|)`.
    coeff: f64,
    /// Cells per axis, `2^{i_k}`.
    cells: [usize; 3],
    /// Strides into [`Term::values`], `x` fastest.
    stride: [usize; 3],
    /// The field at this term's own grid nodes.
    values: Vec<f64>,
}

/// A level-`L` Smolyak reconstruction of a field over a box.
#[derive(Clone, Debug)]
struct Smolyak {
    /// Domain minimum.
    lo: [f64; 3],
    /// Reciprocal of the per-axis span, so `sample` divides nothing.
    inv_span: [f64; 3],
    /// The combination.
    terms: Vec<Term>,
    /// Distinct field evaluations the whole reconstruction cost.
    points: usize,
}

impl Smolyak {
    /// Build the level-`level` combination, evaluating the field once per
    /// distinct sparse-grid point.
    ///
    /// Node coordinates are computed from the **finest**-grid integer index, so
    /// two terms that share a point compute bit-identical coordinates and the
    /// shared cache returns one value rather than two roundings of it.
    fn build<S: Sdf<Scalar = f64>>(field: &S, lo: [f64; 3], hi: [f64; 3], level: u32) -> Self {
        assert!(
            level >= 2,
            "the d = 3 combination reads levels L−2 … L, so L < 2 is not a grid"
        );
        let fine = 1u32 << level;
        let fine_f = f64::from(fine);
        let span = [hi[0] - lo[0], hi[1] - lo[1], hi[2] - lo[2]];

        let mut cache: BTreeMap<[u32; 3], f64> = BTreeMap::new();
        let mut terms: Vec<Term> = Vec::new();
        for s in (level - 2)..=level {
            let coeff = COMBINATION[(level - s) as usize];
            for i0 in 0..=s {
                for i1 in 0..=(s - i0) {
                    let i2 = s - i0 - i1;
                    let cells = [1usize << i0, 1usize << i1, 1usize << i2];
                    let step = [fine >> i0, fine >> i1, fine >> i2];
                    let stride = [1, cells[0] + 1, (cells[0] + 1) * (cells[1] + 1)];
                    let mut values = vec![0.0f64; (cells[0] + 1) * (cells[1] + 1) * (cells[2] + 1)];
                    for c2 in 0..=cells[2] {
                        for c1 in 0..=cells[1] {
                            for c0 in 0..=cells[0] {
                                let key = [
                                    c0 as u32 * step[0],
                                    c1 as u32 * step[1],
                                    c2 as u32 * step[2],
                                ];
                                let v = *cache.entry(key).or_insert_with(|| {
                                    field.sample([
                                        lo[0] + span[0] * f64::from(key[0]) / fine_f,
                                        lo[1] + span[1] * f64::from(key[1]) / fine_f,
                                        lo[2] + span[2] * f64::from(key[2]) / fine_f,
                                    ])
                                });
                                values[c0 * stride[0] + c1 * stride[1] + c2 * stride[2]] = v;
                            }
                        }
                    }
                    terms.push(Term {
                        coeff,
                        cells,
                        stride,
                        values,
                    });
                }
            }
        }

        Self {
            lo,
            inv_span: [1.0 / span[0], 1.0 / span[1], 1.0 / span[2]],
            points: cache.len(),
            terms,
        }
    }

    /// Array reads one [`Sdf::sample`] costs: eight per term.
    fn reads_per_sample(&self) -> u64 {
        8 * self.terms.len() as u64
    }
}

impl Sdf for Smolyak {
    type Scalar = f64;

    fn sample(&self, p: [f64; 3]) -> f64 {
        let u = [
            (p[0] - self.lo[0]) * self.inv_span[0],
            (p[1] - self.lo[1]) * self.inv_span[1],
            (p[2] - self.lo[2]) * self.inv_span[2],
        ];
        let mut acc = 0.0f64;
        for t in &self.terms {
            let mut base = 0usize;
            let mut f = [0.0f64; 3];
            for a in 0..3 {
                let scaled = u[a] * t.cells[a] as f64;
                let floored = scaled.floor();
                // Clamped, not wrapped: the extraction grid touches both
                // boundaries of the domain and the upper one lands exactly on
                // `cells`, which belongs to the last cell with `f = 1`.
                let ci = if floored <= 0.0 {
                    0
                } else if floored >= t.cells[a] as f64 {
                    t.cells[a] - 1
                } else {
                    floored as usize
                };
                f[a] = scaled - ci as f64;
                base += ci * t.stride[a];
            }
            let (sx, sy, sz) = (t.stride[0], t.stride[1], t.stride[2]);
            let v = &t.values;
            let x00 = v[base] + f[0] * (v[base + sx] - v[base]);
            let x10 = v[base + sy] + f[0] * (v[base + sy + sx] - v[base + sy]);
            let x01 = v[base + sz] + f[0] * (v[base + sz + sx] - v[base + sz]);
            let x11 = v[base + sz + sy] + f[0] * (v[base + sz + sy + sx] - v[base + sz + sy]);
            let y0 = x00 + f[1] * (x10 - x00);
            let y1 = x01 + f[1] * (x11 - x01);
            acc += t.coeff * (y0 + f[2] * (y1 - y0));
        }
        acc
    }
}

// ─── one extraction arm ─────────────────────────────────────────────────────

/// What one Marching Cubes extraction, scored against the true field, produced.
#[derive(Clone, Copy, Debug)]
struct Arm {
    /// `max(mesh→field, field→mesh)`, or `inf` when there was nothing to score.
    hausdorff: f64,
    /// Both directions produced samples.
    coverage: bool,
    /// Triangles emitted.
    triangles: u64,
    /// Mesh samples whose projection onto the true surface did not converge and
    /// are therefore **excluded** from `mesh→field`.
    unconverged: u64,
}

/// Extract `surface` on the given grid, measure the result against `truth`.
///
/// The two are the same field for the full arm and the extraction-floor control,
/// and differ only for the sparse arm — where `surface` is the reconstruction
/// and `truth` is what the reconstruction is being judged against. Measuring the
/// reconstruction against itself would compare it with itself and could only
/// ever look good, which is `experiment_p46.rs:129-132`'s rule.
fn arm<S, T>(surface: &S, truth: &T, shape: &RuntimeShape3, origin: [f64; 3], cell_size: f64) -> Arm
where
    S: Sdf<Scalar = f64>,
    T: Sdf<Scalar = f64>,
{
    let mut mesh = MeshBuffer::<f64>::new();
    let mut mc = MarchingCubes::<f64>::new();
    mc.extract(surface, shape, origin, cell_size, &mut mesh)
        .expect("marching cubes on a benchmark grid");
    let cfg = AccuracyConfig::from_cell_size(cell_size).expect("a positive finite cell size");
    let report = accuracy(&mesh.positions, &mesh.indices, truth, shape, origin, &cfg)
        .expect("accuracy on the grid the mesh was extracted on");
    Arm {
        // An unmeasured distance is not a zero. With no coverage
        // `symmetric_hausdorff` reads `0.0`, so a reconstruction that lost the
        // surface entirely would score better than one that found it.
        hausdorff: if report.has_coverage() {
            report.symmetric_hausdorff()
        } else {
            f64::INFINITY
        },
        coverage: report.has_coverage(),
        triangles: mesh.triangle_count() as u64,
        unconverged: report.unconverged_mesh_samples,
    }
}

// ─── per-field measurement ──────────────────────────────────────────────────

/// One rung of the sparse ladder.
#[derive(Clone, Copy, Debug)]
struct SparseRung {
    /// Smolyak level.
    level: u32,
    /// Distinct field evaluations the reconstruction cost.
    points: usize,
    /// Terms in the combination.
    terms: usize,
    /// Array reads one reconstruction sample costs.
    reads_per_sample: u64,
    /// What the reconstruction extracted to.
    arm: Arm,
}

/// Everything one field contributes to the CSV.
#[derive(Clone, Debug)]
struct Measured {
    /// The `for_each_reference_field!` name.
    name: &'static str,
    /// `field.bound().is_exact()`: whether the Hausdorff is a certified distance
    /// or a residual under the same instrument.
    bound_exact: bool,
    /// C1's ladder.
    rungs: Vec<Rung>,
    /// C1's verdict.
    c1: NormVerdict,
    /// C2's baseline.
    full: Arm,
    /// The true field on the sparse arm's own grid: how much of
    /// `hausdorff_sparse` the extraction grid could be responsible for.
    extraction_floor: Arm,
    /// C2's ladder, ascending in level.
    ladder: Vec<SparseRung>,
    /// Wall clock for this field. Read by no clause.
    wall_seconds: f64,
}

impl Measured {
    /// The rung C2's registered columns are read at.
    fn headline(&self) -> &SparseRung {
        self.ladder
            .iter()
            .find(|r| r.level == HEADLINE_LEVEL)
            .expect("the headline level is one of SPARSE_LEVELS")
    }

    /// Field evaluations the full arm spends.
    fn full_points(&self) -> u64 {
        u64::from(FULL_SAMPLES) * u64::from(FULL_SAMPLES) * u64::from(FULL_SAMPLES)
    }

    /// `full_grid_points / sparse_grid_points`: the sample advantage C2 is
    /// offered.
    fn point_ratio(&self) -> f64 {
        self.full_points() as f64 / self.headline().points as f64
    }

    /// Whether *"matched Hausdorff error"* is arithmetically reachable at all: a
    /// zero baseline cannot be matched from above.
    fn c2_reachable(&self) -> bool {
        self.full.hausdorff > 0.0 && self.full.hausdorff.is_finite()
    }

    /// Reads the reconstruction spends over one extraction pass, against the
    /// field evaluations it saved: SHARE's break-even, per field.
    fn break_even_reads_per_eval(&self) -> f64 {
        let head = self.headline();
        let grid = u64::from(SPARSE_EXTRACTION_SAMPLES).pow(3);
        let saved = self.full_points().saturating_sub(head.points as u64);
        if saved == 0 {
            f64::INFINITY
        } else {
            (head.reads_per_sample * grid) as f64 / saved as f64
        }
    }
}

/// Measure one reference field: C1's ladder, then C2's two arms.
fn measure<F>(name: &'static str, field: &F) -> Measured
where
    F: ReferenceField + Sdf<Scalar = f64>,
{
    let started = Instant::now();
    let (lo, hi) = field.domain();

    let rungs: Vec<Rung> = C1_RESOLUTIONS
        .iter()
        .map(|&n| mixed_rung(field, lo, hi, n))
        .collect();
    let c1 = norm_verdict(&rungs);

    // C2's baseline, on the house grid for this field at `FULL_SAMPLES`.
    let (full_shape, full_origin, full_cell) = common::grid::<f64, _>(field, FULL_SAMPLES);
    let full = arm(field, field, &full_shape, full_origin, full_cell);

    // The sparse arm's grid, by the same convention `common::grid` uses: the
    // spacing is read off axis 0 and every reference field's domain is a cube.
    let sparse_shape =
        RuntimeShape3::new([SPARSE_EXTRACTION_SAMPLES; 3]).expect("65³ fits u32 indices");
    let sparse_cell = (hi[0] - lo[0]) / f64::from(SPARSE_EXTRACTION_SAMPLES - 1);
    let extraction_floor = arm(field, field, &sparse_shape, lo, sparse_cell);

    let mut ladder: Vec<SparseRung> = Vec::new();
    for &level in &SPARSE_LEVELS {
        let smolyak = Smolyak::build(field, lo, hi, level);
        let measured = arm(&smolyak, field, &sparse_shape, lo, sparse_cell);
        ladder.push(SparseRung {
            level,
            points: smolyak.points,
            terms: smolyak.terms.len(),
            reads_per_sample: smolyak.reads_per_sample(),
            arm: measured,
        });
    }

    Measured {
        name,
        bound_exact: field.bound().is_exact(),
        rungs,
        c1,
        full,
        extraction_floor,
        ladder,
        wall_seconds: started.elapsed().as_secs_f64(),
    }
}

/// Least-squares slope of `log2(hausdorff)` against the Smolyak level.
///
/// `H²_mix` regularity would put this at [`SMOLYAK_THEORY_RATE`]. `NaN` when
/// fewer than two rungs produced a finite positive distance, because a rate
/// through one point is not a rate.
fn ladder_rate(ladder: &[SparseRung]) -> f64 {
    let pts: Vec<(f64, f64)> = ladder
        .iter()
        .filter(|r| r.arm.hausdorff.is_finite() && r.arm.hausdorff > 0.0)
        .map(|r| (f64::from(r.level), r.arm.hausdorff.log2()))
        .collect();
    if pts.len() < 2 {
        return f64::NAN;
    }
    let count = pts.len() as f64;
    let mx = pts.iter().map(|p| p.0).sum::<f64>() / count;
    let my = pts.iter().map(|p| p.1).sum::<f64>() / count;
    let mut sxy = 0.0f64;
    let mut sxx = 0.0f64;
    for &(x, y) in &pts {
        sxy += (x - mx) * (y - my);
        sxx += (x - mx) * (x - mx);
    }
    if sxx <= 0.0 { 0.0 } else { sxy / sxx }
}

/// `n33=4.872e1|n49=…`: the C1 ladder as one CSV-safe token.
fn norm_ladder_token(rungs: &[Rung]) -> String {
    let mut out = String::new();
    for r in rungs {
        if !out.is_empty() {
            out.push('|');
        }
        out.push_str(&format!("n{}={:.3e}", r.samples, r.band_norm));
    }
    out
}

/// `L5=1.392e-1|L6=…`: one field's sparse Hausdorff ladder.
fn hausdorff_ladder_token(ladder: &[SparseRung]) -> String {
    let mut out = String::new();
    for r in ladder {
        if !out.is_empty() {
            out.push('|');
        }
        out.push_str(&format!("L{}={:.3e}", r.level, r.arm.hausdorff));
    }
    out
}

/// `L5=705|L6=1649|…`: the sparse point counts, which are geometry and identical
/// across fields.
fn points_ladder_token(ladder: &[SparseRung]) -> String {
    let mut out = String::new();
    for r in ladder {
        if !out.is_empty() {
            out.push('|');
        }
        out.push_str(&format!("L{}={}", r.level, r.points));
    }
    out
}

fn main() {
    if !std::env::args().any(|arg| arg == "--bench") {
        return;
    }
    let prereg = isomesh::experiment!("P-158");

    common::experiment::run(prereg, |run| {
        // ── the two calibration arms, which decide whether anything below is a
        //    measurement of a field rather than of the instrument ─────────────
        let smooth: Vec<Rung> = C1_RESOLUTIONS
            .iter()
            .map(|&n| mixed_rung(&SmoothControl, [-2.0; 3], [2.0; 3], n))
            .collect();
        let kink: Vec<Rung> = C1_RESOLUTIONS
            .iter()
            .map(|&n| mixed_rung(&KinkControl, [-2.0; 3], [2.0; 3], n))
            .collect();

        let mut worst_smooth = 0.0f64;
        for r in &smooth {
            let relative = (r.domain_norm - CONTROL_SMOOTH_EXACT).abs() / CONTROL_SMOOTH_EXACT;
            worst_smooth = worst_smooth.max(relative);
            assert!(
                relative <= CONTROL_SMOOTH_TOLERANCE,
                "VOID: the sixth mixed difference reads {:.12e} on x^2y^2z^2 at \
                 {}^3, and the exact value is {CONTROL_SMOOTH_EXACT} at every \
                 spacing because the stencil factorises on a separable input. \
                 Relative error {relative:.3e} exceeds \
                 {CONTROL_SMOOTH_TOLERANCE:e}, so the instrument is wrong and \
                 every norm this bench reports is measuring it and not the field",
                r.domain_norm,
                r.samples
            );
        }

        let mut worst_kink = 0.0f64;
        for r in &kink {
            let measured = r.domain_norm * r.h.powi(5) / CONTROL_KINK_CONSTANT;
            let relative = (measured - 1.0).abs();
            worst_kink = worst_kink.max(relative);
            assert!(
                relative <= CONTROL_KINK_TOLERANCE,
                "VOID: the stencil reads {:.12e} on max(x+y+z, 0) at {}^3 and the \
                 exact value is {CONTROL_KINK_CONSTANT}/h^5 = {:.12e}. A Lipschitz \
                 kink is the weakest singularity C1 has to detect, so an \
                 instrument that mis-scales it cannot distinguish a sharp field \
                 from a smooth one",
                r.domain_norm,
                r.samples,
                CONTROL_KINK_CONSTANT / r.h.powi(5)
            );
        }
        let kink_exponent = growth_exponent(&kink, |r| r.domain_norm);
        assert!(
            (kink_exponent - CONTROL_KINK_EXPONENT).abs() <= CONTROL_EXPONENT_TOLERANCE,
            "VOID: the fitted growth exponent of a Lipschitz kink is \
             {kink_exponent:.6} and the analytic value is \
             {CONTROL_KINK_EXPONENT}. Every C1 verdict is read off this same \
             least-squares fit, so a fit that cannot recover an exact power of \
             h is not evidence about any field"
        );
        let smooth_value = smooth
            .last()
            .expect("the calibration ladder is not empty")
            .domain_norm;

        // ── C1 and C2 on all eight fields ──────────────────────────────────
        let mut measured: Vec<Measured> = Vec::new();
        isomesh::for_each_reference_field!(f64, |name, field| {
            let m = measure(name, &field);
            println!(
                "{:>15}  norm {:.4e} exp {:+.3} ({}) | domain {:.4e} exp {:+.3} | \
                 H full {:.4e} sparse {:.4e} floor {:.4e} | {:.1}s",
                m.name,
                m.c1.norm,
                m.c1.exponent,
                m.c1.regime(),
                m.c1.domain_norm,
                m.c1.domain_exponent,
                m.full.hausdorff,
                m.headline().arm.hausdorff,
                m.extraction_floor.hausdorff,
                m.wall_seconds,
            );
            measured.push(m);
        });

        // ── vacuity controls over the roster ───────────────────────────────
        assert!(
            measured.len() == FIELD_COUNT,
            "VOID: {} fields were measured and the roster has {FIELD_COUNT}, so \
             C1's claim that the norm is measured on all eight is not what this \
             CSV reports",
            measured.len()
        );
        for m in &measured {
            for r in &m.rungs {
                assert!(
                    r.band_points > 0,
                    "VOID: the |u| <= {BAND_CELLS}h band is empty on {} at {}^3, \
                     so its mixed-derivative norm is a maximum over nothing \
                     rather than a norm (M-44)",
                    m.name,
                    r.samples
                );
            }
            assert!(
                m.full.triangles > 0 && m.full.coverage,
                "VOID: the full grid on {} emitted {} triangles with coverage \
                 {}, so hausdorff_full is a distance to nothing and C2 has no \
                 baseline to match on this field",
                m.name,
                m.full.triangles,
                m.full.coverage
            );
        }

        let sphere = measured
            .iter()
            .find(|m| m.name == "sphere")
            .expect("sphere is the first name in the roster");
        assert!(
            sphere.c1.finite,
            "VOID: sphere reports an unbounded mixed-derivative norm \
             ({:.6e} at the finest rung, growth exponent {:+.6}, regime {}), so \
             the measurement is broken rather than the fields being rough — the \
             registered vacuity control for this row",
            sphere.c1.norm,
            sphere.c1.exponent,
            sphere.c1.regime()
        );

        let qualifying: Vec<&'static str> = measured
            .iter()
            .filter(|m| m.c1.finite)
            .map(|m| m.name)
            .collect();
        assert!(
            !qualifying.is_empty(),
            "VOID: no field shows a bounded norm, so C2 is scoped to the empty \
             set and its verdict is vacuous whichever way it comes out"
        );
        let fields_qualifying = qualifying.join("|");

        let head_ratio = sphere.point_ratio();
        assert!(
            head_ratio >= POINT_RATIO_FLOOR,
            "VOID: the headline sparse level offers {head_ratio:.4}x fewer \
             samples and C2 registers at least {POINT_RATIO_FLOOR}x, so a \
             negative verdict would be a statement about this harness's budget \
             rather than about the method"
        );

        // ── verdicts ───────────────────────────────────────────────────────
        let box_exact = measured
            .iter()
            .find(|m| m.name == "box_exact")
            .expect("box_exact is in the roster");
        let csg = measured
            .iter()
            .find(|m| m.name == "csg_difference")
            .expect("csg_difference is in the roster");
        let c1_holds = measured.len() == FIELD_COUNT && !box_exact.c1.finite && !csg.c1.finite;

        let c2_of = |m: &Measured| {
            let head = m.headline();
            m.c1.finite
                && head.arm.coverage
                && m.c2_reachable()
                && m.point_ratio() >= POINT_RATIO_FLOOR
                && head.arm.hausdorff <= m.full.hausdorff
        };
        let c2_global = measured.iter().filter(|m| m.c1.finite).all(c2_of);

        println!(
            "\ncalibration: x^2y^2z^2 -> {smooth_value:.12} (exact \
             {CONTROL_SMOOTH_EXACT}, worst relative {worst_smooth:.3e}); \
             max(x+y+z,0) -> 6/h^5 to {worst_kink:.3e}, exponent \
             {kink_exponent:.6}"
        );
        println!("qualifying (C1's output): {fields_qualifying}");
        println!("C1 {c1_holds}  C2 {c2_global}");

        for m in &measured {
            let head = m.headline();
            let rate = ladder_rate(&m.ladder);
            run.record(&[
                ("field", m.name.to_string()),
                ("mixed_derivative_norm", format!("{:.6e}", m.c1.norm)),
                ("norm_finite", m.c1.finite.to_string()),
                ("sparse_grid_points", head.points.to_string()),
                ("full_grid_points", m.full_points().to_string()),
                ("point_ratio", format!("{:.6}", m.point_ratio())),
                ("hausdorff_sparse", format!("{:.6e}", head.arm.hausdorff)),
                ("hausdorff_full", format!("{:.6e}", m.full.hausdorff)),
                ("fields_qualifying", fields_qualifying.clone()),
                ("c1_holds", c1_holds.to_string()),
                ("c2_holds", c2_of(m).to_string()),
                // ── extras (M-273) ──
                (
                    "band_points_finest",
                    m.rungs
                        .last()
                        .expect("the ladder is not empty")
                        .band_points
                        .to_string(),
                ),
                ("bound_exact", m.bound_exact.to_string()),
                (
                    "break_even_reads_per_eval",
                    format!("{:.1}", m.break_even_reads_per_eval()),
                ),
                ("c2_global", c2_global.to_string()),
                ("c2_reachable", m.c2_reachable().to_string()),
                ("control_kink_constant", format!("{worst_kink:.6e}")),
                ("control_kink_exponent", format!("{kink_exponent:.6}")),
                ("control_smooth_value", format!("{smooth_value:.9}")),
                ("field_wall_seconds", format!("{:.3}", m.wall_seconds)),
                ("full_triangles", m.full.triangles.to_string()),
                ("hausdorff_ladder", hausdorff_ladder_token(&m.ladder)),
                (
                    "hausdorff_sparse_floor",
                    format!("{:.6e}", m.extraction_floor.hausdorff),
                ),
                (
                    "norm_domain_exponent",
                    format!("{:.6}", m.c1.domain_exponent),
                ),
                ("norm_domain_finest", format!("{:.6e}", m.c1.domain_norm)),
                ("norm_growth_exponent", format!("{:.6}", m.c1.exponent)),
                ("norm_ladder", norm_ladder_token(&m.rungs)),
                (
                    "norm_over_floor_min",
                    format!("{:.6e}", m.c1.over_floor_min),
                ),
                ("norm_regime", m.c1.regime().to_string()),
                ("sparse_coverage", head.arm.coverage.to_string()),
                ("sparse_level", head.level.to_string()),
                ("sparse_points_ladder", points_ladder_token(&m.ladder)),
                ("sparse_rate_exponent", format!("{rate:.6}")),
                ("sparse_rate_theory", format!("{SMOLYAK_THEORY_RATE:.6}")),
                ("sparse_reads_per_sample", head.reads_per_sample.to_string()),
                ("sparse_terms", head.terms.to_string()),
                ("sparse_triangles", head.arm.triangles.to_string()),
                ("sparse_unconverged", head.arm.unconverged.to_string()),
                (
                    "extraction_floor_triangles",
                    m.extraction_floor.triangles.to_string(),
                ),
            ]);
        }
    });
}
