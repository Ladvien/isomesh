//! **P-134 — `|Delta|` normalised as a continuous ambiguity magnitude, instead
//! of a binary test.**
//!
//! Ticket: R-134. Pre-registered before this harness existed.
//!
//! ```bash
//! cargo bench --bench experiment_p134
//! ```
//!
//! Writes `docs/experiments/p-134.csv`.
//!
//! # What was missing
//!
//! The crate computes the body-saddle discriminant on every cell it has ever
//! asked about an interior ambiguity, and then throws away everything except its
//! sign. `BodySaddles::roots`
//! (`crates/isomesh/src/marching_cubes/trilinear.rs:236-267`) reads
//! `discriminant < 0` (`:247`), `discriminant == 0` (`:250`) and otherwise, and
//! the number itself never leaves the function: `roots` is private, and the
//! public surface exposes `inside_mask()` (`:288`), `inside_count()` (`:294`),
//! `has_inner_hexagon()` (`:306`) and `interior_vertex()` (`:815`) — four
//! booleans and a point. A `grep` for `discriminant` over `src/` finds the one
//! line at `:246` and nothing that records it.
//!
//! So the shipped classifier is `sign(Delta)`, and the *distance from the
//! decision boundary* — the thing every other adaptive criterion in the crate is
//! built out of — has never been read. `✗12`'s removal of the rank branch from
//! the QEF (`dual_contouring/solve.rs:20-34`) is the same shape of gap one stage
//! later: a threshold with no measurement of how close each cell sits to it.
//!
//! **P-127 is what licenses treating `Delta` as a classical object rather than a
//! transcribed quadratic.** `docs/experiments/p-127.csv` records
//! `symbolic_difference_is_zero=true`, `terms_disc=12`, `terms_cayley=12`,
//! `total_degree=4` and `pencil_matches=3`: `b*b - 4*a*c` from
//! `BodySaddles::coefficients` **is** Cayley's `2x2x2` hyperdeterminant of the
//! eight corner values under `f[u + 2v + 4w]`, with no sign correction. Cayley's
//! form is homogeneous of degree 4 and `GL(2)^3`-relatively invariant with weight
//! `(det g1 det g2 det g3)^2`, so `|Delta| / (max|f_i|)^4` is scale-free by
//! construction rather than by hope, and this row is that construction measured.
//!
//! **The wave-1 number this row is built on.** P-127's C3 reports
//! `f32_sign_disagreements = 14` over `random_rational_trials = 3481` —
//! **non-zero**, so R-134 is *not* retired by its own registered falsifier
//! (*"C3 by zero f32 disagreements, which would retire P-134 before it is
//! written"*). Two further columns of that CSV decide this file's headline
//! threshold: all fourteen disagreements sit in the near-zero stratum
//! (`f32_disagreements_near_zero = 14`, `f32_disagreements_general = 0`, out of
//! `near_zero_constructed = 481` tuples bracketed within `1e-6` of a sign
//! change), and `f64_sign_disagreements = 0`. The failure of the sign is
//! therefore *entirely* a small-magnitude phenomenon — which is exactly the claim
//! a magnitude is supposed to make legible.
//!
//! # The normalisation, and where its threshold comes from
//!
//! ```text
//! m(cell) = |Delta(f0..f7)| / (max_i |f_i|)^4          the `normalisation` column
//! ```
//!
//! Cayley's twelve coefficients are `1` four times, `-2` six times and `+4`
//! twice, so their `l1` norm is `4 + 12 + 8 = 24` — read off
//! `common::poly::cayley_2x2x2()` on every run rather than quoted
//! (`cayley_l1_norm`). Since every `|f_i| / max|f_i| <= 1`, the normalised
//! term-magnitude sum
//!
//! ```text
//! S(cell) = sum_j |c_j| * prod_i (|f_i| / max|f_i|)^(e_ji)      in [0, 24]
//! ```
//!
//! bounds `m` above and bounds its rounding below. `Poly::eval_f64` spends at
//! most four multiplications per monomial and eleven additions over twelve of
//! them, so after normalisation one evaluation errs by at most `15u * S` with
//! `u = f64::EPSILON / 2`; two *different* routes to the same cell therefore
//! disagree by at most `30u * S`. Non-dyadic rescaling adds up to `4u * S` for
//! the perturbed inputs (degree four) and up to `6u * S` for the denominator's
//! fourth power and the division, so `40u * S` is the whole envelope and this file
//! records the disagreement in units of
//!
//! ```text
//! ROUNDING_ENVELOPE * u * S,        ROUNDING_ENVELOPE = 64
//! ```
//!
//! — a factor of 1.6 of headroom over an already worst-case bound.
//! `scale_error_rounding_units` is that ratio and C1's numeric half is that it
//! stays at or below `1`.
//!
//! The same arithmetic in `f32` gives `30 * (f32::EPSILON / 2) * 24 = 4.29e-5`, so
//! a cell whose `m` is below about `5e-5` has an `f32` discriminant whose **sign**
//! is not determined by its inputs. Rounded up to the decade,
//!
//! ```text
//! HEADLINE_THRESHOLD = 1e-4
//! ```
//!
//! is what the registered `cells_above_threshold` counts above: above it the `f32`
//! sign is reliable, below it the cell is in the stratum where P-127's fourteen
//! disagreements live. `f32_disagreements_above_threshold` is that prediction
//! measured on real field data and is expected to be **0**, and
//! `f32_rounding_floor` re-derives `30 * u_f32 * 24` on every row so the threshold
//! is arithmetic rather than memory.
//!
//! # Per-cell symmetric Hausdorff, built here because the shipped one is global
//!
//! `isomesh::validate::accuracy` reports `symmetric_hausdorff()` over a **whole
//! mesh** (`validate/accuracy.rs:236`), and C2 is a statement about the ordering
//! of *cells*. So the per-cell figure is constructed here, in two directions that
//! are deliberately not the same mechanism:
//!
//! - **mesh -> true surface.** Seven points per triangle — three vertices, three
//!   edge midpoints, the centroid — and at each the **first-order** distance to
//!   the zero set, `|f(p)| / ||grad f(p)||`. First order, and this file says so
//!   rather than calling it a distance: it is the leading Taylor term, exact only
//!   in the limit, and its second-order error grows with the curvature of `f`.
//!   Points where `||grad f||` falls below `GRAD_FLOOR = 1e-12` carry no direction
//!   and are skipped and counted (`gradient_floor_points`).
//! - **true surface -> mesh.** Bisected zero crossings on the cell's **28 straight
//!   segments**: the 12 cube edges from the public `EDGE_CORNERS`, the 12 face
//!   diagonals from the public `face_corners`, and the 4 body diagonals. Every
//!   segment whose endpoints straddle `is_inside` gets `BISECTIONS = 16`
//!   halvings, which pins the crossing to `h / 65536`, and the resulting point is
//!   **on the true zero set and inside the closed cell by construction**. Its
//!   distance to the cell's triangles is then exact point-triangle geometry.
//!
//! The per-cell symmetric Hausdorff is the maximum of the two directions, and
//! both directions' medians are recorded beside it.
//!
//! **Why the second direction is bisection and not a projected sub-grid.** A
//! Newton projection of interior seed points needs a convergence test, an in-cell
//! test and an acceptance rate — three decisions about the instrument that would
//! land inside C2's number. Bisection on a straddling segment needs none: the
//! bracket exists before the first sample, halving cannot fail, the point cannot
//! leave the cell, and no gradient is read. And the twelve cube edges are the
//! *same twelve* whose linear interpolation the extractor placed its vertices on,
//! so direction B's floor is precisely the interpolation error this row is about.
//! The cube's 12 edges connect all 8 corners, so a cell with a sign change always
//! has at least one straddling segment: direction B is never empty on a surface
//! cell, and that is asserted rather than assumed
//! (`cells_without_true_surface_points`).
//!
//! A segment that crosses three times contributes the one crossing bisection
//! converges to. Stated because it is a real restriction: direction B is a lower
//! bound on the true one-sided Hausdorff over a non-monotone segment, never an
//! over-report.
//!
//! **The `bound()` warning, honoured explicitly.** Only direction A reads field
//! values as a length, and gotcha 15 of the API cheat sheet is that `|sample|` is
//! meaningless where `field.bound()` is not `Exact` — `csg_difference` is
//! `Underestimate { q: 0.5 }`, `gyroid` is `Lipschitz { l: 3.4641... }`, and
//! `fbm_terrain` and `noise_cavity` are `Unbounded`. Dividing by `||grad f||` is
//! what makes direction A a *first-order distance* rather than a value reading,
//! and that estimate is correct for any implicit function with a non-degenerate
//! gradient — but it degrades exactly where the gradient varies fastest, which is
//! where the bound is weakest. So every row carries `field_bound` and
//! `distance_reading_is_exact`, and `gradient_norm_median` is the calibration: on
//! an `Exact` field it must read `1` to rounding, and how far from `1` it reads on
//! the other four is the size of the caveat, in the CSV, per field. Direction B
//! carries no such caveat at all — it is a distance between two points in space.
//!
//! # `rank_correlation_with_defects`
//!
//! The defect count is the shipped validator's, per cell, and not a second
//! opinion: `validate::validate_features` returns `NonManifoldFeatures`
//! (`validate.rs:597-612`) — non-manifold edges, non-manifold vertices,
//! inconsistently-oriented edges and boundary edges, each as indices — and each
//! feature is attributed to the cell containing its midpoint (edges) or its
//! position (vertices). That is the `T-001` family of defect, per cell, from the
//! pass that defines it.
//!
//! Boundary edges are counted in and reported separately as `boundary_features`,
//! because `fbm_terrain` has `closed_in_domain() == false` and its boundary edges
//! are the field leaving the box rather than a defect. A reader subtracting them
//! has the non-boundary counts in the same row (`nonmanifold_features`,
//! `misoriented_features`).
//!
//! Two further correlations sit beside it, because a clean mesh gives a dead
//! defect column and a dead column makes a `0.0` unreadable:
//! `rank_correlation_with_ambiguous_faces` against the popcount of
//! `AMBIGUOUS_FACES[case]` (the shipped face-ambiguity classifier — a pure
//! function of the sign pattern, and therefore **not** a function of `|Delta|`),
//! and `rank_correlation_with_uncertified` against the negation of
//! `validate::cell_is_certified` (Plantinga & Vegter's certificate,
//! `validate/isotopy.rs:126`). `defect_variance` and `defect_cells` say on every
//! row whether the primary column could have moved at all.
//!
//! # Arms
//!
//! One extraction, one grid, one shared surface-cell population per row (M-281).
//! The arms are the transform families C1 tests and the two directions C2's error
//! is built from.
//!
//! | arm | what it varies | is_control |
//! |---|---|---|
//! | magnitude | `m = abs(Delta) / max(abs f)^4` on every surface cell | no — it is the row's subject |
//! | dyadic scaling | 4 scale factors that are exact powers of two | **yes** — `scale_error_dyadic`, predicted exactly `0` |
//! | non-dyadic scaling | 8 scale factors spanning six decades | no — C1's rounding arm |
//! | octahedral relabelling | all 48 of `common::poly::octahedral_relabellings()` | no — C1's second invariance |
//! | Hausdorff direction A | mesh points, first-order distance to the zero set | no |
//! | Hausdorff direction B | bisected true-surface points, exact distance to the mesh | **yes** — `cells_where_b_dominates` proves the symmetric figure is two-sided |
//! | shipped evaluation route | `b*b - 4*a*c` from `BodySaddles::coefficients` beside `Poly::eval_f64` | **yes** — `shipped_route_sign_disagreements` |
//! | `f32` route | `Poly::eval_f32` on `f32`-rounded corners | **yes** — `f32_disagreements_above_threshold`, predicted `0` |
//! | gradient calibration | `norm(grad f)` on the `Exact` fields | **yes** — `gradient_norm_median` must read `1` |
//!
//! `resolution` counts **samples**, so `n` samples span `n - 1` cells
//! (`benches/common/mod.rs:40-43`) and `cells = (n - 1)^3` is that arithmetic and
//! nothing else. **Two resolutions, 33 and 65**, the authoring contract's default
//! ladder; the registration names none. `129` is left out because the per-cell
//! Hausdorff spends about 340 field samples per surface cell and the surface-cell
//! count is quadratic in the resolution, so a third rung would quadruple the run
//! for a verdict that turns on a rank correlation rather than on a rate.
//!
//! # Which reading carries each clause
//!
//! - **C1 is per-row.** `scale_error_dyadic == 0` exactly, **and**
//!   `scale_error_rounding_units <= 1`, over that row's whole surface-cell
//!   population and all 60 transforms. The two halves answer different questions
//!   and that is why both are recorded: dyadic rescaling of a binary float is
//!   *exact* — the mantissas never change — so a non-zero there is genuine scale
//!   dependence, the registered falsifier, while the non-dyadic and octahedral
//!   arms can only ever be rounding and the envelope is what says so.
//! - **C2 is global**, because the clause counts fields: *at least four of eight*
//!   above `0.5`. The same boolean is written on every row, decided at
//!   `resolution = 65`, the finer grid. `c2_fields_above_bar` records the count
//!   that decided it and `rank_correlation_with_hausdorff` is per-row throughout,
//!   so the file settles the coarser grid too.
//!
//! # SHARE, recomputed before the numbers
//!
//! The registration's SHARE is *"C2 offers a criterion for the adaptivity stage;
//! it does not by itself move any cost."* Discharged, and the arithmetic is why it
//! is worth writing down: the magnitude costs **no field samples at all**. From
//! eight corner values already in hand, `BodySaddles::coefficients` is
//! `3 + 3 + 4 + 3 + 7 + 3 = 23` arithmetic operations
//! (`trilinear.rs:202-213`) and `b*b - 4*a*c` is 5 more (`:246`); `max_i |f_i|` is
//! 8 absolute values and 7 maxima, its fourth power 2 multiplications, and the
//! normalisation 1 division. **46 operations, recorded as
//! `normalised_magnitude_flops`**, against the eight `Sdf::sample` calls the cell
//! has already paid for. So there is no share to move and none is claimed; what
//! C2 decides is whether those 46 operations *order cells* the way geometric error
//! does, which is a criterion and not a saving.
//!
//! # Determinism
//!
//! One thread, no PRNG, `f64` throughout, `x`-innermost sweep order fixed. The
//! only sorts are `f64::total_cmp` for the quantiles and an integer sort on
//! `(cell, triangle)` keys. `Poly` keeps its coefficients in a `BTreeMap`, so its
//! summation order is fixed by the exponent vectors and identical on every host
//! (`common/poly.rs:270-277`).
//!
//! # Vacuity controls
//!
//! `M-44`: a zero that could not have been non-zero is not a measurement. Each of
//! these runs before the first `run.record` and panics with a message starting
//! `VOID: `.
//!
//! - **The registered one.** `hausdorff_variance > 0` on **every** field and
//!   resolution — the registration's own control, verbatim: *"the per-cell
//!   Hausdorff column must have non-zero variance on every field, or the
//!   correlation is against a constant"*. Recorded as a column.
//! - **`delta_magnitude_variance > 0` on every row.** The same failure on the
//!   other side of the pairing: a constant magnitude column would make every
//!   correlation in the row an artefact of tie-breaking.
//! - **`paired_cells >= 3` on every row**, because
//!   `common::beta::rank_correlation` *defines* a coefficient over fewer than
//!   three pairs to be `0.0` rather than measuring one (`beta.rs:495-497`).
//! - **`min_denominator > 0` on every row.** `m` divides by `(max|f_i|)^4`; a
//!   surface cell has a strictly negative corner so the maximum cannot be zero,
//!   and this is that arithmetic checked rather than argued.
//! - **`cells_without_true_surface_points == 0` on every row.** The cube's 12
//!   edges connect all 8 corners, so a cell with a sign change must have a
//!   straddling segment; a non-zero here means direction B silently contributed
//!   nothing and the "symmetric" Hausdorff is one-sided.
//! - **Cayley is the polynomial this row thinks it is**: 12 terms, total degree 4,
//!   `l1` norm 24, all three read off the module. Every rounding bound in the file
//!   is written in units of that `l1` norm, so a change to the polynomial would
//!   otherwise silently loosen C1's envelope.
//! - **`cells_where_b_dominates > 0` across the run** (global). Direction B losing
//!   on a given field is that field's geometry; never winning anywhere means the
//!   symmetric figure is direction A wearing a longer name and the second
//!   mechanism was never at risk.
//! - **`nonzero_transform_deviations > 0` across the run** (global). If not one of
//!   the 60 transforms on any cell moved `m` at all, C1's invariance held without
//!   the arithmetic ever being exercised.
//! - **`defect_cells > 0` across the run** (global, not per row). A clean mesh is
//!   the answer on most fields and a per-row assertion would abort on it;
//!   `fbm_terrain` is open in its domain, so its boundary edges guarantee the
//!   column is live somewhere — which is what licenses reading a `0.0` elsewhere
//!   as that row's geometry rather than a dead instrument.

#![allow(
    clippy::cast_possible_truncation,
    clippy::cast_precision_loss,
    clippy::float_cmp,
    clippy::too_many_lines
)]

mod common;

use common::poly::Poly;
use isomesh::fields::{FieldBound, ReferenceField};
use isomesh::marching_cubes::MarchingCubes;
use isomesh::marching_cubes::table::{AMBIGUOUS_FACES, EDGE_CORNERS, face_corners, is_inside};
use isomesh::marching_cubes::trilinear::BodySaddles;
use isomesh::validate::{self, ValidateConfig};
use isomesh::{MeshBuffer, Sdf, Shape3};

/// Samples per axis. Two rungs, and the header says why `129` is not a third.
const RESOLUTIONS: [u32; 2] = [33, 65];

/// `for_each_reference_field!` yields eight (`fields/mod.rs:195`).
const FIELDS: usize = 8;

/// C2's bar on the Spearman coefficient against per-cell symmetric Hausdorff.
const C2_BAR: f64 = 0.5;

/// C2's count: at least this many of the eight fields must clear [`C2_BAR`].
const C2_MIN_FIELDS: usize = 4;

/// The resolution C2's global verdict is read at — the finer of the two.
const C2_RESOLUTION: u32 = 65;

/// The normalised-magnitude thresholds `threshold_sweep` reports counts at.
const SWEEP: [f64; 5] = [1e-8, 1e-6, 1e-4, 1e-2, 1e-1];

/// The threshold the registered `cells_above_threshold` counts above.
///
/// Derived, not chosen: `30 * (f32::EPSILON / 2) * 24 = 4.29e-5` is the `f32`
/// rounding floor of `m`, so below about `5e-5` the `f32` sign of `Delta` is not
/// determined by its inputs. Rounded up to the decade. `f32_rounding_floor`
/// re-derives it on every row.
const HEADLINE_THRESHOLD: f64 = 1e-4;

/// Scale factors that are exact powers of two.
///
/// Multiplying a binary float by `2^k` changes only its exponent, so every
/// mantissa operation downstream is bit-identical and `m` is **exactly**
/// invariant. A non-zero `scale_error_dyadic` is therefore genuine scale
/// dependence rather than rounding, which is the discriminator C1's falsifier
/// needs.
const DYADIC_SCALES: [f64; 4] = [1.0 / 1024.0, 0.25, 4.0, 1024.0];

/// Scale factors that are not powers of two, spanning six decades.
///
/// The arm that actually rounds, and the one `scale_error_rounding_units` is
/// about.
const ROUGH_SCALES: [f64; 8] = [1e-3, 0.1, 1.0 / 3.0, 0.7, 3.0, 7.0, 10.0, 1e3];

/// The `l1` norm of Cayley's twelve coefficients: `4*1 + 6*2 + 2*4`.
///
/// Read off `cayley_2x2x2()` and asserted rather than trusted, because the whole
/// rounding envelope is written in units of it.
const CAYLEY_L1: i128 = 24;

/// Units of `u * S` the evaluation routes may disagree by before C1 fails.
///
/// One `Poly::eval_f64` errs by at most `15u * S` after normalisation, two routes
/// by `30u * S`, plus `4u * S` for the perturbed inputs of a non-dyadic rescaling
/// and `6u * S` for the denominator's fourth power and the division: `40u * S`
/// worst case, and `64` leaves a factor of 1.6 over it.
const ROUNDING_ENVELOPE: f64 = 64.0;

/// `f64`'s unit roundoff, the half-epsilon the error bounds are written in.
const UNIT_ROUNDOFF: f64 = f64::EPSILON * 0.5;

/// `f32`'s unit roundoff, for the `f32` sign floor.
const UNIT_ROUNDOFF_F32: f64 = f32::EPSILON as f64 * 0.5;

/// The forward-error constant of a twelve-term degree-four evaluation, doubled
/// for the two routes being compared.
const TWO_ROUTE_UNITS: f64 = 30.0;

/// Bisection halvings per straddling segment: pins the crossing to `h / 65536`.
const BISECTIONS: u32 = 16;

/// Below this a gradient norm carries no direction and the first-order distance
/// is not defined.
const GRAD_FLOOR: f64 = 1e-12;

/// Straight segments per cell whose crossings are direction B's true-surface
/// points: 12 cube edges, 12 face diagonals, 4 body diagonals.
const SEGMENT_COUNT: usize = 28;

/// Arithmetic operations the shipped route spends on `m` from eight corner values
/// already in hand. Derived in the header's SHARE section.
const NORMALISED_MAGNITUDE_FLOPS: u32 = 46;

/// The local offset of a corner, as grid steps.
///
/// `crate::cube::corner_offset` is `pub(crate)` and not re-exported
/// (`marching_cubes/table.rs:88-91` omits it), so these are the three lines
/// `cube.rs:149-155` spells out, in the crate's own `f[u + 2v + 4w]` order.
fn corner_offset(corner: u8) -> [usize; 3] {
    [
        usize::from(corner & 1),
        usize::from((corner >> 1) & 1),
        usize::from((corner >> 2) & 1),
    ]
}

/// `a - b`, componentwise.
fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

/// `a . b`.
fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

/// `a x b`.
fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

/// `k * a`, componentwise.
fn scaled(a: [f64; 3], k: f64) -> [f64; 3] {
    [a[0] * k, a[1] * k, a[2] * k]
}

/// `||a||`.
fn norm(a: [f64; 3]) -> f64 {
    dot(a, a).sqrt()
}

/// Componentwise midpoint of two points.
fn mid(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[0].midpoint(b[0]),
        a[1].midpoint(b[1]),
        a[2].midpoint(b[2]),
    ]
}

/// Arithmetic mean of three points.
fn centroid(a: [f64; 3], b: [f64; 3], c: [f64; 3]) -> [f64; 3] {
    [
        (a[0] + b[0] + c[0]) / 3.0,
        (a[1] + b[1] + c[1]) / 3.0,
        (a[2] + b[2] + c[2]) / 3.0,
    ]
}

/// `(1 - t) * a + t * b`.
fn lerp(a: [f64; 3], b: [f64; 3], t: f64) -> [f64; 3] {
    [
        a[0] + (b[0] - a[0]) * t,
        a[1] + (b[1] - a[1]) * t,
        a[2] + (b[2] - a[2]) * t,
    ]
}

/// Distance from `p` to the **closed** triangle `abc`.
///
/// The minimum over three vertices, three clamped edge projections and — only
/// when the plane projection lands inside — the plane projection. One path and no
/// region case analysis: a degenerate triangle simply contributes no interior
/// candidate, and the clamped edge projections already cover its whole extent, so
/// there is nothing to fall back to.
fn point_triangle_distance(p: [f64; 3], a: [f64; 3], b: [f64; 3], c: [f64; 3]) -> f64 {
    let mut best = f64::INFINITY;
    for v in [a, b, c] {
        best = best.min(norm(sub(p, v)));
    }
    for (u, w) in [(a, b), (b, c), (c, a)] {
        let e = sub(w, u);
        let ee = dot(e, e);
        if ee > 0.0 {
            let t = (dot(sub(p, u), e) / ee).clamp(0.0, 1.0);
            best = best.min(norm(sub(sub(p, u), scaled(e, t))));
        }
    }
    let n = cross(sub(b, a), sub(c, a));
    let nn = dot(n, n);
    if nn > 0.0 {
        let q = sub(p, scaled(n, dot(n, sub(p, a)) / nn));
        let inside = dot(cross(sub(b, a), sub(q, a)), n) >= 0.0
            && dot(cross(sub(c, b), sub(q, b)), n) >= 0.0
            && dot(cross(sub(a, c), sub(q, c)), n) >= 0.0;
        if inside {
            best = best.min(norm(sub(p, q)));
        }
    }
    best
}

/// The 28 straight segments of a cell, as corner-index pairs.
///
/// Generated from the public tables rather than transcribed: the 12 edges are
/// `EDGE_CORNERS`, the 12 face diagonals are entries `0,2` and `1,3` of each
/// `face_corners(axis, side)` quad — which `cube.rs:88-101` returns in cyclic
/// order, so those two pairs are the diagonals — and the 4 body diagonals are the
/// antipodal pairs `k, 7 - k`.
fn segments() -> [[u8; 2]; SEGMENT_COUNT] {
    let mut out = [[0u8; 2]; SEGMENT_COUNT];
    let mut n = 0;
    for e in EDGE_CORNERS {
        out[n] = e;
        n += 1;
    }
    for axis in 0..3usize {
        for side in 0..2u8 {
            let q = face_corners(axis, side);
            out[n] = [q[0], q[2]];
            n += 1;
            out[n] = [q[1], q[3]];
            n += 1;
        }
    }
    for k in 0..4u8 {
        out[n] = [k, 7 - k];
        n += 1;
    }
    assert!(
        n == SEGMENT_COUNT,
        "12 edges + 12 face diagonals + 4 body diagonals is {SEGMENT_COUNT}, built {n}"
    );
    out
}

/// Cayley's form with every coefficient replaced by its magnitude.
///
/// The term-magnitude sum `S` is this polynomial evaluated at `|f_i| / max|f_i|`.
/// Because Cayley is octahedrally invariant the group permutes its monomials
/// among themselves with coefficients preserved, so this is invariant too and `S`
/// is one number per cell rather than one per relabelling.
fn abs_cayley(cayley: &Poly) -> Poly {
    let mut out = Poly::zero();
    for (exp, c) in cayley.monomials() {
        out = out.add(&Poly::monomial(exp, c.abs()));
    }
    out
}

/// `max_i |f_i|`. Exact: `abs` and `max` are both exact on binary floats.
fn max_abs(f: &[f64; 8]) -> f64 {
    f.iter().fold(0.0f64, |acc, v| acc.max(v.abs()))
}

/// `Delta`, its normalisation `m = |Delta| / (max|f_i|)^4`, and that denominator's
/// base `max_i |f_i|`.
///
/// `(max|f|)^4` is formed as `(M*M)*(M*M)` so that a dyadic rescaling of `f`
/// leaves every mantissa in the computation untouched — which is what makes
/// `scale_error_dyadic` an exact-zero prediction rather than a tolerance.
fn magnitude(cayley: &Poly, f: &[f64; 8]) -> (f64, f64, f64) {
    let base = max_abs(f);
    let square = base * base;
    let delta = cayley.eval_f64(f);
    (delta, delta.abs() / (square * square), base)
}

/// `b*b - 4*a*c` from the shipped `BodySaddles::coefficients`, as
/// `trilinear.rs:246` forms it.
fn shipped_discriminant(f: &[f64; 8]) -> f64 {
    let [a, b, c] = BodySaddles::<f64>::coefficients(f);
    b * b - 2.0 * 2.0 * a * c
}

/// Worst absolute deviation of `m` over one family of scale factors, and how many
/// of them moved it at all.
fn scale_family(cayley: &Poly, f: &[f64; 8], m0: f64, scales: &[f64]) -> (f64, usize) {
    let mut worst = 0.0f64;
    let mut moved = 0;
    for s in scales {
        let rescaled: [f64; 8] = std::array::from_fn(|i| f[i] * *s);
        let (_, m, _) = magnitude(cayley, &rescaled);
        let dev = (m - m0).abs();
        if dev > 0.0 {
            moved += 1;
        }
        worst = worst.max(dev);
    }
    (worst, moved)
}

/// Worst absolute deviation of `m` over the 48 octahedral relabellings, and how
/// many of them moved it at all.
fn relabel_family(
    cayley: &Poly,
    f: &[f64; 8],
    m0: f64,
    perms: &[[u8; 8]; 48],
) -> (f64, usize) {
    let mut worst = 0.0f64;
    let mut moved = 0;
    for perm in perms {
        let relabelled = common::poly::relabel(perm, f);
        let (_, m, _) = magnitude(cayley, &relabelled);
        let dev = (m - m0).abs();
        if dev > 0.0 {
            moved += 1;
        }
        worst = worst.max(dev);
    }
    (worst, moved)
}

/// C1's whole measurement for one cell, over all 60 transforms.
#[derive(Clone, Copy, Debug, Default)]
struct Invariance {
    /// Worst absolute deviation of `m` under the dyadic scale factors.
    dyadic: f64,
    /// Worst absolute deviation under the non-dyadic scale factors.
    rough: f64,
    /// Worst absolute deviation under the 48 octahedral relabellings.
    octahedral: f64,
    /// The worst of the three in units of `ROUNDING_ENVELOPE * u * S`.
    rounding_units: f64,
    /// How many of the 60 transforms moved `m` at all.
    moved: usize,
}

/// Run every transform on one cell.
///
/// A cell whose normalised term-magnitude sum `S` is exactly zero has at most one
/// non-zero corner, hence `Delta == 0` identically under every transform; it
/// contributes zero deviations and no rounding-unit ratio, which is the honest
/// reading rather than a `0/0`.
fn invariance_of(
    cayley: &Poly,
    abs_poly: &Poly,
    perms: &[[u8; 8]; 48],
    f: &[f64; 8],
    m0: f64,
    base: f64,
) -> Invariance {
    let unit: [f64; 8] = std::array::from_fn(|i| f[i].abs() / base);
    let s_norm = abs_poly.eval_f64(&unit);

    let (dyadic, moved_dyadic) = scale_family(cayley, f, m0, &DYADIC_SCALES);
    let (rough, moved_rough) = scale_family(cayley, f, m0, &ROUGH_SCALES);
    let (octahedral, moved_perms) = relabel_family(cayley, f, m0, perms);

    let worst = dyadic.max(rough).max(octahedral);
    Invariance {
        dyadic,
        rough,
        octahedral,
        rounding_units: if s_norm > 0.0 {
            worst / (ROUNDING_ENVELOPE * UNIT_ROUNDOFF * s_norm)
        } else {
            0.0
        },
        moved: moved_dyadic + moved_rough + moved_perms,
    }
}

/// Nearest-rank quantile of an already-`total_cmp`-sorted slice.
///
/// `index = floor(q * (len - 1))`, no interpolation, so the value returned is one
/// some cell actually took. `0.0` for an empty slice, which cannot reach a
/// recorded row: `paired_cells >= 3` is asserted first.
fn quantile(sorted: &[f64], q: f64) -> f64 {
    if sorted.is_empty() {
        return 0.0;
    }
    let last = sorted.len() - 1;
    sorted[((q * last as f64).floor() as usize).min(last)]
}

/// Arithmetic mean; `0.0` for an empty slice.
fn mean(values: &[f64]) -> f64 {
    if values.is_empty() {
        return 0.0;
    }
    values.iter().sum::<f64>() / values.len() as f64
}

/// Population variance. This is the registered vacuity control's quantity.
fn variance(values: &[f64]) -> f64 {
    if values.is_empty() {
        return 0.0;
    }
    let mu = mean(values);
    values.iter().map(|v| (v - mu) * (v - mu)).sum::<f64>() / values.len() as f64
}

/// Median of a slice, by a `total_cmp` sort of a copy.
fn median(values: &[f64]) -> f64 {
    let mut sorted = values.to_vec();
    sorted.sort_unstable_by(f64::total_cmp);
    quantile(&sorted, 0.5)
}

/// The `bound()` variant as a CSV-safe token.
fn bound_token(bound: FieldBound) -> &'static str {
    match bound {
        FieldBound::Exact => "exact",
        FieldBound::Lipschitz { .. } => "lipschitz",
        FieldBound::Underestimate { .. } => "underestimate",
        FieldBound::Unbounded => "unbounded",
    }
}

/// One measured `(field, resolution)` pair, before the cross-row verdict exists.
#[derive(Clone, Debug)]
struct Row {
    /// The reference field's `name` literal.
    field: &'static str,
    /// Samples per axis.
    samples: u32,
    /// The grid step this row was measured at.
    cell_size: f64,
    /// `field.bound()` as a token.
    bound: &'static str,
    /// Whether `|f|` may be read as a distance at all.
    bound_is_exact: bool,
    /// `(samples - 1)^3`.
    cells: usize,
    /// Cells whose eight corners do not share a sign.
    surface_cells: usize,
    /// Surface cells that also carry at least one triangle: C2's population.
    paired_cells: usize,
    /// Surface cells with no triangle assigned, excluded from the pairing.
    cells_without_triangles: usize,
    /// Surface cells where no segment straddled. Must be zero.
    cells_without_true_surface_points: usize,
    /// Median normalised magnitude over the surface cells.
    magnitude_median: f64,
    /// Mean of the same column.
    magnitude_mean: f64,
    /// Maximum of the same column.
    magnitude_max: f64,
    /// Minimum of the same column.
    magnitude_min: f64,
    /// Variance of the same column.
    magnitude_variance: f64,
    /// Smallest `max_i |f_i|` over the surface cells: the denominator's floor.
    min_denominator: f64,
    /// Counts above each entry of [`SWEEP`].
    sweep_counts: [usize; SWEEP.len()],
    /// Count above [`HEADLINE_THRESHOLD`].
    cells_above_threshold: usize,
    /// Worst deviation under the dyadic scale factors: predicted exactly `0`.
    scale_error_dyadic: f64,
    /// Worst deviation under the non-dyadic scale factors.
    scale_error_rough: f64,
    /// Worst deviation under the 48 octahedral relabellings.
    scale_error_octahedral: f64,
    /// Worst of the three, absolute: the registered `scale_invariance_error`.
    scale_error_worst: f64,
    /// That worst deviation in units of `ROUNDING_ENVELOPE * u * S`.
    scale_error_rounding_units: f64,
    /// How many `(cell, transform)` pairs moved `m` at all.
    nonzero_transform_deviations: usize,
    /// Median per-cell symmetric Hausdorff, in world units.
    hausdorff_median: f64,
    /// Maximum per-cell symmetric Hausdorff, in world units.
    hausdorff_max: f64,
    /// Variance of the same column. The registered vacuity control's quantity.
    hausdorff_variance: f64,
    /// Median of direction A alone, in world units.
    hausdorff_mesh_to_field_median: f64,
    /// Median of direction B alone, in world units.
    hausdorff_field_to_mesh_median: f64,
    /// Paired cells where direction B exceeded direction A.
    cells_where_b_dominates: usize,
    /// Direction A points skipped for a gradient below [`GRAD_FLOOR`].
    gradient_floor_points: usize,
    /// Median `||grad f||` over the direction A points: the bound calibration.
    gradient_norm_median: f64,
    /// Paired cells carrying at least one validity feature.
    defect_cells: usize,
    /// Total validity features attributed to paired cells.
    defect_total: usize,
    /// Boundary edges over the whole mesh — not a defect on an open field.
    boundary_features: usize,
    /// Non-manifold edges plus non-manifold vertices, whole mesh.
    nonmanifold_features: usize,
    /// Inconsistently-oriented edges, whole mesh.
    misoriented_features: usize,
    /// Variance of the per-cell defect column.
    defect_variance: f64,
    /// Paired cells failing `validate::cell_is_certified`.
    uncertified_cells: usize,
    /// Paired cells with at least one ambiguous face.
    ambiguous_face_cells: usize,
    /// Spearman of `m` against the per-cell symmetric Hausdorff: C2's number.
    rho_hausdorff: f64,
    /// Spearman of `m` against the per-cell validity-feature count.
    rho_defects: f64,
    /// Spearman of `m` against `AMBIGUOUS_FACES[case]`'s popcount.
    rho_ambiguous: f64,
    /// Spearman of `m` against the uncertified indicator.
    rho_uncertified: f64,
    /// Surface cells where `eval_f32` disagrees in sign with `eval_f64`.
    f32_sign_disagreements: usize,
    /// Of those, how many sit above [`HEADLINE_THRESHOLD`]. Predicted `0`.
    f32_disagreements_above_threshold: usize,
    /// Surface cells where the shipped route disagrees in sign with the
    /// polynomial route.
    shipped_sign_disagreements: usize,
    /// Worst `|shipped - polynomial| / (max|f_i|)^4` over the surface cells.
    shipped_max_normalised_gap: f64,
    /// The extracted mesh's vertex count.
    vertices: usize,
    /// The extracted mesh's triangle count.
    triangles: usize,
}

/// Measure one `(field, resolution)`.
///
/// Order: sample the grid once, extract once, run the shipped validator once,
/// bucket the triangles once, then one pass over the cells that does C1's sixty
/// transforms, both Hausdorff directions and the `f32`/shipped-route controls out
/// of the same eight corner values. Nothing here is timed — every clause in this
/// registration is a count, a correlation or an exact-zero prediction, which is
/// the authoring contract's preference and not an omission.
fn measure<F>(
    name: &'static str,
    field: &F,
    samples: u32,
    cayley: &Poly,
    abs_poly: &Poly,
    perms: &[[u8; 8]; 48],
    segs: &[[u8; 2]; SEGMENT_COUNT],
) -> Row
where
    F: ReferenceField + Sdf<Scalar = f64>,
{
    let (shape, origin, h) = common::grid::<f64, _>(field, samples);
    let dim = shape.size();
    let sx = dim[0] as usize;
    let sy = dim[1] as usize;
    let cells = [sx - 1, sy - 1, dim[2] as usize - 1];
    let cell_count = cells[0] * cells[1] * cells[2];

    let mut values = vec![0.0f64; shape.element_count()];
    for z in 0..dim[2] as usize {
        let pz = origin[2] + h * z as f64;
        for y in 0..sy {
            let py = origin[1] + h * y as f64;
            for x in 0..sx {
                let px = origin[0] + h * x as f64;
                values[x + y * sx + z * sx * sy] = field.sample([px, py, pz]);
            }
        }
    }

    let mut mesher = MarchingCubes::<f64>::new();
    let mut mesh = MeshBuffer::<f64>::new();
    mesher
        .extract(field, &shape, origin, h, &mut mesh)
        .expect("marching cubes extraction on a reference field grid");

    let cfg = ValidateConfig::from_cell_size(h).expect("a benchmark cell size is positive");
    let (report, features) = validate::validate_features(&mesh.positions, &mesh.indices, &cfg);

    // Which cell contains a point. Clamped: a mesh vertex sits exactly on a cell
    // boundary by construction, and the boundary belongs to the lower cell.
    let locate = |p: [f64; 3]| -> usize {
        let c: [usize; 3] = std::array::from_fn(|k| {
            let t = ((p[k] - origin[k]) / h).floor();
            if t <= 0.0 {
                0
            } else {
                (t as usize).min(cells[k] - 1)
            }
        });
        c[0] + c[1] * cells[0] + c[2] * cells[0] * cells[1]
    };

    let mut defects = vec![0u32; cell_count];
    for group in [
        features.edges.as_slice(),
        features.boundary_edges.as_slice(),
        features.inconsistently_oriented_edges.as_slice(),
    ] {
        for e in group {
            let a = mesh.positions[e[0] as usize];
            let b = mesh.positions[e[1] as usize];
            defects[locate(mid(a, b))] += 1;
        }
    }
    for v in &features.vertices {
        defects[locate(mesh.positions[*v as usize])] += 1;
    }

    let mut buckets: Vec<(u32, u32)> = Vec::with_capacity(mesh.triangle_count());
    for (t, tri) in mesh.indices.chunks_exact(3).enumerate() {
        let a = mesh.positions[tri[0] as usize];
        let b = mesh.positions[tri[1] as usize];
        let c = mesh.positions[tri[2] as usize];
        buckets.push((locate(centroid(a, b, c)) as u32, t as u32));
    }
    buckets.sort_unstable();

    let mut magnitudes: Vec<f64> = Vec::new();
    let mut paired_mag: Vec<f64> = Vec::new();
    let mut paired_haus: Vec<f64> = Vec::new();
    let mut paired_defect: Vec<f64> = Vec::new();
    let mut paired_ambiguous: Vec<f64> = Vec::new();
    let mut paired_uncertified: Vec<f64> = Vec::new();
    let mut direction_a: Vec<f64> = Vec::new();
    let mut direction_b: Vec<f64> = Vec::new();
    let mut gradient_norms: Vec<f64> = Vec::new();

    let mut row = Row {
        field: name,
        samples,
        cell_size: h,
        bound: bound_token(field.bound()),
        bound_is_exact: field.bound().is_exact(),
        cells: cell_count,
        surface_cells: 0,
        paired_cells: 0,
        cells_without_triangles: 0,
        cells_without_true_surface_points: 0,
        magnitude_median: 0.0,
        magnitude_mean: 0.0,
        magnitude_max: 0.0,
        magnitude_min: 0.0,
        magnitude_variance: 0.0,
        min_denominator: f64::INFINITY,
        sweep_counts: [0; SWEEP.len()],
        cells_above_threshold: 0,
        scale_error_dyadic: 0.0,
        scale_error_rough: 0.0,
        scale_error_octahedral: 0.0,
        scale_error_worst: 0.0,
        scale_error_rounding_units: 0.0,
        nonzero_transform_deviations: 0,
        hausdorff_median: 0.0,
        hausdorff_max: 0.0,
        hausdorff_variance: 0.0,
        hausdorff_mesh_to_field_median: 0.0,
        hausdorff_field_to_mesh_median: 0.0,
        cells_where_b_dominates: 0,
        gradient_floor_points: 0,
        gradient_norm_median: 0.0,
        defect_cells: 0,
        defect_total: 0,
        boundary_features: features.boundary_edges.len(),
        nonmanifold_features: (report.non_manifold_edges + report.non_manifold_vertices) as usize,
        misoriented_features: report.inconsistently_oriented_edges as usize,
        defect_variance: 0.0,
        uncertified_cells: 0,
        ambiguous_face_cells: 0,
        rho_hausdorff: 0.0,
        rho_defects: 0.0,
        rho_ambiguous: 0.0,
        rho_uncertified: 0.0,
        f32_sign_disagreements: 0,
        f32_disagreements_above_threshold: 0,
        shipped_sign_disagreements: 0,
        shipped_max_normalised_gap: 0.0,
        vertices: mesh.vertex_count(),
        triangles: mesh.triangle_count(),
    };

    for cz in 0..cells[2] {
        for cy in 0..cells[1] {
            for cx in 0..cells[0] {
                let mut corner = [0.0f64; 8];
                let mut position = [[0.0f64; 3]; 8];
                let mut case = 0u8;
                for (k, (value, point)) in
                    corner.iter_mut().zip(position.iter_mut()).enumerate()
                {
                    let o = corner_offset(k as u8);
                    let (gx, gy, gz) = (cx + o[0], cy + o[1], cz + o[2]);
                    let v = values[gx + gy * sx + gz * sx * sy];
                    *value = v;
                    *point = [
                        origin[0] + h * gx as f64,
                        origin[1] + h * gy as f64,
                        origin[2] + h * gz as f64,
                    ];
                    if is_inside(v) {
                        case |= 1u8 << k;
                    }
                }
                if case == 0 || case == u8::MAX {
                    continue;
                }
                row.surface_cells += 1;

                let (delta, m0, base) = magnitude(cayley, &corner);
                magnitudes.push(m0);
                row.min_denominator = row.min_denominator.min(base);
                for (slot, threshold) in SWEEP.iter().enumerate() {
                    if m0 > *threshold {
                        row.sweep_counts[slot] += 1;
                    }
                }
                if m0 > HEADLINE_THRESHOLD {
                    row.cells_above_threshold += 1;
                }

                // ── C1: the sixty transforms on this cell ──
                let inv = invariance_of(cayley, abs_poly, perms, &corner, m0, base);
                row.scale_error_dyadic = row.scale_error_dyadic.max(inv.dyadic);
                row.scale_error_rough = row.scale_error_rough.max(inv.rough);
                row.scale_error_octahedral = row.scale_error_octahedral.max(inv.octahedral);
                row.scale_error_rounding_units =
                    row.scale_error_rounding_units.max(inv.rounding_units);
                row.nonzero_transform_deviations += inv.moved;

                // ── the two evaluation-route controls ──
                let f32_corners: [f32; 8] = std::array::from_fn(|i| corner[i] as f32);
                let f32_route = f64::from(cayley.eval_f32(&f32_corners));
                if (f32_route < 0.0) != (delta < 0.0) {
                    row.f32_sign_disagreements += 1;
                    if m0 > HEADLINE_THRESHOLD {
                        row.f32_disagreements_above_threshold += 1;
                    }
                }
                let shipped = shipped_discriminant(&corner);
                if (shipped < 0.0) != (delta < 0.0) {
                    row.shipped_sign_disagreements += 1;
                }
                let square = base * base;
                row.shipped_max_normalised_gap = row
                    .shipped_max_normalised_gap
                    .max((shipped - delta).abs() / (square * square));

                // ── the paired population: cells the mesh actually reached ──
                let key = locate(mid(position[0], position[7])) as u32;
                let lo = buckets.partition_point(|(c, _)| *c < key);
                let hi = buckets.partition_point(|(c, _)| *c <= key);
                if lo == hi {
                    row.cells_without_triangles += 1;
                    continue;
                }

                // Direction A: mesh points, first-order distance to the zero set.
                let mut worst_a = 0.0f64;
                for (_, t) in &buckets[lo..hi] {
                    let tri = &mesh.indices[*t as usize * 3..*t as usize * 3 + 3];
                    let a = mesh.positions[tri[0] as usize];
                    let b = mesh.positions[tri[1] as usize];
                    let c = mesh.positions[tri[2] as usize];
                    for p in [
                        a,
                        b,
                        c,
                        mid(a, b),
                        mid(b, c),
                        mid(c, a),
                        centroid(a, b, c),
                    ] {
                        let gradient_norm = norm(field.gradient(p));
                        if gradient_norm < GRAD_FLOOR {
                            row.gradient_floor_points += 1;
                            continue;
                        }
                        gradient_norms.push(gradient_norm);
                        worst_a = worst_a.max(field.sample(p).abs() / gradient_norm);
                    }
                }

                // Direction B: bisected true-surface points, exact distance to the
                // mesh.
                let mut worst_b = 0.0f64;
                let mut straddling = 0usize;
                for seg in segs {
                    let (ia, ib) = (usize::from(seg[0]), usize::from(seg[1]));
                    if is_inside(corner[ia]) == is_inside(corner[ib]) {
                        continue;
                    }
                    straddling += 1;
                    let (pa, pb) = (position[ia], position[ib]);
                    let mut lo_t = 0.0f64;
                    let mut hi_t = 1.0f64;
                    let mut lo_inside = is_inside(corner[ia]);
                    for _ in 0..BISECTIONS {
                        let tm = lo_t.midpoint(hi_t);
                        let vm = field.sample(lerp(pa, pb, tm));
                        if is_inside(vm) == lo_inside {
                            lo_t = tm;
                            lo_inside = is_inside(vm);
                        } else {
                            hi_t = tm;
                        }
                    }
                    let q = lerp(pa, pb, lo_t.midpoint(hi_t));
                    let mut nearest = f64::INFINITY;
                    for (_, t) in &buckets[lo..hi] {
                        let tri = &mesh.indices[*t as usize * 3..*t as usize * 3 + 3];
                        nearest = nearest.min(point_triangle_distance(
                            q,
                            mesh.positions[tri[0] as usize],
                            mesh.positions[tri[1] as usize],
                            mesh.positions[tri[2] as usize],
                        ));
                    }
                    worst_b = worst_b.max(nearest);
                }
                if straddling == 0 {
                    row.cells_without_true_surface_points += 1;
                    continue;
                }

                row.paired_cells += 1;
                direction_a.push(worst_a);
                direction_b.push(worst_b);
                if worst_b > worst_a {
                    row.cells_where_b_dominates += 1;
                }
                paired_mag.push(m0);
                paired_haus.push(worst_a.max(worst_b));

                let feature_count = defects[key as usize];
                if feature_count > 0 {
                    row.defect_cells += 1;
                    row.defect_total += feature_count as usize;
                }
                paired_defect.push(f64::from(feature_count));

                let ambiguous = AMBIGUOUS_FACES[usize::from(case)].count_ones();
                if ambiguous > 0 {
                    row.ambiguous_face_cells += 1;
                }
                paired_ambiguous.push(f64::from(ambiguous));

                let uncertified = !validate::cell_is_certified(&corner);
                if uncertified {
                    row.uncertified_cells += 1;
                }
                paired_uncertified.push(f64::from(u8::from(uncertified)));
            }
        }
    }

    row.scale_error_worst = row
        .scale_error_dyadic
        .max(row.scale_error_rough)
        .max(row.scale_error_octahedral);

    let mut sorted = magnitudes.clone();
    sorted.sort_unstable_by(f64::total_cmp);
    row.magnitude_median = quantile(&sorted, 0.5);
    row.magnitude_min = sorted.first().copied().unwrap_or(0.0);
    row.magnitude_max = sorted.last().copied().unwrap_or(0.0);
    row.magnitude_mean = mean(&magnitudes);
    row.magnitude_variance = variance(&magnitudes);
    if !row.min_denominator.is_finite() {
        row.min_denominator = 0.0;
    }

    let mut haus_sorted = paired_haus.clone();
    haus_sorted.sort_unstable_by(f64::total_cmp);
    row.hausdorff_median = quantile(&haus_sorted, 0.5);
    row.hausdorff_max = haus_sorted.last().copied().unwrap_or(0.0);
    row.hausdorff_variance = variance(&paired_haus);
    row.hausdorff_mesh_to_field_median = median(&direction_a);
    row.hausdorff_field_to_mesh_median = median(&direction_b);
    row.gradient_norm_median = median(&gradient_norms);

    row.defect_variance = variance(&paired_defect);
    row.rho_hausdorff = common::beta::rank_correlation(&paired_mag, &paired_haus);
    row.rho_defects = common::beta::rank_correlation(&paired_mag, &paired_defect);
    row.rho_ambiguous = common::beta::rank_correlation(&paired_mag, &paired_ambiguous);
    row.rho_uncertified = common::beta::rank_correlation(&paired_mag, &paired_uncertified);
    row
}

fn main() {
    if !std::env::args().any(|arg| arg == "--bench") {
        return;
    }
    let prereg = isomesh::experiment!("P-134");

    common::experiment::run(prereg, |run| {
        let cayley = common::poly::cayley_2x2x2();
        let abs_poly = abs_cayley(&cayley);
        let perms = common::poly::octahedral_relabellings();
        let segs = segments();
        let transforms = DYADIC_SCALES.len() + ROUGH_SCALES.len() + perms.len();

        // Every rounding bound in this file is written in units of Cayley's l1
        // norm, and the denominator is its total degree, so both are read off the
        // polynomial rather than quoted from the header.
        let l1: i128 = cayley.monomials().map(|(_, c)| c.abs()).sum();
        assert_eq!(
            l1, CAYLEY_L1,
            "VOID: Cayley's coefficients have l1 norm {l1}, and every rounding bound in this \
             file is written in units of {CAYLEY_L1}, so the envelope no longer bounds anything"
        );
        assert_eq!(
            cayley.terms(),
            12,
            "VOID: `cayley_2x2x2` has {} terms and not the twelve P-127 recorded, so the \
             polynomial this row normalises is not the one the crate ships",
            cayley.terms()
        );
        assert_eq!(
            cayley.total_degree(),
            4,
            "VOID: `cayley_2x2x2` has total degree {} and not 4, so `(max|f_i|)^4` is not the \
             homogeneous denominator and the normalisation is not scale-free by construction",
            cayley.total_degree()
        );

        let f32_floor = TWO_ROUTE_UNITS * UNIT_ROUNDOFF_F32 * CAYLEY_L1 as f64;
        println!(
            "  Cayley: {} terms, total degree {}, l1 {l1}; f32 sign floor {f32_floor:.4e} \
             -> headline threshold {HEADLINE_THRESHOLD:e}; {transforms} transforms per cell",
            cayley.terms(),
            cayley.total_degree()
        );
        println!(
            "  P-127 measured f32_sign_disagreements = 14 of 3481 trials, all fourteen in the \
             near-zero stratum, so this row is not retired by its own falsifier.\n"
        );

        let mut rows: Vec<Row> = Vec::with_capacity(FIELDS * RESOLUTIONS.len());
        for samples in RESOLUTIONS {
            // Inline block per field, not a closure, so no `return` in here
            // (M-253).
            isomesh::for_each_reference_field!(f64, |name, field| {
                let row = measure(name, &field, samples, &cayley, &abs_poly, &perms, &segs);
                println!(
                    "  {:<14} {:>4}^3 {:<14} surface {:>6} paired {:>6}  m~ {:>10.3e} \
                     above {:>6}  d/r/o {:>8.1e}/{:>8.1e}/{:>8.1e}  units {:>7.4}  \
                     H~ {:>9.3e}  rho_H {:>7.4} rho_D {:>7.4} rho_A {:>7.4}  B>A {:>6} \
                     def {:>5}",
                    row.field,
                    row.samples,
                    row.bound,
                    row.surface_cells,
                    row.paired_cells,
                    row.magnitude_median,
                    row.cells_above_threshold,
                    row.scale_error_dyadic,
                    row.scale_error_rough,
                    row.scale_error_octahedral,
                    row.scale_error_rounding_units,
                    row.hausdorff_median,
                    row.rho_hausdorff,
                    row.rho_defects,
                    row.rho_ambiguous,
                    row.cells_where_b_dominates,
                    row.defect_total,
                );
                rows.push(row);
            });
        }

        // ── vacuity controls, all of them, before the first record ──

        assert_eq!(
            rows.len(),
            FIELDS * RESOLUTIONS.len(),
            "VOID: expected {} rows ({FIELDS} fields x {} resolutions) and built {}, so some \
             field or resolution was never measured and its clauses were never at risk",
            FIELDS * RESOLUTIONS.len(),
            RESOLUTIONS.len(),
            rows.len()
        );

        for row in &rows {
            assert!(
                row.paired_cells >= 3,
                "VOID: {} at {}^3 pairs {} cells, and `common::beta::rank_correlation` defines a \
                 coefficient over fewer than three pairs to be 0.0 rather than measuring one",
                row.field,
                row.samples,
                row.paired_cells
            );
            // The registration's own vacuity control, verbatim: "the per-cell
            // Hausdorff column must have non-zero variance on every field, or the
            // correlation is against a constant".
            assert!(
                row.hausdorff_variance > 0.0,
                "VOID: {} at {}^3 has a per-cell symmetric Hausdorff column of zero variance \
                 over {} cells, so C2's rank correlation is taken against a constant and its \
                 value is an accident of tie-breaking rather than a measurement",
                row.field,
                row.samples,
                row.paired_cells
            );
            assert!(
                row.magnitude_variance > 0.0,
                "VOID: {} at {}^3 has a normalised-magnitude column of zero variance over {} \
                 surface cells, so every correlation in this row is against a constant on the \
                 other side of the pairing",
                row.field,
                row.samples,
                row.surface_cells
            );
            assert!(
                row.min_denominator > 0.0,
                "VOID: {} at {}^3 has a surface cell whose eight corners are all zero, so \
                 `(max|f_i|)^4` is a zero denominator and `delta_magnitude` is not a number",
                row.field,
                row.samples
            );
            assert_eq!(
                row.cells_without_true_surface_points, 0,
                "VOID: {} at {}^3 has {} surface cells where no segment straddled, which is \
                 impossible on a connected cube graph — direction B contributed nothing there \
                 and the symmetric Hausdorff is one-sided",
                row.field, row.samples, row.cells_without_true_surface_points
            );
        }

        // Global, not per row: direction B losing on a given field is that
        // field's geometry, but never winning anywhere means the second mechanism
        // was never at risk.
        let b_dominates_total: usize = rows.iter().map(|r| r.cells_where_b_dominates).sum();
        assert!(
            b_dominates_total > 0,
            "VOID: direction B (bisected true-surface points to the mesh) never exceeds \
             direction A on any cell of any field, so the symmetric Hausdorff is direction A \
             wearing a longer name and the second mechanism was never measured"
        );

        let moved_total: usize = rows.iter().map(|r| r.nonzero_transform_deviations).sum();
        assert!(
            moved_total > 0,
            "VOID: not one of the {transforms} transforms moved the normalised magnitude on any \
             cell of any field, so C1's invariance held without the arithmetic ever being \
             exercised and every zero is a zero that could not have been non-zero (M-44)"
        );

        // Global for the same reason: a clean mesh is the answer on most fields,
        // and `fbm_terrain` is open in its domain so the column is live somewhere
        // — which is what licenses reading a 0.0 elsewhere as that row's geometry
        // rather than a dead instrument.
        let defect_cells_total: usize = rows.iter().map(|r| r.defect_cells).sum();
        assert!(
            defect_cells_total > 0,
            "VOID: not one cell in the whole run carries a validity feature, so every \
             `rank_correlation_with_defects` is against a constant column and none of them is a \
             measurement (M-44)"
        );

        // ── C2's global verdict ──
        let c2_fields_above_bar = rows
            .iter()
            .filter(|r| r.samples == C2_RESOLUTION && r.rho_hausdorff > C2_BAR)
            .count();
        let c2_holds = c2_fields_above_bar >= C2_MIN_FIELDS;
        println!(
            "\n  C2 at {C2_RESOLUTION}^3: {c2_fields_above_bar} of {FIELDS} fields above \
             {C2_BAR} — needs {C2_MIN_FIELDS}, so C2 {}",
            if c2_holds { "HOLDS" } else { "is FALSIFIED" }
        );

        for row in &rows {
            let dyadic_is_exact = row.scale_error_dyadic == 0.0;
            let within_envelope = row.scale_error_rounding_units <= 1.0;
            let c1_holds = dyadic_is_exact && within_envelope;
            let sweep = SWEEP
                .iter()
                .zip(row.sweep_counts.iter())
                .map(|(t, n)| format!("{t:e}:{n}"))
                .collect::<Vec<_>>()
                .join("|");

            run.record(&[
                ("field", row.field.to_string()),
                ("resolution", row.samples.to_string()),
                (
                    "normalisation",
                    "abs(Delta)/max(abs(f_i))^4".to_string(),
                ),
                ("delta_magnitude", format!("{:.9e}", row.magnitude_median)),
                (
                    "scale_invariance_error",
                    format!("{:.9e}", row.scale_error_worst),
                ),
                (
                    "rank_correlation_with_defects",
                    format!("{:.6}", row.rho_defects),
                ),
                (
                    "rank_correlation_with_hausdorff",
                    format!("{:.6}", row.rho_hausdorff),
                ),
                ("threshold_sweep", sweep),
                (
                    "cells_above_threshold",
                    row.cells_above_threshold.to_string(),
                ),
                ("c1_holds", c1_holds.to_string()),
                ("c2_holds", c2_holds.to_string()),
                // ── extras (M-273) ──
                ("ambiguous_face_cells", row.ambiguous_face_cells.to_string()),
                ("bisections_per_segment", BISECTIONS.to_string()),
                ("boundary_features", row.boundary_features.to_string()),
                ("c1_dyadic_is_exact", dyadic_is_exact.to_string()),
                ("c1_within_rounding_envelope", within_envelope.to_string()),
                ("c2_bar", format!("{C2_BAR:.2}")),
                ("c2_fields_above_bar", c2_fields_above_bar.to_string()),
                ("c2_min_fields", C2_MIN_FIELDS.to_string()),
                ("c2_resolution", C2_RESOLUTION.to_string()),
                ("cayley_l1_norm", l1.to_string()),
                ("cayley_terms", cayley.terms().to_string()),
                ("cayley_total_degree", cayley.total_degree().to_string()),
                ("cell_size", format!("{:.9}", row.cell_size)),
                ("cells", row.cells.to_string()),
                (
                    "cells_where_b_dominates",
                    row.cells_where_b_dominates.to_string(),
                ),
                (
                    "cells_without_triangles",
                    row.cells_without_triangles.to_string(),
                ),
                (
                    "cells_without_true_surface_points",
                    row.cells_without_true_surface_points.to_string(),
                ),
                ("defect_cells", row.defect_cells.to_string()),
                ("defect_features", row.defect_total.to_string()),
                ("defect_variance", format!("{:.9e}", row.defect_variance)),
                ("delta_magnitude_max", format!("{:.9e}", row.magnitude_max)),
                ("delta_magnitude_mean", format!("{:.9e}", row.magnitude_mean)),
                ("delta_magnitude_min", format!("{:.9e}", row.magnitude_min)),
                (
                    "delta_magnitude_variance",
                    format!("{:.9e}", row.magnitude_variance),
                ),
                ("distance_reading_is_exact", row.bound_is_exact.to_string()),
                (
                    "f32_disagreements_above_threshold",
                    row.f32_disagreements_above_threshold.to_string(),
                ),
                ("f32_rounding_floor", format!("{f32_floor:.6e}")),
                (
                    "f32_sign_disagreements",
                    row.f32_sign_disagreements.to_string(),
                ),
                ("field_bound", row.bound.to_string()),
                (
                    "gradient_floor_points",
                    row.gradient_floor_points.to_string(),
                ),
                (
                    "gradient_norm_median",
                    format!("{:.9}", row.gradient_norm_median),
                ),
                (
                    "hausdorff_field_to_mesh_median",
                    format!("{:.9e}", row.hausdorff_field_to_mesh_median),
                ),
                ("hausdorff_max", format!("{:.9e}", row.hausdorff_max)),
                ("hausdorff_median", format!("{:.9e}", row.hausdorff_median)),
                (
                    "hausdorff_median_in_cells",
                    format!("{:.6}", row.hausdorff_median / row.cell_size),
                ),
                (
                    "hausdorff_mesh_to_field_median",
                    format!("{:.9e}", row.hausdorff_mesh_to_field_median),
                ),
                (
                    "hausdorff_variance",
                    format!("{:.9e}", row.hausdorff_variance),
                ),
                ("headline_threshold", format!("{HEADLINE_THRESHOLD:e}")),
                ("min_denominator", format!("{:.9e}", row.min_denominator)),
                ("misoriented_features", row.misoriented_features.to_string()),
                ("nonmanifold_features", row.nonmanifold_features.to_string()),
                (
                    "nonzero_transform_deviations",
                    row.nonzero_transform_deviations.to_string(),
                ),
                (
                    "normalised_magnitude_flops",
                    NORMALISED_MAGNITUDE_FLOPS.to_string(),
                ),
                ("octahedral_elements", perms.len().to_string()),
                ("paired_cells", row.paired_cells.to_string()),
                (
                    "rank_correlation_with_ambiguous_faces",
                    format!("{:.6}", row.rho_ambiguous),
                ),
                (
                    "rank_correlation_with_uncertified",
                    format!("{:.6}", row.rho_uncertified),
                ),
                ("rounding_envelope", format!("{ROUNDING_ENVELOPE:.0}")),
                (
                    "scale_error_dyadic",
                    format!("{:.9e}", row.scale_error_dyadic),
                ),
                (
                    "scale_error_octahedral",
                    format!("{:.9e}", row.scale_error_octahedral),
                ),
                ("scale_error_rough", format!("{:.9e}", row.scale_error_rough)),
                (
                    "scale_error_rounding_units",
                    format!("{:.6}", row.scale_error_rounding_units),
                ),
                ("segments_per_cell", SEGMENT_COUNT.to_string()),
                (
                    "shipped_max_normalised_gap",
                    format!("{:.6e}", row.shipped_max_normalised_gap),
                ),
                (
                    "shipped_route_sign_disagreements",
                    row.shipped_sign_disagreements.to_string(),
                ),
                ("surface_cells", row.surface_cells.to_string()),
                ("transforms_per_cell", transforms.to_string()),
                ("triangles", row.triangles.to_string()),
                ("uncertified_cells", row.uncertified_cells.to_string()),
                ("vertices", row.vertices.to_string()),
            ]);
        }
    });
}
