//! **P-131 — `Delta = 0` is not the degenerate case, and the branch that assumes
//! it is has never been audited.**
//!
//! Ticket: R-131. Pre-registered before this harness existed.
//!
//! ```bash
//! cargo bench --bench experiment_p131
//! ```
//!
//! Writes `docs/experiments/p-131.csv`.
//!
//! # What was missing
//!
//! `crates/isomesh/src/marching_cubes/trilinear.rs:250-257` is four lines and a
//! comment:
//!
//! ```text
//! if discriminant == R::ZERO {
//!     // A double root is **one** intersection point, not two: the two
//!     // hyperbolas touch rather than cross, and Proposition 1 counts points.
//!     // Reporting two here would let a degenerate zero-area "hexagon" claim
//!     // six saddles.
//!     roots[0] = -b / (R::TWO * a);
//!     return (roots, 1);
//! }
//! ```
//!
//! Every prior id around it is about the **solver** and none is about the
//! **configuration**:
//!
//! - `M-207` defends the sibling branch at `:238-244` — `a == 0` is a root
//!   *count* and not an absence, because the textbook formula divides by `2a`.
//!   That is a finding about arithmetic.
//! - `M-206` records that `interior::SweptFaces` and `trilinear::BodySaddles`
//!   locate the same body saddles to `1.1e-12` from disjoint parametrisations.
//!   P-131's committed sibling `P-127` explains why: `docs/experiments/p-127.csv`
//!   reads `symbolic_difference_is_zero=true`, `terms_disc=12`, `terms_cayley=12`,
//!   `total_degree=4`, `pencil_matches=3` — `b*b - 4*a*c` at `trilinear.rs:246`
//!   **is** Cayley's `2x2x2` hyperdeterminant of the eight corner values, on the
//!   nose and with no sign flip.
//! - `M-214`, `M-215`, `M-216`, `M-217`, `M-221`, `M-228`, `M-229`, `M-230`,
//!   `✗43` and `✗50` all build on the saddle *count*. Not one asks what a cell
//!   **on** the hyperdeterminant hypersurface is.
//!
//! That last question is the row. Because of `P-127`, `discriminant == 0` is
//! exactly `Det_2,2,2(A) = 0`, so de Silva & Lim's classification applies
//! verbatim rather than by analogy. Their §6 gives the two open strata — *"the
//! rank of a tensor is 2 on the set `{A | Det > 0}` and 3 on the set
//! `{A | Det < 0}`"* — and their Prop 7.3 gives the wall between them: a tensor
//! on `Det = 0` is **generically rank 3 with border rank 2**, the W-state being
//! the canonical witness. A rank-3 tensor that is a limit of rank-2 tensors is
//! not two branches touching. It is a Jordan block in the pencil, and the
//! comment's phrase *"the two hyperbolas touch rather than cross"* is a
//! description of the wrong object.
//!
//! Nothing in the repository has ever counted how often that branch runs.
//!
//! # `roots()` is private, so it is mirrored — and the mirror is checked
//!
//! `BodySaddles::roots` (`trilinear.rs:236-267`) is private. This bench
//! transcribes it, `level_crossing` (`:275-280`), the inside-mask loop
//! (`:180-187`), `inner_hexagon` (`:338-351`), the `is_some()` decision of
//! `interior_vertex` (`:815-897`), `detached_ring` (`:1072-1105`) and
//! `fan_tunnel`'s emission count (`:1120-1259`). **Nothing here is trusted
//! because it was copied carefully.** On every cell of every arm the mirror is
//! compared against the shipped public surface:
//!
//! | mirrored | checked against | column |
//! |---|---|---|
//! | the three coordinate pairs | `BodySaddles::axis(0..3)`, **bitwise** | `mirror_u_disagreements`, `mirror_coordinate_disagreements` |
//! | the inside mask | `BodySaddles::inside_mask()` | `mirror_mask_disagreements` |
//! | `interior_vertex().is_some()` | the shipped `Option` | `mirror_interior_disagreements` |
//! | the six hexagon vertices | `BodySaddles::inner_hexagon()`, **bitwise** | `mirror_hexagon_disagreements` |
//! | `fan_tunnel`'s triangle and unresolved counts | the shipped `fan_tunnel` | `fan_tunnel_disagreements` |
//!
//! The comparison is on raw bits (`f64::to_bits`) and not on `==`, because
//! `level_crossing` divides by a difference that can vanish and `NaN == NaN` is
//! false — an `==` comparison would silently pass over exactly the cells this row
//! is about. All five counters are asserted zero before any row is written; a
//! disagreeing mirror makes every number in the file a number about a different
//! function.
//!
//! # The three branch counters, and why they are three and not one
//!
//! `roots()` has four exits and they are not nested the way the clause names
//! read, so the registered columns are pinned to exits rather than to
//! predicates:
//!
//! - **`a_zero_hits`** — the `a == 0` test at `:238` wins first, so a cell with
//!   `a == 0` never reaches the discriminant test at all. Split further in the
//!   extras as `a_zero_no_root_hits` (`b == 0` too, no root) and
//!   `a_zero_linear_hits` (one root, `-c/b`).
//! - **`discriminant_zero_hits`** — the branch at `:250` **firing**: `a != 0`
//!   and `b*b - 2*2*a*c == 0.0`. This is C1's number and it is per field.
//! - **`double_root_hits`** — the subset of those whose reported double root
//!   `-b/(2*a)` lies **strictly** inside `(0, 1)`, so it sets a mask bit and can
//!   reach a triangle. This is the only subset C3 can move, which discharges the
//!   registration's `SHARE` line.
//!
//! The extras carry `delta_zero_cells`, the count of cells whose discriminant is
//! zero *whatever branch runs*. The identity
//! `delta_zero_cells == discriminant_zero_hits + a_zero_no_root_hits` is asserted
//! per arm: with `a == 0` the discriminant is `b*b`, which is zero exactly when
//! `b` is, so the two exits partition the hypersurface and nothing falls between
//! them.
//!
//! # Arms
//!
//! | arm | what it varies | is_control |
//! |---|---|---|
//! | eight reference fields x `{33, 65, 129}` samples | the field and the grid; 24 rows | no |
//! | `synthetic_w_state` | the canonical W-state cell, one cell | **yes**, the registered vacuity fixture |
//! | `synthetic_w_orbit` | the `GL(2, Z)^3` orbit of that cell | no |
//! | `synthetic_rank_one` | product tensors `x_i*y_j*z_k` | **yes**, the `Det = 0` stratum of rank 1 |
//! | `synthetic_biseparable` | `x_i*N[j][k]` and its two relabellings | **yes**, the `Det = 0` stratum of rank 2 |
//! | `synthetic_generic` | uniform integer tensors | **yes**, the arm that must mostly **not** hit |
//! | hexagon search | random `f64` cells with six inside saddles | **yes**, an assert rather than a row |
//!
//! Three resolutions and these three: `33` and `65` are the house pair, and
//! `129` is added rather than substituted because C1's falsifier is *"the branch
//! never firing on any field at any resolution"* and an exact float coincidence
//! is hunted with cells, of which `129^3` has sixty-four times as many as `33^3`.
//! `65` and `129` are also the `u64`-word-boundary sizes gotcha 4 asks for.
//!
//! The hexagon search is a control and not a row on purpose: its only job is to
//! guarantee that `fan_tunnel_cross_checks` is non-zero, so the `fan_tunnel`
//! mirror C3 leans on has been compared against the shipped function on real
//! six-saddle cells. A control that produced a CSV row would invite that row
//! being quoted as a measurement of a field.
//!
//! # What "correcting the classification" is taken to mean, and why
//!
//! There are exactly three models of a double root and only two of them are
//! arithmetic:
//!
//! - **M1, shipped.** One root. `roots = 1`, so only the `k = 0` slots are
//!   filled and at most three of the six mask bits can be set — an inner hexagon
//!   is unreachable by construction (`:180-187` with `:306-308`).
//! - **M2, multiplicity.** The double root reported twice: `roots = 2` with
//!   `u[1] == u[0]`, hence `v[1] == v[0]` and `w[1] == w[0]` because `:174-177`
//!   derives both from `u`. This is the model the comment rejects, and this
//!   bench prices that rejection instead of asserting it.
//! - **M3, refusal.** The critical point is degenerate, Grosso's Proposition 1
//!   does not apply, and the cell has no rule — the position `SeparateDisks`
//!   already occupies at `mod.rs:429-435`. M3 is a refusal and not a count, so
//!   it cannot appear in `mesh_delta_triangles`; it is named here so nobody
//!   reads its absence as an oversight.
//!
//! So `roots_true` is the exit count with **algebraic multiplicity** (M2) and
//! `roots_reported` is the shipped count (M1); their difference is
//! `discriminant_zero_hits` by construction and that identity is asserted rather
//! than assumed. `mesh_delta_triangles` is `M2 - M1` in triangles, summed over
//! the hit cells that `mod.rs:292-305` would actually route through the
//! trilinear path (`AMBIGUOUS_FACES[case] != 0` and a non-empty contour set).
//!
//! **And the delta is provably zero unless all three coordinates are inside.**
//! Under M2 the mask bits arrive in pairs, so `interior_vertex`'s three line
//! counts are each `0` or `2` and its `total` is `0`, `2` or `6`; `2` is always
//! equal to one of them and returns `None`, and `6` is the hexagon. Under M1
//! only the `k = 0` bits exist, so `fc_u` is identically zero and `total` is at
//! most `2`, reaching `Some` exactly when `u`, `v` and `w` are all inside. Both
//! models therefore differ **only** on cells where the double root and its two
//! crossings all lie in `(0, 1)` — column `hits_all_coords_inside`. Both counts
//! are computed anyway, on one path; the derivation is why the number is what it
//! is, not a shortcut around computing it.
//!
//! # SHARE, recomputed before the numbers
//!
//! *"C3 moves only the cells C1 counts, which is why C1 is reported first and
//! separately."* Discharged structurally: `mesh_delta_triangles` is accumulated
//! inside the `Branch::DiscZero` arm and nowhere else, so a cell that is not one
//! of `discriminant_zero_hits` cannot contribute to it. `hits_all_coords_inside`
//! narrows that further, and `double_root_hits` is the registered column that
//! carries the narrowing.
//!
//! # Exact arithmetic, and where the claim of exactness stops
//!
//! The `exact_arithmetic` column is the disclosure and it is `false` on all
//! twenty-four reference-field rows.
//!
//! - **Synthetic arms are exact.** Every corner value is a small integer with
//!   `|f_i| <= 32` (asserted, `entry_bound_violations`), so `a`, `b`, `c`,
//!   `b*b`, `4*a*c` and every `2x2` minor are integers below `2^28` and are
//!   represented in `f64` without rounding. Independently, Cayley's form is
//!   evaluated over `i128` by `common::poly::cayley_2x2x2().eval_i128` and by
//!   `common::poly::repo_discriminant().eval_i128`, and the sign of each is
//!   compared with the `f64` discriminant on every cell
//!   (`exact_sign_disagreements`, asserted zero).
//! - **Reference-field arms are `f64`.** Corner values are arbitrary doubles, so
//!   Cayley's degree-4 form in them needs about 212 bits and `i128` cannot hold
//!   it. The hit test is therefore the crate's own `f64` expression, which is
//!   correct: C1 asks whether *the shipped branch* fires *at `f64`*, not whether
//!   an exact hypersurface is met. What is still exact on those rows is the
//!   classification's linear algebra: every `2x2` minor is evaluated by Kahan's
//!   error-free determinant `a*d - b*c` through `f64::mul_add`, one rounding of
//!   the exact value, so each minor's **sign** is exact. `f64::mul_add` is
//!   reached directly and not through `Real`, which withholds it on purpose
//!   (`real.rs:44-48`); a true fused multiply-add is the whole mechanism.
//!
//! To keep a zero from being a zero that could not have been non-zero (M-44),
//! every arm also reports `closest_nonzero_disc_rel` — the smallest
//! `|disc| / max(|b*b|, |4*a*c|)` over cells with `a != 0` — and
//! `disc_within_one_ulp`, how many cells came within one `f64` epsilon of the
//! branch without taking it. A reference field with `discriminant_zero_hits = 0`
//! and `disc_within_one_ulp = 0` is a different fact from one with
//! `disc_within_one_ulp = 900`.
//!
//! # Rank 1 is arithmetically unreachable on this branch, and that is recorded
//!
//! C2 asks for a classification into rank 1, 2 or 3. Rank 1 cannot occur:
//!
//! - A rank-1 tensor is `f[i + 2j + 4k] = x_i*y_j*z_k`, so every flattening has
//!   rank 1; in particular the `w` flattening does, i.e. `f4..f7 = t*(f0..f3)`.
//! - Then `twist_hi = t*twist_lo` and `du_hi = t*du_lo`, so
//!   `a = du_hi*twist_lo - du_lo*twist_hi = t*du_lo*twist_lo - du_lo*t*twist_lo`
//!   is **identically zero** — and `:238` wins before `:250` is reached.
//!
//! So the branch's rank population is `{2, 3}`, `branch_rank_one_hits` is
//! predicted `0`, and the column is recorded rather than the clause quietly
//! narrowed (the rule P-70's C1 set). `rank_one_reachable_on_branch` carries the
//! verdict. The same algebra explains a sharper asymmetry the extras carry:
//! `a` is built out of `u`-direction differences, so of the three biseparable
//! strata only the one whose **`u`** flattening is rank 1 reaches the branch —
//! there `a = -(x0 - x1)^2 * det N`, non-zero exactly when the tensor is rank 2.
//! `bisep_u_branch_hits`, `bisep_v_branch_hits` and `bisep_w_branch_hits` are the
//! measurement.
//!
//! # The classifier is two independent routes and their agreement is the check
//!
//! On `Det = 0`:
//!
//! - **Route A, flattenings.** `local_rank(axis)` is the rank of the `2x4`
//!   flattening from its six exact-sign minors. All three `1` is rank 1; exactly
//!   one `1` is rank 2 (biseparable, the remaining `2x2` block being invertible);
//!   no `1` is Prop 7.3's rank 3 with border rank 2. `border_rank_two` and
//!   `true_rank_three` come from here.
//! - **Route B, pencils.** For each of the three axis pairings
//!   (`0123|4567`, `0145|2367`, `0246|1357`, `common::poly`'s own convention)
//!   the pencil `det(A0 + lambda*A1)` is regular when `(c0, c1, c2)` is not the
//!   zero triple. A rank-1 flattening in either of the *other* two axes forces
//!   all three coefficients to vanish, so the number of regular pairings is `0`
//!   for rank `<= 1`, `1` for a rank-2 biseparable and `3` for a W-state.
//!   `w_state_like` comes from here.
//!
//! The two routes share no arithmetic. `class_disagreements` is asserted zero,
//! `local_rank_anomalies` counts the pattern "two flattenings of rank 1 and one
//! of rank 2" that cannot exist, and `regular_anomalies` counts a regular-pairing
//! count of `2`, which likewise cannot. On the `Det = 0` branch `w_state_like`
//! and `true_rank_three` must come out equal — that equality *is* the check, and
//! `delta_negative_cells` in the extras is what keeps the distinction visible: a
//! `Det < 0` cell is rank 3 and is **not** W-like, and there are plenty of them.
//!
//! # Vacuity controls
//!
//! Each runs before the first `run.record` and each panic begins `VOID: `.
//!
//! - **The registered one.** A synthetic W-state cell is in the fixture and must
//!   reach the branch. `f = [0, 1, 1, 0, 1, 0, 0, 0]` under `f[u + 2v + 4w]` is
//!   `|001> + |010> + |100>`; it gives `a = 1`, `b = -2`, `c = 1`,
//!   `disc = 4 - 4 = 0`, so `a != 0` and the branch at `:250` fires.
//!   Columns: `w_state_branch`, `w_state_a/b/c`, `w_state_disc_f64`.
//! - **Its `Det` is exactly zero, twice over.** `w_state_delta_i128` and
//!   `w_state_repo_delta_i128` are both `0` in `i128`.
//! - **It is really rank 3.** `w_state_rank = 3`, `w_state_border_rank = 2`,
//!   `w_state_local_ranks = 2|2|2`, `w_state_regular_pairings = 3`. Without this
//!   C2 could not distinguish "no such cells exist" from "the fixture has none".
//! - **The identity this row's whole framing rests on.**
//!   `repo_discriminant().sub(&cayley_2x2x2()).is_zero()` — if that were false,
//!   `discriminant == 0` would not be `Det = 0` and de Silva & Lim would not
//!   apply. Column `symbolic_identity_holds`.
//! - **The degenerate strata are populated**, so `border_rank_two` and
//!   `true_rank_three` are not one column wearing two names.
//! - **The generic arm mostly misses.** If every cell were a hit the branch test
//!   would be a tautology.
//! - **The mirror is the shipped function**, on every cell of every arm.
//! - **The `fan_tunnel` mirror was exercised** on real six-saddle cells.
//!
//! # Timing
//!
//! No clause here is a cost clause, so nothing is a wall clock. `wall_ns` is
//! recorded beside the verdicts and read by nothing (P-126's rule). Every count
//! is an integer over a deterministic fixture: the reference grids are fixed, the
//! `GL(2, Z)^3` orbit is an exhaustive enumeration rather than a sample, and the
//! three sampled strata use `common::poly::Rng`, a SplitMix64 seeded from the
//! `seed` column.

#![allow(
    clippy::cast_possible_truncation,
    clippy::cast_precision_loss,
    clippy::float_cmp,
    clippy::too_many_lines
)]

mod common;

use std::collections::BTreeSet;
use std::time::Instant;

use common::poly;
use isomesh::fields::ReferenceField;
use isomesh::marching_cubes::ambiguity::joined_mask;
use isomesh::marching_cubes::table::{AMBIGUOUS_FACES, is_inside};
use isomesh::marching_cubes::trilinear::{
    BodySaddles, Contours, MAX_PATCH_TRIANGLES, MAX_TUNNEL_CONTOUR, SADDLE_COUNT, local_crossing,
};
use isomesh::{Sdf, Shape3, for_each_reference_field};

/// The three grids. See the header for why `129` is added and not substituted.
const RESOLUTIONS: [u32; 3] = [33, 65, 129];

/// The seed for the three sampled synthetic strata. SplitMix64, so the same
/// stream on every host and every re-run.
const SEED: u64 = 0x0000_0000_0000_0131;

/// Cells drawn for each sampled synthetic stratum.
const SAMPLED_CELLS: usize = 2048;

/// The magnitude every synthetic corner value must respect, so that `a`, `b`,
/// `c`, `b*b`, `4*a*c` and every `2x2` minor are exact in `f64`. With `32` the
/// widest of those is below `2^29`.
const ENTRY_BOUND: i128 = 32;

/// The integer range the sampled strata draw their factors from.
const FACTOR_RANGE: i64 = 3;

/// Six-saddle cells the hexagon search must find before the `fan_tunnel` mirror
/// is considered exercised.
const HEXAGON_TARGET: usize = 64;

/// Draws the hexagon search is allowed.
const HEXAGON_ATTEMPTS: usize = 4_000_000;

/// Which exit of `roots` a cell took.
///
/// The names are the exits at `trilinear.rs:238-266`, in source order, because a
/// count pinned to a predicate rather than to an exit is a count about a
/// different function: `a == 0` is tested first and wins.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Branch {
    /// `a == 0` and `b == 0` — no root (`:239-241`).
    AZeroNoRoot,
    /// `a == 0`, `b != 0` — the linear root `-c/b` (`:242-243`).
    AZeroLinear,
    /// `a != 0`, `discriminant < 0` — no root (`:247-249`).
    DiscNegative,
    /// `a != 0`, `discriminant == 0` — **the branch this row audits**
    /// (`:250-257`).
    DiscZero,
    /// `a != 0`, `discriminant > 0` — Kahan's two roots (`:263-266`).
    TwoRoots,
}

/// The mirror of `BodySaddles::of` and the private `roots` it calls.
///
/// Transcribed from `trilinear.rs:165-190` and `:236-280`. Every arithmetic step
/// is in the crate's own order and association, including
/// `b * b - R::TWO * R::TWO * a * c`, which is `((2.0 * 2.0) * a) * c` and not
/// `4.0 * (a * c)`. The mirror exists only because `roots` is private; it is
/// never believed, it is compared.
#[derive(Clone, Copy, Debug)]
struct Mirror {
    /// The quadratic's coefficients, from the **public** `coefficients`.
    a: f64,
    /// See [`Mirror::a`].
    b: f64,
    /// See [`Mirror::a`].
    c: f64,
    /// `b*b - 2*2*a*c`, computed on every cell whatever exit was taken.
    disc: f64,
    /// Which exit ran.
    branch: Branch,
    /// How many roots that exit reported.
    roots: usize,
    /// The three coordinate pairs, zero-filled above `roots` exactly as the
    /// crate leaves them.
    coordinate: [[f64; SADDLE_COUNT]; 3],
    /// The inside mask, bit `2*axis + k`.
    mask: u8,
}

impl Mirror {
    /// Solve one cell.
    fn of(corner: &[f64; 8]) -> Self {
        let [a, b, c] = BodySaddles::<f64>::coefficients(corner);
        let disc = b * b - 2.0 * 2.0 * a * c;

        let mut u = [0.0_f64; SADDLE_COUNT];
        let (branch, roots) = if a == 0.0 {
            if b == 0.0 {
                (Branch::AZeroNoRoot, 0)
            } else {
                u[0] = -c / b;
                (Branch::AZeroLinear, 1)
            }
        } else if disc < 0.0 {
            (Branch::DiscNegative, 0)
        } else if disc == 0.0 {
            u[0] = -b / (2.0 * a);
            (Branch::DiscZero, 1)
        } else {
            let q = -(b + b.signum() * disc.sqrt()) * 0.5;
            u[0] = q / a;
            u[1] = c / q;
            (Branch::TwoRoots, SADDLE_COUNT)
        };

        let mut v = [0.0_f64; SADDLE_COUNT];
        let mut w = [0.0_f64; SADDLE_COUNT];
        for ((vk, wk), &uk) in v.iter_mut().zip(w.iter_mut()).zip(u.iter()).take(roots) {
            *vk = level_crossing(corner[0], corner[1], corner[2], corner[3], uk);
            *wk = level_crossing(corner[0], corner[1], corner[4], corner[5], uk);
        }
        let coordinate = [u, v, w];
        Self {
            a,
            b,
            c,
            disc,
            branch,
            roots,
            coordinate,
            mask: inside_mask(&coordinate, roots),
        }
    }
}

/// Where the level set crosses the segment from `lo` to `hi` at parameter `u`.
///
/// `trilinear.rs:275-280`, unguarded there and unguarded here: the mask is the
/// authority on which numbers mean anything, and an epsilon would move a saddle
/// rather than reject it.
fn level_crossing(lo0: f64, lo1: f64, hi0: f64, hi1: f64, u: f64) -> f64 {
    let s = 1.0 - u;
    let lo = lo0 * s + lo1 * u;
    let hi = hi0 * s + hi1 * u;
    -lo / (hi - lo)
}

/// Which of the six coordinates lie strictly inside the cell.
///
/// `trilinear.rs:180-187`, with `coordinate_bit(axis, k) = 1 << (2*axis + k)`
/// from `:134-136`. Strict on both ends: a coordinate of exactly `0` or `1` is
/// on a face and is A-002i's configuration rather than a body saddle.
fn inside_mask(coordinate: &[[f64; SADDLE_COUNT]; 3], roots: usize) -> u8 {
    let mut mask = 0u8;
    for (axis, values) in coordinate.iter().enumerate() {
        for (k, &value) in values.iter().enumerate().take(roots) {
            if value > 0.0 && value < 1.0 {
                mask |= 1 << (2 * axis + k);
            }
        }
    }
    mask
}

/// Whether `BodySaddles::interior_vertex` would return `Some`.
///
/// `trilinear.rs:815-897` depends on the mask alone for that decision: `None`
/// when every in-range line lies in one pair of opposite faces (which also
/// covers `total == 0`), `None` again for the five- and six-line cases the
/// hexagon owns, and `Some` for `total` in `2..=4`. The `u` pair's index pairs
/// are crossed, and that is the crate's derivation at `:834` and not a guess.
fn interior_some(mask: u8) -> bool {
    let set = |axis: usize, k: usize| mask & (1 << (2 * axis + k)) != 0;
    let count = |a: (usize, usize), b: (usize, usize)| u8::from(set(a.0, a.1) && set(b.0, b.1));
    let fc_w = count((0, 0), (1, 0)) + count((0, 1), (1, 1));
    let fc_v = count((0, 0), (2, 0)) + count((0, 1), (2, 1));
    let fc_u = count((1, 0), (2, 1)) + count((1, 1), (2, 0));
    let total = fc_w + fc_v + fc_u;
    if total == fc_w || total == fc_v || total == fc_u {
        return false;
    }
    (2..=4).contains(&total)
}

/// The six vertices of the inner hexagon, in the reference implementation's
/// order.
///
/// `trilinear.rs:343-350`. Written as a function of the three coordinate pairs
/// so that M2's degenerate case — all six vertices at one point, because
/// `u[1] == u[0]` — falls out of the same expression rather than out of a second
/// branch.
fn hexagon(coordinate: &[[f64; SADDLE_COUNT]; 3]) -> [[f64; 3]; 6] {
    let [u, v, w] = *coordinate;
    [
        [u[0], v[0], w[0]],
        [u[0], v[0], w[1]],
        [u[1], v[0], w[1]],
        [u[1], v[1], w[1]],
        [u[1], v[1], w[0]],
        [u[0], v[1], w[0]],
    ]
}

/// Squared distance in the cell's local coordinates. `trilinear.rs:1034-1037`.
fn distance_squared(a: [f64; 3], b: [f64; 3]) -> f64 {
    let d = [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
    d[0] * d[0] + d[1] * d[1] + d[2] * d[2]
}

/// Steps between two hexagon indices around the six-ring. `trilinear.rs:1040`.
const fn ring_distance(a: usize, b: usize) -> usize {
    let r = a.abs_diff(b);
    if r > 2 { 6 - r } else { r }
}

/// Which ring, if any, is not part of the tunnel. `trilinear.rs:1072-1105`,
/// including the transcription fix that module documents.
fn detached_ring(u: [f64; SADDLE_COUNT], contours: &Contours, corner: &[f64; 8]) -> Option<usize> {
    if contours.count() != 3 {
        return None;
    }
    let (u_lo, u_hi) = if u[0] < u[1] {
        (u[0], u[1])
    } else {
        (u[1], u[0])
    };
    (0..contours.count()).find(|&r| {
        let ring = contours.ring(r);
        if ring.len() != 3 {
            return false;
        }
        let mut lo = f64::INFINITY;
        let mut hi = f64::NEG_INFINITY;
        for &e in ring {
            let x = local_crossing(e, corner)[0];
            lo = lo.min(x);
            hi = hi.max(x);
        }
        u_lo > hi || u_hi < lo
    })
}

/// How many triangles `Contours::fan_tunnel` emits, and how many ring edges it
/// leaves unresolved.
///
/// `trilinear.rs:1120-1259`, counted rather than emitted. The `spin` pass at
/// `:1227-1254` reorients the closing fan without changing its size, so it is
/// four triangles here as it is there. Cross-checked against the shipped
/// function on every six-saddle cell of every arm.
fn fan_tunnel_count(
    hex: &[[f64; 3]; 6],
    u: [f64; SADDLE_COUNT],
    contours: &Contours,
    corner: &[f64; 8],
) -> (usize, usize) {
    let detached = detached_ring(u, contours, corner);
    let nearest = |edge: u8| -> usize {
        let p = local_crossing(edge, corner);
        let mut best = 0usize;
        let mut best_d = f64::INFINITY;
        for (k, &h) in hex.iter().enumerate() {
            let d = distance_squared(p, h);
            if d < best_d {
                best_d = d;
                best = k;
            }
        }
        best
    };

    let mut n = 0usize;
    let mut unresolved = 0usize;
    for r in 0..contours.count() {
        let ring = contours.ring(r);
        if Some(r) == detached {
            n += 1;
            continue;
        }
        for k in 0..ring.len() {
            let (a, b) = (ring[k], ring[(k + 1) % ring.len()]);
            match ring_distance(nearest(a), nearest(b)) {
                0 => n += 1,
                1 => n += 2,
                2 => n += 3,
                _ => unresolved += 1,
            }
        }
    }
    if contours.count() == 1 {
        n += 4;
    }
    (n, unresolved)
}

/// Kahan's error-free `2x2` determinant `a*d - b*c`.
///
/// `w` is the rounded product, `e` recovers its rounding error exactly and `f`
/// is `a*d - w` exactly; both fused multiply-adds are exact, so `f + e` is a
/// **single** rounding of the exact determinant. A single rounding cannot turn a
/// non-zero value into a zero of the wrong sign, so the sign returned is exact —
/// which is the only property the rank tests read. `f64::mul_add` is used
/// directly rather than through `Real`, which withholds `mul_add` deliberately
/// (`real.rs:44-48`).
fn det2(a: f64, b: f64, c: f64, d: f64) -> f64 {
    let w = b * c;
    let e = (-b).mul_add(c, w);
    let f = a.mul_add(d, -w);
    f + e
}

/// The two rows of the `2x4` flattening along one axis.
///
/// Axis `0` fixes `i` (`u`), axis `1` fixes `j` (`v`), axis `2` fixes `k` (`w`),
/// and the remaining two bits count up — the same reading `common::poly`'s
/// pairings use, so a pairing and its flattening are the same two rows.
fn flattening(f: &[f64; 8], axis: usize) -> ([f64; 4], [f64; 4]) {
    let idx: [[usize; 4]; 2] = match axis {
        0 => [[0, 2, 4, 6], [1, 3, 5, 7]],
        1 => [[0, 1, 4, 5], [2, 3, 6, 7]],
        _ => [[0, 1, 2, 3], [4, 5, 6, 7]],
    };
    (
        std::array::from_fn(|c| f[idx[0][c]]),
        std::array::from_fn(|c| f[idx[1][c]]),
    )
}

/// The rank of one `2x4` flattening: `0`, `1` or `2`, from exact-sign minors.
fn local_rank(f: &[f64; 8], axis: usize) -> u8 {
    let (r, s) = flattening(f, axis);
    if r.iter().chain(s.iter()).all(|&x| x == 0.0) {
        return 0;
    }
    for i in 0..4 {
        for j in i + 1..4 {
            if det2(r[i], r[j], s[i], s[j]) != 0.0 {
                return 2;
            }
        }
    }
    1
}

/// The two opposite-face corner sets of each axis pairing, `common::poly`'s
/// convention: `0` splits along `w`, `1` along `v`, `2` along `u`, each set read
/// row-major into a `2x2`.
const PAIRINGS: [([usize; 4], [usize; 4]); 3] = [
    ([0, 1, 2, 3], [4, 5, 6, 7]),
    ([0, 1, 4, 5], [2, 3, 6, 7]),
    ([0, 2, 4, 6], [1, 3, 5, 7]),
];

/// Is the pencil `det(A0 + lambda*A1)` of this pairing regular — not the zero
/// quadratic?
///
/// `c0 = det A0`, `c2 = det A1` and `c1` is the mixed term
/// `A0[0][0]*A1[1][1] + A1[0][0]*A0[1][1] - A0[0][1]*A1[1][0] - A1[0][1]*A0[1][0]`,
/// written as two exact-sign `2x2` determinants. A rank-1 flattening in either
/// of the other two axes kills all three coefficients, which is what makes the
/// regular-pairing count a rank certificate independent of Route A.
fn pencil_regular(f: &[f64; 8], pairing: usize) -> bool {
    let (lo, hi) = PAIRINGS[pairing];
    let a0 = [f[lo[0]], f[lo[1]], f[lo[2]], f[lo[3]]];
    let a1 = [f[hi[0]], f[hi[1]], f[hi[2]], f[hi[3]]];
    let c0 = det2(a0[0], a0[1], a0[2], a0[3]);
    let c2 = det2(a1[0], a1[1], a1[2], a1[3]);
    let c1 = det2(a0[0], a0[1], a1[2], a1[3]) + det2(a1[0], a1[1], a0[2], a0[3]);
    !(c0 == 0.0 && c1 == 0.0 && c2 == 0.0)
}

/// What one cell on `Det = 0` actually is.
#[derive(Clone, Copy, Debug)]
struct Class {
    /// The three flattening ranks, Route A's input.
    local_ranks: [u8; 3],
    /// Real tensor rank, Route A.
    rank: u8,
    /// Real border rank, Route A.
    border_rank: u8,
    /// How many of the three pairings are regular, Route B's whole content.
    regular_pairings: u8,
    /// Route B's rank, from `regular_pairings` alone.
    rank_b: u8,
    /// Two flattenings of rank 1 and one of rank 2 — impossible, counted.
    local_rank_anomaly: bool,
    /// A regular-pairing count of `2` — impossible, counted.
    regular_anomaly: bool,
}

impl Class {
    /// Classify a cell **on the `Det = 0` hypersurface**.
    ///
    /// Off that hypersurface the rank follows the sign of `Det` (de Silva & Lim
    /// §6) and neither route below applies; this is only ever called from the
    /// `Branch::DiscZero` arm.
    fn of(f: &[f64; 8]) -> Self {
        let local_ranks = [local_rank(f, 0), local_rank(f, 1), local_rank(f, 2)];
        let ones = local_ranks.iter().filter(|&&r| r == 1).count();
        let zeros = local_ranks.iter().filter(|&&r| r == 0).count();

        let mut local_rank_anomaly = false;
        let (rank, border_rank) = if zeros > 0 {
            (0, 0)
        } else {
            match ones {
                3 => (1, 1),
                1 => (2, 2),
                0 => (3, 2),
                _ => {
                    local_rank_anomaly = true;
                    (u8::MAX, u8::MAX)
                }
            }
        };

        let regular_pairings = (0..3).filter(|&p| pencil_regular(f, p)).count() as u8;
        let mut regular_anomaly = false;
        let rank_b = if zeros > 0 {
            0
        } else {
            match regular_pairings {
                0 => 1,
                1 => 2,
                3 => 3,
                _ => {
                    regular_anomaly = true;
                    u8::MAX
                }
            }
        };

        Self {
            local_ranks,
            rank,
            border_rank,
            regular_pairings,
            rank_b,
            local_rank_anomaly,
            regular_anomaly,
        }
    }
}

/// Everything one arm counts.
#[derive(Clone, Debug, Default)]
struct Tally {
    /// Cells swept.
    cells: u64,
    /// `a == 0`, `b == 0`.
    a_zero_no_root: u64,
    /// `a == 0`, `b != 0`.
    a_zero_linear: u64,
    /// `a != 0`, `disc < 0`.
    disc_negative: u64,
    /// `a != 0`, `disc == 0` — C1's number.
    disc_zero: u64,
    /// `a != 0`, `disc > 0`.
    two_roots: u64,
    /// `disc == 0.0` whatever exit ran.
    delta_zero: u64,
    /// `disc < 0.0` whatever exit ran, i.e. the rank-3 stratum that is not W.
    delta_negative: u64,
    /// Hits whose double root is strictly inside `(0, 1)` — `double_root_hits`.
    root_inside: u64,
    /// Hits whose three coordinates are all strictly inside; the only cells the
    /// two models can disagree about.
    all_coords_inside: u64,
    /// Route A's rank census over the hits.
    rank_hits: [u64; 4],
    /// Hits with border rank exactly `2`.
    border_rank_two: u64,
    /// Hits with three regular pairings — Route B's W-state count.
    w_class: u64,
    /// Hits where Route A and Route B disagreed about the rank.
    class_disagreements: u64,
    /// Impossible flattening patterns.
    local_rank_anomalies: u64,
    /// Impossible regular-pairing counts.
    regular_anomalies: u64,
    /// Roots the shipped exits reported, over all cells.
    roots_reported: u64,
    /// Roots with algebraic multiplicity, over all cells.
    roots_true: u64,
    /// Hits that `mod.rs:292-305` would route through the trilinear path.
    trilinear_eligible: u64,
    /// Hits that gain an inner hexagon under M2.
    hexagon_gained: u64,
    /// Hits whose M2 topology is `SeparateDisks`, which has no rule and so no
    /// triangle count.
    separate_disks: u64,
    /// M1's triangles over the eligible hits.
    m1_triangles: u64,
    /// M2's triangles over the same cells.
    m2_triangles: u64,
    /// `m2_triangles - m1_triangles`.
    delta_triangles: i64,
    /// Mirror against `BodySaddles::axis(0)`, bitwise.
    mirror_u: u64,
    /// Mirror against `BodySaddles::axis(1..3)`, bitwise.
    mirror_coordinate: u64,
    /// Mirror against `BodySaddles::inside_mask()`.
    mirror_mask: u64,
    /// Mirror against `BodySaddles::interior_vertex().is_some()`.
    mirror_interior: u64,
    /// Mirror against `BodySaddles::inner_hexagon()`, bitwise.
    mirror_hexagon: u64,
    /// Six-saddle cells on which the shipped `fan_tunnel` was called.
    fan_tunnel_checks: u64,
    /// Of those, how many the mirror got wrong.
    fan_tunnel_disagreements: u64,
    /// The widest patch either model asked for.
    max_patch: u64,
    /// The smallest `|disc| / max(|b*b|, |4*a*c|)` over cells with `a != 0` and a
    /// non-zero discriminant.
    closest_nonzero: Option<f64>,
    /// Cells that came within one `f64` epsilon of the branch without taking it.
    within_one_ulp: u64,
    /// Synthetic arms only: `i128` Cayley disagreed in sign with the `f64`
    /// discriminant.
    exact_sign_disagreements: u64,
    /// Synthetic arms only: a corner value outside [`ENTRY_BOUND`].
    entry_bound_violations: u64,
}

impl Tally {
    /// Sweep one cell: mirror it, check the mirror, classify it if it is a hit,
    /// and price both models of its double root.
    fn sweep(&mut self, corner: &[f64; 8]) -> Mirror {
        self.cells += 1;
        let m = Mirror::of(corner);
        let saddles = BodySaddles::of(corner);

        // ── the mirror is compared, never believed ──────────────────────────
        for (axis, mine) in m.coordinate.iter().enumerate() {
            let theirs = saddles.axis(axis);
            for (a, b) in mine.iter().zip(theirs.iter()) {
                if a.to_bits() != b.to_bits() {
                    if axis == 0 {
                        self.mirror_u += 1;
                    } else {
                        self.mirror_coordinate += 1;
                    }
                }
            }
        }
        if m.mask != saddles.inside_mask() {
            self.mirror_mask += 1;
        }
        if interior_some(m.mask) != saddles.interior_vertex().is_some() {
            self.mirror_interior += 1;
        }
        if let Some(theirs) = saddles.inner_hexagon() {
            let mine = hexagon(&m.coordinate);
            for (p, q) in mine.iter().flatten().zip(theirs.iter().flatten()) {
                if p.to_bits() != q.to_bits() {
                    self.mirror_hexagon += 1;
                }
            }
        }

        // ── the branch census ───────────────────────────────────────────────
        match m.branch {
            Branch::AZeroNoRoot => self.a_zero_no_root += 1,
            Branch::AZeroLinear => self.a_zero_linear += 1,
            Branch::DiscNegative => self.disc_negative += 1,
            Branch::DiscZero => self.disc_zero += 1,
            Branch::TwoRoots => self.two_roots += 1,
        }
        if m.disc == 0.0 {
            self.delta_zero += 1;
        }
        if m.disc < 0.0 {
            self.delta_negative += 1;
        }
        self.roots_reported += m.roots as u64;
        self.roots_true += m.roots as u64 + u64::from(m.branch == Branch::DiscZero);

        // ── how close the misses came (M-44) ────────────────────────────────
        if m.a != 0.0 && m.disc != 0.0 {
            let scale = (m.b * m.b).abs().max((2.0 * 2.0 * m.a * m.c).abs());
            if scale > 0.0 && scale.is_finite() {
                let rel = m.disc.abs() / scale;
                self.closest_nonzero = Some(self.closest_nonzero.map_or(rel, |c| c.min(rel)));
                if rel <= f64::EPSILON {
                    self.within_one_ulp += 1;
                }
            }
        }

        // ── the cell's own contour set, for the fan_tunnel cross-check and
        //    for both models' triangle counts ────────────────────────────────
        let mut case = 0u8;
        for (c, &value) in corner.iter().enumerate() {
            if is_inside(value) {
                case |= 1 << c;
            }
        }
        let ambiguous = AMBIGUOUS_FACES[case as usize];
        let eligible = ambiguous != 0;
        let mask = if eligible {
            joined_mask(corner, ambiguous)
        } else {
            0
        };
        let contours = if eligible {
            Contours::of(case, mask)
        } else {
            Contours::of(0, 0)
        };

        // The mirror of `fan_tunnel`, checked against the shipped one wherever
        // the shipped path would actually call it.
        if eligible && contours.count() > 0 && saddles.has_inner_hexagon() {
            self.fan_tunnel_checks += 1;
            let mut emitted = 0usize;
            let unresolved = contours.fan_tunnel(&saddles, corner, |_| emitted += 1);
            let hex = hexagon(&m.coordinate);
            let (mine, mine_unresolved) =
                fan_tunnel_count(&hex, m.coordinate[0], &contours, corner);
            if mine != emitted || mine_unresolved != unresolved {
                self.fan_tunnel_disagreements += 1;
            }
            self.max_patch = self.max_patch.max(emitted as u64);
        }

        if m.branch != Branch::DiscZero {
            return m;
        }

        // ── C2: what the hit actually is ────────────────────────────────────
        if m.coordinate[0][0] > 0.0 && m.coordinate[0][0] < 1.0 {
            self.root_inside += 1;
        }
        let class = Class::of(corner);
        if class.local_rank_anomaly {
            self.local_rank_anomalies += 1;
        }
        if class.regular_anomaly {
            self.regular_anomalies += 1;
        }
        if class.rank != u8::MAX && (class.rank as usize) < self.rank_hits.len() {
            self.rank_hits[class.rank as usize] += 1;
        }
        if class.border_rank == 2 {
            self.border_rank_two += 1;
        }
        if class.regular_pairings == 3 {
            self.w_class += 1;
        }
        if class.rank != class.rank_b {
            self.class_disagreements += 1;
        }

        // ── C3: M1 against M2, on this cell ─────────────────────────────────
        //
        // M2 is the double root reported with its multiplicity: `roots = 2` with
        // `u[1] == u[0]`, so `:174-177` gives `v[1] == v[0]` and `w[1] == w[0]`
        // and the six hexagon vertices coincide. Built through the same
        // `inside_mask` and the same `hexagon` the shipped path uses.
        let corrected = [
            [m.coordinate[0][0]; SADDLE_COUNT],
            [m.coordinate[1][0]; SADDLE_COUNT],
            [m.coordinate[2][0]; SADDLE_COUNT],
        ];
        let corrected_mask = inside_mask(&corrected, SADDLE_COUNT);
        let all_inside = corrected_mask == 0b0011_1111;
        if all_inside {
            self.all_coords_inside += 1;
        }
        if !eligible || contours.count() == 0 {
            return m;
        }
        self.trilinear_eligible += 1;

        let m1 = contours.triangle_count(interior_some(m.mask)) as u64;
        if all_inside {
            self.hexagon_gained += 1;
            // `Contours::topology` at `:562-571`: a hexagon with two or more
            // rings and one longer than Corollary 6's bound is `SeparateDisks`,
            // which `mod.rs:429-435` refuses. A refusal is not a triangle count,
            // so such a cell is reported and left out of the delta.
            if contours.count() >= 2 && contours.longest() > MAX_TUNNEL_CONTOUR {
                self.separate_disks += 1;
                return m;
            }
            let hex = hexagon(&corrected);
            let (m2, unresolved) = fan_tunnel_count(&hex, corrected[0], &contours, corner);
            // A ring edge three steps apart on the hexagon has no rule
            // (`:1211`, `mod.rs:529-535`). On a hexagon whose six vertices
            // coincide every ring distance is zero, so this cannot fire; it is
            // folded into `separate_disks` rather than silently dropped.
            if unresolved != 0 {
                self.separate_disks += 1;
                return m;
            }
            self.max_patch = self.max_patch.max(m2 as u64);
            self.m1_triangles += m1;
            self.m2_triangles += m2 as u64;
            self.delta_triangles += m2 as i64 - m1 as i64;
        } else {
            let m2 = contours.triangle_count(interior_some(corrected_mask)) as u64;
            self.m1_triangles += m1;
            self.m2_triangles += m2;
            self.delta_triangles += m2 as i64 - m1 as i64;
        }
        m
    }

    /// `a_zero_no_root + a_zero_linear` — the `:238` exit, C1's sibling count.
    fn a_zero(&self) -> u64 {
        self.a_zero_no_root + self.a_zero_linear
    }
}

/// One CSV row's worth of work.
#[derive(Clone, Debug)]
struct Arm {
    /// The `field` column.
    field: String,
    /// The `resolution` column: samples per axis. A synthetic cell **is** a
    /// two-samples-per-axis grid with one cell in it, so `2` is literal there
    /// rather than a sentinel.
    resolution: u32,
    /// Whether this arm's discriminants and minors are exact integers.
    exact: bool,
    /// What it counted.
    tally: Tally,
}

/// Sweep one reference field at one resolution.
///
/// The value grid is built with the same expression `crate::sdf::sample_grid`
/// uses (`sdf.rs:183-187`) — `origin + cell_size * n` with `x` innermost — so a
/// corner value here is bit-identical to the one `MarchingCubes::extract` reads.
fn sweep_field<F>(field: &F, samples: u32) -> Tally
where
    F: ReferenceField + Sdf<Scalar = f64>,
{
    let (shape, origin, h) = common::grid::<f64, _>(field, samples);
    let size = shape.size();
    let mut values = Vec::with_capacity(shape.element_count());
    for z in 0..size[2] {
        for y in 0..size[1] {
            for x in 0..size[0] {
                values.push(field.sample([
                    origin[0] + h * f64::from(x),
                    origin[1] + h * f64::from(y),
                    origin[2] + h * f64::from(z),
                ]));
            }
        }
    }

    let mut tally = Tally::default();
    for z in 0..size[2] - 1 {
        for y in 0..size[1] - 1 {
            for x in 0..size[0] - 1 {
                let corner: [f64; 8] = std::array::from_fn(|c| {
                    let i = shape.linearize([
                        x + (c as u32 & 1),
                        y + ((c as u32 >> 1) & 1),
                        z + ((c as u32 >> 2) & 1),
                    ]);
                    values[i as usize]
                });
                tally.sweep(&corner);
            }
        }
    }
    tally
}

/// The canonical real W-state as a `2x2x2` tensor.
///
/// `|001> + |010> + |100>` is corner values `[0, 1, 1, 0, 1, 0, 0, 0]` under
/// `f[u + 2v + 4w]`: the three tensor entries with exactly one index set are
/// `a(1,0,0) = f1`, `a(0,1,0) = f2` and `a(0,0,1) = f4`.
const W_STATE: [i128; 8] = [0, 1, 1, 0, 1, 0, 0, 0];

/// Widen an integer cell to the `f64` the crate's arithmetic actually sees.
fn to_f64(f: &[i128; 8]) -> [f64; 8] {
    std::array::from_fn(|i| f[i] as f64)
}

/// Every `2x2` integer matrix with entries in `{-1, 0, 1}` and determinant of
/// magnitude one.
///
/// Enumerated rather than listed, in a fixed order. `Det` is a relative
/// invariant of `GL(2)^3` with weight `(det g1 * det g2 * det g3)^2`
/// (`common::poly`'s own header), so a unit determinant leaves `Det` **equal**
/// and not merely zero — the orbit of a `Det = 0` tensor stays exactly on the
/// hypersurface with no rounding anywhere.
fn small_gl2() -> Vec<[[i128; 2]; 2]> {
    let mut out = Vec::new();
    for a in -1..=1i128 {
        for b in -1..=1i128 {
            for c in -1..=1i128 {
                for d in -1..=1i128 {
                    if (a * d - b * c).abs() == 1 {
                        out.push([[a, b], [c, d]]);
                    }
                }
            }
        }
    }
    out
}

/// The `GL(2, Z)^3` orbit of the W-state under [`small_gl2`], deduplicated.
///
/// Exhaustive over the triples, so this fixture is a set and not a sample. Every
/// element is rank 3 with border rank 2 — rank and border rank are invariants of
/// the group action — and every element has `Det = 0` exactly.
fn w_orbit(generators: &[[[i128; 2]; 2]]) -> Vec<[i128; 8]> {
    let mut seen = BTreeSet::new();
    for g1 in generators {
        for g2 in generators {
            for g3 in generators {
                seen.insert(poly::act_gl2_cubed(*g1, *g2, *g3, &W_STATE));
            }
        }
    }
    seen.into_iter().collect()
}

/// Product tensors `f[i + 2j + 4k] = x_i * y_j * z_k`.
///
/// The rank-`<= 1` stratum of `Det = 0`. Every flattening has rank one, so `a`
/// is identically zero and none of these can reach the audited branch — which is
/// the point of having them: the arm measures that unreachability rather than
/// assuming it.
fn rank_one_cells(rng: &mut poly::Rng) -> Vec<[i128; 8]> {
    let mut out = Vec::with_capacity(SAMPLED_CELLS);
    for _ in 0..SAMPLED_CELLS {
        let draw = |r: &mut poly::Rng| -> [i128; 2] {
            [
                i128::from(r.next_i64_in(-FACTOR_RANGE, FACTOR_RANGE + 1)),
                i128::from(r.next_i64_in(-FACTOR_RANGE, FACTOR_RANGE + 1)),
            ]
        };
        let x = draw(rng);
        let y = draw(rng);
        let z = draw(rng);
        out.push(std::array::from_fn(|c| {
            x[c & 1] * y[(c >> 1) & 1] * z[(c >> 2) & 1]
        }));
    }
    out
}

/// Biseparable tensors, in all three orientations, tagged by the axis whose
/// flattening is rank one.
///
/// Orientation `0` is `f[i + 2j + 4k] = x_i * N[j][k]`, the one whose `a` is
/// `-(x0 - x1)^2 * det N` and which therefore **does** reach the branch when the
/// tensor is genuinely rank 2. Orientations `1` and `2` are the `v` and `w`
/// versions, whose `a` cancels identically. That asymmetry is a property of the
/// crate's coefficient construction and is measured, not argued.
fn biseparable_cells(rng: &mut poly::Rng) -> Vec<([i128; 8], usize)> {
    let mut out = Vec::with_capacity(SAMPLED_CELLS * 3);
    for _ in 0..SAMPLED_CELLS {
        let mut pair = || {
            [
                i128::from(rng.next_i64_in(-FACTOR_RANGE, FACTOR_RANGE + 1)),
                i128::from(rng.next_i64_in(-FACTOR_RANGE, FACTOR_RANGE + 1)),
            ]
        };
        let x = pair();
        let n = [pair(), pair()];
        // `f = x_i * N[j][k]`, `f = y_j * N[i][k]`, `f = z_k * N[i][j]`.
        out.push((
            std::array::from_fn(|c| x[c & 1] * n[(c >> 1) & 1][(c >> 2) & 1]),
            0usize,
        ));
        out.push((
            std::array::from_fn(|c| x[(c >> 1) & 1] * n[c & 1][(c >> 2) & 1]),
            1usize,
        ));
        out.push((
            std::array::from_fn(|c| x[(c >> 2) & 1] * n[c & 1][(c >> 1) & 1]),
            2usize,
        ));
    }
    out
}

/// Uniform integer tensors — the arm that must mostly **miss**.
fn generic_cells(rng: &mut poly::Rng) -> Vec<[i128; 8]> {
    (0..SAMPLED_CELLS)
        .map(|_| {
            std::array::from_fn(|_| i128::from(rng.next_i64_in(-FACTOR_RANGE, FACTOR_RANGE + 1)))
        })
        .collect()
}

/// Sweep a list of exact integer cells, cross-checking each against `i128`.
fn sweep_exact(cells: &[[i128; 8]], cayley: &poly::Poly, repo: &poly::Poly) -> Tally {
    let mut tally = Tally::default();
    for f in cells {
        exact_cell(f, cayley, repo, &mut tally);
    }
    tally
}

/// One exact cell: the bound, the `i128` cross-check, and the sweep.
fn exact_cell(f: &[i128; 8], cayley: &poly::Poly, repo: &poly::Poly, tally: &mut Tally) -> Mirror {
    if f.iter().any(|v| v.abs() > ENTRY_BOUND) {
        tally.entry_bound_violations += 1;
    }
    let corner = to_f64(f);
    let m = tally.sweep(&corner);
    let exact = cayley.eval_i128(f);
    let exact_repo = repo.eval_i128(f);
    if exact != exact_repo || exact.signum() as i32 != sign_of(m.disc) {
        tally.exact_sign_disagreements += 1;
    }
    m
}

/// The sign of a float as `-1`, `0` or `1`, with `-0.0` reading `0`.
fn sign_of(x: f64) -> i32 {
    if x > 0.0 {
        1
    } else if x < 0.0 {
        -1
    } else {
        0
    }
}

/// Find six-saddle cells so the `fan_tunnel` mirror is exercised.
///
/// Uniform `f64` cells in `[-1, 1]^8` from the shared SplitMix64, accepted when
/// the shipped `BodySaddles` reports an inner hexagon on a cell the trilinear
/// path would actually take. Returns the cells found and the draws it took.
fn hexagon_search(rng: &mut poly::Rng) -> (Vec<[f64; 8]>, usize) {
    let mut found = Vec::with_capacity(HEXAGON_TARGET);
    let mut attempts = 0usize;
    while found.len() < HEXAGON_TARGET && attempts < HEXAGON_ATTEMPTS {
        attempts += 1;
        let corner: [f64; 8] = std::array::from_fn(|_| rng.next_f64_unit() * 2.0 - 1.0);
        if !BodySaddles::of(&corner).has_inner_hexagon() {
            continue;
        }
        let mut case = 0u8;
        for (c, &value) in corner.iter().enumerate() {
            if is_inside(value) {
                case |= 1 << c;
            }
        }
        let ambiguous = AMBIGUOUS_FACES[case as usize];
        if ambiguous == 0 {
            continue;
        }
        if Contours::of(case, joined_mask(&corner, ambiguous)).count() == 0 {
            continue;
        }
        found.push(corner);
    }
    (found, attempts)
}

fn main() {
    if !std::env::args().any(|arg| arg == "--bench") {
        return;
    }
    let prereg = isomesh::experiment!("P-131");

    common::experiment::run(prereg, |run| {
        let started = Instant::now();

        // ── the identity the whole framing rests on ─────────────────────────
        let cayley = poly::cayley_2x2x2();
        let repo = poly::repo_discriminant();
        let symbolic_identity = repo.sub(&cayley).is_zero();

        // ── the synthetic fixture ───────────────────────────────────────────
        let generators = small_gl2();
        let orbit = w_orbit(&generators);
        let mut rng = poly::Rng::new(SEED);
        let rank_one = rank_one_cells(&mut rng);
        let biseparable = biseparable_cells(&mut rng);
        let generic = generic_cells(&mut rng);

        // The registered vacuity fixture, on its own, so its numbers can be
        // read out of the CSV rather than inferred from an aggregate.
        let mut w_tally = Tally::default();
        let w_mirror = exact_cell(&W_STATE, &cayley, &repo, &mut w_tally);
        let w_class = Class::of(&to_f64(&W_STATE));
        let w_delta = cayley.eval_i128(&W_STATE);
        let w_repo_delta = repo.eval_i128(&W_STATE);

        let orbit_tally = sweep_exact(&orbit, &cayley, &repo);
        let rank_one_tally = sweep_exact(&rank_one, &cayley, &repo);
        let generic_tally = sweep_exact(&generic, &cayley, &repo);

        let mut bisep_tally = Tally::default();
        let mut bisep_branch = [0u64; 3];
        for (f, orientation) in &biseparable {
            let m = exact_cell(f, &cayley, &repo, &mut bisep_tally);
            if m.branch == Branch::DiscZero {
                bisep_branch[*orientation] += 1;
            }
        }

        // ── the fan_tunnel mirror's own control ─────────────────────────────
        let (hexagon_cells, hexagon_attempts) = hexagon_search(&mut rng);
        let mut hexagon_tally = Tally::default();
        for corner in &hexagon_cells {
            hexagon_tally.sweep(corner);
        }

        // ── the reference fields ────────────────────────────────────────────
        let mut arms: Vec<Arm> = Vec::new();
        for_each_reference_field!(f64, |name, field| {
            for samples in RESOLUTIONS {
                arms.push(Arm {
                    field: name.to_string(),
                    resolution: samples,
                    exact: false,
                    tally: sweep_field(&field, samples),
                });
            }
        });
        let reference_arms = arms.len();

        arms.push(Arm {
            field: "synthetic_w_state".to_string(),
            resolution: 2,
            exact: true,
            tally: w_tally.clone(),
        });
        arms.push(Arm {
            field: "synthetic_w_orbit".to_string(),
            resolution: 2,
            exact: true,
            tally: orbit_tally,
        });
        arms.push(Arm {
            field: "synthetic_rank_one".to_string(),
            resolution: 2,
            exact: true,
            tally: rank_one_tally,
        });
        arms.push(Arm {
            field: "synthetic_biseparable".to_string(),
            resolution: 2,
            exact: true,
            tally: bisep_tally,
        });
        arms.push(Arm {
            field: "synthetic_generic".to_string(),
            resolution: 2,
            exact: true,
            tally: generic_tally,
        });

        // ── vacuity controls, every one before the first row ────────────────
        assert!(
            symbolic_identity,
            "VOID: repo_discriminant() - cayley_2x2x2() is not identically zero, so \
             `discriminant == 0` at trilinear.rs:250 is not `Det_2,2,2 = 0` and de Silva \
             & Lim's Prop 7.3 does not apply to any cell this bench counts -- every \
             rank column would be a claim about the wrong hypersurface"
        );
        assert_eq!(
            w_mirror.branch,
            Branch::DiscZero,
            "VOID: the synthetic W-state cell {W_STATE:?} did not reach the \
             `discriminant == 0` branch (a={}, b={}, c={}, disc={}, exit={:?}), so C2 \
             cannot distinguish `no such cells exist` from `the fixture has none` -- \
             which is the registration's own vacuity control, verbatim",
            w_mirror.a,
            w_mirror.b,
            w_mirror.c,
            w_mirror.disc,
            w_mirror.branch
        );
        assert!(
            w_delta == 0 && w_repo_delta == 0,
            "VOID: the W-state's hyperdeterminant is not exactly zero in i128 \
             (cayley={w_delta}, repo={w_repo_delta}), so the fixture is not on the \
             hypersurface and the rank-3 claim it is supposed to witness is unwitnessed"
        );
        assert!(
            w_class.rank == 3
                && w_class.border_rank == 2
                && w_class.regular_pairings == 3
                && w_class.local_ranks == [2, 2, 2],
            "VOID: the W-state did not classify as rank 3 with border rank 2 \
             (rank={}, border_rank={}, regular_pairings={}, local_ranks={:?}), so the \
             classifier cannot recognise Prop 7.3's canonical example and `not the \
             tangential-touch case` is undecidable on every other cell too",
            w_class.rank,
            w_class.border_rank,
            w_class.regular_pairings,
            w_class.local_ranks
        );

        let branch_hits: u64 = arms.iter().map(|a| a.tally.disc_zero).sum();
        let border_two: u64 = arms.iter().map(|a| a.tally.border_rank_two).sum();
        let rank_three: u64 = arms.iter().map(|a| a.tally.rank_hits[3]).sum();
        let rank_two: u64 = arms.iter().map(|a| a.tally.rank_hits[2]).sum();
        assert!(
            rank_two > 0 && rank_three > 0,
            "VOID: the branch's hits are all one rank (rank2={rank_two}, \
             rank3={rank_three}) across {branch_hits} hits, so `border_rank_two` \
             ({border_two}) and `true_rank_three` are one column wearing two names and \
             C2's `at least one is not the tangential-touch case` is decided by the \
             fixture rather than by the arithmetic"
        );
        assert!(
            border_two > rank_three,
            "VOID: every border-rank-2 hit is also rank 3 (border_rank_two={border_two}, \
             true_rank_three={rank_three}), so the rank-2 stratum of the hypersurface is \
             unpopulated and the two columns cannot be read apart"
        );

        let generic_arm = arms
            .iter()
            .find(|a| a.field == "synthetic_generic")
            .expect("the generic control arm was pushed above");
        let generic_misses = generic_arm.tally.two_roots + generic_arm.tally.disc_negative;
        assert!(
            generic_misses > 0,
            "VOID: every cell of the generic control arm reached the discriminant-zero \
             branch, so `the branch fires` is a tautology on this instrument and C1's \
             count measures the fixture rather than the crate"
        );

        let mirror_u: u64 =
            arms.iter().map(|a| a.tally.mirror_u).sum::<u64>() + hexagon_tally.mirror_u;
        let mirror_coordinate: u64 = arms.iter().map(|a| a.tally.mirror_coordinate).sum::<u64>()
            + hexagon_tally.mirror_coordinate;
        let mirror_mask: u64 =
            arms.iter().map(|a| a.tally.mirror_mask).sum::<u64>() + hexagon_tally.mirror_mask;
        let mirror_interior: u64 = arms.iter().map(|a| a.tally.mirror_interior).sum::<u64>()
            + hexagon_tally.mirror_interior;
        let mirror_hexagon: u64 =
            arms.iter().map(|a| a.tally.mirror_hexagon).sum::<u64>() + hexagon_tally.mirror_hexagon;
        assert!(
            mirror_u == 0
                && mirror_coordinate == 0
                && mirror_mask == 0
                && mirror_interior == 0
                && mirror_hexagon == 0,
            "VOID: the bench's mirror of the private `roots` disagrees with the shipped \
             BodySaddles (u={mirror_u}, coordinates={mirror_coordinate}, \
             mask={mirror_mask}, interior_vertex={mirror_interior}, \
             inner_hexagon={mirror_hexagon}), so every branch count, rank and triangle \
             in this file is about a different function than the one at \
             trilinear.rs:236-267"
        );

        let fan_checks: u64 = arms.iter().map(|a| a.tally.fan_tunnel_checks).sum::<u64>()
            + hexagon_tally.fan_tunnel_checks;
        let fan_disagreements: u64 = arms
            .iter()
            .map(|a| a.tally.fan_tunnel_disagreements)
            .sum::<u64>()
            + hexagon_tally.fan_tunnel_disagreements;
        assert!(
            fan_checks > 0,
            "VOID: no six-saddle cell was found in any arm nor in {hexagon_attempts} \
             hexagon-search draws, so the fan_tunnel mirror C3's triangle counts are \
             built on was never compared against the shipped function"
        );
        assert_eq!(
            fan_disagreements, 0,
            "VOID: the fan_tunnel mirror disagreed with the shipped fan_tunnel on \
             {fan_disagreements} of {fan_checks} six-saddle cells, so \
             `mesh_delta_triangles` is a difference of two numbers only one of which is \
             the crate's"
        );

        let max_patch: u64 = arms
            .iter()
            .map(|a| a.tally.max_patch)
            .max()
            .unwrap_or(0)
            .max(hexagon_tally.max_patch);
        assert!(
            max_patch <= MAX_PATCH_TRIANGLES as u64,
            "VOID: a model asked for {max_patch} triangles in one cell against \
             MAX_PATCH_TRIANGLES = {MAX_PATCH_TRIANGLES}, which the shipped path could \
             not emit -- ✗50's defect in the other direction, and a count the crate \
             cannot produce is not a mesh delta"
        );

        let class_disagreements: u64 = arms.iter().map(|a| a.tally.class_disagreements).sum();
        let local_rank_anomalies: u64 = arms.iter().map(|a| a.tally.local_rank_anomalies).sum();
        let regular_anomalies: u64 = arms.iter().map(|a| a.tally.regular_anomalies).sum();
        assert!(
            class_disagreements == 0 && local_rank_anomalies == 0 && regular_anomalies == 0,
            "VOID: the two classification routes disagree (rank {class_disagreements}, \
             impossible flattening patterns {local_rank_anomalies}, impossible \
             regular-pairing counts {regular_anomalies}), so `border_rank_two`, \
             `true_rank_three` and `w_state_like` are three names for an unresolved \
             disagreement rather than for one classification checked twice"
        );

        let exact_sign_disagreements: u64 =
            arms.iter().map(|a| a.tally.exact_sign_disagreements).sum();
        let entry_bound_violations: u64 = arms.iter().map(|a| a.tally.entry_bound_violations).sum();
        assert!(
            exact_sign_disagreements == 0 && entry_bound_violations == 0,
            "VOID: the synthetic arms are not exact after all \
             (i128-vs-f64 sign disagreements {exact_sign_disagreements}, corner values \
             outside +-{ENTRY_BOUND}: {entry_bound_violations}), so the `exact \
             arithmetic where you can` half of C2 is unmet and the classification of \
             every synthetic hit is a claim about rounded numbers"
        );

        for arm in &arms {
            let t = &arm.tally;
            assert_eq!(
                t.delta_zero,
                t.disc_zero + t.a_zero_no_root,
                "VOID: on `{}` at {} the discriminant-zero cells ({}) are not the \
                 branch hits ({}) plus the no-root `a == 0` cells ({}), so the two \
                 exits do not partition the hypersurface and one of the three branch \
                 counters is counting something else",
                arm.field,
                arm.resolution,
                t.delta_zero,
                t.disc_zero,
                t.a_zero_no_root
            );
            assert_eq!(
                t.roots_true - t.roots_reported,
                t.disc_zero,
                "VOID: on `{}` at {} `roots_true - roots_reported` is {} against \
                 {} branch hits, so the multiplicity model is not the shipped model \
                 plus one root per hit and C3 is comparing two unrelated counts",
                arm.field,
                arm.resolution,
                t.roots_true - t.roots_reported,
                t.disc_zero
            );
            assert_eq!(
                t.cells,
                t.a_zero_no_root + t.a_zero_linear + t.disc_negative + t.disc_zero + t.two_roots,
                "VOID: on `{}` at {} the five exits of `roots` do not account for all \
                 {} cells, so the branch census has a hole in it",
                arm.field,
                arm.resolution,
                t.cells
            );
        }

        // ── the global verdicts ─────────────────────────────────────────────
        //
        // C1 and C2 are global by construction and carry the same value on
        // every row; C3 is per row. C1 reads the reference-field arms only,
        // because the clause says `on at least one reference field at f64` and
        // the synthetic cells are a fixture rather than a field.
        let c1 = arms[..reference_arms].iter().any(|a| a.tally.disc_zero > 0);
        let c2 = branch_hits > 0 && rank_three > 0;
        let c3_any = arms.iter().any(|a| a.tally.delta_triangles != 0);
        let c1_hits: u64 = arms[..reference_arms]
            .iter()
            .map(|a| a.tally.disc_zero)
            .sum();
        let rank_one_hits: u64 = arms.iter().map(|a| a.tally.rank_hits[1]).sum();
        let rank_zero_hits: u64 = arms.iter().map(|a| a.tally.rank_hits[0]).sum();
        let wall_ns = started.elapsed().as_secs_f64() * 1e9;

        for arm in &arms {
            let t = &arm.tally;
            // C3's clause is a disjunction -- `changes at least one triangle ...
            // or provably changes none, and which of those is true is stated as
            // the result` -- so the row's verdict is that its mesh consequence
            // is *established*, which fails exactly when some hit's corrected
            // model has no triangulation rule to count. The falsifier's
            // stricter reading, `C3 by no mesh change`, is the separate
            // `c3_mesh_changed` column; the two are recorded side by side
            // rather than one being chosen quietly.
            let c3 = t.separate_disks == 0;
            run.record(&[
                // ── the registered columns, in registration order ──────────
                ("field", arm.field.clone()),
                ("resolution", arm.resolution.to_string()),
                ("discriminant_zero_hits", t.disc_zero.to_string()),
                ("a_zero_hits", t.a_zero().to_string()),
                ("double_root_hits", t.root_inside.to_string()),
                ("border_rank_two", t.border_rank_two.to_string()),
                ("true_rank_three", t.rank_hits[3].to_string()),
                ("w_state_like", t.w_class.to_string()),
                ("roots_reported", t.roots_reported.to_string()),
                ("roots_true", t.roots_true.to_string()),
                ("mesh_delta_triangles", t.delta_triangles.to_string()),
                ("c1_holds", c1.to_string()),
                ("c2_holds", c2.to_string()),
                ("c3_holds", c3.to_string()),
                // ── extras (M-273) ─────────────────────────────────────────
                //
                // The population, and whether its arithmetic is exact.
                ("cells", t.cells.to_string()),
                ("exact_arithmetic", arm.exact.to_string()),
                // The branch census, split by exit rather than by predicate.
                ("a_zero_no_root_hits", t.a_zero_no_root.to_string()),
                ("a_zero_linear_hits", t.a_zero_linear.to_string()),
                ("disc_negative_cells", t.disc_negative.to_string()),
                ("two_root_cells", t.two_roots.to_string()),
                ("delta_zero_cells", t.delta_zero.to_string()),
                ("delta_negative_cells", t.delta_negative.to_string()),
                // How close the misses came, so a zero is a measurement.
                (
                    "closest_nonzero_disc_rel",
                    t.closest_nonzero
                        .map_or_else(|| "none".to_string(), |v| format!("{v:.6e}")),
                ),
                ("disc_within_one_ulp", t.within_one_ulp.to_string()),
                // C2, both routes and their cross-check.
                ("branch_rank_zero_hits", t.rank_hits[0].to_string()),
                ("branch_rank_one_hits", t.rank_hits[1].to_string()),
                ("branch_rank_two_hits", t.rank_hits[2].to_string()),
                ("branch_rank_three_hits", t.rank_hits[3].to_string()),
                ("class_disagreements", t.class_disagreements.to_string()),
                ("local_rank_anomalies", t.local_rank_anomalies.to_string()),
                ("regular_anomalies", t.regular_anomalies.to_string()),
                (
                    "exact_sign_disagreements",
                    t.exact_sign_disagreements.to_string(),
                ),
                (
                    "entry_bound_violations",
                    t.entry_bound_violations.to_string(),
                ),
                // C3, and the arithmetic behind `provably none`.
                ("hits_all_coords_inside", t.all_coords_inside.to_string()),
                ("hits_trilinear_eligible", t.trilinear_eligible.to_string()),
                ("hexagon_gained", t.hexagon_gained.to_string()),
                ("m2_separate_disks", t.separate_disks.to_string()),
                ("m1_triangles", t.m1_triangles.to_string()),
                ("m2_triangles", t.m2_triangles.to_string()),
                ("c3_mesh_changed", (t.delta_triangles != 0).to_string()),
                ("c3_provably_none", (t.all_coords_inside == 0).to_string()),
                // The mirror's own report, per arm.
                ("mirror_u_disagreements", t.mirror_u.to_string()),
                (
                    "mirror_coordinate_disagreements",
                    t.mirror_coordinate.to_string(),
                ),
                ("mirror_mask_disagreements", t.mirror_mask.to_string()),
                (
                    "mirror_interior_disagreements",
                    t.mirror_interior.to_string(),
                ),
                ("mirror_hexagon_disagreements", t.mirror_hexagon.to_string()),
                ("fan_tunnel_cross_checks", t.fan_tunnel_checks.to_string()),
                (
                    "fan_tunnel_disagreements",
                    t.fan_tunnel_disagreements.to_string(),
                ),
                ("max_patch_triangles_seen", t.max_patch.to_string()),
                // Global quantities, identical on every row and labelled so in
                // the header.
                ("c1_any_reference_field", c1.to_string()),
                ("c1_reference_field_hits_total", c1_hits.to_string()),
                ("c2_global", c2.to_string()),
                ("c3_any_arm_changed", c3_any.to_string()),
                ("branch_hits_total", branch_hits.to_string()),
                ("border_rank_two_total", border_two.to_string()),
                ("true_rank_three_total", rank_three.to_string()),
                ("branch_rank_two_total", rank_two.to_string()),
                ("branch_rank_one_total", rank_one_hits.to_string()),
                ("branch_rank_zero_total", rank_zero_hits.to_string()),
                (
                    "rank_one_reachable_on_branch",
                    (rank_one_hits > 0).to_string(),
                ),
                ("bisep_u_branch_hits", bisep_branch[0].to_string()),
                ("bisep_v_branch_hits", bisep_branch[1].to_string()),
                ("bisep_w_branch_hits", bisep_branch[2].to_string()),
                // The registered vacuity fixture, readable from the file.
                (
                    "w_state_corner",
                    W_STATE
                        .iter()
                        .map(ToString::to_string)
                        .collect::<Vec<_>>()
                        .join("|"),
                ),
                ("w_state_a", w_mirror.a.to_string()),
                ("w_state_b", w_mirror.b.to_string()),
                ("w_state_c", w_mirror.c.to_string()),
                ("w_state_disc_f64", w_mirror.disc.to_string()),
                ("w_state_delta_i128", w_delta.to_string()),
                ("w_state_repo_delta_i128", w_repo_delta.to_string()),
                ("w_state_branch", format!("{:?}", w_mirror.branch)),
                ("w_state_double_root", w_mirror.coordinate[0][0].to_string()),
                ("w_state_rank", w_class.rank.to_string()),
                ("w_state_border_rank", w_class.border_rank.to_string()),
                (
                    "w_state_regular_pairings",
                    w_class.regular_pairings.to_string(),
                ),
                (
                    "w_state_local_ranks",
                    w_class
                        .local_ranks
                        .iter()
                        .map(ToString::to_string)
                        .collect::<Vec<_>>()
                        .join("|"),
                ),
                // The identity, and the fixture's own size.
                ("symbolic_identity_holds", symbolic_identity.to_string()),
                ("cayley_terms", cayley.terms().to_string()),
                ("cayley_total_degree", cayley.total_degree().to_string()),
                ("gl2_generators", generators.len().to_string()),
                ("w_orbit_cells", orbit.len().to_string()),
                ("sampled_cells_per_stratum", SAMPLED_CELLS.to_string()),
                ("hexagon_control_cells", hexagon_cells.len().to_string()),
                ("hexagon_control_attempts", hexagon_attempts.to_string()),
                (
                    "hexagon_control_cross_checks",
                    hexagon_tally.fan_tunnel_checks.to_string(),
                ),
                ("fan_tunnel_cross_checks_total", fan_checks.to_string()),
                ("seed", format!("{SEED:#018x}")),
                (
                    "resolutions",
                    RESOLUTIONS
                        .iter()
                        .map(ToString::to_string)
                        .collect::<Vec<_>>()
                        .join("|"),
                ),
                // Time, recorded beside the verdicts and read by nothing.
                ("wall_ns", format!("{wall_ns:.0}")),
            ]);
        }
    });
}
