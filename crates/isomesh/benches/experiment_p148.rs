//! **P-148 — metric interpolation across a chunk seam, where averaging matrices
//! is not averaging.**
//!
//! Ticket: R-148. Pre-registered before this harness existed.
//!
//! ```bash
//! cargo bench --bench experiment_p148
//! ```
//!
//! Writes `docs/experiments/p-148.csv`.
//!
//! # What was missing
//!
//! `common::metric` (owned by R-146) ships both interpolants and its own doc
//! comment names this row as their referee: *"interpolating between two metrics
//! whose anisotropies point in different directions swells `det`, which in mesh
//! terms means the seam quietly asks for coarser elements than either side did.
//! Measuring that swell is P-148 C1"* (`benches/common/metric.rs:542-548`).
//! Nothing in the repository had measured it, and nothing had asked the second
//! question either: a matrix logarithm is three transcendentals deep, and
//! `M-32` says a chunk seam is bit-exact only under a special cell size.
//!
//! **What `M-32` actually says**, quoted rather than paraphrased
//! (`FINDINGS.md:1145`):
//!
//! > *Chunk seams are bit-exact only when the cell size is a power of two.* Two
//! > adjacent chunks meshed independently agree on **16 of 16** shared-plane
//! > vertices bit-for-bit at `h = 0.125`, and on **0 of 14** at `h = 4/35` —
//! > worst gap `1.57e-16` world units. Cause: an extractor computes
//! > `origin + h·local`, so chunk `c`'s last plane is `(o + h·cn) + h·n` while
//! > `c+1`'s first is `o + h·(c+1)n` — equal by algebra, not by IEEE. **22% of
//! > 200,000 random `(origin, h, cells, chunk)` combinations disagree**, by one
//! > or two ulp.
//!
//! Two later findings sharpen it and both are load-bearing here. `M-49` finds
//! the same defect in `ChunkLayout::cell_of`. And `R-004`'s sweep
//! (`FINDINGS.md:2606-2612`) found that **the seam plane is the wrong unit of
//! analysis**: at `h = 0.1` and `1/12` the seam plane's own coordinate agreed
//! *bit for bit* while only 24 of 92 vertices did, because the disagreement is
//! in *every* coordinate two blocks reach by different groupings, not in the one
//! along the seam normal. This fixture reproduces that too — `seam_plane_ulps`
//! is **0 on all six geometries**, including the three that disagree by 7–30 ulp
//! about the cell *centres* where the metric is actually sampled. A harness that
//! had tested the plane coordinate alone would have reported bit-exactness and
//! been wrong.
//!
//! `M-32`'s own methodological lesson is quoted in the rules table at
//! `FINDINGS.md:1825`: *"the non-power-of-two seam test first used `h = 4/33`,
//! which looks irregular and lands in the 78% of cases that happen to agree
//! exactly. It passed while proving nothing."* So this fixture's spacings were
//! **searched** and every one of the six is asserted below.
//!
//! # The seam fixture
//!
//! Two adjacent chunks over a reference field, meshed and metricised
//! independently, meeting on one shared face.
//!
//! * The **transverse** extent of both chunks is the field's own
//!   `ReferenceField::domain()`, so the seam plane carries the whole cross
//!   section of the surface. Along the **seam axis** each chunk is half the
//!   domain, so the two together span it exactly and the shared face sits at the
//!   domain centre — where every one of the eight reference fields has surface.
//! * Chunk origins come from `ChunkLayout::world_of_sample`, which
//!   `chunk.rs:205-208` calls *"the single place a sample's world position is
//!   defined"*: `origin + cell_size · sample`. Chunk A's base sample is
//!   `CHUNK_BASE = 16` along the axis and chunk B's is `CHUNK_BASE +
//!   chunk_cells`, so the two origins are two different multiplications of one
//!   `h` — `M-32`'s premise, set up through the crate's own accessor rather than
//!   by hand.
//! * The metric is piecewise per **cell**, so the two values being interpolated
//!   across the shared face are the centres of A's last cell and B's first cell,
//!   one `h` apart. Each chunk reaches the other's cell through a **ghost
//!   index** in its own local frame, which is how a chunked mesher actually gets
//!   neighbour data:
//!
//!   | quantity | chunk A's expression | chunk B's expression |
//!   |---|---|---|
//!   | A's own cell centre | `A_o + h·(n − 0.5)` | `B_o + h·(−0.5)` |
//!   | B's own cell centre | `A_o + h·(n + 0.5)` | `B_o + h·0.5` |
//!
//!   with `A_o = o + h·16`, `B_o = o + h·(16 + n)`. Equal by algebra, not by
//!   IEEE — exactly `M-32`, one cell off the plane.
//!
//! Measured before this file was written, over the three domain extents the
//! roster has (`4`, `14`, `16`) and the six geometries below. `ghost` and `own`
//! are the ulp gaps between the two columns of that table; `plane` is the gap in
//! the shared sample plane's own coordinate:
//!
//! | extent | arm | cells/axis | `h` | power of two | `ghost` | `own` | `plane` |
//! |---|---|---|---|---|---|---|---|
//! | 4 | pow2 | 32 | `0.125` | yes | 0 | 0 | 0 |
//! | 4 | non | 34 | `2/17` | no | **30** (2.08e-16) | **2** | 0 |
//! | 14 | pow2 | 56 | `0.25` | yes | 0 | 0 | 0 |
//! | 14 | non | 58 | `7/29` | no | **7** (9.71e-17) | **7** | 0 |
//! | 16 | pow2 | 32 | `0.5` | yes | 0 | 0 | 0 |
//! | 16 | non | 34 | `8/17` | no | **30** (8.33e-16) | **2** | 0 |
//!
//! The power-of-two arm takes `h = 2^⌊log₂(extent/32)⌋`, which is a genuine
//! power of two on all three extents rather than merely a dyadic rational — the
//! literal condition `M-32` states, and the reason the cell counts are 32/56/32
//! rather than a single number. The non-power-of-two arm is that cell count plus
//! two, which keeps the transverse window identical and moves `h` off the
//! dyadics; it is 6% coarser and the header says so rather than hiding it,
//! because neither clause is a resolution comparison.
//!
//! # Arms
//!
//! | arm | what it varies | `is_control` |
//! |---|---|---|
//! | `componentwise` | `(1 − t)A + tB`, entry by entry | **yes** — the naive scheme C1 indicts |
//! | `log_euclidean` | `exp((1 − t)log A + t log B)` | no — the intrinsic scheme |
//! | `cell_size_power_of_two = true` | `h ∈ {0.125, 0.25, 0.5}` | **yes** — `M-32`'s bit-exact regime |
//! | `cell_size_power_of_two = false` | `h ∈ {2/17, 7/29, 8/17}` | no |
//! | *the A ≡ A control* | one metric interpolated with itself | **yes** — `swell_self_max`, proves the swell formula is not structurally positive |
//! | *the inconsistent-metric control* | each chunk displaces by its **own** metric, no interpolation | **yes** — `control_seam_opens`, proves the seam counter can go non-zero |
//! | *the golden control* | one vertex displaced by the metric's own step | **yes** — `golden_control_hash_moved`, proves `mesh_hash` is live |
//!
//! Eight fields × three seam axes × two cell sizes × two schemes = **96 rows**.
//!
//! # How the swell is measured, and why not through `Sym3::det`
//!
//! Log-Euclidean interpolation is determinant-monotone because `log det` is
//! *linear* along its path: `log det exp(X) = tr X`, so
//! `log det M(t) = (1 − t) log det A + t log det B` **exactly**
//! (`common/metric.rs:557-566`). That geometric path is therefore the reference
//! both schemes are measured against, one definition for both arms:
//!
//! ```text
//!     log_swell(t) = ln det M_scheme(t) − [(1 − t) ln det A + t ln det B]
//!     swell(t)     = exp(log_swell(t)) − 1
//! ```
//!
//! Component-wise interpolation can only ever sit above it — by Minkowski's
//! determinant inequality `det((1−t)A + tB)^{1/3} ≥ (1−t) det A^{1/3} + t det
//! B^{1/3}`, and then AM–GM on the right — with equality only when `A ∝ B`. So
//! `swell ≥ 0` is a theorem for the control arm and `swell ≡ 0` is a theorem for
//! the intrinsic one, and C1 is the question of whether the fields reach the 5%
//! bar in practice.
//!
//! **`ln det` is taken from the spectrum, not from `Sym3::det`.** The cofactor
//! expansion carries absolute error `≈ ε·|λ|³max`, and `metric_lp` floors each
//! `|Hessian eigenvalue|` at `H_FLOOR = 1e-9`, so a cell with one real curvature
//! and two flat directions gets a metric with eigenvalues around
//! `(1.9e-7, 1.9e-7, 1.9e4)` whose determinant is `7e-10` against an
//! `ε·|λ|³max` of `1.6e-3` — six orders of noise on top of the answer. Rotate
//! that matrix off the world axes and `det` is meaningless. Jacobi's eigenvalues
//! keep their relative accuracy on a graded positive-definite matrix, so
//! `ln det = Σ ln λ` is the only form of this quantity worth recording, and
//! `log_determinant_swell_max` is reported beside `determinant_swell_max`
//! because the ratio itself can run to many orders.
//!
//! # C2, and the two ways it could fail
//!
//! C2 asks whether log-Euclidean interpolation *preserves* whatever `M-32`
//! grants. Each chunk computes the interpolant at the same physical point from
//! its own side: A evaluates `interp(A_own, A_ghost, t)`, B evaluates
//! `interp(B_own, B_ghost, 1 − t)`. `bit_exact_seam` is true only when all six
//! entries agree bit-for-bit on **every** seam-face cell and **every** rung of
//! the ladder. Two independent things can break it:
//!
//! 1. **The positions**, which is `M-32` and belongs to the cell size. Predicted
//!    `true` on the power-of-two arm and `false` on the other, identically under
//!    both schemes.
//! 2. **The weight**, which is new and belongs to the ladder. The reversal
//!    computes `1 − (1 − t)`, and that returns `t` exactly only when `t` is
//!    dyadic: `1 − (1 − 0.1) = 0.09999999999999998`. So the registered
//!    `bit_exact_seam` uses a **dyadic** ladder `t = k/8`, `k = 0…8`, holding the
//!    weight exact so that the column measures the cell size and nothing else —
//!    and the non-dyadic ladder `t = k/10` is run beside it and reported as
//!    `bit_exact_seam_nondyadic` with `weight_reversal_exact`. Predicted
//!    **false even at `h = 0.125`**, which is a condition `M-32` does not state
//!    and a consumer of a metric field would have to meet.
//!
//! Neither `Sym3::log`, `Sym3::exp` nor `Sym3::eigen` is a source of divergence
//! on its own: all three are deterministic pure functions of their input bits,
//! and `interp_componentwise`'s and `interp_log_euclidean`'s reversals differ
//! only in the order of an IEEE addition, which is commutative. So C2's answer is
//! predicted to be that the transcendentals add nothing to `M-32`'s loss and take
//! nothing from its guarantee — but that is a prediction, and the column is the
//! measurement.
//!
//! # C3, and what makes it more than a tautology
//!
//! A metric field is on no shipped extraction path, so nothing about the mesh
//! changes when one is built — which would make "the seam remains closed" a zero
//! that could not have been non-zero (`M-44`). The metric is therefore given the
//! job it exists for: **it moves the seam vertices.** Each chunk displaces its
//! own copy of every shared seam vertex by `DISPLACE_CELLS · h` along the
//! smallest-eigenvalue direction of its own interpolant at `t = 0.5` — the
//! direction the metric calls cheap. Then the two chunks are concatenated,
//! welded at `weld::epsilon_for(h) = h·1e-4`, and the boundary edges whose *both*
//! endpoints are shared seam vertices are counted against the same count taken
//! without any displacement.
//!
//! `M-129` is why the count is localised by vertex identity rather than by
//! plane: *"the seam counter's own exclusion list was missing an axis, and it
//! accused Marching Cubes of a defect it does not have"*. The surface leaves this
//! two-block slab through its transverse walls on the open fields, and those
//! boundary edges are the slab ending rather than the seam failing.
//!
//! The counter is `isomesh::validate::validate`'s `MeshReport::boundary_edges`;
//! `validate_features` supplies the same number as a list to localise
//! (`validate.rs:606-608` states the lengths are equal), and the two are asserted
//! equal on every arm so the localisation is known to be reading the counter and
//! not a second one.
//!
//! **The prediction, with the margin rather than just the verdict** (`M-44`'s
//! rule): the seam stays closed under both schemes on both arms, because even on
//! the non-power-of-two arm the two chunks' displacements disagree by at most the
//! ~1e-16 world units the metric inherits from `M-32`, against a weld epsilon of
//! `h·1e-4 ≈ 1.2e-5` — eleven orders of margin. `seam_worst_displacement_gap`
//! and `gap_over_weld_epsilon` are that margin. The inconsistent-metric control
//! opens the same seam by construction, because two *different* metrics choose
//! two different eigenvectors and the gap becomes `O(DISPLACE_CELLS·h)`, a
//! hundred times the epsilon.
//!
//! # `hashes_moved`, predicted zero by construction and measured anyway
//!
//! Nothing here touches `crates/isomesh/src/**`, and a metric field is not on any
//! shipped extraction path, so **`hashes_moved` is predicted 0 by
//! construction**. It is measured on the fixture that defines the phrase: the 24
//! committed `marching_cubes` rows of `crates/isomesh/golden_hashes.json`, read
//! with this file's own one-line scanner because `src/golden.rs` is
//! `#[cfg(test)]`. For each row the shipped extractor is run on the fixture's own
//! grid (`golden.rs:163-165`: `cell_size = (hi[0] − lo[0]) / (samples − 1)`,
//! `shape = [samples; 3]`, `origin = lo`), the metric-interpolation pipeline is
//! then run over that grid's own mid-plane seam under this row's scheme, and the
//! extraction is repeated and re-hashed. Two controls stop the zero from being
//! free: `golden_fixture_matches_shipped` (the committed hash must equal the
//! shipped extractor's, or `hashes_moved` is measured against a stale baseline)
//! and `golden_control_hash_moved` (displacing one vertex by the metric's own
//! step **must** move the hash, or `mesh_hash` is not watching).
//!
//! # SHARE, recomputed before the numbers
//!
//! The registration's SHARE line is *"C1 moves the metric-construction stage,
//! whose share `P-146` C3 reports."* Recomputed from the committed
//! `docs/experiments/p-146.csv` rather than from memory: `metric_share` there is
//! `metric_ms / extract_ms`, and over its 40 rows it runs from **0.332862**
//! (`thin_plate` at 65³) to **8.747515** (`csg_difference` at 17³, with
//! `metric_share_hi` reaching **8.930033**), median about **2.0**. P-146 C3
//! asked for under 15% and was **falsified on all 40 rows**: building the metric
//! already costs between a third of the extractor and nine times it.
//!
//! So the stage C1 moves is not a rounding error in a frame, and the cost C1's
//! own falsifier mentions — *"which would mean the intrinsic scheme is
//! unnecessary here and saves the transcendental cost"* — is priced here rather
//! than asserted: `transcendental_cost_ratio` is the median over five repeats of
//! log-Euclidean's per-interpolation cost against component-wise's, and
//! `metric_ms` and `interp_ms` are the per-fixture build and interpolation costs
//! so the seam's share of that already-expensive stage can be read off directly.
//! No clause here has a timing threshold, so no verdict rests on a clock
//! (`M-280`: this host's governor swings the same binary 1.45×); five repeats,
//! median headline, min and max recorded.
//!
//! # Vacuity controls
//!
//! * **The registered one — the two chunks' metrics genuinely differ.**
//!   `metric_distance_max` is the affine-invariant distance
//!   `‖log A − log B‖_F` over the seam face, and `pairs_with_distinct_metrics`
//!   counts the cells where it is non-zero. Asserted globally and **per field**;
//!   a single degenerate `(field, axis, h)` geometry is *recorded* with
//!   `pairs_with_distinct_metrics = 0` rather than aborting the sweep, which is
//!   P-70's precedent for a clause that is unreachable on part of its population.
//! * **The non-power-of-two arm really does disagree**, and the power-of-two arm
//!   really does agree — `seam_ghost_ulps`, asserted `> 0` and `== 0`
//!   respectively on every geometry. This is `M-32`'s own fixture trap
//!   (`FINDINGS.md:1825`) and the six numbers in the table above are the search
//!   that chose the spacings.
//! * **The swell formula can report zero.** `swell_self_max` interpolates a
//!   metric with itself under both schemes and must come back below
//!   `SELF_SWELL_TOLERANCE`, or a positive swell would be an artefact of the
//!   formula rather than of the scheme (`M-44`'s converse).
//! * **Every interpolant is a metric.** `min_metric_eigenvalue` over every
//!   scheme, every cell and every rung, asserted `> 0` — `ln det` of anything
//!   else is not a number.
//! * **The seam carries vertices.** `seam_pairs` per row, asserted non-zero on at
//!   least one axis per field. `thin_plate` is expected to fail this on its `y`
//!   axis alone: the plate is `0.4` cells thick and centred on `y = 0`
//!   (`fields/mod.rs:606-617`), so the `y` seam plane lies *inside* it and the
//!   surface never crosses. Such a row carries `c3_row_decisive = false`.
//! * **The seam counter can go non-zero.** `control_seam_opens` — two chunks
//!   displacing by their own un-interpolated metrics must open the seam.
//! * **`mesh_hash` can move**, and the golden baseline is the committed fixture:
//!   `golden_control_hash_moved`, `golden_fixture_matches_shipped`,
//!   `golden_rows == 24`.
//! * **The two boundary-edge instruments agree.** `MeshReport::boundary_edges`
//!   and `NonManifoldFeatures::boundary_edges.len()` asserted equal on every arm.
//!
//! # Verdict granularity
//!
//! * **`c1_holds` is global** and identical on all 96 rows — the clause compares
//!   two schemes across the whole field roster. The per-row fact is
//!   `c1_row_swell_above_bar`.
//! * **`c2_holds` is per geometry** — both scheme rows of one
//!   `(field, axis, cell size)` carry the same value, because the clause is a
//!   comparison *between* the schemes on that geometry. `c2_global_holds` is the
//!   conjunction.
//! * **`c3_holds` is per row**, and `c3_row_decisive` says whether the row had a
//!   seam to close. `c3_global_holds` is the conjunction over decisive rows.

#![allow(
    clippy::cast_possible_truncation,
    clippy::cast_precision_loss,
    clippy::float_cmp,
    clippy::too_many_lines
)]

mod common;

use std::collections::BTreeMap;
use std::time::Instant;

use common::metric::{
    H_FLOOR, Sym3, hessian, interp_componentwise, interp_log_euclidean, metric_lp,
};
use isomesh::chunk::ChunkLayout;
use isomesh::fields::ReferenceField;
use isomesh::marching_cubes::MarchingCubes;
use isomesh::validate::{ValidateConfig, mesh_hash, validate, validate_features};
use isomesh::weld::{Welder, epsilon_for};
use isomesh::{MeshBuffer, RuntimeShape3, Sdf, for_each_reference_field};

// ════════════════════════════════════════════════════════════════════════════
// constants
// ════════════════════════════════════════════════════════════════════════════

/// The `L^p` norm the metric is built for. `p = 2` is P-146's choice and the
/// value every Group D row is measured at, so the exponent `−1/(2p + 3)` is
/// `−1/7` here as it is there.
const P_NORM: f64 = 2.0;

/// Global sample index of chunk A's origin along the seam axis.
///
/// Non-zero on purpose. `M-32`'s mechanism needs chunk `c` with `c ≥ 1`: at
/// `c = 0` the origin is the layout origin itself and the two expressions
/// coincide trivially. Sixteen is the chunk base the repo's other seam
/// measurements use (`FINDINGS.md:10739`).
const CHUNK_BASE: i64 = 16;

/// Rungs of the **dyadic** interpolation ladder: `t = k / 8` for `k = 0…8`.
///
/// A power of two, so `1 − t` is exact and `1 − (1 − t) == t`. That is what lets
/// `bit_exact_seam` measure the cell size rather than the weight; see the header.
const T_STEPS: u32 = 8;

/// Rungs of the deliberately **non-dyadic** ladder: `t = k / 10`, interior only.
///
/// Ten is not a power of two, so `1 − (1 − t) ≠ t` and the reversal a neighbour
/// chunk performs perturbs the weight itself. Reported as
/// `bit_exact_seam_nondyadic`.
const NON_DYADIC_STEPS: u32 = 10;

/// C1's bar, from the registration: swell "above 5% on at least one field".
///
/// Applied to **both** halves of the clause — component-wise must exceed it and
/// log-Euclidean must not. Using the registration's own number on both sides
/// rather than inventing a tolerance for the intrinsic scheme is deliberate; the
/// measured value is recorded so a reader sees whether the second half passed by
/// a hair or by fourteen orders.
const C1_BAR: f64 = 0.05;

/// A metric interpolated with itself must swell by less than this.
///
/// Component-wise reaches exactly zero (`0.5x + 0.5x` is `x` in IEEE).
/// Log-Euclidean goes through `ln`, a sum and `exp`, so it lands at round-off.
///
/// **`1e-5`, and the first draft's `1e-9` was a claim about a
/// well-conditioned eigenproblem this row does not have.** The premise of
/// that number was "six orders above the `f64` round-off of a
/// well-conditioned eigenproblem" — but `common::metric` floors Hessian
/// eigenvalues at `H_FLOOR = 1e-9` precisely so that a flat direction stays
/// representable, which makes the metric's condition number about `1e9` on
/// exactly the fields this row cares about. A `log`/`exp` round trip on a
/// matrix with condition number `κ` loses about `κ · ε`, i.e. `1e9 · 2.2e-16
/// ≈ 2e-7` — so `1e-9` was **unreachable by construction** and the control
/// was gating on the eigensolver's accuracy rather than on the formula.
/// Measured on the first run: `box_exact/x/log_euclidean` swells by
/// `5.04e-7`, which is that bound to within a factor of three.
///
/// `1e-5` sits two orders above the measured round-off and **four orders
/// below `C1_BAR = 0.05`**, so the control can still fail without ever
/// touching the clause — which is the property that made `1e-9` attractive
/// and is preserved. The measured self-swell is recorded per row
/// (`swell_self_max`) so the margin is auditable rather than asserted.
const SELF_SWELL_TOLERANCE: f64 = 1e-5;

/// Metric-driven displacement of a seam vertex, in cells.
///
/// A hundred times `ValidateConfig::WELD_EPSILON_REL = 1e-4`, so two chunks that
/// choose *different* directions land far outside the weld and the seam opens —
/// that is the inconsistent-metric control. Two chunks that agree land on the
/// same point whatever the magnitude. One percent of a cell is also small enough
/// that the displaced mesh is still the same mesh.
const DISPLACE_CELLS: f64 = 1e-2;

/// A vertex is in the seam plane when its seam-axis coordinate is within this
/// many cells of the plane.
///
/// Marching Cubes interpolates along a cut edge, and an edge lying *in* the
/// plane has both endpoints at the plane's own coordinate, so the interpolation
/// returns that coordinate exactly. `1e-9` cells is therefore seven orders wider
/// than needed and still twelve orders inside any real geometry.
const SEAM_TOL_CELLS: f64 = 1e-9;

/// Transverse quantum for matching chunk A's seam vertices to chunk B's, in
/// cells.
///
/// `M-377`'s choice, for `M-377`'s reason: four orders finer than any real
/// surface movement and four orders coarser than an ulp (`FINDINGS.md:11510`).
const MATCH_QUANTUM_CELLS: f64 = 1e-6;

/// Largest fraction of a fixture's seam vertices that may collide on a
/// transverse key before the key is judged too coarse rather than the mesh
/// judged coincident.
///
/// **A collision is a real coincident vertex pair.** Every seam vertex lies on
/// the seam plane, so two sharing a transverse key to within
/// [`MATCH_QUANTUM_CELLS`] are one geometric point reached from two different
/// grid edges — M-48's degenerate crossing. Measured on the first run: exactly
/// 1 on `gyroid` at `h = 0.25`, and 0 everywhere else. `0.05` is two orders
/// above that and still low enough that a genuinely coarse quantum, which
/// would collide a large share of a dense seam, fails the control.
const KEY_COLLISION_FRACTION: f64 = 0.05;

/// Repeats of the timed interpolation comparison. No clause reads a clock; five
/// repeats and a median are the house floor for reporting one anyway (`M-280`).
const TIMING_REPEATS: usize = 5;

/// Interpolations per timed repeat.
const TIMING_PAIRS: usize = 512;

/// The golden roster row this bench compares against.
const GOLDEN_ALGORITHM: &str = "marching_cubes";

/// `8` fields × `3` resolutions of `RESOLUTIONS = [17, 25, 33]` (`golden.rs:73`).
const GOLDEN_ROWS: usize = 24;

/// Seam-axis names, for the `seam_axis` column.
const AXIS_NAMES: [&str; 3] = ["x", "y", "z"];

// ════════════════════════════════════════════════════════════════════════════
// the two schemes
// ════════════════════════════════════════════════════════════════════════════

/// Which interpolant a row is about.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Scheme {
    /// `(1 − t)A + tB`, entry by entry. C1's control arm.
    Componentwise,
    /// `exp((1 − t)log A + t log B)`. Determinant-monotone.
    LogEuclidean,
}

impl Scheme {
    /// Both, in CSV order.
    const ALL: [Self; 2] = [Self::Componentwise, Self::LogEuclidean];

    /// The `interpolation_scheme` column's value.
    fn name(self) -> &'static str {
        match self {
            Self::Componentwise => "componentwise",
            Self::LogEuclidean => "log_euclidean",
        }
    }

    /// Index into a two-element per-scheme array.
    fn slot(self) -> usize {
        match self {
            Self::Componentwise => 0,
            Self::LogEuclidean => 1,
        }
    }

    /// Interpolate, through `common::metric` unchanged.
    fn interp(self, a: &Sym3, b: &Sym3, t: f64) -> Sym3 {
        match self {
            Self::Componentwise => interp_componentwise(a, b, t),
            Self::LogEuclidean => interp_log_euclidean(a, b, t),
        }
    }
}

// ════════════════════════════════════════════════════════════════════════════
// spectral helpers
// ════════════════════════════════════════════════════════════════════════════

/// `ln det M` from the spectrum, **not** from `Sym3::det`.
///
/// See the header: the cofactor expansion's `ε·|λ|³max` absolute error swamps
/// the determinant of a metric with a floored eigenvalue, and the whole of C1 is
/// a determinant ratio.
fn ln_det(m: &Sym3) -> f64 {
    let (values, _) = m.eigen();
    values[0].ln() + values[1].ln() + values[2].ln()
}

/// Smallest eigenvalue. `eigen` returns ascending, so this is `values[0]`.
fn min_eigenvalue(m: &Sym3) -> f64 {
    m.eigen().0[0]
}

/// The eigenvector of the smallest eigenvalue — the direction the metric calls
/// cheap, and therefore the direction a metric-driven mesher moves a vertex
/// along. `eigen`'s sign convention makes it deterministic.
fn cheap_direction(m: &Sym3) -> [f64; 3] {
    let (_, vectors) = m.eigen();
    [vectors[0][0], vectors[1][0], vectors[2][0]]
}

/// Frobenius norm of a symmetric matrix from its six stored entries.
fn frobenius(m: &Sym3) -> f64 {
    let [xx, xy, xz, yy, yz, zz] = m.0;
    (xx * xx + yy * yy + zz * zz + 2.0 * (xy * xy + xz * xz + yz * yz)).sqrt()
}

/// The affine-invariant (log-Euclidean) distance `‖log A − log B‖_F`.
///
/// This is the registered vacuity control's column. It is zero **only** when the
/// two metrics are equal, and unlike a Euclidean entry difference it does not
/// call two wildly different anisotropies "close" because their large entries
/// happen to match.
fn log_distance(a: &Sym3, b: &Sym3) -> f64 {
    frobenius(&a.log().add(&b.log().scale(-1.0)))
}

/// Whether either of the two metrics' Hessians has an eigenvalue at
/// `metric_lp`'s floor.
///
/// `common::metric`'s header is explicit that a number derived from a floored
/// direction is the constant talking rather than a measurement, so the swell is
/// reported both over the whole population and over the off-floor part of it.
fn at_floor(h: &Sym3) -> bool {
    let (values, _) = h.eigen();
    values.iter().any(|value| value.abs() <= H_FLOOR)
}

/// Bitwise equality of all six entries.
fn same_bits(a: &Sym3, b: &Sym3) -> bool {
    a.0.iter()
        .zip(b.0.iter())
        .all(|(left, right)| left.to_bits() == right.to_bits())
}

/// Ulp distance between two `f64`s of the same sign, `0` when they are equal.
fn ulps_between(a: f64, b: f64) -> u64 {
    if a.to_bits() == b.to_bits() {
        0
    } else {
        let (lo, hi) = (a.to_bits().min(b.to_bits()), a.to_bits().max(b.to_bits()));
        hi - lo
    }
}

/// Whether `h` is a power of two: a finite positive `f64` is one exactly when
/// its mantissa field is zero.
fn is_power_of_two(h: f64) -> bool {
    h > 0.0 && h.is_finite() && h.to_bits() & ((1u64 << 52) - 1) == 0
}

// ════════════════════════════════════════════════════════════════════════════
// the seam fixture's geometry
// ════════════════════════════════════════════════════════════════════════════

/// The transverse axes of a seam axis, in ascending order.
fn transverse(axis: usize) -> [usize; 2] {
    match axis {
        0 => [1, 2],
        1 => [0, 2],
        _ => [0, 1],
    }
}

/// Two adjacent chunks meeting on one shared face: one seam axis, one cell size.
#[derive(Clone, Copy, Debug)]
struct Geometry {
    /// Seam normal, `0`/`1`/`2`.
    axis: usize,
    /// Cells across the field's domain along every axis.
    total: u32,
    /// Cells per chunk along the seam axis; `total / 2`.
    chunk_cells: u32,
    /// Cell size.
    h: f64,
    /// The layout origin both chunks are indexed from.
    origin: [f64; 3],
    /// Chunk A's extraction origin, `world_of_sample([CHUNK_BASE, 0, 0])`.
    a_origin: [f64; 3],
    /// Chunk B's extraction origin, one chunk further along the seam axis.
    b_origin: [f64; 3],
    /// Whether `h` is a power of two — `M-32`'s stated condition.
    power_of_two: bool,
}

impl Geometry {
    /// Build the fixture for one field domain, one seam axis and one cell count.
    ///
    /// The transverse window is the domain exactly and the shared face sits at
    /// the domain centre. Both chunk origins come from
    /// `ChunkLayout::world_of_sample`, so the two are two different
    /// multiplications of one `h` — see the header's table.
    fn new(lo: [f64; 3], extent: f64, axis: usize, total: u32) -> Self {
        let chunk_cells = total / 2;
        let h = extent / f64::from(total);
        let mut origin = lo;
        origin[axis] = lo[axis] - h * CHUNK_BASE as f64;

        let layout = ChunkLayout::<f64>::new(chunk_cells, h, origin)
            .expect("the seam fixture's chunk layout is valid");
        let mut base_a = [0i64; 3];
        base_a[axis] = CHUNK_BASE;
        let mut base_b = [0i64; 3];
        base_b[axis] = CHUNK_BASE + i64::from(chunk_cells);

        Self {
            axis,
            total,
            chunk_cells,
            h,
            origin,
            a_origin: layout.world_of_sample(base_a),
            b_origin: layout.world_of_sample(base_b),
            power_of_two: is_power_of_two(h),
        }
    }

    /// The world point `along` cells from `base` along the seam axis, at the
    /// centre of seam-face cell `(u, v)` transversely.
    ///
    /// `base` is one chunk's own extraction origin and `along` is that chunk's
    /// own local offset, so the expression is the extractor's
    /// `origin + h·local` — which is exactly why the two chunks disagree.
    fn point(&self, base: [f64; 3], along: f64, u: u32, v: u32) -> [f64; 3] {
        let [first, second] = transverse(self.axis);
        let mut p = base;
        p[self.axis] = base[self.axis] + self.h * along;
        p[first] = self.origin[first] + self.h * (f64::from(u) + 0.5);
        p[second] = self.origin[second] + self.h * (f64::from(v) + 0.5);
        p
    }

    /// Chunk A's local offset to its own last cell centre.
    fn a_own_offset(&self) -> f64 {
        f64::from(self.chunk_cells) - 0.5
    }

    /// Chunk A's local offset to its ghost of chunk B's first cell centre.
    fn a_ghost_offset(&self) -> f64 {
        f64::from(self.chunk_cells) + 0.5
    }

    /// The shared sample plane's coordinate as chunk A computes it.
    fn plane_from_a(&self) -> f64 {
        self.a_origin[self.axis] + self.h * f64::from(self.chunk_cells)
    }

    /// The shared sample plane's coordinate as chunk B computes it: its own
    /// origin.
    fn plane_from_b(&self) -> f64 {
        self.b_origin[self.axis]
    }

    /// The extraction shape of one chunk: half the domain along the seam axis,
    /// the whole of it transversely. `Shape3::size` counts **samples**.
    fn shape(&self) -> RuntimeShape3 {
        let mut size = [self.total + 1; 3];
        size[self.axis] = self.chunk_cells + 1;
        RuntimeShape3::new(size).expect("the seam fixture's shape fits the index space")
    }
}

// ════════════════════════════════════════════════════════════════════════════
// the metric field over the shared face
// ════════════════════════════════════════════════════════════════════════════

/// The four metrics each seam-face cell pair has: two cells, each seen from both
/// chunks.
#[derive(Debug)]
struct SeamMetrics {
    /// Chunk A's last cell, as chunk A computes it.
    a_own: Vec<Sym3>,
    /// Chunk B's first cell, as chunk A's ghost index computes it.
    a_ghost: Vec<Sym3>,
    /// Chunk B's first cell, as chunk B computes it.
    b_own: Vec<Sym3>,
    /// Chunk A's last cell, as chunk B's ghost index computes it.
    b_ghost: Vec<Sym3>,
    /// Whether either Hessian of the pair has an eigenvalue at `H_FLOOR`.
    at_floor: Vec<bool>,
    /// `‖log A − log B‖_F` per pair, chunk A's view.
    distance: Vec<f64>,
    /// Wall time of the build, in milliseconds.
    build_ms: f64,
}

impl SeamMetrics {
    /// Number of seam-face cell pairs.
    fn len(&self) -> usize {
        self.a_own.len()
    }
}

/// Build the four metric fields over the shared face.
fn seam_metrics<F: Sdf<Scalar = f64>>(field: &F, g: &Geometry) -> SeamMetrics {
    let started = Instant::now();
    let capacity = (g.total as usize) * (g.total as usize);
    let mut out = SeamMetrics {
        a_own: Vec::with_capacity(capacity),
        a_ghost: Vec::with_capacity(capacity),
        b_own: Vec::with_capacity(capacity),
        b_ghost: Vec::with_capacity(capacity),
        at_floor: Vec::with_capacity(capacity),
        distance: Vec::with_capacity(capacity),
        build_ms: 0.0,
    };

    for v in 0..g.total {
        for u in 0..g.total {
            let h_a_own = hessian(field, g.point(g.a_origin, g.a_own_offset(), u, v), g.h);
            let h_a_ghost = hessian(field, g.point(g.a_origin, g.a_ghost_offset(), u, v), g.h);
            let h_b_own = hessian(field, g.point(g.b_origin, 0.5, u, v), g.h);
            let h_b_ghost = hessian(field, g.point(g.b_origin, -0.5, u, v), g.h);

            let a_own = metric_lp(&h_a_own, P_NORM);
            let a_ghost = metric_lp(&h_a_ghost, P_NORM);

            out.at_floor
                .push(at_floor(&h_a_own) || at_floor(&h_a_ghost));
            out.distance.push(log_distance(&a_own, &a_ghost));
            out.a_own.push(a_own);
            out.a_ghost.push(a_ghost);
            out.b_own.push(metric_lp(&h_b_own, P_NORM));
            out.b_ghost.push(metric_lp(&h_b_ghost, P_NORM));
        }
    }

    out.build_ms = started.elapsed().as_secs_f64() * 1e3;
    out
}

// ════════════════════════════════════════════════════════════════════════════
// C1 and C2: swell and bit-exactness, per scheme
// ════════════════════════════════════════════════════════════════════════════

/// Everything one scheme produced on one geometry's metric field.
#[derive(Clone, Debug)]
struct Interp {
    /// C1's headline: `max (det M_scheme(t) / geometric path) − 1`.
    swell_max: f64,
    /// The same, averaged over every cell and every interior rung.
    swell_mean: f64,
    /// `max log_swell`, the numerically robust form of `swell_max`.
    log_swell_max: f64,
    /// `swell_max` over pairs where neither Hessian is at `H_FLOOR`, or `None`
    /// when every pair is.
    swell_max_off_floor: Option<f64>,
    /// Pairs contributing to `swell_max_off_floor`.
    off_floor_pairs: u64,
    /// The A ≡ A control: interpolating a metric with itself must not swell.
    swell_self_max: f64,
    /// Smallest eigenvalue seen on any interpolant — asserted positive.
    min_eigenvalue: f64,
    /// C2's registered column, on the dyadic ladder.
    bit_exact: bool,
    /// Cells on which every dyadic rung agreed bit-for-bit.
    bit_exact_pairs: u64,
    /// The smallest dyadic `t` at which some cell disagreed.
    first_failure_t: Option<f64>,
    /// The same test on the deliberately non-dyadic ladder.
    bit_exact_nondyadic: bool,
    /// Wall time of the whole interpolation sweep, in milliseconds.
    interp_ms: f64,
}

/// Run one scheme over one geometry's seam metrics.
fn measure_interp(scheme: Scheme, m: &SeamMetrics) -> Interp {
    let started = Instant::now();

    let mut swell_max = f64::NEG_INFINITY;
    let mut swell_sum = 0.0f64;
    let mut swell_count = 0u64;
    let mut log_swell_max = f64::NEG_INFINITY;
    let mut off_floor_max = f64::NEG_INFINITY;
    let mut off_floor_pairs = 0u64;
    let mut swell_self_max = 0.0f64;
    let mut min_eig = f64::INFINITY;

    let mut bit_exact_pairs = 0u64;
    let mut bit_exact = true;
    let mut first_failure = f64::INFINITY;
    let mut bit_exact_nondyadic = true;

    for index in 0..m.len() {
        let (a, a_ghost) = (&m.a_own[index], &m.a_ghost[index]);
        let (b, b_ghost) = (&m.b_own[index], &m.b_ghost[index]);

        // ── C1: the swell against the determinant-monotone path ──
        let ln_a = ln_det(a);
        let ln_b = ln_det(a_ghost);
        for k in 1..T_STEPS {
            let t = f64::from(k) / f64::from(T_STEPS);
            let mixed = scheme.interp(a, a_ghost, t);
            min_eig = min_eig.min(min_eigenvalue(&mixed));
            let log_swell = ln_det(&mixed) - ((1.0 - t) * ln_a + t * ln_b);
            let swell = log_swell.exp() - 1.0;
            swell_max = swell_max.max(swell);
            swell_sum += swell;
            swell_count += 1;
            log_swell_max = log_swell_max.max(log_swell);
            if !m.at_floor[index] {
                off_floor_max = off_floor_max.max(swell);
                off_floor_pairs += 1;
            }
        }

        // ── the A ≡ A control ──
        let same = scheme.interp(a, a, 0.5);
        swell_self_max = swell_self_max.max(((ln_det(&same) - ln_a).exp() - 1.0).abs());

        // ── C2: the two chunks' interpolants, dyadic ladder ──
        let mut pair_exact = true;
        for k in 0..=T_STEPS {
            let t = f64::from(k) / f64::from(T_STEPS);
            let from_a = scheme.interp(a, a_ghost, t);
            let from_b = scheme.interp(b, b_ghost, 1.0 - t);
            if !same_bits(&from_a, &from_b) {
                pair_exact = false;
                first_failure = first_failure.min(t);
            }
        }
        if pair_exact {
            bit_exact_pairs += 1;
        } else {
            bit_exact = false;
        }

        // ── the non-dyadic ladder, where the weight itself is perturbed ──
        for k in 1..NON_DYADIC_STEPS {
            let t = f64::from(k) / f64::from(NON_DYADIC_STEPS);
            let from_a = scheme.interp(a, a_ghost, t);
            let from_b = scheme.interp(b, b_ghost, 1.0 - t);
            if !same_bits(&from_a, &from_b) {
                bit_exact_nondyadic = false;
            }
        }
    }

    Interp {
        swell_max,
        swell_mean: swell_sum / swell_count as f64,
        log_swell_max,
        swell_max_off_floor: if off_floor_pairs == 0 {
            None
        } else {
            Some(off_floor_max)
        },
        off_floor_pairs,
        swell_self_max,
        min_eigenvalue: min_eig,
        bit_exact,
        bit_exact_pairs,
        first_failure_t: if first_failure.is_finite() {
            Some(first_failure)
        } else {
            None
        },
        bit_exact_nondyadic,
        interp_ms: started.elapsed().as_secs_f64() * 1e3,
    }
}

/// Whether `1 − (1 − t) == t` on every rung of the non-dyadic ladder.
fn weight_reversal_exact() -> bool {
    (1..NON_DYADIC_STEPS).all(|k| {
        let t = f64::from(k) / f64::from(NON_DYADIC_STEPS);
        1.0 - (1.0 - t) == t
    })
}

// ════════════════════════════════════════════════════════════════════════════
// C3: the meshed seam
// ════════════════════════════════════════════════════════════════════════════

/// The two chunks' meshes and the shared seam vertices matching them up.
#[derive(Debug)]
struct SeamMesh {
    /// Chunk A's mesh, as the shipped extractor produced it.
    a: MeshBuffer<f64>,
    /// Chunk B's mesh.
    b: MeshBuffer<f64>,
    /// `(A's vertex, B's vertex, seam-face cell index)` per shared seam vertex.
    pairs: Vec<(u32, u32, usize)>,
    /// Chunk A's vertices lying in the shared plane.
    seam_a: u64,
    /// Chunk B's vertices lying in the shared plane.
    seam_b: u64,
    /// Seam vertices of A with no partner in B, or vice versa.
    unmatched: u64,
    /// Transverse keys claimed by more than one seam vertex.
    key_collisions: u64,
}

/// Extract both chunks and match their shared seam vertices.
fn seam_mesh<F: Sdf<Scalar = f64>>(field: &F, g: &Geometry) -> SeamMesh {
    let shape = g.shape();
    let mut a = MeshBuffer::<f64>::new();
    MarchingCubes::<f64>::new()
        .extract(field, &shape, g.a_origin, g.h, &mut a)
        .expect("chunk A extraction");
    let mut b = MeshBuffer::<f64>::new();
    MarchingCubes::<f64>::new()
        .extract(field, &shape, g.b_origin, g.h, &mut b)
        .expect("chunk B extraction");

    let tolerance = g.h * SEAM_TOL_CELLS;
    let quantum = g.h * MATCH_QUANTUM_CELLS;
    let [first, second] = transverse(g.axis);

    let key_of = |p: [f64; 3]| {
        (
            ((p[first] - g.origin[first]) / quantum).round() as i64,
            ((p[second] - g.origin[second]) / quantum).round() as i64,
        )
    };
    let cell_of = |p: [f64; 3]| {
        let clamp = |x: f64, base: f64| {
            let raw = ((x - base) / g.h).floor();
            let capped = raw.max(0.0).min(f64::from(g.total - 1));
            capped as usize
        };
        clamp(p[second], g.origin[second]) * g.total as usize + clamp(p[first], g.origin[first])
    };

    let mut key_collisions = 0u64;
    let mut from_b: BTreeMap<(i64, i64), u32> = BTreeMap::new();
    let plane_b = g.plane_from_b();
    let mut seam_b = 0u64;
    for (index, p) in b.positions.iter().enumerate() {
        if (p[g.axis] - plane_b).abs() > tolerance {
            continue;
        }
        seam_b += 1;
        if from_b.insert(key_of(*p), index as u32).is_some() {
            key_collisions += 1;
        }
    }

    let plane_a = g.plane_from_a();
    let mut seam_a = 0u64;
    let mut pairs = Vec::new();
    let mut unmatched = 0u64;
    for (index, p) in a.positions.iter().enumerate() {
        if (p[g.axis] - plane_a).abs() > tolerance {
            continue;
        }
        seam_a += 1;
        match from_b.get(&key_of(*p)) {
            Some(&partner) => pairs.push((index as u32, partner, cell_of(*p))),
            None => unmatched += 1,
        }
    }
    unmatched += seam_b - pairs.len() as u64;

    SeamMesh {
        a,
        b,
        pairs,
        seam_a,
        seam_b,
        unmatched,
        key_collisions,
    }
}

/// Which metric each chunk displaces its seam vertices by.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Displacement {
    /// No displacement at all — C3's baseline.
    None,
    /// Each chunk uses its own interpolant at the shared face. This is the row.
    Interpolated(Scheme),
    /// Each chunk uses its **own** metric with no interpolation. The control:
    /// two different metrics choose two different directions, so the two copies
    /// of a seam vertex land `O(DISPLACE_CELLS·h)` apart and the weld cannot
    /// close them.
    Inconsistent,
}

/// One welded two-chunk mesh's seam census.
#[derive(Clone, Copy, Debug)]
struct SeamCensus {
    /// Boundary edges with both endpoints a shared seam vertex.
    open_edges: u64,
    /// `MeshReport::boundary_edges` over the whole welded slab.
    boundary_total: u64,
    /// Shared seam vertices the two chunks disagree about after displacement.
    vertices_moved: u64,
    /// Shared seam vertices that were displaced at all.
    vertices_displaced: u64,
    /// Largest `|A's landing − B's landing|`, in world units.
    worst_gap: f64,
    /// Vertices the weld removed.
    welded_away: u64,
}

/// Displace, concatenate, weld, validate, and count the seam's open edges.
fn seam_census(g: &Geometry, m: &SeamMetrics, mesh: &SeamMesh, mode: Displacement) -> SeamCensus {
    let step = DISPLACE_CELLS * g.h;
    let offsets = |cell: usize| -> ([f64; 3], [f64; 3]) {
        match mode {
            Displacement::None => ([0.0; 3], [0.0; 3]),
            Displacement::Interpolated(scheme) => (
                cheap_direction(&scheme.interp(&m.a_own[cell], &m.a_ghost[cell], 0.5)),
                cheap_direction(&scheme.interp(&m.b_own[cell], &m.b_ghost[cell], 0.5)),
            ),
            Displacement::Inconsistent => (
                cheap_direction(&m.a_own[cell]),
                cheap_direction(&m.b_own[cell]),
            ),
        }
    };

    let mut a = mesh.a.clone();
    let mut b = mesh.b.clone();
    let mut seam_flag_a = vec![false; a.positions.len()];
    let mut seam_flag_b = vec![false; b.positions.len()];
    let mut vertices_moved = 0u64;
    let mut vertices_displaced = 0u64;
    let mut worst_gap = 0.0f64;

    for &(va, vb, cell) in &mesh.pairs {
        seam_flag_a[va as usize] = true;
        seam_flag_b[vb as usize] = true;
        let (dir_a, dir_b) = offsets(cell);
        let mut landing_a = [0.0f64; 3];
        let mut landing_b = [0.0f64; 3];
        let mut moved = false;
        let mut displaced = false;
        for axis in 0..3 {
            landing_a[axis] = a.positions[va as usize][axis] + step * dir_a[axis];
            landing_b[axis] = b.positions[vb as usize][axis] + step * dir_b[axis];
            if landing_a[axis].to_bits() != landing_b[axis].to_bits() {
                moved = true;
            }
            if step * dir_a[axis] != 0.0 || step * dir_b[axis] != 0.0 {
                displaced = true;
            }
            worst_gap = worst_gap.max((landing_a[axis] - landing_b[axis]).abs());
        }
        a.positions[va as usize] = landing_a;
        b.positions[vb as usize] = landing_b;
        if moved {
            vertices_moved += 1;
        }
        if displaced {
            vertices_displaced += 1;
        }
    }

    let mut joined = a;
    let offset = joined.positions.len();
    joined
        .append(&b)
        .expect("the two-chunk slab fits the index space");
    let mut seam_flag = seam_flag_a;
    seam_flag.extend_from_slice(&seam_flag_b);
    debug_assert_eq!(seam_flag.len(), joined.positions.len());
    debug_assert_eq!(offset + b.positions.len(), joined.positions.len());

    let before = joined.positions.len();
    let mut welder = Welder::<f64>::new();
    welder
        .weld(&mut joined, epsilon_for(g.h))
        .expect("the welder accepts an extractor's own buffer");
    let mut welded_seam = vec![false; joined.positions.len()];
    for (input, &flag) in welder.remap().iter().enumerate() {
        if seam_flag[input] {
            welded_seam[flag as usize] = true;
        }
    }

    let config = ValidateConfig::from_cell_size(g.h).expect("a positive finite cell size");
    let report = validate(&joined, &config);
    let (mirror, features) = validate_features(&joined.positions, &joined.indices, &config);
    assert_eq!(
        report.boundary_edges,
        features.boundary_edges.len() as u64,
        "the localised seam counter and MeshReport::boundary_edges must be one counter"
    );
    assert_eq!(
        report.boundary_edges, mirror.boundary_edges,
        "validate and validate_features must agree on the boundary-edge count"
    );

    let open_edges = features
        .boundary_edges
        .iter()
        .filter(|edge| welded_seam[edge[0] as usize] && welded_seam[edge[1] as usize])
        .count() as u64;

    SeamCensus {
        open_edges,
        boundary_total: report.boundary_edges,
        vertices_moved,
        vertices_displaced,
        worst_gap,
        welded_away: (before - joined.positions.len()) as u64,
    }
}

// ════════════════════════════════════════════════════════════════════════════
// the golden fixture
// ════════════════════════════════════════════════════════════════════════════

/// One committed row of `golden_hashes.json`.
#[derive(Clone, Debug)]
struct Golden {
    /// Reference-field name.
    field: String,
    /// Samples per axis.
    samples: u32,
    /// The committed `mesh_hash`.
    hash: u64,
}

/// `field_of`'s shape (`golden.rs:245`), which is not public: find `"key":`,
/// strip an optional quote, cut at the first of `"`, `,` or `}`.
fn json_field(chunk: &str, key: &str) -> String {
    let needle = format!("\"{key}\":");
    let rest = chunk
        .split_once(&needle)
        .map(|(_, tail)| tail)
        .unwrap_or_else(|| panic!("golden_hashes.json object has no `{key}`: {chunk}"));
    if let Some(stripped) = rest.strip_prefix('"') {
        let end = stripped.find('"').expect("a closing quote");
        stripped[..end].to_string()
    } else {
        let end = rest
            .find(['"', ',', '}'])
            .unwrap_or_else(|| panic!("golden_hashes.json value for `{key}` never ends"));
        rest[..end].trim().to_string()
    }
}

/// The committed `marching_cubes` rows of `M-31`'s 216.
///
/// Read from `crates/isomesh/golden_hashes.json`, the file
/// `golden_hashes_are_unchanged` (`golden/tests.rs:59`) gates against, so
/// `hashes_moved` is movement in **the** fixture rather than in a re-derivation
/// of it.
fn golden_marching_cubes() -> Vec<Golden> {
    let path = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("golden_hashes.json");
    let text = std::fs::read_to_string(&path)
        .unwrap_or_else(|e| panic!("golden_hashes.json is the baseline of `hashes_moved`: {e}"));
    let mut out = Vec::new();
    for chunk in text.split('{').skip(1) {
        if json_field(chunk, "algorithm") != GOLDEN_ALGORITHM {
            continue;
        }
        out.push(Golden {
            field: json_field(chunk, "field"),
            samples: json_field(chunk, "samples")
                .parse()
                .expect("a golden resolution"),
            hash: u64::from_str_radix(&json_field(chunk, "hash"), 16).expect("a hex mesh hash"),
        });
    }
    out
}

/// What one field's three golden rows said.
#[derive(Clone, Copy, Debug, Default)]
struct GoldenField {
    /// Committed rows found for this field.
    rows: u64,
    /// Committed hashes the shipped extractor failed to reproduce.
    stale: u64,
    /// Rows whose hash moved with the metric pipeline running, per scheme.
    moved: [u64; 2],
    /// Rows where displacing one vertex by the metric's step failed to move the
    /// hash — the control that makes a zero meaningful.
    control_silent: u64,
    /// Interpolations the pipeline actually performed.
    interpolations: u64,
}

/// Run the metric-interpolation pipeline over a golden grid's own mid-plane, so
/// that `hashes_moved` is measured with the pipeline running rather than
/// asserted from the absence of a call site.
///
/// Returns the accumulated `ln det` of every interpolant and the number of
/// interpolations, both consumed by the caller so the sweep cannot be elided.
fn golden_pipeline<F: Sdf<Scalar = f64>>(
    field: &F,
    lo: [f64; 3],
    cell: f64,
    cells: u32,
    scheme: Scheme,
) -> (f64, u64) {
    let mid = f64::from(cells / 2);
    let mut accumulator = 0.0f64;
    let mut count = 0u64;
    for v in 0..cells {
        for u in 0..cells {
            let transverse_y = lo[1] + cell * (f64::from(u) + 0.5);
            let transverse_z = lo[2] + cell * (f64::from(v) + 0.5);
            let near = [lo[0] + cell * (mid - 0.5), transverse_y, transverse_z];
            let far = [lo[0] + cell * (mid + 0.5), transverse_y, transverse_z];
            let a = metric_lp(&hessian(field, near, cell), P_NORM);
            let b = metric_lp(&hessian(field, far, cell), P_NORM);
            for k in 1..T_STEPS {
                let t = f64::from(k) / f64::from(T_STEPS);
                accumulator += ln_det(&scheme.interp(&a, &b, t));
                count += 1;
            }
        }
    }
    (accumulator, count)
}

/// Extract one golden row through the shipped path, on the fixture's own grid.
fn golden_extract<F>(field: &F, samples: u32) -> (MeshBuffer<f64>, f64, [f64; 3])
where
    F: ReferenceField + Sdf<Scalar = f64>,
{
    let (lo, hi) = field.domain();
    let cell = (hi[0] - lo[0]) / f64::from(samples - 1);
    let shape = RuntimeShape3::new([samples; 3]).expect("a golden shape");
    let mut mesh = MeshBuffer::<f64>::new();
    MarchingCubes::<f64>::new()
        .extract(field, &shape, lo, cell, &mut mesh)
        .expect("golden extraction");
    (mesh, cell, lo)
}

/// Measure `hashes_moved` for one field over its three committed rows.
fn golden_field<F>(field: &F, name: &str, committed: &[Golden]) -> GoldenField
where
    F: ReferenceField + Sdf<Scalar = f64>,
{
    let mut out = GoldenField::default();
    for row in committed.iter().filter(|row| row.field == name) {
        out.rows += 1;
        let (mesh, cell, lo) = golden_extract(field, row.samples);
        if mesh_hash(&mesh) != row.hash {
            out.stale += 1;
        }

        for scheme in Scheme::ALL {
            let (accumulator, count) = golden_pipeline(field, lo, cell, row.samples - 1, scheme);
            out.interpolations += count;
            assert!(
                accumulator.is_finite(),
                "the golden metric pipeline produced a non-finite ln det on {name} at {}",
                row.samples
            );
            let (again, _, _) = golden_extract(field, row.samples);
            if mesh_hash(&again) != row.hash {
                out.moved[scheme.slot()] += 1;
            }
        }

        // The control: a metric-driven displacement of the size this pipeline
        // would apply **must** move the hash, or the zero above is free.
        let mut nudged = mesh;
        if let Some(position) = nudged.positions.first_mut() {
            let metric = metric_lp(&hessian(field, *position, cell), P_NORM);
            let direction = cheap_direction(&metric);
            for axis in 0..3 {
                position[axis] += DISPLACE_CELLS * cell * direction[axis];
            }
        }
        if mesh_hash(&nudged) == row.hash {
            out.control_silent += 1;
        }
    }
    out
}

// ════════════════════════════════════════════════════════════════════════════
// one row of the CSV, before the verdicts are known
// ════════════════════════════════════════════════════════════════════════════

/// One `(field, seam axis, cell size)` geometry, fully measured.
#[derive(Clone, Debug)]
struct Fixture {
    /// Reference-field name.
    field: &'static str,
    /// The geometry.
    g: Geometry,
    /// Seam-face cell pairs.
    seam_cells: u64,
    /// Pairs whose two metrics differ at all — the registered vacuity control.
    distinct_pairs: u64,
    /// `max ‖log A − log B‖_F` over the shared face.
    distance_max: f64,
    /// The same, averaged.
    distance_mean: f64,
    /// The same, minimised.
    distance_min: f64,
    /// Pairs with a floored Hessian eigenvalue.
    at_floor_pairs: u64,
    /// Ulp gap between the two chunks' expressions for B's first cell centre.
    ghost_ulps: u64,
    /// The same for A's last cell centre.
    own_ulps: u64,
    /// The same for the shared sample plane's own coordinate.
    plane_ulps: u64,
    /// `|A's expression − B's expression|` for B's first cell centre, in world
    /// units.
    ghost_delta: f64,
    /// Metric build time over the shared face, in milliseconds.
    metric_ms: f64,
    /// Per-scheme interpolation results.
    interp: [Interp; 2],
    /// Chunk A's vertices in the shared plane.
    seam_a: u64,
    /// Chunk B's vertices in the shared plane.
    seam_b: u64,
    /// Matched shared seam vertices.
    seam_pairs: u64,
    /// Seam vertices with no partner.
    unmatched: u64,
    /// Transverse keys claimed twice.
    key_collisions: u64,
    /// C3's baseline: the seam's open edges with no displacement.
    baseline: SeamCensus,
    /// The inconsistent-metric control.
    control: SeamCensus,
    /// Per-scheme seam censuses.
    census: [SeamCensus; 2],
}

/// Measure one geometry end to end.
fn measure<F>(name: &'static str, field: &F, axis: usize, total: u32) -> Fixture
where
    F: ReferenceField + Sdf<Scalar = f64>,
{
    let (lo, hi) = field.domain();
    let g = Geometry::new(lo, hi[0] - lo[0], axis, total);
    let metrics = seam_metrics(field, &g);
    let mesh = seam_mesh(field, &g);

    let interp = [
        measure_interp(Scheme::Componentwise, &metrics),
        measure_interp(Scheme::LogEuclidean, &metrics),
    ];
    let baseline = seam_census(&g, &metrics, &mesh, Displacement::None);
    let control = seam_census(&g, &metrics, &mesh, Displacement::Inconsistent);
    let census = [
        seam_census(
            &g,
            &metrics,
            &mesh,
            Displacement::Interpolated(Scheme::Componentwise),
        ),
        seam_census(
            &g,
            &metrics,
            &mesh,
            Displacement::Interpolated(Scheme::LogEuclidean),
        ),
    ];

    let distinct_pairs = metrics.distance.iter().filter(|d| **d > 0.0).count() as u64;
    let distance_max = metrics.distance.iter().copied().fold(0.0f64, f64::max);
    let distance_min = metrics
        .distance
        .iter()
        .copied()
        .fold(f64::INFINITY, f64::min);
    let distance_mean = metrics.distance.iter().sum::<f64>() / metrics.len() as f64;

    let ghost_a = g.a_origin[axis] + g.h * g.a_ghost_offset();
    let ghost_b = g.b_origin[axis] + g.h * 0.5;
    let own_a = g.a_origin[axis] + g.h * g.a_own_offset();
    let own_b = g.b_origin[axis] + g.h * -0.5;

    Fixture {
        field: name,
        g,
        seam_cells: metrics.len() as u64,
        distinct_pairs,
        distance_max,
        distance_mean,
        distance_min,
        at_floor_pairs: metrics.at_floor.iter().filter(|f| **f).count() as u64,
        ghost_ulps: ulps_between(ghost_a, ghost_b),
        own_ulps: ulps_between(own_a, own_b),
        plane_ulps: ulps_between(g.plane_from_a(), g.plane_from_b()),
        ghost_delta: (ghost_a - ghost_b).abs(),
        metric_ms: metrics.build_ms,
        interp,
        seam_a: mesh.seam_a,
        seam_b: mesh.seam_b,
        seam_pairs: mesh.pairs.len() as u64,
        unmatched: mesh.unmatched,
        key_collisions: mesh.key_collisions,
        baseline,
        control,
        census,
    }
}

/// The power-of-two cell count for one domain extent: the largest `h = 2^k` with
/// at most `1/32` of the extent, which is a genuine power of two on all three
/// extents the roster has. See the header's table.
fn power_of_two_cells(extent: f64) -> u32 {
    let mut h = 1.0f64;
    while h > extent / 32.0 {
        h *= 0.5;
    }
    let total = (extent / h).round() as u32;
    assert!(
        total.is_multiple_of(2) && is_power_of_two(extent / f64::from(total)),
        "the power-of-two arm needs an even cell count and a power-of-two cell size; \
         extent {extent} gave {total} cells of {}",
        extent / f64::from(total)
    );
    total
}

// ════════════════════════════════════════════════════════════════════════════
// main
// ════════════════════════════════════════════════════════════════════════════

fn main() {
    if !std::env::args().any(|arg| arg == "--bench") {
        return;
    }
    let prereg = isomesh::experiment!("P-148");

    common::experiment::run(prereg, |run| {
        let committed = golden_marching_cubes();
        let mut fixtures: Vec<Fixture> = Vec::new();
        let mut golden: Vec<(&'static str, GoldenField)> = Vec::new();
        let mut timing: Vec<(Sym3, Sym3)> = Vec::new();

        for_each_reference_field!(f64, |name, field| {
            golden.push((name, golden_field(&field, name, &committed)));
            let (lo, hi) = field.domain();
            let pow2 = power_of_two_cells(hi[0] - lo[0]);
            for axis in 0..3 {
                for total in [pow2, pow2 + 2] {
                    let fixture = measure(name, &field, axis, total);
                    if timing.is_empty() {
                        let metrics = seam_metrics(&field, &fixture.g);
                        for index in 0..metrics.len().min(TIMING_PAIRS) {
                            timing.push((metrics.a_own[index], metrics.a_ghost[index]));
                        }
                    }
                    fixtures.push(fixture);
                }
            }
        });

        // ── the cost of the transcendentals, five repeats, median headline ──
        let mut per_scheme_ns = [0.0f64; 2];
        let mut per_scheme_lo = [0.0f64; 2];
        let mut per_scheme_hi = [0.0f64; 2];
        for scheme in Scheme::ALL {
            // Warm up once, then time.
            let mut sink = 0.0f64;
            for (a, b) in &timing {
                sink += scheme.interp(a, b, 0.5).trace();
            }
            let mut samples = Vec::with_capacity(TIMING_REPEATS);
            for _ in 0..TIMING_REPEATS {
                let started = Instant::now();
                for (a, b) in &timing {
                    sink += scheme.interp(a, b, 0.5).trace();
                }
                samples.push(started.elapsed().as_secs_f64() * 1e9 / timing.len() as f64);
            }
            assert!(
                sink.is_finite(),
                "the timed interpolation loop produced a non-finite trace and was therefore \
                 measuring nothing"
            );
            samples.sort_by(f64::total_cmp);
            per_scheme_ns[scheme.slot()] = samples[TIMING_REPEATS / 2];
            per_scheme_lo[scheme.slot()] = samples[0];
            per_scheme_hi[scheme.slot()] = samples[TIMING_REPEATS - 1];
        }
        let transcendental_cost_ratio = per_scheme_ns[1] / per_scheme_ns[0];

        // ════════════════════════════════════════════════════════════════════
        // vacuity controls — every one before the first `run.record`
        // ════════════════════════════════════════════════════════════════════

        assert_eq!(
            committed.len(),
            GOLDEN_ROWS,
            "VOID: golden_hashes.json yielded {} `{GOLDEN_ALGORITHM}` rows rather than \
             {GOLDEN_ROWS}, so the scanner matched the wrong thing and `hashes_moved` would be \
             measured against a baseline that is not the committed fixture",
            committed.len()
        );
        for (name, g) in &golden {
            assert_eq!(
                g.rows, 3,
                "VOID: {name} has {} committed `{GOLDEN_ALGORITHM}` rows rather than 3, so its \
                 `hashes_moved` is measured against nothing",
                g.rows
            );
            assert_eq!(
                g.stale, 0,
                "VOID: the shipped extractor fails to reproduce {} of {name}'s committed golden \
                 hashes, so `hashes_moved` would be measured against a stale fixture and could \
                 not be read as a cost",
                g.stale
            );
            assert_eq!(
                g.control_silent, 0,
                "VOID: on {} of {name}'s golden rows, displacing a vertex by the metric's own \
                 step did NOT move the mesh hash -- so `hashes_moved = 0` is a zero that could \
                 not have been non-zero (M-44), and is equally consistent with a mesh_hash blind \
                 to positions",
                g.control_silent
            );
        }

        // `M-32`'s own fixture trap (FINDINGS.md:1825): a non-power-of-two
        // spacing that *looks* irregular may land in the 78% that agree
        // exactly. All six geometries were searched before this file existed
        // and every one is asserted here.
        for f in &fixtures {
            if f.g.power_of_two {
                assert_eq!(
                    f.ghost_ulps, 0,
                    "VOID: {}/{}/h={} is the power-of-two arm and the two chunks' expressions \
                     for the neighbour's cell centre already disagree by {} ulp, so C2 has no \
                     bit-exact regime to preserve",
                    f.field, AXIS_NAMES[f.g.axis], f.g.h, f.ghost_ulps
                );
            } else {
                assert!(
                    f.ghost_ulps > 0,
                    "VOID: {}/{}/h={} is the non-power-of-two arm and the two chunks' \
                     expressions for the neighbour's cell centre agree bit for bit, so this arm \
                     lands in M-32's 78% and proves nothing about a non-power-of-two seam",
                    f.field,
                    AXIS_NAMES[f.g.axis],
                    f.g.h
                );
            }
        }

        // The registered vacuity control. Global, and per field, so that one
        // degenerate geometry is *recorded* rather than aborting the sweep
        // (P-70's precedent); every row carries `pairs_with_distinct_metrics`.
        assert!(
            fixtures.iter().any(|f| f.distinct_pairs > 0),
            "VOID: not one geometry in the sweep has a seam-face cell where the two chunks' \
             metrics differ, so both schemes interpolate between equal endpoints everywhere and \
             C1's swell is measured against nothing"
        );
        for (name, _) in &golden {
            assert!(
                fixtures
                    .iter()
                    .any(|f| f.field == *name && f.distinct_pairs > 0),
                "VOID: every seam-face cell of every {name} geometry has ||log A - log B||_F = 0, \
                 so both schemes interpolate between equal endpoints on this field and its swell \
                 columns are unmeasured"
            );
            assert!(
                fixtures
                    .iter()
                    .any(|f| f.field == *name && f.seam_pairs > 0),
                "VOID: no {name} geometry puts a shared seam vertex on its seam plane, so \
                 `seam_open_edges = 0` on every one of its rows is a zero that could not have \
                 been non-zero (M-44) and C3 says nothing about this field"
            );
        }

        // The swell formula must be able to report zero, and every interpolant
        // must be a metric.
        for f in &fixtures {
            for scheme in Scheme::ALL {
                let arm = &f.interp[scheme.slot()];
                assert!(
                    arm.swell_self_max < SELF_SWELL_TOLERANCE,
                    "VOID: {}/{}/{} interpolating a metric with ITSELF swells its determinant by \
                     {:e}, above {SELF_SWELL_TOLERANCE:e} -- so a positive swell is an artefact \
                     of the formula rather than of the scheme",
                    f.field,
                    AXIS_NAMES[f.g.axis],
                    scheme.name(),
                    arm.swell_self_max
                );
                assert!(
                    arm.min_eigenvalue > 0.0,
                    "VOID: {}/{}/{} produced an interpolant whose smallest eigenvalue is {:e}, \
                     so it is not a metric and `ln det` of it is not a number",
                    f.field,
                    AXIS_NAMES[f.g.axis],
                    scheme.name(),
                    arm.min_eigenvalue
                );
            }
        }

        // The seam counter must be able to go non-zero — demonstrated on the
        // fields where it CAN, with the others named.
        //
        // **The control needs the two chunks' own metrics to actually differ.**
        // The inconsistent arm displaces each chunk's shared vertex by its own
        // un-interpolated metric, so the weld only fails to merge them when
        // those two metrics point somewhere different. On `box_exact` they do
        // not: an axis-aligned seam through a polyhedron has the identical
        // Hessian on both sides, so `vertices_moved` is 0, the weld succeeds,
        // and `open_edges` cannot rise. That is a property of the field, not a
        // broken control, and the first run measured it: 60 shared seam
        // vertices, 0 moved, 0 open edges against a baseline of 0.
        //
        // The control is therefore asserted where it is reachable, the
        // unreachable fields are named, and at least one field must
        // demonstrate the counter rising — which is the assertion that
        // licenses reading a 0 elsewhere as geometry rather than a dead
        // instrument.
        let mut counter_demonstrated_on: Vec<&'static str> = Vec::new();
        let mut counter_unreachable_on: Vec<&'static str> = Vec::new();
        for f in &fixtures {
            if f.seam_pairs == 0 {
                continue;
            }
            // **A key collision is a coincident vertex pair, not a broken
            // key, and the distinction is measured.** The quantum is `1e-6`
            // cells (M-377's), and every seam vertex lies on the seam plane,
            // so two vertices sharing a transverse key are the SAME geometric
            // point reached from two different grid edges — M-48's degenerate
            // crossing, the identical mechanism P-145's Euler cross-check ran
            // into. Measured: 1 collision on gyroid at h=0.25 out of its seam
            // vertices. The correspondence is therefore a bijection on
            // POINTS, which is what C3 counts over, and the collision count
            // is a column.
            //
            // A large collision fraction WOULD mean the key is too coarse, so
            // that remains gated.
            let collision_fraction = if f.seam_pairs > 0 {
                f.key_collisions as f64 / f.seam_pairs as f64
            } else {
                0.0
            };
            assert!(
                collision_fraction <= KEY_COLLISION_FRACTION,
                "VOID: {}/{}/h={} matched {} of {} seam vertices to an already-claimed \
                 transverse key ({:.3} of them, above {KEY_COLLISION_FRACTION}), which is \
                 too many to be M-48 coincidences and means the {MATCH_QUANTUM_CELLS}-cell \
                 quantum is too coarse to separate distinct seam vertices",
                f.field,
                AXIS_NAMES[f.g.axis],
                f.g.h,
                f.key_collisions,
                f.seam_pairs,
                collision_fraction
            );
            if f.control.vertices_moved == 0 {
                counter_unreachable_on.push(f.field);
                continue;
            }
            counter_demonstrated_on.push(f.field);
            assert!(
                f.control.open_edges > f.baseline.open_edges,
                "VOID: {}/{}/h={} has {} shared seam vertices, and displacing them by each \
                 chunk's OWN un-interpolated metric left the seam's open edges at {} against a \
                 baseline of {} -- so `seam_open_edges` could not have risen and C3 is unmeasured",
                f.field,
                AXIS_NAMES[f.g.axis],
                f.g.h,
                f.seam_pairs,
                f.control.open_edges,
                f.baseline.open_edges
            );
            assert_eq!(
                f.key_collisions, 0,
                "VOID: {}/{}/h={} matched two seam vertices to one transverse key, so the \
                 A-to-B correspondence C3 is counted over is not a bijection",
                f.field, AXIS_NAMES[f.g.axis], f.g.h
            );
        }
        assert!(
            !counter_demonstrated_on.is_empty(),
            "VOID: no fixture demonstrated `seam_open_edges` rising under an inconsistent \
             displacement ({} field(s) unreachable because both chunks' own metrics agree on \
             the seam: {}), so every C3 zero in this run is a zero that could not have been \
             non-zero (M-44)",
            counter_unreachable_on.len(),
            counter_unreachable_on.join("|")
        );

        assert!(
            !weight_reversal_exact(),
            "VOID: 1 - (1 - t) == t on every rung of the non-dyadic ladder, so \
             `bit_exact_seam_nondyadic` is testing the same thing as `bit_exact_seam` and the \
             weight mechanism the header claims is unreachable"
        );

        // ════════════════════════════════════════════════════════════════════
        // verdicts
        // ════════════════════════════════════════════════════════════════════

        let scheme_swell_max = |slot: usize| {
            fixtures
                .iter()
                .map(|f| f.interp[slot].swell_max)
                .fold(f64::NEG_INFINITY, f64::max)
        };
        let componentwise_swell_max = scheme_swell_max(0);
        let log_euclidean_swell_max = scheme_swell_max(1);
        let c1_holds = componentwise_swell_max > C1_BAR && log_euclidean_swell_max <= C1_BAR;

        let c2_of = |f: &Fixture| {
            let cw = f.interp[0].bit_exact;
            let le = f.interp[1].bit_exact;
            le == cw && (!f.g.power_of_two || le)
        };
        let c2_global_holds = fixtures.iter().all(c2_of);

        let hashes_moved = |name: &str, slot: usize| {
            golden
                .iter()
                .find(|(field, _)| *field == name)
                .map_or(0, |(_, g)| g.moved[slot])
        };
        let golden_interpolations: u64 = golden.iter().map(|(_, g)| g.interpolations).sum();

        let mut c3_global_holds = true;
        let mut decisive_rows = 0u64;
        for f in &fixtures {
            for scheme in Scheme::ALL {
                if f.seam_pairs == 0 {
                    continue;
                }
                decisive_rows += 1;
                if f.census[scheme.slot()].open_edges > f.baseline.open_edges {
                    c3_global_holds = false;
                }
            }
        }

        println!("\n-- P-148 aggregates --");
        println!(
            "  C1  componentwise swell max {componentwise_swell_max:e}, log-Euclidean \
             {log_euclidean_swell_max:e}, bar {C1_BAR} -> {c1_holds}"
        );
        println!(
            "  C2  {} of {} geometries preserve M-32's bit-exactness -> {c2_global_holds}",
            fixtures.iter().filter(|f| c2_of(f)).count(),
            fixtures.len()
        );
        println!(
            "  C3  {decisive_rows} decisive rows, seam closed on all of them: {c3_global_holds}"
        );
        println!(
            "      interpolation cost/pair: componentwise {:.1} ns, log-Euclidean {:.1} ns, \
             ratio {transcendental_cost_ratio:.2}x",
            per_scheme_ns[0], per_scheme_ns[1]
        );
        println!("      golden interpolations run: {golden_interpolations}");

        // ════════════════════════════════════════════════════════════════════
        // rows
        // ════════════════════════════════════════════════════════════════════

        for f in &fixtures {
            let c2_holds = c2_of(f);
            for scheme in Scheme::ALL {
                let slot = scheme.slot();
                let arm = &f.interp[slot];
                let census = &f.census[slot];
                let decisive = f.seam_pairs > 0;
                let c3_holds = census.open_edges <= f.baseline.open_edges;

                run.record(&[
                    // ── the registration's twelve, in registration order ──
                    ("interpolation_scheme", scheme.name().to_string()),
                    ("seam_axis", AXIS_NAMES[f.g.axis].to_string()),
                    ("determinant_swell_max", format!("{:.6e}", arm.swell_max)),
                    ("determinant_swell_mean", format!("{:.6e}", arm.swell_mean)),
                    ("seam_vertices_moved", census.vertices_moved.to_string()),
                    ("seam_open_edges", census.open_edges.to_string()),
                    ("hashes_moved", hashes_moved(f.field, slot).to_string()),
                    ("bit_exact_seam", arm.bit_exact.to_string()),
                    ("cell_size_power_of_two", f.g.power_of_two.to_string()),
                    ("c1_holds", c1_holds.to_string()),
                    ("c2_holds", c2_holds.to_string()),
                    ("c3_holds", c3_holds.to_string()),
                    // ── extras (M-273) ──
                    ("field", f.field.to_string()),
                    ("cell_size", format!("{:.12e}", f.g.h)),
                    ("cells_per_axis", f.g.total.to_string()),
                    ("cells_per_chunk", f.g.chunk_cells.to_string()),
                    ("p_norm", format!("{P_NORM:.1}")),
                    ("seam_cells", f.seam_cells.to_string()),
                    // the registered vacuity control's columns
                    ("metric_distance_max", format!("{:.6e}", f.distance_max)),
                    ("metric_distance_mean", format!("{:.6e}", f.distance_mean)),
                    ("metric_distance_min", format!("{:.6e}", f.distance_min)),
                    ("pairs_with_distinct_metrics", f.distinct_pairs.to_string()),
                    // C1's arithmetic
                    (
                        "log_determinant_swell_max",
                        format!("{:.6e}", arm.log_swell_max),
                    ),
                    (
                        "determinant_swell_max_off_floor",
                        arm.swell_max_off_floor
                            .map_or_else(|| "none".to_string(), |v| format!("{v:.6e}")),
                    ),
                    ("off_floor_pairs", arm.off_floor_pairs.to_string()),
                    ("at_floor_pairs", f.at_floor_pairs.to_string()),
                    (
                        "at_floor_fraction",
                        format!("{:.6}", f.at_floor_pairs as f64 / f.seam_cells as f64),
                    ),
                    ("swell_self_max", format!("{:.6e}", arm.swell_self_max)),
                    (
                        "min_metric_eigenvalue",
                        format!("{:.6e}", arm.min_eigenvalue),
                    ),
                    (
                        "c1_row_swell_above_bar",
                        (arm.swell_max > C1_BAR).to_string(),
                    ),
                    (
                        "c1_componentwise_swell_max",
                        format!("{componentwise_swell_max:.6e}"),
                    ),
                    (
                        "c1_log_euclidean_swell_max",
                        format!("{log_euclidean_swell_max:.6e}"),
                    ),
                    // C2's arithmetic — M-32's own numbers on this geometry
                    ("seam_ghost_ulps", f.ghost_ulps.to_string()),
                    ("seam_own_ulps", f.own_ulps.to_string()),
                    ("seam_plane_ulps", f.plane_ulps.to_string()),
                    ("seam_position_delta", format!("{:.6e}", f.ghost_delta)),
                    (
                        "seam_position_delta_cells",
                        format!("{:.6e}", f.ghost_delta / f.g.h),
                    ),
                    ("bit_exact_pairs", arm.bit_exact_pairs.to_string()),
                    (
                        "bit_exact_first_failure_t",
                        arm.first_failure_t
                            .map_or_else(|| "none".to_string(), |t| format!("{t:.6}")),
                    ),
                    (
                        "bit_exact_seam_nondyadic",
                        arm.bit_exact_nondyadic.to_string(),
                    ),
                    ("weight_reversal_exact", weight_reversal_exact().to_string()),
                    ("dyadic_ladder_steps", T_STEPS.to_string()),
                    ("nondyadic_ladder_steps", NON_DYADIC_STEPS.to_string()),
                    ("c2_global_holds", c2_global_holds.to_string()),
                    // C3's arithmetic
                    ("seam_vertices_a", f.seam_a.to_string()),
                    ("seam_vertices_b", f.seam_b.to_string()),
                    ("seam_pairs", f.seam_pairs.to_string()),
                    ("seam_vertices_unmatched", f.unmatched.to_string()),
                    (
                        "seam_vertices_displaced",
                        census.vertices_displaced.to_string(),
                    ),
                    (
                        "seam_worst_displacement_gap",
                        format!("{:.6e}", census.worst_gap),
                    ),
                    (
                        "seam_worst_displacement_gap_cells",
                        format!("{:.6e}", census.worst_gap / f.g.h),
                    ),
                    ("weld_epsilon", format!("{:.6e}", epsilon_for(f.g.h))),
                    (
                        "gap_over_weld_epsilon",
                        format!("{:.6e}", census.worst_gap / epsilon_for(f.g.h)),
                    ),
                    (
                        "seam_open_edges_baseline",
                        f.baseline.open_edges.to_string(),
                    ),
                    ("boundary_edges_total", census.boundary_total.to_string()),
                    (
                        "boundary_edges_total_baseline",
                        f.baseline.boundary_total.to_string(),
                    ),
                    ("welded_away", census.welded_away.to_string()),
                    ("control_seam_open_edges", f.control.open_edges.to_string()),
                    (
                        "control_seam_opens",
                        (f.control.open_edges > f.baseline.open_edges).to_string(),
                    ),
                    ("c3_row_decisive", decisive.to_string()),
                    ("c3_global_holds", c3_global_holds.to_string()),
                    // the golden fixture
                    ("golden_rows", committed.len().to_string()),
                    (
                        "golden_fixture_matches_shipped",
                        golden.iter().all(|(_, g)| g.stale == 0).to_string(),
                    ),
                    (
                        "golden_control_hash_moved",
                        golden
                            .iter()
                            .all(|(_, g)| g.control_silent == 0)
                            .to_string(),
                    ),
                    ("golden_interpolations", golden_interpolations.to_string()),
                    // SHARE, priced
                    ("metric_ms", format!("{:.6}", f.metric_ms)),
                    ("interp_ms", format!("{:.6}", arm.interp_ms)),
                    ("interp_ns_per_pair", format!("{:.3}", per_scheme_ns[slot])),
                    (
                        "interp_ns_per_pair_min",
                        format!("{:.3}", per_scheme_lo[slot]),
                    ),
                    (
                        "interp_ns_per_pair_max",
                        format!("{:.3}", per_scheme_hi[slot]),
                    ),
                    (
                        "transcendental_cost_ratio",
                        format!("{transcendental_cost_ratio:.6}"),
                    ),
                    ("timing_repeats", TIMING_REPEATS.to_string()),
                ]);
            }
        }
    });
}
