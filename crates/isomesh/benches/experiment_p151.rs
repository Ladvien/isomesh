//! **P-151 — refining where the player is looking, with forty years of theory
//! behind it.**
//!
//! Ticket: R-151. Pre-registered before this harness existed.
//!
//! ```bash
//! cargo bench --bench experiment_p151
//! ```
//!
//! Writes `docs/experiments/p-151.csv`.
//!
//! # What was missing
//!
//! v1's Axis 7 named goal-oriented adaptivity and nothing in this repository has
//! ever used it. The dual-weighted-residual method (Becker & Rannacher) refines
//! to reduce error **in a chosen output functional** rather than globally, and
//! for a game the functional is screen-space error under the current camera.
//! Yano & Darmofal (`10.2514/6.2012-79`, in corpus) give the framing the
//! registration quotes: *"both the sizing decision and the anisotropy decision
//! are driven directly by the behavior of the a posteriori error estimate."*
//!
//! This row consumes `common::metric` (`benches/common/metric.rs`), which `R-146`
//! owns, for one thing only: [`principal_curvatures`], the shape operator of the
//! level set, which supplies the **primal residual** of the estimate. It does
//! **not** build an `M_Lp` metric — the adjoint is not a metric and pretending
//! otherwise would be two answers to one question.
//!
//! **The cost of that shared machinery is already measured and is quoted here
//! rather than re-derived.** `docs/experiments/p-146.csv` reports `metric_share`
//! — the Hessian-plus-`metric_lp` build as a fraction of extraction — over 40
//! rows spanning eight fields and five resolutions: **min `0.332862`, median
//! `2.008345`, max `8.747515`**, and `c3_holds` is `false` on **every one of the
//! 40 rows** against P-146's own 15% bar. In plain words: *evaluating a
//! curvature metric at every band point of the extraction grid costs one to nine
//! times the extraction it is meant to steer.* That number is the reason this
//! harness evaluates its criterion on a **fixed coarse proxy grid** and not on
//! the extraction grid, and C3 is the measurement of whether that is enough.
//!
//! `M-121` is the other prior this row is built on: *"A level change moves the
//! surface by up to 3.14 cells, which is the pop nobody measures"* —
//! `FINDINGS.md:1338`, measured as the worst vertex-to-nearest-vertex distance
//! between a block meshed at its old and its new level, **worst 3.136 cells,
//! typically 0.6–1.6**, over a full flight of an LOD ladder. `M-121`'s method
//! note is carried here verbatim in spirit: the pop *is only measurable at the
//! instant of the switch*, so a block is extracted twice rather than stored at
//! every level. `pop_magnitude_cells` below is that same statistic on this
//! harness's own block ladder.
//!
//! # The construction, and the one thing it is not
//!
//! `crates/isomesh/src/` is frozen for Phase 27, and every inherent `extract`
//! takes a single **scalar** `cell_size` (`marching_cubes/mod.rs:193` and the six
//! others). A per-cell adaptive extractor is therefore a source change by
//! construction. What is built instead is the construction the crate's LOD story
//! is actually about, and the one `M-121` measured:
//!
//! > **Block LOD.** The domain is cut into [`BLOCKS`]³ = 64 disjoint blocks. Each
//! > block is extracted as its own uniform grid at its own level, and the pieces
//! > are `MeshBuffer::append`ed. Levels are [`LEVEL_CELLS`] = 2, 4, 8, 16 cells
//! > per block axis, i.e. 3, 5, 9 and 17 **samples** per axis — all odd, and the
//! > block walls sit on the domain's own integer planes, so `y = 0` is a sample
//! > plane at every level. `M-266` makes that load-bearing: `thin_plate` is
//! > centred on `y = 0` and an even count loses its surface entirely.
//!
//! The blocks are disjoint in **cells**, so no triangle is counted twice at a
//! seam. Where two neighbours sit at different levels the seam cracks, and that
//! crack is not swept up: it is the real cost of block LOD and it is inside every
//! error number below.
//!
//! A **criterion** is a rule for turning a triangle budget into a level per
//! block. Both arms here run the *same* greedy allocator over the *same*
//! occupancy, and differ only in the error model they rank with.
//!
//! # Screen-space error: the definition, and why this one
//!
//! The registration allows two readings — project both surfaces into the camera
//! and measure pixel displacement, or measure world-space error weighted by the
//! projected pixel size `f/z`. **This harness takes the second, and takes it in
//! the surface-to-mesh direction.** Precisely:
//!
//! ```text
//!     E  =  sqrt( (1/|V|) * sum over q in V of ( f_px * d(q, M) / z(q) )^2 )
//! ```
//!
//! - `V` is a fixed set of points **on the true surface** `Γ`, generated once per
//!   field by Newton-projecting a [`PROBE_SEEDS`]³ lattice along `∇f` and kept
//!   only where the projection converged, then filtered to the **visible** ones
//!   (in frustum **and** front-facing). `V` is a property of the field and the
//!   camera and is **identical for both arms and every budget** — one instrument,
//!   or it is not a comparison.
//! - `d(q, M)` is the exact Euclidean point-to-triangle distance from `q` to the
//!   extracted mesh, over a uniform triangle grid ([`TriGrid`]).
//! - `f_px = (height_px / 2) / tan(fov_y / 2)` is the focal length in pixels, and
//!   `z(q)` is the **view-axis depth** of `q`. `f_px · δ / z` is the standard
//!   first-order pixel displacement of a world displacement `δ` transverse to the
//!   view axis; taken isotropically it is the conservative reading (a
//!   displacement *along* the ray moves no pixels), which is the same choice
//!   Lindstrom–Turk screen-space error makes.
//!
//! Three consequences, stated rather than discovered later:
//!
//! 1. **The direction is surface → mesh on purpose.** Mesh → surface (`|f(v)|`
//!    over mesh vertices) is cheaper and is what the crate's `validate::accuracy`
//!    reaches for, but it is blind to a feature the mesh *omits*: every vertex of
//!    a mesh that missed a bump is still close to the surface. `V` is fixed on
//!    the truth, so an omission is counted.
//! 2. **Visibility is `in frustum ∧ front-facing`, an occlusion proxy.** It is
//!    exact for a convex body (`sphere`, `box_exact`) and optimistic for a
//!    non-convex one (`torus`, `thin_plate`'s rim), where a front-facing patch
//!    can still be hidden behind nearer geometry. A depth buffer would settle it;
//!    a depth buffer is a renderer, and this is a bench.
//! 3. **The population is the four `FieldBound::Exact` fields.** Not because the
//!    functional needs `|f|` — it does not, it measures a real distance to real
//!    triangles — but because probe *generation* Newton-projects along `∇f` and
//!    only an exact distance field makes `p − f ∇f/‖∇f‖²` a one-step projection
//!    with a convergence test worth asserting (`fields/mod.rs:83-84`). The same
//!    narrowing P-146 documented, for a different reason, on the same roster:
//!    `sphere`, `torus`, `box_exact`, `thin_plate`.
//!
//! # The adjoint, derived
//!
//! Write the output functional as the **unnormalised, squared** form of the same
//! quantity, which is what makes it differentiable in the surface:
//!
//! ```text
//!     J(Γ_h)  =  integral over Γ of  w(q)^2 * d(q, Γ_h)^2  dA(q),
//!     w(q)    =  vis(q) * f_px / z(q).
//! ```
//!
//! DWR needs the sensitivity of `J` to a **local perturbation of the state**.
//! The state here is the extracted surface, so perturb `Γ_h` inside one block `K`
//! by `δ` along its normal. Then
//!
//! ```text
//!     dJ/dδ  =  -2 * integral over (Γ ∩ K) of w^2 * d dA,
//! ```
//!
//! and linearising `d` over the block by its mean value `ρ_K` gives the **adjoint
//! weight**
//!
//! ```text
//!     z_K  =  |dJ/dδ|  =  2 * ρ_K * integral over (Γ ∩ K) of w^2 dA
//!          =  2 * ρ_K * Z_K,          Z_K := integral over (Γ ∩ K) of w^2 dA.
//! ```
//!
//! `Z_K` is the whole of the camera's contribution and is the object this file
//! calls *the adjoint*: a per-block scalar, assembled by [`adjoint_weights`], and
//! the only thing `adjoint_ms` times. The block's share of the functional is then
//! the DWR pairing of residual against adjoint,
//!
//! ```text
//!     eta_K  =  integral over (Γ ∩ K) of w^2 d^2 dA  =  Z_K * ρ_K^2  =  (1/2) z_K ρ_K.
//! ```
//!
//! The **primal residual** `ρ_K` is the a-priori one, which is the ordinary DWR
//! practice: marching cubes reconstructs the level set by a linear interpolant,
//! so over a cell of size `h` a surface of principal curvature `κ` is missed by
//! the chord height
//!
//! ```text
//!     ρ_K(h)  =  κ_K * h^2 / 8,
//! ```
//!
//! with `κ_K` the largest `|principal curvature|` over the block's proxy probes,
//! from `common::metric::principal_curvatures` — R-146's shape operator, at the
//! proxy cell size. Promoting a block halves `h`, so `ρ` falls by 4 and `eta` by
//! **16**: the gain of one promotion is `(15/16) eta_K`.
//!
//! Discretely, with `N_K` proxy probes in the block and proxy cell size `h_p`,
//!
//! ```text
//!     Z_K  =  h_p^2 * sum over q in K of vis(q) * (f_px / z(q))^2,
//!     A_K  =  h_p^2 * N_K                        (the block's surface area),
//!     eta_goal(K, L)  =  Z_K * (κ_K * h_L^2 / 8)^2.
//! ```
//!
//! # The two criteria, and what separates them
//!
//! Both arms rank blocks by **gain per triangle** of one promotion, under the
//! same triangle model `tri(K, L) = 2 A_K / h_L^2`, and both know the same
//! occupancy `A_K`. They differ in the error model alone:
//!
//! | arm | error model | knows | `is_control` |
//! |---|---|---|---|
//! | `camera_distance` | `A_K · (f_px/r_K)² · (h_L²/8)²` | occupancy, **radial distance** `r_K` from the eye | **yes** — this is the LOD every shipped system ships |
//! | `goal_oriented` | `Z_K · (κ_K h_L²/8)²` | occupancy, **per-probe visibility and depth**, **curvature** | no — this is the row |
//!
//! The camera-distance arm is deliberately **not** a strawman. It is given the
//! same second-order residual model (`ρ ∝ h²`, i.e. `κ ≡ 1`) and the same
//! occupancy, so the only two things it lacks are **visibility** and
//! **curvature**. It is defined on radial distance rather than view depth because
//! that is what a shipped distance LOD uses — it needs no camera orientation.
//! Its ranking reduces to `(f_px/r_K)² h_L⁶`, which at equal level is exactly
//! "refine the nearest block first".
//!
//! Because the goal arm's ranking is `⟨vis · (f_px/z)²⟩_K · κ_K² · h_L⁶`, the
//! block-area factor cancels out of both rankings identically, so the crude
//! `A_K` estimator below cannot tilt the comparison.
//!
//! **Curvature is the confound and C2 is the instrument that names it.** The two
//! arms differ in *two* things, not one, and the registration's C2 asks precisely
//! whether the difference that matters is visibility. That is what the second
//! camera is for.
//!
//! # Arms
//!
//! | arm | what it varies | `is_control` |
//! |---|---|---|
//! | camera `near` | eye at [`NEAR_EYE`], inside the domain, most of the scene outside the frustum | no — this is the registered scene |
//! | camera `wide` | same direction, eye at [`WIDE_RANGE`] units, **the whole domain inside the frustum** | **yes** — C2's control |
//! | criterion `camera_distance` | distance LOD | **yes** — C1's baseline |
//! | criterion `goal_oriented` | the DWR estimate above | no |
//! | budget rung 0..3 | [`BUDGET_FRACTIONS`] of the all-at-finest triangle estimate | — |
//!
//! Four fields x two cameras x two criteria x four rungs = **64 rows**.
//!
//! # The three clauses, and how each is decided
//!
//! **C1 — at matched triangle count, goal-oriented beats camera-distance on
//! screen-space error by at least 20%.** The matched count is
//! `T* = triangles(camera_distance, rung `[`REFERENCE_RUNG`]`)` on the **near**
//! camera. Each arm's `screen_error_at_matched_triangles` is read off its own
//! measured `(triangles, error)` curve by linear interpolation in
//! `(ln T, ln E)`, **clamped to the curve's endpoints**. Clamping matters and its
//! direction is deliberate: the goal arm can *saturate* — a block whose surface
//! is exactly flat has `κ_K = 0`, an honest zero gain, and the greedy will not
//! spend triangles on it — so its curve may stop short of `T*`. Clamping then
//! charges it the error it had at **fewer** triangles, which is larger. The
//! clamp is therefore conservative *against* the hypothesis, and every clamped
//! read-off is flagged on the row.
//!
//! `c1_holds` is the field's near-camera verdict, `1 − E_goal(T*)/E_cam(T*) ≥
//! 0.20`, and it is written identically on all sixteen of that field's rows.
//!
//! **C2 — the saving comes from not refining what is off-screen or back-facing,
//! rather than from a general improvement.** Its falsifier is *"a saving that
//! persists with everything on-screen"*. So the same field is measured a second
//! time from the `wide` camera, from which **no part of the domain is outside the
//! frustum**, and
//!
//! ```text
//!     c2_holds  =  gain_near >= 0.20  AND  gain_wide < 0.5 * gain_near.
//! ```
//!
//! `off_screen_triangles` — the registered column — is the count of extracted
//! triangles whose centroid is outside the frustum **or** whose geometric normal
//! faces away from the eye, which is the union C2 names. Its two halves are
//! recorded separately, because they behave differently and pooling them hides
//! the mechanism: **moving the camera back removes frustum invisibility and
//! cannot remove back-facing.** Roughly half of a closed surface faces away from
//! any eye, at any distance. So the `wide` control isolates the **frustum** half
//! of the mechanism and leaves the back-face half in both arms. A reader must
//! take C2's verdict as a statement about frustum visibility; `back_facing_
//! fraction` is on every row so the residual is visible rather than implied.
//!
//! **C3 — the adjoint cost is under 10% of extraction, or the method is a
//! non-starter for a frame budget and should be said so.** What is timed is
//! [`adjoint_weights`] and nothing else: per proxy probe, one `Sdf::gradient`,
//! one normalisation, one back-face test, one frustum test, one `f_px/z`, one
//! accumulate. Against it is timed the block-LOD extraction of the goal arm at
//! the reference rung. Seven repeats, **interleaved** (adjoint, extract, adjoint,
//! extract, …) after one warm-up, median as the headline, min and max and the
//! max/min scatter as extras — `M-280` measured this host's `amd-pstate-epp`
//! governor swinging the same binary 1.45x between runs, and an interleaved
//! schedule is the only way a *ratio* survives a drifting clock.
//!
//! **The adjoint is timed with nothing cached, and that is a choice.** A cached
//! surface normal would make the per-frame adjoint a handful of flops and C3 a
//! foregone conclusion; in a destructible-terrain game the field is edited
//! between frames and a cached normal is stale. So the gradient is recomputed at
//! every probe, every repeat.
//!
//! **C3 is about the adjoint, but the criterion is not only the adjoint, and the
//! difference is quoted from P-146 rather than hidden.** Two further stages are
//! timed and reported beside it: `proxy_ms` (sampling the proxy grid, finding the
//! surface-proximate cells, projecting one probe into each) and `curvature_ms`
//! (`principal_curvatures` — nineteen samples and a Jacobi sweep per probe, the
//! machinery whose share on the *extraction* grid P-146 measured at 0.33x–8.75x
//! extraction). `criterion_ms` is their sum with the adjoint and
//! `criterion_share` its share of extraction. A reader who takes `adjoint_share`
//! as the cost of the method has been warned on the row.
//!
//! # SHARE, recomputed before the numbers
//!
//! The registration's SHARE is *"C1 moves the LOD selection stage; `M-121`
//! measured a level change moving the surface by up to 3.14 cells, which is the
//! pop this criterion is meant to spend better."*
//!
//! Discharged, and the discharge is a column rather than a sentence.
//! `pop_magnitude_cells` is `M-121`'s statistic on this ladder: for every block
//! the criterion **promoted**, the block is extracted at both its old and its new
//! level and the two vertex sets' two-sided nearest-vertex Hausdorff distance is
//! taken, in units of the **coarse** level's cell size; the row reports the worst
//! over the blocks it promoted. That is the pop a player would see at the instant
//! of the switch.
//!
//! But "spend the pop better" is not "make the pop smaller" — a coarse mesh
//! genuinely *is* a different surface, and `M-121` says so. The quantity that
//! decides whether a pop is *visible* is its size **in pixels**, so
//! `pop_magnitude_pixels` = `pop_world · f_px / z_K`, and **zero for a block with
//! no visible probe**. A goal-oriented criterion that has done its job has a
//! comparable `pop_magnitude_cells` and a smaller `pop_magnitude_pixels`: the
//! same pops, moved off the screen. That is the SHARE, and it is falsifiable on
//! the row rather than asserted in prose.
//!
//! # Vacuity controls
//!
//! Every one runs before the first `run.record` and every panic starts
//! `"VOID: "`. `M-44`: a zero that could not have been non-zero is not a
//! measurement.
//!
//! - **The registered control.** *"The camera must be positioned so that a
//!   non-trivial fraction of the scene is off-screen, reported as a fraction, or
//!   C2 cannot fire."* Reported as `off_screen_fraction` — the fraction of the
//!   field's true-surface probes **outside the frustum** — and asserted above
//!   [`OFF_SCREEN_FLOOR`] on the `near` camera for every field. It is asserted on
//!   the *frustum* fraction and not on the registered `off_screen_triangles`
//!   union on purpose: back-facing alone puts about half of any closed surface
//!   off-screen for free, so a union above 20% is a fact about closed surfaces
//!   and not about where the camera was put.
//! - **The C2 control must be a control.** The `wide` camera's
//!   `off_screen_fraction` must be under [`WIDE_ON_SCREEN_CEIL`] on every field,
//!   or "everything on-screen" is not what was run and C2's falsifier has no
//!   instrument. Column: `off_screen_fraction`.
//! - **The functional must have a population.** `visible_probes > 0` on every
//!   (field, camera), or `screen_space_error` is an RMS over nothing.
//! - **The two criteria must disagree.** At least one (field, camera, rung) must
//!   assign different levels under the two arms, or every `ratio` in the file is
//!   a measurement of the extractor's determinism. Column: `levels_differ`.
//! - **Something must be promoted.** At least one row must promote at least one
//!   block, or `pop_magnitude_cells` is a zero that could not have been
//!   non-zero. Column: `blocks_promoted`.
//!
//! # Determinism
//!
//! One thread, no PRNG, no map iteration, `f64` throughout. The greedy scans
//! blocks in index order and takes a strictly-greater comparison, so ties go to
//! the lowest block index. Probes are generated by sweeping `z`, `y`, `x` with
//! `x` innermost, the crate's order. Sorting is [`f64::total_cmp`]. Wall-clock
//! columns are the only machine-dependent numbers and only C3 reads them.

#![allow(
    clippy::cast_possible_truncation,
    clippy::cast_precision_loss,
    clippy::cast_sign_loss,
    clippy::too_many_arguments,
    clippy::too_many_lines
)]

mod common;

use std::hint::black_box;
use std::time::Instant;

use isomesh::fields::{BoxExact, ReferenceField, Sphere, ThinPlate, Torus};
use isomesh::marching_cubes::MarchingCubes;
use isomesh::{MeshBuffer, RuntimeShape3, Sdf};

use crate::common::metric::principal_curvatures;

// ─── the fixture ─────────────────────────────────────────────────────────────

/// Blocks per axis. Sixty-four blocks over a `[-2, 2]³` domain, so a block is a
/// unit cube and every block wall lands on an integer plane — which is what puts
/// `y = 0` on a sample plane at every level (`M-266`).
const BLOCKS: usize = 4;

/// Cells per block axis at each level. Samples per axis are these plus one:
/// **3, 5, 9, 17** — all odd, all the time.
const LEVEL_CELLS: [u32; 4] = [2, 4, 8, 16];

/// The finest level's index.
const MAX_LEVEL: usize = LEVEL_CELLS.len() - 1;

/// Proxy cells per block axis. Five samples per block axis, a `16³` grid over the
/// whole domain — coarse on purpose: `docs/experiments/p-146.csv` measured the
/// curvature machinery at 0.33x–8.75x extraction when it is evaluated at every
/// band point of the *extraction* grid, and C3 is the question of whether a
/// criterion can be afforded at all.
const PROXY_CELLS: u32 = 4;

/// Triangle budgets, as fractions of the all-blocks-at-[`MAX_LEVEL`] estimate.
///
/// The all-at-level-0 floor is `(h_max/h_0)² = 1/64 = 0.015625` of that estimate,
/// so every rung is reachable from the starting assignment.
const BUDGET_FRACTIONS: [f64; 4] = [0.06, 0.15, 0.35, 0.75];

/// Which rung fixes the matched triangle count `T*`.
const REFERENCE_RUNG: usize = 2;

/// Estimated triangles per surface cell. Cancels from both rankings; it is here
/// only so the budget is denominated in something a reader recognises.
const TRI_PER_CELL: f64 = 2.0;

/// Seeds per axis for the true-surface probe lattice.
const PROBE_SEEDS: u32 = 33;

/// Newton iterations for the projection onto the level set.
const NEWTON_STEPS: usize = 8;

/// A seed is projected only if it already lies this close, in units of the seed
/// lattice spacing. Half a cell: the surface passes within half a spacing of
/// exactly the seeds that bracket it.
const SEED_BAND_CELLS: f64 = 0.5;

/// A projected seed counts as a surface probe when `|f|` falls below this.
const NEWTON_TOL: f64 = 1e-9;

/// Below this gradient magnitude the Newton step has no direction and the seed is
/// dropped.
const GRAD_FLOOR: f64 = 1e-12;

// ─── the camera ──────────────────────────────────────────────────────────────

/// Frame width in pixels.
const PIXELS_X: f64 = 1920.0;

/// Frame height in pixels.
const PIXELS_Y: f64 = 1080.0;

/// Vertical field of view, degrees.
const FOV_Y_DEG: f64 = 38.0;

/// The registered camera's eye. Inside the domain and close to the surface, so a
/// large part of every field falls outside the frustum — which is the
/// registration's own vacuity control.
const NEAR_EYE: [f64; 3] = [0.30, 0.62, 1.30];

/// Both cameras look here.
const TARGET: [f64; 3] = [0.0, 0.0, 0.0];

/// The C2 control camera's range from [`TARGET`], along the same direction. The
/// domain's corner radius is `2√3 = 3.46`, and the vertical half-window at the
/// nearest domain depth is `(20 − 3.46) · tan(19°) = 5.70`, so the whole domain
/// is inside the frustum with room to spare.
const WIDE_RANGE: f64 = 20.0;

/// Near plane. A probe closer than this has no defined pixel size.
const NEAR_PLANE: f64 = 0.01;

/// Far plane. Large enough that the `wide` camera clips nothing.
const FAR_PLANE: f64 = 1000.0;

// ─── the bars ────────────────────────────────────────────────────────────────

/// C1's bar: "beats camera-distance LOD on screen-space error by at least 20%".
const C1_BAR: f64 = 0.20;

/// C2's bar: the saving must not *persist* with everything on-screen. The wide
/// camera's gain must be under this fraction of the near camera's.
const C2_PERSIST: f64 = 0.5;

/// C3's bar: "the adjoint cost is under 10% of extraction".
const C3_BAR: f64 = 0.10;

/// The registered vacuity control's bar on the off-screen fraction.
const OFF_SCREEN_FLOOR: f64 = 0.20;

/// The C2 control camera must have essentially nothing outside the frustum.
const WIDE_ON_SCREEN_CEIL: f64 = 0.02;

/// Timed repeats per stage, interleaved, after one warm-up.
const REPEATS: usize = 7;

/// Floor applied to an error before it is logged for the matched read-off. One
/// ten-thousandth of a pixel is below any display's ability to show it.
const ERROR_FLOOR_PX: f64 = 1e-4;

/// `metric_share` over all 40 rows of `docs/experiments/p-146.csv` — the measured
/// cost of the shared curvature machinery on the *extraction* grid, quoted rather
/// than re-derived. P-146's `c3_holds` is `false` on every one of those rows.
const P146_METRIC_SHARE: [f64; 3] = [0.332_862, 2.008_345, 8.747_515];

/// The name of the output functional, for the `functional` column.
const FUNCTIONAL: &str = "rms_visible_pixel_displacement";

// ─── small vector algebra ────────────────────────────────────────────────────

fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn add(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

fn scale(a: [f64; 3], k: f64) -> [f64; 3] {
    [a[0] * k, a[1] * k, a[2] * k]
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn norm(a: [f64; 3]) -> f64 {
    dot(a, a).sqrt()
}

/// Unit vector. Panics on a zero input: every call site here is fed a difference
/// of two distinct fixture constants or a gradient already tested against
/// [`GRAD_FLOOR`], so a zero means the caller is wrong rather than the input.
fn unit(a: [f64; 3]) -> [f64; 3] {
    let len = norm(a);
    assert!(len > 0.0, "unit vector of a zero vector: {a:?}");
    scale(a, 1.0 / len)
}

fn dist_sq(a: [f64; 3], b: [f64; 3]) -> f64 {
    let d = sub(a, b);
    dot(d, d)
}

// ─── the camera ──────────────────────────────────────────────────────────────

/// A pinhole camera: an eye, an orthonormal frame, and the two half-angle
/// tangents of its frustum.
#[derive(Clone, Copy)]
struct Camera {
    name: &'static str,
    eye: [f64; 3],
    forward: [f64; 3],
    right: [f64; 3],
    up: [f64; 3],
    /// `(height_px / 2) / tan(fov_y / 2)`: pixels per radian at the centre, and
    /// therefore pixels per world unit at unit depth.
    f_px: f64,
    tan_x: f64,
    tan_y: f64,
}

impl Camera {
    fn new(name: &'static str, eye: [f64; 3], target: [f64; 3]) -> Self {
        let forward = unit(sub(target, eye));
        let right = unit(cross(forward, [0.0, 1.0, 0.0]));
        let up = cross(right, forward);
        let tan_y = (FOV_Y_DEG.to_radians() * 0.5).tan();
        Self {
            name,
            eye,
            forward,
            right,
            up,
            f_px: PIXELS_Y * 0.5 / tan_y,
            tan_x: tan_y * (PIXELS_X / PIXELS_Y),
            tan_y,
        }
    }

    /// View-axis depth of a world point.
    fn depth(&self, p: [f64; 3]) -> f64 {
        dot(sub(p, self.eye), self.forward)
    }

    /// Radial distance from the eye — what a shipped distance LOD uses, because
    /// it needs no camera orientation.
    fn range(&self, p: [f64; 3]) -> f64 {
        norm(sub(p, self.eye))
    }

    fn in_frustum(&self, p: [f64; 3]) -> bool {
        let d = sub(p, self.eye);
        let z = dot(d, self.forward);
        (NEAR_PLANE..=FAR_PLANE).contains(&z)
            && dot(d, self.right).abs() <= z * self.tan_x
            && dot(d, self.up).abs() <= z * self.tan_y
    }

    /// `true` when a patch at `p` with outward normal `n` faces the eye.
    fn front_facing(&self, p: [f64; 3], n: [f64; 3]) -> bool {
        dot(n, sub(p, self.eye)) < 0.0
    }

    /// Pixels per world unit of transverse displacement at `p`.
    fn pixels_per_world(&self, p: [f64; 3]) -> f64 {
        self.f_px / self.depth(p).max(NEAR_PLANE)
    }
}

// ─── true-surface probes ─────────────────────────────────────────────────────

/// One point of the true surface, with its outward unit normal.
#[derive(Clone, Copy)]
struct Probe {
    p: [f64; 3],
    n: [f64; 3],
}

/// A probe the camera can see, carrying its own pixel scale.
#[derive(Clone, Copy)]
struct VisibleProbe {
    p: [f64; 3],
    px_per_world: f64,
}

/// Newton-project a [`PROBE_SEEDS`]³ lattice onto the level set.
///
/// The step is `p ← p − f ∇f / ‖∇f‖²`, which for an exact distance field lands on
/// the surface in one iteration on a flat patch and converges quadratically
/// elsewhere. A seed is admitted only if it starts within [`SEED_BAND_CELLS`] of
/// the surface — by the 1-Lipschitz property of an exact distance field that is
/// exactly the seeds whose cell the surface can reach — and is kept only if it
/// converges to [`NEWTON_TOL`] and stays inside the domain.
///
/// Swept `z`, `y`, `x` with `x` innermost, the crate's order, so the probe list
/// is the same list every run.
fn surface_probes<F>(field: &F, lo: [f64; 3], hi: [f64; 3]) -> Vec<Probe>
where
    F: Sdf<Scalar = f64>,
{
    let spacing = (hi[0] - lo[0]) / f64::from(PROBE_SEEDS - 1);
    let band = SEED_BAND_CELLS * spacing * 3.0f64.sqrt();
    let mut probes = Vec::new();

    for kz in 0..PROBE_SEEDS {
        for ky in 0..PROBE_SEEDS {
            for kx in 0..PROBE_SEEDS {
                let seed = [
                    lo[0] + f64::from(kx) * spacing,
                    lo[1] + f64::from(ky) * spacing,
                    lo[2] + f64::from(kz) * spacing,
                ];
                if field.sample(seed).abs() > band {
                    continue;
                }
                let mut p = seed;
                let mut grad = [0.0f64; 3];
                for _ in 0..NEWTON_STEPS {
                    let value = field.sample(p);
                    grad = field.gradient(p);
                    let gg = dot(grad, grad);
                    if gg <= GRAD_FLOOR {
                        break;
                    }
                    p = sub(p, scale(grad, value / gg));
                }
                if field.sample(p).abs() > NEWTON_TOL {
                    continue;
                }
                if (0..3).any(|axis| p[axis] < lo[axis] || p[axis] > hi[axis]) {
                    continue;
                }
                if norm(grad) <= GRAD_FLOOR {
                    continue;
                }
                probes.push(Probe {
                    p,
                    n: unit(field.gradient(p)),
                });
            }
        }
    }
    probes
}

// ─── the proxy geometry ──────────────────────────────────────────────────────

/// One block's camera-independent geometry, read off the coarse proxy grid.
struct BlockProxy {
    /// The block's minimum corner.
    lo: [f64; 3],
    /// One probe per surface-proximate proxy cell, at the cell centre projected
    /// one Newton step onto the level set.
    probes: Vec<[f64; 3]>,
    /// `h_p² · probes.len()` — the block's surface area, up to the constant
    /// factor by which a proximity count overestimates a crossing count. That
    /// factor cancels out of both rankings (see the header), so it is left where
    /// it is rather than fitted.
    area: f64,
    /// Largest `|principal curvature|` over the block's probes.
    kappa: f64,
}

/// The whole proxy: sixty-four blocks and the cell size they were read at.
struct Proxy {
    blocks: Vec<BlockProxy>,
    /// Block edge length.
    block_size: f64,
    /// Proxy cell size.
    h: f64,
}

impl Proxy {
    /// The cell size of one block at one level.
    fn cell_size(&self, level: usize) -> f64 {
        self.block_size / f64::from(LEVEL_CELLS[level])
    }

    /// The block indices that carry surface, in index order.
    fn occupied(&self) -> Vec<usize> {
        (0..self.blocks.len())
            .filter(|&k| !self.blocks[k].probes.is_empty())
            .collect()
    }
}

/// Sample the proxy grid, find the surface-proximate cells, project one probe
/// into each. Timed as `proxy_ms`; carries no camera and no curvature.
///
/// A cell is surface-proximate when the smallest `|f|` over its eight corners is
/// at most half the cell diagonal. For an exact distance field that test is
/// **provably conservative**: if the surface reaches any point of the cell, the
/// nearest corner is within half a diagonal of it and `|f|` there is at most that
/// distance. It is also generous — a shell about two cells thick rather than one
/// — and the header says why that does not tilt anything.
fn build_proxy<F>(field: &F, lo: [f64; 3], hi: [f64; 3]) -> Proxy
where
    F: Sdf<Scalar = f64>,
{
    let block_size = (hi[0] - lo[0]) / BLOCKS as f64;
    let h = block_size / f64::from(PROXY_CELLS);
    let half_diagonal = h * 3.0f64.sqrt() * 0.5;

    let mut blocks = Vec::with_capacity(BLOCKS * BLOCKS * BLOCKS);
    for bz in 0..BLOCKS {
        for by in 0..BLOCKS {
            for bx in 0..BLOCKS {
                let block_lo = [
                    lo[0] + bx as f64 * block_size,
                    lo[1] + by as f64 * block_size,
                    lo[2] + bz as f64 * block_size,
                ];
                let mut probes = Vec::new();
                for cz in 0..PROXY_CELLS {
                    for cy in 0..PROXY_CELLS {
                        for cx in 0..PROXY_CELLS {
                            let cell_lo = [
                                block_lo[0] + f64::from(cx) * h,
                                block_lo[1] + f64::from(cy) * h,
                                block_lo[2] + f64::from(cz) * h,
                            ];
                            let mut nearest = f64::INFINITY;
                            for corner in 0..8u32 {
                                let c = [
                                    cell_lo[0] + f64::from(corner & 1) * h,
                                    cell_lo[1] + f64::from((corner >> 1) & 1) * h,
                                    cell_lo[2] + f64::from((corner >> 2) & 1) * h,
                                ];
                                let value = field.sample(c).abs();
                                if value < nearest {
                                    nearest = value;
                                }
                            }
                            if nearest > half_diagonal {
                                continue;
                            }
                            let centre = [
                                cell_lo[0] + 0.5 * h,
                                cell_lo[1] + 0.5 * h,
                                cell_lo[2] + 0.5 * h,
                            ];
                            let value = field.sample(centre);
                            let grad = field.gradient(centre);
                            let gg = dot(grad, grad);
                            if gg <= GRAD_FLOOR {
                                continue;
                            }
                            probes.push(sub(centre, scale(grad, value / gg)));
                        }
                    }
                }
                let area = h * h * probes.len() as f64;
                blocks.push(BlockProxy {
                    lo: block_lo,
                    probes,
                    area,
                    kappa: 0.0,
                });
            }
        }
    }
    Proxy {
        blocks,
        block_size,
        h,
    }
}

/// The primal residual's curvature, per block: the largest `|principal
/// curvature|` over the block's probes, from `common::metric`. Timed as
/// `curvature_ms`.
///
/// Nineteen field samples and one Jacobi sweep per probe — this is the stage
/// whose share on the *extraction* grid P-146 measured at 0.33x–8.75x extraction
/// (`docs/experiments/p-146.csv`, `metric_share`, 40 rows).
fn block_curvatures<F>(field: &F, proxy: &Proxy) -> Vec<f64>
where
    F: Sdf<Scalar = f64>,
{
    proxy
        .blocks
        .iter()
        .map(|block| {
            let mut worst = 0.0f64;
            for &q in &block.probes {
                if let Some(kappas) = principal_curvatures(field, q, proxy.h) {
                    let magnitude = kappas[0].abs().max(kappas[1].abs());
                    if magnitude > worst {
                        worst = magnitude;
                    }
                }
            }
            worst
        })
        .collect()
}

// ─── the adjoint ─────────────────────────────────────────────────────────────

/// **The adjoint.** `Z_K = h_p² Σ_q vis(q) (f_px / z(q))²` per block — the
/// sensitivity of the screen-space functional to a unit normal perturbation of
/// the surface inside the block, up to the factor `2 ρ_K` derived in the header.
///
/// This function and nothing else is what `adjoint_ms` times. Nothing is cached:
/// the surface normal is recomputed from the field at every probe on every call,
/// because in a destructible-terrain game a cached normal is stale by the next
/// frame, and C3 is a claim about a frame budget.
fn adjoint_weights<F>(field: &F, camera: &Camera, proxy: &Proxy, out: &mut Vec<f64>)
where
    F: Sdf<Scalar = f64>,
{
    out.clear();
    let cell_area = proxy.h * proxy.h;
    for block in &proxy.blocks {
        let mut acc = 0.0f64;
        for &q in &block.probes {
            let grad = field.gradient(q);
            let length = norm(grad);
            if length <= GRAD_FLOOR {
                continue;
            }
            let n = scale(grad, 1.0 / length);
            if !camera.front_facing(q, n) || !camera.in_frustum(q) {
                continue;
            }
            let w = camera.f_px / camera.depth(q);
            acc += w * w;
        }
        out.push(acc * cell_area);
    }
}

/// Mean radial distance from the eye, per block — the camera-distance arm's only
/// camera input. Blocks with no probes report `INFINITY`, which the allocator
/// never reaches because it never considers them.
fn block_ranges(camera: &Camera, proxy: &Proxy) -> Vec<f64> {
    proxy
        .blocks
        .iter()
        .map(|block| {
            if block.probes.is_empty() {
                return f64::INFINITY;
            }
            let total: f64 = block.probes.iter().map(|&q| camera.range(q)).sum();
            total / block.probes.len() as f64
        })
        .collect()
}

// ─── the allocator ───────────────────────────────────────────────────────────

/// Which error model a criterion ranks with.
#[derive(Clone, Copy, PartialEq, Eq)]
enum Criterion {
    /// Occupancy and radial distance only: `A_K (f_px/r_K)² (h_L²/8)²`.
    CameraDistance,
    /// The DWR estimate: `Z_K (κ_K h_L²/8)²`.
    GoalOriented,
}

impl Criterion {
    fn name(self) -> &'static str {
        match self {
            Self::CameraDistance => "camera_distance",
            Self::GoalOriented => "goal_oriented",
        }
    }
}

/// `eta[K][L]`: the block's modelled share of the functional at each level, under
/// one criterion.
///
/// Both models are `weight · ρ(h_L)²` with `ρ = κ h²/8`. They differ in exactly
/// two places, which is the whole comparison: the camera-distance model replaces
/// the per-probe visible weight `Z_K` by the blind `A_K (f_px/r_K)²`, and
/// replaces the measured curvature by `κ ≡ 1`.
fn energies(
    criterion: Criterion,
    proxy: &Proxy,
    camera: &Camera,
    adjoint: &[f64],
    ranges: &[f64],
) -> Vec<[f64; 4]> {
    (0..proxy.blocks.len())
        .map(|k| {
            let block = &proxy.blocks[k];
            let (weight, kappa) = match criterion {
                Criterion::GoalOriented => (adjoint[k], block.kappa),
                Criterion::CameraDistance => {
                    let px = camera.f_px / ranges[k];
                    (block.area * px * px, 1.0)
                }
            };
            let mut out = [0.0f64; 4];
            for (level, slot) in out.iter_mut().enumerate() {
                let h = proxy.cell_size(level);
                let residual = kappa * h * h / 8.0;
                *slot = weight * residual * residual;
            }
            out
        })
        .collect()
}

/// What one greedy allocation produced.
struct Allocation {
    levels: Vec<usize>,
    /// Modelled triangle count at the end.
    estimated: f64,
    /// `true` when the budget stopped the greedy; `false` when it ran out of
    /// positive-gain promotions first. A goal-oriented arm saturates on a field
    /// with flat blocks and that is an honest zero, not a defect.
    budget_bound: bool,
}

/// Greedy: promote the block with the largest gain per triangle whose promotion
/// still fits the budget, until none does.
///
/// Deterministic: the scan is in block index order and the comparison is strictly
/// greater, so a tie goes to the lowest index.
fn allocate(proxy: &Proxy, eta: &[[f64; 4]], occupied: &[usize], target: f64) -> Allocation {
    let mut levels = vec![0usize; proxy.blocks.len()];
    let tri = |k: usize, level: usize| {
        TRI_PER_CELL * proxy.blocks[k].area / (proxy.cell_size(level) * proxy.cell_size(level))
    };
    let mut total: f64 = occupied.iter().map(|&k| tri(k, 0)).sum();

    let budget_bound = loop {
        let mut best: Option<(f64, usize, f64)> = None;
        let mut blocked_by_budget = false;
        for &k in occupied {
            let level = levels[k];
            if level == MAX_LEVEL {
                continue;
            }
            let gain = eta[k][level] - eta[k][level + 1];
            let cost = tri(k, level + 1) - tri(k, level);
            if gain <= 0.0 || cost <= 0.0 {
                continue;
            }
            if total + cost > target {
                blocked_by_budget = true;
                continue;
            }
            let score = gain / cost;
            if best.is_none_or(|(top, _, _)| score > top) {
                best = Some((score, k, cost));
            }
        }
        match best {
            Some((_, k, cost)) => {
                levels[k] += 1;
                total += cost;
            }
            None => break blocked_by_budget,
        }
    };

    Allocation {
        levels,
        estimated: total,
        budget_bound,
    }
}

// ─── extraction ──────────────────────────────────────────────────────────────

/// Extract every occupied block at its own level and append the pieces.
///
/// The blocks are disjoint in cells, so nothing is counted twice at a seam; where
/// two neighbours differ in level the seam cracks, and that crack is inside every
/// error number this file reports.
fn extract_levels<F>(
    field: &F,
    proxy: &Proxy,
    occupied: &[usize],
    levels: &[usize],
    shapes: &[RuntimeShape3; 4],
    mc: &mut MarchingCubes<f64>,
    scratch: &mut MeshBuffer<f64>,
    out: &mut MeshBuffer<f64>,
) where
    F: Sdf<Scalar = f64>,
{
    out.reset();
    for &k in occupied {
        let level = levels[k];
        scratch.reset();
        mc.extract(
            field,
            &shapes[level],
            proxy.blocks[k].lo,
            proxy.cell_size(level),
            scratch,
        )
        .expect("P-151: marching cubes on a block grid");
        out.append(scratch)
            .expect("P-151: block mesh fits the index space");
    }
}

// ─── the error instrument ────────────────────────────────────────────────────

/// A uniform grid over the extracted triangles, for exact point-to-triangle
/// nearest-distance queries.
///
/// Degenerate (zero-area) triangles are excluded and counted: a triangle with no
/// area is not a surface, and the barycentric branch of the closest-point
/// routine divides by its area.
struct TriGrid {
    lo: [f64; 3],
    cell: f64,
    dim: [i64; 3],
    bins: Vec<Vec<u32>>,
    tris: Vec<[[f64; 3]; 3]>,
    degenerate: usize,
}

impl TriGrid {
    fn new(mesh: &MeshBuffer<f64>, lo: [f64; 3], hi: [f64; 3], cell: f64) -> Self {
        let dim = [
            (((hi[0] - lo[0]) / cell).ceil() as i64).max(1),
            (((hi[1] - lo[1]) / cell).ceil() as i64).max(1),
            (((hi[2] - lo[2]) / cell).ceil() as i64).max(1),
        ];
        let mut grid = Self {
            lo,
            cell,
            dim,
            bins: vec![Vec::new(); (dim[0] * dim[1] * dim[2]) as usize],
            tris: Vec::with_capacity(mesh.triangle_count()),
            degenerate: 0,
        };
        for face in mesh.indices.as_chunks::<3>().0 {
            let t = [
                mesh.positions[face[0] as usize],
                mesh.positions[face[1] as usize],
                mesh.positions[face[2] as usize],
            ];
            if norm(cross(sub(t[1], t[0]), sub(t[2], t[0]))) <= 0.0 {
                grid.degenerate += 1;
                continue;
            }
            let index = grid.tris.len() as u32;
            grid.tris.push(t);
            let mut bin_lo = [0i64; 3];
            let mut bin_hi = [0i64; 3];
            for axis in 0..3 {
                let a = t[0][axis].min(t[1][axis]).min(t[2][axis]);
                let b = t[0][axis].max(t[1][axis]).max(t[2][axis]);
                bin_lo[axis] = grid.clamp_index(a, axis);
                bin_hi[axis] = grid.clamp_index(b, axis);
            }
            for z in bin_lo[2]..=bin_hi[2] {
                for y in bin_lo[1]..=bin_hi[1] {
                    for x in bin_lo[0]..=bin_hi[0] {
                        let slot = (x + grid.dim[0] * (y + grid.dim[1] * z)) as usize;
                        grid.bins[slot].push(index);
                    }
                }
            }
        }
        grid
    }

    fn clamp_index(&self, value: f64, axis: usize) -> i64 {
        let raw = ((value - self.lo[axis]) / self.cell).floor() as i64;
        raw.clamp(0, self.dim[axis] - 1)
    }

    /// Distance from `p` to the nearest triangle.
    ///
    /// Rings of bins are searched outward. After all rings up to `r` are done,
    /// every unsearched bin starts at least `r · cell` from `p`, so the search
    /// stops as soon as the best distance is at or under that. Triangles are
    /// binned by their bounding box, so a triangle that comes near `p` is in a
    /// near bin — the bound is sound. The loop terminates unconditionally at the
    /// grid's own extent, where every bin has been visited.
    fn nearest(&self, p: [f64; 3]) -> f64 {
        let base = [
            self.clamp_index(p[0], 0),
            self.clamp_index(p[1], 1),
            self.clamp_index(p[2], 2),
        ];
        let max_ring = self.dim[0].max(self.dim[1]).max(self.dim[2]);
        let mut best = f64::INFINITY;
        for ring in 0..=max_ring {
            if best.is_finite() && best.sqrt() <= ring as f64 * self.cell {
                break;
            }
            for z in (base[2] - ring)..=(base[2] + ring) {
                if z < 0 || z >= self.dim[2] {
                    continue;
                }
                for y in (base[1] - ring)..=(base[1] + ring) {
                    if y < 0 || y >= self.dim[1] {
                        continue;
                    }
                    for x in (base[0] - ring)..=(base[0] + ring) {
                        if x < 0 || x >= self.dim[0] {
                            continue;
                        }
                        let on_shell = (x - base[0]).abs() == ring
                            || (y - base[1]).abs() == ring
                            || (z - base[2]).abs() == ring;
                        if !on_shell {
                            continue;
                        }
                        let slot = (x + self.dim[0] * (y + self.dim[1] * z)) as usize;
                        for &index in &self.bins[slot] {
                            let t = self.tris[index as usize];
                            let d = point_triangle_distance_sq(p, t[0], t[1], t[2]);
                            if d < best {
                                best = d;
                            }
                        }
                    }
                }
            }
        }
        best.sqrt()
    }
}

/// Squared distance from `p` to the triangle `abc`, by the Voronoi-region
/// classification of Ericson, *Real-Time Collision Detection* §5.1.5.
///
/// The final barycentric branch divides by `va + vb + vc`, which is positive
/// whenever the triangle has area and `p` projects inside it; zero-area triangles
/// never reach it because they are excluded from [`TriGrid`].
fn point_triangle_distance_sq(p: [f64; 3], a: [f64; 3], b: [f64; 3], c: [f64; 3]) -> f64 {
    let ab = sub(b, a);
    let ac = sub(c, a);

    let ap = sub(p, a);
    let d1 = dot(ab, ap);
    let d2 = dot(ac, ap);
    if d1 <= 0.0 && d2 <= 0.0 {
        return dist_sq(p, a);
    }

    let bp = sub(p, b);
    let d3 = dot(ab, bp);
    let d4 = dot(ac, bp);
    if d3 >= 0.0 && d4 <= d3 {
        return dist_sq(p, b);
    }

    let vc = d1 * d4 - d3 * d2;
    if vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0 {
        let v = d1 / (d1 - d3);
        return dist_sq(p, add(a, scale(ab, v)));
    }

    let cp = sub(p, c);
    let d5 = dot(ab, cp);
    let d6 = dot(ac, cp);
    if d6 >= 0.0 && d5 <= d6 {
        return dist_sq(p, c);
    }

    let vb = d5 * d2 - d1 * d6;
    if vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0 {
        let w = d2 / (d2 - d6);
        return dist_sq(p, add(a, scale(ac, w)));
    }

    let va = d3 * d6 - d5 * d4;
    if va <= 0.0 && (d4 - d3) >= 0.0 && (d5 - d6) >= 0.0 {
        let w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
        return dist_sq(p, add(b, scale(sub(c, b), w)));
    }

    let denom = 1.0 / (va + vb + vc);
    let v = vb * denom;
    let w = vc * denom;
    dist_sq(p, add(a, add(scale(ab, v), scale(ac, w))))
}

/// What the functional read off one mesh.
struct ErrorReading {
    /// `E` — RMS visible pixel displacement.
    rms_px: f64,
    /// The worst single visible probe, in pixels.
    max_px: f64,
    /// The same population in world units, unweighted — the quantity the pixel
    /// figure is `f/z` times, kept so a reader can separate geometry from camera.
    rms_world: f64,
}

fn screen_error(grid: &TriGrid, probes: &[VisibleProbe]) -> ErrorReading {
    let mut sum_px = 0.0f64;
    let mut sum_world = 0.0f64;
    let mut max_px = 0.0f64;
    for probe in probes {
        let d = grid.nearest(probe.p);
        let px = d * probe.px_per_world;
        sum_px += px * px;
        sum_world += d * d;
        if px > max_px {
            max_px = px;
        }
    }
    let count = probes.len() as f64;
    ErrorReading {
        rms_px: (sum_px / count).sqrt(),
        max_px,
        rms_world: (sum_world / count).sqrt(),
    }
}

/// Triangles outside the frustum, back-facing, and the union of the two.
///
/// `off_screen_triangles` is the union, which is what C2 names. The two halves
/// are reported separately because moving the camera back removes the first and
/// cannot remove the second.
fn cull_counts(mesh: &MeshBuffer<f64>, camera: &Camera) -> (usize, usize, usize) {
    let mut frustum_off = 0;
    let mut back_facing = 0;
    let mut either = 0;
    for face in mesh.indices.as_chunks::<3>().0 {
        let a = mesh.positions[face[0] as usize];
        let b = mesh.positions[face[1] as usize];
        let c = mesh.positions[face[2] as usize];
        let centroid = scale(add(add(a, b), c), 1.0 / 3.0);
        let geometric = cross(sub(b, a), sub(c, a));
        let outside = !camera.in_frustum(centroid);
        let away = dot(geometric, sub(centroid, camera.eye)) > 0.0;
        if outside {
            frustum_off += 1;
        }
        if away {
            back_facing += 1;
        }
        if outside || away {
            either += 1;
        }
    }
    (frustum_off, back_facing, either)
}

// ─── the pop ─────────────────────────────────────────────────────────────────

/// `M-121`'s statistic for one block and one level transition: the two-sided
/// nearest-vertex Hausdorff distance between the block meshed at `level` and at
/// `level + 1`, in units of the **coarse** level's cell size.
///
/// `M-121`'s method note is honoured — the block is extracted twice at the
/// instant of the switch rather than stored at every level — and its direction is
/// widened to two-sided, because a fine vertex with no coarse neighbour is as
/// visible a pop as the reverse. A transition in which either mesh is empty
/// reports `0`: nothing moved because nothing was there.
fn block_pop<F>(
    field: &F,
    proxy: &Proxy,
    block: usize,
    level: usize,
    shapes: &[RuntimeShape3; 4],
    mc: &mut MarchingCubes<f64>,
    coarse: &mut MeshBuffer<f64>,
    fine: &mut MeshBuffer<f64>,
) -> f64
where
    F: Sdf<Scalar = f64>,
{
    let lo = proxy.blocks[block].lo;
    for (target, at) in [(&mut *coarse, level), (&mut *fine, level + 1)] {
        target.reset();
        mc.extract(field, &shapes[at], lo, proxy.cell_size(at), target)
            .expect("P-151: marching cubes on a block grid");
    }
    if coarse.positions.is_empty() || fine.positions.is_empty() {
        return 0.0;
    }
    let one_way = |from: &[[f64; 3]], to: &[[f64; 3]]| {
        from.iter()
            .map(|&v| {
                to.iter()
                    .map(|&w| dist_sq(v, w))
                    .fold(f64::INFINITY, f64::min)
            })
            .fold(0.0f64, f64::max)
    };
    let worst = one_way(&coarse.positions, &fine.positions)
        .max(one_way(&fine.positions, &coarse.positions));
    worst.sqrt() / proxy.cell_size(level)
}

// ─── timing ──────────────────────────────────────────────────────────────────

/// Median, min and max of one stage's interleaved repeats, in milliseconds.
#[derive(Clone, Copy)]
struct Timing {
    median: f64,
    min: f64,
    max: f64,
}

impl Timing {
    /// [`REPEATS`] is odd, so the median is an observation and not the average of
    /// two.
    fn of(mut ms: Vec<f64>) -> Self {
        assert!(!ms.is_empty(), "P-151: timing over no repeats");
        ms.sort_by(f64::total_cmp);
        Self {
            median: ms[ms.len() / 2],
            min: ms[0],
            max: ms[ms.len() - 1],
        }
    }

    /// max / min — how far the governor moved under this stage (`M-280`).
    fn scatter(self) -> f64 {
        self.max / self.min
    }
}

// ─── the matched read-off ────────────────────────────────────────────────────

/// One arm's `(triangles, error)` curve, sorted and deduplicated by triangle
/// count.
fn curve_of(rungs: &[(usize, f64)]) -> Vec<(f64, f64)> {
    let mut points: Vec<(f64, f64)> = rungs
        .iter()
        .map(|&(tri, err)| (tri as f64, err.max(ERROR_FLOOR_PX)))
        .collect();
    points.sort_by(|a, b| a.0.total_cmp(&b.0).then_with(|| a.1.total_cmp(&b.1)));
    points.dedup_by(|a, b| a.0.total_cmp(&b.0).is_eq());
    points
}

/// Read an arm's error off its own curve at `t_star`, by linear interpolation in
/// `(ln T, ln E)`, clamped to the curve's endpoints.
///
/// Returns the error and whether the read-off was clamped. Clamping above the
/// curve's top charges the arm the error it had at **fewer** triangles, which is
/// the larger number — conservative against whichever arm saturated.
fn read_off(curve: &[(f64, f64)], t_star: f64) -> (f64, bool) {
    assert!(!curve.is_empty(), "P-151: read-off from an empty curve");
    if t_star <= curve[0].0 {
        return (curve[0].1, t_star < curve[0].0);
    }
    let last = curve[curve.len() - 1];
    if t_star >= last.0 {
        return (last.1, t_star > last.0);
    }
    for pair in curve.windows(2) {
        let (t0, e0) = pair[0];
        let (t1, e1) = pair[1];
        if t_star <= t1 {
            let f = (t_star.ln() - t0.ln()) / (t1.ln() - t0.ln());
            return ((e0.ln() + f * (e1.ln() - e0.ln())).exp(), false);
        }
    }
    unreachable!("P-151: t_star is inside the curve but matched no interval")
}

// ─── one measured configuration ──────────────────────────────────────────────

/// Everything one (field, camera, criterion, rung) produced.
struct Row {
    criterion: Criterion,
    rung: usize,
    triangles: usize,
    error: ErrorReading,
    frustum_off: usize,
    back_facing: usize,
    off_screen: usize,
    blocks_promoted: usize,
    pop_cells: f64,
    pop_pixels: f64,
    level_histogram: [usize; 4],
    budget_bound: bool,
    estimated_triangles: f64,
    degenerate: usize,
}

/// Everything one (field, camera) pair produced.
struct Pane {
    camera: &'static str,
    probes: usize,
    visible: usize,
    /// Fraction of the field's true-surface probes outside the frustum — the
    /// registered vacuity control's number.
    off_screen_fraction: f64,
    rows: Vec<Row>,
    /// `T*`, the matched triangle count.
    matched_triangles: f64,
    /// Per criterion: the error read off at `T*`, and whether it was clamped.
    matched: [(f64, bool); 2],
    gain: f64,
    levels_differ: bool,
    adjoint: Timing,
    extract: Timing,
    proxy_ms: f64,
    curvature_ms: f64,
}

impl Pane {
    fn matched_for(&self, criterion: Criterion) -> (f64, bool) {
        self.matched[usize::from(criterion == Criterion::GoalOriented)]
    }
}

/// Everything one field produced.
struct FieldResult {
    name: &'static str,
    panes: Vec<Pane>,
}

/// Measure one reference field from both cameras under both criteria.
fn measure_field<F>(field: &F, cameras: &[Camera; 2]) -> FieldResult
where
    F: ReferenceField + Sdf<Scalar = f64>,
{
    let (lo, hi) = field.domain();
    assert!(
        field.bound().is_exact(),
        "P-151: {} declares {:?}, and probe generation Newton-projects along the gradient, \
         which only an exact distance field licenses (fields/mod.rs:83-84)",
        F::NAME,
        field.bound()
    );

    let shapes: [RuntimeShape3; 4] = [0, 1, 2, 3].map(|level| {
        RuntimeShape3::new([LEVEL_CELLS[level] + 1; 3]).expect("P-151: block grid fits u32")
    });

    let probes = surface_probes(field, lo, hi);

    let proxy_start = Instant::now();
    let mut proxy = build_proxy(field, lo, hi);
    let proxy_ms = proxy_start.elapsed().as_secs_f64() * 1e3;

    let curvature_start = Instant::now();
    let kappas = block_curvatures(field, &proxy);
    let curvature_ms = curvature_start.elapsed().as_secs_f64() * 1e3;
    for (block, kappa) in proxy.blocks.iter_mut().zip(kappas) {
        block.kappa = kappa;
    }

    let occupied = proxy.occupied();
    let tri_at_max: f64 = occupied
        .iter()
        .map(|&k| {
            TRI_PER_CELL * proxy.blocks[k].area
                / (proxy.cell_size(MAX_LEVEL) * proxy.cell_size(MAX_LEVEL))
        })
        .sum();

    // The pop is a property of (field, block, level) and of nothing else, so it
    // is measured once and read by every row that promoted that block.
    let mut mc = MarchingCubes::<f64>::new();
    let mut scratch = MeshBuffer::<f64>::new();
    let mut second = MeshBuffer::<f64>::new();
    let mut pops = vec![[0.0f64; MAX_LEVEL]; proxy.blocks.len()];
    for &k in &occupied {
        for (level, slot) in pops[k].iter_mut().enumerate() {
            *slot = block_pop(
                field,
                &proxy,
                k,
                level,
                &shapes,
                &mut mc,
                &mut scratch,
                &mut second,
            );
        }
    }

    let mut mesh = MeshBuffer::<f64>::new();
    let mut adjoint = Vec::with_capacity(proxy.blocks.len());
    let mut panes = Vec::with_capacity(cameras.len());

    for camera in cameras {
        let outside = probes.iter().filter(|q| !camera.in_frustum(q.p)).count();
        let visible: Vec<VisibleProbe> = probes
            .iter()
            .filter(|q| camera.in_frustum(q.p) && camera.front_facing(q.p, q.n))
            .map(|q| VisibleProbe {
                p: q.p,
                px_per_world: camera.pixels_per_world(q.p),
            })
            .collect();

        adjoint_weights(field, camera, &proxy, &mut adjoint);
        let ranges = block_ranges(camera, &proxy);

        let mut rows: Vec<Row> = Vec::with_capacity(2 * BUDGET_FRACTIONS.len());
        let mut allocations: Vec<(Criterion, usize, Vec<usize>)> = Vec::new();

        for criterion in [Criterion::CameraDistance, Criterion::GoalOriented] {
            let eta = energies(criterion, &proxy, camera, &adjoint, &ranges);
            for (rung, &fraction) in BUDGET_FRACTIONS.iter().enumerate() {
                let plan = allocate(&proxy, &eta, &occupied, fraction * tri_at_max);
                extract_levels(
                    field,
                    &proxy,
                    &occupied,
                    &plan.levels,
                    &shapes,
                    &mut mc,
                    &mut scratch,
                    &mut mesh,
                );
                let grid = TriGrid::new(&mesh, lo, hi, proxy.block_size * 0.25);
                let error = screen_error(&grid, &visible);
                let (frustum_off, back_facing, off_screen) = cull_counts(&mesh, camera);

                let mut promoted = 0usize;
                let mut pop_cells = 0.0f64;
                let mut pop_pixels = 0.0f64;
                let mut histogram = [0usize; 4];
                for &k in &occupied {
                    let level = plan.levels[k];
                    histogram[level] += 1;
                    if level == 0 {
                        continue;
                    }
                    promoted += 1;
                    let visible_block = proxy.blocks[k]
                        .probes
                        .iter()
                        .any(|&q| camera.in_frustum(q));
                    for (step, &cells) in pops[k].iter().enumerate().take(level) {
                        if cells > pop_cells {
                            pop_cells = cells;
                        }
                        if visible_block {
                            let pixels = cells
                                * proxy.cell_size(step)
                                * camera.f_px
                                / ranges[k].max(NEAR_PLANE);
                            if pixels > pop_pixels {
                                pop_pixels = pixels;
                            }
                        }
                    }
                }

                rows.push(Row {
                    criterion,
                    rung,
                    triangles: mesh.triangle_count(),
                    error,
                    frustum_off,
                    back_facing,
                    off_screen,
                    blocks_promoted: promoted,
                    pop_cells,
                    pop_pixels,
                    level_histogram: histogram,
                    budget_bound: plan.budget_bound,
                    estimated_triangles: plan.estimated,
                    degenerate: grid.degenerate,
                });
                allocations.push((criterion, rung, plan.levels));
            }
        }

        let levels_differ = (0..BUDGET_FRACTIONS.len()).any(|rung| {
            let of = |criterion: Criterion| {
                allocations
                    .iter()
                    .find(|(c, r, _)| *c == criterion && *r == rung)
                    .map(|(_, _, levels)| levels)
                    .expect("P-151: every criterion allocated every rung")
            };
            of(Criterion::CameraDistance) != of(Criterion::GoalOriented)
        });

        let curve = |criterion: Criterion| {
            curve_of(
                &rows
                    .iter()
                    .filter(|row| row.criterion == criterion)
                    .map(|row| (row.triangles, row.error.rms_px))
                    .collect::<Vec<_>>(),
            )
        };
        let cam_curve = curve(Criterion::CameraDistance);
        let goal_curve = curve(Criterion::GoalOriented);
        let matched_triangles = rows
            .iter()
            .find(|row| row.criterion == Criterion::CameraDistance && row.rung == REFERENCE_RUNG)
            .map(|row| row.triangles as f64)
            .expect("P-151: the camera-distance arm ran the reference rung");
        let cam_matched = read_off(&cam_curve, matched_triangles);
        let goal_matched = read_off(&goal_curve, matched_triangles);
        let gain = 1.0 - goal_matched.0 / cam_matched.0;

        // The timed pair, interleaved after one warm-up. The extraction that is
        // timed is the goal arm at the reference rung — the mesh whose criterion
        // the adjoint paid for.
        let goal_eta = energies(Criterion::GoalOriented, &proxy, camera, &adjoint, &ranges);
        let reference = allocate(
            &proxy,
            &goal_eta,
            &occupied,
            BUDGET_FRACTIONS[REFERENCE_RUNG] * tri_at_max,
        );
        adjoint_weights(field, camera, &proxy, &mut adjoint);
        extract_levels(
            field,
            &proxy,
            &occupied,
            &reference.levels,
            &shapes,
            &mut mc,
            &mut scratch,
            &mut mesh,
        );
        let mut adjoint_ms = Vec::with_capacity(REPEATS);
        let mut extract_ms = Vec::with_capacity(REPEATS);
        for _ in 0..REPEATS {
            let start = Instant::now();
            adjoint_weights(field, camera, &proxy, &mut adjoint);
            adjoint_ms.push(start.elapsed().as_secs_f64() * 1e3);
            black_box(adjoint.len());

            let start = Instant::now();
            extract_levels(
                field,
                &proxy,
                &occupied,
                &reference.levels,
                &shapes,
                &mut mc,
                &mut scratch,
                &mut mesh,
            );
            extract_ms.push(start.elapsed().as_secs_f64() * 1e3);
            black_box(mesh.triangle_count());
        }

        panes.push(Pane {
            camera: camera.name,
            probes: probes.len(),
            visible: visible.len(),
            off_screen_fraction: outside as f64 / probes.len() as f64,
            rows,
            matched_triangles,
            matched: [cam_matched, goal_matched],
            gain,
            levels_differ,
            adjoint: Timing::of(adjoint_ms),
            extract: Timing::of(extract_ms),
            proxy_ms,
            curvature_ms,
        });
    }

    FieldResult {
        name: F::NAME,
        panes,
    }
}

// ─── the run ─────────────────────────────────────────────────────────────────

fn main() {
    if !std::env::args().any(|arg| arg == "--bench") {
        return;
    }
    let prereg = isomesh::experiment!("P-151");

    common::experiment::run(prereg, |run| {
        let wide_eye = scale(unit(sub(NEAR_EYE, TARGET)), WIDE_RANGE);
        let cameras = [
            Camera::new("near", NEAR_EYE, TARGET),
            Camera::new("wide", wide_eye, TARGET),
        ];

        println!(
            "construction: BLOCK LOD, {BLOCKS}^3 blocks, levels {LEVEL_CELLS:?} cells/axis \
             (3|5|9|17 samples), proxy {PROXY_CELLS} cells/block-axis\n  \
             functional {FUNCTIONAL}: RMS over VISIBLE true-surface probes of \
             f_px * d(q, mesh) / z(q), surface->mesh direction\n  \
             camera near eye {:?}, wide eye {:?}, fov_y {FOV_Y_DEG} deg, {PIXELS_X}x{PIXELS_Y} px, \
             f_px {:.3}\n  adjoint = Z_K = h_p^2 sum_q vis(q) (f_px/z(q))^2, nothing cached; \
             {REPEATS} interleaved repeats after one warm-up\n  \
             p-146.csv metric_share over 40 rows: min {:.6} median {:.6} max {:.6}, c3_holds \
             false on every row\n",
            NEAR_EYE,
            wide_eye,
            cameras[0].f_px,
            P146_METRIC_SHARE[0],
            P146_METRIC_SHARE[1],
            P146_METRIC_SHARE[2],
        );

        let results = vec![
            measure_field(&Sphere::<f64>::canonical(), &cameras),
            measure_field(&Torus::<f64>::canonical(), &cameras),
            measure_field(&BoxExact::<f64>::canonical(), &cameras),
            measure_field(&ThinPlate::<f64>::canonical(), &cameras),
        ];

        for result in &results {
            for pane in &result.panes {
                println!(
                    "{:>11} {:>4}  probes {:>5} visible {:>5} off-frustum {:.4}  T* {:>6}  \
                     gain {:+.4}  adjoint {:.4} ms / extract {:.4} ms = {:.5}",
                    result.name,
                    pane.camera,
                    pane.probes,
                    pane.visible,
                    pane.off_screen_fraction,
                    pane.matched_triangles as u64,
                    pane.gain,
                    pane.adjoint.median,
                    pane.extract.median,
                    pane.adjoint.median / pane.extract.median,
                );
                for row in &pane.rows {
                    println!(
                        "              {:>15} rung {} tri {:>6} E {:>11.5} px  off-screen {:>6} \
                         ({:>6} frustum / {:>6} back)  pop {:.4} cells / {:.4} px  levels {:?}{}",
                        row.criterion.name(),
                        row.rung,
                        row.triangles,
                        row.error.rms_px,
                        row.off_screen,
                        row.frustum_off,
                        row.back_facing,
                        row.pop_cells,
                        row.pop_pixels,
                        row.level_histogram,
                        if row.budget_bound { "" } else { " SATURATED" },
                    );
                }
            }
        }
        println!();

        // ── vacuity controls, all before the first record ────────────────────
        //
        // M-44: a zero that could not have been non-zero is not a measurement.
        for result in &results {
            for pane in &result.panes {
                assert!(
                    pane.probes > 0,
                    "VOID: {} produced no true-surface probes at all, so the functional has no \
                     population and every error in the file is an RMS over nothing",
                    result.name
                );
                assert!(
                    pane.visible > 0,
                    "VOID: {} from the {} camera has {} probes and none of them visible, so \
                     screen_space_error is an RMS over an empty set",
                    result.name,
                    pane.camera,
                    pane.probes
                );
            }

            let near = &result.panes[0];
            assert!(
                near.off_screen_fraction > OFF_SCREEN_FLOOR,
                "VOID: the registered control. From the near camera at {NEAR_EYE:?} only \
                 {:.4} of {}'s true surface falls outside the frustum, at or under \
                 {OFF_SCREEN_FLOOR}. C2 asks whether the saving comes from not refining what is \
                 off-screen, and with the scene on-screen it cannot fire",
                near.off_screen_fraction,
                result.name
            );

            let wide = &result.panes[1];
            assert!(
                wide.off_screen_fraction < WIDE_ON_SCREEN_CEIL,
                "VOID: the C2 control is not a control. {:.4} of {}'s true surface is outside \
                 the wide camera's frustum, at or above {WIDE_ON_SCREEN_CEIL}, so 'a saving that \
                 persists with everything on-screen' has no instrument to be measured against",
                wide.off_screen_fraction,
                result.name
            );
        }

        assert!(
            results
                .iter()
                .flat_map(|r| r.panes.iter())
                .any(|pane| pane.levels_differ),
            "VOID: the two criteria assigned identical levels in every configuration, so both \
             arms extracted the same mesh and every gain in the file is a measurement of the \
             extractor's run-to-run determinism"
        );

        assert!(
            results
                .iter()
                .flat_map(|r| r.panes.iter())
                .flat_map(|pane| pane.rows.iter())
                .any(|row| row.blocks_promoted > 0),
            "VOID: no configuration promoted a single block above level 0, so \
             pop_magnitude_cells is a zero that could not have been non-zero and M-121's \
             statistic has nothing to measure"
        );

        // ── the verdicts ─────────────────────────────────────────────────────
        let c3_holds = results
            .iter()
            .flat_map(|r| r.panes.iter())
            .all(|pane| pane.adjoint.median / pane.extract.median < C3_BAR);

        let mut c1_by_field: Vec<bool> = Vec::with_capacity(results.len());
        let mut c2_by_field: Vec<bool> = Vec::with_capacity(results.len());
        for result in &results {
            let gain_near = result.panes[0].gain;
            let gain_wide = result.panes[1].gain;
            c1_by_field.push(gain_near >= C1_BAR);
            c2_by_field.push(gain_near >= C1_BAR && gain_wide < C2_PERSIST * gain_near);
        }
        let c1_fields = c1_by_field.iter().filter(|held| **held).count();
        let c2_fields = c2_by_field.iter().filter(|held| **held).count();

        println!(
            "C1: {c1_fields} of {} fields beat camera-distance LOD by >= {C1_BAR} at matched \
             triangles, near camera",
            results.len()
        );
        println!(
            "C2: {c2_fields} of {} fields kept a >= {C1_BAR} near-camera gain AND lost more than \
             {:.0}% of it with everything on-screen",
            results.len(),
            (1.0 - C2_PERSIST) * 100.0
        );
        println!("C3: every (field, camera) adjoint under {C3_BAR} of extraction -> {c3_holds}\n");

        // ── the rows ─────────────────────────────────────────────────────────
        for (index, result) in results.iter().enumerate() {
            let c1 = c1_by_field[index];
            let c2 = c2_by_field[index];
            let gain_near = result.panes[0].gain;
            let gain_wide = result.panes[1].gain;

            for pane in &result.panes {
                let share = pane.adjoint.median / pane.extract.median;
                let share_lo = pane.adjoint.min / pane.extract.max;
                let share_hi = pane.adjoint.max / pane.extract.min;
                let criterion_ms = pane.adjoint.median + pane.proxy_ms + pane.curvature_ms;

                for row in &pane.rows {
                    let (matched_error, clamped) = pane.matched_for(row.criterion);
                    run.record(&[
                        ("criterion", row.criterion.name().to_string()),
                        ("functional", FUNCTIONAL.to_string()),
                        ("triangles", row.triangles.to_string()),
                        ("screen_space_error", format!("{:.6e}", row.error.rms_px)),
                        (
                            "screen_error_at_matched_triangles",
                            format!("{matched_error:.6e}"),
                        ),
                        ("off_screen_triangles", row.off_screen.to_string()),
                        ("adjoint_ms", format!("{:.6}", pane.adjoint.median)),
                        ("adjoint_share", format!("{share:.6}")),
                        ("pop_magnitude_cells", format!("{:.6}", row.pop_cells)),
                        ("c1_holds", c1.to_string()),
                        ("c2_holds", c2.to_string()),
                        ("c3_holds", c3_holds.to_string()),
                        // ── extras (M-273) ──
                        ("adjoint_ms_max", format!("{:.6}", pane.adjoint.max)),
                        ("adjoint_ms_min", format!("{:.6}", pane.adjoint.min)),
                        ("adjoint_scatter", format!("{:.6}", pane.adjoint.scatter())),
                        ("adjoint_share_hi", format!("{share_hi:.6}")),
                        ("adjoint_share_lo", format!("{share_lo:.6}")),
                        (
                            "back_facing_fraction",
                            format!(
                                "{:.6}",
                                row.back_facing as f64 / (row.triangles.max(1)) as f64
                            ),
                        ),
                        ("back_facing_triangles", row.back_facing.to_string()),
                        ("blocks_occupied", {
                            let occupied: usize = row.level_histogram.iter().sum();
                            occupied.to_string()
                        }),
                        ("blocks_promoted", row.blocks_promoted.to_string()),
                        ("budget_bound", row.budget_bound.to_string()),
                        (
                            "budget_fraction",
                            format!("{:.4}", BUDGET_FRACTIONS[row.rung]),
                        ),
                        ("budget_rung", row.rung.to_string()),
                        ("c1_bar", format!("{C1_BAR:.2}")),
                        ("c1_fields_held", c1_fields.to_string()),
                        ("c2_fields_held", c2_fields.to_string()),
                        ("c3_row_holds", (share < C3_BAR).to_string()),
                        (
                            "c3_row_decisive",
                            ((share_lo < C3_BAR) == (share_hi < C3_BAR)).to_string(),
                        ),
                        ("camera", pane.camera.to_string()),
                        (
                            "camera_eye",
                            {
                                let eye = if pane.camera == "near" {
                                    NEAR_EYE
                                } else {
                                    wide_eye
                                };
                                format!("{:.4}x{:.4}x{:.4}", eye[0], eye[1], eye[2])
                            },
                        ),
                        ("criterion_ms", format!("{criterion_ms:.6}")),
                        (
                            "criterion_share",
                            format!("{:.6}", criterion_ms / pane.extract.median),
                        ),
                        ("curvature_ms", format!("{:.6}", pane.curvature_ms)),
                        ("degenerate_triangles", row.degenerate.to_string()),
                        (
                            "estimated_triangles",
                            format!("{:.1}", row.estimated_triangles),
                        ),
                        ("extract_ms", format!("{:.6}", pane.extract.median)),
                        ("extract_ms_max", format!("{:.6}", pane.extract.max)),
                        ("extract_ms_min", format!("{:.6}", pane.extract.min)),
                        ("field", result.name.to_string()),
                        ("fov_y_deg", format!("{FOV_Y_DEG:.1}")),
                        (
                            "frustum_off_fraction",
                            format!(
                                "{:.6}",
                                row.frustum_off as f64 / (row.triangles.max(1)) as f64
                            ),
                        ),
                        ("frustum_off_triangles", row.frustum_off.to_string()),
                        ("gain_near", format!("{gain_near:.6}")),
                        ("gain_wide", format!("{gain_wide:.6}")),
                        (
                            "level_histogram",
                            format!(
                                "{}|{}|{}|{}",
                                row.level_histogram[0],
                                row.level_histogram[1],
                                row.level_histogram[2],
                                row.level_histogram[3]
                            ),
                        ),
                        ("levels_differ", pane.levels_differ.to_string()),
                        ("matched_clamped", clamped.to_string()),
                        (
                            "matched_triangles",
                            format!("{:.1}", pane.matched_triangles),
                        ),
                        (
                            "off_screen_fraction",
                            format!("{:.6}", pane.off_screen_fraction),
                        ),
                        (
                            "p146_metric_share_max",
                            format!("{:.6}", P146_METRIC_SHARE[2]),
                        ),
                        (
                            "p146_metric_share_median",
                            format!("{:.6}", P146_METRIC_SHARE[1]),
                        ),
                        (
                            "p146_metric_share_min",
                            format!("{:.6}", P146_METRIC_SHARE[0]),
                        ),
                        ("pixels", format!("{PIXELS_X:.0}x{PIXELS_Y:.0}")),
                        ("pop_magnitude_pixels", format!("{:.6}", row.pop_pixels)),
                        ("probes", pane.probes.to_string()),
                        ("proxy_ms", format!("{:.6}", pane.proxy_ms)),
                        ("repeats", REPEATS.to_string()),
                        ("screen_error_max_px", format!("{:.6e}", row.error.max_px)),
                        (
                            "screen_error_world_rms",
                            format!("{:.6e}", row.error.rms_world),
                        ),
                        ("visible_probes", pane.visible.to_string()),
                    ]);
                }
            }
        }
    });
}
