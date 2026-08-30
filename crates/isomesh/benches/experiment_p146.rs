//! **P-146 — a metric built from the trilinear's own Hessian, and what it does
//! to the triangle count at fixed error.**
//!
//! Ticket: R-146. Pre-registered before this harness existed.
//!
//! ```bash
//! cargo bench --bench experiment_p146
//! ```
//!
//! Writes `docs/experiments/p-146.csv`.
//!
//! # What was missing
//!
//! The registration calls this **the largest paid-for gap the sweep found**:
//! five papers on metric-based anisotropic adaptation sit in the corpus at
//! 0.70–0.71 — Yano & Darmofal, Cao, Bawin et al., ParMmg, Mirebeau — and
//! before this phase the strings `Loseille`, `Alauzet`, `metric tensor` and
//! `log-Euclidean` appeared in **no** file of this repository. The optimal
//! `L^p` metric
//!
//! ```text
//!     M_Lp = D_Lp · det(|H_u|)^(−1/(2p + d)) · |H_u|,        d = 3,
//! ```
//!
//! with complexity `C(M) = ∫ √det M` standing in for the vertex count, is the
//! whole of the mechanism, and its primary here is **NASA NTRS 20200003084**,
//! which restates Loseille & Alauzet verbatim; the two SIAM originals are
//! paywalled and the corpus holds only their landing pages. A finding must cite
//! the restatement, not the papers it could not open.
//!
//! Nothing new was needed from the crate to get the Hessian. `isomesh` already
//! samples the trilinear and already differences at cell size (`M-65`), and
//! `benches/common/metric.rs:443` is that same nineteen-point stencil at that
//! same step. This row **owns** `common::metric` and is the first consumer of
//! it; R-147, R-148, R-149, R-150 and R-151 read the same eigensolver.
//!
//! # The anisotropic arm, and the one thing it is not
//!
//! **`crates/isomesh/src/` is frozen for Phase 27, so there is no anisotropic
//! extractor to call and none is written.** Every inherent `extract` in the
//! crate takes a single **scalar** `cell_size`
//! (`crates/isomesh/src/marching_cubes/mod.rs:193` and the six others), so a
//! per-cell anisotropic mesher is a source change by construction, not by
//! choice. What is built here instead is the honest bench-local reading of the
//! same mechanism:
//!
//! > **A metric-driven anisotropic sampling *grid*.** The metric prescribes a
//! > point density `√(e_aᵀ M e_a)` along each world axis `a`; averaging that
//! > over the surface band gives one weight per axis, and the anisotropic arm
//! > spends the **same total sample budget** `N³` on per-axis counts
//! > `(n_x, n_y, n_z)` proportional to those weights, against a uniform
//! > `(N, N, N)`.
//!
//! It is therefore **per-axis global anisotropy, not per-cell anisotropy**, and
//! the distinction is the single most important sentence in this file: a field
//! whose flat direction *rotates* over the surface has per-cell anisotropy this
//! arm structurally cannot spend, and it will read as a null here while a real
//! anisotropic mesher would win on it. Every verdict below is a verdict about
//! *this* construction. A per-cell result is a different row and needs a source
//! change to ask for.
//!
//! The grid stays rectilinear, so it is driven through the shipped extractor by
//! a coordinate warp: [`Stretched`] samples `f(lo + q ⊙ s)` on a **cubic** grid
//! in `q` with `s_a = h_a / h_iso`, and the emitted vertices are mapped back to
//! world space before anything measures them. Sample counts are rounded to the
//! nearest **odd** integer, because `M-266` proved the canonical grids' odd
//! counts are load-bearing — `thin_plate` is centred on `y = 0` and loses its
//! surface entirely on an even count, which would be a parity result wearing a
//! metric's clothes.
//!
//! # Arms
//!
//! | arm | what it varies | is_control |
//! |---|---|---|
//! | **isotropic** | `(N, N, N)` samples, `h` equal on every axis | **yes** — this is uniform refinement, C1's baseline |
//! | **anisotropic** | `(n_x, n_y, n_z)` from the metric, `∏ n_a ≈ N³` | no — this is the row |
//! | **matched-error read-off** | both arms fitted `ln T` against `ln E`, evaluated at one `E*` | — |
//!
//! Both arms are extracted by the same `MarchingCubes::new()` at its shipped
//! defaults, and both are measured by **one instrument**: `validate::accuracy`
//! with the *isotropic* rung's `shape`, `origin` and `cell_size` for the seed
//! lattice in **both** calls (`validate/accuracy.rs:332-337` explicitly licenses
//! a seed lattice that is not the extraction grid). A comparison in which each
//! arm grades its own homework on its own lattice is not a comparison.
//!
//! # C1's population, and the arithmetic that narrows it
//!
//! C1 asks for "at least 25% fewer triangles on **at least three of eight**
//! reference fields". Only four of the eight can carry a Hausdorff number at
//! all. `validate::accuracy` projects by Newton along `∇f` and compares against
//! `|f|`, and `crates/isomesh/src/fields/mod.rs:83-84` states the rule on
//! `FieldBound::Exact` in as many words: *"Only this admits a Hausdorff
//! measurement against the field's own values."* The roster splits
//!
//! | field | `bound()` | in C1's population |
//! |---|---|---|
//! | `sphere`, `torus`, `box_exact`, `thin_plate` | `Exact` | **yes** |
//! | `gyroid` | `Lipschitz { l: 3.464… }` | no |
//! | `csg_difference` | `Underestimate { q: 0.5 }` | no |
//! | `fbm_terrain`, `noise_cavity` | `Unbounded` | no |
//!
//! So **C1's population is 4, not 8, and its bar of 3 is 75% of what can be
//! measured rather than the 37.5% it reads as.** That is a narrowing and not an
//! impossibility — `3 ≤ 4`, so C1 is *reachable* and is run. Should a run find
//! fewer than three measurable fields, `c1_holds` is written
//! `unreachable:population=N<3` with the arithmetic on the row, which is P-70's
//! precedent and a registered outcome rather than a dropped clause. The four
//! skipped fields still emit rows: every metric, aspect-ratio, triangle and
//! cost column is measured on all eight, and only the four Hausdorff-derived
//! columns read `unmeasurable:bound=…`.
//!
//! # C2, predicted per field before the run — as the registration demands
//!
//! C2 says the win lands on the fields with a **flat direction**. For a
//! *per-axis global* grid the operative property is sharper than "has a flat
//! direction": it is **"has a flat direction that is the same world axis at
//! every point of the band"**, because a flat direction that rotates averages
//! back to isotropy in the per-axis weights. The mechanism is measured directly
//! by `flat_axis_aligned_fraction` and `flat_axis` — the fraction of band
//! points whose Hessian has a floored eigenvalue **and** at least one
//! un-floored one, with the floored eigenvector within 5° of a world axis. The
//! second condition matters: on a box face *all three* eigenvalues are floored,
//! the point is flat in every direction, and it prescribes nothing.
//!
//! | field | predicted `flat_axis` | why | predicted C1 win |
//! |---|---|---|---|
//! | `sphere` | `none` | `H(‖x‖−r)` has its zero eigenvalue along the **radial** direction, which rotates over the sphere; the two curvatures are equal | no |
//! | `torus` | `none` | axis-symmetric about `y`, but the tube's flat direction is azimuthal and rotates | no |
//! | `box_exact` | `mixed` | every face is exactly flat, but the three face families are permutation-symmetric on a cube | no |
//! | `csg_difference` | `mixed` | the subtracted sphere sits on the `(0.6, 0.6, 0.6)` diagonal (`fields/mod.rs:917-923`), so the field is permutation-symmetric too | no |
//! | `thin_plate` | `mixed` | a slab is a box; its two large faces have `H ≡ 0` and carry no direction, and the four rim edges spread over all three axes | no |
//! | `gyroid` | `none` | cubic symmetry, and not a distance field, so nothing is exactly flat | no |
//! | **`fbm_terrain`** | **`y`** | `sample` is `p[1] − (base + amp·n(x, 0, z))` (`fields/mod.rs:1352-1362`) — **exactly linear in `y`**, so `∂²f/∂y² = ∂²f/∂x∂y = ∂²f/∂y∂z ≡ 0` and `ŷ` is an exact null eigenvector at every point | **yes** |
//! | `noise_cavity` | `none` | isotropic 3-D value noise capped by a sphere | no |
//!
//! **The prediction is therefore that exactly one of eight fields wins, and
//! that it is the one whose `bound()` is `Unbounded`.** If that holds, C1 is
//! falsified — 0 of 4 against a bar of 3 — *and* C2 holds, and the finding is
//! that this crate's reference-field roster is adversarial to global
//! anisotropy: every field that could exploit it is one the Hausdorff
//! instrument cannot grade.
//!
//! `c2_holds` is per field — `predicted_win == observed_win` — and reads
//! `unmeasurable` where C1 has no verdict; `c2_mechanism_holds` compares the
//! predicted `flat_axis` against the measured one and is available on all
//! eight.
//!
//! **`flat_axis` must be read together with `flat_axis_aligned_fraction`, and
//! that pairing is not decoration.** The label is a two-thirds majority *of the
//! exploitable population*, and a two-thirds majority of a population of one is
//! still a two-thirds majority: a field with six exploitable band points out of
//! eighteen thousand can be labelled `y` while having no flat structure
//! whatever. The fraction is the column that says so, and it sits on the same
//! row. The definition is deliberately left as it stands rather than given a
//! minimum population *after* the predictions above were written, because a
//! measurement rule retuned once its answer is visible is no longer a
//! measurement — `crates/isomesh/src/experiment.rs:26-31`. A reader wanting the
//! mechanism should filter on the fraction, not on the label.
//!
//! # C3, and the arithmetic that says which way it will fall
//!
//! C3: "computing the metric costs under 15% of extraction". What is timed is
//! exactly `hessian` + `metric_lp` over the band points, five repeats after one
//! warm-up, median as the headline with min and max beside it. The band
//! *selection* is deliberately **not** timed: the extractor already visits every
//! grid sample and the band test is one comparison on a value it already has, so
//! charging the metric for that scan would be charging it twice.
//!
//! The scaling decides the verdict before the clock does. The band is a shell,
//! so `band ∝ N²`; extraction is a volume, so `extract ∝ N³`; therefore
//! `metric_share ∝ 1/N` and `metric_share × N` is the scale-free constant —
//! recorded as `share_times_resolution`, from which the crossing resolution is
//! `share_times_resolution / 0.15`. One band point costs nineteen field samples
//! plus one Jacobi eigendecomposition; one extraction cell costs about one field
//! sample. At the ladder's resolutions the band is a few percent of the volume
//! and the per-point cost ratio is one to two orders larger than that, so **C3
//! is predicted FALSIFIED at every rung**, with the crossing recorded rather
//! than guessed. `c3_holds` is the **global** verdict (every row under 15%) and
//! is the same on every row; `c3_row_holds` is that row's own, and
//! `c3_row_decisive` says whether the min/max spread straddles the bar — this
//! machine's `amd-pstate-epp` governor swings the same binary 1.45× between runs
//! (`M-280`), so a share is reported with its scatter rather than averaged into
//! a pass.
//!
//! # SHARE, recomputed before the numbers
//!
//! The registration's SHARE is *"C1 moves the triangle budget, whose share of
//! frame cost `M-135` puts at 29% for the contour and 45% for the collider
//! check."* Discharged as an Amdahl ceiling: a 25% triangle reduction — C1's own
//! bar, and the most C1 can claim while merely holding — bounds the saving at
//! `0.25 × 0.29 = 7.25%` of the contour stage and `0.25 × 0.45 = 11.25%` of the
//! collider check, i.e. **at most 18.5% of a frame**, and only if both stages are
//! exactly linear in triangle count. Against the prediction above the realised
//! share is **zero**: no measurable field is predicted to win, so there is no
//! triangle budget to move and this row proposes no landing. The number that
//! would justify one is `ratio` on a field with `bound() == Exact`, and Phase 28
//! is where such a landing would be registered.
//!
//! # Vacuity controls
//!
//! Every one runs before the first `run.record`, and every panic starts
//! `"VOID: "`. `M-44`: a zero that could not have been non-zero is not a
//! measurement.
//!
//! - **The registered control** — `aspect_ratio_max > 3` on at least one field,
//!   or the "anisotropic" arm is isotropic and C1 is measuring noise. Column:
//!   `aspect_ratio_max`.
//! - **The module author's amendment to it, and it is load-bearing.**
//!   `benches/common/metric.rs:67-74` warns that at a flat direction
//!   `aspect_ratio` is `|λ|max / H_FLOOR ≈ 1e9`–`1e11`, *"a number of order 1e11
//!   that is the floor talking, not a measured anisotropy"*. So the maximum is
//!   asserted **again** over the cells that are **not** at the floor
//!   (`aspect_ratio_max_off_floor > 3` somewhere), and every row carries
//!   `at_floor_cells` beside `aspect_ratio_max` so a reader can tell which
//!   number they are looking at. The module's own measured at-floor counts are
//!   `box_exact` 1686/1790, `fbm_terrain` 1156/1156 and `gyroid` 1/2945 with a
//!   genuine `5.11e3`; the expected witness for this control is therefore
//!   `gyroid`, which is Lipschitz rather than a distance field and so has no
//!   exactly-flat direction to floor.
//! - **The anisotropic arm must be anisotropic.** `axis_ratio = max n_a / min n_a`
//!   must exceed 1.5 on at least one row, or the two arms are the same grid under
//!   two names. Expected witness: `fbm_terrain`.
//! - **The two arms must be on one budget.** `budget_ratio = ∏n_a / N³` in
//!   `[0.25, 4]` on every row, or "same total sample budget" is not what was run.
//! - **The Hausdorff instrument must respond to resolution.** A field enters C1's
//!   population only if both arms' errors span at least 1.2× across the ladder;
//!   otherwise `ratio` would be read off a fit through a horizontal line. At
//!   least one field must survive that test.
//! - **The band must be non-empty** on every row, or `complexity_target` and both
//!   aspect ratios are statistics of nothing.
//!
//! # A residual, recorded so the four skipped fields are not silent
//!
//! `mesh_residual_max_iso` / `mesh_residual_max_aniso` are `max |f(v)| / h` over
//! mesh vertices. That is **not** a Hausdorff distance and it grades nothing —
//! for a non-`Exact` field `|f|` is not a distance, which is the whole reason
//! those fields are out of C1's population. It is recorded because a field-value
//! residual is still a real number about a real mesh, and a row that reports
//! only `unmeasurable` four times over throws away the observation that the
//! anisotropic arm placed its vertices somewhere.
//!
//! # Determinism
//!
//! One thread, no PRNG, no map iteration, `f64` throughout. Sorting is
//! [`f64::total_cmp`]. The band is swept `z`, `y`, `x` with `x` innermost, the
//! crate's order. The one seeded object anywhere near this row is
//! `FbmTerrain`'s committed `0x5EED_1234`, which is the field's and not the
//! harness's. Wall-clock columns are the only machine-dependent numbers, and
//! only C3 reads them.

#![allow(
    clippy::cast_possible_truncation,
    clippy::cast_precision_loss,
    clippy::cast_sign_loss,
    clippy::too_many_lines
)]

mod common;

use std::hint::black_box;
use std::time::Instant;

use isomesh::fields::{FieldBound, ReferenceField};
use isomesh::marching_cubes::MarchingCubes;
use isomesh::validate::{AccuracyConfig, accuracy};
use isomesh::{MeshBuffer, RuntimeShape3, Sdf};

use crate::common::metric::{H_FLOOR, Sym3, complexity, hessian, metric_lp};

/// Samples per axis for the isotropic arm — the ladder both arms are budgeted
/// against, and the ladder the matched-error read-off is fitted through.
///
/// Five rungs spanning `17 → 65`, a factor of `4` in `h` and so roughly `16×`
/// in a second-order error. Odd throughout (`M-266`).
const RUNGS: [u32; 5] = [17, 25, 33, 49, 65];

/// The norm `M_Lp` is optimised for.
///
/// One value, and `L²`. `R-150` owns the norm sweep — *"R-146's metric under
/// ≥ 3 norms"* — and running it here too would be two answers to one question.
const P_NORM: f64 = 2.0;

/// The `metric` column's value: which metric field was built.
const METRIC_NAME: &str = "M_Lp_hessian";

/// A grid sample joins the surface band when `|f| <= BAND_CELLS · h`.
///
/// One cell, not three: the metric is wanted where the mesh is, and a wider
/// shell would charge C3 for points no triangle is near.
const BAND_CELLS: f64 = 1.0;

/// Timed repeats per stage, after one warm-up. Five is the authoring contract's
/// floor for a clause with a ratio threshold.
const REPEATS: usize = 5;

/// Fewest samples any axis of the anisotropic grid may carry.
///
/// Odd, and five rather than two: `Error::GridTooSmall` starts at two, and a
/// two-sample axis has one cell and can only place a crossing against the wall.
const MIN_SAMPLES: u32 = 5;

/// C1's bar: "at least 25% fewer triangles".
const WIN_RATIO: f64 = 0.75;

/// C1's bar on the field count.
const C1_MIN_WINNERS: usize = 3;

/// C3's bar.
const C3_MAX_SHARE: f64 = 0.15;

/// The registered vacuity control's bar on `aspect_ratio_max`.
const ASPECT_FLOOR: f64 = 3.0;

/// The anisotropic arm must be anisotropic somewhere, by at least this factor.
const AXIS_RATIO_FLOOR: f64 = 1.5;

/// `cos(5°)`. A floored eigenvector counts as axis-aligned above this.
const AXIS_ALIGNED_COS: f64 = 0.996_194_698_091_745_5;

/// The ladder must move the error instrument by at least this factor, or a fit
/// through it has no slope to read.
const LADDER_SPAN_FLOOR: f64 = 1.2;

/// `∏ n_a / N³` must land inside this, or the arms are not budget-matched.
const BUDGET_BAND: [f64; 2] = [0.25, 4.0];

/// `for_each_reference_field!` yields eight (`fields/mod.rs:211-255`).
const FIELDS: usize = 8;

/// C2's per-field prediction, written into this file before the harness had
/// ever run, which is the whole point of the clause.
///
/// `(field, predicted C1 win, predicted flat axis)`. The reasoning for every row
/// is the table in the header; the short version is that a *per-axis global*
/// grid can only spend a flat direction that is the **same world axis
/// everywhere**, and exactly one reference field has one.
const PREDICTED: [(&str, bool, &str); FIELDS] = [
    ("sphere", false, "none"),
    ("torus", false, "none"),
    ("box_exact", false, "mixed"),
    ("csg_difference", false, "mixed"),
    ("thin_plate", false, "mixed"),
    ("gyroid", false, "none"),
    ("fbm_terrain", true, "y"),
    ("noise_cavity", false, "none"),
];

/// The registered prediction for one field.
///
/// Panics for a name that was not predicted: an unpredicted field is a C2 clause
/// with nothing to compare against, and defaulting it would answer C2 by
/// choosing the answer.
fn predicted(name: &str) -> (bool, &'static str) {
    for (field, win, axis) in PREDICTED {
        if field == name {
            return (win, axis);
        }
    }
    panic!("P-146: no C2 prediction was registered for field `{name}`");
}

/// The name of a `FieldBound` variant, without its parameters — the CSV writer
/// refuses a `,` inside a value and `Lipschitz { l: 3.46 }` has one.
fn bound_name(bound: FieldBound) -> &'static str {
    match bound {
        FieldBound::Exact => "Exact",
        FieldBound::Lipschitz { .. } => "Lipschitz",
        FieldBound::Underestimate { .. } => "Underestimate",
        FieldBound::Unbounded => "Unbounded",
    }
}

/// Why a field's Hausdorff cannot be measured, from its declared bound.
fn unmeasurable_reason(bound: FieldBound) -> &'static str {
    match bound {
        FieldBound::Exact => "measurable",
        FieldBound::Lipschitz { .. } => "unmeasurable:bound=Lipschitz",
        FieldBound::Underestimate { .. } => "unmeasurable:bound=Underestimate",
        FieldBound::Unbounded => "unmeasurable:bound=Unbounded",
    }
}

// ─── timing ──────────────────────────────────────────────────────────────────

/// Median, min and max of one stage's repeats, in milliseconds.
///
/// The median rather than the mean because `M-280` measured a 1.45× swing on
/// this host's governor and one slow repeat should not move the headline; the
/// min and max are carried so a reader sees the swing rather than taking the
/// median on trust.
#[derive(Clone, Copy)]
struct Timing {
    median: f64,
    min: f64,
    max: f64,
}

/// [`REPEATS`] is odd, so the median is an observation rather than an average of
/// two.
fn timing(mut ms: Vec<f64>) -> Timing {
    ms.sort_by(f64::total_cmp);
    Timing {
        median: ms[ms.len() / 2],
        min: ms[0],
        max: ms[ms.len() - 1],
    }
}

// ─── the metric field ────────────────────────────────────────────────────────

/// The grid samples within one cell of the surface, swept `z`, `y`, `x` with
/// `x` innermost.
fn band_points<F>(field: &F, origin: [f64; 3], h: f64, samples: u32) -> Vec<[f64; 3]>
where
    F: Sdf<Scalar = f64>,
{
    let band = BAND_CELLS * h;
    let mut out = Vec::new();
    for k in 0..samples {
        for j in 0..samples {
            for i in 0..samples {
                let p = [
                    origin[0] + f64::from(i) * h,
                    origin[1] + f64::from(j) * h,
                    origin[2] + f64::from(k) * h,
                ];
                if field.sample(p).abs() <= band {
                    out.push(p);
                }
            }
        }
    }
    out
}

/// The timed stage: `hessian` then `metric_lp`, once per band point.
///
/// This and nothing else is what C3 charges the metric for. `out` is reused
/// across repeats so the allocator is not in the measurement.
fn build_metrics<F>(field: &F, points: &[[f64; 3]], h: f64, out: &mut Vec<Sym3>)
where
    F: Sdf<Scalar = f64>,
{
    out.clear();
    for &p in points {
        let hess = hessian(field, p, h);
        out.push(metric_lp(&hess, P_NORM));
    }
    let _ = black_box(&*out);
}

/// What the metric field says about one `(field, resolution)` — all of it
/// untimed, because none of it is part of building the metric and none of it may
/// be charged to C3.
struct Census {
    points: usize,
    at_floor: usize,
    off_floor: usize,
    /// Band points with a floored eigenvalue, at least one un-floored one, and
    /// the floored eigenvector within 5° of a world axis.
    exploitable_flat: usize,
    axis_hits: [usize; 3],
    aspect_max: f64,
    aspect_mean: f64,
    aspect_max_off_floor: f64,
    aspect_mean_off_floor: f64,
    /// Mean `√(e_aᵀ M e_a)` per axis: the metric's own point density along each
    /// world axis, and the only thing the anisotropic split is derived from.
    weights: [f64; 3],
    complexity: f64,
}

/// Census the metric field.
///
/// Recomputes the Hessian because the at-floor test is a question about `H`, not
/// about `M`: `metric_lp` scales every eigenvalue by one common factor, so the
/// floor is invisible in the metric alone.
fn census_of<F>(field: &F, points: &[[f64; 3]], metrics: &[Sym3], h: f64) -> Census
where
    F: Sdf<Scalar = f64>,
{
    let mut census = Census {
        points: points.len(),
        at_floor: 0,
        off_floor: 0,
        exploitable_flat: 0,
        axis_hits: [0; 3],
        aspect_max: 0.0,
        aspect_mean: 0.0,
        aspect_max_off_floor: 0.0,
        aspect_mean_off_floor: 0.0,
        weights: [0.0; 3],
        complexity: 0.0,
    };
    let mut aspect_sum = 0.0f64;
    let mut aspect_sum_off_floor = 0.0f64;

    for (&p, m) in points.iter().zip(metrics) {
        let (values, vectors) = hessian(field, p, h).eigen();

        // The flat direction is the smallest |eigenvalue|; `eigen` sorts by
        // value and not by magnitude, so a saddle's smallest magnitude is not
        // `values[0]`.
        let mut flat = 0usize;
        let mut stiff = 0usize;
        for (index, value) in values.iter().enumerate() {
            if value.abs() < values[flat].abs() {
                flat = index;
            }
            if value.abs() > H_FLOOR {
                stiff += 1;
            }
        }
        let floored = values[flat].abs() <= H_FLOOR;

        let aspect = m.aspect_ratio();
        aspect_sum += aspect;
        if aspect > census.aspect_max {
            census.aspect_max = aspect;
        }
        if floored {
            census.at_floor += 1;
        } else {
            census.off_floor += 1;
            aspect_sum_off_floor += aspect;
            if aspect > census.aspect_max_off_floor {
                census.aspect_max_off_floor = aspect;
            }
        }

        // `stiff > 0` is what makes a flat direction *exploitable*: a point flat
        // in every direction — a box face — prescribes nothing, and counting it
        // would make every box look like a heightfield.
        if floored && stiff > 0 {
            let vector = [vectors[0][flat], vectors[1][flat], vectors[2][flat]];
            let mut axis = 0usize;
            for (index, component) in vector.iter().enumerate() {
                if component.abs() > vector[axis].abs() {
                    axis = index;
                }
            }
            if vector[axis].abs() >= AXIS_ALIGNED_COS {
                census.exploitable_flat += 1;
                census.axis_hits[axis] += 1;
            }
        }

        for (axis, weight) in census.weights.iter_mut().enumerate() {
            *weight += m.get(axis, axis).sqrt();
        }
    }

    let n = points.len() as f64;
    census.aspect_mean = aspect_sum / n;
    if census.off_floor > 0 {
        census.aspect_mean_off_floor = aspect_sum_off_floor / census.off_floor as f64;
    }
    for weight in &mut census.weights {
        *weight /= n;
    }
    census.complexity = complexity(metrics, h * h * h);
    census
}

/// The dominant axis of the exploitable flat directions, or `mixed` / `none`.
///
/// "Dominant" is a two-thirds majority: a field whose flat directions spread
/// over two or three axes has no axis for the grid to spend a budget on, and
/// calling the largest bucket the answer would launder that into a prediction
/// match.
fn flat_axis(census: &Census) -> &'static str {
    if census.exploitable_flat == 0 {
        return "none";
    }
    let names = ["x", "y", "z"];
    for (axis, &hits) in census.axis_hits.iter().enumerate() {
        if hits * 3 >= census.exploitable_flat * 2 {
            return names[axis];
        }
    }
    "mixed"
}

// ─── the anisotropic grid ────────────────────────────────────────────────────

/// Nearest odd integer, at least one. Ties go up, deterministically.
fn round_odd(x: f64) -> u32 {
    let half = ((x - 1.0) * 0.5).round();
    (2.0f64.mul_add(half, 1.0)).max(1.0) as u32
}

/// Per-axis sample counts from the metric's per-axis point densities, at the
/// isotropic arm's total budget.
///
/// `n_a ∝ weights[a]` with `∏ n_a = N³`, so the constant is
/// `k = (N³)^{1/3} / geomean(weights)` and the two arms differ in *shape* only.
///
/// **There is exactly one clamp and it is a lower one.** An axis that would
/// fall below [`MIN_SAMPLES`] is pinned there and the remaining budget is
/// *re-solved* over the axes still free, which is what keeps `∏ n_a` on the
/// budget instead of merely near it. At most three rounds, because every round
/// that changes anything pins at least one more axis.
///
/// No upper clamp, and the first version of this had one — a ceiling of
/// `budget / MIN_SAMPLES²`, the most a single axis can carry with the other two
/// at the floor. It was wrong for a reason worth writing down: it can bind on
/// **two** axes in the same round, before the lower pins have been resolved,
/// and two axes at a ceiling each sized for the other two being at the floor
/// double-counts the budget. Measured on `fbm_terrain`, whose `ŷ` weight is
/// `~1e-4` of the other two: `29x5x29` against a budget of `9³`, **5.77× over**,
/// rising to `39.5×` at `17³` — caught by the budget vacuity control rather
/// than by reading the code. Re-solving alone is enough: at least one free axis
/// always has weight at or above the free set's geometric mean, so it receives
/// at least `target^{1/count}`, and the last free axis receives the whole
/// remaining budget — `17³/5² = 196` samples at the ladder's coarsest rung, two
/// orders above [`MIN_SAMPLES`]. So no axis can end below the floor and no
/// ceiling is needed to bound one from above.
///
/// Returns the counts and how many axes were pinned at the floor.
fn anisotropic_grid(weights: [f64; 3], samples: u32) -> ([u32; 3], usize) {
    let budget = f64::from(samples).powi(3);
    let mut pinned = [false; 3];
    let mut n = [samples; 3];

    for _round in 0..3 {
        let free: Vec<usize> = (0..3).filter(|&axis| !pinned[axis]).collect();
        if free.is_empty() {
            break;
        }
        let mut held = 1.0f64;
        for (axis, &fixed) in pinned.iter().enumerate() {
            if fixed {
                held *= f64::from(n[axis]);
            }
        }
        let count = free.len() as f64;
        let target = budget / held;
        let mut logsum = 0.0f64;
        for &axis in &free {
            logsum += weights[axis].ln();
        }
        let geometric_mean = (logsum / count).exp();
        let scale = target.powf(1.0 / count) / geometric_mean;

        let mut newly_pinned = false;
        for &axis in &free {
            let raw = round_odd(scale * weights[axis]);
            if raw < MIN_SAMPLES {
                n[axis] = MIN_SAMPLES;
                pinned[axis] = true;
                newly_pinned = true;
            } else {
                n[axis] = raw;
            }
        }
        if !newly_pinned {
            break;
        }
    }

    (n, pinned.iter().filter(|fixed| **fixed).count())
}

/// The field seen through a per-axis coordinate stretch.
///
/// `sample(q) = f(lo + q ⊙ s)`. Extracting this on a **cubic** grid of
/// `cell_size = h` in `q` is exactly extracting `f` on a rectilinear grid whose
/// physical spacings are `h · s`, which is the only way to reach a rectilinear
/// anisotropic grid through an `extract` that takes a scalar `cell_size`
/// (`marching_cubes/mod.rs:193`).
///
/// `gradient` is deliberately left as `Sdf`'s central-difference default: it is
/// read only for the emitted normals, and nothing downstream of this harness
/// reads a normal. Positions are mapped back to world space by the caller before
/// any measurement touches them.
struct Stretched<'a, F> {
    field: &'a F,
    lo: [f64; 3],
    s: [f64; 3],
}

impl<F> Sdf for Stretched<'_, F>
where
    F: Sdf<Scalar = f64>,
{
    type Scalar = f64;

    fn sample(&self, q: [f64; 3]) -> f64 {
        self.field.sample([
            q[0].mul_add(self.s[0], self.lo[0]),
            q[1].mul_add(self.s[1], self.lo[1]),
            q[2].mul_add(self.s[2], self.lo[2]),
        ])
    }
}

/// `max |f(v)| / h` over the mesh's vertices. A residual, **not** a distance —
/// see the header.
fn residual_max<F>(field: &F, mesh: &MeshBuffer<f64>, h: f64) -> f64
where
    F: Sdf<Scalar = f64>,
{
    let mut worst = 0.0f64;
    for &p in &mesh.positions {
        let value = field.sample(p).abs();
        if value > worst {
            worst = value;
        }
    }
    worst / h
}

// ─── the matched-error read-off ──────────────────────────────────────────────

/// Least squares of `ln(triangles)` against `ln(error)`.
struct Fit {
    intercept: f64,
    slope: f64,
    r2: f64,
}

impl Fit {
    /// The fitted triangle count at error `e`.
    fn eval(&self, e: f64) -> f64 {
        self.slope.mul_add(e.ln(), self.intercept).exp()
    }
}

/// Fit one arm's ladder.
///
/// The caller has already established that the errors span
/// [`LADDER_SPAN_FLOOR`], so the slope's denominator is non-zero.
fn fit_log(points: &[(f64, f64)]) -> Fit {
    let n = points.len() as f64;
    let mut mean_x = 0.0f64;
    let mut mean_y = 0.0f64;
    for &(error, triangles) in points {
        mean_x += error.ln();
        mean_y += triangles.ln();
    }
    mean_x /= n;
    mean_y /= n;

    let mut sxx = 0.0f64;
    let mut sxy = 0.0f64;
    for &(error, triangles) in points {
        let dx = error.ln() - mean_x;
        sxx = dx.mul_add(dx, sxx);
        sxy = dx.mul_add(triangles.ln() - mean_y, sxy);
    }
    let slope = sxy / sxx;
    let intercept = slope.mul_add(-mean_x, mean_y);

    let mut ssres = 0.0f64;
    let mut sstot = 0.0f64;
    for &(error, triangles) in points {
        let y = triangles.ln();
        let residual = y - slope.mul_add(error.ln(), intercept);
        ssres = residual.mul_add(residual, ssres);
        sstot = (y - mean_y).mul_add(y - mean_y, sstot);
    }
    // `sstot == 0` means every rung emitted the same triangle count. The
    // ladder-span control has excluded that on the error axis but not on this
    // one, and `0` is the honest reading of "the fit explains none of a variance
    // that does not exist".
    let r2 = if sstot > 0.0 {
        1.0 - ssres / sstot
    } else {
        0.0
    };

    Fit {
        intercept,
        slope,
        r2,
    }
}

/// One field's C1 answer, once both ladders exist.
struct Matched {
    error: f64,
    triangles_iso: f64,
    triangles_aniso: f64,
    ratio: f64,
    win: bool,
    interpolated: bool,
    r2_iso: f64,
    r2_aniso: f64,
}

/// Read both arms off at one matched error, or say why not.
///
/// The matched point `E*` is the **isotropic arm's finest error** — a concrete,
/// measured quantity, and the natural question: "at the error uniform refinement
/// buys at 65³, what does each arm cost?" Both arms are read from their fits
/// rather than one from a fit and one from a raw point, so neither gets a
/// smoothing the other did not.
fn match_arms(rungs: &[Rung]) -> Result<Matched, String> {
    let mut iso: Vec<(f64, f64)> = Vec::with_capacity(rungs.len());
    let mut aniso: Vec<(f64, f64)> = Vec::with_capacity(rungs.len());
    for rung in rungs {
        let pair = rung.hausdorff.map_err(String::from)?;
        if rung.triangles_iso == 0 || rung.triangles_aniso == 0 {
            return Err(String::from("unmeasurable:empty_mesh"));
        }
        if pair[0] <= 0.0 || pair[1] <= 0.0 {
            return Err(String::from("unmeasurable:zero_error"));
        }
        iso.push((pair[0], rung.triangles_iso as f64));
        aniso.push((pair[1], rung.triangles_aniso as f64));
    }

    let span = |points: &[(f64, f64)]| {
        let mut lo = f64::INFINITY;
        let mut hi = 0.0f64;
        for &(error, _) in points {
            lo = lo.min(error);
            hi = hi.max(error);
        }
        (lo, hi)
    };
    let (iso_lo, iso_hi) = span(&iso);
    let (aniso_lo, aniso_hi) = span(&aniso);
    if iso_hi < iso_lo * LADDER_SPAN_FLOOR || aniso_hi < aniso_lo * LADDER_SPAN_FLOOR {
        return Err(String::from("unmeasurable:flat_ladder"));
    }

    let fit_iso = fit_log(&iso);
    let fit_aniso = fit_log(&aniso);
    let error = iso_lo;
    let triangles_iso = fit_iso.eval(error);
    let triangles_aniso = fit_aniso.eval(error);
    let ratio = triangles_aniso / triangles_iso;

    Ok(Matched {
        error,
        triangles_iso,
        triangles_aniso,
        ratio,
        win: ratio <= WIN_RATIO,
        interpolated: error >= aniso_lo && error <= aniso_hi,
        r2_iso: fit_iso.r2,
        r2_aniso: fit_aniso.r2,
    })
}

// ─── one rung, one field ─────────────────────────────────────────────────────

/// Everything one `(field, resolution)` produced.
struct Rung {
    samples: u32,
    census: Census,
    grid: [u32; 3],
    pinned: usize,
    budget_ratio: f64,
    axis_ratio: f64,
    triangles_iso: u64,
    triangles_aniso: u64,
    /// `[isotropic, anisotropic]` symmetric Hausdorff, or the reason there is
    /// none.
    hausdorff: Result<[f64; 2], &'static str>,
    residual_iso: f64,
    residual_aniso: f64,
    metric_ms: Timing,
    extract_ms: Timing,
}

/// One field's five rungs and the verdicts they decide together.
struct FieldResult {
    name: &'static str,
    bound: &'static str,
    predicted_win: bool,
    predicted_axis: &'static str,
    rungs: Vec<Rung>,
    matched: Result<Matched, String>,
}

/// Measure one reference field across the whole ladder.
fn measure_field<F>(field: &F, name: &'static str) -> FieldResult
where
    F: ReferenceField + Sdf<Scalar = f64>,
{
    let bound = field.bound();
    let (predicted_win, predicted_axis) = predicted(name);
    let (lo, hi) = field.domain();
    let extent = hi[0] - lo[0];

    let mut mc = MarchingCubes::<f64>::new();
    let mut mesh = MeshBuffer::<f64>::new();
    let mut aniso_mesh = MeshBuffer::<f64>::new();
    let mut metrics: Vec<Sym3> = Vec::new();
    let mut rungs: Vec<Rung> = Vec::with_capacity(RUNGS.len());

    for samples in RUNGS {
        let (shape, origin, h) = common::grid::<f64, _>(field, samples);

        // ── the metric field, timed exactly as C3 charges it ────────────────
        let points = band_points(field, origin, h, samples);
        assert!(
            !points.is_empty(),
            "VOID: {name} at {samples}^3 put no grid sample within {BAND_CELLS} cell of its \
             surface, so complexity_target, aspect_ratio_max and aspect_ratio_mean would be \
             statistics of an empty population and every zero among them a zero that could not \
             have been non-zero (M-44)"
        );
        build_metrics(field, &points, h, &mut metrics);
        let mut metric_ms = Vec::with_capacity(REPEATS);
        for _ in 0..REPEATS {
            let start = Instant::now();
            build_metrics(field, &points, h, &mut metrics);
            metric_ms.push(start.elapsed().as_secs_f64() * 1e3);
        }
        let metric_ms = timing(metric_ms);
        let census = census_of(field, &points, &metrics, h);

        // ── the isotropic arm, timed ────────────────────────────────────────
        mesh.reset();
        mc.extract(field, &shape, origin, h, &mut mesh)
            .expect("isotropic extraction over the reference grid");
        let mut extract_ms = Vec::with_capacity(REPEATS);
        for _ in 0..REPEATS {
            mesh.reset();
            let start = Instant::now();
            mc.extract(field, &shape, origin, h, &mut mesh)
                .expect("isotropic extraction over the reference grid");
            extract_ms.push(start.elapsed().as_secs_f64() * 1e3);
        }
        let extract_ms = timing(extract_ms);
        let triangles_iso = mesh.triangle_count() as u64;

        // ── the anisotropic arm, at the same total sample budget ────────────
        let (grid, pinned) = anisotropic_grid(census.weights, samples);
        let stretch = [
            extent / f64::from(grid[0] - 1) / h,
            extent / f64::from(grid[1] - 1) / h,
            extent / f64::from(grid[2] - 1) / h,
        ];
        let stretched = Stretched {
            field,
            lo,
            s: stretch,
        };
        let aniso_shape = RuntimeShape3::new(grid).expect("anisotropic grid fits u32");
        aniso_mesh.reset();
        mc.extract(&stretched, &aniso_shape, [0.0; 3], h, &mut aniso_mesh)
            .expect("anisotropic extraction over the stretched grid");
        for p in &mut aniso_mesh.positions {
            p[0] = p[0].mul_add(stretch[0], lo[0]);
            p[1] = p[1].mul_add(stretch[1], lo[1]);
            p[2] = p[2].mul_add(stretch[2], lo[2]);
        }
        let triangles_aniso = aniso_mesh.triangle_count() as u64;

        // ── one instrument, both arms ───────────────────────────────────────
        let hausdorff = if bound.is_exact() {
            let cfg = AccuracyConfig::from_cell_size(h).expect("positive cell size");
            let iso_report = accuracy(&mesh.positions, &mesh.indices, field, &shape, origin, &cfg)
                .expect("accuracy over the isotropic arm");
            let aniso_report = accuracy(
                &aniso_mesh.positions,
                &aniso_mesh.indices,
                field,
                &shape,
                origin,
                &cfg,
            )
            .expect("accuracy over the anisotropic arm on the isotropic seed lattice");
            if iso_report.has_coverage() && aniso_report.has_coverage() {
                Ok([
                    iso_report.symmetric_hausdorff(),
                    aniso_report.symmetric_hausdorff(),
                ])
            } else {
                Err("unmeasurable:no_coverage")
            }
        } else {
            Err(unmeasurable_reason(bound))
        };

        let budget_ratio = (f64::from(grid[0]) * f64::from(grid[1]) * f64::from(grid[2]))
            / f64::from(samples).powi(3);
        let axis_hi = f64::from(grid.iter().copied().max().unwrap_or(samples));
        let axis_lo = f64::from(grid.iter().copied().min().unwrap_or(samples));

        rungs.push(Rung {
            samples,
            census,
            grid,
            pinned,
            budget_ratio,
            axis_ratio: axis_hi / axis_lo,
            triangles_iso,
            triangles_aniso,
            hausdorff,
            residual_iso: residual_max(field, &mesh, h),
            residual_aniso: residual_max(field, &aniso_mesh, h),
            metric_ms,
            extract_ms,
        });
    }

    let matched = match_arms(&rungs);
    FieldResult {
        name,
        bound: bound_name(bound),
        predicted_win,
        predicted_axis,
        rungs,
        matched,
    }
}

/// One rung's line on the console.
fn report(result: &FieldResult, rung: &Rung) {
    let show = |index: usize| {
        rung.hausdorff
            .map_or_else(str::to_string, |pair| format!("{:.9}", pair[index]))
    };
    println!(
        "{:>14} {:>4}^3  aniso {:>5}x{:>5}x{:>5} (ratio {:>9.3}  pinned {}  budget {:>6.3})  \
         tri {:>8} / {:>8}  haus {:>28} / {:>28}  aspect {:>10.3e} (off-floor {:>10.3e}, \
         at-floor {:>7}/{:>7})  flat {:>5} {:>6.3}  metric {:>9.3} ms  extract {:>9.3} ms  \
         share {:>9.4}",
        result.name,
        rung.samples,
        rung.grid[0],
        rung.grid[1],
        rung.grid[2],
        rung.axis_ratio,
        rung.pinned,
        rung.budget_ratio,
        rung.triangles_iso,
        rung.triangles_aniso,
        show(0),
        show(1),
        rung.census.aspect_max,
        rung.census.aspect_max_off_floor,
        rung.census.at_floor,
        rung.census.points,
        flat_axis(&rung.census),
        rung.census.exploitable_flat as f64 / rung.census.points as f64,
        rung.metric_ms.median,
        rung.extract_ms.median,
        rung.metric_ms.median / rung.extract_ms.median,
    );
}

// ─── the run ─────────────────────────────────────────────────────────────────

fn main() {
    if !std::env::args().any(|arg| arg == "--bench") {
        return;
    }
    let prereg = isomesh::experiment!("P-146");

    common::experiment::run(prereg, |run| {
        println!(
            "construction: metric-driven anisotropic sampling GRID, per-axis global, NOT \
             per-cell.\n  metric {METRIC_NAME}, p = {P_NORM}, band |f| <= {BAND_CELLS}h, \
             H_FLOOR = {H_FLOOR:e}\n  ladder {RUNGS:?} samples/axis, {REPEATS} timed repeats \
             after one warm-up\n  Hausdorff: validate::accuracy, ONE seed lattice per rung, \
             shared by both arms\n"
        );

        let mut results: Vec<FieldResult> = Vec::with_capacity(FIELDS);
        isomesh::for_each_reference_field!(f64, |name, field| {
            results.push(measure_field(&field, name));
        });
        assert_eq!(
            results.len(),
            FIELDS,
            "P-146: for_each_reference_field! must yield {FIELDS} fields"
        );

        for result in &results {
            for rung in &result.rungs {
                report(result, rung);
            }
        }
        println!();

        // ── vacuity controls, all before the first record ────────────────────
        //
        // M-44: a zero that could not have been non-zero is not a measurement.
        let aspect_max = results
            .iter()
            .flat_map(|result| result.rungs.iter())
            .map(|rung| rung.census.aspect_max)
            .fold(0.0f64, f64::max);
        assert!(
            aspect_max > ASPECT_FLOOR,
            "VOID: the largest aspect ratio anywhere in the sweep is {aspect_max:e}, which does \
             not exceed {ASPECT_FLOOR}. That is the registration's own vacuity control: the \
             'anisotropic' arm is then isotropic and C1 is measuring noise"
        );

        let aspect_max_off_floor = results
            .iter()
            .flat_map(|result| result.rungs.iter())
            .map(|rung| rung.census.aspect_max_off_floor)
            .fold(0.0f64, f64::max);
        assert!(
            aspect_max_off_floor > ASPECT_FLOOR,
            "VOID: every aspect ratio above {ASPECT_FLOOR} came from a cell sitting at H_FLOOR = \
             {H_FLOOR:e}, where aspect_ratio is |lambda|max / H_FLOOR and is a restatement of that \
             constant rather than a measurement (benches/common/metric.rs:67-74). The largest \
             off-floor ratio is {aspect_max_off_floor:e}"
        );

        let axis_ratio = results
            .iter()
            .flat_map(|result| result.rungs.iter())
            .map(|rung| rung.axis_ratio)
            .fold(0.0f64, f64::max);
        assert!(
            axis_ratio > AXIS_RATIO_FLOOR,
            "VOID: the most anisotropic grid this metric asked for anywhere is {axis_ratio:.4}:1, \
             below {AXIS_RATIO_FLOOR}:1. The two arms are then the same grid under two names and \
             every `ratio` in the file is a measurement of the extractor's run-to-run determinism"
        );

        for result in &results {
            for rung in &result.rungs {
                assert!(
                    rung.budget_ratio >= BUDGET_BAND[0] && rung.budget_ratio <= BUDGET_BAND[1],
                    "VOID: {} at {}^3 spent {:.4}x the isotropic arm's sample budget on \
                     {}x{}x{}, outside [{}, {}]. The arms are not budget-matched and a \
                     triangle-count comparison between them means nothing",
                    result.name,
                    rung.samples,
                    rung.budget_ratio,
                    rung.grid[0],
                    rung.grid[1],
                    rung.grid[2],
                    BUDGET_BAND[0],
                    BUDGET_BAND[1]
                );
            }
        }

        let population = results.iter().filter(|r| r.matched.is_ok()).count();
        assert!(
            population > 0,
            "VOID: no reference field carried a symmetric Hausdorff number that responded to \
             resolution, so C1 has no instrument at all rather than a narrow one. Reasons: {}",
            results
                .iter()
                .filter_map(|r| r
                    .matched
                    .as_ref()
                    .err()
                    .map(|why| format!("{}={why}", r.name)))
                .collect::<Vec<_>>()
                .join(" ")
        );

        // ── the global verdicts ──────────────────────────────────────────────
        let winners = results
            .iter()
            .filter(|r| r.matched.as_ref().is_ok_and(|m| m.win))
            .count();
        let c1_holds = if population >= C1_MIN_WINNERS {
            (winners >= C1_MIN_WINNERS).to_string()
        } else {
            format!("unreachable:population={population}<{C1_MIN_WINNERS}")
        };

        let c3_holds = results
            .iter()
            .flat_map(|result| result.rungs.iter())
            .all(|rung| rung.metric_ms.median / rung.extract_ms.median < C3_MAX_SHARE);

        let measurable: Vec<&FieldResult> = results.iter().filter(|r| r.matched.is_ok()).collect();
        let c2_matches = measurable
            .iter()
            .filter(|r| r.matched.as_ref().is_ok_and(|m| m.win) == r.predicted_win)
            .count();
        let c2_global = c2_matches == measurable.len();

        println!(
            "C1: population {population} of {FIELDS} -- only FieldBound::Exact admits a Hausdorff \
             measurement (fields/mod.rs:83-84); winners {winners}; bar {C1_MIN_WINNERS} -> \
             {c1_holds}"
        );
        println!(
            "C2: {c2_matches} of {} measurable fields matched their pre-registered prediction -> \
             {c2_global}",
            measurable.len()
        );
        println!("C3: every row under {C3_MAX_SHARE} of extraction -> {c3_holds}\n");

        // ── the rows ─────────────────────────────────────────────────────────
        for result in &results {
            let field_columns = match &result.matched {
                Ok(m) => [
                    format!("{:.3}", m.triangles_aniso),
                    format!("{:.6}", m.ratio),
                    m.win.to_string(),
                    (m.win == result.predicted_win).to_string(),
                    String::new(),
                    format!("{:.3}", m.triangles_iso),
                    format!("{:.9}", m.error),
                    format!("{:.6}", m.r2_iso),
                    format!("{:.6}", m.r2_aniso),
                    m.interpolated.to_string(),
                ],
                Err(why) => [
                    why.clone(),
                    why.clone(),
                    why.clone(),
                    String::from("unmeasurable"),
                    why.clone(),
                    why.clone(),
                    why.clone(),
                    why.clone(),
                    why.clone(),
                    why.clone(),
                ],
            };
            let [
                matched_triangles,
                ratio,
                c1_field_win,
                c2_holds,
                c1_skip,
                iso_at_matched,
                matched_error,
                r2_iso,
                r2_aniso,
                interpolated,
            ] = field_columns;

            for rung in &result.rungs {
                let share = rung.metric_ms.median / rung.extract_ms.median;
                let share_lo = rung.metric_ms.min / rung.extract_ms.max;
                let share_hi = rung.metric_ms.max / rung.extract_ms.min;
                let axis = flat_axis(&rung.census);
                let haus = |index: usize| {
                    rung.hausdorff
                        .map_or_else(str::to_string, |pair| format!("{:.9}", pair[index]))
                };

                run.record(&[
                    ("field", result.name.to_string()),
                    ("resolution", rung.samples.to_string()),
                    ("metric", METRIC_NAME.to_string()),
                    ("p_norm", format!("{P_NORM:.1}")),
                    (
                        "complexity_target",
                        format!("{:.6e}", rung.census.complexity),
                    ),
                    ("triangles_isotropic", rung.triangles_iso.to_string()),
                    ("triangles_anisotropic", rung.triangles_aniso.to_string()),
                    ("hausdorff_isotropic", haus(0)),
                    ("hausdorff_anisotropic", haus(1)),
                    ("triangles_at_matched_hausdorff", matched_triangles.clone()),
                    ("ratio", ratio.clone()),
                    (
                        "aspect_ratio_max",
                        format!("{:.6e}", rung.census.aspect_max),
                    ),
                    (
                        "aspect_ratio_mean",
                        format!("{:.6e}", rung.census.aspect_mean),
                    ),
                    ("metric_ms", format!("{:.6}", rung.metric_ms.median)),
                    ("metric_share", format!("{share:.6}")),
                    ("c1_holds", c1_holds.clone()),
                    ("c2_holds", c2_holds.clone()),
                    ("c3_holds", c3_holds.to_string()),
                    // ── extras (M-273) ──
                    (
                        "aspect_ratio_max_off_floor",
                        format!("{:.6e}", rung.census.aspect_max_off_floor),
                    ),
                    (
                        "aspect_ratio_mean_off_floor",
                        format!("{:.6e}", rung.census.aspect_mean_off_floor),
                    ),
                    ("at_floor_cells", rung.census.at_floor.to_string()),
                    (
                        "at_floor_fraction",
                        format!(
                            "{:.6}",
                            rung.census.at_floor as f64 / rung.census.points as f64
                        ),
                    ),
                    ("axes_pinned", rung.pinned.to_string()),
                    ("axis_ratio", format!("{:.6}", rung.axis_ratio)),
                    ("band_points", rung.census.points.to_string()),
                    ("budget_ratio", format!("{:.6}", rung.budget_ratio)),
                    ("c1_field_win", c1_field_win.clone()),
                    ("c1_population", population.to_string()),
                    ("c1_skip_reason", c1_skip.clone()),
                    ("c1_winners", winners.to_string()),
                    ("c2_global_holds", c2_global.to_string()),
                    (
                        "c2_mechanism_holds",
                        (axis == result.predicted_axis).to_string(),
                    ),
                    ("c2_predicted_axis", result.predicted_axis.to_string()),
                    ("c2_predicted_win", result.predicted_win.to_string()),
                    (
                        "c3_row_decisive",
                        ((share_lo < C3_MAX_SHARE) == (share_hi < C3_MAX_SHARE)).to_string(),
                    ),
                    ("c3_row_holds", (share < C3_MAX_SHARE).to_string()),
                    ("extract_ms", format!("{:.6}", rung.extract_ms.median)),
                    ("extract_ms_max", format!("{:.6}", rung.extract_ms.max)),
                    ("extract_ms_min", format!("{:.6}", rung.extract_ms.min)),
                    ("field_bound", result.bound.to_string()),
                    ("flat_axis", axis.to_string()),
                    (
                        "flat_axis_aligned_fraction",
                        format!(
                            "{:.6}",
                            rung.census.exploitable_flat as f64 / rung.census.points as f64
                        ),
                    ),
                    ("fit_r2_anisotropic", r2_aniso.clone()),
                    ("fit_r2_isotropic", r2_iso.clone()),
                    (
                        "grid_anisotropic",
                        format!("{}x{}x{}", rung.grid[0], rung.grid[1], rung.grid[2]),
                    ),
                    ("grid_isotropic", format!("{0}x{0}x{0}", rung.samples)),
                    ("hausdorff_matched", matched_error.clone()),
                    ("matched_interpolated", interpolated.clone()),
                    (
                        "mesh_residual_max_aniso",
                        format!("{:.6e}", rung.residual_aniso),
                    ),
                    (
                        "mesh_residual_max_iso",
                        format!("{:.6e}", rung.residual_iso),
                    ),
                    ("metric_ms_max", format!("{:.6}", rung.metric_ms.max)),
                    ("metric_ms_min", format!("{:.6}", rung.metric_ms.min)),
                    ("metric_share_hi", format!("{share_hi:.6}")),
                    ("metric_share_lo", format!("{share_lo:.6}")),
                    ("off_floor_cells", rung.census.off_floor.to_string()),
                    (
                        "share_times_resolution",
                        format!("{:.6}", share * f64::from(rung.samples)),
                    ),
                    ("triangles_iso_at_matched", iso_at_matched.clone()),
                ]);
            }
        }
    });
}
