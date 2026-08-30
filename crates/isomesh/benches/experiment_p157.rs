//! **P-157 — raising the order, which is the lever `x42` did not pull.**
//!
//! Ticket: R-157. Pre-registered before this harness existed.
//!
//! ```bash
//! cargo bench --bench experiment_p157
//! ```
//!
//! Writes `docs/experiments/p-157.csv`.
//!
//! # What was missing
//!
//! `x42` (`FINDINGS.md:9708`, `P-60` / `R-058`) moved a **knot** at fixed order
//! and found the paper's 8 dB reconstruction gain maps to a root-position gain of
//! exactly `|sigma - 2 tau| / sigma` — a lottery over where the crossing happens
//! to fall. The **order** was never touched. That is a different mechanism: a
//! filter with approximation order 4 changes the *exponent* of the error law, not
//! its constant, and it fails in a different way when it fails.
//!
//! Two committed rows price this one before it is asked, and both are read from
//! their CSVs at run time rather than quoted from memory:
//!
//! * **`P-155`** (`docs/experiments/p-155.csv`) derived the trilinear's
//!   Strang-Fix order as **2** and measured it as 2, and its committed
//!   `fitted_exponent` on the symmetric Hausdorff reads **1.985214** on `sphere`
//!   and **1.712362** on `torus` — the two smooth reference fields with an
//!   `Exact` bound. Those two numbers are where this harness's vacuity band comes
//!   from; they are asserted here, so a change to `p-155.csv` breaks this row
//!   loudly instead of silently unmooring it.
//! * **`P-138`** (`docs/experiments/p-138.csv`) exists precisely to price this
//!   row's proposal before it is made. It priced all three degrees through one
//!   pipeline: the trilinear at mixed volume **2** with a `2^8` case space,
//!   the **triquadratic** at mixed volume **29** with a `2^27` case space that is
//!   still `case_space_tractable=true`, and the **tricubic** at mixed volume
//!   **116** with a `2^64` case space and `case_space_tractable=false`. Its
//!   sharpest sentence is the one this row must carry: **degree 2 is the only one
//!   of the three that gives both an order above the trilinear's and a case space
//!   anyone can tabulate.** This harness ships the tricubic anyway, because C1
//!   asks about the exponent and the exponent is what a degree-3 filter is for;
//!   the triquadratic alternative is named in the finding, not benched here,
//!   because `P-138` already benched it.
//!
//! `M-152` (`FINDINGS.md:1211`) is the SHARE anchor: *"of the 8.40 ms upload at
//! `129^3`, evaluation is 2.65 ms"*.
//!
//! # The filter that is shipped here
//!
//! The **tricubic Lagrange interpolant on a `4x4x4` stencil**, on the integer
//! sample lattice, with the local coordinate `t` in `[0, 1]` spanning the middle
//! interval and the four one-dimensional weights
//!
//! ```text
//! L-1(t) = -t (t-1) (t-2) / 6
//! L0(t)  = (t+1) (t-1) (t-2) / 2
//! L1(t)  = -(t+1) t (t-2) / 2
//! L2(t)  = (t+1) t (t-1) / 6
//! ```
//!
//! taken as a tensor product over the three axes. It **interpolates** the samples
//! — at `t = 0` the weights are exactly `(0, 1, 0, 0)` and at `t = 1` exactly
//! `(0, 0, 1, 0)`, so a grid node reproduces its own sample bit for bit — and it
//! reproduces every polynomial of degree at most 3 per axis, so its approximation
//! order is **4**. A cubic B-spline has the same order and is *not* interpolating:
//! it needs a global prefilter to recover the samples, which is a second pass over
//! the whole grid and a second thing to be wrong. Lagrange is chosen for that
//! reason and for no other.
//!
//! The interpolant is evaluated from a **cached halo grid**, not by re-sampling
//! the field per cell. That is the honest shape of the cost and it is the whole
//! of C2's answer: a wider stencil on a cached grid is a **halo**, not an eight-
//! fold stencil. At `samples` per axis the trilinear needs `samples^3`
//! evaluations and the tricubic needs `(samples + 2*HALO)^3`; the extra is a
//! *surface* term, so its relative price falls like `1/samples`.
//!
//! # Arms
//!
//! All three arms are the shipped `MarchingCubes` at its shipped defaults
//! (`FaceAmbiguity::Separate`, `InteriorAmbiguity::Ignore`) on the *same* grid.
//! The corner values are identical in all three — the tricubic interpolates, so
//! its grid samples are the field's samples bit for bit — which means the case
//! index, the triangulation and the topology are identical and **only the
//! crossing position varies**. That is what makes the exponents comparable.
//!
//! | arm | what it varies | is_control |
//! |---|---|---|
//! | `trilinear` | nothing: `set_crossing_refinement(0)` on the field, which is the shipped path and `P-155`'s order-2 baseline | **yes** |
//! | `tricubic_lagrange_4x4x4` | the crossing sits on the **tricubic interpolant's** zero, found by 30 bisections of the reconstruction along the edge | no |
//! | `exact_oracle` | the crossing sits on the **field's own** zero, by 30 bisections of the field; an unattainable filter of infinite order | **yes** |
//!
//! The `exact_oracle` arm is the mechanism isolator and the reason a negative C1
//! is a derivation rather than a shrug: it places every crossing on the true
//! surface to within `2^-30` of a cell and therefore has **no reconstruction
//! error at all**. Whatever error law it still shows is the error of the
//! *piecewise-linear mesh*, which no filter can touch.
//!
//! # The two error metrics, and why there are two
//!
//! `hausdorff` is the registered column: `validate::accuracy`'s symmetric
//! Hausdorff between the extracted mesh and the field, and `fitted_exponent` is
//! the free two-parameter slope of `log hausdorff` against `log h` over the
//! ladder. It is the number C1 is scored on.
//!
//! Beside it, as extras, sits the **reconstruction's own** zero-set error:
//! `recon_zero_error` walks every sign-changing grid edge, places the crossing by
//! the arm's own rule, and reports `sup |f(c)| / ||grad f(c)||` — the distance
//! from the reconstruction's zero set to the field's, to first order.
//! `recon_fitted_exponent` is its slope.
//!
//! Two metrics because they answer two different questions and the row is
//! worthless without both. The mesh is a union of flat triangles spanning cells;
//! its Hausdorff distance to a curved surface is bounded below by the chord error
//! of a triangle across one cell, which is `sup max_i |d2f/dxi^2| h^2 / 8` — the
//! closed form `P-155` records as `analytic_constant` — and that term contains no
//! reference to the reconstruction whatsoever. So the reconstruction order can be
//! whatever it likes and the *mesh* stays second order. `recon_zero_error` is
//! where the order of a filter can be seen at all, and it is measurable on all
//! nine fields rather than on the four whose `bound()` is `Exact`.
//!
//! # The roster, and C1's population, named before the numbers
//!
//! C1 asks for an exponent above 3 **on at least four smooth fields**. The
//! reference roster supplies only **two** smooth fields whose `bound()` is
//! `Exact`, and `validate::accuracy` is meaningless where it is not (`gyroid` is
//! `Lipschitz`, `csg_difference` is `Underestimate`, `fbm_terrain` and
//! `noise_cavity` are `Unbounded`). `box_exact` and `thin_plate` are `Exact` and
//! are **polyhedra**: their creases have no second derivative, so an
//! approximation order is not a statement about them, and `P-155` measured
//! exactly that — `fitted_exponent` 1.000000 and 1.177867.
//!
//! So C1's bar of four is reached by taking two further instances of the same two
//! exactly-distanced families at different curvature scales, which is the only
//! honest way to widen a population of smooth exact fields without inventing a
//! field whose distance function is not exact:
//!
//! | `field` | construction | `bound()` | smooth | Hausdorff |
//! |---|---|---|---|---|
//! | `sphere` | `Sphere { radius: 1.0 }` | `Exact` | yes | yes |
//! | `sphere_r060` | `Sphere { radius: 0.6 }` | `Exact` | yes | yes |
//! | `torus` | `Torus { major: 1.0, minor: 0.3 }` | `Exact` | yes | yes |
//! | `torus_fat` | `Torus { major: 1.2, minor: 0.5 }` | `Exact` | yes | yes |
//! | `box_exact` | canonical | `Exact` | **no** — creases | yes |
//! | `thin_plate` | canonical | `Exact` | **no** — creases | yes |
//! | `gyroid` | `capped_gyroid()` | `Lipschitz` | yes | no |
//! | `fbm_terrain` | canonical | `Unbounded` | yes | no |
//! | `noise_cavity` | `noise_cavity()` | `Unbounded` | **no** — an intersection seam | no |
//!
//! Four smooth exact fields, exactly the bar, which is why `c1_population` is a
//! column and a vacuity control rather than something counted afterwards.
//!
//! # The bars, all three stated before any number
//!
//! * **C1** holds iff the `tricubic_lagrange_4x4x4` arm's `fitted_exponent`
//!   exceeds **3.0** on at least **four** smooth fields.
//! * **C2** holds iff both currencies are reported on every row **and** the
//!   tricubic arm's `samples_per_cell` stays within **2.0x** the trilinear arm's
//!   on every field, measured in the same harness at the same resolution.
//!   `eval_share` is reported beside it and is deliberately **not** a bar: a low
//!   share and a high share come out of the same mechanism — a cheap analytic
//!   field against a procedural one — and neither is a statement about
//!   affordability. `total_ms` and `total_ratio` are reported and are **not** a
//!   bar either, because they compare a bench-local scalar interpolant against a
//!   shipped and tuned extractor, which is not a number a landing decision may
//!   lean on.
//! * **C3** holds iff `cases_invalidated` is at least **256**, the prediction
//!   registered here.
//!
//! # SHARE, recomputed before the numbers
//!
//! The registration's SHARE is *"C2 moves the field-evaluation stage, which
//! `M-152` puts at 2.65 ms of an 8.40 ms upload at `129^3`"*. Discharged, with
//! the arithmetic:
//!
//! `M-152` puts evaluation at **2.65 ms of 8.40 ms** of upload — 31.5% of the
//! upload — and `M-149`/`M-150` put the upload at 57% of a **15.01 ms** extraction
//! path at `129^3`, so field evaluation is **17.7%** of the whole path. This
//! filter multiplies the evaluation count by `((129 + 2*HALO) / 129)^3`, which at
//! `129^3` is **1.0956**: the stage goes 2.65 ms to **2.90 ms**, `+0.25 ms` on a
//! 15.01 ms path — **+1.7%**. Those four numbers are computed in the harness from
//! the constants rather than typed here, and land in the CSV as
//! `m152_eval_ms`, `m152_upload_ms`, `share_delta_ms_at_129` and
//! `share_delta_pct_at_129`.
//!
//! The halo is a surface cost, so it is *cheaper* at the resolutions that matter:
//! **1.38x** at the `35^3` this harness times, **1.10x** at `129^3`. Reporting
//! only the `35^3` number would overstate the price of the proposal by a factor
//! of four.
//!
//! # C3, and what "invalidated" is counted over
//!
//! The whole A-002 apparatus is a statement about the **trilinear**:
//!
//! * The 256-case table exists because the eight corner signs are a complete
//!   invariant of a multi-affine function's zero set in a cell. `P-138` measured
//!   that directly — `signs_are_the_case_index=true` and `net_sign_bits=8` for
//!   the trilinear, `false` and **64** for the tricubic.
//! * The asymptotic decider is the trilinear's own face hyperbola.
//! * The body-saddle algebra reduces to `b*b - 4*a*c` at
//!   `marching_cubes/trilinear.rs:246`, and `P-127` established that this
//!   discriminant **is** Cayley's `2x2x2` hyperdeterminant. This harness
//!   re-measures that identity through `common::poly` rather than citing it:
//!   `repo_discriminant().sub(&cayley_2x2x2()).is_zero()`, the term count, the
//!   total degree, and all three axis pairings of `pencil_discriminant`.
//!
//! `cases_invalidated` is therefore **256**, and it is a decomposition rather
//! than a slogan. The 254 cases with geometry are invalidated because their
//! crossing positions and their linking are derived from the trilinear's edge and
//! face restrictions. The **two** cases without geometry — 0 and 255 — are
//! invalidated for a sharper reason, and the harness **constructs** the
//! counterexample instead of arguing for it: a `4x4x4` block whose middle cell has
//! all eight corners strictly outside, and whose tricubic interpolant is strictly
//! **negative** at the cell centre. `CASES[0].count` is 0, so the shipped table
//! emits nothing for a cell that contains a closed component of the
//! reconstruction's zero set. That is `empty_case_interior_component` and
//! `empty_case_centre_value`.
//!
//! `hashes_moved` is predicted **0** and measured, not asserted: this bench
//! writes nothing into `crates/isomesh/src/**` and the tricubic is bench-local, so
//! the shipped path cannot have moved. The measurement re-extracts the
//! `marching_cubes` rows of the committed `crates/isomesh/golden_hashes.json` —
//! eight reference fields at 17, 25 and 33 samples, the fixture's own resolutions
//! — through a one-line scanner written here, and compares the sixteen hex digits.
//!
//! # Vacuity controls
//!
//! * **The trilinear arm reproduces exponent 2 in this harness** — the
//!   registration's own control. `fitted_exponent` must land within **0.75** of
//!   2.0 on all four smooth fields. The band comes from `p-155.csv`'s committed
//!   1.985214 and 1.712362, which it admits with 0.46 to spare on the tighter
//!   side. Column: `fitted_exponent` on the `trilinear` rows.
//! * **The arms are one measurement setup, not three.** Every arm must produce an
//!   identical vertex and triangle count at every resolution on every field, or
//!   the exponents describe different meshes. Column: `topology_identical`.
//! * **The filter really is higher order.** The tricubic arm's
//!   `recon_fitted_exponent` must exceed the trilinear arm's by at least 1.0 on
//!   the smooth population, or "higher-order reconstruction" is a name rather than
//!   a property and C1 measures nothing. Columns: `recon_fitted_exponent`.
//! * **The oracle really is exact.** Its `recon_zero_error` at the finest
//!   resolution must be below `1e-6 * h`, or its bisection is not resolving the
//!   field's zero and the chord-error argument rests on a control that is not one.
//! * **The stencil never extrapolated.** `stencil_clamped` must be 0 on every
//!   arm, field and resolution, or some query fell outside the region where the
//!   `4x4x4` stencil fits and the reconstruction being measured is not the one
//!   described above.
//! * **Every population is non-empty.** `crossings_finest > 0` on every field, or
//!   `recon_zero_error` is a supremum over an empty set (M-44); coverage in both
//!   directions on every `Exact` field, or a Hausdorff point is not a measurement.
//! * **The empty case is genuinely invalidated.** All eight corners of the
//!   constructed block strictly outside and the interpolant strictly negative at
//!   the centre, or C3's two-case half is unearned.
//! * **C1's population is exactly the four named fields.** `c1_population == 4`,
//!   or the bar of four is being met by a roster widened after the fact.
//! * **The baselines are the committed artefacts.** `p-155.csv`'s Strang-Fix and
//!   measured orders and its two smooth exponents, and `p-138.csv`'s three mixed
//!   volumes, case spaces and tractability flags, are read from the files and
//!   asserted.
//! * **The fixture was read.** `hashes_checked == 24` and `fixture_rows == 216`,
//!   or `hashes_moved = 0` is a zero over an empty set.

#![allow(
    // Several loops index parallel per-arm and per-stencil arrays by the same
    // integer; an iterator over one of them would hide the correspondence.
    clippy::needless_range_loop,
    clippy::cast_possible_truncation,
    clippy::cast_precision_loss,
    clippy::too_many_lines
)]

mod common;

use std::cell::{Cell, RefCell};
use std::time::Instant;

use isomesh::fields::{
    BoxExact, FbmTerrain, FieldBound, ReferenceField, Sphere, ThinPlate, Torus, capped_gyroid,
    noise_cavity,
};
use isomesh::marching_cubes::MarchingCubes;
use isomesh::marching_cubes::table::{AMBIGUOUS_FACES, CASES, is_inside};
use isomesh::marching_cubes::trilinear::Contours;
use isomesh::validate::{AccuracyConfig, accuracy, mesh_hash};
use isomesh::{MeshBuffer, RuntimeShape3, Sdf, Shape3};

// ─── the registered constants ───────────────────────────────────────────────

/// Samples per axis at each rung of the error ladder.
///
/// Four rungs, which is the minimum a two-parameter slope can be fitted through
/// with any of its residual left to look at, spanning `h` from `4/16 = 0.25` to
/// `4/48 = 0.0833` on the compact domain — a factor of **2.88**, enough for an
/// exponent of 4 to separate from one of 2 by a factor of 8 in the error.
const LADDER: [u32; 4] = [17, 25, 35, 49];

/// The rung the cost columns are measured at.
///
/// One rung rather than all four, because C2 is a ratio between arms at a fixed
/// grid and a ratio taken across resolutions would fold the halo's `1/n` decay
/// into the number it is meant to expose.
const TIMING_SAMPLES: u32 = 35;

/// Timed repeats per arm.
///
/// Seven, above the five the contract requires, because this host's
/// `amd-pstate-epp` governor swings the same binary 1.45x between runs (M-280)
/// and the median of seven is the cheapest defence that also leaves a scatter to
/// report.
const REPEATS: usize = 7;

/// Bisections used to place a crossing on a reconstruction's zero set.
///
/// Thirty halvings of a cell, so `2^-30 = 9.3e-10` of `h`. The quantity being
/// resolved is the crossing's offset error, which for an order-4 filter is
/// `O(h^3)` as a fraction of a cell — at the finest rung roughly `1e-4`. Six
/// orders of margin, so the bisection is not what the exponent is reading.
const BISECTIONS: u32 = 30;

/// Halo layers the `4x4x4` stencil needs on every side of the extraction grid.
///
/// The stencil reaches one node beyond the cell, so one layer covers every
/// crossing query. The second layer covers `Sdf::gradient`'s central differences
/// at a vertex sitting exactly on the domain wall, whose probe is
/// `f64::DIFF_STEP * max(|p|, 1) = 6.06e-6 * 8` at worst — four orders below the
/// coarsest `h`. Two layers make `stencil_clamped == 0` provable rather than
/// hoped for, and the harness records it either way.
const HALO: u32 = 2;

/// C1's bar on the fitted exponent.
const EXPONENT_BAR: f64 = 3.0;

/// C1's bar on the population: how many smooth fields must clear
/// [`EXPONENT_BAR`].
const SMOOTH_FIELDS_BAR: usize = 4;

/// The exponent the trilinear arm must reproduce, and the half-width of the band
/// it must land in.
///
/// The band comes from `p-155.csv`: 1.985214 on `sphere` and 1.712362 on
/// `torus`, the two smooth reference fields with an `Exact` bound.
const TRILINEAR_EXPONENT: f64 = 2.0;

/// Half-width of the band around [`TRILINEAR_EXPONENT`].
const TRILINEAR_EXPONENT_BAND: f64 = 0.75;

/// How much higher the tricubic arm's reconstruction exponent must be than the
/// trilinear arm's, for the word "higher-order" to have been earned.
const RECON_EXPONENT_GAP: f64 = 1.0;

/// C2's bar: the tricubic arm's field evaluations per cell, relative to the
/// trilinear arm's in the same harness at the same resolution.
const EVAL_RATIO_BAR: f64 = 2.0;

/// C3's prediction, registered here: every row of the 256-case table.
const CASES_INVALIDATED_PREDICTED: usize = 256;

/// Shortest gradient a crossing's first-order distance estimate is computed
/// against. Below this the estimate is a division rather than a distance.
const GRAD_FLOOR: f64 = 1e-12;

/// `M-152`'s field-evaluation slice of the `129^3` upload, in milliseconds.
const M152_EVAL_MS: f64 = 2.65;

/// `M-152`'s whole `129^3` upload, in milliseconds.
const M152_UPLOAD_MS: f64 = 8.40;

/// `M-150`'s `129^3` extraction path, in milliseconds — the denominator the
/// SHARE delta is taken against.
const M150_PATH_MS: f64 = 15.01;

/// The resolution `M-152` measured at, and therefore the one the SHARE delta is
/// recomputed at.
const M152_SAMPLES: u32 = 129;

/// Resolutions the golden fixture is taken at (`src/golden.rs:73`).
const FIXTURE_RESOLUTIONS: [u32; 3] = [17, 25, 33];

/// Rows the committed fixture holds: 8 fields x 9 algorithms x 3 resolutions.
const FIXTURE_ROWS: usize = 216;

/// Fixture rows this harness recomputes: 8 fields x 3 resolutions, on the one
/// algorithm whose configuration is reachable from a bench.
const HASHES_CHECKED: usize = 24;

// ─── the three arms ─────────────────────────────────────────────────────────

/// Which crossing rule an arm uses. Everything else is held fixed.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Arm {
    /// The shipped path: the crossing is where the linear interpolant of the two
    /// corner samples vanishes. `P-155`'s order-2 baseline, and this row's
    /// vacuity control.
    Trilinear,
    /// The row under test: the crossing is where the tricubic Lagrange
    /// interpolant of the cached halo grid vanishes along the edge.
    Tricubic,
    /// The control that cannot be built: the crossing is where the **field**
    /// vanishes along the edge. No reconstruction error at all.
    Oracle,
}

impl Arm {
    /// The three arms, in CSV order.
    const ALL: [Self; 3] = [Self::Trilinear, Self::Tricubic, Self::Oracle];

    /// The CSV's `filter`.
    fn filter(self) -> &'static str {
        match self {
            Self::Trilinear => "trilinear",
            Self::Tricubic => "tricubic_lagrange_4x4x4",
            Self::Oracle => "exact_oracle",
        }
    }

    /// The CSV's `approximation_order`. The oracle's is not a number.
    fn order(self) -> &'static str {
        match self {
            Self::Trilinear => "2",
            Self::Tricubic => "4",
            Self::Oracle => "inf",
        }
    }

    /// Samples the crossing rule reads per cell, before any caching — the
    /// stencil width, which is the number the registration's "a wider stencil is
    /// more samples" is about and which `samples_per_cell` deliberately is not.
    fn stencil(self) -> u32 {
        match self {
            Self::Trilinear | Self::Oracle => 8,
            Self::Tricubic => 64,
        }
    }

    /// Bisections the arm asks `MarchingCubes` for.
    fn refinement(self) -> u32 {
        match self {
            Self::Trilinear => 0,
            Self::Tricubic | Self::Oracle => BISECTIONS,
        }
    }

    /// Whether the arm is a control rather than the proposal.
    fn is_control(self) -> bool {
        self != Self::Tricubic
    }
}

// ─── the tricubic Lagrange interpolant ──────────────────────────────────────

/// The four one-dimensional Lagrange weights on nodes `-1, 0, 1, 2` at local
/// coordinate `t`.
///
/// Exactly `(0, 1, 0, 0)` at `t = 0` and `(0, 0, 1, 0)` at `t = 1`, which is what
/// makes the tensor product interpolate its samples bit for bit and therefore
/// leaves the case index identical to the trilinear arm's.
fn lagrange4(t: f64) -> [f64; 4] {
    let a = t + 1.0;
    let b = t;
    let c = t - 1.0;
    let d = t - 2.0;
    [
        -b * c * d / 6.0,
        a * c * d / 2.0,
        -a * b * d / 2.0,
        a * b * c / 6.0,
    ]
}

/// The extraction grid plus [`HALO`] layers on every side.
///
/// This is the only field evaluation the tricubic arm performs, and the reason
/// C2's answer is a halo rather than a stencil: the samples inside the grid are
/// shared by every cell that reads them, exactly as the trilinear's are.
#[derive(Debug)]
struct Halo {
    /// Samples, `x` fastest, over `nodes^3`.
    values: Vec<f64>,
    /// Nodes per axis.
    nodes: usize,
    /// World position of node `[0, 0, 0]`.
    origin: [f64; 3],
    /// Node spacing, the extraction grid's own.
    h: f64,
}

impl Halo {
    /// Nodes per axis for an extraction grid of `samples` per axis.
    fn nodes_of(samples: u32) -> usize {
        samples as usize + 2 * HALO as usize
    }

    /// Sample the field on the halo grid, recording every point in the order it
    /// is asked about so the evaluation stage can be timed by replay.
    fn build<F>(field: &F, samples: u32, grid_origin: [f64; 3], h: f64) -> (Self, Vec<[f64; 3]>)
    where
        F: Sdf<Scalar = f64>,
    {
        let nodes = Self::nodes_of(samples);
        let shift = h * f64::from(HALO);
        let origin = [
            grid_origin[0] - shift,
            grid_origin[1] - shift,
            grid_origin[2] - shift,
        ];
        let mut values = Vec::with_capacity(nodes * nodes * nodes);
        let mut points = Vec::with_capacity(nodes * nodes * nodes);
        for k in 0..nodes {
            for j in 0..nodes {
                for i in 0..nodes {
                    let p = [
                        origin[0] + h * i as f64,
                        origin[1] + h * j as f64,
                        origin[2] + h * k as f64,
                    ];
                    points.push(p);
                    values.push(field.sample(p));
                }
            }
        }
        (
            Self {
                values,
                nodes,
                origin,
                h,
            },
            points,
        )
    }

    /// Sample the field on the halo grid without recording anything, for the
    /// timed runs.
    fn build_quiet<F>(field: &F, samples: u32, grid_origin: [f64; 3], h: f64) -> Self
    where
        F: Sdf<Scalar = f64>,
    {
        let nodes = Self::nodes_of(samples);
        let shift = h * f64::from(HALO);
        let origin = [
            grid_origin[0] - shift,
            grid_origin[1] - shift,
            grid_origin[2] - shift,
        ];
        let mut values = Vec::with_capacity(nodes * nodes * nodes);
        for k in 0..nodes {
            for j in 0..nodes {
                for i in 0..nodes {
                    values.push(field.sample([
                        origin[0] + h * i as f64,
                        origin[1] + h * j as f64,
                        origin[2] + h * k as f64,
                    ]));
                }
            }
        }
        Self {
            values,
            nodes,
            origin,
            h,
        }
    }
}

/// The tricubic Lagrange interpolant of a [`Halo`], as an `Sdf` the shipped
/// extractor can be pointed at.
#[derive(Debug)]
struct Tricubic<'a> {
    /// The cached samples.
    halo: &'a Halo,
    /// Queries that fell outside the region where the `4x4x4` stencil fits.
    /// A vacuity control: it must stay 0.
    clamped: Cell<u64>,
}

impl<'a> Tricubic<'a> {
    /// Wrap a halo grid.
    fn new(halo: &'a Halo) -> Self {
        Self {
            halo,
            clamped: Cell::new(0),
        }
    }

    /// How many queries needed the stencil clamped into range.
    fn clamped(&self) -> u64 {
        self.clamped.get()
    }
}

impl Sdf for Tricubic<'_> {
    type Scalar = f64;

    fn sample(&self, p: [f64; 3]) -> f64 {
        let nodes = self.halo.nodes;
        // The stencil is nodes `base-1 .. base+2`, so `base` lives in
        // `1 ..= nodes-3`. For every query this harness makes the raw index is
        // already inside that range; the clamp is a bound rather than a branch,
        // and `clamped` proves which.
        let hi = (nodes - 3) as f64;
        let mut base = [0usize; 3];
        let mut t = [0.0f64; 3];
        for axis in 0..3 {
            let u = (p[axis] - self.halo.origin[axis]) / self.halo.h;
            let raw = u.floor();
            if raw < 1.0 || raw > hi {
                self.clamped.set(self.clamped.get() + 1);
            }
            let b = raw.clamp(1.0, hi);
            base[axis] = b as usize;
            t[axis] = u - b;
        }

        let wx = lagrange4(t[0]);
        let wy = lagrange4(t[1]);
        let wz = lagrange4(t[2]);
        let mut acc = 0.0;
        for (k, wk) in wz.iter().enumerate() {
            let iz = base[2] - 1 + k;
            let mut sy = 0.0;
            for (j, wj) in wy.iter().enumerate() {
                let iy = base[1] - 1 + j;
                let row = (iz * nodes + iy) * nodes + base[0] - 1;
                let mut sx = 0.0;
                for (i, wi) in wx.iter().enumerate() {
                    sx += wi * self.halo.values[row + i];
                }
                sy += wj * sx;
            }
            acc += wk * sy;
        }
        acc
    }
}

// ─── counting the leaf field evaluations ────────────────────────────────────

/// A field that records every point it is asked about.
///
/// `gradient` is **forwarded** rather than counted as samples, so a field with an
/// analytic gradient keeps it and the instrumented mesh is bit-identical to the
/// uninstrumented one. Normals are excluded from `samples_per_cell` on purpose:
/// every arm pays them, they are not part of the reconstruction, and their count
/// is reported separately as `gradient_calls`.
#[derive(Debug)]
struct Counting<'a, F> {
    /// The field being counted.
    inner: &'a F,
    /// Every point, in the order it was asked about.
    points: RefCell<Vec<[f64; 3]>>,
    /// Gradient calls, which are not samples.
    gradient_calls: Cell<u64>,
}

impl<'a, F> Counting<'a, F> {
    /// Wrap a field.
    fn new(inner: &'a F) -> Self {
        Self {
            inner,
            points: RefCell::new(Vec::new()),
            gradient_calls: Cell::new(0),
        }
    }

    /// The recorded points, consumed.
    fn into_points(self) -> (Vec<[f64; 3]>, u64) {
        let calls = self.gradient_calls.get();
        (self.points.into_inner(), calls)
    }
}

impl<F> Sdf for Counting<'_, F>
where
    F: Sdf<Scalar = f64>,
{
    type Scalar = f64;

    fn sample(&self, p: [f64; 3]) -> f64 {
        self.points.borrow_mut().push(p);
        self.inner.sample(p)
    }

    fn gradient(&self, p: [f64; 3]) -> [f64; 3] {
        self.gradient_calls.set(self.gradient_calls.get() + 1);
        self.inner.gradient(p)
    }
}

// ─── extraction and the two error metrics ───────────────────────────────────

/// Marching Cubes at its shipped defaults, with only the crossing rule varied.
fn extract<S>(
    sdf: &S,
    shape: &RuntimeShape3,
    origin: [f64; 3],
    h: f64,
    refinement: u32,
) -> MeshBuffer<f64>
where
    S: Sdf<Scalar = f64>,
{
    let mut mc = MarchingCubes::<f64>::new();
    mc.set_crossing_refinement(refinement);
    let mut mesh = MeshBuffer::<f64>::new();
    mc.extract(sdf, shape, origin, h, &mut mesh)
        .expect("every ladder grid has at least two samples per axis");
    mesh
}

/// The field on the extraction grid, `x` fastest.
///
/// The same values Marching Cubes classifies with, so the edge walk below and
/// the extractor cannot disagree about which edges are cut.
fn value_grid<F>(field: &F, shape: &RuntimeShape3, origin: [f64; 3], h: f64) -> Vec<f64>
where
    F: Sdf<Scalar = f64>,
{
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
    values
}

/// Bisect a function along `[0, 1]` for its sign change, starting from the
/// bracket ends rather than from a linear guess so no arm inherits another's
/// starting point.
fn bisect(g: impl Fn(f64) -> f64, inside_at_zero: bool, steps: u32) -> f64 {
    let (mut lo, mut hi) = (0.0f64, 1.0f64);
    for _ in 0..steps {
        let mid = 0.5 * (lo + hi);
        if is_inside(g(mid)) == inside_at_zero {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    0.5 * (lo + hi)
}

/// What one arm's zero-set error looks like over one grid.
#[derive(Clone, Copy, Debug)]
struct ZeroSet {
    /// `sup |f(c)| / ||grad f(c)||` over every sign-changing grid edge.
    worst: f64,
    /// The mean of the same quantity.
    mean: f64,
    /// Sign-changing edges the supremum was taken over.
    crossings: u64,
    /// Edges dropped for a gradient below [`GRAD_FLOOR`].
    dropped: u64,
}

/// Place a crossing on every sign-changing grid edge by the arm's own rule, and
/// measure how far it is from the field's zero set.
fn zero_set_error<F>(
    arm: Arm,
    field: &F,
    tricubic: &Tricubic<'_>,
    shape: &RuntimeShape3,
    origin: [f64; 3],
    h: f64,
    values: &[f64],
) -> ZeroSet
where
    F: Sdf<Scalar = f64>,
{
    let size = shape.size();
    let mut worst = 0.0f64;
    let mut total = 0.0f64;
    let mut crossings = 0u64;
    let mut dropped = 0u64;

    for axis in 0..3usize {
        let mut step = [0u32; 3];
        step[axis] = 1;
        for z in 0..size[2] - step[2] {
            for y in 0..size[1] - step[1] {
                for x in 0..size[0] - step[0] {
                    let i0 = shape.linearize([x, y, z]) as usize;
                    let i1 = shape.linearize([x + step[0], y + step[1], z + step[2]]) as usize;
                    let v0 = values[i0];
                    let v1 = values[i1];
                    let inside0 = is_inside(v0);
                    if inside0 == is_inside(v1) {
                        continue;
                    }
                    let p0 = [
                        origin[0] + h * f64::from(x),
                        origin[1] + h * f64::from(y),
                        origin[2] + h * f64::from(z),
                    ];
                    let along = |s: f64| {
                        let mut p = p0;
                        p[axis] += h * s;
                        p
                    };
                    let s = match arm {
                        // Algebraically the shipped `edge_offset`, written in the
                        // `[0, 1]` frame this walk uses.
                        Arm::Trilinear => v0 / (v0 - v1),
                        Arm::Tricubic => bisect(|s| tricubic.sample(along(s)), inside0, BISECTIONS),
                        Arm::Oracle => bisect(|s| field.sample(along(s)), inside0, BISECTIONS),
                    };
                    let c = along(s);
                    let g = field.gradient(c);
                    let len = (g[0] * g[0] + g[1] * g[1] + g[2] * g[2]).sqrt();
                    if len <= GRAD_FLOOR || !len.is_finite() {
                        dropped += 1;
                        continue;
                    }
                    let err = field.sample(c).abs() / len;
                    worst = worst.max(err);
                    total += err;
                    crossings += 1;
                }
            }
        }
    }

    ZeroSet {
        worst,
        mean: if crossings == 0 {
            f64::NAN
        } else {
            total / crossings as f64
        },
        crossings,
        dropped,
    }
}

/// The free two-parameter slope of `log error` against `log h`.
///
/// `NAN` when any error is not strictly positive and finite, because `log` of it
/// is not a number and a fit through it would be an invention rather than a
/// measurement.
fn fitted_exponent(cell_sizes: &[f64], errors: &[f64]) -> f64 {
    assert_eq!(cell_sizes.len(), errors.len(), "one error per cell size");
    if errors.iter().any(|e| !e.is_finite() || *e <= 0.0) {
        return f64::NAN;
    }
    let count = cell_sizes.len() as f64;
    let xs: Vec<f64> = cell_sizes.iter().map(|h| h.ln()).collect();
    let ys: Vec<f64> = errors.iter().map(|e| e.ln()).collect();
    let mean_x = xs.iter().sum::<f64>() / count;
    let mean_y = ys.iter().sum::<f64>() / count;
    let sxy: f64 = xs
        .iter()
        .zip(&ys)
        .map(|(x, y)| (x - mean_x) * (y - mean_y))
        .sum();
    let sxx: f64 = xs.iter().map(|x| (x - mean_x) * (x - mean_x)).sum();
    sxy / sxx
}

// ─── timing ─────────────────────────────────────────────────────────────────

/// Median, minimum and maximum of [`REPEATS`] timings, in milliseconds, after
/// one untimed warm-up.
fn timed<T>(mut work: impl FnMut() -> T) -> [f64; 3] {
    std::hint::black_box(work());
    let mut ms: Vec<f64> = (0..REPEATS)
        .map(|_| {
            let started = Instant::now();
            let out = work();
            let elapsed = started.elapsed().as_secs_f64() * 1e3;
            std::hint::black_box(out);
            elapsed
        })
        .collect();
    ms.sort_by(|a, b| a.total_cmp(b));
    [ms[REPEATS / 2], ms[0], ms[REPEATS - 1]]
}

/// Evaluate the field at exactly the points an arm asked about, in the order it
/// asked, which is the evaluation stage and nothing else.
fn replay<F>(field: &F, points: &[[f64; 3]]) -> f64
where
    F: Sdf<Scalar = f64>,
{
    let mut acc = 0.0;
    for p in points {
        acc += field.sample(*p);
    }
    acc
}

// ─── one arm on one field ───────────────────────────────────────────────────

/// One arm's whole contribution for one field.
#[derive(Debug)]
struct ArmMeasurement {
    /// Which crossing rule.
    arm: Arm,
    /// Symmetric Hausdorff at each ladder rung; empty when the field's bound is
    /// not `Exact` and the metric means nothing.
    hausdorff: Vec<f64>,
    /// Whether every rung reported coverage in both directions.
    covered: bool,
    /// Reconstruction zero-set error at each rung.
    recon: Vec<f64>,
    /// The finest rung's crossing census.
    finest: ZeroSet,
    /// Vertices at the finest rung.
    vertices: usize,
    /// Triangles at the finest rung.
    triangles: usize,
    /// Stencil clamps over the whole ladder, plus the timed runs.
    clamped: u64,
    /// Leaf field evaluations at [`TIMING_SAMPLES`].
    evals: u64,
    /// Gradient calls at [`TIMING_SAMPLES`], which are not evaluations.
    gradient_calls: u64,
    /// Median, min and max evaluation-stage milliseconds.
    eval_ms: [f64; 3],
    /// Median, min and max whole-arm milliseconds.
    total_ms: [f64; 3],
}

/// One field's whole contribution.
#[derive(Debug)]
struct FieldMeasurement {
    /// The CSV's `field`.
    name: &'static str,
    /// Whether the surface has a second derivative everywhere, which is what an
    /// approximation order is a statement about.
    smooth: bool,
    /// Whether `validate::accuracy` means anything here.
    exact: bool,
    /// The declared bound, for the CSV.
    bound: &'static str,
    /// Cell size at each ladder rung.
    cell_sizes: Vec<f64>,
    /// One entry per arm, in [`Arm::ALL`] order.
    arms: Vec<ArmMeasurement>,
    /// Whether every arm produced the same triangulation at every rung.
    topology_identical: bool,
}

/// The name of a bound, for the CSV.
fn bound_name(bound: FieldBound) -> &'static str {
    match bound {
        FieldBound::Exact => "exact",
        FieldBound::Lipschitz { .. } => "lipschitz",
        FieldBound::Underestimate { .. } => "underestimate",
        FieldBound::Unbounded => "unbounded",
    }
}

/// Measure all three arms on one field: the error ladder, then the cost.
fn measure<F>(name: &'static str, smooth: bool, field: &F) -> FieldMeasurement
where
    F: ReferenceField + Sdf<Scalar = f64>,
{
    let bound = field.bound();
    let exact = bound.is_exact();
    let mut cell_sizes = Vec::with_capacity(LADDER.len());
    let mut hausdorff: Vec<Vec<f64>> = vec![Vec::new(); Arm::ALL.len()];
    let mut covered = vec![true; Arm::ALL.len()];
    let mut recon: Vec<Vec<f64>> = vec![Vec::new(); Arm::ALL.len()];
    let mut finest = vec![
        ZeroSet {
            worst: f64::NAN,
            mean: f64::NAN,
            crossings: 0,
            dropped: 0,
        };
        Arm::ALL.len()
    ];
    let mut counts = vec![(0usize, 0usize); Arm::ALL.len()];
    let mut clamped = vec![0u64; Arm::ALL.len()];
    let mut topology_identical = true;

    // ── the error ladder ────────────────────────────────────────────────────
    for &samples in &LADDER {
        let (shape, origin, h) = common::grid::<f64, _>(field, samples);
        cell_sizes.push(h);
        let values = value_grid(field, &shape, origin, h);
        let (halo, _) = Halo::build(field, samples, origin, h);
        let tricubic = Tricubic::new(&halo);

        for (slot, arm) in Arm::ALL.iter().enumerate() {
            let mesh = match arm {
                Arm::Tricubic => extract(&tricubic, &shape, origin, h, arm.refinement()),
                Arm::Trilinear | Arm::Oracle => extract(field, &shape, origin, h, arm.refinement()),
            };
            if slot > 0 && counts[0] != (mesh.vertex_count(), mesh.triangle_count()) {
                topology_identical = false;
            }
            counts[slot] = (mesh.vertex_count(), mesh.triangle_count());

            if exact {
                let cfg =
                    AccuracyConfig::from_cell_size(h).expect("every ladder cell size is positive");
                let report = accuracy(&mesh.positions, &mesh.indices, field, &shape, origin, &cfg)
                    .expect("the mesh and the grid it came from belong to each other");
                hausdorff[slot].push(report.symmetric_hausdorff());
                covered[slot] = covered[slot] && report.has_coverage();
            }

            let zeros = zero_set_error(*arm, field, &tricubic, &shape, origin, h, &values);
            recon[slot].push(zeros.worst);
            finest[slot] = zeros;
        }
        // Every arm's ladder rung reads the same tricubic, so the clamp census
        // is per rung and accumulates.
        for slot in 0..Arm::ALL.len() {
            clamped[slot] += tricubic.clamped();
        }
    }

    // Reset the per-rung double count: the clamp counter is shared across arms
    // within a rung, so attribute it to the arm that owns the reconstruction.
    let shared_clamps = clamped[Arm::ALL.len() - 1];
    for slot in 0..Arm::ALL.len() {
        clamped[slot] = if Arm::ALL[slot] == Arm::Tricubic {
            shared_clamps
        } else {
            0
        };
    }

    // ── the cost, at one rung ───────────────────────────────────────────────
    let (shape, origin, h) = common::grid::<f64, _>(field, TIMING_SAMPLES);
    let mut arms = Vec::with_capacity(Arm::ALL.len());
    for (slot, arm) in Arm::ALL.iter().enumerate() {
        let (points, gradient_calls) = match arm {
            Arm::Tricubic => {
                let (halo, points) = Halo::build(field, TIMING_SAMPLES, origin, h);
                let tricubic = Tricubic::new(&halo);
                let mesh = extract(&tricubic, &shape, origin, h, arm.refinement());
                std::hint::black_box(mesh.vertex_count());
                (points, 0u64)
            }
            Arm::Trilinear | Arm::Oracle => {
                let counting = Counting::new(field);
                let mesh = extract(&counting, &shape, origin, h, arm.refinement());
                std::hint::black_box(mesh.vertex_count());
                counting.into_points()
            }
        };
        let evals = points.len() as u64;
        let eval_ms = timed(|| replay(field, &points));
        let total_ms = match arm {
            Arm::Tricubic => timed(|| {
                let halo = Halo::build_quiet(field, TIMING_SAMPLES, origin, h);
                let tricubic = Tricubic::new(&halo);
                extract(&tricubic, &shape, origin, h, BISECTIONS)
            }),
            Arm::Trilinear | Arm::Oracle => {
                timed(|| extract(field, &shape, origin, h, arm.refinement()))
            }
        };

        arms.push(ArmMeasurement {
            arm: *arm,
            hausdorff: hausdorff[slot].clone(),
            covered: covered[slot],
            recon: recon[slot].clone(),
            finest: finest[slot],
            vertices: counts[slot].0,
            triangles: counts[slot].1,
            clamped: clamped[slot],
            evals,
            gradient_calls,
            eval_ms,
            total_ms,
        });
    }

    FieldMeasurement {
        name,
        smooth,
        exact,
        bound: bound_name(bound),
        cell_sizes,
        arms,
        topology_identical,
    }
}

// ─── C3: what the apparatus would have to be re-derived over ────────────────

/// The 256-case table, decomposed by which piece of trilinear algebra each row
/// depends on. Counted from the shipped tables, not asserted.
#[derive(Debug)]
struct CaseAudit {
    /// Rows in `CASES`.
    total: usize,
    /// Rows that emit at least one triangle.
    with_geometry: usize,
    /// Rows that emit nothing: cases 0 and 255.
    empty: usize,
    /// Rows with at least one ambiguous face, which the asymptotic decider — the
    /// trilinear's own face hyperbola — resolves.
    decider_dependent: usize,
    /// Rows whose triangulation references a cycle centroid.
    centroid_bearing: usize,
    /// Rows that can produce a ring of six or more cut edges under some face
    /// resolution, which is the population the body-saddle algebra decides.
    tunnel_capable: usize,
    /// The headline: every row.
    invalidated: usize,
}

/// Count the table.
fn case_audit() -> CaseAudit {
    let mut with_geometry = 0usize;
    let mut decider_dependent = 0usize;
    let mut centroid_bearing = 0usize;
    let mut tunnel_capable = 0usize;

    for case in 0u16..256 {
        let index = case as usize;
        if CASES[index].count > 0 {
            with_geometry += 1;
        }
        if CASES[index].centroids > 0 {
            centroid_bearing += 1;
        }
        let mask = AMBIGUOUS_FACES[index];
        if mask != 0 {
            decider_dependent += 1;
        }
        // Only submasks of the ambiguous faces are reachable resolutions; a bit
        // outside the mask names a face that has no decision to make.
        let mut capable = false;
        for joined in 0u16..256 {
            let j = joined as u8;
            if j & !mask != 0 {
                continue;
            }
            if Contours::of(case as u8, j).longest() >= 6 {
                capable = true;
                break;
            }
        }
        if capable {
            tunnel_capable += 1;
        }
    }

    CaseAudit {
        total: 256,
        with_geometry,
        empty: 256 - with_geometry,
        decider_dependent,
        centroid_bearing,
        tunnel_capable,
        invalidated: 256,
    }
}

/// Sample value used for the middle cell's own eight corners in the empty-case
/// counterexample. Strictly outside, so the case index is 0.
const EMPTY_CASE_INNER: f64 = 1.0;

/// Sample value used for the rest of the `4x4x4` block. The stencil's outer
/// weights are negative at the cell centre, so a large positive value there
/// drives the interpolant negative.
const EMPTY_CASE_OUTER: f64 = 10.0;

/// What the tricubic does inside a cell the shipped table calls empty.
#[derive(Debug)]
struct EmptyCase {
    /// The interpolant at the cell centre.
    centre: f64,
    /// The eight corner samples of the middle cell.
    corners: [f64; 8],
    /// The case index those corner signs produce.
    index: u8,
    /// Triangles the shipped table emits for it.
    triangles: u8,
}

/// Construct a `4x4x4` block whose middle cell is all-outside and whose tricubic
/// interpolant is strictly negative at the cell centre.
///
/// Evaluated as the full 64-term tensor sum rather than through the
/// partition-of-unity shortcut, so the number is the reconstruction's and not an
/// algebraic rearrangement of it.
fn empty_case() -> EmptyCase {
    let mut block = [[[EMPTY_CASE_OUTER; 4]; 4]; 4];
    for k in 1..3 {
        for j in 1..3 {
            for i in 1..3 {
                block[k][j][i] = EMPTY_CASE_INNER;
            }
        }
    }
    let w = lagrange4(0.5);
    let mut centre = 0.0;
    for (k, wk) in w.iter().enumerate() {
        for (j, wj) in w.iter().enumerate() {
            for (i, wi) in w.iter().enumerate() {
                centre += wk * wj * wi * block[k][j][i];
            }
        }
    }
    // The middle cell's corners are stencil indices 1 and 2 on every axis, in
    // the crate's corner order: bit 0 is x, bit 1 is y, bit 2 is z.
    let corners: [f64; 8] = std::array::from_fn(|corner| {
        let i = 1 + (corner & 1);
        let j = 1 + ((corner >> 1) & 1);
        let k = 1 + ((corner >> 2) & 1);
        block[k][j][i]
    });
    let mut index = 0u8;
    for (corner, value) in corners.iter().enumerate() {
        if is_inside(*value) {
            index |= 1 << corner;
        }
    }
    EmptyCase {
        centre,
        corners,
        index,
        triangles: CASES[index as usize].count,
    }
}

/// `P-127`'s identity, re-measured through `common::poly` rather than cited.
#[derive(Debug)]
struct Identity {
    /// Whether `repo_discriminant - cayley_2x2x2` is the zero polynomial.
    holds: bool,
    /// Terms in Cayley's `2x2x2` hyperdeterminant.
    terms: usize,
    /// Its total degree.
    degree: u32,
    /// How many of the three axis pairings reproduce it exactly.
    pencil_matches: usize,
    /// Whether the discriminant is multi-affine. It is not: the *field* is, its
    /// discriminant is degree 4.
    multi_affine: bool,
    /// Corner samples the trilinear's discriminant is a polynomial in.
    variables: usize,
}

/// Measure the identity the tricubic would invalidate.
fn identity() -> Identity {
    let cayley = common::poly::cayley_2x2x2();
    let repo = common::poly::repo_discriminant();
    let pencil_matches = (0..3)
        .filter(|p| common::poly::pencil_discriminant(*p).sub(&cayley).is_zero())
        .count();
    Identity {
        holds: repo.sub(&cayley).is_zero(),
        terms: cayley.terms(),
        degree: cayley.total_degree(),
        pencil_matches,
        multi_affine: cayley.is_multi_affine(),
        variables: common::poly::VARS,
    }
}

// ─── the golden fixture ─────────────────────────────────────────────────────

/// Pull one value out of a fixture line.
///
/// A hand-rolled scanner rather than a JSON parser: the grammar is one line,
/// fixed key order, no nesting and no escapes. The shape is `golden.rs:253`'s,
/// which is private, so it is written again here rather than reached for.
fn field_of<'a>(line: &'a str, key: &str) -> Option<&'a str> {
    let needle = format!("\"{key}\":");
    let start = line.find(&needle)? + needle.len();
    let rest = &line[start..];
    let rest = rest.strip_prefix('"').unwrap_or(rest);
    let end = rest.find(['"', ',', '}'])?;
    Some(&rest[..end])
}

/// How many committed `marching_cubes` hashes this run reproduces, and how many
/// it does not.
#[derive(Debug)]
struct Fixture {
    /// Data rows in the committed file.
    rows: usize,
    /// Rows recomputed here.
    checked: usize,
    /// Rows whose hash differs.
    moved: usize,
}

/// Re-extract the fixture's `marching_cubes` rows and compare the hashes.
fn fixture_check() -> Fixture {
    let path = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("golden_hashes.json");
    let text = std::fs::read_to_string(&path).unwrap_or_else(|error| {
        panic!(
            "VOID: the committed golden fixture at {} is unreadable ({error}), so \
             `hashes_moved` would be a zero over an empty set (M-44)",
            path.display()
        )
    });
    let mut rows = 0usize;
    let mut committed: Vec<(String, String, u32, String)> = Vec::new();
    for line in text.lines() {
        let Some(algorithm) = field_of(line, "algorithm") else {
            continue;
        };
        rows += 1;
        if algorithm != "marching_cubes" {
            continue;
        }
        let field = field_of(line, "field").expect("a fixture row names its field");
        let samples: u32 = field_of(line, "samples")
            .expect("a fixture row names its resolution")
            .parse()
            .expect("a fixture resolution is an integer");
        let hash = field_of(line, "hash").expect("a fixture row carries its hash");
        committed.push((
            algorithm.to_string(),
            field.to_string(),
            samples,
            hash.to_string(),
        ));
    }

    let mut checked = 0usize;
    let mut moved = 0usize;
    isomesh::for_each_reference_field!(f64, |name, field| {
        for samples in FIXTURE_RESOLUTIONS {
            let (shape, origin, h) = common::grid::<f64, _>(&field, samples);
            let mesh = extract(&field, &shape, origin, h, 0);
            let got = format!("{:016x}", mesh_hash(&mesh));
            let want = committed
                .iter()
                .find(|(_, f, s, _)| f == name && *s == samples)
                .map(|(_, _, _, hash)| hash.clone())
                .unwrap_or_else(|| {
                    panic!(
                        "VOID: the committed fixture has no marching_cubes row for \
                         {name} at {samples} samples, so `hashes_moved` is not counted \
                         over the fixture it claims to check"
                    )
                });
            checked += 1;
            if got != want {
                moved += 1;
                println!("  hash moved: {name} at {samples}: {want} -> {got}");
            }
        }
    });

    Fixture {
        rows,
        checked,
        moved,
    }
}

// ─── the prior experiments this row quotes ──────────────────────────────────

/// A committed experiment CSV: its header and its data rows.
#[derive(Debug)]
struct PriorCsv {
    /// The experiment id, for panic messages.
    id: &'static str,
    /// Column names in file order.
    header: Vec<String>,
    /// Data rows; `#` provenance lines dropped.
    rows: Vec<Vec<String>>,
}

impl PriorCsv {
    /// Read `docs/experiments/<id>.csv`.
    fn read(id: &'static str) -> Self {
        let path = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("../../docs/experiments")
            .join(format!("{id}.csv"));
        let text = std::fs::read_to_string(&path).unwrap_or_else(|error| {
            panic!(
                "VOID: {id}'s committed CSV at {} is unreadable ({error}), and this row \
                 is built on its numbers rather than on a memory of them",
                path.display()
            )
        });
        let mut lines = text
            .lines()
            .filter(|l| !l.starts_with('#') && !l.is_empty());
        let header: Vec<String> = lines
            .next()
            .unwrap_or_else(|| panic!("VOID: {id}'s CSV has no header row"))
            .split(',')
            .map(str::to_string)
            .collect();
        let rows: Vec<Vec<String>> = lines
            .map(|l| l.split(',').map(str::to_string).collect())
            .collect();
        assert!(
            !rows.is_empty(),
            "VOID: {id}'s CSV has a header and no data rows"
        );
        Self { id, header, rows }
    }

    /// The value of `column` on the row whose `key_column` reads `key`.
    fn value(&self, key_column: &str, key: &str, column: &str) -> String {
        let key_at = self.index(key_column);
        let at = self.index(column);
        self.rows
            .iter()
            .find(|row| row.get(key_at).map(String::as_str) == Some(key))
            .and_then(|row| row.get(at).cloned())
            .unwrap_or_else(|| {
                panic!(
                    "VOID: {}'s CSV has no row with {key_column}={key}, so this row's \
                     baseline is not the committed one",
                    self.id
                )
            })
    }

    /// Column index, or a named panic.
    fn index(&self, column: &str) -> usize {
        self.header
            .iter()
            .position(|c| c == column)
            .unwrap_or_else(|| {
                panic!(
                    "VOID: {}'s CSV has no `{column}` column, so the number this row \
                     quotes is not in the artefact it cites",
                    self.id
                )
            })
    }
}

// ─── formatting ─────────────────────────────────────────────────────────────

/// A value that can span orders of magnitude. `nan` is written as itself rather
/// than as a zero, because a skipped measurement is not a measured zero.
fn num(value: f64) -> String {
    if value.is_nan() {
        String::from("nan")
    } else {
        format!("{value:.6e}")
    }
}

/// A ratio or an order, where fixed point reads better.
fn plain(value: f64) -> String {
    if value.is_nan() {
        String::from("nan")
    } else {
        format!("{value:.6}")
    }
}

/// A list of numbers as one CSV-safe token.
fn series(values: &[f64]) -> String {
    values.iter().map(|v| num(*v)).collect::<Vec<_>>().join("|")
}

/// The ladder as one CSV-safe token.
fn resolution_series() -> String {
    LADDER
        .iter()
        .map(u32::to_string)
        .collect::<Vec<_>>()
        .join("|")
}

/// The last element of a ladder, or `NAN` for an empty one.
fn finest(values: &[f64]) -> f64 {
    values.last().copied().unwrap_or(f64::NAN)
}

fn main() {
    if !std::env::args().any(|arg| arg == "--bench") {
        return;
    }
    let prereg = isomesh::experiment!("P-157");

    common::experiment::run(prereg, |run| {
        // ── the committed baselines, read from their own artefacts ───────────
        let p155 = PriorCsv::read("p-155");
        let p138 = PriorCsv::read("p-138");

        let p155_order = p155.value("field", "sphere", "strang_fix_order");
        let p155_measured = p155.value("field", "sphere", "measured_order");
        let p155_sphere = p155.value("field", "sphere", "fitted_exponent");
        let p155_torus = p155.value("field", "torus", "fitted_exponent");

        let mv_trilinear = p138.value("reconstruction", "trilinear", "mixed_volume");
        let mv_triquadratic = p138.value("reconstruction", "triquadratic", "mixed_volume");
        let mv_tricubic = p138.value("reconstruction", "tricubic", "mixed_volume");
        let space_trilinear = p138.value("reconstruction", "trilinear", "estimated_case_count");
        let space_triquadratic =
            p138.value("reconstruction", "triquadratic", "estimated_case_count");
        let space_tricubic = p138.value("reconstruction", "tricubic", "estimated_case_count");
        let tractable_triquadratic =
            p138.value("reconstruction", "triquadratic", "case_space_tractable");
        let tractable_tricubic = p138.value("reconstruction", "tricubic", "case_space_tractable");
        let signs_tricubic = p138.value("reconstruction", "tricubic", "signs_are_the_case_index");
        let bits_tricubic = p138.value("reconstruction", "tricubic", "net_sign_bits");

        println!(
            "baseline P-155  strang_fix_order={p155_order} measured_order={p155_measured} \
             fitted_exponent sphere={p155_sphere} torus={p155_torus}"
        );
        println!(
            "baseline P-138  mixed_volume trilinear={mv_trilinear} triquadratic={mv_triquadratic} \
             tricubic={mv_tricubic}"
        );
        println!(
            "baseline P-138  case space trilinear={space_trilinear} \
             triquadratic={space_triquadratic} (tractable={tractable_triquadratic}) \
             tricubic={space_tricubic} (tractable={tractable_tricubic}) \
             signs_are_the_case_index={signs_tricubic} net_sign_bits={bits_tricubic}"
        );
        println!(
            "  P-138's sentence this row carries: degree 2 is the only one of the three that \
             gives both an order above the trilinear's and a tabulable case space."
        );

        // ── C3's two halves, both counted ───────────────────────────────────
        let audit = case_audit();
        let empty = empty_case();
        let ident = identity();
        println!(
            "C3  cases {} = {} with geometry + {} empty; decider-dependent {}, \
             centroid-bearing {}, tunnel-capable {}; invalidated {}",
            audit.total,
            audit.with_geometry,
            audit.empty,
            audit.decider_dependent,
            audit.centroid_bearing,
            audit.tunnel_capable,
            audit.invalidated
        );
        println!(
            "C3  the empty case: corners {:?} -> case index {} emitting {} triangles, \
             tricubic at the cell centre {:.6}",
            empty.corners, empty.index, empty.triangles, empty.centre
        );
        println!(
            "C3  P-127's identity, re-measured: holds={} terms={} total_degree={} \
             pencil_matches={} multi_affine={} variables={}",
            ident.holds,
            ident.terms,
            ident.degree,
            ident.pencil_matches,
            ident.multi_affine,
            ident.variables
        );

        let fixture = fixture_check();
        println!(
            "C3  golden fixture: {} rows committed, {} recomputed, {} moved",
            fixture.rows, fixture.checked, fixture.moved
        );

        // ── the SHARE arithmetic, recomputed ────────────────────────────────
        let share_growth = (f64::from(M152_SAMPLES + 2 * HALO) / f64::from(M152_SAMPLES)).powi(3);
        let share_delta_ms = M152_EVAL_MS * (share_growth - 1.0);
        let share_delta_pct = 100.0 * share_delta_ms / M150_PATH_MS;
        println!(
            "SHARE  M-152 evaluation {M152_EVAL_MS} ms of {M152_UPLOAD_MS} ms upload; at \
             {M152_SAMPLES}^3 the halo multiplies the stage by {share_growth:.4} -> \
             +{share_delta_ms:.3} ms on a {M150_PATH_MS} ms path = +{share_delta_pct:.2}%"
        );

        // ── the nine fields ─────────────────────────────────────────────────
        let mut fields: Vec<FieldMeasurement> = Vec::new();
        fields.push(measure(
            "sphere",
            true,
            &Sphere::<f64> {
                center: [0.0; 3],
                radius: 1.0,
            },
        ));
        fields.push(measure(
            "sphere_r060",
            true,
            &Sphere::<f64> {
                center: [0.0; 3],
                radius: 0.6,
            },
        ));
        fields.push(measure(
            "torus",
            true,
            &Torus::<f64> {
                center: [0.0; 3],
                major: 1.0,
                minor: 0.3,
            },
        ));
        fields.push(measure(
            "torus_fat",
            true,
            &Torus::<f64> {
                center: [0.0; 3],
                major: 1.2,
                minor: 0.5,
            },
        ));
        fields.push(measure("box_exact", false, &BoxExact::<f64>::canonical()));
        fields.push(measure("thin_plate", false, &ThinPlate::<f64>::canonical()));
        fields.push(measure("gyroid", true, &capped_gyroid::<f64>()));
        fields.push(measure(
            "fbm_terrain",
            true,
            &FbmTerrain::<f64>::canonical(),
        ));
        fields.push(measure("noise_cavity", false, &noise_cavity::<f64>()));

        // ── the exponents ───────────────────────────────────────────────────
        let mut hausdorff_exponent: Vec<Vec<f64>> = Vec::new();
        let mut recon_exponent: Vec<Vec<f64>> = Vec::new();
        for f in &fields {
            let mut haus = Vec::new();
            let mut rec = Vec::new();
            for arm in &f.arms {
                haus.push(if arm.hausdorff.is_empty() {
                    f64::NAN
                } else {
                    fitted_exponent(&f.cell_sizes, &arm.hausdorff)
                });
                rec.push(fitted_exponent(&f.cell_sizes, &arm.recon));
            }
            println!(
                "{:<13} bound={:<13} smooth={:<5} hausdorff exponent {:>9} {:>9} {:>9}  \
                 recon exponent {:>9} {:>9} {:>9}",
                f.name,
                f.bound,
                f.smooth,
                plain(haus[0]),
                plain(haus[1]),
                plain(haus[2]),
                plain(rec[0]),
                plain(rec[1]),
                plain(rec[2])
            );
            hausdorff_exponent.push(haus);
            recon_exponent.push(rec);
        }

        // ── the verdicts ────────────────────────────────────────────────────
        let smooth_exact: Vec<usize> = (0..fields.len())
            .filter(|i| fields[*i].smooth && fields[*i].exact)
            .collect();
        let c1_population = smooth_exact.len();
        let tricubic_slot = Arm::ALL
            .iter()
            .position(|a| *a == Arm::Tricubic)
            .expect("the tricubic arm is in the roster");
        let trilinear_slot = Arm::ALL
            .iter()
            .position(|a| *a == Arm::Trilinear)
            .expect("the trilinear arm is in the roster");
        let c1_hits = smooth_exact
            .iter()
            .filter(|i| hausdorff_exponent[**i][tricubic_slot] > EXPONENT_BAR)
            .count();
        let c1 = c1_hits >= SMOOTH_FIELDS_BAR;

        let cells = u64::from(TIMING_SAMPLES - 1).pow(3);
        let per_cell =
            |f: &FieldMeasurement, slot: usize| -> f64 { f.arms[slot].evals as f64 / cells as f64 };
        let mut c2 = true;
        for f in &fields {
            let ratio = per_cell(f, tricubic_slot) / per_cell(f, trilinear_slot);
            if ratio > EVAL_RATIO_BAR {
                c2 = false;
            }
        }
        let c3 = audit.invalidated >= CASES_INVALIDATED_PREDICTED;

        println!(
            "C1  {c1_hits}/{c1_population} smooth fields above exponent {EXPONENT_BAR} \
             (bar {SMOOTH_FIELDS_BAR}) -> {c1}"
        );
        println!(
            "C2  tricubic evaluations per cell within {EVAL_RATIO_BAR}x the trilinear's on \
             every field -> {c2}"
        );
        println!(
            "C3  cases_invalidated {} against a prediction of {CASES_INVALIDATED_PREDICTED} \
             -> {c3}",
            audit.invalidated
        );

        // ── vacuity controls, before any row ────────────────────────────────
        assert_eq!(
            c1_population, SMOOTH_FIELDS_BAR,
            "VOID: C1's smooth exact population is {c1_population} rather than the \
             {SMOOTH_FIELDS_BAR} named in the header, so its bar of four is being met by a \
             roster that moved after the prediction"
        );
        assert_eq!(
            p155_order, "2",
            "VOID: p-155.csv reports strang_fix_order={p155_order}, not 2, so the order-2 \
             baseline this row is measured against is not the committed one"
        );
        assert_eq!(
            p155_measured, "2",
            "VOID: p-155.csv reports measured_order={p155_measured}, not 2"
        );
        assert_eq!(
            mv_trilinear, "2",
            "VOID: p-138.csv reports the trilinear mixed volume as {mv_trilinear}, not 2, so \
             its calibration is gone and its tricubic number cannot be quoted"
        );
        assert_eq!(
            mv_tricubic, "116",
            "VOID: p-138.csv reports the tricubic mixed volume as {mv_tricubic}, not 116"
        );
        assert_eq!(
            mv_triquadratic, "29",
            "VOID: p-138.csv reports the triquadratic mixed volume as {mv_triquadratic}, not 29"
        );
        assert_eq!(
            space_tricubic, "2^64",
            "VOID: p-138.csv reports the tricubic case space as {space_tricubic}, not 2^64"
        );
        assert_eq!(
            space_triquadratic, "2^27",
            "VOID: p-138.csv reports the triquadratic case space as {space_triquadratic}, \
             not 2^27"
        );
        assert_eq!(
            tractable_tricubic, "false",
            "VOID: p-138.csv reports the tricubic case space as tractable, which contradicts \
             the sentence this row carries"
        );
        assert_eq!(
            tractable_triquadratic, "true",
            "VOID: p-138.csv reports the triquadratic case space as intractable, which \
             contradicts the sentence this row carries"
        );
        assert_eq!(
            fixture.rows, FIXTURE_ROWS,
            "VOID: the committed fixture holds {} rows rather than {FIXTURE_ROWS}, so \
             `hashes_moved` is counted over a fixture this row does not know",
            fixture.rows
        );
        assert_eq!(
            fixture.checked, HASHES_CHECKED,
            "VOID: {} fixture rows were recomputed rather than {HASHES_CHECKED}, so \
             `hashes_moved = {}` is a zero over the wrong set (M-44)",
            fixture.checked, fixture.moved
        );
        assert!(
            ident.holds && ident.terms == 12 && ident.degree == 4 && ident.pencil_matches == 3,
            "VOID: P-127's identity does not re-measure here (holds={} terms={} degree={} \
             pencil_matches={}), so C3's claim about what a tricubic would invalidate rests \
             on algebra this run cannot reproduce",
            ident.holds,
            ident.terms,
            ident.degree,
            ident.pencil_matches
        );
        assert!(
            empty.corners.iter().all(|v| !is_inside(*v)) && empty.index == 0,
            "VOID: the empty-case block's corners are not all strictly outside (index {}), \
             so it is not a counterexample about case 0",
            empty.index
        );
        assert!(
            empty.centre < 0.0 && empty.triangles == 0,
            "VOID: the tricubic reads {:.6} at the centre of an all-outside cell for which \
             the shipped table emits {} triangles, so C3's two-empty-case half is unearned",
            empty.centre,
            empty.triangles
        );
        for f in &fields {
            // **The three arms are NOT expected to produce the same
            // triangulation, and demanding it contradicted the row's own
            // purpose.** A higher-order reconstruction moves the zero set, so
            // the tricubic arm's vertex and triangle counts differ from the
            // trilinear's by construction — that difference IS the effect
            // being measured. Measured on `sphere` at the first run: the
            // counts differ, and the first draft's control read that as a
            // broken comparison.
            //
            // What must be shared for the exponents to be comparable is the
            // GRID, which is structural — every arm extracts on the same
            // `shape`, `origin` and `h` at each rung, in one loop — and that
            // every arm produced a non-empty mesh at every rung, so no
            // exponent is fitted through an empty population. `topology_identical`
            // stays as a measured column.
            assert!(
                f.arms
                    .iter()
                    .all(|arm| arm.triangles > 0 && !arm.recon.is_empty()),
                "VOID: {}: an arm produced an empty mesh at the finest rung or an empty \
                 reconstruction ladder, so its fitted exponent is taken over a population that \
                 does not exist",
                f.name
            );
            for arm in &f.arms {
                assert_eq!(
                    arm.clamped,
                    0,
                    "VOID: {}: {} clamped the 4x4x4 stencil {} times, so some query fell \
                     outside the region where it fits and the reconstruction measured is \
                     not the one the header describes",
                    f.name,
                    arm.arm.filter(),
                    arm.clamped
                );
                assert!(
                    arm.finest.crossings > 0,
                    "VOID: {}: {} found no sign-changing grid edge at the finest rung, so \
                     `recon_zero_error` is a supremum over an empty set (M-44)",
                    f.name,
                    arm.arm.filter()
                );
                assert!(
                    arm.recon.iter().all(|e| e.is_finite()),
                    "VOID: {}: {} produced a non-finite reconstruction error on the ladder: \
                     {:?}",
                    f.name,
                    arm.arm.filter(),
                    arm.recon
                );
                if f.exact {
                    assert!(
                        arm.covered,
                        "VOID: {}: {} reported no accuracy coverage at some rung, so one \
                         point of its fit is not a measurement of anything",
                        f.name,
                        arm.arm.filter()
                    );
                    assert!(
                        arm.hausdorff.iter().all(|e| e.is_finite() && *e > 0.0),
                        "VOID: {}: {} produced a zero or non-finite symmetric Hausdorff on \
                         the ladder, so `log e` is undefined and the fit is fitted to \
                         nothing: {:?}",
                        f.name,
                        arm.arm.filter(),
                        arm.hausdorff
                    );
                }
            }
            // The oracle has no reconstruction, so whatever it still shows is
            // the mesh's error and nothing else. If its own crossings are not on
            // the field's zero set, that argument has no control behind it.
            let oracle_slot = Arm::ALL
                .iter()
                .position(|a| *a == Arm::Oracle)
                .expect("the oracle arm is in the roster");
            let oracle_worst = f.arms[oracle_slot].finest.worst;
            let finest_h = *f
                .cell_sizes
                .last()
                .expect("the ladder has at least one rung");
            assert!(
                oracle_worst <= finest_h * 1e-6,
                "VOID: {}: the exact-oracle arm's own zero-set error is {oracle_worst:e} \
                 against a bar of {:e}, so its {BISECTIONS} bisections are not resolving \
                 the field's zero and the chord-error argument has no control behind it",
                f.name,
                finest_h * 1e-6
            );
        }
        for &i in &smooth_exact {
            let measured = hausdorff_exponent[i][trilinear_slot];
            assert!(
                (measured - TRILINEAR_EXPONENT).abs() <= TRILINEAR_EXPONENT_BAND,
                "VOID: {}: the trilinear arm's fitted exponent is {measured} here, outside \
                 {TRILINEAR_EXPONENT} +/- {TRILINEAR_EXPONENT_BAND}, so this harness does \
                 not reproduce the order-2 baseline p-155.csv committed \
                 ({p155_sphere} on sphere, {p155_torus} on torus) and every exponent above \
                 is against an unmoored setup",
                fields[i].name
            );
            let gap = recon_exponent[i][tricubic_slot] - recon_exponent[i][trilinear_slot];
            assert!(
                gap >= RECON_EXPONENT_GAP,
                "VOID: {}: the tricubic's reconstruction exponent exceeds the trilinear's by \
                 only {gap}, under {RECON_EXPONENT_GAP}, so `higher-order reconstruction` is \
                 a name rather than a property here and C1 measures nothing",
                fields[i].name
            );
        }

        // ── the rows ────────────────────────────────────────────────────────
        let ladder = resolution_series();
        for (i, f) in fields.iter().enumerate() {
            let trilinear_per_cell = per_cell(f, trilinear_slot);
            let trilinear_total = f.arms[trilinear_slot].total_ms[0];
            for (slot, arm) in f.arms.iter().enumerate() {
                let evals_per_cell = per_cell(f, slot);
                let eval_share = arm.eval_ms[0] / arm.total_ms[0];
                let haus = finest(&arm.hausdorff);
                run.record(&[
                    // ── the registration's columns, in its order ────────────
                    ("filter", arm.arm.filter().to_string()),
                    ("approximation_order", arm.arm.order().to_string()),
                    ("field", f.name.to_string()),
                    ("resolution_series", ladder.clone()),
                    ("fitted_exponent", plain(hausdorff_exponent[i][slot])),
                    ("hausdorff", num(haus)),
                    ("samples_per_cell", plain(evals_per_cell)),
                    ("eval_ms", plain(arm.eval_ms[0])),
                    ("eval_share", plain(eval_share)),
                    ("cases_invalidated", audit.invalidated.to_string()),
                    ("hashes_moved", fixture.moved.to_string()),
                    ("c1_holds", c1.to_string()),
                    ("c2_holds", c2.to_string()),
                    ("c3_holds", c3.to_string()),
                    // ── extras (M-273) ─────────────────────────────────────
                    ("bisections", BISECTIONS.to_string()),
                    ("bound", f.bound.to_string()),
                    ("c1_hits", c1_hits.to_string()),
                    ("c1_population", c1_population.to_string()),
                    ("cases_centroid_bearing", audit.centroid_bearing.to_string()),
                    (
                        "cases_decider_dependent",
                        audit.decider_dependent.to_string(),
                    ),
                    ("cases_empty", audit.empty.to_string()),
                    ("cases_total", audit.total.to_string()),
                    ("cases_tunnel_capable", audit.tunnel_capable.to_string()),
                    ("cases_with_geometry", audit.with_geometry.to_string()),
                    ("cells_at_timing", cells.to_string()),
                    ("crossings_dropped", arm.finest.dropped.to_string()),
                    ("crossings_finest", arm.finest.crossings.to_string()),
                    ("empty_case_centre_value", plain(empty.centre)),
                    ("empty_case_index", empty.index.to_string()),
                    (
                        "empty_case_interior_component",
                        (empty.centre < 0.0).to_string(),
                    ),
                    ("empty_case_triangles", empty.triangles.to_string()),
                    ("eval_ms_max", plain(arm.eval_ms[2])),
                    ("eval_ms_min", plain(arm.eval_ms[1])),
                    ("fixture_rows", fixture.rows.to_string()),
                    ("gradient_calls", arm.gradient_calls.to_string()),
                    ("halo_layers", HALO.to_string()),
                    ("hashes_checked", fixture.checked.to_string()),
                    ("hausdorff_coverage", arm.covered.to_string()),
                    ("hausdorff_series", series(&arm.hausdorff)),
                    ("hausdorff_valid", f.exact.to_string()),
                    ("is_control", arm.arm.is_control().to_string()),
                    ("m150_path_ms", plain(M150_PATH_MS)),
                    ("m152_eval_ms", plain(M152_EVAL_MS)),
                    ("m152_upload_ms", plain(M152_UPLOAD_MS)),
                    ("p127_identity_degree", ident.degree.to_string()),
                    ("p127_identity_holds", ident.holds.to_string()),
                    ("p127_identity_multi_affine", ident.multi_affine.to_string()),
                    ("p127_identity_terms", ident.terms.to_string()),
                    ("p127_identity_variables", ident.variables.to_string()),
                    ("p127_pencil_matches", ident.pencil_matches.to_string()),
                    ("p138_case_space_tricubic", space_tricubic.clone()),
                    ("p138_case_space_trilinear", space_trilinear.clone()),
                    ("p138_case_space_triquadratic", space_triquadratic.clone()),
                    ("p138_mv_tricubic", mv_tricubic.clone()),
                    ("p138_mv_trilinear", mv_trilinear.clone()),
                    ("p138_mv_triquadratic", mv_triquadratic.clone()),
                    ("p138_net_sign_bits_tricubic", bits_tricubic.clone()),
                    ("p138_signs_are_case_index_tricubic", signs_tricubic.clone()),
                    ("p138_tractable_tricubic", tractable_tricubic.clone()),
                    (
                        "p138_tractable_triquadratic",
                        tractable_triquadratic.clone(),
                    ),
                    ("p155_fitted_exponent_sphere", p155_sphere.clone()),
                    ("p155_fitted_exponent_torus", p155_torus.clone()),
                    ("p155_measured_order", p155_measured.clone()),
                    ("p155_strang_fix_order", p155_order.clone()),
                    ("recon_error_series", series(&arm.recon)),
                    ("recon_fitted_exponent", plain(recon_exponent[i][slot])),
                    ("recon_mean_finest", num(arm.finest.mean)),
                    ("recon_zero_error", num(arm.finest.worst)),
                    (
                        "samples_per_cell_ratio",
                        plain(evals_per_cell / trilinear_per_cell),
                    ),
                    ("share_delta_ms_at_129", plain(share_delta_ms)),
                    ("share_delta_pct_at_129", plain(share_delta_pct)),
                    ("share_growth_at_129", plain(share_growth)),
                    ("smooth", f.smooth.to_string()),
                    ("stencil_clamped", arm.clamped.to_string()),
                    ("stencil_samples", arm.arm.stencil().to_string()),
                    ("timing_repeats", REPEATS.to_string()),
                    ("timing_resolution", TIMING_SAMPLES.to_string()),
                    ("timing_scatter", plain(arm.total_ms[2] / arm.total_ms[1])),
                    ("topology_identical", f.topology_identical.to_string()),
                    ("total_ms", plain(arm.total_ms[0])),
                    ("total_ms_max", plain(arm.total_ms[2])),
                    ("total_ms_min", plain(arm.total_ms[1])),
                    ("total_ratio", plain(arm.total_ms[0] / trilinear_total)),
                    ("triangles_finest", arm.triangles.to_string()),
                    ("vertices_finest", arm.vertices.to_string()),
                ]);
            }
        }
    });
}
