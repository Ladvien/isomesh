//! **P-162 — `0.2571 dB` is a prediction, not a hope: BCC sampling against the
//! cubic grid at matched point density.**
//!
//! Ticket: R-162. Pre-registered before this harness existed.
//!
//! ```bash
//! cargo bench --bench experiment_p162
//! ```
//!
//! Writes `docs/experiments/p-162.csv`.
//!
//! # What was missing
//!
//! **This crate has sampled on `Z³` for its entire life and has never measured
//! the alternative.** Every extractor's `extract` takes
//! `(sdf, shape: &impl Shape3, origin, cell_size, out)` — seven of them, the
//! signature byte-identical in each (`marching_cubes/mod.rs:193`,
//! `marching_tetrahedra.rs:125`, `surface_nets.rs:193`, `dual_contouring.rs:277`,
//! `manifold_dual_contouring.rs:317`, `greedy_quads.rs:155`,
//! `subgrid/extract.rs:564`) — and `Shape3`'s strides are `[1, sx, sx·sy]`
//! (`shape.rs:11-22`). That is a cubic lattice with no seam at which another one
//! could enter. Grepping `lattice` over `crates/isomesh/src` returns the *chunk*
//! lattice (`chunk.rs:40`), the noise generator's own integer lattice
//! (`fields/noise.rs:71`, `fields/mod.rs:1142`) and the sample lattice's symmetry
//! group (`equivariant.rs:9`, `dual_contouring/solve.rs:95`). Not one line offers
//! a second sampling geometry, and no registration `P-8`–`P-161` asks for one.
//!
//! What is new here is not the idea. Barnes & Sloane (`10.1137/0604005`) proved
//! `A₃*` — the body-centred cubic lattice — optimal **among three-dimensional
//! lattices** in 1983, and this crate has been able to read that for as long as it
//! has existed. What is new is that the claim is **arithmetically closed**:
//! `G(Z³) = 1/12` and `G(A₃*) = 19/(192·∛2)` are exact, so
//!
//! ```text
//!   10·log₁₀(G(Z³) / G(A₃*)) = 0.257097 dB,     G(Z³)/G(A₃*) − 1 = 5.7481%
//! ```
//!
//! with nothing fitted and nothing measured. `common::lattice::Lattice::g`
//! asserts both decimals against those closed forms on **every call**
//! (`benches/common/lattice.rs:239-250`, residuals `2.17e-10` and `3.82e-10`), so
//! a transposed digit cannot become a published prediction. A number that precise
//! is either met or it is wrong, and until this row ran nothing in this repository
//! knew which — nor even whether `G`, which is a *quantisation* figure of merit,
//! predicts anything at all about a *geometric* error. C2 exists because the
//! honest answer to that was "nobody has checked".
//!
//! The machinery is `common::lattice`, which this ticket owns and which R-163 and
//! R-164 consume unchanged. Everything below drives it and the public `isomesh`
//! API and writes nothing to `crates/isomesh/src/**`.
//!
//! # Arms
//!
//! Two per field, eight fields, sixteen rows. `(field, lattice)` is the primary
//! key, and every per-field number is recorded on both of its rows so that P-164
//! can read this CSV as its baseline without joining anything.
//!
//! | arm | what it varies | is_control |
//! |---|---|---|
//! | `Z3` — cubic sites, trilinear reconstruction | nothing: the lattice and the filter the crate ships | **yes** |
//! | `A3*` — BCC sites, linear box-spline reconstruction | the sampling lattice, and the reconstruction filter that belongs to it | no |
//!
//! The two arms differ in the lattice and its filter **and in nothing else**: the
//! same eight reference fields, the same box, matched site counts, the same
//! contouring grid, the same refinement, the same probe count and the same
//! `zero_set_hausdorff`. `common::lattice::trilinear_reconstruct` looks its eight
//! corners up by exactly the binary search `bcc_reconstruct` looks its four up by
//! (`lattice.rs:1026-1031`), so even the lookup cost is a shared constant rather
//! than an asymmetry.
//!
//! # Method, and the four choices that decide what the numbers mean
//!
//! **1. Matched point density, in one direction only.** The cubic lattice anchored
//! on the box centre can realise only an *odd* number of sites per axis, so its
//! attainable totals are `15³`, `17³`, … — gaps of 30% and more — while BCC
//! interleaves two cubic sub-lattices and is far finer grained
//! (`lattice.rs:339-362`). So the cubic arm is built first, at
//! `TARGET_POINTS = 49³ = 117,649`, and the BCC arm is then asked for *the count
//! the cubic arm realised*. Asking both for the same round number matches neither.
//! `samples` is recorded per row from `LatticeGrid::sites.len()`, never from
//! `TARGET_POINTS`, and `density_mismatch` is the residual gap — the vacuity
//! control is exactly this number.
//!
//! `49³` rather than `33³`: `noise_cavity`'s features are about `1/3.45 ≈ 0.29`
//! across (`fields/mod.rs:1152-1153`) and its domain spans `4`, so `49` sites per
//! axis is `3.5` samples per feature. At `25³` it would be `1.7` — below Nyquist,
//! and a comparison of two aliased reconstructions measures aliasing.
//!
//! **2. Where the reconstructed zero set comes from.** Each arm samples its field
//! at every site, then *contours its own reconstruction* on a shared uniform
//! `49³` grid over the box inset by `2.5` lattice scales: every sign change on an
//! axis-aligned grid edge is bisected `14` times against the reconstruction, and
//! the crossing is kept. The inset is not cosmetic — `bcc_reconstruct` panics
//! unless the box spline's weights sum to `1`, which is the exact condition for
//! "no needed site was clipped away" (`lattice.rs:957-964`), and the support
//! radius is `2` lattice units, `2·scale/∛4 = 1.26·scale` in world distance.
//! `2.5·scale` covers that and the trilinear cell's `√3·scale` together.
//!
//! **Why every crossing is guaranteed to project onto the true zero set.** Both
//! filters are **convex combinations**: `bcc_box_spline` returns
//! `(high − low).max(0.0)`, non-negative (`lattice.rs:939-947`), and is a
//! partition of unity, and the trilinear weights are products of `t`/`1−t`. So a
//! reconstruction value is bounded between the minimum and maximum of the sample
//! values in its stencil, and a *zero crossing of the reconstruction implies a
//! genuine sign change among sites within one stencil radius*. There is therefore
//! a true zero within `~2.5·scale` of every point handed to
//! `zero_set_hausdorff`, which is why its "a point of the reconstructed zero set
//! does not project" panic (`lattice.rs:1141-1146`) is unreachable here rather
//! than merely unlikely. A filter with negative lobes would not have this
//! property and this harness would not be sound with one.
//!
//! **3. The dB convention, and it is a choice.** `G` is a **second moment** — a
//! mean *squared* error — so the registered `0.2571 dB` lives in the power
//! convention, `10·log₁₀` of a ratio of squared errors. Hausdorff distance is a
//! *linear* distance, an amplitude. To be commensurable with the prediction it
//! must therefore enter as
//!
//! ```text
//!   measured_gain_db = 20·log₁₀(h_cubic / h_bcc) = 10·log₁₀((h_cubic / h_bcc)²)
//! ```
//!
//! — the two are the same number, and the second form is the one that makes the
//! commensurability obvious. `20·log₁₀` it is, on every row, consistently. The
//! sign convention is the module's: positive means *gain available by moving from
//! the cubic grid to BCC*, which is the direction
//! `Lattice::Cubic.gain_db_over(Lattice::Bcc)` reads in and the direction the
//! registration's `+0.2571` is written in (`lattice.rs:140-159`). A cubic row
//! carries `hausdorff_ratio = 1` and `measured_gain_db = 0` by construction; that
//! is the control saying so, not a missing measurement.
//!
//! **4. C2's escape clause is measured, not asserted.** The registration allows
//! C2 to hold if the deviation *is explained*, "because `G` predicts quantization
//! error and Hausdorff error is not that, and the gap between them is the
//! interesting part". Hausdorff is a **maximum**; `G` predicts a **mean square**.
//! So the same crossing set is also reduced to an RMS, using the linearised
//! distance `|f(p)| / ‖∇f(p)‖` at each crossing — a first-order distance that
//! needs no iteration and whose common second-order bias cancels in the *ratio*
//! between two arms measured on the identical point geometry. `rms_gain_db` is
//! that ratio in the same `20·log₁₀` convention. If the max-form misses the
//! prediction and the mean-square form lands on it, the deviation is not merely
//! explained in prose — it is *located*, and the location is "max versus mean".
//! That is the only reading of "explained" this harness will accept, and it is
//! stated here before the run.
//!
//! Determinism: no RNG of this harness's own. The only stochastic element is
//! `zero_set_hausdorff`'s probe stream, a `SplitMix64` on the fixed seed
//! `0x1362_A3B5_D1E7_9F11` (`lattice.rs:1096-1098`), seeded into the bounding box
//! of the *point set*, so the two arms' probes differ by exactly the difference
//! between their two bounding boxes and by nothing else. `PROBES = 5000`,
//! recorded, because a Hausdorff distance quoted without its probe count is not
//! reproducible.
//!
//! **What `extraction_ms` is and is not.** SHARE says C3 is a complexity cost and
//! not a runtime one, and this column is the reason that sentence is in the
//! registration. This harness's per-evaluation cost is dominated by
//! `LatticeGrid::find`, an `O(log sites)` binary search per stencil tap
//! (`lattice.rs:298-315`), which a shipped extractor with a linear index space
//! would not pay at all — and the BCC arm pays it four times per evaluation
//! against the cubic arm's eight, so BCC is expected to be the *faster* arm here
//! for a reason that has nothing to do with either lattice. `extraction_ms` is
//! recorded because the registration names it; it is not evidence about either
//! lattice's runtime. Five timed repeats, median as the headline, min and max as
//! extras, warmed up once first — M-280 measured this host's `amd-pstate-epp`
//! governor swinging the same binary `1.45×` between runs, so a single wall clock
//! here would be a number about the governor.
//!
//! # SHARE, recomputed before the numbers
//!
//! *"C1 moves the whole sampling stage; C3 is a complexity cost, not a runtime
//! one."*
//!
//! **C1 moves the sampling stage, and this row does not move it.** A BCC arm in
//! the shipped crate needs a site enumerator, a BCC reconstruction filter, a
//! tetrahedral extractor over the BCC Delaunay complex, and every consumer of
//! `Shape3`'s `[1, sx, sx·sy]` strides to stop assuming a cubic index space.
//! Nothing here proposes any of it: `crates/isomesh/src/**` is untouched, no
//! reference field is added, and no golden hash can move. A positive C1 is a
//! landing **ticket**, registered in advance with its golden-hash and doc-facts
//! ripple priced into it — a landing that happens quietly inside a measurement
//! commit is `V-45`'s failure mode.
//!
//! **C3 is a complexity cost, and here is the arithmetic.** The BCC case table is
//! `16` entries against the cubic `256`, generated rather than transcribed and
//! calibrated entry-by-entry against `isomesh::marching_cubes::table::CASES`
//! (`lattice.rs:598-607`). But a case count alone flatters BCC, so the cost is
//! stated in full: BCC's natural cell is its Delaunay tetrahedron, of which there
//! are **6 per lattice site** against the cube's **1** (`lattice.rs:609-625`), and
//! each cell reads **4** corners against the cube's **8**. So per lattice site:
//! `6 × 4 = 24` corner reads against `1 × 8 = 8`, three times as many, against a
//! table sixteen times smaller. That is the trade, and it is a complexity trade —
//! not the wall clock in `extraction_ms`.
//!
//! # Vacuity controls
//!
//! All five run before the first `run.record`, and every panic message starts
//! `VOID: `.
//!
//! - **Matched point density.** `density_mismatch = |samples_bcc − samples_cubic|
//!   / samples_cubic` must be at most `5%` on every field, and both counts are
//!   recorded. This is the registration's own vacuity control, verbatim: without
//!   it "the comparison is a resolution change wearing a lattice's name".
//!   Columns: `samples`, `density_mismatch`.
//! - **The prediction under test is the registered one.**
//!   `Lattice::Cubic.gain_db_over(Lattice::Bcc)` must equal the registered
//!   `0.2571 dB` to `5e-5`, or C2 would be scored against a number nobody
//!   registered. Column: `predicted_gain_db`.
//! - **Both arms measured a real surface.** Each arm's crossing count must reach
//!   `64` and its Hausdorff and RMS must be strictly positive and finite: a
//!   maximum over a handful of points is not a surface, and a ratio of two zeros
//!   is not a gain (M-44). Columns: `points`, `hausdorff`, `rms_error`.
//! - **The error is not this harness's own refinement noise.** Each arm's
//!   Hausdorff must exceed the bisection residual — one eval cell over `2¹⁴` — by
//!   at least `100×`, or `hausdorff_ratio` would be reporting the bisection
//!   rather than the lattice. Columns: `hausdorff`, `eval_samples`,
//!   `refine_steps`.
//! - **The case-table enumeration understands the cube.** The generated cubic
//!   table's total triangle count must equal the shipped
//!   `isomesh::marching_cubes::table::CASES` total, entry-by-entry agreement
//!   being asserted inside `case_table` itself. Two independent enumerations
//!   reaching the same 256 numbers is what licenses believing the 16 beside them.
//!   Columns: `case_table_size`, `total_triangles`.

mod common;

use std::time::Instant;

use common::lattice::{
    BCC_BOX_SPLINE_STENCIL, CaseTable, Lattice, LatticeGrid, TRILINEAR_STENCIL, bcc_reconstruct,
    case_table, lattice_grid, shipped_cubic_triangle_counts, trilinear_reconstruct,
    zero_set_hausdorff,
};
use isomesh::Sdf;
use isomesh::fields::ReferenceField;

// ─── the configuration, all of it derived in the header ─────────────────────

/// Sites the **cubic** arm is asked for: `49³`.
///
/// The cubic lattice anchored on the box centre realises only odd counts per
/// axis, and `49³` is one it can hit exactly. `49` also puts `noise_cavity` at
/// `3.5` samples per feature rather than `25³`'s aliased `1.7`.
const TARGET_POINTS: usize = 117_649;

/// Samples per axis of the shared contouring grid.
///
/// Over the inset box this is a slightly *finer* spacing than the sample
/// lattice's, so the crossing set resolves the reconstruction rather than
/// band-limiting it.
const EVAL_SAMPLES: usize = 49;

/// Bisection steps used to place a crossing on the reconstruction's zero set.
///
/// `2⁻¹⁴` of an eval cell. The residual must be far under the error being
/// measured or the ratio reports the bisection; the fourth vacuity control is
/// that inequality.
const REFINE_STEPS: u32 = 14;

/// Probe seeds for `zero_set_hausdorff`'s truth-to-reconstruction direction.
const PROBES: usize = 5_000;

/// Timed repeats of each arm's extraction. Five, plus one warm-up.
const TIMED_REPEATS: usize = 5;

/// Box inset, in lattice scales.
///
/// The box spline's support radius is `2` lattice units, `2·scale/∛4 =
/// 1.26·scale` in world distance; the trilinear cell reaches `√3·scale`. This
/// covers both.
const INSET_SCALES: f64 = 2.5;

/// Largest site-count gap between the two arms this comparison will accept.
const DENSITY_TOLERANCE: f64 = 0.05;

/// Fewest crossings an arm may report and still be describing a surface.
const MIN_CROSSINGS: usize = 64;

/// How far above the bisection residual a Hausdorff must sit.
const RESIDUAL_HEADROOM: f64 = 100.0;

/// The dB figure the registration predicts, quoted from it.
const REGISTERED_GAIN_DB: f64 = 0.2571;

/// Agreement required between the registration's rounded dB and the module's
/// computed one.
const GAIN_DB_TOLERANCE: f64 = 5e-5;

/// C2's neighbourhood: "within a factor of 2".
const NEIGHBOURHOOD_FACTOR: f64 = 2.0;

/// C1's threshold: "at least five of eight fields".
const FIELDS_REQUIRED: usize = 5;

/// The amplitude dB constant, `20`.
///
/// Not `10`, and the choice is argued in the header: `G` is a mean *squared*
/// error, Hausdorff is a linear distance, so
/// `20·log₁₀(h₁/h₂) = 10·log₁₀((h₁/h₂)²)` is the commensurable form.
const AMPLITUDE_DB: f64 = 20.0;

/// Delaunay tetrahedra per BCC lattice site: **6**.
///
/// Not a mechanism reimplemented here but a constant cited from the module that
/// owns the decomposition: the representative tetrahedron has volume `2/3`
/// against a volume-per-site of `4` in `bcc_box_spline`'s integer coordinates,
/// so `4 / (2/3) = 6` (`benches/common/lattice.rs:609-625`).
const BCC_CELLS_PER_SITE: usize = 6;

/// Cubic cells per cubic lattice site: **1**.
const CUBIC_CELLS_PER_SITE: usize = 1;

// ─── the measurement ────────────────────────────────────────────────────────

/// One arm: one lattice, one field.
#[derive(Debug)]
struct Arm {
    /// Which lattice this arm sampled on.
    lattice: Lattice,
    /// Sites the lattice actually realised in the box — the matched-density
    /// number, read from `LatticeGrid::sites`, never from [`TARGET_POINTS`].
    samples: usize,
    /// The factor the unit-volume generator rows were multiplied by.
    scale: f64,
    /// Crossings of this arm's reconstruction on the shared contouring grid.
    points: usize,
    /// Symmetric Hausdorff between those crossings and the field's true zero
    /// set.
    hausdorff: f64,
    /// RMS of the linearised distance `|f| / ‖∇f‖` over the crossings — the
    /// mean-square form `G` actually predicts.
    rms: f64,
    /// Mean of the same linearised distance.
    mean: f64,
    /// Largest linearised distance over the crossings.
    worst: f64,
    /// Median of [`TIMED_REPEATS`] extractions, in milliseconds.
    ms_median: f64,
    /// Fastest of them.
    ms_min: f64,
    /// Slowest of them.
    ms_max: f64,
}

/// One field: both arms, and the comparison between them.
#[derive(Debug)]
struct FieldRow {
    /// The reference field's name.
    field: &'static str,
    /// World distance the contouring box was inset by.
    inset: f64,
    /// Spacing of the contouring grid, needed by the refinement-noise control.
    eval_cell: f64,
    /// `|samples_bcc − samples_cubic| / samples_cubic`.
    mismatch: f64,
    /// `h_cubic / h_bcc`. Above 1 means BCC is the better lattice here.
    ratio: f64,
    /// `20·log₁₀(ratio)`.
    gain_db: f64,
    /// The same, on the RMS linearised distance instead of the Hausdorff max.
    rms_gain_db: f64,
    /// C1's per-field question: did BCC improve the Hausdorff at all.
    improved: bool,
    /// C2's per-field question, on the Hausdorff.
    prediction_holds: bool,
    /// C2's escape clause, per field: the mean-square form of the same
    /// measurement.
    rms_prediction_holds: bool,
    /// The cubic control arm.
    cubic: Arm,
    /// The BCC arm.
    bcc: Arm,
}

/// The reconstruction filter that belongs to a lattice.
///
/// One filter per lattice and no configuration: the trilinear is what the crate
/// ships on `Z³`, the linear box spline is the simplex element of BCC's own
/// Delaunay complex, and both have approximation order 2 — so the arms differ in
/// the lattice and its natural filter, not in the filter's order.
fn reconstruct(grid: &LatticeGrid, values: &[f64], p: [f64; 3]) -> f64 {
    match grid.lattice {
        Lattice::Cubic => trilinear_reconstruct(grid, values, p),
        Lattice::Bcc => bcc_reconstruct(grid, values, p),
        Lattice::Fcc => unreachable!(
            "P-162 compares Z3 against A3* only; D3 is P-163's row and is built by its own harness"
        ),
    }
}

/// Sample the field at every site of the grid, in site order.
fn sample_sites<F>(field: &F, grid: &LatticeGrid) -> Vec<f64>
where
    F: Sdf<Scalar = f64>,
{
    grid.sites.iter().map(|s| field.sample(*s)).collect()
}

/// Place a crossing on the reconstruction's zero set by bisection.
///
/// `inside_at_a` is the sign class of the reconstruction at `a`; the invariant
/// maintained is that `lo` keeps that class and `hi` does not.
fn refine(
    grid: &LatticeGrid,
    values: &[f64],
    a: [f64; 3],
    b: [f64; 3],
    inside_at_a: bool,
) -> [f64; 3] {
    let mut lo = a;
    let mut hi = b;
    for _ in 0..REFINE_STEPS {
        let m = [
            f64::midpoint(lo[0], hi[0]),
            f64::midpoint(lo[1], hi[1]),
            f64::midpoint(lo[2], hi[2]),
        ];
        if (reconstruct(grid, values, m) < 0.0) == inside_at_a {
            lo = m;
        } else {
            hi = m;
        }
    }
    [
        f64::midpoint(lo[0], hi[0]),
        f64::midpoint(lo[1], hi[1]),
        f64::midpoint(lo[2], hi[2]),
    ]
}

/// Contour one arm's reconstruction on the shared `EVAL_SAMPLES³` grid.
///
/// The grid, the sign rule and the loop order are identical for both arms, so
/// the two crossing sets differ only in the reconstruction that produced them.
/// Axis outermost then `k`, `j`, `i` — a fixed order, so the point list is
/// deterministic and the `Instant` repeats are comparable.
fn zero_crossings(
    grid: &LatticeGrid,
    values: &[f64],
    elo: [f64; 3],
    ehi: [f64; 3],
) -> Vec<[f64; 3]> {
    let n = EVAL_SAMPLES;
    let last = (n - 1) as f64;
    let step = [
        (ehi[0] - elo[0]) / last,
        (ehi[1] - elo[1]) / last,
        (ehi[2] - elo[2]) / last,
    ];
    let at = |c: [usize; 3]| {
        [
            elo[0] + step[0] * c[0] as f64,
            elo[1] + step[1] * c[1] as f64,
            elo[2] + step[2] * c[2] as f64,
        ]
    };
    let index = |c: [usize; 3]| c[0] + n * (c[1] + n * c[2]);

    let mut sampled = vec![0.0f64; n * n * n];
    for k in 0..n {
        for j in 0..n {
            for i in 0..n {
                sampled[index([i, j, k])] = reconstruct(grid, values, at([i, j, k]));
            }
        }
    }

    let mut out = Vec::new();
    for axis in 0..3usize {
        let mut delta = [0usize; 3];
        delta[axis] = 1;
        for k in 0..n - delta[2] {
            for j in 0..n - delta[1] {
                for i in 0..n - delta[0] {
                    let a = [i, j, k];
                    let b = [i + delta[0], j + delta[1], k + delta[2]];
                    let va = sampled[index(a)];
                    let vb = sampled[index(b)];
                    if (va < 0.0) == (vb < 0.0) {
                        continue;
                    }
                    out.push(refine(grid, values, at(a), at(b), va < 0.0));
                }
            }
        }
    }
    out
}

/// One extraction: sample the lattice, then contour the reconstruction.
///
/// This is the unit `extraction_ms` times, and it is deliberately the whole of
/// it — sampling the field is part of what a sampling-lattice change costs.
fn extract<F>(field: &F, grid: &LatticeGrid, elo: [f64; 3], ehi: [f64; 3]) -> Vec<[f64; 3]>
where
    F: Sdf<Scalar = f64>,
{
    let values = sample_sites(field, grid);
    zero_crossings(grid, &values, elo, ehi)
}

/// Linearised distance from `p` to the field's zero set: `|f(p)| / ‖∇f(p)‖`.
///
/// First order, and that is enough: it is used only inside a *ratio* between two
/// arms whose crossings sit at the same distance scale, where the common
/// second-order bias cancels. Iterating instead would put a second convergence
/// criterion inside an error measurement.
fn linear_residual<F>(field: &F, p: [f64; 3]) -> f64
where
    F: Sdf<Scalar = f64>,
{
    let value = field.sample(p);
    let g = field.gradient(p);
    let norm = (g[0] * g[0] + g[1] * g[1] + g[2] * g[2]).sqrt();
    assert!(
        norm > 0.0 && norm.is_finite(),
        "the field's gradient is {g:?} at {p:?}, so the linearised distance to its zero set \
         is undefined there"
    );
    value.abs() / norm
}

/// Measure one arm: extract, time, then measure the two error forms.
fn arm<F>(field: &F, grid: &LatticeGrid, elo: [f64; 3], ehi: [f64; 3]) -> Arm
where
    F: Sdf<Scalar = f64>,
{
    // The warm-up run, whose point set every number below is computed from — so
    // the errors and the timings describe the same extraction rather than two.
    let points = extract(field, grid, elo, ehi);

    let mut spans: Vec<u128> = Vec::with_capacity(TIMED_REPEATS);
    for _ in 0..TIMED_REPEATS {
        let started = Instant::now();
        let again = extract(field, grid, elo, ehi);
        spans.push(started.elapsed().as_nanos());
        assert_eq!(
            again.len(),
            points.len(),
            "the extraction is not deterministic on {}: {} crossings against {}",
            grid.lattice.name(),
            again.len(),
            points.len()
        );
    }
    spans.sort_unstable();

    let hausdorff = zero_set_hausdorff(field, &points, PROBES);

    let mut sum = 0.0f64;
    let mut sum_sq = 0.0f64;
    let mut worst = 0.0f64;
    for p in &points {
        let r = linear_residual(field, *p);
        sum += r;
        sum_sq += r * r;
        worst = worst.max(r);
    }
    let count = points.len() as f64;

    Arm {
        lattice: grid.lattice,
        samples: grid.sites.len(),
        scale: grid.scale,
        points: points.len(),
        hausdorff,
        rms: (sum_sq / count).sqrt(),
        mean: sum / count,
        worst,
        ms_median: spans[TIMED_REPEATS / 2] as f64 / 1e6,
        ms_min: spans[0] as f64 / 1e6,
        ms_max: spans[TIMED_REPEATS - 1] as f64 / 1e6,
    }
}

/// Is `measured` within [`NEIGHBOURHOOD_FACTOR`] of `predicted`, both ways.
fn in_neighbourhood(measured: f64, predicted: f64) -> bool {
    measured >= predicted / NEIGHBOURHOOD_FACTOR && measured <= predicted * NEIGHBOURHOOD_FACTOR
}

/// Measure both arms on one field and derive the comparison.
fn measure<F>(name: &'static str, field: &F) -> FieldRow
where
    F: ReferenceField + Sdf<Scalar = f64>,
{
    let (lo, hi) = field.domain();

    // Cubic first, then BCC anchored on the count the cubic arm realised. The
    // asymmetry is the module's documented protocol and the reason the vacuity
    // control can pass at all (`benches/common/lattice.rs:339-362`).
    let cubic = lattice_grid(Lattice::Cubic, lo, hi, TARGET_POINTS);
    let bcc = lattice_grid(Lattice::Bcc, lo, hi, cubic.sites.len());

    let inset = INSET_SCALES * cubic.scale.max(bcc.scale);
    let elo = [lo[0] + inset, lo[1] + inset, lo[2] + inset];
    let ehi = [hi[0] - inset, hi[1] - inset, hi[2] - inset];
    assert!(
        ehi[0] > elo[0] && ehi[1] > elo[1] && ehi[2] > elo[2],
        "{name}: an inset of {inset} leaves no interior in {lo:?}..{hi:?}"
    );

    let cubic_arm = arm(field, &cubic, elo, ehi);
    let bcc_arm = arm(field, &bcc, elo, ehi);

    let mismatch =
        (bcc_arm.samples as f64 - cubic_arm.samples as f64).abs() / cubic_arm.samples as f64;
    let ratio = cubic_arm.hausdorff / bcc_arm.hausdorff;
    let gain_db = AMPLITUDE_DB * ratio.log10();
    let rms_gain_db = AMPLITUDE_DB * (cubic_arm.rms / bcc_arm.rms).log10();
    let predicted = Lattice::Cubic.gain_db_over(Lattice::Bcc);

    FieldRow {
        field: name,
        inset,
        eval_cell: (ehi[0] - elo[0]) / (EVAL_SAMPLES - 1) as f64,
        mismatch,
        ratio,
        gain_db,
        rms_gain_db,
        improved: bcc_arm.hausdorff < cubic_arm.hausdorff,
        prediction_holds: in_neighbourhood(gain_db, predicted),
        rms_prediction_holds: in_neighbourhood(rms_gain_db, predicted),
        cubic: cubic_arm,
        bcc: bcc_arm,
    }
}

/// The complexity numbers one lattice's extractor would carry: cells per site,
/// reconstruction stencil, and corner reads per site.
///
/// Returned together because C3 is falsified by "a case table large enough to be
/// impractical", and a case count read without its cells-per-site multiplier is
/// the flattering half of the trade.
fn cost(lattice: Lattice, table: &CaseTable) -> (usize, usize, usize) {
    let cells = match lattice {
        Lattice::Cubic => CUBIC_CELLS_PER_SITE,
        Lattice::Bcc | Lattice::Fcc => BCC_CELLS_PER_SITE,
    };
    let stencil = match lattice {
        Lattice::Cubic => TRILINEAR_STENCIL,
        Lattice::Bcc | Lattice::Fcc => BCC_BOX_SPLINE_STENCIL,
    };
    (cells, stencil, cells * table.corners_per_cell)
}

fn main() {
    if !std::env::args().any(|arg| arg == "--bench") {
        return;
    }
    let prereg = isomesh::experiment!("P-162");

    common::experiment::run(prereg, |run| {
        let cubic_table = case_table(Lattice::Cubic);
        let bcc_table = case_table(Lattice::Bcc);
        let predicted = Lattice::Cubic.gain_db_over(Lattice::Bcc);

        println!(
            "predicted gain  G(Z3) = {:.9}  G(A3*) = {:.9}  ->  {predicted:.6} dB  \
             (MSE reduction {:.4}%)",
            Lattice::Cubic.g(),
            Lattice::Bcc.g(),
            // `1 - G(A3*)/G(Z3)`, which is the 5.7481% the registration quotes.
            // The reciprocal reading, `G(Z3)/G(A3*) - 1`, is 6.0986% and is a
            // different quantity: the excess of the worse lattice over the
            // better, not the reduction the better one delivers.
            100.0 * (1.0 - Lattice::Bcc.g() / Lattice::Cubic.g())
        );
        println!(
            "case tables     Z3: {} cases, {} distinct, max {} tris, {} tris total  |  \
             A3*: {} cases, {} distinct, max {} tris, {} tris total",
            cubic_table.cases,
            cubic_table.distinct_up_to_symmetry,
            cubic_table.max_triangles_per_case,
            cubic_table.total_triangles,
            bcc_table.cases,
            bcc_table.distinct_up_to_symmetry,
            bcc_table.max_triangles_per_case,
            bcc_table.total_triangles
        );
        println!();

        let mut rows: Vec<FieldRow> = Vec::new();
        isomesh::for_each_reference_field!(f64, |name, field| {
            // An inline block per field, not a closure: a `return` in here would
            // return from `main` and the run would stop at the first field
            // (M-253).
            let row = measure(name, &field);
            println!(
                "{:>15}  samples {:>7} / {:>7} ({:.3}% apart)  h {:.6e} / {:.6e}  \
                 ratio {:.5}  {:+.4} dB",
                row.field,
                row.cubic.samples,
                row.bcc.samples,
                row.mismatch * 100.0,
                row.cubic.hausdorff,
                row.bcc.hausdorff,
                row.ratio,
                row.gain_db
            );
            println!(
                "{:>15}  rms {:.6e} / {:.6e}  {:+.4} dB   points {:>6} / {:>6}   \
                 ms {:.3} / {:.3}",
                "",
                row.cubic.rms,
                row.bcc.rms,
                row.rms_gain_db,
                row.cubic.points,
                row.bcc.points,
                row.cubic.ms_median,
                row.bcc.ms_median
            );
            rows.push(row);
        });
        println!();

        // ── the vacuity controls, before any verdict is reported ─────────────

        // 1. The registration's own control: genuinely matched point density,
        //    reported as counts.
        for row in &rows {
            assert!(
                row.mismatch <= DENSITY_TOLERANCE,
                "VOID: {}: the arms hold {} cubic sites against {} BCC sites, a {:.3}% gap \
                 against the {:.1}% this comparison allows — at that gap the row is a \
                 resolution change wearing a lattice's name, which is the one thing P-162's \
                 vacuity control forbids",
                row.field,
                row.cubic.samples,
                row.bcc.samples,
                row.mismatch * 100.0,
                DENSITY_TOLERANCE * 100.0
            );
        }

        // 2. The number C2 is scored against is the number that was registered.
        assert!(
            (predicted - REGISTERED_GAIN_DB).abs() < GAIN_DB_TOLERANCE,
            "VOID: the module computes {predicted:.9} dB from G(Z3)/G(A3*) while the \
             registration predicts {REGISTERED_GAIN_DB} dB — C2 would be scored against a \
             prediction nobody registered"
        );

        // 3. Both arms measured a real surface, so no ratio here is a ratio of
        //    two zeros (M-44).
        for row in &rows {
            for side in [&row.cubic, &row.bcc] {
                assert!(
                    side.points >= MIN_CROSSINGS,
                    "VOID: {} on {}: only {} crossings, under the {MIN_CROSSINGS} this \
                     harness will call a surface — a Hausdorff maximum over a handful of \
                     points is not an error measurement",
                    row.field,
                    side.lattice.name(),
                    side.points
                );
                assert!(
                    side.hausdorff > 0.0 && side.hausdorff.is_finite(),
                    "VOID: {} on {}: Hausdorff {} — a zero or non-finite error makes \
                     hausdorff_ratio meaningless and C1's comparison a comparison of two \
                     zeros",
                    row.field,
                    side.lattice.name(),
                    side.hausdorff
                );
                assert!(
                    side.rms > 0.0 && side.rms.is_finite(),
                    "VOID: {} on {}: RMS linearised distance {} — C2's escape clause would \
                     be evaluated on a zero",
                    row.field,
                    side.lattice.name(),
                    side.rms
                );
            }
        }

        // 4. The errors are the lattice's, not this harness's own bisection
        //    residual.
        for row in &rows {
            let residual = row.eval_cell * 2f64.powi(-(REFINE_STEPS as i32));
            for side in [&row.cubic, &row.bcc] {
                assert!(
                    side.hausdorff > RESIDUAL_HEADROOM * residual,
                    "VOID: {} on {}: Hausdorff {:.6e} against a bisection residual of \
                     {residual:.6e} (one {:.6e} eval cell over 2^{REFINE_STEPS}) — under \
                     {RESIDUAL_HEADROOM}x headroom the ratio reports this harness's \
                     refinement, not the sampling lattice",
                    row.field,
                    side.lattice.name(),
                    side.hausdorff,
                    row.eval_cell
                );
            }
        }

        // 5. The enumeration C3's BCC number comes from understands the cube.
        //    `case_table` asserts all 256 entries individually; this is the
        //    aggregate, so a reader of the CSV can see it was asked.
        let shipped: usize = shipped_cubic_triangle_counts()
            .iter()
            .map(|c| *c as usize)
            .sum();
        assert_eq!(
            cubic_table.total_triangles, shipped,
            "VOID: the generated cubic table emits {} triangles against the shipped \
             isomesh::marching_cubes::table::CASES total of {shipped}, so the enumeration \
             that produced the BCC 16 does not reproduce the cube's 256 and C3's number is \
             uncalibrated",
            cubic_table.total_triangles
        );

        // ── the verdicts, global, with their arithmetic ──────────────────────

        let improved = rows.iter().filter(|r| r.improved).count();
        let near = rows.iter().filter(|r| r.prediction_holds).count();
        let near_rms = rows.iter().filter(|r| r.rms_prediction_holds).count();

        // C1: at matched sample count, BCC improves symmetric Hausdorff on at
        // least five of eight fields.
        let c1 = improved >= FIELDS_REQUIRED;

        // C2: the improvement is within a factor of two of the predicted
        // 0.2571 dB — or the deviation is EXPLAINED, and the only reading of
        // "explained" this harness accepts is the one it can measure: `G`
        // predicts a mean square and Hausdorff is a maximum, so the mean-square
        // form of the identical measurement landing in the neighbourhood
        // locates the deviation at max-versus-mean rather than at `G`.
        let c2 = near >= FIELDS_REQUIRED || near_rms >= FIELDS_REQUIRED;

        // C3: the cost is stated, and the BCC table is not "large enough to be
        // impractical" — the threshold being the cubic table this crate already
        // ships and calls practical.
        let c3 = bcc_table.cases <= cubic_table.cases;

        println!(
            "C1  {improved} of {} fields improved (needs {FIELDS_REQUIRED}) -> {c1}",
            rows.len()
        );
        println!(
            "C2  hausdorff within 2x of {predicted:.6} dB on {near} of {}; rms within 2x on \
             {near_rms} of {} -> {c2}",
            rows.len(),
            rows.len()
        );
        println!(
            "C3  A3* {} cases against Z3 {} ({} tris against {}); per site {} cells x {} \
             corners = {} reads against {} x {} = {} -> {c3}",
            bcc_table.cases,
            cubic_table.cases,
            bcc_table.total_triangles,
            cubic_table.total_triangles,
            BCC_CELLS_PER_SITE,
            bcc_table.corners_per_cell,
            BCC_CELLS_PER_SITE * bcc_table.corners_per_cell,
            CUBIC_CELLS_PER_SITE,
            cubic_table.corners_per_cell,
            CUBIC_CELLS_PER_SITE * cubic_table.corners_per_cell
        );
        println!();

        // ── the rows ────────────────────────────────────────────────────────

        for row in &rows {
            for side in [&row.cubic, &row.bcc] {
                let is_control = matches!(side.lattice, Lattice::Cubic);
                let table = if is_control { &cubic_table } else { &bcc_table };
                let (cells, stencil, reads) = cost(side.lattice, table);
                // Relative to the cubic control, so the control's own row reads
                // 1 and 0 dB by construction rather than by omission.
                let ratio = row.cubic.hausdorff / side.hausdorff;
                let gain = AMPLITUDE_DB * ratio.log10();
                let rms_gain = AMPLITUDE_DB * (row.cubic.rms / side.rms).log10();

                run.record(&[
                    ("lattice", side.lattice.name().to_string()),
                    ("G", format!("{:.9}", side.lattice.g())),
                    ("samples", side.samples.to_string()),
                    ("field", row.field.to_string()),
                    ("hausdorff", format!("{:.9}", side.hausdorff)),
                    ("hausdorff_ratio", format!("{ratio:.6}")),
                    ("predicted_gain_db", format!("{predicted:.6}")),
                    ("measured_gain_db", format!("{gain:.6}")),
                    ("prediction_holds", row.prediction_holds.to_string()),
                    ("extraction_ms", format!("{:.4}", side.ms_median)),
                    ("case_table_size", table.cases.to_string()),
                    ("c1_holds", c1.to_string()),
                    ("c2_holds", c2.to_string()),
                    ("c3_holds", c3.to_string()),
                    // ── extras (M-273) ──
                    ("cells_per_site", cells.to_string()),
                    ("corner_reads_per_site", reads.to_string()),
                    ("corners_per_cell", table.corners_per_cell.to_string()),
                    ("density_mismatch", format!("{:.6}", row.mismatch)),
                    (
                        "distinct_up_to_symmetry",
                        table.distinct_up_to_symmetry.to_string(),
                    ),
                    ("eval_samples", EVAL_SAMPLES.to_string()),
                    ("extraction_ms_max", format!("{:.4}", side.ms_max)),
                    ("extraction_ms_min", format!("{:.4}", side.ms_min)),
                    ("extraction_repeats", TIMED_REPEATS.to_string()),
                    ("fields_improved", improved.to_string()),
                    ("fields_in_neighbourhood", near.to_string()),
                    ("fields_rms_in_neighbourhood", near_rms.to_string()),
                    ("inset", format!("{:.6}", row.inset)),
                    ("is_control", is_control.to_string()),
                    ("lattice_scale", format!("{:.9}", side.scale)),
                    (
                        "max_triangles_per_case",
                        table.max_triangles_per_case.to_string(),
                    ),
                    ("mean_linear_error", format!("{:.9}", side.mean)),
                    ("points", side.points.to_string()),
                    ("probes", PROBES.to_string()),
                    ("reconstruction_stencil", stencil.to_string()),
                    ("refine_steps", REFINE_STEPS.to_string()),
                    ("rms_error", format!("{:.9}", side.rms)),
                    ("rms_gain_db", format!("{rms_gain:.6}")),
                    ("rms_prediction_holds", row.rms_prediction_holds.to_string()),
                    ("total_triangles", table.total_triangles.to_string()),
                    ("worst_linear_error", format!("{:.9}", side.worst)),
                ]);
            }
        }
    });
}
