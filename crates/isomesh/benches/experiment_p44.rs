//! **P-44 — the mean centre residual, tested out of sample.**
//!
//! Ticket: R-042a. Pre-registered before this harness existed.
//!
//! ```bash
//! cargo bench --bench experiment_p44
//! ```
//!
//! Writes `docs/experiments/p-44.csv`.
//!
//! # Why this experiment exists at all
//!
//! P-43 registered the **maximum** normalised centre residual and was falsified
//! (✗29): `r(gyroid) = −0.847` against a required `+0.7`. The post-mortem in that
//! run's own extra columns found the **mean** correlating at `r = 0.983`, `0.984`
//! and `0.9998` on the same three fields, with log-log decay exponents agreeing
//! to within `0.05`.
//!
//! **That number was read off the data that killed its predecessor**, which makes
//! it a description of three fields and not yet a claim about a witness. This
//! harness moves it somewhere it can fail: the four reference fields P-43 never
//! touched — `sphere`, `torus`, `box_exact`, `csg_difference` — at the same four
//! resolutions. Nothing here is tuned on the P-43 rows; the statistic, the
//! normalisation and the grid convention are all carried over unchanged so that
//! a hit cannot be manufactured by redefining the quantity.
//!
//! # The witness, normalised exactly as P-43 normalised it
//!
//! Unchanged from P-43, deliberately, so the two CSVs are comparable row for row:
//!
//! ```text
//! residual(cell) = |f(centre) − mean(eight corners)| / h
//! ```
//!
//! At the cell centre all eight trilinear weights are `1/8`, so the trilinear
//! interpolant of the corners *is* their arithmetic mean and no interpolation
//! code is involved. The division by `h` is the registration's "normalised by
//! cell size" and it is the same divisor P-43 used. The registered statistic is
//! the per-grid **mean** of that population.
//!
//! # What the four new fields test that P-43's three could not
//!
//! P-43's explanation for its own failure was structural: a `C¹` **crease**
//! makes `|f(centre) − mean(corners)|` an `O(h)` quantity rather than `O(h²)`,
//! so dividing by `h` leaves a constant and the maximum stops falling with
//! resolution. All three of its fields had one — `gyroid` is
//! `max(Gyroid, Sphere)`, a CSG seam; `thin_plate` is an exact box; `noise_cavity`
//! has an unbounded gradient.
//!
//! This set splits on exactly that property, which makes it an out-of-sample test
//! of the *mechanism* and not only of the statistic:
//!
//! - `sphere` and `torus` are smooth on their domains — no crease anywhere the
//!   surface is. If the crease account is right, their **maximum** residual
//!   should decay with an exponent near `1` where P-43's fields sat near `0`.
//! - `box_exact` and `csg_difference` have creases: box edges and corners, and a
//!   concave `max` seam. Their maximum should behave like P-43's did.
//!
//! The registered clauses are all about the mean. `pearson_r_max`,
//! `decay_exponent_max` and `centre_residual_max` are extra columns carried for
//! this comparison, and they are named here so it is clear they were computed
//! alongside rather than chosen afterwards.
//!
//! # Clause three is wall clock, and how it is measured
//!
//! P-43's cost clause was arithmetically wrong — it priced the witness at `1/8`
//! of the corner evaluations, which is the *per-cell* accounting, while this
//! crate prefills one shared sample grid and evaluates every corner exactly once.
//! The true evaluation ratio is `(n − 1)³ / n³`, which rises toward one. The
//! registration therefore restates the claim on time: the witness must cost under
//! `0.5×` a Marching Cubes extraction of the same grid.
//!
//! Two readings of "computing the witness" exist and both are reported, with the
//! **conservative** one registered:
//!
//! - `witness_ns` — the whole thing from a cold sample grid: fill `n³` corners,
//!   then one streaming pass over `(n − 1)³` cells. This is the registered
//!   column. It is the standalone cost, and it double-counts the corner grid
//!   against a caller who is about to mesh anyway.
//! - `witness_incremental_ns` — the centre pass alone, with the corner grid
//!   already resident. This is what a caller integrated into `extract` would
//!   actually pay, since `MarchingCubes` has already filled that grid.
//!
//! Registering the larger number can only make clause three harder to pass, which
//! is the direction an honest cost claim should err in.
//!
//! The witness allocates nothing and sorts nothing: the mean is a streaming sum
//! and the maximum is one `fmax` per cell, so there is no `O(n log n)` term
//! hiding in the timing. The corner buffer is resized once per grid, outside
//! every timed region.
//!
//! # Timing discipline
//!
//! Buffers — the mesher, the sink and the corner grid — are allocated once and
//! reused. One untimed warm-up per measurement, so no timed run pays the first
//! run's page faults, and that warm-up is also the run whose *value* is reported.
//! Then five timed runs and the **median**; never a mean, which one slow run
//! would drag. [`std::hint::black_box`] consumes the witness's result so it
//! cannot be optimised out of existence, and the extraction's output lands in a
//! buffer that is read afterwards.
//!
//! # Four points, and what the correlations are for
//!
//! Each `pearson_r_mean` is four paired observations of a monotone refinement
//! sweep. It is the registered decision threshold and nothing more: no p-value,
//! no confidence interval, and none would mean anything at that sample size.
//! `exponent_gap` is the stronger of the two statistical clauses precisely
//! because it compares *rates* rather than levels — two quantities can correlate
//! at `r > 0.99` while converging at different orders, and an exponent gap
//! notices that where `r` does not.
//!
//! # A mean statistic against a max metric, carried as a diagnostic
//!
//! The registered correlate is the symmetric Hausdorff distance, which is a
//! **supremum** over both sampled directions. The registered witness statistic
//! is a **mean** over `(n − 1)³` cells. Those are different kinds of quantity,
//! and a mean has no obligation to converge at the order of a maximum: a feature
//! occupying an `O(h)` fraction of the domain contributes `O(h)` to a max and
//! `O(h)·O(h) = O(h²)` to a mean.
//!
//! So `mean_absolute_error` — [`AccuracyReport::mean_absolute_error`], the
//! mesh→field mean the crate already computes on the same report — is correlated
//! against the same witness with the same code path, and its own decay exponent
//! and gap are reported. This is a **diagnostic for reading clause two**, not an
//! alternative clause: clause two is decided against the Hausdorff exponent and
//! nothing else. Naming it here rather than reaching for it afterwards is the
//! same discipline the registration itself applies to P-43's mean.

mod common;

use std::hint::black_box;
use std::time::Instant;

use isomesh::fields::ReferenceField;
use isomesh::marching_cubes::MarchingCubes;
use isomesh::validate::{AccuracyConfig, accuracy};
use isomesh::{MeshBuffer, RuntimeShape3, Sdf};

/// The resolutions named by the registration.
const RESOLUTIONS: [u32; 4] = [17, 33, 65, 129];

/// The four fields P-43 never touched. Every clause is decided on all four.
const FIELDS: [&str; 4] = ["sphere", "torus", "box_exact", "csg_difference"];

/// Clause one's floor on Pearson `r`.
const R_FLOOR: f64 = 0.9;

/// Clause two's ceiling on `|exponent_mean − exponent_hausdorff|`.
const EXPONENT_GAP_CEILING: f64 = 0.15;

/// Clause three's ceiling on `witness_ns / extract_ns`.
const COST_RATIO_CEILING: f64 = 0.5;

/// Timed repetitions per measurement. Odd, so the median is an observation.
const REPS: usize = 5;

/// What one streaming witness pass found.
///
/// Both statistics come out of the same pass. The mean is the registered one; the
/// maximum rides along for one `fmax` per cell and is what tests P-43's crease
/// explanation on fields that do and do not have a crease.
#[derive(Clone, Copy)]
struct Stats {
    /// Per-grid mean normalised residual — the registered statistic.
    mean: f64,
    /// Per-grid maximum normalised residual — P-43's falsified statistic,
    /// carried for the crease comparison.
    max: f64,
}

/// Buffers allocated once and reused, so no timed region contains an allocation.
struct Rig {
    /// The mesher, which keeps its own value and edge caches between calls.
    mc: MarchingCubes<f64>,
    /// The extraction sink.
    out: MeshBuffer<f64>,
    /// The `n³` corner grid the witness fills and then reads.
    corners: Vec<f64>,
}

impl Rig {
    /// Nothing allocated yet.
    fn new() -> Self {
        Self {
            mc: MarchingCubes::new(),
            out: MeshBuffer::new(),
            corners: Vec::new(),
        }
    }

    /// Size the corner grid for `samples` per axis. Never inside a timed region.
    fn reserve(&mut self, samples: u32) {
        let n = samples as usize;
        self.corners.resize(n * n * n, 0.0);
    }

    /// Fill the shared corner grid: one field evaluation per **sample**.
    ///
    /// `i = x + y·sx + z·sx·sy`, the crate's index order, so this grid is laid
    /// out exactly as `MarchingCubes`' own prefill.
    fn fill_corners<F>(&mut self, field: &F, samples: u32, lo: [f64; 3], cell: f64)
    where
        F: Sdf<Scalar = f64>,
    {
        let n = samples as usize;
        for k in 0..n {
            for j in 0..n {
                let z = lo[2] + cell * k as f64;
                let y = lo[1] + cell * j as f64;
                let row = j * n + k * n * n;
                for i in 0..n {
                    self.corners[row + i] = field.sample([lo[0] + cell * i as f64, y, z]);
                }
            }
        }
    }

    /// One field evaluation per **cell**, reduced against the corners already in
    /// `self.corners`.
    ///
    /// Streaming: no allocation, no sort, one running sum and one running max.
    fn centre_pass<F>(&self, field: &F, samples: u32, lo: [f64; 3], cell: f64) -> Stats
    where
        F: Sdf<Scalar = f64>,
    {
        let n = samples as usize;
        let cells = n - 1;
        let half = cell * 0.5;
        let plane = n * n;
        let mut sum = 0.0f64;
        let mut max = 0.0f64;
        for k in 0..cells {
            let z = lo[2] + cell * k as f64 + half;
            for j in 0..cells {
                let y = lo[1] + cell * j as f64 + half;
                let base = j * n + k * plane;
                for i in 0..cells {
                    let c = base + i;
                    let corners = self.corners[c]
                        + self.corners[c + 1]
                        + self.corners[c + n]
                        + self.corners[c + n + 1]
                        + self.corners[c + plane]
                        + self.corners[c + plane + 1]
                        + self.corners[c + plane + n]
                        + self.corners[c + plane + n + 1];
                    let centre = field.sample([lo[0] + cell * i as f64 + half, y, z]);
                    let r = (centre - corners / 8.0).abs() / cell;
                    sum += r;
                    max = max.max(r);
                }
            }
        }
        let count = cells * cells * cells;
        Stats {
            mean: sum / count as f64,
            max,
        }
    }

    /// The standalone witness: corner grid from cold, then the centre pass.
    fn witness<F>(&mut self, field: &F, samples: u32, lo: [f64; 3], cell: f64) -> Stats
    where
        F: Sdf<Scalar = f64>,
    {
        self.fill_corners(field, samples, lo, cell);
        self.centre_pass(field, samples, lo, cell)
    }

    /// One Marching Cubes extraction of the same grid, into the reused sink.
    fn extract<F>(&mut self, field: &F, shape: &RuntimeShape3, lo: [f64; 3], cell: f64)
    where
        F: Sdf<Scalar = f64>,
    {
        self.out.reset();
        self.mc
            .extract(field, shape, lo, cell, &mut self.out)
            .expect("marching cubes extraction");
    }
}

/// Median nanoseconds over [`REPS`] timed runs.
///
/// The median and never the mean: one descheduled run would drag an average and
/// there is no reason to let it.
fn median_ns(mut run: impl FnMut()) -> f64 {
    let mut runs = Vec::with_capacity(REPS);
    for _ in 0..REPS {
        let t = Instant::now();
        run();
        runs.push(t.elapsed().as_secs_f64() * 1e9);
    }
    runs.sort_by(f64::total_cmp);
    runs[runs.len() / 2]
}

/// One measured `(field, resolution)` cell of the matrix.
struct Measured {
    /// Field name, as `for_each_reference_field!` spells it.
    field: &'static str,
    /// Samples per axis.
    samples: u32,
    /// Grid spacing `h`.
    cell: f64,
    /// The witness, from the untimed warm-up run.
    stats: Stats,
    /// `max(mesh→field, field→mesh)` from the crate's own accuracy report.
    hausdorff: f64,
    /// Mesh → field maximum: misplaced geometry.
    forward_max: f64,
    /// Field → mesh maximum: missing geometry.
    reverse_max: f64,
    /// Mesh → field **mean**, the crate's `mean_absolute_error`. A mean-type
    /// error metric, to sit beside the max-type Hausdorff.
    mesh_mean: f64,
    /// Both accuracy directions produced samples.
    has_coverage: bool,
    /// Median nanoseconds for the standalone witness. Registered.
    witness_ns: f64,
    /// Median nanoseconds for the centre pass alone, corners resident.
    incremental_ns: f64,
    /// Median nanoseconds for the corner fill alone.
    corner_fill_ns: f64,
    /// Median nanoseconds for one Marching Cubes extraction. Registered.
    extract_ns: f64,
    /// Cells, `(n − 1)³`.
    cells: u64,
    /// Samples, `n³`.
    corner_evals: u64,
    /// Triangles produced.
    triangles: usize,
    /// Vertices produced.
    vertices: usize,
}

impl Measured {
    /// Clause three's quantity: standalone witness against extraction.
    fn cost_ratio(&self) -> f64 {
        self.witness_ns / self.extract_ns
    }

    /// The same ratio for a witness folded into an extraction that has already
    /// filled the sample grid.
    fn incremental_cost_ratio(&self) -> f64 {
        self.incremental_ns / self.extract_ns
    }

    /// `(n − 1)³ / n³`, the evaluation ratio P-43's cost clause got wrong.
    fn extra_eval_fraction(&self) -> f64 {
        self.cells as f64 / self.corner_evals as f64
    }
}

/// Measure one field at one resolution: value first, then time.
fn measure<F>(field: &F, name: &'static str, samples: u32, rig: &mut Rig) -> Measured
where
    F: ReferenceField + Sdf<Scalar = f64>,
{
    let (shape, lo, cell) = common::grid(field, samples);
    rig.reserve(samples);

    // Untimed warm-up. Every buffer is resident afterwards, and this is the run
    // whose value is reported — the timed runs recompute the same thing.
    let stats = rig.witness(field, samples, lo, cell);
    rig.extract(field, &shape, lo, cell);

    let witness_ns = median_ns(|| {
        black_box(rig.witness(field, samples, lo, cell));
    });
    let corner_fill_ns = median_ns(|| rig.fill_corners(field, samples, lo, cell));
    let incremental_ns = median_ns(|| {
        black_box(rig.centre_pass(field, samples, lo, cell));
    });
    let extract_ns = median_ns(|| rig.extract(field, &shape, lo, cell));

    let cfg = AccuracyConfig::from_cell_size(cell).expect("positive cell size");
    let report = accuracy(
        &rig.out.positions,
        &rig.out.indices,
        field,
        &shape,
        lo,
        &cfg,
    )
    .expect("accuracy over the extraction grid");

    let n = u64::from(samples);
    Measured {
        field: name,
        samples,
        cell,
        stats,
        hausdorff: report.symmetric_hausdorff(),
        forward_max: report.mesh_to_field.max,
        reverse_max: report.field_to_mesh.max,
        mesh_mean: report.mean_absolute_error(),
        has_coverage: report.has_coverage(),
        witness_ns,
        incremental_ns,
        corner_fill_ns,
        extract_ns,
        cells: (n - 1) * (n - 1) * (n - 1),
        corner_evals: n * n * n,
        triangles: rig.out.triangle_count(),
        vertices: rig.out.vertex_count(),
    }
}

/// A statistic of one measured grid, as a plain function pointer.
///
/// One projector type for the correlation and the decay exponent, so the two
/// cannot end up describing different columns.
type Stat = fn(&Measured) -> f64;

/// The registered statistic: the per-grid mean residual.
fn residual_mean(m: &Measured) -> f64 {
    m.stats.mean
}

/// P-43's falsified statistic, carried for the crease comparison.
fn residual_max(m: &Measured) -> f64 {
    m.stats.max
}

/// The registered correlate: a **max**-type error metric.
fn hausdorff(m: &Measured) -> f64 {
    m.hausdorff
}

/// A **mean**-type error metric, for the diagnostic in the module docs: the
/// witness statistic is a mean, and a mean cannot be expected to converge at the
/// order of a maximum.
fn mesh_mean(m: &Measured) -> f64 {
    m.mesh_mean
}

/// Pearson's `r` over paired samples.
///
/// Zero when either sample has no spread: a constant column carries no
/// covariance, and a NaN in the CSV is a column nobody can compare to a floor.
fn pearson(xs: &[f64], ys: &[f64]) -> f64 {
    assert_eq!(xs.len(), ys.len(), "paired samples");
    let n = xs.len() as f64;
    let mx = xs.iter().sum::<f64>() / n;
    let my = ys.iter().sum::<f64>() / n;
    let mut sxy = 0.0f64;
    let mut sxx = 0.0f64;
    let mut syy = 0.0f64;
    for (&x, &y) in xs.iter().zip(ys) {
        let (dx, dy) = (x - mx, y - my);
        sxy += dx * dy;
        sxx += dx * dx;
        syy += dy * dy;
    }
    if sxx <= 0.0 || syy <= 0.0 {
        return 0.0;
    }
    sxy / (sxx * syy).sqrt()
}

/// `r` between one statistic and one error metric, over the rows of `field`.
fn correlate(rows: &[Measured], field: &str, pick: Stat, against: Stat) -> f64 {
    let mut xs = Vec::new();
    let mut ys = Vec::new();
    for row in rows {
        if row.field == field {
            xs.push(pick(row));
            ys.push(against(row));
        }
    }
    pearson(&xs, &ys)
}

/// Least-squares slope of `ln(statistic)` against `ln(h)` over one field's
/// resolutions — the **observed convergence order** of that statistic.
///
/// Zero when any observation is non-positive, so a `ln(0)` never reaches the CSV
/// as `-inf`.
fn decay_exponent(rows: &[Measured], field: &str, pick: Stat) -> f64 {
    let mut xs = Vec::new();
    let mut ys = Vec::new();
    for row in rows {
        if row.field == field {
            let v = pick(row);
            if v <= 0.0 {
                return 0.0;
            }
            xs.push(row.cell.ln());
            ys.push(v.ln());
        }
    }
    let n = xs.len() as f64;
    let mx = xs.iter().sum::<f64>() / n;
    let my = ys.iter().sum::<f64>() / n;
    let mut sxy = 0.0f64;
    let mut sxx = 0.0f64;
    for (&x, &y) in xs.iter().zip(&ys) {
        sxy += (x - mx) * (y - my);
        sxx += (x - mx) * (x - mx);
    }
    if sxx <= 0.0 {
        return 0.0;
    }
    sxy / sxx
}

/// Gap between the mean residual's decay exponent and one metric's.
///
/// Clause two's quantity when `against` is [`hausdorff`]; the diagnostic when it
/// is [`mesh_mean`].
fn exponent_gap(rows: &[Measured], field: &str, against: Stat) -> f64 {
    (decay_exponent(rows, field, residual_mean) - decay_exponent(rows, field, against)).abs()
}

fn main() {
    let prereg = isomesh::experiment!("P-44");

    let mut rig = Rig::new();
    let mut rows: Vec<Measured> = Vec::new();

    // `for_each_reference_field!` inlines its body once per field, so a `return`
    // in here would truncate the sweep silently. Selection is by name.
    isomesh::for_each_reference_field!(f64, |name, field| {
        if FIELDS.contains(&name) {
            for samples in RESOLUTIONS {
                let row = measure(&field, name, samples, &mut rig);
                println!(
                    "{name:>15} {samples:>4}³  mean {:.6e}  max {:.6e}  hausdorff {:.6e}  \
                     witness {:.3} ms  extract {:.3} ms  ratio {:.3}",
                    row.stats.mean,
                    row.stats.max,
                    row.hausdorff,
                    row.witness_ns / 1e6,
                    row.extract_ns / 1e6,
                    row.cost_ratio(),
                );
                rows.push(row);
            }
        }
    });

    println!();
    let mut worst_r = f64::INFINITY;
    let mut worst_gap = 0.0f64;
    let mut worst_gap_mae = 0.0f64;
    for field in FIELDS {
        let r_mean = correlate(&rows, field, residual_mean, hausdorff);
        let gap = exponent_gap(&rows, field, hausdorff);
        let gap_mae = exponent_gap(&rows, field, mesh_mean);
        worst_r = worst_r.min(r_mean);
        worst_gap = worst_gap.max(gap);
        worst_gap_mae = worst_gap_mae.max(gap_mae);
        println!(
            "{field:>15}  r(mean) {r_mean:+.4}  exponent: witness {:+.3}  hausdorff {:+.3} \
             (gap {gap:.3})  mae {:+.3} (gap {gap_mae:.3})",
            decay_exponent(&rows, field, residual_mean),
            decay_exponent(&rows, field, hausdorff),
            decay_exponent(&rows, field, mesh_mean),
        );
    }
    println!("\nP-43's falsified max statistic, on fields that do and do not have a crease:");
    for field in FIELDS {
        println!(
            "{field:>15}  r(max) {:+.4}  exponent(max) {:+.3}",
            correlate(&rows, field, residual_max, hausdorff),
            decay_exponent(&rows, field, residual_max),
        );
    }

    let worst_ratio = rows.iter().map(Measured::cost_ratio).fold(0.0f64, f64::max);
    let worst_incremental = rows
        .iter()
        .map(Measured::incremental_cost_ratio)
        .fold(0.0f64, f64::max);

    let verdict = |held: bool| if held { "HELD" } else { "FALSIFIED" };
    println!(
        "\nC1 r(mean) >= {R_FLOOR} on all four: {}  (worst {worst_r:+.4})",
        verdict(worst_r >= R_FLOOR)
    );
    println!(
        "C2 exponent gap <= {EXPONENT_GAP_CEILING} on all four: {}  (worst {worst_gap:.4}; \
         against the mean-type metric the worst gap is {worst_gap_mae:.4})",
        verdict(worst_gap <= EXPONENT_GAP_CEILING)
    );
    println!(
        "C3 witness_cost_ratio < {COST_RATIO_CEILING}: {}  (worst {worst_ratio:.4}; \
         incremental reading worst {worst_incremental:.4})\n",
        verdict(worst_ratio < COST_RATIO_CEILING)
    );

    common::experiment::run(prereg, |run| {
        for row in &rows {
            let field = row.field;
            run.record(&[
                ("field", field.to_string()),
                ("samples_per_axis", row.samples.to_string()),
                ("centre_residual_mean", format!("{:.9e}", row.stats.mean)),
                ("symmetric_hausdorff", format!("{:.9e}", row.hausdorff)),
                (
                    "pearson_r_mean",
                    format!("{:.6}", correlate(&rows, field, residual_mean, hausdorff)),
                ),
                (
                    "decay_exponent_mean",
                    format!("{:.6}", decay_exponent(&rows, field, residual_mean)),
                ),
                (
                    "decay_exponent_hausdorff",
                    format!("{:.6}", decay_exponent(&rows, field, hausdorff)),
                ),
                (
                    "exponent_gap",
                    format!("{:.6}", exponent_gap(&rows, field, hausdorff)),
                ),
                ("witness_ns", format!("{:.1}", row.witness_ns)),
                ("extract_ns", format!("{:.1}", row.extract_ns)),
                ("witness_cost_ratio", format!("{:.6}", row.cost_ratio())),
                ("centre_residual_max", format!("{:.9e}", row.stats.max)),
                (
                    "pearson_r_max",
                    format!("{:.6}", correlate(&rows, field, residual_max, hausdorff)),
                ),
                ("mean_absolute_error", format!("{:.9e}", row.mesh_mean)),
                (
                    "pearson_r_mean_vs_mae",
                    format!("{:.6}", correlate(&rows, field, residual_mean, mesh_mean)),
                ),
                (
                    "decay_exponent_mae",
                    format!("{:.6}", decay_exponent(&rows, field, mesh_mean)),
                ),
                (
                    "exponent_gap_mae",
                    format!("{:.6}", exponent_gap(&rows, field, mesh_mean)),
                ),
                (
                    "decay_exponent_max",
                    format!("{:.6}", decay_exponent(&rows, field, residual_max)),
                ),
                (
                    "witness_incremental_ns",
                    format!("{:.1}", row.incremental_ns),
                ),
                (
                    "witness_incremental_cost_ratio",
                    format!("{:.6}", row.incremental_cost_ratio()),
                ),
                ("corner_fill_ns", format!("{:.1}", row.corner_fill_ns)),
                (
                    "extra_eval_fraction",
                    format!("{:.6}", row.extra_eval_fraction()),
                ),
                (
                    "witness_ns_per_cell",
                    format!("{:.4}", row.witness_ns / row.cells as f64),
                ),
                (
                    "extract_ns_per_sample",
                    format!("{:.4}", row.extract_ns / row.corner_evals as f64),
                ),
                ("hausdorff_forward_max", format!("{:.9e}", row.forward_max)),
                ("hausdorff_reverse_max", format!("{:.9e}", row.reverse_max)),
                ("has_coverage", row.has_coverage.to_string()),
                ("cell_size", format!("{:.9e}", row.cell)),
                ("cells", row.cells.to_string()),
                ("corner_evals", row.corner_evals.to_string()),
                ("triangles", row.triangles.to_string()),
                ("vertices", row.vertices.to_string()),
                ("timed_reps", REPS.to_string()),
            ]);
        }
    });
}
