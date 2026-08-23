//! **P-43 — one evaluation at the cell centre as an under-sampling witness.**
//!
//! Ticket: R-042. Pre-registered before this harness existed.
//!
//! ```bash
//! cargo bench --bench experiment_p43
//! ```
//!
//! Writes `docs/experiments/p-43.csv`.
//!
//! # The witness, and exactly how it is normalised
//!
//! For every cell of the grid, the field is evaluated once at the cell centre
//! and compared against the trilinear interpolant of the eight corners
//! *evaluated at that same centre*. At the centre the eight trilinear weights
//! are all `1/8`, so the interpolant collapses to the plain mean of the eight
//! corner values and no interpolation code is needed — the arithmetic mean *is*
//! the interpolant there.
//!
//! ```text
//! residual(cell) = |f(centre) − mean(eight corners)| / h
//! ```
//!
//! The division by `h` is the "normalised by cell size" of the registration and
//! it is fixed for the whole experiment. The reason to divide rather than not:
//! the raw difference has the units of `f`, which for these fields is a length,
//! and the quantity wanted is *"how badly does one cell mis-model the field"* —
//! a fraction of a cell, not an absolute distance. For a smooth field the raw
//! difference is `O(h²·‖∇²f‖)`, so the normalised residual is `O(h·‖∇²f‖)`: it
//! still falls with resolution, which is what makes it comparable against a
//! Hausdorff distance that also falls with resolution.
//!
//! # One-sided, in the safe direction
//!
//! This is [`validate::field_bound`](isomesh::validate::field_bound_report)'s
//! discipline pointed the other way. That module samples `‖∇f‖` and says plainly
//! that **a sampled maximum is a lower bound on a supremum**, so it can prove a
//! declared bound wrong and can never prove one right. The same asymmetry holds
//! here and is the whole point: a large centre residual is a *witness* that the
//! eight corners do not determine the field inside the cell, and that is
//! conclusive. A residual of zero is not a certificate of anything — the field
//! may do whatever it likes everywhere the single centre sample did not look.
//! The witness can prove a chunk inadequate; it can never prove one adequate.
//!
//! # Four points, and what `pearson_r` is being used for
//!
//! Clause one is a correlation between the per-grid maximum residual and the
//! symmetric Hausdorff distance [`validate::accuracy`](isomesh::validate::accuracy)
//! already reports, over the four registered resolutions of one field.
//!
//! **That is four points.** `pearson_r` here is the registered decision
//! threshold and nothing more: it is not an inferential statistic, there is no
//! p-value attached to it, and no confidence interval over four paired
//! observations of a monotone sweep would mean anything. Both quantities are
//! expected to fall with `h`, so a high `r` mostly says "they fall together" —
//! which is exactly the claim being tested, and exactly why the adversary below
//! is the informative half of the run. `pearson_r_pooled` repeats the
//! computation over the eight points of both registered fields at once; pooling
//! across fields mixes two different `‖∇²f‖` scales, so it is reported as an
//! extra column and is not the clause.
//!
//! The registration fixes the clause on the **maximum**, so `pearson_r` is the
//! maximum's correlation and nothing else decides C1. Two other statistics of
//! the same population are correlated the same way and reported as extra
//! columns, because a maximum over `(n − 1)³` cells is a single order statistic
//! and a field with a crease in it — a `max`-capped gyroid, a slab's midplane —
//! puts an `O(1)` residual in one cell at every resolution: `pearson_r_mean`
//! over `centre_residual_mean` and `pearson_r_p99` over `centre_residual_p99`.
//! They are diagnostics for reading the C1 result, not alternative clauses, and
//! they are named here so it is clear they were computed alongside rather than
//! chosen afterwards.
//!
//! # The adversary
//!
//! `thin_plate` is measured beside the two registered fields because
//! `falsified_by` names its failure mode in advance: a feature that passes
//! cleanly through the cell centre without perturbing it gives a residual near
//! zero while the reconstruction is wrong. It is a slab, so its field is
//! piecewise linear away from the faces and a linear function is reproduced
//! *exactly* by its own trilinear interpolant — the residual is identically zero
//! on the overwhelming majority of cells. Its rows are reported in full and its
//! own `pearson_r` is computed the same way, but it is not part of clause one;
//! the `registered_field` column marks which rows the clause is decided on.
//!
//! # Which evaluation accounting `extra_eval_fraction` uses
//!
//! Clause two is a claim about the **count** of extra evaluations, not their
//! cost, and the count depends on an accounting choice the registration does not
//! make explicit. It says "the extra evaluations are under 15% of the corner
//! evaluations, which is the structural 1/8 plus slack", and 1/8 is the *per
//! cell* figure: one centre against a cell's eight corners.
//!
//! This crate does not evaluate corners per cell. `MarchingCubes::extract`
//! prefills one shared `n³` value grid before any cell work, so each corner is
//! evaluated exactly once and shared between the up-to-eight cells touching it.
//! Against *that* denominator the extra work is `(n − 1)³ / n³`, which rises
//! towards 1 rather than sitting at 1/8. The registered column follows the
//! crate's own accounting, because that is the number a caller would pay:
//!
//! ```text
//! extra_eval_fraction          = (n − 1)³ / n³        shared corners, what isomesh does
//! extra_eval_fraction_unshared = (n − 1)³ / 8(n − 1)³ = 1/8   per-cell corners
//! ```
//!
//! Both are reported, so the clause can be read either way from the artefact
//! instead of from an argument about which one was meant.
//!
//! # Nothing here is timed
//!
//! Since clause two counts rather than clocks, no wall clock is read anywhere in
//! this file and there is no median-of-runs to report. Every counter in the CSV
//! is one the sweep actually incremented, not a formula re-evaluated.

mod common;

use isomesh::fields::ReferenceField;
use isomesh::marching_cubes::MarchingCubes;
use isomesh::validate::{AccuracyConfig, accuracy};
use isomesh::{MeshBuffer, Sdf};

/// The resolutions named by the registration.
const RESOLUTIONS: [u32; 4] = [17, 33, 65, 129];

/// The two fields clause one is decided on.
const REGISTERED_FIELDS: [&str; 2] = ["noise_cavity", "gyroid"];

/// The adversary named in `falsified_by`, measured beside them.
const ADVERSARY: &str = "thin_plate";

/// Clause one's threshold on Pearson `r`.
const R_FLOOR: f64 = 0.7;

/// Clause two's ceiling on the extra evaluation fraction.
const EXTRA_EVAL_CEILING: f64 = 0.15;

/// A cell is called inadequate when its residual exceeds this multiple of the
/// grid's median residual.
///
/// A relative threshold rather than an absolute one, because the residual's
/// scale is a property of the field's second derivative and no fixed number
/// works for both a gyroid and a noise cavity. It degenerates honestly: when the
/// median is zero — `thin_plate`, where most cells are exactly linear — the
/// threshold is zero and the count becomes "cells with any residual at all",
/// which is the useful reading there. `centre_residual_median` is reported so
/// that degeneracy is visible in the CSV rather than hidden in this comment.
const INADEQUATE_MULTIPLE: f64 = 10.0;

/// One grid's centre-residual population, reduced to what is reported.
struct Residuals {
    /// Largest normalised residual on the grid.
    max: f64,
    /// Mean normalised residual.
    mean: f64,
    /// 99th percentile, nearest-rank.
    p99: f64,
    /// Median. The mean of the two central order statistics when the population
    /// is even, which it always is here: `n` is odd, so `(n − 1)³` is even.
    median: f64,
    /// Cells whose residual exceeds [`INADEQUATE_MULTIPLE`] times the median.
    inadequate: u64,
    /// The threshold that count was taken at.
    threshold: f64,
    /// Cell-centre evaluations, one per cell: `(n − 1)³`.
    centre_evals: u64,
    /// Corner evaluations: `n³`.
    corner_evals: u64,
}

/// One measured `(field, resolution)` cell of the matrix.
struct Measured {
    /// Field name, as `for_each_reference_field!` spells it.
    field: &'static str,
    /// Samples per axis.
    samples: u32,
    /// Grid spacing.
    cell: f64,
    /// The witness.
    residuals: Residuals,
    /// `max(mesh→field, field→mesh)`, from the crate's own accuracy report.
    hausdorff: f64,
    /// Mesh → field maximum, the direction that notices misplaced geometry.
    forward_max: f64,
    /// Field → mesh maximum, the direction that notices missing geometry.
    reverse_max: f64,
    /// Mesh → field mean.
    forward_mean: f64,
    /// Triangles the extractor produced.
    triangles: usize,
    /// Vertices the extractor produced.
    vertices: usize,
}

impl Measured {
    /// `(n − 1)³ / n³`, from the counters rather than the formula.
    ///
    /// The crate's own accounting: one shared `n³` corner grid, one centre per
    /// cell. This is the registered column.
    fn extra_eval_fraction(&self) -> f64 {
        let r = &self.residuals;
        r.centre_evals as f64 / r.corner_evals as f64
    }

    /// The same ratio against **eight corners per cell**, which is the `1/8` the
    /// registration's parenthetical describes. Structurally constant; computed
    /// from the counter anyway so it cannot drift from the sweep.
    fn extra_eval_fraction_unshared(&self) -> f64 {
        let centres = self.residuals.centre_evals as f64;
        centres / (8.0 * centres)
    }
}

/// Every cell's normalised centre residual, reduced.
///
/// The corner values are computed once into a grid and shared between the eight
/// cells that touch each of them, so this costs `n³` corner evaluations and
/// `(n − 1)³` centre evaluations and not `8·(n − 1)³` of the former.
fn centre_residuals<F>(field: &F, samples: u32, lo: [f64; 3], cell: f64) -> Residuals
where
    F: Sdf<Scalar = f64>,
{
    let n = samples as usize;
    let mut values = vec![0.0f64; n * n * n];
    for k in 0..n {
        for j in 0..n {
            for i in 0..n {
                // `i = x + y·sx + z·sx·sy`, the crate's index order.
                values[i + j * n + k * n * n] = field.sample([
                    lo[0] + cell * i as f64,
                    lo[1] + cell * j as f64,
                    lo[2] + cell * k as f64,
                ]);
            }
        }
    }

    let cells = n - 1;
    let half = cell * 0.5;
    let mut residuals = Vec::with_capacity(cells * cells * cells);
    let mut sum = 0.0f64;
    for k in 0..cells {
        for j in 0..cells {
            for i in 0..cells {
                let base = i + j * n + k * n * n;
                let corners = values[base]
                    + values[base + 1]
                    + values[base + n]
                    + values[base + n + 1]
                    + values[base + n * n]
                    + values[base + n * n + 1]
                    + values[base + n * n + n]
                    + values[base + n * n + n + 1];
                let interpolated = corners / 8.0;
                let centre = field.sample([
                    lo[0] + cell * i as f64 + half,
                    lo[1] + cell * j as f64 + half,
                    lo[2] + cell * k as f64 + half,
                ]);
                let r = (centre - interpolated).abs() / cell;
                sum += r;
                residuals.push(r);
            }
        }
    }

    // `total_cmp` is a total order on `f64`, so this is deterministic even if a
    // field hands back a NaN — which would sort to the end and be visible in
    // `max` rather than silently dropped by a partial comparison.
    residuals.sort_unstable_by(f64::total_cmp);
    let len = residuals.len();
    let max = residuals[len - 1];
    let mean = sum / len as f64;
    let p99 = residuals[(len - 1) * 99 / 100];
    let median = if len % 2 == 0 {
        (residuals[len / 2 - 1] + residuals[len / 2]) / 2.0
    } else {
        residuals[len / 2]
    };
    let threshold = INADEQUATE_MULTIPLE * median;
    let inadequate = residuals.iter().filter(|&&r| r > threshold).count() as u64;

    Residuals {
        max,
        mean,
        p99,
        median,
        inadequate,
        threshold,
        centre_evals: len as u64,
        corner_evals: (n * n * n) as u64,
    }
}

/// The witness and the correlate, for one field at one resolution.
fn measure<F>(
    field: &F,
    name: &'static str,
    samples: u32,
    mc: &mut MarchingCubes<f64>,
    out: &mut MeshBuffer<f64>,
) -> Measured
where
    F: ReferenceField + Sdf<Scalar = f64>,
{
    let (shape, lo, cell) = common::grid(field, samples);
    let residuals = centre_residuals(field, samples, lo, cell);

    out.reset();
    mc.extract(field, &shape, lo, cell, out)
        .expect("marching cubes extraction");

    let cfg = AccuracyConfig::from_cell_size(cell).expect("positive cell size");
    let report = accuracy(&out.positions, &out.indices, field, &shape, lo, &cfg)
        .expect("accuracy over the extraction grid");

    Measured {
        field: name,
        samples,
        cell,
        residuals,
        hausdorff: report.symmetric_hausdorff(),
        forward_max: report.mesh_to_field.max,
        reverse_max: report.field_to_mesh.max,
        forward_mean: report.mesh_to_field.mean,
        triangles: out.triangle_count(),
        vertices: out.vertex_count(),
    }
}

/// Pearson's `r` over paired samples.
///
/// Zero when either sample has no spread, which is the honest answer: a constant
/// column carries no covariance to detect, and returning a NaN would propagate
/// into the CSV as a column nobody can compare against `R_FLOOR`.
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

/// A statistic of one measured grid, as a plain function pointer.
///
/// One projector type for both the correlation and the decay exponent, so the
/// two cannot end up describing different columns.
type Stat = fn(&Measured) -> f64;

/// `r` between one statistic and the symmetric Hausdorff distance, over every
/// row belonging to `fields`.
fn correlate(rows: &[Measured], fields: &[&str], pick: Stat) -> f64 {
    let mut xs = Vec::new();
    let mut ys = Vec::new();
    for row in rows {
        if fields.contains(&row.field) {
            xs.push(pick(row));
            ys.push(row.hausdorff);
        }
    }
    pearson(&xs, &ys)
}

/// The registered statistic: the per-grid maximum residual.
fn residual_max(m: &Measured) -> f64 {
    m.residuals.max
}

/// Diagnostic statistic: the per-grid mean residual.
fn residual_mean(m: &Measured) -> f64 {
    m.residuals.mean
}

/// Diagnostic statistic: the per-grid 99th-percentile residual.
fn residual_p99(m: &Measured) -> f64 {
    m.residuals.p99
}

/// The correlate itself, so its own decay can be measured on the same footing.
fn hausdorff(m: &Measured) -> f64 {
    m.hausdorff
}

/// Least-squares slope of `ln(statistic)` against `ln(h)` over one field's
/// resolutions — the **observed convergence order** of that statistic.
///
/// This is the number that explains a correlation rather than merely reporting
/// it. A statistic whose exponent is near zero does not fall with resolution at
/// all, so it cannot track a Hausdorff distance that does, whatever `r` happens
/// to come out at over four points.
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

fn main() {
    let prereg = isomesh::experiment!("P-43");

    let mut mc = MarchingCubes::<f64>::new();
    let mut out = MeshBuffer::<f64>::new();
    let mut rows: Vec<Measured> = Vec::new();

    // `for_each_reference_field!` inlines its body once per field, so a `return`
    // in here would leave the sweep silently truncated. Selection is by name.
    isomesh::for_each_reference_field!(f64, |name, field| {
        if REGISTERED_FIELDS.contains(&name) || name == ADVERSARY {
            for samples in RESOLUTIONS {
                let row = measure(&field, name, samples, &mut mc, &mut out);
                println!(
                    "{name:>13} {samples:>4}³  residual max {:.6e} mean {:.6e}  \
                     hausdorff {:.6e}  ({} tri)",
                    row.residuals.max, row.residuals.mean, row.hausdorff, row.triangles
                );
                rows.push(row);
            }
        }
    });

    let pooled = correlate(&rows, &REGISTERED_FIELDS, residual_max);
    let pooled_mean = correlate(&rows, &REGISTERED_FIELDS, residual_mean);
    let pooled_p99 = correlate(&rows, &REGISTERED_FIELDS, residual_p99);
    let per_field = |f: &str| correlate(&rows, &[f], residual_max);
    let r_noise = per_field("noise_cavity");
    let r_gyroid = per_field("gyroid");
    let r_adversary = per_field(ADVERSARY);
    let worst_registered = r_noise.min(r_gyroid);
    let worst_extra_fraction = rows
        .iter()
        .map(Measured::extra_eval_fraction)
        .fold(0.0f64, f64::max);

    println!(
        "\nmax  r(noise_cavity) {r_noise:.4}  r(gyroid) {r_gyroid:.4}  \
         pooled {pooled:.4}  r({ADVERSARY}) {r_adversary:.4}"
    );
    for stat in [
        ("mean", residual_mean as Stat, pooled_mean),
        ("p99", residual_p99 as Stat, pooled_p99),
    ] {
        let (label, pick, pooled_stat) = stat;
        println!(
            "{label:>4} r(noise_cavity) {:.4}  r(gyroid) {:.4}  pooled {pooled_stat:.4}  \
             r({ADVERSARY}) {:.4}",
            correlate(&rows, &["noise_cavity"], pick),
            correlate(&rows, &["gyroid"], pick),
            correlate(&rows, &[ADVERSARY], pick),
        );
    }
    println!("\nobserved decay exponents in h (log-log slope over the four resolutions):");
    for field in [REGISTERED_FIELDS[0], REGISTERED_FIELDS[1], ADVERSARY] {
        println!(
            "{field:>13}  max {:+.3}  mean {:+.3}  p99 {:+.3}  hausdorff {:+.3}",
            decay_exponent(&rows, field, residual_max),
            decay_exponent(&rows, field, residual_mean),
            decay_exponent(&rows, field, residual_p99),
            decay_exponent(&rows, field, hausdorff),
        );
    }
    println!(
        "\nC1 r >= {R_FLOOR} on both registered fields: {}  (worst {worst_registered:.4})",
        if worst_registered >= R_FLOOR {
            "HELD"
        } else {
            "FALSIFIED"
        }
    );
    println!(
        "C2 extra_eval_fraction < {EXTRA_EVAL_CEILING}: {}  (worst {worst_extra_fraction:.6}, \
         unshared-corner accounting would read 0.125)\n",
        if worst_extra_fraction < EXTRA_EVAL_CEILING {
            "HELD"
        } else {
            "FALSIFIED"
        }
    );

    common::experiment::run(prereg, |run| {
        for row in &rows {
            let r = &row.residuals;
            let registered = REGISTERED_FIELDS.contains(&row.field);
            run.record(&[
                ("field", row.field.to_string()),
                ("samples_per_axis", row.samples.to_string()),
                ("centre_residual_max", format!("{:.9e}", r.max)),
                ("centre_residual_mean", format!("{:.9e}", r.mean)),
                ("symmetric_hausdorff", format!("{:.9e}", row.hausdorff)),
                ("pearson_r", format!("{:.6}", per_field(row.field))),
                (
                    "extra_eval_fraction",
                    format!("{:.6}", row.extra_eval_fraction()),
                ),
                ("registered_field", registered.to_string()),
                (
                    "extra_eval_fraction_unshared",
                    format!("{:.6}", row.extra_eval_fraction_unshared()),
                ),
                (
                    "pearson_r_mean",
                    format!("{:.6}", correlate(&rows, &[row.field], residual_mean)),
                ),
                (
                    "pearson_r_p99",
                    format!("{:.6}", correlate(&rows, &[row.field], residual_p99)),
                ),
                ("pearson_r_mean_pooled", format!("{pooled_mean:.6}")),
                ("pearson_r_p99_pooled", format!("{pooled_p99:.6}")),
                (
                    "decay_exponent_max",
                    format!("{:.6}", decay_exponent(&rows, row.field, residual_max)),
                ),
                (
                    "decay_exponent_mean",
                    format!("{:.6}", decay_exponent(&rows, row.field, residual_mean)),
                ),
                (
                    "decay_exponent_p99",
                    format!("{:.6}", decay_exponent(&rows, row.field, residual_p99)),
                ),
                (
                    "decay_exponent_hausdorff",
                    format!("{:.6}", decay_exponent(&rows, row.field, hausdorff)),
                ),
                ("pearson_r_pooled", format!("{pooled:.6}")),
                ("centre_residual_p99", format!("{:.9e}", r.p99)),
                ("centre_residual_median", format!("{:.9e}", r.median)),
                ("inadequate_threshold", format!("{:.9e}", r.threshold)),
                ("inadequate_cells", r.inadequate.to_string()),
                (
                    "inadequate_fraction",
                    format!("{:.9}", r.inadequate as f64 / r.centre_evals as f64),
                ),
                ("centre_evals", r.centre_evals.to_string()),
                ("corner_evals", r.corner_evals.to_string()),
                ("hausdorff_forward_max", format!("{:.9e}", row.forward_max)),
                ("hausdorff_reverse_max", format!("{:.9e}", row.reverse_max)),
                (
                    "hausdorff_forward_mean",
                    format!("{:.9e}", row.forward_mean),
                ),
                ("cell_size", format!("{:.9e}", row.cell)),
                ("triangles", row.triangles.to_string()),
                ("vertices", row.vertices.to_string()),
            ]);
        }
    });
}
