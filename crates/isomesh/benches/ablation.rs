//! **One algorithm, two vertex rules, eight fields, one process.**
//!
//! Ticket: X-002. This is the seam's acceptance and its first use.
//!
//! # What is being compared, and why it is a clean comparison
//!
//! [`DualContouring`] is generic over its [`VertexRule`], defaulting to its QEF rule.
//! Swapping in [`Centroid`] — Surface Nets' rule — changes **one thing**: where
//! the vertex goes inside a cell. The cell classification, the quad walk, the
//! sign conventions, the buffer reuse and the winding are the same lines of code
//! in both arms, because they *are* the same lines of code.
//!
//! That is what a comparison between two extractors cannot give you. Surface
//! Nets and Dual Contouring already differ by this rule, but they are also two
//! structs with two `extract` methods, so a difference between them is a
//! difference between two implementations. Here the difference is the rule and
//! nothing else.
//!
//! # No branch in the hot loop, and that is a type-system property
//!
//! The rule is a type parameter, not a field to match on, so each arm is
//! monomorphised and the compiler emits one placement per instantiation. There
//! is no value to test and therefore no test to emit. `the_ablation_arms_are_not_branches`
//! in `src/dual_contouring/tests.rs` holds the structural half of that claim;
//! this bench is the measured half.
//!
//! # Running it
//!
//! ```bash
//! cargo bench --bench ablation
//! ```
//!
//! Writes `docs/measurements/ablation.csv`. Guarded on `--bench` for the same
//! reason every bench here is: `--all-targets` re-selects bench targets even
//! with `test = false`, and a debug-build run would overwrite committed evidence
//! with numbers from whichever machine ran the tests.

use std::fmt::Write as _;
use std::fs;
use std::hint::black_box;
use std::path::PathBuf;
use std::time::Instant;

mod common;

use isomesh::dual::VertexRule;
use isomesh::dual_contouring::DualContouring;
#[cfg(feature = "experimental")]
use isomesh::experimental::ProbabilisticQuadric;
use isomesh::fields::ReferenceField;
use isomesh::surface_nets::Centroid;
use isomesh::validate::{
    AccuracyConfig, ValidateConfig, accuracy, self_intersections, validate_indexed,
};
use isomesh::{MeshBuffer, Sdf};

/// `f64`, so the accuracy column is measuring the rule rather than the scalar.
type Scalar = f64;

/// The resolutions each arm runs at.
const RESOLUTIONS: [u32; 3] = [17, 33, 65];

/// Untimed runs first, so the median measures steady-state re-meshing.
const WARMUP_RUNS: u32 = 2;

/// Timed runs per cell of the table.
const TIMED_RUNS: u32 = 5;

/// Self-intersection counting is `O(n²)` in triangles, so it runs at one
/// resolution — the same one the shootout uses, so the two tables compare.
const SELF_INTERSECTION_AT: u32 = 33;

struct Row {
    field: &'static str,
    rule: &'static str,
    samples: u32,
    median_ms: f64,
    vertices: usize,
    triangles: usize,
    non_manifold_edges: u64,
    self_intersections_per_1k: Option<f64>,
    hausdorff: Option<f64>,
}

fn main() {
    if !std::env::args().any(|arg| arg == "--bench") {
        return;
    }

    println!("ablation — one algorithm, two vertex rules, f64, one process\n");
    println!(
        "{:<16} {:<10} {:>7} {:>11} {:>9} {:>10} {:>6} {:>8} {:>11}",
        "field",
        "rule",
        "samples",
        "median_ms",
        "vertices",
        "triangles",
        "nm",
        "si/1k",
        "hausdorff"
    );

    let mut rows = Vec::new();
    isomesh::for_each_reference_field!(f64, |name, field| {
        for samples in RESOLUTIONS {
            // The two arms. Written out rather than looped because they are two
            // *types*: that is the point of the seam, and a loop would need a
            // runtime value to iterate, which is the thing being avoided.
            if let Some(row) = measure(
                name,
                &field,
                "qef",
                &mut DualContouring::<Scalar>::new(),
                samples,
            ) {
                print_row(&row);
                rows.push(row);
            }
            if let Some(row) = measure(
                name,
                &field,
                "centroid",
                &mut DualContouring::<Scalar, Centroid>::with_rule(Centroid),
                samples,
            ) {
                print_row(&row);
                rows.push(row);
            }
            // X-004's arm. Behind the feature because the rule is, and absent
            // rather than stubbed when it is off — a row of zeros would be worse
            // than a missing row, because it would average into the summary.
            #[cfg(feature = "experimental")]
            if let Some(row) = measure(
                name,
                &field,
                "pq_scaled",
                &mut DualContouring::<Scalar, ProbabilisticQuadric>::with_rule(
                    ProbabilisticQuadric::default(),
                ),
                samples,
            ) {
                print_row(&row);
                rows.push(row);
            }
        }
    });

    let path = write_csv(&rows);
    println!("\nwrote {}", path.display());
    report(&rows);
}

fn measure<F, V>(
    field_name: &'static str,
    field: &F,
    rule: &'static str,
    mesher: &mut DualContouring<Scalar, V>,
    samples: u32,
) -> Option<Row>
where
    F: Sdf<Scalar = Scalar> + ReferenceField,
    V: VertexRule<Scalar>,
{
    let (shape, origin, h) = common::grid(field, samples);

    let mut mesh = MeshBuffer::<Scalar>::new();
    let mut times = Vec::with_capacity((WARMUP_RUNS + TIMED_RUNS) as usize);
    for _ in 0..(WARMUP_RUNS + TIMED_RUNS) {
        mesh.reset();
        let start = Instant::now();
        mesher
            .extract(field, &shape, origin, h, &mut mesh)
            .expect("extraction");
        times.push(start.elapsed());
        black_box(mesh.triangle_count());
    }
    if mesh.triangle_count() == 0 {
        return None;
    }
    let mut timed: Vec<f64> = times[WARMUP_RUNS as usize..]
        .iter()
        .map(|d| d.as_secs_f64() * 1000.0)
        .collect();
    timed.sort_by(f64::total_cmp);
    let median_ms = timed[timed.len() / 2];

    let cfg = ValidateConfig::from_cell_size(h).ok()?;
    let report = validate_indexed(&mesh.positions, &mesh.indices, &cfg);

    let self_intersections_per_1k = if samples == SELF_INTERSECTION_AT {
        self_intersections(&mesh.positions, &mesh.indices, h)
            .ok()
            .map(|r| r.per_thousand_triangles())
    } else {
        None
    };

    // Only where the field's own metadata says the quantity is a distance.
    let hausdorff = if field.is_exact_distance() {
        AccuracyConfig::from_cell_size(h)
            .ok()
            .and_then(|acfg| {
                accuracy(&mesh.positions, &mesh.indices, field, &shape, origin, &acfg).ok()
            })
            .map(|a| a.symmetric_hausdorff())
    } else {
        None
    };

    Some(Row {
        field: field_name,
        rule,
        samples,
        median_ms,
        vertices: mesh.vertex_count(),
        triangles: mesh.triangle_count(),
        non_manifold_edges: report.non_manifold_edges,
        self_intersections_per_1k,
        hausdorff,
    })
}

fn print_row(row: &Row) {
    let si = row
        .self_intersections_per_1k
        .map_or_else(|| "       -".to_string(), |v| format!("{v:8.3}"));
    let hausdorff = row
        .hausdorff
        .map_or_else(|| "          -".to_string(), |v| format!("{v:11.3e}"));
    println!(
        "{:<16} {:<10} {:>7} {:>11.3} {:>9} {:>10} {:>6} {si} {hausdorff}",
        row.field,
        row.rule,
        row.samples,
        row.median_ms,
        row.vertices,
        row.triangles,
        row.non_manifold_edges,
    );
}

fn write_csv(rows: &[Row]) -> PathBuf {
    let mut csv = String::from(
        "field,rule,samples,median_ms,vertices,triangles,non_manifold_edges,\
         self_intersections_per_1k,hausdorff\n",
    );
    for row in rows {
        let si = row
            .self_intersections_per_1k
            .map_or_else(String::new, |v| format!("{v}"));
        let hausdorff = row.hausdorff.map_or_else(String::new, |v| format!("{v:e}"));
        let _ = writeln!(
            csv,
            "{},{},{},{:.6},{},{},{},{si},{hausdorff}",
            row.field,
            row.rule,
            row.samples,
            row.median_ms,
            row.vertices,
            row.triangles,
            row.non_manifold_edges,
        );
    }
    let path = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../docs/measurements")
        .join("ablation.csv");
    if let Some(dir) = path.parent() {
        let _ = fs::create_dir_all(dir);
    }
    let _ = fs::write(&path, csv);
    path
}

/// The comparison the table is for, stated rather than left to the reader.
///
/// The prediction, written before the run: the QEF arm should be **more
/// accurate** on fields with sharp features — that is the entire reason A-007
/// exists — and **less well behaved** on self-intersections, because a solved
/// vertex can leave its cell where a centroid cannot. If the Hausdorff columns
/// come out equal, the seam is not swapping what it claims to.
fn report(rows: &[Row]) {
    println!("\nthe rule, isolated — same topology, same walk, same buffers\n");
    println!(
        "{:<16} {:>13} {:>13} {:>9}",
        "field", "qef hausdorff", "centroid", "ratio"
    );
    let mut compared = 0usize;
    isomesh::for_each_reference_field!(f64, |name, _field| {
        let at = |rule: &str| {
            rows.iter()
                .find(|r| r.field == name && r.rule == rule && r.samples == 65)
                .and_then(|r| r.hausdorff)
        };
        if let (Some(qef), Some(centroid)) = (at("qef"), at("centroid")) {
            println!(
                "{name:<16} {qef:>13.3e} {centroid:>13.3e} {:>9.3}",
                qef / centroid
            );
            compared += 1;
        }
    });
    if compared == 0 {
        println!("(no field reported an exact distance at 65^3)");
    }
}
