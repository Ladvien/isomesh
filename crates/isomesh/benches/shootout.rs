//! M-001 — five extractors, identical fields and grids, one process, one run.
//!
//! ```bash
//! cargo bench --bench shootout
//! ```
//!
//! Writes `docs/measurements/shootout.csv`.
//!
//! # Why this is worth doing at all
//!
//! **The comparison does not exist in the literature for post-2020 hardware.**
//! The published figures this project has collected come from different papers,
//! different decades, different machines, different fields and different
//! definitions of what counts as a triangle — and several of them have already
//! failed verification here (`✗1`, `✗14`, `M-51`). One process, one run, one
//! grid removes every one of those degrees of freedom.
//!
//! # What is measured, and what is deliberately not
//!
//! Per (field, algorithm, resolution): wall-clock milliseconds, vertices,
//! triangles, non-manifold edges, self-intersections per thousand triangles, and
//! symmetric Hausdorff distance.
//!
//! Two columns are conditional, and the conditions are the fields' own metadata
//! rather than this file's opinion:
//!
//! - **Hausdorff is only reported where `is_exact_distance()`.** `gyroid` and
//!   `fbm_terrain` are not distance fields, so the quantity T-003 computes
//!   against them is not a distance and printing it would invent a number.
//! - **Self-intersection is only run at the smaller grid.** It is a
//!   broadphase-accelerated all-pairs test and the largest meshes here carry
//!   ~90,000 triangles; the ticket wants λ per algorithm, which the smaller grid
//!   already answers.
//!
//! # Timing, and why the numbers here are not `resolution_sweep`'s
//!
//! The sweep exists to fit `t = a + b·n³` and reports a median of many runs at
//! nine resolutions on one field. This runs two resolutions across seven fields
//! and five algorithms, so it takes fewer samples per point on purpose — its job
//! is the *comparison*, and a column that is 5% noisy does not change which
//! algorithm is fastest. Where the two disagree beyond that, the sweep is the
//! timing authority and this is the topology authority.
//!
//! # Reading the triangle-count column
//!
//! Marching Tetrahedra's cost against Marching Cubes is **not one number**. It
//! is `4.0` where the surface normal lies in one octant and `2.0` where it
//! changes sign, averaging `2.992` over the sphere — M-52. So `box_exact` reads
//! near 3.9 and `sphere` near 3.0, and quoting either alone would be
//! misleading. The CSV carries every field so the spread is visible rather than
//! averaged away.

mod common;

use std::fmt::Write as _;
use std::fs;
use std::hint::black_box;
use std::path::PathBuf;
use std::time::Instant;

use isomesh::dual_contouring::DualContouring;
use isomesh::fields::ReferenceField;
use isomesh::manifold_dual_contouring::ManifoldDualContouring;
use isomesh::marching_cubes::{FaceAmbiguity, MarchingCubes};
use isomesh::marching_tetrahedra::MarchingTetrahedra;
use isomesh::surface_nets::SurfaceNets;
use isomesh::validate::{
    AccuracyConfig, ValidateConfig, accuracy, self_intersections, validate_indexed,
};
use isomesh::{MeshBuffer, RuntimeShape3, Sdf};

type Scalar = f64;

/// The two grids. Coarse enough that the accuracy pass is affordable on seven
/// fields, fine enough that the surface is resolved on all of them.
const RESOLUTIONS: [u32; 2] = [33, 65];

/// Self-intersection counting runs at this resolution only. See the module docs.
const SELF_INTERSECTION_AT: u32 = 33;

const WARMUP_RUNS: u32 = 1;
const TIMED_RUNS: u32 = 3;

/// Every algorithm in the comparison.
///
/// `marching_cubes+decider` is a *row*, not a fifth algorithm: it is
/// `MarchingCubes` with A-002's asymptotic decider switched on. It earns a row
/// because the question "does the decider cost anything" is exactly the kind of
/// thing this table exists to answer, and on five of the seven fields it is
/// measurably identical (M-40), which is itself the answer.
#[derive(Clone, Copy, PartialEq, Eq)]
enum Algorithm {
    MarchingCubes,
    MarchingCubesDecider,
    MarchingTetrahedra,
    SurfaceNets,
    DualContouring,
    ManifoldDualContouring,
}

impl Algorithm {
    const ALL: [Self; 6] = [
        Self::MarchingCubes,
        Self::MarchingCubesDecider,
        Self::MarchingTetrahedra,
        Self::SurfaceNets,
        Self::DualContouring,
        Self::ManifoldDualContouring,
    ];

    fn name(self) -> &'static str {
        match self {
            Self::MarchingCubes => "marching_cubes",
            Self::MarchingCubesDecider => "marching_cubes+decider",
            Self::MarchingTetrahedra => "marching_tetrahedra",
            Self::SurfaceNets => "surface_nets",
            Self::DualContouring => "dual_contouring",
            Self::ManifoldDualContouring => "manifold_dual_contouring",
        }
    }

    /// Extract into `out`, reusing it. Every arm takes the same reused-buffer
    /// path a real workload runs, per rule 6.
    fn extract<F: Sdf<Scalar = Scalar>>(
        self,
        field: &F,
        shape: &RuntimeShape3,
        origin: [Scalar; 3],
        cell_size: Scalar,
        out: &mut MeshBuffer<Scalar>,
    ) {
        match self {
            Self::MarchingCubes => MarchingCubes::<Scalar>::new()
                .extract(field, shape, origin, cell_size, out)
                .expect("extraction"),
            Self::MarchingCubesDecider => {
                let mut mesher = MarchingCubes::<Scalar>::new();
                mesher.set_face_ambiguity(FaceAmbiguity::AsymptoticDecider);
                mesher
                    .extract(field, shape, origin, cell_size, out)
                    .expect("extraction");
            }
            Self::MarchingTetrahedra => MarchingTetrahedra::<Scalar>::new()
                .extract(field, shape, origin, cell_size, out)
                .expect("extraction"),
            Self::SurfaceNets => SurfaceNets::<Scalar>::new()
                .extract(field, shape, origin, cell_size, out)
                .expect("extraction"),
            Self::DualContouring => DualContouring::<Scalar>::new()
                .extract(field, shape, origin, cell_size, out)
                .expect("extraction"),
            Self::ManifoldDualContouring => ManifoldDualContouring::<Scalar>::new()
                .extract(field, shape, origin, cell_size, out)
                .expect("extraction"),
        }
    }
}

/// One row of the table.
struct Row {
    field: &'static str,
    algorithm: &'static str,
    samples: u32,
    median_ms: f64,
    vertices: usize,
    triangles: usize,
    non_manifold_edges: u64,
    /// `None` where the pass was not run at this resolution.
    self_intersections_per_1k: Option<f64>,
    /// `None` where the field is not an exact distance field.
    hausdorff: Option<f64>,
}

fn main() {
    // Same guard as `resolution_sweep`: `cargo test --all-targets` re-selects
    // bench targets even with `test = false`, and a debug-build shootout would
    // take minutes and overwrite the committed CSV with numbers from whichever
    // machine ran the tests.
    if !std::env::args().any(|arg| arg == "--bench") {
        return;
    }

    println!("shootout — seven reference fields, five algorithms, f64, one process\n");
    let mut rows = Vec::new();
    isomesh::for_each_reference_field!(f64, |name, field| {
        for samples in RESOLUTIONS {
            for algorithm in Algorithm::ALL {
                if let Some(row) = measure(name, &field, algorithm, samples) {
                    print_row(&row);
                    rows.push(row);
                }
            }
        }
    });

    let path = write_csv(&rows);
    println!("\nwrote {}", path.display());
    report(&rows);
}

fn measure<F: Sdf<Scalar = Scalar> + ReferenceField>(
    field_name: &'static str,
    field: &F,
    algorithm: Algorithm,
    samples: u32,
) -> Option<Row> {
    let (shape, origin, h) = common::grid(field, samples);

    let mut mesh = MeshBuffer::<Scalar>::new();
    let mut times = Vec::with_capacity((WARMUP_RUNS + TIMED_RUNS) as usize);
    for _ in 0..(WARMUP_RUNS + TIMED_RUNS) {
        mesh.reset();
        let start = Instant::now();
        algorithm.extract(field, &shape, origin, h, &mut mesh);
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
        algorithm: algorithm.name(),
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
    let lambda = row
        .self_intersections_per_1k
        .map_or_else(|| "     -".to_string(), |v| format!("{v:6.3}"));
    let hausdorff = row
        .hausdorff
        .map_or_else(|| "        -".to_string(), |v| format!("{v:9.3e}"));
    println!(
        "{:<14} {:<24} {:>4} {:>9.3} {:>8} {:>8} {:>5} {lambda} {hausdorff}",
        row.field,
        row.algorithm,
        row.samples,
        row.median_ms,
        row.vertices,
        row.triangles,
        row.non_manifold_edges,
    );
}

fn write_csv(rows: &[Row]) -> PathBuf {
    let mut csv = String::from(
        "field,algorithm,samples,median_ms,vertices,triangles,non_manifold_edges,\
         self_intersections_per_1k,hausdorff\n",
    );
    for row in rows {
        let lambda = row
            .self_intersections_per_1k
            .map_or_else(String::new, |v| format!("{v}"));
        let hausdorff = row.hausdorff.map_or_else(String::new, |v| format!("{v:e}"));
        let _ = writeln!(
            csv,
            "{},{},{},{:.6},{},{},{},{lambda},{hausdorff}",
            row.field,
            row.algorithm,
            row.samples,
            row.median_ms,
            row.vertices,
            row.triangles,
            row.non_manifold_edges,
        );
    }
    let path = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../docs/measurements")
        .join("shootout.csv");
    if let Some(dir) = path.parent() {
        let _ = fs::create_dir_all(dir);
    }
    let _ = fs::write(&path, csv);
    path
}

/// The comparison the table is for, stated rather than left to the reader.
fn report(rows: &[Row]) {
    let at = |field: &str, algorithm: &str, samples: u32| -> Option<&Row> {
        rows.iter()
            .find(|r| r.field == field && r.algorithm == algorithm && r.samples == samples)
    };

    println!("\n--- ratios against marching_cubes, per field, at 65^3 ---");
    println!(
        "{:<14} {:>10} {:>10} {:>10} {:>10} {:>10}",
        "field", "mc33 tris", "mt tris", "sn tris", "dc tris", "mdc tris"
    );
    isomesh::for_each_reference_field!(f64, |name, _field| {
        if let Some(mc) = at(name, "marching_cubes", 65) {
            let ratio = |alg: &str| {
                at(name, alg, 65).map_or_else(
                    || "     -".to_string(),
                    |r| format!("{:10.3}", r.triangles as f64 / mc.triangles as f64),
                )
            };
            println!(
                "{name:<14} {} {} {} {} {}",
                ratio("marching_cubes+decider"),
                ratio("marching_tetrahedra"),
                ratio("surface_nets"),
                ratio("dual_contouring"),
                ratio("manifold_dual_contouring")
            );
        }
    });

    println!("\n--- non-manifold edges, total over every field and resolution ---");
    for algorithm in Algorithm::ALL {
        let total: u64 = rows
            .iter()
            .filter(|r| r.algorithm == algorithm.name())
            .map(|r| r.non_manifold_edges)
            .sum();
        println!("{:<24} {total:>6}", algorithm.name());
    }

    println!("\n--- worst self-intersections per 1k triangles, at 33^3 ---");
    for algorithm in Algorithm::ALL {
        let worst = rows
            .iter()
            .filter(|r| r.algorithm == algorithm.name())
            .filter_map(|r| r.self_intersections_per_1k.map(|v| (v, r.field)))
            .fold(None::<(f64, &str)>, |acc, x| match acc {
                Some((best, _)) if best >= x.0 => acc,
                _ => Some(x),
            });
        match worst {
            Some((v, field)) => println!("{:<24} {v:8.3}  on {field}", algorithm.name()),
            None => println!("{:<24}        -", algorithm.name()),
        }
    }
}
