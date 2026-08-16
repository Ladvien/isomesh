//! **Three ways to build a distance field, measured against the analytic one.**
//!
//! Tickets: S-001, S-002, S-003. The accuracy-and-wall-clock comparison S-003
//! asks for, and the table a consumer needs to choose between them.
//!
//! ```bash
//! cargo bench --bench distance_construct
//! ```
//!
//! Writes `docs/measurements/distance_construct.csv`.
//!
//! # What is being compared
//!
//! - **exact** — Felzenszwalb & Huttenlocher's separable transform. Exact in the
//!   sense that it returns the true distance to the nearest opposite-signed
//!   *sample*, which is not the same as the true distance to the surface.
//! - **swept** — Zhao's eight-orthant Gauss–Seidel solve of `|∇d| = 1`.
//! - **marched** — Sethian's front, finalising the smallest tentative value at
//!   each step.
//!
//! Sweeping and marching share the same Godunov update and the same seeding, so
//! the difference between *those two* is purely the order they visit samples.
//! The exact transform differs from both in what it is measuring.
//!
//! # Only fields whose analytic value is a distance
//!
//! Error is measured against `field.sample(p)`, which is only a distance where
//! [`FieldBound::Exact`](isomesh::fields::FieldBound::Exact) says so. The gyroid
//! would produce a column of nonsense, so it is skipped rather than annotated.

use std::fmt::Write as _;
use std::fs;
use std::path::PathBuf;
use std::time::Instant;

use isomesh::construct::{
    signed_distance_field, signed_distance_field_marched, signed_distance_field_swept,
};
use isomesh::fields::ReferenceField;
use isomesh::{RuntimeShape3, Sdf, Shape3};

type Scalar = f64;

/// Grid sizes the comparison runs at. `O(N log N)` against `O(N)` only separates
/// with size, so a single resolution would answer nothing about scaling.
const RESOLUTIONS: [u32; 3] = [17, 33, 65];

/// Timed runs; the best is taken, since the fastest run is the one with the
/// least interference.
const RUNS: u32 = 3;

struct Row {
    field: &'static str,
    method: &'static str,
    samples: u32,
    ms: f64,
    worst: f64,
    worst_near: f64,
}

fn main() {
    if !std::env::args().any(|arg| arg == "--bench") {
        return;
    }

    println!("distance_construct — three constructors against the analytic field\n");
    println!(
        "{:<16} {:<9} {:>7} {:>10} {:>12} {:>12}",
        "field", "method", "samples", "ms", "worst", "worst_near"
    );

    let mut rows: Vec<Row> = Vec::new();

    isomesh::for_each_reference_field!(f64, |name, field| {
        // **`if` rather than `return`.** `for_each_reference_field!` expands an
        // inline block per field, not a closure, so a `return` here exits `main`
        // and the run stops at the first skipped field — which is exactly what it
        // did, silently and with exit code 0 (M-253).
        if field.bound().is_exact() {
            for samples in RESOLUTIONS {
                let (lo, hi) = field.domain();
                let h = (hi[0] - lo[0]) / f64::from(samples - 1);
                let shape = RuntimeShape3::new([samples; 3]).expect("valid shape");
                let size = shape.size();

                let mut grid = Vec::with_capacity(shape.element_count());
                for z in 0..size[2] {
                    for y in 0..size[1] {
                        for x in 0..size[0] {
                            grid.push(field.sample([
                                lo[0] + h * f64::from(x),
                                lo[1] + h * f64::from(y),
                                lo[2] + h * f64::from(z),
                            ]));
                        }
                    }
                }

                for (method, build) in [
                    (
                        "exact",
                        &signed_distance_field::<Scalar>
                            as &dyn Fn(
                                &[Scalar],
                                &RuntimeShape3,
                                Scalar,
                            ) -> isomesh::Result<Vec<Scalar>>,
                    ),
                    ("swept", &signed_distance_field_swept::<Scalar>),
                    ("marched", &signed_distance_field_marched::<Scalar>),
                ] {
                    let mut ms = f64::INFINITY;
                    let mut built = Vec::new();
                    for _ in 0..RUNS {
                        let start = Instant::now();
                        built = build(&grid, &shape, h).expect("construction");
                        ms = ms.min(start.elapsed().as_secs_f64() * 1e3);
                    }

                    // Two error columns, because they answer different questions: a
                    // consumer meshing the surface cares only about the near band,
                    // and a consumer sphere-tracing cares about all of it.
                    let mut worst = 0.0f64;
                    let mut worst_near = 0.0f64;
                    for z in 0..size[2] {
                        for y in 0..size[1] {
                            for x in 0..size[0] {
                                let p = [
                                    lo[0] + h * f64::from(x),
                                    lo[1] + h * f64::from(y),
                                    lo[2] + h * f64::from(z),
                                ];
                                let truth = field.sample(p);
                                let i = ((z * size[1] + y) * size[0] + x) as usize;
                                let err = (built[i] - truth).abs();
                                worst = worst.max(err);
                                if truth.abs() <= 2.0 * h {
                                    worst_near = worst_near.max(err);
                                }
                            }
                        }
                    }

                    println!(
                        "{name:<16} {method:<9} {samples:>7} {ms:>10.3} {worst:>12.5} {worst_near:>12.5}"
                    );
                    rows.push(Row {
                        field: name,
                        method,
                        samples,
                        ms,
                        worst,
                        worst_near,
                    });
                }
            }
        }
    });

    let mut csv = String::from("field,method,samples,ms,worst,worst_near\n");
    for r in &rows {
        let _ = writeln!(
            csv,
            "{},{},{},{:.6},{:.8},{:.8}",
            r.field, r.method, r.samples, r.ms, r.worst, r.worst_near
        );
    }
    let path = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../docs/measurements")
        .join("distance_construct.csv");
    if let Some(dir) = path.parent() {
        let _ = fs::create_dir_all(dir);
    }
    let _ = fs::write(&path, csv);
    println!("\nwrote {}", path.display());

    report(&rows);
}

/// The recommendation, stated rather than left to the reader.
fn report(rows: &[Row]) {
    let mean = |method: &str, pick: &dyn Fn(&Row) -> f64| -> f64 {
        let vals: Vec<f64> = rows
            .iter()
            .filter(|r| r.method == method)
            .map(pick)
            .collect();
        vals.iter().sum::<f64>() / vals.len().max(1) as f64
    };
    println!("\nmean over every field and resolution:");
    println!(
        "{:<9} {:>10} {:>12} {:>12}",
        "method", "ms", "worst", "worst_near"
    );
    for method in ["exact", "swept", "marched"] {
        println!(
            "{method:<9} {:>10.3} {:>12.5} {:>12.5}",
            mean(method, &|r| r.ms),
            mean(method, &|r| r.worst),
            mean(method, &|r| r.worst_near)
        );
    }
    // Scaling, which is the whole question between O(N) and O(N log N).
    println!("\nms at 65³ against 17³, per method:");
    for method in ["exact", "swept", "marched"] {
        let at = |n: u32| -> f64 {
            let vals: Vec<f64> = rows
                .iter()
                .filter(|r| r.method == method && r.samples == n)
                .map(|r| r.ms)
                .collect();
            vals.iter().sum::<f64>() / vals.len().max(1) as f64
        };
        let (small, large) = (at(17), at(65));
        // 65³/17³ is about 56, so a linear method should land near that.
        println!(
            "{method:<9} {small:>8.3} → {large:>8.3} ms  ({:>5.1}× for 56× the samples)",
            large / small
        );
    }
}
