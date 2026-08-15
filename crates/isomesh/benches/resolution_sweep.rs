//! T-006's resolution sweep: fit `t = a + b·n³` and report `a` separately.
//!
//! Run with `cargo bench --bench resolution_sweep`. Writes
//! `docs/measurements/resolution_sweep.csv`.
//!
//! # Why this is not a criterion benchmark
//!
//! Criterion measures one configuration well. This is a *study* across nine of
//! them whose output is a fitted coefficient and a CSV, so it needs the raw
//! per-resolution timings rather than a rendered report. `extract.rs` keeps the
//! criterion role — noticing regressions on a fixed grid — and the two share
//! [`common::grid`] so they cannot disagree about what "64³" means.
//!
//! # The prediction, written down before the first run
//!
//! The ticket's rationale is V-6: **73% of a published 64³ figure was fixed
//! launch overhead**, hence "stop trusting single-grid numbers". That figure is
//! *GPU dispatch* overhead. There is no dispatch on this path, so the honest
//! prediction is that **`a` comes out near zero here and the warning applies to
//! Phase 6, not to the CPU path.** If it does, that is the finding, and it means
//! the rule was imported from a context that does not apply yet.
//!
//! # The misspecification the two-term model invites
//!
//! `extract` does `O(n³)` work sampling and marching, and `O(n²)` work creating
//! vertices, because the surface is two-dimensional. Fitting `a + b·n³` to data
//! that is really `a + c·n² + b·n³` pushes the `n²` term into `a`, which would
//! look exactly like fixed overhead and would not be.
//!
//! So the fit is reported twice: over the whole sweep, and over the large-`n`
//! tail alone. `n²/n³ = 1/n`, so the surface term's share shrinks as `n` grows —
//! if the two `a` values disagree materially, the two-term model is absorbing
//! the surface term and `a` must not be read as fixed cost. That check costs one
//! extra fit and is the difference between a number and an understood number.

mod common;

use std::fmt::Write as _;
use std::fs;
use std::hint::black_box;
use std::path::PathBuf;
use std::time::Instant;

use isomesh::MeshBuffer;
use isomesh::dual_contouring::DualContouring;
use isomesh::extractor::Extractor;
use isomesh::fields::Sphere;
use isomesh::marching_cubes::MarchingCubes;
use isomesh::surface_nets::SurfaceNets;

/// Samples per axis. Spans the 16³→256³ range the ticket asks for, with enough
/// points either side of the tail cut for both fits to be meaningful.
const RESOLUTIONS: [u32; 9] = [16, 24, 32, 48, 64, 96, 128, 192, 256];

/// Untimed runs first, so the fit measures steady-state re-meshing rather than
/// first-touch page faults on freshly grown scratch buffers.
const WARMUP_RUNS: u32 = 2;

/// Timed runs per resolution. The median is taken, not the mean: one scheduler
/// preemption should not move the number the fit sees.
const TIMED_RUNS: u32 = 5;

/// The tail fit uses resolutions at or above this.
const TAIL_FROM: u32 = 64;

/// `f32`, because that is what a game passes and because 256³ in `f64` would
/// hold roughly 400 MB of scratch. The precision comparison lives in
/// `extract.rs`, at a resolution where both fit comfortably in cache-friendly
/// memory.
type Scalar = f32;

struct Row {
    algorithm: &'static str,
    samples: u32,
    n_cubed: f64,
    median_ms: f64,
    vertices: usize,
    triangles: usize,
}

// **Marching Tetrahedra is deliberately not swept**, and the reason outlived
// this file's own `Extractor` trait, deleted at X-001 in favour of the crate's.
//
// Adding a fourth row would rewrite `docs/measurements/resolution_sweep.csv`,
// and that file is committed evidence: ✗14, M-19, M-20, M-21, M-22 and O-11's
// cross-machine comparison all quote exact figures from it. Re-running the
// sweep to add Marching Tetrahedra would move every one of those by measurement
// noise, so the row waits for M-001, which is the ticket that re-measures the
// whole family in one process and one run on purpose.
//
// A registry to enumerate from now exists — `isomesh::for_each_extractor!` —
// which makes the three swept here a **choice** rather than an omission. It is
// still a choice.

fn main() {
    // Cargo passes `--bench` under `cargo bench` and passes no arguments at all
    // under `cargo test`. CI runs `cargo test --workspace --all-targets`, which
    // re-selects bench targets even though the manifest sets `test = false`, so
    // without this guard a full debug-build sweep runs on every test invocation
    // -- minutes of work, and it rewrites the committed CSV with numbers from
    // whichever machine happened to run it. This is the discriminator criterion
    // uses for the same reason.
    if !std::env::args().any(|arg| arg == "--bench") {
        return;
    }

    println!("resolution sweep — sphere, f32, {TIMED_RUNS} timed runs per point\n");
    println!(
        "{:<20} {:>8} {:>12} {:>12} {:>10} {:>10}",
        "alg", "samples", "n^3", "median ms", "verts", "tris"
    );

    let mut rows = sweep("marching_cubes", MarchingCubes::<Scalar>::new);
    rows.extend(sweep("surface_nets", SurfaceNets::<Scalar>::new));
    rows.extend(sweep("dual_contouring", DualContouring::<Scalar>::new));

    let path = write_csv(&rows);
    println!("\nwrote {}", path.display());
    report(&rows);
}

/// Time one extractor across every resolution.
///
/// The extractor and its scratch live inside the per-resolution loop, so 256³'s
/// several hundred megabytes are released before the next resolution rather than
/// held for the whole sweep.
fn sweep<E: Extractor<Scalar>>(name: &'static str, make: impl Fn() -> E) -> Vec<Row> {
    let field = Sphere::<Scalar>::canonical();
    let mut rows = Vec::with_capacity(RESOLUTIONS.len());
    for n in RESOLUTIONS {
        let (shape, origin, h) = common::grid(&field, n);
        let mut extractor = make();
        let mut mesh = MeshBuffer::<Scalar>::new();
        let mut times = Vec::with_capacity((WARMUP_RUNS + TIMED_RUNS) as usize);

        for _ in 0..(WARMUP_RUNS + TIMED_RUNS) {
            // Reused buffer, per rule 6 — this is the re-mesh path a real
            // workload runs, not a cold first call.
            mesh.reset();
            let start = Instant::now();
            extractor
                .extract_into(&field, &shape, origin, h, &mut mesh)
                .expect("extraction");
            times.push(start.elapsed());
            black_box(mesh.triangle_count());
        }

        // The warm-up runs go through the identical path and are discarded here,
        // rather than being measured by different code.
        times.drain(..WARMUP_RUNS as usize);
        times.sort_unstable();

        let median_ms = times[times.len() / 2].as_secs_f64() * 1e3;
        let n_cubed = f64::from(n).powi(3);
        println!(
            "{:<20} {n:>8} {n_cubed:>12.0} {median_ms:>12.3} {:>10} {:>10}",
            name,
            mesh.vertex_count(),
            mesh.triangle_count()
        );
        rows.push(Row {
            algorithm: name,
            samples: n,
            n_cubed,
            median_ms,
            vertices: mesh.vertex_count(),
            triangles: mesh.triangle_count(),
        });
    }
    rows
}

/// Ordinary least squares for `t = a + b·x`. Returns `(a, b, r²)`.
fn fit(xs: &[f64], ts: &[f64]) -> (f64, f64, f64) {
    let n = xs.len() as f64;
    let sx: f64 = xs.iter().sum();
    let st: f64 = ts.iter().sum();
    let sxx: f64 = xs.iter().map(|x| x * x).sum();
    let sxt: f64 = xs.iter().zip(ts).map(|(x, t)| x * t).sum();
    let denom = n * sxx - sx * sx;
    let b = (n * sxt - sx * st) / denom;
    let a = (st - b * sx) / n;

    let mean = st / n;
    let ss_tot: f64 = ts.iter().map(|t| (t - mean) * (t - mean)).sum();
    let ss_res: f64 = xs
        .iter()
        .zip(ts)
        .map(|(x, t)| {
            let e = t - (a + b * x);
            e * e
        })
        .sum();
    let r2 = if ss_tot > 0.0 {
        1.0 - ss_res / ss_tot
    } else {
        1.0
    };
    (a, b, r2)
}

/// Report every algorithm actually present in `rows`, in first-seen order.
///
/// **Derived from the rows rather than listed**, and that is the fix for a live
/// defect rather than a style preference: this used to filter on a hard-coded
/// `["mc", "sn", "dc"]` while `Extractor::NAME` had been spelled out to
/// `marching_cubes` / `surface_nets` / `dual_contouring`. Every selection came
/// back empty, `fit` divided by zero, and the whole `t = a + b·n³` block printed
/// `NaN` — silently, because nothing asserts on a benchmark's stdout. A list that
/// can drift from the thing it names will eventually drift from it.
fn algorithms_in(rows: &[Row]) -> Vec<&'static str> {
    let mut seen: Vec<&'static str> = Vec::new();
    for row in rows {
        if !seen.contains(&row.algorithm) {
            seen.push(row.algorithm);
        }
    }
    seen
}

fn report(rows: &[Row]) {
    for algorithm in algorithms_in(rows) {
        let all: Vec<&Row> = rows.iter().filter(|r| r.algorithm == algorithm).collect();
        // A fit needs two points to have a slope at all; below that `fit`
        // divides by zero and prints NaN, which reads like a measurement.
        if all.len() < 2 {
            println!("\n{algorithm}:  {} point(s) — too few to fit", all.len());
            continue;
        }
        let tail: Vec<&Row> = all
            .iter()
            .copied()
            .filter(|r| r.samples >= TAIL_FROM)
            .collect();

        let xs: Vec<f64> = all.iter().map(|r| r.n_cubed).collect();
        let ts: Vec<f64> = all.iter().map(|r| r.median_ms).collect();
        let (a, b, r2) = fit(&xs, &ts);

        let xs_t: Vec<f64> = tail.iter().map(|r| r.n_cubed).collect();
        let ts_t: Vec<f64> = tail.iter().map(|r| r.median_ms).collect();
        let (a_t, b_t, r2_t) = fit(&xs_t, &ts_t);

        println!("\n{algorithm}:  t = a + b·n³");
        println!(
            "  full sweep  a = {:>10.4} ms   b = {:>8.4} ns/sample   r² = {r2:.6}",
            a,
            b * 1e6
        );
        println!(
            "  n >= {TAIL_FROM}     a = {:>10.4} ms   b = {:>8.4} ns/sample   r² = {r2_t:.6}",
            a_t,
            b_t * 1e6
        );

        println!("  marginal throughput {:.1} M samples/s", 1e3 / (b * 1e6));

        // `a` is only meaningful against the times actually measured, so it is
        // compared to both ends of the sweep rather than to the other fit's `a`.
        // Comparing the two estimates of `a` to each other is useless when both
        // are near zero: any wobble is a large relative change.
        let largest = ts.iter().copied().fold(f64::MIN, f64::max);
        let smallest = ts.iter().copied().fold(f64::MAX, f64::min);
        println!(
            "  a is {:.2}% of the largest run ({largest:.3} ms) and {:.0}% of the smallest ({smallest:.3} ms)",
            a / largest * 100.0,
            a / smallest * 100.0
        );

        // Three separate verdicts, because they answer different questions.
        if a < 0.0 {
            println!("  a < 0 is not physically possible, so `t = a + b·n³` does not");
            println!("  describe this algorithm: its cost grows faster than n³ over this range");
        } else if a.abs() < 0.01 * largest {
            println!("  a is negligible at scale: there is no meaningful fixed cost on this path");
        } else {
            println!(
                "  a is a material share of a large run — worth attributing before optimising"
            );
        }

        if a.abs() > smallest {
            println!(
                "  but a exceeds the smallest measured run, so do not extrapolate the fit below",
            );
            println!(
                "  n = {}: at that end the O(n²) surface term dominates",
                all.first().map_or(0, |r| r.samples)
            );
        }
        if r2 < 0.999 {
            println!("  r² = {r2:.6} — the two-term model is a poor description; see the CSV");
        }
    }
}

fn write_csv(rows: &[Row]) -> PathBuf {
    let dir = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../docs/measurements");
    fs::create_dir_all(&dir).expect("create docs/measurements");
    let path = dir.join("resolution_sweep.csv");

    let mut csv =
        String::from("algorithm,scalar,field,samples,cells,n_cubed,median_ms,vertices,triangles\n");
    for r in rows {
        let cells = u64::from(r.samples - 1).pow(3);
        writeln!(
            csv,
            "{},f32,sphere,{},{},{:.0},{:.6},{},{}",
            r.algorithm, r.samples, cells, r.n_cubed, r.median_ms, r.vertices, r.triangles
        )
        .expect("format csv");
    }
    fs::write(&path, csv).expect("write csv");
    path.canonicalize().unwrap_or(path)
}
