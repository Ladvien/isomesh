//! **A timing here is a property of the binary, not only of the code (M-281).**
//!
//! Ticket: M-001. Run with `cargo bench --bench layout_bias`; it prints and
//! asserts rather than writing a CSV, because its output is a *ratio that must
//! be one*, not a measurement to track.
//!
//! # The incident
//!
//! `benches/family` and `benches/resolution_sweep` measure Marching Cubes on
//! the same field, at the same resolutions, with the same warmup and median
//! rule — and disagreed by a **uniform 1.24–1.36×** at every resolution from
//! 16³ to 256³, on the same machine with the clock held at 4.18 GHz. A uniform
//! ratio at 16³, where the whole run is 40 µs and nothing is under memory
//! pressure, rules out cache, TLB and allocation.
//!
//! # What this file rules out, and what it leaves
//!
//! It puts **both loop shapes in one binary**: `family`'s (warmups in their own
//! loop, then timed runs, `black_box(&mesh)`) and `resolution_sweep`'s (one
//! loop of warmups plus timed runs, drained and sorted,
//! `black_box(mesh.triangle_count())`). They come out **identical** — the
//! assertion below is that their ratio is within 5% — so the loop shape is not
//! it.
//!
//! What is left is the binary. The confirming step cannot live in one file and
//! is recorded here instead of being lost: adding **one unrelated function** to
//! `resolution_sweep.rs` moved its own Marching Cubes 256³ row from
//! **152.5 ms to 130.8 ms**, a 17% change with no change to the measured code.
//! That is classic layout bias — Mytkowicz et al., *Producing Wrong Data
//! Without Doing Anything Obviously Wrong!* (`10.1145/1508284.1508275`) — and
//! the consequence for this repository is stated in M-281: **a millisecond
//! figure is comparable only against other figures from the same binary and the
//! same build.** Ratios measured side by side in one run survive; absolute
//! numbers quoted across benches do not.
//!
//! Both orders are run, so M-197's *"whichever runs second pays"* cannot be it
//! either.

mod common;

use std::hint::black_box;
use std::time::{Duration, Instant};

use isomesh::MeshBuffer;
use isomesh::extractor::Extractor;
use isomesh::fields::Sphere;
use isomesh::marching_cubes::MarchingCubes;

type Scalar = f32;
const WARMUP: usize = 2;
const TIMED: usize = 5;

/// `family`'s shape: warmups in their own loop, then timed runs.
fn family_shape(samples: u32) -> f64 {
    let field = Sphere::<Scalar>::canonical();
    let (shape, origin, h) = common::grid(&field, samples);
    let mut extractor = MarchingCubes::<Scalar>::new();
    let mut mesh = MeshBuffer::<Scalar>::new();
    for _ in 0..WARMUP {
        mesh.reset();
        extractor
            .extract_into(&field, &shape, origin, h, &mut mesh)
            .expect("extraction");
        black_box(&mesh);
    }
    let mut runs: Vec<u128> = Vec::with_capacity(TIMED);
    for _ in 0..TIMED {
        mesh.reset();
        let started = Instant::now();
        extractor
            .extract_into(&field, &shape, origin, h, &mut mesh)
            .expect("extraction");
        runs.push(started.elapsed().as_nanos());
        black_box(&mesh);
    }
    runs.sort_unstable();
    runs[TIMED / 2] as f64 / 1e6
}

/// `resolution_sweep`'s shape: one loop of warmups + timed, drain, sort.
fn sweep_shape(samples: u32) -> f64 {
    let field = Sphere::<Scalar>::canonical();
    let (shape, origin, h) = common::grid(&field, samples);
    let mut extractor = MarchingCubes::<Scalar>::new();
    let mut mesh = MeshBuffer::<Scalar>::new();
    let mut times = Vec::with_capacity(WARMUP + TIMED);
    for _ in 0..(WARMUP + TIMED) {
        mesh.reset();
        let start = Instant::now();
        extractor
            .extract_into(&field, &shape, origin, h, &mut mesh)
            .expect("extraction");
        times.push(start.elapsed());
        black_box(mesh.triangle_count());
    }
    times.drain(..WARMUP);
    times.sort_unstable();
    times[times.len() / 2].as_secs_f64() * 1e3
}

/// The sweep's shape but with `f64` instead of `f32`, in case the scalar is it.
fn sweep_shape_f64(samples: u32) -> f64 {
    let field = Sphere::<f64>::canonical();
    let (shape, origin, h) = common::grid(&field, samples);
    let mut extractor = MarchingCubes::<f64>::new();
    let mut mesh = MeshBuffer::<f64>::new();
    let mut times: Vec<Duration> = Vec::with_capacity(WARMUP + TIMED);
    for _ in 0..(WARMUP + TIMED) {
        mesh.reset();
        let start = Instant::now();
        extractor
            .extract_into(&field, &shape, origin, h, &mut mesh)
            .expect("extraction");
        times.push(start.elapsed());
        black_box(mesh.triangle_count());
    }
    times.drain(..WARMUP);
    times.sort_unstable();
    times[times.len() / 2].as_secs_f64() * 1e3
}

fn main() {
    if !std::env::args().any(|arg| arg == "--bench") {
        return;
    }
    println!(
        "{:>5} {:>12} {:>12} {:>8} {:>14}",
        "n", "family_ms", "sweep_ms", "ratio", "sweep_f64_ms"
    );
    for n in [16u32, 32, 64, 128, 256] {
        // Both orders, so "whichever runs second pays" (M-197) cannot be it.
        let a = family_shape(n);
        let b = sweep_shape(n);
        let b2 = sweep_shape(n);
        let a2 = family_shape(n);
        let f64_ms = sweep_shape_f64(n);
        let ratio = b / a;
        println!(
            "{n:>5} {a:>12.4} {b:>12.4} {ratio:>8.3} {f64_ms:>14.4}   (2nd pass: family {a2:.4}, \
             sweep {b2:.4})"
        );
        assert!(
            (0.95..=1.05).contains(&ratio),
            "the two loop shapes differ by {:.1}% at {n}³ inside one binary, which would make the \
             shape the cause after all and this file's conclusion wrong",
            (ratio - 1.0) * 100.0
        );
    }
}
