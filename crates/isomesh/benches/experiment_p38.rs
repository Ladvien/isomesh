//! **P-38 — A-024's remedy reached one of the two engines.**
//!
//! Ticket: R-037. Pre-registered in the commit before this one.
//!
//! ```bash
//! cargo bench --bench experiment_p38
//! ```
//!
//! Writes `docs/experiments/p-38.csv`.
//!
//! # What is being compared, and why the caller arranges it
//!
//! The measurement predates the fix, and it has to, or there is nothing to
//! compare against. So the two arms are arranged by the **caller** rather than
//! by the mesher, which is exactly how A-024 diagnosed the dual path (M-287):
//!
//! - `before` — the caller's own cube, `[n, n, n]`. At `n = 128` that is a
//!   512-byte row stride and a 65,536-byte plane stride.
//! - `after` — `[n | 1, n, n]`, the layout the fix produces internally. It adds
//!   one sample per row, under 1% more work, and moves both strides.
//!
//! Both are reported per sample of their own grid, so the extra row is paid for
//! rather than hidden. The `pad_z` column is the control that makes this a
//! measurement rather than a story: `[n, n, n | 1]` adds the same work and
//! touches neither stride, so it must keep the whole penalty.
//!
//! `surface_nets` is the second control. It already carries `size[0] | 1`, so
//! its `before` and `after` are the same layout and must agree.
//!
//! # Two fields, and the registered one is the harder case
//!
//! A-024 measured on a field with no surface, isolating the scaffolding. That is
//! `scaffolding_*` here — the canonical sphere sampled a long way from itself, so
//! every corner is outside and no vertex is ever created. The **registered**
//! columns use the canonical sphere on its own domain, surface and all, because
//! that is the workload and it is the case that can dilute the effect.

mod common;

use std::time::Instant;

use isomesh::fields::{ReferenceField, Sphere};
use isomesh::marching_cubes::MarchingCubes;
use isomesh::surface_nets::SurfaceNets;
use isomesh::validate::mesh_hash;
use isomesh::{MeshBuffer, RuntimeShape3};

/// Sizes whose cube is an aliasing period, and the neighbours that bracket them.
const SPIKES: [u32; 2] = [128, 256];

/// Timed repetitions. The median is reported; a mean would be dragged by the
/// first run's page faults.
const REPS: usize = 5;

#[derive(Clone, Copy, PartialEq, Eq)]
enum Which {
    MarchingCubes,
    SurfaceNets,
}

impl Which {
    const fn name(self) -> &'static str {
        match self {
            Self::MarchingCubes => "marching_cubes",
            Self::SurfaceNets => "surface_nets",
        }
    }
}

/// One extraction into a reused buffer, so allocation is not in the timing.
fn extract_once(
    which: Which,
    mc: &mut MarchingCubes<f64>,
    sn: &mut SurfaceNets<f64>,
    field: &Sphere<f64>,
    shape: &RuntimeShape3,
    origin: [f64; 3],
    cell: f64,
    out: &mut MeshBuffer<f64>,
) {
    out.reset();
    match which {
        Which::MarchingCubes => mc
            .extract(field, shape, origin, cell, out)
            .expect("extraction"),
        Which::SurfaceNets => sn
            .extract(field, shape, origin, cell, out)
            .expect("extraction"),
    }
}

/// Median nanoseconds per sample of *this* grid.
fn ns_per_sample(which: Which, size: [u32; 3], origin: [f64; 3], cell: f64) -> f64 {
    let field = Sphere::<f64>::canonical();
    let shape = RuntimeShape3::new(size).expect("valid shape");
    let samples = f64::from(size[0]) * f64::from(size[1]) * f64::from(size[2]);
    let mut mc = MarchingCubes::<f64>::new();
    let mut sn = SurfaceNets::<f64>::new();
    let mut out = MeshBuffer::<f64>::new();

    // One untimed run so every buffer is allocated and every page is resident.
    extract_once(
        which, &mut mc, &mut sn, &field, &shape, origin, cell, &mut out,
    );

    let mut runs = Vec::with_capacity(REPS);
    for _ in 0..REPS {
        let t = Instant::now();
        extract_once(
            which, &mut mc, &mut sn, &field, &shape, origin, cell, &mut out,
        );
        runs.push(t.elapsed().as_secs_f64() * 1e9 / samples);
    }
    runs.sort_by(f64::total_cmp);
    runs[runs.len() / 2]
}

/// The canonical sphere's own domain: `[-2, 2]³`, surface included.
fn surfaced(n: u32) -> ([f64; 3], f64) {
    let (lo, hi) = Sphere::<f64>::canonical().domain();
    ([lo[0], lo[1], lo[2]], (hi[0] - lo[0]) / f64::from(n - 1))
}

/// The same sphere, sampled a long way from itself. Every corner is outside, so
/// no edge is ever cut and what is left is the scaffolding.
fn surface_free(n: u32) -> ([f64; 3], f64) {
    ([10.0; 3], 4.0 / f64::from(n - 1))
}

/// Every `marching_cubes` row of the committed fixture, recomputed.
///
/// Clause three is *"the golden hashes do not move"*, and this is that check in
/// the artefact rather than in someone's memory: the same
/// [`isomesh::validate::mesh_hash`] `golden.rs` commits with, over the same
/// field-and-resolution matrix, compared against `golden_hashes.json` on disk.
fn golden_marching_cubes_unchanged() -> (usize, usize) {
    let path = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("golden_hashes.json");
    let Ok(text) = std::fs::read_to_string(&path) else {
        return (0, 0);
    };

    let mut checked = 0usize;
    let mut matched = 0usize;
    isomesh::for_each_reference_field!(f64, |name, field| {
        for samples in [17u32, 25, 33] {
            let needle = format!(
                "\"algorithm\":\"marching_cubes\",\"field\":\"{name}\",\"samples\":{samples},"
            );
            let Some(line) = text.lines().find(|l| l.contains(&needle)) else {
                continue;
            };
            let Some(committed) = line
                .split("\"hash\":\"")
                .nth(1)
                .and_then(|rest| rest.split('"').next())
                .and_then(|hex| u64::from_str_radix(hex, 16).ok())
            else {
                continue;
            };

            let (lo, hi) = field.domain();
            let cell = (hi[0] - lo[0]) / f64::from(samples - 1);
            let shape = RuntimeShape3::new([samples; 3]).expect("valid shape");
            let mut out = MeshBuffer::<f64>::new();
            MarchingCubes::<f64>::new()
                .extract(&field, &shape, lo, cell, &mut out)
                .expect("extraction");

            checked += 1;
            if mesh_hash(&out) == committed {
                matched += 1;
            }
        }
    });
    (checked, matched)
}

fn main() {
    let prereg = isomesh::experiment!("P-38");

    let (checked, matched) = golden_marching_cubes_unchanged();
    let golden_unchanged = checked > 0 && checked == matched;
    println!("golden marching_cubes rows: {matched}/{checked} match\n");

    common::experiment::run(prereg, |run| {
        for which in [Which::MarchingCubes, Which::SurfaceNets] {
            for n in SPIKES {
                let (origin, cell) = surfaced(n);
                let (free_origin, free_cell) = surface_free(n);

                let before = ns_per_sample(which, [n; 3], origin, cell);
                let after = ns_per_sample(which, [n | 1, n, n], origin, cell);
                let pad_z = ns_per_sample(which, [n, n, n | 1], origin, cell);

                let (lo_origin, lo_cell) = surfaced(n - 1);
                let (hi_origin, hi_cell) = surfaced(n + 1);
                let lo = ns_per_sample(which, [n - 1; 3], lo_origin, lo_cell);
                let hi = ns_per_sample(which, [n + 1; 3], hi_origin, hi_cell);
                let neighbour_mean = (lo + hi) / 2.0;

                let free_before = ns_per_sample(which, [n; 3], free_origin, free_cell);
                let free_after = ns_per_sample(which, [n | 1, n, n], free_origin, free_cell);
                let (free_lo_origin, free_lo_cell) = surface_free(n - 1);
                let (free_hi_origin, free_hi_cell) = surface_free(n + 1);
                let free_lo = ns_per_sample(which, [n - 1; 3], free_lo_origin, free_lo_cell);
                let free_hi = ns_per_sample(which, [n + 1; 3], free_hi_origin, free_hi_cell);
                let free_neighbour_mean = (free_lo + free_hi) / 2.0;

                println!(
                    "{:>14} {n:>4}³  surfaced {before:7.3} → {after:7.3} ns/sample \
                     (neighbours {neighbour_mean:7.3}, pad_z {pad_z:7.3})",
                    which.name()
                );
                println!(
                    "{:>14}        scaffold {free_before:7.3} → {free_after:7.3} ns/sample \
                     (neighbours {free_neighbour_mean:7.3})",
                    ""
                );

                run.record(&[
                    ("extractor", which.name().to_string()),
                    ("samples_per_axis", n.to_string()),
                    ("ns_per_sample_before", format!("{before:.4}")),
                    ("ns_per_sample_after", format!("{after:.4}")),
                    (
                        "neighbour_ratio_before",
                        format!("{:.4}", before / neighbour_mean),
                    ),
                    (
                        "neighbour_ratio_after",
                        format!("{:.4}", after / neighbour_mean),
                    ),
                    ("golden_unchanged", golden_unchanged.to_string()),
                    ("ns_neighbour_lo", format!("{lo:.4}")),
                    ("ns_neighbour_hi", format!("{hi:.4}")),
                    ("ns_pad_z_control", format!("{pad_z:.4}")),
                    ("pad_z_ratio", format!("{:.4}", pad_z / neighbour_mean)),
                    ("scaffold_ns_before", format!("{free_before:.4}")),
                    ("scaffold_ns_after", format!("{free_after:.4}")),
                    (
                        "scaffold_ratio_before",
                        format!("{:.4}", free_before / free_neighbour_mean),
                    ),
                    (
                        "scaffold_ratio_after",
                        format!("{:.4}", free_after / free_neighbour_mean),
                    ),
                    ("golden_rows_checked", checked.to_string()),
                    ("golden_rows_matched", matched.to_string()),
                ]);
            }
        }
    });
}
