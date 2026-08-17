//! **What a local field edit actually costs Marching Cubes, counted.**
//!
//! Ticket: R-020. **Exploratory. Nothing is registered against this run** — this
//! phase has twice registered a clause whose arithmetic died on contact (P-23
//! clause 3, P-24 entirely), and once got it right by measuring first (M-313).
//! This measures first.
//!
//! ```bash
//! cargo bench --bench edit_trace
//! ```
//!
//! Writes `docs/measurements/edit_trace.csv`.
//!
//! # The two questions, and they have different answers
//!
//! R-020 asks whether re-meshing after a local edit has **computation distance
//! `O(|edit|)`** — whether the recorded trace changes in proportion to the cells
//! touched rather than to the grid. For Marching Cubes the trace splits in two,
//! and the split is the point:
//!
//! - **The per-cell work.** A cell's case index is a function of its eight
//!   corners and nothing else, so an edit to `k` samples can change the case of
//!   at most the cells incident to them — **at most `8k`**, since each sample is
//!   a corner of at most eight cells. That is arithmetic, and this run checks it.
//! - **The output.** Vertices are appended in scan order and indices refer to
//!   positions in that buffer, so a cell that emits a different *number* of
//!   triangles shifts **every index after it**. A sequential counter is the
//!   classic instability in Acar's sense, and nothing about the field's locality
//!   protects against it.
//!
//! **So "the computation is edit-proportional" and "the output is
//! edit-proportional" are different claims and this measures both**, on the same
//! edit, at three grid sizes with the brush held fixed.

mod common;

use std::fmt::Write as _;

use isomesh::construct::SampledField;
use isomesh::extractor::Extractor;
use isomesh::marching_cubes::MarchingCubes;
use isomesh::{MeshBuffer, RuntimeShape3, Shape3};

/// `cube::is_inside`, which is private. Negative is inside the solid, and a
/// sample of exactly zero is outside — the convention every extractor here uses.
fn is_inside(v: f64) -> bool {
    v < 0.0
}

/// Samples per axis. The brush is the same at all three.
const RESOLUTIONS: [u32; 3] = [33, 65, 129];

/// Brush radius, in samples. Fixed across resolutions — that is the experiment.
const BRUSH_RADIUS: f64 = 5.0;

/// Corner offsets of a cell, in `cube.rs`'s numbering order.
const CORNERS: [[u32; 3]; 8] = [
    [0, 0, 0],
    [1, 0, 0],
    [0, 1, 0],
    [1, 1, 0],
    [0, 0, 1],
    [1, 0, 1],
    [0, 1, 1],
    [1, 1, 1],
];

/// The eight-bit sign pattern of a cell, which is exactly what the case table
/// is indexed by — so two cells with the same mask emit the same topology.
fn case_index(values: &[f64], size: [u32; 3], cell: [u32; 3]) -> u8 {
    let mut mask = 0u8;
    for (bit, off) in CORNERS.iter().enumerate() {
        let p = [cell[0] + off[0], cell[1] + off[1], cell[2] + off[2]];
        let i =
            (p[2] as usize * size[1] as usize + p[1] as usize) * size[0] as usize + p[0] as usize;
        if values.get(i).copied().is_some_and(is_inside) {
            mask |= 1 << bit;
        }
    }
    mask
}

fn mesh(values: &[f64], shape: &RuntimeShape3, h: f64) -> MeshBuffer<f64> {
    let mut out = MeshBuffer::<f64>::new();
    let Ok(field) = SampledField::new(values, shape, [0.0; 3], h) else {
        return out;
    };
    let _ = MarchingCubes::<f64>::new().extract_into(&field, shape, [0.0; 3], h, &mut out);
    out
}

fn main() {
    if !std::env::args().any(|arg| arg == "--bench") {
        return;
    }

    println!(
        "{:>5} {:>10} {:>8} {:>7} {:>7} {:>9} {:>10} {:>11} {:>12}",
        "n",
        "cells",
        "dirty s",
        "dirty c",
        "case ch",
        "verts",
        "buf moved",
        "geom moved",
        "first moved"
    );
    let mut csv = String::from(
        "samples_per_axis,cells,dirty_samples,dirty_cells,case_changed,eight_k,\
         vertices_before,vertices_after,buffer_moved,geometric_moved,first_moved_index\n",
    );

    for n in RESOLUTIONS {
        let Ok(shape) = RuntimeShape3::new([n; 3]) else {
            continue;
        };
        let count = shape.element_count();
        let h = 4.0 / f64::from(n - 1);

        // A sphere, so there is a surface to disturb. Radius is a third of the
        // domain, well inside it, and the same shape at every resolution.
        let centre = [2.0_f64; 3];
        let radius = 1.2_f64;
        let mut before = Vec::with_capacity(count);
        for z in 0..n {
            for y in 0..n {
                for x in 0..n {
                    let p = [f64::from(x) * h, f64::from(y) * h, f64::from(z) * h];
                    let d = [p[0] - centre[0], p[1] - centre[1], p[2] - centre[2]];
                    before.push((d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt() - radius);
                }
            }
        }

        // The edit: carve a spherical brush at the sphere's equator, so it
        // straddles the surface and actually changes topology rather than
        // rewriting values nobody reads.
        let brush_centre = [n / 2 + (radius / h) as u32, n / 2, n / 2];
        let mut after = before.clone();
        let mut dirty_samples = 0u64;
        let r = BRUSH_RADIUS.ceil() as i64;
        for dz in -r..=r {
            for dy in -r..=r {
                for dx in -r..=r {
                    if (dx * dx + dy * dy + dz * dz) as f64 > BRUSH_RADIUS * BRUSH_RADIUS {
                        continue;
                    }
                    let p = [
                        i64::from(brush_centre[0]) + dx,
                        i64::from(brush_centre[1]) + dy,
                        i64::from(brush_centre[2]) + dz,
                    ];
                    if p.iter().any(|c| *c < 0 || *c >= i64::from(n)) {
                        continue;
                    }
                    let i =
                        (p[2] as usize * n as usize + p[1] as usize) * n as usize + p[0] as usize;
                    if let Some(v) = after.get_mut(i) {
                        // Carve: force it outside.
                        if *v < 1.0 {
                            *v = 1.0;
                            dirty_samples += 1;
                        }
                    }
                }
            }
        }

        // ── the computation side: which cells can possibly re-run ────────────
        let cells = [n - 1, n - 1, n - 1];
        let mut dirty_cells = 0u64;
        let mut case_changed = 0u64;
        for z in 0..cells[2] {
            for y in 0..cells[1] {
                for x in 0..cells[0] {
                    let cell = [x, y, z];
                    let mut touched = false;
                    for off in &CORNERS {
                        let p = [cell[0] + off[0], cell[1] + off[1], cell[2] + off[2]];
                        let i = (p[2] as usize * n as usize + p[1] as usize) * n as usize
                            + p[0] as usize;
                        if before.get(i) != after.get(i) {
                            touched = true;
                            break;
                        }
                    }
                    if touched {
                        dirty_cells += 1;
                        if case_index(&before, [n; 3], cell) != case_index(&after, [n; 3], cell) {
                            case_changed += 1;
                        }
                    }
                }
            }
        }

        // ── the output side: how much of the buffer actually moved ───────────
        let a = mesh(&before, &shape, h);
        let b = mesh(&after, &shape, h);
        let mut moved = 0u64;
        let mut first_moved = -1i64;
        let shared = a.positions.len().min(b.positions.len());
        for i in 0..shared {
            if a.positions.get(i) != b.positions.get(i) {
                moved += 1;
                if first_moved < 0 {
                    first_moved = i as i64;
                }
            }
        }
        moved += (a.positions.len() as i64 - b.positions.len() as i64).unsigned_abs();

        // **The control, and without it neither number means anything (E-208).**
        // The diff above compares by *index*, so a buffer holding the same
        // surface in a different order reads as entirely moved. This compares
        // the two as *sets of positions*: it counts vertices that genuinely
        // appeared or vanished, which is the surface changing rather than the
        // indices shifting. If the two columns agree, the sequential counter is
        // innocent and the surface really did move that much.
        let key = |p: &[f64; 3]| [p[0].to_bits(), p[1].to_bits(), p[2].to_bits()];
        let mut set_a: Vec<[u64; 3]> = a.positions.iter().map(key).collect();
        let mut set_b: Vec<[u64; 3]> = b.positions.iter().map(key).collect();
        set_a.sort_unstable();
        set_b.sort_unstable();
        let (mut i, mut j, mut only_a, mut only_b) = (0usize, 0usize, 0u64, 0u64);
        while i < set_a.len() && j < set_b.len() {
            match set_a[i].cmp(&set_b[j]) {
                core::cmp::Ordering::Less => {
                    only_a += 1;
                    i += 1;
                }
                core::cmp::Ordering::Greater => {
                    only_b += 1;
                    j += 1;
                }
                core::cmp::Ordering::Equal => {
                    i += 1;
                    j += 1;
                }
            }
        }
        only_a += (set_a.len() - i) as u64;
        only_b += (set_b.len() - j) as u64;
        let geometric = only_a + only_b;

        let total_cells = u64::from(cells[0]) * u64::from(cells[1]) * u64::from(cells[2]);
        println!(
            "{n:>5} {total_cells:>10} {dirty_samples:>8} {dirty_cells:>7} {case_changed:>7} \
             {:>9} {moved:>10} {geometric:>11} {first_moved:>12}",
            a.positions.len()
        );
        let _ = writeln!(
            csv,
            "{n},{total_cells},{dirty_samples},{dirty_cells},{case_changed},{},{},{},{moved},{geometric},{first_moved}",
            8 * dirty_samples,
            a.positions.len(),
            b.positions.len()
        );
    }

    let path = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../docs/measurements/edit_trace.csv");
    match std::fs::write(&path, &csv) {
        Ok(()) => println!("\nwrote {}", path.display()),
        Err(e) => println!("\n::error:: {}: {e}", path.display()),
    }
}
