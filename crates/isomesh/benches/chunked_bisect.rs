//! **What does chunking buy on the edit that beat the single grid?**
//!
//! Ticket: R-028. **Comparative, and deliberately without a `P-` id.** The bound
//! itself is not a hypothesis — a lockstep search runs inside one [`Air`] and
//! therefore cannot visit more samples than that `Air` holds, which is a property
//! of the construction and is asserted by
//! `connectivity::world::tests::a_bisect_visits_no_more_than_one_chunk`. What is
//! *not* settled by construction, and is what this measures, is the size of the
//! gap against the unchunked structure on a real world.
//!
//! ```bash
//! cargo bench --bench chunked_bisect
//! ```
//!
//! Writes `docs/measurements/chunked_bisect.csv`.
//!
//! # The edit
//!
//! M-321's bisect: one tunnel spanning the world, one voxel filled at its
//! midpoint. Both frontiers are then half the component, so lockstep stops only
//! after walking half of it — which on a single grid cost **1.1× a full
//! rebuild**, i.e. the incremental structure bought essentially nothing.
//!
//! Two arms over the same world, at several world sizes:
//!
//! - **`single`** — one [`Air`] spanning the whole world.
//! - **`chunked`** — an [`AirWorld`] of `cells`-sized chunks over the same
//!   samples.
//!
//! The question is whether `single` tracks the world while `chunked` tracks the
//! chunk, and by how much they differ at the largest size run here.

mod common;

use std::fmt::Write as _;

use isomesh::chunk::{ChunkId, ChunkLayout};
use isomesh::connectivity::{Air, AirWorld};
use isomesh::{RuntimeShape3, Shape3};

/// Cells per chunk, so a chunk holds `33³` samples — the resolution every other
/// connectivity measurement in this repo uses.
const CELLS: u32 = 32;

/// Chunks along x. The world grows; the chunk does not.
const WIDTHS: [i32; 4] = [2, 4, 8, 16];

/// The tunnel's height and depth within every chunk.
const TUNNEL: u32 = CELLS / 2;

fn main() {
    if !std::env::args().any(|arg| arg == "--bench") {
        return;
    }

    println!(
        "{:>8} {:>7} {:>12} {:>12} {:>12} {:>12} {:>10}",
        "arm", "chunks", "world air", "visited", "rebuild", "vs rebuild", "per chunk"
    );
    let mut csv = String::from(
        "arm,chunks,cells_per_chunk,world_samples,world_air,visited,\
         rebuild_visited,visited_over_rebuild\n",
    );

    let n = CELLS as usize + 1;
    let per_chunk = (n * n * n) as u64;

    for width in WIDTHS {
        let span = width as u32 * CELLS + 1;
        let pinch = span / 2;

        // One cavern spanning the world, pinched to a SINGLE air voxel at the
        // midplane. Filling that voxel severs two halves that are each about
        // half the world -- M-321's fixture, which is the shape the
        // sealed-volume mechanic actually produces, rather than a thin tunnel
        // whose halves are trivially small.
        let single = (|| {
            let shape = RuntimeShape3::new([span, CELLS + 1, CELLS + 1]).ok()?;
            let mut values = vec![1.0_f64; shape.element_count()];
            for z in 0..=CELLS {
                for y in 0..=CELLS {
                    if y == TUNNEL && z == TUNNEL {
                        continue;
                    }
                    let i = ((z as usize * (CELLS as usize + 1)) + y as usize) * span as usize
                        + pinch as usize;
                    if let Some(v) = values.get_mut(i) {
                        *v = -1.0;
                    }
                }
            }
            let (mut air, _) = Air::build(&values, &shape).ok()?;
            if air.components() != 1 {
                return None;
            }
            let f = air.fill(&[[pinch, TUNNEL, TUNNEL]], || true);
            Some((shape.element_count() as u64, air.air_samples(), f.visited))
        })();

        // The same world, one Air per chunk. The pinch sits inside the middle
        // chunk, at that chunk's own midplane.
        let mid_chunk = width / 2;
        let local_pinch = CELLS / 2;
        let chunked = (|| {
            let layout = ChunkLayout::<f64>::new(CELLS, 1.0, [0.0; 3]).ok()?;
            let mut world = AirWorld::new(layout);
            for c in 0..width {
                let mut values = vec![1.0_f64; n * n * n];
                if c == mid_chunk {
                    for z in 0..=CELLS {
                        for y in 0..=CELLS {
                            if y == TUNNEL && z == TUNNEL {
                                continue;
                            }
                            let i = ((z as usize * n) + y as usize) * n + local_pinch as usize;
                            if let Some(v) = values.get_mut(i) {
                                *v = -1.0;
                            }
                        }
                    }
                }
                world.load(ChunkId::new([c, 0, 0]), &values).ok()?;
            }
            if world.components() != 1 {
                return None;
            }
            let f = world.fill(
                ChunkId::new([mid_chunk, 0, 0]),
                &[[local_pinch, TUNNEL, TUNNEL]],
                || true,
            )?;
            Some(f.visited)
        })();

        let (Some((samples, air, single_visited)), Some(chunk_visited)) = (single, chunked) else {
            println!("::error:: width {width}: a fixture did not build");
            continue;
        };

        // A rebuild scans every sample of whatever it rebuilds: the world for
        // the single arm, one chunk for the chunked one.
        for (arm, visited, rebuild) in [
            ("single", single_visited, samples),
            ("chunked", chunk_visited, per_chunk),
        ] {
            let ratio = if rebuild == 0 {
                0.0
            } else {
                visited as f64 / rebuild as f64
            };
            println!(
                "{arm:>8} {width:>7} {air:>12} {visited:>12} {rebuild:>12} {ratio:>11.3}x {:>9.3}",
                visited as f64 / per_chunk as f64
            );
            let _ = writeln!(
                csv,
                "{arm},{width},{CELLS},{samples},{air},{visited},{rebuild},{ratio:.4}"
            );
        }
    }

    let path = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../docs/measurements/chunked_bisect.csv");
    match std::fs::write(&path, &csv) {
        Ok(()) => println!("\nwrote {}", path.display()),
        Err(e) => println!("\n::error:: {}: {e}", path.display()),
    }
}
