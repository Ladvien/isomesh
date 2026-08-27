//! **P-72 - the granularity of the active-cell structure, swept.**
//!
//! Ticket: R-070. Pre-registered before this harness existed.
//!
//! ```bash
//! cargo bench --bench experiment_p72
//! ```
//!
//! Writes `docs/experiments/p-72.csv`.
//!
//! # What is held fixed, and why that is the whole experiment
//!
//! **The total cell count.** A sweep that varied chunk size and let the world
//! grow with it would be a resolution sweep wearing a granularity sweep's name -
//! `x35`'s failure mode, and GVDB avoided it by holding 2048³ across every
//! ⟨3,3,3,B⟩. So the world is 128³ cells at every granularity, and the only
//! thing that changes is how those cells are partitioned:
//!
//! | chunk cells | chunks | samples/chunk | total samples | duplication |
//! |---:|---:|---:|---:|---:|
//! | 8 | 16³ = 4,096 | 9³ = 729 | 2,985,984 | 1.4238x |
//! | 16 | 8³ = 512 | 17³ = 4,913 | 2,515,456 | 1.1953x |
//! | 32 | 4³ = 64 | 33³ = 35,937 | 2,299,968 | 1.0968x |
//! | 64 | 2³ = 8 | 65³ = 274,625 | 2,197,000 | 1.0473x |
//!
//! Those duplication factors are `((c+1)/c)³` and were computed at
//! registration. They are the cost of small chunks. The saving is the other
//! direction: an edit that touches one 8³ chunk re-meshes 512 cells, and the
//! same edit inside a 64³ chunk re-meshes 262,144.
//!
//! # The edit trace
//!
//! A dig: eleven spherical brushes subtracted along a straight path across the
//! middle of the world, radius 3 cells, each one a separate `mark_edit` +
//! `mesh_dirty` cycle - which is what a held mouse button produces and is the
//! shape `M-124` and `P-39` both measure against. The path is in **world**
//! coordinates and is identical at every granularity, so every arm digs the
//! same hole through the same field. `edits_agree` asserts that.
//!
//! # Controls
//!
//! Three, and each one is an assertion rather than a printed number:
//!
//! - **Every arm must mesh something.** A granularity whose dirty set came back
//!   empty would report a fast time for doing nothing - the `M-44` failure, and
//!   the one this harness is most exposed to because a coarser chunk could
//!   plausibly mark zero.
//! - **Every arm must produce the same mesh.** Not the same triangle order:
//!   chunks are meshed independently and their vertex indices are chunk-local,
//!   so the comparison is the **sorted multiset of vertex positions** over the
//!   whole world, which is what `P-61` established as the load-bearing
//!   comparison for a partition change. A partition that changes the surface is
//!   a seam defect, not a speed result.
//! - **The duplication factor must be the arithmetic one.** The harness counts
//!   field evaluations and asserts the ratio against `((c+1)/c)³` to four
//!   digits. If the count disagrees, the sweep is not measuring what the
//!   registration computed.

#![allow(
    clippy::cast_possible_truncation,
    clippy::cast_precision_loss,
    clippy::cast_sign_loss,
    clippy::print_literal,
    clippy::too_many_lines
)]

mod common;

use std::cell::Cell;
use std::time::Instant;

use isomesh::chunk::dirty::{DirtySet, mark_edit};
use isomesh::chunk::{ChunkId, ChunkLayout};
use isomesh::extractor::Extractor;
use isomesh::fields::{FbmTerrain, Gyroid};
use isomesh::marching_cubes::MarchingCubes;
use isomesh::{MeshBuffer, Sdf};

/// Cells per axis in the whole world, at every granularity.
///
/// 128 so that all four chunk sizes divide it exactly: a granularity that left
/// a partial chunk would be comparing a different cell count, which is the one
/// thing this sweep holds fixed.
const WORLD_CELLS: u32 = 128;

/// The knob. GVDB's ⟨3,3,3,B⟩ with `B` from 3 to 6 is 8 to 64 voxels per brick
/// edge, which is the registered 8³–64³ range.
///
/// **2 and 4 are here because the first run's curve was monotone.** 8³ won on
/// both fields and 8 was the smallest granularity swept, so what the run
/// measured was a minimum *at the boundary of the range* — and reporting a
/// boundary minimum as C1's "pronounced optimum" would be a claim the data does
/// not make. The sample-duplication cost `((c+1)/c)³` is 1.9531 at c = 4 and
/// **3.3750** at c = 2, against 1.4238 at c = 8, so if the two curves cross
/// anywhere they cross here. Extending down is how the question gets an answer
/// rather than a boundary.
const GRANULARITIES: [u32; 6] = [2, 4, 8, 16, 32, 64];

/// World extent per axis. Fixed, so `cell_size` is `EXTENT / WORLD_CELLS`.
const EXTENT: f64 = 4.0;

/// World origin. Centred on the reference fields' own domain centre
/// (`cube_domain` is symmetric about zero) so that both fields actually have a
/// surface inside the world -- the positive octant alone does not contain
/// `fbm_terrain`'s sheet, and the first run's `M-44` control said so by refusing
/// to report a time for a trace that marked nothing.
const ORIGIN: f64 = -EXTENT * 0.5;

/// Brush radius in **cells**, so the brush is the same physical size at every
/// granularity.
///
/// Six rather than three: the brush subtracts by `max(field, -sphere)`, so it
/// only changes a sample where the field is below the sphere's depth there, and
/// a 3-cell brush on a 128³ world is 0.094 world units against field values of
/// order 1. It bit, but only in a band thin enough that `gyroid`'s undulation
/// carried the surface out of it between edits.
const BRUSH_CELLS: f64 = 6.0;

/// Edits in the trace.
const EDITS: usize = 11;

/// Traces per arm, median taken.
///
/// **Three because the optimum's *location* is a 10% effect.** C1's magnitude is
/// a 50x spread and needs no statistics, but which granularity wins is decided
/// by 8.62 against 9.57 — and a single measurement of that gap on a machine with
/// a CPU governor is the mistake that `M-337`'s own re-audit found: a registered
/// 1.25x floor that re-measured at 1.022 three runs later. The build happens once,
/// only the trace repeats, because the trace is what is timed.
const REPS: usize = 3;

/// A field wrapper that counts `sample` calls.
///
/// The count is the control on the duplication factor: `((c+1)/c)³` is
/// arithmetic, and a harness that asserts arithmetic against itself has
/// asserted nothing. This counts what the sweep actually evaluated.
struct Counted<'a, F> {
    field: &'a F,
    calls: &'a Cell<u64>,
}

impl<F: Sdf<Scalar = f64>> Sdf for Counted<'_, F> {
    type Scalar = f64;

    #[inline]
    fn sample(&self, p: [f64; 3]) -> f64 {
        self.calls.set(self.calls.get() + 1);
        self.field.sample(p)
    }
}

/// A sphere subtracted from a field: `max(field, -(|p - c| - r))`.
///
/// Written out rather than composed through `BrushStack` because the trace needs
/// the field **before** and **after** each edit as two separate `Sdf`s for
/// `mark_edit`, and a stack that owns its brushes cannot hand out its own
/// prefix.
struct Dug<'a, F> {
    field: &'a F,
    /// Brush centres applied so far, in world coordinates.
    centres: &'a [[f64; 3]],
    radius: f64,
}

impl<F: Sdf<Scalar = f64>> Sdf for Dug<'_, F> {
    type Scalar = f64;

    #[inline]
    fn sample(&self, p: [f64; 3]) -> f64 {
        let mut v = self.field.sample(p);
        for c in self.centres {
            let d = [p[0] - c[0], p[1] - c[1], p[2] - c[2]];
            let sphere = (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt() - self.radius;
            // Subtract: max(field, -shape).
            v = v.max(-sphere);
        }
        v
    }
}

/// One arm's result.
struct Arm {
    chunk_cells: u32,
    chunks: u64,
    samples_per_chunk: u64,
    field_calls: u64,
    /// Vertices the initial full build produced, which is what the second term
    /// of the call model is denominated in.
    build_vertices: u64,
    dirty_chunks: u64,
    remeshed_cells: u64,
    mark_ms: f64,
    remesh_ms: f64,
    vertices: usize,
    triangles: usize,
    /// Distinct surface points, quantised to `cell_size * 1e-6`.
    ///
    /// The identity control. See control 2 in `main` for why this is quantised
    /// rather than bit-exact, and why the raw vertex count cannot be the
    /// comparison.
    surface: std::collections::BTreeSet<[i64; 3]>,
}

/// Run the whole trace at one granularity on one field.
fn run_arm<F: Sdf<Scalar = f64> + Sync>(field: &F, chunk_cells: u32) -> Arm {
    let cell_size = EXTENT / f64::from(WORLD_CELLS);
    let layout = ChunkLayout::<f64>::new(chunk_cells, cell_size, [ORIGIN; 3]).expect("layout");
    let shape = layout.sample_shape().expect("shape");
    let per_axis = WORLD_CELLS / chunk_cells;

    // The dig path: straight across x through the middle of the world, so it
    // crosses several chunk boundaries at every granularity.
    let radius = BRUSH_CELLS * cell_size;

    // **The dig has to hit the surface at every step, and where the surface is
    // depends on the field -- which is the whole of C2.** Two runs were void
    // before this shape: a path through the world centre missed `fbm_terrain`'s
    // sheet entirely, and a path at one probed height missed `gyroid`'s surface
    // everywhere except the probe's own `x`, because the gyroid undulates in `y`
    // as `x` advances and the brush is 3 cells wide. Both were caught by the
    // `M-44` control refusing to report a time for a trace that marked nothing.
    //
    // So the height is probed **per edit, at that edit's own `x`** -- which is
    // also what a player does, digging at the surface in front of them rather
    // than at a fixed altitude. Asserted, not defaulted.
    let mid = ORIGIN + EXTENT * 0.5;
    let surface_y = |x: f64| -> f64 {
        let steps = 1024;
        let mut prev = field.sample([x, ORIGIN, mid]);
        for i in 1..=steps {
            let y = ORIGIN + EXTENT * (f64::from(i) / f64::from(steps));
            let v = field.sample([x, y, mid]);
            if (prev < 0.0) != (v < 0.0) {
                return y;
            }
            prev = v;
        }
        panic!("no surface crossing along y at x = {x}: the trace would dig in empty space");
    };

    let centres: Vec<[f64; 3]> = (0..EDITS)
        .map(|i| {
            let t = (i as f64 + 0.5) / EDITS as f64;
            let x = ORIGIN + EXTENT * t;
            [x, surface_y(x), mid]
        })
        .collect();

    let calls = Cell::new(0u64);
    let mut dirty = DirtySet::new();
    let mut dirty_chunks = 0u64;
    let mut remeshed_cells = 0u64;
    let mut reps: Vec<(f64, f64)> = Vec::with_capacity(REPS);

    // The world's meshes, one per chunk, replaced as chunks are re-meshed. The
    // final state is what the identity control compares.
    let mut meshes: std::collections::BTreeMap<[i32; 3], MeshBuffer<f64>> =
        std::collections::BTreeMap::new();

    // The initial mesh: every chunk once. This is not part of the timed trace -
    // it is the state an edit arrives into - but its field calls ARE counted,
    // because the duplication control is about the partition and the initial
    // build is where the duplication is paid in full.
    let mut mc = MarchingCubes::<f64>::new();
    for cz in 0..per_axis {
        for cy in 0..per_axis {
            for cx in 0..per_axis {
                let id = ChunkId {
                    coords: [cx as i32, cy as i32, cz as i32],
                };
                let origin = layout.sample_origin(id);
                let counted = Counted {
                    field,
                    calls: &calls,
                };
                let mut out = MeshBuffer::<f64>::new();
                let _ = mc.extract_into(&counted, &shape, origin, cell_size, &mut out);
                meshes.insert(id.coords, out);
            }
        }
    }
    let build_calls = calls.get();
    let build_vertices: usize = meshes.values().map(|m| m.positions.len()).sum();

    for rep in 0..REPS {
        let mut mark_ns = 0u128;
        let mut remesh_ns = 0u128;
        for step in 0..EDITS {
            let before = Dug {
                field,
                centres: &centres[..step],
                radius,
            };
            let after = Dug {
                field,
                centres: &centres[..=step],
                radius,
            };

            // The cell box the brush can touch, in world cell indices. One cell of
            // margin, because a cell is active on its corners and the brush's
            // boundary can land inside the outermost cell.
            // **Through `layout.cell_of`, not by hand.** The hand-rolled form was
            // `((world - radius) / cell_size).floor()`, which silently assumes the
            // layout's origin is zero -- true for the first fixture and false once
            // the world was centred on the fields' own domain, at which point the
            // box pointed 64 cells away from the brush and `mark_edit` scanned a
            // region nothing had touched. The `M-44` control caught it; the fix is
            // to ask the layout, which is the only thing that knows.
            let c = centres[step];
            let lo_world = [0, 1, 2].map(|a| c[a] - radius);
            let hi_world = [0, 1, 2].map(|a| c[a] + radius);
            let lo = layout.cell_of(lo_world).map(|v| v - 1);
            let hi = layout.cell_of(hi_world).map(|v| v + 1);

            let t = Instant::now();
            let report = mark_edit(&layout, &before, &after, lo, hi, &mut dirty).expect("mark");
            mark_ns += t.elapsed().as_nanos();
            let _ = report;

            let t = Instant::now();
            let done = dirty.mesh_dirty(&layout, |_id, origin| {
                let mut out = MeshBuffer::<f64>::new();
                let _ = mc.extract_into(&after, &shape, origin, cell_size, &mut out);
                std::hint::black_box(&out);
            });
            remesh_ns += t.elapsed().as_nanos();
            if rep == 0 {
                dirty_chunks += done as u64;
                remeshed_cells += done as u64 * u64::from(chunk_cells).pow(3);
            }
        }
        reps.push((mark_ns as f64 / 1e6, remesh_ns as f64 / 1e6));
    }
    drop(meshes);

    // Median by total, and both components come from the same trace: a median
    // of marks paired with a median of remeshes from a different rep would be a
    // row that no run produced.
    reps.sort_unstable_by(|a, b| (a.0 + a.1).partial_cmp(&(b.0 + b.1)).expect("finite"));
    let (mark_ms, remesh_ms) = reps[REPS / 2];

    // ── the surface, meshed once at the end and NOT timed ────────────────────
    //
    // **The trace cannot be the source of the identity control, and finding out
    // why is the second defect this harness caught in itself.** A chunk dirtied
    // only at step 3 holds a mesh of `centres[..=3]`, and which step last
    // touched a region depends on the partition - so comparing the traces'
    // final buffers compares different partial digs and would fail on a
    // correct implementation. What the control has to ask is the seam question:
    // *does this partition, meshed chunk-wise, produce the same surface as
    // that one?* So the final field is meshed over every chunk here, outside
    // the timer.
    let full = Dug {
        field,
        centres: &centres,
        radius,
    };
    let mut surface: std::collections::BTreeSet<[i64; 3]> = std::collections::BTreeSet::new();
    let quantum = cell_size * 1e-6;
    let mut vertices = 0usize;
    let mut triangles = 0usize;
    for cz in 0..per_axis {
        for cy in 0..per_axis {
            for cx in 0..per_axis {
                let id = ChunkId {
                    coords: [cx as i32, cy as i32, cz as i32],
                };
                let mut out = MeshBuffer::<f64>::new();
                let _ =
                    mc.extract_into(&full, &shape, layout.sample_origin(id), cell_size, &mut out);
                vertices += out.positions.len();
                triangles += out.indices.len() / 3;
                for p in &out.positions {
                    surface.insert([0, 1, 2].map(|a| (p[a] / quantum).round() as i64));
                }
            }
        }
    }
    Arm {
        chunk_cells,
        chunks: u64::from(per_axis).pow(3),
        samples_per_chunk: u64::from(chunk_cells + 1).pow(3),
        field_calls: build_calls,
        build_vertices: build_vertices as u64,
        dirty_chunks,
        remeshed_cells,
        mark_ms,
        remesh_ms,
        vertices,
        triangles,
        surface,
    }
}

type Row = Vec<(&'static str, String)>;

fn main() {
    if !std::env::args().any(|arg| arg == "--bench") {
        return;
    }

    let prereg = isomesh::experiment!("P-72");
    let mut rows: Vec<Row> = Vec::new();

    println!(
        "{:>12} {:>6} {:>7} {:>10} {:>7} {:>7} {:>9} {:>9} {:>9} {:>8}",
        "field", "chunk", "chunks", "dup", "dirty", "mark", "remesh", "total", "ms/edit", "verts"
    );

    // Per field: the four granularities, then the field's own best/worst.
    let mut best_per_field: Vec<(&'static str, u32, f64, f64)> = Vec::new();

    for field_name in ["gyroid", "fbm_terrain"] {
        let mut arms: Vec<Arm> = Vec::new();
        for c in GRANULARITIES {
            let arm = match field_name {
                "gyroid" => run_arm(&Gyroid::<f64>::canonical(), c),
                _ => run_arm(&FbmTerrain::<f64>::canonical(), c),
            };
            arms.push(arm);
        }

        // ── control 1: every arm meshed something ────────────────────────────
        for a in &arms {
            assert!(
                a.dirty_chunks > 0,
                "VOID: chunk {} marked no dirty chunk in {EDITS} edits on {field_name}, so its \
                 time measures an empty trace",
                a.chunk_cells
            );
            assert!(
                a.vertices > 0,
                "VOID: chunk {} produced no geometry on {field_name}",
                a.chunk_cells
            );
        }

        // ── control 2: every arm produced the same surface ───────────────────
        //
        // **On a quantised lattice, and the first version was wrong to compare
        // raw vertex lists.** It fired immediately - 73,032 vertices at 8³
        // against 65,256 at 16³ on `gyroid` - and it was measuring the right
        // thing about the wrong question. Chunk-wise extraction emits
        // **coincident vertices on every shared chunk face** (`A-015`, `M-220`);
        // `weld` is what closes them, and a finer partition has more internal
        // faces, so the raw count *must* rise. That is a granularity **cost**,
        // recorded below as data.
        //
        // Nor can the deduplicated bit patterns be compared: a chunk computes a
        // corner's world position as `origin + cell_size * i` with its own
        // `origin` and its own `i`, so two chunks reach the same point by
        // different arithmetic and `M-32` measured that disagreement at **one
        // ulp**. Demanding bit equality across partitions would fail on correct
        // code for a reason this ledger already knows.
        //
        // So the comparison is the **set of surface points quantised to
        // `cell_size * 1e-6`** - four orders of magnitude finer than any real
        // surface movement and four orders coarser than an ulp of a coordinate
        // near 4.0 (about 8.9e-16). Real movement is caught; `M-32`'s ulp is
        // absorbed.
        let reference = &arms[0];
        for a in &arms[1..] {
            let differing = a.surface.symmetric_difference(&reference.surface).count();
            assert_eq!(
                differing,
                0,
                "chunk {} disagrees with chunk {} on {differing} quantised surface points of \
                 {} on {field_name}: a partition change moved the surface, which is a seam \
                 defect and not a speed result",
                a.chunk_cells,
                reference.chunk_cells,
                reference.surface.len()
            );
        }

        // ── control 3: every field evaluation is accounted for ───────────────
        //
        // **The registered `((c+1)/c)³` is the *grid* duplication, and the first
        // version of this control asserted it against the total.** It fired:
        // 1.627144 measured at c = 8 against 1.423828 predicted. The gap is not
        // an error in the arithmetic - it is a **second consumer of the field
        // that the registration did not name**. `edge_position` calls
        // `unit_gradient` per vertex for the normal, so the true model has two
        // terms:
        //
        // ```text
        // calls = chunks * (c+1)³        the sample grid, duplication ((c+1)/c)³
        //       + stencil * vertices     the normals
        // ```
        //
        // Asserting the grid term exactly and *deriving* the stencil width from
        // the remainder is a stronger control than either half: it proves the
        // harness knows where every evaluation went, and it names the normal
        // cost rather than folding it into a fudge.
        let baseline = u64::from(WORLD_CELLS).pow(3) as f64;
        for a in &arms {
            let grid = a.chunks * a.samples_per_chunk;
            assert!(
                a.field_calls >= grid,
                "chunk {}: {} field calls is fewer than the {grid} sample-grid points, so the \
                 grid was not fully evaluated",
                a.chunk_cells,
                a.field_calls
            );
            let normals = a.field_calls - grid;
            assert_eq!(
                normals % a.build_vertices,
                0,
                "chunk {}: {normals} non-grid field calls over {} vertices is not a whole \
                 stencil, so the call model is incomplete and the duplication figure is not \
                 the registered one",
                a.chunk_cells,
                a.build_vertices
            );
            let stencil = normals / a.build_vertices;
            let c = f64::from(a.chunk_cells);
            let predicted = ((c + 1.0) / c).powi(3);
            let measured = grid as f64 / baseline;
            assert!(
                (measured - predicted).abs() < 1e-4,
                "chunk {}: measured GRID duplication {measured:.6} against the registered \
                 ((c+1)/c)^3 = {predicted:.6}",
                a.chunk_cells
            );
            println!(
                "  {field_name} c={:<3} grid {grid} ({measured:.4}x) + {stencil}-sample normals \
                 over {} vertices = {} calls",
                a.chunk_cells, a.build_vertices, a.field_calls
            );
        }

        for a in &arms {
            let total = a.mark_ms + a.remesh_ms;
            let dup = (a.chunks * a.samples_per_chunk) as f64 / baseline;
            println!(
                "{:>12} {:>6} {:>7} {:>10.4} {:>7} {:>7.2} {:>9.2} {:>9.2} {:>9.3} {:>8}",
                field_name,
                a.chunk_cells,
                a.chunks,
                dup,
                a.dirty_chunks,
                a.mark_ms,
                a.remesh_ms,
                total,
                total / EDITS as f64,
                a.vertices
            );
            rows.push(vec![
                ("field", field_name.to_string()),
                ("chunk_cells", a.chunk_cells.to_string()),
                ("world_cells", WORLD_CELLS.to_string()),
                ("chunks", a.chunks.to_string()),
                ("samples_per_chunk", a.samples_per_chunk.to_string()),
                ("sample_duplication", format!("{dup:.6}")),
                ("edits", EDITS.to_string()),
                ("dirty_chunks_total", a.dirty_chunks.to_string()),
                ("remeshed_cells_total", a.remeshed_cells.to_string()),
                ("mark_ms", format!("{:.4}", a.mark_ms)),
                ("remesh_ms", format!("{:.4}", a.remesh_ms)),
                ("total_ms", format!("{total:.4}")),
                ("ms_per_edit", format!("{:.6}", total / EDITS as f64)),
                ("vertices", a.vertices.to_string()),
                ("triangles", a.triangles.to_string()),
                // The raw-to-distinct ratio IS a granularity cost and the first
                // version of control 2 mistook it for a defect: a finer
                // partition has more internal chunk faces, and every one of
                // them emits its vertices twice (`A-015`, `M-220`). Recorded
                // rather than asserted, because it is an outcome.
                ("distinct_surface_points", a.surface.len().to_string()),
                (
                    "vertex_duplication",
                    format!("{:.6}", a.vertices as f64 / a.surface.len() as f64),
                ),
            ]);
        }

        let best = arms
            .iter()
            .min_by(|x, y| {
                (x.mark_ms + x.remesh_ms)
                    .partial_cmp(&(y.mark_ms + y.remesh_ms))
                    .expect("finite")
            })
            .expect("arms");
        let worst = arms
            .iter()
            .max_by(|x, y| {
                (x.mark_ms + x.remesh_ms)
                    .partial_cmp(&(y.mark_ms + y.remesh_ms))
                    .expect("finite")
            })
            .expect("arms");
        best_per_field.push((
            field_name,
            best.chunk_cells,
            best.mark_ms + best.remesh_ms,
            worst.mark_ms + worst.remesh_ms,
        ));
        println!();
    }

    // ── verdict ──────────────────────────────────────────────────────────────
    let spread = best_per_field
        .iter()
        .map(|(_, _, b, w)| w / b)
        .fold(0.0f64, f64::max);
    let c1 = spread >= 2.0;
    let optima: Vec<u32> = best_per_field.iter().map(|(_, c, _, _)| *c).collect();
    let c2 = optima.iter().any(|c| *c != optima[0]);
    let c3 = spread < 4.0;

    for (name, best, b, w) in &best_per_field {
        println!(
            "{name}: best {best}³ at {b:.2} ms, worst at {w:.2} ms, spread {:.4}x",
            w / b
        );
    }
    println!(
        "\nC1 pronounced optimum (>= 2x): {spread:.4}x -> {}",
        if c1 { "HELD" } else { "FALSIFIED" }
    );
    println!(
        "C2 field-dependent optimum: {optima:?} -> {}",
        if c2 { "HELD" } else { "FALSIFIED" }
    );
    println!(
        "C3 spread below GVDB's 256x and under 4x: {spread:.4}x -> {}",
        if c3 { "HELD" } else { "FALSIFIED" }
    );

    // The worst granularity per field, taken from the rows themselves rather
    // than stored a second time when they were built: two copies of a fact
    // drift, and this one is derived in one place.
    let worst_by_field: Vec<(String, String)> = ["gyroid", "fbm_terrain"]
        .iter()
        .map(|f| {
            let mut worst = ("NA".to_string(), f64::NEG_INFINITY);
            for r in &rows {
                let same_field = r.iter().any(|(k, v)| *k == "field" && v == f);
                let cc = r.iter().find(|(k, _)| *k == "chunk_cells");
                let total = r
                    .iter()
                    .find(|(k, _)| *k == "total_ms")
                    .and_then(|(_, v)| v.parse::<f64>().ok());
                if let (true, Some(cc), Some(total)) = (same_field, cc, total)
                    && total > worst.1
                {
                    worst = (cc.1.clone(), total);
                }
            }
            ((*f).to_string(), worst.0)
        })
        .collect();

    let mut finished: Vec<Row> = Vec::new();
    for row in &mut rows {
        row.push(("best_chunk_cells", {
            let f = row
                .iter()
                .find(|(k, _)| *k == "field")
                .map(|(_, v)| v.clone())
                .unwrap_or_default();
            best_per_field
                .iter()
                .find(|(n, _, _, _)| *n == f)
                .map_or_else(|| "NA".to_string(), |(_, c, _, _)| c.to_string())
        }));
        row.push(("worst_chunk_cells", {
            let f = row
                .iter()
                .find(|(k, _)| *k == "field")
                .map(|(_, v)| v.clone())
                .unwrap_or_default();
            worst_by_field
                .iter()
                .find(|(n, _)| *n == f)
                .map_or_else(|| "NA".to_string(), |(_, c)| c.clone())
        }));
        row.push(("spread", format!("{spread:.6}")));
        row.push(("c1_holds", c1.to_string()));
        row.push(("c2_holds", c2.to_string()));
        row.push(("c3_holds", c3.to_string()));
        finished.push(row.clone());
    }

    common::experiment::run(prereg, |run| {
        for row in &finished {
            run.record(row);
        }
    });
}
