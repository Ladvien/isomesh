//! **P-89 - the granularity curve extended to 1³, against `M-377`'s own model.**
//!
//! Ticket: R-089. Pre-registered before this harness existed.
//!
//! ```bash
//! cargo bench --bench experiment_p89
//! ```
//!
//! Writes `docs/experiments/p-89.csv`.
//!
//! # What this is
//!
//! `M-377` swept chunk granularity over 2³-64³ at a fixed 128³ world and found
//! an *interior* optimum at 4³ - interior only because 2³ was measured and lost.
//! The curve is still **rising** at 2³ on both fields, so 1³ should be worse
//! still, and `M-377`'s own "would be shown wrong by" says so in as many words:
//! *"an optimum below 2³ - the curve is rising at 2³ on both fields, so 1³
//! should be worse still, but it is untested."* This is that test.
//!
//! The fixture is `experiment_p72.rs`. Same 128³ world, same `EXTENT`, same
//! `ORIGIN`, same six-cell brush, same eleven-edit probed-per-edit dig path,
//! same three-trace median, same three controls. The only changes are the
//! granularity list - `[1, 2, 4]` - and the model machinery C2 needs.
//!
//! | chunk cells | chunks | samples/chunk | total samples | duplication |
//! |---:|---:|---:|---:|---:|
//! | **1** | **128³ = 2,097,152** | **2³ = 8** | **16,777,216** | **8.0000x** |
//! | 2 | 64³ = 262,144 | 3³ = 27 | 7,077,888 | 3.3750x |
//! | 4 | 32³ = 32,768 | 5³ = 125 | 4,096,000 | 1.9531x |
//!
//! # C2's model, and why it is read from the committed CSV
//!
//! `M-377`'s two-term model is a **field-evaluation count**:
//!
//! ```text
//! calls = units x (c+1)³        the sample grid
//!       + stencil x vertices    the per-vertex normal, stencil == 6
//! ```
//!
//! It came back with a stencil of exactly 6 on all twelve of `M-377`'s arms,
//! which is what makes it a model to extrapolate against rather than a curve to
//! fit. Applied to the **timed trace** the units are the *dirty* chunks and the
//! vertices are the ones the remesh emitted, so the trace's remesh work is
//!
//! ```text
//! W(c) = dirty_chunks(c) x (c+1)³ + 6 x remesh_vertices(c)
//! ```
//!
//! and this harness **asserts that identity against a counted field** rather
//! than assuming it - control 3 below, run on the remesh as well as on the
//! build. `W` is therefore measured exactly, not modelled.
//!
//! The time model is `M-377`'s two terms, in the same order that entry put them:
//!
//! ```text
//! total_ms(c) = mark_ms          flat in c   ("mark_edit scans the same world
//!                                              region whatever the partition")
//!             + k_field x W(c)   the remesh  ("everything is remesh_ms")
//! ```
//!
//! Both parameters come from **`docs/experiments/p-72.csv` as committed**, read
//! off disk at run time - `P-70`'s discipline, because a number in a comment is
//! a number that drifts:
//!
//! - `mark_ms` is the committed value at c = 2, the nearest arm. The model says
//!   the term does not move; taking the finest committed arm is the reading most
//!   favourable to the model, and the run reports what it cost.
//! - `k_field` is a least-squares-through-origin fit of the committed
//!   `remesh_ms` at c = 2 **and** c = 4 against this run's own `W(2)` and
//!   `W(4)`. Two anchors rather than one so that the fit has a residual: a
//!   single anchor would make `model_error` identically zero at that arm, which
//!   is `P-70`'s C3 - a clause whose predicate is implied by the construction.
//!   The residuals are emitted as `k_residual` and are the model's *in-sample*
//!   error; c = 1 is out of sample.
//!
//! Anchoring on the committed arms is only legitimate if this run's c = 2 and
//! c = 4 arms **are** those arms, so `dirty_chunks_total` is asserted equal to
//! the committed integers (`M-279`: the new run must agree with the old one on
//! everything that is not a clock).
//!
//! # The share line, recomputed before the run
//!
//! This experiment moves nothing - it measures a curve - so there is no runtime
//! share to bound. C2's 10% bar, though, has a reachability question that is
//! arithmetic and answerable from the committed file alone, and it is about the
//! **flat-mark term**, not about the remesh term C2 is nominally testing.
//!
//! `mark_ms` in `p-72.csv` is 3.6734 / 3.6840 / 3.7205 / 3.8009 / 4.0234 /
//! 5.0104 at c = 64 / 32 / 16 / 8 / 4 / 2 on `gyroid`. That is flat to 3% over
//! 64-8 and then **rises 24.5% from c = 4 to c = 2**. Mark is 52.1% of
//! `total_ms` at c = 2 (5.0104 / 9.6252). If the same 1.245x step repeats from
//! c = 2 to c = 1, a *flat* mark term under-predicts the total by
//! 0.245 x 0.521 = **12.8%** before the remesh term is even consulted - above
//! the 10% bar, from the half of the model that is a constant. On
//! `fbm_terrain` the same step is 14.8476 / 13.8520 = 1.072x and mark is 47.9%
//! of the total, so the same failure mode costs only 3.4% and C2 is reachable
//! there. **So C2 is at risk on `gyroid` for a reason that has nothing to do
//! with c = 1's remesh cost, and the run should be read with
//! `predicted_mark_ms` beside `mark_ms`.** Said before running, per `x51`.
//!
//! # Controls
//!
//! All three of `p-72`'s, inherited:
//!
//! - **Every arm must mesh something**, and the assertion is `M-377`'s verbatim
//!   `VOID: chunk N marked no dirty chunk in 11 edits`. That one string caught
//!   four distinct fixture defects in `p-72` - a dig path that missed
//!   `fbm_terrain`'s sheet, a dig path that missed `gyroid`'s undulation, an
//!   edit box computed against an origin of zero, and a world in the positive
//!   octant. It is the registered vacuity control for this row.
//!   The registration's own wording is stronger than the total `p-72`
//!   asserted - *"the 1³ arm must mark dirty chunks on **every edit** or the
//!   row is void"* - so `min_dirty_chunks_per_edit` is a column and is asserted
//!   too. Ten productive edits and one that dug air pass a total, and that is
//!   the exact shape of `p-72`'s third defect.
//! - **Every arm must produce the same surface**: the sorted multiset of vertex
//!   positions quantised to `cell_size x 1e-6`, compared across 1³, 2³ and 4³.
//!   Raw vertex counts cannot be compared (chunk faces duplicate vertices,
//!   `A-015`, `M-220`) and bit patterns cannot be compared (`M-32`, one ulp).
//! - **Every field evaluation is accounted for**, on the build *and* on the
//!   remesh, against `units x (c+1)³ + stencil x vertices` with the stencil
//!   derived from the remainder and asserted to be 6.
//!
//! # Deviations from `p-72`'s fixture, both recorded
//!
//! 1. `p-72` stores every chunk's initial `MeshBuffer` in a `BTreeMap` and uses
//!    it only to sum `build_vertices` before dropping it. At c = 1 that map is
//!    2,097,152 live `MeshBuffer`s for a number that is a running total, so the
//!    sum is accumulated directly. No measured quantity changes.
//! 2. The field-call accounting for the trace runs as a **fourth, untimed pass**
//!    over the same eleven edits. `p-72` counts only the untimed build, because
//!    a `Cell` increment inside the timed region is a perturbation. The pass is
//!    deterministic and its dirty-chunk count is asserted equal to the timed
//!    rep's.

#![allow(
    clippy::cast_possible_truncation,
    clippy::cast_precision_loss,
    clippy::cast_sign_loss,
    clippy::print_literal,
    clippy::too_many_lines
)]

mod common;

use std::cell::Cell;
use std::collections::{BTreeMap, BTreeSet};
use std::path::Path;
use std::time::Instant;

use isomesh::chunk::dirty::{DirtySet, mark_edit};
use isomesh::chunk::{ChunkId, ChunkLayout};
use isomesh::extractor::Extractor;
use isomesh::fields::{FbmTerrain, Gyroid};
use isomesh::marching_cubes::MarchingCubes;
use isomesh::{MeshBuffer, Sdf};

/// Cells per axis in the whole world, at every granularity. `p-72`'s, unchanged
/// - the total cell count is the thing a granularity sweep holds fixed.
const WORLD_CELLS: u32 = 128;

/// The knob, extended below `p-72`'s range.
///
/// `M-377` reported an interior optimum at 4³ whose interiority rests entirely
/// on 2³ losing, and 2³ was itself added because the registered 8³-64³ range
/// gave a boundary minimum. 1³ is the last granularity there is: a chunk of one
/// cell re-samples its corner planes at `((1+1)/1)³ = 8.0`, which is 2.37x the
/// penalty at 2³ and 4.1x the penalty at 4³.
const GRANULARITIES: [u32; 3] = [1, 2, 4];

/// The granularities this run shares with `p-72`, and therefore the arms whose
/// committed `remesh_ms` can anchor the model.
const ANCHORS: [u32; 2] = [2, 4];

/// World extent per axis. `p-72`'s.
const EXTENT: f64 = 4.0;

/// World origin. `p-72`'s - centred on the reference fields' own domain centre,
/// because the positive octant does not contain `fbm_terrain`'s sheet.
const ORIGIN: f64 = -EXTENT * 0.5;

/// Brush radius in **cells**. `p-72`'s six.
const BRUSH_CELLS: f64 = 6.0;

/// Edits in the trace. `p-72`'s eleven, and the number in the vacuity control's
/// message.
const EDITS: usize = 11;

/// Traces per arm, median taken. `p-72`'s three.
///
/// Three is practical at 1³ and was budgeted before the run: what is expensive
/// at 1³ is the 2,097,152-chunk initial build and the 2,097,152-chunk identity
/// pass, and neither of those is inside the repeated trace. The trace itself
/// touches only the dirty set, which at 1³ is a few thousand chunks.
const REPS: usize = 3;

/// C2's bar.
const MODEL_BAR: f64 = 0.10;

/// The stencil `M-377` derived on all twelve of its arms.
const STENCIL: u64 = 6;

/// A field wrapper that counts `sample` calls. `p-72`'s.
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

/// A sphere subtracted from a field: `max(field, -(|p - c| - r))`. `p-72`'s.
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
    /// Field calls made by the untimed initial full build.
    build_calls: u64,
    build_vertices: u64,
    /// Derived from `build_calls`, not assumed. `M-377` got 6 on twelve arms.
    build_stencil: u64,
    dirty_chunks: u64,
    remeshed_cells: u64,
    /// Field calls `mark_edit` made across the trace, counted in the untimed
    /// accounting pass. Two per corner sample, before and after.
    mark_calls: u64,
    /// Field calls the remesh made across the trace. This is `W(c)`.
    remesh_calls: u64,
    remesh_vertices: u64,
    remesh_stencil: u64,
    /// Cells `mark_edit` classified across the trace, from its own
    /// `EditReport`. Constant in `c` — the edit box is a world-cell box.
    region_cells: u64,
    /// Chunks `mark_edit` had to *consider* across the trace, from its own
    /// `EditReport`. This is the term `M-377`'s flat-mark premise does not
    /// have: it is the region's cell count divided by `c³`, so it explodes as
    /// the partition gets finer even though the region and its field work do
    /// not move at all.
    region_chunks: u64,
    /// The fewest chunks any single edit marked.
    ///
    /// The registered vacuity control is *"the 1³ arm must mark dirty chunks on
    /// **every edit** or the row is void"*, which is strictly stronger than a
    /// non-zero total: ten productive edits and one that missed would pass a
    /// total and is exactly the shape of `M-377`'s third fixture defect, where
    /// a fixed dig height found the surface at the probe's own `x` and nowhere
    /// else.
    min_dirty_per_edit: u64,
    mark_ms: f64,
    remesh_ms: f64,
    vertices: usize,
    triangles: usize,
    /// Distinct surface points, quantised to `cell_size * 1e-6`.
    surface: BTreeSet<[i64; 3]>,
    /// Wall clock for the whole arm, build and identity pass included. Recorded
    /// because 1³ was registered as a budget risk and the answer to "was a
    /// three-repeat median practical" has to be a number.
    wall_s: f64,
}

/// Run the whole trace at one granularity on one field.
fn run_arm<F: Sdf<Scalar = f64>>(field: &F, chunk_cells: u32) -> Arm {
    let started = Instant::now();
    let cell_size = EXTENT / f64::from(WORLD_CELLS);
    let layout = ChunkLayout::<f64>::new(chunk_cells, cell_size, [ORIGIN; 3]).expect("layout");
    let shape = layout.sample_shape().expect("shape");
    let per_axis = WORLD_CELLS / chunk_cells;
    let radius = BRUSH_CELLS * cell_size;

    // The dig path, `p-72`'s: the height is probed per edit at that edit's own
    // `x`, because a fixed height missed `gyroid`'s undulation and a path
    // through the world centre missed `fbm_terrain`'s sheet. Both were caught by
    // the vacuity control, not by review.
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

    // The cell box one edit can touch, through `layout.cell_of` because the
    // hand-rolled form assumed an origin of zero (`M-377`'s fourth defect).
    let edit_box = |step: usize| -> ([i64; 3], [i64; 3]) {
        let c = centres[step];
        let lo = layout.cell_of([0, 1, 2].map(|a| c[a] - radius)).map(|v| v - 1);
        let hi = layout.cell_of([0, 1, 2].map(|a| c[a] + radius)).map(|v| v + 1);
        (lo, hi)
    };

    let calls = Cell::new(0u64);
    let mut dirty = DirtySet::new();
    let mut mc = MarchingCubes::<f64>::new();

    // ── the initial full build: every chunk once, untimed, counted ───────────
    //
    // Not part of the trace - it is the state an edit arrives into - but its
    // field calls are what control 3 asserts the grid duplication against.
    let mut build_vertices = 0u64;
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
                build_vertices += out.positions.len() as u64;
            }
        }
    }
    let build_calls = calls.get();

    // ── the timed trace, REPS times, uncounted ───────────────────────────────
    let mut dirty_chunks = 0u64;
    let mut remeshed_cells = 0u64;
    let mut reps: Vec<(f64, f64)> = Vec::with_capacity(REPS);
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
            let (lo, hi) = edit_box(step);

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

    // Median by total; both components come from the same trace.
    reps.sort_unstable_by(|a, b| (a.0 + a.1).partial_cmp(&(b.0 + b.1)).expect("finite"));
    let (mark_ms, remesh_ms) = reps[REPS / 2];

    // ── the accounting pass: the same eleven edits, untimed, counted ─────────
    //
    // `W(c)` is what C2's model is denominated in and it has to be the *trace's*
    // field calls, not the build's. Counting inside the timed region would put a
    // `Cell` increment on every sample of the thing being timed.
    let acct = Cell::new(0u64);
    let mut mark_calls = 0u64;
    let mut remesh_calls = 0u64;
    let mut remesh_vertices = 0u64;
    let mut acct_dirty = 0u64;
    let mut region_cells = 0u64;
    let mut region_chunks = 0u64;
    let mut min_dirty_per_edit = u64::MAX;
    {
        let counted = Counted {
            field,
            calls: &acct,
        };
        for step in 0..EDITS {
            let before = Dug {
                field: &counted,
                centres: &centres[..step],
                radius,
            };
            let after = Dug {
                field: &counted,
                centres: &centres[..=step],
                radius,
            };
            let (lo, hi) = edit_box(step);

            let at_mark = acct.get();
            let report =
                mark_edit(&layout, &before, &after, lo, hi, &mut dirty).expect("mark");
            mark_calls += acct.get() - at_mark;
            region_cells += report.region_cells;
            region_chunks += report.region_chunks;

            let at_remesh = acct.get();
            let verts = &mut remesh_vertices;
            let done = dirty.mesh_dirty(&layout, |_id, origin| {
                let mut out = MeshBuffer::<f64>::new();
                let _ = mc.extract_into(&after, &shape, origin, cell_size, &mut out);
                *verts += out.positions.len() as u64;
            });
            remesh_calls += acct.get() - at_remesh;
            acct_dirty += done as u64;
            min_dirty_per_edit = min_dirty_per_edit.min(done as u64);
        }
    }
    assert_eq!(
        acct_dirty, dirty_chunks,
        "chunk {chunk_cells}: the accounting pass marked {acct_dirty} dirty chunks against the \
         timed trace's {dirty_chunks}, so the trace is not deterministic and W is not the \
         work that was timed"
    );

    // ── the surface, meshed once at the end and NOT timed ────────────────────
    //
    // The identity control cannot come from the trace: a chunk dirtied only at
    // step 3 holds a mesh of `centres[..=3]`, and which step last touched a
    // region depends on the partition. What the control asks is the seam
    // question, so the final field is meshed over every chunk here.
    let full = Dug {
        field,
        centres: &centres,
        radius,
    };
    let mut surface: BTreeSet<[i64; 3]> = BTreeSet::new();
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

    let chunks = u64::from(per_axis).pow(3);
    let samples_per_chunk = u64::from(chunk_cells + 1).pow(3);
    let build_stencil = stencil_of(
        chunk_cells,
        "build",
        build_calls,
        chunks * samples_per_chunk,
        build_vertices,
    );
    let remesh_stencil = stencil_of(
        chunk_cells,
        "remesh",
        remesh_calls,
        dirty_chunks * samples_per_chunk,
        remesh_vertices,
    );

    Arm {
        chunk_cells,
        chunks,
        samples_per_chunk,
        build_calls,
        build_vertices,
        build_stencil,
        dirty_chunks,
        remeshed_cells,
        mark_calls,
        remesh_calls,
        remesh_vertices,
        remesh_stencil,
        region_cells,
        region_chunks,
        min_dirty_per_edit,
        mark_ms,
        remesh_ms,
        vertices,
        triangles,
        surface,
        wall_s: started.elapsed().as_secs_f64(),
    }
}

/// Control 3, as a function because it now runs on two populations.
///
/// `M-377`'s two-term model asserts the grid term **exactly** and *derives* the
/// stencil from the remainder, which is a stronger control than asserting 6:
/// it proves the harness knows where every field evaluation went. The
/// registration names `stencil == 6`, so the derived value is checked against
/// it rather than merely printed.
fn stencil_of(chunk_cells: u32, what: &str, calls: u64, grid: u64, vertices: u64) -> u64 {
    assert!(
        vertices > 0,
        "chunk {chunk_cells}: the {what} emitted no vertices, so its stencil is a division by \
         zero and the two-term model has no second term to check"
    );
    assert!(
        calls >= grid,
        "chunk {chunk_cells}: the {what} made {calls} field calls, fewer than its {grid} \
         sample-grid points, so the grid was not fully evaluated"
    );
    let normals = calls - grid;
    assert_eq!(
        normals % vertices,
        0,
        "chunk {chunk_cells}: the {what}'s {normals} non-grid field calls over {vertices} \
         vertices is not a whole stencil, so the call model is incomplete"
    );
    let stencil = normals / vertices;
    assert_eq!(
        stencil, STENCIL,
        "chunk {chunk_cells}: the {what}'s derived stencil is {stencil}, not the {STENCIL} \
         M-377 got on all twelve of its arms, so this is not the same call model"
    );
    stencil
}

/// One committed `p-72` arm, read off disk.
struct Committed {
    mark_ms: f64,
    remesh_ms: f64,
    total_ms: f64,
    dirty_chunks: u64,
}

/// Read `docs/experiments/p-72.csv` as committed.
///
/// `P-70`'s discipline: C2's model parameters are an artefact in the repository,
/// not numbers typed into this file. Keyed by `(field, chunk_cells)`.
fn committed_p72(root: &Path) -> BTreeMap<(String, u32), Committed> {
    let text = std::fs::read_to_string(root.join("docs/experiments/p-72.csv"))
        .expect("docs/experiments/p-72.csv, which is where C2's model parameters live");
    let mut header: Vec<&str> = Vec::new();
    let mut out = BTreeMap::new();
    for line in text.lines() {
        if line.starts_with('#') {
            continue;
        }
        let cells: Vec<&str> = line.split(',').collect();
        if header.is_empty() {
            header = cells;
            continue;
        }
        let get = |name: &str| -> Option<&str> {
            header
                .iter()
                .position(|h| *h == name)
                .and_then(|i| cells.get(i))
                .copied()
        };
        let (Some(field), Some(c)) = (get("field"), get("chunk_cells")) else {
            continue;
        };
        let num = |name: &str| -> f64 {
            get(name)
                .and_then(|v| v.parse().ok())
                .unwrap_or(f64::NAN)
        };
        let Ok(chunk_cells) = c.parse::<u32>() else {
            continue;
        };
        out.insert(
            (field.to_string(), chunk_cells),
            Committed {
                mark_ms: num("mark_ms"),
                remesh_ms: num("remesh_ms"),
                total_ms: num("total_ms"),
                dirty_chunks: num("dirty_chunks_total") as u64,
            },
        );
    }
    out
}

type Row = Vec<(&'static str, String)>;

fn main() {
    if !std::env::args().any(|arg| arg == "--bench") {
        return;
    }

    let prereg = isomesh::experiment!("P-89");
    let root = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../..");
    let committed = committed_p72(&root);

    // ── the share line, printed before anything is measured ──────────────────
    //
    // Recomputed from the file rather than from the doc comment, for the same
    // reason C2's parameters are.
    println!("SHARE, recomputed from docs/experiments/p-72.csv before the run:");
    println!(
        "  this experiment moves nothing - it measures a curve - so no clause is \
         bounded by a share of a runtime."
    );
    for field_name in ["gyroid", "fbm_terrain"] {
        let a2 = committed
            .get(&(field_name.to_string(), 2))
            .expect("committed p-72 arm at c = 2");
        let a4 = committed
            .get(&(field_name.to_string(), 4))
            .expect("committed p-72 arm at c = 4");
        let step = a2.mark_ms / a4.mark_ms;
        let share = a2.mark_ms / a2.total_ms;
        let cost = (step - 1.0) * share;
        println!(
            "  {field_name}: mark 4->2 rises {step:.4}x and is {:.1}% of total at c=2, so a FLAT \
             mark term costs the model {:.1}% at c=1 if the step repeats -- C2's 10% bar is {}",
            100.0 * share,
            100.0 * cost,
            if cost >= MODEL_BAR {
                "UNREACHABLE from the mark term alone"
            } else {
                "reachable"
            }
        );
    }
    println!();

    let mut rows: Vec<Row> = Vec::new();
    let mut verdicts: Vec<(&'static str, bool, bool)> = Vec::new();

    for field_name in ["gyroid", "fbm_terrain"] {
        let mut arms: Vec<Arm> = Vec::new();
        for c in GRANULARITIES {
            let arm = match field_name {
                "gyroid" => run_arm(&Gyroid::<f64>::canonical(), c),
                _ => run_arm(&FbmTerrain::<f64>::canonical(), c),
            };
            println!(
                "  {field_name} c={:<3} built and traced in {:.1} s",
                arm.chunk_cells, arm.wall_s
            );
            arms.push(arm);
        }

        // ── control 1: every arm meshed something ────────────────────────────
        //
        // `M-377`'s string, verbatim, because it is the registered vacuity
        // control for this row and four of `p-72`'s fixture defects announced
        // themselves through it.
        for a in &arms {
            assert!(
                a.dirty_chunks > 0,
                "VOID: chunk {} marked no dirty chunk in {EDITS} edits on {field_name}, so its \
                 time measures an empty trace",
                a.chunk_cells
            );
            // The registration's own wording, which is stronger than the total
            // `p-72` asserted: *"the 1³ arm must mark dirty chunks on every
            // edit or the row is void"*. A trace where one edit missed still
            // reports a healthy total.
            assert!(
                a.min_dirty_per_edit > 0,
                "VOID: chunk {} has an edit that marked no dirty chunk on {field_name}, so at \
                 least one of the {EDITS} edits dug in empty space",
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
        // Quantised to `cell_size * 1e-6`: four orders finer than any real
        // surface movement, four orders coarser than an ulp near 4.0 (~8.9e-16,
        // `M-32`). Raw vertex counts must differ across partitions - chunk faces
        // duplicate (`A-015`, `M-220`) - and at 1³ every internal face is a
        // chunk face, so this is the granularity where a seam defect would be
        // loudest.
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
        // Asserted inside `run_arm` by `stencil_of`, on the build and on the
        // remesh. What is left here is the registered grid duplication, which at
        // c = 1 is exactly 8.
        let baseline = u64::from(WORLD_CELLS).pow(3) as f64;
        for a in &arms {
            let grid = a.chunks * a.samples_per_chunk;
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
                "  {field_name} c={:<3} build {} calls = {grid} grid ({measured:.4}x) + {} x {} \
                 verts; trace {} mark + {} remesh = {} grid + {} x {} verts",
                a.chunk_cells,
                a.build_calls,
                a.build_stencil,
                a.build_vertices,
                a.mark_calls,
                a.remesh_calls,
                a.dirty_chunks * a.samples_per_chunk,
                a.remesh_stencil,
                a.remesh_vertices
            );
        }

        // ── the anchors must be the committed arms ───────────────────────────
        //
        // `k_field` is fitted against committed `remesh_ms` values, which is only
        // legitimate if this run's c = 2 and c = 4 arms are the arms that
        // produced them. `M-279`: everything that is not a clock must agree.
        for c in ANCHORS {
            let arm = arms
                .iter()
                .find(|a| a.chunk_cells == c)
                .expect("an anchor arm was swept");
            let old = committed
                .get(&(field_name.to_string(), c))
                .expect("committed p-72 anchor arm");
            assert_eq!(
                arm.dirty_chunks, old.dirty_chunks,
                "chunk {c} on {field_name} marked {} dirty chunks against p-72's committed {}, \
                 so this is not the same trace and its committed remesh_ms cannot anchor the \
                 model",
                arm.dirty_chunks, old.dirty_chunks
            );
        }

        // ── C2's model, fitted from the committed file ───────────────────────
        //
        // Least squares through the origin over the two anchors:
        //   k = sum(remesh_committed * W) / sum(W^2)
        // Two anchors so the fit has a residual at both; one would make
        // `model_error` identically zero there, which is a HELD with no
        // instrument.
        let mut num = 0.0f64;
        let mut den = 0.0f64;
        for c in ANCHORS {
            let arm = arms
                .iter()
                .find(|a| a.chunk_cells == c)
                .expect("an anchor arm was swept");
            let old = committed
                .get(&(field_name.to_string(), c))
                .expect("committed p-72 anchor arm");
            let w = arm.remesh_calls as f64;
            num += old.remesh_ms * w;
            den += w * w;
        }
        let k = num / den;
        let flat_mark = committed
            .get(&(field_name.to_string(), 2))
            .expect("committed p-72 arm at c = 2")
            .mark_ms;
        println!(
            "\n  {field_name}: k = {:.6} ns per remesh field call, flat mark term {flat_mark:.4} \
             ms, both from docs/experiments/p-72.csv",
            k * 1e6
        );
        for c in ANCHORS {
            let arm = arms
                .iter()
                .find(|a| a.chunk_cells == c)
                .expect("an anchor arm was swept");
            let old = committed
                .get(&(field_name.to_string(), c))
                .expect("committed p-72 anchor arm");
            let fitted = k * arm.remesh_calls as f64;
            println!(
                "    in-sample at c={c}: model remesh {fitted:.4} ms against committed \
                 {:.4} ms, residual {:.2}%",
                old.remesh_ms,
                100.0 * (fitted - old.remesh_ms).abs() / old.remesh_ms
            );
        }

        // ── the verdicts ─────────────────────────────────────────────────────
        let at = |c: u32| -> &Arm {
            arms.iter()
                .find(|a| a.chunk_cells == c)
                .expect("a swept arm")
        };
        let one = at(1);
        let two = at(2);
        let c1 = one.mark_ms + one.remesh_ms > two.mark_ms + two.remesh_ms;
        let predicted_one = flat_mark + k * one.remesh_calls as f64;
        let measured_one = one.mark_ms + one.remesh_ms;
        let error_one = (measured_one - predicted_one).abs() / predicted_one;
        let c2 = error_one < MODEL_BAR;
        verdicts.push((field_name, c1, c2));

        for a in &arms {
            let total = a.mark_ms + a.remesh_ms;
            let dup = (a.chunks * a.samples_per_chunk) as f64 / baseline;
            let predicted_mark = flat_mark;
            let predicted_remesh = k * a.remesh_calls as f64;
            let predicted = predicted_mark + predicted_remesh;
            let error = (total - predicted).abs() / predicted;
            let k_residual = committed.get(&(field_name.to_string(), a.chunk_cells)).map_or(
                f64::NAN,
                |old| (predicted_remesh - old.remesh_ms).abs() / old.remesh_ms,
            );
            println!(
                "{:>12} c={:<3} dup {dup:>7.4}x dirty {:>6} mark {:>8.3} remesh {:>8.3} total \
                 {:>9.3} predicted {:>9.3} err {:>7.2}%",
                field_name,
                a.chunk_cells,
                a.dirty_chunks,
                a.mark_ms,
                a.remesh_ms,
                total,
                predicted,
                100.0 * error
            );
            // Where a model error lives is not derivable from `model_error`,
            // and the two halves of a two-term model can be wrong separately.
            println!(
                "{:>12}      mark {:>8.3} vs flat {:>8.3} ({:>7.2}%) over {} region cells in \
                 {} region chunks | remesh {:>8.3} vs model {:>8.3} ({:>6.2}%)",
                "",
                a.mark_ms,
                predicted_mark,
                100.0 * (a.mark_ms - predicted_mark).abs() / predicted_mark,
                a.region_cells,
                a.region_chunks,
                a.remesh_ms,
                predicted_remesh,
                100.0 * (a.remesh_ms - predicted_remesh).abs() / predicted_remesh
            );
            rows.push(vec![
                ("field", field_name.to_string()),
                ("chunk_cells", a.chunk_cells.to_string()),
                ("world_cells", WORLD_CELLS.to_string()),
                ("chunks", a.chunks.to_string()),
                ("sample_duplication", format!("{dup:.6}")),
                ("predicted_ms", format!("{predicted:.4}")),
                ("measured_ms", format!("{total:.4}")),
                ("model_error", format!("{error:.6}")),
                ("mark_ms", format!("{:.4}", a.mark_ms)),
                ("remesh_ms", format!("{:.4}", a.remesh_ms)),
                ("total_ms", format!("{total:.4}")),
                ("dirty_chunks_total", a.dirty_chunks.to_string()),
                // Extras. The registered set names no instrument for *where* a
                // model error lives, and the whole point of a two-term model is
                // that the two terms can be wrong separately.
                ("predicted_mark_ms", format!("{predicted_mark:.4}")),
                ("predicted_remesh_ms", format!("{predicted_remesh:.4}")),
                (
                    "mark_model_error",
                    format!("{:.6}", (a.mark_ms - predicted_mark).abs() / predicted_mark),
                ),
                (
                    "remesh_model_error",
                    format!(
                        "{:.6}",
                        (a.remesh_ms - predicted_remesh).abs() / predicted_remesh
                    ),
                ),
                ("k_ns_per_call", format!("{:.6}", k * 1e6)),
                ("k_residual", format!("{k_residual:.6}")),
                ("samples_per_chunk", a.samples_per_chunk.to_string()),
                ("edits", EDITS.to_string()),
                ("reps", REPS.to_string()),
                ("remeshed_cells_total", a.remeshed_cells.to_string()),
                ("build_field_calls", a.build_calls.to_string()),
                ("build_vertices", a.build_vertices.to_string()),
                ("build_stencil", a.build_stencil.to_string()),
                ("mark_field_calls", a.mark_calls.to_string()),
                ("remesh_field_calls", a.remesh_calls.to_string()),
                ("remesh_vertices", a.remesh_vertices.to_string()),
                ("remesh_stencil", a.remesh_stencil.to_string()),
                ("model_work_calls", a.remesh_calls.to_string()),
                // The bookkeeping term the flat-mark premise does not have,
                // read from `mark_edit`'s own `EditReport` rather than derived.
                ("region_cells_total", a.region_cells.to_string()),
                ("region_chunks_total", a.region_chunks.to_string()),
                (
                    "min_dirty_chunks_per_edit",
                    a.min_dirty_per_edit.to_string(),
                ),
                ("vertices", a.vertices.to_string()),
                ("triangles", a.triangles.to_string()),
                ("distinct_surface_points", a.surface.len().to_string()),
                (
                    "vertex_duplication",
                    format!("{:.6}", a.vertices as f64 / a.surface.len() as f64),
                ),
                ("ms_per_edit", format!("{:.6}", total / EDITS as f64)),
                ("arm_wall_s", format!("{:.3}", a.wall_s)),
            ]);
        }
        println!(
            "  {field_name}: C1 1³ {:.3} ms against 2³ {:.3} ms -> {}",
            measured_one,
            two.mark_ms + two.remesh_ms,
            if c1 {
                "1³ is WORSE, interior optimum confirmed"
            } else {
                "1³ WINS, the duplication model is wrong"
            }
        );
        println!(
            "  {field_name}: C2 measured {measured_one:.3} ms against predicted \
             {predicted_one:.3} ms -> {:.2}% ({})\n",
            100.0 * error_one,
            if c2 { "within 10%" } else { "above 10%" }
        );
    }

    let c1_all = verdicts.iter().all(|(_, c1, _)| *c1);
    let c2_all = verdicts.iter().all(|(_, _, c2)| *c2);
    println!(
        "C1 1³ is worse than 2³ on both fields: {:?} -> {}",
        verdicts.iter().map(|(f, c, _)| (*f, *c)).collect::<Vec<_>>(),
        if c1_all { "HELD" } else { "FALSIFIED" }
    );
    println!(
        "C2 the two-term model extrapolates to c = 1 within 10%: {:?} -> {}",
        verdicts.iter().map(|(f, _, c)| (*f, *c)).collect::<Vec<_>>(),
        if c2_all { "HELD" } else { "FALSIFIED" }
    );

    let mut finished: Vec<Row> = Vec::new();
    for row in &mut rows {
        row.push(("c1_holds", c1_all.to_string()));
        row.push(("c2_holds", c2_all.to_string()));
        finished.push(row.clone());
    }

    common::experiment::run(prereg, |run| {
        for row in &finished {
            run.record(row);
        }
    });
}
