//! B-003's tests, run headlessly.
//!
//! The ticket's acceptance is *"meshing a large volume does not stall the render
//! loop — show it in the frame-time graph"*, and a frame-time graph is not a
//! test. What is testable is the property the graph would show: **the work is
//! not on the main thread, and what remains there is bounded per frame**. Both
//! are asserted here, and both have a companion that proves they could fail.

use super::*;

use bevy_app::App;
use bevy_app::TaskPoolPlugin;
use bevy_asset::AssetPlugin;
use isomesh::chunk::ChunkLayout;
use isomesh::fields::Sphere;

/// A field that is slow on purpose, so "off the main thread" is measurable.
///
/// A sphere at these grid sizes meshes in microseconds, which is too fast to
/// tell a stalled frame from a quick one. This one burns a fixed amount of work
/// per sample, so an extraction takes long enough that a main-thread
/// implementation would be visibly slower than an off-thread one.
struct SlowSphere {
    inner: Sphere<f32>,
    spin: u32,
}

impl Sdf for SlowSphere {
    type Scalar = f32;
    fn sample(&self, p: [f32; 3]) -> f32 {
        let mut accumulator = 0.0f32;
        for i in 0..self.spin {
            accumulator += (i as f32).sin();
        }
        // `accumulator` is multiplied by zero rather than discarded, so the loop
        // cannot be optimised away and the cost is real.
        self.inner.sample(p) + accumulator * 0.0
    }
}

/// The deliberately slow field's spin count: tens of milliseconds per extraction
/// in a debug build, so twelve of them are still running when a frame returns.
const SPIN: u32 = 4_000;

fn app(spin: u32, chunks: i32) -> App {
    let mut app = App::new();
    app.add_plugins((TaskPoolPlugin::default(), AssetPlugin::default()))
        .add_plugins(IsomeshPlugin);

    let layout = ChunkLayout::new(8, 0.25, [0.0; 3]).expect("a valid layout");
    let volume = app
        .world_mut()
        .spawn(VoxelVolume::new(
            layout,
            SlowSphere {
                inner: Sphere::<f32>::canonical(),
                spin,
            },
        ))
        .id();

    for x in 0..chunks {
        app.world_mut().spawn((
            VoxelChunk {
                id: ChunkId::new([x, 0, 0]),
                volume,
            },
            NeedsRemesh,
        ));
    }
    app
}

fn meshed(app: &mut App) -> usize {
    app.world_mut()
        .query::<&ChunkMesh>()
        .iter(app.world())
        .count()
}

/// Run `app` until `done` holds, or give up.
///
/// **Spinning `app.update()` in a tight loop is not waiting.** The extractions
/// run on [`AsyncComputeTaskPool`] threads, and a loop that never yields grants
/// them almost no wall clock: B-005 measured 500 iterations completing in ~20 ms
/// against work that needs ~40 ms, draining 2 chunks of 12 and failing 8 runs out
/// of 8. The same test passed in **39** iterations once each one slept a
/// millisecond — twelve times fewer iterations and twice the wall time — which is
/// what identified the loop rather than the plugin as the defect.
///
/// So the bound is a deadline, not an iteration count. Returns the frames used,
/// because two tests assert on that.
fn drain_until(app: &mut App, mut done: impl FnMut(&mut App) -> bool) -> usize {
    const DEADLINE: Duration = Duration::from_secs(10);
    const SLICE: Duration = Duration::from_micros(250);

    let started = Instant::now();
    let mut frames = 0;
    while started.elapsed() < DEADLINE {
        app.update();
        frames += 1;
        if done(app) {
            break;
        }
        std::thread::sleep(SLICE);
    }
    frames
}

#[test]
fn every_chunk_eventually_gets_a_mesh() {
    let mut app = app(0, 12);
    drain_until(&mut app, |app| meshed(app) == 12);
    assert_eq!(meshed(&mut app), 12, "the queue never drained");

    // And nothing is left marked or in flight, so the drain is complete rather
    // than merely far along.
    let waiting = app
        .world_mut()
        .query_filtered::<(), With<NeedsRemesh>>()
        .iter(app.world())
        .count();
    let running = app
        .world_mut()
        .query_filtered::<(), With<MeshingTask>>()
        .iter(app.world())
        .count();
    assert_eq!((waiting, running), (0, 0));
}

#[test]
fn spawning_the_work_does_not_do_the_work() {
    // The acceptance property, stated without a clock. One update queues twelve
    // extractions of a deliberately slow field, and `apply_finished_meshes` runs
    // in that same update, right after `spawn_meshing_tasks`. If any extraction
    // ran on the main thread -- inline, or blocked on at the spawn -- it would
    // be finished and applied by the time `update` returns. So an `update` that
    // returns with every task still running and nothing meshed did its work
    // somewhere else.
    //
    // This replaces a comparison of the frame's wall time against one
    // extraction timed on the main thread. On a debug build the frame's own
    // bookkeeping is of the same order as the extraction it was supposed to be
    // a quarter of, and on a busy CI runner the ratio is noise: it failed three
    // runs in a row (11.6-21.9 ms of frame against 25-31 ms of extraction) on a
    // change that touched neither system.
    let mut app = app(SPIN, 12);
    app.update();

    let running = app
        .world_mut()
        .query_filtered::<(), With<MeshingTask>>()
        .iter(app.world())
        .count();
    assert_eq!(
        (meshed(&mut app), running),
        (0, 12),
        "the update that queued twelve extractions returned with some of them \
         already done: that is work that ran on the main thread"
    );

    // And the queue does drain, so the twelve tasks were real work, off-thread.
    drain_until(&mut app, |app| meshed(app) == 12);
    assert_eq!(meshed(&mut app), 12, "the queue never drained");
}

#[test]
fn the_budget_bounds_what_lands_per_frame() {
    // A zero budget still applies one mesh per frame -- never zero, or the queue
    // livelocks -- so twelve chunks take at least twelve frames to land.
    let mut app = app(0, 12);
    app.insert_resource(MeshBudget {
        per_frame: Duration::ZERO,
        max_in_flight: 64,
    });

    let mut worst = 0usize;
    let mut before = 0usize;
    let frames = drain_until(&mut app, |app| {
        let after = meshed(app);
        worst = worst.max(after - before);
        before = after;
        after == 12
    });
    assert_eq!(meshed(&mut app), 12);
    assert_eq!(
        worst, 1,
        "a zero budget applied {worst} meshes in one frame"
    );
    assert!(frames >= 12, "twelve meshes landed in {frames} frames");
}

#[test]
fn a_generous_budget_is_not_bounded_the_same_way() {
    // M-44 again: the bound above has to prove it could have been looser. With
    // room to work, more than one mesh lands per frame.
    let mut app = app(0, 12);
    app.insert_resource(MeshBudget {
        per_frame: Duration::from_secs(1),
        max_in_flight: 64,
    });

    let mut best = 0usize;
    let mut before = 0usize;
    drain_until(&mut app, |app| {
        let after = meshed(app);
        best = best.max(after - before);
        before = after;
        after == 12
    });
    assert!(
        best > 1,
        "a one-second budget still applied only {best} mesh per frame, so the \
         zero-budget test is not measuring the budget"
    );
}

#[test]
fn queueing_is_capped_so_memory_does_not_run_ahead_of_the_frame() {
    // Every task in flight holds a MeshBuilder. Unbounded queueing is its own
    // stall, and one the frame budget never sees.
    let mut app = app(2_000, 40);
    app.insert_resource(MeshBudget {
        per_frame: Duration::from_micros(1),
        max_in_flight: 3,
    });

    app.update();
    let running = app
        .world_mut()
        .query_filtered::<(), With<MeshingTask>>()
        .iter(app.world())
        .count();
    assert!(running <= 3, "{running} tasks in flight against a cap of 3");

    let stats = *app.world().resource::<MeshStats>();
    assert!(stats.waiting > 0, "nothing was held back by the cap");
}

/// The vertex positions of the mesh `chunk` is showing, if it has one.
fn shown(app: &App, chunk: Entity) -> Option<Vec<[f32; 3]>> {
    let handle = &app.world().entity(chunk).get::<ChunkMesh>()?.0;
    let mesh = app.world().resource::<Assets<Mesh>>().get(handle)?;
    let positions = mesh
        .attribute(Mesh::ATTRIBUTE_POSITION)
        .and_then(|values| values.as_float3())?;
    Some(positions.to_vec())
}

/// Issue #12. The first version asserted the end state (a `ChunkMesh`, no
/// `NeedsRemesh`, no `MeshingTask`), which the swallowing implementation also
/// leaves, so it passed with the edit dropped. This one observes the two things
/// a swallow changes: whether a second extraction ran, and which field it read.
#[test]
fn an_edit_during_extraction_is_requeued_rather_than_swallowed() {
    // NeedsRemesh comes off when the task is *spawned*, so re-marking a chunk
    // mid-extraction queues a second pass. If it came off when the mesh landed,
    // the second mark would be cleared by the first task's completion and the
    // edit would vanish.
    let mut app = app(2_000, 1);
    app.update();

    let chunk = app
        .world_mut()
        .query_filtered::<Entity, With<VoxelChunk>>()
        .iter(app.world())
        .next()
        .expect("one chunk");
    let volume = app
        .world_mut()
        .query_filtered::<Entity, With<VoxelVolume>>()
        .iter(app.world())
        .next()
        .expect("one volume");
    assert!(
        app.world().entity(chunk).contains::<MeshingTask>(),
        "the chunk should be extracting"
    );
    assert!(
        !app.world().entity(chunk).contains::<NeedsRemesh>(),
        "NeedsRemesh should come off at spawn"
    );

    // What each field should make of this chunk, extracted directly. The two
    // must differ here, or the comparisons below could not tell a stale mesh
    // from a fresh one.
    let layout = ChunkLayout::new(8, 0.25, [0.0; 3]).expect("a valid layout");
    let id = ChunkId::new([0, 0, 0]);
    let old_field = SlowSphere {
        inner: Sphere::<f32>::canonical(),
        spin: 2_000,
    };
    let new_field = Sphere::<f32> {
        center: [0.0; 3],
        radius: 1.5,
    };
    let positions_of = |field: &dyn VolumeField| {
        extract_chunk(field, &layout, id, Extractor::default())
            .expect("the reference extraction")
            .positions()
            .to_vec()
    };
    let (old_mesh, new_mesh) = (positions_of(&old_field), positions_of(&new_field));
    assert!(
        !old_mesh.is_empty() && !new_mesh.is_empty() && old_mesh != new_mesh,
        "the two fields must mesh this chunk differently ({} and {} vertices)",
        old_mesh.len(),
        new_mesh.len()
    );

    // The edit arrives while the first extraction is still running: the field
    // is swapped and the chunk marked the way a game does it.
    edit(&mut app, volume, new_field);
    assert_eq!(marked(&mut app), vec![id], "the edit should mark the chunk");

    let mut respawned = 0;
    let mut first_shown = None;
    drain_until(&mut app, |app| {
        // `spawned` resets every frame, so the extractions after the mark are
        // the sum over the frames driven here.
        respawned += app.world().resource::<MeshStats>().spawned;
        if first_shown.is_none() {
            first_shown = shown(app, chunk);
        }
        !app.world().entity(chunk).contains::<NeedsRemesh>()
            && !app.world().entity(chunk).contains::<MeshingTask>()
    });
    assert!(
        app.world().entity(chunk).contains::<ChunkMesh>(),
        "the chunk never finished"
    );
    assert!(
        !app.world().entity(chunk).contains::<NeedsRemesh>(),
        "the re-mark was never serviced"
    );
    assert_eq!(
        respawned, 1,
        "the edit that arrived mid-extraction was extracted {respawned} times after the mark, not once"
    );

    // `replace_field` promises that a running task finishes against the field
    // it was spawned with; the re-extraction must read the one swapped in.
    let last = shown(&app, chunk);
    assert!(
        first_shown.as_deref() == Some(old_mesh.as_slice()),
        "the extraction already running should have finished against the old field \
         ({} vertices shown, {} expected)",
        first_shown.as_ref().map_or(0, Vec::len),
        old_mesh.len()
    );
    assert!(
        last.as_deref() == Some(new_mesh.as_slice()),
        "the chunk ended up showing a mesh the edited field did not make \
         ({} vertices shown, {} expected)",
        last.as_ref().map_or(0, Vec::len),
        new_mesh.len()
    );
}

#[test]
fn a_chunk_whose_volume_is_gone_stops_asking() {
    let mut app = app(0, 4);
    let volume = app
        .world_mut()
        .query_filtered::<Entity, With<VoxelVolume>>()
        .iter(app.world())
        .next()
        .expect("one volume");
    app.world_mut().entity_mut(volume).despawn();

    app.update();
    let still_marked = app
        .world_mut()
        .query_filtered::<(), With<NeedsRemesh>>()
        .iter(app.world())
        .count();
    assert_eq!(still_marked, 0, "orphaned chunks kept requesting a mesh");
}

#[test]
fn the_subgrid_extractor_is_reachable_through_the_component() {
    // The plugin's whole point is that the choice is data. This is the one
    // extractor whose configuration is part of that choice (M-95).
    let mut app = App::new();
    app.add_plugins((TaskPoolPlugin::default(), AssetPlugin::default()))
        .add_plugins(IsomeshPlugin);

    let layout = ChunkLayout::new(6, 0.25, [0.0; 3]).expect("a valid layout");
    let volume = app
        .world_mut()
        .spawn(
            VoxelVolume::new(layout, Sphere::<f32>::canonical())
                .with_extractor(Extractor::Subgrid { samples: 4 }),
        )
        .id();
    app.world_mut().spawn((
        VoxelChunk {
            id: ChunkId::new([0, 0, 0]),
            volume,
        },
        NeedsRemesh,
    ));

    drain_until(&mut app, |app| meshed(app) == 1);
    assert_eq!(meshed(&mut app), 1);
}

/// Two adjacent chunks, meshed independently, welded, counting boundary edges
/// that lie in their shared plane.
///
/// Excludes the two-chunk block's own six walls, which are the block ending
/// rather than a seam failing. **The first version of this omitted `y`**, and a
/// field extending past the chunk's vertical range clipped at the top and put one
/// of those edges in the seam plane -- reported as a Marching Cubes seam failure
/// that was nothing of the sort.
fn seam_open_edges<S: VolumeField>(field: &S, extractor: Extractor) -> usize {
    let layout = ChunkLayout::new(16, 0.25, [0.0; 3]).expect("a valid layout");
    let mut all = isomesh::MeshBuffer::<f32>::new();
    for id in [ChunkId::new([0, 0, 0]), ChunkId::new([1, 0, 0])] {
        let built =
            extract_chunk(field, &layout, id, extractor).expect("a chunk extraction failed");
        let mut buf = isomesh::MeshBuffer::<f32>::new();
        buf.positions.extend_from_slice(built.positions());
        buf.normals.extend_from_slice(built.normals());
        buf.indices.extend_from_slice(built.indices());
        all.append(&buf)
            .expect("two chunks fit the u32 index space");
    }
    let h = layout.cell_size();
    isomesh::weld::Welder::<f32>::new()
        .weld(&mut all, isomesh::weld::epsilon_for(h))
        .expect("weld");
    let cfg = isomesh::validate::ValidateConfig::from_cell_size(f64::from(h)).expect("cfg");
    let (_report, features) =
        isomesh::validate::validate_features(&all.positions, &all.indices, &cfg);

    // The block's own six walls, derived from the layout rather than written
    // out. **Both hand-written versions of this were wrong**: the first omitted
    // `y` entirely, and the second gave `z` an upper bound of 8.0 when a chunk
    // is only 4.0 deep -- so the `z` wall was never excluded and its edges were
    // counted as seam failures. Deriving the bounds removes the opportunity.
    let span = layout.cell_size() * layout.cells() as f32;
    let hi = [span * 2.0, span, span];
    let (seam, tol) = (span, h * 0.25);
    let outer = |v: [f32; 3]| (0..3).any(|a| v[a] < tol || v[a] > hi[a] - tol);
    features
        .boundary_edges
        .iter()
        .filter(|e| {
            let (p, q) = (all.positions[e[0] as usize], all.positions[e[1] as usize]);
            !(outer(p) && outer(q)) && (p[0] - seam).abs() < tol && (q[0] - seam).abs() < tol
        })
        .count()
}

/// A field that crosses the seam plane along its whole length.
struct Waves;
impl Sdf for Waves {
    type Scalar = f32;
    fn sample(&self, p: [f32; 3]) -> f32 {
        p[1] - 0.9 * (p[0] * 0.7).sin() * (p[2] * 0.6).cos() - 2.0
    }
}

#[test]
fn the_seam_counts_are_pinned() {
    // B-006. These are not aspirations: Marching Cubes closes a chunk seam
    // because its vertices sit on grid *edges*, which both chunks compute from
    // identical corner values; the dual methods place one vertex per cell
    // *interior* and a boundary quad needs the neighbour's, so they cannot.
    //
    // Pinned as exact non-zero numbers rather than `> 0`, per M-4's rule: a
    // known defect with a number moves only when someone means it to.
    assert_eq!(
        seam_open_edges(&Waves, Extractor::MarchingCubes),
        0,
        "marching cubes stopped tiling across a chunk boundary"
    );
    assert_eq!(
        seam_open_edges(&Waves, Extractor::SurfaceNets),
        5,
        "surface nets' chunk-seam gap changed"
    );
    assert_eq!(
        seam_open_edges(&Waves, Extractor::DualContouring),
        4,
        "dual contouring's chunk-seam gap changed"
    );
    // B-007: subgrid tiles, and its zero is pinned alongside the others so it
    // cannot quietly stop.
    assert_eq!(
        seam_open_edges(&Waves, Extractor::Subgrid { samples: 4 }),
        0,
        "subgrid stopped tiling across a chunk boundary"
    );
}

#[test]
fn chunk_seams_reports_what_was_measured() {
    // The predicate and the measurement have to agree, or the API is decoration.
    for (extractor, expected) in [
        (Extractor::MarchingCubes, ChunkSeams::Closed),
        (Extractor::SurfaceNets, ChunkSeams::Gapped),
        (Extractor::DualContouring, ChunkSeams::Gapped),
        (Extractor::Subgrid { samples: 4 }, ChunkSeams::Closed),
    ] {
        let open = seam_open_edges(&Waves, extractor);
        let says = extractor.chunk_seams();
        assert_eq!(says, expected, "{extractor:?} reported {says:?}");
        match says {
            ChunkSeams::Closed => assert_eq!(open, 0, "{extractor:?} says Closed with {open} open"),
            ChunkSeams::Gapped => assert!(open > 0, "{extractor:?} says Gapped with none open"),
            ChunkSeams::Unverified => {}
        }
    }
    assert_eq!(
        Extractor::Subgrid { samples: 4 }.chunk_seams(),
        ChunkSeams::Closed
    );
}

/// **The quickstart's own configuration meshes, and every one of its eight chunks
/// carries triangles (B-009).**
///
/// `examples/quickstart.rs` is the file a new user copies, so "it compiles" is
/// the wrong bar — a layout that puts the sphere outside the chunks it spawns
/// compiles perfectly and shows an empty window. That mistake was made and caught
/// while writing it: `ChunkLayout::new`'s first argument is **cells**, not
/// samples, so the original `new(16, 0.25, …)` gave each chunk a 4.0-unit span
/// for a radius-1 sphere and left seven of the eight empty.
///
/// This asserts the numbers the example actually ships: 16 cells at `1/16` of a
/// unit, so one chunk spans exactly 1.0 and the eight tile `[-1, 1]³` — the cube
/// `Sphere::canonical()` fills. **Every chunk must carry triangles**, which is
/// what says the sphere's surface crosses all eight rather than sitting in one.
#[test]
fn the_quickstart_layout_meshes_every_chunk() {
    let mut app = App::new();
    app.add_plugins((TaskPoolPlugin::default(), AssetPlugin::default()))
        .add_plugins(IsomeshPlugin);

    // Exactly what examples/quickstart.rs passes.
    let layout = ChunkLayout::new(16, 1.0 / 16.0, [0.0; 3]).expect("a valid layout");
    let volume = app
        .world_mut()
        .spawn(VoxelVolume::new(layout, Sphere::<f32>::canonical()))
        .id();
    for z in -1..1 {
        for y in -1..1 {
            for x in -1..1 {
                app.world_mut().spawn((
                    VoxelChunk {
                        id: ChunkId::new([x, y, z]),
                        volume,
                    },
                    NeedsRemesh,
                ));
            }
        }
    }

    drain_until(&mut app, |app| meshed(app) == 8);
    assert_eq!(meshed(&mut app), 8, "the quickstart's queue never drained");

    // Not just "a mesh exists" -- a mesh with geometry in it, per chunk. An empty
    // handle would satisfy the count above and show nothing on screen.
    let handles: Vec<_> = app
        .world_mut()
        .query::<&ChunkMesh>()
        .iter(app.world())
        .map(|m| m.0.clone())
        .collect();
    let meshes = app.world().resource::<Assets<Mesh>>();
    let counts: Vec<usize> = handles
        .iter()
        .map(|h| {
            meshes
                .get(h)
                .and_then(|m| m.indices().map(|i| i.len() / 3))
                .unwrap_or(0)
        })
        .collect();
    assert!(
        counts.iter().all(|&n| n > 0),
        "a quickstart chunk meshed to nothing, so the sphere does not cross all \
         eight: {counts:?}"
    );
    std::println!("measured: quickstart chunks meshed to {counts:?} triangles");
}

// ── `IsomeshSystems`, `replace_field`, `mark_changed_chunks` ────────────────
//
// Added for downstream plugins that edit a volume's field (bevy_elements' snow
// does) and need to order against the meshing systems and learn which chunks
// an edit invalidated.

use bevy_ecs::system::RunSystemOnce;

/// Four chunks in a row along x, none marked, over a small sphere that sits
/// wholly inside chunk 0: chunk `[x, 0, 0]` covers world x in `[2x, 2x + 2)`.
fn edit_app() -> (App, Entity) {
    let mut app = App::new();
    app.add_plugins((TaskPoolPlugin::default(), AssetPlugin::default()))
        .add_plugins(IsomeshPlugin);
    let layout = ChunkLayout::new(8, 0.25, [0.0; 3]).expect("a valid layout");
    let volume = app
        .world_mut()
        .spawn(VoxelVolume::new(
            layout,
            Sphere::<f32> {
                center: [1.0, 1.0, 1.0],
                radius: 0.5,
            },
        ))
        .id();
    for x in 0..4 {
        app.world_mut().spawn(VoxelChunk {
            id: ChunkId::new([x, 0, 0]),
            volume,
        });
    }
    (app, volume)
}

fn marked(app: &mut App) -> Vec<ChunkId> {
    let mut ids: Vec<ChunkId> = app
        .world_mut()
        .query_filtered::<&VoxelChunk, With<NeedsRemesh>>()
        .iter(app.world())
        .map(|c| c.id)
        .collect();
    ids.sort();
    ids
}

/// Swap `volume`'s field for `field` and run `mark_changed_chunks` over the
/// change, as a downstream system would.
fn edit(
    app: &mut App,
    volume: Entity,
    field: impl VolumeField,
) -> isomesh::chunk::dirty::EditReport {
    let old = app
        .world_mut()
        .get_mut::<VoxelVolume>(volume)
        .expect("the volume exists")
        .replace_field(field);
    app.world_mut()
        .run_system_once(
            move |mut commands: Commands,
                  volumes: Query<&VoxelVolume>,
                  chunks: Query<(Entity, &VoxelChunk)>| {
                let v = volumes.get(volume).expect("the volume exists");
                mark_changed_chunks(&mut commands, volume, v, &chunks, &*old, &*v.field)
            },
        )
        .expect("the system ran")
        .expect("the region is far below the sample-space bound")
}

#[test]
fn replace_field_returns_the_old_field_and_installs_the_new_one() {
    let layout = ChunkLayout::new(8, 0.25, [0.0; 3]).expect("a valid layout");
    let mut volume = VoxelVolume::new(layout, Sphere::<f32>::canonical());
    // The canonical sphere has radius 1 about the origin, so (0, 0, 0.5) is
    // inside it and outside the replacement's tiny radius.
    let probe = [0.0, 0.0, 0.5];
    let before = volume.field.sample(probe);
    assert!(before < 0.0, "the probe should be inside the old field");

    let old = volume.replace_field(Sphere::<f32> {
        center: [0.0; 3],
        radius: 0.1,
    });
    assert_eq!(old.sample(probe), before, "the old field came back intact");
    assert!(
        volume.field.sample(probe) > 0.0,
        "the volume now holds the new field"
    );
    // Layout and extractor are untouched.
    assert_eq!(volume.layout, layout);
    assert_eq!(volume.extractor, Extractor::default());
}

#[test]
fn mark_changed_chunks_marks_exactly_the_chunks_the_edit_touched() {
    let (mut app, volume) = edit_app();
    assert!(marked(&mut app).is_empty());

    // Grow the sphere in place. Still wholly inside chunk 0, so chunks 1..4
    // saw no change at all and must not be re-meshed.
    let report = edit(
        &mut app,
        volume,
        Sphere::<f32> {
            center: [1.0, 1.0, 1.0],
            radius: 0.6,
        },
    );
    assert_eq!(marked(&mut app), vec![ChunkId::new([0, 0, 0])]);
    assert_eq!(report.dirty_chunks, 1);
    assert!(report.output_changed_cells > 0);
}

#[test]
fn mark_changed_chunks_marks_the_chunk_a_surface_left_and_the_one_it_entered() {
    let (mut app, volume) = edit_app();
    // Move the sphere from chunk 0 into chunk 2: both gained or lost surface.
    edit(
        &mut app,
        volume,
        Sphere::<f32> {
            center: [5.0, 1.0, 1.0],
            radius: 0.5,
        },
    );
    assert_eq!(
        marked(&mut app),
        vec![ChunkId::new([0, 0, 0]), ChunkId::new([2, 0, 0])]
    );
}

#[test]
fn mark_changed_chunks_marks_nothing_for_an_identical_field() {
    let (mut app, volume) = edit_app();
    let report = edit(
        &mut app,
        volume,
        Sphere::<f32> {
            center: [1.0, 1.0, 1.0],
            radius: 0.5,
        },
    );
    assert!(marked(&mut app).is_empty());
    assert_eq!(report.dirty_chunks, 0);
    assert_eq!(report.value_changed_cells, 0);
}

/// `edit_app`'s sphere plus a plate 0.6 cells thick inside cell x = 16 (chunk
/// `[2, 0, 0]`). Every corner of that cell is outside, before and after.
struct SphereAndPlate;

impl Sdf for SphereAndPlate {
    type Scalar = f32;
    fn sample(&self, p: [f32; 3]) -> f32 {
        let sphere = Sphere::<f32> {
            center: [1.0, 1.0, 1.0],
            radius: 0.5,
        }
        .sample(p);
        sphere.min((p[0] - 4.125).abs() - 0.075)
    }
}

fn the_edit_app_with(extractor: Extractor) -> (App, Entity) {
    let (mut app, volume) = edit_app();
    app.world_mut()
        .get_mut::<VoxelVolume>(volume)
        .expect("the volume exists")
        .extractor = extractor;
    (app, volume)
}

/// Issue #11. Subgrid meshes a plate thinner than a cell that no corner sees, so
/// the chunk that holds it must be re-meshed when the plate appears, and again
/// when it goes.
#[test]
fn mark_changed_chunks_marks_a_subgrid_feature_thinner_than_a_cell() {
    let extractor = Extractor::Subgrid { samples: 4 };
    let (mut app, volume) = the_edit_app_with(extractor);
    let layout = app
        .world()
        .get::<VoxelVolume>(volume)
        .expect("the volume exists")
        .layout;
    let chunk_2 = ChunkId::new([2, 0, 0]);

    // Precondition: Subgrid really does mesh the plate in chunk 2.
    let built = extract_chunk(&SphereAndPlate, &layout, chunk_2, extractor).expect("extracts");
    assert!(built.triangle_count() > 0, "the plate is not in the mesh");

    edit(&mut app, volume, SphereAndPlate);
    assert!(
        marked(&mut app).contains(&chunk_2),
        "chunk 2's Subgrid mesh changed and the chunk was not marked"
    );

    // The mirror: from the plate back to the sphere alone.
    for chunk in marked(&mut app) {
        let entity = app
            .world_mut()
            .query::<(Entity, &VoxelChunk)>()
            .iter(app.world())
            .find(|(_, c)| c.id == chunk)
            .map(|(e, _)| e)
            .expect("the chunk exists");
        app.world_mut().entity_mut(entity).remove::<NeedsRemesh>();
    }
    assert!(marked(&mut app).is_empty());
    edit(
        &mut app,
        volume,
        Sphere::<f32> {
            center: [1.0, 1.0, 1.0],
            radius: 0.5,
        },
    );
    assert!(
        marked(&mut app).contains(&chunk_2),
        "the plate left chunk 2's mesh and the chunk was not marked"
    );
}

/// The rule belongs to the extractor: Marching Cubes cannot mesh that plate, so
/// its volume has nothing to re-mesh, and the corner rule says so.
#[test]
fn mark_changed_chunks_leaves_a_marching_cubes_volume_on_the_corner_rule() {
    let (mut app, volume) = edit_app();
    edit(&mut app, volume, SphereAndPlate);
    assert!(
        marked(&mut app).is_empty(),
        "marching cubes meshes nothing of the plate: {:?}",
        marked(&mut app)
    );
}

/// `edit`, over a box the caller chose.
fn edit_in(
    app: &mut App,
    volume: Entity,
    field: impl VolumeField,
    region: ([i64; 3], [i64; 3]),
) -> isomesh::chunk::dirty::EditReport {
    let old = app
        .world_mut()
        .get_mut::<VoxelVolume>(volume)
        .expect("the volume exists")
        .replace_field(field);
    app.world_mut()
        .run_system_once(
            move |mut commands: Commands,
                  volumes: Query<&VoxelVolume>,
                  chunks: Query<(Entity, &VoxelChunk)>| {
                let v = volumes.get(volume).expect("the volume exists");
                mark_changed_chunks_in(&mut commands, volume, v, &chunks, &*old, &*v.field, region)
            },
        )
        .expect("the system ran")
        .expect("the region is far below the sample-space bound")
}

/// The box is what is compared, and so what is paid for: a sphere moves from
/// chunk 0 to chunk 2, and a box over chunk 2's cells sees chunk 2 only.
#[test]
fn mark_changed_chunks_in_compares_only_its_box() {
    let moved = || Sphere::<f32> {
        center: [5.0, 1.0, 1.0],
        radius: 0.5,
    };
    let (mut app, volume) = edit_app();
    // Chunk 2 owns cells x in 16..24.
    let report = edit_in(&mut app, volume, moved(), ([16, 0, 0], [23, 7, 7]));
    assert_eq!(marked(&mut app), vec![ChunkId::new([2, 0, 0])]);
    assert_eq!(report.region_cells, 8 * 8 * 8, "the cost follows the box");

    // A box the edit never reached marks nothing.
    let (mut app, volume) = edit_app();
    let report = edit_in(&mut app, volume, moved(), ([24, 0, 0], [31, 7, 7]));
    assert!(marked(&mut app).is_empty());
    assert_eq!(report.dirty_chunks, 0);

    // An inverted box is nothing to compare.
    let (mut app, volume) = edit_app();
    let report = edit_in(&mut app, volume, moved(), ([5, 5, 5], [4, 5, 5]));
    assert!(marked(&mut app).is_empty());
    assert_eq!(report, isomesh::chunk::dirty::EditReport::default());
}

/// The same bound for a Subgrid volume, where it matters most: the plate in
/// chunk 2 is found by a box over chunk 2 alone, at 4 samples.
#[test]
fn mark_changed_chunks_in_bounds_a_subgrid_volumes_cost() {
    let (mut app, volume) = the_edit_app_with(Extractor::Subgrid { samples: 4 });
    let report = edit_in(&mut app, volume, SphereAndPlate, ([16, 0, 0], [23, 7, 7]));
    assert_eq!(marked(&mut app), vec![ChunkId::new([2, 0, 0])]);
    assert_eq!(report.region_cells, 8 * 8 * 8);
}

#[test]
fn mark_changed_chunks_leaves_another_volumes_chunks_alone() {
    let (mut app, volume) = edit_app();
    // A second volume with a chunk at the same id as the one being edited.
    let layout = ChunkLayout::new(8, 0.25, [0.0; 3]).expect("a valid layout");
    let other = app
        .world_mut()
        .spawn(VoxelVolume::new(layout, Sphere::<f32>::canonical()))
        .id();
    let bystander = app
        .world_mut()
        .spawn(VoxelChunk {
            id: ChunkId::new([0, 0, 0]),
            volume: other,
        })
        .id();

    edit(
        &mut app,
        volume,
        Sphere::<f32> {
            center: [1.0, 1.0, 1.0],
            radius: 0.6,
        },
    );
    assert_eq!(marked(&mut app), vec![ChunkId::new([0, 0, 0])]);
    assert!(
        app.world().get::<NeedsRemesh>(bystander).is_none(),
        "the other volume's chunk was marked"
    );
}

#[test]
fn mark_changed_chunks_on_a_volume_with_no_chunks_reports_nothing() {
    let mut app = App::new();
    app.add_plugins((TaskPoolPlugin::default(), AssetPlugin::default()))
        .add_plugins(IsomeshPlugin);
    let layout = ChunkLayout::new(8, 0.25, [0.0; 3]).expect("a valid layout");
    let volume = app
        .world_mut()
        .spawn(VoxelVolume::new(layout, Sphere::<f32>::canonical()))
        .id();
    let report = edit(&mut app, volume, Sphere::<f32>::default());
    assert_eq!(report, isomesh::chunk::dirty::EditReport::default());
}

/// Marks every chunk. Ordered by the caller relative to [`IsomeshSystems`].
fn mark_every_chunk(mut commands: Commands, chunks: Query<Entity, With<VoxelChunk>>) {
    for entity in &chunks {
        commands.entity(entity).insert(NeedsRemesh);
    }
}

#[test]
fn a_system_before_isomesh_systems_is_meshed_the_same_frame() {
    let (mut app, _) = edit_app();
    app.add_systems(Update, mark_every_chunk.before(IsomeshSystems));
    app.update();
    let stats = *app.world().resource::<MeshStats>();
    assert_eq!(stats.spawned, 4, "the marks were not seen this frame");
    assert!(marked(&mut app).is_empty(), "NeedsRemesh was not consumed");
}

#[test]
fn a_system_after_isomesh_systems_waits_for_the_next_frame() {
    // The companion: proves the test above is sensitive to which side of the
    // set the marking system is on, rather than passing whatever the order.
    let (mut app, _) = edit_app();
    app.add_systems(Update, mark_every_chunk.after(IsomeshSystems));
    app.update();
    let stats = *app.world().resource::<MeshStats>();
    assert_eq!(stats.spawned, 0, "the marks were seen the same frame");
    assert_eq!(marked(&mut app).len(), 4);
    app.update();
    let stats = *app.world().resource::<MeshStats>();
    assert_eq!(stats.spawned, 4, "and the next frame picks them up");
}

// ── Issue #10: an empty chunk is published, not rendered ────────────────────
//
// Sanctuary's wall: 102 of 856 chunks mesh to nothing, and each was a
// `contains no vertices` warning from Bevy's mesh allocator, because the empty
// asset was extracted to the render world. The chunk still finishes -- it keeps
// its `ChunkMesh` -- but its mesh stays in the main world.

use bevy_asset::{AssetId, RenderAssetUsages};

/// One published chunk: its id, its mesh's asset id, vertex count and usages.
type Published = (ChunkId, AssetId<Mesh>, usize, RenderAssetUsages);

/// Every chunk's published mesh, sorted by chunk id.
///
/// Returns asset *ids*, not handles: a test that kept a strong handle would
/// keep the asset alive and could not see it freed.
fn published(app: &mut App) -> Vec<Published> {
    let handles: Vec<(ChunkId, Handle<Mesh>)> = app
        .world_mut()
        .query::<(&VoxelChunk, &ChunkMesh)>()
        .iter(app.world())
        .map(|(chunk, mesh)| (chunk.id, mesh.0.clone()))
        .collect();
    let meshes = app.world().resource::<Assets<Mesh>>();
    let mut rows: Vec<Published> = handles
        .iter()
        .map(|(id, handle)| {
            let mesh = meshes
                .get(handle)
                .expect("a published handle names a live asset");
            (*id, handle.id(), mesh.count_vertices(), mesh.asset_usage)
        })
        .collect();
    rows.sort_by_key(|row| row.0);
    rows
}

/// [`edit_app`] with every chunk meshed once: chunk 0 holds the sphere, chunks
/// 1, 2 and 3 hold nothing.
fn meshed_edit_app() -> (App, Entity) {
    let (mut app, volume) = edit_app();
    let chunks: Vec<Entity> = app
        .world_mut()
        .query_filtered::<Entity, With<VoxelChunk>>()
        .iter(app.world())
        .collect();
    for chunk in chunks {
        app.world_mut().entity_mut(chunk).insert(NeedsRemesh);
    }
    drain_until(&mut app, |app| meshed(app) == 4);
    assert_eq!(meshed(&mut app), 4, "the queue never drained");
    (app, volume)
}

#[test]
fn an_empty_chunk_keeps_its_chunk_mesh_but_not_the_render_world() {
    let (mut app, _) = meshed_edit_app();
    let rows = published(&mut app);

    let ids: Vec<ChunkId> = rows.iter().map(|row| row.0).collect();
    assert_eq!(
        ids,
        (0..4).map(|x| ChunkId::new([x, 0, 0])).collect::<Vec<_>>(),
        "every finished chunk has a ChunkMesh, the empty ones included"
    );

    let (_, _, vertices, usage) = rows[0];
    assert!(vertices > 0, "chunk 0 holds the sphere");
    assert_eq!(
        usage,
        RenderAssetUsages::default(),
        "a chunk with a surface is still drawn and still readable"
    );
    for &(id, _, vertices, usage) in &rows[1..] {
        assert_eq!(vertices, 0, "{id:?} should be empty");
        assert_eq!(
            usage,
            RenderAssetUsages::MAIN_WORLD,
            "{id:?} has no vertices, so the render world must not get it: Bevy's \
             mesh allocator warns `contains no vertices` for every one it is handed"
        );
    }
}

#[test]
fn a_chunk_that_gains_or_loses_its_surface_gets_a_fresh_mesh_with_the_matching_usage() {
    let (mut app, volume) = meshed_edit_app();
    let before = published(&mut app);

    // Move the sphere from chunk 0 into chunk 2: 0 loses its surface, 2 gains one.
    edit(
        &mut app,
        volume,
        Sphere::<f32> {
            center: [5.0, 1.0, 1.0],
            radius: 0.5,
        },
    );
    drain_until(&mut app, |app| {
        let rows = published(app);
        rows[0].1 != before[0].1 && rows[2].1 != before[2].1
    });
    let after = published(&mut app);

    assert_ne!(after[0].1, before[0].1, "chunk 0 was never re-published");
    assert_ne!(after[2].1, before[2].1, "chunk 2 was never re-published");
    assert_eq!(
        (after[0].2, after[0].3),
        (0, RenderAssetUsages::MAIN_WORLD),
        "the chunk the surface left"
    );
    assert!(after[2].2 > 0, "the chunk the surface entered has vertices");
    assert_eq!(
        after[2].3,
        RenderAssetUsages::default(),
        "the chunk the surface entered"
    );
    assert_eq!(
        (after[1].1, after[3].1),
        (before[1].1, before[3].1),
        "chunks the edit never touched keep their assets"
    );

    // The replaced assets go once their last handle drops: the old ChunkMesh
    // components were the only holders.
    let replaced = [before[0].1, before[2].1];
    let live = |app: &App| {
        let meshes = app.world().resource::<Assets<Mesh>>();
        replaced
            .iter()
            .filter(|&&id| meshes.get(id).is_some())
            .count()
    };
    drain_until(&mut app, |app| live(app) == 0);
    assert_eq!(live(&app), 0, "a replaced mesh outlived its handle");
}
