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

#[test]
fn every_chunk_eventually_gets_a_mesh() {
    let mut app = app(0, 12);
    for _ in 0..200 {
        app.update();
        if meshed(&mut app) == 12 {
            break;
        }
    }
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
    // The acceptance property. One update spawns twelve extractions of a
    // deliberately slow field; if any of that ran on the main thread, the update
    // could not come back in a fraction of the time the extractions take.
    let mut app = app(4_000, 12);

    let spawn_frame = Instant::now();
    app.update();
    let spawn_cost = spawn_frame.elapsed();

    // Drain, and time the whole thing. The total is what the work costs.
    let drain = Instant::now();
    for _ in 0..500 {
        app.update();
        if meshed(&mut app) == 12 {
            break;
        }
    }
    let total = drain.elapsed() + spawn_cost;
    assert_eq!(meshed(&mut app), 12, "the queue never drained");

    assert!(
        spawn_cost * 4 < total,
        "the frame that queued the work took {spawn_cost:?} of a total {total:?}, \
         which is not the profile of work that happened elsewhere"
    );
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

    let mut frames = 0;
    let mut worst = 0usize;
    for _ in 0..200 {
        let before = meshed(&mut app);
        app.update();
        let after = meshed(&mut app);
        worst = worst.max(after - before);
        frames += 1;
        if after == 12 {
            break;
        }
    }
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
    for _ in 0..200 {
        let before = meshed(&mut app);
        app.update();
        best = best.max(meshed(&mut app) - before);
        if meshed(&mut app) == 12 {
            break;
        }
    }
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
    assert!(
        app.world().entity(chunk).contains::<MeshingTask>(),
        "the chunk should be extracting"
    );
    assert!(
        !app.world().entity(chunk).contains::<NeedsRemesh>(),
        "NeedsRemesh should come off at spawn"
    );

    // The edit arrives while the first extraction is still running.
    app.world_mut().entity_mut(chunk).insert(NeedsRemesh);
    for _ in 0..500 {
        app.update();
        if !app.world().entity(chunk).contains::<NeedsRemesh>()
            && !app.world().entity(chunk).contains::<MeshingTask>()
        {
            break;
        }
    }
    assert!(
        app.world().entity(chunk).contains::<ChunkMesh>(),
        "the chunk never finished"
    );
    assert!(
        !app.world().entity(chunk).contains::<NeedsRemesh>(),
        "the re-mark was never serviced"
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

    for _ in 0..500 {
        app.update();
        if meshed(&mut app) == 1 {
            break;
        }
    }
    assert_eq!(meshed(&mut app), 1);
}
