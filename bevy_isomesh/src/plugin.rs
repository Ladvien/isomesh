//! B-003 — the plugin, and meshing that does not run on the main thread.
//!
//! [`mesh`](crate::mesh) turns one extraction into one [`Mesh`](bevy_mesh::Mesh).
//! This turns a
//! field and a chunk layout into a stream of them, on the task pool, under a
//! frame budget.
//!
//! # Where the work happens
//!
//! Extraction runs on [`AsyncComputeTaskPool`](bevy_tasks::AsyncComputeTaskPool),
//! never in a system. That is the
//! whole point of the ticket: a 64³ chunk through Marching Cubes is a quarter of
//! a million field evaluations, and M-98 measured subgrid Marching Tetrahedra at
//! 22 ms for a single 33³ grid — a third of a frame for *one chunk*. Anything
//! that size on the main thread is a stall you can see.
//!
//! What remains on the main thread is inserting finished meshes as assets, and
//! that is what the budget bounds.
//!
//! # The budget, and why it takes a `Duration` here when G-006 refused to
//!
//! G-006's [`mesh_within_budget`](isomesh::chunk::dirty::DirtySet::mesh_within_budget)
//! takes a **predicate**, not milliseconds, because `core` has no clock and a
//! `std` feature would have meant two paths through the core crate. Its doc
//! comment says what to do instead: *"a caller with a `std` clock writes
//! `|| start.elapsed() < budget`"*.
//!
//! This crate is that caller. It has `std`, it has a clock, and taking a
//! `Duration` here is where the decision always pointed — the constraint was
//! never "durations are wrong", it was "`no_std` cannot honour one".
//!
//! # It always applies at least one mesh
//!
//! The budget is consulted **after** each mesh is applied, never before, for the
//! same reason G-006 gives: a budget too small for a single chunk would
//! otherwise drain nothing forever while the queue grew, which is a livelock
//! that looks like a leak. Overshooting by at most one mesh is the price.
//!
//! # What this does not do
//!
//! It produces a [`Handle`](bevy_asset::Handle)`<Mesh>` on [`ChunkMesh`] and
//! stops there. Attaching a
//! render component is the application's job, because `Mesh3d` lives behind
//! `bevy_render` and this crate is built on leaf crates so that a consumer doing
//! CPU-only meshing never compiles the renderer. One line in the consumer buys
//! that, and the examples show it.

use std::sync::Arc;
use std::time::{Duration, Instant};

use bevy_app::{App, Plugin, Update};
use bevy_asset::{Assets, Handle};
use bevy_ecs::prelude::*;
use bevy_mesh::Mesh;
use bevy_tasks::{AsyncComputeTaskPool, Task, block_on, futures_lite::future};

use isomesh::chunk::{ChunkId, ChunkLayout};
use isomesh::{RuntimeShape3, Sdf};

use crate::MeshBuilder;

/// A field a [`VoxelVolume`] can be meshed from.
///
/// The bounds are what crossing a thread boundary costs: the task pool needs
/// `Send + Sync + 'static`, and `Sdf` is dyn-compatible once its scalar is
/// named, so a volume can hold a runtime-chosen field without the plugin being
/// generic over it.
pub trait VolumeField: Sdf<Scalar = f32> + Send + Sync + 'static {}

impl<T> VolumeField for T where T: Sdf<Scalar = f32> + Send + Sync + 'static {}

/// How a volume's chunks are extracted.
///
/// A plain enum rather than a boxed strategy: the set is closed, the choice has
/// to cross to a worker thread, and naming it in a component keeps the decision
/// visible in the ECS rather than hidden in a closure.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum Extractor {
    /// [`MarchingCubes`](isomesh::marching_cubes::MarchingCubes).
    #[default]
    MarchingCubes,
    /// [`SurfaceNets`](isomesh::surface_nets::SurfaceNets).
    SurfaceNets,
    /// [`DualContouring`](isomesh::dual_contouring::DualContouring).
    DualContouring,
    /// [`SubgridMarchingTetrahedra`](isomesh::subgrid::extract::SubgridMarchingTetrahedra),
    /// at the given 1D sampling resolution.
    ///
    /// The resolution is part of the choice because it decides *which features
    /// exist* rather than how well they are approximated — M-95 measured that
    /// changing it leaves the topology alone and moves positions by `~1e-12`,
    /// and M-98 that it is the whole of the cost.
    Subgrid {
        /// 1D samples per tetrahedron edge.
        samples: u32,
    },
}

/// Whether independently meshed chunks of an [`Extractor`] meet at their shared
/// face.
///
/// This exists because the plugin cannot mesh a volume any way *but* chunked,
/// and choosing the wrong extractor for that produces a world full of gaps with
/// nothing said about it. It reports rather than refuses, for the same reason
/// [`collider::readiness`](isomesh::collider::readiness) and
/// [`BrushOp::commutes_with`](isomesh::brush::BrushOp::commutes_with) report: a
/// consumer whose gaps are hidden, or who meshes one chunk, is entitled to make
/// that call.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ChunkSeams {
    /// Neighbouring chunks close. Vertices sit on grid **edges**, which two
    /// chunks compute from identical corner values and therefore agree on.
    Closed,
    /// Neighbouring chunks are **not guaranteed** to close, and are measured
    /// open on at least one field. One vertex per cell **interior** means a
    /// boundary quad needs the neighbour's vertex, which the chunk does not
    /// have, so it stops short.
    ///
    /// Not *always* open: whether a given surface leaves a gap depends on how it
    /// meets the seam, and both dual methods measure 0 on some fields. The
    /// guarantee is what is missing, not the geometry on any one frame.
    Gapped,
    /// Not established. Do not assume either way.
    Unverified,
}

impl Extractor {
    /// What this extractor does at a chunk boundary.
    ///
    /// Measured on two adjacent chunks, meshed independently, welded, counting
    /// boundary edges lying in the shared plane and excluding the two-chunk
    /// block's own six walls. Two fields, `waves` and `blobs`:
    ///
    /// | extractor | waves | blobs | |
    /// |---|---:|---:|---|
    /// | [`MarchingCubes`](Extractor::MarchingCubes) | 0 | 0 | [`Closed`](ChunkSeams::Closed) |
    /// | [`SurfaceNets`](Extractor::SurfaceNets) | 5 | 0 | [`Gapped`](ChunkSeams::Gapped) |
    /// | [`DualContouring`](Extractor::DualContouring) | 4 | 1 | [`Gapped`](ChunkSeams::Gapped) |
    /// | [`Subgrid`](Extractor::Subgrid) | 0 | 0 | [`Closed`](ChunkSeams::Closed) |
    ///
    /// `the_seam_counts_are_pinned` holds all four, so they cannot drift, and it
    /// derives the block's bounds from the layout rather than writing them out —
    /// **both hand-written versions were wrong**, one omitting `y` and one giving
    /// `z` a bound twice the chunk's depth, and each inflated a count (M-132).
    ///
    /// # Why `Subgrid` is closed
    ///
    /// Measured 0 across 20 configurations — two fields here, plus six field
    /// phases at three sampling resolutions each — and its edge orientation is a
    /// property of the grid rather than of a tetrahedron: `TETS[t]` is ordered by
    /// inclusion, so a tet edge always runs from the lower cube-corner index to
    /// the higher. **M-79's warning is about a different renumbering** — it says
    /// a mesh that renumbered vertices *per tet* would crack along every shared
    /// face, and chunking renumbers per chunk while leaving the within-cell
    /// corner order untouched.
    #[must_use]
    pub fn chunk_seams(self) -> ChunkSeams {
        match self {
            Self::MarchingCubes => ChunkSeams::Closed,
            Self::SurfaceNets | Self::DualContouring => ChunkSeams::Gapped,
            Self::Subgrid { .. } => ChunkSeams::Closed,
        }
    }
}

/// A field, a chunk layout, and how to mesh it.
#[derive(Component, Clone)]
pub struct VoxelVolume {
    /// Where the chunks are.
    pub layout: ChunkLayout<f32>,
    /// What to mesh.
    ///
    /// Behind an [`Arc`] because every queued chunk task needs its own handle to
    /// it and the tasks outlive the system that spawned them.
    pub field: Arc<dyn VolumeField>,
    /// Which extractor to run.
    pub extractor: Extractor,
}

impl VoxelVolume {
    /// A volume over `field`, meshed with Marching Cubes.
    pub fn new(layout: ChunkLayout<f32>, field: impl VolumeField) -> Self {
        Self {
            layout,
            field: Arc::new(field),
            extractor: Extractor::default(),
        }
    }

    /// The same volume, meshed with `extractor`.
    ///
    /// # Check [`Extractor::chunk_seams`] before choosing
    ///
    /// A volume is always meshed **chunked**, and the dual methods are not
    /// guaranteed to tile: [`SurfaceNets`](Extractor::SurfaceNets) and
    /// [`DualContouring`](Extractor::DualContouring) measure up to 5 open edges
    /// on a single seam, and 0 on other fields. That is structural -- a boundary
    /// quad needs the neighbour cell's vertex -- and not something a future
    /// release fixes. [`MarchingCubes`](Extractor::MarchingCubes) and
    /// [`Subgrid`](Extractor::Subgrid) measure 0 everywhere tried.
    ///
    /// It is not refused here, because a consumer meshing a single chunk, or one
    /// whose gaps are never seen, is entitled to the sharper extractor. It is
    /// said out loud because saying nothing put a cracked world in this
    /// repository's own README for a commit (M-128).
    #[must_use]
    pub fn with_extractor(mut self, extractor: Extractor) -> Self {
        self.extractor = extractor;
        self
    }
}

/// A chunk of a volume. Lives on its own entity, with the volume as its parent.
#[derive(Component, Clone, Copy, Debug, PartialEq, Eq)]
pub struct VoxelChunk {
    /// Which chunk of the volume this is.
    pub id: ChunkId,
    /// The volume entity it belongs to.
    pub volume: Entity,
}

/// Marks a chunk as needing to be meshed.
///
/// Removed when its task is *spawned*, not when the mesh arrives — so a chunk
/// edited while meshing is re-marked and re-queued, rather than having its edit
/// silently dropped into a task that started before it.
#[derive(Component, Clone, Copy, Debug, Default)]
pub struct NeedsRemesh;

/// The finished mesh for a chunk.
///
/// Attaching a render component to this is the application's job — see the
/// module docs.
#[derive(Component, Clone, Debug)]
pub struct ChunkMesh(pub Handle<Mesh>);

/// An extraction in flight on the task pool.
#[derive(Component)]
struct MeshingTask(Task<MeshBuilder>);

/// How much main-thread time per frame the plugin may spend turning finished
/// extractions into assets.
#[derive(Resource, Clone, Copy, Debug)]
pub struct MeshBudget {
    /// Wall-clock ceiling per frame. At least one mesh is always applied.
    pub per_frame: Duration,
    /// How many extractions may be in flight at once.
    ///
    /// Unbounded queueing is its own stall: every task holds a `MeshBuilder`,
    /// and a thousand 64³ chunks in flight is memory the frame budget never
    /// sees. Defaults to twice the pool's thread count.
    pub max_in_flight: usize,
}

impl Default for MeshBudget {
    fn default() -> Self {
        Self {
            // A quarter of a 60 Hz frame, leaving the rest for everything that
            // is not meshing.
            per_frame: Duration::from_micros(4_000),
            max_in_flight: 2 * bevy_tasks::available_parallelism(),
        }
    }
}

/// What the plugin did last frame.
#[derive(Resource, Clone, Copy, Debug, Default)]
pub struct MeshStats {
    /// Extractions spawned onto the pool.
    pub spawned: usize,
    /// Finished meshes turned into assets.
    pub applied: usize,
    /// Extractions still running.
    pub in_flight: usize,
    /// Chunks still waiting for a task slot.
    pub waiting: usize,
}

/// Chunks marked for meshing that are not already meshing.
///
/// A named alias because the tuple is at clippy's complexity ceiling and the
/// filter is the interesting part: `Without<MeshingTask>` is what stops a chunk
/// re-marked mid-extraction from being queued twice.
type WaitingChunks<'w, 's> =
    Query<'w, 's, (Entity, &'static VoxelChunk), (With<NeedsRemesh>, Without<MeshingTask>)>;

/// Registers the meshing systems.
#[derive(Clone, Copy, Debug, Default)]
pub struct IsomeshPlugin;

impl Plugin for IsomeshPlugin {
    fn build(&self, app: &mut App) {
        // The plugin inserts `Mesh` assets, so it needs `Assets<Mesh>` to
        // exist. `MeshPlugin` is where that normally comes from, and adding it
        // here rather than requiring it means the plugin works in a headless
        // app with no renderer -- which is what its own tests run in, and what
        // a server meshing for collision would want.
        if !app.is_plugin_added::<bevy_mesh::MeshPlugin>() {
            app.add_plugins(bevy_mesh::MeshPlugin);
        }
        app.init_resource::<MeshBudget>()
            .init_resource::<MeshStats>()
            // `spawn` then `apply`, chained: a task spawned this frame cannot
            // have finished, so running them the other way round would simply
            // add a frame of latency to every chunk.
            .add_systems(Update, (spawn_meshing_tasks, apply_finished_meshes).chain());
    }
}

/// Queue extractions for chunks that need them.
fn spawn_meshing_tasks(
    mut commands: Commands,
    budget: Res<MeshBudget>,
    mut stats: ResMut<MeshStats>,
    volumes: Query<&VoxelVolume>,
    waiting: WaitingChunks,
    in_flight: Query<(), With<MeshingTask>>,
) {
    let pool = AsyncComputeTaskPool::get();
    let mut running = in_flight.iter().count();
    stats.spawned = 0;

    let mut pending = 0usize;
    for (entity, chunk) in &waiting {
        if running >= budget.max_in_flight {
            pending += 1;
            continue;
        }
        let Ok(volume) = volumes.get(chunk.volume) else {
            // The volume is gone; so is the reason to mesh its chunk.
            commands.entity(entity).remove::<NeedsRemesh>();
            continue;
        };

        let field = Arc::clone(&volume.field);
        let layout = volume.layout;
        let extractor = volume.extractor;
        let id = chunk.id;
        let task = pool.spawn(async move { extract_chunk(&*field, &layout, id, extractor) });

        // `NeedsRemesh` comes off now rather than when the mesh lands, so an
        // edit arriving mid-extraction re-marks the chunk and is re-queued
        // instead of being swallowed by a task that started before it.
        commands
            .entity(entity)
            .remove::<NeedsRemesh>()
            .insert(MeshingTask(task));
        running += 1;
        stats.spawned += 1;
    }

    stats.in_flight = running;
    stats.waiting = pending;
}

/// Turn finished extractions into assets, within the frame budget.
fn apply_finished_meshes(
    mut commands: Commands,
    budget: Res<MeshBudget>,
    mut stats: ResMut<MeshStats>,
    mut meshes: ResMut<Assets<Mesh>>,
    mut tasks: Query<(Entity, &mut MeshingTask)>,
) {
    let started = Instant::now();
    stats.applied = 0;

    for (entity, mut task) in &mut tasks {
        let Some(builder) = block_on(future::poll_once(&mut task.0)) else {
            continue;
        };
        commands
            .entity(entity)
            .remove::<MeshingTask>()
            .insert(ChunkMesh(meshes.add(builder.into_mesh())));
        stats.applied += 1;

        // Checked *after* the work, so a budget too small for one mesh still
        // makes progress. See the module docs.
        if started.elapsed() >= budget.per_frame {
            break;
        }
    }

    stats.in_flight = stats.in_flight.saturating_sub(stats.applied);
}

/// Extract one chunk. Runs on the task pool, so it touches no ECS state.
fn extract_chunk(
    field: &dyn VolumeField,
    layout: &ChunkLayout<f32>,
    id: ChunkId,
    extractor: Extractor,
) -> MeshBuilder {
    let mut out = MeshBuilder::new();
    let Ok(shape) = layout.sample_shape() else {
        return out;
    };
    let origin = layout.sample_origin(id);
    let cell = layout.cell_size();
    // A failed extraction yields an empty chunk rather than a panic: this runs
    // on a worker thread, where a panic is a silent lost task rather than a
    // crash, which is the worst of both.
    let _ = extract_with(extractor, field, &shape, origin, cell, &mut out);
    out
}

fn extract_with(
    extractor: Extractor,
    field: &dyn VolumeField,
    shape: &RuntimeShape3,
    origin: [f32; 3],
    cell: f32,
    out: &mut MeshBuilder,
) -> isomesh::Result<()> {
    match extractor {
        Extractor::MarchingCubes => isomesh::marching_cubes::MarchingCubes::<f32>::new()
            .extract(&field, shape, origin, cell, out),
        Extractor::SurfaceNets => isomesh::surface_nets::SurfaceNets::<f32>::new()
            .extract(&field, shape, origin, cell, out),
        Extractor::DualContouring => isomesh::dual_contouring::DualContouring::<f32>::new()
            .extract(&field, shape, origin, cell, out),
        Extractor::Subgrid { samples } => {
            isomesh::subgrid::extract::SubgridMarchingTetrahedra::<f32>::new(samples)?
                .extract(&field, shape, origin, cell, out)
        }
    }
}

#[cfg(test)]
mod tests;
