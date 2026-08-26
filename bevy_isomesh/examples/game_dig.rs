//! E-202 — carving tunnels, the way a game does it.
//!
//! ```bash
//! cd bevy_isomesh && cargo run --example game_dig --release
//! ```
//!
//! You start **walking** on the rock: `WASD` moves, the mouse looks, `Space`
//! jumps, gravity holds you down and the field stops you walking through stone.
//! `F` switches to the old fly mode, where `Q`/`E` go up and down and nothing
//! falls.
//!
//! **The numbers panel starts hidden**, because this is a game first and the
//! panel sits on top of the rock you are digging. A coloured headline names the
//! mesher — a different colour for each of the eight — and one line of keys sits
//! under it; `H` brings the numbers back, and `ISOMESH_VIEW=hud` opens with them
//! for a capture.
//!
//! A translucent sphere on the rock under the crosshair is the brush that a click
//! would push — orange to carve, cyan to fill. **Hold** the left button to keep
//! carving along the sweep, the right to keep filling; edits are paced at
//! **12.5 a second**, which is below the rate a hand can tell apart and well
//! below the rate that used to bury the frame. The wheel or `[`/`]` resizes the
//! brush, `1`–`8` swap the mesher, `Z` undoes one brush, `X` clears the log, `C`
//! outlines the chunks the last edit re-meshed, and `Tab` releases the cursor.
//!
//! # The rock is textured by the shader, not by the mesh
//!
//! An isosurface has no natural parameterisation, so `MeshBuilder`'s UVs are a
//! dominant-axis planar projection that seams wherever the dominant axis flips —
//! its own doc comment says so. This example never reads them: the terrain is an
//! `ExtendedMaterial` whose fragment shader samples all three planes and
//! interpolates by the normal, which is a function of world position alone and so
//! cannot see a chunk boundary. It samples a **four-layer texture array** and
//! blends grass, leafy surface dirt and deep dirt by slope and world height, so
//! a fresh tunnel reads as its own depth; the fourth layer is concrete and the
//! five slabs lining the sandbox force it. The shader and both `512x2048` packed
//! arrays are **compiled in** (`include_bytes!`, `load_internal_asset!`),
//! because nothing copies an `assets/` tree into `web/dist` and a run-time load
//! path would work natively and 404 in the browser.
//!
//! # What this is actually testing
//!
//! Not "can it mesh a field" — three examples already do that. This is the first
//! one where the mesh is **re-built while someone is holding the mouse down**,
//! and it exists to put two numbers on screen that nothing else in the repo
//! measures under load:
//!
//! - **Chunks touched per edit.** A brush changes the field everywhere, because
//!   an SDF is global; what it changes *visibly* is a shell. G-002's `mark_edit`
//!   compares the field either side of one edit and marks only the chunks whose
//!   cells actually moved.
//! - **E1 — the fraction of the brush's own bounding box that changed.** M-33
//!   measured 15–36% offline. This shows it live, per edit, and it is the number
//!   the entire incremental story rests on: if it were 100%, re-meshing the
//!   bounding box would be as cheap as being clever about it.
//!
//! # The spacing is a power of two, deliberately
//!
//! `h = 0.125`. **M-32** measured that two chunks agree on their shared sample
//! plane bit-for-bit only at a power-of-two cell size; anywhere else they differ
//! by an ulp and the seam needs A-013's weld to close. This example does not
//! weld — each chunk is its own `Mesh3d`, exactly as an engine would keep them —
//! so it uses the spacing where the seam is exact and the surface is continuous
//! without one. `chunk_seam_weld` is the example that shows the other case.
//!
//! # The edit log grows, and the cost grows with it — sub-linearly
//!
//! Edits compose rather than mutate: the field is a `BrushStack` over the base
//! terrain, and carving pushes a brush. That is what makes undo a re-fold of the
//! log rather than a snapshot (E-207's premise), and it means **every field
//! sample walks every brush**.
//!
//! So the cost grows, and it is worth being precise about how much rather than
//! waving at it. Measured over a 60-carve scripted run (`ISOMESH_AUTOCARVE=60`,
//! which prints one line per edit), median milliseconds per re-meshed chunk:
//!
//! | edits in the log | 1–15 | 16–30 | 31–45 | 46–60 |
//! |---|---|---|---|---|
//! | ms per chunk | 0.158 | 0.354 | 0.525 | 0.589 |
//!
//! **3.7× for 7× the log, and flattening** — not proportional, even though every
//! sample really does walk every brush. So the stack walk is a real cost and not
//! the dominant one at these lengths; what else is in there has not been
//! measured and is not asserted here. Press `X` to clear the log and watch it
//! drop back.

mod common;

use std::cell::Cell;
use std::collections::HashMap;
use std::time::Duration;

use bevy::asset::{RenderAssetUsages, load_internal_asset, uuid_handle};
use bevy::image::{
    CompressedImageFormats, ImageAddressMode, ImageSampler, ImageSamplerDescriptor, ImageType,
};
use bevy::input::mouse::{AccumulatedMouseMotion, AccumulatedMouseScroll, MouseScrollUnit};
use bevy::pbr::{ExtendedMaterial, MaterialExtension, MaterialPlugin};
use bevy::platform::time::Instant;
use bevy::prelude::*;
use bevy::render::mesh::{Indices, PrimitiveTopology};
use bevy::render::render_resource::AsBindGroup;
use bevy::render::renderer::{RenderDevice, RenderQueue};
use bevy::shader::{Shader, ShaderRef};
use bevy::window::{CursorGrabMode, CursorOptions, PrimaryWindow};
use bevy_isomesh::MeshBuilder;
use common::{Capture, CommonPlugin, DemoStats, OrbitCamera, ViewFlags};
use isomesh::brush::{Brush, BrushOp, BrushStack};
use isomesh::chunk::dirty::{DirtySet, EditReport, mark_edit};
use isomesh::chunk::{ChunkId, ChunkLayout};
use isomesh::dual_contouring::DualContouring;
use isomesh::fields::{BoundedSdf, FieldBound, Sphere};
use isomesh::greedy_quads::GreedyQuads;
use isomesh::manifold_dual_contouring::ManifoldDualContouring;
use isomesh::marching_cubes::{FaceAmbiguity, MarchingCubes};
use isomesh::marching_tetrahedra::MarchingTetrahedra;
use isomesh::surface_nets::SurfaceNets;
use isomesh::{Real, RuntimeShape3, Sdf};
use isomesh_gpu::{
    ExtractTimings, FieldBuffer, FieldSampler, GpuBrush, GpuOp, GpuShape, GridParams,
    MarchingCubesGpu, Readback, read_bytes_many_deferred,
};

/// Chunk edge, in cells.
const CHUNK_CELLS: u32 = 16;
/// See the module docs: a power of two, so the seam is bit-exact without a weld.
const CELL_SIZE: f32 = 0.125;
/// Chunks along x and z, and up in y.
///
/// A chunk spans `CHUNK_CELLS * CELL_SIZE = 16 * 0.125 = 2.0` units, so this is
/// a **16x8x16-unit sandbox of 256 chunks**. Still small enough to mesh at
/// startup in one go, and worth saying in numbers: 256 chunks at the ~0.16 ms
/// per chunk the table above measures is ~40 ms of Marching Cubes, and only the
/// one or two chunk layers straddling the surface produce a mesh at all --
/// `rebuild` skips a chunk with no crossing rather than spawning an empty draw
/// call over it.
const EXTENT: [i32; 3] = [8, 4, 8];

/// Thickness of the five cuboids that line the sandbox.
///
/// Each slab occupies exactly the thickness immediately *outside* one face, so
/// its inner surface is the boundary plane the walk clamp enforces and no part
/// of it is inside the box. Outside rather than coplanar because `Ground` is
/// solid right up to the wall, and a slab sharing a plane with terrain z-fights
/// it along the whole seam. It also puts the geometry where `aim` refuses to
/// carve -- outside `sandbox` -- which is why the walls can never be dug
/// through.
const WALL_THICKNESS: f32 = 0.5;

/// Gizmos for the re-meshed-chunk outline.
#[derive(Default, Reflect, GizmoConfigGroup)]
struct ChunkGizmos;

/// The terrain before any edit: a slab with a rolling top.
///
/// Hand-rolled rather than `FbmTerrain`, because this needs a floor a player can
/// stand on and a ceiling to dig into, and it must be cheap — it is sampled
/// inside the edit loop.
///
/// **A pure function of position, and that is load-bearing.** [`GpuFields`]
/// samples this once per chunk and folds every later edit over the cached
/// buffer on the device, so the cache has no invalidation because nothing can
/// invalidate it. Give this state -- an animation clock, a seed that is edited,
/// a dependence on the camera -- and the GPU path goes silently stale.
#[derive(Clone, Copy)]
struct Ground;

impl Sdf for Ground {
    type Scalar = f32;

    fn sample(&self, p: [f32; 3]) -> f32 {
        // Distance to a wavy height field, negative below it.
        let height = 0.35 * (p[0] * 0.9).sin() * (p[2] * 0.7).cos() + 0.15 * (p[0] * 2.1).sin();
        p[1] - height
    }
}

/// Lipschitz bound for the whole editable field, which is what makes sphere
/// tracing sound here.
///
/// `Ground` is **not** 1-Lipschitz. With
/// `height = 0.35·sin(0.9x)·cos(0.7z) + 0.15·sin(2.1x)`,
/// `|∂height/∂x| ≤ 0.35·0.9 + 0.15·2.1 = 0.63` and `|∂height/∂z| ≤ 0.35·0.7 =
/// 0.245`, so for `f = y − height`, `|∇f| ≤ √(1 + 0.63² + 0.245²) = 1.207`.
/// `min`/`max` of Lipschitz functions keeps the constant and the brushes are
/// exact spheres, so this one number bounds the whole `BrushStack` no matter how
/// long the edit log grows -- the property `FieldBound::Lipschitz` documents and
/// M-354 exploits. Stepping by `|f|` rather than `|f| / L` would overshoot a
/// slope and tunnel through the surface.
const LIPSCHITZ: f32 = 1.25;

/// `Ground` is not a distance, so it must say what its values are worth before a
/// pruner can bound it over a box.
///
/// [`LIPSCHITZ`] is the number, derived above: `|grad f| <= 1.207` for this
/// height field, and `1.25` is the constant the whole [`BrushStack`] is bounded
/// by, because min/max of Lipschitz functions preserves it and every brush is an
/// exact sphere.
impl BoundedSdf for Ground {
    fn value_bound(&self) -> FieldBound {
        FieldBound::Lipschitz {
            l: LIPSCHITZ as f64,
        }
    }
}

/// ULP of slack added to every enclosure bound.
///
/// Copied with the pruner it belongs to, from `tape_pruning.rs`: the evaluation
/// error of a sphere distance is a few ULP of a magnitude near this sandbox's
/// far corner, and 64 ULP of that is far below the Lipschitz reach it is added
/// to. It buys a wide margin and costs a pruning decision only exactly on the
/// boundary -- where keeping the brush is the safe answer anyway.
const PAD_ULPS: f32 = 64.0;
/// Nearest a brush may be placed, so you cannot dig inside your own eye.
const AIM_NEAR: f32 = 0.30;
/// Furthest the trace looks: the sandbox's full 16x8x16 diagonal is
/// `sqrt(256 + 64 + 256) = 24`, and the camera can stand at any corner of it, so
/// anything shorter would leave the far corner unreachable.
const AIM_FAR: f32 = 25.0;
/// Iteration cap. A trace that has not converged in this many steps is grazing
/// the surface almost tangentially, where the answer is not useful anyway.
const AIM_STEPS: u32 = 128;
/// Surface tolerance, and the minimum step, which is what stops a stall at a
/// near-zero sample.
const AIM_HIT: f32 = 0.01;

/// Extraction time one frame may spend, and nothing else.
///
/// The same quarter-of-a-60-Hz-frame `MeshBudget::default` uses, which leaves
/// 12 ms for everything that is not meshing. `spend` is consulted *after* each
/// chunk, so a budget too small for one chunk still makes progress -- the
/// livelock `DirtySet::mesh_within_budget`'s docs warn about cannot happen.
///
/// `Duration` is a pure type with no clock, so `std` is right here; the *clock*
/// is `bevy::platform::time::Instant`, which is `std` natively and
/// `performance.now()` in a browser.
const MESH_BUDGET: Duration = Duration::from_micros(4_000);

/// Triangles one GPU chunk may emit before it is truncated.
///
/// `extract_indirect` sizes its buffers from a budget rather than reading the
/// scan total back, which is exactly what makes it the only non-blocking entry
/// point -- knowing the real count would cost the wait this design exists to
/// avoid. 8,192 is an order of magnitude above the ~700 triangles a chunk of
/// this sandbox produces under Marching Cubes, and costs `8192 * 9 * 4` = 295 KB
/// per buffer.
const GPU_TRIANGLE_BUDGET: u32 = 8_192;

/// GPU chunks allowed in flight at once.
///
/// Backpressure, not a second budget: each job holds two 295 KB geometry buffers
/// plus a staging copy, so an unbounded queue against a slow adapter is a
/// megabyte a chunk with nothing to stop it. When the list is full the drain
/// stops and its chunks stay in the dirty set, so the work is deferred rather
/// than dropped -- there is one queue and it is still the only one.
const GPU_JOBS_MAX: usize = 16;

// ── the body ────────────────────────────────────────────────────────────────
//
// Collision samples the field directly. No example in this directory does that
// -- `game_walk` and `game_capsule_walk` cast against a parry3d `TriMesh` built
// from the chunk meshes -- and the reason is the one `trace`'s doc comment
// already gives: this demo edits the field on every frame of a held button, so a
// collider cache is invalid before it is built. M-116 puts a convex
// decomposition at 241-272 ms per fragment and M-135 puts the collider check at
// 45% of a usable mesh, the largest single stage. Sampling the field is also the
// only option that works on both targets, and `crates/isomesh` offers no
// raycast, sphere cast or closest-point query to build on.

/// Radius of each of the four spheres the body is made of, so the body is
/// `0.50` wide.
///
/// Was `0.4`, which made a body `0.8` wide -- human height and 2.7 times human
/// width, and wider than the `0.50` cavity the default brush (radius `0.25`)
/// carves, so the player could not follow their own tunnel.
const BODY_RADIUS: f32 = 0.25;
/// Sphere centres below the eye.
///
/// Four spheres at `0.40` spacing against a `0.50` diameter, so consecutive
/// spheres **overlap** and the body is a continuous capsule. The old pair sat
/// exactly `2 * BODY_RADIUS` apart and touched at a single point, which is a
/// pinched waist -- a lip of rock could pass between them and the resolver
/// would see nothing to push out of.
///
/// The eye sits at the top of the topmost sphere and the lowest one rests on the
/// ground, so standing on flat terrain puts the eye at `1.45 + 0.25 = 1.70`:
/// a person, and where `setup` places the camera.
const BODY_OFFSETS: [f32; 4] = [0.25, 0.65, 1.05, 1.45];
/// Downward acceleration. Roughly twice Earth's, which is the usual game figure:
/// real gravity over a 1.7-unit body reads as floating.
const GRAVITY: f32 = 18.0;
/// Launch speed. `8.5^2 / (2 * 18)` is a **2.0-unit apex**, which is a whole
/// chunk: high enough to jump onto the lip of a pit dug with a brush at the
/// large end of the wheel's range, rather than the 1.0 the first version had.
const JUMP_SPEED: f32 = 8.5;
/// Resolution passes per frame. Three, because a body wedged in a corner is
/// pushed out of one sphere into the next and needs another look for each -- and
/// the body went from two spheres to four. Twelve sample-and-gradient pairs a
/// frame, against the eight the two-sphere pair cost at two passes.
const RESOLVE_PASSES: u32 = 3;
/// How far below the lowest sphere to look for ground.
const GROUND_PROBE: f32 = 0.06;
/// Below this, the gradient carries no direction. See [`resolve_body`].
const GRADIENT_EPS: f32 = 1e-4;

/// Where the next edit lands, recomputed every frame and read by both the ghost
/// and the edit. One point, two consumers: a preview that can disagree with the
/// brush it previews is worse than no preview.
#[derive(Resource, Default)]
struct Aim {
    point: Vec3,
    /// `false` when the ray leaves the sandbox without crossing the surface.
    hit: bool,
}

/// Seconds between edits while a button is held.
///
/// **12.5 edits a second.** Was 20 (`0.05`), and this is the plan's own next
/// step for the case a hold still costs too much: every brush is a term in the
/// `(L + 1)` factor on every field sample the demo takes afterwards, so the
/// cheapest way to make a long stroke cheap is to lay down fewer of them. 80 ms
/// is still under the ~100 ms at which a hand reads its own action as
/// instantaneous, and a sweep still cuts one continuous tunnel because the
/// spheres overlap at half a radius.
///
/// The distance gate stays: it is what stops a stationary brush pushing an
/// idempotent duplicate, and this bounds the case the distance gate cannot -- a
/// fast sweep, where the aim point clears half a radius every frame.
const EDIT_PERIOD: f32 = 0.08;

/// Marker for the translucent brush preview.
#[derive(Component)]
struct Ghost;

/// The dim-and-wait overlay, shown while an extractor switch drains.
///
/// There is no full-screen overlay anywhere else under `examples/`, so this is
/// new furniture rather than a reuse: the tree's banner rows are bounded strips
/// at `GlobalZIndex(4)`, and this has to cover the HUD.
#[derive(Component)]
struct LoadingModal;

/// The two ghost skins, built once. Carve is warm, fill is cold, and the colour
/// is the only thing that says which button is about to do what.
#[derive(Resource)]
struct GhostMaterials {
    carve: Handle<StandardMaterial>,
    fill: Handle<StandardMaterial>,
}

/// Gizmo group for the ghost's outline.
#[derive(Default, Reflect, GizmoConfigGroup)]
struct GhostGizmos;

#[derive(Resource)]
struct World {
    layout: ChunkLayout<f32>,
    brushes: Vec<Brush<Sphere<f32>>>,
    /// Scratch for [`prune_into`]: the brushes that can still change the fold
    /// over the chunk being meshed. One `Vec`, reused across every chunk of
    /// every frame, so it stops reallocating after the first few.
    survivors: Vec<Brush<Sphere<f32>>>,
    /// The one re-mesh queue. `setup`, `dig` and `switch_algorithm` only ever
    /// insert into it; [`drain_dirty`] is the only thing that meshes.
    dirty: DirtySet,
    /// Chunks the last budgeted drain could not reach.
    backlog: usize,
    /// Whether the current backlog came from an extractor switch, which is what
    /// the modal reports. A dig's one-to-eight chunks drain in a single frame and
    /// must not flash a modal.
    switching: bool,
    radius: f32,
    /// Chunks re-meshed by the most recent edit, for the outline.
    last_touched: Vec<ChunkId>,
    last_edit: Option<EditReport>,
    last_edit_ms: f64,
    last_chunks: usize,
    show_chunks: bool,
    grabbed: bool,
    /// Which extractor meshes a chunk, and what the last switch to it cost.
    algorithm: Algorithm,
    switch_ms: f64,
    switch_chunks: usize,
    /// When the switch being drained started, so its cost is wall-clock time
    /// across frames rather than one blocking call.
    switch_started: Instant,
    /// Centre of the last brush pushed by a held button, or `None` when no
    /// button is down. This is what makes a hold a *stroke* rather than a burst.
    stroke_last: Option<Vec3>,
    stroke_edits: u32,
    stroke_ms: f64,
    /// Seconds since the last edit, against [`EDIT_PERIOD`]. This is what makes
    /// a hold a fixed *rate* rather than one edit per frame.
    stroke_clock: f32,
    /// Velocity of the walking body. Only `y` is integrated -- horizontal motion
    /// is direct, as it is in every first-person game that is not a vehicle.
    velocity: Vec3,
    /// `true` in walk mode, `false` in fly mode. `F` toggles.
    walking: bool,
    /// Whether the body was standing on something at the end of the last
    /// resolve, which is what makes a jump possible.
    grounded: bool,
    /// Where the last GPU chunk's time went, or `None` before the first one.
    gpu_cost: Option<GpuCost>,
    /// Largest **unclamped** triangle count any GPU chunk has reported since the
    /// last switch.
    ///
    /// The peak rather than the last, because most of the 256 chunks are empty
    /// air and the last one collected is usually one of them -- a HUD reading
    /// `0 tris` after meshing the whole sandbox says nothing. And the peak is the
    /// figure that matters: this pinning at exactly [`GPU_TRIANGLE_BUDGET`] is
    /// the only tell that a chunk was truncated.
    gpu_peak_triangles: u32,
    /// Readbacks the device refused. Zero is the only acceptable value; it is on
    /// the HUD so a non-zero one cannot hide.
    gpu_failures: u32,
}

#[derive(Component)]
struct Chunk(ChunkId);

/// The GPU mesher, built once. Bevy's own device, reached through
/// `RenderDevice::wgpu_device()` -- `isomesh-gpu` has never heard of Bevy, and
/// `gpu_compute_mc.rs` is the example that exists to prove that abstraction does
/// not leak.
#[derive(Resource)]
struct GpuMesher {
    mc: MarchingCubesGpu,
}

/// Per-chunk base field samples, and the sampler that folds edits over them.
///
/// [`Ground`] is a pure function of position, so a chunk's base samples are
/// correct for the life of the process no matter how much is carved -- only the
/// *log* changes. Sampling it once per chunk and folding the log on the device
/// is the whole of M-155 applied here: after the first mesh of a chunk, a
/// re-mesh moves 64 bytes per surviving brush across the bus and no samples at
/// all.
///
/// Bounded at [`EXTENT`]'s 256 chunks by construction, so there is no eviction
/// policy: 256 x 17^3 x 4 bytes is 5.0 MB.
#[derive(Resource)]
struct GpuFields {
    sampler: FieldSampler,
    bases: HashMap<ChunkId, FieldBuffer>,
    /// The pruned survivors in the shader's own form, rebuilt per dispatch and
    /// **reused across every chunk of every frame** -- the same reason
    /// `World::survivors` is a field rather than a local, and it stops
    /// reallocating after the first few chunks.
    log: Vec<GpuBrush>,
    /// Bytes of *samples* uploaded since startup. On the HUD because the claim
    /// this change makes is that it stops growing.
    sample_bytes: u64,
}

/// Chunks whose geometry is on its way back from the GPU.
#[derive(Resource, Default)]
struct GpuPending {
    jobs: Vec<GpuJob>,
}

/// One chunk's geometry, in flight.
struct GpuJob {
    id: ChunkId,
    readback: Readback,
}

/// Where one GPU chunk's time went.
#[derive(Clone, Copy)]
struct GpuCost {
    /// CPU field evaluation plus the upload of this chunk's [`Ground`] base,
    /// measured around `FieldBuffer::sampled`. **Zero on every mesh after the
    /// first**, which is the claim: the base is cached and only the log moves.
    base_ms: f64,
    /// Uploading the pruned brush log and folding it over the base on the
    /// device, measured around `FieldSampler::fold_into`.
    fold_ms: f64,
    /// The three dispatches, as `isomesh-gpu` measures them.
    timings: ExtractTimings,
}

#[derive(Resource)]
struct Look {
    yaw: f32,
    pitch: f32,
}

/// Touch intent for one frame, summed into the same movement and edit paths the
/// keyboard and mouse drive.
///
/// A phone has no pointer lock -- iOS Safari does not implement it -- and it has
/// no keys, so without this the demo is watchable on a phone and not playable.
/// It is intent rather than a second control scheme: [`touch_input`] reduces
/// gestures to the same numbers `WASD` and the mouse produce, and
/// [`move_camera`] and [`dig`] read them beside the keys rather than instead of
/// them.
#[derive(Resource, Default)]
struct TouchIntent {
    /// Screen-space look delta this frame, in pixels.
    look: Vec2,
    /// Walk axes, `x` strafe, `y` forward, each in `-1.0..=1.0`.
    move_axis: Vec2,
    dig: bool,
    fill: bool,
    jump: bool,
    /// Set once any touch has been seen, and never cleared: it reveals the
    /// on-screen buttons and suppresses the pointer-lock request.
    ///
    /// Sticky rather than per-frame because both consumers are about the
    /// *device*, not the gesture: buttons that appeared and vanished between
    /// taps would be unusable, and `grab` asking a phone for pointer lock logs a
    /// console error every frame a finger is up.
    seen: bool,
}

/// Which on-screen button an entity is, so one query serves all three.
#[derive(Component)]
enum TouchButton {
    Jump,
    Dig,
    Fill,
}

/// Virtual-stick displacement, as a testable free function.
///
/// Dead zone first, so a resting thumb does not creep; then scaled by the
/// travel that remains, so the stick reaches full speed at [`STICK_RADIUS`]
/// rather than at one pixel past the dead zone; then clamped to the unit
/// **disc**, so a corner drag is not `sqrt(2)` times faster than a straight one;
/// then `y` is negated, because screen `y` grows downward and forward is up.
fn touch_axes(start: Vec2, current: Vec2) -> Vec2 {
    /// Thumb jitter, in pixels. Below this the stick is centred.
    const DEAD_ZONE: f32 = 8.0;
    /// Drag length that reads as full deflection, in pixels.
    const STICK_RADIUS: f32 = 90.0;
    let d = current - start;
    if d.length() <= DEAD_ZONE {
        return Vec2::ZERO;
    }
    let v = (d - d.normalize() * DEAD_ZONE) / (STICK_RADIUS - DEAD_ZONE);
    let v = v.clamp_length_max(1.0);
    Vec2::new(v.x, -v.y)
}

/// The terrain's material. Aliased because the full name appears in four places
/// -- the resource, the asset collection, the plugin and the spawn -- and a
/// mismatch between any two of them is a type error a hundred lines from its
/// cause.
type TerrainMaterial = ExtendedMaterial<StandardMaterial, TriplanarExtension>;

/// `load_internal_asset!` needs a handle that exists before the asset does, and
/// `MaterialExtension::fragment_shader` is a static function with no access to
/// the world -- so the handle cannot come from a resource. A fixed UUID is the
/// mechanism Bevy uses for its own built-in shaders.
const TRIPLANAR_SHADER: Handle<Shader> = uuid_handle!("6f1c9a5e-4f2b-4c7a-9d3e-2b8c5a71f0d4");

/// Layers in the stacked terrain array: grass, surface dirt, deep dirt, concrete.
///
/// The order is the array slice order and `textures/PROVENANCE.md` tables it.
/// `reinterpret_stacked_2d_as_array` slices a stacked PNG top-down, so a stack
/// built the other way up compiles, renders, and paints the walls with grass.
const TERRAIN_LAYERS: u32 = 4;
/// `settings.z` for the terrain: blend the first three layers by slope and depth.
const LAYER_BLEND: f32 = -1.0;
/// `settings.z` for the walls: array layer 3, the concrete, and nothing else.
const LAYER_CONCRETE: f32 = 3.0;

/// The triplanar half of the terrain material: two packed texture **arrays** and
/// the three numbers that place them.
#[derive(Asset, AsBindGroup, Reflect, Debug, Clone)]
struct TriplanarExtension {
    /// `x` world units per texture tile, `y` blend sharpness, `z` the forced
    /// array layer ([`LAYER_BLEND`] or [`LAYER_CONCRETE`]). See the WGSL for why
    /// this is one `vec4` and not three `f32`s.
    #[uniform(100)]
    settings: Vec4,
    /// RGB colour, A roughness, one layer per terrain material.
    ///
    /// `dimension = "2d_array"` is not cosmetic: the bind group layout it
    /// declares has to match `texture_2d_array<f32>` in the WGSL, and a
    /// mismatch is a pipeline creation failure at the first draw rather than a
    /// compile error here.
    #[texture(101, dimension = "2d_array")]
    #[sampler(102)]
    albedo_roughness: Handle<Image>,
    /// RGB OpenGL-convention normal, A ambient occlusion, same layers.
    #[texture(103, dimension = "2d_array")]
    #[sampler(104)]
    normal_ao: Handle<Image>,
}

impl MaterialExtension for TriplanarExtension {
    fn fragment_shader() -> ShaderRef {
        ShaderRef::Handle(TRIPLANAR_SHADER)
    }
    // `prepass_fragment_shader` and `deferred_fragment_shader` stay `Default` on
    // purpose: this material is opaque, so the depth and shadow passes need the
    // geometry and nothing the texture provides, and the demo runs forward --
    // there is no `DeferredPrepass` on its camera.
}

/// The terrain's material, handed to every chunk `rebuild` spawns.
///
/// The walls get a *second* `TerrainMaterial` instance rather than a second
/// material type -- the forced layer is a uniform, so the two share one
/// pipeline and one pair of texture arrays and differ only in 16 bytes of
/// `settings`. That handle is not a resource: the five `MeshMaterial3d`
/// components hold it, nothing spawns a wall after `setup`, and a resource
/// nothing reads is a field `dead_code` is right about.
#[derive(Resource)]
struct SurfaceMaterial(Handle<TerrainMaterial>);

/// One of the five cuboids that line the sandbox, so the count is assertable.
#[derive(Component)]
struct Wall;

/// A scripted sequence of carves, one per frame, driven by `ISOMESH_AUTOCARVE`.
///
/// The acceptance criterion for this ticket is about what happens *while
/// someone is clicking*, and a screenshot cannot click. Without this the example
/// could be committed compiling, rendering, and silently not carving at all —
/// so the loop runs itself, through exactly the same code path a click takes,
/// and the committed screenshot is of a tunnel that was actually dug.
#[derive(Resource, Default)]
struct AutoCarve {
    remaining: u32,
    step: u32,
    every: u32,
}

impl AutoCarve {
    fn from_env() -> Self {
        Self {
            remaining: std::env::var("ISOMESH_AUTOCARVE")
                .ok()
                .and_then(|v| v.parse().ok())
                .unwrap_or(0),
            // Captured frames per carve. One carve per frame digs the whole
            // visible tunnel in the first half-second of a clip and leaves the
            // rest of it static, which reads as a jump cut rather than as
            // digging.
            every: std::env::var("ISOMESH_AUTOCARVE_EVERY")
                .ok()
                .and_then(|v| v.parse().ok())
                .filter(|n| *n >= 1)
                .unwrap_or(1),
            step: 0,
        }
    }

    /// Where the `n`th scripted carve goes: a tunnel boring into the hill.
    fn centre(n: u32) -> Vec3 {
        let t = n as f32;
        Vec3::new(-0.9 + t * 0.30, 0.55 - t * 0.045, 2.2 - t * 0.34)
    }
}

fn main() {
    let mut app = App::new();
    app.add_plugins(DefaultPlugins.set(WindowPlugin {
        primary_window: Some(Window {
            title: "isomesh - E-202 game dig".into(),
            // Web only, inert on native: bind to the canvas the page supplies
            // rather than letting Bevy append its own. `fit_canvas_to_parent`
            // stays at its `false` default, but note what that does *not* buy:
            // winit observes the element either way and reports its CSS box as
            // the window size, so `web/style.css`'s cap on `.stage` really does
            // move the render resolution -- measured 993x558 in a 1400x900
            // window. The HUD is laid out in pixels and still fits; that is the
            // thing to re-check if the cap ever gets tighter.
            canvas: Some("#isomesh-canvas".into()),
            ..default()
        }),
        ..default()
    }))
    .add_plugins(CommonPlugin)
    // The material is a type, so its plugin is per-type: without this the
    // handle resolves and nothing is ever drawn with it.
    .add_plugins(MaterialPlugin::<TerrainMaterial>::default())
    .init_gizmo_group::<ChunkGizmos>()
    .init_gizmo_group::<GhostGizmos>()
    .init_resource::<Aim>()
    .init_resource::<GpuPending>()
    .init_resource::<TouchIntent>()
    .insert_resource(Look {
        yaw: 0.0,
        pitch: -0.15,
    })
    .insert_resource(AutoCarve::from_env())
    .add_systems(Startup, setup)
    .add_systems(
        Update,
        // Chained, because the order is load-bearing and a tuple does not
        // impose one: `touch_input` reduces the frame's touches to the same
        // numbers the keys produce and must run before anything reads them,
        // `switch_algorithm` marks the sandbox before the frame's
        // aim and edit rather than racing them, `move_camera` moves the camera
        // and resolves the body against the field, `aim` traces from where it now
        // is, `dig` marks at that point, `drain_dirty` meshes what the frame's
        // budget allows -- after `dig`, so an edit's one to eight chunks clear in
        // the same frame -- `gpu_collect` finishes whatever the GPU returned, and
        // `ghost` draws the same point. Unchained, `dig` read a camera transform
        // one frame stale -- which it did before this, invisibly.
        // `loading_modal` is last so it reads the post-drain backlog.
        (
            touch_input,
            grab,
            switch_algorithm,
            move_camera,
            aim,
            dig,
            drain_dirty,
            gpu_collect,
            ghost,
            report,
            outline_chunks,
            loading_modal,
        )
            .chain(),
    );
    // After `DefaultPlugins`, which is what inserts `Assets<Shader>`.
    load_internal_asset!(app, TRIPLANAR_SHADER, "triplanar.wgsl", Shader::from_wgsl);
    app.run();
}

/// RGB colour, A roughness: four 512-square layers stacked top-down.
///
/// One `include_bytes!` per pack, at module scope, so `setup` and the test that
/// gates the stack read the same bytes rather than two paths that can diverge.
const TERRAIN_ALBEDO_ROUGHNESS: &[u8] =
    include_bytes!("textures/terrain_albedo_roughness_array.png");
/// RGB OpenGL-convention normal, A ambient occlusion, same four layers.
const TERRAIN_NORMAL_AO: &[u8] = include_bytes!("textures/terrain_normal_ao_array.png");

/// Decode one compiled-in PNG.
///
/// Compiled in rather than loaded, for the reason the module docs give: nothing
/// copies an `assets/` tree into `web/dist`, so an `AssetServer::load` path would
/// work natively and 404 in the browser. `Image::from_buffer` decodes on the
/// calling thread, which for these two is 512x2048 of PNG twice.
fn embedded_texture(bytes: &[u8], is_srgb: bool) -> Image {
    Image::from_buffer(
        bytes,
        ImageType::Extension("png"),
        // No compressed formats: these are plain 8-bit PNGs, and WebGL2's
        // downlevel limits offer no transcoding target worth the dependency.
        CompressedImageFormats::NONE,
        is_srgb,
        // The default is `ClampToEdge`, which on a tiling triplanar texture
        // smears one row of pixels across the whole world.
        ImageSampler::Descriptor(ImageSamplerDescriptor {
            address_mode_u: ImageAddressMode::Repeat,
            address_mode_v: ImageAddressMode::Repeat,
            address_mode_w: ImageAddressMode::Repeat,
            ..ImageSamplerDescriptor::linear()
        }),
        RenderAssetUsages::RENDER_WORLD,
    )
    .expect("the packed terrain arrays are committed beside this example")
}

/// The same PNG, relabelled as a [`TERRAIN_LAYERS`]-layer 2D array.
///
/// Reinterpreted, not resized: the file is one 512x2048 image and this declares
/// its four stacked 512-squares to be four array layers, which is why
/// [`embedded_texture`] hands back an owned `Image` rather than a handle. The
/// `expect` is the only mis-pack this call can catch -- a height that does not
/// divide by four. A stack built bottom-up divides perfectly and paints the
/// walls with grass, which is why
/// `the_terrain_array_is_four_square_layers_stacked_top_down` asserts the
/// resulting extent on the committed bytes.
fn terrain_array(bytes: &[u8], is_srgb: bool) -> Image {
    let mut image = embedded_texture(bytes, is_srgb);
    image
        .reinterpret_stacked_2d_as_array(TERRAIN_LAYERS)
        .expect("the packed terrain arrays are four 512x512 layers stacked vertically");
    image
}

/// The five slabs that line the sandbox, as `(centre, size)` in world units.
///
/// A free function rather than five literals in `setup`, so the boxes can be
/// asserted without a renderer and so `EXTENT` moves the walls with the chunks.
/// Every number comes from [`sandbox`] and [`WALL_THICKNESS`].
///
/// The side walls overhang by one thickness in the other horizontal axis, so the
/// four vertical corners close rather than showing a `WALL_THICKNESS`-square gap
/// of sky, and the floor spans the overhung footprint for the same reason.
///
/// **No ceiling**, and the floor is not decoration. Fly mode exists to leave the
/// box and look down into it, so a lid would defeat it; but `Ground` is solid to
/// `-inf` while the chunks stop at `y = -5.4`, so a shaft dug to the bottom
/// reaches a depth where the mesh ends and collision does not. The floor slab is
/// what the player sees there instead of a hole into nothing.
fn walls(layout: &ChunkLayout<f32>) -> [(Vec3, Vec3); 5] {
    let (lo, hi) = sandbox(layout);
    let half = WALL_THICKNESS * 0.5;
    let mid = (lo + hi) * 0.5;
    let height = hi.y - lo.y;
    let span_x = hi.x - lo.x + 2.0 * WALL_THICKNESS;
    let span_z = hi.z - lo.z + 2.0 * WALL_THICKNESS;
    [
        (
            Vec3::new(lo.x - half, mid.y, mid.z),
            Vec3::new(WALL_THICKNESS, height, span_z),
        ),
        (
            Vec3::new(hi.x + half, mid.y, mid.z),
            Vec3::new(WALL_THICKNESS, height, span_z),
        ),
        (
            Vec3::new(mid.x, mid.y, lo.z - half),
            Vec3::new(span_x, height, WALL_THICKNESS),
        ),
        (
            Vec3::new(mid.x, mid.y, hi.z + half),
            Vec3::new(span_x, height, WALL_THICKNESS),
        ),
        (
            Vec3::new(mid.x, lo.y - half, mid.z),
            Vec3::new(span_x, WALL_THICKNESS, span_z),
        ),
    ]
}

#[allow(clippy::too_many_arguments)]
fn setup(
    mut commands: Commands,
    mut materials: ResMut<Assets<StandardMaterial>>,
    mut terrain_materials: ResMut<Assets<TerrainMaterial>>,
    mut images: ResMut<Assets<Image>>,
    mut meshes: ResMut<Assets<Mesh>>,
    mut config: ResMut<GizmoConfigStore>,
    camera: Query<Entity, With<OrbitCamera>>,
    auto: Res<AutoCarve>,
    mut flags: ResMut<ViewFlags>,
    device: Res<RenderDevice>,
    queue: Res<RenderQueue>,
) {
    // **This demo opens with the panel hidden.** It is a game first: the numbers
    // are the point of a *screenshot*, and on a screen they sit on top of the
    // rock you are trying to dig. `H` brings them back and the hint line under
    // the headline says so, so nothing is hidden without a way back.
    //
    // `ISOMESH_VIEW=hud` overrides it, which is what the committed still is
    // taken with; `nohud` still produces an empty frame for the GIFs.
    if !flags.hud_requested {
        flags.hud = false;
    }
    // The shared harness spawns an orbit camera. Take its `OrbitCamera` off
    // rather than despawning the entity: the orbit system then skips it, this
    // example drives the same camera directly, and everything else in the
    // harness that expects a camera to exist still finds one.
    for entity in &camera {
        commands
            .entity(entity)
            .remove::<OrbitCamera>()
            .insert(Transform::from_xyz(0.0, 1.70, 6.0));
    }
    let (chunk_gizmos, _) = config.config_mut::<ChunkGizmos>();
    chunk_gizmos.line.width = 2.0;
    let (ghost_gizmos, _) = config.config_mut::<GhostGizmos>();
    ghost_gizmos.line.width = 1.5;
    ghost_gizmos.depth_bias = -1.0;

    let layout = ChunkLayout::<f32>::new(
        CHUNK_CELLS,
        CELL_SIZE,
        [
            -(EXTENT[0] as f32) * CHUNK_CELLS as f32 * CELL_SIZE * 0.5,
            // The ceiling stays at +2.6, so this is 5.9 units of rock to dig
            // down through instead of 0.9. `AutoCarve::centre` descends to
            // y ~= -2.1 over its 60 steps, which was below the old floor and is
            // inside this box -- so those constants are left alone and the
            // committed capture sequence is unchanged.
            -5.4,
            -(EXTENT[2] as f32) * CHUNK_CELLS as f32 * CELL_SIZE * 0.5,
        ],
    )
    .expect("valid layout");

    let albedo_roughness = images.add(terrain_array(TERRAIN_ALBEDO_ROUGHNESS, true));
    let normal_ao = images.add(terrain_array(TERRAIN_NORMAL_AO, false));
    // One closure, two instances: the terrain blends layers 0-2 and the walls
    // force layer 3. Everything else about them -- the pipeline, the two
    // textures, the tiling and the sharpness -- is the same, so a difference
    // anywhere else here would be a difference nothing asked for.
    let mut terrain_skin = |layer: f32| {
        terrain_materials.add(TerrainMaterial {
            base: StandardMaterial {
                // White, not the old 0.62/0.58/0.52 tint: the shader multiplies
                // the sampled colour by this, so a tint here would darken the
                // texture. The field stays as the knob it is, set to neutral.
                base_color: Color::WHITE,
                perceptual_roughness: 0.85,
                ..default()
            },
            extension: TriplanarExtension {
                // 1.5 world units per tile against a 2.0-unit chunk, so a tile
                // is a little smaller than a chunk and the cracks read at the
                // scale a brush of radius 0.25 cuts at. `y` is the blend
                // sharpness, `z` the forced layer.
                settings: Vec4::new(1.5, 4.0, layer, 0.0),
                albedo_roughness: albedo_roughness.clone(),
                normal_ao: normal_ao.clone(),
            },
        })
    };
    let material = terrain_skin(LAYER_BLEND);
    let wall_material = terrain_skin(LAYER_CONCRETE);

    // The sandbox, made visible. `sandbox` is the box the 256 chunks cover and
    // until now only `aim` and `trace` consumed it -- so the boundary existed,
    // stopped the ghost, and nothing on screen said where it was. [`walls`]
    // derives the five boxes; this only spawns them.
    for (centre, size) in walls(&layout) {
        commands.spawn((
            Mesh3d(meshes.add(Cuboid::new(size.x, size.y, size.z))),
            MeshMaterial3d(wall_material.clone()),
            Transform::from_translation(centre),
            Wall,
        ));
    }

    // `bevy::math::primitives::Sphere` in full: `Sphere` in this file is
    // `isomesh::fields::Sphere`, the SDF. `ico(3)` is 1,280 triangles; the
    // default `Ico { subdivisions: 5 }` is 20,480, which is a lot of
    // alpha-blended fill for a cursor.
    let ghost_mesh = meshes.add(
        bevy::math::primitives::Sphere::new(1.0)
            .mesh()
            .ico(3)
            .expect("3 subdivisions is far below the 80-subdivision cap"),
    );
    let skin = |rgba: [f32; 4]| StandardMaterial {
        base_color: Color::srgba(rgba[0], rgba[1], rgba[2], rgba[3]),
        alpha_mode: AlphaMode::Blend,
        // Both faces, so the far side of the ghost shows through the near one
        // and it reads as a volume rather than as a disc.
        cull_mode: None,
        // Unlit, so the cursor does not change colour as you turn: it is a
        // readout, not a thing in the world.
        unlit: true,
        ..default()
    };
    let skins = GhostMaterials {
        carve: materials.add(skin([0.98, 0.42, 0.28, 0.28])),
        fill: materials.add(skin([0.36, 0.86, 1.00, 0.28])),
    };
    commands.spawn((
        Mesh3d(ghost_mesh),
        MeshMaterial3d(skins.carve.clone()),
        Transform::default(),
        Ghost,
    ));
    commands.insert_resource(skins);

    let mut world = World {
        layout,
        brushes: Vec::new(),
        survivors: Vec::new(),
        dirty: DirtySet::new(),
        backlog: 0,
        // The startup fill is a drain like any other, so it shows the same modal
        // rather than blocking the first frame for most of a second.
        switching: true,
        radius: 0.25,
        last_touched: Vec::new(),
        last_edit: None,
        last_edit_ms: 0.0,
        last_chunks: 0,
        show_chunks: true,
        grabbed: false,
        algorithm: Algorithm::from_env(),
        switch_ms: 0.0,
        switch_chunks: 0,
        switch_started: Instant::now(),
        stroke_last: None,
        stroke_edits: 0,
        stroke_ms: 0.0,
        // Starts due, so the first click of a hold lands on the frame it is made
        // rather than a twentieth of a second later.
        stroke_clock: EDIT_PERIOD,
        velocity: Vec3::ZERO,
        // Walk by default, fly when a scripted capture is running: the committed
        // clip is a flight along a hardcoded tunnel, and a body that falls to the
        // floor would film a different thing.
        walking: auto.remaining == 0,
        grounded: false,
        gpu_cost: None,
        gpu_peak_triangles: 0,
        gpu_failures: 0,
    };

    // Mark every chunk once. `drain_dirty` meshes them under its frame budget,
    // and after this only edited chunks are ever marked.
    world.switch_chunks = mark_sandbox(&mut world.dirty);

    // Above the tree's existing `GlobalZIndex(4)` banner layer, because this is
    // meant to cover everything including the HUD. Spawned hidden;
    // `loading_modal` is the only thing that changes that.
    commands
        .spawn((
            Node {
                position_type: PositionType::Absolute,
                left: Val::Px(0.0),
                right: Val::Px(0.0),
                top: Val::Px(0.0),
                bottom: Val::Px(0.0),
                justify_content: JustifyContent::Center,
                align_items: AlignItems::Center,
                ..default()
            },
            BackgroundColor(Color::srgba(0.02, 0.03, 0.05, 0.72)),
            GlobalZIndex(10),
            Visibility::Hidden,
            LoadingModal,
        ))
        .with_child((
            Text::new("Loading… please wait."),
            TextFont {
                font_size: FontSize::Px(22.0),
                ..default()
            },
            TextColor(Color::srgb(0.92, 0.94, 0.98)),
        ));

    // The touch controls, spawned hidden. `touch_input` flips them to
    // `Display::Flex` the first time a finger lands, so a desktop run looks
    // exactly as it did and a phone is not asked to guess that the left third of
    // the screen walks.
    //
    // `Interaction` is the hit test: `bevy_ui`'s focus system already reads
    // `Touches`, so a tap on one of these is reported without this file doing
    // any rectangle arithmetic. These are the first `Button`s anywhere under
    // `examples/`; the HUD palette is reused so they read as part of it.
    //
    // `GlobalZIndex(5)` is above the tree's `4` banner layer and below the
    // loading modal's `10`, which is right: a control you cannot see the effect
    // of is worse than a hidden one.
    //
    // `JUMP` sits 120 px up on the left, clear of the stick's own resting zone;
    // `DIG` and `FILL` are bottom-right under the thumb that is already looking.
    for (button, label, left, right, bottom) in [
        (TouchButton::Jump, "JUMP", Some(24.0), None, 120.0),
        (TouchButton::Dig, "DIG", None, Some(112.0), 24.0),
        (TouchButton::Fill, "FILL", None, Some(24.0), 24.0),
    ] {
        commands
            .spawn((
                Button,
                Node {
                    display: Display::None,
                    position_type: PositionType::Absolute,
                    left: left.map_or(Val::Auto, Val::Px),
                    right: right.map_or(Val::Auto, Val::Px),
                    bottom: Val::Px(bottom),
                    width: Val::Px(72.0),
                    height: Val::Px(72.0),
                    justify_content: JustifyContent::Center,
                    align_items: AlignItems::Center,
                    ..default()
                },
                BackgroundColor(Color::srgba(0.02, 0.03, 0.05, 0.72)),
                GlobalZIndex(5),
                button,
            ))
            .with_child((
                Text::new(label),
                TextFont {
                    font_size: FontSize::Px(16.0),
                    ..default()
                },
                TextColor(Color::srgb(0.92, 0.94, 0.98)),
            ));
    }

    // Bevy's own device, handed to a crate that has never heard of Bevy --
    // `RenderDevice` and `RenderQueue` are inserted into the *main* world by
    // `bevy_render`, so a plain `Startup` system can reach them and there is no
    // `RenderApp` in this file. `gpu_compute_mc.rs` is the example that exists to
    // prove that.
    //
    // `expect`, not a fallback to the CPU: a browser without WebGPU cannot run
    // this demo at all, which `web/play.html` now says before it imports the
    // module, so there is nothing to degrade to and a silent CPU path would make
    // key `8` a lie.
    commands.insert_resource(GpuMesher {
        mc: MarchingCubesGpu::new(device.wgpu_device(), &queue)
            .expect("the compute pipeline for key 8; WebGPU or a native backend is required"),
    });
    // The same device and the same `expect` shape: without this pipeline key `8`
    // could still extract, but only by uploading 19.6 KB of samples per chunk
    // per edit -- the configuration `isomesh-gpu`'s own docs call the wrong one.
    commands.insert_resource(GpuFields {
        sampler: FieldSampler::new(device.wgpu_device())
            .expect("the field-fold pipeline for key 8; WebGPU or a native backend is required"),
        bases: HashMap::new(),
        log: Vec::new(),
        sample_bytes: 0,
    });

    commands.insert_resource(SurfaceMaterial(material));
    commands.insert_resource(world);
}

// ── the bound ───────────────────────────────────────────────────────────────
//
// Copied out of `bevy_isomesh/examples/tape_pruning.rs` -- its `f32` variant,
// specifically, the one that returns an infinite interval for a field declaring
// no constant rather than panicking, which is the right behaviour in a demo a
// stranger runs. Transcribed rather than imported, the same house convention
// `game_edit_tape_trim.rs` states: examples in this directory do not `use` one
// another, and `Brush`/`BrushStack` have no `BoundedSdf` impl to import a
// pruner against.
//
// What it buys, measured in `tape_pruning.rs`: median surviving fraction 0.2969,
// median per-chunk speedup 3.365x, world aggregate 2.473x, and the mesh
// byte-identical on 64 of 64 chunks. That last clause is the whole licence for
// doing this in the example that reports E1: a dropped brush provably cannot
// change the fold anywhere in the box, so the numbers on screen are unchanged.

/// The enclosure of a scalar field over a chunk.
#[derive(Clone, Copy, Debug)]
struct Interval {
    /// Lower bound.
    lo: f32,
    /// Upper bound.
    hi: f32,
}

/// The box a chunk's field evaluations can touch.
#[derive(Clone, Copy, Debug)]
struct ChunkBox {
    /// Centre of the sampled extent.
    centre: [f32; 3],
    /// Circumradius: half the diagonal of the sampled extent, plus the margin
    /// [`Sdf::gradient`]'s central differences reach outside it.
    radius: f32,
}

impl ChunkBox {
    /// The box for a chunk whose sample grid starts at `origin` and spans `span`
    /// on every axis.
    fn new(origin: [f32; 3], span: f32) -> Self {
        let centre = [
            origin[0] + span * 0.5,
            origin[1] + span * 0.5,
            origin[2] + span * 0.5,
        ];
        // `h = DIFF_STEP * max(|p|, 1)` at the furthest corner bounds the
        // differencing reach anywhere in the box.
        let mut far = 1.0f32;
        for lo in origin {
            far = far.max(lo.abs()).max((lo + span).abs());
        }
        let half = span * 0.5 + <f32 as Real>::DIFF_STEP * far;
        Self {
            centre,
            radius: half * 3.0f32.sqrt(),
        }
    }
}

/// Slack for one bound, in absolute units.
fn pad(value: f32, reach: f32) -> f32 {
    PAD_ULPS * f32::EPSILON * (value.abs() + reach)
}

/// `f(centre) ± l·r`, widened so it is an enclosure and not an estimate.
///
/// A field that declares no Lipschitz constant gets an infinite enclosure, which
/// makes it unprunable rather than wrongly prunable. Nothing here takes that
/// path — `Ground` declares [`LIPSCHITZ`] and every brush is an exact sphere —
/// but the alternative would be an `expect` in a demo a stranger runs.
fn enclose<S: BoundedSdf<Scalar = f32>>(field: &S, chunk: ChunkBox) -> Interval {
    let Some(l) = field.value_bound().lipschitz() else {
        return Interval {
            lo: f32::NEG_INFINITY,
            hi: f32::INFINITY,
        };
    };
    let value = field.sample(chunk.centre);
    let reach = l as f32 * chunk.radius;
    let slack = reach + pad(value, reach);
    Interval {
        lo: value - slack,
        hi: value + slack,
    }
}

/// Select the brushes that can still change the fold anywhere in `chunk`.
///
/// Order is preserved, because `Add` and `Subtract` do not commute with each
/// other — [`BrushOp::commutes_with`] is the crate's own statement of that.
/// Returns the number of survivors, which is `out.len()`.
///
/// `out` is reused across chunks, which is why this takes `&mut Vec` and calls
/// `clear()`: after the first few chunks it never reallocates.
fn prune_into(
    tape: &[Brush<Sphere<f32>>],
    base: &Ground,
    chunk: ChunkBox,
    out: &mut Vec<Brush<Sphere<f32>>>,
) -> usize {
    out.clear();
    let mut v = enclose(base, chunk);
    for brush in tape {
        let s = enclose(&brush.shape, chunk);
        match brush.op {
            BrushOp::Add => {
                // Strictly above the running value everywhere, so `min` is
                // forced to select the running value at every point.
                if s.lo > v.hi {
                    continue;
                }
                v = Interval {
                    lo: v.lo.min(s.lo),
                    hi: v.hi.min(s.hi),
                };
            }
            BrushOp::Subtract => {
                // `max(v, -s)`. Negation is exact, so negating the enclosure is
                // exact too.
                let n = Interval {
                    lo: -s.hi,
                    hi: -s.lo,
                };
                if n.hi < v.lo {
                    continue;
                }
                v = Interval {
                    lo: v.lo.max(n.lo),
                    hi: v.hi.max(n.hi),
                };
            }
            BrushOp::SmoothAdd { k } => {
                // **Never pruned in the losing direction**, and that is the
                // registered asymmetry rather than an omission: `smooth_min` at
                // `h == 1` returns `b + (a - b)`, which is not bit-identical to
                // `a`. The enclosure still has to track it, and its floor sags
                // by at most `k/4` below the plain minimum.
                //
                // This example pushes only `Brush::add` and `Brush::subtract`,
                // so the arm is unreachable here -- it stays for `BrushOp`'s
                // exhaustiveness, and because a wrong arm added later would be
                // silent.
                let k = k as f32;
                let lo = v.lo.min(s.lo) - 0.25 * k;
                v = Interval {
                    lo: lo - pad(lo, 0.25 * k),
                    hi: v.hi.min(s.hi),
                };
            }
        }
        out.push(*brush);
    }
    out.len()
}

/// The editable region, in world units: the box the 256 chunks cover.
///
/// The field is unbounded — `Ground` has a height at every `x` and `z` — but only
/// these chunks are ever meshed. A brush outside them edits rock nobody can see
/// and spawns an island chunk beside the sandbox, so the aim stops at the wall
/// and the ghost vanishing is how the demo says so.
fn sandbox(layout: &ChunkLayout<f32>) -> (Vec3, Vec3) {
    let span = layout.cell_size() * CHUNK_CELLS as f32;
    let lo = Vec3::from_array(layout.sample_origin(ChunkId::new([0, 0, 0])));
    let size = Vec3::new(EXTENT[0] as f32, EXTENT[1] as f32, EXTENT[2] as f32) * span;
    (lo, lo + size)
}

/// First surface crossing along a ray inside the sandbox, as a distance, or
/// `None` if there is none within [`AIM_FAR`].
///
/// Sphere tracing rather than a triangle raycast, and that is the whole point:
/// the field is the thing being edited, so tracing the field cannot go stale.
/// The only ray code in this directory is `game_walk`'s parry3d `TriMesh` cast,
/// which would need a per-chunk collider cache that every edit invalidates.
///
/// `f <= AIM_HIT` on the first sample means the camera is already inside rock,
/// and returning `AIM_NEAR` there is deliberate: it is how you dig yourself out.
fn trace(
    field: &impl Sdf<Scalar = f32>,
    origin: Vec3,
    direction: Vec3,
    bounds: (Vec3, Vec3),
) -> Option<f32> {
    let mut t = AIM_NEAR;
    for _ in 0..AIM_STEPS {
        let p = origin + direction * t;
        let f = field.sample([p.x, p.y, p.z]);
        // The box test rides along with the surface test rather than gating the
        // march, because the camera may be *outside* the sandbox looking in: a
        // ray that starts out keeps stepping and can still hit once it is inside.
        // A box is convex, so one that has left cannot come back, and `t` runs
        // out at `AIM_FAR`.
        if f <= AIM_HIT && p.cmpge(bounds.0).all() && p.cmple(bounds.1).all() {
            return Some(t);
        }
        t += (f / LIPSCHITZ).max(AIM_HIT);
        if t > AIM_FAR {
            return None;
        }
    }
    None
}

fn aim(camera: Query<&Transform, With<Camera3d>>, world: Res<World>, mut target: ResMut<Aim>) {
    let Ok(view) = camera.single() else {
        return;
    };
    let field = BrushStack {
        base: Ground,
        brushes: &world.brushes,
    };
    let origin = view.translation;
    let direction = *view.forward();
    match trace(&field, origin, direction, sandbox(&world.layout)) {
        Some(t) => {
            target.point = origin + direction * t;
            target.hit = true;
        }
        None => {
            target.point = origin + direction * AIM_FAR;
            target.hit = false;
        }
    }
}

/// Draw the ghost: a translucent sphere of exactly the brush that a click would
/// push, at exactly where it would go.
fn ghost(
    world: Res<World>,
    target: Res<Aim>,
    buttons: Res<ButtonInput<MouseButton>>,
    skins: Res<GhostMaterials>,
    mut gizmos: Gizmos<GhostGizmos>,
    mut query: Query<
        (
            &mut Transform,
            &mut MeshMaterial3d<StandardMaterial>,
            &mut Visibility,
        ),
        With<Ghost>,
    >,
) {
    let filling = buttons.pressed(MouseButton::Right);
    for (mut transform, mut material, mut visibility) in &mut query {
        if !target.hit {
            *visibility = Visibility::Hidden;
            continue;
        }
        *visibility = Visibility::Inherited;
        *transform =
            Transform::from_translation(target.point).with_scale(Vec3::splat(world.radius));
        let skin = if filling { &skins.fill } else { &skins.carve };
        if material.0 != *skin {
            material.0 = skin.clone();
        }
    }
    if target.hit {
        // The ghost is centred on the surface, so half of it is inside the rock
        // and depth-tested away. The outline is drawn with `depth_bias = -1.0`
        // -- in front of everything -- so the buried half still reads.
        let colour = if filling {
            Color::srgb(0.45, 0.90, 1.00)
        } else {
            Color::srgb(1.00, 0.55, 0.35)
        };
        gizmos.sphere(
            Isometry3d::from_translation(target.point),
            world.radius,
            colour,
        );
    }
}

/// Which mesher builds a chunk: seven CPU extractors and Marching Cubes on the
/// GPU.
///
/// The eighth is the point of having eight. The rendering was always on the GPU
/// -- this is a triplanar `ExtendedMaterial` -- and what key `8` moves there is
/// the **extraction and the field**: a compute shader folds the pruned brush log
/// over a base sampled once per chunk, classifies every cell, prefix-scans, and
/// writes the triangles. After a chunk's first mesh the only CPU work is 64 bytes
/// per surviving brush -- see [`gpu_dispatch`]. Same algorithm as key `1`, so the
/// triangle counts are comparable and the HUD prints where the time went.
///
/// `subgrid_marching_tetrahedra` is **not** offered: it is ~196x Marching Cubes,
/// and this demo re-meshes on every frame of a held mouse button. A key that
/// makes the page stop responding for seconds is not a choice a reader can use.
/// `isomesh_web`'s `NOT_OFFERED` excludes it for the same reason.
#[derive(Clone, Copy, PartialEq, Eq)]
enum Algorithm {
    MarchingCubes,
    MarchingCubesDecider,
    MarchingTetrahedra,
    SurfaceNets,
    DualContouring,
    ManifoldDualContouring,
    GreedyQuads,
    MarchingCubesGpu,
}

impl Algorithm {
    /// In key order: `1` is the reference and the default.
    const ALL: [Self; 8] = [
        Self::MarchingCubes,
        Self::MarchingCubesDecider,
        Self::MarchingTetrahedra,
        Self::SurfaceNets,
        Self::DualContouring,
        Self::ManifoldDualContouring,
        Self::GreedyQuads,
        Self::MarchingCubesGpu,
    ];

    /// The name on the HUD.
    ///
    /// The first six are `isomesh::extractor::ALL_EXTRACTORS`'s own spellings, so
    /// a name here and a name in the registry are the same string. `greedy_quads`
    /// is deliberately **not** one of them: it is
    /// `isomesh::extractor::UNREGISTERED`, because it classifies whole cells
    /// solid or empty and emits the axis-aligned faces between them -- a
    /// Minecraft surface rather than an isosurface, which is why the crate keeps
    /// it out of every sweep. It gets a key here because a digging sandbox is
    /// exactly the case the crate's own comment allows for: *"a caller may
    /// legitimately want to drive it through the same shape."*
    fn name(self) -> &'static str {
        match self {
            Self::MarchingCubes => "marching_cubes",
            Self::MarchingCubesDecider => "marching_cubes+decider",
            Self::MarchingTetrahedra => "marching_tetrahedra",
            Self::SurfaceNets => "surface_nets",
            Self::DualContouring => "dual_contouring",
            Self::ManifoldDualContouring => "manifold_dual_contouring",
            Self::GreedyQuads => "greedy_quads",
            Self::MarchingCubesGpu => "marching_cubes (gpu)",
        }
    }

    /// The HUD headline's colour for this mesher.
    ///
    /// Eight hues, evenly spaced around the wheel and all light enough to read
    /// on dark rock, so which mesher is running is answerable from the corner of
    /// an eye rather than by reading a word. Assigned in **key order** -- `1` is
    /// amber and `8` is pink -- so the mapping is the same thing the key row
    /// already teaches.
    ///
    /// Saturated rather than pastel on purpose: the panel below is a near-white
    /// grey, and a headline that shares its colour is not a headline.
    fn colour(self) -> Color {
        match self {
            Self::MarchingCubes => Color::srgb(1.00, 0.62, 0.16),
            Self::MarchingCubesDecider => Color::srgb(0.96, 0.90, 0.25),
            Self::MarchingTetrahedra => Color::srgb(0.60, 0.95, 0.30),
            Self::SurfaceNets => Color::srgb(0.25, 0.92, 0.55),
            Self::DualContouring => Color::srgb(0.25, 0.88, 0.95),
            Self::ManifoldDualContouring => Color::srgb(0.45, 0.62, 1.00),
            Self::GreedyQuads => Color::srgb(0.80, 0.50, 1.00),
            Self::MarchingCubesGpu => Color::srgb(1.00, 0.40, 0.70),
        }
    }

    /// The digit that selects this mesher, `1`-`8`.
    ///
    /// Read out of [`Self::ALL`] rather than written a second time, because
    /// `switch_algorithm`'s `DIGITS` already indexes that array: two hand-kept
    /// lists would let the HUD advertise a key that selects something else.
    fn key(self) -> usize {
        Self::ALL
            .iter()
            .position(|a| *a == self)
            .expect("every variant is in ALL")
            + 1
    }

    /// `ISOMESH_ALGORITHM=<name>`, so a capture needs no keyboard -- the same
    /// reason `ISOMESH_FIELD` and `ISOMESH_VIEW` exist. Anything else, including
    /// unset, is Marching Cubes.
    fn from_env() -> Self {
        let want = std::env::var("ISOMESH_ALGORITHM").unwrap_or_default();
        Self::ALL
            .into_iter()
            .find(|a| a.name() == want)
            .unwrap_or(Self::MarchingCubes)
    }

    /// Extra cells to borrow from the positive neighbours, beyond the chunk's
    /// own.
    ///
    /// A dual method emits one quad per crossed grid edge and that quad needs
    /// all four cells around the edge, so `dual.rs`'s walk skips the outermost
    /// quad plane on every face -- a chunk given only its own cells stops one
    /// cell short and neither neighbour emits the bridging quad. One borrowed
    /// cell layer makes the planes tile exactly once: this chunk emits global
    /// planes `base+1 ..= base+cells`, the positive neighbour emits
    /// `base+cells+1 ..= base+2*cells`, and plane `base` comes from the negative
    /// neighbour's last. No gap and no duplicate.
    ///
    /// Two preconditions, both true here and both load-bearing. `CELL_SIZE` is a
    /// power of two, which is what makes the borrowed cell's samples -- and
    /// therefore its vertex -- bit-identical in the two chunks that compute it
    /// (M-32). And `SurfaceNets::set_smoothing_passes` stays at its default 0:
    /// `DualMesher::smooth` averages a vertex with its face-adjacent cells, and
    /// those differ between two chunks. This example never calls the setter.
    ///
    /// Zero for the edge-based families. Marching Cubes and Marching Tetrahedra
    /// march *every* cell of the grid and put vertices on grid edges both chunks
    /// compute from identical corner samples, so they already tile; a borrowed
    /// layer would only make them mesh the neighbour's first cell twice.
    ///
    /// `greedy_quads` is zero for a third reason: it caps the domain boundary on
    /// purpose, so each chunk is a closed box and a borrowed layer would only
    /// move the wall outward.
    fn halo_cells(self) -> u32 {
        match self {
            Self::SurfaceNets | Self::DualContouring | Self::ManifoldDualContouring => 1,
            _ => 0,
        }
    }

    /// Extract one chunk on the CPU, or `None` when it holds no surface -- so
    /// empty air costs a sample loop rather than an entity and a draw call over
    /// nothing.
    ///
    /// [`Algorithm::MarchingCubesGpu`] returns `None` here and is *not* a chunk
    /// with no surface: its geometry is asynchronous and arrives through
    /// [`gpu_dispatch`] and [`gpu_collect`]. [`drain_dirty`] routes on the
    /// variant and never calls this for it.
    fn mesh<F: Sdf<Scalar = f32>>(
        self,
        layout: &ChunkLayout<f32>,
        field: &F,
        origin: [f32; 3],
    ) -> Option<Mesh> {
        // `+ 1` is the positive-face sample overlap `sample_shape()` provides;
        // the halo is a further *cell*. The origin is unchanged -- it must come
        // from `sample_origin(id)` and never from arithmetic on it, or the two
        // chunks disagree about where their shared plane is (M-32).
        let shape = RuntimeShape3::new([layout.cells() + 1 + self.halo_cells(); 3]).ok()?;
        let cell = layout.cell_size();
        let mut builder = MeshBuilder::new();
        let extracted = match self {
            Self::MarchingCubes => {
                MarchingCubes::<f32>::new().extract(field, &shape, origin, cell, &mut builder)
            }
            Self::MarchingCubesDecider => {
                let mut mesher = MarchingCubes::<f32>::new();
                mesher.set_face_ambiguity(FaceAmbiguity::AsymptoticDecider);
                mesher.extract(field, &shape, origin, cell, &mut builder)
            }
            Self::MarchingTetrahedra => {
                MarchingTetrahedra::<f32>::new().extract(field, &shape, origin, cell, &mut builder)
            }
            Self::SurfaceNets => {
                SurfaceNets::<f32>::new().extract(field, &shape, origin, cell, &mut builder)
            }
            Self::DualContouring => {
                DualContouring::<f32>::new().extract(field, &shape, origin, cell, &mut builder)
            }
            Self::ManifoldDualContouring => ManifoldDualContouring::<f32>::new().extract(
                field,
                &shape,
                origin,
                cell,
                &mut builder,
            ),
            Self::GreedyQuads => {
                GreedyQuads::<f32>::new().extract(field, &shape, origin, cell, &mut builder)
            }
            // Not reachable: `drain_dirty` sends this variant to `gpu_dispatch`.
            // Spelled out rather than caught by a `_` arm, so adding a ninth
            // mesher is a compile error here rather than a silently empty chunk.
            Self::MarchingCubesGpu => return None,
        };
        extracted.ok()?;
        if builder.indices().is_empty() {
            return None;
        }
        Some(builder.into_mesh())
    }
}

/// Re-mesh dirty chunks nearest the camera first, stopping when the frame budget
/// is gone and keeping the rest for next frame.
///
/// **This is the only place a chunk is meshed.** `setup`, `dig` and
/// `switch_algorithm` mark and return. A 256-chunk switch that used to block for
/// most of a second now drains over a few dozen frames with the frame time flat,
/// which is what makes [`loading_modal`]'s overlay a truthful readout rather
/// than a splash screen -- and it is what stops a held mouse button from
/// collapsing to single-figure FPS, because a frame can no longer be charged for
/// an unbounded queue.
///
/// The four reconcile arms are the whole state machine: a chunk that had a mesh
/// and still has one is replaced, one that lost its surface is despawned rather
/// than left as an empty draw call, one that gained a surface is spawned, and one
/// that never had either costs nothing.
#[allow(clippy::too_many_arguments)]
fn drain_dirty(
    mut world: ResMut<World>,
    mut pending: ResMut<GpuPending>,
    mut fields: ResMut<GpuFields>,
    gpu: Res<GpuMesher>,
    device: Res<RenderDevice>,
    queue: Res<RenderQueue>,
    mut commands: Commands,
    mut meshes: ResMut<Assets<Mesh>>,
    material: Res<SurfaceMaterial>,
    chunks: Query<(Entity, &Chunk)>,
    camera: Query<&Transform, With<Camera3d>>,
) {
    // One place where a switch is declared finished, and it is here rather than
    // after the drain because the last GPU chunk lands a frame or more after the
    // queue empties: the modal must stay up until the geometry is on screen. The
    // cost is a frame of latency in `switch_ms`, which is a wall-clock figure
    // across dozens of frames.
    if world.switching && world.dirty.is_empty() && pending.jobs.is_empty() {
        world.switch_ms = world.switch_started.elapsed().as_secs_f64() * 1000.0;
        world.switching = false;
    }
    if world.dirty.is_empty() {
        world.backlog = 0;
        return;
    }
    let eye = camera.single().map_or(Vec3::ZERO, |view| view.translation);
    // Copied out before the closure borrows `world.brushes`. `ChunkLayout<f32>`
    // is `Copy`, which is what makes that free.
    let layout = world.layout;
    let algorithm = world.algorithm;
    let started = Instant::now();
    let mut built: Vec<(ChunkId, Option<Mesh>)> = Vec::new();
    let mut gpu_cost = None;
    // A `Cell` because two closures need it: the mesher pushes jobs and the
    // budget predicate reads the count. They cannot both hold `&mut`, and the
    // count is the backpressure -- see `GPU_JOBS_MAX`.
    let in_flight = Cell::new(pending.jobs.len());
    // One mutable borrow for the closure: the base cache is read on every GPU
    // chunk and written on the first mesh of each.
    let fields = &mut *fields;
    // Disjoint field borrows: the closure reads `brushes` and writes `survivors`
    // while `dirty` is borrowed mutably. Destructuring is what makes that legal.
    let World {
        dirty,
        brushes,
        survivors,
        ..
    } = &mut *world;
    let report = dirty.mesh_within_budget(
        &layout,
        [eye.x, eye.y, eye.z],
        |id, origin| {
            // The `(L + 1)` factor on every one of this chunk's 4,913 samples,
            // cut to the brushes whose spheres can still reach it. Bit-exact:
            // a dropped brush provably cannot change the fold anywhere in the
            // box, so this is speed with no change to the geometry.
            //
            // Both paths prune against the same box and are handed the same
            // survivors; what each does with them is where the keys part
            // company. The CPU walks them once per sample, the GPU uploads 64
            // bytes per brush and folds them on the device -- see
            // [`gpu_dispatch`].
            let span = layout.cell_size() * layout.cells() as f32;
            let box_ = ChunkBox::new(origin, span);
            prune_into(brushes, &Ground, box_, survivors);
            if algorithm == Algorithm::MarchingCubesGpu {
                if let Some((job, cost)) = gpu_dispatch(
                    &gpu.mc,
                    fields,
                    device.wgpu_device(),
                    &queue,
                    &layout,
                    survivors,
                    id,
                    origin,
                ) {
                    pending.jobs.push(job);
                    in_flight.set(in_flight.get() + 1);
                    gpu_cost = Some(cost);
                }
                return;
            }
            let field = BrushStack {
                base: Ground,
                brushes: survivors,
            };
            built.push((id, algorithm.mesh(&layout, &field, origin)));
        },
        // The predicate the crate asks for: a `no_std` crate cannot read a
        // clock, so the caller owns it. The second clause is backpressure, not a
        // second budget -- a stopped drain leaves its chunks in the set, so the
        // work is deferred rather than dropped.
        || started.elapsed() < MESH_BUDGET && in_flight.get() < GPU_JOBS_MAX,
    );
    // Asset creation is deliberately outside the timed region: the budget bounds
    // *extraction*, and charging it for `Assets::add` would measure Bevy.
    for (id, mesh) in built {
        reconcile(
            &mut commands,
            &chunks,
            &material,
            id,
            mesh.map(|m| meshes.add(m)),
        );
    }
    if let Some(cost) = gpu_cost {
        world.gpu_cost = Some(cost);
    }
    // The next frame's head is where a drained switch is noticed, so a GPU chunk
    // still in flight keeps the modal up.
    world.backlog = report.remaining;
}

/// Attach, replace or drop the mesh of one chunk.
///
/// The four arms are the whole state machine: a chunk that had a mesh and still
/// has one is replaced, one that lost its surface is despawned rather than left
/// as an empty draw call, one that gained a surface is spawned, and one that
/// never had either costs nothing.
///
/// One function, because the CPU drain and the GPU collection both end here and
/// two copies of a four-arm reconcile is two places for the despawn to go
/// missing.
fn reconcile(
    commands: &mut Commands,
    chunks: &Query<(Entity, &Chunk)>,
    material: &SurfaceMaterial,
    id: ChunkId,
    handle: Option<Handle<Mesh>>,
) {
    let existing = chunks.iter().find(|(_, c)| c.0 == id).map(|(e, _)| e);
    match (existing, handle) {
        (Some(entity), Some(handle)) => {
            commands.entity(entity).insert(Mesh3d(handle));
        }
        (Some(entity), None) => {
            commands.entity(entity).despawn();
        }
        (None, Some(handle)) => {
            commands.spawn((
                Mesh3d(handle),
                MeshMaterial3d(material.0.clone()),
                Chunk(id),
            ));
        }
        (None, None) => {}
    }
}

/// Fold the pruned brush log over this chunk's cached base **on the device**,
/// extract on the GPU, and start the readback.
///
/// # The field is evaluated on the device, and that is the point of key `8`
///
/// `isomesh-gpu`'s headline finding is that where the field is evaluated decides
/// everything else: sample on the CPU and hand over a `FieldBuffer` and the
/// **upload is 87% of the path**, which does not take field evaluation off the
/// CPU's budget -- it adds a copy to it. So the eight keys deliberately do
/// **not** pay the same field cost any more. Seven of them walk the survivors
/// once per sample on the CPU; this one uploads 64 bytes per surviving brush and
/// folds them in a compute shader.
///
/// `FieldBuffer::sampled` still appears here exactly once per chunk, because
/// [`Ground`] is not one of `field.wgsl`'s four base fields: the base is sampled
/// on the CPU and uploaded the first time a chunk is meshed, and every later
/// re-mesh of that chunk is `FieldSampler::fold_into` over the cached buffer.
/// [`Ground`] is a pure function of position, so that buffer can never be stale.
/// The HUD prints `base` and `fold` separately and marks the base `(cached)`,
/// which is how a reader sees the per-edit sample upload reach zero.
///
/// The output is a fresh `FieldBuffer` per dispatch rather than one reused
/// scratch buffer: `fold_into` reads sample *positions* from `out.params()`, so
/// the output must carry this chunk's own origin. That is a 19.6 KB allocation
/// and correctness beats saving it.
///
/// `extract_indirect`, not `extract` or `extract_buffers`: it is the **only**
/// entry point that reads nothing back. The others end in
/// `read_buffer_u32` -> `read_bytes_many` -> `device.poll(PollType::Wait)` +
/// `mpsc::recv()`, and on WebGPU `Device::poll` is a documented no-op while
/// `map_async` completes only from the browser's event loop -- so that wait
/// wedges the tab rather than blocking a thread.
#[allow(clippy::too_many_arguments)]
fn gpu_dispatch(
    mc: &MarchingCubesGpu,
    fields: &mut GpuFields,
    device: &wgpu::Device,
    queue: &wgpu::Queue,
    layout: &ChunkLayout<f32>,
    survivors: &[Brush<Sphere<f32>>],
    id: ChunkId,
    origin: [f32; 3],
) -> Option<(GpuJob, GpuCost)> {
    // No halo: `halo_cells` is zero for Marching Cubes, which marches every cell
    // and puts vertices on grid edges both chunks compute identically.
    let params = GridParams::new([layout.cells() + 1; 3], origin, layout.cell_size()).ok()?;

    // The bare terrain, with **no tape**: the log is what the device folds, so
    // baking it into the base would be the stale cache this design does not have.
    let started = Instant::now();
    let mut uploaded = 0u64;
    let base = match fields.bases.entry(id) {
        std::collections::hash_map::Entry::Occupied(slot) => slot.into_mut(),
        std::collections::hash_map::Entry::Vacant(slot) => {
            let sampled = FieldBuffer::sampled(device, queue, params, &Ground).ok()?;
            uploaded = params.field_buffer_size();
            slot.insert(sampled)
        }
    };
    fields.sample_bytes += uploaded;
    // Zero on a hit, and the HUD says `(cached)` when it is.
    let base_ms = if uploaded == 0 {
        0.0
    } else {
        started.elapsed().as_secs_f64() * 1000.0
    };

    // The pruned survivors as the shader's 64-byte records, into the reused
    // scratch. The `SmoothAdd` arm is unreachable here -- this demo pushes only
    // `Brush::add` and `Brush::subtract` -- and is spelled out for `BrushOp`'s
    // exhaustiveness rather than swallowed by a wildcard that would silently
    // mistranslate a brush kind added later.
    fields.log.clear();
    fields.log.extend(survivors.iter().map(|b| GpuBrush {
        shape: GpuShape::Sphere {
            center: b.shape.center,
            radius: b.shape.radius,
        },
        op: match b.op {
            BrushOp::Add => GpuOp::Add,
            BrushOp::Subtract => GpuOp::Subtract,
            BrushOp::SmoothAdd { k } => GpuOp::SmoothAdd { k: k as f32 },
        },
    }));

    let started = Instant::now();
    let buffer = FieldBuffer::new(device, params);
    fields
        .sampler
        .fold_into(device, queue, base, &buffer, &fields.log)
        .ok()?;
    let fold_ms = started.elapsed().as_secs_f64() * 1000.0;

    let geometry = mc
        .extract_indirect(device, queue, &buffer, GPU_TRIANGLE_BUDGET)
        .ok()?;
    // The staging copy is recorded in the same submission, so the geometry
    // buffers may be dropped at the end of this function: wgpu keeps a
    // submission's resources alive until it completes.
    let bytes = u64::from(GPU_TRIANGLE_BUDGET) * 9 * 4;
    let readback = read_bytes_many_deferred(
        device,
        queue,
        &[
            (&geometry.total, 4),
            (&geometry.positions, bytes),
            (&geometry.normals, bytes),
        ],
    )
    .ok()?;
    Some((
        GpuJob { id, readback },
        GpuCost {
            base_ms,
            fold_ms,
            timings: geometry.timings,
        },
    ))
}

/// Finish any GPU chunk whose bytes have arrived.
///
/// One frame of latency at best, more under load, and that is the honest shape of
/// a GPU readback -- the alternative is the blocking wait that deadlocks a
/// browser tab. Nothing else in the frame waits on it: the chunk keeps its
/// previous mesh until the new one lands.
fn gpu_collect(
    mut world: ResMut<World>,
    mut pending: ResMut<GpuPending>,
    device: Res<RenderDevice>,
    mut commands: Commands,
    mut meshes: ResMut<Assets<Mesh>>,
    material: Res<SurfaceMaterial>,
    chunks: Query<(Entity, &Chunk)>,
) {
    if pending.jobs.is_empty() {
        return;
    }
    let wgpu_device = device.wgpu_device();
    // Drain the ready jobs and keep the rest, in order. `retain` cannot do it:
    // the ready ones have to be *moved* out, because `Readback::take` consumes.
    let held = std::mem::take(&mut pending.jobs);
    for job in held {
        if !job.readback.ready(wgpu_device) {
            pending.jobs.push(job);
            continue;
        }
        let id = job.id;
        let Ok(parts) = job.readback.take() else {
            // A refused map is a device problem, not a chunk with no surface, so
            // the chunk keeps the mesh it has and the error is on screen.
            world.gpu_failures += 1;
            continue;
        };
        let Some((total, mesh)) = gpu_mesh(&parts) else {
            // Three buffers were asked for; anything else is a bug in this file,
            // not a device problem.
            world.gpu_failures += 1;
            continue;
        };
        world.gpu_peak_triangles = world.gpu_peak_triangles.max(total);
        reconcile(
            &mut commands,
            &chunks,
            &material,
            id,
            mesh.map(|m| meshes.add(m)),
        );
    }
}

/// Turn one finished readback into its triangle count and a `Mesh`, or `None`
/// when the readback was not the three buffers [`gpu_dispatch`] asked for.
///
/// The count is the **unclamped** total and comes back even when the mesh does
/// not, so the HUD can show a chunk that hit [`GPU_TRIANGLE_BUDGET`]. A count of
/// zero is a chunk with no surface and its `Mesh` is `None`, which
/// [`reconcile`] turns into a despawn -- the same as a CPU `None`.
///
/// `extract_indirect` writes a **triangle soup with no index buffer** -- three
/// `f32` positions per vertex, three vertices per triangle, in cell order -- so
/// the index buffer here is the identity, exactly as `gpu_compute_mc.rs` builds
/// it.
fn gpu_mesh(parts: &[Vec<u8>]) -> Option<(u32, Option<Mesh>)> {
    let [total, positions, normals] = parts else {
        return None;
    };
    let total = u32::from_le_bytes(total.get(..4)?.try_into().ok()?);
    // Clamped, because the buffers were sized from `GPU_TRIANGLE_BUDGET` and the
    // total is the *unclamped* count -- that is the contract `extract_indirect`
    // documents, and reading past the budget would read uninitialised bytes. The
    // HUD prints the raw total beside the budget, so truncation is visible rather
    // than silent.
    let triangles = total.min(GPU_TRIANGLE_BUDGET) as usize;
    if triangles == 0 {
        return Some((total, None));
    }
    let vertices = triangles * 3;
    // Three `f32` per vertex, so twelve bytes each. `as_chunks` drops a trailing
    // partial vertex the same way `chunks_exact` did; `take` bounds it to the
    // triangles the total actually claims, because the buffer is sized to the
    // budget and the rest of it is uninitialised.
    let floats = |bytes: &[u8]| -> Vec<[f32; 3]> {
        bytes
            .as_chunks::<12>()
            .0
            .iter()
            .take(vertices)
            .map(|v| {
                let f = |at: usize| f32::from_le_bytes([v[at], v[at + 1], v[at + 2], v[at + 3]]);
                [f(0), f(4), f(8)]
            })
            .collect()
    };
    let mut mesh = Mesh::new(
        PrimitiveTopology::TriangleList,
        RenderAssetUsages::default(),
    );
    mesh.insert_attribute(Mesh::ATTRIBUTE_POSITION, floats(positions));
    mesh.insert_attribute(Mesh::ATTRIBUTE_NORMAL, floats(normals));
    // The shader never reads `in.uv` -- the terrain is triplanar, projected from
    // world position -- but `pbr_input_from_standard_material` is compiled with
    // `VERTEX_UVS` when the mesh carries them, and the CPU chunks carry them
    // because `MeshBuilder` emits them. Matching the vertex layout is what lets
    // one pipeline serve both.
    mesh.insert_attribute(Mesh::ATTRIBUTE_UV_0, vec![[0.0f32, 0.0]; vertices]);
    mesh.insert_indices(Indices::U32((0..vertices as u32).collect()));
    Some((total, Some(mesh)))
}

/// Mark every chunk in the sandbox dirty, and return how many that is.
///
/// Two callers -- `setup`'s startup fill and `switch_algorithm` -- and both need
/// the count for the HUD. One triple loop, because a second one is a second
/// place for [`EXTENT`] to go stale.
fn mark_sandbox(dirty: &mut DirtySet) -> usize {
    for z in 0..EXTENT[2] {
        for y in 0..EXTENT[1] {
            for x in 0..EXTENT[0] {
                dirty.insert(ChunkId::new([x, y, z]));
            }
        }
    }
    (EXTENT[0] * EXTENT[1] * EXTENT[2]) as usize
}

/// Show the modal exactly while an extractor switch has a backlog.
///
/// Not while a dig drains: one to eight chunks clear inside a single frame's
/// budget, and a modal that flashes on every carve is worse than none. Written
/// only when it differs -- Bevy's UI extraction is change-driven, so an
/// unconditional write turns a static overlay into per-frame work.
fn loading_modal(world: Res<World>, mut modal: Query<&mut Visibility, With<LoadingModal>>) {
    let wanted = if world.switching {
        Visibility::Inherited
    } else {
        Visibility::Hidden
    };
    for mut visibility in &mut modal {
        if *visibility != wanted {
            *visibility = wanted;
        }
    }
}

/// `1`-`8` re-mesh the whole sandbox with a different mesher.
///
/// The whole sandbox, not the dirty set: nothing about the *field* changed, so
/// `mark_edit` would correctly report zero dirty chunks. And re-meshing 256
/// chunks with each algorithm in turn is the comparison -- the HUD prints what
/// it cost, which is the only honest way to show that Marching Tetrahedra is
/// three times the triangles.
///
/// Marking, not meshing: [`drain_dirty`] does the work over as many frames as
/// the budget needs, and `world.switching` is what puts the modal on screen
/// while it does. That is the difference between a switch that locks the demo
/// and one you can watch fill in.
///
/// `handle_keys` in the shared harness also reads `Digit1`-`Digit7`, into
/// `ViewFlags::field`. This example never reads that field, so there is nothing
/// to collide with, and `Digit8` is bound nowhere else at all.
fn switch_algorithm(keys: Res<ButtonInput<KeyCode>>, mut world: ResMut<World>) {
    const DIGITS: [KeyCode; 8] = [
        KeyCode::Digit1,
        KeyCode::Digit2,
        KeyCode::Digit3,
        KeyCode::Digit4,
        KeyCode::Digit5,
        KeyCode::Digit6,
        KeyCode::Digit7,
        KeyCode::Digit8,
    ];
    let Some(index) = DIGITS.iter().position(|k| keys.just_pressed(*k)) else {
        return;
    };
    let chosen = Algorithm::ALL[index];
    if chosen == world.algorithm {
        return;
    }
    world.algorithm = chosen;
    world.switch_chunks = mark_sandbox(&mut world.dirty);
    world.switching = true;
    world.switch_started = Instant::now();
    world.last_touched.clear();
    // The peak belongs to the mesher that produced it, so a switch starts it
    // again rather than showing key 8's figure while key 1 is selected.
    world.gpu_peak_triangles = 0;
}

/// Reduce this frame's touches to the same numbers the keyboard and mouse
/// produce, and reveal the on-screen buttons the first time a finger appears.
///
/// Ordered before [`move_camera`] and [`dig`], which are where those numbers are
/// consumed. **Each touch is classified by where it started**, never by where it
/// is now, so a gesture keeps its role for its whole life: a look sweep that
/// wandered into the left third would otherwise start walking mid-drag.
fn touch_input(
    touches: Res<Touches>,
    windows: Query<&Window, With<PrimaryWindow>>,
    presses: Query<(Ref<Interaction>, &TouchButton)>,
    mut nodes: Query<&mut Node, With<TouchButton>>,
    mut touch: ResMut<TouchIntent>,
    mut world: ResMut<World>,
) {
    let Ok(window) = windows.single() else {
        return;
    };
    // Cleared every frame, so a lifted finger stops the body rather than leaving
    // it walking into a wall for ever.
    touch.look = Vec2::ZERO;
    touch.move_axis = Vec2::ZERO;
    touch.dig = false;
    touch.fill = false;
    touch.jump = false;

    /// Fraction of the screen width that is the movement stick.
    ///
    /// The left 40%: a thumb's whole arc on a phone held in two hands, and it
    /// leaves the majority of the screen -- the part actually being looked at --
    /// to the look drag.
    const STICK_ZONE: f32 = 0.4;
    let split = window.width() * STICK_ZONE;
    for t in touches.iter() {
        if t.start_position().x < split {
            touch.move_axis += touch_axes(t.start_position(), t.position());
        } else {
            // `delta`, this frame's movement, not `distance` from the start:
            // this stands in for a mouse, whose motion is also per-frame. It is
            // also what keeps a thumb parked on `DIG` from spinning the camera --
            // a stationary finger has a delta of zero, so the buttons need no
            // exclusion rectangle.
            touch.look += t.delta();
        }
    }
    for (interaction, button) in &presses {
        if *interaction != Interaction::Pressed {
            continue;
        }
        match button {
            // Level-triggered, exactly like a held mouse button: `dig`'s own
            // pacing decides the rate.
            TouchButton::Dig => touch.dig = true,
            TouchButton::Fill => touch.fill = true,
            // Edge-triggered. `bevy_ui`'s focus system assigns `Pressed` only on
            // the transition into it, so `Ref::is_changed` is the press edge --
            // and without it a thumb left on the button bunny-hops, which
            // `Space` held does not do.
            TouchButton::Jump => touch.jump |= interaction.is_changed(),
        }
    }
    if touches.any_just_pressed() && !touch.seen {
        touch.seen = true;
        // `dig`'s edit gate requires `grabbed`, which on a desktop is what the
        // first left click sets. A phone has no pointer to lock, so the first
        // touch is the same event.
        world.grabbed = true;
        // Flipped once, not every frame: `Node` is change-detected and
        // re-writing it would re-lay-out the UI on every frame of the demo.
        for mut node in &mut nodes {
            node.display = Display::Flex;
        }
    }
}

fn grab(
    keys: Res<ButtonInput<KeyCode>>,
    buttons: Res<ButtonInput<MouseButton>>,
    // `CursorOptions` is its own component on the window entity in Bevy 0.19,
    // not a field of `Window`.
    mut cursor: Query<&mut CursorOptions, With<PrimaryWindow>>,
    mut world: ResMut<World>,
    touch: Res<TouchIntent>,
) {
    let Ok(mut cursor) = cursor.single_mut() else {
        return;
    };
    if buttons.just_pressed(MouseButton::Left) && !world.grabbed {
        world.grabbed = true;
    }
    if keys.just_pressed(KeyCode::Tab) {
        world.grabbed = !world.grabbed;
    }
    // A touch device gets no pointer lock: iOS Safari does not implement the
    // API, and asking for it produces a console error and nothing else. `dig`
    // still needs `world.grabbed`, which `touch_input` sets on the first touch,
    // so the gate here is only about the cursor.
    let (mode, visible) = if world.grabbed && !touch.seen {
        (CursorGrabMode::Locked, false)
    } else {
        (CursorGrabMode::None, true)
    };
    cursor.grab_mode = mode;
    cursor.visible = visible;
}

/// One frame of the body's vertical motion.
///
/// **Gravity is not integrated while the body is standing**, and that is not a
/// shortcut. Integrating it makes a resting body sink `g·dt²` -- five
/// millimetres a frame -- which [`resolve_body`] then pushes back out along the
/// surface *normal*, and on a slope that normal has a horizontal component: a
/// body standing still slides downhill for ever, accumulating speed. Horizontal
/// motion in this demo is direct rather than integrated, so `velocity` carries
/// only the fall, and a body on the ground has none to carry.
///
/// The `<= 0.0` guard is what keeps a jump: `Space` writes `+JUMP_SPEED` on a
/// frame where `grounded` is still true, and this must not erase it.
///
/// Gravity resumes the moment the ground probe stops answering, which is what
/// walking off a ledge is.
fn gravity_step(grounded: bool, velocity: &mut Vec3, dt: f32) {
    if grounded && velocity.y <= 0.0 {
        velocity.y = 0.0;
    } else {
        velocity.y -= GRAVITY * dt;
    }
}

/// Push the body out of the rock, and report whether it is standing on it.
///
/// Four spheres against the field, [`RESOLVE_PASSES`] times. `f` is not an exact
/// distance -- `Ground` has `|grad f|` between 1 and 1.207 -- and it errs the
/// safe way for a height field: `f >= d`, so the test triggers only on a real
/// overlap and at worst lets the body sink 18% of a radius into the steepest
/// slope, which is 7 cm and invisible. Inside a carved region the field *is* the
/// brush's exact distance, so there is no error at all where it matters.
///
/// The push direction is the normalised gradient, which points away from the
/// solid. **M-172 measured `BrushStack::gradient` returning exactly
/// `[0, 0, 0]` on the medial axis**, so below [`GRADIENT_EPS`] the direction is
/// `Vec3::Y`: on the medial axis inside rock, up is the way out.
///
/// **No friction and no slope limit**, stated because both are choices. The push
/// is pure depenetration along the normal, so a body walking into a wall slides
/// along it and a body walking into a slope climbs it however steep it is. In a
/// digging sandbox that is the behaviour you want: dig a pit with vertical walls
/// and you can still walk out of it. A resting body does not slide, because
/// [`gravity_step`] gives it nothing to slide with.
///
/// A free function rather than a closure inside the system, so the field is
/// built once and so this can be tested without an `App`.
fn resolve_body(field: &impl Sdf<Scalar = f32>, eye: &mut Vec3, velocity: &mut Vec3) -> bool {
    // The **deepest** overlap, not the last one looked at. This was
    // `contact = Some(n)` inside the loop, so the velocity cancellation below
    // used whichever sphere happened to come last in `BODY_OFFSETS` -- a body
    // wedged with its head in a lip and its feet on the floor cancelled against
    // the head. Deeper is more of the reason the body is stuck, so it is the
    // direction worth cancelling along.
    let mut contact: Option<(f32, Vec3)> = None;
    for _ in 0..RESOLVE_PASSES {
        for offset in BODY_OFFSETS {
            let c = *eye - Vec3::Y * offset;
            let f = field.sample([c.x, c.y, c.z]);
            if f >= BODY_RADIUS {
                continue;
            }
            let g = Vec3::from_array(field.gradient([c.x, c.y, c.z]));
            let n = if g.length() > GRADIENT_EPS {
                g.normalize()
            } else {
                Vec3::Y
            };
            let depth = BODY_RADIUS - f;
            *eye += n * depth;
            if contact.is_none_or(|(deepest, _)| depth > deepest) {
                contact = Some((depth, n));
            }
        }
    }
    // Otherwise gravity accumulates against the floor and the first step off a
    // ledge is a plummet. Only the component *into* the surface goes: sliding
    // along it is movement the player asked for.
    if let Some((_, n)) = contact {
        let into = velocity.dot(n);
        if into < 0.0 {
            *velocity -= n * into;
        }
    }
    // A separate downward probe, not `contact.is_some()`: a body pressed against
    // a wall is in contact and is not standing on anything. The lowest sphere is
    // the one with the largest offset, found rather than indexed so the
    // contingency of adding another offset cannot silently probe the wrong one.
    //
    // **A cross across the foot sphere's own lower surface**, and each of the
    // three geometries this could be is a different bug:
    //
    // * *One point* straight down reads air whenever the body is wedged with its
    //   centre over a void and part of its foot resting on rock to one side --
    //   `move_camera`'s `world.grounded && Space` then refuses the jump, and the
    //   hole the player dug is a trap. That was the shipped behaviour.
    // * A *flat* cross, five samples at one depth, reads solid whenever the
    //   terrain beside the foot is higher than the terrain under it -- which on
    //   any slope is ground the body is not touching. `gravity_step` cuts gravity
    //   and the body hangs in mid-air. Measured on this field at the origin,
    //   where the height rises 0.63 per unit: **12.5 cm of hover**.
    // * This one: the lateral samples sit at the sphere's own boundary for their
    //   offset, `sqrt(R² - r²)` below the centre, and only then drop by
    //   [`GROUND_PROBE`]. So the tolerance is purely vertical, which is what its
    //   doc comment says it is, and every sample asks the same question -- "is
    //   there rock just under my foot, here?" -- at a point the foot really does
    //   hang over.
    //
    // `0.7` of the radius laterally, so the four outer samples cover most of the
    // footprint while staying `0.71 R` below the centre rather than degenerating
    // to the sphere's equator, where "below" stops being the question.
    const FOOTPRINT: f32 = 0.7;
    let lowest = BODY_OFFSETS.into_iter().fold(f32::MIN, f32::max);
    let foot = *eye - Vec3::Y * lowest;
    let r = BODY_RADIUS * FOOTPRINT;
    let edge = BODY_RADIUS * (1.0 - FOOTPRINT * FOOTPRINT).sqrt() + GROUND_PROBE;
    let centre = BODY_RADIUS + GROUND_PROBE;
    [
        [0.0, -centre, 0.0],
        [r, -edge, 0.0],
        [-r, -edge, 0.0],
        [0.0, -edge, r],
        [0.0, -edge, -r],
    ]
    .into_iter()
    .any(|[dx, dy, dz]| field.sample([foot.x + dx, foot.y + dy, foot.z + dz]) <= 0.0)
}

/// Move the camera. Flying is a direct write; walking integrates gravity and is
/// stopped by the field.
///
/// One system, not two, because both write the same `Transform` and Bevy does not
/// order systems that merely touch the same component -- two of them would race
/// for the camera every frame. `F` switches; the mouse look is shared.
fn move_camera(
    keys: Res<ButtonInput<KeyCode>>,
    motion: Res<AccumulatedMouseMotion>,
    time: Res<Time>,
    mut world: ResMut<World>,
    mut look: ResMut<Look>,
    mut camera: Query<&mut Transform, With<Camera3d>>,
    touch: Res<TouchIntent>,
) {
    let Ok(mut transform) = camera.single_mut() else {
        return;
    };
    if keys.just_pressed(KeyCode::KeyF) {
        world.walking = !world.walking;
        // Zeroed on the switch, so a fall in progress does not carry into fly
        // mode and shoot the camera through the floor when walking resumes.
        world.velocity = Vec3::ZERO;
    }
    if world.grabbed {
        let sensitivity = 0.0022;
        look.yaw -= motion.delta.x * sensitivity;
        look.pitch = (look.pitch - motion.delta.y * sensitivity).clamp(-1.5, 1.5);
    }
    // Touch look is applied unconditionally, outside the `grabbed` gate above:
    // `grabbed` means the pointer is locked, and there is no pointer to lock on
    // a phone. Its own sensitivity, and higher, because a thumb drag across a
    // handset is a fraction of a mouse sweep.
    if touch.look != Vec2::ZERO {
        let sensitivity = 0.0035;
        look.yaw -= touch.look.x * sensitivity;
        look.pitch = (look.pitch - touch.look.y * sensitivity).clamp(-1.5, 1.5);
    }
    transform.rotation = Quat::from_euler(EulerRot::YXZ, look.yaw, look.pitch, 0.0);

    let mut direction = Vec3::ZERO;
    let forward = *transform.forward();
    let right = *transform.right();
    for (key, delta) in [
        (KeyCode::KeyW, forward),
        (KeyCode::KeyS, -forward),
        (KeyCode::KeyA, -right),
        (KeyCode::KeyD, right),
    ] {
        if keys.pressed(key) {
            direction += delta;
        }
    }
    // The virtual stick, in the same basis the keys use, so the flatten and the
    // renormalise below apply to it too and a phone walks at exactly the speed a
    // keyboard does.
    direction += right * touch.move_axis.x + forward * touch.move_axis.y;
    if world.walking {
        // Flattened and renormalised, so looking at your feet does not slow you
        // down. `Vec3::normalize_or_zero` covers looking straight up, where the
        // flattened forward is the zero vector.
        direction.y = 0.0;
        direction = direction.normalize_or_zero();
    } else {
        // Vertical flight, and only in fly mode: walking has `V` for up and
        // gravity for down.
        for (key, delta) in [(KeyCode::KeyQ, -Vec3::Y), (KeyCode::KeyE, Vec3::Y)] {
            if keys.pressed(key) {
                direction += delta;
            }
        }
        direction = direction.normalize_or_zero();
    }
    // 9.0, not the old 6.0: a 16-unit box takes 6.4 s to cross at the walk speed
    // and 1.8 s at this one.
    let speed = if keys.pressed(KeyCode::ShiftLeft) {
        9.0
    } else {
        2.5
    };
    let dt = time.delta_secs();
    if !world.walking {
        transform.translation += direction * speed * dt;
        return;
    }

    gravity_step(world.grounded, &mut world.velocity, dt);
    // `Space`, which is what a hand reaches for, and the on-screen `JUMP` beside
    // it -- `touch_input` has already reduced that to one edge per tap. The
    // shared harness also reads `Space` into `ViewFlags::paused`, and that is
    // harmless here: `paused` is read only by `orbit_camera`, and `setup` takes
    // `OrbitCamera` off this camera so that system's query is empty in this demo.
    if world.grounded && (keys.just_pressed(KeyCode::Space) || touch.jump) {
        world.velocity.y = JUMP_SPEED;
    }
    transform.translation += (direction * speed + Vec3::Y * world.velocity.y) * dt;
    // Sampled against the field rather than a collider, which is the only option
    // that cannot go stale when the player digs: M-116 puts a convex
    // decomposition at 241-272 ms per fragment, and this demo edits the field on
    // every frame of a held button, so a cache is invalid before it is built.
    let World {
        brushes,
        velocity,
        grounded,
        ..
    } = &mut *world;
    let field = BrushStack {
        base: Ground,
        brushes,
    };
    *grounded = resolve_body(&field, &mut transform.translation, velocity);
    // The walls are the boundary and this is what makes them one. `Ground` has a
    // height at every `x` and `z`, so without this a walk past `x = 8` stands on
    // invisible ground outside the box, in front of the wall the player just
    // walked through.
    //
    // Walk mode only: fly mode is how the sandbox is inspected from outside, and
    // `aim` already refuses to carve past `sandbox`. No `y` clamp either -- the
    // field stops the descent and the floor slab is scenery.
    let (lo, hi) = sandbox(&world.layout);
    transform.translation.x = transform
        .translation
        .x
        .clamp(lo.x + BODY_RADIUS, hi.x - BODY_RADIUS);
    transform.translation.z = transform
        .translation
        .z
        .clamp(lo.z + BODY_RADIUS, hi.z - BODY_RADIUS);
}

/// The loop this example exists for: one brush, one incremental re-mesh, for as
/// long as the button is held.
#[allow(clippy::too_many_arguments)]
fn dig(
    buttons: Res<ButtonInput<MouseButton>>,
    scroll: Res<AccumulatedMouseScroll>,
    keys: Res<ButtonInput<KeyCode>>,
    time: Res<Time>,
    mut world: ResMut<World>,
    target: Res<Aim>,
    mut auto: ResMut<AutoCarve>,
    capture: Res<Capture>,
    touch: Res<TouchIntent>,
) {
    // Multiplicative, because radius is a scale: a fixed 0.1 step is a third of
    // the smallest brush and a twentieth of the largest. The wheel is free in
    // this example specifically -- `setup` strips `OrbitCamera` off the harness
    // camera, so `orbit_camera`'s query is empty here and nothing else reads the
    // scroll.
    //
    // `AccumulatedMouseScroll` carries a `unit` because the number's meaning
    // depends on it, and ignoring that is not a simplification but a bug on the
    // platform this demo actually runs on: a browser reports `Pixel` with a
    // `deltaY` around 100 per notch, and 1.12^100 clamps the brush to its
    // maximum on the first flick. Convert to notches once and there is a single
    // scale law for every input device.
    const PIXELS_PER_NOTCH: f32 = 100.0;
    let notches = match scroll.unit {
        MouseScrollUnit::Line => scroll.delta.y,
        MouseScrollUnit::Pixel => scroll.delta.y / PIXELS_PER_NOTCH,
    };
    let mut scale = 1.0;
    if keys.just_pressed(KeyCode::BracketLeft) {
        scale *= 0.85;
    }
    if keys.just_pressed(KeyCode::BracketRight) {
        scale *= 1.18;
    }
    if notches != 0.0 {
        scale *= 1.12_f32.powf(notches);
    }
    if scale != 1.0 {
        world.radius = (world.radius * scale).clamp(0.10, 2.00);
    }
    if keys.just_pressed(KeyCode::KeyC) {
        world.show_chunks = !world.show_chunks;
    }

    // **One carve per *captured* frame, not per rendered frame.** The tunnel
    // bores away from the camera at 0.34 units a step, so at 60 Hz it is out of
    // frame in a quarter of a second -- and a recording cannot begin until the
    // window has stopped resizing, which is 30 frames in. A clip made that way
    // photographs the aftermath while all the action happens off-screen. Every
    // other capture-driven example in this directory paces itself off
    // `capture.taken` for the same reason.
    let ready = !capture.is_active() || capture.taken > 0;
    let due = !capture.is_active() || auto.step <= capture.taken / auto.every.max(1);
    let scripted = if auto.remaining > 0 && ready && due {
        auto.remaining -= 1;
        let centre = AutoCarve::centre(auto.step);
        auto.step += 1;
        Some(centre)
    } else {
        None
    };
    // Both gated on a non-empty log, so the union of the doomed brushes' boxes
    // below is never the union of nothing, and both keys are a documented no-op
    // when there is nothing to undo.
    let clear = keys.just_pressed(KeyCode::KeyX) && !world.brushes.is_empty();
    let undo = keys.just_pressed(KeyCode::KeyZ) && !world.brushes.is_empty();
    let shrink = clear || undo;
    // The on-screen buttons feed the same two booleans the mouse does, so a
    // phone carves through the identical pacing, distance gate and edit path --
    // there is no second edit route to keep in step with this one.
    let left = buttons.pressed(MouseButton::Left) || touch.dig;
    let right = buttons.pressed(MouseButton::Right) || touch.fill;
    if !left && !right {
        world.stroke_last = None;
    }
    // Hold to keep editing, but only where the brush has actually moved. A
    // second subtract sphere at the same centre is idempotent -- the mesh does
    // not change -- so a timer would push a brush every frame for no visible
    // reason, and this demo's own HUD advertises that every field sample walks
    // the whole log. Half a radius leaves no unbroken rock between consecutive
    // spheres, so a sweep cuts one continuous tunnel.
    // Wall-clock pacing on top of the distance gate. The distance gate admits
    // one edit per *frame* once the aim point has moved half a radius, so the
    // log grew at frame rate and a four-second hold pushed ~200 brushes -- which
    // every later field sample then walks. This bounds the rate instead.
    //
    // No catch-up loop, and no subtracting `EDIT_PERIOD`: catching up is exactly
    // the burst this exists to prevent.
    world.stroke_clock += time.delta_secs();
    let due = world.stroke_clock >= EDIT_PERIOD;
    let moved = world
        .stroke_last
        .is_none_or(|last| target.point.distance(last) >= world.radius * 0.5);
    let editable = world.grabbed && target.hit && moved && due;
    // One edit per frame, and a key that shrinks the log wins over a held
    // button. They drive the same `split` index from opposite ends, so letting
    // both run in one frame would truncate the log *and* count a stroke edit for
    // a brush that was never pushed. `!left` on the fill arm makes left win when
    // both buttons are held, so a stuck button cannot alternate.
    let carve = !shrink && (scripted.is_some() || (editable && left));
    let fill = !shrink && editable && right && !left;
    if !(carve || fill || shrink) {
        return;
    }

    let started = Instant::now();
    let layout = world.layout;
    let algorithm = world.algorithm;
    // One aim point, two consumers: the ghost drew this exact sphere this exact
    // frame. The scripted path keeps its own hardcoded centres so the committed
    // capture sequence is unchanged.
    let centre = match scripted {
        Some(centre) => centre,
        None => target.point,
    };

    // How much of the log survives. A push grows it by one and compares against
    // the log without that one; `Z` drops the last brush and `X` drops all of
    // them. Three keys, one index.
    let split = if clear {
        0
    } else {
        world.brushes.len() - usize::from(undo)
    };

    // The region to re-check.
    let (min_cell, max_cell) = if shrink {
        // The doomed brushes' own padded boxes, unioned -- read off the
        // *brushes*, never off `world.radius`, because the wheel may have
        // resized the brush since any of them was pushed and re-checking the
        // wrong box leaves a stale chunk on screen.
        //
        // Not the whole sandbox: outside this box `before` and `after` agree by
        // construction, and the sandbox is 14x the cells for the same answer --
        // with every one of those samples walking the whole log twice.
        let cell = layout.cell_size();
        let mut lo = [f32::INFINITY; 3];
        let mut hi = [f32::NEG_INFINITY; 3];
        for brush in &world.brushes[split..] {
            let reach = brush.shape.radius + cell;
            for axis in 0..3 {
                lo[axis] = lo[axis].min(brush.shape.center[axis] - reach);
                hi[axis] = hi[axis].max(brush.shape.center[axis] + reach);
            }
        }
        (layout.cell_of(lo), layout.cell_of(hi))
    } else {
        let shape = Sphere {
            center: [centre.x, centre.y, centre.z],
            radius: world.radius,
        };
        let brush = if carve {
            Brush::subtract(shape)
        } else {
            Brush::add(shape)
        };
        world.brushes.push(brush);

        // Padded by a cell, which is not tidiness: `cell_of` inverts
        // `world_of_sample` in a cell's interior and not reliably on its corner,
        // for the same power-of-two reason as M-32. A padded range cannot lose a
        // crossing to that; an exact one can.
        let reach = world.radius + layout.cell_size();
        (
            layout.cell_of([centre.x - reach, centre.y - reach, centre.z - reach]),
            layout.cell_of([centre.x + reach, centre.y + reach, centre.z + reach]),
        )
    };

    // `before` and `after` are two slices of the one log, which is what makes
    // this exact: the two fields differ by precisely the brushes between them.
    // A push reads `[..split]` then `[..]`; a shrink reads them the other way
    // round and drops the tail *after* the re-mesh has sampled it.
    //
    // Truncating first would be the bug: `before` and `after` would both be bare
    // `Ground`, `mark_edit` would correctly report nothing changed, and the
    // carve would stay on screen beside an empty log.
    let (report, touched) = {
        // `world.dirty` is the one queue -- marking into it is all `dig` does
        // now, and `drain_dirty` meshes it later in the same frame under its
        // budget. Destructured because `dirty` is written while `brushes` is
        // read, and the borrow checker needs to see the two fields as disjoint.
        let World { dirty, brushes, .. } = &mut *world;
        let (before_log, after_log) = if shrink {
            (&brushes[..], &brushes[..split])
        } else {
            (&brushes[..split], &brushes[..])
        };
        // The same bound as `drain_dirty`'s, over the edit box instead of a
        // chunk. `mark_edit` samples both fields on every cell corner of
        // `[min_cell, max_cell]`, which is 1,024 evaluations of a tape of length
        // `L` -- 2-4 ms on its own at `L = 400`. Pruning is bit-exact, so E1 and
        // the dirty set are unchanged, which matters because they are the numbers
        // this example exists to report.
        //
        // The longest axis, not axis 0: the box is a cube for a push (`reach` is
        // the same on every axis) but a shrink unions several brushes' boxes and
        // need not be. `ChunkBox` encloses a cube, and a cube that covers the
        // longest axis covers the box -- conservative, so it can only keep a
        // brush it could have dropped.
        let widest = (0..3)
            .map(|axis| max_cell[axis] - min_cell[axis] + 1)
            .max()
            .expect("three axes");
        let edit_box = ChunkBox::new(
            layout.world_of_sample(min_cell),
            widest as f32 * layout.cell_size(),
        );
        // Two `Vec`s and not `world.survivors`: this runs once per edit, not once
        // per chunk, so the allocation is not on the hot path -- and the two folds
        // have to exist at the same time, because differing by one brush is the
        // whole mechanism.
        let mut before_kept = Vec::new();
        let mut after_kept = Vec::new();
        prune_into(before_log, &Ground, edit_box, &mut before_kept);
        prune_into(after_log, &Ground, edit_box, &mut after_kept);
        let before = BrushStack {
            base: Ground,
            brushes: &before_kept,
        };
        let after = BrushStack {
            base: Ground,
            brushes: &after_kept,
        };
        let report = mark_edit(&layout, &before, &after, min_cell, max_cell, dirty)
            .expect("a dig brush spans a few cells, far inside the u32 sample space");
        // `mark_edit` names the owner of each changed cell; with a halo, the
        // chunks one step *negative* on any axis whose first cell plane the box
        // reached borrow that plane too, so their meshes went stale as well.
        //
        // The dependency is a product over axes, not three faces: a cell at
        // global `c` is borrowed by chunk `c/16 - 1` on every axis where
        // `c % 16 == 0`, so a cell on a chunk corner is borrowed by all seven
        // chunks in the negative octant. Face-only expansion would leave a
        // diagonal neighbour stale -- a one-cell notch at a chunk corner, which
        // is the exact defect the halo exists to remove.
        //
        // Conservative in one direction only: the box is the *padded* brush box,
        // so this can over-dirty by a chunk and can never under-dirty.
        //
        // Collected before the expansion, because `DirtySet::insert` cannot run
        // inside an iteration over the same set.
        let touched_owners: Vec<ChunkId> = dirty.iter().collect();
        if algorithm.halo_cells() > 0 {
            let cells = i64::from(CHUNK_CELLS);
            for id in &touched_owners {
                let reaches =
                    [0, 1, 2].map(|axis| min_cell[axis] <= i64::from(id.coords[axis]) * cells);
                for dz in 0..=1 {
                    for dy in 0..=1 {
                        for dx in 0..=1 {
                            let step = [dx, dy, dz];
                            if step == [0, 0, 0]
                                || (0..3).any(|axis| step[axis] == 1 && !reaches[axis])
                            {
                                continue;
                            }
                            let mut coords = id.coords;
                            for axis in 0..3 {
                                coords[axis] -= step[axis];
                            }
                            dirty.insert(ChunkId::new(coords));
                        }
                    }
                }
            }
        }
        // From the *expanded* set, so the outline gizmo shows what was actually
        // re-meshed. `EditReport::dirty_chunks` is computed inside `mark_edit`
        // and is unaffected, so the HUD's E1 line keeps reporting the owners it
        // always did.
        (report, dirty.iter().collect::<Vec<ChunkId>>())
    };
    if shrink {
        world.brushes.truncate(split);
        // The next hold starts a fresh stroke rather than measuring against a
        // brush that no longer exists.
        world.stroke_last = None;
    }

    // Marking, so this is the mark cost -- `mark_edit` plus the halo expansion,
    // which is the E1 measurement this example exists to report. The extraction
    // it queues is charged to `drain_dirty`'s budget instead, which is the whole
    // point: an edit can no longer make a frame arbitrarily long.
    world.last_edit_ms = started.elapsed().as_secs_f64() * 1000.0;
    // One line per scripted edit, so the log-growth claim in the module docs can
    // be checked from a terminal rather than by reading a HUD off a screenshot.
    if scripted.is_some() {
        info!(
            "edit {:>3}: {} chunks in {:.3} ms, E1 {:.1}%",
            world.brushes.len(),
            touched.len(),
            world.last_edit_ms,
            100.0 * report.changed_fraction()
        );
    }
    world.last_edit = Some(report);
    world.last_chunks = touched.len();
    world.last_touched = touched;
    // The stroke is the "how quick is it" number the held button exists to show:
    // one hold, N edits, milliseconds each. Scripted carves are excluded because
    // they are paced off the capture, not off a hand.
    if (carve || fill) && scripted.is_none() {
        if world.stroke_last.is_none() {
            world.stroke_edits = 0;
            world.stroke_ms = 0.0;
        }
        world.stroke_edits += 1;
        world.stroke_ms += world.last_edit_ms;
        world.stroke_last = Some(centre);
        world.stroke_clock = 0.0;
    }
}

fn report(
    world: Res<World>,
    fields: Res<GpuFields>,
    mut stats: ResMut<DemoStats>,
    meshes: Res<Assets<Mesh>>,
    chunks: Query<(&Chunk, &Mesh3d)>,
) {
    // Totals across every resident chunk, read back from the assets rather than
    // tracked alongside them -- a running counter would be one more thing that
    // can disagree with what is actually on screen.
    let mut resident = 0usize;
    let mut vertices = 0usize;
    let mut triangles = 0usize;
    for (_, handle) in &chunks {
        resident += 1;
        if let Some(mesh) = meshes.get(&handle.0) {
            vertices += mesh.count_vertices();
            triangles += mesh.indices().map_or(0, |i| i.len() / 3);
        }
    }
    stats.title = format!("E-202 game dig - {resident} chunks resident");
    stats.vertices = vertices;
    stats.triangles = triangles;
    stats.extract_ms = world.last_edit_ms;
    stats.extra = vec![
        format!(
            "edit log {:>4} brushes    brush radius {:.2}",
            world.brushes.len(),
            world.radius
        ),
        String::new(),
        format!(
            "last edit: {:>3} chunks re-meshed in {:.2} ms",
            world.last_chunks, world.last_edit_ms
        ),
        match world.last_edit {
            // `output_changed_cells`, not `value_changed_cells`. M-34: counting
            // cells whose *samples* moved reads 100% and says incremental
            // meshing is pointless; counting cells whose *triangles* move is
            // 15-36% and says the opposite. E1 is the second one.
            Some(r) => format!(
                "           {} of {} cells in the box re-mesh = E1 {:.1}%  ({} moved a sample)",
                r.output_changed_cells,
                r.region_cells,
                100.0 * r.changed_fraction(),
                r.value_changed_cells
            ),
            None => "           (click to carve)".to_string(),
        },
        match world.last_edit {
            Some(r) => format!(
                "           {} of {} chunks in the box were dirty = {:.1}%",
                r.dirty_chunks,
                r.region_chunks,
                100.0 * r.dirty_chunk_fraction()
            ),
            None => String::new(),
        },
        String::new(),
        format!(
            "mesher     last switch {} chunks in {:.0} ms across frames   \
             (rendering is always GPU: triplanar ExtendedMaterial)",
            world.switch_chunks, world.switch_ms
        ),
        match world.gpu_cost {
            // `base` first, because it is the term this demo exists to show
            // going away: it is the CPU field evaluation plus the upload, and
            // after a chunk's first mesh it is zero and says `(cached)`. What is
            // left is a fold on the device and three dispatches.
            Some(cost) => format!(
                "           base {:.2} ms{}   fold {:.2}   count {:.2}   scan {:.2}   emit {:.2}   peak {} tris/chunk of {}{}",
                cost.base_ms,
                if cost.base_ms > 0.0 { "" } else { " (cached)" },
                cost.fold_ms,
                cost.timings.count_ms,
                cost.timings.scan_ms,
                cost.timings.emit_ms,
                world.gpu_peak_triangles,
                GPU_TRIANGLE_BUDGET,
                if world.gpu_failures > 0 {
                    format!("   {} FAILED READBACKS", world.gpu_failures)
                } else {
                    String::new()
                }
            ),
            None => "           (press 8 to extract on the GPU instead)".to_string(),
        },
        // The verdict, because a demo that offers a GPU option and says nothing
        // invites the assumption that GPU means faster. The `0.06`/`0.22` pair is
        // `docs/measurements/gpu_vs_cpu.csv`'s 17^3 row -- the CPU column and the
        // device-field GPU column -- and a chunk here is exactly 17^3. It is a
        // citation rather than a live measurement, which is why it names M-296.
        format!(
            "           {} KB of samples uploaded since start; measured 17^3: CPU 0.06 ms, GPU 0.22 (M-296)",
            fields.sample_bytes / 1024
        ),
        format!(
            "backlog    {:>4} chunks queued, {:.1} ms/frame budget{}",
            world.backlog,
            MESH_BUDGET.as_secs_f64() * 1000.0,
            if world.switching {
                "   (switching)"
            } else {
                ""
            }
        ),
        format!(
            "stroke     {:>3} edits, {:.1} ms total, {:.2} ms/edit",
            world.stroke_edits,
            world.stroke_ms,
            if world.stroke_edits == 0 {
                0.0
            } else {
                world.stroke_ms / f64::from(world.stroke_edits)
            }
        ),
        format!(
            "body       {}   {}",
            if world.walking { "walking" } else { "flying " },
            if !world.walking {
                "[F] walk".to_string()
            } else if world.grounded {
                "on ground   [Space] jump   [F] fly".to_string()
            } else {
                format!("airborne {:+.1} u/s   [F] fly", world.velocity.y)
            }
        ),
        "ground     grass/dirt/deep-dirt array triplanar, 1.5 u/tile, sharpness 4; walls concrete"
            .to_string(),
        "touch      left third drags to walk, right drags to look, [JUMP] [DIG] [FILL] on screen"
            .to_string(),
        String::new(),
        "every field sample walks the log: measured 3.7x ms/chunk for 7x the log".to_string(),
    ];
    // The headline, above everything and in this mesher's own colour. It is here
    // rather than eight lines down the panel because "which mesher am I looking
    // at" is the question this demo exists to let a reader ask, and it survives
    // the panel being hidden -- which is how the demo opens.
    stats.banner = Some((
        format!("[{}] {}", world.algorithm.key(), world.algorithm.name()),
        world.algorithm.colour(),
    ));
    // What is left on screen with the panel off. The full key list, not just
    // `[H] HUD`: this demo starts hidden, so a reader who never presses `H`
    // would otherwise have no way to learn that `Space` jumps.
    stats.hint = Some(
        "[H] numbers   [LMB] carve   [RMB] fill   [WASD] walk   [Space] jump   [F] fly   [1-8] mesher   [Tab] cursor"
            .to_string(),
    );
    // The harness's shared footer advertises `[W] wire`, `[N] normals`,
    // `[G] grid` and `[R] re-mesh`. Every one of those is a lie here: `W` walks
    // forward, and the chunk entities carry `Chunk(ChunkId)` rather than
    // `DemoMesh`/`DemoDomain`, so the harness's wireframe, normal and domain
    // systems never see them. `[H] HUD` is repeated from the shared footer
    // because replacing it is what `DemoStats::keys` does, and a key that hides
    // the panel it is printed on must be findable on that panel.
    stats.keys = Some(
        "[LMB] hold to carve   [RMB] hold to fill   [WASD] move   [Space] jump   [F] walk/fly   [Shift] fast\n\
         [wheel] or [ ] brush   [1-8] mesher   [Z] undo   [X] clear log   [C] chunks   [H] HUD   [Tab] cursor   [F12] shot"
            .to_string(),
    );
}

/// Outline the chunks the last edit re-meshed.
///
/// This is the "chunks-touched count on screen" the ticket asks for, made
/// spatial: the count says how much work an edit cost, the boxes say *where*,
/// and only the second one shows you that a brush straddling a corner costs
/// eight chunks rather than one.
fn outline_chunks(world: Res<World>, mut gizmos: Gizmos<ChunkGizmos>) {
    if !world.show_chunks {
        return;
    }
    let span = world.layout.cell_size() * CHUNK_CELLS as f32;
    for id in &world.last_touched {
        let origin = world.layout.sample_origin(*id);
        let centre = Vec3::new(
            origin[0] + span * 0.5,
            origin[1] + span * 0.5,
            origin[2] + span * 0.5,
        );
        gizmos.cube(
            Transform::from_translation(centre).with_scale(Vec3::splat(span)),
            Color::srgb(0.20, 0.85, 1.0),
        );
    }
}

#[cfg(test)]
mod tests {
    use std::sync::{Arc, LazyLock};

    use bevy::asset::AssetApp;
    use bevy::ecs::system::RunSystemOnce;
    use bevy::mesh::VertexAttributeValues;
    use bevy::render::renderer::WgpuWrapper;
    use bevy::time::TimeUpdateStrategy;

    use super::*;

    /// A frame, fixed. `TimeUpdateStrategy::ManualDuration` is what makes
    /// `time.delta_secs()` this rather than however long the test machine took,
    /// which is what makes an edit count assertable.
    const FRAME: Duration = Duration::from_millis(16);

    /// **One** device and queue for this whole test binary, cloned per test.
    ///
    /// `isomesh_gpu::headless::Gpu` lends its device and `RenderDevice` needs to
    /// own one, so this opens a second. No software fallback, the same policy
    /// that module states: a test that quietly ran on a CPU reference driver
    /// would report a GPU comparison that is not one.
    ///
    /// **A `LazyLock`, and it is a fix rather than a tidy-up.** This used to
    /// build a fresh `wgpu::Instance` per call, so every `harness` test opened
    /// its own Vulkan instance and its own full-limits device. Under
    /// `cargo test`'s default parallelism -- 21 tests, 24 threads -- several
    /// `request_adapter`/`request_device` calls were in flight at once and the
    /// run **deadlocked**: reproduced at `--test-threads` 12 and 24 with the
    /// same five `harness` tests never reporting, and clean at 8 and 16, which is
    /// a race rather than a threshold. `wgpu::Device` and `wgpu::Queue` are
    /// `Arc`-backed handles, so one instance shared and cloned is both correct
    /// and what a real application does -- an engine does not open a device per
    /// system.
    ///
    /// The static holds the instance for the process's life, which is exactly the
    /// lifetime the device needs and no longer.
    fn wgpu_pair() -> (wgpu::Device, wgpu::Queue) {
        static SHARED: LazyLock<(wgpu::Device, wgpu::Queue)> = LazyLock::new(|| {
            let instance = wgpu::Instance::new(
                wgpu::InstanceDescriptor::new_without_display_handle_from_env(),
            );
            let adapter =
                bevy::tasks::block_on(instance.request_adapter(&wgpu::RequestAdapterOptions {
                    power_preference: wgpu::PowerPreference::HighPerformance,
                    force_fallback_adapter: false,
                    compatible_surface: None,
                }))
                .expect("a GPU adapter; there is no software fallback, by design");
            let limits = adapter.limits();
            bevy::tasks::block_on(adapter.request_device(&wgpu::DeviceDescriptor {
                label: Some("game_dig test"),
                required_features: wgpu::Features::empty(),
                required_limits: limits,
                memory_hints: wgpu::MemoryHints::Performance,
                experimental_features: wgpu::ExperimentalFeatures::disabled(),
                trace: wgpu::Trace::Off,
            }))
            .expect("a device on that adapter")
        });
        (*SHARED).clone()
    }

    /// The demo's own systems, in an `App` with no window and no renderer.
    ///
    /// This is the closest thing to running the demo that a machine with no
    /// display can do, and it exercises the wiring that the unit tests above
    /// cannot: the dirty set surviving the frame, the budget, the modal, the edit
    /// pacing and the keys. `RenderDevice` and `RenderQueue` are inserted by hand
    /// because `bevy_render` is what normally does it -- into the **main** world,
    /// which is exactly why `drain_dirty` can take them without a `RenderApp`.
    ///
    /// `grab`, `aim`, `ghost`, `report` and `outline_chunks` are left out on
    /// purpose: they want a window, gizmos or a `DemoStats`, and the tests below
    /// drive [`Aim`] directly so the aim point is a fixture rather than a
    /// consequence of where a camera happens to be looking.
    fn harness(walking: bool) -> App {
        let (device, queue) = wgpu_pair();
        let mc = MarchingCubesGpu::new(&device, &queue).expect("the compute pipeline");
        let fields = GpuFields {
            sampler: FieldSampler::new(&device).expect("the field-fold pipeline"),
            bases: HashMap::new(),
            log: Vec::new(),
            sample_bytes: 0,
        };
        let layout = test_layout();
        let mut world = World {
            layout,
            brushes: Vec::new(),
            survivors: Vec::new(),
            dirty: DirtySet::new(),
            backlog: 0,
            switching: true,
            radius: 0.25,
            last_touched: Vec::new(),
            last_edit: None,
            last_edit_ms: 0.0,
            last_chunks: 0,
            show_chunks: true,
            grabbed: true,
            algorithm: Algorithm::MarchingCubes,
            switch_ms: 0.0,
            switch_chunks: 0,
            switch_started: Instant::now(),
            stroke_last: None,
            stroke_edits: 0,
            stroke_ms: 0.0,
            stroke_clock: EDIT_PERIOD,
            velocity: Vec3::ZERO,
            walking,
            grounded: false,
            gpu_cost: None,
            gpu_peak_triangles: 0,
            gpu_failures: 0,
        };
        world.switch_chunks = mark_sandbox(&mut world.dirty);

        let mut app = App::new();
        app.add_plugins(MinimalPlugins)
            .add_plugins(bevy::asset::AssetPlugin::default())
            .init_asset::<Mesh>()
            .insert_resource(TimeUpdateStrategy::ManualDuration(FRAME))
            .init_resource::<ButtonInput<KeyCode>>()
            .init_resource::<ButtonInput<MouseButton>>()
            .init_resource::<AccumulatedMouseMotion>()
            .init_resource::<AccumulatedMouseScroll>()
            .init_resource::<Aim>()
            .init_resource::<GpuPending>()
            .init_resource::<TouchIntent>()
            .insert_resource(AutoCarve::default())
            .insert_resource(Capture::default())
            .insert_resource(Look {
                yaw: 0.0,
                pitch: 0.0,
            })
            .insert_resource(SurfaceMaterial(Handle::default()))
            .insert_resource(RenderDevice::from(device))
            .insert_resource(RenderQueue(Arc::new(WgpuWrapper::new(queue))))
            .insert_resource(GpuMesher { mc })
            .insert_resource(fields)
            .insert_resource(world)
            .add_systems(
                Update,
                (
                    switch_algorithm,
                    move_camera,
                    dig,
                    drain_dirty,
                    gpu_collect,
                    loading_modal,
                )
                    .chain(),
            );
        app.world_mut()
            .spawn((Camera3d::default(), Transform::from_xyz(0.0, 1.70, 6.0)));
        app.world_mut().spawn((Visibility::Hidden, LoadingModal));
        app
    }

    /// One frame, with the input clearing `InputPlugin` would have done.
    ///
    /// Without it `just_pressed` stays true for ever, and a jump would fire on
    /// every frame of the test rather than on the one the key was pressed.
    fn step(app: &mut App) {
        app.update();
        app.world_mut()
            .resource_mut::<ButtonInput<KeyCode>>()
            .clear();
        app.world_mut()
            .resource_mut::<ButtonInput<MouseButton>>()
            .clear();
    }

    fn eye(app: &mut App) -> Vec3 {
        let mut cameras = app
            .world_mut()
            .query_filtered::<&Transform, With<Camera3d>>();
        cameras
            .iter(app.world())
            .next()
            .expect("the camera")
            .translation
    }

    fn modal(app: &mut App) -> Visibility {
        let mut modals = app
            .world_mut()
            .query_filtered::<&Visibility, With<LoadingModal>>();
        *modals.iter(app.world()).next().expect("the modal")
    }

    /// The 256-chunk fill drains across frames, and the modal is up exactly while
    /// it does.
    ///
    /// The switch that used to block for most of a second is the reason both step
    /// 3 and step 4 exist, and this is the assertion that ties them together: on
    /// every frame with a backlog the overlay is visible, on the frame the queue
    /// empties it is not, and no chunk is lost on the way -- the entity count
    /// only grows.
    #[test]
    fn the_modal_is_up_exactly_while_a_backlog_drains() {
        let mut app = harness(false);
        let mut frames = 0;
        let mut resident = 0usize;
        loop {
            step(&mut app);
            frames += 1;
            let switching = app.world().resource::<World>().switching;
            let backlog = app.world().resource::<World>().backlog;
            let visibility = modal(&mut app);
            if switching {
                assert_eq!(
                    visibility,
                    Visibility::Inherited,
                    "frame {frames}: a backlog of {backlog} with no modal"
                );
            } else {
                assert_eq!(
                    visibility,
                    Visibility::Hidden,
                    "frame {frames}: the modal outlived the backlog"
                );
            }
            let now = app
                .world_mut()
                .query_filtered::<(), With<Chunk>>()
                .iter(app.world())
                .count();
            assert!(
                now >= resident,
                "frame {frames}: chunk entities went from {resident} to {now}"
            );
            resident = now;
            if !switching && backlog == 0 {
                break;
            }
            assert!(frames < 600, "the fill never drained");
        }
        assert!(
            frames >= 2,
            "256 chunks drained inside one frame, so the budget did nothing"
        );
        assert!(
            resident > 32,
            "only {resident} chunks got a mesh out of 256 queued"
        );
        let world = app.world().resource::<World>();
        assert!(
            world.switch_ms > 0.0,
            "the drain never reported what it cost"
        );
        println!("the startup fill drained in {frames} frames, {resident} chunks resident");
    }

    /// Key `8` meshes the same sandbox on the GPU, through the real systems.
    ///
    /// The unit test above compares one chunk's triangle count; this drives the
    /// whole path -- `switch_algorithm` marks, `drain_dirty` dispatches under the
    /// same budget, `gpu_collect` finishes whatever came back -- and asserts the
    /// end state a reader would see: **the same number of chunks on screen as key
    /// `1` produced**, no refused readbacks, and the timing line populated.
    ///
    /// # And then the count gate on GPU-014: the per-edit sample upload is zero
    ///
    /// Every chunk's [`Ground`] base is uploaded exactly once, so after the
    /// switch has drained `GpuFields::sample_bytes` is `EXTENT`'s 256 chunks x
    /// 17^3 samples x 4 bytes = **5,030,912**, and a carve that re-meshes chunks
    /// already meshed adds **nothing** to it: what crosses the bus per edit is
    /// the 64-byte-per-brush log.
    ///
    /// A byte-count equality rather than a timing ratio, deliberately. A
    /// wall-clock assertion on a machine with a CPU governor swings by tens of
    /// percent run to run and cannot gate anything (M-304's sibling lesson);
    /// bytes are exact, machine-independent, and name the mechanism directly.
    ///
    /// Needs a GPU adapter, as [`wgpu_pair`] says.
    #[test]
    fn key_eight_meshes_the_sandbox_on_the_gpu() {
        let mut app = harness(false);
        for _ in 0..600 {
            step(&mut app);
            if app.world().resource::<World>().backlog == 0 {
                break;
            }
        }
        let cpu_resident = app
            .world_mut()
            .query_filtered::<(), With<Chunk>>()
            .iter(app.world())
            .count();
        assert!(cpu_resident > 32, "the CPU fill produced {cpu_resident}");

        app.world_mut()
            .resource_mut::<ButtonInput<KeyCode>>()
            .press(KeyCode::Digit8);
        let mut frames = 0;
        loop {
            step(&mut app);
            frames += 1;
            let world = app.world().resource::<World>();
            let pending = app.world().resource::<GpuPending>().jobs.len();
            assert!(
                pending <= GPU_JOBS_MAX,
                "frame {frames}: {pending} jobs in flight against a cap of {GPU_JOBS_MAX}"
            );
            if !world.switching && world.backlog == 0 && pending == 0 {
                break;
            }
            assert!(frames < 900, "the GPU switch never drained");
        }
        let world = app.world().resource::<World>();
        assert_eq!(
            world.algorithm.name(),
            "marching_cubes (gpu)",
            "Digit8 did not select the GPU mesher"
        );
        assert_eq!(
            world.gpu_failures, 0,
            "{} readbacks were refused",
            world.gpu_failures
        );
        let cost = world.gpu_cost.expect("a GPU chunk was measured");
        assert!(
            cost.fold_ms > 0.0,
            "the device fold was not measured, so nothing was dispatched"
        );
        assert!(
            cost.timings.count_ms > 0.0 && cost.timings.emit_ms > 0.0,
            "the dispatches were not measured: {:?}",
            cost.timings
        );
        assert!(
            world.gpu_peak_triangles > 0 && world.gpu_peak_triangles < GPU_TRIANGLE_BUDGET,
            "the last GPU chunk reported {} triangles against a budget of {GPU_TRIANGLE_BUDGET}",
            world.gpu_peak_triangles
        );
        let gpu_resident = app
            .world_mut()
            .query_filtered::<(), With<Chunk>>()
            .iter(app.world())
            .count();
        assert_eq!(
            gpu_resident, cpu_resident,
            "the GPU left {gpu_resident} chunks on screen where the CPU left {cpu_resident}"
        );
        // Every dispatch in the switch drain was a chunk's first, so this split
        // is the cache-*miss* cost -- the one the cache exists to stop paying.
        println!(
            "cache-miss split: base {:.3} ms   fold {:.3}   count {:.3}   scan {:.3}   emit {:.3}",
            cost.base_ms,
            cost.fold_ms,
            cost.timings.count_ms,
            cost.timings.scan_ms,
            cost.timings.emit_ms
        );

        // The count gate. One base per chunk, uploaded once.
        let per_chunk = u64::from(CHUNK_CELLS + 1).pow(3) * 4;
        let sandbox_chunks = (EXTENT[0] * EXTENT[1] * EXTENT[2]) as u64;
        let after_switch = app.world().resource::<GpuFields>().sample_bytes;
        assert_eq!(
            after_switch,
            sandbox_chunks * per_chunk,
            "every chunk should have uploaded its base exactly once: \
             {sandbox_chunks} chunks x {per_chunk} bytes"
        );
        assert_eq!(
            app.world().resource::<GpuFields>().bases.len(),
            sandbox_chunks as usize,
            "the base cache holds a buffer per chunk and nothing else"
        );

        // Now carve, through the same path `a_held_button_edits_at_the_paced_rate`
        // drives, and drain the chunks it dirtied.
        //
        // **On the surface, not at a fixed height.** `Ground`'s top sits near
        // `y = -0.34` at this `x`/`z`, so a brush at `y = 0.3` subtracts air:
        // `mark_edit` correctly reports zero dirty chunks and the assertion below
        // would pass a test that re-meshed nothing. `Ground.sample([x, 0, z])` is
        // `-height`, so this places the brush centre exactly on the crossing.
        let before = app.world().resource::<World>().brushes.len();
        app.world_mut()
            .resource_mut::<ButtonInput<MouseButton>>()
            .press(MouseButton::Left);
        let mut carve_frames = 0;
        while app.world().resource::<World>().brushes.len() == before {
            let x = -1.0 + carve_frames as f32 * 0.2;
            let z = 1.0;
            let mut target = app.world_mut().resource_mut::<Aim>();
            target.hit = true;
            target.point = Vec3::new(x, -Ground.sample([x, 0.0, z]), z);
            drop(target);
            app.update();
            carve_frames += 1;
            assert!(carve_frames < 60, "the carve never landed");
        }
        app.world_mut()
            .resource_mut::<ButtonInput<MouseButton>>()
            .release(MouseButton::Left);
        for _ in 0..600 {
            step(&mut app);
            let world = app.world().resource::<World>();
            let pending = app.world().resource::<GpuPending>().jobs.len();
            if world.backlog == 0 && pending == 0 {
                break;
            }
        }
        let world = app.world().resource::<World>();
        assert!(
            world.last_chunks > 0,
            "the carve re-meshed no chunks, so it proved nothing about the upload"
        );
        assert_eq!(
            app.world().resource::<GpuFields>().sample_bytes,
            after_switch,
            "an edit uploaded samples: the base cache was missed on {} re-meshed chunks",
            world.last_chunks
        );
        println!(
            "{after_switch} bytes of samples uploaded for {sandbox_chunks} chunks, \
             and {} re-meshed chunks added none",
            world.last_chunks
        );
        // Reported, never asserted: this host's per-dispatch setup cost is not
        // the RTX 3090 Vulkan figure in `docs/measurements/gpu_vs_cpu.csv`, and a
        // bound invented from one machine's wall clock is a gate that cries
        // wolf. The split is what the HUD shows for a cache hit.
        let hit = world.gpu_cost.expect("the carve dispatched a chunk");
        println!(
            "cache-hit split: base {:.3} ms   fold {:.3}   count {:.3}   scan {:.3}   emit {:.3}",
            hit.base_ms,
            hit.fold_ms,
            hit.timings.count_ms,
            hit.timings.scan_ms,
            hit.timings.emit_ms
        );
        println!("the GPU switch drained in {frames} frames, {gpu_resident} chunks resident");
    }

    /// The HUD says what the new field path did, and marks the cached base.
    ///
    /// **This is the only way to see this demo's screen on a machine with no
    /// display**, and the two lines it checks are the deliverable of the HUD half
    /// of GPU-014: `base … (cached)` is how a reader watches the per-edit sample
    /// upload reach zero, and the verdict line is what stops a GPU key from
    /// implying that GPU means faster at this grid size. Both are strings a
    /// reader is meant to read, so a test that only checked `GpuCost`'s fields
    /// would pass with either of them missing.
    ///
    /// `report` is left out of [`harness`] because it wants a `DemoStats`; this
    /// runs it as a one-shot system with one inserted, which is the same system
    /// the demo runs every frame.
    #[test]
    fn the_hud_reports_the_cached_base_and_the_measured_verdict() {
        let mut app = harness(false);
        app.init_resource::<DemoStats>();
        app.world_mut()
            .resource_mut::<ButtonInput<KeyCode>>()
            .press(KeyCode::Digit8);
        for _ in 0..900 {
            step(&mut app);
            let world = app.world().resource::<World>();
            if !world.switching
                && world.backlog == 0
                && app.world().resource::<GpuPending>().jobs.is_empty()
            {
                break;
            }
        }
        // A second pass over an already-meshed chunk, so `base_ms` is a hit.
        let cached = app
            .world()
            .resource::<GpuFields>()
            .bases
            .keys()
            .copied()
            .next()
            .expect("a cached chunk");
        app.world_mut().resource_mut::<World>().dirty.insert(cached);
        step(&mut app);
        app.world_mut()
            .run_system_once(report)
            .expect("the HUD system");

        let lines = app.world().resource::<DemoStats>().extra.clone();
        for line in &lines {
            println!("{line}");
        }
        let gpu = lines
            .iter()
            .find(|l| l.contains("base "))
            .expect("the GPU timing line");
        assert!(
            gpu.contains("(cached)"),
            "a re-meshed chunk did not report a cached base: {gpu}"
        );
        assert!(
            gpu.contains("fold "),
            "the HUD dropped the device fold term: {gpu}"
        );
        let verdict = lines
            .iter()
            .find(|l| l.contains("M-296"))
            .expect("the measured-verdict line");
        assert!(
            verdict.contains("CPU 0.06 ms") && verdict.contains("GPU 0.22"),
            "the verdict line stopped citing the 17^3 row: {verdict}"
        );
        assert!(
            verdict.contains(&format!("{} KB", 256 * 4913 * 4 / 1024)),
            "the uploaded-sample total is not the 256-chunk figure: {verdict}"
        );
    }

    /// A held button edits at [`EDIT_PERIOD`], not once a frame.
    ///
    /// The arithmetic is exact and that is the point: the clock starts due, so an
    /// edit lands on frame 1 and then every **fifth** frame at a 16 ms frame and
    /// an 80 ms period -- frames 1, 6, 11, 16, 21, 26, 31, 36. **Eight edits in
    /// forty frames.** Before any of this the answer was forty, and every one of
    /// those brushes is a term in the `(L + 1)` factor on every field sample the
    /// demo takes afterwards.
    #[test]
    fn a_held_button_edits_at_the_paced_rate() {
        let mut app = harness(false);
        // Drain the startup fill first, so the frames counted below are the
        // stroke's own.
        for _ in 0..600 {
            step(&mut app);
            if app.world().resource::<World>().backlog == 0 {
                break;
            }
        }
        app.world_mut()
            .resource_mut::<ButtonInput<MouseButton>>()
            .press(MouseButton::Left);
        for frame in 0..40 {
            // A moving aim point, because the distance gate is still there and a
            // stationary brush is an idempotent duplicate it correctly refuses.
            let mut target = app.world_mut().resource_mut::<Aim>();
            target.hit = true;
            target.point = Vec3::new(-1.0 + frame as f32 * 0.2, 0.3, 1.0);
            drop(target);
            app.update();
            app.world_mut()
                .resource_mut::<ButtonInput<KeyCode>>()
                .clear();
        }
        let edits = app.world().resource::<World>().brushes.len();
        assert_eq!(
            edits, 8,
            "forty 16 ms frames at an 80 ms period is eight edits, not {edits}"
        );
        // And the modal never appeared: a dig's chunks clear inside one frame's
        // budget, so a reader carving does not get a "Loading" panel.
        assert_eq!(modal(&mut app), Visibility::Hidden);
    }

    /// `Space` jumps only from the ground, and `F` swaps the mode.
    ///
    /// The physics is asserted directly further up; this is the wiring -- that
    /// the keys reach it and that fly mode really stops gravity.
    #[test]
    fn the_jump_and_the_mode_keys_reach_the_body() {
        let mut app = harness(true);
        // Let the body land, and the fill drain with it.
        for _ in 0..600 {
            step(&mut app);
            if app.world().resource::<World>().backlog == 0
                && app.world().resource::<World>().grounded
            {
                break;
            }
        }
        assert!(
            app.world().resource::<World>().grounded,
            "the body never landed; eye {:?}",
            eye(&mut app)
        );
        let standing = eye(&mut app);

        app.world_mut()
            .resource_mut::<ButtonInput<KeyCode>>()
            .press(KeyCode::Space);
        step(&mut app);
        assert!(
            app.world().resource::<World>().velocity.y > 0.0,
            "Space did not launch the body"
        );
        let mut peak = standing.y;
        for _ in 0..40 {
            step(&mut app);
            peak = peak.max(eye(&mut app).y);
        }
        // `8.5^2 / (2 * 18)` is a 2.0-unit apex -- a whole chunk -- and a 16 ms
        // frame samples it coarsely, so this asserts most of the number rather
        // than all of it. At the old `JUMP_SPEED` of 6.0 the apex was 1.0 and
        // this would fail, which is what makes it a gate on the change.
        assert!(
            peak > standing.y + 1.6,
            "the jump reached {peak} from {}, so it is not a chunk high",
            standing.y
        );
        for _ in 0..120 {
            step(&mut app);
        }
        assert!(
            (eye(&mut app).y - standing.y).abs() < 0.05,
            "the body did not come back down to {}: {:?}",
            standing.y,
            eye(&mut app)
        );

        // `F` is fly mode: gravity stops and `E` climbs.
        app.world_mut()
            .resource_mut::<ButtonInput<KeyCode>>()
            .press(KeyCode::KeyF);
        step(&mut app);
        assert!(!app.world().resource::<World>().walking, "F did not toggle");
        let flying_from = eye(&mut app);
        app.world_mut()
            .resource_mut::<ButtonInput<KeyCode>>()
            .press(KeyCode::KeyE);
        for _ in 0..30 {
            app.update();
            app.world_mut()
                .resource_mut::<ButtonInput<KeyCode>>()
                .clear_just_pressed(KeyCode::KeyE);
        }
        assert!(
            eye(&mut app).y > flying_from.y + 0.5,
            "E did not climb in fly mode: {:?} from {flying_from:?}",
            eye(&mut app)
        );
        assert_eq!(
            app.world().resource::<World>().velocity,
            Vec3::ZERO,
            "fly mode is still integrating gravity"
        );
    }

    /// The layout the demo runs on, so a test measures the shipped geometry
    /// rather than a convenient one.
    fn test_layout() -> ChunkLayout<f32> {
        ChunkLayout::<f32>::new(
            CHUNK_CELLS,
            CELL_SIZE,
            [
                -(EXTENT[0] as f32) * CHUNK_CELLS as f32 * CELL_SIZE * 0.5,
                -5.4,
                -(EXTENT[2] as f32) * CHUNK_CELLS as f32 * CELL_SIZE * 0.5,
            ],
        )
        .expect("the layout the demo ships with")
    }

    /// An edit log spread across the sandbox, near the terrain surface.
    ///
    /// A lattice walk rather than an RNG: no dependency, and a failure is
    /// reproducible from the brush index alone. Spread matters -- the whole
    /// point of pruning is that most of the tape is nowhere near any one chunk,
    /// so a tape clustered around the origin would make the test vacuous.
    fn tape(n: usize) -> Vec<Brush<Sphere<f32>>> {
        (0..n)
            .map(|i| {
                let t = i as f32;
                let shape = Sphere {
                    center: [
                        -7.0 + (t * 0.937).rem_euclid(14.0),
                        0.2 * (t * 0.61).sin(),
                        -7.0 + (t * 2.113).rem_euclid(14.0),
                    ],
                    radius: 0.30 + 0.20 * (t * 0.37).cos().abs(),
                };
                // Both ops, because `Add` and `Subtract` take different arms of
                // the pruner and do not commute with each other.
                if i % 3 == 0 {
                    Brush::add(shape)
                } else {
                    Brush::subtract(shape)
                }
            })
            .collect()
    }

    /// Every float of a mesh, as bits, plus its indices.
    ///
    /// Bits, not values: `f32 == f32` calls `-0.0` equal to `0.0`, and `-0.0` is
    /// exactly what the pruner's one soft spot would produce.
    fn signature(mesh: Option<&Mesh>) -> Vec<u32> {
        let Some(mesh) = mesh else {
            return Vec::new();
        };
        let mut out = Vec::new();
        if let Some(VertexAttributeValues::Float32x3(values)) =
            mesh.attribute(Mesh::ATTRIBUTE_POSITION)
        {
            out.extend(values.iter().flatten().map(|f| f.to_bits()));
        }
        if let Some(VertexAttributeValues::Float32x3(values)) =
            mesh.attribute(Mesh::ATTRIBUTE_NORMAL)
        {
            out.extend(values.iter().flatten().map(|f| f.to_bits()));
        }
        if let Some(Indices::U32(indices)) = mesh.indices() {
            out.extend(indices.iter().copied());
        }
        out
    }

    /// Pruning the tape must not move a single bit of any chunk's mesh.
    ///
    /// This is the gate on `drain_dirty`'s speed-up. The pruner drops a brush
    /// only when its enclosure proves the brush cannot change the fold anywhere
    /// in the chunk's box, so "faster" has to mean *byte-identical*, on every
    /// extractor, or the demo's whole premise -- that the numbers on screen are
    /// the numbers the algorithm produced -- is gone.
    ///
    /// Both halves are asserted: the meshes agree, **and** the pruner actually
    /// dropped something. A pruner that returns the whole tape would pass the
    /// first half and buy nothing, which is the failure a timing measurement
    /// cannot distinguish from a slow machine.
    #[test]
    fn pruning_leaves_every_mesh_bit_identical() {
        let layout = test_layout();
        let log = tape(32);
        let span = layout.cell_size() * layout.cells() as f32;
        let mut survivors = Vec::new();
        let mut worst = 0usize;
        let mut best = usize::MAX;
        for algorithm in Algorithm::ALL {
            if algorithm == Algorithm::MarchingCubesGpu {
                // Its geometry comes back from the device; there is nothing to
                // extract on this thread. `gpu_dispatch` samples the same pruned
                // field, so the bound it relies on is the one asserted here.
                continue;
            }
            // A slab at the terrain surface: chunk y = 2 spans y in [-1.4, 0.6]
            // and the height field lives inside that.
            for x in 0..4 {
                for z in 0..2 {
                    let id = ChunkId::new([x, 2, z]);
                    let origin = layout.sample_origin(id);
                    let whole = algorithm.mesh(
                        &layout,
                        &BrushStack {
                            base: Ground,
                            brushes: &log,
                        },
                        origin,
                    );
                    let kept =
                        prune_into(&log, &Ground, ChunkBox::new(origin, span), &mut survivors);
                    worst = worst.max(kept);
                    best = best.min(kept);
                    let pruned = algorithm.mesh(
                        &layout,
                        &BrushStack {
                            base: Ground,
                            brushes: &survivors,
                        },
                        origin,
                    );
                    assert_eq!(
                        signature(whole.as_ref()),
                        signature(pruned.as_ref()),
                        "{} disagreed on chunk {:?} with {kept} of {} brushes kept",
                        algorithm.name(),
                        id.coords,
                        log.len()
                    );
                }
            }
        }
        assert!(
            worst < log.len(),
            "the pruner kept the whole {}-brush tape on some chunk, so it is not firing",
            log.len()
        );
        assert!(
            best <= log.len() / 2,
            "the best chunk kept {best} of {}, which is not the several-fold cut the mechanism is \
             claimed to be",
            log.len()
        );
    }

    /// Pruning both sides of an edit must not move E1 or the dirty set.
    ///
    /// `dig` prunes the pair it hands `mark_edit`, and `mark_edit`'s report *is*
    /// the measurement this example exists to publish. Identical means identical:
    /// every counter, and the same chunks in the same order.
    #[test]
    fn pruning_leaves_the_edit_report_identical() {
        let layout = test_layout();
        let log = tape(32);
        let split = log.len() - 1;
        let pushed = log[split].shape;
        let reach = pushed.radius + layout.cell_size();
        let min_cell = layout.cell_of([
            pushed.center[0] - reach,
            pushed.center[1] - reach,
            pushed.center[2] - reach,
        ]);
        let max_cell = layout.cell_of([
            pushed.center[0] + reach,
            pushed.center[1] + reach,
            pushed.center[2] + reach,
        ]);

        let mut plain = DirtySet::new();
        let plain_report = mark_edit(
            &layout,
            &BrushStack {
                base: Ground,
                brushes: &log[..split],
            },
            &BrushStack {
                base: Ground,
                brushes: &log,
            },
            min_cell,
            max_cell,
            &mut plain,
        )
        .expect("a brush of this radius spans a few cells");
        assert!(
            plain_report.output_changed_cells > 0,
            "the fixture's last brush changed nothing, so this test proves nothing"
        );

        // The same box `dig` builds: the widest axis, so a non-cubic union is
        // still enclosed.
        let widest = (0..3)
            .map(|axis| max_cell[axis] - min_cell[axis] + 1)
            .max()
            .expect("three axes");
        let edit_box = ChunkBox::new(
            layout.world_of_sample(min_cell),
            widest as f32 * layout.cell_size(),
        );
        let mut before_kept = Vec::new();
        let mut after_kept = Vec::new();
        prune_into(&log[..split], &Ground, edit_box, &mut before_kept);
        let kept = prune_into(&log, &Ground, edit_box, &mut after_kept);
        assert!(
            kept < log.len(),
            "the pruner kept the whole tape over the edit box, so it is not firing"
        );

        let mut pruned = DirtySet::new();
        let pruned_report = mark_edit(
            &layout,
            &BrushStack {
                base: Ground,
                brushes: &before_kept,
            },
            &BrushStack {
                base: Ground,
                brushes: &after_kept,
            },
            min_cell,
            max_cell,
            &mut pruned,
        )
        .expect("the same region, one brush shorter");
        assert_eq!(plain_report, pruned_report, "E1 moved under pruning");
        assert_eq!(
            plain.iter().collect::<Vec<_>>(),
            pruned.iter().collect::<Vec<_>>(),
            "the dirty set moved under pruning"
        );
    }

    /// A halo'd chunk must mesh into the borrowed cell layer, and a chunk without
    /// one must not.
    ///
    /// This is the gate on the gaps between areas that keys `4`, `5` and `6` used
    /// to show. The three dual methods emit one quad per crossed grid edge and
    /// `dual.rs`'s walk skips the outermost quad plane on every face, so a chunk
    /// given only its own cells stops one cell short of its positive faces and
    /// neither neighbour bridges it. With one borrowed layer the last quad plane
    /// lands on the shared plane instead, and the tell is geometric and exact:
    /// **some vertex must sit beyond it**, in the cell the neighbour owns.
    ///
    /// The edge-based families must stay on the near side of that plane, because
    /// they already tile and a borrowed layer would make them mesh the
    /// neighbour's first cell a second time.
    #[test]
    fn the_halo_reaches_the_borrowed_layer_and_nothing_else_does() {
        let layout = test_layout();
        // Bare terrain: the height field crosses every column of this chunk, so
        // no brush is needed to put a surface in the borrowed layer.
        let empty: [Brush<Sphere<f32>>; 0] = [];
        let field = BrushStack {
            base: Ground,
            brushes: &empty,
        };
        let id = ChunkId::new([3, 2, 3]);
        let origin = layout.sample_origin(id);
        // The plane this chunk shares with its positive-x neighbour.
        let shared = origin[0] + layout.cells() as f32 * layout.cell_size();
        for algorithm in Algorithm::ALL {
            if algorithm == Algorithm::MarchingCubesGpu {
                continue;
            }
            let mesh = algorithm
                .mesh(&layout, &field, origin)
                .expect("the terrain crosses this chunk");
            let Some(VertexAttributeValues::Float32x3(positions)) =
                mesh.attribute(Mesh::ATTRIBUTE_POSITION)
            else {
                panic!("{} emitted no positions", algorithm.name());
            };
            let max_x = positions
                .iter()
                .map(|p| p[0])
                .fold(f32::NEG_INFINITY, f32::max);
            let min_x = positions.iter().map(|p| p[0]).fold(f32::INFINITY, f32::min);
            assert!(
                min_x >= origin[0] - 1e-6,
                "{} meshed below its own origin",
                algorithm.name()
            );
            if algorithm.halo_cells() > 0 {
                assert!(
                    max_x > shared,
                    "{} stopped at {max_x} and never reached the borrowed layer past {shared}, so \
                     the seam is still open",
                    algorithm.name()
                );
            } else {
                assert!(
                    max_x <= shared + 1e-6,
                    "{} meshed past {shared} to {max_x}, so it is duplicating its neighbour's \
                     first cell",
                    algorithm.name()
                );
            }
        }
    }

    /// A budgeted drain must keep every chunk it could not reach.
    ///
    /// The predicate is a **counter**, not a clock: `MESH_BUDGET` on a fast
    /// machine and a slow one meshes a different number of chunks, and a test
    /// that asserted a frame count would be measuring the machine. What
    /// `drain_dirty` actually depends on is the contract -- the chunks a call
    /// could not reach stay in the set -- and stopping after exactly `K` makes
    /// that exact: the union over calls is the original set, once each. A
    /// dropped chunk is a permanent hole in the terrain, which is the failure
    /// this rules out.
    #[test]
    fn a_budgeted_drain_keeps_what_it_could_not_reach() {
        let layout = test_layout();
        let mut dirty = DirtySet::new();
        let total = mark_sandbox(&mut dirty);
        assert_eq!(total, 256, "the sandbox is 8x4x8 chunks");

        let mut seen: Vec<ChunkId> = Vec::new();
        let mut calls = 0;
        while !dirty.is_empty() {
            let mut this_call = 0;
            let report = dirty.mesh_within_budget(
                &layout,
                [0.0, 1.6, 6.0],
                |id, _| seen.push(id),
                || {
                    this_call += 1;
                    this_call < 7
                },
            );
            assert_eq!(report.meshed, 7.min(seen.len()).min(report.meshed));
            assert_eq!(
                report.remaining,
                dirty.iter().count(),
                "the report and the set disagree about what is left"
            );
            calls += 1;
            assert!(calls <= total, "the drain is not making progress");
        }
        assert_eq!(calls, total.div_ceil(7), "seven chunks a call, 256 chunks");
        seen.sort_unstable();
        let mut expected = DirtySet::new();
        mark_sandbox(&mut expected);
        assert_eq!(
            seen,
            expected.iter().collect::<Vec<_>>(),
            "the union of the calls is not the set that was queued, once each"
        );
    }

    /// The body falls, lands on the terrain, and then stands still.
    ///
    /// `resolve_body` and [`gravity_step`] are the whole of ask 6's physics that
    /// is not a keypress, so this drives them directly, in the order
    /// `move_camera` does: gravity, integrate, resolve. The eye must settle at
    /// `1.70` above the surface -- the height the four overlapping spheres put
    /// it, which is exactly where `setup` places the camera.
    ///
    /// **And then not move.** The drift assertion is the one that matters: the
    /// push is along the surface normal, which on a slope has a horizontal
    /// component, so a body that keeps being pushed out keeps sliding downhill.
    /// A slope is therefore in the fixture rather than only flat ground.
    #[test]
    fn a_dropped_body_lands_on_the_terrain_and_stands_still() {
        let empty: [Brush<Sphere<f32>>; 0] = [];
        let field = BrushStack {
            base: Ground,
            brushes: &empty,
        };
        let dt = 1.0 / 60.0;
        let lowest = BODY_OFFSETS.into_iter().fold(f32::MIN, f32::max);
        // Two columns: one where the height field is flat and one on its
        // steepest slope, so this is not a test of a single lucky point.
        for (x, z) in [(0.0f32, 0.0f32), (1.3, -2.1)] {
            let mut eye = Vec3::new(x, 4.0, z);
            let mut velocity = Vec3::ZERO;
            let mut grounded = false;
            for _ in 0..180 {
                gravity_step(grounded, &mut velocity, dt);
                eye.y += velocity.y * dt;
                grounded = resolve_body(&field, &mut eye, &mut velocity);
            }
            // `Ground` is `y - height`, so a sample at `y = 0` is `-height`.
            let height = -field.sample([eye.x, 0.0, eye.z]);
            let want = height + lowest + BODY_RADIUS;
            assert!(
                grounded,
                "the body never registered ground at ({x}, {z}); eye {eye:?}"
            );
            assert!(
                (eye.y - want).abs() < 0.05,
                "the body settled at {} instead of {want} at ({x}, {z})",
                eye.y
            );
            // Another three seconds of standing. A sliding body covers metres in
            // that time; a standing one must not move at all.
            let settled = eye;
            for _ in 0..180 {
                gravity_step(grounded, &mut velocity, dt);
                eye.y += velocity.y * dt;
                grounded = resolve_body(&field, &mut eye, &mut velocity);
            }
            assert!(
                eye.distance(settled) < 1e-3,
                "the body drifted {} units while standing at ({x}, {z})",
                eye.distance(settled)
            );
            assert_eq!(
                velocity.y, 0.0,
                "a standing body is carrying vertical velocity at ({x}, {z})"
            );
        }
    }

    /// Rock stops the body, and a carved pit lets it fall.
    ///
    /// Two halves of one claim, because "cannot walk through rock" is only
    /// interesting beside "and the hole you dug is a hole".
    ///
    /// The wall is a **large-radius** sphere, and the radius is the whole
    /// fixture: `resolve_body` has no slope limit by design, so a small sphere is
    /// a dome the body walks straight over -- which it did, on the first version
    /// of this test. At radius 40 the face is vertical to within
    /// `1.70² / 80 = 36 mm` over the body's height, which is a wall.
    #[test]
    fn rock_stops_the_body_and_a_hole_does_not() {
        let solid = [Brush::add(Sphere {
            center: [0.0, 0.4, -40.8],
            radius: 40.0,
        })];
        let field = BrushStack {
            base: Ground,
            brushes: &solid,
        };
        let mut eye = Vec3::new(0.0, 1.70, 0.0);
        let mut velocity = Vec3::ZERO;
        let dt = 1.0 / 60.0;
        let mut grounded = false;
        // Walking at the demo's own speed straight at the wall for a second and a
        // half, which is 2.5 units -- well past it if nothing stops the body.
        for _ in 0..90 {
            gravity_step(grounded, &mut velocity, dt);
            eye += Vec3::new(0.0, velocity.y, -2.5) * dt;
            grounded = resolve_body(&field, &mut eye, &mut velocity);
        }
        // The wall's face is at `z = -0.8` and every one of the body's spheres
        // has radius `BODY_RADIUS`, so no centre -- and the eye shares their `z`
        // -- can pass `-0.8 + BODY_RADIUS` by more than the field's own slack.
        // Derived rather than written as a literal: the body got narrower, and a
        // hardcoded bound would have had to be loosened by hand, which is how a
        // geometry test stops testing geometry.
        let closest = -0.8 + BODY_RADIUS;
        assert!(
            eye.z > closest - 0.05,
            "the body walked into the rock: z {} (a free walk reaches -2.25; \
             the wall allows {closest})",
            eye.z
        );

        let carved = [Brush::subtract(Sphere {
            center: [0.0, 0.0, -2.0],
            radius: 1.2,
        })];
        let pit = BrushStack {
            base: Ground,
            brushes: &carved,
        };
        let mut eye = Vec3::new(0.0, 1.70, -2.0);
        let mut velocity = Vec3::ZERO;
        let mut grounded = false;
        for _ in 0..240 {
            gravity_step(grounded, &mut velocity, dt);
            eye.y += velocity.y * dt;
            grounded = resolve_body(&pit, &mut eye, &mut velocity);
        }
        // The pit's floor is the bottom of the subtracted sphere, y = -1.2, so
        // the eye lands well below the 1.70 it stands at on the surface.
        assert!(
            eye.y < 1.0,
            "the body did not fall into the carved pit: eye y {}",
            eye.y
        );
        assert!(grounded, "the body did not land in the pit: eye {eye:?}");
    }

    /// On the medial axis the gradient is zero, and the body still gets out.
    ///
    /// **M-172** measured `BrushStack::gradient` returning exactly `[0, 0, 0]`
    /// there. A resolver that normalised that would produce `NaN` and put the
    /// camera nowhere; `resolve_body` takes `Vec3::Y` instead, because inside
    /// rock, up is the way out.
    #[test]
    fn the_medial_axis_pushes_up_rather_than_to_nowhere() {
        // A sphere added around the origin, and a body exactly at its centre --
        // the one point where the distance function's gradient vanishes.
        let solid = [Brush::add(Sphere {
            center: [0.0, 0.0, 0.0],
            radius: 1.0,
        })];
        let field = BrushStack {
            base: Ground,
            brushes: &solid,
        };
        let centre = Vec3::new(0.0, 0.0, 0.0);
        let gradient = Vec3::from_array(field.gradient([centre.x, centre.y, centre.z]));
        assert!(
            gradient.length() <= GRADIENT_EPS,
            "the fixture is not on a medial axis: gradient {gradient:?}"
        );
        let mut eye = centre + Vec3::Y * BODY_OFFSETS[0];
        let mut velocity = Vec3::ZERO;
        resolve_body(&field, &mut eye, &mut velocity);
        assert!(
            eye.is_finite(),
            "a zero gradient produced a non-finite position: {eye:?}"
        );
        assert!(
            eye.y > centre.y + BODY_OFFSETS[0],
            "the body was not pushed up out of the rock: eye {eye:?}"
        );
    }

    /// The GPU mesher must produce the same surface as key `1`.
    ///
    /// Marching Cubes is Marching Cubes, and after GPU-014 this is the stronger
    /// claim of the two: the CPU walks `BrushStack` per sample while the GPU
    /// folds the *same* log over an uploaded [`Ground`] base in a compute
    /// shader, so the triangle count agreeing means the two field paths agree
    /// too. `isomesh-gpu`'s own `a_fold_over_an_uploaded_base_matches_the_cpu`
    /// measures that fold sample-for-sample; this one asks whether it moved a
    /// crossing. The check a screenshot of key `8` can only suggest.
    ///
    /// **Needs a GPU adapter**, and fails rather than skips without one -- the
    /// same policy `isomesh-gpu`'s own tests hold, and for the reason its
    /// `headless` module states: a silent software fallback reports numbers three
    /// orders of magnitude off and looks exactly like a slow GPU.
    #[test]
    fn the_gpu_mesher_agrees_with_the_cpu_one() {
        let gpu = isomesh_gpu::headless::Gpu::new()
            .expect("a GPU adapter; there is no software fallback, by design");
        let mc = MarchingCubesGpu::new(gpu.device(), gpu.queue()).expect("the compute pipeline");
        let mut fields = GpuFields {
            sampler: FieldSampler::new(gpu.device()).expect("the field-fold pipeline"),
            bases: HashMap::new(),
            log: Vec::new(),
            sample_bytes: 0,
        };
        let layout = test_layout();
        let log = tape(8);
        let field = BrushStack {
            base: Ground,
            brushes: &log,
        };
        let mut compared = 0;
        for x in 0..4 {
            let id = ChunkId::new([x, 2, 1]);
            let origin = layout.sample_origin(id);
            let cpu = Algorithm::MarchingCubes
                .mesh(&layout, &field, origin)
                .map_or(0, |mesh| mesh.indices().map_or(0, |i| i.len() / 3));
            let (job, _) = gpu_dispatch(
                &mc,
                &mut fields,
                gpu.device(),
                gpu.queue(),
                &layout,
                &log,
                id,
                origin,
            )
            .expect("the dispatch");
            // The same poll `gpu_collect` does, with a bound so a device that
            // never completes fails the test rather than hanging the suite.
            let mut spins = 0;
            while !job.readback.ready(gpu.device()) {
                spins += 1;
                assert!(spins < 100_000, "the readback never completed");
            }
            let parts = job.readback.take().expect("the bytes");
            let (total, mesh) = gpu_mesh(&parts).expect("three buffers");
            assert!(
                total < GPU_TRIANGLE_BUDGET,
                "chunk {:?} hit the budget at {total} triangles, so it was truncated",
                id.coords
            );
            assert_eq!(
                mesh.is_some(),
                total > 0,
                "a mesh and a zero count disagree on chunk {:?}",
                id.coords
            );
            let slack = (cpu as f64 * 0.01).max(1.0);
            assert!(
                (f64::from(total) - cpu as f64).abs() <= slack,
                "chunk {:?}: GPU {total} triangles against CPU {cpu} on {}",
                id.coords,
                gpu.report().name
            );
            if cpu > 0 {
                compared += 1;
            }
        }
        assert!(
            compared >= 2,
            "only {compared} chunks had a surface, so this proved almost nothing"
        );
    }

    /// The committed PNGs really are four square layers, stacked.
    ///
    /// **The one silent failure in the pack.** The layer index in
    /// `triplanar.wgsl` is an array slice index, so a stack of three or five
    /// layers, a non-square layer, or a stack assembled bottom-up all compile,
    /// render, and put the wrong material on the surface -- bottom-up puts grass
    /// on the concrete walls. `reinterpret_stacked_2d_as_array` can only catch a
    /// height that does not divide, so the shape is asserted here, on the same
    /// bytes `setup` decodes.
    ///
    /// No GPU: this is a PNG decode and an extent.
    #[test]
    fn the_terrain_array_is_four_square_layers_stacked_top_down() {
        for (name, bytes, is_srgb) in [
            ("albedo/roughness", TERRAIN_ALBEDO_ROUGHNESS, true),
            ("normal/AO", TERRAIN_NORMAL_AO, false),
        ] {
            let flat = embedded_texture(bytes, is_srgb);
            assert_eq!(flat.width(), 512, "{name}: layer edge");
            assert_eq!(
                flat.height(),
                TERRAIN_LAYERS * flat.width(),
                "{name}: not {TERRAIN_LAYERS} square layers stacked"
            );
            assert_eq!(
                flat.texture_descriptor.size.depth_or_array_layers, 1,
                "{name}: the file is meant to be a flat stack before reinterpretation"
            );

            let array = terrain_array(bytes, is_srgb);
            assert_eq!(
                array.texture_descriptor.size.depth_or_array_layers, TERRAIN_LAYERS,
                "{name}: reinterpretation did not produce {TERRAIN_LAYERS} layers"
            );
            assert_eq!(array.width(), 512, "{name}: layer width after reinterpret");
            assert_eq!(
                array.height(),
                512,
                "{name}: layer height after reinterpret"
            );
            // The shader indexes layer 3 for the walls, so the top index has to
            // exist. Stated as the constant rather than as `3` so the two cannot
            // drift.
            assert!(
                LAYER_CONCRETE >= 0.0 && (LAYER_CONCRETE as u32) < TERRAIN_LAYERS,
                "the wall layer {LAYER_CONCRETE} is outside the array"
            );
            assert!(
                LAYER_BLEND < 0.0,
                "the terrain's sentinel {LAYER_BLEND} must be negative or it selects a layer"
            );
        }
    }

    /// The WGSL and this file agree about the array, or nothing renders right.
    ///
    /// The layer indices and the two bindings are stated in both places and
    /// **cannot be checked by running the demo on this host** -- there is no
    /// display here, and both failure modes are quiet where it counts:
    ///
    /// * a `texture_2d<f32>` in the WGSL against `dimension = "2d_array"` in the
    ///   `AsBindGroup` derive is a *pipeline creation* failure at the first draw,
    ///   a long way from either declaration;
    /// * a layer index that disagrees -- concrete listed as a terrain layer, or
    ///   `LAYER_CONCRETE` pointing at the dirt -- compiles, renders, and paints
    ///   the wrong material on the walls.
    ///
    /// So the shader source is read and cross-checked. `include_str!` of the same
    /// path `load_internal_asset!` compiles in, so there is one file and no copy.
    #[test]
    fn the_shader_and_this_file_agree_about_the_texture_array() {
        const WGSL: &str = include_str!("triplanar.wgsl");
        for binding in [101, 103] {
            let line = WGSL
                .lines()
                .find(|l| l.contains(&format!("@binding({binding})")))
                .unwrap_or_else(|| panic!("binding {binding} is not declared in triplanar.wgsl"));
            assert!(
                line.contains("texture_2d_array<f32>"),
                "binding {binding} is `{}`, which cannot match `dimension = \"2d_array\"`",
                line.trim()
            );
        }
        // The three the fragment shader blends, in the order `PROVENANCE.md`
        // stacked them.
        for (name, want) in [
            ("LAYER_GRASS", 0),
            ("LAYER_DIRT_SURFACE", 1),
            ("LAYER_DIRT_DEEP", 2),
        ] {
            let decl = format!("const {name}: i32 = {want};");
            assert!(
                WGSL.contains(&decl),
                "triplanar.wgsl does not declare `{decl}`, so the blend is sampling \
                 a layer this file did not pack there"
            );
            assert_ne!(
                want, LAYER_CONCRETE as i32,
                "{name} and LAYER_CONCRETE are the same slice, so the walls and the \
                 terrain share a material"
            );
        }
        // Every index in range, and the sentinel out of it.
        assert!(
            LAYER_CONCRETE >= 0.0 && (LAYER_CONCRETE as u32) < TERRAIN_LAYERS,
            "the wall layer {LAYER_CONCRETE} is outside a {TERRAIN_LAYERS}-layer array"
        );
        assert!(
            LAYER_BLEND < 0.0,
            "the terrain's sentinel {LAYER_BLEND} must be negative, or `settings.z` \
             selects it as a layer and the blend never runs"
        );
        // `settings.z` is the switch, so the shader has to branch on it. A shader
        // that ignored it would render the blend on the walls and look almost
        // right.
        assert!(
            WGSL.contains("triplanar.settings.z"),
            "triplanar.wgsl never reads `settings.z`, so the forced wall layer does \
             nothing"
        );
    }

    /// The body is a 1.70 x 0.50 capsule, and it hangs above nothing.
    ///
    /// Three claims. The **shape** is an identity on the constants, not an
    /// outcome of a fall, so it is asserted as one: `1.70` tall, `0.50` wide, and
    /// consecutive spheres overlapping rather than touching at a point, which is
    /// the difference between a capsule and the old pinched pair.
    ///
    /// The **settle** is bounded by [`GROUND_PROBE`] rather than by a hand-picked
    /// epsilon: `gravity_step` cuts gravity the moment the probe answers, so a
    /// body may come to rest up to one probe above the terrain and no further.
    /// That bound is what the *flat* footprint violated -- it read solid on
    /// terrain beside the foot that was higher than the terrain under it, and
    /// left the body hanging **12.5 cm** up at the origin.
    #[test]
    fn the_body_is_a_capsule_that_never_hangs_above_the_terrain() {
        let lowest = BODY_OFFSETS.into_iter().fold(f32::MIN, f32::max);
        assert!(
            ((lowest + BODY_RADIUS) - 1.70).abs() < 1e-6,
            "the body is {} tall, not 1.70",
            lowest + BODY_RADIUS
        );
        assert!(
            (2.0 * BODY_RADIUS - 0.50).abs() < 1e-6,
            "the body is {} wide, not 0.50",
            2.0 * BODY_RADIUS
        );
        for pair in BODY_OFFSETS.windows(2) {
            assert!(
                pair[1] - pair[0] < 2.0 * BODY_RADIUS,
                "spheres {} and {} are {} apart against a {} diameter, so the body \
                 has a pinched waist a lip of rock can pass through",
                pair[0],
                pair[1],
                pair[1] - pair[0],
                2.0 * BODY_RADIUS
            );
        }

        let empty: [Brush<Sphere<f32>>; 0] = [];
        let bare = BrushStack {
            base: Ground,
            brushes: &empty,
        };
        let dt = 1.0 / 60.0;
        // Five columns across the height field, including the origin where its
        // gradient is steepest (`0.63` per unit in `x`) -- that is the column the
        // flat footprint hovered over.
        for (x, z) in [
            (0.0f32, 0.0f32),
            (1.3, -2.1),
            (0.6, 0.9),
            (-2.4, 1.1),
            (3.1, -0.4),
        ] {
            let mut eye = Vec3::new(x, 4.0, z);
            let mut velocity = Vec3::ZERO;
            let mut grounded = false;
            for _ in 0..300 {
                gravity_step(grounded, &mut velocity, dt);
                eye.y += velocity.y * dt;
                grounded = resolve_body(&bare, &mut eye, &mut velocity);
            }
            assert!(grounded, "the body never landed at ({x}, {z}): {eye:?}");
            // `Ground` is `y - height`, so a sample at `y = 0` is `-height`. Taken
            // at the *settled* column: the push is along the surface normal, which
            // on a slope has a horizontal component, so the body does not land on
            // the column it was dropped down.
            let want = -bare.sample([eye.x, 0.0, eye.z]) + lowest + BODY_RADIUS;
            assert!(
                eye.y - want <= GROUND_PROBE,
                "the body is hanging {} above the terrain at ({x}, {z}), which is \
                 more than the {GROUND_PROBE} probe that stopped its fall",
                eye.y - want
            );
        }
    }

    /// A pit dug with the demo's own brush is not a trap.
    ///
    /// **The bug this fixes.** A body perched on the rim of a hole it dug has
    /// rock under the *outside* of its foot sphere and open air under the
    /// centre. `resolve_body` probed the centre and nothing else, so `grounded`
    /// went false, `move_camera`'s `world.grounded && Space` refused the jump,
    /// and the hole was inescapable.
    ///
    /// The fixture is asserted to *be* that case rather than assumed to be: the
    /// point the old probe sampled is checked to be air while the body is
    /// standing. Without that first assertion this test would pass on a body
    /// resting on flat ground and prove nothing.
    #[test]
    fn a_body_on_the_rim_of_its_own_pit_is_still_grounded() {
        let lowest = BODY_OFFSETS.into_iter().fold(f32::MIN, f32::max);
        let dt = 1.0 / 60.0;
        // A bowl half a unit across at the surface -- two passes of the default
        // brush -- and the body walking onto its edge rather than its middle.
        let carved = [Brush::subtract(Sphere {
            center: [0.0, 0.0, 0.0],
            radius: 0.5,
        })];
        let pit = BrushStack {
            base: Ground,
            brushes: &carved,
        };
        let mut eye = Vec3::new(0.3, 3.0, 0.0);
        let mut velocity = Vec3::ZERO;
        let mut grounded = false;
        for _ in 0..400 {
            gravity_step(grounded, &mut velocity, dt);
            eye.y += velocity.y * dt;
            grounded = resolve_body(&pit, &mut eye, &mut velocity);
        }
        let foot = eye - Vec3::Y * lowest;
        assert!(
            pit.sample([foot.x, foot.y - BODY_RADIUS - GROUND_PROBE, foot.z]) > 0.0,
            "the fixture is not the interesting case: there is rock directly under \
             the foot's centre, so one sample would already have found it"
        );
        assert!(
            grounded,
            "the body on the rim of its own pit read as falling, so the jump is \
             refused and the hole is a trap: {eye:?}"
        );
    }

    /// The five slabs line the sandbox, outside it, with the corners closed.
    ///
    /// Geometry rather than an entity count, because the count is the part that
    /// cannot be wrong by a little: a slab one thickness off still spawns five
    /// entities and leaves a gap of sky at every corner. Each assertion is a
    /// property of the box against [`sandbox`], so an `EXTENT` change moves both
    /// sides of it.
    #[test]
    fn five_slabs_line_the_sandbox_from_outside_with_the_corners_closed() {
        let layout = test_layout();
        let (lo, hi) = sandbox(&layout);
        let boxes = walls(&layout);
        assert_eq!(boxes.len(), 5, "four sides and a floor, and no ceiling");
        for (centre, size) in boxes {
            let slab_lo = centre - size * 0.5;
            let slab_hi = centre + size * 0.5;
            // Exactly one axis is a thickness, and it is the axis the slab
            // faces along.
            let thin = [size.x, size.y, size.z]
                .into_iter()
                .filter(|s| *s == WALL_THICKNESS)
                .count();
            assert_eq!(thin, 1, "a slab must be thin in exactly one axis: {size:?}");
            // Outside, not overlapping: the slab and the sandbox interior are
            // disjoint, so nothing z-fights the terrain and `aim` -- which
            // refuses to carve outside `sandbox` -- can never dig it.
            let disjoint = slab_hi.x <= lo.x
                || slab_lo.x >= hi.x
                || slab_hi.y <= lo.y
                || slab_lo.y >= hi.y
                || slab_hi.z <= lo.z
                || slab_lo.z >= hi.z;
            assert!(
                disjoint,
                "the slab at {centre:?} of {size:?} reaches inside the sandbox"
            );
            // Flush: the inner face is the boundary plane, so there is no gap
            // between the terrain's edge and the wall the player is stopped by.
            let flush = slab_hi.x == lo.x
                || slab_lo.x == hi.x
                || slab_hi.y == lo.y
                || slab_lo.y == hi.y
                || slab_hi.z == lo.z
                || slab_lo.z == hi.z;
            assert!(
                flush,
                "the slab at {centre:?} of {size:?} does not touch the boundary"
            );
        }
        // The corners close: every vertical slab spans the full height, and the
        // two facing each axis overhang far enough in the other to meet their
        // neighbours' outer faces.
        for (centre, size) in &boxes[..4] {
            assert_eq!(
                (centre.y - size.y * 0.5, centre.y + size.y * 0.5),
                (lo.y, hi.y),
                "a side wall does not span the sandbox's full height"
            );
            let (span, half) = if size.x == WALL_THICKNESS {
                (size.z, hi.z - lo.z)
            } else {
                (size.x, hi.x - lo.x)
            };
            assert_eq!(
                span,
                half + 2.0 * WALL_THICKNESS,
                "a side wall does not overhang, so its corners show sky"
            );
        }
        let (floor_centre, floor_size) = boxes[4];
        assert_eq!(
            floor_centre.y + floor_size.y * 0.5,
            lo.y,
            "the floor's top is not the sandbox's bottom"
        );
    }

    /// `setup` really spawns the five slabs, and paints them with **concrete**.
    ///
    /// [`walls`] proves the boxes are right; this proves they reach the screen
    /// wearing the other material. Both handles in `setup` are a
    /// `Handle<TerrainMaterial>` and their names differ by one word, so handing
    /// the walls `material` instead of `wall_material` compiles, spawns five
    /// entities, and renders the grass-and-dirt blend on the boundary. Nothing
    /// but reading the asset back catches that, and there is no display on this
    /// host to look at it with.
    ///
    /// `setup` runs in full here rather than being partly re-implemented, which
    /// is why the app registers the four asset collections and the two gizmo
    /// groups: `init_gizmo_group` only touches resources, so it needs no
    /// `GizmoPlugin` and no renderer.
    #[test]
    fn setup_dresses_the_five_walls_in_the_concrete_layer() {
        let (device, queue) = wgpu_pair();
        let mut app = App::new();
        app.add_plugins(MinimalPlugins)
            .add_plugins(bevy::asset::AssetPlugin::default())
            .init_asset::<Mesh>()
            .init_asset::<Image>()
            .init_asset::<StandardMaterial>()
            .init_asset::<TerrainMaterial>()
            .init_gizmo_group::<ChunkGizmos>()
            .init_gizmo_group::<GhostGizmos>()
            .init_resource::<ViewFlags>()
            .insert_resource(AutoCarve::default())
            .insert_resource(RenderDevice::from(device))
            .insert_resource(RenderQueue(Arc::new(WgpuWrapper::new(queue))));
        app.world_mut()
            .run_system_once(setup)
            .expect("setup runs with no window and no renderer");

        let mut slabs = app
            .world_mut()
            .query_filtered::<(&MeshMaterial3d<TerrainMaterial>, &Transform), With<Wall>>();
        let spawned: Vec<_> = slabs
            .iter(app.world())
            .map(|(m, t)| (m.0.clone(), t.translation))
            .collect();
        assert_eq!(
            spawned.len(),
            5,
            "four sides and a floor, and no ceiling; got {}",
            spawned.len()
        );

        let terrain = app.world().resource::<SurfaceMaterial>().0.clone();
        let materials = app.world().resource::<Assets<TerrainMaterial>>();
        assert_eq!(
            materials
                .get(&terrain)
                .expect("the terrain material")
                .extension
                .settings
                .z,
            LAYER_BLEND,
            "the terrain is not blending its layers"
        );
        for (handle, centre) in &spawned {
            assert_ne!(
                *handle, terrain,
                "the wall at {centre:?} wears the terrain material, so the boundary \
                 is grass and dirt instead of concrete"
            );
            assert_eq!(
                materials
                    .get(handle)
                    .expect("the wall material")
                    .extension
                    .settings
                    .z,
                LAYER_CONCRETE,
                "the wall at {centre:?} is not forced to the concrete layer"
            );
        }
        // The same five centres [`walls`] computes, so the spawn cannot drift
        // from the geometry the other test gates.
        let layout = app.world().resource::<World>().layout;
        let mut want: Vec<Vec3> = walls(&layout).into_iter().map(|(c, _)| c).collect();
        let mut got: Vec<Vec3> = spawned.into_iter().map(|(_, c)| c).collect();
        let key = |v: &Vec3| (v.x.to_bits(), v.y.to_bits(), v.z.to_bits());
        want.sort_unstable_by_key(key);
        got.sort_unstable_by_key(key);
        assert_eq!(
            got, want,
            "the spawned walls are not the boxes `walls` names"
        );
    }

    /// Walking hard at the `+x` wall stops at it.
    ///
    /// `Ground` has a height at every `x` and `z`, so before the clamp this walk
    /// carried on over invisible ground outside the box -- past the wall the
    /// player had just watched themselves walk through. `Shift` and 300 frames
    /// is 43 units of intent against an 8-unit box, so a body that is not
    /// clamped ends up nowhere near the assertion.
    #[test]
    fn walking_into_the_wall_stops_at_the_sandbox() {
        let mut app = harness(true);
        for _ in 0..600 {
            step(&mut app);
            if app.world().resource::<World>().backlog == 0
                && app.world().resource::<World>().grounded
            {
                break;
            }
        }
        // Yaw and pitch are zero in the harness, so `right` is `+x` and `D` is a
        // walk straight at the wall.
        {
            let mut keys = app.world_mut().resource_mut::<ButtonInput<KeyCode>>();
            keys.press(KeyCode::KeyD);
            keys.press(KeyCode::ShiftLeft);
        }
        for _ in 0..300 {
            app.update();
        }
        let layout = app.world().resource::<World>().layout;
        let (_, hi) = sandbox(&layout);
        let eye = eye(&mut app);
        assert!(
            eye.x <= hi.x - BODY_RADIUS + 1e-4,
            "the body walked out of the sandbox: x {} against a wall at {}",
            eye.x,
            hi.x
        );
        // And it really did reach the wall, so this is a clamp rather than a
        // body that never moved.
        assert!(
            eye.x > hi.x - BODY_RADIUS - 0.05,
            "the body never reached the wall: x {}",
            eye.x
        );
    }

    /// The virtual stick's shape, without a phone.
    ///
    /// Four properties, each a bug someone ships: creep from a resting thumb, a
    /// stick that never reaches full speed, inverted forward, and a corner drag
    /// that is `sqrt(2)` times faster than a straight one.
    #[test]
    fn the_virtual_stick_has_a_dead_zone_a_unit_disc_and_screen_y_flipped() {
        let origin = Vec2::new(100.0, 400.0);
        assert_eq!(
            touch_axes(origin, origin + Vec2::new(4.0, 0.0)),
            Vec2::ZERO,
            "a 4 px twitch moved the body"
        );
        assert_eq!(
            touch_axes(origin, origin),
            Vec2::ZERO,
            "a stationary thumb is not centred"
        );
        let far = touch_axes(origin, origin + Vec2::new(200.0, 0.0));
        assert!(
            (far.length() - 1.0).abs() < 1e-6,
            "a 200 px drag gave {far:?}, not full deflection"
        );
        // Screen `y` grows downward, so dragging the thumb *down* must walk
        // backwards.
        let down = touch_axes(origin, origin + Vec2::new(0.0, 200.0));
        assert!(down.y < 0.0, "dragging down walked forwards: {down:?}");
        let up = touch_axes(origin, origin + Vec2::new(0.0, -200.0));
        assert!(up.y > 0.0, "dragging up walked backwards: {up:?}");
        let corner = touch_axes(origin, origin + Vec2::new(200.0, 200.0));
        assert!(
            (corner.length() - 1.0).abs() < 1e-6,
            "a corner drag gave {} rather than 1.0, so diagonals are faster",
            corner.length()
        );
    }
}
