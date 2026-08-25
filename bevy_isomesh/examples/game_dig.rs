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
//! cannot see a chunk boundary. The shader and both packed 1k textures are
//! **compiled in** (`include_bytes!`, `load_internal_asset!`), because nothing
//! copies an `assets/` tree into `web/dist` and a run-time load path would work
//! natively and 404 in the browser.
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
    ExtractTimings, FieldBuffer, GridParams, MarchingCubesGpu, Readback, read_bytes_many_deferred,
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

/// Gizmos for the re-meshed-chunk outline.
#[derive(Default, Reflect, GizmoConfigGroup)]
struct ChunkGizmos;

/// The terrain before any edit: a slab with a rolling top.
///
/// Hand-rolled rather than `FbmTerrain`, because this needs a floor a player can
/// stand on and a ceiling to dig into, and it must be cheap — it is sampled
/// inside the edit loop.
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

/// Radius of each of the two spheres the body is made of.
const BODY_RADIUS: f32 = 0.4;
/// Sphere centres below the eye.
///
/// The eye sits at the top of the upper sphere and the lower one rests on the
/// ground, so standing on flat terrain puts the eye at `1.2 + 0.4 = 1.6` --
/// exactly where `setup` already places the camera, so walk mode opens on the
/// viewpoint this demo has always had.
const BODY_OFFSETS: [f32; 2] = [0.4, 1.2];
/// Downward acceleration. Roughly twice Earth's, which is the usual game figure:
/// real gravity over a 1.6-unit body reads as floating.
const GRAVITY: f32 = 18.0;
/// Launch speed. `8.5^2 / (2 * 18)` is a **2.0-unit apex**, which is a whole
/// chunk: high enough to jump onto the lip of a pit dug with a brush at the
/// large end of the wheel's range, rather than the 1.0 the first version had.
const JUMP_SPEED: f32 = 8.5;
/// Resolution passes per frame. Two, because a body wedged in a corner is pushed
/// out of one sphere into the other and needs a second look.
const RESOLVE_PASSES: u32 = 2;
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
    /// CPU field evaluation plus the upload, measured around
    /// `FieldBuffer::sampled`.
    upload_ms: f64,
    /// The three dispatches, as `isomesh-gpu` measures them.
    timings: ExtractTimings,
}

#[derive(Resource)]
struct Look {
    yaw: f32,
    pitch: f32,
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

/// The triplanar half of the terrain material: two packed textures and the two
/// numbers that place them.
#[derive(Asset, AsBindGroup, Reflect, Debug, Clone)]
struct TriplanarExtension {
    /// `x` world units per texture tile, `y` blend sharpness. See the WGSL for
    /// why this is one `vec4` and not two `f32`s.
    #[uniform(100)]
    settings: Vec4,
    /// RGB colour, A roughness.
    #[texture(101)]
    #[sampler(102)]
    albedo_roughness: Handle<Image>,
    /// RGB OpenGL-convention normal, A ambient occlusion.
    #[texture(103)]
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

#[derive(Resource)]
struct SurfaceMaterial(Handle<TerrainMaterial>);

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
    .insert_resource(Look {
        yaw: 0.0,
        pitch: -0.15,
    })
    .insert_resource(AutoCarve::from_env())
    .add_systems(Startup, setup)
    .add_systems(
        Update,
        // Chained, because the order is load-bearing and a tuple does not
        // impose one: `switch_algorithm` marks the sandbox before the frame's
        // aim and edit rather than racing them, `move_camera` moves the camera
        // and resolves the body against the field, `aim` traces from where it now
        // is, `dig` marks at that point, `drain_dirty` meshes what the frame's
        // budget allows -- after `dig`, so an edit's one to eight chunks clear in
        // the same frame -- `gpu_collect` finishes whatever the GPU returned, and
        // `ghost` draws the same point. Unchained, `dig` read a camera transform
        // one frame stale -- which it did before this, invisibly.
        // `loading_modal` is last so it reads the post-drain backlog.
        (
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
            .insert(Transform::from_xyz(0.0, 1.6, 6.0));
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

    // Compiled in rather than loaded, for the reason the module docs give:
    // nothing copies an `assets/` tree into `web/dist`, so an
    // `AssetServer::load` path would work natively and 404 in the browser.
    // `Image::from_buffer` decodes here on the main thread once, which is
    // 1024x1024 of PNG twice.
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
        .expect("the packed ground textures are committed beside this example")
    }

    let albedo_roughness = images.add(embedded_texture(
        include_bytes!("textures/ground_0046_albedo_roughness.png"),
        true,
    ));
    let normal_ao = images.add(embedded_texture(
        include_bytes!("textures/ground_0046_normal_ao.png"),
        false,
    ));
    let material = terrain_materials.add(TerrainMaterial {
        base: StandardMaterial {
            // White, not the old 0.62/0.58/0.52 tint: the shader multiplies the
            // sampled colour by this, so a tint here would darken the texture.
            // The field stays as the knob it is, set to neutral.
            base_color: Color::WHITE,
            perceptual_roughness: 0.85,
            ..default()
        },
        extension: TriplanarExtension {
            // 1.5 world units per tile against a 2.0-unit chunk, so a tile is a
            // little smaller than a chunk and the cracks read at the scale a
            // brush of radius 0.25 cuts at. `y` is the blend sharpness.
            settings: Vec4::new(1.5, 4.0, 0.0, 0.0),
            albedo_roughness,
            normal_ao,
        },
    });

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
/// the **extraction**: a compute shader classifies every cell, a prefix scan
/// allocates, and a second dispatch writes the triangles, with the only CPU work
/// being the field sample and the upload. Same algorithm as key `1`, so the
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
            // Both paths pay it, which is what makes the eight keys comparable:
            // the field cost is the same and what differs is the extraction.
            let span = layout.cell_size() * layout.cells() as f32;
            let box_ = ChunkBox::new(origin, span);
            prune_into(brushes, &Ground, box_, survivors);
            let field = BrushStack {
                base: Ground,
                brushes: survivors,
            };
            if algorithm == Algorithm::MarchingCubesGpu {
                if let Some((job, cost)) = gpu_dispatch(
                    &gpu.mc,
                    device.wgpu_device(),
                    &queue,
                    &layout,
                    &field,
                    id,
                    origin,
                ) {
                    pending.jobs.push(job);
                    in_flight.set(in_flight.get() + 1);
                    gpu_cost = Some(cost);
                }
                return;
            }
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

/// Sample the chunk's field on the CPU, upload it, extract on the GPU, and start
/// the readback.
///
/// `FieldBuffer::sampled`, not a device-side field: `field.wgsl`'s base is a
/// `switch` over four hard-coded reference fields and this demo's base is a
/// custom terrain plus an edit tape the shader cannot evaluate. Sampling on the
/// CPU is also what makes this a fair comparison -- all eight keys pay the same
/// field cost against the same pruned tape, and what differs is the extraction.
/// The HUD prints the upload share, so the crate's own "the upload dominates"
/// finding is on screen rather than in a docstring.
///
/// `extract_indirect`, not `extract` or `extract_buffers`: it is the **only**
/// entry point that reads nothing back. The others end in
/// `read_buffer_u32` -> `read_bytes_many` -> `device.poll(PollType::Wait)` +
/// `mpsc::recv()`, and on WebGPU `Device::poll` is a documented no-op while
/// `map_async` completes only from the browser's event loop -- so that wait
/// wedges the tab rather than blocking a thread.
fn gpu_dispatch(
    mc: &MarchingCubesGpu,
    device: &wgpu::Device,
    queue: &wgpu::Queue,
    layout: &ChunkLayout<f32>,
    field: &impl Sdf<Scalar = f32>,
    id: ChunkId,
    origin: [f32; 3],
) -> Option<(GpuJob, GpuCost)> {
    // No halo: `halo_cells` is zero for Marching Cubes, which marches every cell
    // and puts vertices on grid edges both chunks compute identically.
    let params = GridParams::new([layout.cells() + 1; 3], origin, layout.cell_size()).ok()?;
    let started = Instant::now();
    let buffer = FieldBuffer::sampled(device, queue, params, field).ok()?;
    let upload_ms = started.elapsed().as_secs_f64() * 1000.0;
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
            upload_ms,
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
    let floats = |bytes: &[u8]| -> Vec<[f32; 3]> {
        bytes
            .chunks_exact(12)
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

fn grab(
    keys: Res<ButtonInput<KeyCode>>,
    buttons: Res<ButtonInput<MouseButton>>,
    // `CursorOptions` is its own component on the window entity in Bevy 0.19,
    // not a field of `Window`.
    mut cursor: Query<&mut CursorOptions, With<PrimaryWindow>>,
    mut world: ResMut<World>,
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
    let (mode, visible) = if world.grabbed {
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
/// Two spheres against the field, [`RESOLVE_PASSES`] times. `f` is not an exact
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
    let mut contact = None;
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
            *eye += n * (BODY_RADIUS - f);
            contact = Some(n);
        }
    }
    // Otherwise gravity accumulates against the floor and the first step off a
    // ledge is a plummet. Only the component *into* the surface goes: sliding
    // along it is movement the player asked for.
    if let Some(n) = contact {
        let into = velocity.dot(n);
        if into < 0.0 {
            *velocity -= n * into;
        }
    }
    // A separate downward probe, not `contact.is_some()`: a body pressed against
    // a wall is in contact and is not standing on anything. The lowest sphere is
    // the one with the largest offset, found rather than indexed so the
    // contingency of adding another offset cannot silently probe the wrong one.
    let lowest = BODY_OFFSETS.into_iter().fold(f32::MIN, f32::max);
    let foot = *eye - Vec3::Y * lowest;
    field.sample([foot.x, foot.y - (BODY_RADIUS + GROUND_PROBE), foot.z]) <= 0.0
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
    // `Space`, which is what a hand reaches for. The shared harness also reads
    // it into `ViewFlags::paused`, and that is harmless here: `paused` is read
    // only by `orbit_camera`, and `setup` takes `OrbitCamera` off this camera so
    // that system's query is empty in this demo.
    if world.grounded && keys.just_pressed(KeyCode::Space) {
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
    let left = buttons.pressed(MouseButton::Left);
    let right = buttons.pressed(MouseButton::Right);
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
            // Upload first, because it is the biggest of the four and that is
            // the crate's own committed finding rather than a surprise: the
            // field is sampled on the CPU and shipped, and the three dispatches
            // are cheap beside it.
            Some(cost) => format!(
                "           upload {:.2} ms   count {:.2}   scan {:.2}   emit {:.2}   peak {} tris/chunk of {}{}",
                cost.upload_ms,
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
        "ground     ground_0046 triplanar, 1.5 u/tile, sharpness 4".to_string(),
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
    use std::sync::Arc;

    use bevy::asset::AssetApp;
    use bevy::mesh::VertexAttributeValues;
    use bevy::render::renderer::WgpuWrapper;
    use bevy::time::TimeUpdateStrategy;

    use super::*;

    /// A frame, fixed. `TimeUpdateStrategy::ManualDuration` is what makes
    /// `time.delta_secs()` this rather than however long the test machine took,
    /// which is what makes an edit count assertable.
    const FRAME: Duration = Duration::from_millis(16);

    /// A device and queue of this test binary's own.
    ///
    /// `isomesh_gpu::headless::Gpu` lends its device and `RenderDevice` needs to
    /// own one, so this opens a second. No software fallback, the same policy
    /// that module states: a test that quietly ran on a CPU reference driver
    /// would report a GPU comparison that is not one.
    fn wgpu_pair() -> (wgpu::Device, wgpu::Queue) {
        let instance =
            wgpu::Instance::new(wgpu::InstanceDescriptor::new_without_display_handle_from_env());
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
            .spawn((Camera3d::default(), Transform::from_xyz(0.0, 1.6, 6.0)));
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
        assert!(cost.upload_ms > 0.0, "the upload was not measured");
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
        println!("the GPU switch drained in {frames} frames, {gpu_resident} chunks resident");
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
    /// `1.6` above the surface -- the height the two spheres put it, which is
    /// exactly where `setup` has always placed the camera.
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
    /// `1.2² / 80 = 18 mm` over the body's height, which is a wall.
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
        let mut eye = Vec3::new(0.0, 1.6, 0.0);
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
        // The wall's face is at z = -0.8 and the body is a sphere of
        // `BODY_RADIUS`, so its centre cannot pass -0.4 by more than the field's
        // own slack.
        assert!(
            eye.z > -0.5,
            "the body walked into the rock: z {} (a free walk reaches -2.25)",
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
        let mut eye = Vec3::new(0.0, 1.6, -2.0);
        let mut velocity = Vec3::ZERO;
        let mut grounded = false;
        for _ in 0..240 {
            gravity_step(grounded, &mut velocity, dt);
            eye.y += velocity.y * dt;
            grounded = resolve_body(&pit, &mut eye, &mut velocity);
        }
        // The pit's floor is the bottom of the subtracted sphere, y = -1.2, so
        // the eye lands well below the 1.6 it stands at on the surface.
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
    /// Marching Cubes is Marching Cubes: the device-side extraction reads the
    /// same samples `FieldBuffer::sampled` uploaded, so the triangle count has to
    /// agree with the CPU extractor's to within float ordering at a crossing.
    /// This is the check a screenshot of key `8` can only suggest.
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
            let (job, _) =
                gpu_dispatch(&mc, gpu.device(), gpu.queue(), &layout, &field, id, origin)
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
}
