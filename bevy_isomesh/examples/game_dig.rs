//! E-202 — carving tunnels, the way a game does it.
//!
//! ```bash
//! cd bevy_isomesh && cargo run --example game_dig --release
//! ```
//!
//! `WASD`/`QE` fly, mouse looks, and a translucent sphere on the rock under the
//! crosshair is the brush that a click would push — orange to carve, cyan to
//! fill. **Hold** the left button to keep carving along the sweep, the right to
//! keep filling. The wheel or `[`/`]` resizes the brush, `1`–`7` swap the
//! extractor and print what the re-mesh cost, `Z` undoes one brush, `X` clears
//! the log, `C` outlines the chunks the last edit re-meshed, and `Tab` releases
//! the cursor.
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

use bevy::asset::{RenderAssetUsages, load_internal_asset, uuid_handle};
use bevy::image::{
    CompressedImageFormats, ImageAddressMode, ImageSampler, ImageSamplerDescriptor, ImageType,
};
use bevy::input::mouse::{AccumulatedMouseMotion, AccumulatedMouseScroll, MouseScrollUnit};
use bevy::pbr::{ExtendedMaterial, MaterialExtension, MaterialPlugin};
use bevy::platform::time::Instant;
use bevy::prelude::*;
use bevy::render::render_resource::AsBindGroup;
use bevy::shader::{Shader, ShaderRef};
use bevy::window::{CursorGrabMode, CursorOptions, PrimaryWindow};
use bevy_isomesh::MeshBuilder;
use common::{Capture, CommonPlugin, DemoStats, OrbitCamera};
use isomesh::Sdf;
use isomesh::brush::{Brush, BrushStack};
use isomesh::chunk::dirty::{DirtySet, EditReport, mark_edit};
use isomesh::chunk::{ChunkId, ChunkLayout};
use isomesh::dual_contouring::DualContouring;
use isomesh::fields::Sphere;
use isomesh::greedy_quads::GreedyQuads;
use isomesh::manifold_dual_contouring::ManifoldDualContouring;
use isomesh::marching_cubes::{FaceAmbiguity, MarchingCubes};
use isomesh::marching_tetrahedra::MarchingTetrahedra;
use isomesh::surface_nets::SurfaceNets;

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

/// Where the next edit lands, recomputed every frame and read by both the ghost
/// and the edit. One point, two consumers: a preview that can disagree with the
/// brush it previews is worse than no preview.
#[derive(Resource, Default)]
struct Aim {
    point: Vec3,
    /// `false` when the ray leaves the sandbox without crossing the surface.
    hit: bool,
}

/// Marker for the translucent brush preview.
#[derive(Component)]
struct Ghost;

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
    dirty: DirtySet,
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
    /// Centre of the last brush pushed by a held button, or `None` when no
    /// button is down. This is what makes a hold a *stroke* rather than a burst.
    stroke_last: Option<Vec3>,
    stroke_edits: u32,
    stroke_ms: f64,
}

#[derive(Component)]
struct Chunk(ChunkId);

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
    .insert_resource(Look {
        yaw: 0.0,
        pitch: -0.15,
    })
    .insert_resource(AutoCarve::from_env())
    .add_systems(Startup, setup)
    .add_systems(
        Update,
        // Chained, because the order is load-bearing and a tuple does not
        // impose one: `switch_algorithm` re-meshes the sandbox before the
        // frame's aim and edit rather than racing them, `fly` moves the
        // camera, `aim` traces from where it now is, `dig` edits at that
        // point and `ghost` draws the same point. Unchained, `dig` read a
        // camera transform one frame stale -- which it did before this,
        // invisibly.
        (
            grab,
            switch_algorithm,
            fly,
            aim,
            dig,
            ghost,
            report,
            outline_chunks,
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
    chunks: Query<(Entity, &Chunk)>,
) {
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
        dirty: DirtySet::new(),
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
        stroke_last: None,
        stroke_edits: 0,
        stroke_ms: 0.0,
    };

    // Mesh every chunk once. After this, only edited chunks are ever re-meshed.
    for z in 0..EXTENT[2] {
        for y in 0..EXTENT[1] {
            for x in 0..EXTENT[0] {
                world.dirty.insert(ChunkId::new([x, y, z]));
            }
        }
    }
    // The query is empty on this first frame -- `Commands` are deferred and
    // nothing has spawned a `Chunk` yet -- so every arm `rebuild` takes is
    // `(None, Some(_))`, which is the spawn this used to do inline.
    let layout = world.layout;
    let algorithm = world.algorithm;
    {
        let field = BrushStack {
            base: Ground,
            brushes: &world.brushes,
        };
        rebuild(
            &mut commands,
            &mut meshes,
            &material,
            &chunks,
            &layout,
            &field,
            algorithm,
            &mut world.dirty,
        );
    }

    commands.insert_resource(SurfaceMaterial(material));
    commands.insert_resource(world);
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

/// Which extractor meshes a chunk. Seven of the crate's eight, and the eighth is
/// named here rather than silently missing.
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
}

impl Algorithm {
    /// In key order: `1` is the reference and the default.
    const ALL: [Self; 7] = [
        Self::MarchingCubes,
        Self::MarchingCubesDecider,
        Self::MarchingTetrahedra,
        Self::SurfaceNets,
        Self::DualContouring,
        Self::ManifoldDualContouring,
        Self::GreedyQuads,
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
        }
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

    /// Extract one chunk, or `None` when it holds no surface -- so empty air
    /// costs a sample loop rather than an entity and a draw call over nothing.
    fn mesh<F: Sdf<Scalar = f32>>(
        self,
        layout: &ChunkLayout<f32>,
        field: &F,
        origin: [f32; 3],
    ) -> Option<Mesh> {
        let shape = layout.sample_shape().ok()?;
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
        };
        extracted.ok()?;
        if builder.indices().is_empty() {
            return None;
        }
        Some(builder.into_mesh())
    }
}

/// Re-mesh every dirty chunk and reconcile the entity for it. Returns how many
/// chunks were rebuilt.
///
/// The four arms are the whole state machine: a chunk that had a mesh and still
/// has one is replaced, one that lost its surface is despawned rather than left
/// as an empty draw call, one that gained a surface is spawned, and one that
/// never had either costs nothing.
///
/// One function, because there were two divergent copies of it -- `setup`'s
/// spawn-only loop and `dig`'s four-way reconcile -- and the algorithm switch
/// would have been a third.
#[allow(clippy::too_many_arguments)]
fn rebuild<F: Sdf<Scalar = f32>>(
    commands: &mut Commands,
    meshes: &mut Assets<Mesh>,
    material: &Handle<TerrainMaterial>,
    chunks: &Query<(Entity, &Chunk)>,
    layout: &ChunkLayout<f32>,
    field: &F,
    algorithm: Algorithm,
    dirty: &mut DirtySet,
) -> usize {
    // `mesh_dirty` returns the count itself and clears the set afterwards, so
    // there is nothing to tally in the closure and nothing to clear after it.
    dirty.mesh_dirty(layout, |id, origin| {
        let mesh = algorithm.mesh(layout, field, origin).map(|m| meshes.add(m));
        let existing = chunks.iter().find(|(_, c)| c.0 == id).map(|(e, _)| e);
        match (existing, mesh) {
            (Some(entity), Some(handle)) => {
                commands.entity(entity).insert(Mesh3d(handle));
            }
            (Some(entity), None) => {
                commands.entity(entity).despawn();
            }
            (None, Some(handle)) => {
                commands.spawn((Mesh3d(handle), MeshMaterial3d(material.clone()), Chunk(id)));
            }
            (None, None) => {}
        }
    })
}

/// `1`-`7` re-mesh the whole sandbox with a different extractor.
///
/// The whole sandbox, not the dirty set: nothing about the *field* changed, so
/// `mark_edit` would correctly report zero dirty chunks. And re-meshing 256
/// chunks with each algorithm in turn is the comparison -- the HUD prints what
/// it cost, which is the only honest way to show that Marching Tetrahedra is
/// three times the triangles.
///
/// `handle_keys` in the shared harness also reads `Digit1`-`Digit7`, into
/// `ViewFlags::field`. This example never reads that field, so there is nothing
/// to collide with.
fn switch_algorithm(
    keys: Res<ButtonInput<KeyCode>>,
    mut world: ResMut<World>,
    mut commands: Commands,
    mut meshes: ResMut<Assets<Mesh>>,
    material: Res<SurfaceMaterial>,
    chunks: Query<(Entity, &Chunk)>,
) {
    const DIGITS: [KeyCode; 7] = [
        KeyCode::Digit1,
        KeyCode::Digit2,
        KeyCode::Digit3,
        KeyCode::Digit4,
        KeyCode::Digit5,
        KeyCode::Digit6,
        KeyCode::Digit7,
    ];
    let Some(index) = DIGITS.iter().position(|k| keys.just_pressed(*k)) else {
        return;
    };
    let chosen = Algorithm::ALL[index];
    if chosen == world.algorithm {
        return;
    }
    world.algorithm = chosen;

    let started = Instant::now();
    // Copied out before `field` borrows `world.brushes`, so the timing writes
    // after the `rebuild` call have nothing to fight over. `ChunkLayout<f32>` is
    // `Copy`, which is what makes that free.
    let layout = world.layout;
    let mut dirty = DirtySet::new();
    for z in 0..EXTENT[2] {
        for y in 0..EXTENT[1] {
            for x in 0..EXTENT[0] {
                dirty.insert(ChunkId::new([x, y, z]));
            }
        }
    }
    let rebuilt = {
        let field = BrushStack {
            base: Ground,
            brushes: &world.brushes,
        };
        rebuild(
            &mut commands,
            &mut meshes,
            &material.0,
            &chunks,
            &layout,
            &field,
            chosen,
            &mut dirty,
        )
    };
    world.switch_chunks = rebuilt;
    world.switch_ms = started.elapsed().as_secs_f64() * 1000.0;
    world.last_touched.clear();
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

fn fly(
    keys: Res<ButtonInput<KeyCode>>,
    motion: Res<AccumulatedMouseMotion>,
    time: Res<Time>,
    world: Res<World>,
    mut look: ResMut<Look>,
    mut camera: Query<&mut Transform, With<Camera3d>>,
) {
    let Ok(mut transform) = camera.single_mut() else {
        return;
    };
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
        (KeyCode::KeyQ, -Vec3::Y),
        (KeyCode::KeyE, Vec3::Y),
    ] {
        if keys.pressed(key) {
            direction += delta;
        }
    }
    // 9.0, not the old 6.0: a 16-unit box takes 6.4 s to cross at the walk speed
    // and 1.8 s at this one.
    let speed = if keys.pressed(KeyCode::ShiftLeft) {
        9.0
    } else {
        2.5
    };
    if direction != Vec3::ZERO {
        transform.translation += direction.normalize() * speed * time.delta_secs();
    }
}

/// The loop this example exists for: one brush, one incremental re-mesh, for as
/// long as the button is held.
#[allow(clippy::too_many_arguments)]
fn dig(
    buttons: Res<ButtonInput<MouseButton>>,
    scroll: Res<AccumulatedMouseScroll>,
    keys: Res<ButtonInput<KeyCode>>,
    mut world: ResMut<World>,
    mut commands: Commands,
    mut meshes: ResMut<Assets<Mesh>>,
    material: Res<SurfaceMaterial>,
    target: Res<Aim>,
    chunks: Query<(Entity, &Chunk)>,
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
    let moved = world
        .stroke_last
        .is_none_or(|last| target.point.distance(last) >= world.radius * 0.5);
    let editable = world.grabbed && target.hit && moved;
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

    let mut dirty = DirtySet::new();
    // `before` and `after` are two slices of the one log, which is what makes
    // this exact: the two fields differ by precisely the brushes between them.
    // A push reads `[..split]` then `[..]`; a shrink reads them the other way
    // round and drops the tail *after* the re-mesh has sampled it.
    //
    // Truncating first would be the bug: `before` and `after` would both be bare
    // `Ground`, `mark_edit` would correctly report nothing changed, and the
    // carve would stay on screen beside an empty log.
    let (report, rebuilt, touched) = {
        let (before_log, after_log) = if shrink {
            (&world.brushes[..], &world.brushes[..split])
        } else {
            (&world.brushes[..split], &world.brushes[..])
        };
        let before = BrushStack {
            base: Ground,
            brushes: before_log,
        };
        let after = BrushStack {
            base: Ground,
            brushes: after_log,
        };
        let report = mark_edit(&layout, &before, &after, min_cell, max_cell, &mut dirty)
            .expect("a dig brush spans a few cells, far inside the u32 sample space");
        let touched: Vec<ChunkId> = dirty.iter().collect();
        let rebuilt = rebuild(
            &mut commands,
            &mut meshes,
            &material.0,
            &chunks,
            &layout,
            &after,
            algorithm,
            &mut dirty,
        );
        (report, rebuilt, touched)
    };
    if shrink {
        world.brushes.truncate(split);
        // The next hold starts a fresh stroke rather than measuring against a
        // brush that no longer exists.
        world.stroke_last = None;
    }

    world.last_edit_ms = started.elapsed().as_secs_f64() * 1000.0;
    // One line per scripted edit, so the log-growth claim in the module docs can
    // be checked from a terminal rather than by reading a HUD off a screenshot.
    if scripted.is_some() {
        info!(
            "edit {:>3}: {} chunks in {:.3} ms, E1 {:.1}%",
            world.brushes.len(),
            rebuilt,
            world.last_edit_ms,
            100.0 * report.changed_fraction()
        );
    }
    world.last_edit = Some(report);
    world.last_chunks = rebuilt;
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
            "algorithm  {}   (last switch: {} chunks in {:.1} ms)",
            world.algorithm.name(),
            world.switch_chunks,
            world.switch_ms
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
        "ground     ground_0046 triplanar, 1.5 u/tile, sharpness 4".to_string(),
        String::new(),
        "every field sample walks the log: measured 3.7x ms/chunk for 7x the log".to_string(),
    ];
    // The harness's shared footer advertises `[W] wire`, `[N] normals`,
    // `[G] grid` and `[R] re-mesh`. Every one of those is a lie here: `W` flies
    // forward, and the chunk entities carry `Chunk(ChunkId)` rather than
    // `DemoMesh`/`DemoDomain`, so the harness's wireframe, normal and domain
    // systems never see them.
    stats.keys = Some(
        "[LMB] hold to carve   [RMB] hold to fill   [WASD/QE] fly   [Shift] fast\n\
         [wheel] or [ ] brush   [1-7] algorithm   [Z] undo   [X] clear log   [C] chunks   [Tab] cursor   [F12] shot"
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
