//! Shared harness for every isomesh example.
//!
//! Orbit camera, HUD, and one set of keybindings, so that adding an example is
//! writing the thing it demonstrates and nothing else.
//!
//! # Keys
//!
//! | key | |
//! |---|---|
//! | `W` | wireframe |
//! | `N` | normals |
//! | `G` | grid / domain box |
//! | `Space` | pause |
//! | `R` | re-mesh |
//! | `F12` | screenshot |
//! | `Esc` | quit |
//!
//! Drag with the left mouse button to orbit, scroll to zoom.
//!
//! # Why the wireframe is drawn with gizmos
//!
//! Bevy ships a `WireframePlugin`, and it needs the wgpu features
//! `POLYGON_MODE_LINE | IMMEDIATES`. wgpu-types 29.0.4 does list Metal as
//! supporting `POLYGON_MODE_LINE`, so it is available here — but the plugin
//! *silently warns and disables itself* when a feature is missing, which turns a
//! toggle that does nothing into the hardest kind of bug to notice. Requesting
//! non-default features through `WgpuSettings` can also fail renderer
//! initialisation outright.
//!
//! Line gizmos need no features, behave identically on every backend, and draw
//! on top where you can actually see them. At demo triangle counts that costs
//! nothing worth measuring. Revisit if an example ever needs a wireframe over a
//! six-figure mesh.

use std::collections::VecDeque;

use bevy::asset::RenderAssetUsages;
use bevy::input::mouse::{AccumulatedMouseMotion, AccumulatedMouseScroll};
use bevy::mesh::{Indices, VertexAttributeValues};
use bevy::prelude::*;
use bevy::render::view::screenshot::{Screenshot, save_to_disk};

/// Everything an example gets for free.
pub struct CommonPlugin;

impl Plugin for CommonPlugin {
    fn build(&self, app: &mut App) {
        app.init_resource::<ViewFlags>()
            .init_resource::<DemoStats>()
            .init_resource::<FrameTimes>()
            .init_resource::<AutoScreenshot>()
            .init_resource::<Capture>()
            .init_resource::<Spin>()
            // PreStartup, not Startup: system order within a schedule is
            // unspecified, so an example's own `Startup` system that wants to
            // adjust the camera would sometimes run before the camera existed
            // and silently lose its settings. Spawning a schedule earlier makes
            // "the camera exists by the time your Startup runs" a guarantee
            // rather than a coin flip.
            .add_systems(PreStartup, (spawn_camera, spawn_light, spawn_hud))
            .add_systems(Update, (auto_screenshot, capture_sequence))
            .add_systems(
                Update,
                (
                    orbit_camera,
                    handle_keys,
                    update_frame_times,
                    update_hud,
                    draw_wireframe,
                    draw_normals,
                    draw_domain,
                ),
            );
    }
}

/// What the view toggles are currently set to.
#[derive(Resource, Debug)]
pub struct ViewFlags {
    pub wireframe: bool,
    pub normals: bool,
    pub grid: bool,
    /// Freezes automatic camera motion. Manual orbit still works, so a paused
    /// view can still be inspected.
    pub paused: bool,
    /// Set for one frame when `R` is pressed.
    pub remesh_requested: bool,
    /// Which reference field the example should show, when it offers a choice.
    pub field: usize,
    /// Hides the HUD entirely.
    ///
    /// For media. The HUD is the point of a screenshot -- the numbers are the
    /// evidence -- but a GIF meant to be looked at rather than read has the text
    /// sitting on top of the geometry, and a reader who wants the numbers can
    /// open the still.
    pub hud: bool,
}

impl Default for ViewFlags {
    /// Toggles can be preset from `ISOMESH_VIEW`, a comma-separated list of
    /// `wire`, `normals`, `nogrid`, `nohud`.
    ///
    /// That exists so a screenshot of a particular view can be captured without
    /// a human pressing a key, which is what makes the visual acceptance tests
    /// in the examples catalog reproducible rather than anecdotal.
    fn default() -> Self {
        let requested = std::env::var("ISOMESH_VIEW").unwrap_or_default();
        let has = |name: &str| requested.split(',').any(|part| part.trim() == name);
        Self {
            wireframe: has("wire"),
            normals: has("normals"),
            grid: !has("nogrid"),
            paused: false,
            remesh_requested: false,
            field: std::env::var("ISOMESH_FIELD")
                .ok()
                .and_then(|v| v.parse().ok())
                .unwrap_or(0),
            hud: !has("nohud"),
        }
    }
}

/// Resolution requested through `ISOMESH_SAMPLES`, in samples per axis.
///
/// The harness already lets a capture pick its view (`ISOMESH_VIEW`) and its
/// field (`ISOMESH_FIELD`) without a keyboard, and says that is what makes the
/// catalog's visual acceptance tests reproducible. **Resolution was the hole in
/// that claim**: the committed `e111` screenshots are of a 19³ grid and there
/// was no way to ask for one, so re-taking them meant pressing `[` by hand and
/// the images could not be regenerated from a command line at all.
///
/// Examples with a resolution read this once, at startup, and fall back to their
/// own default.
#[allow(dead_code)] // Each example compiles its own copy of this module.
#[must_use]
pub fn samples_override() -> Option<u32> {
    std::env::var("ISOMESH_SAMPLES")
        .ok()
        .and_then(|v| v.parse().ok())
}

/// Numbers the example wants on screen. Filled in by the example, rendered here.
#[derive(Resource, Default, Debug)]
pub struct DemoStats {
    pub title: String,
    pub vertices: usize,
    pub triangles: usize,
    /// Time the last extraction took.
    pub extract_ms: f64,
    /// Extra lines, one per entry.
    pub extra: Vec<String>,
}

/// A rolling window of frame times.
///
/// The HUD reports the **median** of the last 30, not the instantaneous value.
/// An instantaneous frame time is mostly noise and invites drawing conclusions
/// from a single sample, which is the habit the speed analysis in this repo
/// spends several pages arguing against.
#[derive(Resource, Default)]
struct FrameTimes(VecDeque<f64>);

const FRAME_WINDOW: usize = 30;

/// Takes a screenshot and quits, when `ISOMESH_SCREENSHOT` names a path.
///
/// A demo that only compiles is not a demo, and a windowed app cannot be
/// eyeballed from a terminal. This is what makes an example checkable: run it
/// headless-ish, get a PNG, look at the PNG. `docs/2026-08-11-bevy-examples-catalog.md`
/// is explicit that several of these examples are *visual acceptance tests* and
/// that the before/after pair belongs in the commit.
#[derive(Resource)]
struct AutoScreenshot {
    path: Option<String>,
    frames: u32,
}

impl Default for AutoScreenshot {
    fn default() -> Self {
        Self {
            path: std::env::var("ISOMESH_SCREENSHOT").ok(),
            frames: 0,
        }
    }
}

fn auto_screenshot(
    mut auto: ResMut<AutoScreenshot>,
    mut commands: Commands,
    mut exit: MessageWriter<AppExit>,
) {
    let Some(path) = auto.path.clone() else {
        return;
    };
    auto.frames += 1;
    // Long enough for the first extraction, asset upload and one full render.
    match auto.frames {
        90 => {
            commands
                .spawn(Screenshot::primary_window())
                .observe(save_to_disk(path.clone()));
            info!("screenshot -> {path}");
        }
        120 => {
            exit.write(AppExit::Success);
        }
        _ => {}
    }
}

/// Captures a numbered frame sequence, for turning an example into a GIF.
///
/// `ISOMESH_CAPTURE=<dir>` writes `frame_0000.png` onward and exits when the
/// sequence is complete. `ISOMESH_CAPTURE_FRAMES` and `ISOMESH_CAPTURE_EVERY`
/// tune length and stride.
///
/// A still shows what a mesh looks like; only a sequence shows what changing a
/// parameter *does*, which is the thing a reader of the README is actually
/// trying to find out.
#[derive(Resource)]
pub struct Capture {
    dir: Option<String>,
    total: u32,
    every: u32,
    /// Frames captured so far. Examples read this to drive a parameter sweep in
    /// step with the capture rather than with wall-clock time, so the sequence
    /// is reproducible.
    pub taken: u32,
    elapsed: u32,
    settle: u32,
}

impl Default for Capture {
    fn default() -> Self {
        let number = |key: &str, fallback: u32| {
            std::env::var(key)
                .ok()
                .and_then(|v| v.parse().ok())
                .unwrap_or(fallback)
        };
        Self {
            dir: std::env::var("ISOMESH_CAPTURE").ok(),
            total: number("ISOMESH_CAPTURE_FRAMES", 60),
            every: number("ISOMESH_CAPTURE_EVERY", 3).max(1),
            taken: 0,
            elapsed: 0,
            // Long enough for the window, the first extraction and the asset
            // upload to settle, so frame zero is not a grey rectangle.
            settle: number("ISOMESH_CAPTURE_SETTLE", 45),
        }
    }
}

impl Capture {
    /// Whether a sequence is being recorded.
    #[allow(dead_code)] // Each example compiles its own copy of this module.
    #[must_use]
    pub fn is_active(&self) -> bool {
        self.dir.is_some()
    }
}

fn capture_sequence(
    mut capture: ResMut<Capture>,
    mut commands: Commands,
    mut exit: MessageWriter<AppExit>,
) {
    let Some(dir) = capture.dir.clone() else {
        return;
    };
    capture.elapsed += 1;
    if capture.elapsed <= capture.settle {
        return;
    }
    if (capture.elapsed - capture.settle) % capture.every != 0 {
        return;
    }
    if capture.taken >= capture.total {
        exit.write(AppExit::Success);
        return;
    }
    let path = format!("{dir}/frame_{:04}.png", capture.taken);
    commands
        .spawn(Screenshot::primary_window())
        .observe(save_to_disk(path));
    capture.taken += 1;
}

/// Marks the mesh an example wants wireframed and normal-drawn.
#[derive(Component)]
pub struct DemoMesh;

/// The axis-aligned box an example is sampling, drawn when the grid is on.
#[derive(Component)]
pub struct DemoDomain {
    pub min: Vec3,
    pub max: Vec3,
}

/// Orbits a focus point.
#[derive(Component)]
pub struct OrbitCamera {
    pub focus: Vec3,
    pub radius: f32,
    pub yaw: f32,
    pub pitch: f32,
}

impl Default for OrbitCamera {
    fn default() -> Self {
        Self {
            focus: Vec3::ZERO,
            radius: 7.0,
            yaw: 0.6,
            pitch: 0.4,
        }
    }
}

/// Radians of automatic yaw per frame, from `ISOMESH_SPIN`.
#[derive(Resource)]
struct Spin(f32);

impl Default for Spin {
    fn default() -> Self {
        Self(
            std::env::var("ISOMESH_SPIN")
                .ok()
                .and_then(|v| v.parse().ok())
                .unwrap_or(0.0),
        )
    }
}

#[derive(Component)]
struct HudText;

fn spawn_camera(mut commands: Commands) {
    commands.spawn((Camera3d::default(), OrbitCamera::default()));
}

fn spawn_light(mut commands: Commands) {
    commands.spawn((
        DirectionalLight {
            illuminance: 8_000.0,
            shadow_maps_enabled: true,
            ..default()
        },
        Transform::from_xyz(6.0, 10.0, 6.0).looking_at(Vec3::ZERO, Vec3::Y),
    ));
    // In 0.19 `AmbientLight` became a per-camera component; the global default
    // is this resource.
    commands.insert_resource(GlobalAmbientLight {
        brightness: 220.0,
        ..default()
    });
}

fn spawn_hud(mut commands: Commands) {
    commands.spawn((
        Text::new(""),
        TextFont {
            // 0.19 made this a unit-carrying enum rather than a bare float.
            font_size: FontSize::Px(13.0),
            ..default()
        },
        TextColor(Color::srgb(0.92, 0.94, 0.98)),
        Node {
            position_type: PositionType::Absolute,
            top: Val::Px(10.0),
            left: Val::Px(12.0),
            ..default()
        },
        HudText,
    ));
}

fn orbit_camera(
    mouse: Res<ButtonInput<MouseButton>>,
    motion: Res<AccumulatedMouseMotion>,
    scroll: Res<AccumulatedMouseScroll>,
    spin: Res<Spin>,
    flags: Res<ViewFlags>,
    mut query: Query<(&mut OrbitCamera, &mut Transform)>,
) {
    for (mut orbit, mut transform) in &mut query {
        if !flags.paused {
            orbit.yaw += spin.0;
        }
        if mouse.pressed(MouseButton::Left) {
            orbit.yaw -= motion.delta.x * 0.005;
            orbit.pitch = (orbit.pitch - motion.delta.y * 0.005).clamp(-1.5, 1.5);
        }
        if scroll.delta.y != 0.0 {
            orbit.radius = (orbit.radius * (1.0 - scroll.delta.y * 0.1)).clamp(1.0, 200.0);
        }

        let direction = Vec3::new(
            orbit.yaw.cos() * orbit.pitch.cos(),
            orbit.pitch.sin(),
            orbit.yaw.sin() * orbit.pitch.cos(),
        );
        transform.translation = orbit.focus + direction * orbit.radius;
        transform.look_at(orbit.focus, Vec3::Y);
    }
}

fn handle_keys(
    keys: Res<ButtonInput<KeyCode>>,
    mut flags: ResMut<ViewFlags>,
    mut commands: Commands,
    mut exit: MessageWriter<AppExit>,
) {
    flags.remesh_requested = keys.just_pressed(KeyCode::KeyR);

    if keys.just_pressed(KeyCode::KeyW) {
        flags.wireframe = !flags.wireframe;
    }
    if keys.just_pressed(KeyCode::KeyN) {
        flags.normals = !flags.normals;
    }
    if keys.just_pressed(KeyCode::KeyG) {
        flags.grid = !flags.grid;
    }
    if keys.just_pressed(KeyCode::Space) {
        flags.paused = !flags.paused;
    }
    if keys.just_pressed(KeyCode::Escape) {
        exit.write(AppExit::Success);
    }

    for (key, index) in [
        (KeyCode::Digit1, 0),
        (KeyCode::Digit2, 1),
        (KeyCode::Digit3, 2),
        (KeyCode::Digit4, 3),
        (KeyCode::Digit5, 4),
        (KeyCode::Digit6, 5),
        (KeyCode::Digit7, 6),
    ] {
        if keys.just_pressed(key) {
            flags.field = index;
        }
    }

    if keys.just_pressed(KeyCode::F12) {
        let path = format!("screenshot-{}.png", std::process::id());
        commands
            .spawn(Screenshot::primary_window())
            .observe(save_to_disk(path.clone()));
        info!("screenshot -> {path}");
    }
}

fn update_frame_times(time: Res<Time>, mut times: ResMut<FrameTimes>) {
    let ms = time.delta_secs_f64() * 1000.0;
    if ms > 0.0 {
        times.0.push_back(ms);
        while times.0.len() > FRAME_WINDOW {
            times.0.pop_front();
        }
    }
}

fn median(values: &VecDeque<f64>) -> f64 {
    if values.is_empty() {
        return 0.0;
    }
    let mut sorted: Vec<f64> = values.iter().copied().collect();
    sorted.sort_by(f64::total_cmp);
    sorted[sorted.len() / 2]
}

fn update_hud(
    stats: Res<DemoStats>,
    flags: Res<ViewFlags>,
    times: Res<FrameTimes>,
    mut query: Query<&mut Text, With<HudText>>,
) {
    if !flags.hud {
        for mut hud in &mut query {
            hud.0.clear();
        }
        return;
    }

    let frame_ms = median(&times.0);
    let fps = if frame_ms > 0.0 {
        1000.0 / frame_ms
    } else {
        0.0
    };

    let mut text = String::new();
    text.push_str(&stats.title);
    text.push('\n');
    text.push_str(&format!(
        "\n{:>9} vertices\n{:>9} triangles\n{:>9.3} ms extract",
        stats.vertices, stats.triangles, stats.extract_ms
    ));
    text.push_str(&format!(
        "\n{frame_ms:>9.2} ms/frame (median of {FRAME_WINDOW})\n{fps:>9.0} fps"
    ));
    for line in &stats.extra {
        text.push('\n');
        text.push_str(line);
    }
    text.push_str(&format!(
        "\n\n[W] wire {}   [N] normals {}   [G] grid {}\n[Space] {}   [R] re-mesh   [F12] shot   [Esc] quit",
        on_off(flags.wireframe),
        on_off(flags.normals),
        on_off(flags.grid),
        if flags.paused { "resume" } else { "pause" },
    ));

    for mut target in &mut query {
        target.0.clone_from(&text);
    }
}

fn on_off(value: bool) -> &'static str {
    if value { "on" } else { "off" }
}

/// Reads the mesh back out of the asset and draws every triangle edge.
///
/// This is why [`MeshBuilder`](bevy_isomesh::MeshBuilder) asks for
/// [`RenderAssetUsages::default()`] rather than `RENDER_WORLD` alone — the mesh
/// has to still be in main memory to be read.
fn draw_wireframe(
    flags: Res<ViewFlags>,
    meshes: Res<Assets<Mesh>>,
    query: Query<(&Mesh3d, &GlobalTransform), With<DemoMesh>>,
    mut gizmos: Gizmos,
) {
    if !flags.wireframe {
        return;
    }
    let colour = Color::srgb(0.15, 0.95, 0.55);
    for (handle, transform) in &query {
        let Some(mesh) = meshes.get(&handle.0) else {
            continue;
        };
        let Some(VertexAttributeValues::Float32x3(positions)) =
            mesh.attribute(Mesh::ATTRIBUTE_POSITION)
        else {
            continue;
        };
        let Some(Indices::U32(indices)) = mesh.indices() else {
            continue;
        };
        for tri in indices.chunks_exact(3) {
            let p: [Vec3; 3] = [
                transform.transform_point(Vec3::from(positions[tri[0] as usize])),
                transform.transform_point(Vec3::from(positions[tri[1] as usize])),
                transform.transform_point(Vec3::from(positions[tri[2] as usize])),
            ];
            gizmos.line(p[0], p[1], colour);
            gizmos.line(p[1], p[2], colour);
            gizmos.line(p[2], p[0], colour);
        }
    }
}

/// One short line per vertex, along its normal.
///
/// Capped, because a dense mesh would otherwise submit hundreds of thousands of
/// gizmo lines and the frame-time readout would be measuring this rather than
/// the thing under test.
fn draw_normals(
    flags: Res<ViewFlags>,
    meshes: Res<Assets<Mesh>>,
    query: Query<(&Mesh3d, &GlobalTransform), With<DemoMesh>>,
    mut gizmos: Gizmos,
) {
    if !flags.normals {
        return;
    }
    const MAX_LINES: usize = 20_000;
    let colour = Color::srgb(0.98, 0.62, 0.16);
    for (handle, transform) in &query {
        let Some(mesh) = meshes.get(&handle.0) else {
            continue;
        };
        let (
            Some(VertexAttributeValues::Float32x3(positions)),
            Some(VertexAttributeValues::Float32x3(normals)),
        ) = (
            mesh.attribute(Mesh::ATTRIBUTE_POSITION),
            mesh.attribute(Mesh::ATTRIBUTE_NORMAL),
        )
        else {
            continue;
        };
        let stride = (positions.len() / MAX_LINES).max(1);
        for i in (0..positions.len()).step_by(stride) {
            let base = transform.transform_point(Vec3::from(positions[i]));
            let tip = base + transform.affine().matrix3 * Vec3::from(normals[i]) * 0.09;
            gizmos.line(base, tip, colour);
        }
    }
}

fn draw_domain(flags: Res<ViewFlags>, query: Query<&DemoDomain>, mut gizmos: Gizmos) {
    if !flags.grid {
        return;
    }
    let colour = Color::srgb(0.35, 0.38, 0.48);
    for domain in &query {
        let (lo, hi) = (domain.min, domain.max);
        // The eight corners, indexed by the same xyz bit pattern the extractor
        // uses for cube corners.
        let corner = |i: usize| {
            Vec3::new(
                if i & 1 == 0 { lo.x } else { hi.x },
                if i & 2 == 0 { lo.y } else { hi.y },
                if i & 4 == 0 { lo.z } else { hi.z },
            )
        };
        for i in 0..8usize {
            for axis in 0..3usize {
                let bit = 1 << axis;
                if i & bit == 0 {
                    gizmos.line(corner(i), corner(i | bit), colour);
                }
            }
        }
    }
}

/// A material that shows shape rather than texture — the thing every one of
/// these examples is actually about.
#[allow(dead_code)] // Each example compiles its own copy of this module.
pub fn surface_material(materials: &mut Assets<StandardMaterial>) -> Handle<StandardMaterial> {
    materials.add(StandardMaterial {
        base_color: Color::srgb(0.72, 0.76, 0.82),
        perceptual_roughness: 0.45,
        metallic: 0.05,
        ..default()
    })
}

/// Marker so `cargo check --all-targets` compiles this module even though it is
/// only ever reached through `mod common;` in an example.
#[allow(dead_code)]
pub const RENDER_ASSET_USAGES_NOTE: RenderAssetUsages = RenderAssetUsages::MAIN_WORLD;
