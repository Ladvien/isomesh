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
//! | `H` | HUD |
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
            .add_systems(Update, size_window)
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
    /// Whether hiding the HUD leaves a line saying how to bring it back.
    ///
    /// `ISOMESH_VIEW=nohud` clears this as well as [`Self::hud`], because a
    /// media capture wants an empty frame -- `scripts/record_all_gifs.sh` records
    /// nine demos that way. A reader who presses `H` gets the hint, because a
    /// panel that vanishes with no way back is a bug that looks like a crash.
    pub hud_hint: bool,
    /// Whether `ISOMESH_VIEW` asked for the HUD *positively*.
    ///
    /// `nohud` turns it off; this is the other direction, and it exists because
    /// a demo may reasonably open with the panel hidden -- `game_dig` does, it
    /// is a game first and the panel covers the rock. Such a demo clears
    /// [`Self::hud`] in its own `setup` **unless** this is set, which is what
    /// keeps `ISOMESH_SCREENSHOT` able to capture the numbers: the committed
    /// still is taken with `ISOMESH_VIEW=hud`.
    // Only `game_dig` reads this, and each example compiles its own copy of this
    // module -- the same reason the free functions here carry the attribute.
    #[allow(dead_code)]
    pub hud_requested: bool,
}

impl ViewFlags {
    /// The flags an `ISOMESH_VIEW` list and an `ISOMESH_FIELD` index ask for.
    ///
    /// Split out of [`Default`] so it can be tested. `std::env::set_var` is
    /// `unsafe` and this crate's `[lints.rust]` says `unsafe_code = "forbid"`,
    /// so a test cannot set the variable -- and a test that could would race
    /// every other test in the process, because the environment is global.
    fn parse(view: &str, field: usize) -> Self {
        let has = |name: &str| view.split(',').any(|part| part.trim() == name);
        Self {
            wireframe: has("wire"),
            normals: has("normals"),
            grid: !has("nogrid"),
            paused: false,
            remesh_requested: false,
            field,
            hud: !has("nohud"),
            hud_hint: !has("nohud"),
            hud_requested: has("hud"),
        }
    }
}

impl Default for ViewFlags {
    /// Toggles can be preset from `ISOMESH_VIEW`, a comma-separated list of
    /// `wire`, `normals`, `nogrid`, `nohud`.
    ///
    /// That exists so a screenshot of a particular view can be captured without
    /// a human pressing a key, which is what makes the visual acceptance tests
    /// in the examples catalog reproducible rather than anecdotal.
    fn default() -> Self {
        Self::parse(
            &std::env::var("ISOMESH_VIEW").unwrap_or_default(),
            std::env::var("ISOMESH_FIELD")
                .ok()
                .and_then(|v| v.parse().ok())
                .unwrap_or(0),
        )
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
    /// A headline above the panel, in its own colour and a larger size.
    ///
    /// One line, and it exists because the thing a reader most wants to see is
    /// usually not a number -- in `game_dig` it is *which mesher is running*,
    /// and that was eight lines down a monochrome panel. `None` leaves the row
    /// empty and the panel where it has always been, so a demo that sets nothing
    /// looks exactly as it did.
    ///
    /// A second `Text` entity rather than a coloured span inside the panel:
    /// `update_hud` assembles one string, and turning that into spans would
    /// rewrite the text path every demo shares for the benefit of one line.
    pub banner: Option<(String, Color)>,
    /// The one line left on screen when the HUD is hidden.
    ///
    /// `None` prints `[H] HUD`, which is the minimum a reader needs to undo a
    /// keypress. A demo that opens *with* the panel hidden wants more than that
    /// -- its whole key list is otherwise invisible -- so it replaces this.
    /// Suppressed entirely by `ISOMESH_VIEW=nohud`, which wants an empty frame.
    pub hint: Option<String>,
    /// The key list, when an example's bindings are not the harness's.
    ///
    /// `None` prints the shared footer below. `Some` replaces it outright rather
    /// than appending, because the point is that the shared line is *wrong* for
    /// that example -- `game_dig` walks on `WASD`, so a footer offering `[W]
    /// wire` documents a key that walks forward.
    pub keys: Option<String>,
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
    // Created here rather than assumed. A missing directory used to produce one
    // "Cannot save screenshot, IO error" per frame and no files, which reads as
    // a broken capture rather than a missing `mkdir` -- and the errors come from
    // deep inside Bevy's screenshot observer, nowhere near the cause.
    if capture.elapsed == 0
        && let Err(e) = std::fs::create_dir_all(&dir)
    {
        error!("ISOMESH_CAPTURE={dir} cannot be created: {e}");
        capture.dir = None;
        return;
    }
    capture.elapsed += 1;
    if capture.elapsed <= capture.settle {
        return;
    }
    if !(capture.elapsed - capture.settle).is_multiple_of(capture.every) {
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

/// How many frames [`size_window`] keeps re-applying the requested size.
///
/// **Re-applied rather than set once, because neither "set it early" nor "set it
/// and check" works.** At `PreStartup` the OS window does not exist yet, so the
/// write lands on an entity and is then overwritten by whatever the window comes
/// back as. And `Window::resolution` reads back the value *this* system wrote
/// rather than what the platform granted, so a system that stops when the two
/// agree stops on its own echo — measured, and it reproduced the same wrong size
/// as the `PreStartup` version.
///
/// So it re-applies across the frames in which the window is created and
/// configured, then stops. This must stay comfortably below
/// `ISOMESH_CAPTURE_SETTLE` (45 by default) so a sequence's first frame is taken
/// after the size has stopped moving.
const SIZE_WINDOW_FRAMES: u32 = 30;

/// Force the window to a given size, for reproducible captures.
///
/// `ISOMESH_WINDOW=1280x720`. Without it the window takes whatever it is given —
/// on a tiling compositor that is whatever slot happens to be free, and with no
/// window manager at all it is whatever the X server picks, so two captures of
/// the same example come back different shapes and cannot be composited side by
/// side.
///
/// # It ran in `PreStartup` and silently did nothing (E-214, M-235 amended)
///
/// The OS window does not exist yet at `PreStartup` — `bevy_winit` creates it
/// later — so writing `Window::resolution` there set a field on an entity and
/// then had it overwritten by whatever the window came back as. Measured: with
/// no window manager, `1280x720` and `1600x900` **both** produced 836×1356, with
/// no error, which is precisely the guarantee this function claims to provide.
///
/// It runs in `Update` instead, once, after the window exists. A window can
/// always resize itself; that never needed a compositor, which is what the first
/// diagnosis got wrong. `Local<bool>` rather than a resource because the state is
/// this system's alone, and re-applying every frame would fight a human dragging
/// the window edge.
fn size_window(
    mut frames: Local<u32>,
    mut windows: Query<&mut Window, With<bevy::window::PrimaryWindow>>,
) {
    if *frames > SIZE_WINDOW_FRAMES {
        return;
    }
    *frames += 1;
    let Ok(spec) = std::env::var("ISOMESH_WINDOW") else {
        *frames = SIZE_WINDOW_FRAMES + 1;
        return;
    };
    let Some((w, h)) = spec.split_once(['x', 'X']) else {
        error!("ISOMESH_WINDOW={spec} is not WIDTHxHEIGHT");
        return;
    };
    let (Ok(w), Ok(h)) = (w.trim().parse::<f32>(), h.trim().parse::<f32>()) else {
        error!("ISOMESH_WINDOW={spec} is not WIDTHxHEIGHT");
        return;
    };
    for mut window in &mut windows {
        window.resolution.set(w, h);
    }
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

/// The panel: title, numbers, extra lines, keys.
#[derive(Component)]
struct HudText;

/// [`DemoStats::banner`]'s own line, above the panel and in its own colour.
#[derive(Component)]
struct HudBanner;

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

/// Two entities, because the banner carries its own colour and size.
///
/// The panel keeps `top: 10` when there is no banner, so every demo that sets
/// none is pixel-identical to before; `update_hud` pushes it down to 32 only
/// when a banner is present.
fn spawn_hud(mut commands: Commands) {
    commands.spawn((
        Text::new(""),
        TextFont {
            // Larger than the panel, which is what "headline" means when the
            // only font available is the default one -- there is no bold face to
            // ask for, so size and colour do the work.
            font_size: FontSize::Px(19.0),
            ..default()
        },
        TextColor(Color::WHITE),
        Node {
            position_type: PositionType::Absolute,
            top: Val::Px(8.0),
            left: Val::Px(12.0),
            ..default()
        },
        HudBanner,
    ));
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
    // Names F12 screenshots. A counter rather than `std::process::id()`, which
    // compiles on wasm and then panics with "no pids on this platform" -- and
    // `Window::prevent_default_event_handling` defaults to `true`, so F12 in a
    // browser reaches the app instead of devtools. `Local<u32>` is the same idiom
    // `size_window` uses for its frame counter.
    mut taken: Local<u32>,
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
    if keys.just_pressed(KeyCode::KeyH) {
        flags.hud = !flags.hud;
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
        *taken += 1;
        let path = format!("screenshot-{}.png", *taken);
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

#[allow(clippy::type_complexity)]
fn update_hud(
    stats: Res<DemoStats>,
    flags: Res<ViewFlags>,
    times: Res<FrameTimes>,
    mut query: Query<(&mut Text, &mut Node), With<HudText>>,
    mut banner: Query<(&mut Text, &mut TextColor), (With<HudBanner>, Without<HudText>)>,
) {
    // The banner outlives a hidden panel: in `game_dig` it names the mesher, and
    // that is exactly the thing worth seeing with the numbers switched off. It
    // goes with `nohud` though, which wants an empty frame.
    let (headline, tint) = match (&stats.banner, flags.hud_hint) {
        (Some((line, colour)), true) => (line.as_str(), *colour),
        _ => ("", Color::WHITE),
    };
    for (mut text, mut colour) in &mut banner {
        // Written only when it differs, the same change-driven-write rule
        // `active_cells.rs` states for its `BackgroundColor`s: an unconditional
        // write marks the text changed every frame, and Bevy's UI extraction is
        // change-driven, so a static line would become per-frame work.
        if text.0 != headline {
            text.0 = headline.to_string();
        }
        if colour.0 != tint {
            colour.0 = tint;
        }
    }
    // 32 clears the 19 px headline plus its 8 px inset; 10 is where the panel has
    // always sat, and a demo that sets no banner keeps it.
    let top = Val::Px(if headline.is_empty() { 10.0 } else { 32.0 });

    let frame_ms = median(&times.0);
    let fps = if frame_ms > 0.0 {
        1000.0 / frame_ms
    } else {
        0.0
    };

    if !flags.hud {
        // **The frame rate outlives the panel too, and for the same reason the
        // banner does**: with the numbers switched off the two things worth
        // seeing are which mesher is running and what it costs, and the cost was
        // the one of the pair that vanished. It joins the hint rather than the
        // banner because the banner already carries the mesher — one line each,
        // no duplication when the panel is open.
        //
        // Still suppressed by `nohud`, which wants an *empty* frame: the hint is
        // `""` there and this appends nothing to it, so the committed captures
        // are unchanged.
        let hint = match (&stats.hint, flags.hud_hint) {
            (_, false) => String::new(),
            (Some(line), true) => format!("{fps:.0} fps   {line}"),
            (None, true) => format!("{fps:.0} fps   [H] HUD"),
        };
        for (mut hud, mut node) in &mut query {
            if hud.0 != hint {
                hud.0.clone_from(&hint);
            }
            if node.top != top {
                node.top = top;
            }
        }
        return;
    }

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
    match &stats.keys {
        Some(keys) => {
            text.push_str("\n\n");
            text.push_str(keys);
        }
        None => text.push_str(&format!(
            "\n\n[W] wire {}   [N] normals {}   [G] grid {}\n[Space] {}   [R] re-mesh   [H] HUD   [F12] shot   [Esc] quit",
            on_off(flags.wireframe),
            on_off(flags.normals),
            on_off(flags.grid),
            if flags.paused { "resume" } else { "pause" },
        )),
    }

    for (mut target, mut node) in &mut query {
        if target.0 != text {
            target.0.clone_from(&text);
        }
        if node.top != top {
            node.top = top;
        }
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
        for tri in indices.as_chunks::<3>().0 {
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

#[cfg(test)]
mod tests {
    use std::collections::VecDeque;

    use bevy::ecs::system::RunSystemOnce;

    use super::*;

    /// `nohud` clears the hint as well as the HUD, and an empty list leaves both
    /// on.
    ///
    /// This is the gate on the nine GIFs `scripts/record_all_gifs.sh` records
    /// with `ISOMESH_VIEW=nohud`: they must be *empty* frames, so the reader's
    /// `[H] HUD` hint must not appear in them. A browser has no environment, so
    /// this half of the flag cannot be checked from a screenshot.
    #[test]
    fn nohud_hides_the_hint_too() {
        let media = ViewFlags::parse("nohud", 0);
        assert!(!media.hud, "nohud must hide the HUD");
        assert!(!media.hud_hint, "nohud must hide the hint as well");

        let reader = ViewFlags::parse("", 0);
        assert!(reader.hud, "the HUD is on by default");
        assert!(reader.hud_hint, "pressing H must leave a way back");

        // `nohud` in a longer list, which is how `record_all_gifs.sh` writes it
        // for the demos that also want a wireframe.
        let mixed = ViewFlags::parse("wire, nohud", 0);
        assert!(mixed.wireframe);
        assert!(!mixed.hud && !mixed.hud_hint);
    }

    /// `hud` is the positive request, and only `hud` sets it.
    ///
    /// `game_dig` opens with the panel hidden and clears `hud` in its own
    /// `setup` *unless* this is set, so this flag is the only thing that lets
    /// `ISOMESH_SCREENSHOT` still capture the numbers. A default that answered
    /// `true` here would make that demo's committed still un-retakeable.
    #[test]
    fn hud_is_the_positive_request() {
        assert!(ViewFlags::parse("hud", 0).hud_requested);
        assert!(ViewFlags::parse("wire, hud", 0).hud_requested);
        assert!(!ViewFlags::parse("", 0).hud_requested);
        assert!(!ViewFlags::parse("nohud", 0).hud_requested);
        // Not a prefix match: `nohud` contains `hud` as a substring and must not
        // be read as asking for it.
        assert!(!ViewFlags::parse("nogrid,nohud", 0).hud_requested);
    }

    /// With the panel hidden, the frame rate and the mesher are both still on
    /// screen — and with `nohud` neither is.
    ///
    /// Two things a reader needs while playing rather than while reading:
    /// **which mesher is running** and **what it costs**. The banner carried the
    /// first through a hidden panel from the start; the frame rate went with the
    /// panel, which is the gap this closes. `game_dig` is the demo that opens
    /// hidden, so this was its whole steady state.
    ///
    /// Driven through `update_hud` itself rather than by inspecting the strings
    /// it is built from: the `!flags.hud` early return is the branch that used to
    /// drop the number, and a test that formats its own line would not cross it.
    /// `nohud` is in the same test because the two requirements pull opposite
    /// ways — always visible, and an *empty* frame for the GIFs — and a fix for
    /// one is the obvious way to break the other.
    fn hud_lines(view: &str, hud: bool) -> (String, String) {
        let mut app = App::new();
        let mut flags = ViewFlags::parse(view, 0);
        flags.hud = hud;
        app.insert_resource(flags)
            .insert_resource(FrameTimes(VecDeque::from([8.0, 8.0, 8.0])))
            .insert_resource(DemoStats {
                banner: Some(("[1] Marching Cubes".to_string(), Color::WHITE)),
                hint: Some("[H] numbers".to_string()),
                ..default()
            });
        app.world_mut()
            .run_system_once(spawn_hud)
            .expect("spawn_hud needs only `Commands`");
        app.world_mut()
            .run_system_once(update_hud)
            .expect("update_hud needs no window and no renderer");
        let mut panel = app
            .world_mut()
            .query_filtered::<&Text, (With<HudText>, Without<HudBanner>)>();
        let mut banner = app.world_mut().query_filtered::<&Text, With<HudBanner>>();
        let panel = panel.iter(app.world()).next().expect("the panel").0.clone();
        let banner = banner
            .iter(app.world())
            .next()
            .expect("the banner")
            .0
            .clone();
        (panel, banner)
    }

    #[test]
    fn a_hidden_panel_still_shows_the_frame_rate_and_the_mesher() {
        // 8 ms a frame is 125 fps, so the number is checked rather than merely
        // the word.
        let (line, banner) = hud_lines("", false);
        assert!(
            line.starts_with("125 fps"),
            "the frame rate is not on the line the hidden panel leaves: {line:?}"
        );
        assert!(
            line.contains("[H] numbers"),
            "the frame rate displaced the hint instead of joining it: {line:?}"
        );
        assert_eq!(
            banner, "[1] Marching Cubes",
            "the mesher is not on screen with the panel hidden"
        );

        // Open, the number belongs to the panel and must not be on both.
        let (panel, banner) = hud_lines("", true);
        assert!(
            panel.contains("125 fps"),
            "the open panel lost the frame rate: {panel:?}"
        );
        assert_eq!(banner, "[1] Marching Cubes");

        // `nohud` wants an empty frame, and that is what the GIFs are recorded
        // with. Both lines, because the frame rate is new and the banner is not.
        let (line, banner) = hud_lines("nohud", false);
        assert_eq!(line, "", "nohud left a frame rate on an empty frame");
        assert_eq!(banner, "", "nohud left the banner on an empty frame");
    }
}
