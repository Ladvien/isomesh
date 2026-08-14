//! E-303 — the whole thing on the GPU: field, extraction and draw.
//!
//! ```bash
//! cd bevy_isomesh && cargo run --example gpu_mesh_shader --release
//! ```
//!
//! `[` `]` resolution, `Space` pause the moving brush. Drag to orbit.
//!
//! # What is actually happening each frame
//!
//! Nothing crosses the bus except a camera matrix and a handful of brushes:
//!
//! 1. **The field is evaluated** by a compute pass (GPU-011a) — a base sphere
//!    with a moving edit log folded over it (GPU-011b). No samples are uploaded.
//! 2. **Marching Cubes runs** and the per-cell prefix sum is scanned on the GPU
//!    (GPU-010a). Four bytes come home: the triangle count, to size the buffer.
//! 3. **A mesh shader draws it** straight out of the position and normal buffers
//!    the extraction wrote (GPU-008a). There is no vertex buffer, no index
//!    buffer, and the geometry is never read back.
//!
//! Measured for the extraction half at 129³: **0.54 ms** (M-155), against 15.01
//! before GPU-010a. The read-back this last step removes is `~0.63 ms` of the
//! `~1.17 ms` a consumer would otherwise pay to get the mesh onto the CPU —
//! about half.
//!
//! # It refuses rather than falling back
//!
//! On an adapter without `EXPERIMENTAL_MESH_SHADER` the pipeline is not built,
//! the HUD says so, and nothing is drawn. That is one path chosen by a
//! measurement with the choice visible — not a library quietly substituting a
//! vertex-buffer pipeline, which is what the one-path rule forbids. Bevy's
//! default `Functionality` priority requests every feature the adapter has
//! (M-147), so on capable hardware this simply works.
//!
//! # Where the drawing happens
//!
//! Bevy 0.19 has no `Node` trait: `Core3d` is a **schedule** with
//! `Prepass → MainPass → EarlyPostProcess → PostProcess` sets, and its own docs
//! invite additional systems into them. So this is a system in
//! `Core3dSystems::MainPass`, after the opaque pass has cleared the target, with
//! `ViewQuery` for the view's components and `RenderContext` for the encoder.

mod common;

use std::sync::{Arc, Mutex};

use bevy::core_pipeline::core_3d::CORE_3D_DEPTH_FORMAT;
use bevy::core_pipeline::schedule::{Core3d, Core3dSystems};
use bevy::prelude::*;
use bevy::render::RenderApp;
use bevy::render::extract_resource::{ExtractResource, ExtractResourcePlugin};
use bevy::render::render_resource::{
    CompareFunction, DepthBiasState, DepthStencilState, StencilState, StoreOp,
};
use bevy::render::renderer::{RenderContext, RenderDevice, RenderQueue, ViewQuery};
use bevy::render::view::{ExtractedView, ViewDepthTexture, ViewTarget};
use common::{CommonPlugin, DemoStats, OrbitCamera};
use isomesh_gpu::{
    FieldSampler, GpuBrush, GpuField, GpuOp, GpuShape, GridParams, MarchingCubesGpu,
    MeshShaderRenderer,
};

/// Samples per axis, unless `ISOMESH_SAMPLES` says otherwise.
const DEFAULT_SAMPLES: u32 = 65;

/// What the render world needs to know, extracted from the main world each
/// frame. Kilobytes, against the megabytes an uploaded field would be.
#[derive(Resource, Clone)]
struct Scene {
    samples: u32,
    /// Drives the moving brushes. Frozen by `Space`.
    time: f32,
    /// Read the triangle count back this frame. Off by default, because not
    /// reading it is the property this demo exists to show.
    probe: bool,
}

impl ExtractResource for Scene {
    type Source = Self;
    fn extract_resource(source: &Self) -> Self {
        source.clone()
    }
}

/// Filled in by the render world, read by the HUD.
///
/// An `Arc` rather than a channel: the render world writes and the main world
/// reads, once per frame, and a lock held for a struct copy is cheaper than the
/// plumbing a channel would need in both directions.
#[derive(Resource, Clone, Default)]
struct Readout(Arc<Mutex<Stats>>);

impl ExtractResource for Readout {
    type Source = Self;
    fn extract_resource(source: &Self) -> Self {
        source.clone()
    }
}

#[derive(Default, Clone, Copy)]
struct Stats {
    triangles: u32,
    budget: u32,
    extract_ms: f64,
    mesh_shaders: bool,
    ran: bool,
    /// Whether `triangles` came from an actual read-back this frame.
    probed: bool,
}

/// Pipelines, built once the render device exists.
#[derive(Resource)]
struct Kit {
    sampler: FieldSampler,
    marching_cubes: MarchingCubesGpu,
    /// `None` until the first draw, because the colour format is a property of
    /// the view rather than of the device.
    renderer: Option<MeshShaderRenderer>,
    /// Whether the device can run a mesh pipeline at all. Checked once.
    supported: bool,
}

fn main() {
    App::new()
        .add_plugins(DefaultPlugins.set(WindowPlugin {
            primary_window: Some(Window {
                title: "isomesh - E-303 gpu mesh shader".into(),
                ..default()
            }),
            ..default()
        }))
        .add_plugins(CommonPlugin)
        .add_plugins(ExtractResourcePlugin::<Scene>::default())
        .add_plugins(ExtractResourcePlugin::<Readout>::default())
        .insert_resource(Scene {
            samples: common::samples_override().unwrap_or(DEFAULT_SAMPLES),
            time: 0.0,
            probe: false,
        })
        .init_resource::<Readout>()
        .add_systems(Startup, setup)
        .add_systems(Update, (drive, report))
        .add_plugins(DrawPlugin)
        .run();
}

fn setup(mut commands: Commands, camera: Query<Entity, With<OrbitCamera>>) {
    // Keep the shared harness's orbit camera; it already frames a [-2, 2] box.
    for entity in &camera {
        commands
            .entity(entity)
            .insert(Transform::from_xyz(0.0, 1.2, 5.0));
    }
}

fn drive(
    time: Res<Time>,
    keys: Res<ButtonInput<KeyCode>>,
    flags: Res<common::ViewFlags>,
    mut scene: ResMut<Scene>,
) {
    if !flags.paused {
        scene.time += time.delta_secs();
    }
    if keys.just_pressed(KeyCode::BracketLeft) {
        scene.samples = (scene.samples - 16).max(17);
    }
    if keys.just_pressed(KeyCode::BracketRight) {
        scene.samples = (scene.samples + 16).min(129);
    }
    // Held, not toggled: the cost is per frame it is on.
    scene.probe = keys.pressed(KeyCode::KeyT);
}

fn report(scene: Res<Scene>, readout: Res<Readout>, mut stats: ResMut<DemoStats>) {
    let snapshot = readout.0.lock().map(|s| *s).unwrap_or_default();
    stats.title = format!("E-303 gpu mesh shader - {}^3", scene.samples);
    stats.triangles = snapshot.triangles as usize;
    stats.vertices = snapshot.triangles as usize * 3;
    stats.extract_ms = snapshot.extract_ms;
    stats.extra = vec![
        if snapshot.mesh_shaders {
            "drawn by a MESH SHADER, straight from the compute output".to_string()
        } else {
            "mesh shaders UNAVAILABLE on this device -- nothing is drawn, and that is\n  the whole behaviour: no vertex-buffer pipeline is quietly substituted.\n  The compute path has its own demos, gpu_compute_mc and gpu_vs_cpu."
                .to_string()
        },
        String::new(),
        "field evaluated on the gpu    (GPU-011a)".to_string(),
        "brushes folded on the gpu     (GPU-011b)".to_string(),
        "prefix sum scanned on the gpu (GPU-010a)".to_string(),
        "geometry never read back      (GPU-008b)".to_string(),
        String::new(),
        format!(
            "triangle budget {:>9}   {}",
            snapshot.budget,
            if snapshot.probed {
                format!(
                    "actual {} -- read back because you held T",
                    snapshot.triangles
                )
            } else {
                "actual: NOT READ BACK. hold [T] to pay for it".to_string()
            }
        ),
        String::new(),
        "per frame the cpu sends a camera matrix and three brushes, and waits".to_string(),
        "for nothing. the triangle count goes straight into the draw's own".to_string(),
        "arguments (GPU-010b) -- so the cpu cannot know it without asking.".to_string(),
        String::new(),
        "[ ] resolution   [Space] pause the brushes   [T] read the count".to_string(),
    ];
}

/// Adds the draw system to the render world.
struct DrawPlugin;

impl Plugin for DrawPlugin {
    fn build(&self, app: &mut App) {
        let Some(render_app) = app.get_sub_app_mut(RenderApp) else {
            return;
        };
        // `Core3d` is a schedule in Bevy 0.19, and its docs invite extra
        // systems into these sets. After the opaque pass, so the target is
        // already cleared and the depth buffer is live.
        // `EarlyPostProcess` rather than `MainPass`: the sets are chained, so
        // this is guaranteed to run after the opaque pass has cleared the
        // target, without naming one of Bevy's systems and depending on it
        // staying public.
        render_app.add_systems(
            Core3d,
            draw_isosurface.in_set(Core3dSystems::EarlyPostProcess),
        );
    }

    fn finish(&self, app: &mut App) {
        let Some(render_app) = app.get_sub_app_mut(RenderApp) else {
            return;
        };
        let device = render_app.world().resource::<RenderDevice>().clone();
        let queue = render_app.world().resource::<RenderQueue>().clone();
        let kit = Kit {
            sampler: FieldSampler::new(device.wgpu_device()).expect("this crate's own shader"),
            marching_cubes: MarchingCubesGpu::new(device.wgpu_device(), &queue)
                .expect("this crate's own shader"),
            renderer: None,
            // `ISOMESH_NO_MESH_SHADER=1` forces the unsupported branch. Without
            // it that branch is unreachable on any machine that can run this
            // demo, and a branch nothing ever executes is the failure mode this
            // project keeps finding — so it gets a switch and a screenshot
            // rather than a hope.
            supported: MeshShaderRenderer::is_supported(device.wgpu_device())
                && std::env::var("ISOMESH_NO_MESH_SHADER").is_err(),
        };
        render_app.insert_resource(kit);
    }
}

/// The brushes, moving. Three is enough to see the log being folded live.
fn brushes(time: f32) -> [GpuBrush; 3] {
    [
        GpuBrush {
            shape: GpuShape::Sphere {
                center: [0.9 * (time * 0.7).cos(), 0.5 * (time * 0.9).sin(), 0.3],
                radius: 0.5,
            },
            op: GpuOp::Subtract,
        },
        GpuBrush {
            shape: GpuShape::Capsule {
                a: [-0.9, -0.2, 0.0],
                b: [0.9, 0.2 + 0.4 * (time * 0.5).sin(), 0.0],
                radius: 0.18,
            },
            op: GpuOp::SmoothAdd { k: 0.2 },
        },
        GpuBrush {
            shape: GpuShape::BoxExact {
                center: [0.0, -0.9, 0.0],
                half_extents: [1.2, 0.25, 1.2],
            },
            op: GpuOp::Add,
        },
    ]
}

/// Evaluate, extract and draw — the whole pipeline, once per frame.
fn draw_isosurface(
    view: ViewQuery<(&ExtractedView, &ViewTarget, &ViewDepthTexture)>,
    kit: Option<ResMut<Kit>>,
    scene: Option<Res<Scene>>,
    readout: Option<Res<Readout>>,
    render_device: Res<RenderDevice>,
    render_queue: Res<RenderQueue>,
    mut ctx: RenderContext,
) {
    let (Some(mut kit), Some(scene), Some(readout)) = (kit, scene, readout) else {
        return;
    };
    let (extracted, target, depth) = view.into_inner();
    let device = render_device.wgpu_device();
    let queue: &wgpu::Queue = &render_queue;

    let mut stats = Stats {
        mesh_shaders: kit.supported,
        ..Stats::default()
    };
    let publish = |stats: Stats| {
        if let Ok(mut slot) = readout.0.lock() {
            *slot = stats;
        }
    };
    if !kit.supported {
        publish(stats);
        return;
    }

    // Built on first draw: the colour format belongs to the view, not the
    // device, so it is not knowable in `finish`.
    if kit.renderer.is_none() {
        kit.renderer = MeshShaderRenderer::new(
            device,
            extracted.target_format,
            Some(DepthStencilState {
                format: CORE_3D_DEPTH_FORMAT,
                depth_write_enabled: Some(true),
                // Bevy 0.19 uses a reversed-Z projection, so nearer is greater.
                depth_compare: Some(CompareFunction::GreaterEqual),
                stencil: StencilState::default(),
                bias: DepthBiasState::default(),
            }),
            // The view's sample count, read from the depth texture rather than
            // assumed: Bevy defaults to 4x MSAA and a pipeline built for 1
            // fails the draw outright rather than looking worse.
            depth.texture.sample_count(),
        )
        .ok();
    }
    let Some(renderer) = kit.renderer.as_ref() else {
        publish(stats);
        return;
    };

    let started = std::time::Instant::now();
    let cell = 4.0 / (scene.samples - 1) as f32;
    let Ok(grid) = GridParams::new([scene.samples; 3], [-2.0; 3], cell) else {
        return;
    };

    // Field and extraction, both on the GPU. Nothing uploaded but the brushes.
    let Ok(field) =
        kit.sampler
            .sample_stack(device, queue, grid, GpuField::Sphere, &brushes(scene.time))
    else {
        return;
    };
    // Zero read-backs: the triangle count is written into the draw's own
    // arguments and never comes home. The budget is the price -- see
    // `extract_indirect` -- and a quarter of the cell count is generous against
    // the ~2% of cells a surface actually crosses.
    let budget = (grid.cell_count() / 4).clamp(4096, 4_000_000) as u32;
    let Ok(geometry) = kit
        .marching_cubes
        .extract_indirect(device, queue, &field, budget)
    else {
        return;
    };
    stats.extract_ms = started.elapsed().as_secs_f64() * 1000.0;
    stats.budget = budget;
    stats.ran = true;

    // The count is only knowable by reading it, which is the whole point -- so
    // it is read on request rather than every frame.
    if scene.probe
        && let Ok(total) = isomesh_gpu::read_buffer_u32(device, queue, &geometry.total, 4)
    {
        stats.triangles = total.first().copied().unwrap_or(0);
        stats.probed = true;
    }
    publish(stats);

    // The one thing the CPU still sends: where the camera is.
    let clip_from_world = extracted.clip_from_world.unwrap_or_else(|| {
        extracted.clip_from_view * extracted.world_from_view.to_matrix().inverse()
    });
    let uniform = |label, bytes: &[u8]| {
        let buffer = device.create_buffer(&wgpu::BufferDescriptor {
            label: Some(label),
            size: bytes.len() as u64,
            usage: wgpu::BufferUsages::UNIFORM | wgpu::BufferUsages::COPY_DST,
            mapped_at_creation: false,
        });
        queue.write_buffer(&buffer, 0, bytes);
        buffer
    };
    let camera = uniform(
        "isomesh camera",
        &MeshShaderRenderer::camera_bytes(clip_from_world.to_cols_array_2d()),
    );

    let bind_group = renderer.bind_group(
        device,
        &camera,
        // The draw's uniform, written by a kernel from the scanned total rather
        // than uploaded. Same sixteen bytes, filled on the other side.
        &geometry.draw_params,
        &geometry.positions,
        &geometry.normals,
    );

    // Load, not clear: the opaque pass already cleared and this draws over it.
    let mut colour = target.get_color_attachment();
    colour.ops.load = wgpu::LoadOp::Load;

    let encoder = ctx.command_encoder();
    let mut pass = encoder.begin_render_pass(&wgpu::RenderPassDescriptor {
        label: Some("isomesh mesh shader draw"),
        color_attachments: &[Some(colour)],
        depth_stencil_attachment: Some(depth.get_attachment(StoreOp::Store)),
        timestamp_writes: None,
        occlusion_query_set: None,
        multiview_mask: None,
    });
    renderer.draw_indirect(&mut pass, &bind_group, &geometry.indirect);
}
