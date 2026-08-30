//! **P-74 - ambient occlusion traced against the field the mesh came from.**
//!
//! Ticket: R-074. Pre-registered before this harness existed.
//!
//! ```bash
//! cargo bench --bench experiment_p74
//! ```
//!
//! Writes `docs/experiments/p-74.csv`.
//!
//! # C1's SHARE, recomputed before a line of this was written
//!
//! C1 is an absolute bar rather than a ratio, and its registered SHARE is *"the
//! whole AO budget, currently paid in full by a mesh-based pass"* - so there is
//! no `x51` denominator to be unreachable behind. What has to be checked is
//! whether **2.0 ms is arithmetically available at all**, and the two published
//! anchors the registration itself cites disagree:
//!
//! | anchor | implies at 1920x1080, 8 rays, 16 steps |
//! |---|---:|
//! | queries the configuration demands | 265,420,800 per frame |
//! | 2.0 ms at that count | **132.7 G queries/s** |
//! | nvblox's best measured rate, 7.3 G/s on an RTX 3090 Ti | 36.4 ms |
//! | nvblox's worst, 0.8 G/s | 331.8 ms |
//! | RTSDF's own 4.60 ms trace at 1024^2, scaled to 2.07 Mpx and to a 3090 | 6.06 ms |
//! | 8 trilinear loads per query, all resident in L1 (~17.8 TB/s aggregate) | 0.48 ms |
//! | the same traffic all the way out to VRAM (936 GB/s) | 9.1 ms |
//!
//! **The bar sits inside the spread of the memory hierarchy, and outside both
//! published anchors.** It is reachable only if three things hold at once: a
//! flat `129^3` `f32` grid is a far cheaper query than nvblox's block-hashed
//! TSDF lookup, the 8.587 MB grid stays in L2 (the RTX 3090 has 6 MB, so it does
//! not fit and has to be carried by locality), and sphere traces terminate early
//! rather than running all 16 steps. That is what makes it a measurement.
//!
//! # The three estimators all compute the same number, and that is the design
//!
//! AO here is **the cosine-weighted fraction of the hemisphere blocked within
//! `AO_RADIUS_CELLS` cells**, and all three arms estimate that one quantity so
//! that a mean absolute difference between them means something:
//!
//! - **reference** - `REF_RAYS` cosine-weighted Hammersley directions, cast
//!   against the extracted mesh with `parry3d`. Offline, on the CPU.
//! - **field** - the registered 8 rays x 16-step march, sphere-tracing the
//!   resident jump-flooded distance grid on the GPU. Sphere tracing's step *is*
//!   the cone (Hart, *Sphere tracing*, The Visual Computer 1996,
//!   `10.1007/s003710050084`): the marched radius is the distance to the nearest
//!   surface, so the swept volume widens exactly as a cone march's does.
//! - **SSAO** - 8 normal-oriented hemisphere samples against the rasterised
//!   depth buffer, with the standard range check (Mittring, *Finding next gen*,
//!   SIGGRAPH 2007 course, `10.1145/1281500.1281671`).
//!
//! Every arm starts its rays at `p + n * LIFT_CELLS` cells and uses the same
//! `onb` tangent frame (Duff et al., *Building an orthonormal basis, revisited*,
//! JCGT 2017) over the same direction table - **which is uploaded from the
//! Rust side rather than transcribed into WGSL**, so the field arm and the mesh
//! arms cannot sample different hemispheres.
//!
//! # The MAE floor, which is the control that makes C2 readable at all
//!
//! An 8-sample estimator of a Bernoulli integral differs from a 512-sample one
//! by discretisation alone, and that residue is common to *both* methods. So the
//! harness also computes `seam_mae_floor` / `halo_mae_floor`: the same 8
//! directions, cast against the **mesh** by `parry3d`, against the same 512-ray
//! reference. That is the error an *ideal* 8-ray method would post. A method's
//! MAE above the floor is its own; at the floor it is indistinguishable from
//! perfect at this ray count. Without this column, "field MAE 0.09" is
//! unreadable.
//!
//! # The vacuity control the registration names, and how it could have failed
//!
//! `seam_pixels` and `silhouette_pixels` must be non-empty counts. Both are
//! asserted, and both come with a **negative arm run by the same classifier on
//! the same image**, because a count that has only ever been positive is
//! indistinguishable from a classifier that counts pixels rather than reading
//! them (`M-44`):
//!
//! - `seam_pixels_one_chunk` reclassifies the identical G-buffer with the
//!   interior chunk-boundary plane list emptied - the world as one chunk. It
//!   must be exactly zero.
//! - `silhouette_pixels_flat` reclassifies a synthetic image in which every
//!   pixel hits at one constant view distance. It must be exactly zero.
//! - `seam_pixels_screen` is an independent second definition of the same set -
//!   a chunk-id discontinuity between screen neighbours - so the geometric band
//!   is cross-checked rather than trusted.
//!
//! # What the depth buffer is
//!
//! A real one. The per-chunk `MeshBuffer`s are concatenated into one vertex and
//! index buffer and **rasterised through a `wgpu` render pipeline** into a
//! `Depth32Float` depth attachment plus two `Rgba32Float` G-buffer targets
//! (world position + chunk id, world normal + hit mask). SSAO reads the depth
//! attachment and the normal target and nothing else; it reconstructs world
//! positions from depth through `inv_view_proj`. Nothing about the depth image
//! is synthesised: it is what a depth prepass over this crate's output produces.
//!
//! # Clocks
//!
//! `M-280`: a nanosecond is not a unit on a governed CPU. `ao_ms_field` and
//! `ao_ms_ssao` are **GPU timestamp spans** over the compute pass, median of
//! `REPS`, and the `clock` column says so. `resident_ms` and `rebuild_ms` are
//! wall clock around submit-and-wait, because they span a CPU-visible round
//! trip that no GPU span can see - which is the whole of C3.

#![allow(
    clippy::cast_possible_truncation,
    clippy::cast_precision_loss,
    clippy::cast_sign_loss,
    clippy::needless_range_loop,
    clippy::too_many_arguments,
    clippy::too_many_lines
)]

mod common;

use std::time::Instant;

use isomesh::chunk::{ChunkId, ChunkLayout};
use isomesh::extractor::Extractor;
use isomesh::fields::{FbmTerrain, Gyroid};
use isomesh::marching_cubes::MarchingCubes;
use isomesh::{MeshBuffer, Sdf};
use isomesh_gpu::wgpu;
use isomesh_gpu::{Composer, FieldBuffer, GridParams, headless, jump_flood::JumpFlood, read_bytes};
use parry3d::math::Vector;
use parry3d::query::{Ray, RayCast};
use parry3d::shape::TriMesh;

/// The registered resolution.
const WIDTH: u32 = 1920;
/// The registered resolution.
const HEIGHT: u32 = 1080;
/// The registered ray count per pixel.
const RAYS: u32 = 8;
/// The registered march length.
const STEPS: u32 = 16;

/// Cells per axis in the world. `p72`'s world exactly, so this row is comparable
/// with `M-377` and with `M-155`'s `129^3` GPU field.
const WORLD_CELLS: u32 = 128;
/// Cells per chunk edge, so 4^3 = 64 chunks and 9 interior boundary planes.
const CHUNK_CELLS: u32 = 32;
/// Samples per axis in the resident grid.
const SAMPLES: u32 = WORLD_CELLS + 1;
/// World extent per axis.
const EXTENT: f32 = 4.0;
/// World origin per axis, centred on the reference fields' own domain centre.
const ORIGIN: f32 = -EXTENT * 0.5;
/// Cell size. A power of two, so `origin + h*i` is exact and no seam arises from
/// the arithmetic itself (`M-70`, `M-73`).
const CELL: f32 = EXTENT / WORLD_CELLS as f32;

/// AO radius, in cells - one world unit at `CELL`.
///
/// **The first run used 8 cells and the 512-ray reference came back 0.0000 on
/// `gyroid` over 1536 seam pixels and 1536 silhouette pixels.** That is not a
/// result, it is `M-44`: `seam_mae_field` was zero because *both* sides were
/// zero, and a fixture in which the ideal method and the broken method score
/// identically cannot refute anything. The reason is geometric - `Gyroid`'s
/// channels at `scale = 1` are about `pi/2` wide, so a hemisphere of radius
/// 0.25 sees no opposite wall and the true occlusion really is nil. 32 cells is
/// 1.0 world unit, comparable to the field's own feature size, and
/// `assert_reference_non_vacuous` below now makes that a gate rather than a
/// thing to notice afterwards.
const AO_RADIUS_CELLS: f32 = 32.0;
/// Surface epsilon for the sphere trace, in cells.
///
/// **Not `AO_RADIUS_CELLS / STEPS` any more.** Tying them made the epsilon 2
/// cells at the wider radius, which stops the march on a surface it has not
/// reached and reports occlusion that is not there. Decoupled, the 16-step
/// budget is a genuine truncation: a sphere trace steps by the *true* distance
/// and covers far more than `16 * EPS` in open space, and any ray that runs out
/// of steps is counted unoccluded and tallied in `march_exhausted`. That number
/// is the honest cost of the registered 16 steps at this radius.
const EPS_CELLS: f32 = 0.5;
/// How far every arm lifts its ray origin off the surface, in cells.
const LIFT_CELLS: f32 = 1.0;
/// SSAO's depth-comparison bias, in cells.
const SSAO_BIAS_CELLS: f32 = 0.5;

/// How close to a chunk boundary plane a primary hit must be to count as a seam
/// pixel, in cells.
const SEAM_BAND_CELLS: f32 = 0.75;
/// How large a view-distance step between screen neighbours counts as a
/// silhouette, in cells.
const SILHOUETTE_STEP_CELLS: f32 = 4.0;

/// Rays in the offline reference.
const REF_RAYS: u32 = 512;
/// Pixels drawn from each set for the reference, spread by a uniform stride.
const SAMPLE_CAP: usize = 1536;

/// Timed repeats per arm, median taken.
const REPS: usize = 7;

/// Field of view, vertical, in radians. Chosen with `EYE_DISTANCE` so the world
/// box subtends 24.4 degrees against a 30-degree half-angle: the surface fills
/// most of the frame and there is background at the edges for silhouettes.
const FOV_Y: f32 = core::f32::consts::PI / 3.0;
/// Eye distance from the world centre, in world units.
const EYE_DISTANCE: f32 = EXTENT * 1.1;
/// Near plane.
const NEAR: f32 = 0.05;
/// Far plane.
const FAR: f32 = 20.0;

// ─── the shader ──────────────────────────────────────────────────────────────

/// The G-buffer pipeline and both AO kernels, in one module.
///
/// One module so the A/B is one submission shape and one set of constants, and
/// so `onb` has exactly one definition on the GPU side. `#include <grid>` pulls
/// in `GridParams` from the shipped crate rather than restating its layout - the
/// byte agreement with `GridParams::to_std140` is then the crate's own problem
/// and not a second copy in a bench.
const P74_WGSL: &str = r#"
#include <grid>

// A distance the march reads as "nothing here": outside the resident grid there
// is no geometry, and a value this large ends the march as a miss.
const OUTSIDE: f32 = 1.0e9;

struct Camera {
    view_proj: mat4x4<f32>,
    inv_view_proj: mat4x4<f32>,
    // xyz: eye position. w: unused.
    eye: vec4<f32>,
    // x: width, y: height, z: near, w: far.
    screen: vec4<f32>,
}

struct AoParams {
    // x: width, y: height, z: rays per pixel, w: march steps.
    dims: vec4<u32>,
    // x: ao radius, y: surface epsilon, z: normal lift, w: ssao bias.
    tune: vec4<f32>,
}

@group(0) @binding(0) var<uniform> camera: Camera;
@group(0) @binding(1) var<uniform> ao: AoParams;
@group(0) @binding(2) var<uniform> grid: GridParams;
@group(0) @binding(3) var<storage, read> distances: array<f32>;
// xyz: a cosine-weighted direction in the tangent frame. w: SSAO's radius scale.
@group(0) @binding(4) var<storage, read> dirs: array<vec4<f32>>;
@group(0) @binding(5) var gbuf_pos: texture_2d<f32>;
@group(0) @binding(6) var gbuf_nrm: texture_2d<f32>;
@group(0) @binding(7) var gbuf_depth: texture_depth_2d;
@group(0) @binding(8) var<storage, read_write> ao_out: array<f32>;
// 0: starts inside the field. 1: starts within eps of it. 2: primary hits.
// 3: marches that used all STEPS without hitting and without passing the radius.
@group(0) @binding(9) var<storage, read_write> diag: array<atomic<u32>>;

// ── the G-buffer ────────────────────────────────────────────────────────────

struct VsOut {
    @builtin(position) clip: vec4<f32>,
    @location(0) world: vec3<f32>,
    @location(1) normal: vec3<f32>,
    @location(2) chunk: f32,
}

@vertex
fn gbuf_vs(
    @location(0) position: vec3<f32>,
    @location(1) normal: vec3<f32>,
    @location(2) chunk: f32,
) -> VsOut {
    var out: VsOut;
    out.clip = camera.view_proj * vec4<f32>(position, 1.0);
    out.world = position;
    out.normal = normal;
    out.chunk = chunk;
    return out;
}

struct GbufOut {
    @location(0) position: vec4<f32>,
    @location(1) normal: vec4<f32>,
}

// Every triangle belongs to exactly one chunk -- chunks are extracted
// independently -- so all three vertices carry the same `chunk` and its
// interpolant is that value.
@fragment
fn gbuf_fs(vs: VsOut) -> GbufOut {
    var out: GbufOut;
    out.position = vec4<f32>(vs.world, vs.chunk);
    out.normal = vec4<f32>(normalize(vs.normal), 1.0);
    return out;
}

// ── the resident field ──────────────────────────────────────────────────────

fn field_distance(p: vec3<f32>) -> f32 {
    let cells = grid_cells(grid);
    let extent = vec3<f32>(cells);
    let g = (p - grid.placement.xyz) / grid.placement.w;
    if (any(g < vec3<f32>(0.0)) || any(g > extent)) {
        return OUTSIDE;
    }
    let base = min(vec3<u32>(g), cells - vec3<u32>(1u));
    let f = g - vec3<f32>(base);
    let sx = 1u;
    let sy = grid.samples.x;
    let sz = grid.samples.x * grid.samples.y;
    let i = grid_index(grid, base);
    let x00 = mix(distances[i], distances[i + sx], f.x);
    let x10 = mix(distances[i + sy], distances[i + sx + sy], f.x);
    let x01 = mix(distances[i + sz], distances[i + sx + sz], f.x);
    let x11 = mix(distances[i + sy + sz], distances[i + sx + sy + sz], f.x);
    return mix(mix(x00, x10, f.y), mix(x01, x11, f.y), f.z);
}

// Duff et al., "Building an orthonormal basis, revisited", JCGT 2017. Columns
// are (tangent, bitangent, normal), so `frame * d` takes a tangent-space
// direction to world.
fn onb(n: vec3<f32>) -> mat3x3<f32> {
    let s = select(-1.0, 1.0, n.z >= 0.0);
    let a = -1.0 / (s + n.z);
    let b = n.x * n.y * a;
    return mat3x3<f32>(
        vec3<f32>(1.0 + s * n.x * n.x * a, s * b, -s * n.x),
        vec3<f32>(b, s + n.y * n.y * a, -n.y),
        n,
    );
}

// One 16-step sphere trace. The step length is the distance to the nearest
// surface, floored at eps so a grazing ray cannot stall: that floored radius is
// the cone.
fn blocked(origin: vec3<f32>, dir: vec3<f32>) -> f32 {
    let radius = ao.tune.x;
    let eps = ao.tune.y;
    var t = 0.0;
    for (var step = 0u; step < ao.dims.w; step = step + 1u) {
        let d = field_distance(origin + dir * t);
        if (d < eps) {
            return 1.0;
        }
        t = t + max(d, eps);
        if (t > radius) {
            return 0.0;
        }
    }
    // Out of steps, inside the radius, no hit: the truncation the registered
    // 16-step budget buys, counted rather than absorbed. It biases the field arm
    // towards "unoccluded", which is the direction that makes it look good, so
    // it has to be on the row.
    atomicAdd(&diag[3], 1u);
    return 0.0;
}

@compute @workgroup_size(8, 8, 1)
fn ao_field(@builtin(global_invocation_id) gid: vec3<u32>) {
    if (gid.x >= ao.dims.x || gid.y >= ao.dims.y) {
        return;
    }
    let at = vec2<i32>(i32(gid.x), i32(gid.y));
    let index = gid.y * ao.dims.x + gid.x;
    let hit = textureLoad(gbuf_pos, at, 0);
    if (hit.w < 0.0) {
        ao_out[index] = 0.0;
        return;
    }
    atomicAdd(&diag[2], 1u);
    let n = normalize(textureLoad(gbuf_nrm, at, 0).xyz);
    let p = hit.xyz + n * ao.tune.z;
    let start = field_distance(p);
    if (start < 0.0) {
        atomicAdd(&diag[0], 1u);
    }
    if (start < ao.tune.y) {
        atomicAdd(&diag[1], 1u);
    }
    let frame = onb(n);
    var sum = 0.0;
    for (var i = 0u; i < ao.dims.z; i = i + 1u) {
        sum = sum + blocked(p, normalize(frame * dirs[i].xyz));
    }
    ao_out[index] = sum / f32(ao.dims.z);
}

// ── the depth-buffer baseline ───────────────────────────────────────────────

fn reconstruct(at: vec2<i32>, depth: f32) -> vec3<f32> {
    let uv = (vec2<f32>(at) + vec2<f32>(0.5)) / camera.screen.xy;
    let ndc = vec3<f32>(uv.x * 2.0 - 1.0, 1.0 - uv.y * 2.0, depth);
    let h = camera.inv_view_proj * vec4<f32>(ndc, 1.0);
    return h.xyz / h.w;
}

@compute @workgroup_size(8, 8, 1)
fn ao_ssao(@builtin(global_invocation_id) gid: vec3<u32>) {
    if (gid.x >= ao.dims.x || gid.y >= ao.dims.y) {
        return;
    }
    let at = vec2<i32>(i32(gid.x), i32(gid.y));
    let index = gid.y * ao.dims.x + gid.x;
    let centre = textureLoad(gbuf_depth, at, 0);
    if (centre >= 1.0) {
        ao_out[index] = 0.0;
        return;
    }
    let p = reconstruct(at, centre);
    let n = normalize(textureLoad(gbuf_nrm, at, 0).xyz);
    let frame = onb(n);
    let radius = ao.tune.x;
    var sum = 0.0;
    for (var i = 0u; i < ao.dims.z; i = i + 1u) {
        let offset = frame * normalize(dirs[i].xyz) * (radius * dirs[i].w);
        let world = p + offset;
        let clip = camera.view_proj * vec4<f32>(world, 1.0);
        if (clip.w <= 0.0) {
            continue;
        }
        let ndc = clip.xyz / clip.w;
        let uv = vec2<f32>(ndc.x * 0.5 + 0.5, 0.5 - ndc.y * 0.5);
        let px = vec2<i32>(i32(floor(uv.x * camera.screen.x)), i32(floor(uv.y * camera.screen.y)));
        if (px.x < 0 || px.y < 0 || px.x >= i32(ao.dims.x) || px.y >= i32(ao.dims.y)) {
            continue;
        }
        let scene_depth = textureLoad(gbuf_depth, px, 0);
        if (scene_depth >= 1.0) {
            continue;
        }
        let scene = reconstruct(px, scene_depth);
        let scene_distance = length(scene - camera.eye.xyz);
        let sample_distance = length(world - camera.eye.xyz);
        if (scene_distance < sample_distance - ao.tune.w) {
            // The range check: an occluder much closer to the camera than the
            // sample is a different surface, not a fold of this one.
            let delta = sample_distance - scene_distance;
            sum = sum + smoothstep(0.0, 1.0, radius / max(delta, 1.0e-6));
        }
    }
    ao_out[index] = sum / f32(ao.dims.z);
}
"#;

// ─── small maths, so nothing here needs a linear-algebra dependency ──────────

/// A 4x4 matrix, column-major: `m[column][row]`, which is what WGSL reads.
type Mat4 = [[f32; 4]; 4];

fn mat_mul(a: Mat4, b: Mat4) -> Mat4 {
    let mut out = [[0.0f32; 4]; 4];
    for c in 0..4 {
        for r in 0..4 {
            let mut sum = 0.0;
            for k in 0..4 {
                sum += a[k][r] * b[c][k];
            }
            out[c][r] = sum;
        }
    }
    out
}

fn mat_bytes(m: Mat4, out: &mut Vec<u8>) {
    for c in 0..4 {
        for r in 0..4 {
            out.extend_from_slice(&m[c][r].to_le_bytes());
        }
    }
}

fn normalise(v: [f32; 3]) -> [f32; 3] {
    let len = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
    [v[0] / len, v[1] / len, v[2] / len]
}

fn cross(a: [f32; 3], b: [f32; 3]) -> [f32; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn dot(a: [f32; 3], b: [f32; 3]) -> f32 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

/// Right-handed look-at, matching `glam::Mat4::look_at_rh`.
fn look_at(eye: [f32; 3], target: [f32; 3], up: [f32; 3]) -> Mat4 {
    let f = normalise([target[0] - eye[0], target[1] - eye[1], target[2] - eye[2]]);
    let s = normalise(cross(f, up));
    let u = cross(s, f);
    [
        [s[0], u[0], -f[0], 0.0],
        [s[1], u[1], -f[1], 0.0],
        [s[2], u[2], -f[2], 0.0],
        [-dot(s, eye), -dot(u, eye), dot(f, eye), 1.0],
    ]
}

/// Right-handed perspective into wgpu's `z` in `[0, 1]`, matching
/// `glam::Mat4::perspective_rh`.
fn perspective(fov_y: f32, aspect: f32, near: f32, far: f32) -> Mat4 {
    let h = (fov_y * 0.5).cos() / (fov_y * 0.5).sin();
    let w = h / aspect;
    let r = far / (near - far);
    [
        [w, 0.0, 0.0, 0.0],
        [0.0, h, 0.0, 0.0],
        [0.0, 0.0, r, -1.0],
        [0.0, 0.0, r * near, 0.0],
    ]
}

/// The inverse of [`perspective`], derived rather than solved numerically:
/// `x_v = x_c/w`, `y_v = y_c/h`, `z_v = -w_c`, `w_v = z_c/(r*near) + w_c/near`.
fn perspective_inverse(fov_y: f32, aspect: f32, near: f32, far: f32) -> Mat4 {
    let h = (fov_y * 0.5).cos() / (fov_y * 0.5).sin();
    let w = h / aspect;
    let r = far / (near - far);
    [
        [1.0 / w, 0.0, 0.0, 0.0],
        [0.0, 1.0 / h, 0.0, 0.0],
        [0.0, 0.0, 0.0, 1.0 / (r * near)],
        [0.0, 0.0, -1.0, 1.0 / near],
    ]
}

/// The inverse of [`look_at`]: the same basis as columns, plus the eye.
fn look_at_inverse(eye: [f32; 3], target: [f32; 3], up: [f32; 3]) -> Mat4 {
    let f = normalise([target[0] - eye[0], target[1] - eye[1], target[2] - eye[2]]);
    let s = normalise(cross(f, up));
    let u = cross(s, f);
    [
        [s[0], s[1], s[2], 0.0],
        [u[0], u[1], u[2], 0.0],
        [-f[0], -f[1], -f[2], 0.0],
        [eye[0], eye[1], eye[2], 1.0],
    ]
}

/// Duff et al. 2017, byte-for-byte the same expression as the WGSL `onb`.
/// Columns are `(tangent, bitangent, normal)`.
fn onb(n: [f32; 3]) -> [[f32; 3]; 3] {
    let s = if n[2] >= 0.0 { 1.0 } else { -1.0 };
    let a = -1.0 / (s + n[2]);
    let b = n[0] * n[1] * a;
    [
        [1.0 + s * n[0] * n[0] * a, s * b, -s * n[0]],
        [b, s + n[1] * n[1] * a, -n[1]],
        n,
    ]
}

fn apply_frame(frame: &[[f32; 3]; 3], d: [f32; 3]) -> [f32; 3] {
    [
        frame[0][0] * d[0] + frame[1][0] * d[1] + frame[2][0] * d[2],
        frame[0][1] * d[0] + frame[1][1] * d[1] + frame[2][1] * d[2],
        frame[0][2] * d[0] + frame[1][2] * d[1] + frame[2][2] * d[2],
    ]
}

/// Van der Corput radical inverse in base 2.
fn radical_inverse(mut bits: u32) -> f32 {
    bits = bits.rotate_right(16);
    bits = ((bits & 0x5555_5555) << 1) | ((bits & 0xAAAA_AAAA) >> 1);
    bits = ((bits & 0x3333_3333) << 2) | ((bits & 0xCCCC_CCCC) >> 2);
    bits = ((bits & 0x0F0F_0F0F) << 4) | ((bits & 0xF0F0_F0F0) >> 4);
    bits = ((bits & 0x00FF_00FF) << 8) | ((bits & 0xFF00_FF00) >> 8);
    bits as f32 * 2.328_306_4e-10
}

/// `n` cosine-weighted directions in the tangent frame, `z` up.
///
/// Hammersley: `u1 = i/n`, `u2 = radical_inverse(i)`. Cosine weighting is in the
/// mapping, so a plain mean of blocked flags over these directions *is* the
/// cosine-weighted occlusion - which is what makes the 8-ray and 512-ray arms
/// estimates of one number rather than two.
fn hemisphere(n: u32) -> Vec<[f32; 3]> {
    (0..n)
        .map(|i| {
            let u1 = (i as f32 + 0.5) / n as f32;
            let u2 = radical_inverse(i).clamp(0.0, 1.0 - f32::EPSILON);
            let phi = core::f32::consts::TAU * u1;
            let r = u2.sqrt();
            [r * phi.cos(), r * phi.sin(), (1.0 - u2).sqrt()]
        })
        .collect()
}

// ─── the scene ───────────────────────────────────────────────────────────────

/// One field, meshed chunk-wise, in the two forms the harness needs.
struct Scene {
    /// Interleaved `[position; 3], [normal; 3], chunk` per vertex.
    vertices: Vec<f32>,
    indices: Vec<u32>,
    trimesh: TriMesh,
    vertex_count: usize,
    triangle_count: usize,
}

/// Mesh the field over 4^3 chunks of `CHUNK_CELLS`, concatenated into one draw.
///
/// Chunk-wise on purpose: the chunk faces are the seams C2 is about, and a
/// single-chunk extraction would have none. Coincident vertices on shared faces
/// (`A-015`, `M-220`) are left in: they are geometrically the same triangles and
/// neither the rasteriser nor `parry3d` cares.
fn build_scene<F: Sdf<Scalar = f32>>(field: &F) -> Scene {
    let layout = ChunkLayout::<f32>::new(CHUNK_CELLS, CELL, [ORIGIN; 3]).expect("layout");
    let shape = layout.sample_shape().expect("shape");
    let per_axis = WORLD_CELLS / CHUNK_CELLS;
    let mut mc = MarchingCubes::<f32>::new();

    let mut vertices: Vec<f32> = Vec::new();
    let mut indices: Vec<u32> = Vec::new();
    let mut points: Vec<Vector> = Vec::new();
    let mut triangles: Vec<[u32; 3]> = Vec::new();

    let mut chunk = 0u32;
    for cz in 0..per_axis {
        for cy in 0..per_axis {
            for cx in 0..per_axis {
                let id = ChunkId {
                    coords: [cx as i32, cy as i32, cz as i32],
                };
                let mut out = MeshBuffer::<f32>::new();
                let _ = mc.extract_into(field, &shape, layout.sample_origin(id), CELL, &mut out);
                let base = (vertices.len() / 7) as u32;
                for (p, n) in out.positions.iter().zip(&out.normals) {
                    vertices.extend_from_slice(&[p[0], p[1], p[2], n[0], n[1], n[2], chunk as f32]);
                    points.push(Vector::new(p[0], p[1], p[2]));
                }
                for t in out.indices.as_chunks::<3>().0 {
                    let tri = [base + t[0], base + t[1], base + t[2]];
                    indices.extend_from_slice(&tri);
                    triangles.push(tri);
                }
                chunk += 1;
            }
        }
    }

    let vertex_count = vertices.len() / 7;
    let triangle_count = triangles.len();
    let trimesh = TriMesh::new(points, triangles).expect("a parry TriMesh over the extracted mesh");
    Scene {
        vertices,
        indices,
        trimesh,
        vertex_count,
        triangle_count,
    }
}

/// AO by casting `dirs` against the mesh. The reference, and the 8-ray floor.
fn ao_by_mesh(
    trimesh: &TriMesh,
    position: [f32; 3],
    normal: [f32; 3],
    dirs: &[[f32; 3]],
    radius: f32,
    lift: f32,
) -> f32 {
    let frame = onb(normal);
    let origin = Vector::new(
        position[0] + normal[0] * lift,
        position[1] + normal[1] * lift,
        position[2] + normal[2] * lift,
    );
    let mut blocked = 0u32;
    for d in dirs {
        let w = normalise(apply_frame(&frame, *d));
        let ray = Ray::new(origin, Vector::new(w[0], w[1], w[2]));
        if trimesh.cast_local_ray(&ray, radius, true).is_some() {
            blocked += 1;
        }
    }
    blocked as f32 / dirs.len() as f32
}

/// Largest `|f(a) - f(b)| / h` over axis-adjacent samples.
///
/// The control on C3's premise. A sphere trace is only sound on a field whose
/// gradient magnitude is at most 1; a value above 1 means the resident field
/// *cannot* be traced as it stands and the flood is a mandatory conversion
/// rather than an optimisation.
fn lipschitz(samples: &[f32], per_axis: u32, h: f32) -> f32 {
    let s = per_axis as usize;
    let mut worst = 0.0f32;
    let strides = [1usize, s, s * s];
    for axis in 0..3 {
        let stride = strides[axis];
        for z in 0..s {
            for y in 0..s {
                for x in 0..s {
                    let at = [x, y, z];
                    if at[axis] + 1 >= s {
                        continue;
                    }
                    let i = x + y * s + z * s * s;
                    let d = (samples[i + stride] - samples[i]).abs() / h;
                    if d.is_finite() && d > worst {
                        worst = d;
                    }
                }
            }
        }
    }
    worst
}

type Row = Vec<(&'static str, String)>;

fn main() {
    if !std::env::args().any(|arg| arg == "--bench") {
        return;
    }

    let prereg = isomesh::experiment!("P-74");

    // ── C1's arithmetic, printed before anything runs ───────────────────────
    let pixels = u64::from(WIDTH) * u64::from(HEIGHT);
    let queries = pixels * u64::from(RAYS) * u64::from(STEPS);
    let needed = queries as f64 / 2.0e-3;
    println!("C1's SHARE is the whole AO budget, so the question is reachability:");
    println!("  {pixels} pixels x {RAYS} rays x {STEPS} steps = {queries} distance queries/frame");
    println!(
        "  2.0 ms demands {:.1} G queries/s; nvblox measured 0.8-7.3 G/s on an RTX 3090 Ti",
        needed / 1e9
    );
    println!(
        "  so nvblox's best rate implies {:.1} ms, and RTSDF's own 4.60 ms trace at 1024^2 \
         scales to 9.10 ms on a 2080 Ti",
        queries as f64 / 7.3e9 * 1e3
    );
    println!(
        "  the flat {SAMPLES}^3 grid is {:.3} MB against 6 MB of L2: 8 loads/query is 8.49 GB \
         of traffic, 0.48 ms out of L1 and 9.1 ms out of VRAM",
        f64::from(SAMPLES).powi(3) * 4.0 / 1e6
    );
    println!("  => reachable only with early termination and locality. Running.\n");

    // ── the device ──────────────────────────────────────────────────────────
    let gpu = headless::Gpu::with_timestamps()
        .expect("a device with TIMESTAMP_QUERY; C1 has no honest clock without it");
    let device = gpu.device();
    let queue = gpu.queue();
    let adapter = gpu.report().name.replace([',', '"'], " ");
    let backend = format!("{:?}", gpu.report().backend);
    println!("adapter {adapter} on {backend}\n");

    let stamps = isomesh_gpu::StageTimestamps::new(device, queue).expect("timestamps");

    // ── the module ──────────────────────────────────────────────────────────
    let mut composer = Composer::with_builtins();
    composer.insert("p74", P74_WGSL);
    let source = composer.compose("p74", &[]).expect("compose p74");
    let module = device.create_shader_module(wgpu::ShaderModuleDescriptor {
        label: Some("p74 ao"),
        source: wgpu::ShaderSource::Wgsl(source.into()),
    });

    // ── layouts ─────────────────────────────────────────────────────────────
    let uniform = |binding: u32, visibility: wgpu::ShaderStages| wgpu::BindGroupLayoutEntry {
        binding,
        visibility,
        ty: wgpu::BindingType::Buffer {
            ty: wgpu::BufferBindingType::Uniform,
            has_dynamic_offset: false,
            min_binding_size: None,
        },
        count: None,
    };
    let storage = |binding: u32, read_only: bool| wgpu::BindGroupLayoutEntry {
        binding,
        visibility: wgpu::ShaderStages::COMPUTE,
        ty: wgpu::BindingType::Buffer {
            ty: wgpu::BufferBindingType::Storage { read_only },
            has_dynamic_offset: false,
            min_binding_size: None,
        },
        count: None,
    };
    let colour_texture = |binding: u32| wgpu::BindGroupLayoutEntry {
        binding,
        visibility: wgpu::ShaderStages::COMPUTE,
        ty: wgpu::BindingType::Texture {
            sample_type: wgpu::TextureSampleType::Float { filterable: false },
            view_dimension: wgpu::TextureViewDimension::D2,
            multisampled: false,
        },
        count: None,
    };

    let render_layout = device.create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
        label: Some("p74 camera"),
        entries: &[uniform(0, wgpu::ShaderStages::VERTEX)],
    });
    let compute_layout = device.create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
        label: Some("p74 ao"),
        entries: &[
            uniform(0, wgpu::ShaderStages::COMPUTE),
            uniform(1, wgpu::ShaderStages::COMPUTE),
            uniform(2, wgpu::ShaderStages::COMPUTE),
            storage(3, true),
            storage(4, true),
            colour_texture(5),
            colour_texture(6),
            wgpu::BindGroupLayoutEntry {
                binding: 7,
                visibility: wgpu::ShaderStages::COMPUTE,
                ty: wgpu::BindingType::Texture {
                    sample_type: wgpu::TextureSampleType::Depth,
                    view_dimension: wgpu::TextureViewDimension::D2,
                    multisampled: false,
                },
                count: None,
            },
            storage(8, false),
            storage(9, false),
        ],
    });

    let render_pipeline_layout = device.create_pipeline_layout(&wgpu::PipelineLayoutDescriptor {
        label: Some("p74 gbuffer"),
        bind_group_layouts: &[Some(&render_layout)],
        immediate_size: 0,
    });
    let compute_pipeline_layout = device.create_pipeline_layout(&wgpu::PipelineLayoutDescriptor {
        label: Some("p74 ao"),
        bind_group_layouts: &[Some(&compute_layout)],
        immediate_size: 0,
    });

    let gbuffer_pipeline = device.create_render_pipeline(&wgpu::RenderPipelineDescriptor {
        label: Some("p74 gbuffer"),
        layout: Some(&render_pipeline_layout),
        vertex: wgpu::VertexState {
            module: &module,
            entry_point: Some("gbuf_vs"),
            compilation_options: wgpu::PipelineCompilationOptions::default(),
            buffers: &[wgpu::VertexBufferLayout {
                array_stride: 28,
                step_mode: wgpu::VertexStepMode::Vertex,
                attributes: &[
                    wgpu::VertexAttribute {
                        format: wgpu::VertexFormat::Float32x3,
                        offset: 0,
                        shader_location: 0,
                    },
                    wgpu::VertexAttribute {
                        format: wgpu::VertexFormat::Float32x3,
                        offset: 12,
                        shader_location: 1,
                    },
                    wgpu::VertexAttribute {
                        format: wgpu::VertexFormat::Float32,
                        offset: 24,
                        shader_location: 2,
                    },
                ],
            }],
        },
        primitive: wgpu::PrimitiveState {
            topology: wgpu::PrimitiveTopology::TriangleList,
            strip_index_format: None,
            front_face: wgpu::FrontFace::Ccw,
            // Marching Cubes emits a consistent winding, but AO does not depend
            // on which side a triangle faces and culling here would only be a
            // way to lose geometry silently.
            cull_mode: None,
            unclipped_depth: false,
            polygon_mode: wgpu::PolygonMode::Fill,
            conservative: false,
        },
        depth_stencil: Some(wgpu::DepthStencilState {
            format: wgpu::TextureFormat::Depth32Float,
            depth_write_enabled: Some(true),
            depth_compare: Some(wgpu::CompareFunction::Less),
            stencil: wgpu::StencilState::default(),
            bias: wgpu::DepthBiasState::default(),
        }),
        multisample: wgpu::MultisampleState::default(),
        fragment: Some(wgpu::FragmentState {
            module: &module,
            entry_point: Some("gbuf_fs"),
            compilation_options: wgpu::PipelineCompilationOptions::default(),
            targets: &[
                Some(wgpu::ColorTargetState {
                    format: wgpu::TextureFormat::Rgba32Float,
                    blend: None,
                    write_mask: wgpu::ColorWrites::ALL,
                }),
                Some(wgpu::ColorTargetState {
                    format: wgpu::TextureFormat::Rgba32Float,
                    blend: None,
                    write_mask: wgpu::ColorWrites::ALL,
                }),
            ],
        }),
        multiview_mask: None,
        cache: None,
    });

    let make_compute = |entry: &str| {
        device.create_compute_pipeline(&wgpu::ComputePipelineDescriptor {
            label: Some(entry),
            layout: Some(&compute_pipeline_layout),
            module: &module,
            entry_point: Some(entry),
            compilation_options: wgpu::PipelineCompilationOptions::default(),
            cache: None,
        })
    };
    let field_pipeline = make_compute("ao_field");
    let ssao_pipeline = make_compute("ao_ssao");

    // ── targets, allocated once for both fields ─────────────────────────────
    let extent = wgpu::Extent3d {
        width: WIDTH,
        height: HEIGHT,
        depth_or_array_layers: 1,
    };
    let target = |label: &str, format: wgpu::TextureFormat, usage: wgpu::TextureUsages| {
        device.create_texture(&wgpu::TextureDescriptor {
            label: Some(label),
            size: extent,
            mip_level_count: 1,
            sample_count: 1,
            dimension: wgpu::TextureDimension::D2,
            format,
            usage,
            view_formats: &[],
        })
    };
    let colour_usage = wgpu::TextureUsages::RENDER_ATTACHMENT
        | wgpu::TextureUsages::TEXTURE_BINDING
        | wgpu::TextureUsages::COPY_SRC;
    let pos_texture = target(
        "p74 position",
        wgpu::TextureFormat::Rgba32Float,
        colour_usage,
    );
    let nrm_texture = target("p74 normal", wgpu::TextureFormat::Rgba32Float, colour_usage);
    let depth_texture = target(
        "p74 depth",
        wgpu::TextureFormat::Depth32Float,
        wgpu::TextureUsages::RENDER_ATTACHMENT | wgpu::TextureUsages::TEXTURE_BINDING,
    );
    let pos_view = pos_texture.create_view(&wgpu::TextureViewDescriptor::default());
    let nrm_view = nrm_texture.create_view(&wgpu::TextureViewDescriptor::default());
    let depth_view = depth_texture.create_view(&wgpu::TextureViewDescriptor::default());

    let gbuffer_bytes = pixels * 16;
    let readback = |label: &str, size: u64| {
        device.create_buffer(&wgpu::BufferDescriptor {
            label: Some(label),
            size,
            usage: wgpu::BufferUsages::COPY_DST | wgpu::BufferUsages::COPY_SRC,
            mapped_at_creation: false,
        })
    };
    let pos_staging = readback("p74 position staging", gbuffer_bytes);
    let nrm_staging = readback("p74 normal staging", gbuffer_bytes);

    let ao_bytes = pixels * 4;
    let ao_buffer = device.create_buffer(&wgpu::BufferDescriptor {
        label: Some("p74 ao"),
        size: ao_bytes,
        usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_SRC,
        mapped_at_creation: false,
    });
    let diag_buffer = device.create_buffer(&wgpu::BufferDescriptor {
        label: Some("p74 diagnostics"),
        size: 16,
        usage: wgpu::BufferUsages::STORAGE
            | wgpu::BufferUsages::COPY_DST
            | wgpu::BufferUsages::COPY_SRC,
        mapped_at_creation: false,
    });

    // ── the shared direction table ──────────────────────────────────────────
    let dirs8 = hemisphere(RAYS);
    let dirs_ref = hemisphere(REF_RAYS);
    let mut dir_bytes = Vec::with_capacity(dirs8.len() * 16);
    for (i, d) in dirs8.iter().enumerate() {
        // The SSAO radius scale: the cube root distributes sample points
        // uniformly through the hemisphere's volume rather than bunching them at
        // the surface, which is what the classic kernel's ad-hoc `lerp(0.1, 1)`
        // approximates.
        let scale = ((i as f32 + 0.5) / RAYS as f32).cbrt();
        for value in [d[0], d[1], d[2], scale] {
            dir_bytes.extend_from_slice(&value.to_le_bytes());
        }
    }
    let dirs_buffer = device.create_buffer(&wgpu::BufferDescriptor {
        label: Some("p74 directions"),
        size: dir_bytes.len() as u64,
        usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_DST,
        mapped_at_creation: false,
    });
    queue.write_buffer(&dirs_buffer, 0, &dir_bytes);

    // ── camera and AO uniforms ──────────────────────────────────────────────
    let centre = [ORIGIN + EXTENT * 0.5; 3];
    let look = normalise([0.6, 0.45, 1.0]);
    let eye = [
        centre[0] + look[0] * EYE_DISTANCE,
        centre[1] + look[1] * EYE_DISTANCE,
        centre[2] + look[2] * EYE_DISTANCE,
    ];
    let up = [0.0, 1.0, 0.0];
    let aspect = WIDTH as f32 / HEIGHT as f32;
    let view_proj = mat_mul(
        perspective(FOV_Y, aspect, NEAR, FAR),
        look_at(eye, centre, up),
    );
    let inv_view_proj = mat_mul(
        look_at_inverse(eye, centre, up),
        perspective_inverse(FOV_Y, aspect, NEAR, FAR),
    );
    let mut camera_bytes = Vec::with_capacity(160);
    mat_bytes(view_proj, &mut camera_bytes);
    mat_bytes(inv_view_proj, &mut camera_bytes);
    for value in [eye[0], eye[1], eye[2], 0.0] {
        camera_bytes.extend_from_slice(&value.to_le_bytes());
    }
    for value in [WIDTH as f32, HEIGHT as f32, NEAR, FAR] {
        camera_bytes.extend_from_slice(&value.to_le_bytes());
    }
    let camera_buffer = device.create_buffer(&wgpu::BufferDescriptor {
        label: Some("p74 camera"),
        size: camera_bytes.len() as u64,
        usage: wgpu::BufferUsages::UNIFORM | wgpu::BufferUsages::COPY_DST,
        mapped_at_creation: false,
    });
    queue.write_buffer(&camera_buffer, 0, &camera_bytes);

    let radius = AO_RADIUS_CELLS * CELL;
    let eps = EPS_CELLS * CELL;
    let lift = LIFT_CELLS * CELL;
    let mut ao_param_bytes = Vec::with_capacity(32);
    for value in [WIDTH, HEIGHT, RAYS, STEPS] {
        ao_param_bytes.extend_from_slice(&value.to_le_bytes());
    }
    for value in [radius, eps, lift, SSAO_BIAS_CELLS * CELL] {
        ao_param_bytes.extend_from_slice(&value.to_le_bytes());
    }
    let ao_param_buffer = device.create_buffer(&wgpu::BufferDescriptor {
        label: Some("p74 ao params"),
        size: ao_param_bytes.len() as u64,
        usage: wgpu::BufferUsages::UNIFORM | wgpu::BufferUsages::COPY_DST,
        mapped_at_creation: false,
    });
    queue.write_buffer(&ao_param_buffer, 0, &ao_param_bytes);

    let render_bind = device.create_bind_group(&wgpu::BindGroupDescriptor {
        label: Some("p74 camera"),
        layout: &render_layout,
        entries: &[wgpu::BindGroupEntry {
            binding: 0,
            resource: camera_buffer.as_entire_binding(),
        }],
    });

    // ── the resident grid ───────────────────────────────────────────────────
    let grid = GridParams::new([SAMPLES; 3], [ORIGIN; 3], CELL).expect("grid");
    let grid_buffer = device.create_buffer(&wgpu::BufferDescriptor {
        label: Some("p74 grid"),
        size: GridParams::UNIFORM_SIZE,
        usage: wgpu::BufferUsages::UNIFORM | wgpu::BufferUsages::COPY_DST,
        mapped_at_creation: false,
    });
    queue.write_buffer(&grid_buffer, 0, &grid.to_std140());
    let distance_bytes = grid.sample_count() * 4;
    let distance_buffer = device.create_buffer(&wgpu::BufferDescriptor {
        label: Some("p74 resident distances"),
        size: distance_bytes,
        usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_DST,
        mapped_at_creation: false,
    });

    let compute_bind = device.create_bind_group(&wgpu::BindGroupDescriptor {
        label: Some("p74 ao"),
        layout: &compute_layout,
        entries: &[
            wgpu::BindGroupEntry {
                binding: 0,
                resource: camera_buffer.as_entire_binding(),
            },
            wgpu::BindGroupEntry {
                binding: 1,
                resource: ao_param_buffer.as_entire_binding(),
            },
            wgpu::BindGroupEntry {
                binding: 2,
                resource: grid_buffer.as_entire_binding(),
            },
            wgpu::BindGroupEntry {
                binding: 3,
                resource: distance_buffer.as_entire_binding(),
            },
            wgpu::BindGroupEntry {
                binding: 4,
                resource: dirs_buffer.as_entire_binding(),
            },
            wgpu::BindGroupEntry {
                binding: 5,
                resource: wgpu::BindingResource::TextureView(&pos_view),
            },
            wgpu::BindGroupEntry {
                binding: 6,
                resource: wgpu::BindingResource::TextureView(&nrm_view),
            },
            wgpu::BindGroupEntry {
                binding: 7,
                resource: wgpu::BindingResource::TextureView(&depth_view),
            },
            wgpu::BindGroupEntry {
                binding: 8,
                resource: ao_buffer.as_entire_binding(),
            },
            wgpu::BindGroupEntry {
                binding: 9,
                resource: diag_buffer.as_entire_binding(),
            },
        ],
    });

    let groups = (WIDTH.div_ceil(8), HEIGHT.div_ceil(8));
    let wait = || {
        let _ = device.poll(wgpu::PollType::Wait {
            submission_index: None,
            timeout: None,
        });
    };

    // ── C3's instrument, and the arm that proves it can say "yes" ───────────
    //
    // `lipschitz` decides whether the resident field is sphere-traceable. Run on
    // the two registered fields it will be asked for a number above 1; run on an
    // exact SDF over the same grid, through the same function, it must come out
    // at 1. Without this arm, "1.41 > 1, therefore a conversion is needed" rests
    // on an instrument that has only ever been shown failing - `M-44` in the
    // other direction.
    let sphere_lipschitz = {
        let field = isomesh::fields::Sphere::<f32>::canonical();
        let mut samples = Vec::with_capacity(grid.sample_count() as usize);
        for z in 0..SAMPLES {
            for y in 0..SAMPLES {
                for x in 0..SAMPLES {
                    samples.push(field.sample(grid.sample_position([x, y, z])));
                }
            }
        }
        lipschitz(&samples, SAMPLES, CELL)
    };
    println!(
        "C3's instrument, checked against an exact SDF on the same grid: sphere Lipschitz \
         {sphere_lipschitz:.6}\n"
    );
    assert!(
        sphere_lipschitz <= 1.01,
        "the Lipschitz instrument reports {sphere_lipschitz:.6} for `sphere`, which IS a \
         distance field, so it is measuring the grid rather than the field and cannot be used \
         to say the registered fields need converting"
    );

    let mut rows: Vec<Row> = Vec::new();

    for field_name in ["gyroid", "fbm_terrain"] {
        println!("── {field_name} ──────────────────────────────────────────────");

        // ── mesh, and the field samples both arms share ──────────────────────
        let (scene, samples) = match field_name {
            "gyroid" => {
                let field = Gyroid::<f32>::canonical();
                let scene = build_scene(&field);
                let mut samples = Vec::with_capacity(grid.sample_count() as usize);
                for z in 0..SAMPLES {
                    for y in 0..SAMPLES {
                        for x in 0..SAMPLES {
                            samples.push(field.sample(grid.sample_position([x, y, z])));
                        }
                    }
                }
                (scene, samples)
            }
            _ => {
                let field = FbmTerrain::<f32>::canonical();
                let scene = build_scene(&field);
                let mut samples = Vec::with_capacity(grid.sample_count() as usize);
                for z in 0..SAMPLES {
                    for y in 0..SAMPLES {
                        for x in 0..SAMPLES {
                            samples.push(field.sample(grid.sample_position([x, y, z])));
                        }
                    }
                }
                (scene, samples)
            }
        };
        assert!(
            scene.triangle_count > 0,
            "VOID: {field_name} produced no geometry over the world box, so nothing would be \
             rasterised and every count below would be zero for a reason that is not a finding"
        );
        let raw_lipschitz = lipschitz(&samples, SAMPLES, CELL);
        println!(
            "  mesh {} vertices, {} triangles; raw field Lipschitz {raw_lipschitz:.4}",
            scene.vertex_count, scene.triangle_count
        );

        let mut vertex_bytes = Vec::with_capacity(scene.vertices.len() * 4);
        for value in &scene.vertices {
            vertex_bytes.extend_from_slice(&value.to_le_bytes());
        }
        let mut index_bytes = Vec::with_capacity(scene.indices.len() * 4);
        for value in &scene.indices {
            index_bytes.extend_from_slice(&value.to_le_bytes());
        }
        let vertex_buffer = device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("p74 vertices"),
            size: vertex_bytes.len() as u64,
            usage: wgpu::BufferUsages::VERTEX | wgpu::BufferUsages::COPY_DST,
            mapped_at_creation: false,
        });
        queue.write_buffer(&vertex_buffer, 0, &vertex_bytes);
        let index_buffer = device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("p74 indices"),
            size: index_bytes.len() as u64,
            usage: wgpu::BufferUsages::INDEX | wgpu::BufferUsages::COPY_DST,
            mapped_at_creation: false,
        });
        queue.write_buffer(&index_buffer, 0, &index_bytes);

        // ── the G-buffer, and the depth image SSAO will read ────────────────
        let render_once = || {
            let mut encoder =
                device.create_command_encoder(&wgpu::CommandEncoderDescriptor { label: None });
            {
                let mut pass = encoder.begin_render_pass(&wgpu::RenderPassDescriptor {
                    label: Some("p74 gbuffer"),
                    color_attachments: &[
                        Some(wgpu::RenderPassColorAttachment {
                            view: &pos_view,
                            depth_slice: None,
                            resolve_target: None,
                            ops: wgpu::Operations {
                                // `a = -1` is the miss marker: a pixel the
                                // rasteriser never wrote carries a negative
                                // chunk id and every arm skips it.
                                load: wgpu::LoadOp::Clear(wgpu::Color {
                                    r: 0.0,
                                    g: 0.0,
                                    b: 0.0,
                                    a: -1.0,
                                }),
                                store: wgpu::StoreOp::Store,
                            },
                        }),
                        Some(wgpu::RenderPassColorAttachment {
                            view: &nrm_view,
                            depth_slice: None,
                            resolve_target: None,
                            ops: wgpu::Operations {
                                load: wgpu::LoadOp::Clear(wgpu::Color::TRANSPARENT),
                                store: wgpu::StoreOp::Store,
                            },
                        }),
                    ],
                    depth_stencil_attachment: Some(wgpu::RenderPassDepthStencilAttachment {
                        view: &depth_view,
                        depth_ops: Some(wgpu::Operations {
                            load: wgpu::LoadOp::Clear(1.0),
                            store: wgpu::StoreOp::Store,
                        }),
                        stencil_ops: None,
                    }),
                    timestamp_writes: None,
                    occlusion_query_set: None,
                    multiview_mask: None,
                });
                pass.set_pipeline(&gbuffer_pipeline);
                pass.set_bind_group(0, &render_bind, &[]);
                pass.set_vertex_buffer(0, vertex_buffer.slice(..));
                pass.set_index_buffer(index_buffer.slice(..), wgpu::IndexFormat::Uint32);
                pass.draw_indexed(0..scene.indices.len() as u32, 0, 0..1);
            }
            queue.submit(Some(encoder.finish()));
            wait();
        };
        render_once();
        let started = Instant::now();
        render_once();
        let gbuffer_ms = started.elapsed().as_nanos() as f64 / 1e6;

        // ── the G-buffer, on the CPU, for classification and the reference ──
        {
            let mut encoder =
                device.create_command_encoder(&wgpu::CommandEncoderDescriptor { label: None });
            for (texture, staging) in [(&pos_texture, &pos_staging), (&nrm_texture, &nrm_staging)] {
                encoder.copy_texture_to_buffer(
                    wgpu::TexelCopyTextureInfo {
                        texture,
                        mip_level: 0,
                        origin: wgpu::Origin3d::ZERO,
                        aspect: wgpu::TextureAspect::All,
                    },
                    wgpu::TexelCopyBufferInfo {
                        buffer: staging,
                        layout: wgpu::TexelCopyBufferLayout {
                            offset: 0,
                            bytes_per_row: Some(WIDTH * 16),
                            rows_per_image: Some(HEIGHT),
                        },
                    },
                    extent,
                );
            }
            queue.submit(Some(encoder.finish()));
            wait();
        }
        let read_f32 = |buffer: &wgpu::Buffer, bytes: u64| -> Vec<f32> {
            let raw = read_bytes(device, queue, buffer, bytes).expect("readback");
            raw.as_chunks::<4>()
                .0
                .iter()
                .map(|w| f32::from_le_bytes(*w))
                .collect()
        };
        let gbuf_pos = read_f32(&pos_staging, gbuffer_bytes);
        let gbuf_nrm = read_f32(&nrm_staging, gbuffer_bytes);

        // ── the two pixel sets, and the negative arms of each classifier ────
        let count = pixels as usize;
        let hit: Vec<bool> = (0..count).map(|i| gbuf_pos[i * 4 + 3] >= 0.0).collect();
        let view_distance: Vec<f32> = (0..count)
            .map(|i| {
                let p = [gbuf_pos[i * 4], gbuf_pos[i * 4 + 1], gbuf_pos[i * 4 + 2]];
                ((p[0] - eye[0]).powi(2) + (p[1] - eye[1]).powi(2) + (p[2] - eye[2]).powi(2)).sqrt()
            })
            .collect();
        let hit_pixels = hit.iter().filter(|h| **h).count();
        assert!(
            hit_pixels > 0,
            "VOID: the camera sees no geometry on {field_name}, so both registered pixel sets \
             would be empty for want of a scene rather than for want of a seam"
        );

        // The interior chunk-boundary planes: 3 per axis at CHUNK_CELLS steps.
        let per_axis = WORLD_CELLS / CHUNK_CELLS;
        let planes: Vec<f32> = (1..per_axis)
            .map(|k| ORIGIN + (k * CHUNK_CELLS) as f32 * CELL)
            .collect();
        let seam_band = SEAM_BAND_CELLS * CELL;
        let classify_seam = |planes: &[f32]| -> Vec<usize> {
            (0..count)
                .filter(|&i| {
                    hit[i]
                        && (0..3).any(|a| {
                            let c = gbuf_pos[i * 4 + a];
                            planes.iter().any(|p| (c - p).abs() < seam_band)
                        })
                })
                .collect()
        };
        let seam_set = classify_seam(&planes);
        // The negative arm: the same classifier, the same image, one chunk.
        let seam_one_chunk = classify_seam(&[]).len();

        let step = SILHOUETTE_STEP_CELLS * CELL;
        let classify_silhouette = |hit: &[bool], distance: &[f32]| -> Vec<usize> {
            let w = WIDTH as usize;
            let h = HEIGHT as usize;
            let mut out = Vec::new();
            for y in 0..h {
                for x in 0..w {
                    let i = y * w + x;
                    if !hit[i] {
                        continue;
                    }
                    let mut edge = false;
                    for (dx, dy) in [(-1i32, 0i32), (1, 0), (0, -1), (0, 1)] {
                        let nx = x as i32 + dx;
                        let ny = y as i32 + dy;
                        if nx < 0 || ny < 0 || nx >= w as i32 || ny >= h as i32 {
                            continue;
                        }
                        let j = ny as usize * w + nx as usize;
                        if !hit[j] || (distance[i] - distance[j]).abs() > step {
                            edge = true;
                        }
                    }
                    if edge {
                        out.push(i);
                    }
                }
            }
            out
        };
        let silhouette_set = classify_silhouette(&hit, &view_distance);
        // The negative arm: the same classifier over an image with no
        // discontinuity anywhere -- every pixel a hit, at one distance.
        let flat_hit = vec![true; count];
        let flat_distance = vec![1.0f32; count];
        let silhouette_flat = classify_silhouette(&flat_hit, &flat_distance).len();

        // The independent second definition of the seam set.
        let seam_screen = {
            let w = WIDTH as usize;
            let h = HEIGHT as usize;
            let mut n = 0usize;
            for y in 0..h {
                for x in 0..w {
                    let i = y * w + x;
                    if !hit[i] {
                        continue;
                    }
                    let mine = gbuf_pos[i * 4 + 3].round() as i32;
                    let mut differs = false;
                    for (dx, dy) in [(-1i32, 0i32), (1, 0), (0, -1), (0, 1)] {
                        let nx = x as i32 + dx;
                        let ny = y as i32 + dy;
                        if nx < 0 || ny < 0 || nx >= w as i32 || ny >= h as i32 {
                            continue;
                        }
                        let j = ny as usize * w + nx as usize;
                        if hit[j] && gbuf_pos[j * 4 + 3].round() as i32 != mine {
                            differs = true;
                        }
                    }
                    if differs {
                        n += 1;
                    }
                }
            }
            n
        };

        println!(
            "  hits {hit_pixels}, seam {} (screen {seam_screen}, one-chunk control \
             {seam_one_chunk}), silhouette {} (flat control {silhouette_flat})",
            seam_set.len(),
            silhouette_set.len()
        );

        // ── the vacuity control, asserted ──────────────────────────────────
        assert!(
            !seam_set.is_empty(),
            "VOID: no seam pixels on {field_name}; a scene with no visible seam cannot \
             distinguish the two methods, which is exactly what the registration's vacuity \
             control names"
        );
        assert!(
            !silhouette_set.is_empty(),
            "VOID: no silhouette pixels on {field_name}; the haloing clause would have no \
             population"
        );
        assert_eq!(
            seam_one_chunk, 0,
            "the seam classifier reports {seam_one_chunk} seam pixels on a world with no \
             interior chunk boundary, so it is counting pixels rather than reading them and \
             the non-zero count above proves nothing"
        );
        assert_eq!(
            silhouette_flat, 0,
            "the silhouette classifier reports {silhouette_flat} pixels on an image with no \
             depth discontinuity anywhere, so its non-zero count above proves nothing"
        );
        assert!(
            seam_set.len() < hit_pixels,
            "every visible pixel is a seam pixel on {field_name}, so the seam set is the whole \
             image and the clause is not about seams"
        );

        // ── the field buffer, the flood, and the conversion between them ────
        let field_upload_started = Instant::now();
        let field_buffer =
            FieldBuffer::uploaded(device, queue, grid, &samples).expect("field buffer");
        queue.submit(std::iter::empty());
        wait();
        let field_upload_ms = field_upload_started.elapsed().as_nanos() as f64 / 1e6;

        let flood = JumpFlood::new(device).expect("jump flood");
        // Warm, so the first build's shader-cache and allocator costs are not
        // charged to C3.
        let _ = flood
            .build(device, queue, &field_buffer)
            .expect("flood warm");
        let flood_started = Instant::now();
        let distances = flood.build(device, queue, &field_buffer).expect("flood");
        let flood_build_ms = flood_started.elapsed().as_nanos() as f64 / 1e6;
        let flood_lipschitz = lipschitz(&distances, SAMPLES, CELL);

        let mut distance_bytes_host = Vec::with_capacity(distances.len() * 4);
        for value in &distances {
            distance_bytes_host.extend_from_slice(&value.to_le_bytes());
        }
        let upload_started = Instant::now();
        queue.write_buffer(&distance_buffer, 0, &distance_bytes_host);
        queue.submit(std::iter::empty());
        wait();
        let distance_upload_ms = upload_started.elapsed().as_nanos() as f64 / 1e6;
        println!(
            "  field upload {field_upload_ms:.3} ms, flood build {flood_build_ms:.3} ms, \
             distance re-upload {distance_upload_ms:.3} ms ({:.3} MB each way); flooded \
             Lipschitz {flood_lipschitz:.4}",
            distance_bytes as f64 / 1e6
        );
        assert!(
            flood_lipschitz.is_finite() && flood_lipschitz > 0.0,
            "the flooded field is constant, so the trace would have nothing to hit"
        );

        // ── the timed arms ─────────────────────────────────────────────────
        let run_ao =
            |pipeline: &wgpu::ComputePipeline, label: &'static str, stamped: bool| -> f64 {
                let mut encoder =
                    device.create_command_encoder(&wgpu::CommandEncoderDescriptor { label: None });
                {
                    let mut pass = encoder.begin_compute_pass(&wgpu::ComputePassDescriptor {
                        label: Some(label),
                        timestamp_writes: if stamped { stamps.writes(label) } else { None },
                    });
                    pass.set_pipeline(pipeline);
                    pass.set_bind_group(0, &compute_bind, &[]);
                    pass.dispatch_workgroups(groups.0, groups.1, 1);
                }
                let started = Instant::now();
                queue.submit(Some(encoder.finish()));
                wait();
                started.elapsed().as_nanos() as f64 / 1e6
            };

        let mut field_gpu: Vec<f64> = Vec::with_capacity(REPS);
        let mut field_wall: Vec<f64> = Vec::with_capacity(REPS);
        let _ = run_ao(&field_pipeline, "ao_field warm", false);
        for _ in 0..REPS {
            let wall = run_ao(&field_pipeline, "ao_field", true);
            let spans = stamps.resolve(device, queue).expect("spans");
            assert!(spans.complete, "the timestamp set overflowed");
            field_gpu.push(spans.spans[0].ms);
            field_wall.push(wall);
        }

        let mut ssao_gpu: Vec<f64> = Vec::with_capacity(REPS);
        let _ = run_ao(&ssao_pipeline, "ao_ssao warm", false);
        for _ in 0..REPS {
            let _ = run_ao(&ssao_pipeline, "ao_ssao", true);
            let spans = stamps.resolve(device, queue).expect("spans");
            assert!(spans.complete, "the timestamp set overflowed");
            ssao_gpu.push(spans.spans[0].ms);
        }

        // C3: the same trace, but the field is rebuilt first. `JumpFlood::build`
        // returns a `Vec<f32>` -- the distances land on the CPU -- so the honest
        // rebuild is flood + download + re-upload + trace, and the download is
        // inside `build`.
        let mut rebuild_wall: Vec<f64> = Vec::with_capacity(REPS);
        for _ in 0..REPS {
            let started = Instant::now();
            let rebuilt = flood.build(device, queue, &field_buffer).expect("rebuild");
            let mut bytes = Vec::with_capacity(rebuilt.len() * 4);
            for value in &rebuilt {
                bytes.extend_from_slice(&value.to_le_bytes());
            }
            queue.write_buffer(&distance_buffer, 0, &bytes);
            let _ = run_ao(&field_pipeline, "ao_field rebuild", false);
            rebuild_wall.push(started.elapsed().as_nanos() as f64 / 1e6);
        }

        let median = |mut v: Vec<f64>| -> f64 {
            v.sort_unstable_by(|a, b| a.partial_cmp(b).expect("finite"));
            v[v.len() / 2]
        };
        let ao_ms_field = median(field_gpu);
        let ao_ms_ssao = median(ssao_gpu);
        let resident_ms = median(field_wall);
        let rebuild_ms = median(rebuild_wall);

        // ── the images the MAEs are computed from, in one canonical run ─────
        queue.write_buffer(&diag_buffer, 0, &[0u8; 16]);
        queue.submit(std::iter::empty());
        let _ = run_ao(&field_pipeline, "ao_field record", false);
        let ao_staging = readback("p74 ao staging", ao_bytes);
        let copy_ao = |staging: &wgpu::Buffer| {
            let mut encoder =
                device.create_command_encoder(&wgpu::CommandEncoderDescriptor { label: None });
            encoder.copy_buffer_to_buffer(&ao_buffer, 0, staging, 0, ao_bytes);
            queue.submit(Some(encoder.finish()));
            wait();
        };
        copy_ao(&ao_staging);
        let ao_field_image = read_f32(&ao_staging, ao_bytes);
        let diag = {
            let staging = readback("p74 diag staging", 16);
            let mut encoder =
                device.create_command_encoder(&wgpu::CommandEncoderDescriptor { label: None });
            encoder.copy_buffer_to_buffer(&diag_buffer, 0, &staging, 0, 16);
            queue.submit(Some(encoder.finish()));
            wait();
            let raw = read_bytes(device, queue, &staging, 16).expect("diag");
            raw.as_chunks::<4>()
                .0
                .iter()
                .map(|w| u32::from_le_bytes(*w))
                .collect::<Vec<u32>>()
        };
        let _ = run_ao(&ssao_pipeline, "ao_ssao record", false);
        copy_ao(&ao_staging);
        let ao_ssao_image = read_f32(&ao_staging, ao_bytes);

        assert_eq!(
            diag[2] as usize, hit_pixels,
            "the shader visited {} hit pixels and the CPU classifier found {hit_pixels}, so the \
             two disagree about the G-buffer and every per-set number below is over the wrong \
             population",
            diag[2]
        );

        // ── the offline reference, on the sampled pixels ────────────────────
        let sample_of = |set: &[usize]| -> Vec<usize> {
            let stride = (set.len() / SAMPLE_CAP).max(1);
            set.iter()
                .step_by(stride)
                .copied()
                .take(SAMPLE_CAP)
                .collect()
        };
        /// `(mae_field, mae_ssao, mae_floor, ref_mean, field_mean, ssao_mean, sampled)`.
        type Set = (f64, f64, f64, f64, f64, f64, usize);
        let mae_for = |set: &[usize]| -> Set {
            let chosen = sample_of(set);
            let mut field_error = 0.0f64;
            let mut ssao_error = 0.0f64;
            let mut floor_error = 0.0f64;
            let mut reference_sum = 0.0f64;
            let mut field_sum = 0.0f64;
            let mut ssao_sum = 0.0f64;
            for &i in &chosen {
                let p = [gbuf_pos[i * 4], gbuf_pos[i * 4 + 1], gbuf_pos[i * 4 + 2]];
                let n = normalise([gbuf_nrm[i * 4], gbuf_nrm[i * 4 + 1], gbuf_nrm[i * 4 + 2]]);
                let reference = ao_by_mesh(&scene.trimesh, p, n, &dirs_ref, radius, lift);
                let floor = ao_by_mesh(&scene.trimesh, p, n, &dirs8, radius, lift);
                reference_sum += f64::from(reference);
                field_sum += f64::from(ao_field_image[i]);
                ssao_sum += f64::from(ao_ssao_image[i]);
                field_error += f64::from((ao_field_image[i] - reference).abs());
                ssao_error += f64::from((ao_ssao_image[i] - reference).abs());
                floor_error += f64::from((floor - reference).abs());
            }
            let n = chosen.len() as f64;
            (
                field_error / n,
                ssao_error / n,
                floor_error / n,
                reference_sum / n,
                field_sum / n,
                ssao_sum / n,
                chosen.len(),
            )
        };
        let reference_started = Instant::now();
        let (
            seam_mae_field,
            seam_mae_ssao,
            seam_mae_floor,
            seam_reference,
            seam_field_mean,
            seam_ssao_mean,
            seam_sampled,
        ) = mae_for(&seam_set);
        let (
            halo_mae_field,
            halo_mae_ssao,
            halo_mae_floor,
            halo_reference,
            halo_field_mean,
            halo_ssao_mean,
            halo_sampled,
        ) = mae_for(&silhouette_set);
        // ── the control set the registered premise needs ────────────────────
        //
        // **"SSAO's seam and halo error is non-zero on both fields" is only a
        // claim about seams if the seam set differs from the rest of the
        // image.** A method with a uniform bias posts a non-zero MAE on any
        // subset, and scoring that as a seam artefact would be reading a global
        // error as a local one. So the same three estimators also run over a
        // stride sample of *every* visible pixel, and the comparison
        // `seam_mae_* against all_mae_*` is what says whether the seam is a
        // discriminator or a coincidence.
        let all_set: Vec<usize> = (0..count).filter(|&i| hit[i]).collect();
        let (
            all_mae_field,
            all_mae_ssao,
            all_mae_floor,
            all_reference,
            all_field_mean,
            all_ssao_mean,
            all_sampled,
        ) = mae_for(&all_set);
        let reference_ms = reference_started.elapsed().as_nanos() as f64 / 1e6;

        // ── the gate the first run did not have ─────────────────────────────
        //
        // **A mean absolute difference against a reference that is identically
        // zero is not a measurement.** The first run at `AO_RADIUS_CELLS = 8`
        // reported `seam_mae_field` and `halo_mae_field` as exactly 0.0000 on
        // `gyroid`, and the reason was that the 512-ray reference was 0.0000
        // too: with no occlusion in the scene, a perfect method and a broken one
        // both score zero and C2 cannot be refuted in either direction. So the
        // reference's own magnitude is a precondition, asserted here, and the
        // bar is one 512-ray quantum per pixel so it cannot be met by rounding.
        let quantum = 1.0 / f64::from(REF_RAYS);
        for (label, mean) in [("seam", seam_reference), ("silhouette", halo_reference)] {
            assert!(
                mean > quantum,
                "VOID: the {REF_RAYS}-ray reference finds only {mean:.6} occlusion on \
                 {field_name}'s {label} pixels, under one ray's quantum of {quantum:.6}. Every \
                 MAE below would be a difference from zero and C2 would be answered by a \
                 fixture with nothing in it (M-44)."
            );
        }

        let c1_holds = ao_ms_field < 2.0;
        // C2 as registered: falsified by SSAO matching the field arm. "Matching"
        // is decided against the floor an ideal 8-ray method would post, because
        // a difference smaller than the discretisation residue is not a
        // difference between the methods.
        let c2_holds = seam_mae_ssao > seam_mae_field && halo_mae_ssao > halo_mae_field;
        // ── C3's predicate is the registered one, and the first draft's was not ──
        //
        // The first version scored C3 on `resident_ms < rebuild_ms * 0.05`, a
        // threshold this harness invented; it made the verdict a function of how
        // much of the frame the surface happened to cover, which is nothing to do
        // with the clause. The registration's falsifier is *"a non-zero resident
        // cost, which would mean the field is not reusable in the form the tracer
        // wants and there is a conversion nobody costed"* - so the predicate is
        // whether the resident field is traceable as it stands.
        //
        // A sphere trace is sound exactly when the field is 1-Lipschitz. The
        // resident field here is the crate's own `Sdf`, sampled; `raw_lipschitz`
        // is its largest per-cell difference over the cell size. Above 1 the
        // march overshoots and the jump flood is a **mandatory conversion**, and
        // `conversion_bytes` and `conversion_ms` are what it costs.
        //
        // `sphere_lipschitz` is the positive arm: an exact SDF on the same grid
        // through the same function must come out at 1, or the instrument is
        // measuring the grid rather than the field.
        let c3_holds = raw_lipschitz <= 1.0;
        let conversion_bytes = distance_bytes * 2;
        let conversion_ms = flood_build_ms + distance_upload_ms;

        // ── coverage, because C1's clock only ran where there was a surface ──
        //
        // A background pixel returns from `ao_field` before it traces anything,
        // so `ao_ms_field` is the cost of `hit_pixels` traces and not of
        // 2,073,600 of them. The scaled figure is what a frame filled edge to
        // edge with surface would cost - an over-estimate, since it charges the
        // fixed dispatch overhead twice - and it is on the row so the measured
        // number cannot be read as a full-coverage frame.
        let hit_fraction = hit_pixels as f64 / pixels as f64;
        let ao_ms_field_scaled = ao_ms_field / hit_fraction;
        let ao_ms_ssao_scaled = ao_ms_ssao / hit_fraction;
        let marches = hit_pixels as u64 * u64::from(RAYS);
        let exhausted_fraction = f64::from(diag[3]) / marches as f64;

        println!(
            "  ao_field {ao_ms_field:.4} ms (gpu span over {hit_pixels} traced pixels, \
             {:.2}% coverage -> {ao_ms_field_scaled:.4} ms at full coverage), ao_ssao \
             {ao_ms_ssao:.4} ms; C1 {}",
            100.0 * hit_fraction,
            if c1_holds { "HOLDS" } else { "FALSIFIED" }
        );
        println!(
            "  seam  field {seam_mae_field:.4}  ssao {seam_mae_ssao:.4}  ideal-8-ray floor \
             {seam_mae_floor:.4}  (means: ref {seam_reference:.4} field {seam_field_mean:.4} \
             ssao {seam_ssao_mean:.4}, {seam_sampled} px)"
        );
        println!(
            "  halo  field {halo_mae_field:.4}  ssao {halo_mae_ssao:.4}  ideal-8-ray floor \
             {halo_mae_floor:.4}  (means: ref {halo_reference:.4} field {halo_field_mean:.4} \
             ssao {halo_ssao_mean:.4}, {halo_sampled} px)"
        );
        println!(
            "  ALL   field {all_mae_field:.4}  ssao {all_mae_ssao:.4}  ideal-8-ray floor \
             {all_mae_floor:.4}  (means: ref {all_reference:.4} field {all_field_mean:.4} \
             ssao {all_ssao_mean:.4}, {all_sampled} px)  <- the control the seam claim needs"
        );
        println!(
            "  resident {resident_ms:.4} ms, rebuild {rebuild_ms:.4} ms, difference {:.4} ms; \
             raw field Lipschitz {raw_lipschitz:.4} (sphere control {sphere_lipschitz:.4}), \
             mandatory conversion {conversion_ms:.3} ms over {conversion_bytes} bytes; C3 {}",
            rebuild_ms - resident_ms,
            if c3_holds { "HOLDS" } else { "FALSIFIED" }
        );
        println!(
            "  field starts inside {} of {hit_pixels}, within eps {}; {} of {marches} marches \
             ran out of steps ({:.2}%); reference took {reference_ms:.0} ms\n",
            diag[0],
            diag[1],
            diag[3],
            100.0 * exhausted_fraction
        );

        rows.push(vec![
            ("field", field_name.to_string()),
            ("resolution", format!("{WIDTH}x{HEIGHT}")),
            ("rays_per_pixel", RAYS.to_string()),
            ("march_steps", STEPS.to_string()),
            ("ao_ms_field", format!("{ao_ms_field:.6}")),
            ("ao_ms_ssao", format!("{ao_ms_ssao:.6}")),
            ("seam_pixels", seam_set.len().to_string()),
            ("silhouette_pixels", silhouette_set.len().to_string()),
            ("seam_mae_field", format!("{seam_mae_field:.6}")),
            ("seam_mae_ssao", format!("{seam_mae_ssao:.6}")),
            ("halo_mae_field", format!("{halo_mae_field:.6}")),
            ("halo_mae_ssao", format!("{halo_mae_ssao:.6}")),
            ("resident_ms", format!("{resident_ms:.6}")),
            ("rebuild_ms", format!("{rebuild_ms:.6}")),
            ("c1_holds", c1_holds.to_string()),
            ("c2_holds", c2_holds.to_string()),
            ("c3_holds", c3_holds.to_string()),
            ("adapter", adapter.clone()),
            // ── extras ─────────────────────────────────────────────────────
            ("backend", backend.clone()),
            ("clock", "gpu-timestamp-ao/wall-c3".to_string()),
            ("samples_per_axis", SAMPLES.to_string()),
            ("chunk_cells", CHUNK_CELLS.to_string()),
            ("vertices", scene.vertex_count.to_string()),
            ("triangles", scene.triangle_count.to_string()),
            ("hit_pixels", hit_pixels.to_string()),
            ("seam_pixels_screen", seam_screen.to_string()),
            ("seam_pixels_one_chunk", seam_one_chunk.to_string()),
            ("silhouette_pixels_flat", silhouette_flat.to_string()),
            ("seam_mae_floor", format!("{seam_mae_floor:.6}")),
            ("halo_mae_floor", format!("{halo_mae_floor:.6}")),
            ("seam_reference_mean", format!("{seam_reference:.6}")),
            ("halo_reference_mean", format!("{halo_reference:.6}")),
            ("seam_sampled", seam_sampled.to_string()),
            ("halo_sampled", halo_sampled.to_string()),
            ("reference_rays", REF_RAYS.to_string()),
            ("field_start_inside", diag[0].to_string()),
            ("field_start_within_eps", diag[1].to_string()),
            ("gbuffer_ms", format!("{gbuffer_ms:.6}")),
            ("field_upload_ms", format!("{field_upload_ms:.6}")),
            ("flood_build_ms", format!("{flood_build_ms:.6}")),
            ("distance_upload_ms", format!("{distance_upload_ms:.6}")),
            ("distance_bytes", distance_bytes.to_string()),
            ("raw_lipschitz", format!("{raw_lipschitz:.6}")),
            ("flood_lipschitz", format!("{flood_lipschitz:.6}")),
            ("queries_per_frame", queries.to_string()),
            ("ao_radius_cells", format!("{AO_RADIUS_CELLS:.3}")),
            ("eps_cells", format!("{EPS_CELLS:.3}")),
            ("lift_cells", format!("{LIFT_CELLS:.3}")),
            ("hit_fraction", format!("{hit_fraction:.6}")),
            (
                "ao_ms_field_full_coverage",
                format!("{ao_ms_field_scaled:.6}"),
            ),
            (
                "ao_ms_ssao_full_coverage",
                format!("{ao_ms_ssao_scaled:.6}"),
            ),
            ("marches", marches.to_string()),
            ("march_exhausted", diag[3].to_string()),
            (
                "march_exhausted_fraction",
                format!("{exhausted_fraction:.6}"),
            ),
            ("seam_field_mean", format!("{seam_field_mean:.6}")),
            ("seam_ssao_mean", format!("{seam_ssao_mean:.6}")),
            ("halo_field_mean", format!("{halo_field_mean:.6}")),
            ("halo_ssao_mean", format!("{halo_ssao_mean:.6}")),
            ("reference_ms", format!("{reference_ms:.3}")),
            ("all_mae_field", format!("{all_mae_field:.6}")),
            ("all_mae_ssao", format!("{all_mae_ssao:.6}")),
            ("all_mae_floor", format!("{all_mae_floor:.6}")),
            ("all_reference_mean", format!("{all_reference:.6}")),
            ("all_field_mean", format!("{all_field_mean:.6}")),
            ("all_ssao_mean", format!("{all_ssao_mean:.6}")),
            ("all_sampled", all_sampled.to_string()),
            ("sphere_lipschitz", format!("{sphere_lipschitz:.6}")),
            ("conversion_bytes", conversion_bytes.to_string()),
            ("conversion_ms", format!("{conversion_ms:.6}")),
            (
                "rebuild_minus_resident_ms",
                format!("{:.6}", rebuild_ms - resident_ms),
            ),
            (
                "ns_per_traced_pixel",
                format!("{:.4}", ao_ms_field * 1e6 / hit_pixels as f64),
            ),
            (
                "c1_holds_full_coverage",
                (ao_ms_field_scaled < 2.0).to_string(),
            ),
        ]);
    }

    // ── cross-field controls ────────────────────────────────────────────────
    assert_eq!(rows.len(), 2, "both registered fields must produce a row");

    common::experiment::run(prereg, |run| {
        for row in rows {
            run.record(&row);
        }
    });
}
