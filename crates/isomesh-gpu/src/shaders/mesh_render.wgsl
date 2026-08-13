// Drawing isomesh's compute output with mesh shaders, without ever reading it
// back.
//
// A mesh shader consumes the very buffers `marching_cubes.wgsl` wrote --
// positions and normals, in place, on the GPU -- so *drawing* needs no
// read-back, no vertex buffer and no index buffer.
//
// Be precise about what that saves, because an earlier version of this comment
// was not. The extraction still reads the per-cell counts back to prefix-sum
// them on the CPU, so the pipeline as a whole is not read-back-free. Measured
// (M-149): the geometry read-back this removes is **6.7% of the GPU path at
// 129^3** and is the *smallest* of the three data-movement costs -- counts
// read-back 1.97 ms, CPU prefix sum 3.27 ms, upload 8.63 ms, geometry
// read-back 1.00 ms. Removing all of it needs a GPU scan and an indirect
// draw, which is GPU-010.
//
// The input is the triangle SOUP those kernels emit: three vertices per
// triangle, in cell order, `array<f32>` with x fastest. No indices, because
// there are none to have.
//
// Requires `EXPERIMENTAL_MESH_SHADER`, which is Vulkan-only for WGSL: wgpu's own
// source says "naga is only supported on vulkan; on other platforms you will
// have to use passthrough shaders" (V-23). A caller must check the capability
// before building a pipeline from this, and refuse loudly if it is absent.

enable wgpu_mesh_shader;

// Triangles emitted per mesh workgroup, and therefore 3x that many vertices.
//
// 32 is chosen to sit inside the smallest limits worth targeting rather than to
// be optimal: `maxMeshOutputVertices` is commonly 256, so 96 vertices and 32
// primitives clears it with room. Tuning this is a measurement nobody has taken
// yet, and picking a number that merely fits is honest until they do.
const BATCH: u32 = 32u;

struct Camera {
    view_proj: mat4x4<f32>,
}

struct Draw {
    // Triangles in the soup. The dispatch is ceil(triangles / BATCH)
    // workgroups, so the last one is usually partial.
    triangle_count: u32,
}

@group(0) @binding(0) var<uniform> camera: Camera;
@group(0) @binding(1) var<uniform> draw: Draw;
@group(0) @binding(2) var<storage, read> positions: array<f32>;
@group(0) @binding(3) var<storage, read> normals: array<f32>;

struct Vertex {
    @builtin(position) clip: vec4<f32>,
    @location(0) normal: vec3<f32>,
}

struct Primitive {
    @builtin(triangle_indices) indices: vec3<u32>,
}

// The mesh output block. Its shape is not a matter of taste: naga requires a
// `workgroup` struct whose members carry exactly these builtins, with the
// vertex and primitive arrays sized by constants
// (`proc/mod.rs:723`, `valid/interface.rs:1531`).
struct MeshOut {
    @builtin(vertex_count) vertex_count: u32,
    @builtin(primitive_count) primitive_count: u32,
    @builtin(vertices) vertices: array<Vertex, 96>,
    @builtin(primitives) primitives: array<Primitive, 32>,
}

var<workgroup> mesh_out: MeshOut;

// Two near-identical readers rather than one taking a pointer.
//
// WGSL does not allow a function parameter in the `storage` address space --
// naga rejects it as `InvalidArgumentPointerSpace`, and the first version of
// this file was written that way and refused to validate. The duplication is
// the language's, not a preference.
fn read_position(index: u32) -> vec3<f32> {
    let base = index * 3u;
    return vec3<f32>(positions[base], positions[base + 1u], positions[base + 2u]);
}

fn read_normal(index: u32) -> vec3<f32> {
    let base = index * 3u;
    return vec3<f32>(normals[base], normals[base + 1u], normals[base + 2u]);
}

@mesh(mesh_out) @workgroup_size(32)
fn draw_mesh(
    @builtin(workgroup_id) group: vec3<u32>,
    @builtin(local_invocation_index) local: u32,
) {
    let first = group.x * BATCH;
    // The last workgroup is partial whenever the triangle count is not a
    // multiple of BATCH, which is the usual case.
    var live = draw.triangle_count - first;
    if (live > BATCH) {
        live = BATCH;
    }

    // One invocation sets the counts. Every thread writing the same value would
    // be a race that happens to agree, and "happens to agree" is not a property
    // worth relying on across drivers.
    if (local == 0u) {
        mesh_out.vertex_count = live * 3u;
        mesh_out.primitive_count = live;
    }
    workgroupBarrier();

    if (local >= live) {
        return;
    }

    let triangle = first + local;
    let out_base = local * 3u;
    for (var k = 0u; k < 3u; k = k + 1u) {
        let vertex = triangle * 3u + k;
        let world = read_position(vertex);
        mesh_out.vertices[out_base + k].clip = camera.view_proj * vec4<f32>(world, 1.0);
        mesh_out.vertices[out_base + k].normal = read_normal(vertex);
    }
    mesh_out.primitives[local].indices =
        vec3<u32>(out_base, out_base + 1u, out_base + 2u);
}

// A deliberately plain lambert, because the point of the demo is the pipeline
// rather than the shading.
@fragment
fn shade(@location(0) normal: vec3<f32>) -> @location(0) vec4<f32> {
    let light = normalize(vec3<f32>(0.4, 0.8, 0.35));
    let lambert = max(dot(normalize(normal), light), 0.0);
    let colour = vec3<f32>(0.72, 0.62, 0.50) * (0.25 + 0.75 * lambert);
    return vec4<f32>(colour, 1.0);
}
