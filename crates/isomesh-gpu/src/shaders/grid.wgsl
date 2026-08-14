// The shader side of `GridParams`. Every kernel in this crate includes this.
//
// The layout must agree with `GridParams::to_std140` byte for byte: two vec4s,
// 32 bytes, x/y/z used and w carrying padding or the cell size. Two vec4s
// rather than a struct of scalars so std140 and std430 agree and there is no
// conditional padding rule to get wrong.

struct GridParams {
    // xyz: samples per axis. w: padding, always zero.
    samples: vec4<u32>,
    // xyz: world position of sample [0,0,0]. w: cell size.
    placement: vec4<f32>,
}

// Samples per axis.
fn grid_samples(g: GridParams) -> vec3<u32> {
    return g.samples.xyz;
}

// Total samples in the grid.
fn grid_sample_count(g: GridParams) -> u32 {
    return g.samples.x * g.samples.y * g.samples.z;
}

// Cells per axis: one fewer than the samples, because a cell spans two.
fn grid_cells(g: GridParams) -> vec3<u32> {
    return g.samples.xyz - vec3<u32>(1u, 1u, 1u);
}

// Flat index of a sample. x varies fastest: i = x + y*sx + z*sx*sy.
//
// This ordering is a convention shared with the CPU side and a mismatch
// produces a mesh that looks plausible and is transposed, so it lives in one
// function rather than being written out at each call site.
fn grid_index(g: GridParams, at: vec3<u32>) -> u32 {
    return at.x + at.y * g.samples.x + at.z * g.samples.x * g.samples.y;
}

// World position of a sample.
//
// The index is MULTIPLIED by the cell size, never accumulated. isomesh's M-70
// and M-73 both record cracks caused by the other choice: at a spacing that is
// not a power of two, (origin + h*i) + h differs from origin + h*(i + 1) in the
// last bit, and two chunks that disagree there leave a hole a weld cannot close.
fn grid_position(g: GridParams, at: vec3<u32>) -> vec3<f32> {
    return g.placement.xyz + g.placement.w * vec3<f32>(at);
}

// The 3-D index of a flat sample index -- the inverse of grid_index.
//
// Shared here rather than rewritten per kernel: a transposition between the
// forward and inverse mapping produces a field that is sampled correctly and
// stored rotated, which looks like a meshing bug.
fn grid_sample_at(g: GridParams, flat: u32) -> vec3<u32> {
    return vec3<u32>(
        flat % g.samples.x,
        (flat / g.samples.x) % g.samples.y,
        flat / (g.samples.x * g.samples.y),
    );
}

// Whether a sample index is inside the grid.
fn grid_contains(g: GridParams, at: vec3<u32>) -> bool {
    return all(at < g.samples.xyz);
}
