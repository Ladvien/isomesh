// Jump flooding, 3D. Rong & Tan, "Jump flooding in GPU with applications to
// Voronoi diagram and distance transform", I3D 2006.
//
// Ticket: S-005. Three kernels: seed the surface, flood at halving strides,
// resolve to a signed distance.
//
// # The seed is a crossing, not a sample
//
// The textbook seeds each boundary *sample* with its own position, which
// quantises every distance to the grid and puts a floor of half a cell on the
// error. This seeds the linearly interpolated **crossing point** on each cut
// edge instead -- the same sub-cell information `construct::signed_distance_
// field_swept` uses on the CPU -- so the comparison against S-001 measures the
// jump flood's approximation and not a seeding handicap.

#include <grid>

// Which stride this dispatch is at. One buffer per step, because a uniform
// cannot be rewritten between dispatches inside a single encoder.
// A single vec4 rather than `u32` plus padding: std140 aligns a trailing
// vec3<u32> to 16, which makes the obvious `{ stride: u32, pad: vec3<u32> }`
// **32** bytes and not 16. wgpu catches the mismatch, but only at dispatch.
struct FloodStep {
    // x: offset in samples, halving each pass: n/2, n/4, ... 1. yzw unused.
    stride: vec4<u32>,
}

@group(0) @binding(0) var<uniform> grid: GridParams;
@group(0) @binding(1) var<uniform> step: FloodStep;
// The scalar field. Read for its signs, never written.
@group(0) @binding(2) var<storage, read> field: array<f32>;
// Nearest known surface point per sample. xyz is a world position; w is 1.0
// when the entry is valid and 0.0 when nothing has reached this sample yet.
@group(0) @binding(3) var<storage, read_write> src: array<vec4<f32>>;
@group(0) @binding(4) var<storage, read_write> dst: array<vec4<f32>>;
// The signed distance, written only by `resolve`.
@group(0) @binding(5) var<storage, read_write> out: array<f32>;

// Threads per workgroup. 64 is one wavefront on AMD and two warps on NVIDIA;
// the kernel is memory-bound either way, so the number is not delicate.
const WORKGROUP: u32 = 64u;

// Recover the 3D index of a flat one. The inverse of `grid_index`, and it must
// stay the inverse -- x fastest.
fn unflatten(g: GridParams, i: u32) -> vec3<u32> {
    let sx = g.samples.x;
    let sy = g.samples.y;
    return vec3<u32>(i % sx, (i / sx) % sy, i / (sx * sy));
}

// Seed pass: every sample adjacent to a sign change adopts the nearest
// interpolated crossing among its six edges. Everything else is left invalid.
@compute @workgroup_size(WORKGROUP)
fn seed(@builtin(global_invocation_id) id: vec3<u32>) {
    let i = id.x;
    if i >= grid_sample_count(grid) {
        return;
    }
    let at = unflatten(grid, i);
    let s = grid_samples(grid);
    let here = field[i];
    let inside = here < 0.0;
    let p = grid_position(grid, at);
    let h = grid.placement.w;

    var best = 1e30;
    var found = vec3<f32>(0.0, 0.0, 0.0);

    // Six axis neighbours. Written out rather than looped because WGSL has no
    // signed offset that stays in bounds without a branch anyway.
    for (var axis = 0u; axis < 3u; axis = axis + 1u) {
        for (var dir = 0u; dir < 2u; dir = dir + 1u) {
            var d = vec3<i32>(0, 0, 0);
            let delta = select(1, -1, dir == 0u);
            if axis == 0u { d.x = delta; } else if axis == 1u { d.y = delta; } else { d.z = delta; }
            let n = vec3<i32>(at) + d;
            if any(n < vec3<i32>(0, 0, 0)) || any(n >= vec3<i32>(s)) {
                continue;
            }
            let j = grid_index(grid, vec3<u32>(n));
            let there = field[j];
            if (there < 0.0) == inside {
                continue;
            }
            let denom = here - there;
            if denom == 0.0 {
                continue;
            }
            // Fraction of the edge from this sample to the crossing.
            let t = abs(here / denom);
            let dist = t * h;
            if dist < best {
                best = dist;
                found = p + vec3<f32>(d) * (t * h);
            }
        }
    }

    if best < 1e30 {
        dst[i] = vec4<f32>(found, 1.0);
    } else {
        dst[i] = vec4<f32>(0.0, 0.0, 0.0, 0.0);
    }
}

// Flood pass: adopt the best seed among this sample's own entry and its 26
// neighbours at the current stride.
//
// **This is where the approximation lives.** A sample only ever sees seeds that
// reached one of 27 fixed offsets, so a Voronoi region whose owner is not on
// that lattice at any stride is missed. Rong & Tan report the error as rare and
// small rather than absent; S-005's acceptance is to measure it here rather
// than repeat that.
@compute @workgroup_size(WORKGROUP)
fn flood(@builtin(global_invocation_id) id: vec3<u32>) {
    let i = id.x;
    if i >= grid_sample_count(grid) {
        return;
    }
    let at = vec3<i32>(unflatten(grid, i));
    let s = vec3<i32>(grid_samples(grid));
    let p = grid_position(grid, vec3<u32>(at));
    let k = i32(step.stride.x);

    var best = src[i];
    var best_d = select(1e30, distance(p, best.xyz), best.w > 0.5);

    for (var dz = -1; dz <= 1; dz = dz + 1) {
        for (var dy = -1; dy <= 1; dy = dy + 1) {
            for (var dx = -1; dx <= 1; dx = dx + 1) {
                let n = at + vec3<i32>(dx, dy, dz) * k;
                if any(n < vec3<i32>(0, 0, 0)) || any(n >= s) {
                    continue;
                }
                let cand = src[grid_index(grid, vec3<u32>(n))];
                if cand.w < 0.5 {
                    continue;
                }
                let d = distance(p, cand.xyz);
                if d < best_d {
                    best_d = d;
                    best = cand;
                }
            }
        }
    }

    dst[i] = best;
}

// Resolve: distance to the adopted seed, signed by the field's own sign.
//
// The sign comes from the field rather than from the flood, because the flood
// carries no inside/outside information -- and taking it from the field is what
// makes this agree with the CPU constructors on which side of the surface a
// sample is, even where they disagree on how far.
@compute @workgroup_size(WORKGROUP)
fn resolve(@builtin(global_invocation_id) id: vec3<u32>) {
    let i = id.x;
    if i >= grid_sample_count(grid) {
        return;
    }
    let at = unflatten(grid, i);
    let p = grid_position(grid, at);
    let seed = src[i];
    // 1e30 rather than infinity, matching `construct::far` on the CPU: an
    // unreached sample must stay orderable and must not turn into a NaN the
    // first time something subtracts two of them.
    var d = 1e30;
    if seed.w > 0.5 {
        d = distance(p, seed.xyz);
    }
    out[i] = select(d, -d, field[i] < 0.0);
}
