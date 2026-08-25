// Evaluating a reference field on the GPU, so the samples are never uploaded.
//
// After GPU-010a and GPU-012 the upload is 86% of the extraction path -- 7.04 ms
// of 8.15 at 129^3 -- and what remains of it is `field.sample()` running on the
// CPU. Writing the samples here removes the upload entirely: they are produced
// where they are consumed and nothing crosses the bus.
//
// # This is the cheapest version, deliberately
//
// Four reference fields, selected by a uniform. That is enough to answer "is the
// 86% actually recoverable" and to measure how far GPU arithmetic drifts from
// libm's, which is the question the general mechanism has to be designed
// around. Interpreting an edit log, or composing a consumer's own WGSL, is
// GPU-011b.
//
// # The formulas are transcribed from isomesh, not recalled
//
// Each is the same expression `crates/isomesh/src/fields/mod.rs` evaluates, in
// the same order, including `length(v)` written out as `sqrt(dot(v, v))` --
// WGSL's `length` builtin is free to compute it differently, and the point of
// this file is to measure the platform's arithmetic rather than a second
// formula's.

#include <grid>

const FIELD_SPHERE: u32 = 0u;
const FIELD_TORUS: u32 = 1u;
const FIELD_BOX_EXACT: u32 = 2u;
const FIELD_GYROID: u32 = 3u;
// The base is already in `base`, put there by the caller. GPU-011b's log
// interpreter with the base field taken out of the shader's vocabulary and
// handed over as data -- which is how `FieldSampler::fold_into` serves a
// consumer whose base is an arbitrary analytic SDF the four above cannot spell.
const FIELD_SAMPLED: u32 = 4u;

struct FieldSelect {
    // Which reference field is the base.
    id: u32,
    // Brushes in the log. Zero means the base field alone.
    brush_count: u32,
}

// One brush in the edit log, 64 bytes.
//
// Four vec4s rather than a packed struct, for the same reason GridParams is two:
// std140 and std430 then agree and there is no padding rule to get wrong. The
// slots mean different things per shape and the mapping is written out once,
// here and in `GpuBrush::to_std140`, because a disagreement between them
// produces a brush of the right kind in the wrong place.
//
//   header : [kind, op, _, _]
//   a      : sphere/box/torus centre, or capsule endpoint A .xyz ; .w = radius or major
//   b      : box half-extents, or capsule endpoint B .xyz ; .w = minor
//   join   : .x = smooth-add width k
struct Brush {
    header: vec4<u32>,
    a: vec4<f32>,
    b: vec4<f32>,
    join: vec4<f32>,
}

const SHAPE_SPHERE: u32 = 0u;
const SHAPE_BOX: u32 = 1u;
const SHAPE_CAPSULE: u32 = 2u;

const OP_ADD: u32 = 0u;
const OP_SUBTRACT: u32 = 1u;
const OP_SMOOTH_ADD: u32 = 2u;

@group(0) @binding(0) var<uniform> params: GridParams;
@group(0) @binding(1) var<uniform> select: FieldSelect;
@group(0) @binding(2) var<storage, read_write> samples: array<f32>;
@group(0) @binding(3) var<storage, read> brushes: array<Brush>;
// Read-only base samples, for `FIELD_SAMPLED`. Bound on every path because one
// kernel has one layout; the other four ids never read it.
@group(0) @binding(4) var<storage, read> base: array<f32>;

// isomesh's `vec3::length`: `dot` then `sqrt`, in that order.
fn iso_length(v: vec3<f32>) -> f32 {
    return sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
}

// Sphere: centre origin, radius 1. `length(p - c) - r`.
fn sphere(p: vec3<f32>) -> f32 {
    return iso_length(p) - 1.0;
}

// Torus: centre origin, major 1, minor 0.3.
fn torus(p: vec3<f32>) -> f32 {
    let s = sqrt(p.x * p.x + p.z * p.z);
    let q = vec2<f32>(s - 1.0, p.y);
    return sqrt(q.x * q.x + q.y * q.y) - 0.3;
}

// BoxExact: centre origin, half-extents [1, 1, 1].
fn box_exact(p: vec3<f32>) -> f32 {
    let q = abs(p) - vec3<f32>(1.0, 1.0, 1.0);
    let outside = vec3<f32>(max(q.x, 0.0), max(q.y, 0.0), max(q.z, 0.0));
    let inside = min(max(max(q.x, q.y), q.z), 0.0);
    return iso_length(outside) + inside;
}

// Gyroid: scale 1, iso 0. Triply periodic, and the only one here with
// transcendentals -- so the one where WGSL's accuracy bounds and libm's are
// expected to part company.
fn gyroid(p: vec3<f32>) -> f32 {
    return sin(p.x) * cos(p.y) + sin(p.y) * cos(p.z) + sin(p.z) * cos(p.x);
}

fn evaluate(id: u32, p: vec3<f32>) -> f32 {
    switch id {
        case 1u: { return torus(p); }
        case 2u: { return box_exact(p); }
        case 3u: { return gyroid(p); }
        // Sphere is the default rather than a listed case, because WGSL
        // requires a default and an unreachable one would be a branch nothing
        // ever takes.
        default: { return sphere(p); }
    }
}

// The three brush shapes, transcribed from `fields/mod.rs` and `brush.rs`.

fn shape_sphere(p: vec3<f32>, centre: vec3<f32>, radius: f32) -> f32 {
    return iso_length(p - centre) - radius;
}

fn shape_box(p: vec3<f32>, centre: vec3<f32>, half: vec3<f32>) -> f32 {
    let q = abs(p - centre) - half;
    let outside = vec3<f32>(max(q.x, 0.0), max(q.y, 0.0), max(q.z, 0.0));
    return iso_length(outside) + min(max(max(q.x, q.y), q.z), 0.0);
}

// Point-to-segment distance. A zero-length capsule is a sphere, which is the
// right answer rather than a degenerate case to reject -- `brush.rs` says so
// and the `denom > 0` guard is what implements it.
fn shape_capsule(p: vec3<f32>, a: vec3<f32>, b: vec3<f32>, radius: f32) -> f32 {
    let ab = b - a;
    let ap = p - a;
    let denom = ab.x * ab.x + ab.y * ab.y + ab.z * ab.z;
    var t = 0.0;
    if (denom > 0.0) {
        t = clamp((ap.x * ab.x + ap.y * ab.y + ap.z * ab.z) / denom, 0.0, 1.0);
    }
    return iso_length(ap - ab * t) - radius;
}

// isomesh's `brush::smooth_min`, expression for expression.
//
// A `k` of zero degenerates to an ordinary `min` rather than dividing by zero,
// and the parenthesisation is the CPU's: `(b + (a - b) * h) - k * h * (1 - h)`.
// Smooth-min is not associative and not even bit-commutative (M-38), so the
// order the CPU folds in is part of the answer and this walks the log the same
// way.
fn smooth_min(a: f32, b: f32, k: f32) -> f32 {
    if (k <= 0.0) {
        return min(a, b);
    }
    let h = clamp(0.5 + 0.5 * (b - a) / k, 0.0, 1.0);
    return (b + (a - b) * h) - k * h * (1.0 - h);
}

// isomesh's `brush::apply`.
fn apply_op(op: u32, field: f32, shape: f32, k: f32) -> f32 {
    switch op {
        case 1u: { return max(field, -shape); }
        case 2u: { return smooth_min(field, shape, k); }
        default: { return min(field, shape); }
    }
}

fn brush_shape(brush: Brush, p: vec3<f32>) -> f32 {
    switch brush.header.x {
        case 1u: { return shape_box(p, brush.a.xyz, brush.b.xyz); }
        case 2u: { return shape_capsule(p, brush.a.xyz, brush.b.xyz, brush.a.w); }
        default: { return shape_sphere(p, brush.a.xyz, brush.a.w); }
    }
}

@compute @workgroup_size(64)
fn sample_field(@builtin(global_invocation_id) gid: vec3<u32>) {
    let flat = gid.x;
    if (flat >= grid_sample_count(params)) {
        return;
    }
    let p = grid_position(params, grid_sample_at(params, flat));

    // Base field, then the log folded over it first to last -- the same order
    // `BrushStack::sample` walks, which is load-bearing because a mixed
    // add/subtract log does not commute and a smooth one is not even
    // associative (M-36..M-38).
    // `FIELD_SAMPLED` reads the base the caller uploaded. Each invocation owns
    // index `flat` and reads only `base[flat]`, so there is no
    // cross-invocation dependency and no barrier. An `if` rather than WGSL's
    // `select` builtin, because the uniform above is already named `select`.
    var value = 0.0;
    if (select.id == FIELD_SAMPLED) {
        value = base[flat];
    } else {
        value = evaluate(select.id, p);
    }
    for (var i = 0u; i < select.brush_count; i = i + 1u) {
        let brush = brushes[i];
        value = apply_op(brush.header.y, value, brush_shape(brush, p), brush.join.x);
    }
    samples[flat] = value;
}
