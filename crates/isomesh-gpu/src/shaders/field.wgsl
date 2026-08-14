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

struct FieldSelect {
    // Which reference field to evaluate.
    id: u32,
}

@group(0) @binding(0) var<uniform> params: GridParams;
@group(0) @binding(1) var<uniform> select: FieldSelect;
@group(0) @binding(2) var<storage, read_write> samples: array<f32>;

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

@compute @workgroup_size(64)
fn sample_field(@builtin(global_invocation_id) gid: vec3<u32>) {
    let flat = gid.x;
    if (flat >= grid_sample_count(params)) {
        return;
    }
    samples[flat] = evaluate(select.id, grid_position(params, grid_sample_at(params, flat)));
}
