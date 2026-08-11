//! Per-algorithm extraction benchmarks.
//!
//! These exist to notice a regression, not to answer "how fast is the
//! algorithm" — that question is answered by
//! [`resolution_sweep`](../resolution_sweep/index.html), because a single-grid
//! timing measures dispatch and setup as much as it measures the algorithm
//! (V-6, and the reason T-006 asks for the fixed cost separately).
//!
//! # What is timed
//!
//! The **re-mesh** path, not a cold first run: the extractor and the output
//! buffer are built once, outside the timed closure, and the buffer is `reset()`
//! between iterations. That is rule 6's contract and it is what the real
//! workload does — thousands of chunks through one buffer. Timing an allocation
//! per iteration would measure the allocator.
//!
//! Both scalar widths appear at one resolution, which is the cheapest available
//! evidence toward O-8: whether `f64` costs enough to matter on the paths that
//! do not involve a QEF.

// `criterion_group!` expands to a `pub fn` this file cannot attach a doc comment
// to, and the workspace denies `missing_docs`. A bench binary exports nothing,
// so the lint has no subject here.
#![allow(missing_docs)]

mod common;

use std::hint::black_box;

use criterion::{Criterion, criterion_group, criterion_main};
use isomesh::dc::DualContouring;
use isomesh::fields::{BoxExact, ReferenceField, Sphere, Torus};
use isomesh::mc::MarchingCubes;
use isomesh::sn::SurfaceNets;
use isomesh::{MeshBuffer, Real, Sdf};

/// Resolutions the per-field comparison runs at, in samples per axis.
const COMPARISON_SAMPLES: [u32; 2] = [33, 65];

fn bench_mc<R, F>(c: &mut Criterion, label: &str, field: F, samples: u32)
where
    R: Real,
    F: ReferenceField + Sdf<Scalar = R>,
{
    let (shape, origin, h) = common::grid(&field, samples);
    let mut mc = MarchingCubes::<R>::new();
    let mut out = MeshBuffer::<R>::new();
    // One extraction up front so the scratch buffers are already sized; the
    // benchmark is the steady state, not the first call.
    mc.extract(&field, &shape, origin, h, &mut out)
        .expect("extraction");

    c.bench_function(label, |b| {
        b.iter(|| {
            out.reset();
            mc.extract(&field, &shape, origin, h, &mut out)
                .expect("extraction");
            black_box(out.triangle_count())
        });
    });
}

fn bench_sn<R, F>(c: &mut Criterion, label: &str, field: F, samples: u32)
where
    R: Real,
    F: ReferenceField + Sdf<Scalar = R>,
{
    let (shape, origin, h) = common::grid(&field, samples);
    let mut sn = SurfaceNets::<R>::new();
    let mut out = MeshBuffer::<R>::new();
    sn.extract(&field, &shape, origin, h, &mut out)
        .expect("extraction");

    c.bench_function(label, |b| {
        b.iter(|| {
            out.reset();
            sn.extract(&field, &shape, origin, h, &mut out)
                .expect("extraction");
            black_box(out.triangle_count())
        });
    });
}

fn bench_dc<R, F>(c: &mut Criterion, label: &str, field: F, samples: u32)
where
    R: Real,
    F: ReferenceField + Sdf<Scalar = R>,
{
    let (shape, origin, h) = common::grid(&field, samples);
    let mut dc = DualContouring::<R>::new();
    let mut out = MeshBuffer::<R>::new();
    dc.extract(&field, &shape, origin, h, &mut out)
        .expect("extraction");

    c.bench_function(label, |b| {
        b.iter(|| {
            out.reset();
            dc.extract(&field, &shape, origin, h, &mut out)
                .expect("extraction");
            black_box(out.triangle_count())
        });
    });
}

fn algorithms(c: &mut Criterion) {
    for n in COMPARISON_SAMPLES {
        bench_mc(
            c,
            &format!("mc/sphere/f32/{n}"),
            Sphere::<f32>::canonical(),
            n,
        );
        bench_sn(
            c,
            &format!("sn/sphere/f32/{n}"),
            Sphere::<f32>::canonical(),
            n,
        );
        bench_mc(
            c,
            &format!("mc/torus/f32/{n}"),
            Torus::<f32>::canonical(),
            n,
        );
        bench_sn(
            c,
            &format!("sn/torus/f32/{n}"),
            Torus::<f32>::canonical(),
            n,
        );
        bench_mc(
            c,
            &format!("mc/box_exact/f32/{n}"),
            BoxExact::<f32>::canonical(),
            n,
        );
        bench_sn(
            c,
            &format!("sn/box_exact/f32/{n}"),
            BoxExact::<f32>::canonical(),
            n,
        );
        bench_dc(
            c,
            &format!("dc/sphere/f32/{n}"),
            Sphere::<f32>::canonical(),
            n,
        );
        bench_dc(
            c,
            &format!("dc/torus/f32/{n}"),
            Torus::<f32>::canonical(),
            n,
        );
        bench_dc(
            c,
            &format!("dc/box_exact/f32/{n}"),
            BoxExact::<f32>::canonical(),
            n,
        );
    }
}

/// The same field and grid at both widths. Any difference is the cost of `f64`
/// on a path with no matrix solve in it.
fn precision(c: &mut Criterion) {
    bench_mc(c, "precision/mc/sphere/f32", Sphere::<f32>::canonical(), 65);
    bench_mc(c, "precision/mc/sphere/f64", Sphere::<f64>::canonical(), 65);
    bench_sn(c, "precision/sn/sphere/f32", Sphere::<f32>::canonical(), 65);
    bench_sn(c, "precision/sn/sphere/f64", Sphere::<f64>::canonical(), 65);
}

criterion_group!(benches, algorithms, precision);
criterion_main!(benches);
