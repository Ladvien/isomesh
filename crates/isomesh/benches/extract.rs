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
use isomesh::dual_contouring::DualContouring;
use isomesh::fields::{BoxExact, ReferenceField, Sphere, Torus, capped_gyroid};
use isomesh::marching_cubes::{FaceAmbiguity, MarchingCubes};
use isomesh::marching_tetrahedra::MarchingTetrahedra;
use isomesh::surface_nets::SurfaceNets;
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

/// Marching Cubes with the asymptotic decider.
///
/// Identical to [`bench_mc`] but for the one setter, deliberately: the pair is
/// the measurement, and anything else that differed between them would be in
/// the difference.
fn bench_mc33<R, F>(c: &mut Criterion, label: &str, field: F, samples: u32)
where
    R: Real,
    F: ReferenceField + Sdf<Scalar = R>,
{
    let (shape, origin, h) = common::grid(&field, samples);
    let mut mc = MarchingCubes::<R>::new();
    mc.set_face_ambiguity(FaceAmbiguity::AsymptoticDecider);
    let mut out = MeshBuffer::<R>::new();
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

/// Marching Tetrahedra, for M-001's shootout row.
fn bench_mt<R, F>(c: &mut Criterion, label: &str, field: F, samples: u32)
where
    R: Real,
    F: ReferenceField + Sdf<Scalar = R>,
{
    let (shape, origin, h) = common::grid(&field, samples);
    let mut mt = MarchingTetrahedra::<R>::new();
    let mut out = MeshBuffer::<R>::new();
    mt.extract(&field, &shape, origin, h, &mut out)
        .expect("extraction");

    c.bench_function(label, |b| {
        b.iter(|| {
            out.reset();
            mt.extract(&field, &shape, origin, h, &mut out)
                .expect("extraction");
            black_box(out.triangle_count())
        });
    });
}

fn algorithms(c: &mut Criterion) {
    for n in COMPARISON_SAMPLES {
        bench_mc(
            c,
            &format!("marching_cubes/sphere/f32/{n}"),
            Sphere::<f32>::canonical(),
            n,
        );
        bench_sn(
            c,
            &format!("surface_nets/sphere/f32/{n}"),
            Sphere::<f32>::canonical(),
            n,
        );
        bench_mc(
            c,
            &format!("marching_cubes/torus/f32/{n}"),
            Torus::<f32>::canonical(),
            n,
        );
        bench_sn(
            c,
            &format!("surface_nets/torus/f32/{n}"),
            Torus::<f32>::canonical(),
            n,
        );
        bench_mc(
            c,
            &format!("marching_cubes/box_exact/f32/{n}"),
            BoxExact::<f32>::canonical(),
            n,
        );
        bench_sn(
            c,
            &format!("surface_nets/box_exact/f32/{n}"),
            BoxExact::<f32>::canonical(),
            n,
        );
        bench_mt(
            c,
            &format!("marching_tetrahedra/sphere/f32/{n}"),
            Sphere::<f32>::canonical(),
            n,
        );
        bench_mt(
            c,
            &format!("marching_tetrahedra/torus/f32/{n}"),
            Torus::<f32>::canonical(),
            n,
        );
        bench_mt(
            c,
            &format!("marching_tetrahedra/box_exact/f32/{n}"),
            BoxExact::<f32>::canonical(),
            n,
        );
        bench_dc(
            c,
            &format!("dual_contouring/sphere/f32/{n}"),
            Sphere::<f32>::canonical(),
            n,
        );
        bench_dc(
            c,
            &format!("dual_contouring/torus/f32/{n}"),
            Torus::<f32>::canonical(),
            n,
        );
        bench_dc(
            c,
            &format!("dual_contouring/box_exact/f32/{n}"),
            BoxExact::<f32>::canonical(),
            n,
        );
    }
}

/// What the asymptotic decider costs, against the identical extraction without
/// it.
///
/// Two fields on purpose. On `sphere` the census in `mc/tests.rs` finds **no
/// ambiguous face at all**, so the difference there is the price of asking —
/// one table lookup and a branch per surface cell. On `gyroid` the rule fires on
/// about half a percent of surface cells, so that pair carries the price of
/// answering as well: building the cell's triangulation at run time instead of
/// reading it from the compile-time table.
fn decider(c: &mut Criterion) {
    for n in COMPARISON_SAMPLES {
        bench_mc(
            c,
            &format!("decider/marching_cubes/sphere/f32/{n}"),
            Sphere::<f32>::canonical(),
            n,
        );
        bench_mc33(
            c,
            &format!("decider/marching_cubes+decider/sphere/f32/{n}"),
            Sphere::<f32>::canonical(),
            n,
        );
        bench_mc(
            c,
            &format!("decider/marching_cubes/gyroid/f32/{n}"),
            capped_gyroid::<f32>(),
            n,
        );
        bench_mc33(
            c,
            &format!("decider/marching_cubes+decider/gyroid/f32/{n}"),
            capped_gyroid::<f32>(),
            n,
        );
    }
}

/// The same field and grid at both widths. Any difference is the cost of `f64`
/// on a path with no matrix solve in it.
fn precision(c: &mut Criterion) {
    bench_mc(
        c,
        "precision/marching_cubes/sphere/f32",
        Sphere::<f32>::canonical(),
        65,
    );
    bench_mc(
        c,
        "precision/marching_cubes/sphere/f64",
        Sphere::<f64>::canonical(),
        65,
    );
    bench_sn(
        c,
        "precision/surface_nets/sphere/f32",
        Sphere::<f32>::canonical(),
        65,
    );
    bench_sn(
        c,
        "precision/surface_nets/sphere/f64",
        Sphere::<f64>::canonical(),
        65,
    );
}

criterion_group!(benches, algorithms, decider, precision);
criterion_main!(benches);
