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
use isomesh::fields::{BoxExact, ReferenceField, Sphere, Torus, capped_gyroid, noise_cavity};
use isomesh::manifold_dual_contouring::ManifoldDualContouring;
use isomesh::marching_cubes::{FaceAmbiguity, MarchingCubes};
use isomesh::marching_tetrahedra::MarchingTetrahedra;
use isomesh::subgrid::extract::SubgridMarchingTetrahedra;
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

/// Manifold Dual Contouring.
///
/// Paired with [`bench_dc`] on the same fields on purpose: the two differ only
/// in the per-cycle split, so the gap is what A-010's manifold guarantee costs.
fn bench_mdc<R, F>(c: &mut Criterion, label: &str, field: F, samples: u32)
where
    R: Real,
    F: ReferenceField + Sdf<Scalar = R>,
{
    let (shape, origin, h) = common::grid(&field, samples);
    let mut mdc = ManifoldDualContouring::<R>::new();
    let mut out = MeshBuffer::<R>::new();
    mdc.extract(&field, &shape, origin, h, &mut out)
        .expect("extraction");

    c.bench_function(label, |b| {
        b.iter(|| {
            out.reset();
            mdc.extract(&field, &shape, origin, h, &mut out)
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

/// Marching Cubes 33 with the interior rule as well as the face one.
///
/// Identical to [`bench_mc33`] but for the second setter, for the same reason:
/// the pair is the measurement.
fn bench_trilinear<R, F>(c: &mut Criterion, label: &str, field: F, samples: u32)
where
    R: Real,
    F: ReferenceField + Sdf<Scalar = R>,
{
    let (shape, origin, h) = common::grid(&field, samples);
    let mut mc = MarchingCubes::<R>::new();
    mc.set_face_ambiguity(FaceAmbiguity::AsymptoticDecider);
    mc.set_interior_ambiguity(isomesh::marching_cubes::InteriorAmbiguity::Trilinear);
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

/// Subgrid Marching Tetrahedra.
///
/// The one algorithm here whose cost is not a function of the grid alone: it
/// finds every root along every tetrahedron edge, so its work scales with
/// `SUBGRID_SAMPLES` as well, and both numbers belong in the label. There is no
/// published figure to compare against — the paper reports none — so this is the
/// first measurement of it anywhere, and it exists to be a baseline rather than
/// a boast.
fn bench_smt<R, F>(c: &mut Criterion, label: &str, field: F, samples: u32)
where
    R: Real,
    F: ReferenceField + Sdf<Scalar = R>,
{
    let (shape, origin, h) = common::grid(&field, samples);
    let mut smt =
        SubgridMarchingTetrahedra::<R>::new(SUBGRID_SAMPLES).expect("a positive sampling");
    let mut out = MeshBuffer::<R>::new();
    smt.extract(&field, &shape, origin, h, &mut out)
        .expect("extraction");

    c.bench_function(label, |b| {
        b.iter(|| {
            out.reset();
            smt.extract(&field, &shape, origin, h, &mut out)
                .expect("extraction");
            black_box(out.triangle_count())
        });
    });
}

/// The 1D sampling resolution the subgrid benchmarks run at.
///
/// Held constant across every subgrid row so the field-to-field comparison is
/// about the field. It is the same value the golden fixture pins, so a timing
/// and a hash refer to the same configuration.
const SUBGRID_SAMPLES: u32 = 16;

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
        bench_smt(
            c,
            &format!("subgrid_marching_tetrahedra/sphere/f32/{n}"),
            Sphere::<f32>::canonical(),
            n,
        );
        bench_smt(
            c,
            &format!("subgrid_marching_tetrahedra/torus/f32/{n}"),
            Torus::<f32>::canonical(),
            n,
        );
        bench_smt(
            c,
            &format!("subgrid_marching_tetrahedra/box_exact/f32/{n}"),
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
        bench_mdc(
            c,
            &format!("manifold_dual_contouring/sphere/f32/{n}"),
            Sphere::<f32>::canonical(),
            n,
        );
        bench_mdc(
            c,
            &format!("manifold_dual_contouring/torus/f32/{n}"),
            Torus::<f32>::canonical(),
            n,
        );
        bench_mdc(
            c,
            &format!("manifold_dual_contouring/box_exact/f32/{n}"),
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
        // The third rung, and the one that carries the interior rule's own cost:
        // `noise_cavity` is the only reference field with a cell the rule can
        // actually do something in (M-208), so this pair is the price of meshing
        // a tunnel against the price of ignoring one.
        bench_mc33(
            c,
            &format!("decider/marching_cubes+decider/noise_cavity/f32/{n}"),
            noise_cavity::<f32>(),
            n,
        );
        bench_trilinear(
            c,
            &format!("decider/marching_cubes+trilinear/noise_cavity/f32/{n}"),
            noise_cavity::<f32>(),
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
