//! A-014c's tests.
//!
//! The one the ticket exists for is `thin_plate_comes_back_where_greedy_quads_
//! returns_nothing`: A-005 measured zero triangles on that field, and this is
//! the same field through this extractor.

use super::*;
use crate::fields::{ReferenceField, ThinPlate};
use crate::mesh::MeshBuffer;
use crate::shape::RuntimeShape3;

/// A grid of `n` samples per axis over a field's own domain.
fn grid<F: ReferenceField<Scalar = f64>>(field: &F, n: u32) -> (RuntimeShape3, [f64; 3], f64) {
    let (lo, hi) = field.domain();
    let shape = RuntimeShape3::new([n; 3]).expect("a cubic grid");
    let cell = (hi[0] - lo[0]) / f64::from(n - 1);
    (shape, lo, cell)
}

#[test]
fn samples_must_be_positive() {
    assert!(SubgridMarchingTetrahedra::<f64>::new(0).is_err());
    assert!(SubgridMarchingTetrahedra::<f64>::new(1).is_ok());
}

#[test]
fn a_grid_with_no_cells_is_rejected() {
    let mut mt = SubgridMarchingTetrahedra::<f64>::new(8).expect("valid");
    let shape = RuntimeShape3::new([1, 4, 4]).expect("shape");
    let mut out = MeshBuffer::<f64>::default();
    assert_eq!(
        mt.extract(
            &ThinPlate::<f64>::canonical(),
            &shape,
            [0.0; 3],
            0.1,
            &mut out
        ),
        Err(crate::Error::GridTooSmall { size: [1, 4, 4] })
    );
}

#[test]
fn thin_plate_comes_back_where_greedy_quads_returns_nothing() {
    // A-014c's acceptance criterion, and the reason the whole subgrid track
    // exists. A-005 measured `thin_plate` -- 0.4 cells thick -- producing
    // **zero** triangles under greedy quads, because no cell centre is inside
    // it. M-72 measured Marching Cubes aliasing it into a resolution-dependent
    // scatter rather than resolving it.
    //
    // The plate is a sheet, so the mesh has to be a sheet: two sides and a rim,
    // not an empty buffer and not a handful of slivers.
    let field = ThinPlate::<f64>::canonical();
    let (shape, origin, cell) = grid(&field, 17);
    let mut mt = SubgridMarchingTetrahedra::<f64>::new(16).expect("valid");
    let mut out = MeshBuffer::<f64>::default();
    mt.extract(&field, &shape, origin, cell, &mut out)
        .expect("thin_plate should extract");

    // Pinned, not just "more than nothing": 896 triangles at 17³, against
    // greedy quads' zero on the same field. See M-95.
    assert_eq!(out.triangle_count(), 896);
    assert_eq!(out.vertex_count(), 2248);

    // Every vertex is on the surface it came from. This is the check that
    // separates "produced a lot of triangles" from "produced the *right*
    // triangles": a wrong root, a mis-mapped edge orientation or a stale buffer
    // would put a vertex somewhere the field is not zero.
    for p in &out.positions {
        let v = field.sample(*p);
        assert!(
            v.abs() < 1e-9,
            "vertex at {p:?} has field value {v}, so it is not on the surface"
        );
    }
}

#[test]
fn the_extractor_resolves_a_feature_the_grid_cannot() {
    // The claim in its sharpest form: hold the grid fixed and raise only the 1D
    // sampling. A slab thinner than a cell appears -- with no more grid
    // resolution at all -- which is the property M-67 quantified from the other
    // side and no sign-based method can reproduce.
    struct Slab {
        half: f64,
    }
    impl crate::Sdf for Slab {
        type Scalar = f64;
        fn sample(&self, p: [f64; 3]) -> f64 {
            // A slab in z, deliberately off the sample lattice so no grid plane
            // lands on its surface (M-94's fixture trap).
            (p[2] - 0.0137).abs() - self.half
        }
    }

    let shape = RuntimeShape3::new([9; 3]).expect("shape");
    let origin = [-1.0; 3];
    let cell = 0.25;

    // A slab 1/20 of a cell thick.
    let field = Slab { half: cell / 40.0 };
    let mut coarse = SubgridMarchingTetrahedra::<f64>::new(2).expect("valid");
    let mut sparse = MeshBuffer::<f64>::default();
    coarse
        .extract(&field, &shape, origin, cell, &mut sparse)
        .expect("extract");

    let mut fine = SubgridMarchingTetrahedra::<f64>::new(256).expect("valid");
    let mut dense = MeshBuffer::<f64>::default();
    fine.extract(&field, &shape, origin, cell, &mut dense)
        .expect("extract");

    assert_eq!(
        sparse.triangle_count(),
        0,
        "2 samples per edge should step over a slab this thin"
    );
    assert!(
        dense.triangle_count() > 0,
        "256 samples per edge should resolve it on the same grid"
    );
}

#[test]
fn extraction_is_deterministic() {
    let field = ThinPlate::<f64>::canonical();
    let (shape, origin, cell) = grid(&field, 9);
    let mut mt = SubgridMarchingTetrahedra::<f64>::new(8).expect("valid");

    let mut a = MeshBuffer::<f64>::default();
    let mut b = MeshBuffer::<f64>::default();
    mt.extract(&field, &shape, origin, cell, &mut a)
        .expect("extract");
    mt.extract(&field, &shape, origin, cell, &mut b)
        .expect("extract");

    assert_eq!(a.positions, b.positions);
    assert_eq!(a.normals, b.normals);
    assert_eq!(a.indices, b.indices);
    assert!(!a.is_empty());
}

#[test]
fn reusing_the_extractor_does_not_leak_the_previous_field() {
    // The buffers are held across calls, so a stale `along` or `patch` would
    // show up as geometry from the wrong field rather than as a crash.
    let plate = ThinPlate::<f64>::canonical();
    let (shape, origin, cell) = grid(&plate, 9);
    let mut mt = SubgridMarchingTetrahedra::<f64>::new(8).expect("valid");

    let mut first = MeshBuffer::<f64>::default();
    mt.extract(&plate, &shape, origin, cell, &mut first)
        .expect("extract");
    assert!(!first.is_empty());

    // A field with no surface in this domain at all.
    struct Empty;
    impl crate::Sdf for Empty {
        type Scalar = f64;
        fn sample(&self, _p: [f64; 3]) -> f64 {
            1.0
        }
    }
    let mut second = MeshBuffer::<f64>::default();
    mt.extract(&Empty, &shape, origin, cell, &mut second)
        .expect("extract");
    assert!(
        second.is_empty(),
        "a field with no surface produced {} triangles",
        second.triangle_count()
    );
}

#[test]
fn every_index_is_in_range_and_no_triangle_is_degenerate_by_index() {
    let field = ThinPlate::<f64>::canonical();
    let (shape, origin, cell) = grid(&field, 13);
    let mut mt = SubgridMarchingTetrahedra::<f64>::new(12).expect("valid");
    let mut out = MeshBuffer::<f64>::default();
    mt.extract(&field, &shape, origin, cell, &mut out)
        .expect("extract");

    let n = out.vertex_count() as u32;
    for tri in out.indices.chunks_exact(3) {
        assert!(tri.iter().all(|i| *i < n), "index past {n}: {tri:?}");
        assert!(
            tri[0] != tri[1] && tri[1] != tri[2] && tri[0] != tri[2],
            "degenerate by index: {tri:?}"
        );
    }
}

#[test]
fn raising_the_sampling_changes_the_topology_it_finds_but_not_by_much_else() {
    // Once the 1D sampling brackets a feature, refining further does not change
    // *what* is found -- same vertex count, same triangle count, same indices.
    //
    // It does move the vertices, and the first version of this test asserted
    // bit-equality and failed. Bisection converges to *an* ulp of the root, and
    // which one depends on the bracket it started from, so a different
    // `samples` gives a different last bit. Worth knowing precisely: the
    // determinism guarantee is "same arguments, same output", not "same field,
    // same output", and golden hashes over this extractor must therefore pin
    // `samples` alongside the grid. See M-95.
    let field = ThinPlate::<f64>::canonical();
    let (shape, origin, cell) = grid(&field, 33);

    let mut coarse = SubgridMarchingTetrahedra::<f64>::new(16).expect("valid");
    let mut a = MeshBuffer::<f64>::default();
    coarse
        .extract(&field, &shape, origin, cell, &mut a)
        .expect("extract");

    let mut fine = SubgridMarchingTetrahedra::<f64>::new(32).expect("valid");
    let mut b = MeshBuffer::<f64>::default();
    fine.extract(&field, &shape, origin, cell, &mut b)
        .expect("extract");

    assert_eq!(a.triangle_count(), 4328);
    assert_eq!(a.vertex_count(), b.vertex_count());
    assert_eq!(a.indices, b.indices, "the topology found should not change");

    let worst = a
        .positions
        .iter()
        .zip(b.positions.iter())
        .map(|(p, q)| (0..3).map(|k| (p[k] - q[k]).abs()).fold(0.0f64, f64::max))
        .fold(0.0f64, f64::max);
    assert!(
        worst < 1e-12,
        "16 and 32 samples disagree by {worst}, which is more than a refinement gap"
    );
    assert!(
        worst > 0.0,
        "if they agreed exactly this test is asserting nothing"
    );
}

#[test]
fn the_thin_plate_census_reproduces_with_the_real_root_finder() {
    // M-89 was measured with a throwaway probe carrying its own ad-hoc root
    // finder, deliberately not committed -- a second root finder in the repo is
    // the two-path failure the crate forbids. This is the same census through
    // `all_roots`, which is now the only one there is.
    //
    // It counts rather than meshes, so it says what the *field* presents to
    // §3.2 rather than what §3.2 does about it. The 620 multi-root edges are the
    // subgrid signal itself: every one is an edge a sign test reads as a single
    // bit.
    use crate::marching_tetrahedra::table::TETS;
    use crate::subgrid::curves::CurveKind;
    use crate::subgrid::surface::cycles;

    let field = ThinPlate::<f64>::canonical();
    let (shape, origin, cell) = grid(&field, 33);
    let size = shape.size();

    let mut tets_with_crossings = 0u32;
    let mut multi_root_edges = 0u32;
    let mut kinds = [0u32; 3];
    let mut along: [Vec<f64>; 6] = Default::default();

    for z in 0..size[2] - 1 {
        for y in 0..size[1] - 1 {
            for x in 0..size[0] - 1 {
                for tet in &TETS {
                    let mut corners = [[0.0f64; 3]; 4];
                    for (c, slot) in corners.iter_mut().enumerate() {
                        let offset = crate::cube::corner_offset(tet[c]);
                        for axis in 0..3 {
                            let index = f64::from([x, y, z][axis]) + f64::from(offset[axis]);
                            slot[axis] = origin[axis] + cell * index;
                        }
                    }

                    let mut count = [0u32; 6];
                    for (e, slot) in along.iter_mut().enumerate() {
                        slot.clear();
                        let [lo, hi] = crate::marching_tetrahedra::table::TET_EDGES[e];
                        crate::subgrid::roots::all_roots(
                            corners[lo as usize],
                            corners[hi as usize],
                            &field,
                            64,
                            slot,
                        );
                        count[e] = slot.len() as u32;
                        if slot.len() > 1 {
                            multi_root_edges += 1;
                        }
                    }
                    if count.iter().all(|n| *n == 0) {
                        continue;
                    }
                    tets_with_crossings += 1;
                    for c in cycles(&crate::subgrid::coordinates::EdgeCoordinates::new(count)) {
                        match c.kind {
                            CurveKind::Open => kinds[0] += 1,
                            CurveKind::Normal => kinds[1] += 1,
                            CurveKind::NonNormal => kinds[2] += 1,
                        }
                    }
                }
            }
        }
    }

    assert_eq!(tets_with_crossings, 3072, "tets carrying crossings");
    assert_eq!(multi_root_edges, 620, "edges with more than one root");
    assert_eq!(kinds, [0, 3060, 192], "curves open / normal / non-normal");
}

#[test]
fn the_welded_output_is_a_closed_consistently_oriented_manifold() {
    // A-014e's acceptance. §3.2 fixes each polygon's winding from its own
    // boundary curve, which knows nothing about which side the field calls
    // inside, so the extractor imposes it per triangle against the gradient --
    // per triangle and not per patch, because `thin_plate`'s two faces are 0.4
    // cells apart and land in the same tetrahedron.
    //
    // **The weld is not a convenience here, it is a precondition.** The
    // extractor emits each tetrahedron's vertices independently, so before
    // welding the surface is 896 disconnected triangles: 2,240 boundary edges,
    // and an orientation check that can only see the edges interior to a single
    // tetrahedron. Welding is what turns it into a surface at all. See M-96.
    use crate::validate::{ValidateConfig, validate_indexed};

    for (name, n, samples, expected_chi) in
        [("thin_plate", 17u32, 16u32, 2i64), ("sphere", 17, 16, 2)]
    {
        let (mut out, cell) = match name {
            "thin_plate" => {
                let field = ThinPlate::<f64>::canonical();
                let (shape, origin, cell) = grid(&field, n);
                let mut mt = SubgridMarchingTetrahedra::<f64>::new(samples).expect("valid");
                let mut out = MeshBuffer::<f64>::default();
                mt.extract(&field, &shape, origin, cell, &mut out)
                    .expect("extract");
                (out, cell)
            }
            _ => {
                let field = crate::fields::Sphere::<f64>::canonical();
                let (shape, origin, cell) = grid(&field, n);
                let mut mt = SubgridMarchingTetrahedra::<f64>::new(samples).expect("valid");
                let mut out = MeshBuffer::<f64>::default();
                mt.extract(&field, &shape, origin, cell, &mut out)
                    .expect("extract");
                (out, cell)
            }
        };

        let mut welder = crate::weld::Welder::<f64>::new();
        welder
            .weld(&mut out, cell * 1e-6)
            .unwrap_or_else(|e| panic!("{name}: weld failed: {e}"));

        let cfg = ValidateConfig::from_cell_size(cell).expect("a valid spacing");
        let report = validate_indexed(&out.positions, &out.indices, &cfg);

        assert_eq!(report.inconsistently_oriented_edges, 0, "{name}: {report}");
        assert_eq!(report.non_manifold_edges, 0, "{name}: {report}");
        assert_eq!(report.non_manifold_vertices, 0, "{name}: {report}");
        assert_eq!(report.boundary_edges, 0, "{name}: {report}");
        assert!(report.is_closed(), "{name}: {report}");
        assert_eq!(
            report.euler_characteristic, expected_chi,
            "{name}: {report}"
        );
    }
}
