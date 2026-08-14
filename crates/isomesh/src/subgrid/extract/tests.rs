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

    // Pinned, not just "more than nothing": 840 triangles at 17³, against
    // greedy quads' zero on the same field. See M-95.
    //
    // **896 before A-014h, 840 after.** The 56 that went are the zero-area
    // slivers between two roots that both sit on one grid point: naming them by
    // the point makes them one vertex, which makes the triangle between them a
    // repeated-index degenerate, and the extractor declines to emit those
    // (M-185). None of them carried area at any point in this crate's history.
    assert_eq!(out.triangle_count(), 840);
    // **2248 before A-014g's shared vertex table, 450 after it, 422 after
    // A-014h** — and 422 is exactly what a positional weld reaches, so the
    // extractor now emits each vertex once and there is nothing left to merge.
    assert_eq!(out.vertex_count(), 422);

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

    assert_eq!(a.triangle_count(), 3964);
    assert_eq!(a.vertex_count(), b.vertex_count());
    assert_eq!(a.indices, b.indices, "the topology found should not change");

    // **A-014h turned this from "close" into "identical" (M-186).** The
    // paragraph above records that the first version of this test asserted
    // bit-equality and failed, because bisection converges to *an* ulp of the
    // root and which one depended on the bracket it started from. That is no
    // longer true here: every vertex whose last bit moved with `samples` was a
    // root sitting on a grid point, and those are now emitted at the grid
    // point's own canonical position, which no amount of refinement can move.
    //
    // Asserted at the bit rather than by a tolerance, because a tolerance would
    // no longer be measuring anything.
    let differing = a
        .positions
        .iter()
        .zip(b.positions.iter())
        .filter(|(p, q)| (0..3).any(|k| p[k].to_bits() != q[k].to_bits()))
        .count();
    assert_eq!(
        differing,
        0,
        "{differing} of {} vertices moved when the 1D sampling doubled -- \
         A-014h's canonical placement has stopped covering them, and the \
         determinism guarantee this test pins is weaker again",
        a.vertex_count()
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
    // One generic body, dispatched by type rather than by a `match` on the
    // field's name -- the ladder this replaces is exactly the shape the
    // testing law forbids ("no `if field == gyroid` anywhere in test code"),
    // and its two arms had already been copied identically once.
    fn check<F: ReferenceField<Scalar = f64>>(field: &F, n: u32, samples: u32, expected_chi: i64) {
        use crate::validate::{ValidateConfig, validate_indexed};
        let name = F::NAME;

        let (shape, origin, cell) = grid(field, n);
        let mut mt = SubgridMarchingTetrahedra::<f64>::new(samples).expect("valid");
        let mut out = MeshBuffer::<f64>::default();
        mt.extract(field, &shape, origin, cell, &mut out)
            .expect("extract");

        let mut welder = crate::weld::Welder::<f64>::new();
        welder
            .weld(&mut out, crate::weld::epsilon_for(cell))
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

    check(&ThinPlate::<f64>::canonical(), 17, 16, 2);
    check(&crate::fields::Sphere::<f64>::canonical(), 17, 16, 2);
}

#[test]
fn the_validity_suite_over_every_reference_field() {
    // The definition of done's T-001 requirement, with the gate chosen by the
    // field rather than by the test -- `gyroid` is triply periodic and
    // `fbm_terrain` is a heightfield, so neither is closed in a finite box and
    // demanding `is_closed()` of them would be asserting something false.
    //
    // Welded first, without exception: M-96 measured that the raw output is a
    // per-tetrahedron triangle soup, so every one of these counters is
    // meaningless before the weld.
    //
    // Counts are **pinned rather than asserted to zero**, following the Phase 1
    // amendment: a known defect with a number and a ticket that owns it
    // satisfies this gate; an unexplained one does not.
    //
    // **The `inconsistently_oriented` column is the extractor's output, before
    // A-014f's pass.** That pass drives every one of these to zero on any row
    // whose edges are all manifold (M-187), so this column measures what §3.2's
    // winding plus A-014e's per-triangle vote produce on their own -- which is
    // the thing this suite is for -- and `orient`'s own test measures the
    // pipeline. Keeping them apart is what makes it visible when one of them
    // regresses.
    use crate::fields::ReferenceField;
    use crate::validate::{ValidateConfig, validate_indexed};

    // (field, samples) -> (non-manifold edges, non-manifold vertices,
    // inconsistently oriented). **Three resolutions, because Phase 1's gate says
    // three and this ran at 17 alone until T-010** -- and 17³ turned out to be
    // the resolution that decided which defects this project believed it had.
    const RESOLUTIONS: [u32; 3] = [17, 25, 33];
    #[rustfmt::skip]
    let expected: [(&str, u32, u64, u64, u64); 21] = [
        ("sphere",         17, 0, 0,   0), ("torus",       17, 0, 0,   0),
        ("box_exact",      17, 0, 0,   0), ("csg_difference", 17, 3, 6, 6),
        ("thin_plate",     17, 0, 0,   0), ("gyroid",      17, 0, 0, 138),
        ("fbm_terrain",    17, 4, 6,  19),

        ("sphere",         25, 0, 0,   0), ("torus",       25, 4, 6,   6),
        ("box_exact",      25, 0, 0,   0), ("csg_difference", 25, 0, 0, 0),
        ("thin_plate",     25, 0, 0,   8), ("gyroid",      25, 0, 0, 150),
        ("fbm_terrain",    25, 2, 3,  29),

        ("sphere",         33, 0, 0,   0), ("torus",       33, 0, 0,   0),
        ("box_exact",      33, 0, 0,   0), ("csg_difference", 33, 0, 0, 36),
        ("thin_plate",     33, 0, 0,   6), ("gyroid",      33, 0, 0, 330),
        ("fbm_terrain",    33, 8, 12, 53),
    ];

    let mut checked = 0;
    for n in RESOLUTIONS {
        crate::for_each_reference_field!(f64, |name, field| {
            let (lo, hi) = field.domain();
            let shape = RuntimeShape3::new([n; 3]).expect("a cubic grid");
            let cell = (hi[0] - lo[0]) / f64::from(n - 1);

            let mut mt = SubgridMarchingTetrahedra::<f64>::new(16).expect("valid");
            let mut out = MeshBuffer::<f64>::default();
            mt.extract(&field, &shape, lo, cell, &mut out)
                .unwrap_or_else(|e| panic!("{name}: {e}"));

            let mut welder = crate::weld::Welder::<f64>::new();
            welder
                .weld(&mut out, crate::weld::epsilon_for(cell))
                .unwrap_or_else(|e| panic!("{name}: weld failed: {e}"));

            let cfg = ValidateConfig::from_cell_size(cell).expect("a valid spacing");
            let report = validate_indexed(&out.positions, &out.indices, &cfg);

            let want = expected
                .iter()
                .find(|(f, s, ..)| *f == name && *s == n)
                .unwrap_or_else(|| panic!("{name} {n} is not in the pinned table"));
            std::println!(
                "{name:<15} {n}³  nm-edges {:>2}  nm-verts {:>2}  flipped {:>3}",
                report.non_manifold_edges,
                report.non_manifold_vertices,
                report.inconsistently_oriented_edges
            );
            assert_eq!(
                (
                    report.non_manifold_edges,
                    report.non_manifold_vertices,
                    report.inconsistently_oriented_edges
                ),
                (want.2, want.3, want.4),
                "{name} {n}³\n{report}"
            );

            // Every field must at least be a surface with no dangling geometry, and
            // the clean ones must meet their own gate.
            assert_eq!(report.out_of_range_indices, 0, "{name}\n{report}");
            // Only the fields with a clean bill on all three counters are held to
            // their own gate; the others are pinned above and owned by A-014d.
            if want.2 == 0 && want.3 == 0 && want.4 == 0 {
                if field.closed_in_domain() {
                    assert!(report.is_closed(), "{name}\n{report}");
                } else {
                    assert!(report.is_manifold(), "{name}\n{report}");
                }
                if let Some(chi) = field.expected_euler() {
                    assert_eq!(report.euler_characteristic, chi, "{name}\n{report}");
                }
            }
            checked += 1;
        });
    }
    assert_eq!(
        checked, 21,
        "the sweep did not reach every field at every size"
    );
}

// ---------------------------------------------------------------------------
// A-014d's unblocking measurement.
//
// The ticket's `BLOCKED:` line asks for one thing and not another: *"instrument
// the extractor to report which polygons actually coincide across a shared face
// on `csg_difference`, and see whether that condition is expressible without
// neighbour information."* Two implementation attempts of §3.2.3's inset are on
// record and both made the measured result worse (M-101), so what follows
// measures the premise instead of trying a third time.
// ---------------------------------------------------------------------------

/// One tetrahedron's slice of the triangle soup.
///
/// The extractor emits vertices per tetrahedron and shares none, so a triangle's
/// provenance is recoverable from nothing but the order it was written in — as
/// long as something records the boundaries while the extraction runs.
struct TetRun {
    cell: [u32; 3],
    tet: usize,
    first: usize,
    count: usize,
}

/// The same extraction [`SubgridMarchingTetrahedra::extract`] performs, with the
/// tetrahedron boundaries recorded.
///
/// Driving `cell_tet` directly rather than re-deriving provenance afterwards is
/// what makes this measure the extractor instead of a lookalike, and the caller
/// asserts the two agree bit-for-bit.
fn soup_with_provenance<F: Sdf<Scalar = f64>>(
    field: &F,
    origin: [f64; 3],
    cell_size: f64,
    n: u32,
    samples: u32,
) -> (MeshBuffer<f64>, Vec<TetRun>) {
    let mut mt = SubgridMarchingTetrahedra::<f64>::new(samples).expect("valid");
    let mut out = MeshBuffer::<f64>::default();
    let mut runs = Vec::new();
    let mut vertices: u64 = 0;
    for z in 0..n - 1 {
        for y in 0..n - 1 {
            for x in 0..n - 1 {
                for t in 0..TETS.len() {
                    let first = out.indices.len() / 3;
                    mt.cell_tet(
                        field,
                        origin,
                        cell_size,
                        [x, y, z],
                        t,
                        &mut vertices,
                        &mut out,
                    )
                    .unwrap_or_else(|e| panic!("cell [{x}, {y}, {z}] tet {t}: {e}"));
                    let count = out.indices.len() / 3 - first;
                    if count > 0 {
                        runs.push(TetRun {
                            cell: [x, y, z],
                            tet: t,
                            first,
                            count,
                        });
                    }
                }
            }
        }
    }
    (out, runs)
}

/// A tetrahedron's four corners, by the expression `cell_tet` itself uses.
///
/// Written the same way for the same reason M-32 gives: two tetrahedra sharing a
/// corner agree bit-for-bit only if the same expression computed it, so this is
/// a copy of that loop rather than an equivalent one.
fn tet_corners(origin: [f64; 3], cell_size: f64, cell: [u32; 3], t: usize) -> [[u64; 3]; 4] {
    let mut corners = [[0u64; 3]; 4];
    for (c, slot) in corners.iter_mut().enumerate() {
        let offset = corner_offset(TETS[t][c]);
        for axis in 0..3 {
            let index = f64::from(cell[axis]) + f64::from(offset[axis]);
            slot[axis] = (origin[axis] + cell_size * index).to_bits();
        }
    }
    corners
}

/// A position as bits, so coincidence is exact rather than within a tolerance.
///
/// The extractor's own guarantee is bit-identical positions from the two sides
/// of a shared face, so an epsilon here would be measuring the epsilon.
fn point_key(p: [f64; 3]) -> [u64; 3] {
    [p[0].to_bits(), p[1].to_bits(), p[2].to_bits()]
}

/// A triangle's three points, sorted — equal for two copies of one polygon
/// however each was wound.
fn triangle_key(mesh: &MeshBuffer<f64>, tri: usize) -> [[u64; 3]; 3] {
    let mut key = triangle_points(mesh, tri);
    key.sort_unstable();
    key
}

/// A triangle's three points, in the order it was wound.
fn triangle_points(mesh: &MeshBuffer<f64>, tri: usize) -> [[u64; 3]; 3] {
    let at = |k: usize| point_key(mesh.positions[mesh.indices[tri * 3 + k] as usize]);
    [at(0), at(1), at(2)]
}

/// The winding, as the rotation starting at the smallest point.
///
/// Two triangles over the same three points wind the same way iff these agree;
/// §3.2.3's immersed pairs are expected to disagree, which is what makes them
/// two copies of a surface rather than one surface counted twice.
fn winding(mesh: &MeshBuffer<f64>, tri: usize) -> [[u64; 3]; 3] {
    let p = triangle_points(mesh, tri);
    let start = (0..3).min_by_key(|k| p[*k]).unwrap_or(0);
    [p[start], p[(start + 1) % 3], p[(start + 2) % 3]]
}

/// A triangle's three edges, each as its two points sorted.
fn triangle_edges(mesh: &MeshBuffer<f64>, tri: usize) -> [[[u64; 3]; 2]; 3] {
    let p = triangle_points(mesh, tri);
    let mut edges = [[p[0], p[1]], [p[1], p[2]], [p[2], p[0]]];
    for e in &mut edges {
        e.sort_unstable();
    }
    edges
}

/// **A-014d's unblocking measurement**, and the answer it produced.
///
/// The question is whether a tetrahedron can tell, on its own, that a polygon it
/// is about to emit will also be emitted by its neighbour — because §3.2.3's
/// inset has to be applied to *both* copies or neither, and this extractor is
/// strictly per-tetrahedron.
///
/// Three things are counted, and only the third decides the ticket:
///
/// 1. **How many polygons coincide, and between how many tetrahedra.** A pair is
///    the immersion §3.2.3 describes. Anything else is a different defect.
/// 2. **Whether the emission is symmetric.** If a tetrahedron ever emits a
///    face polygon its neighbour does not, no per-tetrahedron rule can predict
///    the duplication, and the inset needs neighbour information for that reason
///    alone.
/// 3. **How many *other* tetrahedra use the edges of a coincident polygon.**
///    This is M-101's failure, stated as a number. §3.2.3 moves the midpoints of
///    polygon edges lying in an edge of the shared face; every other tetrahedron
///    meeting at that edge keeps its vertices where they were, so each one is a
///    seam the inset tears open unless it moves too.
#[test]
fn which_polygons_coincide_across_a_shared_face() {
    use alloc::collections::{BTreeMap, BTreeSet};

    use crate::fields::ReferenceField;

    // field -> (coincident, pairs inside one cell, polygons with foreign edge
    // users, foreign edge users, of those in another cell)
    let expected: [(&str, u64, u64, u64, u64, u64); 7] = [
        ("sphere", 0, 0, 0, 0, 0),
        ("torus", 0, 0, 0, 0, 0),
        ("box_exact", 0, 0, 0, 0, 0),
        ("csg_difference", 3, 3, 0, 0, 0),
        ("thin_plate", 0, 0, 0, 0, 0),
        ("gyroid", 0, 0, 0, 0, 0),
        ("fbm_terrain", 4, 4, 0, 0, 0),
    ];

    let mut rows = 0;
    crate::for_each_reference_field!(f64, |name, field| {
        let (lo, hi) = field.domain();
        let n = 17u32;
        let cell = (hi[0] - lo[0]) / f64::from(n - 1);
        let (mesh, runs) = soup_with_provenance(&field, lo, cell, n, 16);

        // The instrument must be measuring the extractor and not something that
        // resembles it. Same loop, same order, same output -- asserted, because
        // a measurement of a lookalike would be worse than no measurement.
        let shape = RuntimeShape3::new([n; 3]).expect("a cubic grid");
        let mut reference = MeshBuffer::<f64>::default();
        SubgridMarchingTetrahedra::<f64>::new(16)
            .expect("valid")
            .extract(&field, &shape, lo, cell, &mut reference)
            .unwrap_or_else(|e| panic!("{name}: {e}"));
        assert_eq!(
            mesh.positions, reference.positions,
            "{name}: instrument drift"
        );
        assert_eq!(mesh.indices, reference.indices, "{name}: instrument drift");

        let triangles = mesh.indices.len() / 3;
        let mut owner = alloc::vec![usize::MAX; triangles];
        for (r, run) in runs.iter().enumerate() {
            for slot in owner.iter_mut().skip(run.first).take(run.count) {
                *slot = r;
            }
        }

        // Every triangle occupying a given three points, and every triangle
        // touching a given two.
        let mut by_shape: BTreeMap<[[u64; 3]; 3], Vec<usize>> = BTreeMap::new();
        let mut by_edge: BTreeMap<[[u64; 3]; 2], Vec<usize>> = BTreeMap::new();
        for tri in 0..triangles {
            by_shape
                .entry(triangle_key(&mesh, tri))
                .or_default()
                .push(tri);
            for edge in triangle_edges(&mesh, tri) {
                by_edge.entry(edge).or_default().push(tri);
            }
        }

        let mut coincident = 0u64;
        let mut over_two_tets = 0u64;
        let mut not_face_adjacent = 0u64;
        let mut same_winding = 0u64;
        let mut with_foreign_edge_users = 0u64;
        let mut worst_foreign = 0usize;
        let mut foreign_total = 0u64;
        // The architectural question, and the only one whose answer changes what
        // A-014d costs: an inset confined to one cell is a change to `cell_tet`'s
        // caller, and one that crosses cells is a change to the extractor.
        let mut pair_in_one_cell = 0u64;
        let mut foreign_in_one_cell = 0u64;
        let mut foreign_across_cells = 0u64;

        for group in by_shape.values() {
            let tets: BTreeSet<usize> = group.iter().map(|tri| owner[*tri]).collect();
            if tets.len() < 2 {
                continue;
            }
            coincident += 1;
            if tets.len() > 2 {
                over_two_tets += 1;
            }

            // Do the two tetrahedra actually share a face? Three common corners
            // is a face; anything less means the coincidence is not §3.2.3's.
            if tets.len() == 2 {
                let mut it = tets.iter();
                let (Some(a), Some(b)) = (it.next(), it.next()) else {
                    continue;
                };
                let ca = tet_corners(lo, cell, runs[*a].cell, runs[*a].tet);
                let cb = tet_corners(lo, cell, runs[*b].cell, runs[*b].tet);
                let shared = ca.iter().filter(|c| cb.contains(c)).count();
                if shared != 3 {
                    not_face_adjacent += 1;
                }
            }

            // Two copies of one polygon must be wound oppositely to be two sides
            // of a sheet rather than the same side twice.
            let windings: BTreeSet<[[u64; 3]; 3]> =
                group.iter().map(|tri| winding(&mesh, *tri)).collect();
            if windings.len() < group.len() {
                same_winding += 1;
            }

            // The decider: who else is standing on this polygon's edges.
            let cells: BTreeSet<[u32; 3]> = tets.iter().map(|t| runs[*t].cell).collect();
            if cells.len() == 1 {
                pair_in_one_cell += 1;
            }

            let mut foreign = 0usize;
            for edge in triangle_edges(&mesh, group[0]) {
                if let Some(users) = by_edge.get(&edge) {
                    for tri in users.iter().filter(|tri| !tets.contains(&owner[**tri])) {
                        foreign += 1;
                        if cells.contains(&runs[owner[*tri]].cell) {
                            foreign_in_one_cell += 1;
                        } else {
                            foreign_across_cells += 1;
                        }
                    }
                }
            }
            if foreign > 0 {
                with_foreign_edge_users += 1;
            }
            foreign_total += foreign as u64;
            worst_foreign = worst_foreign.max(foreign);
        }

        std::println!(
            "{name:<15} {triangles:>7} tris {:>5} tets | coincident {coincident:>4} \
             (>2 tets {over_two_tets}, not face-adjacent {not_face_adjacent}, \
             same winding {same_winding}) | with foreign edge users {with_foreign_edge_users} \
             (total {foreign_total}, worst {worst_foreign}) | pair in one cell \
             {pair_in_one_cell} | foreign same-cell {foreign_in_one_cell} \
             cross-cell {foreign_across_cells}",
            runs.len(),
        );

        // Two copies of one polygon must be wound oppositely to be two sides of
        // a sheet. Every one of them is wound the *same* way, because A-014e
        // imposes winding from the gradient at the centroid and two coincident
        // triangles share a centroid -- so the duplication survives the weld as
        // a doubled face rather than as a two-sided sheet.
        assert_eq!(same_winding, coincident, "{name}: a pair wound oppositely");
        let want = expected
            .iter()
            .find(|(f, ..)| *f == name)
            .unwrap_or_else(|| panic!("{name} is not in the pinned table"));
        assert_eq!(
            (
                coincident,
                pair_in_one_cell,
                with_foreign_edge_users,
                foreign_total,
                foreign_across_cells
            ),
            (want.1, want.2, want.3, want.4, want.5),
            "{name}"
        );
        rows += 1;
    });
    assert_eq!(rows, 7, "the sweep did not reach every field");

    // **The answer to A-014d's blocking question, as the table above now records
    // it.** An earlier version of this comment described 33 coincident polygons
    // on `csg_difference`, 27 with foreign edge users, 312 such triangles and
    // 150 of them in another cell -- and concluded the inset was an
    // architectural change no per-cell rule could carry. A-014h falsified every
    // one of those numbers (M-186) by giving each crossing a complete identity:
    // the row is `(3, 3, 0, 0, 0)`, which the assertion above enforces. Most of
    // what M-162 measured was one grid point wearing many names, not §3.2.3's
    // immersion. The stale prose outlived the numbers it described by one
    // ticket, which is exactly the failure X3 names.
}

/// **Which of Figure 15's three polygon types A-014d actually has to inset.**
///
/// §3.2.3 names three — *"There are three polygon types: a quad, a hexagon, and
/// a pentagon made at a corner"* — and gives their triangulation patterns only
/// as a picture, *"Figure 15, right"*. That figure has now been read at
/// 1200–2400 DPI: the quad's and the hexagon's patterns transcribe
/// unambiguously, and the pentagon's second inserted point **cannot** be told
/// apart from the drawing, so rule 5 still binds on that one case.
///
/// It only blocks the ticket if a pentagon is reachable, and this is the
/// measurement that says. It also says which patterns need writing at all —
/// implementing a type no field produces would be unverifiable code.
///
/// The polygon is a [`Region`](crate::subgrid::surface::Region) of one face, and
/// its boundary already carries the distinction Figure 15 draws in colour:
/// [`Arc::Chord`] is a segment of `γ` — the thick blue edges, along which the
/// two copies stay glued, which is why the pair is a *tube* rather than a sphere
/// — and [`Arc::Edge`] is a piece of an edge of the face, which is exactly what
/// §3.2.3 inserts a midpoint into and pushes inward.
#[test]
fn which_polygon_types_coincide_across_a_shared_face() {
    use alloc::collections::BTreeMap;

    use crate::fields::ReferenceField;
    use crate::subgrid::coordinates::TET_FACE_COUNT;
    use crate::subgrid::curves::CurveKind;
    use crate::subgrid::surface::{Arc, Node, NonNormalKind, cycles, face_regions};

    /// A polygon's `(loop kind, chords, edges)` — the figure's own vocabulary.
    type Shape = (&'static str, usize, usize);

    let mut rows = 0;
    let mut totals: BTreeMap<Shape, u64> = BTreeMap::new();
    crate::for_each_reference_field!(f64, |name, field| {
        let (lo, hi) = field.domain();
        let n = 17u32;
        let cell = (hi[0] - lo[0]) / f64::from(n - 1);
        let samples = 16u32;

        let mut along: [Vec<f64>; TET_EDGE_COUNT] = core::array::from_fn(|_| Vec::new());
        // The polygon's node ring, sorted, against every (loop kind, chords,
        // edges) seen on it. Two tetrahedra sharing a face build the same ring
        // from the same crossings, so a key with two entries *is* a coincident
        // pair. The loop kind rides along because Figure 15's third type is a
        // pentagon **made at a corner** — a 5-arc region from a contractible
        // loop would not be that polygon, and its pattern would not apply.
        let mut by_ring: BTreeMap<Vec<[u64; 3]>, Vec<Shape>> = BTreeMap::new();

        for z in 0..n - 1 {
            for y in 0..n - 1 {
                for x in 0..n - 1 {
                    for tet in &TETS {
                        let mut corners = [[0.0f64; 3]; 4];
                        for (c, slot) in corners.iter_mut().enumerate() {
                            let offset = corner_offset(tet[c]);
                            for axis in 0..3 {
                                let index = f64::from([x, y, z][axis]) + f64::from(offset[axis]);
                                slot[axis] = lo[axis] + cell * index;
                            }
                        }

                        let mut total = 0usize;
                        for (e, slot) in along.iter_mut().enumerate() {
                            slot.clear();
                            let [a, b] = TET_EDGES[e];
                            all_roots(
                                corners[a as usize],
                                corners[b as usize],
                                &field,
                                samples,
                                slot,
                            );
                            total += slot.len();
                        }
                        if total == 0 {
                            continue;
                        }

                        let mut borrowed: [&[f64]; TET_EDGE_COUNT] = [&[]; TET_EDGE_COUNT];
                        for (slot, v) in borrowed.iter_mut().zip(along.iter()) {
                            *slot = v.as_slice();
                        }
                        let crossings = TetCrossings {
                            corners,
                            along: borrowed,
                        };
                        if crossings.check().is_err() {
                            continue;
                        }
                        let coords = crossings.coordinates();

                        for cycle in cycles(&coords)
                            .iter()
                            .filter(|c| c.kind == CurveKind::NonNormal)
                        {
                            // Only the two kinds whose disk is built in the tet
                            // boundary can coincide with the neighbour's copy;
                            // a diagonal loop fans over an interior centroid.
                            let (kind, inside_corner) = match cycle.non_normal_kind() {
                                Some(NonNormalKind::Contractible) => ("contractible", None),
                                Some(NonNormalKind::Corner) => {
                                    ("corner", cycle.distinguished_corner())
                                }
                                _ => continue,
                            };
                            for face in 0..TET_FACE_COUNT as u8 {
                                let Some(regions) = face_regions(face, &coords, cycle) else {
                                    continue;
                                };
                                for region in regions.iter().filter(|r| r.is_inside()) {
                                    let mut ring = Vec::new();
                                    let mut whole = true;
                                    for node in &region.node {
                                        match node {
                                            Node::Crossing(p) => match crossings.position(*p) {
                                                Some(q) => ring.push(point_key(q)),
                                                None => whole = false,
                                            },
                                            Node::Corner(c) if Some(*c) == inside_corner => {}
                                            Node::Corner(c) => {
                                                ring.push(point_key(corners[*c as usize]));
                                            }
                                        }
                                    }
                                    if !whole {
                                        continue;
                                    }
                                    let chords = region
                                        .arc
                                        .iter()
                                        .filter(|a| matches!(a, Arc::Chord))
                                        .count();
                                    let edges = region.arc.len() - chords;
                                    ring.sort_unstable();
                                    ring.dedup();
                                    by_ring.entry(ring).or_default().push((kind, chords, edges));
                                }
                            }
                        }
                    }
                }
            }
        }

        let mut shape: BTreeMap<(&'static str, usize, usize), u64> = BTreeMap::new();
        for seen in by_ring.values().filter(|v| v.len() >= 2) {
            let entry = shape.entry(seen[0]).or_default();
            *entry += 1;
            *totals.entry(seen[0]).or_default() += 1;
        }
        std::println!("{name:<15} coincident by (kind, chords, edges): {shape:?}");
        rows += 1;
    });
    assert_eq!(rows, 7, "the sweep did not reach every field");
    std::println!("all fields: {totals:?}");
}

/// **Where A-014d's three defects actually come from**, traced to the
/// tetrahedra that emitted them.
///
/// The first measurement above answered the ticket's question and raised a
/// louder one: `box_exact` carries 30 coincident polygons and validates
/// **clean**, while `csg_difference` carries 33 and reports only 3 non-manifold
/// edges. If duplication were the defect those two numbers would track each
/// other, and they do not — so this asks which edges are actually bad and what
/// is standing on them.
///
/// The census is taken on the soup by **position**, because that is what welding
/// merges by, and it is checked against the welded mesh's own report rather than
/// trusted: a provenance census that disagreed with `validate_indexed` would be
/// describing a different mesh.
#[test]
fn the_defects_traced_back_to_the_tetrahedra_that_made_them() {
    use alloc::collections::{BTreeMap, BTreeSet};

    use crate::fields::ReferenceField;
    use crate::validate::{ValidateConfig, validate_indexed};

    // field -> (bad edges, of those duplication-only, three distinct polygons,
    // collapsed triangles)
    let expected: [(&str, u64, u64, u64, u64); 7] = [
        ("sphere", 0, 0, 0, 0),
        ("torus", 0, 0, 0, 0),
        ("box_exact", 0, 0, 0, 0),
        ("csg_difference", 3, 0, 3, 0),
        ("thin_plate", 0, 0, 0, 0),
        ("gyroid", 0, 0, 0, 0),
        ("fbm_terrain", 4, 2, 2, 0),
    ];

    let mut rows = 0;
    crate::for_each_reference_field!(f64, |name, field| {
        let (lo, hi) = field.domain();
        let n = 17u32;
        let cell = (hi[0] - lo[0]) / f64::from(n - 1);
        let (mesh, runs) = soup_with_provenance(&field, lo, cell, n, 16);
        let triangles = mesh.indices.len() / 3;
        let mut owner = alloc::vec![usize::MAX; triangles];
        for (r, run) in runs.iter().enumerate() {
            for slot in owner.iter_mut().skip(run.first).take(run.count) {
                *slot = r;
            }
        }

        // Triangles with a repeated point collapse to nothing under the weld, so
        // they are counted and set aside rather than allowed to invent edges.
        let mut collapsed = 0u64;
        let mut by_edge: BTreeMap<[[u64; 3]; 2], Vec<usize>> = BTreeMap::new();
        for tri in 0..triangles {
            let p = triangle_points(&mesh, tri);
            if p[0] == p[1] || p[1] == p[2] || p[2] == p[0] {
                collapsed += 1;
                continue;
            }
            for edge in triangle_edges(&mesh, tri) {
                by_edge.entry(edge).or_default().push(tri);
            }
        }

        // The same mesh through the real gate, so the census can be checked
        // rather than believed.
        let shape = RuntimeShape3::new([n; 3]).expect("a cubic grid");
        let mut welded = MeshBuffer::<f64>::default();
        SubgridMarchingTetrahedra::<f64>::new(16)
            .expect("valid")
            .extract(&field, &shape, lo, cell, &mut welded)
            .unwrap_or_else(|e| panic!("{name}: {e}"));
        crate::weld::Welder::<f64>::new()
            .weld(&mut welded, crate::weld::epsilon_for(cell))
            .unwrap_or_else(|e| panic!("{name}: weld failed: {e}"));
        let cfg = ValidateConfig::from_cell_size(cell).expect("a valid spacing");
        let report = validate_indexed(&welded.positions, &welded.indices, &cfg);

        // Of the bad edges, how many are bad *only* because one polygon was
        // emitted twice? Collapse each edge's faces to distinct shapes: if
        // fewer than three remain, duplication is the whole story and §3.2.3's
        // inset is the fix. If three or more distinct polygons genuinely meet,
        // it is not.
        let mut bad = 0u64;
        let mut duplication_only = 0u64;
        let mut genuinely_three = 0u64;
        let mut tets_on_worst = 0usize;
        for tris in by_edge.values().filter(|t| t.len() >= 3) {
            bad += 1;
            let distinct: BTreeSet<[[u64; 3]; 3]> =
                tris.iter().map(|t| triangle_key(&mesh, *t)).collect();
            if distinct.len() < 3 {
                duplication_only += 1;
            } else {
                genuinely_three += 1;
            }
            let tets: BTreeSet<usize> = tris.iter().map(|t| owner[*t]).collect();
            tets_on_worst = tets_on_worst.max(tets.len());
        }
        let soup_boundary = by_edge.values().filter(|t| t.len() == 1).count();

        std::println!(
            "{name:<15} welded nm-edges {:>3} nm-verts {:>3} flipped {:>3} boundary {:>4} \
             | soup: bad edges {bad:>3} (duplication-only {duplication_only}, \
             three-distinct {genuinely_three}, worst tets {tets_on_worst}) \
             collapsed tris {collapsed:>3} soup-boundary {soup_boundary}",
            report.non_manifold_edges,
            report.non_manifold_vertices,
            report.inconsistently_oriented_edges,
            report.boundary_edges,
        );
        let want = expected
            .iter()
            .find(|(f, ..)| *f == name)
            .unwrap_or_else(|| panic!("{name} is not in the pinned table"));
        assert_eq!(
            (bad, duplication_only, genuinely_three, collapsed),
            (want.1, want.2, want.3, want.4),
            "{name}"
        );
        rows += 1;
    });
    assert_eq!(rows, 7, "the sweep did not reach every field");

    // **`box_exact` is the result that re-aims the ticket.** It carries the most
    // coincident polygons of any field -- 30, with 348 triangles standing on
    // their edges -- and validates completely clean. Coincidence is therefore
    // not the defect, and 468 of its triangles are zero-area, which §3.2's
    // boundary disks emit by construction (V-21) and the weld removes.
}

/// The three surviving non-manifold edges on `csg_difference`, named.
///
/// The census above found six bad edges in the soup where the welded mesh
/// reports three, so which three survive is a matter of record rather than
/// inference. This matches them by position — the welder merges by position, so
/// a welded vertex carries a soup vertex's exact bits — and prints what is
/// standing on each one.
#[test]
fn the_surviving_non_manifold_edges_are_not_duplicated_polygons() {
    use alloc::collections::{BTreeMap, BTreeSet};

    use crate::fields::{ReferenceField, csg_difference};
    use crate::validate::{ValidateConfig, validate_features};

    let field = csg_difference::<f64>();
    let (lo, hi) = field.domain();
    let n = 17u32;
    let cell = (hi[0] - lo[0]) / f64::from(n - 1);
    let (mesh, runs) = soup_with_provenance(&field, lo, cell, n, 16);
    let triangles = mesh.indices.len() / 3;
    let mut owner = alloc::vec![usize::MAX; triangles];
    for (r, run) in runs.iter().enumerate() {
        for slot in owner.iter_mut().skip(run.first).take(run.count) {
            *slot = r;
        }
    }
    let mut by_edge: BTreeMap<[[u64; 3]; 2], Vec<usize>> = BTreeMap::new();
    for tri in 0..triangles {
        let p = triangle_points(&mesh, tri);
        if p[0] == p[1] || p[1] == p[2] || p[2] == p[0] {
            continue;
        }
        for edge in triangle_edges(&mesh, tri) {
            by_edge.entry(edge).or_default().push(tri);
        }
    }

    let shape = RuntimeShape3::new([n; 3]).expect("a cubic grid");
    let mut welded = MeshBuffer::<f64>::default();
    SubgridMarchingTetrahedra::<f64>::new(16)
        .expect("valid")
        .extract(&field, &shape, lo, cell, &mut welded)
        .expect("extraction");
    crate::weld::Welder::<f64>::new()
        .weld(&mut welded, crate::weld::epsilon_for(cell))
        .expect("weld");
    let cfg = ValidateConfig::from_cell_size(cell).expect("a valid spacing");
    let (report, features) = validate_features(&welded.positions, &welded.indices, &cfg);
    assert_eq!(report.non_manifold_edges, 3, "the pinned count moved");

    let mut duplication_only = 0;
    let mut three_distinct = 0;
    for edge in &features.edges {
        let mut key = [
            point_key(welded.positions[edge[0] as usize]),
            point_key(welded.positions[edge[1] as usize]),
        ];
        key.sort_unstable();
        let Some(tris) = by_edge.get(&key) else {
            panic!("a welded non-manifold edge has no soup edge with the same bits");
        };
        let distinct: BTreeSet<[[u64; 3]; 3]> =
            tris.iter().map(|t| triangle_key(&mesh, *t)).collect();
        let tets: BTreeSet<usize> = tris.iter().map(|t| owner[*t]).collect();
        if distinct.len() < 3 {
            duplication_only += 1;
        } else {
            three_distinct += 1;
        }
        std::println!(
            "  edge: {} soup faces, {} distinct polygons, {} tetrahedra {:?}",
            tris.len(),
            distinct.len(),
            tets.len(),
            tets.iter()
                .map(|t| (runs[*t].cell, runs[*t].tet))
                .collect::<Vec<_>>(),
        );
    }

    // The result A-014d turns on. §3.2.3's inset separates two *coincident
    // copies of one polygon*; not one of the surviving defects is that.
    assert_eq!(
        (duplication_only, three_distinct),
        (0, 3),
        "every surviving non-manifold edge should be three distinct polygons meeting"
    );
}

/// Whether `gyroid`'s 138 inconsistently-oriented edges are A-014d's to fix.
///
/// The ticket says they are — *"the same cause wearing a different face"* — on
/// the grounds that 12 triangles have zero area and 24 have the gradient lying
/// in the triangle's plane, so §3.2.3's inset would give them area and a normal
/// transverse to the gradient. The first measurement above found **no
/// §3.2.3 pair on `gyroid` at all**, which makes "the same cause" checkable
/// rather than assumed: this asks what the flipped edges are actually touching.
///
/// A-014e decides winding by `dot(face_normal, gradient(centroid))`, so a
/// triangle whose normal is perpendicular to the gradient has no answer, and a
/// zero-area triangle has no normal. Both would show up here as a `|cos|` at or
/// near zero.
#[test]
fn what_gyroids_flipped_edges_are_standing_on() {
    use alloc::collections::BTreeSet;

    use crate::fields::{ReferenceField, capped_gyroid};
    use crate::validate::{ValidateConfig, validate_features};

    let field = capped_gyroid::<f64>();
    let (lo, hi) = field.domain();
    let n = 17u32;
    let cell = (hi[0] - lo[0]) / f64::from(n - 1);
    let shape = RuntimeShape3::new([n; 3]).expect("a cubic grid");
    let mut mesh = MeshBuffer::<f64>::default();
    SubgridMarchingTetrahedra::<f64>::new(16)
        .expect("valid")
        .extract(&field, &shape, lo, cell, &mut mesh)
        .expect("extraction");
    crate::weld::Welder::<f64>::new()
        .weld(&mut mesh, crate::weld::epsilon_for(cell))
        .expect("weld");
    let cfg = ValidateConfig::from_cell_size(cell).expect("a valid spacing");
    let (report, features) = validate_features(&mesh.positions, &mesh.indices, &cfg);
    assert_eq!(
        report.inconsistently_oriented_edges, 138,
        "the pinned count moved"
    );

    // Every triangle on a flipped edge, once.
    let mut on_flipped = BTreeSet::new();
    for tri in 0..mesh.indices.len() / 3 {
        let v = [
            mesh.indices[tri * 3],
            mesh.indices[tri * 3 + 1],
            mesh.indices[tri * 3 + 2],
        ];
        for k in 0..3 {
            let mut edge = [v[k], v[(k + 1) % 3]];
            edge.sort_unstable();
            if features.inconsistently_oriented_edges.contains(&edge) {
                on_flipped.insert(tri);
            }
        }
    }

    // The two quantities A-014e's vote depends on. A triangle it cannot decide
    // is one with no area or with its normal perpendicular to the gradient.
    let mut undecidable = 0u64;
    let mut zero_area = 0u64;
    let mut decisive = 0u64;
    for tri in &on_flipped {
        let p = [
            mesh.positions[mesh.indices[tri * 3] as usize],
            mesh.positions[mesh.indices[tri * 3 + 1] as usize],
            mesh.positions[mesh.indices[tri * 3 + 2] as usize],
        ];
        let face = crate::vec3::cross(crate::vec3::sub(p[1], p[0]), crate::vec3::sub(p[2], p[0]));
        let area = crate::vec3::length(face) * 0.5;
        let third = 1.0 / 3.0;
        let centroid = [
            (p[0][0] + p[1][0] + p[2][0]) * third,
            (p[0][1] + p[1][1] + p[2][1]) * third,
            (p[0][2] + p[1][2] + p[2][2]) * third,
        ];
        let g = field.gradient(centroid);
        let scale = crate::vec3::length(face) * crate::vec3::length(g);
        if area <= f64::EPSILON * cell * cell {
            zero_area += 1;
        } else if scale > 0.0 && (crate::vec3::dot(face, g) / scale).abs() < 1e-6 {
            undecidable += 1;
        } else {
            decisive += 1;
        }
    }

    std::println!(
        "gyroid: {} triangles on the 138 flipped edges -- zero-area {zero_area}, \
         gradient-in-plane {undecidable}, decisive {decisive}",
        on_flipped.len(),
    );

    // **The ticket's mechanism is real and accounts for 8% of the defect.** 15
    // of the 186 triangles have the gradient in their own plane, so A-014e's
    // vote genuinely had no answer for them and giving them a transverse normal
    // would settle it. The other **171 have a decisive vote and disagree with
    // their neighbour anyway** -- two triangles of one sheet, each correctly
    // oriented against the gradient at its own centroid, wound opposite ways.
    // That is the per-triangle vote of A-014e, not §3.2.3's immersion, and no
    // amount of insetting reaches it.
    //
    // Zero-area is 0 rather than the soup's 12 because the weld removes a
    // triangle with a repeated vertex before this ever sees it.
    assert_eq!(
        (zero_area, undecidable, decisive),
        (0, 15, 171),
        "gyroid's flipped edges moved"
    );
}

/// Which of §3.2's five cases the reference fields actually reach.
///
/// A-014d's shared vertex table has to give every vertex a key two tetrahedra
/// can compute independently, and the only vertices that admits are **top-level
/// edge crossings**: a crossing is the `index`-th root along a tetrahedron edge,
/// the edge runs low-corner to high-corner so its direction is a property of the
/// grid, and `all_roots` is deterministic on bit-identical endpoints. Steiner
/// points have no such key — a centroid is a property of one tetrahedron.
///
/// Subdivision is the case that decides how much of the patch is keyable, since
/// it re-emits the parent's face crossings as *child* crossings with child-local
/// labels. So this counts which cases fire before any of that is designed.
#[test]
fn which_fill_cases_the_reference_fields_reach() {
    use crate::fields::ReferenceField;
    use crate::subgrid::curves::CurveKind;
    use crate::subgrid::surface::{Pattern, cycles, residual};

    let mut rows = 0;
    crate::for_each_reference_field!(f64, |name, field| {
        let (lo, hi) = field.domain();
        let n = 17u32;
        let cell = (hi[0] - lo[0]) / f64::from(n - 1);
        let mt = SubgridMarchingTetrahedra::<f64>::new(16).expect("valid");

        // [corner cuts, non-normal, quads, octagons, single loop, subdivision]
        let mut cases = [0u64; 6];
        let mut tets = 0u64;
        for z in 0..n - 1 {
            for y in 0..n - 1 {
                for x in 0..n - 1 {
                    for tet in TETS {
                        // The same crossings `cell_tet` computes, by the same
                        // expression -- see `tet_corners` for why that matters.
                        let mut corners = [[0.0f64; 3]; 4];
                        for (c, slot) in corners.iter_mut().enumerate() {
                            let offset = corner_offset(tet[c]);
                            for axis in 0..3 {
                                let index = f64::from([x, y, z][axis]) + f64::from(offset[axis]);
                                slot[axis] = lo[axis] + cell * index;
                            }
                        }
                        let mut along: [Vec<f64>; TET_EDGE_COUNT] = Default::default();
                        let mut total = 0usize;
                        for (e, slot) in along.iter_mut().enumerate() {
                            let [a, b] = TET_EDGES[e];
                            super::super::roots::all_roots(
                                corners[a as usize],
                                corners[b as usize],
                                &field,
                                mt.samples(),
                                slot,
                            );
                            total += slot.len();
                        }
                        if total == 0 {
                            continue;
                        }
                        tets += 1;

                        let mut borrowed: [&[f64]; TET_EDGE_COUNT] = [&[]; TET_EDGE_COUNT];
                        for (slot, v) in borrowed.iter_mut().zip(along.iter()) {
                            *slot = v.as_slice();
                        }
                        let crossings = TetCrossings {
                            corners,
                            along: borrowed,
                        };
                        let coords = crossings.coordinates();
                        let all = cycles(&coords);
                        if all.iter().any(super::super::surface::Cycle::is_corner_cut) {
                            cases[0] += 1;
                        }
                        if all.iter().any(|c| c.kind == CurveKind::NonNormal) {
                            cases[1] += 1;
                        }
                        let Some(pattern) = Pattern::of(&residual(&all)) else {
                            continue;
                        };
                        match pattern.loop_length() {
                            Some(4) => cases[2] += 1,
                            Some(8) => cases[3] += 1,
                            Some(_) if pattern.loop_count() == 1 => cases[4] += 1,
                            Some(_) => cases[5] += 1,
                            None => {}
                        }
                    }
                }
            }
        }

        std::println!(
            "{name:<15} {tets:>5} tets | corner cuts {:>5} non-normal {:>4} quads {:>5} \
             octagons {:>4} single loop {:>4} SUBDIVISION {:>4}",
            cases[0],
            cases[1],
            cases[2],
            cases[3],
            cases[4],
            cases[5],
        );
        // **Subdivision never fires on any reference field.** That is what makes
        // a shared vertex table tractable: every vertex on a tetrahedron's
        // boundary is a top-level crossing with a globally computable key, and
        // the only unkeyable positions are Steiner points, which are interior to
        // one tetrahedron and are not shared with anything.
        assert_eq!(cases[5], 0, "{name}: subdivision fired");
        rows += 1;
    });
    assert_eq!(rows, 7, "the sweep did not reach every field");
}

/// **How complete the shared vertex table is, against a positional weld — and
/// the one thing that stops it being complete.**
///
/// M-96 measured that the raw output used to be a per-tetrahedron triangle soup:
/// every edge a boundary edge, no topology to check until a weld by position had
/// run. A-014g gives each crossing a global identity so that stops being true,
/// and this measures how far it got by validating the raw output **directly**
/// and against the same output welded.
///
/// The answer is that identity-based sharing is complete **exactly when no root
/// lands on a grid sample point**, and the correlation is not approximate:
/// `torus` and `fbm_terrain` have zero vertices on a grid point and the weld
/// removes zero, while `box_exact` has *every* vertex on one and the weld removes
/// 924 of 1262. A root at a tetrahedron edge's endpoint has one position and a
/// different `(edge, index)` on each of the up-to-24 tet edges meeting there, so
/// it is one point wearing many names — which no amount of correct sharing under
/// this key can merge.
///
/// That is a property of how the field sits on the grid rather than of the
/// algorithm: `box_exact`'s faces are axis-aligned and land on sample planes,
/// which is M-94's fixture trap showing up structurally rather than by accident.
/// **The fix is still an identity and not a weld** — a crossing at parameter 0
/// or 1 should be named by the grid point it sits on rather than by the edge it
/// was found along — and it is named as the next step rather than done here.
#[test]
fn how_complete_the_shared_table_is_against_a_positional_weld() {
    use crate::fields::ReferenceField;
    use crate::validate::{ValidateConfig, validate_indexed};

    // field -> (raw vertices, welded vertices, vertices on a grid sample point)
    let expected: [(&str, usize, usize, usize); 7] = [
        ("sphere", 812, 812, 6),
        ("torus", 912, 912, 0),
        ("box_exact", 338, 338, 338),
        ("csg_difference", 482, 482, 289),
        ("thin_plate", 422, 422, 30),
        ("gyroid", 4014, 4014, 1),
        ("fbm_terrain", 1758, 1758, 0),
    ];

    let mut rows = 0;
    crate::for_each_reference_field!(f64, |name, field| {
        let (lo, hi) = field.domain();
        let n = 17u32;
        let cell = (hi[0] - lo[0]) / f64::from(n - 1);
        let shape = RuntimeShape3::new([n; 3]).expect("a cubic grid");
        let cfg = ValidateConfig::from_cell_size(cell).expect("a valid spacing");

        let mut raw = MeshBuffer::<f64>::default();
        SubgridMarchingTetrahedra::<f64>::new(16)
            .expect("valid")
            .extract(&field, &shape, lo, cell, &mut raw)
            .unwrap_or_else(|e| panic!("{name}: {e}"));
        let direct = validate_indexed(&raw.positions, &raw.indices, &cfg);

        let mut welded = raw.clone();
        crate::weld::Welder::<f64>::new()
            .weld(&mut welded, crate::weld::epsilon_for(cell))
            .unwrap_or_else(|e| panic!("{name}: weld failed: {e}"));
        let after = validate_indexed(&welded.positions, &welded.indices, &cfg);

        // A root at a tetrahedron edge's *endpoint* sits on a grid sample point,
        // which up to 24 tet edges meet at -- so it has one position and a
        // different key on each of them. Identity-based sharing cannot merge
        // those; a positional weld can. This counts them, to test that
        // mechanism rather than assume it.
        let on_lattice = raw
            .positions
            .iter()
            .filter(|p| {
                (0..3).all(|k| {
                    let g = (p[k] - lo[k]) / cell;
                    (g - g.round()).abs() < 1e-9
                })
            })
            .count();

        std::println!(
            "{name:<15} {:>6} verts -> {:>6} welded ({} removed, {on_lattice} on a grid point) \
             | nm-edges {} vs {} \
             nm-verts {} vs {} flipped {} vs {} boundary {} vs {}",
            raw.positions.len(),
            welded.positions.len(),
            raw.positions.len() - welded.positions.len(),
            direct.non_manifold_edges,
            after.non_manifold_edges,
            direct.non_manifold_vertices,
            after.non_manifold_vertices,
            direct.inconsistently_oriented_edges,
            after.inconsistently_oriented_edges,
            direct.boundary_edges,
            after.boundary_edges,
        );

        let want = expected
            .iter()
            .find(|(f, ..)| *f == name)
            .unwrap_or_else(|| panic!("{name} is not in the pinned table"));
        assert_eq!(
            (raw.positions.len(), welded.positions.len(), on_lattice),
            (want.1, want.2, want.3),
            "{name}"
        );

        // The weld never *adds* a defect: sharing by identity is a subset of
        // sharing by position, so every counter must be at least as good after
        // welding as before. A table that merged two vertices that are not the
        // same point would break this rather than merely under-merge.
        assert!(
            direct.non_manifold_edges >= after.non_manifold_edges,
            "{name}"
        );
        assert!(
            direct.non_manifold_vertices >= after.non_manifold_vertices,
            "{name}"
        );
        assert!(
            direct.inconsistently_oriented_edges >= after.inconsistently_oriented_edges,
            "{name}"
        );
        assert_eq!(direct.boundary_edges, after.boundary_edges, "{name}");

        // **Complete exactly where nothing lands on a grid point.**
        if on_lattice == 0 {
            assert_eq!(
                raw.positions.len(),
                welded.positions.len(),
                "{name}: no vertex is on a grid point, so the weld should find nothing to do"
            );
        }
        rows += 1;
    });
    assert_eq!(rows, 7, "the sweep did not reach every field");
}

/// **What A-014h can actually reach, measured before designing it (M-180).**
///
/// M-169 established that identity-based sharing is incomplete exactly where a
/// root lands on a grid sample point, and named the remedy: name such a crossing
/// by the point rather than by the edge. That remedy assumes the copies *are the
/// same point*, and this measures whether they are.
///
/// They are mostly not. Three counts per field: raw vertices, distinct raw
/// vertices **by bit pattern**, and vertices after a positional weld. The middle
/// column is the ceiling on any rule that keeps positions where the extractor put
/// them — and on `box_exact` it is 1028 against the weld's 338.
///
/// So the gap M-169 measured is not one population but two, and only the smaller
/// one is an identity problem. The rest are positions that a human would call the
/// same grid point and IEEE calls different numbers, because two tetrahedra reach
/// it along different edges and `a + (b − a)·t` rounds differently on each.
/// Merging those is not naming, it is *moving* — which is a decision about
/// geometry and belongs to the ticket, not to a test.
#[test]
fn how_much_of_the_positional_weld_an_exact_identity_could_ever_reach() {
    use std::collections::HashSet;

    // field -> (raw, distinct by bit pattern, welded at crate::weld::epsilon_for(cell))
    let expected: [(&str, usize, usize, usize); 7] = [
        ("sphere", 812, 812, 812),
        ("torus", 912, 912, 912),
        ("box_exact", 338, 338, 338),
        ("csg_difference", 482, 482, 482),
        ("thin_plate", 422, 422, 422),
        ("gyroid", 4014, 4014, 4014),
        ("fbm_terrain", 1758, 1758, 1758),
    ];

    let mut rows = 0;
    let mut got_rows: Vec<(&str, (usize, usize, usize))> = Vec::new();
    crate::for_each_reference_field!(f64, |name, field| {
        let (shape, lo, cell) = grid(&field, 17);
        let mut raw = MeshBuffer::<f64>::default();
        SubgridMarchingTetrahedra::<f64>::new(16)
            .expect("valid")
            .extract(&field, &shape, lo, cell, &mut raw)
            .unwrap_or_else(|e| panic!("{name}: {e}"));

        let distinct: HashSet<[u64; 3]> = raw
            .positions
            .iter()
            .map(|p| [p[0].to_bits(), p[1].to_bits(), p[2].to_bits()])
            .collect();

        let mut welded = raw.clone();
        crate::weld::Welder::<f64>::new()
            .weld(&mut welded, crate::weld::epsilon_for(cell))
            .unwrap_or_else(|e| panic!("{name}: weld failed: {e}"));

        let got = (raw.positions.len(), distinct.len(), welded.positions.len());
        let want = expected
            .iter()
            .find(|(f, ..)| *f == name)
            .map(|&(_, a, b, c)| (a, b, c))
            .unwrap_or_else(|| panic!("{name} is not in the table"));

        std::println!(
            "{name:<15} raw {:>5}  exact-distinct {:>5}  welded {:>5}  \
             | an exact rule could remove {:>4} of the weld's {:>4}",
            got.0,
            got.1,
            got.2,
            got.0 - got.1,
            got.0 - got.2,
        );
        assert_eq!(got, want, "{name}");
        got_rows.push((name, got));

        // An exact rule can never merge more than a tolerance rule, on any
        // field. This is the direction that must hold everywhere.
        assert!(
            got.1 >= got.2,
            "{name}: exact identity out-merged a positional weld, which is \
             impossible -- one of the two is not doing what it says"
        );
        rows += 1;
    });
    assert_eq!(rows, 7);

    // The split, stated as an assertion rather than left in the table. On
    // `gyroid` every merge the weld makes is exact, so an identity rule would
    // close it completely; on the two CSG-shaped fields it would close a quarter
    // at most. Same defect name, two different problems underneath.
    // **A-014h's acceptance, in the strongest form available**: on every field
    // the three columns agree, so the extractor emits each vertex exactly once
    // and a positional weld has nothing left to do. The gap M-180 measured --
    // 690 of `box_exact`'s duplicates unreachable by any keying rule -- closed
    // by canonicalising an endpoint root onto its grid point rather than by
    // finding a cleverer key (M-184).
    for &(name, (raw, exact, welded)) in &got_rows {
        assert_eq!(
            (raw, raw),
            (exact, welded),
            "{name}: identity-based sharing is no longer complete -- \
             raw {raw}, distinct {exact}, welded {welded}"
        );
    }
}

/// **Where an endpoint root actually lands, and why `t == 0` is the wrong test
/// (M-179).**
///
/// A-014h's stated mechanism is that *"a crossing at parameter 0 or 1 should be
/// named by the grid point it lies on"*. Measured over every tetrahedron edge of
/// every cell at 17³, **no root anywhere reports `t == 0`**, and only 36 report
/// `t == 1` — all of them on `gyroid`. The rule as written would find almost
/// nothing.
///
/// The cause is [`refine`](crate::subgrid::roots)'s deliberate choice to return
/// the **upper** end of its final bracket, so it can keep the ascending-and-
/// distinct contract when a root sits on a sample. A root at an edge's lower
/// endpoint therefore comes back as a tiny positive parameter, never `0`.
///
/// What that tiny parameter does to the *position* is the other half, and it is
/// also not what the ticket assumed: `a + (b − a)·t` never rounds back onto
/// `corners[a]` on any field, while it does land exactly on `corners[b]` for a
/// root at the far end. Both directions were expected to behave alike; only one
/// does.
#[test]
fn no_root_reports_parameter_zero_and_almost_none_reports_one() {
    use crate::subgrid::roots::all_roots;

    // field -> (roots, t == 0, t == 1, position == corners[a], == corners[b])
    let expected: [(&str, usize, usize, usize, usize, usize); 7] = [
        ("sphere", 4200, 0, 0, 0, 18),
        ("torus", 4632, 0, 0, 0, 0),
        ("box_exact", 6312, 0, 0, 0, 1692),
        ("csg_difference", 6408, 0, 0, 0, 1284),
        ("thin_plate", 2248, 0, 0, 0, 112),
        ("gyroid", 20352, 0, 36, 0, 36),
        ("fbm_terrain", 8336, 0, 0, 0, 0),
    ];

    let mut rows = 0;
    crate::for_each_reference_field!(f64, |name, field| {
        let (shape, lo, cell) = grid(&field, 17);
        let dims = shape.size();
        let (mut total, mut t0, mut t1, mut at_a, mut at_b) = (0, 0, 0, 0, 0);

        for z in 0..dims[2] - 1 {
            for y in 0..dims[1] - 1 {
                for x in 0..dims[0] - 1 {
                    for tet in &TETS {
                        // The same corner expression `cell_tet` uses -- M-32's
                        // caveat is that equal by algebra is not equal by IEEE.
                        let mut corners = [[0.0f64; 3]; 4];
                        for (c, slot) in corners.iter_mut().enumerate() {
                            let off = corner_offset(tet[c]);
                            for axis in 0..3 {
                                let index = f64::from([x, y, z][axis]) + f64::from(off[axis]);
                                slot[axis] = lo[axis] + cell * index;
                            }
                        }
                        for &[ia, ib] in &TET_EDGES {
                            let (a, b) = (corners[ia as usize], corners[ib as usize]);
                            let mut ts = Vec::new();
                            all_roots(a, b, &field, 16, &mut ts);
                            for &tt in &ts {
                                total += 1;
                                // Bit comparison, not `==`: "exactly the
                                // endpoint parameter" is the whole question, so
                                // a tolerance would answer a different one.
                                if tt.to_bits() == 0.0f64.to_bits() {
                                    t0 += 1;
                                }
                                if tt.to_bits() == 1.0f64.to_bits() {
                                    t1 += 1;
                                }
                                let p = [0, 1, 2].map(|k| a[k] + (b[k] - a[k]) * tt);
                                let same = |u: [f64; 3], v: [f64; 3]| {
                                    (0..3).all(|k| u[k].to_bits() == v[k].to_bits())
                                };
                                if same(p, a) {
                                    at_a += 1;
                                }
                                if same(p, b) {
                                    at_b += 1;
                                }
                            }
                        }
                    }
                }
            }
        }

        let got = (total, t0, t1, at_a, at_b);
        let want = expected
            .iter()
            .find(|(f, ..)| *f == name)
            .map(|&(_, a, b, c, d, e)| (a, b, c, d, e))
            .unwrap_or_else(|| panic!("{name} is not in the table"));
        std::println!(
            "{name:<15} roots {:>6}  t==0 {:>3}  t==1 {:>3}  \
             at corner a {:>3}  at corner b {:>5}",
            got.0,
            got.1,
            got.2,
            got.3,
            got.4
        );
        assert_eq!(got, want, "{name}");
        rows += 1;
    });
    assert_eq!(rows, 7);

    // The headline, asserted rather than left to the table: the ticket's stated
    // test finds nothing on six of seven fields, and nothing at all at the lower
    // endpoint.
    assert!(
        expected
            .iter()
            .all(|&(_, _, t0, _, at_a, _)| t0 == 0 && at_a == 0),
        "a root reported parameter 0, or landed on the lower corner -- \
         M-179's mechanism has changed and A-014h can be re-scoped"
    );
}

/// **How wide the weld tolerance's plateau is (T-009, M-181).**
///
/// T-009 routed four scattered epsilons — `h·1e-4` via the constant, `h·1e-4` as
/// a literal, `h·1e-6` and `h·1e-5` — onto one policy, and expected to re-pin the
/// censuses that moved. **Nothing moved.** This measures why rather than leaving
/// that as luck: the weld's answer is *flat* across the whole range the four
/// policies spanned, because the duplicate clusters sit orders of magnitude
/// closer together than any genuine pair of neighbours.
///
/// That is the fact that made the scatter survivable, and it is also the fact
/// that made it dangerous: four different numbers that happen to agree today
/// disagree the moment a field puts two real vertices `1e-5·h` apart.
///
/// It doubles as **P-7 limb (a)**, whose falsifier is *"a component count that
/// changes with the epsilon"*.
#[test]
fn the_weld_answer_is_flat_across_the_range_the_four_policies_spanned() {
    // The four policies T-009 replaced were 1e-4, 1e-5 and 1e-6; the ends here
    // are a decade beyond each, so the plateau's edges are visible if it has any.
    let factors: [f64; 10] = [1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 0.5];
    let spanned = 3..=5; // the indices covering 1e-6 ..= 1e-4

    let mut rows = 0;
    crate::for_each_reference_field!(f64, |name, field| {
        let (shape, lo, cell) = grid(&field, 17);
        let mut raw = MeshBuffer::<f64>::default();
        SubgridMarchingTetrahedra::<f64>::new(16)
            .expect("valid")
            .extract(&field, &shape, lo, cell, &mut raw)
            .unwrap_or_else(|e| panic!("{name}: {e}"));

        let counts: Vec<usize> = factors
            .iter()
            .map(|f| {
                let mut m = raw.clone();
                crate::weld::Welder::<f64>::new()
                    .weld(&mut m, cell * f)
                    .unwrap_or_else(|e| panic!("{name} at {f:e}: {e}"));
                m.positions.len()
            })
            .collect();

        std::println!("{name:<15} raw {:>5} | {counts:?}", raw.positions.len());

        let over_the_span = &counts[spanned.clone()];
        assert!(
            over_the_span.iter().all(|c| *c == over_the_span[0]),
            "{name}: the weld is NOT flat across the policies T-009 merged \
             ({over_the_span:?}) -- every census that used one of them has to be \
             re-read, and P-7 limb (a) is falsified"
        );

        // The margin, not just the outcome: the first factor whose answer
        // differs from the policy's. P-7 asks for this explicitly, because it is
        // the quantity that would decide a coarser weld.
        let policy = counts[5]; // 1e-4, the one policy
        let first_change = factors
            .iter()
            .zip(&counts)
            .find(|(_, c)| **c != policy)
            .map(|(f, _)| *f);
        match first_change {
            Some(f) => {
                let margin = f / crate::validate::ValidateConfig::WELD_EPSILON_REL;
                std::println!("{name:<15}   first change at {f:e}  ({margin:.0}x the policy)");
                assert!(
                    margin >= 100.0,
                    "{name}: the policy sits within {margin}x of a tolerance that \
                     changes the answer -- too close to be a plateau"
                );
            }
            None => std::println!("{name:<15}   no change anywhere up to 0.5h"),
        }

        rows += 1;
    });
    assert_eq!(rows, 7);
}

/// **P-7's remaining limbs, measured (M-182).**
///
/// P-7 was registered before A-014f was attempted, because A-014f's proposed
/// remedy — orient each connected component from its most confident triangle and
/// propagate — is a *per-patch* decision, and M-96 is on record that
/// `thin_plate`'s two faces sit 0.4 cells apart and land in one tetrahedron
/// facing opposite ways. The worry was that a weld merges the two sheets through
/// their thickness and propagation then silently inverts one of them, on the one
/// field that justifies the whole subgrid track (M-95).
///
/// Limb (a) — the margin — was closed at T-009: the plate's welded vertex count
/// is flat from `h·10⁻⁹` to `h·10⁻³` and first moves at `h·10⁻¹`, 1000× the
/// policy. This is limb (b) and both remaining falsifiers.
///
/// **Limb (b)** is the stronger claim and the one that actually decides A-014f:
/// `thin_plate` is *closed in its domain*, so its top and bottom are joined at a
/// rim and are one component **by topology rather than by tolerance**. A closed
/// orientable surface carries a coherent orientation across that rim, which makes
/// component-wise propagation the right answer there rather than a risk to it.
///
/// Note what this does **not** say. It says propagation is safe on the field the
/// concern was raised about. Whether propagation is the right remedy for
/// `gyroid`'s 138 flipped edges is A-014f's question, and this does not answer it.
#[test]
fn the_plate_is_one_closed_orientable_component_at_every_tolerance() {
    use crate::fields::ThinPlate;
    use crate::validate::{ValidateConfig, validate_indexed};

    let mut flipped = Vec::new();
    for n in [17u32, 25, 33] {
        let field = ThinPlate::<f64>::canonical();
        let (shape, lo, cell) = grid(&field, n);
        let cfg = ValidateConfig::from_cell_size(cell).expect("a valid spacing");

        let mut raw = MeshBuffer::<f64>::default();
        SubgridMarchingTetrahedra::<f64>::new(16)
            .expect("valid")
            .extract(&field, &shape, lo, cell, &mut raw)
            .unwrap_or_else(|e| panic!("{n}: {e}"));

        let mut welded = raw.clone();
        crate::weld::Welder::<f64>::new()
            .weld(&mut welded, crate::weld::epsilon_for(cell))
            .expect("weld");
        let r = validate_indexed(&welded.positions, &welded.indices, &cfg);

        std::println!(
            "thin_plate {n}³  components {}  boundary {}  flipped {}  \
             nm-edges {}  nm-verts {}  chi {}  genus {:?}",
            r.components,
            r.boundary_edges,
            r.inconsistently_oriented_edges,
            r.non_manifold_edges,
            r.non_manifold_vertices,
            r.euler_characteristic,
            r.genus,
        );

        // Falsifier 1: more than one component after the weld.
        assert_eq!(
            r.components, 1,
            "thin_plate {n}³ welds into {} components -- P-7 falsified, and \
             component-wise propagation is unsafe on the field it most needs to \
             be safe on",
            r.components
        );
        // Limb (b): joined at a rim, so one component by topology and not by
        // tolerance. This is the part that makes the margin measured at T-009
        // beside the point rather than merely comfortable.
        assert_eq!(
            r.boundary_edges, 0,
            "thin_plate {n}³ is not closed in its domain -- limb (b)'s premise \
             is wrong and the two sheets are one component only by tolerance"
        );
        assert_eq!(r.non_manifold_edges, 0);
        assert_eq!(r.non_manifold_vertices, 0);
        assert_eq!(r.euler_characteristic, 2, "a closed genus-0 plate");

        // **Recorded, not gated (M-182).** The surface is orientable -- one
        // closed component of genus 0 -- but §3.2's winding, flipped per triangle
        // by A-014e's gradient vote, is not *coherent* on it at every resolution.
        // That is A-014f's defect appearing on the field A-014e's per-triangle
        // rule was written to protect, and it is an argument *for* propagation
        // rather than against it: there is a coherent orientation to propagate to,
        // and the local vote is not finding it.
        flipped.push((n, r.inconsistently_oriented_edges));
    }
    assert_eq!(
        flipped,
        alloc::vec![(17, 0), (25, 8), (33, 6)],
        "thin_plate's flipped-edge census moved"
    );

    // Falsifier 2: a component count that changes with the epsilon. T-009
    // measured the vertex count across the plateau; this is the count P-7 names.
    let field = ThinPlate::<f64>::canonical();
    let (shape, lo, cell) = grid(&field, 17);
    let cfg = ValidateConfig::from_cell_size(cell).expect("a valid spacing");
    let mut raw = MeshBuffer::<f64>::default();
    SubgridMarchingTetrahedra::<f64>::new(16)
        .expect("valid")
        .extract(&field, &shape, lo, cell, &mut raw)
        .expect("extract");

    for factor in [1e-9f64, 1e-7, 1e-5, 1e-4, 1e-3, 1e-2] {
        let mut m = raw.clone();
        crate::weld::Welder::<f64>::new()
            .weld(&mut m, cell * factor)
            .expect("weld");
        let r = validate_indexed(&m.positions, &m.indices, &cfg);
        assert_eq!(
            r.components, 1,
            "thin_plate at {factor:e}·h welds into {} components -- P-7's second \
             falsifier fires and the remedy is epsilon-dependent",
            r.components
        );
    }
}
