//! Tests for A-003.
//!
//! The structural ones matter more here than usual: there is no published table
//! to cross-check against, so the construction has to be checked against its own
//! stated properties rather than against anyone else's numbers.

use alloc::vec::Vec;

use super::MarchingTetrahedra;
use super::table::{TET_CASES, TET_COUNT, TET_EDGES, TETS, tet_edge_corners};
use crate::cube::corner_offset;
use crate::fields::{ReferenceField, Sphere, Torus};
use crate::validate::{ValidateConfig, check_determinism, validate_indexed};
use crate::{MeshBuffer, RuntimeShape3, Sdf, vec3};

// ─── the decomposition ──────────────────────────────────────────────────────

/// The six tetrahedra are the six axis orderings, generated independently of
/// `build_tets`'s hard-coded list.
#[test]
fn the_six_tetrahedra_are_the_six_axis_orderings() {
    let mut expected: Vec<[u8; 4]> = Vec::new();
    for a in 0..3u8 {
        for b in 0..3u8 {
            for c in 0..3u8 {
                if a == b || b == c || a == c {
                    continue;
                }
                let p1 = 1 << a;
                let p2 = p1 | (1 << b);
                expected.push([0, p1, p2, 7]);
            }
        }
    }
    let mut got: Vec<[u8; 4]> = TETS.to_vec();
    expected.sort_unstable();
    got.sort_unstable();
    assert_eq!(got, expected);
    assert_eq!(got.len(), TET_COUNT);
}

/// They tile the cube: six tetrahedra of volume 1/6 each, disjoint interiors.
///
/// Checked by integer volume rather than by geometry — the signed volume of a
/// tetrahedron on `0/1` corners is `±1/6`, so six of them summing to one unit
/// cube with no overlap means every determinant has magnitude one.
#[test]
fn the_six_tetrahedra_tile_the_cube() {
    let mut total = 0i32;
    for (t, tet) in TETS.iter().enumerate() {
        let p: Vec<[i32; 3]> = tet
            .iter()
            .map(|&c| corner_offset(c).map(|v| v as i32))
            .collect();
        let u = [p[1][0] - p[0][0], p[1][1] - p[0][1], p[1][2] - p[0][2]];
        let v = [p[2][0] - p[0][0], p[2][1] - p[0][1], p[2][2] - p[0][2]];
        let w = [p[3][0] - p[0][0], p[3][1] - p[0][1], p[3][2] - p[0][2]];
        let det = u[0] * (v[1] * w[2] - v[2] * w[1]) - u[1] * (v[0] * w[2] - v[2] * w[0])
            + u[2] * (v[0] * w[1] - v[1] * w[0]);
        assert_eq!(det.abs(), 1, "tetrahedron {t} is degenerate or overlapping");
        total += det.abs();
    }
    // Six unit-determinant tetrahedra = 6 * (1/6) = one cube.
    assert_eq!(total, 6);
}

/// **The crack-free property, and the reason this decomposition rather than the
/// five-tetrahedron one.**
///
/// Two cells adjacent along an axis split their shared face by a diagonal each.
/// Those two diagonals must be the same *world* segment, or the surfaces do not
/// meet. Checked on all three axes by mapping both cells' corners into a common
/// frame.
#[test]
fn every_shared_face_is_split_the_same_way_by_both_cells() {
    // The diagonal of each cube face, as the pair of corners the tetrahedra
    // share on it.
    let face_diagonal = |axis: usize, side: u8| -> [u8; 2] {
        let mut found = Vec::new();
        for t in 0..TET_COUNT {
            for e in 0..TET_EDGES.len() {
                let [a, b] = tet_edge_corners(t, e);
                let pa = corner_offset(a);
                let pb = corner_offset(b);
                // On this face, and a diagonal of it rather than an edge.
                if pa[axis] as u8 == side
                    && pb[axis] as u8 == side
                    && (a ^ b).count_ones() == 2
                    && !found.contains(&[a.min(b), a.max(b)])
                {
                    found.push([a.min(b), a.max(b)]);
                }
            }
        }
        assert_eq!(found.len(), 1, "axis {axis} side {side}: {found:?}");
        found[0]
    };

    for axis in 0..3usize {
        // This cell's far face, and the neighbour's near face along `axis`.
        let near = face_diagonal(axis, 0);
        let far = face_diagonal(axis, 1);

        // Put both in world coordinates. The neighbour sits one step along
        // `axis`, so its corner offsets gain one there.
        let mut step = [0i32; 3];
        step[axis] = 1;
        let world_far: Vec<[i32; 3]> = far
            .iter()
            .map(|&c| corner_offset(c).map(|v| v as i32))
            .collect();
        let mut world_near: Vec<[i32; 3]> = near
            .iter()
            .map(|&c| {
                let o = corner_offset(c).map(|v| v as i32);
                [o[0] + step[0], o[1] + step[1], o[2] + step[2]]
            })
            .collect();
        let mut world_far = world_far;
        world_near.sort_unstable();
        world_far.sort_unstable();
        assert_eq!(
            world_near, world_far,
            "axis {axis}: this cell splits its far face on {far:?} and the next cell splits \
             its near face on {near:?}, and those are different segments -- a crack"
        );
    }
}

// ─── the case table ─────────────────────────────────────────────────────────

/// Every case is a triangle, a quad, or nothing. **A tetrahedron cannot be
/// ambiguous**, and this is that claim as an assertion: three cut edges or four,
/// never anything else.
#[test]
fn every_case_cuts_three_edges_or_four_or_none() {
    let mut histogram = [0usize; 5];
    for (t, cases) in TET_CASES.iter().enumerate() {
        for (case, entry) in cases.iter().enumerate() {
            let cut = (0..TET_EDGES.len())
                .filter(|&e| {
                    let [a, b] = TET_EDGES[e];
                    (case >> a) & 1 != (case >> b) & 1
                })
                .count();
            assert!(
                matches!(cut, 0 | 3 | 4),
                "tet {t} case {case} cuts {cut} edges"
            );
            let expected = match cut {
                0 => 0,
                3 => 1,
                _ => 2,
            };
            assert_eq!(entry.count as usize, expected, "tet {t} case {case}");
            histogram[cut.min(4)] += 1;
        }
    }
    std::println!("measured: cut-edge histogram over 96 (tet, case) pairs -> {histogram:?}");
}

/// Triangles only ever name cut edges, and never repeat one.
#[test]
fn triangles_name_only_cut_edges() {
    for (t, cases) in TET_CASES.iter().enumerate() {
        for (case, entry) in cases.iter().enumerate() {
            for tri in &entry.triangles[..entry.count as usize] {
                assert!(
                    tri[0] != tri[1] && tri[1] != tri[2] && tri[0] != tri[2],
                    "tet {t} case {case} repeats an edge"
                );
                for &e in tri {
                    let [a, b] = TET_EDGES[e as usize];
                    assert_ne!(
                        (case >> a) & 1,
                        (case >> b) & 1,
                        "tet {t} case {case} names uncut edge {e}"
                    );
                }
            }
        }
    }
}

/// **Winding, which no manifold or Euler check can see.**
///
/// Every triangle's normal must point away from the tetrahedron's inside
/// corners. Recomputed here from the corner geometry, independently of the
/// `orient` that built the table.
#[test]
fn every_triangle_faces_away_from_the_solid() {
    for (t, cases) in TET_CASES.iter().enumerate() {
        for (case, entry) in cases.iter().enumerate().take(15).skip(1) {
            for tri in &entry.triangles[..entry.count as usize] {
                let p: Vec<[f64; 3]> = tri
                    .iter()
                    .map(|&e| {
                        let [a, b] = tet_edge_corners(t, e as usize);
                        let pa = corner_offset(a);
                        let pb = corner_offset(b);
                        [
                            f64::from(pa[0] + pb[0]) * 0.5,
                            f64::from(pa[1] + pb[1]) * 0.5,
                            f64::from(pa[2] + pb[2]) * 0.5,
                        ]
                    })
                    .collect();
                let normal = vec3::cross(vec3::sub(p[1], p[0]), vec3::sub(p[2], p[0]));
                // Centroid of the inside corners: the solid side.
                let mut inside = [0.0f64; 3];
                let mut n = 0.0;
                for (i, &corner) in TETS[t].iter().enumerate() {
                    if case & (1 << i) != 0 {
                        let c = corner_offset(corner);
                        for (slot, v) in inside.iter_mut().zip(c) {
                            *slot += f64::from(v);
                        }
                        n += 1.0;
                    }
                }
                for slot in &mut inside {
                    *slot /= n;
                }
                let away = vec3::sub(p[0], inside);
                assert!(
                    vec3::dot(normal, away) > 0.0,
                    "tet {t} case {case} triangle {tri:?} is wound inward"
                );
            }
        }
    }
}

// ─── meshes ─────────────────────────────────────────────────────────────────

fn mesh<F: Sdf<Scalar = f64> + ReferenceField>(field: &F, samples: u32) -> (MeshBuffer<f64>, f64) {
    let (lo, hi) = field.domain();
    let cell_size = (hi[0] - lo[0]) / f64::from(samples - 1);
    let mut mt = MarchingTetrahedra::<f64>::new();
    let mut out = MeshBuffer::<f64>::new();
    let shape = RuntimeShape3::new([samples; 3]).expect("valid shape");
    mt.extract(field, &shape, lo, cell_size, &mut out)
        .expect("extraction");
    (out, cell_size)
}

#[test]
fn a_meshed_sphere_is_closed() {
    for samples in [17u32, 25, 33] {
        let field = Sphere::<f64>::canonical();
        let (out, h) = mesh(&field, samples);
        let report = validate_indexed(
            &out.positions,
            &out.indices,
            &ValidateConfig::from_cell_size(h).expect("valid cell size"),
        );
        assert!(out.triangle_count() > 0);
        assert!(report.is_closed(), "{samples} samples:\n{report}");
        assert_eq!(report.euler_characteristic, 2, "{samples}:\n{report}");
        assert_eq!(report.non_manifold_edges, 0, "{samples}:\n{report}");
        assert_eq!(report.boundary_edges, 0, "{samples}:\n{report}");
        assert_eq!(
            report.inconsistently_oriented_edges, 0,
            "{samples}:\n{report}"
        );
    }
}

/// Signed volume catches a global inversion, which no topology check can.
#[test]
fn a_meshed_sphere_has_positive_signed_volume() {
    let (out, _) = mesh(&Sphere::<f64>::canonical(), 33);
    let mut total = 0.0;
    for t in out.indices.chunks_exact(3) {
        let a = out.positions[t[0] as usize];
        let b = out.positions[t[1] as usize];
        let c = out.positions[t[2] as usize];
        total += vec3::dot(a, vec3::cross(b, c));
    }
    let volume = total / 6.0;
    let exact = 4.0 / 3.0 * core::f64::consts::PI;
    assert!(volume > 0.0, "inside out: {volume}");
    assert!((volume - exact).abs() / exact < 0.02, "{volume} vs {exact}");
}

#[test]
fn a_meshed_torus_has_genus_one() {
    let (out, h) = mesh(&Torus::<f64>::canonical(), 49);
    let report = validate_indexed(
        &out.positions,
        &out.indices,
        &ValidateConfig::from_cell_size(h).expect("valid cell size"),
    );
    assert!(report.is_closed(), "{report}");
    assert_eq!(report.euler_characteristic, 0, "{report}");
    assert_eq!(report.genus, Some(1), "{report}");
}

#[test]
fn extraction_is_deterministic() {
    let field = Sphere::<f64>::canonical();
    let shape = RuntimeShape3::new([25; 3]).expect("valid shape");
    let report = check_determinism(|out: &mut MeshBuffer<f64>| {
        MarchingTetrahedra::<f64>::new()
            .extract(&field, &shape, [-2.0; 3], 4.0 / 24.0, out)
            .expect("extraction");
    });
    assert!(report.is_deterministic(), "{report}");
}

#[test]
fn every_closed_reference_field_meshes_cleanly() {
    crate::for_each_reference_field!(f64, |name, field| {
        for samples in [17u32, 33] {
            let (lo, hi) = field.domain();
            let h = (hi[0] - lo[0]) / f64::from(samples - 1);
            let shape = RuntimeShape3::new([samples; 3]).expect("valid shape");
            let mut out = MeshBuffer::<f64>::new();
            MarchingTetrahedra::<f64>::new()
                .extract(&field, &shape, lo, h, &mut out)
                .expect("extraction");
            if out.triangle_count() == 0 {
                continue;
            }
            let report = validate_indexed(
                &out.positions,
                &out.indices,
                &ValidateConfig::from_cell_size(h).expect("valid cell size"),
            );
            if field.closed_in_domain() {
                assert!(report.is_closed(), "{name} at {samples}^3:\n{report}");
            } else {
                assert!(report.is_manifold(), "{name} at {samples}^3:\n{report}");
            }
            if let Some(chi) = field.expected_euler() {
                assert_eq!(
                    report.euler_characteristic, chi,
                    "{name} at {samples}^3:\n{report}"
                );
            }
        }
    });
}

// ─── the two claims this ticket exists to check ─────────────────────────────

/// **Claim 1 (tier R): tetrahedral methods give "2–3× more vertices".**
///
/// From the v1 catalog's Family A table, sourced to `10.1109/2945.485620`. The
/// ticket's own wording is softer — "expect a large triangle count; record it" —
/// so this records it and pins the range, in both directions: it fails if the
/// cost vanishes as well as if it grows.
#[test]
fn the_triangle_count_against_marching_cubes_is_measured() {
    use crate::marching_cubes::MarchingCubes;

    let mut worst_v = 0.0f64;
    let mut worst_t = 0.0f64;
    crate::for_each_reference_field!(f64, |name, field| {
        for samples in [33u32, 49] {
            let (lo, hi) = field.domain();
            let h = (hi[0] - lo[0]) / f64::from(samples - 1);
            let shape = RuntimeShape3::new([samples; 3]).expect("valid shape");

            let mut mc = MeshBuffer::<f64>::new();
            MarchingCubes::<f64>::new()
                .extract(&field, &shape, lo, h, &mut mc)
                .expect("extraction");
            let mut mt = MeshBuffer::<f64>::new();
            MarchingTetrahedra::<f64>::new()
                .extract(&field, &shape, lo, h, &mut mt)
                .expect("extraction");
            if mc.vertex_count() == 0 {
                continue;
            }
            let v = mt.vertex_count() as f64 / mc.vertex_count() as f64;
            let t = mt.triangle_count() as f64 / mc.triangle_count() as f64;
            worst_v = worst_v.max(v);
            worst_t = worst_t.max(t);
            std::println!(
                "measured: {name} at {samples}^3 -- vertices {} vs {} ({v:.2}x), \
                 triangles {} vs {} ({t:.2}x)",
                mt.vertex_count(),
                mc.vertex_count(),
                mt.triangle_count(),
                mc.triangle_count()
            );
        }
    });
    std::println!("measured: worst ratios -- vertices {worst_v:.2}x, triangles {worst_t:.2}x");
    assert!(worst_v > 1.0, "tetrahedra must cost something");
    assert!(
        worst_v < 4.0,
        "vertex ratio {worst_v} is outside anything reported"
    );
}

/// **Claim 2 (primary source): tetrahedral methods are *geometrically* worse,
/// not merely bulkier.**
///
/// Lewiner et al. 2003 (`10.1080/10867651.2003.10487582`), the Marching Cubes 33
/// paper already in the corpus: *"They generate many more triangles, with a
/// weaker geometrical accuracy of the result: the tetrahedra's tilings are
/// segmented even in obvious configuration, and the vertex position cannot be
/// adjusted to fit the geometrical trilinear approximation as we do with
/// cubes."*
///
/// That is a stronger and more interesting claim than the triangle count, and it
/// is testable with T-003's harness on the field M-10 already measured Marching
/// Cubes and Surface Nets on. Recorded rather than asserted as a bound, because
/// the number is the finding.
#[test]
fn the_accuracy_against_marching_cubes_is_measured() {
    use crate::marching_cubes::MarchingCubes;
    use crate::validate::{AccuracyConfig, accuracy};

    let field = Sphere::<f64>::canonical();
    let samples = 64u32;
    let (lo, hi) = field.domain();
    let h = (hi[0] - lo[0]) / f64::from(samples - 1);
    let shape = RuntimeShape3::new([samples; 3]).expect("valid shape");
    let cfg = AccuracyConfig::from_cell_size(h).expect("valid cell size");

    let mut mc = MeshBuffer::<f64>::new();
    MarchingCubes::<f64>::new()
        .extract(&field, &shape, lo, h, &mut mc)
        .expect("extraction");
    let mut mt = MeshBuffer::<f64>::new();
    MarchingTetrahedra::<f64>::new()
        .extract(&field, &shape, lo, h, &mut mt)
        .expect("extraction");

    let mc_report =
        accuracy(&mc.positions, &mc.indices, &field, &shape, lo, &cfg).expect("accuracy");
    let mt_report =
        accuracy(&mt.positions, &mt.indices, &field, &shape, lo, &cfg).expect("accuracy");
    let mc_h = mc_report.symmetric_hausdorff();
    let mt_h = mt_report.symmetric_hausdorff();

    std::println!(
        "measured: unit sphere at {samples}^3 -- symmetric Hausdorff, marching cubes {mc_h:.4e}, \
         marching tetrahedra {mt_h:.4e}, ratio {:.3}x",
        mt_h / mc_h
    );
    assert!(mt_h.is_finite() && mc_h.is_finite());
}

/// **Orientation does not move the ratio, and curvature does not either.**
///
/// P-1 predicts `2.992` by weighting the seven edge families by `E[|n·e|]` over
/// uniformly random surface orientations. The seven reference fields measure
/// `2.86–3.91×`, so something separates them, and the two obvious candidates are
/// both wrong — recorded here because a falsified hypothesis is worth more than
/// an unexamined one.
///
/// **Orientation**: a plane at four orientations, from axis-aligned to generic,
/// gives `3.919 / 3.939 / 3.945 / 3.943`. Flat within a percent.
///
/// **Curvature**: a sphere swept from radius `0.3` to `1.8` on the same grid
/// gives `3.036 / 3.046 / 3.026 / 2.995 / 2.981` — converging *down* onto P-1's
/// `2.99` as it flattens, when a locally flatter sphere ought to approach the
/// plane's `3.94` if flatness were the variable.
///
/// So a plane sits at `3.94` and a sphere at `3.0` while being locally the same
/// shape at cell scale, and **nothing measured here explains the gap**. What can
/// be said is empirical: flat-faced fields (`box_exact` 3.91, `thin_plate` 3.84,
/// `csg_difference` 3.83) sit at the plane end, smooth closed ones (`sphere` and
/// `torus` 3.04) at P-1's value, and the rough high-genus ones (`gyroid` and
/// `fbm_terrain` 2.87) below it. See O-15.
#[test]
fn orientation_does_not_move_the_vertex_ratio() {
    use crate::marching_cubes::MarchingCubes;

    /// A plane through the origin with the given unit normal.
    struct Plane([f64; 3]);
    impl Sdf for Plane {
        type Scalar = f64;
        fn sample(&self, p: [f64; 3]) -> f64 {
            p[0] * self.0[0] + p[1] * self.0[1] + p[2] * self.0[2]
        }
        fn gradient(&self, _p: [f64; 3]) -> [f64; 3] {
            self.0
        }
    }

    let samples = 49u32;
    let h = 4.0 / f64::from(samples - 1);
    let shape = RuntimeShape3::new([samples; 3]).expect("valid shape");
    let inv3 = 1.0 / 3.0f64.sqrt();
    let inv2 = 1.0 / 2.0f64.sqrt();

    let mut ratios = Vec::new();
    for (name, normal) in [
        ("axis-aligned  (0,0,1)", [0.0, 0.0, 1.0]),
        ("face diagonal (0,1,1)", [0.0, inv2, inv2]),
        ("body diagonal (1,1,1)", [inv3, inv3, inv3]),
        ("generic", [0.3714, 0.5571, 0.7428]),
    ] {
        let field = Plane(normal);
        let mut mc = MeshBuffer::<f64>::new();
        MarchingCubes::<f64>::new()
            .extract(&field, &shape, [-2.0; 3], h, &mut mc)
            .expect("extraction");
        let mut mt = MeshBuffer::<f64>::new();
        MarchingTetrahedra::<f64>::new()
            .extract(&field, &shape, [-2.0; 3], h, &mut mt)
            .expect("extraction");
        let ratio = mt.vertex_count() as f64 / mc.vertex_count() as f64;
        ratios.push((name, ratio));
        std::println!(
            "measured: plane {name} -- marching cubes {} verts, marching tetrahedra {} ({ratio:.3}x)",
            mc.vertex_count(),
            mt.vertex_count()
        );
    }

    // Pinned as *insensitivity*: every orientation within 1% of every other.
    // This fails if orientation ever starts mattering, which would be the thing
    // that explains the gap.
    let lo = ratios.iter().map(|(_, r)| *r).fold(f64::INFINITY, f64::min);
    let hi = ratios.iter().map(|(_, r)| *r).fold(0.0f64, f64::max);
    assert!(
        hi / lo < 1.01,
        "orientation moved the ratio: {lo:.3} to {hi:.3}"
    );
    // And all of them sit at the flat-field end, well above P-1's 2.99.
    assert!(lo > 3.8, "a plane should cost ~3.9x, got {lo:.3}");
}
