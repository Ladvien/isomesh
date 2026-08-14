//! Tests for Marching Cubes.
//!
//! Three independent lines of defence on the case table, because a wrong entry
//! produces a mesh that looks fine:
//!
//! 1. [`super::validate_table`] — structural, consulting no reference at all.
//! 2. `matches_the_published_table` — against Bourke's independently parsed
//!    copy of the classic table, with the numbering correspondence *derived*
//!    rather than assumed.
//! 3. The meshes themselves, through the T-001/T-002/T-004 harnesses.

use alloc::vec::Vec;
use alloc::{format, vec};

use super::reference::{BOURKE_EDGE_TABLE, BOURKE_TRI_TABLE};
use super::table::{
    AMBIGUOUS_FACES, CASES, EDGE_CORNERS, EDGE_COUNT, MAX_TRIANGLES, NO_EDGE, corner_inside,
    face_bit, face_corners, is_inside, segment_links, triangulate,
};
use super::{FaceAmbiguity, MarchingCubes, validate_decider_table, validate_table};
use crate::fields::{BoxExact, CappedGyroid, CsgDifference, ReferenceField, Sphere, Torus};
use crate::validate::{ValidateConfig, check_determinism, self_intersections, validate_indexed};
use crate::{MeshBuffer, RuntimeShape3, Sdf, vec3};

// ─── the derived table, on its own terms ────────────────────────────────────

#[test]
fn the_derived_table_is_structurally_sound() {
    let report = validate_table();
    assert!(report.is_sound(), "{report:?}");
    assert_eq!(report.face_disagreements, 0, "cracks: {report:?}");
}

/// Recorded, and it is what sets [`MAX_TRIANGLES`]: the construction's real
/// worst case is five, which matches the classic table's well-known bound.
#[test]
fn no_case_exceeds_the_triangle_bound() {
    let report = validate_table();
    assert_eq!(report.max_triangles, 5);
    assert!(report.max_triangles as usize <= MAX_TRIANGLES);
}

#[test]
fn empty_and_full_cases_produce_nothing() {
    assert_eq!(CASES[0].count, 0);
    assert_eq!(CASES[255].count, 0);
}

/// A single inside corner must give one triangle whose normal points away from
/// that corner. This is the smallest case where a global winding flip is
/// visible, and no manifold or Euler check can see one.
#[test]
fn a_single_inside_corner_faces_outward() {
    let entry = &CASES[1];
    assert_eq!(entry.count, 1);

    let midpoint = |e: u8| {
        let [a, b] = EDGE_CORNERS[e as usize];
        let pa = crate::cube::corner_offset(a);
        let pb = crate::cube::corner_offset(b);
        [
            f64::from(pa[0] + pb[0]) * 0.5,
            f64::from(pa[1] + pb[1]) * 0.5,
            f64::from(pa[2] + pb[2]) * 0.5,
        ]
    };
    let v: Vec<[f64; 3]> = entry.triangles[0].iter().map(|&e| midpoint(e)).collect();
    let normal = vec3::cross(vec3::sub(v[1], v[0]), vec3::sub(v[2], v[0]));
    let centroid = [
        (v[0][0] + v[1][0] + v[2][0]) / 3.0,
        (v[0][1] + v[1][1] + v[2][1]) / 3.0,
        (v[0][2] + v[1][2] + v[2][2]) / 3.0,
    ];
    // Corner 0 is at the origin and is the only solid one, so "away from the
    // solid" is the direction of the centroid.
    assert!(
        vec3::dot(normal, centroid) > 0.0,
        "winding is inward: normal {normal:?}, centroid {centroid:?}"
    );
}

// ─── against the published table ────────────────────────────────────────────

/// Recover which corner pair each of Bourke's edges joins, from the XOR
/// structure of his edge table alone.
///
/// Edge `j` is cut exactly when its two corners are classified differently, so
/// `bit_j(edgeTable[c]) == bit_a(c) ^ bit_b(c)` for every one of the 256 cases.
/// Only one pair satisfies that for all of them. Nothing here relies on reading
/// his diagram, which is the step that would reintroduce transcription risk.
fn derive_bourke_edge_corners() -> [[u8; 2]; EDGE_COUNT] {
    let mut out = [[0u8; 2]; EDGE_COUNT];
    for (j, slot) in out.iter_mut().enumerate() {
        let mut found = None;
        for a in 0..8u8 {
            for b in (a + 1)..8u8 {
                let ok = (0..256usize).all(|c| {
                    let cut = (BOURKE_EDGE_TABLE[c] >> j) & 1 == 1;
                    let differ = ((c >> a) & 1) != ((c >> b) & 1);
                    cut == differ
                });
                if ok {
                    assert!(
                        found.is_none(),
                        "edge {j} matched more than one corner pair"
                    );
                    found = Some([a, b]);
                }
            }
        }
        *slot = found.unwrap_or_else(|| panic!("edge {j} matched no corner pair"));
    }
    out
}

/// The permutation carrying Bourke's corner labels onto ours, found by search
/// rather than assumed. The cube graph has 48 automorphisms, so several
/// permutations qualify; any one of them is a valid relabelling.
fn derive_corner_permutation(bourke: &[[u8; 2]; EDGE_COUNT]) -> [u8; 8] {
    let ours: Vec<(u8, u8)> = EDGE_CORNERS.iter().map(|e| (e[0], e[1])).collect();
    let mut perm = [0u8; 8];
    let mut used = [false; 8];
    fn search(
        depth: usize,
        perm: &mut [u8; 8],
        used: &mut [bool; 8],
        bourke: &[[u8; 2]; EDGE_COUNT],
        ours: &[(u8, u8)],
    ) -> bool {
        if depth == 8 {
            return bourke.iter().all(|e| {
                let (a, b) = (perm[e[0] as usize], perm[e[1] as usize]);
                ours.contains(&(a.min(b), a.max(b)))
            });
        }
        for candidate in 0..8u8 {
            if used[candidate as usize] {
                continue;
            }
            used[candidate as usize] = true;
            perm[depth] = candidate;
            if search(depth + 1, perm, used, bourke, ours) {
                return true;
            }
            used[candidate as usize] = false;
        }
        false
    }
    assert!(
        search(0, &mut perm, &mut used, bourke, &ours),
        "no relabelling maps Bourke's cube onto ours"
    );
    perm
}

/// The derived table against the published one, case by case.
///
/// Comparison is on the **boundary of the triangle set** rather than the
/// triangles themselves, because the same polygon can be fanned several ways and
/// a fan difference is not a disagreement about the surface. Interior edges
/// cancel; what remains is the loop structure, which is the surface.
#[test]
fn matches_the_published_table() {
    let bourke_edges = derive_bourke_edge_corners();
    let perm = derive_corner_permutation(&bourke_edges);

    let our_edge = |a: u8, b: u8| -> u8 {
        let (lo, hi) = if a < b { (a, b) } else { (b, a) };
        EDGE_CORNERS
            .iter()
            .position(|e| e[0] == lo && e[1] == hi)
            .expect("adjacent corners") as u8
    };

    let boundary = |tris: &[[u8; 3]]| -> Vec<(u8, u8)> {
        let mut directed: Vec<(u8, u8)> = Vec::new();
        for t in tris {
            directed.push((t[0], t[1]));
            directed.push((t[1], t[2]));
            directed.push((t[2], t[0]));
        }
        let mut out: Vec<(u8, u8)> = directed
            .iter()
            .filter(|(a, b)| !directed.contains(&(*b, *a)))
            .copied()
            .collect();
        out.sort_unstable();
        out
    };

    let mut agree = 0usize;
    for (bourke_case, row) in BOURKE_TRI_TABLE.iter().enumerate() {
        // Bourke's case index, relabelled into ours.
        let mut our_case = 0u8;
        for (corner, &mapped_corner) in perm.iter().enumerate() {
            if (bourke_case >> corner) & 1 == 1 {
                our_case |= 1 << mapped_corner;
            }
        }

        let mut theirs: Vec<[u8; 3]> = Vec::new();
        let mut i = 0;
        while i + 2 < 16 && row[i] >= 0 {
            let mapped = |k: usize| {
                let e = bourke_edges[row[i + k] as usize];
                our_edge(perm[e[0] as usize], perm[e[1] as usize])
            };
            theirs.push([mapped(0), mapped(1), mapped(2)]);
            i += 3;
        }

        let entry = &CASES[our_case as usize];
        let ours: Vec<[u8; 3]> = entry.triangles[..entry.count as usize].to_vec();

        let ours_boundary = boundary(&ours);
        let theirs_reversed: Vec<(u8, u8)> = {
            let mut v: Vec<(u8, u8)> = boundary(&theirs).into_iter().map(|(a, b)| (b, a)).collect();
            v.sort_unstable();
            v
        };
        assert_eq!(
            ours_boundary, theirs_reversed,
            "case {bourke_case} (ours {our_case}) disagrees with the published table"
        );
        agree += 1;
    }
    assert_eq!(agree, 256);
}

/// Our winding is the reverse of Bourke's, and that is a convention rather than
/// a defect — `a_single_inside_corner_faces_outward` and the signed-volume test
/// establish that ours is the outward one.
#[test]
fn our_winding_is_the_reverse_of_the_published_tables() {
    let bourke_edges = derive_bourke_edge_corners();
    let perm = derive_corner_permutation(&bourke_edges);
    // Case 1 in Bourke's numbering: his corner 0 alone is inside.
    let mut our_case = 0u8;
    our_case |= 1 << perm[0];
    let entry = &CASES[our_case as usize];
    assert_eq!(entry.count, 1);
}

// ─── meshes ─────────────────────────────────────────────────────────────────

fn mesh<F: Sdf<Scalar = f64> + ReferenceField>(field: &F, samples: u32) -> (MeshBuffer<f64>, f64) {
    let (lo, hi) = field.domain();
    let cell_size = (hi[0] - lo[0]) / f64::from(samples - 1);
    let mut mc = MarchingCubes::<f64>::new();
    let mut out = MeshBuffer::<f64>::new();
    let shape = RuntimeShape3::new([samples; 3]).expect("valid shape");
    mc.extract(field, &shape, lo, cell_size, &mut out)
        .expect("extraction");
    (out, cell_size)
}

/// Signed volume by the divergence theorem. Positive for an outward-wound closed
/// surface in a right-handed system, so this is the check that catches a global
/// inversion — which every manifold and Euler test passes happily.
fn signed_volume(mesh: &MeshBuffer<f64>) -> f64 {
    let mut total = 0.0;
    for t in mesh.indices.chunks_exact(3) {
        let a = mesh.positions[t[0] as usize];
        let b = mesh.positions[t[1] as usize];
        let c = mesh.positions[t[2] as usize];
        total += vec3::dot(a, vec3::cross(b, c));
    }
    total / 6.0
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
        assert_eq!(
            report.euler_characteristic, 2,
            "{samples} samples:\n{report}"
        );
        assert_eq!(report.genus, Some(0));
        assert_eq!(report.boundary_edges, 0);
        assert_eq!(report.non_manifold_edges, 0);
        assert_eq!(report.non_manifold_vertices, 0);
        assert_eq!(report.inconsistently_oriented_edges, 0);
    }
}

#[test]
fn meshed_sphere_has_positive_signed_volume() {
    let field = Sphere::<f64>::canonical();
    let (out, _) = mesh(&field, 33);
    let volume = signed_volume(&out);
    let exact = 4.0 / 3.0 * core::f64::consts::PI;
    assert!(volume > 0.0, "inside out: signed volume {volume}");
    // A polyhedron inscribed in the sphere is a little smaller than it.
    assert!(
        (volume - exact).abs() / exact < 0.02,
        "volume {volume} vs exact {exact}"
    );
}

/// A torus is genus 1, so this is what proves the extractor reproduces topology
/// rather than merely closing up.
#[test]
fn a_meshed_torus_has_genus_one() {
    let field = Torus::<f64>::canonical();
    let (out, h) = mesh(&field, 49);
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
fn sharp_and_concave_fields_stay_manifold() {
    let box_field = BoxExact::<f64>::canonical();
    let (out, h) = mesh(&box_field, 33);
    let report = validate_indexed(
        &out.positions,
        &out.indices,
        &ValidateConfig::from_cell_size(h).expect("valid cell size"),
    );
    assert!(report.is_closed(), "box_exact:\n{report}");

    let csg: CsgDifference<f64> = crate::fields::csg_difference();
    let (out, h) = mesh(&csg, 41);
    let report = validate_indexed(
        &out.positions,
        &out.indices,
        &ValidateConfig::from_cell_size(h).expect("valid cell size"),
    );
    assert!(report.is_closed(), "csg_difference:\n{report}");
    assert_eq!(report.euler_characteristic, 2, "csg_difference:\n{report}");
}

/// The high-genus case. Its Euler characteristic is not known analytically, so
/// it is recorded rather than asserted — exactly what `expected_euler() == None`
/// is telling the harness to do.
#[test]
fn the_capped_gyroid_is_closed_and_its_genus_is_recorded() {
    let field: CappedGyroid<f64> = crate::fields::capped_gyroid();
    assert_eq!(field.expected_euler(), None);

    let (out, h) = mesh(&field, 49);
    let report = validate_indexed(
        &out.positions,
        &out.indices,
        &ValidateConfig::from_cell_size(h).expect("valid cell size"),
    );
    assert!(report.is_closed(), "{report}");
    assert_eq!(
        report.euler_characteristic % 2,
        0,
        "chi must be even for any closed orientable surface\n{report}"
    );
    std::println!(
        "measured: capped gyroid at 49^3 -> chi = {}, genus = {:?}, {} triangles",
        report.euler_characteristic,
        report.genus,
        out.triangle_count()
    );
}

#[test]
fn marching_cubes_is_deterministic() {
    let field = Sphere::<f64>::canonical();
    let (lo, hi) = field.domain();
    let cell_size = (hi[0] - lo[0]) / 24.0;
    let shape = RuntimeShape3::new([25; 3]).expect("valid shape");
    let mut mc = MarchingCubes::<f64>::new();

    let report = check_determinism(|out: &mut MeshBuffer<f64>| {
        mc.extract(&field, &shape, lo, cell_size, out)
            .expect("extraction");
    });
    assert!(report.is_deterministic(), "{report}");
    assert!(report.triangles > 0);
}

/// Marching Cubes does not self-intersect when no grid sample lands *exactly*
/// on the isosurface.
///
/// 27 samples over `[-2, 2]` gives `h = 2/13`, and no lattice point can satisfy
/// `x² + y² + z² = 1` there: every coordinate is an even multiple of `1/13`, so
/// the squares would have to be even integers summing to the odd number 169.
/// The rate is recorded either way, since it is the baseline A-009 compares
/// dual contouring's clamp against.
#[test]
fn a_meshed_sphere_does_not_self_intersect() {
    let field = Sphere::<f64>::canonical();
    let (out, h) = mesh(&field, 27);
    let si = self_intersections(&out.positions, &out.indices, h).expect("self intersections");
    let report = validate_indexed(
        &out.positions,
        &out.indices,
        &ValidateConfig::from_cell_size(h).expect("valid cell size"),
    );
    std::println!(
        "measured: marching cubes on sphere at 27^3 -> {:.3} intersecting pairs per 1000 triangles, \
         {} degenerate triangles",
        si.per_thousand_triangles(),
        report.degenerate_triangles
    );
    assert_eq!(report.degenerate_triangles, 0, "{report}");
    assert!(si.is_intersection_free(), "{si}");
}

/// A grid sample landing exactly on the isosurface produces degenerate triangles
/// and a small number of self-intersections — and the topology survives both.
///
/// This is a property of the sampling, not a defect in the extractor. Zero is
/// classified outside (see [`is_inside`]), so an edge running into an
/// exactly-zero corner interpolates to `t = 1` and places its vertex precisely
/// on that corner. Several grid edges meet there, so several distinct vertices
/// land on the same point and the triangles between them have no area.
///
/// At 25 samples over `[-2, 2]` the spacing is `1/6` and 30 lattice points lie
/// exactly on the unit sphere. The mesh is still closed with `χ = 2`, which is
/// why degenerate triangles are a recorded metric rather than a gate.
#[test]
fn samples_exactly_on_the_surface_produce_slivers_but_not_holes() {
    let field = Sphere::<f64>::canonical();
    let (out, h) = mesh(&field, 25);
    let report = validate_indexed(
        &out.positions,
        &out.indices,
        &ValidateConfig::from_cell_size(h).expect("valid cell size"),
    );
    let si = self_intersections(&out.positions, &out.indices, h).expect("self intersections");

    std::println!(
        "measured: sphere at 25^3 (30 lattice points exactly on the surface) -> \
         {} degenerate triangles, {:.3} intersecting pairs per 1000, chi = {}",
        report.degenerate_triangles,
        si.per_thousand_triangles(),
        report.euler_characteristic
    );

    assert!(report.degenerate_triangles > 0, "premise: slivers appear");
    // The topology is untouched by them.
    assert!(report.is_closed(), "{report}");
    assert_eq!(report.euler_characteristic, 2, "{report}");
}

/// Every reference field, through the gate its own metadata selects. This is the
/// composition the whole of Phase 0 was building toward.
///
/// The sweep is the macro, not a hand-kept list — a hand enumeration silently
/// omitted `thin_plate`, the one closed field this extractor was suspected of
/// failing, and it passes. The resolutions are the decider twin's, and all odd,
/// which is load-bearing for `thin_plate`: its only inside samples lie on the
/// y = 0 midplane, which is a sample plane exactly when `n` is odd. At even `n`
/// plain MC correctly produces an empty mesh (the field's own resolution
/// caveat), so an even resolution here would trip the produced-nothing assert
/// rather than silently passing. Higher-resolution coverage the hand list
/// carried lives on in the per-field tests (`a_meshed_torus_has_genus_one`,
/// `sharp_and_concave_fields_stay_manifold`,
/// `the_capped_gyroid_is_closed_and_its_genus_is_recorded`).
#[test]
fn every_closed_reference_field_meshes_cleanly() {
    crate::for_each_reference_field!(f64, |name, field| {
        for samples in [17u32, 25, 33] {
            let (out, h) = mesh(&field, samples);
            let report = validate_indexed(
                &out.positions,
                &out.indices,
                &ValidateConfig::from_cell_size(h).expect("valid cell size"),
            );
            assert!(
                out.triangle_count() > 0,
                "{name} at {samples}^3 produced nothing"
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
            std::println!(
                "measured: {name} at {samples}^3 -> {} tris, chi {}, {} degenerate",
                out.triangle_count(),
                report.euler_characteristic,
                report.degenerate_triangles
            );
        }
    });
}

#[test]
fn f32_and_f64_both_extract() {
    let mut mc = MarchingCubes::<f32>::new();
    let mut out = MeshBuffer::<f32>::new();
    let shape = RuntimeShape3::new([17; 3]).expect("valid shape");
    mc.extract(
        &Sphere::<f32>::canonical(),
        &shape,
        [-2.0; 3],
        4.0 / 16.0,
        &mut out,
    )
    .expect("extraction");
    assert!(out.triangle_count() > 0);
}

#[test]
fn a_field_with_no_surface_produces_no_triangles() {
    // A sphere far outside the sampled box.
    let field = Sphere::<f64> {
        center: [100.0; 3],
        radius: 1.0,
    };
    let mut mc = MarchingCubes::<f64>::new();
    let mut out = MeshBuffer::<f64>::new();
    let shape = RuntimeShape3::new([9; 3]).expect("valid shape");
    mc.extract(&field, &shape, [-2.0; 3], 0.5, &mut out)
        .expect("extraction");
    assert_eq!(out.triangle_count(), 0);
    assert_eq!(out.vertex_count(), 0);
}

#[test]
fn a_degenerate_grid_is_rejected() {
    let mut mc = MarchingCubes::<f64>::new();
    let mut out = MeshBuffer::<f64>::new();
    let shape = RuntimeShape3::new([1, 4, 4]).expect("valid shape");
    let error = mc
        .extract(
            &Sphere::<f64>::canonical(),
            &shape,
            [-2.0; 3],
            0.5,
            &mut out,
        )
        .expect_err("a one-sample axis contains no cell");
    assert_eq!(error, crate::Error::GridTooSmall { size: [1, 4, 4] });
    assert!(out.is_empty(), "nothing should have been written");
}

/// Vertices are shared between neighbouring cells; without that the output would
/// be a triangle soup and every edge would be a boundary edge.
#[test]
fn vertices_are_shared_between_cells() {
    let field = Sphere::<f64>::canonical();
    let (out, _) = mesh(&field, 25);
    let soup = out.triangle_count() * 3;
    assert!(
        out.vertex_count() < soup / 2,
        "{} vertices for {} triangles looks like a soup",
        out.vertex_count(),
        soup
    );
}

#[test]
fn the_table_reports_render() {
    let report = validate_table();
    let text = format!("{report:?}");
    assert!(text.contains("max_triangles"));
    let _ = vec![NO_EDGE, MAX_TRIANGLES as u8];
    assert!(is_inside(-1.0f64));
    assert!(!is_inside(0.0f64));
    assert!(corner_inside(0b0000_0001, 0));
}

// ─── A-002: the asymptotic decider ──────────────────────────────────────────
//
// The face rule and its own arithmetic are tested in `marching_cubes::ambiguity::tests`.
// What is tested here is the table it feeds and the meshes that come out.

/// Every case at every resolution mask — 16,384 combinations, against the same
/// five structural properties the 256-entry table is held to.
///
/// The face-locality check is the load-bearing one and it is now keyed on the
/// face's decision bit as well as its corner signs. That is still the crack-free
/// property, because the decision is a function of the four sample values the
/// two cells share.
#[test]
fn the_decider_table_is_structurally_sound() {
    let report = validate_decider_table();
    assert!(report.is_sound(), "{report:?}");
    assert_eq!(report.face_disagreements, 0, "cracks: {report:?}");
}

/// The shipped table is the mask-zero construction. This is one of the two
/// assertions that make the table lookup in `extract` a memo rather than a
/// second rule.
#[test]
fn the_separate_mask_reproduces_the_derived_table() {
    for case in 0..=255u8 {
        let built = triangulate(segment_links(case, 0));
        assert_eq!(built.count, CASES[case as usize].count, "case {case}");
        assert_eq!(
            built.triangles, CASES[case as usize].triangles,
            "case {case}"
        );
    }
}

/// And the other: on a case with no ambiguous face, no mask changes anything.
/// `validate_decider_table` proves this per face; this proves it per case, which
/// is the form `extract` actually relies on.
#[test]
fn masks_are_ignored_on_unambiguous_faces() {
    let mut unambiguous = 0usize;
    for case in 0..=255u8 {
        if AMBIGUOUS_FACES[case as usize] != 0 {
            continue;
        }
        unambiguous += 1;
        for mask in 0..64u8 {
            assert_eq!(
                segment_links(case, mask),
                segment_links(case, 0),
                "case {case} mask {mask}"
            );
        }
    }
    // Recorded: how much of the table the memo covers.
    std::println!("measured: {unambiguous} of 256 cases have no ambiguous face");
    assert!(unambiguous > 0);
}

/// `AMBIGUOUS_FACES` is derived from the whole case; the validator's own test
/// derives ambiguity from a single face's 4-bit pattern. They must agree, or one
/// of the two is describing a different set of faces than it claims.
#[test]
fn ambiguous_faces_agrees_with_the_face_pattern() {
    for case in 0..=255u8 {
        let mut expected = 0u8;
        for axis in 0..3usize {
            for side in 0..2u8 {
                let c = face_corners(axis, side);
                let mut pattern = 0usize;
                for (k, &corner) in c.iter().enumerate() {
                    if corner_inside(case, corner) {
                        pattern |= 1 << k;
                    }
                }
                if pattern == 0b0101 || pattern == 0b1010 {
                    expected |= face_bit(axis, side);
                }
            }
        }
        assert_eq!(AMBIGUOUS_FACES[case as usize], expected, "case {case}");
    }
}

/// Recorded, not gated, exactly as `no_case_exceeds_the_triangle_bound` records
/// the separated table's five — which A-015 left untouched, because plain
/// Marching Cubes never produces a cycle long enough to need a centroid.
///
/// **Twelve.** A-002 predicted ten from the longest possible cycle, all twelve
/// cut edges, fanned from one of its own vertices. A-015 replaced that fan with
/// a centroid one precisely because such a cycle has no chord-safe apex, and a
/// centroid fan emits one triangle per cycle edge rather than `k − 2`. So the
/// bound rose by exactly the two triangles the centroid costs, and it now sits
/// on `MAX_TRIANGLES` itself.
#[test]
fn the_decider_does_not_exceed_the_triangle_bound() {
    let report = validate_decider_table();
    std::println!(
        "measured: the decider's worst case is {} triangles (separated table: {})",
        report.max_triangles,
        validate_table().max_triangles
    );
    assert!(
        report.max_triangles as usize <= MAX_TRIANGLES,
        "{} exceeds MAX_TRIANGLES = {MAX_TRIANGLES}",
        report.max_triangles
    );
    assert_eq!(report.max_triangles, 12);
}

/// **A-002's acceptance criterion.** A cell configuration with an ambiguous face
/// on which the two rules report different Euler characteristics.
///
/// The fixture is *searched*, not chosen: this repo has twice had a fixture
/// picked for looking like it exercised a property sit in the region where the
/// property does not apply (M-32, M-38). Every one of the 256 cases is tried,
/// and the census of which ones differ is printed alongside the assertion.
#[test]
fn the_decider_and_marching_cubes_disagree_about_chi() {
    let mut differing = Vec::new();
    let mut first: Option<(u8, i64, i64)> = None;

    for case in 0..=255u8 {
        if AMBIGUOUS_FACES[case as usize] == 0 {
            continue;
        }
        // Deep inside, shallow outside: d_in = 4 against d_out = 1, so every
        // ambiguous face of this case comes out joined.
        let field = one_cell(case, -2.0, 1.0);
        let separate = chi_of(&field, FaceAmbiguity::Separate);
        let decided = chi_of(&field, FaceAmbiguity::AsymptoticDecider);
        if separate != decided {
            differing.push(case);
            if first.is_none() {
                first = Some((case, separate, decided));
            }
        }
    }

    std::println!(
        "measured: {} of 256 cases change chi under the decider at d_in > d_out; \
         first is case {:?}",
        differing.len(),
        first
    );

    let (case, separate, decided) = first.expect("no case changed chi — the decider does nothing");
    // Pinned in both directions: this fails if the difference disappears and if
    // either value moves. Case 6 is corners 1 and 2 inside — diagonally opposite
    // on the z = 0 face, and on no other face together, so it has exactly one
    // ambiguous face. Separated, that face cuts off each corner with its own
    // triangle: two discs, chi = 2. Joined, the two become one disc, chi = 1.
    assert_eq!(case, 0b0000_0110);
    assert_eq!(separate, 2);
    assert_eq!(decided, 1);
    assert_eq!(differing.len(), 88);
}

/// Every reference field under the decider, through the gate its own metadata
/// selects — with one property held out.
///
/// **`is_closed()` is deliberately not the gate here, and manifoldness is not
/// asserted in this test.** `is_closed()` folds in manifoldness, and
/// manifoldness under the decider is owned by the *fan*, not by the ambiguity
/// rule — see `the_fan_lets_adjacent_cells_share_an_interior_chord`, which shows
/// plain Marching Cubes has the identical defect at the identical rate. Asserting
/// it here would attribute A-001's problem to A-002. It is pinned instead, in
/// both directions and per field, by `the_decider_non_manifold_census_is_pinned`.
///
/// What is asserted is everything the ambiguity rule *is* responsible for:
/// crack-freeness on closed fields, consistent orientation, and the published
/// Euler characteristic. Following M-16, no even-χ parity check is made here,
/// because parity is a corollary of manifoldness and manifoldness is held out.
#[test]
fn every_closed_reference_field_meshes_cleanly_under_the_decider() {
    crate::for_each_reference_field!(f64, |name, field| {
        for samples in [17u32, 25, 33] {
            let Some((decided, separate, tris, mc_tris)) = decider_pair(&field, samples) else {
                continue;
            };

            // The gate is the field's own metadata, never the test's opinion.
            // `fbm_terrain` leaves through the sides by construction, so its
            // boundary edges are the field, not a crack.
            if field.closed_in_domain() {
                assert_eq!(
                    decided.boundary_edges, 0,
                    "{name} at {samples}^3 has cracks:\n{decided}"
                );
            }
            assert_eq!(
                decided.inconsistently_oriented_edges, 0,
                "{name} at {samples}^3:\n{decided}"
            );
            if let Some(chi) = field.expected_euler() {
                assert_eq!(
                    decided.euler_characteristic, chi,
                    "{name} at {samples}^3:\n{decided}"
                );
            }
            std::println!(
                "measured: {name} at {samples}^3 -> chi {} (mc {}), tris {tris} (mc {mc_tris}), \
                 non-manifold edges {} (mc {}), boundary {} (mc {})",
                decided.euler_characteristic,
                separate.euler_characteristic,
                decided.non_manifold_edges,
                separate.non_manifold_edges,
                decided.boundary_edges,
                separate.boundary_edges,
            );
        }
    });
}

/// The whole non-manifold census, pinned in both directions.
///
/// **Empty, since A-015.** It was `[("gyroid", 25, 2, 0)]` — the fan chord
/// collision, reached at one resolution of one field — and the chord-safe apex
/// rule removed it. Kept as an assertion rather than deleted, following M-4: the
/// interesting failure now is a non-manifold mesh reappearing, and this is what
/// would see it.
#[test]
fn the_decider_non_manifold_census_is_pinned() {
    let mut offenders: Vec<(&str, u32, u64, u64)> = Vec::new();
    crate::for_each_reference_field!(f64, |name, field| {
        for samples in [17u32, 25, 33] {
            let Some((decided, separate, _, _)) = decider_pair(&field, samples) else {
                continue;
            };
            if decided.non_manifold_edges != 0 || separate.non_manifold_edges != 0 {
                offenders.push((
                    name,
                    samples,
                    decided.non_manifold_edges,
                    separate.non_manifold_edges,
                ));
            }
        }
    });
    std::println!("measured: non-manifold census (field, n, decider, mc) = {offenders:?}");
    assert_eq!(offenders, Vec::new());
}

/// How often the rule even fires. Custodio et al. (2013) report that real-world
/// data is dominated by the unambiguous configurations; this is that claim
/// against this crate's own seven fields, which is the only version of it that
/// can justify a decision here.
#[test]
fn the_ambiguous_face_census_is_recorded() {
    crate::for_each_reference_field!(f64, |name, field| {
        let samples = 33u32;
        let (lo, hi) = field.domain();
        let h = (hi[0] - lo[0]) / f64::from(samples - 1);
        let census = census(&field, lo, h, samples);
        std::println!(
            "measured: {name} at {samples}^3 -> {} surface cells, {} with an ambiguous face \
             ({:.3}%), {} faces joined by the decider",
            census.surface_cells,
            census.ambiguous_cells,
            100.0 * census.ambiguous_cells as f64 / census.surface_cells.max(1) as f64,
            census.joined_faces,
        );
    });
}

#[test]
fn the_decider_is_deterministic() {
    let field = crate::fields::capped_gyroid::<f64>();
    let (lo, hi) = field.domain();
    let samples = 33u32;
    let h = (hi[0] - lo[0]) / f64::from(samples - 1);
    let shape = RuntimeShape3::new([samples; 3]).expect("valid shape");
    let report = check_determinism(|out: &mut MeshBuffer<f64>| {
        let mut mc = MarchingCubes::<f64>::new();
        mc.set_face_ambiguity(FaceAmbiguity::AsymptoticDecider);
        mc.extract(&field, &shape, lo, h, out).expect("extract");
    });
    assert!(report.is_deterministic(), "{report}");
}

/// The decider joins where plain Marching Cubes separates, so the surfaces it
/// produces on ambiguous faces are different ones — self-intersections are a
/// recorded metric, not a gate, and this records them.
#[test]
fn the_decider_self_intersection_count_is_recorded() {
    let field = crate::fields::capped_gyroid::<f64>();
    let (lo, hi) = field.domain();
    let samples = 33u32;
    let h = (hi[0] - lo[0]) / f64::from(samples - 1);
    let shape = RuntimeShape3::new([samples; 3]).expect("valid shape");

    let mut out = MeshBuffer::<f64>::new();
    let mut mc = MarchingCubes::<f64>::new();
    mc.set_face_ambiguity(FaceAmbiguity::AsymptoticDecider);
    mc.extract(&field, &shape, lo, h, &mut out)
        .expect("extract");

    let report = self_intersections(&out.positions, &out.indices, h).expect("grid fits");
    std::println!(
        "measured: gyroid at {samples}^3 under the decider -> {:.4} self-intersections per 1k \
         triangles",
        report.per_thousand_triangles()
    );
}

// ─── A-002 helpers ──────────────────────────────────────────────────────────

/// The trilinear interpolant of a small value grid, as a field.
///
/// The decider's whole claim is agreement with the trilinear interpolant, so a
/// fixture built to test it should *be* one. Sampling at a grid point returns
/// that point's stored value exactly, and the gradient is analytic — nothing
/// approximate sits between the fixture and the result.
struct Trilinear {
    size: [usize; 3],
    values: Vec<f64>,
}

impl Trilinear {
    /// Cell index and local coordinate for one axis, clamped so the interpolant
    /// extends smoothly past the grid rather than folding.
    fn split(&self, p: f64, axis: usize) -> (usize, f64) {
        let last = self.size[axis] - 2;
        let i = libm::floor(p);
        let i = if i < 0.0 {
            0
        } else if i as usize > last {
            last
        } else {
            i as usize
        };
        (i, p - i as f64)
    }

    fn corner(&self, base: [usize; 3], d: [usize; 3]) -> f64 {
        let x = base[0] + d[0];
        let y = base[1] + d[1];
        let z = base[2] + d[2];
        self.values[x + y * self.size[0] + z * self.size[0] * self.size[1]]
    }
}

impl Sdf for Trilinear {
    type Scalar = f64;

    fn sample(&self, p: [f64; 3]) -> f64 {
        let (i, u) = self.split(p[0], 0);
        let (j, v) = self.split(p[1], 1);
        let (k, w) = self.split(p[2], 2);
        let mut total = 0.0;
        for (a, wu) in [(0usize, 1.0 - u), (1, u)] {
            for (b, wv) in [(0usize, 1.0 - v), (1, v)] {
                for (c, ww) in [(0usize, 1.0 - w), (1, w)] {
                    total += self.corner([i, j, k], [a, b, c]) * wu * wv * ww;
                }
            }
        }
        total
    }

    fn gradient(&self, p: [f64; 3]) -> [f64; 3] {
        let (i, u) = self.split(p[0], 0);
        let (j, v) = self.split(p[1], 1);
        let (k, w) = self.split(p[2], 2);
        let mut g = [0.0f64; 3];
        for (a, wu, du) in [(0usize, 1.0 - u, -1.0), (1, u, 1.0)] {
            for (b, wv, dv) in [(0usize, 1.0 - v, -1.0), (1, v, 1.0)] {
                for (c, ww, dw) in [(0usize, 1.0 - w, -1.0), (1, w, 1.0)] {
                    let value = self.corner([i, j, k], [a, b, c]);
                    g[0] += value * du * wv * ww;
                    g[1] += value * wu * dv * ww;
                    g[2] += value * wu * wv * dw;
                }
            }
        }
        g
    }
}

/// A single cell whose corner signs are `case`, inside corners at `inside` and
/// outside ones at `outside`.
fn one_cell(case: u8, inside: f64, outside: f64) -> Trilinear {
    let mut values = vec![0.0f64; 8];
    for (c, slot) in values.iter_mut().enumerate() {
        *slot = if corner_inside(case, c as u8) {
            inside
        } else {
            outside
        };
    }
    Trilinear {
        size: [2, 2, 2],
        values,
    }
}

/// The Euler characteristic of one cell's patch under a given rule.
fn chi_of(field: &Trilinear, rule: FaceAmbiguity) -> i64 {
    let shape = RuntimeShape3::new([
        field.size[0] as u32,
        field.size[1] as u32,
        field.size[2] as u32,
    ])
    .expect("valid shape");
    let mut out = MeshBuffer::<f64>::new();
    let mut mc = MarchingCubes::<f64>::new();
    mc.set_face_ambiguity(rule);
    mc.extract(field, &shape, [0.0; 3], 1.0, &mut out)
        .expect("extract");
    validate_indexed(
        &out.positions,
        &out.indices,
        &ValidateConfig::from_cell_size(1.0).expect("valid cell size"),
    )
    .euler_characteristic
}

/// One field at one resolution, meshed both ways: `(decider report, plain MC
/// report, decider triangles, plain MC triangles)`, or `None` where the field
/// produces no surface at that resolution.
fn decider_pair<F: Sdf<Scalar = f64> + ReferenceField>(
    field: &F,
    samples: u32,
) -> Option<(
    crate::validate::MeshReport,
    crate::validate::MeshReport,
    usize,
    usize,
)> {
    let (lo, hi) = field.domain();
    let h = (hi[0] - lo[0]) / f64::from(samples - 1);
    let shape = RuntimeShape3::new([samples; 3]).expect("valid shape");
    let cfg = ValidateConfig::from_cell_size(h).expect("valid cell size");

    let mut out = MeshBuffer::<f64>::new();
    let mut mc = MarchingCubes::<f64>::new();
    mc.set_face_ambiguity(FaceAmbiguity::AsymptoticDecider);
    mc.extract(field, &shape, lo, h, &mut out).expect("extract");
    if out.triangle_count() == 0 {
        return None;
    }

    let mut plain = MeshBuffer::<f64>::new();
    MarchingCubes::<f64>::new()
        .extract(field, &shape, lo, h, &mut plain)
        .expect("extract");

    Some((
        validate_indexed(&out.positions, &out.indices, &cfg),
        validate_indexed(&plain.positions, &plain.indices, &cfg),
        out.triangle_count(),
        plain.triangle_count(),
    ))
}

#[derive(Default)]
struct Census {
    surface_cells: u64,
    ambiguous_cells: u64,
    joined_faces: u64,
}

/// Walk the same grid `extract` walks and count how often the decider has
/// anything to decide, and how often it answers "joined".
fn census<F: Sdf<Scalar = f64>>(field: &F, origin: [f64; 3], h: f64, samples: u32) -> Census {
    let mut out = Census::default();
    for z in 0..samples - 1 {
        for y in 0..samples - 1 {
            for x in 0..samples - 1 {
                let mut corner_value = [0.0f64; 8];
                let mut case = 0u8;
                for (c, slot) in corner_value.iter_mut().enumerate() {
                    let o = crate::cube::corner_offset(c as u8);
                    *slot = field.sample([
                        origin[0] + h * f64::from(x + o[0]),
                        origin[1] + h * f64::from(y + o[1]),
                        origin[2] + h * f64::from(z + o[2]),
                    ]);
                    if is_inside(*slot) {
                        case |= 1 << c;
                    }
                }
                if case == 0 || case == 255 {
                    continue;
                }
                out.surface_cells += 1;
                let ambiguous = AMBIGUOUS_FACES[case as usize];
                if ambiguous == 0 {
                    continue;
                }
                out.ambiguous_cells += 1;
                out.joined_faces += u64::from(
                    crate::marching_cubes::ambiguity::joined_mask(&corner_value, ambiguous)
                        .count_ones(),
                );
            }
        }
    }
    out
}

/// Every undirected mesh edge one cell contributes, in **global** terms.
///
/// A cut cube edge is identified by its lower sample and its axis, so two cells
/// naming the same cube edge produce the same key. A **cycle centroid** is keyed
/// by the cell itself, because that is exactly what it is: cell-local, nameable
/// by no other cell, which is why it cannot collide (A-015).
fn global_mesh_edges(base: [i64; 3], case: u8, mask: u8) -> Vec<[(i64, i64, i64, u8); 2]> {
    let entry = triangulate(segment_links(case, mask));
    let key = |code: u8| -> (i64, i64, i64, u8) {
        if crate::marching_cubes::table::is_centroid(code) {
            return (base[0], base[1], base[2], code);
        }
        let [c, _] = EDGE_CORNERS[code as usize];
        let o = crate::cube::corner_offset(c);
        (
            base[0] + i64::from(o[0]),
            base[1] + i64::from(o[1]),
            base[2] + i64::from(o[2]),
            super::table::EDGE_AXIS[code as usize],
        )
    };
    let mut out = Vec::new();
    for t in &entry.triangles[..entry.count as usize] {
        for k in 0..3 {
            let a = key(t[k]);
            let b = key(t[(k + 1) % 3]);
            out.push(if a <= b { [a, b] } else { [b, a] });
        }
    }
    out
}

/// **No two adjacent cells can put more than two triangles on one mesh edge.**
///
/// This is A-015's acceptance criterion, and it is exhaustive rather than
/// sampled: two cells stacked along z share the face `z = 1` and have twelve
/// samples between them, so all 4,096 sign patterns fit in a loop, and every
/// canonical resolution mask is tried on each.
///
/// Before A-015 this test recorded the defect instead — **12 of 4,096 patterns
/// putting four triangles on one mesh edge, identically under Marching Cubes'
/// separated rule and under the decider**, which is what established that the
/// fan owned it rather than the ambiguity rule (✗17). It now records its
/// absence, in both directions: `worst == 2` fails if a collision comes back and
/// if the search stops reaching two-cell configurations at all.
#[test]
fn adjacent_cells_never_share_a_mesh_edge_beyond_two_faces() {
    for (label, joined) in [("separate (A-001)", false), ("decider (A-002)", true)] {
        let mut worst = 0usize;
        let mut occurrences = 0usize;
        let mut first: Option<(u32, u8, u8)> = None;

        for bits in 0..(1u32 << 12) {
            // Sample (x, y, z), z in 0..3, at bit x + 2y + 4z.
            let inside = |x: u32, y: u32, z: u32| bits >> (x + 2 * y + 4 * z) & 1 == 1;
            let mut all = Vec::new();
            let mut cases = [0u8; 2];
            for z0 in 0..2u32 {
                let mut case = 0u8;
                for k in 0..8u8 {
                    let o = crate::cube::corner_offset(k);
                    if inside(o[0], o[1], o[2] + z0) {
                        case |= 1 << k;
                    }
                }
                cases[z0 as usize] = case;
                if case == 0 || case == 255 {
                    continue;
                }
                let mask = if joined {
                    AMBIGUOUS_FACES[case as usize]
                } else {
                    0
                };
                all.extend(global_mesh_edges([0, 0, i64::from(z0)], case, mask));
            }

            all.sort_unstable();
            let mut i = 0usize;
            while i < all.len() {
                let mut j = i;
                while j < all.len() && all[j] == all[i] {
                    j += 1;
                }
                worst = worst.max(j - i);
                if j - i > 2 {
                    occurrences += 1;
                    if first.is_none() {
                        first = Some((bits, cases[0], cases[1]));
                    }
                }
                i = j;
            }
        }

        std::println!(
            "measured: {label}: worst faces on one mesh edge = {worst}, {occurrences} of 4096 \
             two-cell sign patterns affected, first at {first:?}"
        );
        assert_eq!(worst, 2, "{label}: {first:?}");
        assert_eq!(occurrences, 0, "{label}");
    }
}

/// **Does the asymptotic decider widen M-32's chunk-seam problem?**
///
/// M-32 measured that two adjacent chunks compute their shared sample plane as
/// `(o + h·cn) + h·n` and `o + h·(c+1)n` — equal by algebra, not by IEEE — and
/// that 22% of random `(origin, h, cells, chunk)` combinations disagree by an
/// ulp. Under plain Marching Cubes that is a rounding error: the sample values
/// shift by an ulp and the vertex on the seam moves by a fraction of a cell.
///
/// The decider raises the stakes, because it turns those values into a
/// **discrete** choice. A flipped choice is not a shifted vertex — it is two
/// chunks building genuinely different surfaces on a face they share, which is a
/// crack. A-002 creates that question and must not leave it to inference.
///
/// The faces examined are the ones that actually straddle a seam: all four
/// corners lie *in* the shared plane, so both chunks read all four through the
/// arithmetic that differs. Planes where the two expressions agree bit for bit
/// are skipped, since they measure nothing — M-32's own method rule.
///
/// Two things are counted, because they fail differently. A **sign
/// disagreement** — the ulp moving a corner across zero — is a crack under
/// *plain* Marching Cubes too, since it changes the case index. A **decision
/// flip** is the one the decider adds. The closest margin any ambiguous seam
/// face came to flipping is recorded, because a count of zero says nothing about
/// how nearly it happened.
#[test]
fn the_decider_at_a_chunk_seam_is_measured() {
    let field = crate::fields::capped_gyroid::<f64>();

    // A plain xorshift, so the sweep is reproducible without a dependency.
    let mut state = 0x2545_f491_4f6c_dd1du64;
    let mut next = || {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        (state >> 11) as f64 / (1u64 << 53) as f64
    };

    let mut divergent_planes = 0usize;
    let mut seam_faces = 0usize;
    let mut sign_disagreements = 0usize;
    let mut ambiguous_faces = 0usize;
    let mut decision_flips = 0usize;
    let mut closest_margin = f64::INFINITY;

    for _ in 0..8_000 {
        let o = -4.0 + 8.0 * next();
        let cells = 4 + (next() * 13.0) as u32;
        // Coarse relative to the gyroid's own period, which is where an
        // ambiguous face is reachable at all — at fine spacings the surface is
        // resolved and the configuration essentially does not occur.
        let h = 0.4 + next();
        let chunk = (next() * 16.0) as i64 - 8;

        let span = h * f64::from(cells);
        let plane_a = (o + span * chunk as f64) + span;
        let plane_b = o + span * (chunk + 1) as f64;
        if plane_a.to_bits() == plane_b.to_bits() || plane_a.abs() > 6.0 {
            continue;
        }
        divergent_planes += 1;

        let corners = |x: f64, y: f64, z: f64| {
            [
                field.sample([x, y, z]),
                field.sample([x, y + h, z]),
                field.sample([x, y + h, z + h]),
                field.sample([x, y, z + h]),
            ]
        };
        let ambiguous = |v: &[f64; 4]| {
            is_inside(v[0]) == is_inside(v[2])
                && is_inside(v[1]) == is_inside(v[3])
                && is_inside(v[0]) != is_inside(v[1])
        };

        let steps = 24i32;
        for iy in -steps..steps {
            for iz in -steps..steps {
                let y = f64::from(iy) * h;
                let z = f64::from(iz) * h;
                let va = corners(plane_a, y, z);
                let vb = corners(plane_b, y, z);
                seam_faces += 1;

                if (0..4).any(|k| is_inside(va[k]) != is_inside(vb[k])) {
                    // Worse than a flipped decider, and not the decider's doing:
                    // the two chunks disagree about the case index itself.
                    sign_disagreements += 1;
                    continue;
                }
                if !ambiguous(&va) {
                    continue;
                }
                ambiguous_faces += 1;

                // How near the decision sat to its own boundary.
                let d02 = va[0] * va[2];
                let d13 = va[1] * va[3];
                let scale = d02.abs().max(d13.abs()).max(f64::MIN_POSITIVE);
                closest_margin = closest_margin.min((d02 - d13).abs() / scale);

                if crate::marching_cubes::ambiguity::face_is_joined(va)
                    != crate::marching_cubes::ambiguity::face_is_joined(vb)
                {
                    decision_flips += 1;
                }
            }
        }
    }

    std::println!(
        "measured: {divergent_planes} seam planes where the two expressions differ bit for bit, \
         {seam_faces} faces on them, {sign_disagreements} where the ulp moved a corner across zero, \
         {ambiguous_faces} ambiguous, {decision_flips} where the decider flipped; \
         closest relative margin {closest_margin:.3e}"
    );

    // The sweep has to actually reach the configuration, or it proves nothing —
    // the trap this repo has fallen into twice (M-32, M-38).
    assert!(
        ambiguous_faces > 100,
        "only {ambiguous_faces} ambiguous seam faces reached; the sweep is not measuring this"
    );
    assert_eq!(sign_disagreements, 0);
    assert_eq!(decision_flips, 0);
}
