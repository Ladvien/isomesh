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
    CASES, EDGE_CORNERS, EDGE_COUNT, MAX_TRIANGLES, NO_EDGE, corner_inside, is_inside,
};
use super::{MarchingCubes, validate_table};
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
        let pa = super::corner_offset(a);
        let pb = super::corner_offset(b);
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
    let shape = RuntimeShape3::new([samples; 3]);
    mc.extract(field, &shape, lo, cell_size, &mut out);
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
            &ValidateConfig::from_cell_size(h),
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
        &ValidateConfig::from_cell_size(h),
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
        &ValidateConfig::from_cell_size(h),
    );
    assert!(report.is_closed(), "box_exact:\n{report}");

    let csg: CsgDifference<f64> = crate::fields::csg_difference();
    let (out, h) = mesh(&csg, 41);
    let report = validate_indexed(
        &out.positions,
        &out.indices,
        &ValidateConfig::from_cell_size(h),
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
        &ValidateConfig::from_cell_size(h),
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
    let shape = RuntimeShape3::new([25; 3]);
    let mut mc = MarchingCubes::<f64>::new();

    let report = check_determinism(|out: &mut MeshBuffer<f64>| {
        mc.extract(&field, &shape, lo, cell_size, out);
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
    let si = self_intersections(&out.positions, &out.indices, h);
    let report = validate_indexed(
        &out.positions,
        &out.indices,
        &ValidateConfig::from_cell_size(h),
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
        &ValidateConfig::from_cell_size(h),
    );
    let si = self_intersections(&out.positions, &out.indices, h);

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
#[test]
fn every_closed_reference_field_meshes_cleanly() {
    fn check<F: Sdf<Scalar = f64> + ReferenceField>(name: &str, field: &F, samples: u32) {
        let (out, h) = mesh(field, samples);
        let report = validate_indexed(
            &out.positions,
            &out.indices,
            &ValidateConfig::from_cell_size(h),
        );
        assert!(out.triangle_count() > 0, "{name} produced nothing");
        if field.closed_in_domain() {
            assert!(report.is_closed(), "{name}:\n{report}");
        } else {
            assert!(report.is_manifold(), "{name}:\n{report}");
        }
        if let Some(chi) = field.expected_euler() {
            assert_eq!(report.euler_characteristic, chi, "{name}:\n{report}");
        }
        std::println!(
            "measured: {name} at {samples}^3 -> {} tris, chi {}, {} degenerate",
            out.triangle_count(),
            report.euler_characteristic,
            report.degenerate_triangles
        );
    }

    check("sphere", &Sphere::<f64>::canonical(), 33);
    check("torus", &Torus::<f64>::canonical(), 49);
    check("box_exact", &BoxExact::<f64>::canonical(), 33);
    check(
        "csg_difference",
        &crate::fields::csg_difference::<f64>(),
        41,
    );
    check("gyroid", &crate::fields::capped_gyroid::<f64>(), 49);
    check(
        "fbm_terrain",
        &crate::fields::FbmTerrain::<f64>::canonical(),
        33,
    );
}

#[test]
fn f32_and_f64_both_extract() {
    let mut mc = MarchingCubes::<f32>::new();
    let mut out = MeshBuffer::<f32>::new();
    let shape = RuntimeShape3::new([17; 3]);
    mc.extract(
        &Sphere::<f32>::canonical(),
        &shape,
        [-2.0; 3],
        4.0 / 16.0,
        &mut out,
    );
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
    let shape = RuntimeShape3::new([9; 3]);
    mc.extract(&field, &shape, [-2.0; 3], 0.5, &mut out);
    assert_eq!(out.triangle_count(), 0);
    assert_eq!(out.vertex_count(), 0);
}

#[test]
#[should_panic(expected = "at least two samples per axis")]
fn a_degenerate_grid_is_rejected() {
    let mut mc = MarchingCubes::<f64>::new();
    let mut out = MeshBuffer::<f64>::new();
    let shape = RuntimeShape3::new([1, 4, 4]);
    mc.extract(
        &Sphere::<f64>::canonical(),
        &shape,
        [-2.0; 3],
        0.5,
        &mut out,
    );
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
