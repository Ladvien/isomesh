//! Tests for A-005.
//!
//! The validity tests **weld first**, because a blocky mesh is closed as a
//! surface and open as an index buffer — every quad carries its own four
//! vertices so the hard edges survive. A-013's welder is what reconciles those
//! two statements, and using it here is the check that it does.

use alloc::vec::Vec;

use super::GreedyQuads;
use crate::fields::{BoxExact, ReferenceField, Sphere};
use crate::marching_cubes::MarchingCubes;
use crate::validate::{ValidateConfig, check_determinism, validate_indexed};
use crate::weld::Welder;
use crate::{MeshBuffer, RuntimeShape3, Sdf, vec3};

fn mesh<F: Sdf<Scalar = f64> + ReferenceField>(field: &F, samples: u32) -> (MeshBuffer<f64>, f64) {
    let (lo, hi) = field.domain();
    let h = (hi[0] - lo[0]) / f64::from(samples - 1);
    let shape = RuntimeShape3::new([samples; 3]).expect("valid shape");
    let mut out = MeshBuffer::<f64>::new();
    GreedyQuads::<f64>::new()
        .extract(field, &shape, lo, h, &mut out)
        .expect("extraction");
    (out, h)
}

/// **Welding closes a greedy mesh only where the merge produced no T-junctions,
/// and a grid-aligned box is the only reference field where that holds.**
///
/// Split vertices are deliberate — a cube corner has three faces at three
/// normals and sharing there would average them away — so the raw index buffer
/// describes an open surface: a merged quad contributes five edges of which four
/// are unshared. [`crate::weld`] fixes that wherever quads meet corner to
/// corner.
///
/// It cannot fix a **T-junction**, and greedy merging manufactures them: a long
/// quad butts against several short ones, and the vertex where the short ones
/// meet simply does not exist on the long one's edge. There is nothing for the
/// weld to merge it *with*. Splitting the long edge would be a different
/// algorithm.
///
/// Measured on a sphere at 33³: `2568 → 848` vertices and boundary edges
/// `2568 → 768`, so the weld closes 70% of the boundary and the rest is
/// T-junctions. On `box_exact`, where every face merges to exactly one quad and
/// no T-junction can arise, it closes completely.
#[test]
fn welding_closes_a_greedy_mesh_only_where_the_merge_left_no_t_junctions() {
    let cfg_for = |h: f64| ValidateConfig::from_cell_size(h).expect("valid cell size");

    // The clean case: six quads, corner to corner, no T-junctions anywhere.
    let (mut boxed, h) = mesh(&BoxExact::<f64>::canonical(), 33);
    let quads = boxed.triangle_count() / 2;
    let before = validate_indexed(&boxed.positions, &boxed.indices, &cfg_for(h));
    assert_eq!(
        before.edges as usize,
        quads * 5,
        "five edges per merged quad"
    );
    assert_eq!(before.boundary_edges as usize, quads * 4);
    Welder::<f64>::new()
        .weld(&mut boxed, crate::weld::epsilon_for(h))
        .expect("valid epsilon");
    let after = validate_indexed(&boxed.positions, &boxed.indices, &cfg_for(h));
    std::println!(
        "measured: greedy box_exact at 33^3 -- boundary edges {} -> {}, chi {}, closed {}",
        before.boundary_edges,
        after.boundary_edges,
        after.euler_characteristic,
        after.is_closed()
    );
    assert_eq!(
        after.boundary_edges, 0,
        "a box has no T-junctions:\n{after}"
    );
    assert!(after.is_closed(), "{after}");
    assert_eq!(after.euler_characteristic, 2, "{after}");

    // The general case: merged runs of differing lengths, so T-junctions.
    let (mut sphere, h) = mesh(&Sphere::<f64>::canonical(), 33);
    let before = validate_indexed(&sphere.positions, &sphere.indices, &cfg_for(h));
    let report = Welder::<f64>::new()
        .weld(&mut sphere, crate::weld::epsilon_for(h))
        .expect("valid epsilon");
    let after = validate_indexed(&sphere.positions, &sphere.indices, &cfg_for(h));
    std::println!(
        "measured: greedy sphere at 33^3 -- {} verts welded to {}, boundary edges {} -> {} \
         ({:.0}% closed), chi {}",
        report.vertices_before,
        report.vertices_after,
        before.boundary_edges,
        after.boundary_edges,
        100.0 * (1.0 - after.boundary_edges as f64 / before.boundary_edges as f64),
        after.euler_characteristic
    );
    assert!(
        after.boundary_edges > 0,
        "greedy merging manufactures T-junctions and no weld removes them"
    );
    assert!(
        after.boundary_edges < before.boundary_edges / 2,
        "the weld should still close most of it: {} -> {}",
        before.boundary_edges,
        after.boundary_edges
    );
}

/// Signed volume, which is the only check that sees a global winding flip.
#[test]
fn a_meshed_box_has_positive_signed_volume() {
    let (out, _) = mesh(&BoxExact::<f64>::canonical(), 33);
    let mut total = 0.0;
    for t in out.indices.chunks_exact(3) {
        let a = out.positions[t[0] as usize];
        let b = out.positions[t[1] as usize];
        let c = out.positions[t[2] as usize];
        total += vec3::dot(a, vec3::cross(b, c));
    }
    assert!(total / 6.0 > 0.0, "inside out: {}", total / 6.0);
}

/// A cube aligned to the grid is the case greedy merging should flatten
/// completely: six faces, two triangles each, whatever the resolution.
#[test]
fn a_grid_aligned_box_merges_to_six_quads() {
    for samples in [17u32, 33, 65] {
        let (out, _) = mesh(&BoxExact::<f64>::canonical(), samples);
        std::println!(
            "measured: greedy box_exact at {samples}^3 -> {} triangles",
            out.triangle_count()
        );
        assert_eq!(
            out.triangle_count(),
            12,
            "{samples}^3 should merge to 6 quads"
        );
    }
}

#[test]
fn extraction_is_deterministic() {
    let field = Sphere::<f64>::canonical();
    let shape = RuntimeShape3::new([25; 3]).expect("valid shape");
    let report = check_determinism(|out: &mut MeshBuffer<f64>| {
        GreedyQuads::<f64>::new()
            .extract(&field, &shape, [-2.0; 3], 4.0 / 24.0, out)
            .expect("extraction");
    });
    assert!(report.is_deterministic(), "{report}");
}

/// **The claim: greedy merging gives `2.76×` fewer triangles than face culling.**
///
/// Tier R, from the UE5 benchmark in the v1 catalog (`6,690` against `18,492`,
/// which is `2.764×`). Face culling here is the same occupancy with the merge
/// switched off — one quad per visible face — so the comparison is exact rather
/// than against someone else's reimplementation.
///
/// **Predicted before running:** the published figure is one number for one
/// scene and will not reproduce as a constant. Merging pays for flat runs, so
/// `box_exact` should collapse almost entirely while a sphere — whose blocky
/// surface is a staircase of short runs — should gain much less. If that holds,
/// `2.76×` is a property of their scene rather than of the algorithm.
#[test]
fn the_saving_over_face_culling_is_measured() {
    let mut ratios: Vec<(&str, f64)> = Vec::new();
    crate::for_each_reference_field!(f64, |name, field| {
        let samples = 33u32;
        let (lo, hi) = field.domain();
        let h = (hi[0] - lo[0]) / f64::from(samples - 1);
        let shape = RuntimeShape3::new([samples; 3]).expect("valid shape");

        let mut merged = MeshBuffer::<f64>::new();
        GreedyQuads::<f64>::new()
            .extract(&field, &shape, lo, h, &mut merged)
            .expect("extraction");
        if merged.triangle_count() > 0 {
            let culled = face_culled_triangles(&field, &shape, lo, h);
            let ratio = culled as f64 / merged.triangle_count() as f64;
            ratios.push((name, ratio));
            std::println!(
                "measured: {name} at {samples}^3 -- face culled {culled} tris, greedy {} ({ratio:.2}x fewer)",
                merged.triangle_count()
            );
        }
    });

    let lo = ratios.iter().map(|(_, r)| *r).fold(f64::INFINITY, f64::min);
    let hi = ratios.iter().map(|(_, r)| *r).fold(0.0f64, f64::max);
    std::println!("measured: greedy saving ranges {lo:.2}x to {hi:.2}x over face culling");
    assert!(lo > 1.0, "merging must never cost triangles");
    // Pinned as a spread, because a single figure is exactly what the published
    // claim gets wrong.
    assert!(
        hi / lo > 3.0,
        "the spread is the finding: {lo:.2} to {hi:.2}"
    );
}

/// Face culling: one quad per visible cell face, no merging. The baseline the
/// published `2.76x` is measured against.
///
/// **This used to be a second implementation** — its own occupancy sampling, its
/// own neighbour test, its own face count — which agreed with the real extractor
/// on the day it was written and had nothing keeping it that way. A-005's own
/// finding is that the unmerged count is the denominator of every ratio in M-56,
/// so a drift here would move numbers this repo has published. E-106 needed the
/// unmerged *mesh* rather than a count, and `Merge::Off` now serves both.
fn face_culled_triangles<F: Sdf<Scalar = f64>>(
    field: &F,
    shape: &RuntimeShape3,
    origin: [f64; 3],
    h: f64,
) -> usize {
    let mut out = MeshBuffer::<f64>::new();
    let mut mesher = GreedyQuads::<f64>::new();
    mesher.set_merge(crate::greedy_quads::Merge::Off);
    mesher
        .extract(field, shape, origin, h, &mut out)
        .expect("extraction");
    out.triangle_count()
}

/// The blocky surface is a different surface; the counts say by how much.
/// Recorded beside Marching Cubes so the tradeoff table has both ends.
#[test]
fn the_triangle_count_against_marching_cubes_is_recorded() {
    crate::for_each_reference_field!(f64, |name, field| {
        let samples = 33u32;
        let (lo, hi) = field.domain();
        let h = (hi[0] - lo[0]) / f64::from(samples - 1);
        let shape = RuntimeShape3::new([samples; 3]).expect("valid shape");

        let mut greedy = MeshBuffer::<f64>::new();
        GreedyQuads::<f64>::new()
            .extract(&field, &shape, lo, h, &mut greedy)
            .expect("extraction");
        let mut mc = MeshBuffer::<f64>::new();
        MarchingCubes::<f64>::new()
            .extract(&field, &shape, lo, h, &mut mc)
            .expect("extraction");
        if mc.triangle_count() > 0 {
            std::println!(
                "measured: {name} at {samples}^3 -- greedy {} tris, marching cubes {} ({:.3}x)",
                greedy.triangle_count(),
                mc.triangle_count(),
                greedy.triangle_count() as f64 / mc.triangle_count() as f64
            );
        }
    });
}
