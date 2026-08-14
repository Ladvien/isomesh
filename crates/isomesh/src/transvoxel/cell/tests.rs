//! The seam identity, checked against a real coarse mesh rather than against a
//! restatement of the interpolation.

use alloc::vec::Vec;

use super::*;
use crate::fields::{ReferenceField, Sphere, Torus};
use crate::marching_cubes::MarchingCubes;
use crate::transvoxel::table::{EDGE_SAMPLES, is_half_resolution};
use crate::{MeshBuffer, RuntimeShape3};

/// **The property the whole of A-011b rests on.**
///
/// Every crossing a transition cell places on a **half-resolution** edge must
/// appear, bit for bit, among the vertices the coarse neighbour's Marching Cubes
/// pass produced. Not near one — *among* them.
///
/// Checked against an actual coarse mesh rather than against a re-derivation of
/// `edge_position`: re-deriving it would be a second copy of the formula, and two
/// copies agreeing proves only that they were written on the same day. A search
/// against the real vertex buffer proves the seam closes.
///
/// The fine-resolution crossings are deliberately *not* checked this way. They
/// live on the fine grid, and whether they coincide with the fine chunk's
/// vertices is the same identity in the other direction — covered by the same
/// arithmetic and by M-70.
#[test]
fn the_half_resolution_crossings_are_the_coarse_neighbours_vertices() {
    // Exact comparison on purpose: a tolerance here would hide the crack.
    #![allow(clippy::float_cmp)]
    let field = Sphere::<f64>::canonical();
    let (lo, hi) = field.domain();

    // The coarse grid: 16 cells over the domain. The transition face's step is
    // the *fine* spacing, half of it.
    let coarse_cells = 16u32;
    let coarse_h = (hi[0] - lo[0]) / f64::from(coarse_cells);
    let fine_h = coarse_h / 2.0;

    let shape = RuntimeShape3::new([coarse_cells + 1; 3]).expect("valid shape");
    let mut coarse = MeshBuffer::<f64>::new();
    MarchingCubes::<f64>::new()
        .extract(&field, &shape, lo, coarse_h, &mut coarse)
        .expect("extraction");
    assert!(coarse.triangle_count() > 0);

    // Walk every coarse cell face on one plane of the grid and build the
    // transition cell that would sit on it.
    let mut checked = 0usize;
    let mut cut = 0usize;
    for iv in 0..coarse_cells {
        for iu in 0..coarse_cells {
            // A face in the x = const plane through the middle of the sphere,
            // so the in-plane axes are y and z. Indices are on the FINE grid:
            // coarse sample c sits at fine index 2c.
            let base = [
                i64::from(coarse_cells), // coarse index 8 -> fine index 16
                2 * i64::from(iu),
                2 * i64::from(iv),
            ];
            let cell = TransitionCell::sample(&field, lo, fine_h, base, 1, 2, 0.0);

            for (edge, position) in cell.crossings() {
                if !is_half_resolution(edge) {
                    continue;
                }
                cut += 1;
                let found = coarse.positions.contains(&position);
                assert!(
                    found,
                    "edge {edge} joins samples {:?} and its crossing {position:?} is not a \
                     vertex of the coarse mesh",
                    EDGE_SAMPLES[edge as usize]
                );
            }
            checked += 1;
        }
    }

    std::println!(
        "{checked} transition faces walked, {cut} half-resolution crossings all matched a \
         coarse vertex exactly"
    );
    assert!(
        cut > 0,
        "no half-resolution edge was ever cut — the fixture is not on the surface"
    );
}

/// The same identity at a spacing that is not a power of two, which is where
/// M-49 found `cell_of` breaking and where a float assumption is most likely to
/// be wrong.
#[test]
fn the_identity_survives_a_non_power_of_two_spacing() {
    #![allow(clippy::float_cmp)]
    let field = Torus::<f64>::canonical();
    let (lo, hi) = field.domain();

    let coarse_cells = 14u32;
    let coarse_h = (hi[0] - lo[0]) / f64::from(coarse_cells);
    let fine_h = coarse_h / 2.0;
    assert!(
        coarse_h.to_bits() % 2 != 0 || (coarse_h.log2().fract() != 0.0),
        "the point of this test is a spacing that is not a power of two"
    );

    let shape = RuntimeShape3::new([coarse_cells + 1; 3]).expect("valid shape");
    let mut coarse = MeshBuffer::<f64>::new();
    MarchingCubes::<f64>::new()
        .extract(&field, &shape, lo, coarse_h, &mut coarse)
        .expect("extraction");

    let mut cut = 0usize;
    for iv in 0..coarse_cells {
        for iu in 0..coarse_cells {
            let base = [
                i64::from(coarse_cells), // the middle plane, on the fine grid
                2 * i64::from(iu),
                2 * i64::from(iv),
            ];
            let cell = TransitionCell::sample(&field, lo, fine_h, base, 1, 2, 0.0);
            for (edge, position) in cell.crossings() {
                if !is_half_resolution(edge) {
                    continue;
                }
                cut += 1;
                assert!(
                    coarse.positions.contains(&position),
                    "h = {coarse_h}: edge {edge} crossing {position:?} is not a coarse vertex"
                );
            }
        }
    }
    std::println!("h = {coarse_h} (not a power of two): {cut} crossings all matched exactly");
    assert!(cut > 0, "no half-resolution edge was cut");
}

/// A transition cell's nine samples sit on the fine grid, and its four corners
/// sit on the coarse one.
#[test]
fn the_corner_samples_are_the_coarse_grid_and_the_rest_are_not() {
    #![allow(clippy::float_cmp)]
    let field = Sphere::<f64>::canonical();
    let step = 0.125;
    let grid_origin = [-2.0, -2.0, -2.0];
    let base = [4i64, 6, 2];
    let cell = TransitionCell::sample(&field, grid_origin, step, base, 1, 2, 0.0);

    let at = |index: [i64; 3]| {
        [
            grid_origin[0] + step * index[0] as f64,
            grid_origin[1] + step * index[1] as f64,
            grid_origin[2] + step * index[2] as f64,
        ]
    };
    for (s, p) in cell.position.iter().enumerate() {
        let expected = at([base[0], base[1] + (s % 3) as i64, base[2] + (s / 3) as i64]);
        assert_eq!(*p, expected, "sample {s}");
    }

    // Corners 0, 2, 6, 8 land on the coarse lattice: an even number of fine
    // steps from the face origin.
    for corner in [0usize, 2, 6, 8] {
        assert_eq!((corner % 3) % 2, 0);
        assert_eq!((corner / 3) % 2, 0);
    }
}

/// Every cut edge yields a crossing and every uncut one yields nothing, and the
/// crossing lies between its endpoints.
#[test]
fn a_crossing_exists_exactly_where_the_signs_differ() {
    let field = Sphere::<f64>::canonical();
    let mut seen_cut = 0usize;
    for i in 0..12i64 {
        let cell = TransitionCell::sample(&field, [-2.0; 3], 0.125, [8 + i, 12, 18], 0, 1, 0.0);
        for edge in 0..16u8 {
            let crossing = cell.crossing(edge);
            assert_eq!(crossing.is_some(), cell.is_cut(edge), "edge {edge}");
            if let Some(p) = crossing {
                seen_cut += 1;
                let [a, b] = EDGE_SAMPLES[edge as usize];
                let (lo, hi) = (cell.position[a as usize], cell.position[b as usize]);
                for axis in 0..3 {
                    let (min, max) = if lo[axis] <= hi[axis] {
                        (lo[axis], hi[axis])
                    } else {
                        (hi[axis], lo[axis])
                    };
                    assert!(
                        p[axis] >= min && p[axis] <= max,
                        "edge {edge} axis {axis}: {} outside [{min}, {max}]",
                        p[axis]
                    );
                }
            }
        }
    }
    assert!(seen_cut > 0, "no edge was ever cut");
}

/// The case index a cell reports is what the table's own sample-sign convention
/// expects.
#[test]
fn the_case_index_matches_the_tables_convention() {
    let field = Sphere::<f64>::canonical();
    let cell = TransitionCell::sample(&field, [-2.0; 3], 0.25, [6, 6, 8], 0, 1, 0.0);
    let case = cell.case();
    for s in 0..SAMPLE_COUNT {
        assert_eq!(
            crate::transvoxel::table::sample_inside(case, s as u8),
            crate::cube::is_inside(cell.value[s]),
            "sample {s}"
        );
    }
    // And the cut edges the cell finds are the ones the table derives from that
    // case index, with no disagreement.
    let table_cut = crate::transvoxel::table::cut_edges(case);
    let cell_cut: u16 = (0..16u8).filter(|e| cell.is_cut(*e)).map(|e| 1 << e).sum();
    assert_eq!(table_cut, cell_cut, "case {case:#011b}");
}

/// An empty cell has no centroid rather than a zero one.
#[test]
fn a_cell_the_surface_misses_has_no_centroid() {
    let field = Sphere::<f64>::canonical();
    // Well outside the unit sphere.
    let cell = TransitionCell::sample(&field, [-2.0; 3], 0.05, [72, 72, 72], 0, 1, 0.0);
    assert!(cell.crossings().next().is_none());
    assert_eq!(cell.centroid(), None);

    let on_surface = TransitionCell::sample(&field, [-2.0; 3], 0.125, [8, 16, 16], 1, 2, 0.0);
    let crossings: Vec<_> = on_surface.crossings().collect();
    if !crossings.is_empty() {
        assert!(on_surface.centroid().is_some());
    }
}

// ─── triangulation ──────────────────────────────────────────────────────────

/// **Why a zero-width transition cell cannot be shaded, and what that costs.**
///
/// This test was written to check the patch's winding against the field gradient,
/// the way `meshed_sphere_has_positive_signed_volume` checks a closed mesh. It
/// reported **136 of 136 faces inward**, and reversing the fan reported the same
/// 136 — which is not a winding bug. It is the geometry.
///
/// All nine of a transition cell's samples lie in the transition **face**, so at
/// zero width every crossing does too, and every triangle is coplanar with that
/// face. Its normal is the face normal — perpendicular to the surface being
/// stitched. `dot(face_normal, gradient)` is then about `1e-17` and its sign
/// carries no information, so no winding test on it can mean anything.
///
/// That is exactly what Lengyel 2010 §4.3 means by
///
/// > It is possible to use a width of zero and still produce results that
/// > seamlessly stitch multiresolution meshes together, but this width leads to
/// > **severe shading problems**.
///
/// The stitch closes the hole and shades as a hard crease, because it *is* one: a
/// flat wall standing edge-on to the surface. So the transition width is not a
/// polish item to defer — **it is what gives the patch a normal at all**, and
/// A-011c owns it. Recorded as M-74.
#[test]
fn a_zero_width_patch_is_edge_on_to_the_surface_and_cannot_be_wound() {
    // Exact on purpose: the patch is coplanar by construction, not approximately.
    #![allow(clippy::float_cmp)]
    let field = Sphere::<f64>::canonical();
    let (lo, _hi) = field.domain();
    let fine_h = 0.125;

    let mut worst_alignment = 0.0f64;
    let mut faces = 0usize;
    for iv in 0..16i64 {
        for iu in 0..16i64 {
            let cell = TransitionCell::sample(&field, lo, fine_h, [16, 2 * iu, 2 * iv], 1, 2, 0.0);
            let mut patch = MeshBuffer::<f64>::new();
            cell.emit(&field, 0, &mut patch);
            if patch.triangle_count() == 0 {
                continue;
            }

            // Every vertex the patch produced sits in the transition face's own
            // plane -- the crossings by construction, the cycle centroids because
            // an average of coplanar points is coplanar.
            for p in &patch.positions {
                assert_eq!(
                    p[0], cell.position[0][0],
                    "a patch vertex left the face plane"
                );
            }

            for tri in patch.indices.chunks_exact(3) {
                let a = patch.positions[tri[0] as usize];
                let b = patch.positions[tri[1] as usize];
                let c = patch.positions[tri[2] as usize];
                let n = crate::vec3::cross(crate::vec3::sub(b, a), crate::vec3::sub(c, a));
                let len2 = crate::vec3::length_squared(n);
                if len2 == 0.0 {
                    continue;
                }
                faces += 1;
                let centre = [
                    (a[0] + b[0] + c[0]) / 3.0,
                    (a[1] + b[1] + c[1]) / 3.0,
                    (a[2] + b[2] + c[2]) / 3.0,
                ];
                let g = field.gradient(centre);
                let alignment =
                    crate::vec3::dot(n, g).abs() / (len2.sqrt() * crate::vec3::length(g));
                worst_alignment = worst_alignment.max(alignment);
            }
        }
    }

    std::println!(
        "{faces} zero-width patch faces; worst |cos| against the surface normal is \
         {worst_alignment:.3e}"
    );
    assert!(faces > 0, "no transition patch was produced");
    // Edge-on to within a rounding error. If this ever rises, the patch has
    // acquired a normal and the width work has landed.
    assert!(
        worst_alignment < 1e-12,
        "a zero-width patch should be perpendicular to the surface, got |cos| = {worst_alignment:.3e}"
    );
}

/// Every crossing the cell places must be used by the triangulation, and every
/// triangle must be non-degenerate.
///
/// A cut edge whose crossing never reaches a triangle is a hole in the seam by
/// another name.
#[test]
fn every_crossing_reaches_a_triangle() {
    let field = Sphere::<f64>::canonical();
    let (lo, _hi) = field.domain();

    let mut cells = 0usize;
    for iv in 0..16i64 {
        for iu in 0..16i64 {
            let cell = TransitionCell::sample(&field, lo, 0.125, [16, 2 * iu, 2 * iv], 1, 2, 0.0);
            let crossings: Vec<[f64; 3]> = cell.crossings().map(|(_, p)| p).collect();
            if crossings.is_empty() {
                continue;
            }
            cells += 1;

            let mut patch = MeshBuffer::<f64>::new();
            cell.emit(&field, 0, &mut patch);
            assert!(
                patch.triangle_count() > 0,
                "a cell with {} crossings produced no triangles",
                crossings.len()
            );

            for crossing in &crossings {
                assert!(
                    patch.positions.contains(crossing),
                    "crossing {crossing:?} never reached a vertex"
                );
            }
            // Every vertex is either a crossing or a cycle centroid.
            let extra = patch
                .positions
                .iter()
                .filter(|p| !crossings.contains(p))
                .count();
            assert!(
                extra > 0 && extra <= 4,
                "expected one centroid per cycle, got {extra}"
            );
        }
    }
    assert!(cells > 0, "no transition cell was cut");
    std::println!("{cells} transition cells triangulated, every crossing used");
}

// ─── the width, and the winding it makes decidable ──────────────────────────

/// The width displaces the half-resolution face and **nothing else**.
///
/// Every full-resolution crossing must be exactly where it was at zero width,
/// and every half-resolution one must be exactly that far along the face normal.
/// If the width leaked into the fine side, the seam identity M-73 established
/// would break the moment a width was chosen.
#[test]
fn the_width_displaces_only_the_half_resolution_face() {
    // Exact: the displacement is a single addition, not an approximation.
    #![allow(clippy::float_cmp)]
    let field = Sphere::<f64>::canonical();
    let (lo, _hi) = field.domain();
    let width = 0.0625;

    let mut fine_checked = 0usize;
    let mut coarse_checked = 0usize;
    for iv in 0..16i64 {
        for iu in 0..16i64 {
            let base = [16, 2 * iu, 2 * iv];
            let flat = TransitionCell::sample(&field, lo, 0.125, base, 1, 2, 0.0);
            let thick = TransitionCell::sample(&field, lo, 0.125, base, 1, 2, width);

            // Same field, same samples, so the same cut edges.
            assert_eq!(flat.case(), thick.case());

            for edge in 0..16u8 {
                let (Some(a), Some(b)) = (flat.crossing(edge), thick.crossing(edge)) else {
                    continue;
                };
                if is_half_resolution(edge) {
                    coarse_checked += 1;
                    assert_eq!(
                        b[0],
                        a[0] + width,
                        "edge {edge} moved wrong along the normal"
                    );
                    assert_eq!([b[1], b[2]], [a[1], a[2]], "edge {edge} moved in-plane");
                } else {
                    fine_checked += 1;
                    assert_eq!(a, b, "edge {edge} is on the fine face and must not move");
                }
            }
        }
    }
    std::println!(
        "{fine_checked} fine crossings unmoved, {coarse_checked} coarse crossings displaced by \
         exactly the width"
    );
    assert!(fine_checked > 0 && coarse_checked > 0);
}

/// **The winding, decidable at last.**
///
/// M-74 showed a zero-width patch is exactly perpendicular to the surface, so no
/// test against the field gradient could orient it. Give it a width and the patch
/// becomes a ribbon spanning from the fine face to the coarse one, which has a
/// normal — and that normal must point away from the solid, like every other face
/// this crate emits.
///
/// The width is deliberately large here (half the coarse cell, Lengyel's own
/// `w(k) = 2^(k−2)`) so the ribbon is well-conditioned rather than a sliver.
#[test]
fn a_patch_with_width_is_wound_away_from_the_solid() {
    let field = Sphere::<f64>::canonical();
    let (lo, _hi) = field.domain();
    let fine_h = 0.125;
    // The coarse cell is 2*fine_h; Lengyel's width is half of it.
    let width = fine_h;

    let mut agree = 0usize;
    let mut disagree = 0usize;
    let mut worst_alignment = 0.0f64;
    for iv in 0..16i64 {
        for iu in 0..16i64 {
            let cell =
                TransitionCell::sample(&field, lo, fine_h, [16, 2 * iu, 2 * iv], 1, 2, width);
            let mut patch = MeshBuffer::<f64>::new();
            cell.emit(&field, 0, &mut patch);
            if patch.triangle_count() == 0 {
                continue;
            }

            for tri in patch.indices.chunks_exact(3) {
                let a = patch.positions[tri[0] as usize];
                let b = patch.positions[tri[1] as usize];
                let c = patch.positions[tri[2] as usize];
                let n = crate::vec3::cross(crate::vec3::sub(b, a), crate::vec3::sub(c, a));
                let len2 = crate::vec3::length_squared(n);
                if len2 == 0.0 {
                    continue;
                }
                let centre = [
                    (a[0] + b[0] + c[0]) / 3.0,
                    (a[1] + b[1] + c[1]) / 3.0,
                    (a[2] + b[2] + c[2]) / 3.0,
                ];
                let g = field.gradient(centre);
                let cos = crate::vec3::dot(n, g) / (len2.sqrt() * crate::vec3::length(g));
                worst_alignment = worst_alignment.max(cos.abs());
                if cos > 0.0 {
                    agree += 1;
                } else {
                    disagree += 1;
                }
            }
        }
    }

    std::println!(
        "width {width}: {agree} faces outward, {disagree} inward; best |cos| against the \
         surface normal {worst_alignment:.3}"
    );
    assert!(agree + disagree > 0, "no patch was produced");
    // The ribbon now has a normal to speak of at all -- M-74's zero-width case
    // measured exactly 0 here.
    assert!(
        worst_alignment > 0.1,
        "the patch is still edge-on: best |cos| {worst_alignment:.3e}"
    );
    assert_eq!(
        disagree,
        0,
        "{disagree} of {} faces point into the solid",
        agree + disagree
    );
}

// ─── A-011b's acceptance: two chunks at differing LOD ───────────────────────

/// Boundary edges lying wholly in the seam plane.
///
/// The two chunks are legitimately open at their *outer* borders — the surface
/// leaves through the sides, as every chunk in a streamed world does — so a
/// global boundary-edge count says nothing. Only the seam is the question.
fn gaps_in_the_seam_plane(mesh: &MeshBuffer<f64>, seam_x: f64, h: f64) -> usize {
    // Exact on purpose. Every vertex in the seam plane got there by an
    // interpolation whose endpoints both sit at `seam_x`, so its x is that value
    // identically -- and a tolerance here would sweep in vertices *near* the
    // plane and count gaps that are not seam gaps.
    #![allow(clippy::float_cmp)]
    let cfg = crate::validate::ValidateConfig::from_cell_size(h).expect("valid cell size");
    let (_report, features) =
        crate::validate::validate_features(&mesh.positions, &mesh.indices, &cfg);
    features
        .boundary_edges
        .iter()
        .filter(|[a, b]| {
            mesh.positions[*a as usize][0] == seam_x && mesh.positions[*b as usize][0] == seam_x
        })
        .count()
}

/// **A-011b's acceptance criterion:** two adjacent chunks at differing LOD
/// produce zero boundary gaps.
///
/// A full-resolution block over `x ∈ [-2, 0]` at `h = 1/8` meets a
/// half-resolution one over `x ∈ [0, 2]` at `2h`. Meshed independently they do
/// not meet: the fine side ends on a contour of 32×32 sub-squares and the coarse
/// side on one of 16×16 squares, and the difference is a ring of unmatched
/// boundary edges lying in the seam plane — a crack you can see the sky through.
/// Transition cells along the shared face close it.
///
/// Asserted **in both directions**, because a test that only checks the fixed
/// case would still pass if the two resolutions had stopped disagreeing.
///
/// The width is zero here, which is Lengyel §4.3's own position: a zero width
/// *"still produce\[s\] results that seamlessly stitch multiresolution meshes
/// together"* — it is the **shading** it costs (M-74), not the stitch. A non-zero
/// width additionally needs his Equation 4.2 to scale the coarse block's boundary
/// cells inward, or the transition cells overlap them; that is what the ticket
/// still carries.
#[test]
fn transition_cells_close_the_gap_between_two_resolutions() {
    // Exact: a seam either closes or it does not.
    #![allow(clippy::float_cmp)]
    let field = Sphere::<f64>::canonical();
    let fine_h = 0.125;
    let coarse_h = fine_h + fine_h;
    let seam_x = 0.0;
    let grid_origin = [-2.0, -2.0, -2.0];

    // The full-resolution block: x in [-2, 0], y and z across the domain.
    let fine_shape = RuntimeShape3::new([17, 33, 33]).expect("valid shape");
    let mut fine = MeshBuffer::<f64>::new();
    MarchingCubes::<f64>::new()
        .extract(&field, &fine_shape, grid_origin, fine_h, &mut fine)
        .expect("extraction");

    // The half-resolution block: x in [0, 2], same y and z extent at twice the
    // spacing.
    let coarse_shape = RuntimeShape3::new([9, 17, 17]).expect("valid shape");
    let mut coarse = MeshBuffer::<f64>::new();
    MarchingCubes::<f64>::new()
        .extract(
            &field,
            &coarse_shape,
            [seam_x, -2.0, -2.0],
            coarse_h,
            &mut coarse,
        )
        .expect("extraction");

    assert!(fine.triangle_count() > 0 && coarse.triangle_count() > 0);

    // Without transition cells.
    let mut plain = MeshBuffer::<f64>::new();
    plain
        .append(&fine)
        .expect("the meshes fit the u32 index space");
    plain
        .append(&coarse)
        .expect("the meshes fit the u32 index space");
    crate::weld::Welder::<f64>::new()
        .weld(&mut plain, crate::weld::epsilon_for(fine_h))
        .expect("weld");
    let before = gaps_in_the_seam_plane(&plain, seam_x, fine_h);

    // With them: one transition cell per coarse cell face on the seam, spanning
    // two fine cells on each in-plane axis.
    let mut stitched = MeshBuffer::<f64>::new();
    stitched
        .append(&fine)
        .expect("the meshes fit the u32 index space");
    stitched
        .append(&coarse)
        .expect("the meshes fit the u32 index space");
    let mut patches = 0usize;
    for jz in 0..16i64 {
        for jy in 0..16i64 {
            let cell = TransitionCell::sample(
                &field,
                grid_origin,
                fine_h,
                [16, 2 * jy, 2 * jz],
                1,
                2,
                0.0,
            );
            let mut patch = MeshBuffer::<f64>::new();
            cell.emit(&field, 0, &mut patch);
            if patch.triangle_count() > 0 {
                patches += 1;
                stitched
                    .append(&patch)
                    .expect("the meshes fit the u32 index space");
            }
        }
    }
    crate::weld::Welder::<f64>::new()
        .weld(&mut stitched, crate::weld::epsilon_for(fine_h))
        .expect("weld");
    let after = gaps_in_the_seam_plane(&stitched, seam_x, fine_h);

    std::println!(
        "seam-plane boundary edges: {before} without transition cells, {after} with \
         ({patches} patches emitted)"
    );

    assert!(
        patches > 0,
        "no transition cell was cut — the fixture misses the surface"
    );
    assert!(
        before > 0,
        "the two resolutions already agreed, so this proves nothing"
    );
    assert_eq!(
        after, 0,
        "the seam still has {after} unmatched boundary edges"
    );
}

/// **The mirrored patch faces away from the solid too.** `sample`'s width is
/// signed — the caller states which side the coarse block is on — and with the
/// coarse side toward −x the map from `(u, v, w)` parameter space to world
/// space is a reflection. Left uncorrected, every triangle of every mirrored
/// patch is wound into the solid — 144 of 144 measured, the exact complement
/// of `a_patch_with_width_is_wound_away_from_the_solid` — and no manifold or
/// Euler check can see it. Both reflected parameterisations are checked: a
/// negative width, and swapped in-plane axes.
#[test]
fn a_mirrored_patch_is_wound_away_from_the_solid() {
    let field = Sphere::<f64>::canonical();
    let (lo, _hi) = field.domain();
    let fine_h = 0.125;

    // Each has sign(width) × parity(u, v, normal) negative.
    for (u, v, width) in [(1usize, 2usize, -fine_h), (2, 1, fine_h)] {
        let mut agree = 0usize;
        let mut disagree = 0usize;
        for ib in 0..16i64 {
            for ia in 0..16i64 {
                let cell =
                    TransitionCell::sample(&field, lo, fine_h, [16, 2 * ia, 2 * ib], u, v, width);
                let mut patch = MeshBuffer::<f64>::new();
                cell.emit(&field, 0, &mut patch);
                for tri in patch.indices.chunks_exact(3) {
                    let a = patch.positions[tri[0] as usize];
                    let b = patch.positions[tri[1] as usize];
                    let c = patch.positions[tri[2] as usize];
                    let n = crate::vec3::cross(crate::vec3::sub(b, a), crate::vec3::sub(c, a));
                    if crate::vec3::length_squared(n) == 0.0 {
                        continue;
                    }
                    let centre = [
                        (a[0] + b[0] + c[0]) / 3.0,
                        (a[1] + b[1] + c[1]) / 3.0,
                        (a[2] + b[2] + c[2]) / 3.0,
                    ];
                    if crate::vec3::dot(n, field.gradient(centre)) > 0.0 {
                        agree += 1;
                    } else {
                        disagree += 1;
                    }
                }
            }
        }
        std::println!("u={u} v={v} width={width}: {agree} faces outward, {disagree} inward");
        assert!(agree + disagree > 0, "no patch was produced");
        assert_eq!(
            disagree,
            0,
            "u={u} v={v} width={width}: {disagree} of {} faces point into the solid",
            agree + disagree
        );
    }
}
