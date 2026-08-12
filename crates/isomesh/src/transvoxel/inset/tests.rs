//! Equation 4.2, and the seam it has to leave closed.

use super::*;
use crate::fields::Sphere;
use crate::marching_cubes::MarchingCubes;
use crate::transvoxel::cell::TransitionCell;
use crate::{MeshBuffer, RuntimeShape3, Sdf};

/// The taper is linear across the first cell, zero in the middle, and linear
/// across the last — checked on hand-placed vertices rather than on a mesh, so a
/// sign error cannot hide in geometry.
#[test]
fn the_taper_is_linear_across_the_boundary_cells_and_zero_between() {
    // Exact: the taper is a multiply and an add.
    #![allow(clippy::float_cmp)]
    let mut mesh = MeshBuffer::<f64>::new();
    for c in [0.0, 0.5, 1.0, 2.0, 4.0, 7.0, 7.5, 8.0] {
        mesh.positions.push([c * 0.25, 0.0, 0.0]);
        mesh.normals.push([1.0, 0.0, 0.0]);
    }
    let width = 0.125;
    inset_boundary(&mut mesh, [0.0; 3], 8, 0.25, width, ALL_FACES).expect("valid");

    let expected = [
        0.0 + width,
        0.125 + 0.5 * width,
        0.25,
        0.5,
        1.0,
        1.75,
        1.875 - 0.5 * width,
        2.0 - width,
    ];
    for (got, want) in mesh.positions.iter().zip(expected) {
        assert_eq!(got[0], want, "positions {:?}", mesh.positions);
    }
}

/// A face with no coarser neighbour must not move, or the block pulls away from a
/// same-resolution neighbour and opens a seam where there was none.
#[test]
fn only_the_selected_faces_are_tapered() {
    #![allow(clippy::float_cmp)]
    let mut mesh = MeshBuffer::<f64>::new();
    mesh.positions.push([0.0, 0.0, 0.0]);
    mesh.normals.push([1.0, 0.0, 0.0]);
    inset_boundary(&mut mesh, [0.0; 3], 8, 0.25, 0.125, face_bit(0, 0)).expect("valid");
    assert_eq!(mesh.positions[0], [0.125, 0.0, 0.0], "only x should move");
}

/// Zero width is a no-op, which is the configuration A-011b shipped.
#[test]
fn a_zero_width_changes_nothing() {
    #![allow(clippy::float_cmp)]
    let mut mesh = MeshBuffer::<f64>::new();
    mesh.positions.push([0.0, 0.5, 1.0]);
    mesh.normals.push([1.0, 0.0, 0.0]);
    let before = mesh.positions.clone();
    inset_boundary(&mut mesh, [0.0; 3], 8, 0.25, 0.0, ALL_FACES).expect("valid");
    assert_eq!(mesh.positions, before);
}

/// Meaningless inputs are refused rather than absorbed.
#[test]
fn a_meaningless_taper_is_rejected() {
    let mut mesh = MeshBuffer::<f64>::new();
    for (cells, cell_size, width) in [
        (8u32, 0.0f64, 0.1f64),
        (8, -0.25, 0.1),
        (8, 0.25, -0.1),
        (8, 0.25, f64::NAN),
        (1, 0.25, 0.1),
    ] {
        assert!(
            inset_boundary(&mut mesh, [0.0; 3], cells, cell_size, width, ALL_FACES).is_err(),
            "cells {cells}, cell_size {cell_size}, width {width} was accepted"
        );
    }
}

/// **A-011c's acceptance.** The seam stays closed at a real width, and the patch
/// stops being edge-on.
#[test]
fn the_seam_stays_closed_at_a_real_width() {
    // Exact: a vertex is in a plane or it is not.
    #![allow(clippy::float_cmp)]
    let field = Sphere::<f64>::canonical();
    let fine_h = 0.125;
    let coarse_h = fine_h + fine_h;
    let width = fine_h;
    let grid_origin = [-2.0, -2.0, -2.0];
    let seam_x = 0.0;

    let fine_shape = RuntimeShape3::new([17, 33, 33]).expect("valid shape");
    let mut fine = MeshBuffer::<f64>::new();
    MarchingCubes::<f64>::new()
        .extract(&field, &fine_shape, grid_origin, fine_h, &mut fine)
        .expect("extraction");

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

    inset_boundary(
        &mut coarse,
        [seam_x, -2.0, -2.0],
        8,
        coarse_h,
        width,
        face_bit(0, 0),
    )
    .expect("valid taper");

    let mut stitched = MeshBuffer::<f64>::new();
    stitched.append(&fine);
    stitched.append(&coarse);
    let mut alignment = 0.0f64;
    for jz in 0..16i64 {
        for jy in 0..16i64 {
            let cell = TransitionCell::sample(
                &field,
                grid_origin,
                fine_h,
                [16, 2 * jy, 2 * jz],
                1,
                2,
                width,
            );
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
                alignment = alignment.max(cos.abs());
            }
            stitched.append(&patch);
        }
    }

    crate::weld::Welder::<f64>::new()
        .weld(&mut stitched, fine_h * 1e-6)
        .expect("weld");

    let cfg = crate::validate::ValidateConfig::from_cell_size(fine_h).expect("valid cell size");
    let (_report, features) =
        crate::validate::validate_features(&stitched.positions, &stitched.indices, &cfg);

    let in_plane = |plane: f64| {
        features
            .boundary_edges
            .iter()
            .filter(|[a, b]| {
                stitched.positions[*a as usize][0] == plane
                    && stitched.positions[*b as usize][0] == plane
            })
            .count()
    };
    let at_fine = in_plane(seam_x);
    let at_coarse = in_plane(seam_x + width);

    std::println!(
        "width {width}: {at_fine} gaps at the fine plane, {at_coarse} at the coarse plane; \
         best |cos| against the surface normal {alignment:.3}"
    );

    assert_eq!(at_fine, 0, "the fine side of the seam is open");
    assert_eq!(at_coarse, 0, "the coarse side of the seam is open");
    assert!(
        alignment > 0.1,
        "the patch is still perpendicular to the surface: {alignment:.3e}"
    );
}
