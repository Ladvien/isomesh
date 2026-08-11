//! Tests for Hermite extraction.
//!
//! The load-bearing one is `normals_agree_on_a_face_and_disagree_across_an_edge`
//! — that is the ticket's acceptance criterion and the property the whole of
//! dual contouring rests on. If normals were constant across a cell, there would
//! be nothing for a vertex solve to solve.

use super::{HermiteCell, HermiteCrossing};
use crate::cube::{EDGE_CORNERS, corner_offset, is_inside};
use crate::fields::{BoxExact, Sphere};
use crate::{Sdf, vec3};

/// Sample the eight corners of a cell, as an extractor would before calling in.
fn corners<S: Sdf<Scalar = f64>>(sdf: &S, origin: [f64; 3], size: f64) -> [f64; 8] {
    let mut values = [0.0; 8];
    for (corner, slot) in values.iter_mut().enumerate() {
        let offset = corner_offset(corner as u8);
        *slot = sdf.sample([
            origin[0] + size * f64::from(offset[0]),
            origin[1] + size * f64::from(offset[1]),
            origin[2] + size * f64::from(offset[2]),
        ]);
    }
    values
}

fn cell<S: Sdf<Scalar = f64>>(sdf: &S, origin: [f64; 3], size: f64) -> HermiteCell<f64> {
    HermiteCell::from_corners(sdf, &corners(sdf, origin, size), origin, size)
}

/// **The acceptance criterion.**
///
/// On a cell straddling one flat face of a box, every crossing normal is the
/// same face normal. On a cell straddling an edge where two faces meet, they are
/// not — and the disagreement is what a vertex solve reads to find the edge.
#[test]
fn normals_agree_on_a_face_and_disagree_across_an_edge() {
    let field = BoxExact::<f64>::canonical(); // [-1, 1]^3
    let size = 0.25;

    // A cell straddling the +x face, well away from any edge of the box.
    let face = cell(&field, [0.9, -0.1, -0.1], size);
    assert!(
        face.len() >= 4,
        "expected a face crossing, got {}",
        face.len()
    );
    let first = face.iter().next().expect("a crossing").normal;
    for crossing in face.iter() {
        let agreement = vec3::dot(crossing.normal, first);
        assert!(
            (agreement - 1.0).abs() < 1e-9,
            "a flat face must give one normal, got {:?} vs {first:?}",
            crossing.normal
        );
    }
    // And that normal is the face's own.
    assert!((first[0] - 1.0).abs() < 1e-9, "{first:?}");

    // A cell straddling the +x/+y edge of the box.
    let edge = cell(&field, [0.9, 0.9, -0.1], size);
    assert!(
        edge.len() >= 4,
        "expected an edge crossing, got {}",
        edge.len()
    );
    let mut minimum_agreement = 1.0f64;
    for a in edge.iter() {
        for b in edge.iter() {
            minimum_agreement = minimum_agreement.min(vec3::dot(a.normal, b.normal));
        }
    }
    // Two perpendicular faces meet here, so some pair of normals is at 90°.
    assert!(
        minimum_agreement < 1e-9,
        "an edge cell must carry disagreeing normals, min dot was {minimum_agreement}"
    );
}

/// Normals are unit length, because the QEF is written in terms of `nᵢ · (x − pᵢ)`
/// and a non-unit normal silently reweights that plane against the others.
#[test]
fn normals_are_unit_length() {
    for field_origin in [[0.9, -0.1, -0.1], [0.9, 0.9, -0.1], [-1.1, -0.1, 0.3]] {
        let cell = cell(&BoxExact::<f64>::canonical(), field_origin, 0.25);
        for crossing in cell.iter() {
            let length = vec3::length(crossing.normal);
            assert!(
                (length - 1.0).abs() < 1e-12,
                "normal {:?} has length {length}",
                crossing.normal
            );
        }
    }
}

/// Every crossing sits on its own edge, between that edge's two corners.
#[test]
fn crossings_lie_on_the_edges_they_belong_to() {
    let field = Sphere::<f64>::canonical();
    let origin = [0.8, -0.2, -0.2];
    let size = 0.4;
    let cell = cell(&field, origin, size);
    assert!(!cell.is_empty());

    for edge in 0..12u8 {
        let Some(crossing) = cell.get(edge) else {
            continue;
        };
        let [lo, hi] = EDGE_CORNERS[edge as usize];
        let (lo_offset, hi_offset) = (corner_offset(lo), corner_offset(hi));
        for axis in 0..3 {
            let a = origin[axis] + size * f64::from(lo_offset[axis]);
            let b = origin[axis] + size * f64::from(hi_offset[axis]);
            let (low, high) = if a < b { (a, b) } else { (b, a) };
            assert!(
                crossing.position[axis] >= low - 1e-12 && crossing.position[axis] <= high + 1e-12,
                "edge {edge} axis {axis}: {} outside [{low}, {high}]",
                crossing.position[axis]
            );
        }
    }
}

/// The mask has to agree with the corner signs, or the solve is fed edges that
/// carry no surface.
#[test]
fn the_mask_matches_the_corner_signs() {
    let field = Sphere::<f64>::canonical();
    let origin = [0.8, -0.2, -0.2];
    let values = corners(&field, origin, 0.4);
    let cell = HermiteCell::from_corners(&field, &values, origin, 0.4);

    let mut expected = 0usize;
    for (edge, [lo, hi]) in EDGE_CORNERS.iter().copied().enumerate() {
        let cut = is_inside(values[lo as usize]) != is_inside(values[hi as usize]);
        assert_eq!(cell.contains(edge as u8), cut, "edge {edge}");
        expected += usize::from(cut);
    }
    assert_eq!(cell.len(), expected);
    assert_eq!(cell.iter().count(), expected);
}

/// A cell the surface misses carries nothing at all.
#[test]
fn a_cell_with_no_crossing_is_empty() {
    let field = Sphere::<f64>::canonical();
    let cell = cell(&field, [1.5, 1.5, 1.5], 0.25);
    assert!(cell.is_empty());
    assert_eq!(cell.len(), 0);
    assert_eq!(cell.iter().count(), 0);
    assert_eq!(cell.centroid(), None);
    assert!(cell.get(0).is_none());
}

/// The centroid is exactly the vertex Surface Nets would place, so dual
/// contouring starts from the same point and the comparison in E-104 isolates
/// the solve.
#[test]
fn the_centroid_is_the_surface_nets_vertex() {
    let field = BoxExact::<f64>::canonical();
    let origin = [0.9, 0.9, -0.1];
    let cell = cell(&field, origin, 0.25);
    let centroid = cell.centroid().expect("a crossing");

    let mut expected = [0.0f64; 3];
    for crossing in cell.iter() {
        for (axis, slot) in expected.iter_mut().enumerate() {
            *slot += crossing.position[axis];
        }
    }
    let inverse = 1.0 / cell.len() as f64;
    for axis in 0..3 {
        assert!((centroid[axis] - expected[axis] * inverse).abs() < 1e-15);
    }
}

/// Extraction reads only the corner values it is handed, so the same cell
/// sampled twice gives the same answer bit for bit.
#[test]
fn extraction_is_deterministic() {
    let field = BoxExact::<f64>::canonical();
    let origin = [0.9, 0.9, -0.1];
    let a = cell(&field, origin, 0.25);
    let b = cell(&field, origin, 0.25);
    assert_eq!(a.len(), b.len());
    for (left, right) in a.iter().zip(b.iter()) {
        assert_eq!(left, right);
    }
}

#[test]
fn f32_extraction_works_too() {
    let field = BoxExact::<f32>::canonical();
    let mut values = [0.0f32; 8];
    for (corner, slot) in values.iter_mut().enumerate() {
        let offset = corner_offset(corner as u8);
        *slot = field.sample([
            0.9 + 0.25 * offset[0] as f32,
            -0.1 + 0.25 * offset[1] as f32,
            -0.1 + 0.25 * offset[2] as f32,
        ]);
    }
    let cell = HermiteCell::from_corners(&field, &values, [0.9, -0.1, -0.1], 0.25);
    assert!(!cell.is_empty());
    for crossing in cell.iter() {
        assert!((vec3::length(crossing.normal) - 1.0).abs() < 1e-5);
    }
}

/// A sphere's normals all point away from its centre, which is the sanity check
/// that the gradient is being read the right way round.
#[test]
fn normals_point_away_from_the_solid() {
    let field = Sphere::<f64>::canonical();
    let cell = cell(&field, [0.8, -0.2, -0.2], 0.4);
    assert!(!cell.is_empty());
    for HermiteCrossing { position, normal } in cell.iter() {
        let outward = vec3::scale(*position, vec3::length(*position).recip());
        assert!(
            vec3::dot(*normal, outward) > 0.999,
            "normal {normal:?} should agree with the radial direction {outward:?}"
        );
    }
}
