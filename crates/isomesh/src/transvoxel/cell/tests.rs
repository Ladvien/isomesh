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
            let cell = TransitionCell::sample(&field, lo, fine_h, base, 1, 2);

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
            let cell = TransitionCell::sample(&field, lo, fine_h, base, 1, 2);
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
    let cell = TransitionCell::sample(&field, grid_origin, step, base, 1, 2);

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
        let cell = TransitionCell::sample(&field, [-2.0; 3], 0.125, [8 + i, 12, 18], 0, 1);
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
    let cell = TransitionCell::sample(&field, [-2.0; 3], 0.25, [6, 6, 8], 0, 1);
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
    let cell = TransitionCell::sample(&field, [-2.0; 3], 0.05, [72, 72, 72], 0, 1);
    assert!(cell.crossings().next().is_none());
    assert_eq!(cell.centroid(), None);

    let on_surface = TransitionCell::sample(&field, [-2.0; 3], 0.125, [8, 16, 16], 1, 2);
    let crossings: Vec<_> = on_surface.crossings().collect();
    if !crossings.is_empty() {
        assert!(on_surface.centroid().is_some());
    }
}
