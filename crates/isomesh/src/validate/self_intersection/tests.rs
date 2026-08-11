//! Fixtures for the self-intersection counter.
//!
//! The load-bearing ones are the negatives. A counter that finds crossings is
//! easy; a counter that does not report every adjacent pair in a well-formed
//! mesh is the whole difficulty, and it is what the closed-surface fixtures at
//! the bottom check.

// The rate assertions compare exact values: zero intersections must be exactly
// zero, and the per-thousand arithmetic on small integers is exact.
#![allow(clippy::float_cmp)]

use super::*;
use alloc::vec;

const H: f64 = 1.0;

/// Two triangles crossing transversely.
///
/// The vertical triangle meets the `z = 0` plane along the segment from
/// `(0.5, 0.5)` to `(1.5, 1.5)`. The first endpoint is inside the flat triangle
/// (`x + y = 1 < 2`) and the second is outside (`x + y = 3 > 2`), so the segment
/// crosses its boundary and the two genuinely pass through each other.
#[test]
fn crossing_triangles_are_detected() {
    let p = vec![
        [0.0, 0.0, 0.0],
        [2.0, 0.0, 0.0],
        [0.0, 2.0, 0.0],
        [0.5, 0.5, -1.0],
        [0.5, 0.5, 1.0],
        [1.5, 1.5, 0.0],
    ];
    let idx = vec![0, 1, 2, 3, 4, 5];
    let r = self_intersections(&p, &idx, H);

    assert_eq!(r.pairs, [[0, 1]]);
    assert_eq!(r.count(), 1);
    assert_eq!(r.triangles, 2);
    assert!((r.per_thousand_triangles() - 500.0).abs() < 1e-9);
    assert!(!r.is_intersection_free());
}

/// **The trap.** Two triangles sharing an edge touch along it by construction.
/// Reporting that would make every closed mesh look catastrophic.
#[test]
fn edge_adjacent_triangles_are_not_counted() {
    let p = vec![
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0],
    ];
    // Both use the edge {0, 1}, folded out of plane so this is not merely the
    // coplanar case in disguise.
    let idx = vec![0, 1, 2, 0, 1, 3];
    let r = self_intersections(&p, &idx, H);

    assert!(r.is_intersection_free(), "{r}");
    assert_eq!(r.adjacent_pairs_skipped, 1);
    assert_eq!(r.tested_pairs, 0);
}

/// The same trap one step further out: a vertex fan. Every pair in it touches at
/// the shared vertex, so a filter that only excluded edge-adjacency would report
/// a defect for each one.
#[test]
fn vertex_adjacent_triangles_are_not_counted() {
    let p = vec![
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [1.0, 1.0, 0.0],
        [-1.0, 0.0, 0.0],
        [-1.0, -1.0, 0.0],
    ];
    let idx = vec![0, 1, 2, 0, 3, 4];
    let r = self_intersections(&p, &idx, H);

    assert!(r.is_intersection_free(), "{r}");
    assert_eq!(r.adjacent_pairs_skipped, 1);
}

#[test]
fn distant_triangles_are_not_counted() {
    let p = vec![
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [50.0, 50.0, 50.0],
        [51.0, 50.0, 50.0],
        [50.0, 51.0, 50.0],
    ];
    let idx = vec![0, 1, 2, 3, 4, 5];
    let r = self_intersections(&p, &idx, H);

    assert!(r.is_intersection_free());
    // The grid separated them, so the exact test never ran.
    assert_eq!(r.tested_pairs, 0);
}

/// Two triangles in the same plane, overlapping. A degenerate but real
/// self-intersection, and the case the transverse interval test cannot see.
#[test]
fn coplanar_overlapping_triangles_are_detected() {
    let p = vec![
        [0.0, 0.0, 0.0],
        [2.0, 0.0, 0.0],
        [0.0, 2.0, 0.0],
        [0.5, 0.5, 0.0],
        [2.5, 0.5, 0.0],
        [0.5, 2.5, 0.0],
    ];
    let idx = vec![0, 1, 2, 3, 4, 5];
    let r = self_intersections(&p, &idx, H);

    assert_eq!(r.count(), 1, "{r}");
}

#[test]
fn coplanar_disjoint_triangles_are_not_counted() {
    let p = vec![
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [3.0, 3.0, 0.0],
        [4.0, 3.0, 0.0],
        [3.0, 4.0, 0.0],
    ];
    let idx = vec![0, 1, 2, 3, 4, 5];
    let r = self_intersections(&p, &idx, H);

    assert!(r.is_intersection_free(), "{r}");
}

/// Coplanar and edge-to-edge. They share a boundary segment and no area, which
/// is a tangential contact rather than an overlap.
#[test]
fn coplanar_touching_triangles_are_not_counted() {
    let p = vec![
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [1.0, 0.0, 0.0],
        [1.0, 1.0, 0.0],
        [0.0, 1.0, 0.0],
    ];
    // Distinct indices, so the adjacency filter does not hide the geometry.
    let idx = vec![0, 1, 2, 3, 4, 5];
    let r = self_intersections(&p, &idx, H);

    assert!(r.is_intersection_free(), "{r}");
    assert_eq!(r.adjacent_pairs_skipped, 0, "these share no index");
    assert_eq!(r.tested_pairs, 1, "the exact test really did run");
}

/// A closed, well-formed surface must report zero. Every pair in it is adjacent
/// or disjoint, so this is the end-to-end version of the two trap fixtures.
#[test]
fn a_tetrahedron_is_intersection_free() {
    let p = vec![
        [1.0, 1.0, 1.0],
        [1.0, -1.0, -1.0],
        [-1.0, 1.0, -1.0],
        [-1.0, -1.0, 1.0],
    ];
    let idx = vec![0, 1, 2, 0, 2, 3, 0, 3, 1, 1, 3, 2];
    let r = self_intersections(&p, &idx, H);

    assert!(r.is_intersection_free(), "{r}");
    assert_eq!(r.triangles, 4);
    // All six pairs share an edge in a tetrahedron.
    assert_eq!(r.adjacent_pairs_skipped, 6);
    assert_eq!(r.per_thousand_triangles(), 0.0);
}

/// A denser closed surface, where the grid actually has work to do and where
/// most pairs are neither adjacent nor trivially far apart.
#[test]
fn a_torus_grid_is_intersection_free() {
    let (p, idx) = torus_grid(12, 8);
    let r = self_intersections(&p, &idx, 0.25);

    assert!(r.is_intersection_free(), "{r}");
    assert_eq!(r.triangles, 192);
    assert!(
        r.tested_pairs > 0,
        "the grid must not have separated everything"
    );
}

fn torus_grid(m: u32, n: u32) -> (alloc::vec::Vec<[f64; 3]>, alloc::vec::Vec<u32>) {
    let tau = 2.0 * core::f64::consts::PI;
    let (major, minor) = (1.0f64, 0.3f64);
    let mut positions = alloc::vec::Vec::new();
    for i in 0..m {
        let theta = tau * f64::from(i) / f64::from(m);
        for j in 0..n {
            let phi = tau * f64::from(j) / f64::from(n);
            let ring = major + minor * Real::cos(phi);
            positions.push([
                ring * Real::cos(theta),
                minor * Real::sin(phi),
                ring * Real::sin(theta),
            ]);
        }
    }
    let at = |i: u32, j: u32| (i % m) * n + (j % n);
    let mut indices = alloc::vec::Vec::new();
    for i in 0..m {
        for j in 0..n {
            let (a, b, c, d) = (at(i, j), at(i + 1, j), at(i + 1, j + 1), at(i, j + 1));
            indices.extend_from_slice(&[a, b, c, a, c, d]);
        }
    }
    (positions, indices)
}

/// A zero-area triangle has no plane to intersect against. Counted, not
/// silently treated as tested.
#[test]
fn degenerate_triangles_are_excluded_and_reported() {
    let p = vec![
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [2.0, 0.0, 0.0], // collinear with the previous two
        [0.5, -1.0, 0.0],
        [0.5, 1.0, 0.0],
        [0.5, 0.0, 1.0],
    ];
    let idx = vec![0, 1, 2, 3, 4, 5];
    let r = self_intersections(&p, &idx, H);

    assert_eq!(r.degenerate_triangles, 1);
    assert_eq!(r.tested_pairs, 0);
    assert!(r.is_intersection_free());
}

/// The rate, not the fraction. `p = 1 − e^{−λT}` saturates with chunk size,
/// which is why this reports `λ`.
#[test]
fn the_rate_is_per_thousand_triangles() {
    let empty: alloc::vec::Vec<[f64; 3]> = alloc::vec::Vec::new();
    let r = self_intersections(&empty, &[], H);
    assert_eq!(r.triangles, 0);
    assert_eq!(r.per_thousand_triangles(), 0.0, "no divide by zero");
    assert!(r.is_intersection_free());

    let (p, idx) = torus_grid(12, 8);
    let r = self_intersections(&p, &idx, 0.25);
    assert_eq!(r.per_thousand_triangles(), 0.0);
}

#[test]
fn malformed_indices_do_not_panic() {
    let p = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0f64]];
    // Out of range, a repeated index, and a trailing partial triangle.
    let idx = vec![0, 1, 9, 0, 1, 1, 0, 1, 2, 0];
    let r = self_intersections(&p, &idx, H);
    assert_eq!(r.triangles, 1, "only the well-formed triangle survives");
}

#[test]
fn results_are_deterministic() {
    let (p, idx) = torus_grid(10, 6);
    let a = self_intersections(&p, &idx, 0.25);
    let b = self_intersections(&p, &idx, 0.25);
    assert_eq!(a, b);
    assert_eq!(alloc::format!("{a}"), alloc::format!("{b}"));
}

#[test]
fn pairs_come_back_sorted_so_they_can_be_bucketed_later() {
    // Three mutually crossing triangles: a triple of pairs, in a fixed order.
    let p = vec![
        [-1.0, -1.0, 0.0],
        [3.0, -1.0, 0.0],
        [-1.0, 3.0, 0.0],
        [0.5, -1.0, -1.0],
        [0.5, 3.0, -1.0],
        [0.5, 0.5, 1.0],
        [-1.0, 0.5, -1.0],
        [3.0, 0.5, -1.0],
        [0.5, 0.5, 1.5],
    ];
    let idx = vec![0, 1, 2, 3, 4, 5, 6, 7, 8];
    let r = self_intersections(&p, &idx, 2.0);

    assert!(r.count() >= 2, "{r}");
    let mut sorted = r.pairs.clone();
    sorted.sort_unstable();
    assert_eq!(r.pairs, sorted);
    for pair in &r.pairs {
        assert!(pair[0] < pair[1], "each pair is [low, high]");
    }
}

#[test]
#[should_panic(expected = "finite positive cell size")]
fn a_meaningless_cell_size_is_rejected() {
    let p = vec![[0.0, 0.0, 0.0f64]];
    let _ = self_intersections(&p, &[], 0.0);
}

/// The blow-up guard. A spacing far finer than the mesh would grow the grid
/// until it exhausted memory, so it says so instead.
#[test]
#[should_panic(expected = "does not describe this mesh")]
fn a_cell_size_that_does_not_describe_the_mesh_is_rejected() {
    let p = vec![[0.0, 0.0, 0.0], [1000.0, 0.0, 0.0], [0.0, 1000.0, 0.0f64]];
    let idx = vec![0, 1, 2];
    let _ = self_intersections(&p, &idx, 1e-3);
}
