//! §3.2's tests.
//!
//! The two worth reading are
//! `the_appendix_formulas_predict_the_curves_section_3_1_actually_finds`, which
//! checks Theorems B.4 and B.6 as *equalities* against the curve reconstruction
//! rather than trusting either in isolation, and
//! `parallel_quads_all_split_the_same_way`, which is the property "the same,
//! arbitrary diagonal" actually asserts.

use super::*;
use crate::marching_tetrahedra::table::TET_EDGES;

/// The reference tetrahedron: the origin and the three unit axes.
const CORNERS: [[f64; 3]; 4] = [
    [0.0, 0.0, 0.0],
    [1.0, 0.0, 0.0],
    [0.0, 1.0, 0.0],
    [0.0, 0.0, 1.0],
];

/// Crossing parameters for an edge carrying `n` of them, evenly spaced and
/// strictly interior, so no two coincide and none lands on a corner.
fn spaced(n: u32) -> Vec<f64> {
    (0..n)
        .map(|i| (f64::from(i) + 1.0) / (f64::from(n) + 1.0))
        .collect()
}

/// A [`TetCrossings`] realising the given edge coordinates on [`CORNERS`].
fn crossings(count: [u32; TET_EDGE_COUNT]) -> Vec<Vec<f64>> {
    count.iter().map(|n| spaced(*n)).collect()
}

fn tet<'a>(owned: &'a [Vec<f64>]) -> TetCrossings<'a, f64> {
    let mut along: [&[f64]; TET_EDGE_COUNT] = [&[]; TET_EDGE_COUNT];
    for (slot, v) in along.iter_mut().zip(owned.iter()) {
        *slot = v.as_slice();
    }
    TetCrossings {
        corners: CORNERS,
        along,
    }
}

/// `(d₁, d₂)` on the three complementary pairs, which is Property II's form.
fn pattern_coords(d1: u32, d2: u32) -> EdgeCoordinates {
    let mut count = [0u32; TET_EDGE_COUNT];
    // Pairs are `e` and `5 - e`: (0,5), (1,4), (2,3).
    count[0] = d1;
    count[5] = d1;
    count[1] = d2;
    count[4] = d2;
    count[2] = d1 + d2;
    count[3] = d1 + d2;
    EdgeCoordinates { count }
}

#[test]
fn edge_coordinates_come_from_the_crossing_lists_not_from_a_second_argument() {
    let owned = crossings([2, 0, 3, 1, 0, 4]);
    let t = tet(&owned);
    assert_eq!(t.coordinates().count, [2, 0, 3, 1, 0, 4]);
}

#[test]
fn a_crossing_sits_where_its_parameter_puts_it_measured_from_the_lower_corner() {
    let owned = crossings([1, 0, 0, 0, 0, 0]);
    let t = tet(&owned);
    // Edge 0 joins corners 0 and 1; one crossing lands at t = 1/2.
    let [lo, hi] = TET_EDGES[0];
    assert_eq!([lo, hi], [0, 1]);
    let p = t.position(FacePoint { edge: 0, index: 0 });
    assert_eq!(p, Some([0.5, 0.0, 0.0]));
}

#[test]
fn unsorted_or_off_edge_crossings_are_rejected_rather_than_silently_repaired() {
    let owned = [
        alloc::vec![0.7, 0.3],
        Vec::new(),
        Vec::new(),
        Vec::new(),
        Vec::new(),
        Vec::new(),
    ];
    assert_eq!(
        tet(&owned).check(),
        Err(NotFillable::Unsorted { edge: 0, at: 1 })
    );

    let owned = [
        alloc::vec![1.5],
        Vec::new(),
        Vec::new(),
        Vec::new(),
        Vec::new(),
        Vec::new(),
    ];
    assert_eq!(
        tet(&owned).check(),
        Err(NotFillable::OffEdge { edge: 0, at: 0 })
    );
}

#[test]
fn a_classic_corner_cut_is_one_triangle() {
    // Corner 0's three incident edges are 01, 02, 03 -- indices 0, 1, 2.
    let coords = EdgeCoordinates::new([1, 1, 1, 0, 0, 0]);
    let cycles = cycles(&coords);
    assert_eq!(cycles.len(), 1);
    assert_eq!(cycles[0].length(), 3);
    assert!(cycles[0].is_corner_cut());

    let owned = crossings(coords.count);
    let mut patch = TetPatch::new();
    assert_eq!(fill(&tet(&owned), &mut patch), Ok(Unfilled::None));
    assert_eq!(patch.triangles.len(), 1);
    assert_eq!(patch.positions.len(), 3);
}

#[test]
fn a_classic_diagonal_cut_is_a_quad_and_becomes_two_triangles() {
    // Separating {0, 1} from {2, 3} cuts edges 02, 03, 12, 13 -- indices 1..4.
    let coords = EdgeCoordinates::new([0, 1, 1, 1, 1, 0]);
    let cycles = cycles(&coords);
    assert_eq!(cycles.len(), 1);
    assert_eq!(cycles[0].length(), 4);
    assert!(!cycles[0].is_corner_cut());

    let owned = crossings(coords.count);
    let mut patch = TetPatch::new();
    assert_eq!(fill(&tet(&owned), &mut patch), Ok(Unfilled::None));
    assert_eq!(patch.triangles.len(), 2);
    assert_eq!(patch.positions.len(), 4);
}

#[test]
fn every_classic_configuration_fills_completely() {
    // The 0/1 corner of the encoding is classic Marching Tetrahedra, and every
    // one of its cases is a corner cut or a quad -- so §3.2's first two cases
    // must already cover all sixteen sign patterns with nothing left over.
    for signs in 0u8..16 {
        let mut count = [0u32; TET_EDGE_COUNT];
        for (e, slot) in count.iter_mut().enumerate() {
            let [lo, hi] = TET_EDGES[e];
            let inside = |c: u8| signs & (1 << c) != 0;
            *slot = u32::from(inside(lo) != inside(hi));
        }
        let coords = EdgeCoordinates::new(count);
        let owned = crossings(count);
        let mut patch = TetPatch::new();
        assert_eq!(
            fill(&tet(&owned), &mut patch),
            Ok(Unfilled::None),
            "sign pattern {signs:04b} left something unfilled"
        );
        // 0 triangles for the two uniform patterns, 1 for a corner cut, 2 for a
        // diagonal -- which is exactly A-003's own table.
        let expected = match coords.count.iter().sum::<u32>() {
            0 => 0,
            3 => 1,
            4 => 2,
            n => panic!("classic configuration with {n} crossings"),
        };
        assert_eq!(patch.triangles.len(), expected, "sign pattern {signs:04b}");
    }
}

#[test]
fn the_appendix_formulas_predict_the_curves_section_3_1_actually_finds() {
    // Theorem B.4: the number of components is gcd(d₁, d₂).
    // Corollary B.6: every component has length 4(d₁ + d₂) / gcd(d₁, d₂).
    //
    // Both are checked against what §3.1's reconstruction independently
    // produces, so agreement is evidence about the implementation and not just
    // about the arithmetic.
    let mut checked = 0;
    for d1 in 1..=6u32 {
        for d2 in 0..=d1 {
            let coords = pattern_coords(d1, d2);
            assert!(
                coords.is_normal(),
                "({d1}, {d2}) should satisfy both normality conditions"
            );

            let pattern = Pattern::of(&coords).expect("Property II holds by construction");
            assert_eq!((pattern.d1, pattern.d2), (d1, d2));

            let cycles = cycles(&coords);
            assert!(
                cycles.iter().all(|c| c.kind == CurveKind::Normal),
                "({d1}, {d2}) produced a non-normal curve"
            );

            let predicted_count = pattern.loop_count();
            let predicted_length = pattern
                .loop_length()
                .expect("a non-empty pattern has a length");

            assert_eq!(
                cycles.len() as u32,
                predicted_count,
                "({d1}, {d2}): Theorem B.4 predicted {predicted_count} loops"
            );
            for c in &cycles {
                assert_eq!(
                    c.length() as u32,
                    predicted_length,
                    "({d1}, {d2}): Corollary B.6 predicted length {predicted_length}"
                );
            }
            checked += 1;
        }
    }
    // A zero here would mean the sweep never ran -- M-44's rule.
    assert_eq!(checked, 27);
}

#[test]
fn property_two_is_checked_and_not_assumed() {
    // A residual that is not of the form (d₁, d₁), (d₂, d₂), (d₁ + d₂, d₁ + d₂)
    // must be refused rather than forced into a pattern.
    assert_eq!(Pattern::of(&EdgeCoordinates::new([1, 0, 0, 0, 0, 2])), None);
    assert_eq!(Pattern::of(&EdgeCoordinates::new([3, 1, 1, 1, 1, 3])), None);
    assert_eq!(
        Pattern::of(&EdgeCoordinates::new([0, 0, 0, 0, 0, 0])),
        Some(Pattern { d1: 0, d2: 0 })
    );
}

#[test]
fn removing_the_corner_cuts_is_what_makes_property_two_hold() {
    // A corner cut on top of a quad pattern: Property II fails on the raw
    // coordinates and holds once the cut is subtracted. This is why §3.2.1
    // emits the triangles first.
    let mut count = pattern_coords(1, 0).count;
    for e in [0usize, 1, 2] {
        count[e] += 1;
    }
    let coords = EdgeCoordinates::new(count);
    assert_eq!(Pattern::of(&coords), None);

    let cycles = cycles(&coords);
    let residual = residual(&coords, &cycles).expect("cuts came from these coordinates");
    assert_eq!(Pattern::of(&residual), Some(Pattern { d1: 1, d2: 0 }));
}

#[test]
fn parallel_quads_all_split_the_same_way() {
    // "we split all quads along the same, arbitrary diagonal". Arbitrary is
    // real -- either diagonal is admissible -- but *the same* is not, and this
    // is what it means: in a family of m parallel quads, every quad's diagonal
    // joins the same two edges of the tet, so no two diagonals can cross.
    let coords = pattern_coords(3, 0);
    let cycles = cycles(&coords);
    assert_eq!(cycles.len(), 3, "gcd(3, 0) = 3 parallel quads");

    let diagonal_edges: Vec<[u8; 2]> = cycles
        .iter()
        .map(|c| {
            assert_eq!(c.length(), 4);
            // fill() splits 0-2, so the diagonal joins points 0 and 2.
            let mut e = [c.points[0].edge, c.points[2].edge];
            e.sort_unstable();
            e
        })
        .collect();

    assert!(
        diagonal_edges.windows(2).all(|w| w[0] == w[1]),
        "quads split along different edge pairs: {diagonal_edges:?}"
    );
}

#[test]
fn a_single_long_loop_fans_around_its_centre_of_mass() {
    // (3, 2): gcd = 1, so one loop of length 4(5)/1 = 20. Case (2) -- one
    // Steiner point at the centre of mass, one triangle per loop edge.
    for (d1, d2, length) in [(3u32, 2u32, 20usize), (3, 1, 16), (5, 3, 32), (4, 3, 28)] {
        let coords = pattern_coords(d1, d2);
        let cycles = cycles(&coords);
        assert_eq!(cycles.len(), 1, "({d1}, {d2}) should be a single loop");
        assert_eq!(cycles[0].length(), length, "({d1}, {d2})");

        let owned = crossings(coords.count);
        let mut patch = TetPatch::new();
        assert_eq!(
            fill(&tet(&owned), &mut patch),
            Ok(Unfilled::None),
            "({d1}, {d2}) should be the single-loop case"
        );
        assert_eq!(patch.triangles.len(), length, "({d1}, {d2})");
        assert_eq!(patch.positions.len(), length + 1, "({d1}, {d2})");

        // The Steiner point is the centre of mass of the loop's vertices, which
        // is what puts it inside their convex hull.
        let steiner = patch.positions[length];
        let mut sum = [0.0f64; 3];
        for p in &patch.positions[..length] {
            for (s, v) in sum.iter_mut().zip(p.iter()) {
                *s += v;
            }
        }
        for axis in 0..3 {
            let mean = sum[axis] / length as f64;
            assert!(
                (steiner[axis] - mean).abs() < 1e-12,
                "({d1}, {d2}): Steiner point is not the centre of mass"
            );
        }
    }
}

#[test]
fn every_implemented_case_emits_an_intersection_free_patch() {
    // The headline claim of §3.2 is that the output is intersection-free, so
    // this asserts **zero** rather than recording a metric -- the opposite of
    // the dual methods, where a non-zero count is a property of the algorithm.
    //
    // It is also the only test that can see a wrong Steiner assignment: a
    // mis-assigned octagon still fans correctly, still indexes valid vertices
    // and still has the right triangle count. It just passes through its
    // neighbours.
    let mut checked = 0;
    for d1 in 1..=5u32 {
        for d2 in 0..=d1 {
            let coords = pattern_coords(d1, d2);
            let owned = crossings(coords.count);
            let mut patch = TetPatch::new();
            let outcome = fill(&tet(&owned), &mut patch).expect("well-formed crossings");
            if outcome != Unfilled::None {
                continue;
            }

            let indices: Vec<u32> = patch.triangles.iter().flatten().copied().collect();
            let report = crate::validate::self_intersections(&patch.positions, &indices, 1.0)
                .expect("a unit tet is a valid spacing");
            assert_eq!(
                report.count(),
                0,
                "({d1}, {d2}): {} self-intersections in {} triangles",
                report.count(),
                patch.triangles.len()
            );
            checked += 1;
        }
    }
    // 20 patterns have 1 ≤ d₁ ≤ 5, 0 ≤ d₂ ≤ d₁. Exactly one of them — (4, 2),
    // with gcd 2 and ℓ = 12 — needs subdivision; every other is a quad, an
    // octagon or a single loop, and fills. A drop below 19 means a case
    // regressed into being unhandled rather than the assertion getting weaker.
    assert_eq!(checked, 19, "the sweep did not reach the implemented cases");
}

#[test]
fn the_single_loop_zero_is_vacuous_and_this_is_what_makes_that_visible() {
    // M-83 taken to its conclusion. A single loop is *one* fan, so every pair
    // of its triangles shares the apex, and `self_intersections` skips every
    // pair sharing a vertex. The zero reported for the single-loop case above
    // is therefore not evidence of anything, and the honest thing is to assert
    // that rather than let the sweep's green tick imply coverage it does not
    // have.
    //
    // The quad and octagon cases are different: those have loops that share no
    // vertex with each other, so their zeros are real. This asserts the
    // difference, so if the counter ever grows the ability to see intra-fan
    // folds, this test fails and the claim can be upgraded.
    let coords = pattern_coords(3, 2);
    let owned = crossings(coords.count);
    let mut patch = TetPatch::new();
    assert_eq!(fill(&tet(&owned), &mut patch), Ok(Unfilled::None));

    let indices: Vec<u32> = patch.triangles.iter().flatten().copied().collect();
    let report = crate::validate::self_intersections(&patch.positions, &indices, 1.0)
        .expect("a unit tet is a valid spacing");

    let n = patch.triangles.len() as u64;
    let pairs = n * (n - 1) / 2;
    assert_eq!(
        report.adjacent_pairs_skipped, pairs,
        "not every pair was skipped, so the counter can see *something* inside \
         a fan and this test's premise needs revisiting"
    );
    assert_eq!(report.count(), 0);
}

#[test]
fn reversing_the_steiner_assignment_is_visibly_wrong() {
    // M-44's rule applied to the test above: a zero has to prove it could have
    // been non-zero. Reversing "innermost to outermost" -- giving loop `i` the
    // Steiner point meant for loop `m - 1 - i` -- is precisely the mistake that
    // rule exists to prevent, and the nested fans then cross.
    //
    // The mutation has to keep the Steiner points **distinct**. Collapsing them
    // onto one apex also produces a wrong mesh, but every triangle would then
    // share that vertex and `self_intersections` skips pairs sharing a vertex
    // by construction (see its module docs, and M-83) -- so that version of the
    // mutation reports zero and proves nothing. Distinct-but-permuted keeps the
    // pairs visible.
    let m = 3u32;
    let coords = pattern_coords(m, m);
    let owned = crossings(coords.count);
    let mut patch = TetPatch::new();
    assert_eq!(fill(&tet(&owned), &mut patch), Ok(Unfilled::None));

    let base = 8 * m;
    let reversed: Vec<u32> = patch
        .triangles
        .iter()
        .flatten()
        .map(|i| {
            if *i >= base {
                base + (m - 1 - (*i - base))
            } else {
                *i
            }
        })
        .collect();

    let report = crate::validate::self_intersections(&patch.positions, &reversed, 1.0)
        .expect("a unit tet is a valid spacing");
    assert!(
        report.count() > 0,
        "reversing the Steiner assignment produced no intersections, so the \
         intersection-free assertion is not measuring anything"
    );
}

#[test]
fn a_cycle_visits_every_segment_once_and_returns_to_its_start() {
    for d1 in 1..=4u32 {
        for d2 in 0..=d1 {
            let coords = pattern_coords(d1, d2);
            for c in cycles(&coords) {
                let mut seen = c.points.clone();
                seen.sort_unstable();
                seen.dedup();
                assert_eq!(
                    seen.len(),
                    c.points.len(),
                    "({d1}, {d2}): a cycle repeated a point"
                );
            }
        }
    }
}

#[test]
fn an_unimplemented_case_is_named_rather_than_guessed_at() {
    // (4, 2): gcd = 2, so two loops of length 4(6)/2 = 12 -- §3.2.1 case (3),
    // subdivision, not implemented. The contract is that fill names which case
    // it reached, and emits nothing for it.
    for (d1, d2, expected) in [
        (4u32, 2u32, Unfilled::Subdivision),
        (6, 3, Unfilled::Subdivision),
        (6, 4, Unfilled::Subdivision),
    ] {
        let coords = pattern_coords(d1, d2);
        let owned = crossings(coords.count);
        let mut patch = TetPatch::new();
        assert_eq!(
            fill(&tet(&owned), &mut patch),
            Ok(expected),
            "({d1}, {d2}) should report {expected:?}"
        );
        assert!(
            patch.is_empty(),
            "({d1}, {d2}): an unhandled case emitted triangles"
        );
    }
}

#[test]
fn an_octagon_is_a_fan_around_its_steiner_point() {
    // (1, 1): gcd = 1, length 4(2)/1 = 8. One octagon, one Steiner point, eight
    // triangles -- one per loop edge.
    let coords = pattern_coords(1, 1);
    let cycles = cycles(&coords);
    assert_eq!(cycles.len(), 1);
    assert_eq!(cycles[0].length(), 8);

    let owned = crossings(coords.count);
    let mut patch = TetPatch::new();
    assert_eq!(fill(&tet(&owned), &mut patch), Ok(Unfilled::None));
    assert_eq!(patch.triangles.len(), 8);
    // Eight crossings plus one Steiner point.
    assert_eq!(patch.positions.len(), 9);

    // Every triangle uses the Steiner point exactly once, and the loop edges it
    // fans over are all eight of them.
    let steiner = patch.positions.len() as u32 - 1;
    for t in &patch.triangles {
        assert_eq!(
            t.iter().filter(|i| **i == steiner).count(),
            1,
            "triangle {t:?} does not fan around the Steiner point"
        );
    }
}

#[test]
fn m_octagons_get_m_steiner_points_ordered_along_the_axis() {
    // (m, m) for m > 1 is the nested-octagon case: m loops, each of length 8,
    // each fanned around its own Steiner point.
    for m in 1..=4u32 {
        let coords = pattern_coords(m, m);
        let owned = crossings(coords.count);
        let mut patch = TetPatch::new();
        assert_eq!(
            fill(&tet(&owned), &mut patch),
            Ok(Unfilled::None),
            "({m}, {m}) should be the octagon case"
        );

        let loops = m as usize;
        assert_eq!(patch.triangles.len(), 8 * loops, "m = {m}");
        // 8 crossings per loop, plus one Steiner point per loop.
        assert_eq!(patch.positions.len(), 8 * loops + loops, "m = {m}");

        // The Steiner points are the last `m` positions and must be strictly
        // ordered along the segment they were placed on -- the property the
        // intersection-free argument rests on, and the one that does not depend
        // on the particular spacing.
        let steiner = &patch.positions[8 * loops..];
        assert_eq!(steiner.len(), loops);
        let axis = [
            steiner[loops - 1][0] - steiner[0][0],
            steiner[loops - 1][1] - steiner[0][1],
            steiner[loops - 1][2] - steiner[0][2],
        ];
        if m > 1 {
            let mut previous = f64::NEG_INFINITY;
            for s in steiner {
                let t = s[0] * axis[0] + s[1] * axis[1] + s[2] * axis[2];
                assert!(t > previous, "m = {m}: Steiner points are not ordered");
                previous = t;
            }
        }
    }
}

#[test]
fn each_octagon_pairs_nested_crossings_on_the_long_edge() {
    // The Steiner assignment reads "innermost to outermost pair", which is only
    // meaningful if each loop's two crossings on the 2m edge are symmetric
    // about that edge's middle: p_j pairs with p_{2m-1-j}. fill() refuses to
    // emit when that fails, so checking it here is checking that the refusal
    // never fires on a well-formed pattern -- and the nesting itself.
    for m in 2..=4u32 {
        let coords = pattern_coords(m, m);
        let cycles = cycles(&coords);
        assert_eq!(cycles.len(), m as usize);

        // The 2m pair is (2, 3) by pattern_coords' construction: d1 + d2 = 2m.
        let long = 2u8;
        assert_eq!(coords.edge(long), 2 * m);

        let mut ranks: Vec<[u32; 2]> = Vec::new();
        for c in &cycles {
            let mut on_long: Vec<u32> = c
                .points
                .iter()
                .filter(|q| q.edge == long)
                .map(|q| q.index)
                .collect();
            on_long.sort_unstable();
            assert_eq!(on_long.len(), 2, "m = {m}: a loop did not cross e twice");
            assert_eq!(
                on_long[0] + on_long[1],
                2 * m - 1,
                "m = {m}: crossings {on_long:?} are not nested about the middle"
            );
            ranks.push([on_long[0], on_long[1]]);
        }
        // Every nesting level is used exactly once.
        ranks.sort_unstable();
        for (i, r) in ranks.iter().enumerate() {
            assert_eq!(r[0], i as u32, "m = {m}: nesting levels are not distinct");
        }
    }
}

#[test]
fn filling_is_deterministic_and_reset_keeps_capacity() {
    let coords = pattern_coords(2, 0);
    let owned = crossings(coords.count);
    let t = tet(&owned);

    let mut a = TetPatch::new();
    let mut b = TetPatch::new();
    assert_eq!(fill(&t, &mut a), Ok(Unfilled::None));
    assert_eq!(fill(&t, &mut b), Ok(Unfilled::None));
    assert_eq!(a.positions, b.positions);
    assert_eq!(a.triangles, b.triangles);

    let capacity = a.positions.capacity();
    a.reset();
    assert!(a.is_empty());
    assert_eq!(a.positions.capacity(), capacity);
}

#[test]
fn no_triangle_indexes_a_vertex_that_does_not_exist_or_repeats_one() {
    for d1 in 1..=4u32 {
        for d2 in 0..=d1 {
            let coords = pattern_coords(d1, d2);
            let owned = crossings(coords.count);
            let mut patch = TetPatch::new();
            let outcome = fill(&tet(&owned), &mut patch).expect("well-formed crossings");
            assert_ne!(outcome, Unfilled::Inconsistent);
            assert_ne!(outcome, Unfilled::NoPattern);

            let n = patch.positions.len() as u32;
            for t in &patch.triangles {
                assert!(t.iter().all(|i| *i < n), "({d1}, {d2}): index out of range");
                assert!(
                    t[0] != t[1] && t[1] != t[2] && t[0] != t[2],
                    "({d1}, {d2}): degenerate triangle {t:?}"
                );
            }
        }
    }
}
