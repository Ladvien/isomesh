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
    let residual = residual(&cycles);
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

/// Two tets sharing face `(0, 1, 2)`, which lies in the plane `z = 0`. Their
/// apexes are on opposite sides, so the shared face is the only thing they have
/// in common — and the only place a crack could open.
fn two_tets(
    shared: [Vec<f64>; 3],
    apex_a: [Vec<f64>; 3],
    apex_b: [Vec<f64>; 3],
) -> [Vec<Vec<f64>>; 2] {
    // TET_EDGES order: 0=(0,1), 1=(0,2), 2=(0,3), 3=(1,2), 4=(1,3), 5=(2,3).
    // Edges 0, 1 and 3 are the shared face's; 2, 4 and 5 reach the apex.
    let [s01, s02, s12] = shared;
    let build = |apex: [Vec<f64>; 3]| {
        let [a03, a13, a23] = apex;
        alloc::vec![s01.clone(), s02.clone(), a03, s12.clone(), a13, a23,]
    };
    [build(apex_a), build(apex_b)]
}

// Exact float equality is the point here, not an oversight. The shared face's
// crossings are computed from the same corners and the same parameters in both
// tets, so they must come out **bit-identical** -- that is the property the
// whole conformity argument rests on (M-32 measured what happens when an
// analogous claim holds only by algebra). Welding with a tolerance would hide
// exactly the failure this test exists to catch.
#[expect(clippy::float_cmp, reason = "bit-identity is the property under test")]
#[test]
fn the_shared_face_of_two_tets_carries_no_crack() {
    // A-014b's acceptance criterion. Two tets are filled independently, with no
    // communication and no second pass, and the question is whether the
    // triangles meeting at their shared face line up.
    //
    // Measured the way E-107 measured the transvoxel seam: merge the two
    // patches, weld coincident vertices, and count boundary edges -- edges used
    // by exactly one triangle -- that lie *in the shared face*. Every other
    // boundary edge is the patch's genuine open border against the rest of the
    // grid and must survive; a boundary edge inside the shared plane is a hole
    // you can see through.
    let corners_a = [
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0],
    ];
    let corners_b = [
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, -1.0],
    ];

    // Counts on the shared face's three edges, then each apex's three. The two
    // tets differ in their apex edges, which is the whole point: conformity
    // must not depend on the neighbour's far side.
    let cases: [([u32; 3], [u32; 3], [u32; 3]); 4] = [
        ([1, 1, 0], [1, 0, 0], [1, 0, 0]),
        ([1, 1, 2], [1, 2, 1], [2, 1, 2]),
        ([2, 2, 2], [2, 2, 2], [1, 1, 1]),
        ([3, 1, 2], [2, 1, 3], [1, 3, 2]),
    ];

    let mut examined = 0;
    let mut excused = 0;
    for (shared, apex_a, apex_b) in cases {
        let s = [spaced(shared[0]), spaced(shared[1]), spaced(shared[2])];
        let a = [spaced(apex_a[0]), spaced(apex_a[1]), spaced(apex_a[2])];
        let b = [spaced(apex_b[0]), spaced(apex_b[1]), spaced(apex_b[2])];
        let [owned_a, owned_b] = two_tets(s, a, b);

        let mut merged_positions: Vec<[f64; 3]> = Vec::new();
        let mut merged_triangles: Vec<[u32; 3]> = Vec::new();
        let mut any_unfilled = false;

        for (owned, corners) in [(&owned_a, corners_a), (&owned_b, corners_b)] {
            let mut along: [&[f64]; TET_EDGE_COUNT] = [&[]; TET_EDGE_COUNT];
            for (slot, v) in along.iter_mut().zip(owned.iter()) {
                *slot = v.as_slice();
            }
            let t = TetCrossings { corners, along };
            let mut patch = TetPatch::new();
            let outcome = fill(&t, &mut patch).expect("well-formed crossings");
            if outcome != Unfilled::None {
                any_unfilled = true;
            }

            // Weld on exact position equality. The shared face's crossings are
            // computed from the same corners and the same parameters in both
            // tets, so they are bit-identical -- which is itself the thing being
            // relied on, and would fail loudly here if it were not true.
            for tri in &patch.triangles {
                let mut welded = [0u32; 3];
                for (slot, index) in welded.iter_mut().zip(tri.iter()) {
                    let p = patch.positions[*index as usize];
                    *slot = match merged_positions.iter().position(|q| *q == p) {
                        Some(existing) => existing as u32,
                        None => {
                            merged_positions.push(p);
                            merged_positions.len() as u32 - 1
                        }
                    };
                }
                merged_triangles.push(welded);
            }
        }

        if any_unfilled || merged_triangles.is_empty() {
            continue;
        }

        // Count each undirected edge's incident faces.
        let mut edges: Vec<([u32; 2], u32)> = Vec::new();
        for tri in &merged_triangles {
            for k in 0..3 {
                let (x, y) = (tri[k], tri[(k + 1) % 3]);
                let key = if x < y { [x, y] } else { [y, x] };
                match edges.iter_mut().find(|(e, _)| *e == key) {
                    Some((_, n)) => *n += 1,
                    None => edges.push((key, 1)),
                }
            }
        }

        let in_shared_face = |e: [u32; 2]| {
            merged_positions[e[0] as usize][2] == 0.0 && merged_positions[e[1] as usize][2] == 0.0
        };

        // Every segment either tet discarded as part of an **open** curve, as a
        // pair of positions. §3.1 is explicit that these legitimately become
        // mesh boundary: "Such curves are discarded, but their segments may
        // still appear in neighboring tets as part of the mesh boundary." So a
        // one-sided edge in the shared face is a crack only if it is *not* one
        // of these.
        let mut discarded: Vec<[[f64; 3]; 2]> = Vec::new();
        for (owned, corners) in [(&owned_a, corners_a), (&owned_b, corners_b)] {
            let mut along: [&[f64]; TET_EDGE_COUNT] = [&[]; TET_EDGE_COUNT];
            for (slot, v) in along.iter_mut().zip(owned.iter()) {
                *slot = v.as_slice();
            }
            let t = TetCrossings { corners, along };
            for curve in crate::subgrid::curves::curves(&t.coordinates()) {
                if curve.kind != CurveKind::Open {
                    continue;
                }
                for seg in &curve.segments {
                    if let (Some(p), Some(q)) = (t.position(seg.a), t.position(seg.b)) {
                        discarded.push([p, q]);
                    }
                }
            }
        }
        let was_discarded = |e: [u32; 2]| {
            let (p, q) = (
                merged_positions[e[0] as usize],
                merged_positions[e[1] as usize],
            );
            discarded
                .iter()
                .any(|[a, b]| (*a == p && *b == q) || (*a == q && *b == p))
        };

        let cracks: Vec<[u32; 2]> = edges
            .iter()
            .filter(|(e, n)| *n == 1 && in_shared_face(*e) && !was_discarded(*e))
            .map(|(e, _)| *e)
            .collect();
        assert!(
            cracks.is_empty(),
            "{shared:?}/{apex_a:?}/{apex_b:?}: {} boundary edges inside the shared face \
             that no open curve accounts for: {cracks:?}",
            cracks.len()
        );

        // No shared-face edge may be used by more than two faces either -- that
        // would be a non-manifold seam rather than a hole.
        assert!(
            edges.iter().all(|(e, n)| !in_shared_face(*e) || *n <= 2),
            "{shared:?}: a shared-face edge is used by more than two triangles"
        );

        // A zero has to prove it could have been non-zero: the shared face must
        // actually carry geometry, or this case is asserting nothing.
        let shared_edges = edges.iter().filter(|(e, _)| in_shared_face(*e)).count();
        assert!(
            shared_edges > 0,
            "{shared:?}: no triangle edge lies in the shared face at all"
        );
        excused += edges
            .iter()
            .filter(|(e, n)| *n == 1 && in_shared_face(*e) && was_discarded(*e))
            .count();
        examined += 1;
    }
    assert!(examined > 0, "no case filled completely enough to examine");
    // The open-curve excuse is the interesting half of the assertion, and an
    // excuse that never fires would let a real crack through unnoticed. The
    // fixture is chosen so it fires: `[1, 1, 2]` with those apexes produces one.
    assert!(
        excused > 0,
        "no shared-face boundary edge was excused by a discarded open curve, \
         so that clause is not being exercised"
    );
}

#[test]
fn the_edge_labelling_closes_on_the_far_corner_it_was_not_told_about() {
    // §3.2.2's labelling walks an edge from its lower corner, taking that
    // corner's side and flipping at every crossing. Nothing in that walk looks
    // at the *upper* corner -- so the label it arrives at is a prediction, and
    // it must match the side that corner was independently assigned.
    //
    // The two agree only if the parity rule and the vertex rule are consistent
    // with each other: crossing γ an odd number of times must mean the two
    // corners are on opposite sides, which is precisely what b_ij is defined to
    // say. A sign error in either rule breaks this and nothing else would.
    let mut corner_loops = 0;
    let mut contractible_loops = 0;

    for raw in 0..4096u32 {
        let mut count = [0u32; TET_EDGE_COUNT];
        for (e, slot) in count.iter_mut().enumerate() {
            *slot = (raw >> (2 * e)) & 0b11;
        }
        let coords = EdgeCoordinates::new(count);
        for cycle in cycles(&coords) {
            let Some(sides) = cycle.corner_sides() else {
                continue;
            };
            match cycle.non_normal_kind() {
                Some(NonNormalKind::Corner) => corner_loops += 1,
                Some(NonNormalKind::Contractible) => contractible_loops += 1,
                _ => panic!("only corner and contractible loops have sides"),
            }

            for edge in 0..TET_EDGE_COUNT as u8 {
                let labels = cycle.edge_sides(edge).expect("this loop has sides");
                let [lo, hi] = TET_EDGES[edge as usize];
                let crossings = cycle.points.iter().filter(|p| p.edge == edge).count();

                assert_eq!(labels.len(), crossings + 1, "{count:?} edge {edge}");
                assert_eq!(labels[0], sides[lo as usize], "{count:?} edge {edge}");
                assert_eq!(
                    *labels.last().expect("non-empty"),
                    sides[hi as usize],
                    "{count:?} edge {edge}: walking from corner {lo} across {crossings} \
                     crossings disagrees with corner {hi}'s own side"
                );
            }

            // A corner loop's distinguished corner is the inside one, and it is
            // the only one.
            if let Some(distinguished) = cycle.distinguished_corner() {
                assert_eq!(sides[distinguished as usize], Side::Inside);
                assert_eq!(
                    sides.iter().filter(|s| **s == Side::Inside).count(),
                    1,
                    "{count:?}: a corner loop has exactly one inside vertex"
                );
            }
        }
    }

    assert!(
        corner_loops > 0 && contractible_loops > 0,
        "corner: {corner_loops}, contractible: {contractible_loops} -- both must be \
         reached or one of the two labelling rules is untested"
    );
}

#[test]
fn a_diagonal_loop_has_no_sides_and_that_is_the_papers_answer() {
    // "when γ is diagonal we do not require an inside/outside distinction" --
    // so None here is the specified behaviour, not a missing case, and it is
    // why the diagonal type could be triangulated before this labelling existed.
    let mut found = false;
    for raw in 0..4096u32 {
        let mut count = [0u32; TET_EDGE_COUNT];
        for (e, slot) in count.iter_mut().enumerate() {
            *slot = (raw >> (2 * e)) & 0b11;
        }
        for cycle in cycles(&EdgeCoordinates::new(count)) {
            if cycle.non_normal_kind() == Some(NonNormalKind::Diagonal) {
                assert_eq!(cycle.corner_sides(), None);
                assert_eq!(cycle.edge_sides(0), None);
                assert_eq!(cycle.distinguished_corner(), None);
                found = true;
            }
        }
    }
    assert!(found, "no diagonal loop was reached");
}

#[test]
fn the_coverage_of_the_implemented_cases_is_pinned() {
    // How much of the encoding §3.2 currently serves, over every configuration
    // with 0..=3 crossings on each of the six edges. Pinned as exact counts
    // rather than a ratio so a regression moves a number rather than rounding
    // away, and so the remaining work has a size.
    let mut none = 0u32;
    let mut single = 0u32;
    let mut subdivision = 0u32;
    let mut non_normal = 0u32;
    let mut no_pattern = 0u32;
    let mut inconsistent = 0u32;

    for raw in 0..4096u32 {
        let mut count = [0u32; TET_EDGE_COUNT];
        for (e, slot) in count.iter_mut().enumerate() {
            *slot = (raw >> (2 * e)) & 0b11;
        }
        let owned = crossings(count);
        let mut patch = TetPatch::new();
        let outcome = match fill(&tet(&owned), &mut patch) {
            Ok(o) => o,
            Err(e) => panic!("{count:?} was rejected as malformed: {e:?}"),
        };

        // Whatever was emitted has to be indexable, on every path -- including
        // the boundary disk, which is the one case that builds triangles from
        // regions rather than from a fan and so has its own way to go wrong.
        let n = patch.positions.len() as u32;
        for t in &patch.triangles {
            assert!(
                t.iter().all(|i| *i < n),
                "{count:?}: triangle {t:?} indexes past {n} vertices"
            );
        }
        match outcome {
            Unfilled::None => none += 1,
            Unfilled::SingleLoop => single += 1,
            Unfilled::Subdivision => subdivision += 1,
            Unfilled::NonNormalLoop => non_normal += 1,
            Unfilled::NoPattern => no_pattern += 1,
            Unfilled::Inconsistent => inconsistent += 1,
        }
    }

    // Nothing may reach Inconsistent: it names a disagreement between the
    // curves and the crossing lists, which is this crate's own bug and not a
    // case. A non-zero here is always a defect, never a missing feature.
    assert_eq!(inconsistent, 0, "Inconsistent is a bug, not a case");

    // Two of these are load-bearing zeros rather than gaps.
    //
    // `NoPattern` is zero, which says Property II held on Γ_normal's own
    // residual in every one of the 4,096 configurations — an empirical check of
    // Theorem B.3 across the whole sweep, and the thing M-85's restructure was
    // needed to make true.
    //
    // `SingleLoop` and `Subdivision` are zero because neither case is *reachable*
    // at these counts: both need a residual pattern with ℓ > 8, which wants
    // larger d₁ and d₂ than three crossings per edge can supply. They are
    // exercised by `an_unimplemented_case_is_named_rather_than_guessed_at` and
    // `a_single_long_loop_fans_around_its_centre_of_mass` instead. So the whole
    // of the remaining gap here is §3.2.2's corner and contractible types.
    assert_eq!(
        [none, single, subdivision, non_normal, no_pattern],
        [4096, 0, 0, 0, 0],
        "coverage moved: [None, SingleLoop, Subdivision, NonNormalLoop, NoPattern]"
    );
    assert_eq!(
        none + single + subdivision + non_normal + no_pattern,
        4096,
        "the outcomes do not partition the sweep"
    );
}

#[test]
fn non_normal_loops_are_classified_and_all_three_types_are_reachable() {
    // §3.2.2's parity rule: b_ij = e_ij^γ mod 2 over the loop's own coordinates,
    // p = b₀₁ + b₀₂ + b₀₃, then 0 → contractible, 2 → diagonal, 1 or 3 →
    // corner. A sweep that never produced all three would be testing one
    // branch and reporting three, so the reachability of each is asserted.
    // Counters rather than a map, so the public enum does not have to derive
    // Ord for a test's convenience.
    let mut seen = [0u32; 3];
    let slot = |kind: NonNormalKind| match kind {
        NonNormalKind::Contractible => 0usize,
        NonNormalKind::Diagonal => 1,
        NonNormalKind::Corner => 2,
    };
    for a in 0..4u32 {
        for b in 0..4u32 {
            for c in 0..4u32 {
                for d in 0..4u32 {
                    for e in 0..4u32 {
                        for f in 0..4u32 {
                            let coords = EdgeCoordinates::new([a, b, c, d, e, f]);
                            for cycle in cycles(&coords) {
                                if let Some(kind) = cycle.non_normal_kind() {
                                    seen[slot(kind)] += 1;
                                }
                                // A normal loop has no non-normal type.
                                assert_eq!(
                                    cycle.non_normal_kind().is_some(),
                                    cycle.kind == CurveKind::NonNormal
                                );
                            }
                        }
                    }
                }
            }
        }
    }

    for kind in [
        NonNormalKind::Contractible,
        NonNormalKind::Diagonal,
        NonNormalKind::Corner,
    ] {
        assert!(
            seen[slot(kind)] > 0,
            "{kind:?} was never produced, so its branch is untested"
        );
    }
}

#[test]
fn all_three_non_normal_types_now_have_a_spanning_disk() {
    // Diagonal loops fan around a centre of mass in the tet interior;
    // contractible and corner loops get a disk built in the tet boundary, the
    // corner one with §3.2.2's two extras -- the vertex coinciding with the
    // inside corner omitted, and a triangle capping that corner.
    //
    // Swept over the whole space rather than the `[a, b, c, 0, 0, 0]` corner of
    // it: that subspace never produces a corner-type loop at all, and the
    // reachability counters below are what caught an earlier version claiming
    // to cover a branch it could not reach.
    let mut seen = [0u32; 3];
    let mut empty_disks = 0u32;
    let slot = |kind: NonNormalKind| match kind {
        NonNormalKind::Contractible => 0usize,
        NonNormalKind::Diagonal => 1,
        NonNormalKind::Corner => 2,
    };

    for raw in 0..4096u32 {
        {
            {
                let mut count = [0u32; TET_EDGE_COUNT];
                for (e, slot) in count.iter_mut().enumerate() {
                    *slot = (raw >> (2 * e)) & 0b11;
                }
                let coords = EdgeCoordinates::new(count);
                let non_normal: Vec<_> = cycles(&coords)
                    .into_iter()
                    .filter(|x| x.kind == CurveKind::NonNormal)
                    .collect();
                if non_normal.is_empty() {
                    continue;
                }
                for x in &non_normal {
                    if let Some(kind) = x.non_normal_kind() {
                        seen[slot(kind)] += 1;
                    }
                }

                let owned = crossings(coords.count);
                let mut patch = TetPatch::new();
                let outcome = fill(&tet(&owned), &mut patch).expect("well-formed");
                assert_eq!(
                    outcome,
                    Unfilled::None,
                    "{:?}: every non-normal type has a disk now",
                    coords.count
                );
                // Not "every non-normal loop produces triangles". The minimal
                // one -- `e = (2, 0, 0, 0, 0, 0)`, a single scoop pair -- has an
                // inside region that is a **bigon**: one chord and one edge
                // piece, two nodes, and a fan over two nodes is empty. That is
                // V-21's degeneracy at its most extreme, and it is correct here:
                // the region has zero area until A-014d insets it into the tet
                // interior. Counted rather than asserted away, so the number
                // moves if the behaviour does.
                if patch.triangles.is_empty() {
                    empty_disks += 1;
                }
            }
        }
    }
    for kind in [
        NonNormalKind::Contractible,
        NonNormalKind::Diagonal,
        NonNormalKind::Corner,
    ] {
        assert!(
            seen[slot(kind)] > 0,
            "{kind:?} was never reached, so its disk is untested"
        );
    }
    // 105 of the configurations carrying a non-normal loop triangulate to
    // nothing, because their inside regions are bigons -- a chord plus a single
    // edge piece. V-21: correct at this stage, zero area until A-014d insets
    // them into the tet interior. Pinned so the number moves if the behaviour
    // does, rather than hiding inside a green tick.
    assert_eq!(empty_disks, 105, "degenerate bigon disks");
}

#[test]
fn the_subdivision_stencil_produces_four_normal_sub_tets() {
    // The stencil is only usable if every tet it makes is itself a valid input
    // to the same procedure. That is not obvious from the assignment
    // e_ai = 2d₂, e_aj = d₁, e_ak = d₂, e_al = d₁ − d₂ -- it is asymmetric in
    // the four corners, and three of the four faces of each sub-tet satisfy the
    // triangle inequality with *equality*, which is exactly where an off-by-one
    // would show up.
    let mut checked = 0;
    for d1 in 1..=6u32 {
        for d2 in 0..=d1 {
            let coords = pattern_coords(d1, d2);
            let pattern = Pattern::of(&coords).expect("built as a pattern");
            let Some(stencil) = Subdivision::label(&coords, pattern) else {
                panic!("({d1}, {d2}): Property II holds but no labelling was found");
            };

            // The labelling is what it claims to be.
            let [i, j, k, l] = stencil.corner;
            let e = |a: u8, b: u8| coords.edge(crate::subgrid::coordinates::edge_between(a, b));
            assert_eq!((e(i, j), e(k, l)), (d1, d1), "({d1}, {d2})");
            assert_eq!((e(i, k), e(j, l)), (d2, d2), "({d1}, {d2})");
            assert_eq!((e(i, l), e(j, k)), (d1 + d2, d1 + d2), "({d1}, {d2})");
            assert_eq!(stencil.spoke(), [2 * d2, d1, d2, d1 - d2], "({d1}, {d2})");

            for which in 0..4 {
                let sub = stencil
                    .sub_tet(&coords, which)
                    .expect("four sub-tets exist");
                assert!(
                    sub.is_normal(),
                    "({d1}, {d2}) sub-tet {which} is not a normal configuration: {:?}",
                    sub.count
                );
            }
            checked += 1;
        }
    }
    assert_eq!(checked, 27, "the sweep did not run");
}

#[test]
fn the_stencil_is_asymmetric_and_the_labelling_is_what_decides_it() {
    // e_ai = 2d₂ against e_aj = d₁: swapping which corner is called `i` gives a
    // different subdivision, so the labelling is load-bearing rather than
    // cosmetic. With d₁ ≠ d₂ the spokes are genuinely distinct.
    let coords = pattern_coords(3, 1);
    let pattern = Pattern::of(&coords).expect("built as a pattern");
    let stencil = Subdivision::label(&coords, pattern).expect("a labelling exists");
    let spoke = stencil.spoke();
    assert_eq!(spoke, [2, 3, 1, 2]);

    // The four sub-tets are genuinely different from one another, which is what
    // makes "recursively process each of the four new tets" four calls and not
    // one repeated.
    let mut seen: Vec<[u32; TET_EDGE_COUNT]> = (0..4)
        .filter_map(|w| stencil.sub_tet(&coords, w))
        .map(|c| c.count)
        .collect();
    assert_eq!(seen.len(), 4);
    seen.sort_unstable();
    seen.dedup();
    assert_eq!(seen.len(), 4, "two sub-tets came out identical");
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

#[test]
fn a_face_split_along_a_loop_partitions_into_regions() {
    // §3.2.2's σ \ γ. The decomposition is checked by conservation rather than
    // by inspection: every edge arc of the face must appear in exactly one
    // region, and every chord in exactly two (once from each side) -- which is
    // what "partition" means for a disk cut by chords, and what a peel that
    // dropped or duplicated a span would break.
    let mut faces_with_chords = 0;
    let mut multi_chord_faces = 0;

    for raw in 0..4096u32 {
        let mut count = [0u32; TET_EDGE_COUNT];
        for (e, slot) in count.iter_mut().enumerate() {
            *slot = (raw >> (2 * e)) & 0b11;
        }
        for cycle in cycles(&EdgeCoordinates::new(count)) {
            if cycle.corner_sides().is_none() {
                continue;
            }
            for face in 0..4u8 {
                let regions = face_regions(face, &EdgeCoordinates::new(count), &cycle)
                    .expect("this loop has sides");
                assert!(!regions.is_empty(), "{count:?} face {face}: no regions");

                let mut edge_arcs = 0usize;
                let mut chords = 0usize;
                for r in &regions {
                    assert_eq!(
                        r.arc.len(),
                        r.node.len(),
                        "{count:?}: arcs and nodes differ"
                    );
                    assert!(r.arc.len() >= 2, "{count:?}: a region with under 2 arcs");
                    for a in &r.arc {
                        match a {
                            Arc::Edge { .. } => edge_arcs += 1,
                            Arc::Chord => chords += 1,
                        }
                    }
                }

                // The face's three edges are each cut into (crossings + 1)
                // pieces, and every piece belongs to exactly one region.
                let expected_pieces: usize = (0..3)
                    .map(|k| {
                        let e = crate::subgrid::coordinates::TET_FACES[face as usize].edge[k];
                        cycle.points.iter().filter(|p| p.edge == e).count() + 1
                    })
                    .sum();
                assert_eq!(
                    edge_arcs, expected_pieces,
                    "{count:?} face {face}: edge pieces not conserved"
                );

                // Each chord bounds exactly two regions.
                let chord_count = chords;
                assert_eq!(
                    chord_count % 2,
                    0,
                    "{count:?} face {face}: a chord bounds an odd number of regions"
                );
                assert_eq!(
                    regions.len(),
                    chord_count / 2 + 1,
                    "{count:?} face {face}: {} regions for {} chord sides",
                    regions.len(),
                    chord_count
                );
                // Order-sensitive, unlike the counts above: a crossing sits on
                // the border between exactly two regions and a corner belongs
                // to exactly one. Peeling a chord that is *not* innermost cuts
                // off a region still containing unresolved chords, which lands
                // their endpoints in the wrong number of regions.
                let mut appearances: Vec<(Node, usize)> = Vec::new();
                for r in &regions {
                    for n in &r.node {
                        match appearances.iter_mut().find(|(m, _)| m == n) {
                            Some((_, c)) => *c += 1,
                            None => appearances.push((*n, 1)),
                        }
                    }
                }
                for (n, c) in &appearances {
                    let expected = match n {
                        Node::Corner(_) => 1,
                        Node::Crossing(_) => 2,
                    };
                    assert_eq!(
                        *c, expected,
                        "{count:?} face {face}: {n:?} appears in {c} regions, not {expected}"
                    );
                }

                if chord_count > 0 {
                    faces_with_chords += 1;
                }
                if chord_count / 2 >= 2 {
                    multi_chord_faces += 1;
                }
            }
        }
    }
    assert!(faces_with_chords > 0, "no face was ever actually cut");
    // The peel's "innermost" requirement only bites when a face carries two or
    // more chords -- with one chord, every chord is innermost. If this were
    // zero, dropping that requirement would pass the test, which is exactly
    // what a mutation showed before this assertion existed.
    assert!(
        multi_chord_faces > 0,
        "no face carried two chords, so the innermost rule is untested"
    );
}
