//! The encoding checked against the paper's own examples, and against A-003's
//! table.

use alloc::vec::Vec;

use super::*;
use crate::marching_tetrahedra::table::TET_CASES;

/// The paper writes `e := (e₀₁, e₀₂, e₀₃, e₂₃, e₁₃, e₁₂)`; this crate uses
/// [`TET_EDGES`]' lexicographic order. Translate one of the paper's vectors.
fn from_paper_order(paper: [u32; 6]) -> EdgeCoordinates {
    let mut count = [0u32; TET_EDGE_COUNT];
    for (value, pair) in paper
        .into_iter()
        .zip([[0, 1], [0, 2], [0, 3], [2, 3], [1, 3], [1, 2]])
    {
        count[edge_between(pair[0], pair[1]) as usize] = value;
    }
    EdgeCoordinates { count }
}

// ─── the geometry the encoding sits on ──────────────────────────────────────

/// The whole reason the ordering can be changed safely: `5 - edge` really is the
/// complementary pairing under this crate's edge order.
///
/// The quad basis vectors are *defined* by which opposite pair they separate, so
/// if this were wrong `decompose` would produce the wrong polygons for every
/// input and nothing else here would notice.
#[test]
fn complementary_edges_share_no_corner() {
    for edge in 0..TET_EDGE_COUNT as u8 {
        let other = complementary(edge);
        assert_ne!(edge, other);
        assert_eq!(
            complementary(other),
            edge,
            "the pairing must be an involution"
        );

        let a = TET_EDGES[edge as usize];
        let b = TET_EDGES[other as usize];
        for corner in a {
            assert!(
                !b.contains(&corner),
                "edges {edge} {a:?} and {other} {b:?} share corner {corner}"
            );
        }
        // Together they must cover all four corners.
        let mut all: Vec<u8> = a.into_iter().chain(b).collect();
        all.sort_unstable();
        assert_eq!(all, alloc::vec![0, 1, 2, 3]);
    }
}

#[test]
fn every_face_names_three_corners_and_their_three_edges() {
    let mut corner_uses = [0u32; 4];
    let mut edge_uses = [0u32; TET_EDGE_COUNT];
    for (f, face) in TET_FACES.iter().enumerate() {
        // Face f omits corner f.
        assert!(
            !face.corner.contains(&(f as u8)),
            "face {f} contains corner {f}"
        );
        for c in face.corner {
            corner_uses[c as usize] += 1;
        }
        for k in 0..3 {
            let expected = edge_between(face.corner[k], face.corner[(k + 1) % 3]);
            assert_eq!(face.edge[k], expected, "face {f} edge {k}");
            edge_uses[face.edge[k] as usize] += 1;
        }
    }
    // Each corner is on three faces; each edge is on two.
    assert_eq!(corner_uses, [3; 4]);
    assert_eq!(edge_uses, [2; TET_EDGE_COUNT]);
}

// ─── section 2.1, the normal-curve conditions ───────────────────────────────

/// `eᵢⱼ = cᵢ + cⱼ` and its inverse must be inverses, on every face and over a
/// wide range of corner coordinates.
#[test]
fn corner_and_edge_coordinates_round_trip() {
    for c0 in 0..6u32 {
        for c1 in 0..6u32 {
            for c2 in 0..6u32 {
                // Build a tet whose face 3 (corners 0, 1, 2) has these corner
                // coordinates, and nothing on the edges to corner 3.
                let mut count = [0u32; TET_EDGE_COUNT];
                count[edge_between(0, 1) as usize] = c0 + c1;
                count[edge_between(1, 2) as usize] = c1 + c2;
                count[edge_between(0, 2) as usize] = c2 + c0;
                let e = EdgeCoordinates { count };

                let recovered = e
                    .corner_coordinates(3)
                    .expect("a constructed face is normal");
                assert_eq!(recovered, [c0, c1, c2], "c = ({c0}, {c1}, {c2})");
            }
        }
    }
}

/// Both conditions reject exactly what the paper says they are for, and neither
/// rejects a case the other should catch.
#[test]
fn the_two_normality_conditions_reject_what_they_are_for() {
    // Odd sum: "what goes in must come out". One arc entering face 3 and none
    // leaving it.
    let mut count = [0u32; TET_EDGE_COUNT];
    count[edge_between(0, 1) as usize] = 1;
    let odd = EdgeCoordinates { count };
    assert!(matches!(odd.validate(), Err(NotNormal::OddSum { .. })));
    assert!(!odd.is_normal());

    // Triangle inequality: "arcs entering one edge must exit a different edge".
    // Four crossings on 01 with only one on each edge that could receive them.
    //
    // The first attempt at this fixture set three edges and left the rest at
    // zero, which made *another* face odd and tripped the sum check before the
    // inequality was ever reached. Every face here sums even on purpose, so the
    // only thing left to fail is the inequality.
    let mut count = [0u32; TET_EDGE_COUNT];
    count[edge_between(0, 1) as usize] = 4;
    count[edge_between(0, 2) as usize] = 1;
    count[edge_between(1, 2) as usize] = 1;
    count[edge_between(0, 3) as usize] = 1;
    count[edge_between(1, 3) as usize] = 1;
    count[edge_between(2, 3) as usize] = 0;
    let pinched = EdgeCoordinates { count };
    for face in 0..TET_FACE_COUNT as u8 {
        let f = TET_FACES[face as usize];
        let sum: u32 = f.edge.iter().map(|e| pinched.count[*e as usize]).sum();
        assert_eq!(
            sum % 2,
            0,
            "face {face} was left odd, so this fixture is wrong"
        );
    }
    assert!(
        matches!(
            pinched.validate(),
            Err(NotNormal::TriangleInequality { .. })
        ),
        "{:?}",
        pinched.validate()
    );

    // And the empty tet is normal, trivially.
    assert!(EdgeCoordinates::empty().is_normal());
}

// ─── section 2.3, the decomposition ─────────────────────────────────────────

/// Every collection of normal polygons must decompose back to itself.
///
/// Generated rather than listed: all triangle counts up to three at each of the
/// four corners, crossed with no quad and with each of the three quads up to
/// three deep. That is 4⁴ × 10 = 2,560 configurations, and each one is a
/// round trip through `e = M·n` and back.
#[test]
fn decompose_round_trips_on_every_normal_surface() {
    let mut checked = 0usize;
    let mut with_quads = 0usize;
    for t0 in 0..4u32 {
        for t1 in 0..4u32 {
            for t2 in 0..4u32 {
                for t3 in 0..4u32 {
                    for (pair, depth) in [(usize::MAX, 0u32)]
                        .into_iter()
                        .chain((0..3).flat_map(|p| (1..4u32).map(move |d| (p, d))))
                    {
                        let mut quad = [0u32; 3];
                        if pair != usize::MAX {
                            quad[pair] = depth;
                        }
                        let n = NormalSurface {
                            triangle: [t0, t1, t2, t3],
                            quad,
                        };
                        let e = n.edge_coordinates();

                        // A normal surface always meets the boundary in normal
                        // curves, so this is also a check on `validate`.
                        assert!(e.is_normal(), "{n:?} gave a non-normal boundary: {e:?}");

                        let back = decompose(&e).unwrap_or_else(|| {
                            panic!("{n:?} -> {e:?} did not decompose");
                        });
                        assert_eq!(back, n, "round trip changed the surface");
                        checked += 1;
                        if pair != usize::MAX {
                            with_quads += 1;
                        }
                    }
                }
            }
        }
    }
    std::println!("{checked} normal surfaces round-tripped, {with_quads} of them with quads");
    assert_eq!(checked, 4 * 4 * 4 * 4 * 10);
}

/// **The paper's own counterexample**, §2.3:
///
/// > Consider for example the edge coordinates `e = (2,1,1,2,1,1)`. Here the
/// > solution to Equation 3 is `n = (0,0,0,0,0,1,1)`, decomposing `e` into
/// > **intersecting** quads.
///
/// Two non-zero quads violate *"only one of the coordinates `qᵢⱼ` can be
/// nonzero"*, so this must not decompose — and it must still be a perfectly valid
/// normal *curve* system on the boundary, which is precisely why the general
/// reconstruction in A-014b is needed rather than a lookup.
#[test]
fn the_papers_intersecting_quad_example_has_no_normal_decomposition() {
    let e = from_paper_order([2, 1, 1, 2, 1, 1]);
    assert!(
        e.is_normal(),
        "the boundary curves are normal even though the surface is not: {:?}",
        e.validate()
    );
    assert_eq!(
        decompose(&e),
        None,
        "two intersecting quads were accepted as a normal surface"
    );

    // And the two quads it *would* be are indeed both non-zero, which is the
    // reason: q02 and q03 in the paper's naming.
    let q02 = NormalSurface {
        triangle: [0; 4],
        quad: {
            let mut q = [0u32; 3];
            q[edge_between(0, 2) as usize] = 1;
            q
        },
    };
    let q03 = NormalSurface {
        triangle: [0; 4],
        quad: {
            let mut q = [0u32; 3];
            q[edge_between(0, 3) as usize] = 1;
            q
        },
    };
    let (a, b) = (q02.edge_coordinates(), q03.edge_coordinates());
    let mut sum = [0u32; TET_EDGE_COUNT];
    for (edge, slot) in sum.iter_mut().enumerate() {
        *slot = a.count[edge] + b.count[edge];
    }
    assert_eq!(
        EdgeCoordinates { count: sum },
        e,
        "the paper's decomposition should reproduce its own example"
    );
}

/// The paper's second example, §2.4:
///
/// > For instance, edge coordinates `e = (1,3,0,3,3,0)` (pictured in Figure 5)
/// > **cannot be encoded by any normal or almost normal surface**.
///
/// Here the boundary curves are not even normal, which is the harder case the
/// general procedure exists for: *"we decompose arbitrary curves on the boundary
/// of each tet into normal and non-normal parts… rather than restrict ourselves
/// to normal curves."*
#[test]
fn the_papers_non_normal_example_is_neither_normal_nor_decomposable() {
    let e = from_paper_order([1, 3, 0, 3, 3, 0]);
    assert!(!e.is_normal(), "expected a non-normal boundary");
    assert!(
        matches!(e.validate(), Err(NotNormal::TriangleInequality { .. })),
        "{:?}",
        e.validate()
    );
    assert_eq!(decompose(&e), None);
}

// ─── the cross-check against something already trusted ──────────────────────

/// **Classic Marching Tetrahedra is the 0/1 case of this encoding**, checked
/// against A-003's own table on all sixteen configurations of all six tets.
///
/// The paper's claim is that marching tetrahedra *"reinvented a small piece of
/// this story"*. If that is true then for any sign configuration, taking `eᵢⱼ = 1`
/// exactly where the signs differ must decompose into normal polygons whose
/// triangle count matches what A-003 emits — counting a quad as the two triangles
/// it is cut into.
///
/// This is the test that connects new machinery to verified machinery. It would
/// fail on a wrong complementary pairing, a wrong incidence vector, or a wrong
/// corner-coordinate formula.
#[test]
fn classic_marching_tetrahedra_is_the_zero_one_case_of_this_encoding() {
    let mut corner_cuts = 0usize;
    let mut diagonal_cuts = 0usize;

    for (t, cases) in TET_CASES.iter().enumerate() {
        for (case, expected) in cases.iter().enumerate() {
            let inside = |corner: u8| case & (1 << corner) != 0;

            let mut count = [0u32; TET_EDGE_COUNT];
            for (edge, [a, b]) in TET_EDGES.iter().copied().enumerate() {
                if inside(a) != inside(b) {
                    count[edge] = 1;
                }
            }
            let e = EdgeCoordinates { count };
            assert!(e.is_classic(), "a sign test produced a count above one");
            assert!(e.is_normal(), "case {case:#06b} has a non-normal boundary");

            let n = decompose(&e)
                .unwrap_or_else(|| panic!("tet {t} case {case:#06b} -> {e:?} did not decompose"));

            // A quad becomes two triangles; a corner triangle stays one.
            let triangles = n.triangle.iter().sum::<u32>() + 2 * n.quad.iter().sum::<u32>();
            assert_eq!(
                triangles,
                u32::from(expected.count),
                "tet {t} case {case:#06b}: encoding says {triangles} triangles, A-003 says {}",
                expected.count
            );

            match (n.triangle.iter().sum::<u32>(), n.quad.iter().sum::<u32>()) {
                (0, 0) => {}
                (1, 0) => corner_cuts += 1,
                (0, 1) => diagonal_cuts += 1,
                other => panic!("tet {t} case {case:#06b} gave {other:?}, not a single cut"),
            }
        }
    }

    // Per tet: 2 trivial cases, 8 corner cuts, 6 diagonal cuts — P-4's arithmetic,
    // arrived at from the opposite direction.
    std::println!("{corner_cuts} corner cuts, {diagonal_cuts} diagonal cuts over six tets");
    assert_eq!(corner_cuts, 6 * 8);
    assert_eq!(diagonal_cuts, 6 * 6);
}

/// What a sign-based method cannot see, stated as a measurement.
///
/// Any edge carrying two or more crossings is invisible to Marching Tetrahedra —
/// it reads the parity and sees none. Over every configuration with counts up to
/// three, this counts how many are indistinguishable from a *classic* one, which
/// is the size of the gap A-014 exists to close.
#[test]
fn the_share_of_configurations_a_sign_test_cannot_distinguish_is_reported() {
    let mut normal = 0usize;
    let mut classic = 0usize;
    let mut aliased = 0usize;

    // Every count in 0..4 on all six edges.
    for raw in 0..4096u32 {
        let mut count = [0u32; TET_EDGE_COUNT];
        for (edge, slot) in count.iter_mut().enumerate() {
            *slot = (raw >> (2 * edge)) & 0b11;
        }
        let e = EdgeCoordinates { count };
        if !e.is_normal() {
            continue;
        }
        normal += 1;
        if e.is_classic() {
            classic += 1;
        } else {
            // A sign test sees only parity, so this configuration is read as the
            // classic one with the same parities.
            aliased += 1;
        }
    }

    std::println!(
        "{normal} normal configurations with counts up to 3: {classic} classic, \
         {aliased} aliased away by a sign test ({:.1}%)",
        100.0 * aliased as f64 / normal as f64
    );
    assert!(
        aliased > classic,
        "the aliased majority is the whole premise"
    );
}
