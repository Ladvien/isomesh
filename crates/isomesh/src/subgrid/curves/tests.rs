//! §3.1's own claim, checked: the same face gives the same curve from either tet.

use alloc::vec::Vec;

use super::*;
use crate::marching_tetrahedra::table::TET_EDGE_COUNT;
use crate::subgrid::coordinates::TET_FACES;

/// Every edge-coordinate vector with counts up to `max`.
fn every_coordinate(max: u32) -> impl Iterator<Item = EdgeCoordinates> {
    let span = max + 1;
    let total = span.pow(TET_EDGE_COUNT as u32);
    (0..total).map(move |raw| {
        let mut count = [0u32; TET_EDGE_COUNT];
        let mut n = raw;
        for slot in &mut count {
            *slot = n % span;
            n /= span;
        }
        EdgeCoordinates::new(count)
    })
}

// ─── the conformity property ────────────────────────────────────────────────

/// **§3.1's central claim**: *"We carefully define this procedure to produce
/// identical curves on triangles shared by neighboring tetrahedra."*
///
/// Two tets sharing a face differ in exactly one thing: their fourth vertex, and
/// the three edges running to it. So the claim is precisely that a face's
/// segments are a function of **that face's own three edge coordinates and
/// nothing else** — vary the other three however you like and the answer must not
/// move. That is what lets each tet run alone and still agree.
///
/// # This was first written as invariance under relabelling, and that is false
///
/// The first version permuted the face's corner labels and demanded the same
/// segments. It fails, correctly, on `e = (3, 0, 0)`: step 3(b) skips *"the first
/// residual point along oriented edge `ij` (assuming a canonical orientation
/// `i < j`)"*, so swapping `i` and `j` skips the point at the other end.
///
/// The construction **depends on vertex order by design**, and that is safe for
/// exactly the reason the relabelling test ignored: in a real tet mesh the two
/// tets share the face's *global vertex indices*, so they agree on `i < j`
/// without communicating. Conformity here is locality plus a shared ordering —
/// not symmetry. See M-79.
#[test]
fn a_shared_face_gives_the_same_segments_from_either_tet() {
    // Face 3 is corners (0, 1, 2); the edges to corner 3 are the ones a
    // neighbouring tet would differ on.
    let off_face: Vec<u8> = (0..TET_EDGE_COUNT as u8)
        .filter(|e| {
            let [lo, hi] = crate::marching_tetrahedra::table::TET_EDGES[*e as usize];
            lo == 3 || hi == 3
        })
        .collect();
    assert_eq!(off_face.len(), 3);

    let mut checked = 0usize;
    for coords in every_coordinate(3) {
        let base = face_segments(3, &coords);

        // Every way the other tet could have filled in its own three edges.
        for other in 0..4u32.pow(3) {
            let mut count = coords.count;
            let mut n = other;
            for e in &off_face {
                count[*e as usize] = n % 4;
                n /= 4;
            }
            assert_eq!(
                face_segments(3, &EdgeCoordinates::new(count)),
                base,
                "the face moved when a neighbour's own edges changed: {:?} -> {count:?}",
                coords.count
            );
        }
        checked += 1;
    }
    std::println!("{checked} faces, each unchanged by all 64 neighbour configurations");
    assert!(checked > 4000);
}

// ─── the three steps ────────────────────────────────────────────────────────

/// Step 1 on a normal face: every crossing is used exactly once, and every
/// segment joins two *different* edges.
#[test]
fn a_normal_face_pairs_every_crossing_across_two_edges() {
    let mut checked = 0usize;
    for coords in every_coordinate(4) {
        let f = TET_FACES[3];
        let e: [u32; 3] = [
            coords.edge(f.edge[0]),
            coords.edge(f.edge[1]),
            coords.edge(f.edge[2]),
        ];
        let even = (e[0] + e[1] + e[2]).is_multiple_of(2);
        let triangle = (0..3).all(|k| e[k] + e[(k + 2) % 3] >= e[(k + 1) % 3]);
        if !(even && triangle) {
            continue;
        }

        let segments = face_segments(3, &coords);
        // Two endpoints per segment, and every crossing on the face used once.
        assert_eq!(
            segments.len() * 2,
            (e[0] + e[1] + e[2]) as usize,
            "{:?} left crossings unpaired",
            coords.count
        );
        assert!(
            segments.iter().all(|s| !s.is_scoop()),
            "{:?} produced a scoop on a normal face",
            coords.count
        );

        let mut used: Vec<FacePoint> = segments.iter().flat_map(|s| [s.a, s.b]).collect();
        used.sort_unstable();
        let before = used.len();
        used.dedup();
        assert_eq!(
            before,
            used.len(),
            "{:?} used a crossing twice",
            coords.count
        );
        checked += 1;
    }
    std::println!("{checked} normal faces paired cleanly");
    assert!(checked > 100);
}

/// **Step 2's forced choice.** An odd sum drops one from each edge coordinate,
/// *"effectively creating three open endpoints"* — and which point is left over
/// is not a choice: the corner cuts claim `cᵢ` from one end and `cⱼ` from the
/// other, and `cᵢ + cⱼ = eᵢⱼ − 1`, so exactly one is unclaimed per edge.
#[test]
fn the_reduced_pass_leaves_exactly_one_point_per_edge() {
    let mut checked = 0usize;
    for coords in every_coordinate(4) {
        let f = TET_FACES[3];
        let e: [u32; 3] = [
            coords.edge(f.edge[0]),
            coords.edge(f.edge[1]),
            coords.edge(f.edge[2]),
        ];
        let odd = !(e[0] + e[1] + e[2]).is_multiple_of(2);
        let triangle = (0..3).all(|k| e[k] + e[(k + 2) % 3] >= e[(k + 1) % 3]);
        if !(odd && triangle) {
            continue;
        }

        let segments = face_segments(3, &coords);
        for (k, edge) in f.edge.iter().enumerate() {
            let used = segments
                .iter()
                .flat_map(|s| [s.a, s.b])
                .filter(|p| p.edge == *edge)
                .count() as u32;
            assert_eq!(
                used,
                e[k].saturating_sub(1),
                "edge {edge} of {:?}: {used} used of {}",
                coords.count,
                e[k]
            );
        }
        checked += 1;
    }
    std::println!("{checked} odd-sum faces each left one endpoint open per edge");
    assert!(checked > 100);
}

/// **Step 3's scoops.** When one edge is longer than the other two together, the
/// excess is paired along that edge — and an odd excess skips one point.
#[test]
fn a_violated_inequality_produces_scoops_on_the_long_edge() {
    let mut even_case = 0usize;
    let mut odd_case = 0usize;
    for coords in every_coordinate(5) {
        let f = TET_FACES[3];
        let e: [u32; 3] = [
            coords.edge(f.edge[0]),
            coords.edge(f.edge[1]),
            coords.edge(f.edge[2]),
        ];
        let Some(k) = (0..3).find(|k| e[*k] > e[(*k + 1) % 3] + e[(*k + 2) % 3]) else {
            continue;
        };
        let r = e[k] - e[(k + 1) % 3] - e[(k + 2) % 3];

        let segments = face_segments(3, &coords);
        let scoops: Vec<&Segment> = segments.iter().filter(|s| s.is_scoop()).collect();
        assert!(
            scoops.iter().all(|s| s.a.edge == f.edge[k]),
            "a scoop landed on an edge that is not the long one: {:?}",
            coords.count
        );
        assert_eq!(
            scoops.len(),
            (r / 2) as usize,
            "{:?}: residual {r} gave {} scoops",
            coords.count,
            scoops.len()
        );
        // Scoops join adjacent crossings -- "least geometric length".
        assert!(
            scoops.iter().all(|s| s.b.index == s.a.index + 1),
            "a scoop skipped a crossing: {:?}",
            coords.count
        );

        if r.is_multiple_of(2) {
            even_case += 1;
        } else {
            odd_case += 1;
        }
    }
    std::println!("{even_case} even residuals, {odd_case} odd");
    assert!(even_case > 0 && odd_case > 0, "both arms must be exercised");
}

// ─── classification ─────────────────────────────────────────────────────────

/// The three kinds are all reachable, and a normal-surface configuration
/// produces only normal curves.
#[test]
fn the_three_curve_kinds_are_all_reachable() {
    let mut seen = [0usize; 3];
    for coords in every_coordinate(3) {
        for curve in curves(&coords) {
            match curve.kind {
                CurveKind::Open => seen[0] += 1,
                CurveKind::Normal => seen[1] += 1,
                CurveKind::NonNormal => seen[2] += 1,
            }
        }
    }
    std::println!(
        "curves seen: {} open, {} normal, {} non-normal",
        seen[0],
        seen[1],
        seen[2]
    );
    assert!(
        seen.iter().all(|n| *n > 0),
        "a kind was never produced: {seen:?}"
    );
}

/// A configuration that **is** a normal surface must give only normal curves —
/// that is the definition, and it ties §3.1 back to A-014a's `decompose`.
#[test]
fn a_normal_surface_gives_only_normal_curves() {
    let mut checked = 0usize;
    for coords in every_coordinate(3) {
        if super::super::coordinates::decompose(&coords).is_none() {
            continue;
        }
        for curve in curves(&coords) {
            assert_eq!(
                curve.kind,
                CurveKind::Normal,
                "{:?} decomposes into normal polygons but gave a {:?} curve",
                coords.count,
                curve.kind
            );
        }
        checked += 1;
    }
    std::println!("{checked} normal-surface configurations gave only normal curves");
    assert!(checked > 50);
}

/// Classic Marching Tetrahedra's configurations produce exactly one closed
/// normal curve each — the single disk per tet that A-014a showed is the 0/1
/// corner of the encoding.
#[test]
fn a_classic_configuration_gives_one_closed_curve() {
    let mut checked = 0usize;
    for case in 0..16u8 {
        let mut count = [0u32; TET_EDGE_COUNT];
        for (edge, [a, b]) in crate::marching_tetrahedra::table::TET_EDGES
            .iter()
            .copied()
            .enumerate()
        {
            let inside = |c: u8| case & (1 << c) != 0;
            if inside(a) != inside(b) {
                count[edge] = 1;
            }
        }
        let coords = EdgeCoordinates::new(count);
        let found = curves(&coords);
        if coords.total() == 0 {
            assert!(found.is_empty());
            continue;
        }
        assert_eq!(
            found.len(),
            1,
            "case {case:#06b} gave {} curves",
            found.len()
        );
        assert_eq!(found[0].kind, CurveKind::Normal, "case {case:#06b}");
        assert_eq!(
            found[0].segments.len(),
            coords.total() as usize,
            "a closed curve has as many segments as crossings"
        );
        checked += 1;
    }
    assert_eq!(checked, 14, "two of the sixteen cases are trivial");
}
