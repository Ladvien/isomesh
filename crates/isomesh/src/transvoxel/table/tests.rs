//! The derivation checked against the dissertation's own count, and against
//! itself.

use alloc::vec::Vec;

use super::*;
use crate::marching_cubes::table::{NO_EDGE, segment_links};

/// Every canonical `(case, joined)` pair — the ones a real cell can present.
fn every_configuration() -> impl Iterator<Item = (u16, u16)> {
    (0..512u16).flat_map(|case| {
        let ambiguous = AMBIGUOUS_FACES[case as usize];
        // Only bits on genuinely ambiguous faces mean anything, so enumerate the
        // subsets of `ambiguous` rather than all 2^9 masks.
        let mut masks = Vec::new();
        let mut sub = ambiguous;
        loop {
            masks.push(sub);
            if sub == 0 {
                break;
            }
            sub = (sub - 1) & ambiguous;
        }
        masks.into_iter().map(move |joined| (case, joined))
    })
}

// ─── the geometry, before any contour is drawn on it ────────────────────────

/// The cell's boundary must be a closed oriented surface, or the links cannot
/// form consistently wound cycles — and nothing downstream would notice.
///
/// Each of the sixteen edges must appear in exactly two faces, traversed in
/// opposite directions. A single transposition in [`FACES`] breaks this and
/// nothing else here would.
#[test]
fn every_edge_is_shared_by_two_faces_in_opposite_directions() {
    // Directed edge (from sample, to sample) -> how many faces traverse it.
    let mut directed: Vec<(u8, u8, u8)> = Vec::new();
    for face in &FACES {
        for k in 0..face.corners {
            let edge = face.edge[k];
            assert_ne!(edge, NO_EDGE, "a walked corner has no edge");
            let from = face.sample[k];
            let to = face.sample[(k + 1) % face.corners];
            directed.push((edge, from, to));
        }
    }
    assert_eq!(directed.len(), 4 * 4 + 4 + 3 * 4, "boundary corner count");

    for edge in 0..EDGE_COUNT as u8 {
        let uses: Vec<&(u8, u8, u8)> = directed.iter().filter(|d| d.0 == edge).collect();
        assert_eq!(uses.len(), 2, "edge {edge} is used {} times", uses.len());
        let (_, a_from, a_to) = *uses[0];
        let (_, b_from, b_to) = *uses[1];
        assert_eq!(
            (a_from, a_to),
            (b_to, b_from),
            "edge {edge} runs {a_from}->{a_to} and {b_from}->{b_to}; it must reverse"
        );
        // And its endpoints must be the samples the edge table names.
        let mut named = EDGE_SAMPLES[edge as usize];
        named.sort_unstable();
        let mut walked = [a_from, a_to];
        walked.sort_unstable();
        assert_eq!(named, walked, "edge {edge} joins the wrong samples");
    }
}

/// The four laterals are triangles, and a triangle has no pairing to choose.
///
/// This is why the transition itself needs no ambiguity rule: the only faces
/// that can be ambiguous are the four quadrants and the half-resolution face,
/// which is exactly what the dissertation names.
#[test]
fn only_quadrants_and_the_back_face_can_be_ambiguous() {
    for (f, face) in FACES.iter().enumerate() {
        if face.corners == 3 {
            for case in 0..512u16 {
                assert_eq!(
                    AMBIGUOUS_FACES[case as usize] & (1 << f),
                    0,
                    "lateral face {f} ambiguous at case {case:#011b}"
                );
            }
        }
    }
    // And at least one case does make a quadrant ambiguous, or the check above is
    // vacuous.
    assert!(
        (0..512u16).any(|c| AMBIGUOUS_FACES[c as usize] != 0),
        "no case is ambiguous anywhere"
    );
}

/// Exactly four edges carry the half resolution, and they are the four that join
/// the corner samples.
#[test]
fn the_half_resolution_edges_are_the_corner_pairs() {
    let half: Vec<u8> = (0..EDGE_COUNT as u8)
        .filter(|e| is_half_resolution(*e))
        .collect();
    assert_eq!(half, alloc::vec![12, 13, 14, 15]);
    for e in half {
        let [a, b] = EDGE_SAMPLES[e as usize];
        assert!(
            matches!(a, 0 | 2 | 6 | 8) && matches!(b, 0 | 2 | 6 | 8),
            "half-resolution edge {e} joins {a} and {b}, which are not both corners"
        );
    }
}

// ─── the contour ────────────────────────────────────────────────────────────

/// The links must cover exactly the cut edges, and form closed cycles.
///
/// Every cut edge is an entry on exactly one of its two faces and an exit on the
/// other, so it has exactly one successor and exactly one predecessor. That is
/// what makes "one vertex per cycle" well defined, and it is the property A-010
/// needed on the cube.
#[test]
fn the_links_cover_the_cut_edges_and_close() {
    for (case, joined) in every_configuration() {
        let next = transition_links(case, joined);
        let cut = cut_edges(case);

        let mut linked = 0u16;
        let mut targeted = 0u16;
        for (e, &t) in next.iter().enumerate() {
            if t != NO_EDGE {
                linked |= 1 << e;
                assert_eq!(
                    targeted & (1 << t),
                    0,
                    "case {case:#011b}/{joined:#b}: edge {t} is the target of two segments"
                );
                targeted |= 1 << t;
            }
        }
        assert_eq!(
            linked, cut,
            "case {case:#011b}/{joined:#b}: linked {linked:#018b}, cut {cut:#018b}"
        );
        assert_eq!(
            targeted, cut,
            "case {case:#011b}/{joined:#b}: targeted {targeted:#018b}, cut {cut:#018b}"
        );
    }
}

/// Walking from any cut edge returns to it, and the cycles partition the cut
/// edges.
#[test]
fn the_cycles_partition_the_cut_edges() {
    let mut longest = 0usize;
    let mut most_cycles = 0usize;
    for (case, joined) in every_configuration() {
        let next = transition_links(case, joined);
        let mut visited = 0u16;
        let mut cycles = 0usize;
        for start in 0..EDGE_COUNT as u8 {
            if next[start as usize] == NO_EDGE || visited & (1 << start) != 0 {
                continue;
            }
            cycles += 1;
            let mut len = 0usize;
            let mut current = start;
            while visited & (1 << current) == 0 {
                visited |= 1 << current;
                len += 1;
                current = next[current as usize];
            }
            assert_eq!(
                current, start,
                "case {case:#011b}/{joined:#b}: a walk closed on {current}, not {start}"
            );
            assert!(len >= 3, "case {case:#011b}: a cycle of length {len}");
            longest = longest.max(len);
        }
        most_cycles = most_cycles.max(cycles);
        assert_eq!(
            visited,
            cut_edges(case),
            "case {case:#011b}: cycles missed an edge"
        );
    }
    std::println!("longest cycle {longest}, most cycles in one cell {most_cycles}");
}

/// The full-resolution face's contour must be exactly what Marching Cubes draws
/// on those squares, or the seam does not close on the fine side.
///
/// Cross-checked against [`segment_links`] rather than against a restatement of
/// it: each quadrant is embedded in a cube case whose near-`z` face carries the
/// same four signs, and the pairing must agree. The two walks were written
/// separately — one over cube faces, one over a nine-face cell — so agreeing on
/// all sixteen sign patterns of all four quadrants is a real check and not a
/// tautology.
///
/// A link is *this face's* when both its ends are edges of this face. Every edge
/// has exactly one successor, set by whichever of its two faces sees it as an
/// entry, and a face shares only one edge with any neighbour — so that rule picks
/// out exactly the face's own pairs on both sides of the comparison.
#[test]
fn each_quadrant_pairs_as_marching_cubes_would() {
    // `cube::face_corners(2, 0)` — the near-z face, counter-clockwise from
    // outside, which is the same winding the quadrants use.
    let cube_face = crate::cube::face_corners(2, 0);
    let cube_edges = [
        crate::cube::edge_index(cube_face[0], cube_face[1]),
        crate::cube::edge_index(cube_face[1], cube_face[2]),
        crate::cube::edge_index(cube_face[2], cube_face[3]),
        crate::cube::edge_index(cube_face[3], cube_face[0]),
    ];

    let mut compared = 0usize;
    for (q, face) in FACES.iter().enumerate().take(4) {
        for signs in 0..16u8 {
            let mut case = 0u16;
            let mut cube_case = 0u8;
            for (k, corner) in cube_face.iter().enumerate() {
                if signs & (1 << k) != 0 {
                    case |= 1 << face.sample[k];
                    cube_case |= 1 << corner;
                }
            }

            let t = transition_links(case, 0);
            let c = segment_links(cube_case, 0);

            // Both sides, reduced to "position k pairs to position m".
            let own = |links: &[u8], edges: &[u8; 4]| -> Vec<(usize, usize)> {
                let mut pairs = Vec::new();
                for (k, e) in edges.iter().enumerate() {
                    let to = links[*e as usize];
                    if let Some(m) = edges.iter().position(|x| *x == to) {
                        pairs.push((k, m));
                    }
                }
                pairs.sort_unstable();
                pairs
            };
            let quadrant_edges = [face.edge[0], face.edge[1], face.edge[2], face.edge[3]];

            assert_eq!(
                own(&t, &quadrant_edges),
                own(&c, &cube_edges),
                "quadrant {q}, signs {signs:#06b}"
            );
            compared += 1;
        }
    }
    assert_eq!(compared, 64);
}

/// A lateral face is where the resolution changes — and where a feature too fine
/// for the coarse grid gets capped off.
///
/// Written first as *"a lateral link always joins a fine edge to a coarse one"*,
/// which is false, and the counterexample is the interesting half. With only the
/// midpoint sample inside, the two fine sub-edges are both cut and the coarse
/// edge is not, so the lateral links **fine to fine** and closes the bump off
/// entirely on the fine side. The coarse neighbour sees both endpoints outside
/// and contributes nothing, which is exactly right: it cannot represent that
/// feature, so the seam must not ask it to.
///
/// The true rule, asserted here, is sharper and follows from the signs: on a
/// lateral with samples `(a, m, b)` the coarse edge is cut precisely when
/// `sign(a) != sign(b)`, and then exactly one fine sub-edge is cut too. So
///
/// - coarse edge cut → the link joins fine to coarse, and the seam is stitched;
/// - coarse edge not cut → the link joins fine to fine, and a sub-coarse feature
///   is capped.
///
/// Both arms are exercised, and the counts are reported because the second is the
/// one nobody expects.
#[test]
fn a_lateral_link_stitches_the_seam_or_caps_a_feature_the_coarse_grid_cannot_see() {
    let mut stitched = 0usize;
    let mut capped = 0usize;
    for (case, joined) in every_configuration() {
        let next = transition_links(case, joined);
        for face in FACES.iter().filter(|f| f.corners == 3) {
            let edges = [face.edge[0], face.edge[1], face.edge[2]];
            // Read the coarse edge off the face rather than off corner
            // positions: the laterals do not all put it in the same place in
            // their walk, which is what the first version of this test got wrong.
            let coarse = edges
                .iter()
                .copied()
                .find(|e| is_half_resolution(*e))
                .expect("every lateral has exactly one half-resolution edge");
            assert_eq!(
                edges.iter().filter(|e| is_half_resolution(**e)).count(),
                1,
                "a lateral must have exactly one half-resolution edge"
            );
            let coarse_cut = cut_edges(case) & (1 << coarse) != 0;

            for e in edges {
                let to = next[e as usize];
                if to == NO_EDGE || !edges.contains(&to) {
                    continue;
                }
                let crosses = is_half_resolution(e) != is_half_resolution(to);
                assert_eq!(
                    crosses,
                    coarse_cut,
                    "case {case:#011b}/{joined:#b}: lateral {:?} links {e} to {to}, \
                     coarse edge {coarse} cut = {coarse_cut}",
                    &face.sample[..3]
                );
                if crosses {
                    stitched += 1;
                } else {
                    capped += 1;
                }
            }
        }
    }
    std::println!("{stitched} lateral links stitch the seam, {capped} cap a sub-coarse feature");
    assert!(stitched > 0, "no lateral ever stitched the seam");
    assert!(
        capped > 0,
        "no lateral ever capped a feature — the counterexample is gone"
    );
}

// ─── the dissertation's number ──────────────────────────────────────────────

/// `D₈` must actually be a group of order eight before it can classify anything.
#[test]
fn the_dihedral_group_is_closed_and_has_order_eight() {
    let mut seen: Vec<Perm> = Vec::new();
    for p in &DIHEDRAL {
        assert!(!seen.contains(p), "D8 has a repeated element");
        seen.push(*p);
    }
    assert_eq!(seen.len(), 8);

    // Every element is a permutation of the nine samples.
    for p in &DIHEDRAL {
        let mut sorted = p.to_vec();
        sorted.sort_unstable();
        assert_eq!(sorted, (0..SAMPLE_COUNT as u8).collect::<Vec<_>>());
    }

    // Closed under composition.
    for a in &DIHEDRAL {
        for b in &DIHEDRAL {
            let mut c = [0u8; SAMPLE_COUNT];
            for s in 0..SAMPLE_COUNT {
                c[s] = b[a[s] as usize];
            }
            assert!(DIHEDRAL.contains(&c), "D8 is not closed under composition");
        }
    }

    // The centre sample is fixed by everything, and the corners map to corners.
    for p in &DIHEDRAL {
        assert_eq!(p[4], 4, "the centre must be fixed");
        for corner in [0u8, 2, 6, 8] {
            assert!(matches!(p[corner as usize], 0 | 2 | 6 | 8));
        }
    }
}

/// **The acceptance criterion for the derivation.** Lengyel 2010 §4.3:
///
/// > Observing the orbits of each of the 512 cases under the action of the group
/// > `D₈` or `D₈ × ℤ₂`, as appropriate, yields exactly **73 distinct equivalence
/// > classes**.
///
/// The group is `D₈ × ℤ₂` — rotations, reflections and inversion — for a case
/// with no ambiguous face, and `D₈` alone otherwise, because inverting an
/// ambiguous case does not give a geometrically equivalent one. Both halves of
/// that rule are the dissertation's.
///
/// Reproducing 73 from geometry defined independently of the paper's numbering is
/// the strongest available evidence that this cell is that cell.
#[test]
fn the_orbit_count_is_lengyels_seventy_three() {
    let canonical = |case: u16| -> u16 {
        let inversion_applies = AMBIGUOUS_FACES[case as usize] == 0;
        let mut best = case;
        for p in &DIHEDRAL {
            let image = permute(case, p);
            best = best.min(image);
            if inversion_applies {
                best = best.min(invert(image));
            }
        }
        best
    };

    // Ambiguity must be a property of the class, not of the representative, or
    // "as appropriate" is not well defined.
    for case in 0..512u16 {
        let ambiguous = AMBIGUOUS_FACES[case as usize] != 0;
        for p in &DIHEDRAL {
            assert_eq!(
                AMBIGUOUS_FACES[permute(case, p) as usize] != 0,
                ambiguous,
                "case {case:#011b} changes ambiguity under a symmetry"
            );
        }
        assert_eq!(
            AMBIGUOUS_FACES[invert(case) as usize] != 0,
            ambiguous,
            "case {case:#011b} changes ambiguity under inversion"
        );
    }

    let mut classes: Vec<u16> = (0..512u16).map(canonical).collect();
    classes.sort_unstable();
    classes.dedup();

    std::println!("{} equivalence classes over 512 cases", classes.len());
    assert_eq!(
        classes.len(),
        73,
        "Lengyel 2010 section 4.3 says exactly 73; got {}",
        classes.len()
    );

    // And the trivial class holds both all-outside and all-inside, which is the
    // dissertation's own description of it.
    assert_eq!(
        canonical(0),
        canonical(0x1ff),
        "the trivial class must be one class"
    );
}

/// Every orbit size divides eight, which is the dissertation's stated corollary
/// of the class size being the index of a stabiliser subgroup in `D₈`.
#[test]
fn every_dihedral_orbit_size_divides_eight() {
    for case in 0..512u16 {
        let mut orbit: Vec<u16> = DIHEDRAL.iter().map(|p| permute(case, p)).collect();
        orbit.sort_unstable();
        orbit.dedup();
        assert_eq!(
            8 % orbit.len(),
            0,
            "case {case:#011b} has orbit size {}",
            orbit.len()
        );
    }
}

/// The 51/22 split the dissertation reports for Table 4.8.
///
/// > Table 4.8 lists the 51 transition cell equivalence classes having four or
/// > fewer sample values classified as inside solid space. For equivalence
/// > classes that do not include geometric inverses due to the presence of an
/// > ambiguous face, the inverse equivalence class is listed, and this accounts
/// > for the remaining 22 classes.
///
/// A second, independent handle on the same classification: 51 + 22 = 73 has to
/// come out of the *same* partition, and a derivation that hit 73 by accident
/// would be unlikely to split it this way as well.
#[test]
fn the_class_split_is_fifty_one_and_twenty_two() {
    let canonical = |case: u16| -> u16 {
        let inversion_applies = AMBIGUOUS_FACES[case as usize] == 0;
        let mut best = case;
        for p in &DIHEDRAL {
            let image = permute(case, p);
            best = best.min(image);
            if inversion_applies {
                best = best.min(invert(image));
            }
        }
        best
    };

    let mut classes: Vec<u16> = (0..512u16).map(canonical).collect();
    classes.sort_unstable();
    classes.dedup();

    // A class "has four or fewer inside" if some member does.
    let mut light = 0usize;
    let mut heavy = 0usize;
    for &rep in &classes {
        let members: Vec<u16> = (0..512u16).filter(|c| canonical(*c) == rep).collect();
        if members.iter().any(|m| m.count_ones() <= 4) {
            light += 1;
        } else {
            heavy += 1;
        }
    }
    std::println!("{light} classes reach four-or-fewer inside, {heavy} do not");
    assert_eq!(light, 51);
    assert_eq!(heavy, 22);
}
