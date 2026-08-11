//! The 256-case Marching Cubes table, **derived at compile time**.
//!
//! # Why this is generated rather than transcribed
//!
//! The implementation brief is explicit that the table must not be typed from
//! memory, because a wrong entry produces a mesh that looks fine and is silently
//! non-manifold. It offers two defences: take the table from a cited source, and
//! write a validator that checks all 256 cases without a reference table.
//!
//! The papers are in the local corpus — Lewiner et al., *Efficient
//! Implementation of Marching Cubes' Cases with Topological Guarantees*, Journal
//! of Graphics Tools 8:2 (2003), `10.1080/10867651.2003.10487582`, and Lengyel's
//! 2010 Transvoxel dissertation — but **their case tables did not survive
//! PDF-to-text conversion**: Lengyel's Tables 3.1 and 3.2 are figures of cube
//! diagrams, and Lewiner's Table 1 converts to scrambled cells. Reading
//! triangulations off diagrams is the guessing the rules forbid.
//!
//! So the table is *constructed* instead, which removes the transcription
//! entirely, and then checked — see [`super::validate_table`]. What is taken
//! from Lengyel (dissertation §3.1.1, p.20) is the structure and the
//! conventions, not the numbers:
//!
//! - the 8-bit case index is the concatenated inside/outside corner bits;
//! - *"a positive value means that a sample point is outside the terrain in
//!   empty space, and a negative value means that it is inside it in solid
//!   space"* — which is this crate's sign convention;
//! - *"A choice must be made globally as to whether a sample value of exactly
//!   zero is considered to be in empty space or in solid space, and a consistent
//!   classification either way is acceptable."* See [`is_inside`].
//!
//! # The construction
//!
//! For one corner-sign configuration:
//!
//! 1. A cube edge is **cut** when its two corners are classified differently.
//!    Surface vertices live on cut edges, one each.
//! 2. On each face, walk the four corners counter-clockwise *as seen from
//!    outside the cube*. Each maximal run of inside corners is entered at one
//!    cut edge and left at another, giving a directed segment across that face.
//! 3. **Every cut edge gets exactly one outgoing and one incoming segment.** A
//!    cut edge lies on exactly two faces, and two faces of a cube with outward
//!    normals induce *opposite* directions on their shared edge — so it is an
//!    entry on one and an exit on the other. The segments therefore link the cut
//!    edges into directed cycles with no choices left over.
//! 4. Fan-triangulate each cycle.
//!
//! Nothing in that is a lookup, so there is nothing to mistype. Step 3 is also
//! what makes the result crack-free: a face's segments are a function of that
//! face's own four corner signs, so two cells sharing a face necessarily agree
//! on it. `face_segments_depend_only_on_that_face` checks exactly that.
//!
//! # The one real choice
//!
//! Step 2 resolves an ambiguous face — one with two diagonally opposite inside
//! corners — by **separating the inside corners**, since each run of inside
//! corners is cut off by its own segment. Both pairings are crack-free as long
//! as the rule is applied consistently, which is why this is a decision rather
//! than a fact. Marching Cubes proper has no way to do better; MC33's asymptotic
//! decider replaces this rule with one that reads the bilinear saddle, and that
//! is A-002's ticket.
//!
//! Interior ambiguity is likewise left alone: where a cell could host either two
//! separate sheets or a tunnel, this construction produces separate sheets.
//! Custodio et al. (2013) put the practical cost in perspective — the vast
//! majority of real-world cells match the unambiguous configurations, and a game
//! needs topological *consistency*, which this has, rather than topological
//! *correctness*, which needs the interior test.

// Cube topology and the sign rule live in `crate::cube`, shared with every
// other extractor. Re-exported here because this module's docs and the table
// construction below are written in terms of them.
pub use crate::cube::{
    CORNER_COUNT, EDGE_AXIS, EDGE_CORNERS, EDGE_COUNT, corner_inside, edge_index, face_corners,
    is_inside,
};

/// Upper bound on triangles per cell.
///
/// The construction cannot exceed this; `no_case_exceeds_the_triangle_bound`
/// records the maximum actually reached, which is five.
pub const MAX_TRIANGLES: usize = 12;

/// Marker for "no edge".
pub const NO_EDGE: u8 = u8::MAX;

/// The triangulation for one corner-sign configuration.
#[derive(Clone, Copy, Debug)]
pub struct McCase {
    /// Number of triangles.
    pub count: u8,
    /// Each triangle as three edge indices. Only the first `count` are valid.
    pub triangles: [[u8; 3]; MAX_TRIANGLES],
}

/// All 256 configurations, built during compilation.
pub static CASES: [McCase; 256] = build_cases();

const fn build_cases() -> [McCase; 256] {
    let mut out = [McCase {
        count: 0,
        triangles: [[0u8; 3]; MAX_TRIANGLES],
    }; 256];
    let mut case = 0usize;
    while case < 256 {
        out[case] = build_case(case as u8);
        case += 1;
    }
    out
}

/// Segment links for one configuration: `next[e]` is the cut edge the segment
/// leaving cut edge `e` arrives at, or [`NO_EDGE`].
///
/// Exposed so the validator can rebuild them independently of the triangles.
pub const fn segment_links(case: u8) -> [u8; EDGE_COUNT] {
    let mut next = [NO_EDGE; EDGE_COUNT];

    let mut axis = 0usize;
    while axis < 3 {
        let mut side = 0u8;
        while side < 2 {
            let c = face_corners(axis, side);

            // Start the walk on an outside corner so the first transition seen
            // is an entry. If every corner is inside, the face is not cut.
            let mut start = 4usize;
            let mut k = 0usize;
            while k < 4 {
                if !corner_inside(case, c[k]) {
                    start = k;
                    break;
                }
                k += 1;
            }

            if start < 4 {
                let mut entry = NO_EDGE;
                let mut j = 0usize;
                while j < 4 {
                    let p = c[(start + j) % 4];
                    let q = c[(start + j + 1) % 4];
                    let p_in = corner_inside(case, p);
                    let q_in = corner_inside(case, q);
                    if !p_in && q_in {
                        entry = edge_index(p, q);
                    } else if p_in && !q_in {
                        next[entry as usize] = edge_index(p, q);
                    }
                    j += 1;
                }
            }

            side += 1;
        }
        axis += 1;
    }

    next
}

const fn build_case(case: u8) -> McCase {
    let next = segment_links(case);

    let mut triangles = [[0u8; 3]; MAX_TRIANGLES];
    let mut count = 0usize;
    let mut visited = [false; EDGE_COUNT];

    let mut e = 0usize;
    while e < EDGE_COUNT {
        if next[e] != NO_EDGE && !visited[e] {
            // Walk the cycle this cut edge belongs to.
            let mut cycle = [0u8; EDGE_COUNT];
            let mut len = 0usize;
            let mut cur = e as u8;
            while !visited[cur as usize] {
                visited[cur as usize] = true;
                cycle[len] = cur;
                len += 1;
                cur = next[cur as usize];
            }

            // Fan from the first vertex. Every cycle has at least three edges:
            // two distinct cube edges share at most one face, so no two-cycle
            // can form.
            let mut i = 1usize;
            while i + 1 < len {
                triangles[count] = [cycle[0], cycle[i], cycle[i + 1]];
                count += 1;
                i += 1;
            }
        }
        e += 1;
    }

    McCase {
        count: count as u8,
        triangles,
    }
}
