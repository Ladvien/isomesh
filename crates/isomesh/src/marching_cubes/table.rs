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
//! 4. Triangulate each cycle — see [`triangulate`] and [`safe_apex`] for which
//!    fan, and A-015 for why the choice of apex is not free.
//!
//! Nothing in that is a lookup, so there is nothing to mistype. Step 3 is also
//! what makes the result crack-free: a face's segments are a function of that
//! face's own four corner signs, so two cells sharing a face necessarily agree
//! on it. `face_segments_depend_only_on_that_face` checks exactly that.
//!
//! # The one real choice
//!
//! Step 2 resolves an ambiguous face — one with two diagonally opposite inside
//! corners — by pairing each entry with the exit that closes its own run, which
//! **separates the inside corners**. Both pairings are crack-free as long as the
//! rule is applied consistently, which is why this is a decision rather than a
//! fact.
//!
//! That decision is now a parameter. [`segment_links`] takes a six-bit mask, one
//! bit per face, selecting the crossed pairing — which joins the inside corners
//! instead — on the faces whose bit is set. [`CASES`] is the all-zero mask, so
//! Marching Cubes proper is unchanged; A-002's asymptotic decider computes the
//! mask per cell from the bilinear saddle, in [`super::ambiguity`]. Bits set on
//! a face that is not ambiguous have no effect, because such a face has at most
//! one entry and so no pairing to choose; `masks_are_ignored_on_unambiguous_faces`
//! is the proof, and it is what makes the table lookup a legitimate memo rather
//! than a second rule.
//!
//! Only a face with **four** cut edges has a choice. A face has 0, 2 or 4 cut
//! edges and nothing else, and four requires alternating signs around the ring —
//! [`AMBIGUOUS_FACES`] records exactly which faces those are, per case.
//!
//! # Interior ambiguity is not resolved here
//!
//! Where a cell could host either two separate sheets or a tunnel, this
//! construction produces separate sheets, under either mask. Resolving it needs
//! Chernyaev's body-saddle test, which is **deliberately not implemented** — see
//! A-002b in the backlog. Custodio et al. (2013) put the practical cost in
//! perspective: the vast majority of real-world cells match the unambiguous
//! configurations, and a game needs topological *consistency*, which this has,
//! rather than topological *correctness*, which needs the interior test.

// Cube topology and the sign rule live in `crate::cube`, shared with every
// other extractor. Re-exported here because this module's docs and the table
// construction below are written in terms of them.
pub use crate::cube::{
    CORNER_COUNT, EDGE_AXIS, EDGE_CORNERS, EDGE_COUNT, corner_inside, edge_index, edge_on_face,
    edges_share_a_face, face_corners, is_inside,
};

/// Upper bound on triangles per cell.
///
/// A centroid fan emits one triangle per cycle edge rather than `k − 2`, so a
/// cell with all twelve edges cut in one centroid-fanned cycle reaches twelve.
/// `no_case_exceeds_the_triangle_bound` and `the_decider_does_not_exceed_the_
/// triangle_bound` record the maxima actually reached.
pub const MAX_TRIANGLES: usize = 12;

/// Marker for "no edge".
pub const NO_EDGE: u8 = u8::MAX;

/// Most cycle centroids one cell can need.
///
/// A cell has at most twelve cut edges and a centroid is only ever created for a
/// cycle with no chord-safe apex, which measurement puts at length eight or more
/// — so one is the real ceiling. Three is kept as the arithmetic bound: four
/// edges is the shortest cycle that could in principle need one, and a bound
/// that does not depend on a measurement cannot go stale. See [`CENTROID_BASE`].
pub const MAX_CENTROIDS: usize = 3;

/// Triangle corner codes at or above this name a **cycle centroid**, not a cube
/// edge: `CENTROID_BASE + c` is the centroid of this cell's cycle `c`.
///
/// # Why a cell ever needs an interior vertex (A-015)
///
/// A polygon triangulated without extra vertices has `k − 3` interior chords,
/// and nothing *in general* stops two cells that share a face from choosing the
/// same chord — measured at A-002 as 12 of 4,096 two-cell sign patterns putting
/// **four** triangles on one mesh edge, identically under both ambiguity rules.
///
/// What makes a chord collidable is specific, though, and local: only a cell
/// containing **both** of its cut edges can emit it, and two cells share a pair
/// of cube edges only when those edges share a cube *face*. So a fan whose every
/// chord joins edges sharing no face is safe, and [`safe_apex`] looks for an apex
/// with that property. Measured over all 256 cases and every canonical mask,
/// **one exists for every cycle of length 3 to 7 and for 48 of the 60 length-8
/// cycles.** Plain Marching Cubes tops out at length 7, so it never reaches this
/// path at all and still places exactly one vertex per crossed grid edge —
/// ✗1/M-2/M-22's `V_mc = C` is intact.
///
/// A cycle with no safe apex fans from a centroid instead, which removes chords
/// entirely: every mesh edge is then either a cycle edge, lying on a face and so
/// belonging to exactly the two cells sharing it, or a spoke to a **cell-local**
/// centroid that no other cell can name. Two faces either way. Only the joined
/// pairing of A-002's decider produces cycles long enough to need this, and the
/// measured cost of it across the seven reference fields is **six vertices and
/// twelve triangles, all on one field at one resolution.**
pub const CENTROID_BASE: u8 = EDGE_COUNT as u8;

/// Is this triangle corner code a cycle centroid rather than a cube edge?
#[inline]
#[must_use]
pub const fn is_centroid(code: u8) -> bool {
    code >= CENTROID_BASE && code != NO_EDGE
}

/// The triangulation for one corner-sign configuration.
#[derive(Clone, Copy, Debug)]
pub struct McCase {
    /// Number of triangles.
    pub count: u8,
    /// How many cycle centroids the triangles reference.
    pub centroids: u8,
    /// Each triangle as three corner codes — a cube edge below
    /// [`CENTROID_BASE`], a cycle centroid at or above it. Only the first
    /// `count` are valid.
    pub triangles: [[u8; 3]; MAX_TRIANGLES],
}

/// A cube has six faces, so a per-face decision is six bits.
pub const FACE_COUNT: usize = 6;

/// The bit a face occupies in a resolution mask: `axis * 2 + side`.
///
/// Stated once here because [`segment_links`], [`AMBIGUOUS_FACES`],
/// [`super::ambiguity`]'s mask builder and the validator all have to agree about
/// it, and a transposition between any two of them would be invisible in the
/// output.
#[inline]
#[must_use]
pub const fn face_bit(axis: usize, side: u8) -> u8 {
    1 << (axis * 2 + side as usize)
}

/// All 256 configurations, built during compilation.
///
/// The all-separate resolution — mask zero — which is Marching Cubes proper.
pub static CASES: [McCase; 256] = build_cases();

const fn build_cases() -> [McCase; 256] {
    let mut out = [McCase {
        count: 0,
        centroids: 0,
        triangles: [[0u8; 3]; MAX_TRIANGLES],
    }; 256];
    let mut case = 0usize;
    while case < 256 {
        out[case] = triangulate(segment_links(case as u8, 0));
        case += 1;
    }
    out
}

/// Which of a case's faces are ambiguous, as a mask over [`face_bit`].
///
/// A face is ambiguous when its four corner signs alternate around the ring, so
/// all four of its edges are cut and the two inside corners can be joined or
/// separated. Those are the only faces on which [`segment_links`]'s mask does
/// anything.
pub static AMBIGUOUS_FACES: [u8; 256] = build_ambiguous_faces();

const fn build_ambiguous_faces() -> [u8; 256] {
    let mut out = [0u8; 256];
    let mut case = 0usize;
    while case < 256 {
        let mut mask = 0u8;
        let mut axis = 0usize;
        while axis < 3 {
            let mut side = 0u8;
            while side < 2 {
                let c = face_corners(axis, side);
                // Opposite corners agree, adjacent ones differ.
                let a = corner_inside(case as u8, c[0]);
                let b = corner_inside(case as u8, c[1]);
                if a == corner_inside(case as u8, c[2])
                    && b == corner_inside(case as u8, c[3])
                    && a != b
                {
                    mask |= face_bit(axis, side);
                }
                side += 1;
            }
            axis += 1;
        }
        out[case] = mask;
        case += 1;
    }
    out
}

/// Segment links for one configuration: `next[e]` is the cut edge the segment
/// leaving cut edge `e` arrives at, or [`NO_EDGE`].
///
/// `joined` selects, per face, which of the two pairings an ambiguous face uses:
/// clear pairs each entry with the exit that closes its own run of inside
/// corners, **separating** them; set crosses the pairing, **joining** them. Bits
/// on faces that are not ambiguous are ignored, since such a face has at most one
/// entry and therefore no pairing to choose. See [`AMBIGUOUS_FACES`] and
/// [`face_bit`].
///
/// Exposed so the validator can rebuild them independently of the triangles.
pub const fn segment_links(case: u8, joined: u8) -> [u8; EDGE_COUNT] {
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
                // Starting outside makes the transitions alternate entry, exit,
                // entry, exit, so entry `n` is always seen before exit `n`. Two
                // of each is the maximum: a face has 0, 2 or 4 cut edges.
                let mut entries = [NO_EDGE; 2];
                let mut exits = [NO_EDGE; 2];
                let mut pairs = 0usize;
                let mut j = 0usize;
                while j < 4 {
                    let p = c[(start + j) % 4];
                    let q = c[(start + j + 1) % 4];
                    let p_in = corner_inside(case, p);
                    let q_in = corner_inside(case, q);
                    if !p_in && q_in {
                        entries[pairs] = edge_index(p, q);
                    } else if p_in && !q_in {
                        exits[pairs] = edge_index(p, q);
                        pairs += 1;
                    }
                    j += 1;
                }

                if pairs == 2 && joined & face_bit(axis, side) != 0 {
                    next[entries[0] as usize] = exits[1];
                    next[entries[1] as usize] = exits[0];
                } else {
                    let mut n = 0usize;
                    while n < pairs {
                        next[entries[n] as usize] = exits[n];
                        n += 1;
                    }
                }
            }

            side += 1;
        }
        axis += 1;
    }

    next
}

/// The first vertex of this cycle that can serve as a fan apex without creating
/// a chord two cells could both name, or `len` if there is none.
///
/// `cycle[..len]` is one directed cycle of cut edges, as [`segment_links`] links
/// them.
///
/// A fan from apex `a` creates the chords `(c[a], c[a+i])` for `2 <= i <= len-2`.
/// Such a chord is safe when its two cut edges share no cube face, because then
/// no *other* cell contains both and the mesh edge can only come from this cell's
/// own fan — which emits it exactly twice. A three-cycle has no chords at all, so
/// its apex is trivially safe.
///
/// **Measured over all 256 cases and every canonical resolution mask:** every
/// cycle of length 3 to 7 has a safe apex, as do 48 of the 60 length-8 cycles;
/// the length-9 and length-12 cycles have none. Plain Marching Cubes tops out at
/// length 7, so it never needs a centroid — the cost of A-015 falls entirely on
/// the long cycles the asymptotic decider's joined pairing can produce.
pub const fn safe_apex(cycle: &[u8; EDGE_COUNT], len: usize) -> usize {
    let mut a = 0usize;
    while a < len {
        let mut safe = true;
        let mut i = 2usize;
        while i + 2 <= len {
            if edges_share_a_face(cycle[a], cycle[(a + i) % len]) {
                safe = false;
                break;
            }
            i += 1;
        }
        if safe {
            return a;
        }
        a += 1;
    }
    len
}

/// Triangulate the directed cycles the segments form.
///
/// Split out of the table construction so the runtime decider path and the
/// compile-time table share one implementation and cannot drift.
///
/// A three-cycle becomes its one triangle. **Anything longer fans from a
/// centroid** rather than from one of its own vertices — see [`CENTROID_BASE`]
/// for why (A-015): a fan from a vertex leaves `k − 3` interior chords, and two
/// cells sharing a face can choose the same chord and put four triangles on one
/// mesh edge. A centroid is cell-local, so its spokes cannot be named by any
/// other cell, and every remaining mesh edge is a cycle edge lying on a face,
/// which exactly two cells share.
///
/// This costs a vertex and two triangles per long cycle. It does **not** change
/// the surface, its Euler characteristic or its genus: both triangulations of a
/// `k`-gon are discs, `χ = 1`.
#[must_use]
pub const fn triangulate(next: [u8; EDGE_COUNT]) -> McCase {
    let mut triangles = [[0u8; 3]; MAX_TRIANGLES];
    let mut count = 0usize;
    let mut centroids = 0usize;
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

            // Every cycle has at least three edges: two distinct cube edges share
            // at most one face, so no two-cycle can form. That holds under either
            // pairing — crossing an ambiguous face's pairing permutes which exit
            // each entry reaches, and both exits are still on that one shared
            // face.
            // Fan from the first chord-safe apex; a cycle with none gets a
            // centroid. See `safe_apex` and [`CENTROID_BASE`].
            let apex = safe_apex(&cycle, len);
            if apex < len {
                let mut i = 2usize;
                while i < len {
                    triangles[count] = [
                        cycle[apex],
                        cycle[(apex + i - 1) % len],
                        cycle[(apex + i) % len],
                    ];
                    count += 1;
                    i += 1;
                }
            } else {
                let centre = CENTROID_BASE + centroids as u8;
                centroids += 1;
                let mut i = 0usize;
                while i < len {
                    // Winding follows the cycle, so the centroid fan keeps the
                    // orientation a vertex fan would have had.
                    triangles[count] = [centre, cycle[i], cycle[(i + 1) % len]];
                    count += 1;
                    i += 1;
                }
            }
        }
        e += 1;
    }

    McCase {
        count: count as u8,
        centroids: centroids as u8,
        triangles,
    }
}
