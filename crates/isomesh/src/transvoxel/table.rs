//! The transition cell, and its contour derived rather than transcribed.
//!
//! # What a transition cell is
//!
//! Two adjacent blocks of voxel data, one sampled at half the other's
//! resolution. Lengyel 2010 (`transvoxel_dissertation_lengyel2010`, §4.3, read
//! this session):
//!
//! > We call a cell inside the half-resolution block that lies along the border
//! > with the full-resolution block a **transition cell**, and we see that a
//! > triangulation compatible with the Marching Cubes algorithm for any such cell
//! > must account for a total of 13 sample values lying on its boundary… nine of
//! > these samples come from the full-resolution data, and the remaining four
//! > come from the half-resolution data.
//!
//! Taken directly that is 2¹³ cases, and 2¹⁷ and 2²⁰ for cells bordered on two
//! and three faces — 1,187,840 configurations, which the dissertation calls *"a
//! monumentally tedious job"* costing *"dozens of megabytes"*. Its fix is to
//! split the cell in two:
//!
//! > The left cell is triangulated using sample values only from the face
//! > bordering the full-resolution block, and the four corner values labeled A,
//! > B, C, and D in the figure are duplicated on the opposite face of the cell.
//! > The transition from full resolution to half resolution takes place entirely
//! > inside this cell, for which we now have a very manageable nine distinct
//! > sample values to consider. The right cell is triangulated in the
//! > conventional manner using only the eight sample values from the
//! > half-resolution data.
//!
//! So **this module is only the left cell**: nine samples, 2⁹ = 512 cases. The
//! right cell is [`MarchingCubes`](crate::marching_cubes::MarchingCubes),
//! unchanged, and that is the whole reason the split is worth making.
//!
//! # Why it is derived
//!
//! Lengyel publishes `transitionCellClass` and `transitionCellData`. Transcribing
//! them is exactly what hard rule 5 exists to prevent — a mistyped case table
//! produces a mesh that looks fine and is subtly non-manifold — and A-001 already
//! established the alternative here: **derive the contour from the cell's own
//! boundary and check the derivation against the paper's published *count*.**
//!
//! The dissertation supplies that count, and it is a demanding one:
//!
//! > Observing the orbits of each of the 512 cases under the action of the group
//! > `D₈` or `D₈ × ℤ₂`, as appropriate, yields exactly **73 distinct equivalence
//! > classes**. One equivalence class is the trivial class containing the two
//! > cases for which the inside/outside state of all nine sample values is the
//! > same.
//!
//! `the_orbit_count_is_lengyels_seventy_three` reproduces it from the derived
//! table. A transposition anywhere in the geometry below moves that number.
//!
//! # The cell as nine faces
//!
//! The insight that makes this the same problem A-001 already solved: the left
//! cell's boundary is a closed oriented surface made of **nine flat faces**, and
//! a contour on it is what you get by pairing crossings within each face and
//! following the pairs from face to face — which is exactly
//! [`segment_links`](crate::marching_cubes::table::segment_links) generalised off
//! the cube.
//!
//! ```text
//!   full-resolution face (w = 0)        half-resolution face (w = 1)
//!        6 --- 7 --- 8                        6 --------- 8
//!        |  q2 |  q3 |                        |           |
//!        3 --- 4 --- 5                        |   back    |
//!        |  q0 |  q1 |                        |           |
//!        0 --- 1 --- 2                        0 --------- 2
//! ```
//!
//! - **four quadrants** on the full-resolution face, each an ordinary four-corner
//!   face and each independently ambiguous,
//! - **one back face**, also four-cornered and also independently ambiguous,
//!   whose corner *values* are samples 0, 2, 8, 6 duplicated,
//! - **four lateral faces**, each of which reduces to a **triangle**. Its two
//!   depth edges join a sample to a duplicate of itself, so neither can ever be
//!   crossed, and what is left is two full-resolution sub-edges and one
//!   half-resolution edge. A three-cornered face carries 0 or 2 crossings and has
//!   no ambiguity to resolve — which is the whole transition, and why it needs no
//!   special case.
//!
//! Sixteen edges can be cut: twelve on the full-resolution face and four on the
//! half-resolution face. The numbering below is this crate's own; no paper's is
//! being matched, which is what makes the table derivable.

use crate::marching_cubes::table::NO_EDGE;

/// Samples on the full-resolution face. The half-resolution face duplicates four
/// of them, so nine values determine the cell.
pub const SAMPLE_COUNT: usize = 9;

/// Cuttable edges: twelve on the full-resolution face, four on the
/// half-resolution face.
///
/// The four depth edges are omitted deliberately rather than counted and always
/// skipped — each joins a sample to its own duplicate, so its endpoints are equal
/// by construction and it cannot carry a sign change.
pub const EDGE_COUNT: usize = 16;

/// The two samples each edge joins.
///
/// `0..6` run along `u` on the full-resolution face, `6..12` along `v`, and
/// `12..16` are the half-resolution face's four edges. A half-resolution edge
/// joins the *same* two samples as the pair of full-resolution sub-edges beside
/// it — that duplication is the transition, and it is why edges 12–15 are
/// separate entries rather than aliases.
pub const EDGE_SAMPLES: [[u8; 2]; EDGE_COUNT] = [
    // full-resolution, along u
    [0, 1],
    [1, 2],
    [3, 4],
    [4, 5],
    [6, 7],
    [7, 8],
    // full-resolution, along v
    [0, 3],
    [3, 6],
    [1, 4],
    [4, 7],
    [2, 5],
    [5, 8],
    // half-resolution
    [0, 2],
    [2, 8],
    [6, 8],
    [0, 6],
];

/// `true` for the four half-resolution edges.
///
/// The only place the two resolutions are told apart, and A-011b needs it: a
/// crossing on one of these sits on the coarse neighbour's grid edge and must
/// land at the position that neighbour's own Marching Cubes pass would give it,
/// or the seam does not close.
#[inline]
#[must_use]
pub const fn is_half_resolution(edge: u8) -> bool {
    edge >= 12
}

/// How many faces bound the cell: four quadrants, one back face, four laterals.
pub const FACE_COUNT: usize = 9;

/// The longest face boundary. Quadrants and the back face are four-cornered;
/// laterals reduce to three.
const MAX_FACE_CORNERS: usize = 4;

/// One face of the cell: its corners in order, and the edge leaving each.
///
/// Corners run **counter-clockwise seen from outside the cell**, which is what
/// makes a crossing an *entry* on one of its two faces and an *exit* on the
/// other, and therefore what makes the links form consistently oriented cycles.
/// `edge[k]` is the edge from `sample[k]` to `sample[k + 1]`.
struct Face {
    sample: [u8; MAX_FACE_CORNERS],
    edge: [u8; MAX_FACE_CORNERS],
    corners: usize,
}

/// The nine faces, wound counter-clockwise seen from outside.
///
/// Each edge appears in exactly two faces and in **opposite directions**, which
/// `every_edge_is_shared_by_two_faces_in_opposite_directions` asserts rather than
/// assumes — an orientation slip here would produce cycles that close and are
/// wound inconsistently, and no Euler or manifold check can see that.
const FACES: [Face; FACE_COUNT] = [
    // The four full-resolution quadrants. Seen from outside the cell the
    // full-resolution face is the *near* face along w, so counter-clockwise runs
    // (i,j) -> (i,j+1) -> (i+1,j+1) -> (i+1,j) — the reverse of the far face's
    // order, exactly as `cube::face_corners` documents for a near face.
    Face {
        sample: [0, 3, 4, 1],
        edge: [6, 2, 8, 0],
        corners: 4,
    },
    Face {
        sample: [1, 4, 5, 2],
        edge: [8, 3, 10, 1],
        corners: 4,
    },
    Face {
        sample: [3, 6, 7, 4],
        edge: [7, 4, 9, 2],
        corners: 4,
    },
    Face {
        sample: [4, 7, 8, 5],
        edge: [9, 5, 11, 3],
        corners: 4,
    },
    // The half-resolution face is the far face along w, so counter-clockwise runs
    // the other way round.
    Face {
        sample: [0, 2, 8, 6],
        edge: [12, 13, 14, 15],
        corners: 4,
    },
    // The four laterals, each a triangle: two full-resolution sub-edges and the
    // half-resolution edge that closes it.
    Face {
        sample: [0, 1, 2, 0],
        edge: [0, 1, 12, NO_EDGE],
        corners: 3,
    },
    Face {
        sample: [2, 5, 8, 0],
        edge: [10, 11, 13, NO_EDGE],
        corners: 3,
    },
    Face {
        sample: [6, 8, 7, 0],
        edge: [14, 5, 4, NO_EDGE],
        corners: 3,
    },
    Face {
        sample: [0, 6, 3, 0],
        edge: [15, 7, 6, NO_EDGE],
        corners: 3,
    },
];

/// Is sample `s` inside the solid, for case index `case`?
#[inline]
#[must_use]
pub const fn sample_inside(case: u16, sample: u8) -> bool {
    case & (1 << sample) != 0
}

/// Which faces of this case are ambiguous, as a mask over face index.
///
/// A face is ambiguous when its four corner signs alternate around the ring, so
/// all four of its edges are cut and the two inside corners can be joined or
/// separated. The laterals are three-cornered and can never be ambiguous, which
/// is why the transition itself needs no decision.
///
/// The dissertation names exactly these: *"does not have an ambiguous
/// half-resolution face and does not have any ambiguous quadrants on its
/// full-resolution face."*
pub static AMBIGUOUS_FACES: [u16; 512] = build_ambiguous_faces();

const fn build_ambiguous_faces() -> [u16; 512] {
    let mut out = [0u16; 512];
    let mut case = 0usize;
    while case < 512 {
        let mut mask = 0u16;
        let mut f = 0usize;
        while f < FACE_COUNT {
            if FACES[f].corners == 4 {
                let c = &FACES[f].sample;
                let a = sample_inside(case as u16, c[0]);
                let b = sample_inside(case as u16, c[1]);
                if a == sample_inside(case as u16, c[2])
                    && b == sample_inside(case as u16, c[3])
                    && a != b
                {
                    mask |= 1 << f;
                }
            }
            f += 1;
        }
        out[case] = mask;
        case += 1;
    }
    out
}

/// Segment links for one configuration: `next[e]` is the cut edge the segment
/// leaving cut edge `e` arrives at, or [`NO_EDGE`].
///
/// The direct analogue of
/// [`segment_links`](crate::marching_cubes::table::segment_links), and it works
/// the same way for the same reason: walking a face from an **outside** corner
/// makes its transitions alternate entry, exit, entry, exit, so entry `n` is
/// always seen before exit `n` and the pairing needs no sorting.
///
/// `joined` selects, per ambiguous face, which of the two pairings it uses —
/// clear separates the two inside corners, set joins them. Bits on faces that are
/// not ambiguous are ignored. See [`AMBIGUOUS_FACES`].
#[must_use]
pub const fn transition_links(case: u16, joined: u16) -> [u8; EDGE_COUNT] {
    let mut next = [NO_EDGE; EDGE_COUNT];

    let mut f = 0usize;
    while f < FACE_COUNT {
        let face = &FACES[f];
        let n = face.corners;

        // Start on an outside corner. If every corner is inside, the face is not
        // cut and there is nothing to pair.
        let mut start = n;
        let mut k = 0usize;
        while k < n {
            if !sample_inside(case, face.sample[k]) {
                start = k;
                break;
            }
            k += 1;
        }

        if start < n {
            let mut entries = [NO_EDGE; 2];
            let mut exits = [NO_EDGE; 2];
            let mut pairs = 0usize;
            let mut j = 0usize;
            while j < n {
                let from = (start + j) % n;
                let to = (start + j + 1) % n;
                let p_in = sample_inside(case, face.sample[from]);
                let q_in = sample_inside(case, face.sample[to]);
                if !p_in && q_in {
                    entries[pairs] = face.edge[from];
                } else if p_in && !q_in {
                    exits[pairs] = face.edge[from];
                    pairs += 1;
                }
                j += 1;
            }

            if pairs == 2 && joined & (1 << f) != 0 {
                next[entries[0] as usize] = exits[1];
                next[entries[1] as usize] = exits[0];
            } else {
                let mut m = 0usize;
                while m < pairs {
                    next[entries[m] as usize] = exits[m];
                    m += 1;
                }
            }
        }

        f += 1;
    }

    next
}

/// Which edges of this case are cut, as a bitmask.
///
/// Derived from the sample signs alone, so it is the independent check that
/// [`transition_links`] covers every cut edge and no others.
#[must_use]
pub const fn cut_edges(case: u16) -> u16 {
    let mut mask = 0u16;
    let mut e = 0usize;
    while e < EDGE_COUNT {
        let [a, b] = EDGE_SAMPLES[e];
        if sample_inside(case, a) != sample_inside(case, b) {
            mask |= 1 << e;
        }
        e += 1;
    }
    mask
}

// ─── the group action, for the equivalence-class count ──────────────────────

/// A sample permutation: `image[s]` is where sample `s` goes.
type Perm = [u8; SAMPLE_COUNT];

/// 90-degree rotation of the 3x3 arrangement, counter-clockwise.
///
/// Sample `s` sits at `(i, j) = (s % 3, s / 3)`; the rotation sends `(i, j)` to
/// `(j, 2 - i)`.
const fn rotate() -> Perm {
    let mut out = [0u8; SAMPLE_COUNT];
    let mut s = 0usize;
    while s < SAMPLE_COUNT {
        let i = s % 3;
        let j = s / 3;
        out[s] = (j + 3 * (2 - i)) as u8;
        s += 1;
    }
    out
}

/// Flip about the horizontal axis: `(i, j)` goes to `(i, 2 - j)`.
const fn flip() -> Perm {
    let mut out = [0u8; SAMPLE_COUNT];
    let mut s = 0usize;
    while s < SAMPLE_COUNT {
        let i = s % 3;
        let j = s / 3;
        out[s] = (i + 3 * (2 - j)) as u8;
        s += 1;
    }
    out
}

/// The eight elements of `D₈`, as sample permutations.
///
/// Lengyel's own generators: *"we have chosen r to be the 90-degree
/// counterclockwise rotation and f to be the flip about the horizontal axis."*
pub static DIHEDRAL: [Perm; 8] = build_dihedral();

const fn build_dihedral() -> [Perm; 8] {
    let r = rotate();
    let f = flip();
    let mut out = [[0u8; SAMPLE_COUNT]; 8];

    // Identity, then successive rotations.
    let mut s = 0usize;
    while s < SAMPLE_COUNT {
        out[0][s] = s as u8;
        s += 1;
    }
    let mut k = 1usize;
    while k < 4 {
        let mut s = 0usize;
        while s < SAMPLE_COUNT {
            out[k][s] = r[out[k - 1][s] as usize];
            s += 1;
        }
        k += 1;
    }
    // Each rotation composed with the flip.
    let mut k = 0usize;
    while k < 4 {
        let mut s = 0usize;
        while s < SAMPLE_COUNT {
            out[4 + k][s] = f[out[k][s] as usize];
            s += 1;
        }
        k += 1;
    }
    out
}

/// Apply a sample permutation to a case index.
#[must_use]
pub const fn permute(case: u16, perm: &Perm) -> u16 {
    let mut out = 0u16;
    let mut s = 0usize;
    while s < SAMPLE_COUNT {
        if case & (1 << s) != 0 {
            out |= 1 << perm[s];
        }
        s += 1;
    }
    out
}

/// Invert every sample's inside/outside state.
#[inline]
#[must_use]
pub const fn invert(case: u16) -> u16 {
    !case & 0x1ff
}

#[cfg(test)]
mod tests;
