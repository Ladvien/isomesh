//! §3.1 — the curves a tet's boundary carries, from arbitrary edge coordinates.
//!
//! [`coordinates`](super::coordinates) says *how many* times the surface crosses
//! each edge. This says which crossing joins which, on every triangular face, for
//! **any** non-negative edge coordinates — including the ones that are not normal
//! curve systems at all and therefore have no lookup-table answer.
//!
//! # The property the whole method rests on
//!
//! Baktash, Gillespie & Crane, `10.48550/arXiv.2606.00454` §3.1:
//!
//! > We carefully define this procedure to produce **identical curves on
//! > triangles shared by neighboring tetrahedra**.
//!
//! That is what makes a per-tet algorithm produce a conforming mesh — no
//! communication between tets, no second pass, no crack.
//!
//! **The mechanism is locality, not symmetry**, and the distinction cost a test
//! to learn. Every decision below reads the face's own three edge coordinates and
//! the grid's **canonical** edge orientation (lower vertex index first) — so two
//! tets, which differ only in their fourth vertex and the edges to it, cannot
//! disagree. `a_shared_face_gives_the_same_segments_from_either_tet` varies those
//! other three edges every way and requires the face not to move.
//!
//! What the procedure is emphatically *not* is invariant under relabelling the
//! face's corners: step 3(b) skips *"the first residual point … assuming a
//! canonical orientation `i < j`"*, so swapping `i` and `j` skips the other end.
//! That is fine, and it is fine for a specific reason — neighbouring tets share
//! the face's **global vertex indices**, so they agree on `i < j` without ever
//! comparing notes. See M-79.
//!
//! # The procedure, verbatim
//!
//! For each triangle `ijk`:
//!
//! > (1) If the edge coordinates satisfy the even sum and triangle inequality
//! > conditions, we compute corner coordinates via Equation 2. For each corner
//! > `i` we then connect up the first `cᵢ` pairs of intersection points along
//! > oriented edges `ij` and `ik` into segments.
//! >
//! > (2) If the triangle inequality is satisfied, but the even sum condition is
//! > violated, then we subtract 1 from each edge coordinate and return to Step
//! > (1). Doing so ensures that the sum is now even, effectively creating three
//! > open endpoints.
//! >
//! > (3) Finally, if the triangle inequality is violated, there will be exactly
//! > one edge `ij` with `r := eᵢⱼ − eⱼₖ − eₖᵢ` residual points… we first construct
//! > as many corner cuts as we can, by applying Step (1) to adjusted edge
//! > coordinates `e′ᵢⱼ := min(eᵢⱼ, eⱼₖ + eₖᵢ)`. Then, to handle residual points:
//! > (a) If `r` is even, we connect consecutive pairs of points along residual
//! > edge `ij` into segments… (b) If `r` is odd, we skip the first residual point
//! > along oriented edge `ij` (assuming a canonical orientation `i < j`), and
//! > connect the remaining consecutive pairs as in the even case.
//!
//! **Nothing here is left to a choice.** Steps 2 and 3 both leave points
//! unclaimed, and *which* points is forced rather than picked: after step 2 the
//! corner cuts claim `cᵢ` from one end and `cⱼ` from the other, and
//! `cᵢ + cⱼ = eᵢⱼ − 1`, so exactly one point — the one at index `cᵢ` — is left
//! over on each edge. `the_reduced_pass_leaves_exactly_one_point_per_edge` checks
//! that rather than trusting the arithmetic.

use alloc::vec::Vec;

use crate::marching_tetrahedra::table::TET_EDGES;

use super::coordinates::{EdgeCoordinates, TET_FACES};

/// One crossing on one edge.
///
/// `index` counts from the edge's **lower-numbered corner**, which is the
/// canonical orientation `TET_EDGES` already stores. Indexing from a *face's*
/// corner would make the same point have two names depending on which tet asked,
/// and the conformity property would be gone.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub struct FacePoint {
    /// Which tet edge the crossing lies on.
    pub edge: u8,
    /// Which crossing along it, from the edge's lower corner.
    pub index: u32,
}

/// A piece of curve joining two crossings on one face.
///
/// Stored with its endpoints in sorted order so two segments are equal exactly
/// when they join the same two points — which is what the conformity test
/// compares.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub struct Segment {
    /// The lower endpoint.
    pub a: FacePoint,
    /// The upper endpoint.
    pub b: FacePoint,
}

impl Segment {
    fn new(a: FacePoint, b: FacePoint) -> Self {
        if a <= b {
            Self { a, b }
        } else {
            Self { a: b, b: a }
        }
    }

    /// Whether both endpoints lie on the same edge.
    ///
    /// The paper's *"scoop"* — a segment that runs along a tet edge rather than
    /// between two of them, and the thing that makes a curve **non-normal**.
    #[must_use]
    pub fn is_scoop(&self) -> bool {
        self.a.edge == self.b.edge
    }
}

/// The crossing `t` places along edge `edge` **from corner `from`**, expressed in
/// the edge's canonical indexing.
///
/// The one place the two orientations are reconciled, so a sign error here shows
/// up as a conformity failure rather than as a subtly wrong curve.
fn from_corner(edge: u8, from: u8, t: u32, count: u32) -> FacePoint {
    let [lo, _hi] = TET_EDGES[edge as usize];
    let index = if from == lo { t } else { count - 1 - t };
    FacePoint { edge, index }
}

/// Every segment one face of the tet contributes.
///
/// `face` indexes [`TET_FACES`]. The result is sorted, so it can be compared
/// directly against the same face seen from a neighbouring tet.
#[must_use]
pub fn face_segments(face: u8, coords: &EdgeCoordinates) -> Vec<Segment> {
    let f = TET_FACES[face as usize];
    // Edge `k` joins corner `k` to corner `k + 1`, so the edge *opposite* corner
    // `k` is edge `k + 1`.
    let count: [u32; 3] = [
        coords.edge(f.edge[0]),
        coords.edge(f.edge[1]),
        coords.edge(f.edge[2]),
    ];

    let mut out = Vec::new();
    let total = count[0] + count[1] + count[2];

    // Which corner, if any, sits opposite an edge too long for the other two.
    // At most one can: two such edges would each exceed the sum of the others.
    let violated = (0..3).find(|k| count[*k] > count[(*k + 1) % 3] + count[(*k + 2) % 3]);

    let (reduced, residual_edge) = match violated {
        // Step 3: clamp the long edge, and remember it for the residual pass.
        Some(k) => {
            let mut r = count;
            r[k] = count[(k + 1) % 3] + count[(k + 2) % 3];
            (r, Some(k))
        }
        // Step 2: an odd sum with the inequality intact -- drop one from each.
        None if !total.is_multiple_of(2) => ([count[0] - 1, count[1] - 1, count[2] - 1], None),
        // Step 1 applies as written.
        None => (count, None),
    };

    // Step 1, on whichever coordinates the branch above produced. Corner `k`'s
    // coordinate counts the arcs cutting it off; each joins edge `k` to edge
    // `k - 1`, both indexed from corner `k`.
    let mut corner = [0u32; 3];
    for (k, slot) in corner.iter_mut().enumerate() {
        let adjacent = reduced[k] + reduced[(k + 2) % 3];
        let opposite = reduced[(k + 1) % 3];
        // Every branch above restores the two conditions, so this is a
        // non-negative even number by construction.
        debug_assert!(adjacent >= opposite && (adjacent - opposite) % 2 == 0);
        *slot = (adjacent - opposite) / 2;
    }

    for k in 0..3usize {
        let at = f.corner[k];
        let (leaving, arriving) = (f.edge[k], f.edge[(k + 2) % 3]);
        for t in 0..corner[k] {
            out.push(Segment::new(
                from_corner(leaving, at, t, count[k]),
                from_corner(arriving, at, t, count[(k + 2) % 3]),
            ));
        }
    }

    // Step 3's residual pass: whatever the corner cuts could not absorb, paired
    // along the long edge itself. These are the scoops.
    if let Some(k) = residual_edge {
        let edge = f.edge[k];
        let n = count[k];
        // Claimed from corner `k`: the first `corner[k]`. Claimed from corner
        // `k + 1`: the last `corner[(k + 1) % 3]`. What is left is the middle
        // run, expressed in canonical indexing.
        let low = from_corner(edge, f.corner[k], corner[k], n).index;
        let high = from_corner(edge, f.corner[(k + 1) % 3], corner[(k + 1) % 3], n).index;
        let (first, last) = if low <= high {
            (low, high)
        } else {
            (high, low)
        };

        let mut point = first;
        let r = last + 1 - first;
        // An odd residual cannot pair up; the paper skips the *first* in
        // canonical order, which is what makes both tets skip the same one.
        if r % 2 != 0 {
            point += 1;
        }
        while point < last {
            out.push(Segment::new(
                FacePoint { edge, index: point },
                FacePoint {
                    edge,
                    index: point + 1,
                },
            ));
            point += 2;
        }
    }

    out.sort_unstable();
    out
}

/// Every segment the whole tet boundary carries.
#[must_use]
pub fn boundary_segments(coords: &EdgeCoordinates) -> Vec<Segment> {
    let mut out = Vec::new();
    for face in 0..TET_FACES.len() as u8 {
        out.extend(face_segments(face, coords));
    }
    out.sort_unstable();
    out
}

/// What kind of curve a connected component of the segment graph is.
///
/// §3.1's own three-way split, and each arm has a different fate in §3.2.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum CurveKind {
    /// *"If `γ` has a degree-1 vertex, it forms an open curve. Such curves are
    /// discarded, but their segments may still appear in neighboring tets as part
    /// of the mesh boundary."*
    Open,
    /// *"If all segments of `γ` connect two distinct edges, then `γ` is a normal
    /// curve, and will ultimately bound a disk on the tet interior."*
    Normal,
    /// *"Otherwise, if any segment of `γ` runs along a tet edge, it is a
    /// non-normal curve, and its spanning disk can include pieces of the tet
    /// boundary."*
    NonNormal,
}

/// One connected curve on the tet boundary.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Curve {
    /// Which of §3.1's three kinds this is.
    pub kind: CurveKind,
    /// Its segments, sorted.
    pub segments: Vec<Segment>,
}

/// Group a tet's boundary segments into curves and classify each.
#[must_use]
pub fn curves(coords: &EdgeCoordinates) -> Vec<Curve> {
    let segments = boundary_segments(coords);

    // Degree of every endpoint, for the open-curve test.
    let mut points: Vec<FacePoint> = segments.iter().flat_map(|s| [s.a, s.b]).collect();
    points.sort_unstable();
    points.dedup();
    let degree = |p: FacePoint| -> usize {
        segments
            .iter()
            .filter(|s| s.a == p || s.b == p)
            .map(|s| usize::from(s.a == p) + usize::from(s.b == p))
            .sum()
    };

    // Connected components over the segment graph.
    let mut unassigned: Vec<usize> = (0..segments.len()).collect();
    let mut out = Vec::new();
    while let Some(seed) = unassigned.pop() {
        let mut group = alloc::vec![seed];
        let mut frontier = alloc::vec![seed];
        while let Some(current) = frontier.pop() {
            let s = segments[current];
            let mut still = Vec::new();
            for other in unassigned.drain(..) {
                let t = segments[other];
                if t.a == s.a || t.a == s.b || t.b == s.a || t.b == s.b {
                    group.push(other);
                    frontier.push(other);
                } else {
                    still.push(other);
                }
            }
            unassigned = still;
        }

        let mut mine: Vec<Segment> = group.into_iter().map(|i| segments[i]).collect();
        mine.sort_unstable();
        let ends: Vec<FacePoint> = {
            let mut p: Vec<FacePoint> = mine.iter().flat_map(|s| [s.a, s.b]).collect();
            p.sort_unstable();
            p.dedup();
            p
        };
        let kind = if ends.iter().any(|p| degree(*p) < 2) {
            CurveKind::Open
        } else if mine.iter().any(Segment::is_scoop) {
            CurveKind::NonNormal
        } else {
            CurveKind::Normal
        };
        out.push(Curve {
            kind,
            segments: mine,
        });
    }
    out.sort_by(|a, b| a.segments.cmp(&b.segments));
    out
}

#[cfg(test)]
mod tests;
