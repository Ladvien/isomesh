//! §3.2 — the surface that fills the curves, from arbitrary edge coordinates.
//!
//! [`curves`](super::curves) says which crossings join which on the tet's
//! boundary. This says what spans them. Baktash, Gillespie & Crane,
//! `10.48550/arXiv.2606.00454` §3.2.1, reduce the problem to four cases:
//!
//! > We first consider the subset of normal curves `Γ_normal ⊂ Γ`, reducing the
//! > general case to four cases we can triangulate directly: **multiple
//! > triangles, parallel quads or octagons, or a single closed loop.**
//!
//! The reduction is what makes the case list finite, and it rests on two
//! properties that hold *only after the corner cuts are removed*:
//!
//! > We first handle corner cuts, i.e., we emit a triangle for each loop
//! > `γ ∈ Γ_normal` of length `ℓ = 3`. We can then establish two properties:
//! >
//! > I. All remaining loops have the same length `ℓ > 3` (Theorem B.6).
//! > II. These loops have edge coordinates `d₁`, `d₂`, and `d₁ + d₂` on opposite
//! > pairs of edges, for some integers `0 ≤ d₂ ≤ d₁` (Theorem B.3).
//!
//! Property II is why [`Pattern`] exists and why it is checked rather than
//! assumed: it is a strong claim about the residual coordinates, and if it ever
//! fails, every case below it is being applied to something it was not derived
//! for.
//!
//! # Why `ℓ > 4` means `ℓ ≥ 8`
//!
//! > the length of any non-triangular normal loop on a tet is a multiple of 4
//! > (Property I, Theorem B.6).
//!
//! The appendix gives the exact form, and it is worth having as an equality
//! rather than a bound — Corollary B.6: the length of every component in a
//! `(d₁, d₂)` pattern is exactly `4(d₁ + d₂) / gcd(d₁, d₂)`, and Theorem B.4:
//! the number of components is exactly `gcd(d₁, d₂)`. Both are computed in
//! [`Pattern::of`] and both are checked against the curves actually found, which
//! is a far stronger test than "the case list did not panic".
//!
//! # What this module does not do
//!
//! It produces **positions and indices, not normals.** The crate's convention is
//! that normals are the field's own gradient (see
//! [`marching_tetrahedra`](crate::marching_tetrahedra)), and a Steiner point
//! placed in the tet interior has no gradient until something samples one there.
//! Rather than invent an averaged normal for it, the fill emits geometry into a
//! [`TetPatch`] and A-014c's extractor — which has the field in scope — attaches
//! gradients as it copies into the sink. One path, and no vertex carries a
//! normal that came from anywhere but the field.

use alloc::vec::Vec;

use crate::marching_tetrahedra::table::{TET_EDGE_COUNT, TET_EDGES};
use crate::real::Real;

use super::coordinates::{EdgeCoordinates, TET_FACE_COUNT, complementary};
use super::curves::{Curve, CurveKind, FacePoint, Segment, curves};

/// A tetrahedron, and every crossing found along its edges.
///
/// The edge coordinates are **derived** from the slice lengths rather than
/// passed alongside them, so the two cannot disagree. That is the whole reason
/// this type exists instead of a pair of arguments.
#[derive(Clone, Copy, Debug)]
pub struct TetCrossings<'a, R: Real> {
    /// The four corner positions, in the corner order [`TET_EDGES`] indexes.
    pub corners: [[R; 3]; 4],
    /// Crossings along each edge, as parameters in `[0, 1]` measured from the
    /// edge's **lower-numbered** corner, in ascending order.
    ///
    /// Ascending order is the caller's contract and is checked by
    /// [`TetCrossings::check`] rather than assumed: the pairing rules in §3.1
    /// index crossings by position along the edge, so an unsorted list silently
    /// joins the wrong points.
    pub along: [&'a [R]; TET_EDGE_COUNT],
}

/// Why a [`TetCrossings`] cannot be filled.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum NotFillable {
    /// An edge's crossing parameters are not in ascending order.
    Unsorted {
        /// Which edge.
        edge: u8,
        /// The index whose parameter is not greater than its predecessor's.
        at: usize,
    },
    /// A crossing parameter is outside `[0, 1]` or is not finite.
    OffEdge {
        /// Which edge.
        edge: u8,
        /// Which crossing along it.
        at: usize,
    },
}

impl<R: Real> TetCrossings<'_, R> {
    /// The edge coordinates these crossings describe.
    #[must_use]
    pub fn coordinates(&self) -> EdgeCoordinates {
        let mut count = [0u32; TET_EDGE_COUNT];
        for (slot, along) in count.iter_mut().zip(self.along.iter()) {
            // A tet edge cannot carry more crossings than a u32 counts; the
            // caller would have exhausted memory first.
            *slot = u32::try_from(along.len()).unwrap_or(u32::MAX);
        }
        EdgeCoordinates { count }
    }

    /// Check the ascending-and-in-range contract.
    ///
    /// # Errors
    ///
    /// [`NotFillable`] naming the edge and the offending crossing.
    pub fn check(&self) -> Result<(), NotFillable> {
        for (e, along) in self.along.iter().enumerate() {
            let edge = e as u8;
            for (i, t) in along.iter().enumerate() {
                if !t.is_finite() || *t < R::ZERO || *t > R::ONE {
                    return Err(NotFillable::OffEdge { edge, at: i });
                }
                // `<=` rather than `!(>)`: the two differ only on NaN, which the
                // finiteness check above has already rejected, and the direct
                // form says "not strictly ascending" without the double
                // negative. Equal parameters are rejected too — two crossings at
                // the same point are one crossing, and the pairing rules index
                // by position.
                if i > 0 && *t <= along[i - 1] {
                    return Err(NotFillable::Unsorted { edge, at: i });
                }
            }
        }
        Ok(())
    }

    /// Where a crossing sits in space.
    ///
    /// [`FacePoint::index`] counts from the edge's lower-numbered corner, which
    /// is the same end `along` measures from — the two conventions are the same
    /// one, and this is the only place that depends on it.
    #[must_use]
    pub fn position(&self, point: FacePoint) -> Option<[R; 3]> {
        let [lo, hi] = TET_EDGES[point.edge as usize];
        let t = *self.along[point.edge as usize].get(point.index as usize)?;
        let (a, b) = (self.corners[lo as usize], self.corners[hi as usize]);
        Some([
            a[0] + (b[0] - a[0]) * t,
            a[1] + (b[1] - a[1]) * t,
            a[2] + (b[2] - a[2]) * t,
        ])
    }
}

/// The triangles one tetrahedron contributes, and the vertices they index.
///
/// Caller-provided and reusable per `CLAUDE.md` rule 6: [`reset`](Self::reset)
/// clears without releasing capacity, because a real workload fills one of these
/// per tet per chunk and re-meshes thousands of chunks per edit.
#[derive(Clone, Debug, Default)]
pub struct TetPatch<R: Real> {
    /// Vertex positions — edge crossings first, then any Steiner points.
    pub positions: Vec<[R; 3]>,
    /// Triangles, indexing [`positions`](Self::positions).
    pub triangles: Vec<[u32; 3]>,
}

impl<R: Real> TetPatch<R> {
    /// An empty patch.
    #[must_use]
    pub const fn new() -> Self {
        Self {
            positions: Vec::new(),
            triangles: Vec::new(),
        }
    }

    /// Clear without releasing capacity.
    pub fn reset(&mut self) {
        self.positions.clear();
        self.triangles.clear();
    }

    /// Whether anything was emitted.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.triangles.is_empty()
    }
}

/// A boundary curve as an ordered cycle, rather than a set of segments.
///
/// [`curves`] returns each component as an unordered [`Segment`] set, which is
/// what the conformity property is stated over. Triangulating one needs the
/// cyclic order, and needs it to be **deterministic** — two runs that walked the
/// same loop from different starts would emit different triangles for identical
/// input and break T-004.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Cycle {
    /// Which of §3.1's three kinds this came from.
    pub kind: CurveKind,
    /// The crossings, in cyclic order. `points[i]` joins `points[i + 1]`, and
    /// the last joins the first.
    pub points: Vec<FacePoint>,
}

impl Cycle {
    /// The loop's length `ℓ` — its segment count, which equals its point count.
    #[must_use]
    pub fn length(&self) -> usize {
        self.points.len()
    }

    /// Whether this is a corner cut: §3.2.1's `ℓ = 3` case.
    #[must_use]
    pub fn is_corner_cut(&self) -> bool {
        self.kind == CurveKind::Normal && self.points.len() == 3
    }
}

/// Order each closed curve into a cycle, dropping the open ones.
///
/// > If `γ` has a degree-1 vertex, it forms an open curve. Such curves are
/// > discarded, but their segments may still appear in neighboring tets as part
/// > of the mesh boundary.
///
/// The walk is over **segment indices**, not points, because a two-segment loop
/// joins the same pair of points twice and a point-keyed walk would see one
/// edge. The start is the lexicographically smallest segment and the first step
/// goes to its smaller endpoint, which makes the traversal a function of the
/// input alone.
#[must_use]
pub fn cycles(coords: &EdgeCoordinates) -> Vec<Cycle> {
    let mut out = Vec::new();
    for curve in curves(coords) {
        if curve.kind == CurveKind::Open {
            continue;
        }
        if let Some(points) = walk(&curve) {
            out.push(Cycle {
                kind: curve.kind,
                points,
            });
        }
    }
    out
}

/// Walk one closed curve's segments into cyclic order.
///
/// Returns `None` if the component is not a single cycle — every point degree 2
/// and every segment used exactly once. That is this crate's own invariant
/// rather than a caller error, and returning `None` keeps it out of the public
/// error type while still refusing to emit a wrong answer.
fn walk(curve: &Curve) -> Option<Vec<FacePoint>> {
    let segments: &[Segment] = &curve.segments;
    if segments.is_empty() {
        return None;
    }

    // Segments are sorted, so index 0 is the canonical start and `a <= b`
    // within it — the walk is determined from here on.
    let mut used = alloc::vec![false; segments.len()];
    let start = segments[0].a;
    let mut points = Vec::with_capacity(segments.len());
    let mut current = start;

    for _ in 0..segments.len() {
        points.push(current);
        let next = segments.iter().enumerate().find_map(|(i, s)| {
            if used[i] {
                return None;
            }
            if s.a == current {
                Some((i, s.b))
            } else if s.b == current {
                Some((i, s.a))
            } else {
                None
            }
        })?;
        used[next.0] = true;
        current = next.1;
    }

    // A closed cycle returns to where it began, having used every segment.
    if current == start && used.iter().all(|u| *u) {
        Some(points)
    } else {
        None
    }
}

/// The `(d₁, d₂)` pattern the loops left after the corner cuts are removed.
///
/// Property II (Theorem B.3) says the residual edge coordinates take the values
/// `d₁`, `d₂` and `d₁ + d₂` on the three **complementary pairs** of edges. This
/// type both extracts that and refuses to invent it: [`of`](Self::of) returns
/// `None` when the residual does not have the form, which is the signal that
/// something upstream is wrong rather than a case to paper over.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Pattern {
    /// The larger of the two quad multiplicities.
    pub d1: u32,
    /// The smaller, `0 ≤ d₂ ≤ d₁`.
    pub d2: u32,
}

impl Pattern {
    /// Read the pattern off residual edge coordinates.
    ///
    /// `None` if they are not `(d₁, d₁), (d₂, d₂), (d₁ + d₂, d₁ + d₂)` over the
    /// three complementary pairs — i.e. if Property II does not hold.
    #[must_use]
    pub fn of(residual: &EdgeCoordinates) -> Option<Self> {
        // The three complementary pairs are `e` and `5 - e` for e in 0..3.
        let mut pair = [0u32; 3];
        for (e, slot) in pair.iter_mut().enumerate() {
            let edge = e as u8;
            let (a, b) = (residual.edge(edge), residual.edge(complementary(edge)));
            if a != b {
                return None;
            }
            *slot = a;
        }
        pair.sort_unstable();
        // Sorted, the three values must be d₂, d₁, d₁ + d₂.
        if pair[0] + pair[1] != pair[2] {
            return None;
        }
        Some(Self {
            d1: pair[1],
            d2: pair[0],
        })
    }

    /// Whether the pattern is empty — no residual loops at all.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.d1 == 0 && self.d2 == 0
    }

    /// How many loops the pattern carries.
    ///
    /// Theorem B.4: *"the number of curves on the boundary of a tetrahedron with
    /// two non-zero diagonal cuts `d₁, d₂` is exactly `gcd(d₁, d₂)`."* With one
    /// of them zero the other counts directly, which `gcd(d, 0) = d` already
    /// gives.
    #[must_use]
    pub fn loop_count(&self) -> u32 {
        gcd(self.d1, self.d2)
    }

    /// The length `ℓ` every loop in the pattern has.
    ///
    /// Corollary B.6: *"the length of every component in a `(d₁, d₂)` pattern is
    /// exactly `4(d₁ + d₂) / gcd(d₁, d₂)`."* `None` for an empty pattern, which
    /// has no loops to have a length.
    #[must_use]
    pub fn loop_length(&self) -> Option<u32> {
        let g = self.loop_count();
        if g == 0 {
            return None;
        }
        Some(4 * (self.d1 + self.d2) / g)
    }
}

/// Greatest common divisor, with `gcd(d, 0) = d`.
const fn gcd(a: u32, b: u32) -> u32 {
    let (mut a, mut b) = (a, b);
    while b != 0 {
        let t = b;
        b = a % b;
        a = t;
    }
    a
}

/// The residual edge coordinates once every corner cut has been removed.
///
/// A corner cut is a length-3 loop, and its three points lie on the three edges
/// incident to one tet corner — so removing it subtracts one from each of those
/// three coordinates. Returns `None` if that would go negative, which cannot
/// happen for cycles derived from these very coordinates and is therefore this
/// crate's bug rather than a caller's.
#[must_use]
pub fn residual(coords: &EdgeCoordinates, cycles: &[Cycle]) -> Option<EdgeCoordinates> {
    let mut count = coords.count;
    for cut in cycles.iter().filter(|c| c.is_corner_cut()) {
        for point in &cut.points {
            let slot = &mut count[point.edge as usize];
            *slot = slot.checked_sub(1)?;
        }
    }
    Some(EdgeCoordinates { count })
}

/// Triangulate one tetrahedron's boundary curves.
///
/// Handles §3.2.1's first two cases — corner cuts (`ℓ = 3`) and quads
/// (`ℓ = 4`). Octagons, the single-loop case, subdivision and §3.2.2's
/// non-normal loops are not here yet, and [`Unfilled`] says which was reached
/// rather than emitting something wrong.
///
/// # Errors
///
/// [`NotFillable`] if the crossings violate [`TetCrossings::check`]'s contract.
pub fn fill<R: Real>(
    tet: &TetCrossings<'_, R>,
    out: &mut TetPatch<R>,
) -> Result<Unfilled, NotFillable> {
    tet.check()?;
    out.reset();

    let coords = tet.coordinates();
    let cycles = cycles(&coords);

    // Every crossing that any cycle uses becomes a vertex, indexed by its
    // FacePoint so the two tets sharing a face agree on which vertex is which.
    let mut keys: Vec<FacePoint> = cycles
        .iter()
        .flat_map(|c| c.points.iter().copied())
        .collect();
    keys.sort_unstable();
    keys.dedup();
    for key in &keys {
        match tet.position(*key) {
            Some(p) => out.positions.push(p),
            // A cycle named a crossing the crossing list does not have, which
            // means `coordinates()` and `along` disagree — impossible, since the
            // first is derived from the second.
            None => return Ok(Unfilled::Inconsistent),
        }
    }
    let index_of = |p: FacePoint| -> u32 {
        // `keys` is sorted and contains every point any cycle uses.
        keys.binary_search(&p).map_or(u32::MAX, |i| i as u32)
    };

    let mut unfilled = Unfilled::None;

    // Case 1 -- corner cuts. One triangle each, and they come first because
    // Properties I and II only hold once they are gone.
    for cut in cycles.iter().filter(|c| c.is_corner_cut()) {
        out.triangles.push([
            index_of(cut.points[0]),
            index_of(cut.points[1]),
            index_of(cut.points[2]),
        ]);
    }

    let Some(residual) = residual(&coords, &cycles) else {
        return Ok(Unfilled::Inconsistent);
    };
    let Some(pattern) = Pattern::of(&residual) else {
        // Property II did not hold. Either a non-normal configuration reached
        // here, or Theorem B.3 does not say what this code thinks it says.
        return Ok(Unfilled::NoPattern);
    };
    if pattern.is_empty() {
        return Ok(unfilled);
    }

    // Case 2 -- quads. "If the remaining loops are quads, then we split all
    // quads along the same, arbitrary diagonal, producing two triangles per
    // loop." The diagonal is arbitrary; *the same* is not, and it is what keeps
    // parallel quads from crossing each other. Splitting every quad from its
    // own first point -- which `walk` chose canonically -- is what makes it the
    // same one for every loop in a parallel family.
    for cycle in cycles.iter().filter(|c| c.kind == CurveKind::Normal) {
        match cycle.length() {
            3 => {}
            4 => {
                let v: Vec<u32> = cycle.points.iter().map(|p| index_of(*p)).collect();
                out.triangles.push([v[0], v[1], v[2]]);
                out.triangles.push([v[0], v[2], v[3]]);
            }
            _ => unfilled = unfilled.worst(Unfilled::NormalLoop),
        }
    }
    if cycles.iter().any(|c| c.kind == CurveKind::NonNormal) {
        unfilled = unfilled.worst(Unfilled::NonNormalLoop);
    }

    Ok(unfilled)
}

/// What [`fill`] met and did not yet handle.
///
/// Not an error and not a fallback: the triangles it *did* emit are correct, and
/// this names what was left on the floor. A caller that requires completeness
/// checks for [`Unfilled::None`].
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub enum Unfilled {
    /// Everything the configuration carried was triangulated.
    None,
    /// A normal loop longer than a quad — §3.2.1's octagon, single-loop and
    /// subdivision cases.
    NormalLoop,
    /// A non-normal loop — §3.2.2.
    NonNormalLoop,
    /// Property II failed on the residual coordinates.
    NoPattern,
    /// The curves and the crossing lists disagree. This crate's own bug.
    Inconsistent,
}

impl Unfilled {
    /// Keep whichever of two outcomes is further from complete.
    fn worst(self, other: Self) -> Self {
        if other > self { other } else { self }
    }
}

/// A tetrahedron has four faces, and every cycle lives on their union.
const _: () = assert!(TET_FACE_COUNT == 4);

#[cfg(test)]
mod tests;
