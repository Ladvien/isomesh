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

use super::coordinates::{EdgeCoordinates, TET_FACE_COUNT, TET_FACES, complementary};
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

/// The edge coordinates carried by the **normal** loops that are not corner
/// cuts — the residual Property II is a statement about.
///
/// Summed from those loops directly rather than by subtracting corner cuts from
/// the tet's own coordinates, and the difference is not cosmetic. Property II
/// (Theorem B.3) is stated over `Γ_normal`, and a tet carrying any *non-normal*
/// loop has points that belong to neither a corner cut nor a residual normal
/// loop. Subtracting only the corner cuts leaves those in, and the result then
/// fails Property II for a reason that has nothing to do with §3.2.1 — which is
/// exactly what happened on `e = (0, 0, 2, 0, 0, 0)`, where a single non-normal
/// loop made the whole configuration look patternless and stopped §3.2.2 from
/// ever running. Summing the residual normal loops cannot make that mistake.
#[must_use]
pub fn residual(cycles: &[Cycle]) -> EdgeCoordinates {
    let mut count = [0u32; TET_EDGE_COUNT];
    for cycle in cycles
        .iter()
        .filter(|c| c.kind == CurveKind::Normal && !c.is_corner_cut())
    {
        for point in &cycle.points {
            count[point.edge as usize] += 1;
        }
    }
    EdgeCoordinates { count }
}

/// Triangulate one tetrahedron's boundary curves.
///
/// Handles §3.2.1's corner cuts (`ℓ = 3`), quads (`ℓ = 4`) and octagons
/// (`ℓ = 8`). The single-loop case, subdivision and §3.2.2's non-normal loops
/// are not here yet, and [`Unfilled`] says which was reached rather than
/// emitting something wrong.
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

    // §3.2.2 -- non-normal loops. Run before §3.2.1's dispatch and independently
    // of it: the two sections partition Γ into Γ_nonnormal and Γ_normal, and a
    // configuration can carry both. Only the diagonal type is implemented -- its
    // spanning disk is in the tet *interior* and is a fan, the same primitive
    // the Steiner cases use. Corner and contractible loops want a disk built in
    // the tet *boundary*, by splitting each face along γ and labelling segments
    // inside/outside, which is a different construction and is not here.
    for cycle in cycles.iter().filter(|c| c.kind == CurveKind::NonNormal) {
        match cycle.non_normal_kind() {
            // "If γ is of diagonal type, we triangulate it by connecting each
            // of its segments to its center of mass a."
            Some(NonNormalKind::Diagonal) => {
                unfilled = unfilled.worst(fill_centroid_fan(cycle, &index_of, out));
            }
            // Contractible: the disk is built in the tet boundary, out of
            // crossings this tet already has.
            Some(NonNormalKind::Contractible) => {
                unfilled = unfilled.worst(fill_boundary_disk(&coords, cycle, &index_of, out));
            }
            // Corner type still wants the same disk plus §3.2.2's two extras --
            // omitting the vertex that coincides with the inside corner, and a
            // triangle at that corner.
            _ => unfilled = unfilled.worst(Unfilled::NonNormalLoop),
        }
    }

    // The dispatch is §3.2.1's own, and it is on the *pattern's* loop length
    // rather than on each cycle's, because Property I says they are the same
    // number and the case is a property of the configuration.
    let residual_loops: Vec<&Cycle> = cycles
        .iter()
        .filter(|c| c.kind == CurveKind::Normal && !c.is_corner_cut())
        .collect();
    let residual = residual(&cycles);
    let Some(pattern) = Pattern::of(&residual) else {
        // Property II did not hold over Γ_normal's own residual, which is the
        // set it is stated for. Either Theorem B.3 does not say what this code
        // thinks it says, or §3.1 produced loops that are not a (d₁, d₂)
        // pattern. Neither is a case to paper over.
        return Ok(unfilled.worst(Unfilled::NoPattern));
    };
    if pattern.is_empty() {
        return Ok(unfilled);
    }

    match pattern.loop_length() {
        // Case 2 -- quads. "If the remaining loops are quads, then we split all
        // quads along the same, arbitrary diagonal, producing two triangles per
        // loop." The diagonal is arbitrary; *the same* is not, and it is what
        // keeps parallel quads from crossing each other. Splitting every quad
        // from its own first point -- which `walk` chose canonically -- is what
        // makes it the same one for every loop in a parallel family.
        Some(4) => {
            for cycle in &residual_loops {
                let v: Vec<u32> = cycle.points.iter().map(|p| index_of(*p)).collect();
                out.triangles.push([v[0], v[1], v[2]]);
                out.triangles.push([v[0], v[2], v[3]]);
            }
        }
        // Case 3 -- octagons.
        Some(8) => {
            unfilled = unfilled.worst(fill_octagons(
                tet,
                &residual,
                pattern,
                &residual_loops,
                &index_of,
                out,
            ));
        }
        // Case 4 -- a single loop, longer than an octagon.
        Some(_) if pattern.loop_count() == 1 => {
            unfilled = unfilled.worst(fill_single_loop(&residual_loops, &index_of, out));
        }
        // Case 5 -- subdivision. Not implemented; named exactly.
        Some(_) => unfilled = unfilled.worst(Unfilled::Subdivision),
        None => {}
    }

    Ok(unfilled)
}

/// §3.2.2's spanning disk for a **contractible** loop, built in the tet
/// boundary.
///
/// > Otherwise, we build a piecewise linear disk contained mostly in the tet
/// > boundary, rather than its interior. … Finally, we emit any polygon `P`
/// > bound by "inside" segments and segments of `γ`.
///
/// Every vertex of that disk is already a crossing this tet has — the disk lies
/// *on* the boundary, so it introduces no new points, which is what makes this
/// case cheaper than any of the Steiner ones despite looking harder.
///
/// A contractible loop marks all four corners outside, so an arc touching a
/// corner is always outside and no qualifying region can contain one. That is
/// asserted rather than assumed: a `Corner` node reaching the emit is
/// [`Unfilled::Inconsistent`], because for a *corner-type* loop it would instead
/// need §3.2.2's omission rule and a corner triangle, and silently fanning over
/// it would produce a plausible wrong disk.
///
/// The triangles are expected to be degenerate wherever a scoop bounds the
/// region — see V-21. A-014d insets them; this stage is graded on connectivity.
fn fill_boundary_disk<R: Real>(
    coords: &EdgeCoordinates,
    cycle: &Cycle,
    index_of: &impl Fn(FacePoint) -> u32,
    out: &mut TetPatch<R>,
) -> Unfilled {
    if cycle.non_normal_kind() != Some(NonNormalKind::Contractible) {
        return Unfilled::NonNormalLoop;
    }

    for face in 0..TET_FACE_COUNT as u8 {
        let Some(regions) = face_regions(face, coords, cycle) else {
            return Unfilled::Inconsistent;
        };
        for region in regions.iter().filter(|r| r.is_inside()) {
            let mut corner = Vec::with_capacity(region.node.len());
            for node in &region.node {
                match node {
                    Node::Crossing(p) => corner.push(index_of(*p)),
                    Node::Corner(_) => return Unfilled::Inconsistent,
                }
            }
            // A fan from the region's first vertex. The region is a polygon of
            // one planar face cut by chords, so it is convex and a fan is a
            // valid triangulation of it.
            for k in 1..corner.len().saturating_sub(1) {
                out.triangles.push([corner[0], corner[k], corner[k + 1]]);
            }
        }
    }
    Unfilled::None
}

/// Fan a loop around the centre of mass of its own vertices.
///
/// Shared by §3.2.1's single-loop case and §3.2.2's diagonal type, which give
/// the same instruction for different reasons — the first because any point of
/// the loop's convex hull works and a centroid is one, the second because the
/// paper names the centre of mass directly.
fn fill_centroid_fan<R: Real>(
    cycle: &Cycle,
    index_of: &impl Fn(FacePoint) -> u32,
    out: &mut TetPatch<R>,
) -> Unfilled {
    let n = cycle.points.len();
    if n < 3 {
        return Unfilled::Inconsistent;
    }

    // Summed in the cycle's canonical order, so the rounding is identical run to
    // run and T-004 has nothing to catch.
    let mut sum = [R::ZERO; 3];
    for point in &cycle.points {
        let Some(index) = usize::try_from(index_of(*point))
            .ok()
            .filter(|i| *i < out.positions.len())
        else {
            return Unfilled::Inconsistent;
        };
        let p = out.positions[index];
        for (slot, value) in sum.iter_mut().zip(p.iter()) {
            *slot += *value;
        }
    }
    let scale = R::ONE / R::from_f64(n as f64);
    let steiner = out.positions.len() as u32;
    out.positions
        .push([sum[0] * scale, sum[1] * scale, sum[2] * scale]);
    fan(cycle, steiner, index_of, out);
    Unfilled::None
}

/// §3.2.2's three kinds of non-normal loop.
///
/// > Viewing a tetrahedron as a topological sphere punctured at its four
/// > vertices, each curve `γ ∈ Γ_nonnormal` has one of three types (no matter
/// > how much it "spirals" around the tet).
///
/// The "no matter how much it spirals" is the useful part: the type is a
/// homotopy class, so it survives any amount of winding and is decidable from
/// parities alone.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum NonNormalKind {
    /// *"Does not separate any vertices."*
    Contractible,
    /// *"Separates two vertices from the other two."*
    Diagonal,
    /// *"Separates one vertex from the other three."*
    Corner,
}

impl Cycle {
    /// Which of §3.2.2's three types a non-normal loop is.
    ///
    /// > We evaluate the loop type using a parity bit `b_ij := mod(e_ij^γ, 2)`
    /// > for each edge `ij`, where `e_ij^γ` are edge coordinates for `γ` alone.
    /// > This value is odd if `γ` separates vertices `i` and `j`, and even if
    /// > they belong to the same connected component of the tet boundary.
    /// > Letting `p := b₀₁ + b₀₂ + b₀₃`, `γ` is then contractible if `p = 0`, is
    /// > of diagonal type if `p = 2`, and is of corner type if `p = 1` or
    /// > `p = 3`.
    ///
    /// `None` for a normal loop, which has no such type.
    ///
    /// Note the asymmetry that makes this well defined: `p` is summed over the
    /// three edges at **corner 0** only, not over all six. A loop separating
    /// corner 0 from the rest gives `p = 3`, one separating some *other* single
    /// corner gives `p = 1`, and both are the corner type — which is why the
    /// two odd values collapse to one answer rather than naming two types.
    #[must_use]
    pub fn non_normal_kind(&self) -> Option<NonNormalKind> {
        if self.kind != CurveKind::NonNormal {
            return None;
        }
        let mut own = [0u32; TET_EDGE_COUNT];
        for point in &self.points {
            own[point.edge as usize] += 1;
        }
        // Edges 0, 1, 2 are (0,1), (0,2) and (0,3) in TET_EDGES' lexicographic
        // order — the three at corner 0.
        let p: u32 = own[0] % 2 + own[1] % 2 + own[2] % 2;
        Some(match p {
            0 => NonNormalKind::Contractible,
            2 => NonNormalKind::Diagonal,
            _ => NonNormalKind::Corner,
        })
    }
}

/// One step around a face's boundary, between two consecutive nodes.
///
/// A face of the tet is a triangle whose boundary carries the loop's crossings.
/// Walking it gives an alternating sequence of nodes and arcs; the arcs are what
/// a region is bounded by.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Arc {
    /// A piece of one of the face's three edges, with the side §3.2.2's
    /// labelling gives it.
    Edge {
        /// Which tet edge it lies on.
        edge: u8,
        /// Which piece along that edge, counting from its lower corner.
        piece: usize,
        /// Inside or outside the loop.
        side: Side,
    },
    /// A segment of `γ` crossing the face's interior. A **scoop** — both
    /// endpoints on one edge — is one of these too, and is exactly the case
    /// V-21 warns realises with zero area until A-014d insets it.
    Chord,
}

/// A point on a face's boundary walk.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Node {
    /// A corner of the tetrahedron.
    Corner(u8),
    /// A crossing of the loop with one of the face's edges.
    Crossing(FacePoint),
}

/// A region of one face, cut out by the loop's segments.
///
/// Its boundary, as arcs in order. Emitted by §3.2.2 only when every [`Arc::Edge`]
/// in it is [`Side::Inside`].
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Region {
    /// The arcs bounding it, in order.
    pub arc: Vec<Arc>,
    /// The nodes between them: `node[i]` precedes `arc[i]`.
    pub node: Vec<Node>,
}

impl Region {
    /// Whether §3.2.2 emits this region: *"we emit any polygon `P` bound by
    /// 'inside' segments and segments of `γ`."*
    ///
    /// A region touching a tet corner always fails this, because every corner is
    /// outside for a contractible loop and all but one for a corner loop — which
    /// is what stops the disk from swallowing a vertex.
    #[must_use]
    pub fn is_inside(&self) -> bool {
        self.arc.iter().all(|a| match a {
            Arc::Edge { side, .. } => *side == Side::Inside,
            Arc::Chord => true,
        })
    }
}

/// Which side of a non-normal loop a piece of the tet boundary is on.
///
/// > Since `γ` is a closed simple curve, it partitions the tet boundary into
/// > two pieces. To define spanning disks, we must classify these pieces as
/// > "inside" or "outside" `γ`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Side {
    /// The piece the spanning disk is built over.
    Inside,
    /// The other one.
    Outside,
}

impl Side {
    /// The other side. Crossing `γ` once swaps them, which is what makes the
    /// labels along an edge alternate.
    #[must_use]
    pub const fn flipped(self) -> Self {
        match self {
            Self::Inside => Self::Outside,
            Self::Outside => Self::Inside,
        }
    }
}

impl Cycle {
    /// The corner a corner-type loop separates from the other three.
    ///
    /// Read off the parity bits the classification already computes: `b_ij` is
    /// odd exactly when `γ` separates `i` from `j`, so the distinguished corner
    /// is the one whose three incident edges are all odd. `None` for any loop
    /// that is not corner type.
    #[must_use]
    pub fn distinguished_corner(&self) -> Option<u8> {
        if self.non_normal_kind()? != NonNormalKind::Corner {
            return None;
        }
        let mut own = [0u32; TET_EDGE_COUNT];
        for point in &self.points {
            own[point.edge as usize] += 1;
        }
        (0..4u8).find(|corner| {
            (0..4u8).filter(|other| other != corner).all(|other| {
                own[super::coordinates::edge_between(*corner, other) as usize] % 2 == 1
            })
        })
    }

    /// Which side of the loop each of the tet's four corners is on.
    ///
    /// > When `γ` is of corner type, we mark the distinguished vertex as inside
    /// > and all other vertices as outside. When `γ` is contractible, all
    /// > vertices are "outside," and when `γ` is diagonal we do not require an
    /// > inside/outside distinction.
    ///
    /// `None` for a diagonal loop — which is not a gap but the paper's own
    /// answer, and the reason the diagonal case can be triangulated without one.
    #[must_use]
    pub fn corner_sides(&self) -> Option<[Side; 4]> {
        match self.non_normal_kind()? {
            NonNormalKind::Diagonal => None,
            NonNormalKind::Contractible => Some([Side::Outside; 4]),
            NonNormalKind::Corner => {
                let distinguished = self.distinguished_corner()?;
                let mut side = [Side::Outside; 4];
                side[distinguished as usize] = Side::Inside;
                Some(side)
            }
        }
    }

    /// The side of each piece an edge is cut into, from its lower corner.
    ///
    /// > Along each edge `ij` of `σ`, the first segment inherits the label of
    /// > vertex `i`; subsequent segments alternate between "inside" and
    /// > "outside" labels.
    ///
    /// An edge carrying `n` of this loop's crossings is cut into `n + 1` pieces,
    /// so that is the length of the result. `i` is the edge's **lower-numbered**
    /// corner, which makes the labelling agree between the two faces sharing the
    /// edge without either consulting the other — the same mechanism that gives
    /// §3.1 its conformity (M-79).
    ///
    /// `None` for a diagonal loop, which has no sides.
    #[must_use]
    pub fn edge_sides(&self, edge: u8) -> Option<Vec<Side>> {
        let sides = self.corner_sides()?;
        let [lo, _hi] = TET_EDGES[edge as usize];
        let crossings = self.points.iter().filter(|p| p.edge == edge).count();

        let mut out = Vec::with_capacity(crossings + 1);
        let mut side = sides[lo as usize];
        out.push(side);
        for _ in 0..crossings {
            side = side.flipped();
            out.push(side);
        }
        Some(out)
    }
}

/// Split one face along a loop's segments — §3.2.2's `σ \ γ`.
///
/// > We first split each triangle `σ` of the tet boundary along the segments of
/// > `γ`, yielding a collection of planar polygons `σ \ γ`.
///
/// Built **combinatorially**, from the boundary walk and the pairing, never from
/// straight-line geometry. V-21 is why: a scoop's chord lies *along* its own
/// edge, so a geometric build collapses it to zero area before A-014d can inset
/// it, and the region it bounds would disappear rather than be emitted
/// degenerate as the paper intends.
///
/// # The peel
///
/// `γ` is simple, so its chords on a face do not cross, and non-crossing chords
/// on a disk nest. So there is always an **innermost** chord — one whose two
/// endpoints have no other chord endpoint between them along one side — and the
/// arcs spanning that gap, plus the chord, bound a region with nothing inside
/// it. Emit it, replace the whole span by the chord itself, and repeat. What
/// remains when no chords are left is the last region.
///
/// That the peel always finds an innermost chord is not assumed: if it cannot,
/// the chords were not non-crossing, which would mean §3.1 produced a
/// self-crossing curve, and this returns `None` rather than inventing a
/// decomposition.
///
/// Returns `None` for a diagonal loop, which has no sides and therefore no
/// inside/outside decomposition to make.
#[must_use]
pub fn face_regions(face: u8, coords: &EdgeCoordinates, cycle: &Cycle) -> Option<Vec<Region>> {
    let f = TET_FACES[face as usize];
    let sides: [Vec<Side>; 3] = [
        cycle.edge_sides(f.edge[0])?,
        cycle.edge_sides(f.edge[1])?,
        cycle.edge_sides(f.edge[2])?,
    ];

    // Walk the face boundary: corner k, then this loop's crossings along
    // edge k in the direction corner k -> corner k+1, then corner k+1, ...
    // `piece` indexes the sub-segments of an edge from its *lower* corner,
    // which is the direction `edge_sides` labels in, so the walk converts once
    // here and nowhere else.
    let mut node: Vec<Node> = Vec::new();
    let mut arc: Vec<Arc> = Vec::new();
    for (k, side_of) in sides.iter().enumerate() {
        let edge = f.edge[k];
        let [lo, _hi] = TET_EDGES[edge as usize];
        let forward = f.corner[k] == lo;
        let mut index: Vec<u32> = cycle
            .points
            .iter()
            .filter(|p| p.edge == edge)
            .map(|p| p.index)
            .collect();
        index.sort_unstable();
        if !forward {
            index.reverse();
        }

        node.push(Node::Corner(f.corner[k]));
        for (step, crossing) in index.iter().enumerate() {
            // Piece `step` of the walk is piece `step` from the lower corner
            // when walking forward, and counts back from the far end otherwise.
            let piece = if forward { step } else { index.len() - step };
            arc.push(Arc::Edge {
                edge,
                piece,
                side: *side_of.get(piece)?,
            });
            node.push(Node::Crossing(FacePoint {
                edge,
                index: *crossing,
            }));
        }
        let piece = if forward { index.len() } else { 0 };
        arc.push(Arc::Edge {
            edge,
            piece,
            side: *side_of.get(piece)?,
        });
    }

    // Which walk positions are joined by a chord.
    //
    // Face assignment comes from `face_segments`, not from "both endpoints lie
    // on an edge this face has". Those differ exactly for a **scoop**: its two
    // endpoints are on one edge, and an edge belongs to *two* faces, so the
    // weaker test hands the same scoop to both of them and a crossing ends up
    // bounding three regions instead of two. §3.1 already knows the answer,
    // because it builds segments per face in the first place.
    let mut partner: Vec<Option<usize>> = alloc::vec![None; node.len()];
    let position = |p: FacePoint| node.iter().position(|n| *n == Node::Crossing(p));
    // This loop's own segments, as endpoint pairs in the same sorted order a
    // `Segment` stores them, so the two can be compared without reaching for
    // `Segment`'s private constructor.
    let mine: Vec<(FacePoint, FacePoint)> = {
        let mut s: Vec<(FacePoint, FacePoint)> = (0..cycle.points.len())
            .map(|k| {
                let (a, b) = (cycle.points[k], cycle.points[(k + 1) % cycle.points.len()]);
                if a <= b { (a, b) } else { (b, a) }
            })
            .collect();
        s.sort_unstable();
        s
    };
    for segment in super::curves::face_segments(face, coords) {
        if mine.binary_search(&(segment.a, segment.b)).is_err() {
            continue;
        }
        if let (Some(i), Some(j)) = (position(segment.a), position(segment.b)) {
            partner[i] = Some(j);
            partner[j] = Some(i);
        }
    }

    Some(peel(node, arc, partner))
}

/// Repeatedly cut off the innermost chord's region. See [`face_regions`].
fn peel(mut node: Vec<Node>, mut arc: Vec<Arc>, mut partner: Vec<Option<usize>>) -> Vec<Region> {
    let mut out = Vec::new();

    loop {
        let n = node.len();
        // An innermost chord: `i -> j` forward with no other chord endpoint
        // strictly between.
        let innermost = (0..n).find_map(|i| {
            let j = partner[i]?;
            let span = (j + n - i) % n;
            if span == 0 {
                return None;
            }
            let clear = (1..span).all(|s| partner[(i + s) % n].is_none());
            clear.then_some((i, j, span))
        });

        let Some((i, j, span)) = innermost else {
            break;
        };

        // The region: the arcs from `i` forward to `j`, closed by the chord.
        let mut region = Region {
            arc: Vec::with_capacity(span + 1),
            node: Vec::with_capacity(span + 1),
        };
        for s in 0..span {
            region.node.push(node[(i + s) % n]);
            region.arc.push(arc[(i + s) % n]);
        }
        region.node.push(node[j]);
        region.arc.push(Arc::Chord);
        out.push(region);

        // Collapse: the span becomes the chord, and its interior nodes go.
        let mut kept_node = Vec::with_capacity(n - span + 1);
        let mut kept_arc = Vec::with_capacity(n - span + 1);
        let mut kept_partner = Vec::with_capacity(n - span + 1);
        kept_node.push(node[i]);
        kept_arc.push(Arc::Chord);
        kept_partner.push(None);
        for s in span..n {
            let at = (i + s) % n;
            kept_node.push(node[at]);
            kept_arc.push(arc[at]);
            kept_partner.push(partner[at]);
        }
        // The two chord endpoints have been consumed; nothing else moved
        // relative to them, so re-derive the pairing by identity.
        let remap: Vec<Option<usize>> = kept_partner
            .iter()
            .map(|p| {
                p.and_then(|old| {
                    let original = |at: usize| (i + at) % n;
                    (0..kept_node.len()).find(|k| *k > 0 && original(span + k - 1) == old)
                })
            })
            .collect();
        node = kept_node;
        arc = kept_arc;
        partner = remap;
        let _ = j;
    }

    // Whatever is left, once no chord remains, is the final region.
    if !arc.is_empty() {
        out.push(Region { arc, node });
    }
    out
}

/// §3.2.1 case (3) — the subdivision stencil, as a labelling and four tets.
///
/// > Noting Property (II) above, let `i, j, k, l` be vertices of the tetrahedron
/// > such that `e_ij = e_kl = d₁`, `e_ik = e_jl = d₂`, and
/// > `e_il = e_jk = d₁ + d₂`. … connecting its vertices to the center of mass
/// > `a` (Figure 13), and assign edge coordinates
/// > `e_ai = 2d₂, e_aj = d₁, e_ak = d₂, e_al = d₁ − d₂`.
/// > We then recursively process each of the four new tets.
///
/// The stencil is **not symmetric in the four corners** — `e_ai` is `2d₂` while
/// `e_aj` is `d₁` — so which corner is called `i` decides the whole
/// subdivision, and [`label`](Subdivision::label) searches for the labelling
/// Property II guarantees rather than assuming corner order supplies it.
///
/// This type is the combinatorics only: the labelling and the four sub-tets'
/// edge coordinates. Placing the crossings on the four new edges and recursing
/// is the geometric half, and is not here yet — see [`Unfilled::Subdivision`].
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Subdivision {
    /// The corners in the stencil's own order: `[i, j, k, l]`.
    pub corner: [u8; 4],
    /// The pattern the stencil was read from.
    pub pattern: Pattern,
}

impl Subdivision {
    /// Find the `i, j, k, l` labelling Property II promises.
    ///
    /// Returns `None` when no labelling of the four corners puts `d₁` on
    /// `(ij, kl)`, `d₂` on `(ik, jl)` and `d₁ + d₂` on `(il, jk)` — which for a
    /// residual satisfying Property II cannot happen, so `None` is a signal that
    /// the residual is not what it was taken to be.
    ///
    /// `d₁ = d₂` makes several labellings equivalent; the search returns the
    /// lexicographically smallest, so the output is deterministic.
    #[must_use]
    pub fn label(residual: &EdgeCoordinates, pattern: Pattern) -> Option<Self> {
        let (d1, d2) = (pattern.d1, pattern.d2);
        let e = |a: u8, b: u8| residual.edge(super::coordinates::edge_between(a, b));
        for i in 0..4u8 {
            for j in 0..4u8 {
                if j == i {
                    continue;
                }
                for k in 0..4u8 {
                    if k == i || k == j {
                        continue;
                    }
                    let l = 6 - i - j - k;
                    if e(i, j) == d1
                        && e(k, l) == d1
                        && e(i, k) == d2
                        && e(j, l) == d2
                        && e(i, l) == d1 + d2
                        && e(j, k) == d1 + d2
                    {
                        return Some(Self {
                            corner: [i, j, k, l],
                            pattern,
                        });
                    }
                }
            }
        }
        None
    }

    /// The four new edge coordinates, `[e_ai, e_aj, e_ak, e_al]`.
    ///
    /// `d₁ − d₂` is non-negative because [`Pattern`] orders them, so the
    /// subtraction cannot wrap.
    #[must_use]
    pub fn spoke(&self) -> [u32; 4] {
        let (d1, d2) = (self.pattern.d1, self.pattern.d2);
        [2 * d2, d1, d2, d1 - d2]
    }

    /// The edge coordinates of the sub-tet opposite corner `which` of `[i,j,k,l]`
    /// — the one that replaces that corner with the centre of mass `a`.
    ///
    /// Corners of the returned tet are in the order `a` first, then the three
    /// surviving corners in their stencil order, so index 0 is always `a`.
    #[must_use]
    pub fn sub_tet(&self, residual: &EdgeCoordinates, which: usize) -> Option<EdgeCoordinates> {
        if which >= 4 {
            return None;
        }
        let spoke = self.spoke();
        // Local corner 0 is `a`; locals 1..3 are the stencil corners that are
        // not `which`, in stencil order.
        let mut kept = [0u8; 3];
        let mut slot = 0;
        for (s, corner) in self.corner.iter().enumerate() {
            if s != which {
                kept[slot] = *corner;
                slot += 1;
            }
        }
        let mut spoke_of = [0u32; 3];
        for (s, value) in spoke_of.iter_mut().enumerate() {
            let stencil = self.corner.iter().position(|c| *c == kept[s])?;
            *value = spoke[stencil];
        }

        let mut count = [0u32; TET_EDGE_COUNT];
        for (local_lo, local_hi) in (0..4u8).flat_map(|a| (a + 1..4u8).map(move |b| (a, b))) {
            let edge = super::coordinates::edge_between(local_lo, local_hi) as usize;
            count[edge] = if local_lo == 0 {
                // A spoke, `a` to a surviving corner.
                spoke_of[local_hi as usize - 1]
            } else {
                // An original edge of the parent tet.
                residual.edge(super::coordinates::edge_between(
                    kept[local_lo as usize - 1],
                    kept[local_hi as usize - 1],
                ))
            };
        }
        Some(EdgeCoordinates { count })
    }
}

/// §3.2.1 case (2) — one loop, one Steiner point.
///
/// > If we have just `m = 1` loop, we insert an additional Steiner point `x` at
/// > any point in the convex hull of the loop vertices, and triangulate the loop
/// > by connecting each of its edges to `x`. In practice we let `x` be the
/// > center of mass.
///
/// *"Any point in the convex hull"* is a genuine degree of freedom rather than
/// an under-specification: the loop is a closed curve on the tet boundary, so
/// every point of its hull sees every edge, and any of them fans without
/// self-overlap. The centre of mass is in the hull for the same reason a
/// centroid always is, so following the paper's practice costs nothing and
/// invents nothing.
///
/// Reached only for `ℓ > 8`, since `ℓ = 4` is the quad case and `ℓ = 8` is the
/// octagon case, which the paper lists first and which therefore takes the
/// `m = 1, ℓ = 8` overlap.
fn fill_single_loop<R: Real>(
    loops: &[&Cycle],
    index_of: &impl Fn(FacePoint) -> u32,
    out: &mut TetPatch<R>,
) -> Unfilled {
    let [cycle] = loops else {
        return Unfilled::Inconsistent;
    };
    fill_centroid_fan(cycle, index_of, out)
}

/// Connect every edge of a loop to one point — §3.2.1's triangulation primitive.
///
/// Used by both Steiner cases. Emits one triangle per loop edge, in cycle order,
/// so the winding is whatever the cycle's own orientation is and is consistent
/// across every loop in a tet.
fn fan<R: Real>(
    cycle: &Cycle,
    steiner: u32,
    index_of: &impl Fn(FacePoint) -> u32,
    out: &mut TetPatch<R>,
) {
    let n = cycle.points.len();
    for k in 0..n {
        let a = index_of(cycle.points[k]);
        let b = index_of(cycle.points[(k + 1) % n]);
        out.triangles.push([a, b, steiner]);
    }
}

/// §3.2.1 case (1) — octagons, and their Steiner points.
///
/// > For `ℓ = 8` (octagons), some pair of opposite edges `e, e'` has `2m`
/// > intersections, and the remaining edges each have `m` (Equation 4). To
/// > triangulate the octagons, we place `m` Steiner points `x₀, …, x_{m-1}`
/// > uniformly along the oriented segment from the midpoint of `e` to the
/// > midpoint of `e'`. Pairs of intersections on `e` are then connected to
/// > consecutive Steiner points, from the innermost to outermost pair. I.e., if
/// > we enumerate intersections along `e` as `p₀, …, p_{2m-1}` (with either
/// > orientation), then all points on the loop passing through `p_{m+i}` are
/// > connected to Steiner point `xᵢ`.
///
/// # Why `ℓ = 8` *is* `d₁ = d₂`
///
/// Corollary B.6 makes the case a statement about the pattern rather than a
/// separate condition to test. With `g = gcd(d₁, d₂)`, `d₁ = ga` and `d₂ = gb`
/// where `gcd(a, b) = 1`, the length `4(d₁ + d₂)/g = 8` forces `a + b = 2`, and
/// coprimality rules out `(2, 0)` — so `a = b = 1` and `d₁ = d₂ = g = m`. The
/// three complementary pairs then carry `m`, `m` and `2m`, which is exactly the
/// paper's *"some pair of opposite edges has 2m … and the remaining edges each
/// have m"*, and the `2m` pair is unique whenever `m > 0`.
///
/// # The one place a choice is made
///
/// The construction is not symmetric in `e` and `e'` — `x₀` sits nearest the
/// midpoint of whichever edge is called `e`. That choice is invisible outside
/// the tet: the octagons' Steiner points are interior, so no neighbouring tet
/// can see them and conformity cannot depend on them. It is fixed to the
/// lower-numbered edge of the pair purely so the output is deterministic.
///
/// Spacing is `(i + 1)/(m + 1)` along the segment, which is the reading of
/// *"uniformly"* that keeps every Steiner point strictly interior. The
/// intersection-free argument depends on the Steiner points being **ordered**
/// along the segment, not on the particular spacing.
fn fill_octagons<R: Real>(
    tet: &TetCrossings<'_, R>,
    residual: &EdgeCoordinates,
    pattern: Pattern,
    loops: &[&Cycle],
    index_of: &impl Fn(FacePoint) -> u32,
    out: &mut TetPatch<R>,
) -> Unfilled {
    let m = pattern.loop_count();
    if m == 0 || pattern.d1 != pattern.d2 || loops.len() as u32 != m {
        return Unfilled::Inconsistent;
    }

    // The unique complementary pair carrying 2m, lower index first.
    let Some(e) = (0..TET_EDGE_COUNT as u8).find(|e| {
        *e < complementary(*e)
            && residual.edge(*e) == 2 * m
            && residual.edge(complementary(*e)) == 2 * m
    }) else {
        return Unfilled::Inconsistent;
    };
    let opposite = complementary(e);

    // p₀ … p_{2m-1}: the residual crossings on `e`, which are the ones these
    // loops actually use. Deriving the enumeration from the loops rather than
    // from an index arithmetic on the full crossing list is what keeps this
    // correct when corner cuts have already claimed points at either end of the
    // edge -- their points are simply not in this list.
    let mut p: Vec<u32> = loops
        .iter()
        .flat_map(|c| c.points.iter())
        .filter(|q| q.edge == e)
        .map(|q| q.index)
        .collect();
    p.sort_unstable();
    p.dedup();
    if p.len() as u32 != 2 * m {
        return Unfilled::Inconsistent;
    }

    let half = R::from_f64(0.5);
    let midpoint = |edge: u8| -> [R; 3] {
        let [lo, hi] = TET_EDGES[edge as usize];
        let (a, b) = (tet.corners[lo as usize], tet.corners[hi as usize]);
        [
            (a[0] + b[0]) * half,
            (a[1] + b[1]) * half,
            (a[2] + b[2]) * half,
        ]
    };
    let (from, to) = (midpoint(e), midpoint(opposite));

    // The Steiner points, appended after every crossing vertex.
    let steiner_base = out.positions.len() as u32;
    let step = R::ONE / R::from_f64(f64::from(m) + 1.0);
    for i in 0..m {
        let t = step * R::from_f64(f64::from(i) + 1.0);
        out.positions.push([
            from[0] + (to[0] - from[0]) * t,
            from[1] + (to[1] - from[1]) * t,
            from[2] + (to[2] - from[2]) * t,
        ]);
    }

    // Assign each loop its Steiner point, and check the nesting the assignment
    // rests on before emitting anything.
    let mut assigned = Vec::with_capacity(loops.len());
    for cycle in loops {
        if cycle.length() != 8 {
            return Unfilled::Inconsistent;
        }
        let mut rank: Vec<u32> = cycle
            .points
            .iter()
            .filter(|q| q.edge == e)
            .filter_map(|q| p.binary_search(&q.index).ok().map(|r| r as u32))
            .collect();
        rank.sort_unstable();
        // Each loop crosses `e` exactly twice, and the pairs are nested about
        // the edge's middle: `p_j` pairs with `p_{2m-1-j}`. If that fails, the
        // "innermost to outermost" assignment has no meaning and guessing one
        // would be exactly the invented convention rule 5 forbids.
        let [low, high] = match rank[..] {
            [low, high] => [low, high],
            _ => return Unfilled::Inconsistent,
        };
        if low + high != 2 * m - 1 || high < m {
            return Unfilled::Inconsistent;
        }
        assigned.push((cycle, steiner_base + (high - m)));
    }

    // Fan each octagon around its Steiner point: one triangle per loop edge.
    for (cycle, steiner) in assigned {
        fan(cycle, steiner, index_of, out);
    }

    Unfilled::None
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
    /// §3.2.1 case (2) — one loop longer than an octagon, wanting a single
    /// Steiner point in the convex hull of its vertices.
    SingleLoop,
    /// §3.2.1 case (3) — several loops longer than octagons, wanting the
    /// Figure-13 subdivision stencil and a recursive call per new tet.
    Subdivision,
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
