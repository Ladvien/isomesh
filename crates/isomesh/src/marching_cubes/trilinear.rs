//! The body saddles of the trilinear interpolant, from one quadratic.
//!
//! Grosso, R., *Construction of Topologically Correct and Manifold Isosurfaces*,
//! Computer Graphics Forum 35(5), pp. 187–196 (`10.1111/cgf.12975`), §4; and
//! Grosso, R., *An Asymptotic Decider for Robust and Topologically Correct
//! Triangulation of Isosurfaces*, CGI '17 (`10.1145/3095140.3095179`), §3.
//!
//! # What this decides, and why it is not [`super::interior`]
//!
//! [`super::ambiguity`] resolves an ambiguous *face*. [`super::interior`] answers
//! Chernyaev's question about one *pair* of opposite faces: are the outside
//! regions those two faces carry joined through the cell? This module answers a
//! different and larger question — **where are the trilinear interpolant's saddle
//! points inside the cell, and how many of them are there** — from which the
//! cell's whole topology follows.
//!
//! The two are independent constructions of the same geometry, which is the point:
//! `the_body_saddle_heights_agree_with_the_swept_saddle_roots` checks one against
//! the other, and neither shares a line of arithmetic with the other.
//!
//! # The construction
//!
//! Restricted to a cell, the interpolant is
//!
//! ```text
//! F(u,v,w) = Σ f_i · (u or 1−u)(v or 1−v)(w or 1−w)
//! ```
//!
//! with `f_i` this crate's corner values and `(u,v,w) ∈ [0,1]³`. **Grosso's corner
//! numbering and this crate's are the same numbering** — his `v0 = (0,0,0)`,
//! `v1 = (1,0,0)`, `v2 = (0,1,0)`, `v3 = (1,1,0)`, `v4 = (0,0,1)` … is exactly
//! `cube.rs`'s "corner `i` sits at `(i&1, (i>>1)&1, (i>>2)&1)`". That is a
//! coincidence rather than a design, so `grosso_corner_numbering_is_ours` pins it;
//! if it ever stops being true every formula below needs a permutation.
//!
//! On a face the interpolant restricts to a bilinear function whose level set is a
//! rectangular hyperbola. For the pair of opposite faces `w = 0` and `w = 1`, the
//! projections of the two hyperbolas onto the `(u,v)` plane meet where a quadratic
//! in `u` vanishes. Where they meet, the line joining the two faces at that
//! `(u,v)` lies **entirely on the level set** (Grosso eq. 5), because the
//! interpolant is linear along it and equal to the isovalue at both ends. Two such
//! lines from different face pairs cross at a **body saddle**.
//!
//! # One quadratic, not three
//!
//! The paper presents "three quadratic equations", one per pair of opposite faces.
//! Solving three is unnecessary and the authors' own implementation does not:
//!
//! > It is enough to compute a pair of solutions for one face. The other solutions
//! > are obtained by evaluating the equations for the common variable.
//!
//! Both saddles share their `u` coordinate set. Given `u`, the `v` coordinate is
//! where the `w = 0` face's hyperbola sits at that `u`, and the `w` coordinate is
//! where the `v = 0` face's hyperbola sits — each a **linear** solve. So the whole
//! classification is one quadratic and four linear interpolations, and the three
//! coordinate sets cannot disagree with each other the way three separate
//! quadratics could.
//!
//! # The coefficients, derived here rather than transcribed
//!
//! Rule 5 forbids guessing a published formula; re-deriving one is how you avoid
//! having to. Writing `g₀(u) = f₀(1−u) + f₁u` and `g₁(u) = f₂(1−u) + f₃u` for the
//! `w = 0` face and `g̃₀`, `g̃₁` from `f₄, f₅` and `f₆, f₇` for `w = 1`, the level
//! set on each face is `v = (i₀ − g₀)/(g₁ − g₀)` and `v = (i₀ − g̃₀)/(g̃₁ − g̃₀)`.
//! Setting them equal and clearing denominators:
//!
//! ```text
//! (i₀ − g₀)(g̃₁ − g̃₀) − (i₀ − g̃₀)(g₁ − g₀) = 0
//! ```
//!
//! Each factor is linear in `u`:
//!
//! ```text
//! g₁ − g₀  = (f₂ − f₀) + u(f₀ + f₃ − f₁ − f₂)
//! g̃₁ − g̃₀ = (f₆ − f₄) + u(f₄ + f₇ − f₅ − f₆)
//! i₀ − g₀  = (i₀ − f₀) − u(f₁ − f₀)
//! i₀ − g̃₀ = (i₀ − f₄) − u(f₅ − f₄)
//! ```
//!
//! so collecting powers of `u` gives `a·u² + b·u + c` with
//!
//! ```text
//! a = (f₅ − f₄)(f₀ + f₃ − f₁ − f₂) − (f₁ − f₀)(f₄ + f₇ − f₅ − f₆)
//! b = (i₀ − f₀)(f₄ + f₇ − f₅ − f₆) − (f₁ − f₀)(f₆ − f₄)
//!   − (i₀ − f₄)(f₀ + f₃ − f₁ − f₂) + (f₅ − f₄)(f₂ − f₀)
//! c = (i₀ − f₀)(f₆ − f₄) − (i₀ − f₄)(f₂ − f₀)
//! ```
//!
//! Those agree term for term with the three the paper prints **and** with the
//! three the authors' implementation computes — a three-way agreement between a
//! derivation, a paper and a program, where V-24 had two (V-30).
//!
//! This crate extracts the **zero** level set and has no isovalue parameter, so
//! `i₀ = 0` throughout and `c` collapses to `f₂f₄ − f₀f₆`. The general form is
//! written out above because it is what the derivation produces; the specialised
//! form is what runs.
//!
//! # What is *not* here
//!
//! Whether six saddles mean a tunnel or a single twelve-vertex contour. That
//! distinction **cannot be made from the saddles alone** and does not need to be:
//! it follows from how many closed contours the cell's cut edges form, which is
//! A-002f's. The paper states an asymptote-side criterion for it (Proposition 1,
//! Corollary 1) whose precise predicate its prose does not pin down; the authors'
//! implementation does not evaluate one, branching on the contour count instead.
//! Following the program rather than inventing a predicate is rule 5 working as
//! intended — see V-31.

#[cfg(test)]
mod tests;

use crate::cube::EDGE_COUNT;
use crate::real::Real;

use super::table::{NO_EDGE, segment_links};

/// The number of candidate body saddles a cell can have. Two, from one quadratic.
pub const SADDLE_COUNT: usize = 2;

/// Bit `k` of [`BodySaddles::inside_mask`], by coordinate and solution.
///
/// The layout matches the authors' implementation so that the two can be compared
/// directly: `u₀, u₁, v₀, v₁, w₀, w₁` in bits `0..6`.
const fn coordinate_bit(axis: usize, solution: usize) -> u8 {
    1 << (2 * axis + solution)
}

/// All six coordinate bits set — the configuration that has an inner hexagon.
pub const ALL_INSIDE: u8 = 0b0011_1111;

/// The trilinear interpolant's body saddles within one cell.
///
/// Constructed by [`BodySaddles::of`]. The two candidates share one `u` set, one
/// `v` set and one `w` set; solution `k` is the point `(u[k], v[k], w[k])`.
///
/// A coordinate is **usable** only when [`BodySaddles::inside_mask`] marks it.
/// Coordinates that are not marked may be non-finite: the linear solve for `v` or
/// `w` divides by a difference that vanishes when the interpolant does not vary
/// along that axis, and no sentinel is substituted for the result. The mask is the
/// single place that says which numbers mean anything, which is why every accessor
/// that returns a *point* is gated on it.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BodySaddles<R> {
    coordinate: [[R; SADDLE_COUNT]; 3],
    inside: u8,
}

impl<R: Real> BodySaddles<R> {
    /// Locate the body saddles of the interpolant on this cell.
    ///
    /// `corner` is the eight corner values in `cube.rs`'s numbering, which
    /// is also Grosso's. The zero level set is assumed, as everywhere in this
    /// crate.
    #[must_use]
    pub fn of(corner: &[R; 8]) -> Self {
        let [a, b, c] = Self::coefficients(corner);
        let (u, roots) = Self::roots(a, b, c);

        // `v` is where the `w = 0` face's hyperbola sits at this `u`, and `w` is
        // where the `v = 0` face's hyperbola sits. Both are linear, which is what
        // makes one quadratic enough for all three coordinates.
        let mut v = [R::ZERO; SADDLE_COUNT];
        let mut w = [R::ZERO; SADDLE_COUNT];
        for ((vk, wk), &u) in v.iter_mut().zip(w.iter_mut()).zip(u.iter()).take(roots) {
            *vk = Self::level_crossing(corner[0], corner[1], corner[2], corner[3], u);
            *wk = Self::level_crossing(corner[0], corner[1], corner[4], corner[5], u);
        }
        let coordinate = [u, v, w];

        let mut inside = 0u8;
        for (axis, values) in coordinate.iter().enumerate() {
            for (k, &value) in values.iter().enumerate().take(roots) {
                if value > R::ZERO && value < R::ONE {
                    inside |= coordinate_bit(axis, k);
                }
            }
        }

        Self { coordinate, inside }
    }

    /// The quadratic's coefficients `[a, b, c]`, at the zero level set.
    ///
    /// Separate from [`BodySaddles::of`] so that
    /// `the_coefficients_reproduce_the_face_hyperbola_difference` can check them
    /// against a direct evaluation of the expression they were derived from,
    /// rather than against themselves.
    #[must_use]
    pub fn coefficients(corner: &[R; 8]) -> [R; 3] {
        // The two faces' bilinear "twist" terms, and the two edge differences the
        // derivation pairs them with.
        let twist_lo = (corner[0] + corner[3]) - (corner[1] + corner[2]);
        let twist_hi = (corner[4] + corner[7]) - (corner[5] + corner[6]);
        let du_lo = corner[1] - corner[0];
        let du_hi = corner[5] - corner[4];
        let dv_lo = corner[2] - corner[0];
        let dv_hi = corner[6] - corner[4];

        let a = du_hi * twist_lo - du_lo * twist_hi;
        // `i₀ = 0`, so `(i₀ − f₀)` is `−f₀` and `(i₀ − f₄)` is `−f₄`.
        let b = (corner[4] * twist_lo - corner[0] * twist_hi) + (du_hi * dv_lo - du_lo * dv_hi);
        let c = corner[2] * corner[4] - corner[0] * corner[6];
        [a, b, c]
    }

    /// The real roots of `a·u² + b·u + c`, and how many there are.
    ///
    /// # Why not the textbook formula
    ///
    /// The authors' implementation uses `(−b ± √d)/2a`, which cancels
    /// catastrophically in whichever branch subtracts `√d` from a near-equal `|b|`
    /// — and `a` is a difference of near-equal products with nothing keeping it
    /// away from zero. [`super::interior`] rejected the same formula for the same
    /// reason at A-002c. Kahan's form adds magnitudes instead, so `q` is accurate
    /// to rounding; `q/a` is the large root and, by Vieta, `c/q` the small one.
    ///
    /// # Why `a == 0` is a root count and not an absence
    ///
    /// With `a` zero the equation is linear and has **one** root, `−c/b`. The
    /// textbook formula loses it — it divides by `2a` and yields infinities — and
    /// the authors' implementation inherits that loss. Solving the smaller
    /// polynomial as a smaller polynomial keeps the count honest, which matters
    /// because Grosso §5.3 selects the cell's interior vertex by *how many* face
    /// pairs have a single solution. The divergence is deliberate and measured;
    /// see M-207.
    fn roots(a: R, b: R, c: R) -> ([R; SADDLE_COUNT], usize) {
        let mut roots = [R::ZERO; SADDLE_COUNT];
        if a == R::ZERO {
            if b == R::ZERO {
                return (roots, 0);
            }
            roots[0] = -c / b;
            return (roots, 1);
        }

        let discriminant = b * b - R::TWO * R::TWO * a * c;
        if discriminant < R::ZERO {
            return (roots, 0);
        }
        if discriminant == R::ZERO {
            // A double root is **one** intersection point, not two: the two
            // hyperbolas touch rather than cross, and Proposition 1 counts points.
            // Reporting two here would let a degenerate zero-area "hexagon" claim
            // six saddles.
            roots[0] = -b / (R::TWO * a);
            return (roots, 1);
        }

        // `signum(0)` is `+1` for `+0.0`; either sign keeps this a sum of
        // magnitudes, because `b` is then zero and adds nothing. `q` cannot be
        // zero here: that would need `|b| = −|b|`, hence `b == 0` *and*
        // `discriminant == 0`, which the branch above has already taken.
        let q = -(b + b.signum() * discriminant.sqrt()) * R::HALF;
        roots[0] = q / a;
        roots[1] = c / q;
        (roots, SADDLE_COUNT)
    }

    /// Where the level set crosses the segment from `lo` to `hi` at parameter `u`.
    ///
    /// `lo` interpolates `lo0 → lo1` across `u` and `hi` likewise; the return value
    /// is the coordinate at which the interpolation between them reaches zero. Not
    /// guarded: the caller's mask records whether the result is usable, and an
    /// epsilon here would move a saddle rather than reject it.
    fn level_crossing(lo0: R, lo1: R, hi0: R, hi1: R, u: R) -> R {
        let s = R::ONE - u;
        let lo = lo0 * s + lo1 * u;
        let hi = hi0 * s + hi1 * u;
        -lo / (hi - lo)
    }

    /// Which of the six coordinates lie strictly inside the cell.
    ///
    /// Bits `0..6` are `u₀, u₁, v₀, v₁, w₀, w₁`. Strict on both ends: a coordinate
    /// of exactly `0` or `1` places the saddle *on* a face rather than inside the
    /// cell, which is the singular configuration A-002i owns and not a body saddle.
    #[must_use]
    pub const fn inside_mask(&self) -> u8 {
        self.inside
    }

    /// How many of the six coordinates lie strictly inside the cell.
    #[must_use]
    pub const fn inside_count(&self) -> u32 {
        self.inside.count_ones()
    }

    /// Does this cell have all six saddle coordinates, and therefore an inner
    /// hexagon?
    ///
    /// **This is not the same as "has a tunnel".** Six saddles mean the cell is
    /// either a tunnel or a single twelve-vertex contour; which one it is follows
    /// from the number of closed contours the cut edges form, not from the saddles.
    /// See this module's header.
    #[must_use]
    pub const fn has_inner_hexagon(&self) -> bool {
        self.inside == ALL_INSIDE
    }

    /// The coordinates along one axis, whether or not they are usable.
    ///
    /// `axis` is `0` for `u`, `1` for `v`, `2` for `w`. Consult
    /// [`inside_mask`](Self::inside_mask) before using a value.
    ///
    /// # Panics
    ///
    /// If `axis` is not less than three.
    #[must_use]
    pub fn axis(&self, axis: usize) -> [R; SADDLE_COUNT] {
        self.coordinate[axis]
    }

    /// The six vertices of the inner hexagon, in order around it.
    ///
    /// `None` unless every coordinate is inside the cell, since the hexagon is not
    /// otherwise defined. Consecutive vertices differ in exactly one coordinate,
    /// so every edge of the hexagon is parallel to an axis — which is what
    /// Grosso's Proposition 2 says, and
    /// `hexagon_edges_are_axis_parallel_and_close` checks.
    ///
    /// # The order is the reference implementation's, not the paper's
    ///
    /// The paper's own listing of these six points is corrupt in the copy this
    /// project holds — its first branch assigns `p₁` and `p₂` the same triple,
    /// which cannot be a hexagon. Guessing the intended point is exactly what
    /// rule 5 forbids, so the order below is the authors' program's (V-31).
    #[must_use]
    pub fn inner_hexagon(&self) -> Option<[[R; 3]; 6]> {
        if !self.has_inner_hexagon() {
            return None;
        }
        let [u, v, w] = self.coordinate;
        Some([
            [u[0], v[0], w[0]],
            [u[0], v[0], w[1]],
            [u[1], v[0], w[1]],
            [u[1], v[1], w[1]],
            [u[1], v[1], w[0]],
            [u[0], v[1], w[0]],
        ])
    }
}

/// Most closed contours one cell's cut edges can form.
///
/// Twelve cut edges and three to a ring is four, and four is reached — by the
/// configuration whose four inside corners are the corners of one tetrahedron,
/// each isolated from the others.
pub const MAX_CONTOURS: usize = 4;

/// The closed rings a cell's cut edges form, in order around each ring.
///
/// Grosso's Algorithm 1 step 2: *"the intersection of the isosurface with the
/// cell is computed at the cell edges. From these intersection points a set of up
/// to four contours is generated **using the asymptotic decider to preserve
/// consistency across cells borders**."*
///
/// # Crack-freeness is inherited here rather than argued
///
/// The rings come from [`super::table::segment_links`], which is the *same*
/// function the table and the face decider already use, driven by the *same*
/// `joined` mask [`super::ambiguity::joined_mask`] computes. So a cell's face
/// connectivity is untouched by anything in this module, and ✗11's face-locality
/// property — a face's segments are a function of that face's own four corner
/// signs and nothing else — keeps covering it, as does A-002's
/// [`validate_decider_table`](super::validate_decider_table) over all 16,384
/// `(case, mask)` pairs. **The interior classification cannot move a face
/// segment**, which is the whole reason this construction needs no grid
/// subdivision pass where Custodio's does (V-29).
///
/// # Why the tunnel test lives here and not in [`BodySaddles`]
///
/// Six body saddles mean the cell is *either* a tunnel *or* one twelve-vertex
/// contour, and the saddles cannot tell you which. The number of contours can,
/// and that is what the authors' own implementation branches on. See
/// [`Contours::topology`].
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Contours {
    /// Every ring's edges, laid end to end.
    edges: [u8; EDGE_COUNT],
    /// Where each ring starts in `edges`, plus a final sentinel at the end.
    start: [u8; MAX_CONTOURS + 1],
    count: u8,
}

impl Contours {
    /// Walk a cell's cut edges into closed rings.
    ///
    /// `case` is the corner-sign index and `joined` the per-face resolution mask,
    /// exactly as [`super::table::segment_links`] takes them.
    #[must_use]
    pub fn of(case: u8, joined: u8) -> Self {
        let next = segment_links(case, joined);
        let mut edges = [NO_EDGE; EDGE_COUNT];
        let mut start = [0u8; MAX_CONTOURS + 1];
        let mut count = 0usize;
        let mut filled = 0usize;
        let mut visited = 0u16;

        for first in 0..EDGE_COUNT as u8 {
            if next[first as usize] == NO_EDGE || visited & (1 << first) != 0 {
                continue;
            }
            debug_assert!(count < MAX_CONTOURS, "more rings than a cell can hold");
            start[count] = filled as u8;
            count += 1;

            let mut current = first;
            while visited & (1 << current) == 0 {
                visited |= 1 << current;
                edges[filled] = current;
                filled += 1;
                current = next[current as usize];
            }
        }
        start[count] = filled as u8;

        Self {
            edges,
            start,
            count: count as u8,
        }
    }

    /// How many closed rings this cell has. Zero for an empty or full cell.
    #[must_use]
    pub const fn count(&self) -> usize {
        self.count as usize
    }

    /// One ring's cut edges, in order around it.
    ///
    /// # Panics
    ///
    /// If `index` is not less than [`count`](Self::count).
    #[must_use]
    pub fn ring(&self, index: usize) -> &[u8] {
        assert!(index < self.count(), "contour {index} does not exist");
        let (lo, hi) = (self.start[index] as usize, self.start[index + 1] as usize);
        &self.edges[lo..hi]
    }

    /// The length of the longest ring, or zero if there are none.
    #[must_use]
    pub fn longest(&self) -> usize {
        (0..self.count())
            .map(|i| self.ring(i).len())
            .max()
            .unwrap_or(0)
    }

    /// What the cell's surface is, once the saddles and the rings are both known.
    ///
    /// **This is the decision six body saddles cannot make alone**, and the one
    /// the paper states through an asymptote-side criterion (Proposition 1,
    /// Corollary 1) whose precise predicate its prose never pins down. The
    /// authors' implementation does not evaluate one — it branches on how many
    /// rings the cell has, and so does this (V-31):
    ///
    /// - fewer than six saddles → every ring is a disk;
    /// - six saddles and **one** ring → that ring has twelve vertices and is
    ///   still a disk, just a long one;
    /// - six saddles and **two or three** rings → two of them are the two ends of
    ///   a tunnel, and a third, if present, is a separate three-vertex disk.
    #[must_use]
    pub fn topology<R: Real>(&self, saddles: &BodySaddles<R>) -> Topology {
        if !saddles.has_inner_hexagon() {
            return Topology::Disks;
        }
        match self.count() {
            0 | 1 => Topology::TwelveVertexContour,
            _ => Topology::Tunnel,
        }
    }
}

/// What the trilinear interpolant does inside one cell.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum Topology {
    /// Every contour bounds a disk. The overwhelming majority of cells.
    Disks,
    /// Two contours are the ends of a cylinder through the cell.
    Tunnel,
    /// A single contour through all twelve cut edges, still a disk. Only
    /// reachable from Marching Cubes' case 13, and rare even there — Grosso
    /// counts **7** in a 512²×641 CT skull, against 2,057 tunnels.
    TwelveVertexContour,
}

/// One of the three pairs of opposite faces, named by the axis its lines run
/// along.
///
/// Pair `U` is the faces `u = 0` and `u = 1`, joined by lines parallel to `u`,
/// and so on.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum Pair {
    /// Faces `u = 0` and `u = 1`.
    U = 0,
    /// Faces `v = 0` and `v = 1`.
    V = 1,
    /// Faces `w = 0` and `w = 1`.
    W = 2,
}

impl Pair {
    /// All three, in `u, v, w` order.
    pub const ALL: [Pair; 3] = [Pair::U, Pair::V, Pair::W];
}

impl<R: Real> BodySaddles<R> {
    /// Which of a pair's two lines lie inside the cell, as a two-bit mask.
    ///
    /// A line joining one pair of opposite faces is fixed by the *other two*
    /// coordinates, and it lies inside the cell exactly when both of them do.
    ///
    /// # The `u` pair's two solutions are crossed, and this is derived not copied
    ///
    /// The obvious reading is that solution `k` contributes the line at
    /// `(v[k], w[k])`, and for the `v` and `w` pairs that is right. **For the `u`
    /// pair it is wrong**, and the inner hexagon is what shows it. Consecutive
    /// hexagon vertices differ in exactly one coordinate, so a pair of vertices
    /// differing *only* in `u` is a `u`-parallel line; reading them off
    /// [`inner_hexagon`](Self::inner_hexagon)'s ring gives
    ///
    /// ```text
    /// (u₀,v₀,w₁) ↔ (u₁,v₀,w₁)   →   the u-line at (v₀, w₁)
    /// (u₀,v₁,w₀) ↔ (u₁,v₁,w₀)   →   the u-line at (v₁, w₀)
    /// ```
    ///
    /// — indices `0,1` and `1,0`, **crossed**. The same reading gives `(u₀,v₀)`
    /// and `(u₁,v₁)` for the `w` pair and `(u₀,w₀)` and `(u₁,w₁)` for the `v`
    /// pair, both uncrossed. The authors' implementation agrees, and it was the
    /// disagreement with the obvious reading that made this worth deriving:
    /// `the_line_inventory_agrees_with_the_inner_hexagon` checks all three pairs
    /// against the ring rather than against the transcription (M-215).
    #[must_use]
    pub fn lines(&self, pair: Pair) -> u8 {
        let inside = self.inside;
        // `(axis, solution)` index pairs, per line, per pair of opposite faces.
        let needed: [[(usize, usize); 2]; 2] = match pair {
            // A u-line is fixed by (v, w) — crossed.
            Pair::U => [[(1, 1), (2, 0)], [(1, 0), (2, 1)]],
            // A v-line is fixed by (u, w).
            Pair::V => [[(0, 0), (2, 0)], [(0, 1), (2, 1)]],
            // A w-line is fixed by (u, v).
            Pair::W => [[(0, 0), (1, 0)], [(0, 1), (1, 1)]],
        };
        let mut mask = 0u8;
        for (line, coordinates) in needed.iter().enumerate() {
            if coordinates
                .iter()
                .all(|&(axis, k)| inside & coordinate_bit(axis, k) != 0)
            {
                mask |= 1 << line;
            }
        }
        mask
    }

    /// How many of a pair's two lines lie inside the cell.
    #[must_use]
    pub fn line_count(&self, pair: Pair) -> u32 {
        self.lines(pair).count_ones()
    }

    /// The point where a pair's line `line` crosses the cell, given the value of
    /// the coordinate it runs along.
    ///
    /// The line is fixed in the other two coordinates; `along` supplies the third.
    ///
    /// Test-only: nothing in the extractor needs a point *on* a line, only the
    /// saddles the lines cross at. It exists so
    /// `the_line_inventory_agrees_with_the_inner_hexagon` can check the pairing
    /// against the hexagon's ring rather than against the transcription it agrees
    /// with. A-002h may promote it when the hexagon is meshed.
    #[cfg(test)]
    pub(super) fn point_on_line_for_test(&self, pair: Pair, line: usize, along: R) -> [R; 3] {
        let [u, v, w] = self.coordinate;
        match pair {
            Pair::U if line == 0 => [along, v[1], w[0]],
            Pair::U => [along, v[0], w[1]],
            Pair::V if line == 0 => [u[0], along, w[0]],
            Pair::V => [u[1], along, w[1]],
            Pair::W if line == 0 => [u[0], v[0], along],
            Pair::W => [u[1], v[1], along],
        }
    }
}

impl<R: Real> BodySaddles<R> {
    /// The one interior vertex an ambiguous disk contour is fanned from, in the
    /// cell's own `[0,1]³` coordinates.
    ///
    /// Grosso §5.3. A contour of six or fewer vertices that touches no ambiguous
    /// face is a flat polygon and needs nothing; a contour of six to nine that
    /// crosses an ambiguous face **twice** cannot be triangulated manifoldly
    /// without one, *"otherwise a third edge in one of the ambiguous faces has to
    /// be introduced, violating manifoldness across cells borders."*
    ///
    /// `None` when every in-range line belongs to a single pair of opposite faces,
    /// because then no two lines cross and there is no interior saddle to use —
    /// and none is needed, since such a contour does not double back.
    ///
    /// # This branch structure is the reference implementation's, and the check
    /// on it is geometric rather than textual
    ///
    /// The paper says only *"if there are only two such lines, the additional
    /// vertex is the intersection point. If there are three such lines, the
    /// additional inner vertex is the midpoint between the two saddle points"* —
    /// which does not determine the four-line case, and does not say which of
    /// three lines' three pairwise intersections the "two saddle points" are. The
    /// index selection below is therefore transcribed from the authors' program
    /// (V-31), and **verified by the one property that a transcription error
    /// cannot survive**: the point it returns must lie *on* the level set.
    /// `the_interior_vertex_lies_on_the_level_set` checks exactly that against an
    /// independent evaluation of the interpolant, which is the same instrument
    /// that validated the inner hexagon.
    #[must_use]
    pub fn interior_vertex(&self) -> Option<[R; 3]> {
        let set = |axis: usize, k: usize| self.inside & coordinate_bit(axis, k) != 0;
        // A 0/1 weight, so every "select the in-range one" below is an exact
        // multiply-and-add rather than a branch.
        let f = |axis: usize, k: usize| {
            if set(axis, k) { R::ONE } else { R::ZERO }
        };
        let count = |a: (usize, usize), b: (usize, usize)| u8::from(set(a.0, a.1) && set(b.0, b.1));

        // Lines per pair of opposite faces. The `u` pair is crossed; see
        // [`lines`](Self::lines) for why that is derived and not copied.
        let fc_w = count((0, 0), (1, 0)) + count((0, 1), (1, 1));
        let fc_v = count((0, 0), (2, 0)) + count((0, 1), (2, 1));
        let fc_u = count((1, 0), (2, 1)) + count((1, 1), (2, 0));
        let total = fc_w + fc_v + fc_u;

        // Every line in one pair — parallel lines never cross, so there is no
        // saddle. Also covers `total == 0`.
        if total == fc_w || total == fc_v || total == fc_u {
            return None;
        }

        let [u, v, w] = self.coordinate;
        let (uc, vc, wc) = match total {
            2 if fc_w == 0 => (
                f(0, 0) * u[0] + f(0, 1) * u[1],
                f(1, 0) * v[0] + f(1, 1) * v[1],
                f(1, 0) * w[1] + f(1, 1) * w[0],
            ),
            2 if fc_v == 0 => (
                f(0, 0) * u[0] + f(0, 1) * u[1],
                f(0, 0) * v[0] + f(0, 1) * v[1],
                f(0, 0) * w[1] + f(0, 1) * w[0],
            ),
            2 => (
                f(1, 0) * u[0] + f(1, 1) * u[1],
                f(1, 0) * v[0] + f(1, 1) * v[1],
                f(1, 0) * w[0] + f(1, 1) * w[1],
            ),
            // Three lines: the mean of the two solutions on each axis, which is
            // the paper's "midpoint between the two saddle points" written so that
            // an axis with only one solution contributes that one unchanged.
            3 => (
                (f(0, 0) * u[0] + f(0, 1) * u[1]) / (f(0, 0) + f(0, 1)),
                (f(1, 0) * v[0] + f(1, 1) * v[1]) / (f(1, 0) + f(1, 1)),
                (f(2, 0) * w[0] + f(2, 1) * w[1]) / (f(2, 0) + f(2, 1)),
            ),
            4 => {
                let on = |axis: usize| f(axis, 0) + f(axis, 1);
                if on(2) == R::ONE {
                    (
                        f(2, 0) * u[0] + f(2, 1) * u[1],
                        f(2, 1) * v[0] + f(2, 0) * v[1],
                        f(2, 0) * w[0] + f(2, 1) * w[1],
                    )
                } else if on(1) == R::ONE {
                    (
                        f(1, 0) * u[0] + f(1, 1) * u[1],
                        f(1, 0) * v[0] + f(1, 1) * v[1],
                        f(1, 1) * w[0] + f(1, 0) * w[1],
                    )
                } else {
                    (
                        f(0, 0) * u[0] + f(0, 1) * u[1],
                        f(0, 0) * v[0] + f(0, 1) * v[1],
                        f(0, 0) * w[0] + f(0, 1) * w[1],
                    )
                }
            }
            // Five or six lines is the hexagon's territory, and a cell there is a
            // tunnel or a twelve-vertex contour rather than a fanned disk. A-002h
            // owns it; `the_disk_path_never_sees_five_or_six_lines` asserts this
            // arm is unreachable from the disk path rather than assuming it.
            _ => return None,
        };
        Some([uc, vc, wc])
    }
}

/// The code a triangle uses to name the cell's one interior vertex, as opposed
/// to one of the twelve cut edges.
///
/// Deliberately the same value as [`super::table::CENTROID_BASE`]: both name a
/// **cell-local** vertex that no other cell can reach, which is exactly why
/// neither is cached and why A-015's index-space budget already accounts for the
/// slot. What differs is where the position comes from — a centroid averages its
/// cycle's edge vertices, this one is a saddle of the interpolant
/// ([`BodySaddles::interior_vertex`]).
pub const INTERIOR: u8 = EDGE_COUNT as u8;

impl Contours {
    /// Fan every ring into triangles, naming cut edges by index and the cell's
    /// interior vertex by [`INTERIOR`].
    ///
    /// `interior` says whether [`BodySaddles::interior_vertex`] produced one. When
    /// it did, rings of more than three vertices fan from it — `k` triangles for a
    /// ring of `k` — because a ring that crosses an ambiguous face twice cannot be
    /// closed manifoldly from one of its own vertices without laying a third edge
    /// in that face. A ring of exactly three is a triangle either way and takes
    /// the direct emission, as the reference implementation does.
    ///
    /// When there is none, every ring fans from its own first vertex: `k − 2`
    /// triangles, and no interior vertex is created at all.
    ///
    /// Winding follows the ring order, which [`super::table::segment_links`]
    /// already orients — so this inherits A-001's counter-clockwise-from-outside
    /// convention rather than choosing one.
    pub fn fan(&self, interior: bool, mut emit: impl FnMut([u8; 3])) {
        for r in 0..self.count() {
            let ring = self.ring(r);
            if ring.len() == 3 {
                emit([ring[0], ring[1], ring[2]]);
            } else if interior {
                for k in 0..ring.len() {
                    emit([ring[k], ring[(k + 1) % ring.len()], INTERIOR]);
                }
            } else {
                for k in 1..ring.len() - 1 {
                    emit([ring[0], ring[k], ring[k + 1]]);
                }
            }
        }
    }

    /// How many triangles [`fan`](Self::fan) will emit.
    #[must_use]
    pub fn triangle_count(&self, interior: bool) -> usize {
        (0..self.count())
            .map(|r| {
                let k = self.ring(r).len();
                if k == 3 {
                    1
                } else if interior {
                    k
                } else {
                    k - 2
                }
            })
            .sum()
    }
}
