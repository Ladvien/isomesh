//! Normal coordinates: how many times the surface crosses each edge.
//!
//! From Baktash, Gillespie & Crane, *Subgrid Marching Tetrahedra*,
//! `10.48550/arXiv.2606.00454` §2, read this session. Definitions and the two
//! validity conditions are quoted where they are used; nothing here is invented.
//!
//! # The encoding
//!
//! A tetrahedron has six edges, and [`EdgeCoordinates`] is one non-negative
//! integer per edge — *"edge coordinates `eᵢⱼ` count the number of intersections
//! of the surface `S` with each of the six edges `ij`."* Classic Marching
//! Tetrahedra is the case where every one of those six numbers is 0 or 1, which
//! is why it *"reinvented a small piece of this story"*.
//!
//! # Edge ordering, and the one place it matters
//!
//! The paper writes `e := (e₀₁, e₀₂, e₀₃, e₂₃, e₁₃, e₁₂)`, ordered so that
//! **complementary pairs are three apart**. This module uses
//! [`crate::marching_tetrahedra::table::TET_EDGES`] instead — the
//! crate's own lexicographic order, already load-bearing in A-003 — and recovers
//! the pairing as [`complementary`], which is `5 - edge`.
//!
//! That is not a cosmetic difference. The quad basis vectors are *defined* by
//! which pair of opposite edges they separate, so a mismatched ordering would
//! decompose configurations into the wrong polygons silently.
//! `complementary_edges_share_no_corner` is the check that the two conventions
//! agree.

use crate::marching_tetrahedra::table::{TET_EDGE_COUNT, TET_EDGES};

/// A tetrahedron has four triangular faces.
pub const TET_FACE_COUNT: usize = 4;

/// The three corners of each face, and the edge joining each consecutive pair.
///
/// Face `f` is the one **opposite corner `f`**, so it is the face corner `f` is
/// *not* on. `edge[k]` joins `corner[k]` to `corner[(k + 1) % 3]`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct TetFace {
    /// The face's three corners.
    pub corner: [u8; 3],
    /// The edge joining each consecutive corner pair.
    pub edge: [u8; 3],
}

/// The four faces, built from [`TET_EDGES`] rather than written out.
pub static TET_FACES: [TetFace; TET_FACE_COUNT] = build_faces();

const fn build_faces() -> [TetFace; TET_FACE_COUNT] {
    let mut out = [TetFace {
        corner: [0; 3],
        edge: [0; 3],
    }; TET_FACE_COUNT];

    let mut opposite = 0usize;
    while opposite < TET_FACE_COUNT {
        // The three corners that are not `opposite`, in increasing order.
        let mut corner = [0u8; 3];
        let mut n = 0usize;
        let mut c = 0u8;
        while c < 4 {
            if c as usize != opposite {
                corner[n] = c;
                n += 1;
            }
            c += 1;
        }

        let mut edge = [0u8; 3];
        let mut k = 0usize;
        while k < 3 {
            edge[k] = edge_between(corner[k], corner[(k + 1) % 3]);
            k += 1;
        }

        out[opposite] = TetFace { corner, edge };
        opposite += 1;
    }
    out
}

/// The index of the edge joining two tet corners.
///
/// # Panics
///
/// At compile time, if the corners are equal or out of range — a tetrahedron's
/// corners are all pairwise adjacent, so every distinct pair names an edge.
pub const fn edge_between(a: u8, b: u8) -> u8 {
    let (lo, hi) = if a < b { (a, b) } else { (b, a) };
    let mut e = 0usize;
    while e < TET_EDGE_COUNT {
        if TET_EDGES[e][0] == lo && TET_EDGES[e][1] == hi {
            return e as u8;
        }
        e += 1;
    }
    panic!("a tetrahedron's corners are pairwise adjacent")
}

/// The edge sharing no corner with this one.
///
/// A tetrahedron's six edges fall into three **complementary pairs**, and under
/// [`TET_EDGES`]' lexicographic order the pairing is exactly `5 - edge`:
/// `01↔23`, `02↔13`, `03↔12`. Asserted, not assumed — see
/// `complementary_edges_share_no_corner`.
///
/// The pairing is what a quad is defined by: *"`qᵢⱼ` for the number of quads
/// separating edge `ij` from the complementary edge."*
#[inline]
#[must_use]
pub const fn complementary(edge: u8) -> u8 {
    (TET_EDGE_COUNT as u8 - 1) - edge
}

/// How many times the surface crosses each of a tetrahedron's six edges.
///
/// In [`TET_EDGES`] order. Classic Marching Tetrahedra is the sub-case where
/// every entry is 0 or 1.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct EdgeCoordinates {
    /// Intersection count per edge, in [`TET_EDGES`] order.
    pub count: [u32; TET_EDGE_COUNT],
}

/// Why a set of edge coordinates does not describe a normal curve system.
///
/// Both conditions are §2.1's, with the paper's own one-line justifications
/// attached — they are the clearest statement of *why* these two and not others.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum NotNormal {
    /// *"The even sum condition captures the idea that **what goes in, must come
    /// out**."*
    ///
    /// A face whose three counts sum to an odd number has an arc with nowhere to
    /// go.
    OddSum {
        /// Which face.
        face: u8,
        /// What its three counts sum to.
        sum: u32,
    },
    /// *"The triangle inequality captures the idea that arcs entering one edge
    /// must **exit a different edge**."*
    ///
    /// `eᵢⱼ + eₖᵢ ≥ eⱼₖ` must hold at every corner of every face.
    TriangleInequality {
        /// Which face.
        face: u8,
        /// The corner of that face where it fails.
        corner: u8,
    },
}

impl EdgeCoordinates {
    /// A tetrahedron the surface misses entirely.
    #[must_use]
    pub const fn empty() -> Self {
        Self {
            count: [0; TET_EDGE_COUNT],
        }
    }

    /// Build from a per-edge count.
    #[must_use]
    pub const fn new(count: [u32; TET_EDGE_COUNT]) -> Self {
        Self { count }
    }

    /// The crossings on one edge.
    #[must_use]
    pub const fn edge(&self, edge: u8) -> u32 {
        self.count[edge as usize]
    }

    /// Total crossings over all six edges.
    ///
    /// The paper's `Σe = n`, which bounds its reconstruction: *"`O(n)`
    /// subdivisions"*, and *"subdivision … requires ≥24 edge intersections"*.
    #[must_use]
    pub fn total(&self) -> u32 {
        self.count.iter().sum()
    }

    /// `true` when every edge carries at most one crossing.
    ///
    /// The classic Marching Tetrahedra regime, and the boundary of what a
    /// sign-based method can represent at all.
    #[must_use]
    pub fn is_classic(&self) -> bool {
        self.count.iter().all(|c| *c <= 1)
    }

    /// The corner coordinates of one face, or `None` if this face's counts do not
    /// describe a normal curve.
    ///
    /// §2.1, verbatim:
    ///
    /// ```text
    /// c₀ = (e₀₁ + e₂₀ − e₁₂) / 2,  c₁ = (e₁₂ + e₀₁ − e₂₀) / 2,  c₂ = (e₂₀ + e₁₂ − e₀₁) / 2
    /// ```
    ///
    /// `cᵢ` is *"the number of segments separating each vertex `i` from the other
    /// two"* — the arcs that cut off that corner. Returned in the face's own
    /// corner order, matching [`TetFace::corner`].
    #[must_use]
    pub fn corner_coordinates(&self, face: u8) -> Option<[u32; 3]> {
        let f = TET_FACES[face as usize];
        // `edge[k]` joins corner k to corner k+1, so the edge *opposite* corner k
        // is `edge[k + 1]`.
        let e = [
            i64::from(self.count[f.edge[0] as usize]),
            i64::from(self.count[f.edge[1] as usize]),
            i64::from(self.count[f.edge[2] as usize]),
        ];
        if (e[0] + e[1] + e[2]) % 2 != 0 {
            return None;
        }
        let mut out = [0u32; 3];
        for (k, slot) in out.iter_mut().enumerate() {
            // Corner k is on edges k and k-1; the opposite edge is k+1.
            let adjacent_a = e[k];
            let adjacent_b = e[(k + 2) % 3];
            let opposite = e[(k + 1) % 3];
            let c = (adjacent_a + adjacent_b - opposite) / 2;
            if c < 0 {
                return None;
            }
            *slot = c as u32;
        }
        Some(out)
    }

    /// Check both of §2.1's conditions on every face.
    ///
    /// # Errors
    ///
    /// [`NotNormal`] naming the first face and corner that fails. Reported rather
    /// than clamped: a configuration that violates these is not a curve system at
    /// all, and inventing the nearest valid one would silently mesh a different
    /// surface.
    pub fn validate(&self) -> Result<(), NotNormal> {
        for face in 0..TET_FACE_COUNT as u8 {
            let f = TET_FACES[face as usize];
            let e: [u32; 3] = [
                self.count[f.edge[0] as usize],
                self.count[f.edge[1] as usize],
                self.count[f.edge[2] as usize],
            ];
            let sum = e[0] + e[1] + e[2];
            if !sum.is_multiple_of(2) {
                return Err(NotNormal::OddSum { face, sum });
            }
            for k in 0..3 {
                // eᵢⱼ + eₖᵢ ≥ eⱼₖ at corner k.
                if e[k] + e[(k + 2) % 3] < e[(k + 1) % 3] {
                    return Err(NotNormal::TriangleInequality {
                        face,
                        corner: f.corner[k],
                    });
                }
            }
        }
        Ok(())
    }

    /// `true` when both §2.1 conditions hold on all four faces.
    #[must_use]
    pub fn is_normal(&self) -> bool {
        self.validate().is_ok()
    }
}

/// Normal surface coordinates: `n := (t₀, t₁, t₂, t₃, q₀₁, q₀₂, q₀₃)`.
///
/// §2.3: *"`tᵢ` to denote the number of triangles at corner `i`, and `qᵢⱼ` for the
/// number of quads separating edge `ij` from the complementary edge."*
///
/// The quads are indexed by the **lower** edge of each complementary pair, so
/// `quad[0]` separates edge 0 (`01`) from edge 5 (`23`), and so on.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct NormalSurface {
    /// Triangles cutting off each corner.
    pub triangle: [u32; 4],
    /// Quads separating each complementary edge pair.
    ///
    /// **At most one may be non-zero.** §2.3: *"for this set of polygons to be
    /// intersection-free, there can be at most one type of diagonal cut, i.e.,
    /// only one of the coordinates `qᵢⱼ` can be nonzero."*
    pub quad: [u32; 3],
}

impl NormalSurface {
    /// The edge coordinates this collection of polygons produces.
    ///
    /// The forward map `e = M·n`. A triangle at corner `i` crosses the three
    /// edges at `i`; a quad on pair `(ab | cd)` crosses the four edges that are
    /// neither `ab` nor `cd`.
    #[must_use]
    pub fn edge_coordinates(&self) -> EdgeCoordinates {
        let mut count = [0u32; TET_EDGE_COUNT];
        for (edge, slot) in count.iter_mut().enumerate() {
            let [a, b] = TET_EDGES[edge];
            *slot += self.triangle[a as usize] + self.triangle[b as usize];
        }
        for (pair, q) in self.quad.iter().enumerate() {
            if *q == 0 {
                continue;
            }
            let lo = pair as u8;
            let hi = complementary(lo);
            for (edge, slot) in count.iter_mut().enumerate() {
                if edge as u8 != lo && edge as u8 != hi {
                    *slot += q;
                }
            }
        }
        EdgeCoordinates { count }
    }

    /// How many polygons this is, in total.
    #[must_use]
    pub fn polygon_count(&self) -> u32 {
        self.triangle.iter().sum::<u32>() + self.quad.iter().sum::<u32>()
    }
}

/// Decompose edge coordinates into normal triangles and quads, if they are the
/// image of any.
///
/// Returns `None` when no such decomposition exists — which is **not an error and
/// not rare**. §2.3 is explicit that the map does not cover every input:
///
/// > In contrast, we cannot always explain a given set of edge coordinates via
/// > intersections with normal triangles and quads: solutions … may yield
/// > negative or fractional coordinates, or describe intersecting polygons.
/// > Consider for example the edge coordinates `e = (2,1,1,2,1,1)`… decomposing
/// > `e` into **intersecting** quads.
///
/// A `None` is therefore the signal that a configuration needs A-014b's general
/// reconstruction rather than a lookup, and counting how often it happens is a
/// measurement worth having before building that.
///
/// Solved by trying each of the four possibilities — no quad, or a quad on one of
/// the three complementary pairs — and **verifying by reconstruction**: a
/// candidate is accepted only if [`NormalSurface::edge_coordinates`] reproduces
/// the input exactly. The algebra is not trusted; its answer is checked.
#[must_use]
pub fn decompose(e: &EdgeCoordinates) -> Option<NormalSurface> {
    for pair in [usize::MAX, 0, 1, 2] {
        let mut quad = [0u32; 3];
        let mut residual = [0i64; TET_EDGE_COUNT];
        for (edge, slot) in residual.iter_mut().enumerate() {
            *slot = i64::from(e.count[edge]);
        }

        if pair != usize::MAX {
            let lo = pair as u8;
            let hi = complementary(lo);
            // Summing the four crossed edges counts every triangle twice and the
            // quad four times; the two uncrossed edges count every triangle twice
            // and the quad not at all. So the difference is 4q.
            let mut crossed = 0i64;
            let mut uncrossed = 0i64;
            for edge in 0..TET_EDGE_COUNT as u8 {
                let v = i64::from(e.count[edge as usize]);
                if edge == lo || edge == hi {
                    uncrossed += v;
                } else {
                    crossed += v;
                }
            }
            let numerator = crossed - 2 * uncrossed;
            if numerator <= 0 || numerator % 4 != 0 {
                continue;
            }
            let q = numerator / 4;
            if q > i64::from(u32::MAX) {
                continue;
            }
            quad[pair] = q as u32;
            for (edge, slot) in residual.iter_mut().enumerate() {
                if edge as u8 != lo && edge as u8 != hi {
                    *slot -= q;
                }
            }
        }

        // With the quad removed, every edge count is `t_a + t_b`, so a corner
        // coordinate on any face gives that corner's triangle count directly.
        let Some(triangle) = triangles_from(&residual) else {
            continue;
        };
        let candidate = NormalSurface { triangle, quad };
        if candidate.edge_coordinates() == *e {
            return Some(candidate);
        }
    }
    None
}

/// Solve `r_ij = t_i + t_j` for the four triangle counts.
fn triangles_from(residual: &[i64; TET_EDGE_COUNT]) -> Option<[u32; 4]> {
    // On face 3 (corners 0, 1, 2) the corner coordinate at corner k is t_k.
    let e01 = residual[edge_between(0, 1) as usize];
    let e02 = residual[edge_between(0, 2) as usize];
    let e12 = residual[edge_between(1, 2) as usize];
    let e03 = residual[edge_between(0, 3) as usize];

    let doubled = [e01 + e02 - e12, e01 + e12 - e02, e02 + e12 - e01];
    let mut triangle = [0u32; 4];
    for (k, slot) in triangle.iter_mut().enumerate().take(3) {
        if doubled[k] < 0 || doubled[k] % 2 != 0 {
            return None;
        }
        let t = doubled[k] / 2;
        if t > i64::from(u32::MAX) {
            return None;
        }
        *slot = t as u32;
    }
    let t3 = e03 - i64::from(triangle[0]);
    if t3 < 0 || t3 > i64::from(u32::MAX) {
        return None;
    }
    triangle[3] = t3 as u32;
    Some(triangle)
}

#[cfg(test)]
mod tests;
