//! Where a transition cell's vertices go, and the one identity the seam rests on.
//!
//! [`table`](super::table) says *which* edges are cut and how they link. This
//! says *where* the crossings are — and the whole of A-011b depends on one
//! property of that placement:
//!
//! > A crossing on a **half-resolution** edge must land exactly where the coarse
//! > neighbour's own Marching Cubes pass puts its vertex.
//!
//! Not "within a tolerance". Exactly. Two vertices a float apart leave a crack
//! the width of a rounding error, and A-013 measured what that costs: a seam of
//! unshared vertices that a renderer draws correctly and a collider reads as a
//! hole (M-69).
//!
//! An earlier version of this paragraph added *"the weld can close a seam it can
//! see; it cannot invent a shared vertex where the two sides disagree in the last
//! bit"*, and **that is ✗18** — measured false at R-004. The welder's rule is
//! first fit within `epsilon_for(h) = h · 1e-4`, not bit-identity, and a
//! disagreement of `1.4e-15` is nine orders of magnitude inside it: the offset
//! arithmetic leaves **0** seam-plane boundary edges once welded, against 63–348
//! under a bit-identity merge (M-278). So the exactness below buys **sharing by
//! construction** — which is what M-69's unwelded consumer gets, and what the
//! weld's own order-dependence (R-002) and ability to *create* a non-manifold
//! edge (M-226) make worth having rather than leaning on.
//!
//! Three facts make the identity hold, and each was measured rather than assumed:
//!
//! - **The endpoints coincide.** A half-resolution edge joins two *corner*
//!   samples of the transition face, and M-70 measured that a level-`k` sample
//!   position is bit-identical to the level-0 position of the sample it sits on,
//!   at every spacing tried including `4/35`.
//! - **The values coincide**, because they are the same field evaluated at the
//!   same point.
//! - **The interpolation coincides**, because [`TransitionCell::crossing`] is
//!   `lo + (hi − lo)·t` with `t = a/(a − b)` — character for character what
//!   `marching_cubes`' `edge_position` computes, and the reason it is written out
//!   here rather than shared is that the two take their operands from different
//!   places. `the_half_resolution_crossings_are_the_coarse_neighbours_vertices`
//!   is what keeps them equal.

use crate::cube::{edge_offset, is_inside, place};
use crate::marching_cubes::table::NO_EDGE;
use crate::vec3;
use crate::{MeshSink, Real, Sdf};

use super::table::{
    AMBIGUOUS_FACES, EDGE_COUNT, EDGE_SAMPLES, SAMPLE_COUNT, is_half_resolution, transition_links,
};

/// One transition cell's nine sample values and their world positions.
///
/// The four half-resolution corners are *not* stored twice: they are samples
/// 0, 2, 6 and 8 of the same nine, which is exactly why the seam closes. See the
/// module docs.
#[derive(Clone, Copy, Debug)]
pub struct TransitionCell<R: Real> {
    /// Field value at each of the nine samples.
    pub value: [R; SAMPLE_COUNT],
    /// World position of each of the nine samples, on the **full-resolution**
    /// face.
    pub position: [[R; 3]; SAMPLE_COUNT],
    /// The **half-resolution** face's four corners, `width` along the face normal
    /// from samples 0, 2, 6 and 8 in that order.
    ///
    /// Their *values* are not stored, because they are those four samples'
    /// values — §4.3's *"the four corner values labeled A, B, C, and D … are
    /// duplicated on the opposite face of the cell."* Only the positions differ,
    /// and only by the width.
    pub back: [[R; 3]; 4],
}

/// Which of [`TransitionCell::back`]'s four slots a sample occupies, if any.
///
/// Samples 0, 2, 6 and 8 are the face's corners and the only ones the
/// half-resolution face duplicates.
#[inline]
#[must_use]
pub const fn back_slot(sample: u8) -> Option<usize> {
    match sample {
        0 => Some(0),
        2 => Some(1),
        6 => Some(2),
        8 => Some(3),
        _ => None,
    }
}

impl<R: Real> TransitionCell<R> {
    /// Sample a transition face.
    ///
    /// `origin` is the **grid's** origin — the world position of global sample
    /// `[0, 0, 0]` — and `base` is the global **fine** sample index of this
    /// face's sample 0. `step` is the fine spacing, so the face spans `2·step`
    /// on each in-plane axis, and `u`/`v` are those axes. The nine samples are
    /// laid out `u`-fastest, matching [`table`](super::table)'s numbering.
    ///
    /// `width` displaces the half-resolution face along the remaining axis, and
    /// **it is signed** — positive runs toward increasing coordinate, so the
    /// caller states which side the coarse block is on rather than the cell
    /// guessing. Lengyel's implementation uses `w(k) = 2^(k−2)`, half the
    /// adjacent full-resolution cell.
    ///
    /// **A zero width is legal and is not a shortcut.** §4.3 says it still
    /// *"seamlessly stitch\[es\] multiresolution meshes together"* — and M-74
    /// measured what it costs: every crossing then lies in the face plane, every
    /// triangle is coplanar with it, and the patch stands **exactly**
    /// perpendicular to the surface it is stitching, with no usable normal. That
    /// is the paper's *"severe shading problems"*, and it is why this parameter
    /// exists rather than being deferred.
    ///
    /// # Why an index and not a face origin
    ///
    /// The first version of this took the face's own world origin and added
    /// local offsets to it. That is wrong, and the test caught it: at
    /// `h = 4/14` a half-resolution crossing came out at `y = -1.11e-16` where
    /// the coarse mesh had it at exactly `0`, because
    ///
    /// ```text
    /// (origin + h·i) + h    ≠    origin + h·(i + 1)
    /// ```
    ///
    /// in IEEE at a spacing that is not a power of two. A hairline difference in
    /// the last bit is not a rounding curiosity here — it is a **crack**, and it
    /// costs the seam its shared vertices: R-004 measured 63–348 unmatched
    /// seam-plane boundary edges under a bit-identity merge at every spacing that
    /// is not a power of two, and, in two of twelve cases, a hole **1.05–2.08
    /// cells** wide where the perturbed sample crossed zero and the two sides
    /// stopped agreeing that an edge was cut at all (M-278).
    ///
    /// The crate's own weld closes the hairline — that much was written here and
    /// is ✗18 — but not the hole, and not for a consumer that never welds
    /// (M-69).
    ///
    /// Indexing from the grid origin is the same expression
    /// [`ChunkLayout::world_of_sample`](crate::chunk::ChunkLayout::world_of_sample)
    /// uses, so the coarse grid's `origin + (2h)·c` and this cell's
    /// `origin + h·(2c)` are bit-identical by M-70. See M-73.
    pub fn sample<S: Sdf<Scalar = R>>(
        sdf: &S,
        origin: [R; 3],
        step: R,
        base: [i64; 3],
        u: usize,
        v: usize,
        width: R,
    ) -> Self {
        debug_assert!(u < 3 && v < 3 && u != v, "u and v must be distinct axes");
        let normal = 3 - u - v;

        let mut value = [R::ZERO; SAMPLE_COUNT];
        let mut position = [[R::ZERO; 3]; SAMPLE_COUNT];
        let mut back = [[R::ZERO; 3]; 4];
        for s in 0..SAMPLE_COUNT {
            let mut index = base;
            index[u] += (s % 3) as i64;
            index[v] += (s / 3) as i64;
            let mut p = [R::ZERO; 3];
            for (axis, slot) in p.iter_mut().enumerate() {
                *slot = origin[axis] + step * R::from_f64(index[axis] as f64);
            }
            position[s] = p;
            value[s] = sdf.sample(p);

            if let Some(k) = back_slot(s as u8) {
                let mut q = p;
                q[normal] += width;
                back[k] = q;
            }
        }
        Self {
            value,
            position,
            back,
        }
    }

    /// Where a sample sits on the face this edge belongs to.
    ///
    /// A half-resolution edge reads its endpoints off the displaced face; every
    /// other edge reads them off the full-resolution one. This is the only place
    /// the two faces are told apart, and it is what turns the patch from a flat
    /// wall into a ribbon with an orientation.
    fn endpoint(&self, edge: u8, sample: u8) -> [R; 3] {
        if is_half_resolution(edge) {
            if let Some(k) = back_slot(sample) {
                return self.back[k];
            }
            debug_assert!(false, "a half-resolution edge joins two corners");
        }
        self.position[sample as usize]
    }

    /// Whether the surface crosses this edge.
    #[must_use]
    pub fn is_cut(&self, edge: u8) -> bool {
        let [a, b] = EDGE_SAMPLES[edge as usize];
        is_inside(self.value[a as usize]) != is_inside(self.value[b as usize])
    }

    /// Where the surface crosses `edge`, or `None` if it does not.
    ///
    /// `mid + (hi − lo)·d` with `d = ((a + b)/2)/(a − b)`, which is Marching
    /// Cubes' own placement (`cube::edge_offset`, R-059). On a cut edge exactly
    /// one endpoint is strictly negative and the other is `>= 0`, so `a − b` is
    /// never zero and no epsilon guard is needed — and an epsilon here would snap
    /// resolvable crossings to the midpoint, which is precisely the sub-voxel
    /// detail A-014 exists to keep.
    #[must_use]
    pub fn crossing(&self, edge: u8) -> Option<[R; 3]> {
        if !self.is_cut(edge) {
            return None;
        }
        let [lo, hi] = EDGE_SAMPLES[edge as usize];
        let (a, b) = (self.value[lo as usize], self.value[hi as usize]);
        let d = edge_offset(a, b);
        let (lo_pos, hi_pos) = (self.endpoint(edge, lo), self.endpoint(edge, hi));
        Some([
            place(lo_pos[0], hi_pos[0], d),
            place(lo_pos[1], hi_pos[1], d),
            place(lo_pos[2], hi_pos[2], d),
        ])
    }

    /// Every crossing, in edge order, with its edge index.
    pub fn crossings(&self) -> impl Iterator<Item = (u8, [R; 3])> + '_ {
        (0..EDGE_COUNT as u8).filter_map(|e| self.crossing(e).map(|p| (e, p)))
    }

    /// The case index this cell presents to [`table`](super::table).
    #[must_use]
    pub fn case(&self) -> u16 {
        let mut case = 0u16;
        for (s, v) in self.value.iter().enumerate() {
            if is_inside(*v) {
                case |= 1 << s;
            }
        }
        case
    }

    /// The centroid of this cell's crossings, or `None` when the surface misses
    /// it.
    ///
    /// Not a vertex rule — a transition cell's vertices are its crossings. This
    /// is here for A-011c, which needs a point inside the cell to fan a long
    /// cycle from, the same way `marching_cubes` fans one from a centroid rather
    /// than from one of its own vertices (A-015).
    #[must_use]
    pub fn centroid(&self) -> Option<[R; 3]> {
        let mut sum = [R::ZERO; 3];
        let mut n = 0u32;
        for (_, p) in self.crossings() {
            for (axis, slot) in sum.iter_mut().enumerate() {
                *slot += p[axis];
            }
            n += 1;
        }
        if n == 0 {
            return None;
        }
        Some(vec3::scale(sum, R::from_f64(f64::from(n)).recip()))
    }

    /// Triangulate this cell into `out`.
    ///
    /// One surface patch per cycle of
    /// [`super::table::transition_links`], fanned from the
    /// cycle's own centroid rather than from one of its vertices.
    ///
    /// **The centroid fan is A-015's rule, and it applies here for the same
    /// reason.** Fanning from a vertex leaves `k − 3` interior chords, and two
    /// cells sharing a face can pick the same chord and put four triangles on one
    /// mesh edge. A centroid is cell-local, so its spokes cannot be named by any
    /// other cell. Transition cells sit *between* two differently-resolved blocks,
    /// where a chord collision is likeliest, so the safe rule is the only one
    /// worth having.
    ///
    /// Winding follows the cycle direction the table produces, which is
    /// counter-clockwise seen from outside **in `(u, v, w)` parameter space**.
    /// The map to world space preserves that orientation only when
    /// `sign(width)` times the handedness of `(u, v, normal)` is positive, so
    /// the emitted order is swapped for the reflected parameterisations —
    /// Lengyel's tables store reversed windings for reflected transition cells
    /// for exactly this reason. `a_patch_with_width_is_wound_away_from_the_solid`
    /// and `a_mirrored_patch_is_wound_away_from_the_solid` establish both ways
    /// round come out facing outward, because no manifold or Euler check can
    /// see a global flip.
    ///
    /// `joined` selects, per ambiguous face, which pairing it uses; pass `0` for
    /// the separating choice, which is Marching Cubes proper.
    pub fn emit<S, M>(&self, sdf: &S, joined: u16, out: &mut M)
    where
        S: Sdf<Scalar = R>,
        M: MeshSink<Scalar = R>,
    {
        let case = self.case();
        debug_assert_eq!(
            joined & !AMBIGUOUS_FACES[case as usize],
            0,
            "a pairing was chosen for a face that is not ambiguous"
        );
        let next = transition_links(case, joined);

        // The FACES cycles are counter-clockwise seen from outside in
        // (u, v, w) *parameter* space. The map to world space preserves that
        // orientation only when sign(width) times the parity of
        // (u, v, normal) is positive -- Lengyel's tables store reversed
        // windings for the reflected transition cells for exactly this
        // reason. e_u, e_v and the back displacement each have exactly one
        // nonzero component, so this triple product is +-step^2*width with no
        // cancellation and its sign is exact. Zero width leaves it exactly
        // zero: no winding is decidable there (M-74), and the shipped order
        // is kept.
        let e_u = vec3::sub(self.position[1], self.position[0]);
        let e_v = vec3::sub(self.position[3], self.position[0]);
        let towards_back = vec3::sub(self.back[0], self.position[0]);
        let reversed = vec3::dot(vec3::cross(e_u, e_v), towards_back) < R::ZERO;

        let mut visited = 0u16;
        for start in 0..EDGE_COUNT as u8 {
            if next[start as usize] == NO_EDGE || visited & (1 << start) != 0 {
                continue;
            }

            // Walk the cycle, collecting its crossings in order.
            let mut cycle = [[R::ZERO; 3]; EDGE_COUNT];
            let mut len = 0usize;
            let mut current = start;
            while visited & (1 << current) == 0 {
                visited |= 1 << current;
                let Some(p) = self.crossing(current) else {
                    debug_assert!(false, "a linked edge must be cut");
                    break;
                };
                cycle[len] = p;
                len += 1;
                current = next[current as usize];
            }
            if len < 3 {
                debug_assert!(
                    false,
                    "a cycle closes across faces and has three or more edges"
                );
                continue;
            }

            let mut sum = [R::ZERO; 3];
            for p in &cycle[..len] {
                for (axis, slot) in sum.iter_mut().enumerate() {
                    *slot += p[axis];
                }
            }
            let centre = vec3::scale(sum, R::from_f64(len as f64).recip());
            let hub = out.vertex(centre, unit_gradient(sdf, centre));

            let mut spoke = [0u32; EDGE_COUNT];
            for (k, p) in cycle[..len].iter().enumerate() {
                spoke[k] = out.vertex(*p, unit_gradient(sdf, *p));
            }
            // `k` before `k + 1`, and the reason is measured — but it took two
            // goes to get a measurement that meant anything. At zero width the
            // patch is coplanar with the transition face and exactly
            // perpendicular to the surface (M-74), so a winding test against the
            // field gradient reported the same count either way round and
            // decided nothing. With Lengyel's width the patch becomes a ribbon
            // with a normal — best `|cos|` against the surface normal 1.000
            // rather than 0.000 — and the orientation is unanimous: this order
            // faces away from the solid on all 144 faces, the other order on
            // none.
            for k in 0..len {
                let (b, c) = (spoke[k], spoke[(k + 1) % len]);
                if reversed {
                    out.triangle(hub, c, b);
                } else {
                    out.triangle(hub, b, c);
                }
            }
        }
    }
}

/// The field's own gradient at a point, normalised — the crate's normal rule.
///
/// Written out rather than shared with `marching_cubes`' copy because that one is
/// private to its module; `normals::recompute` is the public path for anything
/// that wants a different rule.
fn unit_gradient<R: Real, S: Sdf<Scalar = R>>(sdf: &S, p: [R; 3]) -> [R; 3] {
    let g = sdf.gradient(p);
    let length = vec3::length(g);
    debug_assert!(length > R::ZERO, "zero gradient at a surface vertex");
    vec3::scale(g, length.recip())
}

#[cfg(test)]
mod tests;
