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
//! hole (M-69). The weld can close a seam it can *see*; it cannot invent a shared
//! vertex where the two sides disagree in the last bit.
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

use crate::cube::{edge_crossing, is_inside};
use crate::marching_cubes::table::NO_EDGE;
use crate::vec3;
use crate::{MeshSink, Real, Sdf};

use super::table::{AMBIGUOUS_FACES, EDGE_COUNT, EDGE_SAMPLES, SAMPLE_COUNT, transition_links};

/// One transition cell's nine sample values and their world positions.
///
/// The four half-resolution corners are *not* stored twice: they are samples
/// 0, 2, 6 and 8 of the same nine, which is exactly why the seam closes. See the
/// module docs.
#[derive(Clone, Copy, Debug)]
pub struct TransitionCell<R: Real> {
    /// Field value at each of the nine samples.
    pub value: [R; SAMPLE_COUNT],
    /// World position of each of the nine samples.
    pub position: [[R; 3]; SAMPLE_COUNT],
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
    /// the last bit is not a rounding curiosity here — it is a **crack**, and one
    /// no weld can close, because the weld can only merge vertices it can see are
    /// the same and these two are not.
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
    ) -> Self {
        let mut value = [R::ZERO; SAMPLE_COUNT];
        let mut position = [[R::ZERO; 3]; SAMPLE_COUNT];
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
        }
        Self { value, position }
    }

    /// Whether the surface crosses this edge.
    #[must_use]
    pub fn is_cut(&self, edge: u8) -> bool {
        let [a, b] = EDGE_SAMPLES[edge as usize];
        is_inside(self.value[a as usize]) != is_inside(self.value[b as usize])
    }

    /// Where the surface crosses `edge`, or `None` if it does not.
    ///
    /// `lo + (hi − lo)·t` with `t = a/(a − b)`, which is Marching Cubes' own
    /// placement. On a cut edge exactly one endpoint is strictly negative and the
    /// other is `>= 0`, so `a − b` is never zero and no epsilon guard is needed —
    /// and an epsilon here would snap resolvable crossings to the midpoint, which
    /// is precisely the sub-voxel detail A-014 exists to keep.
    #[must_use]
    pub fn crossing(&self, edge: u8) -> Option<[R; 3]> {
        if !self.is_cut(edge) {
            return None;
        }
        let [lo, hi] = EDGE_SAMPLES[edge as usize];
        let (a, b) = (self.value[lo as usize], self.value[hi as usize]);
        let t = edge_crossing(a, b);
        let (lo_pos, hi_pos) = (self.position[lo as usize], self.position[hi as usize]);
        Some([
            lo_pos[0] + (hi_pos[0] - lo_pos[0]) * t,
            lo_pos[1] + (hi_pos[1] - lo_pos[1]) * t,
            lo_pos[2] + (hi_pos[2] - lo_pos[2]) * t,
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
    /// Winding follows the cycle direction the table produces, which is oriented
    /// by walking each face counter-clockwise seen from outside the solid — the
    /// same construction Marching Cubes' own table uses. `a_transition_patch_is_wound_like_marching_cubes`
    /// is what establishes it comes out the right way round rather than inside
    /// out, because no manifold or Euler check can see a global flip.
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
            // `k + 1` before `k`, and the reason is measured rather than
            // guessed: emitted the other way round,
            // `a_transition_patch_is_wound_like_marching_cubes` reported **136 of
            // 136** faces looking into the solid. Unanimity is what makes this a
            // convention and not a bug — a mixed result would have meant the
            // table's cycle orientation was itself inconsistent, which would need
            // fixing rather than flipping.
            for k in 0..len {
                out.triangle(hub, spoke[(k + 1) % len], spoke[k]);
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
