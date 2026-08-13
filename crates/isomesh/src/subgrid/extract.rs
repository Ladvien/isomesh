//! The extractor — [`roots`](super::roots) into [`surface`](super::surface),
//! over a grid.
//!
//! This is where the two halves of subgrid Marching Tetrahedra meet: 1D root
//! finding replaces the sign test, and §3.2's reconstruction triangulates
//! whatever that finds, however many crossings an edge carries.
//!
//! # What it does that the others cannot
//!
//! [`MarchingTetrahedra`](crate::marching_tetrahedra) asks one question per edge
//! and gets one bit. This asks for every zero along that edge. The difference is
//! not incremental: M-67 measured that a sign test cannot distinguish **95.6%**
//! of the configurations a tetrahedron can be in, and A-005 measured
//! `thin_plate` — a plate 0.4 cells thick — returning **zero triangles** from
//! greedy quads, because no cell centre is inside it.
//!
//! # Conventions
//!
//! Identical to [`marching_cubes`](crate::marching_cubes) and
//! [`marching_tetrahedra`](crate::marching_tetrahedra), deliberately: sign
//! negative-inside, zero counts as outside, normals the field's own gradient,
//! winding counter-clockwise seen from outside the solid. The last of those is
//! imposed here rather than inherited — see
//! [`extract`](SubgridMarchingTetrahedra::extract).
//!
//! # Why every cell recomputes its neighbours' edges
//!
//! A grid edge shared by two cells has its roots found twice, and a tetrahedron
//! edge shared inside a cell likewise. That is deliberate for now and it is what
//! makes the result correct without a cache: both calls pass **bit-identical**
//! endpoints, because a corner's position is always `origin + cell_size · index`
//! computed the same way, and [`super::roots::all_roots`] is
//! deterministic for identical arguments. A cache keyed on the grid edge would
//! be faster and is the obvious optimisation, but it is an optimisation with a
//! correctness precondition, and this ticket owes a working extractor before a
//! fast one.

use alloc::vec::Vec;

use crate::cube::corner_offset;
use crate::marching_tetrahedra::table::{TET_EDGE_COUNT, TET_EDGES, TETS};
use crate::mesh::MeshSink;
use crate::real::Real;
use crate::sdf::Sdf;
use crate::shape::Shape3;
use crate::vec3;

use super::roots::all_roots;
use super::surface::{TetCrossings, TetPatch, Unfilled, fill};

/// Subgrid Marching Tetrahedra — Baktash, Gillespie & Crane,
/// `10.48550/arXiv.2606.00454`.
///
/// Holds its working buffers so a repeated extraction allocates nothing new,
/// per `CLAUDE.md` rule 6.
#[derive(Clone, Debug)]
pub struct SubgridMarchingTetrahedra<R: Real> {
    samples: u32,
    along: [Vec<R>; TET_EDGE_COUNT],
    patch: TetPatch<R>,
    index: Vec<u32>,
}

impl<R: Real> SubgridMarchingTetrahedra<R> {
    /// A new extractor sampling each tetrahedron edge `samples` times.
    ///
    /// `samples` is the 1D marching resolution, and it is the knob that decides
    /// which features exist: a pair of crossings closer together than
    /// `1 / samples` of an edge is invisible, exactly as a pair closer than the
    /// grid spacing is invisible to a sign-based method (§1.3). It is **not** a
    /// quality setting on the same axis as grid resolution — raising it resolves
    /// thinner features at the same triangle count, which is the entire point of
    /// the method.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidCellSize`](crate::Error::InvalidCellSize) if `samples` is
    /// zero, which would make every edge unsampled and every mesh empty.
    pub fn new(samples: u32) -> crate::Result<Self> {
        if samples == 0 {
            return Err(crate::Error::InvalidCellSize { value: 0.0 });
        }
        Ok(Self {
            samples,
            along: Default::default(),
            patch: TetPatch::new(),
            index: Vec::new(),
        })
    }

    /// The 1D sampling resolution this extractor was built with.
    #[must_use]
    pub const fn samples(&self) -> u32 {
        self.samples
    }

    /// Extract the zero level set into `out`.
    ///
    /// `shape` counts **samples**, so `[n; 3]` spans `n - 1` cells per axis, and
    /// `origin` is the world position of sample `[0, 0, 0]`. Note that unlike
    /// every other extractor here, the grid is used only for its *geometry* —
    /// the field is never sampled at the grid nodes, because node values are
    /// exactly the information this method replaces.
    ///
    /// # Winding, and why the output must be welded before it is judged
    ///
    /// Counter-clockwise seen from outside the solid, matching every other
    /// extractor here. That is **imposed rather than inherited**: §3.2 fixes
    /// each polygon's vertex order from its own boundary curve, which is
    /// consistent within a tetrahedron and carries no relationship to which side
    /// the field calls inside. Each triangle is therefore flipped, if needed, to
    /// agree with the gradient at its own centroid — per triangle and not per
    /// patch, because a sheet thinner than a cell puts two oppositely-facing
    /// surfaces inside one tetrahedron.
    ///
    /// **Vertices are emitted per tetrahedron and are not shared.** Before
    /// welding, the output is a triangle soup: every edge looks like a boundary
    /// edge and the mesh has no topology to check. Weld with
    /// [`Welder`](crate::weld::Welder) before applying
    /// [`validate_indexed`](crate::validate::validate_indexed) or measuring
    /// self-intersections — M-93 and M-96 are both consequences of forgetting
    /// that.
    ///
    /// # Errors
    ///
    /// [`Error::GridTooSmall`](crate::Error::GridTooSmall) if any axis has fewer
    /// than two samples. [`Error::SubgridUnfilled`](crate::Error::SubgridUnfilled)
    /// if a tetrahedron could not be triangulated — a defect rather than an
    /// unsupported input, since every case §3.2 defines is implemented.
    pub fn extract<S, M>(
        &mut self,
        sdf: &S,
        shape: &impl Shape3,
        origin: [R; 3],
        cell_size: R,
        out: &mut M,
    ) -> crate::Result<()>
    where
        S: Sdf<Scalar = R>,
        M: MeshSink<Scalar = R>,
    {
        let size = shape.size();
        if size[0] < 2 || size[1] < 2 || size[2] < 2 {
            return Err(crate::Error::GridTooSmall { size });
        }

        for z in 0..size[2] - 1 {
            for y in 0..size[1] - 1 {
                for x in 0..size[0] - 1 {
                    let cell = [x, y, z];
                    for t in 0..TETS.len() {
                        self.cell_tet(sdf, origin, cell_size, cell, t, out)?;
                    }
                }
            }
        }
        Ok(())
    }

    /// One tetrahedron of one cell.
    fn cell_tet<S, M>(
        &mut self,
        sdf: &S,
        origin: [R; 3],
        cell_size: R,
        cell: [u32; 3],
        t: usize,
        out: &mut M,
    ) -> crate::Result<()>
    where
        S: Sdf<Scalar = R>,
        M: MeshSink<Scalar = R>,
    {
        // Corner positions. Written as `origin + cell_size · index` with the
        // index formed once, so the two cells sharing a face compute the same
        // world position by the same expression — M-32's caveat is that equal
        // by algebra is not equal by IEEE, and only the same expression is safe.
        let mut corners = [[R::ZERO; 3]; 4];
        for (c, slot) in corners.iter_mut().enumerate() {
            let offset = corner_offset(TETS[t][c]);
            for axis in 0..3 {
                let index = f64::from(cell[axis]) + f64::from(offset[axis]);
                slot[axis] = origin[axis] + cell_size * R::from_f64(index);
            }
        }

        // Every zero along every edge. `TETS[t]` is ordered by inclusion, so a
        // tet edge always runs from the lower cube-corner index to the higher —
        // which makes the traversal direction a property of the grid rather than
        // of the tetrahedron, and is what lets two tetrahedra sharing an edge
        // agree bit-for-bit without consulting each other.
        let mut total = 0usize;
        for (e, slot) in self.along.iter_mut().enumerate() {
            slot.clear();
            let [lo, hi] = TET_EDGES[e];
            all_roots(
                corners[lo as usize],
                corners[hi as usize],
                sdf,
                self.samples,
                slot,
            );
            total += slot.len();
        }
        if total == 0 {
            return Ok(());
        }

        let mut borrowed: [&[R]; TET_EDGE_COUNT] = [&[]; TET_EDGE_COUNT];
        for (slot, v) in borrowed.iter_mut().zip(self.along.iter()) {
            *slot = v.as_slice();
        }
        let crossings = TetCrossings {
            corners,
            along: borrowed,
        };

        let unfilled = fill(&crossings, &mut self.patch).map_err(|_| {
            // `check` can only fail on unsorted or out-of-range parameters, and
            // `all_roots` produces neither — so this is this crate's bug, and it
            // is reported as one rather than swallowed.
            crate::Error::SubgridUnfilled {
                cell,
                tet: t as u8,
                reason: "crossings rejected as malformed",
            }
        })?;
        if unfilled != Unfilled::None {
            return Err(crate::Error::SubgridUnfilled {
                cell,
                tet: t as u8,
                reason: match unfilled {
                    Unfilled::SingleLoop => "single loop",
                    Unfilled::Subdivision => "subdivision",
                    Unfilled::NonNormalLoop => "non-normal loop",
                    Unfilled::NoPattern => "residual is not a (d1, d2) pattern",
                    Unfilled::Inconsistent => "curves disagree with the crossings",
                    Unfilled::None => unreachable!(),
                },
            });
        }

        // Into the sink, with the field's own gradient as each normal.
        self.index.clear();
        self.index.reserve(self.patch.positions.len());
        for (at, position) in self.patch.positions.iter().enumerate() {
            let g = sdf.gradient(*position);
            let length = vec3::length(g);
            // `!is_finite() || <= 0` rather than a negated `>`: NaN is excluded by
            // the finiteness test first, so the comparison is total by the time
            // it runs.
            if !length.is_finite() || length <= R::ZERO {
                return Err(crate::Error::DegenerateNormal { vertex: at as u64 });
            }
            let normal = vec3::scale(g, length.recip());
            self.index.push(out.vertex(*position, normal));
        }
        for tri in &self.patch.triangles {
            let (a, b, c) = (
                self.index.get(tri[0] as usize),
                self.index.get(tri[1] as usize),
                self.index.get(tri[2] as usize),
            );
            let (Some(a), Some(b), Some(c)) = (a, b, c) else {
                return Err(crate::Error::SubgridUnfilled {
                    cell,
                    tet: t as u8,
                    reason: "a triangle indexed a vertex the patch does not have",
                });
            };

            // Orientation. §3.2 fixes each polygon's vertex order from its own
            // boundary curve, which is consistent within a tetrahedron and
            // carries no relation to which side the field calls inside — so the
            // winding has to be imposed here, against the only thing that knows:
            // the gradient, which points away from the solid.
            //
            // Per triangle rather than per patch, because one tetrahedron can
            // carry sheets facing opposite ways — `thin_plate`'s two faces are
            // 0.4 cells apart and routinely land in the same cell.
            let (pa, pb, pc) = (
                self.patch.positions[tri[0] as usize],
                self.patch.positions[tri[1] as usize],
                self.patch.positions[tri[2] as usize],
            );
            let face = vec3::cross(vec3::sub(pb, pa), vec3::sub(pc, pa));
            let third = R::ONE / R::from_f64(3.0);
            let centroid = [
                (pa[0] + pb[0] + pc[0]) * third,
                (pa[1] + pb[1] + pc[1]) * third,
                (pa[2] + pb[2] + pc[2]) * third,
            ];
            let outward = vec3::dot(face, sdf.gradient(centroid));

            if outward < R::ZERO {
                out.triangle(*a, *c, *b);
            } else {
                // Includes the exactly-zero case, which is a triangle with no
                // area — §3.2's boundary disks emit those by construction
                // (V-21) and there is no orientation to choose for one. Left in
                // its original order rather than dropped, so the connectivity
                // §3.2 built stays intact.
                out.triangle(*a, *b, *c);
            }
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests;
