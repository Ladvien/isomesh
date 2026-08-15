//! Marching Tetrahedra.
//!
//! Ticket: A-003.
//!
//! Splits each cell into six tetrahedra and contours each one independently.
//! Both the decomposition and the case table are **derived at compile time** —
//! see [`table`] for the construction and why there was nothing to transcribe
//! from even if transcription were wanted.
//!
//! # What this is for
//!
//! **It is the topological reference.** A tetrahedron's four values determine its
//! linear interpolant completely: there is no face saddle, no body saddle, and
//! no choice to make. Marching Cubes has an ambiguous face and A-002 spends a
//! whole ticket deciding it; Marching Tetrahedra cannot have one. So when two
//! extractors disagree about the topology of a field, this is the one that says
//! which answer the *sampling* supports, independent of any disambiguation rule.
//!
//! That is worth more here than it usually is, because this session found two
//! topology defects in Marching Cubes — a fan-chord collision (`✗17`) and a `χ`
//! difference under the decider — and both were caught by reasoning about one
//! algorithm against itself.
//!
//! # What it costs
//!
//! More triangles, and the amount is measured rather than repeated: see A-003's
//! archive entry. The literature's figure is **2–3× more vertices**
//! (`10.1109/2945.485620`, tier R), and Lewiner et al. (2003) make the stronger
//! claim that tetrahedral methods are *geometrically* worse and not merely
//! bulkier — *"the vertex position cannot be adjusted to fit the geometrical
//! trilinear approximation as we do with cubes"*. Both are checked against this
//! implementation rather than taken on.
//!
//! # Vertex sharing
//!
//! Marching Cubes places every vertex on a cube edge, so a `(sample, axis)` key
//! shares them. A tetrahedron's edges are the cube's twelve **plus** the six
//! face diagonals and the main diagonal — nineteen per cell. The key generalises
//! rather than changing shape: every tetrahedron edge runs from a corner to one
//! with a superset of its bits, so it is `(lower global sample, step)` where the
//! step is the `0/1` offset on each axis. Seven steps, so seven slots per
//! sample.
//!
//! The twelve axis-aligned edges are shared by the four cells around them and
//! the six face diagonals by the two cells across them, exactly as they should
//! be; the main diagonal is interior to its cell and shared by nothing. Nothing
//! special is done to arrange that — it falls out of the key being a statement
//! about world geometry rather than about a cell.

pub mod table;

#[cfg(test)]
mod tests;

use alloc::vec::Vec;

use crate::cube::{corner_offset, is_inside};
use crate::vec3;
use crate::{MeshSink, Real, Sdf, Shape3};

use table::{TET_CASES, TET_COUNT, TETS, tet_edge_corners};

/// How many distinct edge directions leave a sample: the seven non-zero `0/1`
/// steps on three axes.
const STEPS: usize = 7;

/// Marching Tetrahedra over a sampled grid.
///
/// Owns its scratch buffers, like every extractor here, so that re-meshing
/// thousands of chunks does not allocate thousands of times.
///
/// # Example
///
/// ```
/// use isomesh::{MeshBuffer, RuntimeShape3};
/// use isomesh::fields::Sphere;
/// use isomesh::marching_tetrahedra::MarchingTetrahedra;
///
/// let mut mt = MarchingTetrahedra::<f32>::new();
/// let mut out = MeshBuffer::<f32>::new();
///
/// let shape = RuntimeShape3::new([33; 3])?;
/// mt.extract(&Sphere::<f32>::canonical(), &shape, [-2.0; 3], 0.125, &mut out)?;
///
/// assert!(out.triangle_count() > 0);
/// # Ok::<(), isomesh::Error>(())
/// ```
#[derive(Debug)]
pub struct MarchingTetrahedra<R: Real> {
    values: Vec<R>,
    /// One slot per (sample, step): the vertex on that edge, or [`u32::MAX`].
    edge_vertices: Vec<u32>,
}

impl<R: Real> MarchingTetrahedra<R> {
    /// A mesher that has allocated nothing yet.
    #[must_use]
    pub const fn new() -> Self {
        Self {
            values: Vec::new(),
            edge_vertices: Vec::new(),
        }
    }

    /// Extract the zero level set into `out`.
    ///
    /// `shape` counts **samples**, so `[n; 3]` spans `n - 1` cells per axis.
    /// `origin` is the world position of sample `[0, 0, 0]`.
    ///
    /// # Conventions
    ///
    /// Identical to [`crate::marching_cubes`]'s, and deliberately so — sign
    /// negative-inside, zero counts as outside, winding counter-clockwise seen
    /// from outside the solid, normals the field's own gradient. Two extractors
    /// that disagreed about any of those would not be comparable, and comparing
    /// them is what this one is for.
    ///
    /// # Errors
    ///
    /// [`Error::GridTooSmall`](crate::Error::GridTooSmall) if any axis has fewer
    /// than two samples. [`Error::IndexSpaceExhausted`](crate::Error::IndexSpaceExhausted)
    /// if the grid could produce more vertices than a `u32` can address — seven
    /// per sample here rather than Marching Cubes' three, because a tetrahedron
    /// edge can be a face or body diagonal as well as a cube edge.
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
        let sample_count = shape.element_count();
        let bound = STEPS as u64 * sample_count as u64;
        if bound > u64::from(u32::MAX) {
            return Err(crate::Error::IndexSpaceExhausted { needed: bound });
        }

        self.values.clear();
        self.values.reserve(sample_count);
        for z in 0..size[2] {
            for y in 0..size[1] {
                for x in 0..size[0] {
                    let p = [
                        origin[0] + cell_size * R::from_f64(f64::from(x)),
                        origin[1] + cell_size * R::from_f64(f64::from(y)),
                        origin[2] + cell_size * R::from_f64(f64::from(z)),
                    ];
                    self.values.push(sdf.sample(p));
                }
            }
        }

        self.edge_vertices.clear();
        self.edge_vertices.resize(sample_count * STEPS, u32::MAX);

        for z in 0..size[2] - 1 {
            for y in 0..size[1] - 1 {
                for x in 0..size[0] - 1 {
                    let base = [x, y, z];

                    let mut corner_value = [R::ZERO; 8];
                    for (c, slot) in corner_value.iter_mut().enumerate() {
                        let s = corner_sample(shape, base, c as u8);
                        *slot = self.values[s as usize];
                    }

                    for t in 0..TET_COUNT {
                        let mut case = 0u8;
                        for i in 0..4 {
                            if is_inside(corner_value[TETS[t][i] as usize]) {
                                case |= 1 << i;
                            }
                        }
                        let entry = &TET_CASES[t][case as usize];
                        for tri in &entry.triangles[..entry.count as usize] {
                            let mut idx = [0u32; 3];
                            for (k, &edge) in tri.iter().enumerate() {
                                idx[k] = self.vertex_on_edge(
                                    sdf,
                                    shape,
                                    base,
                                    t,
                                    edge as usize,
                                    &corner_value,
                                    origin,
                                    cell_size,
                                    out,
                                );
                            }
                            out.triangle(idx[0], idx[1], idx[2]);
                        }
                    }
                }
            }
        }

        Ok(())
    }

    /// The vertex on one cut edge of one tetrahedron, creating it if this is the
    /// first cell to ask.
    ///
    /// Keyed on `(lower global sample, step)` rather than on anything local, so
    /// every cell that contains this edge agrees about it and the traversal
    /// order cannot affect the result.
    #[allow(clippy::too_many_arguments)]
    fn vertex_on_edge<S, M>(
        &mut self,
        sdf: &S,
        shape: &impl Shape3,
        base: [u32; 3],
        tet: usize,
        edge: usize,
        corner_value: &[R; 8],
        origin: [R; 3],
        cell_size: R,
        out: &mut M,
    ) -> u32
    where
        S: Sdf<Scalar = R>,
        M: MeshSink<Scalar = R>,
    {
        let [lo_corner, hi_corner] = tet_edge_corners(tet, edge);
        let lo = corner_offset(lo_corner);
        let hi = corner_offset(hi_corner);
        // Every tetrahedron edge runs from fewer bits to more, so this is a 0/1
        // step on each axis and the code is in 1..=7.
        let step = (hi[0] - lo[0]) + 2 * (hi[1] - lo[1]) + 4 * (hi[2] - lo[2]);
        debug_assert!(
            (1..=7).contains(&step),
            "tetrahedron edge is not a 0/1 step"
        );
        let lo_sample = corner_sample(shape, base, lo_corner);
        let key = lo_sample as usize * STEPS + (step - 1) as usize;

        let cached = self.edge_vertices[key];
        if cached != u32::MAX {
            return cached;
        }

        let a = corner_value[lo_corner as usize];
        let b = corner_value[hi_corner as usize];
        // As in Marching Cubes: on a cut edge one endpoint is strictly negative
        // and the other is >= 0, so `a - b` is never zero and no epsilon guard
        // is wanted.
        debug_assert!(is_inside(a) != is_inside(b));
        let t = a / (a - b);

        let lo_pos = corner_position(base, lo_corner, origin, cell_size);
        let hi_pos = corner_position(base, hi_corner, origin, cell_size);
        let position = [
            lo_pos[0] + (hi_pos[0] - lo_pos[0]) * t,
            lo_pos[1] + (hi_pos[1] - lo_pos[1]) * t,
            lo_pos[2] + (hi_pos[2] - lo_pos[2]) * t,
        ];

        let g = sdf.gradient(position);
        let len = vec3::length(g);
        debug_assert!(len > R::ZERO, "zero gradient at a surface crossing");
        let normal = vec3::scale(g, len.recip());

        let index = out.vertex(position, normal);
        self.edge_vertices[key] = index;
        index
    }
}

impl<R: Real> Default for MarchingTetrahedra<R> {
    fn default() -> Self {
        Self::new()
    }
}

#[inline]
fn corner_sample(shape: &impl Shape3, base: [u32; 3], corner: u8) -> u32 {
    let o = corner_offset(corner);
    shape.linearize([base[0] + o[0], base[1] + o[1], base[2] + o[2]])
}

#[inline]
fn corner_position<R: Real>(base: [u32; 3], corner: u8, origin: [R; 3], cell_size: R) -> [R; 3] {
    let o = corner_offset(corner);
    [
        origin[0] + cell_size * R::from_f64(f64::from(base[0] + o[0])),
        origin[1] + cell_size * R::from_f64(f64::from(base[1] + o[1])),
        origin[2] + cell_size * R::from_f64(f64::from(base[2] + o[2])),
    ]
}
