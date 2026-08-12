//! The machinery both dual methods share.
//!
//! A dual method places **one vertex per cell** the surface passes through, and
//! **one quad per crossed grid edge**, joining the four cells around that edge.
//! Surface Nets and Dual Contouring are both of these, and V-19 is explicit that
//! they are the *same* method twice over:
//!
//! > Dual Contouring's topology is Surface Nets' topology. The paper's algorithm
//! > is literally: vertex at the QEF minimizer for each sign-changing cube, quad
//! > joining the four cubes of each sign-changing edge. Only vertex *placement*
//! > differs.
//! >
//! > — Ju, Losasso, Schaefer & Warren 2002, `10.1145/566570.566586` §2.2, read
//! > and recorded as V-19.
//!
//! So the sampling, the smoothing, the vertex emission and the quad walk live
//! here once, and the two algorithms are a [`VertexRule`] each. Writing the quad
//! walk twice would mean two places for a winding bug to hide and two things to
//! fix when G-001 adds chunk seams.
//!
//! # What this does not decide
//!
//! Nothing here knows what a sharp feature is. The entire difference between a
//! rounded corner and a crisp one is which [`VertexRule`] is passed in, which is
//! the property that makes E-104 an honest comparison: same grid, same
//! crossings, same topology, same winding — one function swapped.

use alloc::vec::Vec;

use crate::cube::is_inside;
use crate::vec3;
use crate::{MeshSink, Real, Sdf, Shape3};

/// Where a dual method puts a cell's vertex.
///
/// Implementors get the cell's eight corner samples and its position on the
/// grid, and return a **world-space** point. Returning `None` means "no vertex
/// here", which the caller already believes — it has checked the corner signs —
/// so it is a safety valve rather than a control path.
///
/// The rule receives `base`, `origin` and `cell_size` separately rather than a
/// pre-computed cell origin so that an implementation can do its arithmetic in
/// whichever space keeps it exact.
pub(crate) trait VertexRule<R: Real> {
    /// Place this cell's vertex, given its corner samples in the crate's corner
    /// order.
    fn place<S: Sdf<Scalar = R>>(
        &self,
        sdf: &S,
        corner: &[R; 8],
        base: [u32; 3],
        origin: [R; 3],
        cell_size: R,
    ) -> Option<[R; 3]>;
}

/// Scratch and stages shared by every dual method.
///
/// Owns its buffers for the same reason
/// [`MarchingCubes`](crate::marching_cubes::MarchingCubes) does: the real workload re-meshes
/// thousands of chunks and allocation dominates.
#[derive(Debug)]
pub(crate) struct DualMesher<R: Real> {
    values: Vec<R>,
    /// Per cell: where its vertex sits, before it is handed to the sink.
    cell_position: Vec<[R; 3]>,
    /// Scratch for a smoothing pass, so smoothing is not biased by its own
    /// partial results.
    smoothed: Vec<[R; 3]>,
    /// Per cell: the index the sink gave that vertex, or [`u32::MAX`].
    cell_vertex: Vec<u32>,
    /// Per cell: whether the surface passes through it at all.
    cell_active: Vec<bool>,
    /// Laplacian passes. Surface Nets exposes this; Dual Contouring does not,
    /// because averaging a vertex with its neighbours is precisely what destroys
    /// the sharp feature the solve just recovered.
    pub(crate) smoothing_passes: u32,
}

impl<R: Real> DualMesher<R> {
    /// A mesher that has allocated nothing yet.
    pub(crate) const fn new() -> Self {
        Self {
            values: Vec::new(),
            cell_position: Vec::new(),
            smoothed: Vec::new(),
            cell_vertex: Vec::new(),
            cell_active: Vec::new(),
            smoothing_passes: 0,
        }
    }

    /// Sample, place, smooth, emit.
    ///
    /// # Errors
    ///
    /// [`Error::GridTooSmall`](crate::Error::GridTooSmall) if any axis has fewer
    /// than two samples: a dual method places at most one vertex per cell, and a
    /// grid with no cells has nothing to place.
    pub(crate) fn extract<S, M, V>(
        &mut self,
        rule: &V,
        sdf: &S,
        shape: &impl Shape3,
        origin: [R; 3],
        cell_size: R,
        out: &mut M,
    ) -> crate::Result<()>
    where
        S: Sdf<Scalar = R>,
        M: MeshSink<Scalar = R>,
        V: VertexRule<R>,
    {
        let size = shape.size();
        if size[0] < 2 || size[1] < 2 || size[2] < 2 {
            return Err(crate::Error::GridTooSmall { size });
        }

        self.sample(sdf, shape, origin, cell_size);

        let cells = [size[0] - 1, size[1] - 1, size[2] - 1];
        let cell_count = cells[0] as usize * cells[1] as usize * cells[2] as usize;

        self.cell_position.clear();
        self.cell_position.resize(cell_count, [R::ZERO; 3]);
        self.cell_active.clear();
        self.cell_active.resize(cell_count, false);
        self.cell_vertex.clear();
        self.cell_vertex.resize(cell_count, u32::MAX);

        self.place_vertices(rule, sdf, shape, cells, origin, cell_size);
        self.smooth(cells);
        self.emit_vertices(sdf, cells, out);
        self.emit_quads(shape, cells, out);

        Ok(())
    }

    fn sample<S: Sdf<Scalar = R>>(
        &mut self,
        sdf: &S,
        shape: &impl Shape3,
        origin: [R; 3],
        cell_size: R,
    ) {
        let size = shape.size();
        self.values.clear();
        self.values.reserve(shape.element_count());
        for z in 0..size[2] {
            for y in 0..size[1] {
                for x in 0..size[0] {
                    self.values.push(sdf.sample([
                        origin[0] + cell_size * R::from_f64(f64::from(x)),
                        origin[1] + cell_size * R::from_f64(f64::from(y)),
                        origin[2] + cell_size * R::from_f64(f64::from(z)),
                    ]));
                }
            }
        }
    }

    /// One vertex per cell the surface passes through, wherever `rule` says.
    fn place_vertices<S, V>(
        &mut self,
        rule: &V,
        sdf: &S,
        shape: &impl Shape3,
        cells: [u32; 3],
        origin: [R; 3],
        cell_size: R,
    ) where
        S: Sdf<Scalar = R>,
        V: VertexRule<R>,
    {
        for z in 0..cells[2] {
            for y in 0..cells[1] {
                for x in 0..cells[0] {
                    let base = [x, y, z];
                    let mut corner = [R::ZERO; 8];
                    let mut inside_count = 0u32;
                    for (c, slot) in corner.iter_mut().enumerate() {
                        let o = crate::cube::corner_offset(c as u8);
                        let s = shape.linearize([base[0] + o[0], base[1] + o[1], base[2] + o[2]]);
                        *slot = self.values[s as usize];
                        if is_inside(*slot) {
                            inside_count += 1;
                        }
                    }
                    if inside_count == 0 || inside_count == 8 {
                        continue;
                    }

                    let Some(position) = rule.place(sdf, &corner, base, origin, cell_size) else {
                        continue;
                    };
                    let index = cell_index(cells, base);
                    self.cell_position[index] = position;
                    self.cell_active[index] = true;
                }
            }
        }
    }

    fn smooth(&mut self, cells: [u32; 3]) {
        for _ in 0..self.smoothing_passes {
            self.smoothed.clear();
            self.smoothed.extend_from_slice(&self.cell_position);
            for z in 0..cells[2] {
                for y in 0..cells[1] {
                    for x in 0..cells[0] {
                        let index = cell_index(cells, [x, y, z]);
                        if !self.cell_active[index] {
                            continue;
                        }
                        let mut sum = self.cell_position[index];
                        let mut count = 1u32;
                        for axis in 0..3usize {
                            for step in [-1i64, 1] {
                                let mut neighbour = [i64::from(x), i64::from(y), i64::from(z)];
                                neighbour[axis] += step;
                                if neighbour[axis] < 0 || neighbour[axis] >= i64::from(cells[axis])
                                {
                                    continue;
                                }
                                let n = cell_index(
                                    cells,
                                    [
                                        neighbour[0] as u32,
                                        neighbour[1] as u32,
                                        neighbour[2] as u32,
                                    ],
                                );
                                if !self.cell_active[n] {
                                    continue;
                                }
                                for (k, slot) in sum.iter_mut().enumerate() {
                                    *slot += self.cell_position[n][k];
                                }
                                count += 1;
                            }
                        }
                        let inv = R::from_f64(f64::from(count)).recip();
                        self.smoothed[index] = [sum[0] * inv, sum[1] * inv, sum[2] * inv];
                    }
                }
            }
            core::mem::swap(&mut self.cell_position, &mut self.smoothed);
        }
    }

    fn emit_vertices<S, M>(&mut self, sdf: &S, cells: [u32; 3], out: &mut M)
    where
        S: Sdf<Scalar = R>,
        M: MeshSink<Scalar = R>,
    {
        let cell_count = cells[0] as usize * cells[1] as usize * cells[2] as usize;
        for index in 0..cell_count {
            if !self.cell_active[index] {
                continue;
            }
            let position = self.cell_position[index];
            let g = sdf.gradient(position);
            let len = vec3::length(g);
            debug_assert!(len > R::ZERO, "zero gradient at a surface vertex");
            let normal = vec3::scale(g, len.recip());
            self.cell_vertex[index] = out.vertex(position, normal);
        }
    }

    /// One quad per crossed grid edge, joining the four cells around it.
    ///
    /// The winding follows the sign direction along the edge, so the quad faces
    /// away from the solid. `meshed_sphere_has_positive_signed_volume` is what
    /// establishes the direction is the right way round — no manifold or Euler
    /// check can see a globally inverted surface.
    fn emit_quads<M: MeshSink<Scalar = R>>(
        &self,
        shape: &impl Shape3,
        cells: [u32; 3],
        out: &mut M,
    ) {
        let size = shape.size();
        for axis in 0..3usize {
            let u = (axis + 1) % 3;
            let v = (axis + 2) % 3;

            // The edge runs from `p` to `p + e_axis`, and all four surrounding
            // cells must exist, which bounds `p` on the other two axes.
            let mut p = [0u32; 3];
            for a in 0..size[axis] - 1 {
                for b in 1..cells[u] {
                    for c in 1..cells[v] {
                        p[axis] = a;
                        p[u] = b;
                        p[v] = c;

                        let s0 = shape.linearize(p);
                        let mut q = p;
                        q[axis] += 1;
                        let s1 = shape.linearize(q);

                        let inside0 = is_inside(self.values[s0 as usize]);
                        let inside1 = is_inside(self.values[s1 as usize]);
                        if inside0 == inside1 {
                            continue;
                        }

                        // Cell origins: p minus du along u, dv along v. Ordered
                        // counter-clockwise seen from +axis.
                        let mut quad = [0u32; 4];
                        for (slot, (du, dv)) in
                            quad.iter_mut().zip([(0, 0), (1, 0), (1, 1), (0, 1)])
                        {
                            let mut origin = p;
                            origin[u] -= du;
                            origin[v] -= dv;
                            let index = cell_index(cells, origin);
                            *slot = self.cell_vertex[index];
                        }
                        // Every one of the four cells contains this crossed edge,
                        // so every one of them is active and has a vertex.
                        debug_assert!(quad.iter().all(|i| *i != u32::MAX));

                        if inside0 {
                            out.triangle(quad[0], quad[1], quad[2]);
                            out.triangle(quad[0], quad[2], quad[3]);
                        } else {
                            out.triangle(quad[0], quad[2], quad[1]);
                            out.triangle(quad[0], quad[3], quad[2]);
                        }
                    }
                }
            }
        }
    }
}

#[inline]
fn cell_index(cells: [u32; 3], p: [u32; 3]) -> usize {
    p[0] as usize + cells[0] as usize * (p[1] as usize + cells[1] as usize * p[2] as usize)
}
