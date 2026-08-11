//! Surface Nets — the dual method.
//!
//! Marching Cubes puts vertices *on grid edges* and triangulates within each
//! cell. Surface Nets does the opposite: one vertex per cell that the surface
//! passes through, and one quad per grid edge that the surface crosses, joining
//! the four cells around that edge. Hence "dual" — vertices and faces swap roles.
//!
//! The consequences are worth stating precisely, because two of the three
//! things usually claimed for Surface Nets turn out to be measurably false:
//!
//! - **It does *not* produce fewer triangles.** Measured against Marching Cubes
//!   on four fields, the counts differ by exactly `2χ` — four triangles on a
//!   sphere, zero on a torus. That is forced by Euler, not luck: Marching Cubes
//!   places one vertex per crossed grid edge and Surface Nets emits two
//!   triangles per crossed grid edge, and any closed triangulated surface has
//!   `F = 2V − 2χ`. See `triangle_counts_track_marching_cubes_up_to_two_chi`.
//! - **Quads, so the connectivity is regular** — but not uniformly degree four
//!   once triangulated. Measured max degree on a sphere: 10, against Marching
//!   Cubes' 9.
//! - **Rounded corners, genuinely.** A vertex placed at the centroid of its
//!   cell's edge crossings cannot sit on a sharp feature, because an average of
//!   points on a corner's two faces lands between them. Measured on `box_exact`:
//!   the nearest vertex to the corner `(1,1,1)` is 1.15 cells away. That is not
//!   a defect to fix here — it is exactly what dual contouring changes, and it
//!   is the whole point of E-104.
//!
//! # The structural limit
//!
//! One vertex per cell means that where two sheets of the surface pass through
//! the same cell, they are forced to share a vertex and the result is
//! **non-manifold**. Measured: 48 non-manifold edges on the capped gyroid at
//! 49³, 15 on `fbm_terrain` at 33³. The literature review calls this out as "DC's
//! actual structural defect", fixed architecturally by vertex splitting rather
//! than by patching — which is A-010.
//!
//! # On the absence of published timings
//!
//! `docs/research/2026-08-11-meshing-speed-analysis.md` records that Surface
//! Nets and greedy meshing **have no credible published timings at all**,
//! despite being the two things game engines actually ship. Every number this
//! crate measures for it is therefore a reference rather than a comparison, and
//! should be reported with its hardware and grid attached.

use alloc::vec::Vec;

use crate::cube::{EDGE_CORNERS, corner_offset, edge_crossing, is_inside};
use crate::vec3;
use crate::{MeshSink, Real, Sdf, Shape3};

/// Surface Nets over a sampled grid.
///
/// Owns its scratch buffers, for the same reason
/// [`MarchingCubes`](crate::mc::MarchingCubes) does: the real workload re-meshes
/// thousands of chunks and allocation dominates.
///
/// # Example
///
/// ```
/// use isomesh::{MeshBuffer, RuntimeShape3};
/// use isomesh::fields::Sphere;
/// use isomesh::sn::SurfaceNets;
///
/// let mut sn = SurfaceNets::<f32>::new();
/// let mut out = MeshBuffer::<f32>::new();
/// let shape = RuntimeShape3::new([33; 3]);
/// sn.extract(&Sphere::<f32>::canonical(), &shape, [-2.0; 3], 0.125, &mut out);
///
/// assert!(out.triangle_count() > 0);
/// ```
#[derive(Debug)]
pub struct SurfaceNets<R: Real> {
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
    smoothing_passes: u32,
}

impl<R: Real> SurfaceNets<R> {
    /// A mesher that has allocated nothing yet.
    #[must_use]
    pub const fn new() -> Self {
        Self {
            values: Vec::new(),
            cell_position: Vec::new(),
            smoothed: Vec::new(),
            cell_vertex: Vec::new(),
            cell_active: Vec::new(),
            smoothing_passes: 0,
        }
    }

    /// How many Laplacian smoothing passes to run. Default zero.
    ///
    /// Each pass replaces every vertex with the average of itself and the
    /// vertices of the face-adjacent active cells. It visibly relaxes the
    /// staircase that a coarse grid produces.
    ///
    /// **It can move a vertex outside its own cell**, and vertices that leave
    /// their cells are how a dual method starts self-intersecting. The
    /// self-intersection counter is the way to see that happening; A-009 is the
    /// ticket that measures the clamp which prevents it.
    pub fn set_smoothing_passes(&mut self, passes: u32) {
        self.smoothing_passes = passes;
    }

    /// Extract the zero level set into `out`.
    ///
    /// `shape` counts **samples**, so `[n; 3]` spans `n - 1` cells per axis.
    /// `origin` is the world position of sample `[0, 0, 0]`.
    ///
    /// # Conventions
    ///
    /// Identical to Marching Cubes, deliberately, so the two can be compared on
    /// the same field without an adapter: negative is inside, zero is outside,
    /// winding is counter-clockwise seen from outside the solid, and normals come
    /// from the field's gradient.
    ///
    /// Output is triangles, two per quad, because [`MeshSink`] takes triangles.
    /// The quad structure is still visible in the output — each pair shares a
    /// diagonal.
    ///
    /// # Panics
    ///
    /// If `shape` has fewer than two samples on any axis.
    pub fn extract<S, M>(
        &mut self,
        sdf: &S,
        shape: &impl Shape3,
        origin: [R; 3],
        cell_size: R,
        out: &mut M,
    ) where
        S: Sdf<Scalar = R>,
        M: MeshSink<Scalar = R>,
    {
        let size = shape.size();
        assert!(
            size[0] >= 2 && size[1] >= 2 && size[2] >= 2,
            "surface nets needs at least two samples per axis, got {size:?}"
        );

        self.sample(sdf, shape, origin, cell_size);

        let cells = [size[0] - 1, size[1] - 1, size[2] - 1];
        let cell_count = cells[0] as usize * cells[1] as usize * cells[2] as usize;

        self.cell_position.clear();
        self.cell_position.resize(cell_count, [R::ZERO; 3]);
        self.cell_active.clear();
        self.cell_active.resize(cell_count, false);
        self.cell_vertex.clear();
        self.cell_vertex.resize(cell_count, u32::MAX);

        self.place_vertices(shape, cells, origin, cell_size);
        self.smooth(cells);
        self.emit_vertices(sdf, cells, out);
        self.emit_quads(shape, cells, out);
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

    /// One vertex per cell the surface passes through, at the centroid of that
    /// cell's edge crossings.
    ///
    /// The centroid is what makes this Surface Nets rather than dual contouring.
    /// It is cheap, it is always inside the cell — so the mesh cannot
    /// self-intersect through bad placement — and it cannot represent a sharp
    /// feature, because an average of points on a corner's two faces lands
    /// between them.
    fn place_vertices(
        &mut self,
        shape: &impl Shape3,
        cells: [u32; 3],
        origin: [R; 3],
        cell_size: R,
    ) {
        for z in 0..cells[2] {
            for y in 0..cells[1] {
                for x in 0..cells[0] {
                    let base = [x, y, z];
                    let mut corner = [R::ZERO; 8];
                    let mut inside_count = 0u32;
                    for (c, slot) in corner.iter_mut().enumerate() {
                        let o = corner_offset(c as u8);
                        let s = shape.linearize([base[0] + o[0], base[1] + o[1], base[2] + o[2]]);
                        *slot = self.values[s as usize];
                        if is_inside(*slot) {
                            inside_count += 1;
                        }
                    }
                    if inside_count == 0 || inside_count == 8 {
                        continue;
                    }

                    let mut sum = [R::ZERO; 3];
                    let mut crossings = 0u32;
                    for [lo, hi] in EDGE_CORNERS {
                        let (a, b) = (corner[lo as usize], corner[hi as usize]);
                        if is_inside(a) == is_inside(b) {
                            continue;
                        }
                        let t = edge_crossing(a, b);
                        let (lo_o, hi_o) = (corner_offset(lo), corner_offset(hi));
                        for axis in 0..3 {
                            let from = R::from_f64(f64::from(lo_o[axis]));
                            let to = R::from_f64(f64::from(hi_o[axis]));
                            sum[axis] += from + (to - from) * t;
                        }
                        crossings += 1;
                    }
                    // A cell with both signs present always has at least one cut
                    // edge: the cube graph is connected, so some edge joins an
                    // inside corner to an outside one.
                    debug_assert!(crossings > 0);

                    let inv = R::from_f64(f64::from(crossings)).recip();
                    let index = cell_index(cells, base);
                    self.cell_position[index] = [
                        origin[0] + cell_size * (R::from_f64(f64::from(x)) + sum[0] * inv),
                        origin[1] + cell_size * (R::from_f64(f64::from(y)) + sum[1] * inv),
                        origin[2] + cell_size * (R::from_f64(f64::from(z)) + sum[2] * inv),
                    ];
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

impl<R: Real> Default for SurfaceNets<R> {
    fn default() -> Self {
        Self::new()
    }
}

#[inline]
fn cell_index(cells: [u32; 3], p: [u32; 3]) -> usize {
    p[0] as usize + cells[0] as usize * (p[1] as usize + cells[1] as usize * p[2] as usize)
}

#[cfg(test)]
mod tests;
