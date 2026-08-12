//! A-005 — blocky meshing with greedy quad merging.
//!
//! The budget end of the tradeoff table, and the comparison baseline for
//! triangle counts. Every other extractor here places vertices *on* the surface;
//! this one classifies whole cells as solid or empty and emits the axis-aligned
//! faces between them, so the output is a Minecraft surface rather than an
//! isosurface.
//!
//! # Source, and the same problem as A-001 and A-003
//!
//! Lysenko, *Meshing in a Minecraft Game* (2012) — a blog post with **no DOI**,
//! not in the local corpus. There is no case table to transcribe here, only a
//! procedure, and the procedure is stated in prose in the research doc: *"sweep
//! each of 6 face directions, per 2D slice walk +X for a run of identical
//! voxels, extend +Y holding that width, emit one quad."* That is implemented
//! directly and then checked, rather than taken on faith.
//!
//! # Occupancy
//!
//! A cell is solid when the field is negative **at its centre**. One sample per
//! cell, not eight corners — a blocky mesher's whole premise is that a cell is
//! one thing, and asking eight corners then reducing them to one bit would be
//! doing Marching Cubes' work and throwing the answer away.
//!
//! # The domain boundary is capped
//!
//! A cell outside the grid counts as empty, so a solid cell at the edge emits a
//! face there. This differs from every other extractor in the crate, which lets
//! the surface leave through the sides — and it is the honest behaviour for this
//! one, because a blocky mesh of a half-filled box *is* a closed box. It means
//! `fbm_terrain` comes back closed here and open under Marching Cubes, which is
//! a difference in what the two algorithms are, not a defect in either.
//!
//! # Vertices are split, deliberately, and the validity suite welds first
//!
//! Each quad carries its own four vertices. A cube corner has three faces
//! meeting at three different normals, so a shared vertex there would have to
//! average them and the result would not be blocky — the hard edge *is* the
//! output. The consequence is that the index buffer describes an open surface
//! even though the geometry is closed: no two quads share an edge, so every edge
//! is a boundary edge.
//!
//! That is reconciled rather than argued away. [`crate::weld`] exists, and the
//! validity tests weld before validating, which is what turns "closed as a
//! surface, open as an index buffer" into a checkable statement. A consumer that
//! wants a manifold mesh does the same thing and accepts the smoothed normals;
//! one that wants the blocky look uses the buffer as it comes.
//!
//! # Normals cost nothing
//!
//! They are `±1` on one axis, exactly, with no gradient evaluation anywhere.
//! Every other extractor here calls [`Sdf::gradient`] once per vertex.

#[cfg(test)]
mod tests;

use alloc::vec::Vec;

use crate::cube::is_inside;
use crate::{MeshSink, Real, Sdf, Shape3};

/// Blocky meshing with greedy quad merging.
///
/// Owns its occupancy grid and mask scratch so re-meshing thousands of chunks
/// does not allocate thousands of times.
///
/// # Example
///
/// ```
/// use isomesh::{MeshBuffer, RuntimeShape3};
/// use isomesh::fields::Sphere;
/// use isomesh::greedy_quads::GreedyQuads;
///
/// let mut greedy = GreedyQuads::<f32>::new();
/// let mut out = MeshBuffer::<f32>::new();
///
/// let shape = RuntimeShape3::new([33; 3])?;
/// greedy.extract(&Sphere::<f32>::canonical(), &shape, [-2.0; 3], 0.125, &mut out)?;
///
/// assert!(out.triangle_count() > 0);
/// # Ok::<(), isomesh::Error>(())
/// ```
#[derive(Debug)]
pub struct GreedyQuads<R: Real> {
    /// One bit per cell, `x` fastest.
    solid: Vec<bool>,
    /// Scratch for one 2D slice's face mask.
    mask: Vec<bool>,
    merge: Merge,
    _scalar: core::marker::PhantomData<R>,
}

/// Whether adjacent coplanar faces are merged into one quad.
///
/// [`Off`](Merge::Off) is **face culling**: one quad per visible cell face, no
/// merging. It is the baseline the published `2.76x` saving is quoted against,
/// and it exists here as a setting rather than as a second implementation for a
/// specific reason — the test that measured M-56 originally re-derived the
/// unmerged count with its own copy of the occupancy and visibility logic, which
/// is two paths that can drift apart while still agreeing on the day they were
/// written. E-106 needed the unmerged *mesh* rather than a count, and one switch
/// serves both.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum Merge {
    /// Merge greedily: widest run along `u`, then as many full-width rows along
    /// `v` as continue it. The default, and the algorithm's whole point.
    #[default]
    Greedy,
    /// One quad per visible face.
    Off,
}

impl<R: Real> GreedyQuads<R> {
    /// A mesher that has allocated nothing yet.
    #[must_use]
    pub const fn new() -> Self {
        Self {
            solid: Vec::new(),
            mask: Vec::new(),
            merge: Merge::Greedy,
            _scalar: core::marker::PhantomData,
        }
    }

    /// Whether coplanar faces are merged. Defaults to [`Merge::Greedy`].
    ///
    /// [`Merge::Off`] is the face-culling baseline; see [`Merge`] for why it is a
    /// setting rather than a separate function.
    pub fn set_merge(&mut self, merge: Merge) {
        self.merge = merge;
    }

    /// Extract the blocky surface into `out`.
    ///
    /// `shape` counts **samples**, so `[n; 3]` spans `n - 1` cells per axis and
    /// the occupancy grid is one smaller than the sample grid on every axis.
    /// `origin` is the world position of sample `[0, 0, 0]`.
    ///
    /// # Conventions
    ///
    /// Sign and winding are the crate's: negative is inside, and quads are wound
    /// counter-clockwise seen from outside the solid. Normals are exact
    /// axis-aligned unit vectors rather than a sampled gradient.
    ///
    /// # Errors
    ///
    /// [`Error::GridTooSmall`](crate::Error::GridTooSmall) if any axis has fewer
    /// than two samples, since then there is no cell to classify.
    /// [`Error::IndexSpaceExhausted`](crate::Error::IndexSpaceExhausted) if the
    /// output could exceed a `u32` index. Each cell can contribute at most six
    /// faces of four vertices, which is the bound checked here — greedy merging
    /// only ever reduces it.
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
        let cells = [size[0] - 1, size[1] - 1, size[2] - 1];
        let cell_count = cells[0] as u64 * cells[1] as u64 * cells[2] as u64;
        let bound = cell_count.saturating_mul(24);
        if bound > u64::from(u32::MAX) {
            return Err(crate::Error::IndexSpaceExhausted { needed: bound });
        }

        // ── occupancy: one sample per cell centre ───────────────────────────
        let half = cell_size * R::HALF;
        self.solid.clear();
        self.solid.reserve(cell_count as usize);
        for z in 0..cells[2] {
            for y in 0..cells[1] {
                for x in 0..cells[0] {
                    let p = [
                        origin[0] + cell_size * R::from_f64(f64::from(x)) + half,
                        origin[1] + cell_size * R::from_f64(f64::from(y)) + half,
                        origin[2] + cell_size * R::from_f64(f64::from(z)) + half,
                    ];
                    self.solid.push(is_inside(sdf.sample(p)));
                }
            }
        }

        let at = |c: [i64; 3]| -> bool {
            // Outside the grid is empty, so a solid cell at the edge is capped.
            for axis in 0..3 {
                if c[axis] < 0 || c[axis] >= i64::from(cells[axis]) {
                    return false;
                }
            }
            let i = c[0] as usize
                + cells[0] as usize * (c[1] as usize + cells[1] as usize * c[2] as usize);
            self.solid[i]
        };

        // ── six sweeps ──────────────────────────────────────────────────────
        for axis in 0..3usize {
            let u = (axis + 1) % 3;
            let v = (axis + 2) % 3;
            let du = cells[u] as usize;
            let dv = cells[v] as usize;

            for step in [1i64, -1] {
                for slice in 0..cells[axis] as i64 {
                    // A face exists where this cell is solid and the one it
                    // faces is not.
                    self.mask.clear();
                    self.mask.resize(du * dv, false);
                    for b in 0..dv {
                        for a in 0..du {
                            let mut c = [0i64; 3];
                            c[axis] = slice;
                            c[u] = a as i64;
                            c[v] = b as i64;
                            let mut n = c;
                            n[axis] = slice + step;
                            self.mask[a + du * b] = at(c) && !at(n);
                        }
                    }

                    // ── greedy merge ────────────────────────────────────────
                    //
                    // Walk the mask; at each set cell run along `u` for the
                    // widest span, then extend along `v` for as many rows as
                    // keep that full width. Clear what is consumed so nothing is
                    // emitted twice.
                    for b in 0..dv {
                        let mut a = 0usize;
                        while a < du {
                            if !self.mask[a + du * b] {
                                a += 1;
                                continue;
                            }
                            let mut width = 1usize;
                            let mut height = 1usize;
                            if self.merge == Merge::Greedy {
                                while a + width < du && self.mask[a + width + du * b] {
                                    width += 1;
                                }
                                'grow: while b + height < dv {
                                    for k in 0..width {
                                        if !self.mask[a + k + du * (b + height)] {
                                            break 'grow;
                                        }
                                    }
                                    height += 1;
                                }
                            }
                            for row in 0..height {
                                for k in 0..width {
                                    self.mask[a + k + du * (b + row)] = false;
                                }
                            }

                            emit_quad(
                                out, origin, cell_size, axis, u, v, slice, step, a, b, width,
                                height,
                            );
                            a += width;
                        }
                    }
                }
            }
        }

        Ok(())
    }
}

impl<R: Real> Default for GreedyQuads<R> {
    fn default() -> Self {
        Self::new()
    }
}

/// Write one merged quad, wound counter-clockwise seen from outside the solid.
///
/// The face sits at `slice` for a `-1` sweep and `slice + 1` for a `+1` one:
/// a cell spans `[slice, slice + 1]` on its axis, and the face is on whichever
/// side the empty neighbour is.
#[allow(clippy::too_many_arguments)]
fn emit_quad<R: Real, M: MeshSink<Scalar = R>>(
    out: &mut M,
    origin: [R; 3],
    cell_size: R,
    axis: usize,
    u: usize,
    v: usize,
    slice: i64,
    step: i64,
    a: usize,
    b: usize,
    width: usize,
    height: usize,
) {
    let plane = if step > 0 { slice + 1 } else { slice };
    let world = |along_u: usize, along_v: usize| -> [R; 3] {
        let mut p = [R::ZERO; 3];
        let mut grid = [0i64; 3];
        grid[axis] = plane;
        grid[u] = (a + along_u) as i64;
        grid[v] = (b + along_v) as i64;
        for k in 0..3 {
            p[k] = origin[k] + cell_size * R::from_f64(grid[k] as f64);
        }
        p
    };

    let mut normal = [R::ZERO; 3];
    normal[axis] = if step > 0 { R::ONE } else { -R::ONE };

    let corners = [
        world(0, 0),
        world(width, 0),
        world(width, height),
        world(0, height),
    ];

    // `(u, v, axis)` is a right-handed triple, so walking u then v is
    // counter-clockwise seen from `+axis`. A `-1` sweep faces the other way and
    // reverses. `a_meshed_box_has_positive_signed_volume` is what actually
    // checks this, since no topology test can see a global flip.
    let order = if step > 0 { [0, 1, 2, 3] } else { [0, 3, 2, 1] };
    let mut index = [0u32; 4];
    for (slot, &corner) in index.iter_mut().zip(order.iter()) {
        *slot = out.vertex(corners[corner], normal);
    }
    out.triangle(index[0], index[1], index[2]);
    out.triangle(index[0], index[2], index[3]);
}
