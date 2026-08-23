//! The machinery both dual methods share.
//!
//! A dual method places vertices **inside cells** and **one quad per crossed
//! grid edge**, joining the four cells around that edge. Surface Nets and Dual
//! Contouring are both of these, and V-19 is explicit that they are the *same*
//! method twice over:
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
//! # One vertex per cell is a property of the rule, not of the method
//!
//! A-010 widened this. The engine no longer assumes a cell holds one vertex: a
//! rule reports a set of **surface components**, each owning some of the cell's
//! twelve edges, and the quad walk asks each of the four cells *which of its
//! vertices owns this particular edge*. Surface Nets and Dual Contouring answer
//! "the only one" for every edge and their output is unchanged bit-for-bit;
//! Manifold Dual Contouring answers with the Marching Cubes cycle that contains
//! it, which is what makes its output manifold.
//!
//! The edge-indexed lookup is the whole mechanism. Keying on the cell alone is
//! precisely what makes two sheets through one cell share a vertex and pinch.
//!
//! # What this does not decide
//!
//! Nothing here knows what a sharp feature is. The entire difference between a
//! rounded corner and a crisp one is which [`VertexRule`] is passed in, which is
//! the property that makes E-104 an honest comparison: same grid, same
//! crossings, same topology, same winding — one function swapped.

use alloc::vec::Vec;

use crate::cube::{EDGE_COUNT, edge_index, is_inside};
use crate::vec3;
use crate::{MeshSink, Real, Sdf, Shape3};

/// The most vertices one cell can ever need.
///
/// A cell's surface components are the cycles of the Marching Cubes table, and
/// the most a cube configuration produces is four — the four corners of one
/// tetrahedron inside, each isolated from the others, each contributing its own
/// triangle. `every_case_fits_the_slot_budget` proves it over all 256 cases and
/// all 64 face-resolution masks rather than trusting the argument.
pub const MAX_CELL_VERTICES: usize = 4;

/// No vertex owns this edge, because the surface does not cross it.
pub(crate) const NO_SLOT: u8 = u8::MAX;

/// One cell's vertices, and which of its edges each one owns.
///
/// A [`VertexRule`] fills this in. The `slot_of_edge` half is what the quad walk
/// reads: for a crossed grid edge it needs *the* vertex of each surrounding cell
/// that lies on the same sheet of surface, and the cell alone does not identify
/// it.
#[derive(Debug)]
pub struct CellVertices<R: Real> {
    position: [[R; 3]; MAX_CELL_VERTICES],
    /// Which vertex owns each of the twelve edges, or [`NO_SLOT`].
    slot_of_edge: [u8; EDGE_COUNT],
    count: u8,
}

impl<R: Real> CellVertices<R> {
    pub(crate) const fn new() -> Self {
        Self {
            position: [[R::ZERO; 3]; MAX_CELL_VERTICES],
            slot_of_edge: [NO_SLOT; EDGE_COUNT],
            count: 0,
        }
    }

    /// Forget the previous cell. Called once per cell, before the rule runs.
    pub(crate) fn clear(&mut self) {
        self.slot_of_edge = [NO_SLOT; EDGE_COUNT];
        self.count = 0;
    }

    /// Place one vertex that owns the whole cell.
    ///
    /// The single-vertex answer, for rules that do not separate sheets. Every
    /// edge maps to it, including edges the surface does not cross — harmless,
    /// because the quad walk only ever asks about edges it has just found a sign
    /// change on.
    pub fn push_whole_cell(&mut self, position: [R; 3]) {
        debug_assert_eq!(self.count, 0, "a whole-cell vertex must be the only one");
        self.position[0] = position;
        self.slot_of_edge = [0; EDGE_COUNT];
        self.count = 1;
    }

    /// Place one vertex owning the edges in `edges`.
    ///
    /// # Panics
    ///
    /// In debug builds, if more than [`MAX_CELL_VERTICES`] components are pushed
    /// or two components claim the same edge.
    pub fn push_component(&mut self, position: [R; 3], edges: u16) {
        debug_assert!((self.count as usize) < MAX_CELL_VERTICES, "slot budget");
        let slot = self.count;
        self.position[slot as usize] = position;
        for (edge, owner) in self.slot_of_edge.iter_mut().enumerate() {
            if edges & (1 << edge) != 0 {
                debug_assert_eq!(*owner, NO_SLOT, "edge {edge} claimed twice");
                *owner = slot;
            }
        }
        self.count += 1;
    }

    fn count(&self) -> usize {
        self.count as usize
    }
}

/// Where a dual method puts a cell's vertices.
///
/// Implementors get the cell's eight corner samples and its position on the
/// grid, and fill `out` with **world-space** points. Leaving `out` empty means
/// "no vertex here", which the caller already believes it has ruled out — it has
/// checked the corner signs — so it is a safety valve rather than a control path.
///
/// The rule receives `base`, `origin` and `cell_size` separately rather than a
/// pre-computed cell origin so that an implementation can do its arithmetic in
/// whichever space keeps it exact.
pub trait VertexRule<R: Real> {
    /// Place this cell's vertices, given its corner samples in the crate's
    /// corner order.
    fn place<S: Sdf<Scalar = R>>(
        &self,
        sdf: &S,
        corner: &[R; 8],
        base: [u32; 3],
        origin: [R; 3],
        cell_size: R,
        out: &mut CellVertices<R>,
    );
}

/// The local edge a crossed grid edge becomes, seen from each of the four cells
/// around it.
///
/// Indexed `[axis][du][dv]`, matching the `(du, dv)` offsets [`DualMesher::emit_quads`]
/// subtracts from the edge's low grid point to reach each cell's origin. Inside
/// that cell the grid point sits at local offset `du` on `u` and `dv` on `v`, so
/// the edge runs from corner `(du << u) | (dv << v)` along `axis`.
const QUAD_EDGE: [[[u8; 2]; 2]; 3] = build_quad_edge();

const fn build_quad_edge() -> [[[u8; 2]; 2]; 3] {
    let mut out = [[[0u8; 2]; 2]; 3];
    let mut axis = 0usize;
    while axis < 3 {
        let u = (axis + 1) % 3;
        let v = (axis + 2) % 3;
        let mut du = 0usize;
        while du < 2 {
            let mut dv = 0usize;
            while dv < 2 {
                let lo = ((du as u8) << u) | ((dv as u8) << v);
                out[axis][du][dv] = edge_index(lo, lo | (1 << axis));
                dv += 1;
            }
            du += 1;
        }
        axis += 1;
    }
    out
}

/// Scratch and stages shared by every dual method.
///
/// Owns its buffers for the same reason
/// [`MarchingCubes`](crate::marching_cubes::MarchingCubes) does: the real workload re-meshes
/// thousands of chunks and allocation dominates.
///
/// The per-cell arrays are dense over the grid; the per-vertex ones are packed
/// over the cells that produced a vertex, because most cells produce none.
#[derive(Debug)]
pub(crate) struct DualMesher<R: Real> {
    /// Field values, on a grid whose **row length is forced odd** — see
    /// [`row_stride`] and A-024.
    values: Vec<R>,
    /// Samples per row in [`values`](Self::values), which is `size[0] | 1` and
    /// therefore **not** `size[0]` whenever the caller's grid is even.
    row: usize,
    /// One bit per sample, packed 64-to-a-word along `x`: is that sample inside?
    ///
    /// R-039. The active-cell test is *"are these eight corner signs not all
    /// equal"*, which is a question about eight bits, and the scalar form asked
    /// it with eight loads and eight comparisons for **every** cell while
    /// roughly 97% of them produce nothing. Packed, the same question is
    /// answered for sixty-four cells at once — see [`active_word`].
    ///
    /// The bit is `is_inside(value)`, the same comparison the scalar path made,
    /// and **not** the IEEE sign bit: `-0.0` has the sign bit set while
    /// `-0.0 < 0.0` is false, and `box_exact` is exactly zero across its whole
    /// boundary, so a sign-bit shortcut would change the mesh on a reference
    /// field.
    inside: Vec<u64>,
    /// Words per bitmap row, `size[0].div_ceil(64)`.
    ///
    /// Independent of [`row`](Self::row): that pads `values` to break A-024's
    /// stride aliasing, and this packs bits, so the two have no reason to agree
    /// and pretending they do would be a bug waiting for an even grid.
    bit_row: usize,
    /// Per cell: the index of its first vertex, or [`u32::MAX`] when the surface
    /// misses it. This doubles as the active flag.
    cell_first: Vec<u32>,
    /// Per cell: which of its vertices owns each of its twelve edges.
    cell_edge_slot: Vec<[u8; EDGE_COUNT]>,
    /// Per vertex: where it sits, before it is handed to the sink.
    slot_position: Vec<[R; 3]>,
    /// Scratch for a smoothing pass, so smoothing is not biased by its own
    /// partial results.
    smoothed: Vec<[R; 3]>,
    /// Per vertex: the index the sink gave it.
    slot_vertex: Vec<u32>,
    /// Scratch handed to the rule, reused across cells.
    scratch: CellVertices<R>,
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
            row: 0,
            inside: Vec::new(),
            bit_row: 0,
            cell_first: Vec::new(),
            cell_edge_slot: Vec::new(),
            slot_position: Vec::new(),
            smoothed: Vec::new(),
            slot_vertex: Vec::new(),
            scratch: CellVertices::new(),
            smoothing_passes: 0,
        }
    }

    /// Sample, place, smooth, emit.
    ///
    /// # Errors
    ///
    /// [`Error::GridTooSmall`](crate::Error::GridTooSmall) if any axis has fewer
    /// than two samples: a dual method places its vertices inside cells, and a
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

        self.cell_first.clear();
        self.cell_first.resize(cell_count, u32::MAX);
        self.cell_edge_slot.clear();
        self.cell_edge_slot
            .resize(cell_count, [NO_SLOT; EDGE_COUNT]);
        self.slot_position.clear();
        self.slot_vertex.clear();

        self.place_vertices(rule, sdf, shape, cells, origin, cell_size);
        self.smooth(cells, origin, cell_size);
        self.emit_vertices(sdf, out);
        self.emit_quads(shape, cells, out);

        Ok(())
    }

    /// Where sample `p` lives in [`values`](Self::values).
    ///
    /// **Not `shape.linearize`, and the difference is one bit (A-024, M-287).**
    /// `values` is `size[0]·size[1]·size[2]` floats laid out by the caller's
    /// shape, so its row stride is `size[0]·4` bytes and its plane stride
    /// `size[0]·size[1]·4`. At `size[0] = size[1] = 128` those are **512 bytes
    /// and exactly 64 KiB**, which are a cache-set aliasing period on ordinary
    /// hardware twice over, and Surface Nets measured **3.37× the cost of 127³
    /// or 129³** there — on a field with no surface, so it is the scaffolding
    /// and not the geometry. 256³ pays 1.39× for the same reason.
    ///
    /// Forcing the row length **odd** removes both periods at once and cannot
    /// reintroduce either: `4·odd` is never a multiple of 512, and
    /// `4·odd·size[1]` is a multiple of 65,536 only if `size[1]` is a multiple
    /// of 16,384, which is a grid of 2.7×10¹² samples.
    ///
    /// It is unconditional on purpose. A pad applied only when the stride looks
    /// bad would be a second layout reachable from the same call, which is the
    /// shape `CLAUDE.md`'s one-path rule exists to forbid — and a *fixed* pad of
    /// one would be worse than nothing, since it maps every `size[0] = 2ᵏ − 1`
    /// onto the stride it is trying to avoid. `| 1` is idempotent, so it has no
    /// such image.
    ///
    /// The cost is one float per row when the row is even — **0.8% of `values`
    /// at 128³** — and nothing at all at run time: the multiply is by a
    /// different constant, not an extra one.
    #[inline]
    fn index(&self, p: [u32; 3], size: [u32; 3]) -> usize {
        p[0] as usize + self.row * (p[1] as usize + size[1] as usize * p[2] as usize)
    }

    /// Samples per row for a grid of this size.
    #[inline]
    fn row_stride(size: [u32; 3]) -> usize {
        size[0] as usize | 1
    }

    fn sample<S: Sdf<Scalar = R>>(
        &mut self,
        sdf: &S,
        shape: &impl Shape3,
        origin: [R; 3],
        cell_size: R,
    ) {
        let size = shape.size();
        self.row = Self::row_stride(size);
        // One slot per row is padding and is never read; it is filled rather
        // than skipped so the buffer has no uninitialised gaps and `push` stays
        // a single sequential write.
        let pad = self.row - size[0] as usize;
        self.values.clear();
        self.values
            .reserve(self.row * size[1] as usize * size[2] as usize);
        for z in 0..size[2] {
            for y in 0..size[1] {
                for x in 0..size[0] {
                    self.values.push(sdf.sample([
                        origin[0] + cell_size * R::from_f64(f64::from(x)),
                        origin[1] + cell_size * R::from_f64(f64::from(y)),
                        origin[2] + cell_size * R::from_f64(f64::from(z)),
                    ]));
                }
                for _ in 0..pad {
                    self.values.push(R::ZERO);
                }
            }
        }
    }

    /// Pack `is_inside` into [`inside`](Self::inside), one bit per sample.
    ///
    /// R-039. One pass of `n³` comparisons, which buys the removal of `8n³`
    /// loads-and-comparisons in [`place_vertices`](Self::place_vertices) on
    /// every cell that produces nothing — and on a sphere at 128³ that is about
    /// 97% of them.
    fn build_inside_bits(&mut self, size: [u32; 3]) {
        let sx = size[0] as usize;
        let rows = size[1] as usize * size[2] as usize;
        self.bit_row = sx.div_ceil(64);
        self.inside.clear();
        self.inside.resize(self.bit_row * rows, 0);

        for row in 0..rows {
            let src = self.row * row;
            let dst = self.bit_row * row;
            for w in 0..self.bit_row {
                let base = w * 64;
                let n = (sx - base).min(64);
                let mut word = 0u64;
                for k in 0..n {
                    // Branchless on purpose: the whole point of the bitmap is
                    // that the sign test stops being a branch per corner.
                    word |= u64::from(is_inside(self.values[src + base + k])) << k;
                }
                self.inside[dst + w] = word;
            }
        }
    }

    /// Bit `k` of the returned word is `is_inside(sample[64w + k, y, z])`.
    #[inline]
    fn inside_word(&self, w: usize, y: usize, z: usize, size: [u32; 3]) -> u64 {
        self.inside[self.bit_row * (y + size[1] as usize * z) + w]
    }

    /// Bit `k` of the returned word is `is_inside(sample[64w + k + 1, y, z])`.
    ///
    /// The high bit has to come from the next word or the cell straddling a
    /// word boundary reads its `+x` corner as outside — which is a hole every
    /// 64 cells, invisible on a small grid and catastrophic on a large one.
    #[inline]
    fn inside_word_shifted(&self, w: usize, y: usize, z: usize, size: [u32; 3]) -> u64 {
        let lo = self.inside_word(w, y, z, size);
        let hi = if w + 1 < self.bit_row {
            self.inside_word(w + 1, y, z, size)
        } else {
            0
        };
        (lo >> 1) | (hi << 63)
    }

    /// Sixty-four active-cell answers, in four fused word operations per row.
    ///
    /// A cell is active exactly when its eight corner signs are **not all
    /// equal**. Reading `a` as "is `x` inside" and `b` as "is `x+1` inside" for
    /// each of the four rows `(y,z)`, `(y+1,z)`, `(y,z+1)`, `(y+1,z+1)`:
    ///
    /// ```text
    /// any    = OR  over the four rows of (a | b)      any corner inside
    /// all    = AND over the four rows of (a & b)      every corner inside
    /// active = any & !all                             mixed, i.e. crossed
    /// ```
    ///
    /// Source for the mechanism: Museth, *VDB*, `10.1145/2487228.2487235`, whose
    /// topology layer is per-node bitmasks under word-level boolean ops. The
    /// transfer is to the active-cell predicate of a CPU isosurface extractor,
    /// with the mesh held bit-identical because the bit is the same comparison
    /// and the set-bit walk below visits cells in the same lexicographic order
    /// the scalar loop did.
    #[inline]
    fn active_word(&self, w: usize, y: usize, z: usize, size: [u32; 3]) -> u64 {
        let mut any = 0u64;
        let mut all = !0u64;
        for dz in 0..2usize {
            for dy in 0..2usize {
                let a = self.inside_word(w, y + dy, z + dz, size);
                let b = self.inside_word_shifted(w, y + dy, z + dz, size);
                any |= a | b;
                all &= a & b;
            }
        }
        any & !all
    }

    /// Which bits of word `w` are cells that exist.
    ///
    /// The last word of a row runs past the final cell whenever `cells` is not a
    /// multiple of 64, and those bits are answers about samples that are not
    /// cells. `1u64 << 64` is undefined, so the full-word case is named rather
    /// than computed.
    #[inline]
    fn cell_mask(w: usize, cells_x: usize) -> u64 {
        let base = w * 64;
        let remaining = cells_x.saturating_sub(base);
        if remaining >= 64 {
            !0
        } else {
            (1u64 << remaining) - 1
        }
    }

    /// A vertex per surface component of every cell the surface passes through,
    /// wherever `rule` says.
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
        let size = shape.size();
        self.build_inside_bits(size);
        let cells_x = cells[0] as usize;
        // Words that carry a **cell**, which is not the same as words that carry
        // a **sample**. `bit_row` is `size[0].div_ceil(64)` because the bitmap is
        // indexed by sample; the cell row is one shorter, so at `size[0] = 64k+1`
        // — 65, 129, every grid this crate benchmarks at — the sample row needs
        // one more word than the cell row does, and iterating `bit_row` fuses a
        // whole `u64` that `cell_mask` then discards in its entirety.
        //
        // Measured by E-307 before this line existed: packed cost went 0.451–0.470
        // ns/cell at 128 samples to 0.599–0.637 at 129, non-overlapping ranges on a
        // field whose active fraction moved by 0.08 of a percentage point, taking
        // the stage ratio from 5.66–5.76× down to 4.13–4.56×. About 30% of the
        // stage, for a word that could not contain a cell.
        let cell_words = cells_x.div_ceil(64);
        debug_assert!(cell_words <= self.bit_row);

        for z in 0..cells[2] as usize {
            for y in 0..cells[1] as usize {
                for w in 0..cell_words {
                    let mut active = self.active_word(w, y, z, size) & Self::cell_mask(w, cells_x);
                    // `w &= w - 1` clears the lowest set bit, so this visits the
                    // active cells of the row in ascending `x` — the same order
                    // the scalar loop did, which is what keeps vertex creation
                    // order and therefore every index unchanged.
                    while active != 0 {
                        let x = (w * 64 + active.trailing_zeros() as usize) as u32;
                        active &= active - 1;

                        let base = [x, y as u32, z as u32];
                        let mut corner = [R::ZERO; 8];
                        for (c, slot) in corner.iter_mut().enumerate() {
                            let o = crate::cube::corner_offset(c as u8);
                            let s =
                                self.index([base[0] + o[0], base[1] + o[1], base[2] + o[2]], size);
                            *slot = self.values[s];
                        }
                        debug_assert!(
                            {
                                let n = corner.iter().filter(|v| is_inside(**v)).count();
                                n != 0 && n != 8
                            },
                            "the bitmap named an inactive cell at {base:?}"
                        );

                        self.scratch.clear();
                        rule.place(sdf, &corner, base, origin, cell_size, &mut self.scratch);
                        if self.scratch.count() == 0 {
                            continue;
                        }
                        // Every cut edge must have an owning vertex, or the quad
                        // walk would have no corner to use for it.
                        #[cfg(debug_assertions)]
                        {
                            for (edge, [lo, hi]) in
                                crate::cube::EDGE_CORNERS.into_iter().enumerate()
                            {
                                if is_inside(corner[lo as usize]) != is_inside(corner[hi as usize])
                                {
                                    assert_ne!(
                                        self.scratch.slot_of_edge[edge], NO_SLOT,
                                        "cut edge {edge} has no owning vertex"
                                    );
                                }
                            }
                        }

                        let index = cell_index(cells, base);
                        // Truncation is impossible for any grid that fits in
                        // memory: the sink's own index space runs out first.
                        self.cell_first[index] = self.slot_position.len() as u32;
                        self.cell_edge_slot[index] = self.scratch.slot_of_edge;
                        self.slot_position
                            .extend_from_slice(&self.scratch.position[..self.scratch.count()]);
                    }
                }
            }
        }
    }

    /// Laplacian relaxation over the face-adjacent cells.
    ///
    /// **Only meaningful for a rule that places one vertex per cell**, which is
    /// the only kind that exposes it: with two sheets in a cell there is no
    /// single answer to "the neighbour's vertex", and averaging across sheets
    /// would drag one into the other. Asserted rather than branched on.
    fn smooth(&mut self, cells: [u32; 3], origin: [R; 3], cell_size: R) {
        if self.smoothing_passes == 0 {
            return;
        }
        debug_assert_eq!(
            self.slot_position.len(),
            self.cell_first.iter().filter(|f| **f != u32::MAX).count(),
            "smoothing requires one vertex per active cell"
        );

        for _ in 0..self.smoothing_passes {
            self.smoothed.clear();
            self.smoothed.extend_from_slice(&self.slot_position);
            for z in 0..cells[2] {
                for y in 0..cells[1] {
                    for x in 0..cells[0] {
                        let index = cell_index(cells, [x, y, z]);
                        let first = self.cell_first[index];
                        if first == u32::MAX {
                            continue;
                        }
                        let mut sum = self.slot_position[first as usize];
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
                                let n_first = self.cell_first[n];
                                if n_first == u32::MAX {
                                    continue;
                                }
                                for (k, slot) in sum.iter_mut().enumerate() {
                                    *slot += self.slot_position[n_first as usize][k];
                                }
                                count += 1;
                            }
                        }
                        let inv = R::from_f64(f64::from(count)).recip();
                        // Gibson's box constraint, the load-bearing half of the
                        // relaxation: "each node clamped inside its original
                        // cube" (10.1007/bfb0056277). Without it the averaged
                        // vertex leaves its cell, the cells stop partitioning
                        // space, and the mesh self-intersects — the catalog's
                        // headline correction is that the mechanism is this box,
                        // not gradients. Same clamp policy as the QEF path, so
                        // the two cannot disagree about what "inside the cell"
                        // means.
                        let cell_origin = [
                            origin[0] + cell_size * R::from_f64(f64::from(x)),
                            origin[1] + cell_size * R::from_f64(f64::from(y)),
                            origin[2] + cell_size * R::from_f64(f64::from(z)),
                        ];
                        self.smoothed[first as usize] = crate::dual_contouring::apply_clamp(
                            crate::dual_contouring::Clamp::ToCell,
                            [sum[0] * inv, sum[1] * inv, sum[2] * inv],
                            cell_origin,
                            cell_size,
                        );
                    }
                }
            }
            core::mem::swap(&mut self.slot_position, &mut self.smoothed);
        }
    }

    fn emit_vertices<S, M>(&mut self, sdf: &S, out: &mut M)
    where
        S: Sdf<Scalar = R>,
        M: MeshSink<Scalar = R>,
    {
        self.slot_vertex.reserve(self.slot_position.len());
        for &position in &self.slot_position {
            let g = sdf.gradient(position);
            let len = vec3::length(g);
            debug_assert!(len > R::ZERO, "zero gradient at a surface vertex");
            let normal = vec3::scale(g, len.recip());
            self.slot_vertex.push(out.vertex(position, normal));
        }
    }

    /// One quad per crossed grid edge, joining the four cells around it.
    ///
    /// Each cell contributes **the vertex that owns this edge**, not simply its
    /// vertex — see the module docs. The winding follows the sign direction along
    /// the edge, so the quad faces away from the solid.
    /// `meshed_sphere_has_positive_signed_volume` is what establishes the
    /// direction is the right way round — no manifold or Euler check can see a
    /// globally inverted surface.
    fn emit_quads<M: MeshSink<Scalar = R>>(
        &self,
        shape: &impl Shape3,
        cells: [u32; 3],
        out: &mut M,
    ) {
        // Three monomorphisations of one function, not three copies of a loop.
        //
        // **The axis has to be a constant, and A-023 measured what it costs when
        // it is not (M-285).** With `axis`, `u` and `v` as runtime values, every
        // `p[axis] = a` is a dynamically indexed store, so `p` cannot live in
        // registers: each iteration writes three coordinates to the stack and
        // `linearize` reads them straight back, a store-to-load chain the
        // scheduler cannot break. This stage was **82% of the dual mesher's
        // cycles at IPC 0.72** (M-284) while the cell loop beside it, doing more
        // work per iteration, ran at 3.83.
        //
        // The emission order is unchanged — same three passes in the same order,
        // same loop bounds, same triangles in the same sequence — which is why
        // T-007's golden hashes are untouched by this.
        self.emit_quad_axis::<0, M>(shape, cells, out);
        self.emit_quad_axis::<1, M>(shape, cells, out);
        self.emit_quad_axis::<2, M>(shape, cells, out);
    }

    /// One axis of [`emit_quads`](Self::emit_quads), with the axis a constant.
    fn emit_quad_axis<const AXIS: usize, M: MeshSink<Scalar = R>>(
        &self,
        shape: &impl Shape3,
        cells: [u32; 3],
        out: &mut M,
    ) {
        let size = shape.size();
        {
            let axis = AXIS;
            let u = (AXIS + 1) % 3;
            let v = (AXIS + 2) % 3;

            // The edge runs from `p` to `p + e_axis`, and all four surrounding
            // cells must exist, which bounds `p` on the other two axes.
            let mut p = [0u32; 3];
            for a in 0..size[axis] - 1 {
                for b in 1..cells[u] {
                    for c in 1..cells[v] {
                        p[axis] = a;
                        p[u] = b;
                        p[v] = c;

                        let s0 = self.index(p, size);
                        let mut q = p;
                        q[axis] += 1;
                        let s1 = self.index(q, size);

                        let inside0 = is_inside(self.values[s0]);
                        let inside1 = is_inside(self.values[s1]);
                        if inside0 == inside1 {
                            continue;
                        }

                        // Cell origins: p minus du along u, dv along v. Ordered
                        // counter-clockwise seen from +axis.
                        let mut quad = [0u32; 4];
                        for (slot, (du, dv)) in
                            quad.iter_mut().zip([(0u32, 0u32), (1, 0), (1, 1), (0, 1)])
                        {
                            let mut cell = p;
                            cell[u] -= du;
                            cell[v] -= dv;
                            let index = cell_index(cells, cell);
                            let local = QUAD_EDGE[axis][du as usize][dv as usize];
                            let first = self.cell_first[index];
                            let owner = self.cell_edge_slot[index][local as usize];
                            // Every one of the four cells contains this crossed
                            // edge, so every one is active and some component of
                            // it owns the edge. A rule that declined to place a
                            // vertex leaves [`u32::MAX`] here — an index no sink
                            // ever returned, which the validator reports as an
                            // out-of-range index rather than reading past the
                            // end of the vertex list. Branching on the sentinel
                            // rather than letting the sum run past the list is
                            // what keeps that true on 32-bit targets, where
                            // `u32::MAX as usize` wraps back into range and
                            // would silently stitch the quad to an unrelated
                            // vertex.
                            debug_assert!(first != u32::MAX);
                            debug_assert!(owner != NO_SLOT);
                            *slot = if first == u32::MAX || owner == NO_SLOT {
                                u32::MAX
                            } else {
                                self.slot_vertex
                                    .get(first as usize + owner as usize)
                                    .copied()
                                    .unwrap_or(u32::MAX)
                            };
                        }

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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cube::{EDGE_AXIS, EDGE_CORNERS};

    /// The four cells around a crossed grid edge must each name *that* edge, or
    /// the quad walk would look up the wrong component's vertex. A transposition
    /// here would be invisible on any single-vertex rule and wrong on every
    /// multi-vertex one.
    #[test]
    fn quad_edge_names_the_shared_edge_from_all_four_cells() {
        for (axis, per_axis) in QUAD_EDGE.iter().enumerate() {
            let u = (axis + 1) % 3;
            let v = (axis + 2) % 3;
            for (du, per_du) in per_axis.iter().enumerate() {
                for (dv, &edge) in per_du.iter().enumerate() {
                    assert_eq!(
                        EDGE_AXIS[edge as usize] as usize, axis,
                        "axis {axis} du {du} dv {dv}"
                    );
                    // The cell sits at -du on u and -dv on v from the grid point,
                    // so within the cell the grid point is at +du, +dv.
                    let lo = EDGE_CORNERS[edge as usize][0];
                    assert_eq!(usize::from((lo >> u) & 1), du);
                    assert_eq!(usize::from((lo >> v) & 1), dv);
                    assert_eq!((lo >> axis) & 1, 0, "the low corner starts the edge");
                }
            }
        }
    }

    #[test]
    fn a_whole_cell_vertex_owns_every_edge() {
        let mut cell = CellVertices::<f64>::new();
        cell.clear();
        cell.push_whole_cell([1.0, 2.0, 3.0]);
        assert_eq!(cell.count(), 1);
        assert!(cell.slot_of_edge.iter().all(|s| *s == 0));
    }

    #[test]
    fn components_own_disjoint_edges() {
        let mut cell = CellVertices::<f64>::new();
        cell.clear();
        cell.push_component([0.0; 3], 0b0000_0000_0111);
        cell.push_component([1.0; 3], 0b0000_0111_0000);
        assert_eq!(cell.count(), 2);
        assert_eq!(cell.slot_of_edge[0], 0);
        assert_eq!(cell.slot_of_edge[2], 0);
        assert_eq!(cell.slot_of_edge[4], 1);
        assert_eq!(cell.slot_of_edge[6], 1);
        assert_eq!(cell.slot_of_edge[7], NO_SLOT);
    }
}
