//! Chunking: the coordinate system a streamed world is meshed in.
//!
//! The real workload is not "mesh one volume". It is "re-mesh the 3% of cells a
//! brush touched", and that means cutting space into chunks that can be meshed
//! independently and must nonetheless agree along their shared faces.
//!
//! # The overlap, and why it is on the positive faces
//!
//! A chunk **owns** `cells³` cells and **samples** `cells + 1` planes per axis.
//! The extra plane is the overlap: it is the same physical plane as the *next*
//! chunk's first one. Cells are what a chunk owns; samples are what it borrows.
//!
//! ```text
//! chunk 0 samples:  0 1 2 3 4          (owns cells 0..4)
//! chunk 1 samples:          4 5 6 7 8  (owns cells 4..8)
//!                           ^
//!                           shared plane
//! ```
//!
//! Putting the overlap on the positive faces rather than both means each sample
//! plane has exactly one owner, so there is no question of which chunk is
//! authoritative and no double-counting when cells are marked dirty.
//!
//! # What breaks here, and it is supposed to
//!
//! ✗1's identity — `F_sn = F_mc + 2χ` — has a recorded break condition:
//! boundary-clipped meshes. **A single chunk's mesh is boundary-clipped by
//! construction**, so it is a manifold *with boundary* and the identity does not
//! apply to it. That was written down before this module existed; it is not a
//! regression and must not be "fixed".
//!
//! The corresponding gate change is that a per-chunk mesh is held to
//! `is_manifold`, never `is_closed`.

use crate::{Real, RuntimeShape3};

/// Integer coordinates of a chunk on the chunk lattice.
///
/// Signed, because a streamed world extends in both directions from wherever the
/// origin happens to be and a chunk index is not a natural number.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct ChunkId {
    /// Chunk coordinates. Chunk `[0, 0, 0]` owns the cells whose minimum corner
    /// is the layout's origin.
    pub coords: [i32; 3],
}

impl ChunkId {
    /// A chunk at these coordinates.
    #[must_use]
    pub const fn new(coords: [i32; 3]) -> Self {
        Self { coords }
    }

    /// The neighbour one step along `axis`.
    #[must_use]
    pub fn neighbour(self, axis: usize, step: i32) -> Self {
        let mut coords = self.coords;
        coords[axis] += step;
        Self { coords }
    }
}

/// The layout every chunk in a volume shares.
///
/// Private fields and one checked constructor, for the same reason
/// [`ValidateConfig`](crate::validate::ValidateConfig) has them: a layout that
/// exists is one whose arithmetic is meaningful.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ChunkLayout<R: Real> {
    cells: u32,
    cell_size: R,
    origin: [R; 3],
}

impl<R: Real> ChunkLayout<R> {
    /// A layout of `cells³`-cell chunks.
    ///
    /// `origin` is the world position of global sample `[0, 0, 0]`, which is the
    /// minimum corner of chunk `[0, 0, 0]`.
    ///
    /// # Errors
    ///
    /// [`Error::GridTooSmall`](crate::Error::GridTooSmall) if `cells` is zero —
    /// a chunk owning no cells has nothing to mesh.
    ///
    /// [`Error::InvalidCellSize`](crate::Error::InvalidCellSize) if `cell_size`
    /// is not finite and positive.
    pub fn new(cells: u32, cell_size: R, origin: [R; 3]) -> crate::Result<Self> {
        if cells == 0 {
            return Err(crate::Error::GridTooSmall {
                size: [cells + 1; 3],
            });
        }
        if !cell_size.is_finite() || cell_size <= R::ZERO {
            return Err(crate::Error::InvalidCellSize {
                value: f64::from(cell_size.as_f32()),
            });
        }
        Ok(Self {
            cells,
            cell_size,
            origin,
        })
    }

    /// Cells owned per axis.
    #[must_use]
    pub fn cells(&self) -> u32 {
        self.cells
    }

    /// Grid spacing.
    #[must_use]
    pub fn cell_size(&self) -> R {
        self.cell_size
    }

    /// The sample grid one chunk is meshed on: `cells + 1` per axis.
    ///
    /// The `+ 1` is the positive-face overlap. Pass this straight to an
    /// extractor — it counts **samples**, which is the same convention every
    /// extractor uses.
    ///
    /// # Errors
    ///
    /// [`Error::ShapeOverflow`](crate::Error::ShapeOverflow) if `cells + 1`
    /// cubed does not fit the index space.
    pub fn sample_shape(&self) -> crate::Result<RuntimeShape3> {
        RuntimeShape3::new([self.cells + 1; 3])
    }

    /// Global sample index of a chunk's minimum corner.
    ///
    /// Exact integer arithmetic, deliberately: this is the quantity every world
    /// position is derived from, and deriving it in floating point is how two
    /// chunks end up disagreeing about where their shared plane is.
    #[must_use]
    pub fn base_sample(&self, id: ChunkId) -> [i64; 3] {
        let n = i64::from(self.cells);
        [
            i64::from(id.coords[0]) * n,
            i64::from(id.coords[1]) * n,
            i64::from(id.coords[2]) * n,
        ]
    }

    /// World position of a **global** sample index.
    ///
    /// The single place a sample's world position is defined. Everything else
    /// routes through it so that a sample shared by two chunks is computed by one
    /// expression rather than two — see
    /// [`sample_origin`](Self::sample_origin) for why that matters.
    #[must_use]
    pub fn world_of_sample(&self, sample: [i64; 3]) -> [R; 3] {
        [
            self.origin[0] + self.cell_size * R::from_f64(sample[0] as f64),
            self.origin[1] + self.cell_size * R::from_f64(sample[1] as f64),
            self.origin[2] + self.cell_size * R::from_f64(sample[2] as f64),
        ]
    }

    /// World position to pass an extractor as its `origin` for this chunk.
    ///
    /// # The one place chunk seams can crack, and it is arithmetic
    ///
    /// An extractor computes its samples as `origin + cell_size · local`. Chunk
    /// `c`'s **last** plane is therefore `(o + h·cn) + h·n`, and chunk `c+1`'s
    /// **first** plane is `o + h·(c+1)n` — the same point, reached by two
    /// different expressions. Floating-point addition is not associative, so
    /// those are **not** bit-equal in general, and two chunks then sample the
    /// same physical plane at very slightly different places, get very slightly
    /// different field values, and interpolate crossings that do not quite line
    /// up. That is a hairline crack, and it is invisible until someone walks a
    /// character over it.
    ///
    /// **They are exact when `cell_size` is a power of two**, because then
    /// `h · k` is exact for every integer `k` in range and both expressions are
    /// exact. G-001 measured the difference for both cases; see M-32.
    ///
    /// This function returns the chunk origin as
    /// `world_of_sample(base_sample(id))`, which is the best available: one
    /// expression, evaluated once, from exact integers.
    #[must_use]
    pub fn sample_origin(&self, id: ChunkId) -> [R; 3] {
        self.world_of_sample(self.base_sample(id))
    }

    /// Which chunk owns the cell containing `point`.
    ///
    /// Ownership is half-open: a point exactly on a shared plane belongs to the
    /// chunk on its **positive** side, matching where the overlap sits.
    #[must_use]
    pub fn chunk_of(&self, point: [R; 3]) -> ChunkId {
        let inv = self.cell_size.recip();
        let n = R::from_f64(f64::from(self.cells));
        let mut coords = [0i32; 3];
        for (axis, slot) in coords.iter_mut().enumerate() {
            let cell = ((point[axis] - self.origin[axis]) * inv / n).floor();
            *slot = if cell.is_finite() {
                cell.as_f32() as i32
            } else {
                0
            };
        }
        ChunkId { coords }
    }

    /// Global sample index of a chunk-local sample.
    #[must_use]
    pub fn global_sample(&self, id: ChunkId, local: [u32; 3]) -> [i64; 3] {
        let base = self.base_sample(id);
        [
            base[0] + i64::from(local[0]),
            base[1] + i64::from(local[1]),
            base[2] + i64::from(local[2]),
        ]
    }

    /// The chunk and local sample a global sample index belongs to.
    ///
    /// The inverse of [`global_sample`](Self::global_sample) on the samples a
    /// chunk **owns**, which is `0..cells` — the overlap plane at `cells` maps to
    /// the next chunk's `0`, because that is where it is owned.
    #[must_use]
    pub fn local_sample(&self, sample: [i64; 3]) -> (ChunkId, [u32; 3]) {
        let n = i64::from(self.cells);
        let mut coords = [0i32; 3];
        let mut local = [0u32; 3];
        for axis in 0..3 {
            // Floor division, so negative coordinates land in the chunk below
            // rather than being truncated toward zero.
            let c = sample[axis].div_euclid(n);
            let l = sample[axis].rem_euclid(n);
            coords[axis] = c as i32;
            local[axis] = l as u32;
        }
        (ChunkId { coords }, local)
    }
}

#[cfg(test)]
mod tests;
