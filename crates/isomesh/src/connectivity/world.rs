//! Several [`Air`] grids stitched into one world, so `connected` answers across
//! a chunk seam.
//!
//! Ticket: R-028, split out of R-022b on M-321.
//!
//! # Why this exists, with the number that made it a ticket
//!
//! [`Air::fill`] repairs connectivity by lockstep search, which costs the
//! **second-largest** piece of a split. On the measured edit distribution that is
//! 3.38 voxels visited per seed (M-321). On a **bisect** — filling one voxel at
//! the midpoint of a tunnel joining two equal caverns — it is 123,039 per seed,
//! and the whole repair costs **1.1× a full rebuild**.
//!
//! **HDT is not the remedy, and M-321 says why.** Its levels bound a
//! *replacement-edge search*; on a genuine bisect there is no replacement edge to
//! find, because the component really did split. The cost is the unavoidable one
//! of discovering which side is smaller.
//!
//! **Decomposition is the remedy.** A search inside one `Air` cannot cost more
//! than that `Air`, so a per-chunk grid bounds a bisect at the chunk rather than
//! at the world. That is the whole idea, and it is why this type holds many small
//! `Air`s rather than one large one.
//!
//! # Why a consumer needed this to use `fill` at all
//!
//! Without it, `fill` forces a choice where neither option is complete:
//!
//! | | cost |
//! |---|---|
//! | one large `Air` for the world | the bisect tail, 1.1× a rebuild |
//! | one `Air` per chunk | no cross-chunk `connected` at all |
//!
//! Every extractor here is driven per chunk, so the second is the natural shape
//! and it silently cannot answer the question the module exists for. This type is
//! the missing half rather than an optimisation.
//!
//! # How the seam is exact rather than approximate
//!
//! [`ChunkLayout::sample_shape`] is `cells + 1` and
//! [`ChunkLayout::base_sample`] is `coords · cells`, so adjacent chunks **share
//! exactly one sample plane**: chunk `c`'s local `cells` plane and chunk `c+1`'s
//! local `0` plane are the same global samples. Two labels are joined where that
//! shared sample is air on both sides. Nothing is interpolated and nothing is
//! matched by position tolerance.
//!
//! [`ChunkLayout::local_sample`] already resolves the ownership question — the
//! overlap plane maps to the **next** chunk's local `0` — so a global sample has
//! exactly one owner and a query never has to pick between two answers.
//!
//! # Why the global graph is rebuilt rather than maintained
//!
//! It is a union-find, and ✗26 is precisely the reason a union-find must not be
//! asked to delete. It is never asked: **every restitch builds it from scratch**,
//! so it only ever unions. That is sound where the incremental version was not,
//! and it is affordable because the graph is tiny — its nodes are *components*,
//! not samples.
//!
//! What is *not* rebuilt from scratch is the expensive part: the per-seam label
//! pairs, which need an `O(cells²)` scan of a shared plane. Those are cached per
//! seam and recomputed only for seams touching a chunk that changed.

use alloc::vec::Vec;

use super::{Air, Fill, Repair};
use crate::Real;
use crate::chunk::{ChunkId, ChunkLayout};

/// No component.
const NONE: u32 = u32::MAX;

/// One loaded chunk.
#[derive(Debug)]
struct Loaded {
    id: ChunkId,
    air: Air,
    /// World component id per local label. Rebuilt by every restitch.
    global: Vec<u32>,
}

/// The distinct label pairs joined across one shared plane.
///
/// Compressed deliberately: a shared plane has `(cells + 1)²` samples and
/// typically a handful of distinct component pairs across it, so caching the
/// pairs rather than the plane is what keeps a restitch proportional to the
/// number of components rather than to the seam's area.
#[derive(Debug)]
struct Seam {
    lo: ChunkId,
    hi: ChunkId,
    pairs: Vec<(u32, u32)>,
}

/// What the last restitch cost.
///
/// Counts, not durations — the same standing [`Repair`] and [`Fill`] have.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Default)]
pub struct Seams {
    /// Seams whose label pairs were rescanned, at `O(cells²)` each.
    pub rescanned: u64,
    /// Shared-plane samples visited by those rescans.
    pub plane_samples: u64,
    /// Nodes in the global graph — components, not samples.
    pub nodes: u64,
    /// Label pairs unioned across all cached seams.
    pub pairs: u64,
}

/// Many [`Air`] grids over one [`ChunkLayout`], answering `connected` across
/// seams.
///
/// Coordinates in this type's API are **global sample indices**, the same space
/// [`ChunkLayout::global_sample`] produces. Chunk-local coordinates stay inside
/// the per-chunk [`Air`], which is reachable through [`chunk`](Self::chunk).
#[derive(Debug)]
pub struct AirWorld<R: Real> {
    layout: ChunkLayout<R>,
    /// Sorted by [`ChunkId`], so every traversal is deterministic (M-36).
    chunks: Vec<Loaded>,
    seams: Vec<Seam>,
    components: u64,
    last: Seams,
}

impl<R: Real> AirWorld<R> {
    /// An empty world over this layout.
    #[must_use]
    pub fn new(layout: ChunkLayout<R>) -> Self {
        Self {
            layout,
            chunks: Vec::new(),
            seams: Vec::new(),
            components: 0,
            last: Seams::default(),
        }
    }

    /// Add or replace a chunk's samples, air being `value >= 0`.
    ///
    /// `values` is one entry per sample of [`ChunkLayout::sample_shape`], which
    /// is `(cells + 1)³` — the plane a chunk shares with each neighbour is
    /// included, because that shared plane is what the stitch reads.
    ///
    /// # Errors
    ///
    /// [`Error::ShapeOverflow`](crate::Error::ShapeOverflow) if the layout's
    /// sample shape does not fit, or if `values` is the wrong length.
    pub fn load(&mut self, id: ChunkId, values: &[R]) -> crate::Result<Repair> {
        let shape = self.layout.sample_shape()?;
        let (air, repair) = Air::build(values, &shape)?;
        let global = alloc::vec![NONE; air.label_count()];
        let entry = Loaded { id, air, global };
        match self.chunks.binary_search_by(|c| c.id.cmp(&id)) {
            Ok(at) => {
                if let Some(slot) = self.chunks.get_mut(at) {
                    *slot = entry;
                }
            }
            Err(at) => self.chunks.insert(at, entry),
        }
        self.rescan_around(id);
        self.restitch();
        Ok(repair)
    }

    /// Dig in one chunk, then restitch.
    ///
    /// `None` if that chunk is not loaded. Coordinates are **chunk-local**,
    /// matching [`Air::dig`].
    pub fn dig<B: FnMut() -> bool>(
        &mut self,
        id: ChunkId,
        samples: &[[u32; 3]],
        spend: B,
    ) -> Option<Repair> {
        let at = self.chunks.binary_search_by(|c| c.id.cmp(&id)).ok()?;
        let repair = self.chunks.get_mut(at)?.air.dig(samples, spend);
        self.resize_global(at);
        self.rescan_around(id);
        self.restitch();
        Some(repair)
    }

    /// Fill in one chunk, then restitch.
    ///
    /// `None` if that chunk is not loaded. Coordinates are **chunk-local**,
    /// matching [`Air::fill`].
    ///
    /// **This is the operation R-028 exists to bound.** The lockstep search runs
    /// inside one `Air`, so it cannot visit more than that chunk however large
    /// the world is.
    pub fn fill<B: FnMut() -> bool>(
        &mut self,
        id: ChunkId,
        samples: &[[u32; 3]],
        spend: B,
    ) -> Option<Fill> {
        let at = self.chunks.binary_search_by(|c| c.id.cmp(&id)).ok()?;
        let out = self.chunks.get_mut(at)?.air.fill(samples, spend);
        self.resize_global(at);
        self.rescan_around(id);
        self.restitch();
        Some(out)
    }

    /// Are these two **global** samples in the same air component?
    ///
    /// `false` if either is solid, or in a chunk that is not loaded — an
    /// unloaded chunk is not evidence of connection.
    #[must_use]
    pub fn connected(&self, a: [i64; 3], b: [i64; 3]) -> bool {
        match (self.component_at(a), self.component_at(b)) {
            (Some(x), Some(y)) => x == y,
            _ => false,
        }
    }

    /// The world component id of a **global** sample, or `None` if it is solid
    /// or unloaded.
    ///
    /// # Which samples a world can name
    ///
    /// Chunks `0..k` along an axis answer for global samples `0..k·cells`. **The
    /// plane at `k·cells` is not among them**, because
    /// [`ChunkLayout::local_sample`] gives an overlap plane to the *next* chunk
    /// and that chunk is not loaded. The loaded chunks' arrays do contain those
    /// samples and do label them — the stitch reads them — they simply have no
    /// name in this coordinate space.
    ///
    /// So a consumer covering a finite region loads one chunk beyond it, exactly
    /// as it would to mesh that region's far face. This is the ownership rule
    /// `local_sample` already documents, surfacing where it is felt.
    #[must_use]
    pub fn component_at(&self, sample: [i64; 3]) -> Option<u32> {
        let (id, local) = self.layout.local_sample(sample);
        let at = self.chunks.binary_search_by(|c| c.id.cmp(&id)).ok()?;
        let entry = self.chunks.get(at)?;
        let label = entry.air.label_of(local)?;
        match entry.global.get(label as usize) {
            Some(&g) if g != NONE => Some(g),
            _ => None,
        }
    }

    /// Air components across the whole world, seams accounted for.
    #[must_use]
    pub fn components(&self) -> u64 {
        self.components
    }

    /// One chunk's grid, for the per-chunk questions this type does not answer.
    #[must_use]
    pub fn chunk(&self, id: ChunkId) -> Option<&Air> {
        let at = self.chunks.binary_search_by(|c| c.id.cmp(&id)).ok()?;
        self.chunks.get(at).map(|c| &c.air)
    }

    /// Loaded chunks.
    #[must_use]
    pub fn loaded(&self) -> usize {
        self.chunks.len()
    }

    /// What the last restitch cost.
    #[must_use]
    pub fn last_seams(&self) -> Seams {
        self.last
    }

    // --- internals ---------------------------------------------------------

    /// Grow a chunk's per-label table after an edit issued new labels.
    fn resize_global(&mut self, at: usize) {
        let Some(entry) = self.chunks.get_mut(at) else {
            return;
        };
        let want = entry.air.label_count();
        if entry.global.len() < want {
            entry.global.resize(want, NONE);
        }
    }

    /// Recompute the cached label pairs for every seam touching `id`.
    ///
    /// This is the `O(cells²)` half, and it is confined to the six seams of one
    /// chunk rather than to the world.
    fn rescan_around(&mut self, id: ChunkId) {
        let mut touched: Vec<(ChunkId, ChunkId)> = Vec::new();
        for axis in 0..3 {
            for step in [-1i32, 1] {
                let Some(other) = neighbour(id, axis, step) else {
                    continue;
                };
                if self.chunks.binary_search_by(|c| c.id.cmp(&other)).is_err() {
                    continue;
                }
                let pair = if step > 0 { (id, other) } else { (other, id) };
                touched.push(pair);
            }
        }

        for (lo, hi) in touched {
            let pairs = self.scan_seam(lo, hi);
            match self.seams.iter_mut().find(|s| s.lo == lo && s.hi == hi) {
                Some(seam) => seam.pairs = pairs,
                None => self.seams.push(Seam { lo, hi, pairs }),
            }
        }
        // A seam whose chunk was replaced can be left stale otherwise.
        self.seams.retain(|s| {
            self.chunks.binary_search_by(|c| c.id.cmp(&s.lo)).is_ok()
                && self.chunks.binary_search_by(|c| c.id.cmp(&s.hi)).is_ok()
        });
    }

    /// The distinct label pairs across the plane `lo` and `hi` share.
    ///
    /// `hi` is `lo`'s neighbour one step along exactly one axis; `lo`'s local
    /// `cells` plane and `hi`'s local `0` plane are the same global samples.
    fn scan_seam(&mut self, lo: ChunkId, hi: ChunkId) -> Vec<(u32, u32)> {
        let mut pairs: Vec<(u32, u32)> = Vec::new();
        let axis = match (0..3).find(|&a| lo.coords[a] != hi.coords[a]) {
            Some(a) => a,
            None => return pairs,
        };
        let (Ok(a_at), Ok(b_at)) = (
            self.chunks.binary_search_by(|c| c.id.cmp(&lo)),
            self.chunks.binary_search_by(|c| c.id.cmp(&hi)),
        ) else {
            return pairs;
        };
        let cells = self.layout.cells();
        let (u, v) = match axis {
            0 => (1, 2),
            1 => (0, 2),
            _ => (0, 1),
        };

        let mut samples = 0u64;
        for a in 0..=cells {
            for b in 0..=cells {
                let mut lo_local = [0u32; 3];
                lo_local[axis] = cells;
                lo_local[u] = a;
                lo_local[v] = b;
                let mut hi_local = [0u32; 3];
                hi_local[axis] = 0;
                hi_local[u] = a;
                hi_local[v] = b;

                samples += 1;
                let la = self.chunks.get(a_at).and_then(|c| c.air.label_of(lo_local));
                let lb = self.chunks.get(b_at).and_then(|c| c.air.label_of(hi_local));
                if let (Some(la), Some(lb)) = (la, lb)
                    && !pairs.contains(&(la, lb))
                {
                    pairs.push((la, lb));
                }
            }
        }
        self.last.rescanned += 1;
        self.last.plane_samples += samples;
        pairs
    }

    /// Rebuild the global component graph from the cached seams.
    ///
    /// **From scratch, every time.** Its nodes are components rather than
    /// samples, so it is small — and building it fresh means it only ever
    /// unions, which is the one thing a union-find is safe to do (✗26).
    fn restitch(&mut self) {
        let (rescanned, plane_samples) = (self.last.rescanned, self.last.plane_samples);
        self.last = Seams {
            rescanned,
            plane_samples,
            ..Seams::default()
        };

        // Node numbering: chunk `k`'s label `l` is `offset[k] + l`.
        let mut offset = Vec::with_capacity(self.chunks.len());
        let mut total = 0usize;
        for c in &self.chunks {
            offset.push(total);
            total += c.air.label_count();
        }
        self.last.nodes = total as u64;

        let mut parent: Vec<u32> = (0..total as u32).collect();
        fn find(parent: &mut [u32], mut i: u32) -> u32 {
            while parent.get(i as usize).copied().unwrap_or(i) != i {
                let p = parent.get(i as usize).copied().unwrap_or(i);
                let grand = parent.get(p as usize).copied().unwrap_or(p);
                if let Some(slot) = parent.get_mut(i as usize) {
                    *slot = grand;
                }
                i = grand;
            }
            i
        }

        let mut unions = 0u64;
        for seam in &self.seams {
            let (Ok(a_at), Ok(b_at)) = (
                self.chunks.binary_search_by(|c| c.id.cmp(&seam.lo)),
                self.chunks.binary_search_by(|c| c.id.cmp(&seam.hi)),
            ) else {
                continue;
            };
            let (Some(&ao), Some(&bo)) = (offset.get(a_at), offset.get(b_at)) else {
                continue;
            };
            for &(la, lb) in &seam.pairs {
                let x = find(&mut parent, ao as u32 + la);
                let y = find(&mut parent, bo as u32 + lb);
                unions += 1;
                if x != y
                    && let Some(slot) = parent.get_mut(y as usize)
                {
                    *slot = x;
                }
            }
        }
        self.last.pairs = unions;

        // Compact the roots of live labels into dense world ids.
        let mut world_of: Vec<u32> = alloc::vec![NONE; total];
        let mut next = 0u32;
        for (k, c) in self.chunks.iter().enumerate() {
            let Some(&off) = offset.get(k) else { continue };
            for l in 0..c.air.label_count() {
                if c.air.component_size(l as u32) == 0 {
                    continue;
                }
                let root = find(&mut parent, (off + l) as u32);
                if world_of.get(root as usize) == Some(&NONE) {
                    if let Some(slot) = world_of.get_mut(root as usize) {
                        *slot = next;
                    }
                    next += 1;
                }
            }
        }
        self.components = u64::from(next);

        for (k, c) in self.chunks.iter_mut().enumerate() {
            let Some(&off) = offset.get(k) else { continue };
            c.global.clear();
            c.global.resize(c.air.label_count(), NONE);
            for l in 0..c.air.label_count() {
                if c.air.component_size(l as u32) == 0 {
                    continue;
                }
                let root = find(&mut parent, (off + l) as u32);
                let g = world_of.get(root as usize).copied().unwrap_or(NONE);
                if let Some(slot) = c.global.get_mut(l) {
                    *slot = g;
                }
            }
        }
    }
}

/// `id`'s neighbour along `axis`, or `None` if it leaves the chunk lattice.
///
/// [`ChunkId::neighbour`] panics on overflow, and this crate does not panic
/// (CLAUDE.md), so the check is here rather than there.
fn neighbour(id: ChunkId, axis: usize, step: i32) -> Option<ChunkId> {
    let mut coords = id.coords;
    let slot = coords.get_mut(axis)?;
    *slot = slot.checked_add(step)?;
    Some(ChunkId::new(coords))
}

#[cfg(test)]
mod tests;
