//! Which chunks an edit invalidated, and how much of the marked region actually
//! changed.
//!
//! # The number this exists to produce
//!
//! The research doc's item E1 asks for **the fraction of cells that actually
//! change per edit**, and records that nobody has published it. It is the ceiling
//! on every incremental idea in the opportunities doc: if a brush stroke really
//! changes 3% of the cells it touches, then an incremental scheme has 97% to win
//! back and is worth building. If it changes 90%, re-meshing the region wholesale
//! is the right answer and the complexity buys nothing.
//!
//! So this module reports two numbers rather than one, because they answer
//! different questions:
//!
//! - **Value-changed cells** — any of the eight corner samples differs at all.
//!   This is the set that genuinely needs re-meshing, since a vertex position
//!   interpolates the values and moves when they move.
//! - **Sign-changed cells** — the *pattern* of inside/outside corners differs.
//!   This is where the topology moved, and it is a subset of the above.
//!
//! The gap between them is the interesting part: cells that shifted a vertex
//! without changing which triangles exist.
//!
//! # What is not decided here
//!
//! Nothing owns sample data, for the reason [`ChunkLayout`]
//! does not: an edit is expressed as *two fields*, before and after, which is
//! what the crate's field-first design already gives. A consumer holding a
//! fixed-point slab and a consumer holding an expression tree both produce a
//! before and an after, and neither has to adopt a storage type from this crate
//! to be measured.

use alloc::vec::Vec;

use super::{ChunkId, ChunkLayout};
use crate::cube::{corner_offset, is_inside};
use crate::{Real, Sdf};

/// Chunks awaiting a re-mesh.
///
/// Sorted and deduplicated on insert, so iteration order is a pure function of
/// the set's contents rather than of the order edits arrived — the same property
/// [`validate`](crate::validate) is built around, and for the same reason: a
/// re-mesh queue that depends on insertion order gives different output for the
/// same world.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct DirtySet {
    chunks: Vec<ChunkId>,
}

impl DirtySet {
    /// An empty set.
    #[must_use]
    pub const fn new() -> Self {
        Self { chunks: Vec::new() }
    }

    /// Mark one chunk.
    pub fn insert(&mut self, id: ChunkId) {
        if let Err(at) = self.chunks.binary_search(&id) {
            self.chunks.insert(at, id);
        }
    }

    /// How many chunks are waiting.
    #[must_use]
    pub fn len(&self) -> usize {
        self.chunks.len()
    }

    /// Whether anything is waiting.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.chunks.is_empty()
    }

    /// The chunks, in ascending order.
    pub fn iter(&self) -> impl Iterator<Item = ChunkId> + '_ {
        self.chunks.iter().copied()
    }

    /// Forget everything.
    pub fn clear(&mut self) {
        self.chunks.clear();
    }

    /// Re-mesh every dirty chunk, then clear the set.
    ///
    /// `mesh` is called once per chunk with that chunk's id and the world origin
    /// to extract at — which comes from
    /// [`ChunkLayout::sample_origin`](super::ChunkLayout::sample_origin), so the
    /// caller never computes it and cannot get the seam arithmetic wrong (M-32).
    ///
    /// Returns how many chunks were re-meshed. The set is cleared **after** the
    /// last call, so a panic mid-way leaves the queue intact rather than losing
    /// the work.
    pub fn mesh_dirty<R, F>(&mut self, layout: &ChunkLayout<R>, mut mesh: F) -> usize
    where
        R: Real,
        F: FnMut(ChunkId, [R; 3]),
    {
        for id in &self.chunks {
            mesh(*id, layout.sample_origin(*id));
        }
        let done = self.chunks.len();
        self.chunks.clear();
        done
    }
}

/// What an edit actually touched.
///
/// Counts are over **cells**, not samples: a cell is the unit of meshing work,
/// and it is what an incremental scheme would avoid.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct EditReport {
    /// Cells inside the marked region — the work a naive "re-mesh the brush's
    /// bounding box" scheme would do.
    pub region_cells: u64,
    /// Cells where any of the eight corner samples changed at all.
    ///
    /// **Not** the re-meshing set, and the difference matters. An SDF brush
    /// changes values throughout its support — a sphere carved deeper moves the
    /// value at every point the sphere term dominates, including deep inside the
    /// solid. Those cells emit no triangles either way, so their *output* is
    /// unchanged however much their values moved.
    pub value_changed_cells: u64,
    /// Cells whose **triangles** change: a surface cell before or after, whose
    /// values moved. **This is the set that genuinely needs re-meshing**, and the
    /// numerator of [`changed_fraction`](Self::changed_fraction).
    ///
    /// A cell with no sign change at all emits nothing before and nothing after,
    /// so re-meshing it produces the same empty result.
    pub output_changed_cells: u64,
    /// Cells where the inside/outside pattern changed.
    ///
    /// **Not** a subset of [`output_changed_cells`](Self::output_changed_cells),
    /// which is worth stating because the containment looks obvious and is false:
    /// a cell can go from *wholly outside* to *wholly inside* in one edit, and
    /// then every corner flipped while the cell emitted nothing before and
    /// nothing after. See [`swept_cells`](Self::swept_cells).
    pub sign_changed_cells: u64,
    /// Cells the surface passed **entirely through** in one edit: every corner
    /// flipped, and the cell was a surface cell neither before nor after.
    ///
    /// These need no re-mesh — nothing before, nothing after — but they are the
    /// signature of an edit that moved further than one cell per step, which is
    /// the regime where a brush stroke starts skipping geometry between frames.
    /// Worth watching rather than acting on.
    pub swept_cells: u64,
    /// Chunks the region overlaps.
    pub region_chunks: u64,
    /// Chunks containing at least one output-changed cell — what actually goes on
    /// the dirty queue.
    pub dirty_chunks: u64,
}

impl EditReport {
    /// **E1**: the fraction of the marked region whose mesh output changes.
    ///
    /// Zero when the region is empty. This is the number the research asks for
    /// and reports as unpublished, and it is over
    /// [`output_changed_cells`](Self::output_changed_cells) rather than
    /// [`value_changed_cells`](Self::value_changed_cells) — a cell that emits no
    /// triangles before or after has not changed, whatever happened to its
    /// samples.
    #[must_use]
    pub fn changed_fraction(&self) -> f64 {
        if self.region_cells == 0 {
            0.0
        } else {
            self.output_changed_cells as f64 / self.region_cells as f64
        }
    }

    /// The fraction whose samples moved at all, output or not.
    ///
    /// Larger than [`changed_fraction`](Self::changed_fraction), and the gap is
    /// the volume an SDF brush disturbs without changing any geometry.
    #[must_use]
    pub fn value_changed_fraction(&self) -> f64 {
        if self.region_cells == 0 {
            0.0
        } else {
            self.value_changed_cells as f64 / self.region_cells as f64
        }
    }

    /// The fraction of *chunks* in the region that need re-meshing.
    ///
    /// Coarser than [`changed_fraction`](Self::changed_fraction) and usually much
    /// larger, because one changed cell dirties its whole chunk. The gap between
    /// the two is what a finer dirty granularity would buy.
    #[must_use]
    pub fn dirty_chunk_fraction(&self) -> f64 {
        if self.region_chunks == 0 {
            0.0
        } else {
            self.dirty_chunks as f64 / self.region_chunks as f64
        }
    }
}

/// Compare two fields over a region, mark the chunks that changed, and report
/// how much of the region actually moved.
///
/// `before` and `after` are the field either side of one edit. The region is an
/// inclusive range of **global cell indices**, which is what a brush's bounding
/// box converts to — cells, not samples, because a cell is the unit of work.
///
/// Costs `(region + 1)³` field evaluations of each field: it samples the corner
/// grid once and reuses each sample across the eight cells that share it, rather
/// than evaluating eight times per cell.
pub fn mark_edit<R, A, B>(
    layout: &ChunkLayout<R>,
    before: &A,
    after: &B,
    min_cell: [i64; 3],
    max_cell: [i64; 3],
    dirty: &mut DirtySet,
) -> EditReport
where
    R: Real,
    A: Sdf<Scalar = R>,
    B: Sdf<Scalar = R>,
{
    let mut report = EditReport::default();
    for axis in 0..3 {
        if max_cell[axis] < min_cell[axis] {
            return report;
        }
    }

    // Sample the corner grid spanning the region: one plane wider than the cells
    // on each axis, since a cell's corners are its own index and the next.
    let extent = [
        (max_cell[0] - min_cell[0] + 2) as usize,
        (max_cell[1] - min_cell[1] + 2) as usize,
        (max_cell[2] - min_cell[2] + 2) as usize,
    ];
    let count = extent[0] * extent[1] * extent[2];
    let mut changed = alloc::vec![false; count];
    let mut inside_before = alloc::vec![false; count];
    let mut inside_after = alloc::vec![false; count];

    let index = |x: usize, y: usize, z: usize| x + extent[0] * (y + extent[1] * z);
    for z in 0..extent[2] {
        for y in 0..extent[1] {
            for x in 0..extent[0] {
                let sample = [
                    min_cell[0] + x as i64,
                    min_cell[1] + y as i64,
                    min_cell[2] + z as i64,
                ];
                let p = layout.world_of_sample(sample);
                let a = before.sample(p);
                let b = after.sample(p);
                let i = index(x, y, z);
                // Bit comparison, not `!=`: `+0.0 == -0.0` is true and a sign
                // flip on a zero sample is a real change to the classification
                // this crate makes, since zero is outside.
                changed[i] = a.total_cmp(&b) != core::cmp::Ordering::Equal;
                inside_before[i] = is_inside(a);
                inside_after[i] = is_inside(b);
            }
        }
    }

    let mut region_chunks = DirtySet::new();
    for cz in 0..extent[2] - 1 {
        for cy in 0..extent[1] - 1 {
            for cx in 0..extent[0] - 1 {
                report.region_cells += 1;

                let mut any_value = false;
                let mut any_sign = false;
                let mut inside_count_before = 0u32;
                let mut inside_count_after = 0u32;
                for corner in 0..8u8 {
                    let o = corner_offset(corner);
                    let i = index(cx + o[0] as usize, cy + o[1] as usize, cz + o[2] as usize);
                    any_value |= changed[i];
                    any_sign |= inside_before[i] != inside_after[i];
                    inside_count_before += u32::from(inside_before[i]);
                    inside_count_after += u32::from(inside_after[i]);
                }
                // A cell emits triangles only when its corners disagree. One that
                // is wholly inside or wholly outside emits nothing, both before
                // and after, however far its samples moved.
                let was_surface = inside_count_before != 0 && inside_count_before != 8;
                let is_surface = inside_count_after != 0 && inside_count_after != 8;

                let cell = [
                    min_cell[0] + cx as i64,
                    min_cell[1] + cy as i64,
                    min_cell[2] + cz as i64,
                ];
                // A cell is meshed by the chunk that owns it, and cell `k` is
                // owned by the chunk owning sample `k` — the overlap plane sits
                // above the last cell, not on it.
                let (owner, _) = layout.local_sample(cell);
                region_chunks.insert(owner);

                if any_value {
                    report.value_changed_cells += 1;
                }
                if any_value && (was_surface || is_surface) {
                    report.output_changed_cells += 1;
                    dirty.insert(owner);
                }
                if any_sign {
                    report.sign_changed_cells += 1;
                    if !was_surface && !is_surface {
                        report.swept_cells += 1;
                    }
                }
            }
        }
    }

    report.region_chunks = region_chunks.len() as u64;
    report.dirty_chunks = region_chunks
        .iter()
        .filter(|id| dirty.iter().any(|d| d == *id))
        .count() as u64;
    report
}

#[cfg(test)]
mod tests;
