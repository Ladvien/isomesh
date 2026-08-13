//! G-007 — which chunks are resident, as the camera moves.
//!
//! [`ChunkLayout`] says where chunks are and [`DirtySet`](super::dirty::DirtySet)
//! says which need re-meshing. This says which should **exist**: a residency set
//! that follows a camera, and the difference between one frame's set and the
//! next as a pair of lists the caller acts on.
//!
//! # Hysteresis is the whole point
//!
//! A single radius thrashes. A chunk whose distance sits within a hair of the
//! threshold loads and unloads on alternate frames as the camera jitters, and
//! each cycle costs a full re-mesh — the most expensive thing a streaming world
//! does. So there are two radii and a chunk's fate depends on what it already
//! is: outside the resident set it must come *within* [`load`](StreamConfig::load)
//! to enter, and once inside it must go *beyond*
//! [`unload`](StreamConfig::unload) to leave. Between the two it keeps whatever
//! state it has.
//!
//! `unload > load` is enforced in the constructor rather than documented as a
//! precondition, because equal radii is exactly the degenerate case the feature
//! exists to prevent, and a config that cannot express it cannot ship it.
//!
//! # Distance is to the chunk, not to its centre
//!
//! [`mesh_within_budget`](super::dirty::DirtySet::mesh_within_budget) ranks by
//! centre distance, and that is right for an *ordering* — any consistent ranking
//! serves. Residency is a **threshold**, and a threshold has a meaning: a chunk
//! with geometry inside the load radius should be resident, whatever its centre
//! is doing. A chunk is `2·cells·h` across the diagonal, so centre distance gets
//! that wrong by up to half a chunk, and at the small radii a test or a tight
//! memory budget uses, half a chunk is most of the radius.
//!
//! So the measure is the distance from the camera to the chunk's axis-aligned
//! box, which is zero inside it and exact outside.

use alloc::vec::Vec;

use crate::real::Real;

use super::{ChunkId, ChunkLayout};

/// The two radii a residency decision is made against.
///
/// Private fields and one checked constructor, the same shape as
/// [`ValidateConfig`](crate::validate::ValidateConfig): the invalid state is not
/// reported, it is unrepresentable.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct StreamConfig<R: Real> {
    load: R,
    unload: R,
}

impl<R: Real> StreamConfig<R> {
    /// Radii for entering and leaving the resident set.
    ///
    /// # Errors
    ///
    /// [`Error::InvalidCellSize`](crate::Error::InvalidCellSize) if either
    /// radius is not finite and positive, or if `unload` is not **strictly**
    /// greater than `load` — equal radii is the thrashing case, and a band of
    /// zero width is not hysteresis.
    pub fn new(load: R, unload: R) -> crate::Result<Self> {
        let bad = |value: R| !value.is_finite() || value <= R::ZERO;
        if bad(load) || bad(unload) || unload <= load {
            return Err(crate::Error::InvalidCellSize {
                value: f64::from(unload.as_f32()),
            });
        }
        Ok(Self { load, unload })
    }

    /// The radius a chunk must come within to be loaded.
    #[must_use]
    pub fn load(&self) -> R {
        self.load
    }

    /// The radius a resident chunk must pass beyond to be unloaded.
    #[must_use]
    pub fn unload(&self) -> R {
        self.unload
    }
}

/// What changed between one [`ChunkStream::update`] and the next.
///
/// Caller-provided and reusable, per `CLAUDE.md` rule 6 — a streaming world
/// calls this every frame.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct StreamUpdate {
    /// Chunks that became resident, ascending by [`ChunkId`].
    pub loaded: Vec<ChunkId>,
    /// Chunks that stopped being resident, ascending by [`ChunkId`].
    pub unloaded: Vec<ChunkId>,
}

impl StreamUpdate {
    /// An empty update.
    #[must_use]
    pub const fn new() -> Self {
        Self {
            loaded: Vec::new(),
            unloaded: Vec::new(),
        }
    }

    /// Clear without releasing capacity.
    pub fn reset(&mut self) {
        self.loaded.clear();
        self.unloaded.clear();
    }

    /// Whether nothing changed.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.loaded.is_empty() && self.unloaded.is_empty()
    }
}

/// The set of chunks currently resident, and the camera-driven rule that
/// maintains it.
#[derive(Clone, Debug, Default)]
pub struct ChunkStream {
    /// Ascending by [`ChunkId`], so membership is a binary search and iteration
    /// order never depends on when a chunk arrived.
    resident: Vec<ChunkId>,
}

impl ChunkStream {
    /// An empty stream.
    #[must_use]
    pub const fn new() -> Self {
        Self {
            resident: Vec::new(),
        }
    }

    /// The resident chunks, ascending by [`ChunkId`].
    #[must_use]
    pub fn resident(&self) -> &[ChunkId] {
        &self.resident
    }

    /// How many chunks are resident.
    #[must_use]
    pub fn len(&self) -> usize {
        self.resident.len()
    }

    /// Whether nothing is resident.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.resident.is_empty()
    }

    /// Drop every chunk, reporting them all as unloaded.
    ///
    /// For a teleport, where "what changed" is the whole set and walking the
    /// hysteresis band would be a lie about locality.
    pub fn clear(&mut self, out: &mut StreamUpdate) {
        out.reset();
        out.unloaded.append(&mut self.resident);
    }

    /// Move the camera and settle the residency set.
    ///
    /// Nothing is meshed here — the result is two lists, and what to do with
    /// them is the caller's. A typical consumer inserts `loaded` into a
    /// [`DirtySet`](super::dirty::DirtySet) and drops the geometry named by
    /// `unloaded`.
    ///
    /// # Errors
    ///
    /// [`Error::IndexSpaceExhausted`](crate::Error::IndexSpaceExhausted) if the
    /// unload radius spans more chunks than `u32` can count. That is a
    /// misconfiguration — a radius of a thousand chunks is a billion candidates
    /// — and it is reported rather than attempted, because the alternative is an
    /// allocation that takes the process down.
    pub fn update<R: Real>(
        &mut self,
        layout: &ChunkLayout<R>,
        camera: [R; 3],
        config: &StreamConfig<R>,
        out: &mut StreamUpdate,
    ) -> crate::Result<()> {
        out.reset();

        // Only chunks within the *unload* radius can be resident afterwards, so
        // that box bounds the whole decision: anything outside it is either
        // already absent or is about to be evicted, and both are handled by the
        // sweep over the current set below.
        let span = self.candidate_span(layout, camera, config.unload)?;
        let mut next: Vec<ChunkId> = Vec::new();

        // `x` outermost, deliberately. `ChunkId` orders lexicographically on
        // `[x, y, z]`, so `x` has to vary *slowest* for this sweep to come out
        // ascending -- the natural z-outer loop produces exactly the wrong
        // order, and the merge below would be silently wrong rather than
        // loudly. The debug assertion after the loop is what caught it.
        for x in span.lo[0]..=span.hi[0] {
            for y in span.lo[1]..=span.hi[1] {
                for z in span.lo[2]..=span.hi[2] {
                    let id = ChunkId::new([x, y, z]);
                    let distance = distance_to_chunk(layout, id, camera);
                    let was = self.resident.binary_search(&id).is_ok();
                    // The band: inside it, a chunk keeps whatever it already is.
                    // That is the hysteresis, and it is why the test reads the
                    // previous state rather than the distance alone.
                    let keep = if was {
                        distance <= config.unload
                    } else {
                        distance <= config.load
                    };
                    if keep {
                        next.push(id);
                    }
                }
            }
        }

        // `next` is already ascending: the loops walk z, then y, then x, and
        // `ChunkId` orders the same way. Asserted rather than assumed, because
        // the diff below is a merge and a merge on unsorted input is silently
        // wrong rather than loudly.
        debug_assert!(next.windows(2).all(|w| w[0] < w[1]));

        // A sorted merge rather than two membership scans: both sides are
        // ascending, so the diff is linear and its output is ascending too.
        let (mut a, mut b) = (0usize, 0usize);
        while a < self.resident.len() || b < next.len() {
            match (self.resident.get(a), next.get(b)) {
                (Some(old), Some(new)) if old == new => {
                    a += 1;
                    b += 1;
                }
                (Some(old), Some(new)) if old < new => {
                    out.unloaded.push(*old);
                    a += 1;
                }
                (Some(_), Some(new)) => {
                    out.loaded.push(*new);
                    b += 1;
                }
                (Some(old), None) => {
                    out.unloaded.push(*old);
                    a += 1;
                }
                (None, Some(new)) => {
                    out.loaded.push(*new);
                    b += 1;
                }
                (None, None) => break,
            }
        }

        self.resident = next;
        Ok(())
    }

    /// The inclusive box of chunk coordinates a radius can reach.
    fn candidate_span<R: Real>(
        &self,
        layout: &ChunkLayout<R>,
        camera: [R; 3],
        radius: R,
    ) -> crate::Result<Span> {
        let extent = R::from_f64(f64::from(layout.cells())) * layout.cell_size();
        let here = layout.chunk_of(camera);
        // How many chunks the radius can span on one side. `+ 1` because the
        // camera sits somewhere inside its own chunk rather than at its corner,
        // so the far face of the last chunk is up to one extent further away.
        let reach = (radius / extent).floor();
        let reach = if reach.is_finite() && reach >= R::ZERO {
            i64::from(reach.as_f32() as i32) + 1
        } else {
            return Err(crate::Error::IndexSpaceExhausted { needed: u64::MAX });
        };

        let side = 2 * reach + 1;
        let total = side.saturating_mul(side).saturating_mul(side);
        if total > i64::from(u32::MAX) {
            return Err(crate::Error::IndexSpaceExhausted {
                needed: total as u64,
            });
        }

        let mut span = Span {
            lo: [0; 3],
            hi: [0; 3],
        };
        for axis in 0..3 {
            let centre = i64::from(here.coords[axis]);
            span.lo[axis] = (centre - reach).clamp(i64::from(i32::MIN), i64::from(i32::MAX)) as i32;
            span.hi[axis] = (centre + reach).clamp(i64::from(i32::MIN), i64::from(i32::MAX)) as i32;
        }
        Ok(span)
    }
}

/// An inclusive box of chunk coordinates.
struct Span {
    lo: [i32; 3],
    hi: [i32; 3],
}

/// Distance from `camera` to the nearest point of chunk `id`.
///
/// Zero when the camera is inside the chunk. This is the axis-aligned box
/// distance rather than the centre distance, for the reason the module docs
/// give: residency is a threshold, and a chunk with geometry inside the radius
/// should be resident whatever its centre is doing.
#[must_use]
pub fn distance_to_chunk<R: Real>(layout: &ChunkLayout<R>, id: ChunkId, camera: [R; 3]) -> R {
    let origin = layout.sample_origin(id);
    let extent = R::from_f64(f64::from(layout.cells())) * layout.cell_size();
    let mut squared = R::ZERO;
    for axis in 0..3 {
        let lo = origin[axis];
        let hi = lo + extent;
        // Outside on the low side, outside on the high side, or inside and
        // contributing nothing.
        let gap = if camera[axis] < lo {
            lo - camera[axis]
        } else if camera[axis] > hi {
            camera[axis] - hi
        } else {
            R::ZERO
        };
        squared += gap * gap;
    }
    squared.sqrt()
}

#[cfg(test)]
mod tests;
