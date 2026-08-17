//! Connectivity of the air sublevel set, repaired incrementally as you dig.
//!
//! Ticket: R-022a, hypothesis P-23.
//!
//! # The question a game actually asks
//!
//! *Is this cave sealed? Did I just break through? Is this a chokepoint?* None
//! of those is an all-thresholds query about the field. Each is a
//! single-threshold question about the **connected components of the air
//! region** — and it is asked after every edit, at interactive rates.
//!
//! # Why digging is the easy direction, and it is not a small difference
//!
//! Durfee, Dhulipala, Kulkarni, Peng, Sawlani & Sun
//! (`10.48550/arXiv.1908.01956`) state the asymmetry that decides this module's
//! whole design:
//!
//! > *"An **insert** can cause at most two trees in `F` to be joined to form a
//! > single tree."*
//! >
//! > *"A **delete** may split a tree into two, but if there exists another edge
//! > between these two resulting trees, they should then be connected together
//! > to ensure that the forest is maximal."*
//!
//! **Digging removes solid, so air samples and air-air edges only ever
//! appear.** That is insertion-only, an insert never needs a replacement-edge
//! search, and a union-find is the entire structure — near-constant amortised,
//! no logarithm, no sketching.
//!
//! **Filling is the other direction and this module does not do it.** Removing
//! air needs deletion, deletion needs the replacement search, and a union-find
//! cannot express it at any price. That is R-022b, and a `fill` method does not
//! exist rather than existing and being slow — a `dig`-only API is the honest
//! shape for a `dig`-only structure.
//!
//! # What is counted, and why it is a count
//!
//! [`Repair`] reports **union calls** and **effective merges**, which are
//! integers and identical on every machine. A wall-clock ratio is not a gate
//! (✗24): the deterministic quantity is the one to assert on, and the timing is
//! a print beside it.
//!
//! # Cost
//!
//! One `u32` parent and one `u32` size per sample, plus a byte of phase. Path
//! halving and union by size, so a sequence of `m` operations over `n` samples
//! is `O(m · α(n))`. The lattice is **6-connected**, so a newly-air sample
//! contributes at most six union calls — which is the bound P-23's second
//! falsifier watches.

use alloc::vec::Vec;
use core::fmt;

use crate::{Real, Shape3};

/// What one [`Air::dig`] cost.
///
/// Counts, not durations. See the module docs.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Default)]
pub struct Repair {
    /// Samples that were solid and are now air.
    ///
    /// The **dirty set**. Everything else here is read against this.
    pub dirty: u64,
    /// Samples the edit touched that were already air.
    ///
    /// Digging where you have already dug. Reported because a brush applied
    /// twice has a dirty set of zero and should cost nothing, and a harness that
    /// cannot see the difference cannot show that.
    pub already_air: u64,
    /// `union` calls made.
    ///
    /// At most `6 · dirty`, because each newly-air sample offers its six
    /// incident lattice edges and nothing else is visited. **P-23's second
    /// falsifier is this exceeding that bound** — it fails when the *instrument*
    /// is wrong rather than when the world is.
    pub unions: u64,
    /// Of those, the ones that actually merged two different components.
    ///
    /// The rest found both ends already in one component, which is the common
    /// case: most digging does not change connectivity, it just widens what is
    /// already connected.
    pub merges: u64,
}

impl Repair {
    /// Union calls per newly-air sample, or `0.0` when nothing was dug.
    #[must_use]
    pub fn unions_per_dirty(&self) -> f64 {
        if self.dirty == 0 {
            0.0
        } else {
            self.unions as f64 / self.dirty as f64
        }
    }
}

impl fmt::Display for Repair {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "dirty {} (+{} already air), {} unions ({:.2}/dirty), {} merges",
            self.dirty,
            self.already_air,
            self.unions,
            self.unions_per_dirty(),
            self.merges
        )
    }
}

/// Connected components of the air sublevel set, maintained under digging.
///
/// Air is `value >= 0`, the complement of the solid convention every extractor
/// here uses. Build once from a sampled field, then [`dig`](Self::dig) and ask
/// [`connected`](Self::connected).
///
/// # Digging only
///
/// There is deliberately no `fill`. See the module docs: removing air is a
/// deletion, and a union-find cannot do deletions. An API that offered one and
/// silently rebuilt would be the second execution path this crate's rules
/// forbid.
#[derive(Clone, Debug)]
pub struct Air {
    /// `true` where the sample is air.
    air: Vec<bool>,
    parent: Vec<u32>,
    size: Vec<u32>,
    dims: [u32; 3],
}

impl Air {
    /// Build from sampled values, air being `value >= 0`.
    ///
    /// `O(n)` union calls, which is the cost this exists to be compared against.
    /// The returned [`Repair`] describes that build, so the same instrument
    /// measures the rebuild and the incremental path.
    ///
    /// # Errors
    ///
    /// [`Error::ShapeOverflow`](crate::Error::ShapeOverflow) if `values` is not
    /// exactly one entry per sample.
    pub fn build<R: Real>(values: &[R], shape: &impl Shape3) -> crate::Result<(Self, Repair)> {
        let dims = shape.size();
        let count = shape.element_count();
        if values.len() != count {
            return Err(crate::Error::ShapeOverflow {
                size: dims,
                product: values.len() as u64,
            });
        }

        let mut me = Self {
            air: values.iter().map(|v| !crate::cube::is_inside(*v)).collect(),
            parent: (0..count as u32).collect(),
            size: alloc::vec![1u32; count],
            dims,
        };

        // Every air sample offers its three forward edges, so each lattice edge
        // is visited once rather than twice. Doing it from both ends would
        // double `unions` without changing a single component, which is exactly
        // the instrument error P-23's second falsifier watches for.
        let mut repair = Repair::default();
        for z in 0..dims[2] {
            for y in 0..dims[1] {
                for x in 0..dims[0] {
                    let i = me.index([x, y, z]);
                    if me.air.get(i) != Some(&true) {
                        continue;
                    }
                    repair.dirty += 1;
                    me.link_forward([x, y, z], &mut repair);
                }
            }
        }
        Ok((me, repair))
    }

    /// Turn a set of samples to air and repair connectivity.
    ///
    /// Samples already air are counted in
    /// [`already_air`](Repair::already_air) and cost nothing, so applying the
    /// same brush twice repairs nothing the second time.
    ///
    /// Out-of-range coordinates are ignored rather than rejected: a brush
    /// straddling the grid edge is ordinary, not an error.
    pub fn dig(&mut self, samples: &[[u32; 3]]) -> Repair {
        let mut repair = Repair::default();

        // Two passes, and the order is load-bearing. Marking every sample air
        // first means the second pass sees the *finished* phase field, so two
        // newly-air neighbours in one batch are joined. Interleaving them would
        // make the result depend on the order `samples` happens to be in.
        for s in samples {
            if !self.in_range(*s) {
                continue;
            }
            let i = self.index(*s);
            match self.air.get_mut(i) {
                Some(a) if *a => repair.already_air += 1,
                Some(a) => {
                    *a = true;
                    repair.dirty += 1;
                }
                None => {}
            }
        }

        for s in samples {
            if !self.in_range(*s) {
                continue;
            }
            // Both directions here, unlike `build`: a newly-air sample's
            // backward neighbours are *not* going to visit it, because they are
            // not in the dirty set.
            self.link_all(*s, &mut repair);
        }
        repair
    }

    /// Are these two samples in the same air component?
    ///
    /// `false` if either is solid or out of range — solid is not air, and a
    /// sample that does not exist is not connected to anything.
    pub fn connected(&mut self, a: [u32; 3], b: [u32; 3]) -> bool {
        if !self.in_range(a) || !self.in_range(b) {
            return false;
        }
        let (ia, ib) = (self.index(a), self.index(b));
        if self.air.get(ia) != Some(&true) || self.air.get(ib) != Some(&true) {
            return false;
        }
        self.find(ia) == self.find(ib)
    }

    /// How many air components there are.
    ///
    /// `O(n)`; this is a diagnostic, not a query to call per frame.
    pub fn components(&mut self) -> u64 {
        let mut n = 0;
        for i in 0..self.air.len() {
            if self.air.get(i) == Some(&true) && self.find(i) == i {
                n += 1;
            }
        }
        n
    }

    /// Air samples.
    #[must_use]
    pub fn air_samples(&self) -> u64 {
        self.air.iter().filter(|a| **a).count() as u64
    }

    fn in_range(&self, p: [u32; 3]) -> bool {
        p[0] < self.dims[0] && p[1] < self.dims[1] && p[2] < self.dims[2]
    }

    fn index(&self, p: [u32; 3]) -> usize {
        (p[2] as usize * self.dims[1] as usize + p[1] as usize) * self.dims[0] as usize
            + p[0] as usize
    }

    /// Join to the three neighbours at `+x`, `+y`, `+z`.
    fn link_forward(&mut self, p: [u32; 3], repair: &mut Repair) {
        for axis in 0..3 {
            let mut q = p;
            q[axis] += 1;
            self.try_link(p, q, repair);
        }
    }

    /// Join to all six neighbours.
    fn link_all(&mut self, p: [u32; 3], repair: &mut Repair) {
        for axis in 0..3 {
            let mut q = p;
            q[axis] += 1;
            self.try_link(p, q, repair);
            let mut q = p;
            let Some(lower) = q[axis].checked_sub(1) else {
                continue;
            };
            q[axis] = lower;
            self.try_link(p, q, repair);
        }
    }

    fn try_link(&mut self, p: [u32; 3], q: [u32; 3], repair: &mut Repair) {
        if !self.in_range(q) {
            return;
        }
        let (i, j) = (self.index(p), self.index(q));
        if self.air.get(j) != Some(&true) {
            return;
        }
        repair.unions += 1;
        if self.union(i, j) {
            repair.merges += 1;
        }
    }

    fn find(&mut self, mut i: usize) -> usize {
        while self.parent.get(i).copied().map(|p| p as usize) != Some(i) {
            let Some(p) = self.parent.get(i).copied() else {
                return i;
            };
            let grand = self.parent.get(p as usize).copied().unwrap_or(p);
            if let Some(slot) = self.parent.get_mut(i) {
                *slot = grand;
            }
            i = grand as usize;
        }
        i
    }

    /// `true` if this merged two different components.
    fn union(&mut self, a: usize, b: usize) -> bool {
        let (mut ra, mut rb) = (self.find(a), self.find(b));
        if ra == rb {
            return false;
        }
        if self.size.get(ra).copied().unwrap_or(0) < self.size.get(rb).copied().unwrap_or(0) {
            core::mem::swap(&mut ra, &mut rb);
        }
        if let Some(slot) = self.parent.get_mut(rb) {
            *slot = ra as u32;
        }
        let grew = self.size.get(rb).copied().unwrap_or(0);
        if let Some(slot) = self.size.get_mut(ra) {
            *slot += grew;
        }
        true
    }
}

#[cfg(test)]
mod tests;
