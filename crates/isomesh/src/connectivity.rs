//! Connectivity of the air sublevel set, repaired incrementally as you dig
//! **and** as you fill.
//!
//! Tickets: R-022a (hypothesis P-23), R-022b (hypotheses P-25 and P-26).
//!
//! # The question a game actually asks
//!
//! *Is this cave sealed? Did I just break through? Is this a chokepoint?* None
//! of those is an all-thresholds query about the field. Each is a
//! single-threshold question about the **connected components of the air
//! region** — and it is asked after every edit, at interactive rates.
//!
//! # Why the two directions are not symmetric
//!
//! Durfee, Dhulipala, Kulkarni, Peng, Sawlani & Sun
//! (`10.48550/arXiv.1908.01956`) state the asymmetry:
//!
//! > *"An **insert** can cause at most two trees in `F` to be joined to form a
//! > single tree."*
//! >
//! > *"A **delete** may split a tree into two, but if there exists another edge
//! > between these two resulting trees, they should then be connected together
//! > to ensure that the forest is maximal."*
//!
//! **Digging removes solid, so air only ever appears** — insertion-only, and no
//! replacement-edge search. **Filling removes air**, which can split a component
//! and needs the search.
//!
//! # Why this is a flat label array and not a union-find
//!
//! It was a union-find, and adding `fill` broke it (✗26). Parent pointers encode
//! **union history, not spatial adjacency**: `Q → P → A` records that `Q` was
//! merged in via `P`, and says nothing about whether `Q` touches `P` in the
//! lattice. So a filled sample can be an articulation point *of the tree* while
//! being nothing of the kind *in the graph*, and re-rooting it severs its
//! descendants from a component they are still genuinely part of. The queries
//! that break are the **descendants'**, not the filled sample's own.
//!
//! A **flat** array — every sample holding its component id directly — has no
//! such indirection. Re-rooting a shed piece is one write per member, and no
//! surviving sample can route through it. [`connected`](Air::connected) becomes
//! `O(1)` and takes `&self`.
//!
//! **Flat labels fix the representation, not the search.** The lockstep
//! replacement search was always required; the union-find merely promised
//! falsely that it could be skipped.
//!
//! # What the search costs, and the case that is not cheap
//!
//! Lockstep search expands every seed's frontier in step and stops when all but
//! one exhausts, so it costs the **second-largest** piece rather than the
//! component. M-320 measured the smaller side of a split at **one voxel at the
//! median** over 200 brush fills, 120 at the maximum against 227,567 air
//! samples — which is why the levelled HDT scheme (`O(log n)` levels,
//! non-replacement edges demoted to amortise a failed search) is unnecessary
//! **for that edit distribution**.
//!
//! **It is a property of the distribution, not of the structure.** Fill one
//! voxel at the midpoint of a tunnel joining two equal caverns and both
//! frontiers are huge: the search walks until one exhausts, `O(half the
//! component)`. That edit is not exotic — sealing a passage between two spaces
//! *is* the sealed-volume mechanic. P-26 predicts the flat curve on the measured
//! distribution and a **growing** one on a deliberate bisect, and says so
//! separately, because a structure that came out flat on both would mean the
//! adversarial fixture is not adversarial.
//!
//! # Repair is budgeted, and deferring it fails safe in both directions
//!
//! Amortised is the wrong statistic for the frame a breakthrough lands on
//! (M-124: 20.62 ms unbudgeted against a 2.10 ms budgeted peak), so
//! [`dig`](Air::dig) and [`fill`](Air::fill) take the same
//! `spend: FnMut() -> bool` predicate [`mesh_within_budget`] uses (M-78) and
//! leave unfinished work in [`Air::pending`], drained by
//! [`repair`](Air::repair).
//!
//! **Both directions of staleness are conservative, which is not obvious:**
//!
//! | unfinished | labels say | caller reads |
//! |---|---|---|
//! | `fill` repair | the pieces still share a label | *"not sealed yet"* |
//! | `dig` repair | the labels are still distinct | *"not connected yet"* |
//!
//! Water leaking for three frames is recoverable; water **not** leaking out of a
//! room the engine wrongly believes is sealed is a broken game rule. The
//! conservative answer is the cheap one in both directions.
//!
//! # Cost
//!
//! One `bool`, one `u32` label and one `u32` visit stamp per sample, plus a
//! dense `u32` size per live label and a free list. Sizes are a **`Vec` indexed
//! by label, never a `HashMap`** — map iteration order is a determinism hazard
//! and determinism is load-bearing here (M-36).
//!
//! [`mesh_within_budget`]: crate::chunk::dirty::DirtySet::mesh_within_budget

use alloc::vec::Vec;
use core::fmt;

use crate::{Real, Shape3};

/// No label: the sample is solid.
const NONE: u32 = u32::MAX;

/// What one [`Air::dig`] cost.
///
/// Counts, not durations — a wall-clock ratio is not a gate (✗24).
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
    /// Label writes performed.
    ///
    /// **This replaces the union count M-311 reported**, which stopped naming a
    /// quantity that exists when the structure went flat (✗26). The claim that
    /// row measured — repair proportional to the edit, not the lattice — is
    /// unaffected; only the unit changed.
    pub relabels: u64,
    /// Components absorbed into another.
    ///
    /// The common case is zero: most digging widens what is already connected
    /// rather than joining two things.
    pub merges: u64,
}

impl Repair {
    /// Label writes per newly-air sample, or `0.0` when nothing was dug.
    #[must_use]
    pub fn relabels_per_dirty(&self) -> f64 {
        if self.dirty == 0 {
            0.0
        } else {
            self.relabels as f64 / self.dirty as f64
        }
    }
}

impl fmt::Display for Repair {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "dirty {} (+{} already air), {} relabels ({:.2}/dirty), {} merges",
            self.dirty,
            self.already_air,
            self.relabels,
            self.relabels_per_dirty(),
            self.merges
        )
    }
}

/// What one [`Air::fill`] cost.
///
/// [`visited`](Self::visited) is the quantity P-26 predicts: the voxels the
/// lockstep replacement search touched. It is the cost, and it is the number
/// that must stay flat as the lattice grows.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Default)]
pub struct Fill {
    /// Samples that were air and are now solid.
    pub dirty: u64,
    /// Samples the edit touched that were already solid.
    pub already_solid: u64,
    /// Surviving air samples adjacent to something this fill removed.
    ///
    /// The search's starting points. A component needs **two or more** before it
    /// can possibly have split — one seed cannot be separated from itself.
    pub seeds: u64,
    /// Voxels the lockstep search touched.
    ///
    /// **P-26's cost quantity.** Bounded by the second-largest piece, not by the
    /// component: the search stops when all but one frontier exhausts, and the
    /// surviving frontier is never walked to completion.
    pub visited: u64,
    /// Components that shed at least one piece.
    pub splits: u64,
    /// New components created by splitting.
    pub shed: u64,
    /// Components that lost their last air sample.
    ///
    /// Consumed outright rather than severed. These need no search at all, which
    /// is why they are counted apart from [`splits`](Self::splits).
    pub vanished: u64,
    /// Components still awaiting a search when the budget ran out.
    ///
    /// Non-zero means [`Air::connected`] may report two severed pieces as
    /// connected — never the reverse. See the module docs.
    pub pending: u64,
}

impl fmt::Display for Fill {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "dirty {} (+{} already solid), {} seeds, {} visited, \
             {} splits shedding {}, {} vanished, {} pending",
            self.dirty,
            self.already_solid,
            self.seeds,
            self.visited,
            self.splits,
            self.shed,
            self.vanished,
            self.pending
        )
    }
}

/// Connected components of the air sublevel set, maintained under digging and
/// filling.
///
/// Air is `value >= 0`, the complement of the solid convention every extractor
/// here uses. Build once from a sampled field, then [`dig`](Self::dig),
/// [`fill`](Self::fill), and ask [`connected`](Self::connected).
#[derive(Clone, Debug)]
pub struct Air {
    /// `true` where the sample is air.
    air: Vec<bool>,
    /// Component id per sample, [`NONE`] where solid.
    label: Vec<u32>,
    /// Air samples per label, indexed by label. Dense, never a map (M-36).
    size: Vec<u32>,
    /// Retired label ids, reusable.
    free: Vec<u32>,
    /// Live component count, maintained rather than counted.
    live: u64,
    /// Labels whose component may have been severed and not yet searched.
    pending: Vec<u32>,
    /// Visit stamps for the search, avoiding an `O(n)` clear per call.
    stamp: Vec<u32>,
    /// Current stamp generation.
    epoch: u32,
    dims: [u32; 3],
    /// Reusable scratch, so a repair allocates nothing.
    queue: Vec<usize>,
    scratch: Vec<u32>,
}

impl Air {
    /// Build from sampled values, air being `value >= 0`.
    ///
    /// `O(n)`, one flood fill. The returned [`Repair`] describes that build, so
    /// the same instrument measures the rebuild and the incremental path.
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
            label: alloc::vec![NONE; count],
            size: Vec::new(),
            free: Vec::new(),
            live: 0,
            pending: Vec::new(),
            stamp: alloc::vec![0u32; count],
            epoch: 0,
            dims,
            queue: Vec::new(),
            scratch: Vec::new(),
        };

        let mut repair = Repair::default();
        for start in 0..count {
            if me.air.get(start) != Some(&true) || me.label.get(start) != Some(&NONE) {
                continue;
            }
            let l = me.take_label();
            let grown = me.flood(start, l);
            me.set_size(l, grown);
            repair.dirty += grown as u64;
            repair.relabels += grown as u64;
        }
        Ok((me, repair))
    }

    /// Turn a set of samples to air and repair connectivity.
    ///
    /// `spend` is polled once per unit of relabelling work and returning `false`
    /// stops the repair, leaving the remainder in [`pending`](Self::pending).
    /// Unfinished digging reads as *"not connected yet"*, which is the safe
    /// direction — see the module docs. Pass `|| true` for a synchronous repair.
    ///
    /// Out-of-range coordinates are ignored rather than rejected: a brush
    /// straddling the grid edge is ordinary, not an error.
    pub fn dig<B: FnMut() -> bool>(&mut self, samples: &[[u32; 3]], mut spend: B) -> Repair {
        let mut repair = Repair::default();

        // Two passes, and the order is load-bearing. Marking every sample air
        // first means the second pass sees the *finished* phase field, so two
        // newly-air neighbours in one batch are joined. Interleaving them would
        // make the result depend on the order `samples` happens to be in.
        let mut fresh = core::mem::take(&mut self.scratch);
        fresh.clear();
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
                    fresh.push(i as u32);
                }
                None => {}
            }
        }

        // Sorted, so the labelling does not depend on the order the caller
        // happened to list the brush in. Determinism is load-bearing (M-36).
        fresh.sort_unstable();
        fresh.dedup();

        for k in 0..fresh.len() {
            let Some(&i) = fresh.get(k) else { continue };
            let i = i as usize;
            if self.label.get(i) != Some(&NONE) {
                continue;
            }
            if !spend() {
                break;
            }
            self.grow_from(i, &mut repair);
        }

        self.scratch = fresh;
        repair
    }

    /// Turn a set of samples to solid and repair connectivity.
    ///
    /// `spend` is polled once per component searched; returning `false` leaves
    /// the rest in [`pending`](Self::pending), reported as
    /// [`Fill::pending`]. Unfinished filling reads as *"not sealed yet"*, which
    /// is the safe direction — see the module docs. Pass `|| true` for a
    /// synchronous repair.
    pub fn fill<B: FnMut() -> bool>(&mut self, samples: &[[u32; 3]], mut spend: B) -> Fill {
        let mut out = Fill::default();

        let mut gone = core::mem::take(&mut self.scratch);
        gone.clear();
        for s in samples {
            if !self.in_range(*s) {
                continue;
            }
            let i = self.index(*s);
            match self.air.get_mut(i) {
                Some(a) if *a => {
                    *a = false;
                    out.dirty += 1;
                    gone.push(i as u32);
                }
                Some(_) => out.already_solid += 1,
                None => {}
            }
        }
        gone.sort_unstable();
        gone.dedup();

        // Retire the removed samples first, so sizes are correct before any
        // search reads them.
        for &g in &gone {
            let i = g as usize;
            let Some(&l) = self.label.get(i) else {
                continue;
            };
            if l == NONE {
                continue;
            }
            if let Some(slot) = self.label.get_mut(i) {
                *slot = NONE;
            }
            let remaining = self
                .size
                .get(l as usize)
                .copied()
                .unwrap_or(0)
                .saturating_sub(1);
            if let Some(slot) = self.size.get_mut(l as usize) {
                *slot = remaining;
            }
            if remaining == 0 {
                self.retire(l);
                out.vanished += 1;
            }
        }

        // Every surviving air sample adjacent to something removed is a seed.
        // A component with fewer than two cannot have been severed: one seed
        // cannot be separated from itself.
        let mut nb = [0usize; 6];
        for &g in &gone {
            let used = self.neighbours(g as usize, &mut nb);
            for &j in nb.iter().take(used) {
                if self.air.get(j) != Some(&true) {
                    continue;
                }
                let Some(&l) = self.label.get(j) else {
                    continue;
                };
                if l == NONE {
                    continue;
                }
                out.seeds += 1;
                if !self.pending.contains(&l) {
                    self.pending.push(l);
                }
            }
        }
        self.scratch = gone;

        let drained = self.repair(&mut spend);
        out.visited += drained.visited;
        out.splits += drained.splits;
        out.shed += drained.shed;
        out.pending = self.pending.len() as u64;
        out
    }

    /// Drain deferred repair work left by a budgeted [`dig`](Self::dig) or
    /// [`fill`](Self::fill).
    ///
    /// Returns what this call accomplished. [`Fill::pending`] on the result says
    /// what is still outstanding.
    pub fn repair<B: FnMut() -> bool>(&mut self, spend: &mut B) -> Fill {
        let mut out = Fill::default();
        while let Some(&l) = self.pending.first() {
            if !spend() {
                break;
            }
            self.pending.remove(0);
            self.search(l, &mut out);
        }
        out.pending = self.pending.len() as u64;
        out
    }

    /// Are these two samples in the same air component?
    ///
    /// `O(1)`. `false` if either is solid or out of range — solid is not air,
    /// and a sample that does not exist is not connected to anything.
    ///
    /// While [`pending`](Self::pending) is non-zero this may answer `true` for
    /// two samples a fill has already severed. It never answers `false` for two
    /// that are connected. See the module docs.
    #[must_use]
    pub fn connected(&self, a: [u32; 3], b: [u32; 3]) -> bool {
        if !self.in_range(a) || !self.in_range(b) {
            return false;
        }
        let (ia, ib) = (self.index(a), self.index(b));
        match (self.label.get(ia), self.label.get(ib)) {
            (Some(&x), Some(&y)) => x != NONE && x == y,
            _ => false,
        }
    }

    /// The component id of a sample, or `None` if it is solid or out of range.
    ///
    /// Exposed so a higher layer can stitch several `Air`s into one world: read
    /// the labels along a shared face and join the pairs that touch. That is the
    /// per-chunk decomposition, and it is what bounds the bisect tail the module
    /// docs describe — a search inside one `Air` cannot cost more than that
    /// `Air`. Building it needs nothing from this type but this accessor.
    #[must_use]
    pub fn label_of(&self, p: [u32; 3]) -> Option<u32> {
        if !self.in_range(p) {
            return None;
        }
        match self.label.get(self.index(p)) {
            Some(&l) if l != NONE => Some(l),
            _ => None,
        }
    }

    /// How many air components there are.
    ///
    /// `O(1)` — maintained, not counted.
    #[must_use]
    pub fn components(&self) -> u64 {
        self.live
    }

    /// Air samples.
    #[must_use]
    pub fn air_samples(&self) -> u64 {
        self.air.iter().filter(|a| **a).count() as u64
    }

    /// Components awaiting a deferred search.
    #[must_use]
    pub fn pending(&self) -> u64 {
        self.pending.len() as u64
    }

    // --- internals ---------------------------------------------------------

    fn in_range(&self, p: [u32; 3]) -> bool {
        p[0] < self.dims[0] && p[1] < self.dims[1] && p[2] < self.dims[2]
    }

    fn index(&self, p: [u32; 3]) -> usize {
        (p[2] as usize * self.dims[1] as usize + p[1] as usize) * self.dims[0] as usize
            + p[0] as usize
    }

    /// The six axis-aligned neighbours of `i`, skipping any off the lattice.
    fn neighbours(&self, i: usize, out: &mut [usize; 6]) -> usize {
        let (nx, ny) = (self.dims[0] as usize, self.dims[1] as usize);
        let nz = self.dims[2] as usize;
        if nx == 0 || ny == 0 || nz == 0 {
            return 0;
        }
        let x = i % nx;
        let y = (i / nx) % ny;
        let z = i / (nx * ny);
        let mut count = 0;
        let mut push = |v: usize| {
            if let Some(slot) = out.get_mut(count) {
                *slot = v;
                count += 1;
            }
        };
        if x > 0 {
            push(i - 1);
        }
        if x + 1 < nx {
            push(i + 1);
        }
        if y > 0 {
            push(i - nx);
        }
        if y + 1 < ny {
            push(i + nx);
        }
        if z > 0 {
            push(i - nx * ny);
        }
        if z + 1 < nz {
            push(i + nx * ny);
        }
        count
    }

    fn take_label(&mut self) -> u32 {
        self.live += 1;
        match self.free.pop() {
            Some(l) => l,
            None => {
                self.size.push(0);
                (self.size.len() - 1) as u32
            }
        }
    }

    fn retire(&mut self, l: u32) {
        self.live = self.live.saturating_sub(1);
        if let Some(slot) = self.size.get_mut(l as usize) {
            *slot = 0;
        }
        self.free.push(l);
        self.pending.retain(|p| *p != l);
    }

    fn set_size(&mut self, l: u32, n: u32) {
        if let Some(slot) = self.size.get_mut(l as usize) {
            *slot = n;
        }
    }

    fn size_of(&self, l: u32) -> u32 {
        self.size.get(l as usize).copied().unwrap_or(0)
    }

    /// Label every air sample reachable from `start` with `l`. Returns the count.
    fn flood(&mut self, start: usize, l: u32) -> u32 {
        let mut queue = core::mem::take(&mut self.queue);
        queue.clear();
        queue.push(start);
        if let Some(slot) = self.label.get_mut(start) {
            *slot = l;
        }
        let mut seen = 0u32;
        let mut head = 0;
        let mut nb = [0usize; 6];
        while let Some(&i) = queue.get(head) {
            head += 1;
            seen += 1;
            let used = self.neighbours(i, &mut nb);
            for &j in nb.iter().take(used) {
                if self.air.get(j) != Some(&true) || self.label.get(j) != Some(&NONE) {
                    continue;
                }
                if let Some(slot) = self.label.get_mut(j) {
                    *slot = l;
                }
                queue.push(j);
            }
        }
        self.queue = queue;
        seen
    }

    /// Relabel every air sample reachable from `start` currently carrying
    /// `from`, writing `to`. Returns the count.
    fn recolour(&mut self, start: usize, from: u32, to: u32) -> u32 {
        if from == to {
            return 0;
        }
        let mut queue = core::mem::take(&mut self.queue);
        queue.clear();
        queue.push(start);
        if let Some(slot) = self.label.get_mut(start) {
            *slot = to;
        }
        let mut seen = 0u32;
        let mut head = 0;
        let mut nb = [0usize; 6];
        while let Some(&i) = queue.get(head) {
            head += 1;
            seen += 1;
            let used = self.neighbours(i, &mut nb);
            for &j in nb.iter().take(used) {
                if self.label.get(j) != Some(&from) {
                    continue;
                }
                if let Some(slot) = self.label.get_mut(j) {
                    *slot = to;
                }
                queue.push(j);
            }
        }
        self.queue = queue;
        seen
    }

    /// Label the newly-air blob containing `i` and absorb the components it
    /// touches, keeping the largest so the relabelling is union-by-size.
    fn grow_from(&mut self, i: usize, repair: &mut Repair) {
        // Claim the blob with a provisional label, collecting the existing
        // labels it touches.
        let provisional = self.take_label();
        let mut touched = Vec::new();
        let mut queue = core::mem::take(&mut self.queue);
        queue.clear();
        queue.push(i);
        if let Some(slot) = self.label.get_mut(i) {
            *slot = provisional;
        }
        let mut grown = 0u32;
        let mut head = 0;
        let mut nb = [0usize; 6];
        while let Some(&at) = queue.get(head) {
            head += 1;
            grown += 1;
            let used = self.neighbours(at, &mut nb);
            for &j in nb.iter().take(used) {
                if self.air.get(j) != Some(&true) {
                    continue;
                }
                match self.label.get(j).copied() {
                    Some(NONE) => {
                        if let Some(slot) = self.label.get_mut(j) {
                            *slot = provisional;
                        }
                        queue.push(j);
                    }
                    Some(l) if l != provisional && !touched.contains(&l) => touched.push(l),
                    _ => {}
                }
            }
        }
        self.queue = queue;
        self.set_size(provisional, grown);
        repair.relabels += u64::from(grown);

        // Union by size: whichever component is largest keeps its label, and
        // every other one is rewritten into it. Ties go to the smaller id so the
        // outcome does not depend on discovery order (M-36).
        let mut best = provisional;
        for &l in &touched {
            let (sl, sb) = (self.size_of(l), self.size_of(best));
            if sl > sb || (sl == sb && l < best) {
                best = l;
            }
        }

        let mut merged = 0u32;
        for &l in &touched {
            if l == best {
                continue;
            }
            let Some(seed) = self.any_member(l) else {
                continue;
            };
            let moved = self.recolour(seed, l, best);
            repair.relabels += u64::from(moved);
            repair.merges += 1;
            merged += moved;
            self.retire(l);
        }
        if best != provisional {
            // The newly-dug blob is itself a component being absorbed, so this
            // is a join and counts as one. Only when the blob is the largest
            // does it keep its label and absorb the others instead.
            let Some(seed) = self.any_member(provisional) else {
                return;
            };
            let moved = self.recolour(seed, provisional, best);
            repair.relabels += u64::from(moved);
            repair.merges += 1;
            merged += moved;
            self.retire(provisional);
        }
        let total = self.size_of(best) + merged;
        self.set_size(best, total);
    }

    /// Any air sample carrying `l`, or `None`. Linear, and used only where the
    /// component is about to be walked anyway.
    fn any_member(&self, l: u32) -> Option<usize> {
        self.label.iter().position(|x| *x == l)
    }

    /// Lockstep replacement search over the component labelled `l`.
    ///
    /// Every seed grows a frontier one voxel per round; frontiers that meet
    /// merge. The search stops when at most one frontier is still active, so the
    /// surviving piece is **never walked to completion** and the cost is the
    /// second-largest piece rather than the component (P-26).
    fn search(&mut self, l: u32, out: &mut Fill) {
        let mut seeds = Vec::new();
        let mut nb = [0usize; 6];
        // A seed is an air sample of `l` adjacent to solid: only those can have
        // been separated by a removal.
        for i in 0..self.label.len() {
            if self.label.get(i) != Some(&l) {
                continue;
            }
            let used = self.neighbours(i, &mut nb);
            let touches_solid = used < 6
                || nb
                    .iter()
                    .take(used)
                    .any(|&j| self.air.get(j) != Some(&true));
            if touches_solid {
                seeds.push(i);
            }
        }
        if seeds.len() < 2 {
            return;
        }

        self.epoch = self.epoch.wrapping_add(1);
        if self.epoch == 0 {
            self.stamp.iter_mut().for_each(|s| *s = 0);
            self.epoch = 1;
        }
        let epoch = self.epoch;

        // One frontier per seed, each with its own queue. `owner` maps a visited
        // voxel to the frontier that claimed it; `into` resolves frontier merges.
        let mut owner: Vec<u32> = alloc::vec![NONE; self.label.len().min(seeds.len().max(1))];
        owner.clear();
        let mut claim: Vec<u32> = alloc::vec![NONE; self.label.len()];
        let mut queues: Vec<Vec<usize>> = Vec::new();
        let mut into: Vec<u32> = Vec::new();
        let mut done: Vec<bool> = Vec::new();
        for (f, &s) in seeds.iter().enumerate() {
            if self.stamp.get(s) == Some(&epoch) {
                // Already claimed by an earlier seed's frontier this round.
                queues.push(Vec::new());
                into.push(f as u32);
                done.push(true);
                continue;
            }
            if let Some(slot) = self.stamp.get_mut(s) {
                *slot = epoch;
            }
            if let Some(slot) = claim.get_mut(s) {
                *slot = f as u32;
            }
            queues.push(alloc::vec![s]);
            into.push(f as u32);
            done.push(false);
        }

        let resolve = |into: &Vec<u32>, mut f: u32| -> u32 {
            while into.get(f as usize).copied().unwrap_or(f) != f {
                f = into.get(f as usize).copied().unwrap_or(f);
            }
            f
        };

        loop {
            let active: Vec<u32> = (0..queues.len() as u32)
                .filter(|&f| resolve(&into, f) == f && done.get(f as usize) != Some(&true))
                .collect();
            if active.len() <= 1 {
                break;
            }
            for &f in &active {
                if resolve(&into, f) != f || done.get(f as usize) == Some(&true) {
                    continue;
                }
                let Some(q) = queues.get_mut(f as usize) else {
                    continue;
                };
                let Some(at) = q.pop() else {
                    if let Some(slot) = done.get_mut(f as usize) {
                        *slot = true;
                    }
                    continue;
                };
                out.visited += 1;
                let used = self.neighbours(at, &mut nb);
                for &j in nb.iter().take(used) {
                    if self.label.get(j) != Some(&l) {
                        continue;
                    }
                    if self.stamp.get(j) == Some(&epoch) {
                        let g = claim.get(j).copied().unwrap_or(NONE);
                        if g == NONE {
                            continue;
                        }
                        let (rf, rg) = (resolve(&into, f), resolve(&into, g));
                        if rf != rg {
                            // Two frontiers met: the pieces are one piece.
                            let (keep, gone) = if rf < rg { (rf, rg) } else { (rg, rf) };
                            if let Some(slot) = into.get_mut(gone as usize) {
                                *slot = keep;
                            }
                            let moved = queues.get(gone as usize).cloned().unwrap_or_default();
                            if let Some(dst) = queues.get_mut(keep as usize) {
                                dst.extend(moved);
                            }
                            if let Some(slot) = queues.get_mut(gone as usize) {
                                slot.clear();
                            }
                            if done.get(gone as usize) == Some(&true)
                                && let Some(slot) = done.get_mut(keep as usize)
                            {
                                *slot = false;
                            }
                        }
                        continue;
                    }
                    if let Some(slot) = self.stamp.get_mut(j) {
                        *slot = epoch;
                    }
                    if let Some(slot) = claim.get_mut(j) {
                        *slot = resolve(&into, f);
                    }
                    if let Some(q) = queues.get_mut(resolve(&into, f) as usize) {
                        q.push(j);
                    }
                }
            }
        }

        // Frontiers that exhausted are complete pieces. Whichever frontier is
        // still active (or, if all finished, the largest) keeps `l`.
        let roots: Vec<u32> = (0..queues.len() as u32)
            .filter(|&f| resolve(&into, f) == f)
            .collect();
        if roots.len() < 2 {
            return;
        }
        let unfinished: Vec<u32> = roots
            .iter()
            .copied()
            .filter(|&f| done.get(f as usize) != Some(&true))
            .collect();

        let mut members: Vec<(u32, Vec<usize>)> = Vec::new();
        for &f in &roots {
            if unfinished.len() == 1 && unfinished.first() == Some(&f) {
                continue;
            }
            let mut m = Vec::new();
            for i in 0..self.label.len() {
                if self.label.get(i) == Some(&l)
                    && self.stamp.get(i) == Some(&epoch)
                    && claim.get(i).map(|c| resolve(&into, *c)) == Some(f)
                {
                    m.push(i);
                }
            }
            if !m.is_empty() {
                members.push((f, m));
            }
        }
        if members.is_empty() {
            return;
        }
        // If every frontier finished, the largest keeps the old label.
        if unfinished.is_empty() {
            let mut biggest = 0;
            for (k, (_, m)) in members.iter().enumerate() {
                if m.len() > members.get(biggest).map_or(0, |b| b.1.len()) {
                    biggest = k;
                }
            }
            if biggest < members.len() {
                members.remove(biggest);
            }
        }
        if members.is_empty() {
            return;
        }

        out.splits += 1;
        for (_, m) in &members {
            let fresh = self.take_label();
            for &i in m {
                if let Some(slot) = self.label.get_mut(i) {
                    *slot = fresh;
                }
            }
            self.set_size(fresh, m.len() as u32);
            let left = self.size_of(l).saturating_sub(m.len() as u32);
            self.set_size(l, left);
            out.shed += 1;
        }
    }
}

#[cfg(test)]
mod tests;
