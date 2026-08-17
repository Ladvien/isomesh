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
    /// Air–solid faces per label, indexed by label — the boundary surface
    /// area a Sabine estimate needs (R-036), in face units (× h² for world
    /// area). Domain-boundary faces count as solid: the sealed-box
    /// convention, because a component touching the grid edge is a stitching
    /// layer's problem, not this grid's. Delta-maintained everywhere `label`
    /// is; under a budget-truncated op it is conservatively stale exactly as
    /// `label` is, and exact again once drained.
    area: Vec<u32>,
    /// Retired label ids, reusable.
    free: Vec<u32>,
    /// Live component count, maintained rather than counted.
    live: u64,
    /// Components that may have been severed, with the surviving samples
    /// adjacent to the removal that caused it.
    ///
    /// The seeds are carried rather than recomputed. Rederiving them as "air
    /// adjacent to solid" collects the whole cave surface instead of the
    /// fill's own neighbourhood, which is thousands of frontiers rather than a
    /// handful and makes the search explore the entire component — measured,
    /// and it falsified P-26 on the first run.
    pending: Vec<(u32, Vec<u32>)>,
    /// Visit stamps for the search, avoiding an `O(n)` clear per call.
    stamp: Vec<u32>,
    /// Current stamp generation.
    epoch: u32,
    /// Which frontier claimed a sample; meaningful only where `stamp == epoch`.
    claim: Vec<u32>,
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
            area: Vec::new(),
            free: Vec::new(),
            live: 0,
            pending: Vec::new(),
            stamp: alloc::vec![0u32; count],
            epoch: 0,
            claim: alloc::vec![NONE; count],
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
        // Face counts in one pass over the finished labelling — O(6n), once.
        for i in 0..count {
            if me.air.get(i) != Some(&true) {
                continue;
            }
            let Some(&l) = me.label.get(i) else { continue };
            if l == NONE {
                continue;
            }
            let solid = me.solid_faces(i);
            me.area_add(l, solid);
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
        let mut nb_area = [0usize; 6];
        for &g in &gone {
            let i = g as usize;
            let Some(&l) = self.label.get(i) else {
                continue;
            };
            if l == NONE {
                continue;
            }
            // Face bookkeeping against the finished phase field (every gone
            // sample is already marked solid). A face to still-air gains a
            // boundary count on that side; a face to originally-solid (or off
            // the lattice) stops being this component's boundary; a face to
            // another gone sample was air–air and is now solid–solid — zero
            // either way (R-036).
            let used = self.neighbours(i, &mut nb_area);
            for &n in nb_area.iter().take(used) {
                if self.air.get(n) == Some(&true) {
                    if let Some(&ln) = self.label.get(n)
                        && ln != NONE
                    {
                        self.area_add(ln, 1);
                    }
                } else if gone.binary_search(&(n as u32)).is_err() {
                    self.area_sub(l, 1);
                }
            }
            self.area_sub(l, 6 - used as u32);
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
                match self.pending.iter_mut().find(|(p, _)| *p == l) {
                    Some((_, seeds)) => {
                        if !seeds.contains(&(j as u32)) {
                            seeds.push(j as u32);
                        }
                    }
                    None => self.pending.push((l, alloc::vec![j as u32])),
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
        while !self.pending.is_empty() {
            if !spend() {
                break;
            }
            let (l, seeds) = self.pending.remove(0);
            self.search(l, &seeds, &mut out);
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

    /// One past the largest component id this grid has ever issued.
    ///
    /// Component ids are dense and reused, so this bounds them without being a
    /// live count — ask [`component_size`](Self::component_size) whether a given
    /// id is in use. A stitching layer needs it to size its own per-label table
    /// (R-028).
    #[must_use]
    pub fn label_count(&self) -> usize {
        self.size.len()
    }

    /// Air samples carrying this component id; `0` if it is retired or was never
    /// issued.
    #[must_use]
    pub fn component_size(&self, label: u32) -> u32 {
        self.size_of(label)
    }

    /// Air–solid faces bounding this component; `0` if it is retired or was
    /// never issued. Face units — multiply by `h²` for world area. Together
    /// with [`component_size`](Self::component_size) this is what a Sabine
    /// RT60 estimate reads (R-036). Domain-boundary faces count as solid (the
    /// sealed-box convention). Exact in the drained state; under a
    /// budget-truncated op it is conservatively stale exactly as labels are.
    #[must_use]
    pub fn component_area(&self, label: u32) -> u32 {
        self.area.get(label as usize).copied().unwrap_or(0)
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
                self.area.push(0);
                (self.size.len() - 1) as u32
            }
        }
    }

    fn retire(&mut self, l: u32) {
        self.live = self.live.saturating_sub(1);
        if let Some(slot) = self.size.get_mut(l as usize) {
            *slot = 0;
        }
        if let Some(slot) = self.area.get_mut(l as usize) {
            *slot = 0;
        }
        self.free.push(l);
        self.pending.retain(|(p, _)| *p != l);
    }

    fn set_size(&mut self, l: u32, n: u32) {
        if let Some(slot) = self.size.get_mut(l as usize) {
            *slot = n;
        }
    }

    fn size_of(&self, l: u32) -> u32 {
        self.size.get(l as usize).copied().unwrap_or(0)
    }

    fn area_add(&mut self, l: u32, n: u32) {
        if let Some(slot) = self.area.get_mut(l as usize) {
            *slot = slot.saturating_add(n);
        }
    }

    fn area_sub(&mut self, l: u32, n: u32) {
        if let Some(slot) = self.area.get_mut(l as usize) {
            *slot = slot.saturating_sub(n);
        }
    }

    /// Faces of `i` that touch solid or leave the lattice — the sample's
    /// contribution to its component's boundary area.
    fn solid_faces(&self, i: usize) -> u32 {
        let mut nb = [0usize; 6];
        let used = self.neighbours(i, &mut nb);
        let air = nb
            .iter()
            .take(used)
            .filter(|&&j| self.air.get(j) == Some(&true))
            .count();
        6 - air as u32
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
        // Label *and* a member of it, recorded as it is discovered. Searching
        // for a member afterwards would be a scan of the whole lattice per
        // merge, which is the cost this whole module exists to avoid.
        let mut touched: Vec<(u32, usize)> = Vec::new();
        let mut queue = core::mem::take(&mut self.queue);
        queue.clear();
        queue.push(i);
        if let Some(slot) = self.label.get_mut(i) {
            *slot = provisional;
        }
        let mut grown = 0u32;
        let mut blob_faces = 0u32;
        let mut head = 0;
        let mut nb = [0usize; 6];
        while let Some(&at) = queue.get(head) {
            head += 1;
            grown += 1;
            // The blob member's own boundary contribution, read off the
            // finished phase field; faces it shares with established air were
            // counted on the established side while `at` was solid, and stop
            // being boundary now (R-036).
            blob_faces += self.solid_faces(at);
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
                    Some(l) if l != provisional => {
                        self.area_sub(l, 1);
                        if !touched.iter().any(|(t, _)| *t == l) {
                            touched.push((l, j));
                        }
                    }
                    _ => {}
                }
            }
        }
        self.queue = queue;
        self.set_size(provisional, grown);
        self.area_add(provisional, blob_faces);
        repair.relabels += u64::from(grown);

        // Union by size: whichever component is largest keeps its label, and
        // every other one is rewritten into it. Ties go to the smaller id so the
        // outcome does not depend on discovery order (M-36).
        let mut best = provisional;
        for &(l, _) in &touched {
            let (sl, sb) = (self.size_of(l), self.size_of(best));
            if sl > sb || (sl == sb && l < best) {
                best = l;
            }
        }

        let mut merged = 0u32;
        for &(l, member) in &touched {
            if l == best {
                continue;
            }
            let moved = self.recolour(member, l, best);
            repair.relabels += u64::from(moved);
            repair.merges += 1;
            merged += moved;
            let transfer = self.component_area(l);
            self.area_sub(l, transfer);
            self.area_add(best, transfer);
            self.retire(l);
        }
        if best != provisional {
            // The newly-dug blob is itself a component being absorbed, so this
            // is a join and counts as one. Only when the blob is the largest
            // does it keep its label and absorb the others instead.
            let moved = self.recolour(i, provisional, best);
            repair.relabels += u64::from(moved);
            repair.merges += 1;
            merged += moved;
            let transfer = self.component_area(provisional);
            self.area_sub(provisional, transfer);
            self.area_add(best, transfer);
            self.retire(provisional);
        }
        let total = self.size_of(best) + merged;
        self.set_size(best, total);
    }

    /// Lockstep replacement search over the component labelled `l`, starting
    /// from `seeds`.
    ///
    /// Every seed grows a frontier one voxel per round; frontiers that meet
    /// merge. The search stops when at most one frontier is still active, so the
    /// surviving piece is **never walked to completion** and the cost is the
    /// second-largest piece rather than the component (P-26).
    ///
    /// **The seeds are given, not rederived.** Rederiving them as "air adjacent
    /// to solid" collects the entire cave surface rather than this fill's own
    /// neighbourhood — thousands of frontiers instead of a handful, which makes
    /// the walk explore the whole component. That is not a slow version of this
    /// function; it is a different one, and it falsified P-26 on the first run.
    fn search(&mut self, l: u32, seeds: &[u32], out: &mut Fill) {
        if seeds.len() < 2 {
            return;
        }

        self.epoch = self.epoch.wrapping_add(1);
        if self.epoch == 0 {
            self.stamp.iter_mut().for_each(|s| *s = 0);
            self.epoch = 1;
        }
        let epoch = self.epoch;

        // One frontier per seed. `queues` holds what is still to expand and
        // `seen` the members found so far, so an exhausted frontier already
        // knows its piece and nothing has to be scanned for afterwards.
        let mut queues: Vec<Vec<usize>> = Vec::new();
        let mut seen: Vec<Vec<usize>> = Vec::new();
        let mut into: Vec<u32> = Vec::new();
        let mut done: Vec<bool> = Vec::new();
        for (f, &s) in seeds.iter().enumerate() {
            let s = s as usize;
            if self.label.get(s) != Some(&l) || self.stamp.get(s) == Some(&epoch) {
                // Solid now, or already claimed by an earlier seed's frontier.
                queues.push(Vec::new());
                seen.push(Vec::new());
                into.push(f as u32);
                done.push(true);
                continue;
            }
            if let Some(slot) = self.stamp.get_mut(s) {
                *slot = epoch;
            }
            if let Some(slot) = self.claim.get_mut(s) {
                *slot = f as u32;
            }
            queues.push(alloc::vec![s]);
            seen.push(alloc::vec![s]);
            into.push(f as u32);
            done.push(false);
        }

        fn resolve(into: &[u32], mut f: u32) -> u32 {
            while into.get(f as usize).copied().unwrap_or(f) != f {
                f = into.get(f as usize).copied().unwrap_or(f);
            }
            f
        }

        let mut nb = [0usize; 6];
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
                        let g = self.claim.get(j).copied().unwrap_or(NONE);
                        if g == NONE {
                            continue;
                        }
                        let (rf, rg) = (resolve(&into, f), resolve(&into, g));
                        if rf == rg {
                            continue;
                        }
                        // Two frontiers met: their pieces are one piece.
                        let (keep, drop_) = if rf < rg { (rf, rg) } else { (rg, rf) };
                        if let Some(slot) = into.get_mut(drop_ as usize) {
                            *slot = keep;
                        }
                        let moved_q = queues.get(drop_ as usize).cloned().unwrap_or_default();
                        let moved_s = seen.get(drop_ as usize).cloned().unwrap_or_default();
                        if let Some(dst) = queues.get_mut(keep as usize) {
                            dst.extend(moved_q);
                        }
                        if let Some(dst) = seen.get_mut(keep as usize) {
                            dst.extend(moved_s);
                        }
                        if let Some(slot) = queues.get_mut(drop_ as usize) {
                            slot.clear();
                        }
                        if let Some(slot) = seen.get_mut(drop_ as usize) {
                            slot.clear();
                        }
                        let revived = done.get(drop_ as usize) == Some(&true);
                        if revived && let Some(slot) = done.get_mut(keep as usize) {
                            *slot = false;
                        }
                        continue;
                    }
                    if let Some(slot) = self.stamp.get_mut(j) {
                        *slot = epoch;
                    }
                    let root = resolve(&into, f);
                    if let Some(slot) = self.claim.get_mut(j) {
                        *slot = root;
                    }
                    if let Some(q) = queues.get_mut(root as usize) {
                        q.push(j);
                    }
                    if let Some(sset) = seen.get_mut(root as usize) {
                        sset.push(j);
                    }
                }
            }
        }

        // Frontiers that exhausted are complete pieces, and each already holds
        // its own members. Whichever frontier is still active keeps `l`; if all
        // of them finished, the largest does.
        let roots: Vec<u32> = (0..queues.len() as u32)
            .filter(|&f| {
                resolve(&into, f) == f && seen.get(f as usize).is_some_and(|m| !m.is_empty())
            })
            .collect();
        if roots.len() < 2 {
            return;
        }
        let unfinished: Vec<u32> = roots
            .iter()
            .copied()
            .filter(|&f| done.get(f as usize) != Some(&true))
            .collect();

        let keeper = match unfinished.first() {
            Some(&f) if unfinished.len() == 1 => f,
            _ => {
                // All finished: the largest piece keeps the label.
                let mut best = *roots.first().unwrap_or(&0);
                for &f in &roots {
                    let (a, b) = (
                        seen.get(f as usize).map_or(0, Vec::len),
                        seen.get(best as usize).map_or(0, Vec::len),
                    );
                    if a > b {
                        best = f;
                    }
                }
                best
            }
        };

        let mut split = false;
        for &f in &roots {
            if f == keeper {
                continue;
            }
            let Some(members) = seen.get(f as usize).cloned() else {
                continue;
            };
            if members.is_empty() {
                continue;
            }
            let fresh = self.take_label();
            let mut moved_faces = 0u32;
            for &i in &members {
                if let Some(slot) = self.label.get_mut(i) {
                    *slot = fresh;
                }
                moved_faces += self.solid_faces(i);
            }
            self.set_size(fresh, members.len() as u32);
            self.area_sub(l, moved_faces);
            self.area_add(fresh, moved_faces);
            let left = self.size_of(l).saturating_sub(members.len() as u32);
            self.set_size(l, left);
            out.shed += 1;
            split = true;
        }
        if split {
            out.splits += 1;
        }
    }
}

mod world;
pub use world::{AirWorld, Seams};

#[cfg(test)]
mod tests;
