//! **P-120 — array-based union-find for per-chunk labelling, answering `✗26` rather than reopening it.**
//!
//! Ticket: R-120. Pre-registered before this harness existed.
//!
//! ```bash
//! cargo bench --bench experiment_p120
//! ```
//!
//! Writes `docs/experiments/p-120.csv`. **Linux only**, for `experiment_p12`'s
//! reason: C1 is a cost ratio, and `M-281` forbids a nanosecond carrying a
//! verdict on a governed CPU. The instrument is `common::counters::Probe`, which
//! is `perf_event_open`, and off Linux this bench refuses rather than recording a
//! zero it did not measure.
//!
//! # The history, because this row could otherwise look like a reopening
//!
//! `connectivity.rs:29-46` **explicitly refuses** a union-find in favour of flat
//! labels, and says why: it *was* a union-find, and adding `fill` broke it
//! (`✗26`). Parent pointers encode **union history, not spatial adjacency** —
//! `Q → P → A` records that `Q` was merged in via `P` and says nothing about
//! whether `Q` touches `P` in the lattice, so a filled sample can be an
//! articulation point *of the tree* while being nothing of the kind *in the
//! graph*, and re-rooting it severs its descendants from a component they are
//! still genuinely part of.
//!
//! The three union-finds that do live in `crates/isomesh/src/` are private and
//! each justified in prose: `validate.rs:541` `Dsu`, `validate/sealing.rs:835`
//! `UnionFind`, and `connectivity/world.rs:387`, whose module docs state the
//! discipline exactly — *"it is a union-find, and ✗26 is precisely the reason a
//! union-find must not be asked to delete. It is never asked: **every restitch
//! builds it from scratch**, so it only ever unions."*
//!
//! **This row is the per-chunk, rebuilt-from-scratch kind.** It never deletes,
//! it never persists across an edit, it is bench-local, and it therefore
//! **answers `✗26` rather than reopening it**: the structure `✗26` condemned is
//! an *incrementally maintained* one, and C2 is the clause that makes the
//! difference checkable instead of asserted — the final labelling must be
//! independent of the order the unions arrive in, which is `✗26`'s defect
//! measured directly rather than argued away.
//!
//! **Wu, Otoo & Shoshani's 5×–100× is on random binary images, not on voxel
//! islands.** It is an upper bound here and not a prediction. Voxel islands are
//! large, smooth and few; random binary images are the regime that makes a
//! two-pass scan look best, and nothing about this fixture resembles them.
//!
//! # What was missing
//!
//! There was no number anywhere in the repository for what per-chunk labelling
//! costs, in either shape. `M-311`'s union count *"stopped naming a quantity that
//! exists when the structure went flat"* (`connectivity.rs:121-124`), and
//! `Repair::relabels` replaced it — but `relabels` counts label writes on the
//! *incremental* path, not the cost of building a labelling from a chunk. The
//! shipped build is `Air::build` (`connectivity.rs:272-322`): one BFS flood fill
//! per component over a `Vec<u32>` of per-sample labels. Nothing measured it
//! against anything.
//!
//! # The two arms, and why one of them is copied verbatim
//!
//! **`flat`** is `Air::build`'s labelling, mirrored bench-local because
//! `crates/isomesh/src/**` is read-only this phase: `label: vec![NONE; count]`,
//! then `for start in 0..count`, then [`flood`] — `Air::flood`
//! (`connectivity.rs:708-734`) and `Air::neighbours` (`:611-646`) transcribed
//! line for line, including the neighbour push order, which fixes BFS discovery
//! order and therefore which label each component gets.
//!
//! Two things are **excluded from both arms**: `Air::build`'s `set_size`
//! bookkeeping and its second `O(6n)` `solid_faces` pass. Neither is labelling,
//! and both are excluded from the flat arm too — which makes the mirror
//! *cheaper* than the shipped structure, so the exclusion is conservative for C1
//! and cannot manufacture the ratio.
//!
//! A mirror is worth nothing on its own, so **every row builds a real
//! `isomesh::connectivity::Air` on the same chunk values and asserts that its
//! partition and the mirror's are identical**, chunk by chunk, through the public
//! `Air::label_of`. `mirror_matches_shipped` is a column and it is asserted;
//! `ns_per_cell_shipped_air` reports what the whole of `Air::build` costs beside
//! it, so the reader can see exactly what the exclusion was worth.
//!
//! **`union_find`** is Wu, Otoo & Shoshani's array-based structure: a flat
//! `Vec<u32>` of parents with path halving and smaller-root-wins, over
//! *provisional labels* rather than over voxels, driven by the classic two-pass
//! raster scan. Pass one walks the chunk in scan order, unions each air sample
//! with the already-visited air neighbours, and writes a provisional label into
//! the same array the final one lands in. Pass two walks **ascending raster
//! order always** and assigns consecutive final ids on first sight of each root.
//!
//! # Why pass two is what makes C2 an equality rather than a hope
//!
//! Pass two's order is fixed and pass one's is not. The first time a root is seen
//! in ascending raster order is the component's minimum linear index, so the
//! final id of a component is a function of the *partition* alone — not of the
//! scan order, not of the provisional ids, not of the order the unions arrived
//! in. So the final label array is expected to be **bit-identical** across scan
//! orders, which is far stronger than "the partition agrees", and it is what
//! `order_independent` records.
//!
//! Three scan orders are exercised, not two: ascending raster (`x` fastest),
//! **descending** raster — which reverses which neighbour of every pair is the
//! visited one — and `z`-fastest raster, which keeps the backward-neighbour set
//! but allocates every provisional label in a different order and unions them in
//! a different sequence. `distinct_final_states` is the number of distinct
//! label arrays over the three, and C2 wants 1.
//!
//! # The two comparators are shown working, not assumed to work
//!
//! `P-70`'s rule. A comparator that cannot see the defect it is looking for
//! reports a pass for the wrong reason, so both are calibrated against a mutant:
//!
//! - **`control_distinct_final_states`** repeats the three scan orders with pass
//!   two deliberately walking in *pass one's* order instead of raster order —
//!   an order-**sensitive** flattening. It must exceed 1, and it is asserted on
//!   the dug arm. On the base arm it is reported and **not** asserted, because a
//!   chunk with one component has one labelling whatever order it is flattened
//!   in — which is the vacuity control's own argument arriving as a measurement.
//! - **`mutant_partition_mismatches`** runs a union-find with the `−z` unions
//!   dropped and counts the samples whose canonical label then differs from the
//!   flat arm's. It is asserted non-zero on every row: a partition comparator
//!   that cannot see a severed component cannot testify that C3 held.
//!
//! # The dug arm, and the vacuity control
//!
//! *A labeller tested on a connected world is not being tested.* Every reference
//! field's air region is one component on most of these chunks, so the fixture
//! carries a **dug** arm built with `isomesh::brush`: a `BrushStack` over a
//! deterministic 5³ lattice of 125 centres, **two brushes each** —
//! `Brush::add(Sphere { radius: wall })` then
//! `Brush::subtract(Sphere { radius: pocket })` at the same centre. `apply` is
//! `min(field, shape)` and then `max(·, −shape)` (`brush.rs:149-155`), so each
//! centre becomes an air pocket inside a solid wall, and outside the wall the
//! sign of the field is untouched: past the wall radius the `min` can only
//! replace one non-negative value with another and the `max` only one negative
//! value with another.
//!
//! **The pocket carries its own wall, and the first version of this arm did
//! not.** That version subtracted pockets only where the base field was already
//! deep solid, and `thin_plate` has no deep solid at 65³ at all: every candidate
//! was rejected, and the dug arm came out with one component — the registered
//! vacuity control firing on the fixture rather than on the labeller. A dug
//! pocket needs walls, and on a thin field there are none to borrow. Supplying
//! them makes the edit a **pure function of the domain**: the same 125 centres
//! on every field at both resolutions, no sampling-dependent filter, and
//! isolation guaranteed by arithmetic — the wall is 0.035 of the domain edge
//! thick, 2.24 cells at 65³, and a 6-connected step changes a sample's distance
//! to the centre by at most one cell, so a two-cell shell cannot be stepped
//! over.
//!
//! `components` and `components_dug` are **maxima over chunks**, not sums. A sum
//! over 64 chunks exceeds 1 in a perfectly connected world and would make the
//! registered control vacuous by arithmetic; a maximum is 1 exactly when every
//! chunk has one component. `components_dug` is carried onto the base rows too,
//! so the assertion `components_dug > 1` is checkable from any row of the CSV.
//!
//! # What a "cell" is here, and what `ratio` is
//!
//! The shipped structure labels **samples** (`connectivity.rs:226`), and a chunk
//! of `c³` cells is meshed on `(c+1)³` samples (`ChunkLayout::sample_shape`). So
//! `chunk_cells` is [`CHUNK_EDGE_CELLS`]`³` = 32,768, `samples_per_chunk` is 33³,
//! and the registered `ns_per_cell_*` columns are denominated in **cells of the
//! world** — `chunks × chunk_cells` — so that the two registered columns agree
//! with each other. `samples_labelled` and `air_samples` are columns, so nothing
//! is hidden by the choice, and the ratio is the same under either denominator.
//!
//! **`ratio` is the cycle ratio**, `cycles_per_cell_flat /
//! cycles_per_cell_union_find`, and C1 is scored on it. It is *not* the
//! nanosecond ratio: `M-280` and `M-281` say a nanosecond is not a unit on a
//! governed CPU. `ns_ratio` is beside it and `ghz` is on the row, so the two can
//! be reconciled by anyone who wants to.
//!
//! `instruction_ratio` is beside both and is the **deterministic** form
//! (`M-279`): instruction counts reproduce to four figures across runs where a
//! cycle ratio band drifted from 0.984 to 1.035 (`R-105`). It does not carry the
//! verdict here, and the reason is mechanical rather than convenient — the
//! mechanism under test is a **locality** mechanism. A BFS flood fill touches
//! `±1`, `±33` and `±1089` in an order the prefetcher cannot follow; a raster
//! scan touches the same array forwards. An instruction count is blind to that
//! by construction, so scoring C1 on instructions would score a different
//! hypothesis. `l1_misses_per_cell_flat` and `l1_misses_per_cell_union_find` are
//! columns for exactly this reason, and `forms_agree` records whether the two
//! ratios land on the same side of the bar.
//!
//! # C1's bar is close, so the band is a column
//!
//! The measured ratio lands between 1.6× and 2.6× and the registered bar is 2×,
//! which is the worst possible place for it: a cycle ratio is not reproducible
//! to the precision the clause is asking for. So the row does not report the
//! ratio alone. `ratio_rep_min`, `ratio_rep_median` and `ratio_rep_max` are the
//! **nine paired repetitions'** own ratios — each repetition's flat window over
//! its own adjacent union-find window — and `ratio_band_decides` is `true`
//! exactly when that whole band lies on one side of 2×. A row whose band
//! straddles the bar has not decided C1, however its `c1_holds` came out, and
//! the column is there so that fact is a cell of the CSV rather than a caveat
//! someone has to remember. `instruction_ratio` is the form that reproduces:
//! across the whole sweep it moves by less than 0.01 between a field's four
//! rows, while the cycle band on the same rows can be twice as wide as the
//! quantity it is measuring.
//!
//! # The build's missing `popcnt`, and why this row does not pay it
//!
//! This phase established that the repository emits **zero `popcnt`
//! instructions**: there is no `.cargo/config.toml` and no `target-cpu`, so the
//! default `x86-64` baseline is in force and `u64::count_ones()` lowers to the
//! SWAR sequence. Every rank/select and bitmap row in Phase 25 has to price its
//! published figures against that. **This row does not, and the claim is
//! checked rather than assumed:**
//!
//! - `grep -c count_ones crates/isomesh/benches/experiment_p120.rs` → **0**,
//! - `grep -c count_ones crates/isomesh/src/connectivity.rs` → **0** — the
//!   *shipped comparand* does not popcount either,
//! - `objdump -d target/release/deps/experiment_p120-* | grep -c popcnt` →
//!   **0**, and `grep -c 3333333333333333` on the same binary → **0**, so the
//!   SWAR fallback is absent as well. Neither instruction sequence is anywhere
//!   in this bench.
//!
//! `count_ones_per_cell_flat` and `count_ones_per_cell_union_find` are both
//! zero columns for that reason, and `target_feature_popcnt` records the build
//! flag beside them. Wu, Otoo & Shoshani is a scan-and-equivalence-table result
//! rather than a bitmap one: the flat arm reads a `Vec<bool>` and a `Vec<u32>`,
//! the union-find arm reads the same two plus a small parent array, and no
//! clause here inherits a figure that assumed a one-cycle population count. C1
//! would not move under `-C target-cpu=native` for this reason and no other.
//!
//! # SHARE
//!
//! Per-chunk island labelling is **not on the extraction path at all**, so
//! `✗51`'s `1/(1 − share/factor)` bar does not apply and **no extraction speedup
//! is claimed or implied by any clause here**. Each clause's reachable share is a
//! column rather than an argument:
//!
//! - **C1's share is `ns_per_cell_flat`, and it is the whole of it.** The
//!   denominator is the named structure — `Air::build`'s flood fill on the same
//!   chunk — and not a fraction of anything larger. `ns_per_cell_shipped_air` is
//!   the column that says what that structure costs in full, including the
//!   `set_size` and `solid_faces` work both arms drop, so a reader can see that
//!   the comparand was not shaved.
//! - **C2's share is `chunks × orders_tested`**, an enumerated population: every
//!   chunk of every row is labelled under all three scan orders and all three
//!   results are compared. `distinct_final_states` is that comparison's whole
//!   output and `control_distinct_final_states` is the proof it can be non-1.
//! - **C3's share is `air_samples`** — the samples that carry a label at all.
//!   Solid samples are `NONE` in both arms and agree trivially, so they are not
//!   the population, and `air_fraction` is a column so the real one is visible.
//!   `mutant_partition_mismatches` is the proof the comparison can fail.

#![allow(
    clippy::cast_possible_truncation,
    clippy::cast_precision_loss,
    clippy::too_many_lines
)]

mod common;

#[cfg(target_os = "linux")]
mod experiment {
    use std::hint::black_box;
    use std::time::Instant;

    use isomesh::brush::{Brush, BrushStack};
    use isomesh::connectivity::Air;
    use isomesh::fields::{ReferenceField, Sphere};
    use isomesh::marching_cubes::table::is_inside;
    use isomesh::{RuntimeShape3, Sdf};

    use crate::common::counters::{MIN_TIME_RATIO, Probe};
    use crate::common::experiment::Run;

    // ─── the fixture ───────────────────────────────────────────────────────

    /// The registered resolutions, in samples per axis of the world grid.
    const RESOLUTIONS: [u32; 2] = [65, 129];
    /// Cells owned by one chunk per axis, `ChunkLayout`'s own unit. 64 cells at
    /// 65³ and 128 at 129³ divide by it exactly, so no chunk is partial.
    const CHUNK_EDGE_CELLS: u32 = 32;
    /// Measured repetitions per arm. The median is taken **per quantity**.
    const REPS: usize = 9;
    /// Untimed passes before a window, so a batch is steady state.
    const WARMUP: usize = 2;
    /// How long one counter window should last, in nanoseconds.
    const TARGET_BATCH_NS: f64 = 30_000_000.0;
    /// Ceiling on the batch, so a cheap row cannot take minutes.
    const MAX_INNER: usize = 512;
    /// C1's bar, from the registration.
    const SPEEDUP_BAR: f64 = 2.0;

    /// Bubble centres per axis for the dug arm: a 5³ lattice, 125 pockets.
    ///
    /// Five and not four, and the reason is chunk arithmetic rather than taste.
    /// The lattice below puts centres at `0.15 + 0.15k` of the domain edge, and
    /// the chunk walls sit at multiples of `1/chunks_per_axis` — halves at 65³,
    /// quarters at 129³. Five points per axis is the smallest count that puts
    /// **two** centres inside one chunk at both resolutions (`0.30` and `0.45`
    /// both land in the second quarter), so at least one chunk has two air
    /// pockets whatever the field does, and the registered vacuity control
    /// cannot be satisfied by luck.
    const BUBBLE_LATTICE: u32 = 5;
    /// The lattice's lowest centre, as a fraction of the domain edge.
    const BUBBLE_T_LO: f64 = 0.15;
    /// The lattice's highest centre, as a fraction of the domain edge.
    const BUBBLE_T_HI: f64 = 0.75;
    /// Radius of the solid wall a pocket is cut into, as a fraction of the edge.
    ///
    /// `0.06` against a lattice spacing of `0.15`, so no two walls touch and
    /// `0.03` of the edge — 1.9 cells at 65³ — of the field is left untouched
    /// between them.
    const BUBBLE_WALL: f64 = 0.06;
    /// Radius of the air pocket cut inside the wall, as a fraction of the edge.
    ///
    /// 1.6 cells at 65³ and 3.2 at 129³. The nearest sample to any point is
    /// within `√3/2 ≈ 0.87` cells, so a ball this size always contains at least
    /// one sample and a pocket can never be invisible to the labeller. The wall
    /// is `0.035` of the edge thick — 2.24 cells at 65³ — and two cells is the
    /// argument rather than a margin: a 6-connected step changes a sample's
    /// distance to the centre by at most one cell, so a solid shell two cells
    /// thick cannot be stepped over and the pocket is an isolated component **by
    /// construction, on every field**.
    const BUBBLE_RADIUS: f64 = 0.025;

    /// No label: the sample is solid. `connectivity.rs:102`.
    const NONE: u32 = u32::MAX;

    // ─── the shipped flat-label labelling, mirrored ────────────────────────

    /// `Air::neighbours` (`connectivity.rs:611-646`), transcribed.
    ///
    /// The push order is load-bearing and is why this is a transcription rather
    /// than a rewrite: it fixes BFS discovery order, and therefore which label
    /// each component receives.
    fn neighbours(dims: [usize; 3], i: usize, out: &mut [usize; 6]) -> usize {
        let (nx, ny, nz) = (dims[0], dims[1], dims[2]);
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

    /// Reusable buffers for the flat arm, so a pass allocates nothing — which is
    /// what `Air` does too (`connectivity.rs:257-259`).
    #[derive(Default)]
    struct FlatScratch {
        label: Vec<u32>,
        queue: Vec<usize>,
    }

    /// `Air::build`'s labelling: `vec![NONE; count]`, then one BFS flood per
    /// component in ascending sample order.
    ///
    /// Returns the component count. Excludes `set_size` and the `solid_faces`
    /// pass, which are bookkeeping rather than labelling — see the module docs
    /// for why that exclusion is conservative for C1.
    fn flat_label(air: &[bool], dims: [usize; 3], scratch: &mut FlatScratch) -> u32 {
        let count = air.len();
        scratch.label.clear();
        scratch.label.resize(count, NONE);
        let mut next = 0u32;
        let mut nb = [0usize; 6];
        for start in 0..count {
            if !air[start] || scratch.label[start] != NONE {
                continue;
            }
            let l = next;
            next += 1;
            // `Air::flood`, verbatim: an explicit queue walked by `head`, never
            // popped, so it grows to the component and is reused across floods.
            scratch.queue.clear();
            scratch.queue.push(start);
            scratch.label[start] = l;
            let mut head = 0;
            while let Some(&i) = scratch.queue.get(head) {
                head += 1;
                let used = neighbours(dims, i, &mut nb);
                for &j in nb.iter().take(used) {
                    if !air[j] || scratch.label[j] != NONE {
                        continue;
                    }
                    scratch.label[j] = l;
                    scratch.queue.push(j);
                }
            }
        }
        next
    }

    // ─── Wu, Otoo & Shoshani: the array-based union-find ───────────────────

    /// A flat array of parents over **provisional labels**, with path halving.
    ///
    /// Wu, Otoo & Shoshani's structure: no nodes, no pointers, one `Vec<u32>`
    /// indexed by label. Smaller root wins, which makes the resulting forest a
    /// function of the label ids rather than of the machine's allocation
    /// pattern, and keeps `find`'s answers monotone in the ids.
    #[derive(Default)]
    struct Dsu {
        parent: Vec<u32>,
    }

    impl Dsu {
        fn clear(&mut self) {
            self.parent.clear();
        }

        fn make(&mut self) -> u32 {
            let l = self.parent.len() as u32;
            self.parent.push(l);
            l
        }

        #[inline]
        fn find(&mut self, mut x: u32) -> u32 {
            while self.parent[x as usize] != x {
                let grandparent = self.parent[self.parent[x as usize] as usize];
                self.parent[x as usize] = grandparent;
                x = grandparent;
            }
            x
        }

        #[inline]
        fn union(&mut self, a: u32, b: u32) {
            let ra = self.find(a);
            let rb = self.find(b);
            if ra == rb {
                return;
            }
            if ra < rb {
                self.parent[rb as usize] = ra;
            } else {
                self.parent[ra as usize] = rb;
            }
        }

        /// Make every entry point directly at its root, in one ascending sweep.
        ///
        /// Smaller-root-wins and path halving both preserve `parent[i] ≤ i`, so
        /// by the time the sweep reaches `i` its parent is already flattened and
        /// one read finishes the job. This is Wu, Otoo & Shoshani's second-pass
        /// shape and the reason their structure is an *array*: the labelling
        /// pass that follows **reads** the equivalence table and never walks it,
        /// so the tree depth stops being a per-voxel cost.
        fn flatten(&mut self) {
            for i in 0..self.parent.len() {
                self.parent[i] = self.parent[self.parent[i] as usize];
            }
        }
    }

    /// Which way pass one walks the chunk.
    ///
    /// Pass **two** never varies: it is ascending raster order in every case,
    /// which is what makes the final labelling a function of the partition.
    #[derive(Clone, Copy, PartialEq, Eq)]
    enum Order {
        /// `x` fastest, ascending. The canonical raster scan, and the one timed.
        Forward,
        /// Descending index, so the visited neighbour of every pair is the other
        /// one.
        Reverse,
        /// `z` fastest, ascending: the same backward-neighbour set, every
        /// provisional label allocated in a different order.
        ZFast,
    }

    impl Order {
        const ALL: [Self; 3] = [Self::Forward, Self::Reverse, Self::ZFast];

        fn name(self) -> &'static str {
            match self {
                Self::Forward => "forward",
                Self::Reverse => "reverse",
                Self::ZFast => "z_fast",
            }
        }
    }

    /// Reusable buffers for the union-find arm.
    #[derive(Default)]
    struct UfScratch {
        /// Provisional labels in pass one, final labels after pass two — one
        /// array, not two, which is the standard shape and the fair one.
        label: Vec<u32>,
        dsu: Dsu,
        remap: Vec<u32>,
    }

    impl UfScratch {
        /// Size the label array for a chunk without initialising it.
        ///
        /// Pass one writes **every** slot — a provisional label on air, `NONE`
        /// on solid — so a `NONE` fill first would be a second write to every
        /// element and no algorithm would ever perform it. The flat arm's fill
        /// is not the same thing and is kept: there `NONE` is the *unlabelled*
        /// marker the flood scan reads, so `Air::build` genuinely pays it
        /// (`connectivity.rs:284`).
        fn size_for(&mut self, count: usize) {
            if self.label.len() != count {
                self.label.resize(count, NONE);
            }
            self.dsu.clear();
        }
    }

    /// Pass one under [`Order::Forward`]: ascending raster, backward neighbours.
    ///
    /// Split out from the general driver and written straight-line because this
    /// is the arm C1 is measured on, and a `match` on the order inside the hot
    /// loop would measure the `match`.
    fn uf_pass1_forward(air: &[bool], dims: [usize; 3], label: &mut [u32], dsu: &mut Dsu) {
        let (nx, ny, nz) = (dims[0], dims[1], dims[2]);
        let plane = nx * ny;
        let mut i = 0usize;
        for z in 0..nz {
            for y in 0..ny {
                for x in 0..nx {
                    if !air[i] {
                        label[i] = NONE;
                        i += 1;
                        continue;
                    }
                    let a = if x > 0 { label[i - 1] } else { NONE };
                    let b = if y > 0 { label[i - nx] } else { NONE };
                    let c = if z > 0 { label[i - plane] } else { NONE };
                    label[i] = resolve(dsu, a, b, c);
                    i += 1;
                }
            }
        }
    }

    /// Pass one under an arbitrary order, used by the C2 arms.
    fn uf_pass1(air: &[bool], dims: [usize; 3], label: &mut [u32], dsu: &mut Dsu, order: Order) {
        let (nx, ny, nz) = (dims[0], dims[1], dims[2]);
        let plane = nx * ny;
        match order {
            Order::Forward => uf_pass1_forward(air, dims, label, dsu),
            Order::Reverse => {
                for i in (0..air.len()).rev() {
                    if !air[i] {
                        label[i] = NONE;
                        continue;
                    }
                    let x = i % nx;
                    let y = (i / nx) % ny;
                    let z = i / plane;
                    let a = if x + 1 < nx { label[i + 1] } else { NONE };
                    let b = if y + 1 < ny { label[i + nx] } else { NONE };
                    let c = if z + 1 < nz { label[i + plane] } else { NONE };
                    label[i] = resolve(dsu, a, b, c);
                }
            }
            Order::ZFast => {
                for x in 0..nx {
                    for y in 0..ny {
                        for z in 0..nz {
                            let i = x + nx * (y + ny * z);
                            if !air[i] {
                                label[i] = NONE;
                                continue;
                            }
                            let a = if x > 0 { label[i - 1] } else { NONE };
                            let b = if y > 0 { label[i - nx] } else { NONE };
                            let c = if z > 0 { label[i - plane] } else { NONE };
                            label[i] = resolve(dsu, a, b, c);
                        }
                    }
                }
            }
        }
    }

    /// The provisional label of a sample whose visited neighbours carry `a`, `b`
    /// and `c`, unioning whatever they disagree about.
    ///
    /// `NONE` is `u32::MAX`, so the minimum of the three *is* the smallest real
    /// label when there is one and `NONE` when there is not — the comparison is
    /// the emptiness test, and there is no separate branch for it.
    #[inline]
    fn resolve(dsu: &mut Dsu, a: u32, b: u32, c: u32) -> u32 {
        let best = a.min(b).min(c);
        if best == NONE {
            return dsu.make();
        }
        if a != NONE && a != best {
            dsu.union(best, a);
        }
        if b != NONE && b != best {
            dsu.union(best, b);
        }
        if c != NONE && c != best {
            dsu.union(best, c);
        }
        best
    }

    /// Pass two: ascending raster **always**, consecutive ids on first sight.
    ///
    /// The first ascending-raster sight of a root is the component's minimum
    /// linear index, so the id a component receives depends on the partition and
    /// on nothing else. Returns the component count.
    ///
    /// The equivalence table is flattened once before the walk, so this pass is
    /// one array read per air sample rather than a tree walk per air sample.
    fn uf_pass2(label: &mut [u32], dsu: &mut Dsu, remap: &mut Vec<u32>) -> u32 {
        dsu.flatten();
        remap.clear();
        remap.resize(dsu.parent.len(), NONE);
        let mut next = 0u32;
        for slot in label.iter_mut() {
            if *slot == NONE {
                continue;
            }
            let root = dsu.parent[*slot as usize] as usize;
            if remap[root] == NONE {
                remap[root] = next;
                next += 1;
            }
            *slot = remap[root];
        }
        next
    }

    /// Pass two's **order-sensitive mutant**: flattening in pass one's order.
    ///
    /// Not an alternative implementation and never used for a result — this is
    /// `control_distinct_final_states`' instrument, and its whole job is to be
    /// seen disagreeing with itself across scan orders.
    fn uf_pass2_mutant(
        label: &mut [u32],
        dims: [usize; 3],
        dsu: &mut Dsu,
        remap: &mut Vec<u32>,
        order: Order,
    ) {
        dsu.flatten();
        remap.clear();
        remap.resize(dsu.parent.len(), NONE);
        let mut next = 0u32;
        let (nx, ny, nz) = (dims[0], dims[1], dims[2]);
        let mut visit = |slot_index: usize, label: &mut [u32], dsu: &Dsu| {
            if label[slot_index] == NONE {
                return;
            }
            let root = dsu.parent[label[slot_index] as usize] as usize;
            if remap[root] == NONE {
                remap[root] = next;
                next += 1;
            }
            label[slot_index] = remap[root];
        };
        match order {
            Order::Forward => {
                for i in 0..label.len() {
                    visit(i, label, dsu);
                }
            }
            Order::Reverse => {
                for i in (0..label.len()).rev() {
                    visit(i, label, dsu);
                }
            }
            Order::ZFast => {
                for x in 0..nx {
                    for y in 0..ny {
                        for z in 0..nz {
                            visit(x + nx * (y + ny * z), label, dsu);
                        }
                    }
                }
            }
        }
    }

    /// One whole union-find labelling of a chunk under [`Order::Forward`].
    ///
    /// This is the arm C1 is measured on. Returns the component count.
    fn uf_label(air: &[bool], dims: [usize; 3], scratch: &mut UfScratch) -> u32 {
        scratch.size_for(air.len());
        uf_pass1_forward(air, dims, &mut scratch.label, &mut scratch.dsu);
        uf_pass2(&mut scratch.label, &mut scratch.dsu, &mut scratch.remap)
    }

    /// The same, with the `−z` unions dropped: `mutant_partition_mismatches`'
    /// instrument, and never used for a result.
    fn uf_label_mutant(air: &[bool], dims: [usize; 3], scratch: &mut UfScratch) {
        let (nx, ny, nz) = (dims[0], dims[1], dims[2]);
        scratch.size_for(air.len());
        let mut i = 0usize;
        for _z in 0..nz {
            for y in 0..ny {
                for x in 0..nx {
                    if !air[i] {
                        scratch.label[i] = NONE;
                        i += 1;
                        continue;
                    }
                    let a = if x > 0 { scratch.label[i - 1] } else { NONE };
                    let b = if y > 0 { scratch.label[i - nx] } else { NONE };
                    scratch.label[i] = resolve(&mut scratch.dsu, a, b, NONE);
                    i += 1;
                }
            }
        }
        uf_pass2(&mut scratch.label, &mut scratch.dsu, &mut scratch.remap);
    }

    // ─── comparing two labellings as partitions ────────────────────────────

    /// Relabel so ids are consecutive in ascending index of first occurrence.
    ///
    /// Two labellings describe the same partition exactly when their canonical
    /// forms are equal, so this is what turns "the label sets are equal" into an
    /// array comparison rather than a set-of-sets one. A `Vec` indexed by label,
    /// never a map: map iteration order is a determinism hazard (`M-36`).
    fn canonicalise(labels: &[u32], scratch: &mut Vec<u32>, out: &mut Vec<u32>) -> u32 {
        let top = labels
            .iter()
            .copied()
            .filter(|&l| l != NONE)
            .max()
            .map_or(0, |m| m as usize + 1);
        scratch.clear();
        scratch.resize(top, NONE);
        out.clear();
        out.reserve(labels.len());
        let mut next = 0u32;
        for &l in labels {
            if l == NONE {
                out.push(NONE);
                continue;
            }
            if scratch[l as usize] == NONE {
                scratch[l as usize] = next;
                next += 1;
            }
            out.push(scratch[l as usize]);
        }
        next
    }

    /// FNV-1a over a label array, so three scan orders can be compared by one
    /// number each instead of by three stored copies of every chunk.
    fn hash_labels(seed: u64, labels: &[u32]) -> u64 {
        let mut h = seed;
        for &l in labels {
            h ^= u64::from(l);
            h = h.wrapping_mul(0x0000_0100_0000_01b3);
        }
        h
    }

    /// Are these labels `0..k` with every id used?
    fn consecutive(labels: &[u32], k: u32) -> bool {
        let mut seen = vec![false; k as usize];
        for &l in labels {
            if l == NONE {
                continue;
            }
            if l >= k {
                return false;
            }
            seen[l as usize] = true;
        }
        seen.iter().all(|s| *s)
    }

    // ─── the fixture: one field at one resolution, base and dug ────────────

    /// The world grid of one row pair, sampled twice.
    struct World {
        /// Samples per axis.
        n: u32,
        /// Chunks per axis.
        chunks_per_axis: u32,
        /// Cell size, in world units.
        cell: f64,
        /// Minimum corner.
        lo: [f64; 3],
        /// Values of the untouched field.
        base: Vec<f64>,
        /// Values of the field under the dug `BrushStack`.
        dug: Vec<f64>,
        /// Bubbles the clearance filter kept.
        bubbles: usize,
    }

    /// Sample an SDF on a cubic grid of `n` samples per axis.
    fn sample_grid<S: Sdf<Scalar = f64>>(sdf: &S, n: u32, lo: [f64; 3], cell: f64) -> Vec<f64> {
        let n = n as usize;
        let mut out = Vec::with_capacity(n * n * n);
        for z in 0..n {
            for y in 0..n {
                for x in 0..n {
                    out.push(sdf.sample([
                        lo[0] + cell * x as f64,
                        lo[1] + cell * y as f64,
                        lo[2] + cell * z as f64,
                    ]));
                }
            }
        }
        out
    }

    /// The dug arm's edit: 125 air pockets, each with its own solid wall.
    ///
    /// Two brushes per pocket and the order is the whole construction —
    /// `Brush::add` a sphere of [`BUBBLE_WALL`], then `Brush::subtract` one of
    /// [`BUBBLE_RADIUS`] at the same centre. `apply` is `min(field, shape)` then
    /// `max(·, −shape)` (`brush.rs:149-155`), so the result is **air inside the
    /// pocket, solid in the wall, and sign-unchanged outside the wall**: beyond
    /// the wall radius the `min` can only replace a non-negative value with
    /// another non-negative one and the `max` can only replace a negative value
    /// with another negative one.
    ///
    /// The first attempt at this arm subtracted pockets **only where the base
    /// field was already deep solid**, and it is worth recording why that was
    /// wrong rather than merely replaced: `thin_plate` has no deep solid at all
    /// at 65³, every candidate was rejected, and its dug arm came out with one
    /// component — the registered vacuity control firing on the fixture instead
    /// of on the labeller. A pocket needs walls, and on a thin field there are
    /// none to borrow, so the fixture supplies them. Carrying its own wall makes
    /// the edit a pure function of the domain: the same 125 centres on every
    /// field and at both resolutions, no sampling-dependent filter, and
    /// isolation guaranteed by arithmetic rather than checked per field.
    fn dug_brushes(lo: [f64; 3], edge: f64) -> Vec<Brush<Sphere<f64>>> {
        let step = (BUBBLE_T_HI - BUBBLE_T_LO) / f64::from(BUBBLE_LATTICE - 1);
        let mut out = Vec::with_capacity(2 * (BUBBLE_LATTICE as usize).pow(3));
        for ki in 0..BUBBLE_LATTICE {
            for kj in 0..BUBBLE_LATTICE {
                for kk in 0..BUBBLE_LATTICE {
                    let at =
                        |k: u32, axis: usize| lo[axis] + edge * (BUBBLE_T_LO + step * f64::from(k));
                    let center = [at(ki, 0), at(kj, 1), at(kk, 2)];
                    out.push(Brush::add(Sphere {
                        center,
                        radius: BUBBLE_WALL * edge,
                    }));
                    out.push(Brush::subtract(Sphere {
                        center,
                        radius: BUBBLE_RADIUS * edge,
                    }));
                }
            }
        }
        out
    }

    impl World {
        fn build<F>(field: &F, n: u32) -> Self
        where
            F: ReferenceField + Sdf<Scalar = f64>,
        {
            let (_, lo, cell) = crate::common::grid::<f64, F>(field, n);
            let base = sample_grid(field, n, lo, cell);
            let brushes = dug_brushes(lo, cell * f64::from(n - 1));
            // The dug field is the crate's own object, evaluated in full: a
            // `BrushStack` over every brush at every sample. No painted-on
            // shortcut to reconcile with it, and no second definition of what
            // "dug" means.
            let stack = BrushStack {
                base: field,
                brushes: &brushes,
            };
            let dug = sample_grid(&stack, n, lo, cell);
            Self {
                n,
                chunks_per_axis: (n - 1) / CHUNK_EDGE_CELLS,
                cell,
                lo,
                base,
                dug,
                bubbles: brushes.len() / 2,
            }
        }

        /// Samples per axis of one chunk: `cells + 1`, `ChunkLayout`'s own
        /// convention, so adjacent chunks share exactly one sample plane.
        const fn chunk_samples(&self) -> usize {
            CHUNK_EDGE_CELLS as usize + 1
        }

        fn chunk_count(&self) -> usize {
            let c = self.chunks_per_axis as usize;
            c * c * c
        }

        /// The values of one chunk, copied out of the world grid.
        fn chunk_values(&self, arm: &str, chunk: usize, out: &mut Vec<f64>) {
            let source = if arm == "dug" { &self.dug } else { &self.base };
            let c = self.chunks_per_axis as usize;
            let (cx, cy, cz) = (chunk % c, (chunk / c) % c, chunk / (c * c));
            let m = self.chunk_samples();
            let n = self.n as usize;
            let step = CHUNK_EDGE_CELLS as usize;
            out.clear();
            out.reserve(m * m * m);
            for z in 0..m {
                for y in 0..m {
                    let row = (cz * step + z) * n * n + (cy * step + y) * n + cx * step;
                    out.extend_from_slice(&source[row..row + m]);
                }
            }
        }
    }

    // ─── counter windows ───────────────────────────────────────────────────

    /// What one counter window read, per unit of work.
    #[derive(Clone, Copy, Default)]
    struct Counted {
        cycles: f64,
        instructions: f64,
        l1_misses: f64,
        ns: f64,
    }

    impl Counted {
        fn scaled(self, k: f64) -> Self {
            Self {
                cycles: self.cycles * k,
                instructions: self.instructions * k,
                l1_misses: self.l1_misses * k,
                ns: self.ns * k,
            }
        }
    }

    /// One counter window over `inner` repetitions of `body`, divided by `inner`.
    ///
    /// The `perf_event` system calls are all outside the counted region.
    /// Multiplexing is refused rather than scaled: this opens one window at a
    /// time — never a nested pair — so on Zen 3's six general-purpose counters
    /// nothing should be scheduled out, and `MIN_TIME_RATIO` is what says so
    /// (`R-121` paid for the nested-window discovery).
    fn window(probe: &mut Probe, inner: usize, mut body: impl FnMut()) -> Counted {
        probe.reset_and_enable();
        let started = Instant::now();
        for _ in 0..inner {
            body();
        }
        let nanos = started.elapsed().as_nanos() as f64;
        probe.disable();
        let counted = probe.read();
        assert!(
            counted.worst_ratio() >= MIN_TIME_RATIO,
            "a counter ran only {:.1}% of the time it was enabled, so its value is an \
             extrapolation rather than a measurement",
            counted.worst_ratio() * 100.0
        );
        Counted {
            cycles: counted.cycles.count as f64,
            instructions: counted.instructions.count as f64,
            l1_misses: counted.l1d_read_misses.count as f64,
            ns: nanos,
        }
        .scaled(1.0 / inner as f64)
    }

    /// The median of one quantity over the repetitions.
    ///
    /// Taken per quantity rather than per repetition, `experiment_p121`'s rule:
    /// one repetition disturbed by another process on the machine should move
    /// one number, not a row.
    fn median(pick: impl Fn(&Counted) -> f64, reps: &[Counted]) -> f64 {
        let mut values: Vec<f64> = reps.iter().map(&pick).collect();
        values.sort_by(f64::total_cmp);
        values[values.len() / 2]
    }

    fn median_counted(reps: &[Counted]) -> Counted {
        Counted {
            cycles: median(|c| c.cycles, reps),
            instructions: median(|c| c.instructions, reps),
            l1_misses: median(|c| c.l1_misses, reps),
            ns: median(|c| c.ns, reps),
        }
    }

    /// Choose a batch size making one window about [`TARGET_BATCH_NS`] long.
    fn calibrate(mut pass: impl FnMut()) -> usize {
        let started = Instant::now();
        pass();
        let one = started.elapsed().as_nanos() as f64;
        ((TARGET_BATCH_NS / one.max(1.0)).ceil() as usize).clamp(1, MAX_INNER)
    }

    // ─── one row ───────────────────────────────────────────────────────────

    /// Everything one `(field, resolution, arm)` row measured.
    struct Measured {
        arm: &'static str,
        chunks: usize,
        samples_per_chunk: usize,
        air_samples: u64,
        components_max: u32,
        components_total: u64,
        flat: Counted,
        uf: Counted,
        shipped_ns_per_cell: f64,
        /// The lowest and highest of the nine per-repetition cycle ratios.
        ///
        /// `R-105`'s lesson as a column: the same binary's cycle ratio band
        /// drifted from 0.984 to 1.035 across three runs while its instruction
        /// counts held to four figures, so a cycle ratio reported without its
        /// own spread is a number whose reproducibility the reader cannot see.
        ratio_rep_min: f64,
        ratio_rep_max: f64,
        /// The median of the nine **paired** cycle ratios.
        ///
        /// `ratio` is the ratio of the two medians, so a reader can divide
        /// `cycles_per_cell_flat` by `cycles_per_cell_union_find` and get it
        /// back. This is the other estimator — each repetition's own flat window
        /// over its own adjacent union-find window, medianed — which discards a
        /// whole disturbed pair rather than half of one. The two should agree,
        /// and where they do not the row was measured on a busy machine.
        ratio_rep_median: f64,
        labels_consecutive: bool,
        partition_identical: bool,
        order_independent: bool,
        distinct_final_states: usize,
        control_distinct_final_states: usize,
        mutant_partition_mismatches: u64,
        mirror_matches_shipped: bool,
        provisional_labels: u64,
        inner_flat: usize,
        inner_uf: usize,
    }

    /// Label every chunk of one arm, both ways, and check both against the
    /// shipped structure before either is timed.
    fn measure(world: &World, arm: &'static str) -> Measured {
        let m = world.chunk_samples();
        let dims = [m, m, m];
        let chunks = world.chunk_count();

        // ─ the chunks, as the labeller sees them ─
        //
        // Extracted once, so that neither the `Air::build` provenance column nor
        // anything else is charged for a memcpy out of the world grid.
        let mut values: Vec<Vec<f64>> = Vec::with_capacity(chunks);
        let mut air: Vec<Vec<bool>> = Vec::with_capacity(chunks);
        for chunk in 0..chunks {
            let mut one = Vec::new();
            world.chunk_values(arm, chunk, &mut one);
            air.push(one.iter().map(|v| !is_inside(*v)).collect());
            values.push(one);
        }
        let air_samples: u64 = air
            .iter()
            .map(|a| a.iter().filter(|s| **s).count() as u64)
            .sum();
        assert!(
            air_samples > 0,
            "{arm}: no air sample anywhere, so there is nothing to label"
        );

        // ─ correctness, before anything is timed ─
        let mut flat_scratch = FlatScratch::default();
        let mut uf_scratch = UfScratch::default();
        let mut canon_scratch = Vec::new();
        let mut flat_canon = Vec::new();
        let mut uf_canon = Vec::new();
        let mut shipped_labels = Vec::new();
        let mut shipped_canon = Vec::new();

        let mut components_max = 0u32;
        let mut components_total = 0u64;
        let mut labels_consecutive = true;
        let mut partition_identical = true;
        let mut mirror_matches_shipped = true;
        let mut mutant_partition_mismatches = 0u64;
        let mut provisional_labels = 0u64;
        let mut order_hash = [0xcbf2_9ce4_8422_2325u64; 3];
        let mut control_hash = [0xcbf2_9ce4_8422_2325u64; 3];
        let shape = RuntimeShape3::new([m as u32; 3]).expect("chunk shape fits u32");

        for (chunk, a) in air.iter().enumerate() {
            // The shipped flat-label structure, mirrored.
            let flat_components = flat_label(a, dims, &mut flat_scratch);
            canonicalise(&flat_scratch.label, &mut canon_scratch, &mut flat_canon);
            components_max = components_max.max(flat_components);
            components_total += u64::from(flat_components);

            // `M-279`: the mirror's first job is to agree with the structure it
            // mirrors. The real `Air`, on the same values, through public API.
            let (shipped, _) =
                Air::build(&values[chunk], &shape).expect("chunk values fit the shape");
            shipped_labels.clear();
            for z in 0..m as u32 {
                for y in 0..m as u32 {
                    for x in 0..m as u32 {
                        shipped_labels.push(shipped.label_of([x, y, z]).unwrap_or(NONE));
                    }
                }
            }
            canonicalise(&shipped_labels, &mut canon_scratch, &mut shipped_canon);
            mirror_matches_shipped &= shipped_canon == flat_canon;
            mirror_matches_shipped &= shipped.components() == u64::from(flat_components);

            // The union-find, in three scan orders. Pass two never varies.
            for (slot, order) in Order::ALL.iter().enumerate() {
                uf_scratch.size_for(a.len());
                uf_pass1(a, dims, &mut uf_scratch.label, &mut uf_scratch.dsu, *order);
                if *order == Order::Forward {
                    provisional_labels += uf_scratch.dsu.parent.len() as u64;
                }
                let components = uf_pass2(
                    &mut uf_scratch.label,
                    &mut uf_scratch.dsu,
                    &mut uf_scratch.remap,
                );
                order_hash[slot] = hash_labels(order_hash[slot], &uf_scratch.label);
                if *order == Order::Forward {
                    labels_consecutive &= consecutive(&uf_scratch.label, components);
                    canonicalise(&uf_scratch.label, &mut canon_scratch, &mut uf_canon);
                    partition_identical &= uf_canon == flat_canon;
                }

                // The order-sensitive mutant, on the same pass-one result.
                uf_scratch.size_for(a.len());
                uf_pass1(a, dims, &mut uf_scratch.label, &mut uf_scratch.dsu, *order);
                uf_pass2_mutant(
                    &mut uf_scratch.label,
                    dims,
                    &mut uf_scratch.dsu,
                    &mut uf_scratch.remap,
                    *order,
                );
                control_hash[slot] = hash_labels(control_hash[slot], &uf_scratch.label);
            }

            // The severed-union mutant, for `mutant_partition_mismatches`.
            uf_label_mutant(a, dims, &mut uf_scratch);
            canonicalise(&uf_scratch.label, &mut canon_scratch, &mut uf_canon);
            mutant_partition_mismatches += uf_canon
                .iter()
                .zip(&flat_canon)
                .filter(|(x, y)| x != y)
                .count() as u64;
        }

        let distinct_final_states = distinct(&order_hash);
        let control_distinct_final_states = distinct(&control_hash);

        // ─ the timed arms: sibling windows, never nested, and interleaved ─
        //
        // Two windows, one arm each, taken **back to back within a repetition**
        // rather than as nine of one followed by nine of the other. A ratio of
        // two measurements a third of a second apart on a machine with other
        // work on it is a ratio of two different machines; adjacent windows
        // share whatever the clock and the other cores were doing. Nested
        // windows are what is forbidden — Zen 3 has six general-purpose
        // counters and `Probe` opens six, so a nested pair multiplexes and
        // `worst_ratio` refuses (`R-121` paid for that discovery). These are
        // siblings, and `ratio_rep_min` and `ratio_rep_max` report what the
        // ratio did across the nine pairs rather than only its median.
        let mut probe = Probe::open();
        let mut flat_reps = Vec::with_capacity(REPS);
        let mut uf_reps = Vec::with_capacity(REPS);
        let (inner_flat, inner_uf) = {
            let mut flat_pass = || {
                let mut total = 0u32;
                for a in &air {
                    total = total.wrapping_add(flat_label(a, dims, &mut flat_scratch));
                }
                black_box(total);
            };
            let mut uf_pass = || {
                let mut total = 0u32;
                for a in &air {
                    total = total.wrapping_add(uf_label(a, dims, &mut uf_scratch));
                }
                black_box(total);
            };
            // Warm first, then size the batch. A cold first pass mis-sizes
            // `inner` by whatever the page faults cost, and a short window is
            // exactly the window a load spike distorts most.
            for _ in 0..WARMUP {
                flat_pass();
                uf_pass();
            }
            let inner_flat = calibrate(&mut flat_pass);
            let inner_uf = calibrate(&mut uf_pass);
            for _ in 0..REPS {
                flat_reps.push(window(&mut probe, inner_flat, &mut flat_pass));
                uf_reps.push(window(&mut probe, inner_uf, &mut uf_pass));
            }
            (inner_flat, inner_uf)
        };
        let mut rep_ratios: Vec<f64> = flat_reps
            .iter()
            .zip(&uf_reps)
            .map(|(f, u)| f.cycles / u.cycles)
            .collect();
        rep_ratios.sort_by(f64::total_cmp);

        // What the whole of `Air::build` costs — sizes, faces, allocations and
        // all. The column that says the comparand was not shaved by measuring
        // the mirror instead. The chunk's values are already extracted, so this
        // is `Air::build` and nothing else. `Instant` and not counters, because
        // it is provenance rather than a clause: `ghz` is on the row so a reader
        // can convert it, and nothing here is scored on it.
        let shipped_ns = {
            let mut shipped_pass = || {
                for one in &values {
                    let (built, _) = Air::build(one, &shape).expect("chunk values fit the shape");
                    black_box(built.components());
                }
            };
            for _ in 0..WARMUP {
                shipped_pass();
            }
            let inner = calibrate(&mut shipped_pass);
            let started = Instant::now();
            for _ in 0..inner {
                shipped_pass();
            }
            started.elapsed().as_nanos() as f64 / inner as f64
        };

        let cells = (chunks * (CHUNK_EDGE_CELLS as usize).pow(3)) as f64;
        Measured {
            arm,
            chunks,
            samples_per_chunk: m * m * m,
            air_samples,
            components_max,
            components_total,
            flat: median_counted(&flat_reps).scaled(1.0 / cells),
            uf: median_counted(&uf_reps).scaled(1.0 / cells),
            shipped_ns_per_cell: shipped_ns / cells,
            ratio_rep_min: rep_ratios[0],
            ratio_rep_max: rep_ratios[REPS - 1],
            ratio_rep_median: rep_ratios[REPS / 2],
            labels_consecutive,
            partition_identical,
            order_independent: distinct_final_states == 1,
            distinct_final_states,
            control_distinct_final_states,
            mutant_partition_mismatches,
            mirror_matches_shipped,
            provisional_labels,
            inner_flat,
            inner_uf,
        }
    }

    fn distinct(hashes: &[u64; 3]) -> usize {
        let mut seen: Vec<u64> = Vec::with_capacity(3);
        for h in hashes {
            if !seen.contains(h) {
                seen.push(*h);
            }
        }
        seen.len()
    }

    // ─── the run ───────────────────────────────────────────────────────────

    /// Emit one row, scoring the three clauses and asserting the controls.
    fn emit(run: &mut Run, field: &'static str, n: u32, m: &Measured, world: &World, dug: u32) {
        let ratio = m.flat.cycles / m.uf.cycles;
        let ns_ratio = m.flat.ns / m.uf.ns;
        let instruction_ratio = m.flat.instructions / m.uf.instructions;
        let c1 = ratio >= SPEEDUP_BAR;
        let c2 = m.labels_consecutive && m.order_independent;
        let c3 = m.partition_identical;

        // VACUITY CONTROL, from the registration: a labeller tested on a
        // connected world is not being tested.
        assert!(
            dug > 1,
            "{field} {n}: the dug arm has {dug} component(s) in its busiest chunk, so the fixture \
             is a connected world and the row measures nothing"
        );
        assert!(
            m.mirror_matches_shipped,
            "{field} {n} {}: the bench-local flat arm disagrees with `Air::build`, so it is not \
             the shipped structure and C1's denominator is not the registered one",
            m.arm
        );
        assert!(
            m.mutant_partition_mismatches > 0,
            "{field} {n} {}: the severed-union mutant produced the same partition, so the \
             comparator behind C3 cannot see a wrong answer",
            m.arm
        );
        if m.arm == "dug" {
            assert!(
                m.control_distinct_final_states > 1,
                "{field} {n} dug: the order-sensitive mutant flattening produced one state, so \
                 the comparator behind C2 cannot see order dependence"
            );
        }

        run.record(&[
            ("field", field.to_string()),
            ("resolution", n.to_string()),
            ("arm", m.arm.to_string()),
            ("chunk_cells", (CHUNK_EDGE_CELLS as u64).pow(3).to_string()),
            ("components", m.components_max.to_string()),
            ("components_dug", dug.to_string()),
            ("ns_per_cell_flat", format!("{:.5}", m.flat.ns)),
            ("ns_per_cell_union_find", format!("{:.5}", m.uf.ns)),
            ("ratio", format!("{ratio:.4}")),
            ("labels_consecutive", m.labels_consecutive.to_string()),
            ("partition_identical", m.partition_identical.to_string()),
            ("order_independent", m.order_independent.to_string()),
            ("c1_holds", c1.to_string()),
            ("c2_holds", c2.to_string()),
            ("c3_holds", c3.to_string()),
            // ─ extras: the deterministic form, the mechanism, the controls ─
            ("ns_ratio", format!("{ns_ratio:.4}")),
            ("cycles_per_cell_flat", format!("{:.5}", m.flat.cycles)),
            ("cycles_per_cell_union_find", format!("{:.5}", m.uf.cycles)),
            (
                "instructions_per_cell_flat",
                format!("{:.5}", m.flat.instructions),
            ),
            (
                "instructions_per_cell_union_find",
                format!("{:.5}", m.uf.instructions),
            ),
            ("instruction_ratio", format!("{instruction_ratio:.4}")),
            (
                "forms_agree",
                ((instruction_ratio >= SPEEDUP_BAR) == c1).to_string(),
            ),
            // Phase 25's build fact, priced for this row: neither arm calls
            // `count_ones` and the binary contains neither a `popcnt` nor the
            // SWAR fallback. See the module docs for the three checks.
            (
                "target_feature_popcnt",
                cfg!(target_feature = "popcnt").to_string(),
            ),
            ("count_ones_per_cell_flat", "0".to_string()),
            ("count_ones_per_cell_union_find", "0".to_string()),
            (
                "l1_misses_per_cell_flat",
                format!("{:.5}", m.flat.l1_misses),
            ),
            (
                "l1_misses_per_cell_union_find",
                format!("{:.5}", m.uf.l1_misses),
            ),
            (
                "l1_miss_ratio",
                format!("{:.4}", m.flat.l1_misses / m.uf.l1_misses.max(1e-12)),
            ),
            (
                "ns_per_cell_shipped_air",
                format!("{:.5}", m.shipped_ns_per_cell),
            ),
            ("ratio_rep_min", format!("{:.4}", m.ratio_rep_min)),
            ("ratio_rep_max", format!("{:.4}", m.ratio_rep_max)),
            ("ratio_rep_median", format!("{:.4}", m.ratio_rep_median)),
            (
                "ratio_band_decides",
                (m.ratio_rep_min >= SPEEDUP_BAR || m.ratio_rep_max < SPEEDUP_BAR).to_string(),
            ),
            ("ghz", format!("{:.4}", m.flat.cycles / m.flat.ns)),
            ("chunks", m.chunks.to_string()),
            ("samples_per_chunk", m.samples_per_chunk.to_string()),
            (
                "samples_labelled",
                (m.chunks * m.samples_per_chunk).to_string(),
            ),
            ("air_samples", m.air_samples.to_string()),
            (
                "air_fraction",
                format!(
                    "{:.5}",
                    m.air_samples as f64 / (m.chunks * m.samples_per_chunk) as f64
                ),
            ),
            ("components_total", m.components_total.to_string()),
            ("provisional_labels", m.provisional_labels.to_string()),
            ("orders_tested", Order::ALL.len().to_string()),
            ("distinct_final_states", m.distinct_final_states.to_string()),
            (
                "control_distinct_final_states",
                m.control_distinct_final_states.to_string(),
            ),
            (
                "mutant_partition_mismatches",
                m.mutant_partition_mismatches.to_string(),
            ),
            (
                "mirror_matches_shipped",
                m.mirror_matches_shipped.to_string(),
            ),
            ("bubbles", world.bubbles.to_string()),
            ("cell_size", format!("{:.6}", world.cell)),
            ("domain_lo", format!("{:.4}", world.lo[0])),
            ("inner_reps_flat", m.inner_flat.to_string()),
            ("inner_reps_union_find", m.inner_uf.to_string()),
            ("reps", REPS.to_string()),
        ]);

        println!(
            "  {field:<15} {n:>4} {:<5} components {:>4} (dug {:>4})  ratio {ratio:.3}× \
             (instr {instruction_ratio:.3}×, ns {ns_ratio:.3}×)  L1 {:.3}×  \
             C1 {} C2 {} C3 {}",
            m.arm,
            m.components_max,
            dug,
            m.flat.l1_misses / m.uf.l1_misses.max(1e-12),
            if c1 { "✓" } else { "✗" },
            if c2 { "✓" } else { "✗" },
            if c3 { "✓" } else { "✗" },
        );
    }

    pub(crate) fn run(run: &mut Run) {
        let orders: Vec<&str> = Order::ALL.iter().map(|o| o.name()).collect();
        println!(
            "P-120: {}³-cell chunks, scan orders [{}], {REPS} counter windows per arm.\n",
            CHUNK_EDGE_CELLS,
            orders.join(" ")
        );

        for n in RESOLUTIONS {
            isomesh::for_each_reference_field!(f64, |name, field| {
                let world = World::build(&field, n);
                let base = measure(&world, "base");
                let dug = measure(&world, "dug");
                let dug_components = dug.components_max;
                emit(run, name, n, &base, &world, dug_components);
                emit(run, name, n, &dug, &world, dug_components);
            });
        }
    }
}

fn main() {
    if !std::env::args().any(|arg| arg == "--bench") {
        return;
    }

    let prereg = isomesh::experiment!("P-120");

    #[cfg(target_os = "linux")]
    common::experiment::run(prereg, experiment::run);

    // `experiment_p12`'s precedent. C1 is a cost ratio and `M-281` forbids a
    // clock carrying it, so off Linux there is nothing to degrade to: a recorded
    // zero would be a number this harness did not measure.
    #[cfg(not(target_os = "linux"))]
    {
        eprintln!(
            "{} scores C1 on a cycle ratio from hardware performance counters, and this platform \
             has no `perf_event_open`. Refusing rather than recording a column it did not measure.",
            prereg.id
        );
        std::process::exit(1);
    }
}
