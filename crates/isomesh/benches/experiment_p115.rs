//! **P-115 — Tree-Encoded Bitmaps for a subblock-empty summary, because WAH and EWAH are foreclosed.**
//!
//! Ticket: R-115. Pre-registered before this harness existed.
//!
//! ```bash
//! cargo bench --bench experiment_p115
//! ```
//!
//! Writes `docs/experiments/p-115.csv`. **Linux only**, for `experiment_p12`'s
//! reason: C1 is a cost ratio, and on a governed CPU a nanosecond is not a unit
//! (`✗24`, `M-280`, `M-281`). The verdict reads cycles, the instruction form is
//! beside it, and off Linux there is nothing to degrade to, so the harness
//! refuses rather than recording a zero.
//!
//! # What was missing
//!
//! Nothing in this repository had ever compressed the active-cell bitmap. The
//! bitmap itself is built and thrown away inside one function
//! (`dual.rs:359-381`, consumed at `:424-497`), and every downstream row that
//! wants a *summary* of it — `P-107`'s rank directory, `P-114`'s hierarchy,
//! `P-112`'s compaction — assumes the summary is stored flat because there was
//! no measurement saying otherwise.
//!
//! The sweep that produced this phase closed the obvious candidates first.
//! **WAH and EWAH are foreclosed, and this row exists because of the
//! foreclosure rather than in spite of it.** Word-aligned hybrid RLE stores a
//! bitmap as a sequence of fill words and literal words, and the position of
//! bit `k` in that sequence is only knowable by decoding every word before it:
//! random access is O(number of words), not O(1). Every structure this phase
//! wants to build on top of the active-cell bitmap needs `rank(k)` — an O(1)
//! random-access primitive — so an encoding that destroys random access is not
//! a candidate however well it compresses. Tree-Encoded Bitmaps are the only
//! RLE-family scheme that keeps it: the runs are tree *nodes*, the tree is a
//! succinct bit sequence, and a lookup is a path walk with O(1) child indexing.
//! That is the whole reason this is a TEB row and not a WAH row, and it is
//! stated here so a later reader does not re-propose WAH.
//!
//! # The structure, as the paper specifies it
//!
//! Lang, Beischl, Leis, Boncz, Neumann & Kemper, *Tree-Encoded Bitmaps*,
//! SIGMOD 2020, `10.1145/3318464.3380588`, read from the corpus. Implemented
//! faithfully, bench-local, in four parts:
//!
//! 1. **A perfect binary tree over the bitmap's index space.** Bitmap length
//!    `n = 2^h`; leaf `j` of the bottom level carries bit `j` as its *label*.
//!    Both registered resolutions give an exact power of two — 65 samples span
//!    64 cells and `64³ = 2¹⁸`, 129 span 128 and `128³ = 2²¹` — so nothing is
//!    padded and [`Tree::build`] asserts it rather than assuming it.
//! 2. **Bottom-up pruning.** Two sibling leaves with the same label are removed
//!    and the label moves to the parent, so a node at depth `d` that became a
//!    leaf represents a run of `2^(h−d)` equal bits. Pruning is *hereditary*: a
//!    node is prunable exactly when its whole subtree is uniform, which is what
//!    [`Tree::uniform`] is.
//! 3. **The level-order binary marked encoding.** Breadth-first, one bit per
//!    node into `T` — 1 for inner, 0 for leaf — and the leaf labels into `L` in
//!    the same order. The paper's navigation is
//!    `right-child(i) = 2·rank(i)`, `left-child(i) = right-child(i) − 1`,
//!    `label(i) = L[i − rank(i)]`, with `rank` the **inclusive** count of
//!    1-bits in `T[0..=i]`. The paper's worked example — `T = 1100100`,
//!    `L = 0101`, which decodes to `11010000` — is checked by
//!    [`the_papers_worked_example_decodes`].
//! 4. **The rank directory.** A `u32` per 512 bits of `T`, the paper's own
//!    granularity, so a rank is one array load plus at most eight `popcount`s
//!    (four on average) and the walk is O(1) per level. The paper puts the
//!    resulting overhead at 6.25% of the tree bits and that is exactly
//!    `32/512`; `rank_bytes` is a column, so it is checked rather than quoted.
//!
//! Both space optimisations are implemented, because without them the
//! comparison is against a strawman:
//!
//! - **Implicit tree nodes.** The leading 1-bits and trailing 0-bits of `T` are
//!   not stored. A node index below `implicit_inner` is an inner node with
//!   `rank(i) = i + 1`; one past the stored range is a leaf whose rank is the
//!   total. In the paper's worst case this removes the tree encoding entirely
//!   and the TEB degenerates *gracefully* into the plain bitmap plus metadata —
//!   which this implementation reaches, because `prune_depth = h` produces
//!   `tree_bits = 0` and `L` identical to the input.
//! - **The smallest instance, not the fully pruned one.** Pruning bottom-up
//!   level by level gives `h + 1` instances; because the leading run of 1-bits
//!   is free, the *fully* pruned tree is not always the smallest one
//!   (the paper's Figure 6). [`Tree::best_prune_depth`] evaluates the encoded
//!   size of every instance and takes the minimum, and the winner is the
//!   `prune_depth` column.
//!
//! # The flat comparand
//!
//! `flat_bytes` is the same information stored the way the crate stores bit
//! answers today: one bit per **cell**, 64 to a `u64`, over the linear cell
//! index. It is built by mirroring the shipped chain — `build_inside_bits`
//! (`dual.rs:359-381`, one bit per **sample**), `inside_word` (`:385`),
//! `inside_word_shifted` (`:395`), the fused `any & !all` of `active_word`
//! (`:424`) and `cell_mask` (`:445`) — with the crate's own `is_inside` called
//! rather than re-spelled. `crates/isomesh/src/**` is read-only this phase, so
//! the chain is copied, not made `pub`.
//!
//! Two notes the mirror makes visible. The bitmap is indexed by **sample**, so
//! `bit_row = size[0].div_ceil(64)` while the cell row needs
//! `cells_x.div_ceil(64)` — 3 words against 2 at 129 samples, the asymmetry
//! `dual.rs:472-484` documents. And at both registered resolutions `cells_x` is
//! 64 or 128, so `cell_mask` returns `!0` on every word: the mask is mirrored
//! because it is part of the chain, and it is inert on this fixture.
//!
//! # C1 is measured on random access, from a seeded generator
//!
//! Sequential access would answer a different question — a TEB has a run
//! iterator for that, and the paper reports it separately. The registered
//! clause is *random* access, so the order is `ACCESSES` indices drawn from
//! SplitMix64 with a fixed seed (`access_seed`, a column, so the run is
//! reproducible) and masked into range; `cells` is a power of two, so the mask
//! is exact and the draw is unbiased.
//!
//! The two arms are **sibling** counter windows over the same index array,
//! never nested: Zen 3 has six general-purpose counters and `Probe` opens six,
//! so a nested window multiplexes and `Counts::worst_ratio` refuses. A third
//! sibling window walks the index array and does nothing else, so the floor
//! both arms pay for the walk is a column (`cycles_per_access_floor`) rather
//! than an unstated confound. Each arm is batched to its **own**
//! `inner_reps`, because they differ by two orders of magnitude in cost and a
//! shared batch left the flat and floor windows a hundred times shorter than
//! the TEB one — `cycles_per_access_floor` then swung 1.24–5.13 across sixteen
//! rows of identical work. `access_ratio` is scored on the **raw** windows,
//! floor included, which is the form most favourable to the clause: the floor
//! is a common addend, so including it can only make the ratio smaller. A
//! falsification is therefore not an artefact of the harness's own overhead,
//! and no fragile floor-subtracted ratio is reported.
//!
//! `mean_node_visits` is why the ratio comes out where it does, and it is
//! taken from the exhaustive C3 sweep at no extra cost: the walk already keeps
//! the level — `direction` is read from it — so it hands the level back and
//! the timed path pays nothing. It is **not** derived from the node index: a
//! node index in a TEB is a position in the level-order sequence of the
//! *pruned* tree, not a heap index of the perfect tree, so `ilog2(i + 1)` is
//! not the depth. The first version of this column read it that way and
//! under-reported the walk threefold.
//!
//! `access_ratio` is the registered column and is the **cycle** ratio;
//! `ns_access_ratio` is the wall-clock form of the same two windows and no
//! clause consults it. `ghz` is on every row because the row carries `ns`
//! columns, per `M-280`. `instruction_access_ratio` is beside the cycle ratio
//! because instruction counts are deterministic and cycle counts are not
//! (`R-105`); `c1_carrier` names which one carries the verdict, and here it is
//! `cycles`, because C1 is a bound on **cost** and the two forms are not close
//! enough to the bar for the choice to move it.
//!
//! # This build has no `POPCNT`, and a rank structure is made of popcounts
//!
//! `target_feature_popcnt` is a column and it is **false**. There is no
//! `.cargo/config.toml` in this repository and no `target-cpu` anywhere, so
//! the default `x86-64` baseline is in force, and the baseline predates
//! SSE4.2. `u64::count_ones()` therefore lowers to the SWAR sequence:
//! `objdump -d` on this bench's own release binary greps **zero** `popcnt`
//! instructions, while `0x5555555555555555`, `0x3333333333333333`,
//! `0x0f0f0f0f0f0f0f0f` and an SSE2 `psadbw` reduction all appear **inside
//! `Teb::find_leaf`**. The 5900X has the instruction; the build does not use
//! it. `R-108` reports the other half of the same asymmetry independently:
//! `trailing_zeros` lowers to `rep bsf`, which *is* the `TZCNT` encoding, so
//! the incumbent bit-walk gets hardware help on this build and every
//! popcount-based challenger does not.
//!
//! That is a fact about the shipped build rather than about TEB, so the row
//! measures it instead of arguing about it. `popcounts_per_access_teb` and
//! `popcounts_per_access_teb64` are exact structural censuses taken in the
//! exhaustive sweep; `popcounts_per_access_flat` is 0, because a masked load
//! counts nothing.
//!
//! And because a clause must not be falsified by a tuning constant, the paper
//! is followed where it says the rank granularity *"offers a space/time
//! trade-off"* determined **empirically**. A second arm — [`WORD_RANK_WORDS`],
//! one `u32` per tree word — is built from the *same tree, the same labels and
//! the same `prune_depth`*, in the same binary and the same run, so `M-281`'s
//! rule that a comparison lives inside one build is kept. At one word per
//! block the head loop's trip count is a compile-time zero: the whole
//! multi-word popcount, its SSE2 prologue and its unpredictable 0–7 branch all
//! disappear, and a rank becomes one `u32` load plus one masked
//! `count_ones`. It costs space — `teb64_rank_bytes` is eight times
//! `rank_bytes` — and `space_ratio_teb64` says whether C2 would survive it.
//! `access_ratio_teb64` is the resulting cost ratio. **If C1 is falsified on
//! both arms it is not falsified by the granularity.**
//!
//! # SHARE
//!
//! No extraction time is claimed, exactly as registered: C1 is a cost bound,
//! C2 a space bound, C3 an equality, and none of the three is a fraction of a
//! total. Each clause's reachable share is a column:
//!
//! - **C1's share is `access_ratio` against `c1_bar` (2.0)**, with
//!   `access_ratio_teb64` beside it for the other rank granularity and
//!   `instruction_access_ratio` for the deterministic form. The population is
//!   `accesses × inner_reps_<arm>` random lookups, a column per arm.
//! - **C2's share is `space_ratio` against `c2_bar` (1.5)**, with the whole
//!   byte breakdown beside it — `tree_bytes`, `label_bytes`, `rank_bytes`,
//!   `metadata_bytes` — so a reader can see which part of the structure the
//!   ratio came from instead of trusting one number, and `space_ratio_teb64`
//!   for what the faster rank costs.
//! - **C3's share is `answers_equal / cells`, and `cells` is the column.** The
//!   population is every cell, not a sample of them.
//!
//! What the row *decides* is whether a subblock-empty summary can be compressed
//! without giving up the O(1) access `P-107`'s directory is built on. If it can
//! be used at all it would sit inside the offset/compaction work, which `P-121`
//! measured at **0.0073–0.0800 of extraction** depending on field and
//! extractor — so the ceiling on any downstream use is `1/(1 − 0.08)` = 1.087×
//! at the very best, and that ceiling is why no clause here is a speedup.
//!
//! # VACUITY CONTROL
//!
//! `answers_equal == cells` is an equality, and an equality asserted by an
//! instrument that cannot see inequality is worth nothing. So every row also
//! encodes a **truncated** tree — the fully pruned tree with its deepest level
//! cut off, each parent of a deleted pair forced to a leaf carrying its left
//! child's label — and sweeps every cell through it too.
//! `truncated_tree_mismatches` must be non-zero, and it is *guaranteed* to be:
//! in a fully pruned tree no pair of sibling leaves shares a label, so every
//! forced leaf answers wrongly over its right child's whole range. The control
//! is built at `prune_depth = 0` for exactly that reason — on a partially
//! pruned instance the guarantee does not hold, and a control that is only
//! probably non-zero is not a control. Both the equality and the non-zero are
//! `assert!`ed, so a broken instrument cannot produce a CSV.

#![allow(
    clippy::cast_possible_truncation,
    clippy::cast_precision_loss,
    clippy::too_many_arguments,
    clippy::too_many_lines
)]

mod common;

#[cfg(target_os = "linux")]
mod experiment {
    use std::hint::black_box;
    use std::time::Instant;

    use isomesh::marching_cubes::table::is_inside;
    use isomesh::{Real, Sdf};

    use crate::common::counters::{MIN_TIME_RATIO, Probe};

    // ─── the fixture ───────────────────────────────────────────────────────

    /// The registered resolutions, in samples per axis. Both give a cell count
    /// that is an exact power of two, which is what a perfect binary tree over
    /// the index space needs: `64³ = 2¹⁸` and `128³ = 2²¹`.
    const RESOLUTIONS: [u32; 2] = [65, 129];

    /// Distinct random cell indices per arm. 64 KiB of `u32`, streamed
    /// identically by all three windows, so the walk cost is common and is
    /// measured on its own as the floor.
    const ACCESSES: usize = 1 << 14;

    /// The seed the access order is drawn from. A column, so the order is
    /// reproducible from the CSV alone.
    const ACCESS_SEED: u64 = 0x0000_0115_5EED_1115;

    /// Measured repetitions per window. Medians are taken per quantity.
    const REPS: usize = 15;

    /// Untimed passes before anything is counted.
    const WARMUP: usize = 2;

    /// How long one counter window should last, in nanoseconds. Calibrated on
    /// the TEB arm — the slowest — and the same `inner_reps` is then used by
    /// all three windows so each counts exactly the same number of accesses.
    const TARGET_BATCH_NS: f64 = 15_000_000.0;

    /// Ceiling on the batch, so a cheap row cannot take minutes.
    const MAX_INNER: usize = 8192;

    /// C1's bar, from the registration: random access at most 2× the flat
    /// bitmap's.
    const ACCESS_RATIO_BAR: f64 = 2.0;

    /// C2's bar, from the registration: at least 1.5× space saving.
    const SPACE_RATIO_BAR: f64 = 1.5;

    /// The paper's rank granularity, in words: one `u32` per 512 tree bits,
    /// which is the 6.25% overhead it quotes. This is the scored structure.
    const PAPER_RANK_WORDS: usize = 8;

    /// One `u32` per **word** of the tree.
    ///
    /// The paper does not fix 512 — it says the granularity *"offers a
    /// space/time trade-off"* and that 512 was determined *empirically*, on a
    /// machine where `popcount` is one instruction. **It is not one here.**
    /// This build emits no `POPCNT` at all: there is no `.cargo/config.toml`
    /// and no `target-cpu`, so the baseline `x86-64` target is in force,
    /// `cfg!(target_feature = "popcnt")` is false, and `u64::count_ones()`
    /// lowers to the SWAR sequence — `objdump -d … | grep -c popcnt` returns
    /// **0** on this bench's own binary while the `0x3333333333333333` and
    /// `0x0f0f0f0f0f0f0f0f` masks appear inside `Teb::find_leaf`, some of them
    /// hoisted into an SSE2 `psadbw` prologue LLVM emits for a loop that runs
    /// nought to seven times.
    ///
    /// So a second arm re-does the paper's empirical step for this target. At
    /// one word per block the head loop has a compile-time trip count of zero
    /// and disappears: a rank is one `u32` load and one masked popcount.
    /// Following the paper's own instruction is not a deviation from it, and
    /// it is what stops C1's falsification from resting on a tuning constant.
    const WORD_RANK_WORDS: usize = 1;

    /// The TEB's fixed header, counted in `teb_bytes` rather than waved away:
    /// bitmap length, tree height, perfect levels, implicit inner count, tree
    /// bit count, label bit count, explicit one count, and the rank entry
    /// count — eight `u32`s.
    const METADATA_BYTES: usize = 32;

    // ─── the flat active-cell bitmap, mirrored from `dual.rs` ──────────────

    /// `sdf::sample_grid`, which is `pub(crate)`. `src/**` is read-only this
    /// phase, so it is copied rather than made `pub`; the padding to the odd
    /// row stride is part of the contract and is kept.
    fn sample_grid<R: Real, S: Sdf<Scalar = R>>(
        sdf: &S,
        n: u32,
        origin: [R; 3],
        cell_size: R,
        row_stride: usize,
    ) -> Vec<R> {
        let nx = n as usize;
        let pad = row_stride - nx;
        let mut out = Vec::with_capacity(row_stride * nx * nx);
        for z in 0..n {
            for y in 0..n {
                for x in 0..n {
                    out.push(sdf.sample([
                        origin[0] + cell_size * R::from_f64(f64::from(x)),
                        origin[1] + cell_size * R::from_f64(f64::from(y)),
                        origin[2] + cell_size * R::from_f64(f64::from(z)),
                    ]));
                }
                out.resize(out.len() + pad, R::ZERO);
            }
        }
        out
    }

    /// `DualMesher::cell_mask` (`dual.rs:445`), verbatim. Inert at both
    /// registered resolutions — `cells_x` is 64 or 128 — and mirrored anyway,
    /// because leaving a step of the chain out is how a mirror stops being one.
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

    /// One bit per cell, 64 to a `u64`, over the linear cell index.
    struct Flat {
        words: Vec<u64>,
        cells: usize,
    }

    impl Flat {
        #[inline]
        fn get(&self, k: usize) -> bool {
            (self.words[k >> 6] >> (k & 63)) & 1 == 1
        }

        fn bytes(&self) -> usize {
            self.words.len() * 8
        }

        fn ones(&self) -> usize {
            self.words.iter().map(|w| w.count_ones() as usize).sum()
        }

        /// 1-runs, which is the paper's clustering denominator: `f` is the mean
        /// 1-run length and it is what says whether the tree has anything to
        /// prune.
        fn one_runs(&self) -> usize {
            let mut runs = 0usize;
            let mut prev = false;
            for k in 0..self.cells {
                let bit = self.get(k);
                if bit && !prev {
                    runs += 1;
                }
                prev = bit;
            }
            runs
        }
    }

    /// OR a 64-bit chunk into the flat bitmap at an arbitrary bit offset.
    ///
    /// `cell_mask` has already cleared the bits past the end of the row, so a
    /// chunk never spills into the next row's cells.
    #[inline]
    fn deposit(words: &mut [u64], offset: usize, chunk: u64) {
        let wi = offset >> 6;
        let sh = offset & 63;
        words[wi] |= chunk << sh;
        if sh != 0 && wi + 1 < words.len() {
            words[wi + 1] |= chunk >> (64 - sh);
        }
    }

    /// The shipped active-cell chain, mirrored: sample, pack signs one bit per
    /// sample, fold `any & !all` over the four rows, mask, deposit.
    fn active_cell_bitmap<R: Real, S: Sdf<Scalar = R>>(
        field: &S,
        n: u32,
        origin: [R; 3],
        cell_size: R,
    ) -> Flat {
        let sx = n as usize;
        // `DualMesher::row_stride` (`dual.rs:333`): `size[0] | 1`, A-024's odd
        // row.
        let row = sx | 1;
        let values = sample_grid(field, n, origin, cell_size, row);

        // `build_inside_bits` (`dual.rs:359-381`): one bit per **sample**, so
        // the row is `size[0].div_ceil(64)` and not `cells_x.div_ceil(64)`.
        let bit_row = sx.div_ceil(64);
        let rows = sx * sx;
        let mut inside = vec![0u64; bit_row * rows];
        for r in 0..rows {
            let src = row * r;
            let dst = bit_row * r;
            for w in 0..bit_row {
                let base = w * 64;
                let count = (sx - base).min(64);
                let mut word = 0u64;
                for k in 0..count {
                    word |= u64::from(is_inside(values[src + base + k])) << k;
                }
                inside[dst + w] = word;
            }
        }

        // `inside_word` (`:385`) and `inside_word_shifted` (`:395`).
        let word = |w: usize, y: usize, z: usize| inside[bit_row * (y + sx * z) + w];
        let shifted = |w: usize, y: usize, z: usize| {
            let lo = word(w, y, z);
            let hi = if w + 1 < bit_row {
                word(w + 1, y, z)
            } else {
                0
            };
            (lo >> 1) | (hi << 63)
        };

        let cx = sx - 1;
        let cells = cx * cx * cx;
        assert!(
            cells.is_power_of_two(),
            "P-115: a perfect binary tree needs a power-of-two index space and \
             {n} samples give {cells} cells"
        );
        // `dual.rs:484`: the **cell** row, one word shorter than the sample row
        // at 65 and 129 samples.
        let cell_words = cx.div_ceil(64);
        debug_assert!(cell_words <= bit_row);

        let mut words = vec![0u64; cells / 64];
        for z in 0..cx {
            for y in 0..cx {
                let row_base = (y + cx * z) * cx;
                for w in 0..cell_words {
                    // `active_word` (`:424`), the fused four-row fold.
                    let mut any = 0u64;
                    let mut all = !0u64;
                    for dz in 0..2usize {
                        for dy in 0..2usize {
                            let a = word(w, y + dy, z + dz);
                            let b = shifted(w, y + dy, z + dz);
                            any |= a | b;
                            all &= a & b;
                        }
                    }
                    let active = (any & !all) & cell_mask(w, cx);
                    deposit(&mut words, row_base + w * 64, active);
                }
            }
        }

        Flat { words, cells }
    }

    // ─── the tree, before it is encoded ────────────────────────────────────

    /// The perfect binary tree over the bitmap, in heap order.
    ///
    /// Heap order **is** level order for a perfect tree — level `d` occupies
    /// the contiguous index range `[2^d − 1, 2^(d+1) − 2]` and levels ascend —
    /// so the paper's breadth-first encoding is a single ascending sweep with
    /// no queue.
    struct Tree {
        /// Whether the subtree rooted here is all one label. Hereditary
        /// downward, which is the property the presence test relies on.
        uniform: Vec<bool>,
        /// For a leaf, its bit; for an inner node, its **left child's** label,
        /// which is the label a pruned node inherits and the label the
        /// truncated control forces.
        label: Vec<bool>,
        height: u32,
    }

    /// What one level-order encoding costs, before it is materialised.
    #[derive(Clone, Copy)]
    struct Shape {
        nodes: usize,
        leading_ones: usize,
        trailing_zeros: usize,
        labels: usize,
        max_depth: u32,
    }

    impl Shape {
        fn tree_bits(self) -> usize {
            self.nodes - self.leading_ones - self.trailing_zeros
        }

        fn rank_entries(self, block_words: usize) -> usize {
            self.tree_bits().div_ceil(block_words * 64)
        }

        fn bytes(self, block_words: usize) -> usize {
            self.tree_bits().div_ceil(8)
                + self.labels.div_ceil(8)
                + self.rank_entries(block_words) * 4
                + METADATA_BYTES
        }
    }

    impl Tree {
        fn build(flat: &Flat) -> Self {
            let leaves = flat.cells;
            assert!(leaves.is_power_of_two() && leaves >= 64);
            let height = leaves.trailing_zeros();
            let nodes = 2 * leaves - 1;
            let mut uniform = vec![false; nodes];
            let mut label = vec![false; nodes];
            for j in 0..leaves {
                let i = leaves - 1 + j;
                uniform[i] = true;
                label[i] = flat.get(j);
            }
            for i in (0..leaves - 1).rev() {
                let l = 2 * i + 1;
                let r = 2 * i + 2;
                label[i] = label[l];
                uniform[i] = uniform[l] && uniform[r] && label[l] == label[r];
            }
            Self {
                uniform,
                label,
                height,
            }
        }

        /// Is node `i` in the instance pruned at levels `>= prune_depth` and
        /// cut below `cut_at`?
        ///
        /// A node is absent exactly when some strict ancestor became a leaf.
        /// Uniformity is hereditary downward, so if the parent is not uniform
        /// no ancestor is, and the only uniform ancestor that can matter is the
        /// parent at `depth − 1`.
        #[inline]
        fn present(&self, i: usize, depth: u32, prune_depth: u32, cut_at: u32) -> bool {
            depth <= cut_at && (i == 0 || !self.uniform[(i - 1) / 2] || depth <= prune_depth)
        }

        /// Is a present node `i` a leaf of that instance?
        #[inline]
        fn is_leaf(&self, i: usize, depth: u32, prune_depth: u32, cut_at: u32) -> bool {
            depth == self.height || depth == cut_at || (self.uniform[i] && depth >= prune_depth)
        }

        /// The instance's nodes, in level order: `(index, depth, inner, label)`.
        fn walk(
            &self,
            prune_depth: u32,
            cut_at: u32,
            mut visit: impl FnMut(usize, u32, bool, bool),
        ) {
            let deepest = cut_at.min(self.height);
            for depth in 0..=deepest {
                let lo = (1usize << depth) - 1;
                let hi = (1usize << (depth + 1)) - 1;
                for i in lo..hi {
                    if !self.present(i, depth, prune_depth, cut_at) {
                        continue;
                    }
                    let leaf = self.is_leaf(i, depth, prune_depth, cut_at);
                    visit(i, depth, !leaf, self.label[i]);
                }
            }
        }

        fn shape(&self, prune_depth: u32, cut_at: u32) -> Shape {
            let mut nodes = 0usize;
            let mut leading_ones = 0usize;
            let mut trailing_zeros = 0usize;
            let mut labels = 0usize;
            let mut max_depth = 0u32;
            let mut seen_leaf = false;
            self.walk(prune_depth, cut_at, |_, depth, inner, _| {
                nodes += 1;
                max_depth = max_depth.max(depth);
                if inner {
                    if !seen_leaf {
                        leading_ones += 1;
                    }
                    trailing_zeros = 0;
                } else {
                    seen_leaf = true;
                    trailing_zeros += 1;
                    labels += 1;
                }
            });
            Shape {
                nodes,
                leading_ones,
                trailing_zeros,
                labels,
                max_depth,
            }
        }

        /// The paper's Figure 6: the fully pruned tree is not always the
        /// smallest TEB, because the leading run of 1-bits is free and pruning
        /// can shorten it. So every bottom-up instance is priced and the
        /// cheapest wins.
        ///
        /// Priced at the **paper's** granularity, because that is the scored
        /// structure; the word-granularity arm is built at the same
        /// `prune_depth` so the two differ in exactly one thing.
        fn best_prune_depth(&self) -> (u32, Shape, Shape) {
            let mut best = (0u32, self.shape(0, self.height));
            let fully_pruned = best.1;
            for d in 1..=self.height {
                let shape = self.shape(d, self.height);
                if shape.bytes(PAPER_RANK_WORDS) < best.1.bytes(PAPER_RANK_WORDS) {
                    best = (d, shape);
                }
            }
            (best.0, best.1, fully_pruned)
        }

        fn encode<const BW: usize>(&self, prune_depth: u32, cut_at: u32) -> Teb<BW> {
            let shape = self.shape(prune_depth, cut_at);
            let tree_bits = shape.tree_bits();
            let mut tree = vec![0u64; tree_bits.div_ceil(64).max(1)];
            let mut labels = vec![0u64; shape.labels.div_ceil(64).max(1)];
            let stop = shape.nodes - shape.trailing_zeros;
            let mut pos = 0usize;
            let mut lpos = 0usize;
            self.walk(prune_depth, cut_at, |_, _, inner, label| {
                if !inner {
                    if label {
                        labels[lpos >> 6] |= 1u64 << (lpos & 63);
                    }
                    lpos += 1;
                }
                if inner && pos >= shape.leading_ones && pos < stop {
                    let j = pos - shape.leading_ones;
                    tree[j >> 6] |= 1u64 << (j & 63);
                }
                pos += 1;
            });
            debug_assert_eq!(lpos, shape.labels);

            let rank_entries = shape.rank_entries(BW);
            let mut rank = vec![0u32; rank_entries.max(1)];
            let mut acc = 0u32;
            for (b, slot) in rank.iter_mut().enumerate().take(rank_entries) {
                *slot = acc;
                let lo = b * BW;
                let hi = ((b + 1) * BW).min(tree.len());
                for w in &tree[lo..hi] {
                    acc += w.count_ones();
                }
            }

            // The point lookup starts at the last level that is entirely inner,
            // which is what the leading run of implicit 1-bits certifies.
            let mut perfect_levels = 0u32;
            while perfect_levels < cut_at.min(self.height)
                && (1usize << (perfect_levels + 1)) - 1 <= shape.leading_ones
            {
                perfect_levels += 1;
            }

            Teb {
                tree,
                tree_bits,
                rank,
                rank_entries,
                labels,
                label_bits: shape.labels,
                implicit_inner: shape.leading_ones,
                explicit_ones: acc as usize,
                perfect_levels,
                height: self.height,
                nodes: shape.nodes,
            }
        }
    }

    // ─── the encoded structure ─────────────────────────────────────────────

    /// A tree-encoded bitmap: `T`, `L`, and the rank directory over `T`.
    ///
    /// `BW` is the rank granularity in 64-bit words. The paper's is 8;
    /// [`WORD_RANK_WORDS`] is the other arm, and at `BW == 1` the head loop's
    /// trip count is a compile-time zero and vanishes.
    struct Teb<const BW: usize> {
        tree: Vec<u64>,
        tree_bits: usize,
        rank: Vec<u32>,
        rank_entries: usize,
        labels: Vec<u64>,
        label_bits: usize,
        implicit_inner: usize,
        explicit_ones: usize,
        perfect_levels: u32,
        height: u32,
        nodes: usize,
    }

    /// One completed point lookup.
    struct Walk {
        /// The leaf's position in `T`.
        leaf: usize,
        /// Its inclusive rank, which is also its label offset base.
        rank: usize,
        /// The tree level it was found at.
        level: u32,
        /// `count_ones` calls the walk made. Dead in [`Teb::get`], which is
        /// what the timed arms call, so it costs the measurement nothing; live
        /// in the exhaustive sweep, which is untimed.
        popcounts: u32,
    }

    impl<const BW: usize> Teb<BW> {
        /// Whether node `i` is inner, its inclusive rank, and the number of
        /// `count_ones` calls it took.
        ///
        /// Inner-ness and rank are wanted together at every level — the rank
        /// picks the children of an inner node and the label of a leaf — and
        /// computing them separately doubles the popcounts.
        ///
        /// The block's words are taken as **one slice**, so the range is
        /// bounds-checked once rather than once per `count_ones`. That is
        /// worth naming because C1 is a cost ratio and the clause deserves a
        /// competent implementation to be falsified against: indexing
        /// `self.tree[w]` inside the loop cost about a third of the walk.
        /// `unsafe_code = "forbid"` is a workspace lint, so one bounds check
        /// per rank is the floor a safe implementation can reach.
        #[inline]
        fn node(&self, i: usize) -> (bool, usize, u32) {
            if i < self.implicit_inner {
                return (true, i + 1, 0);
            }
            let j = i - self.implicit_inner;
            if j >= self.tree_bits {
                return (false, self.implicit_inner + self.explicit_ones, 0);
            }
            let block = j / (BW * 64);
            let bit = j & 63;
            let window = &self.tree[block * BW..=(j >> 6)];
            let (word, head) = window.split_last().expect("the window holds `j`");
            let mut r = self.rank[block] as usize;
            for w in head {
                r += w.count_ones() as usize;
            }
            // `!0 >> (63 - bit)` is the low `bit + 1` bits; `1 << 64` is not.
            r += (*word & (!0u64 >> (63 - bit))).count_ones() as usize;
            (
                (*word >> bit) & 1 == 1,
                self.implicit_inner + r,
                head.len() as u32 + 1,
            )
        }

        #[inline]
        fn label_bit(&self, idx: usize) -> bool {
            (self.labels[idx >> 6] >> (idx & 63)) & 1 == 1
        }

        /// The paper's Algorithm 1: start at the last perfect level, then walk
        /// down taking the bit of `k` that names the direction.
        ///
        /// Returns the leaf's position in `T`, its inclusive rank, and the tree
        /// **level** it was found at. The level is already a live variable —
        /// `direction` needs it — so handing it back costs the timed path
        /// nothing, and it is the only honest way to report the path length: a
        /// node index here is a *position in the level-order sequence of the
        /// pruned tree*, not a heap index of the perfect tree, so
        /// `ilog2(i + 1)` is not the depth. It was, briefly, and it
        /// under-reported the walk by a factor of three.
        #[inline]
        fn find_leaf(&self, k: usize) -> Walk {
            let mut level = self.perfect_levels;
            let mut i = ((1usize << level) - 1) + (k >> (self.height - level));
            let mut popcounts = 0u32;
            loop {
                let (inner, rank, counted) = self.node(i);
                popcounts += counted;
                if !inner {
                    return Walk {
                        leaf: i,
                        rank,
                        level,
                        popcounts,
                    };
                }
                let direction = (k >> (self.height - level - 1)) & 1;
                // right-child(i) = 2·rank(i); left-child(i) = right-child − 1.
                i = 2 * rank - 1 + direction;
                level += 1;
            }
        }

        #[inline]
        fn get(&self, k: usize) -> bool {
            let walk = self.find_leaf(k);
            self.label_bit(walk.leaf - walk.rank)
        }

        fn tree_bytes(&self) -> usize {
            self.tree_bits.div_ceil(8)
        }

        fn label_bytes(&self) -> usize {
            self.label_bits.div_ceil(8)
        }

        fn rank_bytes(&self) -> usize {
            self.rank_entries * 4
        }

        fn bytes(&self) -> usize {
            self.tree_bytes() + self.label_bytes() + self.rank_bytes() + METADATA_BYTES
        }
    }

    /// The paper's worked example, decoded by this implementation.
    ///
    /// Figure 3b encodes as `T = 1100100`, `L = 0101`, and Figure 6 names the
    /// bitmap it came from: `11010000`. Anything that decodes those seven bits
    /// to something else has the child formula, the rank convention or the
    /// label offset wrong, and every number in this file would then be a
    /// measurement of a different data structure.
    fn the_papers_worked_example_decodes() {
        let bits = [true, true, false, true, false, false, false, false];
        // Eight bits is below `Flat`'s 64-bit floor, so the tree is built from
        // the bits directly rather than through `Tree::build`.
        let leaves = 8usize;
        let nodes = 2 * leaves - 1;
        let mut uniform = vec![false; nodes];
        let mut label = vec![false; nodes];
        for (j, bit) in bits.iter().enumerate() {
            uniform[leaves - 1 + j] = true;
            label[leaves - 1 + j] = *bit;
        }
        for i in (0..leaves - 1).rev() {
            label[i] = label[2 * i + 1];
            uniform[i] =
                uniform[2 * i + 1] && uniform[2 * i + 2] && label[2 * i + 1] == label[2 * i + 2];
        }
        let tree = Tree {
            uniform,
            label,
            height: 3,
        };

        // The unoptimised encoding is the paper's: seven nodes, `1100100`, and
        // four labels, `0101`.
        let shape = tree.shape(0, 3);
        assert_eq!(shape.nodes, 7, "P-115: the paper's tree has seven nodes");
        assert_eq!(shape.labels, 4, "P-115: the paper's tree has four labels");
        assert_eq!(
            shape.leading_ones, 2,
            "P-115: `1100100` opens with two 1-bits"
        );
        assert_eq!(
            shape.trailing_zeros, 2,
            "P-115: `1100100` closes with two 0-bits"
        );

        // Both granularities decode the same bitmap, or the rank arm and the
        // scored arm are not the same structure.
        let teb = tree.encode::<PAPER_RANK_WORDS>(0, 3);
        let teb64 = tree.encode::<WORD_RANK_WORDS>(0, 3);
        for (k, bit) in bits.iter().enumerate() {
            assert_eq!(
                teb.get(k),
                *bit,
                "P-115: the paper's worked example decodes wrongly at bit {k}"
            );
            assert_eq!(
                teb64.get(k),
                *bit,
                "P-115: the word-granularity rank decodes wrongly at bit {k}"
            );
        }
    }

    // ─── the random access order ───────────────────────────────────────────

    #[inline]
    fn splitmix(state: &mut u64) -> u64 {
        *state = state.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = *state;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }

    /// `ACCESSES` uniform random cell indices. `cells` is a power of two, so
    /// the mask is exact and there is no modulo bias to explain away.
    fn access_order(cells: usize, seed: u64) -> Vec<u32> {
        let mask = (cells - 1) as u64;
        let mut state = seed;
        (0..ACCESSES)
            .map(|_| (splitmix(&mut state) & mask) as u32)
            .collect()
    }

    // ─── measurement ───────────────────────────────────────────────────────

    #[derive(Clone, Copy)]
    struct Window {
        cycles: f64,
        instructions: f64,
        ns: f64,
    }

    /// One sibling counter window over `inner` passes of `body`.
    ///
    /// Never nested with another: `Probe` opens all six of Zen 3's
    /// general-purpose counters, so a second window inside this one would
    /// multiplex and `worst_ratio` would refuse. That refusal is asserted here
    /// rather than hoped for.
    fn window(probe: &mut Probe, inner: usize, mut body: impl FnMut() -> u64) -> Window {
        probe.reset_and_enable();
        let start = Instant::now();
        let mut acc = 0u64;
        for _ in 0..inner {
            acc = acc.wrapping_add(body());
        }
        let ns = start.elapsed().as_secs_f64() * 1e9;
        probe.disable();
        black_box(acc);
        let counts = probe.read();
        assert!(
            counts.worst_ratio() >= MIN_TIME_RATIO,
            "P-115: a counter was multiplexed ({:.4}), so its value is an \
             extrapolation and not a measurement",
            counts.worst_ratio()
        );
        Window {
            cycles: counts.cycles.count as f64,
            instructions: counts.instructions.count as f64,
            ns,
        }
    }

    fn median(mut v: Vec<f64>) -> f64 {
        v.sort_by(f64::total_cmp);
        v[v.len() / 2]
    }

    /// Every window normalised to one access.
    fn per_access(windows: &[Window], accesses: f64) -> Vec<Window> {
        windows
            .iter()
            .map(|w| Window {
                cycles: w.cycles / accesses,
                instructions: w.instructions / accesses,
                ns: w.ns / accesses,
            })
            .collect()
    }

    /// Per-quantity medians, which is `P-121`'s discipline: one repetition
    /// disturbed by another process on the machine should move one quantity,
    /// not all three.
    fn medians(windows: &[Window]) -> Window {
        Window {
            cycles: median(windows.iter().map(|w| w.cycles).collect()),
            instructions: median(windows.iter().map(|w| w.instructions).collect()),
            ns: median(windows.iter().map(|w| w.ns).collect()),
        }
    }

    /// The median of the **per-repetition** ratios, not the ratio of the
    /// medians.
    ///
    /// The two arms are measured inside the same repetition, so a disturbance
    /// that hit both cancels in the pairing and does not in the division of two
    /// separately-medianed numbers. It matters here: with unpaired ratios the
    /// same binary and the same seed read 106.7× and 63.6× for `sphere` at 65³
    /// on two consecutive runs, while its **instruction** ratio held at 48.4 to
    /// four figures both times — `R-105`'s observation, on this row's own data.
    fn paired(num: &[Window], den: &[Window], of: impl Fn(&Window) -> f64) -> f64 {
        median(
            num.iter()
                .zip(den.iter())
                .map(|(a, b)| of(a) / of(b))
                .collect(),
        )
    }

    /// How many passes of `body` make one window about [`TARGET_BATCH_NS`].
    ///
    /// **Per arm, not shared.** The first version of this calibrated on the
    /// slowest arm and reused the count, which left the flat and floor windows
    /// a hundred times shorter than the TEB one: `cycles_per_access_floor` then
    /// swung 1.24–5.13 across sixteen rows of *identical* work, and the
    /// floor-subtracted ratio came out **negative** on five of them. A rate is
    /// a rate whatever the batch, so each arm gets a window long enough to
    /// measure one.
    fn calibrate(mut body: impl FnMut() -> u64) -> usize {
        let start = Instant::now();
        black_box(body());
        let one = start.elapsed().as_secs_f64() * 1e9;
        ((TARGET_BATCH_NS / one.max(1.0)) as usize).clamp(1, MAX_INNER)
    }

    // ─── one row ───────────────────────────────────────────────────────────

    struct Row {
        field: &'static str,
        resolution: u32,
        cells: usize,
        active_cells: usize,
        one_runs: usize,
        flat_bytes: usize,
        teb_bytes: usize,
        teb_fully_pruned_bytes: usize,
        teb_unoptimised_bytes: usize,
        teb64_bytes: usize,
        teb64_rank_bytes: usize,
        /// `count_ones` calls per lookup, exact and structural, from the
        /// exhaustive sweep. The build emits no `POPCNT`, so each of these is
        /// a SWAR sequence, and this column is what lets a reader price that.
        popcounts_per_access: f64,
        popcounts_per_access_teb64: f64,
        tree_bytes: usize,
        label_bytes: usize,
        rank_bytes: usize,
        tree_nodes: usize,
        perfect_tree_nodes: usize,
        prune_depth: u32,
        perfect_levels: u32,
        tree_height: u32,
        truncated_cut_depth: u32,
        answers_equal: usize,
        truncated_tree_mismatches: usize,
        /// Mean `node()` visits per lookup, over the exhaustive C3 sweep. The
        /// cost driver behind `cycles_per_access_teb`, and the reason a path
        /// walk cannot approach a single masked load.
        mean_node_visits: f64,
        max_node_visits: u32,
        flat: Window,
        teb: Window,
        teb64: Window,
        floor: Window,
        /// The paired per-repetition medians: TEB over flat, one per quantity.
        ratio: Window,
        /// The same, for the word-granularity rank arm.
        ratio64: Window,
        inner_reps_flat: usize,
        inner_reps_teb: usize,
        inner_reps_teb64: usize,
        inner_reps_floor: usize,
    }

    impl Row {
        fn space_ratio(&self) -> f64 {
            self.flat_bytes as f64 / self.teb_bytes as f64
        }

        /// C2 on the word-granularity arm, whose rank directory is eight times
        /// the paper's. Reported so the space price of the access arm is
        /// visible rather than implied.
        fn space_ratio64(&self) -> f64 {
            self.flat_bytes as f64 / self.teb64_bytes as f64
        }

        fn density(&self) -> f64 {
            self.active_cells as f64 / self.cells as f64
        }

        /// The paper's clustering factor `f`: the mean 1-run length.
        fn clustering(&self) -> f64 {
            if self.one_runs == 0 {
                0.0
            } else {
                self.active_cells as f64 / self.one_runs as f64
            }
        }

        /// The fraction of the perfect tree's nodes that pruning removed. This
        /// is the number that answers "does this regime give the tree anything
        /// to prune".
        fn prune_fraction(&self) -> f64 {
            1.0 - self.tree_nodes as f64 / self.perfect_tree_nodes as f64
        }

        fn cycle_ratio(&self) -> f64 {
            self.ratio.cycles
        }

        fn instruction_ratio(&self) -> f64 {
            self.ratio.instructions
        }

        /// The wall-clock form of the same two windows. No clause consults it
        /// (`M-281`); it is on the row so a reader can see that the clock and
        /// the counter agree.
        fn ns_ratio(&self) -> f64 {
            self.ratio.ns
        }

        fn ghz(&self) -> f64 {
            self.teb.cycles / self.teb.ns
        }

        /// C1 as registered: a bound on **cost**, so the cycle form carries it.
        fn c1(&self) -> bool {
            self.cycle_ratio() <= ACCESS_RATIO_BAR
        }

        /// The same clause on the deterministic form. It exists because
        /// `R-105` measured a cycle ratio band drift from 0.984 to 1.035 on an
        /// identical binary while its instruction counts held to four figures,
        /// and a reader is entitled to see that the two forms agree here rather
        /// than take it on trust.
        fn c1_instructions(&self) -> bool {
            self.instruction_ratio() <= ACCESS_RATIO_BAR
        }

        fn c2(&self) -> bool {
            self.space_ratio() >= SPACE_RATIO_BAR
        }

        fn c3(&self) -> bool {
            self.answers_equal == self.cells
        }
    }

    /// One row.
    ///
    /// **Deliberately not generic over the field.** It was, and the timing
    /// closures were then monomorphised once per field type: eight separate
    /// compilations of the same `walk_flat`, each with its own inlining
    /// decisions, and `instructions_per_access_flat` came out 10.0 on `torus`
    /// and 11.5 on `box_exact` for byte-identical work. The field type is only
    /// needed to *sample*, so sampling happens in [`sweep`] and every row is
    /// measured by one compilation of one function against a plain [`Flat`].
    fn measure(probe: &mut Probe, field_name: &'static str, n: u32, flat: &Flat) -> Row {
        let cells = flat.cells;
        let active_cells = flat.ones();
        // A bitmap with no cell of one kind has a one-node tree, nothing to
        // truncate, and no control. Both would be fixture defects, not results.
        assert!(
            active_cells > 0 && active_cells < cells,
            "P-115: {field_name} at {n}³ has {active_cells} active cells of \
             {cells}, so the bitmap is uniform and the control cannot fire"
        );

        let tree = Tree::build(flat);
        let (prune_depth, shape, fully_pruned) = tree.best_prune_depth();
        let teb = tree.encode::<PAPER_RANK_WORDS>(prune_depth, tree.height);
        // The same tree and the same labels, differing in one thing: the rank
        // granularity. Same binary, same run, same fixture — `M-281`'s rule
        // that a comparison lives inside one build is why this is a second
        // structure rather than a second `RUSTFLAGS`.
        let teb64 = tree.encode::<WORD_RANK_WORDS>(prune_depth, tree.height);
        assert_eq!(teb.nodes, shape.nodes);

        // The vacuity control: the **fully pruned** tree with its deepest level
        // removed. Fully pruned is the instance in which no pair of sibling
        // leaves shares a label, so every forced leaf is guaranteed to answer
        // wrongly over its right child's range — which is what makes this a
        // control rather than a hope.
        let cut = fully_pruned.max_depth - 1;
        let truncated = tree.encode::<PAPER_RANK_WORDS>(0, cut);

        // C3, over every cell rather than a sample of them, and the control
        // over the same population in the same sweep. The walk hands back the
        // level it stopped at and the `count_ones` calls it made, so the mean
        // path length and the popcount census fall out of the sweep that was
        // already exhaustive. Both are dead in `Teb::get`, which is what the
        // timed arms call, so neither costs the measurement anything.
        let mut answers_equal = 0usize;
        let mut truncated_tree_mismatches = 0usize;
        let mut visits_total = 0u64;
        let mut max_node_visits = 0u32;
        let mut popcounts_total = 0u64;
        let mut popcounts64_total = 0u64;
        for k in 0..cells {
            let want = flat.get(k);
            let walk = teb.find_leaf(k);
            if teb.label_bit(walk.leaf - walk.rank) == want {
                answers_equal += 1;
            }
            let visits = walk.level - teb.perfect_levels + 1;
            visits_total += u64::from(visits);
            max_node_visits = max_node_visits.max(visits);
            popcounts_total += u64::from(walk.popcounts);
            popcounts64_total += u64::from(teb64.find_leaf(k).popcounts);
            if truncated.get(k) != want {
                truncated_tree_mismatches += 1;
            }
        }
        assert_eq!(
            teb64.get(cells / 3),
            flat.get(cells / 3),
            "P-115: the two rank granularities are not the same structure"
        );
        assert_eq!(
            answers_equal,
            cells,
            "P-115: {field_name} at {n}³ decoded {} cells wrongly",
            cells - answers_equal
        );
        assert!(
            truncated_tree_mismatches > 0,
            "P-115: {field_name} at {n}³ truncated its tree by a level and every \
             cell still answered the same, so the comparator cannot see \
             inequality and C3's equality is vacuous"
        );

        // C1. Four sibling windows over one index array: the array walk alone,
        // the flat bit test, the TEB path walk at the paper's rank
        // granularity, and the same walk at one `u32` per word.
        let order = access_order(cells, ACCESS_SEED);
        let flat_ref = &flat;
        let teb_ref = &teb;
        let teb64_ref = &teb64;
        let order_ref = &order;
        let mut walk_floor = || {
            order_ref
                .iter()
                .fold(0u64, |a, &k| a.wrapping_add(u64::from(black_box(k))))
        };
        let mut walk_flat = || {
            order_ref.iter().fold(0u64, |a, &k| {
                a.wrapping_add(u64::from(flat_ref.get(k as usize)))
            })
        };
        let mut walk_teb = || {
            order_ref.iter().fold(0u64, |a, &k| {
                a.wrapping_add(u64::from(teb_ref.get(k as usize)))
            })
        };
        let mut walk_teb64 = || {
            order_ref.iter().fold(0u64, |a, &k| {
                a.wrapping_add(u64::from(teb64_ref.get(k as usize)))
            })
        };

        // Each arm gets its own batch, so all four windows are the same
        // *duration* rather than the same iteration count.
        let inner_reps_floor = calibrate(&mut walk_floor);
        let inner_reps_flat = calibrate(&mut walk_flat);
        let inner_reps_teb = calibrate(&mut walk_teb);
        let inner_reps_teb64 = calibrate(&mut walk_teb64);

        for _ in 0..WARMUP {
            black_box(walk_floor());
            black_box(walk_flat());
            black_box(walk_teb());
            black_box(walk_teb64());
        }

        let mut floors = Vec::with_capacity(REPS);
        let mut flats = Vec::with_capacity(REPS);
        let mut tebs = Vec::with_capacity(REPS);
        let mut teb64s = Vec::with_capacity(REPS);
        for _ in 0..REPS {
            floors.push(window(probe, inner_reps_floor, &mut walk_floor));
            flats.push(window(probe, inner_reps_flat, &mut walk_flat));
            tebs.push(window(probe, inner_reps_teb, &mut walk_teb));
            teb64s.push(window(probe, inner_reps_teb64, &mut walk_teb64));
        }
        let floors = per_access(&floors, (ACCESSES * inner_reps_floor) as f64);
        let flats = per_access(&flats, (ACCESSES * inner_reps_flat) as f64);
        let tebs = per_access(&tebs, (ACCESSES * inner_reps_teb) as f64);
        let teb64s = per_access(&teb64s, (ACCESSES * inner_reps_teb64) as f64);

        Row {
            field: field_name,
            resolution: n,
            cells,
            active_cells,
            one_runs: flat.one_runs(),
            flat_bytes: flat.bytes(),
            teb_bytes: teb.bytes(),
            teb_fully_pruned_bytes: fully_pruned.bytes(PAPER_RANK_WORDS),
            // The paper's *basic* TEB: no implicit nodes, no smallest-instance
            // search. It is what the optimisations are worth.
            teb_unoptimised_bytes: fully_pruned.nodes.div_ceil(8)
                + fully_pruned.labels.div_ceil(8)
                + fully_pruned.nodes.div_ceil(PAPER_RANK_WORDS * 64) * 4
                + METADATA_BYTES,
            teb64_bytes: teb64.bytes(),
            teb64_rank_bytes: teb64.rank_bytes(),
            popcounts_per_access: popcounts_total as f64 / cells as f64,
            popcounts_per_access_teb64: popcounts64_total as f64 / cells as f64,
            tree_bytes: teb.tree_bytes(),
            label_bytes: teb.label_bytes(),
            rank_bytes: teb.rank_bytes(),
            tree_nodes: teb.nodes,
            perfect_tree_nodes: 2 * cells - 1,
            prune_depth,
            perfect_levels: teb.perfect_levels,
            tree_height: teb.height,
            truncated_cut_depth: cut,
            answers_equal,
            truncated_tree_mismatches,
            mean_node_visits: visits_total as f64 / cells as f64,
            max_node_visits,
            flat: medians(&flats),
            teb: medians(&tebs),
            teb64: medians(&teb64s),
            floor: medians(&floors),
            ratio: Window {
                cycles: paired(&tebs, &flats, |w| w.cycles),
                instructions: paired(&tebs, &flats, |w| w.instructions),
                ns: paired(&tebs, &flats, |w| w.ns),
            },
            ratio64: Window {
                cycles: paired(&teb64s, &flats, |w| w.cycles),
                instructions: paired(&teb64s, &flats, |w| w.instructions),
                ns: paired(&teb64s, &flats, |w| w.ns),
            },
            inner_reps_flat,
            inner_reps_teb,
            inner_reps_teb64,
            inner_reps_floor,
        }
    }

    fn sweep(probe: &mut Probe) -> Vec<Row> {
        let mut rows = Vec::new();
        for n in RESOLUTIONS {
            isomesh::for_each_reference_field!(f32, |name, field| {
                // The only generic step: turning a field into a bitmap. Every
                // measured arm below sees a `Flat` and nothing else.
                let (_, origin, cell_size) = crate::common::grid(&field, n);
                let flat = active_cell_bitmap(&field, n, origin, cell_size);
                rows.push(measure(probe, name, n, &flat));
            });
        }
        rows
    }

    pub(crate) fn run(run: &mut crate::common::experiment::Run) {
        the_papers_worked_example_decodes();

        let mut probe = Probe::open();
        let rows = sweep(&mut probe);

        println!(
            "\n  popcnt enabled in this build: {}   (`cfg!(target_feature = \"popcnt\")`)",
            cfg!(target_feature = "popcnt")
        );
        println!(
            "\n{:<15} {:>4} {:>8} {:>6} {:>7} {:>8} {:>6} {:>5} {:>7} {:>7} {:>7} {:>7} {:>7} {:>6}",
            "field",
            "n",
            "density",
            "clust",
            "pruned",
            "space",
            "visits",
            "pcnt",
            "cyc_flat",
            "cyc_teb",
            "cyc_t64",
            "ratio",
            "ins_rat",
            "c1c2c3"
        );
        for r in &rows {
            println!(
                "{:<15} {:>4} {:>8.5} {:>6.2} {:>7.4} {:>7.2}× {:>6.2} {:>5.2} {:>7.2} \
                 {:>7.2} {:>7.2} {:>6.1}× {:>6.1}× {:>2}{:>2}{:>2}",
                r.field,
                r.resolution,
                r.density(),
                r.clustering(),
                r.prune_fraction(),
                r.space_ratio(),
                r.mean_node_visits,
                r.popcounts_per_access,
                r.flat.cycles,
                r.teb.cycles,
                r.teb64.cycles,
                r.cycle_ratio(),
                r.instruction_ratio(),
                u8::from(r.c1()),
                u8::from(r.c2()),
                u8::from(r.c3()),
            );
        }

        for r in &rows {
            run.record(&[
                ("field", r.field.to_string()),
                ("resolution", r.resolution.to_string()),
                ("scalar", String::from("f32")),
                ("cells", r.cells.to_string()),
                ("active_cells", r.active_cells.to_string()),
                ("density", format!("{:.8}", r.density())),
                ("one_runs", r.one_runs.to_string()),
                ("clustering_factor", format!("{:.4}", r.clustering())),
                ("flat_bytes", r.flat_bytes.to_string()),
                ("teb_bytes", r.teb_bytes.to_string()),
                ("space_ratio", format!("{:.4}", r.space_ratio())),
                (
                    "teb_fully_pruned_bytes",
                    r.teb_fully_pruned_bytes.to_string(),
                ),
                ("teb_unoptimised_bytes", r.teb_unoptimised_bytes.to_string()),
                ("tree_bytes", r.tree_bytes.to_string()),
                ("label_bytes", r.label_bytes.to_string()),
                ("rank_bytes", r.rank_bytes.to_string()),
                ("metadata_bytes", METADATA_BYTES.to_string()),
                ("tree_nodes", r.tree_nodes.to_string()),
                ("perfect_tree_nodes", r.perfect_tree_nodes.to_string()),
                ("prune_fraction", format!("{:.6}", r.prune_fraction())),
                ("prune_depth", r.prune_depth.to_string()),
                ("perfect_levels", r.perfect_levels.to_string()),
                ("tree_height", r.tree_height.to_string()),
                ("ns_per_access_flat", format!("{:.4}", r.flat.ns)),
                ("ns_per_access_teb", format!("{:.4}", r.teb.ns)),
                ("ns_per_access_floor", format!("{:.4}", r.floor.ns)),
                ("access_ratio", format!("{:.4}", r.cycle_ratio())),
                ("ns_access_ratio", format!("{:.4}", r.ns_ratio())),
                ("mean_node_visits", format!("{:.4}", r.mean_node_visits)),
                ("max_node_visits", r.max_node_visits.to_string()),
                ("cycles_per_access_flat", format!("{:.4}", r.flat.cycles)),
                ("cycles_per_access_teb", format!("{:.4}", r.teb.cycles)),
                ("cycles_per_access_floor", format!("{:.4}", r.floor.cycles)),
                (
                    "instructions_per_access_flat",
                    format!("{:.4}", r.flat.instructions),
                ),
                (
                    "instructions_per_access_teb",
                    format!("{:.4}", r.teb.instructions),
                ),
                (
                    "instructions_per_access_floor",
                    format!("{:.4}", r.floor.instructions),
                ),
                (
                    "instruction_access_ratio",
                    format!("{:.4}", r.instruction_ratio()),
                ),
                // The word-granularity rank arm: the same tree, the same
                // labels, the same run, differing only in the rank block. It
                // is here so C1's falsification cannot be blamed on a tuning
                // constant the paper itself says to determine empirically.
                ("teb64_bytes", r.teb64_bytes.to_string()),
                ("teb64_rank_bytes", r.teb64_rank_bytes.to_string()),
                ("space_ratio_teb64", format!("{:.4}", r.space_ratio64())),
                ("ns_per_access_teb64", format!("{:.4}", r.teb64.ns)),
                ("cycles_per_access_teb64", format!("{:.4}", r.teb64.cycles)),
                (
                    "instructions_per_access_teb64",
                    format!("{:.4}", r.teb64.instructions),
                ),
                ("access_ratio_teb64", format!("{:.4}", r.ratio64.cycles)),
                (
                    "instruction_access_ratio_teb64",
                    format!("{:.4}", r.ratio64.instructions),
                ),
                // `M-281`: this is a fact about the build, not the machine.
                // The 5900X has POPCNT; the default `x86-64` target does not
                // enable it, so every `count_ones` below is a SWAR sequence.
                (
                    "target_feature_popcnt",
                    cfg!(target_feature = "popcnt").to_string(),
                ),
                (
                    "popcounts_per_access_teb",
                    format!("{:.4}", r.popcounts_per_access),
                ),
                (
                    "popcounts_per_access_teb64",
                    format!("{:.4}", r.popcounts_per_access_teb64),
                ),
                ("popcounts_per_access_flat", String::from("0")),
                ("ghz", format!("{:.4}", r.ghz())),
                ("accesses", ACCESSES.to_string()),
                ("inner_reps_flat", r.inner_reps_flat.to_string()),
                ("inner_reps_teb", r.inner_reps_teb.to_string()),
                ("inner_reps_teb64", r.inner_reps_teb64.to_string()),
                ("inner_reps_floor", r.inner_reps_floor.to_string()),
                ("reps", REPS.to_string()),
                ("access_seed", format!("{ACCESS_SEED:#018x}")),
                ("answers_equal", r.answers_equal.to_string()),
                (
                    "truncated_tree_mismatches",
                    r.truncated_tree_mismatches.to_string(),
                ),
                ("truncated_cut_depth", r.truncated_cut_depth.to_string()),
                ("c1_bar", format!("{ACCESS_RATIO_BAR:.1}")),
                ("c2_bar", format!("{SPACE_RATIO_BAR:.1}")),
                ("c1_carrier", String::from("cycles")),
                ("c1_holds", r.c1().to_string()),
                ("c1_holds_instructions", r.c1_instructions().to_string()),
                ("c2_holds", r.c2().to_string()),
                ("c3_holds", r.c3().to_string()),
            ]);
        }
    }
}

fn main() {
    if !std::env::args().any(|arg| arg == "--bench") {
        return;
    }

    let prereg = isomesh::experiment!("P-115");

    #[cfg(target_os = "linux")]
    common::experiment::run(prereg, experiment::run);

    // `experiment_p12`'s precedent. C1 is a cost ratio and `M-281` forbids a
    // nanosecond carrying that verdict, so off Linux there is no instrument and
    // no degraded path — a recorded zero would be a fabricated cost.
    #[cfg(not(target_os = "linux"))]
    {
        eprintln!(
            "{} scores a random-access cost ratio on hardware performance counters, and this \
             platform has no `perf_event_open`. M-281 forbids a nanosecond carrying that verdict.",
            prereg.id
        );
        std::process::exit(1);
    }
}
