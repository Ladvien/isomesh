//! **P-110 — mutable rank/select for a structure written during extraction, scored against the determinism gate.**
//!
//! Ticket: R-110. Pre-registered before this harness existed.
//!
//! ```bash
//! cargo bench --bench experiment_p110
//! ```
//!
//! Writes `docs/experiments/p-110.csv`. **Linux only**, for `experiment_p12`'s
//! reason: C1 is a cost ratio and its verdict reads `perf_event_open`'s
//! instruction counts, not a clock (`✗24`, `M-280`, `M-281`). Off Linux there is
//! nothing to degrade to, so the harness refuses and exits 1 rather than record
//! a fabricated zero.
//!
//! # The gate this row is scored AGAINST rather than around
//!
//! `crates/isomesh/src/validate/determinism.rs:273-303` runs the extractor
//! **three** times — two fresh buffers and a third into a **reused** one — and
//! its own doc comment (`:268-272`) says why: *"two fresh buffers to catch
//! non-determinism in the extractor, and one reused buffer to catch output that
//! depends on the buffer's prior state"*, because *"every algorithm in this
//! crate is meant to be driven by resetting one buffer across thousands of
//! chunks, and nothing else checks that it survives being driven that way"*.
//!
//! That gate is the reason `R-027` was declined and split (`V-45`). A mutable
//! rank/select is a structure **written during** a sweep rather than built once,
//! so it is admissible here only if its final state is a function of the flip
//! **SET** and not of the flip **ORDER**. Order-independence is therefore the
//! clause (C2), not the caveat — and a row that passed C1 and failed C2 would be
//! worthless, which is why C2 is the falsifier that closes the row.
//!
//! # What was missing
//!
//! Nothing in `crates/isomesh/src/` maintains a rank structure across mutations.
//! `dual.rs`'s bitmap is built once per extraction and thrown away, and the only
//! rank-shaped object in the repository is `P-107`'s, which is static. So there
//! was no number for what a *mutable* rank costs against a static one on this
//! machine, and no evidence that a mutable one could survive the determinism
//! gate — only the prior argument, from `V-45`, that a structure carried across
//! calls generally cannot.
//!
//! # The three arms, held identical everywhere except the mechanism
//!
//! All three answer the same question — *how many cell bits are set before
//! `bit`* — over the **same** `Vec<u64>` bitmap, which the caller owns and
//! passes to every query. All three end in the **same call**, `block_tail(words,
//! bit)`, which folds the whole words of `bit`'s own 512-bit block and then its
//! masked word. Neither arm makes any other `count_ones` call, so the popcount
//! work is identical across arms by construction and the columns
//! `count_ones_per_rank_static` and `count_ones_per_rank_mutable` are equal on
//! every row. What differs is only the summary structure above the words.
//!
//! - **`static`** — the comparand, and it is `experiment_p107.rs`'s directory in
//!   shape: a `u16` block rank every `BLOCK_WORDS` = 8 words (512 bits) plus a
//!   `u32` superblock rank every `BLOCKS_PER_SUPER` = 64 blocks (32,768 bits),
//!   **3.2227%** overhead, answering in **≤ 7 words scanned at every
//!   resolution**. `R-107` measured that structure this session (commit
//!   `256b654`); this row reuses it rather than inventing a second meaning for
//!   "rank".
//! - **`mutable`** — Pibiri & Kanda's shape (`10.1016/j.is.2021.101756`): a
//!   b-ary tree of **counts**, not prefix sums, with branching factor `BRANCH` =
//!   64 over the same 512-bit leaf blocks. A `flip(i)` toggles one word bit and
//!   adds ±1 to exactly one slot per level, so an update is `O(log₆₄ n)` writes.
//!   A query scans the slots *before* its own child index at each level, which
//!   is where the paper reaches for SIMD and this crate cannot
//!   (`crates/isomesh/src/` has zero `unsafe`, zero `core::simd` and zero
//!   `#[cfg(target_arch)]`), so the scan here is scalar. **That scan is the
//!   entire cost difference C1 measures.**
//! - **`journal`** — the vacuity control, and a *deliberately* order-sensitive
//!   design rather than a strawman. It is the obvious other way to make a static
//!   directory mutable: keep the static counters, append every flip to a pending
//!   log, rebuild the counters when the log fills (`JOURNAL_CAP` = 1000), and
//!   correct a query by the log entries below its block. Its **answers are
//!   correct** — asserted, `control_answers_correct` — and its **state is not a
//!   function of the flip set**, because after `flips % JOURNAL_CAP ≠ 0` flips
//!   the residual log holds whichever flips arrived last, in the order they
//!   arrived, and the counters hold the state after whichever flips arrived
//!   first. That is `V-45`'s failure mode exactly: output-equal,
//!   state-dependent, invisible to everything except a reused buffer.
//!
//! # The bitmap is mirrored from `dual.rs`, and the mirror is checked
//!
//! `Bitmaps::build` mirrors `dual.rs:359-381` (`build_inside_bits`, one bit per
//! **sample**, packed along X only, `bit_row = size[0].div_ceil(64)`), then
//! `dual.rs:424`'s fused `any & !all` word test masked by `dual.rs:445`'s
//! `cell_mask` into a **cell** bitmap whose row is one word shorter
//! (`cell_words = cells_x.div_ceil(64)`). `M-279`: a mirror licenses nothing
//! until it is shown to reproduce the shipped structure, so every row runs an
//! eight-corner scalar classification of every cell against it
//! (`bitmap_matches_scalar`) and counts the bits past the end of a short row
//! (`pad_bits_set`, asserted 0, and asserted again after the flips). `R-120` and
//! `R-121` both caught real defects that way.
//!
//! The flip set is drawn over **cells**, never over raw bit indices, so a flip
//! cannot land on a pad bit and the padding invariant survives every
//! permutation. What the flips stand for is a sweep editing the field: the
//! post-flip bitmap is no longer the field's own active set, and it is not meant
//! to be — `active_cells_before` and `active_cells_after` are both columns.
//!
//! # The popcount this build does not have, and which way it moves C1
//!
//! There is no `.cargo/config.toml` and no `target-cpu` in the repository, so
//! the default `x86-64` baseline is in force and `u64::count_ones()` **does not
//! lower to `popcnt`** — it lowers to the ~12-instruction SWAR sequence.
//! `cfg!(target_feature = "popcnt")` is false and is the column
//! `target_feature_popcnt`.
//!
//! The contingency is stated exactly rather than modelled. Both arms make the
//! same number of `count_ones` calls per query, because both delegate the same
//! `block_tail` and neither calls it anywhere else — so popcount cost is a
//! **common additive term** in numerator and denominator. Write `rank_static =
//! P + S` and `rank_mutable = P + M`, where `P` is that term. Making `P` cheaper
//! moves `(P + M)/(P + S)` **further from 1** in whichever direction it already
//! is. So:
//!
//! - **`flip` makes zero `count_ones` calls** (`count_ones_per_flip` = 0): it is
//!   one XOR and one ±1 per level. `ns_per_flip` is popcount-independent.
//! - **If `rank_ratio > 1` here, a build with hardware `popcnt` would report a
//!   LARGER ratio**, so a C1 falsified on this build is falsified on that one
//!   too. A C1 that *holds* here is contingent on the absent `popcnt`, and this
//!   row says so rather than leaving it to be discovered.
//!
//! Measuring the other build would mean comparing across binaries, which
//! `M-281` forbids. The call counts are the honest substitute.
//!
//! # SHARE
//!
//! **No time is claimed of extraction**, and no column pretends otherwise. This
//! row has no Amdahl denominator: C1 is a cost bound against a named baseline
//! measured in the same run on the same bitmap, and C2 is the property the
//! crate's determinism gate demands.
//!
//! Every clause's reachable share is a column.
//!
//! - **C1's quantity is `rank_ratio`, and the bar is 2.** The registration's
//!   hypothesis says *"flip(i) plus rank stays within 2x"* while its falsifier
//!   says *"rank costing more than 2x"*; both readings are columns
//!   (`rank_ratio`, `flip_plus_rank_ratio`) and `c1_holds` is their
//!   **conjunction**, which is the stricter of the two and softens neither.
//!   `c1_holds` reads the **instruction** form; `c1_holds_ns` and
//!   `c1_holds_cycles` are beside it. `R-105` watched an identical binary's
//!   cycle ratio move 0.984 → 1.035 across three runs while its instruction
//!   counts held to four figures, so the deterministic form carries the verdict
//!   (`M-280`, `M-281`) and `ghz` is on every row with a nanosecond column.
//! - **C2's quantity is `distinct_final_states`, and the bar is 1.** Its
//!   population is `permutations` = 128 orderings of one flip set of `flips` =
//!   4,096 distinct cell bits, per field. It is an equality over an enumerated
//!   set rather than a ratio, so it has no ceiling.
//! - **C3's quantity is `rank_answers_equal / rank_answers_checked`, and the bar
//!   is 1** — 128 orderings × 1,024 random cell bits against a static directory
//!   rebuilt from *that ordering's own words*, plus one exhaustive sweep of
//!   **every** bit against a running prefix counter (`full_sweep_bits`).
//!
//! **VACUITY CONTROL, asserted rather than recorded.**
//! `control_distinct_final_states` must be **greater than 1** — a deliberately
//! order-sensitive variant the comparator can see. Without it, C2 is a claim
//! that a hash set has one element, proved by a comparator never shown able to
//! hold two. Three further asserts back it up:
//!
//! - `control_answers_correct` — the control is a *working* mechanism whose
//!   **state** is order-sensitive, not broken code. A control that answered
//!   wrongly would be testing the wrong failure.
//! - `flips % JOURNAL_CAP != 0` — with a zero residual the journal's last flush
//!   leaves the counters a pure function of the words, the control would report
//!   `1`, and the vacuity control would itself be vacuous.
//! - `distinct_final_states == 1` ⟺ `final_structures_identical`, where the
//!   second is an **exact bit-for-bit** comparison of every level and every word
//!   against ordering 0's structure. A 64-bit hash collision therefore cannot
//!   report agreement that is not there.
//!
//! **`distinct_final_states` is recorded, not asserted.** It is the clause
//! itself, and `✗51`'s discipline is that a falsified clause with a number is an
//! output; a panic would be a missing row.
//!
//! # Fixture
//!
//! The registration's fixture is **eight reference fields × 65³, 128 seeded
//! permutations of the same flip set per field**. Ordering 0 is the identity —
//! the flip set in its own order — and serves as the reference every other
//! ordering is compared against bit for bit; orderings 1…127 are seeded
//! Fisher–Yates shuffles.
//!
//! 33³ and 129³ are run beside the registered arm because the counter tree is
//! two levels deep at all three (`tree_levels` is a column and reads 2
//! everywhere) while **the width of the level scans is not**: the top level
//! holds 2 slots at 33³, 8 at 65³ and 64 at 129³, and the mean number of `u32`
//! adds a query performs there goes 0.25 → 1.75 → 15.75. C1's whole cost
//! difference is that scan, so a ratio quoted at one width without the
//! neighbouring widths is a number with no error bar. `resolution` and
//! `registered_arm` are columns.

#![allow(
    clippy::cast_possible_truncation,
    clippy::cast_precision_loss,
    clippy::too_many_lines
)]

mod common;

#[cfg(target_os = "linux")]
mod experiment {
    use std::collections::HashSet;
    use std::hint::black_box;
    use std::time::Instant;

    use isomesh::Sdf;
    use isomesh::marching_cubes::table::is_inside;

    use crate::common::counters::{MIN_TIME_RATIO, Probe};

    /// The registered fixture is 65³; the neighbours change the scan width.
    const RESOLUTIONS: [u32; 3] = [33, 65, 129];

    /// The registered arm.
    const REGISTERED_RESOLUTION: u32 = 65;

    /// Orderings of the one flip set, per field. Registered as `permutations`.
    /// Ordering 0 is the identity and is the bit-for-bit reference.
    const PERMUTATIONS: usize = 128;

    /// Distinct cell bits in the flip set. Registered as `flips`.
    ///
    /// Constant across resolutions on purpose: `ns_per_flip` is a per-flip cost,
    /// and a flip count that grew with the grid would make the three
    /// resolutions' figures incomparable.
    const FLIPS: usize = 4096;

    /// Random cell bits per rank window.
    const QUERIES: usize = 2048;

    /// Random cell bits compared per ordering for C3.
    const C3_QUERIES: usize = 1024;

    /// Repetitions per window, medianed per quantity.
    const REPS: usize = 5;

    /// Passes discarded before any window opens.
    const WARMUP: usize = 2;

    /// About this long per counter window, so the ~28 `perf_event` system calls
    /// a window costs land outside it and cannot inflate anything.
    const TARGET_BATCH_NS: f64 = 20_000_000.0;

    /// Ceiling on the batch, so a cheap pass cannot run for a minute.
    const MAX_INNER: usize = 8192;

    /// C1's bar.
    const RATIO_BAR: f64 = 2.0;

    /// Words per leaf block. 512 bits, `experiment_p107.rs`'s geometry.
    const BLOCK_WORDS: usize = 8;

    /// Blocks per superblock in the **static** arm. 32,768 bits, so a rank
    /// within a superblock fits `u16`. `experiment_p107.rs`'s geometry.
    const BLOCKS_PER_SUPER: usize = 64;

    /// Children per internal node of the **mutable** arm's counter tree.
    ///
    /// Equal to [`BLOCKS_PER_SUPER`] so the two arms summarise the same spans
    /// and the comparison is of *mechanism*, not of layout.
    const BRANCH: usize = 64;

    /// `log2(BRANCH)`, the shift from a child index to its parent's.
    const BRANCH_SHIFT: u32 = 6;

    /// Pending flips the **control** buffers before rebuilding its counters.
    ///
    /// [`FLIPS`] `% JOURNAL_CAP != 0` is asserted: with a zero residual the last
    /// flush would leave the counters a pure function of the words, and the
    /// control would report `1` while proving nothing.
    const JOURNAL_CAP: usize = 1000;

    // ─── the bitmaps, mirrored from `dual.rs` ──────────────────────────────

    /// The cell-active bitmap, and the sample-sign bitmap it is folded from.
    ///
    /// Mirrors `dual.rs:359-381`, `:424` and `:445`. The sample bitmap is one
    /// bit per **sample** packed along X only; the cell bitmap's row is one word
    /// shorter, which at 33 samples per axis leaves 32 pad bits per row that
    /// [`Bitmaps::is_cell`] excludes and the caller asserts are clear.
    struct Bitmaps {
        /// `size[0].div_ceil(64)`, the **sample** row.
        bit_row: usize,
        /// One bit per **cell**, `cell_words` per row, `cells_x²` rows.
        active: Vec<u64>,
        cell_words: usize,
        /// Cells per axis, `n − 1`.
        cells_x: usize,
        /// `cell_words * 64 * cells_x²` — the bit space a query may name.
        cell_bits: usize,
    }

    impl Bitmaps {
        /// Pack signs, then fold cells.
        fn build(values: &[f32], n: u32) -> Self {
            let size = [n; 3];
            let sx = n as usize;
            let rows = sx * sx;
            let bit_row = sx.div_ceil(64);
            let mut inside = vec![0u64; bit_row * rows];
            for row in 0..rows {
                let src = sx * row;
                let dst = bit_row * row;
                for w in 0..bit_row {
                    let base = w * 64;
                    let n_bits = (sx - base).min(64);
                    let mut word = 0u64;
                    for k in 0..n_bits {
                        word |= u64::from(is_inside(values[src + base + k])) << k;
                    }
                    inside[dst + w] = word;
                }
            }

            let cells_x = sx - 1;
            let cell_words = cells_x.div_ceil(64);
            let cell_rows = cells_x * cells_x;
            let mut active = vec![0u64; cell_words * cell_rows];
            for z in 0..cells_x {
                for y in 0..cells_x {
                    let dst = cell_words * (z * cells_x + y);
                    for w in 0..cell_words {
                        active[dst + w] =
                            active_word(&inside, bit_row, w, y, z, size) & cell_mask(w, cells_x);
                    }
                }
            }

            Self {
                bit_row,
                active,
                cell_words,
                cells_x,
                cell_bits: cell_words * 64 * cell_rows,
            }
        }

        /// Is cell bit `bit` set?
        #[inline]
        fn get(&self, bit: usize) -> bool {
            self.active[bit >> 6] >> (bit & 63) & 1 == 1
        }

        /// Is `bit` a cell, or a pad bit past the end of a short row?
        #[inline]
        fn is_cell(&self, bit: usize) -> bool {
            bit % (self.cell_words * 64) < self.cells_x
        }

        /// The bit index of the cell at linear index `cell`.
        ///
        /// Every flip and every query is drawn through this, so neither can
        /// name a pad bit.
        #[inline]
        fn bit_of_cell(&self, cell: usize) -> usize {
            let row = cell / self.cells_x;
            let x = cell % self.cells_x;
            row * self.cell_words * 64 + x
        }
    }

    /// Bit `k` is `is_inside(sample[64w + k, y, z])`. `dual.rs:385`.
    #[inline]
    fn inside_word(
        inside: &[u64],
        bit_row: usize,
        w: usize,
        y: usize,
        z: usize,
        size: [u32; 3],
    ) -> u64 {
        inside[bit_row * (y + size[1] as usize * z) + w]
    }

    /// Bit `k` is `is_inside(sample[64w + k + 1, y, z])`. `dual.rs:395`.
    #[inline]
    fn inside_word_shifted(
        inside: &[u64],
        bit_row: usize,
        w: usize,
        y: usize,
        z: usize,
        size: [u32; 3],
    ) -> u64 {
        let lo = inside_word(inside, bit_row, w, y, z, size);
        let hi = if w + 1 < bit_row {
            inside_word(inside, bit_row, w + 1, y, z, size)
        } else {
            0
        };
        (lo >> 1) | (hi << 63)
    }

    /// Sixty-four active-cell answers in four fused word operations. `dual.rs:424`.
    #[inline]
    fn active_word(
        inside: &[u64],
        bit_row: usize,
        w: usize,
        y: usize,
        z: usize,
        size: [u32; 3],
    ) -> u64 {
        let mut any = 0u64;
        let mut all = !0u64;
        for dz in 0..2usize {
            for dy in 0..2usize {
                let a = inside_word(inside, bit_row, w, y + dy, z + dz, size);
                let b = inside_word_shifted(inside, bit_row, w, y + dy, z + dz, size);
                any |= a | b;
                all &= a & b;
            }
        }
        any & !all
    }

    /// Which bits of word `w` are cells that exist. `dual.rs:445`.
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

    // ─── the word-level tail, shared verbatim by all three arms ────────────

    /// Whole words a query at `bit` folds below its block summary.
    ///
    /// The same count in every arm, because every arm calls the same
    /// [`block_tail`]. This is the `count_ones` arithmetic the popcount
    /// contingency is stated in.
    #[inline]
    fn words_folded(bit: usize) -> usize {
        let word = bit >> 6;
        word - (word / BLOCK_WORDS) * BLOCK_WORDS
    }

    /// Fold the words of `bit`'s own block up to, and including, its masked
    /// word. The only `count_ones` any query makes.
    #[inline]
    fn block_tail(words: &[u64], bit: usize) -> u32 {
        let word = bit >> 6;
        let block = word / BLOCK_WORDS;
        let mut r = 0u32;
        for w in &words[(block * BLOCK_WORDS)..word] {
            r += w.count_ones();
        }
        r + (words[word] & ((1u64 << (bit & 63)) - 1)).count_ones()
    }

    // ─── the static arm: `experiment_p107.rs`'s directory ──────────────────

    /// Two levels over the bitmap words: `u16` block ranks within a superblock,
    /// `u32` superblock ranks. Built once, never updated.
    #[derive(Clone, PartialEq, Eq)]
    struct RankDirectory {
        supers: Vec<u32>,
        blocks: Vec<u16>,
    }

    impl RankDirectory {
        /// One pass over the words.
        fn new(words: &[u64]) -> Self {
            let block_count = words.len().div_ceil(BLOCK_WORDS);
            let super_count = block_count.div_ceil(BLOCKS_PER_SUPER);
            let mut supers = Vec::with_capacity(super_count);
            let mut blocks = Vec::with_capacity(block_count);
            let mut total = 0u32;
            for block in 0..block_count {
                if block % BLOCKS_PER_SUPER == 0 {
                    supers.push(total);
                }
                let within = total - supers[block / BLOCKS_PER_SUPER];
                blocks.push(u16::try_from(within).expect(
                    "a superblock spans 32768 bits, so a rank within it fits u16 by construction",
                ));
                let lo = block * BLOCK_WORDS;
                let hi = (lo + BLOCK_WORDS).min(words.len());
                for word in &words[lo..hi] {
                    total += word.count_ones();
                }
            }
            Self { supers, blocks }
        }

        /// Set bits before `bit`. Two directory reads and at most seven folds.
        #[inline]
        fn rank(&self, words: &[u64], bit: usize) -> u32 {
            let block = (bit >> 6) / BLOCK_WORDS;
            self.supers[block / BLOCKS_PER_SUPER]
                + u32::from(self.blocks[block])
                + block_tail(words, bit)
        }

        /// Bytes the directory occupies.
        fn bytes(&self) -> usize {
            self.supers.len() * size_of::<u32>() + self.blocks.len() * size_of::<u16>()
        }
    }

    // ─── the mutable arm: a b-ary tree of counts ───────────────────────────

    /// Pibiri & Kanda's shape: counts rather than prefix sums, so an update
    /// writes one slot per level and a query scans the slots before its own
    /// child index at each level.
    ///
    /// `leaf[b]` is the population of block `b`. `levels[0][i]` is the sum of
    /// leaves `64i … 64i+63`, `levels[k][i]` the sum of `levels[k−1]`'s
    /// children, up to a root of length 1.
    #[derive(Clone, PartialEq, Eq)]
    struct MutableRank {
        leaf: Vec<u16>,
        levels: Vec<Vec<u32>>,
    }

    impl MutableRank {
        fn new(words: &[u64]) -> Self {
            let block_count = words.len().div_ceil(BLOCK_WORDS);
            let mut leaf = Vec::with_capacity(block_count);
            for block in 0..block_count {
                let lo = block * BLOCK_WORDS;
                let hi = (lo + BLOCK_WORDS).min(words.len());
                let population: u32 = words[lo..hi].iter().map(|w| w.count_ones()).sum();
                leaf.push(u16::try_from(population).expect("a block spans 512 bits"));
            }

            let mut levels: Vec<Vec<u32>> = Vec::new();
            let mut below: Vec<u32> = leaf.iter().map(|&c| u32::from(c)).collect();
            loop {
                let len = below.len().div_ceil(BRANCH);
                let mut up = Vec::with_capacity(len);
                for i in 0..len {
                    let lo = i * BRANCH;
                    let hi = (lo + BRANCH).min(below.len());
                    up.push(below[lo..hi].iter().sum::<u32>());
                }
                let done = up.len() == 1;
                below.clone_from(&up);
                levels.push(up);
                if done {
                    break;
                }
            }

            Self { leaf, levels }
        }

        /// Toggle `bit` and repair every counter above it.
        ///
        /// Integer `+1` / `−1` into disjoint slots. That is the whole reason the
        /// final state is a function of the flip **set**: addition on `u16` and
        /// `u32` is commutative and associative, so no slot's value can depend
        /// on the order its increments arrived in, and no slot is ever read
        /// during an update.
        #[inline]
        fn flip(&mut self, words: &mut [u64], bit: usize) {
            let word = bit >> 6;
            words[word] ^= 1u64 << (bit & 63);
            let set = words[word] >> (bit & 63) & 1 == 1;
            let block = word / BLOCK_WORDS;
            if set {
                self.leaf[block] += 1;
            } else {
                self.leaf[block] -= 1;
            }
            let mut idx = block;
            for level in &mut self.levels {
                idx >>= BRANCH_SHIFT;
                if set {
                    level[idx] += 1;
                } else {
                    level[idx] -= 1;
                }
            }
        }

        /// Set bits before `bit`.
        #[inline]
        fn rank(&self, words: &[u64], bit: usize) -> u32 {
            let block = (bit >> 6) / BLOCK_WORDS;
            let mut r = 0u32;
            for (depth, level) in self.levels.iter().enumerate().rev() {
                let idx = block >> (BRANCH_SHIFT * (depth as u32 + 1));
                let start = idx & !(BRANCH - 1);
                for &c in &level[start..idx] {
                    r += c;
                }
            }
            let start = block & !(BRANCH - 1);
            for &c in &self.leaf[start..block] {
                r += u32::from(c);
            }
            r + block_tail(words, bit)
        }

        /// Bytes the tree occupies.
        fn bytes(&self) -> usize {
            self.leaf.len() * size_of::<u16>()
                + self
                    .levels
                    .iter()
                    .map(|level| level.len() * size_of::<u32>())
                    .sum::<usize>()
        }

        /// Every level and every word, so C2 is a claim about the whole
        /// structure rather than about its root.
        fn hash(&self, words: &[u64]) -> u64 {
            let mut h = Fnv::new();
            h.feed(words.len() as u64);
            for &w in words {
                h.feed(w);
            }
            h.feed(self.leaf.len() as u64);
            for &c in &self.leaf {
                h.feed(u64::from(c));
            }
            h.feed(self.levels.len() as u64);
            for level in &self.levels {
                h.feed(level.len() as u64);
                for &c in level {
                    h.feed(u64::from(c));
                }
            }
            h.0
        }
    }

    // ─── the control: a deferred-rebuild journal ───────────────────────────

    /// The static directory plus a pending log, rebuilt when the log fills.
    ///
    /// The obvious other way to make a static structure mutable, and the reason
    /// C2 exists. Its answers are correct; its **state** carries the arrival
    /// order of the last `flips % JOURNAL_CAP` flips and the identity of the
    /// flips that made it into the counters, so two orderings of the same flip
    /// set leave two different structures behind.
    #[derive(Clone, PartialEq, Eq)]
    struct JournalRank {
        directory: RankDirectory,
        /// `(bit, set)` in **arrival order**.
        pending: Vec<(u32, bool)>,
    }

    impl JournalRank {
        fn new(words: &[u64]) -> Self {
            Self {
                directory: RankDirectory::new(words),
                pending: Vec::with_capacity(JOURNAL_CAP),
            }
        }

        #[inline]
        fn flip(&mut self, words: &mut [u64], bit: usize) {
            let word = bit >> 6;
            words[word] ^= 1u64 << (bit & 63);
            let set = words[word] >> (bit & 63) & 1 == 1;
            self.pending.push((bit as u32, set));
            if self.pending.len() == JOURNAL_CAP {
                self.directory = RankDirectory::new(words);
                self.pending.clear();
            }
        }

        /// Set bits before `bit`.
        ///
        /// The counters are as of the last flush and the **words are current**,
        /// so the journal correction applies only to the part of the answer the
        /// counters supply — the rank at the start of `bit`'s own block. The
        /// word fold below that already reflects every pending flip.
        #[inline]
        fn rank(&self, words: &[u64], bit: usize) -> u32 {
            let block = (bit >> 6) / BLOCK_WORDS;
            let block_start_bit = (block * BLOCK_WORDS * 64) as u32;
            let mut base = self.directory.supers[block / BLOCKS_PER_SUPER]
                + u32::from(self.directory.blocks[block]);
            let mut cleared = 0u32;
            for &(idx, set) in &self.pending {
                if idx < block_start_bit {
                    if set {
                        base += 1;
                    } else {
                        cleared += 1;
                    }
                }
            }
            base - cleared + block_tail(words, bit)
        }

        fn hash(&self, words: &[u64]) -> u64 {
            let mut h = Fnv::new();
            h.feed(words.len() as u64);
            for &w in words {
                h.feed(w);
            }
            h.feed(self.directory.supers.len() as u64);
            for &c in &self.directory.supers {
                h.feed(u64::from(c));
            }
            h.feed(self.directory.blocks.len() as u64);
            for &c in &self.directory.blocks {
                h.feed(u64::from(c));
            }
            h.feed(self.pending.len() as u64);
            for &(idx, set) in &self.pending {
                h.feed((u64::from(idx) << 1) | u64::from(set));
            }
            h.0
        }
    }

    // ─── hashing ───────────────────────────────────────────────────────────

    /// FNV-1a on 64-bit lanes with an xor-shift finaliser.
    ///
    /// Only ever used to *count distinct* structures. A collision could report
    /// agreement that is not there, so the authoritative equality is the exact
    /// bit-for-bit comparison in [`measure`], and the two are cross-checked on
    /// every row.
    #[derive(Clone, Copy)]
    struct Fnv(u64);

    impl Fnv {
        const fn new() -> Self {
            Self(0xcbf2_9ce4_8422_2325)
        }

        #[inline]
        fn feed(&mut self, v: u64) {
            self.0 ^= v;
            self.0 = self.0.wrapping_mul(0x0000_0100_0000_01b3);
            self.0 ^= self.0 >> 29;
        }
    }

    // ─── deterministic randomness ──────────────────────────────────────────

    /// splitmix64, so every set and every ordering is the same in every run and
    /// every build.
    struct SplitMix(u64);

    impl SplitMix {
        fn new(seed: u64) -> Self {
            Self(seed)
        }

        fn next_u64(&mut self) -> u64 {
            self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
            let mut z = self.0;
            z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
            z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
            z ^ (z >> 31)
        }

        fn below(&mut self, bound: usize) -> usize {
            (self.next_u64() % bound as u64) as usize
        }
    }

    /// `count` deterministic cell-bit indices, uniform over the cells.
    fn query_set(bits: &Bitmaps, count: usize, seed: u64) -> Vec<usize> {
        let cells = bits.cells_x.pow(3);
        let mut rng = SplitMix::new(seed);
        (0..count)
            .map(|_| bits.bit_of_cell(rng.below(cells)))
            .collect()
    }

    /// [`FLIPS`] **distinct** cell bits.
    ///
    /// Distinctness is what makes "one flip set" a set: a repeated index would
    /// cancel, and the permutation argument would then be about a multiset.
    fn flip_set(bits: &Bitmaps) -> Vec<usize> {
        let cells = bits.cells_x.pow(3);
        assert!(
            cells > FLIPS,
            "the flip set must be a strict subset of the cells"
        );
        let mut rng = SplitMix::new(0x5150_1AD5_0110_0110);
        let mut seen = HashSet::with_capacity(FLIPS);
        let mut out = Vec::with_capacity(FLIPS);
        while out.len() < FLIPS {
            let cell = rng.below(cells);
            if seen.insert(cell) {
                out.push(bits.bit_of_cell(cell));
            }
        }
        out
    }

    /// A seeded ordering of `base`. Fisher–Yates, so every ordering is
    /// reachable and the same seed always gives the same one.
    fn permute(base: &[usize], seed: u64) -> Vec<usize> {
        let mut out = base.to_vec();
        let mut rng = SplitMix::new(seed);
        for i in (1..out.len()).rev() {
            out.swap(i, rng.below(i + 1));
        }
        out
    }

    /// One forward pass and one back, so the window is idempotent.
    ///
    /// A flip is a toggle, so applying the same order twice restores the bitmap
    /// exactly. That is what lets the rank windows measure the same post-flip
    /// population C3 checked, rather than whatever parity the batch size left.
    fn flip_pass(tree: &mut MutableRank, words: &mut [u64], flips: &[usize]) {
        for &bit in flips {
            tree.flip(words, bit);
        }
        for &bit in flips {
            tree.flip(words, bit);
        }
        black_box(&*words);
    }

    // ─── counting ──────────────────────────────────────────────────────────

    /// Cycles, instructions and nanoseconds from one window.
    #[derive(Clone, Copy, Default)]
    struct Counted {
        cycles: f64,
        instructions: f64,
        nanos: f64,
    }

    /// One counter window over `inner` repetitions, divided by `inner`.
    ///
    /// Windows are **siblings, never nested**: Zen 3 has six general-purpose
    /// counters and `Probe` opens six, so a nested window multiplexes and
    /// `Counts::worst_ratio` refuses. `R-121` paid for that discovery.
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
        let scale = 1.0 / inner as f64;
        Counted {
            cycles: counted.cycles.count as f64 * scale,
            instructions: counted.instructions.count as f64 * scale,
            nanos: nanos * scale,
        }
    }

    /// Repetitions to batch into one window, from one timed pass.
    fn choose_inner(mut body: impl FnMut()) -> usize {
        body();
        let started = Instant::now();
        body();
        let pass_ns = started.elapsed().as_nanos() as f64;
        ((TARGET_BATCH_NS / pass_ns.max(1.0)).ceil() as usize).clamp(1, MAX_INNER)
    }

    fn median(values: &mut [f64]) -> f64 {
        values.sort_by(f64::total_cmp);
        values[values.len() / 2]
    }

    fn median_of(reps: &[Counted]) -> Counted {
        let mut cycles: Vec<f64> = reps.iter().map(|c| c.cycles).collect();
        let mut instructions: Vec<f64> = reps.iter().map(|c| c.instructions).collect();
        let mut nanos: Vec<f64> = reps.iter().map(|c| c.nanos).collect();
        Counted {
            cycles: median(&mut cycles),
            instructions: median(&mut instructions),
            nanos: median(&mut nanos),
        }
    }

    // ─── one row ───────────────────────────────────────────────────────────

    struct Row {
        field: &'static str,
        resolution: u32,
        registered_arm: bool,
        cells: usize,
        active_cells_before: u32,
        active_cells_after: u32,
        active_fraction: f64,
        bitmap_words: usize,
        bitmap_bytes: usize,
        cell_words: usize,
        sample_bit_row: usize,
        cell_bits: usize,
        tree_levels: usize,
        top_level_slots: usize,
        static_directory_bytes: usize,
        mutable_directory_bytes: usize,
        static_overhead_fraction: f64,
        mutable_overhead_fraction: f64,
        // C1
        flips: usize,
        ns_per_flip: f64,
        ns_per_rank_mutable: f64,
        ns_per_rank_static: f64,
        ns_per_rank_control: f64,
        rank_ratio: f64,
        flip_plus_rank_ratio: f64,
        cycles_per_flip: f64,
        cycles_per_rank_mutable: f64,
        cycles_per_rank_static: f64,
        rank_ratio_cycles: f64,
        flip_plus_rank_ratio_cycles: f64,
        instructions_per_flip: f64,
        instructions_per_rank_mutable: f64,
        instructions_per_rank_static: f64,
        instructions_per_rank_control: f64,
        rank_ratio_instructions: f64,
        flip_plus_rank_ratio_instructions: f64,
        max_words_scanned: usize,
        count_ones_per_rank_static: f64,
        count_ones_per_rank_mutable: f64,
        count_ones_per_flip: f64,
        target_feature_popcnt: bool,
        // C2
        permutations: usize,
        distinct_final_states: usize,
        final_structures_identical: bool,
        control_distinct_final_states: usize,
        control_structures_identical: bool,
        control_answers_correct: bool,
        // C3
        rank_answers_equal: usize,
        rank_answers_checked: usize,
        full_sweep_bits: usize,
        mutable_equals_prefix_all_bits: usize,
        static_equals_prefix_all_bits: usize,
        // the mirror
        bitmap_matches_scalar: bool,
        pad_bits_set: usize,
        // provenance
        ghz: f64,
        inner_flip: usize,
        inner_rank: usize,
        inner_rank_control: usize,
        c1_holds: bool,
        c1_holds_ns: bool,
        c1_holds_cycles: bool,
        c2_holds: bool,
        c3_holds: bool,
    }

    fn measure<S: Sdf<Scalar = f32>>(
        field: &'static str,
        n: u32,
        sdf: &S,
        origin: [f32; 3],
        cell_size: f32,
    ) -> Row {
        let samples = (n as usize).pow(3);
        let cells_x = n as usize - 1;
        let cells = cells_x.pow(3);

        // ── the values, sampled the way `MarchingCubes::extract` samples ─────
        let mut values = Vec::with_capacity(samples);
        for z in 0..n {
            for y in 0..n {
                for x in 0..n {
                    values.push(sdf.sample([
                        origin[0] + cell_size * x as f32,
                        origin[1] + cell_size * y as f32,
                        origin[2] + cell_size * z as f32,
                    ]));
                }
            }
        }

        let bits = Bitmaps::build(&values, n);
        let active_cells_before: u32 = bits.active.iter().map(|w| w.count_ones()).sum();

        // ── M-279: the mirror has to reproduce the shipped structure ─────────
        //
        // `active_word` is four fused word operations and `cell_mask` discards a
        // whole word at 33 samples per axis. If either is wrong, every clause
        // below is scored over the wrong population.
        let mut bitmap_matches_scalar = true;
        for z in 0..cells_x {
            for y in 0..cells_x {
                for x in 0..cells_x {
                    let mut inside_count = 0u32;
                    for corner in 0..8u8 {
                        let ox = (corner & 1) as usize;
                        let oy = ((corner >> 1) & 1) as usize;
                        let oz = ((corner >> 2) & 1) as usize;
                        let s = (z + oz) * (n as usize) * (n as usize)
                            + (y + oy) * (n as usize)
                            + (x + ox);
                        inside_count += u32::from(is_inside(values[s]));
                    }
                    let scalar_active = inside_count != 0 && inside_count != 8;
                    if bits.get(bits.bit_of_cell((z * cells_x + y) * cells_x + x)) != scalar_active
                    {
                        bitmap_matches_scalar = false;
                    }
                }
            }
        }
        assert!(
            bitmap_matches_scalar,
            "{field} {n}^3: the word-parallel active-cell bitmap disagrees with the eight-corner \
             scalar classification, so every clause below would be scored over the wrong cells"
        );
        assert!(
            active_cells_before > 0,
            "{field} {n}^3: no active cells, so the flip set would be applied to an empty bitmap"
        );

        let pad_bits_set = (0..bits.cell_bits)
            .filter(|&bit| !bits.is_cell(bit) && bits.get(bit))
            .count();
        assert_eq!(
            pad_bits_set, 0,
            "{field} {n}^3: {pad_bits_set} bits past the end of a cell row are set, so \
             `cell_mask` is wrong and every rank after the first short row is too"
        );

        // ── the one flip set, and the control's precondition ─────────────────
        let flips = flip_set(&bits);
        assert_eq!(flips.len(), FLIPS);
        assert_ne!(
            FLIPS % JOURNAL_CAP,
            0,
            "the journal control needs a non-zero residual, or its last flush leaves the \
             counters a pure function of the words, the control reports 1, and the vacuity \
             control is itself vacuous"
        );

        // ── ordering 0: the identity, and the bit-for-bit reference ──────────
        let mut words = bits.active.clone();
        let mut tree = MutableRank::new(&words);
        for &bit in &flips {
            tree.flip(&mut words, bit);
        }
        let mut control_words = bits.active.clone();
        let mut journal = JournalRank::new(&control_words);
        for &bit in &flips {
            journal.flip(&mut control_words, bit);
        }
        assert_eq!(
            words, control_words,
            "{field} {n}^3: the two arms' bitmaps diverged under the same flip set, so a \
             per-ordering rank comparison would be comparing two different populations"
        );
        let directory = RankDirectory::new(&words);
        let active_cells_after: u32 = words.iter().map(|w| w.count_ones()).sum();
        assert_eq!(
            (0..bits.cell_bits)
                .filter(|&bit| !bits.is_cell(bit) && words[bit >> 6] >> (bit & 63) & 1 == 1)
                .count(),
            0,
            "{field} {n}^3: a flip landed on a pad bit"
        );

        // ── C2 and C3 over all `PERMUTATIONS` orderings ──────────────────────
        let c3_queries = query_set(&bits, C3_QUERIES, 0x00C3_0000_0000_00C3);
        let mut hashes = HashSet::with_capacity(PERMUTATIONS);
        let mut control_hashes = HashSet::with_capacity(PERMUTATIONS);
        let mut final_structures_identical = true;
        let mut control_structures_identical = true;
        let mut control_answers_correct = true;
        let mut rank_answers_equal = 0usize;

        hashes.insert(tree.hash(&words));
        control_hashes.insert(journal.hash(&control_words));
        for &bit in &c3_queries {
            if tree.rank(&words, bit) == directory.rank(&words, bit) {
                rank_answers_equal += 1;
            }
            if journal.rank(&control_words, bit) != directory.rank(&control_words, bit) {
                control_answers_correct = false;
            }
        }

        for ordering in 1..PERMUTATIONS {
            let order = permute(&flips, 0xDEC0_0000_0000_0000 ^ ordering as u64);

            let mut other_words = bits.active.clone();
            let mut other_tree = MutableRank::new(&other_words);
            for &bit in &order {
                other_tree.flip(&mut other_words, bit);
            }
            hashes.insert(other_tree.hash(&other_words));
            if other_words != words || other_tree != tree {
                final_structures_identical = false;
            }

            let mut other_control_words = bits.active.clone();
            let mut other_journal = JournalRank::new(&other_control_words);
            for &bit in &order {
                other_journal.flip(&mut other_control_words, bit);
            }
            control_hashes.insert(other_journal.hash(&other_control_words));
            if other_control_words != control_words || other_journal != journal {
                control_structures_identical = false;
            }

            // C3: against a static directory rebuilt from *this* ordering's own
            // words, so an ordering whose bitmap diverged could not hide behind
            // a directory built from someone else's.
            let other_directory = RankDirectory::new(&other_words);
            for &bit in &c3_queries {
                if other_tree.rank(&other_words, bit) == other_directory.rank(&other_words, bit) {
                    rank_answers_equal += 1;
                }
                if other_journal.rank(&other_control_words, bit)
                    != other_directory.rank(&other_control_words, bit)
                {
                    control_answers_correct = false;
                }
            }
        }

        let distinct_final_states = hashes.len();
        let control_distinct_final_states = control_hashes.len();
        assert_eq!(
            distinct_final_states == 1,
            final_structures_identical,
            "{field} {n}^3: the 64-bit structure hash and the exact bit-for-bit comparison \
             disagree, so one of them is lying about C2"
        );
        assert_eq!(
            control_distinct_final_states == 1,
            control_structures_identical,
            "{field} {n}^3: the control's hash and its exact comparison disagree"
        );
        assert!(
            control_distinct_final_states > 1,
            "{field} {n}^3: the deliberately order-sensitive control left the same structure \
             behind under all {PERMUTATIONS} orderings, so the comparator has never been shown \
             able to see disorder and C2's `distinct_final_states == 1` proves nothing"
        );
        assert!(
            control_answers_correct,
            "{field} {n}^3: the control's rank answers are wrong, so it is broken code rather \
             than a working mechanism whose STATE is order-sensitive, and it controls nothing"
        );

        // ── C3's exhaustive half: every bit, against a running prefix ────────
        let mut prefix = 0u32;
        let mut mutable_equals_prefix_all_bits = 0usize;
        let mut static_equals_prefix_all_bits = 0usize;
        for bit in 0..bits.cell_bits {
            if tree.rank(&words, bit) == prefix {
                mutable_equals_prefix_all_bits += 1;
            }
            if directory.rank(&words, bit) == prefix {
                static_equals_prefix_all_bits += 1;
            }
            if words[bit >> 6] >> (bit & 63) & 1 == 1 {
                prefix += 1;
            }
        }
        assert_eq!(
            mutable_equals_prefix_all_bits, bits.cell_bits,
            "{field} {n}^3: the mutable tree disagrees with the prefix sum"
        );
        assert_eq!(
            static_equals_prefix_all_bits, bits.cell_bits,
            "{field} {n}^3: the static directory disagrees with the prefix sum"
        );
        assert_eq!(prefix, active_cells_after);

        // ── the query set, and the popcount arithmetic ───────────────────────
        let queries = query_set(&bits, QUERIES, 0x0A11_0000_0000_0A11);
        let max_words_scanned = queries
            .iter()
            .map(|&bit| words_folded(bit))
            .max()
            .expect("QUERIES is non-zero");
        assert!(
            max_words_scanned < BLOCK_WORDS,
            "{field} {n}^3: a query folded {max_words_scanned} words, which is not bounded by \
             the block width and therefore not the structure either arm claims to be"
        );
        // Both arms call the same `block_tail` and neither calls `count_ones`
        // anywhere else, so this one figure is both arms' per-query popcount
        // count: `words_folded(bit)` whole words plus the masked word.
        let count_ones_per_rank = queries
            .iter()
            .map(|&bit| words_folded(bit) as f64)
            .sum::<f64>()
            / QUERIES as f64
            + 1.0;

        // ── the windows ──────────────────────────────────────────────────────
        let mut flip_words = words.clone();
        let mut flip_tree = tree.clone();

        for _ in 0..WARMUP {
            flip_pass(&mut flip_tree, &mut flip_words, &flips);
            let mut acc = 0u32;
            for &bit in &queries {
                acc = acc.wrapping_add(tree.rank(&words, bit));
                acc = acc.wrapping_add(directory.rank(&words, bit));
            }
            black_box(acc);
        }

        let inner_flip = choose_inner(|| flip_pass(&mut flip_tree, &mut flip_words, &flips));
        let inner_rank = choose_inner(|| {
            let mut acc = 0u32;
            for &bit in &queries {
                acc = acc.wrapping_add(tree.rank(&words, bit));
            }
            black_box(acc);
        });
        // The journal's rank scans its pending log, so it is far dearer and gets
        // its own batch size rather than dragging the other two down. It is a
        // control, not a competitor, and its cost is reported for that reason.
        let inner_rank_control = choose_inner(|| {
            let mut acc = 0u32;
            for &bit in &queries {
                acc = acc.wrapping_add(journal.rank(&control_words, bit));
            }
            black_box(acc);
        });

        let mut probe = Probe::open();
        let mut flip_reps = Vec::with_capacity(REPS);
        let mut rank_mutable = Vec::with_capacity(REPS);
        let mut rank_static = Vec::with_capacity(REPS);
        let mut rank_control = Vec::with_capacity(REPS);

        for _ in 0..REPS {
            flip_reps.push(window(&mut probe, inner_flip, || {
                flip_pass(&mut flip_tree, &mut flip_words, &flips);
            }));
            rank_mutable.push(window(&mut probe, inner_rank, || {
                let mut acc = 0u32;
                for &bit in &queries {
                    acc = acc.wrapping_add(tree.rank(&words, bit));
                }
                black_box(acc);
            }));
            rank_static.push(window(&mut probe, inner_rank, || {
                let mut acc = 0u32;
                for &bit in &queries {
                    acc = acc.wrapping_add(directory.rank(&words, bit));
                }
                black_box(acc);
            }));
            rank_control.push(window(&mut probe, inner_rank_control, || {
                let mut acc = 0u32;
                for &bit in &queries {
                    acc = acc.wrapping_add(journal.rank(&control_words, bit));
                }
                black_box(acc);
            }));
        }
        assert!(
            flip_words == words && flip_tree == tree,
            "{field} {n}^3: the flip windows did not restore the structure, so the rank windows \
             measured a different population from the one C3 checked"
        );

        let flip_counts = median_of(&flip_reps);
        let mutable_counts = median_of(&rank_mutable);
        let static_counts = median_of(&rank_static);
        let control_counts = median_of(&rank_control);

        // Two passes over the flip set per body.
        let per_flip = 1.0 / (2 * FLIPS) as f64;
        let per_query = 1.0 / QUERIES as f64;

        let ns_per_flip = flip_counts.nanos * per_flip;
        let ns_per_rank_mutable = mutable_counts.nanos * per_query;
        let ns_per_rank_static = static_counts.nanos * per_query;
        let instructions_per_flip = flip_counts.instructions * per_flip;
        let instructions_per_rank_mutable = mutable_counts.instructions * per_query;
        let instructions_per_rank_static = static_counts.instructions * per_query;
        let cycles_per_flip = flip_counts.cycles * per_flip;
        let cycles_per_rank_mutable = mutable_counts.cycles * per_query;
        let cycles_per_rank_static = static_counts.cycles * per_query;

        let rank_ratio = ns_per_rank_mutable / ns_per_rank_static;
        let flip_plus_rank_ratio = (ns_per_flip + ns_per_rank_mutable) / ns_per_rank_static;
        let rank_ratio_cycles = cycles_per_rank_mutable / cycles_per_rank_static;
        let flip_plus_rank_ratio_cycles =
            (cycles_per_flip + cycles_per_rank_mutable) / cycles_per_rank_static;
        let rank_ratio_instructions = instructions_per_rank_mutable / instructions_per_rank_static;
        let flip_plus_rank_ratio_instructions =
            (instructions_per_flip + instructions_per_rank_mutable) / instructions_per_rank_static;

        let bitmap_bytes = words.len() * size_of::<u64>();

        // The verdict reads the instruction form: `R-105` watched an identical
        // binary's cycle ratio move 0.984 -> 1.035 across three runs while its
        // instruction counts held to four figures (`M-280`, `M-281`). Both
        // readings of the registered clause must clear the bar, which is the
        // conjunction of the hypothesis's wording and its falsifier's.
        let c1_holds =
            rank_ratio_instructions <= RATIO_BAR && flip_plus_rank_ratio_instructions <= RATIO_BAR;
        let c2_holds = distinct_final_states == 1 && final_structures_identical;
        let rank_answers_checked = PERMUTATIONS * C3_QUERIES;
        let c3_holds = rank_answers_equal == rank_answers_checked
            && mutable_equals_prefix_all_bits == bits.cell_bits
            && static_equals_prefix_all_bits == bits.cell_bits;

        Row {
            field,
            resolution: n,
            registered_arm: n == REGISTERED_RESOLUTION,
            cells,
            active_cells_before,
            active_cells_after,
            active_fraction: f64::from(active_cells_after) / cells as f64,
            bitmap_words: words.len(),
            bitmap_bytes,
            cell_words: bits.cell_words,
            sample_bit_row: bits.bit_row,
            cell_bits: bits.cell_bits,
            tree_levels: tree.levels.len(),
            top_level_slots: tree.levels.first().map_or(0, Vec::len),
            static_directory_bytes: directory.bytes(),
            mutable_directory_bytes: tree.bytes(),
            static_overhead_fraction: directory.bytes() as f64 / bitmap_bytes as f64,
            mutable_overhead_fraction: tree.bytes() as f64 / bitmap_bytes as f64,
            flips: FLIPS,
            ns_per_flip,
            ns_per_rank_mutable,
            ns_per_rank_static,
            ns_per_rank_control: control_counts.nanos * per_query,
            rank_ratio,
            flip_plus_rank_ratio,
            cycles_per_flip,
            cycles_per_rank_mutable,
            cycles_per_rank_static,
            rank_ratio_cycles,
            flip_plus_rank_ratio_cycles,
            instructions_per_flip,
            instructions_per_rank_mutable,
            instructions_per_rank_static,
            instructions_per_rank_control: control_counts.instructions * per_query,
            rank_ratio_instructions,
            flip_plus_rank_ratio_instructions,
            max_words_scanned,
            count_ones_per_rank_static: count_ones_per_rank,
            count_ones_per_rank_mutable: count_ones_per_rank,
            count_ones_per_flip: 0.0,
            target_feature_popcnt: cfg!(target_feature = "popcnt"),
            permutations: PERMUTATIONS,
            distinct_final_states,
            final_structures_identical,
            control_distinct_final_states,
            control_structures_identical,
            control_answers_correct,
            rank_answers_equal,
            rank_answers_checked,
            full_sweep_bits: bits.cell_bits,
            mutable_equals_prefix_all_bits,
            static_equals_prefix_all_bits,
            bitmap_matches_scalar,
            pad_bits_set,
            ghz: mutable_counts.cycles / mutable_counts.nanos,
            inner_flip,
            inner_rank,
            inner_rank_control,
            c1_holds,
            c1_holds_ns: rank_ratio <= RATIO_BAR && flip_plus_rank_ratio <= RATIO_BAR,
            c1_holds_cycles: rank_ratio_cycles <= RATIO_BAR
                && flip_plus_rank_ratio_cycles <= RATIO_BAR,
            c2_holds,
            c3_holds,
        }
    }

    /// Eight reference fields × three resolutions, `f32`.
    ///
    /// No `scalar` column is registered and none is added: every arm is integer
    /// work over a bitmap, and the field's precision changes only which cells
    /// start active — which `active_fraction` already reports.
    fn sweep() -> Vec<Row> {
        let mut rows = Vec::new();
        for n in RESOLUTIONS {
            isomesh::for_each_reference_field!(f32, |name, field| {
                let (_, origin, cell_size) = crate::common::grid(&field, n);
                rows.push(measure(name, n, &field, origin, cell_size));
            });
        }
        rows
    }

    pub(crate) fn run(run: &mut crate::common::experiment::Run) {
        let rows = sweep();

        let registered: Vec<&Row> = rows.iter().filter(|r| r.registered_arm).collect();
        let worst_instructions = registered
            .iter()
            .map(|r| r.rank_ratio_instructions)
            .fold(0.0f64, f64::max);
        let worst_ns = registered
            .iter()
            .map(|r| r.rank_ratio)
            .fold(0.0f64, f64::max);
        let best_ns = registered
            .iter()
            .map(|r| r.rank_ratio)
            .fold(f64::INFINITY, f64::min);

        println!(
            "P-110: C1 on the registered 65^3 arm - rank_ratio_instructions at most \
             {worst_instructions:.4} against the bar of {RATIO_BAR:.1}; the nanosecond form spans \
             {best_ns:.4}-{worst_ns:.4} and carries no verdict"
        );
        println!(
            "P-110: C1 clears the bar on {} of {} rows on instructions and {} on nanoseconds",
            rows.iter().filter(|r| r.c1_holds).count(),
            rows.len(),
            rows.iter().filter(|r| r.c1_holds_ns).count()
        );
        println!(
            "P-110: C2 over {PERMUTATIONS} orderings of {FLIPS} flips - distinct_final_states is \
             at most {} on any row, and the order-sensitive control reports {} to {}",
            rows.iter()
                .map(|r| r.distinct_final_states)
                .max()
                .unwrap_or(0),
            rows.iter()
                .map(|r| r.control_distinct_final_states)
                .min()
                .unwrap_or(0),
            rows.iter()
                .map(|r| r.control_distinct_final_states)
                .max()
                .unwrap_or(0)
        );
        println!(
            "P-110: target_feature_popcnt is {}, both arms make {:.3} count_ones calls per rank \
             and flip makes none, so a hardware popcount would move the ratio further from 1",
            cfg!(target_feature = "popcnt"),
            rows[0].count_ones_per_rank_static
        );

        for row in &rows {
            run.record(&[
                ("field", row.field.to_string()),
                ("resolution", row.resolution.to_string()),
                ("registered_arm", row.registered_arm.to_string()),
                ("cells", row.cells.to_string()),
                ("active_cells_before", row.active_cells_before.to_string()),
                ("active_cells_after", row.active_cells_after.to_string()),
                ("active_fraction", format!("{:.6}", row.active_fraction)),
                ("bitmap_words", row.bitmap_words.to_string()),
                ("bitmap_bytes", row.bitmap_bytes.to_string()),
                ("cell_words_per_row", row.cell_words.to_string()),
                ("sample_bit_row", row.sample_bit_row.to_string()),
                ("cell_bits", row.cell_bits.to_string()),
                ("tree_levels", row.tree_levels.to_string()),
                ("top_level_slots", row.top_level_slots.to_string()),
                (
                    "static_directory_bytes",
                    row.static_directory_bytes.to_string(),
                ),
                (
                    "mutable_directory_bytes",
                    row.mutable_directory_bytes.to_string(),
                ),
                (
                    "static_overhead_fraction",
                    format!("{:.6}", row.static_overhead_fraction),
                ),
                (
                    "mutable_overhead_fraction",
                    format!("{:.6}", row.mutable_overhead_fraction),
                ),
                ("flips", row.flips.to_string()),
                ("permutations", row.permutations.to_string()),
                ("ns_per_flip", format!("{:.4}", row.ns_per_flip)),
                (
                    "ns_per_rank_mutable",
                    format!("{:.4}", row.ns_per_rank_mutable),
                ),
                (
                    "ns_per_rank_static",
                    format!("{:.4}", row.ns_per_rank_static),
                ),
                (
                    "ns_per_rank_control",
                    format!("{:.4}", row.ns_per_rank_control),
                ),
                ("rank_ratio", format!("{:.4}", row.rank_ratio)),
                (
                    "flip_plus_rank_ratio",
                    format!("{:.4}", row.flip_plus_rank_ratio),
                ),
                ("cycles_per_flip", format!("{:.3}", row.cycles_per_flip)),
                (
                    "cycles_per_rank_mutable",
                    format!("{:.3}", row.cycles_per_rank_mutable),
                ),
                (
                    "cycles_per_rank_static",
                    format!("{:.3}", row.cycles_per_rank_static),
                ),
                ("rank_ratio_cycles", format!("{:.4}", row.rank_ratio_cycles)),
                (
                    "flip_plus_rank_ratio_cycles",
                    format!("{:.4}", row.flip_plus_rank_ratio_cycles),
                ),
                (
                    "instructions_per_flip",
                    format!("{:.3}", row.instructions_per_flip),
                ),
                (
                    "instructions_per_rank_mutable",
                    format!("{:.3}", row.instructions_per_rank_mutable),
                ),
                (
                    "instructions_per_rank_static",
                    format!("{:.3}", row.instructions_per_rank_static),
                ),
                (
                    "instructions_per_rank_control",
                    format!("{:.3}", row.instructions_per_rank_control),
                ),
                (
                    "rank_ratio_instructions",
                    format!("{:.4}", row.rank_ratio_instructions),
                ),
                (
                    "flip_plus_rank_ratio_instructions",
                    format!("{:.4}", row.flip_plus_rank_ratio_instructions),
                ),
                ("max_words_scanned", row.max_words_scanned.to_string()),
                (
                    "count_ones_per_rank_static",
                    format!("{:.3}", row.count_ones_per_rank_static),
                ),
                (
                    "count_ones_per_rank_mutable",
                    format!("{:.3}", row.count_ones_per_rank_mutable),
                ),
                (
                    "count_ones_per_flip",
                    format!("{:.3}", row.count_ones_per_flip),
                ),
                (
                    "target_feature_popcnt",
                    row.target_feature_popcnt.to_string(),
                ),
                (
                    "distinct_final_states",
                    row.distinct_final_states.to_string(),
                ),
                (
                    "final_structures_identical",
                    row.final_structures_identical.to_string(),
                ),
                (
                    "control_distinct_final_states",
                    row.control_distinct_final_states.to_string(),
                ),
                (
                    "control_structures_identical",
                    row.control_structures_identical.to_string(),
                ),
                (
                    "control_answers_correct",
                    row.control_answers_correct.to_string(),
                ),
                ("rank_answers_equal", row.rank_answers_equal.to_string()),
                ("rank_answers_checked", row.rank_answers_checked.to_string()),
                ("full_sweep_bits", row.full_sweep_bits.to_string()),
                (
                    "mutable_equals_prefix_all_bits",
                    row.mutable_equals_prefix_all_bits.to_string(),
                ),
                (
                    "static_equals_prefix_all_bits",
                    row.static_equals_prefix_all_bits.to_string(),
                ),
                (
                    "bitmap_matches_scalar",
                    row.bitmap_matches_scalar.to_string(),
                ),
                ("pad_bits_set", row.pad_bits_set.to_string()),
                ("ghz", format!("{:.4}", row.ghz)),
                ("inner_flip", row.inner_flip.to_string()),
                ("inner_rank", row.inner_rank.to_string()),
                ("inner_rank_control", row.inner_rank_control.to_string()),
                ("c1_holds", row.c1_holds.to_string()),
                ("c1_holds_ns", row.c1_holds_ns.to_string()),
                ("c1_holds_cycles", row.c1_holds_cycles.to_string()),
                ("c2_holds", row.c2_holds.to_string()),
                ("c3_holds", row.c3_holds.to_string()),
            ]);
        }
    }
}

fn main() {
    if !std::env::args().any(|arg| arg == "--bench") {
        return;
    }

    let prereg = isomesh::experiment!("P-110");

    #[cfg(target_os = "linux")]
    common::experiment::run(prereg, experiment::run);

    // `experiment_p12`'s precedent. C1 is a cost ratio whose verdict reads
    // instruction counts, and a nanosecond on a governed CPU cannot carry it
    // (`M-281`). A recorded zero would be a fabricated ratio.
    #[cfg(not(target_os = "linux"))]
    {
        eprintln!(
            "{} scores a cost ratio on hardware performance counters, and this platform has no \
             `perf_event_open`. There is no clock substitute.",
            prereg.id
        );
        std::process::exit(1);
    }
}
