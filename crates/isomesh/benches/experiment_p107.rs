//! **P-107 — a rank directory over the active-cell bitmap gives the output slot index in O(1).**
//!
//! Ticket: R-107. Pre-registered before this harness existed.
//!
//! ```bash
//! cargo bench --bench experiment_p107
//! ```
//!
//! Writes `docs/experiments/p-107.csv`. **Linux only**, for `experiment_p12`'s
//! reason: C1 is a share of extraction and C2 is a query cost, and on a governed
//! CPU a nanosecond is not a unit (`✗24`, `M-280`, `M-281`). The verdicts read
//! `perf_event_open`; off Linux there is nothing to degrade to, so the harness
//! refuses and exits 1 rather than record a fabricated zero.
//!
//! # What was missing
//!
//! Flying Edges spends a whole prefix-sum pass turning per-row counts into
//! output offsets. A rank directory turns that pass into a lookup. Before this
//! harness there was no number in the repository for **how much of extraction
//! the offset work is**, so the mechanism could not be scored — which is exactly
//! the shape of `✗51`, where `P-69`'s C1 asked for 2× on a stage that is 11.6%
//! of the marginal cost and had a ceiling of `1/(1 − 0.116/2) = 1.061×`.
//!
//! Two bounds existed and neither is the quantity. `M-135`'s
//! `docs/measurements/stage_breakdown.csv` decomposes the **pipeline**
//! (`contour / normals / weld / collider`, mean shares 29.0 / 0.4 / 25.5 / 45.0),
//! so the whole contour stage is 29% and the offset work is some fraction of
//! that. `✗54` measured the **GPU** scan at 4.37% of `gpu_total_ms` with upload
//! at 87.50% — the same quantity on the other path, and already *below* this
//! clause's bar, but not transferable: on the GPU the scan competes with a PCIe
//! upload that dominates the frame, while on the CPU there is no upload at all.
//!
//! `P-121` has since measured the closest CPU quantity there is, and it is this
//! row's **prior, not its measurement**: `(cycles_emit_prepare +
//! cycles_emit_walk) / cycles_total` is **0.0800** at `sphere 65³ f32
//! marching_cubes` and **0.0073–0.0112** on every `dual_contouring` row. That
//! numerator is a strict *superset* of the offset work — it also carries the
//! `CASES` lookup, the triangle push and the quad emission — so the offset share
//! measured here must come out below it, and 5% was already known to be tight.
//! It is measured here anyway, because a share inherited from a neighbouring
//! definition is how `✗51` happened.
//!
//! # What "offset/compaction" means here, stated before it is measured
//!
//! The registration's C1 is a share, and a share is meaningless until its
//! numerator is pinned. The rank directory replaces exactly one thing: *given a
//! cell, what is its output slot index?* So the numerator is the work that turns
//! the per-cell active predicate into a dense slot index, and nothing else:
//!
//! - **`offset_prepare`** — clearing and refilling the dense slot table
//!   (`slots.clear(); slots.resize(cells, u32::MAX)`). This exists **only**
//!   because slots are addressed densely instead of compacted; a directory
//!   removes it outright, because the slot *is* `rank(cell)`. It is the same
//!   object `marching_cubes/mod.rs:250-251` resizes on every extract — the crate
//!   keys it by `(sample, axis)` and so pays `3·n³` slots (3.3 MB at 65³), while
//!   this harness keys it by cell and pays `(n−1)³` (1.0 MB at 65³). The
//!   harness's numerator is therefore about **3× smaller than the crate's own**,
//!   which biases C1 towards falsification. Both figures are columns
//!   (`slot_table_bytes`, `crate_edge_table_bytes`).
//! - **`offset_scan`** — the prefix-sum pass itself: `dual.rs:489-497`'s set-bit
//!   walk (`x = trailing_zeros`, `active &= active − 1`, ascending `x`) with a
//!   running counter, writing `slots[cell] = next++`.
//!
//! **Excluded, deliberately:** building the bitmap. `P-121` assigned the
//! bit-packing (`dual.rs:359-381`), the fused `any & !all` word test (`:424`)
//! and the sign scan to `classify`, and this harness keeps that boundary so the
//! two rows' shares are comparable. Counting the bitmap build as offset work
//! would inflate C1 with sign tests the directory does not remove.
//!
//! `offset_share_scan_only` is beside `offset_share` so a reader who rejects the
//! memset can score the clause on the pass alone without re-deriving anything.
//!
//! # The denominator, and both of them
//!
//! `cycles_extract_mc` and `cycles_extract_dc` are counter windows over the
//! **shipped** `MarchingCubes::extract` and `DualContouring::extract` on the
//! same grid, in the same run. Nothing is mirrored: the denominator is the real
//! extractor, so there is no mirror-agreement debt to pay.
//!
//! `offset_share` is scored against **`marching_cubes`**, which is the path with
//! a dense slot table and the larger prior. `offset_share_dc` and `c1_holds_dc`
//! are beside it. Choosing the *more favourable* denominator for the registered
//! column is deliberate: a C1 that fails there fails everywhere, and reporting
//! the other one in the same row means the choice cannot hide anything.
//!
//! Numerator and denominator are both batched windows over resident buffers, so
//! the comparison is like-for-like. What the numerator does **not** have is the
//! competing loads an in-situ offset pass would face, so the reading is the
//! offset pass at its cheapest — stated here rather than discovered later.
//!
//! # The directory, and why its O(1) is a count and not a clock
//!
//! Two levels over the **cell** bitmap words, sized to land inside pasta's
//! 3.51%: a `u16` block rank every [`BLOCK_WORDS`] words (512 bits) and a `u32`
//! superblock rank every [`BLOCKS_PER_SUPER`] blocks (32,768 bits). The
//! arithmetic is uniform in the fixture — `2/64 + 4/4096 = 3.2227%` — and at 65³
//! that is 1,056 bytes over a 32,768-byte bitmap, which is the registration's
//! own *"about 1.1 KiB on a 64³ chunk's 32 KiB bitmap"* recovered from the
//! layout rather than quoted at it.
//!
//! A query reads two directory entries and then folds **at most
//! `BLOCK_WORDS − 1` = 7 words**. That bound is what "O(1)" means, and it is
//! reported as `max_words_scanned_rank` — a deterministic integer, asserted
//! `≤ 7`, identical at 33³, 65³ and 129³ — beside `max_words_scanned_scan`,
//! which grows with the bitmap. A ratio of two nanosecond figures on a governed
//! CPU cannot carry that clause; a count can, and the count is the same in every
//! build.
//!
//! # What the comparands are, including the one that wins
//!
//! Three answers to the same question are measured:
//!
//! - **`rank`** — the directory. Random access, O(1).
//! - **`scan`** — a full prefix scan from word 0 for each query. This is the
//!   registered `ns_per_query_scan`: the honest cost of answering a *random*
//!   slot query without a directory.
//! - **`sequential`** — `ns_per_query_sequential`, an extra column and the one
//!   that matters: the shipped path never asks a random question. It walks the
//!   cells in order with a running counter, which is O(1) *amortised* per active
//!   cell and cannot be beaten by a directory. **The directory buys random
//!   access, not throughput**, and a row that reported only `rank` against
//!   `scan` would be quoting a comparand the crate does not use.
//!
//! # The popcount this build does not have
//!
//! There is no `.cargo/config.toml` and no `target-cpu` anywhere in the
//! repository, so the default `x86-64` baseline is in force and
//! `u64::count_ones()` **does not lower to `popcnt`** — it lowers to the
//! twelve-instruction SWAR sequence, and LLVM vectorises the multi-word loops
//! with an SSE2 `psadbw` prologue. Checked on this bench's own binary rather
//! than assumed: `objdump -d` greps **0** `popcnt` and **91** occurrences of the
//! SWAR constant `0x3333333333333333`. `cfg!(target_feature = "popcnt")` is
//! false and is a column, `target_feature_popcnt`.
//!
//! That matters because pasta's rank block, Vigna's broadword rank/select and
//! VDB's popcount child offset are all priced against a one-cycle popcount, and
//! this build executes something an order of magnitude dearer. So the call
//! counts are columns too, and the arithmetic is stated rather than left to be
//! re-derived:
//!
//! - **The `rank` arm makes `count_ones_per_query_rank` ≈ 4.5 calls** — the
//!   mean of `words_folded` over a uniform query set is `(BLOCK_WORDS − 1)/2 =
//!   3.5`, plus the final masked word.
//! - **The `scan` arm makes `count_ones_per_query_scan` calls**, which is half
//!   the bitmap in words and grows with resolution — 2,048 at 33³ and 16,384 at
//!   129³.
//! - **The offset stage makes none.** `Offset::scan` is `trailing_zeros` and
//!   `word &= word − 1`, and the bitmap build is shifts and boolean folds. **C1
//!   is therefore popcount-independent**, which is worth saying because C1 is
//!   the clause that gates the row.
//!
//! So a hardware popcount would move `rank` and `scan` in the *same* direction
//! and leave C1 untouched. No verdict here is contingent on it. Measuring it
//! would mean comparing across binaries, which `M-281` forbids, so it is not
//! measured — the call counts are the honest substitute.
//!
//! # The API `R-112` consumes
//!
//! `P-112` replaces the middle phase of count/scan/scatter with this
//! directory, so the mechanism has to be usable and not only measured. Every
//! Phase 25 row is bench-local, so `R-112` copies these three items rather than
//! importing them, and the shape is fixed here so both rows mean the same thing
//! by "rank":
//!
//! ```ignore
//! struct Bitmaps { bit_row, active: Vec<u64>, cell_words, cells_x, cell_bits }
//! Bitmaps::build(values: &[f32], n: u32) -> Bitmaps   // dual.rs:359-381 + :424 + :445
//! Bitmaps::get(&self, bit: usize) -> bool
//!
//! RankDirectory::new(words: &[u64]) -> RankDirectory
//! RankDirectory::rank(&self, words: &[u64], bit: usize) -> u32   // the output slot
//! RankDirectory::total(&self, words: &[u64]) -> u32
//! RankDirectory::bytes(&self) -> usize
//! ```
//!
//! Three properties of that shape are deliberate. **The words stay with the
//! caller** and are passed to every query, so the directory is a side table
//! `R-112`'s scatter phase can hold beside its own bitmap without moving it.
//! **`rank(bit)` is exclusive** — it counts set bits *before* `bit` — so the
//! slot of an active cell is `rank` of its own bit and no ±1 correction is
//! needed anywhere. **The bit index is the packed cell index**, `(z·cells_x +
//! y)·cell_words·64 + x`, not the linear cell index: the pad bits at the end of
//! a short row are asserted clear (`pad_bits_set`), which is what makes the two
//! indexings agree on every rank.
//!
//! # SHARE
//!
//! Every clause's reachable share is a column.
//!
//! - **C1's share is `offset_share` itself** — the clause *is* the share, and it
//!   gates the row. The bar is 0.05, per row, against `cycles_extract_mc`.
//!   `offset_share_dc`, `offset_share_scan_only` and
//!   `offset_share_instructions` are beside it, and the run prints how many rows
//!   clear the bar under each denominator before it writes anything.
//!   A row under the bar is a **falsified C1 on that row**, recorded as such —
//!   `✗51`'s precedent, and the cheapest useful output this row has.
//! - **C2's share is a space fraction and a word count, not a time**:
//!   `overhead_fraction ≤ 0.0351` **and** `max_words_scanned_rank ≤ 7` with no
//!   dependence on resolution. Neither is a ratio of totals, so `✗51`'s
//!   `1/(1 − share/factor)` bar does not apply to it.
//! - **C3's share is `slots_equal / active_cells`, and the bar is 1.** It is an
//!   equality over an enumerated population — every cell bit, `cells_checked` of
//!   them — so its denominator is exact by construction. It moves no time and
//!   has no ceiling.
//!
//! **VACUITY CONTROL, asserted rather than recorded.** `slots_equal` must equal
//! `active_cells`; `active_cells` must be non-zero; and
//! `short_directory_mismatches` — a directory built **one level short**, the
//! superblock array alone with the block level missing and the query still
//! jumping to the block boundary — must be greater than zero. Without that last
//! one, C3 is an equality between two names for the same computation: the
//! comparand *is* the running prefix counter, so a comparator that cannot fail
//! proves nothing. `bitmap_matches_scalar` is the fourth: the word-parallel
//! bitmap is checked against an eight-corner scalar classification of every
//! cell, so the population C3 is asserted over is the right one.
//!
//! # `ghz` is provenance, and this harness has already needed it
//!
//! `M-280` and `M-281`. `ns_per_query_rank` and `ns_per_query_scan` are
//! registered columns and are reported; `ghz` is on the row so a later reader
//! can see what clock they were taken at. C1 reads cycles, C2 reads bytes and a
//! word count, C3 reads an equality. No clause consults a nanosecond.
//!
//! It is not decoration. In an **earlier, uncommitted run of this harness on
//! this machine**, `sphere 129³` read `ghz = 3.11` against 4.18 on all
//! twenty-three other rows, and its `ns_per_query_rank` read **11.89 ns**
//! against 6.7–6.9 on its neighbours — while `instructions_per_query_rank` read
//! **127.262**, digit for digit the same as `torus 129³`. Same binary, same
//! fixture: the instruction form reproduced and the nanosecond did not. That is
//! why C2's O(1) is scored on `max_words_scanned_rank`, and why the clock
//! figures on this row are reported and never consulted.

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

    use isomesh::dual_contouring::DualContouring;
    use isomesh::marching_cubes::MarchingCubes;
    use isomesh::marching_cubes::table::is_inside;
    use isomesh::{MeshBuffer, RuntimeShape3, Sdf};

    use crate::common::counters::{MIN_TIME_RATIO, Probe};

    /// The registered fixture: eight reference fields at each of these.
    const RESOLUTIONS: [u32; 3] = [33, 65, 129];

    /// Repetitions per window, medianed per quantity.
    ///
    /// Five rather than `experiment_p121`'s nine because every repetition here
    /// runs **two whole shipped extractions** at up to 129³, and the quantities
    /// that carry verdicts are a share of cycles, a byte count and an integer
    /// equality — only the first is noisy at all.
    const REPS: usize = 5;

    /// Passes discarded before any window opens.
    const WARMUP: usize = 2;

    /// About this long per counter window, so the ~28 `perf_event` system calls
    /// a window costs land outside it and cannot inflate anything.
    const TARGET_BATCH_NS: f64 = 20_000_000.0;

    /// Ceiling on the batch, so a cheap pass cannot run for a minute.
    const MAX_INNER: usize = 8192;

    /// Random slot queries per window.
    const QUERIES: usize = 2048;

    /// C2's space bar: pasta's figure, cross-confirmed in two acquired papers.
    const OVERHEAD_BAR: f64 = 0.0351;

    /// C1's share bar.
    const SHARE_BAR: f64 = 0.05;

    /// Words per block. 512 bits, the flat-rank convention.
    pub(crate) const BLOCK_WORDS: usize = 8;

    /// Blocks per superblock. 32,768 bits, so a within-superblock rank fits
    /// `u16` with room to spare.
    pub(crate) const BLOCKS_PER_SUPER: usize = 64;

    // ─── the bitmaps, mirrored from `dual.rs` ──────────────────────────────

    /// The cell-active bitmap, and the sample-sign bitmap it is folded from.
    ///
    /// The sample bitmap mirrors `dual.rs:359-381` `build_inside_bits`: **one
    /// bit per sample**, 64 to a `u64`, packed **along X only**, `bit_row =
    /// size[0].div_ceil(64)`. It is consumed here and not kept, because the
    /// directory sits over cells. `active` is `dual.rs:424`'s fused `any & !all`
    /// word test masked by `dual.rs:445`'s `cell_mask`, and its row is one word
    /// shorter — `cell_words = cells_x.div_ceil(64)` — which at 65 and 129
    /// samples per axis is a real asymmetry and not a rounding detail.
    struct Bitmaps {
        /// `size[0].div_ceil(64)`, the **sample** row, reported so the
        /// cell/sample word asymmetry is visible in the CSV.
        bit_row: usize,
        /// One bit per **cell**, `cell_words` per row, `cells²` rows.
        active: Vec<u64>,
        cell_words: usize,
        /// Cells per axis, `n − 1`.
        cells_x: usize,
        /// `cell_words * 64 * cells²` — the bit space a query may name.
        cell_bits: usize,
    }

    impl Bitmaps {
        /// Pack signs, then fold cells. `values` is row-major with stride `n`,
        /// which is the stride `MarchingCubes::extract` samples on.
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

        /// Is bit `bit` a cell, or one of the pad bits past the end of a row?
        ///
        /// Padding bits — the ones `cell_mask` cleared, at `x ≥ cells_x` in a
        /// row whose width is not a multiple of 64 — must be **zero**, or the
        /// rank of every later cell is wrong. Asserted per row rather than
        /// assumed, because `cell_mask` is the one line of the mirror with no
        /// counterpart in the scalar cross-check.
        #[inline]
        fn is_cell(&self, bit: usize) -> bool {
            bit % (self.cell_words * 64) < self.cells_x
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
    ///
    /// The high bit has to come from the next word or the cell straddling a word
    /// boundary reads its `+x` corner as outside — a hole every 64 cells.
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

    // ─── the rank directory ────────────────────────────────────────────────

    /// Two levels over the bitmap words: `u16` block ranks, `u32` superblock
    /// ranks. The bitmap stays with the caller and is passed to the query, which
    /// is what lets `R-112` own the words and borrow the directory.
    pub(crate) struct RankDirectory {
        /// Absolute rank at the start of each superblock.
        supers: Vec<u32>,
        /// Rank **within** the superblock at the start of each block.
        blocks: Vec<u16>,
    }

    impl RankDirectory {
        /// One pass over the words. `total` is returned by [`Self::total`].
        pub(crate) fn new(words: &[u64]) -> Self {
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

        /// The output slot of the cell at `bit`: how many bits are set before it.
        ///
        /// Two directory reads and at most `BLOCK_WORDS − 1` word folds, so the
        /// cost is bounded by a constant that does not depend on the bitmap size.
        #[inline]
        pub(crate) fn rank(&self, words: &[u64], bit: usize) -> u32 {
            let word = bit >> 6;
            let block = word / BLOCK_WORDS;
            let mut r = self.supers[block / BLOCKS_PER_SUPER] + u32::from(self.blocks[block]);
            for w in &words[(block * BLOCK_WORDS)..word] {
                r += w.count_ones();
            }
            r + (words[word] & ((1u64 << (bit & 63)) - 1)).count_ones()
        }

        /// How many bits the bitmap has set.
        pub(crate) fn total(&self, words: &[u64]) -> u32 {
            let last = self.blocks.len() - 1;
            let mut r = self.supers[last / BLOCKS_PER_SUPER] + u32::from(self.blocks[last]);
            for w in &words[(last * BLOCK_WORDS)..] {
                r += w.count_ones();
            }
            r
        }

        /// Bytes the directory occupies.
        pub(crate) fn bytes(&self) -> usize {
            self.supers.len() * size_of::<u32>() + self.blocks.len() * size_of::<u16>()
        }

        /// Words a query at `bit` folds. Deterministic, and the O(1) evidence.
        #[inline]
        fn words_folded(bit: usize) -> usize {
            let word = bit >> 6;
            word - (word / BLOCK_WORDS) * BLOCK_WORDS
        }
    }

    /// The same directory **one level short**: superblocks only.
    ///
    /// The vacuity control. C3's comparand is the running prefix counter, so a
    /// comparator that cannot fail proves nothing; this is a directory with the
    /// block level missing whose query still jumps to the block boundary, and it
    /// is wrong for every cell that has set bits earlier in its own superblock
    /// but outside its own block.
    struct ShortDirectory {
        supers: Vec<u32>,
    }

    impl ShortDirectory {
        fn new(words: &[u64]) -> Self {
            let block_count = words.len().div_ceil(BLOCK_WORDS);
            let super_count = block_count.div_ceil(BLOCKS_PER_SUPER);
            let mut supers = Vec::with_capacity(super_count);
            let mut total = 0u32;
            for block in 0..block_count {
                if block % BLOCKS_PER_SUPER == 0 {
                    supers.push(total);
                }
                let lo = block * BLOCK_WORDS;
                let hi = (lo + BLOCK_WORDS).min(words.len());
                for word in &words[lo..hi] {
                    total += word.count_ones();
                }
            }
            Self { supers }
        }

        #[inline]
        fn rank(&self, words: &[u64], bit: usize) -> u32 {
            let word = bit >> 6;
            let block = word / BLOCK_WORDS;
            let mut r = self.supers[block / BLOCKS_PER_SUPER];
            for w in &words[(block * BLOCK_WORDS)..word] {
                r += w.count_ones();
            }
            r + (words[word] & ((1u64 << (bit & 63)) - 1)).count_ones()
        }
    }

    /// The registered comparand: a full prefix scan from word zero, per query.
    #[inline]
    fn rank_scan(words: &[u64], bit: usize) -> u32 {
        let word = bit >> 6;
        let mut r = 0u32;
        for w in &words[..word] {
            r += w.count_ones();
        }
        r + (words[word] & ((1u64 << (bit & 63)) - 1)).count_ones()
    }

    // ─── the offset/compaction stage ───────────────────────────────────────

    /// The two passes the directory replaces, and nothing else.
    #[derive(Default)]
    struct Offset {
        /// One `u32` slot per cell, `u32::MAX` for "no slot".
        slots: Vec<u32>,
    }

    impl Offset {
        /// The dense slot table, cleared and refilled — the shape
        /// `marching_cubes/mod.rs:250-251` runs on every extract.
        fn prepare(&mut self, cells: usize) {
            self.slots.clear();
            self.slots.resize(cells, u32::MAX);
            black_box(&self.slots);
        }

        /// The prefix-sum pass: `dual.rs:489-497`'s set-bit walk with a running
        /// counter, ascending `x` so the slot order is the crate's own.
        fn scan(&mut self, bits: &Bitmaps) -> u32 {
            let mut next = 0u32;
            let rows = bits.cells_x * bits.cells_x;
            for row in 0..rows {
                let lin = row * bits.cells_x;
                let base = row * bits.cell_words;
                for w in 0..bits.cell_words {
                    let mut word = bits.active[base + w];
                    while word != 0 {
                        let x = w * 64 + word.trailing_zeros() as usize;
                        word &= word - 1;
                        self.slots[lin + x] = next;
                        next += 1;
                    }
                }
            }
            black_box(&self.slots);
            next
        }
    }

    // ─── counting ──────────────────────────────────────────────────────────

    /// Cycles and instructions from one window.
    #[derive(Clone, Copy, Default)]
    struct Counted {
        cycles: f64,
        instructions: f64,
        nanos: f64,
    }

    /// One counter window over `inner` repetitions, divided by `inner`.
    ///
    /// Every `perf_event` system call is outside the counted region. Windows are
    /// **siblings, never nested**: Zen 3 has six general-purpose counters and
    /// `Probe` opens exactly six, so a nested window multiplexes and
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

    /// The median of a set of readings, taken **per quantity**.
    fn median(values: &mut [f64]) -> f64 {
        values.sort_by(f64::total_cmp);
        values[values.len() / 2]
    }

    /// Median cycles, instructions and nanoseconds of one quantity.
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

    // ─── the query set ─────────────────────────────────────────────────────

    /// `QUERIES` deterministic cell-bit indices, splitmix64 so the set is the
    /// same in every run and every build.
    fn query_set(cell_bits: usize) -> Vec<usize> {
        let mut state = 0x9E37_79B9_7F4A_7C15u64;
        (0..QUERIES)
            .map(|_| {
                state = state.wrapping_add(0x9E37_79B9_7F4A_7C15);
                let mut z = state;
                z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
                z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
                z ^= z >> 31;
                (z % cell_bits as u64) as usize
            })
            .collect()
    }

    // ─── one row ───────────────────────────────────────────────────────────

    struct Row {
        field: &'static str,
        resolution: u32,
        cells: usize,
        active_cells: u32,
        active_fraction: f64,
        bitmap_words: usize,
        bitmap_bytes: usize,
        cell_words: usize,
        sample_bit_row: usize,
        directory_bytes: usize,
        overhead_fraction: f64,
        slot_table_bytes: usize,
        crate_edge_table_bytes: usize,
        // C1
        offset_share: f64,
        offset_share_dc: f64,
        offset_share_scan_only: f64,
        offset_share_instructions: f64,
        cycles_extract_mc: f64,
        cycles_extract_dc: f64,
        cycles_offset: f64,
        cycles_offset_prepare: f64,
        cycles_offset_scan: f64,
        instructions_extract_mc: f64,
        instructions_offset: f64,
        // C2
        max_words_scanned_rank: usize,
        max_words_scanned_scan: usize,
        count_ones_per_query_rank: f64,
        count_ones_per_query_scan: f64,
        target_feature_popcnt: bool,
        ns_per_query_rank: f64,
        ns_per_query_scan: f64,
        ns_per_query_sequential: f64,
        instructions_per_query_rank: f64,
        instructions_per_query_scan: f64,
        instructions_per_query_sequential: f64,
        query_ratio_ns: f64,
        query_ratio_instructions: f64,
        // C3
        slots_equal: u32,
        cells_checked: usize,
        rank_equal_all_bits: usize,
        short_directory_mismatches: usize,
        bitmap_matches_scalar: bool,
        pad_bits_set: usize,
        directory_total_equals_popcount: bool,
        // provenance
        ghz: f64,
        inner_extract: usize,
        inner_offset: usize,
        inner_query: usize,
        c1_holds: bool,
        c1_holds_dc: bool,
        c1_holds_instructions: bool,
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
        let shape = RuntimeShape3::new([n; 3]).expect("the fixture fits u32");
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
        let directory = RankDirectory::new(&bits.active);
        let short = ShortDirectory::new(&bits.active);
        let active_cells: u32 = bits.active.iter().map(|w| w.count_ones()).sum();

        // ── the bitmap is the right population, checked scalar-wise ──────────
        //
        // `active_word` is four fused word operations and `cell_mask` discards a
        // whole word at 65 and 129 samples per axis. If either is wrong, C3 is
        // an equality over the wrong set of cells, so the eight-corner scalar
        // classification is run once per row and compared bit for bit.
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
                    let bit = (z * cells_x + y) * bits.cell_words * 64 + x;
                    if bits.get(bit) != scalar_active {
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
            active_cells > 0,
            "{field} {n}^3: no active cells, so C3 would be an equality over an empty set"
        );
        assert!(
            directory.total(&bits.active) == active_cells,
            "{field} {n}^3: the directory's total disagrees with the bitmap's popcount"
        );

        // `cell_mask`'s own check: every bit past the end of a row must be
        // clear, or the rank of every later cell counts a cell that does not
        // exist. The scalar sweep above cannot see these bits, because it
        // enumerates cells and they are not cells.
        let pad_bits_set = (0..bits.cell_bits)
            .filter(|&bit| !bits.is_cell(bit) && bits.get(bit))
            .count();
        assert_eq!(
            pad_bits_set, 0,
            "{field} {n}^3: {pad_bits_set} bits past the end of a cell row are set, so \
             `cell_mask` is wrong and every rank after the first short row is too"
        );

        // ── C3, and its control ──────────────────────────────────────────────
        //
        // The comparand is the running prefix counter itself, which is what makes
        // `short_directory_mismatches` load-bearing rather than decorative.
        let mut prefix = 0u32;
        let mut slots_equal = 0u32;
        let mut rank_equal_all_bits = 0usize;
        let mut short_directory_mismatches = 0usize;
        for bit in 0..bits.cell_bits {
            let answer = directory.rank(&bits.active, bit);
            if answer == prefix {
                rank_equal_all_bits += 1;
            }
            if short.rank(&bits.active, bit) != prefix {
                short_directory_mismatches += 1;
            }
            if bits.get(bit) {
                if answer == prefix {
                    slots_equal += 1;
                }
                prefix += 1;
            }
        }
        assert_eq!(
            slots_equal,
            active_cells,
            "{field} {n}^3: the directory's slot differs from the prefix sum's on \
             {} of {active_cells} active cells",
            active_cells - slots_equal
        );
        assert!(
            short_directory_mismatches > 0,
            "{field} {n}^3: a directory built one level short still answered every cell \
             correctly, so C3's comparator cannot fail and the equality is vacuous"
        );

        // ── the query set, and the word counts that carry C2's O(1) ──────────
        let queries = query_set(bits.cell_bits);
        let max_words_scanned_rank = queries
            .iter()
            .map(|&bit| RankDirectory::words_folded(bit))
            .max()
            .expect("QUERIES is non-zero");
        let max_words_scanned_scan = queries.iter().map(|&bit| bit >> 6).max().unwrap_or(0);
        assert!(
            max_words_scanned_rank < BLOCK_WORDS,
            "{field} {n}^3: a rank query folded {max_words_scanned_rank} words, which is not \
             bounded by the block width and therefore not O(1)"
        );

        // How many `count_ones` calls each arm makes per query. This build emits
        // **no `popcnt`** — there is no `.cargo/config.toml` and no `target-cpu`,
        // so the `x86-64` baseline is in force and `u64::count_ones()` lowers to
        // the twelve-instruction SWAR sequence. `objdump -d` on this bench's own
        // binary greps zero `popcnt` and 91 occurrences of the SWAR constant
        // `0x3333333333333333`. So a paper figure that assumes a one-cycle
        // popcount — pasta's rank block among them — is pricing something an
        // order of magnitude cheaper than what runs here, and the per-query call
        // count is the arithmetic that says by how much.
        let count_ones_per_query_rank = queries
            .iter()
            .map(|&bit| RankDirectory::words_folded(bit) as f64)
            .sum::<f64>()
            / QUERIES as f64
            + 1.0;
        let count_ones_per_query_scan =
            queries.iter().map(|&bit| (bit >> 6) as f64).sum::<f64>() / QUERIES as f64 + 1.0;

        // ── the windows: shipped extraction, the offset pass, the queries ────
        let mut mc = MarchingCubes::<f32>::new();
        let mut dc = DualContouring::<f32>::new();
        let mut out = MeshBuffer::<f32>::new();
        let mut offset = Offset::default();
        offset.prepare(cells);

        for _ in 0..WARMUP {
            out.reset();
            mc.extract(sdf, &shape, origin, cell_size, &mut out)
                .expect("extraction");
            black_box(&out);
            out.reset();
            dc.extract(sdf, &shape, origin, cell_size, &mut out)
                .expect("extraction");
            black_box(&out);
        }

        let inner_extract = {
            let mut probe_free = || {
                out.reset();
                mc.extract(sdf, &shape, origin, cell_size, &mut out)
                    .expect("extraction");
                black_box(&out);
            };
            choose_inner(&mut probe_free)
        };
        let inner_offset = choose_inner(|| {
            offset.prepare(cells);
            offset.scan(&bits);
        });
        let inner_query = choose_inner(|| {
            let mut acc = 0u32;
            for &bit in &queries {
                acc = acc.wrapping_add(directory.rank(&bits.active, bit));
            }
            black_box(acc);
        });

        let mut probe = Probe::open();
        let mut extract_mc = Vec::with_capacity(REPS);
        let mut extract_dc = Vec::with_capacity(REPS);
        let mut prepare = Vec::with_capacity(REPS);
        let mut scan = Vec::with_capacity(REPS);
        let mut query_rank = Vec::with_capacity(REPS);
        let mut query_scan = Vec::with_capacity(REPS);
        let mut query_seq = Vec::with_capacity(REPS);

        for _ in 0..REPS {
            extract_mc.push(window(&mut probe, inner_extract, || {
                out.reset();
                mc.extract(sdf, &shape, origin, cell_size, &mut out)
                    .expect("extraction");
                black_box(&out);
            }));
            extract_dc.push(window(&mut probe, inner_extract, || {
                out.reset();
                dc.extract(sdf, &shape, origin, cell_size, &mut out)
                    .expect("extraction");
                black_box(&out);
            }));
            prepare.push(window(&mut probe, inner_offset, || offset.prepare(cells)));
            scan.push(window(&mut probe, inner_offset, || {
                black_box(offset.scan(&bits));
            }));
            query_rank.push(window(&mut probe, inner_query, || {
                let mut acc = 0u32;
                for &bit in &queries {
                    acc = acc.wrapping_add(directory.rank(&bits.active, bit));
                }
                black_box(acc);
            }));
            query_scan.push(window(&mut probe, inner_query, || {
                let mut acc = 0u32;
                for &bit in &queries {
                    acc = acc.wrapping_add(rank_scan(&bits.active, bit));
                }
                black_box(acc);
            }));
            query_seq.push(window(&mut probe, inner_offset, || {
                black_box(offset.scan(&bits));
            }));
        }

        let mc_counts = median_of(&extract_mc);
        let dc_counts = median_of(&extract_dc);
        let prepare_counts = median_of(&prepare);
        let scan_counts = median_of(&scan);
        let rank_counts = median_of(&query_rank);
        let scan_q_counts = median_of(&query_scan);
        let seq_counts = median_of(&query_seq);

        let cycles_offset = prepare_counts.cycles + scan_counts.cycles;
        let instructions_offset = prepare_counts.instructions + scan_counts.instructions;
        let offset_share = cycles_offset / mc_counts.cycles;
        let offset_share_dc = cycles_offset / dc_counts.cycles;
        let offset_share_scan_only = scan_counts.cycles / mc_counts.cycles;
        let offset_share_instructions = instructions_offset / mc_counts.instructions;

        let per_query = 1.0 / QUERIES as f64;
        // The sequential comparand's "query" is a cell: the running-counter pass
        // answers every cell once, so its per-query cost is amortised over all of
        // them. That is the comparand the shipped path actually uses.
        let per_cell = 1.0 / cells as f64;

        let bitmap_bytes = bits.active.len() * size_of::<u64>();
        let directory_bytes = directory.bytes();
        let overhead_fraction = directory_bytes as f64 / bitmap_bytes as f64;

        let c1_holds = offset_share >= SHARE_BAR;
        let c2_holds = overhead_fraction <= OVERHEAD_BAR && max_words_scanned_rank < BLOCK_WORDS;
        let c3_holds = slots_equal == active_cells
            && rank_equal_all_bits == bits.cell_bits
            && short_directory_mismatches > 0
            && bitmap_matches_scalar;

        Row {
            field,
            resolution: n,
            cells,
            active_cells,
            active_fraction: f64::from(active_cells) / cells as f64,
            bitmap_words: bits.active.len(),
            bitmap_bytes,
            cell_words: bits.cell_words,
            sample_bit_row: bits.bit_row,
            directory_bytes,
            overhead_fraction,
            slot_table_bytes: cells * size_of::<u32>(),
            crate_edge_table_bytes: samples * 3 * size_of::<u32>(),
            offset_share,
            offset_share_dc,
            offset_share_scan_only,
            offset_share_instructions,
            cycles_extract_mc: mc_counts.cycles,
            cycles_extract_dc: dc_counts.cycles,
            cycles_offset,
            cycles_offset_prepare: prepare_counts.cycles,
            cycles_offset_scan: scan_counts.cycles,
            instructions_extract_mc: mc_counts.instructions,
            instructions_offset,
            max_words_scanned_rank,
            max_words_scanned_scan,
            count_ones_per_query_rank,
            count_ones_per_query_scan,
            target_feature_popcnt: cfg!(target_feature = "popcnt"),
            ns_per_query_rank: rank_counts.nanos * per_query,
            ns_per_query_scan: scan_q_counts.nanos * per_query,
            ns_per_query_sequential: seq_counts.nanos * per_cell,
            instructions_per_query_rank: rank_counts.instructions * per_query,
            instructions_per_query_scan: scan_q_counts.instructions * per_query,
            instructions_per_query_sequential: seq_counts.instructions * per_cell,
            query_ratio_ns: scan_q_counts.nanos / rank_counts.nanos.max(1.0),
            query_ratio_instructions: scan_q_counts.instructions
                / rank_counts.instructions.max(1.0),
            slots_equal,
            cells_checked: bits.cell_bits,
            rank_equal_all_bits,
            short_directory_mismatches,
            bitmap_matches_scalar,
            pad_bits_set,
            directory_total_equals_popcount: true,
            ghz: mc_counts.cycles / mc_counts.nanos,
            inner_extract,
            inner_offset,
            inner_query,
            c1_holds,
            c1_holds_dc: offset_share_dc >= SHARE_BAR,
            c1_holds_instructions: offset_share_instructions >= SHARE_BAR,
            c2_holds,
            c3_holds,
        }
    }

    /// Eight reference fields × three resolutions, `f32`.
    ///
    /// No `scalar` column is registered and none is added: the offset pass is
    /// integer work over a bitmap and does not change with the field's precision,
    /// while the *denominator* does — so an `f64` arm would only move C1's
    /// denominator and would double the fixture to say so.
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

        let clearing = rows.iter().filter(|r| r.c1_holds).count();
        let clearing_dc = rows.iter().filter(|r| r.c1_holds_dc).count();
        let worst_overhead = rows
            .iter()
            .map(|r| r.overhead_fraction)
            .fold(0.0f64, f64::max);
        let best_share = rows.iter().map(|r| r.offset_share).fold(0.0f64, f64::max);

        println!(
            "P-107: offset share clears {SHARE_BAR:.2} on {clearing} of {} rows against \
             marching_cubes and {clearing_dc} against dual_contouring; the largest share \
             measured is {best_share:.4}",
            rows.len()
        );
        println!(
            "P-107: directory overhead is at most {:.4} against C2's bar of {OVERHEAD_BAR:.4}; \
             a rank query folds at most {} words at every resolution",
            worst_overhead,
            rows.iter()
                .map(|r| r.max_words_scanned_rank)
                .max()
                .unwrap_or(0)
        );
        println!(
            "P-107: target_feature_popcnt is {}, so a `count_ones` is the SWAR sequence and not \
             an instruction; a rank query makes {:.1} of them and the offset stage makes none",
            cfg!(target_feature = "popcnt"),
            rows[0].count_ones_per_query_rank
        );

        for row in &rows {
            run.record(&[
                ("field", row.field.to_string()),
                ("resolution", row.resolution.to_string()),
                ("cells", row.cells.to_string()),
                ("active_cells", row.active_cells.to_string()),
                ("active_fraction", format!("{:.6}", row.active_fraction)),
                ("bitmap_words", row.bitmap_words.to_string()),
                ("bitmap_bytes", row.bitmap_bytes.to_string()),
                ("cell_words_per_row", row.cell_words.to_string()),
                ("sample_bit_row", row.sample_bit_row.to_string()),
                ("directory_bytes", row.directory_bytes.to_string()),
                ("overhead_fraction", format!("{:.6}", row.overhead_fraction)),
                ("slot_table_bytes", row.slot_table_bytes.to_string()),
                (
                    "crate_edge_table_bytes",
                    row.crate_edge_table_bytes.to_string(),
                ),
                ("offset_share", format!("{:.6}", row.offset_share)),
                ("offset_share_dc", format!("{:.6}", row.offset_share_dc)),
                (
                    "offset_share_scan_only",
                    format!("{:.6}", row.offset_share_scan_only),
                ),
                (
                    "offset_share_instructions",
                    format!("{:.6}", row.offset_share_instructions),
                ),
                ("cycles_extract_mc", format!("{:.1}", row.cycles_extract_mc)),
                ("cycles_extract_dc", format!("{:.1}", row.cycles_extract_dc)),
                ("cycles_offset", format!("{:.1}", row.cycles_offset)),
                (
                    "cycles_offset_prepare",
                    format!("{:.1}", row.cycles_offset_prepare),
                ),
                (
                    "cycles_offset_scan",
                    format!("{:.1}", row.cycles_offset_scan),
                ),
                (
                    "instructions_extract_mc",
                    format!("{:.1}", row.instructions_extract_mc),
                ),
                (
                    "instructions_offset",
                    format!("{:.1}", row.instructions_offset),
                ),
                (
                    "max_words_scanned_rank",
                    row.max_words_scanned_rank.to_string(),
                ),
                (
                    "max_words_scanned_scan",
                    row.max_words_scanned_scan.to_string(),
                ),
                (
                    "count_ones_per_query_rank",
                    format!("{:.3}", row.count_ones_per_query_rank),
                ),
                (
                    "count_ones_per_query_scan",
                    format!("{:.3}", row.count_ones_per_query_scan),
                ),
                (
                    "target_feature_popcnt",
                    row.target_feature_popcnt.to_string(),
                ),
                ("ns_per_query_rank", format!("{:.4}", row.ns_per_query_rank)),
                ("ns_per_query_scan", format!("{:.4}", row.ns_per_query_scan)),
                (
                    "ns_per_query_sequential",
                    format!("{:.4}", row.ns_per_query_sequential),
                ),
                (
                    "instructions_per_query_rank",
                    format!("{:.3}", row.instructions_per_query_rank),
                ),
                (
                    "instructions_per_query_scan",
                    format!("{:.3}", row.instructions_per_query_scan),
                ),
                (
                    "instructions_per_query_sequential",
                    format!("{:.3}", row.instructions_per_query_sequential),
                ),
                ("query_ratio_ns", format!("{:.3}", row.query_ratio_ns)),
                (
                    "query_ratio_instructions",
                    format!("{:.3}", row.query_ratio_instructions),
                ),
                ("slots_equal", row.slots_equal.to_string()),
                ("cells_checked", row.cells_checked.to_string()),
                ("rank_equal_all_bits", row.rank_equal_all_bits.to_string()),
                (
                    "short_directory_mismatches",
                    row.short_directory_mismatches.to_string(),
                ),
                (
                    "bitmap_matches_scalar",
                    row.bitmap_matches_scalar.to_string(),
                ),
                ("pad_bits_set", row.pad_bits_set.to_string()),
                (
                    "directory_total_equals_popcount",
                    row.directory_total_equals_popcount.to_string(),
                ),
                ("ghz", format!("{:.4}", row.ghz)),
                ("inner_extract", row.inner_extract.to_string()),
                ("inner_offset", row.inner_offset.to_string()),
                ("inner_query", row.inner_query.to_string()),
                ("c1_holds", row.c1_holds.to_string()),
                ("c1_holds_dc", row.c1_holds_dc.to_string()),
                (
                    "c1_holds_instructions",
                    row.c1_holds_instructions.to_string(),
                ),
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

    let prereg = isomesh::experiment!("P-107");

    #[cfg(target_os = "linux")]
    common::experiment::run(prereg, experiment::run);

    // `experiment_p12`'s precedent. C1 is a share of extraction cycles and C2 is
    // a query cost; a nanosecond on a governed CPU cannot carry either (`M-281`),
    // and a recorded zero would be a fabricated share — which is what `✗51` cost.
    #[cfg(not(target_os = "linux"))]
    {
        eprintln!(
            "{} measures a share of extraction with hardware performance counters, and this \
             platform has no `perf_event_open`. There is no clock substitute.",
            prereg.id
        );
        std::process::exit(1);
    }
}
