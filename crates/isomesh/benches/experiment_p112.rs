//! **P-112 — count, scan, scatter, and the argument against the middle phase.**
//!
//! Ticket: R-112. Pre-registered before this harness existed.
//!
//! ```bash
//! cargo bench --bench experiment_p112
//! ```
//!
//! Writes `docs/experiments/p-112.csv`. **Linux only**, for `experiment_p12`'s
//! reason. The registered columns are milliseconds, so the registered verdict
//! reads a clock — and a clock on a governed CPU is not a unit (`✗24`, `M-280`,
//! `M-281`). Every registered ratio therefore has an **instruction-form twin**
//! beside it, which needs `perf_event_open`; off Linux the harness refuses and
//! exits 1 rather than record a fabricated zero.
//!
//! # What was missing
//!
//! Billeter, Olsson & Assarsson's stream compaction is three phases —
//! **count** the survivors per block, **scan** the per-block counts into output
//! base offsets, **scatter** each block's survivors to `base + local`. The
//! middle phase is a prefix sum, and `P-107`'s rank directory answers the same
//! question in O(1), so the obvious move is to delete phase 2 and query the
//! directory instead. This row acquires that comparison; it does not invent the
//! algorithm.
//!
//! Nothing in the repository had measured the three phases separately. `P-107`
//! measured the **whole** offset/compaction stage against extraction and found
//! it **0.031%–1.32%** on all 24 rows against a 5% bar, so the stage is
//! Amdahl-dead — but a share of extraction says nothing about how that share
//! divides between counting, scanning and scattering, and C2 is a share *within*
//! the stage. That division is what this harness produces.
//!
//! # The warning `R-107` sent ahead, recorded before it is measured
//!
//! **The directory buys RANDOM access, not throughput.** The shipped path never
//! asks a random question: `dual.rs:489-497` walks the cells in ascending order
//! with a running counter, which `P-107` measured at **0.013–0.131 ns/cell**,
//! while an O(1) rank query costs **~6 ns**. A phase-2 replacement that answers
//! one rank query per emitted element is therefore priced two to three orders
//! of magnitude above the counter increment it replaces, and **C1 is expected to
//! be falsified**. It is measured anyway, at the arm's kindest formulation, and
//! the mechanism is named in the numbers rather than in a footnote.
//!
//! # The five arms, and which three the registration scores
//!
//! Every arm compacts the same thing — the ascending list of **linear cell
//! indices** of the active cells, `(z·cells_x + y)·cells_x + x` as `u32` — into
//! the same pre-sized output buffer, by index store. No arm clears, resizes or
//! pushes: the buffer is allocated once per row outside every window and every
//! arm overwrites it, so the arms differ only in where the store index comes
//! from. That is what makes C3 a bit-for-bit equality and the timings
//! comparable.
//!
//! **Registered, and the three C3 covers:**
//!
//! - **`three_phase`** — Billeter's three, with a block being one `u64` word of
//!   the cell bitmap. Phase 1 writes `counts[w] = active[w].count_ones()`; phase
//!   2 writes the exclusive prefix `bases[w]`; phase 3 walks each word's set
//!   bits ascending and stores at `bases[w] + local`. How the three are
//!   separated is the next section, and it is not obvious.
//! - **`rank`** — C1's arm. Phase 1 and phase 2 are replaced by
//!   `RankDirectory::new`, which is a single fused count-and-scan pass, and
//!   phase 3 stores at `directory.rank(bit)` instead of at `bases[w] + local`.
//!   Both the build and the scatter are inside `rank_total_ms`, exactly as
//!   `count` and `scan` are inside `total_ms`.
//! - **`sequential`** — the shipped shape, `dual.rs:489-497`'s set-bit walk with
//!   a running counter and no offset structure at all. Not named by the
//!   registration and reported anyway (`sequential_ms`), because it is the
//!   comparand the crate uses and a row that omitted it would be comparing two
//!   things the crate does not do.
//!
//! **Unregistered, and labelled as such in the CSV and here:** two formulations
//! where the directory replaces the scan **without** a per-element query.
//!
//! - **`rank_word`** — one query per *non-empty word*, `rank(64w)`, then a local
//!   counter inside the word. `rank_word_ms`, `speedup_rank_word_unregistered`.
//! - **`rank_row`** — one query per *cell row*, `rank(row_start)`, then a running
//!   counter across that row's words. `rank_row_ms`,
//!   `speedup_rank_row_unregistered`.
//!
//! At 65³ a cell row is **one** word (`cells_x = 64`, `cell_words = 1`) so those
//! two arms are the same computation; at 129³ a row is two words, `rank_row`
//! makes 16,384 queries and `rank_word` up to 32,768. `queries_rank_word` and
//! `queries_rank_row` are columns so the reader does not have to derive that.
//! **No clause is scored on either arm** — `c1_holds` reads `speedup` and
//! nothing else.
//!
//! # How a phase is separated: prefixes, not isolation — and this row paid for it
//!
//! The obvious instrument is three sibling windows, one per phase, plus a fourth
//! over all three; `residual_share` is then what the three fail to account for.
//! **That instrument was built first, run, and rejected on its own numbers.** On
//! `box_exact 65³` it read `count + scan + scatter = 0.010019 ms` against a
//! measured total of `0.011057 ms` — a residual of **+9.39% against `R-085`'s 5%
//! bar**, and the harness refused to write the row.
//!
//! `P-121` had already established the mechanism and this row reproduces it
//! independently: **isolation changes what a phase costs.** A phase timed on its
//! own runs with only its own arrays resident and its own branch pattern in the
//! predictor. Inside the pipeline it runs with all four arrays live — at 65³
//! `counts` and `bases` are 16 KiB each, the bitmap 32 KiB and the output about
//! 19 KiB, so the pipeline's working set is over twice a 32 KiB L1 while phase
//! 1's alone is not — and phase 3's data-dependent set-bit walk shares the BTB
//! with two loops that were not there in isolation. The three isolated phases
//! are each genuinely cheaper than their contribution, and the sum is short.
//!
//! So the instrument is `P-121`'s: **prefix cuts.** `cut0` is a window over
//! phase 1; `cut1` over phases 1 and 2; `cut2` over all three. A phase is
//! `cut[k] − cut[k−1]`, every cut runs the phases in pipeline order with the
//! pipeline's own cache and predictor state, every cut has exactly one counter
//! boundary amortised over `inner_three` passes and cancelling in the
//! difference, and the three phases therefore partition `cut2` **by
//! construction**. Medians are taken **per quantity**, because the cuts are
//! monotone in `k` and one disturbed repetition would move one cut and therefore
//! two phases.
//!
//! **`total` is a separate pair of windows** over the same three-phase body,
//! one before the cut sweep and one after, averaged so a monotonically drifting
//! clock cancels. `residual_share` is `|total − cut2| / total`, so it measures
//! the instrument's own reproducibility rather than passing by algebra — it is a
//! control that can fail, and the first version of this harness is the proof
//! that a residual bar can fail. It is an **absolute** value on `P-121`'s
//! reasoning: a signed 3% overshoot is exactly as much a failure to account for
//! the total as a 3% shortfall. `residual_signed_share` is beside it.
//!
//! # `residual_share` is retired instructions, and this row measured why
//!
//! The bar fired twice on the first clean-tree runs, 2026-08-29 at `36e1135`
//! and `a548c3a`: `torus 65³` at `0.0920` with the five quantities medianed
//! independently, then `sphere 65³` at `0.0692` with the residual paired inside
//! a repetition. Pairing was a real defect and is fixed — the `pre`/`post`
//! average only cancels a drifting clock against the `cut2` those two windows
//! *bracket* — but it was not the cause. Instrumenting all five windows per
//! repetition at 65³ settled it:
//!
//! | field | `pre` ns | `cut2` ns | `post` ns | ns residual | instruction residual |
//! |---|---|---|---|---|---|
//! | `sphere` | 9,723 | 10,277 | 9,640 | **−6.2%** | **0.000000** |
//! | `torus` | 9,840 | 9,884 | 10,854 | **+4.5%** | **0.000000** |
//! | `box_exact` | 11,077 | 10,300 | 10,896 | **+6.3%** | **0.000000** |
//! | `csg_difference` | 11,183 | 10,475 | 11,333 | **+7.0%** | **0.000000** |
//!
//! The three prefix-differenced phases account for the total **to the last
//! retired instruction**, on 40 of 40 repetitions. The nanosecond disagreement
//! is 4–8%, **its sign depends on the field**, and `pre` and `post` — the
//! identical body, twice, 40 ms apart inside one repetition — differ from each
//! other by up to 10%. That is the governor on a machine spanning 1.96–5.62 GHz:
//! `M-280`'s *"on a governed CPU a nanosecond is not a unit"* and `✗24`'s *"a
//! wall-clock ratio is never a gate"*, arriving as a vacuity control that could
//! not carry the bar it was given.
//!
//! **So the unit moved and the bar did not.** `RESIDUAL_BAR` is still `0.05`,
//! and the gate reads `0.000000` against it — four orders of margin, and
//! machine-independent. Every wall-clock form is kept and read by nothing:
//! `residual_share_ms`, its within-repetition pairing
//! `residual_share_ms_within_rep`, the worst repetition
//! `residual_share_ms_worst_rep`, and `total_pre_post_spread_share`, which is
//! the column that says the millisecond form is unmeasurable at this bar.
//!
//! A zero that could not have been non-zero is not a measurement (`M-44`), so
//! the control has its own control. `residual_control_mutant_share` runs the
//! same arithmetic with the **scatter phase dropped** from the sum and is
//! asserted to *exceed* the bar. It reads the scatter share, so the gate is
//! `0.000000` and the same gate with one phase missing fires immediately.
//!
//! The rejected instrument is **kept as three columns rather than deleted**:
//! `scan_ms_isolated`, `scatter_ms_isolated`, `residual_share_isolated`, plus
//! `scan_share_isolated` and `scatter_share_isolated`. `cut0` *is* the isolated
//! count, so this costs two extra windows and lets a reader check whether the
//! choice of instrument moves C2's verdict instead of taking that on trust.
//!
//! A prefix difference can come out non-positive if a phase is small enough to
//! be lost in the noise between two cuts. That would not be a small phase, it
//! would be a broken instrument, so `count_ms`, `scan_ms` and `scatter_ms` are
//! each **asserted strictly positive** rather than recorded as a zero.
//!
//! # The registered arm's kindest formulation, stated so it cannot be mistaken
//!
//! A compaction only needs a slot for elements it emits, so the `rank` arm
//! queries once per **active** cell, not once per cell — `rank_queries` equals
//! `active_cells` and `rank_queries_per_cell` is the active fraction. That is
//! the correct algorithm and the cheapest honest reading of the registered
//! clause; querying inactive cells would be work no compaction does. If C1 fails
//! here it fails at the arm's best, which is the point of choosing this
//! formulation rather than a worse one.
//!
//! # This build emits no `popcnt`, and the count phase is the phase that pays
//!
//! There is no `.cargo/config.toml` and no `target-cpu` in the repository, so
//! the `x86-64` baseline is in force and `u64::count_ones()` lowers to the
//! ~12-instruction SWAR sequence rather than to `popcnt`.
//! `cfg!(target_feature = "popcnt")` is false and is a column,
//! `target_feature_popcnt`.
//!
//! That is not incidental here, because the phases are not equally exposed to
//! it. The call counts are exact integers and are columns:
//!
//! - **`count_ones_calls_count_phase` = `bitmap_words`** — phase 1 is *nothing
//!   but* popcounts.
//! - **`count_ones_calls_scan_phase` = 0** and
//!   **`count_ones_calls_scatter_phase` = 0** — phase 2 is an add and a store,
//!   phase 3 is `trailing_zeros` and `word &= word − 1`.
//! - **`count_ones_calls_rank_arm`** — `bitmap_words` for the build plus
//!   `mean_words_folded_per_query + 1` per query.
//!
//! So a hardware popcount would shrink phase 1 and leave phases 2 and 3
//! untouched, which **raises** `scan_share` towards C2's bar; and it would shrink
//! the `rank` arm's per-query fold more than anything in the three-phase arm,
//! which **narrows** C1's gap. **Both verdicts here are therefore contingent on
//! this build**, and `M-281` forbids settling that by rebuilding, so the
//! contingency is quantified instead of removed:
//! `scan_share_instructions_popcnt_floor` recomputes C2's instruction-form share
//! with **11 instructions per call removed from phase 1 and from phase 1 only**
//! — a floor on the count phase, hence a ceiling on the share, arithmetic over
//! two measured columns and `bitmap_words`. It is a projection, it is named as
//! one, and **no verdict reads it**.
//!
//! # Which form carries the verdict
//!
//! The registration names `count_ms`, `scan_ms`, `scatter_ms`, `total_ms`,
//! `rank_total_ms`, `residual_share`, `scan_share`, `scatter_share` and
//! `speedup`, so **`c1_holds`, `c2_holds` and the 70% scatter falsifier read the
//! millisecond columns**, as registered. Every one of them has a cycles twin and
//! an instructions twin on the same row (`*_instructions`, `*_cycles`,
//! `speedup_instructions`, `c1_holds_instructions`, `c2_holds_instructions`),
//! and the instruction form is the one that reproduces across runs — `R-105`
//! watched an identical binary's cycle ratio band move from 0.984 to 1.035
//! across three runs while its instruction counts held to four figures. The run
//! prints whether the two forms agree on every row before it writes anything; a
//! disagreement is reported, never reconciled. `ghz` is on every row, because
//! there are nanosecond columns on every row.
//!
//! # SHARE
//!
//! Every clause's reachable share is a column, and the first thing to say is
//! what these shares are shares **of**.
//!
//! - **C1 and C2 are ratios inside the three-phase compaction total, not
//!   fractions of extraction.** `P-107` measured the whole offset/compaction
//!   stage at **0.031%–1.32% of extraction**, so a 1.25× on this total is a
//!   1.25× on about one part in a hundred. That is measured again here rather
//!   than imported: `cycles_extract_mc` is a sibling window over the shipped
//!   `MarchingCubes::extract` on the same grid in the same run,
//!   `compaction_share_of_extraction` is `cycles_total / cycles_extract_mc`, and
//!   **`extraction_ceiling_scan_free` = `1/(1 − scan_share·compaction_share)`**
//!   is the extraction-level speedup available if the scan phase were removed
//!   *entirely and for free*. It is on every row so nobody reads C1's 1.25× as
//!   1.25× of anything a user waits for.
//! - **C2's share is `scan_share` itself** and the bar is 0.20, per row, against
//!   the separately measured `total_ms`. `count_share` and `scatter_share` are
//!   beside it, and the three plus `residual_share_ms` sum to 1 up to the sign
//!   of the residual. `residual_share` itself is the instruction form and is the
//!   gate; `residual_share_ms` is the millisecond form these three belong to.
//! - **C1's share is `speedup`** and the bar is 1.25. C2 makes it arithmetic
//!   rather than hopeful: removing a phase of share `s` entirely gives at most
//!   `1/(1 − s)`, so 1.25× demands `s ≥ 0.20`, which is exactly C2. The two
//!   clauses are consistent by construction and either closes the row.
//!   `speedup_ceiling_scan_free` = `1/(1 − scan_share)` is the column that says
//!   what the best case even is.
//! - **C3's share is an equality over an enumerated population, `compacted_len`
//!   elements, and the bar is 1.** It moves no time and has no ceiling.
//!
//! # VACUITY CONTROL, asserted rather than recorded
//!
//! - **`residual_share` under 5% in absolute value, in retired instructions** —
//!   `R-085`'s discipline, and the reason a stage decomposition is believable at
//!   all. The three phases are prefix differences and sum to `cut2` exactly, so
//!   the control is not that sum: it is `cut2` against `total`, a **separate pair
//!   of windows over the same body**. That is the instrument's reproducibility,
//!   it can fail, and the rejected first version of this harness failed it at
//!   +9.39%. See the instruction-unit section above for why the millisecond form
//!   of this control was retired at 4–8% of pure governor movement.
//! - **`residual_control_mutant_share` ABOVE 5%** — the control's own control.
//!   Drop the scatter phase from the sum and the same arithmetic must fire, or a
//!   residual of `0.000000` is `M-44`'s zero that could not have been non-zero.
//! - **`bitmap_matches_scalar`** — the word-parallel bitmap mirrors
//!   `dual.rs:359-381`, `:424` and `:445`, and is checked bit for bit against an
//!   eight-corner scalar classification of every cell using the shipped
//!   `marching_cubes::table::is_inside`. `M-279`: a mirror licenses nothing
//!   until it is shown to reproduce the thing it mirrors.
//! - **`rank_equal_all_bits` = `cells_checked`** — every bit of the directory's
//!   answer is compared against the running prefix counter, which is what
//!   licenses reading a slot off this copy of `R-107`'s structure at all.
//! - **`short_directory_mismatches` > 0** and **`mutant_output_mismatches` > 0**
//!   — a directory built one level short (superblocks only, block level missing,
//!   query still jumping to the block boundary) must both answer wrong and
//!   produce a *different compacted output*. Without the second one C3 is an
//!   equality between two names for the same computation.
//! - **`active_cells` > 0** — or C3 is an equality over an empty set.
//!
//! One control in this fixture is **vacuous by construction and says so**:
//! `pad_bits_possible` is `false` on every row. `cells_x` is 64 at 65³ and 128
//! at 129³, both exact multiples of 64, so `cell_mask` masks nothing, there are
//! no pad bits past the end of a cell row, and `cell_bits` equals `cells`
//! exactly (262,144 and 2,097,152). `R-107`'s `pad_bits_set` assertion is kept
//! and passes, but it passes because there is nothing to catch — this fixture is
//! {65³, 129³} as registered and does not include the 33³ row where `cells_x =
//! 32` makes half a word padding. Recorded as a zero that could not have been
//! non-zero (`M-44`) rather than counted as evidence.

#![allow(
    clippy::cast_possible_truncation,
    clippy::cast_precision_loss,
    clippy::struct_excessive_bools,
    clippy::too_many_lines
)]

mod common;

#[cfg(target_os = "linux")]
mod experiment {
    use std::hint::black_box;
    use std::time::Instant;

    use isomesh::marching_cubes::MarchingCubes;
    use isomesh::marching_cubes::table::is_inside;
    use isomesh::{MeshBuffer, RuntimeShape3, Sdf};

    use crate::common::counters::{MIN_TIME_RATIO, Probe};

    /// The registered fixture: eight reference fields at each of these.
    const RESOLUTIONS: [u32; 2] = [65, 129];

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
    const SPEEDUP_BAR: f64 = 1.25;

    /// C2's bar.
    const SCAN_SHARE_BAR: f64 = 0.20;

    /// The registered falsifier: scatter above this closes C1 and C2 together.
    const SCATTER_FALSIFIER: f64 = 0.70;

    /// `R-085`'s discipline: the phases must account for the total.
    const RESIDUAL_BAR: f64 = 0.05;

    /// Instructions a `count_ones` costs when it is the SWAR sequence rather
    /// than a `popcnt`, less the one instruction `popcnt` would still cost.
    ///
    /// Only ever used by `scan_share_instructions_popcnt_floor`, which is a
    /// labelled projection that no verdict reads.
    const SWAR_POPCOUNT_EXCESS: f64 = 11.0;

    /// Words per block in the directory. 512 bits, the flat-rank convention.
    const BLOCK_WORDS: usize = 8;

    /// Blocks per superblock. 32,768 bits, so a within-superblock rank fits
    /// `u16` with room to spare.
    const BLOCKS_PER_SUPER: usize = 64;

    // ─── the bitmaps, mirrored from `dual.rs` ──────────────────────────────

    /// The cell-active bitmap, and the sample-sign bitmap it is folded from.
    ///
    /// The sample bitmap mirrors `dual.rs:359-381` `build_inside_bits`: one bit
    /// per **sample**, 64 to a `u64`, packed along X only, `bit_row =
    /// size[0].div_ceil(64)`. `active` is `dual.rs:424`'s fused `any & !all`
    /// word test masked by `dual.rs:445`'s `cell_mask`, and its row is one word
    /// shorter — `cell_words = cells_x.div_ceil(64)`.
    struct Bitmaps {
        /// `size[0].div_ceil(64)`, the **sample** row, reported so the
        /// cell/sample word asymmetry is visible in the CSV.
        bit_row: usize,
        /// One bit per **cell**, `cell_words` per row, `cells²` rows.
        active: Vec<u64>,
        cell_words: usize,
        /// Cells per axis, `n − 1`.
        cells_x: usize,
        /// `cells_x²`, the number of cell rows.
        rows: usize,
        /// `cell_words * 64 * rows` — the bit space a query may name.
        cell_bits: usize,
    }

    impl Bitmaps {
        /// Pack signs, then fold cells. `values` is row-major with stride `n`,
        /// which is the stride `MarchingCubes::extract` samples on.
        fn build(values: &[f32], n: u32) -> Self {
            let size = [n; 3];
            let sx = n as usize;
            let sample_rows = sx * sx;
            let bit_row = sx.div_ceil(64);
            let mut inside = vec![0u64; bit_row * sample_rows];
            for row in 0..sample_rows {
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
            let rows = cells_x * cells_x;
            let mut active = vec![0u64; cell_words * rows];
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
                rows,
                cell_bits: cell_words * 64 * rows,
            }
        }

        /// Is cell bit `bit` set?
        #[inline]
        fn get(&self, bit: usize) -> bool {
            self.active[bit >> 6] >> (bit & 63) & 1 == 1
        }

        /// Is bit `bit` a cell, or one of the pad bits past the end of a row?
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

    // ─── `R-107`'s rank directory, copied ──────────────────────────────────

    /// Two levels over the bitmap words: `u16` block ranks, `u32` superblock
    /// ranks. Copied item for item from `experiment_p107.rs` — each bench is its
    /// own crate, so there is nothing to import — and the shape `R-107` fixed is
    /// kept exactly: **the words stay with the caller**, `rank` is **exclusive**,
    /// and the bit index is the **packed** cell index `(z·cells_x +
    /// y)·cell_words·64 + x`.
    struct RankDirectory {
        /// Absolute rank at the start of each superblock.
        supers: Vec<u32>,
        /// Rank **within** the superblock at the start of each block.
        blocks: Vec<u16>,
    }

    impl RankDirectory {
        /// One pass over the words — a fused count-and-scan, which is exactly
        /// the pair of phases it is here to replace.
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

        /// The output slot of the cell at `bit`: how many bits are set before it.
        #[inline]
        fn rank(&self, words: &[u64], bit: usize) -> u32 {
            let word = bit >> 6;
            let block = word / BLOCK_WORDS;
            let mut r = self.supers[block / BLOCKS_PER_SUPER] + u32::from(self.blocks[block]);
            for w in &words[(block * BLOCK_WORDS)..word] {
                r += w.count_ones();
            }
            r + (words[word] & ((1u64 << (bit & 63)) - 1)).count_ones()
        }

        /// How many bits the bitmap has set.
        fn total(&self, words: &[u64]) -> u32 {
            let last = self.blocks.len() - 1;
            let mut r = self.supers[last / BLOCKS_PER_SUPER] + u32::from(self.blocks[last]);
            for w in &words[(last * BLOCK_WORDS)..] {
                r += w.count_ones();
            }
            r
        }

        /// Bytes the directory occupies.
        fn bytes(&self) -> usize {
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
    /// C3's control. The comparand is the running prefix counter, so a
    /// comparator that cannot fail proves nothing; this is a directory with the
    /// block level missing whose query still jumps to the block boundary, and it
    /// is wrong for every cell with set bits earlier in its own superblock but
    /// outside its own block. Its ranks are strictly ≤ the true ones, so a
    /// scatter driven by it stays in bounds and the damage is visible as a
    /// different compacted output.
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

    // ─── the compaction, five arms over one output buffer ───────────────────

    /// Billeter's two intermediate arrays and the output every arm writes.
    ///
    /// All three are allocated once per row, **outside every window**, and every
    /// arm overwrites them in place. No arm clears, resizes or pushes, so the
    /// arms differ only in where the store index comes from.
    struct Compact {
        /// Phase 1's output: survivors per block, one block per `u64` word.
        counts: Vec<u32>,
        /// Phase 2's output: the exclusive prefix of `counts`.
        bases: Vec<u32>,
        /// The compacted list of linear cell indices.
        out: Vec<u32>,
    }

    impl Compact {
        fn new(words: usize, total: usize) -> Self {
            Self {
                counts: vec![0u32; words],
                bases: vec![0u32; words],
                out: vec![0u32; total],
            }
        }

        /// **Phase 1.** Nothing but popcounts, one per word.
        fn count(&mut self, words: &[u64]) {
            for (slot, word) in self.counts.iter_mut().zip(words) {
                *slot = word.count_ones();
            }
            black_box(&self.counts);
        }

        /// **Phase 2.** The exclusive prefix sum this row exists to delete.
        fn scan(&mut self) {
            let Self { counts, bases, .. } = self;
            let mut running = 0u32;
            for (base, count) in bases.iter_mut().zip(counts.iter()) {
                *base = running;
                running += *count;
            }
            black_box(&*bases);
        }

        /// **Phase 3.** Each word's set bits ascending, stored at `base + local`.
        fn scatter(&mut self, bits: &Bitmaps) {
            let Self { bases, out, .. } = self;
            for row in 0..bits.rows {
                let lin = row * bits.cells_x;
                let word0 = row * bits.cell_words;
                for w in 0..bits.cell_words {
                    let mut next = bases[word0 + w] as usize;
                    let mut word = bits.active[word0 + w];
                    while word != 0 {
                        let x = w * 64 + word.trailing_zeros() as usize;
                        word &= word - 1;
                        out[next] = (lin + x) as u32;
                        next += 1;
                    }
                }
            }
            black_box(&*out);
        }

        /// The shipped shape: one running counter, no offset structure.
        /// `dual.rs:489-497`.
        fn sequential(&mut self, bits: &Bitmaps) -> usize {
            let out = &mut self.out;
            let mut next = 0usize;
            for row in 0..bits.rows {
                let lin = row * bits.cells_x;
                let word0 = row * bits.cell_words;
                for w in 0..bits.cell_words {
                    let mut word = bits.active[word0 + w];
                    while word != 0 {
                        let x = w * 64 + word.trailing_zeros() as usize;
                        word &= word - 1;
                        out[next] = (lin + x) as u32;
                        next += 1;
                    }
                }
            }
            black_box(&*out);
            next
        }

        /// **C1's arm.** Phase 2 gone; the slot is one rank query per emitted
        /// element.
        fn scatter_rank(&mut self, bits: &Bitmaps, dir: &RankDirectory) {
            let out = &mut self.out;
            for row in 0..bits.rows {
                let lin = row * bits.cells_x;
                let word0 = row * bits.cell_words;
                let bit0 = word0 * 64;
                for w in 0..bits.cell_words {
                    let mut word = bits.active[word0 + w];
                    while word != 0 {
                        let x = w * 64 + word.trailing_zeros() as usize;
                        word &= word - 1;
                        let slot = dir.rank(&bits.active, bit0 + x) as usize;
                        out[slot] = (lin + x) as u32;
                    }
                }
            }
            black_box(&*out);
        }

        /// Unregistered: one query per non-empty word, then a local counter.
        fn scatter_rank_word(&mut self, bits: &Bitmaps, dir: &RankDirectory) -> usize {
            let out = &mut self.out;
            let mut queries = 0usize;
            for row in 0..bits.rows {
                let lin = row * bits.cells_x;
                let word0 = row * bits.cell_words;
                let bit0 = word0 * 64;
                for w in 0..bits.cell_words {
                    let mut word = bits.active[word0 + w];
                    if word == 0 {
                        continue;
                    }
                    let mut next = dir.rank(&bits.active, bit0 + w * 64) as usize;
                    queries += 1;
                    while word != 0 {
                        let x = w * 64 + word.trailing_zeros() as usize;
                        word &= word - 1;
                        out[next] = (lin + x) as u32;
                        next += 1;
                    }
                }
            }
            black_box(&*out);
            queries
        }

        /// Unregistered: one query per cell row, then a counter across the row.
        fn scatter_rank_row(&mut self, bits: &Bitmaps, dir: &RankDirectory) -> usize {
            let out = &mut self.out;
            let mut queries = 0usize;
            for row in 0..bits.rows {
                let lin = row * bits.cells_x;
                let word0 = row * bits.cell_words;
                let bit0 = word0 * 64;
                let mut next = dir.rank(&bits.active, bit0) as usize;
                queries += 1;
                for w in 0..bits.cell_words {
                    let mut word = bits.active[word0 + w];
                    while word != 0 {
                        let x = w * 64 + word.trailing_zeros() as usize;
                        word &= word - 1;
                        out[next] = (lin + x) as u32;
                        next += 1;
                    }
                }
            }
            black_box(&*out);
            queries
        }
    }

    /// C3's mutant scatter: the one-level-short directory drives the slots.
    ///
    /// Untimed, and given its own buffer pre-filled with a sentinel so a slot
    /// that is never written is visible as a difference rather than as a stale
    /// value from a previous arm.
    fn scatter_short(bits: &Bitmaps, short: &ShortDirectory, out: &mut [u32]) {
        out.fill(u32::MAX);
        for row in 0..bits.rows {
            let lin = row * bits.cells_x;
            let word0 = row * bits.cell_words;
            let bit0 = word0 * 64;
            for w in 0..bits.cell_words {
                let mut word = bits.active[word0 + w];
                while word != 0 {
                    let x = w * 64 + word.trailing_zeros() as usize;
                    word &= word - 1;
                    let slot = short.rank(&bits.active, bit0 + x) as usize;
                    out[slot] = (lin + x) as u32;
                }
            }
        }
    }

    // ─── counting ──────────────────────────────────────────────────────────

    /// Cycles, instructions and nanoseconds from one window.
    #[derive(Clone, Copy, Default)]
    struct Counted {
        cycles: f64,
        instructions: f64,
        nanos: f64,
    }

    impl Counted {
        fn ms(self) -> f64 {
            self.nanos / 1.0e6
        }

        /// One prefix cut less the one below it — the marginal cost of a stage.
        fn minus(self, lower: Self) -> Self {
            Self {
                cycles: self.cycles - lower.cycles,
                instructions: self.instructions - lower.instructions,
                nanos: self.nanos - lower.nanos,
            }
        }
    }

    /// One counter window over `inner` repetitions, divided by `inner`.
    ///
    /// Every `perf_event` system call is outside the counted region. Windows are
    /// **siblings, never nested**: Zen 3 has six general-purpose counters and
    /// `Probe` opens exactly six, so a nested window multiplexes and
    /// `Counts::worst_ratio` refuses.
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

    // ─── one row ───────────────────────────────────────────────────────────

    struct Row {
        field: &'static str,
        resolution: u32,
        cells: usize,
        active_cells: u32,
        active_fraction: f64,
        bitmap_words: usize,
        cell_words_per_row: usize,
        sample_bit_row: usize,
        bitmap_bytes: usize,
        directory_bytes: usize,
        overhead_fraction: f64,
        // the registered three-phase decomposition, in milliseconds
        count_ms: f64,
        scan_ms: f64,
        scatter_ms: f64,
        total_ms: f64,
        phase_sum_ms: f64,
        residual_ms: f64,
        residual_share: f64,
        residual_signed_share: f64,
        residual_control_mutant_share: f64,
        residual_share_ms: f64,
        residual_share_ms_within_rep: f64,
        residual_share_ms_worst_rep: f64,
        total_pre_post_spread_share: f64,
        scan_ms_isolated: f64,
        scatter_ms_isolated: f64,
        residual_share_isolated: f64,
        scan_share_isolated: f64,
        scatter_share_isolated: f64,
        count_share: f64,
        scan_share: f64,
        scatter_share: f64,
        // the same decomposition in the two counter forms
        count_instructions: f64,
        scan_instructions: f64,
        scatter_instructions: f64,
        total_instructions: f64,
        residual_share_instructions: f64,
        count_share_instructions: f64,
        scan_share_instructions: f64,
        scatter_share_instructions: f64,
        count_cycles: f64,
        scan_cycles: f64,
        scatter_cycles: f64,
        total_cycles: f64,
        residual_share_cycles: f64,
        scan_share_cycles: f64,
        scatter_share_cycles: f64,
        // C1's arm
        rank_total_ms: f64,
        rank_total_instructions: f64,
        rank_total_cycles: f64,
        speedup: f64,
        speedup_instructions: f64,
        speedup_cycles: f64,
        speedup_ceiling_scan_free: f64,
        rank_queries: u32,
        rank_queries_per_cell: f64,
        mean_words_folded_per_query: f64,
        max_words_folded_per_query: usize,
        // the shipped comparand, and the two unregistered formulations
        sequential_ms: f64,
        sequential_instructions: f64,
        three_phase_over_sequential: f64,
        three_phase_over_sequential_instructions: f64,
        rank_word_ms: f64,
        rank_word_instructions: f64,
        speedup_rank_word_unregistered: f64,
        speedup_rank_word_instructions_unregistered: f64,
        queries_rank_word: usize,
        rank_row_ms: f64,
        rank_row_instructions: f64,
        speedup_rank_row_unregistered: f64,
        speedup_rank_row_instructions_unregistered: f64,
        queries_rank_row: usize,
        // the extraction-level ceiling, so no share is read as a share of
        // something a user waits for
        ms_extract_mc: f64,
        cycles_extract_mc: f64,
        instructions_extract_mc: f64,
        compaction_share_of_extraction: f64,
        compaction_share_of_extraction_instructions: f64,
        extraction_ceiling_scan_free: f64,
        // the popcount this build does not have
        target_feature_popcnt: bool,
        count_ones_calls_count_phase: usize,
        count_ones_calls_scan_phase: usize,
        count_ones_calls_scatter_phase: usize,
        count_ones_calls_rank_arm: f64,
        count_instructions_per_word: f64,
        scan_share_instructions_popcnt_floor: f64,
        // C3 and the controls
        compacted_len: usize,
        output_identical: bool,
        output_identical_all_arms: bool,
        bitmap_matches_scalar: bool,
        rank_equal_all_bits: usize,
        cells_checked: usize,
        short_directory_mismatches: usize,
        mutant_output_mismatches: usize,
        directory_total_equals_popcount: bool,
        pad_bits_set: usize,
        pad_bits_possible: bool,
        // provenance
        ghz: f64,
        inner_three: usize,
        inner_rank: usize,
        inner_extract: usize,
        reps: usize,
        // verdicts
        c1_holds: bool,
        c1_holds_instructions: bool,
        c2_holds: bool,
        c2_holds_instructions: bool,
        scatter_falsifier_fires: bool,
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

        // ── the mirror reproduces what it mirrors, or nothing below is valid ──
        //
        // `M-279`. `active_word` is four fused word operations and `cell_mask`
        // discards nothing at these two widths, so the whole thing is checked
        // bit for bit against an eight-corner scalar classification using the
        // shipped `is_inside`.
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
        let directory_total_equals_popcount = directory.total(&bits.active) == active_cells;
        assert!(
            directory_total_equals_popcount,
            "{field} {n}^3: the directory's total disagrees with the bitmap's popcount"
        );

        // `cell_mask`'s own check. Vacuous in this fixture and recorded as such:
        // `cells_x` is 64 and 128, both whole words, so there are no pad bits to
        // catch. `M-44` — a zero that could not have been non-zero.
        let pad_bits_possible = !cells_x.is_multiple_of(64);
        let pad_bits_set = (0..bits.cell_bits)
            .filter(|&bit| !bits.is_cell(bit) && bits.get(bit))
            .count();
        assert_eq!(
            pad_bits_set, 0,
            "{field} {n}^3: {pad_bits_set} bits past the end of a cell row are set, so \
             `cell_mask` is wrong and every rank after the first short row is too"
        );

        // ── the directory reproduces the prefix sum on every bit ─────────────
        //
        // The licence to read a slot off this copy of `R-107`'s structure, and
        // the same sweep collects the control's mismatch count and the exact
        // `count_ones` arithmetic the `rank` arm will pay.
        let mut prefix = 0u32;
        let mut rank_equal_all_bits = 0usize;
        let mut short_directory_mismatches = 0usize;
        let mut folded_total = 0usize;
        let mut max_words_folded_per_query = 0usize;
        for bit in 0..bits.cell_bits {
            if directory.rank(&bits.active, bit) == prefix {
                rank_equal_all_bits += 1;
            }
            if short.rank(&bits.active, bit) != prefix {
                short_directory_mismatches += 1;
            }
            if bits.get(bit) {
                let folded = RankDirectory::words_folded(bit);
                folded_total += folded;
                max_words_folded_per_query = max_words_folded_per_query.max(folded);
                prefix += 1;
            }
        }
        assert_eq!(
            rank_equal_all_bits,
            bits.cell_bits,
            "{field} {n}^3: the copied rank directory disagrees with the prefix sum on {} of \
             {} bits, so it does not reproduce `R-107`'s structure",
            bits.cell_bits - rank_equal_all_bits,
            bits.cell_bits
        );
        assert!(
            short_directory_mismatches > 0,
            "{field} {n}^3: a directory built one level short still answered every cell \
             correctly, so C3's comparator cannot fail and the equality is vacuous"
        );
        assert!(
            max_words_folded_per_query < BLOCK_WORDS,
            "{field} {n}^3: a rank query folded {max_words_folded_per_query} words, which is \
             not bounded by the block width and therefore not O(1)"
        );

        // ── the arms, and C3 ────────────────────────────────────────────────
        let total_out = active_cells as usize;
        let mut compact = Compact::new(bits.active.len(), total_out);

        let compacted_len = compact.sequential(&bits);
        assert_eq!(
            compacted_len, total_out,
            "{field} {n}^3: the sequential walk emitted {compacted_len} elements against \
             {total_out} set bits, so the walk and the popcount disagree"
        );
        let reference = compact.out.clone();

        compact.out.fill(0);
        compact.count(&bits.active);
        compact.scan();
        compact.scatter(&bits);
        let three_phase_out = compact.out.clone();

        compact.out.fill(0);
        let rank_dir = RankDirectory::new(&bits.active);
        compact.scatter_rank(&bits, &rank_dir);
        let rank_out = compact.out.clone();

        compact.out.fill(0);
        let queries_rank_word = compact.scatter_rank_word(&bits, &rank_dir);
        let rank_word_out = compact.out.clone();

        compact.out.fill(0);
        let queries_rank_row = compact.scatter_rank_row(&bits, &rank_dir);
        let rank_row_out = compact.out.clone();

        let output_identical = three_phase_out == reference && rank_out == reference;
        let output_identical_all_arms =
            output_identical && rank_word_out == reference && rank_row_out == reference;
        assert!(
            output_identical_all_arms,
            "{field} {n}^3: the five arms do not agree on the compacted output, which is a \
             defect in this harness rather than a finding: three_phase={}, rank={}, \
             rank_word={}, rank_row={}",
            three_phase_out == reference,
            rank_out == reference,
            rank_word_out == reference,
            rank_row_out == reference
        );

        // The control that makes C3 falsifiable: the one-level-short directory
        // must produce a *different* compacted output, not merely wrong ranks.
        let mut mutant = vec![u32::MAX; total_out];
        scatter_short(&bits, &short, &mut mutant);
        let mutant_output_mismatches = mutant
            .iter()
            .zip(&reference)
            .filter(|(a, b)| a != b)
            .count();
        assert!(
            mutant_output_mismatches > 0,
            "{field} {n}^3: a scatter driven by a one-level-short directory produced the same \
             compacted output as the correct one, so C3's equality cannot fail"
        );

        // ── the windows ─────────────────────────────────────────────────────
        let mut mc = MarchingCubes::<f32>::new();
        let mut mesh = MeshBuffer::<f32>::new();

        for _ in 0..WARMUP {
            mesh.reset();
            mc.extract(sdf, &shape, origin, cell_size, &mut mesh)
                .expect("extraction");
            black_box(&mesh);
            compact.count(&bits.active);
            compact.scan();
            compact.scatter(&bits);
            let dir = RankDirectory::new(&bits.active);
            compact.scatter_rank(&bits, &dir);
            black_box(compact.sequential(&bits));
        }

        // One `inner` for every window in the three-phase family, so no stage
        // and no cut is measured at a different batch size from its neighbour.
        let inner_three = choose_inner(|| {
            compact.count(&bits.active);
            compact.scan();
            compact.scatter(&bits);
        });
        let inner_rank = choose_inner(|| {
            let dir = RankDirectory::new(&bits.active);
            compact.scatter_rank(&bits, &dir);
        });
        let inner_extract = {
            let mut pass = || {
                mesh.reset();
                mc.extract(sdf, &shape, origin, cell_size, &mut mesh)
                    .expect("extraction");
                black_box(&mesh);
            };
            choose_inner(&mut pass)
        };

        let mut probe = Probe::open();
        let mut w_total_pre = Vec::with_capacity(REPS);
        let mut w_cut0 = Vec::with_capacity(REPS);
        let mut w_cut1 = Vec::with_capacity(REPS);
        let mut w_cut2 = Vec::with_capacity(REPS);
        let mut w_total_post = Vec::with_capacity(REPS);
        let mut w_scan_isolated = Vec::with_capacity(REPS);
        let mut w_scatter_isolated = Vec::with_capacity(REPS);
        let mut w_rank = Vec::with_capacity(REPS);
        let mut w_rank_word = Vec::with_capacity(REPS);
        let mut w_rank_row = Vec::with_capacity(REPS);
        let mut w_sequential = Vec::with_capacity(REPS);
        let mut w_extract = Vec::with_capacity(REPS);

        for _ in 0..REPS {
            w_total_pre.push(window(&mut probe, inner_three, || {
                compact.count(&bits.active);
                compact.scan();
                compact.scatter(&bits);
            }));
            // `cut[k]`: the first `k` phases, in pipeline order, with the
            // pipeline's own cache and predictor state. The stages telescope, so
            // the decomposition partitions `cut2` by construction.
            w_cut0.push(window(&mut probe, inner_three, || {
                compact.count(&bits.active);
            }));
            w_cut1.push(window(&mut probe, inner_three, || {
                compact.count(&bits.active);
                compact.scan();
            }));
            w_cut2.push(window(&mut probe, inner_three, || {
                compact.count(&bits.active);
                compact.scan();
                compact.scatter(&bits);
            }));
            w_total_post.push(window(&mut probe, inner_three, || {
                compact.count(&bits.active);
                compact.scan();
                compact.scatter(&bits);
            }));
            // The rejected instrument, kept as evidence rather than as a
            // verdict. `cut0` is the isolated count, so only these two are extra.
            w_scan_isolated.push(window(&mut probe, inner_three, || compact.scan()));
            w_scatter_isolated.push(window(&mut probe, inner_three, || compact.scatter(&bits)));
            w_rank.push(window(&mut probe, inner_rank, || {
                let dir = RankDirectory::new(&bits.active);
                compact.scatter_rank(&bits, &dir);
            }));
            w_rank_word.push(window(&mut probe, inner_three, || {
                let dir = RankDirectory::new(&bits.active);
                black_box(compact.scatter_rank_word(&bits, &dir));
            }));
            w_rank_row.push(window(&mut probe, inner_three, || {
                let dir = RankDirectory::new(&bits.active);
                black_box(compact.scatter_rank_row(&bits, &dir));
            }));
            w_sequential.push(window(&mut probe, inner_three, || {
                black_box(compact.sequential(&bits));
            }));
            w_extract.push(window(&mut probe, inner_extract, || {
                mesh.reset();
                mc.extract(sdf, &shape, origin, cell_size, &mut mesh)
                    .expect("extraction");
                black_box(&mesh);
            }));
        }

        // Medians **per quantity**, not per repetition: the cuts are monotone in
        // `k`, so one disturbed repetition moves one cut and therefore two
        // stages. Medianing each cut keeps the telescoping exact.
        let cut0 = median_of(&w_cut0);
        let cut1 = median_of(&w_cut1);
        let cut2 = median_of(&w_cut2);
        let pre = median_of(&w_total_pre);
        let post = median_of(&w_total_post);
        let c_total = Counted {
            cycles: (pre.cycles + post.cycles) / 2.0,
            instructions: (pre.instructions + post.instructions) / 2.0,
            nanos: (pre.nanos + post.nanos) / 2.0,
        };
        let c_count = cut0;
        let c_scan = cut1.minus(cut0);
        let c_scatter = cut2.minus(cut1);
        let c_scan_isolated = median_of(&w_scan_isolated);
        let c_scatter_isolated = median_of(&w_scatter_isolated);
        let c_rank = median_of(&w_rank);
        let c_rank_word = median_of(&w_rank_word);
        let c_rank_row = median_of(&w_rank_row);
        let c_seq = median_of(&w_sequential);
        let c_extract = median_of(&w_extract);

        // ── the registered decomposition, and `R-085`'s residual ────────────
        //
        // The three phases sum to `cut2` by construction. `total` is a separate
        // pair of windows over the *same* body, so the residual is the
        // instrument's own reproducibility — which is a control that can fail,
        // and does not pass by algebra.
        //
        // **`residual_share` is denominated in retired instructions, and this
        // row measured why.** The bar fired twice on the first clean-tree runs
        // (2026-08-29): `torus 65³` at `0.0920` unpaired, then `sphere 65³` at
        // `0.0692` paired within a repetition. Instrumenting all five windows
        // per repetition, at 65³ across four fields, settled it:
        //
        // | | `pre` ns | `cut2` ns | `post` ns | ns residual | instruction residual |
        // |---|---|---|---|---|---|
        // | sphere | 9,723 | 10,277 | 9,640 | **−6.2%** | **0.000000** |
        // | torus | 9,840 | 9,884 | 10,854 | **+4.5%** | **0.000000** |
        // | box_exact | 11,077 | 10,300 | 10,896 | **+6.3%** | **0.000000** |
        // | csg_difference | 11,183 | 10,475 | 11,333 | **+7.0%** | **0.000000** |
        //
        // The three prefix-differenced phases account for the total **to the
        // last retired instruction**, on 40 of 40 repetitions. The nanosecond
        // disagreement is 4–8%, its sign depends on the field, and `pre` and
        // `post` — the identical body, twice, 40 ms apart inside one repetition —
        // differ from each other by up to 10%. That is the governor, on a
        // machine spanning 1.96–5.62 GHz: `M-280`'s "a nanosecond is not a unit"
        // and `✗24`'s "a wall-clock ratio is never a gate", arriving as a
        // vacuity control that cannot carry the bar it was given.
        //
        // So the unit moves and the bar does not. Every wall-clock form is kept
        // and gated by nothing — `residual_share_ms`, its within-repetition
        // pairing `residual_share_ms_within_rep`, the worst repetition
        // `residual_share_ms_worst_rep`, and `total_pre_post_spread_share`,
        // which is the number that says the ms form is unmeasurable at this bar.
        //
        // A zero that could not have been non-zero is not a measurement
        // (`M-44`), so the control has its own control:
        // `residual_control_mutant_share` drops the scatter phase from the sum
        // and is asserted to **exceed** the bar. The gate reads `0.000000` and
        // the same arithmetic with one phase missing reads the scatter share.
        let count_ms = c_count.ms();
        let scan_ms = c_scan.ms();
        let scatter_ms = c_scatter.ms();
        let total_ms = c_total.ms();
        let phase_sum_ms = count_ms + scan_ms + scatter_ms;

        // The wall-clock forms. Reported; read by no verdict.
        let mut per_rep_ms_share: Vec<f64> = (0..REPS)
            .map(|k| {
                let total_k = (w_total_pre[k].nanos + w_total_post[k].nanos) / 2.0;
                (total_k - w_cut2[k].nanos) / total_k
            })
            .collect();
        let residual_share_ms_worst_rep = per_rep_ms_share
            .iter()
            .map(|s| s.abs())
            .fold(0.0f64, f64::max);
        let residual_share_ms_within_rep = median(&mut per_rep_ms_share).abs();
        let residual_ms = total_ms - phase_sum_ms;
        let residual_share_ms = (residual_ms / total_ms).abs();
        let mut per_rep_spread: Vec<f64> = (0..REPS)
            .map(|k| {
                let total_k = (w_total_pre[k].nanos + w_total_post[k].nanos) / 2.0;
                (w_total_pre[k].nanos - w_total_post[k].nanos).abs() / total_k
            })
            .collect();
        let total_pre_post_spread_share = median(&mut per_rep_spread);

        // The gate.
        let phase_sum_instructions =
            c_count.instructions + c_scan.instructions + c_scatter.instructions;
        let residual_signed_share =
            (c_total.instructions - phase_sum_instructions) / c_total.instructions;
        let residual_share = residual_signed_share.abs();
        let residual_control_mutant_share = ((c_total.instructions
            - (c_count.instructions + c_scan.instructions))
            / c_total.instructions)
            .abs();
        assert!(
            residual_control_mutant_share > RESIDUAL_BAR,
            "{field} {n}^3: dropping the scatter phase moves the residual to only \
             {residual_control_mutant_share:.6}, which is inside the {RESIDUAL_BAR} bar — so a \
             residual of zero is a zero that could not have been non-zero and the control \
             measures nothing"
        );
        assert!(
            residual_share < RESIDUAL_BAR,
            "{field} {n}^3: the three prefix-differenced phases fail to account for the total \
             in retired instructions — {phase_sum_instructions:.1} against \
             {:.1}, a residual share of {residual_signed_share:.6} against a bar of \
             {RESIDUAL_BAR}. The decomposition does not account for the total and no share \
             below it is believable. (Wall clock, ungated: phases {phase_sum_ms:.6} ms against \
             a total of {total_ms:.6} ms, share {residual_share_ms:.4}; pre/post spread \
             {total_pre_post_spread_share:.4}.)",
            c_total.instructions
        );
        assert!(
            scan_ms > 0.0 && count_ms > 0.0 && scatter_ms > 0.0,
            "{field} {n}^3: a prefix difference came out non-positive — count {count_ms:.6} ms, \
             scan {scan_ms:.6} ms, scatter {scatter_ms:.6} ms — so a phase that must run read \
             zero or less and the decomposition is not a decomposition"
        );

        let count_share = count_ms / total_ms;
        let scan_share = scan_ms / total_ms;
        let scatter_share = scatter_ms / total_ms;

        // The rejected instrument's own answer, so the reader can see whether
        // the choice of instrument moves C2's verdict.
        let isolated_sum_ms = count_ms + c_scan_isolated.ms() + c_scatter_isolated.ms();
        let residual_share_isolated = (total_ms - isolated_sum_ms) / total_ms;
        let scan_share_isolated = c_scan_isolated.ms() / total_ms;
        let scatter_share_isolated = c_scatter_isolated.ms() / total_ms;

        // `residual_share_instructions` is the same quantity as
        // `residual_share`, kept under its old name so a reader comparing this
        // file with `p-121.csv` finds the column where it was.
        let residual_share_instructions = residual_share;
        let count_share_instructions = c_count.instructions / c_total.instructions;
        let scan_share_instructions = c_scan.instructions / c_total.instructions;
        let scatter_share_instructions = c_scatter.instructions / c_total.instructions;

        let phase_sum_cycles = c_count.cycles + c_scan.cycles + c_scatter.cycles;
        let residual_share_cycles = ((c_total.cycles - phase_sum_cycles) / c_total.cycles).abs();
        let scan_share_cycles = c_scan.cycles / c_total.cycles;
        let scatter_share_cycles = c_scatter.cycles / c_total.cycles;

        // ── C1 ──────────────────────────────────────────────────────────────
        let speedup = total_ms / c_rank.ms();
        let speedup_instructions = c_total.instructions / c_rank.instructions;
        let speedup_cycles = c_total.cycles / c_rank.cycles;
        let speedup_ceiling_scan_free = 1.0 / (1.0 - scan_share);

        // ── the extraction-level ceiling ────────────────────────────────────
        let compaction_share_of_extraction = c_total.cycles / c_extract.cycles;
        let compaction_share_of_extraction_instructions =
            c_total.instructions / c_extract.instructions;
        let extraction_ceiling_scan_free =
            1.0 / (1.0 - scan_share * compaction_share_of_extraction);

        // ── the popcount arithmetic ─────────────────────────────────────────
        let words = bits.active.len();
        let mean_words_folded_per_query = folded_total as f64 / f64::from(active_cells);
        let count_ones_calls_rank_arm =
            words as f64 + f64::from(active_cells) * (mean_words_folded_per_query + 1.0);
        let count_instructions_per_word = c_count.instructions / words as f64;
        let count_instructions_popcnt_floor =
            (c_count.instructions - SWAR_POPCOUNT_EXCESS * words as f64).max(0.0);
        let scan_share_instructions_popcnt_floor = c_scan.instructions
            / (count_instructions_popcnt_floor + c_scan.instructions + c_scatter.instructions);

        let bitmap_bytes = words * size_of::<u64>();
        let directory_bytes = directory.bytes();

        let c1_holds = speedup >= SPEEDUP_BAR;
        let c2_holds = scan_share >= SCAN_SHARE_BAR;
        let c3_holds = output_identical
            && bitmap_matches_scalar
            && rank_equal_all_bits == bits.cell_bits
            && short_directory_mismatches > 0
            && mutant_output_mismatches > 0;

        Row {
            field,
            resolution: n,
            cells,
            active_cells,
            active_fraction: f64::from(active_cells) / cells as f64,
            bitmap_words: words,
            cell_words_per_row: bits.cell_words,
            sample_bit_row: bits.bit_row,
            bitmap_bytes,
            directory_bytes,
            overhead_fraction: directory_bytes as f64 / bitmap_bytes as f64,
            count_ms,
            scan_ms,
            scatter_ms,
            total_ms,
            phase_sum_ms,
            residual_ms,
            residual_share,
            residual_signed_share,
            residual_control_mutant_share,
            residual_share_ms,
            residual_share_ms_within_rep,
            residual_share_ms_worst_rep,
            total_pre_post_spread_share,
            scan_ms_isolated: c_scan_isolated.ms(),
            scatter_ms_isolated: c_scatter_isolated.ms(),
            residual_share_isolated,
            scan_share_isolated,
            scatter_share_isolated,
            count_share,
            scan_share,
            scatter_share,
            count_instructions: c_count.instructions,
            scan_instructions: c_scan.instructions,
            scatter_instructions: c_scatter.instructions,
            total_instructions: c_total.instructions,
            residual_share_instructions,
            count_share_instructions,
            scan_share_instructions,
            scatter_share_instructions,
            count_cycles: c_count.cycles,
            scan_cycles: c_scan.cycles,
            scatter_cycles: c_scatter.cycles,
            total_cycles: c_total.cycles,
            residual_share_cycles,
            scan_share_cycles,
            scatter_share_cycles,
            rank_total_ms: c_rank.ms(),
            rank_total_instructions: c_rank.instructions,
            rank_total_cycles: c_rank.cycles,
            speedup,
            speedup_instructions,
            speedup_cycles,
            speedup_ceiling_scan_free,
            rank_queries: active_cells,
            rank_queries_per_cell: f64::from(active_cells) / cells as f64,
            mean_words_folded_per_query,
            max_words_folded_per_query,
            sequential_ms: c_seq.ms(),
            sequential_instructions: c_seq.instructions,
            three_phase_over_sequential: total_ms / c_seq.ms(),
            three_phase_over_sequential_instructions: c_total.instructions / c_seq.instructions,
            rank_word_ms: c_rank_word.ms(),
            rank_word_instructions: c_rank_word.instructions,
            speedup_rank_word_unregistered: total_ms / c_rank_word.ms(),
            speedup_rank_word_instructions_unregistered: c_total.instructions
                / c_rank_word.instructions,
            queries_rank_word,
            rank_row_ms: c_rank_row.ms(),
            rank_row_instructions: c_rank_row.instructions,
            speedup_rank_row_unregistered: total_ms / c_rank_row.ms(),
            speedup_rank_row_instructions_unregistered: c_total.instructions
                / c_rank_row.instructions,
            queries_rank_row,
            ms_extract_mc: c_extract.ms(),
            cycles_extract_mc: c_extract.cycles,
            instructions_extract_mc: c_extract.instructions,
            compaction_share_of_extraction,
            compaction_share_of_extraction_instructions,
            extraction_ceiling_scan_free,
            target_feature_popcnt: cfg!(target_feature = "popcnt"),
            count_ones_calls_count_phase: words,
            count_ones_calls_scan_phase: 0,
            count_ones_calls_scatter_phase: 0,
            count_ones_calls_rank_arm,
            count_instructions_per_word,
            scan_share_instructions_popcnt_floor,
            compacted_len,
            output_identical,
            output_identical_all_arms,
            bitmap_matches_scalar,
            rank_equal_all_bits,
            cells_checked: bits.cell_bits,
            short_directory_mismatches,
            mutant_output_mismatches,
            directory_total_equals_popcount,
            pad_bits_set,
            pad_bits_possible,
            ghz: c_extract.cycles / c_extract.nanos,
            inner_three,
            inner_rank,
            inner_extract,
            reps: REPS,
            c1_holds,
            c1_holds_instructions: speedup_instructions >= SPEEDUP_BAR,
            c2_holds,
            c2_holds_instructions: scan_share_instructions >= SCAN_SHARE_BAR,
            scatter_falsifier_fires: scatter_share > SCATTER_FALSIFIER,
            c3_holds,
        }
    }

    /// Eight reference fields × {65³, 129³}, `f32`.
    ///
    /// No `scalar` column is registered and none is added: every arm is integer
    /// work over a bitmap and does not change with the field's precision.
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

        let c1 = rows.iter().filter(|r| r.c1_holds).count();
        let c1_i = rows.iter().filter(|r| r.c1_holds_instructions).count();
        let c2 = rows.iter().filter(|r| r.c2_holds).count();
        let c2_i = rows.iter().filter(|r| r.c2_holds_instructions).count();
        let falsifier = rows.iter().filter(|r| r.scatter_falsifier_fires).count();
        let best_speedup = rows.iter().map(|r| r.speedup).fold(0.0f64, f64::max);
        let best_scan_share = rows.iter().map(|r| r.scan_share).fold(0.0f64, f64::max);
        let worst_residual = rows
            .iter()
            .map(|r| r.residual_share.abs())
            .fold(0.0f64, f64::max);
        let agree = rows
            .iter()
            .filter(|r| {
                r.c1_holds == r.c1_holds_instructions && r.c2_holds == r.c2_holds_instructions
            })
            .count();
        let best_rank_word = rows
            .iter()
            .map(|r| r.speedup_rank_word_unregistered)
            .fold(0.0f64, f64::max);
        let best_rank_row = rows
            .iter()
            .map(|r| r.speedup_rank_row_unregistered)
            .fold(0.0f64, f64::max);
        let worst_ceiling = rows
            .iter()
            .map(|r| r.extraction_ceiling_scan_free)
            .fold(0.0f64, f64::max);

        println!(
            "P-112: C1 (speedup >= {SPEEDUP_BAR}) holds on {c1} of {} rows on milliseconds and \
             {c1_i} on instructions; the largest speedup measured is {best_speedup:.4}",
            rows.len()
        );
        println!(
            "P-112: C2 (scan_share >= {SCAN_SHARE_BAR}) holds on {c2} of {} rows on \
             milliseconds and {c2_i} on instructions; the largest scan share measured is \
             {best_scan_share:.4}, whose scan-free ceiling is {:.4}x",
            rows.len(),
            1.0 / (1.0 - best_scan_share)
        );
        let worst_residual_ms = rows
            .iter()
            .map(|r| r.residual_share_ms)
            .fold(0.0f64, f64::max);
        let worst_spread = rows
            .iter()
            .map(|r| r.total_pre_post_spread_share)
            .fold(0.0f64, f64::max);
        let least_mutant = rows
            .iter()
            .map(|r| r.residual_control_mutant_share)
            .fold(f64::INFINITY, f64::min);
        println!(
            "P-112: the 70% scatter falsifier fires on {falsifier} of {} rows; the two verdict \
             forms agree on {agree} of {} rows",
            rows.len(),
            rows.len()
        );
        println!(
            "P-112: the worst residual share is {worst_residual:.6} in RETIRED INSTRUCTIONS \
             against a bar of {RESIDUAL_BAR}, and the same control with the scatter phase \
             dropped reads at least {least_mutant:.4}, so the zero could have been non-zero. \
             The millisecond form of the same residual reaches {worst_residual_ms:.4} and the \
             two identical `total` windows inside one repetition differ by up to \
             {worst_spread:.4} — which is why it is reported and gated by nothing (M-280)"
        );
        println!(
            "P-112: unregistered arms, reported and scored by nothing — one query per word tops \
             out at {best_rank_word:.4}x and one query per row at {best_rank_row:.4}x"
        );
        println!(
            "P-112: the whole three-phase compaction is {:.4}%-{:.4}% of a marching_cubes \
             extraction, so the best extraction-level speedup a free scan could buy anywhere in \
             this fixture is {worst_ceiling:.5}x",
            rows.iter()
                .map(|r| r.compaction_share_of_extraction)
                .fold(f64::INFINITY, f64::min)
                * 100.0,
            rows.iter()
                .map(|r| r.compaction_share_of_extraction)
                .fold(0.0f64, f64::max)
                * 100.0
        );
        println!(
            "P-112: target_feature_popcnt is {}, so phase 1 makes {} SWAR popcounts and phases \
             2 and 3 make none; both verdicts are contingent on that and the projection column \
             says by how much",
            cfg!(target_feature = "popcnt"),
            rows[0].count_ones_calls_count_phase
        );

        for row in &rows {
            run.record(&[
                ("field", row.field.to_string()),
                ("resolution", row.resolution.to_string()),
                ("cells", row.cells.to_string()),
                ("active_cells", row.active_cells.to_string()),
                ("active_fraction", format!("{:.6}", row.active_fraction)),
                ("bitmap_words", row.bitmap_words.to_string()),
                ("cell_words_per_row", row.cell_words_per_row.to_string()),
                ("sample_bit_row", row.sample_bit_row.to_string()),
                ("bitmap_bytes", row.bitmap_bytes.to_string()),
                ("directory_bytes", row.directory_bytes.to_string()),
                ("overhead_fraction", format!("{:.6}", row.overhead_fraction)),
                ("count_ms", format!("{:.6}", row.count_ms)),
                ("scan_ms", format!("{:.6}", row.scan_ms)),
                ("scatter_ms", format!("{:.6}", row.scatter_ms)),
                ("total_ms", format!("{:.6}", row.total_ms)),
                ("phase_sum_ms", format!("{:.6}", row.phase_sum_ms)),
                ("residual_ms", format!("{:.6}", row.residual_ms)),
                ("residual_share", format!("{:.6}", row.residual_share)),
                (
                    "residual_signed_share",
                    format!("{:.6}", row.residual_signed_share),
                ),
                (
                    "residual_control_mutant_share",
                    format!("{:.6}", row.residual_control_mutant_share),
                ),
                ("residual_share_ms", format!("{:.6}", row.residual_share_ms)),
                (
                    "residual_share_ms_within_rep",
                    format!("{:.6}", row.residual_share_ms_within_rep),
                ),
                (
                    "residual_share_ms_worst_rep",
                    format!("{:.6}", row.residual_share_ms_worst_rep),
                ),
                (
                    "total_pre_post_spread_share",
                    format!("{:.6}", row.total_pre_post_spread_share),
                ),
                ("scan_ms_isolated", format!("{:.6}", row.scan_ms_isolated)),
                (
                    "scatter_ms_isolated",
                    format!("{:.6}", row.scatter_ms_isolated),
                ),
                (
                    "residual_share_isolated",
                    format!("{:.6}", row.residual_share_isolated),
                ),
                (
                    "scan_share_isolated",
                    format!("{:.6}", row.scan_share_isolated),
                ),
                (
                    "scatter_share_isolated",
                    format!("{:.6}", row.scatter_share_isolated),
                ),
                ("count_share", format!("{:.6}", row.count_share)),
                ("scan_share", format!("{:.6}", row.scan_share)),
                ("scatter_share", format!("{:.6}", row.scatter_share)),
                (
                    "count_instructions",
                    format!("{:.1}", row.count_instructions),
                ),
                ("scan_instructions", format!("{:.1}", row.scan_instructions)),
                (
                    "scatter_instructions",
                    format!("{:.1}", row.scatter_instructions),
                ),
                (
                    "total_instructions",
                    format!("{:.1}", row.total_instructions),
                ),
                (
                    "residual_share_instructions",
                    format!("{:.6}", row.residual_share_instructions),
                ),
                (
                    "count_share_instructions",
                    format!("{:.6}", row.count_share_instructions),
                ),
                (
                    "scan_share_instructions",
                    format!("{:.6}", row.scan_share_instructions),
                ),
                (
                    "scatter_share_instructions",
                    format!("{:.6}", row.scatter_share_instructions),
                ),
                ("count_cycles", format!("{:.1}", row.count_cycles)),
                ("scan_cycles", format!("{:.1}", row.scan_cycles)),
                ("scatter_cycles", format!("{:.1}", row.scatter_cycles)),
                ("total_cycles", format!("{:.1}", row.total_cycles)),
                (
                    "residual_share_cycles",
                    format!("{:.6}", row.residual_share_cycles),
                ),
                ("scan_share_cycles", format!("{:.6}", row.scan_share_cycles)),
                (
                    "scatter_share_cycles",
                    format!("{:.6}", row.scatter_share_cycles),
                ),
                ("rank_total_ms", format!("{:.6}", row.rank_total_ms)),
                (
                    "rank_total_instructions",
                    format!("{:.1}", row.rank_total_instructions),
                ),
                ("rank_total_cycles", format!("{:.1}", row.rank_total_cycles)),
                ("speedup", format!("{:.6}", row.speedup)),
                (
                    "speedup_instructions",
                    format!("{:.6}", row.speedup_instructions),
                ),
                ("speedup_cycles", format!("{:.6}", row.speedup_cycles)),
                (
                    "speedup_ceiling_scan_free",
                    format!("{:.6}", row.speedup_ceiling_scan_free),
                ),
                ("rank_queries", row.rank_queries.to_string()),
                (
                    "rank_queries_per_cell",
                    format!("{:.6}", row.rank_queries_per_cell),
                ),
                (
                    "mean_words_folded_per_query",
                    format!("{:.4}", row.mean_words_folded_per_query),
                ),
                (
                    "max_words_folded_per_query",
                    row.max_words_folded_per_query.to_string(),
                ),
                ("sequential_ms", format!("{:.6}", row.sequential_ms)),
                (
                    "sequential_instructions",
                    format!("{:.1}", row.sequential_instructions),
                ),
                (
                    "three_phase_over_sequential",
                    format!("{:.6}", row.three_phase_over_sequential),
                ),
                (
                    "three_phase_over_sequential_instructions",
                    format!("{:.6}", row.three_phase_over_sequential_instructions),
                ),
                ("rank_word_ms", format!("{:.6}", row.rank_word_ms)),
                (
                    "rank_word_instructions",
                    format!("{:.1}", row.rank_word_instructions),
                ),
                (
                    "speedup_rank_word_unregistered",
                    format!("{:.6}", row.speedup_rank_word_unregistered),
                ),
                (
                    "speedup_rank_word_instructions_unregistered",
                    format!("{:.6}", row.speedup_rank_word_instructions_unregistered),
                ),
                ("queries_rank_word", row.queries_rank_word.to_string()),
                ("rank_row_ms", format!("{:.6}", row.rank_row_ms)),
                (
                    "rank_row_instructions",
                    format!("{:.1}", row.rank_row_instructions),
                ),
                (
                    "speedup_rank_row_unregistered",
                    format!("{:.6}", row.speedup_rank_row_unregistered),
                ),
                (
                    "speedup_rank_row_instructions_unregistered",
                    format!("{:.6}", row.speedup_rank_row_instructions_unregistered),
                ),
                ("queries_rank_row", row.queries_rank_row.to_string()),
                ("ms_extract_mc", format!("{:.6}", row.ms_extract_mc)),
                ("cycles_extract_mc", format!("{:.1}", row.cycles_extract_mc)),
                (
                    "instructions_extract_mc",
                    format!("{:.1}", row.instructions_extract_mc),
                ),
                (
                    "compaction_share_of_extraction",
                    format!("{:.8}", row.compaction_share_of_extraction),
                ),
                (
                    "compaction_share_of_extraction_instructions",
                    format!("{:.8}", row.compaction_share_of_extraction_instructions),
                ),
                (
                    "extraction_ceiling_scan_free",
                    format!("{:.8}", row.extraction_ceiling_scan_free),
                ),
                (
                    "target_feature_popcnt",
                    row.target_feature_popcnt.to_string(),
                ),
                (
                    "count_ones_calls_count_phase",
                    row.count_ones_calls_count_phase.to_string(),
                ),
                (
                    "count_ones_calls_scan_phase",
                    row.count_ones_calls_scan_phase.to_string(),
                ),
                (
                    "count_ones_calls_scatter_phase",
                    row.count_ones_calls_scatter_phase.to_string(),
                ),
                (
                    "count_ones_calls_rank_arm",
                    format!("{:.1}", row.count_ones_calls_rank_arm),
                ),
                (
                    "count_instructions_per_word",
                    format!("{:.4}", row.count_instructions_per_word),
                ),
                (
                    "scan_share_instructions_popcnt_floor",
                    format!("{:.6}", row.scan_share_instructions_popcnt_floor),
                ),
                ("compacted_len", row.compacted_len.to_string()),
                ("output_identical", row.output_identical.to_string()),
                (
                    "output_identical_all_arms",
                    row.output_identical_all_arms.to_string(),
                ),
                (
                    "bitmap_matches_scalar",
                    row.bitmap_matches_scalar.to_string(),
                ),
                ("rank_equal_all_bits", row.rank_equal_all_bits.to_string()),
                ("cells_checked", row.cells_checked.to_string()),
                (
                    "short_directory_mismatches",
                    row.short_directory_mismatches.to_string(),
                ),
                (
                    "mutant_output_mismatches",
                    row.mutant_output_mismatches.to_string(),
                ),
                (
                    "directory_total_equals_popcount",
                    row.directory_total_equals_popcount.to_string(),
                ),
                ("pad_bits_set", row.pad_bits_set.to_string()),
                ("pad_bits_possible", row.pad_bits_possible.to_string()),
                ("ghz", format!("{:.4}", row.ghz)),
                ("inner_three", row.inner_three.to_string()),
                ("inner_rank", row.inner_rank.to_string()),
                ("inner_extract", row.inner_extract.to_string()),
                ("reps", row.reps.to_string()),
                ("c1_holds", row.c1_holds.to_string()),
                (
                    "c1_holds_instructions",
                    row.c1_holds_instructions.to_string(),
                ),
                ("c2_holds", row.c2_holds.to_string()),
                (
                    "c2_holds_instructions",
                    row.c2_holds_instructions.to_string(),
                ),
                (
                    "scatter_falsifier_fires",
                    row.scatter_falsifier_fires.to_string(),
                ),
                ("c3_holds", row.c3_holds.to_string()),
            ]);
        }
    }
}

fn main() {
    if !std::env::args().any(|arg| arg == "--bench") {
        return;
    }

    let prereg = isomesh::experiment!("P-112");

    #[cfg(target_os = "linux")]
    common::experiment::run(prereg, experiment::run);

    // `experiment_p12`'s precedent. Every registered ratio here has an
    // instruction-form twin, because a millisecond ratio on a governed CPU does
    // not reproduce (`M-281`) — and an instruction count needs
    // `perf_event_open`. A recorded zero would be a fabricated measurement.
    #[cfg(not(target_os = "linux"))]
    {
        eprintln!(
            "{} scores a stage decomposition and needs the instruction-count form of every \
             ratio to do it honestly, and this platform has no `perf_event_open`.",
            prereg.id
        );
        std::process::exit(1);
    }
}
