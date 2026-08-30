//! **P-108 — broadword select to walk the set bits, with no `PEXT` and no table.**
//!
//! Ticket: R-108. Pre-registered before this harness existed.
//!
//! ```bash
//! cargo bench --bench experiment_p108
//! ```
//!
//! Writes `docs/experiments/p-108.csv`. **Linux only**, `experiment_p12`'s
//! precedent: every verdict here is a cycle count or an instruction count, and
//! `M-281` forbids a nanosecond on a governed CPU from carrying one. Off Linux
//! this refuses and exits 1 rather than recording a zero.
//!
//! # C1 is a registered expected null, and it came out null
//!
//! The 2026-08-23 dossier already argued that the cost of `dual.rs:489-497` is
//! dominated by the words neither arm touches, so *"broadword select beats the
//! walk"* was registered in the expectation that it would not. The deliverable
//! is therefore **the number**, not the verdict — and the numbers beside it that
//! say *why*, because a null with a mechanism is worth more than a null.
//!
//! What the sixteen rows say. **The instruction-derived figures below are the
//! reproducible ones**; the cycle-derived ones are from the harness-time run on
//! a machine shared with six sibling benches, and the clean-tree serial re-run
//! is the authoritative dataset for them.
//!
//! - **C1 falsified on every row of every run, by 3.2× to 4.8× in
//!   instructions.** `instruction_ratio` (walk over select, per set bit) is
//!   **0.209 to 0.317**, median 0.248, and that band held on all sixteen rows of
//!   five runs under two different estimators. `cycle_ratio` agrees on the
//!   verdict on every row of every run and spanned **0.098 to 0.456** across six
//!   runs of the same binary — which is `R-105`'s point, and is why the
//!   instruction form is the one quoted.
//! - **Baseline-free, the price is +48.9 to +61.0 instructions per set bit**
//!   (median +53.0), and that is the number to carry:
//!   `extra_instructions_per_set_bit_select` reproduced that band in every run
//!   under both estimators, because a difference cancels the shared baseline and
//!   a ratio does not. Its spread across sixteen fields, resolutions and
//!   densities is only **1.25×**, so the price of broadword select is a property
//!   of the mechanism rather than of the fixture. The *ratio*'s row-to-row
//!   variation is the baseline's: the walk marginal moves 13.1→28.3
//!   instructions per set bit across rows while the extra holds near 53, so
//!   `instruction_ratio` must not be read as a ranking of fields.
//! - **The clause was Amdahl-dead by three orders of magnitude.**
//!   `walk_share_of_extraction` is **0.00008 to 0.00127**, so `c1_ceiling` —
//!   what deleting the walk *entirely* could give a whole dual-contouring
//!   extraction — is **1.00008× to 1.00127×**. `✗51`'s ceiling was 1.061×, and
//!   it was reached by a clause nobody thought was dead.
//! - **The candidate's rationale does not survive measurement either.** Its
//!   theoretical advantage is a countable trip count in place of a
//!   data-dependent loop exit, and the two arms' branch-miss rates are
//!   **indistinguishable**: 0.0083–0.1314 misses per set bit on the walk arm
//!   against 0.0083–0.1303 on the select arm, row for row, agreeing to three
//!   figures. The incumbent's loop-exit branch is already free; what mispredicts
//!   is the shared outer word loop, which both arms pay identically.
//! - **C3's discipline turned out to be load-bearing.** `ratio` (per set bit)
//!   has median 0.210 while `whole_pass_ratio`, which is identically the
//!   *per-word* ratio, has median 0.621. Quoting per-word would have reported a
//!   1.6× loss where the mechanism costs 4.8×.
//! - **The words neither arm touches are most of the pass, as the dossier said.**
//!   `scan_share_of_walk_pass` is **0.492 to 0.954**, median 0.862.
//!
//! # What was missing
//!
//! **Nobody had ever measured the walk against another walk.** `R-039`/`M-337`
//! measured the *bitmap* against the eight-corner scalar gather — the predicate,
//! not the enumeration — and got 5.5× on the stage. The enumeration inside it,
//! `trailing_zeros` plus `active &= active - 1`, has had no comparand since it
//! was written, so there was no number for what a set bit costs to produce.
//!
//! **The `PEXT` foreclosure had no price.** The crate forbids `unsafe`,
//! `core::simd` and `#[cfg(target_arch)]`, so `_pext_u64` is unreachable, and
//! nothing in the repository said what that costs. It costs a known amount:
//! Pandey, Bender and Johnson measured the `PDEP`+`TZCNT` select at **2–4× the
//! broadword one on Haswell**, so the foreclosure is being paid at a published
//! rate rather than an unknown one. That is on every row as
//! `published_pdep_penalty` rather than in a footnote.
//!
//! **The regime was posited rather than cited.** It does not have to be:
//! `M-337`'s measured active fractions are **1.89% at 64³, 0.93% at 128³ and
//! 0.46% at 256³** on `sphere` (`docs/experiments/p-40.csv`, `active_fraction`
//! `0.018916 / 0.009280 / 0.004630`). A 97%-zero bitmap is the measured regime.
//! This harness re-measures it per field anyway, as `active_fraction`, and
//! reproduces `M-337` on the shared field: **0.018158 at 65³ and 0.009171 at
//! 129³** on `sphere`, within 4% of `M-337`'s numbers one sample away. Fourteen
//! of sixteen rows are at or below 5.7% active; the two that are not are
//! `gyroid` (8.2%) and `noise_cavity` (10.8%) at 65³, and both are named rather
//! than dropped.
//!
//! **And one thing nobody could have registered, because it is a property of
//! this crate's build rather than of either algorithm.** `count_ones` is *not an
//! instruction here*. The workspace sets no `target-cpu` and there is no
//! `.cargo/config.toml`, so `isomesh` compiles for baseline `x86-64`, where
//! `POPCNT` (SSE4.2) is absent — `u64::count_ones` lowers to a **19-instruction**
//! SWAR sequence, while `u64::trailing_zeros` lowers to `rep bsfq`, which *is*
//! the `TZCNT` encoding and so gets the hardware instruction on any CPU with
//! BMI1 while degrading to `BSF` on any without. So the incumbent gets hardware
//! help for free and the challenger's popcount does not exist. Column:
//! `target_feature_popcnt`, read from `cfg!(target_feature = "popcnt")` rather
//! than asserted.
//!
//! It is **priced** rather than only named, twice over and the two agree. From
//! disassembly, a standalone `u64::count_ones` on this target is 19
//! instructions: three `movabsq` constant materialisations and sixteen
//! shift/mask/add/multiply. From measurement, a fifth arm is the scan loop plus
//! one `count_ones` per non-zero word and nothing else, so
//! `instructions_per_count_ones` and `cycles_per_count_ones` are its cost **in
//! this binary in this run**, which is the only comparison `M-281` permits. On
//! the **nine** rows where that marginal is at least
//! [`COUNT_ONES_PRICE_FLOOR`] of the pass it reads **22.9 to 30.6 instructions
//! per call**, median 25.3, at 4.0–11.2 cycles — consistent with 19 plus the
//! loop's accumulate and constant reloads.
//!
//! The other seven rows are flagged, not averaged in.
//! `count_ones_price_share` is the fraction of the pass the price rests on and
//! `count_ones_price_readable` is the gate: on `thin_plate` only 64 of 4,096
//! words are non-zero, so the marginal is 1–2% of the pass and the reading —
//! **150** instructions per call — is the two loops' compilation difference
//! rather than a price. A number that small a difference cannot support is
//! marked rather than quoted.
//!
//! And the calls are counted per arm, because a price with no call count settles
//! nothing:
//!
//! | arm | `count_ones` per set bit | column |
//! |---|---|---|
//! | `walk` (incumbent) | **0** | `count_ones_per_set_bit_walk` |
//! | `select_broadword` | **0** | `count_ones_per_set_bit_select_broadword` |
//! | `select_literal` | 1 per rank + 1 per non-zero word = **1.031–1.347** | `count_ones_per_set_bit_select_literal` |
//!
//! **So C1's verdict cannot flip under a hardware popcount, and the arithmetic
//! is on the row.** The verdict is taken on the arm that makes **zero**
//! `count_ones` calls (`count_ones_per_set_bit_select` is `0` on all sixteen
//! rows), so `-C target-cpu=native` cannot make it cheaper by a single
//! instruction. The literal arm makes 1.031–1.347 calls per set bit and would
//! save about `1.15 × (25.4 − 1) ≈ 28` instructions per set bit, landing at
//! roughly 48–72 — still above the broadword arm's measured 62–89 at the bottom
//! of the band and nowhere near the incumbent's 13.1–28.3. The incumbent would
//! itself get cheaper, since `active &= active - 1` becomes one `BLSR` instead
//! of a `lea`/`sub` pair. So a hardware popcount narrows the gap from about 4×
//! to about 3× and cannot close it. That arithmetic is *stated* rather than
//! measured on purpose: `M-281` forbids carrying a verdict across two binaries.
//!
//! That is a fact about the build, not about Vigna, and it would be a handicap
//! to score the candidate on it. So the select arm is measured **twice**:
//!
//! - `literal` — Vigna's published routine verbatim, `place` from
//!   `count_ones`, one `count_ones` per word for the trip count. This is
//!   *"multiply, shift, mask and `count_ones`"* as the ticket words it.
//! - `broadword` — the same routine with **no popcount at all**. Phase 1
//!   already computes the byte-cumulative sums, whose top byte *is* the word's
//!   popcount, so the trip count is `byte_sums(x) >> 56`; and `place` is a
//!   multiply-and-shift over the eight lane bits instead of a popcount over
//!   them. Phase 1 is hoisted out of the rank loop and computed once per
//!   non-zero word.
//!
//! The registered `ns_per_set_bit_select` and `ratio` score **whichever of the
//! two is cheaper**, named in `select_variant`; both are on the row in full. A
//! registered expected null must be given the candidate's best form or the null
//! is the harness's, not the mechanism's.
//!
//! The choice is made on **instructions**, not cycles, and it was made on cycles
//! first: a run taken while six sibling benches were on the machine flipped
//! three of sixteen rows to `literal` on a cycle comparison, which then dragged
//! those rows' *instruction* columns to the literal arm's — a column selected by
//! noise, which is worse than either column. On the deterministic quantity
//! `broadword` wins **every row of every run**, by 12–21 instructions per set
//! bit, so `select_variant` is stable and the columns are comparable between
//! runs. The build's missing `POPCNT` is therefore not what falsified C1; the
//! popcount-free transcription is what the verdict is taken on.
//!
//! # The mirror
//!
//! `crates/isomesh/src/**` is read-only for Phase 25 and `DualMesher` is
//! `pub(crate)`, so the bitmap is rebuilt bench-local from `dual.rs`: the
//! sample-plane packing at `359-381`, `inside_word` at `385`,
//! `inside_word_shifted` at `395`, the fused `any & !all` at `424`, `cell_mask`
//! at `445`, the `cell_words` row length at `484`, and `cell_index` at `767`.
//! Only `is_inside` is public and it is *called*, not re-derived, so the bit is
//! the crate's bit.
//!
//! The one-word asymmetry Phase 25 corrected is reproduced rather than smoothed
//! over: the bitmap is one bit per **sample**, so `bit_row = n.div_ceil(64)` is
//! **3** at 129³ while `cell_words = (n − 1).div_ceil(64)` is **2**. Both are
//! columns (`bitmap_words`, `words`).
//!
//! The bitmap is built once per row and both arms walk that same buffer. This is
//! a walk benchmark; rebuilding the bitmap per pass would measure the packing.
//!
//! # Sibling windows, never nested
//!
//! `R-121` paid for this: Zen 3 has six general-purpose counters and
//! `common::counters::Probe` opens exactly six, so a probe window inside another
//! probe window multiplexes and `Counts::worst_ratio` refuses. Everything here
//! is a **sibling** window and every derived quantity is a difference of two of
//! them.
//!
//! Four arms are measured per row, in the same order every repetition:
//!
//! - `scan` — the whole triple loop, computing `active_word & cell_mask` for
//!   every cell word and then doing **nothing but the zero test the three walk
//!   arms all open with**: `if active == 0 { continue }`, plus one increment on
//!   the ~3–20% of words that are not zero. Its return value is
//!   `words_nonzero` and is asserted against an independent count, so a
//!   baseline that scanned a different word set could not go unnoticed. This is
//!   the cost of the words neither arm touches.
//! - `walk` — `scan`'s loop plus `dual.rs:489-497`'s enumeration.
//! - `select_literal`, `select_broadword` — `scan`'s loop plus the two select
//!   enumerations.
//! - `popcount` — `scan`'s loop plus one `count_ones` per non-zero word and
//!   nothing else, so `(popcount − scan) / words_nonzero` prices
//!   `u64::count_ones` in this build. Its return value is `set_bits` and is
//!   asserted, for the same reason the baseline's is.
//!
//! The registered per-set-bit columns are **marginal**:
//! `(arm − scan) / set_bits`. That subtraction is the whole content of C3. Both
//! select arms take `dual.rs:495`'s zero test before doing any rank work, so all
//! five arms share one per-word control-flow shape and the difference is
//! genuinely per set bit.
//!
//! **What a difference of two windows can and cannot say.** `P-121` established
//! that a stage measured in isolation is not the stage: the outer loop of an arm
//! with an empty body may compile differently from the same loop with a body.
//! That residual is why the baseline-free columns exist. Any error in the
//! baseline shifts *both* marginals by the same absolute amount, so it moves the
//! ratio **toward 1** — which is the direction that would *rescue* C1, not the
//! direction that condemns it. The difference
//! `extra_instructions_per_set_bit_select` cancels the baseline exactly, and its
//! band across sixteen rows (48.9–61.0) is tighter than the ratio's, which is
//! what says the mechanism's price is a property of the mechanism.
//!
//! The emission — one `Vec<u32>` push of `cell_index` per set bit, capacity
//! reserved outside every window — is inside all three walk arms and outside
//! `scan`, so it is charged to the walk. It is about two instructions against
//! the incumbent's four and the challenger's fifty, so including it *narrows*
//! the gap; the choice is named here because it is the conservative one against
//! C1's expected null rather than the flattering one.
//!
//! # SHARE
//!
//! Each clause's reachable share is a column.
//!
//! - **C1's share is `walk_share_of_extraction`** — the marginal cycles of the
//!   set-bit walk over the shipped `DualContouring::extract` on the same field
//!   and grid, measured in the same run as a fifth sibling window. `c1_ceiling`
//!   is `1/(1 − share)`: the most that deleting the walk *entirely* could give a
//!   whole extraction. `✗51`'s rule applied to this row rather than remembered
//!   about it. `extraction_cycles_per_cell` is the denominator, on the row.
//!   Measured: **0.00008 to 0.00127**, so the ceiling is **1.00008× to
//!   1.00127×**. The walk is between one part in eight hundred and one part in
//!   twelve thousand of an extraction, and the extremes are set by the *field*,
//!   not by the walk: `fbm_terrain` costs about 1,140 cycles per cell to
//!   evaluate against `thin_plate`'s 40.
//! - **C2 has no share.** It is an equality: `order_identical` is true only when
//!   all three arms emit element-for-element identical `Vec<u32>`s — same order
//!   *and* same indices. `dual.rs:491-494` says in the source why that matters:
//!   the ascending-`x` order is what keeps vertex creation order and therefore
//!   every index. It is asserted as well as recorded, and it held on all
//!   sixteen rows.
//! - **C3's share is `scan_share_of_walk_pass`** — the fraction of the walk
//!   arm's whole pass that is the shared word scan, which is exactly the
//!   quantity C3 exists to keep out of the numerator. Measured: **0.492 to
//!   0.954**, median 0.862 — so on the typical row six sevenths of the pass is
//!   words neither arm touches, which is the dossier's argument turned into a
//!   number. Beside `ratio` is `whole_pass_ratio`, which is identically the
//!   **per-word** ratio, because a per-word denominator cancels: quoting
//!   per-word means quoting a number diluted by the words neither arm touches,
//!   and the dilution is a factor of **3.0 at the median** (0.210 against
//!   0.621). Having both on the row makes it visible instead of arguable.
//!
//! # Which form carries the verdict
//!
//! **`c1_holds` reads the instruction form, and that decision was forced.** It
//! read the cycle form first. `R-105` had already watched an identical binary's
//! cycle ratio band drift from 0.984 to 1.035 across three runs while its
//! instruction counts held to four figures; this harness, run six times on a
//! machine shared with six sibling benches, saw worse. `cycle_ratio` spanned
//! **0.098 to 0.456** across those runs, and on one of them the cycle
//! *marginal* — a difference of two windows — came out **negative on two of
//! sixteen rows** and flipped `box_exact 65³` to a 1.376× "win". That is a
//! measurement failure wearing the clause's clothes, and no clause may be
//! decided by a quantity that can do it. Over the same six runs
//! `instruction_ratio` was **0.209–0.317** in five and 0.212–0.346 in the first,
//! before the scan baseline was given the walk arms' control-flow shape; and
//! consecutive runs under one estimator agree row for row to four figures.
//!
//! The precision claim has to be stated exactly, because the marginal is a
//! difference of two large window totals and the baseline is 49–95% of the pass.
//! With the estimator *fixed*, per-row instruction ratios reproduce to four
//! figures between runs. **Changing the estimator from median to minimum left
//! ten of sixteen rows bit-identical and moved the other six by up to 10%** —
//! per-window instruction totals vary by a fraction of a percent, and a 15%
//! marginal amplifies that. So the *band* and the *baseline-free difference* are
//! the reproducible quantities and the per-row ratio is not, which is why
//! `extra_instructions_per_set_bit_select` is the number this row asks a reader
//! to carry: 48.9–61.0 in every run under both estimators.
//!
//! The instruction form is a *sound* proxy here as well as a stable one, and the
//! argument runs in the conservative direction. Both arms are pure
//! ALU-and-L1 work over one buffer, with branch-miss rates that are
//! indistinguishable row for row, so cycles differ from instructions only
//! through IPC — and broadword select is a long **dependent** chain
//! (`byte_sums → geq → place → byte_rank → spread → bit_sums → leq`) where the
//! incumbent is three independent instructions. If the two IPCs differ, the
//! challenger's is the worse one. So a 3.2×-to-4.8× instruction gap is a *lower*
//! bound on the cycle gap, which is what the rows show: cycle median 0.211
//! against instruction median 0.248.
//!
//! `c1_holds_cycles` records the cycle verdict beside it and
//! `cycle_marginal_valid` says whether that row's cycle columns are readable at
//! all — both marginals positive. Two independently contaminated readings can
//! give a negative difference where the true difference is positive, so the
//! estimator over repetitions is the **minimum** rather than the median:
//! contention only adds cycles, so the minimum converges on the cost while a
//! median is robust to one disturbed repetition and not to persistent
//! contention. Switching to it took `c1_holds_cycles` back to `false` on all
//! sixteen rows; the run that shipped this file has
//! `cycle_marginal_valid` true on all sixteen, and the flag exists because an
//! earlier one did not.
//!
//! `ns_*` columns are provenance and carry no clause; `ghz` is on
//! every row so a later reader can see the clock they were taken at —
//! 4.182–4.190 GHz across the shipping run, stable to 0.2%, which is *why* the
//! ns and cycle forms agree there and is not a reason to trust the ns
//! form (`M-280`, `M-281`).
//!
//! # Vacuity control
//!
//! Registered: `words_nonzero` must be non-zero and `set_bits` must **exceed**
//! `words_nonzero`, or the per-set-bit denominator is one bit per word and C3
//! measures nothing. Both are asserted per row, and
//! `set_bits_per_nonzero_word` is the margin. Measured: `words_nonzero` runs
//! **64 to 16,979** and `set_bits_per_nonzero_word` **2.886 to 32.000**, median
//! 6.96 — so the tightest row still carries nearly three set bits per non-zero
//! word and the denominator is a set-bit denominator on all sixteen.
//!
//! Not registered, and added because the comparator would otherwise be checked
//! only against the fixture it is measuring: both select routines are verified
//! exhaustively against a linear scan before any timing, over every single-bit
//! word, every two-bit word, a seeded pseudo-random population and the sparse
//! shapes the fixture itself produces — **149,083** `select_ranks_verified`
//! `(word, rank)` pairs, plus `byte_sums(x) >> 56 == x.count_ones()` on every
//! one of those words. A select that answered wrongly would fail C2, but it
//! would fail it *as a list disagreement*, which is a worse diagnosis than a
//! wrong rank.

mod common;

#[cfg(target_os = "linux")]
mod experiment {
    use std::hint::black_box;
    use std::time::Instant;

    use isomesh::dual_contouring::DualContouring;
    use isomesh::fields::ReferenceField;
    use isomesh::marching_cubes::table::is_inside;
    use isomesh::{MeshBuffer, Sdf};

    use crate::common::counters::{MIN_TIME_RATIO, Probe};

    // ─── the fixture ───────────────────────────────────────────────────────

    /// The registered resolutions, in samples per axis.
    const RESOLUTIONS: [u32; 2] = [65, 129];
    /// Measured repetitions per arm; the median of each quantity is reported.
    const REPS: usize = 7;
    /// Repetitions of the whole-extraction share denominator, which is three
    /// orders of magnitude more expensive per pass than a bitmap walk.
    const EXTRACT_REPS: usize = 3;
    /// How long one counter window should last, in nanoseconds.
    const TARGET_BATCH_NS: f64 = 30_000_000.0;
    /// Ceiling on the batch, so a cheap row cannot take minutes.
    const MAX_INNER: usize = 16_384;
    /// Pandey, Bender and Johnson's published price of the `PEXT` foreclosure.
    const PDEP_PENALTY: &str = "2x-4x";
    /// The popcount arm's marginal must be at least this fraction of the pass
    /// for `instructions_per_count_ones` to be a price rather than the two
    /// loops' compilation difference. Fixed here; `count_ones_price_readable`
    /// is the column.
    const COUNT_ONES_PRICE_FLOOR: f64 = 0.05;

    // ─── broadword primitives ──────────────────────────────────────────────
    //
    // Vigna, *Broadword Implementation of Rank/Select Queries*. The paper has no
    // DOI and `paper_download` cannot reach it (`M-415`), so these are
    // transcribed from the published `select_in_word` and then verified
    // exhaustively by [`verify_select`] rather than trusted.
    //
    // Every multiply and every subtraction here is deliberately modular: the
    // products are exactly the ones that would overflow if the lanes were not
    // bounded, and the subtractions are the ones whose borrow is arranged not to
    // cross a byte. `wrapping_*` states that rather than relying on the release
    // profile, so the routines behave identically under `cargo clippy`'s debug
    // build with overflow checks on.

    const ONES_STEP_4: u64 = 0x1111_1111_1111_1111;
    const ONES_STEP_8: u64 = 0x0101_0101_0101_0101;
    const MSBS_STEP_8: u64 = 0x80 * ONES_STEP_8;
    const INCR_STEP_8: u64 = (0x80 << 56)
        | (0x40 << 48)
        | (0x20 << 40)
        | (0x10 << 32)
        | (0x08 << 24)
        | (0x04 << 16)
        | (0x02 << 8)
        | 0x01;

    /// Per byte: `1` at the byte's low bit when `x`'s byte is `≤` `y`'s.
    #[inline]
    fn leq_step_8(x: u64, y: u64) -> u64 {
        ((((y | MSBS_STEP_8).wrapping_sub(x & !MSBS_STEP_8)) ^ x ^ y) & MSBS_STEP_8) >> 7
    }

    /// Per byte: `1` at the byte's low bit when `x`'s byte is zero.
    #[inline]
    fn zcompare_step_8(x: u64) -> u64 {
        ((x | (x | MSBS_STEP_8).wrapping_sub(ONES_STEP_8)) & MSBS_STEP_8) >> 7
    }

    /// Vigna's phase 1: byte `j` of the result is the popcount of bytes `0..=j`.
    ///
    /// So the **top** byte is `x.count_ones()`, which is why the broadword
    /// variant needs no `POPCNT` for its trip count either.
    #[inline]
    fn byte_sums(x: u64) -> u64 {
        let mut s = x.wrapping_sub((x & (0xa * ONES_STEP_4)) >> 1);
        s = (s & (3 * ONES_STEP_4)) + ((s >> 2) & (3 * ONES_STEP_4));
        s = (s + (s >> 4)) & (0x0f * ONES_STEP_8);
        s.wrapping_mul(ONES_STEP_8)
    }

    /// Phase 2's lane comparison: the MSB of byte `j` is set when the cumulative
    /// popcount through byte `j` is `≤ k`, i.e. the target bit is further on.
    #[inline]
    fn geq_lanes(sums: u64, k: u32) -> u64 {
        (u64::from(k).wrapping_mul(ONES_STEP_8) | MSBS_STEP_8).wrapping_sub(sums) & MSBS_STEP_8
    }

    /// Phase 3: the rank-`k` set bit of `x`, given the byte it lives in.
    #[inline]
    fn select_tail(x: u64, sums: u64, k: u32, place: u32) -> u32 {
        let byte_rank = u64::from(k) - (((sums << 8) >> place) & 0xFF);
        let spread = ((x >> place) & 0xFF).wrapping_mul(ONES_STEP_8) & INCR_STEP_8;
        let bit_sums = zcompare_step_8(spread).wrapping_mul(ONES_STEP_8);
        let wanted = byte_rank.wrapping_mul(ONES_STEP_8);
        place + (leq_step_8(bit_sums, wanted).wrapping_mul(ONES_STEP_8) >> 56) as u32
    }

    /// Vigna verbatim: `place` from `count_ones`, phase 1 per call.
    ///
    /// # Panics
    ///
    /// Debug-asserts `k < x.count_ones()`; without that `place` can reach 64 and
    /// the shift in [`select_tail`] is undefined.
    #[inline]
    fn select_in_word(x: u64, k: u32) -> u32 {
        debug_assert!(k < x.count_ones(), "rank past the end of the word");
        let sums = byte_sums(x);
        let place = geq_lanes(sums, k).count_ones() * 8;
        select_tail(x, sums, k, place)
    }

    /// The same routine with **no popcount**: `place` by multiply-and-shift over
    /// the eight lane bits, and phase 1 hoisted to the caller.
    #[inline]
    fn select_from_sums(x: u64, sums: u64, k: u32) -> u32 {
        debug_assert!(k < x.count_ones(), "rank past the end of the word");
        let lanes = geq_lanes(sums, k) >> 7;
        let place = ((lanes.wrapping_mul(ONES_STEP_8) >> 56) as u32) << 3;
        select_tail(x, sums, k, place)
    }

    /// `splitmix64`, so the verification population is reproducible without a
    /// dependency.
    struct Rng(u64);

    impl Rng {
        fn next(&mut self) -> u64 {
            self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
            let mut z = self.0;
            z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
            z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
            z ^ (z >> 31)
        }
    }

    /// Check both select forms against a linear scan, exhaustively over the
    /// structured words and thoroughly over the random ones.
    ///
    /// Returns the number of `(word, rank)` pairs checked.
    ///
    /// # Panics
    ///
    /// On any disagreement, or if the broadword popcount disagrees with
    /// `count_ones`.
    fn verify_select() -> u64 {
        let mut words: Vec<u64> = (0..64).map(|i| 1u64 << i).collect();
        for i in 0..64u32 {
            for j in (i + 1)..64u32 {
                words.push((1u64 << i) | (1u64 << j));
            }
        }
        words.push(u64::MAX);
        words.push(0x8000_0000_0000_0001);
        let mut rng = Rng(0x5EED_1108);
        for _ in 0..4096 {
            words.push(rng.next());
        }
        // The shapes the fixture actually produces: one to six set bits in a
        // word, which is where a select routine's edge cases live.
        for _ in 0..4096 {
            let mut w = 0u64;
            for _ in 0..=(rng.next() % 6) {
                w |= 1u64 << (rng.next() % 64);
            }
            words.push(w);
        }

        let mut pairs = 0u64;
        for x in words {
            if x == 0 {
                continue;
            }
            let sums = byte_sums(x);
            assert_eq!(
                (sums >> 56) as u32,
                x.count_ones(),
                "broadword popcount disagrees with count_ones on {x:#018x}"
            );
            let mut rank = 0u32;
            for bit in 0..64u32 {
                if (x >> bit) & 1 == 0 {
                    continue;
                }
                assert_eq!(
                    select_in_word(x, rank),
                    bit,
                    "literal select({x:#018x}, {rank})"
                );
                assert_eq!(
                    select_from_sums(x, sums, rank),
                    bit,
                    "broadword select({x:#018x}, {rank})"
                );
                rank += 1;
                pairs += 1;
            }
        }
        pairs
    }

    // ─── the bitmap, mirroring `dual.rs` ───────────────────────────────────

    /// The active-cell bitmap, one bit per **sample**, 64 to a `u64`, along X
    /// only — `dual.rs:359-381`.
    struct Bits {
        words: Vec<u64>,
        /// `n.div_ceil(64)`, the **sample** row — `dual.rs:362`.
        bit_row: usize,
        /// Samples per axis.
        n: usize,
        /// Cells per axis, `n − 1`.
        cells: usize,
        /// `cells.div_ceil(64)`, one shorter than [`bit_row`](Self::bit_row) at
        /// both registered resolutions — `dual.rs:484`.
        cell_words: usize,
    }

    impl Bits {
        /// Sample the field on `DualMesher`'s layout and pack the sign bits.
        fn build<S: Sdf<Scalar = f32>>(sdf: &S, n: u32, origin: [f32; 3], cell_size: f32) -> Self {
            let sx = n as usize;
            // `dual.rs:333`: `size[0] | 1`, A-024's odd row stride.
            let row = sx | 1;
            let rows = sx * sx;
            let mut values = vec![0.0f32; row * rows];
            for z in 0..sx {
                for y in 0..sx {
                    let base = row * (y + sx * z);
                    for x in 0..sx {
                        values[base + x] = sdf.sample([
                            origin[0] + cell_size * x as f32,
                            origin[1] + cell_size * y as f32,
                            origin[2] + cell_size * z as f32,
                        ]);
                    }
                }
            }

            let bit_row = sx.div_ceil(64);
            let mut words = vec![0u64; bit_row * rows];
            for r in 0..rows {
                let src = row * r;
                let dst = bit_row * r;
                for w in 0..bit_row {
                    let base = w * 64;
                    let m = (sx - base).min(64);
                    let mut word = 0u64;
                    for k in 0..m {
                        word |= u64::from(is_inside(values[src + base + k])) << k;
                    }
                    words[dst + w] = word;
                }
            }

            let cells = sx - 1;
            Self {
                words,
                bit_row,
                n: sx,
                cells,
                cell_words: cells.div_ceil(64),
            }
        }

        /// `dual.rs:385`.
        #[inline]
        fn word(&self, w: usize, y: usize, z: usize) -> u64 {
            self.words[self.bit_row * (y + self.n * z) + w]
        }

        /// `dual.rs:395`. The high bit comes from the next word or the cell
        /// straddling a word boundary reads its `+x` corner as outside.
        #[inline]
        fn shifted(&self, w: usize, y: usize, z: usize) -> u64 {
            let lo = self.word(w, y, z);
            let hi = if w + 1 < self.bit_row {
                self.word(w + 1, y, z)
            } else {
                0
            };
            (lo >> 1) | (hi << 63)
        }

        /// `dual.rs:424`'s fused `any & !all`, masked by `dual.rs:445`'s
        /// `cell_mask`.
        #[inline]
        fn active(&self, w: usize, y: usize, z: usize) -> u64 {
            let mut any = 0u64;
            let mut all = !0u64;
            for dz in 0..2usize {
                for dy in 0..2usize {
                    let a = self.word(w, y + dy, z + dz);
                    let b = self.shifted(w, y + dy, z + dz);
                    any |= a | b;
                    all &= a & b;
                }
            }
            let remaining = self.cells.saturating_sub(w * 64);
            let mask = if remaining >= 64 {
                !0
            } else {
                (1u64 << remaining) - 1
            };
            (any & !all) & mask
        }

        /// `dual.rs:767`'s `cell_index`.
        #[inline]
        fn index(&self, x: usize, y: usize, z: usize) -> u32 {
            (x + self.cells * (y + self.cells * z)) as u32
        }

        /// Cell words iterated by the walk — `cell_words · cells · cells`.
        fn word_count(&self) -> usize {
            self.cell_words * self.cells * self.cells
        }

        fn cell_count(&self) -> usize {
            self.cells * self.cells * self.cells
        }
    }

    // ─── the four arms ─────────────────────────────────────────────────────

    /// The shared word scan: every cell word computed, **no set bit walked**.
    ///
    /// The per-word body is the zero test all three walk arms open with, and
    /// nothing else — `if active == 0 { continue }`, then one increment on the
    /// ~3% of words that are not zero. So the control flow and the branch
    /// structure are the walk arms', which is what makes the difference a
    /// per-set-bit cost rather than a per-word compilation difference. The
    /// returned count is `words_nonzero`, asserted against [`population`], so
    /// the baseline cannot silently walk a different set of words.
    fn scan_arm(b: &Bits) -> u64 {
        let mut nonzero = 0u64;
        for z in 0..b.cells {
            for y in 0..b.cells {
                for w in 0..b.cell_words {
                    let active = b.active(w, y, z);
                    if active == 0 {
                        continue;
                    }
                    nonzero += 1;
                }
            }
        }
        nonzero
    }

    /// [`scan_arm`] plus **one `count_ones` per non-zero word**, and nothing
    /// else. Subtracting `scan_arm` prices `u64::count_ones` in this build, in
    /// this binary, in this run — which is the only comparison `M-281` permits.
    ///
    /// The workspace sets no `target-cpu`, so `count_ones` is not an
    /// instruction here; `target_feature_popcnt` is the column that says so and
    /// `instructions_per_count_ones` is the price.
    fn popcount_arm(b: &Bits) -> u64 {
        let mut bits = 0u64;
        for z in 0..b.cells {
            for y in 0..b.cells {
                for w in 0..b.cell_words {
                    let active = b.active(w, y, z);
                    if active == 0 {
                        continue;
                    }
                    bits += u64::from(active.count_ones());
                }
            }
        }
        bits
    }

    /// The incumbent — `dual.rs:489-497` exactly.
    fn walk_arm(b: &Bits, out: &mut Vec<u32>) {
        out.clear();
        for z in 0..b.cells {
            for y in 0..b.cells {
                for w in 0..b.cell_words {
                    let mut active = b.active(w, y, z);
                    // `active &= active - 1` clears the lowest set bit, so this
                    // visits the row's active cells in ascending `x`.
                    while active != 0 {
                        let x = w * 64 + active.trailing_zeros() as usize;
                        active &= active - 1;
                        out.push(b.index(x, y, z));
                    }
                }
            }
        }
        black_box(out.as_slice());
    }

    /// Vigna's published routine, once per rank, `count_ones` and all.
    fn select_literal_arm(b: &Bits, out: &mut Vec<u32>) {
        out.clear();
        for z in 0..b.cells {
            for y in 0..b.cells {
                for w in 0..b.cell_words {
                    let active = b.active(w, y, z);
                    if active == 0 {
                        continue;
                    }
                    let n = active.count_ones();
                    for k in 0..n {
                        let x = w * 64 + select_in_word(active, k) as usize;
                        out.push(b.index(x, y, z));
                    }
                }
            }
        }
        black_box(out.as_slice());
    }

    /// The popcount-free form: phase 1 once per non-zero word, its top byte as
    /// the trip count, `place` by multiply-and-shift.
    fn select_broadword_arm(b: &Bits, out: &mut Vec<u32>) {
        out.clear();
        for z in 0..b.cells {
            for y in 0..b.cells {
                for w in 0..b.cell_words {
                    let active = b.active(w, y, z);
                    if active == 0 {
                        continue;
                    }
                    let sums = byte_sums(active);
                    let n = (sums >> 56) as u32;
                    for k in 0..n {
                        let x = w * 64 + select_from_sums(active, sums, k) as usize;
                        out.push(b.index(x, y, z));
                    }
                }
            }
        }
        black_box(out.as_slice());
    }

    /// How many words carry at least one active cell, and how many set bits
    /// they carry — the vacuity control's two numbers.
    fn population(b: &Bits) -> (usize, usize) {
        let mut nonzero = 0;
        let mut set_bits = 0;
        for z in 0..b.cells {
            for y in 0..b.cells {
                for w in 0..b.cell_words {
                    let active = b.active(w, y, z);
                    if active != 0 {
                        nonzero += 1;
                        set_bits += active.count_ones() as usize;
                    }
                }
            }
        }
        (nonzero, set_bits)
    }

    // ─── measurement ───────────────────────────────────────────────────────

    /// One counter window's four quantities, per pass.
    #[derive(Clone, Copy, Default)]
    struct Counted {
        cycles: f64,
        instructions: f64,
        branch_misses: f64,
        ns: f64,
    }

    impl Counted {
        fn minus(self, o: Self) -> Self {
            Self {
                cycles: self.cycles - o.cycles,
                instructions: self.instructions - o.instructions,
                branch_misses: self.branch_misses - o.branch_misses,
                ns: self.ns - o.ns,
            }
        }

        fn per(self, denominator: f64) -> Self {
            let inv = denominator.recip();
            Self {
                cycles: self.cycles * inv,
                instructions: self.instructions * inv,
                branch_misses: self.branch_misses * inv,
                ns: self.ns * inv,
            }
        }
    }

    /// One counter window, undivided. The `perf_event` system calls are all
    /// outside the counted region.
    fn raw_window(probe: &mut Probe, body: impl FnOnce()) -> Counted {
        probe.reset_and_enable();
        let started = Instant::now();
        body();
        let ns = started.elapsed().as_nanos() as f64;
        probe.disable();
        let counts = probe.read();
        assert!(
            counts.worst_ratio() >= MIN_TIME_RATIO,
            "a counter ran only {:.1}% of the time it was enabled, so its value is an \
             extrapolation rather than a measurement",
            counts.worst_ratio() * 100.0
        );
        Counted {
            cycles: counts.cycles.count as f64,
            instructions: counts.instructions.count as f64,
            branch_misses: counts.branch_misses.count as f64,
            ns,
        }
    }

    /// [`raw_window`] over `inner` passes, divided by `inner`.
    fn window(probe: &mut Probe, inner: usize, mut body: impl FnMut()) -> Counted {
        raw_window(probe, || {
            for _ in 0..inner {
                body();
            }
        })
        .per(inner as f64)
    }

    /// The least-contaminated reading of a quantity: its **minimum**.
    ///
    /// Not a median, and the departure from `experiment_p40`'s median-run is
    /// deliberate. Contention only *adds* cycles, cache misses and nanoseconds
    /// to a fixed piece of work, so the minimum over repetitions is the estimate
    /// that converges on the cost; a median is robust to one disturbed
    /// repetition and not to *persistent* contention, which is the actual
    /// condition when seven benches share a machine. It matters here because
    /// every reported cost is a **difference** of two arms: the true scan cost
    /// is strictly below the true walk cost, so estimators that converge give a
    /// positive difference, while two independently inflated medians can give a
    /// negative one. A run of this harness under six sibling benches did exactly
    /// that on two of sixteen rows with medians. Instruction counts are
    /// deterministic, so for them the minimum is also the median and every other
    /// order statistic.
    fn best(v: impl Iterator<Item = f64>) -> f64 {
        v.fold(f64::INFINITY, f64::min)
    }

    /// [`best`] applied to each quantity independently.
    fn best_counted(reps: &[Counted]) -> Counted {
        Counted {
            cycles: best(reps.iter().map(|c| c.cycles)),
            instructions: best(reps.iter().map(|c| c.instructions)),
            branch_misses: best(reps.iter().map(|c| c.branch_misses)),
            ns: best(reps.iter().map(|c| c.ns)),
        }
    }

    // ─── one row ───────────────────────────────────────────────────────────

    /// Everything one `(field, resolution)` measured.
    struct Row {
        field: &'static str,
        resolution: u32,
        cells: usize,
        words: usize,
        bitmap_words: usize,
        words_nonzero: usize,
        set_bits: usize,
        order_identical: bool,
        inner: usize,
        /// Per pass.
        scan: Counted,
        walk: Counted,
        literal: Counted,
        broadword: Counted,
        /// [`scan_arm`] plus one `count_ones` per non-zero word.
        popcount: Counted,
        /// Per whole extraction.
        extraction: Counted,
    }

    impl Row {
        /// The marginal cost of one set bit, over the shared word scan.
        fn per_set_bit(&self, arm: Counted) -> Counted {
            arm.minus(self.scan).per(self.set_bits as f64)
        }

        /// Whichever select variant is cheaper — the candidate's best form,
        /// which is what an expected null owes it. Returns the variant's name,
        /// its marginal cost per set bit, and its whole-pass cost.
        ///
        /// **Chosen on instructions, not cycles**, and that is not a
        /// preference. A run of this harness taken while six sibling benches
        /// were on the machine flipped three of sixteen rows to `literal` on a
        /// cycle comparison, which then moved the *instruction* columns of those
        /// rows too — a noise-driven column choice, which is worse than either
        /// choice. Instruction counts are deterministic (`R-105`, `M-279`), and
        /// on this quantity `broadword` wins every row of every run by 12–21
        /// instructions per set bit, so the choice is stable and the columns are
        /// comparable across runs.
        fn best_select(&self) -> (&'static str, Counted, Counted) {
            let literal = self.per_set_bit(self.literal);
            let broadword = self.per_set_bit(self.broadword);
            if broadword.instructions <= literal.instructions {
                ("broadword", broadword, self.broadword)
            } else {
                ("literal", literal, self.literal)
            }
        }
    }

    /// Measure one `(field, resolution)`.
    ///
    /// # Panics
    ///
    /// On the registered vacuity control, on any arm emitting a different list,
    /// or on a multiplexed counter.
    fn measure<S>(field: &'static str, n: u32, sdf: &S, probe: &mut Probe) -> Row
    where
        S: Sdf<Scalar = f32> + ReferenceField,
    {
        let (shape, origin, cell_size) = crate::common::grid(sdf, n);
        let bits = Bits::build(sdf, n, origin, cell_size);
        let (words_nonzero, set_bits) = population(&bits);

        // ── the registered vacuity control ───────────────────────────────
        assert!(
            words_nonzero > 0,
            "{field} {n}³: no word carries an active cell, so there is nothing to walk"
        );
        assert!(
            set_bits > words_nonzero,
            "{field} {n}³: {set_bits} set bits in {words_nonzero} non-zero words, so the \
             per-set-bit denominator is one bit per word and C3 measures nothing"
        );

        // ── C2: identical order and identical indices ────────────────────
        let mut walk_out = Vec::with_capacity(set_bits);
        let mut literal_out = Vec::with_capacity(set_bits);
        let mut broadword_out = Vec::with_capacity(set_bits);
        walk_arm(&bits, &mut walk_out);
        select_literal_arm(&bits, &mut literal_out);
        select_broadword_arm(&bits, &mut broadword_out);
        assert_eq!(
            walk_out.len(),
            set_bits,
            "{field} {n}³: the walk emitted a different number of cells than the bitmap has \
             set bits"
        );
        let order_identical = walk_out == literal_out && walk_out == broadword_out;
        assert!(
            order_identical,
            "{field} {n}³: a select arm disagreed with `dual.rs:489-497`'s visitation order or \
             its indices, which would change every vertex index"
        );

        // ── the batch size, so one window is about `TARGET_BATCH_NS` ─────
        let started = Instant::now();
        walk_arm(&bits, &mut walk_out);
        let pass_ns = started.elapsed().as_nanos() as f64;
        let inner = ((TARGET_BATCH_NS / pass_ns.max(1.0)).round() as usize).clamp(1, MAX_INNER);

        let mut scan = Vec::with_capacity(REPS);
        let mut walk = Vec::with_capacity(REPS);
        let mut literal = Vec::with_capacity(REPS);
        let mut broadword = Vec::with_capacity(REPS);
        let mut popcount = Vec::with_capacity(REPS);
        // One untimed pass of each arm, so no window pays for a cold i-cache —
        // and the baseline's own count checked against `population`, so a
        // baseline that scanned a different word set could not go unnoticed.
        assert_eq!(
            scan_arm(&bits) as usize,
            words_nonzero,
            "{field} {n}³: the scan baseline and the population count disagree about how many \
             words carry an active cell"
        );
        walk_arm(&bits, &mut walk_out);
        select_literal_arm(&bits, &mut literal_out);
        select_broadword_arm(&bits, &mut broadword_out);
        assert_eq!(
            popcount_arm(&bits) as usize,
            set_bits,
            "{field} {n}³: the count_ones arm and the population count disagree about how many \
             set bits the bitmap carries"
        );
        for _ in 0..REPS {
            scan.push(window(probe, inner, || {
                black_box(scan_arm(&bits));
            }));
            walk.push(window(probe, inner, || walk_arm(&bits, &mut walk_out)));
            literal.push(window(probe, inner, || {
                select_literal_arm(&bits, &mut literal_out);
            }));
            broadword.push(window(probe, inner, || {
                select_broadword_arm(&bits, &mut broadword_out);
            }));
            popcount.push(window(probe, inner, || {
                black_box(popcount_arm(&bits));
            }));
        }

        // ── C1's SHARE denominator: the shipped extraction, same grid ────
        let mut dc = DualContouring::<f32>::new();
        let mut mesh = MeshBuffer::<f32>::new();
        let extract = |dc: &mut DualContouring<f32>, mesh: &mut MeshBuffer<f32>| {
            mesh.reset();
            dc.extract(sdf, &shape, origin, cell_size, mesh)
                .expect("extraction");
        };
        extract(&mut dc, &mut mesh);
        let mut extraction = Vec::with_capacity(EXTRACT_REPS);
        for _ in 0..EXTRACT_REPS {
            extraction.push(raw_window(probe, || extract(&mut dc, &mut mesh)));
        }

        Row {
            field,
            resolution: n,
            cells: bits.cell_count(),
            words: bits.word_count(),
            bitmap_words: bits.words.len(),
            words_nonzero,
            set_bits,
            order_identical,
            inner,
            scan: best_counted(&scan),
            walk: best_counted(&walk),
            literal: best_counted(&literal),
            broadword: best_counted(&broadword),
            popcount: best_counted(&popcount),
            extraction: best_counted(&extraction),
        }
    }

    /// Every `(field, resolution)` the registration names: eight reference
    /// fields × {65³, 129³}, `f32`.
    fn sweep(probe: &mut Probe) -> Vec<Row> {
        let mut rows = Vec::new();
        for n in RESOLUTIONS {
            isomesh::for_each_reference_field!(f32, |name, field| {
                rows.push(measure(name, n, &field, probe));
            });
        }
        rows
    }

    /// `POPCNT` is SSE4.2 and this workspace sets no `target-cpu`, so
    /// `count_ones` is a twelve-instruction SWAR sequence rather than an
    /// instruction. Read from the compiler rather than asserted in prose, and
    /// priced per row by [`popcount_arm`].
    fn target_feature_popcnt() -> bool {
        cfg!(target_feature = "popcnt")
    }

    pub(crate) fn run(run: &mut crate::common::experiment::Run) {
        let verified = verify_select();
        println!(
            "both select forms agree with a linear scan on {verified} (word, rank) pairs; \
             target_feature_popcnt: {}",
            target_feature_popcnt()
        );

        let mut probe = Probe::open();
        let rows = sweep(&mut probe);

        for row in &rows {
            let walk = row.per_set_bit(row.walk);
            let (variant, select, select_pass) = row.best_select();
            let literal = row.per_set_bit(row.literal);
            let broadword = row.per_set_bit(row.broadword);
            let ratio = walk.ns / select.ns;
            let cycle_ratio = walk.cycles / select.cycles;
            let instruction_ratio = walk.instructions / select.instructions;
            let whole_pass_ratio = row.walk.ns / select_pass.ns;
            let scan_share = row.scan.cycles / row.walk.cycles;
            let walk_share = walk.cycles * row.set_bits as f64 / row.extraction.cycles;
            // The price of `u64::count_ones` in this build, measured in this
            // binary and this run: one call per non-zero word, over the same
            // scan baseline. `M-281` forbids reading it off a second build.
            let count_ones_price = row.popcount.minus(row.scan).per(row.words_nonzero as f64);
            // How much of the pass that price rests on. On `thin_plate` only 64
            // of 4,096 words are non-zero, so the popcount arm's marginal is
            // ~1–2% of the pass and the reading (150 instructions per call) is
            // the two loops' compilation difference rather than a price. Flagged
            // rather than dropped or averaged in.
            let count_ones_price_share =
                row.popcount.minus(row.scan).instructions / row.scan.instructions;
            let count_ones_price_readable = count_ones_price_share >= COUNT_ONES_PRICE_FLOOR;
            // How many `count_ones` calls each arm makes per set bit. The
            // incumbent and the broadword select make none; the literal select
            // makes one per rank for `place` plus one per non-zero word for the
            // trip count.
            let literal_count_ones_calls = 1.0 + row.words_nonzero as f64 / row.set_bits as f64;
            // The instruction form carries C1. Not a preference: on a machine
            // shared with six sibling benches the cycle marginal — a difference
            // of two windows — came out NEGATIVE on two of sixteen rows and
            // flipped `box_exact 65³` to a 1.376× "win", which is a
            // measurement failure wearing the clause's clothes. Instruction
            // counts are deterministic (`R-105`, `M-279`) and reproduced to four
            // figures across three runs. The proxy is sound as well as stable:
            // both arms are pure ALU-and-L1 work over one buffer with
            // indistinguishable branch-miss rates, and broadword select is a
            // long *dependent* chain where the incumbent's is three
            // instructions, so if the two IPCs differ it is the challenger's
            // that is worse — which makes a 3.2×-to-4.8× instruction gap the
            // conservative bound on the cycle gap, not an optimistic one.
            let c1_holds = instruction_ratio > 1.0;
            let c1_holds_cycles = cycle_ratio > 1.0;
            // Both marginals must be positive for the cycle columns to mean
            // anything: the true scan cost is strictly below the true walk cost.
            let cycle_marginal_valid = walk.cycles > 0.0 && select.cycles > 0.0;
            let c3_holds = row.words_nonzero > 0 && row.set_bits > row.words_nonzero;
            let ghz = row.walk.cycles / row.walk.ns;

            println!(
                "{:>14} {:>4}³  set_bits {:>7}  words {:>7} ({:>6} nonzero, {:>4.1} bits each)  \
                 walk {:>6.2} → select {:>6.2} ins/bit  ins_ratio {instruction_ratio:5.3} \
                 (cyc {cycle_ratio:6.3}{})  extra {:+6.1} ins/bit  walk share of extraction \
                 {:.5}  [{variant}]",
                row.field,
                row.resolution,
                row.set_bits,
                row.words,
                row.words_nonzero,
                row.set_bits as f64 / row.words_nonzero as f64,
                walk.instructions,
                select.instructions,
                if cycle_marginal_valid { "" } else { " INVALID" },
                select.instructions - walk.instructions,
                walk_share,
            );

            run.record(&[
                ("field", row.field.to_string()),
                ("resolution", row.resolution.to_string()),
                ("scalar", String::from("f32")),
                ("cells", row.cells.to_string()),
                (
                    "active_fraction",
                    format!("{:.6}", row.set_bits as f64 / row.cells as f64),
                ),
                ("set_bits", row.set_bits.to_string()),
                ("words", row.words.to_string()),
                ("words_nonzero", row.words_nonzero.to_string()),
                ("bitmap_words", row.bitmap_words.to_string()),
                (
                    "set_bits_per_nonzero_word",
                    format!("{:.4}", row.set_bits as f64 / row.words_nonzero as f64),
                ),
                ("inner", row.inner.to_string()),
                ("ghz", format!("{ghz:.4}")),
                ("ns_per_set_bit_walk", format!("{:.4}", walk.ns)),
                ("ns_per_set_bit_select", format!("{:.4}", select.ns)),
                ("ratio", format!("{ratio:.4}")),
                ("cycles_per_set_bit_walk", format!("{:.4}", walk.cycles)),
                ("cycles_per_set_bit_select", format!("{:.4}", select.cycles)),
                ("cycle_ratio", format!("{cycle_ratio:.4}")),
                (
                    "instructions_per_set_bit_walk",
                    format!("{:.4}", walk.instructions),
                ),
                (
                    "instructions_per_set_bit_select",
                    format!("{:.4}", select.instructions),
                ),
                ("instruction_ratio", format!("{instruction_ratio:.4}")),
                // Baseline-free: the two marginals share the scan baseline, so
                // it cancels exactly in their difference. A ratio does not have
                // that property, which is why both shapes are on the row.
                (
                    "extra_cycles_per_set_bit_select",
                    format!("{:.4}", select.cycles - walk.cycles),
                ),
                (
                    "extra_instructions_per_set_bit_select",
                    format!("{:.4}", select.instructions - walk.instructions),
                ),
                (
                    "branch_misses_per_set_bit_walk",
                    format!("{:.5}", walk.branch_misses),
                ),
                (
                    "branch_misses_per_set_bit_select",
                    format!("{:.5}", select.branch_misses),
                ),
                ("select_variant", variant.to_string()),
                (
                    "ns_per_set_bit_select_literal",
                    format!("{:.4}", literal.ns),
                ),
                (
                    "cycles_per_set_bit_select_literal",
                    format!("{:.4}", literal.cycles),
                ),
                (
                    "instructions_per_set_bit_select_literal",
                    format!("{:.4}", literal.instructions),
                ),
                (
                    "ns_per_set_bit_select_broadword",
                    format!("{:.4}", broadword.ns),
                ),
                (
                    "cycles_per_set_bit_select_broadword",
                    format!("{:.4}", broadword.cycles),
                ),
                (
                    "instructions_per_set_bit_select_broadword",
                    format!("{:.4}", broadword.instructions),
                ),
                (
                    "cycles_per_word_scan",
                    format!("{:.4}", row.scan.cycles / row.words as f64),
                ),
                (
                    "instructions_per_word_scan",
                    format!("{:.4}", row.scan.instructions / row.words as f64),
                ),
                ("scan_share_of_walk_pass", format!("{scan_share:.4}")),
                (
                    "ns_per_set_bit_walk_whole_pass",
                    format!("{:.4}", row.walk.ns / row.set_bits as f64),
                ),
                (
                    "ns_per_set_bit_select_whole_pass",
                    format!("{:.4}", select_pass.ns / row.set_bits as f64),
                ),
                ("whole_pass_ratio", format!("{whole_pass_ratio:.4}")),
                ("walk_share_of_extraction", format!("{walk_share:.6}")),
                ("c1_ceiling", format!("{:.6}", (1.0 - walk_share).recip())),
                (
                    "extraction_cycles_per_cell",
                    format!("{:.4}", row.extraction.cycles / row.cells as f64),
                ),
                ("target_feature_popcnt", target_feature_popcnt().to_string()),
                (
                    "instructions_per_count_ones",
                    format!("{:.4}", count_ones_price.instructions),
                ),
                (
                    "cycles_per_count_ones",
                    format!("{:.4}", count_ones_price.cycles),
                ),
                (
                    "count_ones_price_share",
                    format!("{count_ones_price_share:.4}"),
                ),
                (
                    "count_ones_price_readable",
                    count_ones_price_readable.to_string(),
                ),
                ("count_ones_per_set_bit_walk", String::from("0")),
                (
                    "count_ones_per_set_bit_select_literal",
                    format!("{literal_count_ones_calls:.4}"),
                ),
                ("count_ones_per_set_bit_select_broadword", String::from("0")),
                (
                    "count_ones_per_set_bit_select",
                    if variant == "broadword" {
                        String::from("0")
                    } else {
                        format!("{literal_count_ones_calls:.4}")
                    },
                ),
                ("published_pdep_penalty", PDEP_PENALTY.to_string()),
                ("select_ranks_verified", verified.to_string()),
                ("order_identical", row.order_identical.to_string()),
                ("c1_holds", c1_holds.to_string()),
                ("c1_holds_cycles", c1_holds_cycles.to_string()),
                ("cycle_marginal_valid", cycle_marginal_valid.to_string()),
                ("c2_holds", row.order_identical.to_string()),
                ("c3_holds", c3_holds.to_string()),
            ]);
        }
    }
}

fn main() {
    if !std::env::args().any(|arg| arg == "--bench") {
        return;
    }

    let prereg = isomesh::experiment!("P-108");

    #[cfg(target_os = "linux")]
    common::experiment::run(prereg, experiment::run);

    // `experiment_p12`'s precedent. Every clause on this row is a cycle or an
    // instruction count, and off Linux there is nothing to degrade to: a
    // recorded zero would be a fabricated measurement.
    #[cfg(not(target_os = "linux"))]
    {
        eprintln!(
            "{} compares two set-bit walks with hardware performance counters, and this \
             platform has no `perf_event_open`. M-281 forbids a nanosecond carrying the verdict.",
            prereg.id
        );
        std::process::exit(1);
    }
}
