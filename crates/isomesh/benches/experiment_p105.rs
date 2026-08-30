//! **P-105 — Harley–Seal carry-save popcount for a per-chunk active count the crate does not yet compute.**
//!
//! Ticket: R-105. Pre-registered before this harness existed.
//!
//! ```bash
//! cargo bench --bench experiment_p105
//! ```
//!
//! Writes `docs/experiments/p-105.csv`. **Linux only.** The verdict is
//! denominated in counters (`✗24`, `M-281`), the counters are
//! `perf_event_open`, and off Linux this bench `exit(1)`s rather than recording
//! a zero for a column it could not measure.
//!
//! # What was missing
//!
//! **Nothing in the shipped crate popcounts the bitmap, and this harness checks
//! that rather than asserting it.** `dual.rs` — the file that owns the
//! active-cell bitmap — contains **zero** occurrences of `count_ones`, which is
//! the column `count_ones_in_dual_rs`, read out of the source file on every run.
//! Its five consumers of the bitmap read whole words
//! ([`dual.rs:385`](../src/dual.rs) `inside_word`), read a word shifted by one
//! (`:395` `inside_word_shifted`), fold four rows into an active-cell word
//! (`:424` `active_word`), mask off the words that carry no cell (`:445`
//! `cell_mask`) and walk set bits (`:490`, `trailing_zeros` + `x &= x - 1`).
//! **There is no counting pass at all.** The eight `count_ones` elsewhere in
//! `crates/isomesh/src/` are over per-cell masks and not over the bitmap:
//! `trilinear.rs:295` counts a `u8` six-coordinate saddle mask, `hermite.rs:170`
//! a `u16` edge mask, `manifold_dual_contouring.rs:245` a `debug_assert!` on an
//! edge mask, and the remaining five are in test modules.
//!
//! So **C1's denominator is a bench-local naive count** and this row prices a
//! stage the crate would have to **add**. Nothing here is an extraction
//! speedup, and no clause claims one.
//!
//! # The build has no `POPCNT`, which is the regime that favours the mechanism
//!
//! `u64::count_ones()` lowers to LLVM's `ctpop.i64`. On a target **with** the
//! `popcnt` feature that is one instruction at four per cycle on Zen 3, and no
//! amount of carry-save arithmetic can beat it. The repository sets no
//! `target-cpu` and no `RUSTFLAGS`, so this bench is built for baseline
//! `x86_64` — **`POPCNT` is SSE4.2-era and is not in the baseline**, and
//! `objdump -d` on the bench binary finds **zero** `popcnt` instructions in it.
//! The column `target_feature_popcnt` records
//! `cfg!(target_feature = "popcnt")` so no reader has to guess which regime a
//! row was measured in: **a row measured under one value of that column is not
//! comparable with a row measured under the other.**
//!
//! This is the regime in which Harley–Seal has its best case, because the
//! popcounts it removes are expensive rather than free. It is worth naming what
//! actually happens to that best case, because the disassembly answers it and
//! the counters agree with the disassembly to three figures:
//!
//! | | per iteration | words | instructions per word |
//! |---|---|---|---|
//! | `count_naive` | 35 SSE2 instructions | 4 | **8.75** |
//! | `count_harley_seal` | 145 scalar instructions | 16 | **9.06** |
//!
//! LLVM **autovectorises the naive fold** — two `movdqu`, then twice the
//! `psrlw`/`pand`/`psubb`/`paddb` byte-population sequence closed by `psadbw`
//! and `paddq`, four `u64` per iteration — and does **not** autovectorise the
//! carry-save chain, whose nine live accumulators keep it in general-purpose
//! registers one word at a time. So the mechanism's saving (one popcount per
//! sixteen words instead of one per word) is competing against a saving LLVM
//! already took (two words per XMM lane pair, four per iteration), and the
//! comparison is between 8.75 and 9.06 instructions per word rather than
//! between five and twelve.
//!
//! `instructions_per_word_naive` and `instructions_per_word_hs` are therefore
//! not decoration: 8.7500 is 35/4 exactly and 9.0625 is 145/16 exactly, so
//! those two columns are the check that the counted window contains the loop it
//! claims to and nothing else.
//!
//! # The arms
//!
//! - **naive:** `words.iter().map(|w| w.count_ones()).sum()`, verbatim from the
//!   registration.
//! - **Harley–Seal:** the carry-save adder of Muła, Kurz & Lemire's AVX2
//!   popcount paper, taken as its **reduction** half and not its kernels —
//!   `csa(a, b, c)` is two XORs, two ANDs and one OR, no intrinsics, no
//!   `unsafe`, no target gating. Sixteen words per popcount is the canonical
//!   fold and the one used here: it is the variant that issues the **fewest**
//!   popcounts per word (one per sixteen against one per eight), so the
//!   mechanism is being given its best arithmetic rather than a handicapped one.
//!
//! Both arms are timed in **interleaved** counted windows — naive, Harley–Seal,
//! naive, … — so the two do not sit on opposite sides of a governor step, and
//! the median of five is taken per arm by wall time. `ghz` is on every row
//! beside the two `ns_per_word_*` columns for the same reason.
//!
//! **`ratio` is denominated in cycles, not in nanoseconds**, and that is the
//! one place this harness had to choose a reading of the registration. `✗24`
//! and `M-281` say a governed clock cannot carry a verdict, so `c1_holds` and
//! `c3_holds` read `cycles_per_word_naive / cycles_per_word_hs`. The wall-clock
//! form of the same quantity is beside it as `ns_ratio`, it is exactly the two
//! registered `ns_per_word_*` columns divided, and `ghz` says what the clock was
//! doing while they were measured. If the two forms ever disagreed about C1's
//! 2× bar the row would say so in `c1_disagreement` rather than leaving a reader
//! to find it by dividing.
//!
//! # The bitmap is the crate's, mirrored rather than approximated
//!
//! `dual.rs:359-381`'s `build_inside_bits` packs one bit per **sample**, 64 to a
//! `u64`, **along X only**, with `bit_row = size[0].div_ceil(64)`. The
//! active-cell bitmap this bench counts is then `active_word` (`:424`) masked by
//! `cell_mask` (`:445`), pushed in the extractor's own `z`, `y`, `w` order.
//!
//! The **cell row is one word shorter than the sample row** and that asymmetry
//! is real and documented at `dual.rs:472-484`; it is mirrored, not fixed. At
//! `resolution = 64` the sample row needs `bit_row = 2` words for 65 samples
//! while `cell_words = 1` covers all 64 cells, and the columns `bit_row` and
//! `cell_words` show it on the row.
//!
//! `resolution` counts **cells**, which is what makes the registration's own
//! arithmetic come out: 64³ cells is `1 × 64 × 64 = 4096` words, exactly the
//! *"64³ chunk's 32 KiB bitmap"* the sibling registration `P-107` prices a rank
//! directory against. Samples per axis are therefore 17, 33, 65, 129 — and 65
//! and 129 are the grids `dual.rs:476` names as *"every grid this crate
//! benchmarks at"*.
//!
//! # SHARE
//!
//! Every clause's reachable share is a **column**, and none of them is a share
//! of extraction, because the stage does not exist yet.
//!
//! - **C1's share of extraction is exactly zero, and `count_ones_in_dual_rs` is
//!   the column that says so.** `✗51`'s `1/(1 − share/factor)` bar therefore
//!   does not apply: there is no share to divide. What C1 *is* denominated in is
//!   `words` at `resolution = 64` — 4096 words, exact by construction — against
//!   a named bench-local baseline, which is why stating it as a ratio is
//!   legitimate. The deciding number is `c1_deciding_ratio` and the row it came
//!   from is `c1_deciding_resolution`.
//! - **C2's share is `words`: every word of every bitmap, on every row.** The
//!   equality is over the whole array, so its population is exact. What could
//!   make it *vacuous* is an all-zero bitmap — both arms return 0 and agree
//!   about nothing — so `active_cells` is a share column and is **asserted
//!   non-zero** on every row.
//! - **C3's share is the span of `words`: 256 → 32768, a factor of 128.** A
//!   monotonicity clause over four points that were all in L1 would be a
//!   statement about one cache level; `bitmap_bytes` runs 2 KiB → 256 KiB, which
//!   crosses out of L1 (32 KiB) and into L2, and both numbers are columns.
//!
//!   **C3 is monotonicity of *the advantage*, and it can be unscoreable rather
//!   than false.** Three conditions are tested from the run's own numbers, and
//!   any one of them makes `c3_vacuous` true and `c3_verdict` read `VACUOUS`:
//!
//!   1. **No advantage to be monotone in.** `advantage_exists` is `ratio_min >
//!      1`, the lowest of the row's five window pairs — an advantage that shows
//!      in every window, not one that shows in the median.
//!   2. **The trend is below the instrument's own resolution.**
//!      `c3_trend_span` is the spread of `ratio` across the four sizes and
//!      `c3_noise_span` is the widest single row's `ratio_spread` over its five
//!      window pairs. A trend narrower than one row's own reproducibility is not
//!      an ordering, and `R-085`'s discipline says an instrument reports the
//!      size of what it cannot resolve rather than ranking inside it.
//!   3. **Some consecutive pair was never ordered.** `c3_steps_separated`
//!      requires every consecutive pair's `[ratio_min, ratio_max]` to be
//!      **non-overlapping**, which is `E-307`'s own bar (`dual.rs:479-482`
//!      earns its conclusion from *"non-overlapping ranges"*) applied to this
//!      row's four points.
//!
//!   `monotone_step_holds` localises which of the three steps moved the wrong
//!   way, for the case where the trend *is* resolvable.
//! - **The vacuity control is `control_counts_differ`**, and it is asserted
//!   rather than reported: the same naive-against-Harley–Seal comparison is run
//!   against a copy of the bitmap with **one bit flipped**, and it must come out
//!   unequal by exactly one (`control_delta`). An equality asserted by an
//!   instrument that cannot see inequality is not evidence.
//! - **`tail_lengths_checked` closes the one hole the fixture leaves.** All four
//!   word counts are multiples of 16, so the fixture never exercises
//!   Harley–Seal's remainder path. Lengths 0..=64 of a SplitMix64 array are
//!   therefore checked against the naive count before the sweep, and
//!   `tail_mismatches` is asserted zero.

mod common;

#[cfg(target_os = "linux")]
mod experiment {
    use std::hint::black_box;
    use std::path::PathBuf;
    use std::time::Instant;

    use isomesh::fields::ReferenceField;
    use isomesh::{Sdf, Shape3, for_each_reference_field};

    use crate::common::counters::{MIN_TIME_RATIO, Probe};

    /// The scalar the bitmap's sign test reads.
    type Scalar = f64;

    /// Untimed windows first, per arm; the counted windows are steady state.
    const WARMUP_RUNS: u32 = 2;
    /// Counted windows per arm. The median by wall time is the one reported.
    const TIMED_RUNS: usize = 5;
    /// Words a counted window aims to process, by repeating the count.
    ///
    /// 2 KiB of bitmap is 256 words and a single fold of it is well under a
    /// microsecond — far too short a window for a counter read to mean
    /// anything. Repeating the count until about 16.7M words have gone past
    /// makes every row's window the same size in work rather than the same size
    /// in bitmaps, so `ns_per_word` and `instructions_per_word` are comparable
    /// across the four sizes, which is what C3 compares.
    const TARGET_WORDS: u64 = 1 << 24;

    /// **Cells** per axis. `resolution = 64` is the registration's *"64³
    /// chunk"*, and it is 4096 words of active-cell bitmap: 32 KiB exactly.
    const RESOLUTIONS: [u32; 4] = [16, 32, 64, 128];

    /// C1's bar.
    const C1_BAR: f64 = 2.0;
    /// The resolution C1 is decided at.
    const C1_RESOLUTION: u32 = 64;

    // ─── the arms ──────────────────────────────────────────────────────────

    /// The registration's own baseline: `iter().map(count_ones).sum()`.
    fn count_naive(words: &[u64]) -> u64 {
        words.iter().map(|w| u64::from(w.count_ones())).sum()
    }

    /// One carry-save adder: three bit-vectors in, sum and carry out.
    ///
    /// `(high, low)` where `low` is the parity of the three inputs and `high` is
    /// their majority. Two XORs, two ANDs, one OR — no intrinsics, no `unsafe`,
    /// no target feature.
    #[inline]
    fn csa(a: u64, b: u64, c: u64) -> (u64, u64) {
        let u = a ^ b;
        ((a & b) | (u & c), u ^ c)
    }

    /// Harley–Seal: sixteen words folded into one popcount.
    ///
    /// Fifteen [`csa`] calls per sixteen words carry the population up through
    /// weights 1, 2, 4, 8 and 16, so the loop issues **one** `count_ones` per
    /// sixteen words instead of sixteen. The four accumulators still holding
    /// weight when the loop ends are counted once each, and the tail — which
    /// this fixture never reaches, every word count being a multiple of 16 — is
    /// the naive fold. [`tail_control`] is what exercises it.
    fn count_harley_seal(words: &[u64]) -> u64 {
        let mut total = 0u64;
        let mut ones = 0u64;
        let mut twos = 0u64;
        let mut fours = 0u64;
        let mut eights = 0u64;

        // `as_chunks` rather than `chunks_exact`: the block is a `[u64; 16]`, so
        // the sixteen reads below are statically in bounds and no bounds check
        // survives into the loop that is being timed.
        let (blocks, remainder) = words.as_chunks::<16>();
        for block in blocks {
            let (twos_a, l) = csa(ones, block[0], block[1]);
            ones = l;
            let (twos_b, l) = csa(ones, block[2], block[3]);
            ones = l;
            let (fours_a, l) = csa(twos, twos_a, twos_b);
            twos = l;
            let (twos_a, l) = csa(ones, block[4], block[5]);
            ones = l;
            let (twos_b, l) = csa(ones, block[6], block[7]);
            ones = l;
            let (fours_b, l) = csa(twos, twos_a, twos_b);
            twos = l;
            let (eights_a, l) = csa(fours, fours_a, fours_b);
            fours = l;
            let (twos_a, l) = csa(ones, block[8], block[9]);
            ones = l;
            let (twos_b, l) = csa(ones, block[10], block[11]);
            ones = l;
            let (fours_a, l) = csa(twos, twos_a, twos_b);
            twos = l;
            let (twos_a, l) = csa(ones, block[12], block[13]);
            ones = l;
            let (twos_b, l) = csa(ones, block[14], block[15]);
            ones = l;
            let (fours_b, l) = csa(twos, twos_a, twos_b);
            twos = l;
            let (eights_b, l) = csa(fours, fours_a, fours_b);
            fours = l;
            let (sixteens, l) = csa(eights, eights_a, eights_b);
            eights = l;
            total += u64::from(sixteens.count_ones());
        }

        total *= 16;
        total += 8 * u64::from(eights.count_ones());
        total += 4 * u64::from(fours.count_ones());
        total += 2 * u64::from(twos.count_ones());
        total += u64::from(ones.count_ones());
        for word in remainder {
            total += u64::from(word.count_ones());
        }
        total
    }

    // ─── the bitmap, mirrored from `dual.rs` ───────────────────────────────

    /// `DualMesher`'s sample buffer: row stride `size[0] | 1` (`dual.rs:333`).
    ///
    /// The excess slot of an odd row is padding, left at zero exactly as
    /// `sample_grid` leaves it, and never read — `build_inside_bits` reads
    /// `min(sx - base, 64)` lanes per word.
    struct Grid {
        values: Vec<Scalar>,
        row: usize,
        size: [u32; 3],
    }

    impl Grid {
        fn sample<F: Sdf<Scalar = Scalar>>(
            field: &F,
            size: [u32; 3],
            origin: [Scalar; 3],
            cell_size: Scalar,
        ) -> Self {
            let row = size[0] as usize | 1;
            let rows = size[1] as usize * size[2] as usize;
            let mut values = vec![0.0; row * rows];
            for z in 0..size[2] {
                for y in 0..size[1] {
                    let base = row * (y as usize + size[1] as usize * z as usize);
                    for x in 0..size[0] {
                        values[base + x as usize] = field.sample([
                            origin[0] + cell_size * Scalar::from(x),
                            origin[1] + cell_size * Scalar::from(y),
                            origin[2] + cell_size * Scalar::from(z),
                        ]);
                    }
                }
            }
            Self { values, row, size }
        }
    }

    /// One bit per **sample**, 64 to a `u64`, along X only.
    ///
    /// `dual.rs:359-381`, copied rather than made `pub`: `dual` is a private
    /// module and `DualMesher` is `pub(crate)`, so there is no public seam and
    /// Phase 25 holds `crates/isomesh/src/**` read-only.
    struct Inside {
        bits: Vec<u64>,
        bit_row: usize,
        size: [u32; 3],
    }

    impl Inside {
        /// `dual.rs:359-381`, with `cube::is_inside`'s strict `value < 0`
        /// inlined — `cube` is a private module, and the predicate is one
        /// comparison whose exact form (zero is **outside**) matters.
        fn build(grid: &Grid) -> Self {
            let size = grid.size;
            let sx = size[0] as usize;
            let rows = size[1] as usize * size[2] as usize;
            let bit_row = sx.div_ceil(64);
            let mut bits = vec![0u64; bit_row * rows];

            for row in 0..rows {
                let src = grid.row * row;
                let dst = bit_row * row;
                for w in 0..bit_row {
                    let base = w * 64;
                    let n = (sx - base).min(64);
                    let mut word = 0u64;
                    for k in 0..n {
                        word |= u64::from(grid.values[src + base + k] < 0.0) << k;
                    }
                    bits[dst + w] = word;
                }
            }
            Self {
                bits,
                bit_row,
                size,
            }
        }

        /// `dual.rs:385`.
        #[inline]
        fn word(&self, w: usize, y: usize, z: usize) -> u64 {
            self.bits[self.bit_row * (y + self.size[1] as usize * z) + w]
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

        /// `dual.rs:424`: `any & !all` over the four `(y, z)` rows.
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
            any & !all
        }
    }

    /// `dual.rs:445`. `1u64 << 64` is undefined, so the full word is named.
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

    /// The active-cell bitmap: the array whose population count this row prices.
    struct Active {
        words: Vec<u64>,
        cell_words: usize,
        cells: usize,
    }

    impl Active {
        /// Built in the extractor's own `z`, `y`, `w` walk (`dual.rs:487-490`),
        /// masked by [`cell_mask`] so no bit is an answer about a sample that is
        /// not a cell.
        fn build(inside: &Inside) -> Self {
            let size = inside.size;
            let cells = [size[0] - 1, size[1] - 1, size[2] - 1];
            let cells_x = cells[0] as usize;
            let cell_words = cells_x.div_ceil(64);
            let mut words = Vec::with_capacity(cell_words * cells[1] as usize * cells[2] as usize);
            for z in 0..cells[2] as usize {
                for y in 0..cells[1] as usize {
                    for w in 0..cell_words {
                        words.push(inside.active(w, y, z) & cell_mask(w, cells_x));
                    }
                }
            }
            Self {
                words,
                cell_words,
                cells: cells[0] as usize * cells[1] as usize * cells[2] as usize,
            }
        }
    }

    // ─── the counted window ────────────────────────────────────────────────

    /// One arm's median counted window, normalised per word.
    struct Arm {
        ns_per_word: f64,
        cycles_per_word: f64,
        instructions_per_word: f64,
        ghz: f64,
        nanos: f64,
        cycles: f64,
    }

    /// One counted window: `reps` folds of the same array.
    ///
    /// `black_box(words)` inside the loop is load-bearing. Without it the fold
    /// is a pure function of an unchanged slice and LLVM hoists all `reps`
    /// iterations into one, which would time the loop counter.
    fn window(
        probe: &mut Probe,
        words: &[u64],
        reps: u64,
        fold: fn(&[u64]) -> u64,
    ) -> (u128, u64, u64) {
        probe.reset_and_enable();
        let started = Instant::now();
        let mut acc = 0u64;
        for _ in 0..reps {
            acc = acc.wrapping_add(fold(black_box(words)));
        }
        let nanos = started.elapsed().as_nanos();
        probe.disable();
        let counted = probe.read();
        black_box(acc);
        assert!(
            counted.worst_ratio() >= MIN_TIME_RATIO,
            "a counter ran only {:.1}% of the time it was enabled, so its value is an \
             extrapolation rather than a measurement",
            counted.worst_ratio() * 100.0
        );
        (nanos, counted.cycles.count, counted.instructions.count)
    }

    fn median(mut runs: Vec<(u128, u64, u64)>, words: u64, reps: u64) -> Arm {
        runs.sort_unstable();
        let (nanos, cycles, instructions) = runs[runs.len() / 2];
        let per = (words * reps) as f64;
        let nanos = nanos as f64;
        let cycles = cycles as f64;
        Arm {
            ns_per_word: nanos / per,
            cycles_per_word: cycles / per,
            instructions_per_word: instructions as f64 / per,
            ghz: cycles / nanos,
            nanos,
            cycles,
        }
    }

    /// Both arms, in interleaved counted windows.
    ///
    /// Interleaved because the two ratios this row reports are ratios *between*
    /// the arms, and a governor step between a block of naive windows and a
    /// block of Harley–Seal windows would land entirely in the ratio. `ghz` is
    /// recorded per arm as well as pooled so a step that happened anyway is
    /// visible rather than absorbed.
    ///
    /// Also returns the **per-window** cycle ratio, one per interleaved pair.
    ///
    /// A median ratio with no spread beside it cannot say whether a 1% trend
    /// across the four sizes is a trend, and C3 is a clause about exactly such a
    /// trend. Both windows of a pair fold the same `words * reps`, so the
    /// per-word normalisation cancels and the raw cycle counts divide directly.
    fn measure(words: &[u64], reps: u64) -> (Arm, Arm, Vec<f64>) {
        for _ in 0..WARMUP_RUNS {
            black_box(count_naive(black_box(words)));
            black_box(count_harley_seal(black_box(words)));
        }
        let mut probe = Probe::open();
        let mut naive = Vec::with_capacity(TIMED_RUNS);
        let mut hs = Vec::with_capacity(TIMED_RUNS);
        let mut pairs = Vec::with_capacity(TIMED_RUNS);
        for _ in 0..TIMED_RUNS {
            let n = window(&mut probe, words, reps, count_naive);
            let h = window(&mut probe, words, reps, count_harley_seal);
            pairs.push(n.1 as f64 / h.1 as f64);
            naive.push(n);
            hs.push(h);
        }
        let n = words.len() as u64;
        (median(naive, n, reps), median(hs, n, reps), pairs)
    }

    // ─── the controls that run before the sweep ────────────────────────────

    /// SplitMix64, so the tail control is byte-identical on every machine.
    fn splitmix(state: &mut u64) -> u64 {
        *state = state.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = *state;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }

    /// Every length 0..=64, because the fixture only has multiples of 16.
    ///
    /// All four registered word counts are 256, 1024, 4096 and 32768, so
    /// `chunks_exact(16)` never leaves a remainder and Harley–Seal's tail is
    /// dead code on every row of the CSV. A fold whose tail is wrong would then
    /// pass C2 on this fixture and corrupt any future caller. Returns the number
    /// of mismatches, which is asserted zero.
    fn tail_control() -> (u32, u32) {
        let mut state = 0x0000_2026u64;
        let buffer: Vec<u64> = (0..64).map(|_| splitmix(&mut state)).collect();
        let mut mismatches = 0;
        for len in 0..=64usize {
            if count_naive(&buffer[..len]) != count_harley_seal(&buffer[..len]) {
                mismatches += 1;
            }
        }
        (65, mismatches)
    }

    /// How many times `count_ones` appears in `dual.rs`.
    ///
    /// The registration's premise, read off the source rather than asserted in
    /// prose. Zero is the claim: the file that owns the bitmap never counts it.
    fn count_ones_in_dual_rs() -> usize {
        let path = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("src/dual.rs");
        let text =
            std::fs::read_to_string(&path).unwrap_or_else(|e| panic!("{}: {e}", path.display()));
        text.matches("count_ones").count()
    }

    // ─── one row ───────────────────────────────────────────────────────────

    struct Row {
        resolution: u32,
        samples: u32,
        words: usize,
        bit_row: usize,
        cell_words: usize,
        cells: usize,
        count: u64,
        count_hs: u64,
        counts_equal: bool,
        control_delta: i64,
        naive: Arm,
        hs: Arm,
        /// One cycle ratio per interleaved window pair, in window order.
        pairs: Vec<f64>,
    }

    impl Row {
        /// The verdict form: cycles, not the clock (`✗24`, `M-281`).
        fn ratio(&self) -> f64 {
            self.naive.cycles_per_word / self.hs.cycles_per_word
        }

        fn ns_ratio(&self) -> f64 {
            self.naive.ns_per_word / self.hs.ns_per_word
        }

        fn instruction_ratio(&self) -> f64 {
            self.naive.instructions_per_word / self.hs.instructions_per_word
        }

        /// The lowest and highest of the per-window ratios.
        ///
        /// The row's own reproducibility, measured rather than assumed. C3 asks
        /// for a trend across four sizes; a trend smaller than this is not a
        /// trend, and `c3_vacuous` is decided on exactly that comparison.
        fn ratio_bounds(&self) -> (f64, f64) {
            self.pairs
                .iter()
                .fold((f64::MAX, f64::MIN), |(lo, hi), &r| (lo.min(r), hi.max(r)))
        }

        fn ratio_spread(&self) -> f64 {
            let (lo, hi) = self.ratio_bounds();
            hi - lo
        }
    }

    fn one<F>(field: &F, resolution: u32) -> Row
    where
        F: ReferenceField + Sdf<Scalar = Scalar>,
    {
        // `resolution` counts cells, so the grid is one sample wider per axis.
        let samples = resolution + 1;
        let (shape, origin, cell_size) = crate::common::grid::<Scalar, F>(field, samples);
        let size = shape.size();

        let grid = Grid::sample(field, size, origin, cell_size);
        let inside = Inside::build(&grid);
        let active = Active::build(&inside);

        let words = active.words.len();
        assert!(words > 0, "the active-cell bitmap has no words");
        let reps = (TARGET_WORDS / words as u64).max(1);

        let count = count_naive(&active.words);
        let count_hs = count_harley_seal(&active.words);
        assert!(
            count > 0,
            "no active cells at {resolution}³, so both arms return 0 and C2's equality is vacuous"
        );

        // VACUITY CONTROL. The same naive-against-Harley–Seal comparison,
        // against a copy with one bit flipped. It must come out unequal, and by
        // exactly one, or the equality above is being asserted by an instrument
        // that cannot see inequality.
        let mut corrupt = active.words.clone();
        let victim = corrupt.len() / 2;
        corrupt[victim] ^= 1;
        let control = count_harley_seal(&corrupt);
        let control_delta = control as i64 - count as i64;
        assert_eq!(
            control_delta.abs(),
            1,
            "flipping one bit of the bitmap moved the compared count by {control_delta}, so the \
             comparison is not reading the array it is given"
        );

        let (naive, hs, pairs) = measure(&active.words, reps);

        Row {
            resolution,
            samples,
            words,
            bit_row: inside.bit_row,
            cell_words: active.cell_words,
            cells: active.cells,
            count,
            count_hs,
            counts_equal: count == count_hs,
            control_delta,
            naive,
            hs,
            pairs,
        }
    }

    // ─── the run ───────────────────────────────────────────────────────────

    pub(crate) fn run(run: &mut crate::common::experiment::Run) {
        let popcnt = cfg!(target_feature = "popcnt");
        let dual_popcounts = count_ones_in_dual_rs();
        let (tail_lengths, tail_mismatches) = tail_control();
        assert_eq!(
            tail_mismatches, 0,
            "Harley–Seal disagrees with the naive count on {tail_mismatches} of {tail_lengths} \
             tail lengths, and the fixture's word counts are all multiples of 16 so the CSV would \
             never have seen it"
        );

        println!(
            "count_ones in crates/isomesh/src/dual.rs: {dual_popcounts}  \
             (the registration's premise: the file that owns the bitmap never counts it)"
        );
        println!(
            "target_feature = \"popcnt\": {popcnt}  \
             (false means count_ones is the SWAR expansion, not one instruction)"
        );
        println!("tail lengths 0..=64 checked: {tail_lengths}, mismatches {tail_mismatches}\n");

        println!(
            "{:>15} {:>5} {:>7} {:>9} {:>8} {:>8} {:>7} {:>7} {:>7} {:>7} {:>6} {:>9}",
            "field",
            "cells",
            "words",
            "active",
            "ns/w n",
            "ns/w hs",
            "ins n",
            "ins hs",
            "ratio",
            "spread",
            "GHz",
            "C1/C3"
        );

        for_each_reference_field!(f64, |name, field| {
            let rows: Vec<Row> = RESOLUTIONS
                .iter()
                .map(|&resolution| one(&field, resolution))
                .collect();

            // C1 is decided at 64³ cells, which is the registration's own
            // "64³ chunk" and this fixture's 4096-word / 32 KiB row.
            let deciding = rows
                .iter()
                .find(|r| r.resolution == C1_RESOLUTION)
                .expect("RESOLUTIONS contains the resolution C1 is decided at");
            let c1_ratio = deciding.ratio();
            let c1_holds = c1_ratio >= C1_BAR;
            // If the clock and the counters disagreed about the bar, the row
            // says so rather than letting the reader find out by dividing the
            // two ns columns.
            let c1_disagreement = (deciding.ns_ratio() >= C1_BAR) != c1_holds;

            // C3: non-decreasing in the word count, exactly, no tolerance —
            // RESOLUTIONS is ascending and so is `words`.
            let mut steps = Vec::with_capacity(rows.len());
            let mut previous = f64::NEG_INFINITY;
            for row in &rows {
                let step = row.ratio() >= previous;
                previous = row.ratio();
                steps.push(step);
            }
            let c3_holds = steps.iter().all(|&s| s);

            // C3 is monotonicity of **the advantage**, and three things can make
            // that unscoreable rather than merely false. All three are decided
            // from this run's own numbers and all three are columns.
            //
            //  1. No row has an advantage that survives its own five windows,
            //     so there is nothing for the clause to be monotone in.
            //  2. The trend the clause asks about is smaller than the row-level
            //     reproducibility of the same quantity, so no ordering of the
            //     four points is a finding. `R-085`'s discipline applied to a
            //     trend instead of a residual: an instrument reports the size of
            //     what it cannot resolve.
            //  3. Some consecutive pair's `[ratio_min, ratio_max]` windows
            //     **overlap**, so which of those two sizes is faster was not
            //     established at all. `E-307` (`dual.rs:479-482`) is the
            //     precedent for requiring non-overlapping ranges before a pair
            //     of measurements is allowed to be ordered.
            let advantage_anywhere = rows.iter().any(|r| r.ratio_bounds().0 > 1.0);
            let trend_lo = rows.iter().map(Row::ratio).fold(f64::MAX, f64::min);
            let trend_hi = rows.iter().map(Row::ratio).fold(f64::MIN, f64::max);
            let c3_trend_span = trend_hi - trend_lo;
            let c3_noise_span = rows.iter().map(Row::ratio_spread).fold(f64::MIN, f64::max);
            let c3_steps_separated = rows.windows(2).all(|pair| {
                let (a_lo, a_hi) = pair[0].ratio_bounds();
                let (b_lo, b_hi) = pair[1].ratio_bounds();
                a_hi < b_lo || b_hi < a_lo
            });
            let c3_vacuous =
                !advantage_anywhere || c3_trend_span <= c3_noise_span || !c3_steps_separated;

            let c1_verdict = if c1_holds { "HELD" } else { "FALSIFIED" };
            let c3_verdict = if c3_vacuous {
                "VACUOUS"
            } else if c3_holds {
                "HELD"
            } else {
                "FALSIFIED"
            };

            for (row, &step) in rows.iter().zip(steps.iter()) {
                let active_fraction = row.count as f64 / row.cells as f64;
                println!(
                    "{name:>15} {:>5} {:>7} {:>9} {:>8.3} {:>8.3} {:>7.2} {:>7.2} {:>7.3} \
                     {:>7.3} {:>6.3} {:>9}",
                    row.resolution,
                    row.words,
                    row.count,
                    row.naive.ns_per_word,
                    row.hs.ns_per_word,
                    row.naive.instructions_per_word,
                    row.hs.instructions_per_word,
                    row.ratio(),
                    row.ratio_spread(),
                    (row.naive.cycles + row.hs.cycles) / (row.naive.nanos + row.hs.nanos),
                    format!("{}/{}", &c1_verdict[..1], &c3_verdict[..1])
                );
                run.record(&[
                    ("field", name.to_string()),
                    ("resolution", row.resolution.to_string()),
                    ("words", row.words.to_string()),
                    ("count_naive", row.count.to_string()),
                    ("count_harley_seal", row.count_hs.to_string()),
                    ("counts_equal", u8::from(row.counts_equal).to_string()),
                    ("ns_per_word_naive", format!("{:.6}", row.naive.ns_per_word)),
                    ("ns_per_word_hs", format!("{:.6}", row.hs.ns_per_word)),
                    ("ratio", format!("{:.6}", row.ratio())),
                    (
                        "instructions_per_word_naive",
                        format!("{:.4}", row.naive.instructions_per_word),
                    ),
                    (
                        "instructions_per_word_hs",
                        format!("{:.4}", row.hs.instructions_per_word),
                    ),
                    (
                        "control_counts_differ",
                        u8::from(row.control_delta != 0).to_string(),
                    ),
                    ("c1_holds", u8::from(c1_holds).to_string()),
                    ("c2_holds", u8::from(row.counts_equal).to_string()),
                    ("c3_holds", u8::from(c3_holds).to_string()),
                    // ─── extras (M-273) ───
                    ("scalar", "f64".to_string()),
                    ("samples_per_axis", row.samples.to_string()),
                    ("cells", row.cells.to_string()),
                    ("active_cells", row.count.to_string()),
                    ("active_fraction", format!("{active_fraction:.6}")),
                    ("bitmap_bytes", (row.words * 8).to_string()),
                    ("bit_row", row.bit_row.to_string()),
                    ("cell_words", row.cell_words.to_string()),
                    ("reps", (TARGET_WORDS / row.words as u64).max(1).to_string()),
                    (
                        "cycles_per_word_naive",
                        format!("{:.4}", row.naive.cycles_per_word),
                    ),
                    (
                        "cycles_per_word_hs",
                        format!("{:.4}", row.hs.cycles_per_word),
                    ),
                    ("ns_ratio", format!("{:.6}", row.ns_ratio())),
                    (
                        "instruction_ratio",
                        format!("{:.6}", row.instruction_ratio()),
                    ),
                    (
                        "ipc_naive",
                        format!(
                            "{:.4}",
                            row.naive.instructions_per_word / row.naive.cycles_per_word
                        ),
                    ),
                    (
                        "ipc_hs",
                        format!(
                            "{:.4}",
                            row.hs.instructions_per_word / row.hs.cycles_per_word
                        ),
                    ),
                    (
                        "ghz",
                        format!(
                            "{:.4}",
                            (row.naive.cycles + row.hs.cycles) / (row.naive.nanos + row.hs.nanos)
                        ),
                    ),
                    ("ghz_naive", format!("{:.4}", row.naive.ghz)),
                    ("ghz_hs", format!("{:.4}", row.hs.ghz)),
                    ("control_delta", row.control_delta.to_string()),
                    ("monotone_step_holds", u8::from(step).to_string()),
                    (
                        "advantage_exists",
                        u8::from(row.ratio_bounds().0 > 1.0).to_string(),
                    ),
                    ("ratio_min", format!("{:.6}", row.ratio_bounds().0)),
                    ("ratio_max", format!("{:.6}", row.ratio_bounds().1)),
                    ("ratio_spread", format!("{:.6}", row.ratio_spread())),
                    ("windows_per_arm", TIMED_RUNS.to_string()),
                    ("c3_trend_span", format!("{c3_trend_span:.6}")),
                    ("c3_noise_span", format!("{c3_noise_span:.6}")),
                    (
                        "c3_steps_separated",
                        u8::from(c3_steps_separated).to_string(),
                    ),
                    ("c3_vacuous", u8::from(c3_vacuous).to_string()),
                    ("c3_verdict", c3_verdict.to_string()),
                    ("c1_verdict", c1_verdict.to_string()),
                    (
                        "c2_verdict",
                        if row.counts_equal {
                            "HELD"
                        } else {
                            "FALSIFIED"
                        }
                        .to_string(),
                    ),
                    ("c1_deciding_ratio", format!("{c1_ratio:.6}")),
                    ("c1_deciding_resolution", C1_RESOLUTION.to_string()),
                    ("c1_disagreement", u8::from(c1_disagreement).to_string()),
                    ("target_feature_popcnt", u8::from(popcnt).to_string()),
                    ("count_ones_in_dual_rs", dual_popcounts.to_string()),
                    ("tail_lengths_checked", tail_lengths.to_string()),
                    ("tail_mismatches", tail_mismatches.to_string()),
                ]);
            }
        });
    }
}

fn main() {
    if !std::env::args().any(|arg| arg == "--bench") {
        return;
    }
    let prereg = isomesh::experiment!("P-105");

    #[cfg(target_os = "linux")]
    common::experiment::run(prereg, experiment::run);

    #[cfg(not(target_os = "linux"))]
    {
        // Refuse rather than record a zero: `instructions_per_word_*` is the
        // verdict form of C1 and C3, and there is no `perf_event_open` here.
        eprintln!(
            "{} needs hardware performance counters, and this platform has no `perf_event_open`.",
            prereg.id
        );
        std::process::exit(1);
    }
}
