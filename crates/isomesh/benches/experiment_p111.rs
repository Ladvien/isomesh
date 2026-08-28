//! **P-111 — table-driven scalar compaction, 8 cells per lookup, branchless.**
//!
//! Ticket: R-111. Pre-registered before this harness existed.
//!
//! ```bash
//! cargo bench --bench experiment_p111
//! ```
//!
//! Writes `docs/experiments/p-111.csv`. **Linux only**, for `experiment_p12`'s
//! reason: C3 is an L1-miss clause, the instrument is `perf_event_open`, and a
//! recorded zero off Linux would be a fabricated cache measurement.
//!
//! # What was missing
//!
//! simdjson's compaction is a 256-entry table mapping an 8-bit mask to a byte
//! permutation, fed to `_mm256_shuffle_epi8`. Reduced to scalar the SIMD half
//! evaporates and what is left is the *table*: `[[u8; 8]; 256]`, the set-bit
//! positions of a byte, 2,048 bytes, `const`-evaluable, safe, deterministic,
//! and reachable with no `_pext_u64`.
//!
//! **That is not a stylistic preference, it is the only admissible form here,
//! and the census is verifiable rather than asserted.** `crates/isomesh/src/`
//! today contains **zero `unsafe` blocks, zero `core::simd` imports and zero
//! `#[cfg(target_arch)]` gates** — the only three matches for those strings
//! anywhere under `src/` are English prose inside `experiment.rs`'s own
//! registrations and one sentence in `subgrid/extract/tests.rs`. `P-69`'s
//! registration records why the `core::simd` door stays shut: *"core::simd is
//! nightly and is staying nightly"*, `rust-lang/portable-simd#364` unresolved
//! and the maintainers' own 2025 summary being *"nightly-only and will remain
//! such"*. `_pext_u64` needs `unsafe` **and** a target gate **and** is a
//! microcoded 18-cycle trap on pre-Zen-3 AMD. A `const fn`-built lookup table
//! needs none of the three. So this row asks the only version of simdjson's
//! question this crate could ever act on.
//!
//! What was missing is the *comparison*. Nothing in the repository had measured
//! a table-driven compaction against `dual.rs:489-497`'s set-bit walk, and
//! nothing had measured the walk's own cost per set bit at all — `M-337`
//! measured the bitmap prepass as a whole, `P-40` measured active-cell
//! enumeration against the eight-corner scalar predicate, and neither
//! decomposes the enumeration from the word scan that feeds it.
//!
//! # The incumbent, mirrored rather than described
//!
//! `crates/isomesh/src/**` is read-only for Phase 25, so the bitmap is rebuilt
//! bench-local from `dual.rs:359-381` (`build_inside_bits`), `:385`
//! (`inside_word`), `:395` (`inside_word_shifted`), `:424` (`active_word`),
//! `:445` (`cell_mask`) and `:484` (`cell_words`, the cell/sample word
//! asymmetry `E-307` paid for). The sign rule is **not** copied: the mirror
//! calls the crate's own public
//! [`isomesh::marching_cubes::table::is_inside`], so the bit can never drift
//! from the extractor's.
//!
//! The incumbent enumeration is `dual.rs:495-497` verbatim:
//!
//! ```text
//! while active != 0 {
//!     let x = w * 64 + active.trailing_zeros() as usize;
//!     active &= active - 1;
//!     ...
//! }
//! ```
//!
//! and `dual.rs:491-494` says out loud why it is shaped that way — clearing the
//! lowest set bit visits the row in ascending `x`, *"the same order the scalar
//! loop did, which is what keeps vertex creation order and therefore every
//! index unchanged"*. That is C2's whole content: an enumeration that produces
//! the same set in a different order is not a drop-in for this loop, it is a
//! different mesh.
//!
//! # The four arms, and why there are four
//!
//! Both mechanisms live *inside* the word loop, so timing "the walk" means
//! timing a pass that also scans every cell word. Denominating per set bit —
//! which the registration does, so the 97% of words neither arm touches cannot
//! inflate the ratio — does not by itself remove the **shared** word scan, and
//! a shared cost pulls any ratio toward 1.0. So the shared part is measured and
//! subtracted rather than argued about:
//!
//! - **`scan`** — [`Bits::sweep_scan`]: `active_word & cell_mask` per word, ORed
//!   into an accumulator, no enumeration. The floor both arms stand on.
//! - **`walk`** — [`Bits::sweep_walk`]: the incumbent, `trailing_zeros` +
//!   `active &= active - 1`.
//! - **`table`** — [`Bits::sweep_table`]: the mechanism as registered. On a
//!   non-zero word, eight unconditional `SET_BITS[byte]` lookups, each writing
//!   **eight** slots and advancing the cursor by `byte.count_ones()`. That is
//!   simdjson's own trick — write more than you need, advance by the count — and
//!   it is what "branchless" means here.
//! - **`table_byteskip`** — [`Bits::sweep_table_byteskip`]: identical but with
//!   `if byte != 0` around the eight stores. Not registered; recorded because it
//!   is the one column that can distinguish *the table is a bad idea* from *the
//!   unconditional stores are a bad idea*, and a null that cannot name its cause
//!   is worth much less than one that can.
//!
//! Every arm writes into a pre-sized `Vec<u32>` through a `len` cursor, never
//! `Vec::push`. The table form needs eight slots of slack past the true length —
//! that is the mechanism — and giving the walk `push` while the table gets a
//! cursor would measure `Vec`'s growth policy instead of the enumeration. The
//! buffers are `cells + 8` long, so `out[len + 7]` is always in bounds and the
//! slack is safe indexing rather than `set_len`. **No `unsafe` anywhere in this
//! file**, which is the point of the row.
//!
//! Windows are **siblings, never nested**: `R-121` paid to discover that Zen 3
//! has six general-purpose counters, `Probe` opens exactly six plus a software
//! event, and two nested windows multiplex until [`MIN_TIME_RATIO`] refuses.
//! `enumeration` is therefore a **prefix difference** — `walk − scan`,
//! `table − scan` — of two independent windows over the same grid in the same
//! repetition, not a window inside a window.
//!
//! # SHARE
//!
//! The stage this row competes in is `cycles_emit_walk`: the set-bit walk plus
//! the per-cell indexing it feeds. `P-121` measured it at 65³ (wave-time
//! numbers, taken on a machine running four concurrent bench agents, so a prior
//! rather than a published measurement — cited as such and re-run on a clean
//! tree by the phase lead). Its share of **`marching_cubes`** extraction, `f32`:
//!
//! ```text
//! gyroid 0.0601   sphere 0.0585   box_exact 0.0469   torus 0.0428
//! csg_difference 0.0426   noise_cavity 0.0383   thin_plate 0.0169
//! fbm_terrain 0.0029
//! ```
//!
//! and of **`dual_contouring`** extraction, `f32`, 0.0002–0.0109. So the
//! Amdahl ceiling on *eliminating this stage entirely* is `1/(1 − 0.0601)` =
//! **1.064×** at its very best, and `1.0003×` on `fbm_terrain` under
//! `dual_contouring`.
//!
//! Each clause's reachable share is therefore a **column**, per field, so no
//! reader can mistake a per-set-bit ratio for an extraction speedup:
//!
//! - **`p121_emit_walk_share_mc_65`** and **`p121_emit_walk_share_dc_65`** —
//!   that field's own measured share.
//! - **`extraction_ceiling_mc`** and **`extraction_ceiling_dc`** —
//!   `1/(1 − share)`, the most a *perfect* replacement could ever return.
//! - **`extraction_speedup_mc`** — what this row's own measured `ratio` would
//!   actually buy: `1/(1 − share·(1 − 1/ratio))`, which is the honest number and
//!   is below 1 whenever the table loses.
//!
//! `fbm_terrain`'s 0.0029 is small because field evaluation is 94% of extraction
//! there, **not** because the walk got cheaper; its absolute per-set-bit cost is
//! in line with the others, which is exactly why this row's clauses are
//! per-set-bit and its ceiling is a separate column.
//!
//! **The clauses themselves are denominated per set bit and remain scoreable
//! whatever the share is.** C1 is a per-set-bit comparison, C2 is an equality
//! and C3 is a per-set-bit cache count. `✗51`'s rule is that a *speedup* clause
//! must know its share before the harness is written; it was read first, it is
//! on every row, and it bounds the *consequence* of the result rather than the
//! result.
//!
//! # Which form carries the verdict
//!
//! `M-280` / `M-281`: on a governed CPU a nanosecond is not a unit. `R-105`
//! measured the identical binary's cycle-ratio band drifting 0.984 → 1.035
//! across three runs while its instruction counts held to four figures. So:
//!
//! - The registered `ratio` is the **ns** form, because
//!   `ns_per_set_bit_{walk,table}` are the registered columns and a registration
//!   is not amended to suit the harness. `c1_holds` is scored on it.
//! - `cycle_ratio` and `instruction_ratio` are beside it, and
//!   `c1_holds_cycles` / `c1_holds_instructions` / `c1_agreement` say whether
//!   the three forms agree. **`instruction_ratio` is the reproducible one and it
//!   is the form the write-up should quote**; a row where `c1_agreement` is
//!   false is a marginal verdict and must be reported as one.
//! - `stores_per_set_bit_walk` and `stores_per_set_bit_table` are **analytic
//!   integers** derived from the bitmap's own `set_bits`, `words_nonzero` and
//!   `bytes_nonzero` — machine-independent, run-independent, and the actual
//!   mechanism behind whatever the ratio turns out to be. They are the column a
//!   reader should check before any timing.
//! - `ghz` is on every row because the row carries ns columns, and no clause
//!   consults it.
//!
//! **`ratio` is `ns_per_set_bit_walk / ns_per_set_bit_table`, so above 1 means
//! the table wins.** The direction is on every row as `ratio_definition`,
//! because an inverted ratio is the easiest way to publish the opposite of what
//! was measured.
//!
//! **This is not a hypothetical, and the harness caught it.** Two consecutive
//! runs of the identical binary, minutes apart, read `thin_plate 129³`'s walk
//! arm at **26.32** and then **43.74** ns per set bit — a 1.66× move on the
//! same code over the same grid — and `c1_holds` therefore read `false` and
//! then `true` for that row. Over the same two runs `instruction_ratio` read
//! **0.9487** and **0.9488**. Six bench agents were running concurrently on
//! this machine, which is precisely the condition Phase 24's clean-tree serial
//! re-run exists for. So the clock's own band is a column too:
//! `ns_per_set_bit_{walk,table}_rep_{min,max}` over [`REPS`] repetitions, and
//! `ratio_from_rep_minima` with `c1_holds_from_rep_minima` beside it. Wall-clock
//! contention can only ever make a pass *slower*, so the fastest repetition of
//! each arm is the least disturbed pair and the ratio of minima is the ns form
//! least able to invent a verdict. [`REPS`] is 15 rather than 9 for the same
//! reason. **None of this replaces the registered median form** — it is
//! recorded beside it so a reader can see when the clock was lying.
//!
//! C3 gets the same treatment for the same reason, and it needs it more. The
//! registered falsifier is *"the table costing more L1 than it saves"*, and
//! `c3_holds` is scored on the literal reading with **zero tolerance**:
//! `l1_misses_per_set_bit_table <= l1_misses_per_set_bit_walk`. That is the
//! registration unamended. But two hardware counter readings differing by half
//! a percent are not a verdict about anything, so the instrument's own
//! resolution is measured and reported beside it:
//! `l1_misses_per_set_bit_walk_rep_min` / `_rep_max` are the walk arm's band
//! over [`REPS`] repetitions, `l1_delta_per_set_bit_table_minus_walk` is the
//! signed difference, and `c3_holds_within_instrument_spread` asks whether the
//! table's median lands inside the walk's own run-to-run band. **Both scorings
//! are columns and neither is hidden**; a write-up must say which it quotes and
//! report the other.
//!
//! # VACUITY CONTROL
//!
//! Asserted, and the assertions cover the **population** rather than the
//! outcome — a guard that asserts the answer cannot report the bad news:
//!
//! - `set_bits > 0` on every row, so a per-set-bit denominator is never a
//!   division by an empty population.
//! - `words_nonzero > 0` and `bytes_nonzero > 0`, so both arms actually enter
//!   their inner loops.
//! - `table_bytes == 2048`, because "2 KiB" is part of the registered claim.
//! - Off Linux, `exit(1)`. `benches/common/mod.rs:26` is the single
//!   `cfg(target_os = "linux")` gate and `experiment_p12` is the
//!   refuse-and-exit precedent.
//!
//! `order_identical`, `c1_holds`, `c2_holds` and `c3_holds` are **recorded and
//! never asserted**: they are the clauses, and a harness that aborts on a
//! falsified clause deletes its own result.

#![allow(
    clippy::cast_precision_loss,
    reason = "counts become f64 ratios; the counts are far below 2^53"
)]

mod common;

#[cfg(target_os = "linux")]
mod experiment {
    use std::hint::black_box;
    use std::time::Instant;

    use crate::common::counters::{MIN_TIME_RATIO, Probe};
    use isomesh::marching_cubes::table::is_inside;

    /// Samples per axis. The registered fixture, and both straddle the word
    /// boundary: 65 samples is 64 cells (one cell word, two sample words) and
    /// 129 is 128 cells (two cell words, three sample words), so the
    /// cell/sample asymmetry at `dual.rs:484` is exercised rather than assumed.
    const SIZES: [u32; 2] = [65, 129];

    /// Repetitions medianed per quantity.
    ///
    /// Fifteen rather than nine: this machine ran six concurrent bench agents
    /// while the row was written, and a median has to survive several disturbed
    /// repetitions rather than one. The `_rep_min` / `_rep_max` columns report
    /// what it survived.
    const REPS: usize = 15;

    /// Target length of one counter window, so the ~28 `perf_event` system
    /// calls a window costs land outside it and cannot inflate an arm.
    const TARGET_BATCH_NS: f64 = 20_000_000.0;

    /// Cap on the batch, so a cheap row cannot spin for minutes.
    const MAX_INNER: usize = 8192;

    /// The mechanism: the set-bit positions of every byte value.
    ///
    /// `SET_BITS[b][i]` is the index of the `i`-th lowest set bit of `b`, and
    /// the slots past `b.count_ones()` are **zero padding that is deliberately
    /// written and never read**. That is the branchless contract: the caller
    /// stores all eight slots at `out[len..len + 8]` and advances `len` by
    /// `b.count_ones()`, so the next byte's stores overwrite the padding and
    /// only `out[..len]` is ever consumed.
    static SET_BITS: [[u8; 8]; 256] = build_set_bits();

    /// The registered "2 KiB", computed rather than claimed.
    const TABLE_BYTES: usize = size_of::<[[u8; 8]; 256]>();

    /// Build [`SET_BITS`] at compile time.
    ///
    /// `const fn`, so the table costs nothing at run time and nothing in the
    /// binary beyond its 2,048 bytes — the whole reason this is the admissible
    /// form of simdjson's mechanism in a crate with no `unsafe`.
    const fn build_set_bits() -> [[u8; 8]; 256] {
        let mut table = [[0u8; 8]; 256];
        let mut byte = 0usize;
        while byte < 256 {
            let mut bit = 0usize;
            let mut slot = 0usize;
            while bit < 8 {
                if byte & (1 << bit) != 0 {
                    table[byte][slot] = bit as u8;
                    slot += 1;
                }
                bit += 1;
            }
            byte += 1;
        }
        table
    }

    /// `DualMesher`'s active-cell bitmap, mirrored bench-local.
    ///
    /// One bit per **sample**, 64 to a `u64`, along `x` only — `dual.rs:359`,
    /// and it is samples rather than cells, which is the correction `P-104`'s
    /// registration carries. `row` is `size[0] | 1`, `DualMesher::row_stride`
    /// (`dual.rs:333`) exactly, so the mirror's addressing is the crate's.
    struct Bits {
        inside: Vec<u64>,
        bit_row: usize,
        size: [u32; 3],
        cells: [usize; 3],
    }

    impl Bits {
        fn build<F>(field: &F, lo: [f64; 3], h: f64, n: u32) -> Self
        where
            F: isomesh::Sdf<Scalar = f64>,
        {
            let sx = n as usize;
            let row = sx | 1;
            let rows = sx * sx;
            let mut values = vec![0.0f64; row * rows];
            for z in 0..sx {
                for y in 0..sx {
                    let base = row * (y + sx * z);
                    for x in 0..sx {
                        values[base + x] = field.sample([
                            lo[0] + h * x as f64,
                            lo[1] + h * y as f64,
                            lo[2] + h * z as f64,
                        ]);
                    }
                }
            }

            let bit_row = sx.div_ceil(64);
            let mut inside = vec![0u64; bit_row * rows];
            for r in 0..rows {
                let src = row * r;
                let dst = bit_row * r;
                for w in 0..bit_row {
                    let base = w * 64;
                    let count = (sx - base).min(64);
                    let mut word = 0u64;
                    for k in 0..count {
                        // Branchless, and the predicate is the crate's own
                        // rather than a copy of `value < 0`.
                        word |= u64::from(is_inside(values[src + base + k])) << k;
                    }
                    inside[dst + w] = word;
                }
            }

            Self {
                inside,
                bit_row,
                size: [n; 3],
                cells: [sx - 1; 3],
            }
        }

        #[inline]
        fn word(&self, w: usize, y: usize, z: usize) -> u64 {
            self.inside[self.bit_row * (y + self.size[1] as usize * z) + w]
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

        /// `dual.rs:424`. Sixty-four active-cell answers, `any & !all`.
        #[inline]
        fn active_word(&self, w: usize, y: usize, z: usize) -> u64 {
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

        /// `dual.rs:445`. `1u64 << 64` is undefined, so the full word is named.
        #[inline]
        fn cell_mask(w: usize, cells_x: usize) -> u64 {
            let remaining = cells_x.saturating_sub(w * 64);
            if remaining >= 64 {
                !0
            } else {
                (1u64 << remaining) - 1
            }
        }

        /// `dual.rs:484`. Words carrying a **cell**, one fewer than `bit_row`
        /// whenever `size[0]` is `64k + 1` — which both fixture sizes are.
        #[inline]
        fn cell_words(&self) -> usize {
            self.cells[0].div_ceil(64)
        }

        fn cell_count(&self) -> usize {
            self.cells[0] * self.cells[1] * self.cells[2]
        }

        /// The floor both arms stand on: the word scan, with no enumeration.
        fn sweep_scan(&self) -> u64 {
            let cells_x = self.cells[0];
            let words = self.cell_words();
            let mut acc = 0u64;
            for z in 0..self.cells[2] {
                for y in 0..self.cells[1] {
                    for w in 0..words {
                        acc |= self.active_word(w, y, z) & Self::cell_mask(w, cells_x);
                    }
                }
            }
            acc
        }

        /// The incumbent, `dual.rs:495-497`.
        fn sweep_walk(&self, out: &mut [u32]) -> usize {
            let cells_x = self.cells[0];
            let words = self.cell_words();
            let mut len = 0usize;
            for z in 0..self.cells[2] {
                for y in 0..self.cells[1] {
                    let row_base = (cells_x * (y + self.cells[1] * z)) as u32;
                    for w in 0..words {
                        let mut active = self.active_word(w, y, z) & Self::cell_mask(w, cells_x);
                        while active != 0 {
                            let x = (w * 64) as u32 + active.trailing_zeros();
                            active &= active - 1;
                            out[len] = row_base + x;
                            len += 1;
                        }
                    }
                }
            }
            len
        }

        /// simdjson's table, reduced to scalar: eight cells per lookup,
        /// branchless inside a non-zero word.
        fn sweep_table(&self, out: &mut [u32]) -> usize {
            let cells_x = self.cells[0];
            let words = self.cell_words();
            let mut len = 0usize;
            for z in 0..self.cells[2] {
                for y in 0..self.cells[1] {
                    let row_base = (cells_x * (y + self.cells[1] * z)) as u32;
                    for w in 0..words {
                        let active = self.active_word(w, y, z) & Self::cell_mask(w, cells_x);
                        // The one test both arms share: the walk's
                        // `while active != 0` is the same branch.
                        if active == 0 {
                            continue;
                        }
                        let word_base = row_base + (w * 64) as u32;
                        let mut rest = active;
                        for b in 0..8usize {
                            let byte = (rest & 0xff) as usize;
                            let base = word_base + (b * 8) as u32;
                            let positions = &SET_BITS[byte];
                            for k in 0..8 {
                                out[len + k] = base + u32::from(positions[k]);
                            }
                            len += byte.count_ones() as usize;
                            rest >>= 8;
                        }
                    }
                }
            }
            len
        }

        /// [`Bits::sweep_table`] with the zero bytes branched over.
        ///
        /// Not registered. It is the column that separates "the table is wrong"
        /// from "the unconditional stores are wrong".
        fn sweep_table_byteskip(&self, out: &mut [u32]) -> usize {
            let cells_x = self.cells[0];
            let words = self.cell_words();
            let mut len = 0usize;
            for z in 0..self.cells[2] {
                for y in 0..self.cells[1] {
                    let row_base = (cells_x * (y + self.cells[1] * z)) as u32;
                    for w in 0..words {
                        let active = self.active_word(w, y, z) & Self::cell_mask(w, cells_x);
                        if active == 0 {
                            continue;
                        }
                        let word_base = row_base + (w * 64) as u32;
                        let mut rest = active;
                        for b in 0..8usize {
                            let byte = (rest & 0xff) as usize;
                            if byte != 0 {
                                let base = word_base + (b * 8) as u32;
                                let positions = &SET_BITS[byte];
                                for k in 0..8 {
                                    out[len + k] = base + u32::from(positions[k]);
                                }
                                len += byte.count_ones() as usize;
                            }
                            rest >>= 8;
                        }
                    }
                }
            }
            len
        }

        /// The analytic census: set bits, non-zero words, non-zero bytes.
        ///
        /// Exact integers, so `stores_per_set_bit_*` is a machine-independent
        /// property of the fixture rather than a timing.
        fn census(&self) -> (u64, u64, u64) {
            let cells_x = self.cells[0];
            let words = self.cell_words();
            let mut set_bits = 0u64;
            let mut nonzero_words = 0u64;
            let mut nonzero_bytes = 0u64;
            for z in 0..self.cells[2] {
                for y in 0..self.cells[1] {
                    for w in 0..words {
                        let active = self.active_word(w, y, z) & Self::cell_mask(w, cells_x);
                        set_bits += u64::from(active.count_ones());
                        if active != 0 {
                            nonzero_words += 1;
                            for b in 0..8u32 {
                                if (active >> (b * 8)) & 0xff != 0 {
                                    nonzero_bytes += 1;
                                }
                            }
                        }
                    }
                }
            }
            (set_bits, nonzero_words, nonzero_bytes)
        }
    }

    /// One counter window's three readings.
    #[derive(Clone, Copy, Default)]
    struct Counted {
        cycles: f64,
        instructions: f64,
        l1d_read_misses: f64,
        nanos: f64,
    }

    impl Counted {
        fn minus(self, other: Self) -> Self {
            Self {
                cycles: self.cycles - other.cycles,
                instructions: self.instructions - other.instructions,
                l1d_read_misses: self.l1d_read_misses - other.l1d_read_misses,
                nanos: self.nanos - other.nanos,
            }
        }
    }

    /// One counter window over `inner` repetitions of `body`, divided by
    /// `inner`.
    ///
    /// The `perf_event` system calls are outside the counted region. Windows
    /// are opened one after another and **never nested** — `R-121` established
    /// that two nested windows multiplex on a six-counter machine and
    /// [`MIN_TIME_RATIO`] then refuses, which is asserted here rather than
    /// hoped for.
    fn window(probe: &mut Probe, inner: usize, mut body: impl FnMut()) -> Counted {
        probe.reset_and_enable();
        let started = Instant::now();
        for _ in 0..inner {
            body();
        }
        let nanos = started.elapsed().as_nanos() as f64;
        probe.disable();
        let counts = probe.read();
        assert!(
            counts.worst_ratio() >= MIN_TIME_RATIO,
            "a counter ran only {:.1}% of the time it was enabled, so its value is an \
             extrapolation rather than a measurement",
            counts.worst_ratio() * 100.0
        );
        let scale = 1.0 / inner as f64;
        Counted {
            cycles: counts.cycles.count as f64 * scale,
            instructions: counts.instructions.count as f64 * scale,
            l1d_read_misses: counts.l1d_read_misses.count as f64 * scale,
            nanos: nanos * scale,
        }
    }

    /// The four arms of one repetition.
    #[derive(Clone, Copy, Default)]
    struct Rep {
        scan: Counted,
        walk: Counted,
        table: Counted,
        byteskip: Counted,
    }

    fn median(pick: &dyn Fn(&Rep) -> f64, reps: &[Rep]) -> f64 {
        let mut values: Vec<f64> = reps.iter().map(pick).collect();
        values.sort_by(f64::total_cmp);
        values[values.len() / 2]
    }

    fn median_counted(pick: &dyn Fn(&Rep) -> Counted, reps: &[Rep]) -> Counted {
        Counted {
            cycles: median(&|r| pick(r).cycles, reps),
            instructions: median(&|r| pick(r).instructions, reps),
            l1d_read_misses: median(&|r| pick(r).l1d_read_misses, reps),
            nanos: median(&|r| pick(r).nanos, reps),
        }
    }
    /// The extreme values of one quantity across the repetitions.
    ///
    /// The instrument's own resolution, measured rather than assumed. A
    /// difference between two arms that lies inside one arm's run-to-run band
    /// is a difference this instrument cannot see, and `R-105` is why that has
    /// to be a column: the identical binary's cycle-ratio band drifted 0.984 →
    /// 1.035 across three runs.
    fn spread(pick: &dyn Fn(&Rep) -> f64, reps: &[Rep]) -> (f64, f64) {
        let values: Vec<f64> = reps.iter().map(pick).collect();
        (
            values.iter().copied().fold(f64::INFINITY, f64::min),
            values.iter().copied().fold(f64::NEG_INFINITY, f64::max),
        )
    }

    /// Everything one `(field, resolution)` produced.
    struct Row {
        field: &'static str,
        resolution: u32,
        cells: usize,
        set_bits: u64,
        nonzero_words: u64,
        nonzero_bytes: u64,
        order_identical: bool,
        byteskip_identical: bool,
        scan: Counted,
        walk: Counted,
        table: Counted,
        byteskip: Counted,
        /// `(min, max)` of the walk arm's L1D read misses per pass, over
        /// [`REPS`] repetitions. C3's resolution.
        walk_l1_spread: (f64, f64),
        /// `(min, max)` of the table arm's, the same way.
        table_l1_spread: (f64, f64),
        /// `(min, max)` of the walk arm's nanoseconds per pass, over [`REPS`]
        /// repetitions.
        ///
        /// The clock's only failure mode on a contended machine is
        /// **inflation** — nothing makes a pass finish faster than it can — so
        /// the minimum is the least contaminated reading and the maximum is the
        /// contamination itself. Both are columns. This is not theoretical: two
        /// consecutive runs of this binary read `thin_plate 129³` at 26.3 and
        /// then 43.7 ns per set bit on the walk arm and **flipped `c1_holds`**,
        /// while `instruction_ratio` held at 0.9487 and 0.9488.
        walk_ns_spread: (f64, f64),
        /// `(min, max)` of the table arm's, the same way.
        table_ns_spread: (f64, f64),
        inner: usize,
    }

    /// `P-121`'s wave-time `cycles_emit_walk / cycles_total` at 65³, `f32`.
    ///
    /// Taken on a machine running four concurrent bench agents, so a **prior**
    /// rather than a published measurement, and cited as such. It bounds the
    /// consequence of this row's result; it is not part of any clause.
    fn p121_shares(field: &str) -> (f64, f64) {
        match field {
            "sphere" => (0.0585, 0.0066),
            "torus" => (0.0428, 0.0058),
            "box_exact" => (0.0469, 0.0109),
            "csg_difference" => (0.0426, 0.0076),
            "thin_plate" => (0.0169, 0.0048),
            "gyroid" => (0.0601, 0.0054),
            "fbm_terrain" => (0.0029, 0.0002),
            "noise_cavity" => (0.0383, 0.0031),
            other => panic!("P-121 published no emit_walk share for `{other}`"),
        }
    }

    /// `1/(1 − share)`: the most eliminating the stage entirely could return.
    fn ceiling(share: f64) -> f64 {
        1.0 / (1.0 - share)
    }

    /// What a mechanism `ratio`× on the stage buys on whole extraction.
    fn realised(share: f64, ratio: f64) -> f64 {
        1.0 / (1.0 - share * (1.0 - 1.0 / ratio))
    }

    fn measure<F>(field_name: &'static str, field: &F, n: u32) -> Row
    where
        F: isomesh::Sdf<Scalar = f64> + isomesh::fields::ReferenceField,
    {
        let (_shape, lo, h) = crate::common::grid::<f64, F>(field, n);
        let bits = Bits::build(field, lo, h, n);
        let cells = bits.cell_count();
        let (set_bits, nonzero_words, nonzero_bytes) = bits.census();

        // ── the population, asserted; the clauses, recorded ──────────────────
        assert!(
            set_bits > 0,
            "{field_name} {n}³: no active cell, so a per-set-bit denominator would be empty"
        );
        assert!(
            nonzero_words > 0 && nonzero_bytes > 0,
            "{field_name} {n}³: neither arm would enter its inner loop"
        );
        assert_eq!(
            TABLE_BYTES, 2048,
            "the registered claim is a 2 KiB table and this one is not"
        );

        // ── C2, once, outside every window ──────────────────────────────────
        let slots = cells + 8;
        let mut walk_out = vec![0u32; slots];
        let mut table_out = vec![0u32; slots];
        let mut skip_out = vec![0u32; slots];
        let walk_len = bits.sweep_walk(&mut walk_out);
        let table_len = bits.sweep_table(&mut table_out);
        let skip_len = bits.sweep_table_byteskip(&mut skip_out);
        assert_eq!(
            walk_len as u64, set_bits,
            "{field_name} {n}³: the walk visited {walk_len} cells and the census counted {set_bits}"
        );
        let order_identical =
            table_len == walk_len && table_out[..table_len] == walk_out[..walk_len];
        let byteskip_identical =
            skip_len == walk_len && skip_out[..skip_len] == walk_out[..walk_len];

        // ── batch size, from an untimed pass ─────────────────────────────────
        let probe_pass = Instant::now();
        for _ in 0..3 {
            black_box(bits.sweep_walk(&mut walk_out));
        }
        let pass_ns = probe_pass.elapsed().as_nanos() as f64 / 3.0;
        let inner = ((TARGET_BATCH_NS / pass_ns.max(1.0)).ceil() as usize).clamp(1, MAX_INNER);

        // ── REPS repetitions of four sibling windows ─────────────────────────
        let mut probe = Probe::open();
        let mut reps: Vec<Rep> = Vec::with_capacity(REPS);
        for _ in 0..REPS {
            let scan = window(&mut probe, inner, || {
                black_box(bits.sweep_scan());
            });
            let walk = window(&mut probe, inner, || {
                let len = bits.sweep_walk(&mut walk_out);
                black_box(len);
                black_box(&walk_out);
            });
            let table = window(&mut probe, inner, || {
                let len = bits.sweep_table(&mut table_out);
                black_box(len);
                black_box(&table_out);
            });
            let byteskip = window(&mut probe, inner, || {
                let len = bits.sweep_table_byteskip(&mut skip_out);
                black_box(len);
                black_box(&skip_out);
            });
            reps.push(Rep {
                scan,
                walk,
                table,
                byteskip,
            });
        }

        Row {
            field: field_name,
            resolution: n,
            cells,
            set_bits,
            nonzero_words,
            nonzero_bytes,
            order_identical,
            byteskip_identical,
            scan: median_counted(&|r| r.scan, &reps),
            walk: median_counted(&|r| r.walk, &reps),
            table: median_counted(&|r| r.table, &reps),
            byteskip: median_counted(&|r| r.byteskip, &reps),
            walk_l1_spread: spread(&|r| r.walk.l1d_read_misses, &reps),
            walk_ns_spread: spread(&|r| r.walk.nanos, &reps),
            table_ns_spread: spread(&|r| r.table.nanos, &reps),
            table_l1_spread: spread(&|r| r.table.l1d_read_misses, &reps),
            inner,
        }
    }

    pub(crate) fn run(run: &mut crate::common::experiment::Run) {
        let mut rows: Vec<Row> = Vec::new();
        for n in SIZES {
            isomesh::for_each_reference_field!(f64, |name, field| {
                // Inline block, so no `return` in here (M-253).
                rows.push(measure(name, &field, n));
            });
        }

        println!(
            "{:<16} {:>4} {:>8} {:>7} {:>8} {:>8} {:>7} {:>7} {:>7} {:>7} {:>6} {:>6} {:>6}",
            "field",
            "n",
            "act_frac",
            "setbits",
            "ns/bit_w",
            "ns/bit_t",
            "ratio",
            "rat_min",
            "cyc_r",
            "instr_r",
            "st/bit",
            "c1",
            "c3"
        );

        for r in &rows {
            let set_bits = r.set_bits as f64;
            let active_fraction = set_bits / r.cells as f64;

            // Whole-pass, per set bit. This is the registered denomination: the
            // enumeration lives inside the word loop and the pass is what the
            // extractor actually runs.
            let ns_walk = r.walk.nanos / set_bits;
            let ns_table = r.table.nanos / set_bits;
            let ns_skip = r.byteskip.nanos / set_bits;
            let ratio = ns_walk / ns_table;

            let cyc_walk = r.walk.cycles / set_bits;
            let cyc_table = r.table.cycles / set_bits;
            let cycle_ratio = cyc_walk / cyc_table;

            let ins_walk = r.walk.instructions / set_bits;
            let ins_table = r.table.instructions / set_bits;
            let ins_skip = r.byteskip.instructions / set_bits;
            let instruction_ratio = ins_walk / ins_table;

            let l1_walk = r.walk.l1d_read_misses / set_bits;
            let l1_table = r.table.l1d_read_misses / set_bits;
            let l1_skip = r.byteskip.l1d_read_misses / set_bits;

            // The shared word scan, subtracted: prefix differences of sibling
            // windows, never a nested window.
            let enum_walk = r.walk.minus(r.scan);
            let enum_table = r.table.minus(r.scan);
            let enum_ns_walk = enum_walk.nanos / set_bits;
            let enum_ns_table = enum_table.nanos / set_bits;
            let enum_ins_walk = enum_walk.instructions / set_bits;
            let enum_ins_table = enum_table.instructions / set_bits;

            // Analytic, exact, machine-independent: the walk stores once per
            // set bit; the table stores eight slots for each of the eight bytes
            // of every non-zero word.
            let stores_walk = 1.0;
            let stores_table = 64.0 * r.nonzero_words as f64 / set_bits;
            let stores_byteskip = 8.0 * r.nonzero_bytes as f64 / set_bits;

            // C3's decisive form: the shared word scan's L1 misses subtracted
            // from both arms, so what is left is the **enumeration's own** cache
            // cost. The 2 KiB table is 32 cache lines in a 32 KiB L1D, so the
            // mechanism's prediction is that this is ~0 and that any whole-pass
            // delta is ambient — on a loaded machine a slower arm simply spends
            // longer being evicted by the other cores, which is a property of
            // the machine and not of the table.
            let l1_scan = r.scan.l1d_read_misses / set_bits;
            let l1_enum_walk = (r.walk.l1d_read_misses - r.scan.l1d_read_misses) / set_bits;
            let l1_enum_table = (r.table.l1d_read_misses - r.scan.l1d_read_misses) / set_bits;

            let (share_mc, share_dc) = p121_shares(r.field);

            let c1_holds = ratio > 1.0;
            let c1_cycles = cycle_ratio > 1.0;
            let c1_instructions = instruction_ratio > 1.0;
            let c2_holds = r.order_identical;
            // "The table's L1 cost does not eat the win", scored **literally**:
            // the table must not read-miss L1 more per set bit than the walk
            // does. Zero tolerance, because a registration is not amended to
            // suit its harness.
            let c3_holds = l1_table <= l1_walk;
            // And the same question at the instrument's own measured
            // resolution, because a 0.5% difference between two hardware
            // counter readings is not a verdict about anything. The walk arm's
            // L1 misses per set bit over REPS repetitions span
            // `walk_l1_spread`; if the table's median lands inside that band,
            // this instrument cannot distinguish the two arms and the strict
            // score above is a coin flip rather than a measurement. **Both are
            // recorded and neither is hidden** — the CSV lets a reader score C3
            // either way, and the write-up must quote whichever it uses.
            let l1_walk_rep_max = r.walk_l1_spread.1 / set_bits;
            let l1_walk_rep_min = r.walk_l1_spread.0 / set_bits;
            let l1_table_rep_max = r.table_l1_spread.1 / set_bits;
            let l1_table_rep_min = r.table_l1_spread.0 / set_bits;
            let c3_within_spread = l1_table <= l1_walk_rep_max;

            // The clock, uncontaminated as far as a contended machine allows.
            // Contention only ever *inflates* a wall-clock reading, so the
            // fastest repetition of each arm is the least disturbed pair, and a
            // ratio of minima is the ns form least able to invent a verdict.
            // Recorded beside the registered median form, never instead of it.
            let ns_walk_rep_min = r.walk_ns_spread.0 / set_bits;
            let ns_walk_rep_max = r.walk_ns_spread.1 / set_bits;
            let ns_table_rep_min = r.table_ns_spread.0 / set_bits;
            let ns_table_rep_max = r.table_ns_spread.1 / set_bits;
            let ratio_minima = ns_walk_rep_min / ns_table_rep_min;

            println!(
                "{:<16} {:>4} {:>8.5} {:>7} {:>8.3} {:>8.3} {:>7.4} {:>7.4} {:>7.4} {:>7.4} \
                 {:>6.2} {:>6} {:>6}",
                r.field,
                r.resolution,
                active_fraction,
                r.set_bits,
                ns_walk,
                ns_table,
                ratio,
                ratio_minima,
                cycle_ratio,
                instruction_ratio,
                stores_table,
                c1_holds,
                c3_holds
            );

            run.record(&[
                // ── registered ──────────────────────────────────────────────
                ("field", r.field.to_string()),
                ("resolution", r.resolution.to_string()),
                ("active_fraction", format!("{active_fraction:.7}")),
                ("set_bits", r.set_bits.to_string()),
                ("ns_per_set_bit_walk", format!("{ns_walk:.5}")),
                ("ns_per_set_bit_table", format!("{ns_table:.5}")),
                ("ratio", format!("{ratio:.5}")),
                ("l1_misses_per_set_bit_walk", format!("{l1_walk:.5}")),
                ("l1_misses_per_set_bit_table", format!("{l1_table:.5}")),
                ("order_identical", r.order_identical.to_string()),
                ("table_bytes", TABLE_BYTES.to_string()),
                ("c1_holds", c1_holds.to_string()),
                ("c2_holds", c2_holds.to_string()),
                ("c3_holds", c3_holds.to_string()),
                // ── the forms a verdict can be read in (M-273) ──────────────
                ("cycles_per_set_bit_walk", format!("{cyc_walk:.5}")),
                ("cycles_per_set_bit_table", format!("{cyc_table:.5}")),
                ("cycle_ratio", format!("{cycle_ratio:.5}")),
                ("instructions_per_set_bit_walk", format!("{ins_walk:.5}")),
                ("instructions_per_set_bit_table", format!("{ins_table:.5}")),
                ("instruction_ratio", format!("{instruction_ratio:.5}")),
                ("c1_holds_cycles", c1_cycles.to_string()),
                ("c1_holds_instructions", c1_instructions.to_string()),
                (
                    "c1_agreement",
                    (c1_holds == c1_cycles && c1_holds == c1_instructions).to_string(),
                ),
                // ── the shared word scan, subtracted ────────────────────────
                (
                    "ns_per_set_bit_scan",
                    format!("{:.5}", r.scan.nanos / set_bits),
                ),
                (
                    "instructions_per_set_bit_scan",
                    format!("{:.5}", r.scan.instructions / set_bits),
                ),
                (
                    "ns_per_set_bit_enumeration_walk",
                    format!("{enum_ns_walk:.5}"),
                ),
                (
                    "ns_per_set_bit_enumeration_table",
                    format!("{enum_ns_table:.5}"),
                ),
                (
                    "enumeration_ratio",
                    format!("{:.5}", enum_ns_walk / enum_ns_table),
                ),
                (
                    "instructions_per_set_bit_enumeration_walk",
                    format!("{enum_ins_walk:.5}"),
                ),
                (
                    "instructions_per_set_bit_enumeration_table",
                    format!("{enum_ins_table:.5}"),
                ),
                (
                    "enumeration_instruction_ratio",
                    format!("{:.5}", enum_ins_walk / enum_ins_table),
                ),
                // ── the byte-skipping diagnostic arm ────────────────────────
                ("ns_per_set_bit_table_byteskip", format!("{ns_skip:.5}")),
                ("ratio_byteskip", format!("{:.5}", ns_walk / ns_skip)),
                (
                    "instructions_per_set_bit_table_byteskip",
                    format!("{ins_skip:.5}"),
                ),
                (
                    "instruction_ratio_byteskip",
                    format!("{:.5}", ins_walk / ins_skip),
                ),
                (
                    "l1_misses_per_set_bit_table_byteskip",
                    format!("{l1_skip:.5}"),
                ),
                ("byteskip_order_identical", r.byteskip_identical.to_string()),
                // ── the analytic mechanism, exact integers ──────────────────
                ("cells", r.cells.to_string()),
                ("words_nonzero", r.nonzero_words.to_string()),
                ("bytes_nonzero", r.nonzero_bytes.to_string()),
                (
                    "cell_words_per_row",
                    r.resolution.saturating_sub(1).div_ceil(64).to_string(),
                ),
                ("stores_per_set_bit_walk", format!("{stores_walk:.5}")),
                ("stores_per_set_bit_table", format!("{stores_table:.5}")),
                (
                    "stores_per_set_bit_table_byteskip",
                    format!("{stores_byteskip:.5}"),
                ),
                // ── the SHARE ceiling, as a column, per field ──────────────
                ("p121_emit_walk_share_mc_65", format!("{share_mc:.4}")),
                ("p121_emit_walk_share_dc_65", format!("{share_dc:.4}")),
                ("extraction_ceiling_mc", format!("{:.5}", ceiling(share_mc))),
                ("extraction_ceiling_dc", format!("{:.5}", ceiling(share_dc))),
                (
                    "extraction_speedup_mc",
                    format!("{:.5}", realised(share_mc, ratio)),
                ),
                (
                    "extraction_speedup_dc",
                    format!("{:.5}", realised(share_dc, ratio)),
                ),
                // ── C3 at the instrument's own resolution ──────────────────
                (
                    "l1_misses_per_set_bit_walk_rep_min",
                    format!("{l1_walk_rep_min:.5}"),
                ),
                (
                    "l1_misses_per_set_bit_walk_rep_max",
                    format!("{l1_walk_rep_max:.5}"),
                ),
                (
                    "l1_misses_per_set_bit_table_rep_min",
                    format!("{l1_table_rep_min:.5}"),
                ),
                (
                    "l1_misses_per_set_bit_table_rep_max",
                    format!("{l1_table_rep_max:.5}"),
                ),
                (
                    "l1_delta_per_set_bit_table_minus_walk",
                    format!("{:.6}", l1_table - l1_walk),
                ),
                (
                    "c3_holds_within_instrument_spread",
                    c3_within_spread.to_string(),
                ),
                ("l1_misses_per_set_bit_scan", format!("{l1_scan:.5}")),
                (
                    "l1_misses_per_set_bit_enumeration_walk",
                    format!("{l1_enum_walk:.5}"),
                ),
                (
                    "l1_misses_per_set_bit_enumeration_table",
                    format!("{l1_enum_table:.5}"),
                ),
                (
                    "l1_enumeration_delta_per_set_bit",
                    format!("{:.6}", l1_enum_table - l1_enum_walk),
                ),
                // ── the clock's own spread, and the ratio it cannot inflate ─
                (
                    "ns_per_set_bit_walk_rep_min",
                    format!("{ns_walk_rep_min:.5}"),
                ),
                (
                    "ns_per_set_bit_walk_rep_max",
                    format!("{ns_walk_rep_max:.5}"),
                ),
                (
                    "ns_per_set_bit_table_rep_min",
                    format!("{ns_table_rep_min:.5}"),
                ),
                (
                    "ns_per_set_bit_table_rep_max",
                    format!("{ns_table_rep_max:.5}"),
                ),
                ("ratio_from_rep_minima", format!("{ratio_minima:.5}")),
                ("c1_holds_from_rep_minima", (ratio_minima > 1.0).to_string()),
                // ── provenance ─────────────────────────────────────────────
                (
                    "ratio_definition",
                    "ns_per_set_bit_walk/ns_per_set_bit_table; above 1 means the table wins"
                        .to_string(),
                ),
                (
                    "l1_miss_ratio_table_over_walk",
                    format!("{:.5}", l1_table / l1_walk.max(f64::MIN_POSITIVE)),
                ),
                ("ghz", format!("{:.4}", r.walk.cycles / r.walk.nanos)),
                ("inner_reps", r.inner.to_string()),
                ("reps", REPS.to_string()),
            ]);
        }
    }
}

fn main() {
    if !std::env::args().any(|arg| arg == "--bench") {
        return;
    }

    let prereg = isomesh::experiment!("P-111");

    #[cfg(target_os = "linux")]
    common::experiment::run(prereg, experiment::run);

    // `experiment_p12`'s precedent, and the registration names it: C3 is an
    // L1-miss clause and the instrument is `perf_event_open`. Off Linux there is
    // nothing to degrade to, and a recorded zero would be a fabricated cache
    // measurement rather than a missing one.
    #[cfg(not(target_os = "linux"))]
    {
        eprintln!(
            "{} scores C3 on L1 read misses per set bit from hardware performance counters, and \
             this platform has no `perf_event_open`. A zero here would be an invention.",
            prereg.id
        );
        std::process::exit(1);
    }
}
