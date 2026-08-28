//! **P-106 — SWAR sign extraction and edge-crossing masks, exhaustively over all 256 patterns.**
//!
//! Ticket: R-106. Pre-registered before this harness existed.
//!
//! ```bash
//! cargo bench --bench experiment_p106
//! ```
//!
//! Writes `docs/experiments/p-106.csv`. **Linux only**, for `experiment_p12`'s
//! reason: C3 is an instructions-per-cell comparison and the only instrument
//! that can read one is `perf_event_open`. A nanosecond on a governed CPU is not
//! a unit (`✗24`, `M-280`, `M-281`), so off Linux there is nothing to degrade
//! to and the harness refuses and exits 1 rather than record a fabricated zero.
//!
//! # Provenance correction this row carries
//!
//! Fujita's *Bitwise Parallel Bulk Computation* — the paper the technique is
//! named after — **did not download** in `M-415`'s acquisition. The technique is
//! therefore registered against the papers that *did* land, FastLanes
//! (`Scalar_T64`: `u64` registers as quasi-SIMD, no intrinsics, no `unsafe`) and
//! the AVX2-popcount paper. And the figure that gets quoted for it — **13.4 ×
//! 10⁹ cell-updates per second** — is **quoted from
//! `docs/research/2026-08-28-bitpacking-simd-acquisition-and-backlog.md`, not
//! from any paper in the corpus.** It is not a comparand, nothing here is scored
//! against it, and no number below is denominated in it.
//!
//! # What was missing
//!
//! An edge of a cell is cut exactly when its two corners differ in sign. That is
//! an XOR, and `cube.rs:26`'s `EDGE_CORNERS` is the pair list — so the twelve
//! cut-edge flags of a cell are twelve XORs of eight bits, and sixty-four cells'
//! worth are twelve XORs of eight **words**. Nothing in the repository had
//! written that circuit down, counted its operations, or checked it against the
//! crate's own table over the whole 256-pattern space.
//!
//! The crate does the equivalent work per cell and per corner. Six independent
//! copies of `case |= 1 << c` exist (`marching_cubes/mod.rs:266`,
//! `manifold_dual_contouring.rs:210`, `property/extraction.rs:438` and `:584`,
//! `marching_cubes/ambiguity/tests.rs:260`, `transvoxel/cell.rs:231` as `u16`),
//! and the *cut* question is then answered downstream from the case index. The
//! one place the crate already thinks in words is `dual.rs:359-424`, which packs
//! one bit per **sample** along X and folds four rows with `any & !all` — the
//! same sign planes this circuit consumes, one question earlier.
//!
//! # The circuit, and how its operations are counted
//!
//! `word_ops` is **not estimated**. The circuit is written once, generic over a
//! four-operation alphabet [`Word`], and instantiated twice: on `u64` for the
//! measured arm, and on [`Tally`] — the same `u64` bits beside a `Cell<u32>`
//! counter — for the count. The number in the CSV is what executing that source
//! increments the counter to, and `tally_agrees_with_u64` asserts the two
//! instantiations compute the same bits on all 256 patterns, which is what makes
//! the count a count *of the measured circuit* rather than of a paraphrase.
//!
//! Two figures come out, and both clear the registered bar:
//!
//! - **`word_ops` = 12** — [`cut_planes`], the registered clause's own reading:
//!   *the twelve cut-edge flags for 64 cells derive from the eight sign planes*.
//!   Twelve XORs, one per edge. It is also **optimal**: twelve outputs each
//!   depending on a distinct pair of inputs cannot share a subexpression, so no
//!   circuit over this alphabet does it in eleven.
//! - **`word_ops_total` = 24** — including [`sign_planes`], which builds the
//!   eight planes from the four sample words of the cell rows plus their
//!   successors: the four `dx = 0` planes are the words themselves, and each
//!   `dx = 1` plane is `(lo >> 1) | (hi << 63)`, three operations, four times.
//!   That is `dual.rs:395`'s `inside_word_shifted` verbatim, and it lands the
//!   stricter reading at **exactly** the bar of 24 — not below it.
//!
//! # The comparand is the cheapest byte path, not a strawman
//!
//! The byte arm reads **only the eight signs** — not the eight values. The
//! shipped loop at `marching_cubes/mod.rs:259-268` additionally materialises
//! `[R; 8]` for the interpolation that follows, and `P-121` measured what that
//! gather costs (the shipped `MarchingCubes` runs **15–24% above** a mirror that
//! defers it). Charging it here would inflate the comparand with work that has
//! nothing to do with cut flags. So the byte arm is: gather eight signs, build
//! the case byte, index a 256-entry `u16` table, store the mask — which is the
//! *tightest* honest form of the path and therefore biases C3 towards
//! falsification.
//!
//! The table itself is derived, never transcribed, and then cross-checked
//! against a second shipped object:
//!
//! - **`cut_mask_from_corners`** — `EDGE_CORNERS` and `corner_inside`: bit `j`
//!   set when edge `j`'s two corners differ in sign. This is the byte arm's
//!   table.
//! - **`cut_mask_from_cases`** — the union of non-centroid edge codes appearing
//!   in `table::CASES[case].triangles`, which is a *different* shipped object
//!   built by a *different* `const fn` (`table.rs:182-195`). Every cut edge
//!   belongs to exactly one cycle and every cycle edge appears in the
//!   triangulation, so the two must agree on all 256 patterns.
//!   `table_matches_cases_patterns` is that count, asserted at 256 before any
//!   measurement runs.
//!
//! `reference.rs:29`'s `BOURKE_EDGE_TABLE` — *"bit `j` is set when Bourke's edge
//! `j` is cut, per case index"*, the independent published listing — is
//! `pub(super)` and so is **unreachable from a bench**. It is not copied here on
//! purpose: a hand-carried 256-entry table is exactly the transcription risk the
//! crate avoided by deriving `CASES`, and `marching_cubes/tests.rs:98` already
//! checks the crate's topology against Bourke's under a *derived* cube symmetry.
//! The derivation above inherits that check rather than duplicating it, and
//! `table_matches_cases_patterns` is the local guard that the derivation is
//! sound.
//!
//! # The twelve flags are strictly less information than the case index
//!
//! Complementing all eight signs preserves "these two differ", so the cut mask
//! is invariant under it: `cut_mask(case) == cut_mask(case ^ 0xFF)` for all 256
//! cases. That is the classic Marching Cubes complement ambiguity, and it is
//! recorded as two columns rather than left as a remark —
//! `complement_pairs_agree` = 256 and `distinct_cut_masks`, the number of
//! distinct `u16` values the 256 cases produce.
//!
//! **The consequence bounds what this mechanism can ever be used for.** A
//! bit-parallel cut mask cannot feed `CASES`, because the map from case to mask
//! is not injective. It answers *which grid edges need a crossing vertex* and
//! *is this cell active at all*; the triangle table still needs the case index,
//! which is `P-103`'s row and not this one. So C3's win is a win on the
//! crossing-allocation half of classification, and any reader who reads it as
//! "extraction gets that much faster" is reading it wrong — which is what the
//! SHARE columns below exist to prevent.
//!
//! # Two layouts, and neither arm pays for the other's
//!
//! The two arms answer the same question in different shapes, and there is no
//! layout-neutral consumer: the byte path's natural output is one `u16` per
//! cell, the circuit's is twelve words per 64 cells. A harness that forced
//! either arm into the other's layout would be measuring a transpose. So all
//! four combinations are on the row, and the registered pair is the one the
//! registration names — *the twelve cut-edge flags for 64 cells*:
//!
//! | producing | per-cell `u16` masks | twelve `u64` planes |
//! |---|---|---|
//! | **byte table** | `instructions_per_cell_table` (registered) | `instructions_per_cell_table_planes` |
//! | **SWAR** | `instructions_per_cell_swar_masks` | `instructions_per_cell_swar` (registered) |
//!
//! `ratio` (registered) is `swar / table`, **each arm in its own native
//! layout**; under 1 is a win. `ratio_planes_layout` and `ratio_masks_layout`
//! compare the two arms *within* one layout, and they are the columns that say
//! how much of the verdict is the circuit and how much is the shape of the
//! answer.
//!
//! One caveat stated rather than buried: the cross-layout arms convert with a
//! straightforward per-cell scatter/gather (twelve shift-mask-shift-or per
//! cell). A bit-matrix transpose (Hacker's Delight 7-3) would be cheaper, so
//! `instructions_per_cell_table_planes` is an **upper bound** on the byte path's
//! cost in the plane layout. It is not used here because a bit-matrix transpose
//! is itself a SWAR mechanism, and putting one inside the arm that is supposed to
//! be the byte-table comparand would make the comparand partly the thing under
//! test.
//!
//! # What is inside each window, and what is not
//!
//! - **Field sampling is outside every window.** The values array is filled once
//!   per row and every arm reads the same one. `P-121` measured field evaluation
//!   at up to **94% of extraction on fbm_terrain**; counting it would drive every
//!   ratio to 1.00 and C3 would hold or fail by dilution.
//! - **The sign packing is inside the SWAR arm.** `instructions_per_cell_swar`
//!   is `pack_signs` (`dual.rs:359-381` mirrored: one bit per **sample**, 64 to
//!   a `u64`, along X only, `bit_row = n.div_ceil(64)`) **plus** the circuit,
//!   because the byte arm starts from the same `values` array and a free bitmap
//!   would be an input one arm was handed. `instructions_per_cell_swar_given_bitmap`
//!   is beside it for the regime `dual.rs` is already in, where the bitmap
//!   exists for other reasons, and `instructions_per_cell_pack` is the third so
//!   the decomposition can be checked: `pack_residual_share` is
//!   `(swar − pack − circuit) / swar`, a **prefix difference over sibling
//!   windows**, and it is reported per row rather than asserted away.
//! - **Windows are siblings, never nested.** Zen 3 has six general-purpose
//!   counters and `Probe` opens exactly six, so two nested windows multiplex and
//!   `Counts::worst_ratio` refuses — `R-121` paid for that discovery. Eight
//!   sibling windows per repetition, and the two **registered** arms alternate
//!   which runs first by repetition parity so neither is permanently the one
//!   that inherits the other's cache state.
//! - **The denominator is the shipped extractor**, `MarchingCubes::extract` on
//!   the same grid in the same run: `instructions_extract_mc` and
//!   `cycles_extract_mc`. Nothing is mirrored into it, so there is no
//!   mirror-agreement debt on the denominator.
//!
//! # The popcount this build does not have
//!
//! There is no `.cargo/config.toml` and no `target-cpu` anywhere in the
//! repository, so the default `x86-64` baseline is in force,
//! `cfg!(target_feature = "popcnt")` is **false**, and `u64::count_ones()`
//! lowers to a ~12-instruction SWAR sequence rather than to the one instruction
//! this CPU has. Several sibling rows are priced against that.
//!
//! **This row is not, on either side.** Neither arm makes a single
//! `count_ones` call: the circuit is XOR, shift and OR, the packer is a
//! shift-and-OR fold, and the byte arm is a table index.
//! `count_ones_calls_per_cell_table` and `count_ones_calls_per_cell_swar` are
//! both **0** and are columns, `target_feature_popcnt` is a column, and **no
//! verdict here is contingent on the popcount lowering** — a
//! `-C target-cpu=native` build could not move `ratio`, and measuring one would
//! mean comparing across binaries, which `M-281` forbids. Checked on this
//! bench's own binary rather than assumed: `objdump -d
//! target/release/deps/experiment_p106-<hash>` greps **0** `popcnt` and **0**
//! occurrences of the SWAR magic constant `0x3333333333333333`, where the
//! sibling that found the fact reports 91 of that constant in its own binary —
//! so the grep does discriminate.
//!
//! # SHARE
//!
//! Every clause's reachable share is a column, and the denominator is measured
//! here rather than inherited.
//!
//! `P-121`'s gate **opened** for this row: `integer_share` — classify plus emit
//! over total — is **0.7260** on `sphere 65³ f32 marching_cubes` against the
//! 0.15 bar, and `classify` **alone** is **50–65%** of `marching_cubes`
//! extraction on the compact fields (0.6459 on `sphere 65³ f32`). Those are
//! **priors, cited, not imported**: they collapse to **1.8–6.7%** on
//! `fbm_terrain` and `noise_cavity`, where field evaluation is up to 94% of
//! extraction, and `P-121`'s `integer_share` on `fbm_terrain` is **0.0215**. A
//! mechanism that halves classification therefore moves extraction by a third on
//! one field and by nothing on another, and quoting one field's share as the
//! row's share is how `✗51` happened.
//!
//! So this row measures its own:
//!
//! - **`case_stream_share_instructions`** — the case-index build alone,
//!   materialised as one byte per cell, over `instructions_extract_mc`. This is
//!   the quantity comparable with `P-121`'s `classify`, measured in this binary,
//!   and it is on every row so the 50–65% prior can be checked rather than
//!   trusted.
//! - **`cut_edge_share_instructions`** — the byte arm (case index plus the cut
//!   lookup and store) over the same denominator. Stated plainly: **the shipped
//!   extractor never materialises a cut mask**, so this is what a cut-mask pass
//!   *would* cost as a fraction of a whole extraction, not a stage share of the
//!   shipped path.
//! - **`swar_saving_share_instructions`** = `(table − swar) / extract`, and
//!   **`extraction_ceiling`** = `1 / (1 − swar_saving_share_instructions)` — the
//!   `✗51` arithmetic, per row. This is C3's reachable share: the ceiling on
//!   what the mechanism could do to an extraction even at zero cost, and on the
//!   sparse fields it is close to 1.00× by construction.
//! - **C1's share is not a share at all.** It is an exact integer derived from
//!   the source (`word_ops`, `word_ops_total`) against an exact bar
//!   (`word_op_bar` = 24), identical in every build and on every machine.
//! - **C2's share is an enumerated population**: `patterns_tested` = 256, the
//!   whole space, so its denominator is exact by construction. `cells_agreeing`
//!   over `cells` extends the same equality to every cell of every field, and
//!   `cells_agreeing_random` repeats it on a 50%-dense pseudorandom sign field
//!   where every word boundary is crossed.
//!
//! **VACUITY CONTROL, asserted rather than recorded.** `patterns_tested` must
//! equal **256** exactly — `✗50` is the incident that earns the word
//! *exhaustively*, where a **sampled** bound became a release-build panic — and
//! the patterns are tested **packed 64 to a word**, each in a different bit
//! lane, so a circuit that is right about the logic and wrong about a lane
//! cannot pass. `mutant_pattern_mismatches` must exceed zero: the same 256-pattern
//! comparison run against a **deliberately broken circuit** (edge 0 computed as
//! `lo & hi`, "both corners inside", instead of `lo ^ hi`, "the corners differ"),
//! which is the shape of the bug a reader would most plausibly write. Three more
//! guards sit beside it, all asserted:
//! `mutant_cell_mismatches_carry` > 0 — the plane build with `hi << 63` dropped,
//! which is `dual.rs:395`'s documented trap and would punch a hole every 64
//! cells; `tally_agrees_with_u64`, without which `word_ops` counts a different
//! circuit from the measured one; and `cut_cells` > 0, because a row where no
//! edge is cut compares two ways of computing zero.
//!
//! The carry mutant is scored on the **pseudorandom** field, deliberately.
//! On every reference field the domain boundary is outside the solid, so the
//! sample at `x = n − 1` is outside, its sign bit is 0, and dropping the carry-in
//! changes **nothing** — the control would read 0 and prove the comparator dead
//! when it is not. `mutant_cell_mismatches_carry` is therefore taken where the
//! defect is visible, and that is the whole reason a random arm exists.
//!
//! # `ns_per_cell` and `ghz` are provenance, not verdicts
//!
//! `M-280` and `M-281`, and `R-105` watched one binary's cycle ratio band move
//! from 0.984 to 1.035 across three runs while its instruction counts held to
//! four figures. C3 reads instructions; `instruction_ratio_rep_spread`
//! *demonstrates* their determinism instead of asserting it, and
//! `cycle_ratio` and the two `ns_per_cell` columns are reported beside `ghz` so
//! a later reader can see what clock they were taken at. No clause consults
//! either.

#![allow(
    clippy::cast_possible_truncation,
    clippy::cast_precision_loss,
    clippy::too_many_lines
)]

mod common;

#[cfg(target_os = "linux")]
mod experiment {
    use std::cell::Cell;
    use std::hint::black_box;
    use std::time::Instant;

    use isomesh::marching_cubes::MarchingCubes;
    use isomesh::marching_cubes::table::{
        CASES, CENTROID_BASE, EDGE_CORNERS, EDGE_COUNT, corner_inside, is_inside,
    };
    use isomesh::{MeshBuffer, RuntimeShape3, Sdf};

    use crate::common::counters::{MIN_TIME_RATIO, Probe};

    // ─── the fixture ───────────────────────────────────────────────────────

    /// The registered fixture's resolution, in samples per axis: eight
    /// reference fields × 65³.
    ///
    /// 65 samples span 64 cells, so `cell_words` is exactly 1 and there are **no
    /// pad lanes** — while the *sample* row needs two words. That asymmetry is
    /// `dual.rs`'s own (`bit_row = n.div_ceil(64)` against `cell_words =
    /// cells.div_ceil(64)`) and it is what makes the `hi << 63` carry-in load
    /// bearing at this exact resolution: word 1 holds one single sample, `x =
    /// 64`, and it is the `+x` corner of every cell at `x = 63`.
    const RESOLUTION: u32 = 65;

    /// Measured repetitions per window per row, and **odd on purpose**: the
    /// median of an odd count is a reading rather than a mean of two.
    ///
    /// Seven rather than more because the quantity that carries C3 is an
    /// instruction count, which `R-105` established is deterministic to four
    /// figures on this machine; `instruction_ratio_rep_spread` reports the band
    /// so the claim is visible rather than assumed.
    const REPS: usize = 7;

    /// Untimed passes of every arm before anything is counted, so buffers are at
    /// final capacity and the pages are faulted in.
    const WARMUP: usize = 2;

    /// About this long per counter window, so the ~28 `perf_event` system calls
    /// a window costs land outside it and cannot inflate anything.
    const TARGET_BATCH_NS: f64 = 20_000_000.0;

    /// Ceiling on the batch, so a cheap pass cannot run for a minute.
    const MAX_INNER: usize = 4096;

    /// C1's bar, from the registration: at most 24 word operations.
    ///
    /// 24 is where the byte-table path wins on instruction count alone, which is
    /// why the clause is an integer bound and not a ratio.
    const WORD_OP_BAR: u32 = 24;

    /// The whole sign-pattern space. C2 is scored over all of it.
    const PATTERNS: usize = 256;

    /// Patterns per word when the sweep is packed one pattern per bit lane.
    const LANES: usize = 64;

    // ─── the circuit's alphabet, in two instantiations ──────────────────────

    /// The four word operations the circuit is allowed, so that counting them is
    /// a property of the source rather than of a reading of it.
    ///
    /// [`u64`] is the measured instantiation; [`Tally`] is the counting one.
    /// `and` exists for the mutant alone — the real circuit never uses it, which
    /// is exactly why substituting it is a defect the comparator must catch.
    trait Word: Copy {
        /// "These two corners differ in sign", per lane.
        fn xor(self, other: Self) -> Self;
        /// The mutant's substitution: "both corners are inside".
        fn and(self, other: Self) -> Self;
        /// Merge two half-planes.
        fn or(self, other: Self) -> Self;
        /// Shift lanes down by `bits`.
        fn shr(self, bits: u32) -> Self;
        /// Shift lanes up by `bits`.
        fn shl(self, bits: u32) -> Self;
    }

    impl Word for u64 {
        #[inline]
        fn xor(self, other: Self) -> Self {
            self ^ other
        }
        #[inline]
        fn and(self, other: Self) -> Self {
            self & other
        }
        #[inline]
        fn or(self, other: Self) -> Self {
            self | other
        }
        #[inline]
        fn shr(self, bits: u32) -> Self {
            self >> bits
        }
        #[inline]
        fn shl(self, bits: u32) -> Self {
            self << bits
        }
    }

    /// The same `u64` bits, beside a counter every operation increments.
    ///
    /// This is how `word_ops` is *derived from the code*: the circuit is written
    /// once and executed twice, and the CSV reports what executing it costs in
    /// this alphabet. `tally_agrees_with_u64` is what says the two
    /// instantiations compute the same thing.
    #[derive(Clone, Copy)]
    struct Tally<'c> {
        bits: u64,
        ops: &'c Cell<u32>,
    }

    impl<'c> Tally<'c> {
        #[inline]
        fn step(self, bits: u64) -> Self {
            self.ops.set(self.ops.get() + 1);
            Self {
                bits,
                ops: self.ops,
            }
        }
    }

    impl Word for Tally<'_> {
        #[inline]
        fn xor(self, other: Self) -> Self {
            self.step(self.bits ^ other.bits)
        }
        #[inline]
        fn and(self, other: Self) -> Self {
            self.step(self.bits & other.bits)
        }
        #[inline]
        fn or(self, other: Self) -> Self {
            self.step(self.bits | other.bits)
        }
        #[inline]
        fn shr(self, bits: u32) -> Self {
            self.step(self.bits >> bits)
        }
        #[inline]
        fn shl(self, bits: u32) -> Self {
            self.step(self.bits << bits)
        }
    }

    /// The eight sign planes of 64 cells, from four sample words and their
    /// successors.
    ///
    /// `lo[i]` is `dual.rs:385`'s `inside_word` for the cell row's `(dy, dz)`
    /// offset `i = dy + 2·dz`; `hi[i]` is the next word of the same sample row,
    /// or zero past the end. Corner `c = (dx, dy, dz)` is `2·i + dx`, so the
    /// `dx = 0` planes are the words themselves and the `dx = 1` planes are
    /// `dual.rs:395`'s `inside_word_shifted`.
    ///
    /// **Twelve word operations**: three per shifted plane, four times.
    #[inline]
    fn sign_planes<W: Word>(lo: [W; 4], hi: [W; 4]) -> [W; 8] {
        let mut planes = [lo[0]; 8];
        for i in 0..4 {
            planes[2 * i] = lo[i];
            planes[2 * i + 1] = lo[i].shr(1).or(hi[i].shl(63));
        }
        planes
    }

    /// The same thing with the carry-in dropped — the plane-build mutant.
    ///
    /// `dual.rs:395` documents the defect this is: without the high bit from the
    /// next word, the cell straddling a word boundary reads its `+x` corner as
    /// outside, which is a hole every 64 cells.
    #[inline]
    fn sign_planes_no_carry<W: Word>(lo: [W; 4], _hi: [W; 4]) -> [W; 8] {
        let mut planes = [lo[0]; 8];
        for i in 0..4 {
            planes[2 * i] = lo[i];
            planes[2 * i + 1] = lo[i].shr(1);
        }
        planes
    }

    /// **The circuit.** Twelve cut-edge planes from eight sign planes.
    ///
    /// An edge is cut exactly when its two corners differ in sign, so each
    /// output is one XOR of the pair `cube.rs:26`'s `EDGE_CORNERS` names.
    ///
    /// **Twelve word operations**, and no fewer is possible over this alphabet:
    /// the twelve edges are twelve distinct corner pairs, so no two outputs share
    /// a subexpression.
    #[inline]
    fn cut_planes<W: Word>(planes: &[W; 8]) -> [W; EDGE_COUNT] {
        let mut cuts = [planes[0]; EDGE_COUNT];
        for (edge, slot) in cuts.iter_mut().enumerate() {
            let [lo, hi] = EDGE_CORNERS[edge];
            *slot = planes[lo as usize].xor(planes[hi as usize]);
        }
        cuts
    }

    /// The registered vacuity control: the circuit with edge 0 computed as
    /// `lo & hi` — "both corners are inside" — instead of `lo ^ hi`.
    ///
    /// The two agree only when both corners are outside, so this must disagree
    /// with the table on three quarters of the 256 patterns. A comparison that
    /// cannot fail is not evidence, and `mutant_pattern_mismatches` is what says
    /// this one can.
    #[inline]
    fn cut_planes_mutant<W: Word>(planes: &[W; 8]) -> [W; EDGE_COUNT] {
        let mut cuts = cut_planes(planes);
        let [lo, hi] = EDGE_CORNERS[0];
        cuts[0] = planes[lo as usize].and(planes[hi as usize]);
        cuts
    }

    /// The twelve flags of one lane, gathered back into a `u16` mask.
    #[inline]
    fn mask_at(cuts: &[u64; EDGE_COUNT], lane: u32) -> u16 {
        let mut mask = 0u16;
        for (edge, word) in cuts.iter().enumerate() {
            mask |= (((word >> lane) & 1) as u16) << edge;
        }
        mask
    }

    /// Word operations in the plane build and in the circuit, by execution.
    fn count_word_ops() -> (u32, u32) {
        let ops = Cell::new(0u32);
        let zero = Tally { bits: 0, ops: &ops };
        let planes = sign_planes([zero; 4], [zero; 4]);
        let build = ops.get();
        let cuts = cut_planes(&planes);
        let circuit = ops.get() - build;
        // Nothing is elided: `Cell::get` reads a value the operations wrote.
        black_box(cuts[EDGE_COUNT - 1].bits);
        (build, circuit)
    }

    // ─── the byte table, derived twice from shipped data ────────────────────

    /// The byte arm's table: bit `j` set when edge `j`'s two corners differ in
    /// sign, from `EDGE_CORNERS` and `corner_inside`.
    fn cut_mask_from_corners(case: u8) -> u16 {
        let mut mask = 0u16;
        for (edge, [lo, hi]) in EDGE_CORNERS.iter().copied().enumerate() {
            if corner_inside(case, lo) != corner_inside(case, hi) {
                mask |= 1 << edge;
            }
        }
        mask
    }

    /// The same twelve flags read off `table::CASES` — a different shipped
    /// object, built by a different `const fn`.
    ///
    /// Every cut edge belongs to exactly one cycle and every cycle edge appears
    /// in the triangulation, so the union of a case's non-centroid triangle
    /// codes is its cut-edge set. Codes at or above `CENTROID_BASE` are cycle
    /// centroids, not cube edges.
    fn cut_mask_from_cases(case: u8) -> u16 {
        let entry = &CASES[case as usize];
        let mut mask = 0u16;
        for triangle in &entry.triangles[..entry.count as usize] {
            for &code in triangle {
                if code < CENTROID_BASE {
                    mask |= 1 << code;
                }
            }
        }
        mask
    }

    /// The 256-entry table the byte arm indexes.
    fn edge_cut_table() -> [u16; PATTERNS] {
        let mut table = [0u16; PATTERNS];
        for (case, slot) in table.iter_mut().enumerate() {
            *slot = cut_mask_from_corners(case as u8);
        }
        table
    }

    // ─── the geometry, mirrored from `dual.rs` ──────────────────────────────

    /// One cubic grid's word layout, `dual.rs`'s own.
    #[derive(Clone, Copy)]
    struct Geom {
        /// Samples per axis.
        n: usize,
        /// Sample words per row, `n.div_ceil(64)` — `dual.rs:363`.
        bit_row: usize,
        /// Cells per axis, `n − 1`.
        cells: usize,
        /// Cell words per row, `cells.div_ceil(64)` — `dual.rs:484`.
        cell_words: usize,
    }

    impl Geom {
        fn new(n: u32) -> Self {
            let n = n as usize;
            let cells = n - 1;
            Self {
                n,
                bit_row: n.div_ceil(64),
                cells,
                cell_words: cells.div_ceil(64),
            }
        }

        /// Samples in the grid.
        #[inline]
        fn samples(self) -> usize {
            self.n * self.n * self.n
        }

        /// Cells in the grid.
        #[inline]
        fn cell_count(self) -> usize {
            self.cells * self.cells * self.cells
        }

        /// Sixty-four-cell groups in the grid.
        #[inline]
        fn groups(self) -> usize {
            self.cells * self.cells * self.cell_words
        }

        /// Lanes of a group that are cells rather than padding.
        ///
        /// At the registered 65³ this is always 64 and `pad_lanes` is 0, so the
        /// circuit needs no `cell_mask` (`dual.rs:445`) and none is charged to
        /// it.
        #[inline]
        fn lanes(self, w: usize) -> usize {
            (self.cells - w * 64).min(64)
        }

        /// Bit `k` of this word is sample `64w + k` of row `(y, z)`.
        /// `dual.rs:385`.
        #[inline]
        fn inside_word(self, inside: &[u64], w: usize, y: usize, z: usize) -> u64 {
            inside[self.bit_row * (y + self.n * z) + w]
        }

        /// The four cell-row words of a group, and their successors.
        #[inline]
        fn fetch(self, inside: &[u64], w: usize, y: usize, z: usize) -> ([u64; 4], [u64; 4]) {
            let mut lo = [0u64; 4];
            let mut hi = [0u64; 4];
            for dz in 0..2 {
                for dy in 0..2 {
                    let i = dy + 2 * dz;
                    lo[i] = self.inside_word(inside, w, y + dy, z + dz);
                    hi[i] = if w + 1 < self.bit_row {
                        self.inside_word(inside, w + 1, y + dy, z + dz)
                    } else {
                        0
                    };
                }
            }
            (lo, hi)
        }

        /// The eight-sign case index of one cell, `mod.rs:259-268`'s shape with
        /// the `[R; 8]` value gather left out — see the module docs.
        #[inline]
        fn case_at(self, values: &[f32], x: usize, y: usize, z: usize) -> u8 {
            let mut case = 0u8;
            for c in 0..8u8 {
                let dx = (c & 1) as usize;
                let dy = ((c >> 1) & 1) as usize;
                let dz = ((c >> 2) & 1) as usize;
                let v = values[(x + dx) + self.n * ((y + dy) + self.n * (z + dz))];
                if is_inside(v) {
                    case |= 1 << c;
                }
            }
            case
        }
    }

    // ─── the arms ──────────────────────────────────────────────────────────

    /// `dual.rs:359-381`'s `build_inside_bits`: one bit per **sample**, 64 to a
    /// `u64`, packed along X only.
    fn pack_signs(values: &[f32], geo: Geom, inside: &mut [u64]) {
        black_box(values.as_ptr());
        for row in 0..geo.n * geo.n {
            let src = geo.n * row;
            let dst = geo.bit_row * row;
            for w in 0..geo.bit_row {
                let base = w * 64;
                let bits = (geo.n - base).min(64);
                let mut word = 0u64;
                for k in 0..bits {
                    word |= u64::from(is_inside(values[src + base + k])) << k;
                }
                inside[dst + w] = word;
            }
        }
        black_box(&*inside);
    }

    /// **The SWAR arm's second half**: the twelve cut-edge planes of every
    /// group, from the sign bitmap.
    ///
    /// `CARRY` false is the plane-build mutant, `MUTANT` true is the circuit
    /// mutant. Const generics rather than closures so the measured
    /// instantiation, `<true, false>`, has no indirection in it.
    fn circuit_pass<const CARRY: bool, const MUTANT: bool>(
        inside: &[u64],
        geo: Geom,
        planes: &mut [u64],
    ) {
        black_box(inside.as_ptr());
        for z in 0..geo.cells {
            for y in 0..geo.cells {
                for w in 0..geo.cell_words {
                    let (lo, hi) = geo.fetch(inside, w, y, z);
                    let signs = if CARRY {
                        sign_planes(lo, hi)
                    } else {
                        sign_planes_no_carry(lo, hi)
                    };
                    let cuts = if MUTANT {
                        cut_planes_mutant(&signs)
                    } else {
                        cut_planes(&signs)
                    };
                    let at = ((z * geo.cells + y) * geo.cell_words + w) * EDGE_COUNT;
                    planes[at..at + EDGE_COUNT].copy_from_slice(&cuts);
                }
            }
        }
        black_box(&*planes);
    }

    /// **The byte arm.** Eight signs, a case byte, a table index, a `u16` store.
    fn table_pass(values: &[f32], geo: Geom, table: &[u16; PATTERNS], masks: &mut [u16]) {
        black_box(values.as_ptr());
        for z in 0..geo.cells {
            for y in 0..geo.cells {
                let row = (z * geo.cells + y) * geo.cells;
                for x in 0..geo.cells {
                    masks[row + x] = table[geo.case_at(values, x, y, z) as usize];
                }
            }
        }
        black_box(&*masks);
    }

    /// The case index alone, materialised one byte per cell.
    ///
    /// `P-121`'s `classify` in this binary, and the SHARE section's own
    /// denominator check.
    fn case_stream_pass(values: &[f32], geo: Geom, cases: &mut [u8]) {
        black_box(values.as_ptr());
        for z in 0..geo.cells {
            for y in 0..geo.cells {
                let row = (z * geo.cells + y) * geo.cells;
                for x in 0..geo.cells {
                    cases[row + x] = geo.case_at(values, x, y, z);
                }
            }
        }
        black_box(&*cases);
    }

    /// The byte arm delivering the **plane** layout: twelve words per group,
    /// scattered a cell at a time.
    ///
    /// The upper bound discussed in the module docs — a bit-matrix transpose
    /// would be cheaper, and is not used here because it is itself the mechanism
    /// under test.
    fn table_planes_pass(values: &[f32], geo: Geom, table: &[u16; PATTERNS], planes: &mut [u64]) {
        black_box(values.as_ptr());
        for z in 0..geo.cells {
            for y in 0..geo.cells {
                for w in 0..geo.cell_words {
                    let mut acc = [0u64; EDGE_COUNT];
                    for k in 0..geo.lanes(w) {
                        let mask = u64::from(table[geo.case_at(values, w * 64 + k, y, z) as usize]);
                        for (edge, word) in acc.iter_mut().enumerate() {
                            *word |= ((mask >> edge) & 1) << k;
                        }
                    }
                    let at = ((z * geo.cells + y) * geo.cell_words + w) * EDGE_COUNT;
                    planes[at..at + EDGE_COUNT].copy_from_slice(&acc);
                }
            }
        }
        black_box(&*planes);
    }

    /// The SWAR arm delivering the **per-cell mask** layout: the circuit, then
    /// one `u16` per lane.
    fn swar_masks_pass(inside: &[u64], geo: Geom, masks: &mut [u16]) {
        black_box(inside.as_ptr());
        for z in 0..geo.cells {
            for y in 0..geo.cells {
                for w in 0..geo.cell_words {
                    let (lo, hi) = geo.fetch(inside, w, y, z);
                    let cuts = cut_planes(&sign_planes(lo, hi));
                    let row = (z * geo.cells + y) * geo.cells + w * 64;
                    for k in 0..geo.lanes(w) {
                        masks[row + k] = mask_at(&cuts, k as u32);
                    }
                }
            }
        }
        black_box(&*masks);
    }

    // ─── counting ──────────────────────────────────────────────────────────

    /// What one counter window read, per repetition of the body.
    #[derive(Clone, Copy, Default)]
    struct Counted {
        cycles: f64,
        instructions: f64,
        nanos: f64,
    }

    /// One counter window over `inner` repetitions, divided by `inner`.
    ///
    /// Every `perf_event` system call is outside the counted region. Windows are
    /// **siblings, never nested**: `Probe` opens exactly the six general-purpose
    /// counters Zen 3 has, so a nested window multiplexes and the ratio assert
    /// below fires.
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

    /// The median of a set of readings.
    fn median(values: &mut [f64]) -> f64 {
        values.sort_by(f64::total_cmp);
        values[values.len() / 2]
    }

    /// The eight sibling windows of one repetition.
    #[derive(Clone, Copy, Default)]
    struct Rep {
        table: Counted,
        swar: Counted,
        pack: Counted,
        circuit: Counted,
        table_planes: Counted,
        swar_masks: Counted,
        case_stream: Counted,
        extract: Counted,
    }

    /// Median cycles, instructions and nanoseconds of one window, per quantity.
    fn median_counted(select: impl Fn(&Rep) -> Counted, reps: &[Rep]) -> Counted {
        let mut cycles: Vec<f64> = reps.iter().map(|r| select(r).cycles).collect();
        let mut instructions: Vec<f64> = reps.iter().map(|r| select(r).instructions).collect();
        let mut nanos: Vec<f64> = reps.iter().map(|r| select(r).nanos).collect();
        Counted {
            cycles: median(&mut cycles),
            instructions: median(&mut instructions),
            nanos: median(&mut nanos),
        }
    }

    // ─── the pattern-level checks, field-independent and run once ───────────

    /// Everything C1 and C2 are scored on, plus the controls that license them.
    #[derive(Clone, Copy)]
    struct Circuit {
        word_ops_plane_build: u32,
        word_ops_circuit: u32,
        patterns_tested: usize,
        /// Patterns whose twelve flags match the table, tested **packed**, one
        /// pattern per bit lane.
        identical: usize,
        /// The same over uniform planes — the weaker test, kept so a lane bug and
        /// a logic bug are distinguishable.
        identical_uniform: usize,
        mutant_mismatches: usize,
        table_matches_cases: usize,
        distinct_cut_masks: usize,
        complement_pairs_agree: usize,
        tally_agrees: bool,
    }

    /// All 256 patterns packed 64 to a word, one pattern per bit lane.
    fn pattern_planes(group: usize) -> [u64; 8] {
        let mut planes = [0u64; 8];
        for lane in 0..LANES {
            let pattern = (group * LANES + lane) as u8;
            for (corner, plane) in planes.iter_mut().enumerate() {
                *plane |= u64::from((pattern >> corner) & 1) << lane;
            }
        }
        planes
    }

    /// The whole sign-pattern space, against the crate's own table.
    fn examine_circuit(table: &[u16; PATTERNS]) -> Circuit {
        let (word_ops_plane_build, word_ops_circuit) = count_word_ops();

        let mut patterns_tested = 0usize;
        let mut identical = 0usize;
        let mut mutant_mismatches = 0usize;
        let mut tally_agrees = true;
        let ops = Cell::new(0u32);
        for group in 0..PATTERNS / LANES {
            let planes = pattern_planes(group);
            let cuts = cut_planes(&planes);
            let mutant = cut_planes_mutant(&planes);

            // The counting instantiation must compute the same bits, or
            // `word_ops` is the operation count of a different circuit.
            let tallied = cut_planes(&planes.map(|bits| Tally { bits, ops: &ops }));
            for (edge, word) in cuts.iter().enumerate() {
                tally_agrees &= tallied[edge].bits == *word;
            }

            for lane in 0..LANES {
                let pattern = group * LANES + lane;
                patterns_tested += 1;
                if mask_at(&cuts, lane as u32) == table[pattern] {
                    identical += 1;
                }
                if mask_at(&mutant, lane as u32) != table[pattern] {
                    mutant_mismatches += 1;
                }
            }
        }

        // The weaker, uniform-plane form of the same sweep: every lane of a
        // group carries the same pattern, so this cannot see a lane bug.
        let mut identical_uniform = 0usize;
        for (pattern, expected) in table.iter().enumerate() {
            let mut planes = [0u64; 8];
            for (corner, plane) in planes.iter_mut().enumerate() {
                *plane = if (pattern >> corner) & 1 == 1 { !0 } else { 0 };
            }
            let cuts = cut_planes(&planes);
            let agrees = cuts
                .iter()
                .enumerate()
                .all(|(edge, word)| *word == if (expected >> edge) & 1 == 1 { !0 } else { 0 });
            if agrees {
                identical_uniform += 1;
            }
        }

        // The reference itself, against a second shipped object.
        let mut table_matches_cases = 0usize;
        let mut complement_pairs_agree = 0usize;
        let mut seen: Vec<u16> = Vec::with_capacity(PATTERNS);
        for (case, mask) in table.iter().enumerate() {
            if cut_mask_from_cases(case as u8) == *mask {
                table_matches_cases += 1;
            }
            if table[case ^ 0xFF] == *mask {
                complement_pairs_agree += 1;
            }
            if !seen.contains(mask) {
                seen.push(*mask);
            }
        }

        Circuit {
            word_ops_plane_build,
            word_ops_circuit,
            patterns_tested,
            identical,
            identical_uniform,
            mutant_mismatches,
            table_matches_cases,
            distinct_cut_masks: seen.len(),
            complement_pairs_agree,
            tally_agrees,
        }
    }

    // ─── the per-cell cross-check, and its two fixtures ─────────────────────

    /// How many cells the two arms agree about, over a whole grid.
    fn agreeing_cells(values: &[f32], geo: Geom, table: &[u16; PATTERNS], carry: bool) -> usize {
        let mut inside = vec![0u64; geo.bit_row * geo.n * geo.n];
        let mut masks = vec![0u16; geo.cell_count()];
        let mut planes = vec![0u64; geo.groups() * EDGE_COUNT];
        pack_signs(values, geo, &mut inside);
        table_pass(values, geo, table, &mut masks);
        if carry {
            circuit_pass::<true, false>(&inside, geo, &mut planes);
        } else {
            circuit_pass::<false, false>(&inside, geo, &mut planes);
        }

        let mut agreeing = 0usize;
        for z in 0..geo.cells {
            for y in 0..geo.cells {
                for w in 0..geo.cell_words {
                    let at = ((z * geo.cells + y) * geo.cell_words + w) * EDGE_COUNT;
                    let cuts: [u64; EDGE_COUNT] = planes[at..at + EDGE_COUNT]
                        .try_into()
                        .expect("twelve words per group");
                    let row = (z * geo.cells + y) * geo.cells + w * 64;
                    for k in 0..geo.lanes(w) {
                        if mask_at(&cuts, k as u32) == masks[row + k] {
                            agreeing += 1;
                        }
                    }
                }
            }
        }
        agreeing
    }

    /// A deterministic pseudorandom sign field, ±1.0 per sample.
    ///
    /// About half the samples are inside, so every word boundary is crossed and
    /// the carry-in mutant is visible — which it is not on any reference field,
    /// where the domain boundary is outside the solid.
    fn random_sign_values(samples: usize) -> Vec<f32> {
        let mut state = 0x243F_6A88_85A3_08D3u64;
        (0..samples)
            .map(|_| {
                state = state.wrapping_add(0x9E37_79B9_7F4A_7C15);
                let mut z = state;
                z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
                z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
                z ^= z >> 31;
                if z & 1 == 0 { 1.0 } else { -1.0 }
            })
            .collect()
    }

    /// The two random-field controls, field-independent and run once.
    #[derive(Clone, Copy)]
    struct Controls {
        cells: usize,
        agreeing: usize,
        carry_mutant_mismatches: usize,
    }

    fn examine_controls(geo: Geom, table: &[u16; PATTERNS]) -> Controls {
        let values = random_sign_values(geo.samples());
        let cells = geo.cell_count();
        Controls {
            cells,
            agreeing: agreeing_cells(&values, geo, table, true),
            carry_mutant_mismatches: cells - agreeing_cells(&values, geo, table, false),
        }
    }

    // ─── one row ───────────────────────────────────────────────────────────

    struct Row {
        field: &'static str,
        resolution: u32,
        cells: usize,
        samples: usize,
        cut_cells: usize,
        cut_fraction: f64,
        cells_agreeing: usize,
        pad_lanes: usize,
        cell_words: usize,
        sample_bit_row: usize,
        output_bytes_table: usize,
        output_bytes_swar: usize,
        table: Counted,
        swar: Counted,
        pack: Counted,
        circuit: Counted,
        table_planes: Counted,
        swar_masks: Counted,
        case_stream: Counted,
        extract: Counted,
        instruction_ratio_rep_spread: f64,
        inner_table: usize,
        inner_swar: usize,
        inner_extract: usize,
        vertices: usize,
    }

    impl Row {
        /// C3's registered reading: each arm in its own native layout. Under 1
        /// is a win.
        fn ratio(&self) -> f64 {
            self.swar.instructions / self.table.instructions
        }

        /// Both arms delivering the twelve planes.
        fn ratio_planes_layout(&self) -> f64 {
            self.swar.instructions / self.table_planes.instructions
        }

        /// Both arms delivering one `u16` per cell.
        fn ratio_masks_layout(&self) -> f64 {
            self.swar_masks.instructions / self.table.instructions
        }

        fn cycle_ratio(&self) -> f64 {
            self.swar.cycles / self.table.cycles
        }

        /// `(swar − pack − circuit) / swar`: how much of the SWAR arm the two
        /// sibling windows fail to account for.
        fn pack_residual_share(&self) -> f64 {
            (self.swar.instructions - self.pack.instructions - self.circuit.instructions)
                / self.swar.instructions
        }

        fn case_stream_share(&self) -> f64 {
            self.case_stream.instructions / self.extract.instructions
        }

        fn cut_edge_share(&self) -> f64 {
            self.table.instructions / self.extract.instructions
        }

        fn swar_saving_share(&self) -> f64 {
            (self.table.instructions - self.swar.instructions) / self.extract.instructions
        }

        /// `✗51`'s arithmetic: the ceiling on extraction if the saving were free.
        fn extraction_ceiling(&self) -> f64 {
            1.0 / (1.0 - self.swar_saving_share())
        }

        fn ghz(&self) -> f64 {
            (self.table.cycles + self.swar.cycles) / (self.table.nanos + self.swar.nanos)
        }
    }

    fn measure<S: Sdf<Scalar = f32>>(
        field: &'static str,
        geo: Geom,
        table: &[u16; PATTERNS],
        sdf: &S,
        origin: [f32; 3],
        cell_size: f32,
    ) -> Row {
        let n = geo.n as u32;
        let shape = RuntimeShape3::new([n; 3]).expect("the fixture fits u32");

        // ── the values, sampled the way `MarchingCubes::extract` samples ─────
        let mut values = Vec::with_capacity(geo.samples());
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

        let mut inside = vec![0u64; geo.bit_row * geo.n * geo.n];
        let mut masks = vec![0u16; geo.cell_count()];
        let mut cases = vec![0u8; geo.cell_count()];
        let mut planes = vec![0u64; geo.groups() * EDGE_COUNT];

        // ── the population, and C2 extended to every cell of this field ──────
        let cells_agreeing = agreeing_cells(&values, geo, table, true);
        table_pass(&values, geo, table, &mut masks);
        let cut_cells = masks.iter().filter(|m| **m != 0).count();
        assert!(
            cut_cells > 0,
            "{field} {n}^3: no cell has a cut edge, so both arms are computing zero and the \
             comparison is between two ways of doing nothing"
        );

        // ── the shipped extractor, the denominator ───────────────────────────
        let mut shipped = MarchingCubes::<f32>::new();
        let mut out = MeshBuffer::<f32>::new();
        shipped
            .extract(sdf, &shape, origin, cell_size, &mut out)
            .expect("extraction");
        let vertices = out.positions.len();

        for _ in 0..WARMUP {
            table_pass(&values, geo, table, &mut masks);
            pack_signs(&values, geo, &mut inside);
            circuit_pass::<true, false>(&inside, geo, &mut planes);
            table_planes_pass(&values, geo, table, &mut planes);
            swar_masks_pass(&inside, geo, &mut masks);
            case_stream_pass(&values, geo, &mut cases);
            shipped
                .extract(sdf, &shape, origin, cell_size, &mut out)
                .expect("extraction");
        }

        let inner_table = choose_inner(|| table_pass(&values, geo, table, &mut masks));
        let inner_swar = choose_inner(|| {
            pack_signs(&values, geo, &mut inside);
            circuit_pass::<true, false>(&inside, geo, &mut planes);
        });
        let inner_pack = choose_inner(|| pack_signs(&values, geo, &mut inside));
        let inner_circuit = choose_inner(|| circuit_pass::<true, false>(&inside, geo, &mut planes));
        let inner_table_planes =
            choose_inner(|| table_planes_pass(&values, geo, table, &mut planes));
        let inner_swar_masks = choose_inner(|| {
            pack_signs(&values, geo, &mut inside);
            swar_masks_pass(&inside, geo, &mut masks);
        });
        let inner_case_stream = choose_inner(|| case_stream_pass(&values, geo, &mut cases));
        let inner_extract = choose_inner(|| {
            shipped
                .extract(sdf, &shape, origin, cell_size, &mut out)
                .expect("extraction");
        });

        let mut probe = Probe::open();
        let mut reps: Vec<Rep> = Vec::with_capacity(REPS);
        for rep in 0..REPS {
            let mut r = Rep::default();
            // The two registered arms alternate which runs first, so neither is
            // permanently the one that inherits the other's cache state.
            if rep % 2 == 0 {
                r.table = window(&mut probe, inner_table, || {
                    table_pass(&values, geo, table, &mut masks);
                });
                r.swar = window(&mut probe, inner_swar, || {
                    pack_signs(&values, geo, &mut inside);
                    circuit_pass::<true, false>(&inside, geo, &mut planes);
                });
            } else {
                r.swar = window(&mut probe, inner_swar, || {
                    pack_signs(&values, geo, &mut inside);
                    circuit_pass::<true, false>(&inside, geo, &mut planes);
                });
                r.table = window(&mut probe, inner_table, || {
                    table_pass(&values, geo, table, &mut masks);
                });
            }
            r.pack = window(&mut probe, inner_pack, || {
                pack_signs(&values, geo, &mut inside);
            });
            r.circuit = window(&mut probe, inner_circuit, || {
                circuit_pass::<true, false>(&inside, geo, &mut planes);
            });
            r.table_planes = window(&mut probe, inner_table_planes, || {
                table_planes_pass(&values, geo, table, &mut planes);
            });
            r.swar_masks = window(&mut probe, inner_swar_masks, || {
                pack_signs(&values, geo, &mut inside);
                swar_masks_pass(&inside, geo, &mut masks);
            });
            r.case_stream = window(&mut probe, inner_case_stream, || {
                case_stream_pass(&values, geo, &mut cases);
            });
            r.extract = window(&mut probe, inner_extract, || {
                shipped
                    .extract(sdf, &shape, origin, cell_size, &mut out)
                    .expect("extraction");
            });
            reps.push(r);
        }

        let ratios: Vec<f64> = reps
            .iter()
            .map(|r| r.swar.instructions / r.table.instructions)
            .collect();
        let spread = ratios.iter().copied().fold(0.0f64, f64::max)
            - ratios.iter().copied().fold(f64::MAX, f64::min);

        let cells = geo.cell_count();
        let per_cell = |c: Counted| Counted {
            cycles: c.cycles / cells as f64,
            instructions: c.instructions / cells as f64,
            nanos: c.nanos / cells as f64,
        };

        Row {
            field,
            resolution: n,
            cells,
            samples: geo.samples(),
            cut_cells,
            cut_fraction: cut_cells as f64 / cells as f64,
            cells_agreeing,
            pad_lanes: geo.cell_words * 64 - geo.cells,
            cell_words: geo.cell_words,
            sample_bit_row: geo.bit_row,
            output_bytes_table: cells * size_of::<u16>(),
            output_bytes_swar: geo.groups() * EDGE_COUNT * size_of::<u64>(),
            table: per_cell(median_counted(|r| r.table, &reps)),
            swar: per_cell(median_counted(|r| r.swar, &reps)),
            pack: per_cell(median_counted(|r| r.pack, &reps)),
            circuit: per_cell(median_counted(|r| r.circuit, &reps)),
            table_planes: per_cell(median_counted(|r| r.table_planes, &reps)),
            swar_masks: per_cell(median_counted(|r| r.swar_masks, &reps)),
            case_stream: per_cell(median_counted(|r| r.case_stream, &reps)),
            extract: per_cell(median_counted(|r| r.extract, &reps)),
            instruction_ratio_rep_spread: spread,
            inner_table,
            inner_swar,
            inner_extract,
            vertices,
        }
    }

    /// The registered fixture: eight reference fields × 65³, `f32`.
    ///
    /// No `scalar` column is registered and none is added: every quantity here
    /// is integer work over signs, and a sign is a sign in either precision. The
    /// *denominator* would move under `f64`, and doubling the fixture to say so
    /// would be measuring `MarchingCubes`, not this circuit.
    fn sweep(geo: Geom, table: &[u16; PATTERNS]) -> Vec<Row> {
        let mut rows = Vec::new();
        isomesh::for_each_reference_field!(f32, |name, field| {
            let (_, origin, cell_size) = crate::common::grid(&field, RESOLUTION);
            rows.push(measure(name, geo, table, &field, origin, cell_size));
        });
        rows
    }

    pub(crate) fn run(run: &mut crate::common::experiment::Run) {
        let geo = Geom::new(RESOLUTION);
        let table = edge_cut_table();

        // Fail fast, before eight rows of measurement: the pattern-level checks
        // are field-independent, so they run once and every row carries them.
        let circuit = examine_circuit(&table);
        let word_ops_total = circuit.word_ops_plane_build + circuit.word_ops_circuit;

        assert_eq!(
            circuit.patterns_tested, PATTERNS,
            "C2 is registered as exhaustive over all 256 sign patterns and ✗50 is what a sampled \
             bound costs"
        );
        assert!(
            circuit.mutant_mismatches > 0,
            "the deliberately broken circuit agreed with the table on all 256 patterns, so the \
             comparison C2 is scored on cannot fail and the 256 is worth nothing"
        );
        assert!(
            circuit.tally_agrees,
            "the counting instantiation of the circuit computes different bits from the measured \
             one, so `word_ops` would be the operation count of a different circuit"
        );
        assert_eq!(
            circuit.table_matches_cases, PATTERNS,
            "the cut-edge table derived from EDGE_CORNERS disagrees with the one read off \
             table::CASES, so the byte arm's comparand is not the crate's own topology"
        );

        let controls = examine_controls(geo, &table);
        assert_eq!(
            controls.agreeing, controls.cells,
            "the two arms disagree on the pseudorandom sign field, where every word boundary is \
             crossed — the plane build is wrong before any reference field is measured"
        );
        assert!(
            controls.carry_mutant_mismatches > 0,
            "dropping the `hi << 63` carry-in changed nothing on the pseudorandom field, so the \
             per-cell comparator cannot see a plane-build defect"
        );

        println!(
            "P-106: word_ops = {} for the circuit and {word_ops_total} including the plane build, \
             against a bar of {WORD_OP_BAR} — derived by executing the same source over a \
             counting word type",
            circuit.word_ops_circuit
        );
        println!(
            "P-106: {} of {} patterns agree with the crate's own table, packed one pattern per \
             bit lane; the broken circuit disagrees on {} of them",
            circuit.identical, circuit.patterns_tested, circuit.mutant_mismatches
        );
        println!(
            "P-106: the 256 cut masks take {} distinct values and every complement pair agrees \
             ({} of 256), so the twelve flags are strictly less information than the case index \
             and cannot feed CASES",
            circuit.distinct_cut_masks, circuit.complement_pairs_agree
        );
        println!(
            "P-106: the pseudorandom control agrees on {} of {} cells and the carry-in mutant \
             breaks {} of them",
            controls.agreeing, controls.cells, controls.carry_mutant_mismatches
        );

        let rows = sweep(geo, &table);

        let worst_ratio = rows.iter().map(Row::ratio).fold(0.0f64, f64::max);
        let best_ratio = rows.iter().map(Row::ratio).fold(f64::MAX, f64::min);
        let worst_ceiling = rows
            .iter()
            .map(Row::extraction_ceiling)
            .fold(0.0f64, f64::max);
        println!(
            "P-106: instructions per cell, each arm in its own layout — ratio {best_ratio:.4} to \
             {worst_ratio:.4} (under 1 is a win); the best extraction ceiling the saving buys is \
             {worst_ceiling:.4}x"
        );
        println!(
            "P-106: target_feature_popcnt is {} and neither arm makes a single count_ones call, \
             so no verdict here is contingent on the popcount lowering",
            cfg!(target_feature = "popcnt")
        );

        let c1_holds = circuit.word_ops_circuit <= WORD_OP_BAR;
        let c2_holds = circuit.identical == PATTERNS && circuit.patterns_tested == PATTERNS;

        for row in &rows {
            run.record(&[
                ("field", row.field.to_string()),
                ("resolution", row.resolution.to_string()),
                ("word_ops", circuit.word_ops_circuit.to_string()),
                ("masks_identical_patterns", circuit.identical.to_string()),
                ("patterns_tested", circuit.patterns_tested.to_string()),
                (
                    "mutant_pattern_mismatches",
                    circuit.mutant_mismatches.to_string(),
                ),
                (
                    "instructions_per_cell_table",
                    format!("{:.4}", row.table.instructions),
                ),
                (
                    "instructions_per_cell_swar",
                    format!("{:.4}", row.swar.instructions),
                ),
                ("ratio", format!("{:.4}", row.ratio())),
                ("c1_holds", c1_holds.to_string()),
                ("c2_holds", c2_holds.to_string()),
                ("c3_holds", (row.ratio() < 1.0).to_string()),
                // ── C1, the rest of the derivation ───────────────────────────
                (
                    "word_ops_plane_build",
                    circuit.word_ops_plane_build.to_string(),
                ),
                ("word_ops_total", word_ops_total.to_string()),
                ("word_op_bar", WORD_OP_BAR.to_string()),
                (
                    "c1_holds_with_plane_build",
                    (word_ops_total <= WORD_OP_BAR).to_string(),
                ),
                ("tally_agrees_with_u64", circuit.tally_agrees.to_string()),
                // ── C2, and the reference that licenses it ────────────────────
                (
                    "masks_identical_patterns_uniform",
                    circuit.identical_uniform.to_string(),
                ),
                (
                    "table_matches_cases_patterns",
                    circuit.table_matches_cases.to_string(),
                ),
                ("distinct_cut_masks", circuit.distinct_cut_masks.to_string()),
                (
                    "complement_pairs_agree",
                    circuit.complement_pairs_agree.to_string(),
                ),
                ("cells_agreeing", row.cells_agreeing.to_string()),
                ("cells_agreeing_random", controls.agreeing.to_string()),
                (
                    "mutant_cell_mismatches_carry",
                    controls.carry_mutant_mismatches.to_string(),
                ),
                // ── the population ───────────────────────────────────────────
                ("cells", row.cells.to_string()),
                ("samples", row.samples.to_string()),
                ("cut_cells", row.cut_cells.to_string()),
                ("cut_fraction", format!("{:.6}", row.cut_fraction)),
                ("vertices", row.vertices.to_string()),
                ("pad_lanes", row.pad_lanes.to_string()),
                ("cell_words_per_row", row.cell_words.to_string()),
                ("sample_bit_row", row.sample_bit_row.to_string()),
                ("output_bytes_table", row.output_bytes_table.to_string()),
                ("output_bytes_swar", row.output_bytes_swar.to_string()),
                // ── C3, the other three layout combinations ──────────────────
                (
                    "instructions_per_cell_table_planes",
                    format!("{:.4}", row.table_planes.instructions),
                ),
                (
                    "instructions_per_cell_swar_masks",
                    format!("{:.4}", row.swar_masks.instructions),
                ),
                (
                    "ratio_planes_layout",
                    format!("{:.4}", row.ratio_planes_layout()),
                ),
                (
                    "ratio_masks_layout",
                    format!("{:.4}", row.ratio_masks_layout()),
                ),
                (
                    "c3_holds_planes_layout",
                    (row.ratio_planes_layout() < 1.0).to_string(),
                ),
                (
                    "c3_holds_masks_layout",
                    (row.ratio_masks_layout() < 1.0).to_string(),
                ),
                // ── the SWAR arm decomposed over sibling windows ─────────────
                (
                    "instructions_per_cell_pack",
                    format!("{:.4}", row.pack.instructions),
                ),
                (
                    "instructions_per_cell_swar_given_bitmap",
                    format!("{:.4}", row.circuit.instructions),
                ),
                (
                    "pack_residual_share",
                    format!("{:.6}", row.pack_residual_share()),
                ),
                (
                    "instructions_per_cell_case_stream",
                    format!("{:.4}", row.case_stream.instructions),
                ),
                // ── SHARE, against this run's own denominator ────────────────
                (
                    "instructions_extract_mc",
                    format!("{:.4}", row.extract.instructions),
                ),
                ("cycles_extract_mc", format!("{:.4}", row.extract.cycles)),
                (
                    "case_stream_share_instructions",
                    format!("{:.6}", row.case_stream_share()),
                ),
                (
                    "cut_edge_share_instructions",
                    format!("{:.6}", row.cut_edge_share()),
                ),
                (
                    "swar_saving_share_instructions",
                    format!("{:.6}", row.swar_saving_share()),
                ),
                (
                    "extraction_ceiling",
                    format!("{:.4}", row.extraction_ceiling()),
                ),
                // ── the popcount this build does not have ────────────────────
                (
                    "target_feature_popcnt",
                    cfg!(target_feature = "popcnt").to_string(),
                ),
                ("count_ones_calls_per_cell_table", String::from("0")),
                ("count_ones_calls_per_cell_swar", String::from("0")),
                // ── provenance: reported, never consulted ────────────────────
                ("cycles_per_cell_table", format!("{:.4}", row.table.cycles)),
                ("cycles_per_cell_swar", format!("{:.4}", row.swar.cycles)),
                ("cycle_ratio", format!("{:.4}", row.cycle_ratio())),
                ("ns_per_cell_table", format!("{:.6}", row.table.nanos)),
                ("ns_per_cell_swar", format!("{:.6}", row.swar.nanos)),
                ("ghz", format!("{:.4}", row.ghz())),
                (
                    "instruction_ratio_rep_spread",
                    format!("{:.6}", row.instruction_ratio_rep_spread),
                ),
                ("reps", REPS.to_string()),
                ("inner_table", row.inner_table.to_string()),
                ("inner_swar", row.inner_swar.to_string()),
                ("inner_extract", row.inner_extract.to_string()),
            ]);
        }
    }
}

fn main() {
    if !std::env::args().any(|arg| arg == "--bench") {
        return;
    }

    let prereg = isomesh::experiment!("P-106");

    #[cfg(target_os = "linux")]
    common::experiment::run(prereg, experiment::run);

    // `experiment_p12`'s precedent. C3 is an instructions-per-cell comparison
    // and the only instrument that can read one is `perf_event_open`; a
    // nanosecond on a governed CPU cannot carry it (`M-281`), and a recorded
    // zero would be a fabricated measurement rather than a missing one.
    #[cfg(not(target_os = "linux"))]
    {
        eprintln!(
            "{} compares instructions per cell with hardware performance counters, and this \
             platform has no `perf_event_open`. There is no clock substitute.",
            prereg.id
        );
        std::process::exit(1);
    }
}
