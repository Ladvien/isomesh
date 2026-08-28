//! **P-114 — a hierarchical bitmap above the active-cell bitmap.**
//!
//! Ticket: R-114. Pre-registered before this harness existed.
//!
//! ```bash
//! cargo bench --bench experiment_p114
//! ```
//!
//! Writes `docs/experiments/p-114.csv`. **Linux only**, `experiment_p12`'s
//! precedent: the registered columns are nanoseconds, and on a governed CPU a
//! nanosecond is not a unit (`M-280`, `M-281`). Every `ns` column here is
//! accompanied by the cycle and the instruction form taken in the same window,
//! and `ghz` is on every row so a later reader can see what clock the cycles
//! were taken at. Off Linux there is nothing to degrade to, so the bench
//! refuses rather than recording a fabricated zero.
//!
//! # What was missing
//!
//! Nothing in this repository had ever put a second level above the active-cell
//! bitmap, and the reason is worth stating precisely: **this is not a new
//! structure.** `dual.rs:417` already names Museth's *VDB*
//! (`10.1145/2487228.2487235`) as the source of the flat bitmap's mechanism —
//! per-node bitmasks under word-level boolean ops — and VDB's node bitmask plus
//! its `popcount` child offset **is the same object one level up**. The paper is
//! in the corpus. So this row is a **re-measurement of a published structure at
//! this crate's granularity**, not a proposal, and the only open question is
//! whether the second level finds anything the first did not at the sizes this
//! crate actually meshes.
//!
//! What the crate has today is exactly one level. `dual.rs:359-381` packs one
//! bit per **sample** along `x`, `bit_row = size[0].div_ceil(64)` words per
//! `(y, z)` row; `dual.rs:405-436` fuses four of those rows into
//! `any & !all` — sixty-four active-cell answers in four word operations;
//! `dual.rs:445` masks the tail of the row, because the cell row is one shorter
//! than the sample row. `place_vertices` (`:487-497`) then walks
//! `cell_words = cells_x.div_ceil(64)` words per `(y, z)` cell row and tests
//! **every one of them**, whether or not the surface is anywhere near. On a
//! sphere at 65³ that is 4,096 word tests to find about 800 words that carry
//! anything.
//!
//! # The mirror, and what is inside each arm
//!
//! `crates/isomesh/src/**` is read-only for Phase 25, so the whole structure is
//! rebuilt bench-local the way `experiment_p40` rebuilds `DualMesher`'s layout
//! in its own `Grid`. [`Mirror`] mirrors `row_stride` (`dual.rs:333`),
//! `build_inside_bits` (`:359-381`), `inside_word` (`:385`),
//! `inside_word_shifted` (`:395`), `active_word` (`:424`) and `cell_mask`
//! (`:445`) line for line, and `is_inside` is inlined as `value < 0.0` because
//! `cube.rs:171-173` is `value < R::ZERO` — **exact zero is outside**.
//!
//! The clause is about the **traversal**, so the arms are drawn around the
//! traversal and nothing else:
//!
//! - **Level 0** is the materialised active-cell bitmap: one word per 64 cells,
//!   `l0_words = cell_words · cells_y · cells_z`, already masked. It is built
//!   **once** and shared by both arms, so its cost is in neither. That is what
//!   "a hierarchical bitmap **above** the active-cell bitmap" means — the
//!   active-cell bitmap is the thing that exists and the hierarchy sits on it.
//!   `ns_per_cell_build_l0` reports what it cost anyway.
//! - **The flat arm** loads all `l0_words` words in `(z, y, w)` order and walks
//!   the set bits of each with `trailing_zeros` + `a &= a - 1`, which is
//!   `dual.rs:489-497` exactly.
//! - **The hierarchical arm** puts one bit per level-0 **word** in level 1 —
//!   one bit per 64 cells, so level 1 is `l0_words.div_ceil(64)` words, a
//!   sixty-fourth the size — and one bit per level-1 word in level 2, a
//!   sixty-fourth of that. `level_1_size_ratio` and `level_2_size_ratio` are
//!   columns, so the hypothesis's first sentence is verified rather than
//!   asserted.
//!
//! Both arms emit the **packed** cell id `word_index · 64 + bit`, which is
//! division-free on both sides: the hierarchical walk reaches a level-0 word by
//! index and would otherwise have to divide by `cell_words` to recover `(w, y,
//! z)`, and paying that on one arm and not the other would be measuring the
//! division. The packed id is converted to the grid id
//! `x + cells_x · (y + cells_y · z)` **outside** every counted window, and the
//! map is monotone, so packed order and grid order are the same order.
//!
//! # Why the hierarchy cannot be built inside the arm it is measured in
//!
//! Building level 1 requires looking at every level-0 word, so a single-shot
//! build-and-traverse can never beat a single-shot traverse — the hierarchy
//! would be a strict superset of the work. That is not a defect of the
//! mechanism and it is not what the clause asks about: the registered clause is
//! *"the active-cell traversal"*, and the real workload traverses one bitmap
//! many times (re-mesh, collider, connectivity, LOD, the quad walk) for each
//! time it is built. So the summary build is measured in its **own** window and
//! reported as `ns_per_cell_build_summary`, with `breakeven_traversals` —
//! build cost divided by the per-traversal saving — saying how many traversals
//! must share one build before the mechanism is ahead. A reader who wants the
//! single-shot verdict reads that column; hiding it would be worse than
//! reporting it.
//!
//! # Counter windows are siblings, never nested
//!
//! Zen 3 has six general-purpose counters and `Probe` opens six plus a software
//! event, so two nested windows multiplex and `Counts::worst_ratio` refuses.
//! `R-121` paid for that discovery. Every window here is a **sibling** window
//! over one arm, batched over `inner_reps` repetitions so the ~28
//! `perf_event` system calls a window costs land outside it, and every quantity
//! is the median of [`REPS`] such windows.
//!
//! # SHARE
//!
//! Each clause's reachable share is a column, measured in this run rather than
//! quoted:
//!
//! - **C1's share is `traversal_share_of_extraction`** — the flat traversal's
//!   cycles per cell over the shipped `DualContouring::extract`'s cycles per
//!   cell, on the same grid, in the same run. C1 is a *stage* clause, so its
//!   own bar of 1.5× has no Amdahl ceiling; what the share bounds is what
//!   winning it would be worth to a whole extraction, and that is
//!   `extraction_speedup_ceiling` = `1/(1 − share · (1 − 1/ratio))`, computed
//!   from this row's own measured `ratio`. `M-337` measured the flat bitmap at
//!   5.5× for this stage on a surface-free field, so the stage is known to be
//!   large enough to matter; this harness measures the share itself and does
//!   not inherit that.
//! - **C1 is scored on the sparsest field as measured.** `sparsest_field` is
//!   the field with the smallest `active_fraction` at the row's own resolution,
//!   derived from this run's data and put on every row so the decision is
//!   readable from one cell. The research doc names `thin_plate`; **this
//!   harness does not inherit that**, and `sparsest_field_active_fraction` and
//!   `sparsest_field_ratio` are columns so the derivation is auditable.
//! - **C2's share is `cells_skipped_per_test`**, and it is the whole clause:
//!   cells never visited because a summary bit was zero, divided by summary
//!   words tested. It is deliberately the **measured** form rather than the
//!   structural constant. The structural form — one zero bit at the top level
//!   covers 64 level-1 bits covering 64 cells each — is `4096` by arithmetic on
//!   every grid in the fixture and therefore cannot fail, and a clause that
//!   cannot fail is worth nothing (`P-70`). It is reported beside it as
//!   `cells_per_top_level_zero_bit` and **asserted wherever it can be**, so the
//!   fan-out is verified without being mistaken for the clause. That
//!   qualification is not hedging: `fbm_terrain` at 65³ has **no** empty
//!   level-1 word at all, because at 65³ one level-1 word covers a whole `z`
//!   slice and a heightfield crosses every slice, so the row has no top-level
//!   zero bit to measure a fan-out with. It still skips cells, via empty
//!   level-0 words inside tested level-1 words, so its
//!   `cells_skipped_per_test` is real. `top_level_fanout_verified` is the
//!   column that says which rows exercised the 4,096 span, and at least one
//!   must — asserted.
//!
//!   The measured form discriminates exactly what C2 claims: the cells skipped
//!   are identical at one level and at two — every empty level-0 word is
//!   skipped either way — while the summary words tested collapse from
//!   `level_1_words` to `level_2_words + level_1_words_nonzero`. C2 is
//!   registered *"at two levels"*, so the phase verdict reads the
//!   `levels == 2` rows; the `levels == 1` rows carry the same column and the
//!   same bar so that discrimination is visible in the data instead of argued
//!   in prose.
//! - **C3 has no share; it is an equality.** Five ordered active-cell lists
//!   must be the same list: the scalar eight-corner gather (the definition),
//!   the flat word walk, the one-level walk, the two-level walk, and — the one
//!   that makes this a statement about the *shipped* mesh — the list the
//!   shipped extractor itself traverses, recovered through a bench-local
//!   [`Recording`] `VertexRule` wrapped around the public `Centroid` and
//!   injected with the public `DualContouring::with_rule`. `place_vertices`
//!   calls `rule.place` exactly once per active cell in traversal order
//!   (`dual.rs:516`), so that list *is* the shipped traversal. On top of that,
//!   `isomesh::validate::mesh_hash` of the wrapped extraction must equal the
//!   unwrapped one, which is what says the wrapper is transparent.
//!
//! # Which form carries the verdict
//!
//! `ratio` is `ns_per_cell_flat / ns_per_cell_hier`, as registered, and
//! `c1_holds` is scored on it. But **instruction counts are deterministic and
//! cycle counts on this machine are not** — `R-105` watched an identical
//! binary's cycle ratio band move from 0.984 to 1.035 across three runs while
//! its instruction counts held to four figures. So `ratio_cycles` and
//! `ratio_instructions` are beside it with their own verdict columns
//! `c1_holds_cycles` and `c1_holds_instructions`, and **the instruction form is
//! the reproducible one**: a re-run that moves `ratio` but not
//! `ratio_instructions` has moved the clock, not the mechanism.

mod common;

/// Samples per axis. `n` samples span `n − 1` cells, `common::grid`'s
/// convention, so 65 and 129 are the two grids whose cell rows are exactly one
/// and two `u64` wide — the sizes `dual.rs:472-484` singles out.
const RESOLUTIONS: [u32; 2] = [65, 129];

/// C1's bar.
const C1_BAR: f64 = 1.5;

/// C2's bar: cells skipped per summary word tested.
const C2_BAR: f64 = 4096.0;

/// Windows per quantity. Odd, so the median is a reading rather than a mean of
/// two.
const REPS: usize = 5;

#[cfg(target_os = "linux")]
mod experiment {
    use std::cell::RefCell;
    use std::hint::black_box;
    use std::rc::Rc;
    use std::time::Instant;

    use isomesh::dual::{CellVertices, VertexRule};
    use isomesh::dual_contouring::DualContouring;
    use isomesh::fields::ReferenceField;
    use isomesh::surface_nets::Centroid;
    use isomesh::validate::mesh_hash;
    use isomesh::{MeshBuffer, Real, Sdf, Shape3};

    use crate::common::counters::{MIN_TIME_RATIO, Probe};
    use crate::{C1_BAR, C2_BAR, REPS, RESOLUTIONS};

    /// About this long per counter window, so the `perf_event` round trip is
    /// negligible against what it measures.
    const TARGET_BATCH_NS: f64 = 20_000_000.0;

    /// Ceiling on the batch factor, so a grid that got faster than expected
    /// cannot turn one window into a minute.
    const MAX_INNER: usize = 1 << 20;

    // ---------------------------------------------------------------- mirror

    /// `DualMesher`'s sample grid and sign bitmap, mirrored bench-local.
    ///
    /// Every line below has a named origin in `crates/isomesh/src/dual.rs`,
    /// which is read-only for this phase. `row` is `size[0] | 1`
    /// (`row_stride`, `:333`, A-024) and the excess slots are left at zero
    /// exactly as `sample_grid` leaves them.
    struct Mirror {
        values: Vec<f64>,
        row: usize,
        size: [u32; 3],
        inside: Vec<u64>,
        bit_row: usize,
    }

    impl Mirror {
        /// Sample `field` on `[n; 3]` over its own domain.
        fn sample<F>(field: &F, n: u32) -> Self
        where
            F: ReferenceField + Sdf<Scalar = f64>,
        {
            let (shape, origin, cell_size) = crate::common::grid::<f64, F>(field, n);
            let size = shape.size();
            let row = size[0] as usize | 1;
            let nx = size[0] as usize;
            let mut values = Vec::with_capacity(row * size[1] as usize * size[2] as usize);
            for z in 0..size[2] {
                let fz = cell_size * (z as f64);
                for y in 0..size[1] {
                    let fy = cell_size * (y as f64);
                    for x in 0..nx {
                        values.push(field.sample([
                            origin[0] + cell_size * (x as f64),
                            origin[1] + fy,
                            origin[2] + fz,
                        ]));
                    }
                    values.resize(values.len() + (row - nx), 0.0);
                }
            }
            let mut mirror = Self {
                values,
                row,
                size,
                inside: Vec::new(),
                bit_row: 0,
            };
            mirror.build_inside_bits();
            mirror
        }

        /// `dual.rs:359-381`, verbatim. One bit per **sample**, 64 to a word,
        /// along `x` only.
        fn build_inside_bits(&mut self) {
            let sx = self.size[0] as usize;
            let rows = self.size[1] as usize * self.size[2] as usize;
            self.bit_row = sx.div_ceil(64);
            self.inside.clear();
            self.inside.resize(self.bit_row * rows, 0);

            for row in 0..rows {
                let src = self.row * row;
                let dst = self.bit_row * row;
                for w in 0..self.bit_row {
                    let base = w * 64;
                    let n = (sx - base).min(64);
                    let mut word = 0u64;
                    for k in 0..n {
                        // `cube.rs:171-173`: `value < R::ZERO`, so exact zero
                        // is **outside**.
                        word |= u64::from(self.values[src + base + k] < 0.0) << k;
                    }
                    self.inside[dst + w] = word;
                }
            }
        }

        /// `dual.rs:385`.
        #[inline]
        fn inside_word(&self, w: usize, y: usize, z: usize) -> u64 {
            self.inside[self.bit_row * (y + self.size[1] as usize * z) + w]
        }

        /// `dual.rs:395`. The high bit comes from the next word or the cell
        /// straddling a word boundary reads its `+x` corner as outside.
        #[inline]
        fn inside_word_shifted(&self, w: usize, y: usize, z: usize) -> u64 {
            let lo = self.inside_word(w, y, z);
            let hi = if w + 1 < self.bit_row {
                self.inside_word(w + 1, y, z)
            } else {
                0
            };
            (lo >> 1) | (hi << 63)
        }

        /// `dual.rs:424`. Sixty-four active-cell answers in four fused word
        /// operations per row.
        #[inline]
        fn active_word(&self, w: usize, y: usize, z: usize) -> u64 {
            let mut any = 0u64;
            let mut all = !0u64;
            for dz in 0..2usize {
                for dy in 0..2usize {
                    let a = self.inside_word(w, y + dy, z + dz);
                    let b = self.inside_word_shifted(w, y + dy, z + dz);
                    any |= a | b;
                    all &= a & b;
                }
            }
            any & !all
        }

        /// Cells per axis: `n` samples span `n − 1` cells.
        fn cells(&self) -> [usize; 3] {
            [
                self.size[0] as usize - 1,
                self.size[1] as usize - 1,
                self.size[2] as usize - 1,
            ]
        }

        /// Level 0: the materialised active-cell bitmap, one word per 64 cells,
        /// laid out `w`-major within `(y, z)` so that ascending word index is
        /// the shipped `(z, y, w)` traversal order.
        fn build_level0(&self) -> Vec<u64> {
            let c = self.cells();
            let cell_words = c[0].div_ceil(64);
            let mut out = vec![0u64; cell_words * c[1] * c[2]];
            for z in 0..c[2] {
                for y in 0..c[1] {
                    let rest = y + c[1] * z;
                    for w in 0..cell_words {
                        out[w + cell_words * rest] = self.active_word(w, y, z) & cell_mask(w, c[0]);
                    }
                }
            }
            out
        }

        /// The definition, for the mirror to be checked against: eight loads
        /// and eight comparisons per cell, `experiment_p40`'s `active_scalar`.
        fn active_scalar(&self) -> Vec<u32> {
            let c = self.cells();
            let mut out = Vec::new();
            for z in 0..c[2] {
                for y in 0..c[1] {
                    for x in 0..c[0] {
                        let mut inside = 0u32;
                        for corner in 0..8usize {
                            let s = (x + (corner & 1))
                                + self.row
                                    * ((y + ((corner >> 1) & 1))
                                        + self.size[1] as usize * (z + ((corner >> 2) & 1)));
                            if self.values[s] < 0.0 {
                                inside += 1;
                            }
                        }
                        if inside != 0 && inside != 8 {
                            out.push((x + c[0] * (y + c[1] * z)) as u32);
                        }
                    }
                }
            }
            out
        }
    }

    /// `dual.rs:445`. `1u64 << 64` is undefined, so the full-word case is named
    /// rather than computed.
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

    // ------------------------------------------------------------- hierarchy

    /// One bit per non-empty word of `below`, so the level above is a
    /// sixty-fourth the size.
    fn build_summary(below: &[u64]) -> Vec<u64> {
        let mut out = vec![0u64; below.len().div_ceil(64)];
        for (i, word) in below.iter().enumerate() {
            if *word != 0 {
                out[i >> 6] |= 1u64 << (i & 63);
            }
        }
        out
    }

    /// `dual.rs:489-497`: every word tested, set bits walked in ascending `x`.
    ///
    /// Returns the number of words tested. Emits the **packed** id
    /// `word · 64 + bit`.
    fn walk_flat(l0: &[u64], out: &mut Vec<u32>) -> usize {
        out.clear();
        for (i, word) in l0.iter().enumerate() {
            let base = (i as u32) << 6;
            let mut active = *word;
            while active != 0 {
                out.push(base + active.trailing_zeros());
                active &= active - 1;
            }
        }
        l0.len()
    }

    /// One summary level: every level-1 word tested, only non-empty level-0
    /// words loaded.
    fn walk_one_level(l0: &[u64], l1: &[u64], out: &mut Vec<u32>) -> usize {
        out.clear();
        let mut tested = l1.len();
        for (j, summary) in l1.iter().enumerate() {
            let mut bits = *summary;
            while bits != 0 {
                let i = (j << 6) + bits.trailing_zeros() as usize;
                bits &= bits - 1;
                tested += 1;
                let base = (i as u32) << 6;
                let mut active = l0[i];
                while active != 0 {
                    out.push(base + active.trailing_zeros());
                    active &= active - 1;
                }
            }
        }
        tested
    }

    /// Two summary levels: one zero bit at the top skips 64 level-1 bits, each
    /// covering 64 cells — 4,096 cells for one word test.
    fn walk_two_levels(l0: &[u64], l1: &[u64], l2: &[u64], out: &mut Vec<u32>) -> usize {
        out.clear();
        let mut tested = l2.len();
        for (k, top) in l2.iter().enumerate() {
            let mut hi = *top;
            while hi != 0 {
                let j = (k << 6) + hi.trailing_zeros() as usize;
                hi &= hi - 1;
                tested += 1;
                let mut bits = l1[j];
                while bits != 0 {
                    let i = (j << 6) + bits.trailing_zeros() as usize;
                    bits &= bits - 1;
                    tested += 1;
                    let base = (i as u32) << 6;
                    let mut active = l0[i];
                    while active != 0 {
                        out.push(base + active.trailing_zeros());
                        active &= active - 1;
                    }
                }
            }
        }
        tested
    }

    /// Packed id → grid id `x + cells_x · (y + cells_y · z)`.
    ///
    /// Monotone, so it does not change the order; run outside every counted
    /// window, so its division is in neither arm.
    fn to_grid_ids(packed: &[u32], cell_words: usize, cells_x: usize) -> Vec<u32> {
        packed
            .iter()
            .map(|p| {
                let i = (*p >> 6) as usize;
                let bit = (*p & 63) as usize;
                let rest = i / cell_words;
                let x = (i % cell_words) * 64 + bit;
                (x + cells_x * rest) as u32
            })
            .collect()
    }

    // ------------------------------------------------------- skip accounting

    /// What the summary levels bought, counted rather than assumed.
    #[derive(Clone, Copy)]
    struct Skip {
        /// Cells in level-0 words the traversal never loaded. **Identical at
        /// one level and at two** — every empty level-0 word is skipped either
        /// way — which is why the clause is denominated per summary word.
        cells_skipped: f64,
        /// Summary words the traversal tested.
        summary_words: usize,
        /// Zero bits at the top level of this arm.
        zero_top_bits: usize,
        /// Cells one zero top-level bit accounts for: 64 at one level, 4,096 at
        /// two. Structural, asserted, and *not* the clause.
        cells_per_top_zero_bit: f64,
    }

    /// Real cells in level-0 word `i`, after `cell_mask`.
    #[inline]
    fn cells_in_word(i: usize, cell_words: usize, cells_x: usize) -> usize {
        cell_mask(i % cell_words, cells_x).count_ones() as usize
    }

    fn account(
        levels: u8,
        l0: &[u64],
        l1: &[u64],
        l2: &[u64],
        cell_words: usize,
        cells_x: usize,
    ) -> Skip {
        let empty_cells = |i: usize| cells_in_word(i, cell_words, cells_x) as f64;
        let cells_skipped: f64 = l0
            .iter()
            .enumerate()
            .filter(|(_, word)| **word == 0)
            .map(|(i, _)| empty_cells(i))
            .sum();

        if levels == 1 {
            let zero_top_bits = l0.iter().filter(|w| **w == 0).count();
            return Skip {
                cells_skipped,
                summary_words: l1.len(),
                zero_top_bits,
                cells_per_top_zero_bit: if zero_top_bits == 0 {
                    0.0
                } else {
                    cells_skipped / zero_top_bits as f64
                },
            };
        }

        // At two levels the top is level 1, and a zero bit there stands for a
        // whole empty level-1 word: 64 level-0 words, 4,096 cells.
        let mut zero_top_bits = 0usize;
        let mut cells_via_top = 0.0f64;
        for (j, summary) in l1.iter().enumerate() {
            if *summary != 0 {
                continue;
            }
            zero_top_bits += 1;
            let block = j << 6;
            for (offset, word) in l0[block..(block + 64).min(l0.len())].iter().enumerate() {
                debug_assert_eq!(*word, 0, "a summary bit was zero over a non-empty word");
                cells_via_top += empty_cells(block + offset);
            }
        }
        Skip {
            cells_skipped,
            summary_words: l2.len() + l1.iter().filter(|w| **w != 0).count(),
            zero_top_bits,
            cells_per_top_zero_bit: if zero_top_bits == 0 {
                0.0
            } else {
                cells_via_top / zero_top_bits as f64
            },
        }
    }

    // ------------------------------------------------------------- the clock

    /// One arm, as measured.
    #[derive(Clone, Copy, Default)]
    struct Arm {
        ns_per_cell: f64,
        cycles_per_cell: f64,
        instructions_per_cell: f64,
        words_tested: usize,
    }

    /// Cycles, instructions and nanoseconds from one counter window.
    #[derive(Clone, Copy, Default)]
    struct Counted {
        cycles: f64,
        instructions: f64,
        nanos: f64,
    }

    /// One counter window over `inner` repetitions of `body`, divided by
    /// `inner`.
    ///
    /// Sibling, never nested: Zen 3 has six general-purpose counters and
    /// `Probe` opens six, so a window inside a window multiplexes and
    /// `worst_ratio` refuses (`R-121`).
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
            nanos: nanos * scale,
        }
    }

    fn median(mut values: Vec<f64>) -> f64 {
        values.sort_by(f64::total_cmp);
        values[values.len() / 2]
    }

    /// Repetitions to make one window about [`TARGET_BATCH_NS`] long.
    fn pick_inner(pass_ns: f64) -> usize {
        ((TARGET_BATCH_NS / pass_ns.max(1.0)).ceil() as usize).clamp(1, MAX_INNER)
    }

    /// Time one pass, to size the batch.
    fn probe_pass_ns(mut body: impl FnMut()) -> f64 {
        body();
        let started = Instant::now();
        body();
        started.elapsed().as_nanos() as f64
    }

    /// The three traversal arms, measured **interleaved**.
    ///
    /// This is the fix for a defect the first run of this harness exposed and
    /// `R-105` had already named. Measuring five windows of the flat arm, then
    /// five of the one-level arm, then five of the two-level arm gives the flat
    /// arm the clock the machine happened to have at the start of the row, and
    /// on a governed CPU that is not the clock the other two get: across two
    /// runs of the identical binary `ratio` moved from 1.215 to 0.681 on
    /// `fbm_terrain` at 129³ while `ratio_instructions` reproduced to four
    /// figures at 1.418 both times. So each repetition takes **one window per
    /// arm, back to back**, and the medians are taken per quantity per arm
    /// afterwards; a clock drifting over a row now moves all three arms
    /// together and largely cancels in the ratio.
    ///
    /// All three arms share one `inner`, sized on the flat arm. That is
    /// deliberate: the batching overhead is then identical across the arms and
    /// cancels in the ratio, at the price of a shorter window on the faster
    /// arm — still milliseconds against microseconds of `perf_event` calls.
    fn measure_arms(
        probe: &mut Probe,
        cells: f64,
        inner: usize,
        arms: &mut [&mut dyn FnMut(); 3],
    ) -> [Arm; 3] {
        let mut counted: [Vec<Counted>; 3] = [
            Vec::with_capacity(REPS),
            Vec::with_capacity(REPS),
            Vec::with_capacity(REPS),
        ];
        for _ in 0..REPS {
            for (slot, body) in counted.iter_mut().zip(arms.iter_mut()) {
                slot.push(window(probe, inner, &mut **body));
            }
        }
        let arm = |c: &Vec<Counted>| Arm {
            ns_per_cell: median(c.iter().map(|w| w.nanos / cells).collect()),
            cycles_per_cell: median(c.iter().map(|w| w.cycles / cells).collect()),
            instructions_per_cell: median(c.iter().map(|w| w.instructions / cells).collect()),
            words_tested: 0,
        };
        [arm(&counted[0]), arm(&counted[1]), arm(&counted[2])]
    }

    // ------------------------------------------------------- the shipped list

    /// The shipped extractor's own ordered active-cell list.
    ///
    /// `place_vertices` calls `rule.place` once per active cell, in traversal
    /// order, **before** it knows whether the rule will produce anything
    /// (`dual.rs:516`), so wrapping the rule recovers the traversal itself
    /// rather than a filtered view of it. `Centroid` is the inner rule because
    /// it is public and because the recorded list is rule-independent — the
    /// walk that produces it is above the placement.
    struct Recording<V> {
        inner: V,
        log: Rc<RefCell<Vec<u32>>>,
        cells: [u32; 3],
    }

    impl<R: Real, V: VertexRule<R>> VertexRule<R> for Recording<V> {
        fn place<S: Sdf<Scalar = R>>(
            &self,
            sdf: &S,
            corner: &[R; 8],
            base: [u32; 3],
            origin: [R; 3],
            cell_size: R,
            out: &mut CellVertices<R>,
        ) {
            self.log
                .borrow_mut()
                .push(base[0] + self.cells[0] * (base[1] + self.cells[1] * base[2]));
            self.inner.place(sdf, corner, base, origin, cell_size, out);
        }
    }

    /// The shipped traversal's ordered list, plus whether the recording wrapper
    /// left the mesh alone.
    fn shipped_list<F>(field: &F, n: u32, cells: [usize; 3]) -> (Vec<u32>, bool)
    where
        F: ReferenceField + Sdf<Scalar = f64>,
    {
        let (shape, origin, cell_size) = crate::common::grid::<f64, F>(field, n);
        let mut plain = MeshBuffer::<f64>::new();
        DualContouring::<f64, _>::with_rule(Centroid)
            .extract(field, &shape, origin, cell_size, &mut plain)
            .expect("extraction");

        let log = Rc::new(RefCell::new(Vec::new()));
        let rule = Recording {
            inner: Centroid,
            log: Rc::clone(&log),
            cells: [cells[0] as u32, cells[1] as u32, cells[2] as u32],
        };
        let mut wrapped = MeshBuffer::<f64>::new();
        DualContouring::<f64, _>::with_rule(rule)
            .extract(field, &shape, origin, cell_size, &mut wrapped)
            .expect("extraction");

        let transparent = mesh_hash(&plain) == mesh_hash(&wrapped);
        (
            Rc::try_unwrap(log).expect("sole owner").into_inner(),
            transparent,
        )
    }

    /// The shipped extraction's cycles per cell, for C1's share.
    fn shipped_cycles_per_cell<F>(probe: &mut Probe, field: &F, n: u32, cells: f64) -> f64
    where
        F: ReferenceField + Sdf<Scalar = f64>,
    {
        let (shape, origin, cell_size) = crate::common::grid::<f64, F>(field, n);
        let mut mesher = DualContouring::<f64>::new();
        let mut out = MeshBuffer::<f64>::new();
        let mut body = || {
            out.reset();
            mesher
                .extract(field, &shape, origin, cell_size, &mut out)
                .expect("extraction");
            black_box(&out);
        };
        body();
        let mut runs = Vec::with_capacity(REPS);
        for _ in 0..REPS {
            runs.push(window(probe, 1, &mut body).cycles / cells);
        }
        median(runs)
    }

    // ------------------------------------------------------------- one fixture

    /// Everything one `(field, resolution)` produced.
    struct Measured {
        field: &'static str,
        resolution: u32,
        cells: usize,
        active_cells: usize,
        active_fraction: f64,
        l0_words: usize,
        l1_words: usize,
        l2_words: usize,
        l0_nonzero: usize,
        l1_nonzero: usize,
        flat: Arm,
        /// Index 0 is one summary level, index 1 is two.
        hier: [Arm; 2],
        skip: [Skip; 2],
        build_summary: Counted,
        build_l0_ns_per_cell: f64,
        shipped_cycles_per_cell: f64,
        mesh_identical: bool,
        inner_reps: usize,
    }

    fn measure<F>(name: &'static str, field: &F, n: u32, probe: &mut Probe) -> Measured
    where
        F: ReferenceField + Sdf<Scalar = f64>,
    {
        let mirror = Mirror::sample(field, n);
        let c = mirror.cells();
        let cells = c[0] * c[1] * c[2];
        let cells_f = cells as f64;
        let cell_words = c[0].div_ceil(64);

        let build_l0_ns_per_cell = {
            let mut runs = Vec::with_capacity(REPS);
            for _ in 0..REPS {
                let started = Instant::now();
                let l0 = mirror.build_level0();
                runs.push(started.elapsed().as_nanos() as f64 / cells_f);
                black_box(&l0);
            }
            median(runs)
        };

        let l0 = mirror.build_level0();
        let l1 = build_summary(&l0);
        let l2 = build_summary(&l1);
        assert_eq!(l0.len(), cell_words * c[1] * c[2], "level-0 word count");
        assert_eq!(
            l1.len(),
            l0.len().div_ceil(64),
            "level 1 is a 64th of level 0"
        );
        assert_eq!(
            l2.len(),
            l1.len().div_ceil(64),
            "level 2 is a 64th of level 1"
        );

        let l0_nonzero = l0.iter().filter(|w| **w != 0).count();
        let l1_nonzero = l1.iter().filter(|w| **w != 0).count();

        // ---- the lists, and C3.
        let mut packed = Vec::with_capacity(l0.iter().map(|w| w.count_ones() as usize).sum());
        let words_flat = walk_flat(&l0, &mut packed);
        let flat_ids = to_grid_ids(&packed, cell_words, c[0]);
        let words_one = walk_one_level(&l0, &l1, &mut packed);
        let one_ids = to_grid_ids(&packed, cell_words, c[0]);
        let words_two = walk_two_levels(&l0, &l1, &l2, &mut packed);
        let two_ids = to_grid_ids(&packed, cell_words, c[0]);

        assert_eq!(
            words_flat,
            l0.len(),
            "the flat arm tests every level-0 word"
        );
        assert_eq!(
            words_one,
            l1.len() + l0_nonzero,
            "one level tests every level-1 word plus the non-empty level-0 words"
        );
        assert_eq!(
            words_two,
            l2.len() + l1_nonzero + l0_nonzero,
            "two levels test every level-2 word plus the non-empty words below"
        );

        let scalar_ids = mirror.active_scalar();
        let (shipped_ids, mesh_transparent) = shipped_list(field, n, c);
        assert_eq!(
            scalar_ids, flat_ids,
            "{name} {n}³: the word mirror disagrees with the eight-corner definition"
        );
        assert_eq!(
            shipped_ids, flat_ids,
            "{name} {n}³: the mirror's traversal is not the shipped extractor's"
        );
        assert!(
            mesh_transparent,
            "{name} {n}³: the recording wrapper moved the mesh"
        );
        let mesh_identical = one_ids == flat_ids && two_ids == flat_ids;
        assert!(
            mesh_identical,
            "{name} {n}³: a summary level changed the active-cell list, so C3's equality \
             is broken rather than merely unmet"
        );

        let active_cells = flat_ids.len();

        // ---- the clock. One sibling window per arm per repetition, the three
        // arms interleaved so a drifting clock cannot land on one of them.
        // Each arm owns its own output buffer, because three closures cannot
        // borrow one buffer at once — and because a per-arm buffer keeps each
        // arm's allocation warm the way a real consumer's would be.
        let mut out_flat = Vec::with_capacity(active_cells);
        let mut out_one = Vec::with_capacity(active_cells);
        let mut out_two = Vec::with_capacity(active_cells);
        let inner = pick_inner(probe_pass_ns(|| {
            black_box(walk_flat(&l0, &mut out_flat));
        }));
        let [mut flat, mut hier_one, mut hier_two] = measure_arms(
            probe,
            cells_f,
            inner,
            &mut [
                &mut || {
                    black_box(walk_flat(&l0, &mut out_flat));
                },
                &mut || {
                    black_box(walk_one_level(&l0, &l1, &mut out_one));
                },
                &mut || {
                    black_box(walk_two_levels(&l0, &l1, &l2, &mut out_two));
                },
            ],
        );
        flat.words_tested = words_flat;
        hier_one.words_tested = words_one;
        hier_two.words_tested = words_two;

        let build_inner = pick_inner(probe_pass_ns(|| {
            black_box(build_summary(&build_summary(&l0)));
        }));
        let build_summary_counted = {
            let mut runs = Vec::with_capacity(REPS);
            for _ in 0..REPS {
                runs.push(window(probe, build_inner, || {
                    let a = build_summary(&l0);
                    black_box(build_summary(&a));
                }));
            }
            Counted {
                cycles: median(runs.iter().map(|c| c.cycles / cells_f).collect()),
                instructions: median(runs.iter().map(|c| c.instructions / cells_f).collect()),
                nanos: median(runs.iter().map(|c| c.nanos / cells_f).collect()),
            }
        };

        let shipped = shipped_cycles_per_cell(probe, field, n, cells_f);

        Measured {
            field: name,
            resolution: n,
            cells,
            active_cells,
            active_fraction: active_cells as f64 / cells_f,
            l0_words: l0.len(),
            l1_words: l1.len(),
            l2_words: l2.len(),
            l0_nonzero,
            l1_nonzero,
            flat,
            hier: [hier_one, hier_two],
            skip: [
                account(1, &l0, &l1, &l2, cell_words, c[0]),
                account(2, &l0, &l1, &l2, cell_words, c[0]),
            ],
            build_summary: build_summary_counted,
            build_l0_ns_per_cell,
            shipped_cycles_per_cell: shipped,
            mesh_identical,
            inner_reps: inner,
        }
    }

    // ------------------------------------------------------------------- run

    pub(crate) fn run(run: &mut crate::common::experiment::Run) {
        let mut probe = Probe::open();
        let mut measured: Vec<Measured> = Vec::new();
        isomesh::for_each_reference_field!(f64, |name, field| {
            for n in RESOLUTIONS {
                let row = measure(name, &field, n, &mut probe);
                println!(
                    "{name:>14} {n:>4}³  active {:>7.4}%  words {:>7} → {:>7}  ratio {:>6.3}× \
                     (instr {:>6.3}×)  skipped/test {:>9.1}",
                    row.active_fraction * 100.0,
                    row.flat.words_tested,
                    row.hier[1].words_tested,
                    row.flat.ns_per_cell / row.hier[1].ns_per_cell,
                    row.flat.instructions_per_cell / row.hier[1].instructions_per_cell,
                    row.skip[1].cells_skipped / row.skip[1].summary_words as f64,
                );
                measured.push(row);
            }
        });

        // VACUITY CONTROL, asserted rather than recorded: a hierarchy that
        // skips nothing is not being tested.
        assert!(
            measured
                .iter()
                .any(|m| m.hier[1].words_tested < m.flat.words_tested),
            "no row tested fewer words with two levels than flat, so the hierarchy skipped \
             nothing and this fixture cannot see the mechanism"
        );
        assert!(
            measured.iter().all(|m| m.skip[1].cells_skipped > 0.0),
            "a row skipped zero cells at two levels"
        );

        // The registration's own fan-out claim, verified structurally rather
        // than assumed: one zero bit at the top of two levels stands for a
        // whole empty level-1 word — 64 level-0 words of 64 cells, 4,096 cells
        // for one word test — and one zero bit at one level stands for 64.
        //
        // Conditional per row, because a field can genuinely have **no** empty
        // level-1 word: at 65³ a level-1 word covers one whole `z` slice, and
        // `fbm_terrain` is a heightfield whose surface crosses every slice. That
        // row still skips cells (via empty level-0 words inside tested level-1
        // words), so its `cells_skipped_per_test` is real; what it cannot do is
        // testify about the top-level fan-out, and reporting 4,096 there would
        // be inventing a measurement. `top_level_fanout_verified` is the column
        // that says which rows did.
        let mut fanout_rows = 0usize;
        for m in &measured {
            if m.skip[1].zero_top_bits > 0 {
                fanout_rows += 1;
                assert!(
                    (m.skip[1].cells_per_top_zero_bit - C2_BAR).abs() < 0.5,
                    "{} {}³: a zero top-level bit covered {} cells, not 4096",
                    m.field,
                    m.resolution,
                    m.skip[1].cells_per_top_zero_bit
                );
            }
            if m.skip[0].zero_top_bits > 0 {
                assert!(
                    (m.skip[0].cells_per_top_zero_bit - 64.0).abs() < 0.5,
                    "{} {}³: a zero one-level bit covered {} cells, not 64",
                    m.field,
                    m.resolution,
                    m.skip[0].cells_per_top_zero_bit
                );
            }
        }
        assert!(
            fanout_rows > 0,
            "no row had an empty level-1 word, so the 4,096-cell fan-out was never exercised"
        );

        // C1's field, from this run's own active_fraction. Per resolution,
        // because "the sparsest field as measured" is measured on a grid.
        let sparsest = |n: u32| -> (&'static str, f64) {
            let pick = measured
                .iter()
                .filter(|m| m.resolution == n)
                .min_by(|a, b| a.active_fraction.total_cmp(&b.active_fraction))
                .expect("every resolution has rows");
            (pick.field, pick.active_fraction)
        };

        for m in &measured {
            let (sparsest_field, sparsest_fraction) = sparsest(m.resolution);
            for (index, levels) in [(0usize, 1u8), (1usize, 2u8)] {
                let hier = m.hier[index];
                let skip = m.skip[index];
                let ratio = m.flat.ns_per_cell / hier.ns_per_cell;
                let ratio_cycles = m.flat.cycles_per_cell / hier.cycles_per_cell;
                let ratio_instructions = m.flat.instructions_per_cell / hier.instructions_per_cell;

                // C1's verdict belongs to the clause, not to the row: it is
                // the sparsest field's ratio at this resolution and this level
                // count, carried on every row so the decision reads from one
                // cell (`P-121`'s idiom).
                let clause_row = measured
                    .iter()
                    .find(|other| other.field == sparsest_field && other.resolution == m.resolution)
                    .expect("the sparsest field is one of the rows");
                let sparsest_ratio =
                    clause_row.flat.ns_per_cell / clause_row.hier[index].ns_per_cell;
                let sparsest_ratio_cycles =
                    clause_row.flat.cycles_per_cell / clause_row.hier[index].cycles_per_cell;
                let sparsest_ratio_instructions = clause_row.flat.instructions_per_cell
                    / clause_row.hier[index].instructions_per_cell;

                let cells_skipped_per_test = skip.cells_skipped / skip.summary_words as f64;
                let share = m.flat.cycles_per_cell / m.shipped_cycles_per_cell;
                let saved = share * (1.0 - 1.0 / ratio);
                let ceiling = if saved < 1.0 {
                    1.0 / (1.0 - saved)
                } else {
                    f64::INFINITY
                };
                let per_traversal_saving = m.flat.ns_per_cell - hier.ns_per_cell;
                let breakeven = if per_traversal_saving > 0.0 {
                    m.build_summary.nanos / per_traversal_saving
                } else {
                    f64::INFINITY
                };

                run.record(&[
                    ("field", m.field.to_string()),
                    ("resolution", m.resolution.to_string()),
                    ("active_fraction", format!("{:.8}", m.active_fraction)),
                    ("levels", levels.to_string()),
                    ("words_tested_flat", m.flat.words_tested.to_string()),
                    ("words_tested_hier", hier.words_tested.to_string()),
                    (
                        "cells_skipped_per_test",
                        format!("{cells_skipped_per_test:.3}"),
                    ),
                    ("ns_per_cell_flat", format!("{:.6}", m.flat.ns_per_cell)),
                    ("ns_per_cell_hier", format!("{:.6}", hier.ns_per_cell)),
                    ("ratio", format!("{ratio:.4}")),
                    ("mesh_identical", m.mesh_identical.to_string()),
                    ("sparsest_field", sparsest_field.to_string()),
                    ("c1_holds", (sparsest_ratio >= C1_BAR).to_string()),
                    ("c2_holds", (cells_skipped_per_test >= C2_BAR).to_string()),
                    ("c3_holds", m.mesh_identical.to_string()),
                    // --- extras (M-273): the cycle and instruction forms, the
                    // shares, the structure's own sizes, and what the summary
                    // cost to build.
                    (
                        "ghz",
                        format!("{:.4}", m.flat.cycles_per_cell / m.flat.ns_per_cell),
                    ),
                    (
                        "ghz_hier",
                        format!("{:.4}", hier.cycles_per_cell / hier.ns_per_cell),
                    ),
                    // `ratio / ratio_cycles`. One means the two windows ran at
                    // the same clock and the nanosecond ratio is the cycle
                    // ratio; anything else is the governor showing up in a
                    // registered column, which is the whole of `M-281`.
                    (
                        "clock_ratio_hier_over_flat",
                        format!("{:.4}", ratio / ratio_cycles),
                    ),
                    ("cells", m.cells.to_string()),
                    ("active_cells", m.active_cells.to_string()),
                    (
                        "word_ratio",
                        format!(
                            "{:.4}",
                            m.flat.words_tested as f64 / hier.words_tested as f64
                        ),
                    ),
                    (
                        "cycles_per_cell_flat",
                        format!("{:.6}", m.flat.cycles_per_cell),
                    ),
                    (
                        "cycles_per_cell_hier",
                        format!("{:.6}", hier.cycles_per_cell),
                    ),
                    ("ratio_cycles", format!("{ratio_cycles:.4}")),
                    (
                        "instructions_per_cell_flat",
                        format!("{:.6}", m.flat.instructions_per_cell),
                    ),
                    (
                        "instructions_per_cell_hier",
                        format!("{:.6}", hier.instructions_per_cell),
                    ),
                    ("ratio_instructions", format!("{ratio_instructions:.4}")),
                    (
                        "c1_holds_cycles",
                        (sparsest_ratio_cycles >= C1_BAR).to_string(),
                    ),
                    (
                        "c1_holds_instructions",
                        (sparsest_ratio_instructions >= C1_BAR).to_string(),
                    ),
                    (
                        "sparsest_field_active_fraction",
                        format!("{sparsest_fraction:.8}"),
                    ),
                    ("sparsest_field_ratio", format!("{sparsest_ratio:.4}")),
                    (
                        "sparsest_field_ratio_instructions",
                        format!("{sparsest_ratio_instructions:.4}"),
                    ),
                    ("cells_skipped_total", format!("{:.0}", skip.cells_skipped)),
                    ("summary_words_tested", skip.summary_words.to_string()),
                    ("zero_top_level_bits", skip.zero_top_bits.to_string()),
                    (
                        "cells_per_top_level_zero_bit",
                        format!("{:.1}", skip.cells_per_top_zero_bit),
                    ),
                    (
                        "top_level_fanout_verified",
                        (skip.zero_top_bits > 0).to_string(),
                    ),
                    ("level_0_words", m.l0_words.to_string()),
                    ("level_1_words", m.l1_words.to_string()),
                    ("level_2_words", m.l2_words.to_string()),
                    ("level_0_words_nonzero", m.l0_nonzero.to_string()),
                    ("level_1_words_nonzero", m.l1_nonzero.to_string()),
                    (
                        "level_1_size_ratio",
                        format!("{:.4}", m.l0_words as f64 / m.l1_words as f64),
                    ),
                    (
                        "level_2_size_ratio",
                        format!("{:.4}", m.l1_words as f64 / m.l2_words as f64),
                    ),
                    ("summary_bytes", ((m.l1_words + m.l2_words) * 8).to_string()),
                    ("level_0_bytes", (m.l0_words * 8).to_string()),
                    (
                        "summary_overhead_fraction",
                        format!(
                            "{:.6}",
                            (m.l1_words + m.l2_words) as f64 / m.l0_words as f64
                        ),
                    ),
                    (
                        "ns_per_cell_build_summary",
                        format!("{:.6}", m.build_summary.nanos),
                    ),
                    (
                        "cycles_per_cell_build_summary",
                        format!("{:.6}", m.build_summary.cycles),
                    ),
                    (
                        "instructions_per_cell_build_summary",
                        format!("{:.6}", m.build_summary.instructions),
                    ),
                    (
                        "ns_per_cell_build_level_0",
                        format!("{:.6}", m.build_l0_ns_per_cell),
                    ),
                    ("breakeven_traversals", format!("{breakeven:.3}")),
                    (
                        "cycles_per_cell_shipped_extraction",
                        format!("{:.4}", m.shipped_cycles_per_cell),
                    ),
                    ("traversal_share_of_extraction", format!("{share:.6}")),
                    ("extraction_speedup_ceiling", format!("{ceiling:.5}")),
                    ("inner_reps", m.inner_reps.to_string()),
                    ("reps", REPS.to_string()),
                ]);
            }
        }
    }
}

fn main() {
    if !std::env::args().any(|arg| arg == "--bench") {
        return;
    }

    let prereg = isomesh::experiment!("P-114");

    #[cfg(target_os = "linux")]
    common::experiment::run(prereg, experiment::run);

    // `experiment_p12`'s precedent. Every registered column here is a
    // nanosecond, and `M-280`/`M-281` forbid a nanosecond carrying a verdict on
    // a governed CPU without the cycle and instruction forms beside it — which
    // come from `perf_event_open`. Off Linux there is nothing to degrade to,
    // and a recorded zero would be a fabricated measurement.
    #[cfg(not(target_os = "linux"))]
    {
        eprintln!(
            "{} scores a traversal ratio, and on a governed CPU that verdict needs the cycle \
             and instruction forms from `perf_event_open`, which this platform does not have.",
            prereg.id
        );
        std::process::exit(1);
    }
}
