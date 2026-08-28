//! **P-113 — Roaring's density thresholds as the chunk-representation decision rule.**
//!
//! Ticket: R-113. Pre-registered before this harness existed.
//!
//! ```bash
//! cargo bench --bench experiment_p113
//! ```
//!
//! Writes `docs/experiments/p-113.csv`.
//!
//! # What was missing
//!
//! **The density.** Roaring switches container type at 4,096 set values per 2¹⁶
//! chunk — 6.25% — and says outright that *"when applications encounter integer
//! sets with lower density (less than 0.1%), a bitmap is unlikely to be the
//! proper data structure."* Applying that rule to this crate needs one number
//! the crate has never measured: **what fraction of a chunk's cells are active.**
//!
//! Two numbers in the ledger look like that number and are not it.
//!
//! - `M-306`'s **16.8%** (gyroid, 688 of 4,096) and **95.1%** (thin_plate, 3,896
//!   of 4,096) are **rejected-brush shares, not active-cell densities**. They
//!   count brushes a chunk's bound rejected, they were measured at **17³ with
//!   `SubgridMarchingTetrahedra` at 16 samples per edge**, and they were measured
//!   **in a debug build**. Nothing about them is a cell count.
//! - The research doc's inference that *"thin_plate's active density is around
//!   4.9%"* is derived from those two shares, so it inherits every one of their
//!   defects. **It must not be quoted as a measurement.** This harness measures
//!   density itself, from the eight-corner sign rule, at 128³, in a release
//!   build, and reports the number it read.
//!
//! What is also missing is the *second* representation. `dual.rs:359-381` builds
//! a bitmap and nothing in the crate builds an array of active-cell indices, so
//! the array arm here is bench-local — as is the bitmap arm, mirrored from
//! `dual.rs:405-436`'s fused `any & !all` fold rather than reached through
//! `DualMesher`, which is `pub(crate)`.
//!
//! # SHARE
//!
//! **No extraction time is claimed, so `✗51`'s `1/(1 − share/factor)` bar does
//! not apply to any clause here.** This row produces a *decision rule*, and each
//! clause's reachable share is a population that is a **column**:
//!
//! - **C1's share is `fields_below_boundary` + `fields_above_boundary`.** The
//!   clause is *"at least one field below and one above"*, so its population is
//!   the eight reference fields at one granularity, and the two counts either
//!   side of 6.25% are exactly the share the clause can see. They sum to 8 unless
//!   a field lands *exactly* on the boundary, which is neither below nor above
//!   and counts for neither side — and one does, so the sum is 7 on that
//!   granularity's rows and the harness prints which field and why.
//!   Its falsifier is a second population: `distinct_winners`, the number of
//!   distinct values the `winner` column takes over those eight fields. One
//!   representation winning on all eight is `distinct_winners == 1`, and that
//!   closes the row.
//! - **C2's share is `rows_below_boundary`.** *"Below the threshold"* is a
//!   sub-population of the fixture's 32 rows, and if it is empty the clause is
//!   **VACUOUS** rather than held — a comparison with no rows to compare is the
//!   `M-44` failure wearing a boolean. The column is there so the reader can see
//!   how many rows carried the verdict.
//! - **C3's share is 8 rows: the eight fields at `chunk_cells = 4`.** At
//!   `M-377`'s optimum a chunk is 64 cells, which is exactly one `u64`, so the
//!   payload is 8 bytes and the overhead is whatever a chunk actually carries
//!   around that word. `bitmap_overhead_bytes_per_active_cell` and
//!   `bitmap_payload_bytes_per_active_cell` are both columns, and C3 is the
//!   comparison between them.
//!
//! # The two representations, and where their crossover comes from
//!
//! Per **occupied** chunk — a chunk with at least one active cell:
//!
//! | arm | payload | overhead |
//! |---|---|---|
//! | bitmap | `c³ / 64` words of 8 bytes | one `Vec<u64>` header |
//! | array | `active` indices of `index_bytes` | one `Vec<W>` header |
//!
//! **Occupied only, and that is a fairness decision rather than a convenience.**
//! An empty chunk is stored as no container at all in Roaring, and in both arms
//! here it is stored as nothing, so it cancels; and both arms need the identical
//! chunk-lookup structure above them, so that cancels too. Counting empty chunks
//! would price a chunk-index structure that neither arm owns. `chunks_total` and
//! `chunks_occupied` are both columns, so the other reading is derivable.
//!
//! The index width is chosen from the chunk size, which is the only thing that
//! bounds a chunk-local index: `u8` at 4³ (64 cells), `u16` at 8³ and 16³ (512
//! and 4,096), `u32` at 64³ (262,144). That choice is where Roaring's own
//! constant comes from, and the harness records it as `byte_crossover_density`:
//! the array is smaller than the bitmap exactly when
//! `active * index_bytes < c³ / 8`, i.e. below a density of `1 / (8 *
//! index_bytes)` — **12.5%** at one byte, **6.25%** at two, **3.125%** at four.
//! Roaring's 4,096-of-2¹⁶ *is* the two-byte case: a 16-bit index against one
//! bit is 16×, and 1/16 is 6.25%. So at 8³ and 16³ this fixture's byte crossover
//! and Roaring's container boundary are the same number, arrived at
//! independently.
//!
//! # The walk, and why it is denominated per active cell
//!
//! Both arms walk every occupied chunk in the same chunk order and every active
//! cell in the same ascending chunk-local order — the bitmap by
//! `trailing_zeros` + `x &= x - 1` (`dual.rs:495-497`'s walk), the array by
//! iteration. The consumer is the same non-vectorisable rotate-and-**add** in
//! both, which makes it a dependency chain rather than a reducible sum and makes
//! the result **order-sensitive**: `walk_checksum` doubles as the control that
//! the two arms visited the same cells in the same order, and it is asserted
//! equal and non-zero.
//!
//! Per **active cell**, not per word or per chunk, because that is the work the
//! consumer wanted done. The bitmap also scans the `c³/64 − popcount` words that
//! hold nothing, and that cost is real and belongs in its per-active-cell
//! figure; `words_scanned` is a column so the reader can see it.
//!
//! Two things about that consumer are findings rather than choices, and both are
//! columns:
//!
//! - **It was `^` first, and the harness refused to run.** The non-zero assertion
//!   on `walk_checksum` fired immediately: an XOR chain over a symmetric fixture
//!   cancelled to *exactly* zero, on both arms, at 4³. XOR is GF(2)-linear, so
//!   two visits to the same index at positions congruent mod 64 annihilate — and
//!   a sphere's 4³ chunks repeat the same 6-bit pattern by the thousand.
//!   `wrapping_add` carries, so it does not. See [`consume`].
//! - **The array arm is saturated at the loop's floor and its column is a
//!   ceiling.** `walk_floor_cycles_per_active_cell` measures the *identical* walk
//!   over one contiguous ramp of the same length — the cheapest this consumer can
//!   be driven at all — and the array arm reads it, on every row, at every index
//!   width. So the array arm's own cost is not resolved. That does **not** void
//!   C2's walk half: the array is the arm that *wins*, and a ceiling that already
//!   loses bounds the true cost from above, so `array < bitmap` holds a fortiori.
//!   The reverse — a saturated arm *winning* — would be void, and
//!   `array_walk_saturated` and `bitmap_walk_saturated` are the columns that let
//!   a reader check which case a row is.
//!
//! # Why this bench refuses to run off Linux
//!
//! `M-280` and `M-281`: on a governed CPU a nanosecond is not a unit. So the
//! walk comparison that decides C2 is scored on **cycles per active cell** from
//! `common::counters::Probe`, with `ns_per_walk_*` reported beside it and a
//! `ghz` column on the same row, taken from the same counted span as the cycles
//! it divides. `perf_event_open` is Linux-only (`benches/common/mod.rs:26`), and
//! the only honest thing to do elsewhere is refuse: a zero in `ghz` is a number
//! nobody measured. `experiment_p12` is the precedent.

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

    use isomesh::fields::ReferenceField;
    use isomesh::{Sdf, for_each_reference_field};

    use crate::common::counters::{MIN_TIME_RATIO, Probe};
    use crate::common::experiment::Run;

    /// `f64`, matching `experiment_p72`, whose fixture shape this reuses. A
    /// density is a count of sign changes and does not depend on the scalar, but
    /// comparing 128³ cell counts against `M-377`'s requires the same grid.
    type Scalar = f64;

    /// Cells per axis in the whole world, at every granularity.
    ///
    /// `experiment_p72`'s constant, for its reason: 128 is divisible by every
    /// granularity in the sweep, so no arm compares a different cell count.
    const WORLD_CELLS: usize = 128;

    /// Samples per axis. `n` samples span `n − 1` cells.
    const SAMPLES: usize = WORLD_CELLS + 1;

    /// Words per sample row in the sign bit-plane. `dual.rs:362`'s
    /// `size[0].div_ceil(64)` — 129 samples need **3** words while 128 cells
    /// need 2, which is the cell/word asymmetry `P-104` is registered against.
    const SAMPLE_WORDS: usize = SAMPLES.div_ceil(64);

    /// Words per cell row. `dual.rs:484`'s `cells_x.div_ceil(64)`.
    const CELL_WORDS: usize = WORLD_CELLS.div_ceil(64);

    /// The registered granularities. 4³ is `M-377`'s optimum and is one `u64`.
    const GRANULARITIES: [usize; 4] = [4, 8, 16, 64];

    /// Roaring's container boundary: 4,096 set values per 2¹⁶ range.
    const ROARING_BOUNDARY: f64 = 4096.0 / 65536.0;

    /// Timed spans per arm, median taken.
    const REPS: usize = 9;

    /// Active cells a timed span must cover before it is believed.
    ///
    /// A span over one field's 19,000 active cells is tens of microseconds, which
    /// is enough for `Instant` but not enough to be quiet. The inner repeat count
    /// is derived per row so every span covers at least this many cells, and the
    /// divisor is `active_cells * inner` rather than `active_cells`.
    const SPAN_CELLS: u64 = 4_000_000;

    // ─── the sign rule, bench-local ────────────────────────────────────────────

    /// `isomesh::cube::is_inside`, which is private (`lib.rs:143` is `mod cube`,
    /// not `pub mod`), copied rather than published.
    ///
    /// **A sample of exactly zero is OUTSIDE.** Lengyel's choice, applied once
    /// everywhere in the crate; a cell is active exactly when its eight corner
    /// answers to this question are not all equal.
    #[inline]
    fn is_inside(value: Scalar) -> bool {
        value < 0.0
    }

    // ─── one field's world ─────────────────────────────────────────────────────

    /// The 128³ active-cell bitmap of one reference field, built once.
    ///
    /// Built **once per field and partitioned per granularity**, never re-sampled:
    /// the cells are the same 128³ cells however they are grouped, so
    /// `active_cells` must be identical across all four granularities. That is an
    /// exact integer equality and it is asserted, not hoped (`M-279`).
    struct World {
        /// Bit `x & 63` of word `CELL_WORDS * (y + WORLD_CELLS * z) + (x >> 6)`.
        cells: Vec<u64>,
        active_cells: u64,
        /// Cells the independent eight-corner scalar loop disagreed with the
        /// word-level `any & !all` fold about. Asserted zero.
        fold_mismatches: u64,
    }

    impl World {
        fn of<F: Sdf<Scalar = Scalar> + ReferenceField<Scalar = Scalar>>(field: &F) -> Self {
            let (lo, hi) = field.domain();
            let cell = (hi[0] - lo[0]) / (SAMPLES - 1) as Scalar;

            // The sample grid, on the field's own canonical domain — the same
            // convention `common::grid` uses, so "128³" here means what it means
            // everywhere else in this crate.
            let mut values = vec![0.0; SAMPLES * SAMPLES * SAMPLES];
            for z in 0..SAMPLES {
                for y in 0..SAMPLES {
                    for x in 0..SAMPLES {
                        values[x + SAMPLES * (y + SAMPLES * z)] = field.sample([
                            lo[0] + cell * x as Scalar,
                            lo[1] + cell * y as Scalar,
                            lo[2] + cell * z as Scalar,
                        ]);
                    }
                }
            }

            // One bit per SAMPLE, 64 to a word, along x. `dual.rs:359-381`.
            let rows = SAMPLES * SAMPLES;
            let mut inside = vec![0u64; SAMPLE_WORDS * rows];
            for row in 0..rows {
                let src = SAMPLES * row;
                for w in 0..SAMPLE_WORDS {
                    let base = w * 64;
                    let n = (SAMPLES - base).min(64);
                    let mut word = 0u64;
                    for k in 0..n {
                        word |= u64::from(is_inside(values[src + base + k])) << k;
                    }
                    inside[SAMPLE_WORDS * row + w] = word;
                }
            }
            let word = |w: usize, y: usize, z: usize| inside[SAMPLE_WORDS * (y + SAMPLES * z) + w];
            let shifted = |w: usize, y: usize, z: usize| {
                let lo = word(w, y, z);
                let hi = if w + 1 < SAMPLE_WORDS {
                    word(w + 1, y, z)
                } else {
                    0
                };
                (lo >> 1) | (hi << 63)
            };

            // The fused fold. `dual.rs:424-436`, plus `cell_mask` (`:445`) to
            // discard the bits of the last word that are samples and not cells.
            let mut cells = vec![0u64; CELL_WORDS * WORLD_CELLS * WORLD_CELLS];
            let mut active_cells = 0u64;
            for z in 0..WORLD_CELLS {
                for y in 0..WORLD_CELLS {
                    for w in 0..CELL_WORDS {
                        let mut any = 0u64;
                        let mut all = !0u64;
                        for dz in 0..2 {
                            for dy in 0..2 {
                                let a = word(w, y + dy, z + dz);
                                let b = shifted(w, y + dy, z + dz);
                                any |= a | b;
                                all &= a & b;
                            }
                        }
                        let remaining = WORLD_CELLS.saturating_sub(w * 64);
                        let mask = if remaining >= 64 {
                            !0u64
                        } else {
                            (1u64 << remaining) - 1
                        };
                        let active = (any & !all) & mask;
                        cells[CELL_WORDS * (y + WORLD_CELLS * z) + w] = active;
                        active_cells += u64::from(active.count_ones());
                    }
                }
            }

            // ── control: the fold against eight independent loads ─────────────
            //
            // The bit extraction above is four lines and trivially right; the
            // fold is the subtle part — the shifted high bit, the four rows, the
            // cell mask — and every density in this file is downstream of it. So
            // the eight corners are read again straight out of the `f64` grid and
            // the two answers are compared cell by cell.
            let mut fold_mismatches = 0u64;
            for z in 0..WORLD_CELLS {
                for y in 0..WORLD_CELLS {
                    for x in 0..WORLD_CELLS {
                        let mut inside_count = 0u32;
                        for corner in 0..8usize {
                            let v = values[(x + (corner & 1))
                                + SAMPLES
                                    * ((y + ((corner >> 1) & 1)) + SAMPLES * (z + (corner >> 2)))];
                            if is_inside(v) {
                                inside_count += 1;
                            }
                        }
                        let scalar_active = inside_count != 0 && inside_count != 8;
                        let word_active =
                            cells[CELL_WORDS * (y + WORLD_CELLS * z) + (x >> 6)] >> (x & 63) & 1
                                == 1;
                        if scalar_active != word_active {
                            fold_mismatches += 1;
                        }
                    }
                }
            }

            Self {
                cells,
                active_cells,
                fold_mismatches,
            }
        }

        #[inline]
        fn active_at(&self, x: usize, y: usize, z: usize) -> bool {
            self.cells[CELL_WORDS * (y + WORLD_CELLS * z) + (x >> 6)] >> (x & 63) & 1 == 1
        }
    }

    // ─── the array arm's index width ───────────────────────────────────────────

    /// A chunk-local cell index, at the narrowest width the chunk size admits.
    trait Idx: Copy {
        /// Distinct indices this width can hold.
        const CAPACITY: usize;
        fn from_local(i: usize) -> Self;
        fn widen(self) -> u64;
    }

    impl Idx for u8 {
        const CAPACITY: usize = 1 << 8;
        #[inline]
        fn from_local(i: usize) -> Self {
            i as Self
        }
        #[inline]
        fn widen(self) -> u64 {
            u64::from(self)
        }
    }

    impl Idx for u16 {
        const CAPACITY: usize = 1 << 16;
        #[inline]
        fn from_local(i: usize) -> Self {
            i as Self
        }
        #[inline]
        fn widen(self) -> u64 {
            u64::from(self)
        }
    }

    impl Idx for u32 {
        const CAPACITY: usize = 1 << 32;
        #[inline]
        fn from_local(i: usize) -> Self {
            i as Self
        }
        #[inline]
        fn widen(self) -> u64 {
            u64::from(self)
        }
    }

    // ─── the two walks ─────────────────────────────────────────────────────────

    /// The consumer, identical in both arms.
    ///
    /// A short serial dependency chain, so neither arm can be autovectorised into
    /// something that is no longer a walk, and **order-sensitive**, so equality of
    /// the two checksums proves the arms visited the same cells in the same order
    /// rather than merely the same number of them.
    ///
    /// # `^` was wrong here, and the harness's own control said so
    ///
    /// The first version was `acc.rotate_left(5) ^ idx`, and the non-zero
    /// assertion below **fired on the first run**: the checksum came out
    /// *exactly* zero at 4³ on `sphere`, on both arms. That is not a walk that
    /// visited nothing — it is XOR's linearity meeting a symmetric fixture. The
    /// final value of an XOR chain is `⊕ᵢ rot^{5(N−1−i)}(idxᵢ)`, so any two
    /// visits with the same index at positions congruent mod 64 **cancel**, and a
    /// sphere's 4³ chunks repeat the same 6-bit active pattern thousands of times
    /// by symmetry. A checksum that can cancel cannot certify a traversal, and one
    /// that cancels to zero is indistinguishable from one that ran on nothing.
    ///
    /// `wrapping_add` carries between bit positions, so it is not GF(2)-linear
    /// and no such pairwise cancellation exists.
    #[inline]
    fn consume(acc: u64, idx: u64) -> u64 {
        acc.rotate_left(7).wrapping_add(idx)
    }

    /// `dual.rs:495-497`'s walk: `trailing_zeros` to find, `x &= x - 1` to clear,
    /// which visits the set bits of a word in ascending order.
    fn walk_bitmap(chunks: &[Vec<u64>]) -> u64 {
        let mut acc = 0u64;
        for words in chunks {
            for (w, &word) in words.iter().enumerate() {
                let mut bits = word;
                while bits != 0 {
                    let idx = (w * 64) as u64 + u64::from(bits.trailing_zeros());
                    bits &= bits - 1;
                    acc = consume(acc, idx);
                }
            }
        }
        acc
    }

    fn walk_array<W: Idx>(chunks: &[Vec<W>]) -> u64 {
        let mut acc = 0u64;
        for indices in chunks {
            for &idx in indices {
                acc = consume(acc, idx.widen());
            }
        }
        acc
    }

    // ─── one arm ───────────────────────────────────────────────────────────────

    struct Arm {
        chunk_cells: usize,
        chunks_total: usize,
        chunks_occupied: usize,
        active_cells: u64,
        /// Active cells over the cells of the **occupied** chunks: the density
        /// Roaring's container rule is evaluated on.
        density: f64,
        /// Active cells over all 128³ world cells, for the reader who wants the
        /// other denominator.
        world_density: f64,
        index_bytes: usize,
        words_per_chunk: usize,
        words_scanned: u64,
        bitmap_bytes: u64,
        bitmap_payload_bytes: u64,
        bitmap_overhead_bytes: u64,
        array_bytes: u64,
        array_payload_bytes: u64,
        array_overhead_bytes: u64,
        byte_crossover_density: f64,
        /// Occupied chunks whose own density is above / below 6.25%: whether
        /// Roaring's rule is a decision *within* one field, not just across
        /// fields.
        chunks_above_boundary: usize,
        chunks_below_boundary: usize,
        chunks_array_smaller: usize,
        chunks_bitmap_smaller: usize,
        ns_per_walk_bitmap: f64,
        ns_per_walk_array: f64,
        cycles_per_active_bitmap: f64,
        cycles_per_active_array: f64,
        instructions_per_active_bitmap: f64,
        instructions_per_active_array: f64,
        /// The same walk over one contiguous ramp of the same length: the cheapest
        /// this consumer can be driven at all.
        ///
        /// Both arms pay `consume`'s serial latency once per active cell, so an arm
        /// sitting at this figure is **saturated** -- its own walk cost is hidden
        /// beneath the chain and its column is a ceiling rather than a reading.
        /// This is the column that says which arm that is.
        walk_floor_cycles: f64,
        walk_floor_ns: f64,
        ghz: f64,
        walk_checksum: u64,
        winner: &'static str,
    }

    /// Build both representations of one field at one granularity, and time both
    /// walks.
    fn measure<W: Idx>(world: &World, chunk_cells: usize, probe: &mut Probe) -> Arm {
        let cells_per_chunk = chunk_cells * chunk_cells * chunk_cells;
        assert!(
            cells_per_chunk <= W::CAPACITY,
            "{cells_per_chunk} chunk-local indices do not fit the chosen array width"
        );
        assert_eq!(
            cells_per_chunk % 64,
            0,
            "a {chunk_cells}³ chunk is {cells_per_chunk} cells, which is not a whole number of \
             words -- the bitmap arm's payload would be a partial word and its bytes a fiction"
        );
        let words_per_chunk = cells_per_chunk / 64;
        let per_axis = WORLD_CELLS / chunk_cells;
        let chunks_total = per_axis * per_axis * per_axis;

        // Both representations of every occupied chunk, built in the same chunk
        // order from the same cell enumeration, so the two walks are the same
        // traversal in two encodings.
        let mut bitmaps: Vec<Vec<u64>> = Vec::new();
        let mut arrays: Vec<Vec<W>> = Vec::new();
        let mut chunks_above_boundary = 0usize;
        let mut chunks_below_boundary = 0usize;
        let mut chunks_array_smaller = 0usize;
        let mut chunks_bitmap_smaller = 0usize;
        let mut active_cells = 0u64;

        let mut words = vec![0u64; words_per_chunk];
        let mut local: Vec<usize> = Vec::with_capacity(cells_per_chunk);
        for cz in 0..per_axis {
            for cy in 0..per_axis {
                for cx in 0..per_axis {
                    words.fill(0);
                    local.clear();
                    // Chunk-local index is `lx + c*(ly + c*lz)`, so at 4³ one
                    // `u64` holds a whole 4×4×4 block. Ascending index order in
                    // both arms.
                    for lz in 0..chunk_cells {
                        for ly in 0..chunk_cells {
                            for lx in 0..chunk_cells {
                                if world.active_at(
                                    cx * chunk_cells + lx,
                                    cy * chunk_cells + ly,
                                    cz * chunk_cells + lz,
                                ) {
                                    let i = lx + chunk_cells * (ly + chunk_cells * lz);
                                    words[i >> 6] |= 1u64 << (i & 63);
                                    local.push(i);
                                }
                            }
                        }
                    }
                    if local.is_empty() {
                        continue;
                    }
                    active_cells += local.len() as u64;
                    // `with_capacity` exactly, then fill: capacity equals length,
                    // so `array_bytes` is the footprint and not a growth curve.
                    // That is also what a real implementation does -- count, then
                    // scatter, which is `P-112`'s three phases.
                    let mut indices: Vec<W> = Vec::with_capacity(local.len());
                    indices.extend(local.iter().map(|&i| W::from_local(i)));
                    let bitmap_bytes = size_of::<Vec<u64>>() + words_per_chunk * 8;
                    let array_bytes = size_of::<Vec<W>>() + indices.len() * size_of::<W>();
                    if array_bytes < bitmap_bytes {
                        chunks_array_smaller += 1;
                    } else if bitmap_bytes < array_bytes {
                        chunks_bitmap_smaller += 1;
                    }
                    if local.len() as f64 / cells_per_chunk as f64 > ROARING_BOUNDARY {
                        chunks_above_boundary += 1;
                    } else {
                        chunks_below_boundary += 1;
                    }
                    bitmaps.push(words.clone());
                    arrays.push(indices);
                }
            }
        }

        let chunks_occupied = bitmaps.len();
        assert!(
            chunks_occupied > 0,
            "no occupied chunk at {chunk_cells}³, so there is nothing to represent"
        );
        assert_eq!(
            active_cells, world.active_cells,
            "partitioning at {chunk_cells}³ changed the active-cell count, which is a property of \
             the 128³ cell set and cannot depend on how it is grouped"
        );

        let bitmap_payload_bytes = (chunks_occupied * words_per_chunk * 8) as u64;
        let bitmap_overhead_bytes = (chunks_occupied * size_of::<Vec<u64>>()) as u64;
        let array_payload_bytes = active_cells * size_of::<W>() as u64;
        let array_overhead_bytes = (chunks_occupied * size_of::<Vec<W>>()) as u64;

        // ── the two walks ────────────────────────────────────────────────────
        let inner = (SPAN_CELLS / active_cells).max(1) as usize;
        let per_span = active_cells * inner as u64;

        let bitmap_span = || {
            let mut acc = 0u64;
            for _ in 0..inner {
                acc = walk_bitmap(&bitmaps);
            }
            acc
        };
        let array_span = || {
            let mut acc = 0u64;
            for _ in 0..inner {
                acc = walk_array(&arrays);
            }
            acc
        };

        // Warm-up, and the checksum comes from here rather than from a timed
        // span: it is the same walk either way, and taking it outside the timer
        // keeps the timed spans identical to each other.
        let bit_acc = black_box(bitmap_span());
        let arr_acc = black_box(array_span());
        assert_eq!(
            bit_acc, arr_acc,
            "the two arms disagree at {chunk_cells}³: the checksum is order-sensitive, so this is \
             either a different cell set or a different visitation order, and either makes the \
             walk comparison meaningless"
        );
        assert_ne!(
            bit_acc, 0,
            "the walk checksum is zero at {chunk_cells}³, so nothing was walked"
        );

        let mut bitmap_ns = Vec::with_capacity(REPS);
        let mut array_ns = Vec::with_capacity(REPS);
        for _ in 0..REPS {
            let t = Instant::now();
            black_box(bitmap_span());
            bitmap_ns.push(t.elapsed().as_secs_f64() * 1e9 / per_span as f64);
            let t = Instant::now();
            black_box(array_span());
            array_ns.push(t.elapsed().as_secs_f64() * 1e9 / per_span as f64);
        }
        bitmap_ns.sort_by(f64::total_cmp);
        array_ns.sort_by(f64::total_cmp);

        // ── the counted spans ────────────────────────────────────────────────
        //
        // `ghz` comes from the same spans as the cycles it divides, so a governor
        // move cannot put a clock and a cycle count from different moments on one
        // row (`M-281`).
        probe.reset_and_enable();
        let t = Instant::now();
        black_box(bitmap_span());
        let bitmap_counted_ns = t.elapsed().as_secs_f64() * 1e9;
        probe.disable();
        let bitmap_counts = probe.read();

        probe.reset_and_enable();
        let t = Instant::now();
        black_box(array_span());
        let array_counted_ns = t.elapsed().as_secs_f64() * 1e9;
        probe.disable();
        let array_counts = probe.read();

        // ── the instrument's own resolution ───────────────────────────────────
        //
        // `consume` is a serial rotate-and-add and both arms pay it once per
        // active cell, so neither arm can read below the chain's latency however
        // cheap its own walk is. An arm sitting **at** that limit is reporting a
        // ceiling rather than a reading, and a harness that does not say which arm
        // that is has published two numbers of different kinds under one heading.
        //
        // Measured, in the identical loop: `walk_array` over ONE contiguous ramp
        // of the same length -- no chunk boundaries, no cache pressure, perfectly
        // sequential. That is the cheapest this consumer can be driven at all, so
        // it is the resolution limit of both columns.
        //
        // **A first attempt used `(0..per_span).fold(0u64, consume)` and was
        // wrong**: it read 3.001 cycles per step while the real array arm read
        // 2.000, i.e. the "floor" was a worse loop than the thing it was meant to
        // bound. A resolution limit has to be measured in the same loop shape or
        // it is just a third arm.
        let ramp: Vec<W> = (0..active_cells as usize)
            .map(|i| W::from_local(i % cells_per_chunk))
            .collect();
        let ramp_chunks = [ramp];
        probe.reset_and_enable();
        let t = Instant::now();
        for _ in 0..inner {
            black_box(walk_array(&ramp_chunks));
        }
        let floor_ns = t.elapsed().as_secs_f64() * 1e9;
        probe.disable();
        let floor_counts = probe.read();

        for (which, counts) in [
            ("bitmap", &bitmap_counts),
            ("array", &array_counts),
            ("floor", &floor_counts),
        ] {
            assert!(
                counts.worst_ratio() >= MIN_TIME_RATIO,
                "a counter was multiplexed at ratio {:.4} on the {which} arm at {chunk_cells}³, \
                 so its cycles are an extrapolation and `ghz` would not be a reading",
                counts.worst_ratio()
            );
        }

        let ghz = (bitmap_counts.cycles.count + array_counts.cycles.count) as f64
            / (bitmap_counted_ns + array_counted_ns);

        let density = active_cells as f64 / (chunks_occupied * cells_per_chunk) as f64;
        let bitmap_bytes = bitmap_payload_bytes + bitmap_overhead_bytes;
        let array_bytes = array_payload_bytes + array_overhead_bytes;
        let cycles_per_active_bitmap = bitmap_counts.cycles.count as f64 / per_span as f64;
        let cycles_per_active_array = array_counts.cycles.count as f64 / per_span as f64;

        // The `winner` is a two-sided verdict, because C2 asks for both: bytes
        // and the walk. A representation that wins one and loses the other is
        // `split`, which is data rather than a tie.
        let array_wins_bytes = array_bytes < bitmap_bytes;
        let array_wins_walk = cycles_per_active_array < cycles_per_active_bitmap;
        let winner = match (array_wins_bytes, array_wins_walk) {
            (true, true) => "array",
            (false, false) => "bitmap",
            _ => "split",
        };

        Arm {
            chunk_cells,
            chunks_total,
            chunks_occupied,
            active_cells,
            density,
            world_density: active_cells as f64 / (WORLD_CELLS * WORLD_CELLS * WORLD_CELLS) as f64,
            index_bytes: size_of::<W>(),
            words_per_chunk,
            words_scanned: (chunks_occupied * words_per_chunk) as u64,
            bitmap_bytes,
            bitmap_payload_bytes,
            bitmap_overhead_bytes,
            array_bytes,
            array_payload_bytes,
            array_overhead_bytes,
            byte_crossover_density: 1.0 / (8.0 * size_of::<W>() as f64),
            chunks_above_boundary,
            chunks_below_boundary,
            chunks_array_smaller,
            chunks_bitmap_smaller,
            ns_per_walk_bitmap: bitmap_ns[REPS / 2],
            ns_per_walk_array: array_ns[REPS / 2],
            cycles_per_active_bitmap,
            cycles_per_active_array,
            instructions_per_active_bitmap: bitmap_counts.instructions.count as f64
                / per_span as f64,
            instructions_per_active_array: array_counts.instructions.count as f64 / per_span as f64,
            walk_floor_cycles: floor_counts.cycles.count as f64 / per_span as f64,
            walk_floor_ns: floor_ns / per_span as f64,
            ghz,
            walk_checksum: bit_acc,
            winner,
        }
    }

    /// The width dispatch. Narrowest type that holds a chunk-local index.
    fn measure_at(world: &World, chunk_cells: usize, probe: &mut Probe) -> Arm {
        let cells = chunk_cells * chunk_cells * chunk_cells;
        if cells <= <u8 as Idx>::CAPACITY {
            measure::<u8>(world, chunk_cells, probe)
        } else if cells <= <u16 as Idx>::CAPACITY {
            measure::<u16>(world, chunk_cells, probe)
        } else {
            measure::<u32>(world, chunk_cells, probe)
        }
    }

    /// C1's two conditions at one granularity, over the eight reference fields.
    struct Gate {
        chunk_cells: usize,
        below: usize,
        above: usize,
        distinct_winners: usize,
        straddles: bool,
        varies: bool,
    }

    pub(crate) fn run(run: &mut Run) {
        let mut probe = Probe::open();

        // ── the fixture ──────────────────────────────────────────────────────
        let mut fields: Vec<(&'static str, Vec<Arm>)> = Vec::new();
        for_each_reference_field!(Scalar, |name, field| {
            let world = World::of(&field);
            assert_eq!(
                world.fold_mismatches, 0,
                "{name}: the word-level `any & !all` fold disagrees with the eight-corner scalar \
                 loop on {} cells, so every density below is computed from a broken predicate",
                world.fold_mismatches
            );
            // VACUITY CONTROL, first half: a field with no active cell has no
            // density and no representation to choose between.
            assert!(
                world.active_cells > 0,
                "{name}: no active cell at {WORLD_CELLS}³, so `density` would be zero and both \
                 arms would be empty"
            );
            let arms: Vec<Arm> = GRANULARITIES
                .iter()
                .map(|&c| measure_at(&world, c, &mut probe))
                .collect();
            fields.push((name, arms));
        });

        // ── VACUITY CONTROL, second half: the density span ───────────────────
        //
        // The whole row turns on 6.25%, and a fixture whose densities all sit on
        // one side of it cannot test the boundary at all. A decade is the
        // registered bar.
        let all_densities: Vec<f64> = fields
            .iter()
            .flat_map(|(_, arms)| arms.iter().map(|a| a.density))
            .collect();
        let min_density = all_densities.iter().copied().fold(f64::INFINITY, f64::min);
        let max_density = all_densities
            .iter()
            .copied()
            .fold(f64::NEG_INFINITY, f64::max);
        let density_span = max_density / min_density;
        assert!(
            density_span >= 10.0,
            "VOID: measured density spans only {density_span:.4}x across the fixture \
             ({min_density:.6} to {max_density:.6}), which is under one decade -- Roaring's 6.25% \
             boundary is then untested by construction and no verdict here means anything"
        );

        // ── C1: a decision, or a constant? ───────────────────────────────────
        //
        // Two registered conditions, and both are evaluated over the eight fields
        // at one granularity, because "across the eight reference fields" is what
        // the clause says and "winning on all eight fields" is what its falsifier
        // says. The per-granularity table is printed so a stricter reading -- all
        // four granularities rather than one -- can be applied by the reader
        // without re-running anything.
        let gates: Vec<Gate> = GRANULARITIES
            .iter()
            .enumerate()
            .map(|(g, &c)| {
                let below = fields
                    .iter()
                    .filter(|(_, arms)| arms[g].density < ROARING_BOUNDARY)
                    .count();
                let above = fields
                    .iter()
                    .filter(|(_, arms)| arms[g].density > ROARING_BOUNDARY)
                    .count();
                let mut winners: Vec<&str> =
                    fields.iter().map(|(_, arms)| arms[g].winner).collect();
                winners.sort_unstable();
                winners.dedup();
                Gate {
                    chunk_cells: c,
                    below,
                    above,
                    distinct_winners: winners.len(),
                    straddles: below >= 1 && above >= 1,
                    varies: winners.len() > 1,
                }
            })
            .collect();
        let c1 = gates.iter().any(|g| g.straddles && g.varies);
        let granularities_straddling = gates.iter().filter(|g| g.straddles).count();
        let granularities_winner_varies = gates.iter().filter(|g| g.varies).count();

        // ── C2: below the threshold, does the array win on both? ─────────────
        let mut rows_below = 0usize;
        let mut c2_failures: Vec<String> = Vec::new();
        for (name, arms) in &fields {
            for a in arms {
                if a.density >= ROARING_BOUNDARY {
                    continue;
                }
                rows_below += 1;
                if a.array_bytes >= a.bitmap_bytes {
                    c2_failures.push(format!(
                        "{name}@{}³ bytes {} >= {}",
                        a.chunk_cells, a.array_bytes, a.bitmap_bytes
                    ));
                }
                if a.cycles_per_active_array >= a.cycles_per_active_bitmap {
                    c2_failures.push(format!(
                        "{name}@{}³ walk {:.4} >= {:.4} cyc/active",
                        a.chunk_cells, a.cycles_per_active_array, a.cycles_per_active_bitmap
                    ));
                }
            }
        }
        // A comparison with no rows to compare is VACUOUS, not held. `M-44`.
        let c2 = if rows_below == 0 {
            "VACUOUS".to_string()
        } else {
            c2_failures.is_empty().to_string()
        };

        // ── C3: at 4³ a chunk is one word, and the header dwarfs it ──────────
        let c3_index = GRANULARITIES
            .iter()
            .position(|&c| c == 4)
            .expect("4³ is in the registered fixture");
        let mut c3 = true;
        for (_, arms) in &fields {
            let a = &arms[c3_index];
            assert_eq!(
                a.words_per_chunk, 1,
                "a 4³ chunk must be exactly one word, and this one is {} -- C3's premise is \
                 arithmetic and a disagreement means the harness partitioned something else",
                a.words_per_chunk
            );
            let overhead = a.bitmap_overhead_bytes as f64 / a.active_cells as f64;
            let payload = a.bitmap_payload_bytes as f64 / a.active_cells as f64;
            if overhead <= payload {
                c3 = false;
            }
        }

        // ── report ───────────────────────────────────────────────────────────
        println!(
            "{:>14} {:>6} {:>8} {:>9} {:>9} {:>9} {:>10} {:>10} {:>9} {:>9} {:>9} {:>7}",
            "field",
            "chunk",
            "occupied",
            "active",
            "density",
            "world_d",
            "bmp_bytes",
            "arr_bytes",
            "cyc/act_b",
            "cyc/act_a",
            "cyc/floor",
            "winner"
        );
        for (name, arms) in &fields {
            for a in arms {
                println!(
                    "{:>14} {:>6} {:>8} {:>9} {:>9.5} {:>9.5} {:>10} {:>10} {:>9.3} {:>9.3} \
                     {:>9.3} {:>7}",
                    name,
                    a.chunk_cells,
                    a.chunks_occupied,
                    a.active_cells,
                    a.density,
                    a.world_density,
                    a.bitmap_bytes,
                    a.array_bytes,
                    a.cycles_per_active_bitmap,
                    a.cycles_per_active_array,
                    a.walk_floor_cycles,
                    a.winner
                );
            }
        }
        println!();
        for g in &gates {
            println!(
                "  {:>3}³: {} field(s) below 6.25%, {} above, {} distinct winner(s) -> straddles \
                 {}, varies {}",
                g.chunk_cells, g.below, g.above, g.distinct_winners, g.straddles, g.varies
            );
        }
        println!(
            "\ndensity span across the fixture: {min_density:.6} to {max_density:.6} = \
             {density_span:.2}x"
        );
        println!(
            "C1 density straddles 6.25% AND the winner varies over the eight fields at some \
             granularity: {granularities_straddling}/4 straddle, {granularities_winner_varies}/4 \
             vary -> {}",
            if c1 { "HELD" } else { "FALSIFIED" }
        );
        println!(
            "C2 below the threshold the array wins on bytes AND walk: {rows_below} row(s) below, \
             {} failure(s) -> {c2}",
            c2_failures.len()
        );
        for f in &c2_failures {
            println!("    {f}");
        }
        // Which arm's walk figure is a ceiling rather than a reading, and whether
        // that threatens the clause it decides.
        let saturated_array = fields
            .iter()
            .flat_map(|(_, arms)| arms.iter())
            .filter(|a| a.cycles_per_active_array <= a.walk_floor_cycles * 1.05)
            .count();
        let saturated_bitmap = fields
            .iter()
            .flat_map(|(_, arms)| arms.iter())
            .filter(|a| a.cycles_per_active_bitmap <= a.walk_floor_cycles * 1.05)
            .count();
        println!(
            "    walk saturation at the instrument floor: array {saturated_array}/32 rows, bitmap \
             {saturated_bitmap}/32 -- a saturated arm reports a CEILING, so `array < bitmap` on a \
             saturated array arm holds a fortiori and `bitmap < array` on a saturated bitmap arm \
             would be void"
        );
        // Roaring's boundary hit exactly, which is worth naming rather than
        // rounding: neither below nor above, so it counts for neither side of C1.
        for (name, arms) in &fields {
            for a in arms {
                if (a.density - ROARING_BOUNDARY).abs() < 1e-12 {
                    println!(
                        "    {name} at {}³ sits EXACTLY on Roaring's 6.25%: {} active over {} \
                         occupied chunks of {} cells is exactly 1/16, so it is neither below nor \
                         above and counts for neither side of C1",
                        a.chunk_cells,
                        a.active_cells,
                        a.chunks_occupied,
                        a.words_per_chunk * 64
                    );
                }
            }
        }
        println!(
            "C3 at 4³ the bitmap's per-chunk overhead ({} B) exceeds its payload ({} B) per active \
             cell: -> {}",
            size_of::<Vec<u64>>(),
            8,
            if c3 { "HELD" } else { "FALSIFIED" }
        );

        for (name, arms) in &fields {
            for (g, a) in arms.iter().enumerate() {
                let gate = &gates[g];
                run.record(&[
                    ("field", (*name).to_string()),
                    ("chunk_cells", a.chunk_cells.to_string()),
                    ("active_cells", a.active_cells.to_string()),
                    ("density", format!("{:.8}", a.density)),
                    ("bitmap_bytes", a.bitmap_bytes.to_string()),
                    ("array_bytes", a.array_bytes.to_string()),
                    (
                        "bytes_per_active_cell_bitmap",
                        format!("{:.6}", a.bitmap_bytes as f64 / a.active_cells as f64),
                    ),
                    (
                        "bytes_per_active_cell_array",
                        format!("{:.6}", a.array_bytes as f64 / a.active_cells as f64),
                    ),
                    ("ns_per_walk_bitmap", format!("{:.6}", a.ns_per_walk_bitmap)),
                    ("ns_per_walk_array", format!("{:.6}", a.ns_per_walk_array)),
                    ("winner", a.winner.to_string()),
                    ("c1_holds", c1.to_string()),
                    ("c2_holds", c2.clone()),
                    ("c3_holds", c3.to_string()),
                    // ── extras (M-273) ───────────────────────────────────────
                    ("world_cells", WORLD_CELLS.to_string()),
                    ("world_density", format!("{:.8}", a.world_density)),
                    ("chunks_total", a.chunks_total.to_string()),
                    ("chunks_occupied", a.chunks_occupied.to_string()),
                    ("cells_per_chunk", (a.words_per_chunk * 64).to_string()),
                    ("words_per_chunk", a.words_per_chunk.to_string()),
                    ("words_scanned", a.words_scanned.to_string()),
                    ("index_bytes", a.index_bytes.to_string()),
                    (
                        "byte_crossover_density",
                        format!("{:.6}", a.byte_crossover_density),
                    ),
                    ("roaring_boundary", format!("{ROARING_BOUNDARY:.6}")),
                    ("bitmap_payload_bytes", a.bitmap_payload_bytes.to_string()),
                    ("bitmap_overhead_bytes", a.bitmap_overhead_bytes.to_string()),
                    ("array_payload_bytes", a.array_payload_bytes.to_string()),
                    ("array_overhead_bytes", a.array_overhead_bytes.to_string()),
                    (
                        "bitmap_payload_bytes_per_active_cell",
                        format!(
                            "{:.6}",
                            a.bitmap_payload_bytes as f64 / a.active_cells as f64
                        ),
                    ),
                    (
                        "bitmap_overhead_bytes_per_active_cell",
                        format!(
                            "{:.6}",
                            a.bitmap_overhead_bytes as f64 / a.active_cells as f64
                        ),
                    ),
                    ("chunks_above_boundary", a.chunks_above_boundary.to_string()),
                    ("chunks_below_boundary", a.chunks_below_boundary.to_string()),
                    ("chunks_array_smaller", a.chunks_array_smaller.to_string()),
                    ("chunks_bitmap_smaller", a.chunks_bitmap_smaller.to_string()),
                    (
                        "cycles_per_active_cell_bitmap",
                        format!("{:.6}", a.cycles_per_active_bitmap),
                    ),
                    (
                        "cycles_per_active_cell_array",
                        format!("{:.6}", a.cycles_per_active_array),
                    ),
                    (
                        "instructions_per_active_cell_bitmap",
                        format!("{:.6}", a.instructions_per_active_bitmap),
                    ),
                    (
                        "instructions_per_active_cell_array",
                        format!("{:.6}", a.instructions_per_active_array),
                    ),
                    (
                        "walk_cycle_ratio",
                        format!(
                            "{:.6}",
                            a.cycles_per_active_bitmap / a.cycles_per_active_array
                        ),
                    ),
                    (
                        "walk_cycle_difference",
                        format!(
                            "{:.6}",
                            a.cycles_per_active_bitmap - a.cycles_per_active_array
                        ),
                    ),
                    (
                        "walk_floor_cycles_per_active_cell",
                        format!("{:.6}", a.walk_floor_cycles),
                    ),
                    (
                        "walk_floor_ns_per_active_cell",
                        format!("{:.6}", a.walk_floor_ns),
                    ),
                    // Which arms are saturated at the chain, i.e. reporting a
                    // ceiling rather than a reading. C2's walk inequality still
                    // holds a fortiori when the WINNER is the saturated arm --
                    // its true cost is bounded above by a figure that already
                    // loses -- and would be void if the LOSER were the saturated
                    // one. This column is how a reader tells those apart.
                    (
                        "bitmap_walk_saturated",
                        (a.cycles_per_active_bitmap <= a.walk_floor_cycles * 1.05).to_string(),
                    ),
                    (
                        "array_walk_saturated",
                        (a.cycles_per_active_array <= a.walk_floor_cycles * 1.05).to_string(),
                    ),
                    ("ghz", format!("{:.4}", a.ghz)),
                    ("walk_checksum", a.walk_checksum.to_string()),
                    ("fold_mismatches", "0".to_string()),
                    ("density_span_over_fixture", format!("{density_span:.6}")),
                    // C1's reachable share, over the eight reference fields at
                    // THIS row's granularity: the two counts either side of 6.25%
                    // must sum to 8 unless a field sits exactly on it, and
                    // `distinct_winners` is the falsifier's own population -- 1
                    // means one representation won on all eight.
                    ("fields_below_boundary", gate.below.to_string()),
                    ("fields_above_boundary", gate.above.to_string()),
                    ("distinct_winners", gate.distinct_winners.to_string()),
                    ("c1_straddles_here", gate.straddles.to_string()),
                    ("c1_winner_varies_here", gate.varies.to_string()),
                    (
                        "granularities_straddling",
                        granularities_straddling.to_string(),
                    ),
                    (
                        "granularities_winner_varies",
                        granularities_winner_varies.to_string(),
                    ),
                    ("rows_below_boundary", rows_below.to_string()),
                    ("c2_failures", c2_failures.len().to_string()),
                ]);
            }
        }
    }
}

fn main() {
    if !std::env::args().any(|arg| arg == "--bench") {
        return;
    }
    let prereg = isomesh::experiment!("P-113");

    #[cfg(target_os = "linux")]
    common::experiment::run(prereg, experiment::run);

    #[cfg(not(target_os = "linux"))]
    {
        eprintln!(
            "{} scores C2's walk comparison on cycles per active cell and puts `ghz` on every row \
             it reports a nanosecond on, and both come from `perf_event_open` -- a Linux system \
             call with no equivalent this bench can reach.\n\
             Refusing rather than writing a zero into a column nobody measured (M-280, M-281).\n\
             Run it on Linux; `perf_event_paranoid = 2` is permissive enough and no root is needed.",
            prereg.id
        );
        std::process::exit(1);
    }
}
