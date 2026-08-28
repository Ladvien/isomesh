//! **P-103 — the case index computed 64 cells at a time from bit-sliced sign planes.**
//!
//! Ticket: R-103. Pre-registered before this harness existed.
//!
//! ```bash
//! cargo bench --bench experiment_p103
//! ```
//!
//! Writes `docs/experiments/p-103.csv`. **Linux only**, for `experiment_p12`'s
//! reason and the registration's own: C1 is a *share of extraction* and C2 is an
//! instruction rate, and the only instrument that can read either is
//! `perf_event_open`. On a governed CPU a nanosecond is not a unit (`✗24`,
//! `M-280`, `M-281`), so off Linux there is nothing to degrade to and the
//! harness refuses and exits 1 rather than record a fabricated zero.
//!
//! # What was missing
//!
//! FastLanes' `Scalar_T64` result — *"uses 64-bit scalar registers as quasi-SIMD
//! and beats naive Scalar up to 8x"*, and *"when incorporating FastLanes in
//! future systems, we recommend just using the Scalar code paths"* — needs a
//! quantity in the crate that is already boolean. There is exactly one: the
//! eight corner signs of a cell.
//!
//! **The comparand is exact, and it is not a shared classifier.**
//! `marching_cubes/mod.rs:258-268` builds the index inline with `case |= 1 << c`
//! guarded by `is_inside` (`cube.rs:171`, `value < 0`, so an exact zero is
//! *outside*), and there are **six independent copies of that loop**:
//!
//! | Site | Index type |
//! |---|---|
//! | `marching_cubes/mod.rs:266` | `u8` |
//! | `manifold_dual_contouring.rs:211` | `u8` |
//! | `property/extraction.rs:439` | `u8` |
//! | `property/extraction.rs:585` | `u8` |
//! | `marching_cubes/ambiguity/tests.rs:264` | `u8` |
//! | `transvoxel/cell.rs:231` | `u16` |
//!
//! There is **no shared classifier anywhere in the crate**. So this row proposes
//! the crate's *first* one rather than modifying an existing one, and it is
//! bench-local: `crates/isomesh/src/**` is read-only for the whole of Phase 25.
//!
//! Nobody had measured what a bit-sliced classifier costs, and — the harder gap
//! — nobody had measured what the *byte* classifier costs as a fraction of an
//! extraction. `P-121` has since built that decomposition, and its numbers are
//! this row's **prior, not its measurement** (`M-281`: a share read off another
//! binary's run is not a share of this one).
//!
//! # C1 is the share gate and it is checked first
//!
//! `✗51`'s rule applied at the top of the row rather than remembered at the
//! bottom: a 2× on a stage under 10% has a ceiling of `1/(1 − 0.10/2) = 1.05×`
//! and is not worth building.
//!
//! `P-121`'s prior for `cycles_classify / cycles_total` on `marching_cubes` at
//! 65³ `f32` is **0.6459** (sphere), 0.5382 (torus), 0.5042 (`box_exact`),
//! 0.4123 (`csg_difference`), 0.5435 (`thin_plate`), 0.1434 (gyroid),
//! **0.0029–0.018** (`fbm_terrain`) and 0.0669 (`noise_cavity`). So the 10% bar
//! was expected to clear comfortably on six of eight fields and to fail on
//! `fbm_terrain`, where field evaluation is up to 94% of extraction. **That
//! expectation is not imported.** `classify_share` is measured here, per row,
//! and `c1_holds` is scored per row.
//!
//! Two things pin the share so it means something:
//!
//! - **The numerator is `P-121`'s `classify` stage, not something adjacent.**
//!   [`Classify::byte_compact`] is the eight-sign loop plus the compaction
//!   `case != 0 && case != u8::MAX`, which is exactly the boundary `P-121` drew
//!   and exactly `mod.rs:258-268` minus the `[R; 8]` corner-value gather — that
//!   gather feeds `edge_position` and belongs to interpolation. Keeping the same
//!   boundary is what makes the two rows' shares comparable at all.
//! - **The denominator is the shipped extractor**, `MarchingCubes::extract` on
//!   the same grid in the same run, sampling included. Nothing is mirrored, so
//!   there is no mirror-agreement debt on the share. It is also the
//!   *conservative* choice: `P-121` measured the shipped `MarchingCubes` at
//!   **15–24% more** than its bit-identical classify-then-compact mirror, so a
//!   share taken against the shipped total is that much **lower** than the same
//!   share taken against `P-121`'s. C1 is scored on the harder denominator.
//!
//! One honest caveat rather than a hidden one: the numerator runs with its input
//! array resident and its output vector already at capacity, so it is the
//! classify stage *at its cheapest*, which biases `classify_share` down and C1
//! toward falsification. `P-121` established that isolation changes what a stage
//! costs — it measured `sample` at 5.19 and 6.74 cycles per sample in two
//! shapes — and this row does not pretend otherwise.
//!
//! `classify_share_instructions` is beside `classify_share` for `M-280`'s
//! reason: instructions are deterministic here and cycles are not. The
//! registered column is the cycle form, because `P-121`'s decomposition is in
//! cycles; both are reported and the printed summary says how many rows clear
//! the bar under each.
//!
//! # C2, and what "produce 64 case indices" is measured as
//!
//! The registration's sentence has two halves that pull in different
//! directions — *"eight sign **bit**-planes"* and *"a whole **u64** of case
//! indices falls out of eight shift-and-OR pairs"* — because a `u64` holds 64
//! one-bit lanes but only 8 byte-wide case indices. Both readings are built and
//! both are columns; C2 is scored on the one the registration names first.
//!
//! - **[`Classify::bitsliced_dense`] — the registered arm, 64 lanes.** One sign
//!   **bit** per sample packed 64 to a `u64` along `x`, which is
//!   `dual.rs:359-381`'s own structure carried onto the Marching Cubes sample
//!   layout. The eight corner planes of a 64-cell word are shifts of that one
//!   plane: corner `c` has offset `[c & 1, (c >> 1) & 1, (c >> 2) & 1]`, so the
//!   four `dx = 0` planes are word loads from four sample rows and the four
//!   `dx = 1` planes are `(lo >> 1) | (hi << 63)` — `dual.rs:395`'s
//!   `inside_word_shifted`, verbatim. **One sign test per sample, not eight per
//!   cell**, and that is the mechanism: the byte path evaluates `is_inside`
//!   `8·(n−1)³` times where this evaluates it `n³` times.
//!
//!   The 64 case indices then come out of [`cases_word`], which is eight
//!   iterations each ending in a shift-and-OR pair, run eight times per word:
//!
//!   ```ignore
//!   let byte   = (plane >> (8 * group)) & 0xff;      // eight cells' sign bits
//!   let masked = (byte * REPLICATE) & LANE_BITS;     // bit j -> lane j, as 1<<j
//!   let ones   = ((masked + HALF) >> 7) & REPLICATE; // lane j -> 0 or 1
//!   out |= ones << c;                                // the shift-and-OR pair
//!   ```
//!
//!   The middle step is carry-free by construction and that is why it is used
//!   rather than the shorter magic-multiply spread: a lane holds at most `0x80`,
//!   `0x80 + 0x7f = 0xff`, so no lane can carry into its neighbour. The
//!   multiply-and-shift form (`b * 0x8040201008040201`) *does* alias — its
//!   partial products overlap and a carry out of bit 14 lands in the bit 15 slot
//!   the query reads — so it would have been a silent one-in-N corruption rather
//!   than a compile error. `cases_identical` would have caught it; it is stated
//!   here because the next reader should not have to rediscover it.
//!
//! - **[`Classify::swar8_dense`] — the other reading, 8 lanes.** One sign
//!   **byte** (`0` or `1`) per sample, read eight at a time as a `u64`. Then the
//!   `dx = 1` plane is a *byte*-offset load rather than a shift, and the whole
//!   assembly is literally `p0 | p1 << 1 | … | p7 << 7` — the registration's
//!   sentence with nothing between it and the machine. Reported as
//!   `instructions_per_cell_swar8` and `instruction_ratio_swar8`, and **not**
//!   what C2 is scored on: the registration says the bit-planes replace *"one
//!   array of bytes"*, and this arm is an array of bytes.
//!
//! - **[`Classify::bitsliced_active`] — what a shared classifier would actually
//!   return.** The same eight planes, folded to `any & !all` (`dual.rs:424`) for
//!   the active-cell word, then a `trailing_zeros` + `x &= x − 1` walk that
//!   materialises a case byte only for the ~1–2% of cells that produce
//!   something. Reported as `instructions_per_cell_bitsliced_active`. It is an
//!   extra column and not C2's arm, because it does not "produce 64 case
//!   indices" — it produces the answer a pipeline needs, which is a different
//!   and smaller claim.
//!
//! **The comparand stores what the arms store.** `mod.rs:262-268` keeps `case`
//! in a register and consumes it in the same loop iteration; a bit-sliced word
//! cannot. So [`Classify::byte_dense`] writes one dense case byte per cell, the
//! same deliverable the two dense arms write, and `instruction_ratio` is
//! deliverable-for-deliverable. `instructions_per_cell_byte_compact` is beside
//! it — the shipped shape, which branches and compacts instead of storing — so a
//! reader can see the comparand in both forms and neither choice can hide
//! anything.
//!
//! **Field sampling is outside every window.** The `values` array is filled once
//! per row before anything is counted, and every arm reads that one array.
//! `P-121` measured field evaluation at up to 94% of extraction on
//! `fbm_terrain`; counting it inside the arms would drive every instruction
//! ratio toward 1.00 and C2 would then hold or fail by dilution rather than by
//! measurement. Sampling *is* inside the C1 denominator, because the denominator
//! is a whole extraction.
//!
//! **Windows are siblings, never nested.** `R-121` paid for that discovery: Zen
//! 3 has six general-purpose counters and `Probe` opens six, so two nested
//! windows multiplex and `Probe::worst_ratio` refuses. There are six flat
//! windows per repetition and the arm order **rotates by repetition index**, so
//! no arm is permanently the one that inherits another's cache and predictor
//! state.
//!
//! # C3: the case index is an integer, and this is how far that is proved
//!
//! `M-31`'s fixture is **216 hashes** = 8 reference fields (`fields/mod.rs:212`)
//! × 9 algorithms (`golden.rs:122`) × 3 resolutions `{17, 25, 33}`
//! (`golden.rs:73`), gated by `golden_hashes_are_unchanged`
//! (`golden/tests.rs:59`) against `crates/isomesh/golden_hashes.json`. The gate
//! is proven able to fire: `P-61` moved 135 of the same 216.
//!
//! `src/**` is read-only, so the bit-sliced classifier cannot be driven through
//! the shipped extractors. C3 is therefore established in three layers, and the
//! registered `hashes_moved` is the strongest of them that is reachable
//! bench-locally:
//!
//! 1. **Mesh hashes against the committed fixture, on the 24 `marching_cubes`
//!    rows.** [`March`] is a bench-local Marching Cubes whose control stream
//!    comes from either classifier and whose payload code
//!    ([`March::emit_cell`]) is the *same instructions* either way. Its byte arm
//!    is required to reproduce the committed `golden_hashes.json` hash for every
//!    one of the 24 `(marching_cubes, field, samples)` rows — asserted, and
//!    `golden_hashes_reproduced` is the column. That is `M-279`'s rule (a new
//!    instrument's first job is to agree with the old one where they overlap)
//!    and it is the licence for reading anything off the mirror. `hashes_moved`
//!    is then how many of those 24 committed hashes the **bit-sliced** arm fails
//!    to reproduce. The grid is the fixture's own: `golden.rs:159-162` computes
//!    `cell_size = (hi[0] − lo[0]) / (samples − 1)` from `field.domain()`, which
//!    is `benches/common/mod.rs:43-45` character for character.
//! 2. **The case index itself, exhaustively, on every golden grid.** For all 8
//!    fields × `{17, 25, 33}` × `{f32, f64}`, every one of the
//!    `golden_cells_checked` cells has its byte and bit-sliced case index
//!    compared; `golden_case_mismatches` must be 0. This is what covers the
//!    other **192** rows: eight of the nine algorithms are not mirrored here,
//!    but an algorithm that reads the same eight-sign integer on every cell of
//!    the same grid cannot produce a different mesh.
//! 3. **The sign predicate, per sample.** `sign_bits_identical` /
//!    `samples_checked` compares every bit of the packed plane against
//!    `is_inside(values[s])` directly. This is the layer that covers the
//!    algorithms which never form a case index at all — `surface_nets`,
//!    `dual_contouring`, `greedy_quads` and the rest read the sign predicate,
//!    and if the plane reproduces it everywhere then so do they.
//!
//! # VACUITY CONTROL, asserted rather than recorded
//!
//! `cases_identical` must equal `cells` on every row, and
//! `corrupt_control_mismatches` — the same comparison with **one sign plane
//! deliberately flipped** — must be non-zero, because a comparator that cannot
//! report a mismatch has not reported a match. Plane 3 is inverted, which flips
//! bit 3 of every cell's case index, so the expected reading is `cells` exactly
//! and it is asserted as such rather than merely as "> 0".
//!
//! That control is deliberately coarse, so a second and tighter one is beside
//! it: `single_bit_control_mismatches` flips **one bit** of plane 3, the bit
//! belonging to the centre cell, and must be exactly **1**. The whole-plane
//! control proves the comparator can see *a* difference; the single-bit control
//! proves it is per-cell and not merely a length or a checksum.
//!
//! # The popcount this build does not have, and why no verdict here needs it
//!
//! There is no `.cargo/config.toml` and no `target-cpu` anywhere in the
//! repository, so the default `x86-64` baseline is in force,
//! `cfg!(target_feature = "popcnt")` is **false**, and `u64::count_ones()`
//! lowers to the ~12-instruction SWAR sequence rather than to the one
//! instruction this CPU actually has. Every published rank/select figure that
//! assumes a hardware popcount is describing something an order of magnitude
//! cheaper than what this build runs.
//!
//! **P-103 does not touch that quantity on either side.** Every arm here makes
//! **zero** `count_ones` calls per cell — the dense arms are shifts, masks, adds
//! and ORs; the active arm is `trailing_zeros` and `x &= x − 1`; the byte path
//! is comparisons and ORs. So `instruction_ratio`, `instruction_ratio_swar8` and
//! `classify_share` are all invariant to the popcount lowering, and no verdict
//! on this row could move under `-C target-cpu=native`. Measuring under such a
//! build would mean comparing across binaries, which `M-281` forbids;
//! `target_feature_popcnt`, `count_ones_calls_per_cell_byte` and
//! `count_ones_calls_per_cell_bitsliced` are the honest substitute and they are
//! columns.
//!
//! # SHARE
//!
//! Each clause's reachable share, as a column.
//!
//! - **C1's share is `classify_share` itself** — the clause *is* the share, and
//!   it gates the row. The bar is 0.10, per row, against
//!   `cycles_extract_shipped`. `classify_share_instructions` is beside it, and
//!   the run prints how many of the 32 rows clear the bar under each form before
//!   it writes anything. A row under the bar is a **falsified C1 on that row**,
//!   recorded as such — `✗51`'s precedent, and the cheapest useful output this
//!   row has.
//! - **C2's share is `instructions_per_cell_byte`** — the byte classifier's own
//!   instruction rate, which is the whole quantity a halving would halve. The
//!   bar is `instruction_ratio < 0.5`, scored on the **worst row**, which is the
//!   tight reading: a cheap field can only push a ratio toward 1.0, never below
//!   the maximum. `worst_instruction_ratio` is a column, and so is
//!   `instruction_ratio_rep_spread`, which *demonstrates* the determinism
//!   instead of asserting it. **C2 is the clause that carries the verdict**
//!   (`M-280`, `M-281`): `R-105` watched one binary's cycle ratio band drift
//!   from 0.984 to 1.035 across three runs while its instruction counts held to
//!   four figures. `cycle_ratio`, `cycle_ratio_rep_best` and
//!   `cycle_ratio_rep_worst` are reported and no clause consults them.
//! - **C3's share is an equality over an enumerated population**, three times
//!   over: 24 committed mesh hashes, `golden_cells_checked` case indices, and
//!   `samples_checked` sign bits. It moves no time and has no ceiling.
//! - **No time claim is made.** `ns_per_cell_byte` and `ghz` are on every row as
//!   provenance so a later reader can see what clock the counts were taken at,
//!   and no clause consults either.

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

    use isomesh::marching_cubes::MarchingCubes;
    use isomesh::marching_cubes::table::{
        CASES, CENTROID_BASE, EDGE_AXIS, EDGE_CORNERS, MAX_CENTROIDS, edge_offset, is_centroid,
        is_inside,
    };
    use isomesh::validate::mesh_hash;
    use isomesh::{MeshBuffer, Real, RuntimeShape3, Sdf};

    use crate::common::counters::{MIN_TIME_RATIO, Probe};

    // ─── the fixture ───────────────────────────────────────────────────────

    /// The registered resolutions, in samples per axis.
    const RESOLUTIONS: [u32; 2] = [33, 65];
    /// `golden.rs:73`'s resolutions, which C3 is established over.
    const GOLDEN_RESOLUTIONS: [u32; 3] = [17, 25, 33];
    /// The 216-hash fixture's eight reference fields (`fields/mod.rs:212`).
    const REFERENCE_FIELDS: usize = 8;
    /// Measured repetitions per row. The median per quantity carries the number;
    /// the per-repetition spread of the instruction ratio is a column, so the
    /// determinism C2 relies on is demonstrated rather than asserted.
    ///
    /// **Six on purpose**, because there are six sibling windows per repetition
    /// and the arm order rotates by repetition index: a multiple of six gives
    /// every arm each position in the rotation exactly once, so no arm's median
    /// is decided by having drawn the favourable slot one extra time.
    const REPS: usize = 6;
    /// Untimed passes of every arm before anything is counted, so the buffers
    /// are at final capacity and the pages are faulted in.
    const WARMUP: usize = 2;
    /// How long one counter window should last, in nanoseconds.
    ///
    /// `experiment_p121`'s figure. The quantities scored here are cycles and
    /// instructions — the two densest events the probe reads — so 30 ms holds
    /// tens of millions of each and the window needs no more than that.
    const TARGET_BATCH_NS: f64 = 30_000_000.0;
    /// Ceiling on the batch, so a cheap row cannot take minutes.
    const MAX_INNER: usize = 4096;
    /// C1's bar, from the registration.
    const CLASSIFY_SHARE_BAR: f64 = 0.10;
    /// C2's bar, from the registration: *under half*.
    const INSTRUCTION_RATIO_BAR: f64 = 0.5;

    // ─── private crate mechanisms, copied rather than made `pub` ────────────

    /// `cube::corner_offset`. Private, and `src/**` is read-only this phase.
    ///
    /// The three bits are `dx`, `dy`, `dz` in that order, which is what makes
    /// the eight corner planes four sample rows and one shift rather than eight
    /// independent gathers.
    #[inline]
    const fn corner_offset(corner: u8) -> [u32; 3] {
        [
            (corner & 1) as u32,
            ((corner >> 1) & 1) as u32,
            ((corner >> 2) & 1) as u32,
        ]
    }

    /// `cube::place`: the centred frame, spelled once in the crate and once
    /// here.
    #[inline]
    fn place<R: Real>(lo: R, hi: R, d: R) -> R {
        (lo + hi) * R::HALF + (hi - lo) * d
    }

    /// `vec3::length`. Left-to-right summation, which is what makes the mirror's
    /// normals bit-identical to the crate's rather than merely close.
    #[inline]
    fn length<R: Real>(a: [R; 3]) -> R {
        (a[0] * a[0] + a[1] * a[1] + a[2] * a[2]).sqrt()
    }

    /// `marching_cubes::unit_gradient`, over `vec3::scale`.
    #[inline]
    fn unit_gradient<R: Real, S: Sdf<Scalar = R>>(sdf: &S, position: [R; 3]) -> [R; 3] {
        let g = sdf.gradient(position);
        let inv = length(g).recip();
        [g[0] * inv, g[1] * inv, g[2] * inv]
    }

    /// `marching_cubes::corner_position`.
    #[inline]
    fn corner_position<R: Real>(
        base: [u32; 3],
        corner: u8,
        origin: [R; 3],
        cell_size: R,
    ) -> [R; 3] {
        let o = corner_offset(corner);
        [
            origin[0] + cell_size * R::from_f64(f64::from(base[0] + o[0])),
            origin[1] + cell_size * R::from_f64(f64::from(base[1] + o[1])),
            origin[2] + cell_size * R::from_f64(f64::from(base[2] + o[2])),
        ]
    }

    /// `marching_cubes::edge_position` with `crossing_refinement == 0`, which is
    /// `MarchingCubes::new`'s default: `refine_crossing` returns `d0` unchanged
    /// at zero steps, so there is no field evaluation on this path.
    #[inline]
    fn edge_position<R: Real>(
        base: [u32; 3],
        edge: u8,
        corner_value: &[R; 8],
        origin: [R; 3],
        cell_size: R,
    ) -> [R; 3] {
        let [lo_corner, hi_corner] = EDGE_CORNERS[edge as usize];
        let d = edge_offset(
            corner_value[lo_corner as usize],
            corner_value[hi_corner as usize],
        );
        let lo_pos = corner_position(base, lo_corner, origin, cell_size);
        let hi_pos = corner_position(base, hi_corner, origin, cell_size);
        [
            place(lo_pos[0], hi_pos[0], d),
            place(lo_pos[1], hi_pos[1], d),
            place(lo_pos[2], hi_pos[2], d),
        ]
    }

    // ─── the grid ──────────────────────────────────────────────────────────

    /// One cubic grid: `n` samples per axis, spanning `n − 1` cells.
    #[derive(Clone, Copy)]
    struct Grid<R: Real> {
        n: u32,
        origin: [R; 3],
        cell_size: R,
    }

    impl<R: Real> Grid<R> {
        #[inline]
        fn samples(self) -> usize {
            let n = self.n as usize;
            n * n * n
        }

        #[inline]
        fn cells(self) -> usize {
            self.n as usize - 1
        }

        #[inline]
        fn cell_count(self) -> usize {
            let c = self.cells();
            c * c * c
        }

        /// `RuntimeShape3::linearize`, which is the layout Marching Cubes
        /// samples into: no row padding, so the stride is the row itself.
        #[inline]
        fn sample_index(self, p: [u32; 3]) -> usize {
            let n = self.n as usize;
            p[0] as usize + n * (p[1] as usize + n * p[2] as usize)
        }

        /// The first sample of row `(y, z)`.
        #[inline]
        fn sample_row(self, y: usize, z: usize) -> usize {
            let n = self.n as usize;
            n * (y + n * z)
        }

        /// Words per **sample** row in the bit plane. `dual.rs:363`'s
        /// `bit_row = size[0].div_ceil(64)` — and the cell row is one shorter,
        /// which is `cell_words`.
        #[inline]
        fn bit_row(self) -> usize {
            (self.n as usize).div_ceil(64)
        }

        /// Words per **cell** row. `dual.rs:484`'s `cell_words`.
        #[inline]
        fn cell_words(self) -> usize {
            self.cells().div_ceil(64)
        }
    }

    // ─── the byte path: the comparand, in one place ─────────────────────────

    /// `mod.rs:259-268`'s eight-sign case index for one cell.
    ///
    /// The single sign loop in this file, called by the byte arms and by
    /// [`March`]'s byte control pass, so no arm has a cheaper way to classify a
    /// cell than another. `is_inside` is `value < 0`, so an exact zero is
    /// **outside**.
    #[inline]
    fn case_of<R: Real>(values: &[R], g: Grid<R>, base: [u32; 3]) -> u8 {
        let mut case = 0u8;
        for c in 0..8u8 {
            let o = corner_offset(c);
            let v = values[g.sample_index([base[0] + o[0], base[1] + o[1], base[2] + o[2]])];
            if is_inside(v) {
                case |= 1 << c;
            }
        }
        case
    }

    /// `mod.rs:259-268` **with** the `[R; 8]` corner-value gather, which
    /// interpolation needs and classification does not.
    #[inline]
    fn case_and_corners<R: Real>(values: &[R], g: Grid<R>, base: [u32; 3]) -> (u8, [R; 8]) {
        let mut case = 0u8;
        let mut corner_value = [R::ZERO; 8];
        for (c, slot) in corner_value.iter_mut().enumerate() {
            let o = corner_offset(c as u8);
            let v = values[g.sample_index([base[0] + o[0], base[1] + o[1], base[2] + o[2]])];
            *slot = v;
            if is_inside(v) {
                case |= 1 << c;
            }
        }
        (case, corner_value)
    }

    // ─── the bit-sliced path ───────────────────────────────────────────────

    /// One byte per lane, all eight lanes.
    const REPLICATE: u64 = 0x0101_0101_0101_0101;
    /// Lane `j` masks bit `j`: `0x01, 0x02, 0x04, …, 0x80` from the low byte up.
    const LANE_BITS: u64 = 0x8040_2010_0804_0201;
    /// `0x7f` per lane. `0x80 + 0x7f = 0xff`, so the add below cannot carry out
    /// of a lane, which is the whole reason this spread is used rather than the
    /// shorter magic-multiply one.
    const HALF: u64 = 0x7f7f_7f7f_7f7f_7f7f;

    /// **Eight case indices, one per byte, from the eight sign bit-planes.**
    ///
    /// `out` byte `j` is the case index of cell `8·group + j`. Eight iterations,
    /// each ending in the shift-and-OR pair the registration names.
    #[inline]
    fn cases_word(plane: &[u64; 8], group: usize) -> u64 {
        let mut out = 0u64;
        for (c, &p) in plane.iter().enumerate() {
            // Eight cells' sign bits for this corner, in the low byte.
            let byte = (p >> (8 * group)) & 0xff;
            // Replicate into eight lanes and keep lane `j`'s own bit, so lane
            // `j` is `1 << j` when set and zero otherwise.
            let masked = (byte * REPLICATE) & LANE_BITS;
            // Normalise every set lane to exactly 1. Carry-free: a lane holds
            // at most `0x80` and `0x80 + 0x7f = 0xff`.
            let ones = ((masked + HALF) >> 7) & REPLICATE;
            out |= ones << c;
        }
        out
    }

    /// What the vacuity control breaks, and how.
    #[derive(Clone, Copy)]
    enum Corrupt {
        /// The mechanism as built.
        None,
        /// Invert one whole sign plane — the registration's own control. Flips
        /// bit `plane` of every cell's case index, so the expected reading is
        /// `cells`.
        Plane(usize),
        /// Flip one bit of one sign plane, at one cell. The tighter control:
        /// the expected reading is exactly 1, which is what says the comparator
        /// is per-cell rather than a length or a checksum.
        Bit(usize, [usize; 3]),
    }

    /// The eight corner planes of cell word `w` in cell row `(y, z)`.
    ///
    /// Four sample rows; the `dx = 1` planes are `dual.rs:395`'s
    /// `inside_word_shifted`, `(lo >> 1) | (hi << 63)`.
    #[inline]
    fn planes_of(
        inside: &[u64],
        bit_row: usize,
        n: usize,
        w: usize,
        y: usize,
        z: usize,
    ) -> [u64; 8] {
        let mut plane = [0u64; 8];
        let rows = [
            y + n * z,
            (y + 1) + n * z,
            y + n * (z + 1),
            (y + 1) + n * (z + 1),
        ];
        for (j, &r) in rows.iter().enumerate() {
            let base = bit_row * r;
            let lo = inside[base + w];
            let hi = if w + 1 < bit_row {
                inside[base + w + 1]
            } else {
                0
            };
            plane[2 * j] = lo;
            plane[2 * j + 1] = (lo >> 1) | (hi << 63);
        }
        plane
    }

    /// Apply the vacuity control to one word's planes, if it targets this word.
    #[inline]
    fn corrupt_planes(plane: &mut [u64; 8], corrupt: Corrupt, w: usize, y: usize, z: usize) {
        match corrupt {
            Corrupt::None => {}
            Corrupt::Plane(p) => plane[p] = !plane[p],
            Corrupt::Bit(p, [tx, ty, tz]) => {
                if y == ty && z == tz && w == tx / 64 {
                    plane[p] ^= 1 << (tx % 64);
                }
            }
        }
    }

    // ─── the thing under test, and its comparands ──────────────────────────

    /// Every arm's state. One `values` array, read by all of them.
    struct Classify<R: Real> {
        /// Filled once per row, **outside** every counter window.
        values: Vec<R>,
        /// The byte path's dense deliverable: one case byte per cell.
        cases_byte: Vec<u8>,
        /// The bit-sliced path's, and the population `cases_identical` counts.
        cases_bitsliced: Vec<u8>,
        /// The eight-lane arm's.
        cases_swar8: Vec<u8>,
        /// The vacuity control's, written outside every window.
        cases_control: Vec<u8>,
        /// `dual.rs:359-381`'s sign bitmap, one bit per sample along `x`.
        inside: Vec<u64>,
        /// One `0`/`1` sign **byte** per sample, padded so an eight-sample load
        /// at the last cell of the last row stays in bounds.
        signs: Vec<u8>,
        /// The byte path's compaction — `P-121`'s `classify` stage, and C1's
        /// numerator.
        active: Vec<[u32; 3]>,
        /// The active arm's output: one case byte per active cell.
        active_cases: Vec<u8>,
    }

    /// Bytes of padding on [`Classify::signs`].
    ///
    /// The eight-lane arm loads eight sign bytes starting at the last cell of
    /// the last cell row of the last plane, which reaches `n³ + 7`.
    const SIGN_PAD: usize = 16;

    impl<R: Real> Classify<R> {
        fn new() -> Self {
            Self {
                values: Vec::new(),
                cases_byte: Vec::new(),
                cases_bitsliced: Vec::new(),
                cases_swar8: Vec::new(),
                cases_control: Vec::new(),
                inside: Vec::new(),
                signs: Vec::new(),
                active: Vec::new(),
                active_cases: Vec::new(),
            }
        }

        /// `sdf::sample_grid` with `row_stride == size[0]`, which is what
        /// `MarchingCubes::extract` passes. **Outside every window.**
        fn sample<S: Sdf<Scalar = R>>(&mut self, sdf: &S, g: Grid<R>) {
            self.values.clear();
            self.values.reserve(g.samples());
            for z in 0..g.n {
                for y in 0..g.n {
                    for x in 0..g.n {
                        self.values.push(sdf.sample([
                            g.origin[0] + g.cell_size * R::from_f64(f64::from(x)),
                            g.origin[1] + g.cell_size * R::from_f64(f64::from(y)),
                            g.origin[2] + g.cell_size * R::from_f64(f64::from(z)),
                        ]));
                    }
                }
            }
            let cells = g.cell_count();
            for buffer in [
                &mut self.cases_byte,
                &mut self.cases_bitsliced,
                &mut self.cases_swar8,
                &mut self.cases_control,
            ] {
                buffer.clear();
                buffer.resize(cells, 0);
            }
            self.inside.clear();
            self.inside
                .resize(g.bit_row() * g.n as usize * g.n as usize, 0);
            self.signs.clear();
            self.signs.resize(g.samples() + SIGN_PAD, 0);
            self.active.clear();
            self.active.reserve(cells / 8);
            self.active_cases.clear();
            self.active_cases.reserve(cells / 8);
        }

        /// **The comparand, dense.** One case byte per cell, so the deliverable
        /// is the one the two dense arms produce.
        fn byte_dense(&mut self, g: Grid<R>) {
            let c = g.cells();
            let mut i = 0usize;
            for z in 0..c as u32 {
                for y in 0..c as u32 {
                    for x in 0..c as u32 {
                        self.cases_byte[i] = case_of(&self.values, g, [x, y, z]);
                        i += 1;
                    }
                }
            }
            black_box(&self.cases_byte);
        }

        /// **The comparand, shipped shape, and C1's numerator.**
        ///
        /// `mod.rs:258-268` plus the compaction: `case == 0` and `case == 255`
        /// are exactly the cells the table gives no triangles for, so dropping
        /// them here is the same decision `extract`'s
        /// `entry.count == 0 { continue }` makes, taken without reading the
        /// table. This is `P-121`'s `classify` stage boundary, kept so the two
        /// rows' shares are comparable.
        fn byte_compact(&mut self, g: Grid<R>) {
            self.active.clear();
            let c = g.cells() as u32;
            for z in 0..c {
                for y in 0..c {
                    for x in 0..c {
                        let base = [x, y, z];
                        let case = case_of(&self.values, g, base);
                        if case != 0 && case != u8::MAX {
                            self.active.push(base);
                        }
                    }
                }
            }
            black_box(&self.active);
        }

        /// `dual.rs:359-381`, on the Marching Cubes sample layout: one sign bit
        /// per sample, 64 to a `u64`, along `x` only.
        ///
        /// **One `is_inside` per sample.** The byte path evaluates it
        /// `8·(n−1)³` times; this evaluates it `n³` times, and that is the
        /// mechanism rather than the assembly trick downstream of it.
        #[inline]
        fn build_inside(values: &[R], inside: &mut [u64], g: Grid<R>) {
            let n = g.n as usize;
            let bit_row = g.bit_row();
            for r in 0..n * n {
                let src = n * r;
                let dst = bit_row * r;
                for w in 0..bit_row {
                    let base = w * 64;
                    let take = (n - base).min(64);
                    let mut word = 0u64;
                    for k in 0..take {
                        word |= u64::from(is_inside(values[src + base + k])) << k;
                    }
                    inside[dst + w] = word;
                }
            }
        }

        /// **The registered arm.** Eight sign bit-planes, then eight
        /// shift-and-OR pairs per `u64` of case indices, eight `u64`s per
        /// 64-cell word.
        fn bitsliced_dense(&mut self, g: Grid<R>) {
            let Self {
                values,
                inside,
                cases_bitsliced,
                ..
            } = self;
            Self::build_inside(values, inside, g);
            bitsliced_cases(inside, g, cases_bitsliced, Corrupt::None);
            black_box(&self.cases_bitsliced);
        }

        /// **What a shared classifier would return**, rather than 64 dense
        /// bytes: the fused `any & !all` active word (`dual.rs:424`) and a case
        /// byte for the active cells only, walked out with `trailing_zeros` and
        /// `x &= x − 1` in ascending `x` (`dual.rs:490-497`).
        fn bitsliced_active(&mut self, g: Grid<R>) {
            let Self {
                values,
                inside,
                active_cases,
                ..
            } = self;
            Self::build_inside(values, inside, g);
            active_cases.clear();
            let n = g.n as usize;
            let bit_row = g.bit_row();
            let cells = g.cells();
            let cell_words = g.cell_words();
            for z in 0..cells {
                for y in 0..cells {
                    for w in 0..cell_words {
                        let plane = planes_of(inside, bit_row, n, w, y, z);
                        let mut any = 0u64;
                        let mut all = !0u64;
                        for &p in &plane {
                            any |= p;
                            all &= p;
                        }
                        let mut bits = (any & !all) & cell_mask(w, cells);
                        while bits != 0 {
                            let k = bits.trailing_zeros() as usize;
                            bits &= bits - 1;
                            let mut case = 0u8;
                            for (c, &p) in plane.iter().enumerate() {
                                case |= (((p >> k) & 1) as u8) << c;
                            }
                            active_cases.push(case);
                        }
                    }
                }
            }
            black_box(&self.active_cases);
        }

        /// **The other reading of the registration's sentence**, eight lanes
        /// wide: one sign *byte* per sample, so the `dx = 1` plane is a
        /// byte-offset load and the assembly is literally
        /// `p0 | p1 << 1 | … | p7 << 7`.
        fn swar8_dense(&mut self, g: Grid<R>) {
            let Self {
                values,
                signs,
                cases_swar8,
                ..
            } = self;
            for (slot, &v) in signs.iter_mut().zip(values.iter()) {
                *slot = u8::from(is_inside(v));
            }
            let cells = g.cells();
            for z in 0..cells {
                for y in 0..cells {
                    let rows = [
                        g.sample_row(y, z),
                        g.sample_row(y + 1, z),
                        g.sample_row(y, z + 1),
                        g.sample_row(y + 1, z + 1),
                    ];
                    let row = &mut cases_swar8[cells * (y + cells * z)..][..cells];
                    let mut off = 0usize;
                    while off < cells {
                        let load = |base: usize| {
                            u64::from_le_bytes(
                                signs[base + off..base + off + 8]
                                    .try_into()
                                    .expect("eight sign bytes"),
                            )
                        };
                        let mut word = 0u64;
                        for (j, &r) in rows.iter().enumerate() {
                            word |= load(r) << (2 * j);
                            word |= load(r + 1) << (2 * j + 1);
                        }
                        let take = (cells - off).min(8);
                        row[off..off + take].copy_from_slice(&word.to_le_bytes()[..take]);
                        off += 8;
                    }
                }
            }
            black_box(&self.cases_swar8);
        }

        /// The vacuity control: the same computation with one sign plane
        /// broken. **Outside every counter window** — it measures nothing, it
        /// establishes that the comparator can fail.
        fn control(&mut self, g: Grid<R>, corrupt: Corrupt) -> usize {
            let Self {
                values,
                inside,
                cases_control,
                cases_byte,
                ..
            } = self;
            Self::build_inside(values, inside, g);
            bitsliced_cases(inside, g, cases_control, corrupt);
            cases_byte
                .iter()
                .zip(cases_control.iter())
                .filter(|(a, b)| a != b)
                .count()
        }

        /// Every bit of the packed plane against `is_inside(values[s])`
        /// directly. C3's third layer, and the one that covers the algorithms
        /// which never form a case index.
        fn sign_bits_identical(&mut self, g: Grid<R>) -> usize {
            let Self { values, inside, .. } = self;
            Self::build_inside(values, inside, g);
            let n = g.n as usize;
            let bit_row = g.bit_row();
            let mut same = 0usize;
            for z in 0..n {
                for y in 0..n {
                    let src = g.sample_row(y, z);
                    let dst = bit_row * (y + n * z);
                    for x in 0..n {
                        let bit = (inside[dst + x / 64] >> (x % 64)) & 1;
                        if (bit == 1) == is_inside(values[src + x]) {
                            same += 1;
                        }
                    }
                }
            }
            same
        }
    }

    /// `dual.rs:445`'s `cell_mask`: the tail word of a cell row is short.
    #[inline]
    const fn cell_mask(w: usize, cells: usize) -> u64 {
        let remaining = cells.saturating_sub(w * 64);
        if remaining >= 64 {
            !0
        } else {
            (1u64 << remaining) - 1
        }
    }

    /// The dense case bytes of every cell, from the sign bit-planes.
    ///
    /// Free rather than a method so the vacuity control can write into a second
    /// buffer while reading the same planes.
    fn bitsliced_cases<R: Real>(inside: &[u64], g: Grid<R>, out: &mut [u8], corrupt: Corrupt) {
        let n = g.n as usize;
        let bit_row = g.bit_row();
        let cells = g.cells();
        let cell_words = g.cell_words();
        for z in 0..cells {
            for y in 0..cells {
                let row = &mut out[cells * (y + cells * z)..][..cells];
                for w in 0..cell_words {
                    let mut plane = planes_of(inside, bit_row, n, w, y, z);
                    corrupt_planes(&mut plane, corrupt, w, y, z);
                    for group in 0..8 {
                        let off = w * 64 + group * 8;
                        if off >= cells {
                            break;
                        }
                        let word = cases_word(&plane, group);
                        let take = (cells - off).min(8);
                        row[off..off + take].copy_from_slice(&word.to_le_bytes()[..take]);
                    }
                }
            }
        }
    }

    // ─── the golden arm: Marching Cubes with a swappable classifier ─────────

    /// Which classifier fills the control stream.
    #[derive(Clone, Copy, PartialEq, Eq)]
    enum Source {
        Byte,
        BitSliced,
    }

    /// A bench-local Marching Cubes, at `f64` because `validate::mesh_hash`
    /// hashes `MeshBuffer<f64>`.
    ///
    /// Two passes: fill a dense case-byte control stream from one classifier,
    /// then emit the payload off it. [`Self::emit_cell`] is the *same*
    /// instructions either way, so the only difference between the two meshes
    /// is where the case byte came from — which is what makes C3 a bit-for-bit
    /// equality rather than a tolerance.
    struct March {
        values: Vec<f64>,
        cases: Vec<u8>,
        inside: Vec<u64>,
        edge_vertices: Vec<u32>,
        mesh: MeshBuffer<f64>,
    }

    impl March {
        fn new() -> Self {
            Self {
                values: Vec::new(),
                cases: Vec::new(),
                inside: Vec::new(),
                edge_vertices: Vec::new(),
                mesh: MeshBuffer::new(),
            }
        }

        fn sample<S: Sdf<Scalar = f64>>(&mut self, sdf: &S, g: Grid<f64>) {
            self.values.clear();
            self.values.reserve(g.samples());
            for z in 0..g.n {
                for y in 0..g.n {
                    for x in 0..g.n {
                        self.values.push(sdf.sample([
                            g.origin[0] + g.cell_size * f64::from(x),
                            g.origin[1] + g.cell_size * f64::from(y),
                            g.origin[2] + g.cell_size * f64::from(z),
                        ]));
                    }
                }
            }
            self.cases.clear();
            self.cases.resize(g.cell_count(), 0);
            self.inside.clear();
            self.inside
                .resize(g.bit_row() * g.n as usize * g.n as usize, 0);
        }

        /// `MeshSink::vertex`.
        #[inline]
        fn vertex<S: Sdf<Scalar = f64>>(&mut self, sdf: &S, position: [f64; 3]) -> u32 {
            let index = self.mesh.positions.len() as u32;
            self.mesh.positions.push(position);
            self.mesh.normals.push(unit_gradient(sdf, position));
            index
        }

        /// `MarchingCubes::vertex_on_edge` — the edge-cache probe.
        #[inline]
        fn vertex_on_edge<S: Sdf<Scalar = f64>>(
            &mut self,
            sdf: &S,
            g: Grid<f64>,
            base: [u32; 3],
            edge: u8,
            corner_value: &[f64; 8],
        ) -> u32 {
            let axis = EDGE_AXIS[edge as usize] as usize;
            let o = corner_offset(EDGE_CORNERS[edge as usize][0]);
            let lo_sample = g.sample_index([base[0] + o[0], base[1] + o[1], base[2] + o[2]]);
            let key = lo_sample * 3 + axis;
            let cached = self.edge_vertices[key];
            if cached != u32::MAX {
                return cached;
            }
            let position = edge_position(base, edge, corner_value, g.origin, g.cell_size);
            let index = self.vertex(sdf, position);
            self.edge_vertices[key] = index;
            index
        }

        /// `mod.rs:307-377` under `MarchingCubes::new`'s defaults —
        /// `FaceAmbiguity::Separate`, `InteriorAmbiguity::Ignore`,
        /// `crossing_refinement == 0` — so `ambiguous` is zero, `mask` is zero
        /// and the entry is the derived table's.
        #[inline]
        fn emit_cell<S: Sdf<Scalar = f64>>(
            &mut self,
            sdf: &S,
            g: Grid<f64>,
            base: [u32; 3],
            case: u8,
            corner_value: &[f64; 8],
        ) {
            let entry = &CASES[case as usize];
            if entry.count == 0 {
                return;
            }
            // Cycle centroids first (A-015): a triangle naming one needs every
            // edge vertex of that cycle averaged before it can be emitted.
            let mut centroid = [0u32; MAX_CENTROIDS];
            for (c, slot) in centroid
                .iter_mut()
                .enumerate()
                .take(entry.centroids as usize)
            {
                let code = CENTROID_BASE + c as u8;
                let mut sum = [0.0f64; 3];
                let mut n = 0u32;
                for tri in &entry.triangles[..entry.count as usize] {
                    if tri[0] != code {
                        continue;
                    }
                    let p = edge_position(base, tri[1], corner_value, g.origin, g.cell_size);
                    sum = [sum[0] + p[0], sum[1] + p[1], sum[2] + p[2]];
                    n += 1;
                }
                let scale = f64::from(n).recip();
                let position = [sum[0] * scale, sum[1] * scale, sum[2] * scale];
                *slot = self.vertex(sdf, position);
            }
            for tri in &entry.triangles[..entry.count as usize] {
                let mut idx = [0u32; 3];
                for (k, &code) in tri.iter().enumerate() {
                    idx[k] = if is_centroid(code) {
                        centroid[(code - CENTROID_BASE) as usize]
                    } else {
                        self.vertex_on_edge(sdf, g, base, code, corner_value)
                    };
                }
                self.mesh.indices.extend_from_slice(&idx);
            }
        }

        /// One extraction, classified by `source`.
        fn march<S: Sdf<Scalar = f64>>(&mut self, sdf: &S, g: Grid<f64>, source: Source) {
            self.edge_vertices.clear();
            self.edge_vertices.resize(g.samples() * 3, u32::MAX);
            self.mesh.reset();

            match source {
                Source::Byte => {
                    let c = g.cells() as u32;
                    let mut i = 0usize;
                    for z in 0..c {
                        for y in 0..c {
                            for x in 0..c {
                                self.cases[i] = case_of(&self.values, g, [x, y, z]);
                                i += 1;
                            }
                        }
                    }
                }
                Source::BitSliced => {
                    let Self {
                        values,
                        inside,
                        cases,
                        ..
                    } = self;
                    Classify::<f64>::build_inside(values, inside, g);
                    bitsliced_cases(inside, g, cases, Corrupt::None);
                }
            }

            let c = g.cells() as u32;
            let mut i = 0usize;
            for z in 0..c {
                for y in 0..c {
                    for x in 0..c {
                        let case = self.cases[i];
                        i += 1;
                        if CASES[case as usize].count == 0 {
                            continue;
                        }
                        let base = [x, y, z];
                        let (byte_case, corner_value) = case_and_corners(&self.values, g, base);
                        debug_assert_eq!(byte_case, case);
                        self.emit_cell(sdf, g, base, case, &corner_value);
                    }
                }
            }
        }
    }

    // ─── the committed golden fixture ──────────────────────────────────────

    /// The value of `"key"` in one `{…}` chunk of `golden_hashes.json`, as a
    /// string.
    fn json_field(chunk: &str, key: &str) -> String {
        let needle = format!("\"{key}\":");
        let at = chunk
            .find(&needle)
            .unwrap_or_else(|| panic!("golden_hashes.json entry has no {key}"))
            + needle.len();
        let rest = &chunk[at..];
        if let Some(stripped) = rest.strip_prefix('"') {
            let end = stripped.find('"').expect("a closed string");
            stripped[..end].to_string()
        } else {
            let end = rest
                .find(|c: char| !c.is_ascii_digit())
                .unwrap_or(rest.len());
            rest[..end].to_string()
        }
    }

    /// The 24 committed `marching_cubes` rows of `M-31`'s 216, as
    /// `(field, samples, hash)`.
    ///
    /// Read from `crates/isomesh/golden_hashes.json`, the file
    /// `golden_hashes_are_unchanged` (`golden/tests.rs:59`) gates against, so
    /// `hashes_moved` is movement in **the** fixture rather than in a
    /// re-derivation of it.
    fn golden_marching_cubes() -> Vec<(String, u32, u64)> {
        let path = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("golden_hashes.json");
        let text = std::fs::read_to_string(&path)
            .unwrap_or_else(|e| panic!("reading {}: {e}", path.display()));
        let mut out = Vec::new();
        for chunk in text.split('{').skip(1) {
            if json_field(chunk, "algorithm") != "marching_cubes" {
                continue;
            }
            let field = json_field(chunk, "field");
            let samples: u32 = json_field(chunk, "samples").parse().expect("a resolution");
            let hash = u64::from_str_radix(&json_field(chunk, "hash"), 16).expect("a hex hash");
            out.push((field, samples, hash));
        }
        out
    }

    // ─── counting ──────────────────────────────────────────────────────────

    /// What one or more counter windows read.
    #[derive(Clone, Copy, Default)]
    struct Counted {
        cycles: f64,
        instructions: f64,
    }

    impl Counted {
        fn scaled(self, by: f64) -> Self {
            Self {
                cycles: self.cycles * by,
                instructions: self.instructions * by,
            }
        }
    }

    /// One counter window, undivided. The `perf_event` system calls are all
    /// outside the counted region.
    fn raw_window(probe: &mut Probe, body: impl FnOnce()) -> (Counted, f64) {
        probe.reset_and_enable();
        let started = Instant::now();
        body();
        let nanos = started.elapsed().as_nanos() as f64;
        probe.disable();
        let counted = probe.read();
        assert!(
            counted.worst_ratio() >= MIN_TIME_RATIO,
            "a counter ran only {:.1}% of the time it was enabled, so its value is an \
             extrapolation rather than a measurement",
            counted.worst_ratio() * 100.0
        );
        (
            Counted {
                cycles: counted.cycles.count as f64,
                instructions: counted.instructions.count as f64,
            },
            nanos,
        )
    }

    /// [`raw_window`] over `inner` repetitions, divided by `inner`.
    fn window(probe: &mut Probe, inner: usize, mut body: impl FnMut()) -> (Counted, f64) {
        let scale = 1.0 / inner as f64;
        let (counted, nanos) = raw_window(probe, || {
            for _ in 0..inner {
                body();
            }
        });
        (counted.scaled(scale), nanos * scale)
    }

    fn median_of(values: &mut [f64]) -> f64 {
        assert!(!values.is_empty(), "a median of nothing");
        values.sort_by(f64::total_cmp);
        values[values.len() / 2]
    }

    fn median(pick: &dyn Fn(&Rep) -> f64, reps: &[Rep]) -> f64 {
        let mut values: Vec<f64> = reps.iter().map(pick).collect();
        median_of(&mut values)
    }

    fn median_counted(pick: &dyn Fn(&Rep) -> Counted, reps: &[Rep]) -> Counted {
        Counted {
            cycles: median(&|r| pick(r).cycles, reps),
            instructions: median(&|r| pick(r).instructions, reps),
        }
    }

    /// The six sibling windows of one repetition.
    #[derive(Clone, Copy)]
    enum Arm {
        ByteDense,
        ByteCompact,
        BitslicedDense,
        BitslicedActive,
        Swar8Dense,
        Extract,
    }

    const ARMS: [Arm; 6] = [
        Arm::ByteDense,
        Arm::ByteCompact,
        Arm::BitslicedDense,
        Arm::BitslicedActive,
        Arm::Swar8Dense,
        Arm::Extract,
    ];

    /// One repetition: six sibling windows, in a rotated order.
    #[derive(Clone, Copy, Default)]
    struct Rep {
        byte_dense: Counted,
        byte_compact: Counted,
        bitsliced_dense: Counted,
        bitsliced_active: Counted,
        swar8_dense: Counted,
        extract: Counted,
        byte_dense_ns: f64,
        extract_ns: f64,
    }

    impl Rep {
        fn instruction_ratio(&self) -> f64 {
            self.bitsliced_dense.instructions / self.byte_dense.instructions
        }

        fn cycle_ratio(&self) -> f64 {
            self.bitsliced_dense.cycles / self.byte_dense.cycles
        }
    }

    // ─── one row ───────────────────────────────────────────────────────────

    /// One measured `(field, resolution, scalar)`.
    struct Row {
        field: &'static str,
        resolution: u32,
        scalar: &'static str,
        cells: usize,
        samples: usize,
        active_cells: usize,
        inner: usize,
        inner_extract: usize,
        byte_dense: Counted,
        byte_compact: Counted,
        bitsliced_dense: Counted,
        bitsliced_active: Counted,
        swar8_dense: Counted,
        extract: Counted,
        byte_dense_ns: f64,
        extract_ns: f64,
        instruction_ratio_spread: f64,
        cycle_ratio_best: f64,
        cycle_ratio_worst: f64,
        cases_identical: usize,
        swar8_identical: usize,
        active_cases_identical: usize,
        corrupt_control_mismatches: usize,
        single_bit_control_mismatches: usize,
        sign_bits_identical: usize,
        inside_bytes: usize,
    }

    impl Row {
        fn per_cell(&self, pick: impl Fn(Counted) -> f64, arm: Counted) -> f64 {
            let _ = &pick;
            pick(arm) / self.cells as f64
        }

        fn instructions_per_cell_byte(&self) -> f64 {
            self.byte_dense.instructions / self.cells as f64
        }

        fn instructions_per_cell_bitsliced(&self) -> f64 {
            self.bitsliced_dense.instructions / self.cells as f64
        }

        /// C2's bar is 0.5 on this.
        fn instruction_ratio(&self) -> f64 {
            self.bitsliced_dense.instructions / self.byte_dense.instructions
        }

        fn instruction_ratio_swar8(&self) -> f64 {
            self.swar8_dense.instructions / self.byte_dense.instructions
        }

        fn instruction_ratio_active(&self) -> f64 {
            self.bitsliced_active.instructions / self.byte_dense.instructions
        }

        fn cycle_ratio(&self) -> f64 {
            self.bitsliced_dense.cycles / self.byte_dense.cycles
        }

        /// C1's bar is 0.10 on this. The registered form, in cycles, because
        /// `P-121`'s decomposition is in cycles.
        fn classify_share(&self) -> f64 {
            self.byte_compact.cycles / self.extract.cycles
        }

        fn classify_share_instructions(&self) -> f64 {
            self.byte_compact.instructions / self.extract.instructions
        }

        fn c1_holds(&self) -> bool {
            self.classify_share() >= CLASSIFY_SHARE_BAR
        }

        fn c2_holds(&self) -> bool {
            self.instruction_ratio() < INSTRUCTION_RATIO_BAR
        }
    }

    /// Measure one `(field, resolution, scalar)`.
    fn measure<R, S>(
        field: &'static str,
        scalar: &'static str,
        n: u32,
        sdf: &S,
        origin: [R; 3],
        cell_size: R,
    ) -> Row
    where
        R: Real,
        S: Sdf<Scalar = R>,
    {
        let g = Grid {
            n,
            origin,
            cell_size,
        };
        let shape = RuntimeShape3::new([n; 3]).expect("the fixture fits u32");

        let mut c = Classify::<R>::new();
        c.sample(sdf, g);

        let mut shipped = MarchingCubes::<R>::new();
        let mut out = MeshBuffer::<R>::new();

        for _ in 0..WARMUP {
            c.byte_dense(g);
            c.byte_compact(g);
            c.bitsliced_dense(g);
            c.bitsliced_active(g);
            c.swar8_dense(g);
            out.reset();
            shipped
                .extract(sdf, &shape, origin, cell_size, &mut out)
                .expect("extraction");
        }

        // ── the two batches, each chosen from its own timed pass ─────────────
        let started = Instant::now();
        c.byte_dense(g);
        let classify_ns = started.elapsed().as_nanos() as f64;
        let inner = ((TARGET_BATCH_NS / classify_ns.max(1.0)).ceil() as usize).clamp(1, MAX_INNER);

        let started = Instant::now();
        out.reset();
        shipped
            .extract(sdf, &shape, origin, cell_size, &mut out)
            .expect("extraction");
        let extract_ns = started.elapsed().as_nanos() as f64;
        let inner_extract =
            ((TARGET_BATCH_NS / extract_ns.max(1.0)).ceil() as usize).clamp(1, MAX_INNER);

        // ── the equalities and the vacuity controls, outside every window ────
        c.byte_dense(g);
        c.bitsliced_dense(g);
        c.bitsliced_active(g);
        c.swar8_dense(g);
        c.byte_compact(g);

        let cells = g.cell_count();
        let cases_identical = c
            .cases_byte
            .iter()
            .zip(c.cases_bitsliced.iter())
            .filter(|(a, b)| a == b)
            .count();
        let swar8_identical = c
            .cases_byte
            .iter()
            .zip(c.cases_swar8.iter())
            .filter(|(a, b)| a == b)
            .count();
        // The active arm's answers, against the byte path's on the same cells:
        // `case != 0 && case != 255` is exactly `any & !all`, so the two lists
        // are the same cells in the same order.
        let active_cases_identical = c
            .active
            .iter()
            .map(|&base| {
                c.cases_byte[base[0] as usize
                    + g.cells() * (base[1] as usize + g.cells() * base[2] as usize)]
            })
            .zip(c.active_cases.iter().copied())
            .filter(|(a, b)| a == b)
            .count();
        assert_eq!(
            c.active.len(),
            c.active_cases.len(),
            "{field} {n}^3 {scalar}: the fused `any & !all` word found {} active cells and the \
             byte path's compaction found {} — the two are supposed to be the same predicate",
            c.active_cases.len(),
            c.active.len()
        );

        let corrupt_control_mismatches = c.control(g, Corrupt::Plane(3));
        let half = g.cells() / 2;
        let single_bit_control_mismatches = c.control(g, Corrupt::Bit(3, [half, half, half]));
        let sign_bits_identical = c.sign_bits_identical(g);

        // The registration's own vacuity control, asserted rather than merely
        // recorded.
        assert_eq!(
            cases_identical,
            cells,
            "{field} {n}^3 {scalar}: the bit-sliced case index differs from the byte path's on \
             {} of {cells} cells — C3 is an integer identity",
            cells - cases_identical
        );
        assert_eq!(
            corrupt_control_mismatches, cells,
            "{field} {n}^3 {scalar}: inverting sign plane 3 flips bit 3 of every cell's case \
             index, so every one of {cells} cells must mismatch; {corrupt_control_mismatches} \
             did. A comparator that cannot report a mismatch has not reported a match"
        );
        assert_eq!(
            single_bit_control_mismatches, 1,
            "{field} {n}^3 {scalar}: flipping one bit of sign plane 3 must move exactly one \
             cell's case index; {single_bit_control_mismatches} moved. The comparator is not \
             per-cell"
        );
        assert_eq!(
            swar8_identical,
            cells,
            "{field} {n}^3 {scalar}: the eight-lane arm's case index differs from the byte \
             path's on {} of {cells} cells",
            cells - swar8_identical
        );
        assert_eq!(
            active_cases_identical,
            c.active.len(),
            "{field} {n}^3 {scalar}: the active arm's case bytes differ from the byte path's on \
             {} of {} active cells",
            c.active.len() - active_cases_identical,
            c.active.len()
        );
        assert_eq!(
            sign_bits_identical,
            g.samples(),
            "{field} {n}^3 {scalar}: the packed sign plane differs from `is_inside` on {} of {} \
             samples",
            g.samples() - sign_bits_identical,
            g.samples()
        );

        // ── REPS repetitions, six sibling windows each ───────────────────────
        let mut probe = Probe::open();
        let mut reps: Vec<Rep> = Vec::with_capacity(REPS);
        for rep in 0..REPS {
            let mut r = Rep::default();
            for k in 0..ARMS.len() {
                match ARMS[(k + rep) % ARMS.len()] {
                    Arm::ByteDense => {
                        let (counted, ns) = window(&mut probe, inner, || c.byte_dense(g));
                        r.byte_dense = counted;
                        r.byte_dense_ns = ns;
                    }
                    Arm::ByteCompact => {
                        r.byte_compact = window(&mut probe, inner, || c.byte_compact(g)).0;
                    }
                    Arm::BitslicedDense => {
                        r.bitsliced_dense = window(&mut probe, inner, || c.bitsliced_dense(g)).0;
                    }
                    Arm::BitslicedActive => {
                        r.bitsliced_active = window(&mut probe, inner, || c.bitsliced_active(g)).0;
                    }
                    Arm::Swar8Dense => {
                        r.swar8_dense = window(&mut probe, inner, || c.swar8_dense(g)).0;
                    }
                    Arm::Extract => {
                        let (counted, ns) = window(&mut probe, inner_extract, || {
                            out.reset();
                            shipped
                                .extract(sdf, &shape, origin, cell_size, &mut out)
                                .expect("extraction");
                            black_box(&out);
                        });
                        r.extract = counted;
                        r.extract_ns = ns;
                    }
                }
            }
            reps.push(r);
        }

        let ratios: Vec<f64> = reps.iter().map(Rep::instruction_ratio).collect();
        let instruction_ratio_spread = ratios.iter().copied().fold(0.0f64, f64::max)
            - ratios.iter().copied().fold(f64::MAX, f64::min);
        let cycle_ratios: Vec<f64> = reps.iter().map(Rep::cycle_ratio).collect();

        Row {
            field,
            resolution: n,
            scalar,
            cells,
            samples: g.samples(),
            active_cells: c.active.len(),
            inner,
            inner_extract,
            byte_dense: median_counted(&|r| r.byte_dense, &reps),
            byte_compact: median_counted(&|r| r.byte_compact, &reps),
            bitsliced_dense: median_counted(&|r| r.bitsliced_dense, &reps),
            bitsliced_active: median_counted(&|r| r.bitsliced_active, &reps),
            swar8_dense: median_counted(&|r| r.swar8_dense, &reps),
            extract: median_counted(&|r| r.extract, &reps),
            byte_dense_ns: median(&|r| r.byte_dense_ns, &reps),
            extract_ns: median(&|r| r.extract_ns, &reps),
            instruction_ratio_spread,
            cycle_ratio_best: cycle_ratios.iter().copied().fold(f64::MAX, f64::min),
            cycle_ratio_worst: cycle_ratios.iter().copied().fold(0.0f64, f64::max),
            cases_identical,
            swar8_identical,
            active_cases_identical,
            corrupt_control_mismatches,
            single_bit_control_mismatches,
            sign_bits_identical,
            inside_bytes: c.inside.len() * 8,
        }
    }

    /// The registered fixture: eight reference fields × {33³, 65³} ×
    /// {`f32`, `f64`}.
    fn sweep() -> Vec<Row> {
        let mut rows = Vec::new();
        for n in RESOLUTIONS {
            isomesh::for_each_reference_field!(f32, |name, field| {
                let (_, origin, cell_size) = crate::common::grid(&field, n);
                rows.push(measure(name, "f32", n, &field, origin, cell_size));
            });
            isomesh::for_each_reference_field!(f64, |name, field| {
                let (_, origin, cell_size) = crate::common::grid(&field, n);
                rows.push(measure(name, "f64", n, &field, origin, cell_size));
            });
        }
        rows
    }

    // ─── C3 ────────────────────────────────────────────────────────────────

    /// What the golden arm established.
    struct Golden {
        /// `(marching_cubes, field, samples)` rows of the committed fixture.
        hash_rows: usize,
        /// How many of those the byte-classified mirror reproduced. Asserted
        /// equal to `hash_rows`: `M-279`'s agreement check, and the licence for
        /// reading anything off the mirror.
        reproduced: usize,
        /// **The registered column.** How many committed hashes the
        /// bit-sliced-classified mirror fails to reproduce.
        moved: usize,
        /// Cells whose case index was compared, over every golden grid.
        cells_checked: usize,
        case_mismatches: usize,
        samples_checked: usize,
        sign_mismatches: usize,
    }

    fn golden() -> Golden {
        let committed = golden_marching_cubes();
        assert_eq!(
            committed.len(),
            REFERENCE_FIELDS * GOLDEN_RESOLUTIONS.len(),
            "golden_hashes.json should carry {} `marching_cubes` rows — eight reference fields \
             times three resolutions — and carries {}",
            REFERENCE_FIELDS * GOLDEN_RESOLUTIONS.len(),
            committed.len()
        );

        let mut g = Golden {
            hash_rows: committed.len(),
            reproduced: 0,
            moved: 0,
            cells_checked: 0,
            case_mismatches: 0,
            samples_checked: 0,
            sign_mismatches: 0,
        };

        let mut march = March::new();
        isomesh::for_each_reference_field!(f64, |name, field| {
            for n in GOLDEN_RESOLUTIONS {
                let (_, origin, cell_size) = crate::common::grid(&field, n);
                let grid = Grid {
                    n,
                    origin,
                    cell_size,
                };
                march.sample(&field, grid);

                march.march(&field, grid, Source::Byte);
                let hash_byte = mesh_hash(&march.mesh);
                march.march(&field, grid, Source::BitSliced);
                let hash_bitsliced = mesh_hash(&march.mesh);

                let (_, _, want) = committed
                    .iter()
                    .find(|(f, s, _)| f == name && *s == n)
                    .unwrap_or_else(|| panic!("no committed marching_cubes hash for {name} {n}"))
                    .clone();
                assert_eq!(
                    hash_byte, want,
                    "{name} {n}^3: the byte-classified mirror hashes {hash_byte:016x} against \
                     the committed golden {want:016x} — then it is not Marching Cubes and the \
                     bit-sliced comparison below is between two things neither of which ships"
                );
                g.reproduced += 1;
                if hash_bitsliced != want {
                    g.moved += 1;
                }
            }
        });

        // Layers 2 and 3: the case index and the sign predicate, exhaustively,
        // on every golden grid in both scalars.
        for n in GOLDEN_RESOLUTIONS {
            isomesh::for_each_reference_field!(f32, |_name, field| {
                let (_, origin, cell_size) = crate::common::grid(&field, n);
                let grid = Grid {
                    n,
                    origin,
                    cell_size,
                };
                let mut c = Classify::<f32>::new();
                c.sample(&field, grid);
                c.byte_dense(grid);
                c.bitsliced_dense(grid);
                g.cells_checked += grid.cell_count();
                g.case_mismatches += c
                    .cases_byte
                    .iter()
                    .zip(c.cases_bitsliced.iter())
                    .filter(|(a, b)| a != b)
                    .count();
                g.samples_checked += grid.samples();
                g.sign_mismatches += grid.samples() - c.sign_bits_identical(grid);
            });
            isomesh::for_each_reference_field!(f64, |_name, field| {
                let (_, origin, cell_size) = crate::common::grid(&field, n);
                let grid = Grid {
                    n,
                    origin,
                    cell_size,
                };
                let mut c = Classify::<f64>::new();
                c.sample(&field, grid);
                c.byte_dense(grid);
                c.bitsliced_dense(grid);
                g.cells_checked += grid.cell_count();
                g.case_mismatches += c
                    .cases_byte
                    .iter()
                    .zip(c.cases_bitsliced.iter())
                    .filter(|(a, b)| a != b)
                    .count();
                g.samples_checked += grid.samples();
                g.sign_mismatches += grid.samples() - c.sign_bits_identical(grid);
            });
        }

        g
    }

    // ─── the run ───────────────────────────────────────────────────────────

    pub(crate) fn run(run: &mut crate::common::experiment::Run) {
        // C3 first: it is an equality over an enumerated population and it
        // costs seconds, so a broken classifier fails before the sweep spends
        // a minute measuring it.
        let gold = golden();
        assert_eq!(
            gold.case_mismatches, 0,
            "the bit-sliced case index differs from the byte path's on {} of {} cells of the \
             golden grids",
            gold.case_mismatches, gold.cells_checked
        );
        assert_eq!(
            gold.sign_mismatches, 0,
            "the packed sign plane differs from `is_inside` on {} of {} samples of the golden \
             grids",
            gold.sign_mismatches, gold.samples_checked
        );
        let c3_holds = gold.moved == 0
            && gold.reproduced == gold.hash_rows
            && gold.case_mismatches == 0
            && gold.sign_mismatches == 0;

        println!(
            "C3 golden: {} of {} committed `marching_cubes` hashes reproduced by the \
             byte-classified mirror; hashes_moved by the bit-sliced classifier = {}",
            gold.reproduced, gold.hash_rows, gold.moved
        );
        println!(
            "    case indices compared over the golden grids: {} cells, {} mismatches; sign \
             bits: {} samples, {} mismatches -> {}",
            gold.cells_checked,
            gold.case_mismatches,
            gold.samples_checked,
            gold.sign_mismatches,
            if c3_holds { "HELD" } else { "FALSIFIED" }
        );
        println!(
            "    the other {} of M-31's 216 rows are covered by those two equalities, not by a \
             mesh hash: eight of the nine algorithms are not mirrored here, and an algorithm \
             reading the same eight-sign integer on every cell of the same grid cannot produce \
             a different mesh\n",
            216 - gold.hash_rows
        );

        let rows = sweep();

        // ── C1: the share gate, checked before C2 is believed ────────────────
        let c1_rows = rows.iter().filter(|r| r.c1_holds()).count();
        let c1_rows_instructions = rows
            .iter()
            .filter(|r| r.classify_share_instructions() >= CLASSIFY_SHARE_BAR)
            .count();
        let c1_all = rows.iter().all(Row::c1_holds);

        // ── C2: the deterministic clause, scored on the worst row ────────────
        let worst_instruction_ratio = rows
            .iter()
            .map(Row::instruction_ratio)
            .fold(0.0f64, f64::max);
        let best_instruction_ratio = rows
            .iter()
            .map(Row::instruction_ratio)
            .fold(f64::MAX, f64::min);
        let worst_ratio_swar8 = rows
            .iter()
            .map(Row::instruction_ratio_swar8)
            .fold(0.0f64, f64::max);
        let c2_all = rows.iter().all(Row::c2_holds);

        let popcnt = cfg!(target_feature = "popcnt");

        println!(
            "{:<16} {:>4} {:>4} {:>7} {:>7} {:>7} {:>9} {:>9} {:>7} {:>7} {:>7} {:>7} {:>5}",
            "field",
            "n",
            "R",
            "clsf%",
            "clsf%i",
            "C1",
            "i/c byte",
            "i/c bsl",
            "i ratio",
            "swar8",
            "active",
            "c ratio",
            "C2"
        );
        for r in &rows {
            println!(
                "{:<16} {:>4} {:>4} {:>7.4} {:>7.4} {:>7} {:>9.3} {:>9.3} {:>7.4} {:>7.4} \
                 {:>7.4} {:>7.4} {:>5}",
                r.field,
                r.resolution,
                r.scalar,
                r.classify_share(),
                r.classify_share_instructions(),
                r.c1_holds(),
                r.instructions_per_cell_byte(),
                r.instructions_per_cell_bitsliced(),
                r.instruction_ratio(),
                r.instruction_ratio_swar8(),
                r.instruction_ratio_active(),
                r.cycle_ratio(),
                r.c2_holds()
            );
        }

        println!(
            "\nC1 share gate: {c1_rows} of {} rows clear {CLASSIFY_SHARE_BAR} in cycles \
             ({c1_rows_instructions} in instructions) -> {}",
            rows.len(),
            if c1_all { "HELD" } else { "FALSIFIED" }
        );
        println!(
            "    scored per row, x51's rule. The denominator is the shipped \
             MarchingCubes::extract, which P-121 measured at 15-24% MORE than its \
             classify-then-compact mirror, so this share is that much lower than P-121's."
        );
        println!(
            "C2 instructions: worst row {worst_instruction_ratio:.4}, best \
             {best_instruction_ratio:.4} (bar < {INSTRUCTION_RATIO_BAR}) -> {}",
            if c2_all { "HELD" } else { "FALSIFIED" }
        );
        println!(
            "    eight-lane variant, reported and not scored: worst row {worst_ratio_swar8:.4}"
        );
        println!(
            "    target_feature_popcnt = {popcnt}; every arm makes zero count_ones calls per \
             cell, so no verdict here is contingent on the popcount lowering"
        );

        for r in &rows {
            let cells = r.cells as f64;
            run.record(&[
                ("field", r.field.to_string()),
                ("resolution", r.resolution.to_string()),
                ("scalar", r.scalar.to_string()),
                ("classify_share", format!("{:.6}", r.classify_share())),
                (
                    "instructions_per_cell_byte",
                    format!("{:.4}", r.instructions_per_cell_byte()),
                ),
                (
                    "instructions_per_cell_bitsliced",
                    format!("{:.4}", r.instructions_per_cell_bitsliced()),
                ),
                ("instruction_ratio", format!("{:.6}", r.instruction_ratio())),
                ("cases_identical", r.cases_identical.to_string()),
                ("cells", r.cells.to_string()),
                (
                    "corrupt_control_mismatches",
                    r.corrupt_control_mismatches.to_string(),
                ),
                ("hashes_moved", gold.moved.to_string()),
                ("c1_holds", r.c1_holds().to_string()),
                ("c2_holds", r.c2_holds().to_string()),
                ("c3_holds", c3_holds.to_string()),
                // ── extra columns (M-273) ──
                (
                    "classify_share_instructions",
                    format!("{:.6}", r.classify_share_instructions()),
                ),
                (
                    "cycles_classify_byte",
                    format!("{:.1}", r.byte_compact.cycles),
                ),
                ("cycles_extract_shipped", format!("{:.1}", r.extract.cycles)),
                (
                    "instructions_classify_byte",
                    format!("{:.1}", r.byte_compact.instructions),
                ),
                (
                    "instructions_extract_shipped",
                    format!("{:.1}", r.extract.instructions),
                ),
                (
                    "instructions_per_cell_byte_compact",
                    format!("{:.4}", r.byte_compact.instructions / cells),
                ),
                (
                    "instructions_per_cell_swar8",
                    format!("{:.4}", r.swar8_dense.instructions / cells),
                ),
                (
                    "instructions_per_cell_bitsliced_active",
                    format!("{:.4}", r.bitsliced_active.instructions / cells),
                ),
                (
                    "instruction_ratio_swar8",
                    format!("{:.6}", r.instruction_ratio_swar8()),
                ),
                (
                    "instruction_ratio_active",
                    format!("{:.6}", r.instruction_ratio_active()),
                ),
                (
                    "instruction_ratio_vs_byte_compact",
                    format!(
                        "{:.6}",
                        r.bitsliced_dense.instructions / r.byte_compact.instructions
                    ),
                ),
                (
                    "instruction_ratio_rep_spread",
                    format!("{:.6}", r.instruction_ratio_spread),
                ),
                (
                    "cycles_per_cell_byte",
                    format!("{:.4}", r.per_cell(|c| c.cycles, r.byte_dense)),
                ),
                (
                    "cycles_per_cell_bitsliced",
                    format!("{:.4}", r.per_cell(|c| c.cycles, r.bitsliced_dense)),
                ),
                ("cycle_ratio", format!("{:.6}", r.cycle_ratio())),
                ("cycle_ratio_rep_best", format!("{:.6}", r.cycle_ratio_best)),
                (
                    "cycle_ratio_rep_worst",
                    format!("{:.6}", r.cycle_ratio_worst),
                ),
                ("swar8_identical", r.swar8_identical.to_string()),
                (
                    "active_cases_identical",
                    r.active_cases_identical.to_string(),
                ),
                (
                    "single_bit_control_mismatches",
                    r.single_bit_control_mismatches.to_string(),
                ),
                ("sign_bits_identical", r.sign_bits_identical.to_string()),
                ("samples_checked", r.samples.to_string()),
                ("active_cells", r.active_cells.to_string()),
                (
                    "active_fraction",
                    format!("{:.6}", r.active_cells as f64 / cells),
                ),
                (
                    "sign_tests_per_cell_byte",
                    format!("{:.4}", 8.0 * cells / cells),
                ),
                (
                    "sign_tests_per_cell_bitsliced",
                    format!("{:.4}", r.samples as f64 / cells),
                ),
                ("inside_bytes", r.inside_bytes.to_string()),
                ("dense_cases_bytes", r.cells.to_string()),
                ("target_feature_popcnt", popcnt.to_string()),
                ("count_ones_calls_per_cell_byte", "0".to_string()),
                ("count_ones_calls_per_cell_bitsliced", "0".to_string()),
                ("inner_reps", r.inner.to_string()),
                ("inner_reps_extract", r.inner_extract.to_string()),
                (
                    "ns_per_cell_byte",
                    format!("{:.6}", r.byte_dense_ns / cells),
                ),
                (
                    "ns_per_cell_extract",
                    format!("{:.6}", r.extract_ns / cells),
                ),
                (
                    "ghz",
                    format!("{:.4}", r.byte_dense.cycles / r.byte_dense_ns),
                ),
                ("golden_hash_rows", gold.hash_rows.to_string()),
                ("golden_hashes_reproduced", gold.reproduced.to_string()),
                ("golden_cells_checked", gold.cells_checked.to_string()),
                ("golden_case_mismatches", gold.case_mismatches.to_string()),
                ("golden_samples_checked", gold.samples_checked.to_string()),
                ("golden_sign_mismatches", gold.sign_mismatches.to_string()),
                ("c1_rows_clearing_bar", c1_rows.to_string()),
                (
                    "c1_rows_clearing_bar_instructions",
                    c1_rows_instructions.to_string(),
                ),
                ("rows", rows.len().to_string()),
                (
                    "worst_instruction_ratio",
                    format!("{worst_instruction_ratio:.6}"),
                ),
                (
                    "best_instruction_ratio",
                    format!("{best_instruction_ratio:.6}"),
                ),
                (
                    "worst_instruction_ratio_swar8",
                    format!("{worst_ratio_swar8:.6}"),
                ),
                ("c1_holds_all_rows", c1_all.to_string()),
                ("c2_holds_all_rows", c2_all.to_string()),
            ]);
        }
    }
}

fn main() {
    if !std::env::args().any(|arg| arg == "--bench") {
        return;
    }

    let prereg = isomesh::experiment!("P-103");

    #[cfg(target_os = "linux")]
    common::experiment::run(prereg, experiment::run);

    // `experiment_p12`'s precedent, and the registration's own shape: C1 is a
    // share of extraction and C2 is an instruction rate. Off Linux there is no
    // `perf_event_open` and nothing to degrade to — a recorded zero would be a
    // fabricated share, and a fabricated share is what `✗51` cost.
    #[cfg(not(target_os = "linux"))]
    {
        eprintln!(
            "{} scores a share of extraction and an instruction rate with hardware performance \
             counters, and this platform has no `perf_event_open`. There is no clock \
             substitute: M-281 forbids a nanosecond carrying this verdict.",
            prereg.id
        );
        std::process::exit(1);
    }
}
