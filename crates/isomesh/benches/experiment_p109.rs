//! **P-109 — Elias–Fano for the edge→vertex structure, which is a dense flat vec and not a map.**
//!
//! Ticket: R-109. Pre-registered before this harness existed.
//!
//! ```bash
//! cargo bench --bench experiment_p109
//! ```
//!
//! Writes `docs/experiments/p-109.csv`. **Linux only**, for `experiment_p12`'s
//! reason: C2 is a cost ratio and its verdict is read off **instruction
//! counts**, which come from `perf_event_open`. On a governed CPU a nanosecond
//! is not a unit (`✗24`, `M-280`, `M-281`), so off Linux there is nothing to
//! degrade to and the harness exits 1 rather than record a fabricated zero.
//!
//! # What was missing, and three corrections against the research doc
//!
//! The sweep proposed Elias–Fano for "the edge→vertex map". Every clause below
//! is scored against the structure the crate actually has, which is not that.
//!
//! **(i) There is no edge→vertex MAP.** `grep -rn
//! 'edge_to_vertex\|edge_cache\|EdgeKey\|edge_id' crates/isomesh/src/` returns
//! nothing. The structure is `edge_vertices: Vec<u32>`, a **dense flat vec** of
//! `sample_count * 3` slots keyed by the arithmetic `lo_sample * 3 + axis`:
//! `marching_cubes/mod.rs:97` declares it, `:250-251` clears and resizes it to
//! `sample_count * 3` slots of `u32::MAX` on **every** `extract`, `:604`
//! computes the key, `:606-609` probes it, `:620-622` writes it, and `u32::MAX`
//! is the sentinel. The same shape appears at `marching_tetrahedra.rs:92` and
//! `property/extraction.rs:425`. So "access" means `edge_vertices[key]` — one
//! bounds-checked `u32` load into a 3.3 MB array at 65³ and a 25.8 MB array at
//! 129³ — and that, not a hash probe, is C2's incumbent.
//!
//! **(ii) The key side is already order-independent; what is discarded is the
//! VALUE.** The key is pure arithmetic on the grid, so it is derivable from the
//! grid alone and needs no encoding to become canonical. The value is a buffer
//! position handed out by the monotonic counter inside `MeshSink::vertex`, and
//! it depends on the order cells are walked. `M-318` already recorded that the
//! encoding is not the obstacle, that index-is-edge-id costs **230× memory**
//! (about 6.4 M slots for 27,822 vertices at 129³), and that the workable shape
//! carries state across extractions. This harness therefore encodes the **key
//! sequence** and carries the values in a parallel `Vec<u32>` — and reports both
//! costs, separately, because a key-side space figure is not the structure.
//!
//! **(iii) `V-45`'s stop applies to a structure that PERSISTS ACROSS CALLS.**
//! One **rebuilt from grid and field alone** is inside `V-45`'s own stated
//! reopening condition, which asks for *"a formulation where the persistent map
//! is derivable from the inputs"*. Everything here is rebuilt from `values` and
//! the grid, once per row, and it is bench-local: `crates/isomesh/src/**` is
//! read-only for Phase 25.
//!
//! # What is encoded, exactly, and what the values cost beside it
//!
//! Vigna, *Quasi-succinct indices* (`10.48550/arXiv.1206.4300`, in the corpus).
//! The encoded object is the **sorted sequence of crossing keys**, `n =
//! crossings` monotone values from the universe `u = 3 · sample_count`:
//!
//! - **Upper bits.** `keys[i] >> l` in a unary bitvector: bit `(keys[i] >> l) +
//!   i` is set, in `m = n + buckets` bits with `buckets = ((u − 1) >> l) + 1`.
//!   Because `l` is chosen so `n · 2^l ≥ u`, `buckets ≤ n` and the bitvector is
//!   **at most 2 bits per crossing** — the `2` in Vigna's `2 + ⌈log₂(u/n)⌉`.
//! - **Lower bits.** `keys[i] & (2^l − 1)`, `l = ⌈log₂(u/n)⌉` bits each, packed
//!   contiguously with no word alignment. Column `low_bits_width`.
//! - **Select hints.** The position of every [`ZERO_SAMPLE`]-th zero and every
//!   [`ONE_SAMPLE`]-th one, `u32` each. About **1 bit per crossing** for the
//!   pair; `select_hint_bytes` is a column and `ef_bytes_no_select` is beside
//!   `ef_bytes` so C1 can be scored on the paper's quantity as well as on the
//!   working structure.
//!
//! `ef_bytes` is upper + low + hints, and `bits_per_crossing` — the registered
//! column C1 is scored on — is `ef_bytes · 8 / crossings`. **That is a bound on
//! the KEY sequence alone.** The values are 32 bits each and are not encoded:
//! `value_bytes = 4 · crossings`, `bits_per_crossing_with_values` is
//! `(ef_bytes + value_bytes) · 8 / crossings`, and it is on every row so a
//! reader cannot mistake a key-side space figure for the whole structure. The
//! incumbent's figure is `dense_bits_per_crossing` = `dense_bytes · 8 /
//! crossings`, which is 32 bits per **slot** and so hundreds of bits per
//! **crossing** at a ~1.8% active fraction — that asymmetry is the whole reason
//! the mechanism was proposed.
//!
//! # The arithmetic C1 runs into, stated before it is measured
//!
//! `l = ⌈log₂(u/n)⌉` and the upper bits are `≥ 1` per crossing, so
//! `bits_per_crossing ≥ 1 + ⌈log₂(u/n)⌉`. A 4-bit budget therefore demands
//! `u/n ≤ 8`: at least one crossing per eight grid edges. The universe is
//! `3 n_axis³` and the crossing count is a surface, `O(n_axis²)`, so `u/n` grows
//! **linearly** in resolution. `universe_per_crossing` is a column and it is the
//! number that decides C1. This is not a reason to skip the clause — it is
//! measured, per field and per resolution, and a falsified C1 with the
//! fixture's own `u/n` beside it is the row's cheapest useful output.
//!
//! # The two different operations, and which one the paper makes O(1)
//!
//! This is the row's load-bearing structural finding, and it is the same shape
//! as `P-107`'s *"the directory buys random access, not throughput"*:
//!
//! - **`access(i)`** — the i-th key. `select1(i) − i` gives the high part,
//!   `low_at(i)` the low. **This is what Vigna's O(1) is about**, and it is
//!   measured as `ns_per_access_ef_positional` /
//!   `instructions_per_access_ef_positional`. The crate never performs it.
//! - **`lookup(key)`** — the value for a key, or the sentinel.
//!   `edge_vertices[key]` does this in one load; Elias–Fano has to *search*: one
//!   `select0(h − 1)` to find where bucket `h = key >> l` starts, a bit-run scan
//!   to where it ends, then a low-bits comparison per candidate. Bucket
//!   occupancy averages `n / buckets ≥ 1`, so the scan is short — but the
//!   `select0` is not free, and this is the operation C2 is denominated in
//!   because it is the operation the crate performs. `ns_per_access_ef` and
//!   `instructions_per_access_ef` are this one.
//!
//! # The query sets, and why the registered one is the crate's own probe order
//!
//! `MarchingCubes::vertex_on_edge` is called once per non-centroid triangle
//! corner, so the probe sequence is fixed by the cell walk and the `CASES`
//! table. Every probed key is a **cut** edge and therefore ends up present, so
//! the steady-state read is a hit in both arms. The registered
//! `ns_per_access_dense` / `ns_per_access_ef` / `access_ratio` are measured over
//! that exact sequence, `probes` accesses of it, replayed against the finished
//! structure.
//!
//! Beside it, [`RANDOM_QUERIES`] uniformly drawn present keys —
//! `*_random` columns — because the probe order has strong locality and the
//! dense array is the arm locality helps. Reporting only the sequential order
//! would flatter the incumbent; reporting only the random order would flatter
//! the encoding. Both are on the row.
//!
//! # Instructions carry C2's verdict
//!
//! `c2_holds` reads `access_ratio_instructions`. `access_ratio` (the registered
//! nanosecond form) is recorded and `c2_holds_ns` is beside it, but `M-280` and
//! `M-281` are why the verdict is not taken from a clock: `R-105` watched the
//! identical binary's cycle ratio band move from 0.984 to 1.035 across three
//! runs while its instruction counts held to four figures. `ghz` is on every
//! row with a nanosecond column.
//!
//! # The popcount this build does not have
//!
//! Elias–Fano's `select` is a popcount-heavy structure and every published
//! figure for it is priced against a one-cycle `popcnt`. **This build emits
//! none.** There is no `.cargo/config.toml` and no `target-cpu` anywhere in the
//! repository, so the `x86-64` baseline is in force and `u64::count_ones()`
//! lowers to the ~12-instruction SWAR sequence. `cfg!(target_feature =
//! "popcnt")` is false and is the column `target_feature_popcnt`.
//!
//! So the call count is a column too. `count_ones_per_access_ef` is **measured,
//! not estimated**: `select0`, `select1` and `select_in_word` are generic over a
//! `const TALLY: bool`, so the counting monomorphisation increments a tally at
//! every `count_ones` site and the measured monomorphisation compiles the tally
//! away entirely. One implementation, two instantiations, nothing added to the
//! hot path.
//!
//! `instructions_per_access_ef_popcnt_credited` subtracts
//! [`SWAR_POPCOUNT_INSTRUCTIONS`] per call — the SWAR sequence's cost less the
//! one instruction a `popcnt` would be — and
//! `access_ratio_instructions_popcnt_credited` is that ratio. `M-281` forbids
//! measuring the alternative by rebuilding with `target-cpu`, so the credit is
//! arithmetic on a measured call count rather than a second binary.
//! `c2_contingent_on_popcnt` is true only if C2 fails as built and would hold
//! under the credit; whether that is so is a result and not an assumption.
//!
//! # The mirror, and why its first job is to agree
//!
//! `edge_vertices` is private, so the map is rebuilt bench-local: the cell walk
//! of `marching_cubes/mod.rs:254-268`, `CASES[case]`, and the same monotonic
//! counter — including the `entry.centroids` vertices, which consume buffer
//! positions without owning a key (`marching_cubes/mod.rs:320-355`).
//!
//! A mirror is worth nothing on its own (`M-279`, and `R-120` and `R-121` both
//! caught real defects this way), so **every row asserts that the mirror
//! reproduces the shipped structure**, twice:
//!
//! - `crossings == shipped_vertices` with `centroid_vertices == 0`. Under the
//!   default `FaceAmbiguity::Separate` the derived table reaches no cycle long
//!   enough to need a centroid, so this is `✗1`/`M-2`/`M-22`'s `V_mc = C`
//!   recovered from the run rather than quoted at it.
//! - `mirror_on_edge == crossings`. For every key the mirror assigned, the
//!   shipped vertex at that buffer position is checked to **lie on the grid edge
//!   the key names**: its two off-axis coordinates are compared *exactly*
//!   against `origin[a] + cell_size · coord`, which is what `place(lo, lo, d) =
//!   (lo + lo)·½ + 0·d` returns bit for bit, and its axis coordinate is checked
//!   to lie inside the segment. Equal vertex *counts* would not identify the
//!   map; this does.
//!
//! # SHARE
//!
//! **This row prices an encoding and moves no extraction time.** C1 is a space
//! bound and C2 a cost bound against a named baseline, so neither is a ratio of
//! a total and `✗51`'s `1/(1 − share/factor)` bar does not apply to either.
//! What the row feeds is `R-027a`'s 45× ceiling, which `M-318` measured and
//! which nothing here re-claims. Each clause's reachable share is a column:
//!
//! - **C1's share is `bits_per_crossing` itself**, bar 4.0, against the key
//!   sequence. `universe_per_crossing` is the fixture quantity that bounds it
//!   from below, `ef_bytes_no_select` scores it without the select hints, and
//!   `bits_per_crossing_with_values` scores the whole structure. Four numbers,
//!   one clause, no re-derivation needed.
//! - **C2's share is `access_ratio_instructions`**, bar 1.5, over the crate's
//!   own probe sequence. `access_ratio` (ns), `access_ratio_random`,
//!   `access_ratio_instructions_random` and the positional pair are beside it.
//! - **C3's share is `values_equal / crossings`, and the bar is 1.** An
//!   equality over an enumerated population, so its denominator is exact by
//!   construction. It moves no time and has no ceiling.
//!
//! **VACUITY CONTROL, asserted rather than recorded.** `crossings` non-zero on
//! every row; `values_equal == crossings`; `positional_values_equal ==
//! crossings`; `absent_sentinel_correct == absent_probed`, so a lookup that
//! answered *something* for every key could not pass; and
//! **`mutant_value_mismatches > 0`**, from a decoder that reassembles keys with
//! the low-bits width **shifted by one** (`mutant_low_bits_width = l + 1`)
//! against a structure encoded at `l`. Without that last one C3 is an equality
//! between two names for the same computation: the comparand *is* the dense
//! array the mirror built, so a comparator that cannot fail proves nothing.

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

    use isomesh::marching_cubes::MarchingCubes;
    use isomesh::marching_cubes::table::{CASES, EDGE_AXIS, EDGE_CORNERS, is_centroid, is_inside};
    use isomesh::{MeshBuffer, RuntimeShape3, Sdf};

    use crate::common::counters::{MIN_TIME_RATIO, Probe};

    /// The registered fixture: eight reference fields at each of these.
    const RESOLUTIONS: [u32; 3] = [33, 65, 129];

    /// Repetitions per window, medianed per quantity.
    const REPS: usize = 5;

    /// Passes discarded before any window opens.
    const WARMUP: usize = 2;

    /// About this long per counter window, so the ~28 `perf_event` system calls
    /// a window costs land outside it and cannot inflate anything.
    const TARGET_BATCH_NS: f64 = 20_000_000.0;

    /// Ceiling on the batch, so a cheap pass cannot run for a minute.
    const MAX_INNER: usize = 4096;

    /// Uniformly drawn present keys for the random-access arm.
    pub(crate) const RANDOM_QUERIES: usize = 4096;

    /// C1's bar: bits per crossing on the encoded key sequence.
    const BITS_PER_CROSSING_BAR: f64 = 4.0;

    /// C2's bar: access cost against direct addressing.
    const ACCESS_RATIO_BAR: f64 = 1.5;

    /// Zeros between select-0 hint entries.
    pub(crate) const ZERO_SAMPLE: usize = 64;

    /// Ones between select-1 hint entries.
    pub(crate) const ONE_SAMPLE: usize = 64;

    /// Instructions a `count_ones` costs here **beyond** what a `popcnt` would.
    ///
    /// The SWAR sequence LLVM emits for the `x86-64` baseline is about twelve
    /// instructions and `popcnt` is one, so a hardware popcount would remove
    /// about eleven per call. Used only for
    /// `instructions_per_access_ef_popcnt_credited`, and it is an arithmetic
    /// credit on a measured call count — never a second binary, which `M-281`
    /// forbids comparing against.
    pub(crate) const SWAR_POPCOUNT_INSTRUCTIONS: f64 = 11.0;

    /// `cube.rs:149`'s `corner_offset`, which is `pub(crate)`.
    ///
    /// Copied verbatim rather than made public: `crates/isomesh/src/**` is
    /// read-only for Phase 25.
    #[inline]
    const fn corner_offset(corner: u8) -> [usize; 3] {
        [
            (corner & 1) as usize,
            ((corner >> 1) & 1) as usize,
            ((corner >> 2) & 1) as usize,
        ]
    }

    // ─── the dense structure, mirrored from `marching_cubes/mod.rs` ─────────

    /// The mirror of `edge_vertices`, plus the probe sequence that reads it.
    struct Dense {
        /// `sample_count * 3` slots, `u32::MAX` where no vertex sits.
        slots: Vec<u32>,
        /// Keys in the order `vertex_on_edge` is called on them.
        probes: Vec<u32>,
        /// Buffer positions handed out in total, edge and centroid alike.
        vertices: u32,
        /// Cycle-centroid vertices, which own no key.
        centroids: u32,
    }

    impl Dense {
        /// Walk cells exactly as `marching_cubes/mod.rs:254-377` does.
        ///
        /// `values` is row-major with stride `n`, which is the stride
        /// `MarchingCubes::extract` samples on (`sdf::sample_grid` with
        /// `size[0]` as the row stride).
        fn build(values: &[f32], n: u32) -> Self {
            let sx = n as usize;
            let cells = sx - 1;
            let mut slots = vec![u32::MAX; sx * sx * sx * 3];
            let mut probes = Vec::new();
            let mut next = 0u32;
            let mut centroids = 0u32;

            for z in 0..cells {
                for y in 0..cells {
                    for x in 0..cells {
                        let mut case = 0u8;
                        for c in 0..8u8 {
                            let o = corner_offset(c);
                            let s = (z + o[2]) * sx * sx + (y + o[1]) * sx + (x + o[0]);
                            if is_inside(values[s]) {
                                case |= 1 << c;
                            }
                        }
                        let entry = CASES[case as usize];
                        if entry.count == 0 {
                            continue;
                        }
                        // Cycle centroids are emitted first and are cell-local:
                        // they consume a buffer position and own no key, so the
                        // counter has to see them or every later value is off.
                        next += u32::from(entry.centroids);
                        centroids += u32::from(entry.centroids);

                        for tri in &entry.triangles[..entry.count as usize] {
                            for &code in tri {
                                if is_centroid(code) {
                                    continue;
                                }
                                let axis = EDGE_AXIS[code as usize] as usize;
                                let lo_corner = EDGE_CORNERS[code as usize][0];
                                let o = corner_offset(lo_corner);
                                let lo_sample = (z + o[2]) * sx * sx + (y + o[1]) * sx + (x + o[0]);
                                let key = lo_sample * 3 + axis;
                                probes.push(key as u32);
                                if slots[key] == u32::MAX {
                                    slots[key] = next;
                                    next += 1;
                                }
                            }
                        }
                    }
                }
            }

            Self {
                slots,
                probes,
                vertices: next,
                centroids,
            }
        }
    }

    /// Crossings whose shipped vertex lies on the grid edge its key names.
    ///
    /// The identification, not a count comparison. `place(lo, hi, d) = (lo +
    /// hi)·½ + (hi − lo)·d`, so on an off-axis coordinate where `lo == hi` it
    /// returns `lo` bit for bit — which makes the two off-axis coordinates of an
    /// edge vertex *exactly* `origin[a] + cell_size · coord`.
    #[allow(
        clippy::float_cmp,
        reason = "an approximate comparison would not identify the map: the off-axis \
                  coordinates of an edge vertex are exact by construction"
    )]
    fn mirror_on_edge(
        dense: &Dense,
        positions: &[[f32; 3]],
        n: u32,
        origin: [f32; 3],
        cell_size: f32,
    ) -> usize {
        let sx = n as usize;
        let mut agreeing = 0usize;
        for (key, &v) in dense.slots.iter().enumerate() {
            if v == u32::MAX {
                continue;
            }
            let position = positions[v as usize];
            let axis = key % 3;
            let lo_sample = key / 3;
            let coord = [lo_sample % sx, (lo_sample / sx) % sx, lo_sample / (sx * sx)];
            let mut ok = true;
            for (a, &c) in coord.iter().enumerate() {
                let lo = origin[a] + cell_size * c as f32;
                if a == axis {
                    let hi = origin[a] + cell_size * (c + 1) as f32;
                    if position[a] < lo || position[a] > hi {
                        ok = false;
                    }
                } else if position[a] != lo {
                    ok = false;
                }
            }
            if ok {
                agreeing += 1;
            }
        }
        agreeing
    }

    // ─── Elias–Fano ────────────────────────────────────────────────────────

    /// Position of the `k`-th set bit of `word`, `k` zero-based.
    ///
    /// Byte at a time, so the cost is bounded at eight `count_ones` and seven
    /// `blsr`-shaped clears rather than the 63 a bare clear loop would take.
    /// Generic over `TALLY` so the `count_ones` sites can be counted without
    /// putting a counter on the measured path.
    #[inline]
    fn select_in_word<const TALLY: bool>(word: u64, k: u32, tally: &mut u64) -> usize {
        let mut w = word;
        let mut left = k;
        let mut base = 0usize;
        for _ in 0..8 {
            if TALLY {
                *tally += 1;
            }
            let count = (w & 0xFF).count_ones();
            if count > left {
                break;
            }
            left -= count;
            w >>= 8;
            base += 8;
        }
        debug_assert!(base < 64, "select_in_word ran off the word");
        let mut byte = w & 0xFF;
        for _ in 0..left {
            byte &= byte - 1;
        }
        base + byte.trailing_zeros() as usize
    }

    /// Vigna's quasi-succinct encoding of a monotone key sequence, with the
    /// buffer positions carried beside it.
    ///
    /// The keys are encoded; the values are not, and `Self::value_bytes` says
    /// what they cost. `Self::lookup` is the operation `edge_vertices[key]`
    /// performs; `Self::key_at` is the O(1) positional access the paper is
    /// about, which the crate never asks for.
    struct EliasFano {
        /// `m = n + buckets` bits, `n` of them set.
        upper: Vec<u64>,
        upper_bits: usize,
        /// `n * low_width` bits, packed with no alignment. One spare word so a
        /// straddling read never indexes past the end.
        low: Vec<u64>,
        low_width: u32,
        low_mask: u64,
        buckets: usize,
        /// Position of every [`ZERO_SAMPLE`]-th zero of `upper`.
        zero_hint: Vec<u32>,
        /// Position of every [`ONE_SAMPLE`]-th one of `upper`.
        one_hint: Vec<u32>,
        /// Buffer positions, in sorted-key order.
        values: Vec<u32>,
    }

    impl EliasFano {
        /// `⌈log₂(u/n)⌉`, as the smallest `l` with `n · 2^l ≥ u`.
        fn low_bits_width(n: usize, universe: usize) -> u32 {
            let mut width = 0u32;
            while (n as u128) << width < universe as u128 {
                width += 1;
            }
            width
        }

        /// Encode `keys`, which must be strictly ascending and below `universe`.
        fn build(keys: &[u32], values: &[u32], universe: usize) -> Self {
            let n = keys.len();
            assert!(n > 0, "an empty key sequence has no encoding to price");
            let low_width = Self::low_bits_width(n, universe);
            let buckets = ((universe - 1) >> low_width) + 1;
            let upper_bits = n + buckets;

            let mut upper = vec![0u64; upper_bits.div_ceil(64)];
            for (i, &key) in keys.iter().enumerate() {
                let position = ((key as usize) >> low_width) + i;
                upper[position >> 6] |= 1u64 << (position & 63);
            }

            let low_mask = if low_width == 0 {
                0
            } else {
                (1u64 << low_width) - 1
            };
            let mut low = vec![0u64; (n * low_width as usize).div_ceil(64) + 1];
            if low_width > 0 {
                for (i, &key) in keys.iter().enumerate() {
                    let value = u64::from(key) & low_mask;
                    let bit = i * low_width as usize;
                    let word = bit >> 6;
                    let offset = bit & 63;
                    low[word] |= value << offset;
                    if offset + low_width as usize > 64 {
                        low[word + 1] |= value >> (64 - offset);
                    }
                }
            }

            let mut zero_hint = Vec::with_capacity(buckets.div_ceil(ZERO_SAMPLE));
            let mut one_hint = Vec::with_capacity(n.div_ceil(ONE_SAMPLE));
            let mut zeros = 0usize;
            let mut ones = 0usize;
            for position in 0..upper_bits {
                if upper[position >> 6] & (1u64 << (position & 63)) == 0 {
                    if zeros.is_multiple_of(ZERO_SAMPLE) {
                        zero_hint.push(position as u32);
                    }
                    zeros += 1;
                } else {
                    if ones.is_multiple_of(ONE_SAMPLE) {
                        one_hint.push(position as u32);
                    }
                    ones += 1;
                }
            }
            assert_eq!(
                zeros, buckets,
                "the upper bitvector has the wrong zero count"
            );
            assert_eq!(ones, n, "the upper bitvector has the wrong one count");

            Self {
                upper,
                upper_bits,
                low,
                low_width,
                low_mask,
                buckets,
                zero_hint,
                one_hint,
                values: values.to_vec(),
            }
        }

        /// The `low_width` bits stored for element `i`.
        #[inline]
        fn low_at(&self, i: usize) -> u64 {
            if self.low_width == 0 {
                return 0;
            }
            let bit = i * self.low_width as usize;
            let word = bit >> 6;
            let offset = bit & 63;
            let mut value = self.low[word] >> offset;
            if offset + self.low_width as usize > 64 {
                value |= self.low[word + 1] << (64 - offset);
            }
            value & self.low_mask
        }

        /// Position of the `h`-th zero of `upper`, `h` zero-based.
        #[inline]
        fn select0<const TALLY: bool>(&self, h: usize, tally: &mut u64) -> usize {
            let sampled = h / ZERO_SAMPLE;
            let mut position = self.zero_hint[sampled] as usize;
            let mut left = h - sampled * ZERO_SAMPLE;
            if left == 0 {
                return position;
            }
            position += 1;
            loop {
                let word = position >> 6;
                let offset = position & 63;
                let tail = !self.upper[word] >> offset;
                if TALLY {
                    *tally += 1;
                }
                let count = tail.count_ones() as usize;
                if count >= left {
                    return position + select_in_word::<TALLY>(tail, left as u32 - 1, tally);
                }
                left -= count;
                position = (word + 1) << 6;
            }
        }

        /// Position of the `i`-th one of `upper`, `i` zero-based.
        #[inline]
        fn select1<const TALLY: bool>(&self, i: usize, tally: &mut u64) -> usize {
            let sampled = i / ONE_SAMPLE;
            let mut position = self.one_hint[sampled] as usize;
            let mut left = i - sampled * ONE_SAMPLE;
            if left == 0 {
                return position;
            }
            position += 1;
            loop {
                let word = position >> 6;
                let offset = position & 63;
                let tail = self.upper[word] >> offset;
                if TALLY {
                    *tally += 1;
                }
                let count = tail.count_ones() as usize;
                if count >= left {
                    return position + select_in_word::<TALLY>(tail, left as u32 - 1, tally);
                }
                left -= count;
                position = (word + 1) << 6;
            }
        }

        /// Index range of the elements whose high part is `h`.
        ///
        /// One `select0` for where the bucket starts, then a bit-run scan for
        /// where it ends. The last bit of `upper` is always a zero — every
        /// bucket is terminated by its own separator — so the scan cannot run
        /// off the end.
        #[inline]
        fn bucket<const TALLY: bool>(&self, h: usize, tally: &mut u64) -> (usize, usize) {
            if h >= self.buckets {
                return (0, 0);
            }
            let start = if h == 0 {
                0
            } else {
                self.select0::<TALLY>(h - 1, tally) + 1
            };
            let mut position = start;
            while position < self.upper_bits {
                let offset = position & 63;
                let word = self.upper[position >> 6] >> offset;
                let run = (!word).trailing_zeros() as usize;
                let available = 64 - offset;
                if run < available {
                    position += run;
                    break;
                }
                position += available;
            }
            (start - h, position - h)
        }

        /// The value stored for `key`, or `u32::MAX`.
        ///
        /// `width` is the low-bits width the *decoder* uses. It is
        /// `self.low_width` on every real call; the mutant passes a shifted one,
        /// which is what makes C3's comparator able to fail.
        #[inline]
        fn lookup<const TALLY: bool>(&self, key: u32, width: u32, tally: &mut u64) -> u32 {
            let h = (key as usize) >> width;
            let mask = if width >= 64 {
                u64::MAX
            } else {
                (1u64 << width) - 1
            };
            let low = u64::from(key) & mask;
            let (first, last) = self.bucket::<TALLY>(h, tally);
            for i in first..last {
                if self.low_at(i) == low {
                    return self.values[i];
                }
            }
            u32::MAX
        }

        /// The `i`-th key. The O(1) positional access the paper is about.
        #[inline]
        fn key_at<const TALLY: bool>(&self, i: usize, tally: &mut u64) -> u32 {
            let high = self.select1::<TALLY>(i, tally) - i;
            ((high << self.low_width) as u64 | self.low_at(i)) as u32
        }

        fn upper_bytes(&self) -> usize {
            self.upper.len() * size_of::<u64>()
        }

        fn low_bytes(&self) -> usize {
            self.low.len() * size_of::<u64>()
        }

        fn hint_bytes(&self) -> usize {
            (self.zero_hint.len() + self.one_hint.len()) * size_of::<u32>()
        }

        fn value_bytes(&self) -> usize {
            self.values.len() * size_of::<u32>()
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

    /// One counter window over `inner` repetitions, divided by `inner`.
    ///
    /// Windows are **siblings, never nested**: Zen 3 has six general-purpose
    /// counters and `Probe` opens exactly six, so a nested window multiplexes
    /// and `Counts::worst_ratio` refuses. `R-121` paid for that discovery.
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

    /// `count` deterministic values below `bound`, splitmix64 so the set is the
    /// same in every run and every build.
    fn draw(count: usize, bound: usize) -> Vec<usize> {
        let mut state = 0x9E37_79B9_7F4A_7C15u64;
        (0..count)
            .map(|_| {
                state = state.wrapping_add(0x9E37_79B9_7F4A_7C15);
                let mut z = state;
                z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
                z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
                z ^= z >> 31;
                (z % bound as u64) as usize
            })
            .collect()
    }

    // ─── one row ───────────────────────────────────────────────────────────

    struct Row {
        field: &'static str,
        resolution: u32,
        samples: usize,
        universe: usize,
        crossings: usize,
        probes: usize,
        universe_per_crossing: f64,
        low_bits_width: u32,
        upper_bits: usize,
        buckets: usize,
        // space
        dense_bytes: usize,
        ef_bytes: usize,
        ef_bytes_no_select: usize,
        upper_bytes: usize,
        low_bytes: usize,
        select_hint_bytes: usize,
        value_bytes: usize,
        bits_per_crossing: f64,
        bits_per_crossing_no_select: f64,
        bits_per_crossing_with_values: f64,
        dense_bits_per_crossing: f64,
        space_ratio_dense_over_ef: f64,
        // access, probe order
        ns_per_access_dense: f64,
        ns_per_access_ef: f64,
        access_ratio: f64,
        instructions_per_access_dense: f64,
        instructions_per_access_ef: f64,
        access_ratio_instructions: f64,
        // access, random present keys
        ns_per_access_dense_random: f64,
        ns_per_access_ef_random: f64,
        access_ratio_random: f64,
        instructions_per_access_dense_random: f64,
        instructions_per_access_ef_random: f64,
        access_ratio_instructions_random: f64,
        // the operation the paper makes O(1)
        ns_per_access_ef_positional: f64,
        instructions_per_access_ef_positional: f64,
        // popcount
        count_ones_per_access_ef: f64,
        count_ones_per_access_ef_positional: f64,
        target_feature_popcnt: bool,
        instructions_per_access_ef_popcnt_credited: f64,
        access_ratio_instructions_popcnt_credited: f64,
        c2_contingent_on_popcnt: bool,
        // correctness
        values_equal: usize,
        positional_values_equal: usize,
        absent_probed: usize,
        absent_sentinel_correct: usize,
        mutant_value_mismatches: usize,
        mutant_low_bits_width: u32,
        // the mirror
        shipped_vertices: usize,
        mirror_vertices: usize,
        centroid_vertices: u32,
        mirror_on_edge: usize,
        // provenance
        ghz: f64,
        inner_probe: usize,
        inner_random: usize,
        // verdicts
        c1_holds: bool,
        c1_holds_no_select: bool,
        c1_holds_with_values: bool,
        c2_holds: bool,
        c2_holds_ns: bool,
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

        // ── the mirror, and the shipped structure it has to reproduce ────────
        let dense = Dense::build(&values, n);
        let mut mc = MarchingCubes::<f32>::new();
        let mut out = MeshBuffer::<f32>::new();
        mc.extract(sdf, &shape, origin, cell_size, &mut out)
            .expect("extraction");

        let mut keys = Vec::new();
        let mut slot_values = Vec::new();
        for (key, &v) in dense.slots.iter().enumerate() {
            if v != u32::MAX {
                keys.push(key as u32);
                slot_values.push(v);
            }
        }
        let crossings = keys.len();
        assert!(
            crossings > 0,
            "{field} {n}^3: no crossings, so every clause would be scored over an empty set"
        );

        let shipped_vertices = out.vertex_count();
        assert_eq!(
            dense.vertices as usize, shipped_vertices,
            "{field} {n}^3: the mirror handed out {} buffer positions and the shipped \
             extractor handed out {shipped_vertices}, so the mirror is not the structure",
            dense.vertices
        );
        assert_eq!(
            crossings, shipped_vertices,
            "{field} {n}^3: {crossings} crossings against {shipped_vertices} shipped \
             vertices, so `V_mc = C` does not hold and the population is wrong"
        );
        let on_edge = mirror_on_edge(&dense, &out.positions, n, origin, cell_size);
        assert_eq!(
            on_edge,
            crossings,
            "{field} {n}^3: {} of {crossings} shipped vertices are not on the grid edge \
             the mirror's key names, so the mirror's key -> value map is not `edge_vertices`",
            crossings - on_edge
        );

        // ── the encoding ─────────────────────────────────────────────────────
        let universe = dense.slots.len();
        let ef = EliasFano::build(&keys, &slot_values, universe);
        let low_width = ef.low_width;
        let mutant_width = low_width + 1;
        assert!(
            keys.windows(2).all(|w| w[0] < w[1]),
            "{field} {n}^3: the key sequence is not strictly ascending, so Elias-Fano \
             does not apply to it"
        );

        // ── C3, and its controls ─────────────────────────────────────────────
        let mut tally = 0u64;
        let mut values_equal = 0usize;
        let mut positional_values_equal = 0usize;
        let mut mutant_value_mismatches = 0usize;
        for (i, &key) in keys.iter().enumerate() {
            let dense_value = dense.slots[key as usize];
            if ef.lookup::<false>(key, low_width, &mut tally) == dense_value {
                values_equal += 1;
            }
            if ef.key_at::<false>(i, &mut tally) == key && ef.values[i] == dense_value {
                positional_values_equal += 1;
            }
            if ef.lookup::<false>(key, mutant_width, &mut tally) != dense_value {
                mutant_value_mismatches += 1;
            }
        }
        assert_eq!(
            values_equal,
            crossings,
            "{field} {n}^3: the decoded map differs from the dense array on {} of \
             {crossings} crossings",
            crossings - values_equal
        );
        assert_eq!(
            positional_values_equal,
            crossings,
            "{field} {n}^3: positional access disagrees on {} of {crossings} crossings",
            crossings - positional_values_equal
        );
        assert!(
            mutant_value_mismatches > 0,
            "{field} {n}^3: a decoder using low-bits width {mutant_width} against a \
             structure encoded at {low_width} still answered every crossing correctly, \
             so C3's comparator cannot fail and the equality is vacuous"
        );

        // A lookup that answered *something* for every key would pass C3, so the
        // sentinel is asserted over keys the structure does not hold.
        let mut absent_probed = 0usize;
        let mut absent_sentinel_correct = 0usize;
        for key in draw(RANDOM_QUERIES, universe) {
            if dense.slots[key] == u32::MAX {
                absent_probed += 1;
                if ef.lookup::<false>(key as u32, low_width, &mut tally) == u32::MAX {
                    absent_sentinel_correct += 1;
                }
            }
        }
        assert!(
            absent_probed > 0,
            "{field} {n}^3: every drawn key was present, so the sentinel is untested"
        );
        assert_eq!(
            absent_sentinel_correct,
            absent_probed,
            "{field} {n}^3: the encoding answered a value for {} of {absent_probed} \
             absent keys",
            absent_probed - absent_sentinel_correct
        );

        // ── the popcount call counts, measured rather than estimated ─────────
        let mut lookup_tally = 0u64;
        for &key in &dense.probes {
            black_box(ef.lookup::<true>(key, low_width, &mut lookup_tally));
        }
        let count_ones_per_access_ef = lookup_tally as f64 / dense.probes.len() as f64;

        let random_indices = draw(RANDOM_QUERIES, crossings);
        let random_keys: Vec<u32> = random_indices.iter().map(|&i| keys[i]).collect();
        let mut positional_tally = 0u64;
        for &i in &random_indices {
            black_box(ef.key_at::<true>(i, &mut positional_tally));
        }
        let count_ones_per_access_ef_positional =
            positional_tally as f64 / random_indices.len() as f64;

        // ── the windows ──────────────────────────────────────────────────────
        let mut sink = 0u64;
        let probe_dense = |sink: &mut u64| {
            let mut acc = 0u32;
            for &key in &dense.probes {
                acc = acc.wrapping_add(dense.slots[key as usize]);
            }
            *sink = sink.wrapping_add(u64::from(acc));
        };
        let probe_ef = |sink: &mut u64| {
            let mut acc = 0u32;
            let mut ignored = 0u64;
            for &key in &dense.probes {
                acc = acc.wrapping_add(ef.lookup::<false>(key, low_width, &mut ignored));
            }
            *sink = sink.wrapping_add(u64::from(acc));
        };
        let random_dense = |sink: &mut u64| {
            let mut acc = 0u32;
            for &key in &random_keys {
                acc = acc.wrapping_add(dense.slots[key as usize]);
            }
            *sink = sink.wrapping_add(u64::from(acc));
        };
        let random_ef = |sink: &mut u64| {
            let mut acc = 0u32;
            let mut ignored = 0u64;
            for &key in &random_keys {
                acc = acc.wrapping_add(ef.lookup::<false>(key, low_width, &mut ignored));
            }
            *sink = sink.wrapping_add(u64::from(acc));
        };
        let positional_ef = |sink: &mut u64| {
            let mut acc = 0u32;
            let mut ignored = 0u64;
            for &i in &random_indices {
                acc = acc.wrapping_add(ef.key_at::<false>(i, &mut ignored));
            }
            *sink = sink.wrapping_add(u64::from(acc));
        };

        for _ in 0..WARMUP {
            probe_dense(&mut sink);
            probe_ef(&mut sink);
            random_dense(&mut sink);
            random_ef(&mut sink);
            positional_ef(&mut sink);
        }

        let inner_probe = choose_inner(|| probe_dense(&mut sink));
        let inner_random = choose_inner(|| random_dense(&mut sink));

        let mut probe = Probe::open();
        let mut dense_probe = Vec::with_capacity(REPS);
        let mut ef_probe = Vec::with_capacity(REPS);
        let mut dense_random = Vec::with_capacity(REPS);
        let mut ef_random = Vec::with_capacity(REPS);
        let mut ef_positional = Vec::with_capacity(REPS);

        for _ in 0..REPS {
            dense_probe.push(window(&mut probe, inner_probe, || probe_dense(&mut sink)));
            ef_probe.push(window(&mut probe, inner_probe, || probe_ef(&mut sink)));
            dense_random.push(window(&mut probe, inner_random, || random_dense(&mut sink)));
            ef_random.push(window(&mut probe, inner_random, || random_ef(&mut sink)));
            ef_positional.push(window(&mut probe, inner_random, || {
                positional_ef(&mut sink);
            }));
        }
        black_box(sink);

        let dense_probe = median_of(&dense_probe);
        let ef_probe = median_of(&ef_probe);
        let dense_random = median_of(&dense_random);
        let ef_random = median_of(&ef_random);
        let ef_positional = median_of(&ef_positional);

        let per_probe = 1.0 / dense.probes.len() as f64;
        let per_random = 1.0 / RANDOM_QUERIES as f64;

        // ── space ────────────────────────────────────────────────────────────
        let dense_bytes = universe * size_of::<u32>();
        let ef_bytes_no_select = ef.upper_bytes() + ef.low_bytes();
        let ef_bytes = ef_bytes_no_select + ef.hint_bytes();
        let value_bytes = ef.value_bytes();
        let bits = |bytes: usize| bytes as f64 * 8.0 / crossings as f64;

        let bits_per_crossing = bits(ef_bytes);
        let instructions_per_access_ef = ef_probe.instructions * per_probe;
        let instructions_per_access_dense = dense_probe.instructions * per_probe;
        let access_ratio_instructions =
            instructions_per_access_ef / instructions_per_access_dense.max(f64::MIN_POSITIVE);
        let credited = (instructions_per_access_ef
            - SWAR_POPCOUNT_INSTRUCTIONS * count_ones_per_access_ef)
            .max(0.0);
        let credited_ratio = credited / instructions_per_access_dense.max(f64::MIN_POSITIVE);

        let c1_holds = bits_per_crossing <= BITS_PER_CROSSING_BAR;
        let c2_holds = access_ratio_instructions <= ACCESS_RATIO_BAR;
        let c3_holds = values_equal == crossings
            && positional_values_equal == crossings
            && absent_sentinel_correct == absent_probed
            && mutant_value_mismatches > 0;

        Row {
            field,
            resolution: n,
            samples,
            universe,
            crossings,
            probes: dense.probes.len(),
            universe_per_crossing: universe as f64 / crossings as f64,
            low_bits_width: low_width,
            upper_bits: ef.upper_bits,
            buckets: ef.buckets,
            dense_bytes,
            ef_bytes,
            ef_bytes_no_select,
            upper_bytes: ef.upper_bytes(),
            low_bytes: ef.low_bytes(),
            select_hint_bytes: ef.hint_bytes(),
            value_bytes,
            bits_per_crossing,
            bits_per_crossing_no_select: bits(ef_bytes_no_select),
            bits_per_crossing_with_values: bits(ef_bytes + value_bytes),
            dense_bits_per_crossing: bits(dense_bytes),
            space_ratio_dense_over_ef: dense_bytes as f64 / (ef_bytes + value_bytes) as f64,
            ns_per_access_dense: dense_probe.nanos * per_probe,
            ns_per_access_ef: ef_probe.nanos * per_probe,
            access_ratio: ef_probe.nanos / dense_probe.nanos.max(f64::MIN_POSITIVE),
            instructions_per_access_dense,
            instructions_per_access_ef,
            access_ratio_instructions,
            ns_per_access_dense_random: dense_random.nanos * per_random,
            ns_per_access_ef_random: ef_random.nanos * per_random,
            access_ratio_random: ef_random.nanos / dense_random.nanos.max(f64::MIN_POSITIVE),
            instructions_per_access_dense_random: dense_random.instructions * per_random,
            instructions_per_access_ef_random: ef_random.instructions * per_random,
            access_ratio_instructions_random: ef_random.instructions
                / dense_random.instructions.max(f64::MIN_POSITIVE),
            ns_per_access_ef_positional: ef_positional.nanos * per_random,
            instructions_per_access_ef_positional: ef_positional.instructions * per_random,
            count_ones_per_access_ef,
            count_ones_per_access_ef_positional,
            target_feature_popcnt: cfg!(target_feature = "popcnt"),
            instructions_per_access_ef_popcnt_credited: credited,
            access_ratio_instructions_popcnt_credited: credited_ratio,
            c2_contingent_on_popcnt: !c2_holds && credited_ratio <= ACCESS_RATIO_BAR,
            values_equal,
            positional_values_equal,
            absent_probed,
            absent_sentinel_correct,
            mutant_value_mismatches,
            mutant_low_bits_width: mutant_width,
            shipped_vertices,
            mirror_vertices: dense.vertices as usize,
            centroid_vertices: dense.centroids,
            mirror_on_edge: on_edge,
            ghz: dense_probe.cycles / dense_probe.nanos.max(f64::MIN_POSITIVE),
            inner_probe,
            inner_random,
            c1_holds,
            c1_holds_no_select: bits(ef_bytes_no_select) <= BITS_PER_CROSSING_BAR,
            c1_holds_with_values: bits(ef_bytes + value_bytes) <= BITS_PER_CROSSING_BAR,
            c2_holds,
            c2_holds_ns: ef_probe.nanos / dense_probe.nanos.max(f64::MIN_POSITIVE)
                <= ACCESS_RATIO_BAR,
            c3_holds,
        }
    }

    /// Eight reference fields × three resolutions, `f32`.
    ///
    /// No `scalar` column is registered and none is added: the key sequence is
    /// integer arithmetic on the grid and the values are `u32` buffer positions,
    /// so nothing in either arm changes with the field's precision. An `f64` arm
    /// would double the fixture to say so.
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
        let c2 = rows.iter().filter(|r| r.c2_holds).count();
        let best_bits = rows
            .iter()
            .map(|r| r.bits_per_crossing)
            .fold(f64::INFINITY, f64::min);
        let best_ratio = rows
            .iter()
            .map(|r| r.access_ratio_instructions)
            .fold(f64::INFINITY, f64::min);
        let worst_u_over_n = rows
            .iter()
            .map(|r| r.universe_per_crossing)
            .fold(0.0f64, f64::max);

        println!(
            "P-109: C1 (<= {BITS_PER_CROSSING_BAR:.1} bits per crossing on the key sequence) \
             holds on {c1} of {} rows; the cheapest encoding measured is {best_bits:.3} bits, \
             and u/n reaches {worst_u_over_n:.1} -- a 4-bit budget needs u/n <= 8",
            rows.len()
        );
        println!(
            "P-109: C2 (<= {ACCESS_RATIO_BAR:.1}x direct addressing, instruction form) holds on \
             {c2} of {} rows; the best instruction ratio measured is {best_ratio:.3}",
            rows.len()
        );
        println!(
            "P-109: target_feature_popcnt is {}, so a `count_ones` is the SWAR sequence and not \
             an instruction; an Elias-Fano lookup makes {:.2} of them and a positional access \
             {:.2}. C2 is contingent on it on {} rows",
            cfg!(target_feature = "popcnt"),
            rows[0].count_ones_per_access_ef,
            rows[0].count_ones_per_access_ef_positional,
            rows.iter().filter(|r| r.c2_contingent_on_popcnt).count()
        );

        for row in &rows {
            run.record(&[
                ("field", row.field.to_string()),
                ("resolution", row.resolution.to_string()),
                ("samples", row.samples.to_string()),
                ("universe", row.universe.to_string()),
                ("crossings", row.crossings.to_string()),
                ("probes", row.probes.to_string()),
                (
                    "universe_per_crossing",
                    format!("{:.4}", row.universe_per_crossing),
                ),
                ("low_bits_width", row.low_bits_width.to_string()),
                ("upper_bits", row.upper_bits.to_string()),
                ("buckets", row.buckets.to_string()),
                ("dense_bytes", row.dense_bytes.to_string()),
                ("ef_bytes", row.ef_bytes.to_string()),
                ("ef_bytes_no_select", row.ef_bytes_no_select.to_string()),
                ("upper_bytes", row.upper_bytes.to_string()),
                ("low_bytes", row.low_bytes.to_string()),
                ("select_hint_bytes", row.select_hint_bytes.to_string()),
                ("value_bytes", row.value_bytes.to_string()),
                ("bits_per_crossing", format!("{:.4}", row.bits_per_crossing)),
                (
                    "bits_per_crossing_no_select",
                    format!("{:.4}", row.bits_per_crossing_no_select),
                ),
                (
                    "bits_per_crossing_with_values",
                    format!("{:.4}", row.bits_per_crossing_with_values),
                ),
                (
                    "dense_bits_per_crossing",
                    format!("{:.4}", row.dense_bits_per_crossing),
                ),
                (
                    "space_ratio_dense_over_ef",
                    format!("{:.4}", row.space_ratio_dense_over_ef),
                ),
                (
                    "ns_per_access_dense",
                    format!("{:.4}", row.ns_per_access_dense),
                ),
                ("ns_per_access_ef", format!("{:.4}", row.ns_per_access_ef)),
                ("access_ratio", format!("{:.4}", row.access_ratio)),
                (
                    "instructions_per_access_dense",
                    format!("{:.4}", row.instructions_per_access_dense),
                ),
                (
                    "instructions_per_access_ef",
                    format!("{:.4}", row.instructions_per_access_ef),
                ),
                (
                    "access_ratio_instructions",
                    format!("{:.4}", row.access_ratio_instructions),
                ),
                (
                    "ns_per_access_dense_random",
                    format!("{:.4}", row.ns_per_access_dense_random),
                ),
                (
                    "ns_per_access_ef_random",
                    format!("{:.4}", row.ns_per_access_ef_random),
                ),
                (
                    "access_ratio_random",
                    format!("{:.4}", row.access_ratio_random),
                ),
                (
                    "instructions_per_access_dense_random",
                    format!("{:.4}", row.instructions_per_access_dense_random),
                ),
                (
                    "instructions_per_access_ef_random",
                    format!("{:.4}", row.instructions_per_access_ef_random),
                ),
                (
                    "access_ratio_instructions_random",
                    format!("{:.4}", row.access_ratio_instructions_random),
                ),
                (
                    "ns_per_access_ef_positional",
                    format!("{:.4}", row.ns_per_access_ef_positional),
                ),
                (
                    "instructions_per_access_ef_positional",
                    format!("{:.4}", row.instructions_per_access_ef_positional),
                ),
                (
                    "count_ones_per_access_ef",
                    format!("{:.4}", row.count_ones_per_access_ef),
                ),
                (
                    "count_ones_per_access_ef_positional",
                    format!("{:.4}", row.count_ones_per_access_ef_positional),
                ),
                (
                    "target_feature_popcnt",
                    row.target_feature_popcnt.to_string(),
                ),
                (
                    "instructions_per_access_ef_popcnt_credited",
                    format!("{:.4}", row.instructions_per_access_ef_popcnt_credited),
                ),
                (
                    "access_ratio_instructions_popcnt_credited",
                    format!("{:.4}", row.access_ratio_instructions_popcnt_credited),
                ),
                (
                    "c2_contingent_on_popcnt",
                    row.c2_contingent_on_popcnt.to_string(),
                ),
                ("values_equal", row.values_equal.to_string()),
                (
                    "positional_values_equal",
                    row.positional_values_equal.to_string(),
                ),
                ("absent_probed", row.absent_probed.to_string()),
                (
                    "absent_sentinel_correct",
                    row.absent_sentinel_correct.to_string(),
                ),
                (
                    "mutant_value_mismatches",
                    row.mutant_value_mismatches.to_string(),
                ),
                (
                    "mutant_low_bits_width",
                    row.mutant_low_bits_width.to_string(),
                ),
                ("shipped_vertices", row.shipped_vertices.to_string()),
                ("mirror_vertices", row.mirror_vertices.to_string()),
                ("centroid_vertices", row.centroid_vertices.to_string()),
                ("mirror_on_edge", row.mirror_on_edge.to_string()),
                ("ghz", format!("{:.4}", row.ghz)),
                ("inner_probe", row.inner_probe.to_string()),
                ("inner_random", row.inner_random.to_string()),
                ("c1_holds", row.c1_holds.to_string()),
                ("c1_holds_no_select", row.c1_holds_no_select.to_string()),
                ("c1_holds_with_values", row.c1_holds_with_values.to_string()),
                ("c2_holds", row.c2_holds.to_string()),
                ("c2_holds_ns", row.c2_holds_ns.to_string()),
                ("c3_holds", row.c3_holds.to_string()),
            ]);
        }
    }
}

fn main() {
    if !std::env::args().any(|arg| arg == "--bench") {
        return;
    }

    let prereg = isomesh::experiment!("P-109");

    #[cfg(target_os = "linux")]
    common::experiment::run(prereg, experiment::run);

    // `experiment_p12`'s precedent. C2's verdict is an instruction ratio and
    // instructions come from `perf_event_open`; a nanosecond on a governed CPU
    // cannot carry it (`M-281`), and a recorded zero would be a fabricated cost.
    #[cfg(not(target_os = "linux"))]
    {
        eprintln!(
            "{} scores its access clause on hardware instruction counts, and this platform \
             has no `perf_event_open`. There is no clock substitute.",
            prereg.id
        );
        std::process::exit(1);
    }
}
