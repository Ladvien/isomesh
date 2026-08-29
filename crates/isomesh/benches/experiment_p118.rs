//! **P-118 — Neal's superaccumulator for cross-cell float accumulation, aimed at `M-177`.**
//!
//! Ticket: R-118. Pre-registered before this harness existed.
//!
//! ```bash
//! cargo bench --bench experiment_p118
//! ```
//!
//! Writes `docs/experiments/p-118.csv`. **Linux only.** C2 is a cost bound and
//! the unit it is scored in is retired instructions (`✗24`, `M-280`), the
//! counters are `perf_event_open`, and off Linux this bench `exit(1)`s rather
//! than record a zero for a column it could not measure.
//!
//! # What was missing
//!
//! **Nothing in the shipped crate accumulates in fixed point, and the
//! distinction the registration turns on is not "ordered" versus "unordered" but
//! "rounded twelve times" versus "rounded once".** `A-016` already made every
//! reduction on the dual vertex path a deterministic function of the *multiset*
//! of its terms: `equivariant::sum` sorts by `(|x|, x)` through `total_cmp`,
//! which is a total order, so `sum_equivariant` cannot be perturbed by the order
//! the terms arrived in. What it still is, is **inexact** — twelve slots summed
//! smallest-magnitude-first is eleven roundings, and eleven roundings are not a
//! function that commutes with negating or permuting the terms.
//!
//! Neal's superaccumulator is a different object: 67 chunks of 64 bits, each
//! holding a 32-bit-wide slice of the double's exponent range with 32 bits of
//! carry headroom, so every addition is an **exact integer** addition and the
//! whole accumulation is exactly the mathematical sum. It is rounded **once**,
//! at read-out. That makes it exactly order-independent for a much stronger
//! reason than a canonical sort: the accumulator's state is the exact sum, and
//! integer addition is associative and commutative with no rounding to lose.
//!
//! That is the second route to `P-101`'s question the registration names. `P-101`
//! needed an *invariant key*; this needs none, because there is nothing to
//! order.
//!
//! # What the committed data already fixes, before this run
//!
//! `docs/experiments/p-101.csv`, `320` rows, is on disk and answers most of C1
//! arithmetically. Counting its `equivariance` block, `32` dual rows per arm:
//!
//! | arm | rows at `elements_vertex_exact = 48` | rows at `pure_permutation_exact = 6` | rows at `pure_sign_flip_exact = 8` |
//! |---|---|---|---|
//! | shipped baseline | **6** of 32 | 8 of 32 | 6 of 32 |
//! | `edge_slot` (transcription) | 6 | 8 | 6 |
//! | `abs_value_abs_offset` (P-101's registered arm) | 6 | 8 | 6 |
//! | `abs_value_abs_offset_equivariant_normal` | 6 | 12 | 6 |
//! | **`edge_slot_equivariant_normal`** | 6 | **32** | 6 |
//!
//! Two things follow, and both were readable before this file existed.
//!
//! 1. **The permutation obstruction is entirely `vec3::length`'s naive dot.**
//!    Changing *only* the normal's normalisation to a magnitude-ordered dot takes
//!    `pure_permutation_exact` from 8 of 32 rows to **32 of 32** — every row
//!    fully equivariant under all six pure axis relabellings. Nothing about the
//!    crossing accumulation moved that column at all.
//! 2. **The sign-flip half does not move for any arm.** `pure_sign_flip_exact`
//!    reads 6 of 32 for the shipped extractor and 6 of 32 for all four of
//!    `P-101`'s arms, including the one that fixed permutations completely. That
//!    is `M-177`, and `M-177` calls it structural.
//!
//! `elements_vertex_exact = 48` requires **both** halves. So **C1 is
//! arithmetically out of reach before this run** unless exactness buys something
//! a magnitude-ordered sum could not — which is a real possibility and the
//! reason the row is run rather than closed on the table above. It is run, the
//! number is produced, and the arms are arranged so that a failure is *located*
//! rather than argued: `superaccumulator_all_twelve_and_normal` moves the one
//! reduction item 1 names as well, so if the sign-flip column still reads 6 of
//! 32 with **every** reduction on the path exact, `M-177` is not an accumulation
//! defect and `P-101`, `P-118` and `R-118` close together.
//!
//! # SHARE
//!
//! The dual vertex solve only. `M-25` puts the sharp-feature solve at 3% over
//! Surface Nets and at **6.5% on Zen 3**, which is the machine this phase runs
//! on, so the ceiling is `1/(1 − 0.065) = 1.070×`. **C1 is a correctness clause
//! and C2 a cost bound; neither is or may become a speedup claim**, and no
//! column in this file is a speedup.
//!
//! # The five arms, one build, one run (`M-281`)
//!
//! | `rule` | centroid's 3 | `solve_with`'s 9 | normal length | replica? |
//! |---|---|---|---|---|
//! | `edge_slot` | `sum_equivariant` | `sum_equivariant` | `vec3::length` | **yes** |
//! | `ordered_naive` | running `f64` sum, ascending edge label | same | `vec3::length` | no |
//! | `superaccumulator_solve_nine` | `sum_equivariant` | **superaccumulator** | `vec3::length` | no |
//! | `superaccumulator_all_twelve` | **superaccumulator** | **superaccumulator** | `vec3::length` | no |
//! | `superaccumulator_all_twelve_and_normal` | **superaccumulator** | **superaccumulator** | **superaccumulator** | no |
//!
//! `edge_slot` is the **instrument check**, and it is what makes a second copy
//! of the crate's arithmetic admissible here at all (`P-61`'s rule, restated at
//! `experiment_p117.rs:53-56`): it is this file's transcription of
//! `HermiteCell::from_corners`, `solve::solve_with` and `apply_clamp`, and it is
//! asserted **bit-identical** to the shipped extractor over every position,
//! normal and index on all 32 equivariance rows and all 48 golden rows. Without
//! it, a difference in a later arm could be a transcription error rather than the
//! change under test.
//!
//! `ordered_naive` is **C2's denominator**, and it is Neal's own comparand:
//! "simple ordered summation" is a running `f64` sum, not the crate's sorted one.
//! The crate's own cost is recorded beside it as
//! `instructions_per_solve_edge_slot`, so a reader can see what an extraction
//! actually pays today without having to infer it.
//!
//! `superaccumulator_solve_nine` exists because the registration's METHOD names
//! the comparand as `solve_with`'s **nine** accumulators, while C1 is about the
//! dual *vertex* and `HermiteCell::centroid`'s three are on the same path. Both
//! readings are measured rather than chosen, and
//! `superaccumulator_all_twelve` — all twelve — is the registered arm, because an
//! accumulator that leaves one inexact reduction in place is not "an exactly
//! order-independent accumulator".
//!
//! Two reductions are **not** varied by any arm, because neither is one of the
//! registration's accumulators and both are already exact-order-independent
//! forms from `M-24`: `dot_equivariant(normal, position − centroid)`, the
//! three-term dot that forms `d`, and `Symmetric3`'s `determinant`, `adjugate`
//! and `mul_vec`.
//!
//! # Neal's superaccumulator, as implemented here
//!
//! Radford Neal, *Fast Exact Summation Using Small and Large Superaccumulators*
//! (2015), the **small** accumulator. The layout is his, and the constants are
//! his arithmetic rather than round numbers:
//!
//! - the exponent's low `5` bits index within a chunk, the high `6` index the
//!   chunk, so a chunk covers `2^5 = 32` bits of value;
//! - `(1 << 6) + 3 = 67` chunks of 64 bits, each carrying a 32-bit slice, which
//!   is where the **32 bits of overlap** the registration names come from — the
//!   headroom a chunk has above the slice it holds;
//! - one `f64` is split across exactly two adjacent chunks, `mantissa << low_exp`
//!   masked to 32 bits into the low chunk and `mantissa >> (32 − low_exp)`,
//!   arithmetic, into the high one. Those two pieces reconstruct
//!   `mantissa · 2^low_exp` exactly, because `(x mod 2^32) + 2^32 · ⌊x / 2^32⌋ = x`.
//!
//! The scale identity, which is what makes the accumulator exact and is asserted
//! by the self-test rather than asserted in prose: with `exp` the raw biased
//! exponent and `mantissa` the signed 53-bit significand,
//! `Σᵢ chunk[i] · 2^(32i) = Σ_terms mantissa · 2^exp = (Σ_terms value) · 2^1075`.
//!
//! **No carry propagation, and the bound is a column rather than a claim.** Neal
//! propagates carries every 63 additions because his accumulator sums arrays.
//! This one sums a dual cell, which has at most `EDGE_COUNT = 12` crossings, so
//! a chunk gains at most `2^52 + 2^32 < 2^53` per addition and at most
//! `12 · 2^53 < 2^57` in total — 64 bits of headroom against 57 bits of use.
//! `worst_chunk_magnitude_bits` is measured over the whole cost corpus and
//! recorded, so the margin is read off the artefact.
//!
//! Read-out is round-to-nearest, ties-to-even, done in integers: carry-propagate
//! the 67 chunks into two's-complement base-`2^32` limbs, take the magnitude,
//! and round `m · 2^-1075` onto the `f64` grid. It is a deterministic function of
//! the exact sum and of nothing else, which is the whole property C1 is
//! denominated in.
//!
//! # C2's unit, stated once
//!
//! `cost_ratio` is the **instruction** ratio, `ordered_naive` in the
//! denominator. On a governed CPU a nanosecond is not a unit (`M-280`; this
//! machine spans 1.96–5.62 GHz under `powersave`/`balance_performance`) and a
//! wall-clock ratio is never a gate (`✗24`). The registered
//! `ns_per_solve_ordered` and `ns_per_solve_superacc` columns are recorded, and
//! `ns_cost_ratio` beside them, and **none of the three is what C2 is scored
//! on**. Both cost arms solve the *same* pre-built corpus of `CellCrossings`, so
//! the crossing construction and the normal normalisation are outside the
//! counted window and the only difference in it is the reduction.
//!
//! # The vacuity control
//!
//! The registration's own words: a harness that cannot reproduce the baseline
//! cannot measure a move away from it. So:
//!
//! - every row carries `baseline_elements_vertex_exact`, the **shipped**
//!   extractor's reading measured in this same run on that same configuration;
//! - on the 32 equivariance rows it is asserted equal to `p-101.csv`'s
//!   `elements_vertex_exact_baseline`, which is the number `✗79 / M-413` quotes,
//!   and the aggregate `baseline_rows_at_48` is asserted to be `6`;
//! - the `edge_slot` arm is asserted bit-identical to the shipped extractor;
//! - `golden_hashes.json` is asserted equal to the shipped extractor's own hash
//!   on all 48 golden dual rows before `hashes_moved` is believed;
//! - and the accumulator itself is tested before any mesh is built:
//!   `superacc_exact_cases`, `superacc_integer_oracle_cases`,
//!   `superacc_permutation_cases` and `superacc_negation_cases` are all asserted
//!   non-zero and all-passing. An accumulator *claimed* exactly
//!   order-independent and not shown to be is `M-44`'s zero that could not have
//!   been non-zero.
//!
//! # Why the golden block carries an equivariance reading too
//!
//! Every registered column has to be populated on every row, and a golden row
//! has no equivariance reading unless one is measured there — so one is. The
//! golden configuration is `origin = lo`, `cell_size = (hi − lo)/(samples − 1)`,
//! at 17, 25 and 33 samples. `2L/16` and `2L/32` are dyadic and their grids
//! mirror bit-exactly; **`2L/24` does not**, which is exactly why `P-57`'s 25³
//! fixture uses `3L/32` instead. The `grid_mirrors` column carries that per row,
//! `false` on the eight golden 25³ configurations, and **C1 is scored on the
//! `equivariance` block only** — a row whose grid does not mirror would be
//! falsified by the fixture rather than by the extractor.

// Exact comparisons on purpose: every clause here is stated in bits.
#![allow(clippy::float_cmp)]

mod common;

use std::cmp::Ordering;
use std::collections::HashMap;
use std::hint::black_box;
use std::path::PathBuf;
use std::time::Instant;

use isomesh::dual::{CellVertices, VertexRule};
use isomesh::dual_contouring::solve::{LAMBDA, dot_equivariant};
use isomesh::dual_contouring::{CLAMP_EPSILON, DualContouring};
use isomesh::fields::ReferenceField;
use isomesh::manifold_dual_contouring::ManifoldDualContouring;
use isomesh::marching_cubes::ambiguity::joined_mask;
use isomesh::marching_cubes::table::{
    EDGE_CORNERS, EDGE_COUNT, NO_EDGE, edge_offset, is_inside, segment_links,
};
use isomesh::validate::{AccuracyConfig, accuracy, mesh_hash};
use isomesh::{MeshBuffer, RuntimeShape3, Sdf};

use crate::common::counters::{MIN_TIME_RATIO, Probe};

// ─── clause constants ───────────────────────────────────────────────────────

/// C2's bar: "cost under 2x simple ordered summation, which is Neal's own
/// figure".
const C2_BAR: f64 = 2.0;

/// `elements_vertex_exact = 48` is C1's target and the octahedral group's order.
const GROUP_ORDER: usize = 48;

/// The shipped extractor's rows at 48, from `p-101.csv`. Asserted, not assumed.
const BASELINE_ROWS_AT_48: usize = 6;

/// Median of this many counted windows per cost arm.
const REPS: usize = 9;

/// About this long per counter window, so the `perf_event` round trip is
/// negligible against what it measures.
const TARGET_BATCH_NS: f64 = 10_000_000.0;

/// Ceiling on the batch factor, so a corpus that solved faster than expected
/// cannot turn one window into a minute.
const MAX_INNER: usize = 1 << 16;

/// The field the cost corpus is built from.
const COST_FIELD: &str = "sphere";

// ─── the crate's own reductions, transcribed ────────────────────────────────
//
// `crate::equivariant` is private, so these are copies rather than calls, and
// each carries the source line it came from. They are copies of eleven lines
// whose behaviour the `edge_slot` arm asserts bit-for-bit against the shipped
// extractor, which is the only thing that makes a second copy of an instrument
// admissible here (`experiment_p117.rs:53-56`).

/// Whether `a` sums before `b`: smaller magnitude first, then smaller value.
/// `equivariant.rs`'s `precedes`.
#[inline]
fn precedes(a: f64, b: f64) -> bool {
    match a.abs().total_cmp(&b.abs()) {
        Ordering::Less => true,
        Ordering::Greater => false,
        Ordering::Equal => a.total_cmp(&b) == Ordering::Less,
    }
}

/// Insertion sort, ascending by magnitude then by signed value.
#[inline]
fn sort_by_magnitude<const N: usize>(t: &mut [f64; N]) {
    let mut i = 1;
    while i < N {
        let mut j = i;
        while j > 0 && precedes(t[j], t[j - 1]) {
            t.swap(j - 1, j);
            j -= 1;
        }
        i += 1;
    }
}

/// Sum smallest-magnitude-first. `equivariant::sum`, which `A-016` landed.
#[inline]
fn sum_equivariant<const N: usize>(mut t: [f64; N]) -> f64 {
    sort_by_magnitude(&mut t);
    let mut acc = 0.0;
    for v in t {
        acc += v;
    }
    acc
}

/// Multiply smallest-magnitude-first. `equivariant::product`.
#[inline]
fn mul_equivariant(mut t: [f64; 3]) -> f64 {
    sort_by_magnitude(&mut t);
    (t[0] * t[1]) * t[2]
}

/// `cube::corner_offset`, which is `pub(crate)`.
#[inline]
fn corner_offset(corner: u8) -> [u32; 3] {
    [
        u32::from(corner & 1),
        u32::from((corner >> 1) & 1),
        u32::from((corner >> 2) & 1),
    ]
}

/// `cube::place`, which is `pub(crate)`.
#[inline]
fn place(lo: f64, hi: f64, d: f64) -> f64 {
    (lo + hi) * 0.5 + (hi - lo) * d
}

/// `vec3::dot`, the naive axis-order sum. **This is `P-101`'s located
/// obstruction**, and the only reduction the fifth arm changes.
#[inline]
fn naive_dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

#[inline]
fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

#[inline]
fn scale(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}

// ─── Neal's small superaccumulator ──────────────────────────────────────────

/// The exponent's low bits, which index within a chunk. Neal's
/// `XSUM_LOW_EXP_BITS`.
const LOW_EXP_BITS: u32 = 5;

/// Mask for those bits.
const LOW_EXP_MASK: usize = (1usize << LOW_EXP_BITS) - 1;

/// The exponent's high bits, which index the chunk. `11 − 5`.
const HIGH_EXP_BITS: u32 = 11 - LOW_EXP_BITS;

/// Neal's `XSUM_SCHUNKS = (1 << XSUM_HIGH_EXP_BITS) + 3`.
const SCHUNKS: usize = (1usize << HIGH_EXP_BITS) + 3;

/// Bits of value one chunk carries. Neal's `XSUM_LOW_MANTISSA_BITS = 1 << 5`.
const LOW_MANTISSA_BITS: u32 = 1 << LOW_EXP_BITS;

/// Mask for one chunk's slice, and for one base-`2^32` limb.
const LIMB_MASK: u64 = (1u64 << LOW_MANTISSA_BITS) - 1;

/// The `f64` significand's stored bits.
const MANTISSA_BITS: u32 = 52;

/// The `f64` exponent field's mask.
const EXP_MASK: u64 = 0x7ff;

/// The implicit bit.
const IMPLICIT: u64 = 1u64 << MANTISSA_BITS;

/// `Σ chunk[i] · 2^(32i)` is the exact sum scaled by `2^SCALE`.
const SCALE: u32 = 1075;

/// Base-`2^32` limbs the read-out canonicalises into: `SCHUNKS` for the chunks
/// themselves plus five for the sign extension and the final carry.
const LIMBS: usize = SCHUNKS + 5;

/// Neal's small superaccumulator: 67 chunks of 64 bits, each holding a 32-bit
/// slice of the value with 32 bits of carry headroom.
///
/// Every addition is exact and integral, so the state **is** the exact sum of
/// the terms and is therefore a function of their multiset alone. Rounding
/// happens once, in [`Superacc::round`].
#[derive(Clone)]
struct Superacc {
    chunk: [i64; SCHUNKS],
}

impl Superacc {
    const fn new() -> Self {
        Self {
            chunk: [0; SCHUNKS],
        }
    }

    /// One exact addition.
    ///
    /// No carry propagation: see the module doc's `12 · 2^53 < 2^57` bound, and
    /// `worst_chunk_magnitude_bits` for the margin actually used.
    ///
    /// # Panics
    ///
    /// On an infinity or a NaN. Neal's accumulator has a path for those; this
    /// one does not, because a dual cell that produced one has a field defect and
    /// recording a number for it would hide that.
    #[inline]
    fn add(&mut self, value: f64) {
        let bits = value.to_bits();
        let raw_exp = ((bits >> MANTISSA_BITS) & EXP_MASK) as usize;
        let frac = bits & (IMPLICIT - 1);
        assert!(
            raw_exp != EXP_MASK as usize,
            "the superaccumulator was handed a non-finite term: {value:?}"
        );
        let (exp, magnitude) = if raw_exp == 0 {
            if frac == 0 {
                // Both encodings of zero contribute nothing, exactly.
                return;
            }
            // Subnormal: exponent 1, no implicit bit.
            (1usize, frac)
        } else {
            (raw_exp, frac | IMPLICIT)
        };

        let mut mantissa = magnitude as i64;
        if bits >> 63 != 0 {
            mantissa = -mantissa;
        }

        let low_exp = (exp & LOW_EXP_MASK) as u32;
        let high_exp = exp >> LOW_EXP_BITS;

        // `(x mod 2^32) + 2^32 · floor(x / 2^32) = x`, with `x = mantissa · 2^low_exp`.
        let low = ((mantissa as u64) << low_exp) & LIMB_MASK;
        let high = mantissa >> (LOW_MANTISSA_BITS - low_exp);
        self.chunk[high_exp] += low as i64;
        self.chunk[high_exp + 1] += high;
    }

    /// The largest `|chunk|`, in bits. The headroom census.
    fn worst_chunk_bits(&self) -> u32 {
        let worst = self
            .chunk
            .iter()
            .map(|c| c.unsigned_abs())
            .max()
            .unwrap_or(0);
        u64::BITS - worst.leading_zeros()
    }

    /// The exact sum, rounded to the nearest `f64`, ties to even.
    ///
    /// A deterministic function of the exact sum and of nothing else.
    ///
    /// # Panics
    ///
    /// If the exact sum is out of `f64` range, or if the accumulator overflowed
    /// its limbs. Both are refusals rather than saturations: a saturated sum is a
    /// number that names no accumulation.
    fn round(&self) -> f64 {
        // ── canonicalise: two's complement, base 2^32 ──
        let mut limb = [0u32; LIMBS];
        let mut carry: i64 = 0;
        for (slot, &c) in limb.iter_mut().zip(self.chunk.iter()) {
            let v = c
                .checked_add(carry)
                .expect("chunk + carry overflows i64: the accumulator overflowed");
            *slot = (v as u64 & LIMB_MASK) as u32;
            carry = v >> LOW_MANTISSA_BITS;
        }
        let mut c = carry;
        for slot in limb[SCHUNKS..].iter_mut() {
            *slot = (c as u64 & LIMB_MASK) as u32;
            c >>= LOW_MANTISSA_BITS;
        }
        assert!(
            limb[LIMBS - 1] == 0 || limb[LIMBS - 1] == u32::MAX,
            "the top limb is not pure sign extension: the accumulator overflowed"
        );

        let negative = limb[LIMBS - 1] & 0x8000_0000 != 0;
        if negative {
            let mut borrow = 1u64;
            for slot in limb.iter_mut() {
                let v = u64::from(!*slot) + borrow;
                *slot = (v & LIMB_MASK) as u32;
                borrow = v >> LOW_MANTISSA_BITS;
            }
        }

        // ── round the magnitude `m · 2^-1075` onto the f64 grid ──
        //
        // Representable values are `q · 2^(g − 1075)` with `g >= 1` and
        // `q < 2^53`, and additionally `q >= 2^52` whenever `g >= 2`. So the
        // grid step for an `m` of `p` bits is `2^g` with `g = max(p − 53, 1)`.
        let p = bit_length(&limb);
        if p == 0 {
            return if negative { -0.0 } else { 0.0 };
        }
        let mut g = (i64::from(p) - 53).max(1) as u32;
        let mut q = shift_right(&limb, g);
        match dropped_versus_half(&limb, g) {
            Ordering::Greater => q += 1,
            Ordering::Equal => q += q & 1,
            Ordering::Less => {}
        }
        if q == 1u64 << 53 {
            q = 1u64 << MANTISSA_BITS;
            g += 1;
        }
        assert!(
            u64::from(g) < EXP_MASK,
            "the exact sum is out of f64 range: exponent field would be {g}"
        );

        let bits = if q < IMPLICIT {
            // Subnormal, which can only happen at `g == 1`.
            debug_assert_eq!(g, 1, "a subnormal read-out at g = {g}");
            q
        } else {
            (u64::from(g) << MANTISSA_BITS) | (q - IMPLICIT)
        };
        f64::from_bits(if negative { bits | (1u64 << 63) } else { bits })
    }
}

/// One past the most significant set bit of a base-`2^32` magnitude, or `0`.
fn bit_length(limb: &[u32; LIMBS]) -> u32 {
    for (i, &l) in limb.iter().enumerate().rev() {
        if l != 0 {
            return i as u32 * LOW_MANTISSA_BITS + (u32::BITS - l.leading_zeros());
        }
    }
    0
}

/// `m >> g`, which fits a `u64` by construction: `g` is chosen so the quotient
/// is 53 or 54 bits.
fn shift_right(limb: &[u32; LIMBS], g: u32) -> u64 {
    let word = (g / LOW_MANTISSA_BITS) as usize;
    let bit = g % LOW_MANTISSA_BITS;
    let mut acc: u128 = 0;
    for k in (0..3).rev() {
        let l = limb.get(word + k).copied().unwrap_or(0);
        acc = (acc << LOW_MANTISSA_BITS) | u128::from(l);
    }
    let q = (acc >> bit) as u64;
    assert!(
        q < 1u64 << 54,
        "the read-out quotient is {q} bits, so `g` was chosen wrong"
    );
    q
}

/// How the dropped low `g` bits of `m` compare with `2^(g − 1)`.
fn dropped_versus_half(limb: &[u32; LIMBS], g: u32) -> Ordering {
    let half = g - 1;
    let word = (half / LOW_MANTISSA_BITS) as usize;
    let bit = half % LOW_MANTISSA_BITS;
    if limb[word] >> bit & 1 == 0 {
        return Ordering::Less;
    }
    let lower_in_word = limb[word] & ((1u32 << bit) - 1) != 0;
    if lower_in_word || limb[..word].iter().any(|l| *l != 0) {
        Ordering::Greater
    } else {
        Ordering::Equal
    }
}

/// Sum a slice exactly, with the headroom census compiled in or out.
///
/// One implementation, two monomorphisations: the counted cost path takes
/// `CENSUS = false` and pays nothing for the census, and the correctness path
/// takes `CENSUS = true`. A second copy of the loop would be a second thing to
/// drift.
#[inline]
fn superacc_sum<const CENSUS: bool>(terms: &[f64], worst_bits: &mut u32) -> f64 {
    let mut acc = Superacc::new();
    for &t in terms {
        acc.add(t);
    }
    if CENSUS {
        *worst_bits = (*worst_bits).max(acc.worst_chunk_bits());
    }
    acc.round()
}

// ─── the symmetric solve, transcribed unchanged ─────────────────────────────
//
// `determinant`, `adjugate` and `mul_vec` are `M-24`'s forms and are **not**
// what this experiment varies. They are here because `solve_with` is one
// function and only its crossing accumulation is under test.

/// A symmetric 3×3 matrix, stored as its six distinct entries.
#[derive(Clone, Copy)]
struct Symmetric3 {
    xx: f64,
    xy: f64,
    xz: f64,
    yy: f64,
    yz: f64,
    zz: f64,
}

impl Symmetric3 {
    #[inline]
    fn outer(n: [f64; 3]) -> [f64; 6] {
        [
            n[0] * n[0],
            n[0] * n[1],
            n[0] * n[2],
            n[1] * n[1],
            n[1] * n[2],
            n[2] * n[2],
        ]
    }

    #[inline]
    fn from_entries(e: [f64; 6]) -> Self {
        Self {
            xx: e[0],
            xy: e[1],
            xz: e[2],
            yy: e[3],
            yz: e[4],
            zz: e[5],
        }
    }

    #[inline]
    fn regularized(mut self, lambda: f64) -> Self {
        self.xx += lambda;
        self.yy += lambda;
        self.zz += lambda;
        self
    }

    #[inline]
    fn adjugate(self) -> Self {
        Self {
            xx: self.yy * self.zz - self.yz * self.yz,
            xy: self.xz * self.yz - self.xy * self.zz,
            xz: self.xy * self.yz - self.xz * self.yy,
            yy: self.xx * self.zz - self.xz * self.xz,
            yz: self.xy * self.xz - self.xx * self.yz,
            zz: self.xx * self.yy - self.xy * self.xy,
        }
    }

    #[inline]
    fn determinant(self) -> f64 {
        sum_equivariant([
            mul_equivariant([self.xx, self.yy, self.zz]),
            2.0 * mul_equivariant([self.xy, self.yz, self.xz]),
            -mul_equivariant([self.xx, self.yz, self.yz]),
            -mul_equivariant([self.yy, self.xz, self.xz]),
            -mul_equivariant([self.zz, self.xy, self.xy]),
        ])
    }

    #[inline]
    fn mul_vec(self, v: [f64; 3]) -> [f64; 3] {
        [
            dot_equivariant([self.xx, self.xy, self.xz], v),
            dot_equivariant([self.xy, self.yy, self.yz], v),
            dot_equivariant([self.xz, self.yz, self.zz], v),
        ]
    }
}

// ─── one cell's crossings ───────────────────────────────────────────────────

/// One edge crossing.
#[derive(Clone, Copy)]
struct Crossing {
    /// World position — `cube::place` in the centred frame.
    position: [f64; 3],
    /// Unit surface normal, normalised by whichever rule this arm uses.
    normal: [f64; 3],
}

const EMPTY_CROSSING: Crossing = Crossing {
    position: [0.0; 3],
    normal: [0.0; 3],
};

/// The crossings on one cell's twelve edges, indexed by edge label.
#[derive(Clone)]
struct CellCrossings {
    edge: [Crossing; EDGE_COUNT],
    mask: u16,
}

/// How a normal's length is computed.
#[derive(Clone, Copy, PartialEq, Eq)]
enum Norm {
    /// `vec3::length`, i.e. `(g0² + g1² + g2²).sqrt()` left to right. Shipped.
    Naive,
    /// The superaccumulator over the three squares, then `sqrt`. Not shipped.
    Superacc,
}

/// How one family of accumulators reduces its terms.
#[derive(Clone, Copy, PartialEq, Eq)]
enum Reduce {
    /// One slot per edge label, reduced by `sum_equivariant`. Shipped (`A-016`).
    EdgeSlots,
    /// A running `f64` sum in ascending edge-label order. Neal's comparand and
    /// C2's denominator.
    Ordered,
    /// Neal's superaccumulator: exact fixed point, one rounding.
    Superacc,
}

/// One arm's configuration: which reduction each family uses.
#[derive(Clone, Copy, PartialEq, Eq)]
struct Config {
    /// `HermiteCell::centroid`'s three.
    centroid: Reduce,
    /// `solve::solve_with`'s nine.
    solve: Reduce,
    /// `vec3::length`'s dot.
    norm: Norm,
}

/// `HermiteCell::from_corners`, transcribed, with the normalisation switchable.
fn crossings_of<S: Sdf<Scalar = f64>>(
    sdf: &S,
    corner_values: &[f64; 8],
    cell_origin: [f64; 3],
    cell_size: f64,
    norm: Norm,
) -> CellCrossings {
    let mut out = CellCrossings {
        edge: [EMPTY_CROSSING; EDGE_COUNT],
        mask: 0,
    };
    for (edge, [lo, hi]) in EDGE_CORNERS.iter().copied().enumerate() {
        let (a, b) = (corner_values[lo as usize], corner_values[hi as usize]);
        if is_inside(a) == is_inside(b) {
            continue;
        }

        let d = edge_offset(a, b);
        let (lo_offset, hi_offset) = (corner_offset(lo), corner_offset(hi));
        let mut position = [0.0f64; 3];
        for (axis, slot) in position.iter_mut().enumerate() {
            let from = f64::from(lo_offset[axis]);
            let to = f64::from(hi_offset[axis]);
            *slot = cell_origin[axis] + cell_size * place(from, to, d);
        }

        let gradient = sdf.gradient(position);
        let length = match norm {
            Norm::Naive => naive_dot(gradient, gradient).sqrt(),
            Norm::Superacc => {
                let mut ignored = 0u32;
                superacc_sum::<false>(
                    &[
                        gradient[0] * gradient[0],
                        gradient[1] * gradient[1],
                        gradient[2] * gradient[2],
                    ],
                    &mut ignored,
                )
                .sqrt()
            }
        };
        let normal = scale(gradient, length.recip());

        out.edge[edge] = Crossing { position, normal };
        out.mask |= 1 << edge;
    }
    out
}

/// Reduce one family's twelve edge slots.
///
/// All three arms are handed the same twelve-slot array and the same mask, so
/// the only thing that varies between them is the reduction. The inactive slots
/// hold `0.0`: `sum_equivariant` sorts them to the front where they contribute
/// nothing, and the other two skip them by mask.
#[inline]
fn reduce_slots(t: &[f64; EDGE_COUNT], mask: u16, how: Reduce, worst_bits: &mut u32) -> f64 {
    match how {
        Reduce::EdgeSlots => sum_equivariant(*t),
        Reduce::Ordered => {
            let mut acc = 0.0f64;
            for (edge, value) in t.iter().enumerate() {
                if mask & (1 << edge) != 0 {
                    acc += *value;
                }
            }
            acc
        }
        Reduce::Superacc => {
            let mut acc = Superacc::new();
            for (edge, value) in t.iter().enumerate() {
                if mask & (1 << edge) != 0 {
                    acc.add(*value);
                }
            }
            *worst_bits = (*worst_bits).max(acc.worst_chunk_bits());
            acc.round()
        }
    }
}

/// The same reduction with the headroom census compiled out, for the counted
/// window.
#[inline]
fn reduce_slots_fast(t: &[f64; EDGE_COUNT], mask: u16, how: Reduce) -> f64 {
    match how {
        Reduce::EdgeSlots => sum_equivariant(*t),
        Reduce::Ordered => {
            let mut acc = 0.0f64;
            for (edge, value) in t.iter().enumerate() {
                if mask & (1 << edge) != 0 {
                    acc += *value;
                }
            }
            acc
        }
        Reduce::Superacc => {
            let mut acc = Superacc::new();
            for (edge, value) in t.iter().enumerate() {
                if mask & (1 << edge) != 0 {
                    acc.add(*value);
                }
            }
            acc.round()
        }
    }
}

/// `solve::solve_with`, transcribed, with the two accumulator families
/// switchable.
///
/// `keep` restricts to a subset of the edges, which is what
/// `HermiteCell::restricted` does for the cycle rule. `census` selects the
/// monomorphisation: the counted cost window takes `false`.
fn solve_cell<const CENSUS: bool>(
    cell: &CellCrossings,
    keep: u16,
    lambda: f64,
    cfg: Config,
    worst_bits: &mut u32,
) -> Option<[f64; 3]> {
    let mask = cell.mask & keep;
    let count = mask.count_ones() as usize;
    if count == 0 {
        return None;
    }
    let reduce = |t: &[f64; EDGE_COUNT], how: Reduce, worst: &mut u32| -> f64 {
        if CENSUS {
            reduce_slots(t, mask, how, worst)
        } else {
            reduce_slots_fast(t, mask, how)
        }
    };

    // ── HermiteCell::centroid's three accumulators ──
    let inverse = (count as f64).recip();
    let mut axes = [[0.0f64; EDGE_COUNT]; 3];
    for edge in 0..EDGE_COUNT {
        if mask & (1 << edge) == 0 {
            continue;
        }
        for (slot, value) in axes.iter_mut().zip(cell.edge[edge].position) {
            slot[edge] = value;
        }
    }
    let mut centroid = [0.0f64; 3];
    for (slot, terms) in centroid.iter_mut().zip(axes.iter()) {
        *slot = reduce(terms, cfg.centroid, worst_bits) * inverse;
    }

    // ── solve_with's nine accumulators ──
    let mut m_terms = [[0.0f64; EDGE_COUNT]; 6];
    let mut g_terms = [[0.0f64; EDGE_COUNT]; 3];
    for edge in 0..EDGE_COUNT {
        if mask & (1 << edge) == 0 {
            continue;
        }
        let c = &cell.edge[edge];
        let normal = c.normal;
        let d = dot_equivariant(normal, sub(c.position, centroid));
        for (slot, value) in m_terms.iter_mut().zip(Symmetric3::outer(normal)) {
            slot[edge] = value;
        }
        for (slot, value) in g_terms
            .iter_mut()
            .zip([normal[0] * d, normal[1] * d, normal[2] * d])
        {
            slot[edge] = value;
        }
    }
    let mut m_sum = [0.0f64; 6];
    for (slot, terms) in m_sum.iter_mut().zip(m_terms.iter()) {
        *slot = reduce(terms, cfg.solve, worst_bits);
    }
    let mut g_sum = [0.0f64; 3];
    for (slot, terms) in g_sum.iter_mut().zip(g_terms.iter()) {
        *slot = reduce(terms, cfg.solve, worst_bits);
    }

    let a = Symmetric3::from_entries(m_sum).regularized(lambda);
    let adj = a.adjugate();
    let det = a.determinant();
    let offset = scale(adj.mul_vec(g_sum), det.recip());
    let x = [
        centroid[0] + offset[0],
        centroid[1] + offset[1],
        centroid[2] + offset[2],
    ];
    if x[0].is_finite() && x[1].is_finite() && x[2].is_finite() {
        Some(x)
    } else {
        None
    }
}

/// `dual_contouring::apply_clamp` with `Clamp::ToCell`, transcribed.
fn clamp_to_cell(x: [f64; 3], cell_origin: [f64; 3], cell_size: f64) -> [f64; 3] {
    let half = cell_size * 0.5;
    let inset = half * (1.0 - CLAMP_EPSILON);
    let mut out = x;
    for (axis, slot) in out.iter_mut().enumerate() {
        let centre = cell_origin[axis] + half;
        *slot = slot.clamp(centre - inset, centre + inset);
    }
    out
}

#[inline]
fn cell_origin_of(base: [u32; 3], origin: [f64; 3], cell_size: f64) -> [f64; 3] {
    [
        origin[0] + cell_size * f64::from(base[0]),
        origin[1] + cell_size * f64::from(base[1]),
        origin[2] + cell_size * f64::from(base[2]),
    ]
}

// ─── the two rules, as this bench's own `VertexRule`s ───────────────────────

/// `dual_contouring::Qef`, reimplemented so the accumulation can be swapped.
struct DcRule {
    cfg: Config,
}

impl VertexRule<f64> for DcRule {
    fn place<S: Sdf<Scalar = f64>>(
        &self,
        sdf: &S,
        corner: &[f64; 8],
        base: [u32; 3],
        origin: [f64; 3],
        cell_size: f64,
        out: &mut CellVertices<f64>,
    ) {
        let cell_origin = cell_origin_of(base, origin, cell_size);
        let cell = crossings_of(sdf, corner, cell_origin, cell_size, self.cfg.norm);
        let mut ignored = 0u32;
        let Some(x) = solve_cell::<false>(&cell, u16::MAX, LAMBDA, self.cfg, &mut ignored) else {
            return;
        };
        out.push_whole_cell(clamp_to_cell(x, cell_origin, cell_size));
    }
}

/// `manifold_dual_contouring::CycleQef`, reimplemented the same way.
///
/// `FaceAmbiguity::Separate` only, which is what `ManifoldDualContouring::new`
/// ships and what the committed golden hashes pin — so the `ambiguous` mask
/// handed to `joined_mask` is `0`, exactly as `CycleQef` computes it under that
/// setting.
struct CycleRule {
    cfg: Config,
}

impl VertexRule<f64> for CycleRule {
    fn place<S: Sdf<Scalar = f64>>(
        &self,
        sdf: &S,
        corner: &[f64; 8],
        base: [u32; 3],
        origin: [f64; 3],
        cell_size: f64,
        out: &mut CellVertices<f64>,
    ) {
        let cell_origin = cell_origin_of(base, origin, cell_size);

        let mut case = 0u8;
        for (c, value) in corner.iter().enumerate() {
            if is_inside(*value) {
                case |= 1 << c;
            }
        }
        let next = segment_links(case, joined_mask(corner, 0));

        let cell = crossings_of(sdf, corner, cell_origin, cell_size, self.cfg.norm);

        let mut visited = 0u16;
        let mut ignored = 0u32;
        for start in 0..EDGE_COUNT as u8 {
            if next[start as usize] == NO_EDGE || visited & (1 << start) != 0 {
                continue;
            }
            let mut edges = 0u16;
            let mut current = start;
            while visited & (1 << current) == 0 {
                visited |= 1 << current;
                edges |= 1 << current;
                current = next[current as usize];
            }
            let Some(x) = solve_cell::<false>(&cell, edges, LAMBDA, self.cfg, &mut ignored) else {
                continue;
            };
            out.push_component(clamp_to_cell(x, cell_origin, cell_size), edges);
        }
    }
}

// ─── one thing that can mesh a field ────────────────────────────────────────

/// Anything this harness can hand a field to.
///
/// A trait rather than a closure because the field type varies — the
/// equivariance arm wraps it in [`Rotated`] — and a closure cannot be generic.
trait Mesher {
    fn mesh<S: Sdf<Scalar = f64>>(
        &mut self,
        field: &S,
        shape: &RuntimeShape3,
        origin: [f64; 3],
        cell_size: f64,
        out: &mut MeshBuffer<f64>,
    );
}

impl<V: VertexRule<f64>> Mesher for DualContouring<f64, V> {
    fn mesh<S: Sdf<Scalar = f64>>(
        &mut self,
        field: &S,
        shape: &RuntimeShape3,
        origin: [f64; 3],
        cell_size: f64,
        out: &mut MeshBuffer<f64>,
    ) {
        out.reset();
        self.extract(field, shape, origin, cell_size, out)
            .expect("extraction");
    }
}

impl Mesher for ManifoldDualContouring<f64> {
    fn mesh<S: Sdf<Scalar = f64>>(
        &mut self,
        field: &S,
        shape: &RuntimeShape3,
        origin: [f64; 3],
        cell_size: f64,
        out: &mut MeshBuffer<f64>,
    ) {
        out.reset();
        self.extract(field, shape, origin, cell_size, out)
            .expect("extraction");
    }
}

// ─── the group, and the comparison keys (P-57's instrument, via P-101) ──────

/// The 6 axis permutations. Crossed with 8 sign patterns this is all 48.
const PERMS: [[usize; 3]; 6] = [
    [0, 1, 2],
    [0, 2, 1],
    [1, 0, 2],
    [1, 2, 0],
    [2, 0, 1],
    [2, 1, 0],
];

/// Probe points for the inverse round-trip check. Both zeros deliberately.
const PROBES: [[f64; 3]; 4] = [
    [0.3, -1.7, 2.9],
    [0.0, -0.0, 1.5],
    [-2.25, 0.656_25, -0.187_5],
    [1.0, 1.0, 1.0],
];

/// One element of the octahedral group, as a signed axis permutation.
#[derive(Clone, Copy, PartialEq, Eq)]
struct Element {
    perm: [usize; 3],
    sign: [i8; 3],
}

impl Element {
    #[inline]
    fn apply(self, p: [f64; 3]) -> [f64; 3] {
        std::array::from_fn(|k| {
            let v = p[self.perm[k]];
            if self.sign[k] < 0 { -v } else { v }
        })
    }

    fn inverse(self) -> Self {
        let mut perm = [0usize; 3];
        let mut sign = [0i8; 3];
        for k in 0..3 {
            let j = self.perm[k];
            perm[j] = k;
            sign[j] = self.sign[k];
        }
        Self { perm, sign }
    }

    fn det(self) -> i32 {
        let mut m = [[0i32; 3]; 3];
        for k in 0..3 {
            m[k][self.perm[k]] = i32::from(self.sign[k]);
        }
        m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
            - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
            + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
    }

    fn label(self) -> String {
        let mut s = String::from("perm=");
        for k in 0..3 {
            s.push(char::from(b'0' + self.perm[k] as u8));
        }
        s.push_str(";sign=");
        for k in 0..3 {
            s.push(if self.sign[k] < 0 { '-' } else { '+' });
        }
        s
    }

    fn short(self) -> String {
        let mut s = String::new();
        for k in 0..3 {
            s.push(char::from(b'0' + self.perm[k] as u8));
        }
        s.push('/');
        for k in 0..3 {
            s.push(if self.sign[k] < 0 { '-' } else { '+' });
        }
        s
    }
}

fn group() -> Vec<Element> {
    let mut out = Vec::with_capacity(GROUP_ORDER);
    for perm in PERMS {
        for bits in 0..8u8 {
            let sign = std::array::from_fn(|k| if bits & (1 << k) == 0 { 1i8 } else { -1i8 });
            out.push(Element { perm, sign });
        }
    }
    out
}

fn verify_group(g: &[Element]) {
    assert_eq!(g.len(), GROUP_ORDER, "the octahedral group has 48 elements");
    assert!(
        g[0] == Element {
            perm: [0, 1, 2],
            sign: [1, 1, 1]
        },
        "element 0 must be the identity: it is the negative control"
    );
    for (i, a) in g.iter().enumerate() {
        for b in &g[i + 1..] {
            assert!(a != b, "duplicate group element at {i}: {}", a.label());
        }
    }
    let (mut rotations, mut reflections) = (0, 0);
    for e in g {
        match e.det() {
            1 => rotations += 1,
            -1 => reflections += 1,
            d => panic!("{} has determinant {d}, not +/-1", e.label()),
        }
        let inv = e.inverse();
        for p in PROBES {
            let round = inv.apply(e.apply(p));
            for k in 0..3 {
                assert_eq!(
                    round[k].to_bits(),
                    p[k].to_bits(),
                    "{} does not round-trip {p:?} bit-exactly",
                    e.label()
                );
            }
        }
    }
    assert_eq!(rotations, 24, "24 elements must have det = +1");
    assert_eq!(reflections, 24, "24 elements must have det = -1");
}

/// `g·f`, i.e. `(g·f)(p) = f(g⁻¹·p)`.
struct Rotated<'a, S> {
    field: &'a S,
    g: Element,
    g_inv: Element,
}

impl<S: Sdf<Scalar = f64>> Sdf for Rotated<'_, S> {
    type Scalar = f64;

    #[inline]
    fn sample(&self, p: [f64; 3]) -> f64 {
        self.field.sample(self.g_inv.apply(p))
    }

    #[inline]
    fn gradient(&self, p: [f64; 3]) -> [f64; 3] {
        self.g.apply(self.field.gradient(self.g_inv.apply(p)))
    }
}

const NEGATIVE_ZERO: u64 = 1u64 << 63;

/// `−0.0` folded onto `+0.0`, identically on both sides of every comparison.
///
/// Both fixtures are centred with an odd sample count, so `0.0` is a grid
/// coordinate and `−(0.0)` is `−0.0`; left raw, every sign-flipping element
/// fails on a disagreement about which *encoding* of zero was written. `✗39`
/// settled this and `P-61` kept it.
#[inline]
fn comparison_key(v: f64) -> u64 {
    let b = v.to_bits();
    if b == NEGATIVE_ZERO { 0 } else { b }
}

fn vertex_keys(positions: &[[f64; 3]], g: Option<Element>) -> Vec<[u64; 3]> {
    let mut out: Vec<[u64; 3]> = positions
        .iter()
        .map(|p| {
            let q = match g {
                Some(e) => e.apply(*p),
                None => *p,
            };
            std::array::from_fn(|k| comparison_key(q[k]))
        })
        .collect();
    out.sort_unstable();
    out
}

fn triangle_keys(
    positions: &[[f64; 3]],
    indices: &[u32],
    g: Option<Element>,
) -> Vec<[[u64; 3]; 3]> {
    let mapped = |i: u32| -> [u64; 3] {
        let p = positions[i as usize];
        let q = match g {
            Some(e) => e.apply(p),
            None => p,
        };
        std::array::from_fn(|k| comparison_key(q[k]))
    };
    let mut out: Vec<[[u64; 3]; 3]> = Vec::with_capacity(indices.len() / 3);
    for tri in indices.as_chunks::<3>().0 {
        let c = [mapped(tri[0]), mapped(tri[1]), mapped(tri[2])];
        let rotations = [[c[0], c[1], c[2]], [c[1], c[2], c[0]], [c[2], c[0], c[1]]];
        out.push(rotations.into_iter().min().expect("three rotations"));
    }
    out.sort_unstable();
    out
}

// ─── fixtures ───────────────────────────────────────────────────────────────

#[derive(Clone, Copy)]
struct Fixture {
    samples: u32,
    cell_size: f64,
    origin: f64,
}

/// `P-57`'s two fixtures, unchanged and reached through `P-61` and `P-101`.
///
/// The 25³ arm uses `3L/32` rather than the crate's `2L/24`, because `L/12` is
/// not dyadic and its grid does not mirror.
fn fixtures(l: f64) -> [Fixture; 2] {
    let h33 = l / 16.0;
    let h25 = 3.0 * l / 32.0;
    [
        Fixture {
            samples: 33,
            cell_size: h33,
            origin: -l,
        },
        Fixture {
            samples: 25,
            cell_size: h25,
            origin: -12.0 * h25,
        },
    ]
}

/// Is the grid bit-exactly closed under a sign flip?
///
/// Without this the relation would be falsified by the fixture rather than by
/// the extractor, which is `P-57`'s own precondition. Recorded per row as
/// `grid_mirrors`, because the golden 25³ configuration does **not** mirror and
/// a reader has to be able to see that from the file.
fn grid_mirrors(samples: u32, origin: f64, cell_size: f64) -> bool {
    let coords: Vec<f64> = (0..samples)
        .map(|i| origin + cell_size * f64::from(i))
        .collect();
    let bits: Vec<u64> = coords.iter().map(|c| comparison_key(*c)).collect();
    coords.iter().all(|c| bits.contains(&comparison_key(-*c)))
}

/// The three resolutions `src/golden.rs` hashes every field at.
const GOLDEN_RESOLUTIONS: [u32; 3] = [17, 25, 33];

// ─── the committed artefacts, as the before-arms ────────────────────────────

/// One equivariance row of `docs/experiments/p-101.csv` — the shipped
/// extractor's reading, which is what `✗79 / M-413` quotes.
struct P101Row {
    elements_vertex_exact: usize,
    pure_permutation_exact: usize,
    pure_sign_flip_exact: usize,
}

/// Read `p-101.csv`'s equivariance block. **Not optional**: it is the vacuity
/// control's source and the check that this harness is `P-101`'s instrument.
fn p101_baseline() -> HashMap<String, P101Row> {
    let path = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../docs/experiments/p-101.csv");
    let text = std::fs::read_to_string(&path)
        .unwrap_or_else(|e| panic!("p-101.csv is the baseline and must be readable: {e}"));
    let mut lines = text.lines().filter(|l| !l.starts_with('#'));
    let header: Vec<&str> = lines
        .next()
        .expect("p-101.csv has a header")
        .split(',')
        .collect();
    let column = |name: &str| -> usize {
        header
            .iter()
            .position(|h| *h == name)
            .unwrap_or_else(|| panic!("p-101.csv has no `{name}` column"))
    };
    let (c_block, c_field, c_extractor, c_resolution) = (
        column("block"),
        column("field"),
        column("extractor"),
        column("resolution"),
    );
    let (c_exact, c_perm, c_sign) = (
        column("elements_vertex_exact_baseline"),
        column("pure_permutation_exact_baseline"),
        column("pure_sign_flip_exact_baseline"),
    );

    let mut out = HashMap::new();
    for line in lines {
        let f: Vec<&str> = line.split(',').collect();
        if f.len() <= c_sign.max(c_exact).max(c_perm) || f[c_block] != "equivariance" {
            continue;
        }
        let parse = |s: &str, what: &str| -> usize {
            s.parse()
                .unwrap_or_else(|_| panic!("p-101.csv `{what}` is not a count: {s:?}"))
        };
        out.insert(
            format!("{}/{}/{}", f[c_field], f[c_extractor], f[c_resolution]),
            P101Row {
                elements_vertex_exact: parse(f[c_exact], "elements_vertex_exact_baseline"),
                pure_permutation_exact: parse(f[c_perm], "pure_permutation_exact_baseline"),
                pure_sign_flip_exact: parse(f[c_sign], "pure_sign_flip_exact_baseline"),
            },
        );
    }
    assert_eq!(
        out.len(),
        32,
        "p-101.csv must carry 32 distinct dual equivariance configurations; found {}",
        out.len()
    );
    out
}

/// Pull one value out of a `golden_hashes.json` line.
///
/// A hand-rolled scanner rather than a JSON parser, for the reason `golden.rs`
/// gives: the grammar is one line, fixed key order, no nesting and no escapes.
fn json_field<'a>(line: &'a str, key: &str) -> Option<&'a str> {
    let needle = format!("\"{key}\":");
    let start = line.find(&needle)? + needle.len();
    let rest = line[start..].trim_start();
    if let Some(q) = rest.strip_prefix('"') {
        q.find('"').map(|end| &q[..end])
    } else {
        let end = rest.find([',', '}']).unwrap_or(rest.len());
        Some(rest[..end].trim())
    }
}

/// The committed golden hashes for the two duals, keyed `field/extractor/samples`.
fn golden_baseline() -> HashMap<String, u64> {
    let path = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("golden_hashes.json");
    let text = std::fs::read_to_string(&path)
        .unwrap_or_else(|e| panic!("golden_hashes.json is C3's baseline: {e}"));
    let mut out = HashMap::new();
    for line in text.lines() {
        let Some(algorithm) = json_field(line, "algorithm") else {
            continue;
        };
        if algorithm != "dual_contouring" && algorithm != "manifold_dual_contouring" {
            continue;
        }
        let field = json_field(line, "field").expect("a golden row names its field");
        let samples = json_field(line, "samples").expect("a golden row names its resolution");
        let hash = json_field(line, "hash").expect("a golden row carries a hash");
        out.insert(
            format!("{field}/{algorithm}/{samples}"),
            u64::from_str_radix(hash, 16).expect("a golden hash is hex"),
        );
    }
    assert_eq!(
        out.len(),
        48,
        "golden_hashes.json must carry 48 dual rows (8 fields x 3 resolutions x 2); found {}",
        out.len()
    );
    out
}

// ─── the equivariance measurement ───────────────────────────────────────────

struct Measured {
    vertices: usize,
    triangles: usize,
    vertex_exact: usize,
    triangle_exact: usize,
    pure_permutation_exact: usize,
    pure_sign_flip_exact: usize,
    first_failing_element: String,
    vertex_failing_labels: String,
    wall_ms: u128,
}

/// The two meshes [`measure`] needs: the reference build, and the one it
/// overwrites for each of the 48 group elements.
///
/// One struct rather than two parameters because the callers below would
/// otherwise carry eight arguments each, and eight arguments is where a caller
/// starts passing them in the wrong order.
struct Bufs {
    reference: MeshBuffer<f64>,
    rotated: MeshBuffer<f64>,
}

fn measure<S: Sdf<Scalar = f64>, M: Mesher>(
    field: &S,
    mesher: &mut M,
    site: &Site,
    elements: &[Element],
    bufs: &mut Bufs,
) -> Measured {
    let shape = RuntimeShape3::new([site.samples; 3]).expect("fixture grid fits u32");
    let started = Instant::now();

    mesher.mesh(
        field,
        &shape,
        site.origin,
        site.cell_size,
        &mut bufs.reference,
    );

    let mut m = Measured {
        vertices: bufs.reference.positions.len(),
        triangles: bufs.reference.triangle_count(),
        vertex_exact: 0,
        triangle_exact: 0,
        pure_permutation_exact: 0,
        pure_sign_flip_exact: 0,
        first_failing_element: String::from("none"),
        vertex_failing_labels: String::new(),
        wall_ms: 0,
    };
    let mut failing: Vec<String> = Vec::new();
    let mut any_failed = false;

    for (index, &g) in elements.iter().enumerate() {
        let wrapped = Rotated {
            field,
            g,
            g_inv: g.inverse(),
        };
        mesher.mesh(
            &wrapped,
            &shape,
            site.origin,
            site.cell_size,
            &mut bufs.rotated,
        );

        let got = vertex_keys(&bufs.rotated.positions, None);
        let want = vertex_keys(&bufs.reference.positions, Some(g));
        let vertex_ok = got == want;

        let got_tri = triangle_keys(&bufs.rotated.positions, &bufs.rotated.indices, None);
        let want_tri = triangle_keys(&bufs.reference.positions, &bufs.reference.indices, Some(g));
        let triangle_ok = got_tri == want_tri;

        if vertex_ok {
            m.vertex_exact += 1;
            if g.sign == [1, 1, 1] {
                m.pure_permutation_exact += 1;
            }
            if g.perm == [0, 1, 2] {
                m.pure_sign_flip_exact += 1;
            }
        } else {
            if !any_failed {
                m.first_failing_element = g.label();
                any_failed = true;
            }
            failing.push(g.short());
        }
        if triangle_ok {
            m.triangle_exact += 1;
        }

        // **The control**: element 0 is the identity and must reproduce the
        // reference exactly, or the extractor is not deterministic and the row
        // measures nothing.
        assert!(
            index != 0 || (vertex_ok && triangle_ok),
            "the identity element is not exact: the extractor is not \
             deterministic, so this row measures nothing"
        );
    }

    m.vertex_failing_labels = failing.join("|");
    m.wall_ms = started.elapsed().as_millis();
    m
}

/// Are two meshes bit-identical in every position, normal and index?
fn bit_identical(a: &MeshBuffer<f64>, b: &MeshBuffer<f64>) -> bool {
    a.indices == b.indices
        && a.positions.len() == b.positions.len()
        && a.normals.len() == b.normals.len()
        && a.positions
            .iter()
            .zip(&b.positions)
            .all(|(p, q)| (0..3).all(|k| p[k].to_bits() == q[k].to_bits()))
        && a.normals
            .iter()
            .zip(&b.normals)
            .all(|(p, q)| (0..3).all(|k| p[k].to_bits() == q[k].to_bits()))
}

/// How many vertex positions differ in bits, index for index.
///
/// A change of reduction cannot move a sign classification, so the two arms walk
/// the same cells in the same order and the vertex arrays are index-aligned.
/// `counts_identical` is recorded beside this so a reader can see the alignment
/// held.
fn positions_moved(a: &MeshBuffer<f64>, b: &MeshBuffer<f64>) -> usize {
    a.positions
        .iter()
        .zip(&b.positions)
        .filter(|(p, q)| (0..3).any(|k| p[k].to_bits() != q[k].to_bits()))
        .count()
}

fn hausdorff<S: Sdf<Scalar = f64>>(
    mesh: &MeshBuffer<f64>,
    field: &S,
    samples: u32,
    origin: [f64; 3],
    cell_size: f64,
) -> f64 {
    let shape = RuntimeShape3::new([samples; 3]).expect("grid fits u32");
    let cfg = AccuracyConfig::from_cell_size(cell_size).expect("valid cell size");
    accuracy(&mesh.positions, &mesh.indices, field, &shape, origin, &cfg)
        .expect("accuracy")
        .symmetric_hausdorff()
}

fn copy_mesh(from: &MeshBuffer<f64>, to: &mut MeshBuffer<f64>) {
    to.reset();
    to.positions.extend_from_slice(&from.positions);
    to.normals.extend_from_slice(&from.normals);
    to.indices.extend_from_slice(&from.indices);
}

// ─── the accumulator's own test, before any mesh is built ───────────────────

/// A seeded xorshift64\*, `experiment_p9`'s, for the permutation battery.
struct Rng(u64);

impl Rng {
    const fn new(seed: u64) -> Self {
        Self(seed | 1)
    }

    fn next(&mut self) -> u64 {
        let mut x = self.0;
        x ^= x >> 12;
        x ^= x << 25;
        x ^= x >> 27;
        self.0 = x;
        x.wrapping_mul(0x2545_F491_4F6C_DD1D)
    }

    fn shuffle<T>(&mut self, slice: &mut [T]) {
        for i in (1..slice.len()).rev() {
            let j = (self.next() % (i as u64 + 1)) as usize;
            slice.swap(i, j);
        }
    }
}

/// Term sets whose exact sum is a representable `f64`, with that sum.
///
/// Each is a case an `f64` running sum gets wrong, so the battery can fail.
#[allow(clippy::type_complexity)]
fn exact_cases() -> Vec<(Vec<f64>, f64)> {
    let two53 = 9_007_199_254_740_992.0f64; // 2^53
    vec![
        // Neal's own motivating shape: a huge term, a small one, the huge one
        // cancelled. A running sum loses the 1.0; so does `sum_equivariant`.
        (vec![1e100, 1.0, -1e100], 1.0),
        (vec![1e100, 1.0, -1e100, 1e-100, -1e-100], 1.0),
        // Ties-to-even at the 53-bit boundary, twice.
        (vec![two53, 1.0, 1.0], two53 + 2.0),
        (
            vec![1.0, f64::EPSILON / 2.0, f64::EPSILON / 2.0],
            1.0 + f64::EPSILON,
        ),
        // Subnormals, and the subnormal-to-normal boundary.
        (
            vec![f64::from_bits(1), f64::from_bits(1)],
            f64::from_bits(2),
        ),
        (vec![f64::from_bits(1); 2], f64::from_bits(2)),
        // An intermediate an `f64` cannot hold at all: `MAX + MAX` overflows,
        // the exact sum does not.
        (vec![f64::MAX, f64::MAX, -f64::MAX], f64::MAX),
        (vec![f64::MAX, f64::MAX, -f64::MAX, -f64::MAX], 0.0),
        // Exact cancellation, and the sign of the zero it produces.
        (vec![1e30, -1e30], 0.0),
        // Twelve terms, which is `EDGE_COUNT` and the depth a dual cell reaches.
        (
            vec![
                1e100, -1e100, 1e50, -1e50, 1e20, -1e20, two53, 1.0, 1.0, -two53, 0.5, 0.5,
            ],
            3.0,
        ),
    ]
}

/// Integer term sets whose exact sum has magnitude under `2^53`, so `as f64` is
/// an **exact** oracle with no rounding of its own.
///
/// Every term is a power of two or a small integer, hence exactly representable,
/// and the sums are chosen small while the terms are not — which is exactly the
/// regime an `f64` running sum fails in.
const INTEGER_CASES: &[&[i64]] = &[
    &[1 << 62, 1, -(1 << 62)],
    &[1 << 62, 1 << 62, 1, -(1 << 62), -(1 << 62)],
    &[1 << 53, 1, 1, -(1 << 53)],
    &[1 << 60, -(1 << 59), -(1 << 59), 7],
    &[
        1 << 62,
        1 << 40,
        1 << 20,
        1,
        -(1 << 62),
        -(1 << 40),
        -(1 << 20),
    ],
];

/// What the accumulator's own battery found.
struct SelfTest {
    exact_cases: usize,
    integer_oracle_cases: usize,
    permutation_cases: usize,
    negation_cases: usize,
    /// Cases where `sum_equivariant` — the shipped reduction — gets a different
    /// answer. **The battery's own vacuity control**: a battery every arm passes
    /// proves nothing about exactness.
    cases_where_shipped_differs: usize,
    worst_chunk_bits: u32,
}

/// Test the accumulator before it is used to measure anything.
fn self_test() -> SelfTest {
    let mut out = SelfTest {
        exact_cases: 0,
        integer_oracle_cases: 0,
        permutation_cases: 0,
        negation_cases: 0,
        cases_where_shipped_differs: 0,
        worst_chunk_bits: 0,
    };
    let mut rng = Rng::new(0x5EED_1118);

    for (terms, want) in exact_cases() {
        let mut worst = 0u32;
        let got = superacc_sum::<true>(&terms, &mut worst);
        out.worst_chunk_bits = out.worst_chunk_bits.max(worst);
        assert_eq!(
            got.to_bits(),
            want.to_bits(),
            "the superaccumulator is not exact on {terms:?}: got {got:?}, want {want:?}"
        );
        out.exact_cases += 1;

        // Order independence, over shuffles of the same multiset.
        let mut shuffled = terms.clone();
        for _ in 0..64 {
            rng.shuffle(&mut shuffled);
            let again = superacc_sum::<true>(&shuffled, &mut worst);
            assert_eq!(
                again.to_bits(),
                got.to_bits(),
                "the superaccumulator moved under a permutation of {terms:?}"
            );
        }
        out.permutation_cases += 1;

        // Negation equivariance, exactly, except that the two encodings of zero
        // are not negations of each other (`✗39`).
        let negated: Vec<f64> = terms.iter().map(|t| -t).collect();
        let neg_got = superacc_sum::<true>(&negated, &mut worst);
        if got == 0.0 {
            assert!(
                neg_got == 0.0,
                "the negated {terms:?} does not sum to zero: {neg_got:?}"
            );
        } else {
            assert_eq!(
                neg_got.to_bits(),
                (-got).to_bits(),
                "the superaccumulator is not exactly negation-equivariant on {terms:?}"
            );
        }
        out.negation_cases += 1;

        // The battery has to be able to tell the arms apart.
        let mut padded = [0.0f64; EDGE_COUNT];
        if terms.len() <= EDGE_COUNT {
            padded[..terms.len()].copy_from_slice(&terms);
            let mask = ((1u32 << terms.len()) - 1) as u16;
            let mut ignored = 0u32;
            let shipped = reduce_slots(&padded, mask, Reduce::EdgeSlots, &mut ignored);
            if shipped.to_bits() != got.to_bits() {
                out.cases_where_shipped_differs += 1;
            }
        }
    }

    for case in INTEGER_CASES {
        let exact: i128 = case.iter().map(|t| i128::from(*t)).sum();
        assert!(
            exact.abs() < 1i128 << 53,
            "the integer oracle only works when the sum is under 2^53: {exact}"
        );
        let terms: Vec<f64> = case
            .iter()
            .map(|t| {
                let f = *t as f64;
                assert_eq!(f as i64, *t, "{t} is not exactly representable as f64");
                f
            })
            .collect();
        let mut worst = 0u32;
        let got = superacc_sum::<true>(&terms, &mut worst);
        out.worst_chunk_bits = out.worst_chunk_bits.max(worst);
        assert_eq!(
            got.to_bits(),
            (exact as f64).to_bits(),
            "the superaccumulator disagrees with the i128 oracle on {case:?}: \
             got {got:?}, want {exact}"
        );
        out.integer_oracle_cases += 1;
    }

    assert!(
        out.exact_cases > 0
            && out.integer_oracle_cases > 0
            && out.permutation_cases > 0
            && out.negation_cases > 0,
        "the accumulator battery is empty, so its zero failures measure nothing"
    );
    assert!(
        out.cases_where_shipped_differs > 0,
        "no case distinguishes the superaccumulator from `sum_equivariant`, so the \
         battery cannot tell the arms apart and its passes are vacuous"
    );
    out
}

// ─── the cost block ─────────────────────────────────────────────────────────

/// Every active cell's crossings on one grid, built once and solved by both
/// cost arms — so the crossing construction and the normal normalisation are
/// outside the counted window.
fn cost_corpus<S: Sdf<Scalar = f64>>(field: &S, fx: &Fixture) -> Vec<CellCrossings> {
    let n = fx.samples as usize;
    let origin = [fx.origin; 3];
    let mut values = vec![0.0f64; n * n * n];
    for z in 0..n {
        for y in 0..n {
            for x in 0..n {
                values[x + n * (y + n * z)] = field.sample([
                    origin[0] + fx.cell_size * x as f64,
                    origin[1] + fx.cell_size * y as f64,
                    origin[2] + fx.cell_size * z as f64,
                ]);
            }
        }
    }
    let mut out = Vec::new();
    for z in 0..n - 1 {
        for y in 0..n - 1 {
            for x in 0..n - 1 {
                let mut corner = [0.0f64; 8];
                for (c, slot) in corner.iter_mut().enumerate() {
                    let o = corner_offset(c as u8);
                    *slot = values
                        [(x + o[0] as usize) + n * ((y + o[1] as usize) + n * (z + o[2] as usize))];
                }
                let inside = corner.iter().filter(|v| is_inside(**v)).count();
                if inside == 0 || inside == 8 {
                    continue;
                }
                let cell_origin =
                    cell_origin_of([x as u32, y as u32, z as u32], origin, fx.cell_size);
                out.push(crossings_of(
                    field,
                    &corner,
                    cell_origin,
                    fx.cell_size,
                    Norm::Naive,
                ));
            }
        }
    }
    out
}

#[derive(Clone, Copy)]
struct Counted {
    cycles: f64,
    instructions: f64,
    nanos: f64,
}

fn median(mut values: Vec<f64>) -> f64 {
    values.sort_by(f64::total_cmp);
    values[values.len() / 2]
}

fn solve_pass(corpus: &[CellCrossings], cfg: Config) -> usize {
    let mut solved = 0usize;
    let mut sink = 0.0f64;
    let mut ignored = 0u32;
    for cell in corpus {
        if let Some(x) = solve_cell::<false>(cell, u16::MAX, LAMBDA, cfg, &mut ignored) {
            solved += 1;
            sink += x[0] + x[1] + x[2];
        }
    }
    black_box(sink);
    solved
}

fn measure_cost(probe: &mut Probe, corpus: &[CellCrossings], cfg: Config) -> (Counted, usize) {
    let started = Instant::now();
    let solved = solve_pass(corpus, cfg);
    let pass_ns = started.elapsed().as_nanos() as f64;
    let inner = ((TARGET_BATCH_NS / pass_ns.max(1.0)).ceil() as usize).clamp(1, MAX_INNER);

    let mut cycles = Vec::with_capacity(REPS);
    let mut instructions = Vec::with_capacity(REPS);
    let mut nanos = Vec::with_capacity(REPS);
    for _ in 0..REPS {
        probe.reset_and_enable();
        let started = Instant::now();
        for _ in 0..inner {
            black_box(solve_pass(corpus, cfg));
        }
        let elapsed = started.elapsed().as_nanos() as f64;
        probe.disable();
        let counts = probe.read();
        assert!(
            counts.worst_ratio() >= MIN_TIME_RATIO,
            "a counter was multiplexed ({:.4}), so its value is an extrapolation",
            counts.worst_ratio()
        );
        let scale = 1.0 / (inner * solved) as f64;
        cycles.push(counts.cycles.count as f64 * scale);
        instructions.push(counts.instructions.count as f64 * scale);
        nanos.push(elapsed * scale);
    }
    (
        Counted {
            cycles: median(cycles),
            instructions: median(instructions),
            nanos: median(nanos),
        },
        solved,
    )
}

// ─── the arms ───────────────────────────────────────────────────────────────

/// One of the five bench-local configurations, plus its CSV labels.
///
/// The column names are spelled out per arm rather than formatted at runtime
/// because `Run::record` takes `&'static str` keys, and manufacturing a
/// `'static` from a `String` is not a thing this repository does.
struct ArmSpec {
    label: &'static str,
    cfg: Config,
    /// Whether this arm must be bit-identical to the shipped extractor.
    is_replica: bool,
    /// Rows reaching `elements_vertex_exact = 48`, over the 32 equivariance rows.
    col_at_48: &'static str,
    /// Rows reaching `pure_permutation_exact = 6`.
    col_perm_6: &'static str,
    /// Rows reaching `pure_sign_flip_exact = 8`.
    col_sign_8: &'static str,
    /// Golden dual hashes moved, of 48.
    col_golden: &'static str,
    /// Golden dual rows whose triangle count moved, of 48.
    col_triangles: &'static str,
}

const ARMS: [ArmSpec; 5] = [
    ArmSpec {
        label: "edge_slot",
        cfg: Config {
            centroid: Reduce::EdgeSlots,
            solve: Reduce::EdgeSlots,
            norm: Norm::Naive,
        },
        is_replica: true,
        col_at_48: "arm_edge_slot_rows_at_48",
        col_perm_6: "arm_edge_slot_rows_at_perm_6",
        col_sign_8: "arm_edge_slot_rows_at_sign_8",
        col_golden: "arm_edge_slot_golden_hashes_moved",
        col_triangles: "arm_edge_slot_golden_triangle_counts_moved",
    },
    ArmSpec {
        label: "ordered_naive",
        cfg: Config {
            centroid: Reduce::Ordered,
            solve: Reduce::Ordered,
            norm: Norm::Naive,
        },
        is_replica: false,
        col_at_48: "arm_ordered_naive_rows_at_48",
        col_perm_6: "arm_ordered_naive_rows_at_perm_6",
        col_sign_8: "arm_ordered_naive_rows_at_sign_8",
        col_golden: "arm_ordered_naive_golden_hashes_moved",
        col_triangles: "arm_ordered_naive_golden_triangle_counts_moved",
    },
    ArmSpec {
        label: "superaccumulator_solve_nine",
        cfg: Config {
            centroid: Reduce::EdgeSlots,
            solve: Reduce::Superacc,
            norm: Norm::Naive,
        },
        is_replica: false,
        col_at_48: "arm_superaccumulator_solve_nine_rows_at_48",
        col_perm_6: "arm_superaccumulator_solve_nine_rows_at_perm_6",
        col_sign_8: "arm_superaccumulator_solve_nine_rows_at_sign_8",
        col_golden: "arm_superaccumulator_solve_nine_golden_hashes_moved",
        col_triangles: "arm_superaccumulator_solve_nine_golden_triangle_counts_moved",
    },
    ArmSpec {
        label: "superaccumulator_all_twelve",
        cfg: Config {
            centroid: Reduce::Superacc,
            solve: Reduce::Superacc,
            norm: Norm::Naive,
        },
        is_replica: false,
        col_at_48: "arm_superaccumulator_all_twelve_rows_at_48",
        col_perm_6: "arm_superaccumulator_all_twelve_rows_at_perm_6",
        col_sign_8: "arm_superaccumulator_all_twelve_rows_at_sign_8",
        col_golden: "arm_superaccumulator_all_twelve_golden_hashes_moved",
        col_triangles: "arm_superaccumulator_all_twelve_golden_triangle_counts_moved",
    },
    ArmSpec {
        label: "superaccumulator_all_twelve_and_normal",
        cfg: Config {
            centroid: Reduce::Superacc,
            solve: Reduce::Superacc,
            norm: Norm::Superacc,
        },
        is_replica: false,
        col_at_48: "arm_superaccumulator_all_twelve_and_normal_rows_at_48",
        col_perm_6: "arm_superaccumulator_all_twelve_and_normal_rows_at_perm_6",
        col_sign_8: "arm_superaccumulator_all_twelve_and_normal_rows_at_sign_8",
        col_golden: "arm_superaccumulator_all_twelve_and_normal_golden_hashes_moved",
        col_triangles: "arm_superaccumulator_all_twelve_and_normal_golden_triangle_counts_moved",
    },
];

/// The arm C1, C2 and C3 are scored on: every reduction the registration names,
/// exact.
const REGISTERED_ARM: &str = "superaccumulator_all_twelve";

/// C2's denominator, and Neal's own comparand.
const ORDERED_ARM: &str = "ordered_naive";

/// The shipped form, priced beside the two cost arms so a reader can see what an
/// extraction pays today.
const SHIPPED_ARM: &str = "edge_slot";

/// The two dual extractors, by the name `p-101.csv` and `golden_hashes.json` use.
const EXTRACTORS: [&str; 2] = ["dual_contouring", "manifold_dual_contouring"];

/// Every row, before the aggregates are known.
type Row = Vec<(&'static str, String)>;

/// One measured configuration, before its arms run.
struct Site {
    field: &'static str,
    block: &'static str,
    samples: u32,
    origin: [f64; 3],
    cell_size: f64,
}

/// One arm's reading on one site, through whichever dual the site names.
fn arm_measure<S: Sdf<Scalar = f64>>(
    field: &S,
    extractor: &str,
    cfg: Config,
    site: &Site,
    elements: &[Element],
    bufs: &mut Bufs,
) -> Measured {
    if extractor == "dual_contouring" {
        let mut m = DualContouring::<f64, DcRule>::with_rule(DcRule { cfg });
        measure(field, &mut m, site, elements, bufs)
    } else {
        let mut m = DualContouring::<f64, CycleRule>::with_rule(CycleRule { cfg });
        measure(field, &mut m, site, elements, bufs)
    }
}

/// The shipped extractor's reading on one site. The vacuity control's source.
fn shipped_measure<S: Sdf<Scalar = f64>>(
    field: &S,
    extractor: &str,
    site: &Site,
    elements: &[Element],
    bufs: &mut Bufs,
) -> Measured {
    if extractor == "dual_contouring" {
        let mut m = DualContouring::<f64>::new();
        measure(field, &mut m, site, elements, bufs)
    } else {
        let mut m = ManifoldDualContouring::<f64>::new();
        measure(field, &mut m, site, elements, bufs)
    }
}

fn main() {
    if !std::env::args().any(|arg| arg == "--bench") {
        return;
    }
    #[cfg(not(target_os = "linux"))]
    {
        eprintln!(
            "P-118 scores C2 in retired instructions from perf_event_open, which is Linux \
             only. Refusing rather than recording a zero for a column that was not measured."
        );
        std::process::exit(1);
    }
    let prereg = isomesh::experiment!("P-118");
    common::experiment::run(prereg, |run| {
        // ── the accumulator's own battery, before anything is measured ───────
        let battery = self_test();
        println!(
            "-- the accumulator, tested before it is used --\n\
             --   exact cases {}, i128-oracle cases {}, permutation cases {} (64 shuffles each),\n\
             --   negation cases {}, cases where `sum_equivariant` disagrees {},\n\
             --   worst |chunk| {} bits of 64 (headroom: the module doc's 12 * 2^53 < 2^57)\n",
            battery.exact_cases,
            battery.integer_oracle_cases,
            battery.permutation_cases,
            battery.negation_cases,
            battery.cases_where_shipped_differs,
            battery.worst_chunk_bits,
        );

        println!(
            "-- what p-101.csv already fixes, quoted before this run --\n\
             --   the shipped duals reach elements_vertex_exact = 48 on {BASELINE_ROWS_AT_48} of 32 rows;\n\
             --   `edge_slot_equivariant_normal` reaches pure_permutation_exact = 6 on 32 of 32,\n\
             --   so the permutation obstruction is vec3::length's naive dot and not the accumulation;\n\
             --   pure_sign_flip_exact = 8 reads 6 of 32 for the shipped arm and for ALL FOUR of\n\
             --   P-101's arms, which is M-177. C1 needs both halves, so C1 is out of reach\n\
             --   before the numbers below are read. It is run anyway, and the fifth arm locates it.\n"
        );

        let elements = group();
        verify_group(&elements);
        let p101 = p101_baseline();
        let golden = golden_baseline();

        let mut rows: Vec<Row> = Vec::new();
        let mut bufs = Bufs {
            reference: MeshBuffer::<f64>::new(),
            rotated: MeshBuffer::<f64>::new(),
        };
        let mut shipped_mesh = MeshBuffer::<f64>::new();
        let mut arm_mesh = MeshBuffer::<f64>::new();

        let mut rows_at_48: HashMap<&str, usize> = HashMap::new();
        let mut rows_at_perm_6: HashMap<&str, usize> = HashMap::new();
        let mut rows_at_sign_8: HashMap<&str, usize> = HashMap::new();
        let mut golden_hashes_moved: HashMap<&str, usize> = HashMap::new();
        let mut golden_hashes_moved_expected: HashMap<&str, usize> = HashMap::new();
        let mut golden_triangles_moved: HashMap<&str, usize> = HashMap::new();
        let mut baseline_rows_at_48 = 0usize;
        let mut equivariance_rows = 0usize;
        let mut golden_rows = 0usize;
        let mut replica_bit_identical = true;
        let mut baseline_matches_p101 = true;
        let mut golden_fixture_matches_shipped = true;
        let mut topology_identical = true;
        let mut worst_hausdorff_delta = 0.0f64;
        let mut worst_chunk_bits = battery.worst_chunk_bits;

        println!(
            "{:<14} {:>4} {:<6} {:<25} {:<40} {:>4} {:>5} {:>3} {:>3} {:>7} {:>6}",
            "field", "n", "block", "extractor", "rule", "vex", "base", "p6", "f8", "moved", "ms"
        );

        isomesh::for_each_reference_field!(f64, |field_name, field| {
            let (lo, hi) = field.domain();
            let l = hi[0];
            for k in 0..3 {
                assert_eq!(
                    lo[k].to_bits(),
                    (-hi[k]).to_bits(),
                    "{field_name}: domain is not the symmetric cube the fixtures assume"
                );
            }

            // The two blocks' configurations, in one list so every row is
            // produced by the same code path.
            let mut sites: Vec<Site> = Vec::new();
            for fx in fixtures(l) {
                assert!(
                    grid_mirrors(fx.samples, fx.origin, fx.cell_size),
                    "{field_name} at {}: P-57's fixture must mirror bit-exactly",
                    fx.samples
                );
                sites.push(Site {
                    field: field_name,
                    block: "equivariance",
                    samples: fx.samples,
                    origin: [fx.origin; 3],
                    cell_size: fx.cell_size,
                });
            }
            for samples in GOLDEN_RESOLUTIONS {
                sites.push(Site {
                    field: field_name,
                    block: "golden",
                    samples,
                    origin: lo,
                    cell_size: (hi[0] - lo[0]) / f64::from(samples - 1),
                });
            }

            for site in &sites {
                let mirrors = grid_mirrors(site.samples, site.origin[0], site.cell_size);
                for extractor in EXTRACTORS {
                    let shipped = shipped_measure(&field, extractor, site, &elements, &mut bufs);
                    copy_mesh(&bufs.reference, &mut shipped_mesh);
                    let shipped_hash = mesh_hash(&shipped_mesh);

                    let key = format!("{}/{extractor}/{}", site.field, site.samples);
                    let committed = if site.block == "golden" {
                        let c = *golden
                            .get(&key)
                            .unwrap_or_else(|| panic!("golden_hashes.json has no {key}"));
                        // Without this, `hashes_moved` would be measured against
                        // a stale fixture and could not be read as a cost.
                        if c != shipped_hash {
                            golden_fixture_matches_shipped = false;
                        }
                        Some(c)
                    } else {
                        None
                    };

                    if site.block == "equivariance" {
                        // **The instrument check.** This harness re-implements
                        // P-57's group, fixtures and comparison keys; if its
                        // shipped arm disagrees with p-101.csv on the same
                        // configuration, it is not P-101's instrument and
                        // nothing below means anything.
                        let b = p101
                            .get(&key)
                            .unwrap_or_else(|| panic!("p-101.csv has no dual row {key}"));
                        let matched = shipped.vertex_exact == b.elements_vertex_exact
                            && shipped.pure_permutation_exact == b.pure_permutation_exact
                            && shipped.pure_sign_flip_exact == b.pure_sign_flip_exact;
                        assert!(
                            matched,
                            "{key}: this harness disagrees with p-101.csv about the SHIPPED \
                             extractor (vertex_exact {} vs {}, pure_perm {} vs {}, pure_sign \
                             {} vs {}) -- the instrument drifted",
                            shipped.vertex_exact,
                            b.elements_vertex_exact,
                            shipped.pure_permutation_exact,
                            b.pure_permutation_exact,
                            shipped.pure_sign_flip_exact,
                            b.pure_sign_flip_exact
                        );
                        baseline_matches_p101 &= matched;
                        equivariance_rows += 1;
                        if shipped.vertex_exact == GROUP_ORDER {
                            baseline_rows_at_48 += 1;
                        }
                    } else {
                        golden_rows += 1;
                    }

                    let shipped_hausdorff = hausdorff(
                        &shipped_mesh,
                        &field,
                        site.samples,
                        site.origin,
                        site.cell_size,
                    );

                    for arm in &ARMS {
                        let armed =
                            arm_measure(&field, extractor, arm.cfg, site, &elements, &mut bufs);
                        copy_mesh(&bufs.reference, &mut arm_mesh);
                        let identical = bit_identical(&arm_mesh, &shipped_mesh);
                        if arm.is_replica {
                            assert!(
                                identical,
                                "{key} / {}: the transcription is not the shipped arithmetic, \
                                 so every other arm's difference is unattributable",
                                arm.label
                            );
                            replica_bit_identical &= identical;
                        }

                        let moved = positions_moved(&arm_mesh, &shipped_mesh);
                        let counts_same = arm_mesh.positions.len() == shipped_mesh.positions.len();
                        let indices_same = arm_mesh.indices == shipped_mesh.indices;
                        let triangles_same =
                            arm_mesh.triangle_count() == shipped_mesh.triangle_count();
                        if !(counts_same && indices_same) {
                            topology_identical = false;
                        }
                        let arm_hash = mesh_hash(&arm_mesh);
                        let hash_moved = arm_hash != shipped_hash;
                        let arm_hausdorff =
                            hausdorff(&arm_mesh, &field, site.samples, site.origin, site.cell_size);
                        let delta = (arm_hausdorff - shipped_hausdorff).abs();
                        worst_hausdorff_delta = worst_hausdorff_delta.max(delta);

                        // C1 is scored on the equivariance block only: a row
                        // whose grid does not mirror would be falsified by the
                        // fixture rather than by the extractor.
                        if site.block == "equivariance" {
                            if armed.vertex_exact == GROUP_ORDER {
                                *rows_at_48.entry(arm.label).or_default() += 1;
                            }
                            if armed.pure_permutation_exact == PERMS.len() {
                                *rows_at_perm_6.entry(arm.label).or_default() += 1;
                            }
                            if armed.pure_sign_flip_exact == 8 {
                                *rows_at_sign_8.entry(arm.label).or_default() += 1;
                            }
                        } else {
                            if hash_moved {
                                *golden_hashes_moved.entry(arm.label).or_default() += 1;
                            }
                            if moved > 0 || !counts_same || !indices_same {
                                *golden_hashes_moved_expected.entry(arm.label).or_default() += 1;
                            }
                            if !triangles_same {
                                *golden_triangles_moved.entry(arm.label).or_default() += 1;
                            }
                        }

                        println!(
                            "{:<14} {:>4} {:<6} {:<25} {:<40} {:>4} {:>5} {:>3} {:>3} {:>7} {:>6}",
                            site.field,
                            site.samples,
                            site.block,
                            extractor,
                            arm.label,
                            armed.vertex_exact,
                            shipped.vertex_exact,
                            armed.pure_permutation_exact,
                            armed.pure_sign_flip_exact,
                            moved,
                            armed.wall_ms
                        );

                        rows.push(vec![
                            ("block", site.block.to_string()),
                            ("field", site.field.to_string()),
                            ("resolution", site.samples.to_string()),
                            ("extractor", extractor.to_string()),
                            ("rule", arm.label.to_string()),
                            ("cell_size", format!("{:.9}", site.cell_size)),
                            ("grid_mirrors", mirrors.to_string()),
                            ("elements_tested", GROUP_ORDER.to_string()),
                            ("elements_vertex_exact", armed.vertex_exact.to_string()),
                            (
                                "baseline_elements_vertex_exact",
                                shipped.vertex_exact.to_string(),
                            ),
                            ("elements_triangle_exact", armed.triangle_exact.to_string()),
                            (
                                "pure_permutation_exact",
                                armed.pure_permutation_exact.to_string(),
                            ),
                            (
                                "baseline_pure_permutation_exact",
                                shipped.pure_permutation_exact.to_string(),
                            ),
                            (
                                "pure_sign_flip_exact",
                                armed.pure_sign_flip_exact.to_string(),
                            ),
                            (
                                "baseline_pure_sign_flip_exact",
                                shipped.pure_sign_flip_exact.to_string(),
                            ),
                            ("first_failing_element", armed.first_failing_element.clone()),
                            ("vertex_failing_labels", armed.vertex_failing_labels.clone()),
                            ("vertices", armed.vertices.to_string()),
                            ("triangles", armed.triangles.to_string()),
                            ("baseline_vertices", shipped.vertices.to_string()),
                            ("baseline_triangles", shipped.triangles.to_string()),
                            ("positions_moved", moved.to_string()),
                            ("counts_identical", counts_same.to_string()),
                            ("indices_identical", indices_same.to_string()),
                            ("triangle_count_identical", triangles_same.to_string()),
                            ("bit_identical_to_shipped", identical.to_string()),
                            ("mesh_hash_shipped", format!("{shipped_hash:016x}")),
                            ("mesh_hash_arm", format!("{arm_hash:016x}")),
                            ("mesh_hash_moved", hash_moved.to_string()),
                            (
                                "golden_hash_committed",
                                committed.map_or_else(String::new, |c| format!("{c:016x}")),
                            ),
                            ("hausdorff_shipped", format!("{shipped_hausdorff:.12e}")),
                            ("hausdorff_arm", format!("{arm_hausdorff:.12e}")),
                            ("hausdorff_delta", format!("{delta:.12e}")),
                            ("wall_ms", armed.wall_ms.to_string()),
                        ]);
                    }
                }
            }
        });

        // ── C2: the cost of the reduction, on one pre-built corpus ───────────
        println!("\n-- cost: the same corpus solved by three reductions --");
        let mut probe = Probe::open();
        let mut cost: HashMap<&str, Counted> = HashMap::new();

        // Named directly rather than through `for_each_reference_field!`: that
        // macro inlines its body once per field, so a `return` in it would
        // return from this whole closure and a write to an outer local in one
        // copy is dead in the next.
        let cost_field = isomesh::fields::Sphere::<f64>::canonical();
        assert_eq!(COST_FIELD, "sphere", "the cost field's name is a column");
        let (_, cost_hi) = cost_field.domain();
        let cost_fixture = fixtures(cost_hi[0])[0];
        let corpus = cost_corpus(&cost_field, &cost_fixture);
        let corpus_cells = corpus.len();
        let corpus_crossings: usize = corpus.iter().map(|c| c.mask.count_ones() as usize).sum();
        assert!(
            corpus_cells > 0 && corpus_crossings > 0,
            "the cost corpus is empty, so its ratio measures nothing"
        );

        // The headroom census, outside the counted window, over exactly the
        // corpus the counted window solves.
        let mut census = 0u32;
        for cell in &corpus {
            let _ = solve_cell::<true>(
                cell,
                u16::MAX,
                LAMBDA,
                Config {
                    centroid: Reduce::Superacc,
                    solve: Reduce::Superacc,
                    norm: Norm::Naive,
                },
                &mut census,
            );
        }
        worst_chunk_bits = worst_chunk_bits.max(census);

        let mut corpus_solved = 0usize;
        for arm in &ARMS {
            if arm.label != REGISTERED_ARM && arm.label != ORDERED_ARM && arm.label != SHIPPED_ARM {
                continue;
            }
            let (counted, solved) = measure_cost(&mut probe, &corpus, arm.cfg);
            corpus_solved = solved;
            println!(
                "  {:<40} {:>10.3} instr  {:>10.3} cyc  {:>10.3} ns  per solve",
                arm.label, counted.instructions, counted.cycles, counted.nanos
            );
            cost.insert(arm.label, counted);
        }

        let ordered = *cost
            .get(ORDERED_ARM)
            .expect("the ordered arm was not priced");
        let superacc = *cost
            .get(REGISTERED_ARM)
            .expect("the registered arm was not priced");
        let shipped_cost = *cost
            .get(SHIPPED_ARM)
            .expect("the shipped arm was not priced");
        let cost_ratio = superacc.instructions / ordered.instructions;
        let ns_cost_ratio = superacc.nanos / ordered.nanos;
        let shipped_over_ordered = shipped_cost.instructions / ordered.instructions;
        // What C2's ratio is *made of*, so the entry can name a mechanism rather
        // than a number: one solve runs twelve reductions (3 centroid + 6 matrix
        // + 3 gradient) and each one pays a whole 67-chunk read-out, over a mean
        // of this many terms.
        let reductions_per_solve = 12usize;
        let mean_terms_per_reduction = corpus_crossings as f64 / corpus_solved.max(1) as f64;

        // ── the aggregates, and the clause verdicts ─────────────────────────
        let at_48 = *rows_at_48.get(REGISTERED_ARM).unwrap_or(&0);
        let hashes_moved = *golden_hashes_moved.get(REGISTERED_ARM).unwrap_or(&0);
        let hashes_moved_expected = *golden_hashes_moved_expected
            .get(REGISTERED_ARM)
            .unwrap_or(&0);
        let triangle_counts_moved = *golden_triangles_moved.get(REGISTERED_ARM).unwrap_or(&0);

        // C1's population is the rows the shipped solve is BELOW 48 on: the
        // registration says "reaches 48 of 48 ... where the shipped solve is
        // below it". Scoring it over all 32 would count the six rows that were
        // already at 48 as the mechanism's work.
        let c1_population = equivariance_rows - baseline_rows_at_48;
        let c1_lifted = at_48.saturating_sub(baseline_rows_at_48);
        let c1_holds = c1_population > 0 && c1_lifted == c1_population;
        let c2_holds = cost_ratio <= C2_BAR;
        let c3_holds = hashes_moved == hashes_moved_expected && triangle_counts_moved == 0;

        // **The vacuity control, asserted rather than merely reported.**
        assert_eq!(
            equivariance_rows, 32,
            "the equivariance block must carry 32 dual configurations"
        );
        assert_eq!(
            golden_rows, 48,
            "the golden block must carry 48 dual configurations"
        );
        assert_eq!(
            baseline_rows_at_48, BASELINE_ROWS_AT_48,
            "the shipped duals must reproduce p-101.csv's {BASELINE_ROWS_AT_48} of 32 rows at 48"
        );
        assert!(
            c1_population > 0,
            "every equivariance row is already at 48, so C1 could not have failed"
        );
        assert!(
            replica_bit_identical,
            "the edge_slot arm must be bit-identical to the shipped extractor"
        );
        assert!(
            baseline_matches_p101,
            "the shipped arm must reproduce p-101.csv row for row"
        );
        assert!(
            golden_fixture_matches_shipped,
            "the committed golden hashes must match the shipped extractor, or \
             `hashes_moved` is measured against a stale fixture"
        );
        assert!(
            worst_chunk_bits > 0,
            "the headroom census never saw a chunk, so its margin measures nothing"
        );
        assert!(
            worst_chunk_bits < 64,
            "a chunk reached {worst_chunk_bits} bits of 64: the no-carry-propagation \
             bound in the module doc is wrong and every superaccumulator number is suspect"
        );

        println!("\n-- aggregates --");
        println!(
            "  C1  rows at elements_vertex_exact = 48, registered arm: {at_48} of 32 \
             (shipped baseline {baseline_rows_at_48}); lifted {c1_lifted} of {c1_population}"
        );
        println!(
            "  C2  instructions per solve: ordered {:.3}, superacc {:.3}, ratio {cost_ratio:.4} \
             against a {C2_BAR:.1}x bar (shipped `sum_equivariant` {:.3})",
            ordered.instructions, superacc.instructions, shipped_cost.instructions
        );
        println!(
            "      ns per solve: ordered {:.3}, superacc {:.3}, ratio {ns_cost_ratio:.4} \
             -- REPORTED, NOT GATED (M-280)",
            ordered.nanos, superacc.nanos
        );
        println!(
            "      what the ratio is made of: {reductions_per_solve} reductions per solve, each \
             paying a whole {SCHUNKS}-chunk read-out over a mean of \
             {mean_terms_per_reduction:.3} terms; `sum_equivariant` itself costs \
             {shipped_over_ordered:.3}x the ordered sum"
        );
        println!(
            "  C3  golden hashes moved {hashes_moved} of 48, predicted from moved positions \
             {hashes_moved_expected}; triangle counts moved {triangle_counts_moved}"
        );
        println!("      topology identical everywhere: {topology_identical}");
        println!("      worst |hausdorff delta|:        {worst_hausdorff_delta:.6e}");
        println!("  verdicts: C1 {c1_holds}, C2 {c2_holds}, C3 {c3_holds}");
        for arm in &ARMS {
            println!(
                "  arm {:<40} 48 {:>2}/32  perm6 {:>2}/32  sign8 {:>2}/32  golden moved {:>2}/48  \
                 triangles moved {:>2}/48",
                arm.label,
                rows_at_48.get(arm.label).unwrap_or(&0),
                rows_at_perm_6.get(arm.label).unwrap_or(&0),
                rows_at_sign_8.get(arm.label).unwrap_or(&0),
                golden_hashes_moved.get(arm.label).unwrap_or(&0),
                golden_triangles_moved.get(arm.label).unwrap_or(&0),
            );
        }

        let mut aggregates: Vec<(&'static str, String)> = vec![
            // ── the registered cost columns ──
            ("ns_per_solve_ordered", format!("{:.4}", ordered.nanos)),
            ("ns_per_solve_superacc", format!("{:.4}", superacc.nanos)),
            ("cost_ratio", format!("{cost_ratio:.6}")),
            ("ns_cost_ratio", format!("{ns_cost_ratio:.6}")),
            ("c2_bar", format!("{C2_BAR:.1}")),
            (
                "instructions_per_solve_ordered",
                format!("{:.4}", ordered.instructions),
            ),
            (
                "instructions_per_solve_superacc",
                format!("{:.4}", superacc.instructions),
            ),
            (
                "instructions_per_solve_edge_slot",
                format!("{:.4}", shipped_cost.instructions),
            ),
            ("cycles_per_solve_ordered", format!("{:.4}", ordered.cycles)),
            (
                "cycles_per_solve_superacc",
                format!("{:.4}", superacc.cycles),
            ),
            (
                "cycles_per_solve_edge_slot",
                format!("{:.4}", shipped_cost.cycles),
            ),
            ("cost_corpus_field", COST_FIELD.to_string()),
            ("cost_corpus_cells", corpus_cells.to_string()),
            ("cost_corpus_solved", corpus_solved.to_string()),
            ("cost_corpus_crossings", corpus_crossings.to_string()),
            ("cost_reps", REPS.to_string()),
            (
                "shipped_over_ordered_instruction_ratio",
                format!("{shipped_over_ordered:.6}"),
            ),
            (
                "superacc_reductions_per_solve",
                reductions_per_solve.to_string(),
            ),
            (
                "superacc_mean_terms_per_reduction",
                format!("{mean_terms_per_reduction:.6}"),
            ),
            // ── the registered C3 columns ──
            ("hashes_moved", hashes_moved.to_string()),
            ("hashes_moved_expected", hashes_moved_expected.to_string()),
            ("triangle_counts_moved", triangle_counts_moved.to_string()),
            ("golden_dual_rows", golden_rows.to_string()),
            // ── the verdicts ──
            ("c1_holds", c1_holds.to_string()),
            ("c2_holds", c2_holds.to_string()),
            ("c3_holds", c3_holds.to_string()),
            ("c1_rows_at_48", at_48.to_string()),
            ("c1_population", c1_population.to_string()),
            ("c1_rows_lifted", c1_lifted.to_string()),
            ("equivariance_rows", equivariance_rows.to_string()),
            ("baseline_rows_at_48", baseline_rows_at_48.to_string()),
            // ── the instrument and the accumulator's own battery ──
            ("replica_bit_identical", replica_bit_identical.to_string()),
            ("baseline_matches_p101", baseline_matches_p101.to_string()),
            (
                "golden_fixture_matches_shipped",
                golden_fixture_matches_shipped.to_string(),
            ),
            ("topology_identical", topology_identical.to_string()),
            (
                "max_hausdorff_delta",
                format!("{worst_hausdorff_delta:.12e}"),
            ),
            ("superacc_chunks", SCHUNKS.to_string()),
            ("superacc_chunk_slice_bits", LOW_MANTISSA_BITS.to_string()),
            ("superacc_scale_exponent", SCALE.to_string()),
            ("worst_chunk_magnitude_bits", worst_chunk_bits.to_string()),
            ("superacc_exact_cases", battery.exact_cases.to_string()),
            (
                "superacc_integer_oracle_cases",
                battery.integer_oracle_cases.to_string(),
            ),
            (
                "superacc_permutation_cases",
                battery.permutation_cases.to_string(),
            ),
            (
                "superacc_negation_cases",
                battery.negation_cases.to_string(),
            ),
            (
                "superacc_cases_where_shipped_differs",
                battery.cases_where_shipped_differs.to_string(),
            ),
            (
                "c1_reachable_before_run",
                "false_p101_locates_the_obstruction_in_vec3_length_and_m177".to_string(),
            ),
        ];
        // Every arm's headline counts on every row, so a reader does not have to
        // filter the file to see whether the mechanism arm moved. That is
        // `p-59.csv`'s shape and `p-61.csv` and `p-101.csv` kept it.
        for arm in &ARMS {
            aggregates.push((
                arm.col_at_48,
                rows_at_48.get(arm.label).unwrap_or(&0).to_string(),
            ));
            aggregates.push((
                arm.col_perm_6,
                rows_at_perm_6.get(arm.label).unwrap_or(&0).to_string(),
            ));
            aggregates.push((
                arm.col_sign_8,
                rows_at_sign_8.get(arm.label).unwrap_or(&0).to_string(),
            ));
            aggregates.push((
                arm.col_golden,
                golden_hashes_moved.get(arm.label).unwrap_or(&0).to_string(),
            ));
            aggregates.push((
                arm.col_triangles,
                golden_triangles_moved
                    .get(arm.label)
                    .unwrap_or(&0)
                    .to_string(),
            ));
        }

        for mut row in rows {
            row.extend(aggregates.iter().cloned());
            run.record(&row);
        }
    });
}
