//! **P-167 — the Fourier spectrum of the 256-case table, and whether it is
//! low-degree.**
//!
//! Ticket: R-167. Pre-registered before this harness existed.
//!
//! ```bash
//! cargo bench --bench experiment_p167
//! ```
//!
//! Writes `docs/experiments/p-167.csv`.
//!
//! # What was missing
//!
//! Marching Cubes is a lookup on a Boolean function of eight variables. The case
//! table is `pub static CASES: [McCase; 256]`
//! (`crates/isomesh/src/marching_cubes/table.rs:180`), built by a `const fn` from
//! the corner signs and nothing else, and `corner_inside(case, corner)` is
//! literally `case & (1 << corner) != 0` — so `{0,1}^8 -> N` is not an analogy,
//! it is the type. The repository has audited that object three ways and never
//! once as a Boolean function:
//!
//! - `validate_table()` (`marching_cubes/mod.rs:836`) checks it **combinatorially**
//!   — cut-edge agreement, face consistency, manifoldness — which is a
//!   correctness gate and says nothing about the function's structure.
//! - `P-116` ran GRAPHGEN's pipeline over the same 256 entries and asked about
//!   *patterns*; `P-121` measured what the classification stage **costs**
//!   (`docs/experiments/p-121.csv`: `cycles_classify / cycles_total` runs from
//!   **1.8%** on `fbm_terrain` to **64.6%** on `sphere` at 65³, `marching_cubes`,
//!   `f32`). Neither asks what the function *is*.
//! - Nothing in the repository has ever taken a Walsh–Hadamard transform. The
//!   machinery is O'Donnell, *Analysis of Boolean Functions*, `arXiv:2105.10386`,
//!   Chapters 1–2 at `n = 8`; the corpus scored it 0.584 and it appears in no
//!   source file.
//!
//! So this row is the first reading of the table as a function rather than as a
//! list. It proposes no source change: the whole mechanism is
//! `benches/common/boolean.rs`, which **this ticket owns** and which R-168 and
//! R-169 consume unchanged.
//!
//! # The vector-valued problem, and how C1's falsifier is answered
//!
//! Walsh–Hadamard is defined for a **scalar** `+/-1`-valued function. The table's
//! output is an integer, so the object is *vector*-valued and C1's falsifier
//! anticipates exactly this: *"C1 by the transform not being well-defined for the
//! vector-valued output, which would mean the framing needs per-bit treatment and
//! should say so."*
//!
//! **This harness says so.** The framing is per output bit — that is what the
//! `output_bit` column is — and the claim that the per-bit family is a faithful
//! decomposition of the vector-valued function is *measured*, not asserted: for
//! every one of the 256 inputs, the integer recomposed from the per-bit
//! **spectral** evaluations is compared against the reading's own output integer,
//! and `reading_reconstruction_exact` reports the result. A per-bit spectrum that
//! did not recompose would be a decomposition of something else.
//!
//! # Three readings, because "the output" is a choice and a choice must be
//! measured
//!
//! `common::boolean` exposes three readings of `CASES` and all three are run
//! here, one row per (reading, bit):
//!
//! - **`shipped_triangle_counts` is the primary, and it is the one to quote.**
//!   It is where the table's triangulation *decisions* live, and it is
//!   octahedrally invariant — the module verifies 0 violations over all 48 × 256
//!   corner relabellings — which is what makes R-169's influence-equality check a
//!   real test. Maximum 5, so four bits.
//! - **`shipped_edge_masks`** is twelve two-corner parities by construction, so
//!   its spectrum is known in closed form *before* it is computed: one Fourier
//!   coefficient at degree 2 per bit. That makes it a third calibration beside
//!   parity and majority, at a degree neither of them occupies. It is a poor
//!   primary — a single edge bit is attached to two named corners and so is not
//!   cube-invariant.
//! - **`shipped_centroid_counts`** is the constant zero: `CASES` is the
//!   all-separate resolution and plain Marching Cubes never reaches a cycle long
//!   enough to need an interior vertex. A constant is a legitimate answer that
//!   says something real — the centroid machinery is dead code for the
//!   unambiguous table — and a bad primary, because a constant's spectrum cannot
//!   falsify C2 either way. Rows whose measured `role` is `constant` carry
//!   `is_degenerate = true` and are excluded from the run-level verdict.
//!
//! Each reading is analysed in `ceil(log2(max + 1)) + 1` bits. The `+ 1` is not
//! padding: the top bit is the **constant-zero witness**, whose measured degree 0
//! proves the reading's output really stops where the arithmetic says it stops.
//! For the triangle counts that is bit 3, which is the reading
//! `common::boolean`'s own documentation describes.
//!
//! # Arms
//!
//! | arm | what it varies | `is_control` |
//! |---|---|---|
//! | `triangle_counts.bit0` … `.bit3` | the primary reading, one row per output bit | no |
//! | `edge_masks.bit0` … `.bit12` | a reading whose spectrum is known before it is computed | no |
//! | `centroid_counts.bit0` | the constant-zero reading | no |
//! | `calibration.parity` | a known degree-8 function through the same transform | **yes** |
//! | `calibration.majority` | a known degree-1-heavy function through the same transform | **yes** |
//!
//! # The three clauses, and the thresholds they are read against
//!
//! Every threshold is a named constant so it cannot drift between the prose and
//! the row.
//!
//! **C1** — the spectrum is computed and its weight distribution by degree
//! reported. Measured per row as: Parseval holds (`weight_sum` within
//! [`NEGLIGIBLE`] of 1 — `Bool8::fourier` also asserts it), the bit's spectral
//! evaluation reproduces its truth table at all 256 inputs, and the reading's
//! integer recomposes from its bits. `fourier_weight_by_degree` is the nine
//! per-degree weights, `|`-joined at nine decimals because the CSV writer refuses
//! a comma.
//!
//! **C2** — the spectrum is **not** concentrated on low degrees.
//! `spectral_concentration` is the weight at or below degree
//! [`CONCENTRATION_DEGREE`] = 2, and "concentrated" means at least
//! [`WEIGHT_TARGET`] = 0.9 of the unit mass is there. Degree 2 is not an
//! arbitrary bar: a degree-2 multilinear form on eight variables already has
//! `1 + 8 + 28 = 37` coefficients, so a degree-3-or-higher "low-degree
//! approximation" costs more multiply-adds than the lookup it would replace and
//! the phrase stops meaning anything. `degree_for_90pct_weight` reports the
//! smallest degree that does reach 0.9, so the negative is quantitative at every
//! threshold rather than only at this one, and `anf_terms` gives the `GF(2)`
//! sparsity the registration asks for. **Note that ANF degree and Fourier degree
//! are different numbers** — an edge bit `x_a xor x_b` has ANF degree 1 and two
//! terms but Fourier degree 2 — so both are recorded.
//!
//! **C3** — *"if C2 is falsified and the spectrum is sparse, a branchless
//! spectral evaluation is benchmarked against the table lookup."* The antecedent
//! is recorded as `c3_reachable = concentrated && sparse`, and **both evaluations
//! are implemented and timed on every row regardless**, because a measured
//! `eval_ms_table` / `eval_ms_spectral` pair costs milliseconds and turns C3's
//! verdict from conditional into quantitative. `c3_holds` is
//! `eval_ms_spectral <= eval_ms_table` — the registration's own falsifier, *"C3
//! by the spectral evaluation losing to a lookup, which is likely — a 256-entry
//! table is already in L1"*. `branchless_feasible` is `sparse` **and** the
//! spectral form winning, which is the only combination under which a caller
//! should consider it.
//!
//! The spectral evaluation is genuinely branchless: `f(x) = 1` iff
//! `sum_S fhat(S) chi_S(x) < 0`, and the character's sign is
//! `1 - 2 * ((S & x).count_ones() & 1)` — arithmetic on a popcount, no test in
//! the inner loop. It sums only the non-negligible coefficients, which is what
//! makes sparsity the thing that decides its cost. `sparse` is `spectral_terms`
//! at or below [`SPARSE_TERM_LIMIT`] = 16: the classification stage
//! already performs eight corner-sign gathers per cell, so twice that many
//! multiply-adds is the outer boundary of "same cost class as one L1 load", and
//! anything past it is not a candidate.
//!
//! # Timing
//!
//! `std::time::Instant`, no criterion. [`WORKLOAD`] = 2²⁰ evaluations per pass
//! over a deterministic `SplitMix64` sequence of case indices seeded at
//! [`WORKLOAD_SEED`]; one warm-up pass of each path, then [`REPEATS`] = 7 timed
//! repeats with the two paths **interleaved** so both see the same governor
//! state — M-280 measured this host's `amd-pstate-epp` swinging the same binary
//! 1.45× between runs. Median is the headline, min and max are extra columns, and
//! `eval_scatter_table` / `eval_scatter_spectral` are max over min so a row whose
//! repeats disagree is visible rather than averaged into a pass. The accumulator
//! goes through `black_box`, and `eval_ms_table > 0` is a vacuity control rather
//! than a comment: a deleted loop would make C3's ratio a division by nothing.
//!
//! # SHARE, recomputed before the numbers
//!
//! Registered: *"C3 moves the classification stage, which the work-graphs result
//! put at a 2.8–3.4× loss when moved to irregular GPU dispatch"* (`V-12`,
//! `10.1145/3675376`, plus the independent profile at
//! `docs/research/2026-08-11-meshing-speed-analysis.md:117`).
//!
//! Discharged, and the arithmetic is the discharge. Classification is `s` of
//! extraction and a replacement `k`× slower makes extraction `1 - s + k*s`.
//! `P-121`'s own CSV puts `s` between **0.018** and **0.646** for
//! `marching_cubes` at 65³. So the SHARE is **do not move this stage**, at every
//! density the spectrum could have: the only question C3 can answer is by how
//! much, and it answers it with `eval_ms_ratio` on every row rather than with an
//! argument. The direction agrees with `V-12`'s independently measured loss,
//! which is worth noting precisely because the two arrive at it from unrelated
//! evidence.
//!
//! # Vacuity controls
//!
//! Every one runs before the first `run.record` and every panic starts `VOID: `.
//!
//! - **The instrument, on two functions with known spectra.**
//!   `common::boolean::self_check()` runs parity and majority-of-8 through the
//!   same code path the table uses; parity's entire weight at degree 8 is
//!   asserted inside it. Here its three returns are asserted against **closed
//!   forms computed in this file** — `1`, `8 * (70/256)^2` and `8 * 35/128` —
//!   rather than against transcribed decimals, so the check compares a
//!   measurement with arithmetic. Columns
//!   `calibration_parity_weight_degree8`, `calibration_majority_weight_degree1`,
//!   `calibration_majority_total_influence`.
//! - **The instrument must separate them.** Parity's weight at or below degree 2
//!   must be zero and majority's must exceed it by more than
//!   [`CALIBRATION_SEPARATION`] = 0.5. A transform that cannot tell a degree-8
//!   function from a degree-1-heavy one would report a plausible-looking spectrum
//!   for both. Both are also emitted as their own rows, so the separation is
//!   visible in the CSV and not only in an assert.
//! - **`sparse` must be a measurement and not a constant column** (M-44). At
//!   least one row of the census must be sparse and at least one must be dense,
//!   or the column could not have come out the other way. Columns
//!   `control_sparse_rows`, `control_dense_rows`.
//! - **The primary reading must be non-constant**, or C2's concentration is
//!   computed against a constant and its verdict is an artefact of the fixture.
//!   At least one bit of `shipped_triangle_counts` must have positive Fourier
//!   degree.
//! - **The transform must be sensitive to the table's contents.**
//!   `common::boolean::corrupt(&triangle_counts, 37, 2)` flips one bit of one
//!   case — the module measures that case and delta as one that moves the
//!   influences — and the L1 distance between the corrupt and shipped
//!   degree-weight vectors must be positive. Column
//!   `corrupt_control_weight_shift`.
//! - **Neither timed loop may be elided.** Both medians must be strictly
//!   positive on every row.
//!
//! # What R-168 and R-169 read out of this CSV
//!
//! R-168 computes noise stability from the spectrum and R-169 reads the
//! influences, so every row names its reading and its bit unambiguously in
//! `output_bit`, and carries `total_influence`, `influence_by_corner` (eight
//! values, `|`-joined), `influences_all_equal` and `corner_symmetry_classes` —
//! the last being the module's *generated* answer that the cube group is
//! transitive on the eight corners, which is what makes R-169's "equal within
//! each octahedral symmetry class" mean "all eight equal".

#![allow(
    clippy::cast_possible_truncation,
    clippy::cast_precision_loss,
    clippy::too_many_lines
)]

mod common;

use std::hint::black_box;
use std::time::Instant;

use crate::common::boolean::{
    Bool8, NEGLIGIBLE, corner_symmetry_classes, corrupt, self_check, shipped_centroid_counts,
    shipped_edge_masks, shipped_triangle_counts,
};

// ─── clause constants ───────────────────────────────────────────────────────

/// The degree at or below which weight counts as "low" for C2.
///
/// Two, and the reason is a cost rather than a convention: a degree-2 multilinear
/// form on eight variables has `1 + 8 + 28 = 37` coefficients, so a degree-3
/// approximation already costs more arithmetic than the single L1 load it would
/// be replacing and "low-degree approximation" stops naming anything a caller
/// could want.
const CONCENTRATION_DEGREE: u32 = 2;

/// The fraction of the unit Fourier mass that has to sit at or below
/// [`CONCENTRATION_DEGREE`] for the spectrum to count as concentrated there.
///
/// Also the target for `degree_for_90pct_weight`, so the two numbers are the same
/// question asked in the two directions: *how much weight is below degree 2*, and
/// *how high does the degree have to go to reach 0.9*.
const WEIGHT_TARGET: f64 = 0.9;

/// At most this many non-negligible Fourier coefficients, out of 256, for the
/// spectrum to count as sparse.
///
/// Sixteen. Classification already gathers eight corner signs per cell, so twice
/// that many multiply-adds is the outer edge of "the same cost class as one
/// lookup"; past it a spectral evaluation is not a candidate whatever its degree.
const SPARSE_TERM_LIMIT: usize = 16;

/// Majority's weight at or below [`CONCENTRATION_DEGREE`] must exceed parity's by
/// more than this, or the transform cannot tell a low-degree function from a
/// high-degree one and no number below it means anything.
const CALIBRATION_SEPARATION: f64 = 0.5;

// ─── timing constants ───────────────────────────────────────────────────────

/// Evaluations per timed pass.
///
/// 2²⁰. Enough that the table path — one L1 load and a shift — takes a quarter of
/// a millisecond rather than a handful of microseconds, which is what makes the
/// ratio against the spectral path a measurement instead of clock granularity.
const WORKLOAD: usize = 1 << 20;

/// Timed repeats per path. Odd, so the median is an observation rather than a
/// mean of two.
const REPEATS: usize = 7;

/// The `SplitMix64` seed for the workload's case indices. Stated so the sequence
/// is reproducible; recorded as `workload_seed` on every row.
const WORKLOAD_SEED: u64 = 0x2545_F491_4F6C_DD1D;

/// The case whose output is flipped for the sensitivity control, and the bits it
/// is flipped by.
///
/// `common::boolean::corrupt`'s own documentation measures `(37, 2)` as a
/// corruption that splits the eight influences — bits 0 and 3 of the triangle
/// count are structurally undetectable by an influence check, so the case and the
/// delta both matter and neither is arbitrary.
const CORRUPT_CASE: usize = 37;

/// See [`CORRUPT_CASE`].
const CORRUPT_DELTA: u32 = 2;

// ─── the workload ───────────────────────────────────────────────────────────

/// `SplitMix64`, ten lines, so the workload is deterministic without a
/// dependency.
struct SplitMix64(u64);

impl SplitMix64 {
    /// The next 64 bits.
    fn next_u64(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }
}

/// [`WORKLOAD`] case indices, uniform over all 256 cases.
///
/// Uniform rather than field-sampled on purpose. A real grid's case histogram is
/// dominated by the two empty cases, which would put both paths in the same
/// branch predictor and the same cache line and measure the *fixture* instead of
/// the two evaluations. C3 is a question about the arithmetic of the two forms.
fn workload() -> Vec<u8> {
    let mut rng = SplitMix64(WORKLOAD_SEED);
    (0..WORKLOAD)
        .map(|_| (rng.next_u64() & 0xFF) as u8)
        .collect()
}

// ─── the two evaluations C3 compares ────────────────────────────────────────

/// The shipped form: one load from a 256-entry table, one shift, one mask.
#[inline]
fn eval_table(values: &[u32; 256], bit: u32, x: u8) -> u8 {
    ((values[usize::from(x)] >> bit) & 1) as u8
}

/// The branchless spectral form: `f(x) = 1` iff
/// `sum_S fhat(S) chi_S(x) < 0`.
///
/// `chi_S(x) = (-1)^|S & x|`, and its sign is taken as
/// `1 - 2 * ((mask & x).count_ones() & 1)` — arithmetic on a popcount, so the
/// inner loop has no test in it and the cost is exactly `terms.len()`
/// multiply-adds. Only the non-negligible coefficients are summed, which is the
/// whole reason sparsity is the quantity that decides whether this form could
/// ever compete.
#[inline]
fn eval_spectral(terms: &[(u8, f64)], x: u8) -> u8 {
    let mut sum = 0.0f64;
    for &(mask, coeff) in terms {
        let parity = (mask & x).count_ones() & 1;
        sum += coeff * (1.0 - 2.0 * f64::from(parity));
    }
    u8::from(sum < 0.0)
}

/// One timed pass over the workload, in milliseconds.
///
/// The accumulator goes through `black_box` so the optimiser cannot notice that
/// nothing reads it and delete the loop. A deleted loop reads as zero
/// milliseconds and would make C3's ratio a division by nothing, which is why the
/// positivity of both medians is a vacuity control.
fn timed_pass(work: &[u8], mut body: impl FnMut(u8) -> u8) -> f64 {
    let started = Instant::now();
    let mut acc = 0u64;
    for &x in work {
        acc += u64::from(body(x));
    }
    black_box(acc);
    started.elapsed().as_secs_f64() * 1e3
}

// ─── the readings ───────────────────────────────────────────────────────────

/// One reading of the table's output, and the bits it is analysed in.
struct Reading {
    /// First half of `output_bit`.
    name: &'static str,
    /// One label per row — second half of `output_bit`.
    labels: Vec<String>,
    /// The per-case output integer.
    values: [u32; 256],
    /// A calibration function rather than a reading of the shipped table.
    is_control: bool,
}

impl Reading {
    /// A reading of the shipped table, in `ceil(log2(max + 1)) + 1` bits.
    ///
    /// The `+ 1` is the constant-zero witness. See the module header.
    fn shipped(name: &'static str, values: [u32; 256]) -> Self {
        let max = values.iter().copied().max().expect("256 cases");
        let bits = (32 - max.leading_zeros()).min(31) + 1;
        Self {
            name,
            labels: (0..bits).map(|b| format!("bit{b}")).collect(),
            values,
            is_control: false,
        }
    }

    /// A calibration function, read as a single `0`/`1` bit through the same
    /// path as a reading of the table — which is the point of it.
    fn calibration(label: &'static str, table: &Bool8) -> Self {
        let mut values = [0u32; 256];
        for (slot, &b) in values.iter_mut().zip(table.0.iter()) {
            *slot = u32::from(b);
        }
        Self {
            name: "calibration",
            labels: vec![String::from(label)],
            values,
            is_control: true,
        }
    }
}

// ─── the analysis, before any clock is read ─────────────────────────────────

/// Everything the transform says about one output bit.
struct Spectrum {
    /// Second half of `output_bit`.
    label: String,
    /// Which bit of the reading's output integer.
    bit: u32,
    /// The non-negligible Fourier coefficients, as `(subset mask, coefficient)`.
    terms: Vec<(u8, f64)>,
    /// Fourier weight by degree, `w[k] = sum_{|S| = k} fhat(S)^2`.
    weights: [f64; 9],
    /// `sum_k w[k]`, which Parseval says is 1.
    weight_sum: f64,
    /// The largest degree carrying a non-negligible coefficient.
    max_degree: u32,
    /// Non-zero coefficients of the `GF(2)` algebraic normal form.
    anf_terms: usize,
    /// The ANF's degree, which is **not** the Fourier degree.
    anf_degree: u32,
    /// How many of the 256 inputs the bit is set on.
    ones: usize,
    /// `Inf_i` for each of the eight corners.
    influences: [f64; 8],
    /// `sum_i Inf_i`.
    total_influence: f64,
    /// Does the spectral evaluation reproduce the truth table at all 256 inputs?
    reconstructs: bool,
}

/// Transform one bit of one reading and record what came out.
fn spectrum(values: &[u32; 256], bit: u32, label: &str) -> Spectrum {
    let table = Bool8::from_values(values, bit);
    let coefficients = table.fourier();
    let terms: Vec<(u8, f64)> = coefficients
        .iter()
        .enumerate()
        .filter(|(_, c)| c.abs() > NEGLIGIBLE)
        .map(|(s, &c)| (s as u8, c))
        .collect();

    let weights = table.weight_by_degree();
    let (anf_coefficients, anf_terms) = table.anf();
    let anf_degree = anf_coefficients
        .iter()
        .enumerate()
        .filter(|&(_, &c)| c != 0)
        .map(|(s, _)| s.count_ones())
        .max()
        .unwrap_or(0);

    let mut influences = [0.0f64; 8];
    for (i, slot) in influences.iter_mut().enumerate() {
        // `Bool8::influence` computes the spectral sum *and* the flip count and
        // asserts they agree, so reaching the next line is itself the check.
        *slot = table.influence(i);
    }

    let reconstructs = table
        .0
        .iter()
        .enumerate()
        .all(|(x, &b)| eval_spectral(&terms, x as u8) == b);

    Spectrum {
        label: String::from(label),
        bit,
        terms,
        weights,
        weight_sum: weights.iter().sum(),
        max_degree: table.max_degree(),
        anf_terms,
        anf_degree,
        ones: table.0.iter().map(|&b| usize::from(b)).sum(),
        influences,
        total_influence: table.total_influence(),
        reconstructs,
    }
}

/// The smallest degree whose cumulative weight reaches `target`.
///
/// Read cumulatively rather than by picking the heaviest degree, because "the
/// spectrum is concentrated below `k`" is a statement about the *sum* and a
/// function can be spread over three low degrees without any one of them being
/// large.
fn degree_for_weight(weights: &[f64; 9], target: f64) -> u32 {
    let mut acc = 0.0f64;
    for (k, w) in weights.iter().enumerate() {
        acc += w;
        if acc >= target {
            return k as u32;
        }
    }
    8
}

/// What the row is, derived from its own measurement rather than transcribed
/// from the module's documentation.
///
/// A constant, a parity, a single-coefficient function or an informative one.
/// Only `constant` is degenerate for C2 — a constant's spectrum cannot falsify
/// the clause either way — and in particular a parity is *maximally* informative
/// here, with zero weight below degree 8.
fn role(s: &Spectrum) -> &'static str {
    if s.max_degree == 0 {
        "constant"
    } else if s.weights[8] > 1.0 - NEGLIGIBLE {
        "parity"
    } else if s.terms.len() == 1 {
        "single_coefficient"
    } else {
        "informative"
    }
}

// ─── one measured row ───────────────────────────────────────────────────────

/// One row: the spectrum, the two clocks and the three verdicts.
struct Row {
    /// `<reading>.<label>`.
    output_bit: String,
    /// The reading this row is a bit of.
    reading: &'static str,
    /// A calibration function rather than the shipped table.
    is_control: bool,
    /// Derived by [`role`].
    role: &'static str,
    /// The transform's output.
    spectrum: Spectrum,
    /// Weight at or below [`CONCENTRATION_DEGREE`].
    concentration: f64,
    /// Cumulative weight at degrees 1, 3 and 4, so the negative is quantitative
    /// at more than one threshold.
    concentration_1: f64,
    /// See `concentration_1`.
    concentration_3: f64,
    /// See `concentration_1`.
    concentration_4: f64,
    /// Smallest degree reaching [`WEIGHT_TARGET`] of the mass.
    degree_for_target: u32,
    /// `concentration >= WEIGHT_TARGET`.
    concentrated: bool,
    /// `spectrum.terms.len() <= SPARSE_TERM_LIMIT`.
    sparse: bool,
    /// Milliseconds per [`WORKLOAD`] evaluations, median of [`REPEATS`].
    table_ms: f64,
    /// Fastest and slowest table repeat.
    table_span: (f64, f64),
    /// Milliseconds per [`WORKLOAD`] evaluations, median of [`REPEATS`].
    spectral_ms: f64,
    /// Fastest and slowest spectral repeat.
    spectral_span: (f64, f64),
    /// Does the reading's integer recompose from its bits' spectral evaluations?
    reading_reconstructs: bool,
}

impl Row {
    /// The registered antecedent of C3: C2 falsified **and** the spectrum sparse.
    fn c3_reachable(&self) -> bool {
        self.concentrated && self.sparse
    }

    /// C1: Parseval, the bit's own reconstruction, and the vector-valued
    /// recomposition.
    fn c1(&self) -> bool {
        (self.spectrum.weight_sum - 1.0).abs() <= NEGLIGIBLE
            && self.spectrum.reconstructs
            && self.reading_reconstructs
    }

    /// C2: the spectrum is **not** concentrated on low degrees.
    fn c2(&self) -> bool {
        !self.concentrated
    }

    /// C3: the branchless spectral evaluation is at least as fast as the lookup.
    fn c3(&self) -> bool {
        self.spectral_ms <= self.table_ms
    }

    /// A branchless spectral evaluation is a candidate: sparse **and** winning.
    fn branchless_feasible(&self) -> bool {
        self.sparse && self.c3()
    }

    /// A constant's spectrum cannot falsify C2 either way, so it is excluded from
    /// the run-level verdict.
    fn is_degenerate(&self) -> bool {
        self.role == "constant"
    }
}

/// Median, min and max of a set of repeats, sorted with `total_cmp`.
fn spread(mut samples: Vec<f64>) -> (f64, f64, f64) {
    samples.sort_by(f64::total_cmp);
    let last = samples.len() - 1;
    (samples[samples.len() / 2], samples[0], samples[last])
}

/// Analyse and time every bit of one reading.
fn measure(reading: &Reading, work: &[u8]) -> Vec<Row> {
    let spectra: Vec<Spectrum> = reading
        .labels
        .iter()
        .enumerate()
        .map(|(bit, label)| spectrum(&reading.values, bit as u32, label))
        .collect();

    // C1's vector-valued half: recompose the output integer from the per-bit
    // spectral evaluations. A per-bit family that does not recompose is a
    // decomposition of some other function.
    let reading_reconstructs = reading.values.iter().enumerate().all(|(x, &want)| {
        let mut got = 0u32;
        for s in &spectra {
            got |= u32::from(eval_spectral(&s.terms, x as u8)) << s.bit;
        }
        got == want
    });

    spectra
        .into_iter()
        .map(|s| {
            let mut table_samples = Vec::with_capacity(REPEATS);
            let mut spectral_samples = Vec::with_capacity(REPEATS);

            // Warm up both paths once before either is timed, then interleave
            // them so a governor step lands on both rather than on one (M-280).
            black_box(timed_pass(work, |x| eval_table(&reading.values, s.bit, x)));
            black_box(timed_pass(work, |x| eval_spectral(&s.terms, x)));
            for _ in 0..REPEATS {
                table_samples.push(timed_pass(work, |x| eval_table(&reading.values, s.bit, x)));
                spectral_samples.push(timed_pass(work, |x| eval_spectral(&s.terms, x)));
            }
            let (table_ms, table_lo, table_hi) = spread(table_samples);
            let (spectral_ms, spectral_lo, spectral_hi) = spread(spectral_samples);

            let concentration = s.weights[..=(CONCENTRATION_DEGREE as usize)].iter().sum();
            Row {
                output_bit: format!("{}.{}", reading.name, s.label),
                reading: reading.name,
                is_control: reading.is_control,
                role: role(&s),
                concentration,
                concentration_1: s.weights[..=1].iter().sum(),
                concentration_3: s.weights[..=3].iter().sum(),
                concentration_4: s.weights[..=4].iter().sum(),
                degree_for_target: degree_for_weight(&s.weights, WEIGHT_TARGET),
                concentrated: concentration >= WEIGHT_TARGET,
                sparse: s.terms.len() <= SPARSE_TERM_LIMIT,
                table_ms,
                table_span: (table_lo, table_hi),
                spectral_ms,
                spectral_span: (spectral_lo, spectral_hi),
                reading_reconstructs,
                spectrum: s,
            }
        })
        .collect()
}

// ─── formatting ─────────────────────────────────────────────────────────────

/// Nine per-degree weights as `w0|w1|…|w8`.
///
/// `|`-joined because `Run::record` refuses a value containing a comma — the
/// writer does not quote, so a comma would shift every later column.
fn joined(values: &[f64]) -> String {
    let parts: Vec<String> = values.iter().map(|v| format!("{v:.9}")).collect();
    parts.join("|")
}

// ─── the run ────────────────────────────────────────────────────────────────

fn main() {
    if !std::env::args().any(|arg| arg == "--bench") {
        return;
    }
    let prereg = isomesh::experiment!("P-167");

    common::experiment::run(prereg, |run| {
        // ── vacuity control 1: the instrument, against closed forms ─────────
        let (parity_weight8, majority_weight1, majority_influence) = self_check();
        let majority_weight1_closed = 8.0 * (70.0 / 256.0f64).powi(2);
        let majority_influence_closed = 8.0 * 35.0 / 128.0;
        assert!(
            (parity_weight8 - 1.0).abs() <= NEGLIGIBLE,
            "VOID: parity's Fourier weight at degree 8 is {parity_weight8}, not 1, so the \
             Walsh-Hadamard transform is not measuring what it claims and no number in this \
             experiment means anything"
        );
        assert!(
            (majority_weight1 - majority_weight1_closed).abs() <= NEGLIGIBLE,
            "VOID: majority-of-8's degree-1 weight is {majority_weight1}, not the closed form \
             8*(70/256)^2 = {majority_weight1_closed}, so the transform is wrong on a function \
             whose spectrum is known exactly"
        );
        assert!(
            (majority_influence - majority_influence_closed).abs() <= NEGLIGIBLE,
            "VOID: majority-of-8's total influence is {majority_influence}, not the closed form \
             8*C(7,3)/2^7 = {majority_influence_closed}, so the spectral and combinatorial \
             readings of influence do not agree"
        );

        // ── vacuity control 2: the instrument must separate the two ─────────
        let parity = Bool8::parity();
        let majority = Bool8::majority();
        let parity_low = parity.concentration_up_to(CONCENTRATION_DEGREE);
        let majority_low = majority.concentration_up_to(CONCENTRATION_DEGREE);
        assert!(
            parity_low <= NEGLIGIBLE,
            "VOID: parity carries {parity_low} of its weight at or below degree \
             {CONCENTRATION_DEGREE}, and a degree-8 character carries none, so the transform is \
             leaking weight downward and every concentration below is inflated"
        );
        assert!(
            majority_low - parity_low > CALIBRATION_SEPARATION,
            "VOID: majority's low-degree weight ({majority_low}) exceeds parity's ({parity_low}) \
             by only {} — the transform cannot tell a degree-1-heavy function from a degree-8 \
             one, so C2's verdict would be an artefact of the instrument",
            majority_low - parity_low
        );

        // ── vacuity control 3: the transform must see the table's contents ──
        let triangle_counts = shipped_triangle_counts();
        let corrupted = corrupt(&triangle_counts, CORRUPT_CASE, CORRUPT_DELTA);
        let shipped_bit1 = Bool8::from_values(&triangle_counts, 1).weight_by_degree();
        let corrupt_bit1 = Bool8::from_values(&corrupted, 1).weight_by_degree();
        let corrupt_shift: f64 = shipped_bit1
            .iter()
            .zip(corrupt_bit1.iter())
            .map(|(a, b)| (a - b).abs())
            .sum();
        assert!(
            corrupt_shift > NEGLIGIBLE,
            "VOID: flipping case {CORRUPT_CASE} by {CORRUPT_DELTA} moved the degree-weight \
             vector by {corrupt_shift}, so the transform is not reading the table it was handed \
             and every spectrum below could belong to any other table"
        );

        // ── the census ──────────────────────────────────────────────────────
        let readings = [
            Reading::shipped("triangle_counts", triangle_counts),
            Reading::shipped("edge_masks", shipped_edge_masks()),
            Reading::shipped("centroid_counts", shipped_centroid_counts()),
            Reading::calibration("parity", &parity),
            Reading::calibration("majority", &majority),
        ];
        let work = workload();
        let rows: Vec<Row> = readings
            .iter()
            .flat_map(|reading| measure(reading, &work))
            .collect();

        // ── vacuity control 4: `sparse` could have come out either way ──────
        let sparse_rows = rows.iter().filter(|r| r.sparse).count();
        let dense_rows = rows.len() - sparse_rows;
        assert!(
            sparse_rows > 0 && dense_rows > 0,
            "VOID: {sparse_rows} sparse rows and {dense_rows} dense ones — `sparse` is a constant \
             column over this census, so it is not a measurement and C3's antecedent could not \
             have been reached or missed (M-44)"
        );

        // ── vacuity control 5: the primary reading is not a constant ────────
        let primary_degrees: Vec<u32> = rows
            .iter()
            .filter(|r| r.reading == "triangle_counts")
            .map(|r| r.spectrum.max_degree)
            .collect();
        assert!(
            primary_degrees.iter().any(|&d| d > 0),
            "VOID: every bit of shipped_triangle_counts has Fourier degree 0, so C2's \
             concentration is computed against a constant and its verdict is a property of the \
             fixture rather than of the table"
        );

        // ── vacuity control 6: neither timed loop was elided ────────────────
        for r in &rows {
            assert!(
                r.table_ms > 0.0 && r.spectral_ms > 0.0,
                "VOID: {} timed at {} ms table and {} ms spectral — a zero means the optimiser \
                 deleted the loop and C3's ratio is a division by nothing",
                r.output_bit,
                r.table_ms,
                r.spectral_ms
            );
        }

        let (classes, class_count) = corner_symmetry_classes();
        println!(
            "the cube group's orbits on the eight corners: {class_count} class(es), \
             assignment {classes:?} — so R-169's 'equal within each octahedral symmetry class' \
             means all eight equal"
        );
        println!(
            "calibration: parity W^8 = {parity_weight8}, majority W^1 = {majority_weight1} \
             (closed form {majority_weight1_closed}), majority total influence = \
             {majority_influence} (closed form {majority_influence_closed})"
        );
        println!(
            "parity's weight at or below degree {CONCENTRATION_DEGREE} is {parity_low}, \
             majority's is {majority_low}; corrupting case {CORRUPT_CASE} by {CORRUPT_DELTA} \
             moved the degree-weight vector by {corrupt_shift}"
        );
        println!(
            "{sparse_rows} of {} rows are sparse at <= {SPARSE_TERM_LIMIT} non-negligible \
             coefficients, {dense_rows} are dense\n",
            rows.len()
        );
        println!(
            "{:<28} {:>4} {:>4} {:>5} {:>10} {:>4} {:>7} {:>7} {:>7}",
            "output_bit", "deg", "anf", "terms", "conc<=2", "d90", "tbl ms", "spc ms", "ratio"
        );

        for r in &rows {
            let ratio = r.spectral_ms / r.table_ms;
            println!(
                "{:<28} {:>4} {:>4} {:>5} {:>10.6} {:>4} {:>7.4} {:>7.4} {:>7.1}  {} {}",
                r.output_bit,
                r.spectrum.max_degree,
                r.spectrum.anf_terms,
                r.spectrum.terms.len(),
                r.concentration,
                r.degree_for_target,
                r.table_ms,
                r.spectral_ms,
                ratio,
                r.role,
                if r.is_control { "control" } else { "" }
            );

            run.record(&[
                ("output_bit", r.output_bit.clone()),
                ("fourier_weight_by_degree", joined(&r.spectrum.weights)),
                ("spectral_concentration", format!("{:.9}", r.concentration)),
                ("anf_terms", r.spectrum.anf_terms.to_string()),
                ("max_degree", r.spectrum.max_degree.to_string()),
                ("sparse", r.sparse.to_string()),
                ("branchless_feasible", r.branchless_feasible().to_string()),
                ("eval_ms_table", format!("{:.6}", r.table_ms)),
                ("eval_ms_spectral", format!("{:.6}", r.spectral_ms)),
                ("c1_holds", r.c1().to_string()),
                ("c2_holds", r.c2().to_string()),
                ("c3_holds", r.c3().to_string()),
                // ── extras (M-273) ──────────────────────────────────────────
                ("reading", String::from(r.reading)),
                ("bit_index", r.spectrum.label.clone()),
                ("role", String::from(r.role)),
                ("is_control", r.is_control.to_string()),
                ("is_degenerate", r.is_degenerate().to_string()),
                ("truth_table_ones", r.spectrum.ones.to_string()),
                ("spectral_terms", r.spectrum.terms.len().to_string()),
                ("anf_degree", r.spectrum.anf_degree.to_string()),
                ("weight_sum", format!("{:.12}", r.spectrum.weight_sum)),
                ("concentration_degree", CONCENTRATION_DEGREE.to_string()),
                ("weight_target", format!("{WEIGHT_TARGET:.2}")),
                ("concentration_up_to_1", format!("{:.9}", r.concentration_1)),
                ("concentration_up_to_3", format!("{:.9}", r.concentration_3)),
                ("concentration_up_to_4", format!("{:.9}", r.concentration_4)),
                ("degree_for_90pct_weight", r.degree_for_target.to_string()),
                ("concentrated", r.concentrated.to_string()),
                ("sparse_term_limit", SPARSE_TERM_LIMIT.to_string()),
                ("c3_reachable", r.c3_reachable().to_string()),
                ("eval_ms_ratio", format!("{ratio:.6}")),
                ("eval_ms_table_min", format!("{:.6}", r.table_span.0)),
                ("eval_ms_table_max", format!("{:.6}", r.table_span.1)),
                ("eval_ms_spectral_min", format!("{:.6}", r.spectral_span.0)),
                ("eval_ms_spectral_max", format!("{:.6}", r.spectral_span.1)),
                (
                    "eval_scatter_table",
                    format!("{:.6}", r.table_span.1 / r.table_span.0),
                ),
                (
                    "eval_scatter_spectral",
                    format!("{:.6}", r.spectral_span.1 / r.spectral_span.0),
                ),
                ("workload_evaluations", WORKLOAD.to_string()),
                ("timed_repeats", REPEATS.to_string()),
                ("workload_seed", format!("{WORKLOAD_SEED:#018x}")),
                (
                    "total_influence",
                    format!("{:.9}", r.spectrum.total_influence),
                ),
                ("influence_by_corner", joined(&r.spectrum.influences)),
                (
                    "influences_all_equal",
                    r.spectrum
                        .influences
                        .iter()
                        .all(|v| (v - r.spectrum.influences[0]).abs() <= NEGLIGIBLE)
                        .to_string(),
                ),
                ("corner_symmetry_classes", class_count.to_string()),
                (
                    "bit_reconstruction_exact",
                    r.spectrum.reconstructs.to_string(),
                ),
                (
                    "reading_reconstruction_exact",
                    r.reading_reconstructs.to_string(),
                ),
                (
                    "calibration_parity_weight_degree8",
                    format!("{parity_weight8:.12}"),
                ),
                (
                    "calibration_majority_weight_degree1",
                    format!("{majority_weight1:.12}"),
                ),
                (
                    "calibration_majority_total_influence",
                    format!("{majority_influence:.12}"),
                ),
                (
                    "corrupt_control_weight_shift",
                    format!("{corrupt_shift:.9}"),
                ),
                ("control_sparse_rows", sparse_rows.to_string()),
                ("control_dense_rows", dense_rows.to_string()),
            ]);
        }

        // ── the run-level verdict ───────────────────────────────────────────
        //
        // Taken over the primary reading's non-constant bits. The two
        // calibration rows are controls and the constant-zero rows cannot
        // falsify C2 either way, so both are excluded and both say so in their
        // own columns.
        let primary: Vec<&Row> = rows
            .iter()
            .filter(|r| r.reading == "triangle_counts" && !r.is_degenerate())
            .collect();
        let held = |clause: fn(&Row) -> bool| primary.iter().filter(|r| clause(r)).count();
        println!();
        println!(
            "over the {} non-constant bits of the primary reading: C1 held on {}, C2 on {}, \
             C3 on {}",
            primary.len(),
            held(Row::c1),
            held(Row::c2),
            held(Row::c3)
        );
        println!(
            "C3's registered antecedent (C2 falsified AND sparse) was reached on {} of {} rows \
             across the whole census; the branchless spectral form was feasible on {}",
            rows.iter().filter(|r| r.c3_reachable()).count(),
            rows.len(),
            rows.iter().filter(|r| r.branchless_feasible()).count()
        );
        let worst = primary
            .iter()
            .map(|r| r.degree_for_target)
            .max()
            .unwrap_or(0);
        let target_pct = WEIGHT_TARGET * 100.0;
        println!(
            "the primary reading's informative bits need degree {worst} to reach \
             {target_pct:.0}% of their Fourier mass, against degree \
             {CONCENTRATION_DEGREE} for a two-corner edge parity"
        );
    });
}
