//! **P-98 - the Plantinga-Vegter certificate, fused into the extraction gather.**
//!
//! Ticket: R-098. Pre-registered before this harness existed.
//!
//! ```bash
//! cargo bench --bench experiment_p98
//! ```
//!
//! Writes `docs/experiments/p-98.csv`.
//!
//! # SHARE, recomputed from `p-62.csv` before this harness ran
//!
//! `M-378` reported the standalone predicate at **0.2107** of extraction
//! (`thin_plate`, 65³) and derived a fused cost of **0.0658** by subtracting the
//! bare eight-corner gather. Recomputing `(predicate_ms - gather_ms) /
//! extract_ms` column by column from `docs/experiments/p-62.csv` at 65³:
//!
//! | field | predicate ms | gather ms | extract ms | standalone | fused | gather/predicate |
//! |---|---:|---:|---:|---:|---:|---:|
//! | `box_exact` | 0.7627 | 0.5116 | 3.8130 | 0.2000 | **0.065854** | 0.671 |
//! | `csg_difference` | 0.7771 | 0.5117 | 4.2717 | 0.1819 | 0.062130 | 0.658 |
//! | `thin_plate` | 0.7209 | 0.5101 | 3.4218 | **0.2107** | 0.061605 | 0.708 |
//! | `torus` | 0.7605 | 0.5120 | 4.0624 | 0.1872 | 0.061171 | 0.673 |
//! | `sphere` | 0.7474 | 0.5075 | 4.2511 | 0.1758 | 0.056432 | 0.679 |
//! | `noise_cavity` | 1.3600 | 0.5081 | 19.2885 | 0.0705 | 0.044166 | 0.374 |
//! | `gyroid` | 1.0266 | 0.5115 | 15.1811 | 0.0676 | 0.033930 | 0.498 |
//! | `fbm_terrain` | 0.8712 | 0.5133 | 51.0242 | 0.0171 | 0.007014 | 0.589 |
//!
//! **The bar reproduces, and not by the arithmetic the prose states.** 0.0658 is
//! the **maximum of the `fused_share` column**, `box_exact`'s 0.065854, and it
//! agrees with `p-62.csv`'s own `fused_share` cell (0.065844) to rounding. It is
//! *not* "21% minus two thirds": that reads `0.2107 - (2/3)(0.2107) = 0.0702`,
//! and `thin_plate`'s own fused share is `0.0616`. The 21% and the 0.0658 come
//! from **different rows** - the worst standalone share is `thin_plate`, the
//! worst fused share is `box_exact` - so the prose is a gloss over a column
//! maximum rather than a derivation. The number is right; the sentence is not a
//! recipe for it.
//!
//! **Reachability, before running anything.** The bar is an *upper bound on
//! added cost*, so what has to fit inside it is the marginal work the fused
//! predicate does per cell. `0.0658 x extract_ms` over 262,144 cells at
//! ~4.2 GHz:
//!
//! | field | budget ms | budget ns/cell | budget cycles/cell | active share |
//! |---|---:|---:|---:|---:|
//! | `thin_plate` | 0.2252 | 0.859 | **3.6** | 0.78% |
//! | `box_exact` | 0.2509 | 0.957 | 4.0 | 2.20% |
//! | `torus` | 0.2673 | 1.020 | 4.3 | 1.61% |
//! | `sphere` | 0.2797 | 1.067 | 4.5 | 1.82% |
//! | `csg_difference` | 0.2811 | 1.072 | 4.5 | 2.29% |
//! | `gyroid` | 0.9989 | 3.810 | 16.0 | 8.18% |
//! | `noise_cavity` | 1.2692 | 4.842 | 20.3 | 10.82% |
//! | `fbm_terrain` | 3.3574 | 12.808 | 53.8 | 3.21% |
//!
//! The tightest field is `thin_plate` at **3.6 cycles per cell**, and the fused
//! predicate's per-cell work on an inactive cell is *one branch on a byte the
//! extractor already built*. Clause two only runs on the active share. So the
//! bar is **reachable with roughly 3x headroom on its own tightest field**, and
//! it is reachable for a reason that is a fact about the predicate rather than
//! about this machine: clause one is `0 in-not box-F(C)`, which for a trilinear
//! is exactly *"all eight corners share a sign"*, which is exactly
//! `case == 0 || case == 255`. **Fusing does not make clause one cheaper; it
//! replaces an eight-value sign scan with one comparison on a byte the extractor
//! already holds - which is a term `M-378`'s subtraction left in the predicate's
//! account, and the measured reason the derived 0.0658 is an over-estimate.**
//!
//! # A share clause needs both halves in the same window, and that is a finding
//!
//! The first instrument timed the numerator and the denominator in separate
//! loops. On this host - twelve cores shared with sibling agents' benches, load
//! average between 3 and 16 for the whole session - four runs of the *identical
//! binary* put the worst 65³ `fused_share` at **0.0457, 0.0471, 0.0527 and
//! 0.0779**: a 1.7x spread straddling the 0.0658 bar, three HELD and one
//! FALSIFIED. Nothing about the predicate changed between those runs.
//!
//! The fix is to time the extraction **inside the same repetition** as the walks
//! and take the minimum of both series ([`Paired::extract`]), so a quiet
//! numerator can never be divided by a loud denominator. Five runs after that
//! change read 0.0475, 0.0484, 0.0501, 0.0511, 0.0529 - a 1.1x spread, all HELD.
//! **`M-281`'s "one build and one run" has to be read as one *window*, not one
//! process, whenever the clause is a ratio.**
//!
//! The other half of the answer is not to depend on a clock at all. `M-337`'s
//! re-audit converted a registered 1.25x timing floor into exact integer
//! equalities; the same move here is [`Counted`], which reports that the fused
//! predicate adds **zero** grid loads to the extractor's `8 x cells` while a
//! standalone pass adds a second `8 x cells`, and that clause two runs on
//! exactly `active_cells`. Those integers are identical on every run.
//!
//! # What "fused" means here, and why it is bench-local
//!
//! `crates/isomesh/src/**` is read-only for this ticket, so the fused walk lives
//! in this file: [`gather_case`] is the extractor's own per-cell prologue
//! (`marching_cubes/mod.rs:259-268`, eight loads and eight sign tests building
//! the case byte), and [`fused_certified`] takes that byte and those eight
//! registers. Nothing is re-read. The marginal cost is measured by difference
//! against [`walk_case`], which is the identical walk with the predicate
//! removed, and **not** by subtracting a separately-timed gather, which is the
//! step `M-378`'s derivation had to take and the one this harness exists to
//! check.
//!
//! # The three clauses
//!
//! - **C1** `fused_predicate_ms / extract_ms <= 0.0658` at 65³ on all eight
//!   fields, where `fused_predicate_ms` is `min(walk_fused) - min(walk_case)`
//!   over 25 short regions timed inside the same repetitions, and every other
//!   time on the row is the same minimum estimator. See [`REPS`] for why a
//!   minimum and not a mean.
//! - **C2** the certified **sets** are identical, `set_difference == 0`, and no
//!   certificate lands on a cell `A-020` calls a tunnel or a twelve-vertex
//!   contour. Sets, not sizes: equal counts over different sets is exactly the
//!   failure this clause exists to catch, so the comparison is a symmetric
//!   difference over a bitset indexed by global cell number.
//! - **C3** the certified fraction is aggregable per chunk with
//!   `extra_passes == 0`.
//!
//! # Controls
//!
//! - **The registered vacuity control: the random arm, inherited verbatim.** The
//!   eight reference fields at 17³-65³ give **seven** tunnel cells in 172,032,
//!   which `M-378` itself called *"a hair from `M-44`'s vacuous zero"*. So
//!   `p-62.csv`'s 400,000-cell arm comes across with the same LCG, the same
//!   seed and the same draw order, and the harness **asserts it reproduces
//!   2,202 tunnel and 180 twelve-vertex cells**. A `unsound_certificates == 0`
//!   over that population is a soundness result; over seven cells it is a
//!   silence.
//! - **The set comparison must be able to report a difference.** A zero
//!   symmetric difference between two implementations of the same predicate is
//!   worthless unless the comparison can see a difference at all. So a third
//!   set is built from [`mutant_certified`] - the canonical transcription error
//!   in a three-axis block, the `z` differences left paired on the `y` corners -
//!   and `mutant_set_difference` is asserted non-zero over the sweep.
//! - **The timer must be able to resolve clause two.** `fused_predicate_ms` is a
//!   difference of two timings, so a small value could mean "free" or
//!   "invisible". [`walk_forced`] runs clause two on *every* cell instead of the
//!   active few, and `forced_margin_ms` is asserted non-zero.
//! - **The rig must be able to resolve a value at the bar, and this is the
//!   control C1's one-sidedness turns on.** An upper bound is passed by
//!   "unresolvably small", so the harness measures its own floor - two timings
//!   of the *identical* baseline loop inside the same repetition - and asserts
//!   `resolution_floor_ms / extract_ms < 0.0658` on every 65³ row. Any row whose
//!   marginal falls at or below that floor is printed as *unresolvably small*
//!   rather than as a cost, and its `fused_predicate_ms` is reported signed:
//!   `thin_plate` read **-0.0180 ms** before the pairing went in, which is not a
//!   negative cost, it is drift.
//! - **A fused marginal must not exceed the forced margin.** The fused walk runs
//!   clause two on the active cells, the forced walk on all of them; the
//!   ordering is arithmetic, and violating it is how this harness's first
//!   baseline was found to be measuring LLVM rather than the predicate.
//! - **Three independent implementations must agree on the count.** The flat
//!   fused walk, the chunk-major fused walk and the shipped
//!   `validate::isotopy_report` are cross-checked against each other on every
//!   row.

#![allow(
    clippy::cast_possible_truncation,
    clippy::cast_precision_loss,
    clippy::cast_sign_loss,
    clippy::too_many_lines
)]

mod common;

use std::time::Instant;

use isomesh::extractor::Extractor;
use isomesh::fields::ReferenceField;
use isomesh::for_each_reference_field;
use isomesh::marching_cubes::MarchingCubes;
use isomesh::marching_cubes::ambiguity::joined_mask;
use isomesh::marching_cubes::table::AMBIGUOUS_FACES;
use isomesh::marching_cubes::trilinear::{BodySaddles, Contours, Topology};
use isomesh::validate::{cell_is_certified, isotopy_report};
use isomesh::{MeshBuffer, RuntimeShape3, Sdf};

/// Samples per axis. `P-62`'s three, because C2 is stated over its rows.
const RESOLUTIONS: [u32; 3] = [17, 33, 65];

/// The registered bar. `M-378`'s derived figure, used as a ceiling rather than a
/// target.
const BAR: f64 = 0.0658;

/// `cube.rs`'s corner numbering: corner `i` sits at `(i&1, (i>>1)&1, (i>>2)&1)`.
const CORNERS: [[usize; 3]; 8] = [
    [0, 0, 0],
    [1, 0, 0],
    [0, 1, 0],
    [1, 1, 0],
    [0, 0, 1],
    [1, 0, 1],
    [0, 1, 1],
    [1, 1, 1],
];

/// Cells per chunk edge for C3's aggregation. `P-72`/`M-377` measured 8³ as the
/// granularity optimum, and 8 divides 16, 32 and 64 exactly, so every arm
/// partitions its cells without a remainder.
const CHUNK_EDGE: usize = 8;

/// Cells drawn for the vacuity control. `P-62`'s arm, verbatim.
const RANDOM_CELLS: u64 = 400_000;

/// Chunks the random arm is partitioned into for C3. 400,000 / 500 = 800 cells
/// each, exactly.
const RANDOM_CHUNKS: u64 = 500;

/// `M-378`'s counts on the random arm. Asserted, because "inherited verbatim" is
/// a claim about the fixture and not a note about intent.
const RANDOM_TUNNELS: u64 = 2_202;
/// `M-378`'s twelve-vertex-contour count on the random arm.
const RANDOM_TWELVE: u64 = 180;

/// Timed regions per measurement.
///
/// **The estimator is the minimum of the series, not the mean or the median, and
/// that is a decision about this host rather than a habit.** Contention and a
/// `powersave` governor can only ever *add* time to a region, so the minimum of
/// many short regions is the closest reading to what the code costs; a mean
/// carries every sibling agent's bench in it. Twenty-five regions because
/// `fused_predicate_ms` is a **difference** of two ~1 ms numbers - the shape
/// `M-337`'s re-audit found swinging 30% when taken once - and because this
/// machine's load average was between 4 and 16 for the whole of this run. The
/// median of the paired differences is reported beside it as
/// `fused_predicate_median_ms`; the two disagreeing is information, not an error.
const REPS: usize = 25;

/// Cells per timed region, so that a 17³ arm and a 65³ arm time comparable work
/// and the timer's own resolution is never the measurement.
///
/// **Short on purpose.** A long region is more likely to be preempted, and a
/// minimum estimator needs at least one region that ran clean: 300,000 cells is
/// one sweep at 65³ (about 1.1 ms), nine at 33³ and seventy-three at 17³.
const TIMED_CELLS: u64 = 300_000;

/// `cube::is_inside`, which is private. Negative is inside; exactly zero is
/// outside.
#[inline(always)]
fn is_inside(v: f64) -> bool {
    v < 0.0
}

/// Sweeps per timed region for a grid of this many cells.
fn sweeps_for(cells: u64) -> u32 {
    (TIMED_CELLS / cells).max(1) as u32
}

/// The clock, on the row. `M-280`: on a governed CPU a nanosecond is not a unit
/// unless the reader can see what the core was doing.
///
/// `NaN` where `/proc/cpuinfo` is not readable - the same honesty as
/// `common::experiment`'s `unknown` provenance, and visible in the CSV.
fn cpu_mhz() -> f64 {
    std::fs::read_to_string("/proc/cpuinfo")
        .ok()
        .and_then(|s| {
            s.lines()
                .find(|l| l.starts_with("cpu MHz"))
                .and_then(|l| l.split(':').nth(1))
                .and_then(|v| v.trim().parse::<f64>().ok())
        })
        .unwrap_or(f64::NAN)
}

/// Minimum milliseconds **per sweep** over [`REPS`] timed regions.
///
/// See [`REPS`] for why the minimum: contention adds time and never removes it.
fn min_ms(sweeps: u32, mut region: impl FnMut()) -> f64 {
    let mut best = f64::INFINITY;
    for _ in 0..REPS {
        let start = Instant::now();
        for _ in 0..sweeps {
            region();
        }
        best = best.min(start.elapsed().as_nanos() as f64 / 1e6 / f64::from(sweeps));
    }
    best
}

/// Median of a collected series.
fn median_of(mut samples: Vec<f64>) -> f64 {
    samples.sort_unstable_by(f64::total_cmp);
    samples[samples.len() / 2]
}

/// Minimum of a collected series.
fn min_of(samples: &[f64]) -> f64 {
    samples.iter().copied().fold(f64::INFINITY, f64::min)
}

// ── the predicate, twice ─────────────────────────────────────────────────────

/// Lower bound of the interval `[lo, hi] * [lo, hi]`, over the four
/// corner differences of one axis.
///
/// The min/max is the exact range of that partial derivative over the cell,
/// because a trilinear's `dF/dx` is bilinear in `(y, z)` and therefore a convex
/// combination of the four `x`-edge differences. The self-product's minimum is
/// `lo²` or `hi²` when the interval keeps a sign and `lo*hi <= 0` when it
/// straddles - and that straddle is the whole test.
#[inline(always)]
fn axis_low(d: [f64; 4]) -> f64 {
    let mut lo = d[0];
    let mut hi = d[0];
    for v in &d[1..] {
        if *v < lo {
            lo = *v;
        }
        if *v > hi {
            hi = *v;
        }
    }
    if lo > 0.0 {
        lo * lo
    } else if hi < 0.0 {
        hi * hi
    } else {
        lo * hi
    }
}

/// Plantinga-Vegter's second clause, `<box-grad-F, box-grad-F> > 0`, from eight
/// corner values already in registers.
///
/// The pairings are `isotopy.rs`'s `X_PAIRS`/`Y_PAIRS`/`Z_PAIRS` written out: an
/// `x` difference pairs corners differing only in bit 0, and so on. `h²` factors
/// out of a sum of three squares, which is why this runs on raw differences.
#[inline(always)]
fn clause_two(c: &[f64; 8]) -> bool {
    let x = axis_low([c[1] - c[0], c[3] - c[2], c[5] - c[4], c[7] - c[6]]);
    let y = axis_low([c[2] - c[0], c[3] - c[1], c[6] - c[4], c[7] - c[5]]);
    let z = axis_low([c[4] - c[0], c[5] - c[1], c[6] - c[2], c[7] - c[3]]);
    x + y + z > 0.0
}

/// The fused certificate: the case byte the extractor already built, plus the
/// eight corners it already loaded.
///
/// **Clause one is not made cheaper by fusion, it is made to vanish.**
/// `0 in-not box-F(C)` is exactly "all eight corners share a sign" for a
/// trilinear, which is exactly `case == 0 || case == 255` - and `case` is the
/// index the extractor built to reach its triangulation table. The standalone
/// predicate re-derives the same fact with a fresh eight-value sign scan
/// (`isotopy.rs:129-132`).
#[inline(always)]
fn fused_certified(case: u8, c: &[f64; 8]) -> bool {
    case == 0 || case == 255 || clause_two(c)
}

/// The fused certificate with the canonical transcription error: the `z` block
/// copy-pasted from `y` and its corner indices left alone.
///
/// This exists so that `set_difference == 0` is a measurement. A symmetric
/// difference that cannot come back non-zero is not a comparison, and the
/// arithmetic-drift this clause is registered against would look exactly like
/// this.
#[inline(always)]
fn mutant_certified(case: u8, c: &[f64; 8]) -> bool {
    if case == 0 || case == 255 {
        return true;
    }
    let x = axis_low([c[1] - c[0], c[3] - c[2], c[5] - c[4], c[7] - c[6]]);
    let y = axis_low([c[2] - c[0], c[3] - c[1], c[6] - c[4], c[7] - c[5]]);
    let z = axis_low([c[2] - c[0], c[3] - c[1], c[6] - c[4], c[7] - c[5]]);
    x + y + z > 0.0
}

/// A certificate that says yes to everything.
///
/// Not a straw man: it is what `M-44`'s rule requires be counted. A run whose
/// `unsound_certificates` is zero and whose `unsound_if_certify_all` is also
/// zero has measured nothing, because the population the counters were gated on
/// was empty. This is the number the gate would have read for the worst possible
/// predicate, and it is counted through the same branch.
#[inline(always)]
const fn always_certified(_case: u8, _c: &[f64; 8]) -> bool {
    true
}

// ── the walks ────────────────────────────────────────────────────────────────

/// The extractor's own per-cell prologue: eight loads, eight sign tests, one
/// case byte (`marching_cubes/mod.rs:259-268`).
#[inline(always)]
fn gather_case(values: &[f64], n: usize, x: usize, y: usize, z: usize) -> ([f64; 8], u8) {
    let base = (z * n + y) * n + x;
    let mut c = [0.0f64; 8];
    let mut case = 0u8;
    for (i, slot) in c.iter_mut().enumerate() {
        let o = CORNERS[i];
        let v = values[base + (o[2] * n + o[1]) * n + o[0]];
        *slot = v;
        if is_inside(v) {
            case |= 1 << i;
        }
    }
    (c, case)
}

/// The same gather with no case byte, which is what a standalone predicate pass
/// has to do.
#[inline(always)]
fn gather(values: &[f64], n: usize, x: usize, y: usize, z: usize) -> [f64; 8] {
    let base = (z * n + y) * n + x;
    let mut c = [0.0f64; 8];
    for (i, slot) in c.iter_mut().enumerate() {
        let o = CORNERS[i];
        *slot = values[base + (o[2] * n + o[1]) * n + o[0]];
    }
    c
}

/// **The first baseline, kept because it is the fixture defect this harness's
/// own control caught.** Gather, build the case byte, discard the values.
///
/// The eight corner values are dead after their sign bits are read, so LLVM
/// drops the loads-into-registers entirely and leaves a sign scan. **The
/// extractor cannot do that**: `corner_value` is live afterwards for
/// `joined_mask` and for every edge interpolation (`mod.rs:285`, `mod.rs:644`).
/// Subtracting this baseline charges the predicate for *materialising the
/// corners*, which is the one cost the registration says fusion removes - and it
/// produced a `fused_predicate_ms` of 0.6920 ms against a `forced_margin_ms` of
/// 0.8328 ms on `sphere` at 65³, i.e. 83% of the cost of clause two on **every**
/// cell for running it on 1.8% of them. That ordering is impossible, which is
/// how the defect was found rather than reported.
fn walk_signs(values: &[f64], n: usize, cells: usize) -> u64 {
    let mut acc = 0u64;
    for z in 0..cells {
        for y in 0..cells {
            for x in 0..cells {
                let (_, case) = gather_case(values, n, x, y, z);
                acc = acc.wrapping_add(u64::from(case));
            }
        }
    }
    acc
}

/// **The baseline.** The extractor's prologue: gather eight corners, build the
/// case byte, and leave the corners live.
///
/// The `black_box` on the array is what makes this the extractor's loop rather
/// than a sign scan, and it is identical in [`walk_fused`] and [`walk_forced`],
/// so it cancels out of every difference taken against it. Its residual bias is
/// conservative: clause two then reads the corners back from a stack slot
/// instead of holding them in registers, which charges the fused predicate a few
/// extra loads on the active cells rather than fewer.
fn walk_case(values: &[f64], n: usize, cells: usize) -> u64 {
    let mut acc = 0u64;
    for z in 0..cells {
        for y in 0..cells {
            for x in 0..cells {
                let (c, case) = gather_case(values, n, x, y, z);
                std::hint::black_box(&c);
                acc = acc.wrapping_add(u64::from(case));
            }
        }
    }
    acc
}

/// **The fused walk.** Identical to [`walk_case`] plus the certificate off the
/// same registers. `walk_fused_ms - walk_case_ms` is the whole of C1.
fn walk_fused(values: &[f64], n: usize, cells: usize) -> (u64, u64) {
    let mut acc = 0u64;
    let mut certified = 0u64;
    for z in 0..cells {
        for y in 0..cells {
            for x in 0..cells {
                let (c, case) = gather_case(values, n, x, y, z);
                std::hint::black_box(&c);
                acc = acc.wrapping_add(u64::from(case));
                if fused_certified(case, &c) {
                    certified += 1;
                }
            }
        }
    }
    (acc, certified)
}

/// Clause one only: the certificate's free half, in isolation.
///
/// **This is the ticket's central mechanistic claim, isolated.** `0 in-not
/// box-F(C)` for a trilinear is exactly `case == 0 || case == 255`, and the case
/// byte is already in a register because the extractor built it to index its
/// table. So `walk_clause_one - walk_case` is what fusing clause one costs, and
/// the prediction is *nothing measurable*. The standalone predicate pays for the
/// same fact with a fresh eight-value sign scan (`isotopy.rs:129-132`).
fn walk_clause_one(values: &[f64], n: usize, cells: usize) -> (u64, u64) {
    let mut acc = 0u64;
    let mut certified = 0u64;
    for z in 0..cells {
        for y in 0..cells {
            for x in 0..cells {
                let (c, case) = gather_case(values, n, x, y, z);
                std::hint::black_box(&c);
                acc = acc.wrapping_add(u64::from(case));
                if case == 0 || case == 255 {
                    certified += 1;
                }
            }
        }
    }
    (acc, certified)
}

/// Clause two on **every** cell rather than the active few.
///
/// The timing control. `fused_predicate_ms` is a difference and could be small
/// because the predicate is nearly free or because the rig cannot see it; this
/// arm decides which, by doing 100% of the clause-two arithmetic instead of the
/// 1-11% an active-cell short circuit leaves.
fn walk_forced(values: &[f64], n: usize, cells: usize) -> (u64, u64) {
    let mut acc = 0u64;
    let mut certified = 0u64;
    for z in 0..cells {
        for y in 0..cells {
            for x in 0..cells {
                let (c, case) = gather_case(values, n, x, y, z);
                std::hint::black_box(&c);
                acc = acc.wrapping_add(u64::from(case));
                if clause_two(&c) {
                    certified += 1;
                }
            }
        }
    }
    (acc, certified)
}

/// The standalone predicate, `P-62`'s loop: gather into an array, call the
/// shipped `validate::cell_is_certified`.
fn walk_standalone(values: &[f64], n: usize, cells: usize) -> u64 {
    let mut certified = 0u64;
    for z in 0..cells {
        for y in 0..cells {
            for x in 0..cells {
                if cell_is_certified(&gather(values, n, x, y, z)) {
                    certified += 1;
                }
            }
        }
    }
    certified
}

/// The bare eight-corner gather, `P-62`'s `gather_ms`: the term `M-378`'s
/// derivation subtracted.
fn walk_gather(values: &[f64], n: usize, cells: usize) {
    for z in 0..cells {
        for y in 0..cells {
            for x in 0..cells {
                std::hint::black_box(gather(values, n, x, y, z));
            }
        }
    }
}

/// One repetition's milliseconds-per-sweep for the walks C1's difference is
/// taken over.
struct Paired {
    signs: f64,
    case: f64,
    case_again: f64,
    clause_one: f64,
    fused: f64,
    forced: f64,
    /// C1's **denominator**, timed in the same repetition as its numerator.
    ///
    /// `M-281` says compare within one build and one run; taken literally on a
    /// machine whose load average moved between 4 and 16 during this run, it
    /// means within one *window*. A numerator from a quiet minute over a
    /// denominator from a loud one is not a share of anything.
    extract: f64,
}

/// Time the six walks **inside one repetition**, [`REPS`] times.
///
/// Paired rather than sequential, because `fused_predicate_ms` is the difference
/// of two numbers near 1 ms on a `powersave` governor: timing all of one walk
/// and then all of the other puts every clock excursion between them straight
/// into the difference, and the first version of this harness read
/// `thin_plate`'s marginal as **-0.0180 ms** that way.
///
/// `case_again` is the **resolution floor**: two timings of the *identical* loop
/// inside the same window, whose spread is the smallest difference this rig can
/// honestly claim to have seen. For an upper-bound clause that floor is the
/// whole question - a cost below it is "at most the bar" only if the floor
/// itself is below the bar, and the harness asserts that it is.
fn paired_walks(
    values: &[f64],
    n: usize,
    cells: usize,
    sweeps: u32,
    extract: &mut dyn FnMut(),
) -> Vec<Paired> {
    let mut out = Vec::with_capacity(REPS);
    for _ in 0..REPS {
        let one = |f: &mut dyn FnMut()| -> f64 {
            let start = Instant::now();
            for _ in 0..sweeps {
                f();
            }
            start.elapsed().as_nanos() as f64 / 1e6 / f64::from(sweeps)
        };
        // `black_box` on the slice, not just on the result: the walks are pure
        // functions of it, so an opaque pointer is what stops LLVM computing one
        // sweep and reusing it for the rest.
        let signs = one(&mut || {
            std::hint::black_box(walk_signs(std::hint::black_box(values), n, cells));
        });
        let case = one(&mut || {
            std::hint::black_box(walk_case(std::hint::black_box(values), n, cells));
        });
        let fused = one(&mut || {
            std::hint::black_box(walk_fused(std::hint::black_box(values), n, cells));
        });
        let clause_one = one(&mut || {
            std::hint::black_box(walk_clause_one(std::hint::black_box(values), n, cells));
        });
        let forced = one(&mut || {
            std::hint::black_box(walk_forced(std::hint::black_box(values), n, cells));
        });
        let case_again = one(&mut || {
            std::hint::black_box(walk_case(std::hint::black_box(values), n, cells));
        });
        // One extraction, not `sweeps` of them: an extraction is 3-60 ms and the
        // sweep count exists to make a 1 ms walk measurable.
        let start = Instant::now();
        extract();
        let extract_ms = start.elapsed().as_nanos() as f64 / 1e6;
        out.push(Paired {
            signs,
            case,
            case_again,
            clause_one,
            fused,
            forced,
            extract: extract_ms,
        });
    }
    out
}

/// Grid loads and clause-two evaluations, **counted rather than timed**.
///
/// This is the machine-independent half of C1, and after four runs of this
/// binary it is the only half that reproduces. `M-337`'s re-audit is the
/// precedent: a registered 1.25x floor that re-measured at 1.022 three runs
/// later, re-registered as exact integer equalities and held on 15 of 15 rows.
/// Here the integers are:
///
/// - the extractor's prologue reads `8 x cells` grid values, and the **fused
///   predicate adds none of them**;
/// - a standalone pass adds exactly `8 x cells` more, which is the whole of
///   `M-378`'s *"two thirds of the standalone cost is re-reading corners"*
///   expressed as a count instead of a clock;
/// - clause two is evaluated on exactly `active_cells`, so the per-field spread
///   in the fused cost is an exact integer and not a timing artefact.
struct Counted {
    /// Grid values read.
    loads: u64,
    /// Cells on which clause two's arithmetic ran.
    clause_two_evals: u64,
}

/// [`gather_case`] with the loads counted.
#[inline(always)]
fn counted_gather_case(
    values: &[f64],
    n: usize,
    x: usize,
    y: usize,
    z: usize,
    counts: &mut Counted,
) -> ([f64; 8], u8) {
    let base = (z * n + y) * n + x;
    let mut c = [0.0f64; 8];
    let mut case = 0u8;
    for (i, slot) in c.iter_mut().enumerate() {
        let o = CORNERS[i];
        let v = values[base + (o[2] * n + o[1]) * n + o[0]];
        counts.loads += 1;
        *slot = v;
        if is_inside(v) {
            case |= 1 << i;
        }
    }
    (c, case)
}

/// The extractor's prologue, counted. No predicate.
fn counted_case(values: &[f64], n: usize, cells: usize) -> Counted {
    let mut counts = Counted {
        loads: 0,
        clause_two_evals: 0,
    };
    for z in 0..cells {
        for y in 0..cells {
            for x in 0..cells {
                let _ = counted_gather_case(values, n, x, y, z, &mut counts);
            }
        }
    }
    counts
}

/// The fused walk, counted. Same loads as [`counted_case`], plus clause two on
/// the cells the case byte says are active.
fn counted_fused(values: &[f64], n: usize, cells: usize) -> Counted {
    let mut counts = Counted {
        loads: 0,
        clause_two_evals: 0,
    };
    for z in 0..cells {
        for y in 0..cells {
            for x in 0..cells {
                let (c, case) = counted_gather_case(values, n, x, y, z, &mut counts);
                if case == 0 || case == 255 {
                    continue;
                }
                counts.clause_two_evals += 1;
                std::hint::black_box(clause_two(&c));
            }
        }
    }
    counts
}

/// The standalone pass, counted: its own gather, on top of the extractor's.
fn counted_standalone(values: &[f64], n: usize, cells: usize) -> Counted {
    let mut counts = Counted {
        loads: 0,
        clause_two_evals: 0,
    };
    for z in 0..cells {
        for y in 0..cells {
            for x in 0..cells {
                let base = (z * n + y) * n + x;
                let mut c = [0.0f64; 8];
                for (i, slot) in c.iter_mut().enumerate() {
                    let o = CORNERS[i];
                    *slot = values[base + (o[2] * n + o[1]) * n + o[0]];
                    counts.loads += 1;
                }
                std::hint::black_box(cell_is_certified(&c));
            }
        }
    }
    counts
}

/// What a chunk-major fused walk can report without a second traversal.
struct ChunkAgg {
    chunks: u64,
    /// Sum over chunks of that chunk's active-cell count.
    active_sum: u64,
    /// Sum over chunks of that chunk's certified **active** cells.
    certified_active_sum: u64,
    /// Certified cells including inactive ones, for cross-checking the flat walk.
    certified_all: u64,
    /// Cells touched. One per cell is `extra_passes == 0`; anything more is C3's
    /// falsifier.
    cells_visited: u64,
    /// Mean over chunks of that chunk's own certified fraction, computed exactly
    /// as `IsotopyReport::certified_fraction` does - **1.0 for a chunk with no
    /// active cells**. Reported because it is what a consumer who stores a
    /// fraction per chunk would read, and it is not the global fraction.
    mean_fraction: f64,
}

/// **C3's walk.** Chunk-major, one visit per cell, per-chunk counters.
fn walk_chunked(values: &[f64], n: usize, cells: usize, edge: usize) -> ChunkAgg {
    let per_axis = cells / edge;
    let mut agg = ChunkAgg {
        chunks: (per_axis as u64).pow(3),
        active_sum: 0,
        certified_active_sum: 0,
        certified_all: 0,
        cells_visited: 0,
        mean_fraction: 0.0,
    };
    let mut fraction_sum = 0.0f64;
    for cz in 0..per_axis {
        for cy in 0..per_axis {
            for cx in 0..per_axis {
                let mut active = 0u64;
                let mut certified = 0u64;
                for z in cz * edge..(cz + 1) * edge {
                    for y in cy * edge..(cy + 1) * edge {
                        for x in cx * edge..(cx + 1) * edge {
                            let (c, case) = gather_case(values, n, x, y, z);
                            std::hint::black_box(&c);
                            agg.cells_visited += 1;
                            let cert = fused_certified(case, &c);
                            if cert {
                                agg.certified_all += 1;
                            }
                            if case != 0 && case != 255 {
                                active += 1;
                                if cert {
                                    certified += 1;
                                }
                            }
                        }
                    }
                }
                agg.active_sum += active;
                agg.certified_active_sum += certified;
                fraction_sum += if active == 0 {
                    1.0
                } else {
                    certified as f64 / active as f64
                };
            }
        }
    }
    agg.mean_fraction = fraction_sum / agg.chunks as f64;
    agg
}

// ── the certified sets ───────────────────────────────────────────────────────

/// A set of cells, one bit per global cell index.
///
/// C2 is registered over the certified **set**, not its size: two predicates
/// that certify the same number of different cells have the same count and
/// nothing else in common, and that is precisely the drift the clause is
/// watching for.
struct CellSet(Vec<u64>);

impl CellSet {
    fn new(cells: u64) -> Self {
        Self(vec![0u64; (cells as usize).div_ceil(64)])
    }

    #[inline(always)]
    fn insert(&mut self, i: u64) {
        self.0[(i / 64) as usize] |= 1u64 << (i % 64);
    }

    fn count(&self) -> u64 {
        self.0.iter().map(|w| u64::from(w.count_ones())).sum()
    }

    /// Cells in exactly one of the two sets.
    fn symmetric_difference(&self, other: &Self) -> u64 {
        self.0
            .iter()
            .zip(&other.0)
            .map(|(a, b)| u64::from((a ^ b).count_ones()))
            .sum()
    }
}

/// One field at one resolution, or the random arm.
struct Row {
    field: &'static str,
    resolution: u32,
    cells: u64,
    active_cells: u64,
    refused_cells: u64,
    certified_fused: u64,
    certified_standalone: u64,
    certified_active_fused: u64,
    set_difference: u64,
    mutant_set_difference: u64,
    unsound_fused: u64,
    unsound_standalone: u64,
    unsound_mutant: u64,
    /// What a predicate that certified **everything** would score on this row's
    /// hidden-topology population, counted through the same `if hidden { if
    /// certified { .. } }` branch the real counters use.
    ///
    /// `M-44`'s rule in one column: a zero in `unsound_certificates` is only a
    /// soundness result if the counter could have been non-zero, and this is the
    /// number it would have read.
    unsound_if_certify_all: u64,
    tunnel_cells: u64,
    twelve_vertex_cells: u64,
    /// The first, defective baseline: gather, case byte, corners discarded.
    /// Reported so the report can name the 30x it moves the verdict by.
    walk_signs_ms: f64,
    walk_case_ms: f64,
    walk_fused_ms: f64,
    walk_forced_ms: f64,
    /// **C1's numerator**: `min(fused) - min(case)` over [`REPS`] short regions
    /// timed inside the same repetitions. See [`REPS`] for the estimator.
    fused_predicate_ms: f64,
    /// The same difference read as the **median of the paired per-repetition
    /// differences** instead.
    ///
    /// Reported next to the registered column because this host is running other
    /// agents' benches on the same twelve cores, and two estimators disagreeing
    /// is information rather than an error.
    fused_predicate_median_ms: f64,
    /// **The ticket's central mechanistic claim, isolated.** What clause one
    /// costs when it is read off the case byte instead of re-scanning eight
    /// values. Predicted: nothing measurable.
    clause_one_margin_ms: f64,
    /// The same numerator against the defective sign-scan baseline. Reported,
    /// never used for a verdict.
    fused_vs_signs_ms: f64,
    /// The timing control: clause two on every cell instead of the active few.
    /// A `fused_predicate_ms` at or above this is impossible and means the
    /// baseline is not the extractor's prologue.
    forced_margin_ms: f64,
    /// The rig's own resolution: two timings of the identical baseline loop in
    /// the same window.
    resolution_floor_ms: f64,
    chunked_fused_ms: f64,
    standalone_predicate_ms: f64,
    bare_gather_ms: f64,
    /// The shipped route to the same fraction: `validate::isotopy_report` over
    /// the values array. A **separate traversal**, which is the contrast C3 is
    /// about.
    isotopy_report_ms: f64,
    extract_ms: f64,
    chunks: u64,
    chunk_cells: u64,
    chunk_active_sum: u64,
    chunk_certified_sum: u64,
    cells_visited: u64,
    mean_chunk_fraction: f64,
    global_fraction: f64,
    /// Grid values the extractor's own prologue reads. Counted, not derived.
    loads_extraction: u64,
    /// Grid values the **fused** predicate adds on top of that. The registered
    /// mechanism says zero, and this is the column that says so in integers.
    loads_fused_additional: u64,
    /// Grid values a **standalone** pass adds. `M-378`'s "two thirds is
    /// re-reading corners", as a count.
    loads_standalone_additional: u64,
    /// Cells on which clause two's arithmetic ran in the fused walk.
    clause_two_evals: u64,
    sweeps: u32,
    cpu_mhz: f64,
}

impl Row {
    fn fused_share(&self) -> f64 {
        self.fused_predicate_ms / self.extract_ms
    }

    /// The floor as a share of extraction. C1 is an upper bound, so a HELD is
    /// only a measurement when this is below the bar.
    fn floor_share(&self) -> f64 {
        self.resolution_floor_ms / self.extract_ms
    }

    /// The registered share, pushed up by the rig's own resolution.
    ///
    /// C1 is one-sided, and six of eight fields at 65³ produce a marginal inside
    /// the noise band - four of them **negative**, which is a code-layout bias
    /// rather than a negative cost. Adding the floor is how an upper bound is
    /// stated honestly when the quantity is smaller than the instrument: if this
    /// is under the bar, no reading the rig could have produced would have
    /// exceeded it.
    fn fused_share_upper(&self) -> f64 {
        (self.fused_predicate_ms.abs() + self.resolution_floor_ms) / self.extract_ms
    }

    /// **The mechanism, as an additive model whose every term is measured where
    /// that term is large.**
    ///
    /// `clause_one_margin_ms` is the branch on the case byte, over all cells.
    /// `forced_margin_ms / cells` is clause two's per-cell cost measured over
    /// 262,144 evaluations - a 40-50% effect on the walk, far above the noise
    /// that swallows the fused marginal itself - and it is charged only on the
    /// active cells, which is an exact integer count. This is the estimate a
    /// reader should trust on the six fields whose direct marginal sits inside
    /// the floor, and it is corroborated by the two where it does not.
    fn modelled_marginal_ms(&self) -> f64 {
        if self.cells == 0 {
            return f64::NAN;
        }
        self.clause_one_margin_ms
            + self.forced_margin_ms * (self.active_cells as f64 / self.cells as f64)
    }

    fn modelled_share(&self) -> f64 {
        self.modelled_marginal_ms() / self.extract_ms
    }

    fn standalone_share(&self) -> f64 {
        self.standalone_predicate_ms / self.extract_ms
    }

    /// Traversals beyond the first. C3's falsifier is any value above zero.
    fn extra_passes(&self) -> u64 {
        self.cells_visited / self.cells - 1
    }
}

/// The same LCG `P-62` drew its arm with, at the same seed, in the same order.
struct Lcg(u64);

impl Lcg {
    /// A value in `[-1, 1)`.
    fn signed(&mut self) -> f64 {
        self.0 = self
            .0
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        f64::from((self.0 >> 40) as u32) / f64::from(1u32 << 23) - 1.0
    }
}

/// **The registered vacuity control.** 400,000 random cells, inherited verbatim
/// from `p-62.csv`.
///
/// Untimed: there is no grid here and no extraction, so a wall time would invite
/// comparison with the sampled rows. What this arm supplies is the *population* -
/// 2,202 tunnels and 180 twelve-vertex contours - without which
/// `unsound_certificates == 0` is a statement about seven cells.
fn random_arm() -> Row {
    let mut rng = Lcg(0x2026_u64 ^ 0x5EED_1234);
    let mut fused_set = CellSet::new(RANDOM_CELLS);
    let mut standalone_set = CellSet::new(RANDOM_CELLS);
    let mut mutant_set = CellSet::new(RANDOM_CELLS);

    let mut row = blank_row("random_cells", 0, RANDOM_CELLS);
    row.chunks = RANDOM_CHUNKS;
    row.chunk_cells = RANDOM_CELLS / RANDOM_CHUNKS;
    row.cells_visited = RANDOM_CELLS;
    row.cpu_mhz = cpu_mhz();

    let mut fraction_sum = 0.0f64;
    let mut chunk_active = 0u64;
    let mut chunk_certified = 0u64;

    for i in 0..RANDOM_CELLS {
        let mut c = [0.0f64; 8];
        let mut case = 0u8;
        for (k, slot) in c.iter_mut().enumerate() {
            let v = rng.signed();
            *slot = v;
            if is_inside(v) {
                case |= 1 << k;
            }
        }
        let active = case != 0 && case != 255;

        let fused = fused_certified(case, &c);
        let standalone = cell_is_certified(&c);
        let mutant = mutant_certified(case, &c);
        if fused {
            fused_set.insert(i);
            row.certified_fused += 1;
        }
        if standalone {
            standalone_set.insert(i);
            row.certified_standalone += 1;
        }
        if mutant {
            mutant_set.insert(i);
        }
        if active {
            row.active_cells += 1;
            if fused {
                row.certified_active_fused += 1;
                chunk_certified += 1;
            } else {
                row.refused_cells += 1;
            }
            chunk_active += 1;
        }

        let mask = joined_mask(&c, AMBIGUOUS_FACES[case as usize]);
        let saddles = BodySaddles::of(&c);
        let hidden = match Contours::of(case, mask).topology(&saddles) {
            Topology::Tunnel => {
                row.tunnel_cells += 1;
                true
            }
            Topology::TwelveVertexContour => {
                row.twelve_vertex_cells += 1;
                true
            }
            _ => false,
        };
        if hidden {
            if fused {
                row.unsound_fused += 1;
            }
            if standalone {
                row.unsound_standalone += 1;
            }
            if mutant {
                row.unsound_mutant += 1;
            }
            // The same branch with the certificate replaced by `true`: what the
            // counter reads for a predicate that says yes to everything, which
            // is the only thing that makes the zeros above measurements.
            if always_certified(case, &c) {
                row.unsound_if_certify_all += 1;
            }
        }

        // C3, on the arm that has the most interesting fraction to aggregate.
        if (i + 1) % row.chunk_cells == 0 {
            fraction_sum += if chunk_active == 0 {
                1.0
            } else {
                chunk_certified as f64 / chunk_active as f64
            };
            chunk_active = 0;
            chunk_certified = 0;
        }
    }

    row.chunk_active_sum = row.active_cells;
    row.chunk_certified_sum = row.certified_active_fused;
    row.mean_chunk_fraction = fraction_sum / RANDOM_CHUNKS as f64;
    row.global_fraction = row.certified_active_fused as f64 / row.active_cells as f64;
    row.set_difference = fused_set.symmetric_difference(&standalone_set);
    row.mutant_set_difference = fused_set.symmetric_difference(&mutant_set);
    assert_eq!(
        fused_set.count(),
        row.certified_fused,
        "random arm: the fused certified set and its counter disagree"
    );
    assert_eq!(
        standalone_set.count(),
        row.certified_standalone,
        "random arm: the standalone certified set and its counter disagree"
    );
    row
}

/// A row with everything zeroed, so no arm can forget a field and land a silent
/// default in the CSV.
fn blank_row(field: &'static str, resolution: u32, cells: u64) -> Row {
    Row {
        field,
        resolution,
        cells,
        active_cells: 0,
        refused_cells: 0,
        certified_fused: 0,
        certified_standalone: 0,
        certified_active_fused: 0,
        set_difference: 0,
        mutant_set_difference: 0,
        unsound_fused: 0,
        unsound_standalone: 0,
        unsound_mutant: 0,
        unsound_if_certify_all: 0,
        tunnel_cells: 0,
        twelve_vertex_cells: 0,
        walk_signs_ms: 0.0,
        walk_case_ms: 0.0,
        walk_fused_ms: 0.0,
        walk_forced_ms: 0.0,
        fused_predicate_ms: 0.0,
        fused_predicate_median_ms: 0.0,
        clause_one_margin_ms: 0.0,
        fused_vs_signs_ms: 0.0,
        forced_margin_ms: 0.0,
        resolution_floor_ms: 0.0,
        chunked_fused_ms: 0.0,
        standalone_predicate_ms: 0.0,
        bare_gather_ms: 0.0,
        isotopy_report_ms: 0.0,
        extract_ms: 0.0,
        chunks: 0,
        chunk_cells: 0,
        chunk_active_sum: 0,
        chunk_certified_sum: 0,
        cells_visited: 0,
        mean_chunk_fraction: f64::NAN,
        global_fraction: f64::NAN,
        loads_extraction: 0,
        loads_fused_additional: 0,
        loads_standalone_additional: 0,
        clause_two_evals: 0,
        sweeps: 0,
        cpu_mhz: f64::NAN,
    }
}

fn measure<F: Sdf<Scalar = f64> + ReferenceField>(
    field: &F,
    field_name: &'static str,
    samples: u32,
) -> Row {
    let shape = RuntimeShape3::new([samples; 3]).expect("a cubic runtime shape");
    let ([lo, hi], _) = ([field.domain().0, field.domain().1], ());
    let h = (hi[0] - lo[0]) / f64::from(samples - 1);

    // The sample grid, once. Every walk and the shipped report read these same
    // values, which is what makes the set comparison about the cell rather than
    // about two samplings.
    let n = samples as usize;
    let mut values = Vec::with_capacity(n * n * n);
    for z in 0..samples {
        for y in 0..samples {
            for x in 0..samples {
                values.push(field.sample([
                    lo[0] + h * f64::from(x),
                    lo[1] + h * f64::from(y),
                    lo[2] + h * f64::from(z),
                ]));
            }
        }
    }

    let cells = n - 1;
    let cell_count = (cells as u64).pow(3);
    let mut row = blank_row(field_name, samples, cell_count);
    row.sweeps = sweeps_for(cell_count);
    assert_eq!(
        cells % CHUNK_EDGE,
        0,
        "{field_name} at {samples}³: {cells} cells per axis is not a whole number of {CHUNK_EDGE}-cell chunks, so the partition would not be exact"
    );

    // ── the timed arms, one build and one run (M-281) ─────────────────────────
    let sweeps = row.sweeps;
    let mut mc = MarchingCubes::<f64>::new();
    let mut out = MeshBuffer::<f64>::new();
    let paired = paired_walks(&values, n, cells, sweeps, &mut || {
        let _ = mc.extract_into(field, &shape, lo, h, &mut out);
        std::hint::black_box(&out);
    });
    let series = |f: fn(&Paired) -> f64| -> Vec<f64> { paired.iter().map(f).collect() };
    let signs = series(|p| p.signs);
    let case = series(|p| p.case);
    let case_again = series(|p| p.case_again);
    let clause_one = series(|p| p.clause_one);
    let fused = series(|p| p.fused);
    let forced = series(|p| p.forced);
    row.walk_signs_ms = min_of(&signs);
    row.walk_case_ms = min_of(&case);
    row.walk_fused_ms = min_of(&fused);
    row.walk_forced_ms = min_of(&forced);
    row.fused_predicate_ms = min_of(&fused) - min_of(&case);
    // C1's denominator, from the same windows as its numerator, on the same
    // minimum estimator: the smallest extraction the run saw, which is the
    // largest share this run can support.
    row.extract_ms = min_of(&series(|p| p.extract));
    row.clause_one_margin_ms = min_of(&clause_one) - min_of(&case);
    row.fused_vs_signs_ms = min_of(&fused) - min_of(&signs);
    row.forced_margin_ms = min_of(&forced) - min_of(&case);
    // **The floor, on the same estimator as the marginal.** Two independent
    // minimum-of-25 readings of the *identical* baseline loop, taken in the same
    // window. Their gap is the smallest difference this rig can claim.
    row.resolution_floor_ms = (min_of(&case_again) - min_of(&case)).abs();
    // The other estimator, reported so a reader can see whether the two agree.
    row.fused_predicate_median_ms =
        median_of(paired.iter().map(|p| p.fused - p.case).collect());
    row.chunked_fused_ms = min_ms(sweeps, || {
        let agg = walk_chunked(
            std::hint::black_box(values.as_slice()),
            n,
            cells,
            CHUNK_EDGE,
        );
        std::hint::black_box(agg.certified_all);
    });
    row.standalone_predicate_ms = min_ms(sweeps, || {
        std::hint::black_box(walk_standalone(
            std::hint::black_box(values.as_slice()),
            n,
            cells,
        ));
    });
    row.bare_gather_ms = min_ms(sweeps, || {
        walk_gather(std::hint::black_box(values.as_slice()), n, cells);
    });
    row.isotopy_report_ms = min_ms(sweeps, || {
        let report = isotopy_report(std::hint::black_box(values.as_slice()), &shape)
            .expect("a report over a grid this harness already sized");
        std::hint::black_box(report.certified);
    });
    row.cpu_mhz = cpu_mhz();

    // ── the counted instrument, machine-independent ──────────────────────────
    let counted_case_only = counted_case(&values, n, cells);
    let counted_fused_walk = counted_fused(&values, n, cells);
    let counted_standalone_walk = counted_standalone(&values, n, cells);
    row.loads_extraction = counted_case_only.loads;
    row.loads_fused_additional = counted_fused_walk.loads - counted_case_only.loads;
    row.loads_standalone_additional = counted_standalone_walk.loads;
    row.clause_two_evals = counted_fused_walk.clause_two_evals;

    // ── the counts and the sets, untimed ─────────────────────────────────────
    let (_, fused_count) = walk_fused(&values, n, cells);
    row.certified_fused = fused_count;
    row.certified_standalone = walk_standalone(&values, n, cells);

    let agg = walk_chunked(&values, n, cells, CHUNK_EDGE);
    row.chunks = agg.chunks;
    row.chunk_cells = (CHUNK_EDGE as u64).pow(3);
    row.chunk_active_sum = agg.active_sum;
    row.chunk_certified_sum = agg.certified_active_sum;
    row.cells_visited = agg.cells_visited;
    row.mean_chunk_fraction = agg.mean_fraction;
    assert_eq!(
        agg.certified_all, fused_count,
        "{field_name} at {samples}³: the chunk-major fused walk certified {} cells and the flat \
         one {fused_count}, so the two disagree about the same predicate on the same values",
        agg.certified_all
    );

    // The shipped report, as a third independent implementation of the same
    // count. It is a **separate traversal**, which is exactly C3's contrast.
    let report = isotopy_report(&values, &shape).expect("a report over this grid");
    assert_eq!(
        report.cells, cell_count,
        "{field_name} at {samples}³: the shipped report walked {} cells, this harness {cell_count}",
        report.cells
    );
    assert_eq!(
        report.certified, agg.certified_active_sum,
        "{field_name} at {samples}³: the shipped isotopy_report certified {} active cells and the \
         fused walk {}, so fusing changed the arithmetic",
        report.certified, agg.certified_active_sum
    );
    row.active_cells = report.active_cells;
    row.certified_active_fused = agg.certified_active_sum;
    row.refused_cells = report.uncertified;
    row.global_fraction = report.certified_fraction();

    // ── the certified sets, and the cross-tabulation ─────────────────────────
    let mut fused_set = CellSet::new(cell_count);
    let mut standalone_set = CellSet::new(cell_count);
    let mut mutant_set = CellSet::new(cell_count);
    let mut index = 0u64;
    for z in 0..cells {
        for y in 0..cells {
            for x in 0..cells {
                let (c, case) = gather_case(&values, n, x, y, z);
                let fused = fused_certified(case, &c);
                let standalone = cell_is_certified(&c);
                let mutant = mutant_certified(case, &c);
                if fused {
                    fused_set.insert(index);
                }
                if standalone {
                    standalone_set.insert(index);
                }
                if mutant {
                    mutant_set.insert(index);
                }

                let mask = joined_mask(&c, AMBIGUOUS_FACES[case as usize]);
                let saddles = BodySaddles::of(&c);
                let hidden = match Contours::of(case, mask).topology(&saddles) {
                    Topology::Tunnel => {
                        row.tunnel_cells += 1;
                        true
                    }
                    Topology::TwelveVertexContour => {
                        row.twelve_vertex_cells += 1;
                        true
                    }
                    _ => false,
                };
                if hidden {
                    if fused {
                        row.unsound_fused += 1;
                    }
                    if standalone {
                        row.unsound_standalone += 1;
                    }
                    if mutant {
                        row.unsound_mutant += 1;
                    }
                    if always_certified(case, &c) {
                        row.unsound_if_certify_all += 1;
                    }
                }
                index += 1;
            }
        }
    }
    row.set_difference = fused_set.symmetric_difference(&standalone_set);
    row.mutant_set_difference = fused_set.symmetric_difference(&mutant_set);
    assert_eq!(
        fused_set.count(),
        row.certified_fused,
        "{field_name} at {samples}³: the fused set and its counter disagree, so the set \
         comparison is not over the population the counts describe"
    );
    assert_eq!(
        standalone_set.count(),
        row.certified_standalone,
        "{field_name} at {samples}³: the standalone set and its counter disagree"
    );

    row
}

type CsvRow = Vec<(&'static str, String)>;

fn main() {
    if !std::env::args().any(|arg| arg == "--bench") {
        return;
    }

    let prereg = isomesh::experiment!("P-98");
    let mut rows: Vec<Row> = Vec::new();

    // Through `for_each_reference_field!`, not a retyped list of eight: the
    // macro is the crate's own definition of "the eight reference fields", and
    // C2 is registered over `p-62.csv`'s rows, which came from the same macro.
    for samples in RESOLUTIONS {
        for_each_reference_field!(f64, |name, field| {
            rows.push(measure(&field, name, samples));
        });
    }
    rows.push(random_arm());

    println!(
        "\n{:>15} {:>4} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>9} {:>7} {:>7} {:>8}",
        "field",
        "n",
        "signs ms",
        "case ms",
        "fused ms",
        "PRED ms",
        "floor",
        "forced",
        "extract",
        "FUSEDSHR",
        "setdiff",
        "mutant",
        "UNSOUND"
    );
    for r in &rows {
        println!(
            "{:>15} {:>4} {:>8.4} {:>8.4} {:>8.4} {:>8.4} {:>8.4} {:>8.4} {:>8.4} {:>9.5} {:>7} \
             {:>7} {:>8}",
            r.field,
            r.resolution,
            r.walk_signs_ms,
            r.walk_case_ms,
            r.walk_fused_ms,
            r.fused_predicate_ms,
            r.resolution_floor_ms,
            r.forced_margin_ms,
            r.extract_ms,
            r.fused_share(),
            r.set_difference,
            r.mutant_set_difference,
            r.unsound_fused + r.unsound_standalone
        );
    }

    // ── controls ─────────────────────────────────────────────────────────────
    //
    // The registered vacuity control, verbatim: the random arm's population is
    // what makes a zero in `unsound_certificates` a soundness result rather than
    // `M-44`'s silence.
    let random = rows
        .iter()
        .find(|r| r.field == "random_cells")
        .expect("the random arm to be in the sweep");
    assert_eq!(
        (random.tunnel_cells, random.twelve_vertex_cells),
        (RANDOM_TUNNELS, RANDOM_TWELVE),
        "VACUITY CONTROL: the random arm was supposed to be inherited verbatim from p-62.csv \
         (2,202 tunnel and 180 twelve-vertex cells) and produced {} and {} instead, so this is \
         not M-378's population and C2's zero is a statement about a different fixture",
        random.tunnel_cells,
        random.twelve_vertex_cells
    );
    let hidden_population: u64 = rows
        .iter()
        .map(|r| r.tunnel_cells + r.twelve_vertex_cells)
        .sum();
    assert!(
        hidden_population > 0,
        "VOID: no tunnel and no twelve-vertex contour anywhere in {} rows",
        rows.len()
    );
    // The set comparison must be able to report a difference. Global rather than
    // per-row for `P-62`'s reason: a field whose every cell is certified under
    // both the real and the mutant arithmetic is a correct outcome, and what has
    // to be shown is that the *instrument* can see a difference.
    let mutant_total: u64 = rows.iter().map(|r| r.mutant_set_difference).sum();
    assert!(
        mutant_total > 0,
        "VOID: the mis-paired-z mutant certified exactly the same set as the fused predicate on \
         every one of {} rows, so `set_difference == 0` is not a comparison and C2 has no \
         instrument",
        rows.len()
    );
    // The unsound counters must be able to increment. `unsound_mutant` came back
    // zero - the mis-paired-z mutant certifies a *different* set without
    // certifying a hidden-topology cell - so the proof that the gate fires is
    // the degenerate certificate, counted through the same branch.
    let certify_all: u64 = rows.iter().map(|r| r.unsound_if_certify_all).sum();
    assert!(
        certify_all > 0,
        "VOID: a predicate that certifies everything would score zero unsound certificates over \
         this fixture, so `unsound_certificates == 0` is M-44's vacuous zero"
    );
    // The timer must be able to resolve clause two.
    let forced_margins: Vec<f64> = rows
        .iter()
        .filter(|r| r.resolution > 0)
        .map(|r| r.forced_margin_ms)
        .collect();
    assert!(
        forced_margins.iter().all(|m| *m > 0.0),
        "VOID: forcing clause two onto every cell did not cost measurable time on some arm, so a \
         small fused_predicate_ms cannot be distinguished from a blind timer: {forced_margins:?}"
    );
    // **The resolution-floor gate, and it is what makes C1's HELD a
    // measurement.** C1 is a one-sided upper bound, so "too small to see" passes
    // it - which is only honest if the rig could have *seen* a cost at the bar.
    // The floor is two timings of the identical baseline loop in the same window;
    // it has to sit below `BAR x extract_ms`, and by a margin, or a HELD is
    // `M-44`'s vacuous zero wearing a stopwatch.
    for r in rows.iter().filter(|r| r.resolution == 65) {
        assert!(
            r.floor_share() < BAR,
            "VOID: {} at 65³ has a timing resolution floor of {:.4} ms ({:.5} of extraction), \
             which is not below the bar {BAR}. A fused_predicate_ms under the bar on this row \
             would be indistinguishable from the rig's own noise.",
            r.field,
            r.resolution_floor_ms,
            r.floor_share()
        );
    }
    // **The baseline-ordering control, and the one that found the first
    // fixture's defect.** The fused walk runs clause two on the active cells
    // only; the forced walk runs it on every cell. So the fused marginal must be
    // *below* the forced margin on every arm. The first version of `walk_case`
    // let LLVM discard the corner values after reading their sign bits, and the
    // subtraction then charged the predicate for materialising them: 0.6920 ms
    // of "predicate" against a 0.8328 ms cost for 55x the clause-two work on
    // `sphere` at 65³. Impossible, and invisible without this ordering.
    for r in rows.iter().filter(|r| r.resolution > 0) {
        assert!(
            r.fused_predicate_ms < r.forced_margin_ms,
            "{} at {}³: the fused predicate ({:.4} ms) cost more than clause two on EVERY cell \
             ({:.4} ms) while running on {} of {} cells. The baseline is not the extractor's \
             prologue -- it is being optimised in a way the extractor cannot be.",
            r.field,
            r.resolution,
            r.fused_predicate_ms,
            r.forced_margin_ms,
            r.active_cells,
            r.cells
        );
    }
    // Every arm must have visited each of its cells exactly once, or C3's
    // `extra_passes` is measuring the wrong loop.
    for r in &rows {
        assert_eq!(
            r.cells_visited, r.cells,
            "{} at {}³: the chunk-major walk visited {} of {} cells, so extra_passes is not a \
             count of traversals",
            r.field, r.resolution, r.cells_visited, r.cells
        );
    }
    // **The counted instrument, and it is the half that reproduces.** Four runs
    // of this binary put the worst 65³ `fused_share` at 0.0457, 0.0471, 0.0527
    // and 0.0779 - a 1.7x spread straddling the 0.0658 bar - so the wall-clock
    // reading of C1 is not decidable on this host. These integers are, and they
    // are the mechanism C1 was derived from.
    for r in rows.iter().filter(|r| r.resolution > 0) {
        assert_eq!(
            r.loads_extraction,
            8 * r.cells,
            "{} at {}³: the extractor's prologue read {} grid values, not 8 x {}",
            r.field,
            r.resolution,
            r.loads_extraction,
            r.cells
        );
        assert_eq!(
            r.loads_fused_additional, 0,
            "{} at {}³: the fused predicate added {} grid loads. Fusing means reading the corners \
             the extractor already read, so anything but zero here is not a fused implementation.",
            r.field, r.resolution, r.loads_fused_additional
        );
        assert_eq!(
            r.loads_standalone_additional,
            8 * r.cells,
            "{} at {}³: a standalone pass added {} grid loads rather than a second 8 x {}, so \
             M-378's 'two thirds is re-reading corners' does not describe this fixture",
            r.field,
            r.resolution,
            r.loads_standalone_additional,
            r.cells
        );
        assert_eq!(
            r.clause_two_evals, r.active_cells,
            "{} at {}³: clause two ran on {} cells while isotopy_report counted {} active. The \
             short circuit is not the case byte.",
            r.field, r.resolution, r.clause_two_evals, r.active_cells
        );
    }

    // ── verdicts ─────────────────────────────────────────────────────────────
    let at65: Vec<&Row> = rows.iter().filter(|r| r.resolution == 65).collect();
    assert_eq!(at65.len(), 8, "C1 is registered over all eight fields at 65³");
    let worst = at65
        .iter()
        .max_by(|a, b| a.fused_share().total_cmp(&b.fused_share()))
        .expect("eight rows at 65³");
    let c1 = at65.iter().all(|r| r.fused_share() <= BAR);
    // The same clause read against the conservative upper bound. Strictly harder
    // than the registration asks, and reported because four of the eight 65³
    // marginals are negative and a negative number passing an upper bound is not
    // an argument.
    let worst_upper = at65
        .iter()
        .max_by(|a, b| a.fused_share_upper().total_cmp(&b.fused_share_upper()))
        .expect("eight rows at 65³");
    let c1_upper = at65.iter().all(|r| r.fused_share_upper() <= BAR);

    let set_difference: u64 = rows.iter().map(|r| r.set_difference).sum();
    let unsound: u64 = rows
        .iter()
        .map(|r| r.unsound_fused + r.unsound_standalone)
        .sum();
    let counts_agree = rows
        .iter()
        .all(|r| r.certified_fused == r.certified_standalone);
    let c2 = set_difference == 0 && unsound == 0 && counts_agree;

    let extra_passes: u64 = rows.iter().map(Row::extra_passes).sum();
    let aggregates = rows.iter().all(|r| {
        r.chunk_active_sum == r.active_cells && r.chunk_certified_sum == r.certified_active_fused
    });
    let c3 = extra_passes == 0 && aggregates;

    println!(
        "\nC1 worst fused share at 65³: {:.5} on {} against the registered bar {BAR} -> {}",
        worst.fused_share(),
        worst.field,
        if c1 { "HELD" } else { "FALSIFIED" }
    );
    println!(
        "   standalone share on the same row {:.5}, bare gather {:.4} ms, fused predicate {:.4} \
         ms, extract {:.4} ms, timing floor {:.4} ms ({:.5} of extraction)",
        worst.standalone_share(),
        worst.bare_gather_ms,
        worst.fused_predicate_ms,
        worst.extract_ms,
        worst.resolution_floor_ms,
        worst.floor_share()
    );
    for r in rows.iter().filter(|r| r.resolution == 65) {
        if r.fused_predicate_ms <= r.resolution_floor_ms {
            println!(
                "   {} at 65³: fused predicate {:.4} ms is AT OR BELOW the rig's {:.4} ms floor \
                 -- read as 'unresolvably small', not as a cost; the bar on this row is {:.4} ms, \
                 which is {:.1}x the floor",
                r.field,
                r.fused_predicate_ms,
                r.resolution_floor_ms,
                BAR * r.extract_ms,
                BAR * r.extract_ms / r.resolution_floor_ms
            );
        }
    }
    println!(
        "   C1 against the conservative upper bound |marginal| + floor: worst {:.5} on {} -> {}. \
         Strictly harder than the registration and the same verdict.",
        worst_upper.fused_share_upper(),
        worst_upper.field,
        if c1_upper { "HELD" } else { "FALSIFIED" }
    );
    println!(
        "   the mechanism, as an additive model measured where each term is large: modelled = \
         clause_one_margin + (clause two ns/cell) x active_cells. Clause two's per-cell cost is \
         taken from the forced arm, a 40-50% effect on the walk."
    );
    for r in rows.iter().filter(|r| r.resolution == 65) {
        println!(
            "     {:>15}: active {:>6}/{:>6} = {:>6.3}%, clause1 {:>7.4} ms, clause2 {:>5.2} \
             ns/cell, MODELLED {:>7.4} ms ({:>8.5}), measured {:>7.4} ms, median-est {:>7.4} ms",
            r.field,
            r.active_cells,
            r.cells,
            100.0 * r.active_cells as f64 / r.cells as f64,
            r.clause_one_margin_ms,
            1e6 * r.forced_margin_ms / r.cells as f64,
            r.modelled_marginal_ms(),
            r.modelled_share(),
            r.fused_predicate_ms,
            r.fused_predicate_median_ms
        );
    }
    let worst_modelled = rows
        .iter()
        .filter(|r| r.resolution == 65)
        .max_by(|a, b| a.modelled_share().total_cmp(&b.modelled_share()))
        .expect("eight rows at 65³");
    println!(
        "   worst MODELLED share at 65³: {:.5} on {} -- a bracket, NOT a second reading of C1. \
         Both of its terms are themselves unpaired timing differences, and it runs {:.2}x the \
         measured marginal on the low-activity fields and {:.2}x on noise_cavity, so it locates \
         the mechanism and does not decide the clause.",
        worst_modelled.modelled_share(),
        worst_modelled.field,
        worst_modelled.modelled_marginal_ms() / worst_modelled.fused_predicate_ms,
        rows.iter()
            .find(|r| r.field == "noise_cavity" && r.resolution == 65)
            .map_or(f64::NAN, |r| r.modelled_marginal_ms() / r.fused_predicate_ms)
    );
    println!(
        "   THE COUNTED DECOMPOSITION, which is the half that reproduces across runs: at 65³ the \
         extractor's prologue reads {} grid values, the FUSED predicate adds {}, and a STANDALONE \
         pass adds {} -- a second copy of the whole gather. Clause two runs on exactly \
         active_cells, which is an exact integer per field, not a clock.",
        worst.loads_extraction, worst.loads_fused_additional, worst.loads_standalone_additional
    );
    println!(
        "C2 set_difference {set_difference} over {} rows, unsound {unsound} over \
         {hidden_population} tunnel/twelve-vertex cells, counts agree {counts_agree} -> {}",
        rows.len(),
        if c2 { "HELD" } else { "FALSIFIED" }
    );
    println!(
        "   the comparison is an instrument: the mis-paired-z mutant differs on {mutant_total} \
         cells and certifies {} of the hidden-topology population",
        rows.iter().map(|r| r.unsound_mutant).sum::<u64>()
    );
    println!(
        "C3 extra_passes {extra_passes}, per-chunk counts sum to the global counts {aggregates} \
         -> {}",
        if c3 { "HELD" } else { "FALSIFIED" }
    );
    for r in &rows {
        if (r.mean_chunk_fraction - r.global_fraction).abs() > 1e-12 {
            println!(
                "   {} at {}³: mean of per-chunk fractions {:.6} vs global {:.6} -- the \
                 aggregable object is the COUNT PAIR, not the fraction",
                r.field, r.resolution, r.mean_chunk_fraction, r.global_fraction
            );
        }
    }
    println!(
        "\nThe shipped route to the same number is a second traversal: isotopy_report costs \
         {:.4} ms at 65³ on {} against a fused marginal of {:.4} ms, i.e. {:.1}x.",
        worst.isotopy_report_ms,
        worst.field,
        worst.fused_predicate_ms,
        worst.isotopy_report_ms / worst.fused_predicate_ms
    );
    println!(
        "The verdict is sensitive to the baseline, and the report has to say so: against the \
         sign-scan baseline that discards the corners, the same fused walk on {} at 65³ reads \
         {:.4} ms ({:.5} of extraction) instead of {:.4} ms ({:.5}). The extractor keeps \
         corner_value live for joined_mask and for every edge interpolation, so the sign scan is \
         not its prologue and that number is not the predicate's cost.",
        worst.field,
        worst.fused_vs_signs_ms,
        worst.fused_vs_signs_ms / worst.extract_ms,
        worst.fused_predicate_ms,
        worst.fused_share()
    );
    println!(
        "The registered caveat, restated because it bounds what a certificate means: a certified \
         cell's patch is a GRAPH over a coordinate plane, not necessarily ONE component."
    );

    common::experiment::run(prereg, |run| {
        for r in &rows {
            let mut csv: CsvRow = vec![
                ("field", r.field.to_string()),
                ("resolution", r.resolution.to_string()),
                ("fused_predicate_ms", format!("{:.4}", r.fused_predicate_ms)),
                (
                    "standalone_predicate_ms",
                    format!("{:.4}", r.standalone_predicate_ms),
                ),
                ("bare_gather_ms", format!("{:.4}", r.bare_gather_ms)),
                ("extract_ms", format!("{:.4}", r.extract_ms)),
                ("fused_share", format!("{:.6}", r.fused_share())),
                ("fused_share_upper", format!("{:.6}", r.fused_share_upper())),
                (
                    "modelled_marginal_ms",
                    format!("{:.4}", r.modelled_marginal_ms()),
                ),
                ("modelled_share", format!("{:.6}", r.modelled_share())),
                (
                    "clause_one_margin_ms",
                    format!("{:.4}", r.clause_one_margin_ms),
                ),
                (
                    "clause_two_ns_per_cell",
                    if r.cells == 0 {
                        format!("{:.4}", f64::NAN)
                    } else {
                        format!("{:.4}", 1e6 * r.forced_margin_ms / r.cells as f64)
                    },
                ),
                ("loads_extraction", r.loads_extraction.to_string()),
                (
                    "loads_fused_additional",
                    r.loads_fused_additional.to_string(),
                ),
                (
                    "loads_standalone_additional",
                    r.loads_standalone_additional.to_string(),
                ),
                ("clause_two_evals", r.clause_two_evals.to_string()),
                (
                    "fused_predicate_median_ms",
                    format!("{:.4}", r.fused_predicate_median_ms),
                ),
                ("c1_holds_upper_bound", c1_upper.to_string()),
                ("standalone_share", format!("{:.6}", r.standalone_share())),
                ("certified_cells_fused", r.certified_fused.to_string()),
                (
                    "certified_cells_standalone",
                    r.certified_standalone.to_string(),
                ),
                ("set_difference", r.set_difference.to_string()),
                (
                    "unsound_certificates",
                    (r.unsound_fused + r.unsound_standalone).to_string(),
                ),
                ("tunnel_cells", r.tunnel_cells.to_string()),
                ("twelve_vertex_cells", r.twelve_vertex_cells.to_string()),
                ("extra_passes", r.extra_passes().to_string()),
                ("c1_holds", c1.to_string()),
                ("c2_holds", c2.to_string()),
                ("c3_holds", c3.to_string()),
                // Extras. The decomposition C1 is a difference of, the controls
                // that make its zeros measurements, and the C3 contrast.
                ("cells", r.cells.to_string()),
                ("active_cells", r.active_cells.to_string()),
                ("refused_cells", r.refused_cells.to_string()),
                ("walk_signs_ms", format!("{:.4}", r.walk_signs_ms)),
                ("fused_vs_signs_ms", format!("{:.4}", r.fused_vs_signs_ms)),
                (
                    "fused_share_vs_signs",
                    format!("{:.6}", r.fused_vs_signs_ms / r.extract_ms),
                ),
                ("walk_case_ms", format!("{:.4}", r.walk_case_ms)),
                ("walk_fused_ms", format!("{:.4}", r.walk_fused_ms)),
                ("walk_forced_ms", format!("{:.4}", r.walk_forced_ms)),
                ("forced_margin_ms", format!("{:.4}", r.forced_margin_ms)),
                (
                    "resolution_floor_ms",
                    format!("{:.4}", r.resolution_floor_ms),
                ),
                ("floor_share", format!("{:.6}", r.floor_share())),
                ("bar_budget_ms", format!("{:.4}", BAR * r.extract_ms)),
                ("chunked_fused_ms", format!("{:.4}", r.chunked_fused_ms)),
                ("isotopy_report_ms", format!("{:.4}", r.isotopy_report_ms)),
                ("extra_passes_shipped_api", 1.to_string()),
                (
                    "certified_active_fused",
                    r.certified_active_fused.to_string(),
                ),
                ("unsound_fused", r.unsound_fused.to_string()),
                ("unsound_standalone", r.unsound_standalone.to_string()),
                ("unsound_mutant", r.unsound_mutant.to_string()),
                (
                    "unsound_if_certify_all",
                    r.unsound_if_certify_all.to_string(),
                ),
                (
                    "mutant_set_difference",
                    r.mutant_set_difference.to_string(),
                ),
                ("hidden_population", hidden_population.to_string()),
                ("chunks", r.chunks.to_string()),
                ("chunk_cells", r.chunk_cells.to_string()),
                ("chunk_active_sum", r.chunk_active_sum.to_string()),
                ("chunk_certified_sum", r.chunk_certified_sum.to_string()),
                ("cells_visited", r.cells_visited.to_string()),
                (
                    "mean_chunk_fraction",
                    format!("{:.6}", r.mean_chunk_fraction),
                ),
                ("global_fraction", format!("{:.6}", r.global_fraction)),
                ("bar", format!("{BAR:.4}")),
                ("sweeps", r.sweeps.to_string()),
                ("reps", REPS.to_string()),
                ("cpu_mhz", format!("{:.1}", r.cpu_mhz)),
            ];
            csv.push(("chunk_edge_cells", CHUNK_EDGE.to_string()));
            run.record(&csv);
        }
    });
}
