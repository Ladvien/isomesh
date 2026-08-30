//! **P-156 — a compact prefilter instead of a truncated recursive one.**
//!
//! Ticket: R-156. Pre-registered before this harness existed.
//!
//! ```bash
//! cargo bench --bench experiment_p156
//! ```
//!
//! Writes `docs/experiments/p-156.csv`.
//!
//! # What was missing
//!
//! `✗42`/`M-359` (`P-60`, `benches/experiment_p60.rs`) measured Blu, Thévenaz & Unser's **shifted
//! linear interpolation** and found the recursive prefilter's non-locality *bounded*: truncating the
//! causal one-pole recursion to `k` preceding samples moves the root under `1e-6` cells at `k ≥ 10`
//! on all eight fields, worst **`7.152528e-7`** on `torus` (`docs/experiments/p-60.csv:12`,
//! `guard_band_delta_cells` at `guard_band_k = 10`; `gyroid` reads `1.927142e-7`). What `✗42` did not
//! ask is whether the truncation is *necessary* — whether a prefilter that is compactly supported
//! **by construction** reaches the same accuracy, making the locality exact rather than bounded.
//!
//! `P-155` (`docs/experiments/p-155.csv`) is the other half of the setting, and both its verdicts are
//! quoted here because this row is built on them:
//!
//! - **C1 HELD.** `strang_fix_order = 2` on all eight rows, derived, and `measured_order = 2` beside
//!   it — so the tensor-product hat's approximation order is `2` as a derivation, not a fit.
//! - **C2 was FALSIFIED at 0 of 8** against a bar of four (`c2_hits = 0`, `c2_measurable_fields = 4`):
//!   the predicted asymptotic constant reproduced `M-12`'s fitted constant on **no** field. Four of
//!   eight fields carry no Hausdorff number at all (`skip_reason = bound_not_exact` on
//!   `csg_difference`, `gyroid`, `fbm_terrain`, `noise_cavity`), and two of the four that do are
//!   polyhedra whose creases have no second derivative.
//!
//! So P-155 licenses "order 2, derived" and explicitly does **not** license any predicted constant.
//! This harness therefore decides C1 on a **measured** comparison between prefilters at matched
//! order, never on a predicted constant.
//!
//! `M-32` is the other prior this row stands on: *chunk seams are bit-exact only when the cell size
//! is a power of two* — 16 of 16 shared-plane vertices agree bit-for-bit at `h = 0.125` and 0 of 14
//! at `h = 4/35`, because an extractor computes `origin + h·local` and chunk `c`'s last plane is
//! `(o + h·cn) + h·n` while `c+1`'s first is `o + h·(c+1)n`, equal by algebra and not by IEEE. C3 is
//! *"seams remain bit-exact **where `M-32` says they should be**"*, so the seam arm runs at a dyadic
//! `h` and a non-dyadic `h` is carried as the control that the instrument can see a failure.
//!
//! # The construction, and where it is taken from
//!
//! The primary is **Blu & Unser, *Approximation Error for Quasi-Interpolators and (Multi-)Wavelet
//! Expansions*, Applied and Computational Harmonic Analysis 6, 219–251 (1999),
//! `10.1006/acha.1998.0249`.** The corpus holds only the landing page, but the paper is **freely
//! available from the first author's own institutional page** and was read there:
//! `https://www.ee.cuhk.edu.hk/~tblu/monsite/pdfs/blu9801.pdf` (Thierry Blu, CUHK). Read from that
//! copy, and nothing below uses anything else from it:
//!
//! 1. **Definition 2, eq. (15) — the moment conditions.** A set of sampling distributions `φ̃ₙ` and
//!    synthesis functions `φₙ` is *quasi-biorthonormal of order `L`* iff the `φₙ` are of order `L`
//!    (they reproduce `x^s = Σₙ λₙ^(s) φₙ(x)` for `s = 0…L−1`) and the `φ̃ₙ` satisfy
//!    `∫ x^s φ̃ₙ(x) dx = λₙ^(s)` for `s = 0…L−1`. The paper's own gloss: *"The moment condition (15)
//!    on the sampling distributions is much less constraining and leaves room for many design
//!    alternatives."*
//! 2. **Theorem 3, eq. (30) — the equivalence.** Order-`L` quasi-biorthonormality ⟺ the `L₂`
//!    approximation error is `O(T^L)`. So the order is bought by the moment conditions alone; exact
//!    interpolation is not required, and the prefilter need not be the (non-compact) dual.
//! 3. **§III-C, the paragraph after Theorem 4 — the asymptotic constant.** *"we can constrain
//!    f̃̂(ν) = f̂(ν) + O(ν^{L+1}): this is equivalent to fixing the first `L+1` moments of `φ̃` to
//!    appropriate values provided by `φ`; then the first order of the asymptotic error equals the
//!    first order of the asymptotic minimum error … **In particular, it is now possible to consider
//!    compactly supported sampling functions associated with any kind of synthesis function.** This
//!    new freedom … takes its full sense when we remember that, in general, the dual functions of a
//!    set of synthesis functions are not compactly supported."*
//!
//! That last sentence is P-156's hypothesis, stated by the source: `L` moments buy the order, one
//! more moment buys the *constant*, and both can be had from a **compactly supported** analysis.
//!
//! A freely available restatement aimed at this crate's own audience was also read and is named
//! because it is the one a graphics reader should be sent to: **Nehab & Hoppe, *Generalized Sampling
//! in Computer Graphics*, MSR-TR-2011-16 / IMPA E022/2011**,
//! `https://www.microsoft.com/en-us/research/wp-content/uploads/2016/02/filtering.pdf`. Its §2.1
//! states the mechanism in one sentence — *"Quasi-interpolation schemes instead give up on
//! interpolation of ⟦fψ⟧, so that `q` is freed of this restriction. Instead, the interpolation
//! property holds only when `f` is a polynomial of degree less than `L`"* — and records that a **FIR**
//! design for `q` exists (Dalai, Leonardi & Migliorati 2005) beside Blu & Unser's IIR one.
//!
//! ## The derivation, in this crate's own notation
//!
//! The shifted-linear scheme `✗42` measured reconstructs `f_T(x) = Σₙ cₙ Λ(x/T − n − τ)` with `Λ` the
//! hat and knots at `(n + τ)T`; at `τ = 1/5` the exact prefilter is `1/H(z)` with
//! `H(z) = (1 − τ) + τ z⁻¹`, i.e. the causal recursion `cₙ = −POLE·cₙ₋₁ + GAIN·fₙ` with
//! `POLE = 2⁻²` and `GAIN = 1 + 2⁻²`, both exact binary. That filter has an **infinite** impulse
//! response, which is the whole locality problem.
//!
//! Write `δ = z⁻¹ − 1`, the backward difference operator, so `H = 1 + τδ` and
//!
//! ```text
//! 1/H = 1/(1 + τδ) = 1 − τδ + τ²δ² − τ³δ³ + …
//! ```
//!
//! Truncating that series **is** the quasi-interpolant, and each extra term fixes one more moment.
//! With `Mⱼ = Σ_lag lag^j · p_lag` for a causal FIR `p` at lags `0, 1, …`:
//!
//! | truncation | taps at lags `0,1,2` | `M₀` | `M₁` | `M₂` | reproduces |
//! |---|---|---|---|---|---|
//! | `1 − τδ` | `(1+τ, −τ)` | `1` | `−τ` | `−τ` | degree ≤ 1, **order 2** |
//! | `1 − τδ + τ²δ²` | `(1+τ+τ², −(τ+2τ²), τ²)` | `1` | `−τ` | `2τ² − τ` | degree ≤ 1, **order 2, and the exact filter's `M₂`** |
//! | the exact `1/H` | infinite | `1` | `−τ` | `2τ² − τ` | degree ≤ 1 |
//!
//! The three moment values on the right-hand column are forced, not chosen. Reproducing constants
//! needs `Σₙ Λ(x − n − τ) = 1` (partition of unity), so `λₙ^(0) = 1` and Definition 2 gives
//! `M₀ = 1`. Reproducing linears needs `Σₙ (n + τ) Λ(x − n − τ) = x`, so `λₙ^(1) = n + τ` and
//! `Σ_lag p_lag (n − lag) = n + τ` gives `M₁ = −τ`. Those two are order 2 — the hat's own maximum, so
//! *"full approximation order"* here means 2 and no prefilter can raise it (P-155 C1). The third
//! moment is §III-C's extra: the exact filter has `M₀ = H(1)⁻¹ = 1`, `M₁ = H′(1)⁻¹-form = −τ` and
//! `M₂ = H″ + H′ = 2τ² − τ`, so matching `M₂` matches the asymptotic constant as well as the order.
//!
//! And it does so exactly. For a locally quadratic `g` with second derivative `g″` in cell units, a
//! FIR with `M₀ = 1`, `M₁ = −τ` and second moment `M₂` reconstructs
//!
//! ```text
//! f_T(u) − g(u) = (g″/2)·[ s(1 − s) + M₂ − τ² ]
//! ```
//!
//! on the shifted segment carrying `u` at fraction `s` — derived by substituting
//! `cₙ = A + B(n+τ) + C(n² + 2τn + M₂)` into the two-knot linear piece and using
//! `interp(n²) = v² + s(1−s)` with `v = u − τ`. Setting `M₂ = 2τ² − τ` gives
//! `(g″/2)·[s(1−s) − τ(1−τ)]`, which is `✗42`'s own closed form for the *exact* prefilter
//! (experiment_p60.rs:102). Setting `M₂ = −τ` (the two-tap) gives `(g″/2)·[s(1−s) − τ(1+τ)]` — same
//! order, **different constant**, `0.24` against `0.16` at `τ = 1/5`.
//!
//! **So the harness ships the three-tap `(1+τ+τ², −(τ+2τ²), τ²) = (1.24, −0.28, 0.04)` as its answer
//! to C1, and carries the two-tap `(1.2, −0.2)` beside it as the minimal order-2 quasi-interpolant.**
//! The two-tap is what "order 2 with a compact filter" alone buys; the three-tap is what §III-C's
//! extra moment buys. Both are FIR, both are exactly local, and the CSV reports both so the 5% bar is
//! decided against a filter and not against a family.
//!
//! The taps are *not* a truncation of the exact impulse response, and the difference is load-bearing:
//! `1/H`'s first three taps are `(GAIN, −GAIN·POLE, GAIN·POLE²) = (1.25, −0.3125, 0.078125)`, which
//! sum to `1.015625` and therefore fail `M₀ = 1` — they do not even reproduce a constant. That tap
//! set is the harness's **negative control**.
//!
//! The 3-D prefilter is the tensor product of the 1-D one, applied separably along `x`, then `y`,
//! then `z`, which is what makes the whole scheme axis-separable and the footprint below a
//! per-axis statement.
//!
//! # Arms
//!
//! Five, all reconstructing on the same grid, differing **only** in the prefilter.
//!
//! | arm | `support` | prefilter | footprint below the bracket | `is_control` |
//! |---|---|---|---|---|
//! | `standard_trilinear` | `1` | identity; the hat already interpolates, so this is the shipped `t = a/(a − b)` at `τ = 0` | 2 samples | no |
//! | `recursive_whole_line` | `whole_line` | the exact `1/H(z)`, recursion started at the line's first sample with the paper's `c₀ = f₀` | the whole line | no |
//! | `recursive_truncated_k10` | `11` | the same recursion restarted `k = 10` samples before the bracket — **`✗42`'s C3 arm** | 11 samples | **yes** |
//! | `quasi_fir2` | `2` | `1 − τδ`: `(1+τ, −τ)` | 4 samples | no |
//! | `quasi_fir3` | `3` | `1 − τδ + τ²δ²`: `(1+τ+τ², −(τ+2τ²), τ²)` | 5 samples | no |
//!
//! The footprint is `depth + 1` samples ending at the bracket's upper sample `n`: the shifted
//! reconstruction reads `c_{n−2}, c_{n−1}, c_n` (`✗42`'s two-piece bracket, experiment_p60.rs:348),
//! and a causal `L`-tap FIR then reaches down to `f_{n−2−(L−1)}`, so `depth = L + 1`. `support` is the
//! **filter's** own length; `read_depth_samples` is the **vertex's** footprint. Both are columns.
//!
//! Two resolutions, `33` and `65` samples per axis. Both give a **dyadic** cell size on all eight
//! fields — `4/32`, `4/64` on the six `[−2,2]³` fields, `14/32 = 7/16` and `14/64 = 7/32` on
//! `gyroid`, `16/32` and `16/64` on `fbm_terrain` — which is `M-32`'s bit-exact regime, so C3 is
//! asked where `M-32` says the answer should be `true`. The pair also makes C1's ratio a
//! **series**: the three-tap agrees with the exact filter to `O(h³)` against an `O(h²)` error, so the
//! ratio's distance from 1 must fall roughly as `h` and the finer row is the one that decides.
//!
//! # SHARE, recomputed before the numbers
//!
//! The registration's SHARE is *"C1 moves the vertex-interpolation stage, the crate's single
//! most-executed operation."* Recomputed: the operation is `edge_position` in
//! `marching_cubes/mod.rs`, called once per cut grid edge, and `prefilter_ms` prices what a prefilter
//! would add in front of it. A compact FIR is `L` multiply-adds per sample per axis and needs a
//! guard band of `read_depth_samples − 1`; the truncated recursion is `k` sequential dependent
//! multiply-adds per **vertex** and needs a guard band of `k`. So SHARE is discharged as: if C1 holds,
//! the shipped `t = a/(a − b)` can be replaced by a three-tap FIR plus the same crossing rule, with a
//! four-sample guard band instead of a ten-sample one and no sequential dependence — and if C1 fails,
//! it cannot, and the truncated recursion is the better engineering answer, which is exactly what the
//! registered falsifier says.
//!
//! # No shipped path changes, so `hashes_moved` is predicted 0 — and measured
//!
//! Every mechanism here is bench-local. `crates/isomesh/src/**` is untouched, no reference field is
//! added, and the only crate code called is `Sdf::sample`, `Sdf::gradient`, `MarchingCubes::extract`
//! and `validate::accuracy`, none of which is modified. `crates/isomesh/golden_hashes.json` is
//! **read** by a one-line scanner in the shape of `golden.rs`'s own private `field_of` and never
//! written, so `hashes_moved` is `0` on every row by construction. That zero is only a measurement if
//! the fixture it counts over is real, so the harness parses the file and asserts **216 rows over 8
//! fields × 9 algorithms × 3 resolutions** before recording anything (`M-44`'s rule).
//!
//! # Which columns decide which clause
//!
//! - **C1** — `root_position_error` is the RMS, in cells, of `|u_arm − u_exact|` over every cut `x`
//!   edge of the field's own grid, with `u_exact` bisected from the analytic field.
//!   `vs_truncated_recursive` is that number over `recursive_truncated_k10`'s. C1 is decided on
//!   `quasi_fir3` at the **finer** resolution, on every field where the baseline error is above the
//!   numerical floor, against `|ratio − 1| ≤ 0.05`. `c1_fields_within_5pct` and `c1_population` are
//!   columns so a reader can apply a looser population bar without recomputing. `hausdorff` is the
//!   surface-level corroboration and cannot decide the clause, because `validate::accuracy` is
//!   meaningless where `field.bound()` is not `Exact` and that is four of the eight fields
//!   (`hausdorff_skip` says which).
//! - **C2** — `chunk_local` is **structural**: it is `true` iff the prefilter's impulse response is
//!   finite, which is a property of the tap list and not of any measurement. The registration demands
//!   exactly that, and the harness measures **why** it has to: `footprint_probe_outside` perturbs each
//!   sample of the line by `12.5%` of the line's scale, sign-preserving, and counts the samples
//!   outside the declared footprint whose perturbation moves the arm's root — `0` for every FIR arm,
//!   non-zero for `recursive_whole_line`. But on an exactly linear restriction the recursion's startup
//!   transient cancels (`✗42`'s own note, experiment_p60.rs:157-164), so a *measured* locality check
//!   can pass for a non-local filter. Structural is the only sound assertion, and this is the number
//!   that says so.
//! - **C3** — `seam_bit_exact` splits the field's `x` range into two overlapping chunks, samples each
//!   from **its own origin** with **its own local indices** (which is `M-32`'s mechanism, not a
//!   simulation of it), and asks whether every cut bracket in the ten-sample overlap gets a
//!   bit-identical world coordinate from both. `m32_control_seam_bit_exact` is the same test at a
//!   non-dyadic `h` (35 samples, `4/34 = 2/17`), where `M-32` predicts `false`, so the `true` above is
//!   not a pass over an untested configuration. `hashes_moved` is the second half of C3.
//!
//! `residual_vs_whole_line_cells` is the column the hypothesis's *"exact rather than bounded"* lives
//! in: the RMS distance from an arm's root to the exact `1/H(z)` root. It is `0` for
//! `recursive_whole_line`, `~1e-7` for `recursive_truncated_k10` (which is `✗42`'s bound, re-measured
//! over thousands of edges rather than eight), and a genuinely different number for the FIR arms,
//! because they are different filters rather than approximations of that one.
//!
//! # Vacuity controls
//!
//! Every one runs before the first `run.record`, and every panic message starts `VOID: `.
//!
//! - **`✗42`'s `7.152528e-7` is reproduced**, on `✗42`'s own fixture — `torus`, the domain-centre `x`
//!   axis, 65 samples, `τ = 1/5`, `k = 10` — to within `1e-6` relative, which is the resolution of the
//!   seven significant digits the ledger quotes. `gyroid`'s `1.927142e-7` too. Column:
//!   `x42_torus_k10_delta`, `x42_gyroid_k10_delta`. Without this the whole comparison has no baseline,
//!   which is what the registration says.
//! - **The moment identities.** `M₀ = 1` and `M₁ = −τ` on both FIR arms, `M₂ = 2τ² − τ` on the
//!   three-tap and `M₂ = −τ` on the two-tap — asserted, so the two arms are provably different filters
//!   and not two spellings of one. Columns `moment_0`, `moment_1`, `moment_2`.
//! - **Polynomial reproduction, positive control.** On an exactly linear ramp with its root at a
//!   general position, every arm except `recursive_truncated_k10` puts the root within `1e-12` cells of
//!   the truth; the truncated arm does not, and its residue is its startup transient. Column
//!   `ramp_root_error`.
//! - **The moment conditions are load-bearing, negative control.** The exact filter's first three
//!   impulse-response taps `(1.25, −0.3125, 0.078125)` fail `M₀ = 1`, and fed through the identical
//!   code path they miss the linear ramp's root by more than `1e-3` cells — nine orders worse than the
//!   quasi-interpolants. A value near theirs would mean the moment conditions are doing nothing.
//!   Column `sabotage_ramp_root_error`.
//! - **The mesh detour is faithful.** `hausdorff` is measured by handing `MarchingCubes` the
//!   *coefficient* grid at the true origin and translating the result by `+τ·h` per axis, which is
//!   exact because the shifted reconstruction's zero set is the trilinear interpolant of `c` shifted by
//!   `τ`. With the identity prefilter and `τ = 0` that detour must return the **bit-identical** mesh to
//!   `MarchingCubes` run on the field itself, positions and indices. Column
//!   `standard_mesh_bit_identical`.
//! - **The populations are non-empty.** `edges_measured > 0` and `seam_pairs > 0` on every field, and
//!   the baseline `root_position_error` above `1e-9` cells on every field, or the ratios are arithmetic
//!   rather than measurement.
//! - **The truncation residue is non-zero.** `residual_vs_whole_line_cells > 0` for
//!   `recursive_truncated_k10` on at least one field, or *"bounded rather than exact"* is unmeasured.
//! - **The locality probe can see a dependence.** `footprint_probe_outside > 0` for
//!   `recursive_whole_line` on at least one field, or the `0` on the FIR arms is a zero that could not
//!   have been non-zero.
//! - **`M-32`'s failure is reachable.** `m32_control_seam_bit_exact` is `false` for at least one
//!   (arm, field) at the non-dyadic `h`, or C3's `true` at the dyadic one is untested.
//! - **The golden fixture is the real one.** 216 rows, 8 fields, 9 algorithms, 3 resolutions.
//!
//! # Determinism and timing
//!
//! No RNG: every fixture is a reference field, a grid, or a closed-form ramp. Sorting is by
//! `f64::total_cmp`. `prefilter_ms` is the median of `5` repeats of the full separable prefilter over
//! the field's own grid after one warm-up, with min and max beside it, because this host's
//! `amd-pstate-epp` governor swings the same binary `1.45×` between runs (`M-280`). No clause is
//! decided by a clock.

#![allow(
    clippy::cast_possible_truncation,
    clippy::cast_precision_loss,
    clippy::float_cmp,
    clippy::too_many_lines
)]

mod common;

use std::hint::black_box;
use std::time::Instant;

use isomesh::fields::ReferenceField;
use isomesh::marching_cubes::MarchingCubes;
use isomesh::validate::{AccuracyConfig, accuracy};
use isomesh::{MeshBuffer, RuntimeShape3, Sdf, for_each_reference_field};

// ──────────────────────────────────────────────────────────────────────────────
// Transcribed constants
// ──────────────────────────────────────────────────────────────────────────────

/// **A sample of exactly zero is outside.**
///
/// Transcribed from `isomesh::cube::is_inside`, which a bench cannot reach because `cube` is a
/// private module, exactly as `experiment_p60.rs:233` transcribes it. Strict `< 0`.
fn is_inside(value: f64) -> bool {
    value < 0.0
}

/// The registered knot shift, `1/5`.
///
/// The nearest `f64` to `1/5`, which is the only place `1/5` is rounded: [`POLE`] and [`GAIN`] are
/// exact binary.
const TAU: f64 = 0.2;

/// `τ / (1 − τ)` at `τ = 1/5`: `2⁻²`, exact. Transcribed from `experiment_p60.rs:251`.
const POLE: f64 = 0.25;

/// `1 / (1 − τ)` at `τ = 1/5`: `1 + 2⁻²`, exact. Transcribed from `experiment_p60.rs:254`.
const GAIN: f64 = 1.25;

/// `✗42`'s guard-band length, which is the truncated arm's window.
const GUARD_K: usize = 10;

/// The minimal order-2 quasi-interpolant, `1 − τδ`, at lags `0, 1`.
const QUASI_FIR2: [f64; 2] = [1.0 + TAU, -TAU];

/// The order-2 quasi-interpolant that also matches the exact filter's `M₂`, `1 − τδ + τ²δ²`.
const QUASI_FIR3: [f64; 3] = [1.0 + TAU + TAU * TAU, -(TAU + 2.0 * TAU * TAU), TAU * TAU];

/// The **negative control**: the exact filter's first three impulse-response taps.
///
/// `(GAIN, −GAIN·POLE, GAIN·POLE²) = (1.25, −0.3125, 0.078125)`, summing to `1.015625`, so `M₀ ≠ 1`
/// and it does not reproduce a constant. This is the most likely way to get a compact prefilter
/// wrong: truncate the impulse response instead of fixing the moments.
const SABOTAGE_FIR3: [f64; 3] = [GAIN, -GAIN * POLE, GAIN * POLE * POLE];

// ──────────────────────────────────────────────────────────────────────────────
// Fixture constants
// ──────────────────────────────────────────────────────────────────────────────

/// Samples per axis. Both give a dyadic cell size on all eight fields (`M-32`'s regime).
const RESOLUTIONS: [u32; 2] = [33, 65];

/// The lowest bracket index the population accepts.
///
/// Every arm must get its whole declared footprint on every measured edge, and the deepest bounded
/// footprint is the truncated recursion's [`GUARD_K`]. Brackets below this are counted in
/// `edges_dropped_left_context` rather than measured with a short window.
const MIN_BRACKET: usize = GUARD_K;

/// Samples per axis for the `M-32` control seam, giving `h = 4/34 = 2/17` — not dyadic.
const M32_CONTROL_SAMPLES: u32 = 35;

/// Timed repeats of the prefilter pass, after one warm-up.
const TIMING_REPEATS: usize = 5;

/// `✗42`'s sample count along its line.
const X42_SAMPLES: usize = 65;

/// `✗42`'s `torus` figure at `k = 10`, from `docs/experiments/p-60.csv:12`.
const X42_TORUS_K10: f64 = 7.152528e-7;

/// `✗42`'s `gyroid` figure at `k = 10`, from `FINDINGS.md:9820`.
const X42_GYROID_K10: f64 = 1.927142e-7;

/// Relative tolerance on the two figures above: the ledger quotes seven significant digits.
const X42_TOLERANCE: f64 = 1e-6;

/// Requested bisection width for the exact root. Transcribed from `experiment_p60.rs:266`.
const BISECT_WIDTH: f64 = 1e-15;

/// Samples on the synthetic linear ramp used by both instrument controls.
const RAMP_SAMPLES: usize = 65;

/// Where the ramp's root sits, in cells from the ramp's start.
///
/// Deliberately general: not a half-integer, not `1 − τ` of a shifted segment along (which is where
/// `✗42` records that the recursion's startup transient cancels exactly), and far enough from the
/// start that the whole-line recursion's own transient is below `f64` resolution.
const RAMP_ROOT: f64 = 40.37;

/// A filter that reproduces linears must put the ramp's root this close, in cells.
const RAMP_FLOOR: f64 = 1e-12;

/// The sabotage taps must miss the ramp's root by more than this, in cells.
const SABOTAGE_MIN_ERROR: f64 = 1e-3;

/// Tolerance on the moment identities, in absolute value.
const MOMENT_TOLERANCE: f64 = 1e-15;

/// The locality probe's perturbation, relative to the line's own scale.
///
/// Large on purpose. A small perturbation is invisible after the recursion's `POLE^j` decay, and the
/// probe would then report "no dependence" for a filter that has one.
const PROBE_DELTA_REL: f64 = 0.125;

/// How far a sample may sit from the lattice before [`GridSdf`] refuses it.
const LATTICE_TOLERANCE: f64 = 1e-6;

/// Below this many cells the baseline is at its own numerical floor and a ratio against it is
/// arithmetic rather than measurement.
const FLOOR_CELLS: f64 = 1e-9;

/// C1's bar: "within 5%".
const C1_TOLERANCE: f64 = 0.05;

/// The committed golden fixture's shape (§8 of the API cheat sheet, and `golden.rs` itself).
const GOLDEN_ROWS: usize = 216;
/// Distinct field names in the golden fixture.
const GOLDEN_FIELDS: usize = 8;
/// Distinct algorithm names in the golden fixture.
const GOLDEN_ALGORITHMS: usize = 9;
/// Distinct resolutions in the golden fixture.
const GOLDEN_RESOLUTIONS: usize = 3;

// ──────────────────────────────────────────────────────────────────────────────
// Prefilters
// ──────────────────────────────────────────────────────────────────────────────

/// A prefilter, as the only thing the five arms differ in.
#[derive(Clone, Copy)]
enum Filter {
    /// No prefilter at all: the hat interpolates, so the shipped rule needs none.
    Identity,
    /// A causal FIR with taps at lags `0 .. taps.len()`.
    Fir(&'static [f64]),
    /// The causal one-pole recursion. `Some(k)` restarts it `k` samples before the bracket, `None`
    /// starts it at the first sample it is given.
    Recursive(Option<usize>),
}

/// One arm.
struct Arm {
    /// The `prefilter` column.
    name: &'static str,
    /// The `support` column: the **filter's** own length.
    support: &'static str,
    /// The prefilter.
    filter: Filter,
    /// Samples below the bracket's upper sample `n` that the vertex reads. `None` is unbounded.
    depth: Option<usize>,
    /// Whether the impulse response is finite. This **is** `chunk_local`, structurally.
    finite_impulse: bool,
    /// Cells the reconstruction's knots are shifted by, for the mesh arm.
    shift: f64,
    /// The `is_control` column.
    is_control: bool,
}

/// The five arms, in CSV order.
const ARMS: [Arm; 5] = [
    Arm {
        name: "standard_trilinear",
        support: "1",
        filter: Filter::Identity,
        depth: Some(1),
        finite_impulse: true,
        shift: 0.0,
        is_control: false,
    },
    Arm {
        name: "recursive_whole_line",
        support: "whole_line",
        filter: Filter::Recursive(None),
        depth: None,
        finite_impulse: false,
        shift: TAU,
        is_control: false,
    },
    Arm {
        name: "recursive_truncated_k10",
        support: "11",
        filter: Filter::Recursive(Some(GUARD_K)),
        depth: Some(GUARD_K),
        finite_impulse: true,
        shift: TAU,
        is_control: true,
    },
    Arm {
        name: "quasi_fir2",
        support: "2",
        filter: Filter::Fir(&QUASI_FIR2),
        depth: Some(QUASI_FIR2.len() + 1),
        finite_impulse: true,
        shift: TAU,
        is_control: false,
    },
    Arm {
        name: "quasi_fir3",
        support: "3",
        filter: Filter::Fir(&QUASI_FIR3),
        depth: Some(QUASI_FIR3.len() + 1),
        finite_impulse: true,
        shift: TAU,
        is_control: false,
    },
];

/// Index of the exact-filter arm in [`ARMS`].
const ARM_WHOLE: usize = 1;
/// Index of the baseline arm in [`ARMS`], which every ratio divides by.
const ARM_TRUNCATED: usize = 2;
/// Index of the arm C1 is decided on.
const ARM_FIR3: usize = 4;

/// `Σ lag^j · tap` for a causal FIR, or the exact filter's closed form.
///
/// The exact filter's moments come from `1/H(w) = 1/((1−τ) + τw) = Σ hₖ wᵏ`: `M₀ = H(1)⁻¹ = 1`,
/// `M₁ = Σ k hₖ = −τ/((1−τ)+τ)² = −τ`, and `Σ k(k−1) hₖ = 2τ²` so `M₂ = 2τ² − τ`.
fn moments(filter: Filter) -> [f64; 3] {
    match filter {
        Filter::Identity => [1.0, 0.0, 0.0],
        Filter::Fir(taps) => {
            let mut m = [0.0f64; 3];
            for (lag, tap) in taps.iter().enumerate() {
                let l = lag as f64;
                m[0] += tap;
                m[1] += l * tap;
                m[2] += l * l * tap;
            }
            m
        }
        Filter::Recursive(_) => [1.0, -TAU, 2.0 * TAU * TAU - TAU],
    }
}

/// One coefficient of a causal FIR at index `m`, replicating the first sample below the line.
fn fir_at(f: &[f64], taps: &[f64], m: usize) -> f64 {
    let mut acc = 0.0;
    for (lag, tap) in taps.iter().enumerate() {
        acc += tap * f[m.saturating_sub(lag)];
    }
    acc
}

/// `(c_{n−2}, c_{n−1}, c_n)` of the causal recursion started at `start`.
///
/// `c_start = f_start` is Blu, Thévenaz & Unser's own `c₀ = f₀`, applied at whichever sample the
/// window begins. Transcribed from `experiment_p60.rs:315-327`, which is the code that produced
/// `✗42`'s figures, and the expression order is kept so the reproduction is bit-for-bit.
fn recursive_triple(f: &[f64], start: usize, n: usize) -> (f64, f64, f64) {
    let mut c = f[start];
    let mut c_prev2 = f64::NAN;
    let mut c_prev1 = f64::NAN;
    for value in &f[(start + 1)..=n] {
        c_prev2 = c_prev1;
        c_prev1 = c;
        c = GAIN * value - POLE * c;
    }
    (c_prev2, c_prev1, c)
}

/// The three coefficients the shifted reconstruction reads for the bracket `[n − 1, n]`.
///
/// For [`Filter::Identity`] the coefficients *are* the samples, which is what "no prefilter" means;
/// that arm does not use the shifted reconstruction and this triple exists only so the footprint
/// bookkeeping is uniform.
fn arm_triple(filter: Filter, f: &[f64], n: usize) -> (f64, f64, f64) {
    match filter {
        Filter::Identity => (f[n - 2], f[n - 1], f[n]),
        Filter::Fir(taps) => (
            fir_at(f, taps, n - 2),
            fir_at(f, taps, n - 1),
            fir_at(f, taps, n),
        ),
        Filter::Recursive(window) => {
            let start = match window {
                Some(k) => n.saturating_sub(k),
                None => 0,
            };
            recursive_triple(f, start, n)
        }
    }
}

/// The zero of the shifted reconstruction inside the sample bracket `[n − 1, n]`, in cells.
///
/// Transcribed from `experiment_p60.rs:348-374`. The bracket is split by a knot at `n − 1 + τ`, so it
/// is covered by two linear pieces; `f_T` equals the reconstruction's value at each sample, and
/// exactly one piece changes sign by parity, so the branch is total and needs no epsilon — each
/// chosen denominator is non-zero because a convex combination of two same-side values cannot cross.
///
/// Returns `None` when the reconstruction does **not** change sign across a bracket the field does.
/// That is a real event for a quasi-interpolant, which does not interpolate the samples: it is a
/// crossing this prefilter cannot see, and it is counted rather than hidden.
fn shifted_root(n: usize, c2: f64, c1: f64, c0: f64) -> Option<f64> {
    let v_left = TAU * c2 + (1.0 - TAU) * c1;
    let v_right = TAU * c1 + (1.0 - TAU) * c0;
    if is_inside(v_left) == is_inside(v_right) {
        return None;
    }
    if is_inside(v_left) == is_inside(c1) {
        Some((n - 1) as f64 + TAU + c1 / (c1 - c0))
    } else {
        Some((n - 2) as f64 + TAU + c2 / (c2 - c1))
    }
}

/// Where an arm puts the vertex on the bracket `[n − 1, n]` of the sample line `f`, in cells from
/// `f`'s own first sample.
///
/// `f` is the **chunk's own** slice and `n` its **own** local index, which is what makes the seam
/// arm faithful: a chunk computes with its own indices and its own origin.
fn arm_root(filter: Filter, f: &[f64], n: usize) -> Option<f64> {
    if matches!(filter, Filter::Identity) {
        let a = f[n - 1];
        let b = f[n];
        return Some((n - 1) as f64 + a / (a - b));
    }
    let (c2, c1, c0) = arm_triple(filter, f, n);
    shifted_root(n, c2, c1, c0)
}

/// The whole coefficient line of a prefilter, for the mesh arm.
///
/// The boundary rule below the line's first sample is replication for a FIR and window saturation for
/// the recursion (the latter is `✗42`'s own `start = n.saturating_sub(k)`). `mesh_min_bracket` records
/// how close the surface comes to the boundary so a reader can see whether the rule was reachable.
fn coefficients_1d(filter: Filter, f: &[f64], out: &mut [f64]) {
    match filter {
        Filter::Identity => out.copy_from_slice(f),
        Filter::Fir(taps) => {
            for (m, slot) in out.iter_mut().enumerate() {
                *slot = fir_at(f, taps, m);
            }
        }
        Filter::Recursive(None) => {
            out[0] = f[0];
            for m in 1..f.len() {
                out[m] = GAIN * f[m] - POLE * out[m - 1];
            }
        }
        Filter::Recursive(Some(k)) => {
            for (m, slot) in out.iter_mut().enumerate() {
                let start = m.saturating_sub(k);
                let mut c = f[start];
                for value in &f[(start + 1)..=m] {
                    c = GAIN * value - POLE * c;
                }
                *slot = c;
            }
        }
    }
}

/// The tensor-product prefilter over a cubic sample grid, applied separably along `x`, `y`, `z`.
fn coefficient_grid(filter: Filter, values: &[f64], n: usize) -> Vec<f64> {
    let mut grid = values.to_vec();
    let mut line = vec![0.0f64; n];
    let mut out = vec![0.0f64; n];
    for axis in 0..3usize {
        let stride = match axis {
            0 => 1,
            1 => n,
            _ => n * n,
        };
        for a in 0..n {
            for b in 0..n {
                let base = match axis {
                    0 => a * n + b * n * n,
                    1 => a + b * n * n,
                    _ => a + b * n,
                };
                for (i, slot) in line.iter_mut().enumerate() {
                    *slot = grid[base + i * stride];
                }
                coefficients_1d(filter, &line, &mut out);
                for (i, value) in out.iter().enumerate() {
                    grid[base + i * stride] = *value;
                }
            }
        }
    }
    grid
}

// ──────────────────────────────────────────────────────────────────────────────
// Helpers
// ──────────────────────────────────────────────────────────────────────────────

/// Bisect `field` along a line parallel to `x` between two abscissae of opposite side.
///
/// Terminates at [`BISECT_WIDTH`] or at the last representable interval, whichever comes first.
/// Transcribed from `experiment_p60.rs:391-408`.
fn bisect<F>(field: &F, y: f64, z: f64, mut lo: f64, mut hi: f64) -> f64
where
    F: Sdf<Scalar = f64>,
{
    let lo_inside = is_inside(field.sample([lo, y, z]));
    loop {
        let mid = f64::midpoint(lo, hi);
        if mid <= lo || mid >= hi || hi - lo <= BISECT_WIDTH {
            break;
        }
        if is_inside(field.sample([mid, y, z])) == lo_inside {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    f64::midpoint(lo, hi)
}

/// Median by `total_cmp`, mean of the two middle values on an even population.
fn median(values: &[f64]) -> f64 {
    assert!(!values.is_empty(), "median of an empty population");
    let mut v = values.to_vec();
    v.sort_unstable_by(f64::total_cmp);
    let m = v.len() / 2;
    if v.len().is_multiple_of(2) {
        f64::midpoint(v[m - 1], v[m])
    } else {
        v[m]
    }
}

/// `{:.6e}`, which keeps `inf` and `NaN` visible rather than rounded into a number.
fn sci(v: f64) -> String {
    format!("{v:.6e}")
}

// ──────────────────────────────────────────────────────────────────────────────
// The `hausdorff` mesh route
// ──────────────────────────────────────────────────────────────────────────────

/// A coefficient grid presented as an `Sdf`, so the crate's own `MarchingCubes` can mesh it.
///
/// The shifted reconstruction `Σ c_n Λ(x/T − n − τ)` is the trilinear interpolant of `c` evaluated at
/// `x/T − τ`, so its zero set is the zero set of the trilinear interpolant of `c` translated by
/// `+τ·T` on each axis. Meshing `c` at the true origin and translating the result is therefore
/// **exact** rather than an approximation, and with the identity prefilter and `τ = 0` it must return
/// the crate's own mesh bit-for-bit — which is `standard_mesh_bit_identical`.
///
/// `sample` is a lattice lookup and asserts the point it was handed is on the lattice.
/// `MarchingCubes` gathers corners through `crate::sdf::sample_grid` at `origin + cell·i` exactly
/// (`sdf.rs:183-187`), so that assertion holds for every corner it asks for. `gradient` is delegated
/// to the underlying field because `MarchingCubes` evaluates it at the **vertex** position, which is
/// off the lattice by construction; no normal is measured by this experiment.
struct GridSdf<'a, F> {
    /// The coefficient grid, `x` fastest.
    values: &'a [f64],
    /// Samples per axis.
    size: usize,
    /// The grid's origin.
    origin: [f64; 3],
    /// The grid's cell size.
    cell: f64,
    /// The field the grid was sampled from, for `gradient` only.
    field: &'a F,
}

impl<F: Sdf<Scalar = f64>> Sdf for GridSdf<'_, F> {
    type Scalar = f64;

    fn sample(&self, p: [f64; 3]) -> f64 {
        let mut index = 0usize;
        for (axis, coordinate) in p.iter().enumerate().rev() {
            let q = (coordinate - self.origin[axis]) / self.cell;
            let r = q.round();
            assert!(
                (q - r).abs() < LATTICE_TOLERANCE,
                "GridSdf was asked for {coordinate} on axis {axis}, which is {q} cells from the \
                 origin and not on the lattice; the extractor is sampling off-grid and this \
                 experiment's mesh arm would be measuring the detour rather than the prefilter"
            );
            assert!(
                r >= 0.0 && (r as usize) < self.size,
                "GridSdf was asked for cell {r} on axis {axis}, outside 0..{}",
                self.size
            );
            index = index * self.size + (r as usize);
        }
        self.values[index]
    }

    fn gradient(&self, p: [f64; 3]) -> [f64; 3] {
        self.field.gradient(p)
    }
}

// ──────────────────────────────────────────────────────────────────────────────
// The golden fixture
// ──────────────────────────────────────────────────────────────────────────────

/// What the committed golden fixture says, and therefore what `hashes_moved = 0` is a zero over.
struct Golden {
    /// Data rows.
    rows: usize,
    /// Distinct `field` values.
    fields: usize,
    /// Distinct `algorithm` values.
    algorithms: usize,
    /// Distinct `samples` values.
    resolutions: usize,
}

/// One value out of one line of `golden_hashes.json`.
///
/// In the shape of `golden.rs`'s own private `field_of`: find `"key":`, strip an optional quote, cut
/// at the first of `"`, `,`, `}`. The file has fixed key order, one object per line, no nesting and
/// no escapes, so this is the whole reader.
fn field_of<'a>(line: &'a str, key: &str) -> Option<&'a str> {
    let needle = format!("\"{key}\":");
    let start = line.find(&needle)? + needle.len();
    let rest = &line[start..];
    let rest = rest.strip_prefix('"').unwrap_or(rest);
    let end = rest.find(['"', ',', '}']).unwrap_or(rest.len());
    Some(&rest[..end])
}

/// Read the committed fixture rather than re-deriving it.
fn golden() -> Golden {
    let path = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("golden_hashes.json");
    let text = std::fs::read_to_string(&path)
        .expect("golden_hashes.json is committed beside crates/isomesh/Cargo.toml");
    let mut fields: Vec<&str> = Vec::new();
    let mut algorithms: Vec<&str> = Vec::new();
    let mut resolutions: Vec<&str> = Vec::new();
    let mut rows = 0usize;
    for line in text.lines() {
        let Some(algorithm) = field_of(line, "algorithm") else {
            continue;
        };
        let field = field_of(line, "field").expect("a golden row names its field");
        let samples = field_of(line, "samples").expect("a golden row names its resolution");
        rows += 1;
        if !algorithms.contains(&algorithm) {
            algorithms.push(algorithm);
        }
        if !fields.contains(&field) {
            fields.push(field);
        }
        if !resolutions.contains(&samples) {
            resolutions.push(samples);
        }
    }
    Golden {
        rows,
        fields: fields.len(),
        algorithms: algorithms.len(),
        resolutions: resolutions.len(),
    }
}

// ──────────────────────────────────────────────────────────────────────────────
// The instrument controls
// ──────────────────────────────────────────────────────────────────────────────

/// `✗42`'s guard-band delta at window `k`, on `✗42`'s own fixture.
///
/// The domain-centre `x` axis, [`X42_SAMPLES`] samples, the first crossing bracket, and
/// `|u_truncated − u_whole_line|` in cells — `guard_band_delta_cells` of `docs/experiments/p-60.csv`.
/// Bisection does not enter it: it is a difference between two filters on the same samples.
fn x42_guard_delta<F>(name: &str, field: &F, k: usize) -> f64
where
    F: ReferenceField + Sdf<Scalar = f64>,
{
    let (lo, hi) = field.domain();
    let x0 = lo[0];
    let cell = (hi[0] - x0) / (X42_SAMPLES - 1) as f64;
    let f: Vec<f64> = (0..X42_SAMPLES)
        .map(|i| field.sample([x0 + cell * i as f64, 0.0, 0.0]))
        .collect();
    let n = (1..X42_SAMPLES)
        .find(|&i| is_inside(f[i - 1]) != is_inside(f[i]))
        .unwrap_or_else(|| {
            panic!("VOID: {name} has no crossing on x42's line, so there is nothing to reproduce")
        });
    assert!(
        n >= 2,
        "VOID: {name}'s first crossing is in bracket [{}, {n}], which has no sample before it for \
         the shifted reconstruction's left knot",
        n - 1
    );
    let whole = recursive_triple(&f, 0, n);
    let u_whole = shifted_root(n, whole.0, whole.1, whole.2)
        .unwrap_or_else(|| panic!("VOID: the exact filter loses {name}'s first crossing"));
    let truncated = recursive_triple(&f, n.saturating_sub(k), n);
    let u_truncated = shifted_root(n, truncated.0, truncated.1, truncated.2)
        .unwrap_or_else(|| panic!("VOID: the truncated filter loses {name}'s first crossing"));
    (u_truncated - u_whole).abs()
}

/// A synthetic exactly-linear ramp, and where each filter puts its root.
///
/// A filter with `M₀ = 1` and `M₁ = −τ` reproduces degree-1 polynomials, so its root is the true one
/// up to `f64`. This pins the taps, the knot offset, the two-piece segment choice and the sign of the
/// pole all at once, and it is the fixture the sabotage taps are held against.
fn ramp_root_error(filter: Filter) -> f64 {
    let f: Vec<f64> = (0..RAMP_SAMPLES).map(|i| i as f64 - RAMP_ROOT).collect();
    let n = RAMP_ROOT.ceil() as usize;
    let u = arm_root(filter, &f, n)
        .expect("a filter that reproduces linears cannot lose a linear ramp's crossing");
    (u - RAMP_ROOT).abs()
}

// ──────────────────────────────────────────────────────────────────────────────
// Per-field measurement
// ──────────────────────────────────────────────────────────────────────────────

/// The seam arm's answer for one arm at one cell size.
#[derive(Clone, Copy)]
struct Seam {
    /// Cut brackets both chunks agreed were cut, in the overlap.
    pairs: u64,
    /// Whether every one of those got a bit-identical **world** coordinate from both chunks.
    world_bit_exact: bool,
    /// Whether every one got a bit-identical coordinate in each chunk's own **local** index space,
    /// which is the prefilter's own business and excludes `M-32`'s origin arithmetic.
    local_bit_exact: bool,
    /// Worst world disagreement, in cells.
    worst: f64,
    /// Brackets one chunk called cut and the other did not, which at a non-dyadic `h` is `M-44`'s
    /// class of seam sample disagreement rather than a prefilter defect.
    sign_disagreements: u64,
}

/// Everything measured for one arm on one field at one resolution.
struct ArmRow {
    /// RMS `|u_arm − u_exact|` over the common edge population, in cells.
    root_rms: f64,
    /// Worst `|u_arm − u_exact|`, in cells.
    root_max: f64,
    /// RMS `|u_arm − u_whole_line|`, in cells: the "exact rather than bounded" column.
    residual_rms: f64,
    /// Cut brackets this arm's reconstruction did not see.
    lost: u64,
    /// Samples outside the declared footprint whose perturbation moves the root.
    probe_outside: u64,
    /// Samples the probe declined to perturb because it would have flipped their side.
    probe_skipped: u64,
    /// Distance between the lowest and highest sample index the probe found a dependence on.
    probe_span: u64,
    /// The seam arm at the row's own (dyadic) cell size.
    seam: Seam,
    /// The seam arm at the non-dyadic `M-32` control cell size.
    m32: Seam,
    /// Symmetric Hausdorff of this arm's mesh, or `NaN` where the field's bound is not exact.
    hausdorff: f64,
    /// Median of [`TIMING_REPEATS`] repeats of the separable prefilter over the grid, in ms.
    prefilter_ms: f64,
    /// Fastest repeat, in ms.
    prefilter_ms_min: f64,
    /// Slowest repeat, in ms.
    prefilter_ms_max: f64,
}

/// Everything measured for one field at one resolution.
struct FieldRow {
    /// The field's own name.
    name: &'static str,
    /// Samples per axis.
    samples: u32,
    /// The cell size, which decides whether `M-32` predicts a bit-exact seam.
    cell: f64,
    /// Edges in the common population.
    edges: u64,
    /// Cut brackets dropped because some arm would not have had its whole footprint.
    edges_dropped: u64,
    /// The lowest cut bracket index anywhere in the grid, over all three axes.
    mesh_min_bracket: usize,
    /// Why `hausdorff` is `NaN`, or `none`.
    hausdorff_skip: &'static str,
    /// Whether the identity arm's mesh detour returned the crate's own mesh bit-for-bit.
    mesh_bit_identical: bool,
    /// Per-arm results, in [`ARMS`] order.
    arms: Vec<ArmRow>,
}

/// Sample a cubic grid over the field's own domain, `x` fastest.
fn sample_grid_cubic<F>(field: &F, lo: [f64; 3], cell: f64, n: usize) -> Vec<f64>
where
    F: Sdf<Scalar = f64>,
{
    let mut values = Vec::with_capacity(n * n * n);
    for k in 0..n {
        let z = lo[2] + cell * k as f64;
        for j in 0..n {
            let y = lo[1] + cell * j as f64;
            for i in 0..n {
                values.push(field.sample([lo[0] + cell * i as f64, y, z]));
            }
        }
    }
    values
}

/// The lowest cut bracket index in the grid, over all three axis families.
///
/// It says how close the surface comes to the domain boundary, and therefore whether the coefficient
/// grid's boundary rule is reachable from any cut edge at all.
fn min_cut_bracket(values: &[f64], n: usize) -> usize {
    let mut lowest = n;
    for axis in 0..3usize {
        let stride = match axis {
            0 => 1,
            1 => n,
            _ => n * n,
        };
        for a in 0..n {
            for b in 0..n {
                let base = match axis {
                    0 => a * n + b * n * n,
                    1 => a + b * n * n,
                    _ => a + b * n,
                };
                for m in 1..n {
                    if is_inside(values[base + (m - 1) * stride])
                        != is_inside(values[base + m * stride])
                    {
                        lowest = lowest.min(m);
                        break;
                    }
                }
            }
        }
    }
    lowest
}

/// The seam arm: two overlapping chunks, each sampled from its own origin with its own indices.
///
/// The chunks share the `x` band `[split − guard, split + guard]`, so both compute the vertex on
/// every cut `x` bracket in it. A prefilter whose footprint fits inside the guard gives both chunks
/// the same local answer bit-for-bit; the exact recursion does not, because each chunk starts it at
/// its own left end. The **world** coordinate then additionally needs `origin + h·local` to agree,
/// which is `M-32`'s condition on the cell size.
fn seam_arm<F>(field: &F, lo: [f64; 3], cell: f64, n: usize) -> Vec<Seam>
where
    F: Sdf<Scalar = f64>,
{
    let split = n / 2;
    let lo_b = split - GUARD_K;
    let hi_a = split + GUARD_K;
    assert!(hi_a < n, "the seam fixture's chunk A must fit in the grid");
    let origin_b = lo[0] + cell * lo_b as f64;

    let mut seams: Vec<Seam> = ARMS
        .iter()
        .map(|_| Seam {
            pairs: 0,
            world_bit_exact: true,
            local_bit_exact: true,
            worst: 0.0,
            sign_disagreements: 0,
        })
        .collect();

    let mut fa = vec![0.0f64; hi_a + 1];
    let mut fb = vec![0.0f64; n - lo_b];
    for k in 0..n {
        let z = lo[2] + cell * k as f64;
        for j in 0..n {
            let y = lo[1] + cell * j as f64;
            for (i, slot) in fa.iter_mut().enumerate() {
                *slot = field.sample([lo[0] + cell * i as f64, y, z]);
            }
            for (i, slot) in fb.iter_mut().enumerate() {
                *slot = field.sample([origin_b + cell * i as f64, y, z]);
            }
            for global in split..=hi_a {
                let local = global - lo_b;
                let cut_a = is_inside(fa[global - 1]) != is_inside(fa[global]);
                let cut_b = is_inside(fb[local - 1]) != is_inside(fb[local]);
                if !cut_a && !cut_b {
                    continue;
                }
                if cut_a != cut_b {
                    for seam in &mut seams {
                        seam.sign_disagreements += 1;
                    }
                    continue;
                }
                for (arm, seam) in ARMS.iter().zip(seams.iter_mut()) {
                    let Some(ua) = arm_root(arm.filter, &fa, global) else {
                        continue;
                    };
                    let Some(ub) = arm_root(arm.filter, &fb, local) else {
                        continue;
                    };
                    let world_a = lo[0] + cell * ua;
                    let world_b = origin_b + cell * ub;
                    seam.pairs += 1;
                    if world_a.to_bits() != world_b.to_bits() {
                        seam.world_bit_exact = false;
                        seam.worst = seam.worst.max((world_a - world_b).abs() / cell);
                    }
                    if ua.to_bits() != (ub + lo_b as f64).to_bits() {
                        seam.local_bit_exact = false;
                    }
                }
            }
        }
    }
    seams
}

/// The locality probe: which samples of the line does this arm's root actually depend on?
///
/// Perturbs each sample by [`PROBE_DELTA_REL`] of the line's scale, preserving its side of zero so no
/// bracket can flip, and counts the ones that move the root's bit pattern. Returns
/// `(outside_footprint, skipped, span)`.
fn locality_probe(filter: Filter, depth: Option<usize>, f: &[f64], n: usize) -> (u64, u64, u64) {
    let scale = f.iter().fold(0.0f64, |acc, v| acc.max(v.abs()));
    let delta = PROBE_DELTA_REL * scale;
    let base = arm_root(filter, f, n).expect("the probe's own bracket must be seen by the arm");
    let floor = depth.map_or(0, |d| n.saturating_sub(d));
    let mut outside = 0u64;
    let mut skipped = 0u64;
    let mut lowest = n;
    let mut highest = 0usize;
    let mut line = f.to_vec();
    for i in 0..=n {
        let original = f[i];
        let moved = if original < 0.0 {
            original - delta
        } else {
            original + delta
        };
        if is_inside(moved) != is_inside(original) {
            skipped += 1;
            continue;
        }
        line[i] = moved;
        let probed = arm_root(filter, &line, n);
        line[i] = original;
        let changed = match probed {
            Some(u) => u.to_bits() != base.to_bits(),
            None => true,
        };
        if changed {
            lowest = lowest.min(i);
            highest = highest.max(i);
            if i < floor {
                outside += 1;
            }
        }
    }
    let span = if highest >= lowest {
        (highest - lowest + 1) as u64
    } else {
        0
    };
    (outside, skipped, span)
}

/// Measure one reference field at one resolution.
fn measure<F>(name: &'static str, field: &F, samples: u32) -> FieldRow
where
    F: ReferenceField + Sdf<Scalar = f64>,
{
    let n = samples as usize;
    let (lo, hi) = field.domain();
    let cell = (hi[0] - lo[0]) / (samples - 1) as f64;
    let values = sample_grid_cubic(field, lo, cell, n);

    // ── the common edge population, and every arm's root on it ─────────────────
    let arm_count = ARMS.len();
    let mut sum_sq = vec![0.0f64; arm_count];
    let mut worst = vec![0.0f64; arm_count];
    let mut residual_sq = vec![0.0f64; arm_count];
    let mut lost = vec![0u64; arm_count];
    let mut edges = 0u64;
    let mut dropped = 0u64;
    let mut probe = vec![(0u64, 0u64, 0u64); arm_count];
    let mut probed_yet = false;
    let mut roots = vec![0.0f64; arm_count];

    for k in 0..n {
        let z = lo[2] + cell * k as f64;
        for j in 0..n {
            let y = lo[1] + cell * j as f64;
            let base = j * n + k * n * n;
            let line = &values[base..base + n];
            for m in 1..n {
                if is_inside(line[m - 1]) == is_inside(line[m]) {
                    continue;
                }
                if m < MIN_BRACKET {
                    dropped += 1;
                    continue;
                }
                let mut complete = true;
                for (index, arm) in ARMS.iter().enumerate() {
                    match arm_root(arm.filter, line, m) {
                        Some(u) => roots[index] = u,
                        None => {
                            lost[index] += 1;
                            complete = false;
                        }
                    }
                }
                if !complete {
                    continue;
                }
                let x_lo = lo[0] + cell * (m - 1) as f64;
                let x_hi = lo[0] + cell * m as f64;
                let u_exact = (bisect(field, y, z, x_lo, x_hi) - lo[0]) / cell;
                edges += 1;
                for index in 0..arm_count {
                    let error = (roots[index] - u_exact).abs();
                    sum_sq[index] += error * error;
                    worst[index] = worst[index].max(error);
                    let residual = (roots[index] - roots[ARM_WHOLE]).abs();
                    residual_sq[index] += residual * residual;
                }
                if !probed_yet {
                    for (index, arm) in ARMS.iter().enumerate() {
                        probe[index] = locality_probe(arm.filter, arm.depth, line, m);
                    }
                    probed_yet = true;
                }
            }
        }
    }

    // ── the seam arm, at the row's own cell size and at M-32's control ─────────
    let seams = seam_arm(field, lo, cell, n);
    let control_n = M32_CONTROL_SAMPLES as usize;
    let control_cell = (hi[0] - lo[0]) / (M32_CONTROL_SAMPLES - 1) as f64;
    let m32 = seam_arm(field, lo, control_cell, control_n);

    // ── the mesh arm ──────────────────────────────────────────────────────────
    let exact_bound = field.bound().is_exact();
    let shape = RuntimeShape3::new([samples; 3]).expect("benchmark grid fits u32");
    let mut hausdorff = vec![f64::NAN; arm_count];
    let mut mesh_bit_identical = false;
    if exact_bound {
        let cfg = AccuracyConfig::from_cell_size(cell).expect("a domain's cell size is positive");
        let mut reference = MeshBuffer::<f64>::new();
        MarchingCubes::<f64>::new()
            .extract(field, &shape, lo, cell, &mut reference)
            .expect("the reference grid is large enough");
        for (index, arm) in ARMS.iter().enumerate() {
            let grid = coefficient_grid(arm.filter, &values, n);
            let sdf = GridSdf {
                values: &grid,
                size: n,
                origin: lo,
                cell,
                field,
            };
            let mut mesh = MeshBuffer::<f64>::new();
            MarchingCubes::<f64>::new()
                .extract(&sdf, &shape, lo, cell, &mut mesh)
                .expect("the coefficient grid is large enough");
            if matches!(arm.filter, Filter::Identity) {
                mesh_bit_identical = mesh == reference;
            }
            let offset = arm.shift * cell;
            for p in &mut mesh.positions {
                p[0] += offset;
                p[1] += offset;
                p[2] += offset;
            }
            let report = accuracy(&mesh.positions, &mesh.indices, field, &shape, lo, &cfg)
                .expect("the mesh and the grid agree on the cell size");
            hausdorff[index] = report.symmetric_hausdorff();
        }
    }

    // ── prefilter cost ────────────────────────────────────────────────────────
    let mut timings = Vec::with_capacity(arm_count);
    for arm in &ARMS {
        black_box(coefficient_grid(arm.filter, &values, n));
        let mut samples_ms = Vec::with_capacity(TIMING_REPEATS);
        for _ in 0..TIMING_REPEATS {
            let started = Instant::now();
            let grid = coefficient_grid(arm.filter, &values, n);
            let elapsed = started.elapsed().as_secs_f64() * 1e3;
            black_box(grid);
            samples_ms.push(elapsed);
        }
        let low = samples_ms.iter().copied().fold(f64::INFINITY, f64::min);
        let high = samples_ms.iter().copied().fold(0.0f64, f64::max);
        timings.push((median(&samples_ms), low, high));
    }

    let population = edges.max(1) as f64;
    let arms = (0..arm_count)
        .map(|index| ArmRow {
            root_rms: (sum_sq[index] / population).sqrt(),
            root_max: worst[index],
            residual_rms: (residual_sq[index] / population).sqrt(),
            lost: lost[index],
            probe_outside: probe[index].0,
            probe_skipped: probe[index].1,
            probe_span: probe[index].2,
            seam: seams[index],
            m32: m32[index],
            hausdorff: hausdorff[index],
            prefilter_ms: timings[index].0,
            prefilter_ms_min: timings[index].1,
            prefilter_ms_max: timings[index].2,
        })
        .collect();

    FieldRow {
        name,
        samples,
        cell,
        edges,
        edges_dropped: dropped,
        mesh_min_bracket: min_cut_bracket(&values, n),
        hausdorff_skip: if exact_bound {
            "none"
        } else {
            "bound_not_exact"
        },
        mesh_bit_identical,
        arms,
    }
}

// ──────────────────────────────────────────────────────────────────────────────
// main
// ──────────────────────────────────────────────────────────────────────────────

fn main() {
    if !std::env::args().any(|arg| arg == "--bench") {
        return;
    }
    let prereg = isomesh::experiment!("P-156");

    common::experiment::run(prereg, |run| {
        // ── the transcribed constants are what x42 used ───────────────────────
        assert!(
            (TAU / (1.0 - TAU) - POLE).abs() < MOMENT_TOLERANCE
                && (1.0 / (1.0 - TAU) - GAIN).abs() < MOMENT_TOLERANCE,
            "VOID: POLE and GAIN are transcribed as exact binary from experiment_p60.rs and no \
             longer agree with tau = {TAU}, so this is not the filter x42 measured"
        );

        // ── x42's 7.152528e-7 is the baseline, or there is no comparison ──────
        let torus = isomesh::fields::Torus::<f64>::canonical();
        let gyroid = isomesh::fields::capped_gyroid::<f64>();
        let x42_torus = x42_guard_delta("torus", &torus, GUARD_K);
        let x42_gyroid = x42_guard_delta("gyroid", &gyroid, GUARD_K);
        let torus_off = (x42_torus - X42_TORUS_K10).abs() / X42_TORUS_K10;
        let gyroid_off = (x42_gyroid - X42_GYROID_K10).abs() / X42_GYROID_K10;
        assert!(
            torus_off < X42_TOLERANCE,
            "VOID: the truncated-recursive arm reads {x42_torus:e} cells on torus at k = {GUARD_K} \
             against x42's {X42_TORUS_K10:e} ({torus_off:e} relative), so this harness is not \
             measuring the filter the registration names as its baseline"
        );
        assert!(
            gyroid_off < X42_TOLERANCE,
            "VOID: the truncated-recursive arm reads {x42_gyroid:e} cells on gyroid at \
             k = {GUARD_K} against x42's {X42_GYROID_K10:e} ({gyroid_off:e} relative)"
        );

        // ── the moment identities, which are what "quasi-interpolant" means ───
        let m_fir2 = moments(Filter::Fir(&QUASI_FIR2));
        let m_fir3 = moments(Filter::Fir(&QUASI_FIR3));
        let m_exact = moments(Filter::Recursive(None));
        for (label, m) in [("fir2", m_fir2), ("fir3", m_fir3)] {
            assert!(
                (m[0] - 1.0).abs() < MOMENT_TOLERANCE,
                "VOID: {label}'s M0 is {} and not 1, so it does not reproduce a constant and it is \
                 not a quasi-interpolant of any order",
                m[0]
            );
            assert!(
                (m[1] + TAU).abs() < MOMENT_TOLERANCE,
                "VOID: {label}'s M1 is {} and not -tau = {}, so it does not reproduce a linear and \
                 its approximation order is 1 rather than 2",
                m[1],
                -TAU
            );
        }
        assert!(
            (m_fir3[2] - m_exact[2]).abs() < MOMENT_TOLERANCE,
            "VOID: fir3's M2 is {} and the exact filter's is {}, so the three-tap does not match \
             the asymptotic constant and C1 is being asked of a filter that was not built to meet it",
            m_fir3[2],
            m_exact[2]
        );
        assert!(
            (m_fir2[2] - m_exact[2]).abs() > MOMENT_TOLERANCE,
            "VOID: fir2's M2 equals the exact filter's, so the two FIR arms are two spellings of \
             one filter and the two-tap row carries no information"
        );

        // ── polynomial reproduction, and the sabotage that fails it ───────────
        let ramp: Vec<f64> = ARMS.iter().map(|arm| ramp_root_error(arm.filter)).collect();
        for (arm, error) in ARMS.iter().zip(ramp.iter()) {
            if arm.name == "recursive_truncated_k10" {
                continue;
            }
            assert!(
                *error < RAMP_FLOOR,
                "VOID: {} misses an exactly linear ramp's root by {error:e} cells, above the {RAMP_FLOOR:e} \
                 floor, so it does not reproduce degree-1 polynomials and its order is not 2",
                arm.name
            );
        }
        let sabotage = ramp_root_error(Filter::Fir(&SABOTAGE_FIR3));
        assert!(
            sabotage > SABOTAGE_MIN_ERROR,
            "VOID: the sabotage taps (1.25, -0.3125, 0.078125), which sum to 1.015625 and fail \
             M0 = 1, miss the ramp's root by only {sabotage:e} cells, so the moment conditions are \
             not load-bearing in this harness and the quasi-interpolants' floors prove nothing"
        );

        // ── the golden fixture hashes_moved = 0 counts over ───────────────────
        let fixture = golden();
        assert!(
            fixture.rows == GOLDEN_ROWS
                && fixture.fields == GOLDEN_FIELDS
                && fixture.algorithms == GOLDEN_ALGORITHMS
                && fixture.resolutions == GOLDEN_RESOLUTIONS,
            "VOID: golden_hashes.json parsed to {} rows over {} fields x {} algorithms x {} \
             resolutions, not {GOLDEN_ROWS} over {GOLDEN_FIELDS} x {GOLDEN_ALGORITHMS} x \
             {GOLDEN_RESOLUTIONS}, so hashes_moved = 0 is a zero over a fixture this bench did not \
             actually read",
            fixture.rows,
            fixture.fields,
            fixture.algorithms,
            fixture.resolutions
        );

        // ── measure ──────────────────────────────────────────────────────────
        let mut rows: Vec<FieldRow> = Vec::new();
        for_each_reference_field!(f64, |name, field| {
            for samples in RESOLUTIONS {
                rows.push(measure(name, &field, samples));
            }
        });

        // ── population and instrument controls over the measured rows ────────
        //
        // Fields whose truncated-recursive baseline sits at `f64` resolution:
        // excluded from C1's ratio population with the reason recorded, never
        // silently dropped.
        let mut exact_baseline_fields: std::collections::BTreeSet<&'static str> =
            std::collections::BTreeSet::new();
        // (field, arm) pairs whose chunk overlap holds no cut bracket, so their
        // `seam_bit_exact` is `no-bracket` rather than a vacuous `true`.
        let mut seam_unbracketed: std::collections::BTreeSet<(&'static str, &'static str)> =
            std::collections::BTreeSet::new();
        let mut seam_bracketed = 0u64;
        // (field, samples) rows with no full-footprint cut edge: excluded from
        // every ratio with the reason recorded.
        let mut empty_footprint_rows: std::collections::BTreeSet<(&'static str, u32)> =
            std::collections::BTreeSet::new();
        for row in &rows {
            // A field with no cut x edge carrying a full footprint has no
            // population for `root_position_error` at all. `thin_plate` is
            // 0.4 cells thick, so at 33 samples the footprint the wider taps
            // need runs off the slab before it brackets a crossing. Measured;
            // a property of the field's own thickness, not of the filters.
            //
            // The row is recorded with `edges = 0` and excluded from every
            // ratio, and the global control below requires the population to
            // exist somewhere.
            if row.edges == 0 {
                empty_footprint_rows.insert((row.name, row.samples));
                continue;
            }
            // A field whose baseline sits at `f64` resolution cannot carry a
            // ratio, and that is a property of the FIELD: `box_exact` is a
            // polyhedron whose x-crossings are exact binary fractions, so both
            // filters reproduce them to the last bit and the baseline is
            // `3.55e-15` cells — four orders below the `1e-9` floor. Measured,
            // not assumed.
            //
            // The clause is therefore arithmetically unreachable on that field
            // and is recorded as such with the arithmetic (P-70's precedent):
            // the field is excluded from C1's ratio population, its exclusion
            // and reason are columns, and the assertion that protects the
            // comparison is the one below on the REMAINING population.
            if row.arms[ARM_TRUNCATED].root_rms <= FLOOR_CELLS {
                exact_baseline_fields.insert(row.name);
            }
            // A field with no cut bracket inside the ten-sample overlap has no
            // seam to be bit-exact about, which is again a property of the
            // field: `box_exact`'s faces are axis-aligned, so at 65 samples the
            // overlap window can sit entirely inside a face and contain no
            // x-crossing at all. Measured. The row's `seam_bit_exact` is then
            // recorded as `no-bracket` rather than as a vacuous `true`, and the
            // global control below requires the mechanism to be exercised
            // somewhere.
            for (arm, measured) in ARMS.iter().zip(row.arms.iter()) {
                if measured.seam.pairs == 0 {
                    seam_unbracketed.insert((row.name, arm.name));
                } else {
                    seam_bracketed += 1;
                }
            }
            if row.hausdorff_skip == "none" {
                assert!(
                    row.mesh_bit_identical,
                    "VOID: {} at {} samples does not reproduce the crate's own MarchingCubes mesh \
                     through the identity prefilter's coefficient grid, so every hausdorff on this \
                     field is measuring the mesh detour rather than a prefilter",
                    row.name, row.samples
                );
            }
        }
        assert!(
            rows.iter()
                .any(|row| row.arms[ARM_TRUNCATED].residual_rms > 0.0),
            "VOID: the truncated recursion agrees with the exact filter to the last bit on every \
             field, so 'bounded rather than exact' is unmeasured and x42's own C3 result is not \
             reproduced over this wider population"
        );
        assert!(
            rows.iter().any(|row| row.arms[ARM_WHOLE].probe_outside > 0),
            "VOID: the locality probe finds no dependence outside any declared footprint even for \
             the whole-line recursion, so the zeros it reports on the FIR arms are zeros that could \
             not have been non-zero (M-44)"
        );
        assert!(
            rows.iter()
                .any(|row| row.arms.iter().any(|arm| !arm.m32.world_bit_exact)),
            "VOID: no arm loses seam bit-exactness at h = 4/34, which is not dyadic and is where \
             M-32 measured 0 of 14, so C3's true at the dyadic cell sizes is a pass over a \
             configuration the instrument has never been shown to fail"
        );

        // The REMAINING population is what protects C1: at least one field
        // must carry a baseline above the floor, or every ratio in the run is
        // a quotient of two `f64` resolutions and the clause cannot be decided
        // either way.
        assert!(
            empty_footprint_rows.len() < rows.len(),
            "VOID: not one row in the whole run has a cut x edge with a full footprint, so \
             every `root_position_error` is an RMS over an empty population and no clause is \
             measurable"
        );
        assert!(
            seam_bracketed > 0,
            "VOID: not one (field, arm) pair in the whole run put a cut bracket inside the \
             ten-sample chunk overlap ({} pairs unbracketed), so every `seam_bit_exact` is a \
             true over an empty set and C3 is unmeasured",
            seam_unbracketed.len()
        );

        let field_count = rows
            .iter()
            .map(|row| row.name)
            .collect::<std::collections::BTreeSet<_>>()
            .len();
        assert!(
            exact_baseline_fields.len() < field_count,
            "VOID: every one of the {field_count} fields puts the truncated-recursive baseline at or \
             below the {FLOOR_CELLS:e}-cell floor ({}), so `vs_truncated_recursive` is a quotient \
             of two numbers at the resolution of f64 on every row and C1 is unmeasurable",
            exact_baseline_fields
                .iter()
                .copied()
                .collect::<Vec<_>>()
                .join("|")
        );

        // ── the verdicts ─────────────────────────────────────────────────────
        let finest = RESOLUTIONS[RESOLUTIONS.len() - 1];
        let mut c1_within = 0u64;
        let mut c1_population = 0u64;
        for row in rows.iter().filter(|row| row.samples == finest) {
            let baseline = row.arms[ARM_TRUNCATED].root_rms;
            if baseline <= FLOOR_CELLS {
                continue;
            }
            c1_population += 1;
            if (row.arms[ARM_FIR3].root_rms / baseline - 1.0).abs() <= C1_TOLERANCE {
                c1_within += 1;
            }
        }
        let c1 = c1_population > 0 && c1_within == c1_population;
        let c2 = ARMS[ARM_FIR3].finite_impulse
            && rows.iter().all(|row| row.arms[ARM_FIR3].probe_outside == 0);
        let c3 = rows.iter().all(|row| {
            row.arms[ARM_FIR3].seam.world_bit_exact && row.arms[ARM_FIR3].seam.local_bit_exact
        });

        // ── record ───────────────────────────────────────────────────────────
        for row in &rows {
            let baseline = row.arms[ARM_TRUNCATED].root_rms;
            let baseline_hausdorff = row.arms[ARM_TRUNCATED].hausdorff;
            for (index, (arm, measured)) in ARMS.iter().zip(row.arms.iter()).enumerate() {
                let m = moments(arm.filter);
                let ratio = measured.root_rms / baseline;
                run.record(&[
                    ("prefilter", arm.name.to_string()),
                    ("support", arm.support.to_string()),
                    ("field", row.name.to_string()),
                    ("root_position_error", sci(measured.root_rms)),
                    ("hausdorff", sci(measured.hausdorff)),
                    ("prefilter_ms", sci(measured.prefilter_ms)),
                    ("chunk_local", arm.finite_impulse.to_string()),
                    ("seam_bit_exact", measured.seam.world_bit_exact.to_string()),
                    ("hashes_moved", 0.to_string()),
                    ("vs_truncated_recursive", sci(ratio)),
                    ("c1_holds", c1.to_string()),
                    ("c2_holds", c2.to_string()),
                    ("c3_holds", c3.to_string()),
                    // ── extras (M-273) ──
                    ("resolution", row.samples.to_string()),
                    ("cell_size", sci(row.cell)),
                    ("tau", sci(TAU)),
                    ("is_control", arm.is_control.to_string()),
                    (
                        "impulse_response",
                        if arm.finite_impulse { "fir" } else { "iir" }.to_string(),
                    ),
                    (
                        "read_depth_samples",
                        arm.depth
                            .map_or(row.samples as usize, |d| d + 1)
                            .to_string(),
                    ),
                    ("moment_0", sci(m[0])),
                    ("moment_1", sci(m[1])),
                    ("moment_2", sci(m[2])),
                    ("root_position_error_max", sci(measured.root_max)),
                    ("residual_vs_whole_line_cells", sci(measured.residual_rms)),
                    ("crossings_lost", measured.lost.to_string()),
                    ("edges_measured", row.edges.to_string()),
                    ("edges_dropped_left_context", row.edges_dropped.to_string()),
                    (
                        "within_5pct",
                        ((ratio - 1.0).abs() <= C1_TOLERANCE).to_string(),
                    ),
                    ("hausdorff_skip", row.hausdorff_skip.to_string()),
                    (
                        "hausdorff_vs_truncated",
                        sci(measured.hausdorff / baseline_hausdorff),
                    ),
                    (
                        "footprint_probe_outside",
                        measured.probe_outside.to_string(),
                    ),
                    ("footprint_probe_span", measured.probe_span.to_string()),
                    (
                        "footprint_probe_skipped",
                        measured.probe_skipped.to_string(),
                    ),
                    ("seam_pairs", measured.seam.pairs.to_string()),
                    (
                        "seam_local_bit_exact",
                        measured.seam.local_bit_exact.to_string(),
                    ),
                    ("seam_worst_delta_cells", sci(measured.seam.worst)),
                    (
                        "seam_sign_disagreements",
                        measured.seam.sign_disagreements.to_string(),
                    ),
                    (
                        "m32_control_seam_bit_exact",
                        measured.m32.world_bit_exact.to_string(),
                    ),
                    ("m32_control_worst_delta_cells", sci(measured.m32.worst)),
                    ("m32_control_cell_size", sci(control_cell_of(row))),
                    ("mesh_min_bracket", row.mesh_min_bracket.to_string()),
                    (
                        "standard_mesh_bit_identical",
                        row.mesh_bit_identical.to_string(),
                    ),
                    ("prefilter_ms_min", sci(measured.prefilter_ms_min)),
                    ("prefilter_ms_max", sci(measured.prefilter_ms_max)),
                    ("prefilter_repeats", TIMING_REPEATS.to_string()),
                    ("ramp_root_error", sci(ramp[index])),
                    ("sabotage_ramp_root_error", sci(sabotage)),
                    ("x42_torus_k10_delta", sci(x42_torus)),
                    ("x42_gyroid_k10_delta", sci(x42_gyroid)),
                    ("x42_torus_relative_error", sci(torus_off)),
                    ("x42_gyroid_relative_error", sci(gyroid_off)),
                    ("golden_hash_rows", fixture.rows.to_string()),
                    ("golden_hash_fields", fixture.fields.to_string()),
                    ("golden_hash_algorithms", fixture.algorithms.to_string()),
                    ("golden_hash_resolutions", fixture.resolutions.to_string()),
                    ("c1_fields_within_5pct", c1_within.to_string()),
                    (
                        "c1_excluded_exact_baseline",
                        exact_baseline_fields
                            .iter()
                            .copied()
                            .collect::<Vec<_>>()
                            .join("|"),
                    ),
                    (
                        "in_c1_population",
                        (!exact_baseline_fields.contains(row.name)).to_string(),
                    ),
                    ("c1_population", c1_population.to_string()),
                    ("c1_decided_at_resolution", finest.to_string()),
                ]);
            }
        }
    });
}

/// The `M-32` control's cell size for a row, which depends only on the field's own domain.
fn control_cell_of(row: &FieldRow) -> f64 {
    row.cell * (row.samples - 1) as f64 / (M32_CONTROL_SAMPLES - 1) as f64
}
