//! **P-76 — triplanar times stochastic filtering, and what a destructible world
//! has left to pay it with.**
//!
//! Ticket: R-076. Pre-registered before this harness existed.
//!
//! ```bash
//! cargo bench --bench experiment_p76
//! ```
//!
//! Writes `docs/experiments/p-76.csv`.
//!
//! # Hypothesis, as registered
//!
//! Triplanar and stochastic filtering compose multiplicatively and nobody states
//! the product. Three planes times three stochastic taps is nine fetches per
//! map; `bevy_isomesh/examples/triplanar.wgsl` already pays the three.
//! Stochastic Texture Filtering (Pharr, Wronski, Kettunen, Bartels, Wolfe,
//! Bako, Salvi & Hasselgren, *"Filtering After Shading With Stochastic Texture
//! Filtering"*, `arXiv:2305.05810`) and Heitz & Neyret's histogram-preserving
//! blending (*"High-Performance By-Example Noise using a Histogram-Preserving
//! Blending Operator"*, `10.1145/3233304`, **"over 20x faster"** than the
//! procedural-noise state of the art — **the hardware is not named in that
//! claim and this harness flags it**) both trade fetches for temporal
//! accumulation. The conflict is the finding: **a destructible world has already
//! spent its temporal budget, because geometry that changes rejects history.**
//!
//! - **C1.** Selecting one triplanar plane per pixel stochastically drops the
//!   fetch count by **exactly 3x** and the fragment cost by at least **2x** at
//!   1920x1080.
//! - **C2.** With TAA resolving the stochastic choice, the ghosting cost under an
//!   active dig is worse than the fetch saving is worth: mean absolute error
//!   against a 3-plane reference in the **8 frames following an edit**, above the
//!   same scene with no digging.
//! - **C3.** Biplanar mapping — two planes, no stochastic term, no temporal debt
//!   — gets at least **half** the saving of C1 with none of C2's cost.
//!
//! VACUITY CONTROL, as registered: **the dig arm must produce a non-zero
//! history-rejection count from `P-77`'s instrument, or C2 cannot fire.**
//! Asserted below, together with its `M-44` other half: a configuration in which
//! that count really is zero, so the non-zero is a measurement rather than a
//! fixture that could not have reported a zero.
//!
//! # The SHARE line, recomputed before a line of this harness was written
//!
//! **C1's fetch half is exact arithmetic and it is reachable — but the
//! registration's own count of the baseline is wrong by 3x, and that is the
//! first finding here.** `triplanar.wgsl` does not pay three fetches per map.
//! `game_dig` builds its terrain material with `settings.z = LAYER_BLEND = -1.0`
//! (`game_dig.rs:577`, `game_dig.rs:998`), which takes the shader's **else**
//! branch (`triplanar.wgsl:142-149`): `layer_ar` is called once per **stratum**
//! and samples all **three planes** inside (`triplanar.wgsl:87-92`), and
//! `layers_na` samples all **three strata** for each of the three planes
//! (`triplanar.wgsl:97-101`). So the terrain fragment issues
//!
//! ```text
//! albedo/roughness : 3 strata x 3 planes = 9
//! normal/AO        : 3 planes x 3 strata = 9
//!                                   total 18 texture fetches per fragment
//! ```
//!
//! **The nine-per-map the registration attributes to triplanar-times-stochastic
//! is already paid by triplanar alone, times the stratum blend.** The product
//! nobody states is not `3 x 3` but `3 planes x 3 strata x N taps`, and at the
//! registration's own N = 3 that is **54 fetches per fragment**, not nine.
//! Measured as [`ARM_STF_3TAP`] below.
//!
//! The 3x itself survives the correction intact, because the plane count factors
//! out of both maps: `2 x 3 planes x 3 strata` over `2 x 1 plane x 3 strata` is
//! `18/6 = 3` **exactly**, in integers, on any machine. That equality is
//! asserted rather than measured.
//!
//! **C1's cost half is reachable if and only if the stochastic arm's own fetch
//! cost is at least the fixed cost, and this harness measures both.** Writing `F`
//! for the per-fragment cost that does not depend on the plane count (the
//! projections, the two `pow` weight sets, the whiteout blend, the lighting) and
//! `S` for the cost of **one** fetch,
//!
//! ```text
//! t(k fetches) = F + k*S      ratio = t(18)/t(6) = (F + 18S)/(F + 6S)
//! ratio >= 2  <=>  F + 18S >= 2F + 12S  <=>  6S >= F
//! ```
//!
//! so the clause is **exactly the claim that the six fetches the stochastic arm
//! still pays cost at least as much as everything else in the fragment put
//! together**, the ceiling is 3 as `F -> 0`, and there is no arithmetic
//! obstruction in either direction. `F` and `S` are recovered here by a
//! least-squares fit of `fragment_ms` against the fetch count over the four
//! bilinear arms at `k = 2, 6, 12, 18` (`fit_intercept_ms`,
//! `fit_slope_ms_per_fetch`, `fit_r2`), so a reader who does not accept this
//! harness's `F` can recompute the ratio for their own.
//!
//! **The fit is an explanation and not the verdict, and the difference between
//! them is itself a result.** Dropping two planes drops more than twelve fetches:
//! it also drops two of the three whiteout blends and two of the three
//! `layers_na` accumulations, so the stochastic arm's `F` is *smaller* than the
//! baseline's and the measured ratio exceeds the pure-fetch prediction. Both
//! numbers are columns (`c1_cost_ratio`, `fit_predicted_ratio_18_over_6`) and the
//! gap between them is the plane-count-dependent ALU the linear model cannot see.
//!
//! **C3's "at least half" is arithmetically forced to be exactly half, and this
//! is the second finding.** Biplanar drops one plane of three, stochastic
//! selection drops two:
//!
//! ```text
//! biplanar saving    18 - 12 = 6 fetches
//! stochastic saving  18 -  6 = 12 fetches
//! ratio              6/12 = 0.5 exactly
//! ```
//!
//! and in *any* cost model linear in the fetch count the time saving is the same
//! `S/(2S) = 0.5`. So C3's threshold is met by **equality**, it holds only
//! because the registered inequality is non-strict, and in wall clock it is
//! decided by whether the residual of the linear fit happens to fall above or
//! below zero — i.e. by noise, not by mechanism. Both halves are reported
//! (`c3_saving_ratio_fetches` is exact, `c3_saving_ratio_ms` is not) and the
//! verdict is stated with that caveat rather than without it.
//!
//! **C2's registered comparison is reachable and its magnitude is the whole
//! question.** `P-77` measured this instrument's steady-state history-rejection
//! rate at **0.0049** in exactly the arm reused here, so the rejected population
//! exists and is small; `P-77` also found the edit-frame rate at **0.0142**, a
//! ratio of 2.9 against its registered 5x, and **falsified its own C1**. The
//! registration for this row says that "no difference between digging and
//! static" *"would mean history rejection is not the bottleneck and stochastic
//! filtering is free here after all"*, and `P-77`'s registration says the
//! symmetric thing — that a single-frame spike closes `R-077` and `P-76`'s C2
//! together. **The linkage is stated explicitly in the report:** `P-77`'s C1 is
//! already falsified on the ratio, so C2 here cannot inherit a 5x temporal
//! catastrophe. What it can still find, and what it does find, is a difference
//! whose *size* is the recommendation.
//!
//! **C2 is also not scoreable as literally registered, and that is said before
//! the run.** *"Worse than the fetch saving is worth"* is a comparison between a
//! mean absolute error and a millisecond, and there is no exchange rate between
//! them anywhere in this repository. The registration's own `falsified_by` names
//! the operative test — *"no difference between digging and static"* — so
//! `c2_holds` encodes exactly that: `mae_vs_reference_digging >
//! mae_vs_reference_static`. The magnitudes, the affected pixel population and
//! the TAA error floor are all reported so a reader can apply their own exchange
//! rate.
//!
//! # `frame_ms` is a registered column and there is no renderer
//!
//! `crates/isomesh` must not depend on Bevy (`CLAUDE.md` hard rule 2), so there
//! is no engine here and **no frame time is invented**. Both time columns are
//! **per-pixel material shading cost over a real visibility buffer**, which is
//! the same choice `P-75` made and labelled (`experiment_p75.rs:78-84`: *"They
//! are shading cost. They are not frame time"*):
//!
//! - **`fragment_ms`** is the median wall clock of **one single-threaded pass of
//!   the material shading stage over the visible fragments of a 1920x1080
//!   visibility buffer**, in the pose `game_dig` opens in. Ray casting, normal
//!   evaluation and the random-seed derivation are outside it; the stage begins
//!   at the world position and normal and ends at a shaded RGB.
//! - **`frame_ms`** is that same cost per fragment scaled to **full coverage**,
//!   `fragment_ns_per_fragment * 1920*1080 / 1e6`: the material shading stage's
//!   cost for a 1080p frame in which every pixel is terrain. It is a shading
//!   stage cost with a stated denominator, **not an engine frame time**, and
//!   nothing else is in it — no geometry pass, no lighting pass, no post.
//!
//! No registered clause depends on a frame time: C1's second half is about *"the
//! fragment cost"*, which is `fragment_ms`. So nothing here is scored VACUOUS on
//! `frame_ms`'s account; the column is reported with the definition above.
//!
//! **The one thing these clocks cannot be is a GPU ratio.** `F` here is P-77's
//! geometric lighting plus a Blinn-Phong term, not `apply_pbr_lighting`; Bevy's
//! PBR is substantially more ALU than that, so the real `F` is **larger** and the
//! real ratio therefore **smaller** than the one measured here. And this
//! harness's 8 MiB of texture arrays fit in this machine's L3, where a GPU's
//! texture cache is kilobytes — so `S` here is on the cheap side of hardware too.
//! `fetches_per_fragment` and `texels_per_fragment` are the machine-independent
//! quantities and they are what C1's first half and C3's exact half are gated on
//! (`M-282`: gate a performance claim on a count, not a wall clock).
//!
//! # The instrument, in enough detail to judge it
//!
//! ## The material shading stage (C1, C3)
//!
//! A transcription of `bevy_isomesh/examples/triplanar.wgsl`'s `LAYER_BLEND`
//! path, fetch for fetch:
//!
//! | quantity | value | source |
//! |---|---|---|
//! | world units per tile | 1.5 | `game_dig.rs:992` (`settings.x`) |
//! | blend sharpness | 4.0 | `game_dig.rs:992` (`settings.y`) |
//! | forced layer, terrain | `LAYER_BLEND = -1.0` | `game_dig.rs:577`, `game_dig.rs:998` |
//! | forced layer, walls | `LAYER_CONCRETE = 3.0` | `game_dig.rs:579`, `game_dig.rs:999` |
//! | array layers | 4 | `game_dig.rs:575` |
//! | plane weights | `pow(abs(n), s) / sum` | `triplanar.wgsl:64-67` |
//! | plane UVs | `p.zy`, `p.xz`, `p.xy` | `triplanar.wgsl:113-115` |
//! | slope ramp | `smoothstep(0.55, 0.82, n.y)` | `triplanar.wgsl:76-77, 122` |
//! | depth ramp | `smoothstep(-1.6, -0.4, p.y)` | `triplanar.wgsl:83-84, 123` |
//! | stratum weights | `up*shallow, (1-up)*shallow, 1-shallow` | `triplanar.wgsl:124` |
//! | roughness clamp | `[0.05, 1.0]` | `triplanar.wgsl:156` |
//! | whiteout normal blend | per-plane swizzle | `triplanar.wgsl:167-172` |
//!
//! **The texture content is synthesised and the layout is not.** `game_dig`
//! compiles in two `512x2048` PNGs reinterpreted as four `512`-square layers
//! (`game_dig.rs:748-755`, `game_dig.rs:786-797`), and this crate cannot decode
//! a PNG without a new dev-dependency — `P-77` hit the same wall
//! (`experiment_p77.rs:312-315`). So the arrays here are **the same shape, the
//! same element type and the same working set**: two `512 x 512 x 4` RGBA8
//! arrays, 4 MiB each, filled with fractal value noise. What the *cost*
//! instrument depends on is the footprint and the access pattern, both of which
//! are the demo's; what the *error* instrument depends on is high-frequency
//! content, which is why the noise is tuned to the demo's `1.5`-unit tile.
//!
//! **A fetch is bilinear except where an arm is stochastic, and the texel count
//! is reported separately, because this is where the registration's arithmetic
//! goes wrong a second time.** One `textureSample` of a filtered array is
//! **four** texel loads. `arXiv:2305.05810`'s whole proposal is to replace that
//! filter with point samples drawn from the filter weights, which is unbiased
//! and costs **one** texel per tap. So three stochastic taps is three texels
//! where one bilinear fetch was four: **at three taps STF is cheaper in texels
//! than the filtered fetch it replaces, not three times more expensive.** The
//! "nine fetches per map" the registration prices as a cost is nine *fetch
//! instructions* and nine texels against triplanar's three instructions and
//! twelve texels. Both counts are columns.
//!
//! **And this instrument cannot convert that texel count into a time, which is
//! stated rather than papered over.** `bilinear_minus_point_ns_per_texel` — the
//! measured cost of the three extra texels a tap quad reads, from two arms that
//! differ in nothing else — comes out **at or below zero**. The four taps are
//! `(x, y)`, `(x+1, y)`, `(x, y+1)`, `(x+1, y+1)` of a 512-wide RGBA8 layer, so
//! the pairs share cache lines and an out-of-order core prefetches the second
//! line for free; the point path meanwhile pays two `Rng::unit` calls to choose
//! its tap. So on **this** machine the *fetch instruction* count is what costs and
//! the texel count is nearly free, and `stf_3tap` at 54 instructions and 54
//! texels measures **slower** than the baseline at 18 instructions and 72 texels.
//! A GPU's texture unit prices these the other way round — one instruction, four
//! taps, one L1 lookup — so **neither machine's wall clock generalises and only
//! the two counts do.** That is the whole reason both are columns.
//!
//! **The fetch tally is inside the timed pass, one increment per fetch, and that
//! is deliberate.** Moving it out would need a second code path. Because it is
//! paid per fetch it is indistinguishable from a slightly larger `S`, and `S` is
//! exactly the quantity the ratio is a function of — so the fit and the ratio are
//! internally consistent. The bias against a real shader is one integer add per
//! fetch, at most one cycle; `tally_ns_bound_per_fetch` states it and
//! `fit_slope_ns_per_fetch` is the number it is a fraction of.
//!
//! ## The temporal instrument (C2)
//!
//! **`P-77`'s, reused rather than re-derived**, because the registered vacuity
//! control names it: a software TAA resolve over a ray-cast depth buffer of
//! `game_dig`'s own field with exact motion vectors. Transcribed from
//! `crates/isomesh/benches/experiment_p77.rs`, with its constants and its two
//! hard-won fixture corrections:
//!
//! 1. The field is `game_dig::Ground` verbatim, brushes are subtracted spheres
//!    of radius `0.25`, and the depth buffer is `game_dig::trace`
//!    (`game_dig.rs:1373-1394`) with `AIM_NEAR = 0.30`, `AIM_FAR = 25.0`,
//!    `AIM_STEPS = 128`, `AIM_HIT = 0.01`, `LIPSCHITZ = 1.25`.
//! 2. Motion vectors are **exact**: the previous frame's view, including the
//!    previous frame's jitter, applied to the current pixel's world hit point.
//! 3. The resolve is Karis 2014's: a 10-entry R2 jitter sequence, bilinear
//!    history fetch at the reprojected position, rectification by **clipping**
//!    against the 3x3 neighbourhood's YCoCg AABB, `history = lerp(clipped,
//!    current, 0.1)`. A sample is **rejected** when the clip moved it. `P-77`'s
//!    k-DOP arms are not re-run: `P-77` proved analytically and then measured
//!    that a k-DOP containing the AABB's axes recovers **exactly zero**, so
//!    there is nothing here for a second volume to add.
//! 4. **The measured arm is `P-77`'s `HEADLINE_ARM`** — standing still, looking
//!    down at the rock under your feet, digging — for `P-77`'s reason: a walking
//!    camera saturates the rejection rate at 86.6% at 7.2 px of reprojection per
//!    frame, and a ratio between two saturated numbers measures nothing. A
//!    walking arm is run anyway, second, because *"a destructible world has
//!    already spent its temporal budget"* is a claim about locomotion too.
//! 5. **Three streams share one G-buffer per world**, so the comparison is
//!    per-pixel rather than between simulations that have drifted: the 3-plane
//!    reference, the stochastic single-plane arm, and biplanar. Each has its own
//!    history buffer and its own resolve. Ray casting happens **once per world
//!    per frame** and is shared by all three, which is what makes 32 frames of
//!    three streams affordable.
//! 6. **`mae_vs_reference_*` is image error, not sample error**: the mean over
//!    hit pixels and RGB channels of `|resolved - reference|`, where the
//!    reference is the deterministic 3-plane shade of *this* frame. Pixels with
//!    no history are included, because a viewer sees them.
//! 7. The registered `history_rejected_*` columns are counted on the
//!    **deterministic** stream, which makes them `P-77`'s own quantity and
//!    cross-checkable against its `0.0049`. The stochastic stream's rejection
//!    counts are separate columns, and they are much larger for a reason that is
//!    the mechanism itself.
//!
//! # Controls, each an assertion rather than a printed number
//!
//! - **THE REGISTERED VACUITY CONTROL.** `history_rejected_digging > 0`.
//! - **`M-44`'s other half.** [`ZERO_PROOF_ARM`] freezes the camera and turns the
//!   jitter off: the reprojection is the identity, the history converges to the
//!   current frame, and the rejection count collapses. If it did not, the
//!   non-zero above would be a property of the harness.
//! - **Exactness.** `fetches(triplanar) == 3 * fetches(stochastic)` as integers,
//!   and the tallied total equals the analytic per-fragment count times the
//!   fragment count, for every arm — so an arm that silently skipped a fetch is a
//!   failure rather than a fast number.
//! - **The stochastic estimator is unbiased in what it selects.** 256 independent
//!   draws per fragment over a 4096-fragment subsample must converge on the
//!   deterministic 3-plane blend, in the **sampled maps** — see [`Shaded`] for why
//!   the shaded colour is explicitly not part of the contract. Gated on a mean,
//!   because unbiasedness is a statement about a mean.
//! - **The population is non-empty**, in both stages, every frame.
//! - **Changed pixels lie inside the brush silhouette**, `P-77`'s control: it
//!   catches a wrong projection, a wrong brush position or a leaky carve.
//! - **One build, one run, interleaved reps, median, and the clock on the row**
//!   (`M-280`, `M-281`). Roughly three sibling agents are compiling while this
//!   runs; every ratio is between arms measured in the same interleaved sweep.
//!
//! # Three fixture defects this harness's own controls found, in order
//!
//! Recorded because they are the most transferable part of the row.
//!
//! 1. **The registered "unbiased" control was asserting the wrong contract, and
//!    the truth is a finding.** The first version compared 256 draws of the
//!    single-plane arm's *shaded colour* against the 3-plane arm's and demanded
//!    they agree; it failed at `0.02229`. It was right to fail: `normalize`,
//!    `max(0, .)`, `powf` and the roughness clamp are nonlinear, so the estimator
//!    is unbiased in the **maps it samples** and biased in the **colour it
//!    shades**. That is `arXiv:2305.05810`'s title — *filtering **after**
//!    shading* — and it means a share of the stochastic error is a **bias no
//!    accumulation removes**. The control now gates the contract
//!    (`estimator_mean_signed_error_sampled`, `1.9e-5` against a `2e-3`
//!    tolerance) and reports the rest (`shading_bias_mean_signed`).
//! 2. **Dropping the clamp's `eps` made the zero-proof arm reject 31,274
//!    samples.** See [`CLAMP_EPS`]: `from_ycocg(to_ycocg(c))` is the identity in
//!    real arithmetic and not in `f32`, and a `1e-7` drift is clipped to `s ~ 0`
//!    at any pixel that is its own 3x3's extremum. With the paper's `1e-5`
//!    restored the same arm rejects **zero** and `M-44`'s control is a control.
//! 3. **`6S >= F`, not `S >= F`.** The reachability line printed the wrong
//!    inequality on the first run: the ratio is `t(18)/t(6)`, so the fixed cost
//!    competes with the **six** fetches the stochastic arm still pays, not with
//!    one. Per fetch `S = 21.7 ms` is well under `F = 71.5 ms` and the wrong form
//!    said "unreachable" while the measured ratio was 2.33. Corrected before the
//!    reported run; the column is `c1_cost_reachable_6s_ge_f`.

// `needless_range_loop` and `too_many_arguments` are allowed for the reason
// `Cargo.toml`'s lint section gives for leaving `pedantic` off: on code that
// indexes three parallel `[f32; 3]`s by one channel index, and on a resolve pass
// whose inputs really are eight distinct buffers, the lint fights the domain.
// `experiment_p74.rs:109-110` set the precedent.
#![allow(
    clippy::cast_possible_truncation,
    clippy::cast_precision_loss,
    clippy::cast_sign_loss,
    clippy::needless_range_loop,
    clippy::too_many_arguments,
    clippy::too_many_lines
)]

mod common;

use std::hint::black_box;
use std::time::Instant;

use isomesh::Sdf;

// ---------------------------------------------------------------------------
// game_dig's constants, restated with the source named.
// ---------------------------------------------------------------------------

/// `game_dig`'s `TriplanarExtension::settings.x`: world units per texture tile
/// (`game_dig.rs:992`).
const TRIPLANAR_TILE: f32 = 1.5;
/// `settings.y`: the blend sharpness (`game_dig.rs:992`).
const TRIPLANAR_SHARPNESS: f32 = 4.0;
/// `game_dig::LAYER_BLEND` (`game_dig.rs:577`): the sentinel that makes
/// `triplanar.wgsl` blend the three terrain strata instead of forcing one layer.
const LAYER_BLEND: f32 = -1.0;
/// `game_dig::LAYER_CONCRETE` (`game_dig.rs:579`): the sandbox walls' forced
/// array layer.
const LAYER_CONCRETE: f32 = 3.0;
/// `game_dig::TERRAIN_LAYERS` (`game_dig.rs:575`).
const TERRAIN_LAYERS: usize = 4;
/// `triplanar.wgsl:76`.
const GRASS_SLOPE_LO: f32 = 0.55;
/// `triplanar.wgsl:77`.
const GRASS_SLOPE_HI: f32 = 0.82;
/// `triplanar.wgsl:83`.
const SHALLOW_Y_LO: f32 = -1.6;
/// `triplanar.wgsl:84`.
const SHALLOW_Y_HI: f32 = -0.4;
/// `game_dig`'s default `World::radius` (`game_dig.rs:1056`).
const BRUSH_RADIUS: f32 = 0.25;
/// `game_dig::LIPSCHITZ` (`game_dig.rs:190`).
const LIPSCHITZ: f32 = 1.25;
/// `game_dig::AIM_NEAR` (`game_dig.rs:216`).
const AIM_NEAR: f32 = 0.30;
/// `game_dig::AIM_FAR` (`game_dig.rs:220`).
const AIM_FAR: f32 = 25.0;
/// `game_dig::AIM_STEPS` (`game_dig.rs:223`).
const AIM_STEPS: u32 = 128;
/// `game_dig::AIM_HIT` (`game_dig.rs:226`).
const AIM_HIT: f32 = 0.01;
/// `game_dig`'s eye at `setup`: `Transform::from_xyz(0.0, 1.70, 6.0)`.
const EYE_START: [f32; 3] = [0.0, 1.70, 6.0];
/// `game_dig`'s opening `Look::pitch` (`game_dig.rs:700`).
const DEMO_PITCH: f32 = -0.15;
/// `game_dig`'s unmodified walk speed.
const WALK_SPEED: f32 = 2.5;
/// Lower corner of the sandbox `game_dig::sandbox` computes.
const SANDBOX_LO: [f32; 3] = [-8.0, -5.4, -8.0];
/// Upper corner of the same box.
const SANDBOX_HI: [f32; 3] = [8.0, 2.6, 8.0];

// ---------------------------------------------------------------------------
// This harness's own knobs.
// ---------------------------------------------------------------------------

/// C1 names 1920x1080 and this is it. Also `P-75`'s buffer, so the two rows'
/// per-pixel shading costs are denominated in the same pixels.
const SCREEN_W: usize = 1920;
/// Height of the same buffer.
const SCREEN_H: usize = 1080;
/// Width of the temporal stage's buffer. `P-77`'s, so its `0.0049` is the number
/// this harness's deterministic stream is checked against.
const TAA_W: usize = 960;
/// Height of the same.
const TAA_H: usize = 540;
/// Timed repetitions of each shading arm. Arms are interleaved inside the rep
/// loop and the **median** is reported, so a sibling agent's build lands on
/// every arm rather than on one (`M-281`). Nine, matching `P-75`'s `shade_reps`,
/// and not five: at five the baseline arm's max-minus-min was **26% of its own
/// median** with three siblings compiling, and C1's cost clause is decided at a
/// threshold only 10% below the measured ratio. `reps_ms_min` and
/// `reps_ms_spread` are both columns so the contamination is visible rather than
/// averaged away, and `c1_cost_ratio_min_based` reports the ratio of the fastest
/// observed pass per arm — the least-contended estimate available.
const SHADE_REPS: usize = 9;
/// Frames of TAA history built before anything is measured. `P-77`'s.
const WARMUP: usize = 24;
/// **The registered `frames_after_edit`.** Frames in the measurement window,
/// starting with the frame the brush lands in.
const WINDOW: usize = 8;
/// TAA blend weight for the current frame. Karis's resolve keeps ~0.9.
const ALPHA: f32 = 0.1;
/// Simulated frame time, 60 Hz.
const DT: f32 = 1.0 / 60.0;
/// Vertical field of view, radians. Bevy's `PerspectiveProjection::default`.
const FOV_Y: f32 = core::f32::consts::FRAC_PI_4;
/// Yaw rate for the walking arm, rad/s. `P-77`'s, and not a `game_dig` constant:
/// a pure forward walk leaves the reprojection stationary where the brush is.
const YAW_RATE: f32 = 0.25;
/// Central-difference half-step for the shading normal. `P-77`'s `1e-3` rather
/// than `game_dig::GRADIENT_EPS = 1e-4`, for `P-77`'s reason: a central
/// difference of an `f32` field at `1e-4` leaves three significant digits and a
/// shading signal built from that noise measures float differencing.
const NORMAL_EPS: f32 = 1e-3;
/// Absolute signal difference above which a pixel counts as changed by the dig.
const CHANGE_TOL: f32 = 1e-4;
/// Slack, in pixels, on the brush-silhouette containment control. `P-77`'s.
const SILHOUETTE_SLACK: f32 = 2.0;
/// Texture array edge, in texels. `game_dig`'s layers are `512`-square
/// (`game_dig.rs:748`). A power of two, which is what lets the wrap be a mask
/// rather than a division — the same address mode the sampler uses.
const TEX_SIZE: usize = 512;
/// Sun direction, `P-77`'s `[0.35, 0.85, 0.40]` normalised.
const SUN_DIR: [f32; 3] = [0.349_128_2, 0.847_882_8, 0.399_003_6];
/// Warm direct term, `P-77`'s.
const ALBEDO_SUN: [f32; 3] = [1.00, 0.92, 0.78];
/// Cool hemispheric term, `P-77`'s.
const ALBEDO_SKY: [f32; 3] = [0.22, 0.30, 0.45];
/// View-dependent rim term, `P-77`'s.
const ALBEDO_RIM: [f32; 3] = [0.10, 0.10, 0.12];
/// Fog colour, `P-77`'s.
const FOG_RGB: [f32; 3] = [0.52, 0.58, 0.66];
/// Fog density per unit, `P-77`'s. In because a hole that reveals a farther
/// surface must move the signal.
const FOG_DENSITY: f32 = 0.06;
/// Draws per fragment in the estimator control.
const UNBIASED_DRAWS: usize = 256;
/// Fragments the estimator control runs over.
const UNBIASED_FRAGMENTS: usize = 4096;
/// Tolerance on the estimator's **mean signed error** in the quantities it
/// selects.
///
/// The per-draw standard deviation is of order the inter-plane albedo spread,
/// ~0.1, so the mean over 256 draws and 4096 fragments has a standard error of
/// `0.1 / (16 * 64) = 1e-4`. This is more than ten sigma, and a wrong CDF —
/// uniform selection in place of weighted — misses by of order the spread
/// itself, 500 times further out. **This is a tolerance on a mean and not on a
/// max**, because unbiasedness is a statement about a mean; the max over 4096
/// fragments is an extreme-value statistic and is reported rather than gated.
const UNBIASED_BIAS_TOL: f32 = 0.002;
/// Name of the C2 arm the registered temporal columns are taken from. `P-77`'s
/// `HEADLINE_ARM`, for `P-77`'s reason: it is the only arm whose steady-state
/// rejection rate (0.49%) leaves the denominator any headroom.
const HEADLINE_TAA_ARM: &str = "dig_at_feet_static";
/// Name of the arm that proves a zero was reachable (`M-44`).
const ZERO_PROOF_ARM: &str = "zero_proof_no_jitter_static";

// ---------------------------------------------------------------------------
// Small vector helpers. `isomesh`'s public API is `[f32; 3]` and the crate's
// whole pitch is that it needs no math library.
// ---------------------------------------------------------------------------

fn add(a: [f32; 3], b: [f32; 3]) -> [f32; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}
fn sub(a: [f32; 3], b: [f32; 3]) -> [f32; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
fn mul(a: [f32; 3], s: f32) -> [f32; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}
fn dot(a: [f32; 3], b: [f32; 3]) -> f32 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
fn norm(a: [f32; 3]) -> [f32; 3] {
    let l = dot(a, a).sqrt();
    if l == 0.0 { a } else { mul(a, 1.0 / l) }
}
/// WGSL's `smoothstep`.
fn smoothstep(e0: f32, e1: f32, x: f32) -> f32 {
    let t = ((x - e0) / (e1 - e0)).clamp(0.0, 1.0);
    t * t * (3.0 - 2.0 * t)
}

// ---------------------------------------------------------------------------
// The field: game_dig's, verbatim.
// ---------------------------------------------------------------------------

/// `game_dig::Ground` (`game_dig.rs:171-175`): distance to a wavy height field,
/// negative below it.
#[derive(Clone, Copy)]
struct Ground;

impl Sdf for Ground {
    type Scalar = f32;

    fn sample(&self, p: [f32; 3]) -> f32 {
        let height = 0.35 * (p[0] * 0.9).sin() * (p[2] * 0.7).cos() + 0.15 * (p[0] * 2.1).sin();
        p[1] - height
    }
}

/// `Ground` with spheres subtracted: `max(f, -(|p - c| - r))`, the carve
/// `game_dig`'s `BrushStack` performs.
struct Dug<'a> {
    brushes: &'a [[f32; 4]],
}

impl Sdf for Dug<'_> {
    type Scalar = f32;

    fn sample(&self, p: [f32; 3]) -> f32 {
        let mut f = Ground.sample(p);
        for b in self.brushes {
            let d = [p[0] - b[0], p[1] - b[1], p[2] - b[2]];
            let sphere = (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt() - b[3];
            f = f.max(-sphere);
        }
        f
    }
}

// ---------------------------------------------------------------------------
// Camera: P-77's, which is game_dig's.
// ---------------------------------------------------------------------------

/// A pinhole camera at one instant, with the frame's jitter baked in.
#[derive(Clone, Copy)]
struct Cam {
    eye: [f32; 3],
    right: [f32; 3],
    up: [f32; 3],
    forward: [f32; 3],
    tan_half: f32,
    aspect: f32,
    width: f32,
    height: f32,
    jitter: [f32; 2],
}

impl Cam {
    /// Bevy's `Quat::from_euler(EulerRot::YXZ, yaw, pitch, 0.0)` written out:
    /// `R = Ry(yaw) * Rx(pitch)`, and the camera looks down local `-Z`.
    fn new(eye: [f32; 3], yaw: f32, pitch: f32, w: usize, h: usize, jitter: [f32; 2]) -> Self {
        let (sy, cy) = (yaw.sin(), yaw.cos());
        let (sp, cp) = (pitch.sin(), pitch.cos());
        Self {
            eye,
            right: [cy, 0.0, -sy],
            up: [sy * sp, cp, cy * sp],
            forward: [-sy * cp, sp, -cy * cp],
            tan_half: (FOV_Y * 0.5).tan(),
            aspect: w as f32 / h as f32,
            width: w as f32,
            height: h as f32,
            jitter,
        }
    }

    /// Horizontal forward, which is the direction `game_dig`'s walk integrates.
    fn walk_dir(yaw: f32) -> [f32; 3] {
        norm([-yaw.sin(), 0.0, -yaw.cos()])
    }

    /// The ray through pixel `(x, y)`, jitter included.
    fn ray(&self, x: usize, y: usize) -> [f32; 3] {
        let px = (x as f32 + 0.5 + self.jitter[0]) / self.width;
        let py = (y as f32 + 0.5 + self.jitter[1]) / self.height;
        let sx = (2.0 * px - 1.0) * self.aspect * self.tan_half;
        let sy = (1.0 - 2.0 * py) * self.tan_half;
        norm(add(
            self.forward,
            add(mul(self.right, sx), mul(self.up, sy)),
        ))
    }

    /// Where a world point lands in this frame's pixel grid, or `None` if it is
    /// behind the near plane. The jitter is subtracted, because the jittered
    /// sample at pixel `p` *is* the ray through `p + jitter`.
    fn project(&self, p: [f32; 3]) -> Option<[f32; 2]> {
        let v = sub(p, self.eye);
        let cz = dot(v, self.forward);
        if cz <= AIM_NEAR {
            return None;
        }
        let cx = dot(v, self.right);
        let cy = dot(v, self.up);
        let ndc_x = cx / (cz * self.tan_half * self.aspect);
        let ndc_y = cy / (cz * self.tan_half);
        Some([
            (ndc_x + 1.0) * 0.5 * self.width - 0.5 - self.jitter[0],
            (1.0 - ndc_y) * 0.5 * self.height - 0.5 - self.jitter[1],
        ])
    }

    /// Screen-space radius, in pixels, of a sphere's silhouette. Exact for the
    /// silhouette cone: its half-angle is `asin(r / d)`.
    fn silhouette_radius(&self, centre: [f32; 3], r: f32) -> f32 {
        let d = sub(centre, self.eye);
        let d = dot(d, d).sqrt();
        if d <= r {
            return f32::INFINITY;
        }
        (r / d).asin().tan() / self.tan_half * 0.5 * self.height
    }
}

/// R2 low-discrepancy sequence, the 10-entry jitter pattern Karis's paper uses.
fn jitter_for(frame: usize) -> [f32; 2] {
    const A1: f32 = 0.754_877_7;
    const A2: f32 = 0.569_840_3;
    let k = (frame % 10 + 1) as f32;
    [(k * A1).fract() - 0.5, (k * A2).fract() - 0.5]
}

// ---------------------------------------------------------------------------
// Ray casting: game_dig::trace, transcribed via P-77.
// ---------------------------------------------------------------------------

/// First surface crossing along a ray inside the sandbox, as a distance.
///
/// A transcription of `game_dig::trace` (`game_dig.rs:1373-1394`), including its
/// box test riding along with the surface test rather than gating the march.
fn trace(field: &Dug<'_>, origin: [f32; 3], direction: [f32; 3]) -> Option<f32> {
    let mut t = AIM_NEAR;
    for _ in 0..AIM_STEPS {
        let p = add(origin, mul(direction, t));
        let f = field.sample(p);
        let inside = (0..3).all(|i| p[i] >= SANDBOX_LO[i] && p[i] <= SANDBOX_HI[i]);
        if f <= AIM_HIT && inside {
            return Some(t);
        }
        t += (f / LIPSCHITZ).max(AIM_HIT);
        if t > AIM_FAR {
            return None;
        }
    }
    None
}

/// Central-difference gradient, normalised.
fn normal_at(field: &Dug<'_>, p: [f32; 3]) -> [f32; 3] {
    let mut g = [0.0f32; 3];
    for i in 0..3 {
        let mut a = p;
        let mut b = p;
        a[i] += NORMAL_EPS;
        b[i] -= NORMAL_EPS;
        g[i] = field.sample(a) - field.sample(b);
    }
    norm(g)
}

/// One fragment of the visibility buffer. Everything the material stage reads
/// and nothing it computes, so the interpolation both stages share is outside
/// every timer (`P-75`'s `interp_ms` split, same reason).
#[derive(Clone, Copy, Default)]
struct Frag {
    hit: bool,
    /// World-space hit point: the triplanar projection's whole input, and the
    /// motion vector's.
    world: [f32; 3],
    /// Geometric normal, normalised.
    n: [f32; 3],
    /// Distance along the ray, for the fog term.
    t: f32,
    /// Ray direction, for the rim and specular terms.
    ray: [f32; 3],
}

/// Ray-cast one frame into a visibility buffer. Rows go to `std::thread`s; the
/// field is pure and shared immutably.
fn render(cam: &Cam, brushes: &[[f32; 4]], w: usize, h: usize, threads: usize) -> Vec<Frag> {
    let mut buf = vec![Frag::default(); w * h];
    let rows_per = h.div_ceil(threads);
    std::thread::scope(|s| {
        for (chunk_index, chunk) in buf.chunks_mut(rows_per * w).enumerate() {
            let y0 = chunk_index * rows_per;
            s.spawn(move || {
                let field = Dug { brushes };
                for (local_y, row) in chunk.chunks_mut(w).enumerate() {
                    let y = y0 + local_y;
                    for (x, frag) in row.iter_mut().enumerate() {
                        let ray = cam.ray(x, y);
                        let Some(t) = trace(&field, cam.eye, ray) else {
                            continue;
                        };
                        let world = add(cam.eye, mul(ray, t));
                        *frag = Frag {
                            hit: true,
                            world,
                            n: normal_at(&field, world),
                            t,
                            ray,
                        };
                    }
                }
            });
        }
    });
    buf
}

// ---------------------------------------------------------------------------
// Noise, for the texture content. P-77's hash and value noise.
// ---------------------------------------------------------------------------

/// One 32-bit integer hash of a lattice cell, in `[0, 1)`. Splitmix64's
/// finaliser truncated to 32 bits.
fn hash_lattice(i: i32, j: i32, k: i32, seed: u32) -> u32 {
    let mut h = (i as u32)
        .wrapping_mul(0x9E37_79B1)
        ^ (j as u32).wrapping_mul(0x85EB_CA6B)
        ^ (k as u32).wrapping_mul(0xC2B2_AE35)
        ^ seed.wrapping_mul(0x27D4_EB2F);
    h ^= h >> 15;
    h = h.wrapping_mul(0x2545_F491);
    h ^= h >> 13;
    h = h.wrapping_mul(0x9E37_79B1);
    h ^= h >> 16;
    h
}

fn hash_unit(i: i32, j: i32, k: i32, seed: u32) -> f32 {
    (hash_lattice(i, j, k, seed) >> 8) as f32 / 16_777_216.0
}

/// Trilinear value noise on a unit lattice, smoothstepped.
fn value_noise(p: [f32; 3], seed: u32) -> f32 {
    let fl = [p[0].floor(), p[1].floor(), p[2].floor()];
    let (i, j, k) = (fl[0] as i32, fl[1] as i32, fl[2] as i32);
    let f = [p[0] - fl[0], p[1] - fl[1], p[2] - fl[2]];
    let s = [
        f[0] * f[0] * (3.0 - 2.0 * f[0]),
        f[1] * f[1] * (3.0 - 2.0 * f[1]),
        f[2] * f[2] * (3.0 - 2.0 * f[2]),
    ];
    let mut acc = 0.0;
    for dz in 0..2 {
        let wz = if dz == 0 { 1.0 - s[2] } else { s[2] };
        for dy in 0..2 {
            let wy = if dy == 0 { 1.0 - s[1] } else { s[1] };
            for dx in 0..2 {
                let wx = if dx == 0 { 1.0 - s[0] } else { s[0] };
                acc += wx * wy * wz * hash_unit(i + dx, j + dy, k + dz, seed);
            }
        }
    }
    acc
}

/// Fractal value noise in `[-1, 1]`, lacunarity 2, gain 0.5, at `base` cells per
/// unit.
fn fbm(p: [f32; 3], octaves: u32, base: f32, seed: u32) -> f32 {
    let mut freq = base;
    let mut amp = 1.0;
    let mut sum = 0.0;
    let mut total = 0.0;
    for o in 0..octaves {
        sum += amp * (value_noise(mul(p, freq), seed.wrapping_add(o * 7919)) * 2.0 - 1.0);
        total += amp;
        freq *= 2.0;
        amp *= 0.5;
    }
    sum / total
}

// ---------------------------------------------------------------------------
// The texture arrays.
// ---------------------------------------------------------------------------

/// A `TEX_SIZE`-square, [`TERRAIN_LAYERS`]-layer RGBA8 array: the same shape,
/// element type and 4 MiB working set as `game_dig`'s packed PNG
/// (`game_dig.rs:748-755`, reinterpreted at `game_dig.rs:786-797`).
///
/// Layer-major, which is how a `texture_2d_array` is laid out and why a
/// per-fragment stratum blend touches three pages 1 MiB apart.
struct TexArray {
    texels: Vec<[u8; 4]>,
}

impl TexArray {
    /// Fill with fractal value noise, one distinct stratum per layer.
    ///
    /// The frequency content matters and is not arbitrary: `game_dig` tiles at
    /// `1.5` world units, so a `512`-square layer carries content down to
    /// `1.5/512 = 0.0029` units, and `P-77`'s first run proved that a *smooth*
    /// signal turns the TAA neighbourhood into a needle and makes the rejection
    /// rate a measurement of the shading model (95.5%, `experiment_p77.rs:159-164`).
    /// Six octaves at 8 cells per tile gives a finest period of `1/256` of the
    /// tile, which is the demo's own detail rate to within a factor of two.
    fn new(seed: u32) -> Self {
        let mut texels = vec![[0u8; 4]; TEX_SIZE * TEX_SIZE * TERRAIN_LAYERS];
        for layer in 0..TERRAIN_LAYERS {
            // A different base tint per layer, so the strata are as far apart in
            // colour as grass, dirt, deep dirt and concrete are. Without that the
            // stratum blend is invisible and the stochastic *layer* arm would
            // report a bias of zero for the wrong reason.
            let tint = [
                [0.30f32, 0.45, 0.20],
                [0.42, 0.32, 0.22],
                [0.28, 0.22, 0.18],
                [0.62, 0.62, 0.60],
            ][layer];
            for y in 0..TEX_SIZE {
                for x in 0..TEX_SIZE {
                    let p = [
                        x as f32 / TEX_SIZE as f32,
                        y as f32 / TEX_SIZE as f32,
                        layer as f32 * 4.0,
                    ];
                    let luma = fbm(p, 6, 8.0, seed);
                    let chroma = fbm(p, 3, 4.0, seed ^ 0x5EED_BEEF);
                    let rough = 0.5 + 0.35 * fbm(p, 2, 2.0, seed ^ 0x1234_5678);
                    let mut t = [0u8; 4];
                    for c in 0..3 {
                        let v = tint[c] * (1.0 + 0.55 * luma) + 0.10 * chroma * (c as f32 - 1.0);
                        t[c] = (v.clamp(0.0, 1.0) * 255.0) as u8;
                    }
                    t[3] = (rough.clamp(0.0, 1.0) * 255.0) as u8;
                    texels[layer * TEX_SIZE * TEX_SIZE + y * TEX_SIZE + x] = t;
                }
            }
        }
        Self { texels }
    }

    /// One texel, with the sampler's repeat address mode. `TEX_SIZE` is a power
    /// of two, so the wrap is a mask; `i32 & (TEX_SIZE - 1)` is correct for
    /// negative indices in two's complement, which is exactly what the hardware
    /// does.
    fn texel(&self, layer: usize, x: i32, y: i32) -> [f32; 4] {
        let xx = (x & (TEX_SIZE as i32 - 1)) as usize;
        let yy = (y & (TEX_SIZE as i32 - 1)) as usize;
        let t = self.texels[layer * TEX_SIZE * TEX_SIZE + yy * TEX_SIZE + xx];
        [
            f32::from(t[0]) * (1.0 / 255.0),
            f32::from(t[1]) * (1.0 / 255.0),
            f32::from(t[2]) * (1.0 / 255.0),
            f32::from(t[3]) * (1.0 / 255.0),
        ]
    }
}

/// How a fetch is filtered.
#[derive(Clone, Copy, PartialEq, Eq)]
enum Filter {
    /// One `textureSample`: **four** texel loads and a bilinear reconstruction.
    Bilinear,
    /// `arXiv:2305.05810`: **one** texel, chosen with probability equal to its
    /// bilinear weight. Unbiased by construction, `1/4` of the bandwidth, and
    /// noisy — which is the debt TAA is supposed to pay off.
    StochasticPoint,
    /// The ALU control: the call, the address arithmetic and the tally, no texel.
    /// Its purpose is to isolate the texel cost from the fetch-call cost, which
    /// is the difference between `S` and the fit's slope.
    NoTexel,
}

impl Filter {
    /// Texel loads per fetch instruction.
    fn texels_per_fetch(self) -> u64 {
        match self {
            Self::Bilinear => 4,
            Self::StochasticPoint => 1,
            Self::NoTexel => 0,
        }
    }
}

/// A 32-bit xorshift, seeded per fragment. Deterministic and reproducible: the
/// stochastic arms are a fixed function of `(pixel, frame)` and rerun identically.
struct Rng(u32);

impl Rng {
    fn new(x: usize, y: usize, frame: usize) -> Self {
        Self(hash_lattice(x as i32, y as i32, frame as i32, 0x9E37_79B9) | 1)
    }
    fn unit(&mut self) -> f32 {
        let mut x = self.0;
        x ^= x << 13;
        x ^= x >> 17;
        x ^= x << 5;
        self.0 = x;
        (x >> 8) as f32 / 16_777_216.0
    }
}

/// One `textureSample`, tallied. `taps` is `arXiv:2305.05810`'s tap count: `1`
/// is a single fetch instruction, `n > 1` averages `n` independent draws.
fn sample(
    tex: &TexArray,
    layer: usize,
    uv: [f32; 2],
    filter: Filter,
    taps: u32,
    rng: &mut Rng,
    tally: &mut u64,
) -> [f32; 4] {
    let x = uv[0] * TEX_SIZE as f32 - 0.5;
    let y = uv[1] * TEX_SIZE as f32 - 0.5;
    let x0 = x.floor();
    let y0 = y.floor();
    let (ix, iy) = (x0 as i32, y0 as i32);
    let (fx, fy) = (x - x0, y - y0);
    let mut out = [0.0f32; 4];
    for _ in 0..taps {
        *tally += 1;
        match filter {
            Filter::Bilinear => {
                for (dy, wy) in [(0i32, 1.0 - fy), (1, fy)] {
                    for (dx, wx) in [(0i32, 1.0 - fx), (1, fx)] {
                        let t = tex.texel(layer, ix + dx, iy + dy);
                        let k = wx * wy;
                        for c in 0..4 {
                            out[c] += k * t[c];
                        }
                    }
                }
            }
            Filter::StochasticPoint => {
                // Pick one of the four bilinear taps with probability equal to
                // its weight: `E[tap] = bilinear`, which is the paper's whole
                // claim and the reason this is a filter rather than a shortcut.
                let dx = i32::from(rng.unit() < fx);
                let dy = i32::from(rng.unit() < fy);
                let t = tex.texel(layer, ix + dx, iy + dy);
                for c in 0..4 {
                    out[c] += t[c];
                }
            }
            Filter::NoTexel => {
                for c in 0..4 {
                    out[c] += 0.5;
                }
            }
        }
    }
    let inv = 1.0 / taps as f32;
    for c in 0..4 {
        out[c] *= inv;
    }
    out
}

// ---------------------------------------------------------------------------
// The arms.
// ---------------------------------------------------------------------------

/// Which triplanar planes a fragment samples.
#[derive(Clone, Copy, PartialEq, Eq)]
enum PlaneMode {
    /// `triplanar.wgsl` as written: all three, weighted by `pow(abs(n), s)`.
    Triplanar,
    /// Two planes, Quilez's *"Biplanar mapping"* weighting: the major and median
    /// axes, weights `clamp((|n| - 1/sqrt3)/(1 - 1/sqrt3))` raised to `k/8`.
    ///
    /// **The UVs are `triplanar.wgsl`'s, not the article's**, and that is a
    /// deviation with a reason: the article writes `vec2(p[ma.y], p[ma.z])`,
    /// which for the x-major case is `p.yz` where the shader uses `p.zy`. On a
    /// tiling texture the difference is a transpose of the sample, which would
    /// put a difference between biplanar and the 3-plane reference that has
    /// nothing to do with the plane count — and C3's error term is exactly that
    /// difference. So the planes are dropped and reweighted; the projections are
    /// left alone.
    Biplanar,
    /// One plane, chosen with probability equal to its triplanar weight. The
    /// estimator is unbiased: `E[T] = sum_j w_j T_j` is the 3-plane blend.
    StochasticOne,
}

impl PlaneMode {
    fn count(self) -> u64 {
        match self {
            Self::Triplanar => 3,
            Self::Biplanar => 2,
            Self::StochasticOne => 1,
        }
    }
}

/// One shading configuration.
#[derive(Clone, Copy)]
struct Arm {
    name: &'static str,
    planes: PlaneMode,
    /// `game_dig`'s `settings.z`: negative blends the three strata, `>= 0.0`
    /// forces that one array layer.
    forced_layer: f32,
    /// Stochastic stratum selection on top of the plane selection: one stratum,
    /// chosen with probability equal to its blend weight. `lw` sums to `1`
    /// identically (`triplanar.wgsl:119-120`), so this estimator is unbiased for
    /// the same reason the plane one is.
    stochastic_layer: bool,
    /// `arXiv:2305.05810`'s tap count.
    taps: u32,
    filter: Filter,
    /// Registered as a shading arm rather than as a control.
    scored: bool,
}

impl Arm {
    fn layer_count(self) -> u64 {
        if self.forced_layer >= 0.0 || self.stochastic_layer {
            1
        } else {
            3
        }
    }

    /// Fetch instructions per fragment. `2` maps, because `triplanar.wgsl`
    /// samples `albedo_roughness_texture` and `normal_ao_texture` over the same
    /// planes and the same strata.
    fn fetches(self) -> u64 {
        if self.filter == Filter::NoTexel {
            return 0;
        }
        2 * self.planes.count() * self.layer_count() * u64::from(self.taps)
    }

    /// Fetch *calls* per fragment, which the `NoTexel` control also pays.
    fn calls(self) -> u64 {
        2 * self.planes.count() * self.layer_count() * u64::from(self.taps)
    }

    fn texels(self) -> u64 {
        self.calls() * self.filter.texels_per_fetch()
    }

    /// Whether the fragment needs a random stream at all. A deterministic arm
    /// must not be charged for a hash it would never issue.
    fn random(self) -> bool {
        self.planes == PlaneMode::StochasticOne
            || self.stochastic_layer
            || self.filter == Filter::StochasticPoint
    }
}

/// `triplanar.wgsl:64-67`.
fn plane_weights(n: [f32; 3], sharpness: f32) -> [f32; 3] {
    let w = [
        n[0].abs().powf(sharpness),
        n[1].abs().powf(sharpness),
        n[2].abs().powf(sharpness),
    ];
    let s = (w[0] + w[1] + w[2]).max(1e-5);
    [w[0] / s, w[1] / s, w[2] / s]
}

/// Quilez's biplanar weights over the major and median axes, normalised.
///
/// `clamp((|n| - 1/sqrt3) / (1 - 1/sqrt3), 0, 1)` then `pow(., k/8)`, with the
/// article's `k` taken from `game_dig`'s own blend sharpness. The sum is clamped
/// off zero for `triplanar.wgsl`'s reason: at `|n| = (1,1,1)/sqrt3` both weights
/// are zero and the article divides by their sum.
fn biplanar_weights(n: [f32; 3], sharpness: f32) -> ([usize; 2], [f32; 2]) {
    let a = [n[0].abs(), n[1].abs(), n[2].abs()];
    let major = if a[0] > a[1] && a[0] > a[2] {
        0
    } else if a[1] > a[2] {
        1
    } else {
        2
    };
    let minor = if a[0] < a[1] && a[0] < a[2] {
        0
    } else if a[1] < a[2] {
        1
    } else {
        2
    };
    // `3 - minor - major` is the article's `me = ivec3(3) - mi - ma`.
    let median = 3 - minor - major;
    const INV_SQRT3: f32 = 0.577_350_3;
    let mut w = [a[major], a[median]];
    for v in &mut w {
        *v = ((*v - INV_SQRT3) / (1.0 - INV_SQRT3))
            .clamp(0.0, 1.0)
            .powf(sharpness / 8.0);
    }
    let s = (w[0] + w[1]).max(1e-5);
    ([major, median], [w[0] / s, w[1] / s])
}

/// The two coordinates plane `j` spans, `triplanar.wgsl:113-115`.
fn plane_uv(p: [f32; 3], j: usize) -> [f32; 2] {
    match j {
        0 => [p[2], p[1]],
        1 => [p[0], p[2]],
        _ => [p[0], p[1]],
    }
}

/// One plane's tangent normal folded into the geometric normal and swizzled back
/// to world space: `triplanar.wgsl:167-172`, written out per plane.
fn whiteout(j: usize, na: [f32; 4], n: [f32; 3]) -> [f32; 3] {
    let t = [na[0] * 2.0 - 1.0, na[1] * 2.0 - 1.0, na[2] * 2.0 - 1.0];
    match j {
        // `t_x = vec3(t_x.xy + n.zy, abs(t_x.z) * n.x)`, contributed as `.zyx`.
        0 => [t[2].abs() * n[0], t[1] + n[1], t[0] + n[2]],
        // `t_y = vec3(t_y.xy + n.xz, abs(t_y.z) * n.y)`, contributed as `.xzy`.
        1 => [t[0] + n[0], t[2].abs() * n[1], t[1] + n[2]],
        // `t_z = vec3(t_z.xy + n.xy, abs(t_z.z) * n.z)`, contributed as `.xyz`.
        _ => [t[0] + n[0], t[1] + n[1], t[2].abs() * n[2]],
    }
}

/// One shaded fragment, plus the two quantities the sampler is *contractually*
/// unbiased in.
///
/// **The split matters and it is the point of `arXiv:2305.05810`'s title.**
/// Stochastic selection is unbiased in what it selects: `E[ar] = sum_l lw_l
/// sum_j w_j T(j, l)`, exactly the 3-plane blend, and the same for the
/// pre-normalisation whiteout sum. It is **not** unbiased in the shaded colour,
/// because `normalize`, `max(0, .)`, `powf` and the roughness clamp are all
/// nonlinear — *filtering after shading is not filtering before shading*. So the
/// estimator control asserts the contract on [`Shaded::ar`] and
/// [`Shaded::n_sum`], and the residual bias in [`Shaded::rgb`] is **reported as a
/// finding** rather than asserted away: it is an error floor no amount of
/// temporal accumulation removes, which is the opposite of the variance TAA is
/// supposed to absorb.
#[derive(Clone, Copy)]
struct Shaded {
    /// The shaded colour. A nonlinear function of the samples.
    rgb: [f32; 3],
    /// `ar` before the roughness clamp: linear in the sampled texels.
    ar: [f32; 4],
    /// The whiteout blend's sum before `normalize`: linear in the sampled texels.
    n_sum: [f32; 3],
}

/// The material shading stage for one fragment: `triplanar.wgsl`'s sampling and
/// blending, then a lighting term.
///
/// **The lighting is `P-77`'s geometry-only model plus a Blinn-Phong lobe, not
/// `apply_pbr_lighting`.** It is identical in every arm, so it is entirely inside
/// `F` and cannot bias the comparison; it is *smaller* than Bevy's, so the ratio
/// measured here is an upper bound on a real pipeline's.
fn shade(
    arm: &Arm,
    ar_tex: &TexArray,
    na_tex: &TexArray,
    f: &Frag,
    rng: &mut Rng,
    tally: &mut u64,
) -> Shaded {
    let p = mul(f.world, 1.0 / TRIPLANAR_TILE);
    let n = f.n;

    // Plane selection and weights.
    //
    // **`plane_weights` is called inside the arms that use it and not before the
    // match, and the first nine-rep run is why.** Hoisted, the biplanar arm paid
    // `triplanar.wgsl`'s three `pow(abs(n), 4)` on top of its own two — a cost it
    // would never issue on hardware, since Quilez's weighting replaces that blend
    // rather than following it. C3's wall-clock saving ratio read **0.4715**
    // against its `0.5` threshold with that charge in, on both the median and the
    // fastest passes, so the defect was large enough to decide the clause on its
    // own.
    let mut planes = [0usize; 3];
    let mut pw = [0.0f32; 3];
    let np = match arm.planes {
        PlaneMode::Triplanar => {
            planes = [0, 1, 2];
            pw = plane_weights(n, TRIPLANAR_SHARPNESS);
            3
        }
        PlaneMode::Biplanar => {
            let (idx, ws) = biplanar_weights(n, TRIPLANAR_SHARPNESS);
            planes[0] = idx[0];
            planes[1] = idx[1];
            pw[0] = ws[0];
            pw[1] = ws[1];
            2
        }
        PlaneMode::StochasticOne => {
            let w3 = plane_weights(n, TRIPLANAR_SHARPNESS);
            let xi = rng.unit();
            let mut acc = 0.0;
            let mut pick = 2;
            for j in 0..3 {
                acc += w3[j];
                if xi < acc {
                    pick = j;
                    break;
                }
            }
            planes[0] = pick;
            pw[0] = 1.0;
            1
        }
    };

    // Stratum selection and weights: `triplanar.wgsl:121-124`.
    let up = smoothstep(GRASS_SLOPE_LO, GRASS_SLOPE_HI, n[1]);
    let shallow = smoothstep(SHALLOW_Y_LO, SHALLOW_Y_HI, f.world[1]);
    let lw3 = [up * shallow, (1.0 - up) * shallow, 1.0 - shallow];
    let mut layers = [0usize; 3];
    let mut lw = [0.0f32; 3];
    let nl = if arm.forced_layer >= 0.0 {
        layers[0] = arm.forced_layer as usize;
        lw[0] = 1.0;
        1
    } else if arm.stochastic_layer {
        let xi = rng.unit();
        let mut acc = 0.0;
        let mut pick = 2;
        for l in 0..3 {
            acc += lw3[l];
            if xi < acc {
                pick = l;
                break;
            }
        }
        layers[0] = pick;
        lw[0] = 1.0;
        1
    } else {
        layers = [0, 1, 2];
        lw = lw3;
        3
    };

    // Albedo/roughness: `layer_ar` per stratum, all selected planes inside.
    let mut ar = [0.0f32; 4];
    for li in 0..nl {
        for ji in 0..np {
            let s = sample(
                ar_tex,
                layers[li],
                plane_uv(p, planes[ji]),
                arm.filter,
                arm.taps,
                rng,
                tally,
            );
            let k = lw[li] * pw[ji];
            for c in 0..4 {
                ar[c] += k * s[c];
            }
        }
    }

    // Normal/AO: `layers_na` per plane, all selected strata inside, then the
    // whiteout blend. Per plane rather than per fragment, because a normal
    // blended across planes has already lost which plane it came from.
    let mut nn = [0.0f32; 3];
    let mut ao = 0.0f32;
    for ji in 0..np {
        let uv = plane_uv(p, planes[ji]);
        let mut na = [0.0f32; 4];
        for li in 0..nl {
            let s = sample(
                na_tex,
                layers[li],
                uv,
                arm.filter,
                arm.taps,
                rng,
                tally,
            );
            for c in 0..4 {
                na[c] += lw[li] * s[c];
            }
        }
        nn = add(nn, mul(whiteout(planes[ji], na, n), pw[ji]));
        ao += pw[ji] * na[3];
    }
    let big_n = norm(nn);

    // Lighting. `perceptual_roughness = clamp(ar.a, 0.05, 1.0)`,
    // `triplanar.wgsl:156`.
    let rough = ar[3].clamp(0.05, 1.0);
    let sun = dot(big_n, SUN_DIR).max(0.0);
    let sky = 0.5 * (big_n[1] + 1.0) * ao;
    let rim = 1.0 - dot(big_n, mul(f.ray, -1.0)).max(0.0);
    let half = norm(sub(SUN_DIR, f.ray));
    let spec = dot(big_n, half).max(0.0).powf(2.0 / (rough * rough)) * (1.0 - rough);
    let fog = 1.0 - (-f.t * FOG_DENSITY).exp();
    let mut out = [0.0f32; 3];
    for c in 0..3 {
        let lit = (ALBEDO_SUN[c] * sun + ALBEDO_SKY[c] * sky + ALBEDO_RIM[c] * rim) * ar[c] + spec;
        out[c] = lit * (1.0 - fog) + FOG_RGB[c] * fog;
    }
    Shaded {
        rgb: out,
        ar,
        n_sum: nn,
    }
}

// ---------------------------------------------------------------------------
// TAA: P-77's resolve, with its baseline volume.
// ---------------------------------------------------------------------------

/// Karis 2014's YCoCg, the space production TAA builds its neighbourhood in.
fn to_ycocg(c: [f32; 3]) -> [f32; 3] {
    [
        0.25 * c[0] + 0.5 * c[1] + 0.25 * c[2],
        0.5 * (c[0] - c[2]),
        -0.25 * c[0] + 0.5 * c[1] - 0.25 * c[2],
    ]
}

/// Inverse of [`to_ycocg`].
fn from_ycocg(c: [f32; 3]) -> [f32; 3] {
    [c[0] + c[1] - c[2], c[0] + c[2], c[0] - c[1] - c[2]]
}

/// The 3x3 neighbourhood's axis-aligned bounding box, and the clip.
///
/// `P-77`'s `aabb_ycocg` — the baseline, because that is what production TAA
/// does. Its k-DOP arms are not repeated: `P-77` showed a k-DOP that contains the
/// AABB's own axes is a subset of it and recovers exactly zero.
struct Aabb {
    lo: [f32; 3],
    hi: [f32; 3],
}

/// The floating-point dilation of every extent, `P-77`'s [`Dop`]-equivalent
/// `DOP_EPS`, which is `10.1145/3681758.3697996`'s own `eps = 1e-5`.
///
/// **This harness's own zero-proof control is what proved it is not decoration.**
/// Dropped, the run rejected **31,274** history samples in the arm whose
/// reprojection is the identity — because `from_ycocg(to_ycocg(c))` is the
/// identity in real arithmetic and not in `f32`, so the history drifts about
/// `1e-7` per frame, and for any pixel that is the extremum of its own 3x3 in
/// some channel the box boundary passes exactly through the current colour and
/// clips that `1e-7` to `s ~ 0`. `1e-5` is two orders above the drift and five
/// below the signal, and with it restored the same arm rejects zero. `P-77` had
/// it and read `0.00000000`; this is the fixture defect this harness's own
/// control caught.
const CLAMP_EPS: f32 = 1e-5;

impl Aabb {
    fn build(colours: &[[f32; 3]]) -> Self {
        let mut lo = [f32::INFINITY; 3];
        let mut hi = [f32::NEG_INFINITY; 3];
        for c in colours {
            for i in 0..3 {
                lo[i] = lo[i].min(c[i]);
                hi[i] = hi[i].max(c[i]);
            }
        }
        for i in 0..3 {
            lo[i] -= CLAMP_EPS;
            hi[i] += CLAMP_EPS;
        }
        Self { lo, hi }
    }

    /// Clip the ray `centre -> history` to the shell. `1.0` means the history was
    /// inside and nothing was rejected. The ray starts inside — the centre
    /// pixel's own colour is one of the neighbourhood colours — and a box is
    /// convex, so this is the nearest slab crossing.
    fn clip(&self, centre: [f32; 3], history: [f32; 3]) -> f32 {
        let d = sub(history, centre);
        let mut s = 1.0f32;
        for i in 0..3 {
            if d[i].abs() < 1e-20 {
                continue;
            }
            let bound = if d[i] > 0.0 { self.hi[i] } else { self.lo[i] };
            let si = (bound - centre[i]) / d[i];
            if si < s {
                s = si;
            }
        }
        s.clamp(0.0, 1.0)
    }
}

/// One TAA history buffer.
struct Taa {
    hist: Vec<[f32; 3]>,
    valid: Vec<bool>,
}

impl Taa {
    fn new(n: usize) -> Self {
        Self {
            hist: vec![[0.0; 3]; n],
            valid: vec![false; n],
        }
    }
}

/// Bilinear fetch from a history buffer, requiring all four taps valid.
fn fetch_hist(
    hist: &[[f32; 3]],
    valid: &[bool],
    w: usize,
    h: usize,
    p: [f32; 2],
) -> Option<[f32; 3]> {
    let x0 = p[0].floor();
    let y0 = p[1].floor();
    if x0 < 0.0 || y0 < 0.0 || x0 + 1.0 >= w as f32 || y0 + 1.0 >= h as f32 {
        return None;
    }
    let (ix, iy) = (x0 as usize, y0 as usize);
    let (fx, fy) = (p[0] - x0, p[1] - y0);
    let mut out = [0.0f32; 3];
    for (dy, wy) in [(0usize, 1.0 - fy), (1, fy)] {
        for (dx, wx) in [(0usize, 1.0 - fx), (1, fx)] {
            let i = (iy + dy) * w + ix + dx;
            if !valid[i] {
                return None;
            }
            let k = wx * wy;
            for c in 0..3 {
                out[c] += hist[i][c] * k;
            }
        }
    }
    Some(out)
}

/// The 3x3 neighbourhood of the current frame's colours, clamped at the border.
fn neighbourhood(
    rgb: &[[f32; 3]],
    gb: &[Frag],
    w: usize,
    h: usize,
    x: usize,
    y: usize,
) -> ([[f32; 3]; 9], usize) {
    let mut out = [[0.0f32; 3]; 9];
    let mut n = 0;
    for dy in -1i32..=1 {
        for dx in -1i32..=1 {
            let sx = (x as i32 + dx).clamp(0, w as i32 - 1) as usize;
            let sy = (y as i32 + dy).clamp(0, h as i32 - 1) as usize;
            let i = sy * w + sx;
            if gb[i].hit {
                out[n] = to_ycocg(rgb[i]);
                n += 1;
            }
        }
    }
    (out, n)
}

/// One frame's temporal counters, for one stream of one world.
#[derive(Clone, Copy, Default)]
struct Stream {
    /// Pixels with a current hit and a valid on-screen reprojection with history.
    population: u64,
    /// Pixels with a current hit and no usable history.
    no_history: u64,
    /// **The registered rejection count**: the clamp moved the reprojected
    /// sample, i.e. it lay outside the neighbourhood volume.
    rejected: u64,
    /// Sum of the clip parameter over the rejected samples. Near 1 means the
    /// history was barely outside (clamp conservatism); near 0 means genuine
    /// disocclusion.
    clip_s: f64,
    /// Sum of reprojection displacement, in pixels, over the population.
    reproj_px: f64,
    /// Hit pixels, the MAE denominator.
    hits: u64,
    /// Sum over hit pixels and channels of `|resolved - reference| / 3`.
    abs_err: f64,
}

impl Stream {
    fn mae(&self) -> f64 {
        if self.hits == 0 {
            0.0
        } else {
            self.abs_err / self.hits as f64
        }
    }
    fn fraction(&self) -> f64 {
        if self.population == 0 {
            0.0
        } else {
            self.rejected as f64 / self.population as f64
        }
    }
    fn accumulate(&mut self, o: &Self) {
        self.population += o.population;
        self.no_history += o.no_history;
        self.rejected += o.rejected;
        self.clip_s += o.clip_s;
        self.reproj_px += o.reproj_px;
        self.hits += o.hits;
        self.abs_err += o.abs_err;
    }
}

/// What one resolve pass hands back: the frame's counters, the fetched history
/// per pixel, and the rectified colour per pixel.
type Resolved = (Stream, Vec<Option<[f32; 3]>>, Vec<[f32; 3]>);

/// Reproject, clip, count. Returns the counters, the fetched history and the
/// clipped colours; the commit is separate so the history buffer is immutable
/// for the whole pass and the rows can be handed to threads.
fn resolve(
    taa: &Taa,
    cur: &[[f32; 3]],
    gb: &[Frag],
    prev_cam: Option<&Cam>,
    w: usize,
    h: usize,
    threads: usize,
) -> Resolved {
    let n = w * h;
    let mut hist_rgb: Vec<Option<[f32; 3]>> = vec![None; n];
    let mut clipped: Vec<[f32; 3]> = vec![[0.0; 3]; n];
    let mut total = Stream::default();
    let rows_per = h.div_ceil(threads);

    std::thread::scope(|s| {
        let mut handles = Vec::new();
        for (ci, (hchunk, cchunk)) in hist_rgb
            .chunks_mut(rows_per * w)
            .zip(clipped.chunks_mut(rows_per * w))
            .enumerate()
        {
            let y0 = ci * rows_per;
            handles.push(s.spawn(move || {
                let mut acc = Stream::default();
                for (ly, (hrow, crow)) in hchunk.chunks_mut(w).zip(cchunk.chunks_mut(w)).enumerate()
                {
                    let y = y0 + ly;
                    for x in 0..w {
                        let i = y * w + x;
                        if !gb[i].hit {
                            continue;
                        }
                        let cy = to_ycocg(cur[i]);
                        let history = prev_cam.and_then(|pc| {
                            pc.project(gb[i].world).and_then(|q| {
                                fetch_hist(&taa.hist, &taa.valid, w, h, q).map(|hrgb| (q, hrgb))
                            })
                        });
                        let Some((q, hrgb)) = history else {
                            acc.no_history += 1;
                            crow[x] = cur[i];
                            continue;
                        };
                        acc.population += 1;
                        let (dx, dy) = (q[0] - x as f32, q[1] - y as f32);
                        acc.reproj_px += f64::from((dx * dx + dy * dy).sqrt());
                        let (nb, k) = neighbourhood(cur, gb, w, h, x, y);
                        let vol = Aabb::build(&nb[..k]);
                        let hy = to_ycocg(hrgb);
                        let sc = vol.clip(cy, hy);
                        crow[x] = from_ycocg(add(cy, mul(sub(hy, cy), sc)));
                        hrow[x] = Some(hrgb);
                        if sc < 1.0 {
                            acc.rejected += 1;
                            acc.clip_s += f64::from(sc);
                        }
                    }
                }
                acc
            }));
        }
        for hd in handles {
            let part = hd.join().expect("a resolve worker panicked");
            total.accumulate(&part);
        }
    });

    (total, hist_rgb, clipped)
}

/// Commit the history: blend the clipped value, seed pixels with no usable
/// history from the current frame, and accumulate the image error against the
/// reference. The committed value **is** the resolved output, which is why the
/// MAE is taken here.
fn commit(
    taa: &mut Taa,
    cur: &[[f32; 3]],
    reference: &[[f32; 3]],
    gb: &[Frag],
    hist_rgb: &[Option<[f32; 3]>],
    clipped: &[[f32; 3]],
    counters: &mut Stream,
    measure: bool,
) {
    for i in 0..gb.len() {
        if !gb[i].hit {
            taa.valid[i] = false;
            continue;
        }
        taa.valid[i] = true;
        taa.hist[i] = match hist_rgb[i] {
            Some(_) => add(mul(clipped[i], 1.0 - ALPHA), mul(cur[i], ALPHA)),
            None => cur[i],
        };
        if measure {
            counters.hits += 1;
            let mut e = 0.0f64;
            for c in 0..3 {
                e += f64::from((taa.hist[i][c] - reference[i][c]).abs());
            }
            counters.abs_err += e / 3.0;
        }
    }
}

/// Shade a whole buffer, multithreaded. **Not the timed pass**: stage A's
/// single-threaded pass is the cost measurement, and stage B only needs the
/// pixels.
fn shade_buffer(
    arm: &Arm,
    ar_tex: &TexArray,
    na_tex: &TexArray,
    gb: &[Frag],
    w: usize,
    h: usize,
    frame: usize,
    threads: usize,
) -> Vec<[f32; 3]> {
    let mut out = vec![[0.0f32; 3]; w * h];
    let rows_per = h.div_ceil(threads);
    std::thread::scope(|s| {
        for (ci, chunk) in out.chunks_mut(rows_per * w).enumerate() {
            let y0 = ci * rows_per;
            s.spawn(move || {
                let mut tally = 0u64;
                for (ly, row) in chunk.chunks_mut(w).enumerate() {
                    let y = y0 + ly;
                    for (x, px) in row.iter_mut().enumerate() {
                        let f = &gb[y * w + x];
                        if !f.hit {
                            continue;
                        }
                        let mut rng = Rng::new(x, y, frame);
                        *px = shade(arm, ar_tex, na_tex, f, &mut rng, &mut tally).rgb;
                    }
                }
                black_box(tally);
            });
        }
    });
    out
}

// ---------------------------------------------------------------------------
// Stage A: the material shading stage's cost. C1 and C3.
// ---------------------------------------------------------------------------

/// What one arm's timed pass produced.
struct ShadeResult {
    /// Median of [`SHADE_REPS`] interleaved passes, ms.
    fragment_ms: f64,
    /// Every rep, in order, so a reader can see the spread on a busy machine.
    reps_ms: Vec<f64>,
    /// Fetch instructions issued, tallied inside the pass.
    fetches: u64,
    /// Mean shaded colour, so the compiler cannot elide the work and a reader
    /// can see the arms agree to within their sampling noise.
    mean_rgb: [f64; 3],
}

/// One single-threaded pass of the material shading stage over the visible
/// fragments. This is the measured quantity.
fn shade_pass(
    arm: &Arm,
    ar_tex: &TexArray,
    na_tex: &TexArray,
    frags: &[Frag],
) -> (f64, u64, [f64; 3]) {
    let mut tally = 0u64;
    let mut acc = [0.0f64; 3];
    let t0 = Instant::now();
    for (i, f) in frags.iter().enumerate() {
        let mut rng = if arm.random() {
            Rng::new(i, 0, 0)
        } else {
            Rng(1)
        };
        let out = shade(arm, ar_tex, na_tex, f, &mut rng, &mut tally);
        for c in 0..3 {
            acc[c] += f64::from(out.rgb[c]);
        }
    }
    let ms = t0.elapsed().as_secs_f64() * 1e3;
    black_box(&acc);
    let inv = 1.0 / frags.len() as f64;
    (ms, tally, [acc[0] * inv, acc[1] * inv, acc[2] * inv])
}

/// Least-squares fit of `t = F + k*S` over `(fetches, ms)` pairs, with `R^2`.
fn fit(points: &[(f64, f64)]) -> (f64, f64, f64) {
    let n = points.len() as f64;
    let sx: f64 = points.iter().map(|p| p.0).sum();
    let sy: f64 = points.iter().map(|p| p.1).sum();
    let sxx: f64 = points.iter().map(|p| p.0 * p.0).sum();
    let sxy: f64 = points.iter().map(|p| p.0 * p.1).sum();
    let denom = n * sxx - sx * sx;
    let slope = (n * sxy - sx * sy) / denom;
    let intercept = (sy - slope * sx) / n;
    let mean = sy / n;
    let ss_tot: f64 = points.iter().map(|p| (p.1 - mean).powi(2)).sum();
    let ss_res: f64 = points
        .iter()
        .map(|p| (p.1 - (intercept + slope * p.0)).powi(2))
        .sum();
    let r2 = if ss_tot == 0.0 {
        1.0
    } else {
        1.0 - ss_res / ss_tot
    };
    (intercept, slope, r2)
}

// ---------------------------------------------------------------------------
// Stage B: the temporal instrument. C2.
// ---------------------------------------------------------------------------

/// One temporal arm.
struct TaaArm {
    name: &'static str,
    /// Camera pitch, radians. `-0.6` is `P-77`'s headline pose — `game_dig`'s
    /// `Look::pitch` is user-driven and clamped to +-1.5, so this is the demo
    /// with the player looking down at the rock under their feet.
    pitch: f32,
    walk: f32,
    yaw_rate: f32,
    jitter: bool,
    /// Whether a brush lands at the start of the window. `false` is the
    /// zero-proof arm.
    edit: bool,
}

/// One arm's output: per-frame counters for both worlds and all three streams.
struct TaaOut {
    /// `[world][stream]`, world 0 = static, world 1 = digging; stream 0 =
    /// reference, 1 = stochastic, 2 = biplanar.
    frames: [[Vec<Stream>; 3]; 2],
    /// Pixels the dig changed, per window frame.
    changed: Vec<u64>,
    /// Changed pixels outside every brush silhouette. Must be zero.
    changed_outside: u64,
    brushes_placed: usize,
}

/// Simulate the static and dug worlds in lockstep over one camera path and one
/// jitter sequence, resolving three shading streams per world.
fn run_taa_arm(
    cfg: &TaaArm,
    arms: [&Arm; 3],
    ar_tex: &TexArray,
    na_tex: &TexArray,
    threads: usize,
) -> TaaOut {
    let (w, h) = (TAA_W, TAA_H);
    let n = w * h;
    let mut taa: Vec<Taa> = (0..6).map(|_| Taa::new(n)).collect();
    let mut prev_cam: Option<Cam> = None;
    let mut brushes: Vec<[f32; 4]> = Vec::new();
    let mut out = TaaOut {
        frames: [
            [Vec::new(), Vec::new(), Vec::new()],
            [Vec::new(), Vec::new(), Vec::new()],
        ],
        changed: Vec::new(),
        changed_outside: 0,
        brushes_placed: 0,
    };

    let mut eye = EYE_START;
    let mut yaw = 0.0f32;

    for frame in 0..WARMUP + WINDOW {
        let jitter = if cfg.jitter {
            jitter_for(frame)
        } else {
            [0.0, 0.0]
        };
        let cam = Cam::new(eye, yaw, cfg.pitch, w, h, jitter);

        // One stroke, through the path a click takes: aim down the camera's
        // forward ray, place a brush of `BRUSH_RADIUS` at the hit.
        if cfg.edit && frame == WARMUP {
            let field = Dug { brushes: &brushes };
            if let Some(t) = trace(&field, cam.eye, cam.forward) {
                let c = add(cam.eye, mul(cam.forward, t));
                brushes.push([c[0], c[1], c[2], BRUSH_RADIUS]);
            }
        }

        let gb_static = render(&cam, &[], w, h, threads);
        let gb_dig = if brushes.is_empty() {
            gb_static.clone()
        } else {
            render(&cam, &brushes, w, h, threads)
        };
        let measure = frame >= WARMUP;

        // Brush silhouettes in this frame's pixel grid.
        let discs: Vec<[f32; 3]> = brushes
            .iter()
            .filter_map(|b| {
                let c = [b[0], b[1], b[2]];
                cam.project(c)
                    .map(|s| [s[0], s[1], cam.silhouette_radius(c, b[3])])
            })
            .collect();

        if measure {
            let mut changed = 0u64;
            for i in 0..n {
                let (a, b) = (&gb_dig[i], &gb_static[i]);
                let differs = a.hit != b.hit
                    || (0..3).any(|c| (a.world[c] - b.world[c]).abs() > CHANGE_TOL);
                if !differs {
                    continue;
                }
                changed += 1;
                let (x, y) = ((i % w) as f32, (i / w) as f32);
                let inside = discs.iter().any(|d| {
                    let dx = x - d[0];
                    let dy = y - d[1];
                    (dx * dx + dy * dy).sqrt() <= d[2] + SILHOUETTE_SLACK
                });
                if !inside {
                    out.changed_outside += 1;
                }
            }
            out.changed.push(changed);
        }

        for (world, gb) in [(0usize, &gb_static), (1, &gb_dig)] {
            // The reference is this world's deterministic 3-plane shade, shared
            // by all three streams as the MAE target.
            let reference = shade_buffer(arms[0], ar_tex, na_tex, gb, w, h, frame, threads);
            for stream in 0..3 {
                let cur = if stream == 0 {
                    reference.clone()
                } else {
                    shade_buffer(arms[stream], ar_tex, na_tex, gb, w, h, frame, threads)
                };
                let slot = world * 3 + stream;
                let (mut counters, hist_rgb, clipped) = resolve(
                    &taa[slot],
                    &cur,
                    gb,
                    prev_cam.as_ref(),
                    w,
                    h,
                    threads,
                );
                commit(
                    &mut taa[slot],
                    &cur,
                    &reference,
                    gb,
                    &hist_rgb,
                    &clipped,
                    &mut counters,
                    measure,
                );
                if measure {
                    out.frames[world][stream].push(counters);
                }
            }
        }

        prev_cam = Some(cam);
        eye = add(eye, mul(Cam::walk_dir(yaw), cfg.walk * DT));
        yaw += cfg.yaw_rate * DT;
    }

    out.brushes_placed = brushes.len();
    out
}

// ---------------------------------------------------------------------------
// Reporting.
// ---------------------------------------------------------------------------

type Row = Vec<(&'static str, String)>;

fn cpu_mhz() -> String {
    std::fs::read_to_string("/sys/devices/system/cpu/cpu0/cpufreq/scaling_cur_freq")
        .ok()
        .and_then(|s| s.trim().parse::<f64>().ok())
        .map_or_else(|| "unknown".to_string(), |khz| format!("{:.0}", khz / 1000.0))
}

/// `game_dig`'s terrain material as written: three planes, three strata, one
/// filtered fetch each. The baseline C1 and C3 are ratios against.
const ARM_TRIPLANAR: Arm = Arm {
    name: "triplanar_3plane_3stratum_bilinear",
    planes: PlaneMode::Triplanar,
    forced_layer: LAYER_BLEND,
    stochastic_layer: false,
    taps: 1,
    filter: Filter::Bilinear,
    scored: true,
};

/// **C1's arm.** One plane, chosen with probability equal to its triplanar
/// weight.
const ARM_STOCHASTIC_PLANE: Arm = Arm {
    name: "stochastic_1plane_3stratum_bilinear",
    planes: PlaneMode::StochasticOne,
    forced_layer: LAYER_BLEND,
    stochastic_layer: false,
    taps: 1,
    filter: Filter::Bilinear,
    scored: true,
};

/// **C3's arm.** Two planes, Quilez's weighting, no stochastic term.
const ARM_BIPLANAR: Arm = Arm {
    name: "biplanar_2plane_3stratum_bilinear",
    planes: PlaneMode::Biplanar,
    forced_layer: LAYER_BLEND,
    stochastic_layer: false,
    taps: 1,
    filter: Filter::Bilinear,
    scored: true,
};

/// **The product the registration says nobody states**, at the registration's own
/// three taps — and it is 54 fetches, not nine, because the stratum blend is a
/// third factor the registration does not name.
const ARM_STF_3TAP: Arm = Arm {
    name: "stf_3tap_3plane_3stratum_point",
    planes: PlaneMode::Triplanar,
    forced_layer: LAYER_BLEND,
    stochastic_layer: false,
    taps: 3,
    filter: Filter::StochasticPoint,
    scored: false,
};

fn main() {
    if !std::env::args().any(|a| a == "--bench") {
        return;
    }
    let prereg = isomesh::experiment!("P-76");
    let threads = std::thread::available_parallelism().map_or(8, |n| n.get());
    let mhz_start = cpu_mhz();

    // Nine shading arms. Four bilinear stratum-blend arms at k = 2, 6, 12, 18
    // are the fit's design points; the rest price the composition the
    // registration is about, and one is the ALU control.
    let arms = [
        ARM_TRIPLANAR,
        ARM_STOCHASTIC_PLANE,
        ARM_BIPLANAR,
        // One plane and one stratum: 2 fetches, a 9x reduction, and the fit's
        // low design point.
        Arm {
            name: "stochastic_1plane_1stratum_bilinear",
            planes: PlaneMode::StochasticOne,
            forced_layer: LAYER_BLEND,
            stochastic_layer: true,
            taps: 1,
            filter: Filter::Bilinear,
            scored: true,
        },
        ARM_STF_3TAP,
        // One stochastic tap in place of the filter, planes untouched: the same
        // 18 fetch instructions as the baseline and a quarter of the texels.
        Arm {
            name: "stf_1tap_3plane_3stratum_point",
            planes: PlaneMode::Triplanar,
            forced_layer: LAYER_BLEND,
            stochastic_layer: false,
            taps: 1,
            filter: Filter::StochasticPoint,
            scored: false,
        },
        // Both stochastic terms composed the way the row is actually about.
        Arm {
            name: "stochastic_1plane_stf_1tap_point",
            planes: PlaneMode::StochasticOne,
            forced_layer: LAYER_BLEND,
            stochastic_layer: false,
            taps: 1,
            filter: Filter::StochasticPoint,
            scored: false,
        },
        // `game_dig`'s sandbox walls: the forced-layer branch, three planes, one
        // stratum. Context, and the fit's middle design point at k = 6 in a
        // second configuration.
        Arm {
            name: "walls_concrete_3plane_1layer_bilinear",
            planes: PlaneMode::Triplanar,
            forced_layer: LAYER_CONCRETE,
            stochastic_layer: false,
            taps: 1,
            filter: Filter::Bilinear,
            scored: false,
        },
        // The ALU control: every projection, weight, blend and lighting term,
        // eighteen fetch calls, zero texels.
        Arm {
            name: "alu_and_call_only_control",
            planes: PlaneMode::Triplanar,
            forced_layer: LAYER_BLEND,
            stochastic_layer: false,
            taps: 1,
            filter: Filter::NoTexel,
            scored: false,
        },
    ];

    // --- the textures -----------------------------------------------------
    let t_tex = Instant::now();
    let ar_tex = TexArray::new(0x0A1B_2C3D);
    let na_tex = TexArray::new(0x4E5F_6071);
    eprintln!(
        "P-76: two {TEX_SIZE}x{TEX_SIZE}x{TERRAIN_LAYERS} RGBA8 arrays ({:.1} MiB total) in {:.2}s",
        2.0 * (TEX_SIZE * TEX_SIZE * TERRAIN_LAYERS * 4) as f64 / (1024.0 * 1024.0),
        t_tex.elapsed().as_secs_f64()
    );

    // --- stage A: the visibility buffer -----------------------------------
    let t_vis = Instant::now();
    let cam = Cam::new(EYE_START, 0.0, DEMO_PITCH, SCREEN_W, SCREEN_H, [0.0, 0.0]);
    let screen = render(&cam, &[], SCREEN_W, SCREEN_H, threads);
    let frags: Vec<Frag> = screen.iter().copied().filter(|f| f.hit).collect();
    let visible = frags.len();
    eprintln!(
        "P-76: visibility buffer {SCREEN_W}x{SCREEN_H}, {visible} visible fragments \
         ({:.1}%) in {:.2}s",
        100.0 * visible as f64 / (SCREEN_W * SCREEN_H) as f64,
        t_vis.elapsed().as_secs_f64()
    );
    assert!(
        visible > 0,
        "P-76: the visibility buffer is empty, so every cost below would be a zero \
         that could not have been non-zero"
    );

    // --- the exactness control, before any timing -------------------------
    assert_eq!(
        ARM_TRIPLANAR.fetches(),
        18,
        "P-76: `game_dig`'s terrain fragment issues 3 planes x 3 strata x 2 maps = 18 \
         fetches; if this is not 18 the transcription of triplanar.wgsl's LAYER_BLEND \
         branch has drifted from the shader"
    );
    assert_eq!(
        ARM_TRIPLANAR.fetches(),
        3 * ARM_STOCHASTIC_PLANE.fetches(),
        "P-76: C1's 3x is an integer identity -- the plane count factors out of both \
         maps -- and it does not hold, so the arms are not the arms C1 names"
    );
    assert_eq!(
        ARM_STF_3TAP.fetches(),
        54,
        "P-76: three planes x three strata x three taps x two maps is 54, not the \
         registration's nine"
    );

    // --- stage A: the timed sweep -----------------------------------------
    let mut reps: Vec<Vec<f64>> = vec![Vec::new(); arms.len()];
    let mut tallies: Vec<u64> = vec![0; arms.len()];
    let mut means: Vec<[f64; 3]> = vec![[0.0; 3]; arms.len()];
    for rep in 0..SHADE_REPS {
        for (ai, arm) in arms.iter().enumerate() {
            let (ms, tally, mean) = shade_pass(arm, &ar_tex, &na_tex, &frags);
            reps[ai].push(ms);
            tallies[ai] = tally;
            means[ai] = mean;
        }
        eprintln!("P-76: shading rep {} of {SHADE_REPS} done", rep + 1);
    }
    let results: Vec<ShadeResult> = (0..arms.len())
        .map(|ai| {
            let mut sorted = reps[ai].clone();
            sorted.sort_unstable_by(f64::total_cmp);
            ShadeResult {
                fragment_ms: sorted[sorted.len() / 2],
                reps_ms: reps[ai].clone(),
                fetches: tallies[ai],
                mean_rgb: means[ai],
            }
        })
        .collect();

    // The tally is the instrument: it must equal the analytic count exactly, for
    // every arm. An arm that skipped a fetch would otherwise be a fast number.
    for (arm, r) in arms.iter().zip(&results) {
        assert_eq!(
            r.fetches,
            arm.calls() * visible as u64,
            "P-76 {}: tallied {} fetch calls over {visible} fragments, analytic {} per \
             fragment -- the kernel did not issue what the arm declares",
            arm.name,
            r.fetches,
            arm.calls()
        );
    }

    // --- the estimator control --------------------------------------------
    //
    // The contract is that stochastic plane selection is unbiased **in what it
    // selects**: `E[ar] = sum_l lw_l sum_j w_j T(j, l)` and likewise for the
    // pre-normalisation whiteout sum, both of which are linear in the sampled
    // texels. That is a statement about a MEAN, so the control is a mean and not
    // a max: the max over 4096 fragments times seven channels of a Monte-Carlo
    // error is an extreme-value statistic and would need a tolerance loose enough
    // to let a real bug through. The mean over the whole subsample has a standard
    // error of `sigma / (sqrt(draws) * sqrt(fragments))`, about `1.4e-4` here, so
    // [`UNBIASED_BIAS_TOL`] is more than ten sigma while a wrong CDF -- uniform
    // selection instead of weighted, say -- misses by of order the inter-plane
    // albedo spread, `0.05` and up. The worst per-fragment error is reported
    // beside it so a reader can see the variance the mean is hiding.
    //
    // The shaded colour is deliberately NOT part of the contract. See [`Shaded`]:
    // `normalize`, `max(0, .)`, `powf` and the roughness clamp are nonlinear, so
    // filtering after shading leaves a bias no accumulation removes. That bias is
    // measured here and reported as a finding.
    let stride = (visible / UNBIASED_FRAGMENTS).max(1);
    let sampled: Vec<&Frag> = frags
        .iter()
        .step_by(stride)
        .take(UNBIASED_FRAGMENTS)
        .collect();
    let mut bias_ar = [0.0f64; 4];
    let mut bias_nsum = [0.0f64; 3];
    let mut bias_rgb = [0.0f64; 3];
    let mut worst_ar = 0.0f32;
    let mut worst_rgb = 0.0f32;
    let mut tally = 0u64;
    for (i, f) in sampled.iter().enumerate() {
        let mut rng = Rng(1);
        let want = shade(&ARM_TRIPLANAR, &ar_tex, &na_tex, f, &mut rng, &mut tally);
        let mut got_ar = [0.0f64; 4];
        let mut got_nsum = [0.0f64; 3];
        let mut got_rgb = [0.0f64; 3];
        for d in 0..UNBIASED_DRAWS {
            let mut rng = Rng::new(i, d, 7);
            let s = shade(&ARM_STOCHASTIC_PLANE, &ar_tex, &na_tex, f, &mut rng, &mut tally);
            for c in 0..4 {
                got_ar[c] += f64::from(s.ar[c]);
            }
            for c in 0..3 {
                got_nsum[c] += f64::from(s.n_sum[c]);
                got_rgb[c] += f64::from(s.rgb[c]);
            }
        }
        let inv = 1.0 / UNBIASED_DRAWS as f64;
        for c in 0..4 {
            let e = got_ar[c] * inv - f64::from(want.ar[c]);
            bias_ar[c] += e;
            worst_ar = worst_ar.max(e.abs() as f32);
        }
        for c in 0..3 {
            bias_nsum[c] += got_nsum[c] * inv - f64::from(want.n_sum[c]);
            let e = got_rgb[c] * inv - f64::from(want.rgb[c]);
            bias_rgb[c] += e;
            worst_rgb = worst_rgb.max(e.abs() as f32);
        }
    }
    let inv_n = 1.0 / sampled.len() as f64;
    for c in 0..4 {
        bias_ar[c] *= inv_n;
    }
    for c in 0..3 {
        bias_nsum[c] *= inv_n;
        bias_rgb[c] *= inv_n;
    }
    let worst_bias_sampled = bias_ar
        .iter()
        .chain(bias_nsum.iter())
        .fold(0.0f64, |m, b| m.max(b.abs()));
    let worst_bias_shaded = bias_rgb.iter().fold(0.0f64, |m, b| m.max(b.abs()));
    assert!(
        worst_bias_sampled <= f64::from(UNBIASED_BIAS_TOL),
        "P-76: the stochastic single-plane estimator is biased in the quantity it \
         SELECTS -- worst mean signed error {worst_bias_sampled:.6} over \
         {UNBIASED_DRAWS} draws and {} fragments against a tolerance of \
         {UNBIASED_BIAS_TOL} -- so the plane CDF is wrong and every MAE below is a \
         bug rather than variance",
        sampled.len()
    );
    println!(
        "P-76 estimator control: mean signed error in the SAMPLED maps \
         {worst_bias_sampled:.6} <= {UNBIASED_BIAS_TOL} (worst single fragment \
         {worst_ar:.5}); residual bias in the SHADED colour \
         {worst_bias_shaded:.6} (worst single fragment {worst_rgb:.5}) -- \
         filtering after shading is not filtering before shading."
    );

    // --- C1 and C3 --------------------------------------------------------
    let idx = |name: &str| arms.iter().position(|a| a.name == name).expect("arm exists");
    let i_tri = idx(ARM_TRIPLANAR.name);
    let i_stoch = idx(ARM_STOCHASTIC_PLANE.name);
    let i_bipl = idx(ARM_BIPLANAR.name);
    let i_alu = idx("alu_and_call_only_control");

    let ms_tri = results[i_tri].fragment_ms;
    let ms_stoch = results[i_stoch].fragment_ms;
    let ms_bipl = results[i_bipl].fragment_ms;
    let ms_alu = results[i_alu].fragment_ms;

    let c1_fetch_ratio_exact =
        ARM_TRIPLANAR.fetches() == 3 * ARM_STOCHASTIC_PLANE.fetches();
    let c1_cost_ratio = ms_tri / ms_stoch;
    let c1 = c1_fetch_ratio_exact && c1_cost_ratio >= 2.0;
    // The same ratio between each arm's *fastest* observed pass. On a machine with
    // three siblings compiling, the minimum is the least-contended estimate of an
    // uncontended cost, and it is reported beside the median rather than instead
    // of it: if the two disagree about C1's threshold, the clause is noise.
    let min_of = |i: usize| results[i].reps_ms.iter().copied().fold(f64::MAX, f64::min);
    let c1_cost_ratio_min = min_of(i_tri) / min_of(i_stoch);
    let c3_ms_saving_ratio_min =
        (min_of(i_tri) - min_of(i_bipl)) / (min_of(i_tri) - min_of(i_stoch));

    let save_bipl_fetches = ARM_TRIPLANAR.fetches() - ARM_BIPLANAR.fetches();
    let save_stoch_fetches = ARM_TRIPLANAR.fetches() - ARM_STOCHASTIC_PLANE.fetches();
    let c3_fetch_saving_ratio = save_bipl_fetches as f64 / save_stoch_fetches as f64;
    let c3_ms_saving_ratio = (ms_tri - ms_bipl) / (ms_tri - ms_stoch);

    // The fit: `t = F + k*S` over the four bilinear stratum-blend arms.
    let fit_points: Vec<(f64, f64)> = arms
        .iter()
        .zip(&results)
        .filter(|(a, _)| a.filter == Filter::Bilinear && a.forced_layer < 0.0)
        .map(|(a, r)| (a.fetches() as f64, r.fragment_ms))
        .collect();
    let (fit_f, fit_s, fit_r2) = fit(&fit_points);
    let predicted_ratio = (fit_f + 3.0 * fit_s * 6.0) / (fit_f + fit_s * 6.0);
    // The cost of the three *extra* texels a bilinear tap quad reads over a point
    // sample, from the two arms that differ in nothing else.
    //
    // **This came out at or below zero and that is the honest answer, not a
    // defect.** The four bilinear taps are `(x, y)`, `(x+1, y)`, `(x, y+1)`,
    // `(x+1, y+1)` of a 512-wide RGBA8 layer: the first pair shares a 64-byte
    // cache line and the second pair shares another, so a bilinear fetch is two
    // lines where a point sample is one, and on an out-of-order core the extra
    // line is prefetched and free. The point path also pays two `Rng::unit` calls
    // to choose its tap. So **this CPU does not price texels additively and a
    // GPU's texture unit does not price bilinear as four loads either** — which
    // is exactly why `texels_per_fragment` is reported as a count and never
    // converted into a time here.
    let i_point = idx("stf_1tap_3plane_3stratum_point");
    let bilinear_minus_point_ns_per_texel = (results[i_tri].fragment_ms
        - results[i_point].fragment_ms)
        * 1e6
        / (3.0 * 18.0 * visible as f64);
    let call_ns = (results[i_point].fragment_ms - ms_alu) * 1e6 / (18.0 * visible as f64);
    // One integer add at ~4.2 GHz is 0.24 ns; the fit's slope is per fetch.
    let tally_bound_ns = 0.24;
    let fit_slope_ns = fit_s * 1e6 / visible as f64;

    println!();
    println!("P-76 stage A, material shading stage at {SCREEN_W}x{SCREEN_H} ({visible} fragments):");
    for (arm, r) in arms.iter().zip(&results) {
        println!(
            "  {:44} {:3} fetches {:4} texels  {:9.3} ms  {:7.2} ns/frag",
            arm.name,
            arm.fetches(),
            arm.texels(),
            r.fragment_ms,
            r.fragment_ms * 1e6 / visible as f64
        );
    }
    println!("  fit t = F + k*S over k = 2, 6, 12, 18 (bilinear, stratum blend):");
    println!("    F = {fit_f:.4} ms   S = {fit_s:.4} ms/fetch   R^2 = {fit_r2:.6}");
    println!(
        "    6S >= F ? {}  (C1's cost half is reachable iff this; 6S = {:.4} ms)",
        6.0 * fit_s >= fit_f,
        6.0 * fit_s
    );
    println!("    predicted ratio (F + 18S)/(F + 6S) = {predicted_ratio:.4}");
    println!("  C1 fetch ratio  18/6 = 3 exactly: {c1_fetch_ratio_exact}");
    println!("  C1 cost ratio   {c1_cost_ratio:.4}  (median reps; needs >= 2)");
    println!("  C1 cost ratio   {c1_cost_ratio_min:.4}  (fastest reps; needs >= 2)");
    println!("  C3 fetch saving ratio {c3_fetch_saving_ratio:.4}  (needs >= 0.5; exactly 0.5 by construction)");
    println!("  C3 ms saving ratio    {c3_ms_saving_ratio:.4}  (median reps; needs >= 0.5)");
    println!("  C3 ms saving ratio    {c3_ms_saving_ratio_min:.4}  (fastest reps; needs >= 0.5)");

    // --- stage B: the temporal instrument ---------------------------------
    let taa_arms = [
        TaaArm {
            name: HEADLINE_TAA_ARM,
            pitch: -0.6,
            walk: 0.0,
            yaw_rate: 0.0,
            jitter: true,
            edit: true,
        },
        TaaArm {
            name: "dig_at_feet_walk",
            pitch: -0.6,
            walk: WALK_SPEED,
            yaw_rate: YAW_RATE,
            jitter: true,
            edit: true,
        },
        TaaArm {
            name: ZERO_PROOF_ARM,
            pitch: -0.6,
            walk: 0.0,
            yaw_rate: 0.0,
            jitter: false,
            edit: false,
        },
    ];

    struct TaaSummary {
        name: &'static str,
        /// `[world][stream]` window aggregate.
        agg: [[Stream; 3]; 2],
        per_frame: [[Vec<Stream>; 3]; 2],
        changed: Vec<u64>,
        brushes: usize,
    }

    let mut taa_summaries: Vec<TaaSummary> = Vec::new();
    for cfg in &taa_arms {
        let t0 = Instant::now();
        let out = run_taa_arm(
            cfg,
            [&ARM_TRIPLANAR, &ARM_STOCHASTIC_PLANE, &ARM_BIPLANAR],
            &ar_tex,
            &na_tex,
            threads,
        );
        eprintln!(
            "P-76: taa arm {} ({TAA_W}x{TAA_H}, {} frames) in {:.1}s, {} brushes",
            cfg.name,
            WARMUP + WINDOW,
            t0.elapsed().as_secs_f64(),
            out.brushes_placed
        );

        // `P-77`'s control: subtracting a sphere moves the zero set only inside
        // the sphere, so a changed pixel outside every projected silhouette is a
        // wrong projection, a wrong brush position or a leaky carve.
        assert_eq!(
            out.changed_outside, 0,
            "P-76 {}: {} changed pixels outside every brush silhouette",
            cfg.name, out.changed_outside
        );
        for world in 0..2 {
            for stream in 0..3 {
                for (i, f) in out.frames[world][stream].iter().enumerate() {
                    assert!(
                        f.hits > 0,
                        "P-76 {}: world {world} stream {stream} frame +{i} has no hit \
                         pixels, so its MAE is a zero that could not have been non-zero",
                        cfg.name
                    );
                    assert!(
                        f.population > 0,
                        "P-76 {}: world {world} stream {stream} frame +{i} has no history \
                         samples at all",
                        cfg.name
                    );
                }
            }
        }

        let mut agg: [[Stream; 3]; 2] = Default::default();
        for world in 0..2 {
            for stream in 0..3 {
                for f in &out.frames[world][stream] {
                    agg[world][stream].accumulate(f);
                }
            }
        }
        taa_summaries.push(TaaSummary {
            name: cfg.name,
            agg,
            per_frame: out.frames,
            changed: out.changed,
            brushes: out.brushes_placed,
        });
    }

    let head = taa_summaries
        .iter()
        .find(|s| s.name == HEADLINE_TAA_ARM)
        .expect("the headline arm always runs");
    let zero = taa_summaries
        .iter()
        .find(|s| s.name == ZERO_PROOF_ARM)
        .expect("the zero-proof arm always runs");

    // Registered columns, from the headline arm. Stream 0 is the deterministic
    // 3-plane shade, which makes these `P-77`'s own quantity.
    let rej_static = head.agg[0][0].rejected;
    let rej_dig = head.agg[1][0].rejected;
    let mae_static = head.agg[0][1].mae();
    let mae_dig = head.agg[1][1].mae();

    // THE REGISTERED VACUITY CONTROL.
    assert!(
        rej_dig > 0,
        "P-76: VACUITY CONTROL FAILED -- the dig arm rejected zero history samples over \
         the {WINDOW} frames after the edit, so C2 cannot fire and every MAE below is \
         measuring something other than a temporal budget"
    );
    // `M-44`'s other half: the zero was reachable. With the jitter off and the
    // camera frozen the reprojection is the identity and the history converges to
    // the current frame.
    let zp = zero.agg[0][0].rejected;
    assert!(
        (zp as f64) * 100.0 < rej_dig as f64,
        "P-76: the vacuity control is not a control -- with jitter off and the camera \
         frozen the resolve still rejected {zp} history samples against the dig arm's \
         {rej_dig}, so nothing in this fixture could have produced the zero the \
         registered control rules out"
    );
    println!();
    println!(
        "P-76 VACUITY CONTROL: dig arm rejected {rej_dig} history samples \
         (static {rej_static}, frozen-and-unjittered {zp}) -- non-zero, and the zero was reachable."
    );

    let c2 = mae_dig > mae_static;

    // The raw MAE is dominated by TAA's OWN ghosting, which is present in the
    // deterministic 3-plane stream too and has nothing to do with stochastic
    // sampling. Stream 0 resolved against its own unresolved reference IS that
    // floor, so subtracting it isolates the sampling-attributable error -- which
    // is what C2 and C3 are actually claims about. Both the registered raw
    // comparison and this decomposition are reported; `c2_holds` is the raw one,
    // because that is what the registration names.
    let floor_static = head.agg[0][0].mae();
    let floor_dig = head.agg[1][0].mae();
    let excess_static = mae_static - floor_static;
    let excess_dig = mae_dig - floor_dig;
    let excess_ratio = if excess_static > 0.0 {
        excess_dig / excess_static
    } else {
        f64::NAN
    };
    let bipl_static = head.agg[0][2].mae();
    let bipl_dig = head.agg[1][2].mae();
    let bipl_excess_static = bipl_static - floor_static;
    let bipl_excess_dig = bipl_dig - floor_dig;
    // The ideal-conditions floor: frozen camera, no jitter, no dig, 24 frames of
    // warmup. The deterministic stream converges to the reference exactly; what
    // the stochastic stream still carries here is the error no accumulation
    // removes, because an exponential blend at alpha = 0.1 has an effective
    // sample count of `(2 - alpha)/alpha = 19` and a residual bias besides.
    let ideal_reference = zero.agg[0][0].mae();
    let ideal_stochastic = zero.agg[0][1].mae();
    let ideal_biplanar = zero.agg[0][2].mae();
    let zp_stochastic = zero.agg[0][1].rejected;

    // C3's second half: biplanar carries no temporal debt, so its dig-minus-static
    // MAE delta must be below the stochastic arm's.
    let bipl_delta = bipl_dig - bipl_static;
    let stoch_delta = mae_dig - mae_static;
    let c3 = c3_fetch_saving_ratio >= 0.5 && c3_ms_saving_ratio >= 0.5 && bipl_delta < stoch_delta;

    println!();
    println!("P-76 stage B, from the `{HEADLINE_TAA_ARM}` arm at {TAA_W}x{TAA_H}:");
    for s in &taa_summaries {
        for (wi, wn) in [(0usize, "static"), (1, "digging")] {
            let pop = s.agg[wi][0].population.max(1) as f64;
            println!(
                "  {:32} {:8}  rejected det {:7} ({:.5}) stoch {:8}  mae ref {:.6} stoch {:.6} bipl {:.6}",
                s.name,
                wn,
                s.agg[wi][0].rejected,
                s.agg[wi][0].rejected as f64 / pop,
                s.agg[wi][1].rejected,
                s.agg[wi][0].mae(),
                s.agg[wi][1].mae(),
                s.agg[wi][2].mae()
            );
        }
    }
    println!("  C2 as registered: mae_digging {mae_dig:.8} vs mae_static {mae_static:.8}, delta {stoch_delta:+.8}");
    println!("  TAA's own ghosting floor (deterministic stream): static {floor_static:.8} digging {floor_dig:.8}");
    println!("  sampling-attributable excess: static {excess_static:.8} digging {excess_dig:.8}, ratio {excess_ratio:.4}");
    println!("  biplanar excess (a BIAS, not variance): static {bipl_excess_static:.8} digging {bipl_excess_dig:.8}");
    println!("  biplanar delta {bipl_delta:+.8} vs stochastic delta {stoch_delta:+.8} (C3's 'none of C2's cost')");
    println!(
        "  ideal conditions (frozen, unjittered, undug): reference {ideal_reference:.8} \
         stochastic {ideal_stochastic:.8} biplanar {ideal_biplanar:.8}; the stochastic \
         stream rejected {zp_stochastic} samples where the deterministic one rejected {zp}"
    );
    println!("  c1 {c1}  c2 {c2}  c3 {c3}");
    println!();

    // --- rows -------------------------------------------------------------
    let mhz_end = cpu_mhz();
    let common_tail = |row: &mut Row| {
        row.push(("history_rejected_static", rej_static.to_string()));
        row.push(("history_rejected_digging", rej_dig.to_string()));
        row.push(("mae_vs_reference_static", format!("{mae_static:.8}")));
        row.push(("mae_vs_reference_digging", format!("{mae_dig:.8}")));
        row.push(("frames_after_edit", WINDOW.to_string()));
        row.push((
            "biplanar_fetches",
            ARM_BIPLANAR.fetches().to_string(),
        ));
        row.push(("biplanar_fragment_ms", format!("{ms_bipl:.6}")));
        row.push(("c1_holds", c1.to_string()));
        row.push(("c2_holds", c2.to_string()));
        row.push(("c3_holds", c3.to_string()));
        row.push(("mae_taa_floor_static", format!("{floor_static:.8}")));
        row.push(("mae_taa_floor_digging", format!("{floor_dig:.8}")));
        row.push(("mae_sampling_excess_static", format!("{excess_static:.8}")));
        row.push(("mae_sampling_excess_digging", format!("{excess_dig:.8}")));
        row.push((
            "mae_sampling_excess_ratio_dig_over_static",
            if excess_ratio.is_nan() {
                "NA".to_string()
            } else {
                format!("{excess_ratio:.4}")
            },
        ));
        row.push(("mae_biplanar_static", format!("{bipl_static:.8}")));
        row.push(("mae_biplanar_digging", format!("{bipl_dig:.8}")));
        row.push((
            "mae_biplanar_excess_static",
            format!("{bipl_excess_static:.8}"),
        ));
        row.push((
            "mae_biplanar_excess_digging",
            format!("{bipl_excess_dig:.8}"),
        ));
        row.push(("mae_ideal_reference", format!("{ideal_reference:.8}")));
        row.push(("mae_ideal_stochastic", format!("{ideal_stochastic:.8}")));
        row.push(("mae_ideal_biplanar", format!("{ideal_biplanar:.8}")));
        row.push((
            "zero_proof_stochastic_rejected",
            zp_stochastic.to_string(),
        ));
        row.push((
            "rejection_fraction_static",
            format!(
                "{:.8}",
                rej_static as f64 / head.agg[0][0].population.max(1) as f64
            ),
        ));
        row.push((
            "rejection_fraction_digging",
            format!(
                "{:.8}",
                rej_dig as f64 / head.agg[1][0].population.max(1) as f64
            ),
        ));
        row.push(("c1_fetch_ratio_is_exactly_3", c1_fetch_ratio_exact.to_string()));
        row.push(("c1_cost_ratio", format!("{c1_cost_ratio:.4}")));
        row.push(("c1_cost_ratio_min_based", format!("{c1_cost_ratio_min:.4}")));
        row.push((
            "c3_saving_ratio_ms_min_based",
            format!("{c3_ms_saving_ratio_min:.4}"),
        ));
        row.push((
            "c3_saving_ratio_fetches",
            format!("{c3_fetch_saving_ratio:.4}"),
        ));
        row.push(("c3_saving_ratio_ms", format!("{c3_ms_saving_ratio:.4}")));
        row.push(("triplanar_fetches", ARM_TRIPLANAR.fetches().to_string()));
        row.push(("triplanar_fragment_ms", format!("{ms_tri:.6}")));
        row.push((
            "stochastic_plane_fetches",
            ARM_STOCHASTIC_PLANE.fetches().to_string(),
        ));
        row.push(("stochastic_plane_fragment_ms", format!("{ms_stoch:.6}")));
        row.push(("stf_3tap_fetches", ARM_STF_3TAP.fetches().to_string()));
        row.push(("fit_intercept_ms", format!("{fit_f:.6}")));
        row.push(("fit_slope_ms_per_fetch", format!("{fit_s:.6}")));
        row.push(("fit_slope_ns_per_fetch", format!("{fit_slope_ns:.4}")));
        row.push(("fit_r2", format!("{fit_r2:.6}")));
        row.push((
            "c1_cost_reachable_6s_ge_f",
            (6.0 * fit_s >= fit_f).to_string(),
        ));
        row.push(("fit_six_fetch_cost_ms", format!("{:.6}", 6.0 * fit_s)));
        row.push(("fit_predicted_ratio_18_over_6", format!("{predicted_ratio:.4}")));
        row.push((
            "bilinear_minus_point_ns_per_texel",
            format!("{bilinear_minus_point_ns_per_texel:.4}"),
        ));
        row.push(("fetch_call_ns_measured", format!("{call_ns:.4}")));
        row.push(("tally_ns_bound_per_fetch", format!("{tally_bound_ns:.2}")));
        row.push(("alu_control_ms", format!("{ms_alu:.6}")));
        row.push(("screen_pixels", (SCREEN_W * SCREEN_H).to_string()));
        row.push(("visible_fragments", visible.to_string()));
        row.push(("shade_reps", SHADE_REPS.to_string()));
        row.push(("taa_width", TAA_W.to_string()));
        row.push(("taa_height", TAA_H.to_string()));
        row.push(("taa_warmup_frames", WARMUP.to_string()));
        row.push(("threads", threads.to_string()));
        row.push(("cpu_mhz_start", mhz_start.clone()));
        row.push(("cpu_mhz_end", mhz_end.clone()));
        row.push((
            "estimator_mean_signed_error_sampled",
            format!("{worst_bias_sampled:.6}"),
        ));
        row.push((
            "estimator_worst_fragment_error_sampled",
            format!("{worst_ar:.6}"),
        ));
        row.push((
            "shading_bias_mean_signed",
            format!("{worst_bias_shaded:.6}"),
        ));
        row.push((
            "shading_bias_worst_fragment",
            format!("{worst_rgb:.6}"),
        ));
        row.push((
            "estimator_bias_tolerance",
            format!("{UNBIASED_BIAS_TOL:.4}"),
        ));
        row.push(("zero_proof_rejected", zp.to_string()));
        row.push((
            "biplanar_temporal_delta_mae",
            format!("{bipl_delta:+.8}"),
        ));
        row.push((
            "stochastic_temporal_delta_mae",
            format!("{stoch_delta:+.8}"),
        ));
    };

    let mut rows: Vec<Row> = Vec::new();

    for (arm, r) in arms.iter().zip(&results) {
        let ns = r.fragment_ms * 1e6 / visible as f64;
        let mut row: Row = vec![
            ("arm", format!("shading/{}", arm.name)),
            ("fetches_per_fragment", arm.fetches().to_string()),
            ("fragment_ms", format!("{:.6}", r.fragment_ms)),
            // Full-coverage 1080p shading-stage cost. See the module docs: this
            // is a shading cost with a stated denominator, not a frame time.
            (
                "frame_ms",
                format!("{:.6}", ns * (SCREEN_W * SCREEN_H) as f64 / 1e6),
            ),
            ("stage", "material_shading".to_string()),
            ("texels_per_fragment", arm.texels().to_string()),
            ("fetch_calls_per_fragment", arm.calls().to_string()),
            ("planes", arm.planes.count().to_string()),
            ("strata", arm.layer_count().to_string()),
            ("taps", arm.taps.to_string()),
            (
                "filter",
                match arm.filter {
                    Filter::Bilinear => "bilinear",
                    Filter::StochasticPoint => "stochastic_point",
                    Filter::NoTexel => "none",
                }
                .to_string(),
            ),
            ("forced_layer", format!("{:.1}", arm.forced_layer)),
            ("scored_arm", arm.scored.to_string()),
            ("fragment_ns_per_fragment", format!("{ns:.4}")),
            (
                "ns_per_fetch_over_alu_control",
                if arm.calls() == 0 {
                    "NA".to_string()
                } else {
                    format!(
                        "{:.4}",
                        (r.fragment_ms - ms_alu) * 1e6 / (arm.calls() * visible as u64) as f64
                    )
                },
            ),
            (
                "reps_ms_spread",
                format!(
                    "{:.4}",
                    r.reps_ms.iter().copied().fold(f64::MIN, f64::max)
                        - r.reps_ms.iter().copied().fold(f64::MAX, f64::min)
                ),
            ),
            (
                "reps_ms_min",
                format!(
                    "{:.4}",
                    r.reps_ms.iter().copied().fold(f64::MAX, f64::min)
                ),
            ),
            ("fetches_tallied", r.fetches.to_string()),
            ("mean_r", format!("{:.6}", r.mean_rgb[0])),
            ("mean_g", format!("{:.6}", r.mean_rgb[1])),
            ("mean_b", format!("{:.6}", r.mean_rgb[2])),
            ("taa_arm", "NA".to_string()),
            ("frame_offset", "NA".to_string()),
            ("world", "NA".to_string()),
            ("stream", "NA".to_string()),
            ("frame_rejected", "NA".to_string()),
            ("frame_rejection_fraction", "NA".to_string()),
            ("frame_population", "NA".to_string()),
            ("frame_no_history", "NA".to_string()),
            ("frame_mae", "NA".to_string()),
            ("frame_mean_reproj_px", "NA".to_string()),
            ("frame_mean_clip_s", "NA".to_string()),
            ("frame_changed_pixels", "NA".to_string()),
            ("taa_brushes_placed", "NA".to_string()),
        ];
        common_tail(&mut row);
        rows.push(row);
    }

    let stream_name = |s: usize| match s {
        0 => "reference_3plane",
        1 => "stochastic_1plane",
        _ => "biplanar_2plane",
    };
    // Stage-B rows carry the arm under test's own cost columns -- the stochastic
    // single-plane arm -- because that is the configuration C2 evaluates.
    for s in &taa_summaries {
        for world in 0..2 {
            for stream in 0..3 {
                let wn = if world == 0 { "static" } else { "digging" };
                let series = &s.per_frame[world][stream];
                let mut emit = |label: String, f: &Stream, offset: Option<usize>| {
                    let changed = offset
                        .and_then(|o| s.changed.get(o))
                        .map_or_else(|| "NA".to_string(), u64::to_string);
                    let ns = ms_stoch * 1e6 / visible as f64;
                    let mut row: Row = vec![
                        ("arm", label),
                        (
                            "fetches_per_fragment",
                            ARM_STOCHASTIC_PLANE.fetches().to_string(),
                        ),
                        ("fragment_ms", format!("{ms_stoch:.6}")),
                        (
                            "frame_ms",
                            format!("{:.6}", ns * (SCREEN_W * SCREEN_H) as f64 / 1e6),
                        ),
                        ("stage", "taa_resolve".to_string()),
                        (
                            "texels_per_fragment",
                            ARM_STOCHASTIC_PLANE.texels().to_string(),
                        ),
                        (
                            "fetch_calls_per_fragment",
                            ARM_STOCHASTIC_PLANE.calls().to_string(),
                        ),
                        ("planes", ARM_STOCHASTIC_PLANE.planes.count().to_string()),
                        ("strata", ARM_STOCHASTIC_PLANE.layer_count().to_string()),
                        ("taps", ARM_STOCHASTIC_PLANE.taps.to_string()),
                        ("filter", "bilinear".to_string()),
                        ("forced_layer", format!("{LAYER_BLEND:.1}")),
                        ("scored_arm", "true".to_string()),
                        ("fragment_ns_per_fragment", format!("{ns:.4}")),
                        ("ns_per_fetch_over_alu_control", "NA".to_string()),
                        ("reps_ms_spread", "NA".to_string()),
                        ("reps_ms_min", "NA".to_string()),
                        ("fetches_tallied", "NA".to_string()),
                        ("mean_r", "NA".to_string()),
                        ("mean_g", "NA".to_string()),
                        ("mean_b", "NA".to_string()),
                        ("taa_arm", s.name.to_string()),
                        (
                            "frame_offset",
                            offset.map_or_else(|| "NA".to_string(), |o| o.to_string()),
                        ),
                        ("world", wn.to_string()),
                        ("stream", stream_name(stream).to_string()),
                        ("frame_rejected", f.rejected.to_string()),
                        (
                            "frame_rejection_fraction",
                            format!("{:.8}", f.fraction()),
                        ),
                        ("frame_population", f.population.to_string()),
                        ("frame_no_history", f.no_history.to_string()),
                        ("frame_mae", format!("{:.8}", f.mae())),
                        (
                            "frame_mean_reproj_px",
                            format!(
                                "{:.4}",
                                if f.population == 0 {
                                    0.0
                                } else {
                                    f.reproj_px / f.population as f64
                                }
                            ),
                        ),
                        (
                            "frame_mean_clip_s",
                            if f.rejected == 0 {
                                "NA".to_string()
                            } else {
                                format!("{:.6}", f.clip_s / f.rejected as f64)
                            },
                        ),
                        ("frame_changed_pixels", changed),
                        ("taa_brushes_placed", s.brushes.to_string()),
                    ];
                    common_tail(&mut row);
                    rows.push(row);
                };
                emit(
                    format!("taa/{}/{wn}/{}/window", s.name, stream_name(stream)),
                    &s.agg[world][stream],
                    None,
                );
                for (i, f) in series.iter().enumerate() {
                    emit(
                        format!("taa/{}/{wn}/{}/frame_plus_{i}", s.name, stream_name(stream)),
                        f,
                        Some(i),
                    );
                }
            }
        }
    }

    common::experiment::run(prereg, |run| {
        for row in &rows {
            run.record(row);
        }
    });
}
