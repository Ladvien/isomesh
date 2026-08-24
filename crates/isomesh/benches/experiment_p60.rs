//! **P-60 — shifted linear interpolation against the crossing rule the mesher
//! actually uses.**
//!
//! Ticket: R-058. Pre-registered before this harness existed.
//!
//! ```bash
//! cargo bench --bench experiment_p60
//! ```
//!
//! Writes `docs/experiments/p-60.csv`.
//!
//! # What is transcribed, and from where
//!
//! Blu, Thévenaz & Unser, *Linear Interpolation Revitalized*
//! (`10.1109/tip.2004.826093`). Four things are taken from the paper and nothing
//! else is:
//!
//! 1. **The shift.** `tau_opt = (1 − sqrt(3)/3)/2 ≈ 0.21`, which the abstract
//!    calls "close to 1/5". The registration commits to `tau = 1/5` and this
//!    harness uses `1/5`, not `tau_opt`.
//! 2. **The prefilter.** `H_tau(z) = 1 / (1 − tau + tau·z^-1)`, expanded in the
//!    paper as `sum_{k>=0} (−1)^k/(1−tau) · (tau/(1−tau))^k · z^-k`, i.e. the
//!    causal recursion `c_n = −(tau/(1−tau))·c_{n−1} + (1/(1−tau))·f_n` with the
//!    paper's own initialisation `c_0 = f_0`. **At `tau = 1/5` the two constants
//!    are `2^-2` and `1 + 2^-2`**, so the recursion is multiplication-free and
//!    both are exact binary — the only rounding of `1/5` anywhere below is in
//!    the knot offset, which has no exact `f64`.
//! 3. **The reconstruction.** `f_T(x) = sum_n c_n·Λ(x/T − n − tau)` with `Λ` the
//!    hat function, so the knots sit at `(n + tau)·T` and the prefilter is what
//!    enforces `f_T(nT) = f_n`.
//! 4. **The two predictions.** "The gain of shifted over standard linear
//!    interpolation is about 8 dB asymptotically, as the sampling step tends to
//!    0", and Fig. 5's "for `w < 3*pi/4`, the optimally shifted method
//!    outperforms the nonshifted one, by up to 8 dB" — against "when `f(x)` is
//!    the step function … the shifted linear interpolation gives rise to a Gibbs
//!    phenomenon—unlike the standard method".
//!
//! Nothing here reconstructs a proof, and no result of the paper beyond those
//! four statements is used.
//!
//! # What this crate owns rather than the source
//!
//! - **The crossing rule.** `t = a/(a − b)` is `edge_position` in
//!   `marching_cubes/mod.rs`, and it is the baseline every number below is
//!   measured against. It is the crate's own definition of where a surface is,
//!   not the paper's.
//! - **"A sample of exactly zero is outside".** `isomesh::cube` is a private
//!   module, so its strict `< 0` predicate is **transcribed** below rather than
//!   re-decided, exactly as P-54 transcribes Hart's ball. It matters more here
//!   than usual: on this line several reference fields sample **exactly zero**,
//!   and which two samples the bracket names decides what the whole row is
//!   about.
//! - **The line and the grid.** `p(t) = origin + t·e_x` through the domain
//!   centre, `y = z = 0`, `x` across the field's own `domain()`, 65 samples, so
//!   64 cells and `T = 2L/64`. Every reference field's domain is the symmetric
//!   cube `[-L, L]^3`, which is asserted per field rather than assumed.
//! - **The exact root.** Bisection of the field itself, which is available only
//!   because every reference field has an analytic `sample`. It bisects to below
//!   `1e-15` **or** to the last representable interval, whichever comes first —
//!   at `x ≈ -8` one ulp is `1.8e-15`, so a fixed width alone would not
//!   terminate. The floor reached is asserted, not assumed.
//!
//! # Where the shifted root lives, and why that is the delicate half
//!
//! The reconstruction is piecewise linear on the **shifted** grid, not on the
//! sample grid. On `u = x/T ∈ [m − 1 + tau, m + tau]` the only two hats that are
//! non-zero are `n = m − 1` and `n = m`, so with `s = u − (m − 1 + tau)`
//!
//! ```text
//! f_T = c_{m-1}·(1 − s) + c_m·s
//! ```
//!
//! a straight line from `c_{m-1}` to `c_m`. Putting `u = n` gives `s = 1 − tau`
//! and therefore `f_T(nT) = tau·c_{n-1} + (1 − tau)·c_n`, which is exactly the
//! inverse of the recursion: **the interpolation property is an algebraic
//! identity of the filter, for every `n >= 1` and independently of the
//! initialisation.** That is why `interpolation_residual` is a real check on
//! this implementation and not a tautology — it is zero only if the pole, the
//! gain and the knot offset agree with each other.
//!
//! The consequence for a root find is that the sample bracket `[n − 1, n]` which
//! the mesher hands over is **split by a knot** at `n − 1 + tau`, so it is
//! covered by two different linear pieces:
//!
//! ```text
//! u:        n-1        n-1+tau                n
//! piece:  ...  m = n-1  |        m = n        ...
//! values:   c_{n-2}→c_{n-1}   c_{n-1}→c_n
//! ```
//!
//! `f_T` equals `f_{n-1}` at the left end and `f_n` at the right end, which have
//! opposite `is_inside`, so exactly one of the two pieces changes sign — a
//! parity argument, so there is one branch and it is total, with no third case
//! to fall back to. `shifted_segment` records which piece it was.
//!
//! # The closed form the measurement is checked against
//!
//! For a signal that is locally quadratic with second derivative `g''` in cell
//! units, the reconstruction error is exactly
//!
//! ```text
//! f_T − g = (g''/2)·[ tau^2 − tau + s(1 − s) ]
//! ```
//!
//! on every shifted segment — independent of which segment, because both sides
//! are translation-invariant. Standard linear interpolation on the sample
//! bracket gives `(g''/2)·sigma(1 − sigma)` with `sigma` the fractional crossing
//! position. Dividing, and writing `P = tau(1 − tau)`, the **root** errors are
//! in the ratio
//!
//! ```text
//! error_ratio = | q(sigma) − P | / ( sigma·(1 − sigma) )
//!
//! q(sigma) = (sigma − tau)·(1 − sigma + tau)   for sigma >= tau
//!          = (sigma + 1 − tau)·(tau − sigma)   for sigma <  tau
//! ```
//!
//! which at `tau = 1/5` is `|sigma − 2/5| / sigma` above `1/5` and
//! `(sigma + 3/5)/(1 − sigma)` below it. It **never exceeds 1**, is exactly `1`
//! at `sigma = tau` (the root sitting on a knot), and is exactly `0` at
//! `sigma = 2·tau`. It clears the registration's 30% bar on 82% of uniformly
//! distributed crossing positions. That is the mechanism C1 is about, so it is a
//! column — `predicted_error_ratio` — and not a paragraph. It predicts nothing
//! where the local quadratic term is absent, which is exactly what
//! `chord_deviation_cells` reports.
//!
//! # The domain-centre axis is a symmetry axis, and that decides most of the CSV
//!
//! The registration says "an axis-aligned line through each reference field",
//! and the line through the domain centre is the only one that needs no second
//! arbitrary constant. It is also the line every symmetric field is *simplest*
//! on, and the restrictions are derivable in advance rather than discovered:
//!
//! ```text
//! sphere          |x| − 1                        linear on each half
//! box_exact       |x| − 1                        the max(d,0)+min(max d,0) form
//! csg_difference  |x| − 1                        the scoop is at +++, not here
//! thin_plate      |x| − 1  outside the kink      f = −half_y inside it
//! torus           | |x| − 1 | − 0.3              piecewise linear, kinks at 0, ±1
//! noise_cavity    |x| − 1.5 where the cap wins   the cap wins up to its own zero
//! gyroid          max(sin x, |x| − 6)            sin, because cos 0 = 1, sin 0 = 0
//! fbm_terrain     −2·fbm(x, 0, 0)                1-D quintic-fade gradient noise
//! ```
//!
//! Six of the eight are **piecewise linear through their first crossing**, and
//! the closed form above is a statement about a quadratic: with `g'' = 0` both
//! reconstructions are exact and both errors collapse to the bisection floor, so
//! `error_ratio` on those rows is the quotient of two numbers at the resolution
//! of `f64` and is not a measurement of the paper's claim. `standard_is_exact`
//! and `chord_deviation_cells` are the two columns that say so per row, and they
//! are why they exist. Only `gyroid` and `fbm_terrain` carry a curvature the
//! standard rule can actually get wrong.
//!
//! **This is reported, not corrected.** Moving the line off the axis to
//! manufacture curvature would be choosing a fixture by its answer.
//!
//! There is one further exactness worth naming, because it is why so many
//! `guard_band_delta_cells` read exactly zero: on an *exactly linear* signal the
//! causal startup transient `eps_n` cancels out of the left-piece root entirely
//! when the root sits `1 − tau` of a shifted segment along, since the recovered
//! position carries `eps_{n-2}·(1 − 1.25·s)` and `s = 1 − tau = 0.8` kills it.
//! The shifted root is then bit-identical for every window length, which is a
//! true zero rather than an untested one — and `guard_band_truncated` still says
//! whether the window was actually shorter than the line.
//!
//! # Which columns decide which clause
//!
//! - **C1** is decided by `median_error_ratio`, the median of `error_ratio` over
//!   `sphere`, `torus`, `gyroid` and `fbm_terrain`, against `<= 0.70`. It is the
//!   same value on every row on purpose: the clause is a statement about a
//!   population and a reader must not have to recompute it to check the verdict.
//! - **C2** is decided by `median_error_ratio_step_like`, the same statistic over
//!   `box_exact` and `csg_difference`, which is what "worse or equal" is about,
//!   with `gibbs_overshoot` beside it as the stated cause.
//! - **C3** is decided by `guard_band_converged` at `guard_band_k = 10`, on all
//!   eight fields, with `guard_band_delta_cells` at `k ∈ {2, 5, 10, 20}` so the
//!   `(tau/(1−tau))^k = (1/4)^k` decay is visible rather than asserted.
//!
//! `is_step_like` is **measured, not assumed**: the second difference of the
//! sampled line exceeding ten times its own median within three samples of the
//! crossing. `median_second_difference` and
//! `max_second_difference_near_crossing` are both columns so the boolean is
//! auditable — and they have to be, because on a line where the restriction is
//! piecewise linear the median second difference is exactly zero and a
//! multiplicative test against zero is a comparison with nothing.
//!
//! # The guard band cannot always be as long as it asks for
//!
//! `guard_band_k = 20` wants twenty samples before the crossing. On this line
//! several fields cross well inside that: the window is then the whole line, its
//! delta is exactly zero, and `guard_band_truncated` records `false` so the row
//! cannot be read as a convergence that was never tested.
//! `guard_band_window_start` gives the first sample index the recursion actually
//! started from.
//!
//! # Two controls, because a ratio of two small numbers proves nothing
//!
//! Both run before any field does, and both are columns on every row.
//!
//! - **`control_quadratic_error_ratio`** is the positive control: a synthetic
//!   exactly-quadratic line on the same 65-sample grid, with its root at a known
//!   `sigma = 1/2` and no bisection involved, where the closed form above says
//!   the ratio is `1/5`. If the pole sign, the gain, the knot offset or the
//!   segment choice were wrong this would not be `0.2`, and every field number
//!   would be measuring the mistake.
//! - **`control_sabotage_error_ratio`** is the negative control, in P-54's shape:
//!   the identical code path fed the **raw samples in place of `c`** — the shift
//!   applied without the prefilter that makes it interpolate, which is the single
//!   most likely way to get this wrong. It biases the root by about `tau` of a
//!   cell and must come out far worse than the standard rule; a value near `0.2`
//!   there would mean the prefilter is doing nothing and the positive control is
//!   passing for the wrong reason.
//!
//! # No extractor is called
//!
//! This reads a 1-D sample line, bisects a field, and filters 65 numbers. It
//! cannot move a golden hash, and nothing in it is reachable from the shipped
//! path.

mod common;

use common::experiment::Run;
use isomesh::fields::ReferenceField;
use isomesh::{Sdf, for_each_reference_field};

/// **A sample of exactly zero is outside.**
///
/// Transcribed from `isomesh::cube::is_inside`, which is not reachable from a
/// bench because `cube` is a private module. Strict `< 0`, verbatim: the crate's
/// own doc records that Lengyel's dissertation §3.1.1 requires the choice be
/// made once and applied everywhere, and that strictness is what makes `a − b`
/// non-zero on a cut edge so `t = a/(a − b)` needs no epsilon guard.
fn is_inside(value: f64) -> bool {
    value < 0.0
}

/// The registered sample count along the line.
const SAMPLES: usize = 65;

/// Cells along the line, so `T = 2L / CELLS`.
const CELLS: usize = SAMPLES - 1;

/// The registered shift, `1/5`.
///
/// The nearest `f64` to `1/5`, which is where the only rounding of the shift
/// lives: [`POLE`] and [`GAIN`] are exact binary and the recursion below is not
/// affected by this constant at all.
const TAU: f64 = 0.2;

/// `tau / (1 − tau)` at `tau = 1/5`: `2^-2`, exact.
const POLE: f64 = 0.25;

/// `1 / (1 − tau)` at `tau = 1/5`: `1 + 2^-2`, exact.
const GAIN: f64 = 1.25;

/// `tau·(1 − tau)`, the constant in the closed-form error ratio.
const TAU_PRODUCT: f64 = TAU * (1.0 - TAU);

/// The registered guard-band window lengths.
const GUARD_K: [usize; 4] = [2, 5, 10, 20];

/// The registered convergence threshold for C3, in cells.
const GUARD_TOLERANCE: f64 = 1e-6;

/// Requested bisection width for the exact root, in field units.
const BISECT_WIDTH: f64 = 1e-15;

/// Points used to sample the crossing interval for `gibbs_overshoot` and
/// `chord_deviation_cells`. The registration asks for at least 1,024.
const DENSE: usize = 1025;

/// A second difference this many times the median counts as step-like.
const STEP_LIKE_FACTOR: f64 = 10.0;

/// "Within three samples of the crossing", in samples.
const STEP_LIKE_RADIUS: usize = 3;

/// Below this many cells the standard rule is at its own numerical floor and a
/// ratio against it is arithmetic rather than measurement.
const EXACT_FLOOR_CELLS: f64 = 1e-12;

/// The four fields C1 is a statement about.
const SMOOTH: [&str; 4] = ["sphere", "torus", "gyroid", "fbm_terrain"];

/// The two fields C2 is a statement about.
const STEP_LIKE: [&str; 2] = ["box_exact", "csg_difference"];

/// Where the positive control's root sits, in cells from the line start.
///
/// A half-integer, so `sigma = 1/2` and the closed form predicts exactly `1/5`.
const CONTROL_ROOT: f64 = 31.5;

/// Half the positive control's second derivative in cell units.
///
/// Small enough that the quadratic's *other* root, at
/// `CONTROL_ROOT − 1/CONTROL_CURVATURE = −18.5`, is off the line, so the first
/// sign change is the intended one.
const CONTROL_CURVATURE: f64 = 0.02;

/// What the closed form says `control_quadratic_error_ratio` is.
const CONTROL_EXPECTED: f64 = 0.2;

/// Tolerance on the positive control. Loose enough to absorb the `O(curvature²)`
/// terms the closed form drops, tight enough that a sign or offset error cannot
/// hide: the sabotage below lands two orders of magnitude away.
const CONTROL_TOLERANCE: f64 = 1e-2;

/// The negative control has to be *visibly* worse than the standard rule.
const SABOTAGE_MIN_RATIO: f64 = 2.0;

/// The `c_n` of the causal prefilter over the samples `start..=end`.
///
/// Index `i` of the result is sample `start + i`. `c_start = f_start` is the
/// paper's `c_0 = f_0`, applied at whichever sample the window begins.
fn prefilter(f: &[f64], start: usize, end: usize, pole: f64) -> Vec<f64> {
    assert!(
        start < end && end < f.len(),
        "prefilter window out of range"
    );
    let mut c = Vec::with_capacity(end - start + 1);
    c.push(f[start]);
    for value in &f[(start + 1)..=end] {
        let prev = c[c.len() - 1];
        c.push(GAIN * value - pole * prev);
    }
    c
}

/// The shifted reconstruction at `u = x/T`, from coefficients based at `start`.
fn reconstruct(c: &[f64], start: usize, u: f64) -> f64 {
    // The segment is the `m` with `m − 1 + tau <= u <= m + tau`.
    let m = (u - TAU).floor() + 1.0;
    let hi = m as usize;
    assert!(
        m >= 1.0 && hi > start && hi - start < c.len(),
        "shifted reconstruction at u = {u} needs a knot outside the filtered window"
    );
    let s = u - (m - 1.0 + TAU);
    c[hi - 1 - start] * (1.0 - s) + c[hi - start] * s
}

/// The zero of the shifted reconstruction inside the sample bracket
/// `[n − 1, n]`, in cells from the line start, and which piece held it.
///
/// Feeding raw samples as `c` is the negative control, so this takes `c` and
/// never looks at `f`: the endpoint values are computed from the reconstruction
/// itself, which is what makes the sabotage traverse the same code.
fn shifted_root(c: &[f64], start: usize, n: usize, label: &str) -> (f64, &'static str) {
    assert!(
        start + 2 <= n,
        "{label}: the shifted root in bracket [{}, {n}] needs c_{}, which is \
         before the window start {start}",
        n - 1,
        n - 2
    );
    let c2 = c[n - 2 - start];
    let c1 = c[n - 1 - start];
    let c0 = c[n - start];
    let v_left = TAU * c2 + (1.0 - TAU) * c1;
    let v_right = TAU * c1 + (1.0 - TAU) * c0;
    assert!(
        is_inside(v_left) != is_inside(v_right),
        "{label}: the shifted reconstruction does not change sign across the \
         bracket [{}, {n}] that the field changes sign across ({v_left:e} to \
         {v_right:e})",
        n - 1
    );
    // Exactly one of the two pieces changes sign, by parity, so this is total.
    if is_inside(v_left) == is_inside(c1) {
        ((n - 1) as f64 + TAU + c1 / (c1 - c0), "right")
    } else {
        ((n - 2) as f64 + TAU + c2 / (c2 - c1), "left")
    }
}

/// The closed-form `error_ratio` for a locally quadratic signal crossing at
/// fractional position `sigma`.
fn predicted_ratio(sigma: f64) -> f64 {
    let q = if sigma >= TAU {
        (sigma - TAU) * (1.0 - sigma + TAU)
    } else {
        (sigma + 1.0 - TAU) * (TAU - sigma)
    };
    (q - TAU_PRODUCT).abs() / (sigma * (1.0 - sigma))
}

/// Bisect `field` along `y = z = 0` between two abscissae of opposite side.
///
/// Terminates at [`BISECT_WIDTH`] or at the last representable interval,
/// whichever comes first; the caller asserts which.
fn bisect<F>(field: &F, mut lo: f64, mut hi: f64) -> (f64, f64)
where
    F: Sdf<Scalar = f64>,
{
    let lo_inside = is_inside(field.sample([lo, 0.0, 0.0]));
    loop {
        let mid = f64::midpoint(lo, hi);
        if mid <= lo || mid >= hi || hi - lo <= BISECT_WIDTH {
            break;
        }
        if is_inside(field.sample([mid, 0.0, 0.0])) == lo_inside {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    (f64::midpoint(lo, hi), hi - lo)
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

/// One guard-band window's answer.
struct Guard {
    /// Samples the window asked for.
    k: usize,
    /// First sample index the recursion actually started from.
    start: usize,
    /// Whether the window really was shorter than the line.
    truncated: bool,
    /// `|truncated root − full-line root|`, in cells.
    delta: f64,
    /// Whether that is below [`GUARD_TOLERANCE`].
    converged: bool,
}

/// Everything measured on one field's line, before the population medians exist.
struct Line {
    name: &'static str,
    /// Line start, in field units.
    x0: f64,
    /// `T`, in field units.
    cell: f64,
    /// The bracket holding the first crossing is `[n − 1, n]`.
    n: usize,
    u_exact: f64,
    u_standard: f64,
    u_shifted: f64,
    error_standard: f64,
    error_shifted: f64,
    ratio: f64,
    /// `sigma`, the standard rule's fractional position in the bracket.
    fraction: f64,
    predicted: f64,
    chord_deviation: f64,
    gibbs: f64,
    median_d2: f64,
    max_d2_near: f64,
    step_like: bool,
    interpolation_residual: f64,
    bisect_floor: f64,
    segment: &'static str,
    guards: [Guard; 4],
}

impl Line {
    /// The exact root, in field units.
    fn exact_x(&self) -> f64 {
        self.x0 + self.cell * self.u_exact
    }

    /// The standard rule's root, in field units.
    fn standard_x(&self) -> f64 {
        self.x0 + self.cell * self.u_standard
    }

    /// The shifted reconstruction's root, in field units.
    fn shifted_x(&self) -> f64 {
        self.x0 + self.cell * self.u_shifted
    }
}

/// Measure one reference field's centre line.
fn measure<F>(name: &'static str, field: &F) -> Line
where
    F: ReferenceField + Sdf<Scalar = f64>,
{
    let (lo, hi) = field.domain();
    assert!(
        (lo[1] + hi[1]).abs() < 1e-12 && (lo[2] + hi[2]).abs() < 1e-12,
        "{name}: domain is not centred on the origin in y and z, so y = z = 0 is \
         not the domain centre"
    );
    let x0 = lo[0];
    let cell = (hi[0] - x0) / CELLS as f64;
    let f: Vec<f64> = (0..SAMPLES)
        .map(|i| field.sample([x0 + cell * i as f64, 0.0, 0.0]))
        .collect();

    let n = (1..SAMPLES)
        .find(|&i| is_inside(f[i - 1]) != is_inside(f[i]))
        .unwrap_or_else(|| {
            panic!(
                "{name}: no sign change on the line y = z = 0 across [{x0}, \
                 {}] at {SAMPLES} samples, so this field has no root to compare \
                 three reconstructions on",
                hi[0]
            )
        });
    assert!(
        n >= 2,
        "{name}: the first crossing is in bracket [{}, {n}], which has no \
         sample before it for the shifted reconstruction's left knot",
        n - 1
    );

    // ── the three roots ────────────────────────────────────────────────────
    let a = f[n - 1];
    let b = f[n];
    let fraction = a / (a - b);
    let u_standard = (n - 1) as f64 + fraction;

    let x_lo = x0 + cell * (n - 1) as f64;
    let x_hi = x0 + cell * n as f64;
    let (x_exact, bisect_floor) = bisect(field, x_lo, x_hi);
    assert!(
        bisect_floor <= BISECT_WIDTH.max(4.0 * f64::EPSILON * x_exact.abs()),
        "{name}: bisection stopped at width {bisect_floor:e}, which is neither \
         below {BISECT_WIDTH:e} nor at the f64 resolution of {x_exact}"
    );
    let u_exact = (x_exact - x0) / cell;

    let c = prefilter(&f, 0, SAMPLES - 1, POLE);
    let (u_shifted, segment) = shifted_root(&c, 0, n, name);

    let error_standard = (u_standard - u_exact).abs();
    let error_shifted = (u_shifted - u_exact).abs();

    // ── the control on the reconstruction itself ───────────────────────────
    let mut interpolation_residual = 0.0f64;
    for i in 1..SAMPLES {
        let residual = (TAU * c[i - 1] + (1.0 - TAU) * c[i] - f[i]).abs();
        interpolation_residual = interpolation_residual.max(residual);
    }

    // ── Gibbs, and the nonlinearity the standard rule has to fight ─────────
    let mut gibbs = 0.0f64;
    let mut deviation = 0.0f64;
    for j in 0..DENSE {
        let s = j as f64 / (DENSE - 1) as f64;
        let u = (n - 1) as f64 + s;
        let x = x_lo + cell * s;
        let truth = field.sample([x, 0.0, 0.0]);
        gibbs = gibbs.max((reconstruct(&c, 0, u) - truth).abs());
        deviation = deviation.max((truth - (a + (b - a) * s)).abs());
    }
    let chord_deviation = deviation / (a - b).abs();

    // ── is_step_like, as a measurement ────────────────────────────────────
    let d2: Vec<f64> = (1..SAMPLES - 1)
        .map(|i| (f[i - 1] - 2.0 * f[i] + f[i + 1]).abs())
        .collect();
    let median_d2 = median(&d2);
    let first = (n - 1).saturating_sub(STEP_LIKE_RADIUS).max(1);
    let last = (n + STEP_LIKE_RADIUS).min(SAMPLES - 2);
    let max_d2_near = (first..=last).fold(0.0f64, |acc, i| acc.max(d2[i - 1]));
    let step_like = max_d2_near > STEP_LIKE_FACTOR * median_d2;

    // ── C3 ────────────────────────────────────────────────────────────────
    let guards = GUARD_K.map(|k| {
        let start = n.saturating_sub(k);
        let window = prefilter(&f, start, n, POLE);
        let (u_truncated, _) = shifted_root(&window, start, n, name);
        let delta = (u_truncated - u_shifted).abs();
        Guard {
            k,
            start,
            truncated: start > 0,
            delta,
            converged: delta < GUARD_TOLERANCE,
        }
    });

    Line {
        name,
        x0,
        cell,
        n,
        u_exact,
        u_standard,
        u_shifted,
        error_standard,
        error_shifted,
        ratio: error_shifted / error_standard,
        fraction,
        predicted: predicted_ratio(fraction),
        chord_deviation,
        gibbs,
        median_d2,
        max_d2_near,
        step_like,
        interpolation_residual,
        bisect_floor,
        segment,
        guards,
    }
}

/// The two instrument controls, on a synthetic line with a root in closed form.
struct Control {
    /// The shifted rule against the standard rule on an exact quadratic.
    ratio: f64,
    /// The same, with the prefilter omitted entirely.
    sabotage: f64,
    standard: f64,
    shifted: f64,
    sabotage_error: f64,
}

/// Run the positive and negative controls.
fn control() -> Control {
    let g: Vec<f64> = (0..SAMPLES)
        .map(|i| {
            let d = i as f64 - CONTROL_ROOT;
            d + CONTROL_CURVATURE * d * d
        })
        .collect();
    let n = (1..SAMPLES)
        .find(|&i| is_inside(g[i - 1]) != is_inside(g[i]))
        .expect("the quadratic control crosses zero by construction");
    assert_eq!(
        n,
        32,
        "the control's first crossing moved to bracket [{}, {n}]; its predicted \
         ratio is only 1/5 for a root at sigma = 1/2 in [31, 32]",
        n - 1
    );

    let a = g[n - 1];
    let b = g[n];
    let u_standard = (n - 1) as f64 + a / (a - b);
    let c = prefilter(&g, 0, SAMPLES - 1, POLE);
    let (u_shifted, _) = shifted_root(&c, 0, n, "control");
    // The sabotage: the shift without the prefilter that makes it interpolate.
    let (u_sabotage, _) = shifted_root(&g, 0, n, "control/sabotage");

    let standard = (u_standard - CONTROL_ROOT).abs();
    let shifted = (u_shifted - CONTROL_ROOT).abs();
    let sabotage_error = (u_sabotage - CONTROL_ROOT).abs();
    Control {
        ratio: shifted / standard,
        sabotage: sabotage_error / standard,
        standard,
        shifted,
        sabotage_error,
    }
}

/// `{:.6e}`, which is enough digits to decide every clause here and keeps `inf`
/// and `NaN` visible rather than rounded into a number.
fn sci(v: f64) -> String {
    format!("{v:.6e}")
}

fn main() {
    if !std::env::args().any(|arg| arg == "--bench") {
        return;
    }
    let prereg = isomesh::experiment!("P-60");

    common::experiment::run(prereg, |run: &mut Run| {
        // ── the instrument, before any field ──────────────────────────────
        let ctl = control();
        println!("instrument control — exact quadratic, root at sigma = 1/2, 65 samples");
        println!(
            "  standard error {:.6e} cells   shifted {:.6e}   ratio {:.9} \
             (closed form {CONTROL_EXPECTED})",
            ctl.standard, ctl.shifted, ctl.ratio
        );
        println!(
            "  sabotage (no prefilter) error {:.6e} cells   ratio {:.6e}",
            ctl.sabotage_error, ctl.sabotage
        );
        assert!(
            (ctl.ratio - CONTROL_EXPECTED).abs() < CONTROL_TOLERANCE,
            "the shifted reconstruction is wrong: on an exact quadratic with \
             sigma = 1/2 the closed form gives {CONTROL_EXPECTED} and this \
             harness measured {:.9}. Every field number below would be \
             measuring that mistake.",
            ctl.ratio
        );
        assert!(
            ctl.sabotage > SABOTAGE_MIN_RATIO,
            "the negative control did not fire: dropping the prefilter entirely \
             gave ratio {:.6e}, which is not worse than the standard rule, so \
             the positive control proves nothing about the prefilter",
            ctl.sabotage
        );
        println!("  both controls pass\n");

        // ── the sweep ─────────────────────────────────────────────────────
        let mut lines: Vec<Line> = Vec::new();
        for_each_reference_field!(f64, |name, field| {
            lines.push(measure(name, &field));
        });
        assert_eq!(
            lines.len(),
            8,
            "the field sweep stopped early, which is M-199's shape"
        );

        println!(
            "field           n   T          sigma      exact root        \
             err_std      err_shift    ratio        pred    chord    seg"
        );
        for line in &lines {
            println!(
                "{:<14} {:>3}  {:.7}  {:.6}  {:>16.12}  {:.4e}  {:.4e}  \
                 {:>10}  {:>6}  {:.2e}  {}",
                line.name,
                line.n,
                line.cell,
                line.fraction,
                line.exact_x(),
                line.error_standard,
                line.error_shifted,
                format!("{:.4e}", line.ratio),
                format!("{:.3}", line.predicted),
                line.chord_deviation,
                line.segment,
            );
        }
        println!();
        println!(
            "field           median|d2|   max|d2| near   step_like  \
             gibbs        interp resid  bisect floor"
        );
        for line in &lines {
            println!(
                "{:<14} {:.4e}   {:.4e}      {:<9}  {:.4e}   {:.4e}     {:.4e}",
                line.name,
                line.median_d2,
                line.max_d2_near,
                line.step_like,
                line.gibbs,
                line.interpolation_residual,
                line.bisect_floor,
            );
            assert!(
                line.interpolation_residual < 1e-12,
                "{}: the shifted reconstruction does not reproduce the samples \
                 it was filtered from (worst residual {:.6e}); the pole, the \
                 gain and the knot offset disagree",
                line.name,
                line.interpolation_residual
            );
        }
        println!();
        println!("field           k    window start  truncated  delta cells  converged");
        for line in &lines {
            for guard in &line.guards {
                println!(
                    "{:<14} {:>2}   {:>11}  {:<9}  {:.6e}  {}",
                    line.name, guard.k, guard.start, guard.truncated, guard.delta, guard.converged,
                );
            }
        }
        println!();

        // ── the two population statistics ─────────────────────────────────
        let smooth: Vec<f64> = lines
            .iter()
            .filter(|l| SMOOTH.contains(&l.name))
            .map(|l| l.ratio)
            .collect();
        assert_eq!(smooth.len(), SMOOTH.len(), "C1's population is incomplete");
        let median_smooth = median(&smooth);

        let stepwise: Vec<f64> = lines
            .iter()
            .filter(|l| STEP_LIKE.contains(&l.name))
            .map(|l| l.ratio)
            .collect();
        assert_eq!(
            stepwise.len(),
            STEP_LIKE.len(),
            "C2's population is incomplete"
        );
        let median_stepwise = median(&stepwise);

        let show = |v: &[f64]| {
            v.iter()
                .map(|r| format!("{r:.6e}"))
                .collect::<Vec<_>>()
                .join(", ")
        };
        println!(
            "C1 population {} ratios [{}] → median {median_smooth:.6e}",
            SMOOTH.join(", "),
            show(&smooth)
        );
        println!(
            "C2 population {} ratios [{}] → median {median_stepwise:.6e}",
            STEP_LIKE.join(", "),
            show(&stepwise)
        );
        println!(
            "rows where the standard rule is already at its numerical floor \
             (< {EXACT_FLOOR_CELLS:e} cells): {}",
            lines
                .iter()
                .filter(|l| l.error_standard < EXACT_FLOOR_CELLS)
                .map(|l| l.name)
                .collect::<Vec<_>>()
                .join(", ")
        );
        println!();

        // ── rows ──────────────────────────────────────────────────────────
        for line in &lines {
            for guard in &line.guards {
                run.record(&[
                    ("field", line.name.to_string()),
                    ("samples", SAMPLES.to_string()),
                    ("tau", format!("{TAU:.9}")),
                    ("root_error_standard", sci(line.error_standard)),
                    ("root_error_shifted", sci(line.error_shifted)),
                    ("error_ratio", sci(line.ratio)),
                    ("median_error_ratio", sci(median_smooth)),
                    ("is_step_like", line.step_like.to_string()),
                    ("gibbs_overshoot", sci(line.gibbs)),
                    ("guard_band_k", guard.k.to_string()),
                    ("guard_band_delta_cells", sci(guard.delta)),
                    ("guard_band_converged", guard.converged.to_string()),
                    // ── extras ──
                    ("median_error_ratio_step_like", sci(median_stepwise)),
                    ("crossing_index", line.n.to_string()),
                    ("cell_size", format!("{:.9}", line.cell)),
                    ("crossing_fraction", format!("{:.9}", line.fraction)),
                    ("exact_root", format!("{:.15}", line.exact_x())),
                    ("root_standard", format!("{:.15}", line.standard_x())),
                    ("root_shifted", format!("{:.15}", line.shifted_x())),
                    ("predicted_error_ratio", sci(line.predicted)),
                    ("chord_deviation_cells", sci(line.chord_deviation)),
                    (
                        "standard_is_exact",
                        (line.error_standard < EXACT_FLOOR_CELLS).to_string(),
                    ),
                    ("median_second_difference", sci(line.median_d2)),
                    ("max_second_difference_near_crossing", sci(line.max_d2_near)),
                    ("interpolation_residual", sci(line.interpolation_residual)),
                    ("shifted_segment", line.segment.to_string()),
                    ("guard_band_window_start", guard.start.to_string()),
                    ("guard_band_truncated", guard.truncated.to_string()),
                    ("control_quadratic_error_ratio", sci(ctl.ratio)),
                    ("control_sabotage_error_ratio", sci(ctl.sabotage)),
                ]);
            }
        }
    });
}
