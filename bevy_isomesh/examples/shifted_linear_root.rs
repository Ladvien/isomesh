//! E-316 — the shifted-linear gain is real, and where the root falls decides
//! whether you get any of it.
//!
//! ```bash
//! cd bevy_isomesh && cargo run --example shifted_linear_root --release
//! ```
//!
//! **Always `--release`.** Nothing here meshes, so a debug build is merely slow
//! rather than unusable — but the sweep re-bisects a field every frame and a
//! debug build makes that visible.
//!
//! It runs itself and loops. `Space` freezes it; `Left`/`Right` drive `sigma` by
//! hand and `R` hands it back to the loop; `Up`/`Down` (or `1`–`9`) change which
//! signal the line is cut through. `ISOMESH_FIELD=0..8` pins one without a
//! keyboard — `0` is the exactly-quadratic control the sweep runs on, `1`–`8`
//! are the eight reference fields on their domain-centre axis, in the order
//! `docs/experiments/p-60.csv` lists them.
//!
//! ```bash
//! DISPLAY=:77 ISOMESH_CAPTURE_FRAMES=110 ISOMESH_CAPTURE_EVERY=2 FPS=14 \
//!   ./scripts/record_gif.sh shifted_linear_root \
//!   docs/gifs/where-the-root-falls-decides-the-gain.gif
//! ```
//!
//! Demonstrates **M-359 / ✗42 (P-60, R-058)** — `docs/experiments/p-60.csv`,
//! 32 rows — as an instrument rather than as a table.
//!
//! # The finding, in one paragraph
//!
//! Every Marching Cubes vertex is `t = a/(a − b)` between two corner samples.
//! Blu, Thévenaz & Unser (`10.1109/tip.2004.826093`) show that shifting the
//! sampling knots by a fixed, signal-independent `tau` and enforcing the
//! interpolation property recovers *"about 8 dB asymptotically"* on the
//! **reconstruction**. P-60 asked whether that transfers to the **root
//! position** a mesher actually uses. It does, and the answer is a closed form
//! rather than a statistic: the two root errors stand in the ratio
//!
//! ```text
//! error_ratio = | sigma - 2 tau | / sigma          for sigma >= tau
//!             = ( sigma + 1 - tau ) / ( 1 - sigma ) for sigma <  tau
//! ```
//!
//! for a root at fraction `sigma` of the cell. At `tau = 1/5` that is **exactly
//! 1 when the root lands on a knot** (`sigma = tau`), **exactly 0 at
//! `sigma = 2 tau = 0.4`**, never above `1` in the quadratic regime, and it
//! clears the registration's 30% bar on `14/17 = 82.35%` of uniformly
//! distributed crossing positions.
//!
//! So C1's registered *"at least 30% lower"* is **FALSIFIED** at a median ratio
//! of `1.486067` — not because the method fails, but because a median over four
//! arbitrary crossing positions is a lottery over `sigma`.
//!
//! # What a studio meets, and why this is not an accuracy win
//!
//! **A mesher cannot choose `sigma`.** The surface falls where it falls; a chunk
//! of terrain hands the extractor sixty-four thousand crossings at sixty-four
//! thousand different fractions. The expected benefit is therefore the **82%**
//! figure and not the 8 dB one, and it is bought at the price of a causal IIR
//! prefilter across the whole grid line. `edge_position` in
//! `marching_cubes/mod.rs` is **unchanged** and this example does not propose
//! changing it. It is on screen so that the next person to read the paper does
//! not have to re-derive why.
//!
//! What *is* shippable is the guard band. The prefilter is a one-pole with
//! `(tau/(1−tau))^k = (1/4)^k` decay, and truncating it to the `k` samples
//! preceding the crossing moves the recovered root by under `1e-6` cells at
//! **`k = 10`** on 8 of 8 fields — while **`k = 5` fails** on both fields that
//! carry a non-zero delta. Ten samples of slab overlap is a number a chunked
//! mesher can price, and `10` is the real threshold rather than a loose one.
//!
//! # What the demo does about it
//!
//! One grid line, three panels, and the numbers that decide each of them.
//!
//! - **THE LINE** shows the 65-sample line around its crossing: the samples
//!   `f_n`, the standard piecewise-linear reconstruction through them, the
//!   prefilter coefficients `c_n` sitting on the **shifted** knot grid at
//!   `(n + tau)·T`, and the exact zero from bisecting the field. The three roots
//!   are a thousandth of a cell apart, so a magnified inset beside it — whose
//!   half-width is printed, never implied — is where they separate.
//! - **THE RATIO** draws the closed form as a curve and puts the live
//!   measurement on top of it. `sigma` is swept by **sliding the field along the
//!   line**, never by moving the grid, so the sample positions stay fixed the
//!   way a mesher's are. The trace nose-dives to zero at `sigma = 2 tau` and
//!   touches `1` at `sigma = tau`, which is the whole finding in one gesture.
//! - **THE GUARD BAND** recomputes the shifted root from only the `k` preceding
//!   samples for `k` in `2, 5, 10, 20` and plots the delta against the
//!   whole-line answer on a log axis, with the `1e-6` bar drawn across it.
//!   `k = 5` is visibly above it and `k = 10` visibly below.
//!
//! # Two controls, because a ratio of two small numbers proves nothing
//!
//! Both are the harness's own, run here on the same code path, and both are
//! printed:
//!
//! - The **positive** control is a synthetic exactly-quadratic line whose
//!   closed-form ratio at its `sigma = 1/2` crossing is exactly `1/5`. It
//!   measures **`1.984127e-1`**, which pins the pole sign, the gain, the knot
//!   offset and the two-piece segment choice at once. It is also the sweep's
//!   own signal, so the clip's opening frame *is* the committed control row.
//! - The **negative** control feeds the identical code path the **raw samples in
//!   place of `c`** — the shift applied without the prefilter that makes it
//!   interpolate, which is the most likely way to get this wrong. It comes back
//!   **39x worse** than the standard rule, a root biased by `0.195` cells, i.e.
//!   `tau`. It is drawn on the wide panel, a visible fifth of a cell to the
//!   side.
//!
//! # Half of C1's population was degenerate, and the demo says so
//!
//! The domain-centre axis is a **symmetry axis** of six of the eight reference
//! fields. `sphere` restricts to `|x| − 1` there, `torus` to `||x| − 1| − 0.3`;
//! both are piecewise linear through the first crossing, both reconstructions
//! are exact on a straight line, and both errors are the **bisection floor**
//! `7.105427e-15`. Their `error_ratio` is a quotient of two numbers at the
//! resolution of `f64` and is not a measurement of the paper's claim. Only
//! `gyroid` and `fbm_terrain` carry a curvature the standard rule can get wrong,
//! and they **split**: `fbm_terrain` is 4.9x better, `gyroid` 2x worse because
//! its crossing sits at `x = −pi` exactly, an inflection of `sin` where the local
//! second derivative vanishes and the error is cubic-dominated.
//!
//! That is why the sweep runs on the quadratic control and not on a reference
//! field: the closed form is a statement about a **quadratic**, and sliding
//! `gyroid` along its line would sweep `sigma` through a shape the formula does
//! not describe. Selecting a reference field pins it at the committed
//! configuration — slide zero — and prints its CSV row; the arrows still slide
//! it, and it still will not follow the curve.
//!
//! # Every number on screen is measured in this process
//!
//! The arithmetic below is transcribed from
//! `crates/isomesh/benches/experiment_p60.rs`, constant for constant: the pole
//! `2^-2`, the gain `1 + 2^-2`, the paper's `c_0 = f_0`, the strict `< 0`
//! inside test, the bisection floor, the two-piece segment choice. The eight
//! reference fields are measured once at startup at slide zero — the committed
//! configuration — and every quoted number is compared against
//! `docs/experiments/p-60.csv` as a `{:.6e}` **string**, so "agrees" means the
//! same digits and not a tolerance somebody chose. The comparison is printed in
//! full on startup and its verdict is on the HUD; a disagreement puts a
//! `CROSS-CHECK FAILED` line up rather than panicking in a stranger's terminal.

mod common;

use bevy::prelude::*;
use common::{Capture, CommonPlugin, DemoStats, OrbitCamera, ViewFlags};
use isomesh::fields::ReferenceField;
use isomesh::{Sdf, for_each_reference_field};

// ─── the arithmetic, transcribed from benches/experiment_p60.rs ─────────────

/// **A sample of exactly zero is outside.**
///
/// Transcribed from `isomesh::cube::is_inside` by way of the P-60 harness,
/// which transcribes it for the same reason: `cube` is a private module. Strict
/// `< 0`, verbatim. It matters more here than usual — on this line several
/// reference fields sample **exactly zero**, and which two samples the bracket
/// names decides what the whole row is about.
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
/// lives: [`POLE`] and [`GAIN`] are exact binary.
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

/// Below this many cells the standard rule is at its own numerical floor and a
/// ratio against it is arithmetic rather than measurement.
const EXACT_FLOOR_CELLS: f64 = 1e-12;

/// The four fields C1 is a statement about.
const SMOOTH: [&str; 4] = ["sphere", "torus", "gyroid", "fbm_terrain"];

/// Half the positive control's second derivative in cell units.
///
/// The harness's value, unchanged, because changing it would move the control
/// off the committed `1.984127e-1` and this demo's opening frame is that row.
const CONTROL_CURVATURE: f64 = 0.02;

/// The sample index the control's crossing bracket starts at.
///
/// The harness pins its root at `31.5`; here the root is `CONTROL_BASE + sigma`,
/// so `sigma = 1/2` reproduces it exactly and every other `sigma` in `(0, 1)`
/// slides the same quadratic along the same fixed grid.
const CONTROL_BASE: f64 = 31.0;

/// The registered bar C1 was measured against: 30% lower is `ratio <= 0.70`.
const REGISTERED_BAR: f64 = 0.70;

/// The `c_n` of the causal prefilter over the samples `start..=end`.
///
/// Index `i` of the result is sample `start + i`. `c_start = f_start` is the
/// paper's `c_0 = f_0`, applied at whichever sample the window begins. At
/// `tau = 1/5` this recursion is `c_n = −2^-2·c_{n−1} + (1 + 2^-2)·f_n` and both
/// constants are exact binary, so it is multiplication-free.
fn prefilter(f: &[f64], start: usize, end: usize) -> Vec<f64> {
    assert!(
        start < end && end < f.len(),
        "prefilter window out of range"
    );
    let mut c = Vec::with_capacity(end - start + 1);
    c.push(f[start]);
    for value in &f[(start + 1)..=end] {
        let prev = c[c.len() - 1];
        c.push(GAIN * value - POLE * prev);
    }
    c
}

/// The shifted reconstruction at `u = x/T`, from coefficients based at `start`.
///
/// The knots sit at `(n + tau)·T`, so on `u ∈ [m − 1 + tau, m + tau]` the only
/// two hats that are non-zero are `n = m − 1` and `n = m` and the value is the
/// straight line between their coefficients.
fn reconstruct(c: &[f64], start: usize, u: f64) -> Option<f64> {
    let m = (u - TAU).floor() + 1.0;
    if m < 1.0 {
        return None;
    }
    let hi = m as usize;
    if hi <= start || hi - start >= c.len() {
        return None;
    }
    let s = u - (m - 1.0 + TAU);
    Some(c[hi - 1 - start] * (1.0 - s) + c[hi - start] * s)
}

/// Which of the two shifted pieces inside the sample bracket held the root.
type Segment = &'static str;

/// The zero of the shifted reconstruction inside the sample bracket
/// `[n − 1, n]`, in cells from the line start, and which piece held it.
///
/// Feeding raw samples as `c` is the negative control, so this takes `c` and
/// never looks at `f`: the endpoint values are computed from the reconstruction
/// itself, which is what makes the sabotage traverse the same code.
///
/// `None` where the reconstruction does not change sign across a bracket the
/// field does — impossible for the prefiltered coefficients, by the
/// interpolation identity, and entirely possible for the sabotage.
fn shifted_root(c: &[f64], start: usize, n: usize) -> Option<(f64, Segment)> {
    if start + 2 > n || n < start || n - start >= c.len() {
        return None;
    }
    let c2 = c[n - 2 - start];
    let c1 = c[n - 1 - start];
    let c0 = c[n - start];
    let v_left = TAU * c2 + (1.0 - TAU) * c1;
    let v_right = TAU * c1 + (1.0 - TAU) * c0;
    if is_inside(v_left) == is_inside(v_right) {
        return None;
    }
    // Exactly one of the two pieces changes sign, by parity, so this is total.
    if is_inside(v_left) == is_inside(c1) {
        Some(((n - 1) as f64 + TAU + c1 / (c1 - c0), "right"))
    } else {
        Some(((n - 2) as f64 + TAU + c2 / (c2 - c1), "left"))
    }
}

/// The closed-form `error_ratio` for a locally quadratic signal crossing at
/// fractional position `sigma`.
///
/// `|q(sigma) − tau(1 − tau)| / (sigma(1 − sigma))`, which at `tau = 1/5` is
/// `|sigma − 2/5| / sigma` above `1/5` and `(sigma + 3/5)/(1 − sigma)` below it.
fn predicted_ratio(sigma: f64) -> f64 {
    let q = if sigma >= TAU {
        (sigma - TAU) * (1.0 - sigma + TAU)
    } else {
        (sigma + 1.0 - TAU) * (TAU - sigma)
    };
    (q - TAU_PRODUCT).abs() / (sigma * (1.0 - sigma))
}

/// The first sample index `i` whose bracket `[i − 1, i]` changes side.
fn first_crossing(f: &[f64]) -> Option<usize> {
    (1..f.len()).find(|&i| is_inside(f[i - 1]) != is_inside(f[i]))
}

/// Bisect a 1-D restriction between two abscissae of opposite side.
///
/// Terminates at [`BISECT_WIDTH`] or at the last representable interval,
/// whichever comes first. Returns the root and the width it stopped at, so the
/// floor is a reported number rather than an assumption.
fn bisect(sample: &dyn Fn(f64) -> f64, mut lo: f64, mut hi: f64) -> (f64, f64) {
    let lo_inside = is_inside(sample(lo));
    loop {
        let mid = f64::midpoint(lo, hi);
        if mid <= lo || mid >= hi || hi - lo <= BISECT_WIDTH {
            break;
        }
        if is_inside(sample(mid)) == lo_inside {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    (f64::midpoint(lo, hi), hi - lo)
}

/// Median by `total_cmp`, mean of the two middle values on an even population.
fn median(values: &[f64]) -> f64 {
    if values.is_empty() {
        return f64::NAN;
    }
    let mut v = values.to_vec();
    v.sort_unstable_by(f64::total_cmp);
    let m = v.len() / 2;
    if v.len().is_multiple_of(2) {
        f64::midpoint(v[m - 1], v[m])
    } else {
        v[m]
    }
}

/// `{:.6e}`, the CSV's own format. Comparisons are made on these strings, so
/// "agrees" means the same digits rather than a tolerance somebody chose.
fn sci(v: f64) -> String {
    format!("{v:.6e}")
}

// ─── one measured line ──────────────────────────────────────────────────────

/// One guard-band window's answer.
#[derive(Clone, Copy, Default)]
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

/// The negative control's answer, where it has one.
///
/// `None` above `sigma = 1 − tau` is not a missing measurement, it is the
/// result: feeding the raw samples to the shifted reconstruction gives the
/// bracket endpoints `tau·f_{n−2} + (1 − tau)·f_{n−1}` and
/// `tau·f_{n−1} + (1 − tau)·f_n`, and on a line of unit slope the second is
/// `(1 − tau) − sigma`. Past `sigma = 0.8` **both endpoints land on the same
/// side** and the un-prefiltered shift does not bracket the crossing at all.
/// Without the prefilter this is not a reconstruction that is merely biased;
/// it is not a reconstruction of the samples.
#[derive(Clone, Copy)]
struct Sabotage {
    /// Its root, in cells from the line start.
    u: f64,
    /// `|u − u_exact|`, in cells. About `tau` on the control.
    error: f64,
    /// `error / error_standard`.
    ratio: f64,
}

/// Everything one 65-sample line measured.
#[derive(Clone)]
struct Line {
    /// The samples, in field value.
    f: Vec<f64>,
    /// The prefilter coefficients over the whole line.
    c: Vec<f64>,
    /// The bracket holding the first crossing is `[n − 1, n]`.
    n: usize,
    /// The exact root, in cells from the line start.
    u_exact: f64,
    /// `t = a/(a − b)`, in cells from the line start.
    u_standard: f64,
    /// The shifted reconstruction's root, in cells from the line start.
    u_shifted: f64,
    /// The negative control, where the un-prefiltered shift finds a root at
    /// all. See [`Sabotage`].
    sabotage: Option<Sabotage>,
    /// `sigma` as the standard rule sees it — `a/(a − b)`, the CSV's
    /// `crossing_fraction`. Not the same number as [`Line::sigma_true`].
    fraction: f64,
    /// `|u_standard − u_exact|`, in cells.
    error_standard: f64,
    /// `|u_shifted − u_exact|`, in cells.
    error_shifted: f64,
    /// `error_shifted / error_standard`, the CSV's `error_ratio`.
    ratio: f64,
    /// Worst `|tau·c_{n−1} + (1 − tau)·c_n − f_n|` over the line.
    interpolation_residual: f64,
    /// Which shifted piece held the root.
    segment: Segment,
    /// The width the bisection stopped at, in field units.
    bisect_floor: f64,
    /// The four guard-band windows.
    guards: [Guard; 4],
}

impl Line {
    /// Measure a line, given its samples and the exact root in cells.
    ///
    /// One measurement path for both signals: the reference-field lines hand it
    /// a bisected root and the quadratic control hands it an analytic one, and
    /// nothing below this point knows which.
    fn measure(f: Vec<f64>, u_exact: f64, bisect_floor: f64) -> Option<Self> {
        let n = first_crossing(&f)?;
        if n < 2 {
            return None;
        }
        let a = f[n - 1];
        let b = f[n];
        let fraction = a / (a - b);
        let u_standard = (n - 1) as f64 + fraction;

        let c = prefilter(&f, 0, SAMPLES - 1);
        let (u_shifted, segment) = shifted_root(&c, 0, n)?;
        // The negative control: the shift without the prefilter that makes it
        // interpolate, through the identical function. `None` past
        // `sigma = 1 − tau`, which is a result rather than a gap — see
        // [`Sabotage`].
        let sabotage = shifted_root(&f, 0, n).map(|(u, _)| u);

        let error_standard = (u_standard - u_exact).abs();
        let error_shifted = (u_shifted - u_exact).abs();
        let sabotage = sabotage.map(|u| {
            let error = (u - u_exact).abs();
            Sabotage {
                u,
                error,
                ratio: error / error_standard,
            }
        });

        let mut interpolation_residual = 0.0f64;
        for i in 1..SAMPLES {
            let residual = (TAU * c[i - 1] + (1.0 - TAU) * c[i] - f[i]).abs();
            interpolation_residual = interpolation_residual.max(residual);
        }

        let mut guards = [Guard::default(); 4];
        for (slot, &k) in guards.iter_mut().zip(GUARD_K.iter()) {
            let start = n.saturating_sub(k);
            let window = prefilter(&f, start, n);
            let (u_truncated, _) = shifted_root(&window, start, n)?;
            let delta = (u_truncated - u_shifted).abs();
            *slot = Guard {
                k,
                start,
                truncated: start > 0,
                delta,
                converged: delta < GUARD_TOLERANCE,
            };
        }

        Some(Self {
            f,
            c,
            n,
            u_exact,
            u_standard,
            u_shifted,
            sabotage,
            fraction,
            error_standard,
            error_shifted,
            ratio: error_shifted / error_standard,
            interpolation_residual,
            segment,
            bisect_floor,
            guards,
        })
    }

    /// Where the root really sits inside its bracket, as a fraction of the cell.
    ///
    /// **Not** `fraction`. `a/(a − b)` is the chord's answer and carries its own
    /// `O(g'')` error — on the control at a true `sigma = 1/2` it reads `0.495`.
    /// The closed form is a statement about where the root *is*, so this is the
    /// abscissa the prediction is evaluated at and the sweep is plotted against.
    fn sigma_true(&self) -> f64 {
        self.u_exact - (self.n - 1) as f64
    }

    /// The closed form at [`Line::sigma_true`], or `None` where there is nothing
    /// to predict: a root sitting exactly on a sample, or a standard rule
    /// already at its numerical floor.
    fn predicted(&self) -> Option<f64> {
        let sigma = self.sigma_true();
        if !(0.0..=1.0).contains(&sigma)
            || sigma <= 0.0
            || sigma >= 1.0
            || self.error_standard < EXACT_FLOOR_CELLS
        {
            return None;
        }
        Some(predicted_ratio(sigma))
    }

    /// Whether the standard rule is already at its own numerical floor, so the
    /// ratio on this line is arithmetic rather than measurement.
    fn standard_is_exact(&self) -> bool {
        self.error_standard < EXACT_FLOOR_CELLS
    }
}

// ─── the signals a line can be cut through ──────────────────────────────────

/// One reference field, restricted to its domain-centre `x` axis.
struct FieldLine {
    /// The field's own name, and the CSV's key.
    name: &'static str,
    /// Line start, in field units.
    x0: f64,
    /// `T`, in field units.
    cell: f64,
    /// `field.sample([x, 0, 0])`, boxed so the eight distinct types can sit in
    /// one list. A virtual call per sample is 65 of them per frame.
    sample: Box<dyn Fn(f64) -> f64 + Send + Sync>,
}

impl FieldLine {
    /// Measure this field's line with the **field** slid `slide` cells along it.
    ///
    /// The grid never moves: sample `i` is always at `x0 + T·i` in the frame the
    /// mesher works in, and sliding the field means asking it for
    /// `f(x − slide·T)`. At `slide = 0` this is the committed configuration, bit
    /// for bit.
    fn line(&self, slide: f64) -> Option<Line> {
        let x_at = |i: f64| self.x0 + self.cell * (i - slide);
        let f: Vec<f64> = (0..SAMPLES).map(|i| (self.sample)(x_at(i as f64))).collect();
        let n = first_crossing(&f)?;
        if n < 2 {
            return None;
        }
        let (x_exact, floor) = bisect(&*self.sample, x_at((n - 1) as f64), x_at(n as f64));
        let u_exact = (x_exact - self.x0) / self.cell + slide;
        Line::measure(f, u_exact, floor)
    }
}

/// The exactly-quadratic control line, with its root at `CONTROL_BASE + sigma`.
///
/// `g(u) = d + C·d²` with `d = u − root`. `C` is small enough that the
/// quadratic's *other* root, at `root − 1/C = root − 50`, is off the line, so
/// the first sign change is the intended one.
fn control_line(sigma: f64) -> Option<Line> {
    let root = CONTROL_BASE + sigma;
    let f: Vec<f64> = (0..SAMPLES)
        .map(|i| {
            let d = i as f64 - root;
            d + CONTROL_CURVATURE * d * d
        })
        .collect();
    Line::measure(f, root, 0.0)
}

/// Every signal the line can be cut through: the control, then the eight
/// reference fields in the order `p-60.csv` lists them.
#[derive(Resource)]
struct Signals {
    /// The eight reference fields, on their domain-centre axis.
    fields: Vec<FieldLine>,
}

impl Signals {
    /// Build the eight restrictions once.
    fn new() -> Self {
        let mut fields: Vec<FieldLine> = Vec::new();
        for_each_reference_field!(f64, |name, field| {
            let (lo, hi) = field.domain();
            let x0 = lo[0];
            let cell = (hi[0] - x0) / CELLS as f64;
            fields.push(FieldLine {
                name,
                x0,
                cell,
                sample: Box::new(move |x: f64| field.sample([x, 0.0, 0.0])),
            });
        });
        Self { fields }
    }

    /// How many signals there are, control included.
    fn len(&self) -> usize {
        self.fields.len() + 1
    }

    /// The name of signal `index`.
    fn name(&self, index: usize) -> &'static str {
        if index == 0 {
            "quadratic control"
        } else {
            self.fields[index - 1].name
        }
    }

    /// `T` for signal `index`, in field units. The control is in cell units, so
    /// its own `T` is exactly one.
    fn cell(&self, index: usize) -> f64 {
        if index == 0 {
            1.0
        } else {
            self.fields[index - 1].cell
        }
    }

    /// Measure signal `index` with the field slid `slide` cells.
    fn line(&self, index: usize, slide: f64) -> Option<Line> {
        if index == 0 {
            control_line(slide)
        } else {
            self.fields[index - 1].line(slide)
        }
    }

    /// The slide range the arrows and the story are allowed to drive.
    ///
    /// The control's slide **is** `sigma`, so it is clamped inside `(0, 1)`:
    /// at either end the root lands exactly on a sample, both errors collapse
    /// to zero and the ratio is `0/0`.
    fn slide_range(index: usize) -> (f64, f64) {
        if index == 0 {
            (0.02, 0.98)
        } else {
            (-1.0, 1.0)
        }
    }
}

// ─── the ledger, compiled in ────────────────────────────────────────────────

/// P-60's committed artefact, embedded at compile time.
///
/// `include_str!` rather than transcribed constants: the path resolves against
/// this source file so no working directory can break it, and a number that
/// lived only here could drift away from the CSV with nothing to say so.
const LEDGER_CSV: &str = include_str!("../../docs/experiments/p-60.csv");

/// The committed CSV, indexed by header name.
#[derive(Resource)]
struct Ledger {
    /// Column names, in file order.
    header: Vec<&'static str>,
    /// Data rows, split on commas, all the same width as the header.
    rows: Vec<Vec<&'static str>>,
}

impl Ledger {
    /// Parse the embedded CSV, skipping its `#` provenance header.
    fn load() -> Option<Self> {
        let mut lines = LEDGER_CSV.lines().filter(|l| !l.starts_with('#'));
        let header: Vec<&str> = lines.next()?.split(',').collect();
        let rows: Vec<Vec<&str>> = lines
            .map(|l| l.split(',').collect::<Vec<&str>>())
            .filter(|r| r.len() == header.len())
            .collect();
        if rows.is_empty() {
            return None;
        }
        Some(Self { header, rows })
    }

    /// One cell, by field name, guard-band width and column name.
    fn cell(&self, field: &str, k: usize, column: &str) -> Option<&'static str> {
        let c_field = self.header.iter().position(|h| *h == "field")?;
        let c_k = self.header.iter().position(|h| *h == "guard_band_k")?;
        let c = self.header.iter().position(|h| *h == column)?;
        let key = k.to_string();
        self.rows
            .iter()
            .find(|r| r[c_field] == field && r[c_k] == key)
            .map(|r| r[c])
    }

    /// How many rows the artefact carries.
    fn rows(&self) -> usize {
        self.rows.len()
    }
}

/// One quoted number, and whether this run reproduced it.
struct Check {
    /// What the number is.
    label: &'static str,
    /// The digits `p-60.csv` carries.
    expected: String,
    /// The digits this run produced, in the same format.
    measured: String,
}

impl Check {
    /// Whether the two strings are identical.
    fn agrees(&self) -> bool {
        self.expected == self.measured
    }
}

// ─── what the eight fields and the two controls measured, once ──────────────

/// The startup measurement: the committed configuration, re-run here.
#[derive(Resource)]
struct Baseline {
    /// One line per reference field, at slide zero.
    fields: Vec<(&'static str, Line)>,
    /// The positive control, at `sigma = 1/2`.
    control: Line,
    /// C1's statistic: the median `error_ratio` over the four smooth fields.
    median_smooth: f64,
    /// How many fields sit at the bisection floor.
    exact_count: usize,
    /// The quoted numbers, checked against the CSV.
    checks: Vec<Check>,
    /// The closed form, checked across the whole `sigma` sweep.
    sweep: SweepCheck,
    /// Rows in the artefact.
    ledger_rows: usize,
}

impl Baseline {
    /// How many quoted numbers agree.
    fn agreed(&self) -> usize {
        self.checks.iter().filter(|c| c.agrees()).count()
    }

    /// The guard row for one field, by name.
    fn guards(&self, name: &str) -> Option<&[Guard; 4]> {
        self.fields
            .iter()
            .find(|(n, _)| *n == name)
            .map(|(_, l)| &l.guards)
    }
}

/// How well the measurement follows the closed form across the whole sweep.
struct SweepCheck {
    /// How many `sigma` positions were tested.
    positions: usize,
    /// Median `|measured − predicted| / predicted`.
    median_rel: f64,
    /// Worst of the same, and where.
    worst_rel: f64,
    /// The `sigma` the worst deviation sits at.
    worst_sigma: f64,
    /// The measured ratio at `sigma = tau`, where the closed form says `1`.
    at_tau: f64,
    /// The measured ratio at `sigma = 2·tau`, where the closed form says `0`.
    at_two_tau: f64,
    /// The largest measured ratio anywhere in the sweep.
    max_measured: f64,
    /// The fraction of tested `sigma` that clears [`REGISTERED_BAR`].
    clear_fraction: f64,
}

/// The closed form's own answer to "how often does this clear the 30% bar".
///
/// `ratio <= 0.7` holds on `sigma <= 1/17` below `tau` and on `sigma >= 4/17`
/// above it, so the measure is `1 − 4/17 + 1/17 = 14/17`.
const CLEAR_FRACTION_EXACT: f64 = 14.0 / 17.0;

/// Run the sweep the demo animates, once, and report how closely it tracks.
fn check_sweep() -> SweepCheck {
    let positions = 199usize;
    let mut rel: Vec<f64> = Vec::with_capacity(positions);
    let mut worst_rel = 0.0f64;
    let mut worst_sigma = f64::NAN;
    let mut max_measured = 0.0f64;
    let mut cleared = 0usize;
    for i in 1..=positions {
        let sigma = i as f64 / (positions + 1) as f64;
        let Some(line) = control_line(sigma) else {
            continue;
        };
        let predicted = predicted_ratio(sigma);
        max_measured = max_measured.max(line.ratio);
        if line.ratio <= REGISTERED_BAR {
            cleared += 1;
        }
        if predicted > 1e-3 {
            let deviation = (line.ratio - predicted).abs() / predicted;
            rel.push(deviation);
            if deviation > worst_rel {
                worst_rel = deviation;
                worst_sigma = sigma;
            }
        }
    }
    let at = |sigma: f64| control_line(sigma).map_or(f64::NAN, |l| l.ratio);
    SweepCheck {
        positions,
        median_rel: median(&rel),
        worst_rel,
        worst_sigma,
        at_tau: at(TAU),
        at_two_tau: at(2.0 * TAU),
        max_measured,
        clear_fraction: cleared as f64 / positions as f64,
    }
}

/// Measure the eight fields, run the controls, and check every quoted number.
fn measure_baseline(signals: &Signals, ledger: &Ledger) -> Baseline {
    let fields: Vec<(&'static str, Line)> = signals
        .fields
        .iter()
        .filter_map(|f| f.line(0.0).map(|l| (f.name, l)))
        .collect();
    let control = control_line(0.5).expect("the quadratic control crosses zero by construction");

    let smooth: Vec<f64> = fields
        .iter()
        .filter(|(n, _)| SMOOTH.contains(n))
        .map(|(_, l)| l.ratio)
        .collect();
    let median_smooth = median(&smooth);
    let exact_count = fields
        .iter()
        .filter(|(_, l)| l.standard_is_exact())
        .count();

    let expect = |field: &str, k: usize, column: &'static str| {
        ledger
            .cell(field, k, column)
            .unwrap_or("MISSING")
            .to_string()
    };
    let ratio_of = |name: &str| {
        fields
            .iter()
            .find(|(n, _)| *n == name)
            .map_or(f64::NAN, |(_, l)| l.ratio)
    };
    let guard_of = |name: &str, k: usize| {
        fields
            .iter()
            .find(|(n, _)| *n == name)
            .and_then(|(_, l)| l.guards.iter().find(|g| g.k == k))
            .map_or(f64::NAN, |g| g.delta)
    };
    let shifted_of = |name: &str| {
        fields
            .iter()
            .find(|(n, _)| *n == name)
            .map_or(f64::NAN, |(_, l)| l.error_shifted)
    };
    let standard_of = |name: &str| {
        fields
            .iter()
            .find(|(n, _)| *n == name)
            .map_or(f64::NAN, |(_, l)| l.error_standard)
    };

    // The CSV's own count of degenerate rows, so both sides of the last check
    // come out of a measurement rather than out of this file.
    let ledger_exact = signals
        .fields
        .iter()
        .filter(|f| expect(f.name, 2, "standard_is_exact") == "true")
        .count();

    let checks = vec![
        Check {
            label: "median_error_ratio, the four smooth fields (C1)",
            expected: expect("sphere", 2, "median_error_ratio"),
            measured: sci(median_smooth),
        },
        Check {
            label: "fbm_terrain  error_ratio        4.9x better",
            expected: expect("fbm_terrain", 2, "error_ratio"),
            measured: sci(ratio_of("fbm_terrain")),
        },
        Check {
            label: "gyroid       error_ratio        2x worse",
            expected: expect("gyroid", 2, "error_ratio"),
            measured: sci(ratio_of("gyroid")),
        },
        Check {
            label: "torus  guard_band_delta_cells   k = 2",
            expected: expect("torus", 2, "guard_band_delta_cells"),
            measured: sci(guard_of("torus", 2)),
        },
        Check {
            label: "torus  guard_band_delta_cells   k = 5   (fails the bar)",
            expected: expect("torus", 5, "guard_band_delta_cells"),
            measured: sci(guard_of("torus", 5)),
        },
        Check {
            label: "torus  guard_band_delta_cells   k = 10  (passes it)",
            expected: expect("torus", 10, "guard_band_delta_cells"),
            measured: sci(guard_of("torus", 10)),
        },
        Check {
            label: "gyroid guard_band_delta_cells   k = 10  (passes it)",
            expected: expect("gyroid", 10, "guard_band_delta_cells"),
            measured: sci(guard_of("gyroid", 10)),
        },
        Check {
            label: "control_quadratic_error_ratio   positive control",
            expected: expect("sphere", 2, "control_quadratic_error_ratio"),
            measured: sci(control.ratio),
        },
        Check {
            label: "control_sabotage_error_ratio    negative control",
            expected: expect("sphere", 2, "control_sabotage_error_ratio"),
            measured: sci(control.sabotage.map_or(f64::NAN, |s| s.ratio)),
        },
        Check {
            label: "sphere root_error_standard      the bisection floor",
            expected: expect("sphere", 2, "root_error_standard"),
            measured: sci(standard_of("sphere")),
        },
        Check {
            label: "torus  root_error_shifted       the causal transient",
            expected: expect("torus", 2, "root_error_shifted"),
            measured: sci(shifted_of("torus")),
        },
        Check {
            label: "fields at the bisection floor   standard_is_exact",
            expected: ledger_exact.to_string(),
            measured: exact_count.to_string(),
        },
    ];

    Baseline {
        fields,
        control,
        median_smooth,
        exact_count,
        checks,
        sweep: check_sweep(),
        ledger_rows: ledger.rows(),
    }
}

/// Print the whole cross-check, once, before the window opens.
fn report_baseline(base: &Baseline, signals: &Signals) {
    info!(
        "E-316 shifted_linear_root -- M-359 / x42, the demo for P-60 (R-058); \
         artefact docs/experiments/p-60.csv, {} rows",
        base.ledger_rows
    );
    info!("cross-check, compared as {{:.6e}} strings -- same digits or not at all:");
    for check in &base.checks {
        info!(
            "  {:<46}  expected (p-60.csv) = {:>13}   measured = {:>13}   {}",
            check.label,
            check.expected,
            check.measured,
            if check.agrees() { "agree" } else { "DISAGREE" }
        );
    }
    info!(
        "  {} of {} quoted numbers agree",
        base.agreed(),
        base.checks.len()
    );

    info!("the eight reference lines, at slide zero -- the committed configuration:");
    for (name, line) in &base.fields {
        info!(
            "  {:<15} T {:.7}  bracket [{:>2}, {:>2}]  sigma {:.6}  err_std {:.6e}  \
             err_shift {:.6e}  ratio {:.6e}  bisected to {:.3e}  {}",
            name,
            signals
                .fields
                .iter()
                .find(|f| f.name == *name)
                .map_or(f64::NAN, |f| f.cell),
            line.n - 1,
            line.n,
            line.fraction,
            line.error_standard,
            line.error_shifted,
            line.ratio,
            line.bisect_floor,
            if line.standard_is_exact() {
                "standard rule at the BISECTION FLOOR -- symmetry axis, piecewise linear"
            } else {
                "a real measurement"
            }
        );
    }
    info!(
        "  {} of {} fields are degenerate on the domain-centre axis, so C1's population \
         was half arithmetic",
        base.exact_count,
        base.fields.len()
    );

    info!(
        "the guard band -- delta from the whole-line root, in cells, bar at 1e-6; \
         `from` is the first sample the recursion started at:"
    );
    for (name, line) in &base.fields {
        let cells: Vec<String> = line
            .guards
            .iter()
            .map(|g| {
                format!(
                    "k={:<2} from {:>2} {:.6e}{}",
                    g.k,
                    g.start,
                    g.delta,
                    if g.truncated { "" } else { " (whole line)" }
                )
            })
            .collect();
        info!("  {:<15} {}", name, cells.join("  "));
    }

    let s = &base.sweep;
    info!(
        "the closed form |sigma - 2 tau| / sigma, swept over {} positions of the exact quadratic:",
        s.positions
    );
    info!(
        "  measured vs predicted: median {:.3}% deviation, worst {:.3}% at sigma {:.4}",
        s.median_rel * 100.0,
        s.worst_rel * 100.0,
        s.worst_sigma
    );
    info!(
        "  at sigma = tau  measured {:.6} against a predicted 1;  at sigma = 2 tau  measured \
         {:.6e} against a predicted 0",
        s.at_tau, s.at_two_tau
    );
    info!(
        "  largest measured ratio anywhere {:.6};  clears the registered {:.2} bar on {:.2}% of \
         sigma against the closed form's 14/17 = {:.2}%",
        s.max_measured,
        REGISTERED_BAR,
        s.clear_fraction * 100.0,
        CLEAR_FRACTION_EXACT * 100.0
    );
    info!(
        "positive control  ratio {:.6e} against a closed-form 2.000000e-1;  negative control \
         (raw samples, no prefilter) ratio {:.6e}, root biased {:.6} cells = tau",
        base.control.ratio,
        base.control.sabotage.map_or(f64::NAN, |s| s.ratio),
        base.control.sabotage.map_or(f64::NAN, |s| s.error)
    );
    info!(
        "  and above sigma = 1 - tau = 0.8 the un-prefiltered shift has no root in the bracket \
         at all: both endpoints land on the same side, so it is not a reconstruction of the \
         samples rather than a biased one"
    );
    info!(
        "NOT SHIPPABLE as an accuracy win: a mesher cannot choose sigma, so the expected benefit \
         is the 82% figure and not the 8 dB one. edge_position is unchanged. The guard band is \
         what ships: 10 samples."
    );
}

// ─── the story ──────────────────────────────────────────────────────────────

/// Where the opening hold on the committed control ends, as a clip fraction.
const HOLD_END: f32 = 0.14;

/// Where the sweep ends and the closing rest begins.
const SWEEP_END: f32 = 0.88;

/// How far into the sweep the descending leg turns around.
///
/// The clip opens on the control at `sigma = 1/2`, walks **down** past the zero
/// at `2 tau` and the peak at `tau`, then turns and climbs the whole way to
/// `0.97`. One continuous motion with no cut, and it crosses both landmarks
/// twice.
const SWEEP_TURN: f32 = 0.34;

/// Seconds for one pass through the story, when nobody is capturing.
const STORY_SECONDS: f32 = 24.0;

/// How fast the arrow keys drive `sigma`, in cells per second.
const MANUAL_RATE: f64 = 0.30;

/// One beat of the story.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum Beat {
    /// Held at `sigma = 1/2`: the committed positive control.
    Control,
    /// `sigma` sliding, the trace filling in the closed form.
    Sweep,
    /// Resting at `sigma = 2·tau`, where the shifted root is exact.
    Zero,
    /// A hand on the arrow keys.
    Manual,
    /// A reference field, pinned at the committed configuration.
    Committed,
}

/// The story's `sigma` at clip fraction `phase`, and which beat that is.
fn story(phase: f32) -> (Beat, f64) {
    let lerp = |a: f32, b: f32, t: f32| f64::from(a + (b - a) * t.clamp(0.0, 1.0));
    if phase < HOLD_END {
        (Beat::Control, 0.5)
    } else if phase < SWEEP_END {
        let local = (phase - HOLD_END) / (SWEEP_END - HOLD_END);
        if local < SWEEP_TURN {
            (Beat::Sweep, lerp(0.5, 0.03, local / SWEEP_TURN))
        } else {
            (
                Beat::Sweep,
                lerp(0.03, 0.97, (local - SWEEP_TURN) / (1.0 - SWEEP_TURN)),
            )
        }
    } else {
        let local = ((phase - SWEEP_END) / (1.0 - SWEEP_END)).clamp(0.0, 1.0);
        // Smoothstep, so the settle onto the zero has no visible stop.
        let t = local * local * (3.0 - 2.0 * local);
        (Beat::Zero, lerp(0.97, 0.40, t))
    }
}

/// `ISOMESH_CAPTURE_FRAMES`, or the harness default.
fn capture_frames() -> u32 {
    std::env::var("ISOMESH_CAPTURE_FRAMES")
        .ok()
        .and_then(|v| v.parse().ok())
        .filter(|n| *n > 0)
        .unwrap_or(60)
}

// ─── layout, in world units ─────────────────────────────────────────────────

/// The wide view of the line: `x`, `y`, width, height.
const PANEL_LINE: [f32; 4] = [0.0, 6.20, 5.60, 3.60];
/// The magnified inset, where the three roots separate.
const PANEL_INSET: [f32; 4] = [6.20, 6.20, 4.00, 3.60];
/// The ratio against `sigma`.
const PANEL_RATIO: [f32; 4] = [0.0, 1.60, 4.60, 4.00];
/// The guard band.
const PANEL_GUARD: [f32; 4] = [5.60, 1.60, 4.60, 4.00];

/// How many cells either side of the crossing the wide view shows.
const LINE_HALF_CELLS: f64 = 3.5;

/// The inset half-width, as a multiple of the larger of the two root errors.
const INSET_MARGIN: f64 = 1.8;

/// The smallest inset half-width, in cells. Below this the two reconstructions
/// are the same `f64` and the picture is honestly empty.
const INSET_FLOOR: f64 = 1e-14;

/// Log-10 bottom and top of the guard-band axis, in cells.
const GUARD_LOG: [f64; 2] = [-8.0, -1.0];

/// Camera distance, chosen so the whole chart clears the HUD panel.
const CAMERA_RADIUS: f32 = 12.8;

/// Where the camera looks, so the chart sits in the right 56% of the frame.
const CAMERA_FOCUS: Vec3 = Vec3::new(1.23, 5.70, 0.0);

/// Width of the HUD backdrop, in logical pixels.
const HUD_PANEL: Vec2 = Vec2::new(548.0, 660.0);

/// How many chart labels are pre-spawned. Unused ones are `Display::None`.
const MAX_LABELS: usize = 44;

// ─── colour ─────────────────────────────────────────────────────────────────

/// The samples and everything the standard rule produces.
const STANDARD: Color = Color::srgb(0.38, 0.72, 1.00);
/// The prefilter coefficients and everything the shifted rule produces.
const SHIFTED: Color = Color::srgb(1.00, 0.62, 0.18);
/// The exact root, from bisecting the field.
const EXACT: Color = Color::srgb(0.97, 0.98, 1.00);
/// Axes, frames and ticks.
const AXIS: Color = Color::srgb(0.42, 0.46, 0.56);
/// A threshold that decides a clause.
const BAR: Color = Color::srgb(0.96, 0.34, 0.34);
/// The closed form.
const CLOSED_FORM: Color = Color::srgb(0.42, 0.94, 0.55);
/// The live measurement, and its trail.
const MEASURED: Color = Color::srgb(1.00, 0.88, 0.24);
/// The negative control.
const SABOTAGE: Color = Color::srgb(0.92, 0.30, 0.62);

/// Line gizmos for the curves.
#[derive(Default, Reflect, GizmoConfigGroup)]
struct ChartGizmos;

/// Thinner gizmos for axes, frames and ticks.
#[derive(Default, Reflect, GizmoConfigGroup)]
struct FrameGizmos;

// ─── resources ──────────────────────────────────────────────────────────────

/// What this frame is showing.
#[derive(Resource)]
struct Shot {
    /// Index into [`Signals`], `0` being the control.
    probe: usize,
    /// Cells the field is slid along the line. For the control this **is**
    /// `sigma`.
    slide: f64,
    /// Which beat, including the two the story never reaches on its own.
    beat: Beat,
    /// Wall-clock progress through the story, kept so a loop can be detected.
    phase: f32,
    /// Whether a hand is on the arrow keys.
    manual: bool,
}

impl Default for Shot {
    fn default() -> Self {
        Self {
            probe: 0,
            slide: 0.5,
            beat: Beat::Control,
            phase: 0.0,
            manual: false,
        }
    }
}

/// The live measurement, and the trail the sweep has drawn so far.
#[derive(Resource, Default)]
struct Live {
    /// This frame's line, or `None` where the signal has no usable crossing.
    line: Option<Line>,
    /// `(sigma, ratio)` for every sweep position visited since the last reset.
    trail: Vec<Vec2>,
}

/// One chart label: where it points and what it says.
#[derive(Clone)]
struct LabelSpec {
    /// The world point the label is anchored to.
    at: Vec3,
    /// What it says.
    text: String,
    /// What colour it says it in.
    colour: Color,
    /// Horizontal anchoring against [`LabelSpec::at`].
    justify: JustifyContent,
    /// Vertical anchoring against [`LabelSpec::at`].
    align: AlignItems,
    /// Point size.
    size: f32,
}

/// The labels this frame wants drawn, rebuilt every frame.
#[derive(Resource, Default)]
struct LabelSpecs(Vec<LabelSpec>);

/// One pre-spawned label, by index into [`LabelSpecs`].
#[derive(Component)]
struct ChartLabel(usize);

/// The bottom caption — the line a viewer reads instead of the HUD.
#[derive(Component)]
struct Caption;

// ─── app ────────────────────────────────────────────────────────────────────

fn main() {
    App::new()
        .add_plugins(DefaultPlugins.set(WindowPlugin {
            primary_window: Some(Window {
                title: "isomesh - E-316 shifted linear root".into(),
                ..default()
            }),
            ..default()
        }))
        .add_plugins(CommonPlugin)
        .init_gizmo_group::<ChartGizmos>()
        .init_gizmo_group::<FrameGizmos>()
        // Not black: the panels are lit slabs and a void behind them reads as a
        // rendering failure rather than as a backdrop.
        .insert_resource(ClearColor(Color::srgb(0.055, 0.065, 0.090)))
        .init_resource::<Shot>()
        .init_resource::<Live>()
        .init_resource::<LabelSpecs>()
        .add_systems(Startup, setup)
        // `PreUpdate`, not `Update`: the harness's `update_hud` lives in
        // `Update` and system order within a schedule is unspecified, so a HUD
        // written there would read the previous frame's numbers while the
        // caption under it carried this frame's. Two numbers on screen
        // disagreeing is worse than either being late.
        .add_systems(PreUpdate, (controls, advance, measure, report).chain())
        .add_systems(Update, (draw, place_labels).chain())
        .run();
}

fn setup(
    mut commands: Commands,
    mut gizmo_config: ResMut<GizmoConfigStore>,
    mut materials: ResMut<Assets<StandardMaterial>>,
    mut meshes: ResMut<Assets<Mesh>>,
    mut camera: Query<&mut OrbitCamera>,
) {
    gizmo_config.config_mut::<ChartGizmos>().0.line.width = 2.4;
    gizmo_config.config_mut::<FrameGizmos>().0.line.width = 1.1;

    for mut orbit in &mut camera {
        // Straight on. This is an instrument, and an instrument in perspective
        // is an instrument you cannot read a value off. `yaw = PI/2`, not 0: the
        // harness builds its direction as `(cos yaw cos pitch, sin pitch,
        // sin yaw cos pitch)`, so `yaw = 0` looks along +X and renders an XY
        // chart edge-on. Orbiting still works — the labels are projected through
        // the camera every frame, so they follow the geometry rather than
        // sitting at fixed pixels and drifting off it.
        orbit.focus = CAMERA_FOCUS;
        orbit.yaw = std::f32::consts::FRAC_PI_2;
        orbit.pitch = 0.0;
        orbit.radius = CAMERA_RADIUS;
    }

    // One lit slab behind each panel. The chart itself is gizmos, which are
    // unlit and would otherwise float against the clear colour with nothing to
    // say where a panel begins and ends.
    let slab = materials.add(StandardMaterial {
        base_color: Color::srgb(0.105, 0.125, 0.165),
        perceptual_roughness: 0.85,
        metallic: 0.0,
        ..default()
    });
    for rect in [PANEL_LINE, PANEL_INSET, PANEL_RATIO, PANEL_GUARD] {
        let mesh = meshes.add(Cuboid::new(rect[2] + 0.28, rect[3] + 0.28, 0.10));
        commands.spawn((
            Mesh3d(mesh),
            MeshMaterial3d(slab.clone()),
            Transform::from_xyz(rect[0] + rect[2] * 0.5, rect[1] + rect[3] * 0.5, -0.09),
        ));
    }

    // Behind the harness HUD, which `CommonPlugin` spawns at the default z and
    // this example leaves alone. `GlobalZIndex(-1)` is the whole mechanism.
    commands.spawn((
        Node {
            position_type: PositionType::Absolute,
            top: Val::Px(6.0),
            left: Val::Px(6.0),
            width: Val::Px(HUD_PANEL.x),
            height: Val::Px(HUD_PANEL.y),
            ..default()
        },
        BackgroundColor(Color::srgba(0.03, 0.04, 0.06, 0.62)),
        GlobalZIndex(-1),
    ));

    commands
        .spawn((
            Node {
                position_type: PositionType::Absolute,
                left: Val::Px(0.0),
                right: Val::Px(0.0),
                bottom: Val::Px(16.0),
                justify_content: JustifyContent::Center,
                ..default()
            },
            GlobalZIndex(4),
        ))
        .with_children(|parent| {
            parent.spawn((
                Text::new(""),
                TextFont {
                    font_size: FontSize::Px(19.0),
                    ..default()
                },
                // `NoWrap`: in a centring flex row the measure is handed the
                // container's whole width but the node's height resolves before
                // the wrap, so a soft wrap pushes the second line off frame.
                TextLayout {
                    linebreak: bevy::text::LineBreak::NoWrap,
                    ..default()
                },
                TextColor(Color::srgb(0.97, 0.95, 0.91)),
                BackgroundColor(Color::srgba(0.03, 0.03, 0.05, 0.84)),
                Node {
                    padding: UiRect::axes(Val::Px(14.0), Val::Px(7.0)),
                    ..default()
                },
                Caption,
            ));
        });

    // A zero-size flex container per label, so the text can be anchored left,
    // centred or right on a projected world point without measuring it: the
    // child simply overflows the container in whichever direction
    // `justify_content` sends it.
    for index in 0..MAX_LABELS {
        commands
            .spawn((
                Node {
                    position_type: PositionType::Absolute,
                    display: Display::None,
                    width: Val::Px(0.0),
                    height: Val::Px(0.0),
                    ..default()
                },
                GlobalZIndex(3),
                ChartLabel(index),
            ))
            .with_children(|parent| {
                parent.spawn((
                    Text::new(""),
                    TextFont {
                        font_size: FontSize::Px(11.0),
                        ..default()
                    },
                    TextLayout {
                        linebreak: bevy::text::LineBreak::NoWrap,
                        ..default()
                    },
                    TextColor(AXIS),
                ));
            });
    }

    let signals = Signals::new();
    let ledger = Ledger::load().expect("p-60.csv parses");
    let baseline = measure_baseline(&signals, &ledger);
    report_baseline(&baseline, &signals);
    commands.insert_resource(signals);
    commands.insert_resource(ledger);
    commands.insert_resource(baseline);
}

/// Keys: which signal, and whether a hand is driving `sigma`.
fn controls(
    keys: Res<ButtonInput<KeyCode>>,
    time: Res<Time>,
    capture: Res<Capture>,
    signals: Res<Signals>,
    mut flags: ResMut<ViewFlags>,
    mut shot: ResMut<Shot>,
) {
    // A capture is a scripted sequence; a stray keypress in it would be a frame
    // nobody can reproduce.
    if capture.is_active() {
        return;
    }
    // The harness's digits cover `1`-`7`, i.e. signals 0-6. These are the two it
    // does not reach, written to the same `flags.field` so there is one source
    // of truth for which signal is up.
    for (key, index) in [(KeyCode::Digit8, 7usize), (KeyCode::Digit9, 8usize)] {
        if keys.just_pressed(key) {
            flags.field = index;
        }
    }
    if keys.just_pressed(KeyCode::ArrowUp) {
        flags.field = (flags.field + 1) % signals.len();
    }
    if keys.just_pressed(KeyCode::ArrowDown) {
        flags.field = (flags.field + signals.len() - 1) % signals.len();
    }

    let right = keys.pressed(KeyCode::ArrowRight);
    let left = keys.pressed(KeyCode::ArrowLeft);
    if right || left {
        let direction = f64::from(i8::from(right) - i8::from(left));
        let (lo, hi) = Signals::slide_range(shot.probe);
        shot.slide =
            (shot.slide + direction * MANUAL_RATE * f64::from(time.delta_secs())).clamp(lo, hi);
        shot.manual = true;
    }
    // `R` is the harness's re-mesh key and there is nothing here to re-mesh, so
    // it hands the sweep back to the story instead.
    if flags.remesh_requested {
        shot.manual = false;
    }
}

/// Decide what this frame is about.
///
/// Under capture the story advances with the captured frame count, so a clip of
/// any length is the whole story rather than a slice of it. Interactively it
/// runs on wall-clock time and loops.
fn advance(
    time: Res<Time>,
    capture: Res<Capture>,
    flags: Res<ViewFlags>,
    signals: Res<Signals>,
    mut shot: ResMut<Shot>,
    mut live: ResMut<Live>,
    mut elapsed: Local<f32>,
) {
    let probe = flags.field.min(signals.len() - 1);
    if probe != shot.probe {
        shot.probe = probe;
        shot.manual = false;
        live.trail.clear();
    }

    if shot.manual {
        shot.beat = Beat::Manual;
        return;
    }

    // A reference field is pinned at the committed configuration. Sliding one
    // sweeps `sigma` through a local shape the closed form does not describe --
    // `gyroid` crosses at an inflection and `sphere` restricts to a straight
    // line -- so the story does not pretend otherwise, and the arrows still do
    // it by hand for anyone who wants to watch it fail.
    if probe != 0 {
        shot.beat = Beat::Committed;
        shot.slide = 0.0;
        live.trail.clear();
        return;
    }

    let phase = if capture.is_active() {
        f32::from(u16::try_from(capture.taken).unwrap_or(u16::MAX))
            / f32::from(u16::try_from(capture_frames()).unwrap_or(1).max(1))
    } else {
        if !flags.paused {
            *elapsed += time.delta_secs();
        }
        (*elapsed / STORY_SECONDS).fract()
    }
    .clamp(0.0, 1.0);

    if phase < shot.phase {
        live.trail.clear();
    }
    shot.phase = phase;
    let (beat, sigma) = story(phase);
    shot.beat = beat;
    shot.slide = sigma;
}

/// Measure the line this frame asks for, and extend the trail.
fn measure(shot: Res<Shot>, signals: Res<Signals>, mut live: ResMut<Live>) {
    live.line = signals.line(shot.probe, shot.slide);
    let Some(line) = &live.line else {
        return;
    };
    if shot.probe != 0 {
        return;
    }
    let sigma = line.sigma_true();
    if !(0.0..=1.0).contains(&sigma) || !line.ratio.is_finite() {
        return;
    }
    let point = Vec2::new(sigma as f32, line.ratio as f32);
    // One sample per frame would draw the same point sixty times while the
    // story holds. Only a moved `sigma` earns a mark.
    if live.trail.last().is_none_or(|p| (p.x - point.x).abs() > 1e-4) {
        live.trail.push(point);
    }
}

// ─── drawing ────────────────────────────────────────────────────────────────

/// A point inside a panel, from fractions of its width and height.
fn at(rect: [f32; 4], fx: f32, fy: f32) -> Vec3 {
    Vec3::new(rect[0] + rect[2] * fx, rect[1] + rect[3] * fy, 0.0)
}

/// Draw a panel's frame.
fn frame(gizmos: &mut Gizmos<FrameGizmos>, rect: [f32; 4]) {
    let corners = [
        at(rect, 0.0, 0.0),
        at(rect, 1.0, 0.0),
        at(rect, 1.0, 1.0),
        at(rect, 0.0, 1.0),
    ];
    for i in 0..4 {
        gizmos.line(corners[i], corners[(i + 1) % 4], AXIS.with_alpha(0.55));
    }
}

/// Everything on screen except the HUD and the caption.
fn draw(
    live: Res<Live>,
    shot: Res<Shot>,
    base: Res<Baseline>,
    signals: Res<Signals>,
    mut chart: Gizmos<ChartGizmos>,
    mut frames: Gizmos<FrameGizmos>,
    mut labels: ResMut<LabelSpecs>,
) {
    labels.0.clear();
    for rect in [PANEL_LINE, PANEL_INSET, PANEL_RATIO, PANEL_GUARD] {
        frame(&mut frames, rect);
    }
    draw_ratio(&live, &base, &mut chart, &mut frames, &mut labels);
    draw_guard(&live, &base, &signals, shot.probe, &mut chart, &mut labels);
    let Some(line) = &live.line else {
        labels.0.push(LabelSpec {
            at: at(PANEL_LINE, 0.5, 0.5),
            text: format!(
                "{} has no usable crossing at slide {:+.3}",
                signals.name(shot.probe),
                shot.slide
            ),
            colour: BAR,
            justify: JustifyContent::Center,
            align: AlignItems::Center,
            size: 13.0,
        });
        return;
    };
    draw_line(line, signals.name(shot.probe), &mut chart, &mut frames, &mut labels);
    draw_inset(line, &mut chart, &mut frames, &mut labels);
}

/// The wide view: the samples, both reconstructions, the shifted knots, and the
/// three roots that are too close together to tell apart here.
fn draw_line(
    line: &Line,
    name: &str,
    chart: &mut Gizmos<ChartGizmos>,
    frames: &mut Gizmos<FrameGizmos>,
    labels: &mut LabelSpecs,
) {
    let rect = PANEL_LINE;
    let u_lo = line.u_exact - LINE_HALF_CELLS;
    let u_hi = line.u_exact + LINE_HALF_CELLS;
    let lo_index = u_lo.ceil().max(0.0) as usize;
    let hi_index = (u_hi.floor().min((SAMPLES - 1) as f64) as usize).max(lo_index);

    // Scale to the largest magnitude anywhere in the window, samples and
    // coefficients alike -- the coefficients overshoot on a steep line and a
    // scale that ignored them would clip the orange curve off the panel.
    let mut vmax = 1e-12f64;
    for i in lo_index..=hi_index {
        vmax = vmax.max(line.f[i].abs()).max(line.c[i].abs());
    }

    let x_of = |u: f64| rect[0] + rect[2] * ((u - u_lo) / (u_hi - u_lo)) as f32;
    let y_zero = rect[1] + rect[3] * 0.52;
    let y_of = |v: f64| y_zero + rect[3] * 0.40 * (v / vmax) as f32;
    let p = |u: f64, v: f64| Vec3::new(x_of(u), y_of(v), 0.0);

    frames.line(
        Vec3::new(rect[0], y_zero, 0.0),
        Vec3::new(rect[0] + rect[2], y_zero, 0.0),
        AXIS,
    );

    // The grid, and the shifted knot grid beside it. The whole point of the
    // panel is that the orange ticks sit tau of a cell to the right of the blue
    // ones and the reconstruction is built on those.
    for i in lo_index..=hi_index {
        let x = x_of(i as f64);
        frames.line(
            Vec3::new(x, y_zero - 0.10, 0.0),
            Vec3::new(x, y_zero + 0.10, 0.0),
            STANDARD.with_alpha(0.45),
        );
        let xk = x_of(i as f64 + TAU);
        frames.line(
            Vec3::new(xk, y_zero - 0.07, 0.0),
            Vec3::new(xk, y_zero + 0.07, 0.0),
            SHIFTED.with_alpha(0.45),
        );
    }

    // The two reconstructions. Both are piecewise linear, so the polylines are
    // the reconstructions rather than a rendering of them.
    for i in lo_index..hi_index {
        chart.line(
            p(i as f64, line.f[i]),
            p((i + 1) as f64, line.f[i + 1]),
            STANDARD,
        );
        chart.line(
            p(i as f64 + TAU, line.c[i]),
            p((i + 1) as f64 + TAU, line.c[i + 1]),
            SHIFTED,
        );
    }
    for i in lo_index..=hi_index {
        chart.circle(
            Isometry3d::from_translation(p(i as f64, line.f[i])),
            0.045,
            STANDARD,
        );
        chart.circle(
            Isometry3d::from_translation(p(i as f64 + TAU, line.c[i])),
            0.035,
            SHIFTED,
        );
    }

    // `tau`, drawn to scale between one sample and its own knot.
    let y_tau = rect[1] + rect[3] * 0.10;
    let (x_a, x_b) = (x_of((line.n - 1) as f64), x_of((line.n - 1) as f64 + TAU));
    frames.line(
        Vec3::new(x_a, y_tau, 0.0),
        Vec3::new(x_b, y_tau, 0.0),
        SHIFTED,
    );
    for x in [x_a, x_b] {
        frames.line(
            Vec3::new(x, y_tau - 0.07, 0.0),
            Vec3::new(x, y_tau + 0.07, 0.0),
            SHIFTED,
        );
    }

    // The roots. Three of them land within a thousandth of a cell of each other
    // and this panel cannot separate them -- which is the reason the inset
    // exists, and is said rather than left for the viewer to wonder about.
    let mark = |g: &mut Gizmos<ChartGizmos>, u: f64, colour: Color, half: f32| {
        let x = x_of(u);
        g.line(
            Vec3::new(x, y_zero - half, 0.0),
            Vec3::new(x, y_zero + half, 0.0),
            colour,
        );
    };
    mark(chart, line.u_exact, EXACT, 0.60);
    match line.sabotage {
        Some(sabotage) if sabotage.u > u_lo && sabotage.u < u_hi => {
            mark(chart, sabotage.u, SABOTAGE, 0.42);
            labels.0.push(LabelSpec {
                at: Vec3::new(x_of(sabotage.u), y_zero - 0.48, 0.0),
                text: format!("no prefilter: {:.3} cells off, i.e. tau", sabotage.error),
                colour: SABOTAGE,
                justify: JustifyContent::Center,
                align: AlignItems::FlexStart,
                size: 11.0,
            });
        }
        Some(_) => {}
        None => {
            labels.0.push(LabelSpec {
                at: at(rect, 0.98, 0.06),
                text: "no prefilter: no root in the bracket at all, past sigma = 1 - tau".into(),
                colour: SABOTAGE,
                justify: JustifyContent::FlexEnd,
                align: AlignItems::FlexStart,
                size: 11.0,
            });
        }
    }

    labels.0.push(LabelSpec {
        at: at(rect, 0.012, 0.965),
        text: format!("THE LINE   {name}   bracket [{}, {}]", line.n - 1, line.n),
        colour: EXACT,
        justify: JustifyContent::FlexStart,
        align: AlignItems::FlexStart,
        size: 12.0,
    });
    labels.0.push(LabelSpec {
        at: at(rect, 0.012, 0.895),
        text: "samples f_n -- t = a/(a-b) reads its root off these".into(),
        colour: STANDARD,
        justify: JustifyContent::FlexStart,
        align: AlignItems::FlexStart,
        size: 11.0,
    });
    labels.0.push(LabelSpec {
        at: at(rect, 0.012, 0.835),
        text: "prefiltered c_n, on the knots at (n + tau) T".into(),
        colour: SHIFTED,
        justify: JustifyContent::FlexStart,
        align: AlignItems::FlexStart,
        size: 11.0,
    });
    labels.0.push(LabelSpec {
        at: Vec3::new(f32::midpoint(x_a, x_b), y_tau - 0.10, 0.0),
        text: "tau = 1/5 cell".into(),
        colour: SHIFTED,
        justify: JustifyContent::Center,
        align: AlignItems::FlexStart,
        size: 11.0,
    });
    labels.0.push(LabelSpec {
        at: Vec3::new(x_of(line.u_exact), y_zero + 0.66, 0.0),
        text: "three roots, a thousandth of a cell apart ->".into(),
        colour: EXACT,
        justify: JustifyContent::Center,
        align: AlignItems::FlexEnd,
        size: 11.0,
    });
}

/// The magnification, where the three roots separate.
fn draw_inset(
    line: &Line,
    chart: &mut Gizmos<ChartGizmos>,
    frames: &mut Gizmos<FrameGizmos>,
    labels: &mut LabelSpecs,
) {
    let rect = PANEL_INSET;
    let half = (line
        .error_standard
        .max(line.error_shifted)
        .max(INSET_FLOOR)
        * INSET_MARGIN)
        .max(INSET_FLOOR);
    let u_lo = line.u_exact - half;
    let u_hi = line.u_exact + half;
    let x_of = |u: f64| rect[0] + rect[2] * ((u - u_lo) / (u_hi - u_lo)) as f32;
    let y_zero = rect[1] + rect[3] * 0.46;

    // Both reconstructions are straight lines across a window this narrow, so
    // the vertical scale is the chord's own slope over the window -- which puts
    // the standard rule's crossing at a readable angle at every magnitude.
    let slope = (line.f[line.n] - line.f[line.n - 1]).abs().max(1e-300);
    let span = slope * half;
    let y_of = |v: f64| y_zero + rect[3] * 0.36 * (v / span).clamp(-1.0, 1.0) as f32;

    frames.line(
        Vec3::new(rect[0], y_zero, 0.0),
        Vec3::new(rect[0] + rect[2], y_zero, 0.0),
        AXIS,
    );

    // The standard rule is the chord across the bracket; the shifted rule is
    // whichever of the two shifted pieces carries the crossing. Both are drawn
    // from the functions themselves, evaluated at the window edges.
    let a = line.f[line.n - 1];
    let b = line.f[line.n];
    let chord = |u: f64| a + (b - a) * (u - (line.n - 1) as f64);
    chart.line(
        Vec3::new(x_of(u_lo), y_of(chord(u_lo)), 0.0),
        Vec3::new(x_of(u_hi), y_of(chord(u_hi)), 0.0),
        STANDARD,
    );
    if let (Some(v_lo), Some(v_hi)) = (
        reconstruct(&line.c, 0, u_lo),
        reconstruct(&line.c, 0, u_hi),
    ) {
        chart.line(
            Vec3::new(x_of(u_lo), y_of(v_lo), 0.0),
            Vec3::new(x_of(u_hi), y_of(v_hi), 0.0),
            SHIFTED,
        );
    }

    let rule = |g: &mut Gizmos<ChartGizmos>, u: f64, colour: Color| {
        let x = x_of(u).clamp(rect[0], rect[0] + rect[2]);
        g.line(
            Vec3::new(x, rect[1] + rect[3] * 0.30, 0.0),
            Vec3::new(x, rect[1] + rect[3] * 0.86, 0.0),
            colour,
        );
    };
    rule(chart, line.u_exact, EXACT);
    rule(chart, line.u_standard, STANDARD);
    rule(chart, line.u_shifted, SHIFTED);

    labels.0.push(LabelSpec {
        at: at(rect, 0.015, 0.965),
        text: format!("MAGNIFIED   +/- {half:.3e} cells"),
        colour: EXACT,
        justify: JustifyContent::FlexStart,
        align: AlignItems::FlexStart,
        size: 12.0,
    });
    // Three rows rather than three labels on one, because two of these roots
    // are routinely within a pixel of each other and stacked text is unreadable
    // exactly when the picture is most interesting.
    for (row, (u, colour, text)) in [
        (
            line.u_standard,
            STANDARD,
            format!("a/(a-b)   error {:.3e}", line.error_standard),
        ),
        (line.u_exact, EXACT, "exact, bisected".to_string()),
        (
            line.u_shifted,
            SHIFTED,
            format!("shifted   error {:.3e}", line.error_shifted),
        ),
    ]
    .into_iter()
    .enumerate()
    {
        let x = x_of(u).clamp(rect[0] + 0.05, rect[0] + rect[2] - 0.05);
        let y = rect[1] + rect[3] * (0.245 - 0.075 * row as f32);
        // A leader from the rule down to its own row, so which caption belongs
        // to which vertical is drawn rather than guessed.
        frames.line(
            Vec3::new(x, y, 0.0),
            Vec3::new(x, rect[1] + rect[3] * 0.30, 0.0),
            colour.with_alpha(0.5),
        );
        labels.0.push(LabelSpec {
            at: Vec3::new(x, y, 0.0),
            text,
            colour,
            justify: JustifyContent::Center,
            align: AlignItems::Center,
            size: 11.0,
        });
    }
}

/// The closed form, with the live measurement riding on it.
fn draw_ratio(
    live: &Live,
    base: &Baseline,
    chart: &mut Gizmos<ChartGizmos>,
    frames: &mut Gizmos<FrameGizmos>,
    labels: &mut LabelSpecs,
) {
    let rect = PANEL_RATIO;
    let top = 1.25f64;
    let x_of = |sigma: f64| rect[0] + rect[2] * (0.10 + 0.86 * sigma) as f32;
    let y_of = |ratio: f64| rect[1] + rect[3] * (0.13 + 0.76 * (ratio / top).clamp(0.0, 1.2)) as f32;
    let p = |sigma: f64, ratio: f64| Vec3::new(x_of(sigma), y_of(ratio), 0.0);

    frames.line(p(0.0, 0.0), p(1.0, 0.0), AXIS);
    frames.line(p(0.0, 0.0), p(0.0, top), AXIS);
    for ratio in [0.5, 1.0] {
        frames.line(
            p(0.0, ratio),
            p(1.0, ratio),
            AXIS.with_alpha(if ratio == 1.0 { 0.55 } else { 0.25 }),
        );
    }
    // The registered bar. C1 asked for 30% lower, which is this line.
    frames.line(p(0.0, REGISTERED_BAR), p(1.0, REGISTERED_BAR), BAR);
    for sigma in [TAU, 2.0 * TAU] {
        frames.line(p(sigma, 0.0), p(sigma, top), CLOSED_FORM.with_alpha(0.35));
    }
    // Where the closed form clears the bar: sigma <= 1/17 and sigma >= 4/17.
    for band in [(0.0, 1.0 / 17.0), (4.0 / 17.0, 1.0)] {
        for step in 0..3 {
            let y = -0.035 * f64::from(step) - 0.02;
            chart.line(p(band.0, y), p(band.1, y), CLOSED_FORM.with_alpha(0.7));
        }
    }

    // The closed form itself, behind everything measured.
    const STEPS: usize = 260;
    let mut previous: Option<Vec3> = None;
    for i in 1..=STEPS {
        let sigma = i as f64 / (STEPS + 1) as f64;
        let point = p(sigma, predicted_ratio(sigma).min(top));
        if let Some(last) = previous {
            chart.line(last, point, CLOSED_FORM);
        }
        previous = Some(point);
    }

    // The trail, then the live point on top of it.
    for pair in live.trail.windows(2) {
        chart.line(
            p(f64::from(pair[0].x), f64::from(pair[0].y).min(top)),
            p(f64::from(pair[1].x), f64::from(pair[1].y).min(top)),
            MEASURED.with_alpha(0.75),
        );
    }
    if let Some(line) = &live.line
        && let Some(predicted) = line.predicted()
    {
        let sigma = line.sigma_true();
        let point = p(sigma, line.ratio.min(top));
        chart.line(p(sigma, 0.0), point, MEASURED.with_alpha(0.45));
        chart.circle(Isometry3d::from_translation(point), 0.085, MEASURED);
        chart.circle(Isometry3d::from_translation(point), 0.125, MEASURED);
        labels.0.push(LabelSpec {
            at: point + Vec3::new(0.0, 0.16, 0.0),
            text: format!(
                "measured {:.6}   closed form {:.6}",
                line.ratio, predicted
            ),
            colour: MEASURED,
            justify: JustifyContent::Center,
            align: AlignItems::FlexEnd,
            size: 11.0,
        });
    }

    labels.0.push(LabelSpec {
        at: at(rect, 0.015, 0.965),
        text: "THE RATIO   shifted error / standard error".into(),
        colour: EXACT,
        justify: JustifyContent::FlexStart,
        align: AlignItems::FlexStart,
        size: 12.0,
    });
    labels.0.push(LabelSpec {
        at: at(rect, 0.015, 0.900),
        text: "closed form  | sigma - 2 tau | / sigma".into(),
        colour: CLOSED_FORM,
        justify: JustifyContent::FlexStart,
        align: AlignItems::FlexStart,
        size: 11.0,
    });
    labels.0.push(LabelSpec {
        at: p(1.0, REGISTERED_BAR),
        text: format!("{REGISTERED_BAR:.2}  C1's registered bar"),
        colour: BAR,
        justify: JustifyContent::FlexEnd,
        align: AlignItems::FlexEnd,
        size: 11.0,
    });
    for (sigma, text) in [(TAU, "sigma = tau"), (2.0 * TAU, "sigma = 2 tau")] {
        labels.0.push(LabelSpec {
            at: p(sigma, top),
            text: text.into(),
            colour: CLOSED_FORM,
            justify: JustifyContent::Center,
            align: AlignItems::FlexEnd,
            size: 11.0,
        });
    }
    for (ratio, text) in [(0.0, "0"), (1.0, "1")] {
        labels.0.push(LabelSpec {
            at: p(0.0, ratio) - Vec3::new(0.06, 0.0, 0.0),
            text: text.into(),
            colour: AXIS,
            justify: JustifyContent::FlexEnd,
            align: AlignItems::Center,
            size: 11.0,
        });
    }
    for (sigma, text) in [(0.0, "sigma 0"), (1.0, "1")] {
        labels.0.push(LabelSpec {
            at: p(sigma, -0.16),
            text: text.into(),
            colour: AXIS,
            justify: JustifyContent::Center,
            align: AlignItems::FlexStart,
            size: 11.0,
        });
    }
    labels.0.push(LabelSpec {
        at: p(0.62, -0.10),
        text: format!(
            "clears the bar on 14/17 = {:.1}%, measured {:.1}%",
            CLEAR_FRACTION_EXACT * 100.0,
            base.sweep.clear_fraction * 100.0
        ),
        colour: CLOSED_FORM,
        justify: JustifyContent::Center,
        align: AlignItems::FlexStart,
        size: 11.0,
    });
}

/// The guard band: the delta from the whole-line root, on a log axis.
fn draw_guard(
    live: &Live,
    base: &Baseline,
    signals: &Signals,
    probe: usize,
    chart: &mut Gizmos<ChartGizmos>,
    labels: &mut LabelSpecs,
) {
    let rect = PANEL_GUARD;
    let y_of = |delta: f64| {
        let log = if delta > 0.0 {
            delta.log10().clamp(GUARD_LOG[0], GUARD_LOG[1])
        } else {
            GUARD_LOG[0]
        };
        rect[1]
            + rect[3]
                * (0.14 + 0.72 * ((log - GUARD_LOG[0]) / (GUARD_LOG[1] - GUARD_LOG[0])) as f32)
    };
    let floor = y_of(0.0);

    for decade in -8..=-1 {
        let y = y_of(10f64.powi(decade));
        chart.line(
            Vec3::new(rect[0] + 0.30, y, 0.0),
            Vec3::new(rect[0] + rect[2] - 0.10, y, 0.0),
            AXIS.with_alpha(0.20),
        );
        if decade % 2 != 0 {
            labels.0.push(LabelSpec {
                at: Vec3::new(rect[0] + 0.26, y, 0.0),
                text: format!("1e{decade}"),
                colour: AXIS,
                justify: JustifyContent::FlexEnd,
                align: AlignItems::Center,
                size: 10.0,
            });
        }
    }
    // C3's threshold.
    let y_bar = y_of(GUARD_TOLERANCE);
    chart.line(
        Vec3::new(rect[0] + 0.30, y_bar, 0.0),
        Vec3::new(rect[0] + rect[2] - 0.10, y_bar, 0.0),
        BAR,
    );
    labels.0.push(LabelSpec {
        at: Vec3::new(rect[0] + rect[2] - 0.10, y_bar, 0.0),
        text: "1e-6".into(),
        colour: BAR,
        justify: JustifyContent::FlexEnd,
        align: AlignItems::FlexEnd,
        size: 11.0,
    });

    let empty: [Guard; 4] = [Guard::default(); 4];
    let live_guards = live.line.as_ref().map_or(&empty, |l| &l.guards);
    let groups: [(&str, &[Guard; 4]); 3] = [
        (signals.name(probe), live_guards),
        ("torus", base.guards("torus").unwrap_or(&empty)),
        ("gyroid", base.guards("gyroid").unwrap_or(&empty)),
    ];

    let left = rect[0] + 0.72;
    let usable = rect[2] - 0.95;
    let group_w = usable / groups.len() as f32;
    for (gi, (name, guards)) in groups.iter().enumerate() {
        let base_x = left + group_w * gi as f32;
        for (bi, guard) in guards.iter().enumerate() {
            let x = base_x + group_w * (0.14 + 0.22 * bi as f32);
            let y = y_of(guard.delta);
            let colour = if !guard.truncated {
                AXIS.with_alpha(0.55)
            } else if guard.converged {
                CLOSED_FORM
            } else {
                SHIFTED
            };
            // Three parallel strokes rather than one, because a gizmo line has
            // no width in world units and a single stroke reads as a hair.
            for offset in [-0.035f32, 0.0, 0.035] {
                chart.line(
                    Vec3::new(x + offset, floor, 0.0),
                    Vec3::new(x + offset, y, 0.0),
                    colour,
                );
            }
            // A cap on the axis, so a delta of exactly zero is a bar that was
            // measured and came back at nothing rather than a bar that is
            // missing. Six of the eight fields read exactly zero at every `k`.
            chart.line(
                Vec3::new(x - 0.055, floor, 0.0),
                Vec3::new(x + 0.055, floor, 0.0),
                colour,
            );
            labels.0.push(LabelSpec {
                at: Vec3::new(x, floor - 0.05, 0.0),
                text: format!("{}", guard.k),
                colour: if guard.truncated {
                    AXIS
                } else {
                    AXIS.with_alpha(0.5)
                },
                justify: JustifyContent::Center,
                align: AlignItems::FlexStart,
                size: 10.0,
            });
        }
        labels.0.push(LabelSpec {
            at: Vec3::new(base_x + group_w * 0.44, floor - 0.30, 0.0),
            text: (*name).to_string(),
            colour: if gi == 0 { MEASURED } else { EXACT },
            justify: JustifyContent::Center,
            align: AlignItems::FlexStart,
            size: 11.0,
        });
    }

    labels.0.push(LabelSpec {
        at: at(rect, 0.015, 0.965),
        text: "THE GUARD BAND   the k preceding samples only".into(),
        colour: EXACT,
        justify: JustifyContent::FlexStart,
        align: AlignItems::FlexStart,
        size: 12.0,
    });
    labels.0.push(LabelSpec {
        at: at(rect, 0.015, 0.900),
        text: "delta from the whole-line root, in cells; grey = whole line".into(),
        colour: AXIS,
        justify: JustifyContent::FlexStart,
        align: AlignItems::FlexStart,
        size: 10.0,
    });
    labels.0.push(LabelSpec {
        at: at(rect, 0.50, 0.035),
        text: "k = 5 fails the bar, k = 10 passes on 8 of 8 -- C3 HELD".into(),
        colour: CLOSED_FORM,
        justify: JustifyContent::Center,
        align: AlignItems::FlexStart,
        size: 11.0,
    });
}

/// Project every label through the camera and move its node there.
///
/// Projected rather than pinned to fixed pixels: the chart is world geometry and
/// a viewer who orbits it should not watch the labels stay behind. A point the
/// camera cannot see hides its label instead of clamping it to an edge, where it
/// would read as data.
fn place_labels(
    specs: Res<LabelSpecs>,
    camera: Query<(&Camera, &GlobalTransform)>,
    mut nodes: Query<(&ChartLabel, &mut Node, &Children)>,
    mut text: Query<(&mut Text, &mut TextColor, &mut TextFont)>,
) {
    let Ok((camera, camera_transform)) = camera.single() else {
        return;
    };
    for (label, mut node, children) in &mut nodes {
        let placed = specs
            .0
            .get(label.0)
            .and_then(|spec| {
                camera
                    .world_to_viewport(camera_transform, spec.at)
                    .ok()
                    .map(|at| (spec, at))
            })
            .map(|(spec, at)| {
                node.display = Display::Flex;
                node.left = Val::Px(at.x);
                node.top = Val::Px(at.y);
                node.justify_content = spec.justify;
                node.align_items = spec.align;
                for child in children {
                    if let Ok((mut value, mut colour, mut font)) = text.get_mut(*child) {
                        if value.0 != spec.text {
                            value.0.clone_from(&spec.text);
                        }
                        colour.0 = spec.colour;
                        font.font_size = FontSize::Px(spec.size);
                    }
                }
            });
        if placed.is_none() {
            node.display = Display::None;
        }
    }
}

// ─── what is on screen in words ─────────────────────────────────────────────

fn report(
    live: Res<Live>,
    shot: Res<Shot>,
    base: Res<Baseline>,
    signals: Res<Signals>,
    mut stats: ResMut<DemoStats>,
    mut caption: Query<&mut Text, With<Caption>>,
) {
    stats.title =
        String::from("E-316 shifted linear root -- M-359 / x42 (P-60, R-058), tau = 1/5");
    // No extractor is called anywhere in P-60, so these three are zero by
    // construction rather than by omission, and the first HUD line says so.
    stats.vertices = 0;
    stats.triangles = 0;
    stats.extract_ms = 0.0;

    let name = signals.name(shot.probe);
    let mut extra = vec![String::from(
        "no extractor is called -- P-60 reads one 65-sample line -- so the three counts",
    )];
    extra.push(String::from(
        "above are 0 by construction, and no golden hash is reachable from here.",
    ));

    if let Some(line) = &live.line {
        let predicted = line.predicted();
        extra.push(format!(
            "\nthe line   {name}   T {:.7}   bracket [{}, {}]   shifted piece: {}",
            signals.cell(shot.probe),
            line.n - 1,
            line.n,
            line.segment
        ));
        extra.push(format!(
            "  exact    bisected    {:>15.9}",
            line.u_exact
        ));
        extra.push(format!(
            "  standard a/(a-b)     {:>15.9}   error {}",
            line.u_standard,
            sci(line.error_standard)
        ));
        extra.push(format!(
            "  shifted  tau = 1/5   {:>15.9}   error {}",
            line.u_shifted,
            sci(line.error_shifted)
        ));
        extra.push(match predicted {
            Some(p) => format!(
                "  ratio {}   closed form at sigma {:.4} = {}   {:.2}% apart",
                sci(line.ratio),
                line.sigma_true(),
                sci(p),
                (line.ratio - p).abs() / p * 100.0
            ),
            None => format!(
                "  ratio {}   DEGENERATE: the root is on a sample, both at the floor",
                sci(line.ratio)
            ),
        });
        extra.push(format!(
            "  tau c(n-1) + (1-tau) c_n == f_n to {} -- the filter interpolates",
            sci(line.interpolation_residual)
        ));
    } else {
        extra.push(format!(
            "\nthe line   {name} has no usable crossing at slide {:+.3}",
            shot.slide
        ));
    }

    extra.push(format!(
        "\nthe sweep  sigma {:.4}   {}",
        shot.slide,
        match shot.beat {
            Beat::Control => "holding on the committed positive control",
            Beat::Sweep => "sliding the FIELD along a fixed grid",
            Beat::Zero => "settling on sigma = 2 tau, where the shifted root is exact",
            Beat::Manual => "driven by hand -- [R] hands it back to the loop",
            Beat::Committed => "a reference field, pinned at the committed configuration",
        }
    ));
    extra.push(String::from(
        "  the CLOSED FORM is exactly 1 at sigma = tau, exactly 0 at 2 tau, never above",
    ));
    extra.push(format!(
        "  1, and clears C1's 0.70 bar on 14/17 = {:.2}% of uniform crossings.",
        CLEAR_FRACTION_EXACT * 100.0
    ));
    extra.push(format!(
        "  MEASURED over {} positions: median {:.2}% off it, peak {:.4} at the kink,",
        base.sweep.positions,
        base.sweep.median_rel * 100.0,
        base.sweep.max_measured
    ));
    extra.push(format!(
        "  clearing on {:.2}%. the gap is the O(g'' squared) term the form drops.",
        base.sweep.clear_fraction * 100.0
    ));

    let guard_row = |name: &str| {
        base.guards(name).map_or_else(
            || String::from("missing"),
            |g| {
                g.iter()
                    .map(|guard| {
                        if guard.truncated {
                            format!("k{:<2} {}", guard.k, sci(guard.delta))
                        } else {
                            format!("k{:<2} whole line", guard.k)
                        }
                    })
                    .collect::<Vec<_>>()
                    .join("  ")
            },
        )
    };
    extra.push(String::from(
        "\nthe guard band   delta from the whole-line root, in cells, bar at 1e-6",
    ));
    extra.push(format!("  torus   {}", guard_row("torus")));
    extra.push(format!("  gyroid  {}", guard_row("gyroid")));
    extra.push(String::from(
        "  k = 5 fails on both fields that carry a delta; k = 10 passes on 8 of 8.",
    ));

    extra.push(format!(
        "\ncontrols  positive  exact quadratic at sigma = 1/2   {} vs a closed-form",
        sci(base.control.ratio)
    ));
    extra.push(format!(
        "                    2.000000e-1.  negative  raw samples, no prefilter: {},",
        sci(base.control.sabotage.map_or(f64::NAN, |s| s.ratio))
    ));
    extra.push(format!(
        "                    {:.0}x worse, root biased {:.6} cells = tau, and above",
        base.control.sabotage.map_or(f64::NAN, |s| s.ratio),
        base.control.sabotage.map_or(f64::NAN, |s| s.error)
    ));
    extra.push(String::from(
        "                    sigma = 1 - tau it has no root in the bracket at all.",
    ));

    let agreed = base.agreed();
    extra.push(format!(
        "\ndocs/experiments/p-60.csv  {} rows  {} of {} quoted numbers agree as {{:.6e}}",
        base.ledger_rows,
        agreed,
        base.checks.len()
    ));
    for (label, index) in [
        ("median_error_ratio", 0usize),
        ("torus guard k = 10", 5),
        ("control quadratic ", 7),
    ] {
        if let Some(check) = base.checks.get(index) {
            extra.push(format!(
                "  {label}   expected {:>12}   measured {:>12}",
                check.expected, check.measured
            ));
        }
    }
    extra.push(format!(
        "  {} of {} fields sit at the BISECTION FLOOR: the domain-centre axis is a",
        base.exact_count,
        base.fields.len()
    ));
    extra.push(format!(
        "  symmetry axis, sphere restricts to |x| - 1, so C1's median came out {}.",
        sci(base.median_smooth)
    ));
    for check in base.checks.iter().filter(|c| !c.agrees()) {
        extra.push(format!(
            "CROSS-CHECK FAILED  {}  expected {}  measured {}",
            check.label, check.expected, check.measured
        ));
    }
    extra.push(String::from(
        "\nNOT SHIPPABLE as accuracy: a mesher cannot choose sigma, so the benefit is the",
    ));
    extra.push(String::from(
        "82% and not the 8 dB; edge_position is unchanged. The guard band ships: 10.",
    ));
    extra.push(String::from(
        "[<-][->] sigma by hand   [up][down] or [1]-[9] which signal   [R] back to the loop",
    ));
    stats.extra = extra;

    let text = caption_for(&live, &shot, &base, name);
    for mut target in &mut caption {
        target.0.clone_from(&text);
    }
}

/// The line a viewer reads instead of the HUD.
fn caption_for(live: &Live, shot: &Shot, base: &Baseline, name: &str) -> String {
    let Some(line) = &live.line else {
        return format!("{name} has no usable crossing here");
    };
    match shot.beat {
        Beat::Control => format!(
            "the committed positive control -- an exact quadratic crossing at sigma = 1/2: \
             ratio {}, and p-60.csv says {}",
            sci(line.ratio),
            base.checks
                .iter()
                .find(|c| c.label.starts_with("control_quadratic"))
                .map_or("?", |c| c.expected.as_str())
        ),
        Beat::Sweep => format!(
            "sliding the field, never the grid: at sigma = {:.3} the shifted root is {} of \
             the standard one",
            line.sigma_true(),
            sci(line.ratio)
        ),
        Beat::Zero => format!(
            "at sigma = 2 tau the shifted root is EXACT -- error {} against {} -- and at \
             sigma = tau it is exactly as wrong",
            sci(line.error_shifted),
            sci(line.error_standard)
        ),
        Beat::Manual => format!(
            "sigma = {:.4} by hand   ratio {}   [R] hands it back to the loop",
            line.sigma_true(),
            sci(line.ratio)
        ),
        Beat::Committed => format!(
            "{name} at the committed configuration: ratio {}, and p-60.csv agrees on {} of {} \
             quoted numbers",
            sci(line.ratio),
            base.agreed(),
            base.checks.len()
        ),
    }
}
