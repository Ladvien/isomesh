//! **P-77 - how much temporal history a dig destroys, and whether 0.2 ms buys it back.**
//!
//! Ticket: R-077. Pre-registered before this harness existed.
//!
//! ```bash
//! cargo bench --bench experiment_p77
//! ```
//!
//! Writes `docs/experiments/p-77.csv`.
//!
//! # Hypothesis, as registered
//!
//! The one measured mitigation in this area is cheap and the problem it
//! mitigates has never been measured under destruction. k-DOP clipping
//! (Ikkala, Lauttia, Jaaskelainen & Makitalo, SIGGRAPH Asia 2024 Technical
//! Communications, `10.1145/3681758.3697996`) replaces TAA's axis-aligned
//! neighbourhood clamp with a tighter k-discrete-oriented polytope at 0.2 ms
//! overhead.
//!
//! - **C1.** The fraction of TAA history samples rejected in the frame after a
//!   brush stroke is at least **5x** the steady-state rejection rate, and the
//!   elevated rate persists for at least **3 frames**.
//! - **C2.** k-DOP clipping recovers at least **half** the rejected samples at
//!   under **0.3 ms**.
//! - **C3.** The rejection is spatially concentrated at the brush: over **80%**
//!   of rejected samples fall within the brush's screen-space bounding box
//!   dilated by one dig radius.
//!
//! VACUITY CONTROL, as registered: the steady-state arm must have a **non-zero**
//! rejection rate, or the 5x ratio is division by a floor. Asserted below.
//!
//! # The SHARE line, recomputed before the code was written
//!
//! **C1 is reachable but bounded.** It is a rate rather than a ratio of a total,
//! so the only arithmetic ceiling is `post <= 1.0`: a 5x ratio needs a
//! steady-state rate under 0.20. The measured steady-state rate is reported and
//! is far below that, so the clause is reachable.
//!
//! **C2's recovery half is arithmetically UNREACHABLE, and the paper says so.**
//! Section 3 of `10.1145/3681758.3697996`: *"As a general-purpose approach, we
//! include the X, Y and Z axes and optimize the rest of the axes to achieve
//! tightest bounds around a unit sphere. This approach **should never be looser
//! than an AABB**."* A k-DOP that contains the AABB's own axes is a subset of
//! that AABB, so every colour the AABB rejects the k-DOP also rejects. The set
//! of samples a k-DOP "recovers" from an AABB in the same colour space is
//! **empty by construction**, and no measurement can make it half. The paper's
//! own framing is the other direction - Figure 1's caption is *"neighborhood
//! clipping with an AABB often **allows** colors which don't fit"* - so k-DOP
//! clipping reduces ghosting by rejecting **more**, never by recovering. C2 as
//! registered inverts its source. This harness runs it anyway (the `x51` rule)
//! and reports the number, plus the quantity the paper actually claims.
//!
//! **C2's cost half is about a GPU TAA resolve that does not exist here.** The
//! source's own numbers, which the registration paraphrases from the abstract:
//! the 0.2 ms is a **GTX 1080 Ti** (Section 4), and on an **RTX 3090** the same
//! delta is 210 - 160 = **50 us**. The GPU is named in the body even though it is
//! not named in the abstract. `kdop_ms` here is the wall-clock of *this
//! harness's own single-threaded CPU clamp pass* over its own pixel population.
//! It is labelled as such and is not a GPU TAA resolve cost.
//!
//! **C3 is not independent of C1, and that is provable rather than measured.**
//! Subtracting a sphere is `f' = max(f, -(|p-c| - r))`. Outside the sphere the
//! second term is negative, so a point that was solid stays solid and a point
//! that was air is untouched: **the zero set changes only inside the brush**.
//! A ray whose depth changed therefore passed through the brush, so every
//! changed pixel lies inside the brush's projected silhouette disc. The harness
//! asserts exactly that. Consequently, writing `b` for the dilated box's area
//! fraction and `R` for C1's ratio, the concentration is `(R - 1 + b)/R` - so
//! C3's 80% is C1's `R >= 5` restated, up to the box's own area. Both numbers
//! are reported; a reader should not treat them as two independent tests.
//!
//! # The instrument, in enough detail to judge it
//!
//! `crates/isomesh` must not depend on Bevy, so there is no renderer here. TAA
//! history rejection is nonetheless a **per-pixel geometric predicate** and the
//! predicate is implemented in full:
//!
//! 1. **The field is `game_dig`'s.** `Ground`'s height field copied verbatim
//!    from `bevy_isomesh/examples/game_dig.rs`:
//!    `0.35*sin(0.9x)*cos(0.7z) + 0.15*sin(2.1x)`, negative below. Brushes are
//!    subtracted spheres of radius **0.25** (`game_dig`'s default `world.radius`).
//! 2. **The camera is `game_dig`'s.** `setup` inserts
//!    `Transform::from_xyz(0.0, 1.70, 6.0)` and `main` inserts
//!    `Look { yaw: 0.0, pitch: -0.15 }`; the walk speed is **2.5** units/s. The
//!    projection is Bevy's default `PerspectiveProjection`: a **45 degree**
//!    vertical field of view. A yaw rate of 0.25 rad/s is added, because a pure
//!    forward walk puts the focus of expansion at the centre of the frame and
//!    leaves the reprojection almost stationary exactly where the brush is - a
//!    fixture that would report an artificially low steady-state rate.
//! 3. **The depth buffer is `game_dig`'s own ray caster.** [`trace`] is a
//!    transcription of `game_dig::trace`: sphere tracing with `AIM_NEAR = 0.30`,
//!    `AIM_FAR = 25.0`, `AIM_STEPS = 128`, `AIM_HIT = 0.01`,
//!    `LIPSCHITZ = 1.25`, and the same sandbox box test. Tracing the field
//!    rather than the extracted mesh is what the demo does and for the demo's
//!    reason: the field is the thing being edited, so it cannot go stale.
//! 4. **Motion vectors are exact, from the camera and the geometry.** For a
//!    pixel with a hit at world point `X`, the previous-frame screen position is
//!    `X` projected through the previous frame's view - including the previous
//!    frame's jitter. The geometry is static between frames except where a brush
//!    landed, so this is the exact motion vector a renderer would write, not an
//!    approximation of one.
//! 5. **The per-pixel signal is a three-channel radiance built only from
//!    geometry** - a warm sun term in `max(0, n.l)`, a cool sky term in
//!    `n.y`, a view-dependent rim term, and distance fog. No textures and no
//!    shadows: the three channels have to be *independent* functions of the
//!    surface, because a grey signal would put every neighbourhood on the
//!    diagonal of colour space, where an AABB is a cube and any oriented volume
//!    wins for free. Distance fog is in because a hole that reveals a farther
//!    surface must move the signal, which is the whole mechanism under test.
//! 6. **The resolve is Karis's.** Jitter is a 10-entry R2 sequence (the paper's
//!    own pattern, Roberts 2018). History is fetched bilinearly at the
//!    reprojected position, rectified by **clipping** - a ray cast from the
//!    current pixel's colour to the history colour, stopped at the bounding
//!    volume's shell, which is the operation Section 3 of the paper describes -
//!    and blended `history = lerp(clipped, current, 0.1)`. A sample is
//!    **rejected** when the clip moved it, i.e. when the reprojected colour lay
//!    outside the neighbourhood volume. This is the registered quantity.
//! 7. **Four bounding volumes are evaluated on the same fetched history**, so
//!    the comparison is per-sample rather than between two pipelines that have
//!    drifted apart: `aabb_rgb` (Lottes 2011), `aabb_ycocg` (Karis 2014, the
//!    baseline, because that is what production TAA does), `dop26_ycocg` and
//!    `dop26_rgb` - Klosowski et al. 1998's 13-axis 26-DOP, three coordinate
//!    axes plus four body diagonals plus six edge diagonals, ranges dilated by
//!    the paper's `eps = 1e-5`. The paper's 32-DOP has 16 axes rather than 13;
//!    adding axes only tightens the volume further, and tightness is what
//!    decides C2, so 13 is a conservative stand-in and not a weaker one.
//!    `dop26_ycocg` contains the YCoCg axes and is therefore a subset of
//!    `aabb_ycocg`; `dop26_rgb` does not, and is the control that proves the
//!    recovery counter can report a non-zero.
//! 8. **Control and treatment run in lockstep.** Every arm simulates two worlds
//!    over one camera path and one jitter sequence: one that is never dug and
//!    one that is. Rejection is therefore compared frame-paired, so a spike
//!    cannot be camera motion in disguise.
//!
//! # Controls, each an assertion rather than a printed number
//!
//! - **The registered vacuity control.** The steady-state arm's rejection count
//!   must be non-zero in the measurement window, or the ratio is division by a
//!   floor.
//! - **The population must be non-empty.** `game_dig`'s opening pitch of -0.15
//!   rad from an eye at 1.70 puts the horizon inside the frame, and this is the
//!   control that caught it: a level camera - the `Transform` alone, before
//!   `Look` is applied - never meets a height field whose maximum is 0.50, and
//!   the first version of this harness had an empty depth buffer.
//! - **Changed pixels must lie inside the brush's silhouette.** The algebra
//!   above says the zero set moves only inside the brush; the harness projects
//!   the brush's silhouette cone and asserts that no pixel outside it changed.
//!   A wrong projection, a wrong brush position or a leaky `Dug` would all show
//!   up here as a large count.
//! - **The clamp timings must be comparable.** One build, one run, one thread
//!   for the clamp passes, and the ratio is reported beside the millisecond
//!   figure (`M-280`, `M-281`).
//!
//! # Two fixture defects this harness's own controls found, in order
//!
//! Both are recorded because they are the most transferable part of the row:
//! `R-076` and `R-091` are gated on this instrument and will hit the same two.
//!
//! 1. **A smooth signal makes the neighbourhood a needle.** Run 1 shaded the
//!    surface with normals and fog alone and measured a steady-state rejection
//!    rate of **95.5%** under a walking camera. The 3x3 AABB over nine samples
//!    of a smooth gradient spans one pixel of that gradient, so any sub-pixel
//!    reprojection error leaves it. Fixed by [`DETAIL_TILE`]'s albedo detail,
//!    which is what `game_dig`'s own textured terrain supplies.
//! 2. **The rejection rate is a function of reprojection displacement and
//!    almost nothing else, and a fast camera saturates it.** With the detail in,
//!    the walking arm still reads **86.6%** - and `frame_mean_reproj_px` says
//!    why: **7.2 pixels per frame** at 960x540. A player walking at
//!    `game_dig`'s 2.5 units/s while looking at ground 3.0 units away moves the
//!    image 9.5 pixels of parallax per 60 Hz frame; a 3x3 neighbourhood can
//!    validate history that moved one. That is a fact about locomotion, not
//!    about digging, and a ratio taken between two saturated rates measures
//!    nothing. So the fixture became a **sweep in reprojection displacement**
//!    (0.61, 2.3, 4.0, 7.2 px) and the registered columns are reported from
//!    [`HEADLINE_ARM`], the only arm where the denominator has headroom.

#![allow(
    clippy::cast_possible_truncation,
    clippy::cast_precision_loss,
    clippy::cast_sign_loss,
    clippy::print_literal,
    clippy::too_many_lines
)]

mod common;

use std::time::Instant;

use isomesh::Sdf;

// ---------------------------------------------------------------------------
// game_dig's constants, restated with the source named.
//
// Every one of these is copied from `bevy_isomesh/examples/game_dig.rs`.
// ---------------------------------------------------------------------------

/// `game_dig::EDIT_PERIOD` is `0.08` s, which the demo's own doc comment calls
/// **12.5 edits a second**. This is that rate.
const STROKE_RATE_HZ: f32 = 12.5;
/// Seconds between strokes while the button is held: `game_dig::EDIT_PERIOD`.
const EDIT_PERIOD: f32 = 0.08;
/// `game_dig`'s default `World::radius`, the brush the demo opens with.
const BRUSH_RADIUS: f32 = 0.25;
/// `game_dig::LIPSCHITZ`. `Ground` is not 1-Lipschitz; this is the bound that
/// makes sphere tracing sound over the whole edited field.
const LIPSCHITZ: f32 = 1.25;
/// `game_dig::AIM_NEAR`.
const AIM_NEAR: f32 = 0.30;
/// `game_dig::AIM_FAR`.
const AIM_FAR: f32 = 25.0;
/// `game_dig::AIM_STEPS`.
const AIM_STEPS: u32 = 128;
/// `game_dig::AIM_HIT`.
const AIM_HIT: f32 = 0.01;
/// `game_dig`'s eye position at `setup`: `Transform::from_xyz(0.0, 1.70, 6.0)`.
const EYE_START: [f32; 3] = [0.0, 1.70, 6.0];
/// `game_dig`'s opening `Look::pitch`.
const PITCH: f32 = -0.15;
/// `game_dig`'s unmodified walk speed (`ShiftLeft` raises it to 9.0).
const WALK_SPEED: f32 = 2.5;
/// The sandbox `game_dig::sandbox` computes: `CHUNK_CELLS * CELL_SIZE = 2.0`
/// units per chunk, `EXTENT = [8, 4, 8]` chunks, origin
/// `(-8.0, -5.4, -8.0)`.
const SANDBOX_LO: [f32; 3] = [-8.0, -5.4, -8.0];
/// Upper corner of the same box: `lo + (16, 8, 16)`.
const SANDBOX_HI: [f32; 3] = [8.0, 2.6, 8.0];

// ---------------------------------------------------------------------------
// This harness's own knobs.
// ---------------------------------------------------------------------------

/// Yaw rate, rad/s. Not a `game_dig` constant: the demo's yaw comes from the
/// mouse and has no rate. See the module docs for why a non-zero one is
/// necessary rather than convenient.
const YAW_RATE: f32 = 0.25;
/// Simulated frame time. 60 Hz, so `EDIT_PERIOD` is one stroke every 4.8 frames.
const DT: f32 = 1.0 / 60.0;
/// Frames of TAA history built before anything is measured. The paper warms 100;
/// convergence here is geometric at `1 - 0.1` per frame, so 24 frames leaves
/// `0.9^24 = 8%` of the initial transient, and the control arm carries whatever
/// remains equally.
const WARMUP: usize = 24;
/// Frames in the measurement window, starting with the frame the first brush
/// lands in.
const WINDOW: usize = 12;
/// TAA blend weight for the current frame. Karis's resolve keeps ~0.9 of the
/// history.
const ALPHA: f32 = 0.1;
/// Vertical field of view, radians. Bevy's `PerspectiveProjection::default`.
const FOV_Y: f32 = core::f32::consts::FRAC_PI_4;
/// Central-difference half-step for the shading normal.
///
/// **Not `game_dig::GRADIENT_EPS = 1e-4`, and that is a deviation with a
/// reason.** A central difference of an `f32` field of magnitude ~1 at `1e-4`
/// leaves about three significant digits in the gradient, and a shading signal
/// built from that noise would make the steady-state rejection rate a
/// measurement of float differencing rather than of geometry. `1e-3` is still
/// 1/125 of `game_dig`'s 0.125 cell.
const NORMAL_EPS: f32 = 1e-3;
/// How much elevation counts as elevated, for `frames_elevated`. Reported
/// per-frame as well, so a reader can apply another threshold.
const ELEVATION_FACTOR: f64 = 2.0;
/// Absolute signal difference above which a pixel counts as changed by the dig.
const CHANGE_TOL: f32 = 1e-4;
/// Slack, in pixels, on the brush-silhouette containment assertion. A ray
/// grazing the silhouette is decided by `f32` arithmetic in two different
/// places (the trace and the projection), and two pixels is far below the
/// 14-pixel silhouette this fixture produces.
const SILHOUETTE_SLACK: f32 = 2.0;
/// The arm the registered columns' headline values are taken from.
///
/// **Standing still, looking down at the rock under your feet, digging.** Not
/// the walking arm, and the reason is not that the walking arm gives a worse
/// number - every arm is falsified on C1 - but that the walking arm's
/// steady-state rejection rate is **86.6%**, so its ratio is the quotient of
/// two saturated numbers and cannot be a measurement of anything. At 7.2 pixels
/// of reprojection displacement per frame no 3x3 neighbourhood can validate a
/// history sample, and that is a fact about locomotion rather than about
/// digging. This arm's steady rate is 0.49%, which leaves 200x of headroom for
/// the registered 5x, so the falsification here is a measurement rather than a
/// division of two ceilings. Standing still to dig into rock is also
/// `game_dig`'s default state: the walk is key-driven and the demo opens with
/// no key held.
const HEADLINE_ARM: &str = "dig_at_feet_static";

/// Sun direction: `[0.35, 0.85, 0.40]` normalised.
const SUN_DIR: [f32; 3] = [0.349_128_2, 0.847_882_8, 0.399_003_6];
/// Warm direct term.
const ALBEDO_SUN: [f32; 3] = [1.00, 0.92, 0.78];
/// Cool hemispheric term.
const ALBEDO_SKY: [f32; 3] = [0.22, 0.30, 0.45];
/// View-dependent rim term, which is what makes the third channel independent.
const ALBEDO_RIM: [f32; 3] = [0.10, 0.10, 0.12];
/// Fog colour.
const FOG_RGB: [f32; 3] = [0.52, 0.58, 0.66];
/// Fog density per unit. Distance fog is in so that a depth change alone moves
/// the signal.
const FOG_DENSITY: f32 = 0.06;

/// World units per albedo tile: `game_dig`'s `TriplanarExtension::settings.x`.
///
/// **The albedo detail is the correction the first run of this harness forced,
/// and it is the most important thing in the fixture.** A three-channel shade of
/// a smooth height field has no content anywhere near the pixel rate, so a 3x3
/// neighbourhood in colour space is a *needle* - the local gradient over one
/// pixel - and any sub-pixel reprojection error leaves it. The first run
/// measured a steady-state rejection rate of **95.5%** on exactly that signal,
/// which is the fixture saying it is measuring the shading model rather than the
/// dig. `game_dig`'s terrain is not smooth: `TERRAIN_ALBEDO_ROUGHNESS` is a
/// 512x512 four-layer array tiled at **1.5 world units** per tile, so real
/// content exists in the demo down to `1.5/512 = 0.0029` units. This harness
/// cannot decode a PNG without a new dev-dependency, so it substitutes value
/// noise on the world hit point with the same 1.5-unit base tile and
/// [`DETAIL_OCTAVES`] octaves. Sampled on the *world* point, not in screen
/// space, so it is a texture on the surface and reprojects exactly.
const DETAIL_TILE: f32 = 1.5;
/// Octaves of albedo detail. Finest period `1.5 / 128 = 0.0117` world units,
/// which is under a pixel footprint out to about 3 units and over it beyond
/// that - so the far field is *less* textured than the demo's, which is the
/// conservative direction: an under-textured neighbourhood is a narrower box
/// and rejects more.
const DETAIL_OCTAVES: u32 = 8;
/// Luminance detail amplitude, as a fraction of albedo.
const DETAIL_LUMA: f32 = 0.30;
/// Chromatic detail amplitude. Four octaves, so chroma varies at a coarser rate
/// than luma - the split the paper's two scenes are chosen to span
/// ("high-frequency chromaticity and luminance variation" in Grass,
/// "high-frequency luminance variation" in Asphalt).
const DETAIL_CHROMA: f32 = 0.10;

/// Klosowski et al. 1998's 26-DOP: three coordinate axes, four body diagonals,
/// six edge diagonals. Unnormalised here; [`Dop`] normalises on construction so
/// that the paper's `eps = 1e-5` dilation is in colour units.
const DOP26_AXES: [[f32; 3]; 13] = [
    [1.0, 0.0, 0.0],
    [0.0, 1.0, 0.0],
    [0.0, 0.0, 1.0],
    [1.0, 1.0, 1.0],
    [1.0, 1.0, -1.0],
    [1.0, -1.0, 1.0],
    [1.0, -1.0, -1.0],
    [1.0, 1.0, 0.0],
    [1.0, -1.0, 0.0],
    [1.0, 0.0, 1.0],
    [1.0, 0.0, -1.0],
    [0.0, 1.0, 1.0],
    [0.0, 1.0, -1.0],
];
/// The paper's floating-point dilation of every extent.
const DOP_EPS: f32 = 1e-5;

// ---------------------------------------------------------------------------
// The field.
// ---------------------------------------------------------------------------

/// `game_dig::Ground`, verbatim: distance to a wavy height field, negative
/// below it.
#[derive(Clone, Copy)]
struct Ground;

impl Sdf for Ground {
    type Scalar = f32;

    fn sample(&self, p: [f32; 3]) -> f32 {
        let height = 0.35 * (p[0] * 0.9).sin() * (p[2] * 0.7).cos() + 0.15 * (p[0] * 2.1).sin();
        p[1] - height
    }
}

/// `Ground` with spheres subtracted: `max(f, -(|p - c| - r))`, which is the
/// carve `game_dig`'s `BrushStack` performs.
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
// Small vector helpers. Deliberately not a linear-algebra crate: `isomesh`'s
// public API is `[f32; 3]` and the whole point of the crate is that it needs no
// math library.
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

// ---------------------------------------------------------------------------
// Camera.
// ---------------------------------------------------------------------------

/// A pinhole camera at one instant, with the frame's jitter baked in.
#[derive(Clone, Copy)]
struct Cam {
    eye: [f32; 3],
    right: [f32; 3],
    up: [f32; 3],
    forward: [f32; 3],
    /// `tan(FOV_Y / 2)`.
    tan_half: f32,
    aspect: f32,
    width: f32,
    height: f32,
    /// Sub-pixel jitter, in pixels, applied to this frame's sample positions.
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
    /// behind the near plane. Fractional; the jitter is subtracted, because the
    /// jittered sample at pixel `p` *is* the ray through `p + jitter`.
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

    /// Screen-space radius, in pixels, of the silhouette of a sphere. Exact for
    /// the silhouette cone: its half-angle is `asin(r / d)`.
    fn silhouette_radius(&self, centre: [f32; 3], r: f32) -> f32 {
        let d = dot(sub(centre, self.eye), sub(centre, self.eye)).sqrt();
        if d <= r {
            return f32::INFINITY;
        }
        let half = (r / d).asin();
        half.tan() / self.tan_half * 0.5 * self.height
    }
}

/// R2 low-discrepancy sequence, the 10-entry jitter pattern the paper uses.
fn jitter_for(frame: usize) -> [f32; 2] {
    const A1: f32 = 0.754_877_7;
    const A2: f32 = 0.569_840_3;
    let k = (frame % 10 + 1) as f32;
    [(k * A1).fract() - 0.5, (k * A2).fract() - 0.5]
}

// ---------------------------------------------------------------------------
// Ray casting: game_dig::trace, transcribed.
// ---------------------------------------------------------------------------

/// First surface crossing along a ray inside the sandbox, as a distance.
///
/// A transcription of `game_dig::trace`, including its box test riding along
/// with the surface test rather than gating the march.
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

/// One pixel of the G-buffer.
#[derive(Clone, Copy, Default)]
struct Px {
    hit: bool,
    /// World-space hit point. The motion vector's whole input.
    world: [f32; 3],
    /// Linear signal, RGB.
    rgb: [f32; 3],
}

/// One 32-bit integer hash of a lattice cell, in `[0, 1)`.
///
/// Splitmix64's finaliser truncated to 32 bits: cheap, and it decorrelates the
/// three coordinate multiplies well enough that the noise has no axis-aligned
/// structure. The multipliers are large odd primes.
fn hash_lattice(i: i32, j: i32, k: i32, seed: u32) -> f32 {
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
    (h >> 8) as f32 / 16_777_216.0
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
                acc += wx * wy * wz * hash_lattice(i + dx, j + dy, k + dz, seed);
            }
        }
    }
    acc
}

/// Fractal value noise in `[-1, 1]`, `octaves` deep, lacunarity 2, gain 0.5.
fn fbm(p: [f32; 3], octaves: u32, seed: u32) -> f32 {
    let mut freq = 1.0 / DETAIL_TILE;
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

/// Per-channel albedo multiplier at a world point: luminance detail down to the
/// pixel rate, chroma at a coarser rate. See [`DETAIL_TILE`] for why this exists.
fn detail(p: [f32; 3], textured: bool) -> [f32; 3] {
    if !textured {
        return [1.0, 1.0, 1.0];
    }
    let luma = fbm(p, DETAIL_OCTAVES, 0x5EED_1234);
    let chroma = fbm(p, 4, 0x1234_5EED);
    [
        1.0 + DETAIL_LUMA * luma + DETAIL_CHROMA * chroma,
        1.0 + DETAIL_LUMA * luma,
        1.0 + DETAIL_LUMA * luma - DETAIL_CHROMA * chroma,
    ]
}

/// Shade a hit: three channels that are independent functions of the surface,
/// modulated by the albedo detail, plus distance fog so that a depth change
/// alone moves the signal.
fn shade(n: [f32; 3], ray: [f32; 3], t: f32, world: [f32; 3], textured: bool) -> [f32; 3] {
    let sun = dot(n, SUN_DIR).max(0.0);
    let sky = 0.5 * (n[1] + 1.0);
    let rim = 1.0 - dot(n, mul(ray, -1.0)).max(0.0);
    let fog = 1.0 - (-t * FOG_DENSITY).exp();
    let tex = detail(world, textured);
    let mut out = [0.0f32; 3];
    for c in 0..3 {
        let lit = ALBEDO_SUN[c] * sun + ALBEDO_SKY[c] * sky + ALBEDO_RIM[c] * rim;
        out[c] = lit * tex[c] * (1.0 - fog) + FOG_RGB[c] * fog;
    }
    out
}

/// Ray-cast one frame. Rows are handed to `std::thread`s; the field is pure and
/// shared immutably.
fn render(
    cam: &Cam,
    brushes: &[[f32; 4]],
    w: usize,
    h: usize,
    threads: usize,
    textured: bool,
) -> Vec<Px> {
    let mut buf = vec![Px::default(); w * h];
    let rows_per = h.div_ceil(threads);
    std::thread::scope(|s| {
        for (chunk_index, chunk) in buf.chunks_mut(rows_per * w).enumerate() {
            let y0 = chunk_index * rows_per;
            s.spawn(move || {
                let field = Dug { brushes };
                for (local_y, row) in chunk.chunks_mut(w).enumerate() {
                    let y = y0 + local_y;
                    for (x, px) in row.iter_mut().enumerate() {
                        let ray = cam.ray(x, y);
                        let Some(t) = trace(&field, cam.eye, ray) else {
                            continue;
                        };
                        let world = add(cam.eye, mul(ray, t));
                        let n = normal_at(&field, world);
                        *px = Px {
                            hit: true,
                            world,
                            rgb: shade(n, ray, t, world, textured),
                        };
                    }
                }
            });
        }
    });
    buf
}

// ---------------------------------------------------------------------------
// Colour spaces and bounding volumes.
// ---------------------------------------------------------------------------

/// Karis 2014's YCoCg, the space production TAA builds its AABB in.
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

/// A k-DOP over a set of colours: min/max extents along each of `N` axes,
/// dilated by the paper's `eps`.
///
/// An AABB is the `N = 3` case with the coordinate axes, which is what
/// `Section 3` means by *"AABBs (which is itself a 6-DOP)"*.
struct Dop<const N: usize> {
    axes: [[f32; 3]; N],
    lo: [f32; N],
    hi: [f32; N],
}

impl<const N: usize> Dop<N> {
    fn build(axes: [[f32; 3]; N], colours: &[[f32; 3]]) -> Self {
        let mut lo = [f32::INFINITY; N];
        let mut hi = [f32::NEG_INFINITY; N];
        for c in colours {
            for i in 0..N {
                let d = dot(*c, axes[i]);
                lo[i] = lo[i].min(d);
                hi[i] = hi[i].max(d);
            }
        }
        for i in 0..N {
            lo[i] -= DOP_EPS;
            hi[i] += DOP_EPS;
        }
        Self { axes, lo, hi }
    }

    /// Clip the ray `centre -> history` to the shell. Returns the parameter in
    /// `[0, 1]`: `1.0` means the history was inside and nothing was rejected.
    ///
    /// The ray starts inside (the centre pixel's own colour is one of the
    /// neighbourhood colours), and a DOP is convex, so this reduces to the
    /// nearest slab crossing - the simplification Section 3 names.
    fn clip(&self, centre: [f32; 3], history: [f32; 3]) -> f32 {
        let d = sub(history, centre);
        let mut s = 1.0f32;
        for i in 0..N {
            let dd = dot(d, self.axes[i]);
            if dd.abs() < 1e-20 {
                continue;
            }
            let c = dot(centre, self.axes[i]);
            let bound = if dd > 0.0 { self.hi[i] } else { self.lo[i] };
            let si = (bound - c) / dd;
            if si < s {
                s = si;
            }
        }
        s.clamp(0.0, 1.0)
    }
}

/// Normalise a DOP axis set once, so `DOP_EPS` is in colour units on every axis.
fn normalised_axes<const N: usize>(axes: [[f32; 3]; N]) -> [[f32; 3]; N] {
    let mut out = axes;
    for a in &mut out {
        *a = norm(*a);
    }
    out
}

/// The identity axes, i.e. an AABB in whatever space it is handed.
const AABB_AXES: [[f32; 3]; 3] = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];

// ---------------------------------------------------------------------------
// One frame's counters.
// ---------------------------------------------------------------------------

#[derive(Clone, Copy, Default)]
struct Frame {
    /// Pixels with a current hit and a valid on-screen reprojection with a hit.
    population: u64,
    /// Pixels with a current hit and no usable history: off-screen, behind the
    /// near plane, or reprojecting onto a miss. Not counted as clamp
    /// rejections; reported so the reader can see how big the excluded set is.
    no_history: u64,
    rejected_aabb_rgb: u64,
    rejected_aabb_ycocg: u64,
    rejected_dop26_rgb: u64,
    rejected_dop26_ycocg: u64,
    /// Rejected by the baseline `aabb_ycocg`, accepted by `dop26_ycocg`. Zero by
    /// the containment property; measured anyway.
    recovered_dop26_ycocg: u64,
    /// Rejected by the baseline, accepted by `dop26_rgb`. The control that the
    /// recovery counter can report a non-zero.
    recovered_dop26_rgb: u64,
    /// Accepted by the baseline, rejected by `dop26_ycocg`. What the paper
    /// actually claims to do.
    extra_dop26_ycocg: u64,
    /// Baseline rejections where the clip discarded at least **half** the
    /// history's deviation from the current colour, i.e. `s < 0.5`. The
    /// registered "rejected" is `s < 1.0` - the clamp moved the sample at all -
    /// and this is the same count under a magnitude threshold, so a reader can
    /// see whether the headline is knife-edge.
    rejected_half: u64,
    /// Sum of the clip parameter `s` over the rejected samples.
    ///
    /// **This is the column that separates the two explanations C2's falsifier
    /// names.** `s` is where on the segment from the current colour to the
    /// history colour the neighbourhood's shell lies, so `s` near 1 means the
    /// history was barely outside - clamp conservatism, which a tighter or
    /// better-shaped volume could plausibly arbitrate - and `s` near 0 means it
    /// was nowhere near the neighbourhood, which is a genuine disocclusion and
    /// no clipping scheme can help.
    clip_s: f64,
    /// Sum of reprojection displacements, in pixels, over the population. The
    /// input `R-076` and `R-091` need: a clamp cannot keep history that has
    /// moved further than the neighbourhood it is compared against.
    reproj_px: f64,
    /// Baseline rejections inside the brush box dilated by one dig radius.
    rejected_in_box: u64,
    /// Pixels whose signal differs from the undug world's.
    changed: u64,
    /// Changed pixels outside every brush silhouette. Must be zero.
    changed_outside: u64,
    /// Single-threaded wall-clock of the baseline clamp pass, ms.
    aabb_ms: f64,
    /// Single-threaded wall-clock of the 26-DOP clamp pass, ms.
    dop_ms: f64,
}

impl Frame {
    fn fraction(&self) -> f64 {
        if self.population == 0 {
            0.0
        } else {
            self.rejected_aabb_ycocg as f64 / self.population as f64
        }
    }
}

/// One TAA history buffer, plus the resolve.
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

/// Bilinear fetch from a buffer, requiring all four taps valid.
fn fetch(hist: &[[f32; 3]], valid: &[bool], w: usize, h: usize, p: [f32; 2]) -> Option<[f32; 3]> {
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

/// The 3x3 neighbourhood of the current frame, clamped at the border.
fn neighbourhood(cur: &[Px], w: usize, h: usize, x: usize, y: usize) -> ([[f32; 3]; 9], usize) {
    let mut out = [[0.0f32; 3]; 9];
    let mut n = 0;
    for dy in -1i32..=1 {
        for dx in -1i32..=1 {
            let sx = (x as i32 + dx).clamp(0, w as i32 - 1) as usize;
            let sy = (y as i32 + dy).clamp(0, h as i32 - 1) as usize;
            let p = &cur[sy * w + sx];
            if p.hit {
                out[n] = p.rgb;
                n += 1;
            }
        }
    }
    (out, n)
}

// ---------------------------------------------------------------------------
// One arm.
// ---------------------------------------------------------------------------

struct ArmCfg {
    name: &'static str,
    w: usize,
    h: usize,
    /// Camera pitch, radians. `game_dig`'s opening `Look::pitch` is `-0.15`,
    /// which aims 11.4 units away; `Look::pitch` is user-driven and clamped to
    /// +-1.5, so `-0.6` is the same demo with the player looking down at the
    /// rock under their feet, 3.0 units away. Both are `game_dig` poses and the
    /// difference is the brush's screen footprint, which is what C1's ceiling is
    /// made of.
    pitch: f32,
    /// Walk speed, units/s. `game_dig`'s unmodified walk is 2.5.
    walk: f32,
    /// Yaw rate, rad/s.
    yaw_rate: f32,
    /// Albedo detail on. Off is the fixture-defect arm: see [`DETAIL_TILE`].
    textured: bool,
    /// TAA jitter on. **Off is the registered vacuity control's proof.** With a
    /// frozen camera and no jitter the reprojection is the identity, the history
    /// converges to the current frame, and the rejection count is *zero* - so
    /// the non-zero steady rate every other arm reports is a measurement rather
    /// than a fixture that could not have produced a zero (`M-44`).
    jitter: bool,
    /// `false` for a single stroke, `true` for a held button at
    /// [`STROKE_RATE_HZ`].
    continuous: bool,
}

struct ArmOut {
    control: Vec<Frame>,
    treatment: Vec<Frame>,
    brushes_placed: usize,
    box_area_fraction: f64,
}

/// Simulate the control (never dug) and treatment (dug) worlds in lockstep over
/// one camera path and one jitter sequence.
fn run_arm(cfg: &ArmCfg, threads: usize) -> ArmOut {
    let (w, h) = (cfg.w, cfg.h);
    let n = w * h;
    let dop26 = normalised_axes(DOP26_AXES);
    let aabb = AABB_AXES;

    let mut taa_c = Taa::new(n);
    let mut taa_t = Taa::new(n);
    let mut prev_cam: Option<Cam> = None;
    let mut brushes: Vec<[f32; 4]> = Vec::new();
    let mut since_edit = f32::INFINITY;

    let mut control = Vec::new();
    let mut treatment = Vec::new();
    let mut box_area_total = 0.0f64;
    let mut box_area_frames = 0usize;

    let total = WARMUP + WINDOW;
    let mut eye = EYE_START;
    let mut yaw = 0.0f32;

    for frame in 0..total {
        let jitter = if cfg.jitter {
            jitter_for(frame)
        } else {
            [0.0, 0.0]
        };
        let cam = Cam::new(eye, yaw, cfg.pitch, w, h, jitter);

        // Edits land at the start of the measurement window, through the same
        // path a click takes: aim down the camera's forward ray, place a brush
        // of `BRUSH_RADIUS` at the hit.
        if frame >= WARMUP {
            let place = if cfg.continuous {
                if since_edit.is_infinite() {
                    since_edit = 0.0;
                    true
                } else {
                    since_edit += DT;
                    if since_edit >= EDIT_PERIOD {
                        since_edit -= EDIT_PERIOD;
                        true
                    } else {
                        false
                    }
                }
            } else {
                frame == WARMUP
            };
            if place {
                let field = Dug { brushes: &brushes };
                if let Some(t) = trace(&field, cam.eye, cam.forward) {
                    let c = add(cam.eye, mul(cam.forward, t));
                    brushes.push([c[0], c[1], c[2], BRUSH_RADIUS]);
                }
            }
        }

        let cur_c = render(&cam, &[], w, h, threads, cfg.textured);
        let cur_t = render(&cam, &brushes, w, h, threads, cfg.textured);

        // Brush silhouettes in this frame's pixel grid: centre and radius.
        let discs: Vec<[f32; 3]> = brushes
            .iter()
            .filter_map(|b| {
                let c = [b[0], b[1], b[2]];
                cam.project(c)
                    .map(|s| [s[0], s[1], cam.silhouette_radius(c, b[3])])
            })
            .collect();
        if !discs.is_empty() {
            // Summed area of the dilated boxes, over-counted where they overlap
            // and clipped to the screen. This is only the `b` in the module
            // docs' `(R - 1 + b)/R`, i.e. the share of the frame C3's container
            // occupies; over-counting it is the conservative direction.
            let area: f64 = discs
                .iter()
                .map(|d| {
                    let sx = (4.0 * d[2]).min(cam.width) as f64;
                    let sy = (4.0 * d[2]).min(cam.height) as f64;
                    sx * sy
                })
                .sum();
            box_area_total += area / (w as f64 * h as f64);
            box_area_frames += 1;
        }

        // The four extra volumes and the cross-tabulation are only needed in the
        // measurement window; a warm-up frame pays for the baseline resolve
        // alone, because that is the only pass that feeds the history buffer.
        let measure = frame >= WARMUP;
        let fc = resolve(
            &mut taa_c,
            &cur_c,
            &cur_c,
            prev_cam.as_ref(),
            w,
            h,
            &aabb,
            &dop26,
            // The control gets the same boxes, so C3's container can be
            // measured on the undug world too. Without that, the edit-
            // attributable concentration would have to be *modelled* from an
            // area fraction instead of counted.
            &discs,
            measure,
        );
        let ft = resolve(
            &mut taa_t,
            &cur_t,
            &cur_c,
            prev_cam.as_ref(),
            w,
            h,
            &aabb,
            &dop26,
            &discs,
            measure,
        );

        if frame >= WARMUP {
            control.push(fc);
            treatment.push(ft);
        }

        prev_cam = Some(cam);
        eye = add(eye, mul(Cam::walk_dir(yaw), cfg.walk * DT));
        yaw += cfg.yaw_rate * DT;
    }

    ArmOut {
        control,
        treatment,
        brushes_placed: brushes.len(),
        box_area_fraction: if box_area_frames == 0 {
            0.0
        } else {
            box_area_total / box_area_frames as f64
        },
    }
}

/// The TAA resolve for one frame, and the frame's counters.
///
/// `reference` is the undug world's G-buffer, used only to count changed pixels;
/// for the control arm it is the same buffer and the count is zero by
/// construction.
#[allow(clippy::too_many_arguments)]
fn resolve(
    taa: &mut Taa,
    cur: &[Px],
    reference: &[Px],
    prev_cam: Option<&Cam>,
    w: usize,
    h: usize,
    aabb_axes: &[[f32; 3]; 3],
    dop_axes: &[[f32; 3]; 13],
    discs: &[[f32; 3]],
    measure: bool,
) -> Frame {
    let mut f = Frame::default();
    let n = w * h;

    // Which pixels have usable history, and where from. One pass, so the four
    // clamp passes below all see the same population.
    let mut hist_rgb: Vec<Option<[f32; 3]>> = vec![None; n];
    if let Some(pc) = prev_cam {
        for y in 0..h {
            for x in 0..w {
                let i = y * w + x;
                if !cur[i].hit {
                    continue;
                }
                let Some(q) = pc.project(cur[i].world) else {
                    f.no_history += 1;
                    continue;
                };
                match fetch(&taa.hist, &taa.valid, w, h, q) {
                    Some(hrgb) => {
                        hist_rgb[i] = Some(hrgb);
                        f.population += 1;
                        let (dx, dy) = (q[0] - x as f32, q[1] - y as f32);
                        f.reproj_px += f64::from((dx * dx + dy * dy).sqrt());
                    }
                    None => f.no_history += 1,
                }
            }
        }
    }

    // Changed-pixel accounting and the silhouette containment control.
    for i in 0..n {
        let (a, b) = (&cur[i], &reference[i]);
        let differs = a.hit != b.hit
            || (0..3).any(|c| (a.rgb[c] - b.rgb[c]).abs() > CHANGE_TOL)
            || (0..3).any(|c| (a.world[c] - b.world[c]).abs() > CHANGE_TOL);
        if !differs {
            continue;
        }
        f.changed += 1;
        let (x, y) = ((i % w) as f32, (i / w) as f32);
        let inside = discs.iter().any(|d| {
            let dx = x - d[0];
            let dy = y - d[1];
            (dx * dx + dy * dy).sqrt() <= d[2] + SILHOUETTE_SLACK
        });
        if !inside {
            f.changed_outside += 1;
        }
    }

    // Four clamp passes over the same population. Separate passes because that
    // is what makes `aabb_ms` and `dop_ms` comparable: each one rebuilds the 3x3
    // neighbourhood exactly as a resolve shader does, and neither pays for the
    // other's work. Single-threaded on purpose: a per-pixel cost divided by an
    // unstated thread count is not a cost.
    let mut clipped: Vec<[f32; 3]> = vec![[0.0; 3]; n];

    // Pass 1: the baseline, Karis 2014's YCoCg AABB. Timed. Drives the history.
    let t0 = Instant::now();
    for y in 0..h {
        for x in 0..w {
            let i = y * w + x;
            let Some(hrgb) = hist_rgb[i] else { continue };
            let (nb, k) = neighbourhood(cur, w, h, x, y);
            let mut cols = [[0.0f32; 3]; 9];
            for j in 0..k {
                cols[j] = to_ycocg(nb[j]);
            }
            let vol = Dop::build(*aabb_axes, &cols[..k]);
            let c = to_ycocg(cur[i].rgb);
            let s = vol.clip(c, to_ycocg(hrgb));
            clipped[i] = from_ycocg(add(c, mul(sub(to_ycocg(hrgb), c), s)));
            if s < 1.0 {
                f.rejected_aabb_ycocg += 1;
                f.clip_s += f64::from(s);
            }
            if s < 0.5 {
                f.rejected_half += 1;
            }
        }
    }
    f.aabb_ms = t0.elapsed().as_secs_f64() * 1e3;

    if !measure {
        commit(taa, cur, &hist_rgb, &clipped);
        return f;
    }

    // Pass 2: the paper's method, in the baseline's colour space. Timed, and
    // the only other timed pass, so `kdop_ms` and `aabb_ms` differ in exactly
    // one thing: the volume.
    let t1 = Instant::now();
    for y in 0..h {
        for x in 0..w {
            let i = y * w + x;
            let Some(hrgb) = hist_rgb[i] else { continue };
            let (nb, k) = neighbourhood(cur, w, h, x, y);
            let mut cols = [[0.0f32; 3]; 9];
            for j in 0..k {
                cols[j] = to_ycocg(nb[j]);
            }
            let vol = Dop::build(*dop_axes, &cols[..k]);
            let s = vol.clip(to_ycocg(cur[i].rgb), to_ycocg(hrgb));
            if s < 1.0 {
                f.rejected_dop26_ycocg += 1;
            }
        }
    }
    f.dop_ms = t1.elapsed().as_secs_f64() * 1e3;

    // Pass 3, untimed: the RGB AABB (Lottes 2011), a 26-DOP in RGB, and the
    // cross-tabulation C2 asks for -- all four verdicts on the *same* fetched
    // history sample, so the comparison is per-sample rather than between two
    // pipelines that have drifted apart. The RGB 26-DOP is the control on the
    // recovery counter: its axes are not the baseline's, so it is not contained
    // in the baseline volume and *can* accept what the baseline rejects.
    for y in 0..h {
        for x in 0..w {
            let i = y * w + x;
            let Some(hrgb) = hist_rgb[i] else { continue };
            let (nb, k) = neighbourhood(cur, w, h, x, y);
            let mut cols = [[0.0f32; 3]; 9];
            for j in 0..k {
                cols[j] = to_ycocg(nb[j]);
            }
            let c = cur[i].rgb;
            let cy = to_ycocg(c);
            let hy = to_ycocg(hrgb);
            let base_rejected = Dop::build(*aabb_axes, &cols[..k]).clip(cy, hy) < 1.0;
            let dop_y_rejected = Dop::build(*dop_axes, &cols[..k]).clip(cy, hy) < 1.0;
            let dop_rgb_rejected = Dop::build(*dop_axes, &nb[..k]).clip(c, hrgb) < 1.0;
            if Dop::build(*aabb_axes, &nb[..k]).clip(c, hrgb) < 1.0 {
                f.rejected_aabb_rgb += 1;
            }
            if dop_rgb_rejected {
                f.rejected_dop26_rgb += 1;
            }
            if base_rejected {
                if !dop_y_rejected {
                    f.recovered_dop26_ycocg += 1;
                }
                if !dop_rgb_rejected {
                    f.recovered_dop26_rgb += 1;
                }
                // C3's container: the brush's screen-space bounding box dilated
                // by one dig radius, i.e. a box of half-side `2 * r_screen`.
                let (px, py) = (x as f32, y as f32);
                let inside = discs
                    .iter()
                    .any(|d| (px - d[0]).abs() <= 2.0 * d[2] && (py - d[1]).abs() <= 2.0 * d[2]);
                if inside {
                    f.rejected_in_box += 1;
                }
            } else if dop_y_rejected {
                f.extra_dop26_ycocg += 1;
            }
        }
    }

    commit(taa, cur, &hist_rgb, &clipped);
    f
}

/// Commit the history: blend the clipped value, and seed pixels with no usable
/// history from the current frame.
fn commit(taa: &mut Taa, cur: &[Px], hist_rgb: &[Option<[f32; 3]>], clipped: &[[f32; 3]]) {
    for i in 0..cur.len() {
        if !cur[i].hit {
            taa.valid[i] = false;
            continue;
        }
        taa.valid[i] = true;
        taa.hist[i] = match hist_rgb[i] {
            Some(_) => add(mul(clipped[i], 1.0 - ALPHA), mul(cur[i].rgb, ALPHA)),
            None => cur[i].rgb,
        };
    }
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

/// The numbers printed at the end, from [`HEADLINE_ARM`].
struct Headline {
    steady: f64,
    post: f64,
    ratio: f64,
    frames_elevated: usize,
    kdop_ms: f64,
    recovered: f64,
    concentration: f64,
    concentration_attributable: f64,
    ceiling: f64,
    reproj_px: f64,
}

fn main() {
    if !std::env::args().any(|a| a == "--bench") {
        return;
    }
    let prereg = isomesh::experiment!("P-77");
    let threads = std::thread::available_parallelism().map_or(8, |n| n.get());
    let mhz = cpu_mhz();

    // Ten arms. The registered columns come from [`HEADLINE_ARM`]; the rest
    // exist because C1's ceiling turned out to be made of two things the
    // registration does not name - the brush's screen footprint and the
    // reprojection displacement - and each arm moves exactly one of them.
    let arms = [
        ArmCfg {
            name: "dig_at_feet_walk",
            w: 960,
            h: 540,
            pitch: -0.6,
            walk: WALK_SPEED,
            yaw_rate: YAW_RATE,
            textured: true,
            jitter: true,
            continuous: false,
        },
        ArmCfg {
            name: "dig_at_feet_static",
            w: 960,
            h: 540,
            pitch: -0.6,
            walk: 0.0,
            yaw_rate: 0.0,
            textured: true,
            jitter: true,
            continuous: false,
        },
        ArmCfg {
            name: "dig_at_feet_walk_held",
            w: 960,
            h: 540,
            pitch: -0.6,
            walk: WALK_SPEED,
            yaw_rate: YAW_RATE,
            textured: true,
            jitter: true,
            continuous: true,
        },
        ArmCfg {
            name: "dig_at_feet_walk_untextured",
            w: 960,
            h: 540,
            pitch: -0.6,
            walk: WALK_SPEED,
            yaw_rate: YAW_RATE,
            textured: false,
            jitter: true,
            continuous: false,
        },
        ArmCfg {
            name: "dig_at_feet_walk_half_res",
            w: 480,
            h: 270,
            pitch: -0.6,
            walk: WALK_SPEED,
            yaw_rate: YAW_RATE,
            textured: true,
            jitter: true,
            continuous: false,
        },
        ArmCfg {
            name: "demo_pose_walk",
            w: 960,
            h: 540,
            pitch: PITCH,
            walk: WALK_SPEED,
            yaw_rate: YAW_RATE,
            textured: true,
            jitter: true,
            continuous: false,
        },
        ArmCfg {
            name: "demo_pose_static",
            w: 960,
            h: 540,
            pitch: PITCH,
            walk: 0.0,
            yaw_rate: 0.0,
            textured: true,
            jitter: true,
            continuous: false,
        },
        // Two intermediate points on the locomotion sweep. The first run showed
        // that the steady-state rejection rate is a function of the
        // *reprojection displacement* and nothing else, so the honest form of
        // C1 is a curve in that variable rather than one number: 0.61 px
        // (jitter only), 1.9 px, 3.6 px, 7.2 px.
        ArmCfg {
            name: "dig_at_feet_walk_quarter",
            w: 960,
            h: 540,
            pitch: -0.6,
            walk: WALK_SPEED * 0.25,
            yaw_rate: YAW_RATE * 0.25,
            textured: true,
            jitter: true,
            continuous: false,
        },
        ArmCfg {
            name: "dig_at_feet_walk_half_speed",
            w: 960,
            h: 540,
            pitch: -0.6,
            walk: WALK_SPEED * 0.5,
            yaw_rate: YAW_RATE * 0.5,
            textured: true,
            jitter: true,
            continuous: false,
        },
        // The vacuity control's other half: the configuration in which the
        // steady-state rejection rate really is zero. Frozen camera, no jitter,
        // so the reprojection is the identity and the history converges to the
        // current frame. See `ArmCfg::jitter`.
        ArmCfg {
            name: "zero_proof_no_jitter_static",
            w: 960,
            h: 540,
            pitch: -0.6,
            walk: 0.0,
            yaw_rate: 0.0,
            textured: true,
            jitter: false,
            continuous: false,
        },
    ];

    let mut rows: Vec<Row> = Vec::new();
    let mut headline: Option<Headline> = None;
    // Steady-state rejection count of the `zero_proof_no_jitter_static` arm.
    let mut zero_proof: Option<u64> = None;
    let mut headline_steady_rejected: u64 = 0;

    for cfg in &arms {
        let started = Instant::now();
        let out = run_arm(cfg, threads);
        eprintln!(
            "P-77: arm {} ({}x{}) in {:.1}s, {} brushes",
            cfg.name,
            cfg.w,
            cfg.h,
            started.elapsed().as_secs_f64(),
            out.brushes_placed
        );

        // --- controls, before any number is reported ---------------------
        for (i, (c, t)) in out.control.iter().zip(&out.treatment).enumerate() {
            assert!(
                c.population > 0,
                "P-77 {}: frame +{i} of the control arm has no history samples at all -- \
                 an empty depth buffer reports a zero that could not have been non-zero",
                cfg.name
            );
            assert_eq!(
                c.changed, 0,
                "P-77 {}: the control world was dug",
                cfg.name
            );
            assert_eq!(
                t.changed_outside, 0,
                "P-77 {}: frame +{i} has {} changed pixels outside every brush silhouette, \
                 so either the projection or the brush placement is wrong -- subtracting a \
                 sphere cannot move the zero set outside the sphere",
                cfg.name, t.changed_outside
            );
        }
        let steady_rejected: u64 = out.control.iter().map(|f| f.rejected_aabb_ycocg).sum();
        if cfg.jitter {
            // THE REGISTERED VACUITY CONTROL.
            assert!(
                steady_rejected > 0,
                "P-77 {}: VACUITY CONTROL FAILED -- the steady-state arm rejected zero \
                 history samples over {WINDOW} frames, so the registered 5x ratio would be \
                 division by a floor",
                cfg.name
            );
        } else {
            // The other half of the same control: this arm exists to show that a
            // zero was reachable. If it is not zero, the vacuity control above is
            // not a control, because nothing in the fixture could have produced
            // the zero it rules out.
            zero_proof = Some(steady_rejected);
        }

        // --- C1 ----------------------------------------------------------
        let steady_pop: u64 = out.control.iter().map(|f| f.population).sum();
        let steady = steady_rejected as f64 / steady_pop as f64;
        let post = out.treatment[0].fraction();
        // A ratio over a zero floor is not a number. The zero-proof arm reports
        // `NaN` here rather than an infinity that a reader might quote.
        let evaluable = steady > 0.0;
        let ratio = if evaluable { post / steady } else { f64::NAN };
        let frames_elevated = if evaluable {
            out.treatment
                .iter()
                .take_while(|f| f.fraction() >= ELEVATION_FACTOR * steady)
                .count()
        } else {
            0
        };

        // --- C2 ----------------------------------------------------------
        let base_rej = out.treatment[0].rejected_aabb_ycocg;
        let recovered = out.treatment[0].recovered_dop26_ycocg as f64 / base_rej.max(1) as f64;
        let recovered_rgb_axes =
            out.treatment[0].recovered_dop26_rgb as f64 / base_rej.max(1) as f64;
        let kdop_ms = out.treatment[0].dop_ms;
        let aabb_ms = out.treatment[0].aabb_ms;

        // --- C3 ----------------------------------------------------------
        let in_box = out.treatment[0].rejected_in_box;
        let concentration = in_box as f64 / base_rej.max(1) as f64;
        // The same fraction over the *edit-attributable* rejections only:
        // treatment minus the frame-paired control, in the box and in total.
        // C3's registered denominator is every rejection in the frame, most of
        // which is the locomotion floor and has nothing to do with the dig, so
        // this is the number that tests the clause's actual mechanism.
        let concentration_attributable = {
            let d_box = out.treatment[0].rejected_in_box as f64
                - out.control[0].rejected_in_box as f64;
            let d_tot = out.treatment[0].rejected_aabb_ycocg as f64
                - out.control[0].rejected_aabb_ycocg as f64;
            if d_tot <= 0.0 { f64::NAN } else { d_box / d_tot }
        };

        // C1's arithmetic ceiling, from the fixture rather than from the clause.
        // The dig can only add rejections at pixels whose signal it changed, and
        // that set is provably inside the brush's silhouette, so the ratio can
        // never exceed `1 + changed / steady_rejected`. Computed on the edit
        // frame, where the changed set is largest relative to the history.
        let ceiling = {
            let c0 = out.control[0].rejected_aabb_ycocg.max(1) as f64;
            (c0 + out.treatment[0].changed as f64) / c0
        };

        let c1 = evaluable && ratio >= 5.0 && frames_elevated >= 3;
        let c2 = recovered >= 0.5 && kdop_ms < 0.3;
        let c3 = concentration > 0.80;

        if cfg.name == HEADLINE_ARM {
            headline_steady_rejected = steady_rejected;
            headline = Some(Headline {
                steady,
                post,
                ratio,
                frames_elevated,
                kdop_ms,
                recovered,
                concentration,
                concentration_attributable,
                ceiling,
                reproj_px: out.control[0].reproj_px
                    / out.control[0].population.max(1) as f64,
            });
        }

        let summary = |arm: String, f: &Frame, offset: Option<usize>, ctrl: Option<&Frame>| -> Row {
            let mut r: Row = vec![
                ("arm", arm),
                ("stroke_rate_hz", format!("{STROKE_RATE_HZ:.1}")),
                ("steady_rejection_fraction", format!("{steady:.8}")),
                ("post_edit_rejection_fraction", format!("{post:.8}")),
                ("ratio", format!("{ratio:.4}")),
                ("frames_elevated", frames_elevated.to_string()),
                ("kdop_ms", format!("{kdop_ms:.4}")),
                ("kdop_recovered_fraction", format!("{recovered:.8}")),
                ("rejected_in_brush_box", in_box.to_string()),
                ("rejected_total", base_rej.to_string()),
                ("concentration_fraction", format!("{concentration:.6}")),
                ("c1_holds", c1.to_string()),
                ("c2_holds", c2.to_string()),
                ("c3_holds", c3.to_string()),
                // Extras.
                ("c1_evaluable", evaluable.to_string()),
                ("taa_jitter", cfg.jitter.to_string()),
                ("fixture", cfg.name.to_string()),
                ("width", cfg.w.to_string()),
                ("height", cfg.h.to_string()),
                ("camera_walk_units_per_s", format!("{:.2}", cfg.walk)),
                ("camera_yaw_rad_per_s", format!("{:.3}", cfg.yaw_rate)),
                ("camera_pitch_rad", format!("{:.3}", cfg.pitch)),
                ("albedo_detail", cfg.textured.to_string()),
                ("held_button", cfg.continuous.to_string()),
                ("ratio_ceiling_from_footprint", format!("{ceiling:.4}")),
                (
                    "concentration_edit_attributable",
                    if concentration_attributable.is_nan() {
                        "NA".to_string()
                    } else {
                        format!("{concentration_attributable:.6}")
                    },
                ),
                (
                    "kdop_ms_per_megapixel",
                    format!("{:.4}", kdop_ms * 1e6 / (cfg.w * cfg.h) as f64),
                ),
                ("brushes_placed", out.brushes_placed.to_string()),
                ("box_area_fraction", format!("{:.8}", out.box_area_fraction)),
                ("threads", threads.to_string()),
                ("cpu_mhz", mhz.clone()),
                ("aabb_ms", format!("{aabb_ms:.4}")),
                (
                    "clamp_ratio_dop_over_aabb",
                    format!("{:.4}", kdop_ms / aabb_ms),
                ),
                (
                    "kdop_recovered_fraction_rgb_axes",
                    format!("{recovered_rgb_axes:.8}"),
                ),
                ("frame_offset", offset.map_or("NA".to_string(), |o| o.to_string())),
                ("population", f.population.to_string()),
                ("no_history", f.no_history.to_string()),
                ("frame_rejected", f.rejected_aabb_ycocg.to_string()),
                ("frame_rejection_fraction", format!("{:.8}", f.fraction())),
                ("frame_rejected_aabb_rgb", f.rejected_aabb_rgb.to_string()),
                ("frame_rejected_dop26_ycocg", f.rejected_dop26_ycocg.to_string()),
                ("frame_rejected_dop26_rgb", f.rejected_dop26_rgb.to_string()),
                ("frame_recovered_dop26_ycocg", f.recovered_dop26_ycocg.to_string()),
                ("frame_recovered_dop26_rgb", f.recovered_dop26_rgb.to_string()),
                ("frame_extra_rejected_dop26_ycocg", f.extra_dop26_ycocg.to_string()),
                ("frame_changed_pixels", f.changed.to_string()),
                ("frame_changed_outside_silhouette", f.changed_outside.to_string()),
                ("frame_rejected_in_box", f.rejected_in_box.to_string()),
                ("frame_rejected_half", f.rejected_half.to_string()),
                (
                    "frame_mean_clip_s_over_rejected",
                    format!(
                        "{:.6}",
                        if f.rejected_aabb_ycocg == 0 {
                            f64::NAN
                        } else {
                            f.clip_s / f.rejected_aabb_ycocg as f64
                        }
                    ),
                ),
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
                ("frame_aabb_ms", format!("{:.4}", f.aabb_ms)),
                ("frame_dop_ms", format!("{:.4}", f.dop_ms)),
            ];
            let (cf, cr) = ctrl.map_or((f64::NAN, 0u64), |c| (c.fraction(), c.rejected_aabb_ycocg));
            r.push(("control_rejection_fraction", if cf.is_nan() { "NA".to_string() } else { format!("{cf:.8}") }));
            r.push(("control_rejected", if ctrl.is_some() { cr.to_string() } else { "NA".to_string() }));
            r.push((
                "control_rejected_in_box",
                ctrl.map_or("NA".to_string(), |c| c.rejected_in_box.to_string()),
            ));
            r.push((
                "paired_frame_ratio",
                if cf.is_nan() || cf == 0.0 {
                    "NA".to_string()
                } else {
                    format!("{:.4}", f.fraction() / cf)
                },
            ));
            // Per-frame edit-attributable concentration: the same
            // treatment-minus-control difference as the summary column, but on
            // this frame alone, so the persistence of the *locality* can be read
            // off the series too.
            r.push((
                "frame_concentration_edit_attributable",
                ctrl.map_or("NA".to_string(), |c| {
                    let d_box = f.rejected_in_box as f64 - c.rejected_in_box as f64;
                    let d_tot =
                        f.rejected_aabb_ycocg as f64 - c.rejected_aabb_ycocg as f64;
                    if d_tot <= 0.0 {
                        "NA".to_string()
                    } else {
                        format!("{:.6}", d_box / d_tot)
                    }
                }),
            ));
            r
        };

        // Window-aggregate rows, then the per-frame series for both worlds.
        let agg = |frames: &[Frame]| -> Frame {
            let mut a = Frame::default();
            for f in frames {
                a.population += f.population;
                a.no_history += f.no_history;
                a.rejected_aabb_rgb += f.rejected_aabb_rgb;
                a.rejected_aabb_ycocg += f.rejected_aabb_ycocg;
                a.rejected_dop26_rgb += f.rejected_dop26_rgb;
                a.rejected_dop26_ycocg += f.rejected_dop26_ycocg;
                a.recovered_dop26_ycocg += f.recovered_dop26_ycocg;
                a.recovered_dop26_rgb += f.recovered_dop26_rgb;
                a.extra_dop26_ycocg += f.extra_dop26_ycocg;
                a.rejected_in_box += f.rejected_in_box;
                a.changed += f.changed;
                a.changed_outside += f.changed_outside;
                a.aabb_ms += f.aabb_ms;
                a.rejected_half += f.rejected_half;
                a.clip_s += f.clip_s;
                a.reproj_px += f.reproj_px;
            }
            a
        };

        rows.push(summary(
            format!("{}/steady_window", cfg.name),
            &agg(&out.control),
            None,
            None,
        ));
        rows.push(summary(
            format!("{}/post_edit_window", cfg.name),
            &agg(&out.treatment),
            None,
            Some(&agg(&out.control)),
        ));
        for (i, (c, t)) in out.control.iter().zip(&out.treatment).enumerate() {
            rows.push(summary(
                format!("{}/frame_plus_{i}", cfg.name),
                t,
                Some(i),
                Some(c),
            ));
        }
    }

    // The vacuity control is only a control if the zero it rules out was
    // reachable. `zero_proof_no_jitter_static` is the same fixture with the
    // jitter off: identity reprojection, history converging to the current
    // frame. Its rejection count must be at least two orders of magnitude below
    // the headline arm's, or the non-zero steady rate is a property of the
    // harness rather than of the resolve.
    let zp = zero_proof.expect("the zero-proof arm always runs");
    assert!(
        (zp as f64) * 100.0 < headline_steady_rejected as f64,
        "P-77: the vacuity control is not a control -- with jitter off and the camera \
         frozen the resolve still rejected {zp} history samples against the headline arm's \
         {headline_steady_rejected}, so nothing in this fixture could have produced the zero \
         the registered control rules out"
    );
    println!(
        "P-77 vacuity control: steady rejections {headline_steady_rejected} with jitter, \
         {zp} without -- the zero was reachable."
    );

    let Headline {
        steady,
        post,
        ratio,
        frames_elevated,
        kdop_ms,
        recovered,
        concentration,
        concentration_attributable,
        ceiling,
        reproj_px,
    } = headline.expect("the headline arm always runs");
    println!();
    println!("P-77 headline, from the `{HEADLINE_ARM}` arm at 960x540:");
    println!("  mean reprojection            {reproj_px:.3} px/frame");
    println!("  steady_rejection_fraction    {steady:.6}  (VACUITY CONTROL: must be > 0)");
    println!("  post_edit_rejection_fraction {post:.6}");
    println!("  ratio                        {ratio:.3}  (C1 needs >= 5)");
    println!("  ratio ceiling from footprint {ceiling:.3}  (arithmetic, from the changed set)");
    println!("  frames_elevated              {frames_elevated}       (C1 needs >= 3)");
    println!("  kdop_recovered_fraction      {recovered:.6}  (C2 needs >= 0.5)");
    println!("  kdop_ms (CPU clamp pass)     {kdop_ms:.4}  (not a GPU resolve cost)");
    println!("  concentration_fraction       {concentration:.4}  (C3 needs > 0.80)");
    println!("  edit-attributable            {concentration_attributable:.4}  (C3's mechanism, paired denominator)");
    println!();

    common::experiment::run(prereg, |run| {
        for row in &rows {
            run.record(row);
        }
    });
}
