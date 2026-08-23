//! E-305 — Lipschitz tape pruning: the frame cost stops tracking the whole
//! edit history.
//!
//! ```bash
//! cd bevy_isomesh && cargo run --example tape_pruning --release
//! ```
//!
//! **Always `--release`.** `P` turns pruning off, `H` the survivor-fraction
//! overlay, `C` the chunk boxes, `X` restarts the stroke.
//!
//! The stroke paces itself off the capture length, so any window works and none
//! of them records a still. The recipe that shows the whole staircase build and
//! then lets the numbers settle on it, measured at **1.52 MB**:
//!
//! ```bash
//! ISOMESH_CAPTURE_FRAMES=60 ISOMESH_CAPTURE_EVERY=8 ISOMESH_WINDOW=1280x720 \
//!   FPS=12 ./scripts/record_gif.sh tape_pruning out.gif
//! ```
//!
//! **Adding `ISOMESH_SPIN` needs the cheap dither and a narrower frame.** The
//! camera is otherwise still, so the only thing moving is the plot and a GIF
//! compresses that beautifully; spinning changes every pixel of every frame.
//! Measured: the recipe above plus `ISOMESH_SPIN=0.004` comes back at **5.15 MB**,
//! over the 4.8 MB the committed clips sit within, and drops to **4.04 MB** with
//! `WIDTH=800 DITHER=bayer:bayer_scale=3`.
//!
//! Demonstrates **M-341 / P-39**. A sculpting stroke lays down 64 brushes over a
//! solid ball split into 64 chunks. Every chunk is re-meshed continuously, and
//! the plot in the corner is the whole point: **the red bars — meshing against
//! the whole tape — climb with every brush; the green bars — meshing against the
//! brushes that can still matter inside that chunk — climb two and a half times
//! slower.**
//!
//! # The green line is not flat, and pretending it were would be the lie
//!
//! The ticket that asked for this demo wanted the pruned cost to be *flat*, and
//! on this fixture it is not, for a reason worth stating rather than tuning
//! away. P-39's tape is **uniformly scattered over the world**: doubling the
//! tape doubles the brush density everywhere, so it doubles the number of
//! brushes overlapping any given chunk too. A constant *fraction* survives —
//! measured, about 0.30 — so the pruned cost still grows linearly, at 0.30 of
//! the slope. What the mechanism buys is the **2.5× divisor**, not a constant.
//!
//! A flat line is what a *moving* stroke would give, because then a chunk's
//! survivor count saturates instead of growing: the brushes laid down after the
//! tool left are all provably losing there. That is the game case and it is
//! strictly kinder to the mechanism. This example keeps P-39's uniform scatter,
//! which is the registered fixture and the harder one, so its numbers can be
//! held against M-341's rather than admired on their own.
//!
//! # The mechanism, in one sample per brush per chunk
//!
//! [`BrushStack::sample`](isomesh::brush::BrushStack) is a linear fold over
//! every brush, and Marching Cubes prefills the entire sample grid before it
//! looks at a single cell, so a 33³ chunk walks the whole edit history 35,937
//! times. Most of that history cannot possibly matter inside any one chunk.
//!
//! A shape with declared Lipschitz constant `l` varies by at most `l·r` over a
//! box of circumradius `r`, so `f(centre) ± l·r` encloses it there. `l` is not
//! hard-coded: it comes from
//! [`BoundedSdf::value_bound`](isomesh::fields::BoundedSdf) via
//! [`FieldBound::lipschitz`](isomesh::fields::FieldBound::lipschitz), which
//! answers `1` for every exact distance field — `Sphere` and `Capsule` both. A
//! brush whose enclosure is strictly clear of the running fold's enclosure loses
//! the `min`/`max` chain at every point in the chunk, so it can be **deleted**
//! from the tape rather than merely skipped, and the pruned tape is a shorter
//! slice of the same `&[Brush<S>]`.
//!
//! # Deleted bit-exactly, which is the claim the speed rests on
//!
//! `apply(Add, f, s)` is IEEE `min(f, s)` and `apply(Subtract, f, s)` is IEEE
//! `max(f, −s)`. Both **select** an operand rather than computing a new value,
//! and negation is exact, so dropping a provably-losing `Add` or `Subtract`
//! moves the result by exactly zero ULP. The HUD's `mesh identical` line is that
//! claim, checked live: every measured chunk is meshed both ways every visit and
//! the two buffers compared **bit for bit** — `to_bits()`, not `==`, because
//! `-0.0 == 0.0` and `-0.0` is precisely the value the one soft spot in the
//! selection argument would produce.
//!
//! The losing test is **strict** (`s.lo > v.hi`) for the same reason. `f32::min`
//! is documented to return *either* operand when they compare equal, which is
//! observable only for `+0.0` against `-0.0`; a strict inequality leaves no tie
//! to resolve. A non-strict test would put a signed-zero hole in a bit-exactness
//! claim for a pruning gain of measure zero.
//!
//! # `SmoothAdd` is not prunable in the losing direction
//!
//! The asymmetry is load-bearing, so [`prune_into`]'s
//! [`SmoothAdd`](isomesh::brush::BrushOp::SmoothAdd) arm exists and never
//! prunes. `smooth_min` at `h == 1` returns `b + (a − b)`, which is *not*
//! bit-identical to `a` — so a smooth brush that provably loses still cannot be
//! deleted without moving the mesh, and the mechanism has nothing to say about
//! it. This example's tape is `Add` and `Subtract` only, which is P-39's
//! registered fixture; the arm is there so the rule is structural rather than a
//! comment.
//!
//! # Two details that make the enclosure an enclosure
//!
//! - **The box is bigger than the sample grid.** Marching Cubes' normals come
//!   from [`Sdf::gradient`], and `BrushStack` does not override it, so the
//!   default central differences sample `DIFF_STEP · max(|coord|, 1)` *outside*
//!   the sampled extent. A brush pruned on a bound that stopped at the grid
//!   could move a normal. [`ChunkBox::new`] inflates by exactly that margin.
//! - **The bound arithmetic is rounded.** The Lipschitz inequality is about the
//!   exact function; `sample(centre)` and `f(c) ± l·r` are both floating point,
//!   and a bound rounded the wrong way is not a bound. Every enclosure is
//!   widened by [`PAD_ULPS`] ULP of the magnitudes involved, which also covers
//!   the few-ULP evaluation error of these closed forms.
//!
//! # What the numbers should say, and what they are compared against
//!
//! M-341 measured this fixture — the same 4×4×4 world of 33-sample chunks, the
//! same 64-brush tape, the same seed — in `f64`, offline, five timed reps per
//! arm per chunk: **median surviving fraction 0.2969, median per-chunk speedup
//! 3.365×, world aggregate 2.473×, range 0.992× to 22.47×, mesh byte-identical
//! on 64 of 64 chunks.** The bound cost 540–1450 ns per chunk for a 64-brush
//! tape, which is 4.4e-5 to 1.3e-3 of the meshing it enables.
//!
//! This runs the same fixture in `f32`, live, one timed run per arm per visit
//! (kept as a best-of over repeat visits at the same tape length, which is how a
//! single-shot timing gets denoised without paying for five). Settled at tape 64
//! after seventeen sweeps: **world aggregate 2.450×, median per-chunk 3.120×,
//! range 0.999× to 25.03×, median survivor fraction 0.3281, 64 of 64 chunks
//! bit-identical, bound cost 0.99 µs per chunk = 1.53e-4 of the meshing it
//! enables, 6.458 ms pruned against 15.824 ms full, 45 fps.** So the aggregate —
//! the one number a game would feel — reproduces M-341 to within 1%.
//!
//! **The survivor fraction is 0.3281 against M-341's 0.2969, and the whole gap
//! is the `f32` gradient margin.** The enclosure has to cover wherever
//! [`Sdf::gradient`] samples, and its central-difference step is
//! `DIFF_STEP · max(|coord|, 1)` — `4.92e-3` in `f32` against `6.06e-6` in
//! `f64`, because `DIFF_STEP` is the cube root of the type's epsilon. At the far
//! corner of the world that inflates a chunk's circumradius from 3.4641 to
//! 3.5324, a 2% looser bound.
//!
//! Measured rather than reasoned: substituting `f64`'s `DIFF_STEP` into
//! [`ChunkBox::new`] and changing nothing else moves the median survivor
//! fraction to **0.2969 — M-341's figure exactly** — and the world aggregate to
//! 2.513×. That substitution is not committed, because in `f32` it is not a
//! bound: the differences really do reach `4.92e-3 · |p|` outside the grid, and
//! a bound that does not cover them can move a normal.
//!
//! # Why the whole world is re-meshed continuously
//!
//! An engine would re-mesh only what an edit touched — `game_dig` is that
//! example, and [`chunk::dirty`](isomesh::chunk::dirty) is that machinery. Here
//! the *cost of meshing a chunk against a tape* is the subject, so a cursor
//! sweeps the chunk grid at a fixed time budget per frame and re-measures. That
//! makes the per-chunk cost visible everywhere at once, and it is why the world
//! aggregate is a **rolling** figure: each chunk contributes its most recent
//! measurement, so the aggregate trails the tape by up to one sweep. It settles
//! once the stroke finishes and the tape holds at 64.
//!
//! # Spacing is a power of two, deliberately
//!
//! `h = 0.125`. **M-32** measured that two chunks agree on their shared sample
//! plane bit-for-bit only at a power-of-two cell size; anywhere else they differ
//! by an ulp and the seam needs a weld. Nothing here welds, so the seam has to
//! be exact on its own.

mod common;

use std::time::Instant;

use bevy::prelude::*;
use bevy_isomesh::to_bevy_mesh;
use common::{Capture, CommonPlugin, DemoDomain, DemoMesh, DemoStats, OrbitCamera, ViewFlags};
use isomesh::brush::{Brush, BrushOp, BrushStack, Capsule};
use isomesh::chunk::{ChunkId, ChunkLayout};
use isomesh::fields::{BoundedSdf, FieldBound, Sphere};
use isomesh::marching_cubes::MarchingCubes;
use isomesh::{MeshBuffer, Real, Sdf};

// ── the fixture, which is P-39's ────────────────────────────────────────────

/// Chunks along each axis. 4³ = 64 chunks, as M-341 measured.
const CHUNKS_PER_AXIS: i32 = 4;

/// Cells per chunk axis, so 33 samples per axis.
const CELLS_PER_CHUNK: u32 = 32;

/// World units per cell. A power of two — see the module docs.
const CELL_SIZE: f32 = 0.125;

/// The world's minimum corner, so the world is `[-8, 8]³`.
const WORLD_ORIGIN: f32 = -8.0;

/// Radius of the solid the brushes carve. Its surface crosses every chunk but
/// the eight corners, so the fixture has interior, surface and empty chunks
/// without being arranged to.
const BASE_RADIUS: f32 = 6.0;

/// Brushes the stroke lays down.
const TAPE_MAX: usize = 64;

/// ULP of slack added to every enclosure bound.
///
/// Four times P-39's sixteen, because this runs in `f32`: the evaluation error
/// of a capsule distance is a few ULP of a magnitude near 30 at the far corner
/// of the world, and 64 ULP of that is 2.4e-4 against a Lipschitz reach of 3.5.
/// It buys a twenty-fold margin and costs a pruning decision only exactly on the
/// boundary.
const PAD_ULPS: f32 = 64.0;

// ── pacing ──────────────────────────────────────────────────────────────────

/// Milliseconds of measurement work the sweep may do per frame.
///
/// A budget rather than a fixed chunk count, so the sweep uses whatever headroom
/// the machine has instead of a number guessed on one of them: on an empty tape
/// it laps the world in a few frames, and at tape 64 it settles to one chunk a
/// frame. Both arms count against it, because both ran.
///
/// Deliberately most of a frame. Meshing a chunk against a 64-brush tape is
/// what the demo is about, and hiding it behind a small budget would make the
/// sweep crawl rather than make the frame cheap.
const MEASURE_BUDGET_MS: f64 = 14.0;

/// Hard cap on chunks per frame, so an empty tape cannot spend the whole sweep
/// in one frame and skip the animation.
const MAX_CHUNKS_PER_FRAME: usize = 6;

/// Seconds between brushes when nobody is capturing.
const STROKE_INTERVAL: f32 = 0.14;

// ── the plot ────────────────────────────────────────────────────────────────

/// One column per tape length, `0..=TAPE_MAX`.
const PLOT_COLS: usize = TAPE_MAX + 1;

/// Bar width and column stride, in logical pixels.
const PLOT_BAR_W: f32 = 3.0;
const PLOT_STRIDE: f32 = 4.0;

/// Plot height, and where its bottom-left corner sits.
const PLOT_H: f32 = 92.0;
const PLOT_LEFT: f32 = 14.0;
const PLOT_BOTTOM: f32 = 26.0;

/// Full-scale of the plot's y axis, in milliseconds per chunk re-mesh.
///
/// Fixed rather than auto-scaled, and that is the honest choice: a scale that
/// grew with the red bars would shrink the green ones on screen while their
/// real cost had not moved, which is the opposite of what the plot is for.
/// Measured on this machine at tape 64: 15.8 ms per chunk against the whole
/// tape and 6.4 ms against the pruned one, so 18 fits both with a little
/// headroom. Bars clamp rather than overflow.
const PLOT_MAX_MS: f64 = 18.0;

// ── shapes ──────────────────────────────────────────────────────────────────

/// A brush shape in this demo.
///
/// One enum so the whole stack is a single `&[Brush<Shape>]` slice, which is
/// what makes a pruned tape a shorter slice of the same type. Both variants are
/// exact distance fields, so both declare `l == 1` — and they declare it through
/// the crate's own [`BoundedSdf`] rather than by a constant written here.
#[derive(Clone, Copy, Debug)]
enum Shape {
    /// A ball.
    Sphere(Sphere<f32>),
    /// A swept segment.
    Capsule(Capsule<f32>),
}

impl Sdf for Shape {
    type Scalar = f32;

    fn sample(&self, p: [f32; 3]) -> f32 {
        match self {
            Self::Sphere(s) => s.sample(p),
            Self::Capsule(c) => c.sample(p),
        }
    }
}

impl BoundedSdf for Shape {
    fn value_bound(&self) -> FieldBound {
        match self {
            Self::Sphere(s) => s.value_bound(),
            Self::Capsule(c) => c.value_bound(),
        }
    }
}

/// A 64-bit LCG, so the 64 brushes are the same 64 brushes on every machine and
/// every run.
///
/// Fixture construction, not output: nothing measured depends on the generator
/// being good, only on it being reproducible. It draws in `f64` and rounds to
/// `f32` at the end, so the brushes are P-39's brushes to the last bit `f32` can
/// hold and the survivor counts are comparable against M-341's.
struct Lcg(u64);

impl Lcg {
    /// Seeded.
    const fn new(seed: u64) -> Self {
        Self(seed)
    }

    /// The next raw draw.
    fn next_u32(&mut self) -> u32 {
        self.0 = self
            .0
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        (self.0 >> 33) as u32
    }

    /// A float in `[0, 1)`, from 24 bits so it is exactly representable.
    fn unit(&mut self) -> f64 {
        f64::from(self.next_u32() >> 8) / 16_777_216.0
    }

    /// A float in `[lo, hi)`.
    fn range(&mut self, lo: f64, hi: f64) -> f64 {
        lo + (hi - lo) * self.unit()
    }
}

/// The 64-brush tape: `Add` and `Subtract` over spheres and capsules, scattered
/// across the whole world.
///
/// Uniform over the world rather than concentrated on the solid's surface. That
/// is the *harder* fixture for the mechanism — surface-concentrated edits leave
/// whole chunks touched by nothing, which prunes trivially — and it is the
/// reading of "scattered over a 4×4×4 chunk world" that cannot be accused of
/// arranging the answer.
fn tape() -> Vec<Brush<Shape>> {
    // P-39's seed, and nothing about the result depends on the value.
    let mut rng = Lcg::new(0x39_5EED_C0DE_1234);
    let mut out = Vec::with_capacity(TAPE_MAX);
    for _ in 0..TAPE_MAX {
        let centre = [
            rng.range(-6.5, 6.5),
            rng.range(-6.5, 6.5),
            rng.range(-6.5, 6.5),
        ];
        let shape = if rng.next_u32() & 1 == 0 {
            Shape::Sphere(Sphere {
                center: [centre[0] as f32, centre[1] as f32, centre[2] as f32],
                radius: rng.range(0.35, 1.1) as f32,
            })
        } else {
            let dir = [
                rng.range(-1.0, 1.0),
                rng.range(-1.0, 1.0),
                rng.range(-1.0, 1.0),
            ];
            let len = (dir[0] * dir[0] + dir[1] * dir[1] + dir[2] * dir[2]).sqrt();
            // A zero direction would be a sphere, which is a fine capsule but
            // not the one this row is meant to contribute.
            let unit = if len > 1e-9 {
                [dir[0] / len, dir[1] / len, dir[2] / len]
            } else {
                [1.0, 0.0, 0.0]
            };
            let half = rng.range(0.25, 1.0);
            let end = |sign: f64| {
                [
                    (centre[0] + sign * unit[0] * half) as f32,
                    (centre[1] + sign * unit[1] * half) as f32,
                    (centre[2] + sign * unit[2] * half) as f32,
                ]
            };
            Shape::Capsule(Capsule {
                a: end(-1.0),
                b: end(1.0),
                radius: rng.range(0.3, 0.8) as f32,
            })
        };
        let op = if rng.next_u32() & 1 == 0 {
            BrushOp::Add
        } else {
            BrushOp::Subtract
        };
        out.push(Brush { shape, op });
    }
    out
}

// ── the bound ───────────────────────────────────────────────────────────────

/// The enclosure of a scalar field over a chunk.
#[derive(Clone, Copy, Debug)]
struct Interval {
    /// Lower bound.
    lo: f32,
    /// Upper bound.
    hi: f32,
}

/// The box a chunk's field evaluations can touch.
#[derive(Clone, Copy, Debug)]
struct ChunkBox {
    /// Centre of the sampled extent.
    centre: [f32; 3],
    /// Circumradius: half the diagonal of the sampled extent, plus the margin
    /// [`Sdf::gradient`]'s central differences reach outside it.
    radius: f32,
}

impl ChunkBox {
    /// The box for a chunk whose sample grid starts at `origin` and spans `span`
    /// on every axis.
    fn new(origin: [f32; 3], span: f32) -> Self {
        let centre = [
            origin[0] + span * 0.5,
            origin[1] + span * 0.5,
            origin[2] + span * 0.5,
        ];
        // `h = DIFF_STEP * max(|p|, 1)` at the furthest corner bounds the
        // differencing reach anywhere in the box.
        let mut far = 1.0f32;
        for lo in origin {
            far = far.max(lo.abs()).max((lo + span).abs());
        }
        let half = span * 0.5 + <f32 as Real>::DIFF_STEP * far;
        Self {
            centre,
            radius: half * 3.0f32.sqrt(),
        }
    }
}

/// Slack for one bound, in absolute units.
fn pad(value: f32, reach: f32) -> f32 {
    PAD_ULPS * f32::EPSILON * (value.abs() + reach)
}

/// `f(centre) ± l·r`, widened so it is an enclosure and not an estimate.
///
/// A field that declares no Lipschitz constant gets an infinite enclosure, which
/// makes it unprunable rather than wrongly prunable. Nothing in this fixture
/// takes that path — both shapes are exact distances — but the alternative would
/// be an `expect` in a demo a stranger runs.
fn enclose<S: BoundedSdf<Scalar = f32>>(field: &S, chunk: ChunkBox) -> Interval {
    let Some(l) = field.value_bound().lipschitz() else {
        return Interval {
            lo: f32::NEG_INFINITY,
            hi: f32::INFINITY,
        };
    };
    let value = field.sample(chunk.centre);
    let reach = l as f32 * chunk.radius;
    let slack = reach + pad(value, reach);
    Interval {
        lo: value - slack,
        hi: value + slack,
    }
}

/// Select the brushes that can still change the fold anywhere in `chunk`.
///
/// Order is preserved, because `Add` and `Subtract` do not commute with each
/// other — [`BrushOp::commutes_with`] is the crate's own statement of that.
/// Returns the number of survivors, which is `out.len()`.
fn prune_into(
    tape: &[Brush<Shape>],
    base: &Sphere<f32>,
    chunk: ChunkBox,
    out: &mut Vec<Brush<Shape>>,
) -> usize {
    out.clear();
    let mut v = enclose(base, chunk);
    for brush in tape {
        let s = enclose(&brush.shape, chunk);
        match brush.op {
            BrushOp::Add => {
                // Strictly above the running value everywhere, so `min` is
                // forced to select the running value at every point.
                if s.lo > v.hi {
                    continue;
                }
                v = Interval {
                    lo: v.lo.min(s.lo),
                    hi: v.hi.min(s.hi),
                };
            }
            BrushOp::Subtract => {
                // `max(v, -s)`. Negation is exact, so negating the enclosure is
                // exact too.
                let n = Interval {
                    lo: -s.hi,
                    hi: -s.lo,
                };
                if n.hi < v.lo {
                    continue;
                }
                v = Interval {
                    lo: v.lo.max(n.lo),
                    hi: v.hi.max(n.hi),
                };
            }
            BrushOp::SmoothAdd { k } => {
                // **Never pruned in the losing direction**, and that is the
                // registered asymmetry rather than an omission: `smooth_min` at
                // `h == 1` returns `b + (a - b)`, which is not bit-identical to
                // `a`. The enclosure still has to track it, and its floor sags
                // by at most `k/4` below the plain minimum.
                let k = k as f32;
                let lo = v.lo.min(s.lo) - 0.25 * k;
                v = Interval {
                    lo: lo - pad(lo, 0.25 * k),
                    hi: v.hi.min(s.hi),
                };
            }
        }
        out.push(*brush);
    }
    out.len()
}

// ── measurement ─────────────────────────────────────────────────────────────

/// Bit-for-bit equality of two meshes.
///
/// `f32 == f32` calls `-0.0` equal to `0.0`, and `-0.0` is exactly the value the
/// selection lemma's one soft spot would produce, so this compares bits.
fn bitwise_identical(a: &MeshBuffer<f32>, b: &MeshBuffer<f32>) -> bool {
    fn same(x: &[f32; 3], y: &[f32; 3]) -> bool {
        x[0].to_bits() == y[0].to_bits()
            && x[1].to_bits() == y[1].to_bits()
            && x[2].to_bits() == y[2].to_bits()
    }
    a.indices == b.indices
        && a.positions.len() == b.positions.len()
        && a.normals.len() == b.normals.len()
        && a.positions
            .iter()
            .zip(&b.positions)
            .all(|(x, y)| same(x, y))
        && a.normals.iter().zip(&b.normals).all(|(x, y)| same(x, y))
}

/// The reusable extraction machinery, kept out of the timed region.
///
/// One extractor and one buffer per arm, reused across chunks, so no allocation
/// lands inside a measurement. Grouped in its own struct purely so a single
/// `&mut Demo` can hand it to [`measure_chunk`] while the tape stays borrowed
/// immutably.
struct Rig {
    /// The extractor, reused so its private grids are not reallocated.
    mc: MarchingCubes<f32>,
    /// Output of the pruned arm.
    pruned: MeshBuffer<f32>,
    /// Output of the full-tape arm.
    full: MeshBuffer<f32>,
    /// The pruned tape, reused.
    survivors: Vec<Brush<Shape>>,
}

/// What one chunk visit found.
#[derive(Clone, Copy)]
struct Measured {
    /// Brushes that survived the bound.
    survivors: usize,
    /// Milliseconds to mesh against the pruned tape.
    ms_pruned: f64,
    /// Milliseconds to mesh against the whole tape.
    ms_full: f64,
    /// Microseconds the whole pruning pass cost.
    bound_us: f64,
    /// Whether the two arms agreed bit for bit.
    identical: bool,
}

/// Mesh one chunk both ways and time both.
///
/// Leaves the pruned arm's mesh in `rig.pruned` and the full arm's in
/// `rig.full`; the caller uploads whichever the toggle says is active. `None`
/// means the extractor refused the chunk, which is logged and skipped rather
/// than unwrapped.
fn measure_chunk(
    rig: &mut Rig,
    layout: &ChunkLayout<f32>,
    base: Sphere<f32>,
    tape: &[Brush<Shape>],
    id: ChunkId,
) -> Option<Measured> {
    let shape = layout.sample_shape().ok()?;
    let origin = layout.sample_origin(id);
    let span = layout.cell_size() * layout.cells() as f32;
    let chunk = ChunkBox::new(origin, span);

    let started = Instant::now();
    let survivors = prune_into(tape, &base, chunk, &mut rig.survivors);
    let bound_us = started.elapsed().as_secs_f64() * 1e6;

    let pruned_field = BrushStack {
        base,
        brushes: &rig.survivors,
    };
    let started = Instant::now();
    rig.pruned.reset();
    rig.mc
        .extract(
            &pruned_field,
            &shape,
            origin,
            layout.cell_size(),
            &mut rig.pruned,
        )
        .ok()?;
    let ms_pruned = started.elapsed().as_secs_f64() * 1e3;

    let full_field = BrushStack { base, brushes: tape };
    let started = Instant::now();
    rig.full.reset();
    rig.mc
        .extract(
            &full_field,
            &shape,
            origin,
            layout.cell_size(),
            &mut rig.full,
        )
        .ok()?;
    let ms_full = started.elapsed().as_secs_f64() * 1e3;

    Some(Measured {
        survivors,
        ms_pruned,
        ms_full,
        bound_us,
        identical: bitwise_identical(&rig.pruned, &rig.full),
    })
}

/// One chunk's most recent measurement.
///
/// `ms_*` are the **best** seen at this tape length, not the latest. A single
/// timed run of a 5 ms extraction carries a scheduler tick's worth of noise;
/// taking the minimum over repeat visits converges on the real cost without
/// paying for five reps up front, and resetting on a tape change keeps it from
/// carrying a shorter tape's number forward.
#[derive(Clone, Copy)]
struct Record {
    /// Whether this chunk has ever been visited.
    measured: bool,
    /// Tape length these numbers were taken at.
    tape_len: usize,
    /// Visits at that tape length.
    reps: u32,
    /// Surviving brushes.
    survivors: usize,
    /// Best pruned-arm milliseconds.
    ms_pruned: f64,
    /// Best full-arm milliseconds.
    ms_full: f64,
    /// Latest bound cost, microseconds.
    bound_us: f64,
    /// Whether the arms agreed on the last visit.
    identical: bool,
    /// Vertices in the active mesh.
    vertices: usize,
    /// Triangles in the active mesh.
    triangles: usize,
    /// Survivor fraction the material was last painted with, so 64 uniform
    /// uploads a frame are not spent repainting the same colour.
    painted: f32,
}

impl Record {
    /// A chunk nothing is known about yet.
    const fn blank() -> Self {
        Self {
            measured: false,
            tape_len: usize::MAX,
            reps: 0,
            survivors: 0,
            ms_pruned: 0.0,
            ms_full: 0.0,
            bound_us: 0.0,
            identical: true,
            vertices: 0,
            triangles: 0,
            painted: f32::NAN,
        }
    }

    /// Fold a fresh visit in.
    fn absorb(&mut self, m: Measured, tape_len: usize) {
        if !self.measured || self.tape_len != tape_len {
            self.tape_len = tape_len;
            self.reps = 1;
            self.ms_pruned = m.ms_pruned;
            self.ms_full = m.ms_full;
        } else {
            self.reps += 1;
            self.ms_pruned = self.ms_pruned.min(m.ms_pruned);
            self.ms_full = self.ms_full.min(m.ms_full);
        }
        self.measured = true;
        self.survivors = m.survivors;
        self.bound_us = m.bound_us;
        self.identical = m.identical;
    }

    /// Milliseconds the active arm cost.
    fn active_ms(&self, prune: bool) -> f64 {
        if prune { self.ms_pruned } else { self.ms_full }
    }
}

// ── state ───────────────────────────────────────────────────────────────────

/// Gizmos for the chunk boxes.
#[derive(Default, Reflect, GizmoConfigGroup)]
struct ChunkGizmos;

/// One bar of the plot, and which arm it draws.
#[derive(Component)]
struct PlotBar {
    /// Tape length this column stands for.
    column: usize,
    /// `true` for the full-tape bar, `false` for the pruned one.
    full: bool,
}

/// The whole demo.
#[derive(Resource)]
struct Demo {
    /// Chunk lattice.
    layout: ChunkLayout<f32>,
    /// The solid the brushes carve.
    base: Sphere<f32>,
    /// All 64 brushes. The live tape is the first [`Demo::tape_len`] of them.
    tape: Vec<Brush<Shape>>,
    /// How much of the tape the stroke has laid down.
    tape_len: usize,
    /// Extraction machinery.
    rig: Rig,
    /// Sweep order.
    ids: Vec<ChunkId>,
    /// One entity per chunk, parallel to [`Demo::ids`].
    entities: Vec<Entity>,
    /// One mesh asset per chunk that currently has a surface, overwritten in
    /// place rather than replaced. `None` for a chunk with no triangles — see
    /// [`sweep`] for why an empty mesh is not an option.
    meshes: Vec<Option<Handle<Mesh>>>,
    /// One material per chunk, so the overlay can tint each independently.
    materials: Vec<Handle<StandardMaterial>>,
    /// One record per chunk.
    records: Vec<Record>,
    /// Where the sweep is.
    cursor: usize,
    /// The chunk the sweep measured most recently.
    last_visited: usize,
    /// Completed sweeps.
    sweeps: u32,
    /// Whether the tape is pruned before meshing.
    prune: bool,
    /// Whether the survivor-fraction tint is on.
    heat: bool,
    /// Whether the chunk boxes are drawn.
    boxes: bool,
    /// The first chunk whose two arms disagreed, if any ever did.
    mismatch: Option<[i32; 3]>,
    /// Per tape length: mean milliseconds per chunk for the full and pruned
    /// arms, as of the moment the stroke reached that length.
    plot: Vec<Option<(f64, f64)>>,
    /// Highest column written, so a stroke that advances by more than one brush
    /// a frame leaves no gaps.
    plot_upto: usize,
    /// Seconds since the last brush, when nobody is capturing.
    stroke_timer: f32,
    /// Brushes per captured frame, so a capture of any length sees the whole
    /// stroke instead of a still.
    per_capture_frame: usize,
}

impl Demo {
    /// Chunks in the world.
    fn chunks(&self) -> usize {
        self.ids.len()
    }

    /// Chunks visited at least once.
    fn measured(&self) -> usize {
        self.records.iter().filter(|r| r.measured).count()
    }

    /// Rolling world totals: milliseconds for the full and the pruned arm, each
    /// chunk contributing its most recent measurement.
    fn world_ms(&self) -> (f64, f64) {
        let mut full = 0.0;
        let mut pruned = 0.0;
        for r in self.records.iter().filter(|r| r.measured) {
            full += r.ms_full;
            pruned += r.ms_pruned;
        }
        (full, pruned)
    }

    /// Mean milliseconds per measured chunk, full arm then pruned arm.
    fn mean_ms(&self) -> (f64, f64) {
        let n = self.measured();
        if n == 0 {
            return (0.0, 0.0);
        }
        let (full, pruned) = self.world_ms();
        (full / n as f64, pruned / n as f64)
    }

    /// Survivor counts across measured chunks, sorted.
    fn survivor_spread(&self) -> Vec<usize> {
        let mut v: Vec<usize> = self
            .records
            .iter()
            .filter(|r| r.measured)
            .map(|r| r.survivors)
            .collect();
        v.sort_unstable();
        v
    }

    /// Per-chunk speedups across measured chunks, sorted.
    fn speedups(&self) -> Vec<f64> {
        let mut v: Vec<f64> = self
            .records
            .iter()
            .filter(|r| r.measured && r.ms_pruned > 0.0)
            .map(|r| r.ms_full / r.ms_pruned)
            .collect();
        v.sort_by(f64::total_cmp);
        v
    }

    /// Forget every measurement, so the next sweep re-takes them all.
    fn forget(&mut self) {
        for r in &mut self.records {
            let painted = r.painted;
            *r = Record::blank();
            r.painted = painted;
        }
    }
}

/// Median of a sorted slice, or zero when it is empty.
fn median<T: Copy + Default>(sorted: &[T]) -> T {
    sorted.get(sorted.len() / 2).copied().unwrap_or_default()
}

// ── app ─────────────────────────────────────────────────────────────────────

fn main() {
    App::new()
        .add_plugins(DefaultPlugins.set(WindowPlugin {
            primary_window: Some(Window {
                title: "isomesh - E-305 Lipschitz tape pruning".into(),
                ..default()
            }),
            ..default()
        }))
        .add_plugins(CommonPlugin)
        .init_gizmo_group::<ChunkGizmos>()
        .add_systems(Startup, setup)
        .add_systems(
            Update,
            (
                controls,
                advance_stroke,
                sweep,
                paint_heat,
                update_plot,
                report,
                draw_chunk_boxes,
            )
                .chain(),
        )
        .run();
}

/// Brushes to add per captured frame, so a six-frame smoke capture still shows
/// the tape grow and a ninety-frame one still shows it grow smoothly.
fn per_capture_frame() -> usize {
    let frames: usize = std::env::var("ISOMESH_CAPTURE_FRAMES")
        .ok()
        .and_then(|v| v.parse().ok())
        // The harness's own default. Read rather than shared, because `Capture`
        // keeps its length private and the alternative is editing the harness.
        .unwrap_or(60);
    // Aim to finish the stroke around halfway through the capture, so the tail
    // of the clip shows the settled numbers rather than cutting on the last
    // brush.
    TAPE_MAX.div_ceil(frames.max(2) / 2 + 1).max(1)
}

fn setup(
    mut commands: Commands,
    mut materials: ResMut<Assets<StandardMaterial>>,
    mut config: ResMut<GizmoConfigStore>,
    camera: Query<Entity, With<OrbitCamera>>,
) {
    let Ok(layout) = ChunkLayout::<f32>::new(
        CELLS_PER_CHUNK,
        CELL_SIZE,
        [WORLD_ORIGIN, WORLD_ORIGIN, WORLD_ORIGIN],
    ) else {
        error!("chunk layout rejected {CELLS_PER_CHUNK} cells at {CELL_SIZE}");
        return;
    };

    for entity in &camera {
        commands.entity(entity).insert(OrbitCamera {
            focus: Vec3::ZERO,
            // Far enough that the 16-unit chunk lattice stays inside the frame,
            // close enough that a brush-sized dent in the ball is visible at the
            // 900px the GIFs are scaled to.
            radius: 32.0,
            yaw: 0.7,
            pitch: 0.34,
        });
    }
    let (chunk_gizmos, _) = config.config_mut::<ChunkGizmos>();
    chunk_gizmos.line.width = 1.0;

    let span = CELL_SIZE * CELLS_PER_CHUNK as f32;
    let extent = span * CHUNKS_PER_AXIS as f32;
    commands.spawn(DemoDomain {
        min: Vec3::splat(WORLD_ORIGIN),
        max: Vec3::splat(WORLD_ORIGIN + extent),
    });

    let mut ids = Vec::new();
    for z in 0..CHUNKS_PER_AXIS {
        for y in 0..CHUNKS_PER_AXIS {
            for x in 0..CHUNKS_PER_AXIS {
                ids.push(ChunkId::new([x, y, z]));
            }
        }
    }

    let mut chunk_entities = Vec::with_capacity(ids.len());
    let mut chunk_materials = Vec::with_capacity(ids.len());
    for _ in 0..ids.len() {
        let material = materials.add(StandardMaterial {
            base_color: Color::srgb(0.72, 0.76, 0.82),
            perceptual_roughness: 0.5,
            metallic: 0.04,
            ..default()
        });
        // `Mesh3d::default()` names no asset, so nothing is drawn and nothing
        // is uploaded until the first sweep gives this chunk a surface. An
        // empty mesh would be worse than nothing here -- see `sweep`.
        chunk_entities.push(
            commands
                .spawn((
                    Mesh3d::default(),
                    MeshMaterial3d(material.clone()),
                    DemoMesh,
                ))
                .id(),
        );
        chunk_materials.push(material);
    }

    spawn_plot(&mut commands);

    let requested = std::env::var("ISOMESH_VIEW").unwrap_or_default();
    let has = |name: &str| requested.split(',').any(|part| part.trim() == name);

    commands.insert_resource(Demo {
        layout,
        base: Sphere {
            center: [0.0; 3],
            radius: BASE_RADIUS,
        },
        tape: tape(),
        tape_len: 0,
        rig: Rig {
            mc: MarchingCubes::<f32>::new(),
            pruned: MeshBuffer::<f32>::new(),
            full: MeshBuffer::<f32>::new(),
            survivors: Vec::with_capacity(TAPE_MAX),
        },
        records: vec![Record::blank(); ids.len()],
        meshes: vec![None; ids.len()],
        ids,
        entities: chunk_entities,
        materials: chunk_materials,
        cursor: 0,
        last_visited: 0,
        sweeps: 0,
        prune: !has("noprune"),
        // Not gated on `flags.hud`: the tint is geometry, not text, and a clip
        // recorded with `nohud` still wants the spatial pattern.
        heat: !has("noheat"),
        boxes: !has("noboxes"),
        mismatch: None,
        plot: vec![None; PLOT_COLS],
        plot_upto: 0,
        stroke_timer: 0.0,
        per_capture_frame: per_capture_frame(),
    });
}

/// The plot: one red bar per tape length for the full tape, one green bar in
/// front of it for the pruned one.
///
/// Root-level absolutely-positioned nodes rather than a flex row, so a column's
/// x is arithmetic rather than a layout outcome, and [`GlobalZIndex`] rather than
/// spawn order, so the stacking is stated rather than inherited.
fn spawn_plot(commands: &mut Commands) {
    commands.spawn((
        Node {
            position_type: PositionType::Absolute,
            left: Val::Px(PLOT_LEFT - 6.0),
            bottom: Val::Px(PLOT_BOTTOM - 6.0),
            width: Val::Px(PLOT_COLS as f32 * PLOT_STRIDE + 12.0),
            height: Val::Px(PLOT_H + 12.0),
            ..default()
        },
        BackgroundColor(Color::srgba(0.05, 0.06, 0.09, 0.62)),
        GlobalZIndex(1),
        PlotPanel,
    ));
    commands.spawn((
        Text::new(format!(
            "ms per chunk re-mesh vs tape length   red: whole tape   green: pruned   (0-{PLOT_MAX_MS:.0} ms)"
        )),
        TextFont {
            font_size: FontSize::Px(11.0),
            ..default()
        },
        TextColor(Color::srgb(0.80, 0.84, 0.90)),
        Node {
            position_type: PositionType::Absolute,
            left: Val::Px(PLOT_LEFT - 6.0),
            bottom: Val::Px(PLOT_BOTTOM + PLOT_H + 9.0),
            ..default()
        },
        GlobalZIndex(2),
        PlotPanel,
    ));
    for column in 0..PLOT_COLS {
        let left = PLOT_LEFT + column as f32 * PLOT_STRIDE;
        for full in [true, false] {
            commands.spawn((
                Node {
                    position_type: PositionType::Absolute,
                    left: Val::Px(left),
                    bottom: Val::Px(PLOT_BOTTOM),
                    width: Val::Px(PLOT_BAR_W),
                    height: Val::Px(0.0),
                    ..default()
                },
                BackgroundColor(if full {
                    Color::srgb(0.92, 0.30, 0.24)
                } else {
                    Color::srgb(0.24, 0.86, 0.42)
                }),
                GlobalZIndex(if full { 2 } else { 3 }),
                PlotBar { column, full },
                PlotPanel,
            ));
        }
    }
}

/// Everything the plot is made of, so `nohud` can hide all of it at once.
#[derive(Component)]
struct PlotPanel;

/// `P` pruning, `H` heat overlay, `C` chunk boxes, `X` restart the stroke.
fn controls(keys: Res<ButtonInput<KeyCode>>, flags: Res<ViewFlags>, demo: Option<ResMut<Demo>>) {
    let Some(mut demo) = demo else { return };
    if keys.just_pressed(KeyCode::KeyP) {
        demo.prune = !demo.prune;
        // The arms cost different amounts, and a record's best-of is per tape
        // length rather than per arm, so the toggle has to invalidate.
        demo.forget();
    }
    if keys.just_pressed(KeyCode::KeyH) {
        demo.heat = !demo.heat;
    }
    if keys.just_pressed(KeyCode::KeyC) {
        demo.boxes = !demo.boxes;
    }
    if keys.just_pressed(KeyCode::KeyX) {
        demo.tape_len = 0;
        demo.stroke_timer = 0.0;
        demo.plot.iter_mut().for_each(|slot| *slot = None);
        demo.plot_upto = 0;
        demo.forget();
    }
    if flags.remesh_requested {
        demo.forget();
    }
}

/// Lay down the next brush, or several.
///
/// Paced off `Capture::taken` when a sequence is being recorded and off the clock
/// otherwise, so the GIF is the stroke rather than a photograph of its aftermath
/// and an interactive run is watchable.
fn advance_stroke(
    time: Res<Time>,
    capture: Res<Capture>,
    flags: Res<ViewFlags>,
    demo: Option<ResMut<Demo>>,
) {
    let Some(mut demo) = demo else { return };
    if flags.paused {
        return;
    }
    let target = if capture.is_active() {
        (capture.taken as usize * demo.per_capture_frame).min(TAPE_MAX)
    } else {
        demo.stroke_timer += time.delta_secs();
        let steps = (demo.stroke_timer / STROKE_INTERVAL) as usize;
        if steps > 0 {
            demo.stroke_timer -= steps as f32 * STROKE_INTERVAL;
        }
        (demo.tape_len + steps).min(TAPE_MAX)
    };
    demo.tape_len = demo.tape_len.max(target);
}

/// The sweep: visit chunks until the frame's measurement budget is spent.
///
/// # Two things about handing Bevy sixty meshes a second
///
/// **The asset is overwritten in place**, not added and dropped. Adding a new
/// mesh and dropping the old one is what every other example here does, and at
/// one edit per click it is fine; at this rate Bevy's mesh slab allocator ends
/// up copying out of slots the asset server has already released.
///
/// **A chunk with no surface gets no asset at all**, rather than an empty one.
/// `bevy_render`'s `MeshAllocator::allocate_meshes` skips any mesh whose vertex
/// buffer is zero bytes and then copies into it unconditionally, so an empty
/// mesh logs `Use-after-free: attempted to copy element data for an unallocated
/// key` twice — once for vertices, once for indices — every time it is
/// extracted. Sixty-four empty chunks at startup is a wall of red that says
/// nothing about this demo. `Mesh3d::default()` names no asset and draws
/// nothing, which is what an empty chunk actually wants.
fn sweep(
    mut demo: Option<ResMut<Demo>>,
    mut meshes: ResMut<Assets<Mesh>>,
    mut commands: Commands,
) {
    let Some(demo) = demo.as_deref_mut() else {
        return;
    };
    let mut spent = 0.0;
    let mut visited = 0;
    while spent < MEASURE_BUDGET_MS && visited < MAX_CHUNKS_PER_FRAME {
        let index = demo.cursor;
        let id = demo.ids[index];
        let tape_len = demo.tape_len;
        let Some(m) = measure_chunk(
            &mut demo.rig,
            &demo.layout,
            demo.base,
            &demo.tape[..tape_len],
            id,
        ) else {
            // A refused chunk is logged once per sweep and skipped, not
            // unwrapped: the rest of the world still meshes.
            warn!("chunk {:?} could not be extracted", id.coords);
            demo.cursor = (index + 1) % demo.ids.len();
            return;
        };
        if !m.identical && demo.mismatch.is_none() {
            demo.mismatch = Some(id.coords);
            error!(
                "chunk {:?}: pruned and full tape disagree bit-for-bit at tape length {}",
                id.coords, tape_len
            );
        }

        let active = if demo.prune {
            &demo.rig.pruned
        } else {
            &demo.rig.full
        };
        let vertices = active.positions.len();
        let triangles = active.indices.len() / 3;
        if triangles == 0 {
            if demo.meshes[index].take().is_some() {
                commands
                    .entity(demo.entities[index])
                    .insert(Mesh3d::default());
            }
        } else if let Some(handle) = &demo.meshes[index] {
            if let Some(mut mesh) = meshes.get_mut(handle) {
                *mesh = to_bevy_mesh(active);
            }
        } else {
            let handle = meshes.add(to_bevy_mesh(active));
            commands
                .entity(demo.entities[index])
                .insert(Mesh3d(handle.clone()));
            demo.meshes[index] = Some(handle);
        }

        let record = &mut demo.records[index];
        record.absorb(m, tape_len);
        record.vertices = vertices;
        record.triangles = triangles;
        // Both arms, because both ran. Budgeting only the active one made the
        // sweep take a second chunk it could not afford and halved the frame
        // rate.
        spent += m.ms_pruned + m.ms_full;
        visited += 1;
        demo.last_visited = index;
        demo.cursor = (index + 1) % demo.ids.len();
        if demo.cursor == 0 {
            demo.sweeps += 1;
        }
    }
}

/// Tint each chunk by how much of the tape survived inside it.
///
/// Green for heavily pruned through red for nothing pruned. The spatial pattern
/// *is* the mechanism: the shell of chunks the stroke never reaches keeps one
/// brush of sixty-four, the chunks in the thick of it keep nearly all.
fn paint_heat(demo: Option<ResMut<Demo>>, mut materials: ResMut<Assets<StandardMaterial>>) {
    let Some(mut demo) = demo else { return };
    let (heat, tape_len) = (demo.heat, demo.tape_len);
    for index in 0..demo.records.len() {
        let record = demo.records[index];
        let fraction = if !heat || !record.measured || tape_len == 0 {
            f32::INFINITY
        } else {
            record.survivors as f32 / tape_len as f32
        };
        if fraction == demo.records[index].painted
            || (fraction - demo.records[index].painted).abs() < 1.0 / 128.0
        {
            continue;
        }
        let Some(mut material) = materials.get_mut(&demo.materials[index]) else {
            continue;
        };
        material.base_color = heat_colour(fraction);
        demo.records[index].painted = fraction;
    }
}

/// Survivor fraction to colour: green at nothing surviving, yellow at half, red
/// at the whole tape. A non-finite fraction means "overlay off".
///
/// Two segments rather than one lerp from green to red, because a single lerp
/// puts almost no separation in `0.0..0.3` and that is where most chunks sit —
/// the median survivor fraction is about a third. Going through yellow spends
/// half the ramp on the half of the range that is actually populated.
fn heat_colour(fraction: f32) -> Color {
    if !fraction.is_finite() {
        return Color::srgb(0.72, 0.76, 0.82);
    }
    let t = fraction.clamp(0.0, 1.0);
    if t < 0.5 {
        Color::srgb(0.16 + 1.48 * t, 0.80, 0.22)
    } else {
        Color::srgb(0.90, 0.80 - 1.16 * (t - 0.5), 0.20)
    }
}

/// Push the rolling per-chunk means into the column for the current tape length,
/// filling any columns the stroke skipped.
fn update_plot(
    demo: Option<ResMut<Demo>>,
    flags: Res<ViewFlags>,
    mut bars: Query<(&PlotBar, &mut Node)>,
    mut panels: Query<&mut Visibility, With<PlotPanel>>,
) {
    let Some(mut demo) = demo else { return };
    for mut visibility in &mut panels {
        *visibility = if flags.hud {
            Visibility::Inherited
        } else {
            Visibility::Hidden
        };
    }
    if !flags.hud {
        return;
    }

    let (full, pruned) = demo.mean_ms();
    if demo.measured() > 0 {
        let column = demo.tape_len.min(PLOT_COLS - 1);
        let from = demo.plot_upto.min(column);
        for slot in &mut demo.plot[from..=column] {
            *slot = Some((full, pruned));
        }
        demo.plot_upto = column;
    }

    for (bar, mut node) in &mut bars {
        let height = demo.plot.get(bar.column).and_then(|slot| *slot).map_or(
            0.0,
            |(full_ms, pruned_ms)| {
                let ms = if bar.full { full_ms } else { pruned_ms };
                ((ms / PLOT_MAX_MS) as f32).clamp(0.0, 1.0) * PLOT_H
            },
        );
        node.height = Val::Px(height);
    }
}

/// The HUD. The numbers are the demo.
fn report(demo: Option<Res<Demo>>, mut stats: ResMut<DemoStats>) {
    let Some(demo) = demo else { return };
    let chunks = demo.chunks();
    let measured = demo.measured();
    let (world_full, world_pruned) = demo.world_ms();
    let (mean_full, mean_pruned) = demo.mean_ms();
    let aggregate = if world_pruned > 0.0 {
        world_full / world_pruned
    } else {
        0.0
    };
    let survivors = demo.survivor_spread();
    let speedups = demo.speedups();
    let fraction = |n: usize| {
        if demo.tape_len == 0 {
            0.0
        } else {
            n as f64 / demo.tape_len as f64
        }
    };
    let bound_us = if measured == 0 {
        0.0
    } else {
        demo.records
            .iter()
            .filter(|r| r.measured)
            .map(|r| r.bound_us)
            .sum::<f64>()
            / measured as f64
    };
    let bound_share = if mean_pruned > 0.0 {
        bound_us / 1000.0 / mean_pruned
    } else {
        0.0
    };
    let agreed = demo.records.iter().filter(|r| r.identical).count();
    let cursor_id = demo.ids.get(demo.last_visited).map_or([0, 0, 0], |c| c.coords);
    let cursor_survivors = demo
        .records
        .get(demo.last_visited)
        .map_or(0, |r| r.survivors);
    // Microseconds of re-mesh cost each brush in the tape is worth, measured
    // against the base-only column of the plot rather than assumed to start at
    // zero. This is "stops tracking the edit history", as a number.
    let (slope_full, slope_pruned) = match (demo.plot.first().and_then(|c| *c), demo.tape_len) {
        (Some((base_full, base_pruned)), n) if n > 0 => (
            (mean_full - base_full) * 1000.0 / n as f64,
            (mean_pruned - base_pruned) * 1000.0 / n as f64,
        ),
        _ => (0.0, 0.0),
    };

    stats.title = format!(
        "E-305 Lipschitz tape pruning - {chunks} chunks of 33 samples, cell {CELL_SIZE}"
    );
    stats.vertices = demo.records.iter().map(|r| r.vertices).sum();
    stats.triangles = demo.records.iter().map(|r| r.triangles).sum();
    stats.extract_ms = demo
        .records
        .get(demo.last_visited)
        .map_or(0.0, |r| r.active_ms(demo.prune));

    stats.extra = vec![
        format!(
            "tape         {:>3} of {TAPE_MAX} brushes{}",
            demo.tape_len,
            if demo.tape_len == TAPE_MAX {
                "   (stroke complete)"
            } else {
                "   (sculpting)"
            }
        ),
        format!(
            "survivors    median {:>3} ({:.4})   min {:>3}   max {:>3}   cursor {cursor_id:?} {cursor_survivors}",
            median(&survivors),
            fraction(median(&survivors)),
            survivors.first().copied().unwrap_or(0),
            survivors.last().copied().unwrap_or(0),
        ),
        format!(
            "bound cost   {bound_us:>7.2} us per chunk = {bound_share:.2e} of the meshing it enables"
        ),
        String::new(),
        format!("pruned       {mean_pruned:>7.3} ms per chunk re-mesh   ({world_pruned:.1} ms for the world)"),
        format!("whole tape   {mean_full:>7.3} ms per chunk re-mesh   ({world_full:.1} ms for the world)"),
        format!(
            "speedup      {aggregate:>7.3}x world aggregate   per-chunk median {:.3}x   range {:.3}x - {:.2}x",
            median(&speedups),
            speedups.first().copied().unwrap_or(0.0),
            speedups.last().copied().unwrap_or(0.0),
        ),
        format!(
            "             rolling over the last sweep - {measured}/{chunks} chunks measured, {} sweeps",
            demo.sweeps
        ),
        format!(
            "cost slope   {slope_full:>7.1} us per brush against the whole tape, {slope_pruned:.1} pruned"
        ),
        String::new(),
        format!(
            "active arm   {}   {:.3} ms per chunk{}",
            if demo.prune { "PRUNED" } else { "WHOLE TAPE" },
            if demo.prune { mean_pruned } else { mean_full },
            if demo.prune {
                String::new()
            } else {
                format!("   (+{:.0}% - press P)", (aggregate - 1.0) * 100.0)
            }
        ),
        match demo.mismatch {
            None => format!(
                "mesh identical: YES - {agreed}/{chunks} chunks bit-exact, positions, normals and indices"
            ),
            Some(id) => format!("mesh identical: NO - chunk {id:?} disagreed; see the log"),
        },
        String::new(),
        "SmoothAdd is never pruned in the losing direction: at h == 1 smooth_min".to_string(),
        "returns b + (a - b), which is not bit-identical to a. This tape is Add/Subtract.".to_string(),
        String::new(),
        format!(
            "[P] prune {}   [H] heat {} (green 0 - yellow 0.5 - red 1.0 survivor fraction)",
            on_off(demo.prune),
            on_off(demo.heat),
        ),
        format!(
            "[C] chunk boxes {}   [X] restart stroke",
            on_off(demo.boxes)
        ),
    ];
}

/// `on`/`off`, for the HUD.
fn on_off(value: bool) -> &'static str {
    if value { "on" } else { "off" }
}

/// One box per chunk in its survivor-fraction colour, plus a bright one on the
/// chunk the sweep is about to measure.
///
/// The tint on the surface only shows chunks that *have* a surface. The boxes are
/// where the empty shell shows up, and the empty shell is where the mechanism
/// wins hardest — one surviving brush of sixty-four.
fn draw_chunk_boxes(demo: Option<Res<Demo>>, mut gizmos: Gizmos<ChunkGizmos>) {
    let Some(demo) = demo else { return };
    if !demo.boxes {
        return;
    }
    let span = demo.layout.cell_size() * demo.layout.cells() as f32;
    for (index, id) in demo.ids.iter().enumerate() {
        let origin = demo.layout.sample_origin(*id);
        let centre = Vec3::new(
            origin[0] + span * 0.5,
            origin[1] + span * 0.5,
            origin[2] + span * 0.5,
        );
        let record = demo.records[index];
        let colour = if index == demo.cursor {
            Color::srgb(0.35, 0.90, 1.0)
        } else if !record.measured || demo.tape_len == 0 {
            Color::srgb(0.28, 0.30, 0.36)
        } else {
            // The same ramp the surface tint uses, so the boxes and the solid
            // are reading off one scale rather than two that look similar.
            heat_colour(record.survivors as f32 / demo.tape_len as f32)
        };
        gizmos.cube(
            Transform::from_translation(centre).with_scale(Vec3::splat(span * 0.995)),
            colour,
        );
    }
}
