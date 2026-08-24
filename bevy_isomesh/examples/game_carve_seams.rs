//! E-312 — the rim of the shaft you just dug is lit wrong, and meshing finer
//! does not fix it.
//!
//! ```bash
//! cd bevy_isomesh && cargo run --example game_carve_seams --release
//! ```
//!
//! **Always `--release`.** The finest step of the ladder is a 129³ dual
//! contouring pass and a debug build meshes 20-50x slower.
//!
//! It runs itself and loops. `Space` freezes it; `2` `3` `4` `5` pin the four
//! beats and `1` hands it back to the loop; `W` puts the wireframe over the
//! rock, `N` draws the normals the mesh is actually carrying.
//! `ISOMESH_FIELD=1|2|3|4` pins a beat for a still without a keyboard — the
//! harness's digit keys are `flags.field + 1`, which is why the keys sit one
//! along from the variable.
//!
//! ```bash
//! FPS=10 ./scripts/record_gif.sh game_carve_seams docs/gifs/e312.gif
//! ```
//!
//! Demonstrates **M-350 / P-56** (`docs/experiments/p-56.csv`) as a player meets
//! it: not as a bound on a swept fixture, but as a ring of wrong shading on the
//! lip of a tunnel, at the exact place the player just dug.
//!
//! # What goes wrong, in one paragraph
//!
//! A carve is `max(f_rock, −f_tool)`. That composition is `C⁰` and not `C¹`
//! along the crease where the two surfaces meet, and
//! [`BrushStack`] does not override
//! [`Sdf::gradient`] — so the normal every extractor asks it for is the crate's
//! six-sample central difference. Within one differencing step of the crease
//! that stencil has samples on **both** branches, and it returns their average.
//! Exactly on the seam surface it returns `(u + v)/2` for the two branch
//! normals `u` and `v`, which is wrong by half the angle between them:
//!
//! ```text
//! angle(central difference, active branch) <= (180 - theta) / 2
//! ```
//!
//! for a crease dihedral `theta`. At a right-angled corner that is **45
//! degrees**, and this example digs a right-angled corner on purpose.
//!
//! # The clip, in four beats
//!
//! 1. **Dig.** A capsule brush pushes into the face. The crease sharpens from
//!    `180°` to `90°` and the bound with it, from `0°` to `45°`, so the ring
//!    grows out of nothing. Game lighting, no tint.
//! 2. **The artifact.** The camera has settled. The per-vertex angular error is
//!    tinted onto the surface, and it lands **exactly** on the band that was
//!    already misshaded in beat 1.
//! 3. **The ladder.** The chunk steps `25³ → 41³ → 65³ → 97³ → 129³`. Thirty
//!    times the triangles; the rim does not improve.
//! 4. **Back at the shipped LOD, then the fix.** The same 2,588 triangles beat 2
//!    showed, still just as wrong — and then the straddling vertices switch to
//!    the analytic gradient and the ring is gone.
//!
//! The last beat is a controlled A/B on purpose. Cutting from a 129³ sweep
//! straight to a repaired 25³ would change the resolution and the normals at
//! once, and a viewer could reasonably read that as "coarsening fixed it", so the
//! resolution changes first and the normals change alone. See [`FIX_AT`].
//!
//! # Why refining cannot help
//!
//! The stencil width is `Real::DIFF_STEP · max(|pₓ|, |p_y|, |p_z|, 1)` — a
//! property of the *scalar*, not of the grid. At `f32` it is `4.9216e-3`, and it
//! is the same `4.9216e-3` at 25³ and at 129³. Halving the cell moves the
//! extracted vertex *closer* to the crease, so the average gets closer to the
//! exact `(u + v)/2` and the tightness climbs toward the bound rather than away
//! from it. Measured, worst rim error against a `45.00°` bound:
//!
//! | chunk | 25³ | 41³ | 65³ | 97³ | 129³ |
//! |---|---|---|---|---|---|
//! | triangles | 2,588 | 7,492 | 19,146 | 43,516 | 78,100 |
//! | rim error | 44.28° | 42.99° | 44.16° | 44.73° | 44.94° |
//! | straddling | 2.77% | 1.60% | 1.21% | 0.92% | 0.73% |
//!
//! `41³` dips, and that is not noise to be smoothed: the number is a **maximum
//! over a finite vertex set**, so whether the bound is attained depends on
//! whether some vertex happens to land on the seam surface. Every rung is inside
//! `4.5%` of the bound and the trend is upward. Nothing here is 30x better for
//! 30x the mesh.
//!
//! The straddling column is the other half: the population is `O(n)` on an
//! `O(n²)` surface — the crease is a curve — so it halves on every doubling and
//! the fix stays cheap at any resolution.
//!
//! What refining does buy is a **thinner** wrong band — the wrong normals live
//! on one row of vertices, so their world width falls with the cell. That is the
//! honest shape of the defect and this demo says so on screen: a finer chunk
//! turns a wide misshaded rim into a hairline that is just as wrong, forever.
//!
//! # `f32`, which is the game-relevant half of the ledger
//!
//! P-56 swept this in `f64`, where the stencil is `6.0555e-6` — **813x
//! narrower**. At its coarsest row the grid had not resolved the crease to
//! within a stencil at all and the worst error reached only `0.61` of the bound.
//! At `f32`, which is what a game meshes in, the band is wide enough that the
//! bound is essentially attained at *every* resolution here, including the
//! coarsest. Different fixtures, so this is not a controlled comparison — but it
//! is the direction a shipped game runs in, and it is why this example is `f32`.
//!
//! # The geometry, and where the crease angle comes from
//!
//! The rock is a half-space: [`RockFace`] is `p·n − d`, with an analytic
//! gradient of `n` everywhere, so the rock branch contributes no error of its
//! own. The tool is one [`Capsule`] pushed along `−n` — a real
//! [`BrushStack`] `Subtract`, the same composition
//! `game_dig` carves with.
//!
//! The bore axis is **parallel to the face normal**, so on the crease the tool's
//! outward normal is radial, radial is perpendicular to `n`, and the dihedral is
//! a right angle by construction. The example does not take that on trust:
//! [`Crease::of`] bisects the tool's own SDF along the face to find the crease
//! radius and reads the dihedral off the two analytic branch gradients there,
//! which is also what makes the number correct *during* the dig — while the
//! capsule's leading cap is what cuts the face, the crease is shallower and the
//! angle sweeps `180° → 90°` as the shaft goes in. The bound sweeps `0° → 45°`
//! with it, and the ring appears out of nothing.
//!
//! # One key light, raking
//!
//! The default harness lighting is a soft three-quarter key plus a lot of
//! ambient, and it hides this. A normal error only shows up as `N·L`, so the
//! light is put **65° off the face normal**, nearly in the plane of the wall,
//! with the ambient cut to a fifth. On the wall `N·L` is then `0.42`; a rim
//! vertex whose normal has been dragged 45° toward the bore reads
//! `(0.42 ∓ 0.91)/√2` depending on which side of the hole it is — clipped to
//! black on the light side and `0.94` on the far side. So the artifact renders
//! as a dark arc and a bright arc on the same rim, more than twice the wall's
//! own brightness. That contrast is a lighting choice; the 45° is not.
//!
//! # What the fix is, and what it does not prove
//!
//! Beat 4 replaces the normal at every straddling vertex with the analytic
//! active-branch gradient — `n` on the rock side, `−∇f_tool` on the bore side.
//! The reported error there then falls to **zero by construction**, so that
//! number is not the evidence and is not offered as any. Three things are:
//!
//! - **The count.** 38 of 1,374 vertices at the shipped LOD, 290 of 39,484 at
//!   129³ — and the share halves on every doubling, because the crease is a
//!   curve on a surface.
//! - **What is left.** The worst error anywhere in the repaired mesh is
//!   `0.0013°`, which is the `f32` round-off floor the smooth sheets already sat
//!   at. Nothing else in the mesh was touched or needed to be.
//! - **The picture.** The ring is gone at the same triangle count. That is the
//!   claim a designer can act on.
//!
//! One thing the fix does **not** do: give a crease a crisp shaded edge. Dual
//! contouring puts one vertex in a crease cell and both sheets share it, so a
//! single normal per vertex cannot represent a normal discontinuity. Taking the
//! active branch makes the rock side exact and leaves the bore's first row of
//! triangles interpolating across the crease — which is what a hard edge looks
//! like in any mesh with shared normals. Splitting the vertex is a different
//! ticket; this one is about the normal being wrong on **both** sides.
//!
//! # Every number on screen is measured in this process
//!
//! The ladder is measured once at startup, at full depth, through the same
//! [`measure`] the per-frame rebuild uses, and logged. The only numbers read
//! from `docs/experiments/p-56.csv` are P-56's own, quoted for comparison and
//! never mixed into a live figure. Two things are re-checked before the window
//! opens, because a wrong fixture voids everything downstream: that [`stencil`]
//! is **bit-identical** to `BrushStack::gradient` on every vertex of a real
//! extraction, and that the crease [`Crease::of`] finds really is a right angle.
//! Either failing puts a `SELF CHECK FAILED` line on the HUD rather than
//! panicking in a stranger's terminal.
//!
//! # Angles use `atan2`, not `acos`
//!
//! The control number here is the error *away* from the rim, and it is
//! `1.4e-3` degrees. `acos` near `1` cannot resolve better than `1.2e-6`
//! degrees at `f64` and roughly `0.02` degrees at `f32` — larger than the
//! control — so [`angle_deg`] uses `2·atan2(|û − v̂|, |û + v̂|)` and widens to
//! `f64` to do it.

mod common;

use std::time::Instant;

use bevy::asset::RenderAssetUsages;
use bevy::mesh::{Indices, PrimitiveTopology};
use bevy::prelude::*;
use common::{Capture, CommonPlugin, DemoMesh, DemoStats, OrbitCamera, ViewFlags};
use isomesh::brush::{Brush, BrushStack, Capsule};
use isomesh::dual_contouring::DualContouring;
use isomesh::{MeshBuffer, Real, RuntimeShape3, Sdf};

// ─── the scene, in world units ──────────────────────────────────────────────

/// Which way the rock face points, before normalising. Leaning back a little:
/// a vertical wall lit from the side is a flat grey rectangle, and the tilt is
/// what makes the raking key light describe a surface.
const FACE_TILT: [f32; 3] = [0.0, 0.45, 1.0];

/// How far the face sits from the origin, along its own normal.
///
/// Not zero, and that is the reason it is a constant rather than a literal
/// `0.0`: a face through the origin puts an exact `f32` zero on any grid corner
/// the plane happens to pass through, and an exactly-zero corner is a coin flip
/// between inside and outside rather than a measurement.
const FACE_OFFSET: f32 = 0.017;

/// Bore radius. `0.55` puts the rim circle at 46% of the sampled window, which
/// fills the frame without the far side of it leaving the grid.
const BORE_RADIUS: f32 = 0.55;

/// How far outside the rock the tool's near end is parked, along the face
/// normal.
///
/// Larger than [`BORE_RADIUS`], so the capsule's *near* cap never reaches the
/// face and the crease is always cut by the cylinder or by the leading cap —
/// never by both at once, which would make "the crease angle" two numbers.
const TOOL_MOUTH: f32 = 0.62;

/// How deep the leading end goes at full depth, below the face.
///
/// Past the far wall of the sampled window on purpose: a shaft that stops
/// inside the grid is a pocket, and the domed end of it reads as a bubble
/// rather than as a tunnel.
const TOOL_DEPTH: f32 = 2.2;

/// Side of the sampled window, in world units. One chunk, centred on the origin.
const WINDOW: f32 = 2.4;

/// The chunk resolutions the demo measures, coarsest first, in samples per axis.
///
/// **Not settable from the environment, unlike every other resolution in this
/// repo.** The ladder *is* the finding — the point is that stepping it changes
/// the triangle count by 30x and the rim error by less than a degree — so a run
/// with one resolution pinned would be a run of a different demo.
const LADDER: [u32; 5] = [25, 41, 65, 97, 129];

// ─── the story ──────────────────────────────────────────────────────────────

/// Where the dig beat ends, as a fraction of the clip.
///
/// The four boundaries are set for `FPS=10` over the harness's default 80
/// frames: a two-second dig, a second on the artifact, half a second on each
/// rung of the ladder, and two seconds on the last beat, which is the one with
/// the payoff in it.
const DIG_END: f32 = 0.26;
/// Where the artifact beat ends.
const ARTIFACT_END: f32 = 0.40;
/// Where the resolution ladder ends and the last beat begins.
const SWEEP_END: f32 = 0.74;

/// How far into the last beat the analytic gradient is switched on.
///
/// **The last beat is a controlled A/B, and this is the only thing that changes
/// across it.** Its first third re-shows the shipped LOD exactly as the sweep
/// left it — same 2,588 triangles, same wrong rim — and then the normals flip.
/// Ending the sweep at 129³ and cutting straight to a repaired 25³ would change
/// the resolution and the normals at once, and a viewer could reasonably read
/// that as "coarsening fixed it".
const FIX_AT: f32 = 0.33;

/// Seconds for one pass through the story, when nobody is capturing.
const STORY_SECONDS: f32 = 26.0;

/// One beat of the story.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum Beat {
    /// The brush sweeps in and the shaft opens. Game lighting, no tint.
    Dig,
    /// The camera has settled on the rim and the error is tinted onto it.
    Artifact,
    /// The chunk resolution steps up the ladder and the rim does not improve.
    Sweep,
    /// Back at the shipped LOD, and part way through it the straddling vertices
    /// take the analytic gradient instead. See [`FIX_AT`].
    Fix,
}

impl Beat {
    /// The beat at `phase`, and how far through it we are.
    fn at(phase: f32) -> (Self, f32) {
        let span = |lo: f32, hi: f32| ((phase - lo) / (hi - lo)).clamp(0.0, 1.0);
        if phase < DIG_END {
            (Self::Dig, span(0.0, DIG_END))
        } else if phase < ARTIFACT_END {
            (Self::Artifact, span(DIG_END, ARTIFACT_END))
        } else if phase < SWEEP_END {
            (Self::Sweep, span(ARTIFACT_END, SWEEP_END))
        } else {
            (Self::Fix, span(SWEEP_END, 1.0))
        }
    }

    /// The state a digit or `ISOMESH_FIELD` pins, as a beat and a progress
    /// through it.
    ///
    /// Each pins the *end* of its beat, which for the last one is after the
    /// normals have been switched. `ISOMESH_FIELD=3` is therefore the finest
    /// chunk with the whole table up, and `4` is the repaired rim.
    fn pinned(index: usize) -> Option<(Self, f32)> {
        match index {
            1 => Some((Self::Dig, 1.0)),
            2 => Some((Self::Artifact, 1.0)),
            3 => Some((Self::Sweep, 1.0)),
            4 => Some((Self::Fix, 1.0)),
            _ => None,
        }
    }
}

// ─── small vector helpers ───────────────────────────────────────────────────

fn dot(u: [f32; 3], v: [f32; 3]) -> f32 {
    u[0] * v[0] + u[1] * v[1] + u[2] * v[2]
}

fn cross(u: [f32; 3], v: [f32; 3]) -> [f32; 3] {
    [
        u[1] * v[2] - u[2] * v[1],
        u[2] * v[0] - u[0] * v[2],
        u[0] * v[1] - u[1] * v[0],
    ]
}

fn scale(u: [f32; 3], k: f32) -> [f32; 3] {
    [u[0] * k, u[1] * k, u[2] * k]
}

/// `u`, scaled to unit length. Falls over loudly rather than returning a
/// direction nobody chose.
fn unit(u: [f32; 3]) -> [f32; 3] {
    let n = dot(u, u).sqrt();
    assert!(n > 0.0 && n.is_finite(), "cannot normalise {u:?}");
    scale(u, n.recip())
}

/// The angle between two directions, in degrees, accurate at both ends.
///
/// `2·atan2(|û − v̂|, |û + v̂|)`, widened to `f64` for the arithmetic. See the
/// module docs for why not `acos`.
fn angle_deg(u: [f32; 3], v: [f32; 3]) -> f64 {
    let (a, b) = (unit(u), unit(v));
    let diff = [
        f64::from(a[0] - b[0]),
        f64::from(a[1] - b[1]),
        f64::from(a[2] - b[2]),
    ];
    let sum = [
        f64::from(a[0] + b[0]),
        f64::from(a[1] + b[1]),
        f64::from(a[2] + b[2]),
    ];
    let dn = (diff[0] * diff[0] + diff[1] * diff[1] + diff[2] * diff[2]).sqrt();
    let sn = (sum[0] * sum[0] + sum[1] * sum[1] + sum[2] * sum[2]).sqrt();
    (2.0 * dn.atan2(sn)).to_degrees()
}

/// `n` with a thousands separator, for a HUD a designer reads.
fn commas(n: usize) -> String {
    let digits = n.to_string();
    let mut out = String::with_capacity(digits.len() + digits.len() / 3);
    for (i, c) in digits.chars().enumerate() {
        if i > 0 && (digits.len() - i).is_multiple_of(3) {
            out.push(',');
        }
        out.push(c);
    }
    out
}

/// Degrees at the precision the magnitude deserves.
///
/// The same field carries `44.28` and `0.0012` — the shipped rim and the
/// repaired one — and a fixed `{:.2}` prints the second as `0.00`, which reads
/// as "not measured" rather than as four orders of magnitude better.
fn deg(v: f64) -> String {
    if v >= 1.0 {
        format!("{v:.2}")
    } else {
        format!("{v:.4}")
    }
}

// ─── the field ──────────────────────────────────────────────────────────────

/// The rock, as a half-space: `p·n − d`, negative inside.
///
/// A plane rather than a heightfield, and the reason is the measurement rather
/// than the art direction: the rock branch's gradient is `n` *exactly*, at every
/// point, so every degree of error the demo reports comes from the composition
/// and none of it from a curved base field. The rock reads as rock because of
/// the light, not because of the field.
#[derive(Clone, Copy)]
struct RockFace {
    /// Unit outward normal — out of the rock, into the air.
    normal: [f32; 3],
    /// Offset along `normal`.
    offset: f32,
}

impl Sdf for RockFace {
    type Scalar = f32;

    fn sample(&self, p: [f32; 3]) -> f32 {
        dot(p, self.normal) - self.offset
    }

    fn gradient(&self, _p: [f32; 3]) -> [f32; 3] {
        self.normal
    }
}

/// The tool's outward normal at `p`, analytically.
///
/// [`Capsule`] does not override [`Sdf::gradient`] either, so this is written
/// here rather than asked for: it is the direction from the nearest point of the
/// segment, which is what the capsule's distance is the length of.
fn tool_normal(tool: &Capsule<f32>, p: [f32; 3]) -> [f32; 3] {
    let ab = [
        tool.b[0] - tool.a[0],
        tool.b[1] - tool.a[1],
        tool.b[2] - tool.a[2],
    ];
    let ap = [p[0] - tool.a[0], p[1] - tool.a[1], p[2] - tool.a[2]];
    let denom = dot(ab, ab);
    let t = if denom > 0.0 {
        (dot(ap, ab) / denom).clamp(0.0, 1.0)
    } else {
        0.0
    };
    unit([ap[0] - ab[0] * t, ap[1] - ab[1] * t, ap[2] - ab[2] * t])
}

/// The composed value at `q`, and whether the rock branch was the active one.
///
/// Derived from the two branch samples rather than read back out of
/// [`BrushStack`], and bit-identical to it: `apply(Subtract, field, shape)` *is*
/// `field.max(-shape)`. Computing both is what makes the straddle test exact
/// instead of a Lipschitz bound on the margin.
fn value(face: &RockFace, tool: &Capsule<f32>, q: [f32; 3]) -> (f32, bool) {
    let rock = face.sample(q);
    let bore = -tool.sample(q);
    if rock >= bore {
        (rock, true)
    } else {
        (bore, false)
    }
}

/// The analytic normal of the active branch at `p` — the right answer.
fn analytic_normal(face: &RockFace, tool: &Capsule<f32>, p: [f32; 3]) -> [f32; 3] {
    if value(face, tool, p).1 {
        face.normal
    } else {
        scale(tool_normal(tool, p), -1.0)
    }
}

/// The differencing step [`Sdf::gradient`] uses at `p`.
fn diff_step(p: [f32; 3]) -> f32 {
    <f32 as Real>::DIFF_STEP * p[0].abs().max(p[1].abs()).max(p[2].abs()).max(1.0)
}

/// What the six-sample stencil at `p` returned, and whether it crossed branches.
struct Stencil {
    gradient: [f32; 3],
    straddles: bool,
}

/// The crate's default central difference, re-implemented so the branch of each
/// of the six samples is visible.
///
/// The gradient it returns is **bit-identical** to `BrushStack::gradient`, and
/// [`self_check`] proves that on every vertex of a real extraction rather than
/// asserting it here: the formula below is `sdf.rs`'s, step for step, and the
/// only thing added is the `branch` array.
fn stencil(face: &RockFace, tool: &Capsule<f32>, p: [f32; 3]) -> Stencil {
    let h = diff_step(p);
    let inv = (2.0 * h).recip();
    let mut gradient = [0.0f32; 3];
    let mut branch = [false; 6];
    for axis in 0..3 {
        let mut lo = p;
        let mut hi = p;
        lo[axis] -= h;
        hi[axis] += h;
        let (v_lo, rock_lo) = value(face, tool, lo);
        let (v_hi, rock_hi) = value(face, tool, hi);
        gradient[axis] = (v_hi - v_lo) * inv;
        branch[2 * axis] = rock_lo;
        branch[2 * axis + 1] = rock_hi;
    }
    Stencil {
        gradient,
        straddles: branch.iter().any(|b| *b != branch[0]),
    }
}

// ─── the crease, measured from the geometry ─────────────────────────────────

/// Where the rock face meets the bore, and at what angle.
#[derive(Clone, Copy, Default)]
struct Crease {
    /// In-plane radius of the rim circle. Zero when the tool has not broken
    /// through yet.
    radius: f32,
    /// The dihedral through the material, in degrees. `180` is flat, `90` a
    /// right angle, `0` a knife edge.
    dihedral_deg: f64,
    /// `(180 − theta)/2` — the most a central difference can be wrong by here.
    bound_deg: f64,
}

impl Crease {
    /// Find the crease by bisecting the tool's own SDF along the rock face.
    ///
    /// The bore axis is parallel to the face normal, so the tool's distance
    /// restricted to the face depends only on the in-plane radius and is
    /// monotone in it — one root, and bisection finds it without a derivative.
    /// The dihedral is then read off the two analytic branch gradients at that
    /// point: with `n₁` the rock's outward normal and `n₂` the tool's,
    /// `cos theta = n₁·n₂`, so `theta` is the angle between them.
    ///
    /// Returns `None` before the tool has broken the surface.
    fn of(face: &RockFace, tool: &Capsule<f32>, in_plane: [f32; 3]) -> Option<Self> {
        let foot = scale(face.normal, face.offset);
        let at = |rho: f32| {
            tool.sample([
                foot[0] + rho * in_plane[0],
                foot[1] + rho * in_plane[1],
                foot[2] + rho * in_plane[2],
            ])
        };
        let mut lo = 0.0f32;
        let mut hi = BORE_RADIUS * 4.0;
        if at(lo) >= 0.0 || at(hi) <= 0.0 {
            return None;
        }
        // 40 halvings of a 2.2-wide bracket lands inside 2e-12, which is well
        // under an f32 ulp at this magnitude.
        for _ in 0..40 {
            let mid = 0.5 * (lo + hi);
            if at(mid) <= 0.0 {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        let radius = 0.5 * (lo + hi);
        let p = [
            foot[0] + radius * in_plane[0],
            foot[1] + radius * in_plane[1],
            foot[2] + radius * in_plane[2],
        ];
        let dihedral_deg = angle_deg(face.normal, tool_normal(tool, p));
        Some(Self {
            radius,
            dihedral_deg,
            bound_deg: 0.5 * (180.0 - dihedral_deg),
        })
    }
}

// ─── the scene ──────────────────────────────────────────────────────────────

/// The rock face and the frame the shot is composed in.
#[derive(Resource, Clone, Copy)]
struct Scene {
    face: RockFace,
    /// A unit vector in the plane of the face. Also the direction the crease
    /// bisection walks.
    in_plane: [f32; 3],
    /// The other one, completing a right-handed frame with the face normal.
    up_plane: [f32; 3],
}

impl Scene {
    fn new() -> Self {
        let normal = unit(FACE_TILT);
        // Gram-Schmidt x against the normal. `FACE_TILT` has no x component, so
        // this is x itself; written out anyway so the frame survives a change to
        // the tilt.
        let k = dot([1.0, 0.0, 0.0], normal);
        let in_plane = unit([1.0 - k * normal[0], -k * normal[1], -k * normal[2]]);
        Self {
            face: RockFace {
                normal,
                offset: FACE_OFFSET,
            },
            in_plane,
            up_plane: cross(normal, in_plane),
        }
    }

    /// The tool at dig progress `u`, from just touching to full depth.
    ///
    /// `u²` rather than `u`, so half the dig beat is spent on the part where the
    /// leading cap is what cuts the face and the crease angle is still sweeping
    /// `180° → 90°`. Linear would spend 80% of the beat deepening a shaft whose
    /// rim had stopped changing.
    fn tool(&self, u: f32) -> Capsule<f32> {
        let height = BORE_RADIUS - (BORE_RADIUS + TOOL_DEPTH) * u * u;
        Capsule {
            a: scale(self.face.normal, TOOL_MOUTH),
            b: scale(self.face.normal, self.face.offset + height),
            radius: BORE_RADIUS,
        }
    }
}

// ─── one measurement ────────────────────────────────────────────────────────

/// Everything one rebuild measured.
#[derive(Clone, Default)]
struct Row {
    samples: u32,
    cell: f32,
    vertices: usize,
    triangles: usize,
    extract_ms: f64,
    crease: Crease,
    /// Vertices whose stencil crossed the seam.
    straddling: usize,
    /// Worst angular error, in degrees, over the normals the mesh is carrying.
    worst_deg: f64,
    /// Mean over the straddling population.
    straddling_mean_deg: f64,
    /// Worst over everything that did not straddle — the control.
    smooth_worst_deg: f64,
    /// Whether the analytic gradient was substituted at straddling vertices.
    fixed: bool,
}

impl Row {
    fn straddling_share(&self) -> f64 {
        if self.vertices == 0 {
            0.0
        } else {
            self.straddling as f64 / self.vertices as f64
        }
    }

    fn ratio(&self) -> f64 {
        if self.crease.bound_deg > 0.0 {
            self.worst_deg / self.crease.bound_deg
        } else {
            0.0
        }
    }
}

/// One rung of the resolution ladder, kept so the viewer can watch the column
/// fill in.
#[derive(Clone, Copy)]
struct Rung {
    samples: u32,
    triangles: usize,
    worst_deg: f64,
    /// Straddling vertices as a fraction of all of them. Halves on every
    /// doubling, which is the `O(n)`-on-an-`O(n²)`-surface claim made visible.
    share: f64,
    /// `worst / bound` for this rung, kept because the fix beat re-measures the
    /// same chunk and needs the *pre-fix* ratio to quote against the ledger.
    ratio: f64,
}

// ─── colour ─────────────────────────────────────────────────────────────────

/// Unpainted rock.
const ROCK_SRGB: [f32; 4] = [0.44, 0.40, 0.35, 1.0];
/// A vertex whose normal is off by [`TINT_SCALE`] degrees or more.
const HOT_SRGB: [f32; 4] = [1.0, 0.13, 0.05, 1.0];
/// Degrees of angular error that saturate the tint.
///
/// The right-angle bound, fixed for the whole clip rather than tracking the
/// live bound: a colour scale that rescales itself is a colour scale that cannot
/// show the rim heating up as the shaft goes in.
const TINT_SCALE: f64 = 45.0;

/// sRGB as a human picks it into the linear RGBA [`Mesh::ATTRIBUTE_COLOR`]
/// wants. Feeding sRGB in raw renders it washed out (E-208).
fn linear(srgb: [f32; 4]) -> [f32; 4] {
    Color::srgba(srgb[0], srgb[1], srgb[2], srgb[3])
        .to_linear()
        .to_f32_array()
}

/// The tint for an error of `error_deg`, in linear RGBA.
fn heat(rock: [f32; 4], hot: [f32; 4], error_deg: f64) -> [f32; 4] {
    let t = (error_deg / TINT_SCALE).clamp(0.0, 1.0) as f32;
    [
        rock[0] + (hot[0] - rock[0]) * t,
        rock[1] + (hot[1] - rock[1]) * t,
        rock[2] + (hot[2] - rock[2]) * t,
        1.0,
    ]
}

// ─── framing ────────────────────────────────────────────────────────────────

/// How the camera direction is mixed out of the face frame: mostly along the
/// face normal, swung right and lifted.
///
/// 33° off the normal, measured. Face-on and the bore is a black disc with no
/// depth; much further round and the rim's far side goes edge-on and half the
/// ring the demo is about stops being visible.
const CAMERA_MIX: [f32; 3] = [0.80, 0.42, 0.30];

/// Orbit radius at the end of the push-in, in world units.
///
/// The window is `2.4` across and Bevy's default vertical FOV is 45°, so this
/// frames `2.0` of height — the rim circle takes a bit over half of it, and the
/// chunk's own boundary sits just outside the frame top and bottom.
const CAMERA_RADIUS: f32 = 2.45;

/// Radius the push-in starts from.
///
/// Only 11% wider than where it lands. Further back and the sampled chunk stops
/// covering the frame, so the opening shot is a lit quadrilateral floating in
/// the sky — which is what 3.85 photographed as.
const CAMERA_RADIUS_WIDE: f32 = 2.72;

/// How far down the face the camera aims, below the bore's own centre, in world
/// units.
///
/// Aiming straight at the bore put the chunk's lower boundary edge inside the
/// frame, and a straight cut across the bottom of a rock face reads as a
/// floating slab. This slides that edge just out of shot. Small: `0.45` was
/// tried first and clipped the top of the bore off the frame.
const FOCUS_DROP: f32 = 0.18;

/// Radians of extra yaw at the start of the push-in, unwound as the shaft goes
/// in. A moving camera on the dig and a dead-still one afterwards, so the beats
/// that ask the viewer to compare frames are frames that can be compared.
const CAMERA_SWING: f32 = 0.26;

/// Degrees the key light sits off the face normal. See the module docs — this is
/// the number that turns a 45° normal error into a black arc and a white one.
const KEY_LIGHT_RAKE: f32 = 65.0;

/// Where in the plane of the face the key light comes from: up and to the left.
const KEY_LIGHT_IN_PLANE: [f32; 2] = [-0.55, 0.84];

/// Width and height of the backdrop the HUD is read against, in logical pixels.
///
/// The wall is a mid-tone and the HUD is light text, so the numbers cross onto
/// it and lose their contrast. Sized for the widest and tallest state the HUD
/// reaches — a 72-character line, and 26 lines of them once the ladder is full
/// and the fix line is up. At the harness's 13 px font the pitch is 15.6 px.
const HUD_PANEL: Vec2 = Vec2::new(600.0, 432.0);

/// Illuminance of the key light, in lux.
///
/// Chosen so the wall itself renders mid-grey rather than white. That is not
/// taste: the artifact's bright side is `2.2x` the wall's brightness, so a wall
/// exposed anywhere near `1.0` clips the very thing the demo is pointing at.
/// Measured on a 1280x720 capture — at 26,000 lux the wall and the bright arc
/// were both white.
const KEY_LIGHT_LUX: f32 = 8_500.0;

/// Ambient fill, against the harness default of 220.
///
/// Ambient is exactly the dial that hides this: it adds light a wrongly-oriented
/// normal receives anyway, so it lifts the dark arc back toward the wall. A
/// fifth of the default keeps the bore from being a featureless hole without
/// filling the artifact in.
const AMBIENT_BRIGHTNESS: f32 = 45.0;

// ─── the ledger, compiled in ────────────────────────────────────────────────

/// P-56's committed artefact, embedded at compile time.
///
/// `include_str!` rather than transcribed constants: the path resolves against
/// this source file so no working directory can break it, and a number that
/// lived only here could drift away from the CSV with nothing to say so.
const LEDGER_CSV: &str = include_str!("../../docs/experiments/p-56.csv");

/// The dihedral this example's geometry builds, and the one row family of the
/// ledger it can honestly be compared against.
const LEDGER_DIHEDRAL: f64 = 90.0;

/// The one committed row this demo quotes: P-56's finest `f64` grid at the same
/// dihedral.
#[derive(Resource, Clone, Copy)]
struct Ledger {
    samples: u32,
    bound_deg: f64,
    worst_deg: f64,
    ratio: f64,
    share: f64,
}

impl Ledger {
    /// Pull the `theta = 90°`, finest-grid row out of the CSV by header name.
    fn load() -> Option<Self> {
        let mut lines = LEDGER_CSV.lines().filter(|l| !l.starts_with('#'));
        let header: Vec<&str> = lines.next()?.split(',').collect();
        let column = |name: &str| header.iter().position(|h| *h == name);
        let (c_dih, c_samples, c_bound, c_worst, c_ratio, c_share) = (
            column("dihedral_deg")?,
            column("samples_per_axis")?,
            column("predicted_bound_deg")?,
            column("straddling_max_error_deg")?,
            column("worst_over_bound_ratio")?,
            column("straddling_share_of_vertices")?,
        );
        let mut best: Option<Self> = None;
        for line in lines {
            let cells: Vec<&str> = line.split(',').collect();
            if cells.len() != header.len() {
                continue;
            }
            let dihedral: f64 = cells[c_dih].parse().ok()?;
            if (dihedral - LEDGER_DIHEDRAL).abs() > 1e-9 {
                continue;
            }
            let row = Self {
                samples: cells[c_samples].parse().ok()?,
                bound_deg: cells[c_bound].parse().ok()?,
                worst_deg: cells[c_worst].parse().ok()?,
                ratio: cells[c_ratio].parse().ok()?,
                share: cells[c_share].parse().ok()?,
            };
            if best.is_none_or(|b| row.samples > b.samples) {
                best = Some(row);
            }
        }
        best
    }

    /// Whether the live measurement reproduces the ledger's claim, and if not,
    /// which of the two reasons it is.
    ///
    /// Not "are the numbers equal" — they cannot be, this is a different fixture
    /// in a different scalar. The claim is that the bound is the same and that
    /// both essentially attain it, so that is what is checked. A mid-dig frame
    /// fails the first half honestly — the crease is not a right angle yet — and
    /// says so rather than reporting a disagreement it was never going to win.
    ///
    /// `shipped_ratio` is the tightness of the mesh **as extracted**, which is
    /// not the same as the live row's during the fix beat: the fix drives the
    /// error to zero by construction, and quoting that against a ledger row
    /// about un-repaired normals would compare two different things.
    fn verdict(&self, bound_deg: f64, shipped_ratio: f64) -> &'static str {
        if (self.bound_deg - bound_deg).abs() >= 0.01 {
            "crease is not 90 deg yet"
        } else if self.ratio > 0.9 && shipped_ratio > 0.9 {
            "agrees"
        } else {
            "no vertex near the corner yet"
        }
    }
}

// ─── resources ──────────────────────────────────────────────────────────────

/// What this frame is showing.
#[derive(Resource, Clone, Copy, PartialEq)]
struct Shot {
    beat: Beat,
    /// Progress through the beat.
    local: f32,
    /// Dig progress, `0` just touching and `1` full depth.
    depth: f32,
    samples: u32,
    fixed: bool,
    tinted: bool,
    /// How many rungs of the ladder the HUD has uncovered. The table is measured
    /// up front; this is only how much of it the viewer has been shown.
    revealed: usize,
}

impl Default for Shot {
    fn default() -> Self {
        Self {
            beat: Beat::Dig,
            local: 0.0,
            depth: 0.0,
            samples: LADDER[0],
            fixed: false,
            tinted: false,
            revealed: 0,
        }
    }
}

/// What the current chunk measured.
#[derive(Resource, Default)]
struct Live {
    row: Row,
    /// The `Shot` the mesh in the asset was built from, so a frame that changes
    /// nothing costs nothing.
    built: Option<Shot>,
}

/// The resolution ladder, measured once at startup.
#[derive(Resource)]
struct Rungs(Vec<Rung>);

/// The extractor and its buffer, kept across frames so the ladder does not
/// allocate a fresh 39,000-vertex buffer per rung.
#[derive(Resource)]
struct Rig {
    dc: DualContouring<f32>,
    buffer: MeshBuffer<f32>,
    /// The entity the rock is drawn on.
    surface: Entity,
    /// The asset the last rebuild wrote, dropped when the next one replaces it.
    mesh: Option<Handle<Mesh>>,
    rock: [f32; 4],
    hot: [f32; 4],
}

/// The startup self-check, so a wrong fixture is visible on screen rather than
/// only in the log.
#[derive(Resource)]
struct SelfCheck {
    /// Vertices where the re-implemented stencil disagreed with
    /// `BrushStack::gradient`. Must be zero.
    gradient_mismatches: usize,
    /// How far the measured right angle sat from 90°, in degrees.
    dihedral_residual_deg: f64,
}

/// The bottom caption — the line a viewer reads instead of the HUD.
#[derive(Component)]
struct Caption;

// ─── app ────────────────────────────────────────────────────────────────────

fn main() {
    App::new()
        .add_plugins(DefaultPlugins.set(WindowPlugin {
            primary_window: Some(Window {
                title: "isomesh - E-312 carve seams".into(),
                ..default()
            }),
            ..default()
        }))
        .add_plugins(CommonPlugin)
        // Not black. A single lit slab against a void reads as a floating
        // polygon; against a dim sky it reads as a rock face, which is what the
        // demo is claiming to be footage of. Dark enough that the HUD stays
        // legible over it.
        .insert_resource(ClearColor(Color::srgb(0.11, 0.13, 0.17)))
        .insert_resource(Scene::new())
        .insert_resource(Shot::default())
        .init_resource::<Live>()
        .add_systems(Startup, setup)
        // `PreUpdate`, not `Update`, and it is a correctness fix rather than a
        // preference. The harness's `update_hud` lives in `Update` and system
        // order within a schedule is unspecified, so it read `DemoStats` from
        // the *previous* frame while the caption -- written here -- was current.
        // Photographed: the HUD said the crease was 150.6 degrees while the
        // caption under it said 135. Two numbers on screen disagreeing is worse
        // than either being late.
        .add_systems(PreUpdate, (advance, rebuild, frame_camera, report).chain())
        .run();
}

fn setup(
    mut commands: Commands,
    mut materials: ResMut<Assets<StandardMaterial>>,
    mut meshes: ResMut<Assets<Mesh>>,
    mut lights: Query<(&mut DirectionalLight, &mut Transform)>,
    mut ambient: ResMut<GlobalAmbientLight>,
    mut camera: Query<&mut OrbitCamera>,
    scene: Res<Scene>,
) {
    // One key light, raking. See the module docs for the 65 degrees; the point
    // is that a normal error is only ever visible as a change in `N·L`, and the
    // harness default lights this scene almost flat.
    let (rake_sin, rake_cos) = KEY_LIGHT_RAKE.to_radians().sin_cos();
    let in_plane = unit([
        KEY_LIGHT_IN_PLANE[0] * scene.in_plane[0] + KEY_LIGHT_IN_PLANE[1] * scene.up_plane[0],
        KEY_LIGHT_IN_PLANE[0] * scene.in_plane[1] + KEY_LIGHT_IN_PLANE[1] * scene.up_plane[1],
        KEY_LIGHT_IN_PLANE[0] * scene.in_plane[2] + KEY_LIGHT_IN_PLANE[1] * scene.up_plane[2],
    ]);
    let toward_light = unit([
        rake_cos * scene.face.normal[0] + rake_sin * in_plane[0],
        rake_cos * scene.face.normal[1] + rake_sin * in_plane[1],
        rake_cos * scene.face.normal[2] + rake_sin * in_plane[2],
    ]);
    for (mut light, mut transform) in &mut lights {
        light.illuminance = KEY_LIGHT_LUX;
        // Off. A light this close to the plane of the wall throws shadow acne
        // across the whole face at any bias that does not also detach the
        // shaft's own shadow -- and a speckled wall beside a real shading
        // artifact is a demo that cannot be read.
        light.shadow_maps_enabled = false;
        *transform = Transform::from_translation(Vec3::from_array(scale(toward_light, 12.0)))
            .looking_at(Vec3::ZERO, Vec3::Y);
    }
    ambient.brightness = AMBIENT_BRIGHTNESS;

    let direction = unit([
        CAMERA_MIX[0] * scene.face.normal[0]
            + CAMERA_MIX[1] * scene.in_plane[0]
            + CAMERA_MIX[2] * scene.up_plane[0],
        CAMERA_MIX[0] * scene.face.normal[1]
            + CAMERA_MIX[1] * scene.in_plane[1]
            + CAMERA_MIX[2] * scene.up_plane[1],
        CAMERA_MIX[0] * scene.face.normal[2]
            + CAMERA_MIX[1] * scene.in_plane[2]
            + CAMERA_MIX[2] * scene.up_plane[2],
    ]);
    for mut orbit in &mut camera {
        orbit.focus = Vec3::from_array(scale(scene.up_plane, -FOCUS_DROP));
        orbit.radius = CAMERA_RADIUS_WIDE;
        orbit.yaw = direction[2].atan2(direction[0]);
        orbit.pitch = direction[1].asin();
    }

    // White base colour, because the rock's colour arrives per vertex: the same
    // attribute carries the heat tint, so there is one path onto a vertex rather
    // than a base colour and an override.
    let material = materials.add(StandardMaterial {
        base_color: Color::WHITE,
        perceptual_roughness: 0.94,
        metallic: 0.0,
        ..default()
    });
    // A placeholder asset, so the surface entity is complete from frame zero
    // and `rebuild` only ever has to swap a handle.
    let handle = meshes.add(Mesh::new(
        PrimitiveTopology::TriangleList,
        RenderAssetUsages::default(),
    ));
    let surface = commands
        .spawn((Mesh3d(handle.clone()), MeshMaterial3d(material), DemoMesh))
        .id();

    // Behind the harness HUD, which is spawned by `CommonPlugin` at the default
    // z and left alone here. `GlobalZIndex(-1)` is the whole mechanism: no
    // reaching into the shared module, and the panel disappears with the HUD
    // only in the sense that it is empty dark pixels when `nohud` clears it.
    commands.spawn((
        Node {
            position_type: PositionType::Absolute,
            top: Val::Px(6.0),
            left: Val::Px(6.0),
            width: Val::Px(HUD_PANEL.x),
            height: Val::Px(HUD_PANEL.y),
            ..default()
        },
        BackgroundColor(Color::srgba(0.04, 0.05, 0.07, 0.58)),
        GlobalZIndex(-1),
    ));

    commands
        .spawn((
            Node {
                position_type: PositionType::Absolute,
                left: Val::Px(0.0),
                right: Val::Px(0.0),
                bottom: Val::Px(20.0),
                justify_content: JustifyContent::Center,
                ..default()
            },
            GlobalZIndex(4),
        ))
        .with_children(|parent| {
            parent.spawn((
                Text::new(""),
                TextFont {
                    font_size: FontSize::Px(20.0),
                    ..default()
                },
                // `NoWrap`: in a centring flex row the measure is handed the
                // container's whole width but the node's height resolves before
                // the wrap, so a soft wrap pushes the second line off frame.
                TextLayout {
                    linebreak: bevy::text::LineBreak::NoWrap,
                    ..default()
                },
                TextColor(Color::srgb(0.97, 0.94, 0.90)),
                BackgroundColor(Color::srgba(0.03, 0.03, 0.05, 0.82)),
                Node {
                    padding: UiRect::axes(Val::Px(14.0), Val::Px(8.0)),
                    ..default()
                },
                Caption,
            ));
        });

    let ledger = Ledger::load().expect("p-56.csv carries a theta=90 row");
    info!(
        "p-56 @ theta={LEDGER_DIHEDRAL}: bound {:.3} deg, worst {:.3} deg at {}^3 (f64), ratio {:.4}, straddling share {:.4}%",
        ledger.bound_deg,
        ledger.worst_deg,
        ledger.samples,
        ledger.ratio,
        ledger.share * 100.0
    );
    info!(
        "stencil width: f32 {:e}, f64 {:e} -- {:.0}x wider at f32, and neither shrinks with the cell",
        <f32 as Real>::DIFF_STEP,
        <f64 as Real>::DIFF_STEP,
        f64::from(<f32 as Real>::DIFF_STEP) / <f64 as Real>::DIFF_STEP
    );

    let mut rig = Rig {
        dc: DualContouring::<f32>::new(),
        buffer: MeshBuffer::<f32>::new(),
        surface,
        mesh: Some(handle),
        rock: linear(ROCK_SRGB),
        hot: linear(HOT_SRGB),
    };
    commands.insert_resource(self_check(&scene, &mut rig));
    commands.insert_resource(Rungs(measure_ladder(&scene, &mut rig)));
    commands.insert_resource(rig);
    commands.insert_resource(ledger);
}

/// Re-measure the fixture before the window opens.
///
/// Two things that void every number downstream if they are wrong, and they are
/// checked rather than asserted in a comment:
///
/// - **The stencil is the crate's.** Every vertex of a real extraction is
///   differenced twice, once through `BrushStack::gradient` and once through
///   [`stencil`], and the two `[f32; 3]`s are compared for **bit** equality. A
///   non-zero count would mean the branch flags belong to a different gradient
///   than the one being graded.
/// - **The crease really is a right angle.** [`Crease::of`] bisects for it; this
///   holds the answer against the 90° the construction claims.
///
/// It logs and reports rather than panicking: a demo a stranger runs is not the
/// place for an assertion, and the HUD says so if either number is wrong.
fn self_check(scene: &Scene, rig: &mut Rig) -> SelfCheck {
    let tool = scene.tool(1.0);
    let brushes = [Brush::subtract(tool)];
    let field = BrushStack {
        base: scene.face,
        brushes: &brushes,
    };
    let samples = LADDER[1];
    let cell = WINDOW / (samples - 1) as f32;
    let shape = RuntimeShape3::new([samples; 3]).expect("window grid fits u32");
    let Rig { dc, buffer, .. } = rig;
    buffer.reset();
    if let Err(e) = dc.extract(&field, &shape, [-WINDOW * 0.5; 3], cell, buffer) {
        error!("self check could not mesh the shaft: {e}");
    }
    let mut gradient_mismatches = 0usize;
    for &p in &buffer.positions {
        if field.gradient(p) != stencil(&scene.face, &tool, p).gradient {
            gradient_mismatches += 1;
        }
    }
    let crease = Crease::of(&scene.face, &tool, scene.in_plane).unwrap_or_default();
    let residual = (crease.dihedral_deg - 90.0).abs();
    if gradient_mismatches == 0 {
        info!(
            "self check: stencil is bit-identical to BrushStack::gradient on all {} vertices",
            buffer.positions.len()
        );
    } else {
        error!(
            "self check: {gradient_mismatches} of {} vertices disagree with BrushStack::gradient",
            buffer.positions.len()
        );
    }
    info!(
        "self check: crease radius {:.6}, dihedral {:.6} deg (residual {:.2e}), bound {:.6} deg",
        crease.radius, crease.dihedral_deg, residual, crease.bound_deg
    );
    SelfCheck {
        gradient_mismatches,
        dihedral_residual_deg: residual,
    }
}

/// Decide what this frame is about.
///
/// Under capture the story advances with the captured frame count, so a clip of
/// any length is the whole story rather than a slice of it. Interactively it
/// runs on wall-clock time and loops, and a digit pins one beat.
fn advance(
    time: Res<Time>,
    capture: Res<Capture>,
    flags: Res<ViewFlags>,
    mut shot: ResMut<Shot>,
    mut elapsed: Local<f32>,
) {
    let (beat, local) = if let Some(pinned) = Beat::pinned(flags.field) {
        pinned
    } else {
        let phase = if capture.is_active() {
            f32::from(u16::try_from(capture.taken).unwrap_or(u16::MAX))
                / f32::from(u16::try_from(capture_frames()).unwrap_or(1).max(1))
        } else {
            if !flags.paused {
                *elapsed += time.delta_secs();
            }
            (*elapsed / STORY_SECONDS).fract()
        };
        Beat::at(phase.clamp(0.0, 1.0))
    };

    // `min` rather than a wrap: `local` reaches exactly 1.0 on the last frame
    // of the beat and that frame belongs to the last rung.
    let rung = ((local * LADDER.len() as f32) as usize).min(LADDER.len() - 1);
    *shot = Shot {
        beat,
        local,
        depth: if beat == Beat::Dig { local } else { 1.0 },
        samples: match beat {
            Beat::Sweep => LADDER[rung],
            // The shipped LOD, for the dig, the artifact and the fix alike.
            _ => LADDER[0],
        },
        fixed: beat == Beat::Fix && local >= FIX_AT,
        tinted: beat != Beat::Dig,
        revealed: match beat {
            // The ladder is the sweep's device; before it there is nothing to
            // compare and the HUD stays four lines long.
            Beat::Dig | Beat::Artifact => 0,
            Beat::Sweep => rung + 1,
            Beat::Fix => LADDER.len(),
        },
    };
}

/// `ISOMESH_CAPTURE_FRAMES`, or the harness default.
fn capture_frames() -> u32 {
    std::env::var("ISOMESH_CAPTURE_FRAMES")
        .ok()
        .and_then(|v| v.parse().ok())
        .filter(|n| *n > 0)
        .unwrap_or(60)
}

/// Extract one chunk and grade every vertex it produced.
///
/// The one measurement path in the example: the startup ladder pass and every
/// frame's rebuild both come through here, so a rung in the table and the row on
/// the HUD cannot be produced by two different pieces of arithmetic. Leaves the
/// mesh in `rig.buffer` and the per-vertex tint in `colors`.
fn measure(
    scene: &Scene,
    tool: &Capsule<f32>,
    samples: u32,
    fixed: bool,
    tinted: bool,
    rig: &mut Rig,
    colors: &mut Vec<[f32; 4]>,
) -> Row {
    let brushes = [Brush::subtract(*tool)];
    let field = BrushStack {
        base: scene.face,
        brushes: &brushes,
    };
    let cell = WINDOW / (samples - 1) as f32;
    let shape = RuntimeShape3::new([samples; 3]).expect("window grid fits u32");

    let Rig {
        dc,
        buffer,
        rock,
        hot,
        ..
    } = rig;
    let (rock, hot) = (*rock, *hot);
    buffer.reset();
    colors.clear();
    let started = Instant::now();
    if let Err(e) = dc.extract(&field, &shape, [-WINDOW * 0.5; 3], cell, buffer) {
        error!("dual contouring failed on the shaft at {samples}^3: {e}");
        return Row::default();
    }
    let extract_ms = started.elapsed().as_secs_f64() * 1000.0;

    let mut straddling = 0usize;
    let mut straddling_sum = 0.0f64;
    let mut worst = 0.0f64;
    let mut smooth_worst = 0.0f64;
    colors.reserve(buffer.positions.len());
    for i in 0..buffer.positions.len() {
        let p = buffer.positions[i];
        let straddles = stencil(&scene.face, tool, p).straddles;
        let analytic = analytic_normal(&scene.face, tool, p);
        if fixed && straddles {
            buffer.normals[i] = analytic;
        }
        // Graded against the normal that is *in the mesh*, so one expression
        // covers both the shipped state and the repaired one: no second code
        // path, and the tint cannot disagree with the number on the HUD.
        let error = angle_deg(buffer.normals[i], analytic);
        if straddles {
            straddling += 1;
            straddling_sum += error;
        } else {
            smooth_worst = smooth_worst.max(error);
        }
        worst = worst.max(error);
        colors.push(if tinted { heat(rock, hot, error) } else { rock });
    }

    Row {
        samples,
        cell,
        vertices: buffer.positions.len(),
        triangles: buffer.indices.len() / 3,
        extract_ms,
        crease: Crease::of(&scene.face, tool, scene.in_plane).unwrap_or_default(),
        straddling,
        worst_deg: worst,
        straddling_mean_deg: if straddling > 0 {
            straddling_sum / straddling as f64
        } else {
            0.0
        },
        smooth_worst_deg: smooth_worst,
        fixed,
    }
}

/// Measure every rung of the ladder once, at full depth, before the window
/// opens.
///
/// **The table is a property of the geometry, not of what the clip happened to
/// visit.** Built by walking as the sweep ran, `ISOMESH_FIELD=3` produced a
/// still with one row in it — the pin never passes through the other four — and
/// the row it produced was whichever resolution the pin landed on rather than
/// the shipped LOD. Measuring up front costs 210 ms of startup and makes the
/// pinned still carry the whole argument.
fn measure_ladder(scene: &Scene, rig: &mut Rig) -> Vec<Rung> {
    let tool = scene.tool(1.0);
    let mut scratch: Vec<[f32; 4]> = Vec::new();
    LADDER
        .iter()
        .map(|&samples| {
            let row = measure(scene, &tool, samples, false, false, rig, &mut scratch);
            info!(
                "ladder {samples}^3: {} tris, rim {:.4} deg of {:.4}, {} of {} straddle ({:.4}%), {:.2} ms",
                row.triangles,
                row.worst_deg,
                row.crease.bound_deg,
                row.straddling,
                row.vertices,
                row.straddling_share() * 100.0,
                row.extract_ms
            );
            Rung {
                samples,
                triangles: row.triangles,
                worst_deg: row.worst_deg,
                share: row.straddling_share(),
                ratio: row.ratio(),
            }
        })
        .collect()
}

/// Mesh the chunk the shot asks for — only when the answer would change.
fn rebuild(
    scene: Res<Scene>,
    shot: Res<Shot>,
    mut live: ResMut<Live>,
    mut rig: ResMut<Rig>,
    mut meshes: ResMut<Assets<Mesh>>,
    mut commands: Commands,
) {
    // The dig moves the tool every frame; the later beats do not move anything,
    // so they extract once and then cost nothing at all. Every beat change that
    // alters what is on screen alters one of these four fields, so the beat
    // itself is not part of the key.
    let stale = live.built.is_none_or(|built| {
        built.samples != shot.samples
            || built.fixed != shot.fixed
            || built.tinted != shot.tinted
            || (shot.beat == Beat::Dig && built.depth != shot.depth)
    });
    if !stale {
        return;
    }

    let tool = scene.tool(shot.depth);
    let mut colors: Vec<[f32; 4]> = Vec::new();
    let row = measure(
        &scene,
        &tool,
        shot.samples,
        shot.fixed,
        shot.tinted,
        &mut rig,
        &mut colors,
    );

    let mut mesh = Mesh::new(
        PrimitiveTopology::TriangleList,
        RenderAssetUsages::default(),
    );
    mesh.insert_attribute(Mesh::ATTRIBUTE_POSITION, rig.buffer.positions.clone());
    mesh.insert_attribute(Mesh::ATTRIBUTE_NORMAL, rig.buffer.normals.clone());
    mesh.insert_attribute(Mesh::ATTRIBUTE_COLOR, colors);
    mesh.insert_indices(Indices::U32(rig.buffer.indices.clone()));
    // A fresh asset and a swapped handle, so the old buffers are released
    // rather than left in `Assets` behind a handle nothing points at.
    let handle = meshes.add(mesh);
    if let Some(old) = rig.mesh.replace(handle.clone()) {
        meshes.remove(&old);
    }
    commands.entity(rig.surface).insert(Mesh3d(handle));

    live.row = row;
    live.built = Some(*shot);
}

/// Push in while the shaft is being dug, then hold dead still.
fn frame_camera(shot: Res<Shot>, mut camera: Query<&mut OrbitCamera>) {
    if shot.beat != Beat::Dig {
        return;
    }
    // Smoothstep, so the push-in has no visible start or stop.
    let t = shot.local * shot.local * (3.0 - 2.0 * shot.local);
    for mut orbit in &mut camera {
        orbit.radius = CAMERA_RADIUS_WIDE + (CAMERA_RADIUS - CAMERA_RADIUS_WIDE) * t;
        orbit.yaw += CAMERA_SWING * (1.0 - t) * 0.02;
    }
}

// ─── what is on screen ──────────────────────────────────────────────────────

fn report(
    live: Res<Live>,
    shot: Res<Shot>,
    ledger: Res<Ledger>,
    rungs: Res<Rungs>,
    check: Res<SelfCheck>,
    mut stats: ResMut<DemoStats>,
    mut caption: Query<&mut Text, With<Caption>>,
) {
    let row = &live.row;
    stats.title = String::from("E-312 game carve seams -- M-350 / P-56, in a tunnel");
    stats.vertices = row.vertices;
    stats.triangles = row.triangles;
    stats.extract_ms = row.extract_ms;

    let mut extra = vec![
        format!(
            "crease             {:>6.1} deg  {}",
            row.crease.dihedral_deg,
            if row.crease.dihedral_deg > 179.0 {
                "the tool has barely broken the face"
            } else if row.crease.dihedral_deg < 90.5 {
                "the wall meets the bore at a right angle"
            } else {
                "sharpening as the shaft goes in"
            }
        ),
        format!(
            "worst normal error {:>7} deg  against a {:.2} deg bound (180 - theta)/2",
            deg(row.worst_deg),
            row.crease.bound_deg
        ),
        format!(
            "away from the rim  {:>7} deg  over {} vertices that do not straddle",
            deg(row.smooth_worst_deg),
            commas(row.vertices - row.straddling)
        ),
        format!(
            "chunk              {:>4}^3      cell {:.4}, {} triangles",
            row.samples,
            row.cell,
            commas(row.triangles)
        ),
    ];

    let shown = &rungs.0[..shot.revealed.min(rungs.0.len())];
    if !shown.is_empty() {
        extra.push(String::from(
            "\nresolution ladder -- the triangles climb, the rim does not improve",
        ));
        for (i, rung) in shown.iter().enumerate() {
            extra.push(format!(
                "  {:>4}^3  {:>8} tris   rim {:>6.2} deg   {:>5.2}% straddle{}",
                rung.samples,
                commas(rung.triangles),
                rung.worst_deg,
                rung.share * 100.0,
                if i == 0 { "  <- shipped LOD" } else { "" }
            ));
        }
        if let (Some(first), Some(last)) = (shown.first(), shown.last())
            && shown.len() > 1
        {
            extra.push(format!(
                "  {:.1}x the triangles, {:+.2} deg on the rim -- refining did not help",
                last.triangles as f64 / first.triangles as f64,
                last.worst_deg - first.worst_deg
            ));
        }
    }

    if row.fixed {
        extra.push(format!(
            "\nfix  analytic gradient on {} of {} vertices ({:.2}%) -- rim now {} deg",
            commas(row.straddling),
            commas(row.vertices),
            row.straddling_share() * 100.0,
            deg(row.worst_deg)
        ));
    } else {
        extra.push(format!(
            "\nrim  {} of {} vertices straddle the seam ({:.2}%), mean error {:.1} deg",
            commas(row.straddling),
            commas(row.vertices),
            row.straddling_share() * 100.0,
            row.straddling_mean_deg
        ));
    }

    // Quoted against the mesh as extracted. During the fix beat that is the
    // ladder's rung for this resolution, not the live row -- see
    // `Ledger::verdict`.
    let shipped_ratio = if row.fixed {
        rungs
            .0
            .iter()
            .find(|r| r.samples == row.samples)
            .map_or_else(|| row.ratio(), |r| r.ratio)
    } else {
        row.ratio()
    };
    extra.push(format!(
        "p-56 f64 @ {:.0} deg  {:.2} of {:.2} deg at {}^3, ratio {:.3}",
        LEDGER_DIHEDRAL, ledger.worst_deg, ledger.bound_deg, ledger.samples, ledger.ratio
    ));
    extra.push(format!(
        "     this run  {:.2} of {:.2} deg at {}^3, ratio {:.3} -- {}",
        shipped_ratio * row.crease.bound_deg,
        row.crease.bound_deg,
        row.samples,
        shipped_ratio,
        ledger.verdict(row.crease.bound_deg, shipped_ratio)
    ));

    if check.gradient_mismatches > 0 {
        extra.push(format!(
            "SELF CHECK FAILED: {} vertices disagree with BrushStack::gradient",
            check.gradient_mismatches
        ));
    }
    if check.dihedral_residual_deg > 1e-3 {
        extra.push(format!(
            "SELF CHECK FAILED: crease dihedral is {:.4} deg off the right angle",
            check.dihedral_residual_deg
        ));
    }
    stats.extra = extra;

    let text = caption_for(&live, &shot, &ledger);
    for mut target in &mut caption {
        target.0.clone_from(&text);
    }
}

/// The line a viewer reads instead of the HUD.
fn caption_for(live: &Live, shot: &Shot, ledger: &Ledger) -> String {
    let row = &live.row;
    match shot.beat {
        Beat::Dig => {
            if row.crease.bound_deg < 1.0 {
                String::from("carving a shaft into the rock face")
            } else {
                format!(
                    "carving -- the rim is a {:.0} deg crease now,\nso its normals can be {:.0} deg wrong",
                    row.crease.dihedral_deg, row.crease.bound_deg
                )
            }
        }
        Beat::Artifact => format!(
            "the rim of your own tunnel is lit wrong: {:.1} deg off,\nand {:.0} deg is the most a right-angled crease can be",
            row.worst_deg, row.crease.bound_deg
        ),
        Beat::Sweep => format!(
            "{} triangles, and the rim is still {:.1} deg off\na finer cell makes the wrong band thinner, never less wrong",
            commas(row.triangles),
            row.worst_deg
        ),
        // Two captions, because this beat is a before and an after. The "before"
        // has to name the return to the shipped LOD or the cut back from 129³
        // reads as a mistake.
        Beat::Fix if !row.fixed => format!(
            "back at {} triangles, and the rim is exactly as wrong: {:.1} deg\nthirty times the mesh had bought nothing",
            commas(row.triangles),
            row.worst_deg
        ),
        Beat::Fix => format!(
            "same {} triangles, analytic gradient on {} of {} vertices ({:.1}%)\nthe rim is fixed -- p-56 says the same in f64 at {}^3: {:.2} of {:.2} deg",
            commas(row.triangles),
            commas(row.straddling),
            commas(row.vertices),
            row.straddling_share() * 100.0,
            ledger.samples,
            ledger.worst_deg,
            ledger.bound_deg
        ),
    }
}
