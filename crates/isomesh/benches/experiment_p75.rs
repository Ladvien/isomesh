//! **P-75 — material weights carried at the vertex, against the edit log walked
//! at the fragment.**
//!
//! Ticket: R-075. Pre-registered before this harness existed.
//!
//! ```bash
//! cargo bench --bench experiment_p75
//! ```
//!
//! Writes `docs/experiments/p-75.csv`.
//!
//! # The SHARE line, recomputed before the harness was written
//!
//! **C1's denominator is `extract_ms`** — one full extraction of `game_dig`'s
//! 256 chunks over the painted world — and its numerator is a stage that is
//! **currently zero**: `game_dig` carries no material attribute at all, its
//! shader derives the three terrain strata per fragment from the interpolated
//! normal and the world height (`bevy_isomesh/examples/triplanar.wgsl:119-126`),
//! and `isomesh::paint::shade` exists but nothing in the example calls it. So
//! the clause is *"a new cost, as a fraction of extraction"* and the arithmetic
//! that bounds it is
//!
//! ```text
//! weights_share ≈ (vertices / samples) · (per-vertex walk / per-sample extraction step)
//! ```
//!
//! `game_dig`'s world is 8×4×8 chunks of 16³ cells, so **1,257,728 samples**
//! (256 · 17³) against a surface of order 4·10⁴ vertices: `vertices / samples`
//! is about **4%**. The second factor is *above* one rather than below it, and
//! that is the reason the clause is not free: `PaintStack::sample` skips sprays
//! **without evaluating their shapes** (`paint.rs:228-241`), while a material
//! weight has to evaluate every one of them. So the predicted share is a few
//! percent and lands within a factor of ~1.5 of the registered 5% — **reachable
//! in both directions, and close enough to the threshold that the measurement
//! decides it.** That is stated here, before the run, per `✗51`'s rule.
//!
//! **C2's ratio is reachable with a wide margin and its second half is nearly
//! implied.** The per-pixel vertex path is twelve fused multiply-adds, a
//! reciprocal and a four-way maximum, and it is *by construction* independent of
//! log length; the per-pixel fragment path is one walk of the whole log and is
//! monotone in it. So `speedup = fragment / vertex` cannot fail to widen unless
//! the fragment path's cost is **flat** in log length — which is exactly the
//! falsifier the registration names ("the fragment path is not paying M-138's
//! cost and the premise is wrong"). The widening clause therefore has content
//! only as a **check on the instrument**, and this harness reports
//! `fragment_ns_per_log_entry` so a reader can see the walk is the cost rather
//! than a constant with a walk-shaped label on it.
//!
//! **C3's 2% is reachable and tight.** Ten sprays land inside the sandbox
//! (below), each a 0.75-radius ball whose boundary ring on the surface is about
//! 4.7 units long, so the *"within one cell"* band is roughly
//! `10 · 4.7 · 0.125 ≈ 5.9` unit² against a surface of order 280 unit² — about
//! **2%** of the surface is within a cell of a material boundary at all. A
//! misclassification rate of 2% would therefore require misclassifying
//! essentially the whole band, and a rate above 2% would require
//! misclassification off the band. Both directions are live.
//!
//! # There is no renderer here, and this is what the instrument actually is
//!
//! `crates/isomesh` must not depend on Bevy, so **C2's numbers are per-pixel
//! shading cost over a real visibility buffer, not engine frame time.** The
//! harness software-rasterises the extracted mesh at 1920×1080 through
//! `game_dig`'s own camera — eye `(0, 1.70, 6.0)` (`game_dig.rs:947`), yaw 0 and
//! pitch −0.15 (`game_dig.rs:698-701`, applied at `game_dig.rs:2450`), Bevy's
//! default 45° vertical field of view — keeping the nearest fragment per pixel
//! with its triangle index and its perspective-correct barycentrics. That is a
//! visibility buffer in the ordinary sense, and the two shading paths are then
//! run over the *same* buffer in the *same* process:
//!
//! - **the varying interpolation both paths pay identically** — world position
//!   and normal — is done once, outside both timers, and reported as
//!   `interp_ms`. On hardware that is the interpolator; charging it to one side
//!   would be charging the comparison.
//! - **the vertex path** interpolates the four-channel weight attribute from the
//!   triangle's three vertices, renormalises, and takes the argmax.
//! - **the fragment path** recomputes the same quantity exactly: the stratum
//!   blend from the interpolated normal and height, then **one walk of the edit
//!   log**, compositing every spray's coverage.
//!
//! `frame_ms_vertex` and `frame_ms_fragment` are milliseconds for one such pass
//! over the whole 1920×1080 buffer. **They are shading cost. They are not frame
//! time and no frame time is invented here.** The rasterisation itself is
//! reported separately (`raster_ms`) and charged to neither: it is the geometry
//! pass both paths share.
//!
//! The five sandbox walls are **not** in the mesh, and that is deliberate rather
//! than an omission: they are Bevy cuboids outside the field (`game_dig.rs:805-841`)
//! carrying the forced concrete layer (`LAYER_CONCRETE`), so they walk no log in
//! either path and adding them would only add pixels the two paths shade
//! identically — which flatters the vertex path.
//!
//! # The scene, restated with its source named
//!
//! Every constant below is `game_dig`'s, quoted with its line, because this
//! bench cannot link `bevy_isomesh`:
//!
//! | quantity | value | source |
//! |---|---|---|
//! | chunk edge | 16 cells | `game_dig.rs:125` |
//! | cell size | 0.125 | `game_dig.rs:127` |
//! | chunks | 8×4×8 | `game_dig.rs:137` |
//! | layout origin | `(−8, −5.4, −8)` | `game_dig.rs:955-969` |
//! | terrain | `y − (0.35·sin(0.9x)·cos(0.7z) + 0.15·sin(2.1x))` | `game_dig.rs:171-175` |
//! | brush radius | 0.25 | `game_dig.rs:1056` |
//! | carve `n` | `(−0.9 + 0.30n, 0.55 − 0.045n, 2.2 − 0.34n)` | `game_dig.rs:665-668` |
//! | stratum blend | `up·shallow, (1−up)·shallow, 1−shallow` | `triplanar.wgsl:119-126` |
//! | slope ramp | `smoothstep(0.55, 0.82, n.y)` | `triplanar.wgsl:76-77, 123` |
//! | depth ramp | `smoothstep(−1.6, −0.4, p.y)` | `triplanar.wgsl:83-84, 124` |
//!
//! # Why the log must contain sprays, and what that says about M-50's trace
//!
//! `M-138`'s cost is `PaintStack`'s, and `PaintStack` is the only thing in this
//! crate for which "the material at this point" is *not* a pure function of
//! position: a splat's coverage carries the factor `ramp(|f_prefix(p)|, depth)`,
//! the field **as it stood at that splat's place in the log**
//! (`paint.rs:182-217`). Strip the sprays and the material is `triplanar.wgsl`'s
//! stratum blend, which needs no walk at all — and then C2's fragment path has
//! nothing to pay and the registered comparison does not exist. So the log here
//! is `M-50`'s 60-carve `ISOMESH_AUTOCARVE=60` trace **with one spray after
//! every third carve**, at the carve's centre, radius 0.75, `softness = 0.10`
//! and `depth = 0.05` — the two widths from `paint.rs`'s own module example
//! (`paint.rs:50-55`) — cycling the four terrain layers so that all four
//! materials appear and every pair can share a boundary. This is a deviation and
//! it is the only one: it is recorded as `sprays` and `log_len` on every row.
//!
//! **The trace has a defect and it is `M-50`'s, not this harness's.** Carve
//! centres leave the sandbox at `n = 30` (`x = 8.1 > 8`, `z = −8.0`) and never
//! come back, and the first four are in the air above the hill. So of the 60
//! brushes `M-50` binned, **26 ever touch solid rock inside the world** — the
//! rest lengthen the log that every sample walks and change no geometry
//! whatsoever. That is reported as `brushes_biting`, and it is very likely the
//! mechanism behind `M-50`'s own "3.7× for 7× the log, and flattening".
//!
//! # Controls, and each one could have failed
//!
//! - **`boundary_vertices` — the registered vacuity control.** Counted over the
//!   **visible** mesh only: a vertex belongs to a triangle that produced at
//!   least one pixel in the visibility buffer, and some mesh edge at that vertex
//!   changes the argmax material. Asserted non-zero, and `boundary_vertices_all`
//!   is the same count without the visibility restriction so a reader can see
//!   the restriction bite. **What it is honestly worth:** the count moves with
//!   the fixture — 294 at fifteen brushes against 526 at thirty — and a variant
//!   with every spray removed still reported 244, because `game_dig`'s stratum
//!   blend already puts a grass/dirt boundary on the tunnel walls. So this
//!   control is *satisfied* rather than *at risk*: what would zero it is a
//!   camera pointed away from the dig, not an absence of materials.
//! - **`materials`** — distinct argmax materials over the visible pixels, by the
//!   **exact** classifier. Asserted at least two, which is the other half of the
//!   registered control. All four are present.
//! - **`field_walk_mismatches`** — the weight walk carries the field through
//!   every carve itself, and that carried value is compared **bit-for-bit**
//!   against `PaintStack::sample` at every vertex. It is the same operation
//!   sequence, so zero is the only correct answer; a non-zero would mean the
//!   walk being timed is not `M-138`'s walk. Asserted zero — and the zero was
//!   **verified by mutation** rather than trusted: swapping `brush::apply`'s two
//!   value arguments in the timed walk gives `left: 24470, right: 0`, every
//!   vertex. `M-44`: a zero that could not have been non-zero is not a
//!   measurement.
//! - **`l1_max_error`** — the weight vector is L1-normalised *by construction*
//!   (`up·shallow + (1−up)·shallow + (1−shallow) = 1` identically, and a splat
//!   composite `w ← w(1−c) + e·c` preserves the sum), so this measures
//!   accumulated rounding and nothing else. Asserted under 1e-4; it reads 6e-8
//!   to 9e-8, so it is a live number rather than a structural zero.
//! - **`clipped_triangles`** — triangles discarded whole because a vertex fell
//!   behind the near plane. The camera sits 1.2 units above the hilltop, so the
//!   only candidates are directly beneath and behind it; the count is on the row
//!   rather than assumed away.
//!
//! # The variant that was built and put back
//!
//! One spray-free arm — the same 60 carves, no `Edit::Spray` at all — was run to
//! find out what the vacuity control is worth and what the weight stage costs
//! without the sprays to evaluate. It reported `weights_share = 0.0203` against
//! `vertices / samples = 25,482 / 1,257,728 = 0.02026`, which is the SHARE
//! arithmetic landing on the nose: with sprays skipped, a per-vertex weight
//! costs **exactly one sample's worth of walk**, so the share *is* the
//! vertex-to-sample ratio and nothing else. It is not committed as an arm
//! because C2's premise requires the sprays; it owes an `E×` row.

#![allow(
    clippy::cast_possible_truncation,
    clippy::cast_precision_loss,
    clippy::cast_sign_loss,
    clippy::too_many_lines
)]

mod common;

use std::time::Instant;

use isomesh::brush::{self, Brush};
use isomesh::chunk::{ChunkId, ChunkLayout};
use isomesh::fields::Sphere;
use isomesh::marching_cubes::MarchingCubes;
use isomesh::paint::{Edit, PaintStack, Splat};
use isomesh::{MeshBuffer, Sdf};

// ── game_dig's scene, restated ──────────────────────────────────────────────

/// Chunk edge in cells. `game_dig.rs:125`.
const CHUNK_CELLS: u32 = 16;
/// Cell size. `game_dig.rs:127`. A power of two, so the seam is bit-exact.
const CELL_SIZE: f32 = 0.125;
/// Chunks along x, y, z. `game_dig.rs:137`: a 16×8×16-unit sandbox of 256 chunks.
const EXTENT: [i32; 3] = [8, 4, 8];
/// Layout origin. `game_dig.rs:955-969`.
const ORIGIN: [f32; 3] = [
    -(EXTENT[0] as f32) * CHUNK_CELLS as f32 * CELL_SIZE * 0.5,
    -5.4,
    -(EXTENT[2] as f32) * CHUNK_CELLS as f32 * CELL_SIZE * 0.5,
];
/// Brush radius. `game_dig.rs:1056`.
const BRUSH_RADIUS: f32 = 0.25;

/// The terrain before any edit. `game_dig.rs:165-176`, transcribed.
#[derive(Clone, Copy, Debug)]
struct Ground;

impl Sdf for Ground {
    type Scalar = f32;

    fn sample(&self, p: [f32; 3]) -> f32 {
        let height = 0.35 * (p[0] * 0.9).sin() * (p[2] * 0.7).cos() + 0.15 * (p[0] * 2.1).sin();
        p[1] - height
    }
}

/// Where the `n`th scripted carve goes. `game_dig.rs:664-668`.
fn carve_centre(n: u32) -> [f32; 3] {
    let t = n as f32;
    [-0.9 + t * 0.30, 0.55 - t * 0.045, 2.2 - t * 0.34]
}

// ── the stratum blend, from triplanar.wgsl ──────────────────────────────────

/// Number of terrain layers, and the width of the weight vector.
/// `game_dig.rs:575`: grass, surface dirt, deep dirt, concrete.
const MATERIALS: usize = 4;

/// `triplanar.wgsl:76-77`.
const GRASS_SLOPE_LO: f32 = 0.55;
/// `triplanar.wgsl:76-77`.
const GRASS_SLOPE_HI: f32 = 0.82;
/// `triplanar.wgsl:83-84`.
const SHALLOW_Y_LO: f32 = -1.6;
/// `triplanar.wgsl:83-84`.
const SHALLOW_Y_HI: f32 = -0.4;

/// WGSL's `smoothstep`.
fn smoothstep(lo: f32, hi: f32, x: f32) -> f32 {
    let t = ((x - lo) / (hi - lo)).clamp(0.0, 1.0);
    t * t * (3.0 - 2.0 * t)
}

/// The stratum weights before any paint. `triplanar.wgsl:119-126`.
///
/// Sums to one identically, for every `up` and every `shallow` — the shader says
/// so in a comment and it is why nothing here needs renormalising.
fn base_weights(world_y: f32, normal_y: f32) -> [f32; MATERIALS] {
    let up = smoothstep(GRASS_SLOPE_LO, GRASS_SLOPE_HI, normal_y);
    let shallow = smoothstep(SHALLOW_Y_LO, SHALLOW_Y_HI, world_y);
    [up * shallow, (1.0 - up) * shallow, 1.0 - shallow, 0.0]
}

// ── the log ─────────────────────────────────────────────────────────────────

/// Carves per spray: one spray after every third carve.
const SPRAY_EVERY: u32 = 3;
/// Spray ball radius, three brush radii.
const SPRAY_RADIUS: f32 = 0.75;
/// Falloff outside the spray shape. `paint.rs:50-55`.
const SPRAY_SOFTNESS: f32 = 0.10;
/// Reach from the spray-time surface. `paint.rs:50-55`.
const SPRAY_DEPTH: f32 = 0.05;

/// The edit type this world's log is made of.
type WorldEdit = Edit<Sphere<f32>, Sphere<f32>, f32>;

/// `M-50`'s 60-carve trace, with a spray after every third carve.
///
/// One log, and the spray's material is a function of its ordinal rather than a
/// parallel array, so there is exactly one place the log is defined.
fn build_log(brushes: u32) -> Vec<WorldEdit> {
    let mut log = Vec::new();
    for n in 0..brushes {
        let centre = carve_centre(n);
        log.push(Edit::Carve(Brush::subtract(Sphere {
            center: centre,
            radius: BRUSH_RADIUS,
        })));
        if n % SPRAY_EVERY == 2 {
            log.push(Edit::Spray(Splat {
                shape: Sphere {
                    center: centre,
                    radius: SPRAY_RADIUS,
                },
                // RGB is unused here: the material index is the spray's ordinal
                // mod four, and alpha is coverage exactly as `paint.rs` uses it.
                color: [0.0, 0.0, 0.0, 1.0],
                softness: SPRAY_SOFTNESS,
                depth: SPRAY_DEPTH,
            }));
        }
    }
    log
}

/// Which layer the `k`th spray paints. Cycles all four, so every pair of
/// materials can share a boundary somewhere on the surface.
const fn spray_material(k: u32) -> usize {
    (k % MATERIALS as u32) as usize
}

/// A linear ramp from 1 at `x <= 0` to 0 at `x >= width`. `paint.rs:161-174`.
fn ramp(x: f32, width: f32) -> f32 {
    if width <= 0.0 {
        if x <= 0.0 { 1.0 } else { 0.0 }
    } else {
        (1.0 - x / width).clamp(0.0, 1.0)
    }
}

/// The exact material weights at `p`, from **one walk of the log**.
///
/// This is `PaintStack::color_at` (`paint.rs:199-217`) with the three-channel
/// colour accumulator replaced by a four-channel material accumulator, and it is
/// the quantity C2's fragment path and C3's reference both evaluate. The carried
/// field is returned too, so the caller can check it against
/// `PaintStack::sample` — the control that says this really is `M-138`'s walk.
///
/// L1 is preserved: `base_weights` sums to one and each composite is
/// `w ← w(1 − c) + e_m·c`.
fn exact_weights(log: &[WorldEdit], p: [f32; 3], normal_y: f32) -> ([f32; MATERIALS], f32) {
    let mut field = Ground.sample(p);
    let mut w = base_weights(p[1], normal_y);
    let mut spray = 0u32;
    for edit in log {
        match edit {
            Edit::Carve(b) => field = brush::apply(b.op, field, b.shape.sample(p)),
            Edit::Spray(s) => {
                let material = spray_material(spray);
                spray += 1;
                let coverage =
                    ramp(s.shape.sample(p), s.softness) * ramp(field.abs(), s.depth) * s.color[3];
                for (i, channel) in w.iter_mut().enumerate() {
                    let target = if i == material { 1.0 } else { 0.0 };
                    *channel += (target - *channel) * coverage;
                }
            }
        }
    }
    (w, field)
}

/// Index of the largest weight. Ties go to the lower index, deterministically.
fn argmax(w: &[f32; MATERIALS]) -> usize {
    let mut best = 0;
    for i in 1..MATERIALS {
        if w[i] > w[best] {
            best = i;
        }
    }
    best
}

// ── the camera and the visibility buffer ────────────────────────────────────

/// Registered resolution.
const SCREEN_W: usize = 1920;
/// Registered resolution.
const SCREEN_H: usize = 1080;
/// Bevy's `PerspectiveProjection::default().near`.
const NEAR: f32 = 0.1;

/// `game_dig`'s camera: `game_dig.rs:947` for the eye, `698-701` for the initial
/// look, `2450` for how the two Euler angles become a rotation.
#[derive(Clone, Copy, Debug)]
struct Camera {
    eye: [f32; 3],
    right: [f32; 3],
    up: [f32; 3],
    fwd: [f32; 3],
    tan_half: f32,
    aspect: f32,
}

impl Camera {
    /// The initial state of `game_dig`'s first-person camera.
    fn game_dig() -> Self {
        // `Quat::from_euler(EulerRot::YXZ, yaw, pitch, 0)` with yaw 0 is a
        // rotation about x by `pitch`, so `-Z` becomes `(0, sin p, -cos p)` and
        // `+Y` becomes `(0, cos p, sin p)`.
        let pitch: f32 = -0.15;
        Self {
            eye: [0.0, 1.70, 6.0],
            right: [1.0, 0.0, 0.0],
            up: [0.0, pitch.cos(), pitch.sin()],
            fwd: [0.0, pitch.sin(), -pitch.cos()],
            // Bevy's `PerspectiveProjection::default().fov` is PI / 4, vertical.
            tan_half: (std::f32::consts::FRAC_PI_8).tan(),
            aspect: SCREEN_W as f32 / SCREEN_H as f32,
        }
    }

    /// Pixel coordinates and view depth, or `None` behind the near plane.
    fn project(&self, p: [f32; 3]) -> Option<(f32, f32, f32)> {
        let v = [p[0] - self.eye[0], p[1] - self.eye[1], p[2] - self.eye[2]];
        let d = dot(v, self.fwd);
        if d < NEAR {
            return None;
        }
        let ndc_x = dot(v, self.right) / (d * self.tan_half * self.aspect);
        let ndc_y = dot(v, self.up) / (d * self.tan_half);
        Some((
            (ndc_x + 1.0) * 0.5 * SCREEN_W as f32,
            (1.0 - ndc_y) * 0.5 * SCREEN_H as f32,
            d,
        ))
    }
}

fn dot(a: [f32; 3], b: [f32; 3]) -> f32 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

/// One surviving fragment: which triangle, and where on it.
///
/// Barycentrics are already perspective-corrected, so they are the weights a
/// hardware interpolator would use.
#[derive(Clone, Copy, Debug)]
struct Fragment {
    tri: u32,
    bary: [f32; 3],
}

/// The whole extracted world, in world space, plus its vertex attributes.
struct World {
    positions: Vec<[f32; 3]>,
    normals: Vec<[f32; 3]>,
    indices: Vec<u32>,
    weights: Vec<[f32; MATERIALS]>,
}

/// Software-rasterise the world into a visibility buffer, and compact it.
///
/// Returns the surviving fragments, a per-triangle "produced at least one pixel"
/// mask, how many triangles were discarded at the near plane, and the
/// rasterisation time.
fn rasterise(world: &World, camera: &Camera) -> (Vec<Fragment>, Vec<bool>, u64, f64) {
    let pixels = SCREEN_W * SCREEN_H;
    let mut depth = vec![0.0f32; pixels];
    let mut tri_of = vec![u32::MAX; pixels];
    let mut bary_of = vec![[0.0f32; 3]; pixels];
    let mut clipped = 0u64;

    let started = Instant::now();
    let tris = world.indices.len() / 3;
    for t in 0..tris {
        let i = [
            world.indices[t * 3] as usize,
            world.indices[t * 3 + 1] as usize,
            world.indices[t * 3 + 2] as usize,
        ];
        let Some(a) = camera.project(world.positions[i[0]]) else {
            clipped += 1;
            continue;
        };
        let Some(b) = camera.project(world.positions[i[1]]) else {
            clipped += 1;
            continue;
        };
        let Some(c) = camera.project(world.positions[i[2]]) else {
            clipped += 1;
            continue;
        };
        let area = (b.0 - a.0) * (c.1 - a.1) - (c.0 - a.0) * (b.1 - a.1);
        if area == 0.0 {
            continue;
        }
        let inv_area = 1.0 / area;
        let lo_x = a.0.min(b.0).min(c.0).floor().max(0.0) as usize;
        let hi_x = (a.0.max(b.0).max(c.0).ceil()).min(SCREEN_W as f32) as usize;
        let lo_y = a.1.min(b.1).min(c.1).floor().max(0.0) as usize;
        let hi_y = (a.1.max(b.1).max(c.1).ceil()).min(SCREEN_H as f32) as usize;
        if lo_x >= hi_x || lo_y >= hi_y {
            continue;
        }
        let inv_d = [1.0 / a.2, 1.0 / b.2, 1.0 / c.2];
        for py in lo_y..hi_y {
            let fy = py as f32 + 0.5;
            let row = py * SCREEN_W;
            for px in lo_x..hi_x {
                let fx = px as f32 + 0.5;
                // Screen-space barycentrics from the three edge functions.
                let w0 = ((c.0 - b.0) * (fy - b.1) - (fx - b.0) * (c.1 - b.1)) * inv_area;
                let w1 = ((a.0 - c.0) * (fy - c.1) - (fx - c.0) * (a.1 - c.1)) * inv_area;
                let w2 = 1.0 - w0 - w1;
                if w0 < 0.0 || w1 < 0.0 || w2 < 0.0 {
                    continue;
                }
                // Perspective correction: interpolate 1/d linearly in screen
                // space, then divide.
                let recip = w0 * inv_d[0] + w1 * inv_d[1] + w2 * inv_d[2];
                let idx = row + px;
                if recip <= depth[idx] {
                    continue;
                }
                depth[idx] = recip;
                tri_of[idx] = t as u32;
                let s = 1.0 / recip;
                bary_of[idx] = [
                    w0 * inv_d[0] * s,
                    w1 * inv_d[1] * s,
                    w2 * inv_d[2] * s,
                ];
            }
        }
    }
    let raster_ms = started.elapsed().as_secs_f64() * 1e3;

    let mut frags = Vec::new();
    let mut seen = vec![false; tris];
    for idx in 0..pixels {
        let t = tri_of[idx];
        if t == u32::MAX {
            continue;
        }
        frags.push(Fragment {
            tri: t,
            bary: bary_of[idx],
        });
        seen[t as usize] = true;
    }
    (frags, seen, clipped, raster_ms)
}

/// World position and normal-y per fragment: the varyings a hardware
/// interpolator produces, which **both** shading paths receive identically.
fn interpolate(world: &World, frags: &[Fragment]) -> (Vec<[f32; 3]>, Vec<f32>, f64) {
    let started = Instant::now();
    let mut pos = Vec::with_capacity(frags.len());
    let mut ny = Vec::with_capacity(frags.len());
    for f in frags {
        let t = f.tri as usize;
        let i = [
            world.indices[t * 3] as usize,
            world.indices[t * 3 + 1] as usize,
            world.indices[t * 3 + 2] as usize,
        ];
        let mut p = [0.0f32; 3];
        let mut n = [0.0f32; 3];
        for (&b, &vi) in f.bary.iter().zip(i.iter()) {
            let vp = world.positions[vi];
            let vn = world.normals[vi];
            for axis in 0..3 {
                p[axis] += b * vp[axis];
                n[axis] += b * vn[axis];
            }
        }
        let len = (n[0] * n[0] + n[1] * n[1] + n[2] * n[2]).sqrt();
        pos.push(p);
        ny.push(if len > 0.0 { n[1] / len } else { 0.0 });
    }
    (pos, ny, started.elapsed().as_secs_f64() * 1e3)
}

/// The vertex path: interpolate the four-channel attribute, renormalise, argmax.
fn shade_vertex(world: &World, frags: &[Fragment], out: &mut [u8]) {
    for (f, slot) in frags.iter().zip(out.iter_mut()) {
        let t = f.tri as usize;
        let i0 = world.indices[t * 3] as usize;
        let i1 = world.indices[t * 3 + 1] as usize;
        let i2 = world.indices[t * 3 + 2] as usize;
        let (w0, w1, w2) = (
            &world.weights[i0],
            &world.weights[i1],
            &world.weights[i2],
        );
        let mut w = [0.0f32; MATERIALS];
        let mut sum = 0.0f32;
        for m in 0..MATERIALS {
            let v = f.bary[0] * w0[m] + f.bary[1] * w1[m] + f.bary[2] * w2[m];
            w[m] = v;
            sum += v;
        }
        let inv = 1.0 / sum;
        for v in &mut w {
            *v *= inv;
        }
        *slot = argmax(&w) as u8;
    }
}

/// The fragment path: the stratum blend, then one walk of the whole edit log.
fn shade_fragment(log: &[WorldEdit], pos: &[[f32; 3]], ny: &[f32], out: &mut [u8]) {
    for ((&p, &n), slot) in pos.iter().zip(ny.iter()).zip(out.iter_mut()) {
        let (w, _) = exact_weights(log, p, n);
        *slot = argmax(&w) as u8;
    }
}

// ── sampling the surface, for C3 ────────────────────────────────────────────

/// Surface points C3 evaluates. The registered 10^6.
const SURFACE_POINTS: usize = 1_000_000;

/// `splitmix64`, so the point set is reproducible without a dependency.
struct Rng(u64);

impl Rng {
    fn next_u64(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9e37_79b9_7f4a_7c15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
        z ^ (z >> 31)
    }

    fn next_f64(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 * (1.0 / (1u64 << 53) as f64)
    }
}

// ── the run ─────────────────────────────────────────────────────────────────

/// Whole-world extractions per bucket. Three, because `M-337`'s re-audit found a
/// registered ratio that moved 30% between runs on a governed core, and one
/// extraction of 256 chunks over an 80-entry log is a quarter of a second.
const REPS: usize = 3;

/// Shading passes per bucket, per path.
///
/// Nine rather than three, and the first run is why: this machine was **not
/// quiet** (twenty-five sibling agents building), and `frame_ms_vertex` — a loop
/// whose work is bit-identical across buckets 2, 3 and 4 because the geometry is
/// frozen — came out 10.24, 10.27 and 6.09 ms. `raster_ms` and `interp_ms` moved
/// with it, so the excursion is machine-wide rather than anything in the loops.
/// Nine reps plus a per-rep ratio is what survives that; the spread is on the
/// row as `speedup_min` and `speedup_max` so a reader can see how much it moved.
const SHADE_REPS: usize = 9;

/// The clock, on the row. `M-280`: on a governed CPU a nanosecond is not a unit
/// unless the reader can see what the core was doing.
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

/// Everything one log bucket produced.
struct Bucket {
    brushes: u32,
    sprays: u32,
    log_len: usize,
    brushes_biting: u32,
    samples: u64,
    vertices: usize,
    triangles: usize,
    visible_triangles: usize,
    extract_ms: f64,
    weights_ms: f64,
    raster_ms: f64,
    interp_ms: f64,
    frame_ms_vertex: f64,
    frame_ms_fragment: f64,
    speedup_min: f64,
    speedup_max: f64,
    visible_pixels: usize,
    pixel_disagreement: u64,
    materials: usize,
    boundary_vertices: usize,
    boundary_vertices_all: usize,
    misclass: Misclassification,
    field_walk_mismatches: u64,
    l1_max_error: f64,
    clipped_triangles: u64,
}

/// Extract the whole world over `log`, timing extraction and the weight stage
/// separately, and return the accumulated mesh with its weights.
///
/// The two timers bracket only their own work: the accumulation into world-space
/// arrays and the control checks happen outside both.
fn extract_world(log: &[WorldEdit]) -> (World, f64, f64, u64, u64, f64) {
    let layout = ChunkLayout::<f32>::new(CHUNK_CELLS, CELL_SIZE, ORIGIN).expect("valid layout");
    let shape = layout.sample_shape().expect("valid sample shape");
    let stack = PaintStack {
        base: Ground,
        edits: log,
        background: [0.5, 0.5, 0.5, 1.0],
    };
    let mut mc = MarchingCubes::<f32>::new();
    let mut buf = MeshBuffer::<f32>::with_capacity(4096, 8192);
    let mut chunk_weights: Vec<[f32; MATERIALS]> = Vec::with_capacity(4096);

    let mut world = World {
        positions: Vec::new(),
        normals: Vec::new(),
        indices: Vec::new(),
        weights: Vec::new(),
    };
    let mut extract_ms = 0.0;
    let mut weights_ms = 0.0;
    let mut samples = 0u64;
    let mut mismatches = 0u64;
    let mut l1_max = 0.0f64;

    for z in 0..EXTENT[2] {
        for y in 0..EXTENT[1] {
            for x in 0..EXTENT[0] {
                let id = ChunkId::new([x, y, z]);
                let origin = layout.sample_origin(id);
                buf.reset();

                let t0 = Instant::now();
                mc.extract(&stack, &shape, origin, CELL_SIZE, &mut buf)
                    .expect("the sandbox grid is large enough to march");
                extract_ms += t0.elapsed().as_secs_f64() * 1e3;
                samples += u64::from(CHUNK_CELLS + 1).pow(3);

                // The vertex-attribute stage: one walk of the log per vertex,
                // taking the extractor's own gradient normal for the stratum
                // blend, exactly as a shader would from the varying.
                let t1 = Instant::now();
                chunk_weights.clear();
                for (p, n) in buf.positions.iter().zip(buf.normals.iter()) {
                    let (w, _) = exact_weights(log, *p, n[1]);
                    chunk_weights.push(w);
                }
                weights_ms += t1.elapsed().as_secs_f64() * 1e3;

                // Controls, outside both timers.
                for (k, p) in buf.positions.iter().enumerate() {
                    let (_, carried) = exact_weights(log, *p, buf.normals[k][1]);
                    if carried.to_bits() != stack.sample(*p).to_bits() {
                        mismatches += 1;
                    }
                    let sum: f64 = chunk_weights[k].iter().map(|&v| f64::from(v)).sum();
                    l1_max = l1_max.max((sum - 1.0).abs());
                }

                let base = world.positions.len() as u32;
                world.positions.extend_from_slice(&buf.positions);
                world.normals.extend_from_slice(&buf.normals);
                world.weights.extend_from_slice(&chunk_weights);
                world.indices.extend(buf.indices.iter().map(|i| i + base));
            }
        }
    }
    (world, extract_ms, weights_ms, samples, mismatches, l1_max)
}

/// How many of the log's carves actually touch solid rock inside the sandbox.
///
/// The control that says `M-50`'s trace is not what it looks like: a brush whose
/// centre has left the world, or is in the air above the hill, lengthens the log
/// every sample walks and changes nothing.
fn brushes_biting(brushes: u32) -> u32 {
    let hi = [
        ORIGIN[0] + EXTENT[0] as f32 * CHUNK_CELLS as f32 * CELL_SIZE,
        ORIGIN[1] + EXTENT[1] as f32 * CHUNK_CELLS as f32 * CELL_SIZE,
        ORIGIN[2] + EXTENT[2] as f32 * CHUNK_CELLS as f32 * CELL_SIZE,
    ];
    (0..brushes)
        .filter(|&n| {
            let c = carve_centre(n);
            let inside = (0..3).all(|a| c[a] >= ORIGIN[a] && c[a] <= hi[a]);
            inside && Ground.sample(c) < BRUSH_RADIUS
        })
        .count() as u32
}

fn run_bucket(brushes: u32) -> Bucket {
    let log = build_log(brushes);
    let sprays = log
        .iter()
        .filter(|e| matches!(e, Edit::Spray(_)))
        .count() as u32;

    // C1. Three reps; the reported pair comes from the rep whose extraction was
    // the median, so the two numbers on the row are from one run (`M-281`).
    let mut reps: Vec<(f64, f64)> = Vec::with_capacity(REPS);
    let mut kept: Option<(World, u64, u64, f64)> = None;
    for _ in 0..REPS {
        let (world, extract_ms, weights_ms, samples, mismatches, l1_max) = extract_world(&log);
        reps.push((extract_ms, weights_ms));
        if kept.is_none() {
            kept = Some((world, samples, mismatches, l1_max));
        }
    }
    let (world, samples, field_walk_mismatches, l1_max_error) = kept.expect("REPS >= 1");
    let mut order: Vec<usize> = (0..reps.len()).collect();
    order.sort_by(|&a, &b| reps[a].0.total_cmp(&reps[b].0));
    let (extract_ms, weights_ms) = reps[order[reps.len() / 2]];

    // C2. One visibility buffer, two shading passes over it.
    let camera = Camera::game_dig();
    let (frags, seen_tri, clipped_triangles, raster_ms) = rasterise(&world, &camera);
    let (pos, ny, interp_ms) = interpolate(&world, &frags);
    let mut class_vertex = vec![0u8; frags.len()];
    let mut class_fragment = vec![0u8; frags.len()];

    // The two paths run back to back inside one rep, and the reported triple is
    // the rep whose **ratio** is the median. `M-281`: the ratio is what the
    // clause is about, so the two absolute times beside it have to come from one
    // run — otherwise a boost-clock excursion that hit only one of them lands in
    // the speedup. `speedup_min`/`speedup_max` put the spread on the row rather
    // than hiding it behind a median.
    let mut shade_reps: Vec<(f64, f64)> = Vec::with_capacity(SHADE_REPS);
    for _ in 0..SHADE_REPS {
        let t = Instant::now();
        shade_vertex(&world, &frags, &mut class_vertex);
        let v = t.elapsed().as_secs_f64() * 1e3;
        let t = Instant::now();
        shade_fragment(&log, &pos, &ny, &mut class_fragment);
        shade_reps.push((v, t.elapsed().as_secs_f64() * 1e3));
    }
    let mut shade_order: Vec<usize> = (0..shade_reps.len()).collect();
    shade_order.sort_by(|&a, &b| {
        (shade_reps[a].1 / shade_reps[a].0).total_cmp(&(shade_reps[b].1 / shade_reps[b].0))
    });
    let (frame_ms_vertex, frame_ms_fragment) = shade_reps[shade_order[shade_reps.len() / 2]];
    let ratios: Vec<f64> = shade_reps.iter().map(|(v, f)| f / v).collect();
    let speedup_min = ratios.iter().copied().fold(f64::INFINITY, f64::min);
    let speedup_max = ratios.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    let pixel_disagreement = class_vertex
        .iter()
        .zip(class_fragment.iter())
        .filter(|(a, b)| a != b)
        .count() as u64;

    // The vacuity control: materials and material boundaries on the VISIBLE
    // surface. `class_fragment` is the exact classifier, so `materials` is
    // ground truth rather than the approximation's opinion of itself.
    let mut present = [false; MATERIALS];
    for &m in &class_fragment {
        present[m as usize] = true;
    }
    let materials = present.iter().filter(|p| **p).count();

    let mut vertex_class = vec![u8::MAX; world.positions.len()];
    for (v, slot) in vertex_class.iter_mut().enumerate() {
        *slot = argmax(&world.weights[v]) as u8;
    }
    // Two counts, and the pair is the point: the registered control is about a
    // boundary crossing the **visible** surface, so a count over the whole mesh
    // beside it shows how much of the world's material boundary the camera
    // actually sees. If the two were equal the restriction would be doing no
    // work and the control would be weaker than it reads.
    let mut is_boundary = vec![false; world.positions.len()];
    let mut is_boundary_all = vec![false; world.positions.len()];
    for (tri, &visible) in world.indices.as_chunks::<3>().0.iter().zip(seen_tri.iter()) {
        let i = [tri[0] as usize, tri[1] as usize, tri[2] as usize];
        for e in 0..3 {
            let (a, b) = (i[e], i[(e + 1) % 3]);
            if vertex_class[a] != vertex_class[b] {
                is_boundary_all[a] = true;
                is_boundary_all[b] = true;
                if visible {
                    is_boundary[a] = true;
                    is_boundary[b] = true;
                }
            }
        }
    }
    let boundary_vertices = is_boundary.iter().filter(|b| **b).count();
    let boundary_vertices_all = is_boundary_all.iter().filter(|b| **b).count();

    // C3. 10^6 area-weighted surface points, against the field's own answer.
    let stack = PaintStack {
        base: Ground,
        edits: &log,
        background: [0.5, 0.5, 0.5, 1.0],
    };
    let misclass = misclassification(&world, &log, &stack);

    Bucket {
        brushes,
        sprays,
        log_len: log.len(),
        brushes_biting: brushes_biting(brushes),
        samples,
        vertices: world.positions.len(),
        triangles: world.indices.len() / 3,
        visible_triangles: seen_tri.iter().filter(|s| **s).count(),
        extract_ms,
        weights_ms,
        raster_ms,
        interp_ms,
        frame_ms_vertex,
        frame_ms_fragment,
        speedup_min,
        speedup_max,
        visible_pixels: frags.len(),
        pixel_disagreement,
        materials,
        boundary_vertices,
        boundary_vertices_all,
        misclass,
        field_walk_mismatches,
        l1_max_error,
        clipped_triangles,
    }
}

/// What one bucket's C3 sweep found.
struct Misclassification {
    points: u64,
    within_one_cell: u64,
    within_one_cell_axis: u64,
    paint: u64,
    blend: u64,
    paint_local: u64,
    blend_local: u64,
}

/// The field's own unit normal at `p`.
///
/// `Sdf::gradient` normalised, which is exactly `marching_cubes::unit_gradient`
/// (`marching_cubes/mod.rs:739-743`) — `PaintStack` does not override
/// `gradient`, so this is the same expression the extractor used for the vertex
/// normals. That identity is what makes the reference *the field's* answer
/// rather than a second opinion of this bench's: at a vertex the reference and
/// the vertex attribute agree by construction, so every disagreement C3 counts
/// is interpolation and nothing else.
fn field_normal<F: Sdf<Scalar = f32>>(field: &F, p: [f32; 3]) -> [f32; 3] {
    let g = field.gradient(p);
    let len = (g[0] * g[0] + g[1] * g[1] + g[2] * g[2]).sqrt();
    if len > 0.0 {
        [g[0] / len, g[1] / len, g[2] / len]
    } else {
        [0.0, 1.0, 0.0]
    }
}

/// Two unit vectors spanning the plane perpendicular to `n`.
fn tangents(n: [f32; 3]) -> ([f32; 3], [f32; 3]) {
    // Cross with whichever axis `n` is least aligned with, so the cross product
    // never degenerates.
    let axis = if n[0].abs() <= n[1].abs() && n[0].abs() <= n[2].abs() {
        [1.0, 0.0, 0.0]
    } else if n[1].abs() <= n[2].abs() {
        [0.0, 1.0, 0.0]
    } else {
        [0.0, 0.0, 1.0]
    };
    let mut t1 = [
        n[1] * axis[2] - n[2] * axis[1],
        n[2] * axis[0] - n[0] * axis[2],
        n[0] * axis[1] - n[1] * axis[0],
    ];
    let len = (t1[0] * t1[0] + t1[1] * t1[1] + t1[2] * t1[2]).sqrt();
    for c in &mut t1 {
        *c /= len;
    }
    let t2 = [
        n[1] * t1[2] - n[2] * t1[1],
        n[2] * t1[0] - n[0] * t1[2],
        n[0] * t1[1] - n[1] * t1[0],
    ];
    (t1, t2)
}

/// C3's instrument: 10^6 points sampled uniformly by area over the extracted
/// mesh, each classified twice.
///
/// **The approximation** is the barycentric blend of the three vertex weights,
/// which is what a rasteriser hands a fragment. **The reference is the field's
/// own answer**: `exact_weights` at the same world position, with the normal
/// taken from `field_normal` — no mesh quantity enters it, which is what
/// *"against the field"* has to mean.
///
/// # Locality, and why the obvious probe is the wrong one
///
/// The first version of this probed the exact classifier six times at one cell
/// along each **world axis**, holding the normal fixed at the interpolated one.
/// It reported 77-87% of misclassifications as boundary-local, and the shortfall
/// was the instrument rather than the weights: `base_weights` depends on the
/// **surface normal** (`triplanar.wgsl:123`), so part of the material boundary
/// lives in normal space and a probe that freezes the normal is blind to it.
/// A probe one cell off the surface along `y` is also blind in the other
/// direction — `SPRAY_DEPTH` is 0.05, under half a cell, so *every* painted
/// point trivially loses its paint one cell out and scores local for free.
///
/// So the reported probe walks the **surface**: four steps of one cell in the
/// tangent plane, each pulled back onto the zero set with one Newton step, each
/// re-classified with the field's own normal *there*. A class change proves the
/// point is within one cell, measured along the surface, of a place where the
/// material class changes. The world-axis reading is kept beside it as
/// `misclassified_within_one_cell_axis` precisely because the two disagreeing is
/// the finding.
///
/// The `paint`/`blend` split says which term disagreed: `paint` where a spray
/// contributed at the point (the reference weights differ bit-for-bit from the
/// unpainted stratum blend), `blend` where none did and the disagreement is the
/// stratum blend's own nonlinearity in normal and height.
fn misclassification<F>(world: &World, log: &[WorldEdit], field: &F) -> Misclassification
where
    F: Sdf<Scalar = f32>,
{
    let tris = world.indices.len() / 3;
    let mut cdf = Vec::with_capacity(tris);
    let mut total = 0.0f64;
    for t in 0..tris {
        let a = world.positions[world.indices[t * 3] as usize];
        let b = world.positions[world.indices[t * 3 + 1] as usize];
        let c = world.positions[world.indices[t * 3 + 2] as usize];
        let u = [b[0] - a[0], b[1] - a[1], b[2] - a[2]];
        let v = [c[0] - a[0], c[1] - a[1], c[2] - a[2]];
        let n = [
            u[1] * v[2] - u[2] * v[1],
            u[2] * v[0] - u[0] * v[2],
            u[0] * v[1] - u[1] * v[0],
        ];
        total += 0.5 * f64::from((n[0] * n[0] + n[1] * n[1] + n[2] * n[2]).sqrt());
        cdf.push(total);
    }
    assert!(total > 0.0, "P-75: the extracted world has no surface area");

    let mut rng = Rng(0x5075_0000_0000_0001);
    let mut out = Misclassification {
        points: 0,
        within_one_cell: 0,
        within_one_cell_axis: 0,
        paint: 0,
        blend: 0,
        paint_local: 0,
        blend_local: 0,
    };
    for _ in 0..SURFACE_POINTS {
        let target = rng.next_f64() * total;
        let t = cdf.partition_point(|&c| c < target).min(tris - 1);
        let mut r0 = rng.next_f64();
        let mut r1 = rng.next_f64();
        if r0 + r1 > 1.0 {
            r0 = 1.0 - r0;
            r1 = 1.0 - r1;
        }
        let bary = [1.0 - r0 - r1, r0, r1];
        let i = [
            world.indices[t * 3] as usize,
            world.indices[t * 3 + 1] as usize,
            world.indices[t * 3 + 2] as usize,
        ];
        let mut p = [0.0f32; 3];
        let mut w = [0.0f32; MATERIALS];
        for (&bk, &vi) in bary.iter().zip(i.iter()) {
            let b = bk as f32;
            for (slot, &pv) in p.iter_mut().zip(world.positions[vi].iter()) {
                *slot += b * pv;
            }
            for (slot, &wv) in w.iter_mut().zip(world.weights[vi].iter()) {
                *slot += b * wv;
            }
        }
        let interpolated = argmax(&w);
        let n = field_normal(field, p);
        let (exact, _) = exact_weights(log, p, n[1]);
        let exact_class = argmax(&exact);
        if interpolated == exact_class {
            continue;
        }
        out.points += 1;

        // Which term disagreed. A spray with zero coverage leaves the weights
        // bit-identical to the unpainted blend, so this is free rather than a
        // second walk.
        let bare = base_weights(p[1], n[1]);
        let painted = exact
            .iter()
            .zip(bare.iter())
            .any(|(a, b)| a.to_bits() != b.to_bits());

        // The world-axis reading, kept for comparison.
        let mut axis_local = false;
        for axis in 0..3 {
            for step in [-CELL_SIZE, CELL_SIZE] {
                let mut q = p;
                q[axis] += step;
                let (probe, _) = exact_weights(log, q, n[1]);
                if argmax(&probe) != exact_class {
                    axis_local = true;
                }
            }
        }
        if axis_local {
            out.within_one_cell_axis += 1;
        }

        // The surface reading.
        let (t1, t2) = tangents(n);
        let mut local = false;
        for dir in [t1, t2] {
            for step in [-CELL_SIZE, CELL_SIZE] {
                let mut q = [
                    p[0] + step * dir[0],
                    p[1] + step * dir[1],
                    p[2] + step * dir[2],
                ];
                // One Newton step back onto the zero set, so the probe stays on
                // the surface where curvature has carried it off.
                let g = field.gradient(q);
                let g2 = g[0] * g[0] + g[1] * g[1] + g[2] * g[2];
                if g2 > 0.0 {
                    let k = field.sample(q) / g2;
                    for axis in 0..3 {
                        q[axis] -= k * g[axis];
                    }
                }
                let nq = field_normal(field, q);
                let (probe, _) = exact_weights(log, q, nq[1]);
                if argmax(&probe) != exact_class {
                    local = true;
                }
            }
        }
        if local {
            out.within_one_cell += 1;
        }
        if painted {
            out.paint += 1;
            if local {
                out.paint_local += 1;
            }
        } else {
            out.blend += 1;
            if local {
                out.blend_local += 1;
            }
        }
    }
    out
}

type Row = Vec<(&'static str, String)>;

fn main() {
    if !std::env::args().any(|arg| arg == "--bench") {
        return;
    }

    let prereg = isomesh::experiment!("P-75");
    let mhz = cpu_mhz();

    // `M-50`'s four log buckets, taken at each bucket's upper bound.
    let buckets: Vec<Bucket> = [15u32, 30, 45, 60].iter().map(|&b| run_bucket(b)).collect();

    println!(
        "{:>7} {:>7} {:>7} {:>8} {:>9} {:>9} {:>8} {:>10} {:>12} {:>8} {:>7} {:>8}",
        "brushes",
        "sprays",
        "biting",
        "verts",
        "extract",
        "weights",
        "share",
        "frame_vtx",
        "frame_frag",
        "speedup",
        "bounds",
        "misclass"
    );

    // C2's widening is a property of the sweep, not of a row, so it is decided
    // after all four buckets exist.
    let mut widens = true;
    for pair in buckets.windows(2) {
        let a = pair[0].frame_ms_fragment / pair[0].frame_ms_vertex;
        let b = pair[1].frame_ms_fragment / pair[1].frame_ms_vertex;
        if b <= a {
            widens = false;
        }
    }

    let mut rows: Vec<Row> = Vec::new();
    for (k, b) in buckets.iter().enumerate() {
        // The registered vacuity control, and the two integrity controls.
        assert!(
            b.boundary_vertices > 0,
            "P-75: no material boundary crosses the visible surface at {} brushes, so the \
             misclassification rate is measured over a single-material scene and every clause \
             is vacuous",
            b.brushes
        );
        assert!(
            b.materials >= 2,
            "P-75: only {} material(s) on the visible surface at {} brushes",
            b.materials,
            b.brushes
        );
        assert_eq!(
            b.field_walk_mismatches, 0,
            "P-75: the timed weight walk does not carry the same field as PaintStack::sample, \
             so it is not M-138's walk"
        );
        assert!(
            b.l1_max_error < 1e-4,
            "P-75: vertex weights are not L1-normalised (max error {})",
            b.l1_max_error
        );
        assert!(
            b.visible_pixels > 0,
            "P-75: nothing is on screen, so both shading paths shade nothing"
        );

        let weights_share = b.weights_ms / b.extract_ms;
        let speedup = b.frame_ms_fragment / b.frame_ms_vertex;
        let misclassified_fraction = b.misclass.points as f64 / SURFACE_POINTS as f64;
        let local_fraction = if b.misclass.points == 0 {
            1.0
        } else {
            b.misclass.within_one_cell as f64 / b.misclass.points as f64
        };
        let c1 = weights_share < 0.05;
        let c2 = speedup >= 4.0 && widens;
        // 0.95 is this harness's reading of "concentrated", stated rather than
        // implied: the registration says boundary-local without a number, and a
        // clause with no threshold cannot be scored.
        let c3 = misclassified_fraction < 0.02 && local_fraction >= 0.95;

        println!(
            "{:>7} {:>7} {:>7} {:>8} {:>9.3} {:>9.3} {:>8.4} {:>10.4} {:>12.4} {:>8.2} {:>7} {:>8}",
            b.brushes,
            b.sprays,
            b.brushes_biting,
            b.vertices,
            b.extract_ms,
            b.weights_ms,
            weights_share,
            b.frame_ms_vertex,
            b.frame_ms_fragment,
            speedup,
            b.boundary_vertices,
            b.misclass.points
        );

        rows.push(vec![
            ("log_bucket", format!("{}", k + 1)),
            ("brushes", b.brushes.to_string()),
            ("sprays", b.sprays.to_string()),
            ("log_len", b.log_len.to_string()),
            ("brushes_biting", b.brushes_biting.to_string()),
            ("materials", b.materials.to_string()),
            ("boundary_vertices", b.boundary_vertices.to_string()),
            (
                "boundary_vertices_all",
                b.boundary_vertices_all.to_string(),
            ),
            ("samples", b.samples.to_string()),
            ("vertices", b.vertices.to_string()),
            ("triangles", b.triangles.to_string()),
            ("visible_triangles", b.visible_triangles.to_string()),
            (
                "weights_ns_per_vertex",
                format!("{:.4}", b.weights_ms * 1e6 / b.vertices as f64),
            ),
            (
                "extract_ns_per_sample",
                format!("{:.4}", b.extract_ms * 1e6 / b.samples as f64),
            ),
            ("weights_ms", format!("{:.6}", b.weights_ms)),
            ("extract_ms", format!("{:.6}", b.extract_ms)),
            ("weights_share", format!("{weights_share:.6}")),
            ("raster_ms", format!("{:.6}", b.raster_ms)),
            ("interp_ms", format!("{:.6}", b.interp_ms)),
            ("frame_ms_vertex", format!("{:.6}", b.frame_ms_vertex)),
            ("frame_ms_fragment", format!("{:.6}", b.frame_ms_fragment)),
            ("speedup", format!("{speedup:.6}")),
            ("speedup_widens", widens.to_string()),
            ("screen_pixels", (SCREEN_W * SCREEN_H).to_string()),
            ("visible_pixels", b.visible_pixels.to_string()),
            ("clipped_triangles", b.clipped_triangles.to_string()),
            (
                "vertex_ns_per_pixel",
                format!("{:.4}", b.frame_ms_vertex * 1e6 / b.visible_pixels as f64),
            ),
            (
                "fragment_ns_per_pixel",
                format!("{:.4}", b.frame_ms_fragment * 1e6 / b.visible_pixels as f64),
            ),
            (
                "fragment_ns_per_log_entry",
                format!(
                    "{:.4}",
                    b.frame_ms_fragment * 1e6 / b.visible_pixels as f64 / b.log_len as f64
                ),
            ),
            ("pixel_disagreement", b.pixel_disagreement.to_string()),
            (
                "pixel_disagreement_fraction",
                format!(
                    "{:.6}",
                    b.pixel_disagreement as f64 / b.visible_pixels as f64
                ),
            ),
            ("surface_points", SURFACE_POINTS.to_string()),
            ("misclassified_points", b.misclass.points.to_string()),
            (
                "misclassified_fraction",
                format!("{misclassified_fraction:.6}"),
            ),
            (
                "misclassified_within_one_cell",
                b.misclass.within_one_cell.to_string(),
            ),
            ("within_one_cell_fraction", format!("{local_fraction:.6}")),
            (
                "misclassified_within_one_cell_axis",
                b.misclass.within_one_cell_axis.to_string(),
            ),
            ("misclassified_paint", b.misclass.paint.to_string()),
            ("misclassified_blend", b.misclass.blend.to_string()),
            (
                "misclassified_paint_local",
                b.misclass.paint_local.to_string(),
            ),
            (
                "misclassified_blend_local",
                b.misclass.blend_local.to_string(),
            ),
            ("speedup_min", format!("{:.6}", b.speedup_min)),
            ("speedup_max", format!("{:.6}", b.speedup_max)),
            ("shade_reps", SHADE_REPS.to_string()),
            ("field_walk_mismatches", b.field_walk_mismatches.to_string()),
            ("l1_max_error", format!("{:.3e}", b.l1_max_error)),
            ("cpu_mhz", format!("{mhz:.3}")),
            ("reps", REPS.to_string()),
            ("c1_holds", c1.to_string()),
            ("c2_holds", c2.to_string()),
            ("c3_holds", c3.to_string()),
        ]);
    }

    common::experiment::run(prereg, |run| {
        for row in &rows {
            run.record(row);
        }
    });
}
