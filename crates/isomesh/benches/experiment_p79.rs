//! **P-79 — how many shadow-atlas pages one dig invalidates.**
//!
//! Ticket: R-079. Pre-registered before this harness existed; the registration
//! is in `crates/isomesh/src/experiment.rs` and is not amended here.
//!
//! ```bash
//! cargo bench --bench experiment_p79
//! ```
//!
//! Writes `docs/experiments/p-79.csv`.
//!
//! # The hypothesis, and what falsifies it
//!
//! **Epic's own answer to deforming geometry was to stop caching.** The Fortnite
//! virtual-shadow-map write-up gives exactly one quantitative figure — a light
//! loop taken from 1.56 ms to 1.08 ms — and describes invalidation only
//! qualitatively: sun movement causes "quite significant shadow page table
//! changes frame to frame", animated deformation invalidates pages, and
//! directional sun shadows were left uncached. For a world being dug away that
//! is the relevant precedent and the number is missing from the literature.
//!
//! - **C1.** One brush stroke invalidates a page count proportional to the
//!   brush's **projected** area from the light, with a constant under 3 — not to
//!   the brush volume and not to the scene. *Falsified by* a constant above 3,
//!   or by scene-proportional invalidation.
//! - **C2.** Invalidating only the pages the brush's light-space bounding volume
//!   touches produces a **pixel-identical** shadow to a full re-render, on all
//!   eight reference fields. *Falsified by* any pixel difference, which would
//!   mean the conservative bound is wrong and localised invalidation is unsound.
//! - **C3.** Cached-with-invalidation beats uncached by at least **2×** in
//!   shadow cost under a continuous dig at `game_dig`'s throttled 12.5 strokes
//!   per second. *Falsified by* under 2×, which vindicates Epic's decision.
//!
//! **VACUITY CONTROL.** The light must be positioned so the brush casts a
//! shadow into the visible frame, asserted by a non-zero `changed_pixels_full`.
//!
//! # What the instrument is
//!
//! There is no Bevy renderer here — `crates/isomesh` must not depend on Bevy —
//! so the shadow map is built from first principles out of geometry this crate
//! produces, and every number below is a property of **this** shadow pass rather
//! than of any engine's frame.
//!
//! A virtual shadow map is a light-space depth image partitioned into pages.
//! So:
//!
//! 1. **The geometry** is the crate's own marching-cubes extraction of a
//!    reference field over a [`CELLS`]³ grid on the field's canonical domain —
//!    the same triangles a renderer would be handed.
//! 2. **The depth image** is an orthographic light-space depth buffer, produced
//!    by a scanline rasteriser over those triangles: for each texel, the
//!    smallest `dot(p − centre, light_dir)` over every triangle covering the
//!    texel centre, and `+∞` where nothing does. That is a shadow map's
//!    definition, not an approximation of one; an orthographic projection makes
//!    depth linear in screen space so the barycentric interpolation is exact.
//! 3. **The pages** are square blocks of [`PAGE_TEXELS`] texels. UE5's virtual
//!    shadow maps use 128×128 physical pages, which is why 128 is the primary
//!    arm and the other three sizes are swept beside it.
//! 4. **The edit** is one `game_dig` brush stroke: a sphere subtracted from the
//!    field, `max(f, −(|p − c| − r))`, centred on the surface point the light
//!    can see nearest the middle of its frame.
//! 5. **`shadow_ms_cached` / `shadow_ms_uncached`** are the cost of *this*
//!    rasteriser, in milliseconds of shadow-pass work per second of continuous
//!    digging at 12.5 strokes/s. They are **not** engine frame times and no
//!    frame rate is invented anywhere in this harness — both arms are charged at
//!    the same 12.5 Hz, so their ratio is the caching benefit alone.
//!
//! The localised arm is a real virtual-shadow-map update and not a shortcut: the
//! invalidated pages are cleared, every triangle is tested against them, and the
//! survivors are rasterised clipped to those pages. The per-triangle cull is
//! O(triangles) in both arms, which is exactly the light loop Epic's one
//! published figure optimised, so it is charged to the cached arm rather than
//! wished away. What this rasteriser does **not** charge is a GPU's draw-call
//! and hierarchical-culling overhead, so the measured ratio is an upper bound on
//! the achievable saving — the direction that favours the hypothesis, stated so
//! the verdict can be read with it in mind.
//!
//! # `game_dig`'s constants, restated
//!
//! The bench must not depend on Bevy, so the numbers are copied with their
//! source named rather than imported: `bevy_isomesh/examples/game_dig.rs` line
//! 127 `CELL_SIZE = 0.125`, line 1056 `World::radius = 0.25`, line 2581
//! `clamp(0.10, 2.00)` on the wheel, line 342 `EDIT_PERIOD = 0.08`. So the
//! default brush is **2 cells** in radius, the wheel reaches **16 cells**, and a
//! held button lays **12.5 strokes a second**.
//!
//! # The SHARE line, recomputed before the code was written
//!
//! C1's constant is `pages_invalidated / (π r² / page_area)`, and for a sphere
//! that is a pure quantisation ratio: a disc of diameter `d` pages has an
//! axis-aligned page footprint of `(⌈d⌉+1)²` at worst and `⌈d⌉²` at best, so the
//! constant is bounded below by `4/π ≈ 1.27` and rises without limit as the
//! brush shrinks below one page. Worked out in advance, in cells, with the
//! one-cell dilation this harness uses for the sound bound:
//!
//! | brush radius | page = 2 cells | 4 cells | 8 cells | 16 cells |
//! |---:|---:|---:|---:|---:|
//! | 2 cells (`game_dig` default) | 2.86 – 5.09 | 5.09 | 5.09 | **20.37** |
//! | 8 cells | 1.61 – 1.99 | 1.99 | 2.86 | 5.09 |
//! | 16 cells (wheel maximum) | 1.44 – 1.61 | 1.61 | 1.99 | 2.86 |
//!
//! **So C1's "under 3" is arithmetically unreachable at `game_dig`'s default
//! brush for any page of 4 cells or more, and unreachable at UE5's 128-texel
//! page for every brush the wheel can reach except the largest.** That is
//! `✗51`'s rule applied before the run: the clause is decided by a
//! brush-diameter-to-page-size ratio the registration never named. The run
//! happens anyway, to produce the number.
//!
//! C3's ratio, by contrast, has no such ceiling — it is the pixel-count ratio
//! between a full frame and a few pages, damped by the O(triangles) cull the
//! cached arm still pays. Both arms are charged at 12.5 Hz precisely so that the
//! ratio cannot be inflated by a frame rate this harness has not measured.
//!
//! # Controls, each an assertion rather than a printed number
//!
//! - **The buffer is calibrated against closed-form geometry before anything
//!   is measured in it.** [`calibrate`] rasterises the unit sphere and checks
//!   the covered-texel count against `π R² / texel²`, which is what its
//!   projected area is from every direction. Nothing else here would catch a
//!   uniform scale error, and every number in the row is a texel count.
//! - **The registered vacuity control.** Every dig row must have
//!   `changed_pixels_full > 0`: a light that cannot see the dig would report a
//!   perfect localisation over an empty change set.
//! - **That zero has to be reachable, or the control is `M-44`.** So the harness
//!   also digs *buried* brushes — one per field, at the deepest interior point
//!   along the same light ray — and asserts that at least one of them produces
//!   `changed_pixels_full == 0`. A buried cavity is behind the first hit, so a
//!   first-hit depth buffer must not move; if the depth buffer were a
//!   last-hit buffer, or the frame were mis-registered, this is the row that
//!   says so.
//! - **Page quantisation must not be allowed to hide an unsound bound.** Whether
//!   the localised composite matches the full re-render depends on where the
//!   brush happens to fall relative to a page boundary — a brush well inside a
//!   page leaks nothing even when the bound is far too tight — so C2 is measured
//!   a second way that page size cannot rescue. Every changed texel's distance
//!   from the projected brush centre is taken, and the harness reports how many
//!   fall outside the undilated disc
//!   (`changed_pixels_outside_tight_disc`), outside the disc grown by
//!   [`DILATE_CELLS`] (`changed_pixels_outside_dilated_disc`), outside the
//!   derived disc of [`SOUND_DILATE_CELLS`]
//!   (`changed_pixels_outside_sound_disc`), and the largest excess seen at all
//!   (`max_change_excess_cells`, and the same divided by the cell diagonal as
//!   `excess_in_cell_diagonals`). Those are statements about the geometry of the
//!   bound and are independent of every page size in the sweep.
//! - **Area against volume, decided by a sweep rather than asserted.** Three
//!   brush radii spanning `game_dig`'s whole wheel, with
//!   `pages_per_projected_area_cell` and `pages_per_volume_cell` beside each
//!   other: whichever is stable across the sweep identifies the law.
//! - **Direction, swept for the same reason.** A sphere's projected area is
//!   direction-invariant, so `pages_invalidated` must be too, up to page
//!   straddling. Four light directions test that instead of assuming it.

#![allow(
    clippy::cast_possible_truncation,
    clippy::cast_precision_loss,
    clippy::cast_sign_loss,
    clippy::too_many_lines
)]

mod common;

use std::time::Instant;

use isomesh::fields::ReferenceField;
use isomesh::for_each_reference_field;
use isomesh::marching_cubes::MarchingCubes;
use isomesh::{MeshBuffer, RuntimeShape3, Sdf};

// ─── game_dig's constants, restated with their source ───────────────────────

/// `bevy_isomesh/examples/game_dig.rs:127` — `const CELL_SIZE: f32 = 0.125`.
const GAME_DIG_CELL: f64 = 0.125;

/// `game_dig.rs:1056` — the `World::radius` the demo starts with.
const GAME_DIG_BRUSH: f64 = 0.25;

/// `game_dig.rs:2581` — the upper end of the mouse-wheel clamp on the radius.
const GAME_DIG_BRUSH_MAX: f64 = 2.00;

/// `game_dig.rs:342` — `const EDIT_PERIOD: f32 = 0.08`, i.e. 12.5 a second.
const GAME_DIG_EDIT_PERIOD: f64 = 0.08;

/// Strokes a second while a button is held. C3's rate.
const STROKE_HZ: f64 = 1.0 / GAME_DIG_EDIT_PERIOD;

// ─── the instrument's own constants ─────────────────────────────────────────

/// Cells per axis in the extraction grid, over each field's canonical domain.
///
/// 64 because that is the crate's canonical grid — `ThinPlate::
/// CANONICAL_CELL_SIZE` is defined as `2 · COMPACT_DOMAIN / 64` — so the
/// geometry this shadow map covers is the geometry every other measurement in
/// the repository is taken on.
const CELLS: u32 = 64;

/// Shadow texels per extraction cell.
///
/// A shadow map coarser than the mesh's own cell cannot resolve the geometry it
/// shadows, and one much finer resolves detail a 64³ extraction does not carry;
/// eight puts the texel at an eighth of a cell, which at `game_dig`'s 12.5 cm
/// cell is 1.6 cm — the order of the screen-pixel footprint a UE5 clipmap level
/// targets in a first-person view, and the reason 128-texel pages are a
/// sensible primary arm.
const TEXELS_PER_CELL: u32 = 8;

/// Shadow-map resolution per axis.
///
/// The light-space frame must circumscribe the domain cube from **every**
/// direction, so its side is at least `64·√3 = 110.85` cells. 896 texels is 112
/// cells at [`TEXELS_PER_CELL`], and 896 is divisible by all four page sizes —
/// a partial page at the frame edge would make `pages_total` a fiction. Holding
/// the frame fixed in cells rather than fitting it per direction is what keeps
/// the texel size, and therefore the projected area in pages, identical across
/// the direction sweep.
const FRAME_TEXELS: u32 = 896;

/// Page sizes swept, in texels. 128 is UE5's physical virtual-shadow-map page.
const PAGE_TEXELS: [u32; 4] = [16, 32, 64, 128];

/// Brush radii swept, in cells: `game_dig`'s default, the middle of its wheel,
/// and the wheel's maximum. `0.25 / 0.125 = 2` and `2.00 / 0.125 = 16`.
const BRUSH_CELLS: [f64; 3] = [
    GAME_DIG_BRUSH / GAME_DIG_CELL,
    8.0,
    GAME_DIG_BRUSH_MAX / GAME_DIG_CELL,
];

/// Cells the brush's bounding volume is grown by before it selects pages.
///
/// One, because that is `game_dig`'s own invalidation reach: `let reach =
/// brush.shape.radius + cell` at `game_dig.rs:2682`, and it is there because
/// marching cubes places a vertex anywhere inside a cell whose corner values
/// changed — so the geometry that moves extends one cell beyond the brush. The
/// undilated bound is measured too, as `pages_invalidated_tight`, precisely so
/// this cell of slack is a number rather than an assumption.
const DILATE_CELLS: f64 = 1.0;

/// Cells the *derived sound* bound grows the brush by: `2·√3`.
///
/// Derived rather than fitted, and the derivation is the mechanism this row
/// turned out to be about. A subtractive brush replaces `f` with
/// `max(f, r − |p − c|)`, so a grid corner `p`'s value changes iff
/// `|p − c| < r − f(p)` — which for a corner **inside** the solid, where
/// `f(p) < 0`, reaches *further than the brush*. A corner deeper than one
/// brush radius still moves, and moving a corner moves every marching-cubes
/// vertex on an edge it terminates, because the crossing parameter
/// `t = f_a / (f_a − f_b)` is a function of the values and not just their signs.
///
/// Only cells that hold a sign change carry geometry, and on such a cell the
/// corners are within one cell diagonal of the surface, so
/// `|f(p)| ≤ √3 · cell · L` for an `L`-Lipschitz field. That puts the changed
/// corner inside `r + √3 · cell · L`, and a marching-cubes triangle is confined
/// to its own cell, whose farthest point is another `√3 · cell` away. So the
/// sound reach is `r + √3 · (1 + L) · cell`, which is `r + 2√3 · cell` for the
/// 1-Lipschitz fields — six of the eight, since `FieldBound::lipschitz`
/// answers `1` for `Exact` and for `Underestimate` alike.
///
/// The two that are not are the prediction this constant makes falsifiable:
/// `gyroid` registers `Lipschitz { l: 3.4641 }`, which needs `√3 · 4.4641 =
/// 7.73` cells, and `fbm_terrain` registers `Unbounded`, for which **no light-
/// space invalidation radius can be derived at all** — only observed.
const SOUND_DILATE_CELLS: f64 = 2.0 * 1.732_050_807_568_877_2;

/// Timed repetitions per arm, median taken.
const REPS: usize = 5;

/// Light directions swept. Not normalised here; [`Light::new`] does that.
///
/// A sphere's projected area does not depend on the direction it is projected
/// along, so `pages_invalidated` must not either. Four directions — overhead,
/// a 45° sun, an off-axis one with no symmetry to hide behind, and a low
/// grazing one — turn that into a measurement.
const LIGHTS: [(&str, [f64; 3]); 4] = [
    ("overhead", [0.0, -1.0, 0.0]),
    ("low_45", [1.0, -1.0, 0.0]),
    ("oblique", [0.6, -0.7, 0.4]),
    ("grazing", [1.0, -0.25, 0.6]),
];

// ─── small vector helpers ───────────────────────────────────────────────────

#[inline]
fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

#[inline]
fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

#[inline]
fn unit(a: [f64; 3]) -> [f64; 3] {
    let l = dot(a, a).sqrt();
    [a[0] / l, a[1] / l, a[2] / l]
}

// ─── the dug field ──────────────────────────────────────────────────────────

/// One `game_dig` brush stroke: a sphere subtracted from a field.
///
/// `max(f, −(|p − c| − r))`, the crate's own difference. Written out here rather
/// than composed through `BrushStack` because the harness needs the field before
/// and after one stroke as two independent `Sdf`s.
struct Dug<'a, F> {
    base: &'a F,
    centre: [f64; 3],
    radius: f64,
}

impl<F: Sdf<Scalar = f64>> Sdf for Dug<'_, F> {
    type Scalar = f64;

    #[inline]
    fn sample(&self, p: [f64; 3]) -> f64 {
        let d = [
            p[0] - self.centre[0],
            p[1] - self.centre[1],
            p[2] - self.centre[2],
        ];
        let inside_sphere = self.radius - dot(d, d).sqrt();
        let f = self.base.sample(p);
        if f > inside_sphere { f } else { inside_sphere }
    }

    /// The gradient of the active branch.
    ///
    /// `∇|p − c|` is undefined at `c`, and this never evaluates it there: the
    /// sphere branch is active only where `r − |p − c| ≥ f`, and an extractor
    /// asks for a gradient only at a vertex, where the composite is zero — so
    /// on that branch `|p − c| = r > 0`.
    #[inline]
    fn gradient(&self, p: [f64; 3]) -> [f64; 3] {
        let d = [
            p[0] - self.centre[0],
            p[1] - self.centre[1],
            p[2] - self.centre[2],
        ];
        let len = dot(d, d).sqrt();
        let inside_sphere = self.radius - len;
        let f = self.base.sample(p);
        if f > inside_sphere {
            self.base.gradient(p)
        } else {
            [-d[0] / len, -d[1] / len, -d[2] / len]
        }
    }
}

// ─── the light-space depth image ────────────────────────────────────────────

/// An orthographic light-space frame: an orthonormal basis, a texel size and a
/// resolution.
struct Light {
    /// Frame's first tangent axis.
    u: [f64; 3],
    /// Frame's second tangent axis.
    v: [f64; 3],
    /// Direction the light travels. Depth increases along it.
    d: [f64; 3],
    /// World point the frame is centred on.
    centre: [f64; 3],
    /// Light-space coordinate of the frame's low corner, on both axes.
    origin: f64,
    /// World size of one texel.
    texel: f64,
    /// Texels per axis.
    res: u32,
}

impl Light {
    fn new(direction: [f64; 3], centre: [f64; 3], texel: f64, res: u32) -> Self {
        let d = unit(direction);
        // The world axis least aligned with `d`, so the cross product is far
        // from degenerate whatever direction was asked for.
        let helper = if d[0].abs() <= d[1].abs() && d[0].abs() <= d[2].abs() {
            [1.0, 0.0, 0.0]
        } else if d[1].abs() <= d[2].abs() {
            [0.0, 1.0, 0.0]
        } else {
            [0.0, 0.0, 1.0]
        };
        let u = unit(cross(d, helper));
        let v = cross(d, u);
        Self {
            u,
            v,
            d,
            centre,
            origin: -0.5 * f64::from(res) * texel,
            texel,
            res,
        }
    }

    /// `(tangent a, tangent b, depth)` of a world point.
    #[inline]
    fn project(&self, p: [f64; 3]) -> (f64, f64, f64) {
        let r = [
            p[0] - self.centre[0],
            p[1] - self.centre[1],
            p[2] - self.centre[2],
        ];
        (dot(r, self.u), dot(r, self.v), dot(r, self.d))
    }

    /// World point of a texel centre at a given depth.
    #[inline]
    fn unproject(&self, x: u32, y: u32, depth: f64) -> [f64; 3] {
        let a = self.origin + (f64::from(x) + 0.5) * self.texel;
        let b = self.origin + (f64::from(y) + 0.5) * self.texel;
        [
            self.centre[0] + a * self.u[0] + b * self.v[0] + depth * self.d[0],
            self.centre[1] + a * self.u[1] + b * self.v[1] + depth * self.d[1],
            self.centre[2] + a * self.u[2] + b * self.v[2] + depth * self.d[2],
        ]
    }

    /// Light-space tangent coordinate of a texel centre, on either axis.
    #[inline]
    fn tangent_of(&self, index: u32) -> f64 {
        self.origin + (f64::from(index) + 0.5) * self.texel
    }

    /// Texel index containing a light-space tangent coordinate, clamped to the
    /// frame. Clamping is the geometry of a bounded atlas, not a fallback: a
    /// page outside the frame does not exist to invalidate.
    #[inline]
    fn texel_of(&self, coord: f64) -> u32 {
        let i = ((coord - self.origin) / self.texel).floor();
        i.clamp(0.0, f64::from(self.res - 1)) as u32
    }

    #[inline]
    fn whole(&self) -> Rect {
        Rect {
            x0: 0,
            y0: 0,
            x1: self.res - 1,
            y1: self.res - 1,
        }
    }
}

/// An inclusive texel rectangle. The full pass uses the whole frame; the
/// localised pass uses the invalidated page block — which is always contiguous,
/// because it is the page footprint of an axis-aligned bounding box.
#[derive(Clone, Copy)]
struct Rect {
    x0: u32,
    y0: u32,
    x1: u32,
    y1: u32,
}

/// Render a mesh into `buf`, restricted to `clip`.
///
/// One function for both arms, which is the point: the full pass is this with
/// `clip` set to the whole frame, and the localised pass is this with `clip` set
/// to the invalidated pages. There is no second code path whose disagreement
/// with the first could be mistaken for a result.
///
/// `clip` is cleared to `+∞` first, because re-rendering a page starts from
/// nothing — a virtual shadow map does not accumulate into a stale page. The
/// per-triangle rejection is the light loop's cull and is charged to whichever
/// arm is being timed.
fn render(mesh: &MeshBuffer<f64>, light: &Light, buf: &mut [f64], clip: Rect) {
    let res = light.res as usize;
    for y in clip.y0..=clip.y1 {
        let row = y as usize * res;
        buf[row + clip.x0 as usize..=row + clip.x1 as usize].fill(f64::INFINITY);
    }

    let clip_x0 = f64::from(clip.x0);
    let clip_x1 = f64::from(clip.x1);
    let clip_y0 = f64::from(clip.y0);
    let clip_y1 = f64::from(clip.y1);
    let inv_texel = 1.0 / light.texel;

    for tri in mesh.indices.as_chunks::<3>().0 {
        let mut px = [0.0f64; 3];
        let mut py = [0.0f64; 3];
        let mut pw = [0.0f64; 3];
        for k in 0..3 {
            let (a, b, w) = light.project(mesh.positions[tri[k] as usize]);
            px[k] = (a - light.origin) * inv_texel;
            py[k] = (b - light.origin) * inv_texel;
            pw[k] = w;
        }

        // Texel centres sit at `i + 0.5`, so texel `i` is covered iff
        // `i + 0.5` is inside the triangle's span.
        let lo_x = (px[0].min(px[1]).min(px[2]) - 0.5).ceil().max(clip_x0);
        let hi_x = (px[0].max(px[1]).max(px[2]) - 0.5).floor().min(clip_x1);
        let lo_y = (py[0].min(py[1]).min(py[2]) - 0.5).ceil().max(clip_y0);
        let hi_y = (py[0].max(py[1]).max(py[2]) - 0.5).floor().min(clip_y1);
        if hi_x < lo_x || hi_y < lo_y {
            continue;
        }

        let area2 = (px[1] - px[0]) * (py[2] - py[0]) - (px[2] - px[0]) * (py[1] - py[0]);
        if area2.abs() <= 0.0 {
            continue;
        }
        let inv = 1.0 / area2;
        // Edge functions written as `A + B·x + C·y`, normalised by the signed
        // doubled area so the inside test is winding-independent.
        let a0 = (px[1] * py[2] - px[2] * py[1]) * inv;
        let b0 = (py[1] - py[2]) * inv;
        let c0 = (px[2] - px[1]) * inv;
        let a1 = (px[2] * py[0] - px[0] * py[2]) * inv;
        let b1 = (py[2] - py[0]) * inv;
        let c1 = (px[0] - px[2]) * inv;
        let dw0 = pw[0] - pw[2];
        let dw1 = pw[1] - pw[2];

        let mut y = lo_y as u32;
        let y_end = hi_y as u32;
        while y <= y_end {
            let yc = f64::from(y) + 0.5;
            let row = y as usize * res;
            let l0_row = a0 + c0 * yc;
            let l1_row = a1 + c1 * yc;
            let mut x = lo_x as u32;
            let x_end = hi_x as u32;
            while x <= x_end {
                let xc = f64::from(x) + 0.5;
                let l0 = l0_row + b0 * xc;
                let l1 = l1_row + b1 * xc;
                if l0 >= 0.0 && l1 >= 0.0 && l0 + l1 <= 1.0 {
                    let depth = pw[2] + l0 * dw0 + l1 * dw1;
                    let slot = &mut buf[row + x as usize];
                    if depth < *slot {
                        *slot = depth;
                    }
                }
                x += 1;
            }
            y += 1;
        }
    }
}

/// Texels whose depth differs, bit for bit.
///
/// Bits rather than `==`: `+∞` must compare equal to `+∞`, and two depths that
/// differ in the last place are a difference this experiment is not allowed to
/// round away.
fn differing(a: &[f64], b: &[f64]) -> usize {
    a.iter()
        .zip(b)
        .filter(|(x, y)| x.to_bits() != y.to_bits())
        .count()
}

/// How far outside its own bounding disc one stroke reached, in the light's own
/// frame.
///
/// This is C2's page-size-independent instrument, and the one that decided the
/// row. Whether the localised composite matches a full re-render can be rescued
/// by luck — a brush that happens to sit well inside a page leaks nothing even
/// when the bound is far too tight. This asks the geometric question instead:
/// how far from the brush did the shadow actually move?
///
/// One pass over the frame, because three passes for three radii would read six
/// megabytes of depth buffer three times to answer one question.
struct Leak {
    /// Changed texels outside the undilated projected disc.
    outside_tight: usize,
    /// Changed texels outside the disc grown by [`DILATE_CELLS`].
    outside_dilated: usize,
    /// Changed texels outside the derived disc of [`SOUND_DILATE_CELLS`].
    outside_sound: usize,
    /// Largest `|projected offset| − radius` over changed texels, in cells.
    max_excess_cells: f64,
}

fn leak(
    light: &Light,
    before: &[f64],
    after: &[f64],
    brush: [f64; 3],
    radius: f64,
    cell: f64,
) -> Leak {
    let (ba, bb, _) = light.project(brush);
    let tight = radius * radius;
    let dilated = (radius + DILATE_CELLS * cell) * (radius + DILATE_CELLS * cell);
    let sound = (radius + SOUND_DILATE_CELLS * cell) * (radius + SOUND_DILATE_CELLS * cell);
    let res = light.res;
    let mut found = Leak {
        outside_tight: 0,
        outside_dilated: 0,
        outside_sound: 0,
        max_excess_cells: f64::NEG_INFINITY,
    };
    for y in 0..res {
        let db = light.tangent_of(y) - bb;
        let row = y as usize * res as usize;
        for x in 0..res {
            let slot = row + x as usize;
            if before[slot].to_bits() == after[slot].to_bits() {
                continue;
            }
            let da = light.tangent_of(x) - ba;
            let d2 = da * da + db * db;
            if d2 > tight {
                found.outside_tight += 1;
            }
            if d2 > dilated {
                found.outside_dilated += 1;
            }
            if d2 > sound {
                found.outside_sound += 1;
            }
            let excess = (d2.sqrt() - radius) / cell;
            if excess > found.max_excess_cells {
                found.max_excess_cells = excess;
            }
        }
    }
    found
}

/// The invalidated page block, as a texel rectangle and a page count.
///
/// The brush's light-space bounding volume is a disc of `radius` about the
/// projected centre — a sphere projects to a disc along any direction — and the
/// pages it touches are the pages of that disc's axis-aligned bounding box.
fn invalidated(light: &Light, brush: [f64; 3], radius: f64, page: u32) -> (Rect, u32) {
    let (a, b, _) = light.project(brush);
    let x_lo = light.texel_of(a - radius) / page;
    let x_hi = light.texel_of(a + radius) / page;
    let y_lo = light.texel_of(b - radius) / page;
    let y_hi = light.texel_of(b + radius) / page;
    let rect = Rect {
        x0: x_lo * page,
        y0: y_lo * page,
        x1: (x_hi + 1) * page - 1,
        y1: (y_hi + 1) * page - 1,
    };
    (rect, (x_hi - x_lo + 1) * (y_hi - y_lo + 1))
}

/// The texel nearest the frame centre that has any geometry behind it, and the
/// world point on the surface it saw.
///
/// This is where the brush goes: the point on the surface the light can see,
/// nearest the middle of its own frame. It is `game_dig`'s aim point in light
/// space instead of camera space, and centring the brush on the surface is what
/// the demo does — `game_dig.rs:4784` notes a brush centred on the floor plane
/// reaches a radius below it.
fn aim(light: &Light, depth: &[f64]) -> Option<[f64; 3]> {
    let res = light.res;
    let mid = f64::from(res) * 0.5;
    let mut best: Option<(f64, [f64; 3])> = None;
    for y in 0..res {
        let dy = f64::from(y) + 0.5 - mid;
        let row = y as usize * res as usize;
        for x in 0..res {
            let d = depth[row + x as usize];
            if !d.is_finite() {
                continue;
            }
            let dx = f64::from(x) + 0.5 - mid;
            let score = dx * dx + dy * dy;
            if best.is_none_or(|(b, _)| score < b) {
                best = Some((score, light.unproject(x, y, d)));
            }
        }
    }
    best.map(|(_, p)| p)
}

/// The deepest interior point along the light ray through `from`, and its
/// clearance in cells.
///
/// The buried arm's placement, and the reachable zero the registered vacuity
/// control needs: a cavity entirely behind the first hit cannot move a
/// first-hit depth buffer.
fn deepest_interior<F: Sdf<Scalar = f64>>(
    field: &F,
    light: &Light,
    from: [f64; 3],
    cell: f64,
    span: f64,
) -> ([f64; 3], f64) {
    let step = cell * 0.25;
    let steps = (span / step).ceil() as u32;
    let mut best = (from, 0.0f64);
    for k in 0..=steps {
        let t = f64::from(k) * step;
        let p = [
            from[0] + t * light.d[0],
            from[1] + t * light.d[1],
            from[2] + t * light.d[2],
        ];
        let clearance = -field.sample(p) / cell;
        if clearance > best.1 {
            best = (p, clearance);
        }
    }
    best
}

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

/// Median of `REPS` timed repetitions. Each repetition does its own untimed
/// setup and returns only the milliseconds inside the region.
fn median_ms(mut rep: impl FnMut() -> f64) -> f64 {
    let mut samples: Vec<f64> = (0..REPS).map(|_| rep()).collect();
    samples.sort_by(f64::total_cmp);
    samples[REPS / 2]
}

type Row = Vec<(&'static str, String)>;

/// Everything one (field, light, brush radius) arm measures before the page
/// sweep splits it into rows.
struct Arm {
    changed_full: usize,
    leak: Leak,
    full_ms: f64,
    triangles: usize,
    clearance_cells: f64,
}

/// One field, all four lights, all three radii, all four page sizes, plus the
/// buried control.
fn run_field<F: ReferenceField<Scalar = f64> + Sync>(
    name: &'static str,
    field: &F,
    mhz: f64,
    out: &mut Vec<Row>,
) {
    let (dmin, dmax) = field.domain();
    let cell = (dmax[0] - dmin[0]) / f64::from(CELLS);
    let centre = [
        (dmin[0] + dmax[0]) * 0.5,
        (dmin[1] + dmax[1]) * 0.5,
        (dmin[2] + dmax[2]) * 0.5,
    ];
    let span = (dmax[0] - dmin[0]) * 3.0f64.sqrt();
    let shape = RuntimeShape3::new([CELLS + 1; 3]).expect("65 cubed samples fit in a u32");

    let mut base = MeshBuffer::<f64>::new();
    MarchingCubes::<f64>::new()
        .extract(field, &shape, dmin, cell, &mut base)
        .expect("the reference field meshes on its own canonical grid");
    let base_triangles = base.indices.len() / 3;

    let texel = cell / f64::from(TEXELS_PER_CELL);
    let texels = FRAME_TEXELS as usize * FRAME_TEXELS as usize;

    for (light_name, direction) in LIGHTS {
        let light = Light::new(direction, centre, texel, FRAME_TEXELS);
        let mut cache = vec![f64::INFINITY; texels];
        render(&base, &light, &mut cache, light.whole());
        let covered = cache.iter().filter(|d| d.is_finite()).count();
        let hit = aim(&light, &cache).expect(
            "the light must see the surface: no texel in the frame has geometry behind it, \
             so this field and direction cannot host the experiment",
        );

        // The dig arms, then the buried control on the same ray.
        let (buried_centre, clearance) = deepest_interior(field, &light, hit, cell, span);
        let arms: [([f64; 3], f64, f64, bool); 4] = [
            (hit, BRUSH_CELLS[0], 0.0, false),
            (hit, BRUSH_CELLS[1], 0.0, false),
            (hit, BRUSH_CELLS[2], 0.0, false),
            (buried_centre, BRUSH_CELLS[0], clearance, true),
        ];

        for (index, (brush, radius_cells, clearance_cells, buried)) in arms.into_iter().enumerate() {
            // The buried control is one row per field, on one light, because it
            // is a reachability proof rather than a sweep.
            if buried && light_name != LIGHTS[0].0 {
                continue;
            }
            let radius = radius_cells * cell;
            let dug_field = Dug {
                base: field,
                centre: brush,
                radius,
            };
            let mut dug = MeshBuffer::<f64>::new();
            MarchingCubes::<f64>::new()
                .extract(&dug_field, &shape, dmin, cell, &mut dug)
                .expect("one subtracted sphere does not stop the field from meshing");

            let mut full = vec![f64::INFINITY; texels];
            let full_ms = median_ms(|| {
                let start = Instant::now();
                render(&dug, &light, &mut full, light.whole());
                start.elapsed().as_secs_f64() * 1e3
            });

            let arm = Arm {
                changed_full: differing(&cache, &full),
                leak: leak(&light, &cache, &full, brush, radius, cell),
                full_ms,
                triangles: dug.indices.len() / 3,
                clearance_cells,
            };

            let mut local = vec![f64::INFINITY; texels];
            for page in PAGE_TEXELS {
                // The buried control needs one page size, not four.
                if buried && page != PAGE_TEXELS[PAGE_TEXELS.len() - 1] {
                    continue;
                }
                // Three bounds, one page grid. `dilated` is the registered
                // arm — `game_dig`'s own invalidation reach — `tight` is the
                // literal brush sphere, and `sound` is the derived
                // `r + 2√3·cell`. All three are measured because which of them
                // is conservative is the question, not an assumption.
                let dilated_radius = radius + DILATE_CELLS * cell;
                let sound_radius = radius + SOUND_DILATE_CELLS * cell;
                let (rect, pages) = invalidated(&light, brush, dilated_radius, page);
                let (_, pages_tight) = invalidated(&light, brush, radius, page);
                let (sound_rect, pages_sound) = invalidated(&light, brush, sound_radius, page);

                let local_ms = median_ms(|| {
                    // Untimed: a real atlas keeps its physical pages resident,
                    // so restoring the cached image is not part of the pass.
                    local.copy_from_slice(&cache);
                    let start = Instant::now();
                    render(&dug, &light, &mut local, rect);
                    start.elapsed().as_secs_f64() * 1e3
                });
                let changed_local = differing(&cache, &local);
                let pixel_diff = differing(&full, &local);

                let sound_ms = median_ms(|| {
                    local.copy_from_slice(&cache);
                    let start = Instant::now();
                    render(&dug, &light, &mut local, sound_rect);
                    start.elapsed().as_secs_f64() * 1e3
                });
                let pixel_diff_sound = differing(&full, &local);

                let page_cells = f64::from(page) / f64::from(TEXELS_PER_CELL);
                let projected_cells = std::f64::consts::PI * radius_cells * radius_cells;
                let volume_cells =
                    4.0 / 3.0 * std::f64::consts::PI * radius_cells * radius_cells * radius_cells;
                let projected_pages = projected_cells / (page_cells * page_cells);
                let constant = f64::from(pages) / projected_pages;
                let constant_tight = f64::from(pages_tight) / projected_pages;
                let constant_sound = f64::from(pages_sound) / projected_pages;
                let pages_total = (FRAME_TEXELS / page) * (FRAME_TEXELS / page);
                let ratio = arm.full_ms / local_ms;
                let ratio_sound = arm.full_ms / sound_ms;

                out.push(vec![
                    ("field", name.to_string()),
                    ("light_direction", light_name.to_string()),
                    ("brush_volume_cells", format!("{volume_cells:.4}")),
                    ("brush_projected_area", format!("{projected_cells:.4}")),
                    ("pages_invalidated", pages.to_string()),
                    ("invalidation_constant", format!("{constant:.4}")),
                    ("changed_pixels_full", arm.changed_full.to_string()),
                    ("changed_pixels_localised", changed_local.to_string()),
                    (
                        "shadow_ms_cached",
                        format!("{:.5}", local_ms * STROKE_HZ),
                    ),
                    (
                        "shadow_ms_uncached",
                        format!("{:.5}", arm.full_ms * STROKE_HZ),
                    ),
                    ("stroke_rate_hz", format!("{STROKE_HZ:.1}")),
                    ("c1_holds", (constant < 3.0).to_string()),
                    ("c2_holds", (pixel_diff == 0).to_string()),
                    ("c3_holds", (ratio >= 2.0).to_string()),
                    // ── extras ──
                    (
                        "arm",
                        if buried { "buried_control" } else { "dig" }.to_string(),
                    ),
                    ("arm_index", index.to_string()),
                    ("brush_radius_cells", format!("{radius_cells:.4}")),
                    ("brush_radius_world", format!("{radius:.6}")),
                    ("buried_clearance_cells", format!("{:.4}", arm.clearance_cells)),
                    ("cell_size_world", format!("{cell:.6}")),
                    ("texel_size_world", format!("{texel:.6}")),
                    ("texels_per_cell", TEXELS_PER_CELL.to_string()),
                    ("shadow_texels", FRAME_TEXELS.to_string()),
                    ("page_texels", page.to_string()),
                    ("page_cells", format!("{page_cells:.4}")),
                    ("pages_total", pages_total.to_string()),
                    (
                        "pages_fraction",
                        format!("{:.8}", f64::from(pages) / f64::from(pages_total)),
                    ),
                    ("pages_invalidated_tight", pages_tight.to_string()),
                    (
                        "invalidation_constant_tight",
                        format!("{constant_tight:.4}"),
                    ),
                    (
                        "brush_projected_area_pages",
                        format!("{projected_pages:.6}"),
                    ),
                    (
                        "pages_per_projected_area_cell",
                        format!("{:.6}", f64::from(pages) / projected_cells),
                    ),
                    (
                        "pages_per_volume_cell",
                        format!("{:.6}", f64::from(pages) / volume_cells),
                    ),
                    (
                        "changed_pixels_outside_tight_disc",
                        arm.leak.outside_tight.to_string(),
                    ),
                    (
                        "changed_pixels_outside_dilated_disc",
                        arm.leak.outside_dilated.to_string(),
                    ),
                    (
                        "changed_pixels_outside_sound_disc",
                        arm.leak.outside_sound.to_string(),
                    ),
                    (
                        "max_change_excess_cells",
                        format!("{:.4}", arm.leak.max_excess_cells),
                    ),
                    (
                        "excess_in_cell_diagonals",
                        format!(
                            "{:.4}",
                            arm.leak.max_excess_cells / 1.732_050_807_568_877_2
                        ),
                    ),
                    ("sound_dilate_cells", format!("{SOUND_DILATE_CELLS:.4}")),
                    ("pages_invalidated_sound", pages_sound.to_string()),
                    (
                        "invalidation_constant_sound",
                        format!("{constant_sound:.4}"),
                    ),
                    ("pixel_diff_sound", pixel_diff_sound.to_string()),
                    ("shadow_ms_sound_pass", format!("{sound_ms:.5}")),
                    ("ratio_per_stroke_sound", format!("{ratio_sound:.4}")),
                    ("pixel_diff_localised_vs_full", pixel_diff.to_string()),
                    ("shadow_ms_full_pass", format!("{:.5}", arm.full_ms)),
                    ("shadow_ms_localised_pass", format!("{local_ms:.5}")),
                    ("ratio_per_stroke", format!("{ratio:.4}")),
                    ("triangles_base", base_triangles.to_string()),
                    ("frame_covered_texels", covered.to_string()),
                    ("triangles_dug", arm.triangles.to_string()),
                    ("light_dx", format!("{:.6}", light.d[0])),
                    ("light_dy", format!("{:.6}", light.d[1])),
                    ("light_dz", format!("{:.6}", light.d[2])),
                    ("cpu_mhz", format!("{mhz:.1}")),
                ]);
            }
        }
    }
}

/// Check the depth buffer against closed-form geometry before trusting it.
///
/// Every number in this experiment is a count of texels in a buffer this file
/// fills itself, so the projection scale, the texel size and the coverage test
/// all have to be right or the whole row is fiction — and none of the other
/// controls would notice a uniform factor. The unit sphere has a projected area
/// of exactly `π` world units² from **any** direction, which is
/// `π / texel²` texels, so this is an absolute calibration with no free
/// parameter.
///
/// The tolerance is 3%, and it is a budget rather than a guess: the marching-
/// cubes silhouette chords the true circle, losing `O(cell²/R)` of radius, and
/// the texel-level boundary of the disc is `2πR / texel ≈ 804` texels against
/// an area of `51,472` — 1.6% on its own. A basis error, a halved texel or a
/// dropped `0.5` texel-centre offset all miss by far more than that.
fn calibrate() {
    let field = isomesh::fields::Sphere::<f64>::canonical();
    let (dmin, dmax) = ReferenceField::domain(&field);
    let cell = (dmax[0] - dmin[0]) / f64::from(CELLS);
    let texel = cell / f64::from(TEXELS_PER_CELL);
    let shape = RuntimeShape3::new([CELLS + 1; 3]).expect("65 cubed samples fit in a u32");
    let mut mesh = MeshBuffer::<f64>::new();
    MarchingCubes::<f64>::new()
        .extract(&field, &shape, dmin, cell, &mut mesh)
        .expect("the unit sphere meshes on its canonical grid");

    let light = Light::new([0.0, -1.0, 0.0], [0.0; 3], texel, FRAME_TEXELS);
    let mut depth = vec![f64::INFINITY; FRAME_TEXELS as usize * FRAME_TEXELS as usize];
    render(&mesh, &light, &mut depth, light.whole());
    let covered = depth.iter().filter(|d| d.is_finite()).count();
    let expected = std::f64::consts::PI * field.radius * field.radius / (texel * texel);
    let error = (covered as f64 - expected) / expected;
    println!(
        "calibration: unit sphere covers {covered} texels, closed form {expected:.0}, \
         error {:+.3}% (budget 3%)",
        error * 100.0
    );
    assert!(
        error.abs() < 0.03,
        "P-79 CALIBRATION: the depth buffer covers {covered} texels where the unit sphere's \
         projected area is {expected:.0}, an error of {:+.2}%. The projection scale, the texel \
         size or the coverage rule is wrong and every texel count in this experiment is \
         measured in the wrong unit.",
        error * 100.0
    );
}

fn main() {
    if !std::env::args().any(|a| a == "--bench") {
        return;
    }

    println!(
        "instrument: orthographic light-space depth buffer, {FRAME_TEXELS}x{FRAME_TEXELS} texels \
         at {TEXELS_PER_CELL} texels per cell over a {CELLS}-cell domain, pages of {PAGE_TEXELS:?} \
         texels; geometry is marching cubes over the field's canonical domain."
    );
    println!(
        "share (recomputed before the run): C1's constant is pages / (pi r^2 / page_area), a \
         quantisation ratio bounded below by 4/pi = 1.273 and unbounded above as the brush \
         shrinks below one page. At game_dig's default 2-cell brush it is 5.09 for pages of 2-8 \
         cells and 20.37 for UE5's 128-texel (16-cell) page, so \"under 3\" is arithmetically \
         unreachable there; it is reachable only for brushes whose dilated diameter exceeds \
         ~2 pages. C3 is charged at {STROKE_HZ} Hz on BOTH arms, so no frame rate enters the \
         ratio."
    );
    println!();
    calibrate();
    println!();

    let mhz = cpu_mhz();
    let mut rows: Vec<Row> = Vec::new();
    for_each_reference_field!(f64, |name, field| {
        run_field(name, &field, mhz, &mut rows);
    });

    let get = |row: &Row, key: &str| -> String {
        row.iter()
            .find(|(k, _)| *k == key)
            .map(|(_, v)| v.clone())
            .expect("every row carries every column this summary reads")
    };
    let num = |row: &Row, key: &str| -> f64 { get(row, key).parse().expect("numeric column") };
    let count = |row: &Row, key: &str| -> u64 { get(row, key).parse().expect("integer column") };

    println!(
        "{:>15} {:>9} {:>6} {:>5} {:>7} {:>8} {:>9} {:>9} {:>8} {:>7} {:>7} {:>6}",
        "field",
        "light",
        "r/cell",
        "page",
        "pages",
        "constant",
        "chg_full",
        "chg_local",
        "diff",
        "full_ms",
        "loc_ms",
        "ratio"
    );
    for row in &rows {
        println!(
            "{:>15} {:>9} {:>6} {:>5} {:>7} {:>8} {:>9} {:>9} {:>8} {:>7.3} {:>7.3} {:>6.1}",
            get(row, "field"),
            get(row, "light_direction"),
            get(row, "brush_radius_cells"),
            get(row, "page_texels"),
            get(row, "pages_invalidated"),
            get(row, "invalidation_constant"),
            get(row, "changed_pixels_full"),
            get(row, "changed_pixels_localised"),
            get(row, "pixel_diff_localised_vs_full"),
            num(row, "shadow_ms_full_pass"),
            num(row, "shadow_ms_localised_pass"),
            num(row, "ratio_per_stroke"),
        );
    }

    let dig: Vec<&Row> = rows.iter().filter(|r| get(r, "arm") == "dig").collect();
    let buried: Vec<&Row> = rows
        .iter()
        .filter(|r| get(r, "arm") == "buried_control")
        .collect();

    let mut summary = String::new();
    for page in PAGE_TEXELS {
        for radius in BRUSH_CELLS {
            let sel: Vec<&&Row> = dig
                .iter()
                .filter(|r| {
                    count(r, "page_texels") == u64::from(page)
                        && (num(r, "brush_radius_cells") - radius).abs() < 1e-9
                })
                .collect();
            let k: Vec<f64> = sel.iter().map(|r| num(r, "invalidation_constant")).collect();
            let ratios: Vec<f64> = sel.iter().map(|r| num(r, "ratio_per_stroke")).collect();
            let pages: Vec<u64> = sel.iter().map(|r| count(r, "pages_invalidated")).collect();
            let diffs: u64 = sel
                .iter()
                .map(|r| count(r, "pixel_diff_localised_vs_full"))
                .sum();
            let sound_diffs: u64 = sel.iter().map(|r| count(r, "pixel_diff_sound")).sum();
            summary.push_str(&format!(
                "page {page:>3} texels, r = {radius:>4} cells: pages {}-{}, constant \
                 {:.2}-{:.2}, C1 holds on {}/{}, pixel diffs {diffs} (sound bound \
                 {sound_diffs}), ratio {:.1}-{:.1}\n",
                pages.iter().copied().min().unwrap_or(0),
                pages.iter().copied().max().unwrap_or(0),
                k.iter().copied().fold(f64::INFINITY, f64::min),
                k.iter().copied().fold(f64::NEG_INFINITY, f64::max),
                k.iter().filter(|v| **v < 3.0).count(),
                k.len(),
                ratios.iter().copied().fold(f64::INFINITY, f64::min),
                ratios.iter().copied().fold(f64::NEG_INFINITY, f64::max),
            ));
        }
    }

    // How far the stroke actually reached, per field. This is the row's
    // mechanism and it is page-size independent, so it gets its own table.
    println!(
        "reach of one stroke beyond its own bounding disc, per field \
         (derived sound dilation {SOUND_DILATE_CELLS:.4} cells = 2 cell diagonals):"
    );
    println!(
        "{:>15} {:>10} {:>10} {:>10} {:>11} {:>12} {:>12}",
        "field", "max cells", "diagonals", "out tight", "out dilated", "out sound", "diff sound"
    );
    for name in [
        "sphere",
        "torus",
        "box_exact",
        "csg_difference",
        "thin_plate",
        "gyroid",
        "fbm_terrain",
        "noise_cavity",
    ] {
        let sel: Vec<&&Row> = dig.iter().filter(|r| get(r, "field") == name).collect();
        println!(
            "{name:>15} {:>10.3} {:>10.3} {:>10} {:>11} {:>12} {:>12}",
            sel.iter()
                .map(|r| num(r, "max_change_excess_cells"))
                .fold(f64::NEG_INFINITY, f64::max),
            sel.iter()
                .map(|r| num(r, "excess_in_cell_diagonals"))
                .fold(f64::NEG_INFINITY, f64::max),
            sel.iter()
                .map(|r| count(r, "changed_pixels_outside_tight_disc"))
                .max()
                .unwrap_or(0),
            sel.iter()
                .map(|r| count(r, "changed_pixels_outside_dilated_disc"))
                .max()
                .unwrap_or(0),
            sel.iter()
                .map(|r| count(r, "changed_pixels_outside_sound_disc"))
                .max()
                .unwrap_or(0),
            sel.iter()
                .map(|r| count(r, "pixel_diff_sound"))
                .max()
                .unwrap_or(0),
        );
    }
    println!("\n{summary}");

    println!("buried control (the reachable zero the vacuity control needs):");
    for row in &buried {
        println!(
            "  {:>15} clearance {:>7} cells, changed_pixels_full {}",
            get(row, "field"),
            get(row, "buried_clearance_cells"),
            get(row, "changed_pixels_full"),
        );
    }
    println!();

    // The CSV is written before the controls fire, so a failed control leaves an
    // artefact to read rather than an empty directory.
    common::experiment::run(isomesh::experiment!("P-79"), |run| {
        for row in &rows {
            run.record(row);
        }
    });

    // ── the registered vacuity control ──
    for row in &dig {
        assert!(
            count(row, "changed_pixels_full") > 0,
            "P-79 VACUITY CONTROL: {} under {} with a {}-cell brush changed no shadow texel, \
             so the light cannot see the dig and every localisation number on this row is a \
             perfect answer over an empty change set",
            get(row, "field"),
            get(row, "light_direction"),
            get(row, "brush_radius_cells"),
        );
    }
    // ── and the proof that zero was reachable (M-44) ──
    let reachable = buried
        .iter()
        .filter(|r| count(r, "changed_pixels_full") == 0)
        .count();
    assert!(
        reachable > 0,
        "P-79: no buried brush produced changed_pixels_full == 0, so the non-zero control \
         above could not have failed and is not a measurement. Either the depth buffer is not \
         a first-hit buffer or no field has {} cells of interior clearance.",
        BRUSH_CELLS[0]
    );
    println!(
        "controls: {} dig rows all changed at least one texel; {}/{} buried rows changed none, \
         so the zero was reachable.",
        dig.len(),
        reachable,
        buried.len()
    );
}
