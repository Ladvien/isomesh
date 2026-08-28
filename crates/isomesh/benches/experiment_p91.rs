//! **P-91 — geomorph against dither, in units of pixels.**
//!
//! Ticket: R-091. Pre-registered before this harness existed; the registration
//! lives in `crates/isomesh/src/experiment.rs` and is not amended here.
//!
//! ```bash
//! cargo bench --bench experiment_p91
//! ```
//!
//! Writes `docs/experiments/p-91.csv`.
//!
//! # Sources, left in comments as required
//!
//! - Lengyel, *Voxel-Based Terrain for Real-Time Virtual Simulations* (PhD
//!   dissertation, UC Davis, 2010), ch. 6 / §6.1 — the one-scalar geomorph, with
//!   the guarantee explicitly withheld: *"(It is unclear whether this point
//!   always exists.)"* Also `10.1080/2151237X.2011.563682` (Transvoxel).
//! - Haydel, Yuksel & Seiler, `10.1145/3618359` — chose morphing over stochastic
//!   LOD because stochastic causes *"significantly increased data movement"*
//!   spikes at transitions. **Their hardware is a cycle-accurate simulation of
//!   TRaX, not a real GPU**, so what ports is the structural argument and not the
//!   number; `resident_bytes_ratio` below is this harness's measured form of it.
//! - Karis, *High-Quality Temporal Supersampling* (SIGGRAPH 2014 Advances in
//!   Real-Time Rendering) — the YCoCg neighbourhood clamp, which is the
//!   rejection predicate `P-77` registered and this harness reuses.
//! - Ikkala, Lauttia, Jääskeläinen & Mäkitalo, `10.1145/3681758.3697996` —
//!   `P-77`'s subject.
//! - Ericson, *Real-Time Collision Detection* (2005) §5.1.5 — closest point on a
//!   triangle, seven-region form.
//! - Möller & Trumbore, `10.1080/10867651.1997.10487468` — ray/triangle test.
//! - StopThePop `arXiv:2402.00525` and FLIP `10.1145/3406183` — the only
//!   published popping metrics, both image-space and both needing a trained
//!   perceptual model. **Item 3.8's caveat is repeated rather than dropped: the
//!   metric below is plausible, not perceptual.** It has not been validated
//!   against `FLIP_t` or a forced-choice study, and must not be quoted as if it
//!   had.
//! - Roberts 2018's R2 low-discrepancy sequence — the jitter pattern, via `P-77`.
//!
//! # Instrument 1 — the pop, on geometry
//!
//! `M-121`'s own construction, taken from
//! `bevy_isomesh/examples/game_lod_flyover.rs` rather than re-derived:
//! `FbmTerrain::<f32>::canonical()`, blocks 4.0 units wide along `x` and ±4.0 in
//! `y` and `z`, `BASE_H = 0.25` doubling per level, meshed with
//! [`MarchingCubes`] at both the old and the new level. `M-121`'s statistic is
//! `worst_gap(before, after) / spacing(min(was, now))` — a one-sided Hausdorff
//! distance over **vertices**, in cells of the **finer** level.
//! [`reproduce_m121`] runs it over all twelve blocks and both level pairs the
//! demo's `MAX_LEVEL = 2` allows, in both directions, and the result is
//! **asserted** against the committed **3.136**: `M-279`'s rule, that a new
//! instrument's first job is to agree with the old one where they overlap.
//!
//! **The pixel conversion is item 3.8's, verbatim:** at distance `d`, vertical
//! FOV `θ` and vertical resolution `H`, a world displacement `δ` subtends
//! `δ·H / (2·d·tan(θ/2))` pixels, with `H = 1080` and `θ = 45°` (Bevy's
//! `PerspectiveProjection::default`) — both on every row as columns rather than
//! only in this comment. That is the **isotropic** form, the maximum over
//! displacement directions, which is what 3.8 states; the across-view component
//! alone rides along as `p99_pixels_across_view`, because 3.8 notes that
//! across-view motion is the visible half.
//!
//! **`d` is the repo's own current switch distance, not a guess.**
//! `game_lod_flyover::level_for` is `|centre − at| / LEVEL_RANGE` truncated with
//! `LEVEL_RANGE = 7.0`, so level `L → L+1` switches at `7·(L+1)` world units and
//! each row's `switch_distance` is that.
//!
//! **Two displacement measures, because they differ by a factor that matters.**
//! `M-121`'s vertex-to-nearest-**vertex** distance contains the coarse mesh's own
//! vertex spacing: a fine vertex lying exactly on the coarse *surface* is still
//! up to half a coarse vertex spacing from the nearest coarse *vertex*. The
//! visible quantity is the distance to the coarse **surface**, so the registered
//! `p99_pixels_of_pop` is taken from the vertex-to-nearest-triangle distance
//! (Ericson §5.1.5, exact) and the `M-121`-denominated version rides along as
//! `p99_pixels_vertex_metric`. Both are on every row; the gap between them is a
//! correction to how a pop should be read, not a choice buried in a comment.
//!
//! # The three methods, as per-vertex residuals
//!
//! - **`none`** — the hard switch. Every vertex's residual is the full
//!   displacement. This is the registered vacuity control's arm.
//! - **`geomorph`** — item 3.6's harness, exactly: for every fine vertex, cast
//!   along its **stored normal** ([`MeshBuffer::normals`], the crate's own),
//!   intersect the coarse mesh **restricted to the containing coarse cell**, and
//!   classify (a) no intersection, (b) intersection whose interpolated normal has
//!   negative dot product with `n`, (c) success. A success morphs the vertex onto
//!   the coarse surface, so its residual is **zero**; a failure has no morph
//!   target, must snap, and keeps the full displacement.
//! - **`dither`** — alpha-to-coverage. A pixel's resolved surface is the
//!   coverage-weighted mixture of the two surfaces, so with `S` coverage slots a
//!   fade transfers coverage in steps of `1/S` and the worst single-frame
//!   resolved displacement is `δ/min(S, F)` over a fade of `F` frames. `S = 8`
//!   (8× MSAA alpha-to-coverage) drives the registered column; `S = 1` and `4`
//!   ride along, and `S = 1` is the honest statement that a dither with no
//!   sub-pixel coverage **is** a hard switch scattered across pixels.
//!
//! # Instrument 2 — the cost, in `P-77`'s rejection predicate
//!
//! A software TAA resolve over a ray-cast depth buffer, transcribed from
//! `crates/isomesh/benches/experiment_p77.rs`: the same R2 jitter, the same
//! bilinear history fetch requiring four valid taps, the same 3×3 YCoCg
//! neighbourhood **AABB** clip (Karis 2014), the same rejection definition — *the
//! clip moved the sample*, `s < 1.0` — the same `α = 0.1` blend, and the same
//! world-space albedo detail. `P-77` measured a **95.5%** steady-state rejection
//! rate on a smooth shade, and that is a fixture defect this harness inherits the
//! fix for rather than rediscovering.
//!
//! `P-77`'s second defect binds the camera: the rejection rate is a function of
//! reprojection displacement and a walking camera **saturates** it at 86.6%, so
//! the camera here is `P-77`'s `HEADLINE_ARM` regime — **static**, jitter on,
//! whose steady rate there was 0.49% — the only regime where a ratio has
//! headroom.
//!
//! What is traced is **the trilinear reconstruction of the field on each level's
//! own sample grid**, not the analytic field: that is the surface the extracted
//! mesh approximates, so both instruments look at the same pair of surfaces. A
//! trilinear interpolant is a convex combination of its eight cell corners, so a
//! cell whose corners share a sign **cannot** contain a crossing — the
//! traversal's rejection test is exact and no Lipschitz constant is declared or
//! needed.
//!
//! # The SHARE line, recomputed before the code was written
//!
//! The registration's own SHARE is *"C1 and C2 move the LOD transition only; C3
//! moves the TAA resolve."* Those are **disjoint budgets**, so nothing here is
//! double-counted: the pop clauses move geometry the transition emits and C3
//! moves history the resolve reuses. Neither clause is a fraction of a total, so
//! there is no `✗51`-style ceiling from a share — but two of the three have hard
//! arithmetic bounds anyway, and both are worth stating before the run.
//!
//! **C1's terrain half is a 1% test on item 3.6's failure rate, and nothing
//! else.** The geomorph residual constructed below is **exactly zero** where the
//! morph target exists (the target is a point *of* the coarse mesh, so the
//! displayed vertex is on the surface it is handing over to) and **the full
//! displacement** where it does not (no target, so the vertex snaps). A p99 over
//! that two-valued distribution is under one pixel only if the failure population
//! is under 1%. The alternative route — failures whose own displacement happens
//! to be sub-pixel — is arithmetically dead: at `d = 7.0`, `H = 1080` and
//! `θ = 45°`, one pixel is **0.00537 world units = 0.0215 fine cells**, and
//! `M-121` already measured the pop at **0.6–3.14 cells**. So C1's terrain half
//! is a direct test of *"does Lengyel's containing-cell rule succeed on more than
//! 99% of fine vertices"*, which is what item 3.6's *"failure rate ≈ 0%"*
//! implicitly claims. `geomorph_fail_fraction_cell` is therefore the column the
//! clause really turns on, and it is on every row.
//!
//! **C3's superadditivity has a ceiling made of the brush's screen footprint.**
//! The combined arm's rejections lie inside the union of the two causes, so
//! `superadditivity ≤ 1 + excess_dig / (excess_lod + excess_dig)`: the dither's
//! flip set is the whole frame and the dig's changed set is a 0.25-radius
//! silhouette at ~3 units, so the dig excess is small and the clause is reachable
//! only just. `superadditivity_ceiling` is on the row so the verdict can be read
//! against what the fixture allowed rather than against 1.0 alone, and
//! `superadditivity_excess_samples` gives the interaction in samples so a reader
//! can see whether a "held" is a mechanism or a rounding.
//!
//! **C2 has no arithmetic bound and needed none.** Its pop half compares two
//! measured distributions and its cost half compares two measured counts; both
//! were reachable in either direction and both are reported.

#![allow(
    clippy::cast_possible_truncation,
    clippy::cast_precision_loss,
    clippy::cast_sign_loss,
    clippy::cast_possible_wrap,
    clippy::needless_range_loop,
    clippy::too_many_arguments,
    clippy::too_many_lines,
    clippy::print_literal,
    dead_code
)]

mod common;

use isomesh::fields::{FbmTerrain, Gyroid};
use isomesh::marching_cubes::MarchingCubes;
use isomesh::{MeshBuffer, RuntimeShape3, Sdf};

// ─── game_lod_flyover's geometry, restated with the source named ────────────
//
// Every constant here is copied from
// `bevy_isomesh/examples/game_lod_flyover.rs`, which is where `M-121` was
// measured. Changing any of them stops this harness reproducing that row.

/// `game_lod_flyover::BASE_H`. Level `L`'s spacing is `BASE_H · 2^L`.
const BASE_H: f32 = 0.25;
/// `game_lod_flyover::BLOCK_W` — one block's extent along `x`.
const BLOCK_W: f32 = 4.0;
/// `game_lod_flyover::CROSS` — half-extent in `y` and `z`.
const CROSS: f32 = 4.0;
/// `game_lod_flyover::BLOCKS`.
const BLOCKS: usize = 12;
/// `game_lod_flyover::LEVEL_RANGE`. `level_for` is `|centre − at| / LEVEL_RANGE`
/// truncated, so level `L → L+1` switches at `LEVEL_RANGE · (L + 1)`.
const LEVEL_RANGE: f32 = 7.0;
/// Level pairs measured: `L → L+1` for these `L`.
///
/// `game_lod_flyover::MAX_LEVEL` is 2, so the demo's own switches are `0→1` and
/// `1→2`. `2→3` is measured as well, because item 3.8's question is *at what
/// distance*, and stopping at the demo's ceiling would truncate the answer.
const PAIRS: [u32; 3] = [0, 1, 2];

/// `M-121`'s committed worst pop, in cells of the finer level.
const M121_WORST_POP_CELLS: f64 = 3.136;
/// Relative slack **above** `M-121`, and the bound is one-sided for a reason.
///
/// [`reproduce_m121`] enumerates every ordered level pair on all twelve blocks,
/// which is a **superset** of the configurations the demo's flight actually
/// visited — the flight's sequence depends on a user-controlled fly speed, so it
/// can skip pairs this sweep cannot. A superset can only push the maximum up, so
/// the honest check is `M-121 <= worst <= M-121 · (1 + tol)`: below is a
/// different measurement, far above is a different fixture. Measured **3.141**
/// against the committed **3.136**, i.e. +0.16%.
const M121_TOLERANCE: f64 = 0.01;

// ─── the camera model, item 3.8's ───────────────────────────────────────────

/// Vertical resolution the pixel conversion is denominated in. On every row.
const VIEW_HEIGHT_PX: f64 = 1080.0;
/// Vertical field of view, radians: Bevy's `PerspectiveProjection::default`.
const FOV_Y: f64 = core::f64::consts::FRAC_PI_4;

/// Item 3.8's conversion, verbatim: `δ·H / (2·d·tan(θ/2))`.
fn pixels_of(delta_world: f64, distance: f64) -> f64 {
    delta_world * VIEW_HEIGHT_PX / (2.0 * distance * (FOV_Y * 0.5).tan())
}

// ─── the three methods ──────────────────────────────────────────────────────

/// Coverage slots per pixel swept for the dither arm.
const DITHER_SLOTS: [u32; 3] = [1, 4, 8];
/// The slot count the registered columns are taken from: 8× MSAA
/// alpha-to-coverage.
const DITHER_SLOTS_HEADLINE: u32 = 8;
/// Frames a transition is spread over: 12 at 60 Hz is 0.2 s, a typical LOD fade.
/// `F ≥ S`, so the coverage quantisation and not the fade length sets dither's
/// worst single-frame step.
const FADE_FRAMES: u32 = 12;

/// A full secondary position: three `f32`.
const BYTES_SECOND_POSITION: u32 = 12;
/// Lengyel's one scalar along a normal the vertex already carries: one `f32`.
///
/// **The docs' arithmetic does not come out at 6 against 12, and this column says
/// so rather than reprinting it.** Item 3.6 says *"store a distance along the
/// normal instead of a full second position, saving 12 bytes per vertex"*; the
/// gameplay doc renders that as *"6 bytes/vertex instead of 12"*. A full
/// secondary position is 12 bytes and one `f32` scalar is 4, so the saving is
/// **3×**, not the 2× the prose implies — stronger than claimed, in the direction
/// that favours geomorph. `bytes_per_vertex_doc_claim` carries the 6 so a reader
/// can redo the arithmetic.
const BYTES_ONE_SCALAR: u32 = 4;
/// The figure the gameplay doc states, carried for comparison only.
const BYTES_DOC_CLAIM: u32 = 6;

// ─── small vector arithmetic ────────────────────────────────────────────────

type V3 = [f32; 3];

fn sub3(a: V3, b: V3) -> V3 {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
fn add3(a: V3, b: V3) -> V3 {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}
fn mul3(a: V3, s: f32) -> V3 {
    [a[0] * s, a[1] * s, a[2] * s]
}
fn dot3(a: V3, b: V3) -> f32 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
fn cross3(a: V3, b: V3) -> V3 {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
fn len3(a: V3) -> f32 {
    dot3(a, a).sqrt()
}

/// Unit vector, or `None` when there is no direction to return.
///
/// `None` rather than a substituted axis: a zero-length normal is an absence of
/// information, and inventing one would put a second execution path under every
/// geomorph verdict in this file.
fn unit3(a: V3) -> Option<V3> {
    let l = len3(a);
    if l > 0.0 { Some(mul3(a, 1.0 / l)) } else { None }
}

// ─── the ladder ─────────────────────────────────────────────────────────────

/// `game_lod_flyover::spacing`.
fn spacing(level: u32) -> f32 {
    BASE_H * (1 << level) as f32
}

/// Cells per axis of one block at `level`: `[along, across, across]`.
fn block_cells(level: u32) -> [u32; 3] {
    let h = spacing(level);
    let along = (BLOCK_W / h).round() as u32;
    let across = ((CROSS * 2.0) / h).round() as u32;
    [along, across, across]
}

/// `game_lod_flyover::mesh_block`, with the origin passed in so one function
/// serves the twelve-block `M-121` reproduction and the centred `gyroid` block.
fn mesh_block<F: Sdf<Scalar = f32>>(field: &F, origin: V3, level: u32) -> MeshBuffer<f32> {
    let c = block_cells(level);
    let shape = RuntimeShape3::new([c[0] + 1, c[1] + 1, c[2] + 1]).expect("valid shape");
    let mut out = MeshBuffer::<f32>::new();
    MarchingCubes::<f32>::new()
        .extract(field, &shape, origin, spacing(level), &mut out)
        .expect("extraction");
    out
}

// ─── nearest-vertex and nearest-triangle queries ────────────────────────────

/// A uniform grid over points, for nearest-vertex queries.
struct VertGrid {
    lo: V3,
    inv_cell: f32,
    dims: [usize; 3],
    cells: Vec<Vec<u32>>,
    points: Vec<V3>,
}

fn grid_bounds(points: &[V3], cell: f32) -> (V3, [usize; 3]) {
    let mut lo = [f32::INFINITY; 3];
    let mut hi = [f32::NEG_INFINITY; 3];
    for p in points {
        for k in 0..3 {
            lo[k] = lo[k].min(p[k]);
            hi[k] = hi[k].max(p[k]);
        }
    }
    let mut dims = [1usize; 3];
    for k in 0..3 {
        dims[k] = (((hi[k] - lo[k]) / cell).ceil() as usize + 1).max(1);
    }
    (lo, dims)
}

fn cell_index(lo: V3, inv_cell: f32, dims: [usize; 3], p: V3) -> [usize; 3] {
    let mut out = [0usize; 3];
    for k in 0..3 {
        let f = ((p[k] - lo[k]) * inv_cell).floor();
        out[k] = (f.max(0.0) as usize).min(dims[k] - 1);
    }
    out
}

impl VertGrid {
    fn new(points: &[V3], cell: f32) -> Self {
        let (lo, dims) = grid_bounds(points, cell);
        let inv_cell = 1.0 / cell;
        let mut cells = vec![Vec::new(); dims[0] * dims[1] * dims[2]];
        for (i, p) in points.iter().enumerate() {
            let c = cell_index(lo, inv_cell, dims, *p);
            cells[c[0] + dims[0] * (c[1] + dims[1] * c[2])].push(i as u32);
        }
        Self {
            lo,
            inv_cell,
            dims,
            cells,
            points: points.to_vec(),
        }
    }

    /// Exact distance to the nearest point, by expanding shells until the shell's
    /// own lower bound exceeds the best found.
    fn nearest(&self, p: V3) -> f32 {
        let c = cell_index(self.lo, self.inv_cell, self.dims, p);
        let cell = 1.0 / self.inv_cell;
        let mut best = f32::INFINITY;
        let max_ring = self.dims.iter().copied().max().unwrap_or(1) as i64;
        for ring in 0..=max_ring {
            if ring > 0 && ((ring - 1) as f32) * cell > best {
                break;
            }
            for_shell(c, ring, self.dims, |idx| {
                for &j in &self.cells[idx] {
                    let d = len3(sub3(self.points[j as usize], p));
                    if d < best {
                        best = d;
                    }
                }
            });
        }
        best
    }
}

/// Visit the flat indices of every in-range cell exactly `ring` cells away in
/// the Chebyshev metric.
fn for_shell(c: [usize; 3], ring: i64, dims: [usize; 3], mut f: impl FnMut(usize)) {
    for dz in -ring..=ring {
        for dy in -ring..=ring {
            for dx in -ring..=ring {
                if dx.abs().max(dy.abs()).max(dz.abs()) != ring {
                    continue;
                }
                let (ix, iy, iz) = (c[0] as i64 + dx, c[1] as i64 + dy, c[2] as i64 + dz);
                if ix < 0
                    || iy < 0
                    || iz < 0
                    || ix >= dims[0] as i64
                    || iy >= dims[1] as i64
                    || iz >= dims[2] as i64
                {
                    continue;
                }
                f(ix as usize + dims[0] * (iy as usize + dims[1] * iz as usize));
            }
        }
    }
}

/// Closest point on a triangle, with its barycentric coordinates.
///
/// Ericson, *Real-Time Collision Detection* (2005) §5.1.5 — the seven-region
/// form, so it is exact on edges and at vertices rather than solving a
/// least-squares system that is singular there.
fn closest_on_tri(p: V3, a: V3, b: V3, c: V3) -> (V3, [f32; 3]) {
    let ab = sub3(b, a);
    let ac = sub3(c, a);
    let ap = sub3(p, a);
    let d1 = dot3(ab, ap);
    let d2 = dot3(ac, ap);
    if d1 <= 0.0 && d2 <= 0.0 {
        return (a, [1.0, 0.0, 0.0]);
    }
    let bp = sub3(p, b);
    let d3 = dot3(ab, bp);
    let d4 = dot3(ac, bp);
    if d3 >= 0.0 && d4 <= d3 {
        return (b, [0.0, 1.0, 0.0]);
    }
    let vc = d1 * d4 - d3 * d2;
    if vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0 {
        let t = d1 / (d1 - d3);
        return (add3(a, mul3(ab, t)), [1.0 - t, t, 0.0]);
    }
    let cp = sub3(p, c);
    let d5 = dot3(ab, cp);
    let d6 = dot3(ac, cp);
    if d6 >= 0.0 && d5 <= d6 {
        return (c, [0.0, 0.0, 1.0]);
    }
    let vb = d5 * d2 - d1 * d6;
    if vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0 {
        let t = d2 / (d2 - d6);
        return (add3(a, mul3(ac, t)), [1.0 - t, 0.0, t]);
    }
    let va = d3 * d6 - d5 * d4;
    if va <= 0.0 && (d4 - d3) >= 0.0 && (d5 - d6) >= 0.0 {
        let t = (d4 - d3) / ((d4 - d3) + (d5 - d6));
        return (add3(b, mul3(sub3(c, b), t)), [0.0, 1.0 - t, t]);
    }
    let inv = 1.0 / (va + vb + vc);
    let v = vb * inv;
    let w = vc * inv;
    (add3(a, add3(mul3(ab, v), mul3(ac, w))), [1.0 - v - w, v, w])
}

/// A uniform grid over a mesh's triangles: nearest-point and box-restricted
/// queries.
struct TriGrid {
    lo: V3,
    inv_cell: f32,
    dims: [usize; 3],
    cells: Vec<Vec<u32>>,
    tris: Vec<[u32; 3]>,
    positions: Vec<V3>,
    normals: Vec<V3>,
    /// Triangles dropped for zero area, counted rather than silently skipped:
    /// `closest_on_tri` divides by `va + vb + vc`, zero for a degenerate.
    degenerate: usize,
}

impl TriGrid {
    fn new(mesh: &MeshBuffer<f32>, cell: f32) -> Self {
        let (lo, dims) = grid_bounds(&mesh.positions, cell);
        let inv_cell = 1.0 / cell;
        let mut cells = vec![Vec::new(); dims[0] * dims[1] * dims[2]];
        let mut tris = Vec::new();
        let mut degenerate = 0usize;
        for t in mesh.indices.as_chunks::<3>().0 {
            let (a, b, c) = (
                mesh.positions[t[0] as usize],
                mesh.positions[t[1] as usize],
                mesh.positions[t[2] as usize],
            );
            if len3(cross3(sub3(b, a), sub3(c, a))) <= 0.0 {
                degenerate += 1;
                continue;
            }
            tris.push([t[0], t[1], t[2]]);
        }
        for (i, t) in tris.iter().enumerate() {
            let mut tlo = [f32::INFINITY; 3];
            let mut thi = [f32::NEG_INFINITY; 3];
            for &vi in t {
                let p = mesh.positions[vi as usize];
                for k in 0..3 {
                    tlo[k] = tlo[k].min(p[k]);
                    thi[k] = thi[k].max(p[k]);
                }
            }
            let c0 = cell_index(lo, inv_cell, dims, tlo);
            let c1 = cell_index(lo, inv_cell, dims, thi);
            for z in c0[2]..=c1[2] {
                for y in c0[1]..=c1[1] {
                    for x in c0[0]..=c1[0] {
                        cells[x + dims[0] * (y + dims[1] * z)].push(i as u32);
                    }
                }
            }
        }
        Self {
            lo,
            inv_cell,
            dims,
            cells,
            tris,
            positions: mesh.positions.clone(),
            normals: mesh.normals.clone(),
            degenerate,
        }
    }

    fn tri_points(&self, j: u32) -> (V3, V3, V3) {
        let t = self.tris[j as usize];
        (
            self.positions[t[0] as usize],
            self.positions[t[1] as usize],
            self.positions[t[2] as usize],
        )
    }

    /// Exact nearest point on the mesh, and the offset to it.
    fn nearest(&self, p: V3) -> (f32, V3) {
        let c = cell_index(self.lo, self.inv_cell, self.dims, p);
        let cell = 1.0 / self.inv_cell;
        let mut best = (f32::INFINITY, [0.0f32; 3]);
        let max_ring = self.dims.iter().copied().max().unwrap_or(1) as i64;
        for ring in 0..=max_ring {
            if ring > 0 && ((ring - 1) as f32) * cell > best.0 {
                break;
            }
            for_shell(c, ring, self.dims, |idx| {
                for &j in &self.cells[idx] {
                    let (a, b, cc) = self.tri_points(j);
                    let (q, _) = closest_on_tri(p, a, b, cc);
                    let off = sub3(q, p);
                    let d = len3(off);
                    if d < best.0 {
                        best = (d, off);
                    }
                }
            });
        }
        best
    }

    /// Triangles overlapping an axis-aligned box, de-duplicated.
    fn tris_in_box(&self, blo: V3, bhi: V3, out: &mut Vec<u32>) {
        out.clear();
        let c0 = cell_index(self.lo, self.inv_cell, self.dims, blo);
        let c1 = cell_index(self.lo, self.inv_cell, self.dims, bhi);
        for z in c0[2]..=c1[2] {
            for y in c0[1]..=c1[1] {
                for x in c0[0]..=c1[0] {
                    out.extend_from_slice(&self.cells[x + self.dims[0] * (y + self.dims[1] * z)]);
                }
            }
        }
        out.sort_unstable();
        out.dedup();
    }
}

/// Möller–Trumbore, `10.1080/10867651.1997.10487468`. Signed `t`: this is a full
/// line, because Lengyel's morph may go either way along the normal.
fn ray_tri(origin: V3, dir: V3, a: V3, b: V3, c: V3) -> Option<(f32, [f32; 3])> {
    let e1 = sub3(b, a);
    let e2 = sub3(c, a);
    let pv = cross3(dir, e2);
    let det = dot3(e1, pv);
    if det.abs() < 1e-12 {
        return None;
    }
    let inv = 1.0 / det;
    let tv = sub3(origin, a);
    let u = dot3(tv, pv) * inv;
    if u < 0.0 || u > 1.0 {
        return None;
    }
    let qv = cross3(tv, e1);
    let v = dot3(dir, qv) * inv;
    if v < 0.0 || u + v > 1.0 {
        return None;
    }
    Some((dot3(e2, qv) * inv, [1.0 - u - v, u, v]))
}

// ─── the geomorph predicate, item 3.6's harness ─────────────────────────────

/// How far the *unrestricted* attribution arm looks along the normal, in coarse
/// cells either way.
///
/// Three, and it is a bound rather than a taste: the largest surface
/// displacement this harness measures anywhere is under one coarse cell, so a
/// ray that finds nothing within three has not been cut short by the reach.
/// `morph_t_max_coarse_cells_unrestricted` reports the largest `|t|` actually
/// used, which is the column that proves the reach was not the binding
/// constraint.
const MORPH_REACH_COARSE_CELLS: f32 = 3.0;

/// One classification of one fine vertex under one restriction rule.
#[derive(Clone, Copy, Default, Debug)]
struct Geomorph {
    /// (c) success: an intersection inside the allowed region whose interpolated
    /// coarse normal agrees with the fine normal.
    success: u64,
    /// (a) no intersection with the coarse mesh inside the allowed region.
    no_hit: u64,
    /// (b) intersection with non-positive normal dot product.
    flipped: u64,
    /// Fine vertices whose stored normal has zero length, so the ray has no
    /// direction. Counted rather than substituted.
    no_normal: u64,
    /// Largest `|t|` over the successes, in coarse cells — the scalar Lengyel
    /// stores, so its range is what a quantisation would have to cover, and for
    /// the unrestricted arm it is also the proof that the reach did not bind.
    morph_t_max_coarse_cells: f64,
}

impl Geomorph {
    fn total(&self) -> u64 {
        self.success + self.no_hit + self.flipped + self.no_normal
    }
    fn fail_fraction(&self) -> f64 {
        let t = self.total();
        if t == 0 {
            0.0
        } else {
            (self.no_hit + self.flipped + self.no_normal) as f64 / t as f64
        }
    }
}

/// Three restriction rules on the same rays, so *"the morph failed"* names one
/// thing rather than three.
///
/// Item 3.6's harness is the **containing coarse cell**, and that is the one the
/// registered columns come from. The other two exist because the first run of
/// this harness read a 40–60% failure rate on *both* fields, and a number that
/// large is either a finding or a fixture defect — and the only way to tell is to
/// ask whether the ray misses the coarse surface or merely meets it next door.
/// `cell` ⊆ `cell27` ⊆ `reach`, so the three fractions are monotone by
/// construction and any inversion is a bug in this function.
struct GeomorphArms {
    /// Item 3.6, verbatim: the intersection must lie in the coarse cell
    /// containing the fine vertex. Registered.
    cell: Geomorph,
    /// The containing cell dilated to its 3×3×3 neighbourhood — the cheapest
    /// relaxation an implementation could actually ship, since a morph target one
    /// cell over is still local.
    cell27: Geomorph,
    /// No cell restriction at all: the nearest intersection within
    /// [`MORPH_REACH_COARSE_CELLS`]. This is the arm that tests 3.6's *mechanism*
    /// — *"a heightfield is star-shaped along the up axis and a triply-periodic
    /// minimal surface is not"* — because star-shapedness is a statement about
    /// the whole ray and not about one cell.
    reach: Geomorph,
}

/// For each fine vertex: cast along its stored normal, intersect the coarse mesh
/// under each of the three restriction rules, and return the per-vertex geomorph
/// residual under the registered rule — zero where the morph lands, the full
/// surface displacement where it does not.
fn geomorph_scan(
    fine: &MeshBuffer<f32>,
    coarse: &TriGrid,
    coarse_origin: V3,
    coarse_h: f32,
    surface_delta: &[f32],
) -> (GeomorphArms, Vec<f32>) {
    let mut arms = GeomorphArms {
        cell: Geomorph::default(),
        cell27: Geomorph::default(),
        reach: Geomorph::default(),
    };
    let mut residual = Vec::with_capacity(fine.positions.len());
    let mut candidates: Vec<u32> = Vec::new();
    let reach = MORPH_REACH_COARSE_CELLS * coarse_h;
    for (i, &v) in fine.positions.iter().enumerate() {
        let Some(n) = unit3(fine.normals[i]) else {
            arms.cell.no_normal += 1;
            arms.cell27.no_normal += 1;
            arms.reach.no_normal += 1;
            residual.push(surface_delta[i]);
            continue;
        };
        // The containing coarse cell, as a box on the coarse lattice, and its
        // 3×3×3 dilation.
        let mut clo = [0.0f32; 3];
        let mut chi = [0.0f32; 3];
        for k in 0..3 {
            let idx = ((v[k] - coarse_origin[k]) / coarse_h).floor();
            clo[k] = coarse_origin[k] + idx * coarse_h;
            chi[k] = clo[k] + coarse_h;
        }
        let dlo = [clo[0] - coarse_h, clo[1] - coarse_h, clo[2] - coarse_h];
        let dhi = [chi[0] + coarse_h, chi[1] + coarse_h, chi[2] + coarse_h];
        // One candidate set, from the largest region any arm can accept, so all
        // three arms see the same triangles and differ only in the containment
        // test. Anything else would let a grid-query difference masquerade as a
        // restriction difference.
        let qlo = [v[0] - reach, v[1] - reach, v[2] - reach];
        let qhi = [v[0] + reach, v[1] + reach, v[2] + reach];
        coarse.tris_in_box(qlo, qhi, &mut candidates);
        // Nearest accepted hit under each rule: (|t|, normal dot product).
        let mut best_cell: Option<(f32, f32)> = None;
        let mut best_27: Option<(f32, f32)> = None;
        let mut best_reach: Option<(f32, f32)> = None;
        for &j in &candidates {
            let (a, b, c) = coarse.tri_points(j);
            let Some((tt, bary)) = ray_tri(v, n, a, b, c) else {
                continue;
            };
            if tt.abs() > reach {
                continue;
            }
            let hit = add3(v, mul3(n, tt));
            let t3 = coarse.tris[j as usize];
            let mut ni = [0.0f32; 3];
            for (w, vi) in bary.iter().zip(t3.iter()) {
                ni = add3(ni, mul3(coarse.normals[*vi as usize], *w));
            }
            let cand = (tt.abs(), dot3(ni, n));
            let closer = |b: &Option<(f32, f32)>| b.is_none_or(|(bt, _)| cand.0 < bt);
            // "Restricted to the containing cell" is a containment test on the
            // intersection *point*: a triangle may straddle a cell wall and be
            // met outside the cell it is listed in.
            let inside = |lo: V3, hi: V3| {
                (0..3).all(|k| hit[k] >= lo[k] - 1e-5 && hit[k] <= hi[k] + 1e-5)
            };
            if inside(clo, chi) && closer(&best_cell) {
                best_cell = Some(cand);
            }
            if inside(dlo, dhi) && closer(&best_27) {
                best_27 = Some(cand);
            }
            if closer(&best_reach) {
                best_reach = Some(cand);
            }
        }
        for (arm, best) in [
            (&mut arms.cell, best_cell),
            (&mut arms.cell27, best_27),
            (&mut arms.reach, best_reach),
        ] {
            match best {
                None => arm.no_hit += 1,
                Some((_, d)) if d <= 0.0 => arm.flipped += 1,
                Some((t, _)) => {
                    arm.success += 1;
                    let cells = f64::from(t) / f64::from(coarse_h);
                    if cells > arm.morph_t_max_coarse_cells {
                        arm.morph_t_max_coarse_cells = cells;
                    }
                }
            }
        }
        // The residual is the registered rule's: zero where item 3.6's morph
        // lands, the full surface displacement where it does not.
        residual.push(match best_cell {
            Some((_, d)) if d > 0.0 => 0.0,
            _ => surface_delta[i],
        });
    }
    debug_assert!(
        arms.cell.fail_fraction() >= arms.cell27.fail_fraction()
            && arms.cell27.fail_fraction() >= arms.reach.fail_fraction(),
        "the three restriction rules are nested, so their failure fractions must be monotone"
    );
    (arms, residual)
}

// ─── statistics ─────────────────────────────────────────────────────────────

/// Nearest-rank percentile on a sorted slice.
fn percentile(sorted: &[f64], q: f64) -> f64 {
    if sorted.is_empty() {
        return 0.0;
    }
    let rank = (q * (sorted.len() - 1) as f64).round() as usize;
    sorted[rank.min(sorted.len() - 1)]
}

fn sorted_f64(v: &[f32]) -> Vec<f64> {
    let mut out: Vec<f64> = v.iter().map(|&x| f64::from(x)).collect();
    out.sort_by(f64::total_cmp);
    out
}

fn max_f32(v: &[f32]) -> f64 {
    f64::from(v.iter().fold(0.0f32, |a, &b| a.max(b)))
}

// ─── the pop, per field per level pair ──────────────────────────────────────

/// Everything the geometry instrument produces for one field and one level pair.
struct Pop {
    field: &'static str,
    level: u32,
    switch_distance: f64,
    fine_h: f64,
    coarse_h: f64,
    origin: V3,
    fine_vertices: usize,
    coarse_vertices: usize,
    fine_triangles: usize,
    coarse_triangles: usize,
    /// `M-121`'s statistic on this pair: worst vertex-to-nearest-vertex, in
    /// cells of the finer level.
    pop_cells_vertex_metric: f64,
    /// Same, but vertex-to-nearest-**triangle**: the visible displacement.
    pop_cells_surface: f64,
    /// Per-vertex displacements and residuals, world units, sorted.
    surface_sorted: Vec<f64>,
    vertex_sorted: Vec<f64>,
    geomorph_sorted: Vec<f64>,
    across_sorted: Vec<f64>,
    geomorph: GeomorphArms,
    degenerate_coarse: usize,
    fine_bytes: u64,
    coarse_bytes: u64,
}

/// Camera position for a row: at `switch_distance` from the block centre along
/// `+z`. The block is 8 units across in `z`, so a camera on that axis has the
/// whole surface at roughly one range and a single `switch_distance` describes
/// the row.
fn eye_for(centre: V3, distance: f64) -> V3 {
    [centre[0], centre[1], centre[2] + distance as f32]
}

fn block_centre(origin: V3) -> V3 {
    [
        origin[0] + BLOCK_W * 0.5,
        origin[1] + CROSS,
        origin[2] + CROSS,
    ]
}

fn measure_pop<F: Sdf<Scalar = f32>>(
    field_name: &'static str,
    field: &F,
    origin: V3,
    level: u32,
) -> Pop {
    let fine = mesh_block(field, origin, level);
    let coarse = mesh_block(field, origin, level + 1);
    assert!(
        !fine.positions.is_empty() && !coarse.positions.is_empty(),
        "{field_name} level {level}: an empty mesh on one side of the switch would make every \
         displacement below a measurement of nothing"
    );
    let fine_h = spacing(level);
    let coarse_h = spacing(level + 1);
    let switch_distance = f64::from(LEVEL_RANGE) * f64::from(level + 1);

    let vgrid = VertGrid::new(&coarse.positions, coarse_h);
    let tgrid = TriGrid::new(&coarse, coarse_h);
    let eye = eye_for(block_centre(origin), switch_distance);

    let mut vertex_delta = Vec::with_capacity(fine.positions.len());
    let mut surface_delta = Vec::with_capacity(fine.positions.len());
    let mut across = Vec::with_capacity(fine.positions.len());
    for &v in &fine.positions {
        vertex_delta.push(vgrid.nearest(v));
        let (d, off) = tgrid.nearest(v);
        surface_delta.push(d);
        // Across-view component: the part of the displacement perpendicular to
        // the view ray, which is the half item 3.8 calls the visible one.
        let view = unit3(sub3(v, eye)).expect("the eye is never on the surface");
        across.push(len3(sub3(off, mul3(view, dot3(off, view)))));
    }

    let (geomorph, geomorph_residual) =
        geomorph_scan(&fine, &tgrid, origin, coarse_h, &surface_delta);
    for (rule, arm) in [
        ("cell", &geomorph.cell),
        ("cell27", &geomorph.cell27),
        ("reach", &geomorph.reach),
    ] {
        assert_eq!(
            arm.total(),
            fine.positions.len() as u64,
            "{field_name} level {level} rule {rule}: the geomorph classification must account \
             for every fine vertex, or a failure class is being dropped rather than counted"
        );
    }
    assert!(
        geomorph.cell.fail_fraction() >= geomorph.cell27.fail_fraction()
            && geomorph.cell27.fail_fraction() >= geomorph.reach.fail_fraction(),
        "{field_name} level {level}: the three restriction rules are nested \
         (cell ⊆ cell27 ⊆ reach), so their failure fractions must be monotone; \
         {:.5} / {:.5} / {:.5} is not",
        geomorph.cell.fail_fraction(),
        geomorph.cell27.fail_fraction(),
        geomorph.reach.fail_fraction(),
    );

    let bytes = |m: &MeshBuffer<f32>| -> u64 {
        (m.positions.len() * 12 + m.normals.len() * 12 + m.indices.len() * 4) as u64
    };

    Pop {
        field: field_name,
        level,
        switch_distance,
        fine_h: f64::from(fine_h),
        coarse_h: f64::from(coarse_h),
        origin,
        fine_vertices: fine.positions.len(),
        coarse_vertices: coarse.positions.len(),
        fine_triangles: fine.indices.len() / 3,
        coarse_triangles: coarse.indices.len() / 3,
        pop_cells_vertex_metric: max_f32(&vertex_delta) / f64::from(fine_h),
        pop_cells_surface: max_f32(&surface_delta) / f64::from(fine_h),
        surface_sorted: sorted_f64(&surface_delta),
        vertex_sorted: sorted_f64(&vertex_delta),
        geomorph_sorted: sorted_f64(&geomorph_residual),
        across_sorted: sorted_f64(&across),
        geomorph,
        degenerate_coarse: tgrid.degenerate,
        fine_bytes: bytes(&fine),
        coarse_bytes: bytes(&coarse),
    }
}

// ─── M-121's reproduction ───────────────────────────────────────────────────

/// Worst distance from a vertex of `a` to the nearest vertex of `b`.
///
/// `game_lod_flyover::worst_gap`, transcribed. Brute force there, brute force
/// here: the point of a reproduction is to run the same arithmetic.
fn worst_gap(a: &MeshBuffer<f32>, b: &MeshBuffer<f32>) -> f32 {
    let mut worst = 0.0f32;
    for p in &a.positions {
        let mut nearest = f32::MAX;
        for q in &b.positions {
            let d = [p[0] - q[0], p[1] - q[1], p[2] - q[2]];
            let dist = d[0] * d[0] + d[1] * d[1] + d[2] * d[2];
            if dist < nearest {
                nearest = dist;
            }
        }
        if nearest.is_finite() && nearest > worst {
            worst = nearest;
        }
    }
    worst.sqrt()
}

/// `M-121` over the demo's whole ladder.
///
/// **Every ordered level pair, not only the adjacent ones, and that is what made
/// the reproduction work.** Restricted to `0↔1` and `1↔2` this reads a worst of
/// **1.618** cells against the committed **3.136** — the first run of this
/// harness did exactly that and disagreed with `M-121` by 1.94×. The demo's
/// levels come from `level_for(centre, at) = ⌊|centre − at| / 7.0⌋` recomputed
/// every frame with a **user-controlled fly speed** (`[` and `]`), so a block can
/// cross two `LEVEL_RANGE` boundaries between frames and switch `0 → 2` in one
/// step. `smooth` bounds the difference between *neighbouring blocks*, not the
/// difference between one block's successive levels. `h = spacing(min(was, now))`
/// then divides a two-level displacement by the *level-0* cell, which is where
/// the factor of two lives.
struct M121 {
    worst_cells: f64,
    /// Which ordered pair produced the worst, so the row can be read without
    /// re-running.
    worst_pair: (u32, u32),
    /// Per-switch values, sorted, so the *"typically 0.6–1.6"* half of the row
    /// can be checked as well as the worst.
    per_switch: Vec<f64>,
    /// The worst restricted to adjacent pairs, i.e. what a reader who assumed
    /// one-level switches would have measured.
    worst_adjacent_cells: f64,
}

fn reproduce_m121() -> M121 {
    let field = FbmTerrain::<f32>::canonical();
    let threads = std::thread::available_parallelism().map_or(8, |n| n.get());
    let pairs: Vec<(u32, u32)> = (0..=2u32)
        .flat_map(|a| (0..=2u32).map(move |b| (a, b)))
        .filter(|(a, b)| a != b)
        .collect();
    let work: Vec<(usize, u32, u32)> = (0..BLOCKS)
        .flat_map(|i| pairs.iter().map(move |&(a, b)| (i, a, b)))
        .collect();
    let chunk = work.len().div_ceil(threads);
    let out: Vec<Vec<(f64, u32, u32)>> = std::thread::scope(|s| {
        let handles: Vec<_> = work
            .chunks(chunk)
            .map(|slice| {
                s.spawn(move || {
                    let mut local = Vec::new();
                    for &(i, was, now) in slice {
                        let origin = [i as f32 * BLOCK_W, -CROSS, -CROSS];
                        let before = mesh_block(&field, origin, was);
                        let after = mesh_block(&field, origin, now);
                        let h = spacing(was.min(now));
                        let moved = f64::from(worst_gap(&before, &after) / h);
                        if moved.is_finite() {
                            local.push((moved, was, now));
                        }
                    }
                    local
                })
            })
            .collect();
        handles
            .into_iter()
            .map(|h| h.join().expect("join"))
            .collect()
    });
    let all: Vec<(f64, u32, u32)> = out.into_iter().flatten().collect();
    let mut per_switch: Vec<f64> = all.iter().map(|&(m, _, _)| m).collect();
    per_switch.sort_by(f64::total_cmp);
    let worst = all
        .iter()
        .copied()
        .fold((0.0f64, 0u32, 0u32), |acc, x| if x.0 > acc.0 { x } else { acc });
    let worst_adjacent = all
        .iter()
        .filter(|&&(_, a, b)| a.abs_diff(b) == 1)
        .fold(0.0f64, |acc, &(m, _, _)| acc.max(m));
    M121 {
        worst_cells: worst.0,
        worst_pair: (worst.1, worst.2),
        per_switch,
        worst_adjacent_cells: worst_adjacent,
    }
}

// ─── instrument 2: P-77's temporal-history rejection ────────────────────────
//
// Transcribed from `crates/isomesh/benches/experiment_p77.rs`, which is the
// instrument C2's cost half and C3 are denominated in. What is reused verbatim:
// `game_dig`'s `Ground` height field, its brush carve `max(f, -(|p-c|-r))`, its
// aim constants, its 12.5 edits/second held-button rate, Roberts 2018's 10-entry
// R2 jitter, the bilinear history fetch requiring four valid taps, the 3x3 YCoCg
// neighbourhood AABB clip (Karis 2014), the rejection definition -- *the clip
// moved the sample*, `s < 1.0` -- the `alpha = 0.1` blend, and the world-space
// albedo detail that `P-77`'s first run had to add after reading a 95.5%
// steady-state rate off a smooth shade.
//
// Three deviations from `P-77`, each with its reason:
//
// 1. **What is traced is the trilinear reconstruction of the field on a level's
//    own sample grid, not the analytic field.** That is the surface an extracted
//    mesh approximates, and it is the only way the two instruments in this file
//    look at the same pair of surfaces. A trilinear interpolant is a convex
//    combination of its eight cell corners, so no Lipschitz constant is declared:
//    the march is a half-cell walk with a sign test and a bisection, and a cell
//    whose corners share a sign cannot contain a crossing.
// 2. **The ray is clipped to `game_dig`'s sandbox box first** rather than
//    carrying the box test along inside the march. The sandbox is convex, so a
//    ray meets it in exactly one interval and the accepted crossing set is
//    identical; it is an optimisation, not a change of question.
// 3. **The geomorph arm gets exact motion vectors and the dither arm cannot.**
//    Under a morph a vertex genuinely moves and its previous position is known,
//    so `prev_world` is one Newton step onto the previous frame's weight surface
//    (`P-80`'s `project_once`). Under dither *nothing moves*: a renderer writes
//    zero motion, and a flipped pixel's history is a content mismatch rather than
//    a reprojection error. That asymmetry is not a thumb on the scale -- it **is**
//    Haydel, Yuksel & Seiler's structural argument, and it is the half of
//    `10.1145/3618359` that survives the port off their TRaX simulator.

/// `game_dig`'s fine cell: `CELL_SIZE = 0.125`.
const DIG_H_FINE: f32 = 0.125;
/// One LOD level coarser, which is what a transition switches to.
const DIG_H_COARSE: f32 = 0.25;
/// `game_dig::sandbox`'s lower corner.
const SANDBOX_LO: V3 = [-8.0, -5.4, -8.0];
/// Upper corner: `lo + (16, 8, 16)`.
const SANDBOX_HI: V3 = [8.0, 2.6, 8.0];
/// `game_dig::AIM_NEAR`.
const AIM_NEAR: f32 = 0.30;
/// `game_dig::AIM_FAR`.
const AIM_FAR: f32 = 25.0;
/// `game_dig`'s eye at `setup`: `Transform::from_xyz(0.0, 1.70, 6.0)`.
const DIG_EYE: V3 = [0.0, 1.70, 6.0];
/// `P-77`'s `HEADLINE_ARM` pose: looking down at the rock under your feet.
/// `Look::pitch` is user-driven and clamped to +-1.5, so this is `game_dig`.
const DIG_PITCH: f32 = -0.6;
/// `game_dig`'s default `World::radius`.
const DIG_BRUSH_RADIUS: f32 = 0.25;
/// `game_dig::EDIT_PERIOD`, i.e. 12.5 edits a second while the button is held.
const DIG_EDIT_PERIOD: f32 = 0.08;
/// 60 Hz.
const DIG_DT: f32 = 1.0 / 60.0;

/// Resolution of the software resolve. `P-77`'s `dig_at_feet_walk_half_res` arm.
const TAA_W: usize = 480;
const TAA_H: usize = 270;
/// Frames of history built before anything is measured, shared by every arm.
///
/// Every arm is identical during warm-up -- fade weight zero, no brushes -- so
/// the warm-up runs **once** and the history is cloned into each arm. That makes
/// the arms exactly frame-paired rather than merely similarly warmed.
const TAA_WARMUP: usize = 16;
/// Frames in the measurement window, and the length of the LOD cross-fade.
const TAA_WINDOW: usize = 8;
/// Karis's blend weight for the current frame.
const TAA_ALPHA: f32 = 0.1;
/// Vertical field of view, radians: Bevy's `PerspectiveProjection::default`.
const TAA_FOV_Y: f32 = core::f32::consts::FRAC_PI_4;
/// Central-difference half-step for the shading normal. `P-77`'s, and for
/// `P-77`'s reason: `game_dig::GRADIENT_EPS = 1e-4` leaves three significant
/// digits in an `f32` gradient and would make the rejection rate a measurement
/// of float differencing.
const TAA_NORMAL_EPS: f32 = 1e-3;
/// Hard ceiling on half-cell march steps. `AIM_FAR / (DIG_H_FINE/2)` is 400.
const TAA_MARCH_STEPS: u32 = 512;
/// Bisections once a bracket is found. Twenty-four halvings of a 0.0625 bracket
/// leaves 4e-9, which is below `game_dig::AIM_HIT = 0.01` by six orders.
const TAA_BISECT: u32 = 24;

/// Sun direction, `[0.35, 0.85, 0.40]` normalised. `P-77`'s.
const SUN_DIR: V3 = [0.349_128_2, 0.847_882_8, 0.399_003_6];
const ALBEDO_SUN: V3 = [1.00, 0.92, 0.78];
const ALBEDO_SKY: V3 = [0.22, 0.30, 0.45];
const ALBEDO_RIM: V3 = [0.10, 0.10, 0.12];
const FOG_RGB: V3 = [0.52, 0.58, 0.66];
const FOG_DENSITY: f32 = 0.06;
/// World units per albedo tile: `game_dig`'s `TriplanarExtension::settings.x`.
/// See `P-77`'s `DETAIL_TILE` for why a smooth shade is a fixture defect.
const DETAIL_TILE: f32 = 1.5;
const DETAIL_OCTAVES: u32 = 8;
const DETAIL_LUMA: f32 = 0.30;
const DETAIL_CHROMA: f32 = 0.10;

/// `game_dig::Ground`, verbatim: distance to a wavy height field, negative below.
#[derive(Clone, Copy)]
struct Ground;

impl Sdf for Ground {
    type Scalar = f32;

    fn sample(&self, p: V3) -> f32 {
        let height = 0.35 * (p[0] * 0.9).sin() * (p[2] * 0.7).cos() + 0.15 * (p[0] * 2.1).sin();
        p[1] - height
    }
}

/// `Ground` with spheres subtracted, which is the carve `game_dig`'s
/// `BrushStack` performs.
fn dug(brushes: &[[f32; 4]], p: V3) -> f32 {
    let mut f = Ground.sample(p);
    for b in brushes {
        let d = [p[0] - b[0], p[1] - b[1], p[2] - b[2]];
        let sphere = (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt() - b[3];
        f = f.max(-sphere);
    }
    f
}

/// One lattice cell's eight corner samples, held so a half-cell march and a
/// central-difference gradient do not re-evaluate them.
struct CellCache {
    cell: [i32; 3],
    v: [f32; 8],
    primed: bool,
}

impl CellCache {
    fn new() -> Self {
        Self {
            cell: [i32::MIN; 3],
            v: [0.0; 8],
            primed: false,
        }
    }

    /// Trilinear reconstruction of [`dug`] on a lattice of spacing `h` anchored
    /// at [`SANDBOX_LO`], which is the surface Marching Cubes approximates on
    /// that grid.
    fn sample(&mut self, brushes: &[[f32; 4]], h: f32, p: V3) -> f32 {
        let mut c = [0i32; 3];
        let mut f = [0f32; 3];
        for k in 0..3 {
            let x = (p[k] - SANDBOX_LO[k]) / h;
            let fl = x.floor();
            c[k] = fl as i32;
            f[k] = x - fl;
        }
        if !self.primed || c != self.cell {
            for i in 0..8 {
                let q = [
                    SANDBOX_LO[0] + (c[0] + (i & 1) as i32) as f32 * h,
                    SANDBOX_LO[1] + (c[1] + ((i >> 1) & 1) as i32) as f32 * h,
                    SANDBOX_LO[2] + (c[2] + ((i >> 2) & 1) as i32) as f32 * h,
                ];
                self.v[i] = dug(brushes, q);
            }
            self.cell = c;
            self.primed = true;
        }
        let mut acc = 0.0;
        for i in 0..8 {
            let wx = if i & 1 == 0 { 1.0 - f[0] } else { f[0] };
            let wy = if (i >> 1) & 1 == 0 { 1.0 - f[1] } else { f[1] };
            let wz = if (i >> 2) & 1 == 0 { 1.0 - f[2] } else { f[2] };
            acc += wx * wy * wz * self.v[i];
        }
        acc
    }
}

/// Which surface a pixel is looking at.
#[derive(Clone, Copy, PartialEq)]
enum Surf {
    /// The level the chunk is at before the switch.
    Fine,
    /// The level it is switching to.
    Coarse,
    /// Lengyel's morph, in field space: the zero set of
    /// `(1-w)*f_fine + w*f_coarse` moves continuously from one surface to the
    /// other, which is the field-space form of a vertex morph and needs no
    /// per-vertex correspondence to exist.
    Morph(f32),
}

/// The whole-transition policy an arm is running.
#[derive(Clone, Copy, PartialEq)]
enum Blend {
    /// No transition: the fine surface, every frame. The steady-state arm.
    Static,
    /// Alpha-to-coverage / ordered dither: a screen-space-fixed threshold per
    /// pixel decides which of the two surfaces that pixel shows, and the fade
    /// weight advances one step per frame. This is what a stochastic LOD blend
    /// does, and a pixel therefore flips surface on exactly one frame.
    Dither,
    /// The geomorph.
    Morph,
}

struct SurfEval<'a> {
    brushes: &'a [[f32; 4]],
    fine: CellCache,
    coarse: CellCache,
}

impl<'a> SurfEval<'a> {
    fn new(brushes: &'a [[f32; 4]]) -> Self {
        Self {
            brushes,
            fine: CellCache::new(),
            coarse: CellCache::new(),
        }
    }

    fn sample(&mut self, s: Surf, p: V3) -> f32 {
        let br = self.brushes;
        match s {
            Surf::Fine => self.fine.sample(br, DIG_H_FINE, p),
            Surf::Coarse => self.coarse.sample(br, DIG_H_COARSE, p),
            Surf::Morph(w) => {
                let a = self.fine.sample(br, DIG_H_FINE, p);
                let b = self.coarse.sample(br, DIG_H_COARSE, p);
                (1.0 - w) * a + w * b
            }
        }
    }

    /// Half a cell of the finest lattice the surface involves.
    fn step_for(s: Surf) -> f32 {
        match s {
            Surf::Coarse => DIG_H_COARSE * 0.5,
            Surf::Fine | Surf::Morph(_) => DIG_H_FINE * 0.5,
        }
    }

    fn gradient(&mut self, s: Surf, p: V3) -> V3 {
        let mut g = [0.0f32; 3];
        for k in 0..3 {
            let mut a = p;
            let mut b = p;
            a[k] += TAA_NORMAL_EPS;
            b[k] -= TAA_NORMAL_EPS;
            g[k] = self.sample(s, a) - self.sample(s, b);
        }
        g
    }

    /// One Newton step onto `s`'s zero set: `P-80`'s `project_once`. Used for the
    /// geomorph arm's exact motion vector and nothing else.
    fn project_once(&mut self, s: Surf, p: V3) -> V3 {
        let f = self.sample(s, p);
        let g = self.gradient(s, p);
        let gg = dot3(g, g);
        if gg <= 0.0 {
            return p;
        }
        // `gradient` is a central difference over `2*eps`, so it is `2*eps*grad f`.
        let scale = f * 2.0 * TAA_NORMAL_EPS / gg;
        sub3(p, mul3(g, scale))
    }
}

/// The ray's interval inside `game_dig`'s sandbox, intersected with the aim
/// range. `None` when the ray never enters.
fn clip_sandbox(o: V3, d: V3) -> Option<(f32, f32)> {
    let mut t0 = AIM_NEAR;
    let mut t1 = AIM_FAR;
    for k in 0..3 {
        if d[k].abs() < 1e-12 {
            if o[k] < SANDBOX_LO[k] || o[k] > SANDBOX_HI[k] {
                return None;
            }
            continue;
        }
        let inv = 1.0 / d[k];
        let mut a = (SANDBOX_LO[k] - o[k]) * inv;
        let mut b = (SANDBOX_HI[k] - o[k]) * inv;
        if a > b {
            core::mem::swap(&mut a, &mut b);
        }
        t0 = t0.max(a);
        t1 = t1.min(b);
        if t0 > t1 {
            return None;
        }
    }
    Some((t0, t1))
}

/// First crossing of `s`'s zero set along the ray, as a distance.
fn trace_surf(ev: &mut SurfEval<'_>, s: Surf, o: V3, d: V3) -> Option<f32> {
    let (t0, t1) = clip_sandbox(o, d)?;
    let step = SurfEval::step_for(s);
    let mut t = t0;
    if ev.sample(s, add3(o, mul3(d, t))) <= 0.0 {
        return Some(t);
    }
    let mut guard = 0u32;
    while t < t1 {
        guard += 1;
        if guard > TAA_MARCH_STEPS {
            break;
        }
        let tn = (t + step).min(t1);
        if ev.sample(s, add3(o, mul3(d, tn))) <= 0.0 {
            let (mut lo, mut hi) = (t, tn);
            for _ in 0..TAA_BISECT {
                let m = 0.5 * (lo + hi);
                if ev.sample(s, add3(o, mul3(d, m))) <= 0.0 {
                    hi = m;
                } else {
                    lo = m;
                }
            }
            return Some(hi);
        }
        if tn >= t1 {
            break;
        }
        t = tn;
    }
    None
}

/// The standard 8x8 ordered-dither (Bayer) matrix, 64 coverage levels.
///
/// Screen-space fixed, which is what alpha-to-coverage's sample mask is: the
/// pattern does not follow the surface, so a pixel's flip frame is decided by
/// the fade weight alone. With [`TAA_WINDOW`] = 8 the weight advances in eighths
/// and 8 of the 64 levels cross per frame -- 12.5% of pixels flip per frame.
const BAYER8: [[u8; 8]; 8] = [
    [0, 32, 8, 40, 2, 34, 10, 42],
    [48, 16, 56, 24, 50, 18, 58, 26],
    [12, 44, 4, 36, 14, 46, 6, 38],
    [60, 28, 52, 20, 62, 30, 54, 22],
    [3, 35, 11, 43, 1, 33, 9, 41],
    [51, 19, 59, 27, 49, 17, 57, 25],
    [15, 47, 7, 39, 13, 45, 5, 37],
    [63, 31, 55, 23, 61, 29, 53, 21],
];

fn bayer(x: usize, y: usize) -> f32 {
    (f32::from(BAYER8[y % 8][x % 8]) + 0.5) / 64.0
}

/// A pinhole camera at one instant, jitter baked in. `P-77`'s `Cam`.
#[derive(Clone, Copy)]
struct TaaCam {
    eye: V3,
    right: V3,
    up: V3,
    forward: V3,
    tan_half: f32,
    aspect: f32,
    width: f32,
    height: f32,
    jitter: [f32; 2],
}

impl TaaCam {
    /// Bevy's `Quat::from_euler(EulerRot::YXZ, yaw, pitch, 0.0)` written out.
    fn new(eye: V3, yaw: f32, pitch: f32, w: usize, h: usize, jitter: [f32; 2]) -> Self {
        let (sy, cy) = (yaw.sin(), yaw.cos());
        let (sp, cp) = (pitch.sin(), pitch.cos());
        Self {
            eye,
            right: [cy, 0.0, -sy],
            up: [sy * sp, cp, cy * sp],
            forward: [-sy * cp, sp, -cy * cp],
            tan_half: (TAA_FOV_Y * 0.5).tan(),
            aspect: w as f32 / h as f32,
            width: w as f32,
            height: h as f32,
            jitter,
        }
    }

    fn ray(&self, x: usize, y: usize) -> V3 {
        let px = (x as f32 + 0.5 + self.jitter[0]) / self.width;
        let py = (y as f32 + 0.5 + self.jitter[1]) / self.height;
        let sx = (2.0 * px - 1.0) * self.aspect * self.tan_half;
        let sy = (1.0 - 2.0 * py) * self.tan_half;
        unit3(add3(
            self.forward,
            add3(mul3(self.right, sx), mul3(self.up, sy)),
        ))
        .expect("the camera basis is orthonormal, so no pixel ray is degenerate")
    }

    /// Where a world point lands in this frame's pixel grid. The jitter is
    /// subtracted because the jittered sample at pixel `p` *is* the ray through
    /// `p + jitter`.
    fn project(&self, p: V3) -> Option<[f32; 2]> {
        let v = sub3(p, self.eye);
        let cz = dot3(v, self.forward);
        if cz <= AIM_NEAR {
            return None;
        }
        let ndc_x = dot3(v, self.right) / (cz * self.tan_half * self.aspect);
        let ndc_y = dot3(v, self.up) / (cz * self.tan_half);
        Some([
            (ndc_x + 1.0) * 0.5 * self.width - 0.5 - self.jitter[0],
            (1.0 - ndc_y) * 0.5 * self.height - 0.5 - self.jitter[1],
        ])
    }
}

/// Roberts 2018's R2 sequence, the 10-entry jitter pattern via `P-77`.
fn jitter_for(frame: usize) -> [f32; 2] {
    const A1: f32 = 0.754_877_7;
    const A2: f32 = 0.569_840_3;
    let k = (frame % 10 + 1) as f32;
    [(k * A1).fract() - 0.5, (k * A2).fract() - 0.5]
}

/// One 32-bit hash of a lattice cell in `[0, 1)`: splitmix64's finaliser
/// truncated. `P-77`'s.
fn hash_lattice(i: i32, j: i32, k: i32, seed: u32) -> f32 {
    let mut h = (i as u32).wrapping_mul(0x9E37_79B1)
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

fn value_noise(p: V3, seed: u32) -> f32 {
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

fn fbm(p: V3, octaves: u32, seed: u32) -> f32 {
    let mut freq = 1.0 / DETAIL_TILE;
    let mut amp = 1.0;
    let mut sum = 0.0;
    let mut total = 0.0;
    for o in 0..octaves {
        sum += amp * (value_noise(mul3(p, freq), seed.wrapping_add(o * 7919)) * 2.0 - 1.0);
        total += amp;
        freq *= 2.0;
        amp *= 0.5;
    }
    sum / total
}

/// Per-channel albedo at a world point, sampled on the **world** point so it is
/// a texture on the surface and reprojects exactly. `P-77`'s.
fn detail(p: V3) -> V3 {
    let luma = fbm(p, DETAIL_OCTAVES, 0x5EED_1234);
    let chroma = fbm(p, 4, 0x1234_5EED);
    [
        1.0 + DETAIL_LUMA * luma + DETAIL_CHROMA * chroma,
        1.0 + DETAIL_LUMA * luma,
        1.0 + DETAIL_LUMA * luma - DETAIL_CHROMA * chroma,
    ]
}

/// Three channels that are independent functions of the surface, plus distance
/// fog so that a depth change alone moves the signal. `P-77`'s.
///
/// `tex_point` is separated from `world` for one reason and it is the reason the
/// geomorph arm is bracketed. `P-77` samples the albedo detail on the **world**
/// hit point, which is exactly right for static geometry and exactly right for
/// `game_dig`'s own material -- its terrain is a **triplanar** projection, so the
/// albedo is a function of world position and nothing else. Under a geomorph the
/// surface slides through that field, so the albedo at a material point genuinely
/// changes; a UV-locked material would ride with the vertex instead. Neither is
/// "the" answer, so both are measured: `tex_point == world` is the triplanar case
/// and `tex_point ==` the vertex's un-morphed position is the UV-locked case.
///
/// **The control came out flat, and it falsified the hypothesis that motivated
/// it.** Triplanar **432,188** rejections against UV-locked **431,578** -- a
/// difference of **0.14%**, where the geomorph's excess over dither is **3.49x**.
/// So the geomorph's temporal cost is not the material sliding under it. The
/// mechanism it leaves standing is the one that cannot be locked away: a morph
/// moves the surface at **every** pixel on **every** frame, so every pixel's
/// normal, depth and fog change every frame, and a UV-locked albedo is sampled at
/// the *pixel's* material point, which moves too. A dither moves nothing and
/// flips 1/S of the pixels once. Under a per-pixel neighbourhood clamp,
/// continuous whole-frame motion destroys more history than a sparse
/// instantaneous flip -- which is the **inverse** of Haydel, Yuksel & Seiler's
/// ordering, because their budget is data movement at the transition and this one
/// is temporal reuse across it. Two different budgets, opposite verdicts.
fn shade(n: V3, ray: V3, t: f32, tex_point: V3) -> V3 {
    let sun = dot3(n, SUN_DIR).max(0.0);
    let sky = 0.5 * (n[1] + 1.0);
    let rim = 1.0 - dot3(n, mul3(ray, -1.0)).max(0.0);
    let fog = 1.0 - (-t * FOG_DENSITY).exp();
    let tex = detail(tex_point);
    let mut out = [0.0f32; 3];
    for c in 0..3 {
        let lit = ALBEDO_SUN[c] * sun + ALBEDO_SKY[c] * sky + ALBEDO_RIM[c] * rim;
        out[c] = lit * tex[c] * (1.0 - fog) + FOG_RGB[c] * fog;
    }
    out
}

/// One pixel of the G-buffer.
#[derive(Clone, Copy, Default)]
struct Px {
    hit: bool,
    world: V3,
    /// Where this surface point was last frame. Equal to `world` for every arm
    /// but the geomorph, where the surface genuinely moved.
    prev_world: V3,
    rgb: V3,
}

/// Which surface pixel `(x, y)` shows under `blend` at fade weight `w`.
fn surf_at(blend: Blend, w: f32, x: usize, y: usize) -> Surf {
    match blend {
        Blend::Static => Surf::Fine,
        Blend::Dither => {
            if bayer(x, y) < w {
                Surf::Coarse
            } else {
                Surf::Fine
            }
        }
        Blend::Morph => Surf::Morph(w),
    }
}

/// Ray-cast one frame. Rows are handed to `std::thread`s; the field is pure.
///
/// `uv_locked` only means anything for [`Blend::Morph`]: see [`shade`].
fn render_taa(
    cam: &TaaCam,
    brushes: &[[f32; 4]],
    blend: Blend,
    w: f32,
    w_prev: f32,
    uv_locked: bool,
    threads: usize,
) -> Vec<Px> {
    let mut buf = vec![Px::default(); TAA_W * TAA_H];
    let rows_per = TAA_H.div_ceil(threads);
    std::thread::scope(|s| {
        for (chunk_index, chunk) in buf.chunks_mut(rows_per * TAA_W).enumerate() {
            let y0 = chunk_index * rows_per;
            s.spawn(move || {
                let mut ev = SurfEval::new(brushes);
                for (local_y, row) in chunk.chunks_mut(TAA_W).enumerate() {
                    let y = y0 + local_y;
                    for (x, px) in row.iter_mut().enumerate() {
                        let surf = surf_at(blend, w, x, y);
                        let ray = cam.ray(x, y);
                        let Some(t) = trace_surf(&mut ev, surf, cam.eye, ray) else {
                            continue;
                        };
                        let world = add3(cam.eye, mul3(ray, t));
                        let n = unit3(ev.gradient(surf, world)).unwrap_or([0.0, 1.0, 0.0]);
                        // Only the morph moves geometry between frames, and only
                        // it can be given an exact motion vector. See deviation 3.
                        let prev_world = match surf {
                            Surf::Morph(_) => ev.project_once(Surf::Morph(w_prev), world),
                            Surf::Fine | Surf::Coarse => world,
                        };
                        // A UV-locked material's texture coordinate is the
                        // vertex's own, i.e. its position on the un-morphed
                        // surface. One Newton step onto `Surf::Fine` is that
                        // point, computed the same way the motion vector is.
                        let tex_point = match surf {
                            Surf::Morph(_) if uv_locked => {
                                ev.project_once(Surf::Fine, world)
                            }
                            _ => world,
                        };
                        *px = Px {
                            hit: true,
                            world,
                            prev_world,
                            rgb: shade(n, ray, t, tex_point),
                        };
                    }
                }
            });
        }
    });
    buf
}

/// Karis 2014's YCoCg, the space production TAA builds its AABB in.
fn to_ycocg(c: V3) -> V3 {
    [
        0.25 * c[0] + 0.5 * c[1] + 0.25 * c[2],
        0.5 * (c[0] - c[2]),
        -0.25 * c[0] + 0.5 * c[1] - 0.25 * c[2],
    ]
}

fn from_ycocg(c: V3) -> V3 {
    [c[0] + c[1] - c[2], c[0] + c[2], c[0] - c[1] - c[2]]
}

/// One TAA history buffer.
#[derive(Clone)]
struct Taa {
    hist: Vec<V3>,
    valid: Vec<bool>,
}

/// Bilinear fetch requiring all four taps valid. `P-77`'s.
fn fetch(taa: &Taa, p: [f32; 2]) -> Option<V3> {
    let x0 = p[0].floor();
    let y0 = p[1].floor();
    if x0 < 0.0 || y0 < 0.0 || x0 + 1.0 >= TAA_W as f32 || y0 + 1.0 >= TAA_H as f32 {
        return None;
    }
    let (ix, iy) = (x0 as usize, y0 as usize);
    let (fx, fy) = (p[0] - x0, p[1] - y0);
    let mut out = [0.0f32; 3];
    for (dy, wy) in [(0usize, 1.0 - fy), (1, fy)] {
        for (dx, wx) in [(0usize, 1.0 - fx), (1, fx)] {
            let i = (iy + dy) * TAA_W + ix + dx;
            if !taa.valid[i] {
                return None;
            }
            let k = wx * wy;
            for c in 0..3 {
                out[c] += taa.hist[i][c] * k;
            }
        }
    }
    Some(out)
}

/// The 3x3 neighbourhood of the current frame, clamped at the border.
fn neighbourhood(cur: &[Px], x: usize, y: usize) -> ([V3; 9], usize) {
    let mut out = [[0.0f32; 3]; 9];
    let mut n = 0;
    for dy in -1i32..=1 {
        for dx in -1i32..=1 {
            let sx = (x as i32 + dx).clamp(0, TAA_W as i32 - 1) as usize;
            let sy = (y as i32 + dy).clamp(0, TAA_H as i32 - 1) as usize;
            let p = &cur[sy * TAA_W + sx];
            if p.hit {
                out[n] = p.rgb;
                n += 1;
            }
        }
    }
    (out, n)
}

/// The axis-aligned neighbourhood clamp, as a ray clipped to the box shell.
/// Returns the parameter in `[0, 1]`; `1.0` means nothing was rejected.
fn clip_aabb(cols: &[V3], centre: V3, history: V3) -> f32 {
    let mut lo = [f32::INFINITY; 3];
    let mut hi = [f32::NEG_INFINITY; 3];
    for c in cols {
        for k in 0..3 {
            lo[k] = lo[k].min(c[k]);
            hi[k] = hi[k].max(c[k]);
        }
    }
    let d = sub3(history, centre);
    let mut s = 1.0f32;
    for k in 0..3 {
        if d[k].abs() < 1e-20 {
            continue;
        }
        let bound = if d[k] > 0.0 { hi[k] } else { lo[k] };
        let si = (bound - centre[k]) / d[k];
        if si < s {
            s = si;
        }
    }
    s.clamp(0.0, 1.0)
}

#[derive(Clone, Copy, Default)]
struct TaaFrame {
    population: u64,
    no_history: u64,
    rejected: u64,
    reproj_px: f64,
}

/// One frame's resolve: reproject, fetch, clamp, blend. The registered quantity
/// is `rejected` -- the clamp moved the sample, `s < 1.0`.
fn resolve_taa(taa: &mut Taa, cur: &[Px], prev_cam: Option<&TaaCam>) -> TaaFrame {
    let mut f = TaaFrame::default();
    let n = TAA_W * TAA_H;
    let mut hist_rgb: Vec<Option<V3>> = vec![None; n];
    if let Some(pc) = prev_cam {
        for y in 0..TAA_H {
            for x in 0..TAA_W {
                let i = y * TAA_W + x;
                if !cur[i].hit {
                    continue;
                }
                let Some(q) = pc.project(cur[i].prev_world) else {
                    f.no_history += 1;
                    continue;
                };
                match fetch(taa, q) {
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

    let mut clipped: Vec<V3> = vec![[0.0; 3]; n];
    for y in 0..TAA_H {
        for x in 0..TAA_W {
            let i = y * TAA_W + x;
            let Some(hrgb) = hist_rgb[i] else { continue };
            let (nb, k) = neighbourhood(cur, x, y);
            let mut cols = [[0.0f32; 3]; 9];
            for j in 0..k {
                cols[j] = to_ycocg(nb[j]);
            }
            let c = to_ycocg(cur[i].rgb);
            let hy = to_ycocg(hrgb);
            let s = clip_aabb(&cols[..k], c, hy);
            clipped[i] = from_ycocg(add3(c, mul3(sub3(hy, c), s)));
            if s < 1.0 {
                f.rejected += 1;
            }
        }
    }

    for i in 0..n {
        if !cur[i].hit {
            taa.valid[i] = false;
            continue;
        }
        taa.valid[i] = true;
        taa.hist[i] = match hist_rgb[i] {
            Some(_) => add3(mul3(clipped[i], 1.0 - TAA_ALPHA), mul3(cur[i].rgb, TAA_ALPHA)),
            None => cur[i].rgb,
        };
    }
    f
}

/// One cost arm.
struct CostArm {
    name: &'static str,
    blend: Blend,
    digging: bool,
    /// The material rides with the vertex rather than with the world position.
    /// Only meaningful under [`Blend::Morph`]; see [`shade`].
    uv_locked: bool,
}

struct CostOut {
    name: &'static str,
    rejected: u64,
    population: u64,
    brushes: usize,
}

/// The six arms C2's cost half and C3 need, over one shared warm-up.
///
/// The camera is `P-77`'s `HEADLINE_ARM` regime -- static, jitter on -- because
/// `P-77`'s second fixture defect is that the rejection rate is a function of
/// reprojection displacement and a walking camera saturates it at 86.6%. A ratio
/// between two saturated rates measures nothing, and C3 is a ratio.
///
/// The geomorph is run **twice**, triplanar and UV-locked, because the first run
/// of this harness read a geomorph rejection 3.5x the dither's and that number is
/// a property of the *material* rather than of the transition. Both are reported;
/// C2's cost half is scored against the arm that favours geomorph.
fn run_cost(threads: usize) -> Vec<CostOut> {
    let n = TAA_W * TAA_H;
    let mut warm = Taa {
        hist: vec![[0.0; 3]; n],
        valid: vec![false; n],
    };
    let no_brushes: Vec<[f32; 4]> = Vec::new();
    let mut prev_cam: Option<TaaCam> = None;
    for frame in 0..TAA_WARMUP {
        let cam = TaaCam::new(DIG_EYE, 0.0, DIG_PITCH, TAA_W, TAA_H, jitter_for(frame));
        let cur = render_taa(&cam, &no_brushes, Blend::Static, 0.0, 0.0, false, threads);
        resolve_taa(&mut warm, &cur, prev_cam.as_ref());
        prev_cam = Some(cam);
    }
    let warm_cam = prev_cam.expect("the warm-up always runs at least one frame");

    let arms = [
        CostArm {
            name: "static",
            blend: Blend::Static,
            digging: false,
            uv_locked: false,
        },
        CostArm {
            name: "dither_lod",
            blend: Blend::Dither,
            digging: false,
            uv_locked: false,
        },
        CostArm {
            name: "geomorph_lod",
            blend: Blend::Morph,
            digging: false,
            uv_locked: false,
        },
        CostArm {
            name: "geomorph_lod_uv_locked",
            blend: Blend::Morph,
            digging: false,
            uv_locked: true,
        },
        CostArm {
            name: "digging",
            blend: Blend::Static,
            digging: true,
            uv_locked: false,
        },
        CostArm {
            name: "dither_lod_and_digging",
            blend: Blend::Dither,
            digging: true,
            uv_locked: false,
        },
    ];

    let mut out = Vec::new();
    for arm in &arms {
        let started = std::time::Instant::now();
        let mut taa = warm.clone();
        let mut prev = Some(warm_cam);
        let mut brushes: Vec<[f32; 4]> = Vec::new();
        let mut since_edit = f32::INFINITY;
        let mut rejected = 0u64;
        let mut population = 0u64;
        for i in 0..TAA_WINDOW {
            let frame = TAA_WARMUP + i;
            let cam = TaaCam::new(DIG_EYE, 0.0, DIG_PITCH, TAA_W, TAA_H, jitter_for(frame));
            if arm.digging {
                // `game_dig`'s held button: a stroke every EDIT_PERIOD, aimed
                // down the camera's forward ray, placed at the hit.
                let place = if since_edit.is_infinite() {
                    since_edit = 0.0;
                    true
                } else {
                    since_edit += DIG_DT;
                    if since_edit >= DIG_EDIT_PERIOD {
                        since_edit -= DIG_EDIT_PERIOD;
                        true
                    } else {
                        false
                    }
                };
                if place {
                    let mut ev = SurfEval::new(&brushes);
                    if let Some(t) = trace_surf(&mut ev, Surf::Fine, cam.eye, cam.forward) {
                        let c = add3(cam.eye, mul3(cam.forward, t));
                        brushes.push([c[0], c[1], c[2], DIG_BRUSH_RADIUS]);
                    }
                }
            }
            // The fade advances one step per frame and completes on the last
            // frame of the window, so a transition is exactly TAA_WINDOW frames.
            let w = (i + 1) as f32 / TAA_WINDOW as f32;
            let w_prev = i as f32 / TAA_WINDOW as f32;
            let cur = render_taa(&cam, &brushes, arm.blend, w, w_prev, arm.uv_locked, threads);
            let f = resolve_taa(&mut taa, &cur, prev.as_ref());
            rejected += f.rejected;
            population += f.population;
            prev = Some(cam);
        }
        eprintln!(
            "P-91: cost arm {} in {:.1}s -- rejected {rejected} of {population}, {} brushes",
            arm.name,
            started.elapsed().as_secs_f64(),
            brushes.len()
        );
        out.push(CostOut {
            name: arm.name,
            rejected,
            population,
            brushes: brushes.len(),
        });
    }
    out
}

fn cost_of(out: &[CostOut], name: &str) -> u64 {
    out.iter()
        .find(|c| c.name == name)
        .map(|c| c.rejected)
        .expect("every cost arm runs")
}

fn cpu_mhz() -> String {
    std::fs::read_to_string("/sys/devices/system/cpu/cpu0/cpufreq/scaling_cur_freq")
        .ok()
        .and_then(|s| s.trim().parse::<f64>().ok())
        .map_or_else(|| "unknown".to_string(), |khz| format!("{:.0}", khz / 1000.0))
}

type Row = Vec<(&'static str, String)>;

/// One method's residual distribution on one (field, level pair), in world units.
struct Method {
    name: &'static str,
    p99_world: f64,
    max_world: f64,
    /// The registered `bytes_per_vertex`.
    bytes: u32,
    /// Coverage slots, 1 for the two non-dither methods.
    slots: u32,
}

/// The three methods on one row's geometry, plus the two extra dither slot
/// counts.
///
/// **Dither's residual is the full displacement divided by its coverage slot
/// count, and that is the whole model.** An alpha-to-coverage cross-fade with
/// `S` samples per pixel resolves a pixel to a coverage-weighted mixture of the
/// two surfaces, so transferring coverage one slot at a time moves the resolved
/// image by `delta/S` per step. `S = 1` is a hard switch scattered across pixels
/// and its row says so by reading exactly `none`'s numbers -- which is also the
/// assertion below that no accidental geometric benefit leaked into the dither
/// arm.
fn methods_for(p: &Pop) -> Vec<Method> {
    let none_p99 = percentile(&p.surface_sorted, 0.99);
    let none_max = p.surface_sorted.last().copied().unwrap_or(0.0);
    let geo_p99 = percentile(&p.geomorph_sorted, 0.99);
    let geo_max = p.geomorph_sorted.last().copied().unwrap_or(0.0);
    let mut out = vec![
        Method {
            name: "none",
            p99_world: none_p99,
            max_world: none_max,
            bytes: 0,
            slots: 1,
        },
        Method {
            name: "geomorph",
            p99_world: geo_p99,
            max_world: geo_max,
            bytes: BYTES_ONE_SCALAR,
            slots: 1,
        },
    ];
    for s in DITHER_SLOTS {
        let q = f64::from(s.min(FADE_FRAMES));
        out.push(Method {
            name: match s {
                1 => "dither_s1",
                4 => "dither_s4",
                _ => "dither",
            },
            p99_world: none_p99 / q,
            max_world: none_max / q,
            bytes: 0,
            slots: s,
        });
    }
    out
}

fn main() {
    if !std::env::args().any(|a| a == "--bench") {
        return;
    }
    let prereg = isomesh::experiment!("P-91");
    let threads = std::thread::available_parallelism().map_or(8, |n| n.get());
    let mhz = cpu_mhz();

    // ─── M-279's rule: agree with the old instrument where they overlap ──────
    let m121 = reproduce_m121();
    println!(
        "M-121 reproduction: worst {:.3} cells (pair {}->{}) over {} switches; \
         p10 {:.3} median {:.3} p90 {:.3}; adjacent-only worst {:.3}",
        m121.worst_cells,
        m121.worst_pair.0,
        m121.worst_pair.1,
        m121.per_switch.len(),
        percentile(&m121.per_switch, 0.1),
        percentile(&m121.per_switch, 0.5),
        percentile(&m121.per_switch, 0.9),
        m121.worst_adjacent_cells,
    );
    assert!(
        m121.worst_cells >= M121_WORST_POP_CELLS
            && m121.worst_cells <= M121_WORST_POP_CELLS * (1.0 + M121_TOLERANCE),
        "P-91: the pop instrument does not reproduce M-121. This sweep is a superset of the \
         demo's flight, so it must land at or just above the committed {M121_WORST_POP_CELLS} \
         cells; it read {:.4}. Below means a different measurement, far above a different \
         fixture.",
        m121.worst_cells
    );

    // ─── the pop, per field per level pair ──────────────────────────────────
    let terrain = FbmTerrain::<f32>::canonical();
    let gyroid = Gyroid::<f32>::canonical();
    let mut pops: Vec<Pop> = Vec::new();
    for level in PAIRS {
        pops.push(measure_pop(
            "fbm_terrain",
            &terrain,
            [0.0, -CROSS, -CROSS],
            level,
        ));
        pops.push(measure_pop(
            "gyroid",
            &gyroid,
            [-BLOCK_W * 0.5, -CROSS, -CROSS],
            level,
        ));
    }
    for p in &pops {
        report(p);
    }

    // ─── the cost, in P-77's rejection predicate ────────────────────────────
    let cost = run_cost(threads);
    let static_rej = cost_of(&cost, "static");
    let dither_rej = cost_of(&cost, "dither_lod");
    let morph_rej = cost_of(&cost, "geomorph_lod");
    let morph_uv_rej = cost_of(&cost, "geomorph_lod_uv_locked");
    let dig_rej = cost_of(&cost, "digging");
    let both_rej = cost_of(&cost, "dither_lod_and_digging");
    // P-77's own vacuity control, inherited because C3 divides by differences
    // taken against this arm.
    assert!(
        static_rej > 0,
        "P-91: the steady-state arm rejected zero history samples, so every excess below is \
         measured against a floor and C3's ratio is not a number"
    );
    let excess_lod = dither_rej as f64 - static_rej as f64;
    let excess_dig = dig_rej as f64 - static_rej as f64;
    let sum_of_parts = static_rej as f64 + excess_lod + excess_dig;
    let superadditivity = if excess_lod + excess_dig == 0.0 {
        f64::NAN
    } else {
        (both_rej as f64 - static_rej as f64) / (excess_lod + excess_dig)
    };
    // C3's arithmetic ceiling, from the fixture rather than from the clause, and
    // it is the `x51` recomputation this row owes.
    //
    // The combined arm's rejections are a subset of the union of the two causes:
    // a pixel is rejected because its history disagrees, and the dither's flip
    // set is the whole frame while the dig's changed set is the brush silhouette.
    // Superadditivity therefore cannot exceed what you get if **every**
    // dig-attributable rejection were also a fresh dither rejection, i.e.
    // `1 + excess_dig / (excess_lod + excess_dig)`. With `excess_dig` a small
    // fraction of `excess_lod` the clause is reachable but only barely, and that
    // ceiling is a property of a 0.25-radius brush at 3 units, not of the
    // hypothesis.
    let superadditivity_ceiling = if excess_lod + excess_dig == 0.0 {
        f64::NAN
    } else {
        1.0 + excess_dig / (excess_lod + excess_dig)
    };
    let c3 = superadditivity > 1.0;

    // ─── the clause verdicts, from the headline rows ────────────────────────
    //
    // The headline is level 0 -> 1 at `LEVEL_RANGE * 1 = 7.0`, which is the
    // switch `game_lod_flyover::level_for` performs first and therefore "the
    // current switch distance".
    let headline = |field: &str, method: &str| -> f64 {
        let p = pops
            .iter()
            .find(|p| p.field == field && p.level == 0)
            .expect("the level 0 row always exists");
        let m = methods_for(p);
        let m = m
            .iter()
            .find(|m| m.name == method)
            .expect("every method runs on every row");
        pixels_of(m.p99_world, p.switch_distance)
    };
    let c1_terrain = headline("fbm_terrain", "geomorph");
    let c1_gyroid = headline("gyroid", "geomorph");
    let c1 = c1_terrain < 1.0 && c1_gyroid >= 1.0;
    let c2_dither_gyroid = headline("gyroid", "dither");
    // C2's cost half is scored against the geomorph arm that favours geomorph --
    // the *lower* of the triplanar and UV-locked rejections -- so "dither costs
    // more" has to beat the best case rather than the worst.
    let morph_rej_best = morph_rej.min(morph_uv_rej);
    let c2 = c2_dither_gyroid < c1_gyroid && dither_rej > morph_rej_best;

    // ─── rows ───────────────────────────────────────────────────────────────
    let mut rows: Vec<Row> = Vec::new();
    for p in &pops {
        let d = p.switch_distance;
        let none_p99_px = pixels_of(percentile(&p.surface_sorted, 0.99), d);
        let none_max_px = pixels_of(p.surface_sorted.last().copied().unwrap_or(0.0), d);
        for m in methods_for(p) {
            let p99_px = pixels_of(m.p99_world, d);
            let max_px = pixels_of(m.max_world, d);
            let this_method_rej = match m.name {
                "geomorph" => morph_rej,
                "none" => static_rej,
                _ => dither_rej,
            };
            // C1 is a **conjunction over two fields with opposite directions**,
            // and `c1_holds` is therefore a single global verdict replicated on
            // every row. A reader who saw `c1_holds = false` on the `gyroid`
            // geomorph row would naturally read it as "gyroid's half failed",
            // which is the exact opposite of what happened. This column carries
            // *this row's own half*: on `fbm_terrain` the clause wants p99 under
            // one pixel, on `gyroid` it wants one pixel or more.
            let c1_half = if m.name == "geomorph" {
                let holds = if p.field == "fbm_terrain" {
                    p99_px < 1.0
                } else {
                    p99_px >= 1.0
                };
                holds.to_string()
            } else {
                "NA".to_string()
            };
            rows.push(vec![
                // ── registered ──────────────────────────────────────────────
                ("field", p.field.to_string()),
                ("switch_distance", format!("{d:.3}")),
                ("method", m.name.to_string()),
                ("p99_pixels_of_pop", format!("{p99_px:.4}")),
                ("max_pixels_of_pop", format!("{max_px:.4}")),
                ("pop_cells", format!("{:.6}", m.max_world / p.fine_h)),
                ("bytes_per_vertex", m.bytes.to_string()),
                ("history_rejected_static", static_rej.to_string()),
                ("history_rejected_lod", dither_rej.to_string()),
                (
                    "history_rejected_lod_and_digging",
                    both_rej.to_string(),
                ),
                ("sum_of_parts", format!("{sum_of_parts:.1}")),
                ("superadditivity", format!("{superadditivity:.6}")),
                (
                    "superadditivity_ceiling",
                    format!("{superadditivity_ceiling:.6}"),
                ),
                ("excess_rejected_lod", format!("{excess_lod:.1}")),
                ("excess_rejected_digging", format!("{excess_dig:.1}")),
                // The interaction in samples, so a "held" can be judged as a
                // mechanism or as a rounding: `both - static - excess_lod -
                // excess_dig`, i.e. what the combined arm rejected that pure
                // additivity does not account for.
                (
                    "superadditivity_excess_samples",
                    format!("{:.0}", both_rej as f64 - sum_of_parts),
                ),
                ("c1_holds", c1.to_string()),
                ("c2_holds", c2.to_string()),
                ("c3_holds", c3.to_string()),
                ("c1_half_this_row_holds", c1_half),
                (
                    "c2_pop_half_holds",
                    (c2_dither_gyroid < c1_gyroid).to_string(),
                ),
                (
                    "c2_cost_half_holds",
                    (dither_rej > morph_rej_best).to_string(),
                ),
                // ── the pixel conversion, on the row rather than in a comment ─
                ("view_height_px", format!("{VIEW_HEIGHT_PX:.0}")),
                ("fov_y_deg", format!("{:.1}", FOV_Y.to_degrees())),
                (
                    "pixels_per_world_unit_at_distance",
                    format!("{:.4}", pixels_of(1.0, d)),
                ),
                ("pop_world_p99", format!("{:.6}", m.p99_world)),
                ("pop_world_max", format!("{:.6}", m.max_world)),
                ("pop_cells_p99", format!("{:.6}", m.p99_world / p.fine_h)),
                // Exact, because `pixels_of` is exactly proportional to 1/d:
                // the distance at which this row's p99 falls under one pixel.
                (
                    "distance_for_p99_under_one_px",
                    format!("{:.2}", d * p99_px),
                ),
                // ── the ladder ───────────────────────────────────────────────
                ("lod_from", p.level.to_string()),
                ("lod_to", (p.level + 1).to_string()),
                ("fine_h", format!("{:.4}", p.fine_h)),
                ("coarse_h", format!("{:.4}", p.coarse_h)),
                ("fine_vertices", p.fine_vertices.to_string()),
                ("coarse_vertices", p.coarse_vertices.to_string()),
                ("fine_triangles", p.fine_triangles.to_string()),
                ("coarse_triangles", p.coarse_triangles.to_string()),
                ("degenerate_coarse_triangles", p.degenerate_coarse.to_string()),
                ("fine_bytes", p.fine_bytes.to_string()),
                ("coarse_bytes", p.coarse_bytes.to_string()),
                // ── the two displacement measures ───────────────────────────
                (
                    "pop_cells_vertex_metric_max",
                    format!("{:.6}", p.pop_cells_vertex_metric),
                ),
                ("pop_cells_surface_max", format!("{:.6}", p.pop_cells_surface)),
                (
                    "p99_pixels_vertex_metric",
                    format!("{:.4}", pixels_of(percentile(&p.vertex_sorted, 0.99), d)),
                ),
                (
                    "p99_pixels_across_view",
                    format!("{:.4}", pixels_of(percentile(&p.across_sorted, 0.99), d)),
                ),
                // ── the dither model ────────────────────────────────────────
                ("dither_slots", m.slots.to_string()),
                ("fade_frames", FADE_FRAMES.to_string()),
                (
                    "dither_flip_fraction_per_frame",
                    format!("{:.6}", 1.0 / f64::from(FADE_FRAMES)),
                ),
                // ── item 3.6's own harness, three restriction rules ─────────
                (
                    "geomorph_fail_fraction_cell",
                    format!("{:.6}", p.geomorph.cell.fail_fraction()),
                ),
                (
                    "geomorph_fail_fraction_cell27",
                    format!("{:.6}", p.geomorph.cell27.fail_fraction()),
                ),
                (
                    "geomorph_fail_fraction_reach",
                    format!("{:.6}", p.geomorph.reach.fail_fraction()),
                ),
                ("geomorph_success_cell", p.geomorph.cell.success.to_string()),
                ("geomorph_no_hit_cell", p.geomorph.cell.no_hit.to_string()),
                ("geomorph_flipped_cell", p.geomorph.cell.flipped.to_string()),
                (
                    "geomorph_no_normal_cell",
                    p.geomorph.cell.no_normal.to_string(),
                ),
                (
                    "geomorph_success_reach",
                    p.geomorph.reach.success.to_string(),
                ),
                ("geomorph_no_hit_reach", p.geomorph.reach.no_hit.to_string()),
                (
                    "geomorph_flipped_reach",
                    p.geomorph.reach.flipped.to_string(),
                ),
                (
                    "morph_t_max_coarse_cells_cell",
                    format!("{:.6}", p.geomorph.cell.morph_t_max_coarse_cells),
                ),
                (
                    "morph_t_max_coarse_cells_reach",
                    format!("{:.6}", p.geomorph.reach.morph_t_max_coarse_cells),
                ),
                ("bytes_per_vertex_doc_claim", BYTES_DOC_CLAIM.to_string()),
                (
                    "bytes_second_position",
                    BYTES_SECOND_POSITION.to_string(),
                ),
                // ── THE REGISTERED VACUITY CONTROL, on every row ────────────
                ("vacuity_none_p99_pixels", format!("{none_p99_px:.4}")),
                ("vacuity_none_max_pixels", format!("{none_max_px:.4}")),
                // ── M-279's agreement check ─────────────────────────────────
                (
                    "m121_committed_worst_pop_cells",
                    format!("{M121_WORST_POP_CELLS:.3}"),
                ),
                (
                    "m121_reproduced_worst_pop_cells",
                    format!("{:.4}", m121.worst_cells),
                ),
                (
                    "m121_reproduced_adjacent_only",
                    format!("{:.4}", m121.worst_adjacent_cells),
                ),
                (
                    "m121_reproduced_median",
                    format!("{:.4}", percentile(&m121.per_switch, 0.5)),
                ),
                // ── the cost instrument ─────────────────────────────────────
                (
                    "history_rejected_lod_geomorph",
                    morph_rej.to_string(),
                ),
                (
                    "history_rejected_lod_geomorph_uv_locked",
                    morph_uv_rej.to_string(),
                ),
                ("history_rejected_digging", dig_rej.to_string()),
                (
                    "history_rejected_this_method",
                    this_method_rej.to_string(),
                ),
                (
                    "taa_population_static",
                    cost
                        .iter()
                        .find(|c| c.name == "static")
                        .map_or(0, |c| c.population)
                        .to_string(),
                ),
                (
                    "taa_brushes_placed",
                    cost
                        .iter()
                        .find(|c| c.name == "digging")
                        .map_or(0, |c| c.brushes)
                        .to_string(),
                ),
                ("taa_width", TAA_W.to_string()),
                ("taa_height", TAA_H.to_string()),
                ("taa_window_frames", TAA_WINDOW.to_string()),
                ("taa_warmup_frames", TAA_WARMUP.to_string()),
                ("dig_h_fine", format!("{DIG_H_FINE:.4}")),
                ("dig_h_coarse", format!("{DIG_H_COARSE:.4}")),
                // ── machine (M-280) ─────────────────────────────────────────
                ("threads", threads.to_string()),
                ("cpu_mhz", mhz.clone()),
            ]);
        }
    }

    // ─── THE REGISTERED VACUITY CONTROL, asserted off the emitted values ────
    //
    // "the no-transition arm must show a measurable pop, reported in pixels, or
    // both methods are being compared against an invisible baseline." Read back
    // off the rows about to be written rather than off an internal accumulator,
    // so what is asserted is what the CSV says (`P-80`'s pattern).
    let column = |row: &Row, key: &str| -> String {
        row.iter()
            .find(|(k, _)| *k == key)
            .map(|(_, v)| v.clone())
            .expect("every row carries every column")
    };
    let mut none_rows = 0usize;
    for row in &rows {
        if column(row, "method") != "none" {
            continue;
        }
        none_rows += 1;
        for key in ["p99_pixels_of_pop", "max_pixels_of_pop", "pop_cells"] {
            let v: f64 = column(row, key).parse().expect("a number");
            assert!(
                v > 0.0,
                "P-91: VACUITY CONTROL FAILED -- the no-transition arm on {} L{}->{} reports \
                 {key} = {v}, so geomorph and dither are being compared against an invisible \
                 baseline",
                column(row, "field"),
                column(row, "lod_from"),
                column(row, "lod_to"),
            );
        }
    }
    assert_eq!(
        none_rows,
        pops.len(),
        "P-91: every (field, level pair) must carry a no-transition row or the vacuity control \
         is not applied to that row's geometry"
    );
    // The other half of the same control: a dither with one coverage slot **is**
    // a hard switch, so its numbers must equal `none`'s exactly. If they do not,
    // a geometric benefit leaked into the dither model and every dither number
    // above is an artefact of this harness rather than of alpha-to-coverage.
    for p in &pops {
        let m = methods_for(p);
        let none = m.iter().find(|m| m.name == "none").expect("none runs");
        let s1 = m
            .iter()
            .find(|m| m.name == "dither_s1")
            .expect("dither_s1 runs");
        assert!(
            (none.p99_world - s1.p99_world).abs() <= f64::EPSILON * none.p99_world.max(1.0)
                && (none.max_world - s1.max_world).abs() <= f64::EPSILON * none.max_world.max(1.0),
            "P-91 {} L{}: a one-slot dither is a hard switch scattered across pixels and must \
             read exactly `none`; {:.9} vs {:.9} says the dither model gained a geometric \
             benefit it has no mechanism for",
            p.field,
            p.level,
            none.p99_world,
            s1.p99_world
        );
    }

    println!();
    println!("P-91 headline, at the repo's own switch distance for level 0 -> 1 (7.0 units),");
    println!("  {VIEW_HEIGHT_PX:.0} px tall, {:.0} deg vertical FOV:", FOV_Y.to_degrees());
    println!(
        "  VACUITY CONTROL   no-transition p99  fbm_terrain {:.2} px   gyroid {:.2} px",
        headline("fbm_terrain", "none"),
        headline("gyroid", "none"),
    );
    println!(
        "  C1  geomorph p99  fbm_terrain {c1_terrain:.2} px (needs < 1.0)   \
         gyroid {c1_gyroid:.2} px (needs >= 1.0)  ->  {}",
        if c1 { "HELD" } else { "FALSIFIED" }
    );
    println!(
        "  C2  gyroid p99    dither {c2_dither_gyroid:.2} px vs geomorph {c1_gyroid:.2} px  \
         -> pop half {}",
        if c2_dither_gyroid < c1_gyroid {
            "HELD"
        } else {
            "FALSIFIED"
        }
    );
    println!(
        "      rejection     dither {dither_rej} vs geomorph {morph_rej} (triplanar) / \
         {morph_uv_rej} (UV-locked)  -> cost half {}  ->  C2 {}",
        if dither_rej > morph_rej_best {
            "HELD"
        } else {
            "FALSIFIED"
        },
        if c2 { "HELD" } else { "FALSIFIED" }
    );
    println!(
        "  C3  rejection     static {static_rej}  dither {dither_rej}  digging {dig_rej}  \
         both {both_rej}  sum_of_parts {sum_of_parts:.0}",
    );
    println!(
        "      superadditivity {superadditivity:.4} against a ceiling of \
         {superadditivity_ceiling:.4}  ->  {}",
        if c3 { "HELD" } else { "FALSIFIED" }
    );
    println!(
        "      interaction     {:.0} samples of {:.0} excess -- a 'held' by sign, not by \
         mechanism: this IS additivity to engineering precision, which is C3's own \
         registered falsifier",
        both_rej as f64 - sum_of_parts,
        excess_lod + excess_dig,
    );
    println!(
        "  item 3.8's rule   a hard switch needs {:.0} units for p99 < 1 px on fbm_terrain, \
         {:.0} on gyroid",
        7.0 * headline("fbm_terrain", "none"),
        7.0 * headline("gyroid", "none"),
    );
    println!();

    common::experiment::run(prereg, |run| {
        for row in &rows {
            run.record(row);
        }
    });
}

fn report(p: &Pop) {
    let d = p.switch_distance;
    println!(
        "{:<12} L{}->{} d={:>4.1} h={:.3}/{:.3} verts {:>6}/{:>6} tris {:>6}/{:>6} \
         pop_cells vtx={:.3} surf={:.3}",
        p.field,
        p.level,
        p.level + 1,
        d,
        p.fine_h,
        p.coarse_h,
        p.fine_vertices,
        p.coarse_vertices,
        p.fine_triangles,
        p.coarse_triangles,
        p.pop_cells_vertex_metric,
        p.pop_cells_surface,
    );
    println!(
        "  none p99 {:>8.2}px max {:>8.2}px | geomorph p99 {:>8.2}px max {:>8.2}px | \
         fail {:.5} (no_hit {} flipped {} no_normal {}) t_max {:.3} coarse cells",
        pixels_of(percentile(&p.surface_sorted, 0.99), d),
        pixels_of(p.surface_sorted.last().copied().unwrap_or(0.0), d),
        pixels_of(percentile(&p.geomorph_sorted, 0.99), d),
        pixels_of(p.geomorph_sorted.last().copied().unwrap_or(0.0), d),
        p.geomorph.cell.fail_fraction(),
        p.geomorph.cell.no_hit,
        p.geomorph.cell.flipped,
        p.geomorph.cell.no_normal,
        p.geomorph.cell.morph_t_max_coarse_cells,
    );
    println!(
        "  across-view p99 {:>8.2}px | vertex-metric p99 {:>8.2}px | degenerate coarse tris {} | \
         bytes fine {} coarse {} | origin {:?}",
        pixels_of(percentile(&p.across_sorted, 0.99), d),
        pixels_of(percentile(&p.vertex_sorted, 0.99), d),
        p.degenerate_coarse,
        p.fine_bytes,
        p.coarse_bytes,
        p.origin,
    );
}
