//! **P-47 — the analytic gradient stops at the eight reference fields.**
//!
//! Ticket: R-043. Pre-registered before this harness existed.
//!
//! ```bash
//! cargo bench --bench experiment_p47
//! ```
//!
//! Writes `docs/experiments/p-47.csv`.
//!
//! # The hole
//!
//! All eight reference fields override [`Sdf::gradient`] analytically, and
//! `fields/mod.rs` says the central-difference default is never used by one.
//! That is true of the eight and stops being true of anything built out of them:
//! [`BrushStack`] does not override `gradient`, and neither does
//! [`Capsule`]. So on a sculpted field every vertex normal
//! `MarchingCubes` emits is the six-sample central difference at `O(h²)`, and
//! each of those six samples evaluates the **whole tape**.
//!
//! Two things are measured, and the accuracy one is the half the crate cannot
//! currently state at all:
//!
//! - **How wrong** the central-difference normal is against the exact one.
//! - **How much slower**, since one traversal carrying derivatives replaces six
//!   traversals that do not.
//!
//! # The instrument: forward-mode dual numbers
//!
//! [`Dual`] is `(value, [dx, dy, dz])` and it lives in this bench — nothing in
//! `crates/isomesh/src` is touched. The shape formulas are re-expressed over
//! `Dual` **operation for operation** against the crate's own `sample`, so the
//! derivative falls out of the chain rule rather than being hand-derived, and the
//! *value* comes out bit-identical. That last part is checked, not assumed:
//! `dual_value_bit_exact` compares `Dual::v` against `Sdf::sample` at every
//! vertex, so a transcription slip in the formulas shows up as a column rather
//! than as a plausible-looking error figure.
//!
//! Three details are decisions rather than mechanics:
//!
//! - **`min`/`max` propagate the selected branch's derivative**, which is the
//!   same selection argument P-39's pruning rested on. A tie keeps the **left**
//!   operand — the lower index in the fold — so a point sitting exactly on a CSG
//!   seam gets one answer and not whichever one the optimiser felt like. The
//!   derivative there is a one-sided limit; that is inherent to a
//!   non-differentiable field, and the choice is at least a pure function of the
//!   inputs.
//! - **`clamp` propagates zero in the clamped region**, which is the true
//!   derivative there, and keeps the interior branch on a tie by the same rule.
//! - **`sqrt` at zero is an error, not a number.** `d/dx √x` is unbounded at the
//!   origin, which is reachable: a point exactly on a capsule's segment. The
//!   crate has [`Error::DegenerateNormal`] for a vertex with no normal to derive
//!   and says it is *"reported rather than substituted"*. So it is reported —
//!   counted in `degenerate_points`, and those vertices are excluded from both
//!   arms rather than being given a fabricated normal in one of them.
//!
//! # Which normal each arm is
//!
//! The central-difference arm is not a lookalike written here: it is
//! `field.gradient(p)` — the default trait implementation, the one
//! `marching_cubes::unit_gradient` calls — normalised by the same `len.recip()`.
//! `central_matches_mesh_normals` checks the result bit-for-bit against the
//! normals `MarchingCubes` actually wrote into the mesh, so the arm under test is
//! demonstrably the crate's own output.
//!
//! # Rows
//!
//! `brush_stack_64` is the registered fixture: P-39's tape, 64 `Add`/`Subtract`
//! brushes over a sphere of radius 6, meshed across a 4×4×4 chunk world. The
//! other three rows are extras and decide nothing:
//!
//! - `brush_stack_64_smooth` — the same tape with every `Add` a `SmoothAdd`,
//!   which is the only thing that exercises the `smooth_min` dual path at all.
//! - `capsule` — one bare capsule, the other field in this crate with no
//!   analytic gradient. Separates "a tape is bad" from "a capsule is bad".
//! - `sphere` — clause three as a row rather than only as a scalar. The sphere
//!   *does* override `gradient`, so its error is the harness's own noise floor.

mod common;

use std::hint::black_box;
use std::ops::{Add, Mul, Neg, Sub};
use std::time::Instant;

use isomesh::brush::{Brush, BrushOp, BrushStack, Capsule};
use isomesh::chunk::{ChunkId, ChunkLayout};
use isomesh::fields::Sphere;
use isomesh::marching_cubes::MarchingCubes;
use isomesh::{Error, MeshBuffer, Real, RuntimeShape3, Sdf};

/// Chunks along each axis of the sculpted world.
const CHUNKS_PER_AXIS: i32 = 4;

/// Cells per chunk axis, so 33 samples per axis.
const CELLS_PER_CHUNK: u32 = 32;

/// World units per cell. The chunk is 4 units across and the world 16.
const CELL_SIZE: f64 = 0.125;

/// The world's minimum corner, so the world is `[-8, 8]³`.
const WORLD_ORIGIN: f64 = -8.0;

/// Brushes in the tape.
const BRUSHES: usize = 64;

/// Radius of the solid the brushes carve.
const BASE_RADIUS: f64 = 6.0;

/// Join width for the smooth extra, in world units — two cells.
const SMOOTH_K: f64 = 0.25;

/// Samples per axis for the two single-grid extras.
const SINGLE_GRID_SAMPLES: u32 = 65;

/// Timed repetitions per arm. The median is reported; a mean would be dragged by
/// whichever run collided with a scheduler tick or a sibling's benchmark.
const REPS: usize = 5;

/// Clause three's threshold: the dual normal must reproduce `p/|p|` this well.
const SPHERE_CONTROL_TOLERANCE: f64 = 1e-12;

/// Angular-error thresholds the vertex counts are taken at, in degrees.
///
/// A ladder rather than one number, because the error turned out to be bimodal:
/// rounding noise almost everywhere and whole degrees on the handful of
/// vertices that landed on a CSG seam. A single count at 1° would hide which of
/// those two the row is describing.
const THRESHOLDS_DEG: [f64; 5] = [1e-3, 1e-2, 1e-1, 1.0, 5.0];

/// Errors at or under this are the smooth-region bulk rather than a seam.
const BULK_DEG: f64 = 1e-3;

/// Tape lengths swept as an extra, to find where one dual traversal starts
/// beating six scalar ones. Prefixes of the same 64-brush tape.
const TAPE_LENGTHS: [usize; 7] = [1, 2, 4, 8, 16, 32, 64];

/// Row names for [`TAPE_LENGTHS`]. `&'static str` because the `fixture` column
/// is one, and a formatted name would have to be leaked to become one.
const TAPE_NAMES: [&str; 7] = [
    "brush_stack_1",
    "brush_stack_2",
    "brush_stack_4",
    "brush_stack_8",
    "brush_stack_16",
    "brush_stack_32",
    "brush_stack_64",
];

// ── forward-mode dual numbers ───────────────────────────────────────────────

/// A value and its gradient with respect to the sample point.
///
/// Deliberately not generic over [`isomesh::Real`]: this experiment is `f64`
/// only, and a type parameter here would buy nothing but a wall of bounds.
#[derive(Clone, Copy, Debug)]
struct Dual {
    v: f64,
    d: [f64; 3],
}

impl Dual {
    /// A quantity that does not depend on the sample point.
    const fn constant(v: f64) -> Self {
        Self { v, d: [0.0; 3] }
    }

    /// The sample point, seeded so every derivative is with respect to it.
    fn seed(p: [f64; 3]) -> [Self; 3] {
        [
            Self {
                v: p[0],
                d: [1.0, 0.0, 0.0],
            },
            Self {
                v: p[1],
                d: [0.0, 1.0, 0.0],
            },
            Self {
                v: p[2],
                d: [0.0, 0.0, 1.0],
            },
        ]
    }

    /// # Errors
    ///
    /// [`Error::DegenerateNormal`] at zero, where the derivative is unbounded.
    /// A capsule's segment is reachable, so this is a real case and not a
    /// formality — and substituting a number here is exactly what the variant's
    /// documentation forbids.
    fn sqrt(self) -> Result<Self, Error> {
        let root = self.v.sqrt();
        // `is_nan` is spelled out rather than folded into a negated comparison:
        // a non-finite gradient is the other half of what the variant covers.
        if root.is_nan() || root <= 0.0 {
            return Err(Error::DegenerateNormal { vertex: 0 });
        }
        let g = 0.5 / root;
        Ok(Self {
            v: root,
            d: [self.d[0] * g, self.d[1] * g, self.d[2] * g],
        })
    }

    /// Divide by a quantity that does not depend on the sample point.
    ///
    /// Kept distinct from a general quotient rule so the value is the crate's
    /// `x / k` and not `x * (1 / k)`, which is a different float.
    fn div_const(self, k: f64) -> Self {
        Self {
            v: self.v / k,
            d: [self.d[0] / k, self.d[1] / k, self.d[2] / k],
        }
    }

    /// `min`, carrying the selected branch's derivative.
    ///
    /// A tie keeps `self`, which is the earlier operand and so the lower index
    /// in the fold. The value matches `f64::min` on every input this experiment
    /// produces, and where the two could differ — `+0.0` against `-0.0` — the
    /// values are equal and only this one is deterministic.
    fn min(self, other: Self) -> Self {
        if other.v < self.v { other } else { self }
    }

    /// `max`, same rule.
    fn max(self, other: Self) -> Self {
        if other.v > self.v { other } else { self }
    }

    /// `clamp` to `[0, 1]`. Zero derivative where clamped, which is the true
    /// derivative there; the interior branch wins a tie.
    fn clamp01(self) -> Self {
        if self.v < 0.0 {
            Self::constant(0.0)
        } else if self.v > 1.0 {
            Self::constant(1.0)
        } else {
            self
        }
    }
}

impl Add for Dual {
    type Output = Self;
    fn add(self, o: Self) -> Self {
        Self {
            v: self.v + o.v,
            d: [self.d[0] + o.d[0], self.d[1] + o.d[1], self.d[2] + o.d[2]],
        }
    }
}

impl Sub for Dual {
    type Output = Self;
    fn sub(self, o: Self) -> Self {
        Self {
            v: self.v - o.v,
            d: [self.d[0] - o.d[0], self.d[1] - o.d[1], self.d[2] - o.d[2]],
        }
    }
}

impl Mul for Dual {
    type Output = Self;
    fn mul(self, o: Self) -> Self {
        Self {
            v: self.v * o.v,
            d: [
                self.d[0] * o.v + self.v * o.d[0],
                self.d[1] * o.v + self.v * o.d[1],
                self.d[2] * o.v + self.v * o.d[2],
            ],
        }
    }
}

impl Neg for Dual {
    type Output = Self;
    fn neg(self) -> Self {
        Self {
            v: -self.v,
            d: [-self.d[0], -self.d[1], -self.d[2]],
        }
    }
}

/// `vec3::dot` over duals, in the crate's left-to-right order so the value is
/// the same float.
fn dot(a: [Dual; 3], b: [Dual; 3]) -> Dual {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

/// `vec3::length` over duals.
///
/// # Errors
///
/// [`Error::DegenerateNormal`] on a zero-length vector.
fn length(a: [Dual; 3]) -> Result<Dual, Error> {
    dot(a, a).sqrt()
}

/// A field that can be differentiated exactly by carrying a dual through its own
/// evaluation, rather than by sampling it six times.
trait DualSdf: Sdf<Scalar = f64> {
    /// The value and its exact gradient at `p`, in one traversal.
    ///
    /// # Errors
    ///
    /// [`Error::DegenerateNormal`] where the derivative does not exist.
    fn dual(&self, p: [Dual; 3]) -> Result<Dual, Error>;
}

impl DualSdf for Sphere<f64> {
    /// `length(p − c) − r`, operation for operation against `Sphere::sample`.
    fn dual(&self, p: [Dual; 3]) -> Result<Dual, Error> {
        let d = [
            p[0] - Dual::constant(self.center[0]),
            p[1] - Dual::constant(self.center[1]),
            p[2] - Dual::constant(self.center[2]),
        ];
        Ok(length(d)? - Dual::constant(self.radius))
    }
}

impl DualSdf for Capsule<f64> {
    /// Point-to-segment distance, operation for operation against
    /// `Capsule::sample`.
    ///
    /// `ab` and `denom` are constants: the segment does not move with the sample
    /// point, so they carry a zero derivative rather than being special-cased
    /// out of the arithmetic.
    fn dual(&self, p: [Dual; 3]) -> Result<Dual, Error> {
        let ab = [
            Dual::constant(self.b[0] - self.a[0]),
            Dual::constant(self.b[1] - self.a[1]),
            Dual::constant(self.b[2] - self.a[2]),
        ];
        let ap = [
            p[0] - Dual::constant(self.a[0]),
            p[1] - Dual::constant(self.a[1]),
            p[2] - Dual::constant(self.a[2]),
        ];
        let denom = dot(ab, ab);
        // A zero-length capsule is a sphere, which is the crate's own comment
        // and the crate's own branch.
        let t = if denom.v > 0.0 {
            dot(ap, ab).div_const(denom.v).clamp01()
        } else {
            Dual::constant(0.0)
        };
        let q = [ap[0] - ab[0] * t, ap[1] - ab[1] * t, ap[2] - ab[2] * t];
        Ok(length(q)? - Dual::constant(self.radius))
    }
}

/// `brush::apply` over duals.
fn apply_dual(op: BrushOp, field: Dual, shape: Dual) -> Dual {
    match op {
        BrushOp::Add => field.min(shape),
        BrushOp::Subtract => field.max(-shape),
        BrushOp::SmoothAdd { k } => smooth_min_dual(field, shape, k),
    }
}

/// `brush::smooth_min` over duals.
///
/// The analytic derivative of the polynomial, obtained by the chain rule over
/// the same expression the crate evaluates — including the `- k·h·(1 − h)`
/// term, whose derivative is the part a hand-written gradient gets wrong.
fn smooth_min_dual(a: Dual, b: Dual, k: f64) -> Dual {
    if k <= 0.0 {
        return a.min(b);
    }
    let half = Dual::constant(0.5);
    let one = Dual::constant(1.0);
    let kd = Dual::constant(k);
    let h = (half + (half * (b - a)).div_const(k)).clamp01();
    (b + (a - b) * h) - kd * h * (one - h)
}

// ── the fixture, copied from experiment_p39 ─────────────────────────────────

/// A brush shape. One enum so the stack is a single `&[Brush<Shape>]`.
#[derive(Clone, Copy, Debug)]
enum Shape {
    Sphere(Sphere<f64>),
    Capsule(Capsule<f64>),
}

impl Sdf for Shape {
    type Scalar = f64;

    fn sample(&self, p: [f64; 3]) -> f64 {
        match self {
            Self::Sphere(s) => s.sample(p),
            Self::Capsule(c) => c.sample(p),
        }
    }
}

impl DualSdf for Shape {
    fn dual(&self, p: [Dual; 3]) -> Result<Dual, Error> {
        match self {
            Self::Sphere(s) => s.dual(p),
            Self::Capsule(c) => c.dual(p),
        }
    }
}

impl<F, S> DualSdf for BrushStack<'_, F, S>
where
    F: DualSdf,
    S: DualSdf,
{
    /// One traversal of the tape, carrying the gradient. This is the whole
    /// mechanism: the same fold `BrushStack::sample` performs, over a type that
    /// happens to know its own derivative.
    fn dual(&self, p: [Dual; 3]) -> Result<Dual, Error> {
        let mut value = self.base.dual(p)?;
        for brush in self.brushes {
            let shape = brush.shape.dual(p)?;
            value = apply_dual(brush.op, value, shape);
        }
        Ok(value)
    }
}

/// A 64-bit LCG, so the 64 brushes are the same 64 brushes on every run.
struct Lcg(u64);

impl Lcg {
    const fn new(seed: u64) -> Self {
        Self(seed)
    }

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

    fn range(&mut self, lo: f64, hi: f64) -> f64 {
        lo + (hi - lo) * self.unit()
    }
}

/// P-39's tape, seed and all, so the two experiments describe the same world.
fn tape() -> Vec<Brush<Shape>> {
    let mut rng = Lcg::new(0x39_5EED_C0DE_1234);
    let mut out = Vec::with_capacity(BRUSHES);
    for _ in 0..BRUSHES {
        let centre = [
            rng.range(-6.5, 6.5),
            rng.range(-6.5, 6.5),
            rng.range(-6.5, 6.5),
        ];
        let shape = if rng.next_u32() & 1 == 0 {
            Shape::Sphere(Sphere {
                center: centre,
                radius: rng.range(0.35, 1.1),
            })
        } else {
            let dir = [
                rng.range(-1.0, 1.0),
                rng.range(-1.0, 1.0),
                rng.range(-1.0, 1.0),
            ];
            let len = (dir[0] * dir[0] + dir[1] * dir[1] + dir[2] * dir[2]).sqrt();
            let unit = if len > 1e-9 {
                [dir[0] / len, dir[1] / len, dir[2] / len]
            } else {
                [1.0, 0.0, 0.0]
            };
            let half = rng.range(0.25, 1.0);
            Shape::Capsule(Capsule {
                a: [
                    centre[0] - unit[0] * half,
                    centre[1] - unit[1] * half,
                    centre[2] - unit[2] * half,
                ],
                b: [
                    centre[0] + unit[0] * half,
                    centre[1] + unit[1] * half,
                    centre[2] + unit[2] * half,
                ],
                radius: rng.range(0.3, 0.8),
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

/// The same tape with every `Add` replaced by a `SmoothAdd`.
fn smooth_tape(hard: &[Brush<Shape>]) -> Vec<Brush<Shape>> {
    hard.iter()
        .map(|b| match b.op {
            BrushOp::Add => Brush::smooth_add(b.shape, SMOOTH_K),
            _ => *b,
        })
        .collect()
}

// ── normals ─────────────────────────────────────────────────────────────────

/// `marching_cubes::unit_gradient`, reproduced exactly — same `len.recip()`.
///
/// # Errors
///
/// [`Error::DegenerateNormal`] where the gradient vanishes.
fn unit(g: [f64; 3]) -> Result<[f64; 3], Error> {
    let len = (g[0] * g[0] + g[1] * g[1] + g[2] * g[2]).sqrt();
    if len.is_nan() || len <= 0.0 {
        return Err(Error::DegenerateNormal { vertex: 0 });
    }
    let r = len.recip();
    Ok([g[0] * r, g[1] * r, g[2] * r])
}

/// The normal the crate produces today: six samples, `O(h²)`.
fn central_normal<F: Sdf<Scalar = f64>>(field: &F, p: [f64; 3]) -> Result<[f64; 3], Error> {
    unit(field.gradient(p))
}

/// The normal one traversal produces, and the value it produced on the way.
fn dual_normal<F: DualSdf>(field: &F, p: [f64; 3]) -> Result<([f64; 3], f64), Error> {
    let d = field.dual(Dual::seed(p))?;
    Ok((unit(d.d)?, d.v))
}

/// Angle between two unit vectors, in degrees.
///
/// `2·asin(|a − b|/2)` rather than `acos(a·b)`: near zero the dot product is
/// `1 − θ²/2` and `acos` throws away half the significant digits, which is
/// precisely the regime a clause about a tenth of a degree lives in.
fn angle_deg(a: [f64; 3], b: [f64; 3]) -> f64 {
    let d = [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
    let chord = (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt();
    2.0 * (chord * 0.5).min(1.0).asin() * (180.0 / std::f64::consts::PI)
}

// ── vertices ────────────────────────────────────────────────────────────────

/// One grid to mesh.
struct Grid {
    shape: RuntimeShape3,
    origin: [f64; 3],
    cell: f64,
}

/// The vertices `MarchingCubes` produces, with the normals it wrote for them.
struct Vertices {
    positions: Vec<[f64; 3]>,
    /// The crate's own emitted normals, kept so the central arm can be checked
    /// against them rather than merely resembling them.
    normals: Vec<[f64; 3]>,
    triangles: usize,
}

/// Mesh every grid and concatenate the vertices.
fn vertices_of<F: Sdf<Scalar = f64>>(field: &F, grids: &[Grid]) -> Vertices {
    let mut mc = MarchingCubes::<f64>::new();
    let mut buf = MeshBuffer::<f64>::new();
    let mut out = Vertices {
        positions: Vec::new(),
        normals: Vec::new(),
        triangles: 0,
    };
    for grid in grids {
        buf.reset();
        mc.extract(field, &grid.shape, grid.origin, grid.cell, &mut buf)
            .expect("extraction");
        out.positions.extend_from_slice(&buf.positions);
        out.normals.extend_from_slice(&buf.normals);
        out.triangles += buf.triangle_count();
    }
    out
}

/// The 64 chunks of the sculpted world.
fn world_grids(layout: &ChunkLayout<f64>) -> Vec<Grid> {
    let mut grids = Vec::with_capacity(64);
    for z in 0..CHUNKS_PER_AXIS {
        for y in 0..CHUNKS_PER_AXIS {
            for x in 0..CHUNKS_PER_AXIS {
                let id = ChunkId::new([x, y, z]);
                grids.push(Grid {
                    shape: layout.sample_shape().expect("chunk grid fits u32"),
                    origin: layout.sample_origin(id),
                    cell: layout.cell_size(),
                });
            }
        }
    }
    grids
}

/// One grid spanning `[-half, half]³` at [`SINGLE_GRID_SAMPLES`] per axis.
fn cube_grid(half: f64) -> Grid {
    Grid {
        shape: RuntimeShape3::new([SINGLE_GRID_SAMPLES; 3]).expect("grid fits u32"),
        origin: [-half; 3],
        cell: 2.0 * half / f64::from(SINGLE_GRID_SAMPLES - 1),
    }
}

// ── timing ──────────────────────────────────────────────────────────────────

/// Time both arms rep by rep, alternating.
///
/// A sibling agent may be benching on the same machine. Alternating means a
/// machine-wide slowdown lands on both arms in the same rep and cancels in the
/// ratio, which a block of five followed by a block of five would not.
///
/// Both arms return the componentwise sum of every normal they computed, and
/// that sum is `black_box`ed, so nothing can be elided and neither arm is
/// charged for a store the other does not make.
fn paired_medians(
    central: &mut impl FnMut(&[[f64; 3]]) -> [f64; 3],
    dual: &mut impl FnMut(&[[f64; 3]]) -> [f64; 3],
    points: &[[f64; 3]],
) -> (f64, f64) {
    // One untimed warm-up per arm.
    black_box(central(black_box(points)));
    black_box(dual(black_box(points)));

    let n = points.len() as f64;
    let mut c = [0.0f64; REPS];
    let mut d = [0.0f64; REPS];
    for (cs, ds) in c.iter_mut().zip(d.iter_mut()) {
        let t = Instant::now();
        let a = central(black_box(points));
        *cs = t.elapsed().as_secs_f64() * 1e9 / n;
        black_box(a);

        let t = Instant::now();
        let b = dual(black_box(points));
        *ds = t.elapsed().as_secs_f64() * 1e9 / n;
        black_box(b);
    }
    c.sort_by(f64::total_cmp);
    d.sort_by(f64::total_cmp);
    (c[REPS / 2], d[REPS / 2])
}

/// The machine's one-minute load average, or `unknown`.
///
/// On the artefact because a concurrent sibling benchmark is a real condition of
/// this measurement and a reader deserves to see it rather than guess.
fn loadavg() -> String {
    std::fs::read_to_string("/proc/loadavg")
        .ok()
        .and_then(|s| s.split_whitespace().next().map(str::to_string))
        .unwrap_or_else(|| String::from("unknown"))
}

// ── one fixture ─────────────────────────────────────────────────────────────

/// Everything one fixture contributed.
struct Row {
    name: &'static str,
    /// What the `central_ns_per_normal` arm actually is for this fixture.
    ///
    /// Load-bearing labelling, not decoration. `Sphere` **does** override
    /// `Sdf::gradient`, so on that row the column is one analytic evaluation and
    /// not six samples; reading it as a central difference would turn the
    /// harness's own control into a false negative on clause two.
    arm: &'static str,
    brushes: usize,
    /// Vertices where both arms produced a normal — the population every number
    /// in the row is over.
    vertices: usize,
    /// Vertices `MarchingCubes` produced, before the degenerate ones were
    /// excluded.
    meshed: usize,
    triangles: usize,
    degenerate: usize,
    mean_deg: f64,
    /// The mean over the vertices that are *not* on a seam.
    ///
    /// Separated because the plain mean turned out to be one outlier divided by
    /// the vertex count, which is a fact about where vertices landed rather
    /// than about the field.
    bulk_mean_deg: f64,
    median_deg: f64,
    p99_deg: f64,
    p999_deg: f64,
    max_deg: f64,
    /// Vertices above 1e-3, 1e-2, 1e-1, 1 and 5 degrees, in that order.
    over: [usize; 5],
    /// Where the worst vertex is, so the claim can be checked by hand.
    worst_at: [f64; 3],
    /// The differencing step `Sdf::gradient` would use at the worst vertex:
    /// `DIFF_STEP · max(|p|, 1)`.
    worst_diff_step: f64,
    /// At the worst vertex, the smallest margin by which any hard `min`/`max`
    /// in the fold was decided.
    ///
    /// `None` where the fixture has no hard seam to be near. Compared against
    /// [`Row::seam_reach`], not against the step itself: the *margin* between
    /// two 1-Lipschitz branches is 2-Lipschitz, so a positional step of `h`
    /// moves it by up to `2h` and the stencil crosses the seam whenever the
    /// margin is under that.
    worst_seam_gap: Option<f64>,
    central_ns: f64,
    dual_ns: f64,
    central_matches_mesh: bool,
    values_bit_exact: bool,
}

impl Row {
    fn speedup(&self) -> f64 {
        self.central_ns / self.dual_ns
    }

    /// How far, in field value, the six-sample stencil can move a seam margin.
    ///
    /// Twice the differencing step: each branch is 1-Lipschitz, so their
    /// difference is 2-Lipschitz. Getting this factor wrong is the difference
    /// between a diagnostic that names the mechanism and one that denies it.
    fn seam_reach(&self) -> f64 {
        2.0 * self.worst_diff_step
    }
}

fn measure<F: DualSdf>(
    name: &'static str,
    arm: &'static str,
    brushes: usize,
    field: &F,
    grids: &[Grid],
) -> Row {
    let meshed = vertices_of(field, grids);

    // Accuracy first, untimed, and it decides the population: a vertex where
    // either arm has no normal is excluded from both rather than given one.
    let mut points = Vec::with_capacity(meshed.positions.len());
    let mut angles = Vec::with_capacity(meshed.positions.len());
    let mut degenerate = 0usize;
    let mut central_matches_mesh = true;
    let mut values_bit_exact = true;
    for (i, p) in meshed.positions.iter().enumerate() {
        let (Ok(c), Ok((g, value))) = (central_normal(field, *p), dual_normal(field, *p)) else {
            degenerate += 1;
            continue;
        };
        // The central arm must be the crate's own output, bit for bit.
        if c.iter()
            .zip(&meshed.normals[i])
            .any(|(a, b)| a.to_bits() != b.to_bits())
        {
            central_matches_mesh = false;
        }
        if value.to_bits() != field.sample(*p).to_bits() {
            values_bit_exact = false;
        }
        angles.push(angle_deg(c, g));
        points.push(*p);
    }

    let mut central_pass = |pts: &[[f64; 3]]| {
        let mut acc = [0.0f64; 3];
        for p in pts {
            let n = central_normal(field, *p).expect("central normal exists on this population");
            acc[0] += n[0];
            acc[1] += n[1];
            acc[2] += n[2];
        }
        acc
    };
    let mut dual_pass = |pts: &[[f64; 3]]| {
        let mut acc = [0.0f64; 3];
        for p in pts {
            let (n, _) = dual_normal(field, *p).expect("dual normal exists on this population");
            acc[0] += n[0];
            acc[1] += n[1];
            acc[2] += n[2];
        }
        acc
    };
    let (central_ns, dual_ns) = paired_medians(&mut central_pass, &mut dual_pass, &points);

    let n = angles.len();
    let mean = angles.iter().sum::<f64>() / n as f64;
    let mut over = [0usize; 5];
    for (slot, threshold) in over.iter_mut().zip(THRESHOLDS_DEG) {
        *slot = angles.iter().filter(|a| **a > threshold).count();
    }
    // The bulk is everything under the smallest threshold a seam can reach.
    let bulk: Vec<f64> = angles.iter().copied().filter(|a| *a <= BULK_DEG).collect();
    let bulk_mean = if bulk.is_empty() {
        0.0
    } else {
        bulk.iter().sum::<f64>() / bulk.len() as f64
    };
    let worst_index = angles
        .iter()
        .enumerate()
        .max_by(|a, b| a.1.total_cmp(b.1))
        .map_or(0, |(i, _)| i);
    let worst_at = points[worst_index];
    let mut sorted = angles;
    sorted.sort_by(f64::total_cmp);
    Row {
        name,
        arm,
        brushes,
        vertices: n,
        meshed: meshed.positions.len(),
        triangles: meshed.triangles,
        degenerate,
        mean_deg: mean,
        bulk_mean_deg: bulk_mean,
        median_deg: sorted[n / 2],
        p99_deg: sorted[(n * 99) / 100],
        p999_deg: sorted[(n * 999) / 1000],
        max_deg: sorted[n - 1],
        over,
        worst_at,
        worst_diff_step: <f64 as Real>::DIFF_STEP
            * worst_at[0]
                .abs()
                .max(worst_at[1].abs())
                .max(worst_at[2].abs())
                .max(1.0),
        worst_seam_gap: None,
        central_ns,
        dual_ns,
        central_matches_mesh,
        values_bit_exact,
    }
}

/// The smallest margin by which any hard `min`/`max` in the fold was decided at
/// `p`.
///
/// The diagnostic that turns "one weird vertex" into a mechanism. A hard
/// `Add` or `Subtract` is a kink: the field is the lower (or upper) of two
/// smooth branches and the gradient jumps across the crossing. A central
/// difference is only wrong there if its stencil *reaches* the crossing, so the
/// quantity that matters is how close `p` is to a tie, measured in field value,
/// against the differencing step. `SmoothAdd` is excluded: it has no kink to be
/// near, which is the entire point of it.
///
/// `f64::INFINITY` when the tape contains no hard operation.
fn seam_gap(base: &Sphere<f64>, tape: &[Brush<Shape>], p: [f64; 3]) -> f64 {
    let mut value = base.sample(p);
    let mut gap = f64::INFINITY;
    for brush in tape {
        let shape = brush.shape.sample(p);
        match brush.op {
            BrushOp::Add => gap = gap.min((value - shape).abs()),
            BrushOp::Subtract => gap = gap.min((value + shape).abs()),
            BrushOp::SmoothAdd { .. } => {}
        }
        value = isomesh::brush::apply(brush.op, value, shape);
    }
    gap
}

/// Clause three: the dual normal against `p/|p|` on a sphere at the origin.
///
/// The largest Euclidean distance between the two unit vectors, over the
/// sphere's own meshed vertices. An instrument that cannot reproduce a closed
/// form has no business measuring a tape.
fn sphere_control(sphere: &Sphere<f64>, grids: &[Grid]) -> f64 {
    let meshed = vertices_of(sphere, grids);
    let mut worst = 0.0f64;
    for p in &meshed.positions {
        let (g, _) = dual_normal(sphere, *p).expect("a sphere's gradient never vanishes");
        let len = (p[0] * p[0] + p[1] * p[1] + p[2] * p[2]).sqrt();
        let exact = [p[0] / len, p[1] / len, p[2] / len];
        let d = [g[0] - exact[0], g[1] - exact[1], g[2] - exact[2]];
        worst = worst.max((d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt());
    }
    worst
}

// ── output ──────────────────────────────────────────────────────────────────

fn row_of(r: &Row, control: f64, load: &str) -> Vec<(&'static str, String)> {
    vec![
        // Registered.
        ("fixture", r.name.to_string()),
        ("brushes", r.brushes.to_string()),
        ("vertices", r.vertices.to_string()),
        // Scientific notation, not fixed point: the bulk error is around 1e-9
        // degrees and a `%.6f` column would print it as `0.000000`, which is a
        // column that hides its own result.
        ("mean_angular_error_deg", format!("{:.6e}", r.mean_deg)),
        ("max_angular_error_deg", format!("{:.6e}", r.max_deg)),
        ("central_ns_per_normal", format!("{:.4}", r.central_ns)),
        ("dual_ns_per_normal", format!("{:.4}", r.dual_ns)),
        ("speedup", format!("{:.4}", r.speedup())),
        ("sphere_control_max_error", format!("{control:.3e}")),
        // Extra: the distribution behind the mean and the max, which is the
        // whole story of this experiment.
        (
            "bulk_mean_angular_error_deg",
            format!("{:.6e}", r.bulk_mean_deg),
        ),
        ("median_angular_error_deg", format!("{:.6e}", r.median_deg)),
        ("p99_angular_error_deg", format!("{:.6e}", r.p99_deg)),
        ("p999_angular_error_deg", format!("{:.6e}", r.p999_deg)),
        ("vertices_over_0p001deg", r.over[0].to_string()),
        ("vertices_over_0p01deg", r.over[1].to_string()),
        ("vertices_over_0p1deg", r.over[2].to_string()),
        ("vertices_over_1deg", r.over[3].to_string()),
        ("vertices_over_5deg", r.over[4].to_string()),
        (
            "fraction_over_0p1deg",
            format!("{:.3e}", r.over[2] as f64 / r.vertices as f64),
        ),
        (
            "worst_vertex",
            format!(
                "{:.6}|{:.6}|{:.6}",
                r.worst_at[0], r.worst_at[1], r.worst_at[2]
            ),
        ),
        // Extra: why the worst vertex is the worst vertex.
        ("gradient_arm", r.arm.to_string()),
        (
            "worst_seam_gap",
            r.worst_seam_gap
                .map_or_else(|| String::from("n/a"), |g| format!("{g:.3e}")),
        ),
        ("worst_diff_step", format!("{:.3e}", r.worst_diff_step)),
        ("worst_seam_reach", format!("{:.3e}", r.seam_reach())),
        (
            "worst_gap_over_reach",
            r.worst_seam_gap.map_or_else(
                || String::from("n/a"),
                |g| format!("{:.3e}", g / r.seam_reach()),
            ),
        ),
        (
            "worst_stencil_straddles_seam",
            r.worst_seam_gap
                .map_or_else(|| String::from("n/a"), |g| (g < r.seam_reach()).to_string()),
        ),
        // Extra: the population, and what was excluded from it.
        ("vertices_meshed", r.meshed.to_string()),
        ("degenerate_points", r.degenerate.to_string()),
        ("triangles", r.triangles.to_string()),
        // Extra: instrument checks. Both must be true for the row to mean
        // anything, and they are on the artefact so that is checkable later.
        (
            "central_matches_mesh_normals",
            r.central_matches_mesh.to_string(),
        ),
        ("dual_value_bit_exact", r.values_bit_exact.to_string()),
        (
            "sphere_control_within_1e-12",
            (control <= SPHERE_CONTROL_TOLERANCE).to_string(),
        ),
        // Extra: conditions.
        ("loadavg_1min", load.to_string()),
        ("reps", REPS.to_string()),
    ]
}

fn main() {
    let prereg = isomesh::experiment!("P-47");

    let base = Sphere::<f64> {
        center: [0.0; 3],
        radius: BASE_RADIUS,
    };
    let hard = tape();
    let smooth = smooth_tape(&hard);
    let layout = ChunkLayout::new(CELLS_PER_CHUNK, CELL_SIZE, [WORLD_ORIGIN; 3])
        .expect("chunk layout is well formed");
    let world = world_grids(&layout);

    // Clause three first: if the arithmetic cannot reproduce a closed form,
    // nothing below is worth printing.
    let unit_sphere = Sphere::<f64>::canonical();
    let single = [cube_grid(2.0)];
    let control = sphere_control(&unit_sphere, &single);
    println!(
        "C3 sphere control: max |n_dual − p/|p|| = {control:.3e} (tolerance {SPHERE_CONTROL_TOLERANCE:.0e}) → {}",
        if control <= SPHERE_CONTROL_TOLERANCE {
            "HELD"
        } else {
            "FALSIFIED"
        }
    );
    assert!(
        control <= SPHERE_CONTROL_TOLERANCE,
        "the dual arithmetic does not reproduce a sphere's normal; every other \
         number in this experiment would be uninterpretable"
    );

    let load = loadavg();
    println!("loadavg (1 min) at start: {load}\n");

    let smooth_stack = BrushStack {
        base,
        brushes: &smooth,
    };
    let lone_capsule = Capsule::<f64> {
        a: [-0.6, 0.0, 0.0],
        b: [0.6, 0.0, 0.0],
        radius: 0.5,
    };

    // The registered fixture is the last rung of a tape-length ladder rather
    // than a lone point: one dual traversal replaces six scalar ones, so where
    // it starts paying is a function of how much there is to traverse, and that
    // crossover is worth having on the artefact.
    let mut rows: Vec<Row> = TAPE_LENGTHS
        .iter()
        .zip(TAPE_NAMES)
        .map(|(len, name)| {
            let stack = BrushStack {
                base,
                brushes: &hard[..*len],
            };
            let mut row = measure(name, "central_difference", *len, &stack, &world);
            row.worst_seam_gap = Some(seam_gap(&base, &hard[..*len], row.worst_at));
            row
        })
        .collect();
    let mut smooth_row = measure(
        "brush_stack_64_smooth",
        "central_difference",
        BRUSHES,
        &smooth_stack,
        &world,
    );
    smooth_row.worst_seam_gap = Some(seam_gap(&base, &smooth, smooth_row.worst_at));
    rows.push(smooth_row);
    // `Capsule` does not override `gradient`, so its arm is the real thing.
    rows.push(measure(
        "capsule",
        "central_difference",
        0,
        &lone_capsule,
        &single,
    ));
    // `Sphere` does. This row is the dual against a hand-written analytic
    // gradient, which is a different and much harder comparison.
    rows.push(measure(
        "sphere",
        "analytic_override",
        0,
        &unit_sphere,
        &single,
    ));

    for r in &rows {
        println!(
            "{:>22}  {:>6} vertices  bulk mean {:9.3e}° p99 {:9.3e}° max {:9.3e}°  \
             over 0.1° {:>3}  {:9.2} → {:8.2} ns/normal  ×{:.3}",
            r.name,
            r.vertices,
            r.bulk_mean_deg,
            r.p99_deg,
            r.max_deg,
            r.over[2],
            r.central_ns,
            r.dual_ns,
            r.speedup(),
        );
        if let Some(gap) = r.worst_seam_gap {
            println!(
                "{:>22}  worst vertex [{:.4}, {:.4}, {:.4}]  seam gap {:.3e} vs stencil reach \
                 {:.3e} (ratio {:.2}) → straddles the seam: {}",
                "",
                r.worst_at[0],
                r.worst_at[1],
                r.worst_at[2],
                gap,
                r.seam_reach(),
                gap / r.seam_reach(),
                gap < r.seam_reach(),
            );
        }
    }

    let registered = rows
        .iter()
        .find(|r| r.name == "brush_stack_64")
        .expect("the registered fixture is in the sweep");
    println!(
        "\nC1 accuracy on {} : mean {:.5}° (> 0.1° required) and max {:.5}° (> 5° required) → {}",
        registered.name,
        registered.mean_deg,
        registered.max_deg,
        if registered.mean_deg > 0.1 && registered.max_deg > 5.0 {
            "HELD"
        } else {
            "FALSIFIED"
        }
    );
    println!(
        "C2 speed on {} : {:.2} → {:.2} ns/normal, ×{:.4} (>= 2x required) → {}",
        registered.name,
        registered.central_ns,
        registered.dual_ns,
        registered.speedup(),
        if registered.speedup() >= 2.0 {
            "HELD"
        } else {
            "FALSIFIED"
        }
    );
    println!("loadavg (1 min) at end: {}", loadavg());

    common::experiment::run(prereg, |run| {
        for r in &rows {
            run.record(&row_of(r, control, &load));
        }
    });
}
