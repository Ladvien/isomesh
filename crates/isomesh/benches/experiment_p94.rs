//! **P-94 — how big is a dig, in bytes.**
//!
//! Ticket: R-094. Pre-registered before this harness existed.
//!
//! ```bash
//! cargo bench --bench experiment_p94
//! ```
//!
//! Writes `docs/experiments/p-94.csv`.
//!
//! # Scope, restated from the registration
//!
//! This measures the **edit log** — what a save file is and what undo rewinds.
//! No sockets, no clock, no session model, and `BACKLOG.md`'s closure of
//! networked editing is not reopened. The only public datum in the world is
//! Gustafsson's *The unlikely story of Teardown Multiplayer*
//! (`blog.voxagon.se`, 13 Mar 2026) — **~1 Mbit per client**, testimony rather
//! than measurement, whose stated reason for replicating *commands* is
//! *"commands are the same regardless of object size"*. The 08-11 doc's
//! **100k operations × 48 B = 4.80 MB** is an estimate. Both are quoted here as
//! comparisons and neither is used as an input.
//!
//! # SHARE, recomputed before this harness was written
//!
//! The registration's SHARE line is *"this measures sizes and moves nothing"*,
//! so there is no fraction-of-a-runtime to be denominated in — but two of the
//! three clauses have an arithmetic reachability question anyway, and one of them
//! is **unreachable**:
//!
//! - **C1 is reachable with margin.** `EDIT_PERIOD` is 0.08 s but the demo's
//!   throttle does not catch up, so at 60 Hz the period quantises to five frames
//!   = 1/12 s and an hour offers **43,200** dig attempts rather than 45,000. The
//!   crate's record is a shape plus an op: serialised as fixed point that is
//!   **17 B**, in memory `size_of::<Brush<Sphere<f32>>>()`. 43,200 × 17 B =
//!   734 kB against a 2 MB bar — 2.7× of headroom. Note what the *estimate*
//!   would have given: 43,200 × 48 B = **2.07 MB**, which fails C1. So C1 is a
//!   claim about the record's width, and the 08-11 doc's width would have lost.
//! - **C1's second half is structural rather than measurable.** `Brush<S>` is a
//!   shape and an op. It contains no chunk reference and no grid index, so per
//!   edit cost *cannot* depend on chunk granularity, and its coordinates are
//!   fixed-width so it cannot depend on world size either. A harness that
//!   asserted this would be asserting the definition of the type. What is done
//!   instead: the whole hour is re-flown in five sandboxes and the log's bytes
//!   are measured from the trace each one produced, with `trace_span_cells` and
//!   `dirty_chunks_per_edit` on the row to show **both knobs moved something
//!   else in the same harness** — which is the only honest way to report a
//!   constant.
//! - **C2 is arithmetically unreachable, and this is written down before the
//!   run.** A capsule is a segment dilated by `radius`. Coaxial capsules
//!   *partition* the axis: capsule `i` is the only brush whose region reaches
//!   full radius over its own axial slab `[p_i, p_{i+1}]`, because a point at
//!   axial offset `u` past `p_i` and radial offset `ρ` is inside the previous
//!   capsule only when `ρ ≤ √(r² − u²)`. Dropping one leaves an uncovered lens
//!   of radial thickness `r − √(r² − (δ/2)²)`; dropping a run of `k` leaves
//!   `r − √(max(0, r² − (kδ/2)²))`, which reaches the **full radius** as soon as
//!   `kδ ≥ 2r`. 200 → 40 means keeping every fifth, so `k = 4`; the demo's own
//!   distance gate is `δ ≥ r/2`, giving `kδ ≥ 2r` exactly — a solid disc of rock
//!   of radius `r` plugging the tunnel. **No sound compressor can reach 40 by
//!   dropping**, at any grid resolution, and the run is here to produce the
//!   number rather than to discover that. What *does* collapse a coaxial run is
//!   a **merge**, and that is measured as well: collinear contiguous segments
//!   have `min_i d(p, seg_i) = d(p, seg_hull)` identically, so the hull capsule
//!   is the same set — exact in ℝ, and only rounding away in `f64`.
//! - **C3 is reachable and is a genuine unknown.** Both sides compress the same
//!   byte stream, so the share is the whole log.
//!
//! # The trace, and why it is a session rather than a sweep
//!
//! The registered vacuity control has two halves and this is the first. The
//! trace is a **replay of `bevy_isomesh/examples/game_dig.rs`'s dig loop**, with
//! its constants and its two gates transcribed rather than paraphrased:
//! `EDIT_PERIOD = 0.08` (game_dig.rs:342), the half-radius distance gate
//! (game_dig.rs:2635), `AIM_NEAR = 0.30`, `AIM_FAR = 25.0`, `AIM_STEPS = 128`,
//! `AIM_HIT = 0.01`, `LIPSCHITZ = 1.25` (game_dig.rs:216–226, 190), the
//! `[0.10, 2.00]` brush clamp (game_dig.rs:2581), `CHUNK_CELLS = 16` and
//! `CELL_SIZE = 0.125` (game_dig.rs:125–127), the sphere-tracing `trace`
//! function (game_dig.rs:1373–1397) including its *"already inside rock, return
//! `AIM_NEAR`, that is how you dig yourself out"* behaviour, and `Ground`
//! itself — `y − (0.35·sin(0.9x)·cos(0.7z) + 0.15·sin(2.1x))`, whose 1.207
//! gradient bound is what makes 1.25 a legal step divisor.
//!
//! **What makes it a session and not `AutoCarve::centre(n)`** is that no brush
//! centre is scripted. The centre is where the ray *stops*, and the ray stops on
//! the surface of the field **as already edited** — so the tunnel advances
//! because the wall receded, and the spacing `δ` between consecutive brushes is
//! an output of the field rather than an input. That is instrumented, not
//! asserted: `closed_loop_edits` counts edits whose hit distance exceeds the hit
//! distance the *unedited* `Ground` would have returned from the same eye and
//! the same direction. An edit can only be counted there if a previous edit
//! moved the surface. It is asserted non-zero, and a synthetic sweep scores
//! **zero on it by construction** — which is why the same number is also
//! reported for an explicit `AutoCarve`-shaped straight-line control arm, so the
//! session/sweep distinction is a measured difference rather than a claim.
//!
//! Two deviations, both declared. The eye is **flown** rather than walked:
//! game_dig itself flies during a scripted capture (game_dig.rs:1074), and the
//! log depends on the eye and the view direction, which fly mode supplies, not
//! on which of four body spheres resolved first. And the aim is traced once per
//! throttle period rather than once per rendered frame, because only the aim at
//! a stroke instant can enter the log; the consequence is that an edit refused
//! by the distance gate waits a whole period instead of a frame, which can only
//! *lower* the edit count — the conservative direction for a size bound.
//!
//! # The compressors
//!
//! The canonical record is **fixed point**, which is Teardown's own account of
//! what they shipped: op in one byte, centre and radius in `i32` at a quantum of
//! 1/4096 of a world unit — 512 steps per cell at `CELL_SIZE`, far finer than
//! any extractor here can express. Fixed point matters for fairness rather than
//! for size: `i32` and `f32` are both four bytes, so `bytes_per_edit` is 17
//! either way, and quantising *before* both compressors is what stops the
//! bespoke coder from winning on precision `zstd` was made to keep.
//!
//! - `bytes_zstd` is `zstd` at **level 19** on that record stream. The strongest
//!   setting, deliberately: C3's bar is easier the weaker the baseline, and
//!   level 3 is on the row as well so the choice is visible.
//! - `bytes_entropy_coded` is an **adaptive binary range coder** (LZMA's, 11-bit
//!   probabilities, 5-bit adaptation) over a **second-order linear predictor** —
//!   `pred_i = 2q_{i-1} − q_{i-2}` per axis, then zigzag, then Elias-gamma with
//!   an adaptive context per bit-length and per bit position. Second order
//!   because the registration's stated mechanism is that *brush parameters are
//!   strongly correlated along a stroke*, and a stroke is a nearly constant
//!   velocity along a ray. It is **decoded back and asserted equal** to the
//!   input, so the number is the length of a real encoding rather than an
//!   entropy figure nobody can invert.
//! - `bytes_zstd_of_deltas` is the column that decides whether a bespoke format
//!   is worth writing at all: `zstd` on the *same* second-order residuals. If it
//!   lands near the range coder, the advantage was the model and not the coder,
//!   and the answer is "keep zstd, transform first".
//!
//! # Controls, and what each one could have said
//!
//! - `overlapping_pairs` — the registered second half of the vacuity control:
//!   exact segment-to-segment distance under `r_i + r_j` over all 19,900 pairs
//!   of the coaxial arm. Asserted non-zero. Zero would mean C2 was compressing
//!   disjoint brushes.
//! - `closed_loop_edits` — asserted non-zero, above.
//! - `entropy_roundtrip_ok` — asserted. The decoder rebuilds every record.
//! - `coaxial_hash_verified` — the greedy drop keeps the **whole tunnel's**
//!   `mesh_hash` bit-identical to the 200-brush fold, checked once over the full
//!   grid after the greedy pass rather than trusted from the per-candidate
//!   tests. This is ✗41's `necessary_only_hash_unchanged` and it is the reason
//!   `surviving_brushes` is a measurement and not a heuristic.
//!
//! No clause's own predicate is asserted anywhere in this file. The controls are
//! the fixture's, not the hypothesis's.

#![allow(
    clippy::cast_possible_truncation,
    clippy::cast_precision_loss,
    clippy::cast_sign_loss,
    clippy::cast_possible_wrap,
    clippy::too_many_lines,
    clippy::similar_names,
    clippy::print_literal
)]

mod common;

use isomesh::brush::{Brush, BrushOp, Capsule, apply};
use isomesh::chunk::ChunkLayout;
use isomesh::fields::Sphere;
use isomesh::marching_cubes::MarchingCubes;
use isomesh::validate::mesh_hash;
use isomesh::{MeshBuffer, RuntimeShape3, Sdf};

// ── constants transcribed from bevy_isomesh/examples/game_dig.rs ─────────────

/// game_dig.rs:342. Seconds between edits while a button is held.
const EDIT_PERIOD: f64 = 0.08;
/// game_dig.rs:216. Nearest a brush may be placed.
const AIM_NEAR: f64 = 0.30;
/// game_dig.rs:220. Furthest the trace looks.
const AIM_FAR: f64 = 25.0;
/// game_dig.rs:223. Iteration cap on the sphere trace.
const AIM_STEPS: u32 = 128;
/// game_dig.rs:226. Surface tolerance, and the minimum step.
const AIM_HIT: f64 = 0.01;
/// game_dig.rs:190. `Ground`'s gradient bound is 1.207; this is the divisor.
const LIPSCHITZ: f64 = 1.25;
/// game_dig.rs:127.
const CELL_SIZE: f64 = 0.125;
/// game_dig.rs:125.
const CHUNK_CELLS: u32 = 16;
/// game_dig.rs:1056. Default brush radius.
const RADIUS_DEFAULT: f64 = 0.25;
/// Frame step. 60 Hz, which is what quantises the throttle to five frames.
const FRAME_DT: f64 = 1.0 / 60.0;
/// One hour.
const SESSION_SECONDS: f64 = 3600.0;

/// Largest radius this session's wheel reaches.
///
/// game_dig clamps to `[0.10, 2.00]`; this session uses `[0.25, 0.60]`, and the
/// reason is the accelerator rather than taste: the nearest-brush search's ring
/// bound is `R_MAX − (k−1)·ACCEL_H`, so `R_MAX` decides how many lattice cells
/// every field sample has to look at. A 2.0 brush is a chamber scoop, not a dig,
/// and it would triple the search for a handful of edits. Declared as a
/// deviation; it cannot affect `bytes_per_edit`, which is the record's width.
const R_MAX: f64 = 1.0;

/// Lattice pitch of the nearest-brush accelerator, in world units.
///
/// `2·R_MAX`, which is what makes the ring bound fire at `k = 2`: a brush in a
/// cell at Chebyshev cell-distance ≥ 2 has `r_i − |p − c_i| ≤ R_MAX − ACCEL_H =
/// −1`, below every threshold this harness queries with. So the search visits
/// exactly the 27 cells of rings 0 and 1 and is *exact*, not approximate.
const ACCEL_H: f64 = 2.0;

/// Quantum of the serialised record, in world units. 1/4096.
const QUANTUM: f64 = 1.0 / 4096.0;

/// Serialised record width: op `u8`, centre `3 × i32`, radius `i32`.
const RECORD_BYTES: usize = 1 + 12 + 4;

/// C1's bar, in bytes. Decimal MB, matching the 08-11 doc's `100k × 48 B =
/// 4.80 MB`.
const C1_BAR_BYTES: f64 = 2_000_000.0;

/// The 08-11 doc's estimate, for comparison only.
const ESTIMATE_BYTES: f64 = 4_800_000.0;

/// Teardown's ~1 Mbit per client, in bytes, for comparison only. Testimony.
const TEARDOWN_BYTES: f64 = 1_000_000.0 / 8.0;

/// Capsules in the coaxial arm.
const COAXIAL_BRUSHES: usize = 200;

// ── vector helpers ──────────────────────────────────────────────────────────

fn add(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn scale(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn length(a: [f64; 3]) -> f64 {
    dot(a, a).sqrt()
}

/// Closest distance between two segments. Ericson, *Real-Time Collision
/// Detection*, §5.1.9, with the parallel case handled by clamping.
fn seg_seg_distance(p1: [f64; 3], q1: [f64; 3], p2: [f64; 3], q2: [f64; 3]) -> f64 {
    const EPS: f64 = 1e-15;
    let d1 = sub(q1, p1);
    let d2 = sub(q2, p2);
    let r = sub(p1, p2);
    let a = dot(d1, d1);
    let e = dot(d2, d2);
    let f = dot(d2, r);
    let (s, t);
    if a <= EPS && e <= EPS {
        return length(r);
    } else if a <= EPS {
        s = 0.0;
        t = (f / e).clamp(0.0, 1.0);
    } else {
        let c = dot(d1, r);
        if e <= EPS {
            t = 0.0;
            s = (-c / a).clamp(0.0, 1.0);
        } else {
            let b = dot(d1, d2);
            let denom = a * e - b * b;
            let s0 = if denom > 0.0 {
                ((b * f - c * e) / denom).clamp(0.0, 1.0)
            } else {
                0.0
            };
            let t0 = (b * s0 + f) / e;
            if t0 < 0.0 {
                t = 0.0;
                s = (-c / a).clamp(0.0, 1.0);
            } else if t0 > 1.0 {
                t = 1.0;
                s = ((b - c) / a).clamp(0.0, 1.0);
            } else {
                t = t0;
                s = s0;
            }
        }
    }
    length(sub(add(p1, scale(d1, s)), add(p2, scale(d2, t))))
}

// ── the field ───────────────────────────────────────────────────────────────

/// game_dig's `Ground`, transcribed. `height = 0.35·sin(0.9x)·cos(0.7z) +
/// 0.15·sin(2.1x)`, so `|∇f| ≤ √(1 + 0.63² + 0.245²) = 1.207` and
/// [`LIPSCHITZ`] = 1.25 is a legal step divisor.
fn ground(p: [f64; 3]) -> f64 {
    p[1] - (0.35 * (0.9 * p[0]).sin() * (0.7 * p[2]).cos() + 0.15 * (2.1 * p[0]).sin())
}

/// The edit log, with the lattice that makes a nearest-brush query exact and
/// cheap.
struct Log {
    centres: Vec<[f64; 3]>,
    radii: Vec<f64>,
    /// Lattice origin, in world units.
    origin: [f64; 3],
    /// Lattice dimensions, in cells.
    dims: [i32; 3],
    buckets: Vec<Vec<u32>>,
}

impl Log {
    fn new(lo: [f64; 3], hi: [f64; 3]) -> Self {
        let dims = [
            (((hi[0] - lo[0]) / ACCEL_H).ceil() as i32).max(1),
            (((hi[1] - lo[1]) / ACCEL_H).ceil() as i32).max(1),
            (((hi[2] - lo[2]) / ACCEL_H).ceil() as i32).max(1),
        ];
        let n = (dims[0] as usize) * (dims[1] as usize) * (dims[2] as usize);
        Self {
            centres: Vec::new(),
            radii: Vec::new(),
            origin: lo,
            dims,
            buckets: vec![Vec::new(); n],
        }
    }

    fn cell_of(&self, p: [f64; 3]) -> [i32; 3] {
        [
            ((p[0] - self.origin[0]) / ACCEL_H).floor() as i32,
            ((p[1] - self.origin[1]) / ACCEL_H).floor() as i32,
            ((p[2] - self.origin[2]) / ACCEL_H).floor() as i32,
        ]
    }

    fn index(&self, c: [i32; 3]) -> Option<usize> {
        if c[0] < 0
            || c[1] < 0
            || c[2] < 0
            || c[0] >= self.dims[0]
            || c[1] >= self.dims[1]
            || c[2] >= self.dims[2]
        {
            return None;
        }
        Some(
            (c[0] as usize)
                + (self.dims[0] as usize) * ((c[1] as usize) + (self.dims[1] as usize) * (c[2] as usize)),
        )
    }

    fn push(&mut self, centre: [f64; 3], radius: f64) {
        assert!(
            radius <= R_MAX,
            "the accelerator's ring bound is derived from R_MAX and this radius exceeds it: {radius}"
        );
        let id = self.centres.len() as u32;
        self.centres.push(centre);
        self.radii.push(radius);
        let cell = self.cell_of(centre);
        if let Some(i) = self.index(cell) {
            self.buckets[i].push(id);
        }
    }

    /// `maxᵢ (rᵢ − |p − cᵢ|)`, computed exactly whenever it exceeds `floor`, and
    /// otherwise some value that does not.
    ///
    /// Exact because a brush in a lattice cell at Chebyshev distance `k ≥ 1` is
    /// at least `(k−1)·ACCEL_H` away, so its contribution is at most
    /// `R_MAX − (k−1)·ACCEL_H`; with `ACCEL_H = 2·R_MAX` that is negative from
    /// `k = 2` on, and every `floor` this harness passes is at least
    /// [`AIM_HIT`]. Rings 0 and 1 are therefore the whole search.
    fn carve_above(&self, p: [f64; 3], floor: f64) -> f64 {
        debug_assert!(floor >= AIM_HIT);
        let cp = self.cell_of(p);
        let mut best = f64::NEG_INFINITY;
        for dz in -1..=1 {
            for dy in -1..=1 {
                for dx in -1..=1 {
                    let Some(i) = self.index([cp[0] + dx, cp[1] + dy, cp[2] + dz]) else {
                        continue;
                    };
                    for &b in &self.buckets[i] {
                        let b = b as usize;
                        let v = self.radii[b] - length(sub(p, self.centres[b]));
                        if v > best {
                            best = v;
                        }
                    }
                }
            }
        }
        best
    }

    /// The fold, to the precision the sphere trace can use.
    ///
    /// `max(ground, carve)`, exact wherever the value is above [`AIM_HIT`]. Below
    /// that the trace's behaviour is identical for every value — it either
    /// reports a hit or takes its minimum step — so the accelerator is allowed to
    /// stop looking, and that is what bounds it to 27 lattice cells.
    fn march_value(&self, p: [f64; 3]) -> f64 {
        let g = ground(p);
        let floor = g.max(AIM_HIT);
        g.max(self.carve_above(p, floor))
    }
}

/// game_dig.rs:1373. First surface crossing along a ray inside the sandbox.
fn trace(
    value: &impl Fn([f64; 3]) -> f64,
    origin: [f64; 3],
    direction: [f64; 3],
    lo: [f64; 3],
    hi: [f64; 3],
) -> Option<f64> {
    let mut t = AIM_NEAR;
    for _ in 0..AIM_STEPS {
        let p = add(origin, scale(direction, t));
        let f = value(p);
        let inside = (0..3).all(|a| p[a] >= lo[a] && p[a] <= hi[a]);
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

// ── the session ─────────────────────────────────────────────────────────────

/// Deterministic input policy. A hand holding the button for an hour.
///
/// Not a sweep: this drives the *eye and the look*, never a brush centre. Where
/// the brush lands is decided by [`trace`] against the field as already dug.
#[derive(Clone, Copy)]
enum Phase {
    /// Bore straight: look fixed, walk along the view.
    Bore,
    /// Sweep the brush sideways while creeping forward.
    Sweep,
    /// Stand and hollow a chamber.
    Chamber,
    /// Fly to a new heading, still carving.
    Relocate,
}

impl Phase {
    fn next(seed: u64) -> Self {
        match seed % 4 {
            0 => Self::Bore,
            1 => Self::Sweep,
            2 => Self::Chamber,
            _ => Self::Relocate,
        }
    }

    /// Radius the wheel sits at during this phase. Every value ≤ [`R_MAX`].
    fn radius(self) -> f64 {
        match self {
            Self::Bore | Self::Relocate => RADIUS_DEFAULT,
            Self::Sweep => 0.40,
            Self::Chamber => 0.60,
        }
    }

    /// Slot in the per-phase counters.
    fn index(self) -> usize {
        match self {
            Self::Bore => 0,
            Self::Sweep => 1,
            Self::Chamber => 2,
            Self::Relocate => 3,
        }
    }
}

/// xorshift64\*. A deterministic input script needs a deterministic hand.
struct Rng(u64);

impl Rng {
    fn next_u64(&mut self) -> u64 {
        let mut x = self.0;
        x ^= x >> 12;
        x ^= x << 25;
        x ^= x >> 27;
        self.0 = x;
        x.wrapping_mul(0x2545_F491_4F6C_DD1D)
    }

    fn unit(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 / (1u64 << 53) as f64
    }
}

/// Everything one hour produced.
struct SessionResult {
    world_cells: u32,
    /// Dig attempts: frames on which the throttle was due and the gates ran.
    strokes: u64,
    /// Brushes pushed.
    centres: Vec<[f64; 3]>,
    radii: Vec<f64>,
    refusals_distance: u64,
    refusals_nohit: u64,
    /// Edits whose hit distance beat the unedited `Ground`'s.
    closed_loop: u64,
    /// Mean spacing between consecutive brush centres.
    mean_step: f64,
    /// Trace bounding box diagonal, in cells.
    span_cells: f64,
    /// Brushes provably contained in another brush, so exactly redundant under
    /// `max`: `|cᵢ − cⱼ| + rᵢ ≤ rⱼ`.
    nested_dropped: usize,
    distinct_radii: usize,
    /// Edits and no-hit refusals split by phase, in [`Phase`]'s declaration
    /// order. The fixture instrument: it names *which* part of the hand policy
    /// cannot find rock, and the first version of this harness needed it.
    edits_by_phase: [u64; 4],
    nohit_by_phase: [u64; 4],
}

/// The sandbox, in game_dig's proportions: `[world_cells, world_cells/2,
/// world_cells]` cells, centred on the origin.
fn sandbox(world_cells: u32) -> ([f64; 3], [f64; 3]) {
    let ex = f64::from(world_cells) * CELL_SIZE;
    let ey = f64::from(world_cells / 2) * CELL_SIZE;
    (
        [-ex * 0.5, -ey * 0.5, -ex * 0.5],
        [ex * 0.5, ey * 0.5, ex * 0.5],
    )
}

fn direction_of(yaw: f64, pitch: f64) -> [f64; 3] {
    let (sy, cy) = yaw.sin_cos();
    let (sp, cp) = pitch.sin_cos();
    [cp * sy, sp, -cp * cy]
}

fn run_session(world_cells: u32) -> SessionResult {
    let (lo, hi) = sandbox(world_cells);
    let margin = 0.75;
    let mut log = Log::new(lo, hi);

    let mut rng = Rng(0x9E37_79B9_7F4A_7C15 ^ u64::from(world_cells));
    let mut eye = [0.0, 1.0, 0.0];
    let mut yaw = 0.4;
    let mut pitch = -0.45;
    let mut phase = Phase::Bore;
    let mut phase_left = 6.0;
    let mut radius = phase.radius();
    let mut heading = yaw;

    let mut clock = 0.0;
    let mut stroke_last: Option<[f64; 3]> = None;
    let mut strokes = 0u64;
    let mut refusals_distance = 0u64;
    let mut refusals_nohit = 0u64;
    let mut closed_loop = 0u64;
    let mut step_sum = 0.0;
    let mut step_count = 0u64;
    let mut edits_by_phase = [0u64; 4];
    let mut nohit_by_phase = [0u64; 4];
    // Frames left of "the hand has lost the wall and is scanning for it".
    let mut recover = 0u32;

    let frames = (SESSION_SECONDS / FRAME_DT).round() as u64;
    let mut phase_time = 0.0;
    for _ in 0..frames {
        // ── the hand ────────────────────────────────────────────────────────
        phase_left -= FRAME_DT;
        phase_time += FRAME_DT;
        if phase_left <= 0.0 {
            phase = Phase::next(rng.next_u64());
            phase_left = 3.0 + rng.unit() * 9.0;
            phase_time = 0.0;
            radius = phase.radius();
            heading = rng.unit() * std::f64::consts::TAU;
        }
        let speed = match phase {
            Phase::Bore => 1.5,
            Phase::Sweep => 0.4,
            Phase::Chamber => 0.1,
            Phase::Relocate => 4.0,
        };
        // While the hand is recovering a lost wall the phase's own look script
        // is suspended. Without that suspension `Chamber` and `Relocate`
        // rewrite `pitch` every frame and the downward scan below cannot
        // accumulate -- measured 604 edits of a possible 2,160 at 128 cells
        // against 2,014 at 1,024 cells, where `Bore` happened to be running.
        if recover == 0 {
            match phase {
                Phase::Bore => {}
                Phase::Sweep => yaw += 0.5 * FRAME_DT,
                Phase::Chamber => {
                    yaw += 1.0 * FRAME_DT;
                    pitch = -0.2 + 0.6 * (2.0 * phase_time).sin();
                }
                Phase::Relocate => {
                    yaw += (heading - yaw).clamp(-1.5 * FRAME_DT, 1.5 * FRAME_DT);
                    pitch += (-0.25 - pitch).clamp(-FRAME_DT, FRAME_DT);
                }
            }
        } else {
            recover -= 1;
        }
        pitch = pitch.clamp(-1.25, 1.2);
        let dir = direction_of(yaw, pitch);
        // Fly mode, which is what game_dig uses for a scripted capture -- but
        // the **feet follow the crosshair**, not the view direction. A digging
        // player walks into the hole they are making, and the first version of
        // this policy flew along the view instead: the eye left the rock, the
        // ray ran out of sandbox before it reached the ground, and the hour
        // produced 1,067 edits of a possible 43,200 with `gate_refusals_nohit`
        // at 210,661. That is the fixture defect `refusals_nohit` exists to
        // catch, and an unsaturated hour would have cleared C1's 2 MB bar by
        // not digging.
        let step_dir = match stroke_last {
            Some(t) if length(sub(t, eye)) > 1e-9 => {
                let d = sub(t, eye);
                scale(d, 1.0 / length(d))
            }
            _ => dir,
        };
        // Never past the face: the eye stops `AIM_NEAR` short of what it is
        // digging, which is the same bound game_dig's `AIM_NEAR` enforces on the
        // brush.
        let room = stroke_last.map_or(f64::INFINITY, |t| {
            (length(sub(t, eye)) - AIM_NEAR).max(0.0)
        });
        let mut moved_eye = add(eye, scale(step_dir, (speed * FRAME_DT).min(room)));
        // The sandbox is the wall the demo lines with five slabs. Here it turns
        // the hand around rather than being five brushes that never enter the
        // log.
        let mut bounced = false;
        for a in 0..3 {
            let (l, h) = (lo[a] + margin, hi[a] - margin);
            if moved_eye[a] < l {
                moved_eye[a] = l;
                bounced = true;
            } else if moved_eye[a] > h {
                moved_eye[a] = h;
                bounced = true;
            }
        }
        eye = moved_eye;
        if bounced {
            yaw += 2.1;
            pitch = -0.3;
        }

        // ── the two gates, from game_dig's `dig` ────────────────────────────
        clock += FRAME_DT;
        if clock < EDIT_PERIOD {
            continue;
        }
        strokes += 1;
        let value = |p: [f64; 3]| log.march_value(p);
        let Some(t) = trace(&value, eye, dir, lo, hi) else {
            refusals_nohit += 1;
            nohit_by_phase[phase.index()] += 1;
            // A hand that has lost the wall **looks down**, and keeps
            // steepening until it finds rock. The world is a heightfield, so
            // there is always rock below: straight down from open air hits the
            // terrain, and straight down inside a cavity hits its floor.
            //
            // Two earlier versions of this are the fixture defects
            // `gate_refusals_nohit` and the per-phase split caught. Flying
            // along the view direction instead of following the crosshair gave
            // **1,040 edits against 5,596 no-hit refusals** in a three-minute
            // session — the ray ran out of sandbox before it reached the
            // ground. Re-aiming at the *last aim point* was worse, **96
            // edits**: the eye stops `AIM_NEAR` short of that point, so the
            // ray starts inside the cavity that point carved, exits through
            // the hole, finds nothing, and the look is pinned back at the same
            // failing direction for the rest of the hour. A monotone scan
            // cannot deadlock, and the yaw kick stops it retracing one meridian.
            pitch = (pitch - 0.35).max(-1.25);
            yaw += 0.37;
            // Half a second of scanning, refreshed on every failure, so the
            // phase script cannot overwrite the scan mid-recovery.
            recover = 30;
            continue;
        };
        let centre = add(eye, scale(dir, t));
        // The distance gate: game_dig.rs:2635.
        if let Some(last) = stroke_last
            && length(sub(centre, last)) < radius * 0.5
        {
            refusals_distance += 1;
            continue;
        }
        // The closed-loop control: how far the ray would have got on the
        // *unedited* field, from the same eye along the same direction.
        let virgin = |p: [f64; 3]| ground(p);
        if let Some(vt) = trace(&virgin, eye, dir, lo, hi) {
            if t > vt + 1e-9 {
                closed_loop += 1;
            }
        } else {
            closed_loop += 1;
        }
        if let Some(last) = stroke_last {
            step_sum += length(sub(centre, last));
            step_count += 1;
        }
        log.push(centre, radius);
        stroke_last = Some(centre);
        edits_by_phase[phase.index()] += 1;
        clock = 0.0;
    }

    // Exactly redundant under `max`: a sphere inside another sphere.
    let n = log.centres.len();
    let mut nested = 0usize;
    for i in 0..n {
        for j in i.saturating_sub(64)..(i + 64).min(n) {
            if i != j
                && length(sub(log.centres[i], log.centres[j])) + log.radii[i] <= log.radii[j]
            {
                nested += 1;
                break;
            }
        }
    }

    let mut blo = [f64::INFINITY; 3];
    let mut bhi = [f64::NEG_INFINITY; 3];
    for c in &log.centres {
        for a in 0..3 {
            blo[a] = blo[a].min(c[a]);
            bhi[a] = bhi[a].max(c[a]);
        }
    }
    let span_cells = if n == 0 {
        0.0
    } else {
        length(sub(bhi, blo)) / CELL_SIZE
    };

    let mut radii_sorted: Vec<u64> = log.radii.iter().map(|r| r.to_bits()).collect();
    radii_sorted.sort_unstable();
    radii_sorted.dedup();

    SessionResult {
        world_cells,
        strokes,
        refusals_distance,
        refusals_nohit,
        closed_loop,
        mean_step: if step_count == 0 {
            0.0
        } else {
            step_sum / step_count as f64
        },
        span_cells,
        nested_dropped: nested,
        distinct_radii: radii_sorted.len(),
        edits_by_phase,
        nohit_by_phase,
        centres: log.centres,
        radii: log.radii,
    }
}

/// The `AutoCarve::centre(n)`-shaped control: a straight line at a fixed
/// cadence, no field feedback at all. Its `closed_loop` count is zero by
/// construction, which is what makes the session's non-zero count mean
/// something.
fn sweep_control(edits: usize) -> (Vec<[f64; 3]>, Vec<f64>) {
    let mut centres = Vec::with_capacity(edits);
    let mut radii = Vec::with_capacity(edits);
    for n in 0..edits {
        let t = n as f64;
        centres.push([
            -0.9 + t * 0.30 * 0.125,
            0.55 - t * 0.045 * 0.125,
            2.2 - t * 0.34 * 0.125,
        ]);
        radii.push(RADIUS_DEFAULT);
    }
    (centres, radii)
}

// ── the record, and the two compressors ─────────────────────────────────────

fn quantise(v: f64) -> i32 {
    (v / QUANTUM).round() as i32
}

/// The canonical save record: op `u8`, centre `3 × i32`, radius `i32`, fixed
/// point at [`QUANTUM`]. Fixed point is Teardown's account of what they shipped,
/// and here it is what makes the two compressors see the same bytes.
fn serialise(centres: &[[f64; 3]], radii: &[f64]) -> Vec<u8> {
    let mut out = Vec::with_capacity(centres.len() * RECORD_BYTES);
    for (c, r) in centres.iter().zip(radii) {
        // Subtract. This session is a dig; nothing fills.
        out.push(1u8);
        for v in c {
            out.extend_from_slice(&quantise(*v).to_le_bytes());
        }
        out.extend_from_slice(&quantise(*r).to_le_bytes());
    }
    out
}

/// Second-order residuals of the same records, in the same 17-byte layout, so
/// `zstd` sees exactly the transform the range coder sees.
fn serialise_residuals(centres: &[[f64; 3]], radii: &[f64]) -> Vec<u8> {
    let res = residuals(centres, radii);
    let mut out = Vec::with_capacity(res.len() * RECORD_BYTES);
    for r in &res {
        out.push(1u8);
        for v in r {
            out.extend_from_slice(&v.to_le_bytes());
        }
    }
    out
}

/// `[dx, dy, dz, dr]` per edit: second order on the centre, first order on the
/// radius.
fn residuals(centres: &[[f64; 3]], radii: &[f64]) -> Vec<[i32; 4]> {
    let mut out = Vec::with_capacity(centres.len());
    let mut prev = [0i32; 3];
    let mut prev2 = [0i32; 3];
    let mut prev_r = 0i32;
    for (i, (c, r)) in centres.iter().zip(radii).enumerate() {
        let q = [quantise(c[0]), quantise(c[1]), quantise(c[2])];
        let qr = quantise(*r);
        let mut row = [0i32; 4];
        for a in 0..3 {
            let pred = match i {
                0 => 0,
                1 => prev[a],
                _ => prev[a].wrapping_mul(2).wrapping_sub(prev2[a]),
            };
            row[a] = q[a].wrapping_sub(pred);
        }
        row[3] = qr.wrapping_sub(prev_r);
        out.push(row);
        prev2 = prev;
        prev = q;
        prev_r = qr;
    }
    out
}

const PROB_BITS: u32 = 11;
const PROB_ONE: u32 = 1 << PROB_BITS;
const PROB_INIT: u16 = (PROB_ONE / 2) as u16;
const MOVE_BITS: u32 = 5;
const RC_TOP: u32 = 1 << 24;
/// Bit-length contexts. Zigzagged residuals of a 1/4096 fixed point over a
/// 128-unit world need at most 21 bits; 40 is slack, and an assert holds it.
const NBITS: usize = 40;

/// LZMA's binary range encoder.
struct RangeEncoder {
    low: u64,
    range: u32,
    cache: u8,
    cache_size: u64,
    out: Vec<u8>,
}

impl RangeEncoder {
    fn new() -> Self {
        Self {
            low: 0,
            range: u32::MAX,
            cache: 0,
            cache_size: 1,
            out: Vec::new(),
        }
    }

    fn shift_low(&mut self) {
        if self.low < 0xFF00_0000 || self.low > 0xFFFF_FFFF {
            let carry = (self.low >> 32) as u8;
            let mut temp = self.cache;
            loop {
                self.out.push(temp.wrapping_add(carry));
                temp = 0xFF;
                self.cache_size -= 1;
                if self.cache_size == 0 {
                    break;
                }
            }
            self.cache = ((self.low >> 24) & 0xFF) as u8;
        }
        self.cache_size += 1;
        self.low = u64::from((self.low as u32) << 8);
    }

    fn bit(&mut self, prob: &mut u16, bit: u32) {
        let bound = (self.range >> PROB_BITS) * u32::from(*prob);
        if bit == 0 {
            self.range = bound;
            *prob += ((PROB_ONE - u32::from(*prob)) >> MOVE_BITS) as u16;
        } else {
            self.low += u64::from(bound);
            self.range -= bound;
            *prob -= (u32::from(*prob) >> MOVE_BITS) as u16;
        }
        while self.range < RC_TOP {
            self.range <<= 8;
            self.shift_low();
        }
    }

    fn finish(mut self) -> Vec<u8> {
        for _ in 0..5 {
            self.shift_low();
        }
        self.out
    }
}

/// The matching decoder. It exists so that `bytes_entropy_coded` is the length
/// of something invertible.
struct RangeDecoder<'a> {
    code: u32,
    range: u32,
    input: &'a [u8],
    pos: usize,
}

impl<'a> RangeDecoder<'a> {
    fn new(input: &'a [u8]) -> Self {
        let mut d = Self {
            code: 0,
            range: u32::MAX,
            input,
            pos: 1,
        };
        for _ in 0..4 {
            let b = d.byte();
            d.code = (d.code << 8) | u32::from(b);
        }
        d
    }

    fn byte(&mut self) -> u8 {
        let b = self.input.get(self.pos).copied().unwrap_or(0);
        self.pos += 1;
        b
    }

    fn bit(&mut self, prob: &mut u16) -> u32 {
        let bound = (self.range >> PROB_BITS) * u32::from(*prob);
        let bit = if self.code < bound {
            self.range = bound;
            *prob += ((PROB_ONE - u32::from(*prob)) >> MOVE_BITS) as u16;
            0
        } else {
            self.code -= bound;
            self.range -= bound;
            *prob -= (u32::from(*prob) >> MOVE_BITS) as u16;
            1
        };
        while self.range < RC_TOP {
            self.range <<= 8;
            let b = self.byte();
            self.code = (self.code << 8) | u32::from(b);
        }
        bit
    }
}

/// Adaptive Elias-gamma over one integer stream.
struct IntModel {
    /// Unary bit-length, one context per position.
    nbits: [u16; NBITS + 1],
    /// Mantissa, one context per (length, position).
    payload: Vec<u16>,
}

impl IntModel {
    fn new() -> Self {
        Self {
            nbits: [PROB_INIT; NBITS + 1],
            payload: vec![PROB_INIT; (NBITS + 1) * NBITS],
        }
    }

    fn zigzag(v: i32) -> u64 {
        let v = i64::from(v);
        ((v << 1) ^ (v >> 63)) as u64
    }

    fn unzigzag(u: u64) -> i32 {
        (((u >> 1) as i64) ^ -((u & 1) as i64)) as i32
    }

    fn encode(&mut self, rc: &mut RangeEncoder, v: i32) {
        let u = Self::zigzag(v);
        let n = (64 - u.leading_zeros()) as usize;
        assert!(n < NBITS, "residual needs {n} bits, over the model's {NBITS}");
        for k in 0..n {
            rc.bit(&mut self.nbits[k], 1);
        }
        rc.bit(&mut self.nbits[n], 0);
        for j in (0..n.saturating_sub(1)).rev() {
            let bit = ((u >> j) & 1) as u32;
            rc.bit(&mut self.payload[n * NBITS + j], bit);
        }
    }

    fn decode(&mut self, rc: &mut RangeDecoder<'_>) -> i32 {
        let mut n = 0usize;
        while n < NBITS && rc.bit(&mut self.nbits[n]) == 1 {
            n += 1;
        }
        if n == 0 {
            return 0;
        }
        let mut u = 1u64;
        for j in (0..n - 1).rev() {
            let bit = u64::from(rc.bit(&mut self.payload[n * NBITS + j]));
            u = (u << 1) | bit;
            let _ = j;
        }
        Self::unzigzag(u)
    }
}

/// Range-code the residual stream, decode it back, and assert equality.
///
/// Returns the encoded length in bytes. The four-byte count prefix is counted:
/// a length nobody can decode without knowing the count is not a length.
fn entropy_code(res: &[[i32; 4]]) -> (usize, bool) {
    let mut rc = RangeEncoder::new();
    let mut models = [
        IntModel::new(),
        IntModel::new(),
        IntModel::new(),
        IntModel::new(),
    ];
    let mut op = PROB_INIT;
    for row in res {
        rc.bit(&mut op, 1);
        for a in 0..4 {
            models[a].encode(&mut rc, row[a]);
        }
    }
    let body = rc.finish();

    let mut models = [
        IntModel::new(),
        IntModel::new(),
        IntModel::new(),
        IntModel::new(),
    ];
    let mut op = PROB_INIT;
    let mut dec = RangeDecoder::new(&body);
    let mut ok = true;
    for row in res {
        if dec.bit(&mut op) != 1 {
            ok = false;
        }
        for a in 0..4 {
            if models[a].decode(&mut dec) != row[a] {
                ok = false;
            }
        }
    }
    (body.len() + 4, ok)
}

/// First-order variant, for the mechanism column: how much of the win is the
/// second-order predictor rather than the coder.
fn entropy_code_first_order(centres: &[[f64; 3]], radii: &[f64]) -> usize {
    let mut res = Vec::with_capacity(centres.len());
    let mut prev = [0i32; 3];
    let mut prev_r = 0i32;
    for (c, r) in centres.iter().zip(radii) {
        let q = [quantise(c[0]), quantise(c[1]), quantise(c[2])];
        let qr = quantise(*r);
        res.push([
            q[0].wrapping_sub(prev[0]),
            q[1].wrapping_sub(prev[1]),
            q[2].wrapping_sub(prev[2]),
            qr.wrapping_sub(prev_r),
        ]);
        prev = q;
        prev_r = qr;
    }
    entropy_code(&res).0
}

fn zstd_len(bytes: &[u8], level: i32) -> usize {
    zstd::bulk::compress(bytes, level)
        .expect("zstd compresses an in-memory buffer")
        .len()
}

// ── the coaxial arm ─────────────────────────────────────────────────────────

/// A base field with a list of capsules subtracted, folded through the crate's
/// own [`apply`].
struct Coaxial<'a> {
    caps: &'a [Capsule<f64>],
}

impl Sdf for Coaxial<'_> {
    type Scalar = f64;

    fn sample(&self, p: [f64; 3]) -> f64 {
        let mut v = ground(p);
        for c in self.caps {
            v = apply(BrushOp::Subtract, v, c.sample(p));
        }
        v
    }
}

struct Grid {
    shape: RuntimeShape3,
    origin: [f64; 3],
    cell: f64,
}

fn hash_of(mc: &mut MarchingCubes<f64>, out: &mut MeshBuffer<f64>, caps: &[Capsule<f64>], g: &Grid) -> u64 {
    out.reset();
    mc.extract(&Coaxial { caps }, &g.shape, g.origin, g.cell, out)
        .expect("tunnel extraction");
    mesh_hash(out)
}

/// What the coaxial arm produced.
struct CoaxialResult {
    /// End-to-end axial length of the bore, in world units. The aim points it
    /// came from are all on one line, captured from the same closed loop.
    axial_length: f64,
    radius: f64,
    delta: f64,
    overlapping_pairs: u64,
    /// ✗41's leave-one-out count: capsules whose individual removal moves
    /// `mesh_hash`.
    necessary: usize,
    /// Survivors of the greedy bit-exact drop.
    surviving: usize,
    /// The greedy result re-verified over the whole tunnel in one go.
    hash_verified: bool,
    /// Capsules after merging maximal collinear contiguous equal-radius runs.
    merged: usize,
    /// Largest `|fold_merged − fold_full|` over the tunnel grid, in cells.
    merge_deviation_cells: f64,
    merge_hash_equal: bool,
    triangles: usize,
    grid_samples: usize,
    /// The wobbled bore: the same closed loop with a hand that drifts, which is
    /// what a real stroke is. `(k, runs, deviation in cells)` for each `k`
    /// capsules folded into one hull capsule.
    wobble: Vec<(usize, usize, f64)>,
    /// Yaw drift of the wobbled bore, in radians per second.
    wobble_yaw_rate: f64,
    /// Largest transverse departure of the wobbled bore's aim points from the
    /// straight line through its ends, in cells.
    wobble_bow_cells: f64,
}

/// Capture a dead-straight bore from the same closed loop: eye buried at
/// `y = −1`, look fixed along `+x`, button held, fly forward.
///
/// The aim points are `eye + [1,0,0]·t` with `eye.y` and `eye.z` never touched,
/// so they are exactly collinear — which is what "along one axis" means, and it
/// is a property of the loop rather than of a hardcoded list.
fn capture_bore(points_wanted: usize, yaw_rate: f64) -> (Vec<[f64; 3]>, f64) {
    let (lo, hi) = sandbox(1024);
    let mut log = Log::new(lo, hi);
    let radius = RADIUS_DEFAULT;
    let mut yaw = 0.0;
    let mut eye = [-50.0, -1.0, 0.0];
    let mut clock = 0.0;
    let mut stroke_last: Option<[f64; 3]> = None;
    let mut points = Vec::with_capacity(points_wanted);
    let frames = 200_000u64;
    for _ in 0..frames {
        if points.len() >= points_wanted {
            break;
        }
        yaw += yaw_rate * FRAME_DT;
        // The look direction is `+x` rotated by `yaw` about `y`, so at
        // `yaw_rate = 0` it is exactly `[1, 0, 0]` and `eye[1]`, `eye[2]` are
        // never touched.
        let dir = [yaw.cos(), 0.0, yaw.sin()];
        // 3 u/s, which is a quarter of a unit per throttle period. The bore is
        // self-regulating at that speed: the eye trails the face by `AIM_NEAR`,
        // the trace is a handful of steps, and the spacing settles on the
        // brush radius. Slower and the face outruns `AIM_STEPS` -- which is
        // game_dig's own reason a player has to walk into their own tunnel.
        eye = add(eye, scale(dir, 3.0 * FRAME_DT));
        clock += FRAME_DT;
        if clock < EDIT_PERIOD {
            continue;
        }
        let value = |p: [f64; 3]| log.march_value(p);
        let Some(t) = trace(&value, eye, dir, lo, hi) else {
            continue;
        };
        let centre = add(eye, scale(dir, t));
        if let Some(last) = stroke_last
            && length(sub(centre, last)) < radius * 0.5
        {
            continue;
        }
        log.push(centre, radius);
        stroke_last = Some(centre);
        points.push(centre);
        clock = 0.0;
    }
    (points, radius)
}

/// `k` capsules folded into one hull capsule, and how far that moves the fold.
///
/// The deviation is measured on the **field**, over every cell corner of the
/// run's padded box, in cells — not on a mesh hash. That is deliberate: a merge
/// is a change of representation rather than a removal, so the question is how
/// far the surface moved, and `mesh_hash` can only answer *whether* it moved.
/// P-96's logic in a different costume: sub-cell is the bar a player can see.
fn merge_deviation(caps: &[Capsule<f64>], k: usize) -> (usize, f64) {
    let mut runs = 0usize;
    let mut worst = 0.0f64;
    let mut i = 0usize;
    while i < caps.len() {
        let j = (i + k).min(caps.len());
        let run = &caps[i..j];
        runs += 1;
        let hull = Capsule {
            a: run[0].a,
            b: run[run.len() - 1].b,
            radius: run[0].radius,
        };
        let mut blo = [f64::INFINITY; 3];
        let mut bhi = [f64::NEG_INFINITY; 3];
        for c in run {
            for a in 0..3 {
                blo[a] = blo[a].min(c.a[a]).min(c.b[a]) - c.radius;
                bhi[a] = bhi[a].max(c.a[a]).max(c.b[a]) + c.radius;
            }
        }
        let steps = [
            (((bhi[0] - blo[0]) / CELL_SIZE).ceil() as u32) + 1,
            (((bhi[1] - blo[1]) / CELL_SIZE).ceil() as u32) + 1,
            (((bhi[2] - blo[2]) / CELL_SIZE).ceil() as u32) + 1,
        ];
        for z in 0..steps[2] {
            for y in 0..steps[1] {
                for x in 0..steps[0] {
                    let p = [
                        blo[0] + f64::from(x) * CELL_SIZE,
                        blo[1] + f64::from(y) * CELL_SIZE,
                        blo[2] + f64::from(z) * CELL_SIZE,
                    ];
                    let mut a = f64::INFINITY;
                    for c in run {
                        a = a.min(c.sample(p));
                    }
                    let d = (a - hull.sample(p)).abs();
                    if d > worst {
                        worst = d;
                    }
                }
            }
        }
        i = j;
    }
    (runs, worst / CELL_SIZE)
}

fn coaxial_arm() -> CoaxialResult {
    let (points, radius) = capture_bore(COAXIAL_BRUSHES + 1, 0.0);
    assert!(
        points.len() == COAXIAL_BRUSHES + 1,
        "the bore capture produced {} aim points, not {}",
        points.len(),
        COAXIAL_BRUSHES + 1
    );
    let caps: Vec<Capsule<f64>> = (0..COAXIAL_BRUSHES)
        .map(|i| Capsule {
            a: points[i],
            b: points[i + 1],
            radius,
        })
        .collect();

    // The registered vacuity control's second half. Exact segment-to-segment
    // distance, all 19,900 pairs.
    let mut overlapping_pairs = 0u64;
    for i in 0..caps.len() {
        for j in (i + 1)..caps.len() {
            let d = seg_seg_distance(caps[i].a, caps[i].b, caps[j].a, caps[j].b);
            if d < caps[i].radius + caps[j].radius {
                overlapping_pairs += 1;
            }
        }
    }

    let mut delta_sum = 0.0;
    for i in 0..COAXIAL_BRUSHES {
        delta_sum += length(sub(points[i + 1], points[i]));
    }
    let delta = delta_sum / COAXIAL_BRUSHES as f64;

    // The tunnel's grid: the whole bore, padded so neither end is clipped and
    // the walls are interior.
    let mut glo = [f64::INFINITY; 3];
    let mut ghi = [f64::NEG_INFINITY; 3];
    for p in &points {
        for a in 0..3 {
            glo[a] = glo[a].min(p[a]);
            ghi[a] = ghi[a].max(p[a]);
        }
    }
    let pad = radius + 4.0 * CELL_SIZE;
    for a in 0..3 {
        glo[a] -= pad;
        ghi[a] += pad;
    }
    let dims = [
        (((ghi[0] - glo[0]) / CELL_SIZE).ceil() as u32) + 1,
        (((ghi[1] - glo[1]) / CELL_SIZE).ceil() as u32) + 1,
        (((ghi[2] - glo[2]) / CELL_SIZE).ceil() as u32) + 1,
    ];
    let grid = Grid {
        shape: RuntimeShape3::new(dims).expect("tunnel grid fits u32"),
        origin: glo,
        cell: CELL_SIZE,
    };
    let grid_samples = (dims[0] as usize) * (dims[1] as usize) * (dims[2] as usize);

    let mut mc = MarchingCubes::<f64>::default();
    let mut mesh = MeshBuffer::<f64>::default();
    let full = hash_of(&mut mc, &mut mesh, &caps, &grid);
    let triangles = mesh.indices.len() / 3;
    assert!(
        triangles > 0,
        "the tunnel meshed to nothing, so no removal could have changed anything"
    );

    // ✗41's leave-one-out necessity.
    let mut necessary = 0usize;
    let mut scratch: Vec<Capsule<f64>> = Vec::with_capacity(caps.len());
    for i in 0..caps.len() {
        scratch.clear();
        scratch.extend(caps.iter().enumerate().filter(|(j, _)| *j != i).map(|(_, c)| *c));
        if hash_of(&mut mc, &mut mesh, &scratch, &grid) != full {
            necessary += 1;
        }
    }

    // The compressor: greedy bit-exact drop, in order.
    let mut kept: Vec<bool> = vec![true; caps.len()];
    for i in 0..caps.len() {
        kept[i] = false;
        scratch.clear();
        scratch.extend(caps.iter().enumerate().filter(|(j, _)| kept[*j]).map(|(_, c)| *c));
        if hash_of(&mut mc, &mut mesh, &scratch, &grid) != full {
            kept[i] = true;
        }
    }
    let surviving = kept.iter().filter(|k| **k).count();
    scratch.clear();
    scratch.extend(caps.iter().enumerate().filter(|(j, _)| kept[*j]).map(|(_, c)| *c));
    let hash_verified = hash_of(&mut mc, &mut mesh, &scratch, &grid) == full;

    // The other compressor: merge maximal collinear contiguous equal-radius
    // runs into their hull. Exact in ℝ — `min_i d(p, seg_i) = d(p, seg_hull)`
    // when the segments are collinear and contiguous — so this is a
    // representation change rather than a redundancy elimination.
    let mut merged_caps: Vec<Capsule<f64>> = Vec::new();
    let mut i = 0usize;
    while i < caps.len() {
        let mut j = i;
        while j + 1 < caps.len()
            && collinear_contiguous(&caps[j], &caps[j + 1])
        {
            j += 1;
        }
        merged_caps.push(Capsule {
            a: caps[i].a,
            b: caps[j].b,
            radius: caps[i].radius,
        });
        i = j + 1;
    }
    let merged = merged_caps.len();
    let merge_hash_equal = hash_of(&mut mc, &mut mesh, &merged_caps, &grid) == full;
    // How far apart the two folds actually are, over every grid sample.
    let a = Coaxial { caps: &caps };
    let b = Coaxial { caps: &merged_caps };
    let mut worst = 0.0f64;
    for z in 0..dims[2] {
        for y in 0..dims[1] {
            for x in 0..dims[0] {
                let p = [
                    glo[0] + f64::from(x) * CELL_SIZE,
                    glo[1] + f64::from(y) * CELL_SIZE,
                    glo[2] + f64::from(z) * CELL_SIZE,
                ];
                let d = (a.sample(p) - b.sample(p)).abs();
                if d > worst {
                    worst = d;
                }
            }
        }
    }

    // The wobbled bore: the same closed loop with a hand that drifts at 0.05
    // rad/s, which over a 50-unit bore bows the tunnel by several cells.
    // Exact collinearity is measure-zero in a real stroke, so this is where the
    // merge has to earn its keep with a tolerance rather than with an identity.
    const WOBBLE_YAW_RATE: f64 = 0.01;
    let (wpoints, wradius) = capture_bore(COAXIAL_BRUSHES + 1, WOBBLE_YAW_RATE);
    assert!(
        wpoints.len() == COAXIAL_BRUSHES + 1,
        "the wobbled bore capture produced {} aim points, not {}",
        wpoints.len(),
        COAXIAL_BRUSHES + 1
    );
    let wcaps: Vec<Capsule<f64>> = (0..COAXIAL_BRUSHES)
        .map(|i| Capsule {
            a: wpoints[i],
            b: wpoints[i + 1],
            radius: wradius,
        })
        .collect();
    // How bent the wobbled bore actually is: the largest distance from an aim
    // point to the straight line through the two ends. Zero would mean the
    // wobble arm is the straight arm again and its tolerance sweep says nothing.
    let chord = sub(wpoints[COAXIAL_BRUSHES], wpoints[0]);
    let chord_len = length(chord);
    let unit = scale(chord, 1.0 / chord_len);
    let mut bow = 0.0f64;
    for p in &wpoints {
        let d = sub(*p, wpoints[0]);
        let along = dot(d, unit);
        let off = length(sub(d, scale(unit, along)));
        if off > bow {
            bow = off;
        }
    }
    assert!(
        bow > CELL_SIZE,
        "the wobbled bore bows by {bow} world units, under one cell, so it is the straight \
         bore under another name and its tolerance sweep would be vacuous"
    );
    let mut wobble = Vec::new();
    for k in [1usize, 2, 5, 10, 20, 50, 200] {
        let (runs, dev) = merge_deviation(&wcaps, k);
        wobble.push((k, runs, dev));
    }

    CoaxialResult {
        axial_length: length(sub(points[COAXIAL_BRUSHES], points[0])),
        radius,
        delta,
        overlapping_pairs,
        necessary,
        surviving,
        hash_verified,
        merged,
        merge_deviation_cells: worst / CELL_SIZE,
        merge_hash_equal,
        triangles,
        grid_samples,
        wobble,
        wobble_yaw_rate: WOBBLE_YAW_RATE,
        wobble_bow_cells: bow / CELL_SIZE,
    }
}

/// Whether two capsules are collinear, contiguous and equal-radius, so their
/// union is exactly the hull capsule.
fn collinear_contiguous(p: &Capsule<f64>, q: &Capsule<f64>) -> bool {
    const TOL: f64 = 1e-12;
    if (p.radius - q.radius).abs() > TOL {
        return false;
    }
    if length(sub(p.b, q.a)) > TOL {
        return false;
    }
    let u = sub(p.b, p.a);
    let v = sub(q.b, q.a);
    let lu = length(u);
    let lv = length(v);
    if lu <= TOL || lv <= TOL {
        return false;
    }
    // Same direction, to within a tolerance on the unit vectors.
    let cu = scale(u, 1.0 / lu);
    let cv = scale(v, 1.0 / lv);
    length(sub(cu, cv)) <= 1e-9
}

/// The lens a sound drop of `k` consecutive coaxial capsules would leave, in
/// world units: `r − √(max(0, r² − (kδ/2)²))`. This is the arithmetic the SHARE
/// paragraph runs and it is on the row so a reader can check it.
fn lens_thickness(radius: f64, delta: f64, k: f64) -> f64 {
    let half = k * delta * 0.5;
    radius - (radius * radius - half * half).max(0.0).sqrt()
}

/// Chunks a brush's padded box overlaps, at this granularity. The knob that has
/// to move something else in the same harness for C1's constancy to be a
/// measurement rather than a definition.
fn dirty_chunks(layout: &ChunkLayout<f64>, centre: [f64; 3], radius: f64) -> u64 {
    let reach = radius + layout.cell_size();
    let a = layout.chunk_of([centre[0] - reach, centre[1] - reach, centre[2] - reach]);
    let b = layout.chunk_of([centre[0] + reach, centre[1] + reach, centre[2] + reach]);
    let mut n = 1u64;
    for axis in 0..3 {
        n *= (b.coords[axis] - a.coords[axis] + 1) as u64;
    }
    n
}

// ── the run ─────────────────────────────────────────────────────────────────

/// Everything one CSV row needs from the compressors.
struct Sizes {
    uncompressed: usize,
    zstd19: usize,
    zstd3: usize,
    zstd_deltas: usize,
    entropy: usize,
    entropy_first: usize,
    roundtrip_ok: bool,
}

fn measure_sizes(centres: &[[f64; 3]], radii: &[f64]) -> Sizes {
    let raw = serialise(centres, radii);
    let res = residuals(centres, radii);
    let (entropy, roundtrip_ok) = entropy_code(&res);
    Sizes {
        uncompressed: raw.len(),
        zstd19: zstd_len(&raw, 19),
        zstd3: zstd_len(&raw, 3),
        zstd_deltas: zstd_len(&serialise_residuals(centres, radii), 19),
        entropy,
        entropy_first: entropy_code_first_order(centres, radii),
        roundtrip_ok,
    }
}

/// Worlds swept. `128` is game_dig's own (8 chunks of 16 cells).
const WORLDS: [u32; 5] = [64, 128, 256, 512, 1024];
/// Granularities swept, over the 128-cell world's one trace.
const GRANULARITIES: [u32; 5] = [2, 4, 8, 32, 64];

fn main() {
    if !std::env::args().any(|a| a == "--bench") {
        return;
    }

    common::experiment::run(isomesh::experiment!("P-94"), |run| {
        println!("SHARE, recomputed before the run:");
        println!(
            "  throttle at 60 Hz quantises to 5 frames = {:.4} s, so an hour offers {} dig \
             attempts, not {}",
            5.0 * FRAME_DT,
            (SESSION_SECONDS / (5.0 * FRAME_DT)) as u64,
            (SESSION_SECONDS * 12.5) as u64
        );
        println!(
            "  record is {RECORD_BYTES} B serialised, {} B in memory \
             (size_of::<Brush<Sphere<f32>>>())",
            core::mem::size_of::<Brush<Sphere<f32>>>()
        );
        println!(
            "  C1 reachable: 43200 x {RECORD_BYTES} B = {:.0} B against a {:.0} B bar; at the \
             08-11 doc's 48 B it would be {:.0} B and FAIL",
            43_200.0 * RECORD_BYTES as f64,
            C1_BAR_BYTES,
            43_200.0 * 48.0
        );
        println!(
            "  C2 UNREACHABLE: 200 -> 40 drops runs of k=4; the lens a sound drop leaves is \
             r - sqrt(r^2 - (k.delta/2)^2), which is the full radius as soon as k.delta >= 2r, \
             and the distance gate guarantees delta >= r/2"
        );
        println!();

        let coax = coaxial_arm();
        println!(
            "coaxial arm: {} capsules, r = {:.6}, mean delta = {:.6} ({:.4} r), \
             {} overlapping pairs, {} triangles over {} samples",
            COAXIAL_BRUSHES,
            coax.radius,
            coax.delta,
            coax.delta / coax.radius,
            coax.overlapping_pairs,
            coax.triangles,
            coax.grid_samples
        );
        println!(
            "  leave-one-out necessary {}, greedy survivors {} (verified {}), \
             merged runs {} (bit-exact {}, deviation {:.3e} cells)",
            coax.necessary,
            coax.surviving,
            coax.hash_verified,
            coax.merged,
            coax.merge_hash_equal,
            coax.merge_deviation_cells
        );
        println!(
            "  lens if one capsule is dropped: {:.6} world units = {:.4} cells; if four: \
             {:.6} = {:.4} cells",
            lens_thickness(coax.radius, coax.delta, 1.0),
            lens_thickness(coax.radius, coax.delta, 1.0) / CELL_SIZE,
            lens_thickness(coax.radius, coax.delta, 4.0),
            lens_thickness(coax.radius, coax.delta, 4.0) / CELL_SIZE
        );
        println!(
            "  wobbled bore: yaw {:.3} rad/s, bow {:.3} cells; k capsules -> one hull:",
            coax.wobble_yaw_rate, coax.wobble_bow_cells
        );
        for (k, runs, dev) in &coax.wobble {
            println!("    k = {k:>3}: {runs:>3} runs, field deviation {dev:.6} cells");
        }

        // The registered vacuity control, second half.
        assert!(
            coax.overlapping_pairs > 0,
            "the coaxial arm's capsules are disjoint, so C2 would be compressing nothing"
        );
        assert!(
            coax.hash_verified,
            "the greedy drop's survivors do not reproduce the 200-capsule mesh bit-exactly, \
             so `surviving_brushes` is not a sound compression"
        );

        let collapse_ratio = COAXIAL_BRUSHES as f64 / coax.surviving as f64;
        let c2_holds = coax.surviving < 40;

        let mut sessions: Vec<SessionResult> = Vec::new();
        for w in WORLDS {
            let s = run_session(w);
            println!(
                "session world {:>4} cells: {} strokes, {} edits, {} closed-loop, \
                 refusals {}/{}, mean step {:.4}, span {:.1} cells, nested {}, radii {}",
                w,
                s.strokes,
                s.centres.len(),
                s.closed_loop,
                s.refusals_distance,
                s.refusals_nohit,
                s.mean_step,
                s.span_cells,
                s.nested_dropped,
                s.distinct_radii
            );
            println!(
                "    by phase (bore/sweep/chamber/relocate): edits {:?}, no-hit {:?}",
                s.edits_by_phase, s.nohit_by_phase
            );
            sessions.push(s);
        }

        // The registered vacuity control, first half.
        for s in &sessions {
            assert!(
                s.closed_loop > 0,
                "world {}: no edit's aim point got past the unedited surface, so this trace \
                 has no feedback from the field and is a sweep rather than a session",
                s.world_cells
            );
            assert!(
                !s.centres.is_empty(),
                "world {}: the session logged nothing",
                s.world_cells
            );
        }

        let sizes: Vec<Sizes> = sessions
            .iter()
            .map(|s| measure_sizes(&s.centres, &s.radii))
            .collect();
        for z in &sizes {
            assert!(
                z.roundtrip_ok,
                "the range coder's output does not decode back to its input, so \
                 `bytes_entropy_coded` is not the length of an encoding"
            );
        }

        // The sweep control, at the 128-cell session's own edit count.
        let idx128 = WORLDS.iter().position(|w| *w == 128).expect("128 is swept");
        let (sc_centres, sc_radii) = sweep_control(sessions[idx128].centres.len());
        let sc_sizes = measure_sizes(&sc_centres, &sc_radii);
        println!(
            "sweep control ({} edits, no field feedback): uncompressed {} B, zstd19 {} B, \
             entropy {} B, advantage {:.4}",
            sc_centres.len(),
            sc_sizes.uncompressed,
            sc_sizes.zstd19,
            sc_sizes.entropy,
            sc_sizes.zstd19 as f64 / sc_sizes.entropy as f64
        );

        // Per-edit widths, and the cross-world comparison C1's second half asks
        // for. `bytes_per_edit_at_other_world_size` is the value from the
        // 1024-cell world on every row, and from the 64-cell world on the
        // 1024-cell row.
        let bpe: Vec<f64> = sizes
            .iter()
            .zip(&sessions)
            .map(|(z, s)| z.uncompressed as f64 / s.centres.len() as f64)
            .collect();
        let bpe_small = bpe[0];
        let bpe_large = bpe[bpe.len() - 1];

        let mut rows: Vec<(usize, u32)> = WORLDS.iter().map(|_| (0, CHUNK_CELLS)).collect();
        for (i, r) in rows.iter_mut().enumerate() {
            r.0 = i;
        }
        for c in GRANULARITIES {
            rows.push((idx128, c));
        }

        let mut c1_size_all = true;
        let mut c1_const_all = true;
        for (i, _) in &rows {
            if sizes[*i].uncompressed as f64 >= C1_BAR_BYTES {
                c1_size_all = false;
            }
            let other = if WORLDS[*i] == 1024 { bpe_small } else { bpe_large };
            if (bpe[*i] - other).abs() > f64::EPSILON {
                c1_const_all = false;
            }
        }
        let c1_holds = c1_size_all && c1_const_all;

        for (i, chunk_cells) in rows {
            let s = &sessions[i];
            let z = &sizes[i];
            let edits = s.centres.len();
            let layout = ChunkLayout::<f64>::new(chunk_cells, CELL_SIZE, sandbox(s.world_cells).0)
                .expect("chunk layout");
            let dirty: u64 = s
                .centres
                .iter()
                .zip(&s.radii)
                .map(|(c, r)| dirty_chunks(&layout, *c, *r))
                .sum();
            let other = if s.world_cells == 1024 {
                bpe_small
            } else {
                bpe_large
            };
            let advantage = z.zstd19 as f64 / z.entropy as f64;
            let c3_holds = advantage >= 2.0;
            let (lo, hi) = sandbox(s.world_cells);

            run.record(&[
                ("session_minutes", format!("{:.1}", SESSION_SECONDS / 60.0)),
                ("strokes", s.strokes.to_string()),
                ("edits", edits.to_string()),
                ("bytes_uncompressed", z.uncompressed.to_string()),
                ("bytes_per_edit", format!("{:.6}", bpe[i])),
                ("world_cells", s.world_cells.to_string()),
                ("chunk_cells", chunk_cells.to_string()),
                (
                    "bytes_per_edit_at_other_world_size",
                    format!("{other:.6}"),
                ),
                ("coaxial_brushes", COAXIAL_BRUSHES.to_string()),
                ("surviving_brushes", coax.surviving.to_string()),
                ("collapse_ratio", format!("{collapse_ratio:.6}")),
                ("overlapping_pairs", coax.overlapping_pairs.to_string()),
                ("bytes_zstd", z.zstd19.to_string()),
                ("bytes_entropy_coded", z.entropy.to_string()),
                ("compression_advantage", format!("{advantage:.6}")),
                ("c1_holds", c1_holds.to_string()),
                ("c2_holds", c2_holds.to_string()),
                ("c3_holds", c3_holds.to_string()),
                // ── extras ──
                ("cell_size", format!("{CELL_SIZE:.6}")),
                ("sandbox_extent", format!("{:.3}", hi[0] - lo[0])),
                ("strokes_nominal_12_5", "45000".to_string()),
                (
                    "bytes_at_43200_edits",
                    format!("{:.0}", 43_200.0 * RECORD_BYTES as f64),
                ),
                (
                    "bytes_per_edit_in_memory",
                    core::mem::size_of::<Brush<Sphere<f32>>>().to_string(),
                ),
                (
                    "hour_bytes_over_estimate",
                    format!("{:.6}", z.uncompressed as f64 / ESTIMATE_BYTES),
                ),
                (
                    "hour_bytes_over_teardown",
                    format!("{:.6}", z.uncompressed as f64 / TEARDOWN_BYTES),
                ),
                (
                    "entropy_bytes_over_teardown",
                    format!("{:.6}", z.entropy as f64 / TEARDOWN_BYTES),
                ),
                (
                    "dirty_chunks_per_edit",
                    format!("{:.6}", dirty as f64 / edits as f64),
                ),
                ("trace_span_cells", format!("{:.3}", s.span_cells)),
                ("quantum", format!("{QUANTUM:.9}")),
                (
                    "position_bits_needed",
                    format!(
                        "{}",
                        ((hi[0] - lo[0]) / QUANTUM).log2().ceil() as u32
                    ),
                ),
                (
                    "entropy_bytes_per_edit",
                    format!("{:.6}", z.entropy as f64 / edits as f64),
                ),
                ("bytes_zstd_level3", z.zstd3.to_string()),
                ("bytes_zstd_of_deltas", z.zstd_deltas.to_string()),
                (
                    "zstd_deltas_over_entropy",
                    format!("{:.6}", z.zstd_deltas as f64 / z.entropy as f64),
                ),
                ("bytes_entropy_first_order", z.entropy_first.to_string()),
                ("entropy_roundtrip_ok", z.roundtrip_ok.to_string()),
                ("closed_loop_edits", s.closed_loop.to_string()),
                ("gate_refusals_distance", s.refusals_distance.to_string()),
                ("gate_refusals_nohit", s.refusals_nohit.to_string()),
                ("edits_bore", s.edits_by_phase[0].to_string()),
                ("edits_sweep", s.edits_by_phase[1].to_string()),
                ("edits_chamber", s.edits_by_phase[2].to_string()),
                ("edits_relocate", s.edits_by_phase[3].to_string()),
                ("nohit_bore", s.nohit_by_phase[0].to_string()),
                ("nohit_sweep", s.nohit_by_phase[1].to_string()),
                ("nohit_chamber", s.nohit_by_phase[2].to_string()),
                ("nohit_relocate", s.nohit_by_phase[3].to_string()),
                ("mean_aim_step", format!("{:.6}", s.mean_step)),
                (
                    "mean_aim_step_over_radius",
                    format!("{:.6}", s.mean_step / RADIUS_DEFAULT),
                ),
                ("distinct_radii", s.distinct_radii.to_string()),
                ("session_nested_dropped", s.nested_dropped.to_string()),
                (
                    "session_nested_fraction",
                    format!("{:.6}", s.nested_dropped as f64 / edits as f64),
                ),
                ("c1_size_holds", c1_size_all.to_string()),
                ("c1_constancy_holds", c1_const_all.to_string()),
                ("coaxial_radius", format!("{:.6}", coax.radius)),
                ("coaxial_delta", format!("{:.6}", coax.delta)),
                (
                    "coaxial_delta_over_radius",
                    format!("{:.6}", coax.delta / coax.radius),
                ),
                (
                    "coaxial_necessary_leave_one_out",
                    coax.necessary.to_string(),
                ),
                ("coaxial_hash_verified", coax.hash_verified.to_string()),
                ("coaxial_triangles", coax.triangles.to_string()),
                ("coaxial_grid_samples", coax.grid_samples.to_string()),
                (
                    "coaxial_axial_length",
                    format!("{:.6}", coax.axial_length),
                ),
                ("coaxial_merged_brushes", coax.merged.to_string()),
                (
                    "coaxial_merge_ratio",
                    format!("{:.6}", COAXIAL_BRUSHES as f64 / coax.merged as f64),
                ),
                (
                    "coaxial_merge_hash_equal",
                    coax.merge_hash_equal.to_string(),
                ),
                (
                    "coaxial_merge_deviation_cells",
                    format!("{:.3e}", coax.merge_deviation_cells),
                ),
                (
                    "lens_cells_drop_1",
                    format!(
                        "{:.6}",
                        lens_thickness(coax.radius, coax.delta, 1.0) / CELL_SIZE
                    ),
                ),
                (
                    "lens_cells_drop_4",
                    format!(
                        "{:.6}",
                        lens_thickness(coax.radius, coax.delta, 4.0) / CELL_SIZE
                    ),
                ),
                (
                    "wobble_yaw_rate",
                    format!("{:.6}", coax.wobble_yaw_rate),
                ),
                (
                    "wobble_bow_cells",
                    format!("{:.6}", coax.wobble_bow_cells),
                ),
                ("wobble_runs_k2", coax.wobble[1].1.to_string()),
                ("wobble_dev_cells_k1", format!("{:.6}", coax.wobble[0].2)),
                ("wobble_dev_cells_k2", format!("{:.6}", coax.wobble[1].2)),
                ("wobble_dev_cells_k5", format!("{:.6}", coax.wobble[2].2)),
                ("wobble_dev_cells_k10", format!("{:.6}", coax.wobble[3].2)),
                ("wobble_dev_cells_k20", format!("{:.6}", coax.wobble[4].2)),
                ("wobble_dev_cells_k50", format!("{:.6}", coax.wobble[5].2)),
                ("wobble_dev_cells_k200", format!("{:.6}", coax.wobble[6].2)),
                ("sweep_control_edits", sc_centres.len().to_string()),
                ("sweep_control_zstd", sc_sizes.zstd19.to_string()),
                ("sweep_control_entropy", sc_sizes.entropy.to_string()),
                (
                    "sweep_control_advantage",
                    format!(
                        "{:.6}",
                        sc_sizes.zstd19 as f64 / sc_sizes.entropy as f64
                    ),
                ),
                ("sweep_control_closed_loop", "0".to_string()),
            ]);
        }
    });
}
