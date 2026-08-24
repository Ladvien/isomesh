//! E-315 — the tape you keep is twenty times too big.
//!
//! ```bash
//! cd bevy_isomesh && cargo run --example game_edit_tape_trim --release
//! ```
//!
//! **Always `--release`.** Startup runs the whole leave-one-out ablation — 1,571
//! chunk re-meshes — before the window opens. In release that is a couple of
//! seconds across every core; a debug build meshes 20-50x slower and you would
//! wait several minutes for a window.
//!
//! It runs itself and tours five chunks, showing each one's tape at three
//! lengths. `[` and `]` step chunk by chunk, `T` cycles the tape length, `B`
//! hides the brush gizmos, `K` the chunk boxes, `X` hands control back to the
//! tour. `Space` freezes the tour. `ISOMESH_FIELD=1|2|3|4|5` pins one of the five
//! stops for a still without a keyboard.
//!
//! ```bash
//! ISOMESH_CAPTURE_FRAMES=150 ISOMESH_CAPTURE_EVERY=2 FPS=10 \
//!   ./scripts/record_gif.sh game_edit_tape_trim \
//!   docs/gifs/the-tape-you-keep-is-twenty-times-too-big.gif
//! ```
//!
//! **150 frames is a floor, not a preference.** The tour is five chunks × three
//! tape lengths, and the capture window is `SETTLE + FRAMES × EVERY` ticks, so
//! this gives each of the fifteen states ten captured frames — one second at
//! `FPS=10`, which is about how long it takes to see that three hashes are the
//! same number. Half of it and the clip is a flicker book. Measured at **0.92 MB
//! over 149 frames**, comfortably inside the 0.7-4.8 MB the committed clips sit
//! within, because the camera is a cut rather than a move and a GIF pays for
//! motion.
//!
//! Demonstrates **M-358 / ✗41 (P-59, R-057)** (`docs/experiments/p-59.csv`) as a
//! studio meets it: not as an ablation sweep, but as a destructible world in
//! which the edit history you are paying to fold is almost entirely dead weight,
//! chunk by chunk, with the mesh not moving by one bit while you delete it.
//!
//! # The finding, in one paragraph
//!
//! A destructible world is a **tape** of brush edits folded over a base field,
//! and every sample of every chunk walks the whole tape. P-39 established that a
//! Lipschitz interval bound over a chunk deletes the brushes that provably
//! cannot win there: 64 brushes prune to a median of **19** survivors per chunk
//! with the mesh byte-identical on 64 of 64 chunks (M-341). P-59 asked how many
//! of those 19 are actually *needed*, by removing each surviving brush on its own
//! and re-meshing. The answer is a factor of twenty:
//!
//! - Removing every non-survivor at once changes no chunk's `mesh_hash` — **64 of
//!   64**, the soundness control, and it is what makes every other number here
//!   readable.
//! - **1,434 of 1,507** surviving brushes across the world (95.16%) can be
//!   dropped individually with the mesh bit-identical.
//! - Re-meshing each chunk from its `necessary` brushes **alone** — every
//!   individually-unnecessary survivor removed at once, order preserved — is
//!   bit-identical on **64 of 64** chunks. So the world's tape cuts from **1,507
//!   survivors to 73 brushes, bit-exactly**: a further **20.6x** on top of P-39's
//!   64 → 19.
//!
//! On this fixture 51 of 64 chunks need **zero** brushes: the base sphere alone
//! meshes to the same bytes as the sphere carved by 64 edits. Chunk `2-1-0` is
//! the stop this tour lingers on for that — 3,227 triangles, 35 survivors, and a
//! `necessary` tape of length zero.
//!
//! # The other half, which is not a footnote
//!
//! Deciding necessity **cost 1,571 re-meshes**. That is what this demo spends its
//! startup on, and it is on the HUD beside every win, because a measurement of
//! headroom is not a pruner. Nothing here is shippable as a runtime test; what is
//! shippable is the knowledge that a tighter bound has 20x of room, which closes
//! no direction and opens one.
//!
//! # C3 failed, and the failure is the interesting part
//!
//! The registration named one cause for the over-keep: a surviving brush is
//! unnecessary because its interval over the chunk is *far from the surface*,
//! more than one cell size clear of zero. That is real and dominant — **1,218 of
//! 1,434** — but it missed 216, and those 216 are concentrated exactly in the
//! chunks with the most surface: 29 in the 7,439-triangle chunk, 24 in the 6,177,
//! 17 in the 4,149. Their enclosures genuinely straddle the surface band; they
//! are unnecessary because another brush **dominates** them there in the
//! `min`/`max` chain. Being distant is sufficient to be droppable, not necessary,
//! and the second route is domination. `0.849372` against the registered `0.90`.
//!
//! The demo draws that split rather than asserting it: a droppable survivor is
//! amber when it is far and magenta when it is near-and-dominated, so the magenta
//! count rising with the triangle count is something you watch happen.
//!
//! **Half of C3's predicate is structurally dead here.** The registered test
//! reads symmetrically — `lo > cell_size || hi < -cell_size` — and
//! `unnecessary_far_by_hi` is `0` on all 64 rows. It cannot fire.
//! `enclose` sets `hi = f(chunk centre) + reach`, the worst chunk's Lipschitz
//! reach is `3.4642`, so `hi < -0.125` would need `f(chunk centre) < -3.5892` —
//! and a brush field's minimum is `-radius`, which over this tape is
//! **`-0.706948`**, five times too shallow. Both numbers are recomputed from the
//! fixture at startup and printed rather than quoted, which is how the second
//! one came out sharper than the ledger's: M-358 bounds the radius by the LCG's
//! draw range, `<= 1.1`, where the largest radius the seed actually produces is
//! `0.706948`. Looser, still true, and the conclusion is unchanged.
//!
//! # What is on screen, and why the mesh is re-extracted rather than re-used
//!
//! The selected chunk is meshed **live from whichever tape is showing**. Pressing
//! `T` does not swap a colour or re-label a buffer: it rebuilds the brush slice,
//! runs Marching Cubes again, hashes the result with the crate's own
//! [`mesh_hash`](isomesh::validate::mesh_hash), and uploads it. All four hashes —
//! the three the ablation recorded and the one on screen — are on the HUD at
//! once. "Bit-identical" is therefore four 64-bit integers a viewer can compare,
//! not an adjective, and if a change ever made them differ this demo would show
//! it rather than keep claiming it.
//!
//! The brush gizmos are the tape itself, one wireframe per brush of the selected
//! chunk, filtered by the showing tape. Stepping `full → survivors → necessary`
//! empties the frame of brushes while the rock does not move a pixel. That is the
//! 20x, as a picture.
//!
//! # The fixture is P-39's, copied brush for brush
//!
//! Same `BRUSHES = 64` and the same LCG seed, the same [`Shape`] enum, the same
//! [`Interval`], [`ChunkBox`] with its central-difference margin, [`pad`],
//! [`enclose`], [`Policy::Sound`] and [`prune_into`], the same 4³ chunks of 32
//! cells at `0.125` over a sphere of radius 6. Copied out of
//! `crates/isomesh/benches/experiment_p59.rs` rather than imported, because
//! benches in this repo do not `use` one another and a shared module would let
//! one experiment's maintenance move another's published numbers.
//!
//! **This example is `f64`, and that is not a preference.** `tape_pruning`
//! (E-305) runs the same fixture in `f32` and measures a median survivor fraction
//! of `0.3281` against M-341's `0.2969`, the whole gap being `f32`'s wider
//! central-difference step inflating the chunk circumradius. P-59's numbers are
//! `f64` numbers, [`mesh_hash`](isomesh::validate::mesh_hash) is defined on
//! `MeshBuffer<f64>` only, and this demo has to reproduce a hash rather than
//! approach it — so it meshes in `f64` and narrows to `f32` on the way to the
//! vertex buffer, which is the only place the renderer is involved.
//!
//! # Every number on screen is measured in this process
//!
//! The ablation runs at startup across every core and nothing is read from the
//! CSV except P-59's own committed values, quoted for comparison. Thirteen of
//! them are cross-checked against this run and printed as
//! `expected (p-59.csv) = X, measured = X`, including the **per-chunk
//! `mesh_hash`**: all 64 of this run's survivors-only hashes are held against the
//! 64 committed ones. A run where any of the thirteen disagrees puts a
//! `CROSS-CHECK FAILED` line on the HUD instead of quietly reporting a different
//! world.

mod common;

use std::collections::BTreeMap;

use bevy::asset::RenderAssetUsages;
use bevy::mesh::{Indices, PrimitiveTopology};
use bevy::platform::time::Instant;
use bevy::prelude::*;
use common::{Capture, CommonPlugin, DemoDomain, DemoMesh, DemoStats, OrbitCamera, ViewFlags};
use isomesh::brush::{Brush, BrushOp, BrushStack, Capsule};
use isomesh::chunk::{ChunkId, ChunkLayout};
use isomesh::fields::{BoundedSdf, FieldBound, Sphere};
use isomesh::marching_cubes::MarchingCubes;
use isomesh::validate::mesh_hash;
use isomesh::{MeshBuffer, Real, RuntimeShape3, Sdf};

// ── the fixture, which is P-39's ────────────────────────────────────────────

/// Chunks along each axis. 4³ = 64 chunks, P-39's layout.
const CHUNKS_PER_AXIS: i32 = 4;

/// Cells per chunk axis, so 33 samples per axis.
const CELLS_PER_CHUNK: u32 = 32;

/// World units per cell. The chunk is 4 units across and the world 16.
///
/// A power of two, which M-32 measured is what makes two neighbouring chunks
/// agree on their shared sample plane bit-for-bit. Nothing here welds, so the
/// seam has to be exact on its own.
const CELL_SIZE: f64 = 0.125;

/// The world's minimum corner, so the world is `[-8, 8]³`.
const WORLD_ORIGIN: f64 = -8.0;

/// Brushes in the tape. **Not** reducible: every number this demo reproduces is
/// for 64, and a shorter tape would be a different fixture.
const BRUSHES: usize = 64;

/// Radius of the solid the brushes carve.
const BASE_RADIUS: f64 = 6.0;

/// ULP of slack added to every enclosure bound. P-39's constant.
const PAD_ULPS: f64 = 16.0;

/// A brush shape in this demo.
///
/// One enum so the whole stack is a single `&[Brush<Shape>]` slice, which is what
/// makes an ablated tape a shorter slice of the same type. Both variants are
/// exact distance fields, so both declare `l == 1` — and they declare it through
/// the crate's own [`BoundedSdf`] rather than by a constant written here.
#[derive(Clone, Copy, Debug)]
enum Shape {
    /// A ball.
    Sphere(Sphere<f64>),
    /// A swept ball.
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

impl BoundedSdf for Shape {
    fn value_bound(&self) -> FieldBound {
        match self {
            Self::Sphere(s) => s.value_bound(),
            Self::Capsule(c) => c.value_bound(),
        }
    }
}

impl Shape {
    /// How far below zero this shape's own field ever reaches.
    ///
    /// One half of the argument that C3's `hi < -cell_size` branch is dead on
    /// this fixture. `enclose` sets `hi = f(chunk centre) + reach`, and
    /// `f >= -depth` everywhere, so `hi >= reach - depth`; with `reach` three
    /// times the largest `depth` here the branch can never fire. Computed rather
    /// than quoted.
    ///
    /// Both variants are **exact** distance fields whose minimum is `-radius`,
    /// attained on the sphere's centre and anywhere on the capsule's segment.
    /// The capsule's half-length does **not** enter: the segment is a set of
    /// points at distance zero, not a point the distance is measured from.
    fn depth(self) -> f64 {
        match self {
            Self::Sphere(s) => s.radius,
            Self::Capsule(c) => c.radius,
        }
    }
}

/// A 64-bit LCG, so the 64 brushes are the same 64 brushes on every machine and
/// every run — and, because the seed below is P-39's, the same 64 brushes P-39
/// and P-59 measured.
struct Lcg(u64);

impl Lcg {
    /// Seed it.
    const fn new(seed: u64) -> Self {
        Self(seed)
    }

    /// The next 32 bits.
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
/// Transcribed from `experiment_p59.rs::tape`, seed included, which transcribed
/// it from `experiment_p39.rs`. Every constant here is load-bearing for
/// comparability, not for the mechanism: change one and the hashes on the HUD
/// stop being the hashes in `p-59.csv`.
fn tape() -> Vec<Brush<Shape>> {
    // P-39's seed, and nothing about the result depends on the value.
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
            // A zero direction would be a sphere, which is a fine capsule but
            // not the one this row is meant to contribute.
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

// ── the bound ───────────────────────────────────────────────────────────────

/// The enclosure of a scalar field over a chunk.
#[derive(Clone, Copy, Debug)]
struct Interval {
    /// Lower bound.
    lo: f64,
    /// Upper bound.
    hi: f64,
}

impl Interval {
    /// Whether the enclosure keeps the field further than `d` from zero
    /// everywhere in the chunk, and on which side.
    ///
    /// C3's predicate, split so the direction is visible rather than asserted —
    /// the second half is `false` everywhere on this fixture and the demo says
    /// why on screen.
    fn far_from_zero(self, d: f64) -> (bool, bool) {
        (self.lo > d, self.hi < -d)
    }
}

/// The box a chunk's field evaluations can touch.
#[derive(Clone, Copy, Debug)]
struct ChunkBox {
    /// Centre of the sampled extent.
    centre: [f64; 3],
    /// Circumradius: half the diagonal of the sampled extent, plus the margin
    /// `Sdf::gradient`'s central differences reach outside it.
    radius: f64,
}

impl ChunkBox {
    /// The box for a chunk whose sample grid starts at `origin` and spans `span`
    /// on every axis.
    fn new(origin: [f64; 3], span: f64) -> Self {
        let centre = [
            origin[0] + span * 0.5,
            origin[1] + span * 0.5,
            origin[2] + span * 0.5,
        ];
        // `h = DIFF_STEP * max(|p|, 1)` at the furthest corner bounds the
        // differencing reach anywhere in the box.
        let mut far = 1.0f64;
        for lo in origin {
            far = far.max(lo.abs()).max((lo + span).abs());
        }
        let margin = <f64 as Real>::DIFF_STEP * far;
        let half = span * 0.5 + margin;
        Self {
            centre,
            radius: half * 3.0f64.sqrt(),
        }
    }
}

/// Slack for one bound, in absolute units.
fn pad(value: f64, reach: f64) -> f64 {
    PAD_ULPS * f64::EPSILON * (value.abs() + reach)
}

/// `f(centre) ± l·r`, widened so it is an enclosure and not an estimate.
fn enclose<S: BoundedSdf<Scalar = f64>>(field: &S, chunk: ChunkBox) -> Interval {
    let l = field
        .value_bound()
        .lipschitz()
        .expect("every field in this fixture declares a Lipschitz constant");
    let value = field.sample(chunk.centre);
    let reach = l * chunk.radius;
    let slack = reach + pad(value, reach);
    Interval {
        lo: value - slack,
        hi: value + slack,
    }
}

/// Which pruning rule to apply.
///
/// P-39 has a second variant, `PruneSmoothLosers`, which is its negative control
/// for the smooth-min asymmetry and prunes in a direction that is *not*
/// bit-exact. It is not this fixture and is not transcribed, so the one variant
/// here is the registered rule.
#[derive(Clone, Copy, Debug)]
enum Policy {
    /// The registered rule. `Add` and `Subtract` prune in the losing direction;
    /// a `SmoothAdd` never prunes, because `b + (a − b)` is not `a`.
    Sound,
}

/// What one pruning pass found.
#[derive(Clone, Copy, Debug, Default)]
struct PruneStats {
    /// Brushes the bound could not rule out.
    survivors: usize,
    /// `Add`s that provably *win* everywhere in the chunk, so the whole tape
    /// prefix and the base field are dead. Counted, not exploited, exactly as in
    /// P-39: `BrushStack` has no way to say "start from this brush".
    dominant_adds: usize,
}

/// Select the brushes that can still change the fold anywhere in `chunk`,
/// recording which index of the original tape each survivor came from.
///
/// Order is preserved, because `Add` and `Subtract` do not commute with each
/// other — `BrushOp::commutes_with` is the crate's own statement of that.
fn prune_into(
    tape: &[Brush<Shape>],
    base: &Sphere<f64>,
    chunk: ChunkBox,
    policy: Policy,
    kept: &mut Vec<usize>,
) -> PruneStats {
    kept.clear();
    let mut stats = PruneStats::default();
    let mut v = enclose(base, chunk);
    for (index, brush) in tape.iter().enumerate() {
        let s = enclose(&brush.shape, chunk);
        match brush.op {
            BrushOp::Add => {
                // Strictly above the running value everywhere, so `min` is
                // forced to select the running value at every point.
                if s.lo > v.hi {
                    continue;
                }
                if s.hi < v.lo {
                    stats.dominant_adds += 1;
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
                // `Policy::Sound` never prunes a `SmoothAdd`, because `smin` at
                // `h == 1` returns `b + (a − b)`, which is not bit-identical to
                // `a`. This binding is the compile-time statement of that: a
                // second `Policy` variant would stop building right here rather
                // than silently inherit the sound rule.
                let Policy::Sound = policy;
                // `smin(a, b) = min(a, b) − k(1 − |a − b|/k)²/4`, so the floor
                // sags by at most `k/4` below the plain minimum.
                let lo = v.lo.min(s.lo) - 0.25 * k;
                v = Interval {
                    lo: lo - pad(lo, 0.25 * k),
                    hi: v.hi.min(s.hi),
                };
            }
        }
        kept.push(index);
        stats.survivors += 1;
    }
    stats
}

// ── what the ablation decided ───────────────────────────────────────────────

/// What leave-one-out decided about one brush of the tape, inside one chunk.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Role {
    /// The interval bound ruled it out before any meshing. P-39's win.
    Pruned,
    /// A survivor whose individual removal changes `mesh_hash`. The only brushes
    /// a chunk actually needs.
    Necessary,
    /// A droppable survivor whose enclosure over the chunk is more than one cell
    /// clear of zero. C3's named cause.
    DroppableFar,
    /// A droppable survivor whose enclosure straddles the surface band, so it is
    /// dominated in the `min`/`max` chain by another brush. The cause the
    /// registration did not name, and the reason C3 failed.
    DroppableDominated,
}

/// Which of the three tape lengths is showing.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Tape {
    /// All 64 brushes, the whole edit history.
    Full,
    /// What P-39's Lipschitz bound keeps. Median 19.
    Survivors,
    /// What leave-one-out proves is needed. Median 0, world maximum 16.
    Necessary,
}

/// The three tape lengths, in the order the tour walks them.
const TAPES: [Tape; 3] = [Tape::Full, Tape::Survivors, Tape::Necessary];

impl Tape {
    /// Whether a brush in this role is in this tape.
    fn keeps(self, role: Role) -> bool {
        match self {
            Self::Full => true,
            Self::Survivors => role != Role::Pruned,
            Self::Necessary => role == Role::Necessary,
        }
    }

    /// The label on the HUD.
    fn name(self) -> &'static str {
        match self {
            Self::Full => "FULL",
            Self::Survivors => "SURVIVORS",
            Self::Necessary => "NECESSARY",
        }
    }

    /// The next one round the cycle, for `T`.
    fn next(self) -> Self {
        match self {
            Self::Full => Self::Survivors,
            Self::Survivors => Self::Necessary,
            Self::Necessary => Self::Full,
        }
    }
}

/// Everything the ablation found in one chunk.
#[derive(Clone)]
struct Ablation {
    /// Chunk coordinates.
    id: [i32; 3],
    /// One role per brush of the 64-brush tape, in tape order.
    roles: Vec<Role>,
    /// Brushes the bound kept.
    survivors: usize,
    /// Survivors whose individual removal changes `mesh_hash`.
    necessary: usize,
    /// Droppable survivors more than a cell clear of zero over the chunk.
    droppable_far: usize,
    /// Droppable survivors that straddle the band and are merely dominated.
    droppable_dominated: usize,
    /// Droppable survivors caught by the `hi < -cell_size` half of C3's
    /// predicate. Zero on every chunk of this fixture, structurally.
    far_by_hi: usize,
    /// `Add`s that provably win everywhere here.
    dominant_adds: usize,
    /// The full 64-brush tape's hash — C1's control.
    hash_full: u64,
    /// The survivors-only reference every other hash is compared against.
    hash_survivors: u64,
    /// The hash of the `necessary` brushes alone, every individually-unnecessary
    /// survivor removed at once.
    hash_necessary: u64,
    /// Vertices of the reference mesh.
    vertices: usize,
    /// Triangles of the reference mesh. Zero makes `necessary == 0` free, which
    /// is why it is a field rather than an argument.
    triangles: usize,
}

impl Ablation {
    /// Survivors that proved droppable.
    fn droppable(&self) -> usize {
        self.survivors - self.necessary
    }

    /// Re-meshes this chunk cost to decide: one per survivor, plus the control.
    fn remeshes(&self) -> usize {
        self.survivors + 1
    }

    /// C1 on this chunk: removing every non-survivor at once changed nothing.
    fn control_unchanged(&self) -> bool {
        self.hash_full == self.hash_survivors
    }

    /// The joint claim: dropping every individually-unnecessary survivor at once
    /// still meshes to the same bytes.
    fn necessary_only_unchanged(&self) -> bool {
        self.hash_necessary == self.hash_survivors
    }

    /// `x-y-z`, the CSV's own chunk label.
    fn label(&self) -> String {
        format!("{}-{}-{}", self.id[0], self.id[1], self.id[2])
    }

    /// How many brushes a tape length keeps here.
    fn count(&self, tape: Tape) -> usize {
        match tape {
            Tape::Full => BRUSHES,
            Tape::Survivors => self.survivors,
            Tape::Necessary => self.necessary,
        }
    }

    /// The hash the ablation recorded for a tape length.
    fn hash(&self, tape: Tape) -> u64 {
        match tape {
            Tape::Full => self.hash_full,
            Tape::Survivors => self.hash_survivors,
            Tape::Necessary => self.hash_necessary,
        }
    }
}

/// Geometry handed back from a worker thread, already narrowed to `f32`.
///
/// The extraction is `f64` because `mesh_hash` is and because P-59's numbers
/// are; the renderer wants `f32`, so the narrowing happens exactly once, here, on
/// the way out of the measurement.
#[derive(Default)]
struct MeshData {
    /// Vertex positions.
    positions: Vec<[f32; 3]>,
    /// Vertex normals.
    normals: Vec<[f32; 3]>,
    /// Triangle indices.
    indices: Vec<u32>,
}

impl MeshData {
    /// Narrow a `f64` extraction into the renderer's `f32`.
    fn of(buffer: &MeshBuffer<f64>) -> Self {
        Self {
            positions: buffer.positions.iter().map(narrow).collect(),
            normals: buffer.normals.iter().map(narrow).collect(),
            indices: buffer.indices.clone(),
        }
    }

    /// A Bevy asset carrying the flat normals the extractor produced.
    fn to_mesh(&self) -> Mesh {
        let mut mesh = Mesh::new(
            PrimitiveTopology::TriangleList,
            RenderAssetUsages::default(),
        );
        mesh.insert_attribute(Mesh::ATTRIBUTE_POSITION, self.positions.clone());
        mesh.insert_attribute(Mesh::ATTRIBUTE_NORMAL, self.normals.clone());
        mesh.insert_indices(Indices::U32(self.indices.clone()));
        mesh
    }
}

/// `[f64; 3]` to `[f32; 3]`.
fn narrow(v: &[f64; 3]) -> [f32; 3] {
    [v[0] as f32, v[1] as f32, v[2] as f32]
}

/// A `f64` world point where the renderer wants a `f32` one.
fn place(p: [f64; 3]) -> Vec3 {
    Vec3::new(p[0] as f32, p[1] as f32, p[2] as f32)
}

// ── the measurement ─────────────────────────────────────────────────────────

/// One chunk's sample grid.
struct Grid {
    /// Sample extents.
    shape: RuntimeShape3,
    /// World position of sample `(0, 0, 0)`.
    origin: [f64; 3],
    /// World units per cell.
    cell: f64,
}

/// Buffers that outlive one worker's slice of the world.
struct Rig {
    /// The extractor.
    mc: MarchingCubes<f64>,
    /// The reference mesh, and the one the renderer reads.
    reference: MeshBuffer<f64>,
    /// Every other extraction, overwritten in place.
    scratch: MeshBuffer<f64>,
    /// Indices into the full tape, one per survivor, in tape order.
    kept: Vec<usize>,
    /// The tape slice handed to `BrushStack`.
    slice: Vec<Brush<Shape>>,
}

impl Rig {
    /// A rig with room for the whole tape, so the ablation does not spend its
    /// time in the allocator.
    fn new() -> Self {
        Self {
            mc: MarchingCubes::new(),
            reference: MeshBuffer::new(),
            scratch: MeshBuffer::new(),
            kept: Vec::with_capacity(BRUSHES),
            slice: Vec::with_capacity(BRUSHES),
        }
    }
}

/// The tape and the solid it carves.
struct Fixture {
    /// The base sphere.
    base: Sphere<f64>,
    /// All 64 brushes.
    hard: Vec<Brush<Shape>>,
}

impl Fixture {
    /// How far below zero any brush in the tape reaches at its own centre.
    fn deepest(&self) -> f64 {
        self.hard
            .iter()
            .map(|b| b.shape.depth())
            .fold(0.0f64, f64::max)
    }
}

/// Mesh `field` over `grid` into `out` and return `mesh_hash` of the result.
fn hash_of<F: Sdf<Scalar = f64>>(
    mc: &mut MarchingCubes<f64>,
    out: &mut MeshBuffer<f64>,
    field: &F,
    grid: &Grid,
) -> u64 {
    out.reset();
    mc.extract(field, &grid.shape, grid.origin, grid.cell, out)
        .expect("chunk extraction");
    mesh_hash(out)
}

/// Leave-one-out over one chunk, and the reference mesh it leaves behind.
///
/// The order of operations is P-59's, because the numbers have to be P-59's:
/// prune, mesh the survivors as the reference, mesh the whole tape as C1's
/// control, then drop each surviving brush in turn *keeping the order of the
/// rest* and compare. Order is preserved because `Add` and `Subtract` do not
/// commute, so a reordered tape would be a different field and the comparison
/// would measure the reorder instead of the removal.
fn measure_chunk(
    rig: &mut Rig,
    fixture: &Fixture,
    layout: &ChunkLayout<f64>,
    id: ChunkId,
) -> (Ablation, MeshData) {
    let Rig {
        mc,
        reference,
        scratch,
        kept,
        slice,
    } = rig;

    let origin = layout.sample_origin(id);
    let span = f64::from(layout.cells()) * layout.cell_size();
    let chunk = ChunkBox::new(origin, span);
    let cell_size = layout.cell_size();
    let grid = Grid {
        shape: layout.sample_shape().expect("chunk sample grid fits u32"),
        origin,
        cell: cell_size,
    };

    let stats = prune_into(&fixture.hard, &fixture.base, chunk, Policy::Sound, kept);

    // The survivors-only reference. Every hash below is compared against it, and
    // it is what the renderer draws.
    slice.clear();
    slice.extend(kept.iter().map(|&i| fixture.hard[i]));
    let hash_survivors = hash_of(
        mc,
        reference,
        &BrushStack {
            base: fixture.base,
            brushes: slice.as_slice(),
        },
        &grid,
    );
    let vertices = reference.vertex_count();
    let triangles = reference.triangle_count();
    let mesh = MeshData::of(reference);

    // C1 first, because if it fails every other number here is void: the whole
    // set of non-survivors removed at once.
    let hash_full = hash_of(
        mc,
        scratch,
        &BrushStack {
            base: fixture.base,
            brushes: fixture.hard.as_slice(),
        },
        &grid,
    );

    let mut roles = vec![Role::Pruned; BRUSHES];
    let mut necessary = 0usize;
    let mut droppable_far = 0usize;
    let mut droppable_dominated = 0usize;
    let mut far_by_hi = 0usize;
    for position in 0..kept.len() {
        slice.clear();
        slice.extend(
            kept.iter()
                .enumerate()
                .filter(|(at, _)| *at != position)
                .map(|(_, &i)| fixture.hard[i]),
        );
        let ablated = hash_of(
            mc,
            scratch,
            &BrushStack {
                base: fixture.base,
                brushes: slice.as_slice(),
            },
            &grid,
        );

        let index = kept[position];
        if ablated == hash_survivors {
            // C3's predicate, on the same enclosure the pruner itself consumed
            // rather than a second bound computed differently.
            let (far_lo, far_hi) =
                enclose(&fixture.hard[index].shape, chunk).far_from_zero(cell_size);
            far_by_hi += usize::from(far_hi);
            if far_lo || far_hi {
                roles[index] = Role::DroppableFar;
                droppable_far += 1;
            } else {
                roles[index] = Role::DroppableDominated;
                droppable_dominated += 1;
            }
        } else {
            roles[index] = Role::Necessary;
            necessary += 1;
        }
    }

    // The joint claim, which leave-one-out cannot support on its own: every
    // individually-unnecessary survivor removed *at once*, order preserved.
    slice.clear();
    slice.extend(
        kept.iter()
            .filter(|&&i| roles[i] == Role::Necessary)
            .map(|&i| fixture.hard[i]),
    );
    let hash_necessary = hash_of(
        mc,
        scratch,
        &BrushStack {
            base: fixture.base,
            brushes: slice.as_slice(),
        },
        &grid,
    );

    (
        Ablation {
            id: id.coords,
            roles,
            survivors: stats.survivors,
            necessary,
            droppable_far,
            droppable_dominated,
            far_by_hi,
            dominant_adds: stats.dominant_adds,
            hash_full,
            hash_survivors,
            hash_necessary,
            vertices,
            triangles,
        },
        mesh,
    )
}

/// Whether the ablation has finished. Every interactive system waits on it.
///
/// The ablation is 1,571 re-meshes and it runs on the app's own thread. Doing it
/// all inside `Startup` freezes the frame for twenty seconds, and fanning it
/// across `std::thread::scope` — which is what this demo used to do — is not
/// available at all on the web build, where thread spawn panics and
/// `SharedArrayBuffer` needs COOP/COEP headers a static host cannot send. So it
/// runs one chunk per frame, and nothing that reads the result exists until it is
/// done.
#[derive(States, Debug, Clone, Copy, PartialEq, Eq, Hash, Default)]
enum Phase {
    /// Sweeping chunks. Only [`measure_next_chunk`] and the progress line run.
    #[default]
    Measuring,
    /// Measured, cross-checked, and the scene built.
    Ready,
}

/// The one line on screen while the ablation sweeps.
#[derive(Component)]
struct Progress;

/// The ablation, one chunk per frame.
///
/// **A count budget rather than a wall-clock one**, because a count is identical
/// on every machine and a stopwatch is not -- M-348's rule. It is also why the
/// progress line quotes chunks and re-meshes rather than a rate.
#[derive(Resource)]
struct Ablating {
    /// Chunk lattice.
    layout: ChunkLayout<f64>,
    /// The base sphere and the 64 brushes.
    fixture: Fixture,
    /// Sweep order, parallel to [`Ablating::rows`] as far as it has got.
    ids: Vec<ChunkId>,
    /// Next index into [`Ablating::ids`], and so the count already measured.
    next: usize,
    /// The extractor and buffers every chunk reuses.
    rig: Rig,
    /// What the ablation has found so far.
    rows: Vec<Ablation>,
    /// The survivors-only geometry, parallel to [`Ablating::rows`].
    geometry: Vec<MeshData>,
    /// Measured seconds, summed over the chunks rather than taken from the wall
    /// clock, so `World::seconds` still reports the cost of the *measurement*
    /// rather than the 64 frames it was spread over.
    elapsed: f64,
}

impl Ablating {
    /// A sweep over `ids` that has measured nothing yet.
    fn new(layout: ChunkLayout<f64>, fixture: Fixture, ids: Vec<ChunkId>) -> Self {
        let chunks = ids.len();
        Self {
            layout,
            fixture,
            ids,
            next: 0,
            rig: Rig::new(),
            rows: Vec::with_capacity(chunks),
            geometry: Vec::with_capacity(chunks),
            elapsed: 0.0,
        }
    }

    /// The progress line, quoting the *committed* re-mesh count.
    ///
    /// From the ledger rather than a literal: the figure on screen before
    /// anything is measured can only honestly be the one P-59 published, and the
    /// measured one is what the cross-check holds it against afterwards.
    fn line(&self, ledger: &Ledger) -> String {
        format!(
            "measuring chunk {} of {} -- {} re-meshes, the cost this finding charges for",
            self.next,
            self.ids.len(),
            commas(ledger.total_remeshes)
        )
    }
}

/// The largest Lipschitz reach any chunk of this world has.
///
/// Recomputed from [`ChunkBox`] rather than quoted, because it is one half of the
/// argument that C3's `hi` branch is structurally dead and a quoted constant
/// could drift away from the geometry it describes.
fn worst_reach() -> f64 {
    let span = f64::from(CELLS_PER_CHUNK) * CELL_SIZE;
    let mut worst = 0.0f64;
    for z in 0..CHUNKS_PER_AXIS {
        for y in 0..CHUNKS_PER_AXIS {
            for x in 0..CHUNKS_PER_AXIS {
                let origin = [
                    WORLD_ORIGIN + f64::from(x) * span,
                    WORLD_ORIGIN + f64::from(y) * span,
                    WORLD_ORIGIN + f64::from(z) * span,
                ];
                worst = worst.max(ChunkBox::new(origin, span).radius);
            }
        }
    }
    worst
}

/// The aggregates this run measured, over all 64 chunks.
struct World {
    /// Chunks measured.
    chunks: usize,
    /// Upper median of the survivor counts — P-39's and P-59's convention.
    median_survivors: usize,
    /// Survivors summed over the world.
    total_survivors: usize,
    /// Necessary brushes summed over the world.
    total_necessary: usize,
    /// Droppable survivors summed over the world.
    total_droppable: usize,
    /// Droppable survivors that are far from the surface.
    total_droppable_far: usize,
    /// Droppable survivors caught by the `hi` half of C3's predicate.
    total_far_by_hi: usize,
    /// Re-meshes the whole ablation cost.
    total_remeshes: usize,
    /// Chunks where C1 held.
    control_unchanged: usize,
    /// Chunks where the joint drop held.
    necessary_only_unchanged: usize,
    /// Chunks whose reference mesh is empty, where `necessary == 0` is free.
    empty_chunks: usize,
    /// Chunks that need no brush at all.
    chunks_needing_none: usize,
    /// Seconds the ablation took.
    seconds: f64,
    /// The worst chunk circumradius in the world, which is the Lipschitz reach.
    reach: f64,
    /// How far below zero the deepest brush in the tape reaches.
    deepest_brush: f64,
}

impl World {
    /// Fold the per-chunk results into the numbers the HUD and the cross-check
    /// read.
    fn of(rows: &[Ablation], fixture: &Fixture, seconds: f64) -> Self {
        // The upper median, `sorted[len / 2]`, which is what `experiment_p59.rs`
        // and `experiment_p39.rs` both use. Copying the convention rather than
        // re-deriving it is what keeps this comparable with M-341.
        let mut survivors: Vec<usize> = rows.iter().map(|r| r.survivors).collect();
        survivors.sort_unstable();
        Self {
            chunks: rows.len(),
            median_survivors: survivors[survivors.len() / 2],
            total_survivors: rows.iter().map(|r| r.survivors).sum(),
            total_necessary: rows.iter().map(|r| r.necessary).sum(),
            total_droppable: rows.iter().map(Ablation::droppable).sum(),
            total_droppable_far: rows.iter().map(|r| r.droppable_far).sum(),
            total_far_by_hi: rows.iter().map(|r| r.far_by_hi).sum(),
            total_remeshes: rows.iter().map(Ablation::remeshes).sum(),
            control_unchanged: rows.iter().filter(|r| r.control_unchanged()).count(),
            necessary_only_unchanged: rows.iter().filter(|r| r.necessary_only_unchanged()).count(),
            empty_chunks: rows.iter().filter(|r| r.triangles == 0).count(),
            chunks_needing_none: rows.iter().filter(|r| r.necessary == 0).count(),
            seconds,
            reach: worst_reach(),
            deepest_brush: fixture.deepest(),
        }
    }

    /// The median survivor *fraction*, which is the figure M-341 published as
    /// `0.2969`.
    fn median_fraction(&self) -> f64 {
        self.median_survivors as f64 / BRUSHES as f64
    }

    /// `total_droppable_far / total_droppable`, the aggregate that decided C3.
    fn far_fraction(&self) -> f64 {
        self.total_droppable_far as f64 / self.total_droppable as f64
    }

    /// How much smaller the world's necessary tape is than its surviving one.
    fn trim(&self) -> f64 {
        self.total_survivors as f64 / self.total_necessary as f64
    }

    /// What `f(chunk centre)` would have to be for C3's `hi` branch to fire.
    fn hi_branch_threshold(&self) -> f64 {
        -(CELL_SIZE + self.reach)
    }
}

// ── the ledger, compiled in ─────────────────────────────────────────────────

/// P-59's committed artefact, embedded at compile time.
///
/// `include_str!` rather than transcribed constants: the path resolves against
/// this source file so no working directory can break it, and a number that lived
/// only here could drift away from the CSV with nothing to say so.
const LEDGER_CSV: &str = include_str!("../../docs/experiments/p-59.csv");

/// The committed values this demo holds itself against.
#[derive(Resource)]
struct Ledger {
    /// Rows in the file.
    chunks: usize,
    /// The `survivors_median` column, constant across rows.
    median_survivors: usize,
    /// `sum(survivors)`.
    total_survivors: usize,
    /// `sum(necessary)`.
    total_necessary: usize,
    /// `sum(unnecessary)`.
    total_droppable: usize,
    /// `sum(unnecessary_far_from_surface)`.
    total_droppable_far: usize,
    /// `sum(unnecessary_far_by_hi)`, which is zero.
    total_far_by_hi: usize,
    /// `sum(remeshes)`.
    total_remeshes: usize,
    /// Rows with `control_hash_unchanged = true`.
    control_unchanged: usize,
    /// Rows with `necessary_only_hash_unchanged = true`.
    necessary_only_unchanged: usize,
    /// The registered C3 bound, from the `c3_bound` column.
    c3_bound: f64,
    /// The aggregate that decided C3, from the `c3_far_fraction` column.
    c3_far_fraction: f64,
    /// The `mesh_hash` of every chunk, by the CSV's own label.
    hash_by_chunk: BTreeMap<String, u64>,
}

impl Ledger {
    /// Read every column this demo quotes out of the committed CSV.
    ///
    /// Column *names*, never positions: the file carries 47 of them and a
    /// positional read would break silently the first time one is inserted.
    fn load() -> Self {
        let mut lines = LEDGER_CSV.lines().filter(|l| !l.starts_with('#'));
        let header: Vec<&str> = lines
            .next()
            .expect("p-59.csv has a header row")
            .split(',')
            .collect();
        let column = |name: &str| {
            header
                .iter()
                .position(|h| *h == name)
                .unwrap_or_else(|| panic!("p-59.csv has no `{name}` column"))
        };
        let c_chunk = column("chunk");
        let c_survivors = column("survivors");
        let c_necessary = column("necessary");
        let c_droppable = column("unnecessary");
        let c_far = column("unnecessary_far_from_surface");
        let c_far_hi = column("unnecessary_far_by_hi");
        let c_remeshes = column("remeshes");
        let c_control = column("control_hash_unchanged");
        let c_joint = column("necessary_only_hash_unchanged");
        let c_median = column("survivors_median");
        let c_bound = column("c3_bound");
        let c_fraction = column("c3_far_fraction");
        let c_hash = column("mesh_hash");

        let mut out = Self {
            chunks: 0,
            median_survivors: 0,
            total_survivors: 0,
            total_necessary: 0,
            total_droppable: 0,
            total_droppable_far: 0,
            total_far_by_hi: 0,
            total_remeshes: 0,
            control_unchanged: 0,
            necessary_only_unchanged: 0,
            c3_bound: 0.0,
            c3_far_fraction: 0.0,
            hash_by_chunk: BTreeMap::new(),
        };
        for line in lines {
            let cells: Vec<&str> = line.split(',').collect();
            assert_eq!(
                cells.len(),
                header.len(),
                "p-59.csv row `{}` is not as wide as the header",
                cells.first().copied().unwrap_or_default()
            );
            let count = |at: usize| -> usize {
                cells[at]
                    .parse()
                    .unwrap_or_else(|_| panic!("p-59.csv column {at} is not an integer"))
            };
            out.chunks += 1;
            out.total_survivors += count(c_survivors);
            out.total_necessary += count(c_necessary);
            out.total_droppable += count(c_droppable);
            out.total_droppable_far += count(c_far);
            out.total_far_by_hi += count(c_far_hi);
            out.total_remeshes += count(c_remeshes);
            out.control_unchanged += usize::from(cells[c_control] == "true");
            out.necessary_only_unchanged += usize::from(cells[c_joint] == "true");
            // Repeated identically on every row by the harness, so the last one
            // read is the value — and a row that disagreed would be a broken file
            // rather than a second opinion.
            out.median_survivors = cells[c_median]
                .parse::<f64>()
                .expect("p-59.csv survivors_median is a float")
                as usize;
            out.c3_bound = cells[c_bound]
                .parse()
                .expect("p-59.csv c3_bound is a float");
            out.c3_far_fraction = cells[c_fraction]
                .parse()
                .expect("p-59.csv c3_far_fraction is a float");
            out.hash_by_chunk.insert(
                cells[c_chunk].to_string(),
                cells[c_hash].parse().expect("p-59.csv mesh_hash is a u64"),
            );
        }
        assert!(out.chunks > 0, "p-59.csv carries no rows");
        out
    }
}

/// One line of the cross-check: what the CSV says and what this run measured.
#[derive(Clone)]
struct Check {
    /// What is being compared.
    name: &'static str,
    /// The committed value, formatted.
    expected: String,
    /// This run's value, formatted through the same expression.
    measured: String,
}

impl Check {
    /// Both sides are rendered by one expression, so the comparison is a string
    /// equality and cannot pass on a rounding difference the printout hides.
    fn holds(&self) -> bool {
        self.expected == self.measured
    }

    /// `name expected -> measured`, for the HUD.
    fn line(&self) -> String {
        format!("{} {} -> {}", self.name, self.expected, self.measured)
    }
}

/// Every number this demo reproduces from `p-59.csv`.
#[derive(Resource)]
struct CrossCheck(Vec<Check>);

impl CrossCheck {
    /// Build every comparison. The first eight are the ones the demo is required
    /// to reproduce; the rest are the C3 aggregates and the cost, which the entry
    /// quotes and which would otherwise sit on the HUD unchecked.
    fn of(ledger: &Ledger, world: &World, hash_matches: usize) -> Self {
        let pair = |name, expected: String, measured: String| Check {
            name,
            expected,
            measured,
        };
        let share = |n: usize, total: usize| format!("{n}/{total}");
        Self(vec![
            pair(
                "chunks",
                ledger.chunks.to_string(),
                world.chunks.to_string(),
            ),
            pair(
                "median survivors",
                ledger.median_survivors.to_string(),
                world.median_survivors.to_string(),
            ),
            pair(
                "median fraction",
                format!("{:.4}", ledger.median_survivors as f64 / BRUSHES as f64),
                format!("{:.4}", world.median_fraction()),
            ),
            pair(
                "sum(survivors)",
                ledger.total_survivors.to_string(),
                world.total_survivors.to_string(),
            ),
            pair(
                "sum(necessary)",
                ledger.total_necessary.to_string(),
                world.total_necessary.to_string(),
            ),
            pair(
                "control",
                share(ledger.control_unchanged, ledger.chunks),
                share(world.control_unchanged, world.chunks),
            ),
            pair(
                "necessary-only",
                share(ledger.necessary_only_unchanged, ledger.chunks),
                share(world.necessary_only_unchanged, world.chunks),
            ),
            pair(
                "per-chunk mesh_hash",
                share(ledger.chunks, ledger.chunks),
                share(hash_matches, world.chunks),
            ),
            pair(
                "sum(unnecessary)",
                ledger.total_droppable.to_string(),
                world.total_droppable.to_string(),
            ),
            pair(
                "sum(unnec_far)",
                ledger.total_droppable_far.to_string(),
                world.total_droppable_far.to_string(),
            ),
            pair(
                "c3_far_fraction",
                format!("{:.6}", ledger.c3_far_fraction),
                format!("{:.6}", world.far_fraction()),
            ),
            pair(
                "far_by_hi",
                ledger.total_far_by_hi.to_string(),
                world.total_far_by_hi.to_string(),
            ),
            pair(
                "sum(remeshes)",
                ledger.total_remeshes.to_string(),
                world.total_remeshes.to_string(),
            ),
        ])
    }

    /// How many of the comparisons reproduced.
    fn held(&self) -> usize {
        self.0.iter().filter(|c| c.holds()).count()
    }

    /// Whether every one of them did.
    fn all_hold(&self) -> bool {
        self.held() == self.0.len()
    }
}

// ── the tour ────────────────────────────────────────────────────────────────

/// The five chunks the tour stops on, and why each one is here.
///
/// Chosen from `p-59.csv` for what they say rather than for how they look, and
/// ordered so the caveat lands **before** the payoff rather than after it:
///
/// 1. `1-1-1` — the busiest chunk in the world, 7,439 triangles. 64 survivors, 15
///    necessary, and 29 of its 49 droppable survivors are the *dominated* kind
///    C3 did not name.
/// 2. `0-1-1` — the world's worst case for necessity: 16 of 64, the largest
///    `necessary_fraction` in the file at `0.25`.
/// 3. `1-2-0` — C3's named cause almost pure: 41 of 42 droppable survivors are
///    more than a cell clear of zero.
/// 4. `0-3-0` — the honest confound, and it comes fourth on purpose. 30
///    survivors, `necessary == 0`, and an **empty** mesh: a chunk with no surface
///    hashes the same after any removal, so its zero is free rather than earned.
///    A viewer who meets this one *after* the payoff would reasonably suspect
///    every zero of being this.
/// 5. `2-1-0` — the payoff. 3,227 triangles, 35 survivors, and a `necessary` tape
///    of length **zero**: the bare sphere meshes to the same bytes as the sphere
///    carved by 64 edits, and this time there is a surface to be wrong about.
const STOPS: [[i32; 3]; 5] = [[1, 1, 1], [0, 1, 1], [1, 2, 0], [0, 3, 0], [2, 1, 0]];

/// Segments of the story: one per stop per tape length.
const SEGMENTS: usize = STOPS.len() * TAPES.len();

/// Seconds for one pass through the tour, when nobody is capturing.
const STORY_SECONDS: f32 = 33.0;

/// What this frame is showing.
#[derive(Resource, Clone, Copy, PartialEq, Eq)]
struct Shot {
    /// Index into the chunk list.
    chunk: usize,
    /// Which tape length.
    tape: Tape,
}

/// Whether the tour is driving, and where it has got to.
#[derive(Resource, Default)]
struct Tour {
    /// Set by any key that picks a chunk or a tape; `X` clears it.
    manual: bool,
    /// Seconds of tour elapsed, when nobody is capturing.
    elapsed: f32,
}

// ── the scene ───────────────────────────────────────────────────────────────

/// Orbit radius, in world units.
///
/// Bevy's default vertical FOV is 45°, so this frames 12.4 units of height: the
/// selected chunk is 4 of them and the ball fills the frame, which is the
/// difference between looking at a world and looking at a marble. Any closer and
/// the brushes that can reach the chunk leave the shot — the Lipschitz reach is
/// 3.46 — and at 20 the chunk's own wedge of surface was six percent of a ball
/// that took half the frame, which is not a thing a viewer can find.
const CAMERA_RADIUS: f32 = 15.0;

/// How far left of the selected chunk the camera aims, as a share of the radius.
///
/// **The HUD owns the left half of a 1280-wide frame** and the orbit camera puts
/// whatever it looks at dead centre, which is exactly behind the panel's right
/// edge. Aiming this far along the camera's own screen-right axis slides the
/// subject to about 74% of the frame width, clear of the text. A share of the
/// radius rather than a distance, so scrolling to zoom keeps the composition
/// instead of walking the subject off the edge.
const FRAME_SHIFT_RATIO: f32 = 0.355;

/// How much of the way to a newly selected chunk the camera moves per frame.
///
/// Per **frame**, not per second, and deliberately: under capture the story is
/// paced by tick count, so a time-based ease would put the camera somewhere
/// slightly different in every recording of the same clip.
const CAMERA_EASE: f32 = 0.11;

/// Ambient fill, against the harness default of 220.
///
/// The subject is a grey ball with grey dents in it, and the harness's ambient
/// lights it almost flat. A third of it leaves the key light describing the
/// carve marks, which is the only reason the rock reads as destructible.
const AMBIENT_BRIGHTNESS: f32 = 75.0;

/// How far out the key light is parked from the chunk it is lighting.
///
/// A directional light's position only sets its direction, so this only has to
/// clear the world; 40 puts it well outside the 16-unit lattice at every stop.
const KEY_LIGHT_DISTANCE: f32 = 40.0;

/// Width and height of the backdrop the HUD is read against, in logical pixels.
///
/// Sized for the widest and tallest the HUD reaches — an 86-character line and 41
/// lines of them. At the harness's 13 px font the pitch is 15.6 px.
const HUD_PANEL: Vec2 = Vec2::new(726.0, 664.0);

/// Colour of a brush the bound kept and leave-one-out proved necessary.
const NEEDED_SRGB: [f32; 3] = [0.25, 0.95, 1.0];
/// Colour of a droppable survivor far from the surface — C3's named cause.
const FAR_SRGB: [f32; 3] = [1.0, 0.72, 0.18];
/// Colour of a droppable survivor that straddles the band and is dominated.
const DOMINATED_SRGB: [f32; 3] = [0.98, 0.32, 0.86];
/// Colour of a brush the interval bound deleted before any meshing.
const PRUNED_SRGB: [f32; 3] = [0.32, 0.34, 0.42];

/// Circle segments per brush gizmo.
///
/// Sixteen rather than the gizmo default of thirty-two: sixty-four brushes at
/// three circles each is 3,072 lines at this setting and twice that at the
/// default, and at the on-screen size of a one-unit ball in an 11.6-unit frame
/// the extra segments are not resolvable.
const BRUSH_SEGMENTS: u32 = 16;

/// Base colour of the chunk the tour is on.
///
/// Warm, not merely brighter. The selected chunk's share of the ball's surface
/// is a wedge of maybe a tenth of it, and a lighter grey wedge against grey rock
/// under a raking key is not a wedge anybody finds. A hue difference is.
const SELECTED_SRGB: [f32; 3] = [0.95, 0.72, 0.30];
/// Base colour of every other chunk.
const CONTEXT_SRGB: [f32; 3] = [0.30, 0.34, 0.42];

/// Gizmos for the chunk boxes and the brushes.
#[derive(Default, Reflect, GizmoConfigGroup)]
struct TapeGizmos;

/// The bottom caption — the line a viewer reads instead of the HUD.
#[derive(Component)]
struct Caption;

/// Everything measured once at startup and read every frame after.
#[derive(Resource)]
struct Demo {
    /// Chunk lattice.
    layout: ChunkLayout<f64>,
    /// The base sphere and the 64 brushes.
    fixture: Fixture,
    /// Sweep order, parallel to [`Demo::rows`].
    ids: Vec<ChunkId>,
    /// What the ablation found, one per chunk.
    rows: Vec<Ablation>,
    /// One entity per chunk.
    entities: Vec<Entity>,
    /// One material per chunk, so the selection can be lit differently.
    materials: Vec<Handle<StandardMaterial>>,
    /// One mesh asset per chunk that has a surface.
    meshes: Vec<Option<Handle<Mesh>>>,
    /// The world aggregates.
    world: World,
    /// Whether the brush gizmos are drawn.
    brushes_shown: bool,
    /// Whether the chunk boxes are drawn.
    boxes_shown: bool,
    /// The `Shot` the selected chunk's mesh was last extracted for.
    built: Option<Shot>,
    /// The hash of the mesh currently on screen, re-extracted rather than
    /// remembered.
    live_hash: u64,
    /// Triangles the live extraction produced.
    live_triangles: usize,
    /// Vertices the live extraction produced.
    live_vertices: usize,
    /// Brushes the live extraction folded.
    live_brushes: usize,
    /// Milliseconds the live extraction took.
    live_ms: f64,
    /// The extractor and buffers the live re-mesh uses.
    rig: Rig,
}

impl Demo {
    /// The chunk the shot is on.
    fn row(&self, shot: Shot) -> &Ablation {
        &self.rows[shot.chunk]
    }

    /// Index of a chunk by its coordinates.
    fn index_of(&self, coords: [i32; 3]) -> usize {
        index_of(&self.ids, coords)
    }

    /// Centre of a chunk, in world units.
    fn centre(&self, index: usize) -> Vec3 {
        chunk_centre(&self.layout, self.ids[index])
    }

    /// Build the tape a shot asks for into the rig's slice.
    fn build_tape(&mut self, shot: Shot) {
        let roles = &self.rows[shot.chunk].roles;
        let hard = &self.fixture.hard;
        self.rig.slice.clear();
        self.rig.slice.extend(
            hard.iter()
                .zip(roles)
                .filter(|&(_, &role)| shot.tape.keeps(role))
                .map(|(brush, _)| *brush),
        );
    }
}

// ── app ─────────────────────────────────────────────────────────────────────

fn main() {
    App::new()
        .add_plugins(DefaultPlugins.set(WindowPlugin {
            primary_window: Some(Window {
                title: "isomesh - E-315 edit tape trim".into(),
                // Web only, inert on native: bind to the 1280x720 canvas the
                // page supplies rather than letting Bevy append its own. The HUD
                // panels are laid out in pixels for that size, so the canvas is
                // fixed and CSS scales it -- `fit_canvas_to_parent` stays at its
                // `false` default for the same reason.
                canvas: Some("#isomesh-canvas".into()),
                ..default()
            }),
            ..default()
        }))
        .add_plugins(CommonPlugin)
        // Not black: a lit ball against a void reads as a floating polygon,
        // against a dim sky it reads as a world. Dark enough that the HUD and the
        // brush gizmos both stay legible over it.
        .insert_resource(ClearColor(Color::srgb(0.09, 0.11, 0.15)))
        .insert_resource(Shot {
            chunk: 0,
            tape: Tape::Full,
        })
        .init_resource::<Tour>()
        .init_gizmo_group::<TapeGizmos>()
        .init_state::<Phase>()
        .add_systems(Startup, setup)
        .add_systems(
            Update,
            measure_next_chunk.run_if(in_state(Phase::Measuring)),
        )
        // `PreUpdate`, not `Update`, for the reason `game_carve_seams` records:
        // the harness's `update_hud` lives in `Update` and system order within a
        // schedule is unspecified, so a caption written in `Update` disagreed
        // with the HUD above it by one frame. Two numbers on screen disagreeing
        // is worse than either being late.
        .add_systems(
            PreUpdate,
            (controls, advance, rebuild, frame_camera, report)
                .chain()
                .run_if(in_state(Phase::Ready)),
        )
        .add_systems(Update, draw_tape.run_if(in_state(Phase::Ready)))
        .run();
}

/// Where the orbit camera should look so that `subject` lands right of centre.
///
/// The orbit camera puts its focus point dead centre of frame, and the HUD owns
/// the left half of that frame. `Transform::look_at` with a `+Y` up builds its
/// screen-right as `Y × D` for a camera-to-focus direction `D`, which for this
/// camera's yaw is `(sin yaw, 0, −cos yaw)`; aiming
/// [`FRAME_SHIFT_RATIO`]`·radius` along it in the negative direction slides the
/// subject the same distance to the right. Recomputed from the live yaw and
/// radius so a human orbiting or zooming keeps the composition.
fn aim_at(subject: Vec3, yaw: f32, radius: f32) -> Vec3 {
    let (sin, cos) = yaw.sin_cos();
    subject - Vec3::new(sin, 0.0, -cos) * (radius * FRAME_SHIFT_RATIO)
}

/// Index of a chunk by its coordinates, in a sweep-ordered id list.
fn index_of(ids: &[ChunkId], coords: [i32; 3]) -> usize {
    ids.iter()
        .position(|id| id.coords == coords)
        .unwrap_or_else(|| panic!("chunk {coords:?} is not in a 4^3 world"))
}

/// Centre of a chunk, in world units.
fn chunk_centre(layout: &ChunkLayout<f64>, id: ChunkId) -> Vec3 {
    let origin = layout.sample_origin(id);
    let half = 0.5 * f64::from(layout.cells()) * layout.cell_size();
    place([origin[0] + half, origin[1] + half, origin[2] + half])
}

/// Every chunk of the world, the tour's five stops first.
///
/// **Order matters now that the sweep is incremental.** The five [`STOPS`] are
/// measured first, in tour order, so the chunks the captions talk about are the
/// ones a viewer sees appear; the remaining 59 keep the `z, y, x` order they were
/// built in. [`Demo::rows`] stays parallel to [`Demo::ids`], which every consumer
/// already assumes, so nothing downstream sees the reorder.
fn sweep_order() -> Vec<ChunkId> {
    let mut ids = Vec::with_capacity(CHUNKS_PER_AXIS.pow(3) as usize);
    ids.extend(STOPS.map(ChunkId::new));
    for z in 0..CHUNKS_PER_AXIS {
        for y in 0..CHUNKS_PER_AXIS {
            for x in 0..CHUNKS_PER_AXIS {
                let coords = [x, y, z];
                if !STOPS.contains(&coords) {
                    ids.push(ChunkId::new(coords));
                }
            }
        }
    }
    ids
}

/// Build the fixture and the sweep order, then hand the ablation to
/// [`measure_next_chunk`].
///
/// It inserts no [`Demo`]: there is nothing to draw until the sweep ends, and a
/// half-measured world on screen would be a picture of a claim nobody has checked
/// yet.
fn setup(
    mut commands: Commands,
    mut config: ResMut<GizmoConfigStore>,
    mut camera: Query<&mut OrbitCamera>,
    mut ambient: ResMut<GlobalAmbientLight>,
) {
    let layout = ChunkLayout::<f64>::new(CELLS_PER_CHUNK, CELL_SIZE, [WORLD_ORIGIN; 3])
        .expect("chunk layout is well formed");
    let fixture = Fixture {
        base: Sphere {
            center: [0.0; 3],
            radius: BASE_RADIUS,
        },
        hard: tape(),
    };
    let ids = sweep_order();
    let ledger = Ledger::load();

    info!(
        "E-315 -- running P-59's leave-one-out ablation over {} chunks, one chunk per frame; \
         this is the {} re-meshes the finding charges for",
        ids.len(),
        commas(ledger.total_remeshes)
    );

    let extent = CELL_SIZE * f64::from(CELLS_PER_CHUNK) * f64::from(CHUNKS_PER_AXIS);
    commands.spawn(DemoDomain {
        min: Vec3::splat(WORLD_ORIGIN as f32),
        max: Vec3::splat((WORLD_ORIGIN + extent) as f32),
    });

    let (tape_gizmos, _) = config.config_mut::<TapeGizmos>();
    tape_gizmos.line.width = 1.6;
    // Always in front. The tape is an overlay on the world, not part of it: most
    // of these sixty-four brushes are buried inside a solid ball six units
    // across, and depth-tested wireframes would show the handful that happen to
    // poke out. An X-ray of the edit history is the picture this demo is for.
    tape_gizmos.depth_bias = -1.0;
    ambient.brightness = AMBIENT_BRIGHTNESS;

    let first = index_of(&ids, STOPS[0]);
    let centre = chunk_centre(&layout, ids[first]);
    let (yaw, pitch) = view_of(centre);
    for mut orbit in &mut camera {
        orbit.radius = CAMERA_RADIUS;
        orbit.yaw = yaw;
        orbit.pitch = pitch;
        orbit.focus = aim_at(centre, yaw, CAMERA_RADIUS);
    }

    spawn_hud_panel(&mut commands);
    spawn_caption(&mut commands);

    let ablating = Ablating::new(layout, fixture, ids);
    spawn_progress(&mut commands, &ablating, &ledger);

    commands.insert_resource(Shot {
        chunk: first,
        tape: Tape::Full,
    });
    commands.insert_resource(ledger);
    commands.insert_resource(ablating);
}

/// The line on screen while the ablation sweeps.
///
/// The only thing on screen: there is no world to look at until the sweep ends,
/// and a blank window for sixty-four frames reads as a hung app.
fn spawn_progress(commands: &mut Commands, ablating: &Ablating, ledger: &Ledger) {
    commands.spawn((
        Text::new(ablating.line(ledger)),
        TextFont {
            font_size: FontSize::Px(20.0),
            ..default()
        },
        TextColor(Color::srgb(0.97, 0.95, 0.91)),
        Node {
            position_type: PositionType::Absolute,
            left: Val::Px(18.0),
            bottom: Val::Px(18.0),
            padding: UiRect::axes(Val::Px(14.0), Val::Px(8.0)),
            ..default()
        },
        BackgroundColor(Color::srgba(0.02, 0.03, 0.05, 0.84)),
        GlobalZIndex(4),
        Progress,
    ));
}

/// Measure one chunk, and on the last one build everything the demo reads.
///
/// One **chunk** per frame rather than one leave-one-out re-mesh: a chunk is
/// about twenty re-meshes, so the sweep finishes in sixty-four frames — long
/// enough for a viewer to watch the count run, short enough that no browser calls
/// the tab unresponsive. [`measure_chunk`] is untouched and each chunk is
/// independent — its own tape slice, its own buffers — so the numbers are the
/// same ones the threaded sweep produced, and the cross-check is what says so.
fn measure_next_chunk(
    mut commands: Commands,
    mut ablating: ResMut<Ablating>,
    ledger: Res<Ledger>,
    mut materials: ResMut<Assets<StandardMaterial>>,
    mut meshes: ResMut<Assets<Mesh>>,
    mut progress: Query<(Entity, &mut Text), With<Progress>>,
    mut next_phase: ResMut<NextState<Phase>>,
) {
    let ablating = ablating.as_mut();
    let id = ablating.ids[ablating.next];
    let started = Instant::now();
    let (row, geometry) = measure_chunk(&mut ablating.rig, &ablating.fixture, &ablating.layout, id);
    ablating.elapsed += started.elapsed().as_secs_f64();
    ablating.rows.push(row);
    ablating.geometry.push(geometry);
    ablating.next += 1;

    let line = ablating.line(&ledger);
    for (_, mut text) in &mut progress {
        text.0.clone_from(&line);
    }
    if ablating.next < ablating.ids.len() {
        return;
    }

    // Moved out field by field rather than by taking the resource whole, because
    // a system cannot take a resource by value. `Ablating` is removed below and
    // the run condition stops this system on the next frame, so nothing can
    // observe the emptied vectors.
    let layout = ablating.layout;
    let fixture = Fixture {
        base: ablating.fixture.base,
        hard: core::mem::take(&mut ablating.fixture.hard),
    };
    let ids = core::mem::take(&mut ablating.ids);
    let rows = core::mem::take(&mut ablating.rows);
    let geometry = core::mem::take(&mut ablating.geometry);
    let world = World::of(&rows, &fixture, ablating.elapsed);
    commands.remove_resource::<Ablating>();

    let hash_matches = rows
        .iter()
        .filter(|r| ledger.hash_by_chunk.get(&r.label()) == Some(&r.hash_survivors))
        .count();
    let cross = CrossCheck::of(&ledger, &world, hash_matches);
    report_to_console(&fixture, &ledger, &world, &cross);

    // The scene: one entity, one material and one mesh per chunk, built from the
    // survivors-only reference the ablation already produced.
    let mut entities = Vec::with_capacity(ids.len());
    let mut chunk_materials = Vec::with_capacity(ids.len());
    let mut chunk_meshes = Vec::with_capacity(ids.len());
    for data in &geometry {
        let material = materials.add(StandardMaterial {
            base_color: Color::srgb(CONTEXT_SRGB[0], CONTEXT_SRGB[1], CONTEXT_SRGB[2]),
            perceptual_roughness: 0.85,
            metallic: 0.02,
            ..default()
        });
        // A chunk with no surface gets no asset at all rather than an empty one:
        // `MeshAllocator` skips a zero-byte vertex buffer and then copies into it
        // anyway, which logs a use-after-free per empty chunk per frame.
        let handle = if data.indices.is_empty() {
            None
        } else {
            Some(meshes.add(data.to_mesh()))
        };
        entities.push(
            commands
                .spawn((
                    handle.clone().map_or_else(Mesh3d::default, Mesh3d),
                    MeshMaterial3d(material.clone()),
                    DemoMesh,
                ))
                .id(),
        );
        chunk_materials.push(material);
        chunk_meshes.push(handle);
    }

    for (entity, _) in &progress {
        commands.entity(entity).despawn();
    }

    commands.insert_resource(Demo {
        layout,
        fixture,
        ids,
        rows,
        entities,
        materials: chunk_materials,
        meshes: chunk_meshes,
        world,
        brushes_shown: true,
        boxes_shown: true,
        built: None,
        live_hash: 0,
        live_triangles: 0,
        live_vertices: 0,
        live_brushes: 0,
        live_ms: 0.0,
        rig: Rig::new(),
    });
    commands.insert_resource(cross);
    next_phase.set(Phase::Ready);
}

/// The report: the fixture, the cross-check, the win and the cost.
///
/// `info!` rather than `println!`, and not only for consistency with the other
/// two demos: on wasm `println!` goes to an unsupported stdout and is discarded,
/// while Bevy's `LogPlugin` routes `tracing` to `console.log`. This is the
/// evidence, and a reader who only ever sees the GIF should still be able to run
/// the example -- native or in a browser -- and read every comparison.
fn report_to_console(fixture: &Fixture, ledger: &Ledger, world: &World, cross: &CrossCheck) {
    let adds = fixture.hard.iter().filter(|b| b.op == BrushOp::Add).count();
    info!("\nE-315 game_edit_tape_trim -- M-358 / P-59 (R-057)");
    info!(
        "fixture:  {BRUSHES} brushes ({adds} Add, {} Subtract) over a sphere of radius \
         {BASE_RADIUS}; {CHUNKS_PER_AXIS}^3 chunks of {CELLS_PER_CHUNK} cells at {CELL_SIZE}",
        BRUSHES - adds
    );
    info!(
        "ablation: {} re-meshes over {} chunks in {:.2} s, one chunk per frame",
        commas(world.total_remeshes),
        world.chunks,
        world.seconds
    );

    info!("\ncross-check against docs/experiments/p-59.csv");
    for check in &cross.0 {
        info!(
            "  {:<20} expected (p-59.csv) = {:<10} measured = {:<10} {}",
            check.name,
            check.expected,
            check.measured,
            if check.holds() { "ok" } else { "MISMATCH" }
        );
    }
    if cross.all_hold() {
        info!("  all {} reproduce", cross.0.len());
    } else {
        info!(
            "  CROSS-CHECK FAILED: {} of {} reproduce -- this is not P-59's fixture",
            cross.held(),
            cross.0.len()
        );
    }

    info!(
        "\nthe win:  {} survivors -> {} necessary brushes world-wide, bit-exact on {}/{} chunks",
        commas(world.total_survivors),
        commas(world.total_necessary),
        world.necessary_only_unchanged,
        world.chunks
    );
    info!(
        "          {:.1}x on top of P-39's {BRUSHES} -> {} (M-341). {} of {} chunks need no brush \
         at all, {} of those mesh empty.",
        world.trim(),
        world.median_survivors,
        world.chunks_needing_none,
        world.chunks,
        world.empty_chunks
    );
    info!(
        "the cost: {} re-meshes to decide it. Headroom, NOT a shippable pruner.",
        commas(world.total_remeshes)
    );
    info!(
        "C3:       {} far of {} = {:.6} against the registered {:.2} -- FALSIFIED",
        commas(world.total_droppable_far),
        commas(world.total_droppable),
        world.far_fraction(),
        ledger.c3_bound
    );
    info!(
        "          the other {} droppable survivors straddle the surface band and are dominated",
        commas(world.total_droppable - world.total_droppable_far)
    );
    info!("          in the min/max chain, which is the cause the registration did not name.");
    info!(
        "          unnecessary_far_by_hi is {} of {} and cannot fire: no brush's field ever",
        world.total_far_by_hi,
        commas(world.total_droppable)
    );
    info!(
        "          goes below {:.4}, while `hi < -{CELL_SIZE}` needs f(chunk centre) < {:.4}",
        -world.deepest_brush,
        world.hi_branch_threshold()
    );
    info!(
        "          at a Lipschitz reach of {:.4}. Half the registered predicate is dead.",
        world.reach
    );
    info!("");
}

/// The dark panel the HUD is read against.
fn spawn_hud_panel(commands: &mut Commands) {
    // Behind the harness HUD, which `CommonPlugin` spawns at the default z.
    // `GlobalZIndex(-1)` is the whole mechanism: no reaching into the shared
    // module, and the panel is empty dark pixels when `nohud` clears the text.
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
}

/// The line a viewer reads instead of the HUD.
fn spawn_caption(commands: &mut Commands) {
    commands
        .spawn((
            Node {
                position_type: PositionType::Absolute,
                left: Val::Px(0.0),
                right: Val::Px(0.0),
                bottom: Val::Px(18.0),
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
                TextColor(Color::srgb(0.97, 0.95, 0.91)),
                BackgroundColor(Color::srgba(0.02, 0.03, 0.05, 0.84)),
                Node {
                    padding: UiRect::axes(Val::Px(14.0), Val::Px(8.0)),
                    ..default()
                },
                Caption,
            ));
        });
}

// ── steering ────────────────────────────────────────────────────────────────

/// `[` and `]` step chunks, `T` cycles the tape, `B` and `K` toggle the overlays,
/// `X` hands control back to the tour.
fn controls(
    keys: Res<ButtonInput<KeyCode>>,
    mut demo: ResMut<Demo>,
    mut shot: ResMut<Shot>,
    mut tour: ResMut<Tour>,
) {
    if keys.just_pressed(KeyCode::KeyB) {
        demo.brushes_shown = !demo.brushes_shown;
    }
    if keys.just_pressed(KeyCode::KeyK) {
        demo.boxes_shown = !demo.boxes_shown;
    }
    if keys.just_pressed(KeyCode::KeyX) {
        tour.manual = false;
    }

    let last = demo.rows.len() - 1;
    if keys.just_pressed(KeyCode::BracketRight) {
        shot.chunk = if shot.chunk == last {
            0
        } else {
            shot.chunk + 1
        };
        tour.manual = true;
    }
    if keys.just_pressed(KeyCode::BracketLeft) {
        shot.chunk = if shot.chunk == 0 {
            last
        } else {
            shot.chunk - 1
        };
        tour.manual = true;
    }
    if keys.just_pressed(KeyCode::KeyT) {
        shot.tape = shot.tape.next();
        tour.manual = true;
    }
}

/// Decide which chunk and which tape length this frame is about.
///
/// Under capture the tour advances with the captured frame count, so a clip of
/// any length is the whole tour rather than a slice of it. Interactively it runs
/// on wall-clock time and loops, and a digit pins one stop.
fn advance(
    time: Res<Time>,
    capture: Res<Capture>,
    flags: Res<ViewFlags>,
    demo: Res<Demo>,
    mut tour: ResMut<Tour>,
    mut shot: ResMut<Shot>,
) {
    // Steering wins over the pin, so `[`, `]` and `T` still work after a digit
    // has been pressed; `X` gives the pin back.
    if tour.manual {
        return;
    }
    // A digit pins a stop at its last tape length, which is the frame with the
    // whole argument in it: same mesh, same hash, almost no brushes.
    if (1..=STOPS.len()).contains(&flags.field) {
        *shot = Shot {
            chunk: demo.index_of(STOPS[flags.field - 1]),
            tape: Tape::Necessary,
        };
        return;
    }

    let phase = if capture.is_active() {
        f32::from(u16::try_from(capture.taken).unwrap_or(u16::MAX))
            / f32::from(u16::try_from(capture_frames()).unwrap_or(1).max(1))
    } else {
        if !flags.paused {
            tour.elapsed += time.delta_secs();
        }
        (tour.elapsed / STORY_SECONDS).fract()
    };
    // `min` rather than a wrap: `phase` reaches exactly 1.0 on the last captured
    // frame and that frame belongs to the last segment.
    let segment = ((phase.clamp(0.0, 1.0) * SEGMENTS as f32) as usize).min(SEGMENTS - 1);
    *shot = Shot {
        chunk: demo.index_of(STOPS[segment / TAPES.len()]),
        tape: TAPES[segment % TAPES.len()],
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

/// Re-extract the selected chunk from the tape that is showing.
///
/// **The mesh on screen comes from the tape on the HUD**, every time either
/// changes. Swapping a stored buffer would make "bit-identical" a claim about
/// something that happened at startup; extracting again and hashing again makes
/// it a claim about the triangles currently being rasterised.
fn rebuild(
    mut demo: ResMut<Demo>,
    shot: Res<Shot>,
    mut meshes: ResMut<Assets<Mesh>>,
    mut materials: ResMut<Assets<StandardMaterial>>,
    mut commands: Commands,
) {
    let shot = *shot;
    if demo.built == Some(shot) {
        return;
    }
    let previous = demo.built.map(|built| built.chunk);
    let demo = demo.as_mut();

    let origin = demo.layout.sample_origin(demo.ids[shot.chunk]);
    let grid = Grid {
        shape: demo
            .layout
            .sample_shape()
            .expect("chunk sample grid fits u32"),
        origin,
        cell: demo.layout.cell_size(),
    };
    demo.build_tape(shot);

    let started = Instant::now();
    demo.live_hash = hash_of(
        &mut demo.rig.mc,
        &mut demo.rig.reference,
        &BrushStack {
            base: demo.fixture.base,
            brushes: demo.rig.slice.as_slice(),
        },
        &grid,
    );
    demo.live_ms = started.elapsed().as_secs_f64() * 1000.0;
    demo.live_vertices = demo.rig.reference.vertex_count();
    demo.live_triangles = demo.rig.reference.triangle_count();
    demo.live_brushes = demo.rig.slice.len();

    let data = MeshData::of(&demo.rig.reference);
    if data.indices.is_empty() {
        if demo.meshes[shot.chunk].take().is_some() {
            commands
                .entity(demo.entities[shot.chunk])
                .insert(Mesh3d::default());
        }
    } else if let Some(handle) = &demo.meshes[shot.chunk] {
        if let Some(mut mesh) = meshes.get_mut(handle) {
            *mesh = data.to_mesh();
        }
    } else {
        let handle = meshes.add(data.to_mesh());
        commands
            .entity(demo.entities[shot.chunk])
            .insert(Mesh3d(handle.clone()));
        demo.meshes[shot.chunk] = Some(handle);
    }

    // Light the selection and hand the previous one back to the background.
    if previous != Some(shot.chunk) {
        if let Some(old) = previous
            && let Some(mut material) = materials.get_mut(&demo.materials[old])
        {
            material.base_color = Color::srgb(CONTEXT_SRGB[0], CONTEXT_SRGB[1], CONTEXT_SRGB[2]);
        }
        if let Some(mut material) = materials.get_mut(&demo.materials[shot.chunk]) {
            material.base_color = Color::srgb(SELECTED_SRGB[0], SELECTED_SRGB[1], SELECTED_SRGB[2]);
        }
    }

    demo.built = Some(shot);
}

/// Swing round to look at the selected chunk from outside the ball, then hold
/// dead still.
///
/// **Still is the point.** The frames a viewer is asked to compare — the same
/// rock under a 64-brush tape and under a two-brush one — have to be comparable,
/// and a camera drifting through the comparison makes that impossible to judge.
/// So the re-frame is a **cut**, taken on the frame the chunk changes, and every
/// frame after it is identical until the next chunk.
///
/// The cut aims along the chunk's own outward radial, so the piece of ball
/// surface the chunk owns faces the camera instead of hiding behind the far side
/// of the world. Between cuts nothing here fights a human dragging the view: the
/// yaw and pitch are written once and then left alone, and the focus follows
/// whatever yaw the mouse leaves behind.
fn frame_camera(
    demo: Res<Demo>,
    shot: Res<Shot>,
    mut aimed: Local<Option<usize>>,
    mut camera: Query<&mut OrbitCamera>,
    mut lights: Query<&mut Transform, With<DirectionalLight>>,
) {
    let centre = demo.centre(shot.chunk);
    let cut = *aimed != Some(shot.chunk);
    for mut orbit in &mut camera {
        if cut {
            let (yaw, pitch) = view_of(centre);
            orbit.yaw = yaw;
            orbit.pitch = pitch;
            orbit.radius = CAMERA_RADIUS;
            orbit.focus = aim_at(centre, yaw, CAMERA_RADIUS);
            // The key comes with the camera. A fixed light and a camera that
            // swings all the way round the world means half the stops are shot
            // against the terminator, and photographing a black ball is not a
            // demonstration of anything. Over the left shoulder and lifted, so
            // the carve marks still cast the short shadows that make them read
            // as dents rather than as texture.
            let key = (camera_direction(yaw, pitch) + Vec3::Y * 0.5
                - Vec3::new(yaw.sin(), 0.0, -yaw.cos()) * 0.32)
                .normalize();
            for mut transform in &mut lights {
                *transform = Transform::from_translation(centre + key * KEY_LIGHT_DISTANCE)
                    .looking_at(centre, Vec3::Y);
            }
        } else {
            let target = aim_at(centre, orbit.yaw, orbit.radius);
            orbit.focus = orbit.focus.lerp(target, CAMERA_EASE);
        }
    }
    *aimed = Some(shot.chunk);
}

/// The unit vector from the orbit focus toward the camera, for a yaw and pitch.
///
/// The same expression the harness's own `orbit_camera` places the camera with,
/// written out here because the key light has to know where the camera is before
/// that system next runs.
fn camera_direction(yaw: f32, pitch: f32) -> Vec3 {
    Vec3::new(
        yaw.cos() * pitch.cos(),
        pitch.sin(),
        yaw.sin() * pitch.cos(),
    )
}

/// Yaw and pitch that look at `centre` from outside the world's own ball.
///
/// The chunk lattice is centred on the origin and so is the solid, so the
/// direction from the origin through a chunk is the direction from which that
/// chunk's surface is nearest the camera and nothing else is in the way. The
/// pitch is flattened toward the horizon and lifted a little: a raw radial
/// through a corner chunk puts the camera 35° under the world, which reads as a
/// mistake rather than as a choice.
fn view_of(centre: Vec3) -> (f32, f32) {
    // No chunk of a 4³ lattice centred on the origin sits at the origin — the
    // centres are at ±2 and ±6 on every axis — so this cannot fail, and if the
    // layout ever changed so that it could, a panic naming the reason beats a
    // silent NaN yaw.
    let direction = centre
        .try_normalize()
        .expect("no chunk of this lattice is centred on the world origin");
    (
        direction.z.atan2(direction.x),
        (direction.y.asin() * 0.55 + 0.22).clamp(-0.9, 1.0),
    )
}

// ── what is on screen ───────────────────────────────────────────────────────

/// The tape itself: one wireframe per brush of the selected chunk, filtered by
/// the tape that is showing.
///
/// This is the 20x as a picture. In `FULL` all sixty-four are drawn and the
/// pruned ones are grey; in `SURVIVORS` the grey ones vanish; in `NECESSARY`
/// everything but the cyan vanishes — and the rock does not move.
fn draw_tape(demo: Res<Demo>, shot: Res<Shot>, mut gizmos: Gizmos<TapeGizmos>) {
    if demo.boxes_shown {
        // The selected chunk only. Sixty-four always-on-top wireframe cubes over
        // one 16-unit world is a cage the subject cannot be seen through, and the
        // harness already draws the world's own extent on `G`.
        let span = (demo.layout.cell_size() * f64::from(demo.layout.cells())) as f32;
        gizmos.cube(
            Transform::from_translation(demo.centre(shot.chunk))
                .with_scale(Vec3::splat(span * 0.995)),
            Color::srgb(0.55, 0.95, 1.0),
        );
    }
    if !demo.brushes_shown {
        return;
    }

    for (brush, &role) in demo.fixture.hard.iter().zip(&demo.row(*shot).roles) {
        if !shot.tape.keeps(role) {
            continue;
        }
        let srgb = match role {
            Role::Necessary => NEEDED_SRGB,
            Role::DroppableFar => FAR_SRGB,
            Role::DroppableDominated => DOMINATED_SRGB,
            Role::Pruned => PRUNED_SRGB,
        };
        let colour = Color::srgb(srgb[0], srgb[1], srgb[2]);
        match brush.shape {
            Shape::Sphere(s) => {
                gizmos
                    .sphere(place(s.center), s.radius as f32, colour)
                    .resolution(BRUSH_SEGMENTS);
            }
            Shape::Capsule(c) => {
                let (a, b) = (place(c.a), place(c.b));
                gizmos
                    .sphere(a, c.radius as f32, colour)
                    .resolution(BRUSH_SEGMENTS);
                gizmos
                    .sphere(b, c.radius as f32, colour)
                    .resolution(BRUSH_SEGMENTS);
                gizmos.line(a, b, colour);
            }
        }
    }
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

/// `on`/`off`, for the HUD.
fn on_off(value: bool) -> &'static str {
    if value { "on" } else { "off" }
}

/// The HUD. The numbers are the demo.
fn report(
    demo: Res<Demo>,
    shot: Res<Shot>,
    tour: Res<Tour>,
    ledger: Res<Ledger>,
    cross: Res<CrossCheck>,
    mut stats: ResMut<DemoStats>,
    mut caption: Query<&mut Text, With<Caption>>,
) {
    let shot = *shot;
    let row = demo.row(shot);
    let world = &demo.world;
    stats.title =
        String::from("E-315 edit tape trim -- M-358 / P-59: the tape you keep is 20x too big");
    stats.vertices = demo.live_vertices;
    stats.triangles = demo.live_triangles;
    stats.extract_ms = demo.live_ms;

    let recorded = |tape: Tape, note: &str| {
        format!(
            "  {:<10} {:>3} brushes  hash {:<20} {note}",
            match tape {
                Tape::Full => "full",
                Tape::Survivors => "survivors",
                Tape::Necessary => "necessary",
            },
            row.count(tape),
            row.hash(tape)
        )
    };
    let identical = row.control_unchanged() && row.necessary_only_unchanged();

    let mut extra = vec![
        format!(
            "chunk      {:<6} {:>2} of {}   {} tris, {} verts, {}^3 samples, cell {CELL_SIZE}",
            row.label(),
            shot.chunk + 1,
            demo.rows.len(),
            commas(row.triangles),
            commas(row.vertices),
            CELLS_PER_CHUNK + 1
        ),
        format!(
            "tape       {:<10} {:>2} of {BRUSHES} brushes drawn{}{}",
            shot.tape.name(),
            row.count(shot.tape),
            if row.triangles == 0 {
                "   (empty chunk: necessary == 0 is free)"
            } else {
                ""
            },
            if tour.manual { "   [steering]" } else { "" }
        ),
        String::new(),
        recorded(Tape::Full, "the whole edit history"),
        recorded(Tape::Survivors, "P-39's Lipschitz bound, M-341"),
        recorded(Tape::Necessary, "leave-one-out"),
        format!(
            "  {:<10} {:>3} brushes  hash {:<20} extracted this frame",
            "on screen", demo.live_brushes, demo.live_hash
        ),
        format!(
            "  bit-identical: {}",
            if identical {
                "YES -- one mesh, four tapes, four equal hashes"
            } else {
                "NO -- the bound is unsound on this chunk; see the log"
            }
        ),
        String::new(),
        format!(
            "over-keep  {} of this chunk's {} survivors drop individually, and all {} at once",
            row.droppable(),
            row.survivors,
            row.droppable()
        ),
        format!(
            "cost       {} re-meshes here, {} for the world -- headroom, NOT a pruner",
            row.remeshes(),
            commas(world.total_remeshes)
        ),
        format!(
            "C3 here    {:>3} far (clear of zero by > 1 cell)   {:>3} dominated (straddle the band)",
            row.droppable_far, row.droppable_dominated
        ),
        format!(
            "C3 world   {} far of {} = {:.6} against the registered {:.2} -- FALSIFIED",
            commas(world.total_droppable_far),
            commas(world.total_droppable),
            world.far_fraction(),
            ledger.c3_bound
        ),
        format!(
            "           far_by_hi {} of {}: `hi < -{CELL_SIZE}` needs f(chunk centre) < {:.4}",
            world.total_far_by_hi,
            commas(world.total_droppable),
            world.hi_branch_threshold()
        ),
        format!(
            "           at reach {:.4}, and no brush's field goes below {:.4}. C3's other half",
            world.reach, -world.deepest_brush
        ),
        String::from("           is structurally dead on this fixture, and the CSV shows it."),
        format!(
            "world      {} survivors -> {} necessary, bit-exact on {}/{} -- {:.1}x, {} need none",
            commas(world.total_survivors),
            commas(world.total_necessary),
            world.necessary_only_unchanged,
            world.chunks,
            world.trim(),
            world.chunks_needing_none
        ),
        format!(
            "           {} dominant Adds here; ablation {:.2} s, one chunk/frame, {} empty chunks",
            row.dominant_adds, world.seconds, world.empty_chunks
        ),
        String::new(),
        String::from("p-59.csv   expected -> measured"),
    ];
    for pair in cross.0.chunks(2) {
        let mut line = format!("  {:<44}", pair[0].line());
        if let Some(second) = pair.get(1) {
            line.push_str(&second.line());
        }
        if pair.iter().any(|c| !c.holds()) {
            line.push_str("   MISMATCH");
        }
        extra.push(line);
    }
    extra.push(if cross.all_hold() {
        format!("  all {} reproduce", cross.0.len())
    } else {
        format!(
            "  CROSS-CHECK FAILED: {} of {} reproduce",
            cross.held(),
            cross.0.len()
        )
    });
    extra.push(String::new());
    extra.push(String::from(
        "brushes    cyan needed   amber droppable-far   magenta dominated   grey pruned",
    ));
    extra.push(format!(
        "[T] tape   [ [ ] and [ ] ] chunk   [B] brushes {}   [K] boxes {}   [X] tour",
        on_off(demo.brushes_shown),
        on_off(demo.boxes_shown)
    ));
    stats.extra = extra;

    let text = caption_for(row, shot, world);
    for mut target in &mut caption {
        target.0.clone_from(&text);
    }
}

/// The line a viewer reads instead of the HUD.
fn caption_for(row: &Ablation, shot: Shot, world: &World) -> String {
    match shot.tape {
        Tape::Full => format!(
            "chunk {}: {} folded out of the whole {BRUSHES}-brush edit history\n{}",
            row.label(),
            if row.triangles == 0 {
                String::from("no triangles at all")
            } else {
                format!("{} triangles", commas(row.triangles))
            },
            stop_caption(row)
        ),
        Tape::Survivors => format!(
            "P-39's Lipschitz bound keeps {} of {BRUSHES} -- the rock has not moved,\nand the hash is still {}",
            row.survivors, row.hash_survivors
        ),
        Tape::Necessary if row.necessary == 0 => format!(
            "leave-one-out needs NONE of them: the bare sphere meshes to the same\n{}. {} re-meshes to find that out.",
            if row.triangles == 0 {
                String::from("nothing")
            } else {
                format!("{} triangles, byte for byte", commas(row.triangles))
            },
            row.remeshes()
        ),
        Tape::Necessary => format!(
            "leave-one-out needs {} of the {} survivors -- same mesh, same hash.\nWorld-wide {} -> {} brushes, bit-exact, for {} re-meshes.",
            row.necessary,
            row.survivors,
            commas(world.total_survivors),
            commas(world.total_necessary),
            commas(world.total_remeshes)
        ),
    }
}

/// The one-line reason this chunk is a stop on the tour.
fn stop_caption(row: &Ablation) -> String {
    if row.triangles == 0 {
        return format!(
            "an empty chunk: {} survivors kept, and zero needed for free",
            row.survivors
        );
    }
    if row.necessary == 0 {
        return String::from("watch this one -- it needs none of them");
    }
    if row.droppable_dominated > row.droppable_far {
        return format!(
            "{} of its {} droppable survivors are dominated, not distant",
            row.droppable_dominated,
            row.droppable()
        );
    }
    format!(
        "{} survivors, {} needed, and {} of the rest are simply far from the surface",
        row.survivors, row.necessary, row.droppable_far
    )
}
