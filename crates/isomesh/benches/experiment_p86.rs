//! **P-86 — how many of a character's stops are geometry artefacts.**
//!
//! Ticket: R-086. Pre-registered before this harness existed.
//!
//! ```bash
//! cargo bench --bench experiment_p86
//! ```
//!
//! Writes `docs/experiments/p-86.csv`.
//!
//! # Hypothesis, as registered
//!
//! `M-115` found that a moving body is stopped harder and more often by ordinary
//! terrain than by a chunk join, and left the harder question open: **how many of
//! those stops are real.** The crate records degenerate near-zero-area triangles
//! as a metric rather than a gate, because Marching Cubes genuinely emits slivers
//! whenever a grid corner sits near zero. But a sliver is exactly what a capsule
//! controller catches on, and `M-185` found that completing the crossing identity
//! turned a sliver into a repeated-index triangle the extractor now declines to
//! emit. Nobody has connected the recorded metric to the gameplay symptom.
//!
//! * **C1** — over `game_capsule_walk`'s 495 seam crossings plus a 10⁴-step
//!   randomised walk on `fbm_terrain`, at least **20%** of controller stops occur
//!   on a triangle in the bottom decile of the aspect-ratio distribution.
//! * **C2** — stops are **not** concentrated at chunk seams: seam-adjacent
//!   triangles carry no more stops per triangle than interior ones.
//! * **C3** — the field-query path removes **all** of C1's stops, because it
//!   never sees a triangle.
//!
//! **Falsified by** C1 under 20% (slivers are not what stops a character, and the
//! recorded metric has no gameplay consequence); C2 by a seam excess (which
//! reopens `M-133` as a gameplay defect rather than a topology one); C3 by any
//! surviving stop (the terrain genuinely stops the character there and C1 was
//! measuring geometry that is correct).
//!
//! **Vacuity control, as registered:** the walk must produce a non-zero stop
//! count on both arms, and the bottom aspect-ratio decile must be non-empty, both
//! reported as counts. Both are `assert!`ed here, and two more are added below
//! because the registered pair does not cover C3.
//!
//! # SHARE — recomputed before the harness was written
//!
//! The registration says *"this is a rate, not a ratio of a total"*, and that is
//! right as far as it goes: `worst_decile_fraction` is `stops_on_worst_decile /
//! stops_total`, both of which this harness produces, so 20% is reachable in the
//! trivial sense that the numerator can be anything from 0 to the denominator.
//!
//! The share that matters is the **null rate**, and the registration does not
//! name it. A decile is 10% of the triangles by construction, so if a stop picked
//! a triangle uniformly at random C1's bar would be exactly **2× chance**. But a
//! stop does not pick uniformly: a swept capsule hits a triangle roughly in
//! proportion to the area it presents, and a bottom-decile triangle is
//! *near-zero-area by definition*. So the honest null is the decile's **area
//! share**, which this harness emits as `worst_decile_area_share`, and C1's bar is
//! `0.20 / worst_decile_area_share` times chance. That number is the whole
//! mechanism of the result and it is computed rather than assumed.
//!
//! # What was lifted, and how it was checked
//!
//! `CLAUDE.md` forbids Bevy under `crates/`, and `game_capsule_walk` is an Avian
//! rigid body. So neither controller is *imported*; both are re-implemented here
//! from the shipped arithmetic, and each is checked against something the shipped
//! path is already known to do:
//!
//! * **The mesh controller** is collide-and-slide over `parry3d` shape casts — a
//!   swept capsule, up to [`SLIDE_ITERATIONS`] plane projections, then a swept
//!   vertical move. `game_capsule_walk` uses a *dynamic* Avian capsule with
//!   `Friction::new(0.4)`, no step offset and no slope limit; this has no friction
//!   and no step offset either, which is the direction that matters — a step
//!   offset would climb over exactly the lips this experiment is about. The stall
//!   definition is `measure`'s verbatim: a step whose horizontal progress falls
//!   short of what was commanded by at least [`STALL_SHORTFALL`].
//! * **The field controller** is `game_dig::resolve_body` transcribed — deepest
//!   overlap of a sphere chain, push along the normalised gradient,
//!   [`RESOLVE_PASSES`] passes, cancel only the velocity component *into* the
//!   surface, and the `✗45`/`M-363` **chord** ground probe. Checked by
//!   reproducing `✗45` itself: [`ground_probe_lift_check`] settles the same body
//!   on the same field twice, once with the chord probe and once with the flat
//!   cross the shipped code used to have, and asserts the flat cross hovers
//!   *higher*. A lift that had quietly picked up the pre-fix geometry cannot pass
//!   that.
//!
//! # Why the field arm is a per-stop replay and not a second walk
//!
//! C3 asks whether the field path removes **C1's** stops — the specific ones. Two
//! independent walks diverge within a few steps and their stop *counts* are then
//! two different questions, so `stops_field_path` is a replay: at each stop the
//! mesh controller records `(position, vertical velocity, commanded step)`, and
//! the field controller is handed that exact state and asked whether it also
//! falls short. An independent field-only walk over the same route runs anyway and
//! is reported as `stops_field_walk`, because it is the number a reader will
//! expect to see.

#![allow(
    clippy::cast_possible_truncation,
    clippy::cast_precision_loss,
    clippy::cast_sign_loss,
    clippy::print_literal,
    clippy::too_many_lines
)]

mod common;

use std::collections::{HashMap, HashSet};
use std::time::Instant;

use isomesh::chunk::ChunkLayout;
use isomesh::extractor::Extractor;
use isomesh::fields::FbmTerrain;
use isomesh::marching_cubes::MarchingCubes;
use isomesh::{MeshBuffer, Sdf};
use parry3d::math::{Pose, Vector};
use parry3d::query::{ShapeCastOptions, cast_shapes, contact};
use parry3d::shape::{Capsule, Triangle};

// ── the world, identical to E-203 and E-206b so the numbers are comparable ───

/// Cells per chunk axis. `game_capsule_walk`'s value.
const CHUNK_CELLS: u32 = 16;
/// World units per cell. `game_capsule_walk`'s value.
const CELL_SIZE: f32 = 0.5;
/// World units per chunk axis.
const CHUNK_SPAN: f32 = CHUNK_CELLS as f32 * CELL_SIZE;
/// The two vertical chunk layers that can hold `fbm_terrain`'s sheet: its height
/// bound is `amplitude * fbm_bound(4, 0.5)` = `2 * 1.875` = 3.75, and the layer
/// boundary is `y = 0`, so the surface genuinely straddles both.
const VERTICAL_LAYERS: [i32; 2] = [-1, 0];

/// How close to a chunk boundary counts as "at a seam", in world units. One cell,
/// which is `game_capsule_walk`'s `SEAM_BAND`.
const SEAM_BAND: f32 = CELL_SIZE;

// ── the body, identical to E-206b's capsule ──────────────────────────────────

/// Capsule radius. `game_capsule_walk`'s value.
const CAPSULE_RADIUS: f32 = 0.4;
/// Capsule segment length (so total height is `0.9 + 2 * 0.4` = 1.7, which is
/// also `game_dig`'s body height).
const CAPSULE_LENGTH: f32 = 0.9;
/// Half the segment.
const CAPSULE_HALF: f32 = CAPSULE_LENGTH * 0.5;

/// Commanded speed. `game_capsule_walk`'s value.
const WALK_SPEED: f32 = 7.0;
/// Fixed step. `game_capsule_walk` runs on real frame times; a fixed step is what
/// makes this reproducible, and is a deviation recorded in the report.
const DT: f32 = 1.0 / 60.0;
/// `game_dig`'s gravity.
const GRAVITY: f32 = 18.0;
/// Fall speed cap, so the swept AABB a step needs has a bound. Reached only if
/// the body leaves the world, which the fixture control forbids.
const MAX_FALL: f32 = 30.0;

/// A step whose horizontal progress falls this far short of what was commanded is
/// a stall. `game_capsule_walk::measure`'s threshold, verbatim.
const STALL_SHORTFALL: f32 = 0.5;

/// Plane projections per swept move. Four is the usual collide-and-slide budget:
/// enough for a corner of three planes plus one spare.
const SLIDE_ITERATIONS: u32 = 4;
/// Back-off along the motion at an impact, so the next cast does not start
/// touching. One thousandth of a cell.
const SKIN: f32 = 5.0e-4;

// ── the field controller, transcribed from game_dig ──────────────────────────

/// Spheres in the field controller's body, spaced along the capsule segment.
///
/// `game_dig` uses four of radius 0.25 on a 1.7-unit column. Seven here, because
/// the union of spheres **under**-approximates the capsule between their centres
/// by `R - sqrt(R² - (spacing/2)²)`, and under-approximating is the direction that
/// would quietly help C3: a thinner body is stopped less. Seven gives a spacing of
/// 0.15 and an inset of 0.0073, which is an eighth of [`GROUND_PROBE`].
const BODY_SPHERES: usize = 7;
/// `game_dig::RESOLVE_PASSES`.
const RESOLVE_PASSES: u32 = 3;
/// `game_dig::GROUND_PROBE`.
const GROUND_PROBE: f32 = 0.06;
/// `game_dig::GRADIENT_EPS`.
const GRADIENT_EPS: f32 = 1.0e-4;
/// `resolve_body`'s footprint fraction for the chord cross.
const FOOTPRINT: f32 = 0.7;

// ── the two fixtures ─────────────────────────────────────────────────────────

/// Seam crossings the path arm runs for. `M-106`'s figure, and the registration's.
const TARGET_CROSSINGS: u64 = 495;
/// Hard cap on the path arm, so a body that stops crossing seams fails loudly
/// rather than running for ever.
const MAX_PATH_STEPS: u64 = 200_000;
/// Measured steps in the randomised arm. The registration's 10⁴.
const RANDOM_STEPS: u64 = 10_000;
/// Side of the randomised arm's square, in world units — eight chunks, so the
/// walk crosses seams inside a region small enough to keep every triangle it
/// could touch in memory at once.
const RANDOM_EXTENT: f32 = 64.0;
/// Steps discarded before measuring, so a body still falling from its spawn is
/// not counted as blocked. `game_capsule_walk::SETTLE_SECONDS` is 0.75 s, which
/// at [`DT`] is 45 steps; 120 is a generous round number.
const SETTLE_STEPS: u64 = 120;
/// Spawn height. Above `fbm_terrain`'s 3.75 height bound and no higher.
const SPAWN_HEIGHT: f32 = 6.0;

/// Broadphase bucket side, in world units. Two cells: a Marching Cubes triangle
/// spans at most one cell, so a triangle lands in one or two buckets.
const BUCKET: f32 = 1.0;
/// How far from the body a chunk must be to stay unmeshed. Capsule radius plus a
/// step plus the ground probe, rounded up hard.
const REACH: f32 = 2.0;

/// `2√3`, from `validate_indexed`'s mean-ratio pass.
const TWO_ROOT_THREE: f64 = 3.464_101_615_137_754_6;

// ── triangles ────────────────────────────────────────────────────────────────

/// One world triangle, with everything the attribution needs.
#[derive(Clone, Copy)]
struct Tri {
    /// The three corners, in world units.
    v: [[f32; 3]; 3],
    /// Axis-aligned bounds, for the broadphase reject.
    lo: [f32; 3],
    /// Axis-aligned bounds, for the broadphase reject.
    hi: [f32; 3],
    /// Mean-ratio quality `q = 2√3·|cross| / Σlᵢ²`, transcribed from
    /// `validate_indexed`: **1 for equilateral, 0 for degenerate**. This is the
    /// crate's own recorded shape metric, so "bottom decile of the aspect-ratio
    /// distribution" is the bottom decile of `q`.
    q: f32,
    /// A corner within [`SEAM_BAND`] of a chunk boundary plane on any axis.
    seam: bool,
    /// The same on `x` and `z` only, which is `game_capsule_walk`'s definition.
    seam_xz: bool,
    /// The same on `y` only — the horizontal chunk plane at `y = 0`.
    seam_y: bool,
    /// `area <= ValidateConfig::AREA_EPSILON_REL · cell_size²` — the crate's
    /// `degenerate_triangles` metric, the one this experiment exists to connect to
    /// a gameplay symptom.
    degenerate: bool,
    /// `|∇h|` at the centroid, from the field's own analytic gradient.
    ///
    /// The confound instrument. Stops happen where the terrain is steep, so a
    /// bucket of triangles that is steeper than its complement will carry more
    /// stops for a reason that is nothing to do with what the bucket is named
    /// after — and the seam band on `y` turns out to be exactly such a bucket.
    slope: f32,
}

/// What survives an arm: the population a decile is taken over.
#[derive(Clone, Copy)]
struct TriMeta {
    /// Mean-ratio quality.
    q: f32,
    /// Twice the area, which is what an area share is weighted by.
    two_area: f32,
    /// Seam-adjacent on any axis.
    seam: bool,
    /// Seam-adjacent on `x` or `z`.
    seam_xz: bool,
    /// Seam-adjacent on `y`.
    seam_y: bool,
    /// Numerically zero area.
    degenerate: bool,
    /// `|∇h|` at the centroid.
    slope: f32,
    /// The controller actually tested this triangle at least once.
    ///
    /// The honest denominator for C2. `ensure` meshes a whole chunk when the body
    /// comes within [`REACH`] of it, so a chunk the walk merely clipped
    /// contributes ~500 triangles to the denominator of which only the ones near
    /// its own boundary were ever reachable — which manufactures a seam excess
    /// out of the meshing policy. A triangle the broadphase handed to the
    /// controller is one the capsule could have hit.
    tested: bool,
}

/// The meshed world, grown as the body walks into it.
struct Terrain {
    /// The field being meshed.
    field: FbmTerrain<f32>,
    /// The chunk partition.
    layout: ChunkLayout<f32>,
    /// The extractor. `game_capsule_walk` uses `Extractor::MarchingCubes`.
    mc: MarchingCubes<f32>,
    /// Every triangle meshed so far.
    tris: Vec<Tri>,
    /// Triangle indices by broadphase bucket.
    buckets: HashMap<[i32; 3], Vec<u32>>,
    /// Chunks already meshed, so a revisited column is not counted twice.
    meshed: HashSet<[i32; 3]>,
    /// Scratch for the mesh sink.
    out: MeshBuffer<f32>,
    /// Scratch for candidate gathering.
    candidates: Vec<u32>,
    /// Per-triangle stamp, so candidate gathering deduplicates without sorting.
    stamp: Vec<u32>,
    /// Current stamp generation.
    generation: u32,
    /// Whether each triangle was ever handed to the controller.
    tested: Vec<bool>,
}

impl Terrain {
    /// An empty world over `field`.
    fn new(field: FbmTerrain<f32>) -> Self {
        Self {
            field,
            layout: ChunkLayout::new(CHUNK_CELLS, CELL_SIZE, [0.0; 3]).expect("a valid layout"),
            mc: MarchingCubes::new(),
            tris: Vec::new(),
            buckets: HashMap::new(),
            meshed: HashSet::new(),
            out: MeshBuffer::new(),
            candidates: Vec::new(),
            stamp: Vec::new(),
            generation: 0,
            tested: Vec::new(),
        }
    }

    /// Mesh every chunk the body at `(x, z)` could touch.
    ///
    /// The 3×3 column neighbourhood filtered by real distance to the column's
    /// box, rather than the neighbourhood entire: the path arm walks 2.3 km and
    /// the unfiltered swath is three columns wide the whole way, which is three
    /// times the triangles for geometry the capsule is never within a metre of.
    fn ensure(&mut self, x: f32, z: f32) {
        let here = self.layout.chunk_of([x, 0.0, z]);
        for dx in -1..=1 {
            for dz in -1..=1 {
                let cx = here.coords[0] + dx;
                let cz = here.coords[2] + dz;
                let near = |v: f32, c: i32| {
                    let lo = c as f32 * CHUNK_SPAN;
                    let hi = lo + CHUNK_SPAN;
                    (lo - v).max(v - hi).max(0.0)
                };
                let gap = near(x, cx).hypot(near(z, cz));
                if gap > REACH {
                    continue;
                }
                for layer in VERTICAL_LAYERS {
                    if self.meshed.insert([cx, layer, cz]) {
                        self.mesh_chunk([cx, layer, cz]);
                    }
                }
            }
        }
    }

    /// Mesh one chunk and file its triangles.
    fn mesh_chunk(&mut self, coords: [i32; 3]) {
        let id = isomesh::chunk::ChunkId { coords };
        let shape = self.layout.sample_shape().expect("a valid sample shape");
        let origin = self.layout.sample_origin(id);
        self.out.reset();
        self.mc
            .extract_into(&self.field, &shape, origin, CELL_SIZE, &mut self.out)
            .expect("marching cubes over a valid chunk");

        // `area <= AREA_EPSILON_REL · cell_size²`, compared on `2A = |cross|` and
        // squared, exactly as `validate_indexed` does it.
        let two_area_limit = 2.0 * 1.0e-6 * f64::from(CELL_SIZE) * f64::from(CELL_SIZE);
        let limit_sq = two_area_limit * two_area_limit;

        for f in self.out.indices.as_chunks::<3>().0 {
            let v = [
                self.out.positions[f[0] as usize],
                self.out.positions[f[1] as usize],
                self.out.positions[f[2] as usize],
            ];
            let a = v[0].map(f64::from);
            let b = v[1].map(f64::from);
            let c = v[2].map(f64::from);
            let u = [b[0] - a[0], b[1] - a[1], b[2] - a[2]];
            let w = [c[0] - a[0], c[1] - a[1], c[2] - a[2]];
            let cross = [
                u[1] * w[2] - u[2] * w[1],
                u[2] * w[0] - u[0] * w[2],
                u[0] * w[1] - u[1] * w[0],
            ];
            let len_sq = cross[0] * cross[0] + cross[1] * cross[1] + cross[2] * cross[2];
            let e = [w[0] - u[0], w[1] - u[1], w[2] - u[2]];
            let edge_sq = u[0] * u[0]
                + u[1] * u[1]
                + u[2] * u[2]
                + (w[0] * w[0] + w[1] * w[1] + w[2] * w[2])
                + (e[0] * e[0] + e[1] * e[1] + e[2] * e[2]);
            let q = if edge_sq > 0.0 {
                TWO_ROOT_THREE * len_sq.sqrt() / edge_sq
            } else {
                0.0
            };

            let mut lo = v[0];
            let mut hi = v[0];
            for p in &v[1..] {
                for axis in 0..3 {
                    lo[axis] = lo[axis].min(p[axis]);
                    hi[axis] = hi[axis].max(p[axis]);
                }
            }

            // Seam adjacency is a proximity band, the same test
            // `game_capsule_walk::measure` applies to the body: distance to the
            // nearest chunk plane on an axis, under one cell.
            let to_plane = |value: f32| {
                let plane = (value / CHUNK_SPAN).round() * CHUNK_SPAN;
                (value - plane).abs()
            };
            let mut seam = false;
            let mut seam_xz = false;
            let mut seam_y = false;
            for p in &v {
                for (axis, coordinate) in p.iter().enumerate() {
                    if to_plane(*coordinate) < SEAM_BAND {
                        seam = true;
                        if axis == 1 {
                            seam_y = true;
                        } else {
                            seam_xz = true;
                        }
                    }
                }
            }

            let centroid = [0, 1, 2].map(|axis| (v[0][axis] + v[1][axis] + v[2][axis]) / 3.0);
            let g = self.field.gradient(centroid);

            let index = self.tris.len() as u32;
            self.tris.push(Tri {
                v,
                lo,
                hi,
                q: q as f32,
                seam,
                seam_xz,
                seam_y,
                degenerate: len_sq <= limit_sq,
                slope: g[0].hypot(g[2]),
            });
            let key_lo = [0, 1, 2].map(|axis| (lo[axis] / BUCKET).floor() as i32);
            let key_hi = [0, 1, 2].map(|axis| (hi[axis] / BUCKET).floor() as i32);
            for bx in key_lo[0]..=key_hi[0] {
                for by in key_lo[1]..=key_hi[1] {
                    for bz in key_lo[2]..=key_hi[2] {
                        self.buckets.entry([bx, by, bz]).or_default().push(index);
                    }
                }
            }
        }
        self.stamp.resize(self.tris.len(), 0);
        self.tested.resize(self.tris.len(), false);
    }

    /// Gather every triangle whose bounds meet the box `lo..hi`, into
    /// [`Self::candidates`].
    fn gather(&mut self, lo: [f32; 3], hi: [f32; 3]) {
        self.candidates.clear();
        self.generation += 1;
        let generation = self.generation;
        let key_lo = [0, 1, 2].map(|axis| (lo[axis] / BUCKET).floor() as i32);
        let key_hi = [0, 1, 2].map(|axis| (hi[axis] / BUCKET).floor() as i32);
        for bx in key_lo[0]..=key_hi[0] {
            for by in key_lo[1]..=key_hi[1] {
                for bz in key_lo[2]..=key_hi[2] {
                    let Some(bucket) = self.buckets.get(&[bx, by, bz]) else {
                        continue;
                    };
                    for &index in bucket {
                        if self.stamp[index as usize] == generation {
                            continue;
                        }
                        self.stamp[index as usize] = generation;
                        let tri = &self.tris[index as usize];
                        if (0..3).all(|axis| tri.lo[axis] <= hi[axis] && tri.hi[axis] >= lo[axis]) {
                            self.tested[index as usize] = true;
                            self.candidates.push(index);
                        }
                    }
                }
            }
        }
    }

    /// The population an arm's decile is taken over.
    fn population(&self) -> Vec<TriMeta> {
        self.tris
            .iter()
            .enumerate()
            .map(|(index, t)| {
                let a = t.v[0].map(f64::from);
                let b = t.v[1].map(f64::from);
                let c = t.v[2].map(f64::from);
                let u = [b[0] - a[0], b[1] - a[1], b[2] - a[2]];
                let w = [c[0] - a[0], c[1] - a[1], c[2] - a[2]];
                let cross = [
                    u[1] * w[2] - u[2] * w[1],
                    u[2] * w[0] - u[0] * w[2],
                    u[0] * w[1] - u[1] * w[0],
                ];
                let two_area =
                    (cross[0] * cross[0] + cross[1] * cross[1] + cross[2] * cross[2]).sqrt();
                TriMeta {
                    q: t.q,
                    two_area: two_area as f32,
                    seam: t.seam,
                    seam_xz: t.seam_xz,
                    seam_y: t.seam_y,
                    degenerate: t.degenerate,
                    slope: t.slope,
                    tested: self.tested[index],
                }
            })
            .collect()
    }
}

// ── the mesh controller ──────────────────────────────────────────────────────

/// A swept capsule cast against one triangle, in the world frame.
fn cast(pos: Vector, motion: Vector, tri: &Tri) -> Option<(f32, Vector)> {
    let capsule = Capsule::new_y(CAPSULE_HALF, CAPSULE_RADIUS);
    let target = Triangle::new(
        Vector::from_array(tri.v[0]),
        Vector::from_array(tri.v[1]),
        Vector::from_array(tri.v[2]),
    );
    let hit = cast_shapes(
        &Pose::from_translation(pos),
        motion,
        &capsule,
        &Pose::IDENTITY,
        Vector::ZERO,
        &target,
        ShapeCastOptions {
            max_time_of_impact: 1.0,
            target_distance: 0.0,
            // A time-zero impact whose relative velocity is *separating* is a
            // body resting in contact and walking away from it, which is not a
            // block. Left `true` the controller would report a stop on every
            // frame it stood on the ground.
            stop_at_penetration: false,
            compute_impact_geometry_on_penetration: true,
        },
    )
    .expect("parry supports capsule-vs-triangle shape casting");
    hit.map(|h| (h.time_of_impact, h.normal1))
}

/// Collide-and-slide one swept move, returning what it achieved and the
/// triangles that blocked it, earliest first.
fn slide(terrain: &Terrain, pos: Vector, motion: Vector, blockers: &mut Vec<u32>) -> Vector {
    let mut moved = Vector::ZERO;
    let mut remaining = motion;
    for _ in 0..SLIDE_ITERATIONS {
        let length = remaining.length();
        if length < 1.0e-7 {
            break;
        }
        let at = pos + moved;
        let mut best: Option<(f32, Vector, u32)> = None;
        for &index in &terrain.candidates {
            let tri = &terrain.tris[index as usize];
            if let Some((toi, normal)) = cast(at, remaining, tri)
                && best.is_none_or(|(t, _, _)| toi < t)
            {
                best = Some((toi, normal, index));
            }
        }
        match best {
            None => {
                moved += remaining;
                break;
            }
            Some((toi, normal, index)) => {
                blockers.push(index);
                let back = (SKIN / length).min(toi);
                let advance = remaining * (toi - back).max(0.0);
                moved += advance;
                remaining -= advance;
                // The plane projection. `normal1` is the outward normal on the
                // capsule at the impact, so it points *into* the triangle, and
                // removing the component along it is the slide.
                remaining -= normal * remaining.dot(normal);
            }
        }
    }
    moved
}

/// Push the capsule out of any triangle it starts a step inside.
///
/// Without this a body that lands a hair inside a face stays inside it, every
/// later cast returns a time-zero impact, and the stop count becomes a count of
/// one wedged frame repeated — a fixture defect that would look exactly like a
/// finding. Two passes, deepest contact first, which is `resolve_body`'s rule in
/// the mesh controller's costume.
fn depenetrate(terrain: &Terrain, pos: &mut Vector) {
    let capsule = Capsule::new_y(CAPSULE_HALF, CAPSULE_RADIUS);
    for _ in 0..2 {
        let mut deepest: Option<(f32, Vector)> = None;
        for &index in &terrain.candidates {
            let tri = &terrain.tris[index as usize];
            let target = Triangle::new(
                Vector::from_array(tri.v[0]),
                Vector::from_array(tri.v[1]),
                Vector::from_array(tri.v[2]),
            );
            if let Ok(Some(c)) = contact(
                &Pose::from_translation(*pos),
                &capsule,
                &Pose::IDENTITY,
                &target,
                0.0,
            ) && c.dist < 0.0
                && deepest.is_none_or(|(d, _)| c.dist < d)
            {
                deepest = Some((c.dist, c.normal1));
            }
        }
        let Some((dist, normal)) = deepest else {
            return;
        };
        *pos -= normal * (dist.abs() + SKIN);
    }
}

/// The mesh controller's state.
struct MeshBody {
    /// Capsule centre.
    pos: Vector,
    /// Vertical velocity. Horizontal motion is commanded rather than integrated,
    /// which is `game_dig`'s shape and the reason a resting body does not slide.
    vy: f32,
    /// Whether the last vertical move was blocked from below.
    grounded: bool,
}

/// What one step of a controller did.
struct Outcome {
    /// Commanded horizontal distance.
    asked: f32,
    /// Achieved horizontal distance.
    moved: f32,
    /// Fraction of the commanded distance not covered, clamped to `0..=1`.
    shortfall: f32,
}

impl MeshBody {
    /// One step: depenetrate, slide horizontally, then fall.
    fn step(
        &mut self,
        terrain: &mut Terrain,
        commanded: Vector,
        blockers: &mut Vec<u32>,
    ) -> Outcome {
        blockers.clear();
        self.vy = (self.vy - GRAVITY * DT).max(-MAX_FALL);
        let fall = Vector::new(0.0, self.vy * DT, 0.0);

        // One gather per step, over a box covering both moves plus the capsule.
        // Correct because the box is the union of every position the step visits,
        // expanded by the skin.
        let span = commanded.length() + fall.length() + CAPSULE_RADIUS + GROUND_PROBE + SKIN;
        let lo = [
            self.pos.x - span,
            self.pos.y - CAPSULE_HALF - span,
            self.pos.z - span,
        ];
        let hi = [
            self.pos.x + span,
            self.pos.y + CAPSULE_HALF + span,
            self.pos.z + span,
        ];
        terrain.gather(lo, hi);

        depenetrate(terrain, &mut self.pos);

        let before = self.pos;
        let horizontal = slide(terrain, self.pos, commanded, blockers);
        self.pos += horizontal;

        let mut vertical_blockers = Vec::new();
        let dropped = slide(terrain, self.pos, fall, &mut vertical_blockers);
        self.pos += dropped;
        self.grounded = !vertical_blockers.is_empty() && self.vy < 0.0;
        if self.grounded {
            self.vy = 0.0;
        }

        let asked = commanded.length();
        let moved = Vector::new(self.pos.x - before.x, 0.0, self.pos.z - before.z).length();
        Outcome {
            asked,
            moved,
            shortfall: ((asked - moved) / asked).clamp(0.0, 1.0),
        }
    }
}

// ── the field controller, transcribed from game_dig::resolve_body ────────────

/// The sphere-chain offsets from the capsule centre, in `-Y`.
fn body_offsets() -> [f32; BODY_SPHERES] {
    let mut offsets = [0.0; BODY_SPHERES];
    for (i, o) in offsets.iter_mut().enumerate() {
        let t = i as f32 / (BODY_SPHERES - 1) as f32;
        *o = -CAPSULE_HALF + t * CAPSULE_LENGTH;
    }
    offsets
}

/// Which ground probe geometry to use. Two arms of one lift check, not two paths
/// in the controller: [`GroundProbe::FlatCross`] is the geometry `✗45` removed and
/// exists here only so [`ground_probe_lift_check`] can prove the transcription
/// took the surviving one.
#[derive(Clone, Copy, PartialEq, Eq)]
enum GroundProbe {
    /// `M-363`'s chord: lateral samples on the foot sphere's own lower surface.
    Chord,
    /// `✗45`'s five samples at one depth, which hovers on a slope.
    FlatCross,
}

/// Push the body out of the rock, and report whether it is standing on it.
///
/// `game_dig::resolve_body`, transcribed: deepest overlap over the sphere chain,
/// [`RESOLVE_PASSES`] passes, push along the normalised gradient, `Vec3::Y` below
/// [`GRADIENT_EPS`] because `M-172` measured an exactly zero gradient on the
/// medial axis, cancel only the velocity component *into* the surface, and a
/// separate downward probe for `grounded` because a body pressed against a wall is
/// in contact and is not standing on anything.
fn resolve_field(
    field: &impl Sdf<Scalar = f32>,
    centre: &mut Vector,
    velocity: &mut Vector,
    probe: GroundProbe,
) -> bool {
    let offsets = body_offsets();
    let mut deepest: Option<(f32, Vector)> = None;
    for _ in 0..RESOLVE_PASSES {
        for offset in offsets {
            let c = *centre + Vector::Y * offset;
            let f = field.sample([c.x, c.y, c.z]);
            if f >= CAPSULE_RADIUS {
                continue;
            }
            let g = Vector::from_array(field.gradient([c.x, c.y, c.z]));
            let n = if g.length() > GRADIENT_EPS {
                g.normalize()
            } else {
                Vector::Y
            };
            let depth = CAPSULE_RADIUS - f;
            *centre += n * depth;
            if deepest.is_none_or(|(d, _)| depth > d) {
                deepest = Some((depth, n));
            }
        }
    }
    if let Some((_, n)) = deepest {
        let into = velocity.dot(n);
        if into < 0.0 {
            *velocity -= n * into;
        }
    }

    let foot = *centre - Vector::Y * CAPSULE_HALF;
    let r = CAPSULE_RADIUS * FOOTPRINT;
    let edge = match probe {
        GroundProbe::Chord => CAPSULE_RADIUS * (1.0 - FOOTPRINT * FOOTPRINT).sqrt() + GROUND_PROBE,
        GroundProbe::FlatCross => CAPSULE_RADIUS + GROUND_PROBE,
    };
    let middle = CAPSULE_RADIUS + GROUND_PROBE;
    [
        [0.0, -middle, 0.0],
        [r, -edge, 0.0],
        [-r, -edge, 0.0],
        [0.0, -edge, r],
        [0.0, -edge, -r],
    ]
    .into_iter()
    .any(|[dx, dy, dz]| field.sample([foot.x + dx, foot.y + dy, foot.z + dz]) <= 0.0)
}

/// The field controller's state.
struct FieldBody {
    /// Capsule centre.
    pos: Vector,
    /// Vertical velocity.
    vy: f32,
    /// Whether the ground probe answered.
    grounded: bool,
}

impl FieldBody {
    /// One step. The horizontal command is a direct write and the fall is
    /// integrated, then the field pushes the body out — `game_dig::move_camera`'s
    /// order, which is the shipped field path's own shape.
    fn step(&mut self, field: &impl Sdf<Scalar = f32>, commanded: Vector) -> Outcome {
        // `gravity_step`: gravity is not integrated while standing, or a resting
        // body sinks `g·dt²` a frame and is pushed along the surface normal for
        // ever.
        if self.grounded && self.vy <= 0.0 {
            self.vy = 0.0;
        } else {
            self.vy = (self.vy - GRAVITY * DT).max(-MAX_FALL);
        }
        let before = self.pos;
        self.pos += commanded + Vector::new(0.0, self.vy * DT, 0.0);
        let mut velocity = Vector::new(0.0, self.vy, 0.0);
        self.grounded = resolve_field(field, &mut self.pos, &mut velocity, GroundProbe::Chord);
        self.vy = velocity.y;

        let asked = commanded.length();
        let moved = Vector::new(self.pos.x - before.x, 0.0, self.pos.z - before.z).length();
        Outcome {
            asked,
            moved,
            shortfall: ((asked - moved) / asked).clamp(0.0, 1.0),
        }
    }
}

/// Replay one mesh-controller stop against the field, from the same state.
///
/// Returns the shortfall and where the field controller put the body, so the
/// penetration probe can ask the mesh what that position means.
fn field_replay(
    field: &impl Sdf<Scalar = f32>,
    pos: Vector,
    vy: f32,
    grounded: bool,
    commanded: Vector,
) -> (f32, Vector) {
    let mut body = FieldBody { pos, vy, grounded };
    let outcome = body.step(field, commanded);
    (outcome.shortfall, body.pos)
}

/// The deepest overlap between the capsule at `pos` and the gathered triangles.
///
/// Zero when the capsule is clear. Reads only [`Terrain::candidates`], so it
/// costs one contact query per candidate and asks nothing new of the broadphase.
fn penetration(terrain: &Terrain, pos: Vector) -> f32 {
    let capsule = Capsule::new_y(CAPSULE_HALF, CAPSULE_RADIUS);
    let mut deepest = 0.0_f32;
    for &index in &terrain.candidates {
        let tri = &terrain.tris[index as usize];
        let target = Triangle::new(
            Vector::from_array(tri.v[0]),
            Vector::from_array(tri.v[1]),
            Vector::from_array(tri.v[2]),
        );
        if let Ok(Some(c)) = contact(
            &Pose::from_translation(pos),
            &capsule,
            &Pose::IDENTITY,
            &target,
            0.0,
        ) && c.dist < 0.0
        {
            deepest = deepest.max(-c.dist);
        }
    }
    deepest
}

// ── the walks ────────────────────────────────────────────────────────────────

/// `game_capsule_walk::path`: where the capsule is asked to be, horizontally,
/// after `d` metres. A diagonal with a wobble, so the walk crosses seams on both
/// axes and meets three-chunk corners.
fn path(d: f32) -> (f32, f32) {
    (d * 0.82, d * 0.57 + 9.0 * (d * 0.05).sin())
}

/// xorshift64. Deterministic, so the randomised arm is a fixture rather than a
/// die roll.
struct Rng(u64);

impl Rng {
    /// Next value in `0.0..1.0`.
    fn next(&mut self) -> f32 {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        ((self.0 >> 40) as f32) / 16_777_216.0
    }
}

/// One recorded stop.
struct Stop {
    /// Quality of the triangle the body ran into first.
    q: f32,
    /// Quality of the most sliver-like of every triangle that blocked the step —
    /// the generous reading of "a stop on a bottom-decile triangle".
    q_worst_blocker: f32,
    /// The first blocker is seam-adjacent on any axis.
    seam: bool,
    /// The first blocker is seam-adjacent on `x` or `z`.
    seam_xz: bool,
    /// The first blocker is one of the crate's `degenerate_triangles`.
    degenerate: bool,
    /// The field controller, handed this stop's own pre-step state, also fell
    /// short.
    field_stalls: bool,
    /// `|∇h|` at the blocking triangle. The confound instrument again: if a
    /// bucket's stops are steep, the terrain is what stopped the body.
    slope: f32,
    /// How far the **field** controller's replayed body ends up inside the
    /// triangle mesh, in world units.
    ///
    /// C3's mechanism, measured rather than argued. `f = y − h` is a *vertical*
    /// gap, so on a slope it overstates the true distance by `1/cos θ` and the
    /// push-out lets the body sit `R·(1 − cos θ)` inside the rock. If the field
    /// path "removes" a stop by standing in the hillside, the number that says so
    /// is this one.
    field_penetration: f32,
    /// The same for the **mesh** controller's own post-step position, which
    /// depenetrates and should therefore be near zero. The control that proves
    /// [`Stop::field_penetration`] is measuring penetration and not the probe.
    mesh_penetration: f32,
    /// The previous step was not a stop.
    ///
    /// **The effective sample size, and it is much smaller than the stop count.** A
    /// body climbing a 60° face stalls for every frame of the climb, so
    /// consecutive stops are the same event and Poisson error bars over the raw
    /// count would be several times too tight. Every rate is also reported over
    /// episodes.
    first_of_episode: bool,
}

/// One arm's result.
struct Arm {
    /// Measured steps.
    steps: u64,
    /// Steps in which the body's chunk column changed.
    seam_crossings: u64,
    /// Horizontal distance covered, world units.
    metres: f32,
    /// Every triangle the body could have touched.
    population: Vec<TriMeta>,
    /// Every stop.
    stops: Vec<Stop>,
    /// Stops of an independent field-only walk over the same route.
    stops_field_walk: u64,
    /// Steps in which a stop was recorded with no blocking triangle. Must be
    /// zero: a stop with nothing to attribute it to is an attribution failure,
    /// not a null.
    unattributed: u64,
}

/// Which fixture an arm runs.
#[derive(Clone, Copy, PartialEq, Eq)]
enum Fixture {
    /// `game_capsule_walk`'s path, to [`TARGET_CROSSINGS`] seam crossings.
    Path,
    /// [`RANDOM_STEPS`] steps of a randomised heading in a bounded region.
    Random,
}

/// Drop the body onto the terrain and let it settle.
fn settle_mesh(terrain: &mut Terrain, x: f32, z: f32) -> MeshBody {
    let mut body = MeshBody {
        pos: Vector::new(x, SPAWN_HEIGHT, z),
        vy: 0.0,
        grounded: false,
    };
    let mut blockers = Vec::new();
    for _ in 0..SETTLE_STEPS {
        terrain.ensure(body.pos.x, body.pos.z);
        body.step(terrain, Vector::ZERO, &mut blockers);
    }
    body
}

/// Run one arm.
fn run_arm(fixture: Fixture, name: &'static str) -> Arm {
    let field = FbmTerrain::<f32>::canonical();
    let mut terrain = Terrain::new(field);

    let (start_x, start_z) = match fixture {
        Fixture::Path => path(0.0),
        Fixture::Random => (RANDOM_EXTENT * 0.5, RANDOM_EXTENT * 0.5),
    };
    let mut body = settle_mesh(&mut terrain, start_x, start_z);
    assert!(
        body.grounded,
        "VOID: the {name} arm's body never reached the ground in {SETTLE_STEPS} settle steps, so \
         every measured step would be a measurement of falling"
    );

    let mut rng = Rng(0x5EED_1234_ABCD_EF01);
    let mut heading_angle = rng.next() * std::f32::consts::TAU;
    let mut commanded_distance = 0.0_f32;
    let mut column = terrain.layout.chunk_of([body.pos.x, 0.0, body.pos.z]);

    let mut arm = Arm {
        steps: 0,
        seam_crossings: 0,
        metres: 0.0,
        population: Vec::new(),
        stops: Vec::new(),
        stops_field_walk: 0,
        unattributed: 0,
    };
    let mut blockers = Vec::new();
    let mut commands: Vec<Vector> = Vec::new();
    let mut stopped_last_step = false;

    loop {
        let done = match fixture {
            Fixture::Path => arm.seam_crossings >= TARGET_CROSSINGS,
            Fixture::Random => arm.steps >= RANDOM_STEPS,
        };
        if done {
            break;
        }
        assert!(
            arm.steps < MAX_PATH_STEPS,
            "VOID: the {name} arm ran {MAX_PATH_STEPS} steps and reached only {} of \
             {TARGET_CROSSINGS} seam crossings, so the fixture never delivered the registered \
             crossing count",
            arm.seam_crossings
        );

        // Steering. The path arm aims at the point on the path a little ahead of
        // the commanded distance, exactly as `drive` does; the randomised arm
        // turns by a bounded random increment and is reflected at the region
        // border so the meshed set stays bounded.
        let heading = match fixture {
            Fixture::Path => {
                commanded_distance += WALK_SPEED * DT;
                let (tx, tz) = path(commanded_distance);
                let to = Vector::new(tx - body.pos.x, 0.0, tz - body.pos.z);
                if to.length() > 1.0e-6 {
                    to.normalize()
                } else {
                    Vector::X
                }
            }
            Fixture::Random => {
                heading_angle += (rng.next() - 0.5) * 0.7;
                let mut h = Vector::new(heading_angle.cos(), 0.0, heading_angle.sin());
                let inside = |v: f32| v > 1.0 && v < RANDOM_EXTENT - 1.0;
                if !inside(body.pos.x + h.x) || !inside(body.pos.z + h.z) {
                    let back = Vector::new(
                        RANDOM_EXTENT * 0.5 - body.pos.x,
                        0.0,
                        RANDOM_EXTENT * 0.5 - body.pos.z,
                    );
                    if back.length() > 1.0e-6 {
                        h = back.normalize();
                        heading_angle = h.z.atan2(h.x);
                    }
                }
                h
            }
        };
        let commanded = heading * (WALK_SPEED * DT);
        commands.push(commanded);

        terrain.ensure(body.pos.x, body.pos.z);
        let before = body.pos;
        let vy_before = body.vy;
        let grounded_before = body.grounded;
        let outcome = body.step(&mut terrain, commanded, &mut blockers);

        arm.steps += 1;
        arm.metres += outcome.moved;
        let now = terrain.layout.chunk_of([body.pos.x, 0.0, body.pos.z]);
        if now.coords[0] != column.coords[0] || now.coords[2] != column.coords[2] {
            arm.seam_crossings += 1;
            column = now;
        }
        assert!(
            body.pos.y > -8.0 && body.pos.y < SPAWN_HEIGHT + 2.0,
            "VOID: the {name} arm's body left the meshed layers at y = {} after {} steps, so its \
             stops are falls rather than blocks",
            body.pos.y,
            arm.steps
        );

        if outcome.shortfall < STALL_SHORTFALL || outcome.asked <= f32::EPSILON {
            stopped_last_step = false;
            continue;
        }
        if blockers.is_empty() {
            arm.unattributed += 1;
            stopped_last_step = false;
            continue;
        }
        let first = terrain.tris[blockers[0] as usize];
        let worst = blockers
            .iter()
            .map(|&i| terrain.tris[i as usize].q)
            .fold(f32::INFINITY, f32::min);
        let (field_shortfall, field_pos) =
            field_replay(&field, before, vy_before, grounded_before, commanded);
        arm.stops.push(Stop {
            q: first.q,
            q_worst_blocker: worst,
            seam: first.seam,
            seam_xz: first.seam_xz,
            degenerate: first.degenerate,
            field_stalls: field_shortfall >= STALL_SHORTFALL,
            slope: first.slope,
            field_penetration: penetration(&terrain, field_pos),
            mesh_penetration: penetration(&terrain, body.pos),
            first_of_episode: !stopped_last_step,
        });
        stopped_last_step = true;
    }

    // The independent field-only walk over the same commanded sequence. Its
    // trajectory diverges from the mesh arm's by construction — that is why C3 is
    // a replay and this is only a reported number.
    let mut field_body = FieldBody {
        pos: Vector::new(start_x, SPAWN_HEIGHT, start_z),
        vy: 0.0,
        grounded: false,
    };
    for _ in 0..SETTLE_STEPS {
        field_body.step(&field, Vector::ZERO);
    }
    for commanded in &commands {
        if field_body.step(&field, *commanded).shortfall >= STALL_SHORTFALL {
            arm.stops_field_walk += 1;
        }
    }

    arm.population = terrain.population();
    arm
}

// ── controls that are not the measurement ────────────────────────────────────

/// A floor with a vertical wall, for the field controller's instrument control.
///
/// Union of two half-spaces, so the value is the smaller of the two and the
/// gradient is the active one's. Solid below `y = 0` and beyond `x = WALL_AT`.
struct FloorAndWall;

/// Where [`FloorAndWall`]'s wall stands.
const WALL_AT: f32 = 4.0;

impl Sdf for FloorAndWall {
    type Scalar = f32;

    fn sample(&self, p: [f32; 3]) -> f32 {
        p[1].min(WALL_AT - p[0])
    }

    fn gradient(&self, p: [f32; 3]) -> [f32; 3] {
        if p[1] <= WALL_AT - p[0] {
            [0.0, 1.0, 0.0]
        } else {
            [-1.0, 0.0, 0.0]
        }
    }
}

/// Prove the field controller can report a stall at all.
///
/// `stops_field_path == 0` is C3's HELD, and a zero from an instrument that
/// cannot produce a non-zero is `M-44`. So the same [`FieldBody::step`], with the
/// same body and the same threshold, is driven into a vertical wall; the count it
/// returns is asserted non-zero.
fn field_instrument_control() -> u64 {
    let mut body = FieldBody {
        pos: Vector::new(0.0, CAPSULE_HALF + CAPSULE_RADIUS, 0.0),
        vy: 0.0,
        grounded: false,
    };
    for _ in 0..SETTLE_STEPS {
        body.step(&FloorAndWall, Vector::ZERO);
    }
    let mut stalls = 0;
    let commanded = Vector::new(WALK_SPEED * DT, 0.0, 0.0);
    for _ in 0..600 {
        if body.step(&FloorAndWall, commanded).shortfall >= STALL_SHORTFALL {
            stalls += 1;
        }
    }
    stalls
}

/// The height of `fbm_terrain` under `(x, z)`.
///
/// Exact rather than sampled: the field is `f(p) = p.y − h(x, z)`, so
/// `h = −f([x, 0, z])`.
fn height(field: &FbmTerrain<f32>, x: f32, z: f32) -> f32 {
    -field.sample([x, 0.0, z])
}

/// The largest hover the **chord** probe can stop a fall at, at this foot.
///
/// Derived from the probe's own geometry rather than fitted. The probe answers
/// when any of five samples is solid, so the fall stops at the *largest* of their
/// individual thresholds, with `hover` measured from the capsule's bottom to the
/// terrain directly under the foot:
///
/// ```text
/// centre sample  : hover <= GROUND_PROBE
/// lateral sample : hover <= GROUND_PROBE − (R − R·sqrt(1 − k²)) + rise
/// ```
///
/// where `k` is [`FOOTPRINT`] and `rise` is the terrain's own rise at that
/// lateral offset — taken from the field, so curvature is included and nothing is
/// linearised. `✗45`'s flat cross is the same expression with the inset term
/// deleted, which is why it hovers higher wherever `rise > 0`, and why that
/// difference is the whole of `M-363`.
fn chord_hover_bound(field: &FbmTerrain<f32>, foot: Vector) -> f32 {
    let r = CAPSULE_RADIUS * FOOTPRINT;
    let inset = CAPSULE_RADIUS - CAPSULE_RADIUS * (1.0 - FOOTPRINT * FOOTPRINT).sqrt();
    let h0 = height(field, foot.x, foot.z);
    let mut bound = GROUND_PROBE;
    for [dx, dz] in [[r, 0.0], [-r, 0.0], [0.0, r], [0.0, -r]] {
        let rise = height(field, foot.x + dx, foot.z + dz) - h0;
        bound = bound.max(GROUND_PROBE - inset + rise);
    }
    bound
}

/// What the ground-probe lift check measured.
struct LiftCheck {
    /// Worst resting hover of the chord probe under gravity, over the steep
    /// columns.
    chord_settle: f32,
    /// The same for `✗45`'s flat cross.
    flat_settle: f32,
    /// Worst resting hover of the chord probe where the column it settled in is
    /// no steeper than the 0.63 gradient `M-363` measured on `Ground`.
    chord_gentle: f32,
    /// Worst amount by which a **chord**-probe body rested above
    /// [`chord_hover_bound`] at its own settled foot. Must be ~0.
    chord_excess: f32,
    /// The same for a **flat-cross** body. Must be large: this is the built-in
    /// mutation test, and it is what makes the chord's ~0 mean something.
    flat_excess: f32,
    /// Mean of `flat_trigger − chord_trigger` over the columns where a lateral
    /// sample is what answers.
    inset_measured: f32,
    /// `R − R·sqrt(1 − k²)`, the chord's own inset. The number
    /// [`LiftCheck::inset_measured`] must reproduce.
    inset_expected: f32,
    /// Steepest `|∇h|` found in the scan.
    steepest: f32,
    /// Columns measured.
    columns: usize,
}

/// Prove the transcribed ground probe is `M-363`'s chord and not `✗45`'s flat
/// cross.
///
/// Three measurements, because the first two versions of this check both measured
/// something else and said so:
///
/// * **A quasi-static trigger height.** The body is lowered in 1 mm steps with no
///   gravity and no push-out until the probe answers. That is the probe's geometry
///   and nothing else, and the difference between the two geometries is then a
///   *prediction*: wherever a lateral sample is what answers, the flat cross
///   triggers exactly [`LiftCheck::inset_expected`] higher, because deleting the
///   `sqrt(R² − r²)` term is the whole of the difference. The caller asserts that
///   equality to a millimetre.
/// * **A resting hover under gravity**, checked against [`chord_hover_bound`] at
///   the column the body settled in — the push-out is horizontal on a slope, so a
///   body slides out of the column it was dropped in, and bucketing by the *spawn*
///   column read a 0.2342 hover against a 0.63 gradient in the first version.
/// * **The same under the flat cross**, whose excess over the chord's bound must
///   be large. Without it, "the chord obeys its bound" is a claim about an
///   assertion nobody has seen fail.
///
/// The spawn height is also a measurement, and the second version got it wrong:
/// dropped from 0.25 above the terrain the probe answers *at the spawn*, gravity
/// is cut on the first step, and every column reports 0.245 — the spawn height
/// minus one step. It is now above the largest trigger height the geometry admits.
fn ground_probe_lift_check(field: &FbmTerrain<f32>) -> LiftCheck {
    let inset_expected = CAPSULE_RADIUS - CAPSULE_RADIUS * (1.0 - FOOTPRINT * FOOTPRINT).sqrt();

    /// Hover the body is lowered from. Above `GROUND_PROBE + max rise`, which on
    /// this field is under 0.7.
    const FROM: f32 = 1.5;
    /// Lowering step, and so the resolution of the trigger measurement.
    const STEP: f32 = 0.001;

    let trigger = |x: f32, z: f32, probe: GroundProbe| -> f32 {
        let base = height(field, x, z) + CAPSULE_HALF + CAPSULE_RADIUS;
        let mut hover = FROM;
        while hover > -CAPSULE_RADIUS {
            let mut pos = Vector::new(x, base + hover, z);
            let mut velocity = Vector::ZERO;
            // The push-out cannot fire while `hover > 0`, so this reads the probe
            // and only the probe.
            if resolve_field(field, &mut pos, &mut velocity, probe) {
                return hover;
            }
            hover -= STEP;
        }
        hover
    };

    let settle = |x: f32, z: f32, probe: GroundProbe| -> (f32, Vector) {
        let start = height(field, x, z) + CAPSULE_HALF + CAPSULE_RADIUS + FROM;
        let mut pos = Vector::new(x, start, z);
        let mut vy = 0.0_f32;
        let mut grounded = false;
        for _ in 0..600 {
            if grounded && vy <= 0.0 {
                vy = 0.0;
            } else {
                vy = (vy - GRAVITY * DT).max(-MAX_FALL);
            }
            pos.y += vy * DT;
            let mut velocity = Vector::new(0.0, vy, 0.0);
            grounded = resolve_field(field, &mut pos, &mut velocity, probe);
            vy = velocity.y;
        }
        let foot = pos - Vector::Y * CAPSULE_HALF;
        let hover = pos.y - CAPSULE_HALF - CAPSULE_RADIUS - height(field, pos.x, pos.z);
        (hover, foot)
    };

    // Two column sets, both selected from the field's own analytic gradient
    // rather than hoped for: steep ones, where the two probe geometries differ,
    // and gentle ones at the gradient `M-363` measured `Ground` at.
    let mut steepest = 0.0_f32;
    let mut steep: Vec<(f32, f32)> = Vec::new();
    let mut gentle: Vec<(f32, f32)> = Vec::new();
    let mut scan = Rng(0xC0FF_EE00_1234_5678);
    while steep.len() < 256 || gentle.len() < 256 {
        let x = scan.next() * 64.0;
        let z = scan.next() * 64.0;
        let g = field.gradient([x, 0.0, z]);
        let slope = g[0].hypot(g[2]);
        steepest = steepest.max(slope);
        if slope > 1.0 {
            if steep.len() < 256 {
                steep.push((x, z));
            }
        } else if slope < 0.4 && gentle.len() < 256 {
            gentle.push((x, z));
        }
    }

    let mut chord_settle = 0.0_f32;
    let mut flat_settle = 0.0_f32;
    let mut chord_gentle = 0.0_f32;
    let mut chord_excess = f32::NEG_INFINITY;
    let mut flat_excess = f32::NEG_INFINITY;
    let mut inset_sum = 0.0_f32;
    let mut inset_n = 0_u32;
    let mut gentle_seen = 0_usize;
    for (set, is_steep) in [(&steep, true), (&gentle, false)] {
        for &(x, z) in set {
            let (hover, foot) = settle(x, z, GroundProbe::Chord);
            chord_excess = chord_excess.max(hover - chord_hover_bound(field, foot));
            let (flat_hover, flat_foot) = settle(x, z, GroundProbe::FlatCross);
            flat_excess = flat_excess.max(flat_hover - chord_hover_bound(field, flat_foot));
            if is_steep {
                chord_settle = chord_settle.max(hover);
                flat_settle = flat_settle.max(flat_hover);
            }
            let g = field.gradient([foot.x, 0.0, foot.z]);
            if g[0].hypot(g[2]) <= 0.63 {
                chord_gentle = chord_gentle.max(hover);
                gentle_seen += 1;
            }

            // Only where a lateral sample is what answers: where the centre
            // sample answers first the two geometries agree by construction, and
            // averaging those in would dilute a prediction into an estimate.
            let chord_trigger = trigger(x, z, GroundProbe::Chord);
            let flat_trigger = trigger(x, z, GroundProbe::FlatCross);
            if chord_trigger > GROUND_PROBE + STEP {
                inset_sum += flat_trigger - chord_trigger;
                inset_n += 1;
            }
        }
    }
    assert!(
        gentle_seen > 0,
        "VOID: no body settled in a column no steeper than 0.63, so M-363's own hover figure has \
         nothing to be compared against"
    );
    assert!(
        inset_n > 0,
        "VOID: no column had a lateral sample answer the ground probe, so the two probe \
         geometries cannot be told apart and this check asserts nothing"
    );
    LiftCheck {
        chord_settle,
        flat_settle,
        chord_gentle,
        chord_excess,
        flat_excess,
        inset_measured: inset_sum / inset_n as f32,
        inset_expected,
        steepest,
        columns: steep.len() + gentle.len(),
    }
}

// ── rows ─────────────────────────────────────────────────────────────────────

/// The tenth percentile of `q`, the count at or below it, and its area share.
fn decile(population: &[TriMeta]) -> (f32, usize, f64) {
    let mut sorted: Vec<f32> = population.iter().map(|t| t.q).collect();
    sorted.sort_unstable_by(f32::total_cmp);
    let p10 = sorted[sorted.len() / 10];
    let count = sorted.partition_point(|&q| q <= p10);
    let total: f64 = population.iter().map(|t| f64::from(t.two_area)).sum();
    let inside: f64 = population
        .iter()
        .filter(|t| t.q <= p10)
        .map(|t| f64::from(t.two_area))
        .sum();
    (p10, count, inside / total)
}

/// One `stops per triangle` comparison: a seam predicate over a population.
struct Rate {
    /// Seam-adjacent triangles in the population.
    seam_tris: u64,
    /// The rest of the population.
    interior_tris: u64,
    /// Stops attributed to a seam-adjacent triangle.
    seam_stops: u64,
    /// Stops attributed to any other triangle.
    interior_stops: u64,
    /// `seam_stops / seam_tris`.
    per_seam: f64,
    /// `interior_stops / interior_tris`.
    per_interior: f64,
    /// `per_seam / per_interior`. C2 holds at or below 1.
    excess: f64,
}

/// Compute one comparison.
///
/// `tested_only` restricts the **denominator** to triangles the broadphase handed
/// to the controller; every numerator is already inside that set, because a
/// blocker was a candidate. `xz_only` selects `game_capsule_walk::measure`'s own
/// seam test over the all-axis one.
fn rate(population: &[TriMeta], stops: &[&Stop], tested_only: bool, xz_only: bool) -> Rate {
    let mut seam_tris = 0_u64;
    let mut interior_tris = 0_u64;
    for t in population {
        if tested_only && !t.tested {
            continue;
        }
        if if xz_only { t.seam_xz } else { t.seam } {
            seam_tris += 1;
        } else {
            interior_tris += 1;
        }
    }
    let seam_stops = stops
        .iter()
        .filter(|s| if xz_only { s.seam_xz } else { s.seam })
        .count() as u64;
    let interior_stops = stops.len() as u64 - seam_stops;
    let per_seam = seam_stops as f64 / seam_tris as f64;
    let per_interior = interior_stops as f64 / interior_tris as f64;
    Rate {
        seam_tris,
        interior_tris,
        seam_stops,
        interior_stops,
        per_seam,
        per_interior,
        excess: if per_interior > 0.0 {
            per_seam / per_interior
        } else {
            f64::INFINITY
        },
    }
}

/// Mean of `f` over the triangles `keep` accepts, or zero over an empty set.
fn mean_slope(population: &[TriMeta], keep: impl Fn(&TriMeta) -> bool) -> f64 {
    let mut sum = 0.0_f64;
    let mut n = 0_u64;
    for t in population.iter().filter(|t| keep(t)) {
        sum += f64::from(t.slope);
        n += 1;
    }
    if n == 0 { 0.0 } else { sum / n as f64 }
}

/// Assemble one CSV row from one or more arms pooled.
fn row(
    label: &str,
    arms: &[&Arm],
    field_control_stalls: u64,
    lift: &LiftCheck,
    wall: f64,
) -> Vec<(&'static str, String)> {
    let population: Vec<TriMeta> = arms
        .iter()
        .flat_map(|a| a.population.iter().copied())
        .collect();
    let stops: Vec<&Stop> = arms.iter().flat_map(|a| a.stops.iter()).collect();
    let steps: u64 = arms.iter().map(|a| a.steps).sum();
    let crossings: u64 = arms.iter().map(|a| a.seam_crossings).sum();
    let metres: f32 = arms.iter().map(|a| a.metres).sum();
    let field_walk: u64 = arms.iter().map(|a| a.stops_field_walk).sum();
    let unattributed: u64 = arms.iter().map(|a| a.unattributed).sum();

    // ── vacuity control 1, registered: the decile must be non-empty ──────────
    //
    // Over the whole meshed population, which is what "the aspect-ratio
    // distribution" names. The tested-only percentile is reported beside it and
    // the two agree to three digits, which is itself worth knowing: the swath the
    // capsule tested is not a shape-biased sample of the mesh.
    let (p10, decile_count, decile_area) = decile(&population);
    assert!(
        decile_count > 0,
        "VACUOUS ({label}): the bottom aspect-ratio decile is empty over {} triangles, so C1 \
         asks about a set that does not exist",
        population.len()
    );
    let tested: Vec<TriMeta> = population.iter().copied().filter(|t| t.tested).collect();
    assert!(
        !tested.is_empty(),
        "VACUOUS ({label}): the controller tested no triangle at all"
    );
    let (p10_tested, decile_tested, _) = decile(&tested);

    // ── vacuity control 2, registered: a non-zero stop count ─────────────────
    let stops_total = stops.len() as u64;
    assert!(
        stops_total > 0,
        "VACUOUS ({label}): {steps} steps produced no stop at all, so every fraction below is \
         0/0 and the walk measured nothing"
    );

    // ── vacuity control 3, added: attribution must be complete ──────────────
    //
    // A stop with no blocking triangle is a stop this harness cannot classify,
    // and silently dropping those would move `worst_decile_fraction` by an
    // unknown amount in an unknown direction.
    assert_eq!(
        unattributed, 0,
        "{label}: {unattributed} stops had no blocking triangle, so the attribution is \
         incomplete and every fraction below is over a subset of unknown size"
    );

    // ── vacuity control 4, added: C3's instrument must be able to fire ───────
    assert!(
        field_control_stalls > 0,
        "VACUOUS ({label}): the field controller reported no stall against a vertical wall, so \
         `stops_field_path == 0` would be a zero that could not have been non-zero"
    );

    // ── C2, four ways ───────────────────────────────────────────────────────
    //
    // The registered columns carry `xz_only` over the tested population, and both
    // choices are corrections this harness's own counts forced:
    //
    // * **The `y` band is not a seam test.** A band one cell either side of
    //   `y = 0` selects the terrain whose *height* is near zero, and on an fbm
    //   heightfield that is the zero-crossing band — where the surface is
    //   steepest. It put 67.5% of all triangles in the "seam" bucket and its
    //   excess is a slope excess. `mean_slope_y_band` against
    //   `mean_slope_off_y_band` is the number that says so. `game_capsule_walk`
    //   tests `x` and `z`, and C2 exists to confirm `M-115` from the controller's
    //   side, so `x`/`z` is the definition that answers the registered clause.
    // * **`ensure` inflates the interior denominator.** It meshes a whole chunk
    //   when the body comes within `REACH` of it, so a clipped chunk contributes
    //   ~500 triangles of which only the ones by its own boundary were reachable
    //   — a seam excess manufactured by the meshing policy. The tested set is the
    //   denominator that cannot do that.
    let primary = rate(&population, &stops, true, true);
    let all_meshed_xz = rate(&population, &stops, false, true);
    let tested_all_axes = rate(&population, &stops, true, false);
    let all_meshed_all_axes = rate(&population, &stops, false, false);
    assert!(
        primary.seam_tris > 0 && primary.interior_tris > 0,
        "VACUOUS ({label}): {} seam and {} interior triangles were tested, so C2's comparison \
         has an empty side",
        primary.seam_tris,
        primary.interior_tris
    );

    // The same comparison over **episodes**, which is the effective sample size:
    // a body climbing one face stalls on every frame of the climb, so the raw
    // stop count is autocorrelated and a Poisson bar over it is several times too
    // tight. `seam_excess_episodes_sigma` is the one-sigma spread of the ratio
    // from the two episode counts alone, which is what says whether a 1.13 excess
    // is a finding or a fixture.
    let episodes: Vec<&Stop> = stops
        .iter()
        .copied()
        .filter(|s| s.first_of_episode)
        .collect();
    let episode_rate = rate(&population, &episodes, true, true);
    let episode_sigma = episode_rate.excess
        * (1.0 / episode_rate.seam_stops.max(1) as f64
            + 1.0 / episode_rate.interior_stops.max(1) as f64)
            .sqrt();

    let on_worst = stops.iter().filter(|s| s.q <= p10).count() as u64;
    let on_worst_tested = stops.iter().filter(|s| s.q <= p10_tested).count() as u64;
    let on_worst_any = stops.iter().filter(|s| s.q_worst_blocker <= p10).count() as u64;
    let on_degenerate = stops.iter().filter(|s| s.degenerate).count() as u64;
    let degenerate_tris = population.iter().filter(|t| t.degenerate).count() as u64;
    // **What makes `stops_on_degenerate == 0` a measurement.** The crate's
    // `degenerate_triangles` metric counts triangles of numerically zero area, and
    // a zero stop count over triangles the capsule was never within reach of would
    // be `M-44` exactly. This is the count the controller actually tested.
    let degenerate_tested = population
        .iter()
        .filter(|t| t.degenerate && t.tested)
        .count() as u64;
    assert!(
        degenerate_tris == 0 || degenerate_tested > 0,
        "VACUOUS ({label}): {degenerate_tris} degenerate triangles exist and the controller was \
         handed none of them, so `stops_on_degenerate` could not have been non-zero"
    );

    let fraction = on_worst as f64 / stops_total as f64;

    // C3 is over C1's stops: the bottom-decile ones.
    let field_path = stops
        .iter()
        .filter(|s| s.q <= p10 && s.field_stalls)
        .count() as u64;
    let field_all = stops.iter().filter(|s| s.field_stalls).count() as u64;
    let field_pen_max = stops
        .iter()
        .map(|s| s.field_penetration)
        .fold(0.0_f32, f32::max);
    let field_pen_mean = stops
        .iter()
        .map(|s| f64::from(s.field_penetration))
        .sum::<f64>()
        / stops_total as f64;
    let mesh_pen_max = stops
        .iter()
        .map(|s| s.mesh_penetration)
        .fold(0.0_f32, f32::max);
    let slope_at_stops = stops.iter().map(|s| f64::from(s.slope)).sum::<f64>() / stops_total as f64;

    // **What C3's zero is worth.** The same replay fires `field_all` times over
    // all `stops_total` stops, so under the null that a bottom-decile stop is an
    // ordinary stop the expected count inside C1's set is this — and if it is of
    // order one, `stops_field_path == 0` is a weak HELD rather than a strong one.
    // Reported rather than argued, and it is the honest reading of the clause.
    let field_path_expected = on_worst as f64 * field_all as f64 / stops_total as f64;

    let c1 = fraction >= 0.20;
    let c2 = primary.excess <= 1.0;
    let c3 = field_path == 0;

    vec![
        ("field", String::from("fbm_terrain")),
        ("arm", String::from(label)),
        ("steps", steps.to_string()),
        ("seam_crossings", crossings.to_string()),
        ("stops_total", stops_total.to_string()),
        ("stops_on_worst_decile", on_worst.to_string()),
        ("worst_decile_fraction", format!("{fraction:.6}")),
        (
            "stops_per_triangle_seam",
            format!("{:.9}", primary.per_seam),
        ),
        (
            "stops_per_triangle_interior",
            format!("{:.9}", primary.per_interior),
        ),
        ("seam_excess_ratio", format!("{:.6}", primary.excess)),
        ("stops_field_path", field_path.to_string()),
        ("aspect_ratio_p10", format!("{p10:.6}")),
        ("c1_holds", c1.to_string()),
        ("c2_holds", c2.to_string()),
        ("c3_holds", c3.to_string()),
        // ── the C1 null rate ────────────────────────────────────────────────
        ("triangles", population.len().to_string()),
        ("triangles_tested", tested.len().to_string()),
        ("decile_triangles", decile_count.to_string()),
        ("decile_triangles_tested", decile_tested.to_string()),
        ("worst_decile_area_share", format!("{decile_area:.8}")),
        ("c1_bar_over_chance", format!("{:.2}", 0.20 / decile_area)),
        ("stops_any_blocker_worst_decile", on_worst_any.to_string()),
        ("aspect_ratio_p10_tested", format!("{p10_tested:.6}")),
        ("stops_on_worst_decile_tested", on_worst_tested.to_string()),
        ("degenerate_triangles", degenerate_tris.to_string()),
        ("degenerate_triangles_tested", degenerate_tested.to_string()),
        ("stops_on_degenerate", on_degenerate.to_string()),
        // ── C2, the three variants the registered columns are not ───────────
        ("seam_triangles_tested_xz", primary.seam_tris.to_string()),
        (
            "interior_triangles_tested_xz",
            primary.interior_tris.to_string(),
        ),
        ("stops_on_seam_xz", primary.seam_stops.to_string()),
        ("stops_on_interior_xz", primary.interior_stops.to_string()),
        (
            "seam_excess_all_meshed_xz",
            format!("{:.6}", all_meshed_xz.excess),
        ),
        (
            "seam_excess_tested_all_axes",
            format!("{:.6}", tested_all_axes.excess),
        ),
        (
            "seam_excess_all_meshed_all_axes",
            format!("{:.6}", all_meshed_all_axes.excess),
        ),
        (
            "seam_triangles_all_axes",
            all_meshed_all_axes.seam_tris.to_string(),
        ),
        (
            "stops_on_seam_all_axes",
            all_meshed_all_axes.seam_stops.to_string(),
        ),
        // ── the effective sample size ───────────────────────────────────────
        ("stop_episodes", episodes.len().to_string()),
        ("episode_seam_stops", episode_rate.seam_stops.to_string()),
        (
            "episode_interior_stops",
            episode_rate.interior_stops.to_string(),
        ),
        (
            "seam_excess_episodes",
            format!("{:.6}", episode_rate.excess),
        ),
        ("seam_excess_episodes_sigma", format!("{episode_sigma:.6}")),
        // ── the slope confound, measured ────────────────────────────────────
        (
            "mean_slope_y_band",
            format!("{:.4}", mean_slope(&population, |t| t.seam_y)),
        ),
        (
            "mean_slope_off_y_band",
            format!("{:.4}", mean_slope(&population, |t| !t.seam_y)),
        ),
        (
            "mean_slope_xz_band",
            format!("{:.4}", mean_slope(&population, |t| t.seam_xz)),
        ),
        (
            "mean_slope_off_xz_band",
            format!("{:.4}", mean_slope(&population, |t| !t.seam_xz)),
        ),
        ("mean_slope_at_stops", format!("{slope_at_stops:.4}")),
        (
            "stops_field_path_expected",
            format!("{field_path_expected:.3}"),
        ),
        // ── C3's mechanism ─────────────────────────────────────────────────
        ("stops_field_replay_all", field_all.to_string()),
        ("stops_field_walk", field_walk.to_string()),
        ("field_control_stalls", field_control_stalls.to_string()),
        ("field_penetration_max", format!("{field_pen_max:.4}")),
        ("field_penetration_mean", format!("{field_pen_mean:.4}")),
        ("mesh_penetration_max", format!("{mesh_pen_max:.4}")),
        // ── the lift check and the clock ────────────────────────────────────
        ("chord_hover_max", format!("{:.6}", lift.chord_settle)),
        (
            "chord_hover_gentle_max",
            format!("{:.6}", lift.chord_gentle),
        ),
        ("flat_cross_hover_max", format!("{:.6}", lift.flat_settle)),
        (
            "chord_excess_over_bound",
            format!("{:.6}", lift.chord_excess),
        ),
        ("flat_excess_over_bound", format!("{:.6}", lift.flat_excess)),
        (
            "chord_inset_measured",
            format!("{:.6}", lift.inset_measured),
        ),
        (
            "chord_inset_expected",
            format!("{:.6}", lift.inset_expected),
        ),
        ("steepest_gradient", format!("{:.4}", lift.steepest)),
        ("lift_check_columns", lift.columns.to_string()),
        ("metres_travelled", format!("{metres:.2}")),
        ("wall_seconds", format!("{wall:.1}")),
        ("dt_seconds", format!("{DT:.6}")),
        ("stall_shortfall", format!("{STALL_SHORTFALL:.2}")),
    ]
}

fn main() {
    if !std::env::args().any(|arg| arg == "--bench") {
        return;
    }

    let started = Instant::now();
    common::experiment::run(isomesh::experiment!("P-86"), |run| {
        let field = FbmTerrain::<f32>::canonical();

        // ── the lift checks, before any measurement ─────────────────────────
        let lift = ground_probe_lift_check(&field);
        println!(
            "  lift check over {} columns (steepest |∇h| {:.3}): chord inset measured {:.4} \
             against the geometric {:.4}; resting hover chord {:.4} flat-cross {:.4}; excess \
             over the chord's own bound: chord {:+.5}, flat-cross {:+.5}; chord hover where the \
             settled column is no steeper than 0.63 is {:.4}, against M-363's 0.042 on Ground",
            lift.columns,
            lift.steepest,
            lift.inset_measured,
            lift.inset_expected,
            lift.chord_settle,
            lift.flat_settle,
            lift.chord_excess,
            lift.flat_excess,
            lift.chord_gentle
        );
        // The prediction: deleting the `sqrt(R² − r²)` term is the whole
        // difference between the two geometries, so wherever a lateral sample is
        // what answers, the flat cross triggers exactly one inset higher.
        assert!(
            (lift.inset_measured - lift.inset_expected).abs() <= 0.002,
            "the transcribed ground probe is not M-363's chord: the flat cross triggers {:.4} \
             higher, and the chord's own inset R(1 − sqrt(1 − k²)) is {:.4}",
            lift.inset_measured,
            lift.inset_expected
        );
        // The body stops at the first position where the probe answers, so it can
        // rest one step of descent above the exact threshold and no further.
        assert!(
            lift.chord_excess <= 0.02,
            "the transcribed chord probe does not obey its own geometry: some settled body rests \
             {:+.5} above the bound its five sample points allow",
            lift.chord_excess
        );
        // And the mutation test for the assertion above: the geometry ✗45 removed
        // must fail it, or "the chord obeys its bound" is a claim about a check
        // nobody has seen fire.
        assert!(
            lift.flat_excess > 0.02,
            "the bound check cannot discriminate: ✗45's flat cross rests only {:+.5} above the \
             chord's bound, so the chord passing it means nothing",
            lift.flat_excess
        );

        let field_control_stalls = field_instrument_control();
        println!(
            "  field instrument control: {field_control_stalls} stalls against a vertical wall"
        );

        let path_arm = run_arm(Fixture::Path, "path_495_crossings");
        println!(
            "  path arm:   {} steps, {} crossings, {} stops, {} triangles",
            path_arm.steps,
            path_arm.seam_crossings,
            path_arm.stops.len(),
            path_arm.population.len()
        );
        let random_arm = run_arm(Fixture::Random, "random_10k_steps");
        println!(
            "  random arm: {} steps, {} crossings, {} stops, {} triangles",
            random_arm.steps,
            random_arm.seam_crossings,
            random_arm.stops.len(),
            random_arm.population.len()
        );

        let wall = started.elapsed().as_secs_f64();
        for (label, arms) in [
            ("path_495_crossings", vec![&path_arm]),
            ("random_10k_steps", vec![&random_arm]),
            ("combined", vec![&path_arm, &random_arm]),
        ] {
            let values = row(label, &arms, field_control_stalls, &lift, wall);
            for (key, value) in &values {
                if [
                    "stops_total",
                    "stops_on_worst_decile",
                    "worst_decile_fraction",
                    "worst_decile_area_share",
                    "seam_excess_ratio",
                    "seam_excess_all_meshed_xz",
                    "seam_excess_tested_all_axes",
                    "stops_field_path",
                    "stops_field_replay_all",
                    "aspect_ratio_p10",
                    "mean_slope_y_band",
                    "mean_slope_off_y_band",
                    "mean_slope_xz_band",
                    "mean_slope_off_xz_band",
                    "mean_slope_at_stops",
                    "field_penetration_max",
                    "field_penetration_mean",
                    "mesh_penetration_max",
                ]
                .contains(key)
                {
                    println!("  {label:>18} {key:>24} = {value}");
                }
            }
            run.record(&values);
        }
    });
}
