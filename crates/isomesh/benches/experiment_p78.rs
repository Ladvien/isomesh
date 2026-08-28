//! **P-78 — how many light probes one dig invalidates.**
//!
//! Ticket: R-078. Pre-registered before this harness existed.
//!
//! ```bash
//! cargo bench --bench experiment_p78
//! ```
//!
//! Writes `docs/experiments/p-78.csv`.
//!
//! # H, and the falsifier
//!
//! Probe-based GI is the only GI family that tolerates geometry changing every
//! frame, and its cost under a destructible world is unmeasured. The best
//! rigorous datum is a bachelor's thesis (Dell'Ova, BTH 2025, **DOI
//! unverified** — flagged at registration and flagged again here) finding
//! **87–96% of radiance-cascade frame time is the gather pass**.
//!
//! - **C1.** The probe invalidation set is edit-proportional, not
//!   volume-proportional: it tracks the brush's dilated support with a constant
//!   factor **under 4**, across M-50's four edit-log buckets and three probe
//!   densities. Falsified by a factor above 4, **or by a count that grows with
//!   world size at fixed edit size**.
//! - **C2.** Tracing the field beats tracing the extracted mesh for probe
//!   updates by at least **3×** at equal ray count.
//! - **C3.** A dig that opens a new air component invalidates **strictly more**
//!   probes than a dig of the same volume that does not.
//!
//! **Registered vacuity control:** the breakthrough arm must actually merge two
//! air components, *asserted from the union-find rather than assumed from the
//! brush*.
//!
//! # SHARE, recomputed before the code was written
//!
//! **C2's share is the gather pass at 87–96%, and the clause is reachable**
//! because C2 is a ratio *on the gather itself* rather than on a frame. A 3×
//! gather win is a 2.38× (at 87%) to 2.78× (at 96%) win on total GI cost by
//! Amdahl. Nothing about the arithmetic makes 3× unreachable.
//!
//! **C1's ratio needs a unit decision, and one of the two readings is
//! unfalsifiable by arithmetic alone.** `probes_invalidated` is a count of
//! probes and `brush_support_cells` is a count of cells, so the literal quotient
//! of the two registered columns is density-dependent by construction. Its
//! *ceiling* — every probe in the whole world invalidated — is
//! `probes_total / brush_support_cells`, and with a radius-4 brush whose support
//! is ~10³ cells that ceiling is
//!
//! | world | s = 2 | s = 4 | s = 8 |
//! |---|---:|---:|---:|
//! | 32³ = 32,768 cells | 4.1 | 0.51 | 0.06 |
//! | 48³ = 110,592 cells | 13.8 | 1.7 | 0.22 |
//! | 64³ = 262,144 cells | 32.8 | 4.1 | 0.51 |
//!
//! so on **seven of the nine** (world, density) pairs the literal ratio *cannot
//! reach 4 even if the harness invalidated every probe in the world*. Scoring C1
//! on that quotient would be M-44's vacuous zero wearing a ratio's clothes.
//!
//! So `invalidation_factor` is recorded in **probe units**: the invalidated probe
//! count over the number of probes the dilated support itself contains,
//! `brush_support_cells / s³`. Its ceiling is
//! `probes_total · s³ / brush_support_cells = world_cells / brush_support_cells`
//! — **density-independent**, and of order 90–920 as measured. The harness
//! asserts that ceiling exceeds 4 on every row, so a factor above 4 was always
//! reachable. `invalidation_factor_literal` is emitted beside it so nothing is
//! hidden, and `factor_ceiling` / `factor_ceiling_literal` put the arithmetic
//! above on every row.
//!
//! **The probe-unit factor is exactly a volume ratio**, which is why it is the
//! right one: `probes_invalidated · s³ / brush_support_cells`, and
//! `invalidated_volume_cells` is the numerator emitted directly. So "a constant
//! factor" is the claim that the invalidated *region* is a bounded multiple of
//! the support, and "constant across three probe densities" is the test of it.
//!
//! **Its conditioning is not the same at all three densities, and that has to be
//! on the record.** The support holds `brush_support_cells / s³` probes — about
//! 45 at `s = 2`, 5.6 at `s = 4` and **0.70 at `s = 8`** — so one invalidated
//! probe is worth 0.02, 0.18 and 1.43 factor units respectively, and the 4.0
//! threshold sits at 179, 22 and **2.8 probes**. At `s = 8` the factor is
//! quantised in steps larger than a third of the bound it is being tested
//! against, so an `s = 8` row cannot distinguish a wide invalidation from a
//! narrow one. `s = 2` is the density at which this clause is decided; `s = 8`
//! is reported because the registration asked for three densities, and read as
//! quantisation rather than as evidence.
//!
//! # The instrument, and why it measures the registered quantity
//!
//! **There is no Bevy renderer here** — `crates/isomesh` must not depend on Bevy
//! — so "probes whose visibility changes" is implemented as the thing a probe
//! update actually does: **a probe's gathered radiance is a function of where its
//! gather rays land**, so a probe is invalidated exactly when one of its rays
//! lands somewhere else.
//!
//! - Probes sit on a regular lattice at spacing `s` cells, offset to cell
//!   centres. A probe is *live* in a state when the field there is air
//!   (`value >= 0`, the crate's own convention), which is what a probe volume
//!   does with probes buried in geometry.
//! - Each live probe casts `DIRS` = 32 rays on a Fibonacci sphere — a fixed,
//!   shared direction set, so before and after are the same gather.
//! - A ray is sphere-traced against the field with a **global Lipschitz
//!   divisor**, because the gyroid is not a distance field and `|g|` must not be
//!   used as one. Outcome per ray: hit distance, `MISS`, `UNRESOLVED` (step
//!   budget exhausted) or `DEAD` (probe not in air).
//! - A probe is **invalidated** when any of its 32 outcomes changes kind, or
//!   moves by more than `HIT_TOL` = 0.5 cells. Where the brush does not
//!   dominate, `max(f, −sphere)` returns `f` *bit-identically*, so an unaffected
//!   ray is bit-identical before and after and the count has no numerical floor.
//!
//! That is a visibility measurement, not a lighting one: it says which probes'
//! gather results moved, which is precisely the cache-invalidation quantity C1
//! and C3 are about. It says nothing about radiance magnitude, and does not
//! claim to.
//!
//! **C2's two arms trace the same rays through the same geometry by two
//! methods**: sphere-tracing the `BrushStack` (the crate's real field, base
//! plus every logged brush, no pruning) against `parry3d`'s BVH ray-cast over
//! the mesh `surface_nets` extracts from that same field. Ray-set equality is
//! asserted, not assumed: both arms are driven from one liveness vector and one
//! direction set, both return a ray counter, and `rays_field == rays_mesh` is an
//! `assert!`. The BVH build and the extraction are reported as separate columns
//! rather than folded into the gather, because the gather is the 87–96%.
//!
//! # The world
//!
//! `cell_size` is 1, so a cell is a unit and a sample index is its own world
//! position. The field is
//!
//! ```text
//! f = max( min( min(gyroid, −box), annulus ), −ball )
//! ```
//!
//! - **gyroid** at period 12 cells: a pervasive, connected cave labyrinth whose
//!   air fraction is world-size-independent. That is load-bearing for C1's world
//!   sweep — a world that grew by adding *solid* would make "does not grow with
//!   world size" true by having nowhere for new probes to live. `probes_air_before`
//!   is on every row so the reader can check the population really grew.
//! - **−box** seals the world at `[2.25, W − 2.25]`, so rays terminate on rock
//!   rather than escaping into an infinite gyroid the extracted mesh does not
//!   contain. The face is deliberately off-lattice so the cap meshes cleanly.
//! - **annulus + ball** is a **sealed pocket**: a ball of air of radius 5 inside
//!   a solid shell 3 cells thick, at a fixed world position `[16, 16, 16]` in
//!   every world size. This is C3's second air component and it exists by
//!   construction; that it *is* a second component is asserted from the
//!   union-find, not from the geometry.
//!
//! The edit log is M-50's four buckets — 8, 23, 38, 53 subtracted spheres of
//! radius 3, at bucket midpoints — placed on a Fibonacci shell 12.5–15.5 cells
//! from the pocket, clamped into the **smallest** world's interior so the log is
//! byte-identical at every world size. Nothing in the log may reach the pocket:
//! every centre is at least 11.5 cells out, against a shell whose outer radius is
//! 8, and the union-find assertion below is what actually checks it.
//!
//! # Controls, each an assertion rather than a printed number
//!
//! - **The registered vacuity control.** `Repair::merges` is the number of
//!   *pre-existing* components the newly-air blob touched, so **1 is an ordinary
//!   widening and 0 is an isolated bubble** — only `>= 2` is M-311's join. Before
//!   the measured edit, `air.connected(pocket, main)` is **false**,
//!   `components() >= 2`, and the pocket's label holds the ball's volume. After
//!   the breakthrough edit: `merges >= 2`, `components()` has **dropped**, and
//!   `connected` is **true**. After the control edit: `merges == 1`,
//!   `components()` is **unchanged**, and `connected` is still **false**. The
//!   same three quantities separate the two arms of every pair, and the arm is
//!   *chosen* by them rather than checked afterwards — a direction that failed to
//!   merge is not the breakthrough arm.
//! - **C1's factor could have exceeded 4**: `factor_ceiling > 4` is asserted per
//!   row.
//! - **Equal ray count**: `rays_field == rays_mesh`, asserted per row.
//! - **The support is not clipped**: no output-changed cell touches the boundary
//!   of the region `mark_edit` was given, so `brush_support_cells` is the whole
//!   support and not a window on it. **This control fired on the first run** and
//!   is why `brush_support_cells` is M-314's `output_changed_cells` rather than
//!   its `value_changed_cells` — see [`Support`].
//! - **The support agrees with the crate's own instrument**: this harness builds
//!   its own per-cell bitmap (it needs one, to say which probes sit inside the
//!   support) and asserts both its popcount and its value-changed count equal
//!   M-314's `output_changed_cells` and `value_changed_cells`. Two instruments,
//!   two numbers, no slack.
//! - **Nothing is measured over an empty set**: `probes_live > 0` and
//!   `brush_support_cells > 0` per row.
//!
//! # Machine discipline
//!
//! M-280: the CPU is governed, so `cpu_mhz` is on every row and C2's verdict is
//! a **ratio**. M-281: `gather_ms_field` and `gather_ms_mesh` are taken in the
//! same loop iteration of the same build, min of `REPS` passes, and the headline
//! figure is the median over all 72 rows rather than any single pair.

#![allow(
    clippy::cast_possible_truncation,
    clippy::cast_precision_loss,
    clippy::cast_sign_loss,
    clippy::cast_possible_wrap,
    clippy::too_many_lines
)]

mod common;

use std::time::Instant;

use isomesh::brush::{Brush, BrushStack};
use isomesh::chunk::ChunkLayout;
use isomesh::chunk::dirty::{DirtySet, mark_edit};
use isomesh::collider::triangle_indices;
use isomesh::connectivity::Air;
use isomesh::fields::Sphere;
use isomesh::surface_nets::SurfaceNets;
use isomesh::{MeshBuffer, RuntimeShape3, Sdf};
use parry3d::math::Vector;
use parry3d::query::{Ray, RayCast};
use parry3d::shape::TriMesh;

/// The scalar. `f32` because `parry3d` is `f32`-only, and a `f64` field arm
/// against an `f32` mesh arm would put a precision difference inside C2's ratio.
type S = f32;

/// Worlds swept, in cells per axis. Growing at **fixed cell size**, so the edit
/// keeps its physical size and the world genuinely grows around it — which is
/// the configuration C1's second falsifier is about.
const WORLDS: [u32; 3] = [32, 48, 64];

/// Probe lattice spacings, in cells. All three divide every world.
const SPACINGS: [u32; 3] = [2, 4, 8];

/// M-50's four edit-log buckets, at their midpoints.
const BUCKETS: [(&str, usize); 4] = [("1-15", 8), ("16-30", 23), ("31-45", 38), ("46-60", 53)];

/// Cave period in cells. The gyroid scale is `2π / period`.
const CAVE_PERIOD: S = 12.0;

/// Gyroid level set. Zero is the balanced surface, and it is what this uses.
///
/// **The first fixture inset the sealing block 6.25 cells on every face, which
/// put only 22.6% of the world inside it and cut the air labyrinth into
/// interleaved fragments.** With three air components at 32³, the harness's own
/// control-arm search reported 0 of 159 candidate digs usable and the run
/// stopped. The fix was the block, not the level set: at an inset of
/// [`BLOCK_INSET`] the air inside is 2–3 components including the pocket, 64 of
/// 64 breakthrough directions merge, and 112 of the control candidates widen
/// exactly one component. So the canonical `iso = 0` is kept and no level-set
/// fudge is in the fixture.
const CAVE_ISO: S = 0.0;

/// How far the sealing block is inset from the sampling domain, in cells.
///
/// Off-lattice on purpose: a face landing exactly on a sample plane would give
/// `f == 0` there, which the crate classifies as air, and the cap would mesh
/// degenerately.
const BLOCK_INSET: S = 2.25;

/// Gather rays per probe, on a Fibonacci sphere.
///
/// A real DDGI probe casts 64–256; 32 is the smallest set that resolves a
/// radius-4 hole at the cave scale, and the direction under-resolution biases
/// `probes_invalidated` **downward** — i.e. toward C1 holding. Stated rather
/// than hidden.
const DIRS: usize = 32;

/// Sphere-trace step budget. Exhaustion is an outcome (`UNRESOLVED`), not an
/// error: it is deterministic and identical before and after wherever the field
/// is.
const MAX_STEPS: u32 = 96;

/// Surface hit threshold, in cells.
const SURF_EPS: S = 1.0e-3;

/// Floor on a trace step, so a stalled ray consumes its budget rather than
/// looping.
const MIN_STEP: S = 1.0e-3;

/// How far a hit may move before the probe is invalidated, in cells.
const HIT_TOL: S = 0.5;

/// Radius of a logged brush, in cells.
const LOG_BRUSH_R: S = 3.0;

/// Radius of the **measured** brush, in cells. The same in both arms, which is
/// what makes them "the same volume of dig" to first order; `dug_samples` on
/// both rows is the second-order check.
const EDIT_BRUSH_R: S = 4.0;

/// Sealed pocket: centre, air radius, shell thickness. Fixed world position at
/// every world size.
const POCKET: [S; 3] = [16.0, 16.0, 16.0];
/// Pocket air radius, in cells.
const POCKET_R: S = 5.0;
/// Pocket shell thickness, in cells. Three cells, so no sample-level leak.
const SHELL_T: S = 3.0;

/// Distance from the pocket centre at which the breakthrough brush is placed.
/// With `EDIT_BRUSH_R` = 4 it spans `d ∈ [4, 12]`, so it removes the whole shell
/// `[5, 8]` along its direction and reaches four cells past it.
const BT_STANDOFF: S = 8.0;

/// Timed passes per gather, min taken.
const REPS: usize = 2;

/// Ray outcome sentinels. Distinct from any hit distance, which is `>= 0`.
const MISS: S = -1.0;
/// Step budget exhausted.
const UNRESOLVED: S = -2.0;
/// The probe was not in air in this state.
const DEAD: S = -3.0;

/// The base field: gyroid caves inside a sealed block, with one sealed pocket.
#[derive(Clone, Copy, Debug)]
struct Rock {
    /// Block centre.
    center: [S; 3],
    /// Block half-extents. Everything outside the block is solid.
    half: [S; 3],
    /// Gyroid spatial frequency: the cave period is `2π / scale`.
    scale: S,
    /// Gyroid level set. Zero is the balanced surface.
    iso: S,
}

impl Rock {
    /// Global Lipschitz divisor for a sphere trace of this field.
    ///
    /// `∂g/∂x = scale·(cos a·cos b − sin c·sin a)` is bounded by `2·scale`, so
    /// `|∇g| <= 2√3·scale`. Every other operand (box, sphere, annulus) is an
    /// exact distance field with constant 1, and `min`/`max` take the max of
    /// their operands' constants.
    fn lipschitz(&self) -> S {
        (2.0 * 3.0_f32.sqrt() * self.scale).max(1.0)
    }
}

impl Sdf for Rock {
    type Scalar = S;

    #[inline]
    fn sample(&self, p: [S; 3]) -> S {
        let (a, b, c) = (self.scale * p[0], self.scale * p[1], self.scale * p[2]);
        let g = a.sin() * b.cos() + b.sin() * c.cos() + c.sin() * a.cos() - self.iso;

        // Exact box distance, negative inside; `-box` is therefore the SDF of
        // "solid everywhere outside the block", and `min` unions that solid in.
        let mut q = [0.0; 3];
        for axis in 0..3 {
            q[axis] = (p[axis] - self.center[axis]).abs() - self.half[axis];
        }
        let outside = {
            let m = [q[0].max(0.0), q[1].max(0.0), q[2].max(0.0)];
            (m[0] * m[0] + m[1] * m[1] + m[2] * m[2]).sqrt()
        };
        let inside = q[0].max(q[1]).max(q[2]).min(0.0);
        let cave = g.min(-(outside + inside));

        // The pocket: union a solid annulus, then subtract the inner ball. The
        // annulus wins over the gyroid wherever they disagree, so the shell is
        // solid whatever the cave network does there.
        let d = {
            let v = [p[0] - POCKET[0], p[1] - POCKET[1], p[2] - POCKET[2]];
            (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt()
        };
        let annulus = (d - (POCKET_R + SHELL_T)).max(POCKET_R - d);
        cave.min(annulus).max(POCKET_R - d)
    }
}

/// `n` directions on a Fibonacci sphere.
fn fibonacci_dirs(n: usize) -> Vec<[S; 3]> {
    let golden = core::f32::consts::PI * (3.0 - 5.0_f32.sqrt());
    (0..n)
        .map(|i| {
            let z = 1.0 - 2.0 * (i as S + 0.5) / n as S;
            let r = (1.0 - z * z).max(0.0).sqrt();
            let th = golden * i as S;
            [r * th.cos(), r * th.sin(), z]
        })
        .collect()
}

/// The edit log for one bucket: `n` subtracted spheres on a shell around the
/// pocket, clamped into the **smallest** world so the log is identical at every
/// world size.
fn edit_log(n: usize) -> Vec<Brush<Sphere<S>>> {
    let lo = 6.0;
    let hi = WORLDS[0] as S - 6.0;
    let mut out = Vec::with_capacity(n);
    for (i, d) in fibonacci_dirs(256).iter().enumerate() {
        if out.len() == n {
            break;
        }
        let rad = 12.5 + 1.5 * (i % 3) as S;
        let mut c = [
            POCKET[0] + d[0] * rad,
            POCKET[1] + d[1] * rad,
            POCKET[2] + d[2] * rad,
        ];
        for v in &mut c {
            *v = v.clamp(lo, hi);
        }
        // Must not reach the pocket shell: outer radius 8, brush radius 3.
        if distance(c, POCKET) < 11.5 {
            continue;
        }
        out.push(Brush::subtract(Sphere {
            center: c,
            radius: LOG_BRUSH_R,
        }));
    }
    assert!(
        out.len() == n,
        "the Fibonacci candidate set did not yield {n} usable log brushes"
    );
    out
}

fn distance(a: [S; 3], b: [S; 3]) -> S {
    let v = [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
    (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt()
}

/// Sphere-trace one ray. Returns a hit distance, [`MISS`] or [`UNRESOLVED`].
#[inline]
fn trace_one<F: Sdf<Scalar = S>>(
    field: &F,
    o: [S; 3],
    d: [S; 3],
    max_t: S,
    lip: S,
    steps: &mut u64,
) -> S {
    let mut t = 0.0;
    for _ in 0..MAX_STEPS {
        *steps += 1;
        let p = [o[0] + d[0] * t, o[1] + d[1] * t, o[2] + d[2] * t];
        let v = field.sample(p);
        if v < SURF_EPS {
            return t;
        }
        t += (v / lip).max(MIN_STEP);
        if t > max_t {
            return MISS;
        }
    }
    UNRESOLVED
}

/// Which probes are in air in this field state. `value >= 0` is the crate's own
/// air convention (`cube::is_inside` is `value < 0`).
fn liveness<F: Sdf<Scalar = S>>(field: &F, probes: &[[S; 3]]) -> Vec<bool> {
    probes.iter().map(|p| field.sample(*p) >= 0.0).collect()
}

/// Gather every live probe against the field. Returns `(rays, steps)`.
fn gather_field<F: Sdf<Scalar = S>>(
    field: &F,
    probes: &[[S; 3]],
    live: &[bool],
    dirs: &[[S; 3]],
    max_t: S,
    lip: S,
    out: &mut [S],
) -> (u64, u64) {
    let (mut rays, mut steps) = (0u64, 0u64);
    for (i, p) in probes.iter().enumerate() {
        let base = i * dirs.len();
        if !live[i] {
            out[base..base + dirs.len()].fill(DEAD);
            continue;
        }
        for (k, d) in dirs.iter().enumerate() {
            out[base + k] = trace_one(field, *p, *d, max_t, lip, &mut steps);
            rays += 1;
        }
    }
    (rays, steps)
}

/// Gather the same ray set against the extracted mesh's BVH. Returns rays cast.
fn gather_mesh(
    tri: &TriMesh,
    probes: &[[S; 3]],
    live: &[bool],
    dirs: &[[S; 3]],
    max_t: S,
    out: &mut [S],
) -> u64 {
    let mut rays = 0u64;
    for (i, p) in probes.iter().enumerate() {
        let base = i * dirs.len();
        if !live[i] {
            out[base..base + dirs.len()].fill(DEAD);
            continue;
        }
        let origin = Vector::new(p[0], p[1], p[2]);
        for (k, d) in dirs.iter().enumerate() {
            let ray = Ray::new(origin, Vector::new(d[0], d[1], d[2]));
            out[base + k] = tri.cast_local_ray(&ray, max_t, true).unwrap_or(MISS);
            rays += 1;
        }
    }
    rays
}

/// Outcome class: 0 hit, 1 miss, 2 unresolved, 3 dead.
#[inline]
fn kind(v: S) -> u8 {
    if v >= 0.0 {
        0
    } else if v > -1.5 {
        1
    } else if v > -2.5 {
        2
    } else {
        3
    }
}

/// Whether two outcomes for the same ray differ enough to invalidate the probe.
#[inline]
fn outcome_differs(a: S, b: S) -> bool {
    let (ka, kb) = (kind(a), kind(b));
    if ka != kb {
        return true;
    }
    ka == 0 && (a - b).abs() > HIT_TOL
}

/// The per-cell **output-changed** bitmap over the region an edit was measured
/// in, plus the value-changed count for the same region.
///
/// # Why `brush_support_cells` is the output-changed set and not the
/// value-changed one
///
/// The first version used `value_changed_cells`, and the harness's own
/// no-clipping control refused it: `max(f, −sphere)` moves the value wherever
/// `−sphere > f`, so a brush's *value* support extends `|f|` beyond the ball —
/// and `f` reaches −30 deep inside the sealed block, so that support is tens of
/// cells across and **grows with the world**, because the block's interior
/// distance does. A quantity that grows with world size cannot be the
/// denominator of a clause about world-proportionality; it would have hidden
/// exactly the effect C1's second falsifier is looking for.
///
/// M-314's own words for the other column: *"cells whose triangles change: this
/// is the set that genuinely needs re-meshing"*. That is the brush's geometric
/// support, it is what M-50's E1 is denominated in, and it is local to the ball
/// by construction — a surface cell twenty cells away has `f ≈ 0` there, so
/// `max(f, −16)` is `f` and nothing about it moves. `value_changed_cells` is
/// still emitted, named `value_changed_cells_in_region` because within this
/// region it is a **window** on an unbounded set rather than a measurement of
/// one.
struct Support {
    min_cell: [i64; 3],
    extent: [usize; 3],
    /// Output-changed cells: value moved **and** the cell was or is a surface
    /// cell.
    changed: Vec<bool>,
    /// Popcount of [`changed`](Self::changed).
    count: u64,
    /// Cells in this region whose corner values moved at all.
    value_changed: u64,
}

impl Support {
    /// Whether the cell containing `p` is in the support.
    fn contains_point(&self, p: [S; 3]) -> bool {
        let mut idx = [0usize; 3];
        for axis in 0..3 {
            let c = p[axis].floor() as i64 - self.min_cell[axis];
            if c < 0 || c >= self.extent[axis] as i64 {
                return false;
            }
            idx[axis] = c as usize;
        }
        self.changed[idx[0] + self.extent[0] * (idx[1] + self.extent[1] * idx[2])]
    }
}

/// Build the support bitmap with the same per-corner comparison and the same
/// surface test `mark_edit` makes. Two instruments, and the caller asserts they
/// agree.
fn support_of<A, B>(before: &A, after: &B, min_cell: [i64; 3], max_cell: [i64; 3]) -> Support
where
    A: Sdf<Scalar = S>,
    B: Sdf<Scalar = S>,
{
    let mut planes = [0usize; 3];
    for axis in 0..3 {
        planes[axis] = (max_cell[axis] - min_cell[axis] + 2) as usize;
    }
    let corners = planes[0] * planes[1] * planes[2];
    let mut moved = vec![false; corners];
    let mut in_before = vec![false; corners];
    let mut in_after = vec![false; corners];
    for z in 0..planes[2] {
        for y in 0..planes[1] {
            for x in 0..planes[0] {
                let p = [
                    (min_cell[0] + x as i64) as S,
                    (min_cell[1] + y as i64) as S,
                    (min_cell[2] + z as i64) as S,
                ];
                let (a, b) = (before.sample(p), after.sample(p));
                let i = x + planes[0] * (y + planes[1] * z);
                // Bit comparison, not `!=`: a sign flip on a zero sample is a
                // real change to the classification the crate makes.
                moved[i] = a.total_cmp(&b).is_ne();
                // `cube::is_inside` is `value < 0`, and zero is outside.
                in_before[i] = a < 0.0;
                in_after[i] = b < 0.0;
            }
        }
    }

    let extent = [planes[0] - 1, planes[1] - 1, planes[2] - 1];
    let mut changed = vec![false; extent[0] * extent[1] * extent[2]];
    let mut count = 0u64;
    let mut value_changed = 0u64;
    for z in 0..extent[2] {
        for y in 0..extent[1] {
            for x in 0..extent[0] {
                let mut any = false;
                let mut n_before = 0u32;
                let mut n_after = 0u32;
                for corner in 0..8usize {
                    let o = [corner & 1, (corner >> 1) & 1, (corner >> 2) & 1];
                    let i = (x + o[0]) + planes[0] * ((y + o[1]) + planes[1] * (z + o[2]));
                    any |= moved[i];
                    n_before += u32::from(in_before[i]);
                    n_after += u32::from(in_after[i]);
                }
                if any {
                    value_changed += 1;
                }
                let was_surface = n_before != 0 && n_before != 8;
                let is_surface = n_after != 0 && n_after != 8;
                if any && (was_surface || is_surface) {
                    changed[x + extent[0] * (y + extent[1] * z)] = true;
                    count += 1;
                }
            }
        }
    }
    Support {
        min_cell,
        extent,
        changed,
        count,
        value_changed,
    }
}

/// The samples one subtracted sphere turns from solid to air.
///
/// `after = max(before, −sphere)` is `>= 0` only where `before >= 0` already or
/// the sphere contains the point, so the flip set is inside the closed ball and
/// a one-cell margin is generous.
fn flipped_samples<A, B>(before: &A, after: &B, center: [S; 3], radius: S, w: u32) -> Vec<[u32; 3]>
where
    A: Sdf<Scalar = S>,
    B: Sdf<Scalar = S>,
{
    let mut out = Vec::new();
    let lo = |v: S| ((v - radius - 1.0).floor().max(0.0)) as u32;
    let hi = |v: S| ((v + radius + 1.0).ceil().min(w as S)) as u32;
    for z in lo(center[2])..=hi(center[2]) {
        for y in lo(center[1])..=hi(center[1]) {
            for x in lo(center[0])..=hi(center[0]) {
                let p = [x as S, y as S, z as S];
                if before.sample(p) < 0.0 && after.sample(p) >= 0.0 {
                    out.push([x, y, z]);
                }
            }
        }
    }
    out
}

/// Mean core clock, so a governed nanosecond has its clock on the row (M-280).
fn cpu_mhz() -> f64 {
    let text = std::fs::read_to_string("/proc/cpuinfo").unwrap_or_default();
    let mut sum = 0.0;
    let mut n = 0.0;
    for line in text.lines() {
        if let Some(rest) = line.strip_prefix("cpu MHz")
            && let Some((_, v)) = rest.split_once(':')
            && let Ok(v) = v.trim().parse::<f64>()
        {
            sum += v;
            n += 1.0;
        }
    }
    if n > 0.0 { sum / n } else { f64::NAN }
}

/// Everything one arm of one (world, bucket) pair needs.
struct Arm {
    breakthrough: bool,
    center: [S; 3],
    brushes: Vec<Brush<Sphere<S>>>,
    support: Support,
    /// M-314's `output_changed_cells` — the registered `brush_support_cells`.
    support_cells: u64,
    /// M-314's `value_changed_cells`, a window on an unbounded set.
    value_changed_cells: u64,
    dug: u64,
    relabels: u64,
    merges: u64,
    components_after: u64,
    connected_after: bool,
    tri: TriMesh,
    tris: usize,
    extract_ms: f64,
    bvh_ms: f64,
}

type Row = Vec<(&'static str, String)>;

fn main() {
    if !std::env::args().any(|arg| arg == "--bench") {
        return;
    }

    let prereg = isomesh::experiment!("P-78");

    let mhz = cpu_mhz();
    let dirs = fibonacci_dirs(DIRS);
    let bt_candidates = fibonacci_dirs(64);
    let ct_candidates = fibonacci_dirs(96);

    println!(
        "{:>4} {:>7} {:>4} {:>8} {:>7} {:>8} {:>7} {:>9} {:>9} {:>7}",
        "W", "bucket", "s", "support", "probes", "inval", "factor", "field_ms", "mesh_ms", "speed"
    );

    let mut rows: Vec<Row> = Vec::new();
    // C1's world clause: probes_invalidated keyed by (bucket, spacing, arm) so
    // the growth across worlds is read off one table rather than eyeballed.
    let mut world_scan: Vec<(&'static str, u32, bool, u32, u64, u64)> = Vec::new();

    for w in WORLDS {
        let rock = Rock {
            center: [w as S * 0.5; 3],
            half: [w as S * 0.5 - BLOCK_INSET; 3],
            scale: core::f32::consts::TAU / CAVE_PERIOD,
            iso: CAVE_ISO,
        };
        let lip = rock.lipschitz();
        let shape = RuntimeShape3::new([w + 1; 3]).expect("sample grid fits u32");
        let layout = ChunkLayout::<S>::new(16, 1.0, [0.0; 3]).expect("valid layout");
        let max_t = 3.0_f32.sqrt() * w as S;

        for (bucket, log_len) in BUCKETS {
            let log = edit_log(log_len);
            let before = BrushStack {
                base: rock,
                brushes: &log,
            };

            // ── the union-find, built once on the pre-edit field ─────────────
            let mut values = Vec::with_capacity(((w + 1) * (w + 1) * (w + 1)) as usize);
            for z in 0..=w {
                for y in 0..=w {
                    for x in 0..=w {
                        values.push(before.sample([x as S, y as S, z as S]));
                    }
                }
            }
            let (air, _) = Air::build(&values, &shape).expect("one value per sample");
            let components_before = air.components();

            let pocket_sample = [POCKET[0] as u32, POCKET[1] as u32, POCKET[2] as u32];
            let pocket_label = air
                .label_of(pocket_sample)
                .expect("the pocket centre is air by construction");
            let pocket_air_samples = air.component_size(pocket_label);

            // The main air network: the largest component that is not the
            // pocket, found from the union-find rather than from the geometry.
            let mut main_label = u32::MAX;
            let mut main_size = 0u32;
            for l in 0..air.label_count() as u32 {
                if l == pocket_label {
                    continue;
                }
                let size = air.component_size(l);
                if size > main_size {
                    main_size = size;
                    main_label = l;
                }
            }
            assert!(
                main_label != u32::MAX && main_size > pocket_air_samples,
                "no main air network larger than the pocket at W={w} bucket {bucket}"
            );
            let mut main_sample = None;
            'find: for z in 0..=w {
                for y in 0..=w {
                    for x in 0..=w {
                        if air.label_of([x, y, z]) == Some(main_label) {
                            main_sample = Some([x, y, z]);
                            break 'find;
                        }
                    }
                }
            }
            let main_sample = main_sample.expect("the main label has a member");

            // ── the registered vacuity control, half one ────────────────────
            assert!(
                components_before >= 2,
                "W={w} bucket {bucket}: only {components_before} air component, so the \
                 breakthrough arm has nothing to merge"
            );
            assert!(
                !air.connected(pocket_sample, main_sample),
                "W={w} bucket {bucket}: the pocket is already connected to the main air \
                 network, so the edit log leaked into it and C3 has no fixture"
            );
            let expected_pocket = 4.0 / 3.0 * core::f32::consts::PI * POCKET_R.powi(3);
            assert!(
                f64::from(pocket_air_samples) > 0.5 * f64::from(expected_pocket),
                "W={w} bucket {bucket}: the pocket holds {pocket_air_samples} air samples \
                 against a ball volume of {expected_pocket:.0}, so it is not a cavity"
            );

            // ── the breakthrough arm ────────────────────────────────────────
            //
            // `Repair::merges` is the number of **pre-existing** components the
            // newly-air blob touched (read from `grow_from`: `touched.len()`,
            // whichever way union-by-size resolves). So `1` is an ordinary
            // widening, `0` is an isolated bubble in bulk rock — which *adds* a
            // component — and only `>= 2` is M-311's join. The first fixture
            // asked for `>= 1` on this arm and `== 0` on the control, which is
            // "any dig" against "a dig that opens a new bubble": both wrong, and
            // the control search returning empty is what exposed it.
            let mut bt: Option<([S; 3], u64, u64, u64, u64, bool)> = None;
            for d in &bt_candidates {
                let c = [
                    POCKET[0] + d[0] * BT_STANDOFF,
                    POCKET[1] + d[1] * BT_STANDOFF,
                    POCKET[2] + d[2] * BT_STANDOFF,
                ];
                let mut brushes = log.clone();
                brushes.push(Brush::subtract(Sphere {
                    center: c,
                    radius: EDIT_BRUSH_R,
                }));
                let after = BrushStack {
                    base: rock,
                    brushes: &brushes,
                };
                let samples = flipped_samples(&before, &after, c, EDIT_BRUSH_R, w);
                let mut probe = air.clone();
                let repair = probe.dig(&samples, || true);
                let joined = probe.connected(pocket_sample, main_sample);
                if repair.merges >= 2 && probe.components() < components_before && joined {
                    bt = Some((
                        c,
                        repair.dirty,
                        repair.relabels,
                        repair.merges,
                        probe.components(),
                        joined,
                    ));
                    break;
                }
            }
            let (bt_center, bt_dug, bt_relabels, bt_merges, bt_components, bt_joined) =
                bt.expect("no candidate direction merged the pocket into the main air network");

            // ── the control arm: same brush radius, matched dug volume, and it
            //    must widen exactly one component — no join, no new bubble ────
            let mut ct: Option<([S; 3], u64, u64, u64, u64, bool)> = None;
            let mut best_gap = u64::MAX;
            for d in &ct_candidates {
                for rad in [14.0_f32, 16.0, 18.0] {
                    let mut c = [
                        POCKET[0] + d[0] * rad,
                        POCKET[1] + d[1] * rad,
                        POCKET[2] + d[2] * rad,
                    ];
                    for v in &mut c {
                        *v = v.clamp(6.0, WORLDS[0] as S - BLOCK_INSET - EDIT_BRUSH_R);
                    }
                    if distance(c, POCKET) < 13.0 {
                        continue;
                    }
                    let mut brushes = log.clone();
                    brushes.push(Brush::subtract(Sphere {
                        center: c,
                        radius: EDIT_BRUSH_R,
                    }));
                    let after = BrushStack {
                        base: rock,
                        brushes: &brushes,
                    };
                    let samples = flipped_samples(&before, &after, c, EDIT_BRUSH_R, w);
                    if samples.is_empty() {
                        continue;
                    }
                    let mut probe = air.clone();
                    let repair = probe.dig(&samples, || true);
                    let joined = probe.connected(pocket_sample, main_sample);
                    if repair.merges != 1 || probe.components() != components_before || joined {
                        continue;
                    }
                    let gap = bt_dug.abs_diff(repair.dirty);
                    if gap < best_gap {
                        best_gap = gap;
                        ct = Some((
                            c,
                            repair.dirty,
                            repair.relabels,
                            repair.merges,
                            probe.components(),
                            joined,
                        ));
                    }
                }
            }
            let (ct_center, ct_dug, ct_relabels, ct_merges, ct_components, ct_joined) = ct.expect(
                "no candidate dig widened exactly one air component at a comparable volume",
            );

            // ── the registered vacuity control, half two ────────────────────
            assert!(
                bt_merges >= 2 && bt_components < components_before && bt_joined,
                "VACUITY: the breakthrough arm reported {bt_merges} union-find merges, \
                 {components_before} -> {bt_components} components and joined={bt_joined}; \
                 it must absorb at least two pre-existing components and drop the count"
            );
            assert!(
                ct_merges == 1 && ct_components == components_before && !ct_joined,
                "the control arm reported {ct_merges} merges and {components_before} -> \
                 {ct_components} components, so it is not a same-volume dig that widens one \
                 component and opens nothing"
            );

            // ── per-arm geometry: support, EditReport, mesh, BVH ─────────────
            let mut arms: Vec<Arm> = Vec::new();
            for (is_bt, center, dug, relabels, merges, components_after, connected_after) in [
                (
                    true,
                    bt_center,
                    bt_dug,
                    bt_relabels,
                    bt_merges,
                    bt_components,
                    bt_joined,
                ),
                (
                    false,
                    ct_center,
                    ct_dug,
                    ct_relabels,
                    ct_merges,
                    ct_components,
                    ct_joined,
                ),
            ] {
                let mut brushes = log.clone();
                brushes.push(Brush::subtract(Sphere {
                    center,
                    radius: EDIT_BRUSH_R,
                }));
                let after = BrushStack {
                    base: rock,
                    brushes: &brushes,
                };

                // The support of `max(f, −sphere)` reaches past the ball by up
                // to |f|, which the gyroid bounds by 3, so radius + 8 is
                // generous — and the assertion below checks it was.
                let pad = EDIT_BRUSH_R + 8.0;
                let min_cell = [
                    (center[0] - pad).floor() as i64,
                    (center[1] - pad).floor() as i64,
                    (center[2] - pad).floor() as i64,
                ];
                let max_cell = [
                    (center[0] + pad).ceil() as i64,
                    (center[1] + pad).ceil() as i64,
                    (center[2] + pad).ceil() as i64,
                ];
                let support = support_of(&before, &after, min_cell, max_cell);

                let mut dirty = DirtySet::new();
                let report = mark_edit(&layout, &before, &after, min_cell, max_cell, &mut dirty)
                    .expect("the region fits the sample space");
                assert!(
                    support.count == report.output_changed_cells
                        && support.value_changed == report.value_changed_cells,
                    "the harness's own bitmap counted {} output-changed and {} value-changed \
                     cells against M-314's {} and {}",
                    support.count,
                    support.value_changed,
                    report.output_changed_cells,
                    report.value_changed_cells
                );
                assert!(
                    support.count > 0,
                    "the measured edit changed no cell's output at all"
                );

                // The support must not touch the region boundary, or
                // `brush_support_cells` is a window rather than the support.
                let mut clipped = false;
                for z in 0..support.extent[2] {
                    for y in 0..support.extent[1] {
                        for x in 0..support.extent[0] {
                            let edge = x == 0
                                || y == 0
                                || z == 0
                                || x + 1 == support.extent[0]
                                || y + 1 == support.extent[1]
                                || z + 1 == support.extent[2];
                            if edge
                                && support.changed
                                    [x + support.extent[0] * (y + support.extent[1] * z)]
                            {
                                clipped = true;
                            }
                        }
                    }
                }
                assert!(
                    !clipped,
                    "the output-changed set reaches the boundary of the measured region, so \
                     brush_support_cells is truncated"
                );

                let mut mesh = MeshBuffer::<S>::new();
                let t0 = Instant::now();
                SurfaceNets::<S>::new()
                    .extract(&after, &shape, [0.0; 3], 1.0, &mut mesh)
                    .expect("extraction");
                let extract_ms = t0.elapsed().as_secs_f64() * 1000.0;

                let vertices: Vec<Vector> = mesh
                    .positions
                    .iter()
                    .map(|p| Vector::new(p[0], p[1], p[2]))
                    .collect();
                let t1 = Instant::now();
                let tri = TriMesh::new(vertices, triangle_indices(&mesh)).expect("trimesh");
                let bvh_ms = t1.elapsed().as_secs_f64() * 1000.0;

                arms.push(Arm {
                    breakthrough: is_bt,
                    center,
                    brushes,
                    support,
                    support_cells: report.output_changed_cells,
                    value_changed_cells: report.value_changed_cells,
                    dug,
                    relabels,
                    merges,
                    components_after,
                    connected_after,
                    tri,
                    tris: mesh.triangle_count(),
                    extract_ms,
                    bvh_ms,
                });
            }

            // ── the probe sweep ─────────────────────────────────────────────
            for s in SPACINGS {
                let per_axis = w / s;
                let mut probes = Vec::with_capacity((per_axis * per_axis * per_axis) as usize);
                for k in 0..per_axis {
                    for j in 0..per_axis {
                        for i in 0..per_axis {
                            probes.push([
                                (i as S + 0.5) * s as S,
                                (j as S + 0.5) * s as S,
                                (k as S + 0.5) * s as S,
                            ]);
                        }
                    }
                }
                let cells = probes.len() * dirs.len();

                let live_before = liveness(&before, &probes);
                let mut out_before = vec![DEAD; cells];
                let mut steps_before = 0u64;
                let (_, sb) = gather_field(
                    &before,
                    &probes,
                    &live_before,
                    &dirs,
                    max_t,
                    lip,
                    &mut out_before,
                );
                steps_before += sb;
                let probes_air_before = live_before.iter().filter(|b| **b).count() as u64;

                let mut pair: Vec<(bool, Row, u64)> = Vec::new();
                for arm in &arms {
                    let after = BrushStack {
                        base: rock,
                        brushes: &arm.brushes,
                    };
                    let live_after = liveness(&after, &probes);

                    let mut out_after = vec![DEAD; cells];
                    let mut rays_field = 0u64;
                    let mut steps_field = 0u64;
                    let mut field_ms = f64::INFINITY;
                    for _ in 0..REPS {
                        let mut steps = 0u64;
                        let t = Instant::now();
                        let (r, st) = gather_field(
                            &after,
                            &probes,
                            &live_after,
                            &dirs,
                            max_t,
                            lip,
                            &mut out_after,
                        );
                        let ms = t.elapsed().as_secs_f64() * 1000.0;
                        steps += st;
                        rays_field = r;
                        steps_field = steps;
                        field_ms = field_ms.min(ms);
                    }

                    let mut out_mesh = vec![DEAD; cells];
                    let mut rays_mesh = 0u64;
                    let mut mesh_ms = f64::INFINITY;
                    for _ in 0..REPS {
                        let t = Instant::now();
                        rays_mesh = gather_mesh(
                            &arm.tri,
                            &probes,
                            &live_after,
                            &dirs,
                            max_t,
                            &mut out_mesh,
                        );
                        mesh_ms = mesh_ms.min(t.elapsed().as_secs_f64() * 1000.0);
                    }

                    // C2's equal-ray-count control, asserted rather than assumed.
                    assert!(
                        rays_field == rays_mesh,
                        "unequal ray counts: {rays_field} against the field, {rays_mesh} \
                         against the mesh"
                    );
                    assert!(rays_field > 0, "no probe was live, so nothing was gathered");

                    // Invalidation: liveness change, or any ray outcome change.
                    let mut invalidated = 0u64;
                    let mut born = 0u64;
                    let mut live_either = 0u64;
                    let mut in_support = 0u64;
                    for (i, p) in probes.iter().enumerate() {
                        let (lb, la) = (live_before[i], live_after[i]);
                        if !lb && !la {
                            continue;
                        }
                        live_either += 1;
                        if !lb && la {
                            born += 1;
                        }
                        if arm.support.contains_point(*p) {
                            in_support += 1;
                        }
                        let base = i * dirs.len();
                        let changed = lb != la
                            || (0..dirs.len()).any(|k| {
                                outcome_differs(out_before[base + k], out_after[base + k])
                            });
                        if changed {
                            invalidated += 1;
                        }
                    }
                    assert!(live_either > 0, "no probe was live in either state");

                    // C2 compares two ways of tracing **one** scene, so the two
                    // arms have to be seeing it. A mesh gather that mostly
                    // missed would be fast for the wrong reason, and a field
                    // gather that mostly ran out of steps would be slow for the
                    // wrong reason. Both are reported and both are gated.
                    let (mut hit_f, mut hit_m, mut unres_f, mut agree) = (0u64, 0u64, 0u64, 0u64);
                    let mut deltas: Vec<S> = Vec::new();
                    for (i, _) in probes.iter().enumerate() {
                        if !live_after[i] {
                            continue;
                        }
                        let base = i * dirs.len();
                        for k in 0..dirs.len() {
                            let (f, m) = (out_after[base + k], out_mesh[base + k]);
                            let (kf, km) = (kind(f) == 0, kind(m) == 0);
                            hit_f += u64::from(kf);
                            hit_m += u64::from(km);
                            unres_f += u64::from(kind(f) == 2);
                            agree += u64::from(kf == km);
                            if kf && km {
                                deltas.push((f - m).abs());
                            }
                        }
                    }
                    deltas.sort_by(S::total_cmp);
                    let median_delta = deltas.get(deltas.len() / 2).copied().unwrap_or(S::NAN);
                    let hit_frac_field = hit_f as f64 / rays_field as f64;
                    let hit_frac_mesh = hit_m as f64 / rays_mesh as f64;
                    let agree_frac = agree as f64 / rays_field as f64;
                    assert!(
                        hit_frac_field > 0.25 && hit_frac_mesh > 0.25,
                        "a gather that mostly misses is not a gather: field hits \
                         {hit_frac_field:.3}, mesh hits {hit_frac_mesh:.3}"
                    );
                    assert!(
                        agree_frac > 0.8,
                        "the two arms disagree on hit-or-miss for {:.1}% of the rays, so they \
                         are not tracing the same scene",
                        100.0 * (1.0 - agree_frac)
                    );

                    let cube = f64::from(s * s * s);
                    let expected_in_support = arm.support_cells as f64 / cube;
                    let factor = invalidated as f64 / expected_in_support;
                    let factor_literal = invalidated as f64 / arm.support_cells as f64;
                    let world_cells = u64::from(w) * u64::from(w) * u64::from(w);
                    let ceiling = world_cells as f64 / arm.support_cells as f64;
                    let ceiling_literal = probes.len() as f64 / arm.support_cells as f64;

                    // C1's instrument must be able to say "above 4".
                    assert!(
                        ceiling > 4.0,
                        "the invalidation factor could not have exceeded 4 at W={w} s={s}: \
                         ceiling {ceiling:.2}"
                    );

                    let speedup = mesh_ms / field_ms;
                    let mean_steps = steps_field as f64 / rays_field as f64;

                    println!(
                        "{w:>4} {bucket:>7} {s:>4} {:>8} {:>7} {:>8} {factor:>7.2} \
                         {field_ms:>9.2} {mesh_ms:>9.2} {speedup:>7.3}",
                        arm.support_cells, live_either, invalidated
                    );

                    let other_dug = if arm.breakthrough { ct_dug } else { bt_dug };
                    let row: Row = vec![
                        ("log_bucket", bucket.to_string()),
                        ("probe_density", format!("{:.6}", 1.0 / cube)),
                        ("world_cells", world_cells.to_string()),
                        ("brush_support_cells", arm.support_cells.to_string()),
                        ("probes_invalidated", invalidated.to_string()),
                        ("invalidation_factor", format!("{factor:.6}")),
                        ("gather_ms_field", format!("{field_ms:.4}")),
                        ("gather_ms_mesh", format!("{mesh_ms:.4}")),
                        ("gather_speedup", format!("{speedup:.6}")),
                        (
                            "breakthrough",
                            if arm.breakthrough { "yes" } else { "no" }.to_string(),
                        ),
                        ("air_components_before", components_before.to_string()),
                        ("air_components_after", arm.components_after.to_string()),
                        // filled in below, once both arms are known
                        ("probes_invalidated_breakthrough", String::new()),
                        ("c1_holds", (factor < 4.0).to_string()),
                        ("c2_holds", (speedup >= 3.0).to_string()),
                        ("c3_holds", String::new()),
                        // ── extras ────────────────────────────────────────────
                        ("world_cells_per_axis", w.to_string()),
                        ("probe_spacing_cells", s.to_string()),
                        ("probes_total", probes.len().to_string()),
                        ("probes_air_before", probes_air_before.to_string()),
                        ("probes_live", live_either.to_string()),
                        ("probes_born", born.to_string()),
                        ("probes_in_support_actual", in_support.to_string()),
                        (
                            "probes_in_support_expected",
                            format!("{expected_in_support:.4}"),
                        ),
                        (
                            "invalidation_factor_literal",
                            format!("{factor_literal:.6}"),
                        ),
                        ("factor_ceiling", format!("{ceiling:.4}")),
                        ("factor_ceiling_literal", format!("{ceiling_literal:.4}")),
                        (
                            "invalidated_volume_cells",
                            (invalidated * u64::from(s) * u64::from(s) * u64::from(s)).to_string(),
                        ),
                        ("hit_frac_field", format!("{hit_frac_field:.4}")),
                        ("hit_frac_mesh", format!("{hit_frac_mesh:.4}")),
                        ("hit_agree_frac", format!("{agree_frac:.4}")),
                        ("hit_median_delta_cells", format!("{median_delta:.4}")),
                        (
                            "unresolved_frac_field",
                            format!("{:.4}", unres_f as f64 / rays_field as f64),
                        ),
                        (
                            "value_changed_cells_in_region",
                            arm.value_changed_cells.to_string(),
                        ),
                        ("dug_samples", arm.dug.to_string()),
                        ("dug_samples_other_arm", other_dug.to_string()),
                        ("relabels", arm.relabels.to_string()),
                        ("merges", arm.merges.to_string()),
                        ("pocket_air_samples", pocket_air_samples.to_string()),
                        ("main_component_samples", main_size.to_string()),
                        ("pocket_connected_before", "false".to_string()),
                        ("pocket_connected_after", arm.connected_after.to_string()),
                        ("rays_field", rays_field.to_string()),
                        ("rays_mesh", rays_mesh.to_string()),
                        ("dirs_per_probe", DIRS.to_string()),
                        ("max_steps", MAX_STEPS.to_string()),
                        ("mean_steps_per_ray", format!("{mean_steps:.2}")),
                        ("steps_before_gather", steps_before.to_string()),
                        ("mesh_tris", arm.tris.to_string()),
                        ("extract_ms", format!("{:.3}", arm.extract_ms)),
                        ("bvh_build_ms", format!("{:.3}", arm.bvh_ms)),
                        ("log_brushes", log_len.to_string()),
                        ("edit_brush_radius_cells", format!("{EDIT_BRUSH_R:.2}")),
                        ("hit_tol_cells", format!("{HIT_TOL:.3}")),
                        (
                            "edit_center",
                            format!(
                                "{:.2}_{:.2}_{:.2}",
                                arm.center[0], arm.center[1], arm.center[2]
                            ),
                        ),
                        ("cpu_mhz", format!("{mhz:.1}")),
                    ];
                    pair.push((arm.breakthrough, row, invalidated));
                }

                let inv_bt = pair
                    .iter()
                    .find(|(b, _, _)| *b)
                    .map(|(_, _, v)| *v)
                    .expect("the breakthrough arm ran");
                let inv_ct = pair
                    .iter()
                    .find(|(b, _, _)| !*b)
                    .map(|(_, _, v)| *v)
                    .expect("the control arm ran");
                let c3 = inv_bt > inv_ct;
                for (is_bt, mut row, _) in pair {
                    for (k, v) in &mut row {
                        if *k == "probes_invalidated_breakthrough" {
                            *v = inv_bt.to_string();
                        } else if *k == "c3_holds" {
                            *v = c3.to_string();
                        }
                    }
                    row.push(("probes_invalidated_control", inv_ct.to_string()));
                    world_scan.push((bucket, s, is_bt, w, if is_bt { inv_bt } else { inv_ct }, 0));
                    rows.push(row);
                }
            }
        }
    }

    // ── C1's world clause, read off one table ────────────────────────────────
    println!("\nprobes_invalidated against world size, at fixed edit size:");
    println!(
        "{:>7} {:>4} {:>14} {:>8} {:>8} {:>8} {:>10}",
        "bucket", "s", "arm", "W=32", "W=48", "W=64", "64/32"
    );
    let mut worst_growth: f64 = 0.0;
    for (bucket, _) in BUCKETS {
        for s in SPACINGS {
            for is_bt in [true, false] {
                let at = |w: u32| {
                    world_scan
                        .iter()
                        .find(|(b, sp, bt, ww, _, _)| {
                            *b == bucket && *sp == s && *bt == is_bt && *ww == w
                        })
                        .map_or(0, |(_, _, _, _, v, _)| *v)
                };
                let (a, b, c) = (at(32), at(48), at(64));
                let ratio = if a == 0 {
                    f64::NAN
                } else {
                    c as f64 / a as f64
                };
                if ratio.is_finite() {
                    worst_growth = worst_growth.max(ratio);
                }
                println!(
                    "{bucket:>7} {s:>4} {:>14} {a:>8} {b:>8} {c:>8} {ratio:>10.3}",
                    if is_bt { "breakthrough" } else { "control" }
                );
            }
        }
    }
    println!(
        "\nworst 64/32 growth in probes_invalidated: {worst_growth:.3}x against a world grown 8x \
         in cells"
    );

    common::experiment::run(prereg, |run| {
        for row in &rows {
            run.record(row);
        }
    });
}
