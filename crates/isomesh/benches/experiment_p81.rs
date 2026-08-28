//! **P-81 - a capsule against the field, instead of a capsule against the triangles.**
//!
//! Ticket: R-081. Pre-registered before this harness existed.
//!
//! ```bash
//! cargo bench --bench experiment_p81
//! ```
//!
//! Writes `docs/experiments/p-81.csv`.
//!
//! # The source, and the citation the docs got wrong
//!
//! The docs propose SDF collision citing Liu et al. `10.1016/j.cagd.2024.102305`,
//! which is SDF-**vs**-SDF for two comparable solids: a strictly harder problem
//! that needs interval arithmetic, and interval arithmetic is a determinism
//! liability. The game's case is a small analytic proxy against one enormous
//! static field, and the paper for that is
//!
//! > Miles Macklin, Kenny Erleben, Matthias Mueller, Nuttapong Chentanez,
//! > Stefan Jeschke and Zach Corse, *Local Optimization for Robust Signed
//! > Distance Field Collision*, PACMCGIT 3(1), 2020, `10.1145/3384538`
//! > (open PDF: `mmacklin.com/sdfcontact.pdf`)
//!
//! - per-element local optimisation between an SDF isosurface and mesh elements,
//! - projected gradient descent vs Frank-Wolfe vs **golden-section search, with
//!   GSS winning on the 1-D edge problem** - which is exactly this fixture,
//!   because a capsule *is* an edge with a radius,
//! - their decisive line on 129k-triangle rigid shells: *"this mesh-based
//!   collision took approximately 15 ms per-step, compared to under 0.5 ms using
//!   SDF-based contact"*, 48-445 microseconds per timestep, CUDA on a GTX 2080 Ti.
//!
//! **Discounted at registration, and the discount is why C1's bar is 3x and not
//! 30x:** their 30x is GPU-parallel against a deep BVH over 129k triangles, and
//! this crate's baseline is a much smaller per-chunk `parry3d` `TriMesh` on one
//! CPU thread.
//!
//! # SHARE, recomputed against `P-85` before any code was written
//!
//! The registration says *"SHARE: C1 moves the query half of the collider's 45%,
//! and `P-85` will have measured how large that half is."* `P-85` has run, and
//! **its answer is that there is no query half.** From `docs/experiments/p-85.csv`,
//! the `fbm_terrain` 33 row:
//!
//! | stage | ms | share of the collider build |
//! |---|---:|---:|
//! | `handoff` (weld + `collider::readiness`) | 3.900064 | 0.812955 |
//! | `bvh` (`Bvh::from_iter`, binned) | 0.859457 | 0.179151 |
//! | `construct` (`TriMesh::new` minus the BVH) | 0.035011 | 0.007298 |
//! | `copy` (positions + `triangle_indices`) | 0.002710 | 0.000565 |
//! | **residual** | **0.000150** | **0.0000313** |
//! | total | 4.797392 | |
//!
//! Every one of those four stages is **construction**. The decomposition closes
//! to a residual of 0.00313% of the total, so the largest slice of that 45% that
//! could be a *query* is the residual: **s <= 3.13e-5**.
//!
//! Amdahl, stated before the run so it cannot be retro-fitted:
//!
//! ```text
//! collider stage : 1 / (1 - s + s/k)  ->  ceiling 1/(1 - s) = 1.0000313x
//! whole pipeline : s' = 0.45 * s      ->  ceiling             1.0000141x
//! ```
//!
//! **So C1's registered SHARE line is falsified as a premise.** Not by this
//! harness - by `P-85`, before this harness existed. `M-135`'s 45% is
//! `collider::readiness`, a build-time validity walk over a `MeshBuffer`; it is
//! paid once per chunk when the chunk is meshed, and a character controller's
//! per-frame collision query is not inside it at all. Any speedup C1 measures
//! moves the **runtime per-frame collision budget of a moving body**, which is a
//! quantity nothing in this repository has ever measured, and which is *not* the
//! 45%.
//!
//! **The clause itself is still reachable, and that is why this runs anyway.**
//! C1 is a ratio between two independently measured per-query costs, not a share
//! of a fixed total, so it is unbounded above and the 3x bar has no arithmetic
//! ceiling; and the 20 microsecond bar is an absolute that a per-query
//! measurement answers directly. What changes is what a HELD *means*: it is a
//! statement about the runtime collision budget, and it must not be quoted as
//! moving the 45%. `X51` is the standing precedent for saying this out loud
//! before running rather than after.
//!
//! # What the three instruments are
//!
//! ## C1 - two per-query costs over one pose set
//!
//! One 33^3 chunk of `fbm_terrain` on the field's own domain (`common::grid`, so
//! the grid is `P-85`'s and `M-135`'s: 32 cells per axis, `cell_size` 0.5),
//! welded exactly as `P-85` welds it, then handed to `parry3d::shape::TriMesh`.
//! [`C1_QUERIES`] capsule poses are generated once and **both arms answer the
//! same poses**:
//!
//! - **`query_us_trimesh`** - `parry3d::query::contact(capsule, trimesh)`, the
//!   shipped narrowphase, which descends the `TriMesh`'s BVH through
//!   `contact_shape_composite_shape`. The capsule is built from its two world
//!   endpoints (`Capsule::new(a, b, r)`) and both poses are `Pose::IDENTITY`, so
//!   no transform work is charged to either arm.
//! - **`query_us_field`** - [`field_contact`]: golden-section search over the
//!   capsule's segment against `FbmTerrain::sample`, then one `gradient` call at
//!   the minimiser for the normal and the first-order Euclidean rectification
//!   `d = phi / |grad phi|`. That rectification is needed because `fbm_terrain` is
//!   `y - h(x, z)`, which overstates the true distance by `1/cos(theta)` on a
//!   slope; without it the two arms would be answering different questions at the
//!   contact boundary.
//!
//! **The poses are tilted, and that is not cosmetic.** An *upright* capsule over
//! a heightfield is a degenerate 1-D problem: with `x` and `z` fixed, `phi` is
//! monotone in `y`, so the minimiser is always the bottom endpoint and GSS is
//! doing twenty evaluations to find something a comparison would have found in
//! zero. Orientations are therefore sampled over the sphere, and
//! `endpoint_minimisers` is recorded and asserted below [`ENDPOINT_BAR`] so a
//! reader can see the optimisation had an interior answer to find. The upright
//! case is a subset of the population, and it is the one C3 uses.
//!
//! ## C2 - a fixed iteration count, and how it is guaranteed structurally
//!
//! Four separate guarantees, and the last two are measurements rather than
//! claims:
//!
//! 1. **A constant trip count.** `for _ in 0..GSS_ITERATIONS` with no `break`,
//!    no `return`, no tolerance test. There is no convergence criterion in
//!    [`gss_min`] at all.
//! 2. **A branchless interval update.** Textbook GSS chooses its next interval
//!    with `if f(c) < f(d)`. Here the comparison becomes a multiplier
//!    `m = f32::from(fc <= fd)` and every one of the six state updates is
//!    `m * left + (1 - m) * right`. The *sequence* of floating-point operations
//!    is therefore identical for every input; only the values differ. `<=` rather
//!    than `<` so a tie is resolved the same way on every machine.
//! 3. **An exact evaluation count, asserted.** Every query returns the number of
//!    `Sdf::sample` calls it spent, and `gss_evals_per_query` is asserted equal
//!    to `GSS_ITERATIONS + 2` over all [`C2_QUERIES`] of them. A data-dependent
//!    early exit would show up here as an inequality, over 10^6 chances to fire.
//! 4. **A branch-miss floor, measured against a branchy control.**
//!    [`gss_min_branchy`] is the same algorithm with the `if` left in; it exists
//!    only for this control and answers no clause. `branch_misses_per_query_field`
//!    and `branch_misses_per_query_branchy` are both recorded, so "no
//!    data-dependent branching" is a number rather than an adjective.
//!
//! The arithmetic underneath is bit-portable for the same reason the crate's
//! noise is: `fields::noise` evaluates only `+ - * /` and `floor`, every one of
//! them IEEE-754 exact, and its module docs forbid `mul_add` precisely because
//! fusion changes results. This harness holds the same rule - **no `mul_add`
//! anywhere in the query path** - and adds `sqrt` and one division, both of which
//! IEEE-754 requires to be correctly rounded. Rust does not enable FP contraction,
//! so `m * a + n * c` is a multiply, a multiply and an add on both targets.
//!
//! **The cross-machine comparison, and what happens when it is absent.** A
//! record is `(contact point x, y, z, depth)` as four raw `f32` bit patterns,
//! folded FNV-1a-64 over explicit little-endian bytes so the digest does not
//! inherit the host's byte order. [`C2_BLOCKS`] blocks of [`C2_BLOCK`] queries,
//! each reseeded, plus one digest over the whole run. The same bench, source
//! unchanged, is built and run on [`PEER_HOST`] (Apple M5,
//! [`PEER_TARGET`]); its block digests are committed here as
//! [`PEER_BLOCK_DIGESTS`] and the Linux run compares its own against them.
//! `differing_contacts` is exactly zero when every block agrees; when a block
//! disagrees it is recorded as the *upper bound* `differing_blocks * block size`
//! and `differing_is_bound` says so, because a 64-bit digest localises a
//! disagreement to its block and no finer.
//!
//! [`PEER_BLOCK_DIGESTS`] is an `Option`, and `None` is scored **BLOCKED** with
//! `c2_blocker` naming what is missing. It is deliberately not a zero-filled
//! table: a stub there would make `contacts_bit_identical` a fabricated
//! measurement, which is `X35`'s failure and `P-70`'s C3 in one costume.
//!
//! **Local self-consistency is reported either way, and asserted.** The C2 arm
//! makes two independent passes over the same pose stream in the same process and
//! `locally_self_consistent` compares their digests. A query that is not
//! reproducible on one machine cannot be compared across two, so this is the
//! floor the cross-machine clause stands on and it is measured rather than
//! assumed.
//!
//! ## C3 - what a ghost contact is, because the definition *is* the clause
//!
//! The fixture is `M-106`'s: the `game_capsule_walk` path over streamed
//! `fbm_terrain`, driven until [`TARGET_CROSSINGS`] = 495 chunk-column crossings,
//! at [`WALK_SPEED`] and [`DT`]. **It is a kinematic sweep, not a controller**,
//! and that is a deviation with a reason: a controller's trajectory diverges
//! between the two arms by construction (`P-86` says so of its own C3 and had to
//! make it a replay), and a comparison of two collision methods needs the two
//! methods to be asked the *same* question. So the capsule is placed upright with
//! its foot [`REST_PENETRATION`] below the terrain height under it - a resting
//! penetration of the same order as `game_dig`'s `GROUND_PROBE` - and stepped
//! along the path. Both arms are queried at every pose.
//!
//! Around the body a 3x3 column window of chunks (two vertical layers) is meshed
//! and rebuilt whenever the body changes column. Edges are keyed by the **raw
//! bits of the two endpoint positions**, not by index, so a seam edge whose two
//! triangles come from two independently meshed chunks is recognised as one edge
//! (`M-69`'s duplicate vertices are bit-identical: both chunks interpolate the
//! same crossing from the same two corner samples). The window is one chunk -
//! eight world units - larger than the body's reach on every side, so every edge
//! the capsule can touch has both of its triangles present.
//!
//! **A contact is a ghost contact iff all three hold:**
//!
//! 1. the witness point on the triangle lies on an **edge or a vertex** of that
//!    triangle - the minimum barycentric coordinate is at or below
//!    [`FEATURE_TOL`] - rather than in its interior;
//! 2. every edge so implicated is **internal**: shared by exactly two triangles
//!    of the window mesh. This is the condition the registration says the field
//!    query cannot meet, *"because there are no internal edges"*;
//! 3. the returned contact normal lies **outside the cone spanned by the faces
//!    that meet at that feature**. Writing `d` for the dihedral between the
//!    incident faces and `a` for the smallest angle between the contact normal
//!    and any of them, a normal that is legitimately "somewhere in between the
//!    connecting face normals" (Jolt's phrase) has `a <= d`; a ghost normal has
//!    `a > d`. The condition is `a - d > NORMAL_TOL_DEG`, and
//!    `mesh_ghost_excess_max_deg` / `field_normal_excess_max_deg` report how far
//!    past it each arm went.
//!
//! **Condition 3 is stated as an excess over the dihedral rather than as a raw
//! angle, and the first version of this harness got that wrong.** Scored on the
//! raw deviation alone - "more than 5 degrees from every incident face" - the
//! rule counts two different things on the two arms. On the mesh arm it finds
//! edge-derived normals, which is the class C3 names. On the field arm it finds
//! *Marching Cubes faceting*: `mean_dihedral_deg` is 16.78 on this fixture, so a
//! smooth normal sitting exactly between two facets is 8.4 degrees from both and
//! trips a 5-degree threshold while being the most correct normal available. The
//! excess form is immune to that by construction, applies identically to both
//! arms, and is what the two `ghost_contacts_*` columns are scored on;
//! `ghost_contacts_trimesh_rawdev` and `ghost_contacts_field_rawdev` keep the
//! naive reading in the file so the difference between the two is visible rather
//! than argued.
//!
//! Condition 3 only ever *reduces* the count; `edge_witness_contacts` records the
//! population before it, so its effect is visible. `ghost_contacts_trimesh_jolt` additionally requires the edge's
//! dihedral to be at or below Jolt's default `mActiveEdgeCosThresholdAngle` of
//! **5 degrees** (`JoltPhysics/Docs/Architecture.md`: *"Whenever a body hits an
//! inactive edge, the contact normal is the face normal"*), which is the strictly
//! smaller class that Jolt v5.0.0's internal-edge-removal work and `avian#612`
//! actually fix.
//!
//! ## C3's field half is VACUOUS, and this is where that was found out
//!
//! **The registered clause cannot be scored on the field arm, and pretending
//! otherwise would be `P-70`'s C3 exactly.** The class is *"collisions against
//! internal edges between adjacent triangles"*, and the registration's own reason
//! for predicting zero is *"because there are no internal edges"*. That is a
//! definitional emptiness, not a measurement: the field query returns
//! `grad phi` at a point on the isosurface and there is no code path by which an
//! edge could enter it.
//!
//! Three instruments were built to give it teeth, and all three failed in the
//! same way - **each one has to import the mesh as its scoring reference, and
//! then what it measures is the mesh's discretisation error rather than anything
//! the field did.** They are all in the CSV because their *disagreement* is the
//! evidence:
//!
//! | column | rule | reading |
//! |---|---|---|
//! | `ghost_contacts_field_projected_rawdev` | project the field's contact point onto the mesh, count normals >5 deg from every incident face | 1,600 |
//! | `ghost_contacts_field_projected` | the same, as an excess over the local dihedral | 373 |
//! | `ghost_contacts_field_at_mesh_witness` | evaluate `grad phi` at the mesh arm's *own* ghost witness point and apply the identical cone test there | 14,090 |
//!
//! Thirty-eight times apart, over one population, from one field. A quantity
//! whose three honest measurements disagree by 38x is not being measured; the
//! mesh's facet normals are, and `mean_dihedral_deg` = 16.78 is why. So
//! `ghost_contacts_field` is recorded as **VACUOUS** with `c3_field_verdict`
//! naming the reason, the three surrogates are kept beside it, and the half of C3
//! that *is* measurable - the `TriMesh` arm - carries its own verdict in
//! `c3_trimesh_holds`.
//!
//! **What a successor id should register instead.** The intrinsic difference
//! between an edge-derived normal and a surface normal is not their angle to any
//! facet, it is their *stability*: an edge normal is a function of where the body
//! is, so it swings by tens of degrees when the body moves a millimetre, while a
//! surface normal moves with the surface's curvature. Degrees of normal change
//! per centimetre of body travel is reference-free, applies identically to both
//! arms, and is the gameplay symptom itself. It is a different question from the
//! one `P-81` registered, so it belongs to a new id rather than to a quiet
//! substitution here.
//!
//! The pairing that produced the 14,090 is still in the harness and is still the
//! tightest of the three, so it is described here. So `ghost_contacts_field` is scored at **the mesh
//! arm's own ghost witness point**: same point, same triangle, same implicated
//! edges, same incident faces, same cone test, and the *only* thing that differs
//! is where the normal came from - parry's from the internal edge, the field's
//! from `Sdf::gradient` at that point. `ghost_poses_trimesh` is the denominator
//! and is asserted non-zero; `fg.internal` is asserted so the two arms cannot
//! drift onto different features; and `field_ghost_excess_max_deg` reports the
//! largest excess the field's gradient produced anywhere in that population, so
//! a zero comes with its margin.
//!
//! **A weaker reading is kept in the file with its mechanism named.** The first
//! version of this instrument scored the field at *its own* contact point after
//! projecting that point onto the mesh, which is `ghost_contacts_field_projected`.
//! That projection displaces the query by up to
//! `field_projection_dist_max`, and over that distance on `fbm_terrain` the true
//! surface normal rotates by more than a facet cone is wide - so the count is
//! dominated by the displacement rather than by anything the field query did. It
//! is reported because the difference between it and the paired count is the
//! measurement of that confound.
//!
//! # The vacuity control
//!
//! Registered: *"the `TriMesh` arm must report a non-zero ghost-contact count, or
//! C3 cannot fire."* `ghost_contacts_trimesh > 0` is asserted. Three further
//! controls exist because that one alone does not prove the fixture reached the
//! configurations it names:
//!
//! - `seam_crossings == 495` - `M-106`'s fixture actually delivered.
//! - `internal_edges > 0` and `boundary_edges` recorded - the position-keyed edge
//!   map really did join triangles across chunk seams. A window whose every edge
//!   read as a boundary would make condition 2 unsatisfiable and C3's `TriMesh`
//!   count a guaranteed zero, which is the *other* vacuity and the registration
//!   does not name it.
//! - `ghost_poses_trimesh > 0` - the field arm's zero is over a non-empty
//!   population.
//!
//! And for C1: `contact_disagreements` and `contact_fraction_*` are recorded, and
//! both arms are asserted to find contacts on a non-trivial fraction of the pose
//! set. Two query paths that both answer "no contact" everywhere would produce a
//! fast, meaningless ratio.
//!
//! # Units
//!
//! `M-280`: on a governed CPU a nanosecond is not a unit. Every row carries `ghz`
//! measured over the same span as the cycles, plus `cycles_per_query_*`, which is
//! clock-independent; `speedup` is a ratio taken inside one repeat of one binary
//! (`M-281`). The absolute microsecond figures move with the governor and C1's
//! 20-microsecond bar is read from them, so the clock is on the row and the
//! per-repeat spread is too.

#![allow(
    clippy::cast_possible_truncation,
    clippy::cast_precision_loss,
    clippy::cast_sign_loss,
    clippy::too_many_lines,
    clippy::similar_names,
    // `INV_PHI` is written to more digits than `f32` can hold on purpose: the
    // literal must be unambiguous about which real number it names, because C2's
    // whole claim is that the same constant lands on the same `f32` on two
    // targets. Truncating it to `0.618_034` happens to round the same way today
    // and says nothing about why.
    clippy::excessive_precision
)]

mod common;

use std::collections::HashMap;
use std::hint::black_box;
use std::time::Instant;

use parry3d::math::{Pose, Vector};
use parry3d::query::{PointQuery, contact};
use parry3d::shape::{Capsule, TriMesh, Triangle};

use isomesh::chunk::{ChunkId, ChunkLayout};
use isomesh::collider;
use isomesh::fields::{FbmTerrain, ReferenceField};
use isomesh::marching_cubes::MarchingCubes;
use isomesh::weld::{self, Welder};
use isomesh::{MeshBuffer, RuntimeShape3, Sdf};

use crate::common::grid;

#[cfg(target_os = "linux")]
use crate::common::counters::{MIN_TIME_RATIO, Probe};

// ── the body ─────────────────────────────────────────────────────────────────

/// Capsule radius. `game_capsule_walk`'s and `game_dig`'s value.
const CAPSULE_RADIUS: f32 = 0.4;
/// Half the capsule's segment, so total height is `0.9 + 2 * 0.4 = 1.7`.
const CAPSULE_HALF: f32 = 0.45;

// ── the field query ──────────────────────────────────────────────────────────

/// Golden-section iterations, fixed and with no convergence test.
///
/// The interval shrinks by `1/phi` per iteration, so twenty gives
/// `0.618^20 = 4.0e-5` of the segment - `3.6e-5` world units on a 0.9-unit
/// segment, which is one fourteen-thousandth of a `cell_size` of 0.5. Precision
/// is not what fixes the count; determinism is. See the module docs' four
/// structural guarantees.
const GSS_ITERATIONS: u32 = 20;

/// `1/phi = (sqrt(5) - 1)/2`, to more digits than `f32` can hold, so the literal
/// rounds to the same `f32` on every target.
const INV_PHI: f32 = 0.618_033_988_749_895;
/// `1 - 1/phi = 1/phi^2`. A const expression rather than a second literal, so the
/// two cannot drift apart.
const INV_PHI2: f32 = 1.0 - INV_PHI;

/// What one field query answered.
#[derive(Clone, Copy)]
struct FieldContact {
    /// The point on the isosurface, world space. C2's payload.
    point: [f32; 3],
    /// Unit outward normal there, `grad phi / |grad phi|`.
    normal: [f32; 3],
    /// `radius - phi/|grad phi|`. Non-negative means contact.
    depth: f32,
    /// The minimiser's parameter along the capsule segment.
    t: f32,
    /// `Sdf::sample` calls spent. Always `GSS_ITERATIONS + 2`.
    evals: u32,
}

/// Golden-section search for the minimum of `phi` along `a..b`, branchless.
///
/// Macklin et al. 2020's 1-D edge problem, which is the case their comparison
/// gives to GSS. Returns the minimiser's parameter and the field value there.
///
/// The interval reuse is the standard one - `1/phi^2 + 1/phi = 1`, so one of the
/// two interior points of the new interval is always one of the old two - and
/// costs exactly one `sample` per iteration. The choice of which is which is a
/// multiplier, not a branch: `m` is 1.0 to keep the left interval and 0.0 to keep
/// the right, and each of the six state updates evaluates both candidates and
/// blends. No `mul_add`, for the reason `fields::noise`'s module docs give.
fn gss_min<F: Sdf<Scalar = f32>>(
    field: &F,
    a: [f32; 3],
    b: [f32; 3],
    evals: &mut u32,
) -> (f32, f32) {
    let d = [b[0] - a[0], b[1] - a[1], b[2] - a[2]];
    let at = |t: f32| [a[0] + t * d[0], a[1] + t * d[1], a[2] + t * d[2]];

    let mut lo = 0.0_f32;
    let mut hi = 1.0_f32;
    let mut h = hi - lo;
    let mut c = lo + INV_PHI2 * h;
    let mut e = lo + INV_PHI * h;
    let mut fc = field.sample(at(c));
    let mut fe = field.sample(at(e));
    *evals += 2;

    for _ in 0..GSS_ITERATIONS {
        let m = f32::from(fc <= fe);
        let n = 1.0 - m;
        let nlo = m * lo + n * c;
        let nhi = m * e + n * hi;
        h *= INV_PHI;
        // Left: the fresh point is the new interval's lower interior point.
        // Right: it is the upper one. Both are computed, one is selected.
        let x = m * (nlo + INV_PHI2 * h) + n * (nhi - INV_PHI2 * h);
        let fx = field.sample(at(x));
        *evals += 1;
        let nc = m * x + n * e;
        let ne = m * c + n * x;
        let nfc = m * fx + n * fe;
        let nfe = m * fc + n * fx;
        lo = nlo;
        hi = nhi;
        c = nc;
        e = ne;
        fc = nfc;
        fe = nfe;
    }

    let m = f32::from(fc <= fe);
    let n = 1.0 - m;
    (m * c + n * e, m * fc + n * fe)
}

/// The same search with the `if` left in, for the branch-miss control only.
///
/// **This answers no clause.** It exists so that "no data-dependent branching"
/// in [`gss_min`] is a measured difference in `branch_misses` rather than a
/// property asserted about source code. Its arithmetic is the textbook update,
/// which is why its results may differ from [`gss_min`]'s in the last bits: a
/// blend `1.0 * a + 0.0 * c` is `a + 0.0`, and that is `a` for every `a` except
/// `-0.0`. `branchy_digest` records what it produced.
fn gss_min_branchy<F: Sdf<Scalar = f32>>(field: &F, a: [f32; 3], b: [f32; 3]) -> (f32, f32) {
    let d = [b[0] - a[0], b[1] - a[1], b[2] - a[2]];
    let at = |t: f32| [a[0] + t * d[0], a[1] + t * d[1], a[2] + t * d[2]];

    let mut lo = 0.0_f32;
    let mut hi = 1.0_f32;
    let mut h = hi - lo;
    let mut c = lo + INV_PHI2 * h;
    let mut e = lo + INV_PHI * h;
    let mut fc = field.sample(at(c));
    let mut fe = field.sample(at(e));

    for _ in 0..GSS_ITERATIONS {
        h *= INV_PHI;
        if fc <= fe {
            hi = e;
            e = c;
            fe = fc;
            c = lo + INV_PHI2 * h;
            fc = field.sample(at(c));
        } else {
            lo = c;
            c = e;
            fc = fe;
            e = hi - INV_PHI2 * h;
            fe = field.sample(at(e));
        }
    }

    if fc <= fe { (c, fc) } else { (e, fe) }
}

/// One capsule-vs-field contact query, by GSS.
///
/// `a` and `b` are the capsule segment's world endpoints. The rectification
/// `d = phi / |grad phi|` is the first-order Euclidean distance: `fbm_terrain` is
/// `y - h(x, z)`, which overstates the true distance by `1/cos(theta)` on a
/// slope, and without it this arm and the `TriMesh` arm disagree about where the
/// surface is. The gradient is one extra `fbm` evaluation and is needed anyway,
/// because Macklin's contact normal *is* the SDF gradient at the minimiser.
///
/// The minimiser of `phi` is not exactly the minimiser of `phi/|grad phi|`; the
/// difference is second order in the slope's variation along the segment and is
/// not corrected, which is stated rather than hidden.
fn field_contact<F: Sdf<Scalar = f32>>(field: &F, a: [f32; 3], b: [f32; 3]) -> FieldContact {
    let mut evals = 0;
    let (t, phi) = gss_min(field, a, b, &mut evals);
    let d = [b[0] - a[0], b[1] - a[1], b[2] - a[2]];
    let p = [a[0] + t * d[0], a[1] + t * d[1], a[2] + t * d[2]];
    let g = field.gradient(p);
    let gl = (g[0] * g[0] + g[1] * g[1] + g[2] * g[2]).sqrt();
    let normal = [g[0] / gl, g[1] / gl, g[2] / gl];
    let euclid = phi / gl;
    FieldContact {
        point: [
            p[0] - normal[0] * euclid,
            p[1] - normal[1] * euclid,
            p[2] - normal[2] * euclid,
        ],
        normal,
        depth: CAPSULE_RADIUS - euclid,
        t,
        evals,
    }
}

// ── poses ────────────────────────────────────────────────────────────────────

/// A capsule, as its two segment endpoints in world space.
#[derive(Clone, Copy)]
struct Segment {
    a: [f32; 3],
    b: [f32; 3],
}

impl Segment {
    fn capsule(&self) -> Capsule {
        Capsule::new(
            Vector::new(self.a[0], self.a[1], self.a[2]),
            Vector::new(self.b[0], self.b[1], self.b[2]),
            CAPSULE_RADIUS,
        )
    }
}

/// xorshift64. Integer state and one IEEE division to `f32`, so the stream and
/// every pose derived from it are bit-identical on both targets.
struct Rng(u64);

impl Rng {
    fn next(&mut self) -> f32 {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        ((self.0 >> 40) as f32) / 16_777_216.0
    }

    /// A vector in `[-1, 1]^3`, by rejection, then normalised.
    ///
    /// Rejection is a data-dependent loop, and that is fine: it is *pose
    /// generation*, not the query, and the fixed-iteration claim is about the
    /// query. The stream is the same on both machines because the comparisons are
    /// on the same IEEE values.
    fn direction(&mut self) -> [f32; 3] {
        loop {
            let v = [
                2.0 * self.next() - 1.0,
                2.0 * self.next() - 1.0,
                2.0 * self.next() - 1.0,
            ];
            let l2 = v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
            if l2 > 0.01 && l2 <= 1.0 {
                let l = l2.sqrt();
                return [v[0] / l, v[1] / l, v[2] / l];
            }
        }
    }
}

/// The terrain height under `(x, z)`.
///
/// Exact rather than sampled: `phi(p) = p.y - h(x, z)`, so `h = -phi([x, 0, z])`.
fn height(field: &FbmTerrain<f32>, x: f32, z: f32) -> f32 {
    -field.sample([x, 0.0, z])
}

/// How far above or below the resting height a generated pose sits, in world
/// units. Symmetric, so roughly half the poses penetrate.
const POSE_BAND: f32 = 0.35;

/// One tilted capsule near the surface, inside `lo..hi` on `x` and `z`.
fn pose(field: &FbmTerrain<f32>, rng: &mut Rng, lo: f32, span: f32) -> Segment {
    let x = lo + rng.next() * span;
    let z = lo + rng.next() * span;
    let axis = rng.direction();
    let y =
        height(field, x, z) + CAPSULE_HALF + CAPSULE_RADIUS + (2.0 * rng.next() - 1.0) * POSE_BAND;
    let c = [x, y, z];
    Segment {
        a: [
            c[0] - axis[0] * CAPSULE_HALF,
            c[1] - axis[1] * CAPSULE_HALF,
            c[2] - axis[2] * CAPSULE_HALF,
        ],
        b: [
            c[0] + axis[0] * CAPSULE_HALF,
            c[1] + axis[1] * CAPSULE_HALF,
            c[2] + axis[2] * CAPSULE_HALF,
        ],
    }
}

// ── digests ──────────────────────────────────────────────────────────────────

/// FNV-1a-64 offset basis.
const FNV_OFFSET: u64 = 0xcbf2_9ce4_8422_2325;
/// FNV-1a-64 prime.
const FNV_PRIME: u64 = 0x0000_0100_0000_01b3;

/// Fold one 32-bit word into an FNV-1a-64 digest, little-endian.
fn fold(h: u64, word: u32) -> u64 {
    let mut h = h;
    for byte in word.to_le_bytes() {
        h ^= u64::from(byte);
        h = h.wrapping_mul(FNV_PRIME);
    }
    h
}

/// Fold a whole contact record: three coordinates and the depth, raw bits.
fn fold_contact(h: u64, c: &FieldContact) -> u64 {
    let mut h = h;
    for v in [c.point[0], c.point[1], c.point[2], c.depth] {
        h = fold(h, v.to_bits());
    }
    h
}

// ── hardware counters ────────────────────────────────────────────────────────

/// Cycles, branch misses and the wall time they were taken over.
struct Hw {
    cycles: u64,
    branch_misses: u64,
    ns: u128,
}

/// Run `f` with the hardware counters open.
///
/// `perf_event_open` is a Linux system call with no macOS equivalent a bench can
/// reach, which is why `common::counters` is `cfg`-gated at its module
/// declaration. The non-Linux arm measures the same span with the same
/// `Instant`; only the counter columns go missing, and they are recorded as
/// `unavailable` rather than invented (the convention `common/mod.rs` names).
#[cfg(target_os = "linux")]
fn counted<T>(f: impl FnOnce() -> T) -> (T, Option<Hw>) {
    let mut probe = Probe::open();
    let start = Instant::now();
    probe.reset_and_enable();
    let out = f();
    probe.disable();
    let ns = start.elapsed().as_nanos();
    let counts = probe.read();
    assert!(
        counts.worst_ratio() >= MIN_TIME_RATIO,
        "a counter was multiplexed ({:.4}), so `ghz` would be an extrapolation",
        counts.worst_ratio()
    );
    (
        out,
        Some(Hw {
            cycles: counts.cycles.count,
            branch_misses: counts.branch_misses.count,
            ns,
        }),
    )
}

#[cfg(not(target_os = "linux"))]
fn counted<T>(f: impl FnOnce() -> T) -> (T, Option<Hw>) {
    let start = Instant::now();
    let out = f();
    black_box(start.elapsed());
    (out, None)
}

// ── C1 ───────────────────────────────────────────────────────────────────────

/// Poses per C1 arm.
///
/// Twenty thousand: enough that one repeat of either arm is tens of
/// milliseconds, so `Instant`'s resolution is irrelevant, and small enough that
/// nine interleaved repeats of both arms cost a couple of seconds.
const C1_QUERIES: usize = 20_000;
/// Repeats discarded before timing.
const WARMUP: usize = 3;
/// Timed repeats per arm, median taken.
///
/// Nine, and interleaved arm-inside-repeat, for `M-281`: a machine excursion then
/// lands on both arms of the same repeat and the ratio survives it. `M-337`'s
/// re-audit - a registered 1.25x floor that re-measured at 1.022 - is the
/// standing warning against a single reading on a governed CPU.
const REPS: usize = 9;
/// C1's ratio bar.
const SPEEDUP_BAR: f64 = 3.0;
/// C1's absolute bar, microseconds per query.
const QUERY_BUDGET_US: f64 = 20.0;
/// Contact prediction distance handed to `parry3d::query::contact`. Zero, so
/// both arms answer "is there an overlap" and nothing wider.
const PREDICTION: f32 = 0.0;
/// The registered resolution, plus one corroboration arm.
const RESOLUTIONS: [u32; 2] = [33, 65];
/// How many of the pose set's minimisers may sit at a segment endpoint before
/// the 1-D optimisation is measuring a comparison rather than a search.
const ENDPOINT_BAR: f64 = 0.90;
/// How close to an endpoint counts as at one.
const ENDPOINT_TOL: f32 = 0.02;
/// The smallest share of poses either arm must find in contact, so a fast ratio
/// cannot come from two arms that both reject everything.
const CONTACT_FLOOR: f64 = 0.05;

/// One C1 arm's numbers.
struct C1 {
    samples: u32,
    triangles: usize,
    cell_size: f32,
    field_us: f64,
    trimesh_us: f64,
    speedup: f64,
    speedup_min: f64,
    speedup_max: f64,
    field_us_median: f64,
    trimesh_us_median: f64,
    contact_share_field: f64,
    contact_share_trimesh: f64,
    disagreements: u64,
    endpoint_share: f64,
    /// `None` where `perf_event_open` does not exist, so the column reads
    /// `unavailable` rather than `0`. A plausible-but-false `0.0000` GHz on a row
    /// is exactly the cell `P-64` is remembered for.
    cycles_field: Option<u64>,
    cycles_trimesh: Option<u64>,
    ghz: Option<f64>,
}

/// Contour, weld and hand one chunk to `parry3d`, exactly as `P-85` does.
fn trimesh_at(field: &FbmTerrain<f32>, samples: u32) -> (TriMesh, f32, usize, [f32; 3]) {
    let (shape, lo, cell_size): (RuntimeShape3, [f32; 3], f32) = grid(field, samples);
    let mut mesh = MeshBuffer::<f32>::new();
    MarchingCubes::<f32>::new()
        .extract(field, &shape, lo, cell_size, &mut mesh)
        .expect("extraction on a reference field's own grid");
    Welder::<f32>::new()
        .weld(&mut mesh, weld::epsilon_for(cell_size))
        .expect("the welder accepts a mesh Marching Cubes just produced");
    let triangles = mesh.triangle_count();
    let vertices: Vec<Vector> = mesh
        .positions
        .iter()
        .map(|p| Vector::new(p[0], p[1], p[2]))
        .collect();
    let indices = collider::triangle_indices(&mesh);
    let trimesh = TriMesh::new(vertices, indices).expect("a non-empty index buffer");
    (trimesh, cell_size, triangles, lo)
}

fn run_c1(field: &FbmTerrain<f32>, samples: u32) -> C1 {
    let (trimesh, cell_size, triangles, lo) = trimesh_at(field, samples);
    let span = cell_size * (samples - 1) as f32;
    // Keep the whole capsule plus its radius inside the meshed box, so no pose
    // is answered by an absent triangle.
    let margin = CAPSULE_HALF + CAPSULE_RADIUS + cell_size;
    let mut rng = Rng(0x0081_5EED_1234_ABCD);
    let poses: Vec<Segment> = (0..C1_QUERIES)
        .map(|_| pose(field, &mut rng, lo[0] + margin, span - 2.0 * margin))
        .collect();
    let capsules: Vec<Capsule> = poses.iter().map(Segment::capsule).collect();

    // Agreement, endpoint and contact-share controls, outside the timed loops.
    let mut contact_field = 0_u64;
    let mut contact_trimesh = 0_u64;
    let mut disagreements = 0_u64;
    let mut endpoints = 0_u64;
    for (seg, capsule) in poses.iter().zip(&capsules) {
        let f = field_contact(field, seg.a, seg.b);
        assert_eq!(
            f.evals,
            GSS_ITERATIONS + 2,
            "the field query spent a data-dependent number of samples"
        );
        let hit_f = f.depth >= 0.0;
        let hit_t = contact(
            &Pose::IDENTITY,
            capsule,
            &Pose::IDENTITY,
            &trimesh,
            PREDICTION,
        )
        .expect("parry supports capsule-vs-trimesh contact")
        .is_some();
        contact_field += u64::from(hit_f);
        contact_trimesh += u64::from(hit_t);
        disagreements += u64::from(hit_f != hit_t);
        endpoints += u64::from(f.t < ENDPOINT_TOL || f.t > 1.0 - ENDPOINT_TOL);
    }

    let n = C1_QUERIES as f64;
    let mut reps: Vec<(f64, f64)> = Vec::with_capacity(REPS);
    for rep in 0..(WARMUP + REPS) {
        let t = Instant::now();
        for seg in &poses {
            black_box(field_contact(field, seg.a, seg.b));
        }
        let field_ns = t.elapsed().as_nanos() as f64;
        let t = Instant::now();
        for capsule in &capsules {
            black_box(
                contact(
                    &Pose::IDENTITY,
                    capsule,
                    &Pose::IDENTITY,
                    &trimesh,
                    PREDICTION,
                )
                .expect("parry supports capsule-vs-trimesh contact"),
            );
        }
        let trimesh_ns = t.elapsed().as_nanos() as f64;
        if rep >= WARMUP {
            reps.push((field_ns / n / 1000.0, trimesh_ns / n / 1000.0));
        }
    }

    // The reported repeat is the one whose *speedup* is the median, so the row's
    // three numbers come from one repeat and `speedup` reproduces from them.
    let mut order: Vec<usize> = (0..reps.len()).collect();
    order.sort_by(|&i, &j| {
        (reps[i].1 / reps[i].0)
            .partial_cmp(&(reps[j].1 / reps[j].0))
            .expect("no NaN in a wall time")
    });
    let mid = order[reps.len() / 2];
    let mut field_sorted: Vec<f64> = reps.iter().map(|r| r.0).collect();
    let mut tri_sorted: Vec<f64> = reps.iter().map(|r| r.1).collect();
    field_sorted.sort_by(|a, b| a.partial_cmp(b).expect("no NaN"));
    tri_sorted.sort_by(|a, b| a.partial_cmp(b).expect("no NaN"));
    let ratios: Vec<f64> = reps.iter().map(|r| r.1 / r.0).collect();

    let (_, hw_field) = counted(|| {
        for seg in &poses {
            black_box(field_contact(field, seg.a, seg.b));
        }
    });
    let (_, hw_trimesh) = counted(|| {
        for capsule in &capsules {
            black_box(
                contact(
                    &Pose::IDENTITY,
                    capsule,
                    &Pose::IDENTITY,
                    &trimesh,
                    PREDICTION,
                )
                .expect("parry supports capsule-vs-trimesh contact"),
            );
        }
    });

    C1 {
        samples,
        triangles,
        cell_size,
        field_us: reps[mid].0,
        trimesh_us: reps[mid].1,
        speedup: reps[mid].1 / reps[mid].0,
        speedup_min: ratios.iter().copied().fold(f64::INFINITY, f64::min),
        speedup_max: ratios.iter().copied().fold(f64::NEG_INFINITY, f64::max),
        field_us_median: field_sorted[field_sorted.len() / 2],
        trimesh_us_median: tri_sorted[tri_sorted.len() / 2],
        contact_share_field: contact_field as f64 / n,
        contact_share_trimesh: contact_trimesh as f64 / n,
        disagreements,
        endpoint_share: endpoints as f64 / n,
        cycles_field: hw_field.as_ref().map(|h| h.cycles),
        cycles_trimesh: hw_trimesh.as_ref().map(|h| h.cycles),
        ghz: hw_field.as_ref().map(|h| h.cycles as f64 / h.ns as f64),
    }
}

// ── C2 ───────────────────────────────────────────────────────────────────────

/// The registered 10^6.
const C2_QUERIES: u64 = 1_000_000;
/// Digest blocks. Sixty-four, so [`PEER_BLOCK_DIGESTS`] is eight readable lines
/// and a disagreement is localised to 15,625 queries rather than to the whole
/// run.
const C2_BLOCKS: usize = 64;
/// Queries per block. `C2_QUERIES / C2_BLOCKS`, exact by construction.
const C2_BLOCK: u64 = C2_QUERIES / C2_BLOCKS as u64;

/// The cross-machine peer.
const PEER_HOST: &str = "mac_air";
/// The peer's target triple.
const PEER_TARGET: &str = "aarch64-apple-darwin";
/// The peer's toolchain, so the comparison names the compiler as well as the CPU.
const PEER_RUSTC: &str = "1.96.1-31fca3adb";

/// The M5's [`C2_BLOCKS`] block digests over the same 10^6 contact records, or
/// `None` if the peer run has not been performed.
///
/// Measured 2026-08-27 on `mac_air`: MacBook Air, Apple M5, `aarch64-apple-darwin`,
/// rustc 1.96.1 (31fca3adb 2026-06-26), cargo 1.96.1, macOS 26.5.2. The peer built
/// this file **unmodified** - no macOS patch was needed, because `counted`'s
/// `cfg(not(target_os = "linux"))` arm and the gated counters import already cover
/// it - and its source md5 was checked equal on both machines immediately before
/// the peer build. The peer printed `evals=22000000`, so the two pose streams
/// agree and the digest comparison is valid rather than a comparison of two
/// different workloads. Raw peer output with its provenance header lives in
/// `docs/experiments/p-81-m5-digests.txt`.
///
/// **`None` is a verdict, not a placeholder.** A zero-filled table standing in
/// for a measurement would make `contacts_bit_identical` a fabricated answer -
/// `X35`'s failure and `P-70`'s C3 in one - so the absent case is carried in the
/// type and scored `BLOCKED` in the CSV, naming what is missing. When the peer
/// run lands this becomes `Some([...])` with its host, target, toolchain and date
/// in this comment, and C2 is scored on it.
///
/// Peer: [`PEER_HOST`] (MacBook Air, Apple M5, [`PEER_TARGET`], rustc
/// [`PEER_RUSTC`]). The peer runs *this* binary's source unchanged and prints the
/// digests on its `P-81 c2 block digests` line; nothing else crosses the wire.
const PEER_BLOCK_DIGESTS: Option<[u64; C2_BLOCKS]> = Some([
    0x60d3_40e2_a801_5da9,
    0xc5c5_fa16_ffc3_33e8,
    0x963d_90b2_1578_1b8a,
    0x92a9_a88e_ccc9_b92e,
    0x4334_2364_a40a_27ab,
    0x4095_bf5d_dbd3_84bb,
    0x97e3_e09b_a596_3e82,
    0xf078_cedc_7f89_ea35,
    0x04ce_8b83_041e_3644,
    0x78bb_a7d4_cfc0_2986,
    0xa649_a268_e9cd_9a54,
    0x9a74_f1cb_2204_1ab0,
    0x5459_32da_d1f9_a55f,
    0xe5ab_781e_8f80_8a2c,
    0x8e9a_d6a8_092a_ffb0,
    0x69bd_eb66_e1f9_774f,
    0x66b2_7ee5_9734_dae9,
    0x72fc_c009_f31f_a845,
    0xfa33_cd50_9ff5_8a57,
    0x99ef_61fe_d2b8_4a42,
    0x9b0c_4fa2_c613_31ae,
    0x92ba_97da_6b20_c2a8,
    0x9ede_b048_1efd_8f90,
    0x37c4_3a83_30ae_7c5c,
    0x0bdd_a647_65f4_2af8,
    0xf0e2_0223_4d83_0656,
    0x7022_3801_1fd9_498d,
    0xb805_3270_dd9b_3483,
    0xd02f_2f82_d4ba_f999,
    0x28b5_d8d9_f175_1f04,
    0xa423_b451_cafa_3c20,
    0xc5c9_b8d8_46f4_6222,
    0x6823_251e_9ad2_eda2,
    0x42fc_44a8_1163_ed39,
    0x9ec8_f2ba_f072_506f,
    0xf49b_bd06_5418_e291,
    0x94bb_242f_20b4_cf4f,
    0xfd55_4fe6_933d_7e91,
    0xbb4a_0566_d840_23d7,
    0x5061_0a00_e569_fc10,
    0x353c_723b_1f04_e2ec,
    0xd269_7ada_f5e3_6609,
    0x8e91_1118_1f2f_f89b,
    0xa184_6c1b_719e_342e,
    0xf632_10fe_cc4b_5124,
    0xd9e2_fe23_a325_768f,
    0x559e_5635_67f8_6275,
    0x15b6_abfa_be11_1c06,
    0x0976_69bd_293b_9b32,
    0x7a80_2bc7_3570_6309,
    0x22c3_8891_da97_c55a,
    0x4b66_03eb_52ba_ef90,
    0x3a88_38c3_430a_c52b,
    0xa30c_530a_eed7_3919,
    0x713b_71e2_48cd_0bca,
    0x44e9_1947_efe9_ba3b,
    0xab10_58c0_3aea_6287,
    0x3eea_1314_228a_83fa,
    0x38c4_1f70_1985_1d7e,
    0xda6b_afe6_6885_b951,
    0x8764_04f4_9693_1fa3,
    0xd961_fb0e_408c_ae3b,
    0xfd42_e863_cce5_f1ca,
    0x688a_93e1_43b8_0c65,
]);

/// What the C2 arm measured.
struct C2 {
    /// Per-block digests over the local run's contact records.
    blocks: Vec<u64>,
    /// One digest over the whole run.
    digest: u64,
    /// The branchy control's digest over the same poses.
    branchy_digest: u64,
    /// `Sdf::sample` calls, summed. Asserted `== C2_QUERIES * (ITER + 2)`.
    evals: u64,
    /// Branch misses per query, branchless and branchy. `None` off Linux.
    branch_misses_field: Option<f64>,
    branch_misses_branchy: Option<f64>,
    /// Contacts found, so the arm is not digesting a million rejections.
    contacts: u64,
    /// A second, independent pass over the same pose stream in the same process.
    ///
    /// This is the *local* self-consistency the cross-machine clause needs
    /// underneath it: if the query is not even reproducible on one machine, a
    /// peer comparison measures nothing. Reported whether or not a peer exists.
    repeat_digest: u64,
}

fn run_c2(field: &FbmTerrain<f32>) -> C2 {
    // A wide box so the pose set is not one hill: 256 world units, which at
    // `frequency` 0.25 is 64 periods of the coarsest octave.
    const SPAN: f32 = 256.0;
    let mut blocks = Vec::with_capacity(C2_BLOCKS);
    let mut digest = FNV_OFFSET;
    let mut evals = 0_u64;
    let mut contacts = 0_u64;
    let mut rng = Rng(0x0081_C2C2_5EED_0001);

    // The local repeat, first, so `digest` and `repeat_digest` are two passes of
    // the same code over the same stream rather than one pass reused.
    let mut repeat_digest = FNV_OFFSET;
    {
        let mut rng = Rng(0x0081_C2C2_5EED_0001);
        for _ in 0..C2_QUERIES {
            let seg = pose(field, &mut rng, -SPAN * 0.5, SPAN);
            repeat_digest = fold_contact(repeat_digest, &field_contact(field, seg.a, seg.b));
        }
    }

    for _ in 0..C2_BLOCKS {
        let mut block = FNV_OFFSET;
        for _ in 0..C2_BLOCK {
            let seg = pose(field, &mut rng, -SPAN * 0.5, SPAN);
            let c = field_contact(field, seg.a, seg.b);
            assert_eq!(
                c.evals,
                GSS_ITERATIONS + 2,
                "the field query spent a data-dependent number of samples"
            );
            evals += u64::from(c.evals);
            contacts += u64::from(c.depth >= 0.0);
            block = fold_contact(block, &c);
            digest = fold_contact(digest, &c);
        }
        blocks.push(block);
    }

    // The branchy control, over the same pose stream, and the two branch-miss
    // readings. Both loops do one field query per pose and differ only in how
    // the interval is chosen.
    let mut branchy_digest = FNV_OFFSET;
    let mut rng = Rng(0x0081_C2C2_5EED_0001);
    let poses: Vec<Segment> = (0..C2_BLOCK)
        .map(|_| pose(field, &mut rng, -SPAN * 0.5, SPAN))
        .collect();
    let (_, hw_field) = counted(|| {
        for seg in &poses {
            black_box(field_contact(field, seg.a, seg.b));
        }
    });
    let (_, hw_branchy) = counted(|| {
        for seg in &poses {
            black_box(gss_min_branchy(field, seg.a, seg.b));
        }
    });
    for seg in &poses {
        let (t, phi) = gss_min_branchy(field, seg.a, seg.b);
        branchy_digest = fold(branchy_digest, t.to_bits());
        branchy_digest = fold(branchy_digest, phi.to_bits());
    }
    let per = C2_BLOCK as f64;

    C2 {
        blocks,
        digest,
        branchy_digest,
        evals,
        branch_misses_field: hw_field.as_ref().map(|h| h.branch_misses as f64 / per),
        branch_misses_branchy: hw_branchy.as_ref().map(|h| h.branch_misses as f64 / per),
        contacts,
        repeat_digest,
    }
}

// ── C3: the world ────────────────────────────────────────────────────────────

/// Cells per chunk axis. `game_capsule_walk`'s value, and `P-86`'s.
const CHUNK_CELLS: u32 = 16;
/// World units per cell.
const CELL_SIZE: f32 = 0.5;
/// The two vertical layers `fbm_terrain`'s sheet straddles: its height bound is
/// `2 * fbm_bound(4, 0.5)` = 3.75 and the layer boundary is `y = 0`.
const VERTICAL_LAYERS: [i32; 2] = [-1, 0];
/// Commanded speed. `game_capsule_walk`'s value.
const WALK_SPEED: f32 = 7.0;
/// Fixed step, so the fixture is reproducible.
const DT: f32 = 1.0 / 60.0;
/// `M-106`'s crossing count, and the registration's.
const TARGET_CROSSINGS: u64 = 495;
/// Hard cap, so a fixture that stops crossing seams fails loudly.
const MAX_STEPS: u64 = 400_000;
/// How far the capsule's foot sits below the terrain height under it. Half of
/// `game_dig`'s `GROUND_PROBE` of 0.06: a resting penetration a ground-probing
/// controller genuinely holds.
const REST_PENETRATION: f32 = 0.03;
/// Broadphase bucket side, world units. Two cells, so a Marching Cubes triangle
/// lands in one or two buckets.
const BUCKET: f32 = 1.0;
/// Barycentric tolerance below which a witness point counts as on the edge
/// opposite that coordinate rather than in the face.
const FEATURE_TOL: f64 = 1.0e-3;
/// How far a contact normal must sit from every incident face normal to be an
/// edge-derived normal, in degrees.
const NORMAL_TOL_DEG: f64 = 5.0;
/// A second, looser reading of the same tolerance, so the count is not an
/// artefact of one threshold.
const NORMAL_TOL_WIDE_DEG: f64 = 15.0;
/// Jolt's default `mActiveEdgeCosThresholdAngle`, in degrees
/// (`JoltPhysics/Docs/Architecture.md`: "By default this is 5 degrees"). An edge
/// whose dihedral is at or below this is *inactive*, and a contact on an inactive
/// edge is the class Jolt v5.0.0 and `avian#612` replace with the face normal.
const JOLT_ACTIVE_EDGE_DEG: f64 = 5.0;

/// One world triangle.
#[derive(Clone, Copy)]
struct Tri {
    v: [[f32; 3]; 3],
    lo: [f32; 3],
    hi: [f32; 3],
    /// Unit outward normal, from the winding. `MeshSink::triangle`'s convention
    /// is counter-clockwise seen from *outside* the solid, so
    /// `cross(v1 - v0, v2 - v0)` points away from the rock.
    n: [f64; 3],
}

/// A vertex key: the raw bits of its three coordinates.
///
/// Two chunks meeting at a seam interpolate the same crossing from the same two
/// corner samples, so their seam vertices are bit-identical and this key joins
/// them. Keying by index would make every seam edge look like a boundary and
/// would make C3's condition 2 unsatisfiable there - which is the vacuity the
/// registration does not name and `internal_edges` exists to rule out.
type VKey = [u32; 3];

fn vkey(p: [f32; 3]) -> VKey {
    [p[0].to_bits(), p[1].to_bits(), p[2].to_bits()]
}

/// The triangles sharing one edge.
#[derive(Clone, Copy, Default)]
struct EdgeRec {
    tris: [u32; 2],
    count: u32,
}

/// The meshed window around the body, and its edge topology.
struct Window {
    tris: Vec<Tri>,
    buckets: HashMap<[i32; 3], Vec<u32>>,
    edges: HashMap<[VKey; 2], EdgeRec>,
}

fn edge_key(a: [f32; 3], b: [f32; 3]) -> [VKey; 2] {
    let (ka, kb) = (vkey(a), vkey(b));
    if ka <= kb { [ka, kb] } else { [kb, ka] }
}

fn cross(u: [f64; 3], w: [f64; 3]) -> [f64; 3] {
    [
        u[1] * w[2] - u[2] * w[1],
        u[2] * w[0] - u[0] * w[2],
        u[0] * w[1] - u[1] * w[0],
    ]
}

fn dot(u: [f64; 3], w: [f64; 3]) -> f64 {
    u[0] * w[0] + u[1] * w[1] + u[2] * w[2]
}

fn sub(u: [f64; 3], w: [f64; 3]) -> [f64; 3] {
    [u[0] - w[0], u[1] - w[1], u[2] - w[2]]
}

impl Window {
    /// Mesh the 3x3 column neighbourhood of `column`, two layers deep.
    ///
    /// One chunk - eight world units - of margin on every side, against a body
    /// whose reach is under two, so every edge the capsule can touch has both of
    /// its triangles in the window and `internal` means internal.
    fn build(field: &FbmTerrain<f32>, layout: &ChunkLayout<f32>, column: [i32; 2]) -> Self {
        let mut w = Self {
            tris: Vec::new(),
            buckets: HashMap::new(),
            edges: HashMap::new(),
        };
        let shape = layout.sample_shape().expect("a valid sample shape");
        let mut mesh = MeshBuffer::<f32>::new();
        let mut mc = MarchingCubes::<f32>::new();
        for dx in -1..=1 {
            for dz in -1..=1 {
                for layer in VERTICAL_LAYERS {
                    let id = ChunkId {
                        coords: [column[0] + dx, layer, column[1] + dz],
                    };
                    let origin = layout.sample_origin(id);
                    mesh.reset();
                    mc.extract(field, &shape, origin, CELL_SIZE, &mut mesh)
                        .expect("marching cubes over a valid chunk");
                    w.absorb(&mesh);
                }
            }
        }
        w
    }

    fn absorb(&mut self, mesh: &MeshBuffer<f32>) {
        for f in mesh.indices.as_chunks::<3>().0 {
            let v = [
                mesh.positions[f[0] as usize],
                mesh.positions[f[1] as usize],
                mesh.positions[f[2] as usize],
            ];
            let a = v[0].map(f64::from);
            let b = v[1].map(f64::from);
            let c = v[2].map(f64::from);
            let raw = cross(sub(b, a), sub(c, a));
            let len = dot(raw, raw).sqrt();
            if len <= 0.0 {
                // A degenerate triangle has no normal and no edges worth a
                // dihedral. Skipping it is not a fallback: it is not a face.
                continue;
            }
            let n = [raw[0] / len, raw[1] / len, raw[2] / len];
            let mut lo = v[0];
            let mut hi = v[0];
            for p in &v[1..] {
                for axis in 0..3 {
                    lo[axis] = lo[axis].min(p[axis]);
                    hi[axis] = hi[axis].max(p[axis]);
                }
            }
            let index = self.tris.len() as u32;
            self.tris.push(Tri { v, lo, hi, n });
            for (i, j) in [(1_usize, 2_usize), (2, 0), (0, 1)] {
                let rec = self.edges.entry(edge_key(v[i], v[j])).or_default();
                if rec.count < 2 {
                    rec.tris[rec.count as usize] = index;
                }
                rec.count += 1;
            }
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
    }

    /// Every triangle whose bounds meet the box `lo..hi`.
    fn gather(&self, lo: [f32; 3], hi: [f32; 3], out: &mut Vec<u32>) {
        out.clear();
        let key_lo = [0, 1, 2].map(|axis| (lo[axis] / BUCKET).floor() as i32);
        let key_hi = [0, 1, 2].map(|axis| (hi[axis] / BUCKET).floor() as i32);
        for bx in key_lo[0]..=key_hi[0] {
            for by in key_lo[1]..=key_hi[1] {
                for bz in key_lo[2]..=key_hi[2] {
                    if let Some(list) = self.buckets.get(&[bx, by, bz]) {
                        out.extend_from_slice(list);
                    }
                }
            }
        }
        out.sort_unstable();
        out.dedup();
        out.retain(|&i| {
            let t = &self.tris[i as usize];
            (0..3).all(|axis| t.lo[axis] <= hi[axis] && t.hi[axis] >= lo[axis])
        });
    }
}

/// What classifying one contact against the mesh found.
struct Feature {
    /// The witness point is on at least one edge of its triangle.
    on_edge: bool,
    /// Every implicated edge is shared by exactly two triangles.
    internal: bool,
    /// The largest dihedral among the implicated edges, degrees.
    dihedral_deg: f64,
    /// The smallest angle between the contact normal and any incident face
    /// normal, degrees.
    deviation_deg: f64,
}

impl Feature {
    /// How far outside the incident faces' cone the contact normal points.
    ///
    /// Non-positive for a normal that is "somewhere in between the connecting
    /// face normals", which is what an interpolated normal on a faceted surface
    /// looks like and is not a ghost. Positive only for a normal no face at the
    /// feature points anywhere near - the edge-derived normal C3 is about.
    fn cone_excess_deg(&self) -> f64 {
        self.deviation_deg - self.dihedral_deg
    }
}

/// Classify one contact: which feature of `tri` the witness lies on, whether it
/// is internal, and how far the normal is from the faces that meet there.
///
/// This function *is* C3's clause. See the module docs' three conditions.
fn classify(window: &Window, tri_index: u32, witness: [f32; 3], normal: [f32; 3]) -> Feature {
    let tri = &window.tris[tri_index as usize];
    let a = tri.v[0].map(f64::from);
    let b = tri.v[1].map(f64::from);
    let c = tri.v[2].map(f64::from);
    let w = witness.map(f64::from);
    let raw = cross(sub(b, a), sub(c, a));
    let denom = dot(raw, raw);
    let bary = [
        dot(cross(sub(b, w), sub(c, w)), raw) / denom,
        dot(cross(sub(c, w), sub(a, w)), raw) / denom,
        dot(cross(sub(a, w), sub(b, w)), raw) / denom,
    ];

    let nrm = [
        f64::from(normal[0]),
        f64::from(normal[1]),
        f64::from(normal[2]),
    ];
    let nl = dot(nrm, nrm).sqrt();
    let nrm = [nrm[0] / nl, nrm[1] / nl, nrm[2] / nl];

    // Start from the triangle the witness is on: its own face normal is always a
    // candidate for "the normal a face would have given".
    let mut best_cos = dot(nrm, tri.n);
    let mut on_edge = false;
    let mut internal = true;
    let mut dihedral: f64 = 0.0;

    for (k, (i, j)) in [(1_usize, 2_usize), (2, 0), (0, 1)].into_iter().enumerate() {
        if bary[k] > FEATURE_TOL {
            continue;
        }
        on_edge = true;
        let rec = window
            .edges
            .get(&edge_key(tri.v[i], tri.v[j]))
            .copied()
            .unwrap_or_default();
        if rec.count != 2 {
            internal = false;
            continue;
        }
        let other = if rec.tris[0] == tri_index {
            rec.tris[1]
        } else {
            rec.tris[0]
        };
        let on = window.tris[other as usize].n;
        best_cos = best_cos.max(dot(nrm, on));
        dihedral = dihedral.max(dot(tri.n, on).clamp(-1.0, 1.0).acos().to_degrees());
    }

    Feature {
        on_edge,
        internal: on_edge && internal,
        dihedral_deg: dihedral,
        deviation_deg: best_cos.clamp(-1.0, 1.0).acos().to_degrees(),
    }
}

/// `game_capsule_walk::path`: where the capsule is asked to be after `d` metres.
/// A diagonal with a wobble, so the walk crosses seams on both axes and meets
/// three-chunk corners. Transcribed so the fixture is `M-106`'s and not a new one.
fn path(d: f32) -> (f32, f32) {
    (d * 0.82, d * 0.57 + 9.0 * (d * 0.05).sin())
}

/// What the C3 arm measured.
struct C3 {
    steps: u64,
    seam_crossings: u64,
    contacts_trimesh: u64,
    edge_witness_contacts: u64,
    internal_edge_contacts: u64,
    ghost_trimesh: u64,
    ghost_trimesh_wide: u64,
    ghost_trimesh_jolt: u64,
    ghost_trimesh_rawdev: u64,
    ghost_poses: u64,
    ghost_field: u64,
    ghost_field_projected: u64,
    ghost_field_projected_rawdev: u64,
    field_at_internal_edge: u64,
    field_ghost_excess_max_deg: f64,
    field_projection_dist_max: f64,
    contacts_field: u64,
    field_dev_max_deg: f64,
    field_excess_max_deg: f64,
    mesh_ghost_dev_max_deg: f64,
    mesh_ghost_excess_max_deg: f64,
    internal_edges: u64,
    boundary_edges: u64,
    nonmanifold_edges: u64,
    inactive_edge_share: f64,
    mean_dihedral_deg: f64,
    windows: u64,
    outward_share: f64,
}

fn run_c3(field: &FbmTerrain<f32>) -> C3 {
    let layout = ChunkLayout::<f32>::new(CHUNK_CELLS, CELL_SIZE, [0.0; 3]).expect("a valid layout");
    let mut out = C3 {
        steps: 0,
        seam_crossings: 0,
        contacts_trimesh: 0,
        edge_witness_contacts: 0,
        internal_edge_contacts: 0,
        ghost_trimesh: 0,
        ghost_trimesh_wide: 0,
        ghost_trimesh_jolt: 0,
        ghost_trimesh_rawdev: 0,
        ghost_poses: 0,
        ghost_field: 0,
        ghost_field_projected: 0,
        ghost_field_projected_rawdev: 0,
        field_at_internal_edge: 0,
        field_ghost_excess_max_deg: f64::NEG_INFINITY,
        field_projection_dist_max: 0.0,
        contacts_field: 0,
        field_dev_max_deg: 0.0,
        field_excess_max_deg: 0.0,
        mesh_ghost_dev_max_deg: 0.0,
        mesh_ghost_excess_max_deg: 0.0,
        internal_edges: 0,
        boundary_edges: 0,
        nonmanifold_edges: 0,
        inactive_edge_share: 0.0,
        mean_dihedral_deg: 0.0,
        windows: 0,
        outward_share: 0.0,
    };

    let mut column = {
        let (x, z) = path(0.0);
        let id = layout.chunk_of([x, 0.0, z]);
        [id.coords[0], id.coords[2]]
    };
    let mut window = Window::build(field, &layout, column);
    let mut edge_stats = (0_u64, 0_u64, 0_u64, 0_u64, 0.0_f64, 0_u64);
    let mut outward = (0_u64, 0_u64);
    let mut candidates: Vec<u32> = Vec::new();
    let mut distance = 0.0_f32;

    let tally = |w: &Window, s: &mut (u64, u64, u64, u64, f64, u64), o: &mut (u64, u64)| {
        for rec in w.edges.values() {
            match rec.count {
                1 => s.1 += 1,
                2 => {
                    s.0 += 1;
                    let d = dot(
                        w.tris[rec.tris[0] as usize].n,
                        w.tris[rec.tris[1] as usize].n,
                    )
                    .clamp(-1.0, 1.0)
                    .acos()
                    .to_degrees();
                    s.4 += d;
                    s.5 += 1;
                    if d <= JOLT_ACTIVE_EDGE_DEG {
                        s.3 += 1;
                    }
                }
                _ => s.2 += 1,
            }
        }
        // The winding control: a face normal from `cross(v1 - v0, v2 - v0)` must
        // point the same way as the field's own gradient, or every dihedral and
        // every deviation in this harness has its sign inverted.
        for t in &w.tris {
            let centroid = [0, 1, 2].map(|k| (t.v[0][k] + t.v[1][k] + t.v[2][k]) / 3.0);
            let g = field.gradient(centroid);
            let gl = (g[0] * g[0] + g[1] * g[1] + g[2] * g[2]).sqrt();
            let agree = t.n[0] * f64::from(g[0] / gl)
                + t.n[1] * f64::from(g[1] / gl)
                + t.n[2] * f64::from(g[2] / gl);
            o.1 += 1;
            o.0 += u64::from(agree > 0.0);
        }
    };
    tally(&window, &mut edge_stats, &mut outward);
    out.windows += 1;

    while out.seam_crossings < TARGET_CROSSINGS {
        assert!(
            out.steps < MAX_STEPS,
            "VOID: {MAX_STEPS} steps reached only {} of {TARGET_CROSSINGS} seam crossings, so the \
             fixture never delivered M-106's crossing count",
            out.seam_crossings
        );
        distance += WALK_SPEED * DT;
        let (x, z) = path(distance);
        let here = layout.chunk_of([x, 0.0, z]);
        if [here.coords[0], here.coords[2]] != column {
            out.seam_crossings += 1;
            column = [here.coords[0], here.coords[2]];
            window = Window::build(field, &layout, column);
            tally(&window, &mut edge_stats, &mut outward);
            out.windows += 1;
        }
        out.steps += 1;

        // The body: upright, foot `REST_PENETRATION` into the rock under it.
        let foot = height(field, x, z) - REST_PENETRATION;
        let centre = [x, foot + CAPSULE_RADIUS + CAPSULE_HALF, z];
        let seg = Segment {
            a: [centre[0], centre[1] - CAPSULE_HALF, centre[2]],
            b: [centre[0], centre[1] + CAPSULE_HALF, centre[2]],
        };
        let capsule = seg.capsule();
        let reach = CAPSULE_RADIUS + 0.01;
        let lo = [
            centre[0] - reach,
            centre[1] - CAPSULE_HALF - reach,
            centre[2] - reach,
        ];
        let hi = [
            centre[0] + reach,
            centre[1] + CAPSULE_HALF + reach,
            centre[2] + reach,
        ];
        window.gather(lo, hi, &mut candidates);

        // The mesh arm's manifold: one contact per overlapping triangle, which is
        // what a trimesh narrowphase hands the solver and what Jolt's
        // internal-edge removal operates on.
        let mut ghost_here: Option<(Feature, u32, [f32; 3])> = None;
        for &index in &candidates {
            let tri = &window.tris[index as usize];
            let target = Triangle::new(
                Vector::new(tri.v[0][0], tri.v[0][1], tri.v[0][2]),
                Vector::new(tri.v[1][0], tri.v[1][1], tri.v[1][2]),
                Vector::new(tri.v[2][0], tri.v[2][1], tri.v[2][2]),
            );
            let Some(c) = contact(
                &Pose::IDENTITY,
                &capsule,
                &Pose::IDENTITY,
                &target,
                PREDICTION,
            )
            .expect("parry supports capsule-vs-triangle contact") else {
                continue;
            };
            out.contacts_trimesh += 1;
            let f = classify(&window, index, c.point2.to_array(), c.normal2.to_array());
            out.edge_witness_contacts += u64::from(f.on_edge);
            out.internal_edge_contacts += u64::from(f.internal);
            if f.internal && f.deviation_deg > NORMAL_TOL_DEG {
                out.ghost_trimesh_rawdev += 1;
            }
            if f.internal && f.cone_excess_deg() > NORMAL_TOL_DEG {
                out.ghost_trimesh += 1;
                out.mesh_ghost_dev_max_deg = out.mesh_ghost_dev_max_deg.max(f.deviation_deg);
                out.mesh_ghost_excess_max_deg =
                    out.mesh_ghost_excess_max_deg.max(f.cone_excess_deg());
                if f.dihedral_deg <= JOLT_ACTIVE_EDGE_DEG {
                    out.ghost_trimesh_jolt += 1;
                }
                // The looser reading of condition 3. `NORMAL_TOL_WIDE_DEG` is
                // strictly larger, so this class is a subset and the two counts
                // together say whether the verdict is threshold-fragile.
                if f.cone_excess_deg() > NORMAL_TOL_WIDE_DEG {
                    out.ghost_trimesh_wide += 1;
                }
                if ghost_here
                    .as_ref()
                    .is_none_or(|(g, _, _)| f.cone_excess_deg() > g.cone_excess_deg())
                {
                    ghost_here = Some((f, index, c.point2.to_array()));
                }
            }
        }

        // The field arm, asked at every pose, and scored by the same rule at the
        // poses where the mesh arm ghosted.
        let fc = field_contact(field, seg.a, seg.b);
        if fc.depth >= 0.0 {
            out.contacts_field += 1;
        }
        if let Some((_, ghost_tri, ghost_witness)) = ghost_here {
            out.ghost_poses += 1;

            // **C3's answer.** Same point, same feature, same incident faces; the
            // only thing that changes is where the normal came from. parry's came
            // from an internal edge and left the faces' cone; this one is the
            // field's own gradient at that exact point. Nothing here can differ
            // between the two arms except the normal, which is the whole clause.
            let g = field.gradient(ghost_witness);
            let gl = (g[0] * g[0] + g[1] * g[1] + g[2] * g[2]).sqrt();
            let fg = classify(
                &window,
                ghost_tri,
                ghost_witness,
                [g[0] / gl, g[1] / gl, g[2] / gl],
            );
            assert!(
                fg.internal,
                "the paired field test was handed a feature the mesh arm did not call internal, \
                 so the two arms are not being scored at the same place"
            );
            out.field_ghost_excess_max_deg =
                out.field_ghost_excess_max_deg.max(fg.cone_excess_deg());
            if fg.cone_excess_deg() > NORMAL_TOL_DEG {
                out.ghost_field += 1;
            }

            // A second, weaker reading kept for the mechanism it exposes rather
            // than for the verdict: score the field's own contact point after
            // projecting it onto the mesh. That projection moves the query up to
            // half a cell, and over half a cell on `fbm_terrain` the true normal
            // rotates by more than a facet cone is wide - so this count is
            // dominated by the displacement, not by anything the field did.
            // `field_projection_dist_max` is the evidence for that sentence.
            let mut best: Option<(f64, u32, [f32; 3])> = None;
            window.gather(
                [
                    fc.point[0] - CELL_SIZE,
                    fc.point[1] - CELL_SIZE,
                    fc.point[2] - CELL_SIZE,
                ],
                [
                    fc.point[0] + CELL_SIZE,
                    fc.point[1] + CELL_SIZE,
                    fc.point[2] + CELL_SIZE,
                ],
                &mut candidates,
            );
            for &index in &candidates {
                let tri = &window.tris[index as usize];
                let target = Triangle::new(
                    Vector::new(tri.v[0][0], tri.v[0][1], tri.v[0][2]),
                    Vector::new(tri.v[1][0], tri.v[1][1], tri.v[1][2]),
                    Vector::new(tri.v[2][0], tri.v[2][1], tri.v[2][2]),
                );
                let probe = Vector::new(fc.point[0], fc.point[1], fc.point[2]);
                let w = target.project_local_point(probe, false);
                let d = f64::from((w.point - probe).length());
                if best.as_ref().is_none_or(|&(bd, _, _)| d < bd) {
                    best = Some((d, index, w.point.to_array()));
                }
            }
            if let Some((dist, index, w)) = best {
                out.field_projection_dist_max = out.field_projection_dist_max.max(dist);
                let f = classify(&window, index, w, fc.normal);
                out.field_dev_max_deg = out.field_dev_max_deg.max(f.deviation_deg);
                out.field_excess_max_deg = out.field_excess_max_deg.max(f.cone_excess_deg());
                out.field_at_internal_edge += u64::from(f.internal);
                if f.internal && f.deviation_deg > NORMAL_TOL_DEG {
                    out.ghost_field_projected_rawdev += 1;
                }
                if f.internal && f.cone_excess_deg() > NORMAL_TOL_DEG {
                    out.ghost_field_projected += 1;
                }
            }
        }
    }

    out.internal_edges = edge_stats.0;
    out.boundary_edges = edge_stats.1;
    out.nonmanifold_edges = edge_stats.2;
    out.inactive_edge_share = if edge_stats.5 == 0 {
        0.0
    } else {
        edge_stats.3 as f64 / edge_stats.5 as f64
    };
    out.mean_dihedral_deg = if edge_stats.5 == 0 {
        0.0
    } else {
        edge_stats.4 / edge_stats.5 as f64
    };
    out.outward_share = if outward.1 == 0 {
        0.0
    } else {
        outward.0 as f64 / outward.1 as f64
    };
    out
}

// ── the run ──────────────────────────────────────────────────────────────────

type Row = Vec<(&'static str, String)>;

fn main() {
    if !std::env::args().any(|a| a == "--bench") {
        return;
    }

    let field = FbmTerrain::<f32>::canonical();

    // C2 first, and printed immediately, because it is the cross-machine clause
    // and the peer runs this same binary.
    let c2 = run_c2(&field);
    assert_eq!(
        c2.evals,
        C2_QUERIES * u64::from(GSS_ITERATIONS + 2),
        "the field query's sample count is not a constant multiple of the query count, so the \
         fixed-iteration property C2 rests on does not hold"
    );
    assert!(
        c2.contacts > 0,
        "VOID: none of {C2_QUERIES} C2 poses touched the field, so the digest is a digest of \
         rejections"
    );
    print!("P-81 c2 block digests:");
    for b in &c2.blocks {
        print!(" {b:#018x}");
    }
    println!();
    println!("P-81 c2 digest: {:#018x}", c2.digest);
    println!("P-81 c2 branchy digest: {:#018x}", c2.branchy_digest);
    println!(
        "P-81 c2 arch: {} rustc-target-independent evals={}",
        std::env::consts::ARCH,
        c2.evals
    );

    // `None` means the peer run has not been performed, and the columns say
    // `BLOCKED` rather than inventing an answer.
    let peer = PEER_BLOCK_DIGESTS.map(|peer| {
        let differing_blocks = c2
            .blocks
            .iter()
            .zip(peer.iter())
            .filter(|(a, b)| a != b)
            .count() as u64;
        (
            differing_blocks == 0,
            differing_blocks,
            differing_blocks * C2_BLOCK,
        )
    });
    let blocked = || String::from("BLOCKED");
    // The word `common/mod.rs` names for a counter this platform does not have.
    // Never a zero: see `C1::cycles_field`.
    let unavailable = || String::from("unavailable");
    let locally_self_consistent = c2.digest == c2.repeat_digest;
    assert!(
        locally_self_consistent,
        "the field query is not reproducible within one process, so a cross-machine \
         comparison would be meaningless: {:#018x} then {:#018x}",
        c2.digest, c2.repeat_digest
    );

    let c3 = run_c3(&field);
    assert_eq!(
        c3.seam_crossings, TARGET_CROSSINGS,
        "the C3 fixture did not deliver M-106's 495 crossings"
    );
    assert!(
        c3.internal_edges > 0,
        "VOID: the window mesh has no internal edges at all, so C3's condition 2 could never be \
         met and the TriMesh count is a guaranteed zero for the wrong reason"
    );
    assert!(
        c3.outward_share > 0.99,
        "the winding normal disagrees with the field gradient on {:.4} of triangles, so every \
         dihedral and deviation in this harness has the wrong sign",
        1.0 - c3.outward_share
    );
    // The registered vacuity control.
    assert!(
        c3.ghost_trimesh > 0,
        "VACUITY: the TriMesh arm reported zero ghost contacts over {} contacts at {} internal \
         edges, so C3 cannot fire",
        c3.contacts_trimesh,
        c3.internal_edge_contacts
    );
    assert!(
        c3.ghost_poses > 0,
        "VACUITY: the field arm's zero would be over an empty population"
    );

    let c1: Vec<C1> = RESOLUTIONS.iter().map(|&s| run_c1(&field, s)).collect();
    for arm in &c1 {
        assert!(
            arm.contact_share_field >= CONTACT_FLOOR && arm.contact_share_trimesh >= CONTACT_FLOOR,
            "VOID: the {}^3 pose set is in contact on {:.4}/{:.4} of poses, below {CONTACT_FLOOR}, \
             so both arms are being timed on rejections",
            arm.samples,
            arm.contact_share_field,
            arm.contact_share_trimesh
        );
        assert!(
            arm.endpoint_share < ENDPOINT_BAR,
            "VOID: {:.4} of the {}^3 pose set's minimisers sit at a segment endpoint, so the 1-D \
             search has nothing interior to find",
            arm.endpoint_share,
            arm.samples
        );
    }

    let c1_holds = c1[0].speedup >= SPEEDUP_BAR && c1[0].field_us < QUERY_BUDGET_US;
    // The TriMesh half is measured and decides its own verdict; the field half is
    // VACUOUS, so the clause as registered cannot be scored HELD or FALSIFIED.
    let c3_trimesh_holds = c3.ghost_trimesh > 0;

    common::experiment::run(isomesh::experiment!("P-81"), |run| {
        for (i, arm) in c1.iter().enumerate() {
            let registered = i == 0;
            let na = || String::from("NA");
            let mut row: Row = vec![
                ("field", String::from(FbmTerrain::<f32>::NAME)),
                ("arm", format!("c1_{}", arm.samples)),
                ("chunk_cells", arm.samples.to_string()),
                ("cells_per_axis", (arm.samples - 1).to_string()),
                ("cell_size", format!("{:.6}", arm.cell_size)),
                ("triangles", arm.triangles.to_string()),
                ("queries", C1_QUERIES.to_string()),
                ("gss_iterations", GSS_ITERATIONS.to_string()),
                ("query_us_field", format!("{:.6}", arm.field_us)),
                ("query_us_trimesh", format!("{:.6}", arm.trimesh_us)),
                ("speedup", format!("{:.6}", arm.speedup)),
                ("speedup_min", format!("{:.6}", arm.speedup_min)),
                ("speedup_max", format!("{:.6}", arm.speedup_max)),
                (
                    "query_us_field_median",
                    format!("{:.6}", arm.field_us_median),
                ),
                (
                    "query_us_trimesh_median",
                    format!("{:.6}", arm.trimesh_us_median),
                ),
                (
                    "cycles_per_query_field",
                    arm.cycles_field.map_or_else(unavailable, |c| {
                        format!("{:.1}", c as f64 / C1_QUERIES as f64)
                    }),
                ),
                (
                    "cycles_per_query_trimesh",
                    arm.cycles_trimesh.map_or_else(unavailable, |c| {
                        format!("{:.1}", c as f64 / C1_QUERIES as f64)
                    }),
                ),
                (
                    "cycles_speedup",
                    match (arm.cycles_field, arm.cycles_trimesh) {
                        (Some(f), Some(t)) => format!("{:.6}", t as f64 / f.max(1) as f64),
                        _ => unavailable(),
                    },
                ),
                (
                    "ghz",
                    arm.ghz.map_or_else(unavailable, |g| format!("{g:.4}")),
                ),
                (
                    "contact_share_field",
                    format!("{:.6}", arm.contact_share_field),
                ),
                (
                    "contact_share_trimesh",
                    format!("{:.6}", arm.contact_share_trimesh),
                ),
                ("contact_disagreements", arm.disagreements.to_string()),
                ("endpoint_minimisers", format!("{:.6}", arm.endpoint_share)),
                (
                    "c1_holds",
                    (arm.speedup >= SPEEDUP_BAR && arm.field_us < QUERY_BUDGET_US).to_string(),
                ),
            ];
            if registered {
                row.extend([
                    ("c2_queries", C2_QUERIES.to_string()),
                    ("c2_block", C2_BLOCK.to_string()),
                    (
                        "contacts_bit_identical",
                        peer.map_or_else(blocked, |p| p.0.to_string()),
                    ),
                    (
                        "differing_contacts",
                        peer.map_or_else(blocked, |p| p.2.to_string()),
                    ),
                    (
                        "differing_blocks",
                        peer.map_or_else(blocked, |p| p.1.to_string()),
                    ),
                    (
                        "differing_is_bound",
                        peer.map_or_else(blocked, |p| (p.1 > 0).to_string()),
                    ),
                    (
                        "c2_blocker",
                        String::from(if peer.is_some() {
                            "none"
                        } else {
                            "no-peer-digest-table-committed"
                        }),
                    ),
                    (
                        "locally_self_consistent",
                        locally_self_consistent.to_string(),
                    ),
                    ("local_digest", format!("{:#018x}", c2.digest)),
                    ("local_repeat_digest", format!("{:#018x}", c2.repeat_digest)),
                    ("branchy_digest", format!("{:#018x}", c2.branchy_digest)),
                    ("peer_host", String::from(PEER_HOST)),
                    ("peer_target", String::from(PEER_TARGET)),
                    ("peer_rustc", String::from(PEER_RUSTC)),
                    ("c2_contacts", c2.contacts.to_string()),
                    (
                        "branch_misses_per_query_field",
                        c2.branch_misses_field
                            .map_or_else(unavailable, |v| format!("{v:.4}")),
                    ),
                    (
                        "branch_misses_per_query_branchy",
                        c2.branch_misses_branchy
                            .map_or_else(unavailable, |v| format!("{v:.4}")),
                    ),
                    ("c2_holds", peer.map_or_else(blocked, |p| p.0.to_string())),
                    ("seam_crossings", c3.seam_crossings.to_string()),
                    ("c3_steps", c3.steps.to_string()),
                    ("c3_windows", c3.windows.to_string()),
                    ("contacts_trimesh", c3.contacts_trimesh.to_string()),
                    ("contacts_field", c3.contacts_field.to_string()),
                    (
                        "edge_witness_contacts",
                        c3.edge_witness_contacts.to_string(),
                    ),
                    (
                        "internal_edge_contacts",
                        c3.internal_edge_contacts.to_string(),
                    ),
                    ("ghost_contacts_trimesh", c3.ghost_trimesh.to_string()),
                    (
                        "ghost_contacts_trimesh_rawdev",
                        c3.ghost_trimesh_rawdev.to_string(),
                    ),
                    (
                        "ghost_contacts_trimesh_15deg",
                        c3.ghost_trimesh_wide.to_string(),
                    ),
                    (
                        "ghost_contacts_trimesh_jolt",
                        c3.ghost_trimesh_jolt.to_string(),
                    ),
                    ("ghost_poses_trimesh", c3.ghost_poses.to_string()),
                    // **VACUOUS, and said out loud rather than scored HELD.** See
                    // the module docs: the registered class is defined on a
                    // structure the field arm does not have, and every instrument
                    // that tries to score it has to import the mesh as the
                    // reference - at which point it measures the mesh's own
                    // discretisation error. The three surrogates below are that
                    // error read three ways, and they disagree by 38x, which is
                    // the evidence for this cell.
                    ("ghost_contacts_field", String::from("VACUOUS")),
                    (
                        "c3_field_verdict",
                        String::from("VACUOUS-field-query-has-no-internal-edges"),
                    ),
                    (
                        "ghost_contacts_field_at_mesh_witness",
                        c3.ghost_field.to_string(),
                    ),
                    (
                        "ghost_contacts_field_projected",
                        c3.ghost_field_projected.to_string(),
                    ),
                    (
                        "ghost_contacts_field_projected_rawdev",
                        c3.ghost_field_projected_rawdev.to_string(),
                    ),
                    (
                        "field_contacts_at_internal_edge",
                        c3.field_at_internal_edge.to_string(),
                    ),
                    (
                        "field_ghost_excess_max_deg",
                        format!("{:.4}", c3.field_ghost_excess_max_deg),
                    ),
                    (
                        "field_projection_dist_max",
                        format!("{:.6}", c3.field_projection_dist_max),
                    ),
                    (
                        "field_normal_dev_max_deg",
                        format!("{:.4}", c3.field_dev_max_deg),
                    ),
                    (
                        "field_normal_excess_max_deg",
                        format!("{:.4}", c3.field_excess_max_deg),
                    ),
                    (
                        "mesh_ghost_dev_max_deg",
                        format!("{:.4}", c3.mesh_ghost_dev_max_deg),
                    ),
                    (
                        "mesh_ghost_excess_max_deg",
                        format!("{:.4}", c3.mesh_ghost_excess_max_deg),
                    ),
                    ("internal_edges", c3.internal_edges.to_string()),
                    ("boundary_edges", c3.boundary_edges.to_string()),
                    ("nonmanifold_edges", c3.nonmanifold_edges.to_string()),
                    (
                        "inactive_edge_share",
                        format!("{:.6}", c3.inactive_edge_share),
                    ),
                    ("mean_dihedral_deg", format!("{:.4}", c3.mean_dihedral_deg)),
                    ("outward_share", format!("{:.6}", c3.outward_share)),
                    ("c3_trimesh_holds", c3_trimesh_holds.to_string()),
                    ("c3_holds", String::from("VACUOUS")),
                ]);
            } else {
                row.extend([
                    ("c2_queries", na()),
                    ("c2_block", na()),
                    ("contacts_bit_identical", na()),
                    ("differing_contacts", na()),
                    ("differing_blocks", na()),
                    ("differing_is_bound", na()),
                    ("c2_blocker", na()),
                    ("locally_self_consistent", na()),
                    ("local_digest", na()),
                    ("local_repeat_digest", na()),
                    ("branchy_digest", na()),
                    ("peer_host", na()),
                    ("peer_target", na()),
                    ("peer_rustc", na()),
                    ("c2_contacts", na()),
                    ("branch_misses_per_query_field", na()),
                    ("branch_misses_per_query_branchy", na()),
                    ("c2_holds", na()),
                    ("seam_crossings", na()),
                    ("c3_steps", na()),
                    ("c3_windows", na()),
                    ("contacts_trimesh", na()),
                    ("contacts_field", na()),
                    ("edge_witness_contacts", na()),
                    ("internal_edge_contacts", na()),
                    ("ghost_contacts_trimesh", na()),
                    ("ghost_contacts_trimesh_rawdev", na()),
                    ("ghost_contacts_trimesh_15deg", na()),
                    ("ghost_contacts_trimesh_jolt", na()),
                    ("ghost_poses_trimesh", na()),
                    ("ghost_contacts_field", na()),
                    ("c3_field_verdict", na()),
                    ("ghost_contacts_field_at_mesh_witness", na()),
                    ("ghost_contacts_field_projected", na()),
                    ("ghost_contacts_field_projected_rawdev", na()),
                    ("field_contacts_at_internal_edge", na()),
                    ("field_ghost_excess_max_deg", na()),
                    ("field_projection_dist_max", na()),
                    ("field_normal_dev_max_deg", na()),
                    ("field_normal_excess_max_deg", na()),
                    ("mesh_ghost_dev_max_deg", na()),
                    ("mesh_ghost_excess_max_deg", na()),
                    ("internal_edges", na()),
                    ("boundary_edges", na()),
                    ("nonmanifold_edges", na()),
                    ("inactive_edge_share", na()),
                    ("mean_dihedral_deg", na()),
                    ("outward_share", na()),
                    ("c3_trimesh_holds", na()),
                    ("c3_holds", na()),
                ]);
            }
            run.record(&row);
        }
    });

    println!(
        "P-81 verdicts: c1 {} (speedup {:.3}x, field {:.3} us) | c2 {} | c3 trimesh-half {} \
         (trimesh {}, field VACUOUS: at-witness {} / raw-dev {} / projected {})",
        c1_holds,
        c1[0].speedup,
        c1[0].field_us,
        peer.map_or_else(blocked, |p| format!("{} ({} blocks differ)", p.0, p.1)),
        c3_trimesh_holds,
        c3.ghost_trimesh,
        c3.ghost_field,
        c3.ghost_trimesh_rawdev,
        c3.ghost_field_projected
    );
}
