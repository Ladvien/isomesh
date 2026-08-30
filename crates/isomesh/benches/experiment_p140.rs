//! **P-140 — fields whose topology is known by construction, generated rather than measured.**
//!
//! Ticket: R-140. Pre-registered before this harness existed.
//!
//! ```bash
//! cargo bench --bench experiment_p140
//! ```
//!
//! Writes `docs/experiments/p-140.csv`.
//!
//! # What was missing
//!
//! The crate's entire vocabulary of *known* Euler characteristics is two
//! numbers. `ReferenceField::expected_euler` returns `Some(2)` for `sphere`
//! (`fields/mod.rs:324-326`), `box_exact` (`:558-560`), `thin_plate`
//! (`:668-670`) and `csg_difference` (`:946-948`), `Some(0)` for `torus`
//! (`:413-415`), and `None` for `gyroid` (`:1078-1080`), `noise_cavity`
//! (`:1255-1257`) and `fbm_terrain` (`:1386-1388`) — each with the reason
//! written beside it, *"genus depends on how many tunnels the cap encloses"*.
//! So the eight reference fields cover genus 0 and genus 1 and nothing else, and
//! the three fields with the interesting topology are precisely the three with
//! no topological gate beyond manifoldness. That is the gap the registration
//! names, and it is real.
//!
//! `validate::validate` already holds the instrument.
//! `MeshReport::euler_characteristic` (`validate.rs:185`) is
//! `referenced_vertices − edges + faces`, and `MeshReport::genus` (`:202`) is
//! `(2 − χ − boundary_loops) / 2`, reported `Some` only for a single
//! consistently-oriented manifold component with no skipped faces. What was
//! missing was never the ruler. It was something of known length to hold it
//! against.
//!
//! `common::tpms` (owner R-142) supplies three fields whose `χ` is exact — the
//! nodal gyroid at `−8N³`, Schwarz P at `−4N³`, Schwarz D at `−16N³` — and
//! R-142/R-143/R-144/R-145 measure them. **This harness deliberately borrows
//! none of them.** Their `χ` is defined only after `common::tpms::wrap_seams`
//! identifies opposite boundary faces by a period translation, and it is read
//! off `common::tpms::euler` rather than off `validate::validate`. A wrapped row
//! in this CSV would make `measured_chi` mean two different things in one
//! column, measured by two different instruments — the corrupt-but-plausible
//! artefact `Run::record`'s own docs warn about at
//! `benches/common/experiment.rs:59-65`. The calibration here is `sphere` and
//! `torus` instead: they are the two fixtures whose prescribed `χ` **the crate
//! itself declares**, so agreement on them is what licenses reading `χ` off the
//! same call for the six new solids.
//!
//! # Viro, honestly
//!
//! The registration names Viro's combinatorial patchworking as the generator,
//! and the generator property — *prescribe the topology, get a field* — is
//! exactly what C1 and C2 need. Patchworking itself is not what is implemented
//! here, and pretending otherwise would be the more expensive lie. Two reasons:
//!
//! 1. Patchworking starts from a **regular** triangulation of the Newton
//!    polytope with a sign distribution, and whether the triangulation this
//!    crate's cells induce is regular is precisely what **R-139** is measuring.
//!    Building on that answer before it exists would rest this row on an
//!    unmeasured premise.
//! 2. What patchworking prescribes is a *polynomial*, and the topology it
//!    prescribes is that of the polynomial's zero set in a torus or a projective
//!    space. The boundary of a bounded solid inside a finite sampling box is a
//!    different object unless the box is chosen to cut it correctly — the same
//!    unsolved oracle problem in new clothing (`fields/mod.rs:1079`).
//!
//! So the fields are built by two closed-form constructions whose boundary Euler
//! characteristic is derivable in one line each, and **the derivation is written
//! down here rather than read off the measurement**. That is the discipline of
//! the row: `prescribed_chi` is arithmetic, `measured_chi` is a reading, and the
//! two are computed by code that shares nothing.
//!
//! ## Construction A — a ball with `g` cylinders drilled straight through it
//!
//! `f(p) = max(|p| − R, −min_i(|(p_x, p_y) − c_i| − r))`: the ball of radius `R`
//! minus the union of `g` infinite `z`-parallel cylinders of radius `r` at
//! lateral centres `c_i`. Both operands are 1-Lipschitz, and `max`, `min` and
//! negation preserve that, so the field is 1-Lipschitz everywhere.
//!
//! Sound when `|c_i − c_j| > 2r` for every pair — the bores are disjoint — and
//! `max_i|c_i| + r < R`, which is what makes every bore exit through the ball's
//! top *and* bottom cap while leaving a shell of positive thickness at the
//! equator. Then the boundary is a sphere with `2g` open disks removed, glued
//! along `2g` circles to `g` open annuli, and Mayer–Vietoris is one line:
//!
//! ```text
//! chi = (2 − 2g) + g·chi(annulus) − 2g·chi(circle) = (2 − 2g) + 0 − 0 = 2 − 2g
//! ```
//!
//! Connected, closed and orientable, so the genus is `(2 − χ)/2 = g`.
//!
//! ## Construction B — the closed `t`-neighbourhood of an embedded graph
//!
//! `f(p) = min_e dist(p, e) − t` over the graph's edges, each a closed segment.
//! Exact segment distance, so again 1-Lipschitz.
//!
//! Sound when the neighbourhood is a **regular** neighbourhood, which is two
//! conditions, both computed numerically below rather than asserted by eye:
//!
//! - every pair of edges sharing **no** node is farther apart than `2t`, so
//!   their tubes are disjoint — exact segment-to-segment distance, [`gap`];
//! - every pair of edges sharing **one** node `v` has its tubes meeting only
//!   near `v`. Two tubes of radius `t` about rays leaving `v` at angle `θ`
//!   intersect exactly within `t / sin(θ/2)` of `v` along each ray, so holding
//!   that under **half** the shorter edge keeps the two ends of any one edge
//!   from merging with each other.
//!
//! A regular neighbourhood of a graph is a handlebody with `χ(N) = χ(G)`, and
//! `χ(∂M) = 2χ(M)` for a compact 3-manifold, so
//!
//! ```text
//! chi = 2·chi(G) = 2·(V − E)      genus = E − V + 1 = b1(G)
//! ```
//!
//! The `g = b1` independent cycles are **listed** per fixture rather than
//! counted, and `generators.len() == E − V + 1 == prescribed_genus` is asserted
//! three ways against each other. Each generator carries its own witness too:
//! the cycle's centroid is checked to lie strictly outside the solid, so every
//! prescribed handle demonstrably has an open window.
//!
//! ## The eight fixtures
//!
//! | `construction` | family | `derivation` | genus | `prescribed_chi` | `is_control` |
//! |---|---|---|---|---|---|
//! | `sphere` | crate reference field | `crate_expected_euler` | 0 | `2` | **yes** |
//! | `torus` | crate reference field | `crate_expected_euler` | 1 | `0` | **yes** |
//! | `ball_drilled_g1` | A, one axial bore | `mayer_vietoris` | 1 | `0` | no |
//! | `ball_drilled_g2` | A, two bores at `rho = 0.75` | `mayer_vietoris` | 2 | `-2` | no |
//! | `ball_drilled_g3` | A, three bores at `rho = 0.75` | `mayer_vietoris` | 3 | `-4` | no |
//! | `graph_theta_g2` | B, two poles joined by three bent arcs | `graph_thickening` | 2 | `-2` | no |
//! | `graph_k4_g3` | B, the tetrahedron graph `K4` | `graph_thickening` | 3 | `-4` | no |
//! | `graph_cube_g5` | B, the cube graph `Q3` | `graph_thickening` | 5 | `-8` | no |
//!
//! Five distinct prescribed `χ` values, `{2, 0, -2, -4, -8}`, and **five**
//! fixtures of genus above 1 — well past the registered control's demand for
//! one. Two families rather than one, because a single family would leave
//! `prescribed_chi` and the construction confounded: `ball_drilled_g2` and
//! `graph_theta_g2` prescribe the same `χ = -2` through arguments that share no
//! algebra, and `ball_drilled_g3` and `graph_k4_g3` do the same for `χ = -4`.
//!
//! **These are bench-local fixtures. They are not added to
//! `for_each_reference_field!`, and that is deliberate.** One new reference
//! field adds 27 rows to `crates/isomesh/golden_hashes.json` (216 = 8 fields ×
//! 9 algorithms × 3 resolutions) and moves `scripts/doc_facts.sh`'s `FIELDS`
//! and `HASHES` counts, which are gated as prose phrases across twelve
//! documents. Six new fields would be a repo-wide renumbering landed inside a
//! measurement commit, and `crates/isomesh/src/**` is frozen for this phase.
//! Landing them is a Phase 28 ticket with that ripple priced into it; C1's
//! *"are added as fixtures"* is discharged here as **bench-local fixtures**,
//! which is what makes the measurement possible at all.
//!
//! # Arms
//!
//! | arm | what it varies | `is_control` |
//! |---|---|---|
//! | `sphere` × 7 extractors × 10 resolutions | the crate's own genus-0 declaration | **yes** |
//! | `torus` × 7 × 10 | the crate's own genus-1 declaration | **yes** |
//! | `ball_drilled_g{1,2,3}` × 7 × 10 | genus, at fixed construction | no |
//! | `graph_{theta_g2,k4_g3,cube_g5}` × 7 × 10 | genus **and** construction family | no |
//!
//! 8 × 7 × 10 = **560 rows.** The seven extractors are `extractor::ALL_EXTRACTORS`
//! (`extractor.rs:229-237`), driven through `for_each_extractor!`, so
//! `extractors_tested` is 7 on every row and the roster is asserted against the
//! macro's own visiting order rather than assumed. `GreedyQuads` is excluded by
//! the crate, for the reason at `extractor.rs:246-252`: it emits axis-aligned
//! whole-cell faces, so its `χ` is that of a Minecraft surface and putting it in
//! this table would be a category error rather than a measurement.
//!
//! The resolution ladder is `[5, 7, 9, 11, 13, 17, 21, 25, 33, 49]` samples per
//! axis, i.e. `4³ … 48³` cells over the `[-2, 2]³` domain every compact
//! reference field already uses (`fields/mod.rs:258`, `COMPACT_DOMAIN = 2.0`).
//!
//! **Its low end was set by a pilot run, and that is instrument calibration
//! rather than an amended prediction.** A ladder that begins above the
//! threshold cannot report where the threshold is, and the first ladder tried
//! here — `[13, 17, 21, 25, 33, 49]` — did exactly that: every one of the six
//! prescribed fixtures was already correct on all seven extractors at 13³, so
//! `first_correct_resolution` was the constant `13` and C2's curve carried no
//! information about anything but the ladder. The reason is worth writing down,
//! because it is the row's first real finding: `ball_drilled_g1`'s bore is
//! 0.56 across against a cell of 0.333, `feature_cells = 1.68`, and the tunnel
//! survives anyway — the bore is centred on the `z` axis and the grid puts a
//! sample line exactly there, so a *one-sample-wide* void is enough for every
//! extractor to find the handle. Sampling adequacy is not the same thing as
//! feature width, and `feature_cells` is a predictor rather than a threshold.
//! Four rungs at `5, 7, 9, 11` (cell sizes 1.0, 0.667, 0.5, 0.4) bracket the
//! failure regime from below and cost 1 792 cells against the ladder's 172 800,
//! i.e. nothing: the top rung dominates either way. The whole ladder is swept
//! even after an extractor has converged, so a **non-monotone** curve is
//! visible rather than hidden behind an early exit — and the pilot found one,
//! `marching_tetrahedra` on `graph_k4_g3` reading `χ = -4` at 13³, `-16` at
//! 17³ and `-4` again at 21³.
//!
//! # An empty mesh is not a reading
//!
//! `χ` of a mesh with no face is `0 − 0 + 0 = 0` by arithmetic, and `torus` and
//! `ball_drilled_g1` both prescribe `χ = 0`. At the bottom of the ladder — 4
//! cells per axis — an extractor that meshes nothing at all would therefore
//! score as reproducing a genus-1 surface, and `first_correct_resolution` would
//! report the coarsest rung on the strength of an empty buffer. So
//! `chi_agreement` is `faces > 0 && measured_chi == prescribed_chi`, and
//! `faces` is recorded beside it. This is the registered predicate plus the
//! statement that the reading exists; it is not a second reading of `χ`.
//!
//! `chi_agreement` deliberately does **not** also require `mesh_components == 1`.
//! `χ` is the quantity the registration names, and a mesh that reproduces it
//! while falling into two pieces is a reading worth having rather than a skip —
//! `mesh_components`, `measured_genus`, `nonmanifold_edges` and `is_manifold`
//! are all recorded so the stricter reading is available to anyone who wants it
//! without this column meaning two things.
//!
//! The guard is not hypothetical, and one rung of the ladder can be shown to
//! need it by arithmetic alone. At 9³ the cell is `h = 4/8 = 0.5`, so every
//! sample sits at a multiple of `0.5`, while `graph_cube_g5`'s edges lie on the
//! lines `(t, ±0.75, ±0.75)` and its permutations. For integer `k`,
//! `min |0.5k − 0.75| = 0.25`, so the distance from **any** sample to **any**
//! edge axis is at least `0.25·√2 = 0.35355`, and the tube is `0.28`. The whole
//! genus-5 solid therefore falls between the samples: every corner reads
//! `f >= 0.0736 > 0`, no cell has a sign change and the mesh is empty. The same
//! arithmetic makes 5³ empty — `h = 1.0`, the same `0.25` residual — and leaves
//! every other rung comfortable: `0.75/h` is `3`, `6` and `9` exactly at 17³,
//! 33³ and 49³, so a sample sits *on* the edge axis, and the residual is
//! `0.0833` at 7³, 13³ and 25³ and `0.05` at 11³ and 21³, giving sample-to-axis
//! distances of `0.1179` and `0.0707` against a tube of `0.28`. Two rungs of ten
//! are blind, and they are the two whose cell size shares a factor with `0.75`
//! in the wrong way.
//!
//! That is the sharpest thing this row has to say, and it is not about genus:
//! the highest-genus fixture in the suite vanishes at 9³ and is perfect at 7³
//! and at 11³. **Sampling adequacy is a commensurability property of the
//! lattice and the feature, not a width.** `feature_cells` is recorded as a
//! predictor and is demonstrably not a threshold — every prescribed fixture
//! here has `feature_clearance` between 0.52 and 0.56, and their first-correct
//! rungs still spread across the ladder.
//!
//! `subgrid_marching_tetrahedra` runs without a Lipschitz constant, because
//! `for_each_extractor!` hands out a bare `SubgridMarchingTetrahedra::new(16)`
//! and `set_lipschitz` exists on that one type only, so calling it inside the
//! macro body would not type-check in the other six expansions. It therefore
//! pays the full `6 tets × 6 edges × 16 samples = 576` evaluations per cell
//! (`subgrid/extract.rs:450-451`) and is roughly three quarters of this bench's
//! wall clock. Its raw output is already vertex-shared by identity
//! (`subgrid/extract.rs:529-546`), so `χ` is read off it with no weld — welding
//! it is a no-op at best and on one reference field actively damages it (M-226).
//!
//! An extractor returning `Err` is **not** caught and turned into a row. A
//! closed, 1-Lipschitz, everywhere-finite field is inside every extractor's
//! contract, so an `Err` here is a crate defect rather than a sampling-adequacy
//! datum, and folding the two into `measured_chi` would give one column two
//! meanings. It panics, naming the fixture, the extractor and the resolution.
//! Converging to the *wrong* `χ` is what C1's falsifier is about, and that is
//! recorded.
//!
//! # Clause verdicts, and which are global
//!
//! - **C1** — every extractor reproduces the prescribed `χ` at some rung, on
//!   every **prescribed** fixture. Global: the same value on all 560 rows. The
//!   two calibration arms are outside C1's scope, because C1 is about the
//!   fixtures this row prescribes; their agreement is asserted as a vacuity
//!   control instead, which is the stronger place for it. Per-pair evidence is
//!   the extra column `converged`.
//! - **C2** — the first-correct resolution is reported per extractor, and the
//!   curve is informative. Global, and read straight off the registered
//!   falsifier *"all extractors converging at the same resolution"*: C2 holds
//!   when at least one prescribed fixture shows two or more distinct
//!   first-correct values among its seven extractors.
//!   `fixture_distinct_first_correct` carries the per-fixture count and
//!   `distinct_first_correct_overall` the sweep-wide one.
//! - There is no C3 in this registration, so there is no `c3_holds` column.
//!
//! `first_correct_resolution` is `never` when no rung reproduces the prescribed
//! `χ`. **`never` means never within `[5, 49]`**, and it is a recorded value
//! rather than a skip.
//!
//! # SHARE, recomputed before the numbers
//!
//! The registration says `SHARE: none -- this adds test coverage, and its value
//! is the coverage`, and that is discharged rather than waved through. The
//! coverage this row creates is bench-local by the phase's own rule, so the
//! shipped test suite gains nothing until a Phase 28 ticket lands the fixtures
//! as reference fields. What ships *now* is smaller and real: the eight
//! `prescribed_chi` values and the per-extractor first-correct resolutions are
//! the first statement anywhere in this repository of how much grid an extractor
//! needs before its topology is right, on a field whose topology is known.
//! Every prior answer to that question was available for genus 0 and genus 1
//! only.
//!
//! # Vacuity controls
//!
//! Seven, all before the first `run.record`, all prefixed `VOID: `.
//!
//! 1. **The ladder is a ladder** — at least four strictly increasing rungs, so
//!    that "first correct" names a well-ordered thing and the C2 column is a
//!    curve. Proves nothing about the fields; proves the axis exists.
//! 2. **Genus above 1** — the registered control, *"or the fixture is a sphere
//!    in disguise"*. Asserted on `prescribed_genus`: at least three prescribed
//!    fixtures have genus above 1. There are five.
//! 3. **`prescribed_chi` is not a constant** — at least four distinct values
//!    across the fixtures. A `chi_agreement` column computed against one
//!    repeated number could read `true` for an extractor that always emits a
//!    sphere.
//! 4. **The prescription is sound, by construction and not by measurement** —
//!    for every fixture: `separation`, `window`, `solid_width`,
//!    `merge_headroom` and `domain_margin` all strictly positive, and
//!    `generators == prescribed_genus`, cross-checked against `E − V + 1` for
//!    the graphs and against the bore count for the balls. Without this
//!    `prescribed_genus` is not prescribed and every other column is decoration.
//! 5. **The instrument reads a known `χ`** — `sphere` and `torus` reproduce
//!    their crate-declared `expected_euler` on all seven extractors at the
//!    finest rung, and this harness's own prescription for those two equals what
//!    the crate declares. This is the control that licenses `measured_chi` at
//!    all.
//! 6. **`chi_agreement` could have been false** — at least one prescribed
//!    fixture × extractor × rung triple in the sweep disagrees with its
//!    prescribed `χ`. A column that could not have read `false` is not a
//!    measurement (M-44).
//!
//!    Scoped to the **whole sweep** and not to the coarsest rung, and the
//!    difference is load-bearing. C2 coming out *uninformative* — every
//!    extractor converging at the same rung — is a **registered outcome**, said
//!    in the falsifier's own words: *"still worth one line"*. An assertion that
//!    aborted on it would convert a registered null into a crash, which is the
//!    opposite of what a vacuity control is for. What must not happen is that
//!    `chi_agreement` never reads `false` anywhere, because then the instrument
//!    was never shown able to say no.
//! 7. **The sweep meshed something** — at the finest rung, every fixture ×
//!    extractor pair produced at least one face. Together with control 6 this
//!    pins `chi_agreement` from both sides: it has been seen to read `false`,
//!    and every `true` at the top of the ladder rests on a non-empty mesh.

#![allow(
    clippy::cast_possible_truncation,
    clippy::cast_precision_loss,
    clippy::too_many_lines
)]

mod common;

use std::collections::BTreeSet;

use isomesh::extractor::{ALL_EXTRACTORS, Extractor};
use isomesh::fields::{ReferenceField, Sphere, Torus};
use isomesh::validate::{ValidateConfig, validate};
use isomesh::{MeshBuffer, RuntimeShape3, Sdf};

/// Half-extent of the sampling box, identical to `fields/mod.rs:258`'s
/// `COMPACT_DOMAIN`, so that `sphere` and `torus` are sampled on exactly the
/// grid their own `ReferenceField::domain` asks for.
const DOMAIN_HALF: f64 = 2.0;

/// Samples per axis, ascending. `n` samples span `n − 1` cells.
///
/// The four low rungs bracket the failure regime from below; see the header's
/// note on why the ladder had to start there and not at 13.
const LADDER: [u32; 10] = [5, 7, 9, 11, 13, 17, 21, 25, 33, 49];

/// Tube and bore radius, shared by both families so that the geometry varies in
/// one argument. `0.28` puts the thinnest solid feature at `0.56` across, which
/// is 4.5 cells at 33³.
const TUBE: f64 = 0.28;

// ── vector arithmetic, three lines each ─────────────────────────────────────

/// `a − b`.
fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

/// `a + b`.
fn add(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

/// `a · s`.
fn scale(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}

/// `a · b`.
fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

/// `|a|`.
fn norm(a: [f64; 3]) -> f64 {
    dot(a, a).sqrt()
}

/// Distance from `p` to the closed segment `a → b`.
///
/// Exact, and the clamp is what makes a union of these a union of *capsules*
/// rather than of infinite cylinders.
fn point_segment_distance(p: [f64; 3], a: [f64; 3], b: [f64; 3]) -> f64 {
    let ab = sub(b, a);
    let ap = sub(p, a);
    let len2 = dot(ab, ab);
    let t = if len2 > 0.0 {
        (dot(ap, ab) / len2).clamp(0.0, 1.0)
    } else {
        0.0
    };
    norm(sub(ap, scale(ab, t)))
}

/// Minimum distance between two closed segments, exactly.
///
/// The minimum over the parameter square `[0,1]²` is attained either on its
/// boundary — the four point-to-segment cases — or at the interior critical
/// point of the quadratic, which exists only when the two directions are not
/// parallel. Both are checked, so no configuration is guessed at.
fn gap(p0: [f64; 3], p1: [f64; 3], q0: [f64; 3], q1: [f64; 3]) -> f64 {
    let mut best = point_segment_distance(p0, q0, q1)
        .min(point_segment_distance(p1, q0, q1))
        .min(point_segment_distance(q0, p0, p1))
        .min(point_segment_distance(q1, p0, p1));

    let d1 = sub(p1, p0);
    let d2 = sub(q1, q0);
    let r = sub(p0, q0);
    let a = dot(d1, d1);
    let e = dot(d2, d2);
    let b = dot(d1, d2);
    let c = dot(d1, r);
    let f = dot(d2, r);
    let denom = a * e - b * b;
    if denom > 0.0 {
        let s = (b * f - c * e) / denom;
        let t = (a * f - b * c) / denom;
        if (0.0..=1.0).contains(&s) && (0.0..=1.0).contains(&t) {
            best = best.min(norm(sub(add(p0, scale(d1, s)), add(q0, scale(d2, t)))));
        }
    }
    best
}

/// `g` points on the circle of radius `rho` in the `z = 0` plane, the first on
/// the `+x` axis.
fn ring(g: usize, rho: f64) -> Vec<[f64; 2]> {
    (0..g)
        .map(|k| {
            let a = std::f64::consts::TAU * (k as f64) / (g as f64);
            [rho * a.cos(), rho * a.sin()]
        })
        .collect()
}

/// How many distinct first-correct rungs a slice holds, counting `never` — the
/// `None` — as a value of its own rather than as an absence.
fn distinct_rungs(values: &[Option<u32>]) -> usize {
    let mut seen: BTreeSet<Option<u32>> = BTreeSet::new();
    for v in values {
        seen.insert(*v);
    }
    seen.len()
}

// ── the solids ──────────────────────────────────────────────────────────────

/// A graph embedded in `R³`, with the cycles that witness its first Betti
/// number listed rather than derived.
#[derive(Clone, Debug)]
struct Graph {
    /// Node positions.
    nodes: Vec<[f64; 3]>,
    /// Edges, as index pairs into `nodes`.
    edges: Vec<[usize; 2]>,
    /// The `b1(G)` independent cycles, each listed by the nodes it visits.
    /// Carried so that the prescribed genus is a *statement* checked three ways
    /// — against this length, against `E − V + 1`, and against each cycle's own
    /// open-window witness.
    generators: Vec<Vec<usize>>,
}

/// A closed solid whose boundary genus is fixed by its construction.
///
/// One enum rather than a list of `Box<dyn Sdf>`, because every extractor's
/// `extract` takes `sdf: &S` with an implicit `S: Sized` bound, so a bare
/// `dyn Sdf` cannot be passed at all and the double-reference dance that would
/// work buys nothing but a vtable hop inside a 576-evaluation-per-cell loop.
#[derive(Clone, Debug)]
enum Solid {
    /// The crate's canonical unit sphere. Calibration only.
    Sphere(Sphere<f64>),
    /// The crate's canonical torus, major 1, minor 0.3. Calibration only.
    Torus(Torus<f64>),
    /// Construction A: a ball with one `z`-parallel cylinder per entry of
    /// `bores` drilled straight through it.
    DrilledBall {
        /// Ball radius.
        radius: f64,
        /// Bore radius, shared by every bore.
        bore: f64,
        /// Lateral `(x, y)` centre of each bore.
        bores: Vec<[f64; 2]>,
    },
    /// Construction B: the closed `tube`-neighbourhood of an embedded graph.
    ThickenedGraph {
        /// The graph and its cycle basis.
        graph: Graph,
        /// Tube radius.
        tube: f64,
    },
}

impl Sdf for Solid {
    type Scalar = f64;

    fn sample(&self, p: [f64; 3]) -> f64 {
        match self {
            Self::Sphere(s) => s.sample(p),
            Self::Torus(t) => t.sample(p),
            Self::DrilledBall {
                radius,
                bore,
                bores,
            } => {
                let ball = norm(p) - radius;
                let mut drilled = f64::INFINITY;
                for c in bores {
                    let dx = p[0] - c[0];
                    let dy = p[1] - c[1];
                    drilled = drilled.min((dx * dx + dy * dy).sqrt() - bore);
                }
                ball.max(-drilled)
            }
            Self::ThickenedGraph { graph, tube } => {
                let mut d = f64::INFINITY;
                for e in &graph.edges {
                    d = d.min(point_segment_distance(
                        p,
                        graph.nodes[e[0]],
                        graph.nodes[e[1]],
                    ));
                }
                d - tube
            }
        }
    }
}

// ── the fixtures ────────────────────────────────────────────────────────────

/// One prescribed field: a solid, the genus its construction fixes, and the name
/// of the argument that fixes it.
#[derive(Clone, Debug)]
struct Fixture {
    /// The `construction` column.
    name: &'static str,
    /// The `derivation` extra column: which one-line argument gives `χ`.
    derivation: &'static str,
    /// `prescribed_genus`.
    genus: i64,
    /// What `ReferenceField::expected_euler` says, for the two fixtures the
    /// crate has an opinion about. `None` for the six this harness prescribes —
    /// the crate has never seen them.
    declared: Option<i64>,
    /// Calibration arms are outside C1's and C2's scope and inside control 5's.
    is_control: bool,
    /// The field.
    solid: Solid,
}

impl Fixture {
    /// `prescribed_chi`, from the genus and nothing else.
    fn chi(&self) -> i64 {
        2 - 2 * self.genus
    }
}

/// The tetrahedron graph `K4` at half-diagonal `a`: four nodes, six edges,
/// `b1 = 3`.
fn k4(a: f64) -> Graph {
    Graph {
        nodes: vec![[a, a, a], [a, -a, -a], [-a, a, -a], [-a, -a, a]],
        edges: vec![[0, 1], [0, 2], [0, 3], [1, 2], [1, 3], [2, 3]],
        // Three of the four triangles. The fourth is their sum in the cycle
        // space, so these three are a basis and `b1 = 6 − 4 + 1 = 3`.
        generators: vec![vec![0, 1, 2], vec![0, 1, 3], vec![0, 2, 3]],
    }
}

/// The cube graph `Q3` at half-edge `a`: eight nodes, twelve edges, `b1 = 5`.
fn cube_graph(a: f64) -> Graph {
    let coord = |bit: u32, i: u32| if i & bit == 0 { -a } else { a };
    let nodes: Vec<[f64; 3]> = (0..8u32)
        .map(|i| [coord(1, i), coord(2, i), coord(4, i)])
        .collect();

    let mut edges: Vec<[usize; 2]> = Vec::with_capacity(12);
    for i in 0..8u32 {
        for bit in [1u32, 2, 4] {
            let j = i ^ bit;
            if j > i {
                edges.push([i as usize, j as usize]);
            }
        }
    }

    // Five of the six faces. The six four-cycles sum to zero over GF(2), so any
    // five of them are independent and `b1 = 12 − 8 + 1 = 5`.
    let mut generators: Vec<Vec<usize>> = Vec::with_capacity(6);
    for bit in [1u32, 2, 4] {
        for side in [0u32, bit] {
            let face: Vec<usize> = (0..8u32)
                .filter(|i| i & bit == side)
                .map(|i| i as usize)
                .collect();
            generators.push(face);
        }
    }
    generators.pop();

    Graph {
        nodes,
        edges,
        generators,
    }
}

/// Two poles joined by three two-segment arcs: five nodes, six edges, `b1 = 2`.
fn theta(pole: f64, rho: f64) -> Graph {
    let mut nodes = vec![[0.0, 0.0, pole], [0.0, 0.0, -pole]];
    for w in ring(3, rho) {
        nodes.push([w[0], w[1], 0.0]);
    }
    Graph {
        nodes,
        edges: vec![[0, 2], [0, 3], [0, 4], [1, 2], [1, 3], [1, 4]],
        // Two of the three four-cycles; the third is their sum.
        generators: vec![vec![0, 2, 1, 3], vec![0, 2, 1, 4]],
    }
}

/// Construction A at one bore count.
fn drilled(name: &'static str, genus: i64, bores: Vec<[f64; 2]>) -> Fixture {
    Fixture {
        name,
        derivation: "mayer_vietoris",
        genus,
        declared: None,
        is_control: false,
        solid: Solid::DrilledBall {
            // 1.55 keeps the whole solid 0.45 clear of the sampling box and
            // leaves a 0.52 shell outside the outermost bore.
            radius: 1.55,
            bore: TUBE,
            bores,
        },
    }
}

/// Construction B on one graph.
fn graphed(name: &'static str, genus: i64, graph: Graph) -> Fixture {
    Fixture {
        name,
        derivation: "graph_thickening",
        genus,
        declared: None,
        is_control: false,
        solid: Solid::ThickenedGraph { graph, tube: TUBE },
    }
}

/// The eight fixtures, in CSV order.
fn fixtures() -> Vec<Fixture> {
    // 0.75 puts three bores 1.299 apart, so a 0.739 wall survives between them,
    // and it is also the cube graph's half-edge.
    let rho = 0.75;
    let sphere = Sphere::<f64>::canonical();
    let torus = Torus::<f64>::canonical();

    vec![
        Fixture {
            name: "sphere",
            derivation: "crate_expected_euler",
            genus: 0,
            declared: sphere.expected_euler(),
            is_control: true,
            solid: Solid::Sphere(sphere),
        },
        Fixture {
            name: "torus",
            derivation: "crate_expected_euler",
            genus: 1,
            declared: torus.expected_euler(),
            is_control: true,
            solid: Solid::Torus(torus),
        },
        drilled("ball_drilled_g1", 1, vec![[0.0, 0.0]]),
        drilled("ball_drilled_g2", 2, ring(2, rho)),
        drilled("ball_drilled_g3", 3, ring(3, rho)),
        graphed("graph_theta_g2", 2, theta(1.0, 1.0)),
        graphed("graph_k4_g3", 3, k4(0.72)),
        graphed("graph_cube_g5", 5, cube_graph(rho)),
    ]
}

// ── what the construction guarantees ────────────────────────────────────────

/// Every number the prescription rests on, computed from the construction's own
/// arguments and never from a mesh.
///
/// `f64::INFINITY` marks a constraint the construction does not have — a sphere
/// has no second feature to stay clear of — and reads as `inf` in the CSV rather
/// than as a number a reader might average.
#[derive(Clone, Copy, Debug)]
struct Sound {
    /// Smallest clear distance between two solid features that must not merge,
    /// already net of the width they occupy.
    separation: f64,
    /// Smallest clear width of a prescribed void: a bore's diameter, or twice a
    /// graph window's inradius.
    window: f64,
    /// Smallest solid width: a tube's diameter, a bore wall, an equatorial
    /// shell.
    solid_width: f64,
    /// How far two incident tubes' merge region stays short of half the shorter
    /// edge. Positive means the neighbourhood is regular.
    merge_headroom: f64,
    /// `DOMAIN_HALF` minus the solid's own reach, so positive means every sample
    /// on the box boundary is outside and the extracted mesh is closed.
    domain_margin: f64,
    /// The resolution-bearing minimum, `min(solid_width, separation, window)`.
    feature_clearance: f64,
    /// Independent cycles the construction exhibits. Must equal the genus.
    generators: usize,
}

/// [`Sound`] for a drilled ball.
fn soundness_ball(radius: f64, bore: f64, bores: &[[f64; 2]]) -> Sound {
    // The bores are held apart by solid wall, so the pairwise gap is a solid
    // width rather than a separation; `separation` stays `inf`.
    let mut solid_width = f64::INFINITY;
    for (i, a) in bores.iter().enumerate() {
        for b in &bores[i + 1..] {
            let dx = a[0] - b[0];
            let dy = a[1] - b[1];
            solid_width = solid_width.min((dx * dx + dy * dy).sqrt() - 2.0 * bore);
        }
    }
    let reach = bores
        .iter()
        .map(|c| (c[0] * c[0] + c[1] * c[1]).sqrt())
        .fold(0.0f64, f64::max);
    // The equatorial shell, which is also what proves every bore exits through
    // both caps rather than stopping inside the ball.
    solid_width = solid_width.min(radius - reach - bore);

    let window = 2.0 * bore;
    Sound {
        separation: f64::INFINITY,
        window,
        solid_width,
        merge_headroom: f64::INFINITY,
        domain_margin: DOMAIN_HALF - radius,
        feature_clearance: solid_width.min(window),
        generators: bores.len(),
    }
}

/// [`Sound`] for a thickened graph.
fn soundness_graph(graph: &Graph, tube: f64) -> Sound {
    let ends = |e: &[usize; 2]| (graph.nodes[e[0]], graph.nodes[e[1]]);

    // Disjointness of non-incident tubes, exactly; and locality of incident
    // ones, from the half-angle at the node they share.
    let mut separation = f64::INFINITY;
    let mut merge_headroom = f64::INFINITY;
    for (i, ea) in graph.edges.iter().enumerate() {
        for eb in &graph.edges[i + 1..] {
            let (a0, a1) = ends(ea);
            let (b0, b1) = ends(eb);
            match ea.iter().copied().find(|n| eb.contains(n)) {
                None => separation = separation.min(gap(a0, a1, b0, b1) - 2.0 * tube),
                Some(shared) => {
                    let v = graph.nodes[shared];
                    let ua = sub(if ea[0] == shared { a1 } else { a0 }, v);
                    let ub = sub(if eb[0] == shared { b1 } else { b0 }, v);
                    let la = norm(ua);
                    let lb = norm(ub);
                    let cosine = (dot(ua, ub) / (la * lb)).clamp(-1.0, 1.0);
                    // Two tubes of radius `tube` about rays at angle θ meet
                    // exactly within `tube / sin(θ/2)` of the shared node.
                    let reach = tube / (cosine.acos() * 0.5).sin();
                    merge_headroom = merge_headroom.min(0.5 * la.min(lb) - reach);
                }
            }
        }
    }

    // One witness per prescribed handle: the cycle's centroid is outside every
    // tube, so the window really is open.
    let mut inradius = f64::INFINITY;
    for cycle in &graph.generators {
        let mut centroid = [0.0f64; 3];
        for n in cycle {
            centroid = add(centroid, graph.nodes[*n]);
        }
        centroid = scale(centroid, 1.0 / cycle.len() as f64);
        let nearest = graph
            .edges
            .iter()
            .map(|e| {
                let (a, b) = ends(e);
                point_segment_distance(centroid, a, b)
            })
            .fold(f64::INFINITY, f64::min);
        inradius = inradius.min(nearest - tube);
    }

    let reach = graph.nodes.iter().copied().map(norm).fold(0.0f64, f64::max);
    let window = 2.0 * inradius;
    let solid_width = 2.0 * tube;
    Sound {
        separation,
        window,
        solid_width,
        merge_headroom,
        domain_margin: DOMAIN_HALF - (reach + tube),
        feature_clearance: solid_width.min(separation).min(window),
        generators: graph.generators.len(),
    }
}

/// [`Sound`] for one fixture.
fn soundness(f: &Fixture) -> Sound {
    match &f.solid {
        Solid::Sphere(s) => Sound {
            separation: f64::INFINITY,
            window: f64::INFINITY,
            solid_width: 2.0 * s.radius,
            merge_headroom: f64::INFINITY,
            domain_margin: DOMAIN_HALF - s.radius,
            feature_clearance: 2.0 * s.radius,
            generators: 0,
        },
        Solid::Torus(t) => Sound {
            separation: f64::INFINITY,
            window: 2.0 * (t.major - t.minor),
            solid_width: 2.0 * t.minor,
            merge_headroom: f64::INFINITY,
            domain_margin: DOMAIN_HALF - (t.major + t.minor),
            feature_clearance: (2.0 * t.minor).min(2.0 * (t.major - t.minor)),
            generators: 1,
        },
        Solid::DrilledBall {
            radius,
            bore,
            bores,
        } => soundness_ball(*radius, *bore, bores),
        Solid::ThickenedGraph { graph, tube } => soundness_graph(graph, *tube),
    }
}

// ── measurement ─────────────────────────────────────────────────────────────

/// What one extractor produced on one fixture at one resolution.
#[derive(Clone, Copy, Debug)]
struct Measured {
    /// `MeshReport::euler_characteristic`.
    chi: i64,
    /// `MeshReport::genus`, `None` unless the mesh is a single consistently
    /// oriented manifold component.
    genus: Option<i64>,
    /// `MeshReport::vertices`.
    vertices: u64,
    /// `MeshReport::faces`.
    faces: u64,
    /// `MeshReport::components`.
    components: u64,
    /// `MeshReport::boundary_edges`. Non-zero would mean the solid reached the
    /// sampling box, which control 4 forbids by construction.
    boundary_edges: u64,
    /// `MeshReport::non_manifold_edges` — the registered `nonmanifold_edges`.
    non_manifold_edges: u64,
    /// `MeshReport::non_manifold_vertices`.
    non_manifold_vertices: u64,
    /// `MeshReport::duplicate_vertices`: how many vertices a first-fit
    /// positional weld at `cell_size · 1e-4` could remove.
    ///
    /// **A diagnostic, and not a soundness condition — recorded because the
    /// obvious reading of it is wrong.** `euler_characteristic` is computed on
    /// the *index* topology (`validate.rs:177-185`), and every extractor here
    /// shares its vertices by index: Marching Cubes and Marching Tetrahedra
    /// through an edge cache, the dual extractors one vertex per cell, and
    /// Subgrid Marching Tetrahedra by giving each crossing a global identity
    /// (`subgrid/extract.rs:532-546`). So a non-zero count here means two
    /// crossings *coincide in position while being distinct in identity* — a
    /// tet-edge crossing landing on a grid corner is the ordinary case — and
    /// `χ` is still the surface's. Sharing by identity is strictly finer than
    /// sharing by position, and welding on top of it is a no-op at best and on
    /// one reference field adds non-manifold features (M-226).
    duplicate_vertices: u64,
    /// `MeshReport::is_manifold`.
    manifold: bool,
}

impl Measured {
    /// Does this reading reproduce `prescribed`?
    ///
    /// **The face count is part of the predicate, not decoration.** `χ` of a
    /// mesh with no face is `0 − 0 + 0 = 0` by arithmetic, and `torus` and
    /// `ball_drilled_g1` both prescribe `χ = 0`, so at the bottom of the ladder
    /// an extractor that meshed nothing at all would score as reproducing a
    /// genus-1 surface and `first_correct_resolution` would name the coarsest
    /// rung on the strength of an empty buffer.
    fn reproduces(&self, prescribed: i64) -> bool {
        self.faces > 0 && self.chi == prescribed
    }
}

/// Extract `solid` at `samples` per axis with all seven extractors and validate
/// each mesh. One entry per entry of `ALL_EXTRACTORS`, in that order, asserted.
fn measure(name: &str, solid: &Solid, samples: u32) -> Vec<Measured> {
    let shape = RuntimeShape3::new([samples; 3]).expect("P-140 grid fits u32");
    let cell_size = 2.0 * DOMAIN_HALF / f64::from(samples - 1);
    let origin = [-DOMAIN_HALF; 3];
    let cfg =
        ValidateConfig::from_cell_size(cell_size).expect("P-140 cell size is positive and finite");

    let mut visited: Vec<&'static str> = Vec::with_capacity(ALL_EXTRACTORS.len());
    let mut rows: Vec<Measured> = Vec::with_capacity(ALL_EXTRACTORS.len());
    let mut mesh = MeshBuffer::<f64>::new();

    isomesh::for_each_extractor!(f64, |entry, extractor| {
        mesh.reset();
        if let Err(e) = extractor.extract_into(solid, &shape, origin, cell_size, &mut mesh) {
            panic!("P-140: {entry} refused {name} at {samples}^3 samples: {e}");
        }
        let report = validate(&mesh, &cfg);
        visited.push(entry);
        rows.push(Measured {
            chi: report.euler_characteristic,
            genus: report.genus,
            vertices: report.vertices,
            faces: report.faces,
            components: report.components,
            boundary_edges: report.boundary_edges,
            non_manifold_edges: report.non_manifold_edges,
            non_manifold_vertices: report.non_manifold_vertices,
            duplicate_vertices: report.duplicate_vertices,
            manifold: report.is_manifold(),
        });
    });

    assert_eq!(
        visited,
        ALL_EXTRACTORS.to_vec(),
        "P-140: for_each_extractor! visited a different roster than ALL_EXTRACTORS \
         names, so every per-extractor column would be mislabelled"
    );
    rows
}

fn main() {
    if !std::env::args().any(|arg| arg == "--bench") {
        return;
    }
    let prereg = isomesh::experiment!("P-140");

    common::experiment::run(prereg, |run| {
        let fixtures = fixtures();
        let sound: Vec<Sound> = fixtures.iter().map(soundness).collect();
        let extractors = ALL_EXTRACTORS.len();
        let ladder: String = LADDER
            .iter()
            .map(u32::to_string)
            .collect::<Vec<String>>()
            .join("|");

        // ── control 1: the ladder is a ladder ───────────────────────────────
        assert!(
            LADDER.len() >= 4 && LADDER.windows(2).all(|w| w[0] < w[1]),
            "VOID: the resolution ladder {ladder} is not four or more strictly \
             increasing rungs, so `first correct` does not name a well-ordered \
             thing and C2 is not a curve"
        );

        // ── control 2: genus above 1 ────────────────────────────────────────
        let above_one = fixtures
            .iter()
            .filter(|f| !f.is_control && f.genus > 1)
            .count();
        assert!(
            above_one >= 3,
            "VOID: only {above_one} prescribed fixtures have genus above 1, so the \
             suite is a sphere and a torus in disguise and neither C1 nor C2 says \
             anything the crate did not already assert"
        );

        // ── control 3: prescribed_chi is not a constant ─────────────────────
        let mut distinct_chi: BTreeSet<i64> = BTreeSet::new();
        for f in &fixtures {
            distinct_chi.insert(f.chi());
        }
        assert!(
            distinct_chi.len() >= 4,
            "VOID: the fixtures prescribe only {} distinct chi values ({:?}), so \
             `chi_agreement` could read true for an extractor that emits the same \
             topology every time",
            distinct_chi.len(),
            distinct_chi
        );

        // ── control 4: the prescription is sound by construction ────────────
        println!(
            "{:>16} {:>5} {:>5} {:>4} {:>9} {:>9} {:>9} {:>9} {:>9} {:>7}",
            "construction",
            "genus",
            "chi",
            "gens",
            "solid",
            "separate",
            "window",
            "merge",
            "margin",
            "clear"
        );
        for (f, s) in fixtures.iter().zip(&sound) {
            println!(
                "{:>16} {:>5} {:>5} {:>4} {:>9.4} {:>9.4} {:>9.4} {:>9.4} {:>9.4} {:>7.4}",
                f.name,
                f.genus,
                f.chi(),
                s.generators,
                s.solid_width,
                s.separation,
                s.window,
                s.merge_headroom,
                s.domain_margin,
                s.feature_clearance
            );
            assert_eq!(
                s.generators as i64, f.genus,
                "VOID: {} exhibits {} independent cycles and prescribes genus {}, \
                 so `prescribed_genus` is not what the construction fixes and \
                 `prescribed_chi` is arithmetic over the wrong number",
                f.name, s.generators, f.genus
            );
            assert!(
                s.solid_width > 0.0,
                "VOID: {} has a solid feature of width {}, so two parts of the \
                 solid that must stay apart touch and the construction's genus is \
                 not the genus it prescribes",
                f.name,
                s.solid_width
            );
            assert!(
                s.separation > 0.0,
                "VOID: {} has two non-incident tubes {} from merging, so the \
                 neighbourhood is not regular and genus = b1(G) does not hold",
                f.name,
                s.separation
            );
            assert!(
                s.window > 0.0,
                "VOID: {} has a prescribed void of clear width {}, so at least one \
                 handle's window is closed and the prescribed genus is an \
                 overcount",
                f.name,
                s.window
            );
            assert!(
                s.merge_headroom > 0.0,
                "VOID: {} has two incident tubes merging {} past half the shorter \
                 edge, so the merge regions at one edge's two ends run into each \
                 other and the neighbourhood is not regular",
                f.name,
                s.merge_headroom
            );
            assert!(
                s.domain_margin > 0.0,
                "VOID: {} reaches to within {} of the [-2, 2]^3 sampling box, so \
                 the box cuts the surface, the mesh is not closed and chi measures \
                 a surface with boundary",
                f.name,
                s.domain_margin
            );
        }

        // The graph fixtures state their genus a second way, and the two must
        // agree: b1 = E − V + 1, counted off the edge list.
        for f in &fixtures {
            match &f.solid {
                Solid::ThickenedGraph { graph, .. } => {
                    let b1 = graph.edges.len() as i64 - graph.nodes.len() as i64 + 1;
                    assert_eq!(
                        b1,
                        f.genus,
                        "VOID: {} has V={} E={}, so b1 = {} while it prescribes \
                         genus {} — the cycle list and the edge list disagree \
                         about what was built",
                        f.name,
                        graph.nodes.len(),
                        graph.edges.len(),
                        b1,
                        f.genus
                    );
                }
                Solid::DrilledBall { bores, .. } => {
                    assert_eq!(
                        bores.len() as i64,
                        f.genus,
                        "VOID: {} drills {} bores and prescribes genus {}, and \
                         chi = 2 - 2g is arithmetic over the bore count",
                        f.name,
                        bores.len(),
                        f.genus
                    );
                }
                Solid::Sphere(_) | Solid::Torus(_) => {}
            }
        }

        // ── the sweep ───────────────────────────────────────────────────────
        // sweep[fixture][rung][extractor]
        let mut sweep: Vec<Vec<Vec<Measured>>> = Vec::with_capacity(fixtures.len());
        for f in &fixtures {
            let mut per_rung: Vec<Vec<Measured>> = Vec::with_capacity(LADDER.len());
            for samples in LADDER {
                let rows = measure(f.name, &f.solid, samples);
                let seen: Vec<String> = rows
                    .iter()
                    .map(|m| match m.genus {
                        Some(g) => format!("{:>4}/g{g}", m.chi),
                        None => format!("{:>4}/--", m.chi),
                    })
                    .collect();
                println!(
                    "{:>16} {:>3}^3 want chi {:>3}  {}",
                    f.name,
                    samples,
                    f.chi(),
                    seen.join(" ")
                );
                per_rung.push(rows);
            }
            sweep.push(per_rung);
        }

        // first[fixture][extractor]: the smallest rung reproducing the
        // prescribed chi, or None for "never within the ladder".
        let first: Vec<Vec<Option<u32>>> = fixtures
            .iter()
            .enumerate()
            .map(|(fi, f)| {
                (0..extractors)
                    .map(|xi| {
                        LADDER
                            .iter()
                            .copied()
                            .enumerate()
                            .find(|(ri, _)| sweep[fi][*ri][xi].reproduces(f.chi()))
                            .map(|(_, samples)| samples)
                    })
                    .collect()
            })
            .collect();

        let fixture_distinct: Vec<usize> = first.iter().map(|per| distinct_rungs(per)).collect();
        let mut overall: Vec<Option<u32>> = Vec::new();
        for per in &first {
            overall.extend_from_slice(per);
        }
        let distinct_overall = distinct_rungs(&overall);

        // ── control 5: the instrument reads a known chi ─────────────────────
        let finest = LADDER.len() - 1;
        for (fi, f) in fixtures.iter().enumerate().filter(|(_, f)| f.is_control) {
            assert_eq!(
                f.declared,
                Some(f.chi()),
                "VOID: this harness prescribes chi {} for the calibration arm {} \
                 and the crate declares {:?}, so the two rulers disagree before a \
                 single mesh has been extracted",
                f.chi(),
                f.name,
                f.declared
            );
            for (xi, entry) in ALL_EXTRACTORS.iter().enumerate() {
                let m = sweep[fi][finest][xi];
                assert!(
                    m.reproduces(f.chi()),
                    "VOID: {entry} reads chi {} over {} faces on {} at {}^3, where \
                     the crate itself declares chi {}. The instrument cannot read \
                     a chi the crate already asserts, so nothing licenses reading \
                     it off the six prescribed solids",
                    m.chi,
                    m.faces,
                    f.name,
                    LADDER[finest],
                    f.chi()
                );
            }
        }

        // ── control 6: chi_agreement could have been false ──────────────────
        // Scoped to the whole sweep. C2 coming out uninformative is a
        // registered outcome ("still worth one line"), so an assert on the
        // coarsest rung alone would turn a registered null into a crash. What
        // must not happen is that the column never reads `false` at all.
        let misses: usize = fixtures
            .iter()
            .enumerate()
            .filter(|(_, f)| !f.is_control)
            .map(|(fi, f)| {
                (0..LADDER.len())
                    .map(|ri| {
                        (0..extractors)
                            .filter(|xi| !sweep[fi][ri][*xi].reproduces(f.chi()))
                            .count()
                    })
                    .sum::<usize>()
            })
            .sum();
        assert!(
            misses > 0,
            "VOID: `chi_agreement` reads true on every one of the {} prescribed \
             fixture x rung x extractor readings in this sweep, so it is a column \
             that could not have read false and nothing here shows the instrument \
             able to say no (M-44)",
            6 * LADDER.len() * extractors
        );

        // ── control 7: the sweep meshed something ───────────────────────────
        for (fi, f) in fixtures.iter().enumerate() {
            for (xi, entry) in ALL_EXTRACTORS.iter().enumerate() {
                assert!(
                    sweep[fi][finest][xi].faces > 0,
                    "VOID: {entry} produced no face at all on {} at {}^3, so every \
                     `chi_agreement` this run reports for it was scored against an \
                     empty buffer whose chi is 0 by arithmetic",
                    f.name,
                    LADDER[finest]
                );
            }
        }

        // ── verdicts ────────────────────────────────────────────────────────
        // C1 is scored on the prescribed fixtures; the calibration arms are
        // control 5's business, asserted above rather than reported here.
        let c1 = fixtures
            .iter()
            .enumerate()
            .filter(|(_, f)| !f.is_control)
            .all(|(fi, _)| first[fi].iter().all(Option::is_some));
        // C2's falsifier is "all extractors converging at the same resolution",
        // so it holds when at least one prescribed fixture discriminates.
        let c2 = fixtures
            .iter()
            .enumerate()
            .filter(|(_, f)| !f.is_control)
            .any(|(fi, _)| fixture_distinct[fi] >= 2);

        println!(
            "\nC2 curve — first rung reproducing the prescribed chi (never = not within {ladder})"
        );
        for (fi, f) in fixtures.iter().enumerate() {
            for (xi, entry) in ALL_EXTRACTORS.iter().enumerate() {
                println!(
                    "{:>16} g{} chi {:>3}  {:>28}  {}",
                    f.name,
                    f.genus,
                    f.chi(),
                    entry,
                    first[fi][xi].map_or_else(|| String::from("never"), |n| format!("{n}^3"))
                );
            }
            println!(
                "{:>16} distinct first-correct rungs across the seven extractors: {}",
                f.name, fixture_distinct[fi]
            );
        }
        println!(
            "\nC1 {c1} (every extractor converged on every prescribed fixture)  \
             C2 {c2} ({distinct_overall} distinct first-correct rungs over the sweep)"
        );

        // ── the rows ────────────────────────────────────────────────────────
        for (fi, f) in fixtures.iter().enumerate() {
            let s = &sound[fi];
            let prescribed = f.chi();
            for (ri, samples) in LADDER.iter().copied().enumerate() {
                let cells = u64::from(samples - 1).pow(3);
                let cell_size = 2.0 * DOMAIN_HALF / f64::from(samples - 1);
                for (xi, entry) in ALL_EXTRACTORS.iter().enumerate() {
                    let m = sweep[fi][ri][xi];
                    run.record(&[
                        ("construction", f.name.to_string()),
                        ("prescribed_genus", f.genus.to_string()),
                        ("prescribed_chi", prescribed.to_string()),
                        ("measured_chi", m.chi.to_string()),
                        ("chi_agreement", m.reproduces(prescribed).to_string()),
                        ("resolution", samples.to_string()),
                        ("cells", cells.to_string()),
                        ("nonmanifold_edges", m.non_manifold_edges.to_string()),
                        ("extractors_tested", extractors.to_string()),
                        ("c1_holds", c1.to_string()),
                        ("c2_holds", c2.to_string()),
                        // ── extras (M-273) ──
                        ("extractor", (*entry).to_string()),
                        ("derivation", f.derivation.to_string()),
                        ("is_control", f.is_control.to_string()),
                        ("cell_size", format!("{cell_size:.6}")),
                        ("chi_error", (m.chi - prescribed).to_string()),
                        ("converged", first[fi][xi].is_some().to_string()),
                        (
                            "first_correct_resolution",
                            first[fi][xi].map_or_else(|| String::from("never"), |n| n.to_string()),
                        ),
                        (
                            "fixture_distinct_first_correct",
                            fixture_distinct[fi].to_string(),
                        ),
                        (
                            "distinct_first_correct_overall",
                            distinct_overall.to_string(),
                        ),
                        ("resolution_ladder", ladder.clone()),
                        ("generators", s.generators.to_string()),
                        ("feature_clearance", format!("{:.6}", s.feature_clearance)),
                        (
                            "feature_cells",
                            format!("{:.6}", s.feature_clearance / cell_size),
                        ),
                        ("solid_width", format!("{:.6}", s.solid_width)),
                        ("separation", format!("{:.6}", s.separation)),
                        ("window", format!("{:.6}", s.window)),
                        ("merge_headroom", format!("{:.6}", s.merge_headroom)),
                        ("domain_margin", format!("{:.6}", s.domain_margin)),
                        (
                            "measured_genus",
                            m.genus
                                .map_or_else(|| String::from("none"), |g| g.to_string()),
                        ),
                        ("mesh_components", m.components.to_string()),
                        ("vertices", m.vertices.to_string()),
                        ("faces", m.faces.to_string()),
                        ("boundary_edges", m.boundary_edges.to_string()),
                        ("non_manifold_vertices", m.non_manifold_vertices.to_string()),
                        ("duplicate_vertices", m.duplicate_vertices.to_string()),
                        ("is_manifold", m.manifold.to_string()),
                    ]);
                }
            }
        }
    });
}
