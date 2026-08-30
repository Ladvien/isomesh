//! E-321 — intrinsic Delaunay flipping improves the connectivity and moves not
//! one vertex.
//!
//! ```bash
//! cd bevy_isomesh && cargo run --example intrinsic_flips --release
//! ```
//!
//! Keys: `F` before/after, `V` the vertex markers, `S` the sliver overlay,
//! `1`-`7` field, `[` `]` resolution, `W` wireframe (**on by default here**),
//! `N` normals, `G` the sampling box, `H` HUD, `F12` screenshot, `Esc` quit.
//!
//! `ISOMESH_FIELD=7` reaches the eighth field, which no digit key does — the
//! harness maps `1`-`7` onto indices 0-6 (`common/mod.rs:634-644`).
//! `ISOMESH_SAMPLES` must name 17, 25 or 33: those are the only resolutions
//! `crates/isomesh/golden_hashes.json` commits a hash at, and with no committed
//! hash there is no `hashes_moved` to report.
//!
//! Under `ISOMESH_CAPTURE` it needs no keyboard. The arm alternates every ten
//! captured frames, so the clip is the connectivity changing while the yellow
//! vertex crosses stay exactly where they are — which is the finding. The period
//! is twenty captured frames and `record_gif.sh` defaults to eighty, so the clip
//! is four whole cycles and loops on itself.
//!
//! ```bash
//! ISOMESH_CAPTURE_FRAMES=80 ISOMESH_CAPTURE_EVERY=2 \
//!   scripts/record_gif.sh intrinsic_flips docs/gifs/e321.gif
//! ```
//!
//! # What is on screen
//!
//! One surface, drawn shaded and never swapped: the mesh
//! `MarchingCubes::<f64>::new()` emitted, which is the surface **both** arms
//! triangulate. Over it, one of two wireframes — the extrinsic connectivity
//! Marching Cubes produced, or the intrinsic Delaunay connectivity a flipping
//! loop reached from it. `F` swaps them.
//!
//! The yellow crosses are the vertex set, drawn from **one** buffer that both
//! arms index. There is no second copy in this file to get out of step: `arm_mesh`
//! takes the same `positions` slice for either arm, and the crosses are drawn
//! from that same slice by a system that never looks at which arm is showing. So
//! the alternation is the edges moving over a vertex set that provably cannot
//! follow them.
//!
//! Red edges are slivers — triangles whose smallest angle is under
//! [`SLIVER_DEGREES`]. Watch them thin out, and watch how many stay.
//!
//! # What the number means
//!
//! `worst-decile min angle` is C1's statistic: the mean of the worst tenth of
//! the per-triangle minimum angles, in degrees, over the triangles with positive
//! Heron area. The 10th percentile and the global minimum sit beside it so the
//! verdict can be re-read under all three instruments rather than resting on one.
//!
//! For the **after** arm every one of those angles is computed from the intrinsic
//! edge lengths, not from the vertex positions. That is the whole mechanism: a
//! flip takes its new edge's length from the two triangles unfolded into the
//! plane along their shared edge (law of cosines on the intrinsic quadrilateral),
//! never from the chord between the two vertices in `R^3`. Take the chord instead
//! and the result is an ordinary extrinsic remesh of a *different* surface. The
//! panel's `chords that are not geodesic` row is that difference, counted: it is
//! `0` before flipping by construction, and non-zero after.
//!
//! `vertex_positions_moved` is a live bit-for-bit comparison of the extraction's
//! position buffer against a copy taken before the flipping loop ran.
//! `hashes_moved` is `isomesh::validate::mesh_hash` over that buffer against the
//! committed hash **read out of `crates/isomesh/golden_hashes.json`** — the file
//! `golden_hashes_are_unchanged` (`golden/tests.rs:59`) gates, not a
//! re-derivation of it.
//!
//! A zero that could not have been non-zero is not a measurement, so the control
//! is on the panel beside it: the same positions with the *flipped* indices
//! written into a `MeshBuffer` **does** move that hash. `mesh_hash` can therefore
//! see connectivity, and the flipping loop demonstrably changed some.
//!
//! `flips` is reported against its budget. The budget is
//! [`FLIP_BUDGET_PER_EDGE`] = 64 flips per interior edge, so `64 *
//! interior_edges`; the loop is Lawson's, seeded with every interior edge and
//! re-queueing all six corner slots of the two triangles a flip modified, which
//! is exactly the set of edges whose cotangent sum can have changed. An empty
//! queue is therefore a proved fixed point, and `non-Delaunay interior edges`
//! reading `0` after is what licenses calling this arm *the* intrinsic Delaunay
//! triangulation rather than a partial run.
//!
//! # What this shows, and what it does not
//!
//! This is P-175, and it is one clause held out of three.
//!
//! **C2 HELD.** Zero vertex positions moved and zero golden hashes moved, on all
//! 24 rows of `docs/experiments/p-175.csv`. That is not an absence of work, it is
//! a prediction confirmed: after `Intrinsic::build` returns, no method in this
//! file reads a position again, so a pinned vertex set is a property of the
//! mechanism rather than a result it happened to get. It matters here
//! specifically because this crate commits 216 golden hashes over exact vertex
//! bits and `P-61` has already moved 135 of them — any remesh that moves a
//! position moves a hash, and an intrinsic one cannot.
//!
//! **C1 FALSIFIED.** The bar was 10 degrees of gain on the worst decile on at
//! least four fields. Nothing came close: the largest gain in the sweep is under
//! two degrees, `p-175.csv` records `c1_fields = 0` on every row, and the panel
//! computes the gain live so a reader can watch it fall short. The flips are real
//! and the direction is right — slivers drop on every field that has any — the
//! size is not what was registered.
//!
//! **C3 FALSIFIED, by the falsifier the registration named in advance.** An
//! intrinsic flip leaves the extrinsic surface bit-identical, so if every
//! downstream consumer reads extrinsic triangles then nothing measurable changes
//! downstream. P-175 surveyed all eight places in `crates/isomesh/src/**` that
//! read triangle connectivity and found `0` that benefit. The one that reads
//! exactly what a flip improves — `validate::MeshReport::mean_ratio` — computes
//! it from the extrinsic cross product and edge lengths at
//! `validate.rs:917-939`, which no intrinsic flip can move. There is no cotangent
//! Laplacian, no geodesic solver and no parameterisation in this crate, and those
//! three are the operators an intrinsic triangulation exists to serve.
//!
//! So: a clean, cheap, hash-safe connectivity improvement with **nothing in this
//! crate that can consume it**. The picture is worth having anyway, because the
//! thing it makes visible — edges moving while vertices do not — is the property,
//! and the property is what a future cotangent operator would be built on.
//!
//! One more honesty note about the picture itself. The after arm's edges are
//! drawn as straight chords between the same vertices. The intrinsic edge is a
//! geodesic *across* the surface and bends wherever it crosses an original edge,
//! so a flipped edge's drawn line is a chord of the thing being measured rather
//! than the thing itself. The panel quantifies the gap live, and on a smooth
//! field it is a small fraction of a cell. **That smallness is not an argument
//! that the two are the same thing.** The
//! difference the ticket is about is where a length comes *from*: a chord is a
//! measurement of `R^3` and moves with the vertices, an unfolded diagonal is a
//! measurement of the surface and does not. That is also why the shaded surface
//! is the **before** mesh at all times and the after arm contributes no shaded
//! geometry: drawing those chords as filled triangles would be drawing a
//! different polyhedron, which is precisely the error this measurement exists to
//! detect.

mod common;

use std::time::Instant;

use bevy::asset::RenderAssetUsages;
use bevy::prelude::*;
use bevy::render::mesh::{Indices, PrimitiveTopology};
use common::{Capture, CommonPlugin, DemoDomain, DemoMesh, DemoStats, OrbitCamera, ViewFlags};
use isomesh::fields::ReferenceField;
use isomesh::marching_cubes::MarchingCubes;
use isomesh::validate::mesh_hash;
use isomesh::{MeshBuffer, RuntimeShape3, Sdf};

// ════════════════════════════════════════════════════════════════════════════
// constants
// ════════════════════════════════════════════════════════════════════════════

/// The resolutions this demo can show, in samples per axis.
///
/// `golden.rs:73`'s `RESOLUTIONS`, and they have to be: `hashes_moved` is a
/// comparison against the committed fixture and the fixture exists at no other
/// resolution. A demo offering 41^3 would have a blank where its headline number
/// goes.
const LADDER: [u32; 3] = [17, 25, 33];

/// How many reference fields `for_each_reference_field!` expands to.
///
/// Asserted against the macro in [`the_field_order_is_the_fixtures_own`], so a
/// ninth field cannot appear without this following it.
const FIELD_COUNT: usize = 8;

/// The name the golden fixture gives `MarchingCubes::<f64>::new()`.
const GOLDEN_ALGORITHM: &str = "marching_cubes";

/// A triangle whose smallest angle is below this many degrees is a sliver.
///
/// P-175's `sliver_threshold_degrees`, and the same bar its registered vacuity
/// control sets for the worst decile.
const SLIVER_DEGREES: f64 = 15.0;

/// A cotangent sum below `-COTAN_TOLERANCE` fails the intrinsic Delaunay test.
///
/// The sum is dimensionless — `(b^2 + c^2 - a^2) / 4A` is a length squared over
/// a length squared — so a plain absolute tolerance is scale-free here.
const COTAN_TOLERANCE: f64 = 1e-12;

/// Flips allowed per interior edge before the loop is declared non-terminating.
const FLIP_BUDGET_PER_EDGE: u64 = 64;

/// A drawn chord this far from its intrinsic length, as a fraction of the cell
/// size, is not the geodesic it stands for.
///
/// Loose enough that the `before` arm reads exactly zero — its lengths came from
/// those very chords — and tight enough that a flipped edge over a curved
/// surface reads as differing.
const CHORD_TOLERANCE: f64 = 1e-9;

/// C1's bar: degrees of gain on the worst decile.
const C1_DEGREES: f64 = 10.0;

/// C1's other bar: how many fields had to clear [`C1_DEGREES`].
const C1_FIELDS: usize = 4;

/// `p-175.csv`'s `c1_fields`, identical on all 24 rows.
///
/// A sweep-level count, so a single row cannot recompute it — quoted as a
/// citation naming P-175, which is this repo's rule for a figure a demo cannot
/// re-measure (`game_dig.rs:2946-2952`). Held against the artefact by
/// [`the_cited_sweep_figures_are_what_p175_committed`].
const CITED_C1_FIELDS: u32 = 0;

/// `p-175.csv`'s `c3_consumers_benefiting`.
const CITED_C3_BENEFITING: u32 = 0;

/// `p-175.csv`'s `c3_consumers_surveyed`.
const CITED_C3_SURVEYED: u32 = 8;

/// The corner slot that means "this halfedge has no twin".
const NO_TWIN: u32 = u32::MAX;

/// Captured frames each arm is held for before the other one is shown.
///
/// Off `Capture::taken` rather than a clock, which is what the field is public
/// for (`common/mod.rs:288-291`): the sequence is then reproducible and its
/// length is the clip's. Ten gives a twenty-frame period against
/// `record_gif.sh`'s eighty, so the clip is four whole cycles and its last frame
/// is the state its first frame was in.
const ALTERNATE_FRAMES: u32 = 10;

/// Half the length of a vertex marker's arms, as a fraction of the cell size.
const MARKER_CELLS: f32 = 0.22;

// ════════════════════════════════════════════════════════════════════════════
// the intrinsic triangulation
// ════════════════════════════════════════════════════════════════════════════

/// The next corner slot around the same triangle.
const fn next_slot(h: u32) -> u32 {
    (h - h % 3) + (h % 3 + 1) % 3
}

/// The previous corner slot around the same triangle.
const fn prev_slot(h: u32) -> u32 {
    (h - h % 3) + (h % 3 + 2) % 3
}

/// Heron's area from three side lengths, in Kahan's stable ordering.
///
/// Zero rather than a `NaN` when the triangle inequality does not hold, which is
/// the value the formula tends to and the one every caller here tests for.
fn area_of(l0: f64, l1: f64, l2: f64) -> f64 {
    let mut side = [l0, l1, l2];
    side.sort_by(f64::total_cmp);
    let (c, b, a) = (side[0], side[1], side[2]);
    let inner = (a + (b + c)) * (c - (a - b)) * (c + (a - b)) * (a + (b - c));
    if inner > 0.0 {
        0.25 * inner.sqrt()
    } else {
        0.0
    }
}

/// One undirected edge, keyed for the pairing sort.
///
/// Field order is the sort order: the two endpoints, then the slot, so a run of
/// equal endpoints is contiguous and the pairing inside it is index-ordered and
/// therefore independent of the order the triangles arrived in.
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
struct EdgeKey {
    /// The lower-numbered endpoint.
    lo: u32,
    /// The higher-numbered endpoint.
    hi: u32,
    /// The corner slot this halfedge is.
    slot: u32,
    /// Whether the halfedge runs `lo` to `hi` rather than the other way.
    forward: bool,
}

/// How the edges were classified when the twins were paired.
///
/// Every one of these is a *precondition* of the flip rather than a defect
/// report: an intrinsic flip is defined only across an interior edge of an
/// oriented manifold, and an edge that is not one is simply not flippable.
#[derive(Clone, Copy, Default)]
struct EdgeCensus {
    /// Edges with exactly two oppositely-oriented halves. The flippable ones.
    interior: u64,
    /// Edges with one half.
    boundary: u64,
    /// Edges with three or more halves.
    non_manifold: u64,
    /// Edges whose two halves run the same way.
    inconsistently_oriented: u64,
    /// Edges whose two endpoints are the same vertex.
    zero_length: u64,
}

/// Angles of one triangulation, in degrees.
#[derive(Clone, Copy, Default)]
struct AngleStats {
    /// Mean of the worst tenth of the per-triangle minimum angles. C1's
    /// statistic, and `p-175.csv`'s `min_angle_statistic` names it.
    worst_decile_mean: f64,
    /// The 10th-percentile per-triangle minimum angle.
    percentile_10: f64,
    /// The smallest angle anywhere.
    global_min: f64,
    /// Triangles whose smallest angle is under [`SLIVER_DEGREES`].
    slivers: u64,
    /// Triangles with no positive Heron area, excluded from the three above.
    degenerate: u64,
}

/// What one flipping run did.
#[derive(Clone, Copy, Default)]
struct FlipReport {
    /// Flips performed.
    flips: u64,
    /// Flips the numerical guard on the convexity theorem refused. Expected 0.
    rejected: u64,
    /// [`FLIP_BUDGET_PER_EDGE`] times the interior edge count.
    budget: u64,
    /// Whether the budget ran out before the queue emptied.
    exhausted: bool,
}

/// How far the drawn chords are from the intrinsic lengths they stand for.
#[derive(Clone, Copy, Default)]
struct ChordGap {
    /// Undirected edges considered, which is every edge the wireframe draws.
    edges: u64,
    /// Edges whose chord differs from its intrinsic length.
    differing: u64,
    /// The largest such difference, in cells.
    worst_cells: f64,
}

/// A triangulation of a fixed vertex set, carried as lengths rather than points.
///
/// Slot `3t+i` is corner `i` of triangle `t`. `corner[3t+i]` is the vertex there,
/// `length[3t+i]` is the length of the halfedge running from corner `i` to corner
/// `i+1`, and `twin[3t+i]` is the slot on the other side of that halfedge.
///
/// **After [`Intrinsic::build`] returns, no method here reads a position.** That
/// is what makes the retriangulation intrinsic, and what makes
/// `vertex_positions_moved` a fact about the mechanism rather than an assertion
/// about it. The two methods that do take positions —
/// [`Intrinsic::chord_gap`] and [`Intrinsic::sliver_edges`] — are for *drawing*
/// and are called after the flipping has finished.
struct Intrinsic {
    /// The vertex at each corner slot.
    corner: Vec<u32>,
    /// The halfedge length leaving each corner slot.
    length: Vec<f64>,
    /// The slot across each halfedge, or [`NO_TWIN`].
    twin: Vec<u32>,
    /// How the edges classified at build time.
    census: EdgeCensus,
}

impl Intrinsic {
    /// The extrinsic mesh, read once, as an intrinsic triangulation.
    fn build(positions: &[[f64; 3]], indices: &[u32]) -> Self {
        let slots = indices.len();
        let mut corner = Vec::with_capacity(slots);
        corner.extend_from_slice(indices);
        let mut length = Vec::with_capacity(slots);
        let mut keyed: Vec<EdgeKey> = Vec::with_capacity(slots);

        for (tri, verts) in indices.as_chunks::<3>().0.iter().enumerate() {
            for i in 0..3 {
                let (u, v) = (verts[i], verts[(i + 1) % 3]);
                let from = positions[u as usize];
                let to = positions[v as usize];
                let d = [to[0] - from[0], to[1] - from[1], to[2] - from[2]];
                length.push((d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt());
                let forward = u < v;
                let (lo, hi) = if forward { (u, v) } else { (v, u) };
                keyed.push(EdgeKey {
                    lo,
                    hi,
                    slot: (tri * 3 + i) as u32,
                    forward,
                });
            }
        }
        keyed.sort_unstable();

        let mut twin = vec![NO_TWIN; slots];
        let mut census = EdgeCensus::default();
        let mut at = 0usize;
        while at < keyed.len() {
            let mut end = at + 1;
            while end < keyed.len()
                && keyed[end].lo == keyed[at].lo
                && keyed[end].hi == keyed[at].hi
            {
                end += 1;
            }
            let run = &keyed[at..end];
            if run[0].lo == run[0].hi {
                census.zero_length += 1;
            } else if run.len() == 1 {
                census.boundary += 1;
            } else if run.len() > 2 {
                census.non_manifold += 1;
            } else if run[0].forward == run[1].forward {
                census.inconsistently_oriented += 1;
            } else {
                twin[run[0].slot as usize] = run[1].slot;
                twin[run[1].slot as usize] = run[0].slot;
                census.interior += 1;
            }
            at = end;
        }

        Self {
            corner,
            length,
            twin,
            census,
        }
    }

    /// How many triangles there are.
    fn triangles(&self) -> usize {
        self.corner.len() / 3
    }

    /// The area of triangle `tri`, from its three lengths.
    fn area(&self, tri: u32) -> f64 {
        let base = (tri * 3) as usize;
        area_of(
            self.length[base],
            self.length[base + 1],
            self.length[base + 2],
        )
    }

    /// The cotangent of the angle at the corner opposite halfedge `h`.
    ///
    /// `cot = (b^2 + c^2 - a^2) / 4A`, so it comes from lengths alone and no
    /// angle is ever evaluated to decide a flip.
    fn cot_opposite(&self, h: u32) -> f64 {
        let a = self.length[h as usize];
        let b = self.length[next_slot(h) as usize];
        let c = self.length[prev_slot(h) as usize];
        (b * b + c * c - a * a) / (4.0 * area_of(a, b, c))
    }

    /// Whether `h` names an interior edge an intrinsic flip is defined across.
    ///
    /// Three preconditions, in the order they can fail: the edge has a twin and
    /// that twin is in a different triangle; the two triangles share no *other*
    /// edge, so the flip cannot fold one onto itself; and both have positive
    /// area, without which the two cotangents are infinities.
    fn flippable(&self, h: u32) -> bool {
        let o = self.twin[h as usize];
        if o == NO_TWIN || o / 3 == h / 3 {
            return false;
        }
        for other in [next_slot(h), prev_slot(h)] {
            let across = self.twin[other as usize];
            if across != NO_TWIN && across / 3 == o / 3 {
                return false;
            }
        }
        self.area(h / 3) > 0.0 && self.area(o / 3) > 0.0
    }

    /// `cot(alpha) + cot(beta)` across the interior edge at `h`.
    fn cot_sum(&self, h: u32) -> f64 {
        self.cot_opposite(h) + self.cot_opposite(self.twin[h as usize])
    }

    /// Flippable interior edges whose cotangent weight is negative.
    ///
    /// Equivalently, edges whose two opposite angles sum past `pi`. Read *after*
    /// the flipping loop this is the proof that the intrinsic Delaunay fixed
    /// point was reached; read before, it is the count of work available.
    fn negative_cotan_edges(&self) -> u64 {
        let mut count = 0;
        for h in 0..self.corner.len() as u32 {
            if self.twin[h as usize] < h || !self.flippable(h) {
                continue;
            }
            if self.cot_sum(h) < -COTAN_TOLERANCE {
                count += 1;
            }
        }
        count
    }

    /// Replace the interior edge at `h` with the other diagonal of its intrinsic
    /// quadrilateral. `false` when the unfolding is degenerate.
    ///
    /// **The length update is intrinsic, and that distinction is the whole
    /// ticket.** The two triangles are unfolded into the plane along their shared
    /// edge `ab`: `a = (0, 0)`, `b = (L, 0)`, apex `c` above at
    /// `x = (L^2 + |ca|^2 - |bc|^2) / 2L`, `y = 2A_abc / L`, apex `d` below by
    /// the same construction with the sign of `y` reversed. The new edge's length
    /// is `|cd|` **in that layout** — the length of a geodesic across the
    /// surface, not the length of the chord between those two vertices in `R^3`.
    ///
    /// The flip is always well defined where it is attempted: if
    /// `alpha + beta > pi` then `d` lies strictly inside the circumcircle of
    /// `abc`, the unfolded quadrilateral is strictly convex, and the new diagonal
    /// always crosses the old one. The guard below is the numerical remainder of
    /// that theorem and `FlipReport::rejected` records how often it fired.
    fn flip(&mut self, h: u32) -> bool {
        let o = self.twin[h as usize];
        let (hn, hp) = (next_slot(h), prev_slot(h));
        let (on, op) = (next_slot(o), prev_slot(o));

        // `h` runs a -> b in triangle (a, b, c); `o` runs b -> a in (b, a, d).
        let shared = self.length[h as usize];
        let l_bc = self.length[hn as usize];
        let l_ca = self.length[hp as usize];
        let l_ad = self.length[on as usize];
        let l_db = self.length[op as usize];

        let xc = (shared * shared + l_ca * l_ca - l_bc * l_bc) / (2.0 * shared);
        let yc = 2.0 * area_of(shared, l_bc, l_ca) / shared;
        let xd = (shared * shared + l_ad * l_ad - l_db * l_db) / (2.0 * shared);
        let yd = -2.0 * area_of(shared, l_ad, l_db) / shared;
        let (dx, dy) = (xc - xd, yc - yd);
        let l_new = (dx * dx + dy * dy).sqrt();

        if l_new <= 0.0
            || !l_new.is_finite()
            || area_of(l_ad, l_new, l_ca) <= 0.0
            || area_of(l_bc, l_new, l_db) <= 0.0
        {
            return false;
        }

        let corner_c = self.corner[hp as usize];
        let corner_d = self.corner[op as usize];
        let tw_hn = self.twin[hn as usize];
        let tw_on = self.twin[on as usize];

        // Triangle `h/3` becomes (a, d, c): slot `h` is a -> d, `hn` is the new
        // diagonal d -> c, and `hp` keeps c -> a untouched.
        self.length[h as usize] = l_ad;
        self.twin[h as usize] = tw_on;
        if tw_on != NO_TWIN {
            self.twin[tw_on as usize] = h;
        }
        self.corner[hn as usize] = corner_d;
        self.length[hn as usize] = l_new;
        self.twin[hn as usize] = on;

        // Triangle `o/3` becomes (b, c, d): slot `o` is b -> c, `on` is the new
        // diagonal c -> d, and `op` keeps d -> b untouched.
        self.length[o as usize] = l_bc;
        self.twin[o as usize] = tw_hn;
        if tw_hn != NO_TWIN {
            self.twin[tw_hn as usize] = o;
        }
        self.corner[on as usize] = corner_c;
        self.length[on as usize] = l_new;
        self.twin[on as usize] = hn;

        true
    }

    /// Lawson's flip loop, run to the intrinsic Delaunay fixed point.
    ///
    /// The queue is seeded with every interior edge, and on each flip all six
    /// corner slots of the two modified triangles are re-queued — which is
    /// exactly the set of edges whose cotangent sum can have changed, so an empty
    /// queue is a proved fixed point rather than a hope.
    ///
    /// The data structure is a Delta-complex, not a simplicial complex: after the
    /// build nothing is looked up by vertex pair again, so a flip producing a
    /// self-edge or a second edge between an already-joined pair is legal and
    /// needs no special case. Refusing those flips is what would make the result
    /// something other than *the* intrinsic Delaunay triangulation.
    fn delaunay_flip(&mut self) -> FlipReport {
        let budget = FLIP_BUDGET_PER_EDGE * self.census.interior.max(1);
        let mut queued = vec![false; self.corner.len()];
        let mut queue: std::collections::VecDeque<u32> =
            std::collections::VecDeque::with_capacity(self.corner.len());
        for h in 0..self.corner.len() as u32 {
            if self.twin[h as usize] != NO_TWIN && self.twin[h as usize] > h {
                queued[h as usize] = true;
                queue.push_back(h);
            }
        }

        let mut flips = 0u64;
        let mut rejected = 0u64;
        while let Some(h) = queue.pop_front() {
            queued[h as usize] = false;
            if flips >= budget {
                break;
            }
            if !self.flippable(h) || self.cot_sum(h) >= -COTAN_TOLERANCE {
                continue;
            }
            let o = self.twin[h as usize];
            // Slot indices, not edges: after the flip these same six slots are
            // the six corners of the two modified triangles.
            let touched = [h, next_slot(h), prev_slot(h), o, next_slot(o), prev_slot(o)];
            if !self.flip(h) {
                rejected += 1;
                continue;
            }
            flips += 1;
            for slot in touched {
                let across = self.twin[slot as usize];
                if across == NO_TWIN {
                    continue;
                }
                let canonical = slot.min(across);
                if !queued[canonical as usize] {
                    queued[canonical as usize] = true;
                    queue.push_back(canonical);
                }
            }
        }

        FlipReport {
            flips,
            rejected,
            budget,
            exhausted: flips >= budget,
        }
    }

    /// The smallest angle of triangle `tri`, in degrees, from its three lengths.
    ///
    /// `None` when Heron's area is not positive: such a triangle has no angles,
    /// and every statistic here excludes it rather than reporting a `NaN`.
    fn min_angle(&self, tri: u32) -> Option<f64> {
        let base = (tri * 3) as usize;
        let (l0, l1, l2) = (
            self.length[base],
            self.length[base + 1],
            self.length[base + 2],
        );
        if area_of(l0, l1, l2) <= 0.0 {
            return None;
        }
        // The angle opposite `opp`, between the sides `x` and `y`. Clamped
        // because a right angle can round its cosine a hair outside [-1, 1].
        let angle = |opp: f64, x: f64, y: f64| {
            ((x * x + y * y - opp * opp) / (2.0 * x * y))
                .clamp(-1.0, 1.0)
                .acos()
                .to_degrees()
        };
        Some(
            angle(l0, l1, l2)
                .min(angle(l1, l2, l0))
                .min(angle(l2, l0, l1)),
        )
    }

    /// Per-triangle minimum angles, in degrees, and what they say about slivers.
    fn angles(&self) -> AngleStats {
        let mut mins: Vec<f64> = Vec::with_capacity(self.triangles());
        let mut degenerate = 0u64;
        for tri in 0..self.triangles() as u32 {
            match self.min_angle(tri) {
                Some(a) => mins.push(a),
                None => degenerate += 1,
            }
        }
        if mins.is_empty() {
            return AngleStats {
                degenerate,
                ..AngleStats::default()
            };
        }
        mins.sort_by(f64::total_cmp);
        let n = mins.len();
        let decile = (n / 10).max(1);
        let sum: f64 = mins[..decile].iter().sum();
        AngleStats {
            worst_decile_mean: sum / decile as f64,
            percentile_10: mins[(n / 10).min(n - 1)],
            global_min: mins[0],
            slivers: mins.iter().filter(|a| **a < SLIVER_DEGREES).count() as u64,
            degenerate,
        }
    }

    /// Which triangles are slivers, agreeing with [`AngleStats::slivers`].
    fn sliver_triangles(&self) -> Vec<u32> {
        (0..self.triangles() as u32)
            .filter(|tri| self.min_angle(*tri).is_some_and(|a| a < SLIVER_DEGREES))
            .collect()
    }

    /// The three drawn chords of every sliver triangle, as a flat line list.
    fn sliver_edges(&self, positions: &[Vec3]) -> Vec<[Vec3; 2]> {
        let mut out = Vec::new();
        for tri in self.sliver_triangles() {
            let base = (tri * 3) as usize;
            let p = [
                positions[self.corner[base] as usize],
                positions[self.corner[base + 1] as usize],
                positions[self.corner[base + 2] as usize],
            ];
            out.push([p[0], p[1]]);
            out.push([p[1], p[2]]);
            out.push([p[2], p[0]]);
        }
        out
    }

    /// One canonical halfedge per undirected edge, in slot order.
    ///
    /// The interior edges are visited from their lower slot; a boundary edge has
    /// only one half and is visited from it. This is exactly the set of lines the
    /// wireframe draws.
    fn canonical_slots(&self) -> impl Iterator<Item = u32> + '_ {
        (0..self.corner.len() as u32).filter(move |h| {
            let o = self.twin[*h as usize];
            o == NO_TWIN || o > *h
        })
    }

    /// How far each drawn chord is from the intrinsic length it stands for.
    ///
    /// Zero `differing` before flipping, because those lengths *were* the chords.
    /// Non-zero after, because a flipped edge's length came from the unfolded
    /// quadrilateral instead — which is the one thing this demo has to be able to
    /// show rather than claim.
    fn chord_gap(&self, positions: &[[f64; 3]], cell_size: f64) -> ChordGap {
        let mut out = ChordGap::default();
        for h in self.canonical_slots() {
            out.edges += 1;
            let from = positions[self.corner[h as usize] as usize];
            let to = positions[self.corner[next_slot(h) as usize] as usize];
            let d = [to[0] - from[0], to[1] - from[1], to[2] - from[2]];
            let chord = (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt();
            let gap = (chord - self.length[h as usize]).abs();
            if gap > CHORD_TOLERANCE * cell_size {
                out.differing += 1;
            }
            out.worst_cells = out.worst_cells.max(gap / cell_size);
        }
        out
    }

    /// Edges a simplicial complex could not carry: self-edges, and second edges
    /// between a vertex pair that already has one.
    ///
    /// Both are legal in a Delta-complex and neither is refused, which is why
    /// they are counted rather than prevented. They are also why two chords can
    /// end up drawn on top of each other in the after arm.
    fn delta_complex(&self) -> (u64, u64) {
        let mut pairs: Vec<(u32, u32)> = Vec::new();
        let mut self_edges = 0u64;
        for h in self.canonical_slots() {
            let u = self.corner[h as usize];
            let v = self.corner[next_slot(h) as usize];
            if u == v {
                self_edges += 1;
                continue;
            }
            pairs.push((u.min(v), u.max(v)));
        }
        pairs.sort_unstable();
        let total = pairs.len();
        pairs.dedup();
        (self_edges, (total - pairs.len()) as u64)
    }
}

// ════════════════════════════════════════════════════════════════════════════
// the committed fixture
// ════════════════════════════════════════════════════════════════════════════

/// `crates/isomesh/golden_hashes.json`, compiled in.
///
/// The committed fixture `golden_hashes_are_unchanged` (`golden/tests.rs:59`)
/// gates against — not a re-derivation of it, which is what makes `hashes_moved`
/// movement in *the* fixture. `include_str!` rather than a runtime read, so the
/// baseline cannot depend on which directory the demo was launched from.
const GOLDEN_HASHES: &str = include_str!("../../crates/isomesh/golden_hashes.json");

/// The committed `marching_cubes` hash for one field at one resolution.
///
/// A hand-rolled scanner rather than a JSON parser, for the reason `golden.rs`
/// gives: the grammar is one line per entry, fixed key order, no nesting and no
/// escapes. The trailing comma after `samples` is part of the needle, so `17`
/// cannot match `170`, and both quotes around the algorithm are, so
/// `marching_cubes` cannot match `marching_cubes+decider`.
fn golden_hash(field: &str, samples: u32) -> Option<u64> {
    let needle = format!(
        "\"algorithm\":\"{GOLDEN_ALGORITHM}\",\"field\":\"{field}\",\"samples\":{samples},"
    );
    let line = GOLDEN_HASHES.lines().find(|l| l.contains(&needle))?;
    let at = line.find("\"hash\":\"")? + "\"hash\":\"".len();
    u64::from_str_radix(line.get(at..at + 16)?, 16).ok()
}

// ════════════════════════════════════════════════════════════════════════════
// the measurement
// ════════════════════════════════════════════════════════════════════════════

/// Everything one `(field, resolution)` produced, geometry included.
struct Measurement {
    /// The reference field's name, as the golden fixture spells it.
    field: &'static str,
    /// Samples per axis.
    samples: u32,
    /// The grid spacing, in world units.
    cell_size: f64,
    /// The one vertex buffer. Both arms index it and the markers are drawn from
    /// it, so there is no second copy in this demo to disagree with it.
    positions: Vec<Vec3>,
    /// The one normal buffer, likewise.
    normals: Vec<Vec3>,
    /// The connectivity Marching Cubes emitted.
    before_indices: Vec<u32>,
    /// The connectivity the flipping loop reached.
    after_indices: Vec<u32>,
    /// The domain the field was sampled over.
    domain: (Vec3, Vec3),
    /// How the edges classified.
    census: EdgeCensus,
    /// Angles before flipping.
    before: AngleStats,
    /// Angles after flipping, from the intrinsic lengths.
    after: AngleStats,
    /// Flippable edges failing the Delaunay test, before.
    non_delaunay_before: u64,
    /// The same after: `0` is the fixed point.
    non_delaunay_after: u64,
    /// What the flipping loop did.
    flips: FlipReport,
    /// Chord-versus-geodesic before flipping. `differing` is `0` here.
    gap_before: ChordGap,
    /// Chord-versus-geodesic after flipping.
    gap_after: ChordGap,
    /// Self-edges and doubled vertex pairs after flipping.
    delta_after: (u64, u64),
    /// Positions differing bit-for-bit from a copy taken before the flips.
    positions_moved: u64,
    /// `1` if `mesh_hash` over the untouched buffer left the committed hash.
    hashes_moved: u64,
    /// Whether positions, normals and indices are all bit-identical to the copy.
    extrinsic_identical: bool,
    /// C2's positive control: the same positions with the flipped indices move
    /// the hash. Without this, `hashes_moved = 0` proves nothing.
    control_hash_moved: bool,
    /// `mesh_hash` over the extraction, live.
    live_hash: u64,
    /// The hash `golden_hashes.json` commits for this row.
    committed_hash: u64,
    /// How long the Marching Cubes extraction took.
    extract_ms: f64,
    /// How long the flipping loop took, on this rebuild.
    flip_ms: f64,
    /// The sliver overlay's line list, before.
    sliver_edges_before: Vec<[Vec3; 2]>,
    /// The sliver overlay's line list, after.
    sliver_edges_after: Vec<[Vec3; 2]>,
}

impl Measurement {
    /// C1's arithmetic for this row: degrees gained on the worst decile.
    fn gain(&self) -> f64 {
        self.after.worst_decile_mean - self.before.worst_decile_mean
    }

    /// Half the length of a vertex marker's arms, in world units.
    fn marker(&self) -> f32 {
        MARKER_CELLS * self.cell_size as f32
    }
}

/// Every reference field's name, in the order the golden fixture holds them.
///
/// Driven off `for_each_reference_field!` rather than a second list, so the
/// index this demo shows and the row `golden_hashes.json` holds cannot drift.
fn field_names() -> Vec<&'static str> {
    let mut names = Vec::with_capacity(FIELD_COUNT);
    isomesh::for_each_reference_field!(f64, |name, field| {
        let _ = &field;
        names.push(name);
    });
    names
}

/// The name of reference field `index`.
fn field_name(index: usize) -> Option<&'static str> {
    field_names().get(index).copied()
}

/// Measure reference field `index` at `samples` per axis.
///
/// Selected by name rather than by a counter threaded through the expansion: the
/// macro expands its body once per field with the field bound to a concrete type,
/// so there is no dynamic dispatch and no second field list. Note the absence of
/// a `return` in the body — the macro is syntax, not a closure, and one there
/// would leave this function on `sphere` (`fields/mod.rs:198-210`).
fn measure(index: usize, samples: u32) -> Option<Measurement> {
    let wanted = field_name(index)?;
    let mut out = None;
    isomesh::for_each_reference_field!(f64, |name, field| {
        if name == wanted {
            out = measure_one(name, &field, samples);
        }
    });
    out
}

/// One `(field, resolution)`, end to end.
///
/// `f64` throughout, and that is load-bearing rather than a preference: the
/// golden hashes are taken over `MarchingCubes::<f64>` (`golden.rs:171`), so an
/// `f32` extraction here would be quoting a fixture it cannot reproduce — the
/// mistake `game_edit_tape_trim.rs:126-131` records. The Bevy mesh is `f32`
/// because a `Mesh` attribute is, and it is for drawing only; every number on the
/// panel comes from the `f64` buffer.
fn measure_one<F>(name: &'static str, field: &F, samples: u32) -> Option<Measurement>
where
    F: ReferenceField + Sdf<Scalar = f64>,
{
    // `samples - 1` in the denominator: the repo's convention is samples per
    // axis, not cells, and `golden.rs:164` uses exactly this expression.
    let (lo, hi) = field.domain();
    let cell_size = (hi[0] - lo[0]) / f64::from(samples - 1);
    let shape = RuntimeShape3::new([samples; 3]).ok()?;
    let mut mesh = MeshBuffer::<f64>::new();
    let started = Instant::now();
    MarchingCubes::<f64>::new()
        .extract(field, &shape, lo, cell_size, &mut mesh)
        .ok()?;
    let extract_ms = started.elapsed().as_secs_f64() * 1e3;
    if mesh.indices.is_empty() {
        return None;
    }
    let committed_hash = golden_hash(name, samples)?;

    // Independent copies, taken before the flipping run, so `positions_moved`
    // compares two buffers rather than one buffer with itself.
    let positions_before = mesh.positions.clone();
    let normals_before = mesh.normals.clone();
    let indices_before = mesh.indices.clone();

    let base = Intrinsic::build(&mesh.positions, &mesh.indices);
    let before = base.angles();
    let non_delaunay_before = base.negative_cotan_edges();
    let gap_before = base.chord_gap(&mesh.positions, cell_size);

    let mut flipped = Intrinsic::build(&mesh.positions, &mesh.indices);
    let started = Instant::now();
    let flips = flipped.delaunay_flip();
    let flip_ms = started.elapsed().as_secs_f64() * 1e3;
    let after = flipped.angles();
    let non_delaunay_after = flipped.negative_cotan_edges();
    let gap_after = flipped.chord_gap(&mesh.positions, cell_size);
    let delta_after = flipped.delta_complex();

    // C2's positive control: the flipped connectivity written back into a
    // `MeshBuffer` beside the untouched positions. If this does not move the
    // hash, neither `mesh_hash` nor the flipper can be believed (M-44).
    let mut control = MeshBuffer::<f64>::new();
    control.positions.clone_from(&mesh.positions);
    control.normals.clone_from(&mesh.normals);
    control.indices.clone_from(&flipped.corner);
    let control_hash_moved = mesh_hash(&control) != committed_hash;

    let positions_moved = moved_positions(&mesh.positions, &positions_before);
    let normals_identical = mesh.normals.len() == normals_before.len()
        && mesh.normals.iter().zip(&normals_before).all(|(now, was)| {
            now.iter()
                .zip(was.iter())
                .all(|(x, y)| x.to_bits() == y.to_bits())
        });
    let extrinsic_identical =
        positions_moved == 0 && normals_identical && mesh.indices == indices_before;
    let live_hash = mesh_hash(&mesh);

    let positions: Vec<Vec3> = mesh
        .positions
        .iter()
        .map(|p| Vec3::new(p[0] as f32, p[1] as f32, p[2] as f32))
        .collect();
    let normals: Vec<Vec3> = mesh
        .normals
        .iter()
        .map(|n| Vec3::new(n[0] as f32, n[1] as f32, n[2] as f32))
        .collect();
    let sliver_edges_before = base.sliver_edges(&positions);
    let sliver_edges_after = flipped.sliver_edges(&positions);

    Some(Measurement {
        field: name,
        samples,
        cell_size,
        positions,
        normals,
        before_indices: mesh.indices.clone(),
        after_indices: flipped.indices(),
        domain: (
            Vec3::new(lo[0] as f32, lo[1] as f32, lo[2] as f32),
            Vec3::new(hi[0] as f32, hi[1] as f32, hi[2] as f32),
        ),
        census: base.census,
        before,
        after,
        non_delaunay_before,
        non_delaunay_after,
        flips,
        gap_before,
        gap_after,
        delta_after,
        positions_moved,
        hashes_moved: u64::from(live_hash != committed_hash),
        extrinsic_identical,
        control_hash_moved,
        live_hash,
        committed_hash,
        extract_ms,
        flip_ms,
        sliver_edges_before,
        sliver_edges_after,
    })
}

impl Intrinsic {
    /// The flipped connectivity, as a flat index buffer.
    fn indices(&self) -> Vec<u32> {
        self.corner.clone()
    }
}

/// Positions differing bit-for-bit, plus any length difference.
///
/// Bits rather than values, because `+0.0 == -0.0` compares equal and hashes
/// differently — the distinction `mesh_hash` is built on
/// (`validate/mesh_hash.rs:66-68`).
fn moved_positions(now: &[[f64; 3]], was: &[[f64; 3]]) -> u64 {
    now.iter()
        .zip(was)
        .filter(|(a, b)| {
            a.iter()
                .zip(b.iter())
                .any(|(x, y)| x.to_bits() != y.to_bits())
        })
        .count() as u64
        + (now.len() as u64).abs_diff(was.len() as u64)
}

// ════════════════════════════════════════════════════════════════════════════
// the app
// ════════════════════════════════════════════════════════════════════════════

/// The vertex markers, in their own group so they are never lost behind the
/// surface they sit exactly on.
#[derive(Default, Reflect, GizmoConfigGroup)]
struct VertexGizmos;

/// The sliver overlay, biased harder still: a sliver's edges lie on the same
/// triangles the wireframe draws, so at the shared bias they z-fight with it and
/// the highlight looks intermittent (`manifold_check.rs:185-191`).
#[derive(Default, Reflect, GizmoConfigGroup)]
struct SliverGizmos;

/// Which of the two mesh entities a handle belongs to.
#[derive(Component, Clone, Copy, Debug, PartialEq, Eq)]
enum Panel {
    /// The surface both arms triangulate: shaded, and always the extrinsic mesh.
    Surface,
    /// The arm being shown: wireframed by the harness, and deliberately given no
    /// material. The after arm's triangles are chords *across* the surface rather
    /// than pieces of it, so shading them would be shading a different
    /// polyhedron.
    Wire,
}

/// What the demo is showing.
#[derive(Resource)]
struct Demo {
    /// Index into `for_each_reference_field!`.
    field: usize,
    /// Index into [`LADDER`].
    rung: usize,
    /// Whether the flipped arm is the one on screen.
    after: bool,
    /// Whether the vertex markers are drawn.
    vertices_shown: bool,
    /// Whether the sliver overlay is drawn.
    slivers_shown: bool,
}

/// The last measurement, or `None` before the first rebuild.
#[derive(Resource, Default)]
struct Measured(Option<Measurement>);

/// The two arms as mesh assets, built from one position buffer.
#[derive(Resource, Default)]
struct Arms {
    /// The extrinsic connectivity.
    before: Handle<Mesh>,
    /// The intrinsic Delaunay connectivity.
    after: Handle<Mesh>,
}

/// The shaded surface's material.
#[derive(Resource)]
struct SurfaceMaterial(Handle<StandardMaterial>);

fn main() {
    App::new()
        .add_plugins(DefaultPlugins.set(WindowPlugin {
            primary_window: Some(Window {
                title: "isomesh - E-321 intrinsic flips".into(),
                ..default()
            }),
            ..default()
        }))
        .add_plugins(CommonPlugin)
        .init_gizmo_group::<VertexGizmos>()
        .init_gizmo_group::<SliverGizmos>()
        .insert_resource(Demo {
            field: 0,
            rung: starting_rung(),
            after: false,
            vertices_shown: true,
            slivers_shown: true,
        })
        .init_resource::<Measured>()
        .init_resource::<Arms>()
        .add_systems(Startup, setup)
        // `PreUpdate`, for E-306's and E-312's reason: the harness's `update_hud`
        // renders `DemoStats` and its `capture_sequence` advances
        // `Capture::taken`, both in `Update` with no ordering against an
        // example's own systems. In `Update` the panel would render a frame-old
        // readout beside a current picture, and `controls` would read `taken` on
        // either side of its increment -- for a demo whose whole claim is "these
        // numbers describe this picture", the one defect that matters.
        //
        // After `InputSystems` so a keypress is seen in the frame it happened.
        .add_systems(
            PreUpdate,
            (controls, rebuild, show, report)
                .chain()
                .after(bevy::input::InputSystems),
        )
        .add_systems(Update, (draw_vertices, draw_slivers))
        .run();
}

/// The rung `ISOMESH_SAMPLES` asks for, defaulting to 25.
///
/// 25 rather than the smallest or the largest: 17 has too few triangles for the
/// alternation to read as a retriangulation and 33 is dense enough that the
/// vertex markers merge into the wireframe. A value off [`LADDER`] is refused
/// rather than rounded, because there is no committed hash beside it and
/// `hashes_moved` is this demo's headline.
fn starting_rung() -> usize {
    let Some(samples) = common::samples_override() else {
        return 1;
    };
    match LADDER.iter().position(|n| *n == samples) {
        Some(rung) => rung,
        None => {
            // `eprintln!` rather than `error!`: this runs before `add_plugins`
            // has installed the log subscriber, so a `tracing` event here would
            // be swallowed.
            eprintln!("ISOMESH_SAMPLES={samples} is not one of {LADDER:?}; using 25");
            1
        }
    }
}

/// Complain if `ISOMESH_FIELD` names a field this demo does not have.
///
/// The harness parses that variable into `ViewFlags::field` with no range check
/// (`common/mod.rs:152-155`), so an out-of-range value would otherwise be quietly
/// ignored and the demo would open on `sphere` with no explanation.
fn check_pinned_field() {
    if let Ok(raw) = std::env::var("ISOMESH_FIELD")
        && !matches!(raw.trim().parse::<usize>(), Ok(index) if index < FIELD_COUNT)
    {
        error!("ISOMESH_FIELD={raw} is not one of 0..{FIELD_COUNT}");
    }
}

fn setup(
    mut commands: Commands,
    mut materials: ResMut<Assets<StandardMaterial>>,
    mut gizmo_config: ResMut<GizmoConfigStore>,
    mut camera: Query<&mut OrbitCamera>,
    mut flags: ResMut<ViewFlags>,
) {
    check_pinned_field();

    // Both arms are the same surface and the same vertices. Only the *edges*
    // differ, and edges are only visible in the wireframe -- the same reason
    // E-105 and E-106 start wireframed.
    flags.wireframe = true;

    for mut orbit in &mut camera {
        // Off both axes: down an axis a triangulation of a symmetric field reads
        // as a plane figure and the flips look like a 2D mesh being edited, which
        // is the one thing this picture must not suggest.
        orbit.yaw = 0.85;
        orbit.pitch = 0.42;
    }

    // Just enough negative bias to beat z-fighting, and no more. A marker sits
    // exactly on the surface so at the shared bias it flickers; at -1.0 every
    // vertex on the far side would show through and the picture would be a
    // hairball rather than a triangulation.
    let (verts, _) = gizmo_config.config_mut::<VertexGizmos>();
    verts.line.width = 2.6;
    verts.depth_bias = -0.25;

    let (slivers, _) = gizmo_config.config_mut::<SliverGizmos>();
    slivers.line.width = 3.4;
    slivers.depth_bias = -0.35;

    // Opaque, and darker than the harness's usual surface grey. The subject is
    // the green wireframe and the yellow markers over it, and it has to occlude
    // its own far side or the alternation is unreadable (`affine_rejection.rs`
    // darkens its surface for the same reason).
    commands.insert_resource(SurfaceMaterial(materials.add(StandardMaterial {
        base_color: Color::srgb(0.26, 0.29, 0.35),
        perceptual_roughness: 0.62,
        metallic: 0.04,
        ..default()
    })));

    // Spawned rather than assumed: `draw_domain` queries for it, and without one
    // the `G` toggle would silently do nothing. Filled in by the first rebuild.
    commands.spawn(DemoDomain {
        min: Vec3::splat(-1.0),
        max: Vec3::splat(1.0),
    });
}

/// Keys, and the self-driving alternation a capture needs.
///
/// Under `ISOMESH_CAPTURE` the keyboard is ignored outright: the recorder presses
/// nothing, and a subject that only changes on a keypress captures as a still.
fn controls(
    keys: Res<ButtonInput<KeyCode>>,
    capture: Res<Capture>,
    flags: Res<ViewFlags>,
    mut demo: ResMut<Demo>,
) {
    if capture.is_active() {
        // Off `Capture::taken`, so the sequence is in step with the frames that
        // reach the GIF rather than with wall-clock time. `taken` is advanced in
        // `Update` and read here in `PreUpdate`, so the swap lands on the frame
        // after the tenth capture -- a phase offset, not a period error.
        demo.after = (capture.taken / ALTERNATE_FRAMES) % 2 == 1;
        return;
    }
    if keys.just_pressed(KeyCode::KeyF) {
        demo.after = !demo.after;
    }
    if keys.just_pressed(KeyCode::KeyV) {
        demo.vertices_shown = !demo.vertices_shown;
    }
    if keys.just_pressed(KeyCode::KeyS) {
        demo.slivers_shown = !demo.slivers_shown;
    }
    if keys.just_pressed(KeyCode::BracketRight) {
        demo.rung = (demo.rung + 1).min(LADDER.len() - 1);
    }
    if keys.just_pressed(KeyCode::BracketLeft) {
        demo.rung = demo.rung.saturating_sub(1);
    }
    // The harness owns the digit keys and `ISOMESH_FIELD`; this demo only has to
    // refuse an index it has no field for.
    if flags.field < FIELD_COUNT {
        demo.field = flags.field;
    }
}

/// Extract, flip, and build both arms — only when the answer would change.
#[allow(clippy::too_many_arguments)]
fn rebuild(
    demo: Res<Demo>,
    mut measured: ResMut<Measured>,
    mut arms: ResMut<Arms>,
    mut meshes: ResMut<Assets<Mesh>>,
    material: Res<SurfaceMaterial>,
    mut commands: Commands,
    mut panels: Query<(&mut Mesh3d, &Panel)>,
    mut domain: Query<&mut DemoDomain>,
    mut camera: Query<&mut OrbitCamera>,
    mut flags: ResMut<ViewFlags>,
    mut last: Local<Option<(usize, usize)>>,
) {
    let key = (demo.field, demo.rung);
    if *last == Some(key) && !flags.remesh_requested {
        return;
    }
    *last = Some(key);
    flags.remesh_requested = false;

    let samples = LADDER[demo.rung];
    let Some(m) = measure(demo.field, samples) else {
        error!(
            "E-321: field {} at {samples} samples did not measure",
            demo.field
        );
        return;
    };

    // The HUD is the evidence and a headless capture has no HUD to read. One line
    // per rebuild, so `ISOMESH_CAPTURE` leaves the measurement in the log where a
    // script can hold it against `docs/experiments/p-175.csv` -- E-203 learned
    // this the hard way.
    info!(
        "{} at {}^3: {} triangles, {} flips, worst decile {:.6} -> {:.6} deg (gain {:.6}), \
         slivers {} -> {}, non-Delaunay {} -> {}, chords off their geodesic {} -> {} of {} \
         (worst {:.6} cell), vertex_positions_moved {}, hashes_moved {}, control_hash_moved {}",
        m.field,
        m.samples,
        m.before_indices.len() / 3,
        m.flips.flips,
        m.before.worst_decile_mean,
        m.after.worst_decile_mean,
        m.gain(),
        m.before.slivers,
        m.after.slivers,
        m.non_delaunay_before,
        m.non_delaunay_after,
        m.gap_before.differing,
        m.gap_after.differing,
        m.gap_after.edges,
        m.gap_after.worst_cells,
        m.positions_moved,
        m.hashes_moved,
        m.control_hash_moved,
    );

    // One `positions` slice, one `normals` slice, two index buffers. The claim
    // this demo makes is that the vertex set is shared, and this is that claim
    // written as code rather than asserted in a comment.
    arms.before = meshes.add(arm_mesh(&m.positions, &m.normals, &m.before_indices));
    arms.after = meshes.add(arm_mesh(&m.positions, &m.normals, &m.after_indices));

    for mut box_ in &mut domain {
        box_.min = m.domain.0;
        box_.max = m.domain.1;
    }
    // Derived from the domain, not hardcoded: `gyroid`'s box is three and a half
    // times the compact fields' and a fixed radius put E-304's camera *inside*
    // its own subject (`BACKLOG_ARCHIVE.md:157`). At Bevy's default 45 degrees a
    // distance of 1.45 widths shows 0.60 of a width either side of the focus
    // against a half-extent of 0.50, so the mesh fills the frame with a margin
    // rather than touching the edges.
    let width = m.domain.1.x - m.domain.0.x;
    for mut orbit in &mut camera {
        orbit.focus = (m.domain.0 + m.domain.1) * 0.5;
        orbit.radius = width * 1.45;
    }

    let showing = if demo.after {
        arms.after.clone()
    } else {
        arms.before.clone()
    };
    if panels.is_empty() {
        commands.spawn((
            Mesh3d(arms.before.clone()),
            MeshMaterial3d(material.0.clone()),
            Transform::IDENTITY,
            Panel::Surface,
        ));
        // No `MeshMaterial3d`, on purpose: `Mesh3d` requires only `Transform`
        // (`bevy_mesh/components.rs:101`), so an entity with no material
        // contributes no shaded geometry while `DemoMesh` still gets it a
        // wireframe and normals from the harness.
        commands.spawn((Mesh3d(showing), Transform::IDENTITY, DemoMesh, Panel::Wire));
    } else {
        for (mut mesh, panel) in &mut panels {
            let wanted = match panel {
                Panel::Surface => arms.before.clone(),
                Panel::Wire => showing.clone(),
            };
            if mesh.0 != wanted {
                mesh.0 = wanted;
            }
        }
    }

    measured.0 = Some(m);
}

/// One arm as a Bevy mesh: the shared positions and normals, that arm's indices.
fn arm_mesh(positions: &[Vec3], normals: &[Vec3], indices: &[u32]) -> Mesh {
    let mut mesh = Mesh::new(
        PrimitiveTopology::TriangleList,
        RenderAssetUsages::default(),
    );
    mesh.insert_attribute(
        Mesh::ATTRIBUTE_POSITION,
        positions
            .iter()
            .map(Vec3::to_array)
            .collect::<Vec<[f32; 3]>>(),
    );
    mesh.insert_attribute(
        Mesh::ATTRIBUTE_NORMAL,
        normals
            .iter()
            .map(Vec3::to_array)
            .collect::<Vec<[f32; 3]>>(),
    );
    mesh.insert_indices(Indices::U32(indices.to_vec()));
    mesh
}

/// Point the wire entity at the arm the toggle selects.
///
/// Written only when it differs, the same change-driven-write rule the harness
/// states for its `Text`: an unconditional `mesh.0 = …` marks the component
/// changed every frame and Bevy's mesh extraction is change-driven, so a toggle
/// nobody pressed would become per-frame work.
fn show(demo: Res<Demo>, arms: Res<Arms>, mut panels: Query<(&mut Mesh3d, &Panel)>) {
    let wanted = if demo.after {
        &arms.after
    } else {
        &arms.before
    };
    for (mut mesh, panel) in &mut panels {
        if *panel == Panel::Wire && mesh.0 != *wanted {
            mesh.0 = wanted.clone();
        }
    }
}

/// The vertex set, drawn from the one buffer both arms index.
///
/// This system does not know which arm is showing, and that is the point: there
/// is no path by which the markers could follow a retriangulation.
///
/// Uncapped, unlike the harness's `draw_normals` and its `MAX_LINES` stride
/// (`common/mod.rs:795-826`). A strided marker set would show *some* vertices
/// pinned, which is not the claim; and the worst case here is `noise_cavity` at
/// 33 samples, whose vertex count (`p-175.csv`'s `vertices` column) comes to
/// three times that many lines — the same order as the cap — on a demo where no
/// clause reads a clock. `V` turns them off.
fn draw_vertices(demo: Res<Demo>, measured: Res<Measured>, mut gizmos: Gizmos<VertexGizmos>) {
    /// The pinned vertex set.
    const YELLOW: Color = Color::srgb(1.0, 0.95, 0.20);

    if !demo.vertices_shown {
        return;
    }
    let Some(m) = &measured.0 else {
        return;
    };
    let r = m.marker();
    for p in &m.positions {
        gizmos.line(*p - Vec3::X * r, *p + Vec3::X * r, YELLOW);
        gizmos.line(*p - Vec3::Y * r, *p + Vec3::Y * r, YELLOW);
        gizmos.line(*p - Vec3::Z * r, *p + Vec3::Z * r, YELLOW);
    }
}

/// The sliver triangles of the arm on screen.
fn draw_slivers(demo: Res<Demo>, measured: Res<Measured>, mut gizmos: Gizmos<SliverGizmos>) {
    /// Triangles under [`SLIVER_DEGREES`].
    const RED: Color = Color::srgb(1.0, 0.30, 0.26);

    if !demo.slivers_shown {
        return;
    }
    let Some(m) = &measured.0 else {
        return;
    };
    let edges = if demo.after {
        &m.sliver_edges_after
    } else {
        &m.sliver_edges_before
    };
    for [a, b] in edges {
        gizmos.line(*a, *b, RED);
    }
}

/// `yes` or a shouted `NO`, for a boolean a reader has to be able to check.
fn yes_no(value: bool) -> &'static str {
    if value { "yes" } else { "NO" }
}

/// The panel.
fn report(demo: Res<Demo>, measured: Res<Measured>, mut stats: ResMut<DemoStats>) {
    /// The extrinsic arm's banner.
    const AMBER: Color = Color::srgb(1.0, 0.78, 0.30);
    /// The intrinsic arm's banner.
    const GREEN: Color = Color::srgb(0.35, 0.95, 0.55);

    let arm = if demo.after { "AFTER" } else { "BEFORE" };
    let Some(m) = &measured.0 else {
        stats.title = format!(
            "E-321  intrinsic flips - field {} at {} samples/axis",
            field_name(demo.field).unwrap_or("?"),
            LADDER[demo.rung]
        );
        stats.extra = vec!["measuring...".into()];
        return;
    };

    stats.title = format!(
        "E-321  intrinsic flips - {} at {} samples/axis - showing {arm}",
        m.field, m.samples
    );
    stats.vertices = m.positions.len();
    stats.triangles = if demo.after {
        m.after_indices.len() / 3
    } else {
        m.before_indices.len() / 3
    };
    stats.extract_ms = m.extract_ms;
    stats.banner = Some(if demo.after {
        (
            "AFTER  -  intrinsic Delaunay connectivity, same vertices".into(),
            GREEN,
        )
    } else {
        (
            "BEFORE -  marching cubes connectivity, same vertices".into(),
            AMBER,
        )
    });
    stats.hint = Some("[F] before/after   [H] HUD".into());
    stats.keys = Some(
        "[F] before/after   [V] vertices   [S] slivers   [1-7] field (ISOMESH_FIELD=7 for the 8th)\n\
         [ ] resolution   [W] wire   [N] normals   [G] box   [H] HUD   [F12] shot   [Esc] quit"
            .into(),
    );

    let census = &m.census;
    stats.extra = vec![
        format!(
            "{} samples/axis   [ and ] step {:?}, the resolutions golden_hashes.json holds",
            m.samples, LADDER
        ),
        String::new(),
        format!("{:<28}{:>10}{:>10}", "", "before", "after"),
        format!(
            "{:<28}{:>10.3}{:>10.3}",
            "worst-decile min angle", m.before.worst_decile_mean, m.after.worst_decile_mean
        ),
        format!(
            "{:<28}{:>10.3}{:>10.3}",
            "10th-percentile min angle", m.before.percentile_10, m.after.percentile_10
        ),
        format!(
            "{:<28}{:>10.3}{:>10.3}",
            "global min angle", m.before.global_min, m.after.global_min
        ),
        format!(
            "{:<28}{:>10}{:>10}",
            format!("slivers below {SLIVER_DEGREES:.0} deg"),
            m.before.slivers,
            m.after.slivers
        ),
        format!(
            "{:<28}{:>10}{:>10}",
            "non-Delaunay interior edges", m.non_delaunay_before, m.non_delaunay_after
        ),
        format!(
            "{:<28}{:>10}{:>10}",
            "chords that are not geodesic", m.gap_before.differing, m.gap_after.differing
        ),
        String::new(),
        format!(
            "flips {} of a {} budget ({FLIP_BUDGET_PER_EDGE} per interior edge x {}), rejected {}, exhausted {}, {:.3} ms",
            m.flips.flips,
            m.flips.budget,
            census.interior,
            m.flips.rejected,
            yes_no(m.flips.exhausted),
            m.flip_ms
        ),
        format!(
            "edges: interior {}  boundary {}  non-manifold {}  mis-oriented {}  zero-length {}",
            census.interior,
            census.boundary,
            census.non_manifold,
            census.inconsistently_oriented,
            census.zero_length
        ),
        format!(
            "after: {} self-edges, {} doubled vertex pairs (a Delta-complex, never refused); {} degenerate",
            m.delta_after.0, m.delta_after.1, m.after.degenerate
        ),
        String::new(),
        format!(
            "vertex_positions_moved {}   bit-for-bit against a copy taken before the flips",
            m.positions_moved
        ),
        format!(
            "hashes_moved           {}   mesh_hash live {:016x} vs committed {:016x}",
            m.hashes_moved, m.live_hash, m.committed_hash
        ),
        format!(
            "                           {GOLDEN_ALGORITHM}/{}/{} in golden_hashes.json, the file golden/tests.rs:59 gates",
            m.field, m.samples
        ),
        format!(
            "control: the same positions with the FLIPPED indices moves that hash: {}",
            yes_no(m.control_hash_moved)
        ),
        format!(
            "         so hashes_moved = 0 could have been 1 (M-44); extrinsic buffer untouched: {}",
            yes_no(m.extrinsic_identical)
        ),
        String::new(),
        "C2 HELD       no position moved, no golden hash moved. after the build this file".into(),
        "              reads no position, so a pinned vertex set is the mechanism (P-175)".into(),
        format!(
            "C1 FALSIFIED  the worst decile gained {:.3} deg here; the bar was {C1_DEGREES:.0} deg on",
            m.gain()
        ),
        format!(
            "              {C1_FIELDS} fields, and p-175.csv records c1_fields = {CITED_C1_FIELDS} of 8 (P-175)"
        ),
        format!(
            "C3 FALSIFIED  {CITED_C3_BENEFITING} of {CITED_C3_SURVEYED} connectivity consumers in crates/isomesh/src benefit:"
        ),
        "              every one reads EXTRINSIC triangles, so a retriangulation they cannot"
            .into(),
        "              see changes nothing downstream. the registered falsifier (P-175)".into(),
        String::new(),
        "the after arm is drawn as CHORDS between the same vertices. every flipped edge took"
            .into(),
        "its length from the unfolded quadrilateral, never from the two positions -- worst".into(),
        format!(
            "|chord - intrinsic length| = {:.6} cell on {} of {} edges. that is the whole ticket.",
            m.gap_after.worst_cells, m.gap_after.differing, m.gap_after.edges
        ),
    ];
}

#[cfg(test)]
mod tests {
    use core::time::Duration;

    use bevy::asset::AssetApp;
    use bevy::ecs::system::RunSystemOnce;
    use bevy::time::TimeUpdateStrategy;

    use super::*;

    /// A frame, fixed, so nothing here depends on how long the test machine took.
    const FRAME: Duration = Duration::from_millis(16);

    /// `p-175.csv`, read at compile time so every figure this demo reproduces or
    /// quotes is held against the artefact rather than against a memory of it.
    const P175_CSV: &str = include_str!("../../docs/experiments/p-175.csv");

    /// One CSV as a header row plus data rows, comment lines dropped.
    ///
    /// The experiment CSVs carry the hypothesis, the falsifier and the provenance
    /// as `#` lines above the header, which is why this cannot be a
    /// `lines().next()`.
    fn table(csv: &str) -> (Vec<&str>, Vec<Vec<&str>>) {
        let mut rows = csv.lines().filter(|l| !l.starts_with('#'));
        let header: Vec<&str> = rows
            .next()
            .expect("the CSV has a header row")
            .split(',')
            .collect();
        let data = rows
            .filter(|l| !l.trim().is_empty())
            .map(|l| l.split(',').collect())
            .collect();
        (header, data)
    }

    /// One cell of `p-175.csv`, by field, resolution and column name.
    fn cell(field: &str, samples: u32, column: &str) -> String {
        let (header, rows) = table(P175_CSV);
        let at = |name: &str| {
            header
                .iter()
                .position(|h| *h == name)
                .unwrap_or_else(|| panic!("{name} is not a column of p-175.csv"))
        };
        let wanted = at(column);
        let res = samples.to_string();
        for row in &rows {
            if row.get(at("field")) == Some(&field)
                && row.get(at("resolution")) == Some(&res.as_str())
            {
                return row
                    .get(wanted)
                    .unwrap_or_else(|| panic!("row is short of {column}"))
                    .to_string();
            }
        }
        panic!("p-175.csv has no {field} row at {samples}");
    }

    /// One integer cell.
    fn integer(field: &str, samples: u32, column: &str) -> u64 {
        cell(field, samples, column).parse().unwrap_or_else(|_| {
            panic!("p-175.csv {column} for {field}/{samples} is not an integer")
        })
    }

    /// One float cell.
    fn number(field: &str, samples: u32, column: &str) -> f64 {
        cell(field, samples, column)
            .parse()
            .unwrap_or_else(|_| panic!("p-175.csv {column} for {field}/{samples} is not a number"))
    }

    /// Whether a live `f64` is the number the CSV holds.
    ///
    /// Not `==`: the artefact renders to six decimals, so a bit-exact comparison
    /// would be a test of two formatters agreeing. Half a unit in the last place
    /// of that rendering is the exact rounding bound; a hair over it absorbs the
    /// decimal-to-binary round trip.
    fn reproduces(live: f64, cited: f64) -> bool {
        (live - cited).abs() <= 6e-7 + 1e-9 * cited.abs()
    }

    /// The demo's headless app: its own systems, no window and no renderer.
    ///
    /// `setup` is left out because it wants an `Assets<StandardMaterial>` and a
    /// `GizmoConfigStore`; the one thing it produces that `rebuild` needs is the
    /// surface material handle, inserted here by hand. `report` is left out too
    /// and run as a one-shot below, which is the same system the demo runs every
    /// frame.
    ///
    /// No stall-detecting drain: the rebuild is synchronous, so one stepped frame
    /// produces the whole measurement. The panel is written from a resource, not
    /// from geometry published through `Commands`, so there is nothing to wait on.
    fn harness(field: usize, rung: usize, after: bool) -> App {
        let mut app = App::new();
        app.add_plugins(MinimalPlugins)
            .add_plugins(bevy::asset::AssetPlugin::default())
            .init_asset::<Mesh>()
            .insert_resource(TimeUpdateStrategy::ManualDuration(FRAME))
            .init_resource::<ButtonInput<KeyCode>>()
            // `field` pinned, the rest left to `ViewFlags::default()`.
            // `Default` reads `ISOMESH_FIELD` from the environment and
            // `unsafe_code = "forbid"` means a test here cannot set that
            // variable (`common/mod.rs:122-126`), so it has to be out-ranked
            // rather than arranged. No system under test reads any other flag,
            // and `ViewFlags::parse` is private to the harness module.
            .insert_resource(ViewFlags {
                field,
                ..ViewFlags::default()
            })
            .insert_resource(Capture::default())
            .insert_resource(Demo {
                field,
                rung,
                after,
                vertices_shown: true,
                slivers_shown: true,
            })
            .init_resource::<Measured>()
            .init_resource::<Arms>()
            .init_resource::<DemoStats>()
            .insert_resource(SurfaceMaterial(Handle::default()))
            .add_systems(Update, (controls, rebuild, show).chain());
        app
    }

    /// One frame, with the input clearing `InputPlugin` would have done.
    ///
    /// Without it `just_pressed` stays true for ever, and `F` would toggle the
    /// arm on every frame of a test rather than on the one the key was pressed.
    fn step(app: &mut App) {
        app.update();
        app.world_mut()
            .resource_mut::<ButtonInput<KeyCode>>()
            .clear();
    }

    /// Run the demo once and hand back the panel a reader would read.
    fn panel(field: usize, rung: usize, after: bool) -> (String, Vec<String>, Option<String>) {
        let mut app = harness(field, rung, after);
        step(&mut app);
        app.world_mut()
            .run_system_once(report)
            .expect("the HUD system");
        let stats = app.world().resource::<DemoStats>();
        (
            stats.title.clone(),
            stats.extra.clone(),
            stats.banner.clone().map(|(line, _)| line),
        )
    }

    /// Find the one panel line containing `needle`.
    fn line<'a>(lines: &'a [String], needle: &str) -> &'a str {
        lines
            .iter()
            .find(|l| l.contains(needle))
            .map(String::as_str)
            .unwrap_or_else(|| panic!("no panel line contains {needle:?}"))
    }

    /// The field indices this demo shows are the fixture's own, in order.
    ///
    /// [`field_name`] and [`measure`] both count through
    /// `for_each_reference_field!`, and `golden_hashes.json` was generated by the
    /// same macro (`golden.rs:213`). This is what makes an index a valid key into
    /// the fixture, so a ninth reference field cannot silently shift the mapping.
    #[test]
    fn the_field_order_is_the_fixtures_own() {
        let mut names = Vec::new();
        for index in 0..FIELD_COUNT {
            names.push(field_name(index).expect("a name per field"));
        }
        assert!(
            field_name(FIELD_COUNT).is_none(),
            "FIELD_COUNT is {FIELD_COUNT} but the macro expands further"
        );

        // The fixture's own order: first appearance of each field name.
        let mut committed: Vec<&str> = Vec::new();
        for chunk in GOLDEN_HASHES.split('{').skip(1) {
            let at = chunk.find("\"field\":\"").expect("a field key") + "\"field\":\"".len();
            let rest = &chunk[at..];
            let end = rest.find('"').expect("a closed string");
            let name = &rest[..end];
            if !committed.contains(&name) {
                committed.push(name);
            }
        }
        assert_eq!(
            names, committed,
            "the demo's field order is not golden_hashes.json's"
        );
    }

    /// Every row this demo can show has a committed hash, and there are 24.
    ///
    /// The bench's own vacuity control: `golden_hashes.json` must yield exactly
    /// 24 `marching_cubes` rows, or `hashes_moved` is measured against a scanner
    /// that matched nothing.
    #[test]
    fn every_row_this_demo_can_show_has_a_committed_hash() {
        let mut found = 0;
        for index in 0..FIELD_COUNT {
            let name = field_name(index).expect("a name per field");
            for samples in LADDER {
                assert!(
                    golden_hash(name, samples).is_some(),
                    "golden_hashes.json has no {GOLDEN_ALGORITHM}/{name}/{samples}"
                );
                found += 1;
            }
        }
        assert_eq!(found, 24, "the fixture rows this demo can move are 8 x 3");

        // And the scanner is selective: the decider variant hashes the same
        // field at the same resolution and must not be picked up.
        let sphere = golden_hash("sphere", 17).expect("the sphere row");
        assert!(
            GOLDEN_HASHES.contains("\"marching_cubes+decider\",\"field\":\"sphere\""),
            "the fixture no longer holds a same-field row for another algorithm"
        );
        assert_ne!(
            sphere, 0,
            "a zero hash would mean the hex parse silently produced nothing"
        );
    }

    /// The two sweep-level figures the panel quotes are what `p-175.csv` holds.
    ///
    /// **This is the test that makes those numbers citations rather than
    /// remembered figures.** They are the two the demo cannot recompute: C1's
    /// field count is a property of the whole sweep, and C3's is a survey of the
    /// crate's source. Re-run the experiment and this goes red before the panel
    /// can quote a figure the artefact no longer holds.
    #[test]
    fn the_cited_sweep_figures_are_what_p175_committed() {
        let (_, rows) = table(P175_CSV);
        assert_eq!(rows.len(), 24, "p-175.csv is 8 fields x 3 resolutions");
        for index in 0..FIELD_COUNT {
            let name = field_name(index).expect("a name per field");
            for samples in LADDER {
                assert_eq!(
                    integer(name, samples, "c1_fields"),
                    u64::from(CITED_C1_FIELDS),
                    "p-175.csv c1_fields for {name}/{samples}"
                );
                assert_eq!(
                    integer(name, samples, "c3_consumers_benefiting"),
                    u64::from(CITED_C3_BENEFITING),
                    "p-175.csv c3_consumers_benefiting for {name}/{samples}"
                );
                assert_eq!(
                    integer(name, samples, "c3_consumers_surveyed"),
                    u64::from(CITED_C3_SURVEYED),
                    "p-175.csv c3_consumers_surveyed for {name}/{samples}"
                );
                assert_eq!(
                    number(name, samples, "sliver_threshold_degrees"),
                    SLIVER_DEGREES,
                    "p-175.csv sliver_threshold_degrees for {name}/{samples}"
                );
                assert_eq!(
                    cell(name, samples, "min_angle_statistic"),
                    "worst_decile_mean_degrees",
                    "p-175.csv min_angle_statistic for {name}/{samples}"
                );
            }
        }
    }

    /// The live measurement is `p-175.csv`'s row, on every row the demo shows.
    ///
    /// **This is what makes the panel a measurement rather than an illustration.**
    /// The flipping loop here is the bench's, written example-locally, and if the
    /// two ever disagree then one of them is wrong and the demo is the one on
    /// screen. Twenty columns per row, twenty-four rows.
    #[test]
    fn the_live_measurement_reproduces_p175s_row() {
        for index in 0..FIELD_COUNT {
            let name = field_name(index).expect("a name per field");
            for samples in LADDER {
                let m = measure(index, samples)
                    .unwrap_or_else(|| panic!("{name} at {samples} did not measure"));
                assert_eq!(m.field, name, "measure returned the wrong field");

                for (column, live) in [
                    ("vertices", m.positions.len() as u64),
                    ("triangles", (m.before_indices.len() / 3) as u64),
                    ("interior_edges", m.census.interior),
                    ("boundary_edges", m.census.boundary),
                    ("non_manifold_edges", m.census.non_manifold),
                    (
                        "inconsistently_oriented_edges",
                        m.census.inconsistently_oriented,
                    ),
                    ("zero_length_edges", m.census.zero_length),
                    ("slivers_before", m.before.slivers),
                    ("slivers_after", m.after.slivers),
                    ("degenerate_triangles", m.before.degenerate),
                    ("non_delaunay_before", m.non_delaunay_before),
                    ("non_delaunay_after", m.non_delaunay_after),
                    ("flips", m.flips.flips),
                    ("flips_rejected", m.flips.rejected),
                    ("flip_budget", m.flips.budget),
                    ("vertex_positions_moved", m.positions_moved),
                    ("hashes_moved", m.hashes_moved),
                ] {
                    assert_eq!(
                        integer(name, samples, column),
                        live,
                        "p-175.csv {column} for {name}/{samples}"
                    );
                }

                for (column, live) in [
                    ("min_angle_before", m.before.worst_decile_mean),
                    ("min_angle_after", m.after.worst_decile_mean),
                    ("p10_min_angle_before", m.before.percentile_10),
                    ("p10_min_angle_after", m.after.percentile_10),
                    ("global_min_angle_before", m.before.global_min),
                    ("global_min_angle_after", m.after.global_min),
                    ("min_angle_gain", m.gain()),
                ] {
                    let cited = number(name, samples, column);
                    assert!(
                        reproduces(live, cited),
                        "p-175.csv {column} for {name}/{samples} is {cited}, live {live}"
                    );
                }

                for (column, live) in [
                    ("extrinsic_geometry_identical", m.extrinsic_identical),
                    ("flip_budget_exhausted", m.flips.exhausted),
                    ("c1_holds", m.gain() >= C1_DEGREES),
                    (
                        "c2_holds",
                        m.positions_moved == 0 && m.hashes_moved == 0 && m.extrinsic_identical,
                    ),
                ] {
                    assert_eq!(
                        cell(name, samples, column),
                        live.to_string(),
                        "p-175.csv {column} for {name}/{samples}"
                    );
                }

                // The control column is only asserted where there was a flip to
                // see: with `flips = 0` the flipped connectivity *is* the
                // committed connectivity and the hash correctly does not move.
                let cited_control = cell(name, samples, "control_hash_moved") == "true";
                assert_eq!(
                    m.control_hash_moved, cited_control,
                    "p-175.csv control_hash_moved for {name}/{samples}"
                );
                if m.flips.flips > 0 {
                    assert!(
                        m.control_hash_moved,
                        "{name}/{samples} flipped {} edges and the hash did not move, so \
                         hashes_moved = 0 would prove nothing (M-44)",
                        m.flips.flips
                    );
                }
            }
        }
    }

    /// The new edge's length comes from the unfolded quadrilateral, not the chord.
    ///
    /// **This is the ticket.** An extrinsic remesh would leave every edge length
    /// equal to the chord between its endpoints, so the falsifier is available and
    /// cheap: before flipping the two agree on every edge to the bit, and after
    /// flipping they must disagree on some. Both halves are asserted, because the
    /// `after` count alone would also be satisfied by an instrument that reported
    /// a difference everywhere.
    #[test]
    fn the_new_edge_length_is_intrinsic_rather_than_the_chord() {
        let mut fields_with_flips = 0;
        for index in 0..FIELD_COUNT {
            let name = field_name(index).expect("a name per field");
            let m = measure(index, 25).unwrap_or_else(|| panic!("{name} at 25 did not measure"));
            assert_eq!(
                m.gap_before.differing, 0,
                "{name}: the unflipped arm's lengths came from those very chords, so none of \
                 its {} edges can differ",
                m.gap_before.edges
            );
            assert_eq!(
                m.gap_before.edges, m.gap_after.edges,
                "{name}: a flip preserves the edge count"
            );
            if m.flips.flips == 0 {
                assert_eq!(
                    m.gap_after.differing, 0,
                    "{name}: nothing was flipped, so nothing can have left its chord"
                );
                continue;
            }
            fields_with_flips += 1;
            assert!(
                m.gap_after.differing > 0,
                "{name}: {} flips and not one edge left its chord -- the length update read \
                 positions instead of the unfolded quadrilateral",
                m.flips.flips
            );
            assert!(
                m.gap_after.differing <= m.gap_after.edges,
                "{name}: more differing edges than edges"
            );
            assert!(
                m.gap_after.worst_cells > 0.0,
                "{name}: a differing edge with a zero worst gap is a broken instrument"
            );
        }
        assert!(
            fields_with_flips >= 4,
            "only {fields_with_flips} fields flipped anything, so this test is nearly vacuous"
        );
    }

    /// The position instrument can report movement.
    ///
    /// `vertex_positions_moved = 0` is the demo's headline, and a zero that could
    /// not have been non-zero is not a measurement (M-44). One ULP on one
    /// coordinate is the smallest movement there is, and it must be seen: the
    /// comparison is over bits, not values.
    #[test]
    fn the_position_instrument_can_report_movement() {
        let untouched = [[1.0f64, 2.0, 3.0], [-4.0, 5.0, 6.0]];
        assert_eq!(moved_positions(&untouched, &untouched), 0);

        let mut nudged = untouched;
        nudged[1][2] = f64::from_bits(nudged[1][2].to_bits() + 1);
        assert_eq!(
            moved_positions(&nudged, &untouched),
            1,
            "one ULP on one coordinate has to count as a moved position"
        );

        // Signed zero is the case `mesh_hash` is built around: the two compare
        // equal and hash differently, so the instrument must follow the bits.
        let plus = [[0.0f64, 0.0, 0.0]];
        let minus = [[-0.0f64, 0.0, 0.0]];
        assert_eq!(
            moved_positions(&minus, &plus),
            1,
            "-0.0 and +0.0 are the same value and different bits"
        );

        // A shorter buffer is movement too, or a truncation would read as none.
        assert_eq!(moved_positions(&untouched[..1], &untouched), 1);
    }

    /// The flipping loop reaches the intrinsic Delaunay fixed point.
    ///
    /// Both halves of the bench's control: the budget was not exhausted, and no
    /// flippable edge fails the Delaunay test afterwards. Without the second, the
    /// `after` arm is a partial run and C1's gain is not the intrinsic Delaunay
    /// gain; without the first, the loop might merely have been stopped.
    #[test]
    fn the_flip_loop_reaches_the_fixed_point() {
        let mut flipped_something = 0;
        for index in 0..FIELD_COUNT {
            let name = field_name(index).expect("a name per field");
            let m = measure(index, 33).unwrap_or_else(|| panic!("{name} at 33 did not measure"));
            assert!(
                !m.flips.exhausted,
                "{name}: the flip budget of {} ran out",
                m.flips.budget
            );
            assert_eq!(
                m.non_delaunay_after, 0,
                "{name}: {} flippable edges still fail the Delaunay test",
                m.non_delaunay_after
            );
            assert_eq!(
                m.flips.rejected, 0,
                "{name}: the convexity guard fired, which the theorem says it cannot"
            );
            assert_eq!(
                m.flips.budget,
                FLIP_BUDGET_PER_EDGE * m.census.interior.max(1),
                "{name}: the budget is not 64 per interior edge"
            );
            if m.flips.flips > 0 {
                flipped_something += 1;
                assert!(
                    m.non_delaunay_before > 0,
                    "{name}: flips happened with nothing failing the criterion"
                );
            }
        }
        assert!(
            flipped_something >= 4,
            "only {flipped_something} fields flipped anything at 33^3"
        );
    }

    /// The panel reports C2 as held, with a control that could have said otherwise.
    ///
    /// **This is the only way to see this demo's screen on a machine with no
    /// display**, and these are the lines the whole example exists to put in front
    /// of a reader: the two zeros, the hash they were compared against, and the
    /// control that makes them mean something. A test that only checked
    /// `Measurement`'s fields would pass with any of them missing from the panel.
    #[test]
    fn the_hud_reports_the_pinned_vertices_and_the_control_that_licenses_them() {
        let (title, lines, banner) = panel(0, 1, false);
        for l in &lines {
            println!("{l}");
        }
        assert!(
            title.contains("E-321") && title.contains("sphere") && title.contains("BEFORE"),
            "the title lost its ticket, field or arm: {title}"
        );
        assert_eq!(
            banner.as_deref(),
            Some("BEFORE -  marching cubes connectivity, same vertices"),
            "the banner must name the arm on screen"
        );

        let moved = line(&lines, "vertex_positions_moved");
        assert!(
            moved.contains("vertex_positions_moved 0"),
            "a position moved, which would mean this is not an intrinsic flip: {moved}"
        );
        let hashes = line(&lines, "hashes_moved");
        let committed = golden_hash("sphere", 25).expect("the sphere row at 25");
        assert!(
            hashes.contains("hashes_moved           0"),
            "a golden hash moved: {hashes}"
        );
        assert!(
            hashes.contains(&format!("{committed:016x}")),
            "the panel stopped naming the committed hash it compared against: {hashes}"
        );
        let fixture = line(&lines, "in golden_hashes.json");
        assert!(
            fixture.contains("marching_cubes/sphere/25"),
            "the panel does not say which fixture row it read: {fixture}"
        );

        let control = line(&lines, "control:");
        assert!(
            control.ends_with("yes"),
            "the control did not move the hash, so the two zeros above prove nothing: {control}"
        );
        assert!(
            line(&lines, "M-44").contains("could have been 1"),
            "the panel dropped the reason the control is there"
        );
    }

    /// The panel says C1 and C3 were falsified, and by how much.
    ///
    /// A demo that implied a win the measurement did not find would be worse than
    /// no demo, so the two negative verdicts are asserted as strings a reader
    /// reads — including the live gain, which has to be under the bar it names.
    #[test]
    fn the_hud_says_c1_and_c3_were_falsified() {
        let (_, lines, _) = panel(0, 1, true);
        for l in &lines {
            println!("{l}");
        }

        assert!(
            line(&lines, "C2 HELD").contains("no golden hash moved"),
            "C2's verdict line changed shape"
        );

        let c1 = line(&lines, "C1 FALSIFIED");
        assert!(
            c1.contains("the bar was 10 deg"),
            "C1's line stopped naming the registered bar: {c1}"
        );
        let gain: f64 = c1
            .split("gained ")
            .nth(1)
            .and_then(|rest| rest.split(' ').next())
            .and_then(|n| n.parse().ok())
            .unwrap_or_else(|| panic!("C1's line does not carry a live gain: {c1}"));
        assert!(
            gain < C1_DEGREES,
            "the panel says C1 was falsified and reports a gain of {gain} deg, which clears \
             the {C1_DEGREES} deg bar"
        );
        assert!(
            line(&lines, "c1_fields").contains(&format!("= {CITED_C1_FIELDS} of 8")),
            "C1's citation stopped naming p-175.csv's field count"
        );

        let c3 = line(&lines, "C3 FALSIFIED");
        assert!(
            c3.contains(&format!(
                "{CITED_C3_BENEFITING} of {CITED_C3_SURVEYED} connectivity consumers"
            )),
            "C3's line stopped naming the survey: {c3}"
        );
        assert!(
            line(&lines, "EXTRINSIC triangles").contains("changes nothing")
                || lines.iter().any(|l| l.contains("changes nothing")),
            "C3's line stopped saying why nothing changes downstream"
        );
        assert_eq!(
            lines.iter().filter(|l| l.contains("(P-175)")).count(),
            3,
            "each of the three verdicts cites P-175 exactly once"
        );
    }

    /// The chord caveat is on the panel, with a live number behind it.
    ///
    /// The after arm's edges are drawn as chords of geodesics, and a reader who
    /// is not told that will read the picture as an extrinsic remesh — the one
    /// misreading this demo must not permit.
    #[test]
    fn the_hud_admits_the_after_arm_is_drawn_as_chords() {
        let (_, lines, _) = panel(0, 1, true);
        let row = line(&lines, "chords that are not geodesic");
        let counts: Vec<u64> = row
            .split_whitespace()
            .filter_map(|w| w.parse().ok())
            .collect();
        assert_eq!(
            counts.len(),
            2,
            "the chord row must carry a before and an after count: {row}"
        );
        assert_eq!(counts[0], 0, "the unflipped arm's chords are its lengths");
        assert!(
            counts[1] > 0,
            "the flipped arm reports no edge off its chord, so the length update is not \
             intrinsic: {row}"
        );

        let caveat = line(&lines, "drawn as CHORDS");
        assert!(
            caveat.contains("between the same vertices"),
            "the caveat stopped saying the vertices are shared: {caveat}"
        );
        let worst = line(&lines, "|chord - intrinsic length|");
        assert!(
            worst.contains(" cell on ") && worst.contains(" edges."),
            "the caveat stopped quantifying the gap: {worst}"
        );
    }

    /// Exactly the arm carries what the harness wireframes.
    ///
    /// `common::draw_wireframe` selects `(&Mesh3d, &GlobalTransform)` filtered by
    /// `With<DemoMesh>` (`common/mod.rs:760-765`). A demo that forgets the marker
    /// draws no wireframe at all, and one that puts it on both entities draws the
    /// extrinsic arm on top of the intrinsic one — so the alternation would show
    /// nothing either way. M-241 records a committed GIF that advertised a sweep
    /// it never performed; a wireframe demo whose wireframe is empty is the same
    /// defect, and neither is visible from a terminal.
    #[test]
    fn exactly_the_arm_is_what_the_harness_wireframes() {
        let mut app = harness(0, 1, false);
        step(&mut app);

        let marked: Vec<(Option<Handle<Mesh>>, bool, Option<Panel>)> = app
            .world()
            .iter_entities()
            .filter(|e| e.contains::<DemoMesh>())
            .map(|e| {
                (
                    e.get::<Mesh3d>().map(|m| m.0.clone()),
                    e.contains::<GlobalTransform>(),
                    e.get::<Panel>().copied(),
                )
            })
            .collect();

        assert_eq!(
            marked.len(),
            1,
            "the harness wireframes {} entities, not the one arm",
            marked.len()
        );
        let (handle, has_transform, panel) = &marked[0];
        assert_eq!(
            *panel,
            Some(Panel::Wire),
            "the wireframe marker is on the wrong entity"
        );
        assert!(
            has_transform,
            "no GlobalTransform, so draw_wireframe's query skips it and the demo draws nothing"
        );
        assert_eq!(
            handle.as_ref(),
            Some(&app.world().resource::<Arms>().before),
            "the wireframed handle is not the arm on screen"
        );

        // And the shaded entity is deliberately *not* wireframed: it holds the
        // extrinsic mesh at all times, so a marker there would draw the before
        // arm underneath whichever arm the toggle selected.
        let shaded = app
            .world()
            .iter_entities()
            .filter(|e| e.get::<Panel>() == Some(&Panel::Surface))
            .count();
        assert_eq!(shaded, 1, "there must be exactly one shaded surface");
        assert!(
            app.world()
                .iter_entities()
                .filter(|e| e.get::<Panel>() == Some(&Panel::Surface))
                .all(|e| !e.contains::<DemoMesh>()),
            "the shaded surface is wireframed too, so both arms are drawn at once"
        );
    }

    /// `F` swaps the arm and the panel follows it; nothing else moves.
    ///
    /// The alternation is the demo, so the toggle has to change the mesh the
    /// harness wireframes **and** the numbers the panel attributes to it — while
    /// the shaded surface and the vertex buffer stay exactly as they were.
    #[test]
    fn the_toggle_swaps_the_wireframe_and_leaves_the_surface_alone() {
        let mut app = harness(0, 1, false);
        step(&mut app);

        let handles = |app: &App| {
            let mut surface = None;
            let mut wire = None;
            for (mesh, panel) in app
                .world()
                .iter_entities()
                .filter_map(|e| Some((e.get::<Mesh3d>()?, e.get::<Panel>()?)))
            {
                match panel {
                    Panel::Surface => surface = Some(mesh.0.clone()),
                    Panel::Wire => wire = Some(mesh.0.clone()),
                }
            }
            (
                surface.expect("a surface panel"),
                wire.expect("a wire panel"),
            )
        };

        let (surface_before, wire_before) = handles(&app);
        let arms_before = {
            let arms = app.world().resource::<Arms>();
            (arms.before.clone(), arms.after.clone())
        };
        assert_eq!(
            wire_before, arms_before.0,
            "the wire panel does not start on the extrinsic arm"
        );
        assert_eq!(
            surface_before, arms_before.0,
            "the shaded surface must be the extrinsic mesh"
        );
        let vertices_before = app
            .world()
            .resource::<Measured>()
            .0
            .as_ref()
            .expect("a measurement")
            .positions
            .clone();

        app.world_mut()
            .resource_mut::<ButtonInput<KeyCode>>()
            .press(KeyCode::KeyF);
        step(&mut app);

        assert!(
            app.world().resource::<Demo>().after,
            "F did not select the flipped arm"
        );
        let (surface_after, wire_after) = handles(&app);
        assert_eq!(
            surface_after, surface_before,
            "the shaded surface changed, and it is the one thing that must not"
        );
        assert_eq!(
            wire_after, arms_before.1,
            "the wireframe did not move to the flipped arm"
        );
        assert_ne!(
            wire_after, wire_before,
            "both arms resolved to one asset, so the toggle shows nothing"
        );
        assert_eq!(
            app.world()
                .resource::<Measured>()
                .0
                .as_ref()
                .expect("a measurement")
                .positions,
            vertices_before,
            "the vertex buffer moved, which is the one thing this demo claims cannot happen"
        );

        // And the panel now describes the arm on screen.
        app.world_mut()
            .run_system_once(report)
            .expect("the HUD system");
        let title = app.world().resource::<DemoStats>().title.clone();
        assert!(
            title.contains("AFTER"),
            "the panel still says BEFORE: {title}"
        );
    }

    /// Under capture the arm alternates with no keyboard, on a period that loops.
    ///
    /// The recorder presses nothing, so a subject that only changes on a keypress
    /// captures as a still — M-241 caught a committed GIF that advertised a sweep
    /// it never performed. Driven off `Capture::taken`, so the sequence is the
    /// captured frames rather than the wall clock.
    #[test]
    fn the_capture_sequence_alternates_and_loops() {
        let arm_at = |taken: u32| (taken / ALTERNATE_FRAMES) % 2 == 1;

        assert!(!arm_at(0), "the clip must open on the extrinsic arm");
        assert!(!arm_at(ALTERNATE_FRAMES - 1), "the first arm is held");
        assert!(arm_at(ALTERNATE_FRAMES), "the arm must swap after ten");
        assert!(arm_at(2 * ALTERNATE_FRAMES - 1), "the second arm is held");
        assert!(!arm_at(2 * ALTERNATE_FRAMES), "the period is twenty");

        // `record_gif.sh:47` defaults to 80 frames, so the clip is a whole number
        // of periods and its last frame is the state its first frame was in.
        const RECORDED: u32 = 80;
        assert_eq!(
            RECORDED % (2 * ALTERNATE_FRAMES),
            0,
            "an 80-frame clip is not a whole number of {}-frame periods, so the GIF jumps \
             where it loops",
            2 * ALTERNATE_FRAMES
        );

        // And the keyboard is refused: `controls` returns before reading a key.
        let mut app = harness(0, 1, false);
        app.insert_resource(Capture::default());
        step(&mut app);
        app.world_mut()
            .resource_mut::<ButtonInput<KeyCode>>()
            .press(KeyCode::KeyF);
        let capturing = app.world().resource::<Capture>().is_active();
        step(&mut app);
        if capturing {
            assert!(
                !app.world().resource::<Demo>().after,
                "a keypress moved the arm during a capture"
            );
        } else {
            assert!(
                app.world().resource::<Demo>().after,
                "F is ignored outside a capture too"
            );
        }
    }

    /// `[` and `]` walk the committed resolutions and stop at the ends.
    ///
    /// The ladder is exactly `golden.rs:73`'s, and stepping off it would leave the
    /// panel's headline number with no fixture to compare against.
    #[test]
    fn the_resolution_ladder_is_the_committed_one() {
        assert_eq!(LADDER, [17, 25, 33], "the ladder left golden.rs:73's");

        let mut app = harness(0, 0, false);
        step(&mut app);
        app.world_mut()
            .resource_mut::<ButtonInput<KeyCode>>()
            .press(KeyCode::BracketLeft);
        step(&mut app);
        assert_eq!(
            app.world().resource::<Demo>().rung,
            0,
            "[ stepped below the ladder"
        );

        for expected in [1usize, 2, 2] {
            {
                // Released before it is pressed again. `step`'s `clear()` is
                // what `InputPlugin` does -- it drops `just_pressed` and
                // *keeps* `pressed` -- so a second `press` of a key still held
                // registers no new edge, and the same key twice in a row would
                // silently move the rung once.
                let mut keys = app.world_mut().resource_mut::<ButtonInput<KeyCode>>();
                keys.release(KeyCode::BracketRight);
                keys.press(KeyCode::BracketRight);
            }
            step(&mut app);
            assert_eq!(
                app.world().resource::<Demo>().rung,
                expected,
                "] did not land on rung {expected}"
            );
        }
        let m = app.world().resource::<Measured>();
        assert_eq!(
            m.0.as_ref().expect("a measurement").samples,
            33,
            "the measurement did not follow the ladder"
        );
    }
}
