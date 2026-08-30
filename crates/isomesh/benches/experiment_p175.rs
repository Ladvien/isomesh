//! **P-175 — fixing sliver triangles without moving a single vertex.**
//!
//! Ticket: R-175. Pre-registered before this harness existed.
//!
//! ```bash
//! cargo bench --bench experiment_p175
//! ```
//!
//! Writes `docs/experiments/p-175.csv`.
//!
//! # What was missing
//!
//! This crate has measured its own slivers three times and has three times
//! declined to do anything about them, always for the same stated reason.
//! `M-48` found that when a grid **sample** lands on the isosurface every cut
//! edge meeting there places its own vertex at the same point and the
//! edge-keyed cache shares none of them. `M-185` found that completing the
//! crossing identity turns such a sliver into a repeated-index triangle.
//! `✗66`/`M-398` walked a capsule over one and found a sliver is hit at
//! **2.30×** its area share — real enrichment, six times too little to matter —
//! and closed with *"a Marching Cubes sliver is a thin triangle lying in the
//! same plane as its neighbours"* and *"do not add sliver collapsing or
//! quality-driven remeshing"*. `T-026`'s `mean_ratio` row states the standing
//! policy outright: both quality metrics are *"recorded, never gated … Marching
//! Cubes emits slivers wherever a grid corner sits near zero, and that is the
//! algorithm rather than a defect."*
//!
//! Every one of those is an argument about whether slivers are worth *fixing*.
//! None of them asked whether they **can** be fixed by connectivity alone, and
//! that is a different question with a hard constraint attached: this crate
//! commits 216 golden hashes over exact vertex bits (`M-31`'s fixture, now
//! `crates/isomesh/golden_hashes.json`, gated by `golden_hashes_are_unchanged`
//! at `golden/tests.rs:59`, and proven able to fire — `P-61` moved 135 of the
//! same 216). Any remesh that moves a position moves a hash. An **intrinsic**
//! triangulation cannot: it retriangulates the same piecewise-linear surface by
//! tracking edge *lengths* rather than positions, so the vertex set is
//! pointwise fixed by construction.
//!
//! `intrinsic Delaunay`, `signpost` and `Soliman` are zero in this repo, so the
//! mechanism is written here, bench-locally, and driven through the public API.
//!
//! # The mechanism, stated before the numbers
//!
//! **The flip criterion.** An interior edge `ab` shared by triangles `abc` and
//! `bad` is intrinsically Delaunay when the two opposite angles sum to at most
//! `pi`, equivalently when the cotangent weight `cot(alpha) + cot(beta) >= 0`.
//! Both cotangents come from lengths alone — `cot = (b^2 + c^2 - a^2) / (4A)`
//! with `A` from Heron in Kahan's stable ordering — so no angle is ever
//! evaluated to decide a flip and no position is ever read after the first
//! build. An edge is flipped when the sum falls below `-COTAN_TOLERANCE`.
//!
//! **The length update is intrinsic, and that distinction is the whole ticket.**
//! The two triangles are unfolded into the plane along their shared edge:
//! `a = (0, 0)`, `b = (L, 0)`, apex `c` above at
//! `x = (L^2 + |ca|^2 - |bc|^2) / 2L`, `y = 2A_abc / L`, apex `d` below by the
//! same construction with the sign of `y` reversed. The new edge's length is
//! `|cd|` **in that unfolded layout** — the length of a geodesic across the
//! surface, not the length of the chord between those two vertices in `R^3`.
//! Take the chord instead and the result is an ordinary extrinsic remesh of a
//! *different* surface, which is exactly the error this column set exists to
//! detect.
//!
//! **The flip is always well defined.** If `alpha + beta > pi` then `d` lies
//! strictly inside the circumcircle of `abc`, and for any such `d` the
//! tangent-chord angle at `a` bounds `angle(bad) < alpha`, so the quadrilateral
//! `acbd` has interior angle `angle(cab) + angle(bad) < angle(cab) + alpha =
//! pi - angle(abc) < pi` at `a`, and symmetrically at `b`. A strictly
//! non-Delaunay edge therefore always unfolds to a strictly convex
//! quadrilateral and the new diagonal always crosses the old one.
//! `flips_rejected` records how often the numerical guard on that theorem
//! fired; it is expected to be **0**.
//!
//! **The data structure is a Δ-complex, not a simplicial complex.** Corner slots
//! `3t+i` carry a vertex, a halfedge length and a twin slot, and after the build
//! nothing is ever looked up by vertex pair again. A flip that would produce a
//! self-edge, or a second edge between a pair that already has one, is therefore
//! legal and needs no special case — that is what the signpost literature means
//! by allowing a Δ-complex, and refusing those flips is what would make the
//! result something other than the intrinsic Delaunay triangulation.
//!
//! **The flip budget.** `FLIP_BUDGET_PER_EDGE = 64` flips per interior edge, so
//! the budget is `64 · interior_edges` and is recorded per row as `flip_budget`.
//! The loop is Lawson's: a work queue seeded with every interior edge, and on
//! each flip all six corner slots of the two modified triangles are re-queued —
//! which is exactly the set of edges whose cotangent sum can have changed, so an
//! empty queue is a proved fixed point rather than a hope.
//! `flip_budget_exhausted` and `non_delaunay_after` are both recorded, and
//! `non_delaunay_after` reading 0 is what licenses calling the "after" arm *the*
//! intrinsic Delaunay triangulation instead of a partial run.
//!
//! # Arms
//!
//! | arm | what it varies | is_control |
//! |---|---|---|
//! | `extrinsic` | the mesh exactly as `MarchingCubes<f64>::new()` emits it — the connectivity every consumer in the crate actually reads | no, it is the baseline |
//! | `intrinsic_delaunay` | the same vertex set and the same surface, flipped to the intrinsic Delaunay fixed point | no |
//! | `extrinsic_writeback` | the flipped connectivity written back into a `MeshBuffer` beside the untouched positions, and hashed | **yes** — C2's positive control, column `control_hash_moved` |
//! | `box_exact`, `thin_plate` | two fields whose Marching Cubes output is *already* intrinsically Delaunay | **yes** — the zero-flip control, `flips = 0` beside `non_delaunay_before = 0` |
//!
//! Eight reference fields at **17, 25 and 33** samples per axis: 24 rows. Those
//! three resolutions are not the house default of `33³`/`65³` — they are
//! `golden.rs:73`'s `RESOLUTIONS`, and they have to be, because `hashes_moved`
//! is a comparison against the committed fixture and the fixture exists at no
//! other resolution. The extractor is plain `marching_cubes` for the same
//! reason: it is row 2 of `Algorithm::ALL`, and its 24 rows are the only ones of
//! the 216 this harness runs. Claiming the other 192 would be claiming a
//! measurement nobody took.
//!
//! # SHARE, recomputed before the numbers
//!
//! The registration's SHARE reads *"C1 moves mesh quality, not cost; C3 decides
//! whether that quality reaches anything."* Discharged: **no clause here is a
//! cost claim and no clause reads a wall clock.** C1 is a difference of two
//! angles in degrees, C2 is three integer equalities, C3 is a count over a
//! survey of the crate's own source. `flip_ms` is recorded because the
//! registration names it and is read by no clause — `P-126`'s rule for
//! `wall_seconds` exactly — so `M-280`/`✗24`'s 1.45× governor scatter has
//! nothing here to bite on. It is still taken properly: one untimed warm-up,
//! then **7** repeats, the median as the headline with `flip_ms_min` and
//! `flip_ms_max` beside it, each repeat rebuilding the triangulation outside the
//! timed region so `flip_ms` is the flipping loop rather than the build. Every
//! repeat's flip count is asserted equal to the first run's, which makes the
//! repeats a determinism check as well as a timing one.
//!
//! # Vacuity controls
//!
//! Six, all before the first `run.record`:
//!
//! * **The registered one.** The worst-decile angle *before* flipping must be
//!   below 15° on at least one field, or there are no slivers to fix. Column:
//!   `min_angle_before`.
//! * **Something was flippable.** `non_delaunay_before` summed over the sweep
//!   must exceed zero, or the Delaunay criterion is being asked of a
//!   triangulation that already satisfies it everywhere and `min_angle_after`
//!   cannot differ from `min_angle_before` for any reason C1 is about.
//! * **Flips happened.** `flips` summed over the sweep must exceed zero.
//! * **The hash instrument can report movement.** On every row with `flips > 0`,
//!   writing the flipped connectivity back into a `MeshBuffer` beside the
//!   untouched positions **must** move the hash. Without this, `hashes_moved = 0`
//!   is a zero that could not have been non-zero (`M-44`) — it would be equally
//!   consistent with a `mesh_hash` blind to connectivity, or with a flipper that
//!   silently did nothing. Column: `control_hash_moved`.
//! * **The baseline is the committed fixture.** `golden_hashes.json` must yield
//!   exactly 24 `marching_cubes` rows, or `hashes_moved` is measured against a
//!   scanner that matched nothing.
//! * **The fixed point was reached.** `flip_budget_exhausted` false and
//!   `non_delaunay_after == 0` on every row, or the "after" arm is a partial run
//!   and C1's gain is not the intrinsic Delaunay gain.
//!
//! # What `min_angle_before` and `min_angle_after` are
//!
//! The **mean of the worst decile** of per-triangle minimum angles, in degrees,
//! over the non-degenerate triangles — the statistic the clause names, carried
//! in the CSV's own `min_angle_statistic` column so a reader of `p-175.csv`
//! never has to come back here to find out. The 10th-percentile value and the
//! global minimum are recorded beside it as `p10_min_angle_*` and
//! `global_min_angle_*`, so C1's verdict can be re-read under all three
//! statistics rather than resting on one choice of instrument. A triangle whose
//! Heron area is not positive has no angles and is excluded from all three;
//! `degenerate_triangles` counts them.
//!
//! A sliver is a triangle whose smallest angle is below `SLIVER_DEGREES = 15°`,
//! the same bar the registered vacuity control sets for the worst decile, and it
//! is carried per row as `sliver_threshold_degrees`.

#![allow(
    clippy::cast_possible_truncation,
    clippy::cast_precision_loss,
    clippy::too_many_lines
)]

mod common;

use std::collections::{BTreeSet, VecDeque};
use std::time::Instant;

use isomesh::fields::ReferenceField;
use isomesh::marching_cubes::MarchingCubes;
use isomesh::validate::mesh_hash;
use isomesh::{MeshBuffer, Sdf};

// ════════════════════════════════════════════════════════════════════════════
// constants
// ════════════════════════════════════════════════════════════════════════════

/// The three resolutions `golden.rs:73` commits hashes at.
const GOLDEN_RESOLUTIONS: [u32; 3] = [17, 25, 33];

/// The fixture rows this harness can move: 8 fields × 3 resolutions.
const GOLDEN_ROWS: usize = 24;

/// The name the fixture gives `MarchingCubes::<f64>::new()`.
const GOLDEN_ALGORITHM: &str = "marching_cubes";

/// A triangle whose smallest angle is below this many degrees is a sliver.
const SLIVER_DEGREES: f64 = 15.0;

/// An edge this short relative to the cell size is a near-duplicate vertex pair.
const SHORT_EDGE_FRACTION: f64 = 0.05;

/// A cotangent sum below `-COTAN_TOLERANCE` fails the intrinsic Delaunay test.
///
/// The sum is dimensionless — `(b² + c² − a²) / 4A` is a length squared over a
/// length squared — so a plain absolute tolerance is scale-free here.
const COTAN_TOLERANCE: f64 = 1e-12;

/// Flips allowed per interior edge before the loop is declared non-terminating.
const FLIP_BUDGET_PER_EDGE: u64 = 64;

/// C1's bar: degrees of gain on the worst decile.
const C1_DEGREES: f64 = 10.0;

/// C1's other bar: how many fields must clear `C1_DEGREES`.
const C1_FIELDS: usize = 4;

/// Timed repeats of the flipping loop, after one untimed warm-up.
const TIMED_REPEATS: usize = 7;

/// The corner slot that means "this halfedge has no twin".
const NO_TWIN: u32 = u32::MAX;

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
/// therefore independent of the order the triangles arrived in (T-004).
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
struct EdgeKey {
    lo: u32,
    hi: u32,
    slot: u32,
    forward: bool,
}

/// How the edges of a triangulation were classified when the twins were paired.
///
/// Every one of these is a *precondition* of the flip rather than a defect
/// report: an intrinsic flip is defined only across an interior edge of an
/// oriented manifold, and an edge that is not one is simply not flippable.
#[derive(Clone, Copy, Default)]
struct EdgeCensus {
    interior: u64,
    boundary: u64,
    non_manifold: u64,
    inconsistently_oriented: u64,
    zero_length: u64,
}

/// Angles of one triangulation, in degrees.
#[derive(Clone, Copy, Default)]
struct AngleStats {
    worst_decile_mean: f64,
    percentile_10: f64,
    global_min: f64,
    slivers: u64,
    degenerate: u64,
    short_edge_slivers: u64,
    shortest_edge: f64,
}

/// What one flipping run did.
#[derive(Clone, Copy)]
struct FlipReport {
    flips: u64,
    rejected: u64,
    budget: u64,
    exhausted: bool,
}

/// A triangulation of a fixed vertex set, carried as lengths rather than points.
///
/// Slot `3t+i` is corner `i` of triangle `t`. `corner[3t+i]` is the vertex there,
/// `length[3t+i]` is the length of the halfedge running from corner `i` to corner
/// `i+1`, and `twin[3t+i]` is the slot on the other side of that halfedge. After
/// `build` returns, no method here reads a position: that is what makes the
/// retriangulation intrinsic, and what makes `vertex_positions_moved` a fact
/// about the mechanism rather than an assertion about it.
struct Intrinsic {
    corner: Vec<u32>,
    length: Vec<f64>,
    twin: Vec<u32>,
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

    fn triangles(&self) -> u64 {
        (self.corner.len() / 3) as u64
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
    /// The count a cotangent Laplacian would carry as negative off-diagonal
    /// weights, which is C3's evidence — and, read *after* the flipping loop,
    /// the proof that the loop reached the intrinsic Delaunay fixed point.
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
    fn flip(&mut self, h: u32) -> bool {
        let o = self.twin[h as usize];
        let (hn, hp) = (next_slot(h), prev_slot(h));
        let (on, op) = (next_slot(o), prev_slot(o));

        // `h` runs a → b in triangle (a, b, c); `o` runs b → a in (b, a, d).
        let shared = self.length[h as usize];
        let l_bc = self.length[hn as usize];
        let l_ca = self.length[hp as usize];
        let l_ad = self.length[on as usize];
        let l_db = self.length[op as usize];

        // Unfold into the plane along `ab`, `c` above and `d` below.
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

        // Triangle `h/3` becomes (a, d, c): slot `h` is a → d, `hn` is the new
        // diagonal d → c, and `hp` keeps c → a untouched.
        self.length[h as usize] = l_ad;
        self.twin[h as usize] = tw_on;
        if tw_on != NO_TWIN {
            self.twin[tw_on as usize] = h;
        }
        self.corner[hn as usize] = corner_d;
        self.length[hn as usize] = l_new;
        self.twin[hn as usize] = on;

        // Triangle `o/3` becomes (b, c, d): slot `o` is b → c, `on` is the new
        // diagonal c → d, and `op` keeps d → b untouched.
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
    fn delaunay_flip(&mut self) -> FlipReport {
        let budget = FLIP_BUDGET_PER_EDGE * self.census.interior.max(1);
        let mut queued = vec![false; self.corner.len()];
        let mut queue: VecDeque<u32> = VecDeque::with_capacity(self.corner.len());
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
            // the six corners of the two modified triangles, which is exactly
            // the set of edges whose cotangent sum can have changed.
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

    /// Per-triangle minimum angles, in degrees, and what they say about slivers.
    fn angles(&self, short_edge: f64) -> AngleStats {
        let mut mins: Vec<f64> = Vec::with_capacity(self.corner.len() / 3);
        let mut degenerate = 0u64;
        let mut short_edge_slivers = 0u64;
        let mut shortest = f64::INFINITY;
        for sides in self.length.as_chunks::<3>().0 {
            let (l0, l1, l2) = (sides[0], sides[1], sides[2]);
            let least = l0.min(l1).min(l2);
            shortest = shortest.min(least);
            if area_of(l0, l1, l2) <= 0.0 {
                degenerate += 1;
                continue;
            }
            // The angle opposite `opp`, between the sides `x` and `y`. Clamped
            // because a right angle can round its cosine a hair outside [−1, 1].
            let angle = |opp: f64, x: f64, y: f64| {
                ((x * x + y * y - opp * opp) / (2.0 * x * y))
                    .clamp(-1.0, 1.0)
                    .acos()
                    .to_degrees()
            };
            let smallest = angle(l0, l1, l2)
                .min(angle(l1, l2, l0))
                .min(angle(l2, l0, l1));
            if smallest < SLIVER_DEGREES && least < short_edge {
                short_edge_slivers += 1;
            }
            mins.push(smallest);
        }
        if mins.is_empty() {
            return AngleStats {
                degenerate,
                shortest_edge: shortest,
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
            short_edge_slivers,
            shortest_edge: shortest,
        }
    }

    /// The flipped connectivity, as a flat index buffer.
    fn indices(&self) -> Vec<u32> {
        self.corner.clone()
    }
}

// ════════════════════════════════════════════════════════════════════════════
// C3 — the consumer survey
// ════════════════════════════════════════════════════════════════════════════

/// One place in `crates/isomesh/src/**` that reads triangle connectivity.
struct Consumer {
    symbol: &'static str,
    site: &'static str,
    /// Whether the consumer's answer depends on triangle *angles* at all.
    reads_angles: bool,
    /// Whether handing it the intrinsic triangulation would improve its answer.
    benefits: bool,
    why: &'static str,
}

/// Every consumer of connectivity the crate has, and what an intrinsic
/// retriangulation of the same vertex set would do for each.
///
/// This list is C3's answer, and it is auditable rather than asserted: the
/// harness prints it in full beside the numbers. `grep -E
/// 'cotangent|cot_|laplacian|geodesic'` over `crates/isomesh/src` returns four
/// doc comments — `dual.rs:229`, `dual.rs:550`, `dual_contouring.rs:57`,
/// `surface_nets.rs:154` — all describing the *uniform* smoothing pass, plus one
/// registration string, and **no implementation**. There is no cotangent
/// Laplacian, no geodesic solver and no parameterisation in this crate, and
/// those three are the operators an intrinsic triangulation exists to serve.
const CONSUMERS: [Consumer; 8] = [
    Consumer {
        symbol: "mass::mass_properties",
        site: "mass.rs:198",
        reads_angles: false,
        benefits: false,
        why: "divergence theorem over triangles, invariant under retriangulation of a fixed \
              polyhedron; and an intrinsic edge is a geodesic rather than a chord, so feeding it \
              this connectivity would integrate over a solid that does not exist",
    },
    Consumer {
        symbol: "normals::recompute/AreaWeightedFaces",
        site: "normals.rs:77 and :114",
        reads_angles: false,
        benefits: false,
        why: "weights face normals by area rather than by angle; the angle-weighted variant that \
              would care is not implemented, and the faces it needs are the extrinsic ones",
    },
    Consumer {
        symbol: "validate::MeshReport::mean_ratio",
        site: "validate.rs:257 and :903-943",
        reads_angles: true,
        benefits: false,
        why: "the crate's only triangle shape-quality number, and the one thing here that reads \
              exactly what a flip improves -- but it is computed from the extrinsic cross product \
              and edge lengths at validate.rs:917-939, which no intrinsic flip can move",
    },
    Consumer {
        symbol: "collider::readiness",
        site: "collider.rs:183 and :168",
        reads_angles: false,
        benefits: false,
        why: "counts boundary, non-manifold and duplicate features; a flip changes none of them, \
              and physics reads the extrinsic triangles regardless",
    },
    Consumer {
        symbol: "surface_nets::set_smoothing_passes",
        site: "surface_nets.rs:154 via dual.rs:550",
        reads_angles: false,
        benefits: false,
        why: "a uniform Laplacian over the face-adjacent CELL graph, not a cotangent Laplacian \
              over the mesh; it has no angle-derived weight to improve, and it runs before any \
              triangle exists",
    },
    Consumer {
        symbol: "orient::orient",
        site: "orient.rs:148",
        reads_angles: false,
        benefits: false,
        why: "propagates winding across shared edges; purely combinatorial and already exact",
    },
    Consumer {
        symbol: "validate::self_intersections",
        site: "validate/self_intersection.rs:153",
        reads_angles: false,
        benefits: false,
        why: "tests extrinsic triangle pairs for intersection; an intrinsic triangle is not a \
              subset of R^3 and cannot be handed to it",
    },
    Consumer {
        symbol: "validate::pinch_census",
        site: "validate/pinch.rs:417",
        reads_angles: false,
        benefits: false,
        why: "coincident-vertex grouping over positions; a flip moves no position, so every count \
              it reports is identical before and after",
    },
];

// ════════════════════════════════════════════════════════════════════════════
// the committed fixture
// ════════════════════════════════════════════════════════════════════════════

/// One committed golden row.
struct Golden {
    field: String,
    samples: u32,
    hash: u64,
}

/// The value of `"key"` in one `{…}` chunk of `golden_hashes.json`.
///
/// A hand-rolled scanner rather than a JSON parser, for the reason `golden.rs`
/// gives: the grammar is one line, fixed key order, no nesting and no escapes.
/// The crate's own reader is `#[cfg(test)]` and out of a bench's reach.
fn json_field(chunk: &str, key: &str) -> String {
    let needle = format!("\"{key}\":");
    let at = chunk
        .find(&needle)
        .unwrap_or_else(|| panic!("golden_hashes.json entry has no {key}"))
        + needle.len();
    let rest = &chunk[at..];
    if let Some(stripped) = rest.strip_prefix('"') {
        let end = stripped.find('"').expect("a closed string");
        stripped[..end].to_string()
    } else {
        let end = rest
            .find(|c: char| !c.is_ascii_digit())
            .unwrap_or(rest.len());
        rest[..end].to_string()
    }
}

/// The 24 committed `marching_cubes` rows of `M-31`'s 216.
///
/// Read from `crates/isomesh/golden_hashes.json`, the file
/// `golden_hashes_are_unchanged` (`golden/tests.rs:59`) gates against, so
/// `hashes_moved` is movement in **the** fixture rather than in a re-derivation
/// of it.
fn golden_marching_cubes() -> Vec<Golden> {
    let path = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("golden_hashes.json");
    let text = std::fs::read_to_string(&path)
        .unwrap_or_else(|e| panic!("golden_hashes.json is C2's baseline: {e}"));
    let mut out = Vec::new();
    for chunk in text.split('{').skip(1) {
        if json_field(chunk, "algorithm") != GOLDEN_ALGORITHM {
            continue;
        }
        out.push(Golden {
            field: json_field(chunk, "field"),
            samples: json_field(chunk, "samples").parse().expect("a resolution"),
            hash: u64::from_str_radix(&json_field(chunk, "hash"), 16).expect("a hex hash"),
        });
    }
    out
}

// ════════════════════════════════════════════════════════════════════════════
// one row
// ════════════════════════════════════════════════════════════════════════════

/// Everything one `(field, resolution)` produced.
struct Row {
    field: &'static str,
    resolution: u32,
    vertices: u64,
    triangles: u64,
    census: EdgeCensus,
    before: AngleStats,
    after: AngleStats,
    non_delaunay_before: u64,
    non_delaunay_after: u64,
    report: FlipReport,
    cell_size: f64,
    flip_ms: f64,
    flip_ms_min: f64,
    flip_ms_max: f64,
    build_ms: f64,
    positions_moved: u64,
    hashes_moved: u64,
    extrinsic_identical: bool,
    control_hash_moved: bool,
}

impl Row {
    /// C1's arithmetic for this row: degrees gained on the worst decile.
    fn gain(&self) -> f64 {
        self.after.worst_decile_mean - self.before.worst_decile_mean
    }

    fn c1(&self) -> bool {
        self.gain() >= C1_DEGREES
    }

    fn c2(&self) -> bool {
        self.positions_moved == 0 && self.hashes_moved == 0 && self.extrinsic_identical
    }
}

fn measure<F>(name: &'static str, field: &F, samples: u32, committed: &[Golden]) -> Row
where
    F: ReferenceField + Sdf<Scalar = f64>,
{
    let (shape, origin, cell_size) = common::grid::<f64, _>(field, samples);
    let mut mesh = MeshBuffer::<f64>::new();
    MarchingCubes::<f64>::new()
        .extract(field, &shape, origin, cell_size, &mut mesh)
        .expect("marching cubes over a reference field");

    let golden = committed
        .iter()
        .find(|row| row.field == name && row.samples == samples)
        .map(|row| row.hash)
        .unwrap_or_else(|| panic!("golden_hashes.json has no {GOLDEN_ALGORITHM}/{name}/{samples}"));

    // Independent copies, taken before the flipping run, so `vertex_positions_moved`
    // compares two buffers rather than one buffer with itself.
    let positions_before = mesh.positions.clone();
    let normals_before = mesh.normals.clone();
    let indices_before = mesh.indices.clone();

    let short_edge = SHORT_EDGE_FRACTION * cell_size;

    let base = Intrinsic::build(&mesh.positions, &mesh.indices);
    let before = base.angles(short_edge);
    let non_delaunay_before = base.negative_cotan_edges();

    let mut flipped = Intrinsic::build(&mesh.positions, &mesh.indices);
    let report = flipped.delaunay_flip();
    let after = flipped.angles(short_edge);
    let non_delaunay_after = flipped.negative_cotan_edges();

    // One untimed warm-up, then the repeats. The build sits outside the timed
    // region, so `flip_ms` is the flipping loop and `build_ms` is the rest.
    let mut warm = Intrinsic::build(&mesh.positions, &mesh.indices);
    warm.delaunay_flip();

    let mut builds = Vec::with_capacity(TIMED_REPEATS);
    let mut timings = Vec::with_capacity(TIMED_REPEATS);
    for _ in 0..TIMED_REPEATS {
        let started = Instant::now();
        let mut scratch = Intrinsic::build(&mesh.positions, &mesh.indices);
        builds.push(started.elapsed().as_secs_f64() * 1e3);
        let started = Instant::now();
        let repeat = scratch.delaunay_flip();
        timings.push(started.elapsed().as_secs_f64() * 1e3);
        assert_eq!(
            repeat.flips, report.flips,
            "{name}/{samples}: the flip loop is not deterministic, so no repeat of it measures \
             the same work"
        );
    }
    builds.sort_by(f64::total_cmp);
    timings.sort_by(f64::total_cmp);

    // C2's positive control: the flipped connectivity written back into a
    // `MeshBuffer` beside the untouched positions. If this does not move the
    // hash, neither `mesh_hash` nor the flipper can be believed.
    let mut control = MeshBuffer::<f64>::new();
    control.positions.clone_from(&mesh.positions);
    control.normals.clone_from(&mesh.normals);
    control.indices = flipped.indices();
    let control_hash_moved = mesh_hash(&control) != golden;

    let positions_moved = mesh
        .positions
        .iter()
        .zip(&positions_before)
        .filter(|(now, was)| {
            now.iter()
                .zip(was.iter())
                .any(|(x, y)| x.to_bits() != y.to_bits())
        })
        .count() as u64
        + (mesh.positions.len() as u64).abs_diff(positions_before.len() as u64);
    let normals_identical = mesh.normals.len() == normals_before.len()
        && mesh.normals.iter().zip(&normals_before).all(|(now, was)| {
            now.iter()
                .zip(was.iter())
                .all(|(x, y)| x.to_bits() == y.to_bits())
        });
    let extrinsic_identical =
        positions_moved == 0 && normals_identical && mesh.indices == indices_before;

    Row {
        field: name,
        resolution: samples,
        vertices: mesh.vertex_count() as u64,
        triangles: base.triangles(),
        census: base.census,
        before,
        after,
        non_delaunay_before,
        non_delaunay_after,
        report,
        cell_size,
        flip_ms: timings[timings.len() / 2],
        flip_ms_min: timings[0],
        flip_ms_max: timings[timings.len() - 1],
        build_ms: builds[builds.len() / 2],
        positions_moved,
        hashes_moved: u64::from(mesh_hash(&mesh) != golden),
        extrinsic_identical,
        control_hash_moved,
    }
}

// ════════════════════════════════════════════════════════════════════════════
// main
// ════════════════════════════════════════════════════════════════════════════

fn main() {
    if !std::env::args().any(|arg| arg == "--bench") {
        return;
    }
    let prereg = isomesh::experiment!("P-175");

    let committed = golden_marching_cubes();
    let mut rows: Vec<Row> = Vec::with_capacity(GOLDEN_ROWS);
    isomesh::for_each_reference_field!(f64, |name, field| {
        for samples in GOLDEN_RESOLUTIONS {
            rows.push(measure(name, &field, samples, &committed));
        }
    });

    // ── C3, from the survey rather than from an opinion ──────────────────────
    let surveyed = CONSUMERS.len() as u64;
    let benefiting = CONSUMERS.iter().filter(|c| c.benefits).count() as u64;
    let c3 = benefiting > 0;

    // ── C1's run-level verdict: fields clearing the bar at every resolution ──
    let fields: BTreeSet<&'static str> = rows.iter().map(|row| row.field).collect();
    let c1_fields = fields
        .iter()
        .filter(|name| rows.iter().filter(|row| row.field == **name).all(Row::c1))
        .count();
    let c1_run = c1_fields >= C1_FIELDS;

    common::experiment::run(prereg, |run| {
        // ── vacuity controls ─────────────────────────────────────────────────
        assert_eq!(
            committed.len(),
            GOLDEN_ROWS,
            "VOID: golden_hashes.json yielded {} `{GOLDEN_ALGORITHM}` rows rather than \
             {GOLDEN_ROWS}, so the scanner matched the wrong thing and `hashes_moved` would be \
             measured against a baseline that is not the committed fixture",
            committed.len()
        );
        assert_eq!(
            rows.len(),
            GOLDEN_ROWS,
            "VOID: the sweep produced {} rows rather than {GOLDEN_ROWS}, so the fixture and the \
             measurement do not cover the same population",
            rows.len()
        );

        let slivered = rows
            .iter()
            .filter(|row| row.before.worst_decile_mean < SLIVER_DEGREES)
            .count();
        assert!(
            slivered > 0,
            "VOID: the worst-decile angle before flipping is at or above {SLIVER_DEGREES} deg on \
             every one of the {GOLDEN_ROWS} rows, so there are no slivers to fix and C1's gain \
             would be measured against a triangulation with nothing wrong with it -- this is the \
             registered vacuity control"
        );

        let non_delaunay: u64 = rows.iter().map(|row| row.non_delaunay_before).sum();
        assert!(
            non_delaunay > 0,
            "VOID: not one edge in the whole sweep fails the intrinsic Delaunay criterion, so the \
             flipper is asked to change nothing and `min_angle_after` cannot differ from \
             `min_angle_before` for any reason C1 is about"
        );

        let total_flips: u64 = rows.iter().map(|row| row.report.flips).sum();
        assert!(
            total_flips > 0,
            "VOID: zero flips across the sweep against {non_delaunay} non-Delaunay edges, so the \
             `after` arm is the `before` arm and every clause here is unmeasured"
        );

        for row in rows.iter().filter(|row| row.report.flips > 0) {
            assert!(
                row.control_hash_moved,
                "VOID: {}/{} flipped {} edges, and writing that connectivity back beside the same \
                 positions did NOT move the golden hash -- so `hashes_moved = 0` is a zero that \
                 could not have been non-zero (M-44), and is equally consistent with a mesh_hash \
                 blind to connectivity or a flipper that changed nothing",
                row.field, row.resolution, row.report.flips
            );
        }

        for row in &rows {
            assert!(
                !row.report.exhausted && row.non_delaunay_after == 0,
                "VOID: {}/{} stopped with {} non-Delaunay edges still standing (budget {}, \
                 exhausted {}), so the `after` arm is a partial run and its angles are not the \
                 intrinsic Delaunay angles C1 asks about",
                row.field,
                row.resolution,
                row.non_delaunay_after,
                row.report.budget,
                row.report.exhausted
            );
        }

        // ── C3's survey, printed in full beside the numbers ──────────────────
        println!("\nP-175 C3 -- every consumer of triangle connectivity in crates/isomesh/src:");
        for consumer in &CONSUMERS {
            println!(
                "  {:<38} {:<34} angles={:<5} benefits={:<5} {}",
                consumer.symbol,
                consumer.site,
                consumer.reads_angles,
                consumer.benefits,
                consumer.why
            );
        }
        println!("  {benefiting} of {surveyed} benefit, so c3_holds = {c3}\n");

        // ── rows ─────────────────────────────────────────────────────────────
        for row in &rows {
            let c1 = row.c1();
            let c2 = row.c2();
            run.record(&[
                ("field", row.field.to_string()),
                ("resolution", row.resolution.to_string()),
                ("slivers_before", row.before.slivers.to_string()),
                ("slivers_after", row.after.slivers.to_string()),
                (
                    "min_angle_before",
                    format!("{:.6}", row.before.worst_decile_mean),
                ),
                (
                    "min_angle_after",
                    format!("{:.6}", row.after.worst_decile_mean),
                ),
                ("flips", row.report.flips.to_string()),
                ("vertex_positions_moved", row.positions_moved.to_string()),
                ("hashes_moved", row.hashes_moved.to_string()),
                (
                    "extrinsic_geometry_identical",
                    row.extrinsic_identical.to_string(),
                ),
                ("flip_ms", format!("{:.6}", row.flip_ms)),
                ("c1_holds", c1.to_string()),
                ("c2_holds", c2.to_string()),
                ("c3_holds", c3.to_string()),
                // ── extras (M-273) ──
                ("boundary_edges", row.census.boundary.to_string()),
                ("build_ms", format!("{:.6}", row.build_ms)),
                ("c1_fields", c1_fields.to_string()),
                ("c1_run", c1_run.to_string()),
                ("c3_consumers_benefiting", benefiting.to_string()),
                ("c3_consumers_surveyed", surveyed.to_string()),
                ("control_hash_moved", row.control_hash_moved.to_string()),
                ("degenerate_triangles", row.before.degenerate.to_string()),
                ("flip_budget", row.report.budget.to_string()),
                ("flip_budget_exhausted", row.report.exhausted.to_string()),
                ("flip_ms_max", format!("{:.6}", row.flip_ms_max)),
                ("flip_ms_min", format!("{:.6}", row.flip_ms_min)),
                ("flips_rejected", row.report.rejected.to_string()),
                (
                    "global_min_angle_after",
                    format!("{:.6}", row.after.global_min),
                ),
                (
                    "global_min_angle_before",
                    format!("{:.6}", row.before.global_min),
                ),
                (
                    "inconsistently_oriented_edges",
                    row.census.inconsistently_oriented.to_string(),
                ),
                ("interior_edges", row.census.interior.to_string()),
                ("min_angle_gain", format!("{:.6}", row.gain())),
                (
                    "min_angle_statistic",
                    "worst_decile_mean_degrees".to_string(),
                ),
                ("non_delaunay_after", row.non_delaunay_after.to_string()),
                ("non_delaunay_before", row.non_delaunay_before.to_string()),
                ("non_manifold_edges", row.census.non_manifold.to_string()),
                (
                    "p10_min_angle_after",
                    format!("{:.6}", row.after.percentile_10),
                ),
                (
                    "p10_min_angle_before",
                    format!("{:.6}", row.before.percentile_10),
                ),
                (
                    "short_edge_slivers_after",
                    row.after.short_edge_slivers.to_string(),
                ),
                (
                    "shortest_edge_over_cell",
                    format!("{:.6}", row.before.shortest_edge / row.cell_size),
                ),
                ("sliver_threshold_degrees", format!("{SLIVER_DEGREES:.6}")),
                ("triangles", row.triangles.to_string()),
                ("vertices", row.vertices.to_string()),
                ("zero_length_edges", row.census.zero_length.to_string()),
            ]);
        }
    });
}
