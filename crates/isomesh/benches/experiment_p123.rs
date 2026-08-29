//! **P-123 — where `M-318`'s 45× goes, decomposed into three terms.**
//!
//! Ticket: R-027a. Pre-registered before this harness existed.
//!
//! ```bash
//! cargo bench --bench experiment_p123
//! ```
//!
//! Writes `docs/experiments/p-123.csv`. **No clock and no hardware counter is
//! read anywhere in this file**, and that is forced by the registration rather
//! than chosen for convenience: every clause here is an integer count of values
//! or a ratio of two such counts, so the whole dataset reproduces bit for bit on
//! any machine and `M-280`'s governed-clock problem cannot reach it. `P-112` is
//! the precedent the registration names — its residual read `0.000000` in
//! retired instructions and 4–8% in nanoseconds **with a field-dependent sign**,
//! on the identical binary. This file records no `ns_` column at all.
//!
//! # What was missing
//!
//! `M-318` measured `R-027`'s **ceiling**, and then `V-45` stopped `R-027`
//! before anybody asked where that ceiling's 45× is spent. What `edit_trace`
//! produced is three totals per resolution (`FINDINGS.md:5971-5975`):
//!
//! | n | `vertices` | `buffer moved` | `geometric moved` | `keyed` |
//! |---:|---:|---:|---:|---:|
//! | 33 | 1,758 | 1,348 | 330 | **318** |
//! | 65 | 6,918 | 4,257 | 322 | **310** |
//! | 129 | 27,822 | **15,706** | 346 | **346** |
//!
//! Every one of those is a *total*. None of them says which of the 15,706 slots
//! at 129³ changed because the surface moved there, which changed because a cell
//! **earlier in traversal order** emitted a different number of vertices, and
//! which changed because of **emission order alone**. That distinction is the
//! whole decision: only the third is recoverable by a canonical reorder at
//! emission, which is the one shape of `R-027` that leaves `extract_into` a pure
//! function of its inputs. The other two are recoverable only by a persistent
//! edge → slot map, which is state carried across extractions.
//!
//! The ceiling is also not a constant 45×, and saying so costs nothing: it is
//! `buffer_moved(n) / keyed(n)`, which reads **4.239× / 13.732× / 45.393×** over
//! the three resolutions because the numerator grows with the `O(n²)` vertex
//! count and the denominator is flat in `n`. `M-318` quotes the 129³ figure, and
//! this harness inherits that row and no other.
//!
//! # What this row must not reopen
//!
//! `R-027`'s only working shape is the one `V-45` refused, and the reason is not
//! an API cost. `validate/determinism.rs:268-272` runs `check_determinism` three
//! times, the third into a **reused** buffer, under a doc comment saying exactly
//! why — *"one reused buffer to catch output that depends on the buffer's prior
//! state … nothing else checks that it survives being driven that way."* A
//! persistent edge → slot map **is** output that depends on the buffer's prior
//! state, so landing it would turn that third run from a gate into the thing it
//! was built to catch. **No clause in this row proposes it, and no column here
//! is evidence for it.** The decomposition can only ever say whether the
//! *reorder* — a pure function of the grid, carried across no call — is worth
//! writing; if it is not, `R-027` closes and the answer is *do not build this*.
//!
//! # SHARE, recomputed before any number below
//!
//! **This row moves nothing and proposes no source change**, so it has no
//! extraction share, no Amdahl ceiling and no clause that is or may become a
//! speedup claim — `✗51`'s bar does not apply. What stands in a share's place is
//! the denominator each clause is a fraction *of*, recomputed from `M-318`'s own
//! 129³ row before this harness ran:
//!
//! - **The irreducible floor is 346 of 15,706 = 2.20298%.** That is what a
//!   grid-edge naming leaves behind, and no mechanism can go below it.
//! - **The share any mechanism could remove is 15,360 of 15,706 = 97.79702%.**
//!   That is the whole prize, and it is the number `M-318` reported as 45×.
//! - **The share a canonical reorder at emission could remove is
//!   `churn_order_only / churn_total`**, which is `C2`'s column and is
//!   *unmeasured* until this run. `M-318`'s own note on its first shape —
//!   *"Stable order … Does not help. A crossing appearing still shifts every
//!   index after it"* — is a prediction about that column, not a measurement of
//!   it, and `C2`'s 50% bar is where the two part company.
//!
//! # The three arms, one build, one run (`M-281`)
//!
//! All three arms sweep the **same** edge population with the **same** crossing
//! positions in one binary, so a difference between them is the emission order
//! and nothing else. Only the order differs; the arms are permutations of one
//! arrangement, asserted so.
//!
//! | arm | emission order | what it is for |
//! |---|---|---|
//! | `scan` | cells in `z,y,x` scan order, and within a cell the case table's own first-reference order — the shipped rule, transcribed from `marching_cubes/mod.rs:254-378` and asserted **bit-identical** to `extract_into` | the subject: `M-318`'s churn, decomposed |
//! | `canonical` | the grid-edge key `lo_sample * 3 + axis`, ascending — `M-318`'s first shape, *"emit in grid-edge order rather than cell-scan order"* | the vacuity control's zero half: it must drive `churn_order_only` to **exactly 0** |
//! | `permuted` | scan order bucketed by the **first-emitting cell's cut-edge count**, a field-dependent key | `M-44`'s other half: the same detector must be shown able to read **non-zero** |
//!
//! The `scan` arm is the instrument check. It is this file's own transcription of
//! the march, the edge cache and `edge_position`, and if it were not bit-identical
//! to the shipped extractor then a difference in another arm could be a
//! transcription error rather than the ordering under test — `P-61`'s rule for a
//! second copy of an instrument, and `experiment_p101.rs`'s `edge_slot` arm is
//! the model. `replica_bit_identical` is a column and an `assert!`.
//!
//! # The decomposition, and why it is keyed on the edge
//!
//! `churn_total` is `M-318`'s `buffer_moved`: the slots at which the two
//! extractions' position buffers differ, plus the length difference, compared by
//! **index** (`edit_trace.rs:240-251`). The *classification* of each differing
//! slot is keyed on the **edge**, never on the position. `M-318`
//! (`FINDINGS.md:5968-5969`) states the reason and it is transcribed here
//! verbatim: *position-keying makes the answer equal `geometric_moved` by
//! construction and would measure nothing.* So "is this value bit-identical to
//! one in the previous extraction" is answered by asking whether the **same grid
//! edge** carried it, not by matching float triples against a set.
//!
//! The three terms, fixed at registration:
//!
//! - **`churn_geometric`** — a value that changed because its crossing moved:
//!   appeared, vanished, or still there with a moved root. This is exactly
//!   `edit_trace.rs:210-235`'s `keyed_moved` predicate, transcribed and split
//!   three ways, and it is a **count of values** rather than of slots. That is
//!   what `C1`'s own text says it is (*"these are counts of values, not times"*)
//!   and it is what the vacuity control forces: the control requires the
//!   geometric term to be *unchanged* by a reorder, and a count of crossings is
//!   a function of the two value arrays alone, so it is order-invariant by
//!   construction. A slot-level reading is not — `churn_geometric_slots` is
//!   recorded beside it and does move between arms, which is the measurement
//!   that shows why the crossing-level reading is the one the control admits.
//! - **`churn_predecessor_shift`** — a differing slot whose occupants on both
//!   sides are unchanged crossings and whose **rank among the survivors** is the
//!   same in both extractions. Its slot moved and its neighbours did not, so the
//!   only thing that changed is how many values were emitted before it: a cell
//!   earlier in traversal order emitted a different number of vertices. A
//!   reorder cannot recover this; only a persistent slot map can.
//! - **`churn_order_only`** — a differing slot whose occupants are unchanged
//!   crossings and whose rank among the survivors **did** change, so some other
//!   surviving crossing crossed it. This is the emission order permuting values
//!   that did not move, and it is exactly what a canonical reorder removes.
//!
//! `churn_geometric_slots + churn_predecessor_shift + churn_order_only ==
//! churn_total` **exactly**, asserted per arm: the cascade is tested in the
//! registration's own order and every differing slot lands in exactly one
//! bucket. `residual_share` is therefore an exact integer quantity and not a
//! timing residual:
//!
//! ```text
//! churn_residual = churn_total - (churn_geometric + predecessor_shift + order_only)
//!                = churn_geometric_slots - churn_geometric
//!                = edges_moved - geometric_slot_collisions
//! ```
//!
//! Both identities are checkable from the recorded columns, and the second names
//! the residual's whole content: a crossing that *moved* claims a slot in each
//! extraction and is one changed value, while a slot that loses a vanished
//! crossing and gains an appeared one hosts two changed values and is one slot.
//! `changed_edge_slots_before` and `changed_edge_slots_after` are recorded so a
//! reader can rebuild the arithmetic without rerunning anything. A residual
//! outside those two effects would be a decomposition defect, which is what
//! `C1`'s 5% bar exists to catch.
//!
//! # Which unit carries each verdict, and why
//!
//! - **`C1`** — `residual_share`, the absolute value of an integer difference
//!   over an integer count, **per row** (`R-085`'s discipline, and the
//!   registration asks for it per row in as many words). Machine-independent.
//! - **`C2`** — `order_only_share`, a ratio of two integer counts of values.
//!   Registered as *"at least 50% … on at least one reference field"*, so
//!   `c2_holds` is that disjunction over the whole run and reads the same on
//!   every row; `order_only_share` is the per-row number it is taken over, and
//!   `order_only_share_max` and `c2_best_field` name where the maximum was.
//! - **`C3`** — an integer equality against `M-318`'s published `keyed` column,
//!   318 / 310 / 346. Scored only on the `sphere` rows, because those are
//!   `M-318`'s fixture; `noise_cavity` has no published comparand and its
//!   `geometric_matches_m318` and `c3_holds` read `vacuous` rather than a
//!   verdict a clause could not have reached (`docs/experiments.md:237-242`).
//!
//! Nothing here is a cost or a wall-clock ratio, so `M-280` and `✗24` are
//! satisfied by there being no such column to gate on.
//!
//! # The fixtures
//!
//! `sphere` is `edit_trace`'s own fixture, transcribed constant for constant —
//! a sphere of radius 1.2 centred at `[2, 2, 2]` over `[0, 4]³` with
//! `h = 4/(n − 1)`, carved by a radius-5-sample brush at the equator that forces
//! every sample under `1.0` to exactly `1.0`. That is what makes `M-318`'s row
//! reproducible here at all, and `m318_buffer_moved_reproduced` and
//! `m318_vertices_reproduced` report whether it did.
//!
//! `noise_cavity` is `M-318`'s own named risk — *"a field whose edits change
//! crossings far from the dirty set … `noise_cavity` is the candidate and was not
//! run"* — and it is run here because `C2` asks for *at least one reference
//! field* and one field cannot support a disjunction. It has no equator to place
//! a brush on, so its brush is centred on the **cut edge nearest the grid
//! centre**, measured in L-infinity over samples and broken by the edge key: a
//! deterministic function of the field and the resolution and of nothing else.
//!
//! # The vacuity controls, as `assert!`s rather than columns
//!
//! A run that cannot fire aborts instead of recording a pass.
//!
//! 1. **`churn_total > 0` on every row, before any share of it is reported.**
//!    The registration's first half. A share of zero churn is not a measurement.
//! 2. **`canonical_control_order_only == 0` exactly**, with the canonical arm's
//!    edge multiset asserted equal to the `scan` arm's so the reorder is
//!    provably a permutation and the geometry provably untouched. The
//!    registration's second half: *a decomposition whose order-only term cannot
//!    reach zero is measuring the traversal rather than the ordering.*
//! 3. **`permuted_control_order_only > 0`**, which is `M-44` read in the other
//!    direction — a term that reads zero under a canonical order has to be shown
//!    able to read non-zero from the *same* detector in the *same* run, or the
//!    zero is not a measurement either.
//! 4. **`replica_bit_identical`** on both extractions of every row, before any
//!    arm's difference is attributed to anything.
//! 5. **`cut_edges == vertices`** on both extractions: `✗1`/`M-2`/`M-22`'s
//!    `V_mc = C`, which is what makes the transcribed edge cache provably
//!    complete rather than merely plausible.
//! 6. **`CASES[case].centroids == 0` for all 256 cases**, checked once up front.
//!    `FaceAmbiguity::Separate` and `InteriorAmbiguity::Ignore` are
//!    `MarchingCubes::new`'s defaults (`mod.rs:110-112`), so the transcription
//!    reads `CASES` directly and never reaches the decider or the trilinear
//!    path; this is the check that no cell can want a cycle centroid the replica
//!    would not emit.
//!
//! # Transcription discipline
//!
//! Every mechanism this file needs out of a private module is **copied here,
//! with the source line it came from on the row that uses it**.
//! `crates/isomesh/src/` is read-only for this row — the registration says
//! *NO NEW FIELD, NO NEW EXTRACTOR, NO SOURCE CHANGE* — and `edit_trace.rs` is
//! read-only too, so its fixture, brush and counters are transcribed rather than
//! called. A copy whose line number is on the row that uses it is auditable in a
//! way a `pub` would not be (`experiment_p117.rs:53-56`).

// Exact equality throughout: every clause here asks whether a value *changed*,
// and a tolerance would call an edit smaller than it was
// (`edit_trace.rs:35-36`, transcribed with its reason).
#![allow(
    clippy::float_cmp,
    reason = "the question is whether a value changed at all"
)]

mod common;

use std::collections::HashMap;

use isomesh::construct::SampledField;
use isomesh::extractor::Extractor;
use isomesh::fields::{ReferenceField, noise_cavity};
use isomesh::marching_cubes::MarchingCubes;
use isomesh::marching_cubes::table::{
    CASES, EDGE_AXIS, EDGE_CORNERS, NO_EDGE, corner_inside, edge_offset, is_centroid, is_inside,
};
use isomesh::{MeshBuffer, RuntimeShape3, Sdf, Shape3};

/// Samples per axis. `edit_trace.rs:58` — the brush is the same at all three,
/// which is what makes the `keyed` column's flatness in `n` legible.
const RESOLUTIONS: [u32; 3] = [33, 65, 129];

/// The two fields. `sphere` is `M-318`'s fixture and carries `C3`;
/// `noise_cavity` is `M-318`'s own named risk and is what lets `C2`'s
/// *"at least one reference field"* be a disjunction rather than a single row.
const FIELDS: [&str; 2] = ["sphere", "noise_cavity"];

/// Brush radius, in samples. `edit_trace.rs:61` — fixed across resolutions.
const BRUSH_RADIUS: f64 = 5.0;

/// `C1`'s bar: the residual must be under 5% of the churn.
const RESIDUAL_BAR: f64 = 0.05;

/// `C2`'s bar: the order-only term must reach half the churn somewhere.
const ORDER_ONLY_BAR: f64 = 0.50;

/// `M-318`'s own table, `FINDINGS.md:5971-5975`, as the comparand `C3` is scored
/// against and the fixture check `m318_*_reproduced` is taken from. Columns:
/// samples per axis, `vertices`, `buffer moved`, `geometric moved`, `keyed`.
const M318: [[u64; 5]; 3] = [
    [33, 1_758, 1_348, 330, 318],
    [65, 6_918, 4_257, 322, 310],
    [129, 27_822, 15_706, 346, 346],
];

/// Not applicable — `experiment_p101.rs:1388`'s spelling, kept so the two files
/// read the same way.
const NA: &str = "";

/// The word for a clause that could not have fired, which is not the same thing
/// as one that held (`docs/experiments.md:237-242`).
const VACUOUS: &str = "vacuous";

// ─── the crate's own arithmetic, transcribed ────────────────────────────────
//
// `crate::cube` is a private module (`lib.rs:143`), so these are copies rather
// than calls. Each carries the line it came from, and the `scan` arm asserts the
// whole assembly bit-identical to the shipped extractor, which is the only thing
// that makes a second copy of an instrument admissible.

/// `cube::corner_offset` (`cube.rs:149`), which is `pub(crate)`.
#[inline]
fn corner_offset(corner: u8) -> [u32; 3] {
    [
        u32::from(corner & 1),
        u32::from((corner >> 1) & 1),
        u32::from((corner >> 2) & 1),
    ]
}

/// `cube::place` (`cube.rs:233`), which is `pub(crate)`.
#[inline]
fn place(lo: f64, hi: f64, d: f64) -> f64 {
    (lo + hi) * 0.5 + (hi - lo) * d
}

/// `marching_cubes::corner_position` (`mod.rs:758-766`), which is private.
#[inline]
fn corner_position(base: [u32; 3], corner: u8, origin: [f64; 3], h: f64) -> [f64; 3] {
    let o = corner_offset(corner);
    [
        origin[0] + h * f64::from(base[0] + o[0]),
        origin[1] + h * f64::from(base[1] + o[1]),
        origin[2] + h * f64::from(base[2] + o[2]),
    ]
}

/// `marching_cubes::edge_position` (`mod.rs:631-663`) at
/// `crossing_refinement == 0`, which is `MarchingCubes::new`'s default
/// (`mod.rs:112`). `refine_crossing` returns `d0` unchanged at zero steps
/// (`mod.rs:707-709`), so the `Sdf` bisection drops out of the transcription
/// entirely and no field evaluation is needed to place a crossing.
#[inline]
fn edge_position(
    base: [u32; 3],
    edge: u8,
    corner_value: &[f64; 8],
    origin: [f64; 3],
    h: f64,
) -> [f64; 3] {
    let [lo_corner, hi_corner] = EDGE_CORNERS[edge as usize];
    let a = corner_value[lo_corner as usize];
    let b = corner_value[hi_corner as usize];
    // `cube::edge_offset` is `pub` (`cube.rs:222`, and `cube.rs:214-220` says
    // why), so this one is a call and not a copy.
    let d = edge_offset(a, b);
    let lo = corner_position(base, lo_corner, origin, h);
    let hi = corner_position(base, hi_corner, origin, h);
    [
        place(lo[0], hi[0], d),
        place(lo[1], hi[1], d),
        place(lo[2], hi[2], d),
    ]
}

/// The crate's own global edge name: `lo_sample * 3 + axis`, which is the key
/// into `edge_vertices` at `marching_cubes/mod.rs:602-604`. `M-318`'s
/// counterfactual is denominated in exactly this name and so is every
/// classification below.
#[inline]
fn edge_key(shape: &RuntimeShape3, base: [u32; 3], edge: u8) -> u32 {
    let axis = u32::from(EDGE_AXIS[edge as usize]);
    let o = corner_offset(EDGE_CORNERS[edge as usize][0]);
    let lo_sample = shape.linearize([base[0] + o[0], base[1] + o[1], base[2] + o[2]]);
    lo_sample * 3 + axis
}

/// The eight-bit sign pattern of a cell, which is what the case table is indexed
/// by. `edit_trace.rs:77-88`, itself the crate's rule at `mod.rs:259-268`.
fn case_index(values: &[f64], shape: &RuntimeShape3, cell: [u32; 3]) -> u8 {
    let mut mask = 0u8;
    for c in 0..8u8 {
        let o = corner_offset(c);
        let i = shape.linearize([cell[0] + o[0], cell[1] + o[1], cell[2] + o[2]]);
        if is_inside(values[i as usize]) {
            mask |= 1 << c;
        }
    }
    mask
}

/// How many of a case's twelve cube edges are cut. The `permuted` control arm's
/// bucket key, and the only field-dependent quantity in any emission order here.
fn cut_edges_per_case() -> [u8; 256] {
    let mut out = [0u8; 256];
    for (case, count) in out.iter_mut().enumerate() {
        let case = case as u8;
        let mut cut = 0u8;
        for &[lo, hi] in &EDGE_CORNERS {
            if corner_inside(case, lo) != corner_inside(case, hi) {
                cut += 1;
            }
        }
        *count = cut;
    }
    out
}

// ─── one transcribed extraction ─────────────────────────────────────────────

/// One extraction of the march, transcribed from `marching_cubes/mod.rs:254-378`
/// for `MarchingCubes::new`'s default configuration.
///
/// Holds what the shipped path throws away: the edge each vertex sits on, and
/// therefore the emission order as a sequence of *edges* rather than of floats.
/// `M-318` is about precisely that discarded name (`FINDINGS.md:5965-5966` —
/// *"the crate already keys its edge cache that way and discards the name when
/// packing"*).
struct Replica {
    /// `lo_sample * 3 + axis` → slot, `u32::MAX` for "no vertex on this edge".
    /// The shipped structure verbatim: `mod.rs:97` declares it, `:251` fills it
    /// with the same sentinel, `:606-609` probes it and `:620-622` writes it.
    slot_of: Vec<u32>,
    /// slot → edge key. The emission order itself.
    order: Vec<u32>,
    /// slot → position, which must be bit-identical to `MeshBuffer::positions`.
    position: Vec<[f64; 3]>,
    /// slot → the case index of the cell that first emitted it. Read by the
    /// `permuted` control arm and by nothing else.
    first_case: Vec<u8>,
    /// Triangles emitted, for the record.
    triangles: u64,
}

/// Sweep the grid in scan order and record what the shipped extractor would
/// emit, in the order it would emit it.
fn march(values: &[f64], shape: &RuntimeShape3, n: u32, origin: [f64; 3], h: f64) -> Replica {
    let mut r = Replica {
        slot_of: vec![u32::MAX; values.len() * 3],
        order: Vec::new(),
        position: Vec::new(),
        first_case: Vec::new(),
        triangles: 0,
    };
    // `mod.rs:254-256`: z, then y, then x, over cells rather than samples.
    for z in 0..n - 1 {
        for y in 0..n - 1 {
            for x in 0..n - 1 {
                let base = [x, y, z];
                // `mod.rs:259-268`: the case and the eight corner values in one
                // pass, because `edge_position` needs the values again.
                let mut case = 0u8;
                let mut corner_value = [0.0f64; 8];
                for (c, slot) in corner_value.iter_mut().enumerate() {
                    let o = corner_offset(c as u8);
                    let s = shape.linearize([base[0] + o[0], base[1] + o[1], base[2] + o[2]]);
                    let v = values[s as usize];
                    *slot = v;
                    if is_inside(v) {
                        case |= 1 << c;
                    }
                }
                // `mod.rs:278-311` collapses to this under the defaults:
                // `FaceAmbiguity::Separate` makes `ambiguous` zero, a zero
                // `ambiguous` makes `mask` zero, and a zero `ambiguous` takes
                // the `CASES` branch and never the decider or the trilinear one.
                let entry = CASES[case as usize];
                if entry.count == 0 {
                    continue;
                }
                for tri in &entry.triangles[..entry.count as usize] {
                    r.triangles += 1;
                    // `mod.rs:357-376`: three codes per triangle, in order, each
                    // resolved through the edge cache. Centroids are ruled out
                    // for every case up front, so this is the whole rule.
                    for &code in tri {
                        assert!(
                            code != NO_EDGE && !is_centroid(code),
                            "case {case} names code {code}, which is not a cube edge — the \
                             transcription at mod.rs:357-376 is incomplete"
                        );
                        let key = edge_key(shape, base, code);
                        if r.slot_of[key as usize] == u32::MAX {
                            r.slot_of[key as usize] = r.order.len() as u32;
                            r.order.push(key);
                            r.position
                                .push(edge_position(base, code, &corner_value, origin, h));
                            r.first_case.push(case);
                        }
                    }
                }
            }
        }
    }
    r
}

impl Replica {
    /// The position on `edge`, which must be present.
    #[inline]
    fn position_of(&self, edge: u32) -> [f64; 3] {
        self.position[self.slot_of[edge as usize] as usize]
    }

    /// Is this transcription the shipped arithmetic, bit for bit?
    fn matches(&self, mesh: &MeshBuffer<f64>) -> bool {
        self.position.len() == mesh.positions.len()
            && self
                .position
                .iter()
                .zip(&mesh.positions)
                .all(|(a, b)| a == b)
    }
}

/// The shipped extraction, `edit_trace.rs:90-97` transcribed.
fn mesh(values: &[f64], shape: &RuntimeShape3, origin: [f64; 3], h: f64) -> MeshBuffer<f64> {
    let mut out = MeshBuffer::<f64>::new();
    let field = SampledField::new(values, shape, origin, h).expect("one value per sample");
    MarchingCubes::<f64>::new()
        .extract_into(&field, shape, origin, h, &mut out)
        .expect("the reference grids fit the index space");
    out
}

// ─── M-318's counterfactual, transcribed and split ──────────────────────────

/// Which grid edges' crossings moved, and how.
///
/// `edit_trace.rs:210-235` verbatim, with the single `keyed_moved` counter split
/// into its three cases. Counted from the value arrays alone with no extractor
/// involved, and **keyed on the edge**: `M-318` (`FINDINGS.md:5968-5969`) is
/// explicit that position-keying makes the answer equal `geometric_moved` by
/// construction and would measure nothing.
struct Census {
    /// Edge key → did this edge's crossing appear, vanish or move.
    changed: Vec<bool>,
    appeared: u64,
    vanished: u64,
    moved: u64,
    cut_before: u64,
    cut_after: u64,
}

impl Census {
    /// The geometric term: values that changed because a crossing moved. One per
    /// crossing, which is `M-318`'s `keyed` column exactly.
    fn geometric(&self) -> u64 {
        self.appeared + self.vanished + self.moved
    }
}

fn census(before: &[f64], after: &[f64], shape: &RuntimeShape3, n: u32) -> Census {
    let mut c = Census {
        changed: vec![false; before.len() * 3],
        appeared: 0,
        vanished: 0,
        moved: 0,
        cut_before: 0,
        cut_after: 0,
    };
    for z in 0..n {
        for y in 0..n {
            for x in 0..n {
                for axis in 0..3usize {
                    let mut to = [x, y, z];
                    to[axis] += 1;
                    if to[axis] >= n {
                        continue;
                    }
                    let i = shape.linearize([x, y, z]) as usize;
                    let j = shape.linearize(to) as usize;
                    let (ai, aj) = (before[i], before[j]);
                    let (bi, bj) = (after[i], after[j]);
                    let had = is_inside(ai) != is_inside(aj);
                    let has = is_inside(bi) != is_inside(bj);
                    if had {
                        c.cut_before += 1;
                    }
                    if has {
                        c.cut_after += 1;
                    }
                    // `edit_trace.rs:227-231`: appeared or vanished, or still
                    // there and the root moved because an endpoint value did.
                    let changed = match (had, has) {
                        (false, true) => {
                            c.appeared += 1;
                            true
                        }
                        (true, false) => {
                            c.vanished += 1;
                            true
                        }
                        (true, true) if ai != bi || aj != bj => {
                            c.moved += 1;
                            true
                        }
                        _ => false,
                    };
                    if changed {
                        c.changed[i * 3 + axis] = true;
                    }
                }
            }
        }
    }
    c
}

// ─── the arms, as permutations of one arrangement ───────────────────────────

/// One emission order over one extraction's edge population.
struct Arrangement {
    /// slot → edge key.
    order: Vec<u32>,
    /// Edge key → its index among the **survivors**, in this order's slot order.
    /// A survivor is an edge whose crossing did not change, so it is present in
    /// both extractions carrying a bit-identical value.
    surv_rank: HashMap<u32, u32>,
}

impl Arrangement {
    fn new(order: Vec<u32>, changed: &[bool]) -> Self {
        let mut surv_rank = HashMap::with_capacity(order.len());
        for &e in &order {
            if !changed[e as usize] {
                let rank = surv_rank.len() as u32;
                surv_rank.insert(e, rank);
            }
        }
        Self { order, surv_rank }
    }
}

/// `M-318`'s first shape: *"emit in grid-edge order rather than cell-scan
/// order"*. The key is intrinsic to the edge, so the survivors' relative order
/// cannot depend on which other edges exist — which is why this arm is the
/// control that must drive `churn_order_only` to exactly zero.
fn canonical_order(r: &Replica) -> Vec<u32> {
    let mut order = r.order.clone();
    order.sort_unstable();
    order
}

/// Scan order re-bucketed by the first-emitting cell's cut-edge count.
///
/// A legitimate deterministic emission rule and a **field-dependent** one: an
/// edit that changes a cell's case changes its cut-edge count, which moves every
/// edge that cell first emitted into a different bucket and so permutes
/// survivors against survivors elsewhere. That is what makes this the control
/// for the other direction — `churn_order_only` reading zero is only a
/// measurement if the same detector can be made to read non-zero.
fn permuted_order(r: &Replica, cut_of_case: &[u8; 256]) -> Vec<u32> {
    let mut slots: Vec<u32> = (0..r.order.len() as u32).collect();
    slots.sort_by_key(|&s| (cut_of_case[r.first_case[s as usize] as usize], s));
    slots.into_iter().map(|s| r.order[s as usize]).collect()
}

// ─── the decomposition ──────────────────────────────────────────────────────

/// One arm's decomposition of one churn.
#[derive(Default)]
struct Terms {
    /// `M-318`'s `buffer_moved`: slots at which the two buffers differ, plus the
    /// length difference (`edit_trace.rs:240-251`).
    total: u64,
    /// Differing slots hosting a changed crossing on either side. Order-
    /// dependent, and recorded to show why the registered geometric term is
    /// counted per crossing instead.
    geometric_slots: u64,
    predecessor_shift: u64,
    order_only: u64,
    /// Differing slots hosting a changed crossing on **both** sides — one slot,
    /// two changed values, and the whole content of a negative residual.
    collisions: u64,
    changed_slots_before: u64,
    changed_slots_after: u64,
}

/// Classify every differing slot into exactly one of the three buckets, in the
/// registration's own cascade order: geometric, then predecessor-shift, then
/// order-only.
///
/// The diff itself is by index over positions, which is `M-318`'s quantity. The
/// classification is by edge, which is `M-318`'s rule.
fn decompose(
    a: &Arrangement,
    b: &Arrangement,
    before: &Replica,
    after: &Replica,
    changed: &[bool],
) -> Terms {
    let mut t = Terms::default();
    for i in 0..a.order.len().max(b.order.len()) {
        let oa = a.order.get(i).copied();
        let ob = b.order.get(i).copied();
        // Present on both sides carrying the same bits is not churn. Beyond one
        // buffer's end always is, which is `edit_trace.rs:251`'s length term
        // arriving slot by slot rather than as an absolute difference.
        if let (Some(ea), Some(eb)) = (oa, ob)
            && before.position_of(ea) == after.position_of(eb)
        {
            continue;
        }
        t.total += 1;
        let ca = oa.is_some_and(|e| changed[e as usize]);
        let cb = ob.is_some_and(|e| changed[e as usize]);
        if ca {
            t.changed_slots_before += 1;
        }
        if cb {
            t.changed_slots_after += 1;
        }
        if ca && cb {
            t.collisions += 1;
        }
        if ca || cb {
            t.geometric_slots += 1;
            continue;
        }
        // Neither occupant's crossing moved, so whatever value is here is
        // bit-identical to one in the other extraction and only its slot moved.
        let edge = ob
            .or(oa)
            .expect("a differing slot has an occupant on at least one side");
        if a.surv_rank.get(&edge) == b.surv_rank.get(&edge) {
            t.predecessor_shift += 1;
        } else {
            t.order_only += 1;
        }
    }
    assert_eq!(
        t.geometric_slots + t.predecessor_shift + t.order_only,
        t.total,
        "the cascade left a differing slot unclassified, so it is not a partition"
    );
    t
}

// ─── the fixtures ───────────────────────────────────────────────────────────

/// One field, one resolution, one edit.
struct Fixture {
    field: &'static str,
    n: u32,
    h: f64,
    origin: [f64; 3],
    shape: RuntimeShape3,
    before: Vec<f64>,
    after: Vec<f64>,
    edits: u64,
    brush_centre: [u32; 3],
}

/// `edit_trace.rs:135-167`: carve a spherical brush of [`BRUSH_RADIUS`] samples,
/// forcing every sample under `1.0` to exactly `1.0` — outside, by
/// `cube::is_inside`'s rule that a sample of exactly zero is outside and
/// anything non-negative with it (`cube.rs:170-173`).
fn carve(before: &[f64], shape: &RuntimeShape3, n: u32, centre: [u32; 3]) -> (Vec<f64>, u64) {
    let mut after = before.to_vec();
    let mut edits = 0u64;
    let r = BRUSH_RADIUS.ceil() as i64;
    for dz in -r..=r {
        for dy in -r..=r {
            for dx in -r..=r {
                if (dx * dx + dy * dy + dz * dz) as f64 > BRUSH_RADIUS * BRUSH_RADIUS {
                    continue;
                }
                let p = [
                    i64::from(centre[0]) + dx,
                    i64::from(centre[1]) + dy,
                    i64::from(centre[2]) + dz,
                ];
                if p.iter().any(|c| *c < 0 || *c >= i64::from(n)) {
                    continue;
                }
                #[allow(
                    clippy::cast_sign_loss,
                    reason = "the bound above proves every component is non-negative"
                )]
                let i = shape.linearize([p[0] as u32, p[1] as u32, p[2] as u32]) as usize;
                if after[i] < 1.0 {
                    after[i] = 1.0;
                    edits += 1;
                }
            }
        }
    }
    (after, edits)
}

/// `edit_trace`'s own fixture, transcribed constant for constant
/// (`edit_trace.rs:113-167`). This is what makes `M-318`'s row reproducible
/// here, and `C3` is scored against it and against nothing else.
fn sphere_fixture(n: u32) -> Fixture {
    let shape = RuntimeShape3::new([n; 3]).expect("a reference grid fits u32");
    let h = 4.0 / f64::from(n - 1);
    let centre = [2.0_f64; 3];
    let radius = 1.2_f64;
    let mut before = Vec::with_capacity(shape.element_count());
    for z in 0..n {
        for y in 0..n {
            for x in 0..n {
                let p = [f64::from(x) * h, f64::from(y) * h, f64::from(z) * h];
                let d = [p[0] - centre[0], p[1] - centre[1], p[2] - centre[2]];
                before.push((d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt() - radius);
            }
        }
    }
    // `edit_trace.rs:138`: the brush sits at the sphere's equator, so it
    // straddles the surface and changes topology rather than rewriting values
    // nobody reads.
    #[allow(
        clippy::cast_possible_truncation,
        clippy::cast_sign_loss,
        reason = "edit_trace.rs:138's own expression, transcribed"
    )]
    let brush_centre = [n / 2 + (radius / h) as u32, n / 2, n / 2];
    let (after, edits) = carve(&before, &shape, n, brush_centre);
    Fixture {
        field: "sphere",
        n,
        h,
        origin: [0.0; 3],
        shape,
        before,
        after,
        edits,
        brush_centre,
    }
}

/// `M-318`'s own named risk, run: *"a field whose edits change crossings far
/// from the dirty set … `noise_cavity` is the candidate and was not run."*
///
/// Sampled over its declared domain with `sample_grid`'s exact expression
/// (`sdf.rs:183-187`), and edited by the same brush. It has no equator, so the
/// brush is centred on the cut edge nearest the grid centre — see
/// [`brush_at_nearest_crossing`].
fn noise_cavity_fixture(n: u32) -> Fixture {
    let field = noise_cavity::<f64>();
    let (lo, hi) = field.domain();
    let shape = RuntimeShape3::new([n; 3]).expect("a reference grid fits u32");
    let h = (hi[0] - lo[0]) / f64::from(n - 1);
    let mut before = Vec::with_capacity(shape.element_count());
    for z in 0..n {
        for y in 0..n {
            for x in 0..n {
                before.push(field.sample([
                    lo[0] + h * f64::from(x),
                    lo[1] + h * f64::from(y),
                    lo[2] + h * f64::from(z),
                ]));
            }
        }
    }
    let brush_centre = brush_at_nearest_crossing(&before, &shape, n);
    let (after, edits) = carve(&before, &shape, n, brush_centre);
    Fixture {
        field: "noise_cavity",
        n,
        h,
        origin: lo,
        shape,
        before,
        after,
        edits,
        brush_centre,
    }
}

/// The lower sample of the cut edge nearest the grid centre, in L-infinity over
/// samples, ties broken by the edge key.
///
/// A deterministic function of the field and the resolution and of nothing else,
/// which is the property `edit_trace`'s sphere gets for free from its own radius
/// and a noise field cannot.
fn brush_at_nearest_crossing(values: &[f64], shape: &RuntimeShape3, n: u32) -> [u32; 3] {
    let mid = n / 2;
    let mut best = (u32::MAX, u32::MAX, [mid; 3]);
    for z in 0..n {
        for y in 0..n {
            for x in 0..n {
                for axis in 0..3usize {
                    let mut to = [x, y, z];
                    to[axis] += 1;
                    if to[axis] >= n {
                        continue;
                    }
                    let i = shape.linearize([x, y, z]);
                    let j = shape.linearize(to);
                    if is_inside(values[i as usize]) == is_inside(values[j as usize]) {
                        continue;
                    }
                    let far = x.abs_diff(mid).max(y.abs_diff(mid)).max(z.abs_diff(mid));
                    let key = i * 3 + axis as u32;
                    if (far, key) < (best.0, best.1) {
                        best = (far, key, [x, y, z]);
                    }
                }
            }
        }
    }
    assert!(
        best.0 != u32::MAX,
        "the field has no crossing at all, so there is no surface to disturb"
    );
    best.2
}

/// `edit_trace.rs:169-195`: the cells an edit can possibly re-run, and the
/// subset whose case index actually moved. Recorded because `case_changed` is
/// the population the `permuted` control arm's bucket key depends on, so a zero
/// there would explain a zero from that control.
fn dirty_cells(before: &[f64], after: &[f64], shape: &RuntimeShape3, n: u32) -> (u64, u64) {
    let mut dirty = 0u64;
    let mut case_changed = 0u64;
    for z in 0..n - 1 {
        for y in 0..n - 1 {
            for x in 0..n - 1 {
                let cell = [x, y, z];
                let touched = (0..8u8).any(|c| {
                    let o = corner_offset(c);
                    let i =
                        shape.linearize([cell[0] + o[0], cell[1] + o[1], cell[2] + o[2]]) as usize;
                    before[i] != after[i]
                });
                if touched {
                    dirty += 1;
                    if case_index(before, shape, cell) != case_index(after, shape, cell) {
                        case_changed += 1;
                    }
                }
            }
        }
    }
    (dirty, case_changed)
}

/// `M-318`'s published row for this resolution, or `None`.
fn m318_row(n: u32) -> Option<&'static [u64; 5]> {
    M318.iter().find(|row| row[0] == u64::from(n))
}

fn main() {
    if !std::env::args().any(|arg| arg == "--bench") {
        return;
    }
    let prereg = isomesh::experiment!("P-123");
    common::experiment::run(prereg, |run| {
        // **Vacuity control 6, and it is the one that makes the transcription
        // complete rather than plausible.** Under `MarchingCubes::new`'s
        // defaults the march reads `CASES` directly (`mod.rs:307-311`), and a
        // case wanting a cycle centroid would emit a vertex the replica has no
        // rule for. `table.rs:126-131` says plain Marching Cubes never reaches
        // that path; this checks it over all 256 rather than believing it.
        for (case, entry) in CASES.iter().enumerate() {
            assert_eq!(
                entry.centroids, 0,
                "CASES[{case}] wants {} cycle centroids, which the transcription of \
                 mod.rs:320-355 does not implement",
                entry.centroids
            );
        }
        let cut_of_case = cut_edges_per_case();

        println!(
            "{:>13} {:>5} {:>6} {:>9} {:>7} {:>9} {:>7} {:>9} {:>9} {:>8} {:>8}",
            "field",
            "n",
            "edits",
            "churn",
            "geom",
            "predshift",
            "order",
            "resid",
            "order/tot",
            "canon",
            "perm"
        );

        let mut rows: Vec<Vec<(&'static str, String)>> = Vec::new();
        let mut best_order_share = 0.0f64;
        let mut best_order_field = "none";

        for field in FIELDS {
            for n in RESOLUTIONS {
                let fx = if field == "sphere" {
                    sphere_fixture(n)
                } else {
                    noise_cavity_fixture(n)
                };

                // ── the shipped path, and the transcription beside it ────────
                let mesh_before = mesh(&fx.before, &fx.shape, fx.origin, fx.h);
                let mesh_after = mesh(&fx.after, &fx.shape, fx.origin, fx.h);
                let before = march(&fx.before, &fx.shape, fx.n, fx.origin, fx.h);
                let after = march(&fx.after, &fx.shape, fx.n, fx.origin, fx.h);
                let replica_ok = before.matches(&mesh_before) && after.matches(&mesh_after);
                assert!(
                    replica_ok,
                    "{field} at {n}³: the transcribed march is not the shipped arithmetic, so \
                     no arm's difference is attributable to its emission order (P-61's rule)"
                );

                // ── M-318's counterfactual, recomputed here ─────────────────
                let c = census(&fx.before, &fx.after, &fx.shape, fx.n);
                assert_eq!(
                    c.cut_before,
                    before.order.len() as u64,
                    "{field} at {n}³: {} cut edges but {} vertices — ✗1/M-2/M-22's V_mc = C \
                     is what makes the transcribed edge cache provably complete",
                    c.cut_before,
                    before.order.len()
                );
                assert_eq!(
                    c.cut_after,
                    after.order.len() as u64,
                    "{field} at {n}³: V_mc = C fails on the edited grid"
                );
                let churn_geometric = c.geometric();

                // `edit_trace.rs:169-195`'s two counters, hoisted because the
                // permuted control's failure message names `case_changed`: a
                // control that could not permute anything would be explained by
                // no cell having changed case, and that explanation belongs in
                // the panic rather than in a later column.
                let (dirty, case_changed) = dirty_cells(&fx.before, &fx.after, &fx.shape, fx.n);

                // ── the three arms over one shared population ──────────────
                let scan_a = Arrangement::new(before.order.clone(), &c.changed);
                let scan_b = Arrangement::new(after.order.clone(), &c.changed);
                let canon_a = Arrangement::new(canonical_order(&before), &c.changed);
                let canon_b = Arrangement::new(canonical_order(&after), &c.changed);
                let perm_a = Arrangement::new(permuted_order(&before, &cut_of_case), &c.changed);
                let perm_b = Arrangement::new(permuted_order(&after, &cut_of_case), &c.changed);

                // The reorders are permutations and nothing more: same edges,
                // same count, so the geometry is provably untouched and the only
                // thing an arm can move is where a value lands.
                for (label, arm) in [
                    ("canonical", &canon_a),
                    ("permuted", &perm_a),
                    ("canonical_after", &canon_b),
                    ("permuted_after", &perm_b),
                ] {
                    let reference = if label.ends_with("after") {
                        &after.order
                    } else {
                        &before.order
                    };
                    let mut mine = arm.order.clone();
                    let mut theirs = reference.clone();
                    mine.sort_unstable();
                    theirs.sort_unstable();
                    assert_eq!(
                        mine, theirs,
                        "{field} at {n}³: the {label} arm is not a permutation of the shipped \
                         emission, so it changes more than the order"
                    );
                }

                let scan = decompose(&scan_a, &scan_b, &before, &after, &c.changed);
                let canon = decompose(&canon_a, &canon_b, &before, &after, &c.changed);
                let perm = decompose(&perm_a, &perm_b, &before, &after, &c.changed);

                // ── the vacuity controls ────────────────────────────────────
                assert!(
                    scan.total > 0,
                    "{field} at {n}³: churn_total is zero, and a share of zero churn is not a \
                     measurement (M-44). The edit changed {} samples and {} crossings",
                    fx.edits,
                    churn_geometric
                );
                assert_eq!(
                    canon.order_only, 0,
                    "{field} at {n}³: the canonical grid-edge order left {} order-only slots, \
                     so this decomposition is measuring the traversal rather than the ordering",
                    canon.order_only
                );
                assert_eq!(
                    canon.geometric_slots + canon.predecessor_shift,
                    canon.total,
                    "{field} at {n}³: the canonical control's churn is not entirely geometric \
                     plus predecessor-shift, which is what an order-only term of zero means"
                );
                assert!(
                    perm.order_only > 0,
                    "{field} at {n}³: the permuted control could not make the order-only \
                     detector read non-zero, so its zero under the canonical order is not a \
                     measurement either (M-44). {} cells changed case",
                    case_changed
                );

                // ── the numbers, and the clause verdicts ────────────────────
                let residual = scan.total as i64
                    - (churn_geometric + scan.predecessor_shift + scan.order_only) as i64;
                assert_eq!(
                    residual,
                    scan.geometric_slots as i64 - churn_geometric as i64,
                    "{field} at {n}³: the residual is not the difference between the slot-level \
                     and crossing-level geometric terms, so one of them is miscounted"
                );
                #[allow(
                    clippy::cast_precision_loss,
                    reason = "counts here are at most a few million and exact in f64"
                )]
                let residual_share = residual.unsigned_abs() as f64 / scan.total as f64;
                #[allow(
                    clippy::cast_precision_loss,
                    reason = "counts here are at most a few million and exact in f64"
                )]
                let order_only_share = scan.order_only as f64 / scan.total as f64;
                let c1 = residual_share < RESIDUAL_BAR;
                if order_only_share > best_order_share {
                    best_order_share = order_only_share;
                    best_order_field = fx.field;
                }

                let published = if fx.field == "sphere" {
                    m318_row(n)
                } else {
                    None
                };
                let (matches_m318, c3) = match published {
                    Some(row) => {
                        let ok = churn_geometric == row[4];
                        (ok.to_string(), ok.to_string())
                    }
                    None => (VACUOUS.to_string(), VACUOUS.to_string()),
                };
                let m318_cols: [String; 6] = match published {
                    Some(row) => [
                        row[1].to_string(),
                        row[2].to_string(),
                        row[3].to_string(),
                        row[4].to_string(),
                        (before.order.len() as u64 == row[1]).to_string(),
                        (scan.total == row[2]).to_string(),
                    ],
                    None => [
                        NA.to_string(),
                        NA.to_string(),
                        NA.to_string(),
                        NA.to_string(),
                        VACUOUS.to_string(),
                        VACUOUS.to_string(),
                    ],
                };

                println!(
                    "{:>13} {:>5} {:>6} {:>9} {:>7} {:>9} {:>7} {:>9} {order_only_share:>9.6} \
                     {:>8} {:>8}",
                    fx.field,
                    n,
                    fx.edits,
                    scan.total,
                    churn_geometric,
                    scan.predecessor_shift,
                    scan.order_only,
                    residual,
                    canon.order_only,
                    perm.order_only
                );

                let cells = u64::from(n - 1).pow(3);
                rows.push(vec![
                    ("field", fx.field.to_string()),
                    ("resolution", n.to_string()),
                    ("edits", fx.edits.to_string()),
                    ("churn_total", scan.total.to_string()),
                    ("churn_geometric", churn_geometric.to_string()),
                    (
                        "churn_predecessor_shift",
                        scan.predecessor_shift.to_string(),
                    ),
                    ("churn_order_only", scan.order_only.to_string()),
                    ("residual_share", format!("{residual_share:.9}")),
                    ("order_only_share", format!("{order_only_share:.9}")),
                    ("geometric_matches_m318", matches_m318),
                    ("canonical_control_order_only", canon.order_only.to_string()),
                    ("c1_holds", c1.to_string()),
                    ("c3_holds", c3),
                    // ── the decomposition, made checkable ──────────────────
                    ("churn_residual", residual.to_string()),
                    ("churn_geometric_slots", scan.geometric_slots.to_string()),
                    ("geometric_slot_collisions", scan.collisions.to_string()),
                    (
                        "changed_edge_slots_before",
                        scan.changed_slots_before.to_string(),
                    ),
                    (
                        "changed_edge_slots_after",
                        scan.changed_slots_after.to_string(),
                    ),
                    (
                        "shifted_survivors",
                        (scan.predecessor_shift + scan.order_only).to_string(),
                    ),
                    ("survivors", scan_a.surv_rank.len().to_string()),
                    ("edges_appeared", c.appeared.to_string()),
                    ("edges_vanished", c.vanished.to_string()),
                    ("edges_moved", c.moved.to_string()),
                    ("cut_edges_before", c.cut_before.to_string()),
                    ("cut_edges_after", c.cut_after.to_string()),
                    // ── the arms ───────────────────────────────────────────
                    ("canonical_control_total", canon.total.to_string()),
                    ("canonical_control_geometric", churn_geometric.to_string()),
                    (
                        "canonical_control_geometric_slots",
                        canon.geometric_slots.to_string(),
                    ),
                    (
                        "canonical_control_predecessor_shift",
                        canon.predecessor_shift.to_string(),
                    ),
                    ("permuted_control_total", perm.total.to_string()),
                    ("permuted_control_order_only", perm.order_only.to_string()),
                    (
                        "permuted_control_predecessor_shift",
                        perm.predecessor_shift.to_string(),
                    ),
                    // ── the fixture ────────────────────────────────────────
                    ("cells", cells.to_string()),
                    ("dirty_cells", dirty.to_string()),
                    ("case_changed", case_changed.to_string()),
                    ("vertices_before", mesh_before.positions.len().to_string()),
                    ("vertices_after", mesh_after.positions.len().to_string()),
                    ("triangles_before", before.triangles.to_string()),
                    ("triangles_after", after.triangles.to_string()),
                    ("brush_radius_samples", format!("{BRUSH_RADIUS:.1}")),
                    ("brush_centre_x", fx.brush_centre[0].to_string()),
                    ("brush_centre_y", fx.brush_centre[1].to_string()),
                    ("brush_centre_z", fx.brush_centre[2].to_string()),
                    // ── the instrument ─────────────────────────────────────
                    ("replica_bit_identical", replica_ok.to_string()),
                    ("m318_vertices", m318_cols[0].clone()),
                    ("m318_buffer_moved", m318_cols[1].clone()),
                    ("m318_geometric_moved", m318_cols[2].clone()),
                    ("m318_keyed_moved", m318_cols[3].clone()),
                    ("m318_vertices_reproduced", m318_cols[4].clone()),
                    ("m318_buffer_moved_reproduced", m318_cols[5].clone()),
                    ("verdict_unit", "integer_value_counts".to_string()),
                    ("c1_bar", format!("{RESIDUAL_BAR:.6}")),
                    ("c2_bar", format!("{ORDER_ONLY_BAR:.6}")),
                ]);
            }
        }

        // **C2 is a disjunction over fields, so it is scored once over the whole
        // run** — the registration says *"at least 50% … on at least one
        // reference field"*, and a per-row reading of that clause would be a
        // different clause.
        let c2 = best_order_share >= ORDER_ONLY_BAR;

        println!("\n-- aggregates --");
        println!(
            "  C1  residual_share, per row, against {RESIDUAL_BAR:.2}: {} of {} rows under it",
            rows.iter()
                .filter(|r| r.iter().any(|(k, v)| *k == "c1_holds" && v == "true"))
                .count(),
            rows.len()
        );
        println!(
            "  C2  max order_only_share {best_order_share:.9} on {best_order_field}, \
             against {ORDER_ONLY_BAR:.2} — {c2}"
        );
        let c3_tally = |want: &str| {
            rows.iter()
                .filter(|r| r.iter().any(|(k, v)| *k == "c3_holds" && v == want))
                .count()
        };
        println!(
            "  C3  churn_geometric against M-318's 318/310/346: {} held, {} falsified, \
             {} {VACUOUS} (noise_cavity has no published comparand)",
            c3_tally("true"),
            c3_tally("false"),
            c3_tally(VACUOUS)
        );
        println!(
            "  vacuity: canonical control drove order-only to 0 on every row, permuted \
             control drove it above 0 on every row, replica bit-identical on every row"
        );
        println!("  no clock and no counter was read: verdict_unit = integer_value_counts");

        for mut row in rows {
            row.push(("c2_holds", c2.to_string()));
            row.push(("order_only_share_max", format!("{best_order_share:.9}")));
            row.push(("c2_best_field", best_order_field.to_string()));
            run.record(&row);
        }
    });
}
