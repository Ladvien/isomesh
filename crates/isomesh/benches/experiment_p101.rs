//! **P-101 — can the duals be made equivariant by canonicalising the accumulation.**
//!
//! Ticket: R-101. Pre-registered before this harness existed.
//!
//! ```bash
//! cargo bench --bench experiment_p101
//! ```
//!
//! Writes `docs/experiments/p-101.csv`.
//!
//! # H, and the falsifier
//!
//! `M-372` says the duals accumulate crossings in an order that axis relabelling
//! permutes, and points at `M-177`. `P-101` asks whether accumulating in a
//! *relabelling-invariant* order — the crossings sorted by `(|value|, |offset|)`
//! with an invariant tie-break — takes `dual_contouring` from 3 of 16 rows at 48
//! to at least 12 (C1), for free (C2), with `manifold_dual_contouring` doing
//! worse (C3).
//!
//! Falsified by: C1 under 12, C2 by any geometric change, C3 by the two matching.
//!
//! # SHARE, recomputed from the source *before* this harness was written
//!
//! **The registered mechanism has zero share, and the reason is that A-016
//! already landed it.** Every reduction on the path from a cell's eight corner
//! samples to its dual vertex, read out of `crates/isomesh/src/`:
//!
//! | step | where | order-dependent on the axis labels? |
//! |---|---|---|
//! | corner sampling | the field | no |
//! | `is_inside`, `cube::edge_offset` | `cube.rs` | no — exactly antisymmetric (`P-61`'s four IEEE guarantees) |
//! | `cube::place`, the cell origin | `cube.rs`, the rule | no — componentwise |
//! | `Sdf::gradient` | the field | no — the rotated wrapper is `g·∇f(g⁻¹p)`, exact |
//! | **`vec3::length(gradient)`** | **`vec3.rs`** | **YES — `dot(a, a)` is `a0² + a1² + a2²` summed left to right in axis-index order** |
//! | `HermiteCell::centroid` | `hermite.rs` | no — twelve slots per edge, `sum_equivariant` (A-016) |
//! | `solve_with`'s `AᵀA` and `g` | `dual_contouring/solve.rs` | no — twelve slots per edge, `sum_equivariant` (A-016) |
//! | `determinant`, `adjugate`, `mul_vec` | same | no — `mul_equivariant` / `dot_equivariant` (M-24) |
//! | `apply_clamp` | `dual_contouring.rs` | no — componentwise `min`/`max` |
//!
//! So the accumulation the registration proposes to reorder is **already**
//! relabelling-invariant, and has been since A-016. The share a *different*
//! invariant order can move is 0 of the 1 remaining axis-index-dependent
//! reduction. Two numbers from the committed `p-61.csv` say the same thing
//! without reading any source:
//!
//! - `surface_nets` — the dual with **no normals at all** and the *same*
//!   `sum_equivariant` centroid — reaches `pure_permutation_exact` **6 of 6** on
//!   `sphere` at 33³, where `dual_contouring` reaches **2 of 6**. The
//!   accumulation is common to both, so the four-element gap cannot be the
//!   accumulation.
//! - `surface_nets` reaches `pure_sign_flip_exact` **1 of 8** on that same row.
//!   The sign-flip failure is therefore already present with the accumulation
//!   held constant, which is `M-177` and not an ordering question at all.
//!
//! **C1 is arithmetically unreachable before the run.** It is run anyway, to
//! produce the number, and with two extra arms that name what *is* in the way.
//!
//! # Five arms, one build, one run (M-281)
//!
//! | `accumulation_key` | crossing accumulation | normal normalisation |
//! |---|---|---|
//! | `shipped` | twelve slots per edge, magnitude-ordered | `vec3::length`, naive |
//! | `edge_slot` | same, **reimplemented here** | same |
//! | `abs_value_abs_offset` | **running sum in key order** | same |
//! | `edge_slot_equivariant_normal` | twelve slots, magnitude-ordered | **`dot_equivariant(g, g).sqrt()`** |
//! | `abs_value_abs_offset_equivariant_normal` | running sum in key order | equivariant |
//!
//! `edge_slot` is the **instrument check**: it is this file's own transcription
//! of `HermiteCell::from_corners`, `solve::solve_with` and `apply_clamp`, and it
//! is asserted **bit-identical** to the shipped extractor over every position,
//! normal and index on all 32 fixtures and all 48 golden rows. Without it, a
//! difference in a later arm could be a transcription error rather than the
//! change under test — which is `P-61`'s rule for a second copy of an instrument.
//!
//! The two `equivariant_normal` arms are not registered. They are the mechanism:
//! they change exactly the one reduction the table above names and nothing else.
//! If `dual_contouring`'s pure permutations go 2 → 6 there and the key arm does
//! not move, the obstruction is *located* rather than argued.
//!
//! # The key, and why it is spelled out here
//!
//! The registration says `(|value|, |offset|)`. An edge crossing has no single
//! "value", so this harness fixes the reading and says so:
//!
//! - `value = a + b`, the sum of the crossing's two corner samples. **Symmetric
//!   in `(a, b)`**, which it has to be: a reflection along the edge's own axis
//!   swaps which corner is low.
//! - `offset = cube::edge_offset(a, b) = ((a + b)/2)/(a − b)`, the crate's own
//!   signed offset from the edge midpoint. `d → −d` under that swap, so `|d|` is
//!   the invariant.
//!
//! `M-175` is why the tie-break is written before the first run rather than after
//! the first surprise. Ties are counted (`exact_key_ties`) and broken by, in
//! order: `|a − b|`; the crossing position's three component *magnitudes* sorted
//! descending; the normal's, likewise. Each of those is invariant under a signed
//! axis permutation, which is the group C1 is denominated in. A pair that ties on
//! all of them is counted separately as `unbreakable_key_ties` and ordered by
//! edge index — the one axis-dependent step left, reported rather than hidden.
//!
//! # The vacuity control
//!
//! `p-61.csv`'s three already-passing dual rows must still pass, and
//! `exact_key_ties` must be non-zero. Both are columns and both are `assert!`ed.
//! The shipped arm's `elements_vertex_exact`, `pure_permutation_exact` and
//! `pure_sign_flip_exact` are additionally asserted equal to the committed
//! `p-61.csv` values on all 32 dual rows before any new number is reported, and
//! the committed `golden_hashes.json` is asserted equal to the shipped arm's own
//! hash on all 48 golden dual rows before `hashes_moved` is believed.

// Exact comparisons on purpose: every clause here is stated in bits.
#![allow(clippy::float_cmp)]

mod common;

use std::cell::Cell;
use std::cmp::Ordering;
use std::collections::HashMap;
use std::path::PathBuf;
use std::rc::Rc;
use std::time::Instant;

use isomesh::dual::{CellVertices, VertexRule};
use isomesh::dual_contouring::solve::{LAMBDA, dot_equivariant};
use isomesh::dual_contouring::{CLAMP_EPSILON, DualContouring};
use isomesh::fields::ReferenceField;
use isomesh::manifold_dual_contouring::ManifoldDualContouring;
use isomesh::marching_cubes::ambiguity::joined_mask;
use isomesh::marching_cubes::table::{
    EDGE_CORNERS, EDGE_COUNT, NO_EDGE, edge_offset, is_inside, segment_links,
};
use isomesh::validate::{AccuracyConfig, accuracy, mesh_hash};
use isomesh::{MeshBuffer, RuntimeShape3, Sdf};

// ─── the crate's own reductions, transcribed ────────────────────────────────
//
// `crate::equivariant` is private, so these are copies rather than calls. They
// are copies of eleven lines whose behaviour the `edge_slot` arm asserts
// bit-for-bit against the shipped extractor, which is the only thing that makes
// a second copy of an instrument admissible here.

/// Whether `a` sums before `b`: smaller magnitude first, then smaller value.
#[inline]
fn precedes(a: f64, b: f64) -> bool {
    match a.abs().total_cmp(&b.abs()) {
        Ordering::Less => true,
        Ordering::Greater => false,
        Ordering::Equal => a.total_cmp(&b) == Ordering::Less,
    }
}

/// Insertion sort, ascending by magnitude then by signed value.
#[inline]
fn sort_by_magnitude<const N: usize>(t: &mut [f64; N]) {
    let mut i = 1;
    while i < N {
        let mut j = i;
        while j > 0 && precedes(t[j], t[j - 1]) {
            t.swap(j - 1, j);
            j -= 1;
        }
        i += 1;
    }
}

/// Sum smallest-magnitude-first.
#[inline]
fn sum_equivariant<const N: usize>(mut t: [f64; N]) -> f64 {
    sort_by_magnitude(&mut t);
    let mut acc = 0.0;
    for v in t {
        acc += v;
    }
    acc
}

/// Multiply smallest-magnitude-first.
#[inline]
fn mul_equivariant(mut t: [f64; 3]) -> f64 {
    sort_by_magnitude(&mut t);
    (t[0] * t[1]) * t[2]
}

/// `cube::corner_offset`, which is `pub(crate)`.
#[inline]
fn corner_offset(corner: u8) -> [u32; 3] {
    [
        u32::from(corner & 1),
        u32::from((corner >> 1) & 1),
        u32::from((corner >> 2) & 1),
    ]
}

/// `cube::place`, which is `pub(crate)`.
#[inline]
fn place(lo: f64, hi: f64, d: f64) -> f64 {
    (lo + hi) * 0.5 + (hi - lo) * d
}

/// `vec3::dot`, the naive axis-order sum. **This is the obstruction.**
#[inline]
fn naive_dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

#[inline]
fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

#[inline]
fn scale(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}

// ─── the symmetric solve, transcribed unchanged ─────────────────────────────
//
// `determinant`, `adjugate` and `mul_vec` are M-24's forms and are **not** what
// this experiment varies. They are here because `solve_with` is one function and
// only its crossing accumulation is under test.

/// A symmetric 3×3 matrix, stored as its six distinct entries.
#[derive(Clone, Copy)]
struct Symmetric3 {
    xx: f64,
    xy: f64,
    xz: f64,
    yy: f64,
    yz: f64,
    zz: f64,
}

impl Symmetric3 {
    #[inline]
    fn outer(n: [f64; 3]) -> [f64; 6] {
        [
            n[0] * n[0],
            n[0] * n[1],
            n[0] * n[2],
            n[1] * n[1],
            n[1] * n[2],
            n[2] * n[2],
        ]
    }

    #[inline]
    fn from_entries(e: [f64; 6]) -> Self {
        Self {
            xx: e[0],
            xy: e[1],
            xz: e[2],
            yy: e[3],
            yz: e[4],
            zz: e[5],
        }
    }

    #[inline]
    fn regularized(mut self, lambda: f64) -> Self {
        self.xx += lambda;
        self.yy += lambda;
        self.zz += lambda;
        self
    }

    #[inline]
    fn adjugate(self) -> Self {
        Self {
            xx: self.yy * self.zz - self.yz * self.yz,
            xy: self.xz * self.yz - self.xy * self.zz,
            xz: self.xy * self.yz - self.xz * self.yy,
            yy: self.xx * self.zz - self.xz * self.xz,
            yz: self.xy * self.xz - self.xx * self.yz,
            zz: self.xx * self.yy - self.xy * self.xy,
        }
    }

    #[inline]
    fn determinant(self) -> f64 {
        sum_equivariant([
            mul_equivariant([self.xx, self.yy, self.zz]),
            2.0 * mul_equivariant([self.xy, self.yz, self.xz]),
            -mul_equivariant([self.xx, self.yz, self.yz]),
            -mul_equivariant([self.yy, self.xz, self.xz]),
            -mul_equivariant([self.zz, self.xy, self.xy]),
        ])
    }

    #[inline]
    fn mul_vec(self, v: [f64; 3]) -> [f64; 3] {
        [
            dot_equivariant([self.xx, self.xy, self.xz], v),
            dot_equivariant([self.xy, self.yy, self.yz], v),
            dot_equivariant([self.xz, self.yz, self.zz], v),
        ]
    }
}

// ─── one cell's crossings, with the registered sort key on each ─────────────

/// One edge crossing, plus the three invariant scalars the key is built from.
#[derive(Clone, Copy)]
struct Crossing {
    /// World position — `cube::place` in the centred frame.
    position: [f64; 3],
    /// Unit surface normal, normalised by whichever rule this arm uses.
    normal: [f64; 3],
    /// `a + b`. Symmetric in the two corner samples, so the reflection along the
    /// edge's own axis leaves it alone.
    value: f64,
    /// `cube::edge_offset(a, b)`. Negated by that reflection, so `|offset|` is
    /// the invariant.
    offset: f64,
    /// `a − b`. First tie-break; `|span|` is invariant for the same reason.
    span: f64,
}

const EMPTY_CROSSING: Crossing = Crossing {
    position: [0.0; 3],
    normal: [0.0; 3],
    value: 0.0,
    offset: 0.0,
    span: 0.0,
};

/// The crossings on one cell's twelve edges, indexed by edge label.
struct CellCrossings {
    edge: [Crossing; EDGE_COUNT],
    mask: u16,
}

/// How a normal's length is computed.
#[derive(Clone, Copy, PartialEq, Eq)]
enum Norm {
    /// `vec3::length`, i.e. `(g0² + g1² + g2²).sqrt()` left to right. Shipped.
    Naive,
    /// `dot_equivariant(g, g).sqrt()`, magnitude-ordered. Not shipped.
    Equivariant,
}

/// Which order the crossings are accumulated in.
#[derive(Clone, Copy, PartialEq, Eq)]
enum Order {
    /// One slot per edge label, reduced by `sum_equivariant`. Shipped (A-016).
    EdgeSlots,
    /// A running sum over the crossings sorted by `(|value|, |offset|)`. The
    /// registered change.
    KeyRunning,
}

/// `HermiteCell::from_corners`, transcribed, with the normalisation switchable.
fn crossings_of<S: Sdf<Scalar = f64>>(
    sdf: &S,
    corner_values: &[f64; 8],
    cell_origin: [f64; 3],
    cell_size: f64,
    norm: Norm,
) -> CellCrossings {
    let mut out = CellCrossings {
        edge: [EMPTY_CROSSING; EDGE_COUNT],
        mask: 0,
    };
    for (edge, [lo, hi]) in EDGE_CORNERS.iter().copied().enumerate() {
        let (a, b) = (corner_values[lo as usize], corner_values[hi as usize]);
        if is_inside(a) == is_inside(b) {
            continue;
        }

        let d = edge_offset(a, b);
        let (lo_offset, hi_offset) = (corner_offset(lo), corner_offset(hi));
        let mut position = [0.0f64; 3];
        for (axis, slot) in position.iter_mut().enumerate() {
            let from = f64::from(lo_offset[axis]);
            let to = f64::from(hi_offset[axis]);
            *slot = cell_origin[axis] + cell_size * place(from, to, d);
        }

        let gradient = sdf.gradient(position);
        let length = match norm {
            Norm::Naive => naive_dot(gradient, gradient).sqrt(),
            Norm::Equivariant => dot_equivariant(gradient, gradient).sqrt(),
        };
        let normal = scale(gradient, length.recip());

        out.edge[edge] = Crossing {
            position,
            normal,
            value: a + b,
            offset: d,
            span: a - b,
        };
        out.mask |= 1 << edge;
    }
    out
}

/// The full sort key, as bit patterns of non-negative doubles.
///
/// Every component is `|x|` for a finite `x`, and for non-negative finite doubles
/// the IEEE bit pattern orders exactly as the value does — so derived `Ord` here
/// *is* the numeric order, with no float comparator in the sort.
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
struct Key {
    /// `|value|`, the registered primary.
    value: u64,
    /// `|offset|`, the registered secondary.
    offset: u64,
    /// `|span|`, the first tie-break.
    span: u64,
    /// The position's component magnitudes, descending. Invariant under a signed
    /// axis permutation, which is the whole group.
    position: [u64; 3],
    /// The normal's, likewise.
    normal: [u64; 3],
}

/// Component magnitudes, descending, as bits.
#[inline]
fn magnitudes(v: [f64; 3]) -> [u64; 3] {
    let mut m = [v[0].abs(), v[1].abs(), v[2].abs()];
    m.sort_unstable_by(|a, b| b.total_cmp(a));
    [m[0].to_bits(), m[1].to_bits(), m[2].to_bits()]
}

#[inline]
fn key_of(c: &Crossing) -> Key {
    Key {
        value: c.value.abs().to_bits(),
        offset: c.offset.abs().to_bits(),
        span: c.span.abs().to_bits(),
        position: magnitudes(c.position),
        normal: magnitudes(c.normal),
    }
}

/// What one arm's cells did, accumulated across a whole extraction.
#[derive(Default)]
struct Tally {
    cells: Cell<u64>,
    crossings: Cell<u64>,
    multi_crossing_cells: Cell<u64>,
    /// Adjacent pairs in the sorted list that tie on `(|value|, |offset|)` — the
    /// registered primary. **The tie-break's own vacuity control.**
    key_ties: Cell<u64>,
    /// Pairs that tie on the *whole* invariant key, and are therefore ordered by
    /// edge index. The residual axis dependence, counted.
    unbreakable_ties: Cell<u64>,
}

/// `solve::solve_with`, transcribed, with the crossing accumulation switchable.
///
/// `keep` restricts to a subset of the edges, which is what
/// `HermiteCell::restricted` does for the cycle rule.
fn solve_cell(
    cell: &CellCrossings,
    keep: u16,
    lambda: f64,
    order: Order,
    tally: &Tally,
) -> Option<[f64; 3]> {
    let mask = cell.mask & keep;
    let count = mask.count_ones() as usize;
    if count == 0 {
        return None;
    }

    // The key order, and the tie census, for **every** arm — so `exact_key_ties`
    // is a property of the fixture rather than of the arm that happened to use
    // it. Fixed-size, because this runs once per surface cell.
    let mut sorted = [(
        Key {
            value: 0,
            offset: 0,
            span: 0,
            position: [0; 3],
            normal: [0; 3],
        },
        0u8,
    ); EDGE_COUNT];
    let mut n = 0usize;
    for edge in 0..EDGE_COUNT {
        if mask & (1 << edge) == 0 {
            continue;
        }
        sorted[n] = (key_of(&cell.edge[edge]), edge as u8);
        n += 1;
    }
    let sorted = &mut sorted[..n];
    sorted.sort_unstable();
    for pair in sorted.windows(2) {
        let (a, b) = (&pair[0].0, &pair[1].0);
        if a.value == b.value && a.offset == b.offset {
            tally.key_ties.set(tally.key_ties.get() + 1);
            if a == b {
                tally.unbreakable_ties.set(tally.unbreakable_ties.get() + 1);
            }
        }
    }
    tally.cells.set(tally.cells.get() + 1);
    tally.crossings.set(tally.crossings.get() + count as u64);
    if count > 1 {
        tally
            .multi_crossing_cells
            .set(tally.multi_crossing_cells.get() + 1);
    }

    let inverse = (count as f64).recip();
    let centroid = match order {
        Order::EdgeSlots => {
            let mut axes = [[0.0f64; EDGE_COUNT]; 3];
            for edge in 0..EDGE_COUNT {
                if mask & (1 << edge) == 0 {
                    continue;
                }
                for (slot, value) in axes.iter_mut().zip(cell.edge[edge].position) {
                    slot[edge] = value;
                }
            }
            let sum = axes.map(sum_equivariant);
            [sum[0] * inverse, sum[1] * inverse, sum[2] * inverse]
        }
        Order::KeyRunning => {
            let mut sum = [0.0f64; 3];
            for &(_, edge) in sorted.iter() {
                for (slot, value) in sum.iter_mut().zip(cell.edge[edge as usize].position) {
                    *slot += value;
                }
            }
            [sum[0] * inverse, sum[1] * inverse, sum[2] * inverse]
        }
    };

    let (m, g) = match order {
        Order::EdgeSlots => {
            let mut m_terms = [[0.0f64; EDGE_COUNT]; 6];
            let mut g_terms = [[0.0f64; EDGE_COUNT]; 3];
            for edge in 0..EDGE_COUNT {
                if mask & (1 << edge) == 0 {
                    continue;
                }
                let c = &cell.edge[edge];
                let normal = c.normal;
                let d = dot_equivariant(normal, sub(c.position, centroid));
                for (slot, value) in m_terms.iter_mut().zip(Symmetric3::outer(normal)) {
                    slot[edge] = value;
                }
                for (slot, value) in g_terms
                    .iter_mut()
                    .zip([normal[0] * d, normal[1] * d, normal[2] * d])
                {
                    slot[edge] = value;
                }
            }
            (
                Symmetric3::from_entries(m_terms.map(sum_equivariant)),
                g_terms.map(sum_equivariant),
            )
        }
        Order::KeyRunning => {
            let mut m_sum = [0.0f64; 6];
            let mut g_sum = [0.0f64; 3];
            for &(_, edge) in sorted.iter() {
                let c = &cell.edge[edge as usize];
                let normal = c.normal;
                let d = dot_equivariant(normal, sub(c.position, centroid));
                for (slot, value) in m_sum.iter_mut().zip(Symmetric3::outer(normal)) {
                    *slot += value;
                }
                for (slot, value) in g_sum
                    .iter_mut()
                    .zip([normal[0] * d, normal[1] * d, normal[2] * d])
                {
                    *slot += value;
                }
            }
            (Symmetric3::from_entries(m_sum), g_sum)
        }
    };

    let a = m.regularized(lambda);
    let adj = a.adjugate();
    let det = a.determinant();
    let offset = scale(adj.mul_vec(g), det.recip());
    let x = [
        centroid[0] + offset[0],
        centroid[1] + offset[1],
        centroid[2] + offset[2],
    ];
    if x[0].is_finite() && x[1].is_finite() && x[2].is_finite() {
        Some(x)
    } else {
        None
    }
}

/// `dual_contouring::apply_clamp` with `Clamp::ToCell`, transcribed.
fn clamp_to_cell(x: [f64; 3], cell_origin: [f64; 3], cell_size: f64) -> [f64; 3] {
    let half = cell_size * 0.5;
    let inset = half * (1.0 - CLAMP_EPSILON);
    let mut out = x;
    for (axis, slot) in out.iter_mut().enumerate() {
        let centre = cell_origin[axis] + half;
        *slot = slot.clamp(centre - inset, centre + inset);
    }
    out
}

#[inline]
fn cell_origin_of(base: [u32; 3], origin: [f64; 3], cell_size: f64) -> [f64; 3] {
    [
        origin[0] + cell_size * f64::from(base[0]),
        origin[1] + cell_size * f64::from(base[1]),
        origin[2] + cell_size * f64::from(base[2]),
    ]
}

// ─── the two rules, as this bench's own `VertexRule`s ───────────────────────

/// `dual_contouring::Qef`, reimplemented so the accumulation can be swapped.
struct DcRule {
    order: Order,
    norm: Norm,
    tally: Rc<Tally>,
}

impl VertexRule<f64> for DcRule {
    fn place<S: Sdf<Scalar = f64>>(
        &self,
        sdf: &S,
        corner: &[f64; 8],
        base: [u32; 3],
        origin: [f64; 3],
        cell_size: f64,
        out: &mut CellVertices<f64>,
    ) {
        let cell_origin = cell_origin_of(base, origin, cell_size);
        let cell = crossings_of(sdf, corner, cell_origin, cell_size, self.norm);
        let Some(x) = solve_cell(&cell, u16::MAX, LAMBDA, self.order, &self.tally) else {
            return;
        };
        out.push_whole_cell(clamp_to_cell(x, cell_origin, cell_size));
    }
}

/// `manifold_dual_contouring::CycleQef`, reimplemented the same way.
///
/// `FaceAmbiguity::Separate` only, which is what `ManifoldDualContouring::new`
/// ships and what the committed golden hashes pin — so the `ambiguous` mask
/// handed to `joined_mask` is `0`, exactly as `CycleQef` computes it under that
/// setting.
struct CycleRule {
    order: Order,
    norm: Norm,
    tally: Rc<Tally>,
}

impl VertexRule<f64> for CycleRule {
    fn place<S: Sdf<Scalar = f64>>(
        &self,
        sdf: &S,
        corner: &[f64; 8],
        base: [u32; 3],
        origin: [f64; 3],
        cell_size: f64,
        out: &mut CellVertices<f64>,
    ) {
        let cell_origin = cell_origin_of(base, origin, cell_size);

        let mut case = 0u8;
        for (c, value) in corner.iter().enumerate() {
            if is_inside(*value) {
                case |= 1 << c;
            }
        }
        let next = segment_links(case, joined_mask(corner, 0));

        let cell = crossings_of(sdf, corner, cell_origin, cell_size, self.norm);

        let mut visited = 0u16;
        for start in 0..EDGE_COUNT as u8 {
            if next[start as usize] == NO_EDGE || visited & (1 << start) != 0 {
                continue;
            }
            let mut edges = 0u16;
            let mut current = start;
            while visited & (1 << current) == 0 {
                visited |= 1 << current;
                edges |= 1 << current;
                current = next[current as usize];
            }
            let Some(x) = solve_cell(&cell, edges, LAMBDA, self.order, &self.tally) else {
                continue;
            };
            out.push_component(clamp_to_cell(x, cell_origin, cell_size), edges);
        }
    }
}

// ─── one thing that can mesh a field ────────────────────────────────────────

/// Anything this harness can hand a field to.
///
/// A trait rather than a closure because the field type varies — the equivariance
/// arm wraps it in [`Rotated`] — and a closure cannot be generic.
trait Mesher {
    fn mesh<S: Sdf<Scalar = f64>>(
        &mut self,
        field: &S,
        shape: &RuntimeShape3,
        origin: [f64; 3],
        cell_size: f64,
        out: &mut MeshBuffer<f64>,
    );
}

impl<V: VertexRule<f64>> Mesher for DualContouring<f64, V> {
    fn mesh<S: Sdf<Scalar = f64>>(
        &mut self,
        field: &S,
        shape: &RuntimeShape3,
        origin: [f64; 3],
        cell_size: f64,
        out: &mut MeshBuffer<f64>,
    ) {
        out.reset();
        self.extract(field, shape, origin, cell_size, out)
            .expect("extraction");
    }
}

impl Mesher for ManifoldDualContouring<f64> {
    fn mesh<S: Sdf<Scalar = f64>>(
        &mut self,
        field: &S,
        shape: &RuntimeShape3,
        origin: [f64; 3],
        cell_size: f64,
        out: &mut MeshBuffer<f64>,
    ) {
        out.reset();
        self.extract(field, shape, origin, cell_size, out)
            .expect("extraction");
    }
}

// ─── the group, and the comparison keys (P-57's instrument, via P-61) ───────

/// The 6 axis permutations. Crossed with 8 sign patterns this is all 48.
const PERMS: [[usize; 3]; 6] = [
    [0, 1, 2],
    [0, 2, 1],
    [1, 0, 2],
    [1, 2, 0],
    [2, 0, 1],
    [2, 1, 0],
];

/// The order of the octahedral group.
const GROUP_ORDER: usize = 48;

/// Probe points for the inverse round-trip check. Both zeros deliberately.
const PROBES: [[f64; 3]; 4] = [
    [0.3, -1.7, 2.9],
    [0.0, -0.0, 1.5],
    [-2.25, 0.656_25, -0.187_5],
    [1.0, 1.0, 1.0],
];

/// One element of the octahedral group, as a signed axis permutation.
#[derive(Clone, Copy, PartialEq, Eq)]
struct Element {
    perm: [usize; 3],
    sign: [i8; 3],
}

impl Element {
    #[inline]
    fn apply(self, p: [f64; 3]) -> [f64; 3] {
        std::array::from_fn(|k| {
            let v = p[self.perm[k]];
            if self.sign[k] < 0 { -v } else { v }
        })
    }

    fn inverse(self) -> Self {
        let mut perm = [0usize; 3];
        let mut sign = [0i8; 3];
        for k in 0..3 {
            let j = self.perm[k];
            perm[j] = k;
            sign[j] = self.sign[k];
        }
        Self { perm, sign }
    }

    fn det(self) -> i32 {
        let mut m = [[0i32; 3]; 3];
        for k in 0..3 {
            m[k][self.perm[k]] = i32::from(self.sign[k]);
        }
        m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
            - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
            + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
    }

    fn label(self) -> String {
        let mut s = String::from("perm=");
        for k in 0..3 {
            s.push(char::from(b'0' + self.perm[k] as u8));
        }
        s.push_str(";sign=");
        for k in 0..3 {
            s.push(if self.sign[k] < 0 { '-' } else { '+' });
        }
        s
    }

    fn short(self) -> String {
        let mut s = String::new();
        for k in 0..3 {
            s.push(char::from(b'0' + self.perm[k] as u8));
        }
        s.push('/');
        for k in 0..3 {
            s.push(if self.sign[k] < 0 { '-' } else { '+' });
        }
        s
    }
}

fn group() -> Vec<Element> {
    let mut out = Vec::with_capacity(GROUP_ORDER);
    for perm in PERMS {
        for bits in 0..8u8 {
            let sign = std::array::from_fn(|k| if bits & (1 << k) == 0 { 1i8 } else { -1i8 });
            out.push(Element { perm, sign });
        }
    }
    out
}

fn verify_group(g: &[Element]) {
    assert_eq!(g.len(), GROUP_ORDER, "the octahedral group has 48 elements");
    assert!(
        g[0] == Element {
            perm: [0, 1, 2],
            sign: [1, 1, 1]
        },
        "element 0 must be the identity: it is the negative control"
    );
    for (i, a) in g.iter().enumerate() {
        for b in &g[i + 1..] {
            assert!(a != b, "duplicate group element at {i}: {}", a.label());
        }
    }
    let (mut rotations, mut reflections) = (0, 0);
    for e in g {
        match e.det() {
            1 => rotations += 1,
            -1 => reflections += 1,
            d => panic!("{} has determinant {d}, not +/-1", e.label()),
        }
        let inv = e.inverse();
        for p in PROBES {
            let round = inv.apply(e.apply(p));
            for k in 0..3 {
                assert_eq!(
                    round[k].to_bits(),
                    p[k].to_bits(),
                    "{} does not round-trip {p:?} bit-exactly",
                    e.label()
                );
            }
        }
    }
    assert_eq!(rotations, 24, "24 elements must have det = +1");
    assert_eq!(reflections, 24, "24 elements must have det = -1");
}

/// `g·f`, i.e. `(g·f)(p) = f(g⁻¹·p)`.
struct Rotated<'a, S> {
    field: &'a S,
    g: Element,
    g_inv: Element,
}

impl<S: Sdf<Scalar = f64>> Sdf for Rotated<'_, S> {
    type Scalar = f64;

    #[inline]
    fn sample(&self, p: [f64; 3]) -> f64 {
        self.field.sample(self.g_inv.apply(p))
    }

    #[inline]
    fn gradient(&self, p: [f64; 3]) -> [f64; 3] {
        self.g.apply(self.field.gradient(self.g_inv.apply(p)))
    }
}

const NEGATIVE_ZERO: u64 = 1u64 << 63;

/// `−0.0` folded onto `+0.0`, identically on both sides of every comparison.
///
/// Both fixtures are centred with an odd sample count, so `0.0` is a grid
/// coordinate and `−(0.0)` is `−0.0`; left raw, every sign-flipping element fails
/// on a disagreement about which *encoding* of zero was written. `✗39` settled
/// this and `P-61` kept it.
#[inline]
fn comparison_key(v: f64) -> u64 {
    let b = v.to_bits();
    if b == NEGATIVE_ZERO { 0 } else { b }
}

#[inline]
fn monotone(bits: u64) -> i128 {
    if bits & NEGATIVE_ZERO == 0 {
        i128::from(bits)
    } else {
        -i128::from(bits & !NEGATIVE_ZERO)
    }
}

#[inline]
fn ulp_distance(a: u64, b: u64) -> u128 {
    (monotone(a) - monotone(b)).unsigned_abs()
}

/// Multiset symmetric difference of two **sorted** key lists.
fn multiset_difference(a: &[[u64; 3]], b: &[[u64; 3]]) -> (Vec<[u64; 3]>, Vec<[u64; 3]>) {
    let (mut only_a, mut only_b) = (Vec::new(), Vec::new());
    let (mut i, mut j) = (0usize, 0usize);
    while i < a.len() && j < b.len() {
        match a[i].cmp(&b[j]) {
            Ordering::Less => {
                only_a.push(a[i]);
                i += 1;
            }
            Ordering::Greater => {
                only_b.push(b[j]);
                j += 1;
            }
            Ordering::Equal => {
                i += 1;
                j += 1;
            }
        }
    }
    only_a.extend_from_slice(&a[i..]);
    only_b.extend_from_slice(&b[j..]);
    (only_a, only_b)
}

fn vertex_keys(positions: &[[f64; 3]], g: Option<Element>) -> Vec<[u64; 3]> {
    let mut out: Vec<[u64; 3]> = positions
        .iter()
        .map(|p| {
            let q = match g {
                Some(e) => e.apply(*p),
                None => *p,
            };
            std::array::from_fn(|k| comparison_key(q[k]))
        })
        .collect();
    out.sort_unstable();
    out
}

fn triangle_keys(positions: &[[f64; 3]], indices: &[u32], g: Option<Element>) -> Vec<[[u64; 3]; 3]> {
    let mapped = |i: u32| -> [u64; 3] {
        let p = positions[i as usize];
        let q = match g {
            Some(e) => e.apply(p),
            None => p,
        };
        std::array::from_fn(|k| comparison_key(q[k]))
    };
    let mut out: Vec<[[u64; 3]; 3]> = Vec::with_capacity(indices.len() / 3);
    for tri in indices.as_chunks::<3>().0 {
        let c = [mapped(tri[0]), mapped(tri[1]), mapped(tri[2])];
        let rotations = [[c[0], c[1], c[2]], [c[1], c[2], c[0]], [c[2], c[0], c[1]]];
        out.push(rotations.into_iter().min().expect("three rotations"));
    }
    out.sort_unstable();
    out
}

// ─── fixtures ───────────────────────────────────────────────────────────────

#[derive(Clone, Copy)]
struct Fixture {
    samples: u32,
    cell_size: f64,
    origin: f64,
}

/// P-57's two fixtures, unchanged and reached through P-61.
///
/// The 25³ arm uses `3L/32` rather than the crate's `2L/24`, because `L/12` is
/// not dyadic and its grid does not mirror.
fn fixtures(l: f64) -> [Fixture; 2] {
    let h33 = l / 16.0;
    let h25 = 3.0 * l / 32.0;
    [
        Fixture {
            samples: 33,
            cell_size: h33,
            origin: -l,
        },
        Fixture {
            samples: 25,
            cell_size: h25,
            origin: -12.0 * h25,
        },
    ]
}

/// Is the grid bit-exactly closed under a sign flip?
///
/// Without this the relation would be falsified by the fixture rather than by the
/// extractor, which is P-57's own precondition.
fn grid_symmetric(fx: &Fixture) -> bool {
    let coords: Vec<f64> = (0..fx.samples)
        .map(|i| fx.origin + fx.cell_size * f64::from(i))
        .collect();
    let bits: Vec<u64> = coords.iter().map(|c| comparison_key(*c)).collect();
    coords.iter().all(|c| bits.contains(&comparison_key(-*c)))
}

/// The three resolutions `src/golden.rs` hashes every field at.
const GOLDEN_RESOLUTIONS: [u32; 3] = [17, 25, 33];

// ─── the committed artefacts, as the before-arms ────────────────────────────

/// One dual row of `docs/experiments/p-61.csv`, keyed `field/extractor/samples`.
struct P61Row {
    elements_vertex_exact: usize,
    pure_permutation_exact: usize,
    pure_sign_flip_exact: usize,
}

/// Read `p-61.csv`'s equivariance block. **Not optional**: it is the vacuity
/// control's source and the check that this harness is P-61's instrument.
fn p61_baseline() -> HashMap<String, P61Row> {
    let path = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../docs/experiments/p-61.csv");
    let text = std::fs::read_to_string(&path)
        .unwrap_or_else(|e| panic!("p-61.csv is the baseline and must be readable: {e}"));
    let mut lines = text.lines().filter(|l| !l.starts_with('#'));
    let header: Vec<&str> = lines
        .next()
        .expect("p-61.csv has a header")
        .split(',')
        .collect();
    let column = |name: &str| -> usize {
        header
            .iter()
            .position(|h| *h == name)
            .unwrap_or_else(|| panic!("p-61.csv has no `{name}` column"))
    };
    let (c_block, c_field, c_extractor, c_samples) = (
        column("block"),
        column("field"),
        column("extractor"),
        column("samples_per_axis"),
    );
    let (c_exact, c_perm, c_sign) = (
        column("elements_vertex_exact"),
        column("pure_permutation_exact"),
        column("pure_sign_flip_exact"),
    );

    let mut out = HashMap::new();
    for line in lines {
        let f: Vec<&str> = line.split(',').collect();
        if f.len() <= c_sign || f[c_block] != "equivariance" {
            continue;
        }
        if f[c_extractor] != "dual_contouring" && f[c_extractor] != "manifold_dual_contouring" {
            continue;
        }
        let parse = |s: &str, what: &str| -> usize {
            s.parse()
                .unwrap_or_else(|_| panic!("p-61.csv `{what}` is not a count: {s:?}"))
        };
        out.insert(
            format!("{}/{}/{}", f[c_field], f[c_extractor], f[c_samples]),
            P61Row {
                elements_vertex_exact: parse(f[c_exact], "elements_vertex_exact"),
                pure_permutation_exact: parse(f[c_perm], "pure_permutation_exact"),
                pure_sign_flip_exact: parse(f[c_sign], "pure_sign_flip_exact"),
            },
        );
    }
    assert_eq!(
        out.len(),
        32,
        "p-61.csv must carry 32 dual equivariance rows; found {}",
        out.len()
    );
    out
}

/// Pull one value out of a `golden_hashes.json` line.
///
/// A hand-rolled scanner rather than a JSON parser, for the reason `golden.rs`
/// gives: the grammar is one line, fixed key order, no nesting and no escapes.
fn json_field<'a>(line: &'a str, key: &str) -> Option<&'a str> {
    let needle = format!("\"{key}\":");
    let start = line.find(&needle)? + needle.len();
    let rest = line[start..].trim_start();
    if let Some(q) = rest.strip_prefix('"') {
        q.find('"').map(|end| &q[..end])
    } else {
        let end = rest.find([',', '}']).unwrap_or(rest.len());
        Some(rest[..end].trim())
    }
}

/// The committed golden hashes for the two duals, keyed `field/extractor/samples`.
fn golden_baseline() -> HashMap<String, u64> {
    let path = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("golden_hashes.json");
    let text = std::fs::read_to_string(&path)
        .unwrap_or_else(|e| panic!("golden_hashes.json is C2's baseline: {e}"));
    let mut out = HashMap::new();
    for line in text.lines() {
        let Some(algorithm) = json_field(line, "algorithm") else {
            continue;
        };
        if algorithm != "dual_contouring" && algorithm != "manifold_dual_contouring" {
            continue;
        }
        let field = json_field(line, "field").expect("a golden row names its field");
        let samples = json_field(line, "samples").expect("a golden row names its resolution");
        let hash = json_field(line, "hash").expect("a golden row carries a hash");
        out.insert(
            format!("{field}/{algorithm}/{samples}"),
            u64::from_str_radix(hash, 16).expect("a golden hash is hex"),
        );
    }
    assert_eq!(
        out.len(),
        48,
        "golden_hashes.json must carry 48 dual rows (8 fields x 3 resolutions x 2); found {}",
        out.len()
    );
    out
}

// ─── the equivariance measurement ───────────────────────────────────────────

struct Measured {
    vertices: usize,
    triangles: usize,
    vertex_exact: usize,
    triangle_exact: usize,
    pure_permutation_exact: usize,
    pure_sign_flip_exact: usize,
    worst_differing_vertices: usize,
    first_failing_element: String,
    worst_component_ulp: u128,
    vertex_failing_labels: String,
    wall_ms: u128,
}

fn measure<S: Sdf<Scalar = f64>, M: Mesher>(
    field: &S,
    mesher: &mut M,
    fx: &Fixture,
    elements: &[Element],
    reference: &mut MeshBuffer<f64>,
    rotated: &mut MeshBuffer<f64>,
) -> Measured {
    let shape = RuntimeShape3::new([fx.samples; 3]).expect("fixture grid fits u32");
    let origin = [fx.origin; 3];
    let started = Instant::now();

    mesher.mesh(field, &shape, origin, fx.cell_size, reference);

    let mut m = Measured {
        vertices: reference.positions.len(),
        triangles: reference.triangle_count(),
        vertex_exact: 0,
        triangle_exact: 0,
        pure_permutation_exact: 0,
        pure_sign_flip_exact: 0,
        worst_differing_vertices: 0,
        first_failing_element: String::from("none"),
        worst_component_ulp: 0,
        vertex_failing_labels: String::new(),
        wall_ms: 0,
    };
    let mut failing: Vec<String> = Vec::new();
    let mut any_failed = false;

    for (index, &g) in elements.iter().enumerate() {
        let wrapped = Rotated {
            field,
            g,
            g_inv: g.inverse(),
        };
        mesher.mesh(&wrapped, &shape, origin, fx.cell_size, rotated);

        let got = vertex_keys(&rotated.positions, None);
        let want = vertex_keys(&reference.positions, Some(g));
        let vertex_ok = got == want;

        let got_tri = triangle_keys(&rotated.positions, &rotated.indices, None);
        let want_tri = triangle_keys(&reference.positions, &reference.indices, Some(g));
        let triangle_ok = got_tri == want_tri;

        if vertex_ok {
            m.vertex_exact += 1;
            if g.sign == [1, 1, 1] {
                m.pure_permutation_exact += 1;
            }
            if g.perm == [0, 1, 2] {
                m.pure_sign_flip_exact += 1;
            }
        } else {
            let (only_got, only_want) = multiset_difference(&got, &want);
            m.worst_differing_vertices = m.worst_differing_vertices.max(only_got.len());
            for (a, b) in only_got.iter().zip(only_want.iter()) {
                for k in 0..3 {
                    if a[k] != b[k] {
                        m.worst_component_ulp = m.worst_component_ulp.max(ulp_distance(a[k], b[k]));
                    }
                }
            }
            if !any_failed {
                m.first_failing_element = g.label();
                any_failed = true;
            }
            failing.push(g.short());
        }
        if triangle_ok {
            m.triangle_exact += 1;
        }

        // **The control**: element 0 is the identity and must reproduce the
        // reference exactly, or the extractor is not deterministic and the row
        // measures nothing.
        assert!(
            index != 0 || (vertex_ok && triangle_ok),
            "the identity element is not exact: the extractor is not \
             deterministic, so this row measures nothing"
        );
    }

    m.vertex_failing_labels = failing.join("|");
    m.wall_ms = started.elapsed().as_millis();
    m
}

/// Are two meshes bit-identical in every position, normal and index?
fn bit_identical(a: &MeshBuffer<f64>, b: &MeshBuffer<f64>) -> bool {
    a.indices == b.indices
        && a.positions.len() == b.positions.len()
        && a.normals.len() == b.normals.len()
        && a.positions
            .iter()
            .zip(&b.positions)
            .all(|(p, q)| (0..3).all(|k| p[k].to_bits() == q[k].to_bits()))
        && a.normals
            .iter()
            .zip(&b.normals)
            .all(|(p, q)| (0..3).all(|k| p[k].to_bits() == q[k].to_bits()))
}

/// How many vertex positions differ in bits, index for index.
///
/// A pure reordering of a sum cannot move a sign classification, so the two arms
/// walk the same cells in the same order and the vertex arrays are index-aligned.
/// `counts_identical` is recorded beside this so a reader can see that the
/// alignment held.
fn positions_moved(a: &MeshBuffer<f64>, b: &MeshBuffer<f64>) -> usize {
    a.positions
        .iter()
        .zip(&b.positions)
        .filter(|(p, q)| (0..3).any(|k| p[k].to_bits() != q[k].to_bits()))
        .count()
}

fn hausdorff<S: Sdf<Scalar = f64>>(
    mesh: &MeshBuffer<f64>,
    field: &S,
    samples: u32,
    origin: [f64; 3],
    cell_size: f64,
) -> f64 {
    let shape = RuntimeShape3::new([samples; 3]).expect("grid fits u32");
    let cfg = AccuracyConfig::from_cell_size(cell_size).expect("valid cell size");
    accuracy(&mesh.positions, &mesh.indices, field, &shape, origin, &cfg)
        .expect("accuracy")
        .symmetric_hausdorff()
}

fn copy_mesh(from: &MeshBuffer<f64>, to: &mut MeshBuffer<f64>) {
    to.reset();
    to.positions.extend_from_slice(&from.positions);
    to.normals.extend_from_slice(&from.normals);
    to.indices.extend_from_slice(&from.indices);
}

// ─── the arms ───────────────────────────────────────────────────────────────

/// One of the four bench-local configurations, plus its CSV labels.
///
/// The three column names are spelled out per arm rather than formatted at
/// runtime because `Run::record` takes `&'static str` keys, and manufacturing a
/// `'static` from a `String` is not a thing this repository does.
struct ArmSpec {
    label: &'static str,
    order: Order,
    norm: Norm,
    /// Whether this arm must be bit-identical to the shipped extractor.
    is_replica: bool,
    /// Column carrying this arm's `dual_contouring` rows-at-48 on every row.
    col_dc: &'static str,
    /// Column carrying its `manifold_dual_contouring` rows-at-48.
    col_mdc: &'static str,
    /// Column carrying how many of the 48 golden dual hashes it moves.
    col_golden: &'static str,
    /// Column carrying how many of the 32 dual rows reach `pure_permutation_exact
    /// = 6`, i.e. full equivariance under the six pure axis relabellings. **This
    /// is the column the mechanism lives in**: 48 of 48 is out of reach for a
    /// magnitude-ordered sum on the sign-flip half (M-177), so a change that
    /// fixes *relabelling* shows up here and nowhere else.
    col_perm: &'static str,
}

const ARMS: [ArmSpec; 4] = [
    ArmSpec {
        label: "edge_slot",
        order: Order::EdgeSlots,
        norm: Norm::Naive,
        is_replica: true,
        col_dc: "arm_edge_slot_rows_at_48_dual_contouring",
        col_mdc: "arm_edge_slot_rows_at_48_manifold",
        col_golden: "arm_edge_slot_golden_hashes_moved",
        col_perm: "arm_edge_slot_dual_rows_at_perm_6",
    },
    ArmSpec {
        label: "abs_value_abs_offset",
        order: Order::KeyRunning,
        norm: Norm::Naive,
        is_replica: false,
        col_dc: "arm_abs_value_abs_offset_rows_at_48_dual_contouring",
        col_mdc: "arm_abs_value_abs_offset_rows_at_48_manifold",
        col_golden: "arm_abs_value_abs_offset_golden_hashes_moved",
        col_perm: "arm_abs_value_abs_offset_dual_rows_at_perm_6",
    },
    ArmSpec {
        label: "edge_slot_equivariant_normal",
        order: Order::EdgeSlots,
        norm: Norm::Equivariant,
        is_replica: false,
        col_dc: "arm_edge_slot_equivariant_normal_rows_at_48_dual_contouring",
        col_mdc: "arm_edge_slot_equivariant_normal_rows_at_48_manifold",
        col_golden: "arm_edge_slot_equivariant_normal_golden_hashes_moved",
        col_perm: "arm_edge_slot_equivariant_normal_dual_rows_at_perm_6",
    },
    ArmSpec {
        label: "abs_value_abs_offset_equivariant_normal",
        order: Order::KeyRunning,
        norm: Norm::Equivariant,
        is_replica: false,
        col_dc: "arm_abs_value_abs_offset_equivariant_normal_rows_at_48_dual_contouring",
        col_mdc: "arm_abs_value_abs_offset_equivariant_normal_rows_at_48_manifold",
        col_golden: "arm_abs_value_abs_offset_equivariant_normal_golden_hashes_moved",
        col_perm: "arm_abs_value_abs_offset_equivariant_normal_dual_rows_at_perm_6",
    },
];

/// The registered arm — the one C1, C2 and C3 are scored on.
const REGISTERED_ARM: &str = "abs_value_abs_offset";

/// The two dual extractors, by the name `p-61.csv` and `golden_hashes.json` use.
const EXTRACTORS: [&str; 2] = ["dual_contouring", "manifold_dual_contouring"];

/// Every row, before the aggregates are known.
type Row = Vec<(&'static str, String)>;

const NA: &str = "";

fn main() {
    if !std::env::args().any(|arg| arg == "--bench") {
        return;
    }
    let prereg = isomesh::experiment!("P-101");
    common::experiment::run(prereg, |run| {
        let elements = group();
        verify_group(&elements);
        let p61 = p61_baseline();
        let golden = golden_baseline();

        println!(
            "-- SHARE, recomputed before the run --\n\
             -- A-016 already made the crossing accumulation relabelling-invariant\n\
             -- (twelve slots per edge, sum_equivariant). The only axis-index-ordered\n\
             -- reduction left on the dual vertex path is vec3::length's naive dot, and\n\
             -- the registered change cannot reach it. Share = 0 of 1.\n\
             -- C1 is UNREACHABLE before the numbers below are read; it is run anyway.\n"
        );

        let mut rows: Vec<Row> = Vec::new();
        let mut golden_rows: Vec<Row> = Vec::new();

        let mut reference = MeshBuffer::<f64>::new();
        let mut rotated = MeshBuffer::<f64>::new();
        let mut shipped_mesh = MeshBuffer::<f64>::new();
        let mut arm_mesh = MeshBuffer::<f64>::new();

        let mut rows_at_48: HashMap<(&str, &str), usize> = HashMap::new();
        // Rows reaching full equivariance under the six *pure axis relabellings*.
        // The sign-flip half is out of reach for a magnitude-ordered sum
        // (M-177), so this is the counter a relabelling fix can actually move.
        let mut rows_at_perm_6: HashMap<&str, usize> = HashMap::new();
        let mut baseline_rows_at_48: HashMap<&str, usize> = HashMap::new();
        let mut hashes_moved: HashMap<&str, usize> = HashMap::new();
        let mut hashes_moved_expected: HashMap<&str, usize> = HashMap::new();
        let mut worst_hausdorff_delta: f64 = 0.0;
        let mut total_key_ties: u64 = 0;
        let mut total_unbreakable_ties: u64 = 0;
        let mut replica_bit_identical = true;
        let mut baseline_matches_p61 = true;
        let mut golden_fixture_matches_shipped = true;
        let mut golden_counts_changed = 0usize;
        let mut topology_identical = true;
        let mut vacuity_rows_still_48 = 0usize;

        // ── block: equivariance ─────────────────────────────────────────────
        println!("-- equivariance: P-57's 48-element instrument, five arms, one build --");
        println!(
            "{:<15} {:>4} {:<25} {:<40} {:>4} {:>5} {:>3} {:>3} {:>6} {:>6}",
            "field", "n", "extractor", "accumulation_key", "vex", "base", "p6", "f8", "ties", "ms"
        );

        isomesh::for_each_reference_field!(f64, |field_name, field| {
            let (lo, hi) = field.domain();
            let l = hi[0];
            for k in 0..3 {
                assert_eq!(
                    lo[k].to_bits(),
                    (-hi[k]).to_bits(),
                    "{field_name}: domain is not the symmetric cube the fixtures assume"
                );
            }

            for fx in fixtures(l) {
                assert!(
                    grid_symmetric(&fx),
                    "{field_name} at {}: the grid is not bit-exactly closed under a sign \
                     flip, so the relation would be falsified by the fixture",
                    fx.samples
                );
                let origin = [fx.origin; 3];

                for extractor in EXTRACTORS {
                    let shipped = if extractor == "dual_contouring" {
                        let mut m = DualContouring::<f64>::new();
                        measure(&field, &mut m, &fx, &elements, &mut reference, &mut rotated)
                    } else {
                        let mut m = ManifoldDualContouring::<f64>::new();
                        measure(&field, &mut m, &fx, &elements, &mut reference, &mut rotated)
                    };
                    copy_mesh(&reference, &mut shipped_mesh);

                    let k = format!("{field_name}/{extractor}/{}", fx.samples);
                    let b = p61
                        .get(&k)
                        .unwrap_or_else(|| panic!("p-61.csv has no dual row {k}"));
                    // **The instrument check.** This harness re-implements P-57's
                    // group, fixtures and comparison keys; if its shipped arm
                    // disagrees with the committed p-61.csv on the same
                    // configuration, it is not P-61's instrument and nothing
                    // below means anything.
                    let matched = shipped.vertex_exact == b.elements_vertex_exact
                        && shipped.pure_permutation_exact == b.pure_permutation_exact
                        && shipped.pure_sign_flip_exact == b.pure_sign_flip_exact;
                    assert!(
                        matched,
                        "{k}: this harness disagrees with p-61.csv about the SHIPPED \
                         extractor (vertex_exact {} vs {}, pure_perm {} vs {}, pure_sign \
                         {} vs {}) -- the instrument drifted",
                        shipped.vertex_exact,
                        b.elements_vertex_exact,
                        shipped.pure_permutation_exact,
                        b.pure_permutation_exact,
                        shipped.pure_sign_flip_exact,
                        b.pure_sign_flip_exact
                    );
                    baseline_matches_p61 &= matched;
                    if shipped.vertex_exact == GROUP_ORDER {
                        *baseline_rows_at_48.entry(extractor).or_default() += 1;
                    }

                    let shipped_hausdorff =
                        hausdorff(&shipped_mesh, &field, fx.samples, origin, fx.cell_size);
                    let shipped_hash = mesh_hash(&shipped_mesh);

                    for arm in &ARMS {
                        let tally = Rc::new(Tally::default());
                        let armed = if extractor == "dual_contouring" {
                            let mut m = DualContouring::<f64, DcRule>::with_rule(DcRule {
                                order: arm.order,
                                norm: arm.norm,
                                tally: Rc::clone(&tally),
                            });
                            measure(&field, &mut m, &fx, &elements, &mut arm_mesh, &mut rotated)
                        } else {
                            let mut m = DualContouring::<f64, CycleRule>::with_rule(CycleRule {
                                order: arm.order,
                                norm: arm.norm,
                                tally: Rc::clone(&tally),
                            });
                            measure(&field, &mut m, &fx, &elements, &mut arm_mesh, &mut rotated)
                        };
                        // `measure` leaves the arm's own reference mesh in
                        // `arm_mesh` and writes every rotated run into `rotated`,
                        // so this comparison is arm-vs-shipped on one grid.
                        let identical = bit_identical(&arm_mesh, &shipped_mesh);
                        if arm.is_replica {
                            assert!(
                                identical,
                                "{k} / {}: the transcription is not the shipped arithmetic, \
                                 so every other arm's difference is unattributable",
                                arm.label
                            );
                            replica_bit_identical &= identical;
                        }

                        let moved = positions_moved(&arm_mesh, &shipped_mesh);
                        let counts_same = arm_mesh.positions.len() == shipped_mesh.positions.len();
                        let indices_same = arm_mesh.indices == shipped_mesh.indices;
                        if !(counts_same && indices_same) {
                            topology_identical = false;
                        }

                        let arm_hausdorff =
                            hausdorff(&arm_mesh, &field, fx.samples, origin, fx.cell_size);
                        let delta = (arm_hausdorff - shipped_hausdorff).abs();
                        worst_hausdorff_delta = worst_hausdorff_delta.max(delta);
                        let arm_hash = mesh_hash(&arm_mesh);

                        if armed.vertex_exact == GROUP_ORDER {
                            *rows_at_48.entry((arm.label, extractor)).or_default() += 1;
                        }
                        if armed.pure_permutation_exact == PERMS.len() {
                            *rows_at_perm_6.entry(arm.label).or_default() += 1;
                        }
                        // The vacuity control: the three p-61 rows that already
                        // pass must still pass under the registered arm.
                        if arm.label == REGISTERED_ARM
                            && b.elements_vertex_exact == GROUP_ORDER
                            && armed.vertex_exact == GROUP_ORDER
                        {
                            vacuity_rows_still_48 += 1;
                        }
                        total_key_ties += tally.key_ties.get();
                        total_unbreakable_ties += tally.unbreakable_ties.get();

                        println!(
                            "{:<15} {:>4} {:<25} {:<40} {:>4} {:>5} {:>3} {:>3} {:>6} {:>6}",
                            field_name,
                            fx.samples,
                            extractor,
                            arm.label,
                            armed.vertex_exact,
                            shipped.vertex_exact,
                            armed.pure_permutation_exact,
                            armed.pure_sign_flip_exact,
                            tally.key_ties.get(),
                            armed.wall_ms
                        );

                        rows.push(vec![
                            ("block", "equivariance".to_string()),
                            ("field", field_name.to_string()),
                            ("resolution", fx.samples.to_string()),
                            ("extractor", extractor.to_string()),
                            ("accumulation_key", arm.label.to_string()),
                            ("cell_size", format!("{:.9}", fx.cell_size)),
                            ("elements_tested", GROUP_ORDER.to_string()),
                            ("elements_vertex_exact", armed.vertex_exact.to_string()),
                            (
                                "elements_vertex_exact_baseline",
                                shipped.vertex_exact.to_string(),
                            ),
                            (
                                "elements_vertex_exact_p61",
                                b.elements_vertex_exact.to_string(),
                            ),
                            ("elements_triangle_exact", armed.triangle_exact.to_string()),
                            (
                                "pure_permutation_exact",
                                armed.pure_permutation_exact.to_string(),
                            ),
                            (
                                "pure_permutation_exact_baseline",
                                shipped.pure_permutation_exact.to_string(),
                            ),
                            (
                                "pure_sign_flip_exact",
                                armed.pure_sign_flip_exact.to_string(),
                            ),
                            (
                                "pure_sign_flip_exact_baseline",
                                shipped.pure_sign_flip_exact.to_string(),
                            ),
                            ("worst_component_ulp", armed.worst_component_ulp.to_string()),
                            (
                                "worst_differing_vertices",
                                armed.worst_differing_vertices.to_string(),
                            ),
                            ("first_failing_element", armed.first_failing_element.clone()),
                            ("vertex_failing_labels", armed.vertex_failing_labels.clone()),
                            ("exact_key_ties", tally.key_ties.get().to_string()),
                            (
                                "unbreakable_key_ties",
                                tally.unbreakable_ties.get().to_string(),
                            ),
                            ("cells_solved", tally.cells.get().to_string()),
                            ("crossings_visited", tally.crossings.get().to_string()),
                            (
                                "multi_crossing_cells",
                                tally.multi_crossing_cells.get().to_string(),
                            ),
                            ("vertices", armed.vertices.to_string()),
                            ("triangles", armed.triangles.to_string()),
                            ("vertices_baseline", shipped.vertices.to_string()),
                            ("positions_moved", moved.to_string()),
                            ("counts_identical", counts_same.to_string()),
                            ("indices_identical", indices_same.to_string()),
                            ("bit_identical_to_shipped", identical.to_string()),
                            ("hausdorff_shipped", format!("{shipped_hausdorff:.12e}")),
                            ("hausdorff_arm", format!("{arm_hausdorff:.12e}")),
                            ("hausdorff_delta", format!("{delta:.12e}")),
                            ("mesh_hash_shipped", format!("{shipped_hash:016x}")),
                            ("mesh_hash_arm", format!("{arm_hash:016x}")),
                            ("mesh_hash_moved", (arm_hash != shipped_hash).to_string()),
                            ("wall_ms", armed.wall_ms.to_string()),
                        ]);
                    }
                }
            }
        });

        // ── block: golden ───────────────────────────────────────────────────
        //
        // C2's cost, over the configuration `src/golden.rs` actually hashes:
        // `origin = lo`, `cell_size = (hi − lo)/(samples − 1)`, three
        // resolutions. The 33³ equivariance fixture coincides with the 33-sample
        // golden row; 17 and 25 do not, which is why this block exists rather
        // than being read off the rows above.
        println!("\n-- golden: would a committed hash move, and by what --");
        println!(
            "{:<15} {:>4} {:<25} {:<40} {:>7} {:>8}",
            "field", "n", "extractor", "accumulation_key", "moved", "vertices"
        );

        isomesh::for_each_reference_field!(f64, |field_name, field| {
            let (lo, hi) = field.domain();
            for samples in GOLDEN_RESOLUTIONS {
                let cell_size = (hi[0] - lo[0]) / f64::from(samples - 1);
                let shape = RuntimeShape3::new([samples; 3]).expect("golden grid fits u32");

                for extractor in EXTRACTORS {
                    if extractor == "dual_contouring" {
                        let mut m = DualContouring::<f64>::new();
                        m.mesh(&field, &shape, lo, cell_size, &mut shipped_mesh);
                    } else {
                        let mut m = ManifoldDualContouring::<f64>::new();
                        m.mesh(&field, &shape, lo, cell_size, &mut shipped_mesh);
                    }
                    let shipped_hash = mesh_hash(&shipped_mesh);
                    let g_key = format!("{field_name}/{extractor}/{samples}");
                    let committed = *golden
                        .get(&g_key)
                        .unwrap_or_else(|| panic!("golden_hashes.json has no {g_key}"));
                    // Without this, `hashes_moved` would be measured against a
                    // stale fixture and could not be read as a cost.
                    if committed != shipped_hash {
                        golden_fixture_matches_shipped = false;
                    }
                    let shipped_hausdorff =
                        hausdorff(&shipped_mesh, &field, samples, lo, cell_size);

                    for arm in &ARMS {
                        let tally = Rc::new(Tally::default());
                        if extractor == "dual_contouring" {
                            let mut m = DualContouring::<f64, DcRule>::with_rule(DcRule {
                                order: arm.order,
                                norm: arm.norm,
                                tally: Rc::clone(&tally),
                            });
                            m.mesh(&field, &shape, lo, cell_size, &mut arm_mesh);
                        } else {
                            let mut m = DualContouring::<f64, CycleRule>::with_rule(CycleRule {
                                order: arm.order,
                                norm: arm.norm,
                                tally: Rc::clone(&tally),
                            });
                            m.mesh(&field, &shape, lo, cell_size, &mut arm_mesh);
                        }
                        let arm_hash = mesh_hash(&arm_mesh);
                        let identical = bit_identical(&arm_mesh, &shipped_mesh);
                        if arm.is_replica {
                            assert!(
                                identical,
                                "{g_key} / {}: the transcription is not the shipped \
                                 arithmetic on the golden configuration either",
                                arm.label
                            );
                        }
                        let moved = positions_moved(&arm_mesh, &shipped_mesh);
                        let counts_same =
                            arm_mesh.positions.len() == shipped_mesh.positions.len();
                        let indices_same = arm_mesh.indices == shipped_mesh.indices;
                        if !(counts_same && indices_same) {
                            topology_identical = false;
                            golden_counts_changed += 1;
                        }
                        let hash_moved = arm_hash != shipped_hash;
                        if hash_moved {
                            *hashes_moved.entry(arm.label).or_default() += 1;
                        }
                        if moved > 0 {
                            *hashes_moved_expected.entry(arm.label).or_default() += 1;
                        }
                        let arm_hausdorff = hausdorff(&arm_mesh, &field, samples, lo, cell_size);
                        let delta = (arm_hausdorff - shipped_hausdorff).abs();
                        worst_hausdorff_delta = worst_hausdorff_delta.max(delta);

                        println!(
                            "{:<15} {:>4} {:<25} {:<40} {:>7} {:>8}",
                            field_name,
                            samples,
                            extractor,
                            arm.label,
                            hash_moved,
                            arm_mesh.positions.len()
                        );

                        golden_rows.push(vec![
                            ("block", "golden".to_string()),
                            ("field", field_name.to_string()),
                            ("resolution", samples.to_string()),
                            ("extractor", extractor.to_string()),
                            ("accumulation_key", arm.label.to_string()),
                            ("cell_size", format!("{cell_size:.9}")),
                            ("elements_vertex_exact", NA.to_string()),
                            ("elements_vertex_exact_baseline", NA.to_string()),
                            ("worst_component_ulp", NA.to_string()),
                            ("exact_key_ties", tally.key_ties.get().to_string()),
                            (
                                "unbreakable_key_ties",
                                tally.unbreakable_ties.get().to_string(),
                            ),
                            ("cells_solved", tally.cells.get().to_string()),
                            ("crossings_visited", tally.crossings.get().to_string()),
                            (
                                "multi_crossing_cells",
                                tally.multi_crossing_cells.get().to_string(),
                            ),
                            ("vertices", arm_mesh.positions.len().to_string()),
                            ("triangles", arm_mesh.triangle_count().to_string()),
                            ("vertices_baseline", shipped_mesh.positions.len().to_string()),
                            ("positions_moved", moved.to_string()),
                            ("counts_identical", counts_same.to_string()),
                            ("indices_identical", indices_same.to_string()),
                            ("bit_identical_to_shipped", identical.to_string()),
                            ("hausdorff_shipped", format!("{shipped_hausdorff:.12e}")),
                            ("hausdorff_arm", format!("{arm_hausdorff:.12e}")),
                            ("hausdorff_delta", format!("{delta:.12e}")),
                            ("golden_hash_committed", format!("{committed:016x}")),
                            ("mesh_hash_shipped", format!("{shipped_hash:016x}")),
                            ("mesh_hash_arm", format!("{arm_hash:016x}")),
                            ("mesh_hash_moved", hash_moved.to_string()),
                        ]);
                    }
                }
            }
        });

        // ── the aggregates, and the clause verdicts ─────────────────────────
        let dc_at_48 = *rows_at_48
            .get(&(REGISTERED_ARM, "dual_contouring"))
            .unwrap_or(&0);
        let mdc_at_48 = *rows_at_48
            .get(&(REGISTERED_ARM, "manifold_dual_contouring"))
            .unwrap_or(&0);
        let dc_base = *baseline_rows_at_48.get("dual_contouring").unwrap_or(&0);
        let mdc_base = *baseline_rows_at_48
            .get("manifold_dual_contouring")
            .unwrap_or(&0);
        let moved = *hashes_moved.get(REGISTERED_ARM).unwrap_or(&0);
        let moved_expected = *hashes_moved_expected.get(REGISTERED_ARM).unwrap_or(&0);

        let c1_holds = dc_at_48 >= 12;
        let c2_holds =
            topology_identical && moved == moved_expected && worst_hausdorff_delta <= 1e-12;
        let c3_holds = mdc_at_48 < dc_at_48;

        // **The vacuity control, asserted rather than merely reported.** A run in
        // which the three already-passing rows regressed, or in which the
        // tie-break was never reached, measures the fixture and not the change
        // (M-44, M-175).
        assert_eq!(
            dc_base, 3,
            "the shipped dual_contouring must reproduce p-61.csv's 3 of 16 rows at 48"
        );
        assert_eq!(
            mdc_base, 3,
            "the shipped manifold_dual_contouring must reproduce p-61.csv's 3 of 16"
        );
        assert!(
            total_key_ties > 0,
            "exact_key_ties is zero: the tie-break was never exercised, which is exactly \
             the M-175 failure this column exists to prevent"
        );
        assert!(
            replica_bit_identical,
            "the edge_slot arm must be bit-identical to the shipped extractor"
        );
        assert!(
            baseline_matches_p61,
            "the shipped arm must reproduce p-61.csv row for row"
        );
        assert!(
            golden_fixture_matches_shipped,
            "the committed golden hashes must match the shipped extractor, or \
             `hashes_moved` is measured against a stale fixture"
        );

        println!("\n-- aggregates --");
        println!(
            "  C1  dual_contouring rows at 48, registered arm: {dc_at_48} of 16 \
             (shipped baseline {dc_base})"
        );
        println!(
            "  C3  manifold rows at 48, registered arm:        {mdc_at_48} of 16 \
             (shipped baseline {mdc_base})"
        );
        println!("  C2  golden dual hashes moved / expected:        {moved} / 48, {moved_expected} predicted from moved positions");
        println!("      topology identical everywhere:              {topology_identical}");
        println!("      worst |hausdorff delta|:                    {worst_hausdorff_delta:.6e}");
        println!(
            "  vacuity: key ties {total_key_ties}, unbreakable {total_unbreakable_ties}, \
             already-passing rows still at 48: {vacuity_rows_still_48} of {} \
             (3 per extractor, both counted)",
            dc_base + mdc_base
        );
        println!("  verdicts: C1 {c1_holds}, C2 {c2_holds}, C3 {c3_holds}");
        for arm in &ARMS {
            println!(
                "  arm {:<40} dc {:>2}/16  mdc {:>2}/16  perm6 {:>2}/32  golden moved {:>2}/48",
                arm.label,
                rows_at_48
                    .get(&(arm.label, "dual_contouring"))
                    .unwrap_or(&0),
                rows_at_48
                    .get(&(arm.label, "manifold_dual_contouring"))
                    .unwrap_or(&0),
                rows_at_perm_6.get(arm.label).unwrap_or(&0),
                hashes_moved.get(arm.label).unwrap_or(&0)
            );
        }

        let mut aggregates: Vec<(&'static str, String)> = vec![
            ("hashes_moved", moved.to_string()),
            ("hashes_moved_expected", moved_expected.to_string()),
            ("c1_holds", c1_holds.to_string()),
            ("c2_holds", c2_holds.to_string()),
            ("c3_holds", c3_holds.to_string()),
            ("c1_rows_at_48", dc_at_48.to_string()),
            ("c1_population", 16.to_string()),
            ("c3_rows_at_48_manifold", mdc_at_48.to_string()),
            ("baseline_rows_at_48_dual_contouring", dc_base.to_string()),
            ("baseline_rows_at_48_manifold", mdc_base.to_string()),
            (
                "vacuity_baseline_rows_still_at_48",
                vacuity_rows_still_48.to_string(),
            ),
            ("total_key_ties", total_key_ties.to_string()),
            (
                "total_unbreakable_key_ties",
                total_unbreakable_ties.to_string(),
            ),
            (
                "max_hausdorff_delta",
                format!("{worst_hausdorff_delta:.12e}"),
            ),
            ("topology_identical", topology_identical.to_string()),
            ("replica_bit_identical", replica_bit_identical.to_string()),
            ("baseline_matches_p61", baseline_matches_p61.to_string()),
            (
                "golden_fixture_matches_shipped",
                golden_fixture_matches_shipped.to_string(),
            ),
            ("golden_counts_changed", golden_counts_changed.to_string()),
            ("golden_dual_rows", 48.to_string()),
            (
                "share_reachable_before_run",
                "false_the_accumulation_is_already_invariant".to_string(),
            ),
        ];
        // Every arm's headline count on every row, so a reader does not have to
        // filter the file to see whether the mechanism arm moved. That is
        // `p-59.csv`'s shape and `p-61.csv` kept it.
        for arm in &ARMS {
            let dc = *rows_at_48
                .get(&(arm.label, "dual_contouring"))
                .unwrap_or(&0);
            let mdc = *rows_at_48
                .get(&(arm.label, "manifold_dual_contouring"))
                .unwrap_or(&0);
            let g = *hashes_moved.get(arm.label).unwrap_or(&0);
            aggregates.push((arm.col_dc, dc.to_string()));
            aggregates.push((arm.col_mdc, mdc.to_string()));
            aggregates.push((arm.col_golden, g.to_string()));
            aggregates.push((
                arm.col_perm,
                rows_at_perm_6.get(arm.label).unwrap_or(&0).to_string(),
            ));
        }

        for mut row in rows.into_iter().chain(golden_rows) {
            row.extend(aggregates.iter().cloned());
            run.record(&row);
        }
    });
}
