//! **P-55 — a 2D monotonicity theorem, ported to 3D and measured on mesh edges.**
//!
//! Ticket: R-050. Pre-registered in the commit before this harness existed.
//!
//! ```bash
//! cargo bench --bench experiment_p55
//! ```
//!
//! Writes `docs/experiments/p-55.csv`.
//!
//! # The predicate, exactly as registered
//!
//! ```text
//! is_monotone(f, a, b, w):
//!     k    = max(2, ceil(‖b − a‖ / w) + 1)
//!     d    = b − a
//!     g(t) = ∇f(a + t·d) · d            (chain rule: this is d/dt f(a + t·d))
//!     non-monotone  ⟺  two sampled g(tᵢ) disagree in sign
//! ```
//!
//! `t₀ = 0` and `t_{k−1} = 1` evaluate at `a` and `b` **exactly** rather than at
//! `a + 0·d` and `a + 1·d`, so an endpoint's directional derivative is the one at
//! the vertex the mesh actually has rather than at a rounded copy of it.
//!
//! ## `w` is the cell size, and that is a choice with a reason
//!
//! `w` is the sampling pitch along the edge: `k` grows so that consecutive
//! samples are no further apart than `w`. Setting it to the grid's cell size
//! makes the rate at which the edge is probed equal the rate at which the field
//! was probed to build the mesh in the first place — the paper's own resolution
//! argument, that a monotonicity certificate is only as fine as the sampling
//! underneath it. Any smaller `w` would test the field at a frequency the
//! extraction never saw; any larger one would test it at less than the
//! extraction saw. Mesh edges are chords inside one cell, so in practice this
//! puts `k` at 2 or 3, and `k_samples_min`/`k_samples`/`k_samples_mean` report
//! the distribution rather than asserting it.
//!
//! # What the source paper says, and what it does not
//!
//! Finken, Li, Wang, Guo & Levine (arXiv:2608.12142) prove **Theorem 1**: a PL
//! function monotonic with respect to a Morse `f` has no spurious critical
//! points. Three things about the transport of that statement to this file, all
//! of them the registration's own admissions and none of them repairable here:
//!
//! - **The theorem is 2D.** Its pigeonhole step is *"since a triangle has only
//!   three edges"*, which has no hexahedral analogue. What runs below is a **3D
//!   port**, and a clause that holds is evidence about this crate's meshes — it
//!   is not the theorem transported to 3D.
//! - **It does not apply to the trilinear interpolant.** Interior critical points
//!   genuinely exist under trilinear interpolation; that is the origin of the
//!   ambiguous-face problem this crate spends `marching_cubes/` on. So even a
//!   perfect result here would certify nothing about the interpolant marching
//!   cubes actually contours.
//! - **The epsilon is isomesh's, not theirs.** The paper gives a bare
//!   sign-disagreement predicate with no epsilon, no relative tolerance and no
//!   flat-region guard. The rule below — discard `|g|` under
//!   `coef · (|f(a)| + |f(b)|)`, with `coef = 1e-12` fixed by the registration
//!   before this harness existed and `1e-14`/`1e-10` recorded beside it — must
//!   not be attributed to Finken et al.
//!
//! ## The gradient substitution, and a place the registration is wrong
//!
//! The paper obtains `∇f` by **autodiff on a neural field**. This harness calls
//! [`Sdf::gradient`], and that substitution changes the noise story in a way the
//! paper never analyses.
//!
//! The registration describes the substitute as *"the crate's central
//! difference"*, and **that is not what runs**: `fields/mod.rs` opens with *"every
//! one of them overrides [`Sdf::gradient`] with an analytic gradient. The
//! central-difference default is never used by a reference field"*, and it is
//! accurate — all eight entries of `for_each_reference_field!` reach an analytic
//! gradient, the three composed ones (`csg_difference`, `gyroid`, `noise_cavity`)
//! by delegating to whichever operand is active. So every `∇f` below is exact
//! rather than `O(h²)`, `gradient_evals` is one field evaluation each rather than
//! six, and the central-difference noise the registration worried about is
//! absent. The `Sdf::gradient` *call* is the registered one; the implementation
//! it lands in is not the registered one.
//!
//! What the analytic gradients do carry instead is **discontinuity**: `box_exact`
//! and `thin_plate` switch face normals across a box edge and across the interior
//! medial axis, and the three CSG composites switch operands across a seam. Those
//! are deterministic selections from a subdifferential, so the sign of `g` next to
//! one of them is a selection rather than a limit.
//!
//! # What this harness found the registration got wrong
//!
//! Not a failure of the mechanism and not a tuning problem: **applied to the
//! edges of the extracted surface, the predicate is saturated by chord geometry
//! and C1's zero is unreachable for any curved field.**
//!
//! A mesh edge joins two vertices that both lie on the extracted zero set, so it
//! is a *chord* of the surface. For a strictly convex surface every chord's
//! interior is strictly inside the solid, so `f` along the edge goes
//! `0 → negative → 0`: it cannot be monotone, and `g` must change sign. On a
//! sphere of radius `r` the two endpoint derivatives are available in closed
//! form — with `a·b = r² cos θ` for central angle `θ`,
//!
//! ```text
//! g(0) = (a·b − r²)/r = r(cos θ − 1) < 0
//! g(1) = (r² − a·b)/r = r(1 − cos θ) > 0
//! ```
//!
//! — for **every** pair of distinct points on the sphere. So the flag is a
//! theorem about chords, not an observation about the mesh, and the count is
//! bounded below by the number of non-degenerate edges.
//!
//! Three measurements say it is that and not something else:
//!
//! - `sphere` at 17³ reaches `non_monotone_per_1k = 1000.000` exactly — all 804
//!   of 804 edges — and stays near 900 at every finer level.
//! - `worst_reversal` on `sphere` falls as `h²`: `7.092e-2`, `1.997e-2`,
//!   `5.420e-3`, `1.408e-3` across the ladder, ratios `3.55`, `3.68`, `3.85`
//!   converging on `4`. That is `r(1 − cos θ) ≈ r θ²/2` with `θ ∝ h`, which is
//!   the sagitta of the chord and nothing else.
//! - `flagged_by_endpoints_at_1e12` settles the registration's own reading of a
//!   falsified C1 — *"k is the problem"*. **22,573 of the sweep's 901,583 flags
//!   (2.5%) needed an interior sample at all**; the other 97.5% were decided by
//!   `g(0)` and `g(1)` before any interior sample existed, so no choice of `k`
//!   and no `w` could have unflagged them. `all_flags_from_endpoints` is `true`
//!   on 17 of the 32 rows outright, `sphere` and `box_exact` at 65³ among them.
//!
//! `box_exact` is the instructive exception in both directions. Its faces are
//! planar and its analytic gradient there is an axis-aligned unit vector, while
//! the chord's component along that axis is an exact `f64` zero — so `g` is
//! *exactly* zero (`zero_g_samples` is 38,520 at 65³) and the flat-region guard
//! the paper omits is not needed at all. What is left is the corner population:
//! an edge straddling a convex box edge picks up one face normal at each end,
//! and its `worst_reversal` is exactly the cell size at every resolution
//! (`0.25`, `0.125`, `0.0625`, `0.03125`), i.e. `worst_reversal_over_w == 1`.
//! Because box edges are `O(n)` against `O(n²)` surface cells, `box_exact` is
//! the one field whose rate genuinely halves per refinement — the behaviour C2
//! predicted for `gyroid` and `fbm_terrain`, on a field C1 predicted zero for.
//!
//! ## The tolerance is inert, for the same structural reason
//!
//! The registered guard scales by `|f(a)| + |f(b)|` — the interpolation residual
//! at the two mesh *vertices*, which is what "both endpoints are on the extracted
//! zero set" means numerically. Measured at 65³ that scale runs from **exactly
//! zero** on `box_exact` and `3.26e-18` on `thin_plate` up to `1.54e-2` on
//! `noise_cavity`, and over all 32 rows `max_abs_tolerance` never exceeds
//! `7.63e-13` — while the smallest *deciding* magnitude anywhere in the sweep,
//! the minimum `worst_reversal`, is `1.41e-3` (`sphere` at 129³). Nine orders of
//! magnitude separate the guard from the thing it is meant to guard.
//!
//! The consequence is measured rather than argued. `nonzero_g_discarded_at_*` is
//! between 0 and 21 per row, `guard_inert` is `true` on 18 of 32 rows, and the
//! `1e-14`/`1e-12`/`1e-10` counts are **identical on 31 of the 32 rows**. The one
//! exception is `gyroid` at 129³, where `1e-10` removes 6 edges of 212,739 —
//! `0.003%`, and in the direction that makes the count smaller. The registration
//! asked for the sensitivity of the answer to be visible; it is visible, and
//! across four orders of magnitude of tolerance it is nil.
//!
//! ## Where the predicate would have to be applied instead
//!
//! The paper's PL function lives on the **ambient** simplicial complex and its
//! monotonicity condition is over that complex's edges. Ported to a hexahedral
//! grid the corresponding objects are the grid's own edges and diagonals, whose
//! endpoints are *not* on the zero set and for which `|f(a)| + |f(b)|` is a real
//! scale. This harness does not measure that — it measures what was registered.
//! Naming it is a finding, not a substitution.
//!
//! # Which edges, and why the count is reported twice
//!
//! [`MeshBuffer`] never welds, so marching cubes emits three fresh vertices per
//! triangle and an edge shared by two triangles appears twice with distinct
//! indices. Welding by proximity would need an epsilon this experiment has no
//! business inventing, so edges are deduplicated by the **exact bit pattern** of
//! their two endpoint positions, ordered canonically with
//! [`f64::total_cmp`] — a total order, so a NaN coordinate would sort into view
//! rather than being dropped by a partial comparison.
//!
//! `edges` is therefore the mesh's set of distinct segments and is what every
//! clause reads; `edge_instances` is `3 · triangles` and `non_monotone_instances`
//! the same count over instances. The predicate is invariant under swapping `a`
//! and `b` — reversing `d` negates every `g(t)`, which leaves *disagreement*
//! alone — so the two readings differ only in the denominator and each is
//! evaluated once per distinct segment.
//!
//! Exact-bit deduplication turns out to be a free check on marching cubes
//! itself, and it passes: on the five **closed** fields — `sphere`, `torus`,
//! `box_exact`, `csg_difference`, `thin_plate` — `edge_instances` is exactly
//! `2 · edges` on all 20 rows, so every mesh edge is shared by exactly two
//! triangles *and* the vertex a shared grid edge produces is bit-identical from
//! both incident cells. The three that deviate deviate for their own documented
//! reasons: `fbm_terrain` is open, so it carries `92`/`190`/`390`/`788`
//! single-triangle boundary edges up the ladder, while `gyroid` (`+10` at every
//! resolution) and `noise_cavity` (`+132` to `+200`) carry a small population of
//! edges shared by more than two triangles — the non-manifold residue those two
//! fields are in the fixture set to produce.
//!
//! # Determinism
//!
//! One thread, no map iteration, no PRNG, `f64` throughout, and every gated
//! quantity an integer, an exact ratio of integers, or a comparison decided by
//! `total_cmp`. The two `ns` columns are recorded because they are interesting
//! and gate nothing.

mod common;

use std::cmp::Ordering;
use std::time::Instant;

use isomesh::fields::ReferenceField;
use isomesh::marching_cubes::MarchingCubes;
use isomesh::{MeshBuffer, Sdf};

/// Samples per axis. The crate's golden ladder, plus the step past it that C2
/// needs in order to be a trend rather than two points.
const GRIDS: [u32; 4] = [17, 33, 65, 129];

/// The three tolerance coefficients, with the registered one in the middle so
/// the CSV reads as a sensitivity strip.
const COEFFS: [f64; 3] = [1e-14, 1e-12, 1e-10];

/// Index of the coefficient the registration fixes.
const REGISTERED: usize = 1;

/// The one extractor this experiment is about.
const EXTRACTOR: &str = "marching_cubes";

#[inline]
fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

#[inline]
fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

/// Lexicographic total order on a point. `total_cmp` rather than `partial_cmp`
/// so a NaN coordinate is ordered rather than silently ungrouped.
fn lex(a: &[f64; 3], b: &[f64; 3]) -> Ordering {
    a[0].total_cmp(&b[0])
        .then_with(|| a[1].total_cmp(&b[1]))
        .then_with(|| a[2].total_cmp(&b[2]))
}

/// One undirected mesh edge, endpoints in canonical order.
#[derive(Clone, Copy)]
struct Segment {
    a: [f64; 3],
    b: [f64; 3],
}

impl Segment {
    fn new(p: [f64; 3], q: [f64; 3]) -> Self {
        if lex(&p, &q) == Ordering::Greater {
            Self { a: q, b: p }
        } else {
            Self { a: p, b: q }
        }
    }

    fn order(&self, other: &Self) -> Ordering {
        lex(&self.a, &other.a).then_with(|| lex(&self.b, &other.b))
    }
}

/// The mesh's distinct segments, and how many triangle sides each carries.
fn segments(mesh: &MeshBuffer<f64>) -> (Vec<Segment>, Vec<u32>) {
    let mut all = Vec::with_capacity(mesh.indices.len());
    for tri in mesh.indices.as_chunks::<3>().0 {
        let p0 = mesh.positions[tri[0] as usize];
        let p1 = mesh.positions[tri[1] as usize];
        let p2 = mesh.positions[tri[2] as usize];
        all.push(Segment::new(p0, p1));
        all.push(Segment::new(p1, p2));
        all.push(Segment::new(p2, p0));
    }
    all.sort_unstable_by(Segment::order);

    let mut unique: Vec<Segment> = Vec::new();
    let mut multiplicity: Vec<u32> = Vec::new();
    for s in all {
        let same = unique
            .last()
            .is_some_and(|last| s.order(last) == Ordering::Equal);
        if same {
            if let Some(m) = multiplicity.last_mut() {
                *m += 1;
            }
        } else {
            unique.push(s);
            multiplicity.push(1);
        }
    }
    (unique, multiplicity)
}

/// What one edge contributes, at all three tolerances from a single pass of
/// gradients — the `g(tᵢ)` do not depend on the coefficient, only the discard
/// rule does.
#[derive(Clone, Copy)]
struct EdgeOutcome {
    k: u32,
    /// One flag per entry of [`COEFFS`].
    non_monotone: [bool; 3],
    /// How many `g` the guard threw away, per coefficient.
    discarded: [u32; 3],
    /// How many of those were **not exactly zero**, per coefficient. This is
    /// the guard's only substantive action: discarding an exact zero changes
    /// nothing, because a zero is neutral in the sign test either way.
    nonzero_discarded: [u32; 3],
    /// Flagged by the two **endpoints alone**, at the registered coefficient —
    /// `g(0)` and `g(1)` already disagree and no interior sample is needed.
    /// When this equals the flag, `k` is not what decided the row.
    flagged_by_endpoints: bool,
    /// Weaker side of the disagreement at the registered coefficient: the
    /// distance from zero the *minority* sign reached. A tiny value means the
    /// flag is noise; a large one means the field genuinely turns around.
    reversal: f64,
    /// `weaker / stronger`, dimensionless, at the registered coefficient.
    reversal_ratio: f64,
    /// `|f(a)| + |f(b)|` — the quantity the registered tolerance scales by.
    endpoint_abs_f: f64,
    /// Smallest strictly positive `|g|` on this edge, or infinity if all are zero.
    min_abs_g: f64,
    zero_g: u32,
    nonfinite_g: u32,
    degenerate: bool,
}

/// The registered predicate, evaluated once for one edge.
fn examine<F>(field: &F, seg: Segment, w: f64, g: &mut Vec<f64>) -> EdgeOutcome
where
    F: Sdf<Scalar = f64>,
{
    let d = sub(seg.b, seg.a);
    let len = dot(d, d).sqrt();
    // A saturating cast, so a non-finite length lands on `k = 2` rather than on
    // an arbitrary index; `nonfinite_g` is what would then report it.
    let k = 2.max((len / w).ceil() as u32 + 1);

    g.clear();
    for i in 0..k {
        let p = if i == 0 {
            seg.a
        } else if i == k - 1 {
            seg.b
        } else {
            let t = f64::from(i) / f64::from(k - 1);
            [
                seg.a[0] + t * d[0],
                seg.a[1] + t * d[1],
                seg.a[2] + t * d[2],
            ]
        };
        g.push(dot(field.gradient(p), d));
    }

    let endpoint_abs_f = field.sample(seg.a).abs() + field.sample(seg.b).abs();

    let mut out = EdgeOutcome {
        k,
        non_monotone: [false; 3],
        discarded: [0; 3],
        nonzero_discarded: [0; 3],
        flagged_by_endpoints: false,
        reversal: 0.0,
        reversal_ratio: 0.0,
        endpoint_abs_f,
        min_abs_g: f64::INFINITY,
        zero_g: 0,
        nonfinite_g: 0,
        // A total comparison, so a NaN endpoint is a degenerate edge rather
        // than an unordered one that quietly counts as ordinary.
        degenerate: !len.is_finite() || len.total_cmp(&0.0) != Ordering::Greater,
    };
    for &v in g.iter() {
        if !v.is_finite() {
            out.nonfinite_g += 1;
            continue;
        }
        let m = v.abs();
        if m > 0.0 {
            if m < out.min_abs_g {
                out.min_abs_g = m;
            }
        } else {
            out.zero_g += 1;
        }
    }

    for (c, &coef) in COEFFS.iter().enumerate() {
        let tol = coef * endpoint_abs_f;
        // Largest kept value of each sign. Zero is neither: an exactly flat
        // sample is not a disagreement with anything.
        let mut pos = 0.0_f64;
        let mut neg = 0.0_f64;
        for &v in g.iter() {
            if v.abs() < tol {
                out.discarded[c] += 1;
                if v.abs() > 0.0 {
                    out.nonzero_discarded[c] += 1;
                }
                continue;
            }
            if v > pos {
                pos = v;
            }
            if -v > neg {
                neg = -v;
            }
        }
        let flagged = pos > 0.0 && neg > 0.0;
        out.non_monotone[c] = flagged;
        if c == REGISTERED {
            let (g0, g1) = (g.first().copied(), g.last().copied());
            out.flagged_by_endpoints = match (g0, g1) {
                (Some(x), Some(y)) => {
                    let keep_x = x.abs() >= tol;
                    let keep_y = y.abs() >= tol;
                    keep_x && keep_y && ((x > 0.0 && y < 0.0) || (x < 0.0 && y > 0.0))
                }
                _ => false,
            };
            if flagged {
                out.reversal = pos.min(neg);
                out.reversal_ratio = pos.min(neg) / pos.max(neg);
            }
        }
    }
    out
}

/// One CSV row's worth of measurement, before the cross-row booleans exist.
struct Measured {
    field: &'static str,
    samples: u32,
    cell: f64,
    triangles: u64,
    edge_instances: u64,
    edges: u64,
    non_monotone: [u64; 3],
    non_monotone_instances: u64,
    per_1k: [f64; 3],
    k_min: u32,
    k_max: u32,
    k_sum: u64,
    gradient_evals: u64,
    worst_reversal: f64,
    worst_reversal_ratio: f64,
    discarded: [u64; 3],
    nonzero_discarded: [u64; 3],
    edges_with_discard: [u64; 3],
    flagged_by_endpoints: u64,
    tol_sum: f64,
    tol_max: f64,
    endpoint_abs_f_sum: f64,
    endpoint_abs_f_max: f64,
    min_abs_g: f64,
    zero_g: u64,
    nonfinite_g: u64,
    degenerate: u64,
    extract_ns: f64,
    predicate_ns: f64,
}

impl Measured {
    fn per_1k_registered(&self) -> f64 {
        self.per_1k[REGISTERED]
    }
}

/// Extract at `samples³` and run the predicate over every distinct mesh edge.
fn measure<F>(name: &'static str, field: &F, samples: u32) -> Measured
where
    F: ReferenceField + Sdf<Scalar = f64>,
{
    let (shape, origin, cell) = common::grid(field, samples);

    let mut mesh = MeshBuffer::<f64>::new();
    let mut mc = MarchingCubes::<f64>::new();
    let started = Instant::now();
    mc.extract(field, &shape, origin, cell, &mut mesh)
        .expect("marching cubes on a reference field");
    let extract_ns = started.elapsed().as_secs_f64() * 1e9;

    let (unique, multiplicity) = segments(&mesh);

    let mut row = Measured {
        field: name,
        samples,
        cell,
        triangles: mesh.triangle_count() as u64,
        edge_instances: multiplicity.iter().copied().map(u64::from).sum(),
        edges: unique.len() as u64,
        non_monotone: [0; 3],
        non_monotone_instances: 0,
        per_1k: [0.0; 3],
        k_min: u32::MAX,
        k_max: 0,
        k_sum: 0,
        gradient_evals: 0,
        worst_reversal: 0.0,
        worst_reversal_ratio: 0.0,
        discarded: [0; 3],
        nonzero_discarded: [0; 3],
        edges_with_discard: [0; 3],
        flagged_by_endpoints: 0,
        tol_sum: 0.0,
        tol_max: 0.0,
        endpoint_abs_f_sum: 0.0,
        endpoint_abs_f_max: 0.0,
        min_abs_g: f64::INFINITY,
        zero_g: 0,
        nonfinite_g: 0,
        degenerate: 0,
        extract_ns,
        predicate_ns: 0.0,
    };

    let mut g: Vec<f64> = Vec::with_capacity(8);
    let started = Instant::now();
    for (seg, &mult) in unique.iter().zip(multiplicity.iter()) {
        let out = examine(field, *seg, cell, &mut g);

        row.k_min = row.k_min.min(out.k);
        row.k_max = row.k_max.max(out.k);
        row.k_sum += u64::from(out.k);
        row.gradient_evals += u64::from(out.k);
        row.zero_g += u64::from(out.zero_g);
        row.nonfinite_g += u64::from(out.nonfinite_g);
        row.degenerate += u64::from(out.degenerate);
        if out.min_abs_g < row.min_abs_g {
            row.min_abs_g = out.min_abs_g;
        }
        row.endpoint_abs_f_sum += out.endpoint_abs_f;
        if out.endpoint_abs_f > row.endpoint_abs_f_max {
            row.endpoint_abs_f_max = out.endpoint_abs_f;
        }
        let tol = COEFFS[REGISTERED] * out.endpoint_abs_f;
        row.tol_sum += tol;
        if tol > row.tol_max {
            row.tol_max = tol;
        }
        for c in 0..COEFFS.len() {
            row.discarded[c] += u64::from(out.discarded[c]);
            row.nonzero_discarded[c] += u64::from(out.nonzero_discarded[c]);
            if out.discarded[c] > 0 {
                row.edges_with_discard[c] += 1;
            }
            if out.non_monotone[c] {
                row.non_monotone[c] += 1;
            }
        }
        if out.flagged_by_endpoints {
            row.flagged_by_endpoints += 1;
        }
        if out.non_monotone[REGISTERED] {
            row.non_monotone_instances += u64::from(mult);
            if out.reversal > row.worst_reversal {
                row.worst_reversal = out.reversal;
                row.worst_reversal_ratio = out.reversal_ratio;
            }
        }
    }
    row.predicate_ns = started.elapsed().as_secs_f64() * 1e9;

    if row.edges == 0 {
        row.k_min = 0;
    }
    if !row.min_abs_g.is_finite() {
        row.min_abs_g = 0.0;
    }
    for c in 0..COEFFS.len() {
        // An empty mesh has no rate. Recorded as zero rather than as a NaN,
        // with `edges = 0` on the same row saying why it is vacuous.
        row.per_1k[c] = if row.edges == 0 {
            0.0
        } else {
            1000.0 * row.non_monotone[c] as f64 / row.edges as f64
        };
    }
    row
}

/// Strictly decreasing across the resolution ladder — C2's boolean, computed
/// from the four measured values rather than asserted.
fn falls(v: &[u64]) -> bool {
    v.windows(2).all(|w| w[1] < w[0])
}

/// The same, for a rate. `total_cmp` so a NaN cannot pass for "less".
fn falls_rate(v: &[f64]) -> bool {
    v.windows(2)
        .all(|w| w[1].total_cmp(&w[0]) == Ordering::Less)
}

fn main() {
    let prereg = isomesh::experiment!("P-55");

    let mut rows: Vec<Measured> = Vec::new();
    isomesh::for_each_reference_field!(f64, |name, field| {
        for samples in GRIDS {
            let row = measure(name, &field, samples);
            println!(
                "{:>14} {:>4}³  tris {:>7}  edges {:>8} ({:>8} inst)  \
                 non-monotone {:>8}  per-1k {:>8.3}  [1e-14 {:>8}  1e-10 {:>8}]  \
                 k {}..{}  worst reversal {:>10.3e}  tol≤{:>9.2e}  \
                 discarded {}  zero-g {}",
                row.field,
                row.samples,
                row.triangles,
                row.edges,
                row.edge_instances,
                row.non_monotone[REGISTERED],
                row.per_1k[REGISTERED],
                row.non_monotone[0],
                row.non_monotone[2],
                row.k_min,
                row.k_max,
                row.worst_reversal,
                row.tol_max,
                row.discarded[REGISTERED],
                row.zero_g,
            );
            rows.push(row);
        }
    });

    // ── the cross-row booleans, all from measured values ────────────────────
    let mut names: Vec<&'static str> = Vec::new();
    for row in &rows {
        if !names.contains(&row.field) {
            names.push(row.field);
        }
    }

    let mut falls_count: Vec<(&'static str, bool)> = Vec::new();
    let mut falls_per_1k: Vec<(&'static str, bool)> = Vec::new();
    for &name in &names {
        let counts: Vec<u64> = GRIDS
            .iter()
            .filter_map(|&n| {
                rows.iter()
                    .find(|r| r.field == name && r.samples == n)
                    .map(|r| r.non_monotone[REGISTERED])
            })
            .collect();
        let rates: Vec<f64> = GRIDS
            .iter()
            .filter_map(|&n| {
                rows.iter()
                    .find(|r| r.field == name && r.samples == n)
                    .map(Measured::per_1k_registered)
            })
            .collect();
        falls_count.push((name, falls(&counts)));
        falls_per_1k.push((name, falls_rate(&rates)));
    }

    // C3's argmax, per resolution, decided by `total_cmp`.
    let mut highest: Vec<(u32, &'static str)> = Vec::new();
    for &n in &GRIDS {
        let mut best: Option<(&'static str, f64)> = None;
        for row in rows.iter().filter(|r| r.samples == n) {
            let v = row.per_1k_registered();
            let better = match best {
                None => true,
                Some((_, b)) => v.total_cmp(&b) == Ordering::Greater,
            };
            if better {
                best = Some((row.field, v));
            }
        }
        highest.push((n, best.map_or("none", |(f, _)| f)));
    }

    common::experiment::run(prereg, |run| {
        for row in &rows {
            let falls_here = falls_count
                .iter()
                .find(|(f, _)| *f == row.field)
                .is_some_and(|&(_, b)| b);
            let falls_rate_here = falls_per_1k
                .iter()
                .find(|(f, _)| *f == row.field)
                .is_some_and(|&(_, b)| b);
            let top = highest
                .iter()
                .find(|(n, _)| *n == row.samples)
                .map_or("none", |&(_, f)| f);
            let edges = row.edges.max(1) as f64;

            run.record(&[
                // ── the twelve registered metrics ──────────────────────────
                ("field", row.field.to_string()),
                ("extractor", EXTRACTOR.to_string()),
                ("samples_per_axis", row.samples.to_string()),
                ("edges", row.edges.to_string()),
                (
                    "non_monotone_edges",
                    row.non_monotone[REGISTERED].to_string(),
                ),
                (
                    "non_monotone_per_1k",
                    format!("{:.6}", row.per_1k[REGISTERED]),
                ),
                // The largest k any edge on this row needed; the spread is in
                // `k_samples_min` and `k_samples_mean` beside it.
                ("k_samples", row.k_max.to_string()),
                ("tolerance", format!("{:e}", COEFFS[REGISTERED])),
                ("non_monotone_at_1e14", row.non_monotone[0].to_string()),
                ("non_monotone_at_1e10", row.non_monotone[2].to_string()),
                ("worst_reversal", format!("{:.6e}", row.worst_reversal)),
                ("falls_with_resolution", falls_here.to_string()),
                // ── the sensitivity strip, side by side ────────────────────
                ("tolerance_rule", "coef*(|f(a)|+|f(b)|)".to_string()),
                ("per_1k_at_1e14", format!("{:.6}", row.per_1k[0])),
                ("per_1k_at_1e12", format!("{:.6}", row.per_1k[REGISTERED])),
                ("per_1k_at_1e10", format!("{:.6}", row.per_1k[2])),
                (
                    "counts_equal_across_tolerances",
                    (row.non_monotone[0] == row.non_monotone[REGISTERED]
                        && row.non_monotone[REGISTERED] == row.non_monotone[2])
                        .to_string(),
                ),
                ("g_discarded_at_1e14", row.discarded[0].to_string()),
                ("g_discarded_at_1e12", row.discarded[REGISTERED].to_string()),
                ("g_discarded_at_1e10", row.discarded[2].to_string()),
                (
                    "edges_with_discarded_g_at_1e10",
                    row.edges_with_discard[2].to_string(),
                ),
                // Discarding an exact zero changes nothing — zero is neutral in
                // the sign test with or without the guard. These three are the
                // only substantive thing the guard ever did.
                (
                    "nonzero_g_discarded_at_1e14",
                    row.nonzero_discarded[0].to_string(),
                ),
                (
                    "nonzero_g_discarded_at_1e12",
                    row.nonzero_discarded[REGISTERED].to_string(),
                ),
                (
                    "nonzero_g_discarded_at_1e10",
                    row.nonzero_discarded[2].to_string(),
                ),
                (
                    "guard_inert",
                    (row.nonzero_discarded[0] == 0
                        && row.nonzero_discarded[REGISTERED] == 0
                        && row.nonzero_discarded[2] == 0)
                        .to_string(),
                ),
                // ── why the guard is or is not inert ───────────────────────
                ("max_abs_tolerance", format!("{:.6e}", row.tol_max)),
                ("mean_abs_tolerance", format!("{:.6e}", row.tol_sum / edges)),
                (
                    "mean_endpoint_abs_f",
                    format!("{:.6e}", row.endpoint_abs_f_sum / edges),
                ),
                (
                    "max_endpoint_abs_f",
                    format!("{:.6e}", row.endpoint_abs_f_max),
                ),
                ("min_abs_g", format!("{:.6e}", row.min_abs_g)),
                ("zero_g_samples", row.zero_g.to_string()),
                ("nonfinite_g_samples", row.nonfinite_g.to_string()),
                // ── the mesh, and the two edge readings ───────────────────
                ("triangles", row.triangles.to_string()),
                ("edge_instances", row.edge_instances.to_string()),
                (
                    "non_monotone_instances",
                    row.non_monotone_instances.to_string(),
                ),
                ("degenerate_edges", row.degenerate.to_string()),
                // ── the predicate's own parameters ────────────────────────
                ("cell_size_w", format!("{:.9}", row.cell)),
                ("k_samples_min", row.k_min.to_string()),
                ("k_samples_mean", format!("{:.4}", row.k_sum as f64 / edges)),
                ("gradient_evals", row.gradient_evals.to_string()),
                (
                    "worst_reversal_ratio",
                    format!("{:.6}", row.worst_reversal_ratio),
                ),
                (
                    "worst_reversal_over_w",
                    format!("{:.6}", row.worst_reversal / row.cell),
                ),
                // The registration's own falsification reading for C1 is "k is
                // the problem". These two settle it: an edge flagged by its two
                // endpoints alone was decided before any interior sample
                // existed, so no choice of k could unflag it.
                (
                    "flagged_by_endpoints_at_1e12",
                    row.flagged_by_endpoints.to_string(),
                ),
                (
                    "flagged_only_by_interior_at_1e12",
                    row.non_monotone[REGISTERED]
                        .saturating_sub(row.flagged_by_endpoints)
                        .to_string(),
                ),
                (
                    "all_flags_from_endpoints",
                    (row.flagged_by_endpoints == row.non_monotone[REGISTERED]).to_string(),
                ),
                // ── C2's other reading, and C3 ────────────────────────────
                ("per_1k_falls_with_resolution", falls_rate_here.to_string()),
                ("highest_per_1k_field", top.to_string()),
                ("is_highest_per_1k", (top == row.field).to_string()),
                // ── speed, recorded beside the verdict, gating nothing ────
                ("extract_ns", format!("{:.0}", row.extract_ns)),
                ("predicate_ns", format!("{:.0}", row.predicate_ns)),
                (
                    "predicate_ns_per_edge",
                    format!("{:.2}", row.predicate_ns / edges),
                ),
            ]);
        }
    });
}
