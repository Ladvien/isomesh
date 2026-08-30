//! **P-161 — which approximation class `A^s` each field is in, which decides whether LOD can help it at all.**
//!
//! Ticket: R-161. Pre-registered before this harness existed.
//!
//! ```bash
//! cargo bench --bench experiment_p161
//! ```
//!
//! Writes `docs/experiments/p-161.csv`.
//!
//! # What was missing
//!
//! Two things, and the second is the harder gap.
//!
//! **The theory is in the corpus and has never been pointed at a field.** The
//! registration names four adaptive-approximation papers held at 0.62–0.64, and
//! one of them — Bonito, Canuto, Nochetto & Veeser, Acta Numerica 2024,
//! `10.1017/s0962492924000011` — is already quoted by a *sibling* registration
//! (`P-146`, `FINDINGS.md:25286`) for exactly the sentence this row measures:
//! for `W^{2,p}` regularity uniform refinement gives `O(N^(−2/3))` in 3D and the
//! optimal graded mesh gives `O(N^(−2/3))`, *the same exponent*, and that order
//! **"cannot be improved upon assuming either higher regularity … or a graded
//! mesh"**. So the theory's own prediction is that adaptivity buys a
//! **constant** on a smooth field and a **rate** only where regularity fails.
//! Nobody has checked which of the eight reference fields is which. `Stevenson`
//! is cited in this repository *only* for the newest-vertex-bisection closure
//! bound — `docs/research/2026-08-13-adjacent-field-acquisition.md:91-94`, which
//! notes the bound is global-over-history rather than per-edit, and
//! `docs/research/2026-08-12-axes-and-vocabulary.md:68,75-77` — never for the
//! approximation-class theory in the same PDFs, and the reason is visible:
//! Stevenson's paper is one of the paywalled rows
//! (`docs/research/2026-08-13-acquisition-gaps.md:19`,
//! `docs/research/2026-08-15-paywalled-backlog.md:180`), 11 MB of HTML reduced to
//! 18 KB of MathML noise.
//!
//! **`M-12` is this row's exponent in the wrong currency.** `M-12`
//! (`FINDINGS.md:1128`) measured Marching Cubes' error falling like `h²` — mean
//! `2.7168e-3` at 32³ against `6.5015e-4` at 64³, ratio **4.179** against an
//! ideal 4.13. `P-159` makes the point that field evaluation is the currency
//! actually spent. Re-denominated: `N ~ h^(−3)`, so `h²` **is** `N^(−2/3)`, and
//! `M-12`'s law is the statement `A^(2/3)` for the uniform family. That number
//! has never been written down as an exponent in `N`, and its adaptive
//! counterpart has never been measured at all.
//!
//! **There is no octree-adaptive sampler in the crate, and that is deliberate.**
//! `isomesh::lod` is *downsampling* — building each level from the level below —
//! and its own module doc separates that from the **re-sampling** family and says
//! re-sampling "needs no API: the caller already knows how to sample a field, so
//! adding a `Downsample::Reevaluate` variant would be a second path to something
//! that is not downsampling at all" (`crates/isomesh/src/lod.rs:6-16`). The
//! adaptive arm here is therefore necessarily bench-local, and
//! `crates/isomesh/src/**` is read-only for this phase.
//!
//! # What is measured, and why it is not a mesh Hausdorff
//!
//! `A^s` is a property of **the field under a refinement family**, not of one
//! extractor. The error functional is the defect of the piecewise-trilinear
//! reconstruction of the field itself, scored on a fixed probe set in a band
//! about the zero set:
//!
//! - **Error** — the discrete `L²` (root-mean-square) of
//!   `|f(p) − reconstruct(p)|` over [`PROBES`] probes drawn once per field by
//!   rejection sampling to `|f(p)| ≤ band`, with `band = span/32`. The probe set
//!   is **identical across every resolution and both families**, so the slope is
//!   not measured against a moving target, and the readings are *paired*, so
//!   sampling noise largely cancels in the slope. The `L∞` series is carried
//!   beside it as a control on the norm.
//! - **The band is defined by field VALUE, not distance.** That is the only
//!   definition available on four of the eight fields — `gyroid` is `Lipschitz`,
//!   `csg_difference` is `Underestimate { q: 0.5 }`, `fbm_terrain` and
//!   `noise_cavity` are `Unbounded`, and `fields/mod.rs:104-109` says plainly
//!   that for `Unbounded` "no accuracy figure against `|sample|` means
//!   anything", while `:81-85` says only `Exact` "admits a Hausdorff
//!   measurement". An interpolation-defect measurement is about *values*, so it
//!   is well posed on all eight. `validate::accuracy` is still run, gated on
//!   `bound().is_exact()`, as an independent second instrument — see the arms
//!   table — and the four skips are recorded rather than papered over.
//! - **A fixed band, not a shrinking one, and this decides the answer.** With
//!   the band held fixed as `h → 0` both families are volumetric inside it, so a
//!   smooth field gives `s_uniform = s_adaptive = 2/3` and adaptivity buys only
//!   a constant — which is precisely the Bonito–Canuto–Nochetto–Veeser sentence.
//!   A band shrinking with `h` would hand the adaptive family the surface's
//!   codimension for free and make every field "gain", measuring the functional
//!   rather than the field.
//!
//! # Degrees of freedom: field evaluations, and nothing else
//!
//! `N` is **the number of distinct points at which the field was evaluated to
//! build the reconstruction**. For the uniform family that is exactly `n³`. For
//! the adaptive family it is the size of the sample cache — every corner *and*
//! every cell-centre the refinement indicator asked for. That is the
//! conservative choice: the adaptive family is charged for its own error
//! estimator, so it gets no free oracle. It costs the adaptive arm a constant
//! (measured beside it as `dof_adaptive_corners`, the reconstruction corners
//! alone) and therefore **cannot move the slope**, because a constant factor in
//! `N` is an intercept shift in `ln N`. Triangles were the other candidate and
//! were rejected: `thin_plate` is 0.4 cells thick, so three of the eight fields
//! have no triangle count worth comparing at coarse `h`, and `P-159`'s point is
//! that evaluations are the currency actually spent.
//!
//! # The adaptive strategy, stated so it can be criticised
//!
//! An octree over the field's own domain, root grid [`BASE_CELLS`]³ cells, at
//! most [`MAX_EXTRA_DEPTH`] levels below that. One leaf is refined at a time,
//! into eight children, and the reconstruction is per-leaf trilinear on that
//! leaf's own eight corners — discontinuous across a level boundary, exactly as
//! an LOD octree is, and the error functional does not care.
//!
//! - **Indicator.** `defect × h`, where `defect = |f(centre) − trilinear(corners,
//!   ½,½,½)|` is the classical second-difference interpolation defect and `h` is
//!   the leaf's world edge length. The weight is derived, not tuned: band probes
//!   are spread through a shell of fixed thickness, a leaf of size `h` meeting
//!   the surface carries a probe share proportional to its surface patch `~h²`,
//!   so its contribution to `Σ err²` is `~defect²·h²` and the greedy that
//!   equidistributes the measured `L²` refines the largest `defect·h`.
//! - **Marking.** Greedy maximum, one leaf per step, through a binary heap keyed
//!   on the indicator's bit pattern (non-negative and finite, so `to_bits` is
//!   monotone) with the leaf's integer origin as a deterministic tie-break.
//! - **Activity.** A leaf is refinable only if it can carry band probes: its
//!   corners straddle zero, **or** `min|f| ≤ band + 2·spread` over its eight
//!   corners and its centre, where `spread` is the max-minus-min of those nine
//!   values. The `2` is derived. Any point of a leaf is within `h√3/2` of the
//!   centre, so a probe with `|f| ≤ band` puts `|f(centre)| ≤ band +
//!   |∇f|·h·0.866`; and for a field locally linear on the leaf `spread ≥
//!   |∇f|·h`; so `min|f| ≤ band + 0.866·spread`, comfortably inside the test.
//!   `probes_in_inactive_leaves` measures that derivation rather than trusting
//!   it, and is asserted zero.
//! - **Budgets.** The adaptive run is a single refinement sweep snapshotted the
//!   moment its evaluation count crosses each uniform budget, so the two
//!   families are fitted over **the same `N` values** — asserted, within 5%.
//!
//! # The one way this instrument can lie, and the two columns that catch it
//!
//! The adaptive family earns a **constant** before it earns a rate. Its root
//! grid covers the whole domain and its refinement concentrates into the band,
//! so at fixed `N` it eventually holds a factor of roughly
//! `band_volume_fraction^(−s)` over uniform — and a constant *being earned*
//! across the fitted window is indistinguishable, in a six-point log-log slope,
//! from a higher rate. The window is 2.03 decades and the transition costs some
//! of it, so `s_adaptive` is an **effective** exponent over the achievable range
//! and not an asymptote. That is exactly C1's registered failure mode — `s` not
//! estimable from the achievable resolution range — arriving as a widened
//! interval rather than as an absence, and it is why C2 is scored against the
//! *combined* half-width instead of against zero: a transient inflates
//! `s_adaptive` and widens its interval at the same time, so the conservative
//! test declines rather than claims. `fit_r2_adaptive` and `s_adaptive_tail`
//! (the slope over the last three budgets alone) are the two columns that make
//! the curvature visible, and neither is consulted by a clause.
//!
//! # Arms
//!
//! | arm | what it varies | is_control |
//! |---|---|---|
//! | `uniform` | one global cell size; `n` ∈ [`UNIFORM_SAMPLES`] | no — the baseline family, `s_uniform` |
//! | `adaptive` | octree depth under `defect × h`, at six matched evaluation budgets | no — the treatment family, `s_adaptive` |
//! | `random` | the same octree, marking a **uniformly random** active leaf instead of the largest indicator | **yes** — isolates the indicator from the octree. If `s_random` matches `s_adaptive`, the gain is the tree and not the criterion |
//! | `sup_norm` | both families rescored in `L∞` instead of `L²` | **yes** — the split must not be an artefact of the norm |
//! | `mesh_uniform` | `MarchingCubes` + `validate::accuracy` symmetric Hausdorff against triangle count, four resolutions, `Exact` fields only | **yes** — a second, wholly independent instrument for `s_uniform`; skipped and recorded on the four non-`Exact` fields |
//!
//! Six resolutions rather than the registered floor of four: the three golden
//! resolutions 17/25/33 (`golden.rs:72`) plus 49, 65 and 81, which span
//! `4913 → 531_441` degrees of freedom — **2.03 decades** — and leave four
//! degrees of freedom in each fit, so the 95% interval carries `t = 2.776`
//! rather than the two-point-slope `12.706` a four-resolution sweep would.
//!
//! # SHARE, recomputed before the numbers
//!
//! The registration says **"SHARE: none — this predicts where a stage helps, it
//! does not change the stage."** Discharged, and it is true in the strong sense
//! here: this harness reads no shipped code path other than `MarchingCubes` and
//! `validate::accuracy` in a control arm, changes nothing, and produces no ratio
//! a later row could denominate a speedup in. What it produces is a per-field
//! exponent, which is the *precondition* for the octree-LOD line rather than a
//! share of it: on a field where `s_adaptive = s_uniform` no octree scheme can
//! beat uniform sampling asymptotically, however well implemented, and the
//! number that says so is a column here.
//!
//! # Vacuity controls
//!
//! Every one runs before the first `record` and every panic starts `VOID: `.
//!
//! - **Eight fields.** `rows == 8`, or the "all eight fields" in C1 is a claim
//!   about a subset. Proven by `field`.
//! - **Six resolutions per family per field**, never fewer than the registered
//!   four, and a `ln N` window of at least one decade. Proven by `resolutions`
//!   and `dof_uniform`; without a window there is no slope to fit and C1's
//!   registered failure mode — `s` not estimable from the achievable range —
//!   could not be distinguished from a short sweep.
//! - **A confidence interval on every fit that exists.** Proven by
//!   `s_uniform_ci` and `s_adaptive_ci`, and the split in C2 is scored against
//!   the *combined* half-width in `s_difference_ci` rather than against zero, so
//!   it cannot be fitted noise. A missing fit is **not** asserted away: C1's own
//!   falsifier is `s` not being estimable, so an unfitted field is a registered
//!   outcome and is recorded as `class_membership = unfitted` with
//!   `c1_field = false`.
//! - **The probe set is full and fixed.** `probes == PROBES` on every field, so
//!   no field is scored on a thin or empty population. `band_volume_fraction` is
//!   the acceptance rate that filled it, and is the mechanism column for C2: a
//!   field whose band is nearly the whole domain has nothing to adapt to.
//! - **The adaptive arm actually varies its density.** `adaptive_max_depth ≥ 3`,
//!   i.e. somewhere at least 8× finer than its own root grid — `P-160`'s control
//!   in this row's currency — and `adaptive_truncated == false`, so the series
//!   is not measured against the depth cap.
//! - **Every probe lands in a leaf the criterion was willing to refine.**
//!   `probes_in_inactive_leaves == 0`, or the adaptive error is partly measured
//!   over cells the strategy declined to touch and `s_adaptive` is a property of
//!   the activity test rather than of the field.
//! - **Matched budgets.** Every adaptive snapshot is at or just past its uniform
//!   budget and within 5% of it, or the two slopes are fitted over different
//!   windows and their difference is an artefact.

#![allow(
    clippy::cast_possible_truncation,
    clippy::cast_precision_loss,
    clippy::cast_sign_loss,
    clippy::too_many_lines
)]

mod common;

use std::cmp::Reverse;
use std::collections::{BTreeMap, BinaryHeap};
use std::time::Instant;

use isomesh::fields::{FieldBound, ReferenceField};
use isomesh::marching_cubes::MarchingCubes;
use isomesh::validate::{AccuracyConfig, accuracy};
use isomesh::{MeshBuffer, Sdf};

// ─── the sweep ───────────────────────────────────────────────────────────────

/// Samples per axis for the uniform family. Six, so each fit has four degrees of
/// freedom; the first three are the golden fixture's own resolutions.
const UNIFORM_SAMPLES: [u32; 6] = [17, 25, 33, 49, 65, 81];

/// Resolutions for the `mesh_uniform` control. Four, the registered floor: at
/// 65³ and above `validate::accuracy`'s reverse pass costs more than the whole
/// primary measurement, and the clause it serves is a cross-check not a verdict.
const MESH_SAMPLES: [u32; 4] = [17, 25, 33, 49];

/// Probes per field, drawn once and reused by every arm and every resolution.
const PROBES: usize = 4096;

/// Rejection-sampling attempt cap while filling the probe set.
const PROBE_ATTEMPTS: usize = 4_000_000;

/// Band half-width as a fraction of the domain span.
const BAND_FRACTION: f64 = 1.0 / 32.0;

/// Cells per axis of the adaptive octree's root grid.
const BASE_CELLS: u32 = 4;

/// Levels the octree may add below the root grid.
const MAX_EXTRA_DEPTH: u32 = 9;

/// Integer lattice resolution, one unit being half the finest cell edge.
const UNITS: u32 = BASE_CELLS << (MAX_EXTRA_DEPTH + 1);

/// [`UNITS`] as a float, written out and checked rather than cast.
const UNITS_F: f64 = 4096.0;

const _UNITS_AGREE: () = assert!(UNITS == 4096);

/// Slack multiplier in the activity test. Derived in the header, not tuned.
const ACTIVITY_SLACK: f64 = 2.0;

/// Depth a leaf must reach before the adaptive arm counts as adaptive.
const MIN_ADAPTIVE_DEPTH: u32 = 3;

/// Ceiling on how far an adaptive snapshot may overshoot its uniform budget.
const BUDGET_OVERSHOOT: f64 = 1.05;

/// Two-sided 97.5% Student-t quantiles, indexed by degrees of freedom.
///
/// Index 0 is unreachable — a two-point fit has no residual — and is `NaN` so
/// that reaching it poisons the interval instead of inventing one.
const T_975: [f64; 10] = [
    f64::NAN,
    12.706,
    4.303,
    3.182,
    2.776,
    2.571,
    2.447,
    2.365,
    2.306,
    2.262,
];

/// Base seed for the probe sets, salted per field by its name.
const SEED: u64 = 0x1_6161_6161;

/// Salt separating the `random` control's marking stream from the probe stream.
const RANDOM_SALT: u64 = 0xA5A5_A5A5_A5A5_A5A5;

// ─── determinism ─────────────────────────────────────────────────────────────

/// `SplitMix64`. Ten lines, seeded, and written here because `rand` is not a
/// dependency and a bench that cannot be replayed is not a measurement.
#[derive(Clone, Copy, Debug)]
struct Rng(u64);

impl Rng {
    const fn new(seed: u64) -> Self {
        Self(seed)
    }

    fn next_u64(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }

    /// Uniform in `[0, 1)`, 53 mantissa bits.
    fn unit(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 * (1.0 / 9_007_199_254_740_992.0)
    }
}

/// FNV-1a over the field name, xored with [`SEED`].
///
/// A per-field seed, so two fields never share a probe pattern, derived from the
/// name rather than from an enumeration order a later field would shift.
fn seed_for(name: &str) -> u64 {
    let mut h = 0xcbf2_9ce4_8422_2325_u64;
    for b in name.as_bytes() {
        h ^= u64::from(*b);
        h = h.wrapping_mul(0x0000_0100_0000_01b3);
    }
    h ^ SEED
}

// ─── the reconstruction ──────────────────────────────────────────────────────

/// The eight-corner trilinear interpolant. Corners are indexed `x + 2y + 4z`,
/// the crate's own corner order (`marching_cubes/table.rs:88-91`'s re-exported
/// `EDGE_CORNERS` convention; `cube::corner_offset` itself is private).
fn trilinear(c: &[f64; 8], u: f64, v: f64, w: f64) -> f64 {
    let mut acc = 0.0;
    for (i, &ci) in c.iter().enumerate() {
        let wx = if i & 1 == 0 { 1.0 - u } else { u };
        let wy = if (i >> 1) & 1 == 0 { 1.0 - v } else { v };
        let wz = if (i >> 2) & 1 == 0 { 1.0 - w } else { w };
        acc += ci * wx * wy * wz;
    }
    acc
}

/// The cubic box a field is measured over, plus the band that defines its probes.
#[derive(Clone, Copy, Debug)]
struct Domain {
    /// Minimum corner.
    lo: [f64; 3],
    /// Edge length. Every reference field's domain is a cube; asserted.
    span: f64,
    /// Probes satisfy `|f| <= band`.
    band: f64,
}

/// The fixed population every error in this bench is scored over.
#[derive(Clone, Debug)]
struct Probes {
    /// Probe positions, in world coordinates.
    points: Vec<[f64; 3]>,
    /// `field.sample(point)`, the ground truth the reconstruction is compared to.
    truth: Vec<f64>,
    /// Rejection-sampling attempts that produced them.
    attempts: usize,
}

impl Probes {
    /// Acceptance rate, which is the band's share of the domain volume.
    fn band_volume_fraction(&self) -> f64 {
        self.points.len() as f64 / self.attempts as f64
    }
}

/// Draw the probe set: uniform in the box, kept when `|f| <= band`.
fn probe_set<F: Sdf<Scalar = f64>>(field: &F, dom: Domain, seed: u64) -> Probes {
    let mut rng = Rng::new(seed);
    let mut points = Vec::with_capacity(PROBES);
    let mut truth = Vec::with_capacity(PROBES);
    let mut attempts = 0usize;
    while points.len() < PROBES && attempts < PROBE_ATTEMPTS {
        attempts += 1;
        let x = rng.unit();
        let y = rng.unit();
        let z = rng.unit();
        let p = [
            dom.lo[0] + dom.span * x,
            dom.lo[1] + dom.span * y,
            dom.lo[2] + dom.span * z,
        ];
        let v = field.sample(p);
        if v.is_finite() && v.abs() <= dom.band {
            points.push(p);
            truth.push(v);
        }
    }
    Probes {
        points,
        truth,
        attempts,
    }
}

// ─── the uniform family ──────────────────────────────────────────────────────

/// Root-mean-square and sup of the trilinear defect on a uniform `n³` grid.
///
/// The grid is sampled in full — `n³` evaluations, which is this arm's own
/// degrees of freedom — and every probe is read out of the cell containing it.
fn uniform_error<F: Sdf<Scalar = f64>>(
    field: &F,
    dom: Domain,
    n: u32,
    probes: &Probes,
) -> (f64, f64) {
    let nn = n as usize;
    let h = dom.span / f64::from(n - 1);
    let mut values = vec![0.0f64; nn * nn * nn];
    for z in 0..nn {
        for y in 0..nn {
            for x in 0..nn {
                let p = [
                    dom.lo[0] + h * x as f64,
                    dom.lo[1] + h * y as f64,
                    dom.lo[2] + h * z as f64,
                ];
                values[x + nn * (y + nn * z)] = field.sample(p);
            }
        }
    }

    let mut sum = 0.0;
    let mut sup = 0.0;
    for (p, &tv) in probes.points.iter().zip(&probes.truth) {
        let cell: [(usize, f64); 3] = std::array::from_fn(|a| {
            let g = (p[a] - dom.lo[a]) / h;
            let c = g.floor().clamp(0.0, f64::from(n - 2));
            (c as usize, (g - c).clamp(0.0, 1.0))
        });
        let base = cell[0].0 + nn * (cell[1].0 + nn * cell[2].0);
        let corners: [f64; 8] = std::array::from_fn(|i| {
            values[base + (i & 1) + nn * (((i >> 1) & 1) + nn * ((i >> 2) & 1))]
        });
        let e = (trilinear(&corners, cell[0].1, cell[1].1, cell[2].1) - tv).abs();
        sum += e * e;
        if e > sup {
            sup = e;
        }
    }
    ((sum / probes.points.len() as f64).sqrt(), sup)
}

// ─── the adaptive family ─────────────────────────────────────────────────────

/// How a leaf is chosen for refinement.
#[derive(Clone, Copy, Debug)]
enum Marking {
    /// Largest `defect × h`. The registered adaptive strategy.
    Greedy,
    /// A uniformly random active leaf, seeded by the contained value. The
    /// control that isolates the indicator from the octree.
    Random(u64),
}

/// Heap entry: indicator bits, then the origin reversed so ties break on the
/// smallest integer origin, then the node index.
type Cand = (u64, Reverse<[u32; 3]>, u32);

/// One octree cell.
#[derive(Clone, Copy, Debug)]
struct Node {
    /// Minimum corner, in integer lattice units.
    origin: [u32; 3],
    /// Edge length, in integer lattice units. Always even.
    size: u32,
    /// Levels below the root grid.
    depth: u32,
    /// Field values at the eight corners, indexed `x + 2y + 4z`.
    corners: [f64; 8],
    /// Index of the first of eight children, or `u32::MAX` for a leaf.
    children: u32,
    /// Whether the cell can carry a band probe, and so may be refined.
    active: bool,
    /// `defect × h`, the refinement indicator.
    indicator: f64,
}

/// The octree, its sample cache, and the field it is built over.
struct Tree<'a, F> {
    field: &'a F,
    dom: Domain,
    nodes: Vec<Node>,
    /// Root-grid cell `bx + BASE_CELLS·(by + BASE_CELLS·bz)` to node index.
    roots: Vec<u32>,
    /// Every field evaluation, keyed on its integer lattice point. The `bool` is
    /// "this point is a corner of some cell", which separates the
    /// reconstruction's own degrees of freedom from the estimator's.
    cache: BTreeMap<[u32; 3], (f64, bool)>,
    corner_count: usize,
}

impl<'a, F: Sdf<Scalar = f64>> Tree<'a, F> {
    /// The root grid: `BASE_CELLS³` cells over the whole domain.
    fn new(field: &'a F, dom: Domain) -> Self {
        let cells = (BASE_CELLS * BASE_CELLS * BASE_CELLS) as usize;
        let mut t = Self {
            field,
            dom,
            nodes: Vec::new(),
            roots: vec![0; cells],
            cache: BTreeMap::new(),
            corner_count: 0,
        };
        let step = UNITS / BASE_CELLS;
        for bz in 0..BASE_CELLS {
            for by in 0..BASE_CELLS {
                for bx in 0..BASE_CELLS {
                    let origin = [bx * step, by * step, bz * step];
                    let corners: [f64; 8] = std::array::from_fn(|i| {
                        let i = i as u32;
                        t.sample(
                            [
                                origin[0] + (i & 1) * step,
                                origin[1] + ((i >> 1) & 1) * step,
                                origin[2] + ((i >> 2) & 1) * step,
                            ],
                            true,
                        )
                    });
                    let idx = t.create(origin, step, 0, corners);
                    t.roots[(bx + BASE_CELLS * (by + BASE_CELLS * bz)) as usize] = idx;
                }
            }
        }
        t
    }

    /// Field evaluations so far. This is the arm's `N`.
    fn dof(&self) -> usize {
        self.cache.len()
    }

    /// Evaluations that are a cell corner, i.e. the reconstruction's own degrees
    /// of freedom without the estimator's.
    fn corner_dof(&self) -> usize {
        self.corner_count
    }

    /// Sample one lattice point, memoised.
    fn sample(&mut self, k: [u32; 3], corner: bool) -> f64 {
        if let Some(entry) = self.cache.get_mut(&k) {
            if corner && !entry.1 {
                entry.1 = true;
                self.corner_count += 1;
            }
            return entry.0;
        }
        let p = [
            self.dom.lo[0] + self.dom.span * f64::from(k[0]) / UNITS_F,
            self.dom.lo[1] + self.dom.span * f64::from(k[1]) / UNITS_F,
            self.dom.lo[2] + self.dom.span * f64::from(k[2]) / UNITS_F,
        ];
        let v = self.field.sample(p);
        self.cache.insert(k, (v, corner));
        if corner {
            self.corner_count += 1;
        }
        v
    }

    /// Push a node, computing its defect, its activity and its indicator.
    fn create(&mut self, origin: [u32; 3], size: u32, depth: u32, corners: [f64; 8]) -> u32 {
        let half = size / 2;
        let centre = self.sample(
            [origin[0] + half, origin[1] + half, origin[2] + half],
            false,
        );
        let defect = (centre - trilinear(&corners, 0.5, 0.5, 0.5)).abs();

        let mut min_v = centre;
        let mut max_v = centre;
        let mut min_abs = centre.abs();
        let mut any_neg = false;
        let mut any_pos = false;
        for &c in &corners {
            if c < min_v {
                min_v = c;
            }
            if c > max_v {
                max_v = c;
            }
            if c.abs() < min_abs {
                min_abs = c.abs();
            }
            if c < 0.0 {
                any_neg = true;
            } else {
                any_pos = true;
            }
        }
        let spread = max_v - min_v;
        let active = (any_neg && any_pos) || min_abs <= self.dom.band + ACTIVITY_SLACK * spread;

        let h = self.dom.span * f64::from(size) / UNITS_F;
        let idx = self.nodes.len() as u32;
        self.nodes.push(Node {
            origin,
            size,
            depth,
            corners,
            children: u32::MAX,
            active,
            indicator: defect * h,
        });
        idx
    }

    /// Offer a leaf to the marking strategy, if it is one the strategy may pick.
    fn enqueue(&self, idx: u32, heap: &mut BinaryHeap<Cand>, pool: &mut Vec<u32>, mark: Marking) {
        let node = self.nodes[idx as usize];
        if !node.active || node.depth >= MAX_EXTRA_DEPTH {
            return;
        }
        match mark {
            Marking::Greedy => heap.push((node.indicator.to_bits(), Reverse(node.origin), idx)),
            Marking::Random(_) => pool.push(idx),
        }
    }

    /// Split one leaf into eight, sampling the new corners and the new centres.
    fn refine(
        &mut self,
        idx: u32,
        heap: &mut BinaryHeap<Cand>,
        pool: &mut Vec<u32>,
        mark: Marking,
    ) {
        let node = self.nodes[idx as usize];
        let half = node.size / 2;
        let first = self.nodes.len() as u32;
        for oct in 0..8u32 {
            let corg = [
                node.origin[0] + (oct & 1) * half,
                node.origin[1] + ((oct >> 1) & 1) * half,
                node.origin[2] + ((oct >> 2) & 1) * half,
            ];
            let corners: [f64; 8] = std::array::from_fn(|i| {
                let i = i as u32;
                self.sample(
                    [
                        corg[0] + (i & 1) * half,
                        corg[1] + ((i >> 1) & 1) * half,
                        corg[2] + ((i >> 2) & 1) * half,
                    ],
                    true,
                )
            });
            let child = self.create(corg, half, node.depth + 1, corners);
            assert_eq!(
                child,
                first + oct,
                "the eight children of a leaf must be consecutive, or `children` is a lie"
            );
        }
        self.nodes[idx as usize].children = first;
        for oct in 0..8u32 {
            self.enqueue(first + oct, heap, pool, mark);
        }
    }

    /// A world point in fractional lattice units, clamped to the domain.
    fn units_of(&self, p: [f64; 3]) -> [f64; 3] {
        std::array::from_fn(|a| {
            ((p[a] - self.dom.lo[a]) / self.dom.span * UNITS_F).clamp(0.0, UNITS_F)
        })
    }

    /// The leaf containing a point given in fractional lattice units.
    fn leaf_of(&self, t: [f64; 3]) -> u32 {
        let base = f64::from(UNITS / BASE_CELLS);
        let bi: [u32; 3] = std::array::from_fn(|a| {
            (t[a] / base).floor().clamp(0.0, f64::from(BASE_CELLS - 1)) as u32
        });
        let mut idx = self.roots[(bi[0] + BASE_CELLS * (bi[1] + BASE_CELLS * bi[2])) as usize];
        loop {
            let node = self.nodes[idx as usize];
            if node.children == u32::MAX {
                return idx;
            }
            let half = f64::from(node.size / 2);
            let oct = u32::from(t[0] - f64::from(node.origin[0]) >= half)
                | (u32::from(t[1] - f64::from(node.origin[1]) >= half) << 1)
                | (u32::from(t[2] - f64::from(node.origin[2]) >= half) << 2);
            idx = node.children + oct;
        }
    }

    /// Root-mean-square defect, sup defect, and probes whose leaf the strategy
    /// declined to refine.
    fn probe_error(&self, probes: &Probes) -> (f64, f64, usize) {
        let mut sum = 0.0;
        let mut sup = 0.0;
        let mut inactive = 0usize;
        for (p, &tv) in probes.points.iter().zip(&probes.truth) {
            let t = self.units_of(*p);
            let node = self.nodes[self.leaf_of(t) as usize];
            if !node.active {
                inactive += 1;
            }
            let s = f64::from(node.size);
            let l: [f64; 3] =
                std::array::from_fn(|a| ((t[a] - f64::from(node.origin[a])) / s).clamp(0.0, 1.0));
            let e = (trilinear(&node.corners, l[0], l[1], l[2]) - tv).abs();
            sum += e * e;
            if e > sup {
                sup = e;
            }
        }
        ((sum / probes.points.len() as f64).sqrt(), sup, inactive)
    }

    /// Leaves, and the deepest level any leaf reached.
    fn shape(&self) -> (usize, u32) {
        let mut leaves = 0usize;
        let mut deepest = 0u32;
        for node in &self.nodes {
            if node.children == u32::MAX {
                leaves += 1;
                if node.depth > deepest {
                    deepest = node.depth;
                }
            }
        }
        (leaves, deepest)
    }
}

/// One adaptive sweep, snapshotted at each evaluation budget.
#[derive(Clone, Debug, Default)]
struct Adaptive {
    /// Field evaluations at each snapshot.
    dof: Vec<u64>,
    /// Reconstruction corners at each snapshot, the estimator excluded.
    corner_dof: Vec<u64>,
    /// Root-mean-square defect at each snapshot.
    rms: Vec<f64>,
    /// Sup defect at each snapshot.
    sup: Vec<f64>,
    /// Leaves at the last snapshot.
    leaves: usize,
    /// Deepest leaf at the last snapshot.
    max_depth: u32,
    /// The marking ran out of refinable leaves before the last budget.
    truncated: bool,
    /// Worst per-snapshot count of probes sitting in a leaf the strategy would
    /// not refine. Asserted zero for the registered arm.
    probes_inactive: usize,
}

/// Refine one octree through the budget ladder.
fn adaptive_run<F: Sdf<Scalar = f64>>(
    field: &F,
    dom: Domain,
    probes: &Probes,
    budgets: &[u64],
    mark: Marking,
) -> Adaptive {
    let mut tree = Tree::new(field, dom);
    let mut heap: BinaryHeap<Cand> = BinaryHeap::new();
    let mut pool: Vec<u32> = Vec::new();
    let mut rng = Rng::new(match mark {
        Marking::Greedy => 0,
        Marking::Random(s) => s,
    });
    for idx in 0..tree.nodes.len() as u32 {
        tree.enqueue(idx, &mut heap, &mut pool, mark);
    }

    let mut out = Adaptive::default();
    for &budget in budgets {
        while (tree.dof() as u64) < budget {
            let pick = match mark {
                Marking::Greedy => heap.pop().map(|c| c.2),
                Marking::Random(_) => {
                    if pool.is_empty() {
                        None
                    } else {
                        let i = (rng.next_u64() % pool.len() as u64) as usize;
                        Some(pool.swap_remove(i))
                    }
                }
            };
            let Some(idx) = pick else {
                out.truncated = true;
                break;
            };
            tree.refine(idx, &mut heap, &mut pool, mark);
        }
        if out.truncated {
            break;
        }
        let (rms, sup, inactive) = tree.probe_error(probes);
        let (leaves, deepest) = tree.shape();
        out.dof.push(tree.dof() as u64);
        out.corner_dof.push(tree.corner_dof() as u64);
        out.rms.push(rms);
        out.sup.push(sup);
        out.leaves = leaves;
        out.max_depth = deepest;
        out.probes_inactive = out.probes_inactive.max(inactive);
    }
    out
}

// ─── the fit ─────────────────────────────────────────────────────────────────

/// One ordinary-least-squares fit of `ln(error)` on `ln(dof)`.
#[derive(Clone, Copy, Debug)]
struct Fit {
    /// `s`, the negated slope. `error ~ N^(-s)`.
    s: f64,
    /// Half-width of the two-sided 95% interval on `s`, from the residual
    /// standard error and Student-t at `points − 2` degrees of freedom.
    half_width: f64,
    /// Coefficient of determination, which is how a transient shows up.
    r2: f64,
    /// Points the fit used.
    points: usize,
}

/// Fit `s` over a series, or `None` when the series cannot carry a slope.
///
/// Refuses fewer than four points — the registration's floor — and refuses a
/// non-positive error, because `ln 0` is not a number and a zero error is a
/// resolution at which the reconstruction is *exact* rather than a data point on
/// a decay law. Both refusals surface as `c1_field = false`, which is C1's own
/// registered failure mode rather than a silence.
fn fit(dof: &[u64], err: &[f64]) -> Option<Fit> {
    let k = dof.len();
    if k < 4 || err.len() != k {
        return None;
    }
    if !err.iter().all(|e| e.is_finite() && *e > 0.0) {
        return None;
    }
    if !dof.iter().all(|d| *d > 0) {
        return None;
    }
    let xs: Vec<f64> = dof.iter().map(|d| (*d as f64).ln()).collect();
    let ys: Vec<f64> = err.iter().map(|e| e.ln()).collect();
    let kf = k as f64;
    let xbar = xs.iter().sum::<f64>() / kf;
    let ybar = ys.iter().sum::<f64>() / kf;
    let sxx: f64 = xs.iter().map(|x| (x - xbar) * (x - xbar)).sum();
    let syy: f64 = ys.iter().map(|y| (y - ybar) * (y - ybar)).sum();
    let sxy: f64 = xs
        .iter()
        .zip(&ys)
        .map(|(x, y)| (x - xbar) * (y - ybar))
        .sum();
    if !sxx.is_finite() || sxx <= 0.0 {
        return None;
    }
    let slope = sxy / sxx;
    let sse = (syy - slope * sxy).max(0.0);
    let df = k - 2;
    let se = (sse / df as f64 / sxx).sqrt();
    let t = T_975[df.min(T_975.len() - 1)];
    Some(Fit {
        s: -slope,
        half_width: t * se,
        r2: if syy > 0.0 { 1.0 - sse / syy } else { 0.0 },
        points: k,
    })
}

/// The slope over the last three points, with no interval.
///
/// A diagnostic, not a verdict: with three points the interval carries
/// `t = 12.706` and says nothing. What it does say is whether the six-point
/// slope is being dragged by the coarse end of the sweep.
fn tail_slope(dof: &[u64], err: &[f64]) -> f64 {
    let k = dof.len();
    if k < 3 || err.len() != k {
        return f64::NAN;
    }
    let d = &dof[k - 3..];
    let e = &err[k - 3..];
    if !e.iter().all(|v| v.is_finite() && *v > 0.0) {
        return f64::NAN;
    }
    let xs: Vec<f64> = d.iter().map(|v| (*v as f64).ln()).collect();
    let ys: Vec<f64> = e.iter().map(|v| v.ln()).collect();
    let n = 3.0;
    let xbar = xs.iter().sum::<f64>() / n;
    let ybar = ys.iter().sum::<f64>() / n;
    let sxx: f64 = xs.iter().map(|x| (x - xbar) * (x - xbar)).sum();
    let sxy: f64 = xs
        .iter()
        .zip(&ys)
        .map(|(x, y)| (x - xbar) * (y - ybar))
        .sum();
    if sxx <= 0.0 {
        return f64::NAN;
    }
    -(sxy / sxx)
}

// ─── the mesh control ────────────────────────────────────────────────────────

/// `validate::accuracy`'s symmetric Hausdorff against triangle count, or the
/// reason the field was skipped.
///
/// **Gated on `bound().is_exact()`**, because `|sample|` is not a distance on
/// the other four and `fields/mod.rs:104-109` says so: for `Unbounded` "no
/// accuracy figure against `|sample|` means anything". The skip is a recorded
/// column, not an omission.
fn mesh_series<F>(field: &F) -> Result<(Vec<u64>, Vec<f64>), &'static str>
where
    F: ReferenceField + Sdf<Scalar = f64>,
{
    match field.bound() {
        FieldBound::Exact => {}
        FieldBound::Lipschitz { .. } => return Err("bound=Lipschitz"),
        FieldBound::Underestimate { .. } => return Err("bound=Underestimate"),
        FieldBound::Unbounded => return Err("bound=Unbounded"),
    }
    let mut dof = Vec::with_capacity(MESH_SAMPLES.len());
    let mut err = Vec::with_capacity(MESH_SAMPLES.len());
    for &n in &MESH_SAMPLES {
        let (shape, origin, h) = common::grid::<f64, _>(field, n);
        let mut mesh = MeshBuffer::<f64>::new();
        let mut mc = MarchingCubes::<f64>::new();
        mc.extract(field, &shape, origin, h, &mut mesh)
            .expect("marching cubes on a reference field grid");
        let cfg = AccuracyConfig::from_cell_size(h).expect("a positive cell size");
        let report = accuracy(&mesh.positions, &mesh.indices, field, &shape, origin, &cfg)
            .expect("accuracy on an exact field");
        dof.push(mesh.triangle_count() as u64);
        err.push(report.symmetric_hausdorff());
    }
    Ok((dof, err))
}

// ─── one field ───────────────────────────────────────────────────────────────

/// Everything measured for one field, before any verdict is taken.
#[derive(Clone, Debug)]
struct Row {
    field: &'static str,
    seed: u64,
    band: f64,
    band_volume_fraction: f64,
    probes: usize,
    dof_uniform: Vec<u64>,
    err_uniform: Vec<f64>,
    sup_uniform: Vec<f64>,
    greedy: Adaptive,
    random: Adaptive,
    fit_uniform: Option<Fit>,
    fit_adaptive: Option<Fit>,
    fit_random: Option<Fit>,
    fit_uniform_sup: Option<Fit>,
    fit_adaptive_sup: Option<Fit>,
    fit_mesh: Option<Fit>,
    mesh_dof: Vec<u64>,
    mesh_err: Vec<f64>,
    mesh_skip: &'static str,
    wall_ns: f64,
}

/// Measure one field: the probe set, both families, all three controls.
fn measure<F>(name: &'static str, field: &F) -> Row
where
    F: ReferenceField + Sdf<Scalar = f64>,
{
    let started = Instant::now();
    let (lo, hi) = field.domain();
    let span = hi[0] - lo[0];
    assert!(
        (hi[1] - lo[1] - span).abs() < 1e-12 && (hi[2] - lo[2] - span).abs() < 1e-12,
        "{name}: P-161 measures over a cubic domain and this one is not cubic"
    );
    let dom = Domain {
        lo,
        span,
        band: span * BAND_FRACTION,
    };
    let seed = seed_for(name);
    let probes = probe_set(field, dom, seed);

    let budgets: Vec<u64> = UNIFORM_SAMPLES
        .iter()
        .map(|&n| u64::from(n).pow(3))
        .collect();

    let mut dof_uniform = Vec::with_capacity(UNIFORM_SAMPLES.len());
    let mut err_uniform = Vec::with_capacity(UNIFORM_SAMPLES.len());
    let mut sup_uniform = Vec::with_capacity(UNIFORM_SAMPLES.len());
    for (&n, &budget) in UNIFORM_SAMPLES.iter().zip(&budgets) {
        let (rms, sup) = uniform_error(field, dom, n, &probes);
        dof_uniform.push(budget);
        err_uniform.push(rms);
        sup_uniform.push(sup);
    }

    let greedy = adaptive_run(field, dom, &probes, &budgets, Marking::Greedy);
    let random = adaptive_run(
        field,
        dom,
        &probes,
        &budgets,
        Marking::Random(seed ^ RANDOM_SALT),
    );

    let (mesh_dof, mesh_err, mesh_skip) = match mesh_series(field) {
        Ok((d, e)) => (d, e, "none"),
        Err(reason) => (Vec::new(), Vec::new(), reason),
    };

    Row {
        field: name,
        seed,
        band: dom.band,
        band_volume_fraction: probes.band_volume_fraction(),
        probes: probes.points.len(),
        fit_uniform: fit(&dof_uniform, &err_uniform),
        fit_uniform_sup: fit(&dof_uniform, &sup_uniform),
        fit_adaptive: fit(&greedy.dof, &greedy.rms),
        fit_adaptive_sup: fit(&greedy.dof, &greedy.sup),
        fit_random: fit(&random.dof, &random.rms),
        fit_mesh: fit(&mesh_dof, &mesh_err),
        dof_uniform,
        err_uniform,
        sup_uniform,
        greedy,
        random,
        mesh_dof,
        mesh_err,
        mesh_skip,
        wall_ns: started.elapsed().as_secs_f64() * 1e9,
    }
}

// ─── verdicts ────────────────────────────────────────────────────────────────

/// What one field's two fits say.
#[derive(Clone, Copy, Debug)]
struct Verdict {
    /// `s` was estimable for both families: both fits exist and neither interval
    /// reaches zero.
    estimable: bool,
    /// `s_adaptive − s_uniform`.
    difference: f64,
    /// Combined half-width, `sqrt(hw_u² + hw_a²)`.
    combined_half_width: f64,
    /// The difference clears its own combined interval, so LOD buys a rate here.
    helps: bool,
}

fn verdict_of(row: &Row) -> Verdict {
    match (row.fit_uniform, row.fit_adaptive) {
        (Some(u), Some(a)) => {
            let combined = (u.half_width * u.half_width + a.half_width * a.half_width).sqrt();
            let difference = a.s - u.s;
            Verdict {
                estimable: u.s - u.half_width > 0.0
                    && a.s - a.half_width > 0.0
                    && u.half_width.is_finite()
                    && a.half_width.is_finite(),
                difference,
                combined_half_width: combined,
                helps: difference > combined && combined.is_finite(),
            }
        }
        _ => Verdict {
            estimable: false,
            difference: f64::NAN,
            combined_half_width: f64::NAN,
            helps: false,
        },
    }
}

// ─── formatting ──────────────────────────────────────────────────────────────

fn num(x: f64) -> String {
    format!("{x:.6}")
}

fn join_u64(v: &[u64]) -> String {
    v.iter().map(u64::to_string).collect::<Vec<_>>().join("|")
}

fn join_sci(v: &[f64]) -> String {
    v.iter()
        .map(|x| format!("{x:.4e}"))
        .collect::<Vec<_>>()
        .join("|")
}

fn s_of(f: Option<Fit>) -> f64 {
    f.map_or(f64::NAN, |f| f.s)
}

fn interval(f: Option<Fit>) -> String {
    f.map_or_else(
        || String::from("NaN|NaN"),
        |f| format!("{:.6}|{:.6}", f.s - f.half_width, f.s + f.half_width),
    )
}

fn r2_of(f: Option<Fit>) -> f64 {
    f.map_or(f64::NAN, |f| f.r2)
}

fn points_of(f: Option<Fit>) -> usize {
    f.map_or(0, |f| f.points)
}

// ─── main ────────────────────────────────────────────────────────────────────

fn main() {
    if !std::env::args().any(|arg| arg == "--bench") {
        return;
    }
    let prereg = isomesh::experiment!("P-161");

    let mut rows: Vec<Row> = Vec::new();
    isomesh::for_each_reference_field!(f64, |name, field| {
        // Inline block, so no `return` in here (M-253).
        rows.push(measure(name, &field));
    });

    common::experiment::run(prereg, |run| {
        // ── vacuity controls ─────────────────────────────────────────────────

        assert_eq!(
            rows.len(),
            8,
            "VOID: C1 claims all eight reference fields and only {} were measured, so the clause \
             would be scored over a subset chosen after the fact",
            rows.len()
        );

        for row in &rows {
            assert_eq!(
                row.probes, PROBES,
                "VOID: {} filled only {} of {PROBES} band probes in {PROBE_ATTEMPTS} attempts, so \
                 every error on this field is a mean over a population the harness could not \
                 assemble (M-44)",
                row.field, row.probes
            );
            assert!(
                row.dof_uniform.len() >= 4,
                "VOID: {} has {} uniform resolutions and the registration's floor is four, so no \
                 slope may be fitted",
                row.field,
                row.dof_uniform.len()
            );
            assert!(
                !row.greedy.truncated && row.greedy.dof.len() >= 4,
                "VOID: {}'s adaptive sweep stopped after {} of {} budgets, so s_adaptive would be \
                 a property of MAX_EXTRA_DEPTH = {MAX_EXTRA_DEPTH} rather than of the field",
                row.field,
                row.greedy.dof.len(),
                UNIFORM_SAMPLES.len()
            );
            let window = (*row.dof_uniform.last().expect("a non-empty sweep") as f64
                / row.dof_uniform[0] as f64)
                .log10();
            assert!(
                window >= 1.0,
                "VOID: {}'s degrees of freedom span only {window:.3} decades, and a slope fitted \
                 over less than one decade cannot distinguish C1's registered failure -- s not \
                 estimable from the achievable range -- from a short sweep",
                row.field
            );
            assert_eq!(
                row.greedy.probes_inactive, 0,
                "VOID: {} put {} probes in leaves the adaptive criterion declined to refine, so \
                 s_adaptive is partly a property of the activity test rather than of the field",
                row.field, row.greedy.probes_inactive
            );
            assert!(
                row.greedy.max_depth >= MIN_ADAPTIVE_DEPTH,
                "VOID: {}'s adaptive octree reached depth {} below its root grid, so it never \
                 varied its sample density and both arms are uniform (P-160's control)",
                row.field,
                row.greedy.max_depth
            );
            for (u, a) in row.dof_uniform.iter().zip(&row.greedy.dof) {
                assert!(
                    *a >= *u && (*a as f64) <= BUDGET_OVERSHOOT * (*u as f64),
                    "VOID: {}'s adaptive snapshot spent {a} evaluations against a uniform budget \
                     of {u}, so the two slopes are fitted over different windows and their \
                     difference is an artefact of the ladder",
                    row.field
                );
            }
            if let Some(f) = row.fit_uniform {
                assert!(
                    f.half_width.is_finite(),
                    "VOID: {}'s s_uniform carries no confidence interval, so the C2 split would be \
                     scored against zero rather than against noise",
                    row.field
                );
            }
            if let Some(f) = row.fit_adaptive {
                assert!(
                    f.half_width.is_finite(),
                    "VOID: {}'s s_adaptive carries no confidence interval, so the C2 split would \
                     be scored against zero rather than against noise",
                    row.field
                );
            }
        }

        assert!(
            rows.iter()
                .any(|r| r.fit_uniform.is_some() && r.fit_adaptive.is_some()),
            "VOID: not one field produced both fits, so neither clause is measured and the file \
             would report eight unfitted rows as if that were a result"
        );

        // ── verdicts: global before per-row ──────────────────────────────────

        let verdicts: Vec<Verdict> = rows.iter().map(verdict_of).collect();
        let c1 = verdicts.iter().all(|v| v.estimable);
        let gaining = verdicts.iter().filter(|v| v.helps).count();
        let not_gaining = verdicts.len() - gaining;
        let c2 = gaining >= 1 && not_gaining >= 1;

        for (row, v) in rows.iter().zip(&verdicts) {
            let su = s_of(row.fit_uniform);
            let sa = s_of(row.fit_adaptive);
            // `A^s` is defined for the adaptive family, which is the object the
            // cited papers name; the uniform class is beside it as an extra.
            let class = row
                .fit_adaptive
                .map_or_else(|| String::from("unfitted"), |f| format!("A^{:.2}", f.s));
            let class_uniform = row
                .fit_uniform
                .map_or_else(|| String::from("unfitted"), |f| format!("A^{:.2}", f.s));
            // The factor the best-N-term error falls by when the evaluation
            // budget doubles.
            let decay = (-sa).exp2();
            // Degrees of freedom the adaptive family saves per decade of error
            // reduction. Exactly 1.0 when s_difference is 0 -- the registration's
            // "octree LOD provably buys nothing asymptotically on that field".
            let gain = 10f64.powf(1.0 / su - 1.0 / sa);
            let diff_lo = v.difference - v.combined_half_width;
            let diff_hi = v.difference + v.combined_half_width;

            run.record(&[
                ("field", row.field.to_string()),
                ("s_uniform", num(su)),
                ("s_adaptive", num(sa)),
                ("s_difference", num(v.difference)),
                ("class_membership", class),
                ("n_term_decay", num(decay)),
                ("lod_asymptotic_gain", num(gain)),
                ("c1_holds", c1.to_string()),
                ("c2_holds", c2.to_string()),
                // ── extras (M-273) ──
                ("adaptive_leaves", row.greedy.leaves.to_string()),
                ("adaptive_max_depth", row.greedy.max_depth.to_string()),
                ("adaptive_truncated", row.greedy.truncated.to_string()),
                ("band_half_width", num(row.band)),
                ("band_volume_fraction", num(row.band_volume_fraction)),
                ("c1_field", v.estimable.to_string()),
                ("class_uniform", class_uniform),
                ("dof_adaptive", join_u64(&row.greedy.dof)),
                ("dof_adaptive_corners", join_u64(&row.greedy.corner_dof)),
                ("dof_random", join_u64(&row.random.dof)),
                ("dof_uniform", join_u64(&row.dof_uniform)),
                ("err_adaptive", join_sci(&row.greedy.rms)),
                ("err_random", join_sci(&row.random.rms)),
                ("err_uniform", join_sci(&row.err_uniform)),
                ("fields_gaining", gaining.to_string()),
                ("fields_not_gaining", not_gaining.to_string()),
                (
                    "fit_points_adaptive",
                    points_of(row.fit_adaptive).to_string(),
                ),
                ("fit_points_uniform", points_of(row.fit_uniform).to_string()),
                ("fit_r2_adaptive", num(r2_of(row.fit_adaptive))),
                ("fit_r2_uniform", num(r2_of(row.fit_uniform))),
                ("lod_helps", v.helps.to_string()),
                (
                    "mesh_accuracy_skipped",
                    (row.mesh_skip != "none").to_string(),
                ),
                ("mesh_accuracy_skip_reason", row.mesh_skip.to_string()),
                ("mesh_dof", join_u64(&row.mesh_dof)),
                ("mesh_hausdorff", join_sci(&row.mesh_err)),
                ("mesh_s_uniform", num(s_of(row.fit_mesh))),
                ("mesh_s_uniform_ci", interval(row.fit_mesh)),
                ("probes", row.probes.to_string()),
                (
                    "probes_in_inactive_leaves",
                    row.greedy.probes_inactive.to_string(),
                ),
                ("random_truncated", row.random.truncated.to_string()),
                ("resolutions", row.dof_uniform.len().to_string()),
                ("s_adaptive_ci", interval(row.fit_adaptive)),
                ("s_adaptive_sup", num(s_of(row.fit_adaptive_sup))),
                (
                    "s_adaptive_tail",
                    num(tail_slope(&row.greedy.dof, &row.greedy.rms)),
                ),
                ("s_difference_ci", format!("{diff_lo:.6}|{diff_hi:.6}")),
                ("s_random", num(s_of(row.fit_random))),
                ("s_random_ci", interval(row.fit_random)),
                ("s_uniform_ci", interval(row.fit_uniform)),
                ("s_uniform_sup", num(s_of(row.fit_uniform_sup))),
                (
                    "s_uniform_tail",
                    num(tail_slope(&row.dof_uniform, &row.err_uniform)),
                ),
                ("seed", format!("{:016x}", row.seed)),
                ("sup_adaptive", join_sci(&row.greedy.sup)),
                ("sup_uniform", join_sci(&row.sup_uniform)),
                ("wall_ns", num(row.wall_ns)),
            ]);
        }
    });
}
