//! **P-124 — the monotone-edge condition on the ambient complex.**
//!
//! Ticket: R-052. Pre-registered before this harness existed.
//!
//! ```bash
//! cargo bench --bench experiment_p124
//! ```
//!
//! Writes `docs/experiments/p-124.csv`.
//!
//! # What was missing
//!
//! **`✗36 / M-351` closed `P-55` by *proof*, and the proof was about the complex
//! rather than about the fixture** (`FINDINGS.md:8954`). A mesh edge joins two
//! vertices that both lie on the extracted zero set, so it is a **chord**; for a
//! strictly convex surface the chord's interior is strictly inside the solid, `f`
//! runs `0 → negative → 0`, and the directional derivative *must* reverse. On a
//! sphere of radius `r` with central angle `θ` the two endpoint derivatives are
//! closed-form — `g(0) = r(cos θ − 1) < 0`, `g(1) = r(1 − cos θ) > 0` — for
//! **every** distinct pair of surface points. `sphere` duly read `1000.000` per
//! 1k at 17³ (804 of 804), and **97.5%** of that sweep's 901,583 flags were
//! decided by `g(0)` and `g(1)` before any interior sample existed. No `k`, no
//! `w` and no tolerance reaches that.
//!
//! Worse, the tolerance was scaled by `|f(a)| + |f(b)|` — the residual at two
//! vertices that are *on the zero set*. `max_abs_tolerance` was **exactly 0.0**
//! on `box_exact` at all four resolutions and never exceeded `7.63e-13`
//! anywhere, against a minimum deciding reversal of `1.408e-3`: **nine orders of
//! magnitude of inert guard**, and identical counts at `1e-14`/`1e-12`/`1e-10` on
//! 31 of 32 rows.
//!
//! `✗36` named the non-vacuous port and did not run it: *"the paper's PL
//! function lives on the **ambient** complex, and monotonicity is a condition on
//! *that* complex's edges — whose endpoints are not on the zero set and for which
//! `|f(a)| + |f(b)|` is a real scale"*. **Nothing in this repository has ever
//! evaluated the predicate there.** That is this row, and it is the only thing
//! this row is: it decides whether the ambient reading of Finken et al. is
//! measurable at all. It proposes no source change and moves no time.
//!
//! # What the theorem actually says, quoted rather than paraphrased
//!
//! Finken, Li, Wang, Guo & Levine, *Topology-Preserving Meshing of Implicit
//! Scalar Fields via Monotonicity Constraints*, arXiv:2608.12142 §3.1, corpus
//! `doc_id` `10.48550_arXiv.2608.12142`. **Theorem 1.** Let `f̂` be a PL function
//! on a mesh monotonic with respect to a Morse `f`. Then **(1)** every critical
//! point of `f̂` coincides with a critical point of `f`, and **(2)** any critical
//! point of `f` either coincides with a critical point of `f̂` **or shares a
//! triangle with at least one other critical point of `f`**.
//!
//! Part 2 splits by critical-point type. The extremum half: expanding closed
//! isocontours *"must eventually intersect an edge. Continuing to increase the
//! value forces multiple intersections with the same edge, violating
//! monotonicity"*. **Part 2b, the saddle half — the one C1 is gated on:** *"the
//! four contour branches expand outward from `x` … Since a triangle has only
//! three edges, at least two branches must intersect the same edge, again
//! violating monotonicity."*
//!
//! So the **contrapositive of Part 2b** is: a simplex whose every edge is
//! monotone contains **no isolated interior critical point** of `f` — any
//! critical point inside it is paired, hence *zero or at least two*. That is C1's
//! sentence, and it is why one non-monotone edge is enough to reject a cell: the
//! certificate is a property of the whole simplex boundary and is destroyed by a
//! single reversal.
//!
//! Two limitations carried over from `✗36` unchanged, because they are the
//! paper's own. **The proof is 2D** and Part 2b's pigeonhole step is literally
//! 2D-combinatorial (`docs/research/2026-08-23-phase-20-source-corrections.md:139-142`);
//! a tetrahedron has four faces and six edges, and four branches into four faces
//! does not pigeonhole. **It does not apply to the trilinear interpolant**, which
//! is what `marching_cubes` contours. Everything below is therefore a **labelled
//! 3D port** on a *simplicial* complex — the one setting where the theorem's
//! hypothesis (`f̂` PL, critical points only at vertices) is actually met — and a
//! clause that holds is evidence about this crate's grids, not a transported
//! proof.
//!
//! # The complex, and why it is the shipped one
//!
//! **Kuhn's / Freudenthal's six tetrahedra per cell**, taken from the crate
//! rather than re-derived: `isomesh::marching_tetrahedra::table::TETS`
//! (`table.rs:87`) and `TET_EDGES` (`table.rs:121`) are `pub`, so the complex
//! this harness sweeps is bit-for-bit the complex `MarchingTetrahedra` marches.
//! Three properties of it are load-bearing here and all three are already
//! proved:
//!
//! - **It tiles without alternation.** Each cube face carries one diagonal fixed
//!   by the corner numbering alone — `−x: 0–6`, `+x: 1–7`, `−y: 0–5`, `+y: 2–7`,
//!   `−z: 0–3`, `+z: 4–7` (`table.rs:36-53`) — and two cells adjacent along an
//!   axis split their shared face on **the same two lattice points**. So "the tet
//!   edge `{v, v+δ}`" is well defined without asking which cell named it.
//! - **`P-100` measured that, not just argued it**: `✗78 / M-412`, `open_edges`
//!   **0 on 80 of 80** rows across a chunk seam, with a mismatched-diagonal
//!   control reading **80–3,765** on the same rows (`FINDINGS.md:21609`,
//!   `:21770`). A zero seam count there is a measurement.
//! - **Every tet edge steps `0` or `1` per axis, never `−1`**: `TETS`' corners
//!   are ordered by inclusion, `TETS[t][0] == 0` and `TETS[t][3] == 7`
//!   (`table.rs:84-86`). `offsets_are_monotone_steps` asserts it on all 36 edge
//!   instances before the sweep starts, so the offset derivation below cannot
//!   silently underflow.
//!
//! From those 36 instances (6 tets × 6 edges) the harness *derives* — does not
//! transcribe — the **19 distinct cell-local edges**: **12** axis edges (all
//! twelve cube edges), **6** face diagonals (one per face) and **1** body
//! diagonal (`0–7`, in every tetrahedron). Reduced to lattice offsets that is
//! exactly **seven** positive `δ`:
//!
//! | class | `δ` | per cell | distinct in an `n³` grid |
//! |---|---|---|---|
//! | axis | `(1,0,0) (0,1,0) (0,0,1)` | 12 | `3n²(n−1)` |
//! | face diagonal | `(0,1,1) (1,0,1) (1,1,0)` | 6 | `3n(n−1)²` |
//! | body diagonal | `(1,1,1)` | 1 | `(n−1)³` |
//!
//! Every distinct tet edge in the lattice is `{v, v+δ}` for exactly one `(v, δ)`,
//! which makes the population **enumerable in closed form with no deduplication
//! structure at all** — `tet_edge_count_matches_closed_form` asserts the swept
//! count against `3n²(n−1) + 3n(n−1)² + (n−1)³` on every row. The vertex figure
//! is the classic Freudenthal **14-neighbourhood** (`±` the seven `δ`), which is
//! what the PL-extremum test below reads.
//!
//! `resolution` counts **samples**, so `n` samples span `n − 1` cells
//! (`benches/common/mod.rs:37-38`); `cells` and `tets = 6 · cells` are that
//! arithmetic and nothing else. `corner_offset` is `pub(crate)`
//! (`crates/isomesh/src/cube.rs:149-155`), so it is **copied into this file with
//! the source line it came from on the row that uses it** —
//! `crates/isomesh/src/**` is read-only for this row and a copy whose line number
//! is in the comment is auditable in a way a `pub` would not be
//! (`experiment_p117.rs:53-56`).
//!
//! # The predicate, and the one number the registration did not fix
//!
//! ```text
//! is_monotone(f, a, b, w):
//!     k    = max(2, ceil(‖b − a‖ / w) + 1)          Finken §4, verbatim
//!     d    = b − a
//!     g(t) = ∇f(a + t·d) · d                        chain rule: d/dt f(a + t·d)
//!     tol  = coef · max(|f(a)|, |f(b)|)             THE REGISTRATION'S SCALE
//!     non-monotone  ⟺  a kept g > 0 and a kept g < 0 both exist
//! ```
//!
//! **`max(|f(a)|, |f(b)|)` is the registration's decision and it is the whole
//! repair of `✗36`.** On the ambient complex the endpoints are grid samples, not
//! surface points, so `|f|` there is a real distance and the guard has a real
//! scale — where on a mesh edge it is the interpolation residual and is the
//! quantity `✗36` proved inert to nine decimal orders. `max` rather than `+` is
//! registered and is also the stronger reading: it is the larger of the two, so
//! it cannot be dragged to zero by one endpoint sitting on the zero set.
//!
//! **The registration fixes the scale and does not fix the coefficient.** That
//! number is not invented here either: it is `P-55`'s, unchanged — `1e-12`, with
//! `1e-14` and `1e-10` recorded beside it as a sensitivity strip, exactly the
//! shape `✗36` asked for and the shape that made the old guard's inertness
//! visible. `nonzero_g_discarded_at_1e12` and `guard_inert` are the columns that
//! say whether it is inert *here*; the registration is not amended to claim
//! otherwise.
//!
//! Endpoint gradients come from the grid cache, so `g(0)` and `g(k−1)` of a
//! shared tet edge are bit-identical from every cell that contains it. Interior
//! samples are fresh `Sdf::gradient` calls. All eight reference fields override
//! `Sdf::gradient` with an **analytic** gradient — `fields/mod.rs` opens with
//! *"the central-difference default is never used by a reference field"* — which
//! `✗36` had to correct in `P-55`'s own registration and which is stated here up
//! front: every `∇f` below is exact, not `O(h²)`.
//!
//! # The arms — one build, one run, one shared population (`M-281`)
//!
//! | arm | complex | edges swept | what it is for |
//! |---|---|---|---|
//! | **ambient** | Kuhn's six tets on the grid | axis **and** face **and** body diagonals | the row |
//! | **mesh-edge control** | the `marching_cubes` surface | chords of the zero set | `✗36`'s reading, reproduced in the same run |
//!
//! Both arms call **the same `examine`** — same `k` rule, same `w`, same
//! tolerance scale, same sign test. The only difference is which segments are fed
//! to it, which is the entire hypothesis.
//!
//! `diagonal_only_failures` is the column that decides whether the ambient
//! reading *bought* anything: **cells rejected by a diagonal with every one of
//! their axis edges monotone**. Face and body diagonals are the part of the
//! complex a chord predicate structurally cannot reach, and axis edges are the
//! ones `marching_cubes` already interpolates along. **If that column reads 0 the
//! ambient reading is the axis-edge reading with extra work**, and the entry says
//! so.
//!
//! # SHARE, recomputed before the numbers
//!
//! **Zero, registered rather than discovered.** This row changes no shipped code
//! path, adds no stage to an extraction and proposes no landing. Nothing in
//! `crates/isomesh/src/` evaluates monotonicity of anything, so there is no total
//! for a fraction of it to be taken from and no Amdahl ceiling to compute. What
//! stands in a share's place is one integer per clause over a denominator that is
//! exact by construction:
//!
//! | clause | quantity | denominator | exact because |
//! |---|---|---|---|
//! | C1 | `non_monotone_cells` | `cells = (n−1)³` | the grid |
//! | C2 | `non_monotone_cells / cells`, ratio per doubling | same | the grid |
//! | C3 | fields with `0 < non_monotone_cells < cells` | 8 | `for_each_reference_field!` |
//!
//! # Which unit carries each verdict, and why none of them is a nanosecond
//!
//! **No clause here is a cost or a wall-clock ratio, so `M-280` and `✗24` do not
//! bite and no `perf_event_open` counter is opened.** C1 and C3 are integer
//! comparisons. C2 is a ratio of two integer *populations* over two exactly
//! enumerated denominators — the good kind: it is computed by one cross-multiplied
//! division of integers-as-`f64`, is identical on any machine at any clock, and
//! cannot be moved by a governor step. This machine spans 1.96–5.62 GHz under
//! `powersave`/`balance_performance`; the two `ns` columns are recorded because
//! they are interesting and **are read by nothing**.
//!
//! # Two places the registration is arithmetically wrong, reported rather than amended
//!
//! `crates/isomesh/src/experiment.rs:27-31` forbids amending a registration to
//! fit the code. Both of these are stated here, measured as registered where that
//! is possible, and recorded in the CSV so the `FINDINGS.md` entry can quote them.
//!
//! **1. C2's literal reading is dead before the run, and it is dead by
//! arithmetic.** The registration says box_exact's non-monotone population *"is
//! `O(n)` and **halves per refinement**, so the count at 33³, 65³ and 129³ must
//! fall by a factor in [1.7, 2.3] per doubling"*. A population that is `O(n)`
//! does not have a **count** that falls — its count *doubles*; it is the
//! *density* that halves. Those are `✗36`'s two separate sentences ("box edges
//! are `O(n)` against `O(n²)` surface cells" and "the one field whose **rate**
//! genuinely halves per refinement") compressed into one, and `✗36`'s own table
//! is the density: **178.1 / 86.1 / 42.3 / 21.0 per 1k**, ratios `2.07 / 2.04 /
//! 2.01`, squarely inside the registered band. The band is a band on the
//! **rate**.
//!
//! Nor does the `O(n)` transfer, and the reason is the interesting half. On the
//! *surface* complex box_exact's non-monotone set is the population of chords
//! straddling a convex **box edge** — a 1D feature in a 2D complex. On the
//! *ambient* complex the field is `‖∇f‖ = 1` and locally affine inside each face's
//! slab, so `g` there is constant and monotone; what is left is the interior
//! **medial axis**, where the exact distance switches face normal — a **2D**
//! feature in a **3D** complex. Different dimensions, *same codimension*, and a
//! codimension-1 feature's density halves per refinement in either. **That is why
//! the band transfers even though the `O(n)` does not**, and it is why the band is
//! two-sided: above `2.3` the feature is thinner than codimension 1 and is
//! vanishing rather than tracking a surface, below `1.7` it is codimension 0 and
//! the predicate is measuring the volume. `box_exact_population_ratio` is the
//! density ratio and carries C2. `box_exact_count_ratio_*` and
//! `box_exact_count_falls` record the literal reading beside it: the count is
//! expected to *rise*, `≈4×` per doubling if the locus is 2D, and the boolean
//! settles it from the file.
//!
//! **2. The vacuity control's sign is inverted, and asserting it as written could
//! never pass.** The registration asks that *"a control arm evaluating the SAME
//! predicate over MESH edges must reproduce `✗36`'s saturation — 0 non-monotone —
//! in the same run"*, one sentence after *"the instrument must be shown able to
//! read a non-monotone edge"*. Those two cannot both be satisfied, and `✗36`'s
//! saturation is **917.6 / 893.7 / 42.3 per 1k non-monotone**, not `0`. The `0`
//! is unreachable for a third, structural reason: on a mesh edge
//! `max(|f(a)|, |f(b)|)` is the interpolation residual, so `tol ≈ 0`, so the guard
//! discards nothing and the chord theorem flags the edge — a `0` there would
//! require the tolerance to swallow a reversal nine orders of magnitude larger
//! than itself.
//!
//! So the harness implements the sentence's **stated purpose** rather than its
//! arithmetic, and both readings are in the file. `mesh_edge_control_non_monotone`
//! is recorded as registered; what is `assert!`ed is that it is **non-zero on
//! every row** and **saturated on `sphere`**. Writing `== 0` instead would abort
//! the first row and produce no CSV, which is not a result — it is a fixture
//! defect, and `✗36` already published the proof that it is one.
//!
//! # The vacuity controls, as asserts rather than as columns
//!
//! `M-44`: a zero that could not have been non-zero is not a measurement. Every
//! zero this harness can report is `assert!`-guarded, so a run that cannot fire
//! aborts instead of recording a pass.
//!
//! | zero at risk | control, asserted | why it licenses the zero |
//! |---|---|---|
//! | `non_monotone_edges`, `non_monotone_cells` | `mesh_edge_control_non_monotone > 0` on **every** row | the same `examine`, in the same run, demonstrably flags a non-monotone edge |
//! | the control itself being the weak one | `2 · control_non_monotone > control_edges` on **`sphere`**, every resolution | `✗36`'s chord theorem is a *theorem*: a strict majority of chords of a strictly convex surface reverse. A control that merely fired would not show the mesh reading is *saturated* rather than *discriminating* |
//! | `diagonal_only_failures` | `tet_edges_face_diagonal + tet_edges_body_diagonal > 0`, and the closed-form equality | a `0` then means the diagonals were **swept and monotone**, not that they were never visited |
//! | `critical_points_lower_bound` | `pl_local_extrema` recorded beside it | separates "no extremum on the grid at all" from "every extremum sits in a rejected star" |
//!
//! `critical_points_lower_bound` is Part **1** of the theorem, not Part 2b: an
//! interior grid vertex that is a **strict** PL local extremum over its
//! 14-neighbour link is a critical point of `f̂` in any dimension (`f̂` is affine on
//! each tet, so a strict link minimum is a strict minimum on every incident tet),
//! and if every one of its fourteen incident tet edges is monotone, Part 1 says it
//! coincides with a critical point of `f`. The column is that count. It **gates
//! nothing** — the 3D transport of Part 1's differentiability step is exactly as
//! unproved as Part 2b's — and it is recorded because a lower bound on `f`'s
//! critical points derived from the census is the thing the census is *for*.
//!
//! # Determinism
//!
//! One thread, no map iteration, no PRNG, `f64` throughout. The sweep order is
//! `offset`-major then `z`, `y`, `x`, fixed. Every gated quantity is an integer,
//! an exact ratio of integers, or a comparison decided by [`f64::total_cmp`] — a
//! total order, so a NaN sorts into view rather than being dropped by a partial
//! comparison.

mod common;

use std::cmp::Ordering;
use std::time::Instant;

use isomesh::fields::ReferenceField;
use isomesh::marching_cubes::MarchingCubes;
use isomesh::marching_tetrahedra::table::{TET_EDGES, TETS};
use isomesh::{MeshBuffer, Sdf};

/// Samples per axis. Exactly the ladder C2 names, and no more: C1 and C3 are
/// evaluated on every one of them.
const RESOLUTIONS: [u32; 3] = [33, 65, 129];

/// The three tolerance coefficients, with `P-55`'s in the middle so the CSV
/// reads as a sensitivity strip. The registration fixes the *scale* the
/// coefficient multiplies and not the coefficient.
const COEFFS: [f64; 3] = [1e-14, 1e-12, 1e-10];

/// Index of the coefficient every clause is read at.
const REGISTERED: usize = 1;

/// C2's band, two-sided as registered.
const BAND: [f64; 2] = [1.7, 2.3];

/// The field C2 names as its resolution witness.
const WITNESS: &str = "box_exact";

/// `for_each_reference_field!` yields eight (`fields/mod.rs:195`).
const FIELDS: usize = 8;

/// C3's bar: at least six of the eight must discriminate.
const DISCRIMINATING_MIN: usize = 6;

/// Which part of the complex a tet edge lives in.
///
/// The three classes are exactly the three the hypothesis distinguishes: axis
/// edges are the grid edges `marching_cubes` already interpolates along, and the
/// other two are the set a chord predicate structurally cannot reach.
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Debug)]
enum Class {
    Axis,
    FaceDiagonal,
    BodyDiagonal,
}

impl Class {
    /// One axis step per set component of the offset.
    fn of(steps: u32) -> Self {
        match steps {
            1 => Self::Axis,
            2 => Self::FaceDiagonal,
            _ => Self::BodyDiagonal,
        }
    }
}

/// The local offset of a cube corner, as grid steps.
///
/// Copied from `crates/isomesh/src/cube.rs:149-155`, which is `pub(crate)`.
/// Bit 0 is `x`, bit 1 is `y`, bit 2 is `z`.
fn corner_offset(corner: u8) -> [usize; 3] {
    [
        usize::from(corner & 1),
        usize::from((corner >> 1) & 1),
        usize::from((corner >> 2) & 1),
    ]
}

/// The distinct positive lattice offsets carrying the tet edges of the shipped
/// Kuhn complex, derived from `TETS` (`table.rs:87`) and `TET_EDGES`
/// (`table.rs:121`).
///
/// Returns `(offset, class)` sorted by class then lexicographically, so the
/// sweep order is fixed by construction rather than by the order the tetrahedra
/// happen to be built in.
fn tet_edge_offsets() -> Vec<([usize; 3], Class)> {
    let mut out: Vec<([usize; 3], Class)> = Vec::new();
    for tet in &TETS {
        for [ea, eb] in TET_EDGES {
            let lo = corner_offset(tet[usize::from(ea)]);
            let hi = corner_offset(tet[usize::from(eb)]);
            // `table.rs:84-86`: TETS' corners are ordered by inclusion, so every
            // edge runs from fewer bits to more and its two offsets differ by a
            // 0/1 step on each axis and never by -1.
            assert!(
                hi[0] >= lo[0] && hi[1] >= lo[1] && hi[2] >= lo[2],
                "P-124: tet edge {lo:?} -> {hi:?} steps backwards, which table.rs:84-86 \
                 says cannot happen"
            );
            let delta = [hi[0] - lo[0], hi[1] - lo[1], hi[2] - lo[2]];
            let steps = delta.iter().sum::<usize>() as u32;
            let entry = (delta, Class::of(steps));
            if !out.contains(&entry) {
                out.push(entry);
            }
        }
    }
    out.sort_unstable_by(|a, b| a.1.cmp(&b.1).then_with(|| a.0.cmp(&b.0)));
    out
}

/// The 19 distinct cell-local tet edges, as `(class, count)` totals.
///
/// A structural check on the shipped table rather than a value the sweep uses:
/// 6 tetrahedra × 6 edges is 36 instances, and they must reduce to 12 cube
/// edges, 6 face diagonals and 1 body diagonal.
fn cell_local_census() -> [usize; 3] {
    let mut pairs: Vec<[u8; 2]> = Vec::new();
    for tet in &TETS {
        for [ea, eb] in TET_EDGES {
            let (a, b) = (tet[usize::from(ea)], tet[usize::from(eb)]);
            let pair = if a <= b { [a, b] } else { [b, a] };
            if !pairs.contains(&pair) {
                pairs.push(pair);
            }
        }
    }
    let mut census = [0usize; 3];
    for [a, b] in pairs {
        let steps = (a ^ b).count_ones() as usize;
        census[steps - 1] += 1;
    }
    census
}

/// Distinct tet edges in an `n³` sample grid: `3n²(n−1) + 3n(n−1)² + (n−1)³`.
fn closed_form_tet_edges(n: u64) -> u64 {
    let m = n - 1;
    3 * n * n * m + 3 * n * m * m + m * m * m
}

#[inline]
fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

/// One edge with its endpoint data already in hand.
///
/// The ambient arm fills `fa`/`fb`/`ga`/`gb` from the grid cache and the control
/// arm computes them on the spot, which is what makes the two arms literally the
/// same predicate rather than two implementations of it.
#[derive(Clone, Copy)]
struct Edge {
    a: [f64; 3],
    b: [f64; 3],
    fa: f64,
    fb: f64,
    ga: [f64; 3],
    gb: [f64; 3],
}

/// What one edge contributes, at all three coefficients from a single pass of
/// gradients — the `g(tᵢ)` do not depend on the coefficient, only the discard
/// rule does.
#[derive(Clone, Copy)]
struct EdgeOutcome {
    k: u32,
    /// Gradients evaluated *afresh* for this edge, i.e. the interior samples;
    /// the two endpoints come from the cache in the ambient arm.
    fresh_gradients: u32,
    /// One flag per entry of [`COEFFS`].
    non_monotone: [bool; 3],
    /// How many `g` the guard threw away, per coefficient.
    discarded: [u32; 3],
    /// How many of those were **not exactly zero**, per coefficient. Discarding
    /// an exact zero changes nothing — zero is neutral in the sign test either
    /// way — so this is the guard's only substantive action.
    nonzero_discarded: [u32; 3],
    /// Flagged by `g(0)` and `g(k−1)` alone at the registered coefficient, so no
    /// choice of `k` could have unflagged it. `✗36`'s 97.5%.
    flagged_by_endpoints: bool,
    /// Distance from zero the *minority* sign reached, at the registered
    /// coefficient. Tiny means the flag is noise; large means the field turns.
    reversal: f64,
    /// `max(|f(a)|, |f(b)|)` — the quantity the registered tolerance scales by.
    endpoint_scale: f64,
    zero_g: u32,
    nonfinite_g: u32,
    degenerate: bool,
}

/// The registered predicate, evaluated once for one edge.
///
/// `k = max(2, ceil(‖e‖/w) + 1)` is Finken §4 verbatim; the tolerance
/// `coef · max(|f(a)|, |f(b)|)` is the registration's, and is the one thing
/// `✗36` asked to be changed.
fn examine<F>(field: &F, e: Edge, w: f64, g: &mut Vec<f64>) -> EdgeOutcome
where
    F: Sdf<Scalar = f64>,
{
    let d = [e.b[0] - e.a[0], e.b[1] - e.a[1], e.b[2] - e.a[2]];
    let len = dot(d, d).sqrt();
    // A saturating cast, so a non-finite length lands on `k = 2` rather than on
    // an arbitrary index; `nonfinite_g` is what would then report it.
    let k = 2.max((len / w).ceil() as u32 + 1);

    g.clear();
    g.push(dot(e.ga, d));
    let mut fresh_gradients = 0u32;
    for i in 1..k - 1 {
        let t = f64::from(i) / f64::from(k - 1);
        let p = [e.a[0] + t * d[0], e.a[1] + t * d[1], e.a[2] + t * d[2]];
        g.push(dot(field.gradient(p), d));
        fresh_gradients += 1;
    }
    g.push(dot(e.gb, d));

    let endpoint_scale = e.fa.abs().max(e.fb.abs());

    let mut out = EdgeOutcome {
        k,
        fresh_gradients,
        non_monotone: [false; 3],
        discarded: [0; 3],
        nonzero_discarded: [0; 3],
        flagged_by_endpoints: false,
        reversal: 0.0,
        endpoint_scale,
        zero_g: 0,
        nonfinite_g: 0,
        // A total comparison, so a NaN endpoint is a degenerate edge rather than
        // an unordered one that quietly counts as ordinary.
        degenerate: !len.is_finite() || len.total_cmp(&0.0) != Ordering::Greater,
    };
    for &v in g.iter() {
        if v.is_finite() {
            if v.abs() > 0.0 {
                continue;
            }
            out.zero_g += 1;
        } else {
            out.nonfinite_g += 1;
        }
    }

    for (c, &coef) in COEFFS.iter().enumerate() {
        let tol = coef * endpoint_scale;
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
            let (g0, g1) = (g[0], g[g.len() - 1]);
            out.flagged_by_endpoints = g0.abs() >= tol
                && g1.abs() >= tol
                && ((g0 > 0.0 && g1 < 0.0) || (g0 < 0.0 && g1 > 0.0));
            if flagged {
                out.reversal = pos.min(neg);
            }
        }
    }
    out
}

/// The ambient arm's measurement for one `(field, resolution)`.
struct Ambient {
    cells: u64,
    tet_edges: u64,
    per_class: [u64; 3],
    non_monotone_edges: [u64; 3],
    non_monotone_per_class: [u64; 3],
    non_monotone_cells: u64,
    axis_non_monotone_cells: u64,
    diagonal_only_cells: u64,
    interior_vertices: u64,
    pl_local_extrema: u64,
    certified_pl_local_extrema: u64,
    k_min: u32,
    k_max: u32,
    sample_evals: u64,
    gradient_evals: u64,
    discarded: [u64; 3],
    nonzero_discarded: [u64; 3],
    zero_g: u64,
    nonfinite_g: u64,
    degenerate: u64,
    tol_max: f64,
    tol_sum: f64,
    scale_min: f64,
    worst_reversal: f64,
    flagged_by_endpoints: u64,
    predicate_ns: f64,
}

impl Ambient {
    fn certified_cells(&self) -> u64 {
        self.cells - self.non_monotone_cells
    }
}

/// Sweep every distinct tet edge of the ambient complex once.
fn measure_ambient<F>(field: &F, samples: u32, offsets: &[([usize; 3], Class)]) -> Ambient
where
    F: ReferenceField + Sdf<Scalar = f64>,
{
    let (_shape, origin, w) = common::grid(field, samples);
    let n = samples as usize;
    let m = n - 1;
    let vid = |v: [usize; 3]| v[0] + n * (v[1] + n * v[2]);
    let cid = |c: [usize; 3]| c[0] + m * (c[1] + m * c[2]);
    let pos = |v: [usize; 3]| {
        [
            origin[0] + w * v[0] as f64,
            origin[1] + w * v[1] as f64,
            origin[2] + w * v[2] as f64,
        ]
    };

    // One evaluation per grid sample, reused by every incident tet edge, so a
    // shared edge's flag cannot depend on which cell asked for it.
    let mut values = vec![0.0_f64; n * n * n];
    let mut grads = vec![[0.0_f64; 3]; n * n * n];
    for k in 0..n {
        for j in 0..n {
            for i in 0..n {
                let p = pos([i, j, k]);
                let at = vid([i, j, k]);
                values[at] = field.sample(p);
                grads[at] = field.gradient(p);
            }
        }
    }

    let mut out = Ambient {
        cells: (m * m * m) as u64,
        tet_edges: 0,
        per_class: [0; 3],
        non_monotone_edges: [0; 3],
        non_monotone_per_class: [0; 3],
        non_monotone_cells: 0,
        axis_non_monotone_cells: 0,
        diagonal_only_cells: 0,
        interior_vertices: 0,
        pl_local_extrema: 0,
        certified_pl_local_extrema: 0,
        k_min: u32::MAX,
        k_max: 0,
        sample_evals: (n * n * n) as u64,
        gradient_evals: (n * n * n) as u64,
        discarded: [0; 3],
        nonzero_discarded: [0; 3],
        zero_g: 0,
        nonfinite_g: 0,
        degenerate: 0,
        tol_max: 0.0,
        tol_sum: 0.0,
        scale_min: f64::INFINITY,
        worst_reversal: 0.0,
        flagged_by_endpoints: 0,
        predicate_ns: 0.0,
    };

    let mut cell_dirty = vec![false; m * m * m];
    let mut cell_axis_dirty = vec![false; m * m * m];
    let mut vertex_dirty = vec![false; n * n * n];
    let mut g: Vec<f64> = Vec::with_capacity(8);

    let started = Instant::now();
    for &(delta, class) in offsets {
        let span = [n - delta[0], n - delta[1], n - delta[2]];
        for vz in 0..span[2] {
            for vy in 0..span[1] {
                for vx in 0..span[0] {
                    let v = [vx, vy, vz];
                    let u = [vx + delta[0], vy + delta[1], vz + delta[2]];
                    let (ia, ib) = (vid(v), vid(u));
                    let outcome = examine(
                        field,
                        Edge {
                            a: pos(v),
                            b: pos(u),
                            fa: values[ia],
                            fb: values[ib],
                            ga: grads[ia],
                            gb: grads[ib],
                        },
                        w,
                        &mut g,
                    );

                    out.tet_edges += 1;
                    out.per_class[class as usize] += 1;
                    out.k_min = out.k_min.min(outcome.k);
                    out.k_max = out.k_max.max(outcome.k);
                    out.gradient_evals += u64::from(outcome.fresh_gradients);
                    out.zero_g += u64::from(outcome.zero_g);
                    out.nonfinite_g += u64::from(outcome.nonfinite_g);
                    out.degenerate += u64::from(outcome.degenerate);
                    if outcome.endpoint_scale < out.scale_min {
                        out.scale_min = outcome.endpoint_scale;
                    }
                    let tol = COEFFS[REGISTERED] * outcome.endpoint_scale;
                    out.tol_sum += tol;
                    if tol > out.tol_max {
                        out.tol_max = tol;
                    }
                    for c in 0..COEFFS.len() {
                        out.discarded[c] += u64::from(outcome.discarded[c]);
                        out.nonzero_discarded[c] += u64::from(outcome.nonzero_discarded[c]);
                        if outcome.non_monotone[c] {
                            out.non_monotone_edges[c] += 1;
                        }
                    }
                    if !outcome.non_monotone[REGISTERED] {
                        continue;
                    }

                    out.non_monotone_per_class[class as usize] += 1;
                    if outcome.flagged_by_endpoints {
                        out.flagged_by_endpoints += 1;
                    }
                    if outcome.reversal > out.worst_reversal {
                        out.worst_reversal = outcome.reversal;
                    }
                    vertex_dirty[ia] = true;
                    vertex_dirty[ib] = true;

                    // Incident cells. A set axis of `delta` pins that cell
                    // coordinate; a clear one leaves the two cells sharing the
                    // lattice line through the edge, clipped at the boundary.
                    let mut lo = [0usize; 3];
                    let mut hi = [0usize; 3];
                    for (a, &step) in delta.iter().enumerate() {
                        if step == 1 {
                            lo[a] = v[a];
                            hi[a] = v[a];
                        } else {
                            lo[a] = v[a].saturating_sub(1);
                            hi[a] = v[a].min(m - 1);
                        }
                    }
                    for cz in lo[2]..=hi[2] {
                        for cy in lo[1]..=hi[1] {
                            for cx in lo[0]..=hi[0] {
                                let at = cid([cx, cy, cz]);
                                cell_dirty[at] = true;
                                if class == Class::Axis {
                                    cell_axis_dirty[at] = true;
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    for (dirty, axis) in cell_dirty.iter().zip(cell_axis_dirty.iter()) {
        if !*dirty {
            continue;
        }
        out.non_monotone_cells += 1;
        if *axis {
            out.axis_non_monotone_cells += 1;
        } else {
            out.diagonal_only_cells += 1;
        }
    }

    // Theorem 1 Part 1: a strict PL local extremum over the 14-neighbour link is
    // a critical point of `f̂`, and in a star of monotone edges it coincides with
    // a critical point of `f`. Interior means all fourteen neighbours exist,
    // which for `±` the seven offsets is exactly `1 ≤ v ≤ n − 2` on every axis.
    for vz in 1..m {
        for vy in 1..m {
            for vx in 1..m {
                let v = [vx, vy, vz];
                let at = vid(v);
                out.interior_vertices += 1;
                let here = values[at];
                // `all_above`: every link neighbour is strictly above `here`, so
                // `here` is a strict local minimum. `all_below` is the maximum
                // case. `f̂` is affine on each tet, so a strict link extremum is
                // a strict extremum on every incident tet and therefore a
                // critical point of `f̂` — in any dimension.
                let mut all_above = true;
                let mut all_below = true;
                'link: for &(delta, _) in offsets {
                    for sign in [1isize, -1] {
                        let nb = [
                            (vx as isize + sign * delta[0] as isize) as usize,
                            (vy as isize + sign * delta[1] as isize) as usize,
                            (vz as isize + sign * delta[2] as isize) as usize,
                        ];
                        match values[vid(nb)].total_cmp(&here) {
                            Ordering::Greater => all_below = false,
                            Ordering::Less => all_above = false,
                            Ordering::Equal => {
                                all_above = false;
                                all_below = false;
                            }
                        }
                        if !(all_above || all_below) {
                            break 'link;
                        }
                    }
                }
                if !(all_above || all_below) {
                    continue;
                }
                out.pl_local_extrema += 1;
                if !vertex_dirty[at] {
                    out.certified_pl_local_extrema += 1;
                }
            }
        }
    }

    out.predicate_ns = started.elapsed().as_secs_f64() * 1e9;
    if out.tet_edges == 0 {
        out.k_min = 0;
    }
    if !out.scale_min.is_finite() {
        out.scale_min = 0.0;
    }
    out
}

/// One undirected mesh edge, endpoints in canonical order.
///
/// Deduplication is by the **exact bit pattern** of the two endpoint positions,
/// transcribed from `experiment_p55.rs:231-289`: `MeshBuffer` never welds, so
/// marching cubes emits three fresh vertices per triangle, and welding by
/// proximity would need an epsilon this experiment has no business inventing.
#[derive(Clone, Copy)]
struct Segment {
    a: [f64; 3],
    b: [f64; 3],
}

/// Lexicographic total order on a point. `total_cmp` rather than `partial_cmp`
/// so a NaN coordinate is ordered rather than silently ungrouped.
fn lex(a: &[f64; 3], b: &[f64; 3]) -> Ordering {
    a[0].total_cmp(&b[0])
        .then_with(|| a[1].total_cmp(&b[1]))
        .then_with(|| a[2].total_cmp(&b[2]))
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

/// The mesh-edge control arm: `✗36`'s reading, in this run, on this population.
struct Control {
    triangles: u64,
    edges: u64,
    non_monotone: u64,
    degenerate: u64,
    tol_max: f64,
    nonzero_discarded: u64,
    extract_ns: f64,
}

impl Control {
    /// Non-monotone per 1k distinct segments, `✗36`'s own unit.
    fn per_1k(&self) -> f64 {
        if self.edges == 0 {
            0.0
        } else {
            1000.0 * self.non_monotone as f64 / self.edges as f64
        }
    }

    /// Saturated rather than discriminating: a strict majority of the chords
    /// reverse, which is what `✗36` proved for a strictly convex surface.
    fn saturated(&self) -> bool {
        2 * self.non_monotone > self.edges
    }
}

/// Extract at `samples³` and run the **same** predicate over every distinct mesh
/// edge.
fn measure_control<F>(field: &F, samples: u32) -> Control
where
    F: ReferenceField + Sdf<Scalar = f64>,
{
    let (shape, origin, w) = common::grid(field, samples);
    let mut mesh = MeshBuffer::<f64>::new();
    let mut mc = MarchingCubes::<f64>::new();
    let started = Instant::now();
    mc.extract(field, &shape, origin, w, &mut mesh)
        .expect("marching cubes on a reference field");
    let extract_ns = started.elapsed().as_secs_f64() * 1e9;

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

    let mut out = Control {
        triangles: mesh.triangle_count() as u64,
        edges: 0,
        non_monotone: 0,
        degenerate: 0,
        tol_max: 0.0,
        nonzero_discarded: 0,
        extract_ns,
    };
    let mut g: Vec<f64> = Vec::with_capacity(8);
    let mut previous: Option<Segment> = None;
    for seg in all {
        if previous.is_some_and(|last| seg.order(&last) == Ordering::Equal) {
            continue;
        }
        previous = Some(seg);
        out.edges += 1;
        let outcome = examine(
            field,
            Edge {
                a: seg.a,
                b: seg.b,
                fa: field.sample(seg.a),
                fb: field.sample(seg.b),
                ga: field.gradient(seg.a),
                gb: field.gradient(seg.b),
            },
            w,
            &mut g,
        );
        if outcome.non_monotone[REGISTERED] {
            out.non_monotone += 1;
        }
        out.degenerate += u64::from(outcome.degenerate);
        out.nonzero_discarded += u64::from(outcome.nonzero_discarded[REGISTERED]);
        let tol = COEFFS[REGISTERED] * outcome.endpoint_scale;
        if tol > out.tol_max {
            out.tol_max = tol;
        }
    }
    out
}

/// One CSV row, before the cross-row verdicts exist.
struct Row {
    field: &'static str,
    samples: u32,
    cell: f64,
    ambient: Ambient,
    control: Control,
}

impl Row {
    /// C1 and C3's per-row reading: the predicate said something other than
    /// "all" or "none".
    fn discriminates(&self) -> bool {
        self.ambient.non_monotone_cells > 0 && self.ambient.non_monotone_cells < self.ambient.cells
    }

    /// Non-monotone cells per cell — a dimensionless density over a denominator
    /// that is exact by construction.
    fn rate(&self) -> f64 {
        self.ambient.non_monotone_cells as f64 / self.ambient.cells as f64
    }
}

/// `rate(coarse) / rate(fine)`, as **one** division of integers-as-`f64`.
///
/// Cross-multiplied rather than taken as a quotient of quotients, so it cannot
/// pick up a second rounding, and it returns `0.0` rather than a NaN when the
/// finer population is empty — a value outside [`BAND`], which is the honest
/// verdict for a witness that vanished.
fn rate_ratio(coarse: &Row, fine: &Row) -> f64 {
    let num = coarse.ambient.non_monotone_cells as f64 * fine.ambient.cells as f64;
    let den = fine.ambient.non_monotone_cells as f64 * coarse.ambient.cells as f64;
    if den > 0.0 { num / den } else { 0.0 }
}

/// `count(fine) / count(coarse)` — the registration's literal reading, recorded
/// so the file settles whether the population falls at all.
fn count_ratio(coarse: &Row, fine: &Row) -> f64 {
    let den = coarse.ambient.non_monotone_cells as f64;
    if den > 0.0 {
        fine.ambient.non_monotone_cells as f64 / den
    } else {
        0.0
    }
}

fn in_band(x: f64) -> bool {
    x.total_cmp(&BAND[0]) != Ordering::Less && x.total_cmp(&BAND[1]) != Ordering::Greater
}

fn main() {
    if !std::env::args().any(|arg| arg == "--bench") {
        return;
    }
    let prereg = isomesh::experiment!("P-124");
    common::experiment::run(prereg, |run| {
        // ── the complex, checked against the shipped table before any field ──
        let offsets = tet_edge_offsets();
        let census = cell_local_census();
        assert_eq!(
            census,
            [12, 6, 1],
            "P-124: the shipped six-tet split must carry 12 cube edges, 6 face \
             diagonals and 1 body diagonal per cell (table.rs:36-53)"
        );
        assert_eq!(
            offsets.len(),
            7,
            "P-124: the Kuhn complex's tet edges must reduce to seven positive \
             lattice offsets, got {offsets:?}"
        );
        let by_class = |c: Class| offsets.iter().filter(|(_, k)| *k == c).count();
        assert_eq!(
            [
                by_class(Class::Axis),
                by_class(Class::FaceDiagonal),
                by_class(Class::BodyDiagonal),
            ],
            [3, 3, 1],
            "P-124: three axis offsets, three face diagonals, one body diagonal"
        );
        println!(
            "complex: Kuhn's six tets, {} distinct cell-local edges (12 axis / 6 face \
             diagonal / 1 body diagonal), {} positive lattice offsets\n",
            census.iter().sum::<usize>(),
            offsets.len()
        );

        // ── the sweep: eight fields × three resolutions, one binary ──────────
        let mut rows: Vec<Row> = Vec::new();
        isomesh::for_each_reference_field!(f64, |name, field| {
            for samples in RESOLUTIONS {
                let ambient = measure_ambient(&field, samples, &offsets);
                let control = measure_control(&field, samples);
                let (_shape, _origin, cell) = common::grid(&field, samples);

                // ── M-44, the registration's named control ──────────────────
                //
                // The same `examine`, in the same run, must demonstrably flag a
                // non-monotone edge — otherwise an ambient zero is a silence and
                // not a measurement. The registration's parenthetical asks for
                // `0` here; `✗36` measured 917.6 / 893.7 / 42.3 per 1k and
                // proved the chord case in closed form, so `0` is unreachable
                // and asserting it would abort before the first row.
                assert!(
                    control.non_monotone > 0,
                    "P-124: {name} at {samples}³ — the mesh-edge control read 0 of \
                     {} non-monotone, so this run cannot show the predicate able to \
                     fire and no ambient count from it is a measurement (M-44)",
                    control.edges
                );
                if name == "sphere" {
                    assert!(
                        control.saturated(),
                        "P-124: sphere at {samples}³ — the mesh-edge control read \
                         {} of {} non-monotone, not a strict majority. ✗36 proves \
                         every chord of a strictly convex surface reverses, so a \
                         non-saturated control means this is not ✗36's predicate \
                         and the ambient reading's discrimination cannot be \
                         attributed to the complex",
                        control.non_monotone,
                        control.edges
                    );
                }
                assert_eq!(
                    ambient.tet_edges,
                    closed_form_tet_edges(u64::from(samples)),
                    "P-124: {name} at {samples}³ — swept {} tet edges against the \
                     closed form; the enumeration is not the complex",
                    ambient.tet_edges
                );
                assert!(
                    ambient.per_class[Class::FaceDiagonal as usize]
                        + ambient.per_class[Class::BodyDiagonal as usize]
                        > 0,
                    "P-124: {name} at {samples}³ — no diagonal was swept, so a zero \
                     in diagonal_only_failures would mean nothing"
                );

                println!(
                    "{name:>14} {samples:>4}³  cells {:>9}  tet-edges {:>10} \
                     (ax {:>9} / fd {:>9} / bd {:>9})  nm-edges {:>9}  nm-cells {:>9} \
                     ({:>8.6})  diag-only {:>8}  crit≥ {:>6}  k {}..{}  \
                     tol≤{:>9.2e}  control {:>8}/{:>8} ({:>8.3}/1k)",
                    ambient.cells,
                    ambient.tet_edges,
                    ambient.per_class[Class::Axis as usize],
                    ambient.per_class[Class::FaceDiagonal as usize],
                    ambient.per_class[Class::BodyDiagonal as usize],
                    ambient.non_monotone_edges[REGISTERED],
                    ambient.non_monotone_cells,
                    ambient.non_monotone_cells as f64 / ambient.cells as f64,
                    ambient.diagonal_only_cells,
                    ambient.certified_pl_local_extrema,
                    ambient.k_min,
                    ambient.k_max,
                    ambient.tol_max,
                    control.non_monotone,
                    control.edges,
                    control.per_1k(),
                );
                rows.push(Row {
                    field: name,
                    samples,
                    cell,
                    ambient,
                    control,
                });
            }
        });

        // The macro inlines its body once per field and a `return` inside it
        // would exit `main` silently (`M-199`, `fields/mod.rs:198-210`), so the
        // population every clause is denominated in is checked before any
        // verdict is computed.
        assert_eq!(
            rows.len(),
            FIELDS * RESOLUTIONS.len(),
            "P-124: the sweep must cover all eight reference fields at all three \
             resolutions; C1 is `over every cell of every field` and C3's \
             denominator is the eight"
        );

        // ── C1: an integer over every cell, both falsifiers as registered ────
        let any_non_zero = rows.iter().any(|r| r.ambient.non_monotone_cells > 0);
        let any_saturated = rows
            .iter()
            .any(|r| r.ambient.non_monotone_cells == r.ambient.cells);
        let c1_holds = any_non_zero && !any_saturated;

        // ── C2: box_exact's density, ratio per doubling, band [1.7, 2.3] ─────
        let witness: Vec<&Row> = RESOLUTIONS
            .iter()
            .filter_map(|&n| rows.iter().find(|r| r.field == WITNESS && r.samples == n))
            .collect();
        assert_eq!(
            witness.len(),
            RESOLUTIONS.len(),
            "P-124: C2's resolution witness `{WITNESS}` is missing a rung of the ladder"
        );
        let rate_ratios = [
            rate_ratio(witness[0], witness[1]),
            rate_ratio(witness[1], witness[2]),
        ];
        let count_ratios = [
            count_ratio(witness[0], witness[1]),
            count_ratio(witness[1], witness[2]),
        ];
        let count_falls = witness[1].ambient.non_monotone_cells
            < witness[0].ambient.non_monotone_cells
            && witness[2].ambient.non_monotone_cells < witness[1].ambient.non_monotone_cells;
        // The gated scalar: whichever of the two doublings sits furthest from a
        // clean halving, decided by `total_cmp`.
        let worst_ratio = if (rate_ratios[0] - 2.0)
            .abs()
            .total_cmp(&(rate_ratios[1] - 2.0).abs())
            == Ordering::Greater
        {
            rate_ratios[0]
        } else {
            rate_ratios[1]
        };
        let c2_holds = in_band(rate_ratios[0]) && in_band(rate_ratios[1]);

        // ── C3: at least six of eight discriminate, at every resolution ──────
        let mut discriminating: Vec<(u32, usize)> = Vec::new();
        for &n in &RESOLUTIONS {
            let here = rows
                .iter()
                .filter(|r| r.samples == n && r.discriminates())
                .count();
            discriminating.push((n, here));
        }
        let c3_holds = discriminating.iter().all(|&(_, k)| k >= DISCRIMINATING_MIN);

        println!(
            "\nC1  non-monotone cells: any non-zero {any_non_zero}, any saturated \
             {any_saturated} -> {c1_holds}\n\
             C2  {WITNESS} density ratios {:.6} / {:.6} (band {:.1}..{:.1}) -> {c2_holds}; \
             count ratios {:.6} / {:.6}, count falls {count_falls}\n\
             C3  discriminating fields {:?} of {FIELDS} (bar {DISCRIMINATING_MIN}) -> \
             {c3_holds}\n",
            rate_ratios[0],
            rate_ratios[1],
            BAND[0],
            BAND[1],
            count_ratios[0],
            count_ratios[1],
            discriminating,
        );

        for row in &rows {
            let a = &row.ambient;
            let previous = RESOLUTIONS
                .iter()
                .position(|&n| n == row.samples)
                .filter(|&i| i > 0)
                .and_then(|i| {
                    rows.iter()
                        .find(|r| r.field == row.field && r.samples == RESOLUTIONS[i - 1])
                });
            let own_ratio = previous.map_or(0.0, |p| rate_ratio(p, row));
            let here = discriminating
                .iter()
                .find(|(n, _)| *n == row.samples)
                .map_or(0, |&(_, k)| k);
            let edges = a.tet_edges.max(1) as f64;

            run.record(&[
                // ── the fourteen registered metrics ──────────────────────────
                ("field", row.field.to_string()),
                ("resolution", row.samples.to_string()),
                ("cells", a.cells.to_string()),
                // Six tetrahedra per cell, `table.rs:72`.
                ("tets", (6 * a.cells).to_string()),
                ("tet_edges", a.tet_edges.to_string()),
                (
                    "non_monotone_edges",
                    a.non_monotone_edges[REGISTERED].to_string(),
                ),
                ("non_monotone_cells", a.non_monotone_cells.to_string()),
                // Cells rejected by a diagonal with every axis edge monotone —
                // the part of the signal a chord predicate cannot reach.
                ("diagonal_only_failures", a.diagonal_only_cells.to_string()),
                // Theorem 1 Part 1, over certified stars. Gates nothing.
                (
                    "critical_points_lower_bound",
                    a.certified_pl_local_extrema.to_string(),
                ),
                ("box_exact_population_ratio", format!("{worst_ratio:.6}")),
                (
                    "mesh_edge_control_non_monotone",
                    row.control.non_monotone.to_string(),
                ),
                ("c1_holds", c1_holds.to_string()),
                ("c2_holds", c2_holds.to_string()),
                ("c3_holds", c3_holds.to_string()),
                // ── the complex, so every denominator is checkable ───────────
                (
                    "tet_edges_axis",
                    a.per_class[Class::Axis as usize].to_string(),
                ),
                (
                    "tet_edges_face_diagonal",
                    a.per_class[Class::FaceDiagonal as usize].to_string(),
                ),
                (
                    "tet_edges_body_diagonal",
                    a.per_class[Class::BodyDiagonal as usize].to_string(),
                ),
                (
                    "tet_edges_closed_form",
                    closed_form_tet_edges(u64::from(row.samples)).to_string(),
                ),
                (
                    "tet_edge_count_matches_closed_form",
                    (a.tet_edges == closed_form_tet_edges(u64::from(row.samples))).to_string(),
                ),
                // ── where the non-monotonicity lives ────────────────────────
                (
                    "non_monotone_axis_edges",
                    a.non_monotone_per_class[Class::Axis as usize].to_string(),
                ),
                (
                    "non_monotone_face_diagonal_edges",
                    a.non_monotone_per_class[Class::FaceDiagonal as usize].to_string(),
                ),
                (
                    "non_monotone_body_diagonal_edges",
                    a.non_monotone_per_class[Class::BodyDiagonal as usize].to_string(),
                ),
                (
                    "non_monotone_edge_rate",
                    format!("{:.9}", a.non_monotone_edges[REGISTERED] as f64 / edges),
                ),
                (
                    "axis_only_reachable_cells",
                    a.axis_non_monotone_cells.to_string(),
                ),
                ("certified_cells", a.certified_cells().to_string()),
                ("non_monotone_cell_rate", format!("{:.9}", row.rate())),
                ("c1_row_discriminates", row.discriminates().to_string()),
                ("degenerate_edges", a.degenerate.to_string()),
                // ── C2's two readings, side by side ────────────────────────
                (
                    "box_exact_rate_ratio_33_to_65",
                    format!("{:.6}", rate_ratios[0]),
                ),
                (
                    "box_exact_rate_ratio_65_to_129",
                    format!("{:.6}", rate_ratios[1]),
                ),
                (
                    "box_exact_count_ratio_33_to_65",
                    format!("{:.6}", count_ratios[0]),
                ),
                (
                    "box_exact_count_ratio_65_to_129",
                    format!("{:.6}", count_ratios[1]),
                ),
                ("box_exact_count_falls", count_falls.to_string()),
                ("rate_ratio_from_previous", format!("{own_ratio:.6}")),
                // ── C3's tally, on the row that contributed to it ──────────
                ("discriminating_fields_here", here.to_string()),
                ("discriminating_fields_bar", DISCRIMINATING_MIN.to_string()),
                // ── Theorem 1 Part 1's population, so a zero is readable ───
                ("interior_vertices", a.interior_vertices.to_string()),
                ("pl_local_extrema", a.pl_local_extrema.to_string()),
                // ── the tolerance, and whether it did anything ─────────────
                //
                // Semicolon rather than comma on purpose: `Run::record` panics
                // on a value containing a separator because the writer does not
                // quote (`benches/common/experiment.rs:52-65`, `P-64`).
                ("tolerance_rule", "coef*max(|f(a)|;|f(b)|)".to_string()),
                ("tolerance_coef", format!("{:e}", COEFFS[REGISTERED])),
                ("max_abs_tolerance", format!("{:.6e}", a.tol_max)),
                ("mean_abs_tolerance", format!("{:.6e}", a.tol_sum / edges)),
                ("min_endpoint_scale", format!("{:.6e}", a.scale_min)),
                (
                    "non_monotone_edges_at_1e14",
                    a.non_monotone_edges[0].to_string(),
                ),
                (
                    "non_monotone_edges_at_1e10",
                    a.non_monotone_edges[2].to_string(),
                ),
                (
                    "counts_equal_across_tolerances",
                    (a.non_monotone_edges[0] == a.non_monotone_edges[REGISTERED]
                        && a.non_monotone_edges[REGISTERED] == a.non_monotone_edges[2])
                        .to_string(),
                ),
                (
                    "nonzero_g_discarded_at_1e12",
                    a.nonzero_discarded[REGISTERED].to_string(),
                ),
                (
                    "guard_inert",
                    (a.nonzero_discarded[0] == 0
                        && a.nonzero_discarded[REGISTERED] == 0
                        && a.nonzero_discarded[2] == 0)
                        .to_string(),
                ),
                ("zero_g_samples", a.zero_g.to_string()),
                ("nonfinite_g_samples", a.nonfinite_g.to_string()),
                // ── the predicate's own parameters ────────────────────────
                ("cell_size_w", format!("{:.9}", row.cell)),
                ("k_samples_min", a.k_min.to_string()),
                ("k_samples_max", a.k_max.to_string()),
                ("sample_evals", a.sample_evals.to_string()),
                ("gradient_evals", a.gradient_evals.to_string()),
                ("worst_reversal", format!("{:.6e}", a.worst_reversal)),
                (
                    "flagged_by_endpoints_at_1e12",
                    a.flagged_by_endpoints.to_string(),
                ),
                (
                    "all_flags_from_endpoints",
                    (a.flagged_by_endpoints == a.non_monotone_edges[REGISTERED]).to_string(),
                ),
                // ── the mesh-edge control arm, ✗36 in the same run ────────
                ("mesh_edge_control_edges", row.control.edges.to_string()),
                (
                    "mesh_edge_control_triangles",
                    row.control.triangles.to_string(),
                ),
                (
                    "mesh_edge_control_per_1k",
                    format!("{:.6}", row.control.per_1k()),
                ),
                (
                    "mesh_edge_control_saturated",
                    row.control.saturated().to_string(),
                ),
                (
                    "mesh_edge_control_degenerate",
                    row.control.degenerate.to_string(),
                ),
                (
                    "mesh_edge_control_max_abs_tolerance",
                    format!("{:.6e}", row.control.tol_max),
                ),
                (
                    "mesh_edge_control_nonzero_g_discarded",
                    row.control.nonzero_discarded.to_string(),
                ),
                // ── time, recorded beside the verdict, gating nothing ─────
                ("predicate_ns", format!("{:.0}", a.predicate_ns)),
                (
                    "predicate_ns_per_edge",
                    format!("{:.2}", a.predicate_ns / edges),
                ),
                (
                    "control_extract_ns",
                    format!("{:.0}", row.control.extract_ns),
                ),
            ]);
        }
    });
}
