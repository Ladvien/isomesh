//! **P-139 — is the Kuhn triangulation regular, so Viro patchworking transfers?**
//!
//! Ticket: R-139. Pre-registered before this harness existed.
//!
//! ```bash
//! cargo bench --bench experiment_p139
//! ```
//!
//! Writes `docs/experiments/p-139.csv`.
//!
//! # What was missing
//!
//! Viro's combinatorial patchworking — the **T-construction** — takes a sign at
//! every lattice point of a Newton polytope, glues a piecewise-linear
//! hypersurface over a triangulation of that polytope, and concludes that the
//! result is *isotopic to the zero set of a real polynomial with that Newton
//! polytope*. Signs at the eight corners of a cube, glued over Kuhn's six
//! tetrahedra, is **marching tetrahedra**. So the theorem is either a free gift
//! to `isomesh::marching_tetrahedra` or it is inapplicable, and **exactly one
//! fact decides which**: the theorem's hypothesis is that the triangulation be
//! **regular** — induced by a convex lifting of the lattice points.
//!
//! Nothing in this repository has ever asked. `marching_tetrahedra/table.rs:15-53`
//! proves the six tetrahedra *tile without alternation* — the property `✗11` is
//! about — and `P-100 / ✗78 / M-412` measured `open_edges = 0` on 80 of 80 rows
//! across a chunk seam. Both are about **matching**, not about **convexity**.
//! `P-124` swept the same complex for Finken et al.'s monotone-edge condition
//! (`experiment_p124.rs:76-117`) and derived its 19 cell-local edges, and still
//! never lifted it. Regularity is a different question with a different
//! instrument, and this row is only that instrument plus the two cross-checks
//! the registration asks for.
//!
//! `2026-08-23`'s reading — *"marching tetrahedra contours the Lovász
//! extension"* — is the third side of the same statement, and C3 closes it.
//!
//! # C1 — regularity, decided exactly rather than recalled
//!
//! The complex is taken from the crate, not re-derived:
//! `isomesh::marching_tetrahedra::table::TETS` (`table.rs:87`) and `TET_EDGES`
//! (`table.rs:121`) are `pub`, so the six tetrahedra swept here are bit-for-bit
//! the six `MarchingTetrahedra` marches — the six *monotone paths* from corner
//! `0` to corner `7`, one per ordering of the axes (`table.rs:17-25`).
//!
//! A triangulation of a point configuration `A` is **regular** iff some height
//! function `w: A -> R` has a lower convex hull projecting to exactly that
//! triangulation. Equivalently, and this is the form that is a *linear
//! program*: for every cell `T` the affine function `g_T` interpolating `w` on
//! `T`'s vertices must satisfy `w(v) > g_T(v)` for **every** `v` in `A` outside
//! `T` — `g_T` supports the hull touching precisely at `T`. For 8 corners and 6
//! tetrahedra that is `6 x 4 = 24` strict inequalities in 8 unknowns.
//!
//! Every coefficient is an **integer**, and the reason is worth stating because
//! it is what makes the whole clause exact in `i64`. The barycentric coordinates
//! come from Cramer's rule on the homogeneous matrix `[[p_i, 1]]`, whose
//! determinant is the simplex's normalised volume — and every Kuhn tetrahedron
//! is **unimodular**, `|det| = 1`, asserted on all six before anything is
//! solved. So `det * lambda_i` is an integer and `lambda_i` is one too. For
//! `T = {0, 1, 3, 7}` (the `x >= y >= z` ordering) the interpolant at `(x, y, z)`
//! is `(1-x)w_0 + (x-y)w_1 + (y-z)w_3 + z*w_7`, coefficients in `{-1, 0, 1}`.
//!
//! The 24 rows reduce, after dividing by the gcd, to **9 distinct** ones, and
//! the identity of those nine is C3's punchline rather than a curiosity — see
//! below. The system is **homogeneous** (every row's coefficients sum to zero,
//! because a row is an affine dependency), so it has no constant column and
//! Fourier--Motzkin elimination stays homogeneous: an all-zero row means
//! `0 > 0` and the system is infeasible.
//!
//! Elimination order `w_7, w_6, ..., w_0`. Measured stage sizes are
//! `9, 3, 2, 1, 1, 0, 0, 0, 0` — no blow-up, no coefficient growth, every entry
//! stays in `{-1, 0, 1}`, so `i128` here is belt-and-braces rather than
//! necessary. Back-substitution then *constructs* a witness: pick each variable
//! strictly between its lower and upper bounds, `lo + 1` / `hi - 1` when one
//! side is open and `0` when both are.
//!
//! **What the witness turns out to be is the interesting part.** The
//! back-substitution is told nothing about set functions and returns
//! `w = (0, 0, 0, -1, 0, -1, -1, -3)`, which read as a function of the corner's
//! popcount is
//!
//! ```text
//! w(S) = -|S|(|S| - 1)/2      -- minus the number of unordered pairs in S
//! ```
//!
//! the negated edge count of the complete graph induced on `S`. Its geometry is
//! the picture: the main diagonal `0-7` is pushed lowest, and the six
//! tetrahedra — all of which contain that diagonal (`table.rs:19-20`) — fan
//! around it. `secondary_polytope_vertex` is that triangulation's GKZ vector,
//! computed as the normalised volume of the tetrahedra at each corner:
//! `6|2|2|2|2|2|2|6`, summing to `4 x 6 = 24`. By Gel'fand--Kapranov--Zelevinsky
//! the regular triangulations are exactly the **vertices** of the secondary
//! polytope, so `is_regular = true` is what makes that vector a vertex.
//!
//! The witness is then verified twice, by two routes that share no code with the
//! solver: against all 24 inequalities directly, and against a **brute-force
//! lower-hull enumeration** — all 70 four-subsets of the eight corners, each
//! tested for "every other lifted corner strictly above my affine hull". That
//! enumeration finds simplicial lower facets only, which is exactly the right
//! instrument: a subdivision with a non-simplicial cell reports *fewer* than six.
//!
//! # The two liftings this row was told to try, and why both are wrong
//!
//! The brief proposed `w(v) = |v|^2` and a staircase order function
//! `w(u,v,w) = -(u + 2v + 4w)`. **Both are affine on `{0,1}^3` and neither
//! induces any triangulation at all**, and the harness measures that rather
//! than asserting it. `|v|^2 = v_0^2 + v_1^2 + v_2^2 = v_0 + v_1 + v_2` when
//! every coordinate is `0` or `1`; the staircase is affine by construction. An
//! affine lifting makes all eight points coplanar in `R^4`, every one of the 24
//! inequalities evaluate to **exactly 0**, and the lower hull a single
//! non-simplicial cell: `hull_facets = 0`. `w(v) = |v|^2` is the *Delaunay*
//! lifting, and the eight corners of a cube are cospherical, which is precisely
//! the degenerate case — so the classical recollection is right about the
//! function and wrong about the point set.
//!
//! Seven candidates are therefore run through the same four instruments, and
//! four of them are controls:
//!
//! | candidate | `w` | expected | why it is in the table |
//! |---|---|---|---|
//! | `neg_pairs` | `-|S|(|S|-1)/2` | accept | the exhibited lifting; also the LP's own witness |
//! | `neg_card_squared` | `-|S|^2` | accept | the family is `-(concave in |S|)`; shows the witness is not unique |
//! | `k3_cut` | `|S|(3-|S|)` | accept | the cut function of `K_3`, the textbook strictly submodular one |
//! | `abs_squared` | `|v|^2` | **reject** | affine, hence modular; the brief's first suggestion |
//! | `staircase_order` | `-(u+2v+4w)` | **reject** | affine, hence modular; the brief's second suggestion |
//! | `pos_pairs` | `+|S|(|S|-1)/2` | **reject** | strictly *super*modular; its extension is concave |
//! | `one_diamond_flipped` | `(0,0,0,0,0,-1,-1,-3)` | **reject** | the sharp control: still **convex**, violates exactly **one** of the nine, and merges two tetrahedra — `hull_facets = 4` |
//!
//! `one_diamond_flipped` is the control that matters. A checker that only
//! detected affineness would pass it.
//!
//! # C2 — the two instruments certify different surfaces, and the corpus says how often that shows
//!
//! T-015's Plantinga--Vegter checker (`validate/isotopy.rs`, `cell_is_certified`
//! at `:126`, `isotopy_report` at `:188`) certifies that **the trilinear
//! interpolant's** zero set inside a cell is one component isotopic to a disc.
//! Viro certifies that **the piecewise-linear T-construction surface** is
//! isotopic to a real algebraic hypersurface with the cube as Newton polytope.
//! Two different surfaces, two different hypotheses. The registered comparison
//! is therefore between two per-cell predicates:
//!
//! - **PW** — the transfer's hypotheses hold for this cell: the triangulation is
//!   regular (global, C1) **and** unimodular (global, asserted) **and** the sign
//!   distribution is *defined*, i.e. no corner value is exactly zero. Viro's
//!   `epsilon_omega` lives in `{+1, -1}`; a zero is not a sign, the T-surface
//!   passes through a lattice point, and the isotopy between midpoint gluing and
//!   the mesher's interpolated crossing degenerates. That is `M-48`'s degenerate
//!   crossing, and it is a genuine per-cell hypothesis failure rather than a
//!   technicality.
//! - **PV** — `cell_is_certified(&corner)`.
//!
//! `isotopy_agreement` is the share of **active** cells where the two verdicts
//! coincide, and active is the only honest denominator: an inactive cell passes
//! PV trivially by clause one (`isotopy.rs:123-124`) and carries an empty
//! T-surface, so agreement there is agreement on a constant. `pv_disagreements`
//! splits into `pw_yes_pv_no` and `pw_no_pv_yes`, and **both directions are
//! populated** — neither instrument dominates the other. `agreement_hypotheses_only`
//! records the softer reading in which PW is the global constant, so a reader can
//! see both denominators without re-running anything.
//!
//! Two further columns exist because they are where the divergence is *sharpest*,
//! and both were built from a constructed witness before the sweep:
//!
//! - `certified_multi_component` — PV-certified cells whose Kuhn PL surface has
//!   **more than one connected component**. `[1, -1, -1, 1, 4, 4, 4, 4]` is such
//!   a cell: its `z`-differences are `3, 5, 5, 3`, so the interval inner product
//!   is `3^2 - 2^2 - 2^2 = +1 > 0` and PV certifies "one component, isotopic to a
//!   disc", while the six tetrahedra emit two disjoint patches around corners `1`
//!   and `2`. With `4` replaced by `3` the `z`-range becomes `[2, 4]` and the sum
//!   is `2^2 - 2^2 - 2^2 = -4 < 0`: the same sign mask, *not* certified. That is
//!   the control proving the certification is a property of the values and not of
//!   the mask. So PV's certificate does not bound the component count of the
//!   surface marching tetrahedra actually emits — the clearest possible statement
//!   that the two instruments are not measuring one thing.
//! - `coeff_signs_match` / `t_components_differ` — Viro's `epsilon_omega` is the
//!   sign of the **coefficient of the monomial** `x^omega`, not the value at the
//!   corner. Those two sign vectors are related by the Möbius transform (the
//!   trilinear interpolant's monomial coefficients are the finite differences of
//!   its corner values, `c_omega = sum over v <= omega of (-1)^(|omega|-|v|) f_v`),
//!   and they are **not** equal. The columns count the cells where the two
//!   readings of "signs at the eight lattice points" give the same sign vector,
//!   and where they give T-surfaces with different component counts. This is the
//!   price of the source document's identification, measured rather than assumed.
//!
//! The component count comes from a **256-case table derived at run time** from
//! `TETS x TET_EDGES` — union-find over the cut cell-local edges, one class per
//! connected patch. Its 19 cell-local edges are classified by
//! `popcount(a ^ b)` and asserted against `P-124`'s census exactly —
//! **12** cube edges, **6** face diagonals, **1** body diagonal — and every
//! tetrahedron is asserted to cut `0`, `3` or `4` of its six edges, which is the
//! statement that a tetrahedron has **no ambiguity** (`table.rs:56-61`).
//!
//! # C3 — the Lovász identity, and the strict/weak distinction that decides it
//!
//! For `f: 2^[3] -> R`, the **Lovász extension** sorts the coordinates
//! `x_{i1} >= x_{i2} >= x_{i3}` and returns
//! `f(empty) + sum_k x_{ik} (f({i1..ik}) - f({i1..i(k-1)}))`. The region where a
//! given ordering holds is `{x : x_{i1} >= x_{i2} >= x_{i3}}`, whose vertices are
//! `empty, {i1}, {i1,i2}, [3]` — **a Kuhn tetrahedron**, and the four barycentric
//! coordinates are `(1 - x_{i1}, x_{i1} - x_{i2}, x_{i2} - x_{i3}, x_{i3})`. So
//! the Lovász extension *is* PL interpolation over exactly this triangulation,
//! and the harness checks it on `13^3 = 2197` exact rational probes per candidate
//! by two routes that share no code: the sorting formula, and barycentric
//! location by Cramer's rule. Everything is integer because a probe at `q/12`
//! with integer heights gives `12 * value` in `Z`.
//!
//! Then the algebra. The 9 distinct regularity inequalities turn out to be, on
//! the nose, **six local strict-submodularity diamonds** —
//! `f(S+i) + f(S+j) > f(S+i+j) + f(S)` for each of the 3 pairs and each of the 2
//! sets in the complement — plus **three complementary-pair conditions** of the
//! form `f(A) + f(B) > f(A cup B) + f(A cap B)` with `A cap B` empty. The six are
//! literally among the nine, and each of the remaining three is found by search
//! to be the **sum of exactly two** of them: nine Farkas certificates over
//! non-negative integer coefficients, so the two systems define the same open
//! cone and the implication is proved in both directions with no second solver.
//!
//! That gives the identity in its correct form, and the correct form is not the
//! one usually quoted:
//!
//! ```text
//! w submodular (weakly)         <=> the Lovasz extension of w is convex     (Lovasz 1983)
//! w STRICTLY submodular         <=> the lower hull of w is exactly Kuhn's six tetrahedra
//! ```
//!
//! The **weak** version is Lovász's theorem and is confirmed on all seven
//! candidates. The **strict** version is Viro's hypothesis. The gap between them
//! is not pedantry: it is exactly `abs_squared`, `staircase_order` and
//! `one_diamond_flipped` — three liftings that are convex and submodular and
//! induce a subdivision that is **not** the Kuhn triangulation. Both equivalences
//! are recorded as `7/7` agreement counts, so a single candidate breaking either
//! one is visible in the CSV.
//!
//! # Arms
//!
//! | arm | what it varies | is_control |
//! |---|---|---|
//! | exact Fourier--Motzkin on the 24-row system | nothing; one algebraic decision | no |
//! | the same solver on `w_0 > 0 && -w_0 > 0` | a system known infeasible | **yes** |
//! | seven named lifting candidates | the closed form of `w` | **four of seven** |
//! | brute-force lower-hull enumeration | all 70 four-subsets | no |
//! | eight reference fields at `33` and `65` | field and resolution | no |
//! | four constructed eight-value cells | the PW/PV verdict pair | **yes, all four** |
//! | `2197` exact rational probes per candidate | the Lovász identity's two routes | no |
//!
//! `33` and `65` are the contract's default pair and also two of the three
//! resolutions `docs/measurements/isotopy.csv` already carries, so the
//! `active_cells` and `uncertified_cells` columns are directly comparable with a
//! measurement that predates this row.
//!
//! # SHARE, recomputed before the numbers
//!
//! The registration's SHARE line is **`none` — this is a theorem acquisition**,
//! and it stays none after the fact. Nothing here proposes a source change, moves
//! a stage or touches a golden hash: the whole row is an algebraic decision plus
//! two read-only sweeps through the public API. If C1 holds, what the crate gains
//! is a *sentence* — that marching tetrahedra's output is, for every
//! non-degenerate sign configuration, isotopic to a genuine real algebraic
//! hypersurface with the cube as Newton polytope — and a citation for it. That
//! sentence is worth having and is worth exactly zero microseconds.
//!
//! # Vacuity controls
//!
//! `M-44`: a zero that could not have been non-zero is not a measurement. Every
//! one runs **before** the first `run.record` and every panic starts `VOID: `.
//!
//! | zero at risk | control, asserted | column that proves the fixture could have failed |
//! |---|---|---|
//! | `is_regular`, and the solver's ability to say no | the same Fourier--Motzkin run on `w_0 > 0 && -w_0 > 0` must report **infeasible** | `solver_refuses_infeasible` |
//! | `liftings_rejected` | at least one candidate accepted **and** at least one rejected; `abs_squared`, `staircase_order` and `one_diamond_flipped` named individually | `liftings_accepted`, `liftings_rejected` |
//! | `pv_disagreements` | the registration's own control: `uncertified_cells` summed over the corpus must exceed zero | `uncertified_cells` |
//! | `isotopy_agreement`'s denominator | `active_cells > 0` on every row, and `certified_cells > 0` somewhere | `active_cells`, `certified_cells` |
//! | `pw_no_pv_yes` | the `degenerate_certified` fixture — active, one corner exactly `0.0`, PV-certified — must classify as `pw = false, pv = true` | `fixtures_fired` |
//! | `pw_yes_pv_no` | the `alternating_uncertified` fixture must classify as `pw = true, pv = false` | `fixtures_fired` |
//! | `certified_multi_component` | the `multi_component` fixture must be certified with **2** components, and its `M = 3` twin must be **uncertified** on the same sign mask | `fixtures_fired` |
//! | the cell enumeration itself | `isotopy_report`'s counts must equal this harness's own walk on every row, or the two are counting different cells | `report_matches_walk` |
//!
//! # Determinism
//!
//! One thread, no PRNG, no map iteration. C1 and C3 are exact integer arithmetic
//! throughout — `i64` for the geometry, `i128` for elimination, exact rationals
//! only in back-substitution. The sweep is `f64`, `z`/`y`/`x` with `x` innermost,
//! byte-identical to `isotopy.rs:203-213` and to `common::grid`. Sign
//! classification is a three-way `<`/`>` test, so an exact zero is its own class
//! and no float equality is compared anywhere. `wall_ms` is recorded for
//! bookkeeping and **gates nothing** — every clause in this row is an integer
//! count or an exact comparison, which is `M-280`'s lesson taken up front rather
//! than after a governor swing.

mod common;

use std::time::Instant;

use isomesh::fields::ReferenceField;
use isomesh::marching_tetrahedra::table::{TET_EDGES, TETS};
use isomesh::validate::{cell_is_certified, isotopy_report};
use isomesh::{Sdf, Shape3};

/// Samples per axis. Two resolutions, which is what C2 asks for; both are also
/// rows of `docs/measurements/isotopy.csv`, so the sweep is checkable against a
/// measurement that predates it.
const RESOLUTIONS: [u32; 2] = [33, 65];

/// Corners of a cube.
const CORNERS: usize = 8;

/// Distinct cell-local edges of the Kuhn complex: 12 cube edges, 6 face
/// diagonals, 1 body diagonal. `P-124`'s census (`experiment_p124.rs:100-110`).
const CELL_EDGES: usize = 19;

/// Sign configurations of a cube.
const CASES: usize = 256;

/// Tetrahedra per cube, and edges per tetrahedron, from the shipped table.
const KUHN_TETS: usize = 6;

/// Denominator of the exact rational lattice the Lovász identity is probed on.
const PROBE_DEN: i64 = 12;

/// Denominator of the convexity lattice: endpoints on `{0, 2, 4}/4`, so every
/// midpoint lands on `{0, 1, 2, 3, 4}/4` and stays integral.
const MIDPOINT_DEN: i64 = 4;

// ─── cube geometry ──────────────────────────────────────────────────────────

/// The local offset of a cube corner, as grid steps.
///
/// `crate::cube::corner_offset` is `pub(crate)`
/// (`crates/isomesh/src/cube.rs:149-155`), so the three lines are copied here
/// with the source they came from rather than paraphrased from a comment.
const fn corner_offset(corner: u8) -> [u32; 3] {
    [
        (corner & 1) as u32,
        ((corner >> 1) & 1) as u32,
        ((corner >> 2) & 1) as u32,
    ]
}

/// The corner's position in the unit cube, as exact integers.
fn corner_point(corner: u8) -> [i64; 3] {
    let o = corner_offset(corner);
    [i64::from(o[0]), i64::from(o[1]), i64::from(o[2])]
}

/// `corner_offset` as `usize` grid steps, built once per sweep.
fn corner_steps() -> [[usize; 3]; CORNERS] {
    let mut out = [[0usize; 3]; CORNERS];
    for (i, slot) in out.iter_mut().enumerate() {
        let corner = u8::try_from(i).expect("eight corners fit u8");
        let o = corner_offset(corner);
        *slot = [
            usize::try_from(o[0]).expect("a corner step is 0 or 1"),
            usize::try_from(o[1]).expect("a corner step is 0 or 1"),
            usize::try_from(o[2]).expect("a corner step is 0 or 1"),
        ];
    }
    out
}

/// Three-way sign class: `-1`, `0` or `+1`.
///
/// Written with `<` and `>` and never `==`, so an exact zero is its own class
/// without a float equality anywhere. Viro's `epsilon` lives in `{+1, -1}`, so
/// the `0` class is precisely "this lattice point has no sign".
fn sign_class(v: f64) -> i8 {
    if v < 0.0 {
        -1
    } else if v > 0.0 {
        1
    } else {
        0
    }
}

// ─── exact integer linear algebra over the eight corners ────────────────────

/// Determinant of a 3x3 integer matrix.
fn det3(m: &[[i64; 3]; 3]) -> i64 {
    m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
        - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
        + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
}

/// Determinant of a 4x4 integer matrix, expanded along the first row.
fn det4(m: &[[i64; 4]; 4]) -> i64 {
    /// Columns kept in the minor, one entry per column of the first row.
    const REST: [[usize; 3]; 4] = [[1, 2, 3], [0, 2, 3], [0, 1, 3], [0, 1, 2]];
    let mut total = 0i64;
    for (col, keep) in REST.iter().enumerate() {
        let minor = [
            [m[1][keep[0]], m[1][keep[1]], m[1][keep[2]]],
            [m[2][keep[0]], m[2][keep[1]], m[2][keep[2]]],
            [m[3][keep[0]], m[3][keep[1]], m[3][keep[2]]],
        ];
        let sign = if col % 2 == 0 { 1 } else { -1 };
        total += sign * m[0][col] * det3(&minor);
    }
    total
}

/// The rows `[x, y, z, 1]` of four cube corners.
fn simplex_matrix(corners: &[u8; 4]) -> [[i64; 4]; 4] {
    let mut m = [[0i64; 4]; 4];
    for (row, corner) in m.iter_mut().zip(corners) {
        let p = corner_point(*corner);
        *row = [p[0], p[1], p[2], 1];
    }
    m
}

/// Cramer's rule for the affine combination of `corners` equal to `q / scale`.
///
/// Returns `(det, num)` with `scale * lambda_i = num[i] / det`, so the affine
/// interpolant of heights `w` at `q / scale` is `sum(num[i] * w_i) / (det * scale)`.
/// Every entry is a `0/1` corner coordinate or a lattice probe, so this is exact
/// in `i64`.
fn cramer(corners: &[u8; 4], q: [i64; 3], scale: i64) -> (i64, [i64; 4]) {
    let m = simplex_matrix(corners);
    let det = det4(&m);
    let mut num = [0i64; 4];
    for (i, slot) in num.iter_mut().enumerate() {
        let mut replaced = m;
        replaced[i] = [q[0], q[1], q[2], scale];
        *slot = det4(&replaced);
    }
    (det, num)
}

// ─── the regularity system ──────────────────────────────────────────────────

/// One strict homogeneous inequality `c . w > 0` over the eight corner heights.
#[derive(Clone, Copy, Debug)]
struct Constraint {
    /// Integer coefficients, one per cube corner.
    c: [i64; CORNERS],
    /// The Kuhn tetrahedron whose supporting affine function this is.
    tet: usize,
    /// The corner outside that tetrahedron which must sit strictly above it.
    outside: u8,
}

/// The 24 strict inequalities that say "the lower hull of `w` is exactly `TETS`".
///
/// For every Kuhn tetrahedron `T` and every cube corner `v` outside it, `w(v)`
/// must sit strictly above the affine function interpolating `w` on `T`.
fn regularity_system() -> Vec<Constraint> {
    let mut out = Vec::with_capacity(KUHN_TETS * 4);
    for (tet_index, tet) in TETS.iter().copied().enumerate() {
        for outside in 0..8u8 {
            if tet.contains(&outside) {
                continue;
            }
            let (det, num) = cramer(&tet, corner_point(outside), 1);
            assert_eq!(
                det.abs(),
                1,
                "P-139: Kuhn tetrahedron {tet:?} has normalised volume {}, not 1 — \
                 the barycentric coordinates are not integers and this whole clause \
                 is built on their being integers",
                det.abs()
            );
            let sgn = det.signum();
            let mut c = [0i64; CORNERS];
            c[usize::from(outside)] += 1;
            for (n, corner) in num.iter().zip(&tet) {
                c[usize::from(*corner)] -= sgn * n;
            }
            out.push(Constraint {
                c,
                tet: tet_index,
                outside,
            });
        }
    }
    out
}

/// Strict submodularity's six local diamonds on `2^[3]`.
///
/// `f(S + i) + f(S + j) > f(S + i + j) + f(S)` for each of the three pairs
/// `{i, j}` and each of the two subsets of the complement. Local submodularity
/// on these six is equivalent to submodularity on all pairs of sets, and C3
/// proves that equivalence here with Farkas certificates rather than citing it.
fn submodular_diamonds() -> Vec<[i64; CORNERS]> {
    let mut out = Vec::with_capacity(6);
    for (i, j) in [(0usize, 1usize), (0, 2), (1, 2)] {
        let (bi, bj) = (1usize << i, 1usize << j);
        let complement = 7 & !(bi | bj);
        for base in [0usize, complement] {
            let mut c = [0i64; CORNERS];
            c[base] -= 1;
            c[base | bi] += 1;
            c[base | bj] += 1;
            c[base | bi | bj] -= 1;
            out.push(c);
        }
    }
    out
}

/// Evaluate `c . w`.
fn dot(c: &[i64; CORNERS], w: &[i64; CORNERS]) -> i64 {
    c.iter().zip(w).map(|(a, b)| a * b).sum()
}

/// The first regularity constraint `w` fails, as `(tetrahedron, outside corner)`.
///
/// `None` means `w` is a convex lifting inducing exactly `TETS`. For a rejected
/// candidate the pair names the **wall that broke**, which is what distinguishes
/// `one_diamond_flipped` — one wall, two tetrahedra merged — from an affine
/// lifting, where the very first constraint already reads `0 > 0`.
fn first_violation(system: &[Constraint], w: &[i64; CORNERS]) -> Option<(usize, u8)> {
    system
        .iter()
        .find(|constraint| dot(&constraint.c, w) <= 0)
        .map(|constraint| (constraint.tet, constraint.outside))
}

// ─── exact Fourier--Motzkin ─────────────────────────────────────────────────

/// A homogeneous strict inequality `row . w > 0`, widened for elimination.
type Row = [i128; CORNERS];

/// Greatest common divisor, non-negative.
fn gcd_i128(a: i128, b: i128) -> i128 {
    let (mut a, mut b) = (a.abs(), b.abs());
    while b != 0 {
        let t = a % b;
        a = b;
        b = t;
    }
    a
}

/// Divide a row by the gcd of its entries, so equal half-spaces compare equal.
fn primitive(row: Row) -> Row {
    let g = row.iter().fold(0i128, |acc, v| gcd_i128(acc, *v));
    if g == 0 {
        return row;
    }
    let mut out = row;
    for slot in &mut out {
        *slot /= g;
    }
    out
}

/// One exact Fourier--Motzkin elimination run over a homogeneous strict system.
#[derive(Clone, Debug)]
struct Elimination {
    /// The system as it stood **before** eliminating `order[k]`, for each `k`,
    /// plus the fully projected system at the end. Back-substitution reads these.
    stages: Vec<Vec<Row>>,
    /// No stage ever produced an all-zero row, i.e. never `0 > 0`.
    feasible: bool,
    /// Largest stage, which is what would make this method unaffordable.
    max_rows: usize,
}

/// Project a homogeneous strict system by eliminating variables in `order`.
///
/// The system carries no constant column, because every regularity row is an
/// affine dependency and so sums to zero. An all-zero row therefore reads
/// `0 > 0` and is exactly the infeasibility certificate.
fn fourier_motzkin(rows: &[Row], order: &[usize; CORNERS]) -> Elimination {
    let mut current: Vec<Row> = rows.iter().copied().map(primitive).collect();
    current.sort_unstable();
    current.dedup();
    let mut stages = Vec::with_capacity(CORNERS + 1);
    let mut max_rows = current.len();

    for v in order {
        stages.push(current.clone());
        if current.iter().any(|r| r.iter().all(|x| *x == 0)) {
            return Elimination {
                stages,
                feasible: false,
                max_rows,
            };
        }
        let mut next: Vec<Row> = current.iter().copied().filter(|r| r[*v] == 0).collect();
        for upper in current.iter().filter(|r| r[*v] > 0) {
            for lower in current.iter().filter(|r| r[*v] < 0) {
                let (a, b) = (upper[*v], -lower[*v]);
                let mut combined = [0i128; CORNERS];
                for (slot, (u, l)) in combined.iter_mut().zip(upper.iter().zip(lower)) {
                    *slot = b * u + a * l;
                }
                if combined.iter().all(|x| *x == 0) {
                    stages.push(vec![combined]);
                    return Elimination {
                        stages,
                        feasible: false,
                        max_rows,
                    };
                }
                next.push(primitive(combined));
            }
        }
        next.sort_unstable();
        next.dedup();
        max_rows = max_rows.max(next.len());
        current = next;
    }

    let feasible = !current.iter().any(|r| r.iter().all(|x| *x == 0));
    stages.push(current);
    Elimination {
        stages,
        feasible,
        max_rows,
    }
}

/// An exact rational over `i128`, in lowest terms with a positive denominator.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
struct Rat {
    /// Numerator.
    num: i128,
    /// Denominator, always positive.
    den: i128,
}

impl Rat {
    /// Reduce and normalise the sign.
    fn new(num: i128, den: i128) -> Self {
        assert!(den != 0, "P-139: a rational with a zero denominator");
        let g = gcd_i128(num, den).max(1);
        let sign = if den < 0 { -1 } else { 1 };
        Self {
            num: sign * num / g,
            den: sign * den / g,
        }
    }

    /// An integer as a rational.
    fn int(n: i128) -> Self {
        Self { num: n, den: 1 }
    }

    /// Sum.
    fn add(self, other: Self) -> Self {
        Self::new(
            self.num * other.den + other.num * self.den,
            self.den * other.den,
        )
    }

    /// Product with an integer.
    fn scale(self, k: i128) -> Self {
        Self::new(self.num * k, self.den)
    }

    /// Quotient by a non-zero integer.
    fn divide(self, k: i128) -> Self {
        Self::new(self.num, self.den * k)
    }

    /// Negation.
    fn negate(self) -> Self {
        Self {
            num: -self.num,
            den: self.den,
        }
    }

    /// Midpoint of two rationals.
    fn midpoint(self, other: Self) -> Self {
        self.add(other).divide(2)
    }

    /// `self` compared with `other`, by cross-multiplication.
    fn is_greater(self, other: Self) -> bool {
        self.num * other.den > other.num * self.den
    }
}

/// Build a strictly feasible point by back-substitution through `elim`'s stages.
///
/// Each variable is placed strictly between the tightest lower and upper bounds
/// the stage imposes: `lo + 1` or `hi - 1` when one side is open, the midpoint
/// when both are, and `0` when the variable is unconstrained. Deterministic —
/// `lo` and `hi` are a max and a min, so the row order cannot change the answer.
fn back_substitute(elim: &Elimination, order: &[usize; CORNERS]) -> Option<[Rat; CORNERS]> {
    let mut w = [Rat::int(0); CORNERS];
    for index in (0..CORNERS).rev() {
        let v = order[index];
        let mut lo: Option<Rat> = None;
        let mut hi: Option<Rat> = None;
        for row in &elim.stages[index] {
            let mut rest = Rat::int(0);
            for (j, coefficient) in row.iter().enumerate() {
                if j != v {
                    rest = rest.add(w[j].scale(*coefficient));
                }
            }
            let c = row[v];
            if c == 0 {
                if !rest.is_greater(Rat::int(0)) {
                    return None;
                }
                continue;
            }
            let bound = rest
                .negate()
                .divide(c.abs())
                .scale(if c > 0 { 1 } else { -1 });
            if c > 0 {
                lo = Some(match lo {
                    Some(current) if current.is_greater(bound) => current,
                    _ => bound,
                });
            } else {
                hi = Some(match hi {
                    Some(current) if bound.is_greater(current) => current,
                    _ => bound,
                });
            }
        }
        w[v] = match (lo, hi) {
            (None, None) => Rat::int(0),
            (Some(l), None) => l.add(Rat::int(1)),
            (None, Some(h)) => h.add(Rat::int(-1)),
            (Some(l), Some(h)) => l.midpoint(h),
        };
    }
    Some(w)
}

/// Clear denominators, so a rational witness can be printed as integers.
fn clear_denominators(w: &[Rat; CORNERS]) -> [i64; CORNERS] {
    let mut lcm = 1i128;
    for r in w {
        lcm = lcm / gcd_i128(lcm, r.den) * r.den;
    }
    let mut out = [0i64; CORNERS];
    for (slot, r) in out.iter_mut().zip(w) {
        let scaled = r.num * (lcm / r.den);
        *slot = i64::try_from(scaled).expect("the witness fits i64");
    }
    out
}

// ─── candidate liftings ─────────────────────────────────────────────────────

/// A named candidate convex lifting, so a rejection can be reported by name.
#[derive(Clone, Copy, Debug)]
struct Lifting {
    /// Short token, used in `rejected_liftings`.
    name: &'static str,
    /// The closed form, for the header and the stanza.
    closed_form: &'static str,
    /// The eight heights.
    w: [i64; CORNERS],
    /// What this candidate is in the table for.
    expected_accept: bool,
}

/// Number of set bits of a corner index, i.e. `|S|`.
fn cardinality(corner: usize) -> i64 {
    i64::from(corner.count_ones())
}

/// The seven candidates: three that should induce Kuhn's triangulation and four
/// controls, two of which are the liftings this row was told to try.
fn candidates() -> Vec<Lifting> {
    let by_cardinality = |f: fn(i64) -> i64| {
        let mut w = [0i64; CORNERS];
        for (corner, slot) in w.iter_mut().enumerate() {
            *slot = f(cardinality(corner));
        }
        w
    };
    let mut abs_squared = [0i64; CORNERS];
    let mut staircase = [0i64; CORNERS];
    for corner in 0..CORNERS {
        let c = u8::try_from(corner).expect("eight corners fit u8");
        let p = corner_point(c);
        abs_squared[corner] = p[0] * p[0] + p[1] * p[1] + p[2] * p[2];
        staircase[corner] = -(p[0] + 2 * p[1] + 4 * p[2]);
    }

    vec![
        Lifting {
            name: "neg_pairs_of_set_bits",
            closed_form: "w(S)=-|S|*(|S|-1)/2",
            w: by_cardinality(|k| -(k * (k - 1)) / 2),
            expected_accept: true,
        },
        Lifting {
            name: "neg_cardinality_squared",
            closed_form: "w(S)=-|S|^2",
            w: by_cardinality(|k| -k * k),
            expected_accept: true,
        },
        Lifting {
            name: "k3_cut_function",
            closed_form: "w(S)=|S|*(3-|S|)",
            w: by_cardinality(|k| k * (3 - k)),
            expected_accept: true,
        },
        Lifting {
            name: "abs_squared",
            closed_form: "w(v)=|v|^2",
            w: abs_squared,
            expected_accept: false,
        },
        Lifting {
            name: "staircase_order",
            closed_form: "w(u.v.w)=-(u+2v+4w)",
            w: staircase,
            expected_accept: false,
        },
        Lifting {
            name: "pos_pairs_of_set_bits",
            closed_form: "w(S)=+|S|*(|S|-1)/2",
            w: by_cardinality(|k| (k * (k - 1)) / 2),
            expected_accept: false,
        },
        Lifting {
            name: "one_diamond_flipped",
            closed_form: "the exhibited lifting with w(3) raised to 0",
            w: [0, 0, 0, 0, 0, -1, -1, -3],
            expected_accept: false,
        },
    ]
}

/// Every simplicial lower facet of the lifted eight corners, and their total
/// normalised volume.
///
/// Brute force over all 70 four-subsets, each tested for "every other lifted
/// corner sits strictly above my affine hull". This finds **simplicial** facets
/// only, which is the right instrument rather than a limitation: a subdivision
/// with a non-simplicial cell reports fewer facets, and an affine lifting — one
/// flat cell — reports none at all.
fn lower_facets(w: &[i64; CORNERS]) -> (Vec<[u8; 4]>, i64) {
    let mut found = Vec::new();
    let mut volume = 0i64;
    for a in 0..8u8 {
        for b in a + 1..8u8 {
            for c in b + 1..8u8 {
                for d in c + 1..8u8 {
                    let subset = [a, b, c, d];
                    let det = det4(&simplex_matrix(&subset));
                    if det == 0 {
                        continue;
                    }
                    let sgn = det.signum();
                    let strict = (0..8u8).filter(|v| !subset.contains(v)).all(|v| {
                        let (_, num) = cramer(&subset, corner_point(v), 1);
                        let sum: i64 = num
                            .iter()
                            .zip(&subset)
                            .map(|(n, corner)| n * w[usize::from(*corner)])
                            .sum();
                        sgn * (det * w[usize::from(v)] - sum) > 0
                    });
                    if strict {
                        found.push(subset);
                        volume += det.abs();
                    }
                }
            }
        }
    }
    (found, volume)
}

/// Is `w` affine on `{0,1}^3`? An affine lifting is *modular*, is weakly convex,
/// and induces no subdivision at all — which is why both of the liftings this row
/// was told to try are rejected.
fn is_affine(w: &[i64; CORNERS]) -> bool {
    let (a, bx, by, bz) = (w[0], w[1] - w[0], w[2] - w[0], w[4] - w[0]);
    (0..CORNERS).all(|corner| {
        let c = u8::try_from(corner).expect("eight corners fit u8");
        let p = corner_point(c);
        w[corner] == a + bx * p[0] + by * p[1] + bz * p[2]
    })
}

// ─── the Lovász extension, two ways ─────────────────────────────────────────

/// `scale * (the Lovász extension of `w` at `q / scale`)`, by the sorting formula.
fn lovasz_scaled(w: &[i64; CORNERS], q: [i64; 3], scale: i64) -> i64 {
    let mut axis = [0usize, 1, 2];
    axis.sort_by_key(|d| (-q[*d], *d));
    let mut total = scale * w[0];
    let mut set = 0usize;
    for d in axis {
        let previous = set;
        set |= 1 << d;
        total += q[d] * (w[set] - w[previous]);
    }
    total
}

/// `scale * (the PL interpolant of `w` over the Kuhn triangulation at `q / scale`)`.
///
/// Located by barycentric coordinates rather than by sorting, so C3's two sides
/// are computed by routes that share no code. `None` means the probe fell in no
/// tetrahedron, which for a point of the unit cube would mean the six do not
/// tile it.
fn kuhn_pl_scaled(w: &[i64; CORNERS], q: [i64; 3], scale: i64) -> Option<i64> {
    for tet in &TETS {
        let (det, num) = cramer(tet, q, scale);
        let sgn = det.signum();
        if num.iter().all(|n| sgn * n >= 0) {
            let sum: i64 = num
                .iter()
                .zip(tet)
                .map(|(n, corner)| n * w[usize::from(*corner)])
                .sum();
            return Some(sgn * sum);
        }
    }
    None
}

/// Probe the Lovász identity on the whole `13^3` rational lattice.
///
/// Returns `(probes, mismatches)`. Exact: a probe at `q / 12` with integer
/// heights makes `12 * value` an integer on both sides.
fn lovasz_identity(w: &[i64; CORNERS]) -> (u64, u64) {
    let mut probes = 0u64;
    let mut mismatches = 0u64;
    for x in 0..=PROBE_DEN {
        for y in 0..=PROBE_DEN {
            for z in 0..=PROBE_DEN {
                let q = [x, y, z];
                let sorted = lovasz_scaled(w, q, PROBE_DEN);
                let barycentric = kuhn_pl_scaled(w, q, PROBE_DEN)
                    .expect("the six Kuhn tetrahedra tile the unit cube");
                probes += 1;
                if sorted != barycentric {
                    mismatches += 1;
                }
            }
        }
    }
    (probes, mismatches)
}

/// Midpoint convexity violations of the Lovász extension of `w`.
///
/// An independent, coarser instrument than "all nine inequalities hold weakly":
/// every unordered pair from the `{0, 1/2, 1}^3` lattice, tested for
/// `2 * f(mid) <= f(a) + f(b)` in integers.
fn midpoint_convexity_violations(w: &[i64; CORNERS]) -> u64 {
    let lattice: Vec<[i64; 3]> = (0..3)
        .flat_map(|x| (0..3).flat_map(move |y| (0..3).map(move |z| [2 * x, 2 * y, 2 * z])))
        .collect();
    let mut violations = 0u64;
    for (i, a) in lattice.iter().enumerate() {
        for b in &lattice[i + 1..] {
            let mid = [(a[0] + b[0]) / 2, (a[1] + b[1]) / 2, (a[2] + b[2]) / 2];
            let va = lovasz_scaled(w, *a, MIDPOINT_DEN);
            let vb = lovasz_scaled(w, *b, MIDPOINT_DEN);
            let vm = lovasz_scaled(w, mid, MIDPOINT_DEN);
            if 2 * vm > va + vb {
                violations += 1;
            }
        }
    }
    violations
}

/// Farkas certificates: how many of the distinct regularity rows are
/// non-negative integer combinations of the six submodularity diamonds.
///
/// Search over coefficients in `{0, 1, 2}^6`, which is `729` combinations. A
/// non-negative combination is a proof of implication for strict inequalities,
/// so finding one per row proves `diamonds => regularity`; the reverse holds
/// because the six diamonds are literally among the distinct rows.
fn farkas_certificates(distinct: &[Row], diamonds: &[[i64; CORNERS]]) -> (usize, usize) {
    let mut certified = vec![false; distinct.len()];
    let mut total_weight = 0usize;
    for code in 0..729u32 {
        let mut coefficients = [0u32; 6];
        let mut rest = code;
        for slot in &mut coefficients {
            *slot = rest % 3;
            rest /= 3;
        }
        if coefficients.iter().all(|c| *c == 0) {
            continue;
        }
        let mut combined = [0i128; CORNERS];
        for (coefficient, diamond) in coefficients.iter().zip(diamonds) {
            for (slot, entry) in combined.iter_mut().zip(diamond) {
                *slot += i128::from(*coefficient) * i128::from(*entry);
            }
        }
        let reduced = primitive(combined);
        if let Some(position) = distinct.iter().position(|r| *r == reduced)
            && !certified[position]
        {
            certified[position] = true;
            total_weight += usize::try_from(coefficients.iter().sum::<u32>())
                .expect("a certificate weight is at most 12");
        }
    }
    (certified.iter().filter(|c| **c).count(), total_weight)
}

// ─── the T-construction inside one cube, as a 256-case table ────────────────

/// The distinct cell-local edges of the Kuhn complex, the corner-pair lookup, and
/// the `[axis, face diagonal, body diagonal]` census.
///
/// Derived from `TETS x TET_EDGES`, never transcribed. An edge's class is
/// `popcount(a ^ b)`, because a cube corner index *is* its `0/1` coordinate, so
/// `1` is an axis step, `2` a face diagonal and `3` the body diagonal.
fn cell_local_edges() -> (usize, [[u8; CORNERS]; CORNERS], [usize; 3]) {
    let mut lookup = [[u8::MAX; CORNERS]; CORNERS];
    let mut count = 0usize;
    let mut census = [0usize; 3];
    for tet in TETS.iter().copied() {
        for pair in &TET_EDGES {
            let a = usize::from(tet[usize::from(pair[0])]);
            let b = usize::from(tet[usize::from(pair[1])]);
            if lookup[a][b] == u8::MAX {
                let id = u8::try_from(count).expect("19 cell-local edges fit u8");
                lookup[a][b] = id;
                lookup[b][a] = id;
                count += 1;
                let class = usize::try_from((a ^ b).count_ones())
                    .expect("a cell-local edge steps 1, 2 or 3 axes");
                census[class - 1] += 1;
            }
        }
    }
    (count, lookup, census)
}

/// Union-find root of `x`, with path halving.
fn uf_find(parent: &mut [u8; CELL_EDGES], mut x: usize) -> usize {
    while usize::from(parent[x]) != x {
        let grandparent = parent[usize::from(parent[x])];
        parent[x] = grandparent;
        x = usize::from(grandparent);
    }
    x
}

/// Join the classes of `a` and `b`.
fn uf_union(parent: &mut [u8; CELL_EDGES], a: usize, b: usize) {
    let (ra, rb) = (uf_find(parent, a), uf_find(parent, b));
    if ra != rb {
        parent[ra] = u8::try_from(rb).expect("19 cell-local edges fit u8");
    }
}

/// What the T-construction emits inside one cube for one sign mask.
#[derive(Clone, Copy, Debug, Default)]
struct KuhnCase {
    /// Connected components of the piecewise-linear surface.
    components: u8,
    /// Cell-local edges the surface cuts.
    cut_edges: u8,
    /// Triangles the six tetrahedra emit.
    triangles: u8,
}

/// The T-construction for one sign mask, by union-find over cut edges.
///
/// A tetrahedron cuts `0`, `3` or `4` of its six edges and nothing else — the
/// statement that a tetrahedron carries **no ambiguity** (`table.rs:56-61`) — and
/// each cut set is one connected patch, a triangle or a fanned quad. So the
/// surface's components are the classes of the cut cell-local edges.
fn kuhn_case(mask: u8, lookup: &[[u8; CORNERS]; CORNERS]) -> KuhnCase {
    let inside = |corner: u8| mask & (1 << corner) != 0;
    let mut parent = [0u8; CELL_EDGES];
    for (i, slot) in parent.iter_mut().enumerate() {
        *slot = u8::try_from(i).expect("19 cell-local edges fit u8");
    }
    let mut touched = [false; CELL_EDGES];
    let mut triangles = 0u8;

    for tet in TETS.iter().copied() {
        let mut cut = [0usize; 4];
        let mut cuts = 0usize;
        for pair in &TET_EDGES {
            let a = tet[usize::from(pair[0])];
            let b = tet[usize::from(pair[1])];
            if inside(a) != inside(b) {
                assert!(
                    cuts < 4,
                    "P-139: tetrahedron {tet:?} cut more than four of its six edges \
                     at mask {mask:#010b}, which no linear interpolant can do"
                );
                cut[cuts] = usize::from(lookup[usize::from(a)][usize::from(b)]);
                cuts += 1;
            }
        }
        assert!(
            cuts == 0 || cuts == 3 || cuts == 4,
            "P-139: tetrahedron {tet:?} cut {cuts} of its six edges at mask \
             {mask:#010b} — a tetrahedron cuts 0, 3 or 4 and nothing else \
             (table.rs:56-61), so the case table is not the shipped complex"
        );
        triangles += match cuts {
            0 => 0,
            3 => 1,
            _ => 2,
        };
        for edge in &cut[..cuts] {
            touched[*edge] = true;
            uf_union(&mut parent, cut[0], *edge);
        }
    }

    let mut seen = [false; CELL_EDGES];
    let mut components = 0u8;
    let mut cut_edges = 0u8;
    for (edge, hit) in touched.iter().enumerate() {
        if !*hit {
            continue;
        }
        cut_edges += 1;
        let root = uf_find(&mut parent, edge);
        if !seen[root] {
            seen[root] = true;
            components += 1;
        }
    }

    KuhnCase {
        components,
        cut_edges,
        triangles,
    }
}

/// The whole 256-case table, derived once.
fn kuhn_case_table(lookup: &[[u8; CORNERS]; CORNERS]) -> [KuhnCase; CASES] {
    let mut table = [KuhnCase::default(); CASES];
    for (mask, slot) in table.iter_mut().enumerate() {
        let mask = u8::try_from(mask).expect("256 cases fit u8");
        *slot = kuhn_case(mask, lookup);
    }
    table
}

/// The eight monomial coefficients of the cell's trilinear interpolant.
///
/// The Möbius transform of the corner values,
/// `c_omega = sum over v <= omega of (-1)^(|omega| - |v|) f_v`, computed by the
/// three-pass butterfly. These are the signs Viro's `epsilon_omega` actually is —
/// the sign of the coefficient of `x^omega` — and they are **not** the corner
/// values' signs.
fn monomial_coefficients(f: &[f64; CORNERS]) -> [f64; CORNERS] {
    let mut c = *f;
    for bit in [1usize, 2, 4] {
        for corner in 0..CORNERS {
            if corner & bit != 0 {
                c[corner] -= c[corner ^ bit];
            }
        }
    }
    c
}

// ─── the corpus sweep ───────────────────────────────────────────────────────

/// One reference field at one resolution.
#[derive(Clone, Debug)]
struct Sweep {
    /// Reference-field name.
    field: &'static str,
    /// Samples per axis.
    samples: u32,
    /// Cells examined, `(samples - 1)^3`.
    cells: u64,
    /// Cells the surface passes through.
    active: u64,
    /// Active cells Plantinga--Vegter certifies.
    certified: u64,
    /// Active cells it does not.
    uncertified: u64,
    /// Active cells with a corner value of exactly zero, where Viro's sign
    /// distribution is undefined.
    exact_zero: u64,
    /// Active cells where patchworking's hypotheses hold and PV does not certify.
    pw_yes_pv_no: u64,
    /// Active cells where PV certifies and patchworking's hypotheses fail.
    pw_no_pv_yes: u64,
    /// Active cells where the two verdicts coincide.
    agree: u64,
    /// PV-certified cells whose Kuhn PL surface has more than one component.
    certified_multi: u64,
    /// Uncertified cells with more than one component — the control for the line
    /// above, so its zero could have been non-zero.
    uncertified_multi: u64,
    /// Active cells where the monomial-coefficient signs equal the corner signs.
    coeff_signs_match: u64,
    /// Active cells with a monomial coefficient of exactly zero.
    coeff_zero: u64,
    /// Active cells whose two sign readings give different component counts.
    components_differ: u64,
    /// `isotopy_report`'s counts equal this harness's own walk.
    report_matches: bool,
    /// Wall clock, bookkeeping only.
    wall_ms: f64,
}

impl Sweep {
    /// Share of active cells where the two instruments agree.
    fn agreement(&self) -> f64 {
        self.agree as f64 / self.active as f64
    }

    /// Share of active cells PV certifies — the softer reading of C2, in which
    /// patchworking's verdict is the global constant `true`.
    fn agreement_hypotheses_only(&self) -> f64 {
        self.certified as f64 / self.active as f64
    }

    /// Disagreements, both directions.
    fn disagreements(&self) -> u64 {
        self.pw_yes_pv_no + self.pw_no_pv_yes
    }
}

/// Sweep one field at one resolution.
fn measure<F>(field: &F, name: &'static str, samples: u32, table: &[KuhnCase; CASES]) -> Sweep
where
    F: ReferenceField + Sdf<Scalar = f64>,
{
    let started = Instant::now();
    let (shape, origin, h) = common::grid::<f64, _>(field, samples);
    let n = usize::try_from(samples).expect("a resolution fits usize");

    let mut values = Vec::with_capacity(shape.element_count());
    for z in 0..samples {
        for y in 0..samples {
            for x in 0..samples {
                values.push(field.sample([
                    origin[0] + h * f64::from(x),
                    origin[1] + h * f64::from(y),
                    origin[2] + h * f64::from(z),
                ]));
            }
        }
    }

    let report = isotopy_report(&values, &shape).expect("T-015 report over a 33 or 65 sample grid");
    let steps = corner_steps();

    let mut out = Sweep {
        field: name,
        samples,
        cells: 0,
        active: 0,
        certified: 0,
        uncertified: 0,
        exact_zero: 0,
        pw_yes_pv_no: 0,
        pw_no_pv_yes: 0,
        agree: 0,
        certified_multi: 0,
        uncertified_multi: 0,
        coeff_signs_match: 0,
        coeff_zero: 0,
        components_differ: 0,
        report_matches: false,
        wall_ms: 0.0,
    };

    for z in 0..n - 1 {
        for y in 0..n - 1 {
            for x in 0..n - 1 {
                let mut corner = [0.0f64; CORNERS];
                for (slot, step) in corner.iter_mut().zip(&steps) {
                    let index = ((z + step[2]) * n + (y + step[1])) * n + (x + step[0]);
                    *slot = values[index];
                }
                out.cells += 1;

                let classes = corner.map(sign_class);
                let inside = classes[0] < 0;
                if classes.iter().all(|c| (*c < 0) == inside) {
                    continue;
                }
                out.active += 1;

                let pv = cell_is_certified(&corner);
                let has_zero = classes.contains(&0);
                let pw = !has_zero;
                if has_zero {
                    out.exact_zero += 1;
                }
                if pv {
                    out.certified += 1;
                } else {
                    out.uncertified += 1;
                }
                if pw == pv {
                    out.agree += 1;
                } else if pw {
                    out.pw_yes_pv_no += 1;
                } else {
                    out.pw_no_pv_yes += 1;
                }

                let mut mask = 0u8;
                for (i, class) in classes.iter().enumerate() {
                    if *class < 0 {
                        mask |= 1 << i;
                    }
                }
                let components = table[usize::from(mask)].components;
                if pv && components > 1 {
                    out.certified_multi += 1;
                } else if !pv && components > 1 {
                    out.uncertified_multi += 1;
                }

                let coefficients = monomial_coefficients(&corner);
                let coefficient_classes = coefficients.map(sign_class);
                if coefficient_classes.contains(&0) {
                    out.coeff_zero += 1;
                }
                if coefficient_classes == classes {
                    out.coeff_signs_match += 1;
                }
                let mut coefficient_mask = 0u8;
                for (i, class) in coefficient_classes.iter().enumerate() {
                    if *class < 0 {
                        coefficient_mask |= 1 << i;
                    }
                }
                if table[usize::from(coefficient_mask)].components != components {
                    out.components_differ += 1;
                }
            }
        }
    }

    out.report_matches = report.cells == out.cells
        && report.active_cells == out.active
        && report.certified == out.certified
        && report.uncertified == out.uncertified;
    out.wall_ms = started.elapsed().as_secs_f64() * 1e3;
    out
}

// ─── the constructed cells that prove each verdict pair can occur ───────────

/// A named eight-value cell built to make one column able to fire.
#[derive(Clone, Copy, Debug)]
struct Fixture {
    /// Short token for the report.
    name: &'static str,
    /// The eight corner values.
    corner: [f64; CORNERS],
    /// Expected `cell_is_certified`.
    certified: bool,
    /// Expected patchworking verdict, i.e. no corner value is exactly zero.
    patchworking: bool,
    /// Expected connected components of the Kuhn PL surface.
    components: u8,
}

/// The four constructed cells, one per verdict pair the corpus columns can hold.
const FIXTURES: [Fixture; 4] = [
    // Corner 0 is exactly zero, so Viro has no sign there, while the four
    // strictly positive z-differences give an interval inner product of 16 and
    // Plantinga--Vegter certifies. This is the only thing that licenses a
    // `pw_no_pv_yes` of zero.
    Fixture {
        name: "degenerate_certified",
        corner: [0.0, -1.0, -1.0, -1.0, 4.0, 4.0, 4.0, 4.0],
        certified: true,
        patchworking: false,
        components: 1,
    },
    // Certified, and the six tetrahedra emit two disjoint patches: the four
    // z-differences are 3, 5, 5, 3, so the interval inner product is
    // 3^2 - 2^2 - 2^2 = +1 > 0 while corners 1 and 2 are separated. PV's
    // certificate does not bound the component count of the surface marching
    // tetrahedra emits.
    Fixture {
        name: "multi_component_certified",
        corner: [1.0, -1.0, -1.0, 1.0, 4.0, 4.0, 4.0, 4.0],
        certified: true,
        patchworking: true,
        components: 2,
    },
    // The same sign mask with 4 replaced by 3: the z-range becomes [2, 4] and the
    // sum is 2^2 - 2^2 - 2^2 = -4 < 0, so it is **not** certified. Proves the line
    // above is a property of the values and not of the mask.
    Fixture {
        name: "multi_component_control",
        corner: [1.0, -1.0, -1.0, 1.0, 3.0, 3.0, 3.0, 3.0],
        certified: false,
        patchworking: true,
        components: 2,
    },
    // Alternating signs, no zero: Viro's hypotheses hold and Plantinga--Vegter
    // refuses. This is what licenses a `pw_yes_pv_no` of zero.
    Fixture {
        name: "alternating_uncertified",
        corner: [1.0, -1.0, 1.0, -1.0, -1.0, 1.0, -1.0, 1.0],
        certified: false,
        patchworking: true,
        components: 2,
    },
];

fn main() {
    if !std::env::args().any(|arg| arg == "--bench") {
        return;
    }
    let prereg = isomesh::experiment!("P-139");

    common::experiment::run(prereg, |run| {
        // ── the complex, checked against the shipped table ───────────────────
        let (edge_count, edge_of, census) = cell_local_edges();
        assert_eq!(
            edge_count, CELL_EDGES,
            "P-139: TETS x TET_EDGES reduced to {edge_count} distinct cell-local \
             edges, not 19 — the complex is not P-124's"
        );
        assert_eq!(
            census,
            [12, 6, 1],
            "P-139: the shipped six-tet split carries {census:?} cell-local edges \
             by class, not P-124's 12 cube edges / 6 face diagonals / 1 body \
             diagonal (experiment_p124.rs:100-110) — the complex is not the one \
             every claim below is about"
        );
        let table = kuhn_case_table(&edge_of);
        let max_components = table.iter().map(|c| c.components).max().unwrap_or(0);
        let max_triangles = table.iter().map(|c| c.triangles).max().unwrap_or(0);
        let max_cut_edges = table.iter().map(|c| c.cut_edges).max().unwrap_or(0);
        // A cut edge carries one surface vertex and a triangle three, but a
        // vertex is shared, so the two counts move together and neither can be
        // zero unless the other is: the table is not a lookup and this is the
        // cheapest check that it was built rather than guessed.
        assert!(
            table
                .iter()
                .all(|c| (c.cut_edges == 0) == (c.triangles == 0)
                    && (c.cut_edges == 0) == (c.components == 0)),
            "P-139: the derived 256-case table has a mask whose cut edges, \
             triangles and components disagree about whether the cell is active"
        );

        // ── unimodularity and the GKZ vector, both computed ──────────────────
        let mut gkz = [0i64; CORNERS];
        let mut normalised_volume = 0i64;
        let mut unimodular = true;
        for tet in TETS.iter().copied() {
            let det = det4(&simplex_matrix(&tet)).abs();
            unimodular &= det == 1;
            normalised_volume += det;
            for corner in tet {
                gkz[usize::from(corner)] += det;
            }
        }
        assert!(
            unimodular && normalised_volume == 6,
            "P-139: the six Kuhn tetrahedra have total normalised volume \
             {normalised_volume} and unimodular = {unimodular} — Viro's \
             T-construction needs a unimodular triangulation of the Newton \
             polytope, so this is a hypothesis of the transfer and not a detail"
        );
        let gkz_sum: i64 = gkz.iter().sum();
        assert_eq!(
            gkz_sum, 24,
            "P-139: the GKZ vector sums to {gkz_sum}, not (dim + 1) x volume = 24"
        );
        let secondary_vertex = gkz
            .iter()
            .map(std::string::ToString::to_string)
            .collect::<Vec<_>>()
            .join("|");

        // ── C1: the exact linear-feasibility decision ────────────────────────
        let system = regularity_system();
        assert_eq!(
            system.len(),
            24,
            "P-139: the regularity system has {} rows, not 6 tetrahedra x 4 \
             outside corners",
            system.len()
        );
        let mut distinct: Vec<Row> = system
            .iter()
            .map(|constraint| primitive(constraint.c.map(i128::from)))
            .collect();
        distinct.sort_unstable();
        distinct.dedup();

        let order = [7usize, 6, 5, 4, 3, 2, 1, 0];
        let elimination = fourier_motzkin(&distinct, &order);

        // ── VOID: the solver must be able to say no ──────────────────────────
        let mut positive = [0i128; CORNERS];
        positive[0] = 1;
        let mut negative = [0i128; CORNERS];
        negative[0] = -1;
        let refutation = fourier_motzkin(&[positive, negative], &order);
        assert!(
            !refutation.feasible,
            "VOID: the same Fourier-Motzkin run called `w_0 > 0 && -w_0 > 0` \
             feasible, so `is_regular = true` is a constant this solver would \
             have printed for any system and C1 is not a measurement (M-44)"
        );
        let solver_refuses_infeasible = !refutation.feasible;

        let witness = back_substitute(&elimination, &order);
        let witness_heights = witness.map(|w| clear_denominators(&w));
        let witness_verified = witness_heights.is_some_and(|heights| {
            system
                .iter()
                .all(|constraint| dot(&constraint.c, &heights) > 0)
        });

        // ── the candidate liftings, all seven through all four instruments ───
        let diamonds = submodular_diamonds();
        assert_eq!(
            diamonds.len(),
            6,
            "P-139: 3 pairs x 2 complements = 6 diamonds"
        );
        let kuhn_sorted = {
            let mut all: Vec<[u8; 4]> = TETS.to_vec();
            all.sort_unstable();
            all
        };

        let mut accepted = 0usize;
        let mut rejected: Vec<&'static str> = Vec::new();
        let mut weak_iff_convex = 0usize;
        let mut strict_iff_kuhn = 0usize;
        let mut lovasz_probes = 0u64;
        let mut lovasz_mismatches = 0u64;
        let mut exhibited: Option<Lifting> = None;
        let mut exhibited_facets = 0usize;
        let mut control_first_violation = String::from("none");

        println!(
            "{:<26} {:>34} {:>7} {:>7} {:>7} {:>7} {:>7} {:>7} {:>10}",
            "lifting",
            "closed form",
            "strict",
            "weak",
            "affine",
            "facets",
            "kuhn",
            "midpt",
            "first fail"
        );
        for lifting in candidates() {
            let strict_sub = diamonds.iter().all(|d| dot(d, &lifting.w) > 0);
            let weak_sub = diamonds.iter().all(|d| dot(d, &lifting.w) >= 0);
            let strict_reg = system.iter().all(|c| dot(&c.c, &lifting.w) > 0);
            let weak_reg = system.iter().all(|c| dot(&c.c, &lifting.w) >= 0);
            let (facets, volume) = lower_facets(&lifting.w);
            let mut facet_set = facets.clone();
            facet_set.sort_unstable();
            let induces_kuhn = facet_set == kuhn_sorted && volume == 6;
            let midpoint_violations = midpoint_convexity_violations(&lifting.w);
            let broke = first_violation(&system, &lifting.w);
            let broke_at = broke.map_or_else(
                || String::from("-"),
                |(tet, outside)| format!("t{tet}/v{outside}"),
            );
            if lifting.name == "one_diamond_flipped" {
                control_first_violation = broke_at.clone();
            }
            let (probes, mismatches) = lovasz_identity(&lifting.w);
            lovasz_probes += probes;
            lovasz_mismatches += mismatches;

            // Lovász 1983, in its weak form: the extension is convex exactly
            // when the set function is submodular. The midpoint probe is the
            // independent witness; `weak_reg` is the exact statement.
            if weak_sub == weak_reg && weak_sub == (midpoint_violations == 0) {
                weak_iff_convex += 1;
            }
            // Viro's hypothesis, in its strict form. This is the equivalence
            // that the two affine candidates break and the weak one does not.
            if strict_sub == strict_reg && strict_sub == induces_kuhn {
                strict_iff_kuhn += 1;
            }

            println!(
                "{:<26} {:>34} {:>7} {:>7} {:>7} {:>7} {:>7} {:>7} {:>10}",
                lifting.name,
                lifting.closed_form,
                strict_sub,
                weak_sub,
                is_affine(&lifting.w),
                facets.len(),
                induces_kuhn,
                midpoint_violations,
                broke_at
            );

            assert_eq!(
                induces_kuhn,
                lifting.expected_accept,
                "P-139: candidate `{}` ({}) was registered as expected_accept = \
                 {} and the lower hull says {} with {} simplicial facets, first \
                 failing constraint {broke_at} — the prediction is written into \
                 this file, so a mismatch is a result and not a bug to paper over",
                lifting.name,
                lifting.closed_form,
                lifting.expected_accept,
                induces_kuhn,
                facets.len()
            );
            assert_eq!(
                broke.is_none(),
                induces_kuhn,
                "P-139: candidate `{}` has first_violation = {broke_at} and \
                 induces_kuhn = {induces_kuhn} — the 24-inequality system and the \
                 brute-force lower hull are two routes to one answer and they \
                 must not disagree",
                lifting.name
            );

            if induces_kuhn {
                accepted += 1;
                if exhibited.is_none() {
                    exhibited = Some(lifting);
                    exhibited_facets = facets.len();
                }
            } else {
                rejected.push(lifting.name);
            }
        }
        println!();

        let candidate_count = accepted + rejected.len();

        // ── VOID: the lifting checker must be able to reject ─────────────────
        assert!(
            accepted > 0 && !rejected.is_empty(),
            "VOID: {accepted} of {candidate_count} candidate liftings were \
             accepted and {} rejected — an instrument that accepts everything or \
             rejects everything makes `is_regular` a constant (M-44)",
            rejected.len()
        );
        for name in ["abs_squared", "staircase_order", "one_diamond_flipped"] {
            assert!(
                rejected.contains(&name),
                "VOID: `{name}` was not rejected. `abs_squared` and \
                 `staircase_order` are affine on the 0/1 cube and induce no \
                 subdivision; `one_diamond_flipped` is convex and merges two \
                 tetrahedra. A checker that passes any of the three would pass \
                 anything convex, and C1's accept would mean nothing (M-44)"
            );
        }

        let exhibited = exhibited.expect("a lifting was accepted, so one was exhibited");
        let witness_matches_exhibited = witness_heights.is_some_and(|heights| {
            let difference = {
                let mut d = [0i64; CORNERS];
                for (slot, (a, b)) in d.iter_mut().zip(heights.iter().zip(&exhibited.w)) {
                    *slot = a - b;
                }
                d
            };
            is_affine(&difference)
        });

        // ── C3: the Farkas certificates that identify the two cones ──────────
        let (certificates, certificate_weight) = farkas_certificates(&distinct, &diamonds);
        let diamonds_are_rows = diamonds
            .iter()
            .all(|d| distinct.contains(&primitive(d.map(i128::from))));

        let c1 = elimination.feasible
            && witness_verified
            && witness_matches_exhibited
            && exhibited_facets == KUHN_TETS
            && unimodular;
        let c3 = lovasz_mismatches == 0
            && certificates == distinct.len()
            && diamonds_are_rows
            && weak_iff_convex == candidate_count
            && strict_iff_kuhn == candidate_count;

        let stage_sizes = elimination
            .stages
            .iter()
            .map(|s| s.len().to_string())
            .collect::<Vec<_>>()
            .join("|");
        // The header claims elimination causes no coefficient growth. That is a
        // measurement, so it goes in the CSV rather than in the prose: growth
        // would be information about the system, not a failure of the harness.
        let fm_max_coefficient = elimination
            .stages
            .iter()
            .flatten()
            .flatten()
            .map(|x| x.abs())
            .max()
            .unwrap_or(0);
        println!(
            "C1  fourier-motzkin over {} distinct rows (of {}), order {order:?}\n    \
             feasible {}   stages {stage_sizes}   max stage {}   witness {:?}\n    \
             exhibited {} = {}   hull facets {}   unimodular {}   GKZ {}\n\
             C3  lovasz probes {lovasz_probes} mismatches {lovasz_mismatches}   \
             certificates {certificates}/{}  (total weight {certificate_weight})\n    \
             weak submodular <=> convex {weak_iff_convex}/{candidate_count}   \
             strict submodular <=> kuhn {strict_iff_kuhn}/{candidate_count}\n",
            distinct.len(),
            system.len(),
            elimination.feasible,
            elimination.max_rows,
            witness_heights,
            exhibited.name,
            exhibited.closed_form,
            exhibited_facets,
            unimodular,
            secondary_vertex,
            distinct.len(),
        );

        // ── VOID: every verdict pair must be reachable in this run ───────────
        let mut fixtures_fired = 0usize;
        for fixture in FIXTURES {
            let classes = fixture.corner.map(sign_class);
            let inside = classes[0] < 0;
            let active = !classes.iter().all(|c| (*c < 0) == inside);
            let pv = cell_is_certified(&fixture.corner);
            let pw = !classes.contains(&0);
            let mut mask = 0u8;
            for (i, class) in classes.iter().enumerate() {
                if *class < 0 {
                    mask |= 1 << i;
                }
            }
            let components = table[usize::from(mask)].components;
            assert!(
                active
                    && pv == fixture.certified
                    && pw == fixture.patchworking
                    && components == fixture.components,
                "VOID: fixture `{}` was constructed to be active with \
                 certified = {}, patchworking = {} and {} components, and this \
                 run reads active = {active}, certified = {pv}, \
                 patchworking = {pw}, components = {components}. Without it the \
                 corpus column it licenses is a zero that could not have been \
                 non-zero (M-44)",
                fixture.name,
                fixture.certified,
                fixture.patchworking,
                fixture.components
            );
            fixtures_fired += 1;
        }
        println!(
            "fixtures: {fixtures_fired} of {} fired as constructed\n",
            FIXTURES.len()
        );

        // ── C2: the corpus sweep ─────────────────────────────────────────────
        let mut sweeps: Vec<Sweep> = Vec::new();
        println!(
            "{:<16} {:>5} {:>9} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>9}",
            "field",
            "n",
            "cells",
            "active",
            "cert",
            "uncert",
            "zero",
            "py/pn",
            "pn/py",
            "cert>1",
            "unc>1",
            "agree"
        );
        isomesh::for_each_reference_field!(f64, |name, field| {
            for samples in RESOLUTIONS {
                let sweep = measure(&field, name, samples, &table);
                println!(
                    "{:<16} {:>5} {:>9} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>9.6}",
                    sweep.field,
                    sweep.samples,
                    sweep.cells,
                    sweep.active,
                    sweep.certified,
                    sweep.uncertified,
                    sweep.exact_zero,
                    sweep.pw_yes_pv_no,
                    sweep.pw_no_pv_yes,
                    sweep.certified_multi,
                    sweep.uncertified_multi,
                    sweep.agreement()
                );
                sweeps.push(sweep);
            }
        });
        println!();

        // ── VOID: the registration's own control, and the denominators ───────
        let corpus_uncertified: u64 = sweeps.iter().map(|s| s.uncertified).sum();
        let corpus_certified: u64 = sweeps.iter().map(|s| s.certified).sum();
        assert!(
            corpus_uncertified > 0,
            "VOID: T-015 reported no isotopy failure anywhere in the corpus, so \
             C2's agreement is agreement on a constant — the registration's own \
             vacuity control"
        );
        assert!(
            corpus_certified > 0,
            "VOID: T-015 certified nothing in the corpus, so the agreement \
             fraction is a constant from the other side"
        );
        for sweep in &sweeps {
            assert!(
                sweep.active > 0,
                "VOID: {} at {}³ has no active cell, so `isotopy_agreement` has \
                 no denominator",
                sweep.field,
                sweep.samples
            );
            assert!(
                sweep.report_matches,
                "VOID: {} at {}³ — this harness walked {} cells / {} active / {} \
                 certified / {} uncertified and `isotopy_report` disagrees, so the \
                 two instruments are not counting the same cells and no comparison \
                 between them means anything",
                sweep.field,
                sweep.samples,
                sweep.cells,
                sweep.active,
                sweep.certified,
                sweep.uncertified
            );
        }

        let rejected_names = rejected.join("|");
        let witness_column = witness_heights.map_or_else(
            || String::from("none"),
            |heights| {
                heights
                    .iter()
                    .map(std::string::ToString::to_string)
                    .collect::<Vec<_>>()
                    .join("|")
            },
        );

        for sweep in &sweeps {
            let c2 = sweep.disagreements() == 0;
            run.record(&[
                // ── the eleven registered metrics ────────────────────────────
                ("triangulation", String::from("kuhn_freudenthal_six_tet")),
                ("lifting_function", exhibited.name.to_string()),
                ("is_regular", c1.to_string()),
                ("secondary_polytope_vertex", secondary_vertex.clone()),
                ("patchworking_applies", (sweep.exact_zero == 0).to_string()),
                ("isotopy_agreement", format!("{:.6}", sweep.agreement())),
                ("pv_disagreements", sweep.disagreements().to_string()),
                ("cells_checked", sweep.cells.to_string()),
                ("c1_holds", c1.to_string()),
                ("c2_holds", c2.to_string()),
                ("c3_holds", c3.to_string()),
                // ── extras (M-273) ───────────────────────────────────────────
                ("field", sweep.field.to_string()),
                ("resolution", sweep.samples.to_string()),
                ("active_cells", sweep.active.to_string()),
                ("certified_cells", sweep.certified.to_string()),
                ("uncertified_cells", sweep.uncertified.to_string()),
                ("exact_zero_corner_cells", sweep.exact_zero.to_string()),
                ("pw_yes_pv_no", sweep.pw_yes_pv_no.to_string()),
                ("pw_no_pv_yes", sweep.pw_no_pv_yes.to_string()),
                (
                    "agreement_hypotheses_only",
                    format!("{:.6}", sweep.agreement_hypotheses_only()),
                ),
                (
                    "certified_multi_component",
                    sweep.certified_multi.to_string(),
                ),
                (
                    "uncertified_multi_component",
                    sweep.uncertified_multi.to_string(),
                ),
                ("coeff_signs_match", sweep.coeff_signs_match.to_string()),
                ("coeff_zero_cells", sweep.coeff_zero.to_string()),
                ("t_components_differ", sweep.components_differ.to_string()),
                ("report_matches_walk", sweep.report_matches.to_string()),
                ("lifting_closed_form", exhibited.closed_form.to_string()),
                ("lifting_heights", witness_column.clone()),
                (
                    "witness_matches_exhibited",
                    witness_matches_exhibited.to_string(),
                ),
                ("witness_verified", witness_verified.to_string()),
                ("kuhn_unimodular", unimodular.to_string()),
                ("kuhn_normalised_volume", normalised_volume.to_string()),
                ("kuhn_cell_local_edges", edge_count.to_string()),
                (
                    "kuhn_edge_census",
                    format!("{}|{}|{}", census[0], census[1], census[2]),
                ),
                ("kuhn_max_components", max_components.to_string()),
                ("kuhn_max_triangles", max_triangles.to_string()),
                ("kuhn_max_cut_edges", max_cut_edges.to_string()),
                ("control_first_violation", control_first_violation.clone()),
                ("hull_facets", exhibited_facets.to_string()),
                ("regularity_inequalities", system.len().to_string()),
                ("distinct_inequalities", distinct.len().to_string()),
                ("submodular_diamonds", diamonds.len().to_string()),
                (
                    "diamonds_are_regularity_rows",
                    diamonds_are_rows.to_string(),
                ),
                ("farkas_certificates", certificates.to_string()),
                ("farkas_total_weight", certificate_weight.to_string()),
                ("fm_feasible", elimination.feasible.to_string()),
                ("fm_max_stage_rows", elimination.max_rows.to_string()),
                ("fm_stage_rows", stage_sizes.clone()),
                ("fm_max_coefficient", fm_max_coefficient.to_string()),
                (
                    "solver_refuses_infeasible",
                    solver_refuses_infeasible.to_string(),
                ),
                ("liftings_tested", candidate_count.to_string()),
                ("liftings_accepted", accepted.to_string()),
                ("liftings_rejected", rejected.len().to_string()),
                ("rejected_liftings", rejected_names.clone()),
                ("lovasz_probes", lovasz_probes.to_string()),
                ("lovasz_mismatches", lovasz_mismatches.to_string()),
                (
                    "weak_submodular_iff_convex",
                    format!("{weak_iff_convex}of{candidate_count}"),
                ),
                (
                    "strict_submodular_iff_kuhn",
                    format!("{strict_iff_kuhn}of{candidate_count}"),
                ),
                ("fixtures_fired", fixtures_fired.to_string()),
                ("wall_ms", format!("{:.3}", sweep.wall_ms)),
            ]);
        }

        let global_c2 = sweeps.iter().all(|s| s.disagreements() == 0);
        println!(
            "C1 {c1}   C2 {global_c2} ({} of {} rows clean)   C3 {c3}\n\
             corpus: {} active, {corpus_certified} certified, {corpus_uncertified} \
             uncertified, {} degenerate, {} certified-with-two-components",
            sweeps.iter().filter(|s| s.disagreements() == 0).count(),
            sweeps.len(),
            sweeps.iter().map(|s| s.active).sum::<u64>(),
            sweeps.iter().map(|s| s.exact_zero).sum::<u64>(),
            sweeps.iter().map(|s| s.certified_multi).sum::<u64>(),
        );
    });
}
