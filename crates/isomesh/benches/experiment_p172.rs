//! **P-172 — a null registered on purpose: Conley index dies on plateaus, and
//! the census says how big they are.**
//!
//! Ticket: R-172. Pre-registered before this harness existed.
//!
//! ```bash
//! cargo bench --bench experiment_p172
//! ```
//!
//! Writes `docs/experiments/p-172.csv`.
//!
//! # What was missing
//!
//! Conley index theory is the standard answer to *"Morse theory needs
//! non-degenerate critical points and my function has not got them"*: it asks
//! only for an **isolating neighbourhood** `N` — a compact set whose maximal
//! invariant set lies in `int N` — and it is computer-assisted by construction.
//! The corpus holds Dey, Haas & Lipiński on Conley–Morse persistence barcodes.
//! Nothing in this repository had ever asked whether the crate's own fields
//! satisfy that one hypothesis.
//!
//! The neighbouring measurements say why the question is live:
//!
//! - `M-352` / `P-53` is the `=`-corner case: this crate's fields tie **exactly**
//!   rather than approximately, because `box_exact`'s zero set is a coordinate
//!   plane and its samples land on it.
//! - `P-58` / `R-056` (`experiment_p58.rs:56-63`) had to solve that for discrete
//!   Morse theory. Robins, Wood & Sheppard's `ProcessLowerStars`
//!   (`10.1109/tpami.2011.95`, Algorithm 1, §3.1) *"wants distinct voxel values,
//!   and this crate's reference fields tie exactly"*, and the paper's own
//!   tie-break — Eq. (8)'s global ramp `η(i + Ij + IJk)/(3IJK)` — is unusable
//!   here because it depends on the whole image dimensions and is therefore
//!   chunk-dependent and hash-breaking. `p-58.csv` records the size of the
//!   problem: `box_exact` at `17³` has **36 distinct values among 4 913
//!   voxels**, `sphere` 116, `torus` 325.
//! - `docs/research/2026-08-12-axes-and-vocabulary.md:117` is the sentence the
//!   whole line descends from: *"Morse theory — the parent theory. Ambiguity is
//!   really 'a critical point falls inside a cell'."*
//!
//! So the crate has a measured tie census (`P-58`) and no plateau census, and no
//! number anywhere saying which of its cells are inside Conley's hypothesis.
//! `crates/isomesh/src/**` is read-only for the phase, so both criteria below are
//! bench-local and driven through the public API.
//!
//! # The criterion for a plateau, stated exactly
//!
//! The unit of measurement is the **grid cell**, and the object measured is the
//! **trilinear interpolant of its eight corner samples** — the piecewise model
//! every extractor in this crate actually consumes
//! (`marching_cubes/trilinear.rs`), not the analytic field. That choice is
//! load-bearing and is the reason the numbers below are not zero; see *"What the
//! registration got slightly wrong"*.
//!
//! > A cell is a **plateau cell** when at least one of its six faces is
//! > numerically constant: the four corner samples on that face agree to within
//! > `FLAT_REL · cell_size`, with `FLAT_REL = 1e-9`.
//!
//! One criterion, not two. The tolerance is *relative to the cell size* because
//! for an exact signed distance field a corner-to-corner change of one
//! `cell_size` is unit slope, so `FLAT_REL` is a dimensionless slope floor and
//! the criterion means the same thing at `17³` and at `65³`. `plateau_cells` is
//! the count and `plateau_fraction` is that over `(n−1)³`.
//!
//! This is the registration's own mechanism transcribed rather than
//! reinterpreted: *"`box_exact` has exactly-equal samples across whole faces"*.
//! `plateau_cells_exact` repeats the scan with the tolerance replaced by
//! `f64::total_cmp` equality, so a reader can see whether `FLAT_REL` is doing any
//! work at all.
//!
//! # The criterion for an isolating neighbourhood, stated exactly
//!
//! For the negative gradient flow `ẋ = −∇f` the invariant sets are the critical
//! points and the orbits connecting them, so `Inv(Q) ⊂ int Q` fails for a closed
//! cell `Q` exactly when `∇f` vanishes somewhere on `∂Q`.
//!
//! > A cell **admits an isolating neighbourhood** when the trilinear interpolant
//! > of its eight corners has **no critical point on `∂Q`**: on none of the six
//! > faces do all three partial derivatives vanish simultaneously.
//!
//! That is exactly computable from the eight corners, and the algebra is short
//! enough to state. On the face `(axis, side)` take the crate's own tangential
//! axes `u = (axis+1) % 3` and `v = (axis+2) % 3`
//! (`cube.rs:88-101`, `face_corners`) and let `q00, q10, q01, q11` be the four
//! corner samples in `(s, t)` order. The restriction is bilinear,
//! `b(s,t) = q00 + B·s + C·t + D·s·t` with `B = q10−q00`, `C = q01−q00`,
//! `D = q00−q10−q01+q11`, and the two tangential partials are `b_s = B + D·t`,
//! `b_t = C + D·s`. The normal partial `∂f/∂u_axis` is independent of `u_axis`
//! and equals the opposite face's bilinear minus this one's. Three cases, and
//! they are exhaustive:
//!
//! 1. **`D` non-flat.** `b_s = b_t = 0` has the unique solution
//!    `t* = −B/D`, `s* = −C/D`. If it lies in `[0,1]²` and the normal partial is
//!    flat there, `∂Q` carries a critical point.
//! 2. **`B`, `C` and `D` all flat** — the constant face, i.e. the plateau. Every
//!    point of the face is tangentially critical, so `∂Q` carries a critical
//!    point unless the normal partial is strictly one-signed over the whole
//!    face. A bilinear on `[0,1]²` has a saddle for its only interior critical
//!    point (Hessian `[[0,D],[D,0]]`, `det = −D² ≤ 0`) and is linear on each
//!    boundary edge, so its **range over the square is the range over its four
//!    corners**. The test is therefore exact and needs no sampling: the four
//!    normal differences must all be `> tol` or all be `< −tol`.
//! 3. **`D` flat and at least one of `B`, `C` not.** Then one tangential partial
//!    is a non-zero constant and no point of the face is tangentially critical.
//!
//! `conley_applicable_fraction` is the fraction of cells admitting one, and the
//! registered boolean `isolating_neighbourhood_exists` is the per-row statement
//! that **every** cell of that grid admits one — `isolating_cells == cells`,
//! compared as integers so no float equality is involved.
//!
//! Case 2 is where the plateau bites, and it bites for a reason worth spelling
//! out because the naive expectation is wrong. A constant face is *not*
//! automatically an obstruction: `box_exact` in its face slab has `f = x − 1`,
//! whose `x`-constant faces are constant with a strictly positive normal
//! derivative, and Conley is perfectly happy there. The obstruction is the
//! **medial axis**, where the argmax in
//! `f = max(|x|,|y|,|z|) − 1` (`fields/mod.rs:438-443`) switches: on the cell
//! `[0.5,0.75] × [0.5,0.75] × [0,0.25]` the face `x = 0.75` is constant at
//! `−0.25` and its four normal differences are `(−0.25, 0, −0.25, 0)` — two
//! exact zeros, so `∇f = 0` at two boundary points and there is no isolating
//! neighbourhood. `boundary_critical_without_plateau` is the column that decides
//! whether the plateau is the mechanism or merely a correlate.
//!
//! # The criterion for discrete Morse theory is *not* invented here
//!
//! It is `P-58`'s, cited rather than restated, because a second criterion for a
//! concept the repository has already registered would be two paths to one
//! answer:
//!
//! > A cell is **discrete-Morse-applicable** when its eight corner samples are
//! > pairwise distinct under `f64::total_cmp`.
//!
//! `experiment_p58.rs:58-63` states the requirement (*"Algorithm 1 wants
//! distinct voxel values"*) and `p-58.csv`'s `distinct_values` / `tied_voxels`
//! columns are its measurement, taken on the grids `[17, 33, 65]`. This harness
//! uses **the same three resolutions for that reason**: the
//! `vs_discrete_morse_applicable` column here sits on the same grids as
//! `p-58.csv`'s tie census and the two can be read against each other without a
//! rerun.
//!
//! Distinctness is tested **exactly**, with no tolerance, while the plateau and
//! isolation tests are tolerant. That asymmetry is deliberate and is the honest
//! shape of the two hypotheses: Algorithm 1's requirement is about **labels**,
//! and a total order either exists or does not; Conley's is about the
//! **geometry** of an interpolant, and every quantity in it is a computed
//! difference.
//!
//! # What the registration got slightly wrong, recorded rather than smoothed
//!
//! The hypothesis says `box_exact` has *"`∇f = 0` on a set of positive
//! measure"*. Measured, that is true in **two** dimensions and false in
//! **three**, and both halves are columns:
//!
//! - `constant_faces` counts numerically constant cell faces, and it is large.
//! - `constant_cells` counts cells whose **eight** corners agree — a genuine
//!   three-dimensional plateau, on which the interpolant is constant and every
//!   point of the closed cell is an equilibrium. It is predicted **zero on every
//!   field at every resolution**, and the arithmetic is short: `box_exact`'s
//!   analytic gradient is *unit length everywhere* (`fields/mod.rs:446-475`
//!   returns a normalised vector in the exterior branch and a signed basis
//!   vector `g[j] = d[j].signum()` in the interior branch), so the field is
//!   nowhere locally constant in 3D and no cell can have eight equal corners.
//!
//! So the mechanism is real and its dimension is one lower than claimed. A
//! plateau-shaped obstruction that lives on `∂Q` rather than in `int Q` is if
//! anything *worse* for Conley, because `∂Q` is exactly where the theory forbids
//! invariant behaviour — but it is not what the sentence said, and saying so is
//! cheaper than a reader rediscovering it.
//!
//! This is also why the criterion is applied to the interpolant of the sampled
//! grid rather than to the analytic field. Applied to the analytic gradient the
//! whole row would be vacuous: `‖∇f‖ = 1` identically on `box_exact`, so a
//! gradient-norm floor would report zero plateau cells and fire this row's own
//! vacuity control on a fixture that manifestly has the defect.
//!
//! # Arms
//!
//! | arm | what it varies | is_control |
//! |---|---|---|
//! | `box_exact` × {17, 33, 65} | the registered positive control — exactly-equal samples across whole faces | **yes**, the registration's own |
//! | `thin_plate` × {17, 33, 65} | the same box formula at sub-cell thickness (`fields/mod.rs:573`) | no |
//! | `csg_difference` × {17, 33, 65} | `BoxExact ∖ Sphere` — a box plateau cut by a smooth surface | no |
//! | `sphere` × {17, 33, 65} | smooth, radially symmetric, massively tied *in value* and not *in face* | **yes**, negative control |
//! | `torus` × {17, 33, 65} | smooth, genus 1 | **yes**, negative control |
//! | `gyroid` × {17, 33, 65} | smooth and non-convex, `Lipschitz` bound | **yes**, negative control |
//! | `fbm_terrain` × {17, 33, 65} | not closed in domain, `Unbounded` | no |
//! | `noise_cavity` × {17, 33, 65} | value noise ∩ sphere, `Unbounded` | no |
//!
//! Twenty-four rows. `17, 33, 65` are `P-58`'s `RESOLUTIONS`, for the
//! comparability reason above, and no clause reads a wall clock.
//!
//! # SHARE, recomputed before the numbers
//!
//! The registration's SHARE line is *"none — this closes a candidate"*, and that
//! discharges without arithmetic: nothing here proposes to move a stage, so
//! there is no share to price and no `1/(1 − share/factor)` ceiling to check.
//! What the row produces is a **refusal**, and the refusal is quantified: after
//! it, a future ticket proposing Conley index theory for ambiguity resolution
//! has to explain the `box_exact` number rather than argue past it.
//!
//! `ns_per_cell` is on every row as provenance and **gates nothing** — `M-280`
//! and `M-348` are the incidents; a census is a count.
//!
//! # Vacuity controls
//!
//! Every one runs before the first `run.record` and every panic starts `VOID: `.
//!
//! - **The registered control.** Every `box_exact` row must report
//!   `plateau_cells > 0`, or the obstruction this whole row is built on does not
//!   exist in the fixture. Column: `plateau_cells`.
//! - **The negative control on the plateau criterion.** At least one row must
//!   report `plateau_cells == 0`. Without it a criterion that fires on every cell
//!   would report `0.5` on `box_exact` and mean nothing — `M-44`'s rule read
//!   backwards: a non-zero that could not have been zero is not a measurement
//!   either. Column: `plateau_cells`.
//! - **The isolation test must be able to fail.** At least one row must report
//!   `conley_applicable_fraction < 1`, or `∂Q`-criticality is being asserted
//!   rather than measured. Column: `conley_applicable_fraction`.
//! - **The isolation test must be able to pass.** At least one row must report
//!   `isolating_neighbourhood_exists = true`, or the fraction is a constant and
//!   C2's ranking is computed against one. Column:
//!   `isolating_neighbourhood_exists`.
//! - **The discrete-Morse column must be non-trivial.** Every row must report
//!   `vs_discrete_morse_applicable > 0`, and at least one row must report it
//!   `< 1`, or C2 compares Conley against a constant. Column:
//!   `vs_discrete_morse_applicable`.
//!
//! # Which column decides which clause
//!
//! - **C1** — *"the plateau fraction is measured per field, quantifying the
//!   obstruction rather than asserting it"*, falsified by *"negligible plateaus
//!   on every field"*. That is a clause about the sweep, not about a row, so
//!   `c1_holds` carries the **same global verdict on every row**: it is `true`
//!   when some row's `plateau_fraction` reaches `NEGLIGIBLE_PERCENT`, which is
//!   **1%** of its cells. Scored as the integer
//!   comparison `plateau_cells · 100 ≥ cells`, so no threshold is decided by a
//!   float. The per-row form is the extra column `row_plateau_is_negligible`.
//! - **C2** — the ranking. `c2_holds` is per-row and is
//!   `conley_applicable_fraction ≤ vs_discrete_morse_applicable`, scored as
//!   `isolating_cells ≤ discrete_morse_applicable_cells` because the two
//!   fractions share the denominator `cells`, so the comparison is exact
//!   integers. Its falsifier is explicit: *"C2 by Conley applying more widely
//!   than discrete Morse, which would reverse the ranking."*
//!
//!   The comparison is `≤`, so a field on which **both** fractions saturate at
//!   `1` satisfies C2 at equality rather than by a strict ranking. That is not a
//!   hypothetical: `fbm_terrain` is generic — a four-octave value-noise
//!   heightfield has no exactly-tied corner pair — so its discrete-Morse
//!   fraction is `1` by construction, and if its Conley fraction is `1` too then
//!   C2 holds there vacuously. `conley_minus_morse_fraction` is the column that
//!   tells the two apart: a strict ranking is a **negative** value, and a
//!   saturated tie is exactly `0`.
//!
//! # The one methodological fork in C2, decided here and stated
//!
//! `P-58` did not merely record discrete Morse theory's tie problem, it
//! **repaired** it: a chunk-local exact total order on `(value, linear_index)`
//! with `f64::total_cmp` on values, which perturbs no sample, needs no `η`, and
//! whose census `P-58` then measured to be invariant under a second tie-break
//! (`census_matches_reverse_order`). So there are two defensible denominators for
//! *"the fraction where discrete Morse theory applies"* and they give opposite
//! verdicts on C2.
//!
//! The registered column takes the **raw requirement**, because the clause says
//! *"the fraction where discrete Morse theory — already registered — applies"*
//! and Algorithm 1's requirement is distinct values; the repair is `P-58`'s
//! finding, not the paper's hypothesis. Ranking Conley's raw requirement against
//! discrete Morse's *repaired* one would compare a method against a
//! method-plus-a-repository-invention.
//!
//! The other reading is not dropped, it is a column:
//! `conley_minus_repaired_morse_fraction` is
//! `conley_applicable_fraction − 1`, the signed gap against the repaired
//! denominator, which is `1` by construction because a total order always
//! exists. Both are on every row, so the
//! finding can state the whole shape — that both requirements fail on this
//! crate's data, that discrete Morse's failure is repairable by relabelling and
//! Conley's is not, because one is a statement about labels and the other about
//! the geometry of an interpolant.

#![allow(
    clippy::cast_possible_truncation,
    clippy::cast_precision_loss,
    clippy::too_many_lines
)]

mod common;

use std::time::Instant;

use isomesh::Sdf;
use isomesh::marching_cubes::table::{EDGE_AXIS, EDGE_CORNERS, face_corners};

type Scalar = f64;

/// Samples per axis. `P-58`'s `RESOLUTIONS`, so `vs_discrete_morse_applicable`
/// sits on the same grids as `p-58.csv`'s `distinct_values`.
const RESOLUTIONS: [u32; 3] = [17, 33, 65];

/// Dimensionless slope floor: a difference of samples is *flat* when it is below
/// `FLAT_REL · cell_size`.
///
/// For an exact signed distance field a corner-to-corner change of one
/// `cell_size` is unit slope, so this is a slope threshold and means the same
/// thing at every resolution.
const FLAT_REL: Scalar = 1e-9;

/// C1's bar for "not negligible", as a whole percent of a grid's cells.
///
/// Applied as the integer comparison `plateau_cells * 100 >= cells`, so the
/// threshold is not decided by a float.
const NEGLIGIBLE_PERCENT: u64 = 1;

/// Corner `c`'s offset in this crate's cube numbering: bit `k` is axis `k`.
///
/// `crate::cube::corner_offset` is private and not re-exported
/// (`marching_cubes/table.rs:88-91`), so the three lines live here and are
/// asserted against the public `EDGE_CORNERS` / `EDGE_AXIS` in
/// [`assert_conventions`].
const fn corner_offset(c: usize) -> [usize; 3] {
    [c & 1, (c >> 1) & 1, (c >> 2) & 1]
}

/// What one cell's eight corner samples say about the two methods.
#[derive(Clone, Copy)]
struct Verdict {
    /// How many of the six faces are numerically constant.
    constant_faces: u32,
    /// At least one face is constant under `f64::total_cmp` equality, with no
    /// tolerance at all.
    plateau_exact: bool,
    /// The trilinear interpolant has no critical point on the cell boundary.
    isolating: bool,
}

impl Verdict {
    /// At least one of the six faces is numerically constant.
    const fn plateau(self) -> bool {
        self.constant_faces > 0
    }

    /// All six faces are constant, so the interpolant is constant on the closed
    /// cell — the three-dimensional plateau the registration claims.
    const fn constant_cell(self) -> bool {
        self.constant_faces == 6
    }
}

/// The two criteria, evaluated on one cell's eight corner samples.
///
/// `tol` is the absolute flatness threshold, `FLAT_REL * cell_size`. See the
/// module header for the three exhaustive cases and why the constant-face case
/// needs only the four corner values of the normal partial.
fn verdict(cv: &[Scalar; 8], tol: Scalar) -> Verdict {
    let mut constant_faces = 0u32;
    let mut plateau_exact = false;
    let mut isolating = true;

    for axis in 0..3usize {
        // The crate's own tangential axes for this face (`cube.rs:89-90`).
        let u = (axis + 1) % 3;
        let v = (axis + 2) % 3;
        for side in 0..2usize {
            let at = |s: usize, t: usize| cv[(side << axis) | (s << u) | (t << v)];
            let opp = |s: usize, t: usize| cv[((1 - side) << axis) | (s << u) | (t << v)];

            let q00 = at(0, 0);
            let q10 = at(1, 0);
            let q01 = at(0, 1);
            let q11 = at(1, 1);
            let b = q10 - q00;
            let c = q01 - q00;
            let d = q00 - q10 - q01 + q11;

            if b.abs() <= tol && c.abs() <= tol && d.abs() <= tol {
                // Case 2: the constant face. Every point is tangentially
                // critical, so the boundary is clean only if the normal partial
                // is strictly one-signed over the whole face — and a bilinear's
                // range over the square is its range over the four corners.
                constant_faces += 1;
                if q10.total_cmp(&q00).is_eq()
                    && q01.total_cmp(&q00).is_eq()
                    && q11.total_cmp(&q00).is_eq()
                {
                    plateau_exact = true;
                }
                let normal = [
                    opp(0, 0) - q00,
                    opp(1, 0) - q10,
                    opp(0, 1) - q01,
                    opp(1, 1) - q11,
                ];
                let all_positive = normal.iter().all(|&x| x > tol);
                let all_negative = normal.iter().all(|&x| x < -tol);
                if !(all_positive || all_negative) {
                    isolating = false;
                }
            } else if d.abs() > tol {
                // Case 1: `b_s = B + D t` and `b_t = C + D s` have one common
                // zero. Only a zero inside the closed face counts.
                let s = -c / d;
                let t = -b / d;
                if (0.0..=1.0).contains(&s) && (0.0..=1.0).contains(&t) {
                    let bilinear = |a00: Scalar, a10: Scalar, a01: Scalar, a11: Scalar| {
                        a00 + (a10 - a00) * s + (a01 - a00) * t + (a00 - a10 - a01 + a11) * s * t
                    };
                    let near = bilinear(q00, q10, q01, q11);
                    let far = bilinear(opp(0, 0), opp(1, 0), opp(0, 1), opp(1, 1));
                    if (far - near).abs() <= tol {
                        isolating = false;
                    }
                }
            }
            // Case 3: `D` flat and one of `B`, `C` not — one tangential partial
            // is a non-zero constant, so no point of the face is critical and
            // there is nothing to record.
        }
    }

    Verdict {
        constant_faces,
        plateau_exact,
        isolating,
    }
}

/// `P-58`'s discrete-Morse requirement, restricted to one cell.
///
/// `experiment_p58.rs:58-63`: Algorithm 1 of `10.1109/tpami.2011.95` wants
/// distinct values at every vertex. Exact, via `total_cmp`, with no tolerance —
/// the requirement is about labels, not about geometry.
fn corners_distinct(cv: &[Scalar; 8]) -> bool {
    let mut sorted = *cv;
    sorted.sort_unstable_by(|a, b| a.total_cmp(b));
    sorted.windows(2).all(|w| w[0].total_cmp(&w[1]).is_ne())
}

/// One `(field, resolution)` arm's census.
struct Row {
    field: &'static str,
    resolution: u32,
    cells: u64,
    plateau_cells: u64,
    plateau_cells_exact: u64,
    constant_cells: u64,
    constant_faces: u64,
    isolating_cells: u64,
    plateau_and_isolating: u64,
    boundary_critical_without_plateau: u64,
    morse_cells: u64,
    tolerance: Scalar,
    ns_per_cell: Scalar,
}

impl Row {
    fn plateau_fraction(&self) -> Scalar {
        self.plateau_cells as Scalar / self.cells as Scalar
    }

    fn conley_fraction(&self) -> Scalar {
        self.isolating_cells as Scalar / self.cells as Scalar
    }

    fn morse_fraction(&self) -> Scalar {
        self.morse_cells as Scalar / self.cells as Scalar
    }

    /// The registered boolean: every cell of this grid admits an isolating
    /// neighbourhood. Integer comparison, so no float equality is involved.
    fn isolating_everywhere(&self) -> bool {
        self.isolating_cells == self.cells
    }

    /// C2's ranking, exact: the two fractions share the denominator `cells`.
    fn c2(&self) -> bool {
        self.isolating_cells <= self.morse_cells
    }

    /// C1's bar, on this row alone. The clause itself is global.
    fn row_is_negligible(&self) -> bool {
        self.plateau_cells * 100 < self.cells * NEGLIGIBLE_PERCENT
    }
}

/// Census one grid.
fn census(field: &'static str, n: u32, values: &[Scalar], cell_size: Scalar) -> Row {
    let nu = n as usize;
    let plane = nu * nu;
    let tol = FLAT_REL * cell_size;

    let mut row = Row {
        field,
        resolution: n,
        cells: (nu - 1).pow(3) as u64,
        plateau_cells: 0,
        plateau_cells_exact: 0,
        constant_cells: 0,
        constant_faces: 0,
        isolating_cells: 0,
        plateau_and_isolating: 0,
        boundary_critical_without_plateau: 0,
        morse_cells: 0,
        tolerance: tol,
        ns_per_cell: 0.0,
    };

    let started = Instant::now();
    for z in 0..nu - 1 {
        for y in 0..nu - 1 {
            for x in 0..nu - 1 {
                let base = z * plane + y * nu + x;
                let mut cv = [0.0; 8];
                for (c, slot) in cv.iter_mut().enumerate() {
                    let o = corner_offset(c);
                    *slot = values[base + o[2] * plane + o[1] * nu + o[0]];
                }

                let v = verdict(&cv, tol);
                row.constant_faces += u64::from(v.constant_faces);
                if v.constant_cell() {
                    row.constant_cells += 1;
                }
                if v.plateau_exact {
                    row.plateau_cells_exact += 1;
                }
                if v.plateau() {
                    row.plateau_cells += 1;
                    if v.isolating {
                        row.plateau_and_isolating += 1;
                    }
                } else if !v.isolating {
                    row.boundary_critical_without_plateau += 1;
                }
                if v.isolating {
                    row.isolating_cells += 1;
                }
                if corners_distinct(&cv) {
                    row.morse_cells += 1;
                }
            }
        }
    }
    row.ns_per_cell = started.elapsed().as_secs_f64() * 1e9 / row.cells as Scalar;

    // Transcription checks, not vacuity controls: each one fires only if the
    // census contradicts its own arithmetic.
    assert!(
        row.plateau_cells <= row.cells
            && row.isolating_cells <= row.cells
            && row.morse_cells <= row.cells,
        "{field} at {n}³: a per-cell count exceeds the cell count"
    );
    assert!(
        row.plateau_cells_exact <= row.plateau_cells,
        "{field} at {n}³: an exactly-constant face is not tolerantly constant"
    );
    assert!(
        row.constant_cells <= row.plateau_cells,
        "{field} at {n}³: a constant cell is not a plateau cell"
    );
    assert!(
        row.plateau_and_isolating <= row.plateau_cells,
        "{field} at {n}³: more isolating plateaus than plateaus"
    );
    assert!(
        row.boundary_critical_without_plateau <= row.cells - row.isolating_cells,
        "{field} at {n}³: more plateau-free boundary-critical cells than boundary-critical cells"
    );
    assert!(
        row.constant_faces <= 6 * row.cells,
        "{field} at {n}³: {} constant faces against {} cell faces",
        row.constant_faces,
        6 * row.cells
    );

    row
}

/// The corner numbering and the face parameterisation, checked against the
/// crate's own public tables rather than assumed.
fn assert_conventions() {
    // Bit `k` is axis `k`, from `EDGE_CORNERS` / `EDGE_AXIS` (`P-58`'s check).
    for (e, corners) in EDGE_CORNERS.iter().enumerate() {
        assert_eq!(
            corners[0] ^ corners[1],
            1u8 << EDGE_AXIS[e],
            "edge {e} does not join corners differing in the bit of its own axis"
        );
    }
    for c in 0..8usize {
        for (k, &v) in corner_offset(c).iter().enumerate() {
            assert_eq!(v, (c >> k) & 1, "corner bit {k} is not axis {k}");
        }
    }
    // The `(side << axis) | (s << u) | (t << v)` face parameterisation [`verdict`]
    // uses must enumerate exactly the crate's `face_corners(axis, side)`. Getting
    // this wrong silently measures a different four samples per face.
    for axis in 0..3usize {
        let u = (axis + 1) % 3;
        let v = (axis + 2) % 3;
        for side in 0..2usize {
            let mut mine = [0u8; 4];
            for (i, slot) in mine.iter_mut().enumerate() {
                let s = i & 1;
                let t = (i >> 1) & 1;
                *slot = ((side << axis) | (s << u) | (t << v)) as u8;
            }
            mine.sort_unstable();
            let mut theirs = face_corners(axis, side as u8);
            theirs.sort_unstable();
            assert_eq!(
                mine, theirs,
                "the face ({axis}, {side}) parameterisation does not enumerate face_corners"
            );
        }
    }
}

fn main() {
    if !std::env::args().any(|arg| arg == "--bench") {
        return;
    }
    assert_conventions();

    let prereg = isomesh::experiment!("P-172");

    common::experiment::run(prereg, |run| {
        let mut rows: Vec<Row> = Vec::new();

        isomesh::for_each_reference_field!(f64, |name, field| {
            // Inlined block, so no `return` in here (M-199 / M-253).
            for &n in &RESOLUTIONS {
                let (_shape, origin, h) = common::grid::<Scalar, _>(&field, n);
                let nu = n as usize;
                let mut values = Vec::with_capacity(nu * nu * nu);
                for z in 0..nu {
                    for y in 0..nu {
                        for x in 0..nu {
                            values.push(field.sample([
                                origin[0] + h * x as Scalar,
                                origin[1] + h * y as Scalar,
                                origin[2] + h * z as Scalar,
                            ]));
                        }
                    }
                }
                rows.push(census(name, n, &values, h));
            }
        });

        // ── vacuity controls, all before the first record ──
        for row in rows.iter().filter(|r| r.field == "box_exact") {
            assert!(
                row.plateau_cells > 0,
                "VOID: box_exact at {}³ reports no plateau cell, so the obstruction this row \
                 is built on does not exist in the fixture and both clauses are unmeasured",
                row.resolution
            );
        }
        assert!(
            rows.iter().any(|r| r.plateau_cells == 0),
            "VOID: every row reports a plateau cell, so the criterion fires everywhere and \
             box_exact's fraction is a non-zero that could not have been zero (M-44 read \
             backwards)"
        );
        assert!(
            rows.iter().any(|r| r.isolating_cells < r.cells),
            "VOID: no row reports a boundary-critical cell, so `∂Q`-criticality is being \
             asserted rather than measured and C2's Conley fraction is the constant 1"
        );
        assert!(
            rows.iter().any(Row::isolating_everywhere),
            "VOID: no row reports an isolating neighbourhood on every cell, so \
             conley_applicable_fraction never separates the fields and C2 ranks against a \
             constant"
        );
        assert!(
            rows.iter().all(|r| r.morse_cells > 0),
            "VOID: a row reports no discrete-Morse-applicable cell at all, so P-58's \
             requirement is the constant false and C2 compares Conley against nothing"
        );
        assert!(
            rows.iter().any(|r| r.morse_cells < r.cells),
            "VOID: every cell of every row satisfies P-58's distinctness requirement, so the \
             tie census p-58.csv recorded has vanished and C2 ranks against the constant 1"
        );

        // C1 is a clause about the sweep, so one verdict goes on every row.
        let c1 = rows
            .iter()
            .any(|r| r.plateau_cells * 100 >= r.cells * NEGLIGIBLE_PERCENT);
        let c2_global = rows.iter().all(Row::c2);

        println!(
            "{:<15} {:>4} {:>8} {:>8} {:>9} {:>8} {:>9} {:>8} {:>9} {:>7} {:>6} {:>6}",
            "field",
            "n",
            "cells",
            "plateau",
            "plateau_f",
            "isolate",
            "conley_f",
            "morse",
            "morse_f",
            "crit-pl",
            "const",
            "c2"
        );

        for row in &rows {
            println!(
                "{:<15} {:>4} {:>8} {:>8} {:>9.6} {:>8} {:>9.6} {:>8} {:>9.6} {:>7} {:>6} {:>6}",
                row.field,
                row.resolution,
                row.cells,
                row.plateau_cells,
                row.plateau_fraction(),
                row.isolating_cells,
                row.conley_fraction(),
                row.morse_cells,
                row.morse_fraction(),
                row.boundary_critical_without_plateau,
                row.constant_cells,
                row.c2(),
            );

            run.record(&[
                ("field", row.field.to_string()),
                ("resolution", row.resolution.to_string()),
                ("plateau_cells", row.plateau_cells.to_string()),
                ("plateau_fraction", format!("{:.6}", row.plateau_fraction())),
                (
                    "isolating_neighbourhood_exists",
                    row.isolating_everywhere().to_string(),
                ),
                (
                    "conley_applicable_fraction",
                    format!("{:.6}", row.conley_fraction()),
                ),
                (
                    "vs_discrete_morse_applicable",
                    format!("{:.6}", row.morse_fraction()),
                ),
                ("c1_holds", c1.to_string()),
                ("c2_holds", row.c2().to_string()),
                // ── extras (M-273) ──
                ("cells", row.cells.to_string()),
                ("isolating_cells", row.isolating_cells.to_string()),
                (
                    "boundary_critical_cells",
                    (row.cells - row.isolating_cells).to_string(),
                ),
                (
                    "boundary_critical_without_plateau",
                    row.boundary_critical_without_plateau.to_string(),
                ),
                (
                    "plateau_and_isolating",
                    row.plateau_and_isolating.to_string(),
                ),
                ("plateau_cells_exact", row.plateau_cells_exact.to_string()),
                ("constant_cells", row.constant_cells.to_string()),
                ("constant_faces", row.constant_faces.to_string()),
                (
                    "discrete_morse_applicable_cells",
                    row.morse_cells.to_string(),
                ),
                (
                    "tied_corner_cells",
                    (row.cells - row.morse_cells).to_string(),
                ),
                (
                    "conley_minus_morse_fraction",
                    format!("{:.6}", row.conley_fraction() - row.morse_fraction()),
                ),
                (
                    "conley_minus_repaired_morse_fraction",
                    format!("{:.6}", row.conley_fraction() - 1.0),
                ),
                (
                    "row_plateau_is_negligible",
                    row.row_is_negligible().to_string(),
                ),
                ("c2_global_holds", c2_global.to_string()),
                ("flat_relative_tolerance", format!("{FLAT_REL:.3e}")),
                ("flat_absolute_tolerance", format!("{:.3e}", row.tolerance)),
                ("negligible_percent", NEGLIGIBLE_PERCENT.to_string()),
                ("ns_per_cell", format!("{:.3}", row.ns_per_cell)),
                (
                    "criterion_source",
                    String::from(
                        "plateau = a numerically constant cell face; isolating = no critical \
                         point of the trilinear on dQ; discrete Morse = P-58/R-056 distinct \
                         vertex values (10.1109/tpami.2011.95 Alg.1 s3.1)",
                    ),
                ),
            ]);
        }
    });
}
