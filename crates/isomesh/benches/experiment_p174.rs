//! **P-174 — a null registered on purpose: varifolds against the normal cycles already benched.**
//!
//! Ticket: R-174. Pre-registered before this harness existed.
//!
//! ```bash
//! cargo bench --bench experiment_p174
//! ```
//!
//! Writes `docs/experiments/p-174.csv`.
//!
//! # What was missing
//!
//! The crate's curvature line is normal cycles. `benches/experiment_p42.rs:13-14`
//! and `benches/experiment_p45.rs:26-27` both cite Sun & Morvan, *Curvature
//! measures, normal cycles and asymptotic cones*, `10.5802/acirm.50`, and both
//! compute the same two **totals**: the Gaussian measure `Sum_v (2pi - alpha_v)`
//! and the mean measure `Sum_e l(e) beta(e)` of Theorem 5 (2)(b).
//!
//! Two prior findings say exactly why this row exists.
//!
//! - **`✗30` / `M-340`** (`FINDINGS.md:8009`, from `P-42`) — *"with `B` the whole
//!   closed surface, the registered residual is discrete Gauss-Bonnet in an
//!   approximation check's costume, and what it measures is one f64 epsilon per
//!   vertex"*. The committed artefact agrees: `docs/experiments/p-42.csv` reports
//!   `gaussian_total = 12.566370614359` against `gaussian_expected =
//!   12.566370614359` with `residual = 2.8e-13` on the 33³ sphere, and
//!   `-0.00000000000` on the torus at every resolution. The whole-surface
//!   Gaussian measure is an **identity**. It carries no approximation content,
//!   so the only way to get information out of it is to **localise** it — and
//!   nothing in the crate or in `P-42`/`P-45` ever has.
//! - **`✗32` / `M-343`** (`FINDINGS.md:8262`, from `P-45`) — the two measures
//!   split: the Gaussian one is chunk-local and does not sum, the mean one sums
//!   and is not chunk-local. Also a statement about totals over regions, not
//!   about a curvature **at a point**.
//!
//! So the crate has never produced a pointwise curvature, and has never put its
//! estimator against a rival formalism. Discrete varifolds are that rival:
//! Buet, Leonardi & Masnou, *A varifold approach to surface approximation*,
//! `10.1007/s00205-017-1141-0` (acquired and converted this session), which
//! carries quantitative convergence bounds for exactly this construction.
//!
//! This is registered **expecting a null**, because a mature estimator is rarely
//! beaten by a newer formalism on the same data.
//!
//! # The two estimators, and the one thing that makes the comparison fair
//!
//! Both are localisations of a measure at **one shared radius `eps`**, over the
//! **same** neighbour search. That is not a convenience — it is what "at matched
//! cost" in `C1` means, and without it the comparison is rigged. A mollified
//! estimator with radius `eps` averages over `O((eps/h)²)` triangles; a
//! vertex-star estimator averages over six. Putting those two against each other
//! measures the radius, not the formalism.
//!
//! ## The mollifier
//!
//! One profile serves both arms, radially symmetric and compactly supported on
//! the ball of radius `eps`:
//!
//! ```text
//! rho_eps(z) = xi(|z| / eps),        xi(t)  = (1 - t²)³ on [0, 1], 0 outside
//! grad rho_eps(z) = xi'(|z|/eps) * z_hat / eps,   xi'(t) = -6 t (1 - t²)²
//! ```
//!
//! `xi(1) = xi'(1) = xi''(1) = 0`, so `rho_eps` is `C²` across the support
//! boundary — the regularity Buet-Leonardi-Masnou's first-variation estimate
//! asks for. `xi'(0) = 0`, so `grad rho_eps` is `0` at the centre and the
//! `z / |z|` factor never divides by zero. The normalising constant is omitted
//! on purpose: every quantity below is a **ratio** of a `rho`-weighted sum to a
//! `rho`-weighted sum, so the constant cancels exactly, and the surviving
//! `1 / eps` in `grad rho_eps` is the `1 / length` that curvature has.
//!
//! `eps = 0.8 * L0 * sqrt(h / L0)` with `L0 = 1` the unit length all three
//! fixtures live at, i.e. `eps = 0.8 * sqrt(h)`. That exponent is
//! Buet-Leonardi-Masnou's own balance: their error is `O(eps) + O(h / eps)`, so
//! `eps ~ sqrt(h)` sends both terms to zero together and is the only scaling
//! under which the varifold estimator converges at its best rate. Across the
//! ladder it gives `eps` = 0.28284, 0.23594, 0.20000, 0.16865, 0.14142 and
//! `eps / h` = 2.263, 2.713, 3.200, 3.795, 4.525 — `eps` falling, `eps / h`
//! rising, which is the regime the bounds are stated in. Both are asserted, not
//! assumed.
//!
//! ## Arm 1 — `varifold` (the challenger)
//!
//! A mesh is the varifold `V = Sum_T area(T) * delta_{(x_T, n_T)}`. The
//! registered form, implemented verbatim:
//!
//! ```text
//! Hvec_eps(x) = - ( Sum_T area(T) * P_{n_T} grad rho_eps(x - x_T) )
//!                 / ( Sum_T area(T) *          rho_eps(x - x_T) )
//! ```
//!
//! with `P_n = I - n (x) n`. Taking `X(y) = rho_eps(x - y) e_j` in
//! `delta V(X) = Sum_T area(T) div_{P_T} X(x_T)` and using
//! `delta V(X) = integral (k1 + k2) (X . n) d||V||` on a closed surface gives
//! `-num_j / den = (k1 + k2)(x) * n_j(x)`. So the registered expression is the
//! mean curvature **vector** `(k1 + k2) n_outward`, and the scalar is
//!
//! ```text
//! H_varifold(x) = 0.5 * ( Hvec_eps(x) . n_eps(x) )
//! ```
//!
//! Checked before this harness was written: on a unit-sphere triangulation the
//! expression returns `Hvec . n` = 1.9788 / 1.9837 / 1.9736 / 1.9701 at
//! `eps` = 0.15 / 0.2 / 0.3 / 0.4 against `k1 + k2 = 2`. The factor `0.5` is
//! that measurement, not a guess, and the calibration control re-asserts it here.
//!
//! `n_eps` is the varifold's **own** mollified normal, not an oracle:
//! `n_eps(x) = N(x) / |N(x)|` with `N(x) = Sum_T area(T) rho_eps(x - x_T) n_T`,
//! outward because the crate guarantees counter-clockwise-seen-from-outside
//! winding (`lib.rs:56-67`, the same guarantee `experiment_p42.rs:153-155`
//! rests on). `normal_alignment` records `n_eps . grad f / |grad f|` so a reader
//! can see the projection is doing no work.
//!
//! Gaussian curvature is not a first-variation quantity, so it comes from the
//! varifold's second-order data — the tangential derivative of that same
//! mollified normal, which is available **analytically in the same pass**:
//!
//! ```text
//! d_b N_a       = Sum_T area(T) (grad rho_eps)_b (n_T)_a
//! d_b (n_eps)_a = d_b N_a / |N| - N_a (N . d_b N) / |N|³
//! S_ij          = t_i . ( grad n_eps . t_j ),   {t_1, t_2} orthonormal ⊥ n_eps
//! K_varifold    = det S,        H_from_shape = 0.5 * tr S
//! ```
//!
//! There is no finite-difference step and therefore no second length scale, and
//! the whole arm is **one gather** over the triangle centroids: 19 scalar
//! accumulators. That is what keeps `cost_ratio` near one honestly instead of by
//! tuning. `det S` and `tr S` are invariant under the choice of orthonormal
//! tangent basis, so the (deterministic) basis choice cannot move a number.
//! `varifold_h_shape_error` is the independent cross-check: `0.5 tr S` and the
//! first-variation `H` are different formulas for the same quantity.
//!
//! ## Arm 2 — `normal_cycles` (the incumbent, and the control)
//!
//! Cohen-Steiner & Morvan's theorem is stated for a **ball**, not a vertex star,
//! so the faithful localisation of the shipped estimator at radius `eps` is its
//! own measure over that ball, divided by the mass of the same ball:
//!
//! ```text
//! den(x)            = Sum_T area(T) rho_eps(x - x_T)                -> integral rho dA
//! H_normalcycles(x) = 0.5 * Sum_e rho_eps(x - m_e) l(e) beta(e) / den(x)
//! K_normalcycles(x) = Sum_v rho_eps(x - p_v) (2pi - alpha_v)   / den(x)
//! ```
//!
//! `m_e` is the edge midpoint; `beta(e)` is `experiment_p42.rs:296`'s signed
//! dihedral `atan2((n1 x n2) . e_hat, n1 . n2)` with `e_hat` oriented as the
//! first face traverses the edge, re-derived here rather than re-tuned. The
//! `0.5` on the mean measure is `experiment_p42.rs:61-66`'s cylinder
//! calibration: a rounded `[-1,1]³` box gives `integral H da = 12 (pi/4) 2 = 6pi`
//! while `Sum l beta = 12 * 2 * (pi/2) = 12pi`, so `integral H da = ½ Sum l beta`.
//! Both ratios tend to `H(x)` and `K(x)` as `eps -> 0` with no further factor.
//!
//! `nc_star_mean_error` and `nc_star_gaussian_error` carry the **tight**
//! localisation for contrast — `0.25 Sum_{e in v} l beta / A_v` and
//! `(2pi - alpha_v) / A_v` with `A_v` the barycentric third — because the gap
//! between the two localisations of one measure is the most useful diagnostic in
//! the file, and it is what justifies not benching the vertex-star form as an arm.
//!
//! # Fixtures, and the closed forms they are scored against
//!
//! | fixture | source | exact `H` | exact `K` | `integral H dA` | `chi` |
//! |---|---|---|---|---|---|
//! | `sphere` | `fields::Sphere::canonical()`, `r = 1` | `1/r` | `1/r²` | `4 pi r` | 2 |
//! | `torus` | `fields::Torus::canonical()`, `R = 1`, `a = 0.3` | `(R + 2a cos v) / (2a(R + a cos v))` | `cos v / (a(R + a cos v))` | `2 pi² R` | 0 |
//! | `capsule` | bench-local, `r = 0.47`, `half = 0.565` | barrel `1/(2r)`, caps `1/r` | barrel `0`, caps `1/r²` | `2 pi * half + 4 pi r` | 2 |
//!
//! The torus is derived from the crate's own SDF, whose axis of revolution is
//! **y**, not z (`fields/mod.rs:389-391`: `s = sqrt(x² + z²)`, `q = (s - major,
//! y)`). So `cos v = (s - R) / |q|`, which is also the projection onto the
//! surface — `|q|` is the point's actual distance from the core circle, so a
//! vertex that is a hair off the zero set is scored at its own foot point rather
//! than at a nominal `a`. Derivation: `k_v = 1/a` along the meridian and
//! `k_u = cos v / (R + a cos v)` along the parallel, so
//! `H = (k_u + k_v)/2 = (R + 2a cos v) / (2a(R + a cos v))` and `K = k_u k_v`.
//! Then `dA = a(R + a cos v) du dv` collapses the mean integral to
//! `integral (R + 2a cos v)/2 du dv = 2 pi² R`, independent of the minor radius.
//! `docs/experiments/p-42.csv` already carries that number as
//! `mean_smooth_integral_h = 19.739208802` for the torus, which is `2 pi² R` to
//! nine digits — the closed form is not new here, only its pointwise version is.
//!
//! The capsule is the third fixture because it is the only one that carries a
//! region of **exactly zero** Gaussian curvature next to a region of large
//! positive Gaussian curvature, which is the one regime `sphere` and `torus`
//! cannot test: an estimator that manufactures curvature out of mesh noise is
//! invisible where the true value is `1` and obvious where it is `0`. Its radius
//! is `0.47` and not `0.5` on purpose. At `0.5`, `(0.5, y, 0)` is a grid point
//! with `f` **exactly** zero at every resolution in the ladder, which is the
//! `=`-corner case of `M-352` / `P-53`; at `0.47` the condition
//! `a² + b² = r²/h²` has no integer solution for any `h` in
//! {1/8, 1/12, 1/16, 1/24, 1/32}, since `r²/h²` = 14.14, 31.81, 56.55, 127.24,
//! 226.20. Its two seam circles are genuine `C⁰` curvature discontinuities and
//! penalise **both** mollified arms; `smooth_only_mean_error` isolates the band
//! `| |y| - half | <= eps` so the penalty is priced instead of hidden, and the
//! headline error stays unfiltered.
//!
//! # Arms
//!
//! | arm | what it varies | `is_control` |
//! |---|---|---|
//! | `normal_cycles` | the shipped Sun & Morvan measure, localised over a `rho_eps` ball | **yes** — the incumbent C1 is scored against |
//! | `varifold` | Buet-Leonardi-Masnou's first variation of `Sum area(T) delta_{(x_T,n_T)}` | no |
//! | `field` | `sphere`, `torus`, `capsule` — constant `K`, sign-changing `K`, zero-and-positive `K` | — |
//! | `resolution` | 33, 47, 65, 91, 129 samples per axis | — |
//!
//! Five resolutions and not the house three, because `C2` fits a slope: three
//! points give one degree of freedom for a standard error and the `c2_holds`
//! test below is a standard-error test. 129³ is affordable here because the
//! measured cost is a gather over `O((eps/h)²)` triangles per vertex, which is
//! ~192 at 129³ against `p-42.csv`'s 19,230 vertices, not a re-extraction.
//!
//! ## Why 47 and 91 and not 49 and 97
//!
//! **The rungs are chosen so the fixtures stay transversal to the grid, and the
//! obvious ladder is wrong.** A grid corner where `f` is *exactly* zero is the
//! `=`-corner case of `M-352` / `P-53`. When that corner is a local extremum of
//! `f` along all three grid lines it is harmless — only one incident edge is
//! cut, so only one vertex lands there. When the surface crosses **transversally**
//! through it, several incident edges are cut, marching cubes caches one vertex
//! per grid edge, and the result is several **coincident** vertices: triangles
//! with a side of length zero, which destroys the per-face `Sum alpha = pi`
//! identity the whole angle-defect measure rests on.
//!
//! The unit sphere on `[-2, 2]³` hits this whenever `a² + b² + c² = ((n-1)/4)²`
//! has a non-axis integer solution, which needs `n = 4t + 1`. Measured, by
//! extracting every odd `n` from 17 to 145 and counting zero-length sides: the
//! sphere is **safe at every `n = 4t + 3`**, and among `n = 4t + 1` only where
//! `t²` has the trivial representation alone — `t` a power of two, i.e.
//! `n = 17, 33, 65, 129`. That is exactly why `P-42` chose 33 / 65 / 129. The
//! naive geometric fill-in, 49 and 97, are both `4t + 1` with `t = 12` and
//! `t = 24`, and `144 = 4² + 8² + 8²` and `576 = 8² + 16² + 16²` put a grid
//! corner on the sphere at `(2/3, 2/3, 1/3)`. Measured at 49³: **40 zero-length
//! sides**, and the localised normal-cycle mean-curvature error jumps from
//! `1.02e-2` at 33³ to `8.78e-1` — an 87-fold blow-up, on a rung between two
//! clean ones. The varifold arm barely notices (`1.34e-2`), because it reads only
//! face centroids, normals and areas and never a vertex angle.
//!
//! So the ladder is `33, 47, 65, 91, 129`: `47 = 4(11) + 3` and
//! `91 = 4(22) + 3` are transversal by the residue alone, they sit at the
//! geometric midpoints of `[33, 65]` and `[65, 129]`, and the resulting
//! `cell_size` sequence 0.125, 0.086957, 0.0625, 0.044444, 0.03125 has
//! successive ratios 1.4375, 1.3913, 1.4063, 1.4222 — near-uniform in `ln h`,
//! which is what a least-squares slope wants. `torus` and `capsule` measured
//! clean at all 65 resolutions tested, the torus because `sqrt(x² + z²)` is
//! irrational at almost every grid point and the capsule because of its radius.
//! The `zero_length_sides` control below is what turns this from a comment into
//! a gate.
//!
//! Marching cubes at its shipped defaults (`FaceAmbiguity::Separate`,
//! `InteriorAmbiguity::Ignore`, no crossing refinement), in `f64`, exactly as
//! `experiment_p42.rs:470` runs it — the incumbent has to be measured on the
//! incumbent's mesh.
//!
//! # SHARE, recomputed before the numbers
//!
//! The registration's SHARE is *"C1 moves the curvature stage only, whose share
//! must be reported."* Discharged as two columns per row: `extract_ms` is the
//! marching-cubes stage and `estimator_ms` is the curvature stage, so
//! `curvature_share = estimator_ms / (extract_ms + estimator_ms)` is the share of
//! the pipeline a change of curvature estimator can move. Nothing else in the
//! pipeline is touched by either arm: both consume a finished `MeshBuffer` and
//! neither re-samples the field.
//!
//! # Timing
//!
//! `std::time::Instant`, no criterion. `M-280` measured this host's
//! `amd-pstate-epp` governor swinging one binary 1.45× between runs, and
//! `cost_ratio` is a ratio, so every arm is warmed once and then run five times;
//! `estimator_ms` is the **median** and `estimator_ms_min` / `_max` /
//! `_scatter` carry the spread. `cost_matched` reports whether `cost_ratio`
//! landed in `[0.5, 2]`, and it is **not** a condition on `c1_holds`.
//!
//! Read `cost_ratio`'s direction before reading its band. It is
//! `varifold_ms / normal_cycles_ms`, so **below one means the challenger is the
//! cheaper arm** — which is the structural expectation here, not an accident: the
//! varifold needs one gather over face centroids, while normal cycles needs
//! three, over face centroids and edge midpoints and vertices, because its three
//! measures live on three different geometric supports. A `cost_ratio` under
//! `0.5` therefore makes the registered null **more** conservative, not less: the
//! varifold is being given at least as much budget as the incumbent and still has
//! to win on error. Only a `cost_ratio` above `COST_BAND` would weaken `C1`, and
//! that direction is the one the band is really guarding.
//!
//! Extraction happens once per `(field, resolution)` and is outside both clocks.
//!
//! Determinism: no RNG anywhere. Every sum runs in a fixed order — vertices and
//! faces in index order, edges in `(lo, hi, face)` sorted order as
//! `experiment_p42.rs:281` does, and neighbours in bin order from a counting
//! sort — so the file is reproducible bit for bit.
//!
//! # Verdicts
//!
//! - `error_ratio = varifold.mean_curvature_error / normal_cycles.mean_curvature_error`,
//!   one value per `(field, resolution)`, written on **both** of that pair's rows.
//! - `c1_holds = error_ratio >= 1`, i.e. the registered null held and the
//!   varifold did not win. Per row. `c1_holds_global` is the conjunction.
//! - `convergence_exponent` is the least-squares slope of `ln(mean_curvature_error)`
//!   on `ln(cell_size)` over the five resolutions — one per `(estimator, field)`,
//!   written on that estimator's five rows for that field. Positive means
//!   converging. `gaussian_exponent` is the same fit for the Gaussian error.
//! - `c2_holds` asks whether the two exponents are **distinguishable**, which is
//!   the clause's own falsifier read forwards: the falsifier is *"C2 by
//!   indistinguishable exponents"*. The instrument is the fit itself rather than
//!   a round number — `c2_holds = |exponent_gap| > 2 sqrt(se_nc² + se_var²)`,
//!   a two-sigma separation of two least-squares slopes. Per field, written on
//!   all six of that field's rows. `exponent_gap_sigma` is the separation in
//!   sigmas so a reader can re-decide at another threshold.
//!
//! # Vacuity controls
//!
//! Every one runs before the first `run.record`, and every panic starts `VOID: `.
//!
//! - **the ladder is the regime the bounds are stated in** — `eps` strictly
//!   decreasing and `eps / h` strictly increasing along the ladder, and all three
//!   fixtures on the identical domain. Proves `convergence_exponent` is a fit
//!   through one regime and not a splice of two. Columns `eps`, `eps_over_h`.
//! - **the mesh is a closed manifold** — `boundary_edges`,
//!   `non_manifold_edges` and `zero_length_sides` all zero at every
//!   `(field, resolution)`. A boundary vertex's defect is `pi - alpha_v` and not
//!   `2pi - alpha_v`, and a zero-length side costs the per-face `Sum alpha = pi`
//!   identity, so either one calibrates the incumbent against the wrong
//!   reference. `p-42.csv` reports all three as `0` for sphere and torus at
//!   33/65/129, and this control is not a formality: it **fired during
//!   development** on a first ladder containing 49³, where the canonical unit
//!   sphere is transversal through a grid corner and marching cubes emitted 40
//!   zero-length sides. See `## Why 47 and 91 and not 49 and 97` above; the
//!   ladder is what it is because this assert refused that rung.
//! - **the registered control, half 1: the incumbent reproduces the analytic
//!   curvature globally** — on `sphere` and `torus`, at every resolution,
//!   `|½ Sum l beta - integral H dA| / |integral H dA| <= 0.02` and
//!   `|Sum defect - 2 pi chi| <= 1e-6`. Columns `global_int_h_rel_error`,
//!   `global_defect_measured`. `p-42.csv`'s `mean_half_relative_error` is
//!   1.80e-3 / 4.49e-4 / 1.12e-4 on the sphere, so 0.02 is a ten-fold margin at
//!   the coarsest rung and still fires on a lost factor.
//! - **the registered control, half 2: both arms reproduce it pointwise** — at
//!   the finest rung on `sphere` and `torus`, both arms'
//!   `mean_curvature_error <= 0.25`. An estimator with a worse than 25%
//!   area-weighted relative RMS error on a **constant**-curvature sphere is not
//!   measuring curvature, and `error_ratio` between two such numbers would be a
//!   ratio of two noises.
//! - **the reference is not a constant zero** — `h_exact_rms > 0` and
//!   `k_exact_rms > 0` for every fixture, and every one of the sixty fitted
//!   error values strictly positive and finite. A relative error normalised by
//!   zero, or a zero error fed to `ln`, is not a measurement (`M-44`).
//! - **the torus really does change the sign of `K`** — `k_positive_vertices > 0`
//!   **and** `k_negative_vertices > 0`. Otherwise the one fixture that is
//!   supposed to have varying, sign-changing curvature is secretly a sphere and
//!   `C1` is decided on constant-curvature data alone.
//! - **the capsule really does carry a zero-`K` region** —
//!   `k_zero_exact_vertices > 0` and `k_positive_vertices > 0` on the capsule, so
//!   the `K = 0` regime is populated rather than nominal.
//! - **the mollifier is a convolution and not a handful of terms** — mean
//!   support `>= 8` triangles at every `(field, resolution)`. Below that
//!   `Sum_T area(T) rho_eps` is a small sum wearing a convolution's name and the
//!   varifold arm is being blamed for its neighbour count. Column
//!   `support_triangles_mean`.

#![allow(
    clippy::cast_possible_truncation,
    clippy::cast_precision_loss,
    clippy::too_many_lines
)]

mod common;

use std::time::Instant;

use isomesh::fields::{ReferenceField, Sphere, Torus};
use isomesh::marching_cubes::MarchingCubes;
use isomesh::{MeshBuffer, RuntimeShape3, Sdf};

/// Samples per axis. Five rungs because `C2` fits a slope and needs a standard
/// error; see the header's `# Arms`.
const LADDER: [u32; 5] = [33, 47, 65, 91, 129];

/// Rungs in the ladder, as a `usize` so `as_chunks` can group `cells` by field.
const RUNGS: usize = LADDER.len();

/// `eps = EPS_C * L0 * sqrt(h / L0)` with `L0 = 1`. Buet-Leonardi-Masnou's own
/// balance of their `O(eps) + O(h/eps)`.
const EPS_C: f64 = 0.8;

/// Unit length the three fixtures live at, so the `sqrt` is dimensionally sound.
const L0: f64 = 1.0;

/// Timed repeats per arm per row. `M-280`: this host's governor swings one
/// binary 1.45×, and `cost_ratio` is a ratio.
const REPEATS: usize = 5;

/// Half-span of the crate's `COMPACT_DOMAIN`, which `Sphere` and `Torus` return
/// from `ReferenceField::domain()`. Asserted against theirs, never assumed.
const DOMAIN_HALF_SPAN: f64 = 2.0;

/// Ceiling on an arm's area-weighted relative RMS mean-curvature error before it
/// counts as calibrated.
const CALIBRATION_CEILING: f64 = 0.25;

/// Ceiling on the incumbent's **global** relative mean-curvature error.
const GLOBAL_MEAN_TOLERANCE: f64 = 0.02;

/// `Sum_v (2pi - alpha_v)` is `2 pi chi` as an identity, so this is a float
/// tolerance and nothing more.
const GAUSS_BONNET_TOLERANCE: f64 = 1e-6;

/// `cost_ratio` inside `[1/COST_BAND, COST_BAND]` counts as matched cost.
const COST_BAND: f64 = 2.0;

/// Fewer triangles than this under the mollifier and it is not a convolution.
const MIN_SUPPORT: f64 = 8.0;

/// Sigmas of slope separation required for `c2_holds`.
const SIGMA_GATE: f64 = 2.0;

/// Capsule tube radius. `0.47` and not `0.5`: see the header.
const CAPSULE_RADIUS: f64 = 0.47;

/// Capsule half-length along `y`, so the seam circles sit at `|y| = 0.565`.
const CAPSULE_HALF_LENGTH: f64 = 0.565;

// ─── small vector helpers ───────────────────────────────────────────────────

fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn scale(a: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}

fn norm(a: [f64; 3]) -> f64 {
    dot(a, a).sqrt()
}

/// Unit vector, or the zero vector for a zero input.
///
/// Zero rather than `NaN` on purpose, and for the same reason
/// `experiment_p42.rs:293-295` gives: a zero-area face then contributes
/// `atan2(0, 0) = 0` and a zero weight to every sum, instead of poisoning the
/// accumulator it lands in.
fn unit(a: [f64; 3]) -> [f64; 3] {
    let l = norm(a);
    if l > 0.0 {
        scale(a, l.recip())
    } else {
        [0.0; 3]
    }
}

/// Interior angle at `apex` in the triangle `apex, u, w`.
fn corner_angle(apex: [f64; 3], u: [f64; 3], w: [f64; 3]) -> f64 {
    let a = sub(u, apex);
    let b = sub(w, apex);
    let denom = norm(a) * norm(b);
    if denom > 0.0 {
        (dot(a, b) / denom).clamp(-1.0, 1.0).acos()
    } else {
        0.0
    }
}

// ─── the mollifier ─────────────────────────────────────────────────────────

/// `xi(t) = (1 - t²)³` on `[0, 1]`, zero outside. `C²` at `t = 1`.
fn xi(t: f64) -> f64 {
    if t >= 1.0 {
        0.0
    } else {
        let s = 1.0 - t * t;
        s * s * s
    }
}

/// `xi'(t) = -6 t (1 - t²)²`. Zero at both ends of the support.
fn xi_prime(t: f64) -> f64 {
    if t >= 1.0 {
        0.0
    } else {
        let s = 1.0 - t * t;
        -6.0 * t * s * s
    }
}

// ─── fixtures ───────────────────────────────────────────────────────────────

/// A round-ended cylinder along `y`. Exact distance field, closed, and the only
/// fixture with a region of exactly zero Gaussian curvature.
struct Capsule {
    radius: f64,
    half_length: f64,
}

impl Sdf for Capsule {
    type Scalar = f64;

    fn sample(&self, p: [f64; 3]) -> f64 {
        let t = p[1].clamp(-self.half_length, self.half_length);
        norm([p[0], p[1] - t, p[2]]) - self.radius
    }

    fn gradient(&self, p: [f64; 3]) -> [f64; 3] {
        let t = p[1].clamp(-self.half_length, self.half_length);
        unit([p[0], p[1] - t, p[2]])
    }
}

/// The three fields with a closed-form curvature.
///
/// An `enum` and not a generic parameter because `ReferenceField` demands a
/// single `const NAME` per type and these three have to be walked as one list.
enum Fixture {
    Sphere(Sphere<f64>),
    Torus(Torus<f64>),
    Capsule(Capsule),
}

impl Sdf for Fixture {
    type Scalar = f64;

    fn sample(&self, p: [f64; 3]) -> f64 {
        match self {
            Self::Sphere(s) => s.sample(p),
            Self::Torus(t) => t.sample(p),
            Self::Capsule(c) => c.sample(p),
        }
    }

    fn gradient(&self, p: [f64; 3]) -> [f64; 3] {
        match self {
            Self::Sphere(s) => s.gradient(p),
            Self::Torus(t) => t.gradient(p),
            Self::Capsule(c) => c.gradient(p),
        }
    }
}

impl Fixture {
    fn name(&self) -> &'static str {
        match self {
            Self::Sphere(_) => "sphere",
            Self::Torus(_) => "torus",
            Self::Capsule(_) => "capsule",
        }
    }

    /// Sampling domain. `Sphere` and `Torus` answer from the crate; the capsule
    /// answers from `DOMAIN_HALF_SPAN`, and the ladder control asserts the three
    /// agree so the literal cannot drift from `fields/mod.rs`.
    fn domain(&self) -> ([f64; 3], [f64; 3]) {
        match self {
            Self::Sphere(s) => s.domain(),
            Self::Torus(t) => t.domain(),
            Self::Capsule(_) => ([-DOMAIN_HALF_SPAN; 3], [DOMAIN_HALF_SPAN; 3]),
        }
    }

    /// `(H, K)` in closed form at the foot point of `p`. Derivations in the
    /// header's fixture table.
    fn exact_curvature(&self, p: [f64; 3]) -> (f64, f64) {
        match self {
            Self::Sphere(s) => (s.radius.recip(), s.radius.recip() * s.radius.recip()),
            Self::Torus(t) => {
                // The crate revolves about y: fields/mod.rs:389-391.
                let s = (p[0] * p[0] + p[2] * p[2]).sqrt();
                let dr = s - t.major;
                let q = (dr * dr + p[1] * p[1]).sqrt();
                let cos_v = dr / q;
                let denom = t.major + t.minor * cos_v;
                let h = (t.major + 2.0 * t.minor * cos_v) / (2.0 * t.minor * denom);
                (h, cos_v / (t.minor * denom))
            }
            Self::Capsule(c) => {
                if p[1].abs() <= c.half_length {
                    (0.5 / c.radius, 0.0)
                } else {
                    (c.radius.recip(), c.radius.recip() * c.radius.recip())
                }
            }
        }
    }

    /// `integral H dA` over the whole closed surface.
    fn integral_mean_curvature(&self) -> f64 {
        match self {
            Self::Sphere(s) => 4.0 * std::f64::consts::PI * s.radius,
            Self::Torus(t) => 2.0 * std::f64::consts::PI * std::f64::consts::PI * t.major,
            Self::Capsule(c) => {
                // Barrel: (1/2r)(2 pi r L) = pi L, with L = 2 * half_length.
                // Caps:   (1/r)(4 pi r²)  = 4 pi r.
                std::f64::consts::TAU * c.half_length + 4.0 * std::f64::consts::PI * c.radius
            }
        }
    }

    /// Euler characteristic, so `Sum defect` has something to be checked against.
    fn euler(&self) -> i64 {
        match self {
            Self::Sphere(_) | Self::Capsule(_) => 2,
            Self::Torus(_) => 0,
        }
    }

    /// `true` where the exact curvature is continuous at scale `eps`. Always
    /// `true` off the capsule; on the capsule it excludes the seam band.
    fn is_smooth_at(&self, p: [f64; 3], eps: f64) -> bool {
        match self {
            Self::Sphere(_) | Self::Torus(_) => true,
            Self::Capsule(c) => (p[1].abs() - c.half_length).abs() > eps,
        }
    }
}

// ─── mesh topology, shared by both arms ─────────────────────────────────────

/// A triangle with three distinct in-range indices.
struct Face {
    verts: [u32; 3],
    /// Outward unit normal from the crate's counter-clockwise-seen-from-outside
    /// winding (`lib.rs:56-67`).
    normal: [f64; 3],
    area: f64,
    centroid: [f64; 3],
}

/// An occurrence of an undirected edge in one face.
struct EdgeRef {
    lo: u32,
    hi: u32,
    face: u32,
    /// `true` when this face traverses the edge `lo -> hi`.
    forward: bool,
}

/// Everything both arms read off the mesh, computed once per `(field, resolution)`.
struct Topology {
    faces: Vec<Face>,
    centroids: Vec<[f64; 3]>,
    /// Barycentric third of the incident face areas.
    vertex_area: Vec<f64>,
    /// `2 pi - alpha_v` for referenced vertices, `0` elsewhere.
    vertex_defect: Vec<f64>,
    referenced: Vec<bool>,
    edge_mid: Vec<[f64; 3]>,
    /// `l(e) * beta(e)`, parallel to `edge_mid`.
    edge_measure: Vec<f64>,
    /// `Sum_{e incident to v} l(e) beta(e)`, for the vertex-star diagnostic.
    vertex_edge_measure: Vec<f64>,
    area_total: f64,
    mean_measure_total: f64,
    defect_total: f64,
    boundary_edges: u64,
    non_manifold_edges: u64,
    zero_length_sides: u64,
}

/// `experiment_p42.rs:296`'s signed dihedral: `+` convex, `-` concave.
///
/// `atan2((n1 x n2) . e_hat, n1 . n2)` with `e_hat` the edge direction as the
/// **first** face traverses it. Swapping the two faces flips `e_hat` and
/// `n1 x n2` together, so the value does not depend on the sort order the group
/// arrived in.
fn signed_dihedral(first: &Face, second: &Face, edge_dir: [f64; 3]) -> f64 {
    dot(cross(first.normal, second.normal), edge_dir).atan2(dot(first.normal, second.normal))
}

impl Topology {
    fn build(mesh: &MeshBuffer<f64>) -> Self {
        let nv = mesh.positions.len();
        let mut faces = Vec::with_capacity(mesh.indices.len() / 3);
        let mut zero_length_sides = 0u64;

        // Same skip predicate as experiment_p42.rs:236-242, so this accumulator
        // and p-42.csv describe the same face set.
        for tri in mesh.indices.as_chunks::<3>().0 {
            let (a, b, c) = (tri[0], tri[1], tri[2]);
            let in_range = (a as usize) < nv && (b as usize) < nv && (c as usize) < nv;
            if !in_range || a == b || b == c || c == a {
                continue;
            }
            let pa = mesh.positions[a as usize];
            let pb = mesh.positions[b as usize];
            let pc = mesh.positions[c as usize];
            let twice = cross(sub(pb, pa), sub(pc, pa));
            let sides = [norm(sub(pb, pc)), norm(sub(pc, pa)), norm(sub(pa, pb))];
            if sides.iter().any(|&s| s <= 0.0) {
                zero_length_sides += 1;
            }
            faces.push(Face {
                verts: [a, b, c],
                normal: unit(twice),
                area: norm(twice) / 2.0,
                centroid: [
                    (pa[0] + pb[0] + pc[0]) / 3.0,
                    (pa[1] + pb[1] + pc[1]) / 3.0,
                    (pa[2] + pb[2] + pc[2]) / 3.0,
                ],
            });
        }

        let mut vertex_area = vec![0.0_f64; nv];
        let mut angle_sum = vec![0.0_f64; nv];
        let mut referenced = vec![false; nv];
        let mut area_total = 0.0_f64;

        for face in &faces {
            let [a, b, c] = face.verts;
            let pa = mesh.positions[a as usize];
            let pb = mesh.positions[b as usize];
            let pc = mesh.positions[c as usize];
            area_total += face.area;
            for (v, apex, u, w) in [(a, pa, pb, pc), (b, pb, pc, pa), (c, pc, pa, pb)] {
                angle_sum[v as usize] += corner_angle(apex, u, w);
                vertex_area[v as usize] += face.area / 3.0;
                referenced[v as usize] = true;
            }
        }

        let mut refs: Vec<EdgeRef> = Vec::with_capacity(faces.len() * 3);
        for (fi, face) in faces.iter().enumerate() {
            let [a, b, c] = face.verts;
            for (u, v) in [(a, b), (b, c), (c, a)] {
                refs.push(EdgeRef {
                    lo: u.min(v),
                    hi: u.max(v),
                    face: fi as u32,
                    forward: u < v,
                });
            }
        }
        refs.sort_unstable_by_key(|e| (e.lo, e.hi, e.face));

        let mut edge_mid = Vec::new();
        let mut edge_measure = Vec::new();
        let mut vertex_edge_measure = vec![0.0_f64; nv];
        let mut boundary_edges = 0u64;
        let mut non_manifold_edges = 0u64;
        let mut mean_measure_total = 0.0_f64;

        let mut i = 0usize;
        while i < refs.len() {
            let mut j = i + 1;
            while j < refs.len() && refs[j].lo == refs[i].lo && refs[j].hi == refs[i].hi {
                j += 1;
            }
            let group = &refs[i..j];
            let (lo, hi) = (group[0].lo as usize, group[0].hi as usize);
            match group.len() {
                1 => boundary_edges += 1,
                2 => {
                    let first = &faces[group[0].face as usize];
                    let second = &faces[group[1].face as usize];
                    let along = sub(mesh.positions[hi], mesh.positions[lo]);
                    let dir = if group[0].forward {
                        unit(along)
                    } else {
                        unit(scale(along, -1.0))
                    };
                    let contribution = norm(along) * signed_dihedral(first, second, dir);
                    mean_measure_total += contribution;
                    vertex_edge_measure[lo] += contribution;
                    vertex_edge_measure[hi] += contribution;
                    let (a, b) = (mesh.positions[lo], mesh.positions[hi]);
                    edge_mid.push([
                        f64::midpoint(a[0], b[0]),
                        f64::midpoint(a[1], b[1]),
                        f64::midpoint(a[2], b[2]),
                    ]);
                    edge_measure.push(contribution);
                }
                _ => non_manifold_edges += 1,
            }
            i = j;
        }

        let mut vertex_defect = vec![0.0_f64; nv];
        let mut defect_total = 0.0_f64;
        for (v, (&sum, &seen)) in angle_sum.iter().zip(&referenced).enumerate() {
            if seen {
                vertex_defect[v] = std::f64::consts::TAU - sum;
                defect_total += vertex_defect[v];
            }
        }

        let centroids = faces.iter().map(|f| f.centroid).collect();

        Self {
            faces,
            centroids,
            vertex_area,
            vertex_defect,
            referenced,
            edge_mid,
            edge_measure,
            vertex_edge_measure,
            area_total,
            mean_measure_total,
            defect_total,
            boundary_edges,
            non_manifold_edges,
            zero_length_sides,
        }
    }
}

// ─── the shared neighbour search ────────────────────────────────────────────

/// A uniform bin grid with cell size `eps`, so a `rho_eps` gather visits 27 bins.
///
/// One counting sort, one `Vec<u32>` of items. Both arms use it, which is half of
/// why `cost_ratio` is near one for structural reasons rather than tuned ones.
struct Bins {
    origin: [f64; 3],
    inv_cell: f64,
    dims: [i64; 3],
    starts: Vec<u32>,
    items: Vec<u32>,
}

impl Bins {
    fn build(points: &[[f64; 3]], cell: f64) -> Self {
        let mut lo = [f64::INFINITY; 3];
        let mut hi = [f64::NEG_INFINITY; 3];
        for p in points {
            for a in 0..3 {
                lo[a] = lo[a].min(p[a]);
                hi[a] = hi[a].max(p[a]);
            }
        }
        if points.is_empty() {
            lo = [0.0; 3];
            hi = [0.0; 3];
        }
        let inv_cell = cell.recip();
        let mut dims = [1i64; 3];
        for a in 0..3 {
            dims[a] = (((hi[a] - lo[a]) * inv_cell).floor() as i64 + 1).max(1);
        }
        let ncells = (dims[0] * dims[1] * dims[2]) as usize;

        let cell_of = |p: &[f64; 3]| -> usize {
            let mut idx = 0i64;
            for a in 0..3 {
                let c = (((p[a] - lo[a]) * inv_cell).floor() as i64).clamp(0, dims[a] - 1);
                idx = idx * dims[a] + c;
            }
            idx as usize
        };

        let mut starts = vec![0u32; ncells + 1];
        for p in points {
            starts[cell_of(p) + 1] += 1;
        }
        for c in 0..ncells {
            starts[c + 1] += starts[c];
        }
        let mut cursor = starts.clone();
        let mut items = vec![0u32; points.len()];
        for (i, p) in points.iter().enumerate() {
            let c = cell_of(p);
            items[cursor[c] as usize] = i as u32;
            cursor[c] += 1;
        }

        Self {
            origin: lo,
            inv_cell,
            dims,
            starts,
            items,
        }
    }

    /// Indices of every point in the 27 bins around `x`, into `out`.
    ///
    /// The query's own bin index is deliberately **not** clamped: a query just
    /// outside the point cloud must still see the boundary bin, and one that is
    /// far outside must see nothing rather than the nearest wall.
    fn gather(&self, x: [f64; 3], out: &mut Vec<u32>) {
        out.clear();
        let mut base = [0i64; 3];
        for a in 0..3 {
            base[a] = ((x[a] - self.origin[a]) * self.inv_cell).floor() as i64;
        }
        for dz in -1..=1 {
            for dy in -1..=1 {
                for dx in -1..=1 {
                    let g = [base[0] + dx, base[1] + dy, base[2] + dz];
                    if (0..3).any(|a| g[a] < 0 || g[a] >= self.dims[a]) {
                        continue;
                    }
                    let idx = ((g[0] * self.dims[1] + g[1]) * self.dims[2] + g[2]) as usize;
                    out.extend_from_slice(
                        &self.items[self.starts[idx] as usize..self.starts[idx + 1] as usize],
                    );
                }
            }
        }
    }
}

// ─── arm 1: the varifold ────────────────────────────────────────────────────

/// Per-vertex output of one varifold pass.
struct VarifoldField {
    mean: Vec<f64>,
    gaussian: Vec<f64>,
    /// `0.5 tr S`, the independent cross-check on `mean`.
    mean_from_shape: Vec<f64>,
    mass: Vec<f64>,
    /// The varifold's own mollified unit normal.
    normal: Vec<[f64; 3]>,
    support: Vec<u32>,
}

/// One gather over the triangle centroids: first variation **and** the analytic
/// tangential derivative of the mollified normal, from 19 accumulators.
fn varifold_pass(
    positions: &[[f64; 3]],
    topo: &Topology,
    bins: &Bins,
    eps: f64,
    scratch: &mut Vec<u32>,
) -> VarifoldField {
    let n = positions.len();
    let mut out = VarifoldField {
        mean: vec![0.0; n],
        gaussian: vec![0.0; n],
        mean_from_shape: vec![0.0; n],
        mass: vec![0.0; n],
        normal: vec![[0.0; 3]; n],
        support: vec![0; n],
    };

    for (vi, &x) in positions.iter().enumerate() {
        bins.gather(x, scratch);
        let mut mass = 0.0_f64;
        // N = Sum A rho n
        let mut nvec = [0.0_f64; 3];
        // grad_of_n[a][b] = d_b N_a = Sum A (grad rho)_b n_a
        let mut grad_of_n = [[0.0_f64; 3]; 3];
        // Sum A * P_{n_T} grad rho
        let mut projected = [0.0_f64; 3];
        let mut support = 0u32;

        for &fi in scratch.iter() {
            let face = &topo.faces[fi as usize];
            let z = sub(x, face.centroid);
            let d = norm(z);
            if d >= eps || face.area <= 0.0 {
                continue;
            }
            support += 1;
            let t = d / eps;
            let rho = xi(t);
            let weight = face.area * rho;
            mass += weight;
            for (slot, &component) in nvec.iter_mut().zip(&face.normal) {
                *slot += weight * component;
            }
            // grad rho_eps(z) = xi'(t) * z_hat / eps
            let grad = scale(unit(z), xi_prime(t) / eps);
            let along_normal = dot(face.normal, grad);
            for a in 0..3 {
                projected[a] += face.area * (grad[a] - along_normal * face.normal[a]);
                for b in 0..3 {
                    grad_of_n[a][b] += face.area * grad[b] * face.normal[a];
                }
            }
        }

        out.support[vi] = support;
        out.mass[vi] = mass;
        let len_n = norm(nvec);
        if mass <= 0.0 || len_n <= 0.0 {
            continue;
        }
        let nhat = scale(nvec, len_n.recip());
        out.normal[vi] = nhat;

        // The registered form: Hvec = -num/den, which is (k1 + k2) n_outward.
        let hvec = scale(projected, -mass.recip());
        out.mean[vi] = 0.5 * dot(hvec, nhat);

        // grad nhat = (I/|N| - N (x) N/|N|³) grad N, exactly.
        let mut dn = [[0.0_f64; 3]; 3];
        let inv = len_n.recip();
        let inv3 = inv * inv * inv;
        for b in 0..3 {
            let mut ndotcol = 0.0_f64;
            for c in 0..3 {
                ndotcol += nvec[c] * grad_of_n[c][b];
            }
            for a in 0..3 {
                dn[a][b] = grad_of_n[a][b] * inv - nvec[a] * ndotcol * inv3;
            }
        }

        // Deterministic orthonormal tangent basis. `det S` and `tr S` are
        // invariant under this choice, so it cannot move a recorded number.
        let axis = {
            let (ax, ay, az) = (nhat[0].abs(), nhat[1].abs(), nhat[2].abs());
            if ax <= ay && ax <= az {
                [1.0, 0.0, 0.0]
            } else if ay <= az {
                [0.0, 1.0, 0.0]
            } else {
                [0.0, 0.0, 1.0]
            }
        };
        let t1 = unit(sub(axis, scale(nhat, dot(axis, nhat))));
        let t2 = cross(nhat, t1);
        let apply = |t: [f64; 3]| -> [f64; 3] {
            let mut r = [0.0_f64; 3];
            for a in 0..3 {
                for b in 0..3 {
                    r[a] += dn[a][b] * t[b];
                }
            }
            r
        };
        let c1 = apply(t1);
        let c2 = apply(t2);
        let s11 = dot(t1, c1);
        let s12 = dot(t1, c2);
        let s21 = dot(t2, c1);
        let s22 = dot(t2, c2);
        out.gaussian[vi] = s11 * s22 - s12 * s21;
        out.mean_from_shape[vi] = 0.5 * (s11 + s22);
    }

    out
}

// ─── arm 2: the normal cycles ───────────────────────────────────────────────

/// Per-vertex output of one normal-cycle pass.
struct NormalCycleField {
    mean: Vec<f64>,
    gaussian: Vec<f64>,
    mass: Vec<f64>,
}

/// The three bin grids the incumbent's three measures live on.
struct NcBins {
    area: Bins,
    edge: Bins,
    vertex: Bins,
}

/// The shipped Sun & Morvan measures over a `rho_eps` ball, divided by the mass
/// of the same ball.
fn normal_cycle_pass(
    positions: &[[f64; 3]],
    topo: &Topology,
    bins: &NcBins,
    eps: f64,
    scratch: &mut Vec<u32>,
) -> NormalCycleField {
    let n = positions.len();
    let mut out = NormalCycleField {
        mean: vec![0.0; n],
        gaussian: vec![0.0; n],
        mass: vec![0.0; n],
    };

    for (vi, &x) in positions.iter().enumerate() {
        // den(x) = Sum_T area(T) rho_eps(x - x_T) -> integral rho dA
        bins.area.gather(x, scratch);
        let mut mass = 0.0_f64;
        for &fi in scratch.iter() {
            let face = &topo.faces[fi as usize];
            let d = norm(sub(x, face.centroid));
            if d < eps {
                mass += face.area * xi(d / eps);
            }
        }
        out.mass[vi] = mass;
        if mass <= 0.0 {
            continue;
        }
        let inv_mass = mass.recip();

        // 0.5 Sum_e rho_eps(x - m_e) l(e) beta(e) / den, the factor from
        // experiment_p42.rs:61-66's cylinder calibration.
        bins.edge.gather(x, scratch);
        let mut mean = 0.0_f64;
        for &ei in scratch.iter() {
            let d = norm(sub(x, topo.edge_mid[ei as usize]));
            if d < eps {
                mean += xi(d / eps) * topo.edge_measure[ei as usize];
            }
        }
        out.mean[vi] = 0.5 * mean * inv_mass;

        // Sum_v rho_eps(x - p_v) (2pi - alpha_v) / den.
        bins.vertex.gather(x, scratch);
        let mut gaussian = 0.0_f64;
        for &wi in scratch.iter() {
            let w = wi as usize;
            if !topo.referenced[w] {
                continue;
            }
            let d = norm(sub(x, positions[w]));
            if d < eps {
                gaussian += xi(d / eps) * topo.vertex_defect[w];
            }
        }
        out.gaussian[vi] = gaussian * inv_mass;
    }

    out
}

// ─── error metric and slope fit ─────────────────────────────────────────────

/// Area-weighted relative RMS error.
///
/// `sqrt(Sum w (est - exact)² / Sum w) / sqrt(Sum w exact² / Sum w)`, which
/// cancels to `sqrt(Sum w (est - exact)² / Sum w exact²)`. Dimensionless, so the
/// three fixtures are comparable, and area-weighted so a cloud of sliver
/// vertices cannot outvote the surface.
fn relative_rms(weight: &[f64], est: &[f64], exact: &[f64], used: &[bool]) -> f64 {
    let mut residual = 0.0_f64;
    let mut reference = 0.0_f64;
    for (i, &u) in used.iter().enumerate() {
        if !u {
            continue;
        }
        let d = est[i] - exact[i];
        residual += weight[i] * d * d;
        reference += weight[i] * exact[i] * exact[i];
    }
    if reference > 0.0 {
        (residual / reference).sqrt()
    } else {
        f64::NAN
    }
}

/// Least-squares `ln(error) = c + exponent * ln(cell_size)`, with the standard
/// error of the slope.
///
/// A positive `exponent` means the error falls as the grid refines.
struct Fit {
    exponent: f64,
    stderr: f64,
}

fn fit_power_law(cell_size: &[f64], error: &[f64]) -> Fit {
    let n = cell_size.len();
    let xs: Vec<f64> = cell_size.iter().map(|h| h.ln()).collect();
    let ys: Vec<f64> = error.iter().map(|e| e.ln()).collect();
    let xm = xs.iter().sum::<f64>() / n as f64;
    let ym = ys.iter().sum::<f64>() / n as f64;
    let sxx: f64 = xs.iter().map(|x| (x - xm) * (x - xm)).sum();
    let sxy: f64 = xs.iter().zip(&ys).map(|(x, y)| (x - xm) * (y - ym)).sum();
    let exponent = sxy / sxx;
    let ssr: f64 = xs
        .iter()
        .zip(&ys)
        .map(|(x, y)| {
            let r = y - (ym + exponent * (x - xm));
            r * r
        })
        .sum();
    Fit {
        exponent,
        stderr: (ssr / (n as f64 - 2.0) / sxx).sqrt(),
    }
}

/// The four slopes one field's ladder supports, and C2's verdict from two of them.
///
/// `convergence_exponent` is a property of an `(estimator, field)` pair and is
/// written on all five of that pair's resolution rows; `c2_holds` is a property
/// of the field and is written on all ten of its rows.
struct FieldFit {
    nc: Fit,
    varifold: Fit,
    nc_gaussian: Fit,
    varifold_gaussian: Fit,
    /// `|exponent_varifold - exponent_nc|`.
    gap: f64,
    /// `sqrt(se_nc² + se_varifold²)`, the combined standard error of the two
    /// slopes, which is the scale `gap` has to be judged against.
    sigma: f64,
    c2: bool,
}

impl FieldFit {
    fn of(window: &[Cell; RUNGS]) -> Self {
        let hs: Vec<f64> = window.iter().map(|c| c.cell_size).collect();
        let pluck = |f: fn(&Cell) -> f64| -> Vec<f64> { window.iter().map(f).collect() };
        let nc = fit_power_law(&hs, &pluck(|c| c.nc.mean_error));
        let varifold = fit_power_law(&hs, &pluck(|c| c.varifold.mean_error));
        let gap = (varifold.exponent - nc.exponent).abs();
        let sigma = (nc.stderr * nc.stderr + varifold.stderr * varifold.stderr).sqrt();
        Self {
            nc_gaussian: fit_power_law(&hs, &pluck(|c| c.nc.gaussian_error)),
            varifold_gaussian: fit_power_law(&hs, &pluck(|c| c.varifold.gaussian_error)),
            nc,
            varifold,
            gap,
            sigma,
            c2: gap > SIGMA_GATE * sigma,
        }
    }
}

/// Median, min and max of the timed repeats, in milliseconds.
struct Clock {
    median: f64,
    min: f64,
    max: f64,
}

fn clock_of(mut samples: Vec<f64>) -> Clock {
    samples.sort_unstable_by(f64::total_cmp);
    Clock {
        median: samples[samples.len() / 2],
        min: samples[0],
        max: samples[samples.len() - 1],
    }
}

// ─── one measurement ────────────────────────────────────────────────────────

/// What one arm produced at one `(field, resolution)`.
struct Arm {
    mean_error: f64,
    gaussian_error: f64,
    smooth_only_mean_error: f64,
    clock: Clock,
}

/// One `(field, resolution)` cell, both arms and every census that scores it.
struct Cell {
    field: &'static str,
    samples: u32,
    cell_size: f64,
    eps: f64,
    vertices: usize,
    triangles: usize,
    used_vertices: usize,
    mesh_area: f64,
    boundary_edges: u64,
    non_manifold_edges: u64,
    zero_length_sides: u64,
    support_mean: f64,
    h_exact_rms: f64,
    k_exact_rms: f64,
    k_positive: u64,
    k_negative: u64,
    k_zero_exact: u64,
    global_int_h: f64,
    global_int_h_exact: f64,
    global_defect: f64,
    global_defect_exact: f64,
    normal_alignment: f64,
    varifold_h_shape_error: f64,
    nc_star_mean_error: f64,
    nc_star_gaussian_error: f64,
    extract_ms: f64,
    nc: Arm,
    varifold: Arm,
}

fn measure(fixture: &Fixture, samples: u32) -> Cell {
    // The same grid arithmetic as benches/common/mod.rs:44-53, spelled out
    // because `Fixture` cannot implement `ReferenceField` (one `const NAME` per
    // type, three variants).
    let (lo, hi) = fixture.domain();
    let cell_size = (hi[0] - lo[0]) / f64::from(samples - 1);
    let shape = RuntimeShape3::new([samples; 3]).expect("benchmark grid fits u32");

    let mut mesh = MeshBuffer::<f64>::new();
    let mut mc = MarchingCubes::<f64>::new();
    let started = Instant::now();
    mc.extract(fixture, &shape, lo, cell_size, &mut mesh)
        .expect("marching cubes extracts the fixture on its own canonical domain");
    let extract_ms = started.elapsed().as_secs_f64() * 1e3;

    let topo = Topology::build(&mesh);
    let eps = EPS_C * L0 * (cell_size / L0).sqrt();

    let bins_area = Bins::build(&topo.centroids, eps);
    let nc_bins = NcBins {
        area: Bins::build(&topo.centroids, eps),
        edge: Bins::build(&topo.edge_mid, eps),
        vertex: Bins::build(&mesh.positions, eps),
    };

    let mut scratch: Vec<u32> = Vec::with_capacity(4096);

    // Warm up once, then REPEATS timed runs. The values are deterministic, so
    // the last run's output is every run's output.
    let mut varifold = varifold_pass(&mesh.positions, &topo, &bins_area, eps, &mut scratch);
    let mut varifold_ms = Vec::with_capacity(REPEATS);
    for _ in 0..REPEATS {
        let t = Instant::now();
        varifold = varifold_pass(&mesh.positions, &topo, &bins_area, eps, &mut scratch);
        varifold_ms.push(t.elapsed().as_secs_f64() * 1e3);
    }

    let mut nc = normal_cycle_pass(&mesh.positions, &topo, &nc_bins, eps, &mut scratch);
    let mut nc_ms = Vec::with_capacity(REPEATS);
    for _ in 0..REPEATS {
        let t = Instant::now();
        nc = normal_cycle_pass(&mesh.positions, &topo, &nc_bins, eps, &mut scratch);
        nc_ms.push(t.elapsed().as_secs_f64() * 1e3);
    }

    // Exact curvature at every vertex, and the masks that score it.
    let n = mesh.positions.len();
    let mut h_exact = vec![0.0_f64; n];
    let mut k_exact = vec![0.0_f64; n];
    let mut used = vec![false; n];
    let mut smooth = vec![false; n];
    let mut nc_star_mean = vec![0.0_f64; n];
    let mut nc_star_gaussian = vec![0.0_f64; n];
    let (mut k_positive, mut k_negative, mut k_zero_exact) = (0u64, 0u64, 0u64);
    let mut align_num = 0.0_f64;
    let mut align_den = 0.0_f64;

    for (vi, &p) in mesh.positions.iter().enumerate() {
        let (h, k) = fixture.exact_curvature(p);
        h_exact[vi] = h;
        k_exact[vi] = k;
        let live = topo.referenced[vi]
            && topo.vertex_area[vi] > 0.0
            && varifold.mass[vi] > 0.0
            && nc.mass[vi] > 0.0
            && h.is_finite()
            && k.is_finite();
        used[vi] = live;
        if !live {
            continue;
        }
        smooth[vi] = fixture.is_smooth_at(p, eps);
        if k > 0.0 {
            k_positive += 1;
        } else if k < 0.0 {
            k_negative += 1;
        } else {
            k_zero_exact += 1;
        }
        // The tight localisation of the same measure, for contrast.
        nc_star_mean[vi] = 0.25 * topo.vertex_edge_measure[vi] / topo.vertex_area[vi];
        nc_star_gaussian[vi] = topo.vertex_defect[vi] / topo.vertex_area[vi];
        align_num += topo.vertex_area[vi] * dot(varifold.normal[vi], unit(fixture.gradient(p)));
        align_den += topo.vertex_area[vi];
    }

    let w = &topo.vertex_area;
    let smooth_used: Vec<bool> = used.iter().zip(&smooth).map(|(&u, &s)| u && s).collect();

    let mut h_ref_sq = 0.0_f64;
    let mut k_ref_sq = 0.0_f64;
    let mut weight_sum = 0.0_f64;
    let mut used_vertices = 0usize;
    for (vi, &u) in used.iter().enumerate() {
        if u {
            h_ref_sq += w[vi] * h_exact[vi] * h_exact[vi];
            k_ref_sq += w[vi] * k_exact[vi] * k_exact[vi];
            weight_sum += w[vi];
            used_vertices += 1;
        }
    }

    let support_mean = if used_vertices > 0 {
        used.iter()
            .enumerate()
            .filter(|&(_, &u)| u)
            .map(|(vi, _)| f64::from(varifold.support[vi]))
            .sum::<f64>()
            / used_vertices as f64
    } else {
        0.0
    };

    Cell {
        field: fixture.name(),
        samples,
        cell_size,
        eps,
        vertices: n,
        triangles: topo.faces.len(),
        used_vertices,
        mesh_area: topo.area_total,
        boundary_edges: topo.boundary_edges,
        non_manifold_edges: topo.non_manifold_edges,
        zero_length_sides: topo.zero_length_sides,
        support_mean,
        h_exact_rms: if weight_sum > 0.0 {
            (h_ref_sq / weight_sum).sqrt()
        } else {
            0.0
        },
        k_exact_rms: if weight_sum > 0.0 {
            (k_ref_sq / weight_sum).sqrt()
        } else {
            0.0
        },
        k_positive,
        k_negative,
        k_zero_exact,
        global_int_h: 0.5 * topo.mean_measure_total,
        global_int_h_exact: fixture.integral_mean_curvature(),
        global_defect: topo.defect_total,
        global_defect_exact: std::f64::consts::TAU * fixture.euler() as f64,
        normal_alignment: if align_den > 0.0 {
            align_num / align_den
        } else {
            f64::NAN
        },
        varifold_h_shape_error: relative_rms(w, &varifold.mean_from_shape, &h_exact, &used),
        nc_star_mean_error: relative_rms(w, &nc_star_mean, &h_exact, &used),
        nc_star_gaussian_error: relative_rms(w, &nc_star_gaussian, &k_exact, &used),
        extract_ms,
        nc: Arm {
            mean_error: relative_rms(w, &nc.mean, &h_exact, &used),
            gaussian_error: relative_rms(w, &nc.gaussian, &k_exact, &used),
            smooth_only_mean_error: relative_rms(w, &nc.mean, &h_exact, &smooth_used),
            clock: clock_of(nc_ms),
        },
        varifold: Arm {
            mean_error: relative_rms(w, &varifold.mean, &h_exact, &used),
            gaussian_error: relative_rms(w, &varifold.gaussian, &k_exact, &used),
            smooth_only_mean_error: relative_rms(w, &varifold.mean, &h_exact, &smooth_used),
            clock: clock_of(varifold_ms),
        },
    }
}

// ─── main ───────────────────────────────────────────────────────────────────

fn main() {
    if !std::env::args().any(|arg| arg == "--bench") {
        return;
    }
    let prereg = isomesh::experiment!("P-174");

    common::experiment::run(prereg, |run| {
        let fixtures = [
            Fixture::Sphere(Sphere::<f64>::canonical()),
            Fixture::Torus(Torus::<f64>::canonical()),
            Fixture::Capsule(Capsule {
                radius: CAPSULE_RADIUS,
                half_length: CAPSULE_HALF_LENGTH,
            }),
        ];

        let mut cells: Vec<Cell> = Vec::with_capacity(fixtures.len() * LADDER.len());
        for fixture in &fixtures {
            for &samples in &LADDER {
                let cell = measure(fixture, samples);
                println!(
                    "  {:8} {:>4}³  V={:6} T={:6} eps={:.4} eps/h={:.2} sup={:5.0} | \
                     nc H={:.4e} K={:.4e} {:7.2}ms | var H={:.4e} K={:.4e} {:7.2}ms",
                    cell.field,
                    cell.samples,
                    cell.vertices,
                    cell.triangles,
                    cell.eps,
                    cell.eps / cell.cell_size,
                    cell.support_mean,
                    cell.nc.mean_error,
                    cell.nc.gaussian_error,
                    cell.nc.clock.median,
                    cell.varifold.mean_error,
                    cell.varifold.gaussian_error,
                    cell.varifold.clock.median,
                );
                cells.push(cell);
            }
        }

        // ── vacuity controls ────────────────────────────────────────────────

        // The ladder is one regime: eps falling, eps/h rising, one domain.
        for fixture in &fixtures {
            let (flo, fhi) = fixture.domain();
            assert!(
                (flo[0] + DOMAIN_HALF_SPAN).abs() < 1e-12
                    && (fhi[0] - DOMAIN_HALF_SPAN).abs() < 1e-12,
                "VOID: {} is sampled on [{}, {}] but the other fixtures use \
                 [-{DOMAIN_HALF_SPAN}, {DOMAIN_HALF_SPAN}], so `cell_size` does not mean \
                 the same thing across the ladder and the fitted exponent splices two \
                 regimes",
                fixture.name(),
                flo[0],
                fhi[0]
            );
        }
        for window in cells.as_chunks::<RUNGS>().0 {
            for pair in window.windows(2) {
                let (a, b) = (&pair[0], &pair[1]);
                assert!(
                    b.eps < a.eps && b.eps / b.cell_size > a.eps / a.cell_size,
                    "VOID: {} goes from {}³ (eps={:.6}, eps/h={:.4}) to {}³ \
                     (eps={:.6}, eps/h={:.4}); Buet-Leonardi-Masnou's bound needs eps -> 0 \
                     with h/eps -> 0, so a rung that breaks either monotonicity makes \
                     `convergence_exponent` a fit across two regimes",
                    a.field,
                    a.samples,
                    a.eps,
                    a.eps / a.cell_size,
                    b.samples,
                    b.eps,
                    b.eps / b.cell_size
                );
            }
        }

        for cell in &cells {
            // A closed manifold, or both arms are calibrated against the wrong
            // reference.
            assert!(
                cell.boundary_edges == 0
                    && cell.non_manifold_edges == 0
                    && cell.zero_length_sides == 0,
                "VOID: {} at {}³ has {} boundary edges, {} non-manifold edges and {} \
                 zero-length sides; a boundary vertex's defect is `pi - alpha_v` not \
                 `2pi - alpha_v` and a zero-length side costs the per-face `Sum alpha = pi` \
                 identity, so the incumbent would be scored against the wrong analytic \
                 reference",
                cell.field,
                cell.samples,
                cell.boundary_edges,
                cell.non_manifold_edges,
                cell.zero_length_sides
            );

            // The reference is not a constant zero, and no error is zero or NaN.
            assert!(
                cell.h_exact_rms > 0.0 && cell.k_exact_rms > 0.0,
                "VOID: {} at {}³ has h_exact_rms={} and k_exact_rms={}; a relative error \
                 normalised by zero is not a measurement (M-44)",
                cell.field,
                cell.samples,
                cell.h_exact_rms,
                cell.k_exact_rms
            );
            for (arm, label) in [(&cell.nc, "normal_cycles"), (&cell.varifold, "varifold")] {
                assert!(
                    arm.mean_error.is_finite()
                        && arm.mean_error > 0.0
                        && arm.gaussian_error.is_finite()
                        && arm.gaussian_error > 0.0,
                    "VOID: {label} on {} at {}³ reports mean_error={} and \
                     gaussian_error={}; `convergence_exponent` takes the log of these, and \
                     an error that is exactly zero on a discretised mesh is a broken \
                     fixture rather than a perfect estimator (M-44)",
                    cell.field,
                    cell.samples,
                    arm.mean_error,
                    arm.gaussian_error
                );
            }

            // A convolution, not a handful of terms.
            assert!(
                cell.support_mean >= MIN_SUPPORT,
                "VOID: {} at {}³ averages only {:.2} triangles under the mollifier \
                 (eps={:.6}, h={:.6}); below {MIN_SUPPORT} `Sum_T area(T) rho_eps` is a \
                 small sum wearing a convolution's name and the varifold arm is being \
                 scored on its neighbour count",
                cell.field,
                cell.samples,
                cell.support_mean,
                cell.eps,
                cell.cell_size
            );

            // The registered control, half 1: the incumbent reproduces the
            // analytic curvature globally, on sphere and torus.
            if cell.field != "capsule" {
                let rel = (cell.global_int_h - cell.global_int_h_exact).abs()
                    / cell.global_int_h_exact.abs();
                assert!(
                    rel <= GLOBAL_MEAN_TOLERANCE,
                    "VOID: the normal-cycle mean measure on {} at {}³ gives \
                     ½ Sum l beta = {:.9} against the analytic integral H dA = {:.9}, a \
                     relative error of {rel:.6} above {GLOBAL_MEAN_TOLERANCE}; the \
                     incumbent is not reproducing the analytic curvature, so C1 would \
                     compare a calibrated arm against an uncalibrated one",
                    cell.field,
                    cell.samples,
                    cell.global_int_h,
                    cell.global_int_h_exact
                );
                assert!(
                    (cell.global_defect - cell.global_defect_exact).abs() <= GAUSS_BONNET_TOLERANCE,
                    "VOID: the normal-cycle Gaussian measure on {} at {}³ sums to {:.12} \
                     against 2 pi chi = {:.12}; that sum is an identity (M-340), so a \
                     miss means the angle accumulator is wrong and its localisation \
                     cannot be scored",
                    cell.field,
                    cell.samples,
                    cell.global_defect,
                    cell.global_defect_exact
                );
            }
        }

        // The registered control, half 2: both arms reproduce it pointwise, at
        // the finest rung, on sphere and torus.
        let finest = LADDER[LADDER.len() - 1];
        for cell in cells
            .iter()
            .filter(|c| c.samples == finest && c.field != "capsule")
        {
            for (arm, label) in [(&cell.nc, "normal_cycles"), (&cell.varifold, "varifold")] {
                assert!(
                    arm.mean_error <= CALIBRATION_CEILING,
                    "VOID: {label} on {} at {}³ has an area-weighted relative RMS \
                     mean-curvature error of {:.6}, above {CALIBRATION_CEILING}; an \
                     estimator that far off the analytic curvature is not calibrated, and \
                     `error_ratio` between two such numbers is a ratio of two noises",
                    cell.field,
                    cell.samples,
                    arm.mean_error
                );
            }
        }

        // The torus really changes the sign of K, and the capsule really carries
        // a zero-K region.
        for cell in &cells {
            match cell.field {
                "torus" => assert!(
                    cell.k_positive > 0 && cell.k_negative > 0,
                    "VOID: the torus at {}³ scores {} vertices with K > 0 and {} with \
                     K < 0; the one fixture that is supposed to carry varying, \
                     sign-changing Gaussian curvature is behaving like a sphere, and C1 \
                     would be decided on constant-curvature data alone",
                    cell.samples,
                    cell.k_positive,
                    cell.k_negative
                ),
                "capsule" => assert!(
                    cell.k_zero_exact > 0 && cell.k_positive > 0,
                    "VOID: the capsule at {}³ scores {} vertices with K == 0 exactly and \
                     {} with K > 0; the only fixture that tests whether an estimator \
                     manufactures Gaussian curvature out of mesh noise has an unpopulated \
                     zero-curvature region",
                    cell.samples,
                    cell.k_zero_exact,
                    cell.k_positive
                ),
                _ => {}
            }
        }

        // ── rows ────────────────────────────────────────────────────────────

        // `cells` was filled field-major, so each chunk of `RUNGS` is one
        // field's whole ladder and the fits are computed exactly once.
        let windows = cells.as_chunks::<RUNGS>().0;
        let fits: Vec<FieldFit> = windows.iter().map(FieldFit::of).collect();
        let c1_global = cells
            .iter()
            .all(|c| c.varifold.mean_error / c.nc.mean_error >= 1.0);
        let c2_global = fits.iter().all(|f| f.c2);

        for (window, fit_set) in windows.iter().zip(&fits) {
            for cell in window {
                let error_ratio = cell.varifold.mean_error / cell.nc.mean_error;
                let cost_ratio = cell.varifold.clock.median / cell.nc.clock.median;
                let c1 = error_ratio >= 1.0;

                for (estimator, arm, fit, gauss_fit) in [
                    ("normal_cycles", &cell.nc, &fit_set.nc, &fit_set.nc_gaussian),
                    (
                        "varifold",
                        &cell.varifold,
                        &fit_set.varifold,
                        &fit_set.varifold_gaussian,
                    ),
                ] {
                    let total_ms = cell.extract_ms + arm.clock.median;
                    run.record(&[
                        // ── registered (11), in registration order ──
                        ("estimator", estimator.to_string()),
                        ("field", cell.field.to_string()),
                        ("resolution", cell.samples.to_string()),
                        ("mean_curvature_error", format!("{:.6e}", arm.mean_error)),
                        (
                            "gaussian_curvature_error",
                            format!("{:.6e}", arm.gaussian_error),
                        ),
                        ("error_ratio", format!("{error_ratio:.6}")),
                        ("estimator_ms", format!("{:.4}", arm.clock.median)),
                        ("cost_ratio", format!("{cost_ratio:.6}")),
                        ("convergence_exponent", format!("{:.6}", fit.exponent)),
                        ("c1_holds", c1.to_string()),
                        ("c2_holds", fit_set.c2.to_string()),
                        // ── extras (M-273) ──
                        ("boundary_edges", cell.boundary_edges.to_string()),
                        ("c1_holds_global", c1_global.to_string()),
                        ("c2_holds_global", c2_global.to_string()),
                        ("cell_size", format!("{:.9}", cell.cell_size)),
                        (
                            "cost_matched",
                            (cost_ratio >= COST_BAND.recip() && cost_ratio <= COST_BAND)
                                .to_string(),
                        ),
                        (
                            "curvature_share",
                            format!("{:.6}", arm.clock.median / total_ms),
                        ),
                        ("eps", format!("{:.9}", cell.eps)),
                        ("eps_over_h", format!("{:.6}", cell.eps / cell.cell_size)),
                        ("estimator_ms_max", format!("{:.4}", arm.clock.max)),
                        ("estimator_ms_min", format!("{:.4}", arm.clock.min)),
                        (
                            "estimator_ms_scatter",
                            format!("{:.6}", arm.clock.max / arm.clock.min),
                        ),
                        ("exponent_gap", format!("{:.6}", fit_set.gap)),
                        (
                            "exponent_gap_sigma",
                            format!("{:.6}", fit_set.gap / fit_set.sigma),
                        ),
                        ("exponent_stderr", format!("{:.6}", fit.stderr)),
                        ("extract_ms", format!("{:.4}", cell.extract_ms)),
                        ("gaussian_exponent", format!("{:.6}", gauss_fit.exponent)),
                        (
                            "global_defect_exact",
                            format!("{:.9}", cell.global_defect_exact),
                        ),
                        (
                            "global_defect_measured",
                            format!("{:.9}", cell.global_defect),
                        ),
                        (
                            "global_int_h_exact",
                            format!("{:.9}", cell.global_int_h_exact),
                        ),
                        ("global_int_h_measured", format!("{:.9}", cell.global_int_h)),
                        (
                            "global_int_h_rel_error",
                            format!(
                                "{:.6e}",
                                (cell.global_int_h - cell.global_int_h_exact).abs()
                                    / cell.global_int_h_exact.abs()
                            ),
                        ),
                        ("h_exact_rms", format!("{:.6}", cell.h_exact_rms)),
                        ("k_exact_rms", format!("{:.6}", cell.k_exact_rms)),
                        ("k_negative_vertices", cell.k_negative.to_string()),
                        ("k_positive_vertices", cell.k_positive.to_string()),
                        ("k_zero_exact_vertices", cell.k_zero_exact.to_string()),
                        ("mesh_area", format!("{:.6}", cell.mesh_area)),
                        (
                            "nc_star_gaussian_error",
                            format!("{:.6e}", cell.nc_star_gaussian_error),
                        ),
                        (
                            "nc_star_mean_error",
                            format!("{:.6e}", cell.nc_star_mean_error),
                        ),
                        ("non_manifold_edges", cell.non_manifold_edges.to_string()),
                        ("normal_alignment", format!("{:.9}", cell.normal_alignment)),
                        ("repeats", REPEATS.to_string()),
                        (
                            "smooth_only_mean_error",
                            format!("{:.6e}", arm.smooth_only_mean_error),
                        ),
                        (
                            "support_triangles_mean",
                            format!("{:.3}", cell.support_mean),
                        ),
                        ("triangles", cell.triangles.to_string()),
                        ("used_vertices", cell.used_vertices.to_string()),
                        (
                            "varifold_h_shape_error",
                            format!("{:.6e}", cell.varifold_h_shape_error),
                        ),
                        ("vertices", cell.vertices.to_string()),
                        ("zero_length_sides", cell.zero_length_sides.to_string()),
                    ]);
                }
            }
        }
    });
}
