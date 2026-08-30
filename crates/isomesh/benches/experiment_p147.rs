//! **P-147 — the rate does not improve and the constant does, stated before the
//! measurement rather than discovered in it.**
//!
//! Ticket: R-147. Pre-registered before this harness existed.
//!
//! ```bash
//! cargo bench --bench experiment_p147
//! ```
//!
//! Writes `docs/experiments/p-147.csv`.
//!
//! # What was missing
//!
//! `M-12` (FINDINGS.md:1128) measured Marching Cubes' error falling like `h²` —
//! mean error `2.7168e-3` at `32³` against `6.5015e-4` at `64³`, a ratio of
//! `4.179` against the ideal `4.13` — and stopped there. `M-113` then found the
//! *fitted constant* does not survive across configurations, and `P-155` (wave 1,
//! `docs/experiments/p-155.csv`) derived the order from Strang–Fix and failed to
//! predict the constant on **0 of 8** fields against a bar of four. So the repo
//! holds a law of the form `E = C · N^(−α)` in which `α` is understood and `C` is
//! not, and it has never once **decomposed** a change to the mesh into its effect
//! on `α` versus its effect on `C`.
//!
//! That decomposition is this row, and the registration states its answer in
//! advance rather than discovering it: Bonito, Canuto, Nochetto & Veeser (Acta
//! Numerica 2024, corpus `10.1017/s0962492924000011`) say flatly that the order
//! *"cannot be improved upon assuming either higher regularity … or a graded
//! mesh"*. Uniform refinement gives `O(N^(−2/3))` in 3-D and the optimal
//! anisotropic mesh gives `O(N^(−2/3))`: the **same exponent**. What anisotropy
//! buys is the constant, `‖√|det H|‖_{L^τ}` in place of `|f|_{W^{2,p}}`, and by
//! AM–GM the former is never larger.
//!
//! # What P-146 actually measured, quoted rather than assumed
//!
//! This row consumes `common::metric` and reproduces R-146's arm. **R-146's own
//! clauses did not survive**, and that is quoted here from
//! `docs/experiments/p-146.csv` rather than paraphrased, because a consumer that
//! assumes a win it did not get is reading its own hopes:
//!
//! - **C1 FALSIFIED.** `c1_holds=false` on all 40 rows; `c1_winners=0` against a
//!   bar of `C1_MIN_WINNERS=3`, over a `c1_population` of **4** — only
//!   `FieldBound::Exact` was admitted to it. Of those four the metric-driven arm
//!   *tied* on `sphere` and `box_exact` (`ratio=1.000000`) and was **worse** on
//!   `torus` (`ratio=1.153211`) and `thin_plate` (`ratio=2.298935`, i.e. 2.3×
//!   *more* triangles at matched error).
//! - **C3 FALSIFIED.** `metric_share` ran from `0.332862` (`thin_plate` at 65³)
//!   to `8.747515` (`csg_difference` at 17³) — between 33% and 875% of
//!   extraction, never once under the 15% bar.
//! - **C2 reads `unmeasurable` on 20 of 40 rows**, because `validate::accuracy`
//!   was gated on `field.bound().is_exact()` and four of the eight reference
//!   fields are `Lipschitz`, `Underestimate` or `Unbounded`.
//!
//! **The one number this row is built on** is `p-146.csv`'s `axis_ratio`, and it
//! is the reason every prediction below reads the way it does. The per-axis
//! metric-driven grid came out *identical to the isotropic grid* on **five of
//! eight** fields — `axis_ratio=1.000000` at every rung on `sphere`,
//! `box_exact`, `csg_difference`, `gyroid` and `noise_cavity` — mildly graded on
//! `torus` (`1.095`–`1.174`), and strongly graded on only two: `thin_plate`
//! (`7.67`–`30.24`) and `fbm_terrain` (`6.20`–`47.40`).
//!
//! So R-147 is not a rerun of R-146's comparison. R-146 asked *"how many
//! triangles at matched error"* and got no win; this asks *"what does the error
//! law's exponent do, and separately what does its constant do"* — a question
//! whose answer is informative even where the arms tie, because a tie is
//! `Δexponent = 0` exactly and that is C1's own prediction arriving by a route
//! nobody wanted.
//!
//! # The construction, reproduced from `experiment_p146.rs` and checked
//!
//! The registration compares two *arms*, and the two rows are commensurable only
//! if the arms are the same object. The two benches cannot share bench-local
//! code, so `experiment_p146.rs`'s construction is **restated here verbatim** —
//! [`band_points`] (`p146:413`), [`round_odd`] (`:587`), [`anisotropic_grid`]
//! (`:620`) with its lower-clamp-and-re-solve loop and no upper clamp, and
//! [`Stretched`] (`:676`) with the same `mul_add` map back to world space. The
//! ladder is P-146's own `RUNGS = [17, 25, 33, 49, 65]` and not `P-155`'s
//! `[19, 27, 35, 47, 63]`, for the same reason.
//!
//! **The reproduction is checked rather than asserted in prose.**
//! [`read_p146_baseline`] parses the committed `docs/experiments/p-146.csv` and
//! [`Reproduction`] compares, per field and per rung:
//!
//! | column | compared against | rows |
//! |---|---|---|
//! | `triangles_isotropic` | this file's isotropic triangle count | 40 |
//! | `triangles_anisotropic` | this file's anisotropic triangle count | 40 |
//! | `grid_anisotropic` | `format!("{nx}x{ny}x{nz}")` from [`anisotropic_grid`] | 40 |
//! | `hausdorff_isotropic` | `AccuracyReport::symmetric_hausdorff` at `{:.9}` | 20 |
//! | `hausdorff_anisotropic` | the same, anisotropic arm | 20 |
//!
//! A mismatch aborts and names the field, the rung and both values. Nothing is
//! hard-coded: if `p-146.csv` is ever re-measured the check re-reads it.
//! `hausdorff_*` is compared only on the 20 rows where P-146 recorded a number
//! rather than `unmeasurable:bound=…`. Triangle counts and grid shapes are exact
//! integers and cover all 40 rows, so the construction is pinned even where the
//! instrument was not.
//!
//! # Arms
//!
//! | arm | what it varies | is_control |
//! |---|---|---|
//! | **isotropic** | `(N, N, N)` samples, `h` equal on every axis | **yes** — this is uniform refinement, the exponent `2/3` belongs to it |
//! | **anisotropic** | `(n_x, n_y, n_z)` from the metric's per-axis point densities, `∏ n_a ≈ N³` | no |
//! | **the ladder** | five rungs, `17 → 65` samples/axis, a factor of `4` in `h` | — |
//! | **the fit** | ordinary least squares of `ln E` on `ln(dof / dof_ref)` | — |
//!
//! Both arms are extracted by one `MarchingCubes::new()` at its shipped
//! defaults and graded by **one instrument on one seed lattice per rung** — the
//! *isotropic* rung's `shape`, `origin` and `cell_size` in both calls, which
//! `validate/accuracy.rs:332-337` explicitly licenses. An arm that grades its
//! own homework on its own lattice is not a comparison.
//!
//! # The instrument, and why it reaches all eight fields where P-146 reached four
//!
//! `E` is **`AccuracyReport::field_to_mesh.max`**: every point of the seed
//! lattice is Newton-projected onto the zero set along `∇f`, and the Euclidean
//! distance from that projected point to the nearest mesh triangle is measured.
//! It is a one-sided Hausdorff distance, `field → mesh`, sampled on a lattice
//! and therefore a **lower bound** on the true one.
//!
//! P-146 read `symmetric_hausdorff()` and gated it on `field.bound().is_exact()`,
//! which is correct for P-146's clause and costs it four fields. Reading
//! `validate/accuracy.rs:565-603`, **`project` compares nothing against `|f|`**:
//! it is Newton on `f` along `∇f` and it returns a Euclidean length. Every one
//! of the eight reference fields ships an *analytic* `gradient` — `Sphere` at
//! `fields/mod.rs:309`, `Torus` `:395`, `BoxExact` `:545`, `ThinPlate` `:655`,
//! `Difference` `:716`, `Intersection` `:750`, `Gyroid` `:1028`, `NoiseVolume`
//! `:1211`, `FbmTerrain` `:1365` — so the projection is exact-gradient on all
//! eight and the measured distances are geometric on all eight.
//!
//! What a non-`Exact` bound actually costs is **seed coverage, not distance**.
//! `project`'s `max_first_step` band test is free precisely because for a true
//! distance field the first Newton step *is* the distance
//! (`accuracy.rs:562-564`); for `Lipschitz`, `Underestimate` and `Unbounded`
//! fields it is not, so `band_radius` admits and rejects slightly the wrong
//! seeds. That biases *which* points are measured. It cannot bias the fitted
//! exponent, because `band_radius = h · BAND_RADIUS_REL` scales with `h`: the
//! seed band is the same number of cells wide at every rung, so the seed
//! population is geometrically similar along the whole ladder. And it cannot
//! bias `constant_ratio`, because both arms are seeded from **one** lattice per
//! rung and the bias is common to the numerator and the denominator. Coverage is
//! recorded per rung and per arm in `seed_coverage_series_*` so a reader can see
//! it rather than take that argument on trust.
//!
//! Two more instruments are fitted and recorded beside it, and neither decides a
//! clause: `field_to_mesh.mean` (`exponent_mean_*`), which says whether an
//! `L^∞` verdict is an artefact of a single worst point, and `symmetric_hausdorff`
//! (`hausdorff_symmetric_series_*`), which is the column P-146's numbers are
//! reproduced against.
//!
//! **Why not `max |f(v)|`, the residual P-146 also recorded.** Because it is
//! *identically zero* on the two grid-aligned polyhedra and therefore cannot be
//! fitted at all. `box_exact` is `[-1, 1]³` inside `[-2, 2]³`, and `x = ±1` is a
//! sample plane at every rung of this ladder (`17, 25, 33, 49, 65` all put a
//! sample at `+1`), so every corner on a face has `f == 0` exactly, every
//! crossing interpolates to `t = 0`, and every vertex lands on the surface:
//! `p-146.csv` reads `mesh_residual_max_iso = 0.000000e0` on all five
//! `box_exact` rungs and `~1e-17` on all five `thin_plate` rungs. The error
//! those meshes make is **tangential** — the box's twelve edges are chamfered —
//! and a residual measured along the normal is structurally blind to it. The
//! reverse direction sees it: P-146 measured `box_exact`'s symmetric Hausdorff
//! at exactly `1.154701 · h` on every rung (`0.288675135` at `h = 0.25` through
//! `0.072168784` at `h = 0.0625`), which is **first order**, and that is the
//! `W^{2,p}` failure this row is trying to see.
//!
//! # The abscissa, which is the one methodological choice that matters
//!
//! The registration's exponent is `N^(−2/3)` *in three dimensions*, so its `N`
//! is a **3-D degree-of-freedom count**, not a triangle count: uniform
//! refinement has `E ∝ h²` and `dof ∝ h^(−3)`, giving `E ∝ dof^(−2/3)` exactly
//! as quoted. The abscissa is therefore `dof = n_x · n_y · n_z`, the arm's own
//! sample count — which P-146's construction budget-matches between the arms by
//! design, so both arms are fitted over the same ladder of the same quantity.
//!
//! Against a **triangle** count the same data gives different exponents, and
//! that is not a defect but the second half of the answer, so it is recorded as
//! `exponent_vs_triangles_*` and `exponent_difference_vs_triangles`. A surface
//! mesh has `T ∝ h^(−2)`, so `E ∝ T^(−1)` under uniform refinement and the
//! exponent against `T` is `1` where the exponent against `dof` is `2/3`. A
//! reader comparing this row to a paper must check which `N` the paper meant.
//!
//! # The constant, and the units trap in it
//!
//! `fitted_constant_*` is **not** the raw fit intercept. Writing
//! `E = C · dof^(−α)`, the prefactor `C` carries units of `length · dof^α`, so
//! when the two arms' `α` differ — which is exactly the case C1 is registered to
//! test — `C_iso` and `C_aniso` are quantities of different dimension and their
//! ratio is not a number. That would make `constant_ratio` unreadable in the one
//! case the row exists to examine.
//!
//! So the fit is taken in the normalised abscissa `x = ln(dof / dof_ref)` with
//! `dof_ref = 65³ = 274625`, the isotropic ladder's finest rung, and
//! `fitted_constant_*` is the fit's value at `x = 0`: **the fitted error, in
//! world length units, at the reference sample budget**. Its ratio is
//! dimensionless at any pair of exponents, and where the exponents agree it is
//! proportional to the classical prefactor, so nothing is lost. The raw
//! prefactors are recorded anyway as `prefactor_*` so the fit can be
//! reconstructed without knowing `dof_ref`.
//!
//! This makes `constant_ratio` the exact dual of P-146's `ratio`: P-146 read the
//! *triangle count at matched error*, this reads the *error at matched budget*.
//!
//! # Which fields are smooth, decided before the run
//!
//! C1 is quantified "on every **smooth** field", so the population is fixed
//! here, by one criterion — **is `f` continuously differentiable on its own zero
//! set** — and not by a number the run produces. `C¹` on the zero set is the
//! property that makes linear interpolation along a grid edge second order; a
//! gradient jump *on the surface* is what drops it to first.
//!
//! | field | in C1's population | why |
//! |---|---|---|
//! | `sphere` | **yes** | `‖x‖ − r` is `C^∞` away from the origin, which is not on the surface |
//! | `torus` | **yes** | `C^∞` away from the core circle, which is not on the surface |
//! | `box_exact` | no | a polyhedron: `∇f` jumps along all twelve edges and at all eight corners, and those are *on* the zero set |
//! | `csg_difference` | no | `max(f_box, −f_sphere)`, and `fields/mod.rs:706-725` documents the gradient as *"discontinuous at the seam"* — **C3's subject** |
//! | `thin_plate` | no | a slab is a box; its four rim edges are creases on the zero set |
//! | `gyroid` | no | `capped_gyroid` is `max(gyroid, sphere_6)` over `[-7, 7]³`, and the gyroid's nodal surface is space-filling, so the cap circle is a crease on the zero set |
//! | `fbm_terrain` | **yes** | `p_y − amp · fbm(x, 0, z)`; Perlin's fade is the quintic `t³(6t² − 15t + 10)` with `u'' = 60t(t−1)(2t−1)` vanishing at both ends (`fields/noise.rs:85-91`), so the noise is `C²` across every lattice cell wall and `∇f` never vanishes — the zero set is a `C²` graph |
//! | `noise_cavity` | no | `max(noise, sphere_1.5)`: the cap circle is a crease |
//!
//! So **C1's population is three**, and `c1_holds` is the conjunction over those
//! three. `c1_smooth` and `c1_row_holds` put each field's own membership and own
//! verdict on its own row.
//!
//! # Predicted verdicts, with the arithmetic, before the harness ran
//!
//! **C1 — predicted FALSIFIED, by `fbm_terrain`, and the falsifier calls that
//! "the more interesting result".** Two of the three smooth fields cannot
//! disagree: `sphere`'s arms are the *same grid* (`p-146.csv` `axis_ratio=1.0`
//! at every rung), so `Δexponent = 0` exactly; `torus`'s grid moves by at most
//! `1.174:1`, a sub-rung perturbation. `fbm_terrain` is the one that can, and
//! the prediction is that it does, for a reason that does not contradict the
//! cited theory but reads its `d`:
//!
//! > `f = p_y − amp·fbm(x, 0, z)` is **exactly affine in `y`**, so every second
//! > difference touching `y` cancels identically and the crossing on a `y`-edge
//! > is placed *exactly* at any spacing. P-146's grid pins `n_y` at
//! > `MIN_SAMPLES = 5` and spends the whole budget laterally
//! > (`grid_anisotropic` runs `31x5x31` → `237x5x233`), so
//! > `h_xz ∝ dof^(−1/2)` and `E ∝ h_xz² ∝ dof^(−1)`, against the isotropic
//! > `h ∝ dof^(−1/3)` and `E ∝ dof^(−2/3)`. Predicted
//! > `Δexponent ≈ 1 − 2/3 = 0.333`, an order of magnitude over the `0.1` bar.
//!
//! The theory is not wrong: `N^(−2/d)` with a genuinely degenerate direction is
//! `N^(−2/2)`, not `N^(−2/3)`. What the registration's clause assumes is that
//! `d = 3` for every field in the roster, and one field's Hessian has an
//! **exact** null direction rather than a small one. Corroboration from wave 1:
//! `p-146.csv`'s `mesh_residual_max_*` for `fbm_terrain`, converted back to
//! `|f|` by multiplying by each rung's `h`, gives `0.2106 → 0.00675` on the
//! anisotropic arm against `0.331 → 0.0728` on the isotropic one — a log-log
//! slope of `0.85` against `0.38` over the same `dof` ladder.
//!
//! **C2 — predicted FALSIFIED, on its first conjunct, arithmetically.**
//! `constant_ratio` is exactly `1.000000` on the five fields whose two arms are
//! the same grid, because `Stretched` with `s = [1, 1, 1]` samples
//! `f(q · 1 + lo)` on the same lattice the isotropic arm samples and emits a
//! bit-identical mesh (`p-146.csv` confirms it: identical triangle counts, and
//! `mesh_residual_max_aniso == mesh_residual_max_iso` to the last digit on
//! `sphere`, `box_exact` and `csg_difference`). Of the remaining three,
//! `thin_plate`'s lateral resolution *collapses* to `9`–`21` samples and P-146
//! already measured its Hausdorff getting **worse** by `2.3×`, and `torus`'s got
//! worse on four of five rungs. So at most one field of eight improves, against
//! a bar of a strict majority. The second conjunct then fails too: a Spearman
//! correlation over eight points of which five are tied at exactly `1.0` is
//! computed against a five-fold tie block, and ties averaged cannot produce
//! `0.7`. Recorded as `c2_improved_fields` and `c2_rank_correlation`, with
//! `τ = 2`- and `τ = ∞`-norm gaps beside the registered one so the verdict's
//! dependence on the one free parameter is visible rather than assumed away.
//!
//! **C3 — predicted FALSIFIED, and the reason is not the one its falsifier
//! offers.** The falsifier reads *"C3 by no exponent difference on the CSG field,
//! which would say our sharp fields are smoother than assumed"*. There will be
//! no exponent difference, but `csg_difference` is not smooth: `BoxExact`'s
//! half-extents are `[1, 1, 1]` and the subtracted sphere sits on the
//! `(0.6, 0.6, 0.6)` diagonal (`fields/mod.rs:917-923`), so the field is
//! **permutation-symmetric in `(x, y, z)`** — its three per-axis metric weights
//! are equal to rounding, `anisotropic_grid` returns `[N, N, N]`, and the two
//! arms are the same grid. `p-146.csv` measured exactly that:
//! `axis_ratio=1.000000` and `axes_pinned=0` on all five `csg_difference` rungs.
//!
//! `exponent_difference` is therefore predicted **identically zero**, and the
//! honest reading is that a *per-axis global* grid has three degrees of freedom
//! and a permutation-symmetric field annihilates all three. This is reported as
//! `c3_holds = false` — the registration's falsifier is explicit and is
//! honoured rather than reinterpreted — with `c3_reason` and
//! `c3_arms_identical` naming the arithmetic on every row, so nobody reads the
//! `false` as evidence that the CSG field is smooth. A per-cell construction is
//! a different row and needs a source change to ask for.
//!
//! Because a `false` reached that way carries little, one **companion** number
//! is recorded beside it and is explicitly *not* C3's verdict:
//! `companion_regularity_deficit` is the median isotropic exponent over C1's
//! smooth population minus `csg_difference`'s own isotropic exponent. That is
//! the observable "lacks `W^{2,p}` regularity" actually predicts, it is
//! measurable, and it is predicted **positive and large**: `P-155` fitted
//! `1.985` in `h` on `sphere` and exactly `1.000` on `box_exact`
//! (`p-155.csv`'s `fitted_exponent`), i.e. `0.662` against `0.333` in `dof`, so
//! a deficit near `0.33` is expected on the field that shares `box_exact`'s
//! creases.
//!
//! # SHARE, recomputed before the numbers
//!
//! The registration's SHARE is *"none — this is `M-12`'s law, decomposed"*, and
//! that is discharged rather than restated. This row proposes **no landing** and
//! moves no stage of the frame: it changes no extractor, adds no field, and its
//! anisotropic arm is a bench-local sampling grid that the shipped `extract` —
//! which takes one **scalar** `cell_size` at `marching_cubes/mod.rs:193` and six
//! other sites — cannot express. What it produces is a statement about the shape
//! of `M-12`'s law: which of `α` and `C` a mesh change is allowed to move. The
//! consumer-facing value is entirely in that, and the number that would justify
//! a landing is a `constant_ratio` below one on a field whose arms genuinely
//! differ. On the prediction above there is exactly one such field, so the
//! realised share is zero and Phase 28 is where any landing would be registered.
//!
//! # Vacuity controls
//!
//! Every one runs before the first `run.record`, and every panic starts
//! `"VOID: "`. `M-44`: a zero that could not have been non-zero is not a
//! measurement.
//!
//! - **The registered control** — the `am_gm_gap` column must vary by at least
//!   `2×` across fields, or C2's correlation is against a constant. Asserted as
//!   `max > 2 · min` so that a `min` of exactly zero still requires a non-zero
//!   `max`. Columns: `am_gm_gap`, `am_gm_gap_spread`.
//! - **The gap must be non-zero somewhere**, or the numerator `‖√|det H|‖` is
//!   zero on every field and the ratio is `0/x` eight times over. Column:
//!   `am_gm_gap_spread`.
//! - **The anisotropic arm must be anisotropic somewhere.** `axis_ratio_max`
//!   over the whole sweep must exceed `1.5`, or `exponent_difference` is zero on
//!   every row by construction and C1, C2 and C3 are all measuring the
//!   extractor's determinism. This is P-146's own control, carried because C3's
//!   verdict is a zero and a zero needs the chance to have been non-zero.
//!   Column: `axis_ratio_max`.
//! - **C1's population must contain a field whose arms differ.** Without one,
//!   `c1_holds` is `true` for the same reason `0 < 0.1` is true and says nothing
//!   about the exponent. Columns: `c1_smooth`, `arms_identical`.
//! - **The ladder must move the instrument** on every field a clause reads —
//!   C1's three smooth fields and C3's `csg_difference` — by at least `1.2×` on
//!   both arms, or a slope is being read off a horizontal line. Columns:
//!   `ladder_span_isotropic`, `ladder_span_anisotropic`.
//! - **Every rung of every arm must have produced a mesh, a seed and a positive
//!   error.** `ln 0` is not a data point. Columns: `triangles_series_*`,
//!   `seed_coverage_series_*`, `error_max_series_*`.
//! - **The band must be non-empty** on every rung, or the metric, the gap and
//!   the flat-direction fraction are statistics of nothing. Column:
//!   `band_points_series`.
//!
//! Beside them, and not a vacuity control but a correctness gate, is the
//! reproduction check against `p-146.csv` described above: 120 exact comparisons
//! on 40 rows, recorded as `reproduction_checks` and `reproduction_mismatches`.
//!
//! # Determinism
//!
//! One thread, no PRNG, no map iteration, `f64` throughout. Sorting is
//! [`f64::total_cmp`]. The band is swept `z`, `y`, `x` with `x` innermost, the
//! crate's order. The only seeded object anywhere near this row is
//! `FbmTerrain`'s and `NoiseVolume`'s committed `0x5EED_1234`, which belongs to
//! the fields. `wall_seconds` is the only machine-dependent column and no clause
//! reads it — the registration names no cost threshold, so this row records
//! counts and distances and never a ratio of two clocks.

#![allow(
    clippy::cast_possible_truncation,
    clippy::cast_precision_loss,
    clippy::cast_sign_loss,
    clippy::too_many_lines
)]

mod common;

use std::path::PathBuf;
use std::time::Instant;

use isomesh::fields::{FieldBound, ReferenceField};
use isomesh::marching_cubes::MarchingCubes;
use isomesh::validate::{AccuracyConfig, accuracy};
use isomesh::{MeshBuffer, RuntimeShape3, Sdf};

use crate::common::beta::rank_correlation;
use crate::common::metric::{H_FLOOR, Sym3, am_gm_gap, hessian, metric_lp};

/// Samples per axis for the isotropic arm.
///
/// **P-146's ladder, unchanged** (`experiment_p146.rs:275`), because this row's
/// whole claim to commensurability is that it runs P-146's arm. Five rungs
/// spanning `17 → 65`, a factor of `4` in `h` and so roughly `16×` in a
/// second-order error, all odd (`M-266`). Five is one more than the registered
/// floor of four, and the extra rung is what buys a slope confidence interval
/// with three degrees of freedom instead of two.
const RUNGS: [u32; 5] = [17, 25, 33, 49, 65];

/// The reference sample budget the fitted constants are read off at: the
/// isotropic ladder's finest rung, `65³`.
///
/// Derived from [`RUNGS`] rather than written down, so the two cannot drift.
const REFERENCE_DOF: u32 = {
    let n = RUNGS[RUNGS.len() - 1];
    n * n * n
};

/// Student's `t` at 97.5% with `RUNGS.len() − 2 = 3` degrees of freedom.
///
/// Hard-coded because `std` has no `t` distribution and this row needs exactly
/// one quantile. [`LADDER_RUNGS_ARE_FIVE`] is the assertion that keeps the two
/// consistent: change [`RUNGS`]'s length and the build stops here.
const T_95_DF3: f64 = 3.182_446_305_284_63;

/// A compile-time guard on [`T_95_DF3`]'s degrees of freedom.
const LADDER_RUNGS_ARE_FIVE: () = assert!(
    RUNGS.len() == 5,
    "T_95_DF3 is the 97.5% quantile for 3 degrees of freedom, i.e. exactly five rungs"
);

/// The norm `M_Lp` is optimised for. P-146's value, unchanged; `R-150` owns the
/// norm sweep.
const P_NORM: f64 = 2.0;

/// The `l^τ` exponent the AM–GM gap is taken in.
///
/// Derived, not chosen: the `L^p` interpolation bound of Loseille & Alauzet
/// (NASA NTRS 20200003084, the restatement `common::metric` was written from)
/// carries `∫ (det|H|)^{p/(2p+d)}`, and with `g = √|det H|` that integrand is
/// `g^{2p/(2p+d)}`. So `τ = 2p/(2p+d) = 4/7` at `p = 2`, `d = 3`. The gap is a
/// ratio of two `l^τ` norms over one population, so the counting measure's
/// uniform cell volume cancels and no Riemann weight is needed.
///
/// `am_gm_gap_tau2` and `am_gm_gap_tau_inf` record the same gap at `τ = 2` and
/// `τ = ∞`, with their own rank correlations, because C2's verdict must not
/// depend on this one number without saying so.
const TAU: f64 = 2.0 * P_NORM / (2.0 * P_NORM + 3.0);

/// A grid sample joins the surface band when `|f| <= BAND_CELLS · h`. P-146's
/// value (`experiment_p146.rs:290`).
const BAND_CELLS: f64 = 1.0;

/// Fewest samples any axis of the anisotropic grid may carry. P-146's value
/// (`experiment_p146.rs:300`).
const MIN_SAMPLES: u32 = 5;

/// C1's and C3's bar: `|Δexponent| < 0.1` is "statistically indistinguishable",
/// `> 0.1` is "an exponent difference".
const EXPONENT_BAR: f64 = 0.1;

/// C2's bar on the Spearman correlation.
const CORRELATION_BAR: f64 = 0.7;

/// The registered vacuity control's bar: the gap must vary by at least `2×`.
const GAP_SPREAD_FLOOR: f64 = 2.0;

/// The anisotropic arm must be anisotropic somewhere, by at least this factor.
/// P-146's value (`experiment_p146.rs:315`).
const AXIS_RATIO_FLOOR: f64 = 1.5;

/// `cos(5°)`. P-146's axis-alignment test, carried so
/// `flat_direction_axis_aligned_fraction` is the same column P-146 reported.
const AXIS_ALIGNED_COS: f64 = 0.996_194_698_091_745_5;

/// The ladder must move the instrument by at least this factor on any field a
/// clause reads. P-146's value (`experiment_p146.rs:322`).
const LADDER_SPAN_FLOOR: f64 = 1.2;

/// `for_each_reference_field!` yields eight (`fields/mod.rs:211-255`).
const FIELDS: usize = 8;

/// C3's field, named by the clause itself.
const C3_FIELD: &str = "csg_difference";

/// The `instrument` column's value.
const INSTRUMENT: &str = "field_to_mesh_max";

/// C1's population, decided before the run: is `f` continuously differentiable
/// **on its own zero set**?
///
/// The reasoning for every row is the table in the header. Short version: a
/// gradient jump *on the surface* drops linear edge interpolation from second
/// order to first, and five of the eight reference fields have one — four from a
/// `min`/`max` combinator and one from being a polyhedron.
const SMOOTH: [(&str, bool); FIELDS] = [
    ("sphere", true),
    ("torus", true),
    ("box_exact", false),
    ("csg_difference", false),
    ("thin_plate", false),
    ("gyroid", false),
    ("fbm_terrain", true),
    ("noise_cavity", false),
];

/// Whether a field is in C1's population.
///
/// Panics for a name that was not classified: an unclassified field is a C1
/// clause with an undefined population, and defaulting it would decide C1 by
/// choosing which fields it is quantified over.
fn is_smooth(name: &str) -> bool {
    for (field, smooth) in SMOOTH {
        if field == name {
            return smooth;
        }
    }
    panic!("P-147: field `{name}` is not classified in SMOOTH, so C1's population is undefined");
}

/// The name of a `FieldBound` variant, without its parameters — the CSV writer
/// refuses a `,` inside a value and `Lipschitz { l: 3.46 }` has one.
fn bound_name(bound: FieldBound) -> &'static str {
    match bound {
        FieldBound::Exact => "Exact",
        FieldBound::Lipschitz { .. } => "Lipschitz",
        FieldBound::Underestimate { .. } => "Underestimate",
        FieldBound::Unbounded => "Unbounded",
    }
}

// ─── series formatting ───────────────────────────────────────────────────────

/// A ladder of reals as one CSV-safe token, `a|b|c|d|e`, in scientific form.
fn series_e(values: &[f64]) -> String {
    let mut out = String::new();
    for (index, value) in values.iter().enumerate() {
        if index > 0 {
            out.push('|');
        }
        out.push_str(&format!("{value:.6e}"));
    }
    out
}

/// A ladder of reals as one CSV-safe token, in fixed form.
fn series_f(values: &[f64]) -> String {
    let mut out = String::new();
    for (index, value) in values.iter().enumerate() {
        if index > 0 {
            out.push('|');
        }
        out.push_str(&format!("{value:.6}"));
    }
    out
}

/// A ladder of counts as one CSV-safe token.
fn series_u64(values: &[u64]) -> String {
    let mut out = String::new();
    for (index, value) in values.iter().enumerate() {
        if index > 0 {
            out.push('|');
        }
        out.push_str(&value.to_string());
    }
    out
}

/// A ladder of already-formatted tokens as one CSV-safe token.
fn series_str(values: &[String]) -> String {
    values.join("|")
}

/// Median by [`f64::total_cmp`]. Panics on an empty sample: a median of nothing
/// is not zero, it is unasked.
fn median(values: &mut [f64]) -> f64 {
    assert!(!values.is_empty(), "P-147: median over an empty sample");
    values.sort_by(f64::total_cmp);
    let n = values.len();
    if n % 2 == 1 {
        values[n / 2]
    } else {
        f64::midpoint(values[n / 2 - 1], values[n / 2])
    }
}

/// `max / min` over a sample, or `INFINITY` when the minimum is zero.
fn spread(values: &[f64]) -> f64 {
    let mut lo = f64::INFINITY;
    let mut hi = 0.0f64;
    for &value in values {
        lo = lo.min(value);
        hi = hi.max(value);
    }
    if lo > 0.0 { hi / lo } else { f64::INFINITY }
}

// ─── the metric field, restated from experiment_p146.rs ──────────────────────

/// The grid samples within one cell of the surface, swept `z`, `y`, `x` with
/// `x` innermost.
///
/// Restated verbatim from `experiment_p146.rs:413`. The sweep order is
/// load-bearing: the per-axis weights are means taken in this order, so a
/// different order gives a different last bit and a different rounded grid.
fn band_points<F>(field: &F, origin: [f64; 3], h: f64, samples: u32) -> Vec<[f64; 3]>
where
    F: Sdf<Scalar = f64>,
{
    let band = BAND_CELLS * h;
    let mut out = Vec::new();
    for k in 0..samples {
        for j in 0..samples {
            for i in 0..samples {
                let p = [
                    origin[0] + f64::from(i) * h,
                    origin[1] + f64::from(j) * h,
                    origin[2] + f64::from(k) * h,
                ];
                if field.sample(p).abs() <= band {
                    out.push(p);
                }
            }
        }
    }
    out
}

/// What the metric field says about one `(field, resolution)`.
struct Census {
    /// Band points.
    points: usize,
    /// Mean `√(e_aᵀ M e_a)` per axis: the metric's own point density along each
    /// world axis, and the only thing the anisotropic split is derived from.
    weights: [f64; 3],
    /// Band points whose smallest `|Hessian eigenvalue|` sits at `H_FLOOR`.
    at_floor: usize,
    /// Band points with a floored eigenvalue **and** at least one un-floored
    /// one — a genuinely exploitable flat direction. This is
    /// `flat_direction_fraction`'s numerator.
    flat_direction: usize,
    /// The same, additionally requiring the floored eigenvector within 5° of a
    /// world axis. P-146's `flat_axis_aligned_fraction`, carried for
    /// commensurability.
    flat_axis_aligned: usize,
    /// `‖√|det H|‖_{l^τ} / ‖tr|H|/3‖_{l^τ}` at `τ = TAU`, `2` and `∞`.
    gap_tau: f64,
    gap_tau2: f64,
    gap_tau_inf: f64,
}

/// Census the Hessian population and the metric field built from it.
///
/// The Hessians are passed in rather than recomputed. P-146 recomputed them
/// (`experiment_p146.rs:499`) because its census and its timed stage had to be
/// separable; there is no timed stage here, and `hessian` is a pure function of
/// `(field, p, h)`, so one evaluation and two readers is the same arithmetic.
///
/// # Why `flat_direction` requires an un-floored eigenvalue
///
/// A point flat in *every* direction — the middle of a box face, where all three
/// second differences vanish — prescribes no direction at all, and counting it
/// would make every polyhedron look like a heightfield. P-146 made the same
/// choice at `experiment_p146.rs:531-534` and its reason is quoted there.
fn census_of(hessians: &[Sym3], metrics: &[Sym3]) -> Census {
    let mut census = Census {
        points: hessians.len(),
        weights: [0.0; 3],
        at_floor: 0,
        flat_direction: 0,
        flat_axis_aligned: 0,
        gap_tau: 0.0,
        gap_tau2: 0.0,
        gap_tau_inf: 0.0,
    };

    for (hess, m) in hessians.iter().zip(metrics) {
        let (values, vectors) = hess.eigen();

        // `eigen` sorts by value and not by magnitude, so a saddle's smallest
        // magnitude is not `values[0]`.
        let mut flat = 0usize;
        let mut stiff = 0usize;
        for (index, value) in values.iter().enumerate() {
            if value.abs() < values[flat].abs() {
                flat = index;
            }
            if value.abs() > H_FLOOR {
                stiff += 1;
            }
        }
        let floored = values[flat].abs() <= H_FLOOR;

        if floored {
            census.at_floor += 1;
        }
        if floored && stiff > 0 {
            census.flat_direction += 1;
            let vector = [vectors[0][flat], vectors[1][flat], vectors[2][flat]];
            let mut axis = 0usize;
            for (index, component) in vector.iter().enumerate() {
                if component.abs() > vector[axis].abs() {
                    axis = index;
                }
            }
            if vector[axis].abs() >= AXIS_ALIGNED_COS {
                census.flat_axis_aligned += 1;
            }
        }

        for (axis, weight) in census.weights.iter_mut().enumerate() {
            *weight += m.get(axis, axis).sqrt();
        }
    }

    let n = hessians.len() as f64;
    for weight in &mut census.weights {
        *weight /= n;
    }
    census.gap_tau = am_gm_gap(hessians, TAU);
    census.gap_tau2 = am_gm_gap(hessians, 2.0);
    census.gap_tau_inf = am_gm_gap(hessians, f64::INFINITY);
    census
}

// ─── the anisotropic grid, restated from experiment_p146.rs ──────────────────

/// Nearest odd integer, at least one. Ties go up, deterministically.
/// Restated from `experiment_p146.rs:587`.
fn round_odd(x: f64) -> u32 {
    let half = ((x - 1.0) * 0.5).round();
    (2.0f64.mul_add(half, 1.0)).max(1.0) as u32
}

/// Per-axis sample counts from the metric's per-axis point densities, at the
/// isotropic arm's total budget.
///
/// Restated verbatim from `experiment_p146.rs:620`, including the two facts its
/// author measured and wrote down there:
///
/// - `n_a ∝ weights[a]` with `∏ n_a = N³`, so the two arms differ in **shape**
///   only and never in budget.
/// - **There is exactly one clamp and it is a lower one.** An axis below
///   [`MIN_SAMPLES`] is pinned there and the remaining budget is *re-solved* over
///   the axes still free. An upper clamp was tried and removed: it can bind on
///   two axes in one round before the lower pins resolve, and P-146 measured
///   `29x5x29` against a budget of `9³` — **5.77× over**, rising to `39.5×` at
///   `17³`.
///
/// Returns the counts and how many axes were pinned at the floor.
fn anisotropic_grid(weights: [f64; 3], samples: u32) -> ([u32; 3], usize) {
    let budget = f64::from(samples).powi(3);
    let mut pinned = [false; 3];
    let mut n = [samples; 3];

    for _round in 0..3 {
        let free: Vec<usize> = (0..3).filter(|&axis| !pinned[axis]).collect();
        if free.is_empty() {
            break;
        }
        let mut held = 1.0f64;
        for (axis, &fixed) in pinned.iter().enumerate() {
            if fixed {
                held *= f64::from(n[axis]);
            }
        }
        let count = free.len() as f64;
        let target = budget / held;
        let mut logsum = 0.0f64;
        for &axis in &free {
            logsum += weights[axis].ln();
        }
        let geometric_mean = (logsum / count).exp();
        let scale = target.powf(1.0 / count) / geometric_mean;

        let mut newly_pinned = false;
        for &axis in &free {
            let raw = round_odd(scale * weights[axis]);
            if raw < MIN_SAMPLES {
                n[axis] = MIN_SAMPLES;
                pinned[axis] = true;
                newly_pinned = true;
            } else {
                n[axis] = raw;
            }
        }
        if !newly_pinned {
            break;
        }
    }

    (n, pinned.iter().filter(|fixed| **fixed).count())
}

/// The field seen through a per-axis coordinate stretch.
///
/// Restated from `experiment_p146.rs:676`. `sample(q) = f(lo + q ⊙ s)`, so
/// extracting this on a **cubic** grid of `cell_size = h` in `q` is exactly
/// extracting `f` on a rectilinear grid whose physical spacings are `h · s` —
/// the only way to reach a rectilinear anisotropic grid through an `extract`
/// that takes a scalar `cell_size` (`marching_cubes/mod.rs:193`).
///
/// `gradient` is deliberately left as `Sdf`'s default: it is read only for the
/// emitted normals, and no clause here reads a normal. Positions are mapped
/// back to world space by the caller before any measurement touches them, and
/// **the accuracy instrument is handed the real `field`, never this wrapper**,
/// so no stretched gradient reaches a distance.
struct Stretched<'a, F> {
    field: &'a F,
    lo: [f64; 3],
    s: [f64; 3],
}

impl<F> Sdf for Stretched<'_, F>
where
    F: Sdf<Scalar = f64>,
{
    type Scalar = f64;

    fn sample(&self, q: [f64; 3]) -> f64 {
        self.field.sample([
            q[0].mul_add(self.s[0], self.lo[0]),
            q[1].mul_add(self.s[1], self.lo[1]),
            q[2].mul_add(self.s[2], self.lo[2]),
        ])
    }
}

// ─── the power-law fit ───────────────────────────────────────────────────────

/// One arm's fitted error law `E = C · (dof / dof_ref)^(−α)`.
struct Fit {
    /// `α`, positive for a converging arm.
    exponent: f64,
    /// `C`: the fitted error, in world length units, at [`REFERENCE_DOF`].
    constant: f64,
    /// The un-normalised prefactor `P` of `E = P · dof^(−α)`, so the fit can be
    /// rebuilt without knowing [`REFERENCE_DOF`].
    prefactor: f64,
    /// Coefficient of determination in log-log space.
    r2: f64,
    /// 95% half-width on `α`: `T_95_DF3 · se(α)`.
    exponent_ci95: f64,
    /// 95% half-width on `ln C`. Multiply out with `exp` for a factor.
    constant_ln_ci95: f64,
}

/// Ordinary least squares of `ln E` on `x = ln(dof / dof_ref)`.
///
/// `ln E = a + b·x` with `α = −b` and `C = exp(a)`, so `C` is the fit's value at
/// `dof = dof_ref` and is a **length** whatever `α` comes out as. The header's
/// "units trap" section is why that matters: when the two arms' exponents
/// differ, their raw prefactors are quantities of different dimension and their
/// ratio is not a number.
///
/// Standard errors are the textbook OLS ones — `s² = SSres/(n−2)`,
/// `se(b) = s/√Sxx`, `se(a) = s·√(1/n + x̄²/Sxx)` — scaled by [`T_95_DF3`].
///
/// # Panics
///
/// If any error is not strictly positive; `ln 0` is not a data point and the
/// caller's vacuity control has already established otherwise.
fn fit_power_law(dof: &[u64], error: &[f64]) -> Fit {
    assert_eq!(
        dof.len(),
        error.len(),
        "P-147: fit_power_law needs paired ladders"
    );
    let n = dof.len() as f64;
    let reference = f64::from(REFERENCE_DOF);

    let mut xs: Vec<f64> = Vec::with_capacity(dof.len());
    let mut ys: Vec<f64> = Vec::with_capacity(error.len());
    for (&d, &e) in dof.iter().zip(error) {
        assert!(
            e > 0.0 && e.is_finite(),
            "P-147: fit_power_law got a non-positive error {e}"
        );
        xs.push((d as f64 / reference).ln());
        ys.push(e.ln());
    }

    let mean_x = xs.iter().sum::<f64>() / n;
    let mean_y = ys.iter().sum::<f64>() / n;
    let mut sxx = 0.0f64;
    let mut sxy = 0.0f64;
    for (&x, &y) in xs.iter().zip(&ys) {
        let dx = x - mean_x;
        sxx = dx.mul_add(dx, sxx);
        sxy = dx.mul_add(y - mean_y, sxy);
    }
    let slope = sxy / sxx;
    let intercept = slope.mul_add(-mean_x, mean_y);

    let mut ssres = 0.0f64;
    let mut sstot = 0.0f64;
    for (&x, &y) in xs.iter().zip(&ys) {
        let residual = y - slope.mul_add(x, intercept);
        ssres = residual.mul_add(residual, ssres);
        sstot = (y - mean_y).mul_add(y - mean_y, sstot);
    }
    // `sstot == 0` means every rung produced the same error. The ladder-span
    // control has excluded that on every field a clause reads, and `0` is the
    // honest reading of "the fit explains none of a variance that does not
    // exist" for the fields it has not.
    let r2 = if sstot > 0.0 {
        1.0 - ssres / sstot
    } else {
        0.0
    };
    let sigma = (ssres / (n - 2.0)).sqrt();
    let se_slope = sigma / sxx.sqrt();
    let se_intercept = sigma * (1.0 / n + mean_x * mean_x / sxx).sqrt();

    let exponent = -slope;
    Fit {
        exponent,
        constant: intercept.exp(),
        prefactor: exponent.mul_add(reference.ln(), intercept).exp(),
        r2,
        exponent_ci95: T_95_DF3 * se_slope,
        constant_ln_ci95: T_95_DF3 * se_intercept,
    }
}

// ─── one arm, one rung ───────────────────────────────────────────────────────

/// One arm at one rung.
struct Arm {
    /// `n_x · n_y · n_z`: the sample budget, the fit's abscissa.
    dof: u64,
    triangles: u64,
    /// `field_to_mesh.max` — the instrument. See the header.
    error_max: f64,
    /// `field_to_mesh.mean`, fitted separately so an `L^∞` verdict can be
    /// checked against an `L^1`-ish one.
    error_mean: f64,
    /// `mesh_to_field.max`, a distance only where `bound()` is `Exact`.
    mesh_to_field_max: f64,
    /// `symmetric_hausdorff()`, the column P-146's numbers are reproduced
    /// against.
    symmetric: f64,
    /// Seeds that produced a distance: the instrument's coverage.
    seed_samples: u64,
}

/// Everything one `(field, resolution)` produced.
struct Rung {
    samples: u32,
    grid: [u32; 3],
    pinned: usize,
    axis_ratio: f64,
    budget_ratio: f64,
    census: Census,
    iso: Arm,
    aniso: Arm,
}

impl Rung {
    /// The two arms are the same grid, so the two meshes are bit-identical and
    /// every difference between them is exactly zero by construction.
    fn arms_identical(&self) -> bool {
        self.grid == [self.samples; 3]
    }
}

/// One field's whole ladder and the fits it decides.
struct FieldRow {
    name: &'static str,
    bound: &'static str,
    smooth: bool,
    rungs: Vec<Rung>,
    fit_iso: Fit,
    fit_aniso: Fit,
    fit_mean_iso: Fit,
    fit_mean_aniso: Fit,
    fit_tri_iso: Fit,
    fit_tri_aniso: Fit,
    wall_seconds: f64,
}

impl FieldRow {
    /// `|α_aniso − α_iso|`, the quantity C1 and C3 both read.
    fn exponent_difference(&self) -> f64 {
        (self.fit_aniso.exponent - self.fit_iso.exponent).abs()
    }

    /// `C_aniso / C_iso`: the anisotropic arm's fitted error at the reference
    /// budget, as a fraction of the isotropic arm's. Below one is an
    /// improvement.
    fn constant_ratio(&self) -> f64 {
        self.fit_aniso.constant / self.fit_iso.constant
    }

    /// The registered `am_gm_gap`: the finest rung's, which is the
    /// best-resolved Hessian population on the ladder.
    fn gap(&self) -> f64 {
        self.last().census.gap_tau
    }

    /// The finest rung, whose census supplies every per-field scalar.
    fn last(&self) -> &Rung {
        self.rungs
            .last()
            .expect("P-147: every field has at least one rung")
    }

    /// Every rung's two arms are the same grid.
    fn arms_identical(&self) -> bool {
        self.rungs.iter().all(Rung::arms_identical)
    }

    /// The most anisotropic grid this field's metric asked for.
    fn axis_ratio_max(&self) -> f64 {
        self.rungs
            .iter()
            .map(|rung| rung.axis_ratio)
            .fold(0.0f64, f64::max)
    }
}

/// Measure one reference field across the whole ladder.
fn measure_field<F>(field: &F, name: &'static str) -> FieldRow
where
    F: ReferenceField + Sdf<Scalar = f64>,
{
    let started = Instant::now();
    let bound = field.bound();
    let (lo, hi) = field.domain();
    let extent = hi[0] - lo[0];

    let mut mc = MarchingCubes::<f64>::new();
    let mut iso_mesh = MeshBuffer::<f64>::new();
    let mut aniso_mesh = MeshBuffer::<f64>::new();
    let mut rungs: Vec<Rung> = Vec::with_capacity(RUNGS.len());

    for samples in RUNGS {
        let (shape, origin, h) = common::grid::<f64, _>(field, samples);

        // ── the metric field ────────────────────────────────────────────────
        let points = band_points(field, origin, h, samples);
        assert!(
            !points.is_empty(),
            "VOID: {name} at {samples}^3 put no grid sample within {BAND_CELLS} cell of its \
             surface, so am_gm_gap and flat_direction_fraction would be statistics of an empty \
             population and every zero among them a zero that could not have been non-zero (M-44)"
        );
        let hessians: Vec<Sym3> = points.iter().map(|&p| hessian(field, p, h)).collect();
        let metrics: Vec<Sym3> = hessians.iter().map(|hess| metric_lp(hess, P_NORM)).collect();
        let census = census_of(&hessians, &metrics);

        // ── the isotropic arm ───────────────────────────────────────────────
        iso_mesh.reset();
        mc.extract(field, &shape, origin, h, &mut iso_mesh)
            .expect("isotropic extraction over the reference grid");

        // ── the anisotropic arm, at the same total sample budget ────────────
        let (grid, pinned) = anisotropic_grid(census.weights, samples);
        let stretch = [
            extent / f64::from(grid[0] - 1) / h,
            extent / f64::from(grid[1] - 1) / h,
            extent / f64::from(grid[2] - 1) / h,
        ];
        let stretched = Stretched {
            field,
            lo,
            s: stretch,
        };
        let aniso_shape = RuntimeShape3::new(grid).expect("anisotropic grid fits u32");
        aniso_mesh.reset();
        mc.extract(&stretched, &aniso_shape, [0.0; 3], h, &mut aniso_mesh)
            .expect("anisotropic extraction over the stretched grid");
        for p in &mut aniso_mesh.positions {
            p[0] = p[0].mul_add(stretch[0], lo[0]);
            p[1] = p[1].mul_add(stretch[1], lo[1]);
            p[2] = p[2].mul_add(stretch[2], lo[2]);
        }

        // ── one instrument, one seed lattice, both arms ─────────────────────
        let cfg = AccuracyConfig::from_cell_size(h).expect("positive cell size");
        let iso_report = accuracy(
            &iso_mesh.positions,
            &iso_mesh.indices,
            field,
            &shape,
            origin,
            &cfg,
        )
        .expect("accuracy over the isotropic arm");
        let aniso_report = accuracy(
            &aniso_mesh.positions,
            &aniso_mesh.indices,
            field,
            &shape,
            origin,
            &cfg,
        )
        .expect("accuracy over the anisotropic arm on the isotropic seed lattice");

        let arm = |mesh: &MeshBuffer<f64>, report: &isomesh::validate::AccuracyReport<f64>,
                   dof: u64| Arm {
            dof,
            triangles: mesh.triangle_count() as u64,
            error_max: report.field_to_mesh.max,
            error_mean: report.field_to_mesh.mean,
            mesh_to_field_max: report.mesh_to_field.max,
            symmetric: report.symmetric_hausdorff(),
            seed_samples: report.field_to_mesh.samples,
        };
        let iso_dof = u64::from(samples) * u64::from(samples) * u64::from(samples);
        let aniso_dof = u64::from(grid[0]) * u64::from(grid[1]) * u64::from(grid[2]);
        let iso = arm(&iso_mesh, &iso_report, iso_dof);
        let aniso = arm(&aniso_mesh, &aniso_report, aniso_dof);

        for (label, side) in [("isotropic", &iso), ("anisotropic", &aniso)] {
            assert!(
                side.triangles > 0,
                "VOID: {name}'s {label} arm emitted no triangle at {samples}^3, so its error is \
                 the distance from the surface to nothing and every exponent fitted through it is \
                 a fit through an absence (M-44)"
            );
            assert!(
                side.seed_samples > 0,
                "VOID: {name}'s {label} arm at {samples}^3 had no seed reach the surface, so \
                 field_to_mesh.max is zero because nothing was measured rather than because the \
                 mesh is exact (M-44)"
            );
            assert!(
                side.error_max > 0.0 && side.error_max.is_finite(),
                "VOID: {name}'s {label} arm at {samples}^3 measured error_max = {} over {} \
                 seeds. `ln 0` is not a data point and a zero error here is the instrument \
                 failing, not the mesh succeeding (M-44)",
                side.error_max,
                side.seed_samples
            );
        }

        let axis_hi = f64::from(grid.iter().copied().max().unwrap_or(samples));
        let axis_lo = f64::from(grid.iter().copied().min().unwrap_or(samples));
        rungs.push(Rung {
            samples,
            grid,
            pinned,
            axis_ratio: axis_hi / axis_lo,
            budget_ratio: aniso_dof as f64 / iso_dof as f64,
            census,
            iso,
            aniso,
        });
    }

    let dof_iso: Vec<u64> = rungs.iter().map(|rung| rung.iso.dof).collect();
    let dof_aniso: Vec<u64> = rungs.iter().map(|rung| rung.aniso.dof).collect();
    let tri_iso: Vec<u64> = rungs.iter().map(|rung| rung.iso.triangles).collect();
    let tri_aniso: Vec<u64> = rungs.iter().map(|rung| rung.aniso.triangles).collect();
    let max_iso: Vec<f64> = rungs.iter().map(|rung| rung.iso.error_max).collect();
    let max_aniso: Vec<f64> = rungs.iter().map(|rung| rung.aniso.error_max).collect();
    let mean_iso: Vec<f64> = rungs.iter().map(|rung| rung.iso.error_mean).collect();
    let mean_aniso: Vec<f64> = rungs.iter().map(|rung| rung.aniso.error_mean).collect();

    FieldRow {
        name,
        bound: bound_name(bound),
        smooth: is_smooth(name),
        fit_iso: fit_power_law(&dof_iso, &max_iso),
        fit_aniso: fit_power_law(&dof_aniso, &max_aniso),
        fit_mean_iso: fit_power_law(&dof_iso, &mean_iso),
        fit_mean_aniso: fit_power_law(&dof_aniso, &mean_aniso),
        fit_tri_iso: fit_power_law(&tri_iso, &max_iso),
        fit_tri_aniso: fit_power_law(&tri_aniso, &max_aniso),
        rungs,
        wall_seconds: started.elapsed().as_secs_f64(),
    }
}

// ─── the reproduction check against p-146.csv ────────────────────────────────

/// One row of the committed `docs/experiments/p-146.csv`, as text.
///
/// Kept as strings so the comparison is byte-for-byte against what P-146 wrote
/// rather than against a re-parse and re-format of it.
struct Baseline {
    field: String,
    resolution: u32,
    triangles_isotropic: String,
    triangles_anisotropic: String,
    grid_anisotropic: String,
    hausdorff_isotropic: String,
    hausdorff_anisotropic: String,
}

/// Parse the committed `p-146.csv`.
///
/// Hand-rolled, in the shape `src/golden.rs:245`'s `field_of` uses: skip the
/// `#` provenance block, read the header, resolve the five column names, split
/// the rest on `,`. `common::experiment` refuses a `,` inside any value
/// (`benches/common/experiment.rs:69-75`), so a positional split is exact and
/// needs no quoting rules.
///
/// # Panics
///
/// If the file is absent, headerless, missing a column this check reads, or has
/// a row whose cell count disagrees with the header. Every one of those means
/// the reproduction cannot be checked, and an unchecked claim of
/// commensurability is the thing this function exists to prevent.
fn read_p146_baseline() -> Vec<Baseline> {
    let path = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../docs/experiments/p-146.csv");
    let text = std::fs::read_to_string(&path).unwrap_or_else(|why| {
        panic!(
            "P-147 reproduces P-146's construction and must read its committed CSV at {}: {why}",
            path.display()
        )
    });

    let mut lines = text.lines().filter(|line| !line.starts_with('#') && !line.is_empty());
    let header: Vec<&str> = lines
        .next()
        .expect("p-146.csv has a header line after its provenance block")
        .split(',')
        .collect();
    let column = |name: &str| -> usize {
        header
            .iter()
            .position(|candidate| *candidate == name)
            .unwrap_or_else(|| {
                panic!("p-146.csv has no `{name}` column, which P-147's reproduction check reads")
            })
    };
    let c_field = column("field");
    let c_resolution = column("resolution");
    let c_tri_iso = column("triangles_isotropic");
    let c_tri_aniso = column("triangles_anisotropic");
    let c_grid = column("grid_anisotropic");
    let c_haus_iso = column("hausdorff_isotropic");
    let c_haus_aniso = column("hausdorff_anisotropic");

    let mut out = Vec::new();
    for line in lines {
        let cells: Vec<&str> = line.split(',').collect();
        assert_eq!(
            cells.len(),
            header.len(),
            "p-146.csv row has {} cells against a header of {}: {line}",
            cells.len(),
            header.len()
        );
        out.push(Baseline {
            field: cells[c_field].to_string(),
            resolution: cells[c_resolution]
                .parse()
                .expect("p-146.csv `resolution` is an integer"),
            triangles_isotropic: cells[c_tri_iso].to_string(),
            triangles_anisotropic: cells[c_tri_aniso].to_string(),
            grid_anisotropic: cells[c_grid].to_string(),
            hausdorff_isotropic: cells[c_haus_iso].to_string(),
            hausdorff_anisotropic: cells[c_haus_aniso].to_string(),
        });
    }
    out
}

/// How many comparisons the reproduction check made and how many disagreed.
struct Reproduction {
    checks: usize,
    mismatches: Vec<String>,
}

/// Compare this run's arms against P-146's committed numbers.
///
/// Triangle counts and grid shapes are exact integers and are compared on every
/// row. `hausdorff_*` is compared only where P-146 recorded a number rather than
/// `unmeasurable:bound=…`, and at the `{:.9}` precision P-146 wrote it at
/// (`experiment_p146.rs:1240`).
fn check_reproduction(rows: &[FieldRow], baseline: &[Baseline]) -> Reproduction {
    let mut out = Reproduction {
        checks: 0,
        mismatches: Vec::new(),
    };
    for row in rows {
        for rung in &row.rungs {
            let found = baseline
                .iter()
                .find(|entry| entry.field == row.name && entry.resolution == rung.samples);
            let Some(entry) = found else {
                out.mismatches.push(format!(
                    "{} at {}^3 has no row in p-146.csv",
                    row.name, rung.samples
                ));
                continue;
            };

            let mut compare = |what: &str, theirs: &str, ours: String| {
                out.checks += 1;
                if theirs != ours {
                    out.mismatches.push(format!(
                        "{} at {}^3 {what}: p-146.csv says `{theirs}`, this run says `{ours}`",
                        row.name, rung.samples
                    ));
                }
            };
            compare(
                "triangles_isotropic",
                &entry.triangles_isotropic,
                rung.iso.triangles.to_string(),
            );
            compare(
                "triangles_anisotropic",
                &entry.triangles_anisotropic,
                rung.aniso.triangles.to_string(),
            );
            compare(
                "grid_anisotropic",
                &entry.grid_anisotropic,
                format!("{}x{}x{}", rung.grid[0], rung.grid[1], rung.grid[2]),
            );
            if !entry.hausdorff_isotropic.starts_with("unmeasurable") {
                compare(
                    "hausdorff_isotropic",
                    &entry.hausdorff_isotropic,
                    format!("{:.9}", rung.iso.symmetric),
                );
            }
            if !entry.hausdorff_anisotropic.starts_with("unmeasurable") {
                compare(
                    "hausdorff_anisotropic",
                    &entry.hausdorff_anisotropic,
                    format!("{:.9}", rung.aniso.symmetric),
                );
            }
        }
    }
    out
}

// ─── console report ──────────────────────────────────────────────────────────

/// One field's line on the console.
fn report(row: &FieldRow) {
    println!(
        "{:<15} {:<14} smooth={:<5} a_iso={:>8.4}+-{:<7.4} a_ani={:>8.4}+-{:<7.4} \
         |d|={:>8.4} C_iso={:>11.4e} C_ani={:>11.4e} ratio={:>10.4} gap={:>11.4e} \
         flat={:>8.4} axis_ratio<={:>9.3} r2={:>7.4}/{:<7.4} identical={}",
        row.name,
        row.bound,
        row.smooth,
        row.fit_iso.exponent,
        row.fit_iso.exponent_ci95,
        row.fit_aniso.exponent,
        row.fit_aniso.exponent_ci95,
        row.exponent_difference(),
        row.fit_iso.constant,
        row.fit_aniso.constant,
        row.constant_ratio(),
        row.gap(),
        row.last().census.flat_direction as f64 / row.last().census.points as f64,
        row.axis_ratio_max(),
        row.fit_iso.r2,
        row.fit_aniso.r2,
        row.arms_identical()
    );
}

// ─── the run ─────────────────────────────────────────────────────────────────

fn main() {
    if !std::env::args().any(|arg| arg == "--bench") {
        return;
    }
    let () = LADDER_RUNGS_ARE_FIVE;
    let prereg = isomesh::experiment!("P-147");

    common::experiment::run(prereg, |run| {
        println!(
            "construction: P-146's metric-driven anisotropic sampling GRID, per-axis global, NOT \
             per-cell -- restated verbatim and checked against docs/experiments/p-146.csv.\n  \
             metric M_Lp_hessian, p = {P_NORM}, band |f| <= {BAND_CELLS}h, H_FLOOR = \
             {H_FLOOR:e}\n  ladder {RUNGS:?} samples/axis; abscissa dof = nx*ny*nz; constants read \
             off at dof_ref = {REFERENCE_DOF}\n  instrument {INSTRUMENT} = \
             AccuracyReport::field_to_mesh.max, ONE seed lattice per rung shared by both arms\n  \
             am_gm_gap at tau = 2p/(2p+d) = {TAU:.6}, over the finest rung's band Hessians\n"
        );

        let mut rows: Vec<FieldRow> = Vec::with_capacity(FIELDS);
        isomesh::for_each_reference_field!(f64, |name, field| {
            rows.push(measure_field(&field, name));
        });
        assert_eq!(
            rows.len(),
            FIELDS,
            "P-147: for_each_reference_field! must yield {FIELDS} fields"
        );

        for row in &rows {
            report(row);
        }
        println!();

        // ── the reproduction check ───────────────────────────────────────────
        //
        // Not a vacuity control: a mismatch here means this file is not running
        // P-146's arm and the two rows are not commensurable, which is a defect
        // rather than an unmeasured clause.
        let baseline = read_p146_baseline();
        let reproduction = check_reproduction(&rows, &baseline);
        assert!(
            reproduction.mismatches.is_empty(),
            "P-147 restates P-146's construction and must reproduce its committed numbers \
             exactly. {} of {} comparisons disagreed:\n  {}",
            reproduction.mismatches.len(),
            reproduction.checks,
            reproduction.mismatches.join("\n  ")
        );
        println!(
            "reproduction: {} comparisons against docs/experiments/p-146.csv, 0 mismatches\n",
            reproduction.checks
        );

        // ── vacuity controls, all before the first record ────────────────────
        //
        // M-44: a zero that could not have been non-zero is not a measurement.
        let gaps: Vec<f64> = rows.iter().map(FieldRow::gap).collect();
        let gap_lo = gaps.iter().copied().fold(f64::INFINITY, f64::min);
        let gap_hi = gaps.iter().copied().fold(0.0f64, f64::max);
        let gap_spread = spread(&gaps);
        assert!(
            gap_hi > 0.0,
            "VOID: the AM-GM gap is zero on every one of the {FIELDS} fields, so its numerator \
             ||sqrt(|det H|)|| vanished everywhere and C2's rank correlation is between a \
             constant zero and the constant ratios"
        );
        assert!(
            gap_hi > GAP_SPREAD_FLOOR * gap_lo,
            "VOID: the AM-GM gap runs {gap_lo:e} to {gap_hi:e} across the {FIELDS} fields, a \
             spread of {gap_spread:.4}x, which does not exceed {GAP_SPREAD_FLOOR}x. That is the \
             registration's own vacuity control: C2's correlation is then against a constant"
        );

        let axis_ratio_max = rows
            .iter()
            .map(FieldRow::axis_ratio_max)
            .fold(0.0f64, f64::max);
        assert!(
            axis_ratio_max > AXIS_RATIO_FLOOR,
            "VOID: the most anisotropic grid this metric asked for anywhere is \
             {axis_ratio_max:.4}:1, below {AXIS_RATIO_FLOOR}:1. Both arms are then the same grid \
             under two names, every exponent_difference is zero by construction, and C1, C2 and \
             C3 are all measuring the extractor's run-to-run determinism"
        );

        let smooth: Vec<&FieldRow> = rows.iter().filter(|row| row.smooth).collect();
        let c1_population = smooth.len();
        assert!(
            smooth.iter().any(|row| !row.arms_identical()),
            "VOID: every one of C1's {c1_population} smooth fields got a byte-identical pair of \
             arms, so |Deltaexponent| is exactly 0 on all of them and c1_holds is true for the \
             same reason 0 < {EXPONENT_BAR} is true. C1 would then be a restatement of `the grid \
             did not change` (M-44)"
        );

        for row in rows.iter().filter(|row| row.smooth || row.name == C3_FIELD) {
            for (label, errors) in [
                (
                    "isotropic",
                    row.rungs
                        .iter()
                        .map(|rung| rung.iso.error_max)
                        .collect::<Vec<f64>>(),
                ),
                (
                    "anisotropic",
                    row.rungs
                        .iter()
                        .map(|rung| rung.aniso.error_max)
                        .collect::<Vec<f64>>(),
                ),
            ] {
                let span = spread(&errors);
                assert!(
                    span >= LADDER_SPAN_FLOOR,
                    "VOID: {}'s {label} arm moved its error by only {span:.4}x across the whole \
                     ladder {RUNGS:?}, under {LADDER_SPAN_FLOOR}x. A clause reads this field's \
                     exponent and an exponent read off a horizontal line is a slope through \
                     noise: {}",
                    row.name,
                    series_e(&errors)
                );
            }
        }

        // ── the global verdicts ──────────────────────────────────────────────
        //
        // All three clauses are quantified over the roster rather than over one
        // row, so all three verdict columns carry the same value on every row
        // and the per-row facts sit in neighbouring columns.
        let c1_holds = smooth
            .iter()
            .all(|row| row.exponent_difference() < EXPONENT_BAR);

        let ratios: Vec<f64> = rows.iter().map(FieldRow::constant_ratio).collect();
        let gaps_tau2: Vec<f64> = rows
            .iter()
            .map(|row| row.last().census.gap_tau2)
            .collect();
        let gaps_tau_inf: Vec<f64> = rows
            .iter()
            .map(|row| row.last().census.gap_tau_inf)
            .collect();
        let correlation = rank_correlation(&gaps, &ratios);
        let correlation_tau2 = rank_correlation(&gaps_tau2, &ratios);
        let correlation_tau_inf = rank_correlation(&gaps_tau_inf, &ratios);
        let improved = ratios.iter().filter(|ratio| **ratio < 1.0).count();
        // "The fitted constant improves" over the population, read as a strict
        // majority of the eight fields. A correlation computed over a population
        // in which nothing improved would be a correlation between two ways of
        // writing "no change", which is why the two conjuncts are both required.
        let c2_holds = improved * 2 > FIELDS && correlation >= CORRELATION_BAR;

        let c3_row = rows
            .iter()
            .find(|row| row.name == C3_FIELD)
            .expect("P-147: C3 names csg_difference and the roster must contain it");
        let c3_difference = c3_row.exponent_difference();
        let c3_holds = c3_difference > EXPONENT_BAR;
        let c3_arms_identical = c3_row.arms_identical();
        let c3_reason = if c3_arms_identical {
            format!(
                "arms_identical:axis_ratio_max={:.6}:per_axis_grid_cannot_grade_a_\
                 permutation_symmetric_field",
                c3_row.axis_ratio_max()
            )
        } else {
            format!("graded:axis_ratio_max={:.6}", c3_row.axis_ratio_max())
        };

        // The companion to C3, and explicitly not its verdict: does the field
        // that lacks W^{2,p} regularity converge more slowly than the smooth
        // ones do, under the *same* uniform refinement? That is what the
        // regularity deficit predicts and it needs no anisotropy to ask.
        let mut smooth_exponents: Vec<f64> =
            smooth.iter().map(|row| row.fit_iso.exponent).collect();
        let smooth_median = median(&mut smooth_exponents);
        let companion_deficit = smooth_median - c3_row.fit_iso.exponent;
        let companion_holds = companion_deficit > EXPONENT_BAR;

        println!(
            "C1: population {c1_population} of {FIELDS} smooth fields; \
             max |Deltaexponent| {:.6}; bar {EXPONENT_BAR} -> {c1_holds}",
            smooth
                .iter()
                .map(|row| row.exponent_difference())
                .fold(0.0f64, f64::max)
        );
        println!(
            "C2: {improved} of {FIELDS} fields improved their constant (bar: a strict majority, \
             {}); Spearman(am_gm_gap, constant_ratio) = {correlation:.6} at tau={TAU:.6}, \
             {correlation_tau2:.6} at tau=2, {correlation_tau_inf:.6} at tau=inf; bar \
             {CORRELATION_BAR} -> {c2_holds}",
            FIELDS / 2 + 1
        );
        println!(
            "C3: {C3_FIELD} |Deltaexponent| = {c3_difference:.6}; bar {EXPONENT_BAR} -> \
             {c3_holds} ({c3_reason})"
        );
        println!(
            "    companion (NOT C3's verdict): smooth median isotropic exponent {smooth_median:.6} \
             - {C3_FIELD}'s {:.6} = {companion_deficit:.6} -> {companion_holds}\n",
            c3_row.fit_iso.exponent
        );

        // ── the rows, one per field ──────────────────────────────────────────
        for row in &rows {
            let last = row.last();
            let points = last.census.points as f64;

            let samples: Vec<f64> = row
                .rungs
                .iter()
                .map(|rung| f64::from(rung.samples))
                .collect();
            let grids: Vec<String> = row
                .rungs
                .iter()
                .map(|rung| format!("{}x{}x{}", rung.grid[0], rung.grid[1], rung.grid[2]))
                .collect();
            let grids_iso: Vec<String> = row
                .rungs
                .iter()
                .map(|rung| format!("{0}x{0}x{0}", rung.samples))
                .collect();
            let max_iso: Vec<f64> = row.rungs.iter().map(|rung| rung.iso.error_max).collect();
            let max_aniso: Vec<f64> = row.rungs.iter().map(|rung| rung.aniso.error_max).collect();

            run.record(&[
                ("field", row.name.to_string()),
                (
                    "resolution_series",
                    series_f(&samples).replace(".000000", ""),
                ),
                (
                    "fitted_exponent_isotropic",
                    format!("{:.6}", row.fit_iso.exponent),
                ),
                (
                    "fitted_exponent_anisotropic",
                    format!("{:.6}", row.fit_aniso.exponent),
                ),
                ("exponent_difference", format!("{:.6}", row.exponent_difference())),
                (
                    "fitted_constant_isotropic",
                    format!("{:.6e}", row.fit_iso.constant),
                ),
                (
                    "fitted_constant_anisotropic",
                    format!("{:.6e}", row.fit_aniso.constant),
                ),
                ("constant_ratio", format!("{:.6}", row.constant_ratio())),
                ("am_gm_gap", format!("{:.6e}", row.gap())),
                (
                    "flat_direction_fraction",
                    format!("{:.6}", last.census.flat_direction as f64 / points),
                ),
                ("c1_holds", c1_holds.to_string()),
                ("c2_holds", c2_holds.to_string()),
                ("c3_holds", c3_holds.to_string()),
                // ── extras (M-273) ──
                (
                    "am_gm_gap_series",
                    series_e(
                        &row.rungs
                            .iter()
                            .map(|rung| rung.census.gap_tau)
                            .collect::<Vec<f64>>(),
                    ),
                ),
                ("am_gm_gap_spread", format!("{gap_spread:.6e}")),
                ("am_gm_gap_tau", format!("{TAU:.6}")),
                ("am_gm_gap_tau2", format!("{:.6e}", last.census.gap_tau2)),
                (
                    "am_gm_gap_tau_inf",
                    format!("{:.6e}", last.census.gap_tau_inf),
                ),
                ("arms_identical", row.arms_identical().to_string()),
                (
                    "at_floor_fraction",
                    format!("{:.6}", last.census.at_floor as f64 / points),
                ),
                ("axes_pinned_series", {
                    let pinned: Vec<u64> = row.rungs.iter().map(|r| r.pinned as u64).collect();
                    series_u64(&pinned)
                }),
                ("axis_ratio_max", format!("{:.6}", row.axis_ratio_max())),
                (
                    "axis_ratio_series",
                    series_f(
                        &row.rungs
                            .iter()
                            .map(|rung| rung.axis_ratio)
                            .collect::<Vec<f64>>(),
                    ),
                ),
                ("band_points_series", {
                    let counts: Vec<u64> =
                        row.rungs.iter().map(|r| r.census.points as u64).collect();
                    series_u64(&counts)
                }),
                (
                    "budget_ratio_series",
                    series_f(
                        &row.rungs
                            .iter()
                            .map(|rung| rung.budget_ratio)
                            .collect::<Vec<f64>>(),
                    ),
                ),
                ("c1_population", c1_population.to_string()),
                (
                    "c1_row_holds",
                    (row.exponent_difference() < EXPONENT_BAR).to_string(),
                ),
                ("c1_smooth", row.smooth.to_string()),
                ("c2_improved_fields", improved.to_string()),
                ("c2_rank_correlation", format!("{correlation:.6}")),
                (
                    "c2_rank_correlation_tau2",
                    format!("{correlation_tau2:.6}"),
                ),
                (
                    "c2_rank_correlation_tau_inf",
                    format!("{correlation_tau_inf:.6}"),
                ),
                (
                    "c2_row_improves",
                    (row.constant_ratio() < 1.0).to_string(),
                ),
                ("c3_arms_identical", c3_arms_identical.to_string()),
                ("c3_exponent_difference", format!("{c3_difference:.6}")),
                ("c3_field", C3_FIELD.to_string()),
                ("c3_reason", c3_reason.clone()),
                (
                    "companion_deficit_holds",
                    companion_holds.to_string(),
                ),
                (
                    "companion_regularity_deficit",
                    format!("{companion_deficit:.6}"),
                ),
                (
                    "constant_ln_ci95_anisotropic",
                    format!("{:.6}", row.fit_aniso.constant_ln_ci95),
                ),
                (
                    "constant_ln_ci95_isotropic",
                    format!("{:.6}", row.fit_iso.constant_ln_ci95),
                ),
                ("dof_series_anisotropic", {
                    let dof: Vec<u64> = row.rungs.iter().map(|r| r.aniso.dof).collect();
                    series_u64(&dof)
                }),
                ("dof_series_isotropic", {
                    let dof: Vec<u64> = row.rungs.iter().map(|r| r.iso.dof).collect();
                    series_u64(&dof)
                }),
                ("error_max_series_anisotropic", series_e(&max_aniso)),
                ("error_max_series_isotropic", series_e(&max_iso)),
                (
                    "error_mean_series_anisotropic",
                    series_e(
                        &row.rungs
                            .iter()
                            .map(|rung| rung.aniso.error_mean)
                            .collect::<Vec<f64>>(),
                    ),
                ),
                (
                    "error_mean_series_isotropic",
                    series_e(
                        &row.rungs
                            .iter()
                            .map(|rung| rung.iso.error_mean)
                            .collect::<Vec<f64>>(),
                    ),
                ),
                (
                    "exponent_ci95_anisotropic",
                    format!("{:.6}", row.fit_aniso.exponent_ci95),
                ),
                (
                    "exponent_ci95_isotropic",
                    format!("{:.6}", row.fit_iso.exponent_ci95),
                ),
                ("exponent_difference_ci95", {
                    let iso = row.fit_iso.exponent_ci95;
                    let aniso = row.fit_aniso.exponent_ci95;
                    format!("{:.6}", iso.hypot(aniso))
                }),
                ("exponent_difference_resolved", {
                    let iso = row.fit_iso.exponent_ci95;
                    let aniso = row.fit_aniso.exponent_ci95;
                    (row.exponent_difference() > iso.hypot(aniso)).to_string()
                }),
                (
                    "exponent_difference_vs_triangles",
                    format!(
                        "{:.6}",
                        (row.fit_tri_aniso.exponent - row.fit_tri_iso.exponent).abs()
                    ),
                ),
                (
                    "exponent_mean_anisotropic",
                    format!("{:.6}", row.fit_mean_aniso.exponent),
                ),
                (
                    "exponent_mean_isotropic",
                    format!("{:.6}", row.fit_mean_iso.exponent),
                ),
                (
                    "exponent_vs_triangles_anisotropic",
                    format!("{:.6}", row.fit_tri_aniso.exponent),
                ),
                (
                    "exponent_vs_triangles_isotropic",
                    format!("{:.6}", row.fit_tri_iso.exponent),
                ),
                ("field_bound", row.bound.to_string()),
                (
                    "fit_r2_anisotropic",
                    format!("{:.6}", row.fit_aniso.r2),
                ),
                ("fit_r2_isotropic", format!("{:.6}", row.fit_iso.r2)),
                (
                    "flat_direction_axis_aligned_fraction",
                    format!("{:.6}", last.census.flat_axis_aligned as f64 / points),
                ),
                ("grid_series_anisotropic", series_str(&grids)),
                ("grid_series_isotropic", series_str(&grids_iso)),
                (
                    "hausdorff_symmetric_series_anisotropic",
                    series_e(
                        &row.rungs
                            .iter()
                            .map(|rung| rung.aniso.symmetric)
                            .collect::<Vec<f64>>(),
                    ),
                ),
                (
                    "hausdorff_symmetric_series_isotropic",
                    series_e(
                        &row.rungs
                            .iter()
                            .map(|rung| rung.iso.symmetric)
                            .collect::<Vec<f64>>(),
                    ),
                ),
                ("instrument", INSTRUMENT.to_string()),
                (
                    "ladder_span_anisotropic",
                    format!("{:.6}", spread(&max_aniso)),
                ),
                ("ladder_span_isotropic", format!("{:.6}", spread(&max_iso))),
                (
                    "mesh_to_field_max_series_anisotropic",
                    series_e(
                        &row.rungs
                            .iter()
                            .map(|rung| rung.aniso.mesh_to_field_max)
                            .collect::<Vec<f64>>(),
                    ),
                ),
                (
                    "mesh_to_field_max_series_isotropic",
                    series_e(
                        &row.rungs
                            .iter()
                            .map(|rung| rung.iso.mesh_to_field_max)
                            .collect::<Vec<f64>>(),
                    ),
                ),
                (
                    "prefactor_anisotropic",
                    format!("{:.6e}", row.fit_aniso.prefactor),
                ),
                (
                    "prefactor_isotropic",
                    format!("{:.6e}", row.fit_iso.prefactor),
                ),
                ("reference_dof", REFERENCE_DOF.to_string()),
                ("reproduction_checks", reproduction.checks.to_string()),
                (
                    "reproduction_mismatches",
                    reproduction.mismatches.len().to_string(),
                ),
                ("seed_coverage_series_anisotropic", {
                    let seeds: Vec<u64> =
                        row.rungs.iter().map(|r| r.aniso.seed_samples).collect();
                    series_u64(&seeds)
                }),
                ("seed_coverage_series_isotropic", {
                    let seeds: Vec<u64> = row.rungs.iter().map(|r| r.iso.seed_samples).collect();
                    series_u64(&seeds)
                }),
                ("triangles_series_anisotropic", {
                    let tris: Vec<u64> = row.rungs.iter().map(|r| r.aniso.triangles).collect();
                    series_u64(&tris)
                }),
                ("triangles_series_isotropic", {
                    let tris: Vec<u64> = row.rungs.iter().map(|r| r.iso.triangles).collect();
                    series_u64(&tris)
                }),
                ("wall_seconds", format!("{:.3}", row.wall_seconds)),
            ]);
        }
    });
}
