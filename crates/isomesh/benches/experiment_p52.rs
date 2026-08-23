//! **P-52 — tangency-aware vertex placement, two iterations, clamped to the cell.**
//!
//! Ticket: R-047. Pre-registered at `db0ca10`, before this harness existed.
//!
//! ```bash
//! cargo bench --bench experiment_p52
//! ```
//!
//! Writes `docs/experiments/p-52.csv`: four fields × three rules = twelve rows.
//!
//! # This is not Sellán, Batty & Stein's algorithm, and the difference is large
//!
//! *Reach For the Spheres: Tangency-Aware Surface Reconstruction of SDFs*,
//! SIGGRAPH Asia 2023, `10.1145/3610548.3618196`, arXiv `2308.09813`. Their
//! method is a **global sparse solve with per-iteration remeshing, run to
//! convergence over a multi-resolution schedule** — their §4.1 halves `h`
//! repeatedly and stops only when the energy has not fallen by `1e-3·ε` in the
//! last 100 iterations, and their Eq. (11) damps every step with a mass-matrix
//! proximal term `‖V − V^{t−1}‖²_M / 2τ` that has no analogue here. Their own
//! Fig. 17 ablation measures **clamping away far spheres as progressive detail
//! loss**, and this harness clamps every vertex to its cell twice.
//!
//! So **none of their reported accuracy is claimed, reproduced, or cited as a
//! bar.** Table 1 is not used. Exactly one thing is borrowed: **Eq. (8)**, the
//! per-sample tangency target
//!
//! ```text
//! t_i = p_i + σ_i·|s_i|·(c_i − p_i) / ‖c_i − p_i‖
//! ```
//!
//! Every threshold in the registration is **this crate's own**: M-315's
//! placement ceiling (1.5%–21.5% of symmetric Hausdorff available to *any*
//! placement rule on `sphere` and `torus`, which is why C1 is asked only of the
//! sharp fields and C2 is a no-harm clause), and M-66's `35.796°` worst normal
//! error on `box_exact`, **identical at every resolution** — a corner does not
//! soften with `h`, so placement is the only thing that can move it.
//!
//! # The three arms run the crate's own pipeline, not a copy of it
//!
//! `isomesh::dual::VertexRule` is **public**, and `DualContouring::with_rule`'s
//! own doc-test names `Centroid` as *"the experiment's entry point"*. So the two
//! non-baseline arms are bench-local `impl VertexRule<f64>` blocks fed to the
//! crate's real [`DualContouring`], and nothing under `crates/isomesh/src` is
//! touched. That is stronger than reimplementing the pipeline here: the cell
//! classification, the bitmap active-cell prepass, the Hermite crossings, the
//! quad walk, the winding and the vertex emission are **the same machine code**
//! for all three arms, so a difference between two rows can only be the rule.
//!
//! - **`qef`** — `DualContouring::<f64>::new()`. Not a reconstruction of the
//!   baseline: it *is* the production constructor, `Qef { clamp: ToCell, lambda:
//!   None }`, so "does the baseline match `dual_contouring`" is not a question
//!   this file can get wrong.
//! - **`centroid`** — `isomesh::surface_nets::Centroid`, the crate's own rule,
//!   through the dual pipeline. The sanity anchor.
//! - **`tangency`** — below.
//!
//! ## The one thing that had to be transcribed, and how it is checked
//!
//! `dual_contouring::apply_clamp` is `pub(crate)`, so [`clamp_to_cell`] restates
//! it: the cell scaled about its own centre by `(1 − CLAMP_EPSILON)`, using the
//! crate's public constant. A transcription that has drifted would silently make
//! the tangency arm play by different rules from the baseline, so it is
//! **measured rather than asserted**: [`verify_clamp_transcription`] extracts
//! once with the crate's `Qef` (crate solve + crate clamp) and once with
//! [`LocalClampQef`] (crate solve + *this file's* clamp) and compares positions
//! and indices **bit for bit** through `f64::to_bits`. The result is on every CSV
//! row as `clamp_transcription_verified`. If that column is ever `false`, no
//! ratio in the file means anything.
//!
//! # The tangency rule, and the three places it had to make a decision
//!
//! Start from the cell's Hermite crossing centroid — `HermiteCell::centroid()`,
//! i.e. the answer the `centroid` arm gives — then apply Eq. (8) twice, clamping
//! after each. The tangency arm is therefore exactly *"the centroid arm plus two
//! Eq. (8) steps"*, which is what makes the delta between those two rows
//! attributable to the operator rather than to a different starting point.
//!
//! **1. `p_i` is the eight corner samples, not the crossings.** This is forced,
//! not chosen. A crossing has `s ≈ 0` by construction, so `|s_i| ≈ 0`, so
//! `t_i ≈ p_i`, so the least-squares point is the crossing centroid and the rule
//! is the `centroid` arm with extra arithmetic. Only the corner samples carry a
//! non-trivial radius, and the registration's own wording is *"the tangent points
//! its cell's samples imply"*. The corner values are already in hand — the engine
//! passes them to `place` — so the rule costs no extra field evaluation for them.
//!
//! **2. The residual is the full 3-vector, so the solve is a mean.** Eq. (9) of
//! the paper is `½ Σ ‖c_i(Ω) − t_i(Ω)‖²`, a Frobenius-norm fit; there is no
//! point-to-plane or dot-product residual anywhere in the paper. Minimising
//! `Σ‖v − t_i‖²` over a free `v` gives `v = mean(t_i)` — rank 3 per sample, no
//! matrix, no solve. A dot-product form would be a different operator with a
//! different null space and is not what was registered.
//!
//! **3. `σ_i` classifies `p_i` against the reconstruction — not `c_i` against the
//! field.** This is the one place the ticket's gloss and the source disagree, and
//! the source wins because the registration says *"using only Eq. (8) of Sellan,
//! Batty & Stein"*. The paper's two bullets under Eq. (8) read, verbatim:
//!
//! > where `σ_i` depends on the orientation of the surface `Ω` at `c_i(Ω)`:
//! > • If `p_i` is inside/outside `Ω` and the sign of `s_i` is negative/positive,
//! > then `σ_i = 1`. • If `p_i` is inside/outside `Ω` and the sign of `s_i` is
//! > positive/negative, then `σ_i = −1`.
//!
//! and, immediately after: *"We use the mesh element's normal vector at `c_i(Ω)`
//! to distinguish between inside and outside."* So the classified point is
//! **`p_i`**, the classifier is **the current reconstruction's normal at the
//! vertex**, and `σ_i = −1` is the *misclassification* case — §5.4 confirms the
//! geometry from the other side: *"merely by always making `σ_i = 1` in (8)…
//! this means we move the surface towards the closer of the sphere's two
//! possible tangent points"*.
//!
//! Transferred per cell: the local reconstruction is the plane through the
//! current vertex `c` with the cell's mean crossing normal `N`, and `p_i` is
//! inside it when `(p_i − c)·N < 0`. `N` is the equivariant mean of the Hermite
//! crossings' unit gradients — which is the second half of the Hermite data, and
//! the reason [`HermiteCell`] is the right type to build here rather than a bare
//! list of positions. It is deliberately **not** normalised: only the sign of the
//! dot product is used, and that is scale-invariant, so the division would be
//! arithmetic nobody reads.
//!
//! Where `N` is the zero vector the orientation carries no information, and the
//! paper's own Footnote 2 covers exactly that case — *"we do not distinguish
//! between inside and outside spheres if `c_i(Ω)` is on the boundary of a mesh
//! element, since normal information is not reliable there — `t_i(Ω)` is simply
//! the closest point on the sphere"* — i.e. `σ_i = +1`. That is a rule the source
//! states, not a fallback invented here.
//!
//! **What the ticket's gloss would have measured instead.** *"σ = +1 if the sign
//! of `s_i` agrees with the side `c_i` is on"* classifies `c`, by the sign of
//! `f(c)`. On a plane through a cell that operator does not have the surface as a
//! fixed point: starting at the crossing centroid with `f(c) = 0` (outside, by
//! the crate's sign rule), the four outside corners agree and the four inside
//! corners disagree, the mean lands at `−0.289h`, and the next iteration — now
//! reading `c` as inside — lands at `+0.258h`. It oscillates about the answer it
//! started on. Under Eq. (8) as the paper writes it the same cell is a **fixed
//! point**: every corner agrees with the plane through `c`, `σ ≡ +1`, the outside
//! targets sit at `+0.211h` and the inside ones at `−0.211h`, and the mean is `0`.
//!
//! ## Cost, and a line in the registration that is not quite right
//!
//! The registration budgets *"one normalize and one fma per sample"*. The
//! normalize is there; the `σ` test adds a three-component dot product per sample
//! per iteration, and building the `HermiteCell` costs one `Sdf::gradient` per
//! crossing — which the `qef` arm also pays, so the arms stay comparable, but the
//! registration's per-sample accounting omits both. `ns_per_sample` is recorded
//! and **gates nothing** (M-348, ✗24): every clause here is a distance ratio or a
//! count.
//!
//! # Determinism
//!
//! No map iteration, no PRNG, no threads, no clock in any gated quantity. The
//! eight per-axis addends of the least-squares mean and the twelve of the mean
//! normal are reduced by [`sum_equivariant`], which is `crate::equivariant`'s
//! algorithm transcribed — private, so it could not be called — summing ascending
//! by magnitude with the signed value as tie-break, so the result is a function
//! of the *set* of terms and a lattice rotation cannot move a vertex.
//!
//! # How error is measured, and why in two pieces
//!
//! `symmetric_hausdorff` is the crate's own `validate::accuracy`, which samples
//! **both** directions — mesh→field over vertices and centroids, field→mesh over
//! projected lattice seeds. That is the registered headline and it is not
//! recomputed here.
//!
//! `vertex_term` and `centroid_term` are the mesh→field half split in two, and
//! they are computed the way M-315 computed them in `placement_ceiling.rs`:
//! against `|f(p)|` directly. All four fields declare
//! [`FieldBound::Exact`](isomesh::fields::FieldBound::Exact), asserted per field
//! and recorded as `field_bound_exact`, so `|f(p)|` **is** the distance to the
//! surface and no projection is needed. `mesh_to_field_max` is carried beside them
//! so a reader can check `max(vertex_term, centroid_term)` against the crate's own
//! forward number in the same row.
//!
//! That split is C4, and C4 is the clause that pays whatever C1 does. M-315
//! measured Dual Contouring's Hausdorff **vertex-dominated on 8 of 8 rows**, and
//! its centroid error already **better than the perfect-placement floor by
//! 2.9–3.6×** on `sphere` — the QEF minimises distance to tangent planes and buys
//! well-centred facets at the cost of badly-placed vertices. A rule that pulls
//! vertices onto spheres spends that trade in reverse, so C4 predicts the vertex
//! term improving and the centroid term worsening. If the vertex term improves and
//! the centroid term does *not* worsen, M-315's trade is not what the QEF is doing
//! and a placement rule is not zero-sum — which is worth more than C1.
//!
//! ## Ratios point the same way everywhere
//!
//! Every `*_ratio_vs_qef` column is **`qef ÷ rule`**, so **above 1 is better than
//! the baseline** for all three of them and the C1 bar reads directly as
//! `hausdorff_ratio_vs_qef ≥ 1.25`. The `qef` rows are `1.0` by construction.
//!
//! ## Two columns the registration did not ask for, and why they are here
//!
//! A falsified clause with no mechanism is worth very little, so the artefact
//! carries the two numbers that decide *why*.
//!
//! `sigma_negative_fraction` answers whether the `σ` half of Eq. (8) fires at
//! all. If it is zero the row measured Eq. (8) with `σ ≡ 1`, which is a smaller
//! claim than the registered one and has to be legible from the file.
//!
//! `offset_from_cell_centre_in_cells` is the mean distance from a vertex to its
//! own cell's centre, in cells, and it names the operator's bias. Write Eq. (8)
//! as `t_i = p_i + λ_i·(c − p_i)` with `λ_i = σ_i|s_i| / ‖c − p_i‖`. If the eight
//! `λ_i` were equal, the least-squares mean would be exactly
//! `(1 − λ)·mean(p_i) + λ·c`, and **`mean(p_i)` over the eight corners of a cube
//! is the cell centre**. So the per-cell transfer of Eq. (8) is a contraction
//! toward the cell centre with a data-dependent strength — the price of
//! substituting one shared `c` for the paper's *per-sample* closest point
//! `c_i(Ω)`, which is the whole reason `t_i` lands on the surface for them and
//! need not here. That substitution is forced: a per-cell rule has one vertex and
//! no mesh to query, so there is no per-sample `c_i` to be had.

mod common;

use std::cell::Cell;
use std::time::Instant;

use common::experiment::Run;
use isomesh::dual::{CellVertices, VertexRule};
use isomesh::dual_contouring::solve::solve;
use isomesh::dual_contouring::{CLAMP_EPSILON, DualContouring};
use isomesh::fields::{BoxExact, FieldBound, ReferenceField, Sphere, ThinPlate, Torus};
use isomesh::hermite::HermiteCell;
use isomesh::surface_nets::Centroid;
use isomesh::validate::{
    AccuracyConfig, ValidateConfig, accuracy, self_intersections, validate_indexed,
};
use isomesh::{MeshBuffer, Real, RuntimeShape3, Sdf};

/// The registered resolution.
const SAMPLES: u32 = 65;

/// Eq. (8) applied exactly twice, as registered.
const ITERATIONS: u32 = 2;

/// Timing repetitions. The median is reported and gates nothing.
const REPS: usize = 5;

/// C1's bar, from M-66 rather than from the paper. See the module docs.
const SHARP_BAR: f64 = 1.25;

/// C2's tolerance, as a fraction of the `qef` value.
const SMOOTH_TOLERANCE: f64 = 0.10;

// ─── transcribed from the crate ─────────────────────────────────────────────

/// `dual_contouring::apply_clamp` with [`Clamp::ToCell`], which is `pub(crate)`.
///
/// The cell scaled about its own centre by `(1 − CLAMP_EPSILON)`, so the vertex
/// stays *strictly* interior. Checked bit for bit against the crate by
/// [`verify_clamp_transcription`] before any row is written.
fn clamp_to_cell<R: Real>(x: [R; 3], cell_origin: [R; 3], cell_size: R) -> [R; 3] {
    let half = cell_size * R::HALF;
    let inset = half * R::from_f64(1.0 - CLAMP_EPSILON);
    let mut out = x;
    for (axis, slot) in out.iter_mut().enumerate() {
        let centre = cell_origin[axis] + half;
        *slot = slot.clamp(centre - inset, centre + inset);
    }
    out
}

/// `crate::equivariant::precedes`: smaller magnitude first, then smaller value.
fn precedes<R: Real>(a: R, b: R) -> bool {
    match a.abs().total_cmp(&b.abs()) {
        core::cmp::Ordering::Less => true,
        core::cmp::Ordering::Greater => false,
        core::cmp::Ordering::Equal => a.total_cmp(&b) == core::cmp::Ordering::Less,
    }
}

/// `crate::equivariant::sum_equivariant`: sum smallest-magnitude-first, so the
/// result is a function of the *set* of terms rather than of the axis or corner
/// labelling a lattice rotation permutes.
fn sum_equivariant<R: Real, const N: usize>(mut t: [R; N]) -> R {
    let mut i = 1;
    while i < N {
        let mut j = i;
        while j > 0 && precedes(t[j], t[j - 1]) {
            t.swap(j - 1, j);
            j -= 1;
        }
        i += 1;
    }
    let mut acc = R::ZERO;
    for value in t {
        acc += value;
    }
    acc
}

/// `crate::cube::corner_offset`: corner `c` sits at `(c & 1, (c >> 1) & 1,
/// (c >> 2) & 1)`, bit 0 being x.
fn corner_offset<R: Real>(c: usize) -> [R; 3] {
    [
        R::from_f64(f64::from((c as u32) & 1)),
        R::from_f64(f64::from(((c as u32) >> 1) & 1)),
        R::from_f64(f64::from(((c as u32) >> 2) & 1)),
    ]
}

/// `crate::cube::is_inside`: the crate's one global decision about the sign.
/// Zero is **outside**.
fn is_inside<R: Real>(value: R) -> bool {
    value < R::ZERO
}

// ─── the rules ──────────────────────────────────────────────────────────────

/// The crate's QEF solve with **this file's** clamp. Verification only; never a
/// recorded arm.
#[derive(Clone, Copy, Debug)]
struct LocalClampQef;

impl<R: Real> VertexRule<R> for LocalClampQef {
    fn place<S: Sdf<Scalar = R>>(
        &self,
        sdf: &S,
        corner: &[R; 8],
        base: [u32; 3],
        origin: [R; 3],
        cell_size: R,
        out: &mut CellVertices<R>,
    ) {
        let cell_origin = cell_origin(base, origin, cell_size);
        let cell = HermiteCell::from_corners(sdf, corner, cell_origin, cell_size);
        let Some(x) = solve(&cell) else {
            return;
        };
        out.push_whole_cell(clamp_to_cell(x, cell_origin, cell_size));
    }
}

/// Eq. (8) of Sellán, Batty & Stein, iterated and clamped. See the module docs
/// for what this is and — at length — for what it is not.
///
/// Carries two counters so the artefact can answer *"was the `σ` half of Eq. (8)
/// exercised at all?"* — a rule whose `σ` is constant is a different, smaller
/// claim than the one that was registered, and that has to be visible in the CSV
/// rather than argued in a comment. They are updated once per cell, not once per
/// sample, so they cost nothing the timing column could see.
#[derive(Clone, Copy, Debug)]
struct Tangency<'a> {
    /// How many times Eq. (8) is applied. Registered as two.
    iterations: u32,
    /// `(corner, iteration)` pairs whose `σ` came out `−1`.
    sigma_negative: &'a Cell<u64>,
    /// `(corner, iteration)` pairs evaluated.
    sigma_samples: &'a Cell<u64>,
}

impl<R: Real> VertexRule<R> for Tangency<'_> {
    fn place<S: Sdf<Scalar = R>>(
        &self,
        sdf: &S,
        corner: &[R; 8],
        base: [u32; 3],
        origin: [R; 3],
        cell_size: R,
        out: &mut CellVertices<R>,
    ) {
        let cell_origin = cell_origin(base, origin, cell_size);
        let cell = HermiteCell::from_corners(sdf, corner, cell_origin, cell_size);
        // The engine only calls a rule on a cell whose corner signs disagree, and
        // the cube graph is connected, so some edge joins an inside corner to an
        // outside one and the cell is never empty. Handled rather than asserted
        // because `centroid` is the only source of the starting point.
        let Some(start) = cell.centroid() else {
            return;
        };

        // The local reconstruction's orientation: the equivariant mean of the
        // crossings' unit gradients, left unnormalised because only the sign of
        // `(p_i − c)·N` is read.
        let mut axes = [[R::ZERO; 12]; 3];
        for (edge, slot) in (0..12u8).zip(0..12usize) {
            let Some(crossing) = cell.get(edge) else {
                continue;
            };
            for (axis, value) in crossing.normal.into_iter().enumerate() {
                axes[axis][slot] = value;
            }
        }
        let normal = [
            sum_equivariant(axes[0]),
            sum_equivariant(axes[1]),
            sum_equivariant(axes[2]),
        ];
        // Footnote 2: with no reliable normal the paper does not distinguish
        // inside from outside and takes the closest point on the sphere.
        let oriented =
            normal[0] * normal[0] + normal[1] * normal[1] + normal[2] * normal[2] > R::ZERO;

        let mut corner_position = [[R::ZERO; 3]; 8];
        for (c, slot) in corner_position.iter_mut().enumerate() {
            let offset = corner_offset::<R>(c);
            for (axis, coordinate) in slot.iter_mut().enumerate() {
                *coordinate = cell_origin[axis] + cell_size * offset[axis];
            }
        }

        let mut vertex = clamp_to_cell(start, cell_origin, cell_size);
        let mut negatives = 0u64;
        for _ in 0..self.iterations {
            let mut target = [[R::ZERO; 8]; 3];
            for (c, position) in corner_position.iter().enumerate() {
                let s = corner[c];
                let d = [
                    vertex[0] - position[0],
                    vertex[1] - position[1],
                    vertex[2] - position[2],
                ];
                let length = (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt();
                // The clamp keeps the vertex strictly inside its cell, so it can
                // never coincide with a corner and the division below cannot be
                // by zero. Written as a positive test anyway: a degenerate cell
                // would otherwise write a NaN into the mesh, where it is
                // invisible until something downstream reads it. A corner the
                // vertex sits on has no direction, so its target is itself.
                if length > R::ZERO {
                    // σ: is `p_i` on the side of the local reconstruction that
                    // its own sign says it should be on?
                    //
                    // `d` runs from the corner to the vertex, so the corner is on
                    // the negative side of the plane through `c` when
                    // `(p_i − c)·N < 0`, i.e. when `d·N > 0`.
                    let side = d[0] * normal[0] + d[1] * normal[1] + d[2] * normal[2];
                    let reconstruction_inside = side > R::ZERO;
                    let agrees = reconstruction_inside == is_inside(s);
                    let sigma = if !oriented || agrees {
                        R::ONE
                    } else {
                        negatives += 1;
                        -R::ONE
                    };
                    let step = sigma * s.abs() / length;
                    for (axis, slot) in target.iter_mut().enumerate() {
                        slot[c] = position[axis] + d[axis] * step;
                    }
                } else {
                    for (axis, slot) in target.iter_mut().enumerate() {
                        slot[c] = position[axis];
                    }
                }
            }
            // `argmin_v Σ‖v − t_i‖²` over a free `v` is `mean(t_i)`.
            let inverse = R::from_f64(8.0).recip();
            let mean = [
                sum_equivariant(target[0]) * inverse,
                sum_equivariant(target[1]) * inverse,
                sum_equivariant(target[2]) * inverse,
            ];
            vertex = clamp_to_cell(mean, cell_origin, cell_size);
        }
        self.sigma_negative
            .set(self.sigma_negative.get() + negatives);
        self.sigma_samples
            .set(self.sigma_samples.get() + 8 * u64::from(self.iterations));

        out.push_whole_cell(vertex);
    }
}

/// World position of a cell's corner 0.
fn cell_origin<R: Real>(base: [u32; 3], origin: [R; 3], cell_size: R) -> [R; 3] {
    [
        origin[0] + cell_size * R::from_f64(f64::from(base[0])),
        origin[1] + cell_size * R::from_f64(f64::from(base[1])),
        origin[2] + cell_size * R::from_f64(f64::from(base[2])),
    ]
}

// ─── measurement ────────────────────────────────────────────────────────────

/// Extract one mesh, or fail loudly. A row built on a failed extraction is a row
/// about nothing.
fn extract<F, V>(
    mesher: &mut DualContouring<f64, V>,
    field: &F,
    shape: &RuntimeShape3,
    origin: [f64; 3],
    cell_size: f64,
) -> MeshBuffer<f64>
where
    F: Sdf<Scalar = f64>,
    V: VertexRule<f64>,
{
    let mut out = MeshBuffer::<f64>::new();
    let Ok(()) = mesher.extract(field, shape, origin, cell_size, &mut out) else {
        panic!("P-52: extraction failed at {SAMPLES}³");
    };
    out
}

/// Median seconds-per-sample of `REPS` extractions, in nanoseconds.
///
/// Recorded because it is interesting. Gates nothing — see the module docs.
fn time_ns_per_sample<F, V>(
    mesher: &mut DualContouring<f64, V>,
    field: &F,
    shape: &RuntimeShape3,
    origin: [f64; 3],
    cell_size: f64,
) -> f64
where
    F: Sdf<Scalar = f64>,
    V: VertexRule<f64>,
{
    let samples = f64::from(SAMPLES) * f64::from(SAMPLES) * f64::from(SAMPLES);
    let mut out = MeshBuffer::<f64>::new();
    let mut runs = Vec::with_capacity(REPS);
    for _ in 0..REPS {
        out.reset();
        let start = Instant::now();
        let ok = mesher
            .extract(field, shape, origin, cell_size, &mut out)
            .is_ok();
        let elapsed = start.elapsed();
        assert!(ok, "P-52: extraction failed while timing");
        std::hint::black_box(&out);
        runs.push(elapsed.as_secs_f64() * 1e9 / samples);
    }
    runs.sort_by(f64::total_cmp);
    runs[REPS / 2]
}

/// Every coordinate of both meshes identical to the last bit.
fn bit_identical(a: &MeshBuffer<f64>, b: &MeshBuffer<f64>) -> bool {
    a.indices == b.indices
        && a.positions.len() == b.positions.len()
        && a.positions
            .iter()
            .zip(&b.positions)
            .all(|(p, q)| p.iter().zip(q).all(|(x, y)| x.to_bits() == y.to_bits()))
}

/// Triangle centroids, skipping faces with an unusable index — the same filter
/// `placement_ceiling.rs` applies, so the two are comparable.
fn centroids(positions: &[[f64; 3]], indices: &[u32]) -> Vec<[f64; 3]> {
    let mut out = Vec::with_capacity(indices.len() / 3);
    for t in indices.chunks_exact(3) {
        let (Some(a), Some(b), Some(c)) = (
            positions.get(t[0] as usize),
            positions.get(t[1] as usize),
            positions.get(t[2] as usize),
        ) else {
            continue;
        };
        out.push([
            (a[0] + b[0] + c[0]) / 3.0,
            (a[1] + b[1] + c[1]) / 3.0,
            (a[2] + b[2] + c[2]) / 3.0,
        ]);
    }
    out
}

/// The worst of a set of distances, or zero for an empty set.
fn worst(values: impl Iterator<Item = f64>) -> f64 {
    values.fold(0.0_f64, f64::max)
}

/// Mean distance from a vertex to the centre of its own cell, in cells.
///
/// The diagnostic that names the mechanism rather than only its size. Eq. (8)
/// reads `t_i = p_i + λ_i(c − p_i)` with `λ_i = σ_i|s_i|/‖c − p_i‖`, so if the
/// eight `λ_i` were equal the least-squares mean would be
/// `(1 − λ)·mean(p_i) + λ·c` — and `mean(p_i)` over the eight corners of a cube
/// **is the cell centre**. The operator is therefore a contraction of the vertex
/// toward the cell centre by a factor it computes from the samples, and this
/// column is what makes that visible instead of inferred.
///
/// Nearest cell centre rather than `floor`: every rule here keeps its vertex
/// inside its own cell, so the nearest centre is that cell's, and a vertex
/// sitting exactly on a shared face — which `box_exact`'s lattice-aligned
/// crossings do produce for the `centroid` arm — gets the same answer from
/// either neighbour instead of an arbitrary one.
fn mean_offset_from_cell_centre(positions: &[[f64; 3]], origin: [f64; 3], cell_size: f64) -> f64 {
    if positions.is_empty() {
        return 0.0;
    }
    let mut sum = 0.0_f64;
    for p in positions {
        let mut squared = 0.0_f64;
        for axis in 0..3 {
            let index = ((p[axis] - origin[axis]) / cell_size - 0.5).round();
            let centre = origin[axis] + cell_size * (index + 0.5);
            let offset = p[axis] - centre;
            squared += offset * offset;
        }
        sum += squared.sqrt();
    }
    sum / positions.len() as f64 / cell_size
}

/// Everything one row reports about one mesh.
struct Metrics {
    symmetric_hausdorff: f64,
    mesh_to_field_max: f64,
    vertex_term: f64,
    centroid_term: f64,
    self_intersections_per_1k: f64,
    vertices: u64,
    triangles: u64,
    non_manifold_edges: u64,
    offset_from_cell_centre: f64,
}

fn measure<F>(
    field: &F,
    mesh: &MeshBuffer<f64>,
    shape: &RuntimeShape3,
    origin: [f64; 3],
    cell_size: f64,
) -> Metrics
where
    F: Sdf<Scalar = f64>,
{
    let Ok(accuracy_config) = AccuracyConfig::from_cell_size(cell_size) else {
        panic!("P-52: {cell_size} is not a usable cell size");
    };
    let Ok(report) = accuracy(
        &mesh.positions,
        &mesh.indices,
        field,
        shape,
        origin,
        &accuracy_config,
    ) else {
        panic!("P-52: the accuracy harness rejected a {SAMPLES}³ mesh");
    };
    assert!(
        report.has_coverage(),
        "P-52: no accuracy coverage — the mesh missed the surface entirely"
    );

    let Ok(validate_config) = ValidateConfig::from_cell_size(cell_size) else {
        panic!("P-52: {cell_size} is not a usable cell size");
    };
    let validity = validate_indexed(&mesh.positions, &mesh.indices, &validate_config);

    let Ok(intersections) = self_intersections(&mesh.positions, &mesh.indices, cell_size) else {
        panic!("P-52: the self-intersection harness rejected a {SAMPLES}³ mesh");
    };

    // `|f(p)|` is the distance itself: every field here declares `FieldBound::Exact`,
    // which `sweep` asserts before this runs.
    let vertex_term = worst(mesh.positions.iter().map(|p| field.sample(*p).abs()));
    let centroid_term = worst(
        centroids(&mesh.positions, &mesh.indices)
            .iter()
            .map(|p| field.sample(*p).abs()),
    );

    let triangles = mesh.indices.len() as u64 / 3;
    let per_1k = if triangles == 0 {
        0.0
    } else {
        intersections.pairs.len() as f64 * 1000.0 / triangles as f64
    };

    Metrics {
        symmetric_hausdorff: report.symmetric_hausdorff(),
        mesh_to_field_max: report.mesh_to_field.max,
        vertex_term,
        centroid_term,
        self_intersections_per_1k: per_1k,
        vertices: mesh.positions.len() as u64,
        triangles,
        non_manifold_edges: validity.non_manifold_edges,
        offset_from_cell_centre: mean_offset_from_cell_centre(&mesh.positions, origin, cell_size),
    }
}

/// `qef ÷ rule`, so above one is better than the baseline.
///
/// A zero denominator would mean a rule landed every sample exactly on the
/// surface, which is not a thing that happens at `65³` and would be a defect in
/// the measurement rather than a result; it is reported as infinity so it cannot
/// be mistaken for a mild win.
fn ratio(qef: f64, rule: f64) -> f64 {
    if rule > 0.0 {
        qef / rule
    } else {
        f64::INFINITY
    }
}

// ─── the experiment ─────────────────────────────────────────────────────────

/// One row's worth of verdict material, kept for the closing summary.
struct Outcome {
    field: &'static str,
    sharp: bool,
    rule: &'static str,
    hausdorff_ratio: f64,
    vertex_ratio: f64,
    centroid_ratio: f64,
    counts_identical: bool,
}

/// The crate's `Qef` against the crate's solve with this file's clamp, bit for
/// bit. See the module docs.
fn verify_clamp_transcription<F>(
    field: &F,
    shape: &RuntimeShape3,
    origin: [f64; 3],
    cell_size: f64,
) -> bool
where
    F: Sdf<Scalar = f64>,
{
    let crate_side = extract(
        &mut DualContouring::<f64>::new(),
        field,
        shape,
        origin,
        cell_size,
    );
    let bench_side = extract(
        &mut DualContouring::<f64, LocalClampQef>::with_rule(LocalClampQef),
        field,
        shape,
        origin,
        cell_size,
    );
    bit_identical(&crate_side, &bench_side)
}

fn sweep<F>(run: &mut Run, field: &F, sharp: bool, outcomes: &mut Vec<Outcome>)
where
    F: ReferenceField + Sdf<Scalar = f64>,
{
    let (shape, origin, cell_size) = common::grid::<f64, F>(field, SAMPLES);
    // The vertex and centroid terms read `|f(p)|` as a distance, which is only
    // true for an exact field. Asserted rather than assumed, and carried onto
    // every row.
    let exact = field.bound() == FieldBound::Exact;
    assert!(
        exact,
        "P-52: {} does not publish an exact distance, so |f| is not a distance",
        F::NAME
    );
    let clamp_ok = verify_clamp_transcription(field, &shape, origin, cell_size);

    let qef_mesh = extract(
        &mut DualContouring::<f64>::new(),
        field,
        &shape,
        origin,
        cell_size,
    );
    let centroid_mesh = extract(
        &mut DualContouring::<f64, Centroid>::with_rule(Centroid),
        field,
        &shape,
        origin,
        cell_size,
    );
    // Fresh counters for the recorded extraction only, so the fraction describes
    // the mesh in the CSV and not that mesh plus five timing repeats.
    let sigma_negative = Cell::new(0u64);
    let sigma_samples = Cell::new(0u64);
    let tangency_mesh = extract(
        &mut DualContouring::<f64, Tangency>::with_rule(Tangency {
            iterations: ITERATIONS,
            sigma_negative: &sigma_negative,
            sigma_samples: &sigma_samples,
        }),
        field,
        &shape,
        origin,
        cell_size,
    );
    let sigma_negative_fraction = if sigma_samples.get() == 0 {
        0.0
    } else {
        sigma_negative.get() as f64 / sigma_samples.get() as f64
    };

    let qef_ns = time_ns_per_sample(
        &mut DualContouring::<f64>::new(),
        field,
        &shape,
        origin,
        cell_size,
    );
    let centroid_ns = time_ns_per_sample(
        &mut DualContouring::<f64, Centroid>::with_rule(Centroid),
        field,
        &shape,
        origin,
        cell_size,
    );
    let timing_negative = Cell::new(0u64);
    let timing_samples = Cell::new(0u64);
    let tangency_ns = time_ns_per_sample(
        &mut DualContouring::<f64, Tangency>::with_rule(Tangency {
            iterations: ITERATIONS,
            sigma_negative: &timing_negative,
            sigma_samples: &timing_samples,
        }),
        field,
        &shape,
        origin,
        cell_size,
    );

    let qef = measure(field, &qef_mesh, &shape, origin, cell_size);
    let arms = [
        ("qef", 0u32, &qef_mesh, qef_ns),
        ("centroid", 0, &centroid_mesh, centroid_ns),
        ("tangency", ITERATIONS, &tangency_mesh, tangency_ns),
    ];

    for (rule, iterations, mesh, ns) in arms {
        let m = measure(field, mesh, &shape, origin, cell_size);
        let counts_identical = m.vertices == qef.vertices
            && m.triangles == qef.triangles
            && m.non_manifold_edges == qef.non_manifold_edges;
        let hausdorff_ratio = ratio(qef.symmetric_hausdorff, m.symmetric_hausdorff);
        let vertex_ratio = ratio(qef.vertex_term, m.vertex_term);
        let centroid_ratio = ratio(qef.centroid_term, m.centroid_term);

        println!(
            "{:<11} {rule:<9} H {:>10.4e} ×{hausdorff_ratio:>6.3}   \
             vtx {:>10.4e} ×{vertex_ratio:>6.3}   cen {:>10.4e} ×{centroid_ratio:>6.3}   \
             v {:>6} t {:>6} nme {:>4} same {counts_identical}  off {:>5.3} cells",
            F::NAME,
            m.symmetric_hausdorff,
            m.vertex_term,
            m.centroid_term,
            m.vertices,
            m.triangles,
            m.non_manifold_edges,
            m.offset_from_cell_centre,
        );

        run.record(&[
            ("field", F::NAME.to_string()),
            ("samples_per_axis", SAMPLES.to_string()),
            ("rule", rule.to_string()),
            ("iterations", iterations.to_string()),
            (
                "symmetric_hausdorff",
                format!("{:.9e}", m.symmetric_hausdorff),
            ),
            ("hausdorff_ratio_vs_qef", format!("{hausdorff_ratio:.6}")),
            ("vertex_term", format!("{:.9e}", m.vertex_term)),
            ("vertex_term_ratio_vs_qef", format!("{vertex_ratio:.6}")),
            ("centroid_term", format!("{:.9e}", m.centroid_term)),
            ("centroid_term_ratio_vs_qef", format!("{centroid_ratio:.6}")),
            (
                "self_intersections_per_1k",
                format!("{:.6}", m.self_intersections_per_1k),
            ),
            ("vertices", m.vertices.to_string()),
            ("triangles", m.triangles.to_string()),
            ("non_manifold_edges", m.non_manifold_edges.to_string()),
            ("counts_identical_to_qef", counts_identical.to_string()),
            ("ns_per_sample", format!("{ns:.4}")),
            // Extra columns. The first four are provenance for things a reader
            // would otherwise have to take on trust; the last two are the
            // mechanism, and are the reason a falsified C1 is still evidence.
            ("clamp_transcription_verified", clamp_ok.to_string()),
            ("field_bound_exact", exact.to_string()),
            ("mesh_to_field_max", format!("{:.9e}", m.mesh_to_field_max)),
            ("cell_size", format!("{cell_size:.9}")),
            ("sharp_field", sharp.to_string()),
            (
                "offset_from_cell_centre_in_cells",
                format!("{:.6}", m.offset_from_cell_centre),
            ),
            (
                "sigma_negative_fraction",
                format!(
                    "{:.6}",
                    if rule == "tangency" {
                        sigma_negative_fraction
                    } else {
                        0.0
                    }
                ),
            ),
        ]);

        outcomes.push(Outcome {
            field: F::NAME,
            sharp,
            rule,
            hausdorff_ratio,
            vertex_ratio,
            centroid_ratio,
            counts_identical,
        });
    }
}

fn main() {
    let prereg = isomesh::experiment!("P-52");

    let mut outcomes: Vec<Outcome> = Vec::new();
    common::experiment::run(prereg, |run| {
        sweep(run, &Sphere::<f64>::canonical(), false, &mut outcomes);
        sweep(run, &Torus::<f64>::canonical(), false, &mut outcomes);
        sweep(run, &BoxExact::<f64>::canonical(), true, &mut outcomes);
        sweep(run, &ThinPlate::<f64>::canonical(), true, &mut outcomes);
    });

    let tangency: Vec<&Outcome> = outcomes.iter().filter(|o| o.rule == "tangency").collect();

    println!("\n--- clauses ---");

    let mut c1 = true;
    for o in tangency.iter().filter(|o| o.sharp) {
        let held = o.hausdorff_ratio >= SHARP_BAR;
        c1 &= held;
        println!(
            "C1 {:<11} hausdorff ×{:.4} vs {SHARP_BAR} — {}",
            o.field,
            o.hausdorff_ratio,
            if held { "HELD" } else { "FALSIFIED" }
        );
    }
    println!("C1 {}", if c1 { "HELD" } else { "FALSIFIED" });

    let mut c2 = true;
    for o in tangency.iter().filter(|o| !o.sharp) {
        // The column is `qef ÷ rule`; C2 is about how far the rule's value moved
        // from the baseline's, which is the other way up.
        let drift = (o.hausdorff_ratio.recip() - 1.0).abs();
        let held = drift <= SMOOTH_TOLERANCE;
        c2 &= held;
        println!(
            "C2 {:<11} |Δ| {:.4} vs {SMOOTH_TOLERANCE} — {}",
            o.field,
            drift,
            if held { "HELD" } else { "FALSIFIED" }
        );
    }
    println!("C2 {}", if c2 { "HELD" } else { "FALSIFIED" });

    let c3_failures = outcomes.iter().filter(|o| !o.counts_identical).count();
    println!(
        "C3 rows with counts differing from qef: {c3_failures} of {} — {}",
        outcomes.len(),
        if c3_failures == 0 {
            "HELD"
        } else {
            "FALSIFIED"
        }
    );

    let mut c4_fields = 0;
    for o in &tangency {
        let predicted = o.vertex_ratio > 1.0 && o.centroid_ratio < 1.0;
        c4_fields += usize::from(predicted);
        println!(
            "C4 {:<11} vertex ×{:.4} ({}), centroid ×{:.4} ({}) — {}",
            o.field,
            o.vertex_ratio,
            if o.vertex_ratio > 1.0 {
                "improved"
            } else {
                "not improved"
            },
            o.centroid_ratio,
            if o.centroid_ratio < 1.0 {
                "worsened"
            } else {
                "not worsened"
            },
            if predicted { "as predicted" } else { "no" }
        );
    }
    println!(
        "C4 {c4_fields} of 4 fields — {}",
        if c4_fields >= 3 { "HELD" } else { "FALSIFIED" }
    );
    println!("--- end clauses ---");
}
