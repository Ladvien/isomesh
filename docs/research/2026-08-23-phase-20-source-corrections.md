# What four reads of the primary literature found, before Phase 20 ran anything

**Date:** 2026-08-23 · **Method:** four parallel full-text reads of the sources
`docs/research/2026-08-23-findings-audit-and-phase-20-registrations.md` proposes to build on, done
*before* the registrations were written rather than after the harnesses failed.

Three of six registrations rested on a claim the source does not make. Each would have cost a harness to
discover — one of them would have produced a HELD verdict against a bar copied from a paper that never
published it, which is worse than a falsification because nothing would have looked wrong.

This is ✗21's family: *a property lifted from a summary, not from the thing itself*. The audit that
proposed these registrations caught four instances of it in the existing ledger. It then committed three
more in its own proposals. That is not a criticism of the audit so much as evidence about how the error
propagates: **every one of these came from a plausible, well-formed sentence about a paper, written by
someone who had not opened it.**

---

## The corrections, ranked by what they would have cost

### 1. The "13× on sharp features" figure does not exist — P-52

**Claimed** (audit §C, P-49): *Dual Contouring of Signed Distance Data* (`arXiv:2604.00157`) reports
*"edge Chamfer 0.0262 against MC 0.417 and DC 0.350 — 13× on sharp features"*, and a tangency residual of
the form `((t_j − q_j)·d_j)²`.

**Found.** The tangency literature in corpus is **Reach For the Spheres: Tangency-aware surface
reconstruction of SDFs** — Sellán, Batty & Stein, SIGGRAPH Asia 2023, `10.1145/3610548.3618196`, indexed
as `10.1145_3610548.3618196`. Read end to end:

- There is **no edge-Chamfer metric**, no sharp-feature metric, and no per-feature error table. The only
  metrics are Hausdorff, Chamfer and their own energy `E_φ`.
- There is **no 13× figure of any kind**. The only integer factors in the paper are the memory factors
  37 and 125, and an unquantified *"often by several integer factors"* for accuracy against Marching
  Cubes.
- Sharp features appear three times and **all three are qualitative** — Fig. 3's caption, a sentence
  about prior work, and a reference title.
- The residual is **not** a point-to-plane dot product. It is the full squared vector
  `‖c_i − t_i‖²` between the closest point on the mesh and the tangent point on the sphere, with
  `t_i = p_i + σ_i|s_i|(c_i − p_i)/‖c_i − p_i‖` (Eq. 8). With `t_i` frozen it is **rank-3** per sample and
  penalises tangential motion; the dot-product form is rank-1 and is a different, strictly weaker
  objective.

**And the reduction is contradicted by the source's own ablation.** The audit proposed a per-cell,
two-inner-iteration, clamped-to-cell variant. The paper's Fig. 2 caption is literally an argument against
the per-voxel formulation — it *"discards much of the available global information"* — Fig. 3 attributes
sharp-feature recovery to *"global information (not just per-voxel data)"*, and Fig. 17 measures clamping
away far spheres as **progressive loss of detail**. Their method is a global sparse solve with
per-iteration remeshing, run to a stopping rule that needs 10 non-improving iterations per coarse stage
and 100 at the finest, taking 1.1–9.4 seconds per shape in Python.

**Consequence.** P-52 keeps Eq. (8) — one normalize and one fma per sample, allocation-free and
`Real`-generic, the one genuinely portable piece — and states in the registration that this is **not
their algorithm**. Every bar is rebuilt from this crate's own measurements: M-315's placement ceiling and
M-66's `35.796°` corner error identical at every resolution. Table 1 is not cited anywhere.

**What this bought.** A registration that would have HELD or FAILED against a fabricated 13×, and either
way produced a ledger row citing a number no reader could find.

---

### 2. The paper reports no degenerate-triangle count — P-53

**Claimed** (audit §C, P-50): Custódio's third label removes MC33's coincident-vertex degeneracy, framed
so that a 10× reduction reads as reproducing a published result.

**Found.** `10.1186/s13173-019-0086-6` reports **no count of degenerate triangles removed, on any
dataset.** It reports radii-ratio histograms, Betti numbers and blocked-cube percentages, and **no
timings at all.**

Two structural findings that change what a bench can honestly do:

- **The label rule itself is a pure pre-pass** over the eight-corner classification — a ternary sign, one
  strict comparison per corner instead of one non-strict. That half reproduces exactly, bench-local, with
  no crate change. Good news, and it is why P-53 is feasible under this phase's boundaries.
- **The triangulator is not a lookup table at all.** The paper says outright the triangulation is
  generated *"without the need of a look up table"* — it is a per-cube convex hull, with a cross-cell
  dedup of already-triangulated grid faces and a grid-slice face-splitting preprocess inherited from
  Custódio 2013. Both are **cross-cell stateful** and cannot be a pre-pass.

**Consequence.** P-53 reproduces the label pre-pass plus coincident-vertex collapse and measures **this
crate's own** degenerate count against **this crate's own** baseline. The registration says plainly that
10× is our bar, not their figure, and that the convex-hull triangulator is out of scope.

---

### 3. No correlation argument, no tightening figure, and no min/max rule — P-54

**Claimed** (audit §C, P-51): affine arithmetic is *correlation-aware*, so `sin(a)·cos(b)` over a box gets
a tighter interval than `|∇|·r`, and the rejected-cell count on `gyroid` should rise from 688 to ≥1,400 —
a number derived from the 2× gap between `gyroid`'s declared `2√3` and M-267's measured `1.731`.

**Found**, reading `10.1016/j.cag.2010.07.003` in full:

- The Revised Affine Form is `x̂ = x₀ + Σᵢ xᵢεᵢ + eₓ[−1,1]` — for 3D, **five stored reals, fixed size,
  never growing**, because every non-affine error accumulates into the single trailing term. That is
  genuinely `no_std`- and no-heap-compatible, which is the best news in this document.
- **There is no correlation argument in the paper and no quantified tightening figure.** No *"x% fewer
  cells subdivided"*, no *"y× tighter bound"*. The only justification is a single sentence that interval
  methods *"are slow because of the interval overestimation"*. All evidence is end-to-end wall clock.
- **`min`, `max` and `abs` are not in the paper at all.** It handles set-theoretic union and intersection
  only through smooth R-function surrogates, and §4.4 concludes that general conditionals are an open
  problem.
- It needs no symbolic expression graph — it is forward abstract interpretation by operator overloading —
  **but** the field body must be generic over the arithmetic type and must contain no data-dependent
  branch, `min`/`max`, or `abs`.

**This last point turns a weakness into the experiment's sharpest clause.** This crate's fields split
exactly along that line. `gyroid` is `sin(x)cos(y) + sin(y)cos(z) + sin(z)cos(x)` — smooth trig with
shared arguments and **no min/max** — and it is precisely the field M-248 measured at 16.8% rejection.
`box_exact` and `csg_difference` are built from `min`/`max`, for which the only sound treatment collapses
the affine form to an interval and destroys every correlation.

**Consequence.** P-54's prediction is deliberately **non-uniform**: ≥1.5× on `gyroid`, **<5%** on the
min/max fields. A uniform prediction could not distinguish correlation from a merely better constant,
which is the defect the audit's own derivation had — it reached a plausible number by the wrong
mechanism. The registration states that the 1.5× is derived here from M-267, not reproduced from the
paper, and that C2's mechanism is our reasoning about the *absence* of a min/max rule.

**Also found:** **Spelunking the Deep** (Sharp & Jacobson, `10.1145/3528223.3530155`) **is** in corpus,
twice — the DOI copy and an arXiv duplicate `10.48550_arXiv.2202.02444`, the latter carrying page numbers
in its chunk metadata. The audit marked it `[acquire]`.

---

### 4. The theorem is real, and it is two-dimensional — P-55

**Claimed** (audit §C, P-52): `arXiv:2608.12142` supplies a monotone-edge topology certificate, with the
audit correctly noting the abstract's hedge — *"correct with regards to critical points in our
experiments"*.

**Found.** The citation **resolves exactly**: Finken, Li, Wang, Guo & Levine, *Topology-Preserving Meshing
of Implicit Scalar Fields via Monotonicity Constraints*, IEEE VIS 2026 short paper, submitted 12 Aug 2026.
One of four checked that needed no correction to its identity — and the audit's hedge was too generous to
itself in one direction and too harsh in another:

- **Theorem 1 is proved**, with a four-case proof, not merely observed. What is empirical is that the
  sampling-based refinement *achieves* monotonicity.
- **The proof is 2D**, and its Part-2b step is literally 2D-combinatorial: *"since a triangle has only
  three edges, at least two branches must intersect the same edge"*. A hexahedral cell has 12 edges and 6
  faces, and a 3D saddle's local level set is a surface, not four curve branches. The paper's own
  Discussion names the 2D restriction as a limitation.
- **It does not apply to the trilinear interpolant.** The theorem is about a PL function on a simplicial
  mesh where critical points can only occur at vertices. Under trilinear interpolation interior critical
  points genuinely exist — that is the origin of the ambiguous-face problem this crate has an
  asymptotic decider for.
- **There is no tolerance rule.** The predicate is a bare sign disagreement between sampled directional
  derivatives at `k = max(2, ⌈‖e‖/w⌉ + 1)` points. No epsilon, no relative tolerance, no flat-region
  guard, no noise model. On a flat region the derivative hovers at zero and the sign test fires on
  rounding.
- Their gradients come from autodiff on a neural field; this crate would use central differences, which
  changes the noise story in a way the paper never analyses.

**Consequence.** P-55 is registered as a **labelled 3D port**, never as a transported proof, and the
tolerance — `1e-12 · (|f(a)| + |f(b)|)`, recorded at `1e-14` and `1e-10` beside it — is stated in the
registration as isomesh's invention. A held C1 is evidence about this crate's meshes and nothing more.

---

## What did not need correcting

**P-51 and P-56 rest on measurements already in this ledger**, which is why they were registered first
and without waiting on a paper. P-51's bars come from M-27 and M-12; P-56's mechanism is M-283's
`(180° − θ)/2` applied a second time, and its expected order of magnitude comes from P-47's own artefact.
Neither borrows a number from outside.

That is the pattern worth keeping: **the two registrations that needed no external claim are the two that
could be written immediately, and the four that needed one all needed correcting.**

---

## The rule this suggests, offered rather than adopted

`scripts/doc_facts.sh` gates counts. F-009 proposes gating a `FINDINGS.md` table against the CSV it
names. Neither catches this class, because this class is a claim about a **document outside the
repository**.

The cheapest defence is procedural and this phase used it: **a registration that cites an external number
as a bar must name the corpus `doc_id` it was read from, and the read must happen before the registration
is written.** Three of six here failed that test and all three were caught in one parallel pass costing
four agent reads. The alternative — discovering it when the harness produces a verdict against a
fabricated threshold — costs a bench, a CSV, a ledger entry and a correction banner, and leaves a wrong
number in the git history forever.
