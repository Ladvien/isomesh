# Fifty registrations from eighteen fields of mathematics, and the identity found in our own source

**Date:** 2026-08-29 · **Repo state:** `466dd07`, 482 findings entries, 52 open tickets
**Corpus:** home-still, 9,537 documents, 291,094 chunks, both scribe instances healthy
**Companion:** `docs/research/2026-08-29-axes-and-vocabulary-v2.md` (the words)

**Method.** Eighteen areas of mathematics were probed against the corpus and against the repository
simultaneously, on the rule the 2026-08-23 sweep earned: *an area with a paper in the library and no
reference anywhere in the repo is unmined with the source already paid for*, and *absence from the
ledger is not evidence of absence from the library.* Five parallel probes ran roughly 200 semantic
queries and 40 word-boundary greps. Fourteen papers were acquired and converted during the session;
six are paywalled and are recorded as such rather than cited.

**The calibration, stated once.** On this corpus a query about a subject genuinely present scores
**0.65–0.72**; the noise floor for an absent subject is **0.54–0.60** and the top hit is visibly
unrelated. Every ABSENT verdict below is that pattern, and every one was cross-checked with a grep.

---

## 0. The headline, and it is not a paper

`crates/isomesh/src/marching_cubes/trilinear.rs:246` computes `b*b - 4*a*c` from coefficients built at
`:200-214`, attributed to Grosso and defended by `M-207`'s divergence from the textbook quadratic.

**That expression is identically Cayley's `2×2×2` hyperdeterminant of the eight corner values.**

Not proportional. Not numerically close. **Symbolically equal** — both are the same 12-term,
degree-4 homogeneous polynomial in `f₀…f₇`, and their difference simplifies to zero under the corner
indexing the crate already uses (`f[u + 2v + 4w]`, read off `twist_lo`/`du_lo`/`dv_lo`). Three further
properties were proved in the same check:

- **It is the pencil determinant.** `disc(det(A₀ + λA₁))` — where `A₀`, `A₁` are the two opposite-face
  `2×2` corner matrices — equals it exactly, and equals it for **all three axis pairings**
  (`0123|4567`, `0145|2367`, `0246|1357`). `M-206`'s "two independently derived constructions locate
  the same body saddles, to `1.1e-12`" is not a coincidence; it is two of the three slicings of one
  pencil.
- **Its sign is a `GL(2)³` invariant.** `Δ(g·A) = (det g₁ · det g₂ · det g₃)² · Δ(A)`, verified over 40
  random invertible triples in exact rationals. The weight is a **square**, so the sign — and therefore
  the body-saddle count — is invariant under every independent per-axis invertible reparametrisation.
  **Non-cubic cells cannot change it.** Nothing in the ledger currently asserts that.
- **It is invariant under all 48 octahedral relabellings and under negating the field**, absolutely
  rather than merely in sign. `✗39`, `✗49` and `M-177` establish equivariance facts of exactly this
  shape by empirical sweeps; here they are algebra.

**What the literature says about the object.** de Silva & Lim, *Tensor rank and the ill-posedness of
the best low-rank approximation problem*, arXiv:math/0607647 — acquired, converted (44 pp., 100 chunks)
and indexed during this session — devotes §5–§7 to real `2×2×2` classification. Verbatim, §6: *"the
rank of a tensor is 2 on the set {A | Det₂,₂,₂(A) > 0} and 3 on the set {A | Det₂,₂,₂(A) < 0}."*
Their Prop 5.9 is `Δ(A) > 0 ⟹ rank_⊗(A) ≤ 2`; Prop 5.10 is `rank_⊗(A) ≤ 2 ⟹ Δ(A) ≥ 0`. Their sign
convention (discriminant of `det(λ₁A₁ + λ₂A₂)`) agrees with the standard explicit Cayley normalisation
and therefore with ours.

**Two corrections that must travel with the identity, or it will be misused.** (i) `Δ = 0` is **not**
the degenerate stratum — de Silva & Lim's Prop 7.3 and remark record that tensors on `Δ = 0` are
*generically rank 3*, and `Δ = 0` carries ranks 0, 1, 2 **and** 3. The crate's `discriminant == 0`
branch therefore models a stratum whose structure it has not checked. (ii) `sign(Δ)` **alone does not
classify the orbit**; §7.2 needs the pair `(sign Δ, rank_⊞)` — multilinear rank as well — to separate
the eight real orbits. Any registration promising "the case taxonomy becomes a corollary" that cites
only the sign is overclaiming.

**And one thing the identity does *not* do.** The repo's quadratic solves for **intersections of the
face level-set hyperbolas**, which lie on the zero level set — *not* for `∇f = 0`. Verified: on a
random trilinear, the repo's roots give `|f| ≤ 5e-16` while the true critical points sit elsewhere. So
the separate BKK result below bounds a *different* two, and `SADDLE_COUNT = 2` is not derived from it.

---

## 1. What the eighteen probes found, in one table

| Area | Corpus | Repo cites | Verdict |
|---|---|---|---|
| Cayley hyperdeterminant / tensor orbits | ABSENT (0.542) → **acquired** | 0 | **Identity proved. Group A.** |
| BKK / mixed volume / Newton polytopes | ABSENT (0.608) | 0 | `MV = 2`, verified two ways. Group B |
| Viro patchworking / toric | ABSENT (0.564) → **acquired** | 0 | Hinges on Kuhn regularity. Group B |
| CAD / real quantifier elimination | ABSENT (0.573) | 1 line, **rejected on the wrong axis** | Group B |
| Minimal-surface topology (TPMS genus) | — | 0 | **χ oracle exists. Group C** |
| Metric-based anisotropic adaptation | **PRESENT, 0.713**, 5 papers | **0** | Biggest paid-for gap. Group D |
| Jones `β`-numbers / rectifiability | ABSENT (0.578) → **acquired** | 0 | Group E |
| Strang–Fix / quasi-interpolation | ABSENT (0.62 = the Blu paper already held) | 0 | Group F |
| Sparse grids / Smolyak | ABSENT (0.597) | 0 | Registered expecting a null. Group F |
| Information-based complexity | ABSENT (0.561) → **acquired** | 0 | The denominator. Group G |
| `N`-term / adaptive approximation classes | **PRESENT, 0.62–0.64**, 4 papers | 0 (Stevenson cited only for NVB closure) | Cheapest. Group G |
| Lattice quantization | ABSENT → **acquired** | 0 | Group H |
| Min rectangle partition | PRESENT but wrong objective → **acquired** Eppstein | ranked #8, never done | Group I |
| Analysis of Boolean functions | ABSENT (0.584) → **acquired** | 0 | Group J |
| PL continuation (Allgower & Schmidt) | root ABSENT, descendant present | 0 for the root | **1985 precedence.** Group K |
| Conley index | PRESENT (0.662) | 0 | Registered expecting death. Group K |
| Discrete Hodge / Helmholtz | PRESENT, never applied to Hermite data | **0 across the whole repo** | Group L |
| Varifolds; intrinsic triangulations | 0.637 / **PRESENT 0.677** | 0 | Group L |
| Nodal percolation (`ℓ_c > 0`, `d ≥ 3`) | ABSENT → **acquired** | 0 | Group M |

**Killed during the probe, so nobody re-searches them.** *Random matrix theory for QEF conditioning*
and *pseudospectra for the truncation threshold* — the QEF matrix `M = Σnᵢnᵢᵀ` is symmetric hence
**normal**, so its `ε`-pseudospectra are `ε`-balls around its eigenvalues and hand back exactly the
singular-value threshold they were meant to replace; RandNLA dies separately on `M` being `3×3`.
*Scheduling theory* — per-cell output is bounded, so span `T_∞ = O(1)` and Blumofe–Leiserson collapses
to the trivial bound. *Geometric separators for chunking* — an axis-aligned cube already achieves the
`Θ(n^{2/3})` isoperimetric optimum. *Ising/dimer/limit shapes* — the 2D results rest on bipartite
determinantal structure the 3D sign grid has not got. *Renormalization group for LOD* — no critical
point to define relevance against. *Compressed sensing* — deterministic grid, no fixed sparsity basis.
*Helly/LP-type certificates for the cell clamp* — the target is a 30-flop `3×3` adjugate; nothing
cheaper is worth writing.

**Corrections owed to earlier documents.** `2026-08-16-research-backlog-sota.md:5` says persistent
homology "returns zero"; that is true of `BACKLOG.md` and **false of the corpus**, which holds Dey, Hou
& Morozov on zigzag updates and Dey, Fan & Wang on simplicial-map persistence. The
`2026-08-26` audit rejected CAD as *"global, doubly exponential, needs arbitrary precision"* — correct
at runtime, and irrelevant to deriving a table **once, offline**. `2026-08-26-audit...md:365` flags
Etiene et al.'s DOI as unverified; it is **`10.1109/tvcg.2011.109`**. And `docs/measurements/resolution_sweep.csv`
is sphere-only **by construction** (`benches/resolution_sweep.rs:128` hard-codes `Sphere::canonical()`).

**Blocked acquisitions — paywalled, no open-access route found this session.** Loseille & Alauzet
Part I `10.1137/090754078` and Part II `10.1137/10078654X` (landing pages only; the framework is
restated verbatim in NASA NTRS 20200003084, which is what Group D cites); Binev–Dahmen–DeVore
`10.1007/s00211-003-0492-7`; Novak `10.1006/jcom.1996.0015`; Blu & Unser Part I `10.1109/78.790659`;
Nemhauser–Wolsey–Fisher `10.1007/BF01588971`; Boissonnat & Wintraecken `10.4230/LIPIcs.SoCG.2020.20`;
Allgower & Schmidt `10.1137/0722020`; Pöthkow et al. `10.1111/j.1467-8659.2011.01942.x`. **DOI trap:**
`10.1137/0722019` is *not* Allgower & Schmidt — it is Nanda on the QR algorithm, same volume and issue.

**Protocol.** Phase 15's protocol applies in full. All fifty `P-` entries are registered in
`crates/isomesh/src/experiment.rs` **before** any harness commit — that registration is the first
commit of the phase and is the only write to `crates/isomesh/src/**` that precedes a measurement.
`✗51`'s rule is applied: every clause stated as a ratio of a total carries a SHARE line, and every
registration carries a VACUITY CONTROL naming the column that proves the fixture could have failed.

**Nine registrations are expected to return nulls, and that is registered rather than hoped:**
`P-132`, `P-146`, `P-154`, `P-155` (whose null is the most valuable result in the phase), `P-159`,
`P-162`, `P-164`, `P-168`, `P-170`.

**Dependencies, stated once.** `P-123` runs before everything else in Group A. `P-127` rides `P-123`.
`P-138` runs before `P-139` and `P-141`. `P-142` runs before `P-143`, `P-145` and `P-147`. `P-151`
runs before `P-152` and `P-153`. `P-148` runs before `P-149` and `P-150`. Everything else is
independent, and eleven rows are `S`.

---

# Group A — the cell is a tensor (real algebraic geometry, tensor geometry, invariant theory)

#### P-123 — registered for R-123, before the harness: The body-saddle discriminant is Cayley's hyperdeterminant, and the proof is a test rather than a comment

**Ticket:** `R-123` (S). **Records:** `expression`, `terms_disc`, `terms_cayley`, `total_degree`, `symbolic_difference_is_zero`, `pencil_axis_pairing`, `pencil_matches`, `random_rational_trials`, `max_abs_ratio_deviation`, `f32_sign_disagreements`, `f64_sign_disagreements`, `c1_holds`, `c2_holds`, `c3_holds`.

**Hypothesis.** THE CRATE HAS SHIPPED A CLASSICAL INVARIANT FOR ITS ENTIRE LIFE AND CALLED IT A TRANSCRIBED QUADRATIC. (C1) `b*b - 4*a*c` from `BodySaddles::coefficients` is **identically** Cayley's `2×2×2` hyperdeterminant of `[f₀..f₇]` under the indexing `f[u + 2v + 4w]`: same 12 terms, total degree 4, symbolic difference exactly zero — asserted by a committed symbolic check, not by a numeric sample. (C2) The same polynomial equals `disc(det(A₀ + λA₁))` for the two opposite-face `2×2` corner matrices, and does so for **all three** axis pairings, which is the mechanism behind `M-206`. (C3) In exact rational arithmetic over at least 3,000 random 8-tuples the two expressions agree with ratio exactly 1 and zero sign disagreements; in `f32` they do **not**, and the disagreement count is the number `P-130` exists to act on. SHARE: C1 and C2 move no runtime cost at all — this row is a renaming plus two theorems, and its whole value is what Group A's other rows can then assume.

**Falsified by.** C1 by any non-zero symbolic difference, which would mean the indexing convention was misread and every downstream row in Group A collapses. C2 by any axis pairing disagreeing, which would make `M-206` a coincidence after all. C3 by zero `f32` disagreements, which would retire `P-130` before it is written. VACUITY CONTROL: the trial set must contain 8-tuples of both discriminant signs and at least 50 within `1e-6` of zero, reported as counts, or C3 is sampling only the easy stratum.

#### P-124 — registered for R-124, before the harness: The saddle count cannot depend on cell aspect ratio, and this is the first thing in the ledger that says so

**Ticket:** `R-124` (S). **Records:** `field`, `resolution`, `cell_aspect_ratio`, `axis_scales`, `saddle_count_isotropic`, `saddle_count_anisotropic`, `count_disagreements`, `sign_delta_isotropic`, `sign_delta_anisotropic`, `sign_disagreements`, `gl2_weight_check`, `c1_holds`, `c2_holds`.

**Hypothesis.** `Δ(g·A) = (det g₁ det g₂ det g₃)²·Δ(A)` — the weight is a perfect square, so `sign(Δ)` is an **absolute** invariant of the `GL(2)³` action, and the body-saddle count is a property of the eight values alone. (C1) On all eight reference fields at `33³` and `65³`, meshing with per-axis cell scales of `(1,1,1)`, `(1,2,4)` and `(1,1,8)` produces **bit-identical** body-saddle counts per cell, with zero disagreements. (C2) The relative-invariant weight holds numerically: for 500 random invertible `g` triples, `Δ(g·A)/Δ(A)` equals `(det g₁ det g₂ det g₃)²` to `f64` rounding. SHARE: C1 is a correctness clause and moves no cost; it removes an untested assumption from every anisotropic-cell path, including all of Group D.

**Falsified by.** C1 by any disagreement, which would mean the extractor's saddle path depends on geometry the invariant says it cannot — i.e. a bug, and a valuable one. C2 by a weight that is not a square, which would mean the indexing or the action is wrong and `P-123` C1 must be re-read. VACUITY CONTROL: the anisotropic arms must produce a non-zero count of cells whose *vertex positions* move, reported as a column, or the fixture is not exercising anisotropy at all.

#### P-125 — registered for R-125, before the harness: Octahedral and negation invariance as algebra, at a fraction of the cost of the 48-element sweep

**Ticket:** `R-125` (S). **Records:** `element`, `sweep_ms`, `invariant_ms`, `speedup`, `cells_checked`, `sweep_violations`, `invariant_violations`, `agreement`, `negation_violations`, `c1_holds`, `c2_holds`, `c3_holds`.

**Hypothesis.** `✗49` established bit-exact octahedral equivariance of plain Marching Cubes by sweeping all 48 elements; `Δ` is invariant under all 48 **and** under negating the field, provably. (C1) A `Δ`-based check agrees with the 48-element sweep on every cell of every reference field — zero disagreements — while costing at least `10×` less. (C2) `Δ` is invariant under field negation on every cell, which `M-177` says reordering cannot buy for vertex placement; the contrast localises `M-177`'s obstruction to the *solve*, not to the *classification*. (C3) The invariant check catches at least one configuration class the sweep is structurally blind to, or reports honestly that it catches none. SHARE: C1 moves the equivariance-test stage of CI, currently the 48-element sweep in full.

**Falsified by.** C1 by any disagreement or under `10×`. C2 by a negation violation, which would contradict the symbolic proof and send `P-123` back. C3 by "none", which is a fine outcome and should be recorded as such rather than hunted. VACUITY CONTROL: the sweep arm must report a non-zero violation count on a deliberately broken table, or both instruments are agreeing on nothing.

#### P-126 — registered for R-126, before the harness: `Δ > 0` means real tensor rank 2, and the census says whether rank predicts anything

**Ticket:** `R-126` (M). **Records:** `field`, `resolution`, `cells`, `delta_positive`, `delta_negative`, `delta_zero`, `rank_two_cells`, `rank_three_cells`, `ambiguous_cells`, `rank_vs_case_index_agreement`, `mutual_information`, `c1_holds`, `c2_holds`, `c3_holds`.

**Hypothesis.** de Silva & Lim §6: rank is 2 on `{Δ > 0}` and 3 on `{Δ < 0}`. A cell's real tensor rank is therefore computable from a sign the crate already computes. (C1) Across eight fields at three resolutions, the `(Δ > 0, Δ < 0, Δ = 0)` partition is stable and `Δ = 0` is rare — under 0.1% of surface cells at `f64`. (C2) Real tensor rank carries information about ambiguity that the 8-bit case index does not: mutual information between rank and "this cell is ambiguous" exceeds that between rank and the case index by a stated margin. (C3) On `sphere`, where no cell is ambiguous, rank is 2 on essentially every surface cell — the control that says the statistic is not measuring noise. SHARE: C1 and C2 are counts, not timings, and move nothing; they decide whether Group A continues past `P-127`.

**Falsified by.** C1 by an unstable partition or `Δ = 0` above 0.1%, the latter meaning the degenerate stratum is common enough that `P-127` becomes urgent rather than optional. C2 by rank being a function of the case index, which would make it a renaming with no new signal — the most likely outcome and worth knowing cheaply. C3 by `sphere` showing rank-3 cells, which would say the statistic is dominated by arithmetic rather than geometry. VACUITY CONTROL: `gyroid` and `csg_difference` must contribute a non-zero ambiguous-cell count, or C2's mutual information is computed against a constant.

#### P-127 — registered for R-127, before the harness: `Δ = 0` is not the degenerate case, and the branch that assumes it is has never been audited

**Ticket:** `R-127` (M). **Records:** `field`, `resolution`, `discriminant_zero_hits`, `a_zero_hits`, `double_root_hits`, `border_rank_two`, `true_rank_three`, `w_state_like`, `roots_reported`, `roots_true`, `mesh_delta_triangles`, `c1_holds`, `c2_holds`, `c3_holds`.

**Hypothesis.** `roots()` treats `discriminant == 0` as "a double root is one intersection point", which is correct arithmetic about the quadratic and **may be the wrong model of the configuration**: de Silva & Lim Prop 7.3 records that tensors on `Δ = 0` are *generically rank 3*, with border rank 2 — the W-state is the canonical example, `Δ = 0` exactly with real rank 3. (C1) The `discriminant == 0` branch fires on at least one reference field at `f64`, and the count is reported per field rather than assumed zero. (C2) On the cells where it fires, the true configuration is classified — rank 1, 2, or 3 — and at least one is **not** the tangential-touch case the comment describes. (C3) Correcting the classification changes at least one triangle on at least one field, or provably changes none, and which of those is true is stated as the result. SHARE: C3 moves only the cells C1 counts, which is why C1 is reported first and separately.

**Falsified by.** C1 by the branch never firing on any field at any resolution, which retires the row and is worth the afternoon. C2 by every hit being a genuine tangency, which would vindicate the comment. C3 by no mesh change, which would make this a documentation fix and should be recorded as one. VACUITY CONTROL: a synthetic W-state cell must be added to the fixture and must reach the branch, or C2 cannot distinguish "no such cells exist" from "the fixture has none".

#### P-128 — registered for R-128, before the harness: The eight real orbits, and whether orbit class predicts the defects the validity suite finds

**Ticket:** `R-128` (M). **Records:** `field`, `resolution`, `orbit_class`, `sign_delta`, `multilinear_rank`, `cells_in_orbit`, `nonmanifold_edges_in_orbit`, `self_intersections_in_orbit`, `orphaned_vertices_in_orbit`, `defect_rate_per_orbit`, `chi_square`, `c1_holds`, `c2_holds`, `c3_holds`.

**Hypothesis.** de Silva & Lim §7.2: orbit class is determined by the pair `(sign Δ, rank_⊞)` — **the sign alone is not enough**, and any claim that the taxonomy follows from `Δ` must carry the multilinear rank too. (C1) All eight real orbits are populated across the reference fields, or the unpopulated ones are named. (C2) T-001 defect rates differ significantly across orbits — chi-square against the null that defects are orbit-independent. (C3) The orbit partition is **not** a relabelling of the 256-case index: at least one case index spans two orbits and at least one orbit spans two case indices. SHARE: C2 is a rate per orbit, not a share of a total; the population counts are reported so the rates have denominators.

**Falsified by.** C1 by orbits that no field reaches, which narrows the taxonomy usefully. C2 by orbit-independent defects, which would mean the classification is mathematically real and operationally inert — the honest and quite likely outcome. C3 by a bijection with the case index, which would make the whole orbit apparatus a renaming. VACUITY CONTROL: the defect columns must be non-zero in total, asserted from T-001's own counts, or C2 is a chi-square on an empty table.

#### P-129 — registered for R-129, before the harness: An exact sign for the hyperdeterminant, and how often `f32` gets it wrong

**Ticket:** `R-129` (M). **Records:** `field`, `resolution`, `scalar`, `cells`, `sign_disagreements_f32`, `sign_disagreements_f64`, `disagreement_rate`, `exact_ms`, `float_ms`, `overhead_ratio`, `filtered_fallback_rate`, `triangles_changed`, `nonmanifold_delta`, `c1_holds`, `c2_holds`, `c3_holds`.

**Hypothesis.** `Δ` is a degree-4 polynomial in 8 inputs — squarely in the range where Shewchuk's adaptive-precision expansions (`10.1007/pl00009321`, already mined at 5 files) give an exactly-correct sign, and where a floating-point filter takes the exact path on a small minority of cells. (C1) At `f32`, the sign of `Δ` disagrees with the exact sign on a non-zero fraction of surface cells on at least two reference fields. (C2) A filtered exact predicate costs under `1.5×` the naive float evaluation in aggregate, because the filter succeeds on the overwhelming majority of cells. (C3) Correcting the sign changes the mesh — a non-zero triangle delta, or a non-zero change in T-001's non-manifold-edge count — on at least one field. SHARE: C2 moves the body-saddle stage only, whose share of extraction must be reported alongside.

**Falsified by.** C1 by zero disagreements at `f32`, which would say the predicate is not near-degenerate on real fields and retire C2 and C3 with it. C2 by above `1.5×`. C3 by no mesh change, which would mean the sign errors fall on cells whose classification does not reach the output — interesting, and a reason to look at where they do fall. VACUITY CONTROL: at least one field must produce cells with `|Δ|` below the `f32` error bound, reported as a count, or C1 is measuring a problem the fixture excludes.

#### P-130 — registered for R-130, before the harness: `|Δ|` normalised as a continuous ambiguity magnitude, instead of a binary test

**Ticket:** `R-130` (S). **Records:** `field`, `resolution`, `normalisation`, `delta_magnitude`, `scale_invariance_error`, `rank_correlation_with_defects`, `rank_correlation_with_hausdorff`, `threshold_sweep`, `cells_above_threshold`, `c1_holds`, `c2_holds`.

**Hypothesis.** `Δ` is homogeneous of degree 4, so `|Δ| / (max|fᵢ|)⁴` is scale-free and is a candidate **magnitude** of ambiguity where the crate currently has only a sign — the same move persistence makes against the interior test, but computable in `O(1)` from data already in hand. (C1) The normalised magnitude is invariant under scaling all eight corner values, to `f64` rounding, and under the 48 octahedral elements. (C2) It rank-correlates with per-cell symmetric Hausdorff error above 0.5 on at least four of eight fields, making it a candidate refinement criterion rather than merely a classifier. SHARE: C2 offers a criterion for the adaptivity stage; it does not by itself move any cost.

**Falsified by.** C1 by scale dependence, which would mean the normalisation is wrong. C2 by correlation below 0.5 on five or more fields, which would say ambiguity magnitude and geometric error are different phenomena — likely, and it would close the row cleanly. VACUITY CONTROL: the per-cell Hausdorff column must have non-zero variance on every field, or the correlation is against a constant.

#### P-131 — registered for R-131, before the harness: Where the identity stops holding, which is the boundary of everything Group A claims

**Ticket:** `R-131` (S). **Records:** `reconstruction`, `is_multi_affine`, `identity_holds`, `symbolic_residual_terms`, `smooth_min_k`, `deviation_at_k`, `tricubic_degree`, `cases_touched`, `c1_holds`, `c2_holds`.

**Hypothesis.** Every claim in Group A rests on the reconstruction being **multi-affine**. That is a hypothesis about the filter, not a fact about the crate, and it must be recorded as a boundary before anything downstream leans on it. (C1) The identity holds for the trilinear and fails for any reconstruction that is not multi-affine — demonstrated symbolically for a tricubic and for `smooth_min(k)` with `k > 0`, with the residual reported rather than merely asserted non-zero. (C2) Under `smooth_min(k)`, the deviation from the multi-affine `Δ` grows with `k` and is bounded in the `O(k)` seam shell the 2026-08-23 memo measured — connecting Group A's boundary to `M-38`'s smoothing result rather than leaving it a separate fact. SHARE: none — this is a scope statement and moves nothing.

**Falsified by.** C1 by the identity surviving a non-multi-affine reconstruction, which would be a strictly better result and would widen Group A. C2 by deviation that does not track `k`, which would mean the smoothing and the invariant are unrelated and `M-38` gains nothing here. VACUITY CONTROL: the tricubic arm must produce a symbolic residual with named non-zero terms, not a numeric non-zero, or C1 has not established the failure structurally.

#### P-132 — registered for R-132, before the harness: A null registered on purpose — the hyperdeterminant does not resolve the interior test

**Ticket:** `R-132` (S). **Records:** `configuration_class`, `delta_sign`, `interior_test_verdict`, `agreement_rate`, `cases_where_delta_insufficient`, `chernyaev_disagreement_overlap`, `c1_holds`, `c2_holds`.

**Hypothesis.** THIS ROW EXISTS TO STOP AN OVERCLAIM BEFORE SOMEBODY MAKES IT. It is tempting to read Group A as "the 730-subcase table becomes a corollary of an invariant". It does not, and the registration says why in advance: `sign(Δ)` counts *hyperbola intersections in the plane of a pencil*, and the interior test asks *which side of the surface a tunnel opens on*, which is a different question with a different answer. (C1) Over the configurations `M-165` identifies (Chernyaev's numerator-only test wrong on 1,966 of 15,6xx), `sign(Δ)` alone agrees with the correct interior verdict on **strictly fewer than all** of them, and the disagreement set is characterised. (C2) The classes where `Δ` is insufficient are named explicitly, so the boundary is in the ledger rather than in somebody's head. SHARE: none — this is a negative result registered before its positive twin can be believed.

**Falsified by.** C1 by full agreement, which would be a major and unexpected win and would justify reopening the A-002 interior line entirely. C2 by an uncharacterisable disagreement set, which would mean the instrument is too coarse. VACUITY CONTROL: the configuration set must reproduce `M-165`'s 12.6% disagreement rate as a control column, or the fixture is not the population `M-165` measured.

---

# Group B — Newton polytopes, patchworking, and machine-derived tables (algebraic geometry, toric geometry, symbolic computation)

#### P-133 — registered for R-133, before the harness: A trilinear has at most two critical points, and that is a theorem rather than an assertion

**Ticket:** `R-133` (S). **Records:** `field`, `resolution`, `mixed_volume`, `critical_points_complex`, `critical_points_real`, `critical_points_in_cell`, `cells_with_zero`, `cells_with_one`, `cells_with_two`, `on_hyperplane_count`, `saddle_count_relationship`, `c1_holds`, `c2_holds`, `c3_holds`.

**Hypothesis.** `∂f/∂xᵢ` of a multi-affine `f` is multi-affine in the other two variables, so its Newton polytope is a unit square; the mixed volume of the three squares is `Vol([0,2]³) − 3·2 + 0 = 8 − 6 = 2`, verified by convex hull and independently by symbolic solution of `∇f = 0`. By Bernstein's theorem `f` has **at most 2** critical points in `(ℂ*)³`. (C1) The mixed volume computation is committed as a check, and a census over the reference fields reports how many cells have 0, 1 or 2 **real** critical points **inside** the unit cell — the bound is on neither. (C2) At least one cell on at least one field has a critical point on a coordinate hyperplane, which Bernstein's count excludes by construction, demonstrating the caveat rather than stating it. (C3) THE RELATIONSHIP TO `SADDLE_COUNT` IS STATED AS A NON-IDENTITY: the repo's quadratic solves for face-hyperbola intersections on the zero level set, not for `∇f = 0`, and the two 2s are shown to be different by exhibiting a cell where the root sets differ. SHARE: none — this replaces an authority citation with a derivation.

**Falsified by.** C1 by any cell with three or more real critical points, which would refute the mixed-volume computation. C2 by no hyperplane cases, which would make the caveat vacuous here and is worth recording. C3 by the two root sets coinciding, which would be a much stronger result and would mean `SADDLE_COUNT` *is* derivable from BKK. VACUITY CONTROL: the census must include cells of all three counts, reported as a distribution, or the bound is being checked against a single stratum.

#### P-134 — registered for R-134, before the harness: What the case explosion would cost at tricubic, priced before anyone proposes it

**Ticket:** `R-134` (M). **Records:** `reconstruction`, `degree_per_variable`, `mixed_volume`, `critical_points_bound`, `bernstein_route_alive`, `subdivision_exactness`, `estimated_case_count`, `c1_holds`, `c2_holds`.

**Hypothesis.** The 2026-08-23 memo killed Bernstein–Bézier form for the trilinear and said explicitly: *revisit the entire hypothesis if the crate ever adopts a tricubic reconstruction filter.* `P-153` may propose exactly that, so the price should be on the table first. (C1) The mixed volume for `∇f = 0` with a tricubic reconstruction is computed and reported, giving the critical-point bound that would replace the trilinear's 2. (C2) The Bernstein-form arguments the memo killed are re-checked at degree 3 per variable and reported as alive or dead, with the specific reason — the memo predicts alive, and a prediction registered before the check is worth more than one after. SHARE: none — this prices a decision, it does not make one.

**Falsified by.** C1 by a mixed volume that makes the tricubic case space computationally hopeless, which would be a strong argument against `P-153` and should be said plainly. C2 by the Bernstein arguments staying dead at degree 3, which would contradict the memo's own prediction. VACUITY CONTROL: the trilinear must be run through the same pipeline and reproduce `MV = 2`, or the tricubic number has no calibration.

#### P-135 — registered for R-135, before the harness: Is the Kuhn triangulation regular, and does marching tetrahedra therefore inherit a 1980s theorem

**Ticket:** `R-135` (M). **Records:** `triangulation`, `lifting_function`, `is_regular`, `secondary_polytope_vertex`, `patchworking_applies`, `isotopy_agreement`, `pv_disagreements`, `cells_checked`, `c1_holds`, `c2_holds`, `c3_holds`.

**Hypothesis.** Viro's combinatorial patchworking builds a real algebraic hypersurface from **signs at the lattice points of a Newton polytope**, glued over a triangulation, and states when the piecewise-linear result is isotopic to the true hypersurface. Signs at the eight corners of a cube, glued over the Kuhn/Freudenthal triangulation, is **marching tetrahedra**. The theorem's hypothesis is that the triangulation be **regular** (induced by a convex lifting). (C1) The Kuhn triangulation of the unit cube is regular, demonstrated by exhibiting an explicit convex lifting function, or shown not to be — this single fact decides the whole transfer. (C2) If regular, the patchworked surface is isotopic to a real algebraic hypersurface, and that isotopy claim agrees with T-015's Plantinga–Vegter checker on every cell of every reference field. (C3) The connection to the existing Lovász-extension framing is made explicit: the Lovász extension is PL interpolation over exactly this triangulation, so `2026-08-23`'s "marching tetrahedra contours the Lovász extension" and Viro's construction are the same statement from two fields. SHARE: none — this is a theorem acquisition.

**Falsified by.** C1 by non-regularity, which kills the transfer in an afternoon and is the outcome to hope for if the alternative is a long shot. C2 by disagreement with T-015, which would mean one of the two instruments is wrong and finding out which is worth more than the transfer. C3 by the two framings differing, which would be genuinely surprising. VACUITY CONTROL: T-015 must report at least one isotopy failure somewhere in the corpus of test cells, or C2's agreement is agreement on a constant.

#### P-136 — registered for R-136, before the harness: Fields whose topology is known by construction, generated rather than measured

**Ticket:** `R-136` (M). **Records:** `construction`, `prescribed_genus`, `prescribed_chi`, `measured_chi`, `chi_agreement`, `resolution`, `cells`, `nonmanifold_edges`, `extractors_tested`, `c1_holds`, `c2_holds`.

**Hypothesis.** `CLAUDE.md` records that `χ` cannot be asserted on `gyroid` or `fbm_terrain`, so those fields have no topological gate beyond manifoldness. Viro's construction is a **generator**: prescribe the topology, get a field. That turns the oracle problem inside out. (C1) At least three patchworked fields with prescribed, non-trivial genus are added as fixtures, and every extractor in the crate reproduces the prescribed `χ` at sufficient resolution. (C2) The resolution at which each extractor first reproduces the prescribed `χ` is reported per extractor, which is a **sampling-adequacy curve** the crate does not currently have for any field with known topology beyond the sphere and torus. SHARE: none — this adds test coverage, and its value is the coverage.

**Falsified by.** C1 by any extractor failing to converge to the prescribed `χ` at any resolution, which would be a defect and the point of the fixture. C2 by all extractors converging at the same resolution, which would make the curve uninformative but is still worth one line. VACUITY CONTROL: at least one prescribed field must have genus above 1, or the fixture is a sphere in disguise.

#### P-137 — registered for R-137, before the harness: Deriving the case table by quantifier elimination, offline, once

**Ticket:** `R-137` (L). **Records:** `solver`, `configuration`, `wall_clock_s`, `terminated`, `cases_derived`, `cases_matching_table`, `cases_disagreeing`, `custodio_disagreement_overlap`, `machine_checkable_certificate`, `c1_holds`, `c2_holds`, `c3_holds`.

**Hypothesis.** THE REPOSITORY REJECTED CAD ON THE WRONG AXIS. `discovery-dossier.md:298` dismisses it as *"global, doubly exponential, needs arbitrary precision"* — all true, and all irrelevant, because deriving a case table is an **offline computation run once in the project's lifetime**, not a runtime cost. The A-002 series' recurring expense is transcription (`M-207`, `M-219`, `M-221`, `✗22`, `M-228`, `M-231`), which is exactly what machine derivation removes. (C1) A real-quantifier-elimination engine — nlsat via Z3, or QEPCAD — terminates on the single-cell topology query over the 8 corner parameters within a stated wall-clock budget, or is reported as not terminating, which is itself the answer. (C2) The derived taxonomy is compared against the crate's committed table, and every disagreement is adjudicated by hand against the trilinear. (C3) Any disagreement is checked against Custodio et al.'s published corrections to Chernyaev, so a machine disagreement is either a known bug or a new one. SHARE: none at runtime; this changes how the table is *produced*, not what it costs.

**Falsified by.** C1 by non-termination within budget, which closes the row honestly and is a real possibility. C2 by zero disagreements, which would be a strong independent validation of `validate_table()` and worth the whole ticket. C3 by disagreements outside Custodio's set that survive adjudication, which would be a genuine new finding. VACUITY CONTROL: the pipeline must first re-derive a **known-correct** sub-case the table already contains, or a clean run proves nothing about the solver's fidelity.

---

# Group C — deterministic topological oracles (minimal surface theory)

#### P-138 — registered for R-138, before the harness: The gyroid's Euler characteristic is `−8` per cubic cell, and the hole `CLAUDE.md` records is closable

**Ticket:** `R-138` (S). **Records:** `field`, `periods_per_axis`, `resolution`, `wrap_mode`, `chi_predicted`, `chi_measured`, `chi_agreement`, `genus_measured`, `boundary_edges`, `nonmanifold_edges`, `cells`, `c1_holds`, `c2_holds`, `c3_holds`.

**Hypothesis.** `CLAUDE.md` records that `χ` cannot be asserted on `gyroid`, which has left that field with no topological gate. It is assertable, and the reason the project could not see it is that the extraction box is not periodic-conforming. Triply periodic minimal surfaces have genus 3 **per primitive translational cell**; the gyroid's primitive lattice is **body-centred cubic**, so the conventional cubic cell holds two primitive cells and `χ = −8`, not `−4`. Verified numerically on the nodal surfaces at 32, 64, 96 and 128 voxels per period: Schwarz P `−4`, gyroid `−8`, Schwarz D `−16` per conventional cubic cell, stable across resolution. (C1) Extracting the gyroid over an integer number of periods with periodic wrap gives `χ = −8N³` exactly, for `N` in `{1, 2, 3}`, on every extractor in the crate that produces a closed surface. (C2) The measured genus is `1 + 4N³`, consistent with C1 through `χ = 2 − 2g`. (C3) The oracle is added to the validity suite as a gate for `gyroid`, closing a hole the ledger currently records as permanent. SHARE: none — this is a correctness gate where there was none.

**Falsified by.** C1 by any `N` disagreeing, which would mean the extraction is not periodic-conforming and the wrap is where the defect is — itself a finding, since chunk seams are the recurring defect class. C2 by an inconsistent genus, which would mean the surface is not closed and orientable and C1's arithmetic does not apply. C3 by the gate being unimplementable within the existing suite. VACUITY CONTROL: the non-wrapped arm must be run and must **fail** the `−8N³` prediction, or the experiment has not shown that periodicity is what matters.

#### P-139 — registered for R-139, before the harness: Schwarz P and D as second and third fields with exactly known topology

**Ticket:** `R-139` (S). **Records:** `field`, `space_group`, `primitive_lattice`, `chi_per_cubic_cell`, `chi_predicted`, `chi_measured`, `periods`, `resolution`, `body_centering_invariance`, `c1_holds`, `c2_holds`.

**Hypothesis.** The crate has exactly two fields with analytically known topology — `sphere` and `torus` — and both are genus `≤ 1`. Schwarz P (`χ = −4` per cubic cell) and Schwarz D (`χ = −16`) add two high-genus fields whose predictions differ from the gyroid's and from each other, so a single arithmetic error cannot pass all three. (C1) Both are added as reference fields and both reproduce `χ = N³·χ_cell` under periodic wrap. (C2) THE MECHANISM IS ASSERTED RATHER THAN ASSUMED: the gyroid's nodal function is invariant under the body-centring shift `(π,π,π)` while P's and D's are negated by it, which is exactly why the gyroid's cubic cell holds two primitive cells and theirs do not — checked as a symbolic identity, not inferred from the `χ` values it explains. SHARE: none.

**Falsified by.** C1 by either field disagreeing. C2 by the invariance failing, which would mean the explanation for the factor of two is wrong even if the numbers are right — and an unexplained correct number is exactly what `M-206` was before `P-123`. VACUITY CONTROL: the three fields must produce three *different* `χ` values at `N = 1`, or the suite cannot distinguish a correct oracle from a constant.

#### P-140 — registered for R-140, before the harness: Whether a periodic value-noise terrain admits an oracle at all

**Ticket:** `R-140` (M). **Records:** `field`, `octaves`, `period`, `wrap_mode`, `chi_measured`, `chi_variance_across_seeds`, `chi_stable`, `oracle_exists`, `resolution_convergence`, `c1_holds`, `c2_holds`.

**Hypothesis.** `P-138` closes `gyroid`; `fbm_terrain` remains open, and the question is whether it can be closed at all rather than whether some formula applies. A periodic lattice-noise field with a fixed seed is a **deterministic** function with a definite `χ`, so the question is empirical: does it converge with resolution, and is it stable enough to gate on? (C1) With periodic wrap and a fixed seed, measured `χ` converges as resolution increases and is stable to at least `65³`. (C2) `χ` varies across seeds, so the gate must be per-seed rather than per-field — stated in advance so nobody writes a single constant into the suite. SHARE: none.

**Falsified by.** C1 by `χ` still moving at the highest tested resolution, which would say the field is not resolved and no oracle is available at game resolutions — the likely outcome and worth recording. C2 by seed-independent `χ`, which would be surprising and would make the gate much cheaper. VACUITY CONTROL: at least three seeds must be run, or C2 cannot fire.

#### P-141 — registered for R-141, before the harness: Stratified Morse ground truth, and the citation this repository has had wrong since Phase 23

**Ticket:** `R-141` (M). **Records:** `field`, `resolution`, `chi_stratified_morse`, `chi_extracted`, `agreement`, `disagreement_cells`, `method_cost_ms`, `applicable_fields`, `doi_verified`, `c1_holds`, `c2_holds`.

**Hypothesis.** `2026-08-26-audit-and-phase-23-registrations.md:365` flags Etiene et al. 2012's DOI as unverified. It is **`10.1109/tvcg.2011.109`**, TVCG 18(6), and the paper computes ground-truth topological invariants per instance using stratified Morse theory and digital topology — a **deterministic** oracle, strictly stronger than any statistical one wherever it applies. (C1) The method computes a ground-truth `χ` on at least two fields where the crate currently has none, and agrees with `P-138`'s analytic prediction on `gyroid`, which is the cross-check that validates both. (C2) The cost is reported and the set of fields where it applies is stated, because an oracle that costs more than the mesh is a test-only instrument and should be labelled one. SHARE: none — verification cost, not runtime cost.

**Falsified by.** C1 by disagreement with `P-138`, which would put one of the two oracles in doubt and is the most informative possible outcome. C2 by inapplicability to every field the crate cares about. VACUITY CONTROL: the method must be run on `sphere` and return `χ = 2`, or it is not calibrated.

---

# Group D — element shape and anisotropy (Riemannian geometry, adaptive FEM)

#### P-142 — registered for R-142, before the harness: A metric built from the trilinear's own Hessian, and what it does to the triangle count at fixed error

**Ticket:** `R-142` (L). **Records:** `field`, `resolution`, `metric`, `p_norm`, `complexity_target`, `triangles_isotropic`, `triangles_anisotropic`, `hausdorff_isotropic`, `hausdorff_anisotropic`, `triangles_at_matched_hausdorff`, `ratio`, `aspect_ratio_max`, `aspect_ratio_mean`, `metric_ms`, `metric_share`, `c1_holds`, `c2_holds`, `c3_holds`.

**Hypothesis.** THE LARGEST PAID-FOR GAP THIS SWEEP FOUND. The corpus holds five papers on metric-based anisotropic adaptation at 0.70–0.71 — Yano & Darmofal, Cao, Bawin et al., ParMmg, Mirebeau — and the repository cites **none** of them; `Loseille`, `Alauzet`, `metric tensor` and `log-Euclidean` are zero across every file. The optimal `L^p` metric is `M_Lp = D_Lp · det(|H_u|)^(−1/(2p+d)) · |H_u|` with complexity `C(M) = ∫√det M` standing in for vertex count (NASA NTRS 20200003084, restating Loseille & Alauzet verbatim; the SIAM originals are paywalled). The Hessian is available: the crate already samples the trilinear and already computes central differences at cell size (`M-65`). (C1) At matched symmetric Hausdorff error, a metric-driven mesh uses at least 25% fewer triangles than uniform refinement on at least three of eight reference fields. (C2) The fields where it wins are the ones with a **flat direction** — ridges, creases, the CSG box edges — and the fields where it does not win are the isotropically curved ones, predicted per field **before** the run. (C3) Computing the metric costs under 15% of extraction. SHARE: C1 moves the triangle budget, whose share of frame cost `M-135` puts at 29% for the contour and 45% for the collider check.

**Falsified by.** C1 by under 25% everywhere, which would mean our fields have no exploitable anisotropy and closes Group D's largest row. C2 by the win landing on the wrong fields, which would mean the mechanism is not the one claimed even if the number is real. C3 by above 15%. VACUITY CONTROL: the reported maximum aspect ratio must exceed 3 on at least one field, or the "anisotropic" arm is isotropic and C1 is measuring noise.

#### P-143 — registered for R-143, before the harness: The rate does not improve and the constant does, stated before the measurement rather than discovered in it

**Ticket:** `R-143` (M). **Records:** `field`, `resolution_series`, `fitted_exponent_isotropic`, `fitted_exponent_anisotropic`, `exponent_difference`, `fitted_constant_isotropic`, `fitted_constant_anisotropic`, `constant_ratio`, `am_gm_gap`, `flat_direction_fraction`, `c1_holds`, `c2_holds`, `c3_holds`.

**Hypothesis.** THE SEDUCTIVE VERSION OF GROUP D IS WRONG AND THE REGISTRATION SAYS SO IN ADVANCE. For `W^{2,p}` regularity, uniform refinement gives `O(N^(−2/3))` in 3D and the optimal anisotropic mesh gives `O(N^(−2/3))` — **the same exponent**; Bonito, Canuto, Nochetto & Veeser (Acta Numerica 2024, in corpus at `10.1017/s0962492924000011`) state flatly that the order *"cannot be improved upon assuming either higher regularity … or a graded mesh"*. What anisotropy buys is the constant: `‖√|det H|‖_{L^τ}` replaces `|f|_{W^{2,p}}`, and by AM–GM the former is never larger. (C1) The fitted error exponent is statistically indistinguishable between isotropic and anisotropic arms — `|Δexponent| < 0.1` — on every smooth field. (C2) The fitted **constant** improves, and its improvement tracks the measured AM–GM gap `‖√|det H|‖ / ‖tr|H|/d‖` per field, rank correlation above 0.7. (C3) A field lacking `W^{2,p}` regularity — `csg_difference`, whose gradient is discontinuous — shows an exponent difference above 0.1, because that is the one regime where grading wins in the exponent. SHARE: none — this is `M-12`'s law, decomposed.

**Falsified by.** C1 by an exponent gain on a smooth field, which would contradict the cited theory and would be the more interesting result. C2 by a constant improvement uncorrelated with the AM–GM gap, which would mean the mechanism is not the one the theory names. C3 by no exponent difference on the CSG field, which would say our sharp fields are smoother than assumed. VACUITY CONTROL: the AM–GM gap column must vary by at least `2×` across fields, or C2's correlation is against a constant.

#### P-144 — registered for R-144, before the harness: Metric interpolation across a chunk seam, where averaging matrices is not averaging

**Ticket:** `R-144` (M). **Records:** `interpolation_scheme`, `seam_axis`, `determinant_swell_max`, `determinant_swell_mean`, `seam_vertices_moved`, `seam_open_edges`, `hashes_moved`, `bit_exact_seam`, `cell_size_power_of_two`, `c1_holds`, `c2_holds`, `c3_holds`.

**Hypothesis.** A metric field must be interpolated, and component-wise interpolation of symmetric positive-definite matrices is not intrinsic — it can swell determinants, which in mesh terms means spuriously coarsening. The log-Euclidean scheme `exp((1−t)log M₁ + t log M₂)` is determinant-monotone. This lands directly on `M-32`: seams are bit-exact only when the cell size is a power of two, and a matrix logarithm is not obviously going to preserve that. (C1) Component-wise interpolation produces measurable determinant swell — above 5% on at least one field — and log-Euclidean does not. (C2) Log-Euclidean interpolation preserves `M-32`'s power-of-two bit-exactness, or it does not, and which is stated as the result rather than assumed. (C3) The seam remains closed under both schemes: zero additional open edges against the existing seam counter. SHARE: C1 moves the metric-construction stage, whose share `P-142` C3 reports.

**Falsified by.** C1 by no swell, which would mean the intrinsic scheme is unnecessary here and saves the transcendental cost. C2 by loss of bit-exactness, which would make log-Euclidean incompatible with the crate's determinism goal and is the outcome that decides the row. C3 by new open edges under either scheme. VACUITY CONTROL: the fixture must contain a seam where the two chunks' metrics genuinely differ, reported as a metric-distance column, or both schemes are interpolating between equal endpoints.

#### P-145 — registered for R-145, before the harness: The features where a flat direction exists, isolated and measured on their own

**Ticket:** `R-145` (M). **Records:** `feature_class`, `principal_curvature_ratio`, `cells`, `triangles_isotropic`, `triangles_anisotropic`, `saving`, `hausdorff_delta`, `saving_vs_curvature_ratio_correlation`, `c1_holds`, `c2_holds`.

**Hypothesis.** `P-142` measures the aggregate; this measures the mechanism. The AM–GM gap collapses toward zero exactly where one principal curvature vanishes, so the saving should be a function of the curvature ratio and nothing else. (C1) Binning surface cells by principal-curvature ratio, the anisotropic triangle saving is monotone in the ratio, with correlation above 0.7. (C2) Cells with ratio near 1 — spherical caps — show **no** saving, which is the control that says the effect is anisotropy rather than a general refinement improvement. SHARE: C1 decomposes `P-142`'s C1 saving by feature class; the bin populations are reported so each saving has a denominator.

**Falsified by.** C1 by a non-monotone relationship, which would mean the saving comes from somewhere other than anisotropy. C2 by a saving on isotropic cells, same conclusion and more damning. VACUITY CONTROL: every curvature-ratio bin must be non-empty on at least one field, reported as counts, or the monotonicity is fitted through gaps.

#### P-146 — registered for R-146, before the harness: A null registered on purpose — the `L¹` anomaly in the source's own table

**Ticket:** `R-146` (S). **Records:** `norm`, `metric`, `element_count`, `error_measured`, `optimal_metric_wins`, `anomaly_reproduced`, `cao_table_agreement`, `c1_holds`, `c2_holds`.

**Hypothesis.** THIS ROW EXISTS BECAUSE THE HEADLINE NUMBER FOR GROUP D IS NOT WHAT IT LOOKS LIKE. Cao's Table 2 (`10.1090/S0025-5718-07-01981-3`, in corpus) does show `5.43e-7` versus `2.79e-6` at `N_e ≈ 16,000` — a real 5.1× spread, and the element counts match to within 1.4%, not exactly. But the two metrics are optimised for **different norms**: in `H¹` the ranking reverses. And Cao writes that his theory holds *"in all the cases except the `L¹`-error of quadratic interpolation"* — the `L¹` column is the one place his own optimal metric does not win. Building on that column is the weakest available ground and the registration says so before anyone quotes 5.1× as a free win. (C1) On our fields, the metric optimised for a given norm wins in that norm — the property Cao's theory asserts and his `L¹` column violates. (C2) Any violation we find is characterised, and the comparison is always reported **within** a norm, never across two. SHARE: none — this is a methodology guard.

**Falsified by.** C1 by a violation, which would reproduce Cao's anomaly on different data and would be a genuine contribution to a question his paper leaves open. C2 by an uncharacterisable violation. VACUITY CONTROL: at least three norms must be run, or "within a norm" is not a constraint.

#### P-147 — registered for R-147, before the harness: Refining where the player is looking, with forty years of theory behind it

**Ticket:** `R-147` (M). **Records:** `criterion`, `functional`, `triangles`, `screen_space_error`, `screen_error_at_matched_triangles`, `off_screen_triangles`, `adjoint_ms`, `adjoint_share`, `pop_magnitude_cells`, `c1_holds`, `c2_holds`, `c3_holds`.

**Hypothesis.** v1's Axis 7 named goal-oriented adaptivity and nothing has used it. The dual-weighted-residual method refines to reduce error **in a chosen output functional** rather than globally, and for a game the functional is screen-space error under the current camera. Yano & Darmofal (in corpus, `10.2514/6.2012-79`) give the optimisation framing: *"both the sizing decision and the anisotropy decision are driven directly by the behavior of the a posteriori error estimate."* (C1) At matched triangle count, a goal-oriented criterion beats camera-distance LOD on screen-space error by at least 20%. (C2) The saving comes from **not refining what is off-screen or back-facing**, reported as an off-screen triangle count, rather than from a general improvement. (C3) The adjoint cost is under 10% of extraction, or the method is a non-starter for a frame budget and should be said so. SHARE: C1 moves the LOD selection stage; `M-121` measured a level change moving the surface by up to 3.14 cells, which is the pop this criterion is meant to spend better.

**Falsified by.** C1 by under 20%. C2 by a saving that persists with everything on-screen, which would mean the mechanism is not visibility. C3 by above 10%. VACUITY CONTROL: the camera must be positioned so that a non-trivial fraction of the scene is off-screen, reported as a fraction, or C2 cannot fire.

---

# Group E — multiscale flatness (harmonic analysis, quantitative rectifiability)

#### P-148 — registered for R-148, before the harness: The `β`-number is the planarity assumption the QEF never checks

**Ticket:** `R-148` (M). **Records:** `field`, `resolution`, `beta_infinity`, `beta_times_diam`, `qef_residual`, `rank_correlation`, `cells`, `beta_ms`, `beta_share`, `singular_cells`, `c1_holds`, `c2_holds`, `c3_holds`.

**Hypothesis.** Dual Contouring minimises squared distance to tangent planes, which silently assumes the cell's surface patch is near-planar. Jones' `β_∞(Q)` — the width of the thinnest slab containing the patch, over `diam Q` — is the formal name for that assumption and comes with theorems; Azzam & Schul's higher-dimensional Traveling Salesman Theorem (`arXiv:1609.02892`, acquired and indexed this session, 79 pp.) is the surface case. `β` is absent from the corpus and from every file in the repo. (C1) `β_∞(Q)·diam(Q)` is computable per cell at under 10% of extraction cost. (C2) It rank-correlates with the QEF residual above 0.7 on at least six of eight fields. (C3) IT CARRIES SIGNAL THE RESIDUAL DOES NOT: there exist cells with low residual and high `β`, and the count is reported — those are the cells where the plane fit is confidently wrong. SHARE: C1 moves the vertex-placement stage, whose cost `M-25` puts at 3% over Surface Nets.

**Falsified by.** C1 by above 10%. C2 by correlation below 0.7 on three or more fields, which would mean `β` and the residual measure different things and C3 becomes the whole row. C3 by an empty set, which would mean `β` is a reparametrisation of the residual and buys a theorem but no new signal — the most likely outcome and the reason C3 exists. VACUITY CONTROL: the residual column must have non-zero variance on every field, or the correlation is against a constant.

#### P-149 — registered for R-149, before the harness: `β` against curvature against camera distance, as a refinement criterion

**Ticket:** `R-149` (M). **Records:** `criterion`, `field`, `triangles`, `hausdorff`, `hausdorff_at_matched_triangles`, `criterion_ms`, `rank_correlation_with_error`, `worst_decile_hausdorff`, `c1_holds`, `c2_holds`.

**Hypothesis.** Three candidate refinement criteria now exist — `β`, estimated curvature, and camera distance — and only the third has ever been used. (C1) At matched triangle count, `β`-driven refinement beats camera-distance refinement on symmetric Hausdorff error on at least five of eight fields. (C2) `β` beats curvature specifically on fields with **sub-cell structure**, where curvature estimated at the cell size is aliased and `β` is not, because `β` is defined by a slab fit rather than by a derivative — `thin_plate` and `noise_cavity` are the named predictions. SHARE: C1 moves the adaptivity stage, and the triangle counts are reported so the comparison is at a fixed budget.

**Falsified by.** C1 by under five fields. C2 by curvature matching `β` on `thin_plate`, which would mean the aliasing argument is wrong and `M-72`'s mechanism does not reach the criterion. VACUITY CONTROL: the worst-decile Hausdorff column must differ across criteria, or all three are refining the same cells.

#### P-150 — registered for R-150, before the harness: A triangle budget computed before meshing, from a convergent sum

**Ticket:** `R-150` (M). **Records:** `field`, `resolution`, `beta_sum`, `predicted_triangles`, `actual_triangles`, `prediction_error`, `prediction_error_ratio`, `sum_converges`, `scales_used`, `c1_holds`, `c2_holds`.

**Hypothesis.** The Traveling Salesman Theorem's content is that `∑_Q β(Q)²·diam(Q)^d` converges precisely for rectifiable sets — which turns multiscale flatness into a **number computed over the whole field before any mesh exists**. That is an a-priori budget, which no criterion in the crate currently provides; `M-13` gives `surface cells ≈ 1.5·A/h²` after the fact. (C1) The `β`-sum predicts final triangle count within 25% on at least five of eight fields at three resolutions. (C2) The sum converges as scales are added — reported as a convergence column — because a divergent sum would mean the field is not rectifiable at the resolutions we mesh at, which would be a more interesting finding than the budget. SHARE: none — this predicts a cost, it does not change one.

**Falsified by.** C1 by above 25% on four or more fields, which would leave `M-13`'s area law the better predictor and close the row. C2 by divergence, which reframes the whole thing. VACUITY CONTROL: at least four dyadic scales must be summed, or convergence cannot be observed.

---

# Group F — reconstruction order (approximation theory)

#### P-151 — registered for R-151, before the harness: Predicting M-12's fitted constant instead of fitting it

**Ticket:** `R-151` (S). **Records:** `field`, `resolution_series`, `fitted_constant`, `fitted_constant_ci`, `predicted_constant`, `ratio`, `within_ci`, `strang_fix_order`, `measured_order`, `c1_holds`, `c2_holds`.

**Hypothesis.** `M-12` measured Marching Cubes' error falling like `h²` with a **fitted** constant, and `M-113` found the fitted constant does not survive across configurations. Strang–Fix theory says the order is a property of the filter and gives a closed form for the asymptotic constant; the trilinear's order-2 is then a derivation rather than a measurement. (C1) The Strang–Fix approximation order of the tensor-product hat is 2, derived and asserted as a test, matching `M-12`'s measured exponent. (C2) The predicted asymptotic constant reproduces `M-12`'s fitted constant within its confidence interval on at least four of eight fields. SHARE: none — this converts a fit into a prediction.

**Falsified by.** C1 by a derived order other than 2, which would mean the filter is not what we think it is. C2 by the prediction missing on five or more fields, which would mean the asymptotic regime is not reached at game resolutions — a genuinely useful negative, because it would say `M-12`'s law is empirical rather than asymptotic and should not be extrapolated. VACUITY CONTROL: the fitted constants must differ across fields by at least `2×`, or the prediction is matching a universal constant by accident.

#### P-152 — registered for R-152, before the harness: A compact prefilter instead of a truncated recursive one

**Ticket:** `R-152` (M). **Records:** `prefilter`, `support`, `field`, `root_position_error`, `hausdorff`, `prefilter_ms`, `chunk_local`, `seam_bit_exact`, `hashes_moved`, `vs_truncated_recursive`, `c1_holds`, `c2_holds`, `c3_holds`.

**Hypothesis.** `✗42`/`M-359` measured shifted linear interpolation and found the recursive prefilter's non-locality **bounded** — truncation at `k ≥ 10` moves the root under `7.15e-7` cells on all eight fields — so the chunk-locality objection is already retired empirically. Quasi-interpolation achieves full approximation order with a **compactly supported** prefilter, making locality exact rather than bounded (Blu & Unser, `10.1006/acha.1998.0249`; the corpus holds only the landing page, so the construction must be taken from a freely available restatement and that source named in the finding). (C1) A compactly supported quasi-interpolant matches the truncated recursive prefilter's accuracy to within 5%. (C2) It is exactly chunk-local — zero dependence outside the chunk, asserted structurally rather than measured. (C3) Chunk seams remain bit-exact where `M-32` says they should be, and the golden hashes that move are counted, since any reconstruction change is a re-baseline. SHARE: C1 moves the vertex-interpolation stage, the crate's single most-executed operation.

**Falsified by.** C1 by above 5%, which would mean compactness costs accuracy and the truncated recursive filter is the better engineering answer. C2 by any external dependence, which would mean the quasi-interpolant was mis-constructed. C3 by loss of seam bit-exactness. VACUITY CONTROL: the truncated-recursive arm must reproduce `✗42`'s `7.15e-7` figure as a control, or the comparison has no baseline.

#### P-153 — registered for R-153, before the harness: Raising the order, which is the lever `✗42` did not pull

**Ticket:** `R-153` (M). **Records:** `filter`, `approximation_order`, `field`, `resolution_series`, `fitted_exponent`, `hausdorff`, `samples_per_cell`, `eval_ms`, `eval_share`, `cases_invalidated`, `hashes_moved`, `c1_holds`, `c2_holds`, `c3_holds`.

**Hypothesis.** `✗42` tested a **knot shift** at fixed order and found the gain maps to root position only as a lottery over the crossing position. The untested lever is the **order** itself: a filter with approximation order 4 changes the exponent, not the constant, and that is a different mechanism with a different failure mode. (C1) A higher-order reconstruction gives a fitted error exponent above 3 on at least four smooth fields, against the trilinear's 2. (C2) The cost is stated honestly in the currency that matters: field evaluations per cell, and the resulting share of extraction time, because a wider stencil is more samples and samples are this crate's dominant cost on procedural fields. (C3) THE COLLATERAL IS PRICED AT REGISTRATION: the entire A-002 apparatus — the asymptotic decider, the interior test, the body-saddle algebra, and `P-123`'s identity — is derived from the **trilinear** and would need re-deriving; the count of invalidated cases is a reported column, not a discovery. SHARE: C2 moves the field-evaluation stage, which `M-152` puts at 2.65 ms of an 8.40 ms upload at `129³`.

**Falsified by.** C1 by an exponent at or below 3, which would mean the order does not materialise on our fields at our resolutions. C2 by a share that makes it unaffordable regardless of accuracy. C3 by a smaller invalidation than predicted, which would be good news and should be recorded as such. VACUITY CONTROL: the trilinear arm must reproduce exponent 2 in the same harness, or the comparison is between two different measurement setups.

#### P-154 — registered for R-154, before the harness: A null registered on purpose — sparse grids need a regularity our sharp fields have not got

**Ticket:** `R-154` (M). **Records:** `field`, `mixed_derivative_norm`, `norm_finite`, `sparse_grid_points`, `full_grid_points`, `point_ratio`, `hausdorff_sparse`, `hausdorff_full`, `fields_qualifying`, `c1_holds`, `c2_holds`.

**Hypothesis.** Sparse grids give `O(h⁻¹·|log h|^(d−1))` points instead of `O(h⁻ᵈ)` — but **only for functions with bounded mixed second derivatives** (`H²_mix`; Bungartz & Griebel, `10.1017/s0962492904000182`). This registration exists to test the hypothesis rather than the method, because the hypothesis is where it dies. (C1) The mixed-derivative norm `‖∂⁶u/∂x²∂y²∂z²‖` is measured on all eight reference fields, and is unbounded on the CSG fields — `box_exact` and `csg_difference` — where the surface has sharp edges. (C2) On the fields where it *is* bounded, a sparse grid reaches matched Hausdorff error with at least `2×` fewer samples. SHARE: C2 moves the field-evaluation stage on qualifying fields only, and the qualifying set is C1's output.

**Falsified by.** C1 by a bounded norm on the CSG fields, which would be surprising and would open the method to everything. C2 by under `2×` on the smooth fields, which would say the asymptotic advantage does not materialise at game resolutions — the expected outcome, and the reason this row is `M` rather than `L`. VACUITY CONTROL: `sphere` must show a bounded norm, or the measurement is broken rather than the fields being rough.

---

# Group G — the denominator (information-based complexity, nonlinear approximation)

#### P-155 — registered for R-155, before the harness: How far from optimal the extractor is, and the null that would be worth the most

**Ticket:** `R-155` (M). **Records:** `field_class`, `regularity`, `n_samples`, `minimal_error_rate`, `measured_error_rate`, `ratio`, `order_optimal`, `constant_gap`, `class_assumptions`, `c1_holds`, `c2_holds`, `c3_holds`.

**Hypothesis.** EVERY NUMBER IN THIS LEDGER IS A NUMERATOR. `M-12`'s constant, `F-007`'s 13–15%, `✗42`'s ratio — each says *did this help* and none says *how much is left*. Information-based complexity answers the second, and field evaluation is this crate's dominant cost, so the bound is denominated in the currency actually spent (Krieg & Ullrich, `arXiv:2602.02066`, acquired and indexed this session). (C1) The `n`-th minimal error rate for a stated regularity class is established from the literature and written down with **all** its hypotheses, and the class our fields plausibly belong to is named. (C2) The measured extractor error is compared against it, and the answer is stated as a ratio. (C3) IF THE MINIMAL RATE IS ALREADY `Θ(h²)`, MARCHING CUBES IS ORDER-OPTIMAL AND ONLY THE CONSTANT IS IN PLAY — which would cap `P-153` and much of Group D's ambition at "a better constant", and is the most valuable single sentence this phase could produce. SHARE: none — this is a denominator, not a saving.

**Falsified by.** C1 by no rate being available for any class our fields belong to, which would say IBC does not reach this problem and should be closed out. C2 by a ratio that is not computable because the class assumptions are unverifiable on procedural fields. C3 is not falsifiable — it is a branch, and both branches are recorded. VACUITY CONTROL: the class membership must be **argued from a measured property of each field** — a regularity estimate, reported as a column — not assumed, or the whole comparison is against a class we merely hope applies.

#### P-156 — registered for R-156, before the harness: Whether adaptivity beats uniform sampling, measured because the famous theorem does not apply

**Ticket:** `R-156` (M). **Records:** `field`, `sample_budget`, `hausdorff_uniform`, `hausdorff_adaptive`, `gain`, `gain_exceeds_two`, `class_convex`, `class_symmetric`, `operator_linear`, `c1_holds`, `c2_holds`, `c3_holds`.

**Hypothesis.** IT IS TEMPTING TO CITE A BOUND HERE AND IT WOULD BE WRONG. Gal & Micchelli and Novak (`10.1006/jcom.1996.0015`) prove adaption improves worst-case error by at most a factor of two — under four hypotheses: the class is **convex**, the class is **symmetric**, the error is **worst-case**, and the solution operator is a **continuous linear mapping**. Level-set extraction is not a continuous linear operator on a convex symmetric class, and the same literature records gains **up to order `n`** for nonlinear problems with restricted information. So the cap does not bind and the question is empirical. (C1) At matched sample budget, octree-adaptive sampling beats uniform on symmetric Hausdorff error, and the gain is reported per field. (C2) The gain **exceeds 2** on at least one field, which is the concrete demonstration that the factor-2 result does not apply here. (C3) The four hypotheses are each checked against our setting and recorded as a row of booleans, so no future reader re-imports the bound. SHARE: C1 moves the field-evaluation stage at fixed budget.

**Falsified by.** C1 by adaptive losing, which would be a strong argument against the whole octree line. C2 by every gain falling under 2, which would be consistent with the theorem applying after all and would send C3 back for a harder look. C3 by any hypothesis being arguably satisfied. VACUITY CONTROL: the adaptive arm must actually vary its sample density, reported as a density ratio, or both arms are uniform.

#### P-157 — registered for R-157, before the harness: Which approximation class each field is in, which decides whether LOD can help it at all

**Ticket:** `R-157` (M). **Records:** `field`, `s_uniform`, `s_adaptive`, `s_difference`, `class_membership`, `n_term_decay`, `lod_asymptotic_gain`, `c1_holds`, `c2_holds`.

**Hypothesis.** The corpus holds four adaptive-approximation papers at 0.62–0.64 — Karkulik & Praetorius, Gaspoz & Morin, Mommer & Stevenson, and the Bonito/Canuto/Nochetto/Veeser survey — and `Stevenson` is cited in this repository **only** for the newest-vertex-bisection closure bound, never for the approximation-class theory in the same PDFs. The class `A^s` is the right object: if a field's `s` is the same under uniform and adaptive refinement, octree LOD provably buys nothing asymptotically on that field. (C1) `s` is estimated for both refinement strategies on all eight fields. (C2) The fields split: at least one has `s_adaptive > s_uniform` and at least one does not, so the answer is field-dependent and the crate can say which fields LOD is for. SHARE: none — this predicts where a stage helps, it does not change the stage.

**Falsified by.** C1 by `s` not being estimable from the achievable resolution range, which is a real risk and should close the row honestly. C2 by every field behaving alike, which would be a simpler and still useful answer. VACUITY CONTROL: at least four resolutions per field, with the fitted `s` carrying a confidence interval, or the split in C2 is fitted noise.

---

# Group H — the sampling lattice (lattice theory, information theory)

#### P-158 — registered for R-158, before the harness: 0.257 dB is a prediction, not a hope

**Ticket:** `R-158` (M). **Records:** `lattice`, `G`, `samples`, `field`, `hausdorff`, `hausdorff_ratio`, `predicted_gain_db`, `measured_gain_db`, `prediction_holds`, `extraction_ms`, `case_table_size`, `c1_holds`, `c2_holds`, `c3_holds`.

**Hypothesis.** `G(Z³) = 1/12 = 0.0833333` and `G(A₃*) = 19/(192·∛2) = 0.078543281`, both verified by direct Voronoi-cell integration this session, giving a **5.748% MSE reduction, 0.2571 dB** at matched point density. Barnes & Sloane (`10.1137/0604005`) prove `A₃*` optimal **among 3D lattices** — not among all quantizers, where the bracketing bounds remain unproven. (C1) At matched sample count, BCC sampling improves symmetric Hausdorff error over the cubic grid on at least five of eight fields. (C2) The improvement is in the neighbourhood of the predicted 0.257 dB — within a factor of 2 — or the deviation is explained, because `G` predicts *quantization* error and Hausdorff error is not that, and the gap between them is the interesting part. (C3) The cost is stated: a BCC extractor needs a different case table, and its size is reported alongside the cubic 256. SHARE: C1 moves the whole sampling stage; C3 is a complexity cost, not a runtime one.

**Falsified by.** C1 by under five fields. C2 by an unexplained deviation beyond `2×`, which would mean `G` is the wrong figure of merit for surface extraction and would be worth writing down. C3 by a case table large enough to be impractical. VACUITY CONTROL: both arms must be at genuinely matched point density, reported as counts, or the comparison is a resolution change wearing a lattice's name.

#### P-159 — registered for R-159, before the harness: A null registered on purpose — FCC versus BCC is 0.011 dB and should be unmeasurable

**Ticket:** `R-159` (S). **Records:** `lattice`, `G`, `predicted_gap_db`, `measured_gap_db`, `measurement_scatter_db`, `gap_below_scatter`, `c1_holds`.

**Hypothesis.** `G(D₃) = 2^(−11/3) = 0.078745066` against `G(A₃*) = 0.078543281` — a gap of **0.011 dB**, about 4% of the size of the cubic-vs-BCC gap. This registration exists so that a null is on the record before somebody spends a week distinguishing them. (C1) The measured FCC–BCC Hausdorff difference is **below the measurement scatter** of the harness itself, reported as a scatter column from repeated runs. SHARE: none.

**Falsified by.** C1 by a resolvable difference, which would mean either the harness is more precise than expected — useful — or that Hausdorff error responds to something `G` does not capture, which is `P-158` C2's question answered from the other side. VACUITY CONTROL: the harness scatter must be estimated from at least five repeated runs of the *same* lattice, or "below scatter" has no denominator.

#### P-160 — registered for R-160, before the harness: The BCC box-spline paper this project has owned and never opened

**Ticket:** `R-160` (M). **Records:** `filter`, `lattice`, `approximation_order`, `support_size`, `evaluations_per_cell`, `hausdorff`, `eval_ms`, `vs_trilinear_on_cubic`, `c1_holds`, `c2_holds`.

**Hypothesis.** Entezari, Van De Ville & Möller, *Practical box splines on the BCC lattice* (`10.1109/tvcg.2007.70429`), is in the corpus and the 2026-08-23 sweep recorded it as *"✅ (today) | 0 files"* — acquired and never cited, now for a second time. It is the reconstruction half of `P-158`'s sampling half, and the two are only meaningful together: a better lattice sampled with the wrong filter is not obviously better at anything. (C1) A BCC box spline achieves at least the trilinear's approximation order on the BCC lattice, derived and then measured. (C2) The combination — BCC lattice plus box spline — beats cubic-plus-trilinear at matched sample count by more than the lattice change alone does in `P-158`, which is the test of whether the two compose. SHARE: C2 moves the sampling and reconstruction stages jointly, and `P-158`'s C1 is the baseline it is measured against.

**Falsified by.** C1 by a lower order, which would mean the lattice gain is paid back at reconstruction. C2 by no improvement over `P-158` alone, which would mean they do not compose and the filter is doing nothing. VACUITY CONTROL: `P-158` must have completed and its per-field numbers must be the reported baseline, or C2 has nothing to exceed.

---

# Group I — optimality of greedy (combinatorial optimisation, graph theory)

#### P-161 — registered for R-161, before the harness: How far greedy meshing is from the optimum, which is computable in polynomial time

**Ticket:** `R-161` (M). **Records:** `field`, `resolution`, `vertices_n`, `holes_h`, `good_diagonals_g`, `optimum_rectangles`, `greedy_rectangles`, `ratio`, `matching_ms`, `optimum_ms`, `c1_holds`, `c2_holds`, `c3_holds`.

**Hypothesis.** `M-56` measured greedy meshing's saving over face culling at 1.70× to 256× and called it field-dependent — but never measured it against the **optimum**, which for a rectilinear region is a classical polynomial-time problem. Eppstein (`arXiv:0908.3916`, acquired and converted this session) gives it exactly: the minimum number of rectangles is **`n/2 + h − g − 1`**, where `g` is the maximum number of disjoint good diagonals — axis-parallel interior segments joining two concave vertices — found as a maximum independent set in a **bipartite** intersection graph of horizontal against vertical chords, hence by König's theorem via maximum matching, in **`O(n^{3/2} log n)`** with range searching. Independently discovered by Lipski et al. (1979), Ohtsuki (1982), and Ferrari, Sankar & Sklansky (1984). (C1) The optimum is computed on all seven fields at three resolutions. (C2) Greedy's ratio to the optimum is reported per field, turning `M-56`'s unbounded range into a distance. (C3) The optimum's cost is under 10× greedy's, or it is a measurement instrument rather than a shippable algorithm and is labelled one. SHARE: C2 moves the greedy-meshing stage, whose saving `M-56` bounds.

**Falsified by.** C1 by the reduction not applying — for instance if the merged regions have holes the formula's `h` does not model. C2 by greedy already achieving the optimum on every field, which would close `M-56`'s question with the best possible answer. C3 by above 10×. VACUITY CONTROL: at least one field must have `g > 0` and at least one `h > 0`, reported as columns, or the formula is being checked in its degenerate case.

#### P-162 — registered for R-162, before the harness: A null registered on purpose — greedy meshing is not a matroid problem, but LOD budgeting might be

**Ticket:** `R-162` (M). **Records:** `problem`, `is_matroid`, `is_submodular`, `marginal_returns_monotone`, `greedy_ratio`, `bound_applies`, `violating_chunks`, `c1_holds`, `c2_holds`.

**Hypothesis.** Greedy is provably optimal on matroids and within `(1−1/e)` for monotone submodular maximisation under a cardinality constraint (Nemhauser, Wolsey & Fisher; paywalled, so the statement must be sourced from a freely available restatement). Greedy **meshing** is a minimum-cardinality cover, not a max-weight independent set, so no matroid — that half is registered as a negative. The LOD **chunk budget** is the other half: choose which chunks to refine under a fixed frame budget to minimise screen-space error. (C1) Greedy meshing is shown not to be a matroid problem, by exhibiting a counterexample to the exchange property. (C2) The LOD budget's marginal screen-space-error reduction is measured per chunk and shown to be diminishing — or a chunk with **increasing** marginal returns is exhibited, which would kill the `(1−1/e)` guarantee outright and is the falsifier that matters. SHARE: C2 governs the LOD selection stage, whose amortised cost `M-124` tracks to within one chunk.

**Falsified by.** C1 by the exchange property holding, which would be a much stronger result. C2 by increasing marginal returns, which closes the guarantee — and finding one such chunk is cheaper than proving none exists, so the experiment is designed to look for it. VACUITY CONTROL: the marginal-return curve must be computed over at least 50 chunks spanning at least two LOD levels, or "diminishing" is a claim about a handful of points.

---

# Group J — the case table as a Boolean function (Fourier analysis on the hypercube)

#### P-163 — registered for R-163, before the harness: The Fourier spectrum of the 256-case table, and whether it is low-degree

**Ticket:** `R-163` (M). **Records:** `output_bit`, `fourier_weight_by_degree`, `spectral_concentration`, `anf_terms`, `max_degree`, `sparse`, `branchless_feasible`, `eval_ms_table`, `eval_ms_spectral`, `c1_holds`, `c2_holds`, `c3_holds`.

**Hypothesis.** Marching Cubes is a lookup on a Boolean function of 8 variables, and the analysis of Boolean functions has exact machinery for such objects — Walsh–Hadamard expansion, degree, spectral concentration, algebraic normal form (O'Donnell, `arXiv:2105.10386`, acquired and indexed this session). Absent from the corpus at 0.584 and from every file in the repo. (C1) The Walsh–Hadamard spectrum of the table, treated as a vector-valued Boolean function, is computed and its weight distribution by degree reported. (C2) The spectrum is **not** concentrated on low degrees — the honest prediction, since a case table encoding topology should need high-degree terms — and the ANF term count is reported so the negative is quantitative. (C3) If C2 is falsified and the spectrum *is* sparse, a branchless spectral evaluation is benchmarked against the table lookup. SHARE: C3 moves the classification stage, which the work-graphs result put at a 2.8–3.4× loss when moved to irregular GPU dispatch.

**Falsified by.** C1 by the transform not being well-defined for the vector-valued output, which would mean the framing needs per-bit treatment and should say so. C2 by concentration on low degrees, which is the good outcome and triggers C3. C3 by the spectral evaluation losing to a lookup, which is likely — a 256-entry table is already in L1. VACUITY CONTROL: a known low-degree function (parity, majority) must be run through the same transform and reproduce its known spectrum, or the instrument is unvalidated.

#### P-164 — registered for R-164, before the harness: A null registered on purpose — noise stability as a robustness number, and whether it predicts anything measured

**Ticket:** `R-164` (M). **Records:** `noise_rate`, `noise_stability`, `predicted_flip_rate`, `measured_flip_rate_f32`, `agreement`, `ulp_perturbation`, `topology_changes`, `field`, `c1_holds`, `c2_holds`.

**Hypothesis.** Noise stability — the probability that a Boolean function's output survives independently flipping each input with probability `ρ` — is a principled measure of how robust the extracted topology is to perturbing sample values, and it is computable directly from `P-163`'s spectrum. The question is whether it describes the perturbations this crate actually suffers. (C1) Noise stability is computed for the case table and its predicted case-flip rate at a given `ρ` is stated. (C2) THE PREDICTION IS TESTED AGAINST THE WRONG NOISE MODEL ON PURPOSE, AND SAID SO: `f32` rounding is not independent per-corner bit flips — it is a deterministic, spatially correlated perturbation — so the predicted and measured flip rates are expected to **disagree**, and the size of the disagreement measures how far the independent-noise model is from this crate's reality. SHARE: none — this is a modelling check.

**Falsified by.** C1 by stability not being computable from the spectrum, which would send `P-163` back. C2 by agreement, which would mean the independent-flip model describes `f32` behaviour after all — surprising, and it would make noise stability a genuinely usable robustness number rather than a mismatched one. VACUITY CONTROL: the measured flip rate must be non-zero at some perturbation magnitude, reported as a sweep, or C2 compares two zeros.

#### P-165 — registered for R-165, before the harness: Which corner matters, as an influence rather than an intuition

**Ticket:** `R-165` (S). **Records:** `corner_index`, `influence`, `total_influence`, `symmetry_class`, `influence_equal_within_class`, `refinement_priority_correlation`, `c1_holds`, `c2_holds`.

**Hypothesis.** The **influence** of a variable is the probability that flipping it changes the output — a per-corner number the crate could use to decide which sample to refine first, and a sanity check on the table's symmetry. (C1) All eight corner influences are computed and are **equal within each octahedral symmetry class**, which they must be if the table is octahedrally correct — making this a cheap independent check on `validate_table()`. (C2) Total influence is reported and compared against the average sensitivity a refinement heuristic would need. SHARE: none.

**Falsified by.** C1 by unequal influences within a symmetry class, which would be a table defect and worth the whole ticket by itself. C2 by total influence being uninformative, which is the expected outcome for a symmetric table. VACUITY CONTROL: a deliberately corrupted table must produce unequal within-class influences, or C1 cannot detect the defect it is designed for.

---

# Group K — foundations and their hypotheses (numerical continuation, differential topology, dynamical systems)

#### P-166 — registered for R-166, before the harness: Marching Cubes has a 1985 ancestor, and it proves less than it is tempting to claim

**Ticket:** `R-166` (S). **Records:** `source`, `year`, `doi`, `doi_verified`, `claim_type`, `hypothesis`, `guarantee`, `is_homeomorphism`, `descendant_citations_in_repo`, `root_citations_in_repo`, `c1_holds`, `c2_holds`, `c3_holds`.

**Hypothesis.** This repository cites Plantinga–Vegter in 18 files and `Freudenthal`/`Kuhn` in 6 including `marching_tetrahedra/table.rs`, and cites the parent theory in **none**. (C1) Allgower & Schmidt, *An algorithm for piecewise linear approximation of an implicitly defined manifold*, SIAM J. Numer. Anal. **22(2):322–346, 1985**, DOI **`10.1137/0722020`**, predates Lorensen & Cline (1987) and describes PL approximation of `H⁻¹(0)` by simplicial subdivision — recorded with the DOI verified, because **`10.1137/0722019` is a different paper** (Nanda on the QR algorithm, same volume and issue) and is the trap. (C2) THE LIMIT IS RECORDED WITH THE PRIORITY: AS85 proves `‖H(x)‖_∞ < ε` along the constructed manifold under a **full-rank** hypothesis at the seed — a residual bound with a regularity condition, **not** a homeomorphism; the isotopy result for isomanifolds is Boissonnat & Wintraecken, SoCG 2020, whose own text says the earlier result holds *"without a homeomorphism with the zero set of f"*. (C3) Dobkin, Levy, Thurston & Wilks (ACM TOG 9(4), 1990, `10.1145/88560.88575`) brought PL continuation into graphics — and traces **curves**, not surfaces, claiming robustness and loop detection rather than topological correctness. SHARE: none — this is a citation repair and a scope statement.

**Falsified by.** C1 by the citation not resolving as stated. C2 by AS85 in fact proving a homeomorphism, which would make the precedence far stronger and should be checked before being claimed. C3 by DLTW covering surfaces, same. VACUITY CONTROL: each DOI must be resolved and its title recorded in the finding, not merely quoted from this document.

#### P-167 — registered for R-167, before the harness: How often the hypothesis every correctness theorem needs actually fails

**Ticket:** `R-167` (M). **Records:** `field`, `resolution`, `cells`, `non_transverse_cells`, `non_transverse_fraction`, `gradient_norm_min`, `regular_value_violations`, `overlap_with_ambiguous`, `overlap_with_defects`, `c1_holds`, `c2_holds`, `c3_holds`.

**Hypothesis.** Every PL-correctness theorem — AS85's, Plantinga–Vegter's, Boissonnat–Wintraecken's — needs the isovalue to be a **regular value** and the surface to meet each simplex transversally. This crate has never measured how often that fails, which means it has never known which of its cells are inside the theory and which are outside it. (C1) The non-transverse fraction is measured per field per resolution and is non-zero on at least two fields. (C2) The non-transverse cells **overlap** the ambiguous cells substantially — rank correlation or set overlap reported — which would say the ambiguity problem and the transversality failure are the same problem seen twice. (C3) They also overlap T-001's defect cells, or they do not, and which is stated. SHARE: none — this partitions the cells into in-theory and out-of-theory, which nothing currently does.

**Falsified by.** C1 by zero non-transverse cells everywhere, which would mean every theorem applies everywhere and is a strong and welcome result. C2 by no overlap, which would mean they are genuinely different populations and both need handling. C3 by no overlap with defects, which would say transversality failure is harmless in practice. VACUITY CONTROL: `box_exact`, whose faces have exactly-equal samples, must be included, or the field most likely to violate transversality has been excluded from the census.

#### P-168 — registered for R-168, before the harness: A null registered on purpose — Conley index dies on plateaus, and the census says how big they are

**Ticket:** `R-168` (M). **Records:** `field`, `resolution`, `plateau_cells`, `plateau_fraction`, `isolating_neighbourhood_exists`, `conley_applicable_fraction`, `vs_discrete_morse_applicable`, `c1_holds`, `c2_holds`.

**Hypothesis.** Conley index theory is a robust, computer-assisted alternative to Morse theory that tolerates degenerate critical points, and the corpus holds Dey, Haas & Lipiński on Conley–Morse persistence barcodes at 0.662. It requires an **isolating neighbourhood** around the invariant set, and `box_exact` has exactly-equal samples across whole faces, so `∇f = 0` on a set of positive measure and no isolating neighbourhood exists — the method is undefined precisely where it would be needed. (C1) The plateau fraction is measured per field, quantifying the obstruction rather than asserting it. (C2) The fraction of cells where an isolating neighbourhood exists is reported, and compared against the fraction where discrete Morse theory — already registered — applies, so the two are ranked on the same data rather than on argument. SHARE: none — this closes a candidate.

**Falsified by.** C1 by negligible plateaus on every field, which would revive the method and make this row a false alarm worth having raised. C2 by Conley applying more widely than discrete Morse, which would reverse the ranking. VACUITY CONTROL: `box_exact` must report a non-zero plateau fraction, or the obstruction this row is built on does not exist in the fixture.

---

# Group L — is the input consistent, and is the output well shaped (discrete exterior calculus, geometric measure theory, discrete differential geometry)

#### P-169 — registered for R-169, before the harness: The curl of the Hermite data is the inconsistency `λ = 0.01` is guessing about

**Ticket:** `R-169` (M). **Records:** `field`, `resolution`, `curl_residual`, `curl_residual_normalised`, `cells`, `second_vertex_cells`, `separation_auc`, `lambda_global`, `lambda_per_cell`, `qef_residual_delta`, `sharpness_delta`, `self_intersections_delta`, `curl_ms`, `curl_share`, `c1_holds`, `c2_holds`, `c3_holds`.

**Hypothesis.** `Hodge`, `Helmholtz`, `integrability` and `irrotational` are **zero across every file in this repository**, and the corpus holds the machinery only as flow visualisation. A `HermiteCell`'s `≤12` edge normals are samples of `∇f`, and a true gradient field is curl-free — so the curl residual over those edges is a **direct, per-cell measure of "these normals cannot come from one smooth sheet"**, computable in a few dozen flops on a struct the crate already builds. It has a reason to fire where the crate already hurts: CSG `min`/`max` makes the gradient discontinuous. (C1) The curl residual costs under 2% of extraction. (C2) It separates the cells `M-60` says need a second vertex from those that do not — AUC above 0.8 — on `gyroid` and `fbm_terrain`, the only two fields where `M-60`'s rate is non-zero, with the other five as a zero control. (C3) A per-cell `λ` derived from the residual beats the global `λ = 0.01` on at least one of sharpness or self-intersection count, without losing on the other. SHARE: C1 moves the vertex-placement stage; C3 moves the same stage's quality, not its cost.

**Falsified by.** C1 by above 2%. C2 by AUC below 0.8, which would mean inconsistency and multi-sheet-ness are different phenomena and the residual predicts nothing `M-60` cares about. C3 by the per-cell `λ` losing on either axis, which would vindicate the global constant and is a perfectly good result. VACUITY CONTROL: the five fields with zero second-vertex cells must be included and must report near-zero residual, or C2's separation is measured on a population with no negatives.

#### P-170 — registered for R-170, before the harness: A null registered on purpose — varifolds against the normal cycles already benched

**Ticket:** `R-170` (M). **Records:** `estimator`, `field`, `resolution`, `mean_curvature_error`, `gaussian_curvature_error`, `error_ratio`, `estimator_ms`, `cost_ratio`, `convergence_exponent`, `c1_holds`, `c2_holds`.

**Hypothesis.** The crate's curvature estimation is the normal-cycles line — Sun & Morvan, cited in `benches/experiment_p42.rs` and `p45.rs`. Discrete varifolds (Buet, Leonardi & Masnou, `10.1007/s00205-017-1141-0`, acquired and converted this session) are the rival, with quantitative convergence bounds. This is registered as a head-to-head expecting a null, because a mature estimator is rarely beaten by a newer formalism on the same data. (C1) Varifold mean curvature is compared against the normal-cycle estimator at matched cost on the analytic fields, where exact curvature is known. (C2) The convergence exponent is fitted for both, since a method that is worse at our resolutions but converges faster is a different recommendation from one that is simply worse. SHARE: C1 moves the curvature stage only, whose share must be reported.

**Falsified by.** C1 by varifolds winning, which would be a real result and would redirect the `P-42`/`P-45` line. C2 by indistinguishable exponents, which would make the choice purely about constants. VACUITY CONTROL: both estimators must be run against the analytic curvature on `sphere` and `torus` and reproduce it, or neither is calibrated.

#### P-171 — registered for R-171, before the harness: Fixing sliver triangles without moving a single vertex

**Ticket:** `R-171` (M). **Records:** `field`, `resolution`, `slivers_before`, `slivers_after`, `min_angle_before`, `min_angle_after`, `flips`, `vertex_positions_moved`, `hashes_moved`, `extrinsic_geometry_identical`, `flip_ms`, `c1_holds`, `c2_holds`, `c3_holds`.

**Hypothesis.** Intrinsic triangulations improve connectivity quality while leaving every vertex where it is — the corpus holds Gillespie, Sharp & Crane on integer coordinates at 0.677 and the intrinsic-parameterisation paper at 0.676, and `intrinsic Delaunay`, `signpost` and `Soliman` are zero in the repo. That property is unusually valuable **here specifically**: this crate commits 216 golden hashes on vertex positions, and a method that cannot move a position cannot move a hash. (C1) Intrinsic Delaunay flipping raises the minimum triangle angle materially — at least 10° on the worst decile — on at least four fields. (C2) Zero vertex positions move and zero position hashes move, asserted rather than assumed. (C3) THE HONEST FALSIFIER IS STATED AT REGISTRATION: intrinsic flips leave the extrinsic surface bit-identical, so if every downstream consumer reads extrinsic triangles, nothing measurable changes — and whether any consumer benefits is C3's question, answered by naming one or reporting none. SHARE: C1 moves mesh quality, not cost; C3 decides whether that quality reaches anything.

**Falsified by.** C1 by under 10°. C2 by any hash movement, which would mean the implementation is not intrinsic. C3 by no consumer benefiting, which closes the row and is a likely outcome worth reaching cheaply. VACUITY CONTROL: the worst-decile angle before flipping must be below 15° on at least one field, or there are no slivers to fix.

---

# Group M — connectivity of the sign grid (probability, random fields)

#### P-172 — registered for R-172, before the harness: Caves connect in three dimensions in a way they cannot in two, and it is a theorem

**Ticket:** `R-172` (M). **Records:** `field`, `isovalue`, `components`, `largest_component_fraction`, `component_size_distribution`, `giant_component_exists`, `percolation_isovalue`, `two_d_slice_comparison`, `air_union_find_agreement`, `c1_holds`, `c2_holds`, `c3_holds`.

**Hypothesis.** Duminil-Copin, Rivera, Rodriguez & Vanneuville (`arXiv:2108.08008` / `10.1214/22-aop1594`, acquired and converted this session, 61 pp.) prove that for smooth Gaussian fields in dimension `d ≥ 3` the critical level `ℓ_c(d)` is **strictly positive** — the nodal set percolates in 3D in a regime where the 2D analogue does not. For a voxel game that is a statement about caves: at isovalue 0 the excavated region should have one giant component rather than many isolated pockets. (C1) Sweeping the isovalue on `fbm_terrain` and `noise_cavity`, a giant component appears — one component holding above 50% of the air volume — and the isovalue at which it appears is reported. (C2) The 3D behaviour differs qualitatively from a 2D slice of the same field, which is the theorem's actual content and the control that says the effect is dimensional. (C3) The measured components agree with the crate's existing `Air` union-find connectivity, which is an independent check on both. SHARE: none — this predicts a property of generated worlds, and `M-311`'s union-find is the instrument.

**Falsified by.** C1 by no giant component at any isovalue, which for a cave generator would be a design-relevant negative — it would mean the field produces isolated pockets and connectivity must be authored rather than emergent. C2 by 2D and 3D behaving alike, which would mean our fields are too far from the theorem's Gaussian hypothesis for it to describe them — and `fbm_terrain` is hash-based lattice noise, not a Gaussian field, so this is a real risk and is registered as one. C3 by disagreement with the union-find. VACUITY CONTROL: the isovalue sweep must include values producing both many-small-components and one-large-component regimes, reported as a distribution, or the transition is outside the swept range.

---

<!-- BEGIN BACKLOG BLOCK -->
## Phase 26 — fifty registrations from eighteen fields of mathematics, each before its harness

**Added 2026-08-29, above Phase 25 for the reason every phase goes on top: rule 1 reads top-down.**
Nothing here supersedes Phase 17's or Phase 18's open rows, and Phase 25's twenty bit-packing rows stay
where they are — they are in flight and this phase does not displace them.

**Source: `docs/research/2026-08-29-fifty-experiments-from-unmined-mathematics.md`**, with the vocabulary in
`docs/research/2026-08-29-axes-and-vocabulary-v2.md`. Eighteen mathematical areas probed against the corpus and
the repository at once; fourteen papers acquired and converted; eight recorded as paywalled rather than cited.

**Renumbered on merge.** These fifty rows were drafted `R-103`–`R-152`; Phase 25's bit-packing sweep had
already taken that range, so the whole block moved up by twenty and lands here as `R-123`–`R-172`. Both
source documents carry the shifted numbers; nothing else about them changed.

**The phase's centre is an identity, not a paper.** `b*b - 4*a*c` in `marching_cubes/trilinear.rs:246` is
**identically Cayley's `2×2×2` hyperdeterminant** of the eight corner values — proved symbolically, not sampled.
It is also the discriminant of the pencil `det(A₀ + λA₁)` for **all three** axis pairings, which is the mechanism
behind `M-206`; its sign is a `GL(2)³` invariant, so the body-saddle count **cannot depend on cell aspect ratio**;
and it is invariant under all 48 octahedral relabellings and under negating the field. `R-123` records it,
`R-124` and `R-125` cash the invariances, and **`R-132` exists to stop the overclaim** that the 730-subcase table
follows from it — it does not, and the registration says why in advance.

**Two holes in `CLAUDE.md` close.** `χ` is asserted as unavailable on `gyroid`; it is `−8N³` per `N³` conventional
cubic cells under periodic wrap, because the gyroid's primitive lattice is body-centred and the cubic cell holds
two primitive cells — verified numerically at four resolutions, with Schwarz P at `−4` and Schwarz D at `−16` as
independent checks (`R-138`, `R-139`). And Etiene et al.'s DOI, flagged unverified at
`2026-08-26-audit-and-phase-23-registrations.md:365`, is `10.1109/tvcg.2011.109` (`R-141`).

**The largest paid-for gap is Group D.** Five papers on metric-based anisotropic mesh adaptation sit in the corpus
at 0.70–0.71 and the repository cites none of them. `R-143` registers the honest version in advance: the
convergence **rate** does not improve — the exponent is fixed by polynomial degree and dimension — and what moves
is the **constant M-12 fits**, concentrated on features with a flat direction.

**Order matters in five places and nowhere else.** `R-123` runs before the rest of Group A and `R-127` rides it.
`R-138` runs before `R-139` and `R-141`. `R-142` runs before `R-143`, `R-145` and `R-147`. `R-151` runs before
`R-152` and `R-153`. `R-148` runs before `R-149` and `R-150`. Everything else is independent and eleven rows are `S`.

**Nine rows are expected to return nulls and that is registered rather than hoped:** `R-132` (which exists to
prevent an overclaim), `R-146`, `R-154`, `R-155` (whose null — that Marching Cubes is already order-optimal and
only the constant is in play — is the most valuable sentence available in this phase), `R-159`, `R-162`, `R-164`,
`R-168`, `R-170`. Phase 23's two most useful rows were `✗51` and `✗54`, both of which said *do not build this*.

**Seven candidates were killed during the probe and are recorded so nobody re-searches them:** random matrix
theory and pseudospectra for the QEF (the matrix is normal, so pseudospectra return the singular-value threshold
they were meant to replace); scheduling theory (bounded per-cell output makes the span `O(1)`); geometric
separators for chunking (an axis-aligned cube already attains the `Θ(n^{2/3})` isoperimetric optimum);
Ising/dimer limit shapes (no determinantal structure in 3D); renormalization-group LOD (no critical point);
compressed sensing (deterministic grid, no sparsity basis); Helly-type certificates for the cell clamp (the target
is a 30-flop `3×3` adjugate).

| | Ticket | Size | Blocked by |
|---|---|---|---|
| ☐ | **R-123** | S | — |
| ☐ | **R-124** | S | `R-123` |
| ☐ | **R-125** | S | `R-123` |
| ☐ | **R-126** | M | `R-123` |
| ☐ | **R-127** | M | `R-123` |
| ☐ | **R-128** | M | `R-123` |
| ☐ | **R-129** | M | `R-123` |
| ☐ | **R-130** | S | `R-123` |
| ☐ | **R-131** | S | `R-123` |
| ☐ | **R-132** | S | `R-123` |
| ☐ | **R-133** | S | — |
| ☐ | **R-134** | M | — |
| ☐ | **R-135** | M | — |
| ☐ | **R-136** | M | — |
| ☐ | **R-137** | L | — |
| ☐ | **R-138** | S | — |
| ☐ | **R-139** | S | `R-138` |
| ☐ | **R-140** | M | — |
| ☐ | **R-141** | M | `R-138` |
| ☐ | **R-142** | L | — |
| ☐ | **R-143** | M | `R-142` |
| ☐ | **R-144** | M | — |
| ☐ | **R-145** | M | `R-142` |
| ☐ | **R-146** | S | — |
| ☐ | **R-147** | M | `R-142` |
| ☐ | **R-148** | M | — |
| ☐ | **R-149** | M | `R-148` |
| ☐ | **R-150** | M | `R-148` |
| ☐ | **R-151** | S | — |
| ☐ | **R-152** | M | `R-151` |
| ☐ | **R-153** | M | `R-151` |
| ☐ | **R-154** | M | — |
| ☐ | **R-155** | M | — |
| ☐ | **R-156** | M | — |
| ☐ | **R-157** | M | — |
| ☐ | **R-158** | M | — |
| ☐ | **R-159** | S | `R-158` |
| ☐ | **R-160** | M | `R-158` |
| ☐ | **R-161** | M | — |
| ☐ | **R-162** | M | — |
| ☐ | **R-163** | M | — |
| ☐ | **R-164** | M | `R-163` |
| ☐ | **R-165** | S | `R-163` |
| ☐ | **R-166** | S | — |
| ☐ | **R-167** | M | — |
| ☐ | **R-168** | M | — |
| ☐ | **R-169** | M | — |
| ☐ | **R-170** | M | — |
| ☐ | **R-171** | M | — |
| ☐ | **R-172** | M | — |

---
<!-- END BACKLOG BLOCK -->
