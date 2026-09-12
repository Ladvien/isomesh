# What governs each axis

One paragraph per axis, saying which theorem or mechanism the ledger currently believes governs it, and the
`P-`/`M-`/`✗` numbers that established that belief. **This file is the loop's memory of what it believes.**
`docs/research/2026-09-12-discovery-loop-prompt.md` §3f requires an experiment whose `Surprise:` line is not
`none` to edit the relevant paragraph in the same commit, citing the entry.

**Provenance.** Built 2026-09-12 from `docs/research/2026-08-12-axes-and-vocabulary.md` (axes 1–10),
`docs/research/2026-08-29-phase-27-axes-and-vocabulary-v2.md` (axes 11–14 and the vocabulary upgrades on 3,
4 and 6), and the landed `FINDINGS.md` entries through `M-489`. The prompt as published named
`2026-09-11-phase-30-axes-and-vocabulary-v4.md` Part 4 as the source; that document does not exist in this
repository and neither does the Phase 30 registrations document, so the two that do exist were used instead.
Every number below is quoted from a `FINDINGS.md` entry; nothing here is a new measurement.

**A paragraph with no citation is a belief with no evidence.** If you cannot name the entry, write
*"unmeasured"* and that is itself the finding.

---

## Axis 1 — Domain decomposition

**Governed by: nothing measured — every cell in this crate is a cube, and that has never been the
variable.** The axis has vocabulary (`✗1`, `✗14`, `✗25` all compare *algorithms* over a fixed cubic
partition) but no entry in which the partition itself moves. Two results bear on it without governing
it. `✗107 / M-459` found the reference-field roster *adversarial to global anisotropy*, which is a
statement about what the fields would reward a different partition for. And `✗127 / M-491` splits the
nine fields **6 bottleneck-bound / 3 curvature-bound** by Theorem 3.4's two cases — the split is
exactly smooth-versus-rough, `fbm_terrain` having **no bottleneck at all** (`bottleneck_pairs` **0**)
where `noise_cavity` has **1,857,461**. A partition experiment that ignores that split would measure
one half of its roster.

## Axis 2 — Sign inference and the underlying field

**Governed by: the field's own regularity, not the sampler's cleverness.** `M-250` refines the edge
crossing on the real field and gets **13–15%** on curved fields and **nothing at all** on the CSG one;
`M-247` finds repeated CSG destroys the worst case and leaves the typical case untouched. `✗42`'s C1
falsification is the sharpest version: a reconstruction gain quoted in dB does not transfer to a mesh
error, because the mesh error is dominated by where the surface is, not by how well the samples are
interpolated between.

## Axis 3 — Ambiguity resolution

**Governed by an identity, and it is proved rather than measured: the cell's discriminant is Cayley's
2×2×2 hyperdeterminant.** `M-440` establishes that `bb − 4ac` at `marching_cubes/trilinear.rs:246` *is* the
hyperdeterminant of the eight corner values, so the body-saddle test is a `GL(2)³` relative invariant with a
published real-orbit classification rather than a transcribed quadratic. `M-206` confirms two independently
derived constructions locate the same body saddles to **1.1e-12**. What the identity does **not** buy is
prediction: `✗123 / M-482` measured `refinement_priority_correlation` **exactly 0.000000** against a 0.50
bar on 200 of 200 rows, and the reason is the other half of the same entry — `symmetry_class_count` is
**1**, all eight corners are a single orbit of the 48-element cube group, so C1 holding *is* C2 failing.
`✗11` already retired the folklore that plain Marching Cubes leaks through ambiguous faces.

## Axis 4 — Vertex placement

**Governed by: a small, field-dependent minority of cells, and the crate does not yet know how to find
them.** `M-60` — only two of seven fields ever need a second vertex in a cell, and the rate *falls* with
resolution. `M-486` tried to predict exactly those cells from a discrete-exterior-calculus curl residual:
the residual is **free** (worst `curl_share` **0.008061907**, 0.81% of an extraction against a 2% ceiling,
24 of 24 rows) and it **does not predict** them — of nine rows with a reachable AUC only `gyroid` at 33³
clears 0.8, and `fbm_terrain` at 25³ reads **0.270443196**, worse than chance. The per-cell `λ` it licenses
is still worth something on a different quantity: `noise_cavity`'s self-intersections fall **388 → 195** at
33³.

## Axis 5 — Connectivity

**Governed by: the extractor's own case table; nothing downstream re-wires it, and one attempt proved it
cannot.** `✗3` retired "every interior Surface Nets vertex has four neighbours". `✗125 / M-488` is the
strong statement: **15,391** intrinsic Delaunay flips reach the fixed point on every row, remove **1,570 of
10,678** slivers — and `vertex_positions_moved` **0**, `hashes_moved` **0**,
`extrinsic_geometry_identical` **true**, with `c3_holds` false at **0 of 8** surveyed consumers. An
intrinsic connectivity change is invisible to everything this crate ships.

## Axis 6 — Guarantees and soundness

**Governed by: which theorem you cite, and whether its hypothesis is checked — and the ledger has caught
itself citing the wrong one.** `V-51 / M-483` establishes that Allgower & Schmidt 1985 is **not** a
homeomorphism theorem (`is_homeomorphism` **false**; it constructs a PL manifold along which
`‖H(x)‖_∞ < ε` under a full-rank hypothesis at the seed), and that the corpus's only `true` is Boissonnat &
Wintraecken 2020, thirty-five years later — with `descendant_citations_in_repo` **21** against
`root_citations_in_repo` **4**. `✗124 / M-484` shows the hypothesis these theorems need fails measurably on
this crate's own fixtures: `box_exact` reads `non_transverse_fraction` **0.191406250 / 0.094238281 /
0.046936035** over 784 / 3,088 / 12,304 cells. `✗15` and `✗19` retired the two unconditional-manifoldness
claims. The standing belief: **soundness here is conditional, the condition is transversality, and the
condition is violated on the one field whose answer is known exactly.**

## Axis 7 — Adaptivity and multiresolution

**Governed by: no approximation guarantee at all, and that is now proved rather than suspected.**
`✗121 / M-479` — greedy meshing is not a matroid (`is_matroid` **false** on all 6 rows, minimal
counterexample a `2×2` block, `min_basis` **1** against `max_basis` **4**), and the LOD chunk budget is not
submodular where it matters: `violating_chunks` **92 of 221** on `capped_gyroid` and **155 of 186** on
`fbm_terrain`, a worst chunk whose marginal gains *increase* at the third refinement
(`2.984502e1 | 9.237806e0 | 1.598680e1`). Any future claim of a `(1 − 1/e)` bound for the LOD queue is
unlicensed; the classical home is Wolsey's submodular set-covering line, **named and not read**.

## Axis 8 — Incrementality

**Governed by: how few chunks an edit touches, and the number is small.** `✗41` — **9.45×** on the 13
chunks where a brush actually matters. The axis's open question is not the speedup but the *settling*: the
ledger has no entry on how many frames a publish takes to become visible to a consumer.

## Axis 9 — Execution model

**Governed by: what the device advertises, which is not what the documentation says.** `M-147` — Bevy's
device already has mesh shaders enabled, with no configuration and no `unsafe` anywhere in this repository;
`V-23` resolved `CLAUDE.md`'s Metal contradiction by finding both sources right about different layers (the
feature reaches Metal, the WGSL compiler does not). `✗28` is the memory-side counterpart: the 128³ penalty
is a property of the **access pattern**, not of the stride.

## Axis 10 — Consistency across chunks

**Governed by: the float backend, and it was verified rather than argued.** `M-31` — `libm` delivers the
bit-identical cross-platform meshes it was chosen for, 216 golden hashes generated on macOS/arm64 passing
unchanged on Linux/x86-64. `✗18` retired "a hairline seam difference is a crack no weld can close".

## Axis 11 — Reconstruction filter and approximation order

**Governed by: Strang–Fix, and the crate is at order 2 with both filters it has tried.** `M-12` measured
`L = 2`; `✗120 / M-477` held C1 on **8 of 8** fields with `order_gap` spanning **−0.072307** to
**+0.128478** against a −0.15 bar, *both* filters order 2, and the derivation executed to `affine_residual`
**4.441e-16**. C2 there is the cautionary half: it was **arithmetically unreachable**, `c2_margin_db`
topping out at **0.000000478** against a required **0.001000**, because the comparand's non-control arm was
already the same lattice and filter — `x > x` is false whatever `x` is.

## Axis 12 — Sampling lattice and quantization efficiency

**Governed by: the table already being in L1, which removes the motive for compressing it.**
`✗122 / M-480` — `eval_ms_ratio` **953.694027** on the dense bits, `spectral_concentration` **0.243408203**
and **0.428710938** against a 0.9 bar, `spectral_terms` **256 of 256**: the two informative bits of the
shipped table are dense, its twelve edge bits are each one degree-2 coefficient, and a 256-entry table costs
nothing to read. `M-481` is the registered null that held — `closed_form_max_gap`
**0.000000000000000**, stability computable a second way to **6e-15**.

## Axis 13 — Element shape and anisotropy

**Governed by: the reference-field roster, which is adversarial to it.** `✗107 / M-459` — at matched
Hausdorff the metric-driven arm is **never** cheaper (`c1_winners` **0 of 4**) and is **2.298935×** dearer
on `thin_plate`, with `metric_share` **0.332862 to 8.747515** against a 0.15 bar on 40 of 40 rows. The
second half of that entry is the standing constraint on every accuracy claim in this file: **C2 was
`unmeasurable` on 20 of its 40 rows, because `validate::accuracy` is meaningless where `field.bound()` is
not `Exact` — four of the eight reference fields are ungradeable.** `✗125 / M-488` adds that flipping does
not rescue the angles: worst-decile minimum angle rises by at most **1.642437°** against a 10° bar.

## Axis 14 — Optimality and lower bounds

**Governed by: one denominator that exists and is closed, and one that does not exist at all.** The closed
one is greedy meshing: `M-478` reproduces Eppstein's optimum `n/2 + h − g − components` against exhaustive
search on **33,267** components with `brute_disagreements` **0**, and the shipped mesher's worst distance
from it is `ratio` **1.027132** — **1.120%** over optimal, exactly optimal on 8 of the 21 rows where the
ratio is defined, while computing the optimum is *cheaper* than running the mesher on 9 of 24 rows
(`cost_ratio` down to **0.424099**). Greedy meshing is finished. The missing one is **accuracy**: no entry
states an `n`-th minimal error for this crate's field class, so every accuracy number in `FINDINGS.md` is a
numerator. Information-based complexity was ranked #5 of the six transfers in axes-v2 Part 4 and **was never
registered here** — Phases 28–30 did not land on this machine (`P-176` is the highest registered `P-`).

**Reach is not the missing denominator for resolution either — `✗127 / M-491`, 2026-09-12.** The
first candidate for Axis 7's and Axis 14's missing denominator was the reach, on the strength of the
sampling theorems `V-51` sorted out. Measured, `h*/τ` — the coarsest spacing whose extracted topology
matches the finest rung's, over the reach — spans **0.250000000** to **617.006450883** across the
roster, a spread of **2468.025804** against a registered bar of **2.0**. Across the four *smooth*
fields alone it is **3.339**, which is the usable half of the result. The mechanism is that a rough
field's reach is minute at the scale it is actually smooth (`gyroid` **0.000709069**, `noise_cavity`
**0.000145884**) while its topology settles hundreds of times coarser. **The denominator for accuracy
is still unclaimed, and reach is not it for resolution.**

---

## Instruments the loop grades with

Not an axis — the oracles themselves, because a `Surprise:` line is only as good as the instrument that
produced it. `M-487` established that at matched cost the **normal-cycles** incumbent beats discrete
varifolds on mean curvature on **28 of 30** rows (`error_ratio` **1.000821 to 3.358471**), with one
exception that is a C⁰ seam artefact and one genuine reversal on the `torus`'s **Gaussian** curvature
(1.59× better at 129³, exponent **1.575152** against **1.127950**). `M-485` found the registered obstruction
to a Conley-index instrument real but one dimension lower than claimed — `plateau_fraction` **0.500000** on
`box_exact` with `constant_cells` **0**, yet `conley_applicable_fraction` never below **0.957031**.
`M-489` gives the loop a percolation instrument: on `noise_cavity` the giant component appears at
`percolation_isovalue` **0.029340**, and at `iso = 0` one component holds **0.814493** of the air against
**65** isolated pockets — with the theorem's qualitative 2D/3D difference present on `noise_cavity` and
**absent** on `fbm_terrain`, whose 2D slices percolate too at **7.272575**.

**A digital-topology χ is ground truth about the sampled grid, not about the field — `✗126 / M-490`,
2026-09-12.** `M-458`'s mesh-free oracle is exact and reproducible: `✗126` re-implemented it, reproduced
`M-458`'s `noise_cavity` surface χ of **82** at 33³ and **−152** at 65³ *on the nose*, and read **2** on
`sphere` and **0** on `torus` from the same code. What does not follow is that any of those numbers is the
field's χ. On one fixed realisation of `noise_cavity` the three rungs read **82 / −152 / −186**, still
moving at 129³ (`chi_convergence_gap` **2.853658537** then **0.223684211**). **Any entry quoting a random
field's χ must quote the resolution it was read at, and must not call it ground truth without a
convergence column.**

**There is now an analytic expectation to grade a random field's χ against, and it works — on a field
that is actually Gaussian.** Taylor's Gaussian kinematic formula, instantiated on the flat 3-torus as
`λ^(3/2)·Vol·(2π)^(−2)(u² − 1)e^(−u²/2)`, matched an exactly-Gaussian fixture on **5 of 5** rungs at
`relative_departure` **0.004342748 / 0.010112782 / 0.004829332**, including the two rungs where it predicts
**exactly zero** for any `λ`. Of the shipped noise fields, `fbm_terrain` lands inside a matched band at
**0.25σ** and `noise_cavity` misses by **6.77σ** with the *sign* of χ wrong (**15.000000000** measured
against **−6.602285203** predicted). So the formula is an oracle for `fbm_terrain` and is not one for
`noise_cavity`, and the open question is whether the difference is Gaussianity or **stationarity** — the
cap `noise_cavity` carries is not something a stationary-field theorem allows.
