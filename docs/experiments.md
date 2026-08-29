# Experiments that held, and what they buy you

**One hundred and twenty-two times**, this project wrote down what it expected to measure — and what result
would prove it wrong — **before running the measurement**. One hundred and fifteen of those are enforced by
the compiler: `crates/isomesh/src/experiment.rs` refuses to build a harness whose id is not registered, so
the prediction cannot be edited after the numbers arrive. The other seven are `P-1`…`P-7`, which predate the
macro and live in `FINDINGS.md` as prose. **These two numbers are counted from `experiment.rs` rather than
remembered**, because nothing gates this paragraph and it has rotted before: it read *"one hundred and
two"* and *"ninety-five"* for a whole phase after they stopped being true.

This page is the scorecard, and the falsifications are the reason it is worth reading. The most useful
thing in the whole record is still a bug that made two true hypotheses look false for a day.

If you only read one section, read [the instrument story](#when-a-measurement-is-impossible-suspect-the-instrument).

---

## Why pre-register at all

A prediction written after the number is known is not evidence, and it reads exactly like one written
before. This project has already caught itself writing expectations into docs that the measurement
then disproved — twice, at ✗1 and ✗3.

So the prediction moved somewhere the compiler can see it:

```rust
use isomesh::experiment;

// This line does not compile if "P-8" is not in `PREREGISTERED`.
let p = experiment!("P-8");
assert!(p.hypothesis.starts_with("A weld gated"));
```

Three properties make it a gate rather than a habit:

- **An unregistered id is a compile error**, not a test failure you can ignore.
- **`Preregistration` is `#[non_exhaustive]`**, so nothing outside the crate can construct one. The
  macro is the only door.
- **`falsified_by` is not optional.** There is no way to express a registration without it, which is
  the point of the field existing: *a hypothesis with no falsifier is not a hypothesis.*

Registering is a commit, so git carries the ordering that prose cannot. `scripts/backlog_gate.sh`
checks that every registered id reaches `FINDINGS.md`, so the prose can elaborate and cannot
contradict.

| era | how it was registered | held | falsified | never ran |
|---|---|---|---|---|
| P-1 … P-7 | prose in `FINDINGS.md` | 5 | 2 | — |
| P-8 … P-17 | `crates/isomesh/src/experiment.rs` | 5 | 4 | 1 |
| P-18 … P-21 | `crates/isomesh/src/experiment.rs` | 3 | 1 | — |
| P-22 … P-37 | `crates/isomesh/src/experiment.rs` | *not tallied here — several split their verdict across clauses, and a column that rounds those to one word would be the kind of summary this page exists to distrust* | | |
| P-38 … P-47 | `crates/isomesh/src/experiment.rs` | 3 | 7 | — |
| P-48 … P-56 | `crates/isomesh/src/experiment.rs` | *split, and tallied honestly below: 4 experiments held every clause, 2 were falsified outright, and 3 split — the phase's most useful result is a clause that failed* | | |
| P-57 … P-60 | `crates/isomesh/src/experiment.rs` | *4 registrations, every one falsified on at least one clause and none on all — `FINDINGS.md` carries them (`✗39`–`✗42`); they are not summarised on this page* | | |
| P-61 … P-72 | `crates/isomesh/src/experiment.rs` | *12 registered, **10 ran and 2 are blocked on a paper nobody can get** — tallied below. 1 held every clause, 2 were falsified on every clause, and 7 split* | | |
| P-73 … P-102 | `crates/isomesh/src/experiment.rs` | *30 registered and **all 30 ran** — tallied below. 6 held every clause, 2 were falsified on every clause, and 22 split, so **24 of the 30 carry at least one falsified clause*** | | |
| P-103 … P-122 | `crates/isomesh/src/experiment.rs` | *20 registered and **all 20 ran** — tallied below. Sixty clauses: **34 held, 24 falsified, 2 vacuous**. 2 rows held every clause, 1 held every clause it could fire, and **17 of the 20 carry at least one falsified clause*** | | |
| P-123 … P-126 | `crates/isomesh/src/experiment.rs` | *4 registered and **all 4 ran** — tallied below. Twelve clauses: **9 held, 2 falsified, 1 vacuous**. 2 rows held every clause, 2 carry a falsification, and **1 landed a source change, registered in advance*** | | |

A page called "experiments that held" which quietly dropped the falsifications would be exactly the
failure this machinery exists to prevent. They are all here.

**Phase 19 is the worst era on this page and the most useful one.** Seven of ten falsified, and four of
those seven were falsified by their *own* registered falsifier naming the mechanism in advance. Two were
falsified by an error in the registration itself — P-43's cost clause was wrong arithmetic and P-42's
Gaussian clause was a tautology — which is the machinery catching the person operating it, and is the
strongest evidence on this page that the practice is doing something.

---

## The scorecard

| prediction | verdict | what it means if you depend on this crate |
|---|---|---|
| **P-11** — one canonical `world_of_sample` closes seam cracks at every cell size, not just powers of two | **held**, both clauses (M-278) | Chunk seams are exact **only at power-of-two cell sizes** today. Canonical reconstruction gives 0 unmatched seam edges at every spacing tried; what the crate can offer now gives 0 only at powers of two and 63–348 at `0.1`, `1/12`, `1/14`. The weld hides it; an unwelded consumer gets it in full. The fix is priced and open as X-005 |
| **P-15** — more than half the dual mesher's cycles are in `emit_quads` | **held**, and by more than it claimed (M-284) | 45% of the instructions and **82% of the time**, at `r² = 0.995`, predicted to +0.12% on a held-out shape. This is the measurement that made the 4.26× possible |
| **P-9** — k-way weld order changes the output | **held** | The weld is **not confluent**, so vertex order is pinned deliberately rather than left to `HashMap` iteration. Your meshes reproduce byte-for-byte |
| **P-13** — the non-convergent angle at a sharp feature is bounded below by the feature's dihedral | **held** — after M-289 | It is a **property of sharp edges**, not a defect with a location. Nothing to fix; something to know |
| **P-16** — every past-90° vertex sits on a cell straddling the crease | **held, at 0%** — after M-289 | Same conclusion from the other side: two faces meeting inside one cell, not a winding bug |
| **P-8** — a weld gated on the link condition reaches zero non-manifold output | **falsified**, both clauses (E×4) | **The gate makes things worse.** The crate does not ship one, and now you know why rather than wondering |
| **P-12** — the dual's superlinearity is a cache-miss effect | **falsified** (M-279) … and later true | The best methodological lesson here. See below |
| **P-14** — the residual non-manifold vertices are one-vertex-per-cell meeting two surface components | **falsified** on clause one (M-275) | The vertex-to-edge ratio came out **exactly 2:1 on every row** — 597/314, 99/48, 38/19, 106/53. There are not 597 defects plus 314 defects; there are 314, each counted twice more by the vertex walk |
| **P-17** — Manifold Dual Contouring's residue is an interior ambiguity | **falsified** (M-291) | `Interior::Joined` fires on **100%** of ambiguous-face pairs, offenders and control alike — no discriminating power at all. Bounded exhaustively instead (M-292), then reduced to a 48-sample fixture (M-294) |
| **P-10** — vertex inflation from gate-plus-split is under 1% | **never ran** | Superseded. P-8 killed the gate P-10 presupposed before P-10 could measure its cost |
| **P-19** — the winding number's on-demand crossover is set by row sharing, so `Q*` is of order `N²` not `N³` | **held**, both arms (M-303) | `Q*/N²` came out **1.13–1.43** at 17³/33³/65³ against a registered 0.5×–4× band. The naive reading overstates the win by a factor of `N`: at 65³ the real crossover is **4,791 queries against 274,625 samples**, so an on-demand field pays below about **1.7% of the grid**. The hole-punched control moved `Q*/N²` to 13–42, and per-point cost tracked boundary-edge count to **×15.7 against ×16.2** — the mechanism confirmed by number, not by direction |
| **P-20** — a weld key moves no topology metric beyond the splits it names | **held on both falsifiers, wording amended** (M-305) | A constant key moves nothing on **0 of 51** configurations, and non-manifold vertices never rise by more than the split count — so E×4's manufactured bowties do not come back. But boundary edges go **0 → 24** on a creased cube, which H said would not happen. That is not a defect: it **is** the split, seen from the edge column. The exact mirror of E×4, where a *pairwise* refusal did its damage in the vertex column |
| **P-21** — a freshly extracted mesh separates exactly the sample pairs the field's own sign separates | **held**, both clauses (M-307) | **Marching Cubes, MC+decider and Marching Tetrahedra seal 24 of 24** — eight fields × three resolutions, 9.15 M probes, zero holes and zero membranes. **All three duals fail `fbm_terrain` with the identical count** (92 / 138 / 190) and **every one of those holes is on a face of the sampled domain**: a dual's quad needs all four cells around a sign-changing grid edge and at the domain face only one or two exist. **For a chunked world that face is the chunk seam**, so a dual chunk's collider is not watertight on its own. The property that splits the family is not primal-versus-dual but *does the method put its crossing on the probed grid edge* — subgrid Marching Tetrahedra is primal and scores worst, 18 of 24 |
| **P-18** — every convex-decomposition precondition is already reported by `ColliderReadiness` | **falsified**, and not where predicted (M-300) | Registered *expecting* to die on self-intersection-freedom. It died on **non-manifold vertices** instead, and self-intersection-freedom turned out not to be an input precondition of anything (✗20). The audit was the product either way |
| **P-61** — a crossing stored as an offset from the edge midpoint is bit-exactly antisymmetric, so a mirrored grid gives a mirrored vertex | **falsified** on C1 and C2, **held** on C3 and C4 (✗49 / M-372) | **Your vertex positions moved in 0.0.10, and in exchange plain Marching Cubes is now equivariant under all 48 octahedral elements instead of 6.** 0 mismatches on 5,800,000 straddling pairs against 1,035,808 for the old lower-corner form. 135 of 216 golden hashes were rebaselined; triangle counts, Hausdorff and self-intersections are unchanged to twelve digits, and 2,285 of 28,124 cut edges — **8.12%** — moved by ≤ 268 ULP. **If you hash meshes, 0.0.9's hashes will not match.** C1 and C2 were cost clauses and both were priced wrong at registration |
| **P-63** — `O-12` is finite at 2¹⁸, and one sweep settles the vertex-link question for Marching Cubes | **held** on C1 and C2; C3 falsified by its own block (M-374) | **The oldest open question in the ledger is half closed, by exhaustion rather than by sampling.** Every face incident to an edge vertex comes from one of the four cells sharing that edge — 18 corners, 2¹⁸ = 262,144 patterns, the *whole* space — and not one produces a vertex with a split link, on both interior-ambiguity settings at four magnitude draws. The pre-fix defect reproduces **5,302 times** in the same walk, so the zero means something. **What is still open is the dual family**, which needs 4³ = 64 corners and is filed as `R-072` |
| **P-64** — bounded model checking proves the combinatorics the property tests only sample | **held**, all three (M-379) | **The case table is proved indexable, not spot-checked.** All 256 sign patterns in 2.04 s and all 16,384 (pattern, mask) pairs in 161 s. The useful half is the failure on the way: pointing Kani straight at the constructor ran CBMC out of **32 GB**, and the proof only exists because the work was split between `const` evaluation over the real constructor and a model checker over the consumer. A `#[kani::should_panic]` control corrupts one triangle and reports the failure, so a vacuous proof cannot pass as a real one |
| **P-68** — a running error bound turns a vertex position into a position *and* a certified interval | **falsified** on C1 by one unit; **held** on C2, C3 and C4 (✗55 / M-382) | **The bound is sound at coefficient 3 and unsound at the registered 2** — 10 violations in 19,415, which is the kind of margin no sample size finds by luck. The ground truth is exact: `f64` samples are dyadic rationals, so the true crossing is an exact rational and the error is computable in `i128` with no floating point in the loop. That is what makes this a soundness statement rather than two approximations compared |
| **P-69** — restructuring the sample loop autovectorises, with bit-identity as the gate | **falsified** on C1, **held** on C2, C3 **vacuous** (✗51 / M-375) | **Nothing widened: total `%ymm` across all eleven monomorphisations is zero**, and the prescribed shape is a 3–7% *regression*, so it was reverted. C1's 2× was arithmetically unreachable from the start — the loop is **11.6%** of the marginal extraction cost, so its ceiling is 1.06×. The durable result was in no clause: there were **three** copies of that loop and they are now one (`sdf::sample_grid`), with the original arithmetic and all 216 hashes unchanged |
| **P-70** — subgroup ballot compaction beats the Hillis-Steele scan | **falsified** on C1 by arithmetic done at registration; **held** on C2 and C3 (✗54 / M-381) | **The mechanism is real — 1.4114× at 2 M elements, bit-identical, reproducing the literature's 1.33–1.49× band — and it is not shipped.** Applied to the path it lives in it moves **1.27%**, because `scan_ms` is 0.3657 of 8.3694 while `upload_ms` is 87.50%. A second WGSL path for 1.27% is what the one-path rule forbids, and naga 29.0.4 rejects `enable subgroups;` outright — and that 87.50% is the sample-grid upload on a path GPU-011a has since deleted (✗73 / M-405) |
| **P-71** — the 83% is a blocking round-trip, and both targets can avoid it | **falsified** on all three (✗52 / M-376) | **The largest piece of a GPU extraction is the geometry copy, not the count wait** — 0.6663 ms of 1.1331 at 129³ against map-wait's 0.3191 — so `extract_indirect` removes **32%**, not the registered 60%, and it removes it by *not delivering the bytes*. For a consumer that only draws, that is 100% of what it needs. For a collider consumer the copy is unavoidable, which is why `DeferredGeometry` ships in 0.0.10: measured under the frame-budget scheduler at **1.41 frames of latency, worst case 2** |
| **P-72** — the granularity of the active-cell structure is a parameter with an optimum | **held** on C1 at 51×, **falsified** on C2 and C3 (M-377) | **Chunk size is the largest single knob in this crate and its optimum is interior at 4³.** Eleven edits cost **4.59 ms at 4³ against 432.73 ms at 64³** at a fixed total cell count, because remesh is `dirty_chunks × chunk_cells³` — 25,344 cells re-meshed against 8,126,464. The registered 8³–64³ range would have given the wrong shape entirely: it is monotone there with the minimum at the smallest granularity swept. **The demos are not retuned on this**, because the fixture excluded per-chunk entity and draw cost — and P-89 has since measured below it, at **4.36×** and **1.76×** against 1³ (✗69 / M-401) |
| **P-62** — the Plantinga–Vegter certificate is sound against this crate's own tunnel classifier | **held** on C1 and C2, **falsified** on C3 (M-378) | **Zero unsound certificates over 2,389 tunnel cells**, from a predicate that refuses **95%** of the population — so the zero is not the trivial one. It is not free: the cost is 21%, not the registered ceiling, and the reason is structural rather than fixable in isolation — a standalone pass re-gathers the eight corners the extractor already has in registers — fused into that gather it costs **0.0553** of extraction and certifies the same set (M-410) |
| **P-66** — a derivative sign test is a usable under-resolution witness | **falsified** on all three (✗53 / M-380) | **The third attempt at this line died too, and it died informatively.** The witness is sound at k = 17 samples and unsound at k = 5, C2's second half was stated backwards, and `thin_plate` ranks **last** on the metric C3 predicted it would top. Two earlier attempts (P-43, P-44) were *value* witnesses; this was a *derivative* witness, and the family is now three for three. Nothing here changes what the crate does |

---

## Phase 19 — three that held, seven that did not, and two the registration itself got wrong

| prediction | verdict | what it means if you depend on this crate |
|---|---|---|
| **P-39** — a brush that provably cannot win the min/max chain inside a chunk can be deleted from the tape, bit-exactly | **held**, all three clauses (M-341) | **Your edit history stops being a tax.** Median 21 of 64 brushes survive per chunk, a chunk meshes **3.4× faster**, and the mesh is byte-identical on 64 of 64 — because IEEE `min`/`max` *select* an operand rather than computing one. The bound costs 1.5 × 10⁻⁴ of what it saves. `SmoothAdd` is excluded from the losing-direction prune and that exclusion is load-bearing: prune it and 48 of 64 chunks change |
| **P-40** — the active-cell test is one bit, so 64 cells decide at once | **demoted, then re-earned as a count** (M-337 → M-348 → M-349) | Registered as a ≥1.25× wall clock and recorded at 1.336×; an external audit found the committed CSV said **1.1925×**, and three quiet-machine re-runs read 1.022× / 1.184× / 1.177×. The mechanism was never in doubt — what was wrong was gating it on a stopwatch, which ✗24 had already forbidden. Re-registered as P-50 and held as **three exact equalities on 15 of 15 rows**: gathers performed equals cells on the scalar path and active cells on the bitmap path, word groups equals `cells_x.div_ceil(64)·cells_y·cells_z`, and the two produce the same **ordered** active list. At 128³ it removes **99.07%** of gathers on `sphere` — an integer, identical on every machine |
| **P-41** — the sign lattice is not well-composed, and that is where the duals go non-manifold | **held**, and the relation is a bijection (M-338) | **You can count the defects before the mesh exists.** `critical cells == non-manifold vertices`, exactly, on every affected field — 602, 141, 58 — and exactly zero on the five clean ones, from both directions. Co-location 2442/2442 against a ~1% chance baseline. The census is a function of the sign bytes alone, so it is available at 13.5 ms against 63 ms to mesh the same grid |
| **P-46** — repairing the sign lattice drives the duals to zero non-manifold output | **falsified**, and the conditional it leaves is the product (M-344) | The repair is **free and total where critical cells are isolated** (58 → 0 in one sweep, 141 → 0 in five, Hausdorff ratio 1.000000) and **provably stuck where they are dense**: an exhaustive walk over all 4,096 face-adjacent sign patterns finds **36 where no corner leaves both cells clean**. Gate on maximum cluster size, not on count |
| **P-47** — composed fields have a real gradient hole | **falsified** in the mean, held on speed (M-345) | The central difference is **7.6 × 10⁻⁵ degrees** off on a 64-brush stack — three orders under the bar, so your normals are fine. But it reaches **4.365°** exactly where a vertex's stencil straddles a CSG seam, and exact gradients are **2.8× cheaper** anyway. Justified by cost, not by accuracy, which is not the reason anyone expected |
| **P-45** — normal-cycle curvature measures are additive, so a per-chunk number composes | **falsified** on two of three clauses (M-343) | The two measures have **opposite** defects. The Gaussian one is chunk-local and does not sum — it overshoots by exactly `π × excess chunk incidence`, the missing `−N(A∩B)` term, because a partition of faces is not a partition of space. The mean one sums to 3 × 10⁻¹³ and needs a one-ring to compute |
| **P-42** — the curvature residual falls inside a bound computed from the mesh | **falsified** on convergence; the clause was a tautology (M-340) | With `B` the whole closed surface, `Σ(2π − α_v) = 2πχ` **identically**, so the residual is one f64 epsilon per vertex and grows ×4 per halving because the vertex count does. Theorem 6's bound does not converge here either: `max circumradius` is not `O(h)` because Marching Cubes' needles get relatively worse with resolution |
| **P-43** — the max cell-centre residual witnesses under-sampling | **falsified**, and the cost clause was wrong arithmetic (M-339) | One `C¹` crease cell in 2.1 million pins the maximum at every resolution, so it does not fall at all while the Hausdorff does. And *"the structural ⅛"* was per-cell corner accounting; this crate prefills one shared grid, so the real cost is `(n−1)³/n³` — a doubling of the sampling work, not an eighth |
| **P-44** — the *mean* residual witnesses it instead | **falsified out of sample** (M-342) | The correlation reproduced on all four untouched fields at `r ≥ 0.98` and the line still dies: the exponent gap **is** the `/h` normalisation. A mean's convergence order is field-independent because it is set by the smooth bulk; a supremum's is field-dependent because it is set by the worst feature. No power of `h` reconciles them, and an `r` of 0.98 over four monotone points is a statement about monotonicity |
| **P-38** — Marching Cubes carries the same 128³ aliasing defect the dual path was cured of | **falsified** at 0.98 against 1.5 (M-336) | **A-024's `size[0] \| 1` is correctly scoped and now that is measured.** The 64 KiB plane stride is not the defect; the *walk* is. Marching Cubes streams eight corners and steps `x`; `emit_quads` walks every grid edge on three axes with endpoints a row and a plane apart. A future extractor inherits the pad only if it inherits the walk |

**Two of these were falsified by the registration rather than by the world**, and that is the part worth
keeping. P-43's cost clause was arithmetically wrong when it was written; P-42's Gaussian clause was
unfalsifiable by construction, because taking `B` to be the whole closed surface turns the measurement
into discrete Gauss–Bonnet. Neither error survived contact with a number, and neither was repaired by
amending the prediction — both were re-registered under new ids (P-44, P-45) and both of *those* then
failed too, honestly, on data they had not seen.


---

## Phase 20 — where three of six were wrong before anything ran

Every registration in this phase was checked against the paper it cites **before** the harness was
written. Three of six rested on a claim the source does not make. That is the cheapest error this project
has ever caught, and it is worth more than any single verdict below.

| what the registration cited | what the source actually says |
|---|---|
| *"edge Chamfer 0.0262 vs MC 0.417 — 13× on sharp features"* | **No such figure exists.** The tangency paper measures Hausdorff, Chamfer and its own energy, treats sharp features qualitatively, and its residual is a rank-3 point-to-point term, not the point-to-plane form that was registered |
| Custódio's `=` label reduces degenerate triangles by *N* | The paper reports **no degenerate-triangle count on any dataset**, and its triangulator is a per-cube convex hull explicitly *"without the need of a look up table"* |
| Affine arithmetic is *correlation-aware*, with a measured tightening | The paper contains **no correlation argument, no tightening figure, and no `min`/`max` rule at all** — and that last gap turned out to be the experiment |
| A monotone-edge certificate for 3D meshes | **The theorem is real and it is two-dimensional.** Its pigeonhole step is *"a triangle has only three edges"*, and it does not apply to the trilinear interpolant |

| prediction | verdict | what it means if you depend on this crate |
|---|---|---|
| **P-54** — a correlation-aware bound rejects empty cells a Lipschitz ball cannot | **held**, all three (M-354) | **A field whose global Lipschitz constant is worthless becomes cullable.** `gyroid` 688 → **2,647** of 4,096 cells rejected, **3.85×**, with 1,959 gained and **zero lost**. On `gyroid_uncapped` Hart's test rejects **nothing at all** and affine rejects 1,098. And where the field is built from `min`/`max` the gain is *exactly zero* — `box_exact`'s two rejected sets are identical cell for cell — because the collapse to an interval destroys the correlation. Mesh byte-identical on 12 of 12, with a deliberately sabotaged bound proving the gate can fail. Costs ~12× the arithmetic: a capability, not a speedup |
| **P-56** — a central difference across a CSG seam is wrong by at most half the crease angle | **held** on C1 and C3 (M-350) | **The error has a closed form and it is attained.** `(180° − θ)/2`, with median tightness **0.9748** over 24 rows and 9 rows above 0.99. It is provable, not fitted: the difference returns `u − Λw` for a diagonal `Λ` with entries `≤ ½`. P-47's lone 4.365° outlier back-solves to a 171° seam against a 4.5° bound — explained, not excused. It does **not** shrink with resolution, because the stencil step does not |
| **P-53** — a third corner label removes the degeneracy integer volumes create | **held**, **split**, **held** (M-352) | **Every** degenerate triangle traces to an exactly-equal corner — 164 of 164 on `fuel`, 58,097 of 58,097 on `bonsai` — and the repair moves **no geometry at all** (`max_snap_distance` is exactly 0). But it is only safe on one of the two volumes: on `bonsai` it welds **520 separate components**, because 516 of 17,201 collapse groups join vertices that share no triangle. **The shippable result is the precondition** — one union-find over the baseline mesh — not the repair |
| **P-51** — this crate's extractors violate the empty ball a distance sample asserts | **held**, **falsified**, **held** (✗38 / M-355) | **They do not, and that closes four papers in one pass.** Marching Cubes pierces 0 of 22,790 vertices and Dual Contouring 0 of 22,798, exhaustively, two margins negative. The *other* half is wide open: **2.9%–80.5% of samples have a sphere no vertex ever touches**, with MC worse than DC on 5 of 5 fields and **16×** worse on `thin_plate`. The worst misses are exact constants — `√2` cells for MC on a box, because MC places no vertex on a box edge at all |
| **P-52** — tangency placement beats the QEF on sharp features | **falsified** on three of four (✗37 / M-353) | **64.8× and 87.2× worse**, and the algebra says why in one line: the least-squares point is `(1 − λ)·mean(pᵢ) + λ·c`, and the mean of a cube's eight corners is the **cell centre**. Confined to a cell, the tangency energy stops being about the surface. This confirms the source's own ablation — globality is its mechanism, not an implementation detail |
| **P-55** — monotone mesh edges certify the critical-point structure | **falsified** on all three (✗36 / M-351) | **The zero was unreachable by geometry, not by sampling.** A mesh edge is a *chord*: on a convex field `f` runs `0 → negative → 0`, so `g(0) = r(cos θ − 1) < 0` and `g(1) = r(1 − cos θ) > 0` for every pair of surface points. `sphere` flags **804 of 804** edges against a registered zero. 97.5% of flags were decided by the endpoints before any interior sample existed. The predicate belongs on the *ambient* complex |

**What this phase says about the practice, which is not the same as what it says about meshing.** Four of
the six registrations contained an error I put there: a bar copied from a figure that does not exist, a
tolerance that scales by a quantity which is identically zero where it is applied, a curve fit where a
scale-free equality was available, and a field list naming five members of a set with four. **Every one
was caught by the artefact rather than by review** — and three of them were caught before the harness
existed, by reading the paper. The registrations that needed no external claim, P-51 and P-56, were the
two that could be written immediately and the two that needed no correction.

---

## Phase 23 — twelve registered, ten run, and two blocked at the door

All twelve came out of a single audit of this project's own practice, and the audit's central finding
turned up again *inside* the phase that exists to apply it: **a clause stated as a ratio of a total must
name the share of that total it can move.** P-69 registered a 2× on a loop that is 11.6% of the cost it was
denominated in, so its ceiling was 1.06× before anything ran. P-70 registered a speedup on a stage that is
4.4% of its own path. Both numbers were available at registration and neither was computed.

**Two registrations are blocked, and they are blocked on acquisition rather than on a ticket.**

| blocked | what it needs | why no ticket can unblock it |
|---|---|---|
| **P-65** — does MCPro's procedural construction resolve `UnresolvedSixSaddle`? | Stahl & Grosso, GRAPP 2025, `10.5220/0013309800003912` | Six acquisition routes were tried. What is in the library is a **383-character SciTePress landing page** that the catalog reported as *"converted, embedded"* — so the pipeline's own status flag is not a presence oracle either (`M-371`). Running it would mean inventing the quadrant subdivision, halfedge assembly and third routing from an abstract |
| **P-67** — does reduced affine arithmetic keep P-54's rejection rate in a fixed-size struct? | Knoll et al., `10.1111/j.1467-8659.2008.01189.x` | C1 and C2 are runnable against this crate's own P-54 baseline. **C3 is the informative clause and it is not** — it exists to reproduce a measured 1.5–2× / 3–4× band *and* a superquadric case where intervals win, and reproducing a band from a summary is exactly the failure `✗21` records. Splitting the registration to run the easy two would measure that the method is cheaper without ever testing where it is not |

Neither is on the scorecard above. **A scorecard row implies a verdict**, and these have none; they are
registered in `crates/isomesh/src/experiment.rs` so the compiler holds the ids, their questions are on
record in `FINDINGS.md`, and the unblocking event is external — a PDF appearing.

**Three of the ten changed what ships.** P-61 moved every vertex position and rebaselined 135 of 216
golden hashes, to buy bit-exact equivariance under all 48 octahedral elements. P-71 landed
`DeferredGeometry` as a third extraction contract. P-72 measured the largest single knob in the crate — a
**51×** spread in edit-plus-remesh across chunk granularity at a fixed cell count — and deliberately
retuned nothing, because its fixture excluded the per-chunk entity and draw cost a real demo pays.

**One found a crash rather than an answer.** P-63's exhaustive sweep panicked inside `marching_cubes`, in
release, on an ordinary trilinear cell: `MAX_PATCH_TRIANGLES` was a **sampled** maximum of 24 where the
triangulator's own buffer was the **derived** 40. Fixed under its own ticket before any number in the phase
was committed, and it is in 0.0.10. An exhaustive harness is a fuzzer that knows when it is finished.

---

## Phase 24 — thirty run, and twenty-four with a falsified clause

The largest phase in the record, and the one where the falsifications are the point. Six registrations held
every clause. Two were falsified on every clause. The other twenty-two split, which is what a clause-level
verdict is for: a row that reads "held" and a row that reads "falsified" would both be summaries of
something more useful than either.

| prediction | verdict | what it means if you depend on this crate |
|---|---|---|
| **P-75** — per-vertex material weights, computed once inside the extraction walk, beat a per-fragment edit-log walk | **held**, all three (M-387) | **Carry the weights on the vertex.** They cost **3.11–3.16%** of extraction and the share is flat across a 4× log; they beat the per-fragment walk by **12.93× to 40.89×**, widening monotonically across `M-50`'s four buckets; and interpolating them across a triangle misclassifies **0.18–0.35%** of the surface with **99.36–99.73%** of that inside one cell of a material boundary. C1 held only because the surface carries **half** the vertices the pre-run estimate assumed — the same estimate had predicted it falsified |
| **P-76** — stochastic triplanar selection trades fetches for temporal accumulation, and a destructible world has already spent its temporal budget | **held** on C1 and C3; C2 held on its own falsifier, **vacuous** as worded (M-388) | **Keep three planes, and every falsifier failed.** The fetch count drops by **18 = 3 × 6** exactly and the fragment cost by **2.1115×**; biplanar gets **0.5000** of the saving exactly in fetches and **0.5429** in wall clock. The recommendation survives all of that because stochastic selection manufactures **39,495** extra history rejections in a scene with no dig in it at all, against the dig's own **8,483**. C2 as written is vacuous: no instrument in the file can price an image error against a millisecond |
| **P-81** — golden-section search against the resident field beats the shipped `TriMesh` collider on query cost, deterministically | **held** on C1 and C2; C3's `TriMesh` half held, its field half **vacuous** (M-393) | **8.067×** over `TriMesh` at **4.766 µs** a query against a 3× / 20 µs bar, and **bit-identical over 10⁶ queries** across Zen 3 and an Apple M5 — so this path is not behind the wall every other float-sensitive one is. The `TriMesh` arm reports **93,298 ghost contacts** over one 495-crossing fixture: if a mesh collider ships, internal-edge removal is not optional. But P-85 falsifies C1's registered *share* premise — the collider's 45% is construction to within **3.12e-5** — so this moves a moving body's per-frame collision budget and **not** the 45% |
| **P-98** — the Plantinga–Vegter certificate, fused into the extractor's own eight-corner gather, costs what M-378 derived | **held**, all three (M-410) | **The certificate is now shippable as a default rather than a debug gate.** **0.0553 worst** of extraction at 65³ against M-378's derived **0.0658**, the **same certified set** on 25 of 25 rows with **0** unsound over 2,389 hidden-topology cells, and per-chunk aggregation at **0** extra passes. The catch is in what aggregates: the **count pair, not the fraction** — average the fractions and you are wrong by **12.9 points** on `noise_cavity`. M-378's prose gloss on its own 0.0658 is not a recipe for it |
| **P-99** — the under-resolution witness works once its denominator is the surface rather than the volume | **held**, all three (M-411) | **The reformulation ✗53 asked for lands, and the two findings around it are worth more than the verdict.** The single-root rate ranks `thin_plate` first at 0.441176 and 0.469697, fails to converge on `thin_plate` and converges on the other seven, and the clause ✗53 stated backwards holds in the direction the predicate allows: false negatives **322 / 94 / 8 / 1 / 0** against monotone non-decreasing false positives, over `k` = {2,3,5,9,17}. But C1's premise is off by two — `rank_all_edges` reads **6**, not eighth — and the metric C1 endorses **needs the oracle too**, so the deployable LOD input is `rate_sign_change`, which reproduces the ranking on **all 120 rows** |
| **P-102** — ✗43's withdrawn rate reproduces under P-63's exhaustive machinery | **held**, all three (M-414) | **The shipped extractor is closed on 8,064 configurations, and the reproduction is worth less than it looks.** **0 of 8,064** unclosed post-fix on **five** independent magnitude fixtures, **2 of 8,064** pre-fix at 6³ at `t7c6o0+t9c6o0`, and the ±1 arm reports VOID rather than a pass. But C2's share is `fan_configurations` = **2**, so *"exactly 2"* was reachable by exactly 2 and *"zero elsewhere"* was **arithmetically forced rather than measured** |
| **P-74** — sphere-tracing the resident field gives ambient occlusion under 2 ms, with no seam or halo error and no construction cost | **held** on C1, **split** on C2, **falsified** on C3 (✗57 / M-386) | **The first lighting measurement in this repository: field AO is cheap and it is not free.** `ao_ms_field` is **0.3996 / 0.1363 ms**, still under the bar at **0.9406 / 1.2215 ms** once the 42.49% / 11.16% coverage is corrected out — though the bar is an absolute millisecond and its own verdict flipped between two runs of one binary. SSAO's error is **1.52–19.09×** the field arm's on all six sets, so the comparison holds, while both absolute halves fail: the field arm's own seam and halo error are 0.021313 / 0.034184 and 0.007323 / 0.028428, not zero. And the construction RTSDF pays is not already paid here — **10.41 / 12.74 ms** of it — because `raw_lipschitz` is **1.4140 / 2.1793** and a sphere trace on the resident field is unsound |
| **P-78** — probe invalidation is edit-proportional, and tracing the field beats tracing the extracted mesh for probe updates | C1's factor **falsified**, C1's world-size falsifier **held**, C2 **falsified inverted**, C3 **held** (✗59 / M-390) | **Probe GI survives the question that would have killed it and loses the one it was expected to win.** The invalidation factor is **3.36–8.38** with `c1_holds` true on **18 of 72** rows — but invalidation does not grow with the world: **1.391×** against a world that grew **8.00×**. C2 is falsified *inverted*: the `parry3d` BVH over the extracted mesh is **1.95–7.27× faster** than sphere-tracing the field. And the topological event does have a lighting signature — a breakthrough dig invalidates strictly more probes on **36 of 36** pairs at **1.25–2.53×**, on two digs matched to **171 samples exactly** |
| **P-79** — shadow-page invalidation is proportional to the brush's projected area, localised invalidation is pixel-identical, and caching beats Epic's uncached choice | **falsified** on C1 and C2; C3 **held** at the registered configuration and **falsified** over the sweep (✗60 / M-391) | **For a world being dug away, page invalidation is a *quantisation* and not an area.** The constant is **20.3718** against a bar of 3. Localised invalidation is not pixel-identical — **4 texels** on one row, with a **66-of-96-arm** geometric leak underneath it that a page-level test would have shipped; `radius + 2√3·cell` closes it at zero leak on 384/384. Caching is worth **3.96×** on `game_dig`'s own brush and loses on **189 of 384** rows, so Epic's uncached choice is right at the wheel's maximum and wrong at the default |
| **P-80** — the fine-to-coarse residual is computable from the field alone, and a normal map baked from it restores what coarsening costs except where the loss is a silhouette | **falsified** on C1 and C2; C3 **held** on `thin_plate`, **falsified** on `gyroid` (✗61 / M-392) | **The map works and the clause missed by two thirds of a degree.** C1 fails on one field of eight — `noise_cavity` at 0.899724 and 0.836800 against 0.95, both halves of its falsifier firing — and holds on the other seven at a worst `residual_p95_cells` of 0.057262 against a 0.1 bar. C2 misses by **0.62°** against a headroom of **0.94°** its own reference left and nobody computed. C3's prediction holds on `thin_plate` at **28.60°** and fails on `gyroid` at **13.22°** — where the *no-map* arm is itself only **14.85°**, so the clause's premise failed rather than the map succeeding: a mean shading-normal angle over the surviving coarse surface cannot see a silhouette |
| **P-82** — a `game_dig` projectile tunnels a thinning wall and continuous collision stops it, at under 25 µs an element | C1 **falsified** by its own falsifier; C2 **held**; C3 **held** as registered and **falsified** off that fixture (✗62 / M-394) | **The discrete path never tunnels, so there is nothing to buy.** The capture window is `2R + t`, not `t`: `game_dig`'s 2-cell body radius makes `2R` **4.0 cells** against a fastest per-frame step of **2.4588 cells**, giving **0 tunnels over 500,000 shots** on all 25 scored thickness×speed combinations — and thinning the wall from 2 cells to a twentieth of a cell changes the zero not at all. Cost is **1.5644 µs/element tested** against a 25 µs bar. The 8-swept-sphere proxy agrees to **0.1629 cells** as registered and misses on a convex rim above aspect 8.0, where the registered fixture is support-exact and cannot tell the two apart |
| **P-83** — the divergence-theorem surface integral gives mass properties to 1e-4, deterministically, at under 2% of extraction | C1 **falsified**; C2 splits — determinism **held**, the asymmetry half **falsified**; C3 **falsified** (✗63 / M-395) | **The whole gap is Marching Cubes, and no integrator could have passed C1.** Off by **44× to 2,275×** at 33³ on all eight fields, while a volume integral over the same triangles agrees with the surface integral to **1.7e-13**. Cross-machine determinism holds on **48 of 48 hashes** across Zen 3 and an Apple M5; the paper's non-symmetric tensor does not reproduce — an exact **`0.000000e0`** on `box_exact` — so the symmetrisation step can be dropped. Cost fails on **20 of 24** rows, and the paper's surface form costs **2.14×** the classical origin-fan it was registered as cheaper than. No swap is recommended |
| **P-84** — a brush clipped against a precomputed convex-cell partition fractures in under 10 ms a fragment, and the partition is a pure function of the edit log | C1 **held** by three to four orders; C2 **falsified** as registered and **held** under one named repair; C3 **falsified** (✗64 / M-396) | **Runtime convex decomposition is not the problem; ordering is.** The fracture costs **10,464–12,726 cycles per fragment** against `M-116`'s 240.7–271.8 ms. The partition is a pure function of the edit log **only** if a crossing vertex is solved from its stable-id plane triple *and* coincident planes are collapsed — **6 / 6 / 2 / 1** distinct partitions over all 40,320 orderings, and one of those four arms is the shippable one. Piece count grows **6.412534 and 6.019270** at 16 cells a chunk against a 4× bar and **1.360286 and 1.381061** at 32, so the clause is resolution-dependent in a way the registration did not allow for |
| **P-85** — one stage of the collider build is over half its cost, and it is the triangle copy rather than the BVH the docs blame | **held** on C1; **falsified** *twice over* on C2; **falsified** on C3 (✗65 / M-397) | **The collider's 45% is this crate's own weld and readiness gate.** One stage is **80.3–81.8%** of the build, and it is neither candidate: the copy is **0.035–0.050%**, the BVH does not dominate either, and all of `parry3d` is **18.2–19.7%**. And the cost *is* predictable from triangle count — `gyroid` costs **1.017–1.047×** `sphere`'s against a 1.5× bar, at matched triangle *and* vertex count — which is the best available outcome for a scheduler and the opposite of what was registered |
| **P-86** — slivers are what a capsule controller catches on | **falsified** on C1; **falsified** by the letter of its column on C2; **held** weakly on C3 (✗66 / M-398) | **The recorded `degenerate_triangles` metric has no gameplay consequence.** **3.01%** of controller stops land on the bottom aspect-ratio decile against a 20% bar — and that bar was **15.29× chance** rather than 2×, because a capsule hits a sliver at **2.30×** its area share, which is real enrichment and still six times too little. The seam column reads **1.134**, and **1.055 ± 0.101** over episodes, which does not reopen `M-133`. C3 holds at 0 of 89 against a null expectation of **1.052** — a held clause with almost nothing behind it |
| **P-87** — an octree navigation graph built on the sign bitmap repairs a dig locally, edit-proportionally, and merges convexly by 5× | **held** on C1 and C2, **falsified** on C3 (✗67 / M-399) | **Navigation on this crate's own sign bitmap is affordable and it never touches a triangle.** P-72's eleven-edit dig repairs in **0.1433–0.4366 ms** for the worst single edit and **0.8793–2.5396 ms** for all eleven, over a repair set **1.50–1.83×** the crate's own dirty-cell count, growing **1.6355×** while the world grows **512×**. The convexity-preserving merge gives **2.5492–3.8597×**, not the registered 5×, with `gyroid` worst as predicted — because the order of magnitude is in the **octree**, not in the merge |
| **P-88** — an octree of free cells gives CALIBRE's clearance without a medial axis, locally, and never overstates it | C1 **falsified**; C2 **held** for the cheap reading and **falsified** for the tight one; C3 **held** on the analytic fixtures and **falsified** on the reference fields (✗68 / M-400) | **The octree does give a clearance number without a medial axis, and there is no version of it that is both accurate enough and local.** Within one voxel on **0.688892** and **0.831008** of samples against a 0.90 bar. λ-membership flips stay at 0.8482 of the changed samples, well inside the 4× bar, but the flip set is a subset of P-87's repair set only for the `region` reading — `ball` puts 149–682 flips outside it at every λ ≥ 1.0. And the one-sided clause dies where it matters: **62 / 649 overstatements** on the reference fields, 0 on the 480,000 analytic samples. A gate that lies in the permissive direction is worse than no gate, so CALIBRE still needs ρ |
| **P-89** — 1³ is worse than 2³, and M-377's two-term model extrapolates there | **held** on C1, **falsified** on C2 (✗69 / M-401) | **M-377's 4³ optimum is now interior with measured neighbours on both sides and nothing left untested below it.** 1³ loses by **4.36×** on `gyroid` and **1.76×** on `fbm_terrain`. The two-term model that fitted a stencil of exactly 6 on all twelve of M-377's arms does *not* extrapolate to c = 1 — **267%** and **46%** error — because at one cell a chunk the cost stops being field evaluations |
| **P-90** — Lipschitz brush pruning still pays at a 4³ brick, the per-brick list is under 8 bytes, and the two culls compose | **held** on C1, **falsified** on C2 and C3 (✗70 / M-402) | **Pruning survives the small brick; the bookkeeping does not, and the second cull is the first one twice.** C1 reads **6.3916×** and **10.4802×** on the registered 46–60 bucket, but only at the brick: P-39's chunk-level cull alone reads **1.5799×** and **1.8825×**, under its own 2× bar and under its own clock-free evaluation ratios of **1.5956** and **1.8993**. The list costs **9.398804** and **10.996338 bytes per brick**, lost to the index and not the payload. C3 is the registered null — `additional_removed` **0 on 8 of 8** — because the brick enclosure is *nested* inside the chunk's and the two culls are one test at two scales |
| **P-91** — geomorph works on terrain and fails on caves, and dither is the expensive option | **falsified** on C1's terrain half and C2's cost half; **held** on C2's pop half and, by sign only, on C3 (✗71 / M-403) | **Item 3.6's prediction is exactly inverted.** Geomorph FAILS on `fbm_terrain` at **25.0771 px** and the `gyroid` half it was predicted to lose is the one that held, at **4.1900 px** — because Lengyel's containing-cell rule is missing a morph target on **39.69%** of terrain fine vertices against the **< 1%** the clause needs. Dither, registered as the expensive option, costs **3.49× less** temporal history (123,769 against geomorph's best case 431,578) and gets 0.9683 px on `gyroid`. C3 holds by sign at 1.002476 against a fixture ceiling of 1.075770, which is additivity to engineering precision |
| **P-92** — re-extracting a chunk is cheaper than decoding an encoded mesh of it | **falsified** on C1 on every reading; **held** on C2 and C3 (✗72 / M-404) | **You do not choose between extracting and decoding.** Marginal extraction is **3.00×** the 7.3 ns bar at the cheapest field and **58.4×** at a dug chunk, and it misses by 3.41× even against a decode measured on this machine rather than quoted. But the field plus edit log is **77.6–246.2×** smaller than meshopt's encoding of the same triangles, and a replay is byte-identical across an Apple M5 and the Zen 3 at **0 differing bytes over 10,000 edits**. So the title inverts: you already extracted because the chunk was dug, and the answer is **send the log, keep the mesh local, regenerate because the extraction was going to happen anyway** |
| **P-93** — there is a crossover rate below which 64³ wins on total frame cost, because 4³ uploads more vertex data | **falsified** on C1 and C3, **held** on C2 by **509×** (✗73 / M-405) | **4³ wins unconditionally, and granularity is settled rather than a trade.** C1 was denominated in a quantity that no longer exists: the registration's 87.50% is the **sample-grid** upload on a path `GPU-011a` deleted, while the **vertex** upload C1 is actually about is **1.71%** of a 12.5 edits/s frame at the 4³ optimum — reaching the window floor of 0.1 edits/s would need a bus **20.7× slower** than this rig's. The weld does cost more than the upload it saves, by **509×**. And a mixed per-chunk granularity beats neither fixed choice: 0.2376× and 0.1171× against a 1.5× bar |
| **P-94** — an hour of digging is under 2 MB, coaxial edits collapse, and a bespoke coder beats a general-purpose one | C1 **held** — its constancy half **with no instrument** — C2 and C3 **falsified** (✗74 / M-406) | **The log is small, and it does not collapse.** An hour of continuous digging is **415,038 B** and a saturated hour is **734,400 B** against a 2 MB bar; the per-edit constancy half is 17.000000 bytes on all ten rows and could not have come out otherwise, which is said out loud rather than scored. All **200 of 200** coaxial capsules are individually necessary bit-exactly, so ✗41's survivor arithmetic does not transfer to a real stroke. And the bespoke coder beats `zstd -19` by only **1.42–1.53×** while `zstd` on the *same* residuals lands within **1.39–1.43×** of it — the whole win is the **model** and not the coder, so the format is not worth writing |
| **P-95** — undo beats a re-fold below some stable edit separation, and M-50's curve gives the checkpoint cadence | **falsified** on C1 and C2, **held** on C3 (✗75 / M-407) | **Undo is always a re-fold, so the checkpoint cadence is the only knob.** Re-folding from a checkpoint is cheaper at *every* separation — **16.3×/15.9×/19.0× at `d = 0`** and still **1.30×/1.29×/1.55× at `d = 127`** — and the fitted crossover lies **40.95 / 38.50 / 73.86 edits past the bottom of the log**, which is no crossover at all. M-50's curve predicts a cadence of **193.8** against a fit-free measured **32**, and the disagreement is a *shape* difference, not a scale factor. Undo is bit-exact: 0 hash mismatches over 3,654 undo/redo round trips |
| **P-96** — smooth union's 40,320 orderings differ by less than a tenth of a cell, and the spread is a function of the blend radius | C1 **vacuous** on M-38's own fixture, then **held** and **falsified** on a fixture that can move; C2 and C3 **held** (✗76 / M-408) | **M-38 is retired rather than confirmed.** Its 40,317 distinct results are meshes **2.2e-16 world units** apart, in a branch where **no smoothing happens** at all. On a fixture that can move, the spread holds to **0.0375** of blend radius and is falsified from **0.075**, reaching 0.656099534 at M-38's own `k = 0.15`; it falls with `k`, recovers M-36's single result at `k = 0` exactly, and stays inside the 10k seam shell. So smooth union is gameplay-commutative below **half a cell** of blend radius, and the ordering authority M-38 seemed to demand is needed only above it |
| **P-97** — a hundred-thousand-edit trace replays byte-identically across machines, and if it does not, bisection names the edit | **falsified** on C1 and C2, **held** on C3 (✗77 / M-409) | **The fold is bit-identical and two of the meshes are not.** C1 fails on **2 of 8** fields by exactly **3 bytes** each, and C2 fails with an **empty divergence population** — `first_differing_edit` is `NONE` on all eight and the folded grid differs by 0 bytes over 287,496 each, so there is no differing edit to find. The `f64` fold agrees across x86-64 and AArch64 at **every one of 100,000 prefixes** on all eight fields; the only divergence in 2,470,824 compared bytes is the **sign bit of a NaN** that `unit_gradient` emits in release from a plateau `M-31`'s zero-edit fixtures cannot have. Replay is linear at **1.0055–1.0173×** against a 1.2 bar |
| **P-100** — a 24-tetrahedron barycentric split reaches full octahedral equivariance where the six-tet split is stuck at 12 | **held** on C1 and C3, **falsified** on both halves of C2 (✗78 / M-412) | **M-372's obstruction is real and it is in the *decomposition*.** The barycentric split reaches **48 of 48** with `worst_component_ulp` **0**, bit-exactly, where both six-tet arms read 12 on 32 of 32 in the same run — and it still tiles across chunk seams at **0 open edges**, against a mismatched-diagonal control reading 80–3,765 on the same 80 rows. Removing the obstruction costs **2.4405× the six-tet split and 7.9353× Marching Cubes** against registered bars of 2× and 6×. Invariance is purchasable; the price is a CAD price, not a game's |
| **P-101** — accumulating dual crossings in a relabelling-invariant order takes the duals from 3 of 16 to 12, for free | **falsified** on C1 and C3; **held** on C2 at a stated cost (✗79 / M-413) | **The registered key reorders a reduction `A-016` already made invariant.** `dual_contouring` stays at **3 of 16**, identical to its own baseline, and is *worse than shipped on 9 of them*. The reordering is free exactly as promised — a stated cost of **40 of 48** golden dual hashes, `max_hausdorff_delta` 3.008704396734e-14, topology identical on 320 of 320 — and `manifold_dual_contouring` does not do worse: the two duals agree on **64 of 64** row-pairs. The one reduction that is not invariant is `vec3::length`, and changing only that takes both duals to `pure_permutation_exact` **6 of 6 on 32 of 32 rows** |
| **P-73** — every connectivity-weighted normal is far worse than a central difference, and the `acos` in it is a golden-hash liability | **falsified** on all three (✗56 / M-385) | **The 2005 result does not transfer to trilinear output, and the determinism objection does not exist.** **28 of 60** testable rows are *under* 3×, and the clause is separately vacuous on **6 of 16** cases where the gradient median is exactly zero. The off-surface canary predicts nothing — **0 of 7** fields clear ρ = 0.7, the best ρ anywhere in the file is 0.584622, and `fbm_terrain` is *negative*. And the connectivity route moves **0 of 18** hashes between an Apple M5 and the Zen 3. Nothing in `src/` changed, and that is the finding: do not build the dossier's ~40-line pseudonormal |
| **P-77** — a dig's TAA history rejection is five times steady state and persists, k-DOP clipping recovers half of it, and it is concentrated at the brush | **falsified** on all three (✗58 / M-389) | **A dig's temporal debt is one frame.** The post-edit rate is **2.8997×** steady against a fixture ceiling of **2.9852** computed before the run, and the spike lasts a **single frame** — both halves of C1's own falsifier. k-DOP recovers exactly **0** on all 140 rows, by the paper's own containment guarantee, and its cost half is vacuous because no GPU TAA resolve pass exists in this crate. Concentration is **67.23%** against an 80% bar, while the mechanism the clause was reaching for holds at **1.000000** on the edit-attributable denominator: under locomotion there is no history left to owe |

**Seven of the thirty were decided by arithmetic that was available before the harness existed**, which is
Phase 23's central finding turning up again in the phase that was supposed to have absorbed it. P-93's C1 was
denominated in an 87.50% belonging to a path `GPU-011a` had already deleted. P-92's own gloss on C1 was
arithmetically unreachable before the run. P-83's C1 could not have been passed by any integrator, because
the gap it measures is Marching Cubes rather than the quadrature. P-82's capture window is `2R + t` and the
clause was written against `t`. P-86's 20% bar was **15.29× chance** rather than the 2× it reads as. P-94's
constancy half is 17.000000 bytes an edit and could not have come out otherwise. And P-102's C2 was
*"exactly 2"* against a population of exactly 2. None of those needed a machine; all of them needed the
number to be written down where a reader could check it, which is the only thing this page has ever claimed
for pre-registration.

**The vacuity machinery earned its keep, and it is the reason six of these rows are readable at all.** C1 of
P-96 is vacuous on `M-38`'s own fixture — 40,317 distinct results in a branch where no smoothing happens —
and the row is only a result because a second fixture that *can* move was run beside it. P-102's ±1 arm
reports VOID rather than a pass. P-73 is vacuous on 6 of 16 cases where the denominator is exactly zero,
P-81's C3 on a field arm with nothing to count, P-77's on a GPU pass this crate does not have, and P-95's on
an `N` that does not exist to be stable. A clause that cannot fire is not a clause that held.

---


## Phase 25 — twenty run, seventeen with a falsified clause, and one mechanism that won and was refused anyway

Twenty registrations drawn from a bit-packing and broadword survey, and **not one of them was allowed to
change `crates/isomesh/src/`**. That was the phase's condition rather than its outcome: every arm is a
bench-local reimplementation, so no golden hash moved and no consumer-visible behaviour changed. A
mechanism that earned landing was to get its own ticket afterwards, because a landing not registered in
advance is `V-45`'s failure mode.

**Sixty clauses: 34 held, 24 falsified, 2 vacuous.** Two rows held every clause (`M-419`, `M-430`), one
held every clause that could fire (`M-432`), and the other seventeen carry at least one falsification.

| prediction | verdict | what it means if you depend on this crate |
|---|---|---|
| **P-106** — the twelve cut-edge flags of a cell are a fixed boolean circuit over its eight sign bits, so sixty-four cells' worth fit in one word each | **held**, all three (M-419) | **The phase's one clean mechanism win, and it still should not be built yet.** The circuit is **12** word operations and **24** including the plane build, against a bar of 24 — derived by executing the same source over a *counting* word type rather than by reading the code — and it agrees with the crate's own table on **256 of 256** patterns where a deliberately broken circuit misses **192**. It costs **0.1204** of the byte table's instructions per cell, a **7.7–8.3×** cut on classify. But the 256 cut masks take only **128 distinct values** and every complement pair agrees, so **the twelve flags carry strictly less information than the case index and cannot feed `CASES`** — and the best extraction-level ceiling the saving buys is **2.4752×** on one field and **1.0231×** on another |
| **P-117** — nothing on the hot path is FMA-contraction sensitive, and if something is, `M-31`'s cross-machine golden hashes would already have moved | **held**, all three (M-430) | **The risk is real, named, and unrealised — which is the only outcome that lets you stop worrying about it.** **Eleven** contraction-sensitive expressions under `crates/isomesh/src/`, with measured divergence under deliberate fusion from **2 ULP** (`vec3::dot`) to two **sign-crossing** divergences where an exact `0.0` becomes `-4.44e-16` (`marching_cubes::corner_position`, `dual::place_vertices`). And **0** of **216** golden hashes differ across `x86_64-unknown-linux-gnu` and `aarch64-apple-darwin`, with `mul_add_call_sites` **0** — `rustc` does not contract `a*b + c` by default and this crate never asks it to. So a future `+fma` build, a `-ffast-math` equivalent, or one hand-written `mul_add` would move golden hashes at eleven known sites. That is a build-flag constraint, not a code change |
| **P-119** — double-buffered welding removes `P-9`'s chunk-order spread, under 1.25× cost, with the mesh invariant under every permutation | C1 **vacuous** on the registered column; C2 and C3 **held** (M-432) | **The mechanism works and the registered perturbation had nothing to move.** `baseline_spread` is **0 on 16 of 16** rows for chunk order, so C1 as worded cannot fire — and the census that explains it did not exist until this phase, because the harness declared it and populated none of it. `buckets_of_three_or_more_spanning_chunks`, the only bucket kind whose member count chunk append order can change, reads **0 on 16 of 16**. Meanwhile on `P-9`'s own within-bucket protocol the single-buffer arm reproduces the published spread on all 16 rows and double-buffering takes it to **0**. Cost is **1.0888–1.1059×** in instructions against a 1.25 bar, and `peak_bytes_ratio` is **2.86–4.28×** — a real price for determinism that no clause bounded |
| **P-114** — one bit per 64 cells gives a second level at a sixty-fourth the size, and two levels skip 4,096 cells per test | **held** on C1 and C3; **falsified** on C2 on 21 of 32 rows (✗90 / M-427) | **This is the `✗51` shape: it wins its microbenchmark by 2.57× and earns no ticket.** C1 holds on the sparsest field as measured in all three units — `thin_plate` at **1.612×** (65³) and **2.568×** (129³) against a 1.5× bar, and **2.129×**/**3.318×** in instructions, so the instruction form is the *larger* margin. Then the share closes it: `traversal_share_of_extraction` is **0.000088–0.000994**, so the extraction-level ceiling peaks at **1.00039×** and dips *below* 1. The registration had inherited a 5.5× stage figure and asserted the stage mattered; the harness measured the share instead. And the win does not generalise — `ratio` is **below 1** on all eight `gyroid` and `noise_cavity` rows, at 4–11% active fraction |
| **P-121** — decompose extraction into float work and integer work, and a bit-work total under 15% closes three rows before anyone writes a harness | **held** on C1 and C2; **falsified** on C3 on 18 of 72 rows (✗95 / M-434) | **The phase's gate, and the number three other rows are denominated in.** Worst `residual_share` **0.005195** against a 0.05 bar, and `max_reference_integer_share_at_65` **0.724533** against 0.15 — so `R-103`, `R-104` and `R-106` were live rather than Amdahl-dead. C3's falsification is a **qualification on every share this phase quotes**: on 2 of 8 groups the surface-free arm does not separate `classify` from `emit`, so `integer_share` is a share of their *sum* and no row may quote `classify` alone as isolated. Several stage columns go negative — prefix differences inverted by clock movement, bounded to half a percent by C1's own bar |
| **P-103** — FastLanes' `Scalar_T64` beats naive scalar by up to 8×, and the one quantity already boolean here is a cell's eight corner signs | **falsified** on C1 as a whole-fixture gate and on C2; **held** on C3 (✗80 / M-416) | **The mechanism did exactly what the paper says and it bought nothing.** Sign tests per cell fall from **8.0000** to **1.0476–1.0967** — a 7.3–7.6× cut in precisely the quantity `Scalar_T64` targets — and `instruction_ratio` still spans **0.579154–1.08354** against a bar of 0.5. The reason is in the comparand: the bit-sliced arm's instructions per cell span only **26.69–29.96** across all 32 rows while the byte arm's span **24.77–50.81**, so every apparent win is the byte path paying an `f64` penalty, and at 65³ `f32` the bit-sliced arm **loses outright** |
| **P-104** — the active-cell bitmap is FastLanes' horizontal layout, and `M-287`'s stride tax says so | C1 **vacuous**; C2 **falsified**; C3 **held** (✗81 / M-417) | **The layout the clause is about is unreachable through the public API, and a coincidence was retracted.** `row_stride(size) = size[0] | 1` since `A-024`, and an odd row greater than one can never divide a power of two, so no caller can reach the taxed layout. The first clean run appeared to reproduce `M-287`'s published 3.37 at **3.3714**; the re-run reads **8.1052** on the same measurement of the same binary, while `c2_instruction_ratio_fold` reads **1.8054** on both. Of the 65 columns the two datasets share, **30 are byte-identical and 35 moved** — instruction columns by at most **0.006%**, cycle columns by up to **221.8%** |
| **P-105** — nothing in the crate popcounts the bitmap, so Harley–Seal's carry-save reduction is worth pricing | **falsified** on C1; **held** on C2; C3 **vacuous** on 28 of 32 rows and **falsified** on the rest (✗82 / M-418) | **LLVM already took the saving Harley–Seal is selling.** The naive fold autovectorises to **35/4 = 8.75** instructions per word exactly; the carry-save chain's nine live accumulators keep it scalar at **145/16 = 9.0625**. So the mechanism issues *more* instructions per word, and the null was measured in its **best case** — `objdump` finds zero `popcnt` in the binary, because the repository sets no `target-cpu`. C3 is vacuous because a four-point monotone trend spans less than one row's own five-window reproducibility, up to 25× less. And `count_ones_in_dual_rs` is **0**, so the row prices a stage the crate would have to *add* |
| **P-107** — Flying Edges spends a prefix-sum pass to turn per-row counts into output offsets, and a two-level rank directory makes it a lookup | **falsified** on C1; **held** on C2 and C3 (✗83 / M-420) | **Amdahl-dead exactly the way `✗51` was, and the registration named it in advance.** `offset_share` is **0.000328–0.009714** against a 0.05 bar, so the offset stage is at most **0.97%** of extraction and the whole extraction-level ceiling is **1.0098×** even if the stage were free. The structure is correct and cheap — `overhead_fraction` **0.032227** against 0.0351, a rank query folding **at most 7 words at every resolution** — and it buys **random access, not throughput**: the shipped sequential walk costs a fraction of an instruction per query against the rank query's ~127 |
| **P-108** — Vigna's broadword select needs neither `PEXT`/`PDEP` nor a table, and the walk it replaces is a registered expected null | **falsified** on C1 as registered; **held** on C2 and C3 (✗84 / M-421) | **The null arrived, and the share was the real answer.** Broadword select costs **+48.9 to +61.0** instructions per set bit over `dual.rs:489-497`, an instruction ratio of **0.2085–0.3165**. But `walk_share_of_extraction` is **0.000103–0.001462** — between one part in ten thousand and one and a half parts in a thousand — so deleting the walk entirely buys at most **1.00146×**. The scored arm makes **zero** `count_ones` calls, and the unscored literal variant that does is priced: a hardware `POPCNT` narrows the gap from about 4× to about 3× and cannot close it |
| **P-109** — Elias–Fano encodes a monotone sequence in about `2 + ⌈log₂(u/n)⌉` bits with O(1) access, and within-chunk edge indices are monotone | **falsified** on C1 and C2; **held** on C3 (✗85 / M-422) | **Right structure, wrong sparsity, and the cost model forecloses it arithmetically.** A 4-bit budget demands `u/n ≤ 8`; `universe_per_crossing` measures **16.5304–786.3330**, so the fixture's *densest* row is already 2.07× past the ceiling. Measured **7.2959–12.6711** bits per crossing, with `low_bits_width` equal to `⌈log₂(u/n)⌉` on 24 of 24 and the total within one bit of Vigna's own model. Access costs **36.8–158.7×** direct addressing in instructions. The encoding still beats the *shipped* structure's space by up to 563×, which is the number a future reader should quote rather than the bar |
| **P-110** — Pibiri & Kanda's contribution is rank/select that supports `flip(i)` while staying fast, which is the shape a structure written during a sweep needs | **falsified** on C1; **held** on C2 and C3 (✗86 / M-423) | **The order-independence proof is the durable half.** Over **128** orderings of 4,096 flips, `distinct_final_states` is **1** on every row while the order-sensitive control reports **128** — so `V-45`'s objection does not arise, and the control is what makes the 1 admissible. C1 fails in every unit: the rank-only conjunct is a **2.1–2.4%** near miss whose sign the unit decides, and the flip-plus-rank conjunct decides the verdict at 2.24–2.38×. The cost is located structurally — static rank is flat across all 24 rows while mutable rank climbs monotonically with the top-level scan |
| **P-111** — simdjson's 256-entry byte-to-permutation table, reduced to scalar, is a 2 KiB const this crate can actually have | **falsified** on C1 in all four units and on C3 on 14 of 16 rows; **held** on C2 (✗87 / M-424) | **Decidable from a bit census with no clock at all.** `stores_per_set_bit_walk` is **1.00000** on 16 of 16; `stores_per_set_bit_table` is **2.00000–22.17637**, exactly `64·words_nonzero/set_bits`. A mechanism issuing 2–22 stores per useful store cannot beat one issuing one. The loss is *located* rather than observed: a byte-skipping variant drops stores to near 1 and still loses on all 16 rows, so the cost is the lookup. And C3's L1 claim survives the instrument on only **11 of 16** rows — one row where the table genuinely reads fewer misses per set bit |
| **P-112** — Billeter's count-scan-scatter has a prefix sum in the middle, and `P-107`'s O(1) rank replaces it | **falsified** on C1 and C2; **held** on C3 (✗88 / M-425) | **Closed by its own registered falsifier, and it left a method rule behind.** The 70% scatter falsifier fires on **7 of 16** rows; `scan_share` is **0.036580–0.158941** against 0.20, whose scan-free ceiling is **1.1890×** — under the 1.25 C1 asked for, so C2 makes C1 arithmetically unreachable. The rank arm is **2.2–8.6× slower**. And the whole three-phase compaction is **0.0162%–0.4825%** of an extraction, so a free scan buys **1.00063×**. The `R-085` residual bar fired twice before the dataset landed and the diagnosis is the transferable part: the decomposition is exact **to the last retired instruction** while its nanosecond form reads 4–8% with a field-dependent sign |
| **P-113** — Roaring switches container type at 6.25% density and says a bitmap is the wrong structure below 0.1% | **falsified** on C1 and C2; **held** on C3 (✗89 / M-426) | **The clause quantified over the wrong axis, and a span control passed anyway.** Density is a function of chunk size, not of field: **0 of 4** granularities have the eight fields straddling 6.25%, and within any one granularity they never cross it. The `density_span_over_fixture` control cleared its one-decade bar and the span is generated *entirely* by chunk size — a span control proves the range is wide, not wide along the axis the clause is about. C2's two failures were decidable from `size_of` alone. And the walk timer is **saturated on the array arm on 32 of 32 rows**, so `bitmap < array` was unobservable and the verdict is carried by instructions |
| **P-115** — Tree-Encoded Bitmaps are the only RLE-family scheme that preserves the random access rank needs | **falsified** on C1; **held** on C2 and C3 (✗91 / M-428) | **One measurement retires TEB, WAH and EWAH for this crate, and the registration pre-committed to that reading.** `access_ratio` is **57.77–115.77×** the flat bitmap's against a bar of **2.0** — 29× to 58× past it — on both rank granularities and in instructions. Space compresses **2.37–563.75×**, and the pairing is the finding: `thin_plate` 129³ is **563.75× smaller and 63.9× slower to read**, against a flat bitmap whose entire size is 256 KiB. Subtracting the measured instrument floor makes the ratios *worse*, and the implementation is not a strawman — `teb_bytes` is below the fully-pruned variant on 6 of 16 rows |
| **P-116** — GRAPHGEN synthesises the optimal decision tree over a decision table, compresses it to a DRAG, and emits code | **held** on C1 and C3; **falsified** on C2 on 8 of 8 rows (✗92 / M-429) | **The synthesis is sound and inert, and that is a fact about the Marching Cubes case table.** `average_path_len` is **8.000000** — the optimal tree is the *complete* eight-level tree — and `min` equals `max` over **6,561** DP states, so **every** ordering of the eight conditions is optimal. The DRAG then merges exactly **one** node of **511**, and it is a leaf merge. With 256 patterns, **0** don't-cares and **255** distinct actions there is no redundancy to exploit. The emitted 1,311 lines compile under `rustc`, Kani discharges **62** checks in 2.04 s with a sabotage that fails as it must, and the generated form costs **1.65%–3.83% more** |
| **P-118** — Neal's superaccumulator is exactly order-independent, so it is a second route to `P-101`'s question that needs no invariant key | **falsified** on C1 and C2; **held** on C3 (✗93 / M-431) | **`M-177` is closed: it is not an accumulation defect.** With *every* reduction on the dual vertex path exact — twelve superaccumulated accumulators plus `vec3::length`'s dot — `pure_sign_flip_exact` still reads **6 of 32**, while the same arm takes `pure_permutation_exact` from 8 of 32 to **32 of 32**. So the permutation half is entirely `vec3::length` and the sign-flip half is structural, as `M-177` said. C1 was out of reach before the run and a column says so. C2 fails at **7.93×** instructions against a 2× bar, with a named mechanism: twelve reductions per solve, each paying a whole 67-chunk read-out over about four terms. Incidentally measured for the first time: `sum_equivariant` costs **4.64×** a naive ordered sum on the dual solve |
| **P-120** — Wu, Otoo & Shoshani replace pointer-based equivalence trees with a flat array and report 5× to 100× | **falsified** on C1 on 4 of 32 rows; **held** on C2 and C3 (✗94 / M-433) | **The `✗26` objection does not reproduce, and that is a claim about *this* union-find.** `order_independent`, `labels_consecutive` and `partition_identical` are all true on **32 of 32** across three scan orders — but every arm here is per-chunk, bench-local, rebuilt from scratch and never asked to delete, and `connectivity.rs:29-46`'s refusal is about an incrementally maintained structure asked to handle `fill`. C1 clears 2.0× on 28 of 32 rows and tops out at **2.3759×**; the four failures are all `fbm_terrain`, the only field whose air fraction is near half, so its *comparand* gets cheap. On instructions the bar clears on only 14 of 32 |
| **P-122** — Stream VByte keeps a dense control stream and a variable-length data stream separate so the decoder never branches on payload | **falsified** on C1; **held** on C2 and C3 (✗96 / M-435) | **The split never reduces mispredictions and is cheaper anyway.** `branch_miss_ratio` spans **1.001214–1.983313** on 16 of 16 rows — it *doubles* them at worst — while instructions per cell fall by **2.9%–14.6%**, and the harness closes the obvious alternative explanation by showing the split makes *more* corner gathers, not fewer. There was very little to remove: `branch_misses_per_cell_single` is **0.008008–0.250699**, under one misprediction per 125 cells at the sparsest. And **one of sixteen rows' verdict depends on which order the two arms were measured in**, so the harness records both orderings and refuses to quote that row as firm |

**Eleven of the twenty were decided by arithmetic available before the harness existed**, which is Phase
24's central finding turning up a third time. `P-107`'s and `P-108`'s and `P-112`'s clauses are share bars
against shares of under 1%. `P-109`'s 4-bit budget needs `u/n ≤ 8` against a measured 16.5–786.3.
`P-111`'s stores-per-set-bit ratio is `64·words_nonzero/set_bits` and needs no clock. `P-113`'s byte
crossover is `1/(8·index_bytes)` and its two failures fall exactly where `size_of` puts them. `P-114`'s
one-level rows cannot reach 4,096 by construction. `P-116`'s DRAG has nothing to compress because the
table has no don't-cares. `P-118`'s C1 was closed by `p-101.csv` before `experiment_p118.rs` was written,
and the file records that in a column. `P-103`'s ceiling and `P-121`'s own bar complete the eleven.

**Five refute a figure this repository or a cited paper had already published.** `M-287`'s 3.37× stride tax
is not reproducible as a cycle ratio — the same measurement of the same binary reads 3.3714 and 8.1052.
`M-337`'s 5.5× stage figure, which `P-114`'s registration inherited, is a share of under 0.1% of
extraction. `M-306`'s rejected-brush shares were read as active-cell densities and are high by 12.5×.
Roaring's 6.25% is the two-byte case only. And `M-177`'s obstruction, which three separate mechanisms have
now failed to move, is confirmed structural rather than arithmetic.

**Two harnesses aborted on their own vacuity controls before their datasets landed, and both aborts were
worth more than the clause they interrupted.** `experiment_p112.rs`'s residual bar fired twice and the
diagnosis retired a wall-clock control in favour of an instruction one, with a mutant control attached to
prove the resulting zero could have been non-zero. `experiment_p116.rs`'s emitted-text cross-check was
comparing the source against its own reverse, and the three-conjunct instrumentation located the defect in
exactly one of them. A third defect — a committed CSV column recording a decimal digit count rather than a
stride — was found by writing this page's row for `P-104` against the dataset rather than against the run.

---

## Phase 26 — four experiments the backlog had left unregistered, and two of them close a question

Not a sweep. These four were already in `BACKLOG.md`, unblocked, and had **never been registered**:
`R-027a` split out of `R-027` on `V-45` in Phase 17, `R-052` and `R-053` left as the successors to
`R-050` and `R-048` when Phase 20 closed, and `R-072` named in `O-12`'s own text and then filed by
nobody. That last one is the reason the phase exists at all — an open question was invisible to
`scripts/backlog_gate.sh` clause 6 for want of a `P-` id to be missing, which is exactly the drift the
clause exists to surface and could not.

**Three of the four are corrections to an experiment that already ran**, and each names in its own
registration the error it is not repeating. **One landed a source change**, registered in advance.

| prediction | verdict | what it means if you depend on this crate |
|---|---|---|
| **P-126** — `O-12`'s remaining half: is the *dual* vertex link a single cycle, over all 2²⁷ sign patterns of a 3×3×3 block? | **held**, all three (M-439) | **`O-12` is closed by exhaustion, and the closure is a standing gate rather than a claim.** `worst_link_components` reads **1** on every non-control arm over **134,217,728** patterns and **5,067,767,808** dual vertices, while all three controls read **2** with 99.9M–221.8M link-defective vertices. C1 is additionally *proved*: two lemmas checked exhaustively at start-up — a cycle's consecutive edges always share a cell face over 256 cases × 64 masks, and with the joined bit agreed the two cells sharing a face induce the identical pairing over 49,152 combinations — so the sweep's job is that no reachable pattern escapes the proof. `control/open_block` reproduces `P-63`'s own vacuity as a *measurement* at the right block size: 5.07 G dual vertices of which only **860,880,896** carry a complete link |
| **P-124** — port Finken et al.'s monotone-edge condition to the **ambient** complex, where `✗36` proved the mesh-edge reading is saturated | **held**, all three (M-437) | **The detector has signal, and it is almost entirely in the edges a chord predicate cannot reach.** At 65³ over 262,144 cells, `diagonal_only_failures` is **17,862 of 17,862** on `sphere`, **10,563 of 10,567** on `box_exact` and **12,669 of 13,693** on `thin_plate` — cells rejected by a tetrahedral diagonal with every axis edge monotone. The mesh-edge control runs the same predicate in the same run and reproduces `✗36`'s saturation, which is what makes the ambient count a property of the *complex* rather than of the tolerance. The tolerance was fixed at registration — scaled by `max(|f(a)|, |f(b)|)` — because a tolerance vanishing at the endpoints is precisely what `✗36` caught |
| **P-123** — decompose `M-318`'s 45× buffer-churn gap; if the order-only term dominates, a canonical reorder at emission recovers it | **held** on C1, **falsified** on C2, C3 held on 3 rows and **vacuous** on 3 (✗97 / M-436) | **The gap is not ordering, and `R-027` is closed for good.** The order-only term is **0.0000%–0.137%** of the churn; the predecessor-shift term is **97.9%–99.4%**. `sphere` 129³ reproduces `M-318` to the digit — `churn_total` **15,706**, `churn_geometric` **346** — and spends 15,370 of it on slots shifted because an earlier cell emitted a different *number* of vertices. That is Acar's sequential-counter instability, and no reorder at emission touches it: the only remaining shape is the persistent edge→slot map `V-45` refused, which makes `extract_into` stop being a pure function of its inputs. The vacuity control fired both ways — the canonical arm drove the order-only term to **exactly 0** while leaving geometry unchanged, and a permuted arm drove the same detector non-zero |
| **P-125** — ship the pinch predicate as a `validate` report, so a caller can ask whether the `=`-corner repair is safe on **its** data | **held** on C1 and C2, **falsified** on C3 on 2 of 6 rows (✗98 / M-438) | **The clearance is sound and the count is not, and the difference matters.** C1 reproduces `M-352` exactly — **516 of 17,201** collapse groups on `bonsai`, **0 of 50** on `fuel` — from the baseline mesh alone with no repair applied, and C2 gives one identical census over every face permutation, which is `✗26`'s objection asked in advance. C3 fails because a pinch count is an **upper bound**: 516 pinches weld **129** components, since several pinches can weld the same pair and a pinch inside one component welds nothing. So *zero pinches is a genuine clearance* — confirmed on all three zero rows — and a non-zero count bounds what a weld could reach rather than predicting it |

**The four registrations were written the same day and three of them contained an error the harness
author found and reported rather than amended.** That is the machinery catching the person operating it,
for the third phase running, and all three were available before any harness existed:

- **`P-124`'s C2 is internally inconsistent.** It says the non-monotone population *"is O(n) and HALVES
  PER REFINEMENT, so the count … must fall by a factor in [1.7, 2.3] per doubling"*. An O(n)
  population's **count grows**; it is the **density** that halves, and `✗36`'s own table is a density
  (178.1 / 86.1 / 42.3 / 21.0 per 1k). The literal reading was arithmetically dead.
- **`P-126` wrote 2²⁷ for a 64-corner block**, which is 2⁶⁴. The harness resolved it with a
  period-3-per-axis identification that makes both registered numbers true, recorded the identification
  as a column, and stated the cost: the sweep is exhaustive over 3-periodic fields and visits each
  configuration 128 times.
- **`P-126`'s C2 names `FaceAmbiguity::Connected`, which does not exist** — `ambiguity.rs` ships only
  `Separate` and `AsymptoticDecider` — so the clause as literally written cannot fire, and the harness
  says why: an always-joined rule is a function of the shared face's four signs, both cells agree, and a
  link provably cannot split when they agree.

**And one harness rejected its own first implementation on a measurement.** `P-125`'s report groups
coincident vertices, and the obvious grouping is equality on the lattice cell key — a genuine
equivalence relation, and `MeshReport::weld_buckets`' own notion. Measured against the shipped welder on
`sphere` 25³, where `M-48` recorded `Welder::weld` removing 48 vertices and 96 triangles at this very
epsilon, it finds **30 of the 48 merges and 60 of the 96 folding faces**, because a coincidence class
straddling a cell face splits between two buckets. For a report whose entire job is issuing clearances,
under-reporting is a **false clearance**. Rewritten as ε-connected components: **48 and 96, exactly the
weld**, and erring on the safe side because the welder's classes refine the closure.

**Three of the four experiments contain no clock at all** — every column is an integer count over an
enumerated population — so `M-280` has nothing to move and their clean re-runs could not have shifted a
digit. `P-126`'s `wall_seconds` sums to about **1,015 s**, inside the ticket's own 4–45 minute estimate,
and no clause reads it.

---

## The 4.26×, and not one triangle changed

Two changes, three days apart, on the same function.

**A-023 made the loop axis a constant.** `DualMesher::emit_quads` took `axis`, `u` and `v` as runtime
values, so `p[axis] = a` was a dynamically indexed store. The coordinate array could not live in
registers, and every iteration wrote three coordinates to the stack that `linearize` read straight
back — a store-to-load forwarding chain. Splitting it into three `const`-generic monomorphisations
left the same three passes in the same order over the same bounds:

| `emit_quads`, per iteration | before | after |
|---|---|---|
| cycles | 43.26 | **3.33** |
| instructions | 31.24 | **12.02** |
| share of the mesher's cycles at 193³ | 82% | **27–32%** |

**A-024 forced the row length odd.** `values` was laid out by the caller's shape, so at 128 samples
per axis the row stride was **512 bytes** and the plane stride **exactly 64 KiB** — a cache-set
aliasing period on this machine, twice over. The remedy is `size[0] | 1`, unconditional and
idempotent. Unconditional because a pad applied only when the stride looks bad is a second layout
reachable from one call; idempotent because a *fixed* pad of one would map every `size[0] = 2ᵏ − 1`
onto the very stride it is avoiding.

Together, on a sphere at 256³, `f32`, one thread, one binary, one run:

| | before | after |
|---|---|---|
| Surface Nets | 693.8 ms | **162.7 ms** — 4.26× |
| Surface Nets IPC | 1.20 | **4.09** |
| `SN / MC` | 5.43× | **1.26×** |
| per-sample cost, 16³ → 256³ | 29.63 → 41.35 ns (+40%) | **8.71 → 9.70 ns (+11%)** |
| cache misses per sample | 3.72 | **1.56** |

**Not one triangle changed.** The golden hashes pass untouched — an optimisation that changes the mesh
is a bug in the optimisation.

And it cost this repository one of its own beliefs. ✗14 read *"Surface Nets never wins on Zen 3, at
any resolution: 2.46× behind even at 16³."* Surface Nets is now **faster** than Marching Cubes at 16³,
24³ and 32³:

| n | 16 | 24 | 32 | 48 | 64 | 96 | 128 | 192 | 256 |
|---|---|---|---|---|---|---|---|---|---|
| Marching Cubes ns/sample | 9.37 | 8.87 | 8.29 | 8.03 | 7.74 | 7.41 | 7.48 | 7.62 | 7.67 |
| Surface Nets ns/sample | **8.71** | **8.59** | **8.05** | 8.24 | 8.42 | 8.53 | 8.75 | 9.53 | 9.70 |
| ratio | **0.93** | **0.97** | **0.97** | 1.03 | 1.09 | 1.15 | 1.17 | 1.25 | 1.26 |

The reason it had been true was one missing `const`.

---

## One bit of a row length, at the chunk size everybody uses

128³ is the canonical chunk size in voxel engines, and it was the worst possible one.

The diagnosis came before any code changed, by letting the *caller* arrange the pad — adding one
sample moves the stride while changing the work by under 1%:

| shape | plane bytes | cycles/sample | vs 127³ and 129³ |
|---|---|---|---|
| 127³ | 64,516 | 33.10 | — |
| **128³** | **65,536 = 2¹⁶** | **108.51** | **3.37×** |
| 129³ | 66,564 | 31.39 | — |
| 129×128×128 — pad `x` | 66,048 | 31.48 | 0.98× |
| 128×129×128 — pad `y` | 66,048 | 36.45 | 1.13× |
| **128×128×129 — pad `z`, the control** | **65,536 = 2¹⁶** | **107.89** | **3.35×** |

**The `z` row is what makes this a measurement rather than a story.** It adds the same 0.8% of work and
touches neither stride, and it keeps the entire penalty. Without it, "adding a sample helped" would be
indistinguishable from "the fixture got slightly different."

Measured on a field with no surface in it at all — a sphere entirely inside the domain — so this is the
scaffolding, not the geometry.

---

## When a measurement is impossible, suspect the instrument

This is the most useful entry in the file, and it is an entry about being wrong.

R-006 and R-008 both compared vertex normals against a *reference* gradient computed from an exact
wedge SDF written for the experiment. Both came back falsified. P-16 missed its threshold by four
times over. Two independent hypotheses about sharp features, both dead, on the same afternoon.

The reference gradient was returning noise.

Its exterior branch normalised `away = q − dir·t` guarded only by `e > 0.0`. For a point lying *on* the
ray, that difference is a catastrophic-cancellation residue of order `ε·|q|` — a random direction with
a plausible magnitude. About half of all Marching Cubes vertices land epsilon-outside the wedge, so
about half were being compared against a random unit vector.

| | as reported | corrected |
|---|---|---|
| past-90° vertices, Marching Cubes | 6,959 | **472** |
| past-90° vertices, Dual Contouring | 4,868 | **232** |
| rows containing one, of 168 | 75 | **8** |
| **P-13** | falsified | **held** |
| **P-16** | falsified at 4× its threshold | **held, at 0%** |

Both falsified entries are still in `FINDINGS.md`, with their original wrong text, under a banner
pointing at the correction. Deleting them would delete the evidence that the instrument can be the
thing that is broken.

The rule it earned, now in the method section: **a reference implementation used as ground truth needs
the same scrutiny as the thing it checks.** And its shorter form — *when a measurement is impossible,
suspect the instrument before the world.*

P-12 is the same lesson wearing different clothes. It predicted the dual's superlinearity was a cache
effect; the miss rates came back flat and it was recorded falsified. That was **a correct measurement
of a machine whose bottleneck was somewhere else.** Once A-023 removed the store-to-load stall that was
covering it, the 128³ spike went from 1.24× — indistinguishable from noise — to 2.6×, and A-024 then
took it to 1.01×. A null measured under a dominant confound is a statement about the confound.

---

## The paper is wrong, and here is the fixture

Schaefer, Ju & Warren, *Manifold Dual Contouring* (IEEE TVCG 13(3), 2007), §3, on the uniform grid:

> *"this surface is always a manifold because the original MC algorithm always constructs a manifold
> and the dual preserves the topology of the surface."*

Over eight reference fields at three resolutions, no chunking and no weld, so there is only one
mechanism in the count:

| | non-manifold edges | non-manifold vertices |
|---|---|---|
| Marching Cubes | **0** | **0** |
| Marching Cubes + asymptotic decider | **0** | **0** |
| Manifold Dual Contouring, this crate's default table | 143 | — |
| Manifold Dual Contouring, the paper's decider-modified table | **114** | 222 |

The premise holds here. The conclusion does not.

Every one of those defects is on `noise_cavity`, and the mechanism now fits in 48 samples. A `2×2×3`
column of `±1` inside a `4×4×3` lattice, trilinear between the samples, everything else outside — its
middle plane reads `out, in, in, out`, the face saddle. On those twelve hand-written values Marching
Cubes emits 40 triangles and **0** non-manifold edges; both duals emit 20 triangles and **1**, carrying
four distinct faces.

Four quads meet on one dual edge because all four cut edges of the shared face land in one cycle on
each side, so all four connect the same pair of vertices. And the manifold construction buys nothing
here: splitting a cell by cycle cannot split a cell that has one cycle.

The sharper half is that **the offending set is not a set of sign configurations.** Scale the shared
face's two inside corners — no sign moves — and the asymptotic decider's saddle crosses zero with the
defect:

| inside corners | saddle `s` | the face | non-manifold edges | triangles |
|---|---|---|---|---|
| `−0.25` | `+0.375` | separated | **1** | 20 |
| `−1` | `0`, an exact tie | separated | **1** | 20 |
| `−4` | `−1.5` | **joined** | **0** | 20 |

One configuration, two answers, chosen by the field rather than by the rule. Which is why a face rule
still leaves 25–49 offending faces per resolution while being combinatorially capable of avoiding all
of them — over all 4,096 two-cell sign patterns, **not one** offends under every consistent mask.

---

## Reading the tiers

Every number in this repository's documentation carries one of these, and the legend matters because
the tiers are not decoration:

| tier | means | the bar |
|---|---|---|
| **M** | measured here | We ran it. Code and numbers in this repo, reproducible by checkout |
| **V** | verified externally | We read the primary source ourselves. DOI or file attached |
| **R** | reported | A credible source asserts it; we have not checked |
| **F** | folklore | Widely repeated, no verified primary source found |
| **✗** | falsified | Tested and found false — including what we believed, and why |

**Never cite an R or F claim as justification for a design decision without saying which tier it is.**

Two standing rules keep the ledger honest. **Re-tier rather than rewrite** — when an R becomes an M,
the R text stays and the measurement goes below it, because the gap between what was reported and what
was measured is itself data. And **record the ones we got right for the wrong reason**: ✗10 sits in the
falsified section even though its outcome was fine, because the reasoning was wrong and the reasoning
is what generalises.

---

## Where the raw results live

| | |
|---|---|
| the registrations | [`crates/isomesh/src/experiment.rs`](../crates/isomesh/src/experiment.rs) |
| the full ledger, 512 entries | [`FINDINGS.md`](../FINDINGS.md) |
| per-experiment result data | [`docs/experiments/`](experiments/) — each CSV carries its `# hypothesis:` and `# falsified by:` header |
| timing and quality measurements | [`docs/measurements/`](measurements/) |
| what falsified beliefs cost | `FINDINGS.md` Part 1 — 19 entries, and 79 falsified beliefs across the whole ledger, each recording where the wrong belief came from, because provenance predicts the next error |
