# Experiments that held, and what they buy you

Seventeen times, this project wrote down what it expected to measure — and what result would prove it
wrong — **before running the measurement**. Ten of those are enforced by the compiler.

This page is the scorecard. Five of the ten compiler-enforced predictions held, four were falsified,
one never ran, and the most useful thing in the whole record is a bug that made two true hypotheses
look false for a day.

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
| **P-70** — subgroup ballot compaction beats the Hillis-Steele scan | **falsified** on C1 by arithmetic done at registration; **held** on C2 and C3 (✗54 / M-381) | **The mechanism is real — 1.4376× at 2 M elements, bit-identical, reproducing the literature's 1.33–1.49× band — and it is not shipped.** Applied to the path it lives in it moves **1.33%**, because `scan_ms` is 0.3657 of 8.3694 while `upload_ms` is 87.50%. A second WGSL path for 1.33% is what the one-path rule forbids, and naga 29.0.4 rejects `enable subgroups;` outright |
| **P-71** — the 83% is a blocking round-trip, and both targets can avoid it | **falsified** on all three (✗52 / M-376) | **The largest piece of a GPU extraction is the geometry copy, not the count wait** — 0.6663 ms of 1.1331 at 129³ against map-wait's 0.3191 — so `extract_indirect` removes **32%**, not the registered 60%, and it removes it by *not delivering the bytes*. For a consumer that only draws, that is 100% of what it needs. For a collider consumer the copy is unavoidable, which is why `DeferredGeometry` ships in 0.0.10: measured under the frame-budget scheduler at **1.41 frames of latency, worst case 2** |
| **P-72** — the granularity of the active-cell structure is a parameter with an optimum | **held** on C1 at 51×, **falsified** on C2 and C3 (M-377) | **Chunk size is the largest single knob in this crate and its optimum is interior at 4³.** Eleven edits cost **4.56 ms at 4³ against 435.81 ms at 64³** at a fixed total cell count, because remesh is `dirty_chunks × chunk_cells³` — 25,344 cells re-meshed against 8,126,464. The registered 8³–64³ range would have given the wrong shape entirely: it is monotone there with the minimum at the smallest granularity swept. **The demos are not retuned on this**, because the fixture excluded per-chunk entity and draw cost |
| **P-62** — the Plantinga–Vegter certificate is sound against this crate's own tunnel classifier | **held** on C1 and C2, **falsified** on C3 (M-378) | **Zero unsound certificates over 2,389 tunnel cells**, from a predicate that refuses **95%** of the population — so the zero is not the trivial one. It is not free: the cost is 21%, not the registered ceiling, and the reason is structural rather than fixable in isolation — a standalone pass re-gathers the eight corners the extractor already has in registers |
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
| the full ledger, 434 entries | [`FINDINGS.md`](../FINDINGS.md) |
| per-experiment result data | [`docs/experiments/`](experiments/) — each CSV carries its `# hypothesis:` and `# falsified by:` header |
| timing and quality measurements | [`docs/measurements/`](measurements/) |
| what falsified beliefs cost | `FINDINGS.md` Part 1 — 24 entries, each recording where the wrong belief came from, because provenance predicts the next error |
