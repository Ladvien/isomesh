# Eight unmined lenses, and the six experiments they earned

**Date:** 2026-08-23
**Source:** eight parallel corpus sweeps over home-still (9,481 documents, 289,866 chunks), deliberately
pointed at lenses the 2026-08-18 novelty table and the 2026-08-17 mechanics dossier never used:
non-cubic sampling lattices, certified root isolation, bit-level and data-layout kernels, machine-checked
combinatorial verification, discrete differential geometry and digital topology, staged computation and
partial evaluation, kinetic and event-driven geometry, and player-visible by-products of extraction.
**Tier:** every row here is **R** — read by a sweep agent, not measured here — except where it cites an
existing **M** row, in which case the M row wins. Six of the candidates have been promoted to `P-`
registrations (P-38 … P-43) and become **M** or **✗** only when their harness runs.

**Excluded by design.** The 2026-08-18 novelty table's fifteen rows and its rejected list, and Phase 18's
banked mechanics, were handed to every sweep as a do-not-re-propose list. Where a candidate below touches
banked ground the overlap is named in its own row rather than hidden.

---

## What the lenses were for

The previous two rounds asked *what else does the isosurface literature do*. That question is close to
exhausted: eight sweeps over primal, dual, adaptive, placement, GPU, incremental, learned and
adjacent-math returned fifteen candidates, and the marginal one was already thin. This round asked a
different question — *what do other fields do that nobody has pointed at a mesher* — and the yield
pattern is the answer to whether that was worth doing. Three of the eight lenses produced a candidate
that is now registered; two produced candidates held in reserve; three produced mostly `OUT` rows, and
the `OUT` rows are the cheap half of the result because they foreclose attractive work.

**The single most valuable finding is not a mechanism.** `abstract_search` is unusable for this subject
area. Across the eight sweeps it returned zero relevant hits on **every** graphics or geometry query —
RANSAC and bipolar-disorder papers for a BCC-lattice query, Chaum's mix networks for a digital-topology
query, *Attention Is All You Need* for a CSG-pruning query — while `distill_search` returned rank-1 hits
on the same subjects. The graphics corpus carries no abstract text, so the abstract index degrades to
title-or-nothing and the known 37.3% null-title rate becomes effectively 100% for this region.
`openalex_search` was **down** for the whole session: 30-second SSE timeouts on every one of ~25 calls
across all eight agents. The mandated third-channel novelty check was substituted with `paper_get`,
`paper_search` and `paper_citations`, and every candidate below says which it got.

That caveat has teeth, and it fired: **Tilove 1984** — the 1984 CACM paper that is the correctness
provenance for the flagship candidate — ranked *third* on its query with `title: null` and would have
been invisible to a title-based reader.

---

## The six registered experiments, and the four the results earned

**Outcomes added 2026-08-23, after the fact and marked as such.** Everything below the table was written
before any of it ran; this table is the only part that knows how it came out.

| P- | Ticket | Bar | One-line claim | Outcome |
|---|---|---|---|---|
| **P-38** | R-037 | speed | Marching Cubes never received A-024's odd-stride fix, and pays the same 128³ aliasing penalty the dual path was cured of | **✗28 / M-336 — FALSIFIED** at 0.98 against 1.5. The stride is not the defect; the *walk* is |
| **P-39** | R-038 | speed | A brush that provably cannot win the min/max chain inside a chunk can be deleted from the tape, bit-exactly, before the chunk is meshed | **M-341 — HELD**, all three. 3.36× median, 64/64 bit-exact |
| **P-40** | R-039 | speed | The active-cell test is one bit per sample, so 64 cells decide at once in four fused word operations, with the same cells visited in the same order | **M-337 — HELD**, all three. Stage 5.5×, Surface Nets 1.34×, 12/12 bit-identical |
| **P-41** | R-040 | capability | The sign lattice this crate meshes is not well-composed, and that is where its dual extractors go non-manifold | **M-338 — HELD**, and the relation is a *bijection* |
| **P-42** | R-041 | capability | Curvature as a normal-cycle *measure* is additive over chunks and carries a computable error bound | **✗30 / M-340 — FALSIFIED.** With `B` the whole surface the clause was discrete Gauss–Bonnet, i.e. a tautology |
| **P-43** | R-042 | capability | One evaluation at the cell centre is a witness that a chunk's resolution is inadequate | **✗29 / M-339 — FALSIFIED**, both clauses, one of them by the registration's own arithmetic |
| **P-44** | R-042a | capability | …the *mean* residual instead, tested on the four fields P-43 never touched | **✗31 / M-342 — FALSIFIED.** Correlation reproduced; the exponent gap *is* the normalisation. Line closed |
| **P-45** | R-041a | capability | …additivity instead, with a real boundary on `B` | **✗32 / M-343 — FALSIFIED.** The two measures have opposite defects |
| **P-46** | R-040a | capability | …repair the sign lattice, now that the bijection gives it a target | **✗33 / M-344 — FALSIFIED.** Free and total when sparse, *provably* stuck when dense |
| **P-47** | R-043 | capability | Composed fields fall back to a six-sample central difference, and nobody has priced it | **✗34 / M-345 — FALSIFIED** in the mean by three orders; the 2.8× speedup held |
| **P-48** | R-044 | capability | Close the gap `isotopy.rs` names in its own header: certify the *analytic* field | **M-347 — HELD**, all three. Zero unsound over 1.07 × 10⁹ evaluations |
| **P-49** | R-045 | capability | Connectivity is a boolean; a game asks for a clearance | **M-346 — HELD**, all three. Exactly zero error where the answer is known |

**Five held, seven falsified.** The two reserve candidates that were promoted late — P-48 and P-49 — both
held, and neither was in the original six. The two that were *re-registrations of falsified experiments*
— P-44 and P-45 — both failed again, on data they had not seen, which is what a re-registration is for.
And two of the seven failures were errors in the registration rather than in the world: P-43's cost
clause was wrong arithmetic and P-42's Gaussian clause was unfalsifiable by construction. Neither was
repaired by amending the prediction.

**The lens that paid was not the one that looked most promising.** Staged computation produced the
flagship (P-39), digital topology produced the bijection (P-41), and bit-level kernels produced the
speedup (P-40) — but bit-level kernels' *own* rank-1 candidate, the aliasing fix, was the phase's first
null. The lenses that produced nothing usable were lattices and formal verification; their value was the
foreclosed list.

---

### P-38 · R-037 — the aliasing fix that reached one of two engines

**Not a literature candidate. An in-repo contradiction, found by reading both engines side by side.**

`dual.rs` carries `row_stride(size) = size[0] | 1` and a doc comment recording what it bought: *"Surface
Nets measured 3.37× the cost of 127³ or 129³ there — on a field with no surface, so it is the scaffolding
and not the geometry. 256³ pays 1.39×"* (A-024, M-287). `MarchingCubes` has no such thing. Its `values`
buffer is indexed through `shape.linearize` (`marching_cubes/mod.rs:687`), so at `size[0] = size[1] = 128`
the row stride is 512 bytes and the plane stride exactly 65,536 — the two aliasing periods A-024 exists to
break. Its `edge_vertices` cache is worse: `3 · sample_count` `u32`, a plane stride of `3 · 65,536` at 128,
and it is the buffer with the scattered access pattern.

**Why it runs first.** P-40 measures a traversal-stage speedup at 128³. If Marching Cubes is sitting in a
3× aliasing hole at exactly that size, any 128³ comparison is measuring the hole. This is the ✗14 lesson
in a new costume: *a null measured under a dominant confound is a statement about the confound.*

**Bit-identity is structural**, not hoped for: the change permutes where floats are stored. It does not
change which floats are computed, the order cells are visited, or the order vertices are created.

---

### P-39 · R-038 — Lipschitz tape pruning, and the lemma that makes it bit-exact

**The flagship, and the one candidate whose payoff grows with how long the player has been playing.**

`BrushStack::sample` is a linear fold over every brush in the stack, and `MarchingCubes::extract` prefills
the whole sample grid before any cell work, so a 33³ chunk evaluates the entire edit history 35,937 times.
In a destructible world the tape only grows. Nothing in the crate attacks this.

Bound each brush's shape over the chunk AABB — with a declared Lipschitz constant `l` the bound over a box
of circumradius `r` centred at `c` is `f(c) ± l·r`, so it is **one sample per brush per chunk** — and delete
every brush that provably cannot affect the result anywhere inside. Keeter measures the reduction at *"two
orders of magnitude"* and calls it *"the optimization that makes our algorithm practical"*; Barbier et al.
generalise it from intervals to the Lipschitz property and name this exact use case as future work:
*"could also lead to significant improvement for SDF discretization or polygonization."*

**This is not the tabled Sharp & Jacobson row.** That one rejects *cells* that contain no surface. This one
makes the cells that *do* contain surface cheaper by shortening the expression they evaluate. Orthogonal,
and composable.

**The lemma, derived from this repo's source rather than from any paper.** No paper in this literature
discusses IEEE bit-identity of pruning, because no paper in this literature has a 216-golden-hash gate —
Keeter and Barbier both prune to preserve a rendered image, which tolerates ULP drift.

- `apply(Add, f, s) = min(f, s)` and `apply(Subtract, f, s) = max(f, −s)`. IEEE `min`/`max` **select** an
  operand; they do not compute a new value, and negation is exact. So deleting a provably-losing `Add` or
  `Subtract` brush changes the result by **exactly zero ULP**. `brush.rs`'s own module docs already assert
  the selection property for the commutativity argument (M-36 … M-38).
- `smooth_min` is **not** bit-exactly prunable in the losing direction. At `h == 1` it returns
  `b + (a − b)`, which is not bit-identical to `a`. It *is* exactly dominant at `h == 0`. The registration
  carries that asymmetry, because discovering it in the benchmark would read as a determinism regression.

Sources: Barbier et al., *Lipschitz Pruning*, CGF 44(2) 2025, `10.1111/cgf.70057` — **not in corpus**,
verified by `paper_get`, downloadable. Keeter, *Massively parallel rendering of complex closed-form implicit
surfaces*, `10.1145/3386569.3392429` — **not in corpus**, double-checked by `catalog_read` and
`distill_exists`. Tilove, *A null-object detection algorithm for constructive solid geometry*, CACM 1984,
`10.1145/358105.358195` — **in corpus and indexed**, and it is the set-theoretic statement of the same
argument forty years early: primitive redundancy under spatial localisation.

---

### P-40 · R-039 — the active-cell test is one bit

`DualMesher::place_vertices` gathers eight corners and counts insides for **every** cell, and skips roughly
97% of them. Pack the corner sign into one bit per sample, 64 to a `u64` along x, and a whole word of 64
cells decides at once: with `s(y,z)` the row word,

```text
any    = OR  over the four rows (y,z),(y+1,z),(y,z+1),(y+1,z+1) of (w | w>>1)
all    = AND over the same four rows of (w & w>>1)
active = any & !all
```

— about twenty word operations for sixty-four cells, against sixty-four scalar gathers. Walk the set bits
with `while w != 0 { let b = w.trailing_zeros(); w &= w - 1; }`, which is safe stable `core`, `no_std`, and
lowers to `tzcnt`/`popcnt` on x86-64 and `rbit`+`clz` on aarch64 with no target-feature gate and no `unsafe`.

Mechanism source: Museth, *VDB*, `10.1145/2487228.2487235` — **in corpus**, 45 chunks, confirmed by
`catalog_read` after `paper_download` wrongly reported it unavailable. VDB's whole topology layer is
per-node bitmasks under word-level boolean ops; the transfer is to the active-cell predicate of a CPU
isosurface extractor with the mesh held bit-identical.

**One trap, pre-registered.** Do **not** extract the IEEE sign bit as the inside bit. `-0.0` has the sign
bit set and `-0.0 < 0.0` is false, and `box_exact` is exactly zero across its whole boundary, so signed
zeros are reachable in this crate's own fixtures. The comparison stays; the win is in the cell test.

**Amdahl is registered as a split prediction, not hidden:** field evaluation already dominates at high
resolution (M-296), so the whole-extractor ratio is claimed on `sphere` and the stage ratio is claimed on a
surface-free field where the active path never runs.

---

### P-41 · R-040 — the sign lattice is not well-composed

A binary digital set is **well-composed** exactly when its boundary is a 2-manifold, and Latecki
characterises that by the absence of small critical configurations: in 3D, a diagonal pair sharing only an
edge, or only a vertex. In a well-composed set 6-, 18- and 26-connectivity coincide, so *"is this cave
sealed"* stops depending on which connectivity was picked.

The move is **not** to repair the mesh. It is to repair the *sign pattern* upstream, before the case table
is ever indexed — a dual extractor emits non-manifold output exactly where the sign lattice is not
well-composed, and if the lattice is well-composed by construction those cases are unreachable.

**The registration is the detector, not the repair**, and that is deliberate. Two numbers decide whether
the repair is worth writing: the critical-configuration census over the eight reference fields, and the
fraction of the dual extractors' non-manifold incidents that land in cells the census flagged. Below 90%
co-location the sign lattice is not the cause and the whole line dies before an extractor is touched.

Sources: Latecki, Eckhardt & Rosenfeld, *Well-Composed Sets*, CVIU 1995, `10.1006/cviu.1995.1006`;
Latecki, *3D Well-Composed Pictures*, `10.1006/gmip.1997.0422`; Boutry, Géraud & Najman, *A Tutorial on
Well-Composedness*, JMIV 2018, `10.1007/s10851-017-0769-6` — all three DOI-verified via `paper_get`, all
three **not in corpus** (`distill_exists` false for the tutorial; four differently-worded corpus queries
returned nothing on target). Two near-miss DOIs were caught and discarded during verification:
`10.1006/cviu.1995.1013` is Brechbühler's SPHARM paper, and `10.1016/j.dam.2015.01.006` is a paper about
Turán numbers.

**Overlap, named:** this is not the banked union-find topological safety gate. That gates a vertex
*reposition*, after extraction, on the mesh. This is a census of the sign field, before extraction.

---

### P-42 · R-041 — curvature as a measure, with a bound

`MeshReport` reports `mean_ratio` and `euler_characteristic` and nothing that measures curvature or bounds
its own error. Normal cycle theory gives both, and gives them cheaply: on a triangle mesh the Gaussian
measure collapses to the vertex angle defect and the mean measure to edge length times signed dihedral
angle.

Angle defect is textbook and is **not** the claim. Two things are:

- **Additivity.** `N(A ∪ B) = N(A) + N(B) − N(A ∩ B)`, so the measure of a chunk is the sum over its cells
  with the shared boundary subtracted once. A per-chunk number composes into a per-world number with no
  global pass — which is the property that distinguishes this from every other curvature estimator and is
  exactly what a chunked mesher needs.
- **A stated error bar.** Cohen-Steiner & Morvan's Theorem 6 bounds the difference from the smooth
  surface's measure by `C·K·ε` with `K = Σ cr(t)² + Σ_{t∩∂B≠∅} cr(t)` and `ε = max cr(t)` — every quantity
  computable from the mesh alone. `validate::accuracy` measures distance and `validate::isotopy` certifies
  topology; neither says anything about a differential quantity, and neither states a bound it derived
  rather than sampled.

The falsifier is aimed at the interesting failure, not the boring one: Theorem 6 needs the polyhedron
closely inscribed in the surface with the projection a bijection, and marching-cubes vertices lie on the
*trilinear interpolant's* zero set, not on the field's. If the measured residual escapes its own bound,
that hypothesis has failed, and that is a better result than the bound holding.

Sources: Cohen-Steiner & Morvan, SoCG 2003, `10.1145/777792.777839` — DOI-verified, **not in corpus**; Sun
& Morvan, `10.5802/acirm.50` — **in corpus**, and it carries the transcribable formulas, the additivity
law and Theorem 6's explicit constants, so rule 5 is satisfiable without the primary.

---

### P-43 · R-042 — a witness that the chunk is under-sampled

isomesh has no way to tell a caller that a chunk's fixed resolution is insufficient for the field it is
carrying, so an edit that carves a feature thinner than a voxel silently produces wrong geometry with no
diagnostic. `validate::field_bound` *samples* the gradient and is explicit that a sampled maximum is a
lower bound on a supremum; it can disprove a declared bound and cannot certify sampling adequacy.

Evaluate the field once at each cell centre and compare against the trilinear interpolant of the eight
corners. The normalised residual is a directly computable witness of inadequacy, in the one-sided direction
that is safe: it can prove a chunk under-sampled and can never prove it adequate — the same discipline
`field_bound.rs` already uses, pointed the other way. Cost is exactly one extra evaluation per cell, 12.5%
of the corner evaluations, and it can run at reduced stride.

The falsifier is a correlation against a number the repo already produces: per-chunk centre residual
against the symmetric Hausdorff `validate::accuracy` reports. Below `r = 0.7` it is a witness of nothing.
The expected failure mode is named in advance — a feature that passes cleanly through the cell centre
without perturbing it, which is exactly what `thin_plate` is constructed to be.

Sources: Petersen & Middleton, `10.1016/s0019-9958(62)90633-2` (DOI-verified) for the sampling-theoretic
root; Chazal, Cohen-Steiner & Lieutier, *A Sampling Theory for Compact Sets in Euclidean Space*,
`10.1007/s00454-009-9144-8` — **in corpus and indexed** — for the version that tolerates non-smooth shapes,
which is what `csg_difference` and `box_exact` need.

---

## Held in reserve, with the reason

Ranked, and each one is a registration away if a run above comes back null.

| Candidate | Bar | Why it is not registered yet |
|---|---|---|
| **Nielson's `DisC[T]`/`G[T]` invariants as a certified surface-free sub-box** (`10.1109/tvcg.2003.1207437`, in corpus, title null) | capability | The strongest reserve. `trilinear.rs` solves one quadratic and derives the rest by linear solve; the cuboid needs the roots of all three, which is real extra arithmetic. Also the corpus chunk's printed 12-term discriminant contains an OCR corruption (`-2bdjf`), so the coefficient must be re-derived, not transcribed |
| **Interval-valued `Sdf`, closing the gap `isotopy.rs` names in its own header** | capability | `smooth_min` needs a sound enclosure and it is not monotone in both arguments; Fryazinov et al. (in corpus) derive exactly this and call conditional operators *"a very complex task"* |
| **Forward-mode dual numbers over the brush tape** | capability | `BrushStack` and `Capsule` do not override `Sdf::gradient`, so any composed field silently falls back to the 6-sample central difference at `O(h²)`. Real hole; but it changes `mesh.normals`, so it is a fourth `NormalStrategy` measured against the third, never a replacement |
| **Deterministic simplicial collapse as a "this chunk is a ball" certificate** | capability | One-sided by nature — `(1,0,0,0)` proves a ball, anything else proves nothing — and Benedetti & Lutz's heuristic succeeds because it is *random*, which this crate cannot be. Falsifier is sharp; payoff is narrower than it first reads |
| **A15 acute-tetrahedra marching stencil** (`10.1145/2504459.2504507`, in corpus) | capability | Would let the crate state a lower bound on output triangle angle for the first time. Two costs: the mesh changes, so a golden rebaseline; and the tile coordinates must come from Eppstein/Sullivan/Üngör, which is not in corpus and has no open PDF |
| **Metamorphic relations over the octahedral group** | capability | 48 elements × 8 fields × 7 extractors, all exact after remapping. Held because the granularity is subtle: `table.rs` picks a `safe_apex` by lowest edge index, which is *not* invariant under axis relabelling, so the relation holds at cut-edge-cycle level and not at triangle level. Stating it wrong manufactures 2,688 false positives |
| **Per-vertex edit provenance** (Tilove spatial localisation) | capability | Free once P-39 lands — the surviving-brush set *is* the provenance. Deliberately not registered separately; it rides P-39 or not at all |

---

## Foreclosed, and why — the cheap half of the result

Each of these looks transferable under its lens and is not. Recording them stops the next sweep paying for
them again.

- **PEXT/PDEP-driven compaction.** `core::arch::x86_64::_pext_u64` is an `unsafe fn`, x86-64 only, gated on
  `target_feature = "bmi2"`, needing either `std`'s runtime detection or a compile-time flag this crate
  cannot impose on consumers. Three independent violations; the `unsafe` alone is disqualifying. The safe
  substitute is already P-40's set-bit walk, and on a bitmap that is 97% zeros the cost is dominated by the
  words skipped entirely, not the bits extracted.
- **Morton / Hilbert layout for the sample grid.** Bit-identity *is* preservable — storage layout can change
  while traversal stays lexicographic — but without PDEP the encode is ~9 shift/AND/OR per coordinate per
  access, strictly more work, to buy locality that P-38's single `| 1` already buys for one bit. The
  surviving form is Morton ordering *across* chunks in the scheduler, which touches no mesh bytes.
- **Polyhedral model, time skewing, diamond tiling.** Every win in that literature comes from fusing
  iterations *over time*. Isosurface extraction is a single-sweep stencil with no `t` loop to skew; at
  `T = 1` the machinery collapses to plain 3-D blocking. The one transferable idea is *finalization* — the
  certificate that a z-plane can be retired — which is what would make a slab window provably safe.
- **Minimal perfect hashing for the edge cache.** `edge_vertices` is already a perfect hash: direct
  addressing on `(sample, axis)`, one probe. It is merely not minimal, and an MPHF needs the key set up
  front, which is only known after marching. A two-plane slab window gets the same footprint reduction with
  zero construction cost.
- **Interval Newton, Krawczyk, α-theory, certified homotopy continuation.** Disqualified by the
  unbounded-iteration rule, and the disqualification is stated in the corpus rather than assumed: Cheng &
  Wen (in corpus) write that *"the termination of the method is not ensured"* and their own algorithm emits
  a *"suspected root box"* set — a second, degraded answer path, which is the one-path rule's exact target.
- **Cylindrical algebraic decomposition.** Global, doubly exponential, needs arbitrary precision. Fails
  chunk-locality, the frame budget and rule 3 simultaneously.
- **Sturm sequences.** For the degree-2 and degree-3 polynomials that actually arise on a cell they are
  strictly more expensive than the discriminant and Bernstein sign counts, and buy nothing. Out on cost.
- **Circle packings, discrete conformal structures, heat/wave kernel signatures.** All in corpus, all
  requiring either a convex optimisation over the whole surface or a Laplace–Beltrami eigendecomposition.
  Not chunk-local. Excluded by the constraint, not by taste.
- **Cup products and cohomology rings.** Dropped on mathematical grounds rather than availability: for a
  surface embedded in ℝ³ the ring adds nothing over `(χ, orientability, β₁)`, so the discriminating power
  the lens hoped for does not exist in this dimension.
- **Runtime codegen / JIT specialisation of the SDF tape** — the literal first Futamura projection. Needs a
  code generator and W^X memory: a dependency and `unsafe`. P-39's reserve companion, the two-stage
  `Prepared` brush, is the same idea reachable inside the constraints.
- **`egg` / equality saturation as a runtime component.** A crate, a `HashMap` with nondeterministic
  iteration order, and an unbounded saturation loop. Four violations. Its offline use to *derive* a rewrite
  is admissible, but an arithmetically equivalent rewrite of float code is generally not bit-equivalent —
  the in-corpus *Sparsity-Specific Code Optimization* paper says so in as many words — so the yield is near
  zero against a golden-hash gate.
- **Superoptimisation of the case kernel.** The observation worth banking is that a 256-entry sign domain
  makes superoptimisation a *proof* rather than a heuristic. The instance is still Amdahl-dead: it optimises
  the part of the pipeline that is provably not the bottleneck.
- **Blue-noise and Poisson-disk sample placement.** Destroys the marching-cell structure and needs a
  neighbour search that is not a fixed-stencil chunk-local operation. Lloyd relaxation is already rejected.
- **Aperiodic Wang / corner tiles to break the diagonal grain.** The transfer — hash of grid coordinates
  selecting a per-cell stencil — *is* the already-banked randomised tetrahedral split. Recorded so the lens
  is not pointed at it a third time.
- **HCP sampling.** Not a lattice (not closed under addition), so it cannot host a periodic `const` stencil.
  Out on structure, independent of literature.

---

## Method notes

**Every "this paper is absent" claim here was checked, and one was wrong before it was written.**
`paper_download` is not a presence oracle: it returned *"No open-access PDF found"* for Museth's VDB, which
`catalog_read` then confirmed is present, converted and embedded with 45 chunks. Absence claims below rest
on `distill_exists` or `catalog_read`, never on a failed download.

**Two DOIs in the literature are traps and are recorded so nobody re-derives them.**
`10.1109/VISUAL.2001.964496` is *not* Theussl's *Optimal Regular Volume Sampling* — it is a tennis
visualisation paper; the correct DOI is `10.1109/visual.2001.964498`. And the *Dual Contouring of Hermite
Data* DOI that three earlier sweeps cited, `10.1145/566654.566586`, does not resolve; the paper is at
`10.1145/566570.566586`, which is what `hermite.rs` already cites correctly.

**Read the paper before the row.** P-41 and P-42 both rest on primaries that are *not* in the corpus, and
both are one `paper_download` away. Rule 5 is satisfied for P-42 by the in-corpus Sun & Morvan companion,
which carries the formulas; P-41's repair — as opposed to its census, which is registered — must not be
written from the tutorial's summary.
