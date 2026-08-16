# Phase 17 — SOTA

**Written 2026-08-16.** Every prior `R-` ticket is archived. These are the directions the research
turned up that **never became work** — verified: *self-adjusting*, *change propagation*, *contour
tree*, *persistent homology*, *dynamic connectivity* and *second fundamental form* return **zero**
occurrences across `BACKLOG.md` and `BACKLOG_ARCHIVE.md`.

Each ticket states **the gap and why it is unoccupied**, a pre-registered hypothesis, the harness,
the falsifier, and **what it is worth if it holds**. Phase 15's protocol applies: `P-` entry in the
commit *before* the measuring commit; `M-`/`✗`/`E×-` in the same commit as the result.

**A note on tier.** These are tier **F** — hypotheses. Several are one negative result away from being
closed, and closing one is a finding. Do not treat a null as a wasted ticket.

---

## 17a — The empty column

### R-020 · Isosurface extraction as a self-adjusting computation

**The gap, and it is the sharpest one found.** Across the ten-axis decomposition, axis 8 —
incrementality — reads **"full re-mesh" for every published algorithm.** An independent sweep of arXiv,
OpenAlex, CrossRef, Semantic Scholar and CORE found **incremental isosurface extraction does not
exist**. What exists is three adjacent things, none of which is it: time-varying *indexing* of a
precomputed series (T-BON), incremental construction toward convergence on a *fixed* field, and
dirty-region invalidation with no output-sensitivity bound — which is what every voxel engine and TSDF
system actually ships, at chunk granularity, with nothing proved about it.

**Axes where everyone differs are crowded. The axis where everyone is identical is where the space is.**

**The theory is in the corpus.** `10.48550/arXiv.2105.06712` (Acar, Anderson, Blelloch, Baweja) bounds
update cost by a **computation distance** between two executions and gives **work *and* span** — the
proof shape for edit-proportional re-meshing, already parallel, which matters because chunk meshing
already is.

**H:** re-meshing after a local field edit has computation distance `O(|edit|)` — i.e. the recorded
trace changes proportionally to the number of cells touched, not to grid size. **E1 is already
measured and supports this**: M-33/M-50 put a brush at **15–36%** of cells in its own bounding box,
reproduced live under a mouse.

**Harness:** instrument the extractor to record a dependency trace; apply a single-cell edit; measure
trace delta against edit size across grid sizes. **Records:** trace delta vs `|edit|` vs `n`.

**Falsified by:** trace delta scaling with `n` rather than `|edit|` — which would mean the extractor is
*unstable* in Acar's sense and edit-proportional repair is provably unavailable, closing the direction
with a reason. That is a real finding, not a failure.

**Worth if it holds:** the first isosurface extractor with a proved incremental bound. **Risk to
flag:** the closest prior art — Acar's own *Dynamic Well-Spaced Point Sets* (`10.1145/1810959.1811011`)
and *Kinetic Mesh Refinement in 2D* (`10.1145/1998196.1998254`) — is **unobtainable through the
pipeline** and does output-sensitive incremental *meshing*. Read before claiming novelty.

### R-021 · Maintain the contour tree, not the triangles

**The reframe from `10.48550/arXiv.1406.4005`:** the maintainable object is scalar-field level-set
topology, `O(log n)` per certificate failure, with certificates failing only on adjacent-vertex value
swaps or saddle collisions — and it handles general update operations, not just continuous motion.

**Caveat that must be carried:** that paper is **2-manifolds only** (`h: ℝ² → ℝ`, a triangulated
terrain). The two 3D results the question hinges on — Tarasov & Vyalyi 1998 and Safa & Wang 2014 — are
**both unobtainable**, and Edelsbrunner's 3-manifold Reeb maintenance is `O(n)` per certificate
failure, asymptotically no better than rebuilding.

**H:** the contour tree of a chunk can be maintained under a brush edit in time proportional to the
dirty set, where full recomputation is not. **Falsified by:** maintenance cost tracking chunk volume.

**Worth if it holds:** it changes axis 8 from "re-mesh fewer triangles" to "maintain topology,
re-derive geometry" — a different algorithm class.

### R-022 · Dynamic connectivity on the air sublevel set

**The cheap half of R-021, and it is buildable now.** The questions a game asks — *is this sealed? did
I break through? is this a chokepoint?* — are **not** all-thresholds queries. They are single-threshold
questions about connected components and bridges of the air sublevel set. That is dynamic
connectivity, which **is** measured and **is** cheap: microsecond queries, millions of updates/sec,
`O(log V)` depth when the spanning forest is unchanged — the common case, since most digging does not
alter connectivity. Bridges (chokepoints) are polylog amortized.

**The unoccupied part: dynamic connectivity has never been run on a voxel lattice.** Every measured
system was benchmarked on social/web graphs — Twitter 81K vertices, Stanford 280K. A voxel air-graph is
a **bounded-degree 6-connected lattice with 10⁶–10⁹ vertices.** Bounded degree should help; sheer V may
not. **And batching is untouched** — games edit thousands of voxels per explosion, not one;
`10.48550/arXiv.2002.05129` (batch-dynamic trees) is the right tool, is in the corpus, and has never
been pointed at this.

**H:** batched dynamic connectivity sustains a brush-sized edit under 1 ms at 128³.
**Falsified by:** per-edit cost scaling with lattice size rather than with the dirty set.

**Worth if it holds:** breakthrough-as-an-engine-event and sealed-volume-as-a-predicate, neither of
which any engine has. The benchmark alone is a contribution, since nobody has published it.

---

## 17b — One disputed table, replaced by one scalar

### R-023 · Persistence-thresholded ambiguity resolution

**arXiv returned literally zero** for persistence applied to Marching Cubes ambiguity. The nearest work
(Kissi & Tierny 2024; Brüel-Gabrielsson 2018) simplifies **the field globally**, then contours. This
hypothesis is a **local decision inside the cell**, and nobody has tried it.

**It is live right now.** A-002i and A-020b are blocked on architecture for exactly the singular and
tunnel cases, and A-002b's own reasoning is that *there is no correct published table to transcribe* —
Custodio et al. proved Chernyaev's interior test tracks a quadratic where the true saddle trajectory is
hyperbolic, and Lewiner's reference implementation omits cases 10 and 12 entirely. **The guaranteed
version is 730 subcases.**

**The move:** stop asking *"is there a tunnel"* and ask *"does this tunnel have persistence above ε."*
Below threshold, mesh closed; above, mesh open. A disputed 730-entry table becomes a **computable
scalar with a stability theorem behind it** — and a knob a game wants.

**Scope it honestly.** Persistence is defined for a filtration of a *space*; the ambiguous cases are
ambiguous precisely because eight corner samples underdetermine the trilinear interpolant, so a
tunnel's persistence depends on which interpolant you assume. **This does not remove the modelling
choice — it relocates it**, from a hardcoded table into one tunable, stability-backed threshold. Claim
exactly that.

**H:** a persistence threshold reproduces MC33's topology on the fields where MC33 is agreed correct,
and differs only on cells where the published algorithms disagree with each other.
**Falsified by:** disagreeing on cells where Chernyaev and Lewiner agree.

**Worth if it holds:** it retires A-002b, A-002i and A-020b together, and it is the one idea here
nobody has published even a negative result on.

---

## 17c — Two things sitting on this crate's own seam

### R-024 · Does field-sealed imply mesh-sealed?

**Nobody has established this, and every paper treats the two as interchangeable.** A cell can be
topologically connected in the field and still produce a watertight surface, or the reverse, depending
on the case table and the interpolant.

**This is a day of work and it is publishable alone:** extract a mesh; compute connected components of
the air sublevel set; compute connected components of the mesh complement; **assert they agree.**

**H:** they agree for Marching Cubes on all eight reference fields, and **disagree for at least one
dual method** — the duals place vertices by solve rather than on the interpolant, which is exactly
where the correspondence could break.

**Falsified by:** universal agreement — a stronger correctness statement than the crate currently
makes, and worth saying so.

**Worth either way:** it is the precondition for every mechanic in R-022, and the gap sits exactly on
the seam this crate occupies.

### R-025 · Second-order vertex placement

**Both ingredients are published separately and have never been composed.** Jet/Hessian fitting
(Cazals & Pouget; Jiao & Zha, in corpus) and QEF placement (Ju et al.) — nobody fits the **second
fundamental form** per cell instead of tangent planes.

**It attacks a measured term.** P-2's error model is `O(|e|²κ)`, and on a true SDF the Hessian's
nonzero eigenvalues at a surface point are `−κᵢ/(1−κᵢd)` — **principal curvatures fall out of samples
already taken**, no medial axis involved.

**H:** curvature-aware placement improves Hausdorff on smooth fields (`sphere`, `torus`, `gyroid`) by
>20% over planar QEF at fixed resolution, and **does not** improve it on `box_exact` — where the
surface is flat and the second-order term is zero.

**Falsified by:** no improvement on smooth fields, meaning curvature estimation is too noisy at game
resolutions — which is Aamari & Levrard's minimax bound biting, and worth recording as such.

---

## 17d — The result you already have and have not claimed

### R-026 · Write up the head-to-head

**M-001 produced the comparison that does not exist in the literature**, and M-004's writeup ticket is
archived while the paper is not written. Verified: **no paper since 2020 benchmarks Marching Cubes vs
Surface Nets vs Dual Contouring against each other**, and Surface Nets — the thing engines actually
ship — **has no credible published timings at all.**

You additionally hold results that **contradict** published figures: M-51 and M-55 falsify the
literature's `2–3×` Marching Tetrahedra ratio (measured `~3×` triangles for `4.3%` worse geometry, not
86%), M-1's `V_sn = V_mc + χ` identity, M-53's four-corner table of manifold × intersection-free, and
M-54's `101×` Dual Contouring accuracy advantage on sharp fields.

**This is the least speculative item in the phase and the only one whose result is already in hand.**
The remaining work is Open SciVis volumes for comparability (H-005), mesh-quality metrics for the
table reviewers expect (H-003), and prose.

---

## Ordering

| | Why |
|---|---|
| **R-024** | One day, publishable alone, and it gates R-022 |
| **R-026** | The result exists; only the writing is missing |
| **R-022** | Buildable now on measured foundations; the benchmark itself is unpublished |
| **R-023** | Retires three blocked tickets at once, and nobody has even a negative result |
| **R-020** | The biggest space, and the one most at risk from unread prior art — get Acar's two papers first |
| **R-025** | Cleanest hypothesis, most likely to null out honestly |
| **R-021** | Highest ceiling, worst evidence base — two load-bearing papers unobtainable |
