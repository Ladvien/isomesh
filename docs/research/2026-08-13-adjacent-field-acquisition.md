# Adjacent-field acquisition sweep — five domains, five novelty questions

**Date:** 2026-08-13
**Method:** five parallel agents over home-still (9,132 docs / 276,918 chunks). Local corpus first via
`distill_search`, absence confirmed by `catalog_read` (never by search), gaps filled by `paper_search`,
then `paper_download` → `scribe_convert` → `distill_index`.
**Result:** ~65 papers acquired. **Five novelty questions asked independently. All five came back "no."**

---

## The headline

Each agent was given one falsifiable novelty question and told to answer plainly, including "yes, here
it is." None found prior art.

| # | Question | Answer | Nearest prior art |
|---|---|---|---|
| 1 | Does **incremental isosurface extraction** exist? | **No** | Time-varying *indexing* of precomputed series (T-BON); incremental construction toward convergence on a *fixed* field; dirty-region invalidation with no bound |
| 2 | Has anyone used **persistence to resolve MC ambiguity**? | **No** — arXiv returned literally zero | Kissi/Tierny 2024 and Brüel-Gabrielsson 2018 both simplify **the field globally**, then contour. The hypothesis is a **local decision inside the cell** |
| 3 | Has anyone driven **voxel LOD by an a posteriori estimator or LFS** rather than camera distance? | **No** | Error-driven octrees exist but **offline only** (Manson & Schaefer, Koschier). Runtime LOD is universally Luebke's `p = εx / (2d·tan(θ/2))`, where ε is an *a priori* bound |
| 4 | Is there a **curvature-aware (second-order) vertex placement** rule for dual methods? | **No** | Both ingredients published separately — jet/Hessian fitting (Cazals & Pouget, Jiao & Zha) and QEF placement (Ju et al.) — never composed |
| 5 | Has anyone modelled **chunk consistency as a sheaf condition**? | **No** — the two literatures do not cite each other | SEAM (2026) does the right *shape* for PDE domains; applied-sheaf work has settled into sensor fusion, signal processing, multi-agent consensus, GNNs, networking. **Geometry is absent from that list** |

Five independent "no"s is not five findings of equal weight — but it does say the axes decomposition
was pointing somewhere real rather than at crowded ground.

---

## The five acquisitions that matter most

**1. `10.48550/arXiv.2105.06712` — Efficient Parallel Self-Adjusting Computation** *(already in corpus)*
Acar, Anderson, Blelloch, Baweja. Bounds update cost by a **computation distance** between two
executions, and gives **work *and* span**. This is the exact proof shape for edit-proportional
re-meshing, and it is already parallel — which matters, because chunk meshing already is.

**2. `10.48550/arXiv.1406.4005` — Maintaining Contour Trees of Dynamic Terrains** *(new)*
The single best structural idea from the sweep: **don't maintain triangles — maintain the contour
tree.** Maintains scalar-field level-set topology under a *changing* field, O(log n) per certificate
failure, certificates fail only on adjacent-vertex value swaps or saddle collisions, and it explicitly
handles general update operations rather than only continuous motion. Nearest existing analogue to
what we want, one dimension down.

**3. `10.1111/cgf.13933` — Probabilistic Quadrics** *(already in corpus, invisible to search)*
Trettner & Kobbelt. **The `catalog_read`-before-claiming-absence rule paid for itself again.** It says
outright that quadric minimisation "is in many cases not robust and requires an SVD or some ad-hoc
regularization," then derives a closed form solvable by a plain linear system — **50× faster than
SVD** — and demonstrates it *on isosurface extraction*. It reframes our problem: rather than reaching
for Tikhonov or truncated SVD, treat the input planes as uncertain and the regularization falls out.
**This supersedes the `λ ≈ 0.01` hack in the audit doc.**

**4. `10.4310/hha.2022.v24.n1.a16` — Cellular Sheaves of Lattices and the Tarski Laplacian** *(new)*
Ghrist & Riess. The most elegant find in the sweep: **lattice-valued sheaves make "global section =
consistent global state" and "lattice = join-semilattice = CRDT merge" the same object.** That fuses
both halves of the chunk-consistency problem — geometric seams and concurrent multiplayer edits — into
one formalism. Riess's thesis (`arXiv 2304.02568`) is the long form.

**5. `10.1017/s0962492924000011` — Adaptive Finite Element Methods, Acta Numerica 2024** *(new)*
One self-contained modern survey covering residual estimators with two-sided bounds, Dörfler marking,
newest vertex bisection, and full rate-optimality. Replaces most of the closed FEM book literature.

Also worth naming: **`10.1016/j.cad.2020.102856` — Attene, Indirect Predicates**. Reframing that
matters — an edge-crossing point is a *construction*, not an input. Keep it symbolic as (line, plane)
and get exact sign tests at near-float speed. Directly answers "what if the field is exactly zero at a
corner."

And **`10.48550/arXiv.2606.04227` — Volk, Incremental Sheaf Cohomology, O(1)-in-n Lazy Edit
Processing** *(already in corpus)*: literally maintains H¹ under a stream of local edits, deferring
global assembly to sync points. Caveats stated up front — 1-dimensional complexes only, evaluated on
Barabási–Albert graphs not grids, single-author 2026 preprint, unverified.

---

## Three cautions to carry, before any of this becomes a claim

Each agent volunteered the objection a reviewer would raise. These are the honest parts.

**On persistence (question 2).** Persistence is defined for a filtration of a *space*. MC ambiguity is
a question about one level set of a trilinear interpolant inside one cell — and the ambiguous cases
are ambiguous *precisely because eight corner samples underdetermine the interpolant*. So a tunnel's
persistence depends on which interpolant you assume. **This does not remove the modelling choice; it
relocates it** — from a hardcoded 730-entry table into a single tunable, stability-backed threshold.
That relocation is the real contribution and should be claimed as exactly that.

**On sheaves (question 5).** For a regular grid decomposition the nerve is simple, and a cover of
contractible boxes will often have **vanishing H¹ for the naive constant sheaf** — the machinery would
then be describing a problem that isn't there. The obstruction becomes non-trivial only once
restriction maps are non-identity: differing LOD across a seam, differing sign conventions, concurrent
edits with different clocks. Which is exactly the regime Volk's locality result does *not* cover.
**Confirm the interesting case is the one we actually have before betting on the formalism.**

**On FEM refinement bounds (question 3).** The newest-vertex-bisection closure bound
(`N − N₀ ≤ C·M`) is **global-over-history, not per-edit**. It bounds cumulative closure cost across
the whole refinement history against cumulative marks. It does *not* say a single local edit
propagates a bounded distance — which is the property a chunked engine actually needs.

---

## Infrastructure problems found, which cost real coverage

**`paper_download` resolves arXiv + Unpaywall only.** Every `10.1145/*` (ACM) and `10.1109/*` (IEEE)
failed, and so did open-access venues like LIPIcs and JoCG. This is not a paywall problem — it is a
resolver problem, and it blocked the single most on-point prior art in existence:

> **Acar's own incremental *meshing* line** — *Dynamic Well-Spaced Point Sets* (`10.1145/1810959.1811011`)
> and *Kinetic Mesh Refinement in 2D* (`10.1145/1998196.1998254`). Output-sensitive incremental meshing,
> by the person who built the theory. No arXiv version. **These need sourcing by other means.**

Also lost the same way: Acar's *A Cost Semantics for SAC* (the trace-distance paper), Carr et al.
*Flexible Isosurfaces* (canonical "which level-set features to keep"), Nielson's *On Marching Cubes*,
the original *Asymptotic Decider*, Smith/Levien/Owens *Decoupled Fallback* (single-pass GPU scan on
**exactly our platform** — Apple M-series, WebGPU), and Sorensen et al. on GPU progress models (which
finds **Apple and ARM GPUs do not support the linear occupancy-bound model** — a direct constraint on
any persistent-threads design on Metal).

**Corpus pollution — needs purging.** Several DOIs resolved to publisher landing pages and were
embedded as navigation boilerplate: `10.1090_s0025-5718-07-01959-x` (7 chunks of nav text),
`10.1007_s00371-007-0163-2` (6), `10.1137_060675666` (1), `10.1006_acha.1997.0238`,
`10.1007_978-3-642-24550-3_29` (the canonical Shapiro CRDT paper — 3.7 KB stub). These will surface as
false hits.

**`scribe_convert` is unreliable.** olmocr returns `workspace/markdown does not exist`; the MCP call
times out at 60 s but **the server usually completes asynchronously via the `glm_ocr` fallback**. Do
not retry on the error — wait and confirm with `catalog_read`. More than ~2 parallel conversions makes
failures much likelier. `pipeline_drift` is 80 against a threshold of 3.

**One whole subtopic came back empty:** lifting scheme / integer wavelets — all four canonical papers
behind Elsevier/SIAM. That is a real hole, not an oversight.

---

## What I'd do next, in order

1. **Read `2105.06712` and `1406.4005` together.** They answer "what is the theorem?" and "what is the
   object to maintain?" respectively. If the contour tree is the right maintainable object, that
   reframes axis 8 entirely — from "re-mesh fewer triangles" to "maintain topology, re-derive geometry."
2. **Replace the audit's `λ ≈ 0.01` with Probabilistic Quadrics.** It is already in the corpus, it is
   closed-form, it is 50× faster than SVD, and it was demonstrated on isosurface extraction. This is a
   concrete improvement to A-007/A-008 available today.
3. **Source Acar's two meshing papers by other means.** They are the closest existing prior art to the
   top-ranked transfer and the resolver cannot reach them.
4. **Purge the five polluted catalog entries** before they contaminate a future search.
5. **Before betting on sheaves, check H¹ is non-trivial for our actual decomposition.** One
   afternoon's calculation, and it decides whether the formalism is describing a real obstruction.
