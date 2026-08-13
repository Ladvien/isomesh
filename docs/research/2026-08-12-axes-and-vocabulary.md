# Meshing algorithms decomposed into axes — and the words for each

**Date:** 2026-08-12
**Purpose:** two things at once. A decomposition of "an isosurface extractor" into independent design
axes, and the **vocabulary** for each — because the single biggest obstacle to raiding another field
is not understanding its math, it's not knowing what it calls your problem.

**How to read the tiers.** Vocabulary is stable and safe. The *transfer opportunities* are mine unless
marked otherwise, and they are hypotheses — tier **F** in `FINDINGS.md` terms until someone checks the
literature. I've marked which is which.

---

## Why axes at all

"Marching Cubes vs Dual Contouring" is the wrong unit of comparison. Every extractor is a **pipeline
of ten roughly independent decisions**, and the famous algorithms are just popular *combinations* of
choices. Surface Nets and Dual Contouring differ on exactly one axis (vertex placement); Marching
Cubes and Marching Tetrahedra differ on exactly one (domain decomposition).

The consequence: **the design space is far larger than the named algorithms cover**, and the unoccupied
regions are where novelty lives. If there are 10 axes with 3–5 published choices each, the named
algorithms occupy maybe 15 points in a space of ~10⁵. Most of the space is empty because nobody
enumerated it, not because it's bad.

---

## The ten axes

| # | Axis | The question it answers | Marching Cubes | Surface Nets | Dual Contouring |
|---|---|---|---|---|---|
| 1 | **Domain decomposition** | How is space partitioned? | cubic grid | cubic grid | cubic grid |
| 2 | **Sign inference** | What decides inside/outside? | corner signs | corner signs | corner signs |
| 3 | **Ambiguity resolution** | What happens when signs are consistent with two topologies? | nothing (MC33: asymptotic decider) | n/a — one vertex per cell | n/a |
| 4 | **Vertex placement** | Where exactly does the vertex go? | linear interp on edge | centroid of crossings | QEF minimizer |
| 5 | **Connectivity** | Which vertices get joined? | primal, per-cell case table | dual, quad per crossed edge | dual, quad per crossed edge |
| 6 | **Feature preservation** | Are sharp edges kept? | no | no | yes (via Hermite data) |
| 7 | **Adaptivity** | Can resolution vary in space? | no (octree MC: yes) | no | yes (octree DC) |
| 8 | **Incrementality** | What happens after an edit? | full re-mesh | full re-mesh | full re-mesh |
| 9 | **Execution model** | How does it map to hardware? | per-cell, variable output | per-cell then per-edge | per-cell then per-edge |
| 10 | **Guarantees** | What is provably true of the output? | manifold (with correct table) | none | neither manifold nor intersection-free |

**Read that table for the empty cells.** Axis 8 is "full re-mesh" for every published algorithm. Axis
10 has a "none" and a "neither." Axis 3 doesn't even apply to two of the three. Those are gaps.

---

# The vocabulary, by axis

---

## Axis 1 — Domain decomposition

**What it is:** the partition of space you march over.

| Word | What it means | Why you want it |
|---|---|---|
| **Cell complex** / **CW complex** | A space built by gluing cells (points, edges, faces, volumes) along their boundaries | The general framing. "Grid" is a special case; saying *cell complex* opens the topology literature |
| **Simplicial complex** | A cell complex made only of simplices (triangles, tetrahedra) with a gluing rule | Tet meshes are these. Half of computational topology assumes them |
| **Conforming** / **non-conforming** | Adjacent cells meet full-face-to-full-face, no partial overlaps | This is what "crack-free" *means* formally. Search this word, not "crack" |
| **Hanging node** / **T-junction** / **T-vertex** | A vertex on one cell's face that is mid-edge for its neighbour | The exact defect Transvoxel exists to fix. FEM has decades on it |
| **2:1 balance** / **restricted octree** | Neighbouring octree cells differ by at most one level | The standard precondition that makes transition cells finite in number |
| **Kuhn** / **Freudenthal triangulation** | The 6-tet decomposition of a cube from the monotone corner paths | The right name for Marching Tetrahedra's decomposition. Predates Doi & Koide by 50 years |
| **BCC / FCC lattice** | Body- and face-centred cubic point lattices | More isotropic than the cubic grid — fewer orientation artefacts for the same sample count |
| **Delone set** / **ε-net** | Points that are well-separated and cover densely | The sampling-theory framing of "resolution" |
| **Voronoi diagram** / **power diagram** | Partition by nearest site (weighted, for power) | The dual of Delaunay. Dual methods are secretly doing this |
| **Restricted Delaunay** | Delaunay triangulation intersected with a surface | The construction behind provably-correct surface reconstruction |
| **Red-green refinement**, **newest vertex bisection**, **longest edge bisection** | Three named strategies for refining a simplicial mesh conformingly | FEM's answer to adaptive subdivision without cracks. Each has convergence theorems |

**Field that owns it:** computational geometry; adaptive finite element methods.
**Search terms that actually retrieve:** `conforming adaptive refinement`, `2:1 balanced octree`,
`newest vertex bisection convergence`, `body-centered cubic isosurface`.

**Transfer opportunity [mine, unverified]:** the game literature reinvented conforming refinement as
"Transvoxel transition cells." FEM solved the general problem with *newest vertex bisection*, which has
a proven bound on the number of refinements propagated per edit. That bound is exactly the
"edit-proportional repair" property the opportunities doc says field-derived LOD needs.

---

## Axis 2 — Sign inference and the underlying field

**What it is:** how a continuous field becomes a discrete inside/outside decision.

| Word | What it means | Why you want it |
|---|---|---|
| **Level set** / **sublevel set** | `{x : f(x) = 0}` and `{x : f(x) ≤ 0}` | The formal name for what you're extracting |
| **Nodal domain** | A connected component where a function has constant sign | The thing your sign grid is sampling |
| **Signed distance function (SDF)** vs **implicit function** | SDF satisfies `‖∇f‖ = 1`; an implicit merely has the right zero set | Enormous practical difference. Most "SDFs" in games are not SDFs, and error bounds that assume `‖∇f‖=1` silently break |
| **Eikonal equation** | `‖∇f‖ = 1` — the PDE an SDF satisfies | Search this to find redistancing / fast marching / fast sweeping |
| **Lipschitz bound** | `\|f(x) − f(y)\| ≤ L‖x−y‖` | The weakest useful assumption. Sphere tracing needs it; so does any "this cell is empty" proof |
| **Inclusion function** / **interval extension** | A function returning a guaranteed range over a box | The formal tool for "prove this cell contains no surface." Underneath `fidget` |
| **Reach** / **local feature size (LFS)** / **medial axis** | How far you can offset a surface before it self-intersects; distance to the medial axis | **The** sampling-theory quantity. Every "if you sample densely enough, the reconstruction is correct" theorem is stated in terms of LFS |
| **ε-sample** | Sample set where every surface point has a sample within `ε · LFS` | The precondition of the Amenta–Bern reconstruction guarantees |
| **Nyquist rate** / **aliasing** | Sampling below twice the highest frequency loses information irrecoverably | Why sub-voxel features are invisible to sign-based methods — and what subgrid MT sidesteps |

**Field that owns it:** surface reconstruction theory; interval analysis; PDE (level set methods).
**Search terms:** `local feature size reconstruction guarantee`, `epsilon-sample surface`,
`interval arithmetic implicit surface`, `Lipschitz implicit sphere tracing`.

**Transfer opportunity [mine, unverified]:** LFS gives a *spatially varying* sampling requirement.
Every voxel engine uses uniform resolution or camera-distance LOD. Nobody drives resolution by
estimated local feature size — which is the theoretically correct criterion and computable from the
field.

---

## Axis 3 — Ambiguity resolution

**What it is:** when corner signs are consistent with more than one topology, which do you pick?

| Word | What it means | Why you want it |
|---|---|---|
| **Trilinear interpolant** | The degree-3 (multilinear) function agreeing with 8 corner values | The *reference* surface MC is approximating. "Correct" means matching this |
| **Asymptotic decider** | Resolve a face ambiguity by the sign of the bilinear saddle value | Nielson & Hamann 1991. The face case, cheap and settled |
| **Body saddle** / **interior ambiguity** | The 3D analogue — a tunnel inside a cell | The hard case. A-002b defers it, correctly |
| **Morse theory** | Study of a function's topology via its critical points | The parent theory. *Ambiguity* is really "a critical point falls inside a cell" |
| **Critical point**, **index**, **saddle** | Where `∇f = 0`; the signature of the Hessian there | Vocabulary for classifying what kind of ambiguity you have |
| **Morse–Smale complex** | Decomposition of a domain by gradient flow between critical points | The full topological skeleton. Overkill for one cell, ideal for a chunk |
| **Persistent homology** / **persistence diagram** | Tracks topological features as a threshold sweeps, recording birth and death | Gives every tunnel and handle a **numeric significance** |
| **Persistence-based simplification** | Cancel features whose persistence is below a threshold | Principled topological denoising |
| **Topologically faithful** / **isotopic** | Output is homeomorphic (or ambient-isotopic) to the true level set | The strong correctness statements are phrased this way |

**Field that owns it:** Morse theory; computational topology / TDA.
**Search terms:** `persistent homology isosurface simplification`, `Morse-Smale complex extraction`,
`topologically correct isosurface trilinear`.

**Transfer opportunity [mine, unverified] — this is the one I'd chase first.** The interior ambiguity
problem is currently "transcribe a 730-subcase table that the literature says is partly wrong"
(Custodio et al. showed Chernyaev's interior test mis-tracks the saddle trajectory). Persistence
reframes it: **don't ask "is there a tunnel," ask "does this tunnel have persistence above ε."** Below
threshold, mesh it closed; above, mesh it open. That replaces a disputed case table with a computable
scalar, gives you a *tunable* knob a game wants, and is stable under perturbation by the stability
theorem for persistence diagrams — a real theorem, not a heuristic.

---

## Axis 4 — Vertex placement

**What it is:** given that a cell contains surface, where do you put the point?

| Word | What it means | Why you want it |
|---|---|---|
| **Hermite data** | Position **and normal** at each edge crossing | The input that makes sharp features recoverable. Without normals you cannot do it |
| **Quadratic error function (QEF)** | `Σ (nᵢ·(x − pᵢ))²` — squared distance to the tangent planes | Dual Contouring's objective |
| **Normal equations** | Solving `AᵀA x = Aᵀb` | The naive QEF solve. **Squares the condition number** — the reason f32 breaks |
| **Condition number** | Amplification factor from input error to output error | The number to instrument in E-112 |
| **Rank-revealing QR** / **truncated SVD** | Factorizations that expose and truncate near-null directions | The numerically sound solve. Truncation is how you detect "this cell is flat, don't move the vertex" |
| **Tikhonov regularization** / **ridge** | Add `λI` to stabilize an ill-posed solve | The `λ ≈ 0.01` in the audit's formula. Now you know its name and its literature |
| **Pseudoinverse** / **Moore–Penrose** | Minimum-norm least-squares solution | What truncated SVD computes |
| **Equivariance** vs **invariance** | Output rotates with input; output unchanged by rotation | Vertex placement should be *equivariant*. Most implementations aren't |
| **Moving frames** / **Cartan** | Systematic construction of equivariant quantities | The general machinery behind the closed-form rule in the audit |
| **First fundamental theorem of invariant theory** | Every invariant of vectors under `O(n)` is a function of their pairwise dot products | Why "build it from dot products and cross products" is *complete*, not a trick |
| **Second fundamental form** / **shape operator** | The curvature tensor of a surface | If you want *curvature-aware* placement instead of planar, this is the object |

**Field that owns it:** numerical linear algebra; invariant theory; differential geometry.
**Search terms:** `rank revealing QR quadratic error function`, `equivariant vertex placement`,
`Tikhonov regularization ill-conditioned least squares`, `moving frames invariant theory vectors`.

**Transfer opportunity [partly done]:** the audit's closed-form three-plane rule is a moving-frames
construction. **Unexplored:** fitting a *quadric* instead of a plane per crossing — i.e. using the
second fundamental form — which would place vertices on curved surfaces with second-order rather than
first-order accuracy. That directly attacks P-2's `O(|e|²κ)` error term.

---

## Axis 5 — Connectivity

**What it is:** which vertices get joined into faces.

| Word | What it means | Why you want it |
|---|---|---|
| **Primal** vs **dual** | Vertices on grid edges + faces per cell (MC), vs vertices in cells + faces per grid edge (SN/DC) | The single cleanest way to classify extractors |
| **Poincaré duality** | The formal correspondence between k-cells and (n−k)-cells | Why the primal/dual pairing is exact, and why `F_sn = F_mc + 2χ` is forced |
| **Star** and **link** of a vertex | All cells containing it; the boundary of that | **The** tool for checking manifoldness. A mesh is a manifold iff every vertex link is a circle |
| **Pseudomanifold** | Every edge in exactly 2 faces, but vertex links may be pinched | The precise name for what plain DC produces — and why edge-manifoldness alone isn't enough |
| **Vertex splitting** | Emit multiple vertices per cell where the link is pinched | Manifold DC's mechanism. A-010 |
| **Non-manifold vertex** | A vertex whose link is two or more circles | The defect A-015 was chasing. Edge counts will not catch it |
| **Orientability** / **consistent orientation** | All faces wound the same way relative to the surface | `inconsistently_oriented_edges` in T-001 is checking this |
| **Euler characteristic** `χ` | `V − E + F`, a topological invariant | Already load-bearing here |
| **Genus** | Number of handles; `χ = 2 − 2g − b` for orientable surfaces with `b` boundaries | |

**Field that owns it:** combinatorial / algebraic topology.
**Search terms:** `vertex link manifold condition mesh`, `pseudomanifold repair`,
`manifold dual contouring vertex splitting`.

---

## Axis 6 — Guarantees and soundness

| Word | What it means |
|---|---|
| **Watertight** | No boundary edges. Weaker than manifold — an informal term, be careful reading it |
| **Embedding** vs **immersion** | Injective (no self-intersection) vs only locally injective | **The exact vocabulary for "manifold but self-intersecting."** DC output is an immersion, not an embedding |
| **Self-intersection-free** | An embedding | |
| **Exact predicates** / **adaptive precision arithmetic** | Geometric tests guaranteed correct despite floating point (Shewchuk) | The standard fix for "my algorithm is correct but crashes on degenerate input" |
| **Simulation of Simplicity (SoS)** | Symbolic perturbation to remove degeneracies | Edelsbrunner & Mücke. Makes "what if the field is exactly zero at a corner" go away by construction |
| **Robustness** vs **stability** | Correct despite arithmetic error; output varies continuously with input | Different properties; papers conflate them |

**Field that owns it:** robust computational geometry.
**Search terms:** `Shewchuk adaptive precision predicates`, `simulation of simplicity degeneracy`,
`immersion embedding mesh self-intersection`.

**Note the framing win:** "manifold or intersection-free, pick one" (already falsified as ✗2) is, in
proper vocabulary, "is this immersion an embedding?" — and *that* question has an actual literature.

---

## Axis 7 — Adaptivity and multiresolution

| Word | What it means | Why you want it |
|---|---|---|
| **A posteriori error estimator** | Computable-after-the-fact bound on local error | FEM's answer to "where should I refine?" Games use camera distance instead |
| **Residual-based estimator** | Error estimated from how badly the solution fails the equation locally | The workhorse variety |
| **Dual weighted residual** / **goal-oriented adaptivity** | Refine to reduce error *in a specific output quantity*, not globally | Refine where the **player is looking**, with theory behind it |
| **Multiresolution analysis** / **wavelet** / **lifting scheme** | Decomposition into coarse + detail at every scale, invertible | The principled version of "mip the field." Lifting gives it in-place and integer-exact |
| **Subdivision surface** | Refinement scheme converging to a smooth limit surface | The other direction: coarse mesh + rules instead of dense samples |
| **Restriction** / **prolongation** | Coarsening and refining operators between levels | Multigrid vocabulary; exactly what LOD needs |

**Field that owns it:** adaptive FEM; multiresolution analysis; multigrid.
**Search terms:** `a posteriori error estimator adaptive refinement`, `goal oriented adaptivity dual
weighted residual`, `lifting scheme integer wavelet`.

**Transfer opportunity [mine, unverified]:** LOD selection in every voxel engine is `distance to
camera`. FEM has forty years of theory on *provably near-optimal* refinement criteria. A field-derived
LOD driven by a residual estimator would refine where the surface is geometrically complex rather than
where it happens to be close — different, and defensible.

---

## Axis 8 — Incrementality

**The most empty column in the table. Every published extractor is full-re-mesh.**

| Word | What it means | Why you want it |
|---|---|---|
| **Self-adjusting computation** / **change propagation** | Umut Acar's framework: run once, record a dependency trace, then update the output in time proportional to the change | **The formal theory of the dirty-set problem**, with theorems |
| **Dynamic dependence graph (DDG)** | The recorded trace that makes the above work | |
| **Trace stability** / **computation distance** | How much the trace changes when the input changes — the quantity the time bound is stated in | If your algorithm is *unstable* in this sense, incremental repair provably can't be cheap. This is why a METIS-partitioned DAG can't be repaired |
| **Batch-dynamic algorithm** | Handles a batch of updates at once, work-efficiently, in parallel | The modern shape. DynHAC is one |
| **Work–span (work–depth) model** | Total operations vs critical-path length | The standard cost model for parallel algorithms. Learn to state your algorithm in it |
| **Kinetic data structure** | Maintains a combinatorial structure as inputs move continuously, with certificate failures | For *animated* fields rather than edited ones |
| **Persistent data structure** | Old versions survive updates | What HashDAG uses for editable SVDAGs, and what the world-as-a-file row needs |
| **Locality of reference** / **write-once** | Each output location written once; updates touch a bounded neighbourhood | The precondition change-propagation frameworks require |

**Field that owns it:** dynamic algorithms; self-adjusting computation; parallel algorithms.
**Search terms:** `self-adjusting computation change propagation`, `batch-dynamic parallel algorithm`,
`trace stability incremental computation`, `kinetic data structure`.

**Transfer opportunity [mine, unverified] — the biggest gap in the table.** Nobody has stated
isosurface extraction as a self-adjusting computation. The pieces fit suspiciously well: the field is
the input, the mesh is the output, edits are input changes, and the extractor is already
cell-local — which is exactly the *stability* property Acar's bounds require. If it works you get
edit-proportional re-meshing with a **proof**, instead of a dirty-set heuristic.

---

## Axis 9 — Execution model

| Word | What it means | Why you want it |
|---|---|---|
| **Arithmetic intensity** / **roofline model** | FLOPs per byte moved; the plot showing whether you're compute- or bandwidth-bound | Settles "is this worth optimizing" before you optimize |
| **Stream compaction** / **prefix sum (scan)** | Turning a sparse flagged array into a dense one | The primitive underneath every GPU MC. Variable output per cell is *the* difficulty |
| **Warp / wavefront divergence** | Threads in a SIMD group taking different branches | Why MC's variable case-per-cell output "makes for asymmetrical code" (Media Molecule) |
| **Occupancy** | Fraction of hardware thread slots in use | |
| **Coalescing** / **bank conflict** | Memory access patterns that are fast vs serialized | The mechanism behind the 14.5× LDS-vs-atomics result |
| **Persistent threads** | Long-lived workers pulling from a queue instead of one-thread-per-item | The pattern for irregular work before work graphs existed |
| **Task / mesh / amplification shader** | Pipeline stages that generate geometry on-GPU without a vertex buffer | The 23.4× result |
| **Indirect dispatch** | GPU decides its own launch dimensions | Lets the GPU size the next stage without a CPU round trip |
| **Work graphs** | GPU-side dynamic task spawning | Two-sided: big wins on irregular work, 2.8–3.4× *losses* on classification |

**Field that owns it:** parallel algorithms; GPU architecture.
**Search terms:** `stream compaction marching cubes GPU`, `mesh shader isosurface`,
`roofline model memory bound kernel`, `persistent threads irregular parallelism`.

---

## Axis 10 — Consistency across chunks

**What it is:** independent pieces agreeing where they meet. The seam problem, properly named.

| Word | What it means | Why you want it |
|---|---|---|
| **Sheaf** / **presheaf** | Data assigned to open sets, with restriction maps | The formal object for "each chunk has data" |
| **Gluing axiom** / **sheaf condition** | Locally-agreeing sections glue to a unique global section | **This is literally "chunks that agree on overlaps form one mesh."** The word exists |
| **Čech cohomology** / **H¹** | The obstruction to gluing local data globally | If `H¹ ≠ 0`, no amount of local fixing produces global consistency. A *theorem*, not advice |
| **Cellular sheaf** / **sheaf Laplacian** | Discrete, computable versions on a cell complex | The tractable form. Active area in applied topology |
| **Join-semilattice** / **commutative monoid** / **idempotent** | Order/algebra structures under which merge order stops mattering | The exact algebra behind conflict-free merging |
| **CRDT** (**CvRDT** / **CmRDT**) | Data types that converge without coordination | The engineering name for the above |
| **CALM theorem** | A program has a coordination-free distributed implementation **iff** it is monotone | An equivalence, not a heuristic. Tells you when dig-only gets multiplayer free |
| **Monotonicity** | Output only grows as input grows | The property you must preserve to stay coordination-free |

**Field that owns it:** applied sheaf theory; distributed systems theory.
**Search terms:** `cellular sheaf cohomology consistency`, `sheaf Laplacian data fusion`,
`CALM theorem monotonicity coordination-free`, `CRDT semilattice convergence`.

---

# How to use this to actually find things

**The move that works:** take your problem, find its word in the table above, then search *that word
plus the source field*, not plus "mesh" or "voxel." Searching `crack-free LOD voxel` returns game blog
posts. Searching `conforming adaptive refinement hanging node` returns forty years of FEM theory with
proofs.

**Three signals for triaging a paper fast:**

- **Does it state a theorem or report a measurement?** Both are fine; conflating them isn't. A paper
  claiming "topologically correct" should say *homeomorphic* or *isotopic* to something specific.
- **What is the error measured against?** "Compared to the trilinear interpolant" and "compared to the
  true surface" are different claims. Papers slide between them.
- **What are the sampling assumptions?** If a guarantee doesn't mention **LFS**, **reach**, or a
  **Lipschitz bound**, it isn't a guarantee about arbitrary inputs.

**Words that signal a paper is worth your time here:** *equivariant*, *a posteriori*, *conforming*,
*isotopic*, *work-efficient*, *batch-dynamic*, *stability* (in the trace sense), *local feature size*.

**Words that signal folklore:** *watertight* (undefined), *high quality*, *robust* (unqualified),
*fast* (no denominator), *artifact-free*.

---

# The five transfers I'd rank first

All tier **F** — hypotheses, not findings. Each names the field to raid and what would falsify it.

| # | Transfer | From | The move | Falsified by |
|---|---|---|---|---|
| 1 | **Isosurface extraction as a self-adjusting computation** | Self-adjusting computation (Acar) | State the extractor with a dynamic dependence graph; get edit-proportional re-meshing with a proven bound instead of a dirty-set heuristic | Showing the trace is unstable under a local field edit — i.e. computation distance is not `O(edit size)` |
| 2 | **Persistence-thresholded ambiguity resolution** | Persistent homology | Replace the disputed 730-subcase interior test with "mesh the tunnel iff its persistence > ε" | Persistence of a cell-local feature not being computable in `O(1)`, or the threshold not being stable |
| 3 | **LFS-driven adaptive resolution** | Surface reconstruction theory | Drive refinement by estimated local feature size — the theoretically correct criterion — instead of camera distance | LFS estimation costing more than the meshing it saves |
| 4 | **A posteriori error estimators for LOD** | Adaptive FEM | Refine where a residual estimator says error is large; goal-oriented variants refine where the player looks | The estimator not correlating with measured Hausdorff error on our fields |
| 5 | **Second-order vertex placement** | Differential geometry | Fit the second fundamental form per cell instead of tangent planes; attack the `O(\|e\|²κ)` error term directly | Curvature estimation from Hermite data being too noisy at game resolutions |

**Why #1 first:** it is the only one that addresses the axis where *every* published algorithm is
identical (axis 8, "full re-mesh"). Axes where everyone already differs are crowded. The axis where
everyone is the same is where the unoccupied space is.
