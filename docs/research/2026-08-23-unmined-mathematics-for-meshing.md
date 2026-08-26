# Areas of mathematics not yet brought to bear, read against the new findings

**Date:** 2026-08-23 · **Repo state:** `61c6201`, 431 entries
**Corpus caveat:** home-still was **not reachable** from this session — the MCP server is not connected,
the proxied local Filesystem MCP fails on every call with a JSON Schema dialect error (draft-07 against
a 2020-12-only validator), and the folder grant was declined. Everything below is web-primary. A
home-still sweep would change the confidence on the *negative* claims — "nobody has done this" — much
more than on the positive ones.

---

## What the ledger has already mined, so the gaps are visible

Combinatorial topology (Euler characteristic, manifoldness, normal surface theory, well-composedness,
cubical critical configurations). Differential geometry (normal cycles, curvature measures, medial
axis, λ-medial). Numerical linear algebra (QEF, Tikhonov, nested-dissection Cholesky, heat operator).
IEEE floating point (equivariance, reproducible ordering, bit-identity). Interval and Lipschitz bounds.
Dynamic graph algorithms (union-find, batch-dynamic connectivity, HDT, replacement search).
Self-adjusting computation. Sampling theory. Rigid-block equilibrium. Sabine acoustics. Dissolution
kinetics. Convex decomposition.

**Four things are conspicuously absent, and three of them have a hook into a live finding.**

---

## 1. Real algebraic geometry of *multi-affine* hypersurfaces — the strongest find

**The observation that makes it relevant.** The trilinear interpolant is not a generic cubic. It is
**multi-affine**: degree 1 in each variable separately, total degree 3. That is a special class with its
own literature, and the visualisation community does not cite it.

**The theorem.** Basu & Perrucci, *Topology of real multi-affine hypersurfaces and a homological
stability property*, **Advances in Mathematics 420 (2023), 108982**, `10.1016/j.aim.2023.108982`,
arXiv:2204.01595. Their **Theorem 2**: the number of semi-algebraically connected components of a real
hypersurface defined by a multi-affine polynomial of total degree `d` is bounded by **2^(d−1)**,
**independent of the ambient dimension `n`**. Compare the classical Petrovskiĭ–Oleĭnik–Thom–Milnor
bound of `d(2d−1)^(n−1)`.

For `d = 3` — the trilinear — that is **at most 4 components**.

**Why that number should make you sit up.** Four is the maximum number of contour components MC33
allows in a cell, and it is the `MAX_CENTROIDS`-adjacent budget A-015 and M-217 argued about
case-by-case. If the box-restricted version of Theorem 2 holds, **the case taxonomy becomes a
corollary of a theorem instead of a table you derived and had to guard with `validate_table()`.** That
is the difference between rule 5 forbidding you to guess and rule 5 becoming unnecessary.

**Two caveats that must not be skipped.** Theorem 2 is for the hypersurface in ℝⁿ, **not restricted to
a box** — intersecting with the unit cube can split components, so the bound does not transfer to a
cell without further work. And their **Theorem 3** shows that for multi-affine hypersurfaces
intersected with a ball, the sum of Betti numbers `b₀…b₅` grows *exponentially* in `n`, so the clean
`b₀` bound emphatically does not generalise upward. The `d = 3` coincidence with MC33's four is
**unverified** — I have not checked the correspondence, only noticed it.

**The experiment shape.** Establish or refute the box-restricted `b₀ ≤ 4` by exhaustive enumeration
over the 256 sign configurations × ambiguity masks — machinery A-002f/A-020 already has — and report
whether any configuration reaches four. If the bound holds and is tight, the ledger gains a *derived*
budget where it currently has a searched one.

---

## 2. Pseudo-Boolean functions and the Lovász extension — the framing that names two of your extractors

This came out of the same probe and is worth more than it first looks.

**The identification.** The **multilinear extension** of a pseudo-Boolean function is exactly the
trilinear interpolant of eight corner values. The **Lovász extension** is the piecewise-linear
extension, affine on each simplex of the Freudenthal/Kuhn order-triangulation of the cube — **which is
precisely the 6-tetrahedron decomposition P-3 verified and `marching_tetrahedra` uses.** So:

> **Marching Tetrahedra contours the Lovász extension. Marching Cubes contours the multilinear
> extension.** Both extractors are in this crate, and the mathematical literature already has a Morse
> theory for one of them.

Jost et al., *Discrete-to-Continuous Extensions: Lovász Extension and Morse Theory*, **Discrete &
Computational Geometry (2022)**, `10.1007/s00454-022-00461-1`, arXiv:2003.06021, gives discrete-to-
continuous Morse theory through the Lovász extension. **That is a rigorous account of why Marching
Tetrahedra has no ambiguity** — which the ledger currently states as an observation (M-67, M-81, P-4's
16-entry table) rather than as a consequence.

**A second identity, cleaner than the algebra it replaces.** The standard multilinear-extension fact

```
∂²F/∂xᵢ∂xⱼ  =  f(S ∪ {i,j}) − f(S ∪ {i}) − f(S ∪ {j}) + f(S)
```

says the mixed second derivative **is** the discrete cross-difference, whose sign is sub- or
supermodularity in the pair `(i, j)`. On a cell face that cross-difference is exactly the
`(A + C) − (B + D)` denominator M-204 spent two parentheses fixing. So **the asymptotic decider is a
submodularity test on the face's two-variable pseudo-Boolean function** — Nielson–Hamann restated in a
basis-free language where the grouping M-204 needed is forced by the definition rather than by a
fixture search. That is not new mathematics; it is a framing that would have made M-204's bug
unrepresentable.

**One clean negative, so nobody re-searches it.** The pseudo-Boolean *optimisation* literature — Boros
& Hammer `10.1016/S0166-218X(01)00341-9`, roof duality, submodular minimisation — is **useless here**.
A multi-affine function attains its optimum at a vertex, so optimising the multilinear extension
reduces to Boolean optimisation and the literature contains nothing about level sets, components,
tunnels or saddles.

---

## 3. Min-plus / polyhedral geometry of the edit tape — and the certificate it hands you for free

**The starting observation.** `Add = min`, `Subtract = max ∘ negate`. A hard-op tape is a min-plus
expression, and M-36/M-37/M-38's whole commutativity result is a statement about that algebra: same-kind
ops commute bit-exactly because min and max *select* an operand; mixed ops do not.

**The confirmation.** For a pure-`Add` tape of spheres, the seam — the locus where the argmin switches
— **is the additively-weighted (Apollonius) Voronoi diagram of the (centre, radius) sites**, by
definition and not by analogy: `δ(x, Pᵢ) = ‖x − cᵢ‖ − wᵢ` is literally the sphere SDF. Each seam sheet
is one sheet of a hyperboloid of revolution with foci `cᵢ, cⱼ` and `2a = |rᵢ − rⱼ|`, degenerating to the
perpendicular bisector when the radii match. Verified numerically at 400k samples. It is computable:
CGAL `Apollonius_graph_2` in 2D, and **Plateau-Holleville et al., *In Search of Empty Spheres: 3D
Apollonius Diagrams on GPU*, ACM TOG 44(4) 2025, `10.1145/3730868`** — 100k sites in 8.1 s on a 4090.

**The refutation, and it is the genuinely unexplored object.** The moment the tape contains a
`Subtract`, the bisector between an Add-sphere `i` and a Subtract-sphere `j` is

```
|x − cᵢ| + |x − cⱼ| = rᵢ + rⱼ
```

— a **prolate spheroid** with foci `cᵢ, cⱼ`, non-empty exactly when the spheres intersect. Verified
numerically. **Sum-of-distances bisectors are produced by no standard Voronoi diagram.** The
minimization diagram of `{+dᵢ} ∪ {−dⱼ}` over ℝ³ has no name in the literature, no complexity bound
better than Sharir's generic `O(n^{3+ε})` (`10.1007/BF02574384`), and no implementation. Lie sphere
geometry (arXiv:2408.09279, 2024) is the only framework that plausibly hosts it.

**That is also the structural explanation for M-37.** Mixed add/subtract gives 11 distinct results
across 40,320 orderings not because of arithmetic but because the two operand families live on two
different bisector geometries.

**Where "tropical" earns its keep, exactly once.** Substitute the **power** pseudo-distance
`|x − c|² − r²`. Then

```
min_i (|x − cᵢ|² − rᵢ²)  =  |x|²  +  min_i (−2cᵢ·x + |cᵢ|² − rᵢ²)
```

— `|x|²` plus a **min of affine functions**, which is a genuine min-plus polynomial with a polyhedral
corner locus (the power diagram), a Newton polytope, and the whole tropical toolkit. Verified to
2.8e-14. **And it comes with an algorithm this crate already owns:** Felzenszwalb & Huttenlocher,
`10.4086/toc.2012.v008a019` — the same paper S-001/M-251 implemented for the exact distance transform —
computes the lower envelope of `n` equal-shape parabolas in **O(n) amortised**. Along any grid line,
`min_i(|x − cᵢ|² − rᵢ²)` is exactly that envelope, so per-scanline argmin costs `O(n + gridsize)`
instead of `O(n · gridsize)`. **You have the machinery and have not pointed it at the tape.**

The price is measured and must be stated: `argmin(power)` agrees with `argmin(true SDF)` on only **84%**
of samples over 8 spheres. It is a sound pruner and accelerator, **not** a substitute for the SDF
argmin.

**The item with the best ratio of value to effort in this whole memo, and it is not a paper.** Carry the
**margin** `f₍₂₎ − f₍₁₎` — the gap between the two smallest primitive values — through the fold you
already run. Then, with `L` the Lipschitz bound and `h` the stencil radius:

> **`f₍₂₎ − f₍₁₎ > 2hL`  ⟹  the stencil cannot straddle a seam.**

One extra comparison. Sound, one-sided, and it is exactly the certificate that would have caught
P-47's single 4.4° vertex *before* it produced a normal — the vertex that carried the entire reported
mean of 7.6e-5° while the bulk sat at 1.9e-8°. Every ingredient is published: Rossignac & Requicha's
**active zones** `10.1145/49155.51123` (the 1986 successor to Tilove that P-39 should cite alongside
1984), Barbier's Lipschitz bounds `10.1111/cgf.70057`, Sharp & Jacobson's range analysis
`10.1145/3528223.3530155`, and Scholtes' **essentially active index set** (*Introduction to Piecewise
Differentiable Equations*, `10.1007/978-1-4614-4340-7`) — which is the analytic name for exactly the
brush set P-39's pruning computes. **Nobody has assembled them into a seam-aware mesher.**

**And the smoothing parameter turns out to be a length.** `smooth_min(k)` is the log-semiring; hard
`min` is its Maslov tropical limit, with `|smin_k − min| ≤ k·ln n` concentrated in an `O(k)` shell
around the seam. Measured: at `k = 0.1`, max deviation 0.135 against the bound 0.208, and outside a
`10k`-wide margin shell the deviation is **1.3e-5**. So `k` **is** a computable seam-neighbourhood
width, and the margin `f₍₂₎ − f₍₁₎` is the coordinate that measures it. That is also why M-38 measured
40,317 distinct orderings: smoothing destroys the operand-selection property that makes min and max
exact.

---

## 4. Certified subdivision instead of hand-derived algebra — the engineering answer to a stated complaint

The A-002 series' recurring cost is not conceptual, it is transcription: M-207's textbook quadratic
losing a root, M-219's typo in the reference implementation's max, M-221's `0 × NaN`, ✗22's
unrecoverable counterexample, M-228's undefined three-step case, M-231's `[9,3]` cells.

Every partial derivative of a trilinear is **multi-affine in the remaining two variables**, so its
exact range over any sub-box is min/max of that sub-box's corner values — free, and *exact* rather than
merely tight. That makes a **Plantinga–Vegter-style certified subdivision** (`10.1145/1057432.1057465`,
which T-015 already implements for isotopy) unusually cheap for solving `∇f = 0`, i.e. for locating
body saddles. **I found nobody doing PV on the trilinear-grid Marching Cubes problem.** It is a
correctness-engineering win against your actual complaint, not a mathematical novelty.

The general machinery, if you want it: Mourrain & Pavone, *Subdivision methods for solving polynomial
equations*, `10.1016/j.jsc.2008.04.016`; Sherbrooke & Patrikalakis, `10.1016/0167-8396(93)90019-Y`;
Alberti, Comte & Mourrain 2005, which is the closest existing thing — Bernstein form plus a multivariate
Descartes rule as a cell-regularity criterion, for general algebraic surfaces.

---

## The kill: Bernstein–Bézier form is novel *and* dead for the trilinear

I chased this hard because it looked like the obvious systematic replacement for the ad-hoc algebra,
and it does not work. Recording it so it is not re-attempted.

For a degree-(1,1,1) tensor-product function the Bernstein basis is `B₀(t) = 1−t`, `B₁(t) = t`, so the
**eight Bernstein coefficients are identically the eight corner values**. De Casteljau subdivision of a
multi-affine function produces sub-cell coefficients that are exactly `f` at the sub-cell corners —
subdivision **is** resampling. Consequences, each fatal:

- The convex-hull enclosure is not "tight after subdivision", it is **exact at every level** — Garloff's
  vertex condition holds unconditionally for multi-affine functions, so overestimation is zero. Which
  also makes it uninformative: it is `[min, max]` of the corner values, i.e. the 1992 branch-on-need
  octree test (`10.1145/130881.130882`).
- The **sign pattern of the coefficients is the Marching Cubes case index** — precisely the object
  already known to be insufficient to resolve the ambiguity. Variation-diminishing gives `V(b) ≥ #roots`
  with equality mod 2, which is a restatement of the ambiguity, not a resolution of it.
- The sharpness theorems are vacuous here. The canonical pair is *linear in 1/degree under degree
  elevation* (Rivlin 1970) and *quadratic in box width under subdivision* (Stahl 1995) — both describe
  overestimation that is already zero.

**Verdict: novel because nobody published it, dead because it provably cannot work.** The one part that
survives — Bernstein/Descartes root isolation along a *general* line through a cell, where the trilinear
restricts to a cubic in `t` — is already standard practice (Loop & Blinn `10.1145/1141911.1141939`;
Reimers & Seland, Eurographics 2008).

**Revisit the entire hypothesis if the crate ever adopts a tricubic reconstruction filter.** At degree ≥ 2
per variable every part of it becomes true, useful, and as far as I can tell unpublished for isosurface
topology.

---

## Areas I judged worth naming but did not probe

- **Discrete Morse theory on the cubical grid** (Forman; Robins–Wood–Sheppard for cubical complexes,
  linear time). A combinatorial gradient vector field computable in one pass, giving critical cells
  directly. This is what P-52's monotone-edge gate is groping toward, and it is the natural way to
  assert topology on `gyroid` and `fbm_terrain` where χ is unavailable. Item 2's Lovász–Morse result is
  the tetrahedral half of the same picture.
- **Box-constrained QP for vertex placement.** The cell clamp is a *projection after an unconstrained
  solve*, which is not the constrained minimiser. A 3-variable box-constrained QP by active set
  terminates in ≤ 3 iterations and gives the true minimiser of the QEF inside the cell. M-28 and M-30
  measured what the clamp *does*; nobody has measured the gap between "solve then clamp" and "solve
  subject to the box". Cheap, classical (KKT), zero topology risk, and the `VertexRule` seam is already
  there.
- **Optimal rectangle partition** as the baseline greedy meshing is measured against. M-56 measured
  greedy at 1.70×–256× over face culling but never against the *optimum*, which for a rectilinear region
  is a classical polynomial-time problem (maximum independent set on a bipartite intersection graph —
  Lipski, Ohtsuki). That would turn "greedy is field-dependent" into "greedy is within X% of optimal".
- **Lie sphere geometry** (arXiv:2408.09279) as the frame for §3's unnamed mixed-sign diagram. It turns
  second-order sphere conditions into first-order linear ones in a higher-dimensional space, reducing
  the diagram to a convex hull — the only route I found by which a non-affine min becomes genuinely
  tropical.

---

## Ranked, with the hook each has into a live finding

| | Area | Hook | Cost |
|---|---|---|---|
| 1 | **Margin `f₍₂₎ − f₍₁₎` + seam certificate** | P-47's single 4.4° vertex; P-39's pruning; M-37's non-commutativity | One comparison in the existing fold |
| 2 | **Multi-affine algebraic geometry** (Basu & Perrucci) | The case-table budget A-015/M-217 derived by search | Enumeration over machinery A-020 already has |
| 3 | **Pseudo-Boolean / Lovász framing** | M-204's grouped denominator; MT's ambiguity-freedom (M-67, M-81) | Reading; then a doc rewrite, not code |
| 4 | **Power-diagram substitution + F&H lower envelope** | S-001's own algorithm, pointed at the tape instead of the grid | Moderate — and it is a pruner, not a replacement |
| 5 | **PV certified subdivision for body saddles** | The A-002 transcription-defect series | Real work, but T-015 already has PV |
| 6 | **Box-constrained QEF** | M-28/M-30's clamp, measured but never compared to the constrained optimum | Small — the seam exists |
| 7 | **Discrete Morse on the cubical grid** | P-52; the fields where χ cannot be asserted | Moderate |
| 8 | **Optimal rectangle partition** | M-56's unbounded greedy ratio | Small, and it is a bound not a feature |
| — | ~~Bernstein–Bézier form~~ | — | **Dead. Do not attempt.** |

**What a home-still sweep would change.** Chiefly the confidence on the negatives. "Tropical geometry
has never been applied to CSG" and "nobody has run Plantinga–Vegter on the trilinear grid" rest on
roughly forty web searches, not on a corpus sweep, and this project's own record (✗4, V-29, V-32, V-43)
is that absence claims made from search results are unfounded three times out of four.
