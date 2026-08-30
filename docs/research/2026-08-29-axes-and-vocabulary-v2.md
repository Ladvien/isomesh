# Meshing decomposed into axes, v2 — the four axes v1 missed, and the words for them

**Date:** 2026-08-29
**Supersedes nothing in** `docs/research/2026-08-12-axes-and-vocabulary.md` — that document's ten axes
and their vocabulary are still correct and still the right first read. This one does three things it
did not: it names **four axes that were missing**, it upgrades the vocabulary on **three axes where
the ledger has since outgrown v1's words**, and it records one identity that was **proved during this
sweep** rather than found in a paper.

**How to read the tiers.** Vocabulary is stable and safe. Claims marked **[proved here]** were checked
symbolically in this session and the check is reproducible. Claims marked **[verified]** were traced to
a primary source and the settling sentence is quoted in
`docs/research/2026-08-29-fifty-experiments-from-unmined-mathematics.md`. Everything marked **[F]** is
a hypothesis in `FINDINGS.md` terms — no better than an opinion until someone measures it.

---

## 0. What v1 got right, and the shape of what it missed

v1's central move was correct and has held up for seventeen days: *an extractor is a pipeline of
roughly independent decisions, the named algorithms are popular combinations, and the empty regions of
the product space are where novelty lives.* Everything below is built on that.

But every one of v1's ten axes is a **mechanism** axis. Each asks *how does the algorithm do the thing*.
Not one asks *how well could the thing possibly be done* — and that turns out to be where the
mathematics with the deepest literature actually lives. Three consequences, each of which cost this
project something:

- **v1 has no vocabulary for element shape.** Axis 1 is about how space is partitioned and Axis 7 about
  how finely — but *never* about the shape of what comes out. Every cell in this crate is a cube and
  every triangle is as isotropic as that cube allows. The entire forty-year literature on **anisotropic
  mesh adaptation** answers a question v1 cannot phrase.
- **v1 has no vocabulary for reconstruction quality.** Axis 2 says "SDF vs implicit" and stops. That
  the trilinear is a *reconstruction filter* with a measurable **approximation order**, and that the
  order is a design choice, is invisible in v1's framing — which is why `✗42`'s shifted-linear
  registration reached for a knot shift rather than for the order.
- **v1 has no denominator.** Every measured number in `FINDINGS.md` is a numerator. M-12's fitted `h²`
  constant, `F-007`/`M-250`'s 13–15%, `✗42`'s ratio — each answers *did this help?* and none can answer
  *how much is left?* There is a mathematics whose entire purpose is to supply that denominator, and
  this project has never used a word of it.

So: **four new axes**, numbered on from v1's ten.

---

# Part 1 — the four missing axes

---

## Axis 11 — Reconstruction filter and approximation order

**What it is:** the trilinear interpolant is not a fact about the world. It is a *choice of
reconstruction filter*, and the field of filter design has a precise theory of what that choice costs.

| Word | What it means | Why you want it |
|---|---|---|
| **Shift-invariant space** | The span of integer translates of one generator `φ` | The formal home of "reconstruct a continuous field from samples on a lattice". Trilinear is `φ = ` tensor-product hat |
| **Approximation order** `L` | Error falls as `O(h^L)` as the sample spacing `h → 0` | M-12 measured `L = 2` for this crate. **`L` is a property of the filter, not of the problem** — that is the sentence v1 could not say |
| **Strang–Fix conditions** | `φ` has approximation order `L` **iff** its Fourier transform has an `L`-fold zero at every non-zero lattice point | The theorem that tells you `L` before you measure it. The right search term, and it retrieves a real literature |
| **Asymptotic constant** | The multiplier `C` in `error ≈ C·h^L·‖f^(L)‖` | **This is the number M-12 fitted.** Blu & Unser give a closed form for it, so the fit has a prediction to be checked against |
| **Quasi-interpolation** | Achieving full approximation order with a **compactly supported** prefilter instead of a recursive infinite-support one | The exact escape from `✗42`'s locality objection — a compact prefilter cannot break chunking by construction |
| **Interpolating vs approximating** | The filter reproduces the samples exactly, vs merely gets close | Interpolation is a *constraint*, and dropping it buys accuracy. Most graphics code assumes it without noticing |
| **Marsden's identity** | The polynomial-reproduction identity that makes order arguments go through | The workhorse lemma; seeing it cited is a signal a paper is doing real approximation theory |
| **Box spline** | The multivariate generalisation of the B-spline, defined by a direction matrix | The natural filter family on non-cubic lattices — the bridge from Axis 11 to Axis 12 |
| **Order vs smoothness** | `C^k` continuity and `O(h^L)` accuracy are *different* properties | Papers conflate them constantly. A smoother filter is not automatically a more accurate one |

**Field that owns it:** approximation theory; signal processing (the Unser/Blu school at EPFL).
**Search terms:** `Strang-Fix conditions approximation order`, `quasi-interpolant compact support
polynomial reproduction`, `asymptotic constant interpolation error Fourier`, `box spline BCC lattice`.

**What this buys immediately.** `✗42` registered shifted linear interpolation and got a nuanced
result: the recursive prefilter's non-locality turned out to be *bounded* (truncation at `k ≥ 10`
moves the root under `7.15e-7` cells), but the reconstruction gain mapped to root-position gain only as
a lottery over the crossing position. Axis 11's vocabulary says why: **a knot shift and an order
increase are different mechanisms**, and `✗42` tested the first. `h² → h⁴` is the untested lever, and
it is a different filter rather than the same filter moved.

---

## Axis 12 — Sampling lattice and quantization efficiency

**What it is:** v1 mentioned BCC and FCC lattices in one table row of Axis 1 and moved on. There is a
quantitative theory of how good a lattice is at sampling, and it produces exact numbers.

| Word | What it means | Why you want it |
|---|---|---|
| **Normalized second moment** `G` | `G = (1/d)·(∫_V ‖x‖² dx) / V^(1+2/d)` over the Voronoi cell — dimensionless, scale-free | **The** figure of merit for a sampling lattice. Lower is better |
| **Voronoi (quantization) cell** | The set of points closer to one lattice site than any other | Cube for `Z³`, truncated octahedron for BCC, rhombic dodecahedron for FCC |
| **`A₃*` / `D₃*`** | The body-centred cubic lattice, in root-system notation | **[verified]** `G(A₃*) = 19/(192·∛2) = 0.078543281`, against `G(Z³) = 1/12 = 0.0833333` — a **5.75% MSE reduction, 0.257 dB** |
| **`D₃` / `A₃`** | The face-centred cubic lattice | **[verified]** `G(D₃) = 2^(−11/3) = 0.078745066` — *worse* than BCC, and the BCC-vs-FCC gap is 0.011 dB, about 4% of the size of the cubic-vs-BCC gap |
| **Dual lattice** | `Λ* = {y : ⟨x,y⟩ ∈ ℤ ∀x ∈ Λ}` | BCC and FCC are each other's duals; the *quantization*-optimal and *packing*-optimal lattices in 3D are not the same one |
| **Covering radius / packing radius** | Largest hole; largest non-overlapping ball | Different optimality criteria pick different lattices. Say which one you mean |
| **Gersho's conjecture** | The optimal quantizer's cells tile space with one polytope — **unproven for `d ≥ 3`** | The gap between "optimal among lattices" and "optimal, full stop" |

**Field that owns it:** lattice theory; information theory (Conway & Sloane; Agrell).
**Search terms:** `normalized second moment lattice quantizer`, `optimal lattice quantizer three
dimensions`, `body-centred cubic sampling reconstruction`, `box splines BCC lattice`.

**The precision that matters, stated once so nobody overclaims it.** Barnes & Sloane (1983,
`10.1137/0604005`) prove `A₃*` is optimal **among three-dimensional lattices** — codebooks closed under
addition, cells congruent translates of one polytope. That is *not* optimality over all 3-D vector
quantizers, where the codebook need not be a lattice. The general bounds bracketing that larger class
remain unproven. **"BCC is proved optimal in 3D" is true only with the word *lattice* in it.**

---

## Axis 13 — Element shape and anisotropy

**What it is:** the axis that says a triangle has a *shape*, that the shape is a free variable, and
that there is an optimal one determined by the surface.

| Word | What it means | Why you want it |
|---|---|---|
| **Metric tensor field** `M(x)` | A symmetric positive-definite matrix at every point, defining local lengths | The object that encodes "how big and which way stretched should an element be, here" |
| **Unit mesh** | A mesh whose every edge has length 1 **as measured in `M`** | The formal statement of "conforming to a metric". Mesh generation becomes: build a unit mesh for this `M` |
| **Continuous mesh** | The metric field `M` treated as the continuous limit of a discrete mesh | Loseille & Alauzet's framing. It makes "optimal mesh" a *calculus of variations* problem instead of a heuristic |
| **Complexity** `C(M) = ∫√det M` | The continuous stand-in for vertex count | **[verified]** The constraint you optimise under. It is what makes "at fixed triangle budget" a well-posed clause |
| **Optimal `L^p` metric** | **[verified]** `M_Lp = D_Lp · det(\|H_u\|)^(−1/(2p+d)) · \|H_u\|` | Two factors doing two jobs: `det(\|H\|)^(−1/(2p+d))` is the **density**, `\|H\|` is the **orientation and aspect ratio** |
| **`\|H\|`** | The Hessian with its eigenvalues replaced by their absolute values | The metric wants a positive-definite object; the Hessian is not one. This is the standard repair |
| **Log-Euclidean interpolation** | Interpolate metrics as `exp((1−t)log M₁ + t log M₂)` | Component-wise averaging of SPD matrices is not intrinsic and can swell determinants. **This is a chunk-seam question in disguise** |
| **Aspect ratio / stretching / alignment** | Three separable properties of an element | Papers that say only "quality" are not doing this. Papers that separate the three are |
| **Anisotropic quality measure** | e.g. Vasilevskii's `\|△\|_M / \|∂△\|_M²` | The measurable version of "is this element right for this metric" |

**Field that owns it:** adaptive FEM and CFD mesh adaptation; Riemannian geometry as applied there.
**Search terms:** `continuous mesh framework interpolation error`, `metric-based anisotropic mesh
adaptation`, `optimal anisotropic meshes L^p norm interpolation`, `log-Euclidean metric interpolation SPD`.

**The honest statement of what anisotropy buys, because the seductive version is wrong.** For a
function with `W^{2,p}` regularity, uniform refinement gives `O(N^(−2/3))` in 3D, and the optimal
anisotropic mesh gives **`O(N^(−2/3))` as well — the same exponent.** The exponent `−n/d` is fixed by
the polynomial degree and the dimension and no amount of grading or stretching improves it. What
anisotropy improves is **the constant**: `‖√|det H|‖_{L^τ}` replaces `|f|_{W^{2,p}}`, and by AM–GM the
former is never larger and collapses toward zero exactly where a principal curvature vanishes —
ridges, cylinders, creases, boundary layers. A better rate appears only against functions too rough
for the uniform estimate's hypothesis. **So the pitch is "the same `h²` law with a smaller fitted
constant, concentrated on features with one flat direction" — which is precisely the constant M-12
fits and currently has no way to move.**

---

## Axis 14 — Optimality and lower bounds

**What it is:** the denominator. Given `n` field evaluations, what is the *best any algorithm could
do*, and how far from it are we?

| Word | What it means | Why you want it |
|---|---|---|
| **`n`-th minimal error** | The smallest worst-case error achievable by *any* algorithm using `n` pieces of information | The denominator. "We are within `X` of it" is a sentence this project cannot currently write |
| **Information-based complexity (IBC)** | The theory of that question — information, algorithms, and the gap | The parent field (Traub, Woźniakowski, Novak). Search this, not "how good is my mesher" |
| **Optimal recovery** | Reconstructing a function (or a functional of it) from limited samples, minimax-optimally | The framing that fits "extract a level set from `n` samples" |
| **Adaptive vs non-adaptive information** | Do you get to choose sample `k+1` after seeing samples `1..k`? | Octree refinement **is** adaptive information. That is what LOD is, in this vocabulary |
| **Kolmogorov `n`-width** `d_n(F)` | The best worst-case error approximating a class `F` from *any* `n`-dimensional linear subspace | Lower bounds for **linear** methods. Marching Cubes is nonlinear, so this bounds a neighbour, not it |
| **Metric entropy / `ε`-entropy** | `log` of the minimum number of `ε`-balls covering a class | The information-theoretic floor on how many **bits** any representation of a surface to accuracy `ε` needs |
| **`N`-term approximation** | Best approximation using `N` terms freely chosen from a dictionary | The nonlinear theory. Adaptive meshing lives here, not in the width theory |
| **Approximation class `A^s`** | The functions for which `N`-term error decays like `N^(−s)` | **The right question about a field**: which class is `gyroid` in? `fbm_terrain`? If uniform and adaptive give the same `s`, octree LOD buys nothing asymptotically on that field |

**Field that owns it:** information-based complexity; nonlinear approximation theory.
**Search terms:** `information-based complexity optimal recovery`, `n-th minimal error adaptive
information`, `nonlinear N-term approximation convergence rate`, `approximation classes adaptive finite
element`.

**One theorem, and the four hypotheses that stop it applying here.** There is a famous result — Gal &
Micchelli; Novak, *On the power of adaption*, `10.1006/jcom.1996.0015` — that **adaption improves the
worst-case error by at most a factor of two**. It is enormously tempting to point that at octree LOD
and declare a cap. **Do not.** The theorem's hypotheses are: the class is **convex**, the class is
**symmetric**, the error is **worst-case**, and the solution operator is a **continuous linear
mapping**. Recovering a level set is not a continuous linear operator on a convex symmetric class, and
the same literature records that for nonlinear problems with restricted information adaption can gain
**up to order `n`**. So the factor-2 result is the wrong tool, correctly identified — and the right
move is to *measure* the adaptive-vs-uniform gap rather than to cite a bound at it.

---

# Part 2 — vocabulary upgrades on three v1 axes

---

## Axis 3 (ambiguity), upgraded — the cell is a tensor, and its discriminant has a name

v1 said: *"Morse theory | the parent theory. Ambiguity is really 'a critical point falls inside a
cell'."* That was right, and the corpus sweep of 2026-08-23 acted on it. This is the **other** parent
theory, and it was sitting in the crate's own source the whole time.

| Word | What it means | Why you want it |
|---|---|---|
| **Multi-affine** | Degree 1 in each variable separately | The trilinear is not a generic cubic. Saying *multi-affine* opens a literature that saying *cubic* closes |
| **`2×2×2` tensor** | The eight corner values, indexed `f[u + 2v + 4w]` | **The eight corners are a tensor and the trilinear is its multilinear form.** Everything below follows from taking that seriously |
| **Pencil** | The one-parameter family `A₀ + λA₁` of the two opposite-face `2×2` matrices | **[proved here]** The body-saddle quadratic **is** `det(A₀ + λA₁)`, and gives the *same* polynomial for all three axis pairings |
| **Cayley's hyperdeterminant** `Δ` | The degree-4 discriminant of a `2×2×2` tensor | **[proved here]** `b² − 4ac` in `marching_cubes/trilinear.rs:246` is **identically** `Δ` — same 12-term degree-4 polynomial, difference symbolically zero |
| **Relative invariant** | Transforms by a character: `Δ(g·A) = (det g₁ det g₂ det g₃)²·Δ(A)` | **[proved here]** The weight is a **square**, so **`sign(Δ)` is an absolute `GL(2)³` invariant** — the body-saddle count cannot depend on cell aspect ratio or on any per-axis affine reparametrisation |
| **Tensor rank** `rank_⊗` | Fewest rank-1 terms summing to the tensor | **[verified]** de Silva & Lim: `Δ > 0 ⟹ rank ≤ 2`; `rank ≤ 2 ⟹ Δ ≥ 0` |
| **Multilinear rank** `rank_⊞` | The triple of slice-space dimensions | **[verified]** `sign(Δ)` **alone does not** classify the orbit — you need the pair `(sign Δ, rank_⊞)` for the eight real orbits |
| **Border rank** | Rank of a limit of rank-`r` tensors | Why `Δ = 0` is *not* the degenerate stratum: **[verified]** tensors on `Δ = 0` are **generically rank 3** |
| **Newton polytope** | Convex hull of a polynomial's exponent vectors | For `∂f/∂x` of a trilinear it is a unit square in the *other two* coordinates |
| **Mixed volume / BKK bound** | Bernstein's count of solutions in `(ℂ*)ⁿ` | **[verified]** `MV = 8 − 6 + 0 = 2`: a trilinear has at most **2** critical points in the complex torus |
| **Patchworking** | Viro's construction of real hypersurfaces from signs at Newton-polytope lattice points | Signs at the eight corners of a cube, glued over a triangulation. **That is marching tetrahedra, described by a real algebraic geometer in 1980** |
| **Regular / coherent triangulation** | One induced by a convex lifting function | Viro's theorem needs it. Whether the Kuhn triangulation qualifies is the hinge of the whole transfer |
| **Cylindrical algebraic decomposition (CAD)** | Collins' algorithm for real quantifier elimination | Doubly exponential and hopeless **at runtime** — and irrelevant to that, because deriving a case table is an **offline, once-ever** computation |

**Field that owns it:** real algebraic geometry; tensor geometry; symbolic computation.
**Search terms:** `hyperdeterminant 2x2x2 real orbit classification`, `multi-affine hypersurface
topology`, `Bernstein mixed volume Newton polytope`, `combinatorial patchworking real algebraic
hypersurface`, `cylindrical algebraic decomposition quantifier elimination`.

**Why this is worth more than a renaming.** The A-002 series' whole cost has been transcription —
`M-207`'s lost root, `M-219`'s typo, `M-221`'s `0 × NaN`, `✗22`, `M-228`, `M-231`. Naming the quantity
converts assertions into theorems: the count is **octahedrally invariant** and **negation-invariant**
`[proved here]`, which `✗39`/`✗49` currently establish by a 48-element empirical sweep; and it is
**aspect-ratio independent**, which nothing currently establishes at all. `M-206`'s "two independently
derived constructions locate the same body saddles to `1.1e-12`" stops being a happy coincidence and
becomes the three axis pairings of one pencil.

---

## Axis 4 (vertex placement), upgraded — is the input even consistent?

v1's Axis 4 vocabulary is about *solving* the QEF. It has no word for the prior question: **can these
normals have come from one smooth surface at all?**

| Word | What it means | Why you want it |
|---|---|---|
| **Integrability** | A vector field is a gradient **iff** its curl vanishes | The `≤12` edge normals of a `HermiteCell` are samples of `∇f`. If their curl residual is large, no single sheet explains them |
| **Helmholtz–Hodge decomposition** | Split a field into curl-free + divergence-free + harmonic | The measurement instrument. The curl part **is** the inconsistency, as a number, per cell |
| **Discrete exterior calculus (DEC)** | Differential forms on a cell complex; `d`, `∧`, `⋆` as matrices | The computable form on a cube. A few dozen flops on data the crate already builds |
| **Jones' `β`-number** | `β_∞(Q)` = width of the thinnest slab containing `S ∩ Q`, over `diam Q` | **The formal name for "how planar is this cell"** — exactly what the QEF silently assumes and never checks |
| **Analyst's Traveling Salesman Theorem** | `∑ β(Q)²·diam(Q)^d < ∞` characterises rectifiability | Turns multiscale flatness into a *convergent sum* — i.e. an a-priori budget, not a per-cell heuristic |
| **Varifold** | A measure on (point, tangent-plane) pairs | Curvature from discrete data **with convergence bounds**. The rival to the normal cycles already benched |
| **Second fundamental form** | The curvature tensor | v1 already named it. Still unexplored, and `β` is how you'd know where it matters |

**Field that owns it:** geometric measure theory; quantitative rectifiability; discrete exterior calculus.
**Search terms:** `Helmholtz-Hodge decomposition discrete vector field`, `Jones beta numbers
quantitative rectifiability`, `analyst traveling salesman theorem higher dimension`, `varifold
approximation mean curvature convergence`.

**The hook.** `λ = 0.01` in `dual_contouring/solve.rs` is a **global guess about how inconsistent the
QEF's input might be.** Hodge measures that inconsistency **directly, per cell**. And `M-60` says only
two of seven fields ever need a second vertex in a cell — with no predictive rule behind the split.
A curl residual is a candidate rule, and it has a reason to fire exactly where the project already
hurts: CSG `min`/`max` makes the gradient discontinuous.

---

## Axis 6 (guarantees), upgraded — say which theorem, and check its hypothesis

| Word | What it means | Why you want it |
|---|---|---|
| **Regular value / transversality** | `∇f ≠ 0` on the level set; the surface meets each simplex cleanly | **The hypothesis every PL-correctness theorem needs.** Nobody in this project has measured how often it fails |
| **PL continuation / simplicial pivoting** | Following an implicitly defined manifold through a simplicial subdivision | **[verified]** Allgower & Schmidt, `10.1137/0722020`, SIAM J. Numer. Anal. 22(2):322–346, **1985** — two years before Lorensen & Cline |
| **Residual bound vs topological guarantee** | `‖H(x)‖_∞ < ε` vs *homeomorphic to the true set* | **[verified]** AS85 proves the **first**, under a full-rank hypothesis. It does **not** prove a homeomorphism — do not overclaim the precedence |
| **Isotopy vs homeomorphism** | Deformable to, within the ambient space, vs merely abstractly the same | The strong statement. v1 named it; the modern isotopy theorem for isomanifolds is Boissonnat & Wintraecken, SoCG 2020 |

**The precedence claim, stated exactly.** PL approximation of an implicitly defined manifold by
simplicial subdivision **predates Marching Cubes by two years and this project has never cited it**.
What it delivers is an `ε`-residual bound plus a regularity hypothesis — genuinely less than the
topological guarantee it is tempting to attribute to it. Both halves belong in the ledger: the
priority *and* the limit.

---

# Part 3 — how to use this

**The move that works, restated.** Take the problem, find its word, then search *the word plus its
source field*. v1's example still holds: `crack-free LOD voxel` returns blog posts; `conforming
adaptive refinement hanging node` returns forty years of FEM with proofs. Three new instances:

- "my triangles are the wrong shape" → **not** `mesh quality voxel` → `metric-based anisotropic mesh
  adaptation continuous mesh`
- "is this cell flat enough for a plane fit" → **not** `planarity test` → `Jones beta number
  quantitative rectifiability`
- "how good could this possibly get" → **not** `optimal meshing` → `n-th minimal error information-based
  complexity`

**Two new triage signals, from what this sweep caught.**

- **Does the paper's headline number compare like with like?** Cao's Table 2 shows a real 5.1× spread
  in `L¹` error at matched element count — and the two metrics being compared are optimised for
  **different norms**, and the `L¹` column is the one place his own theory does not hold. The numbers
  are right; the story is not.
- **Does an optimality claim name its class?** "BCC is the optimal 3D quantizer" is true among
  *lattices* and unproven in general. "Adaption buys at most 2×" is true for *linear operators on
  convex symmetric classes* and false in general. **An optimality claim without its class is folklore
  wearing a theorem's clothes.**

**Words that now signal a paper is worth your time here**, extending v1's list: *approximation order*,
*asymptotic constant*, *complexity* (in the `∫√det M` sense), *unit mesh*, *relative invariant*,
*mixed volume*, *regular triangulation*, *`n`-th minimal error*, *approximation class*, *curl-free*,
*regular value*.

**Words that still signal folklore**, extending v1's: *optimal* (no class), *provably better* (no
denominator), *quality* (unseparated from shape, size and alignment), *fractal* (no measured exponent),
*information-theoretic* (used as an intensifier).

---

# Part 4 — the six transfers I'd rank first

All are registered as `P-103`…`P-152` in
`docs/research/2026-08-29-fifty-experiments-from-unmined-mathematics.md`. Ranked by *what they change
if they land*, not by how likely they are to land.

| # | Transfer | From | The move | Falsified by |
|---|---|---|---|---|
| 1 | **The cell's discriminant is Cayley's hyperdeterminant** `[proved here]` | Tensor geometry / real algebraic geometry | The body-saddle test becomes a `GL(2)³` relative invariant with a published real-orbit classification, instead of a transcribed quadratic | Nothing — the identity is symbolic. What can fail is the *consequences*: that orbit class predicts defects, and that exact-sign evaluation changes any mesh |
| 2 | **The gyroid has an exact Euler characteristic after all** `[verified]` | Minimal surface theory | `χ = −8N³` per `N³` conventional cubic cells under periodic wrap, closing a hole `CLAUDE.md` records as unclosable | The extraction not being periodic-conforming, which would mean the oracle exists and this crate cannot reach it |
| 3 | **Anisotropy moves the constant M-12 fits** `[verified]` | Continuous-mesh theory | Drive element shape from `det\|H\|^(−1/(2p+d))·\|H\|` instead of leaving it isotropic by default | The AM–GM gap being small on our fields — i.e. no feature having a flat direction, which would mean there is no constant to recover |
| 4 | **Curl residual replaces a global `λ`** `[F]` | Discrete exterior calculus | Measure per-cell Hermite inconsistency directly instead of guessing it once, globally, at `0.01` | The residual failing to separate `M-60`'s second-vertex cells from the rest |
| 5 | **A denominator for the whole crate** `[F]` | Information-based complexity | State the `n`-th minimal error for the field class, and report the extractor's distance from it | The minimal-error rate already being `Θ(h²)`, which would say Marching Cubes is order-optimal and only the constant is in play — **which is itself the most valuable negative available** |
| 6 | **Marching tetrahedra is combinatorial patchworking** `[F]` | Toric / real algebraic geometry | A 1980s theorem about when sign-interpolation is isotopic to the true hypersurface | The Kuhn triangulation not being regular, which kills the transfer outright in one afternoon |

**Why #1 first.** It is the only entry on this list that is already *true* rather than plausible, it
was found inside this crate's own source rather than in a paper, and it converts three empirical facts
(`M-206`'s coincidence, `✗49`'s 48-element sweep, the untested aspect-ratio question) into corollaries
of one invariant. v1 ranked self-adjusting computation first on the argument that *the axis where
everyone is identical is where the empty space is*. That argument still holds — and Axis 3 turned out
to have a second empty region nobody had looked at, because the field that owns it does not read
graphics papers and the graphics papers do not use its words.
