# Meshing decomposed into axes, v3 — four more axes, and the words for grading a mesh without a mesh

**Date:** 2026-09-04
**Supersedes nothing in** `2026-08-12-axes-and-vocabulary.md` (ten mechanism axes) or
`2026-08-29-phase-27-axes-and-vocabulary-v2.md` (four "how well could it be done" axes). Both are still the
right first reads. This document names **four axes that v1 and v2 both missed**, upgrades the vocabulary on
**five axes where Phase 27's results outgrew v2's words**, and records the reading rule that produced them.

**How to read the tiers.** Vocabulary is stable and safe. **[verified]** means the settling sentence was read
in the primary source during this session and is quoted in
`docs/research/2026-09-04-phase-29-fifty-experiments-from-twenty-six-fields.md`. **[abstract]** means the
claim was checked against the paper's abstract only — the PDF is in home-still and still converting, and the
registration that leans on it says so. **[F]** is a hypothesis in `FINDINGS.md` terms.

---

## 0. What v2 got right, and the shape of what it missed

v2's move was to add a *denominator* — Axes 11–14 asked how well the thing could possibly be done. Phase 27
then cashed that in, and the four most consequential results were all denominators:

- **`M-472`**: against `W_inf^2`, Marching Cubes is order-optimal; only the constant is in play, and the
  constant gap is bounded at **6.3–7.3×**.
- **`✗116`**: a fourth-order filter reconstructs at order 4 and the mesh still converges at order 2 — *"the
  second order belongs to the piecewise-linear mesh and no filter can reach it."*
- **`✗107`/`✗108`**: anisotropy failed, and the finding was about the *roster*: five of eight reference fields
  make the two arms one grid, and the one field where the mechanism works is `Unbounded` and therefore
  ungradeable.
- **`M-474`**: the approximation class `s` is field-dependent — LOD helps on three fields and is measurably
  *worse* than uniform on `fbm_terrain`.

Read together, those four say something v2 could not phrase: **every grade this project can give a mesh is
computed from another mesh or from a bound the field must volunteer.** `validate::accuracy` needs
`field.bound() == Exact`; the χ oracle of `P-145` is the first instrument that reads the field's own signs
and it arrived in Phase 27. There is a mathematics whose whole business is measuring area, curvature and
topology of a set from a lattice of samples *without constructing a surface* — stereology and integral
geometry — and this project has never used a word of it. That is Axis 15. The other three missing axes follow
from the same reading: the field as a *random object* with computable expectations (Axis 16), the output's
*representation order* as a design choice distinct from the filter's (Axis 17 — the lever `✗116` named and
nobody has pulled), and *time* (Axis 18 — v1's Axis 8 asked what happens after an edit, never what happens
during one).

---

# Part 1 — the four missing axes

---

## Axis 15 — Mesh-free measurement (integral geometry, stereology, digital geometry)

**What it is:** the instruments that grade an extraction from the sign grid alone, so that a field with no
`Exact` bound can still be graded, and so that the mesh has a reference that is not another mesh.

| Word | What it means | Why you want it |
|---|---|---|
| **Intrinsic volumes** `V₀…V₃` / **Minkowski functionals** / **Lipschitz–Killing curvatures** | The four additive, motion-invariant, continuous functionals of a 3-D body: volume, surface area, integrated mean curvature, Euler characteristic — three names for one family, up to normalisation | **[verified]** Taylor 2006 states the normalisation: `M_{k−j}/(k−j)! = L_j · ω_{k−j}`. The whole of what a mesh is *for*, as four numbers a lattice can estimate directly |
| **Hadwiger's theorem** | Every continuous, additive, motion-invariant functional on convex bodies is a linear combination of the intrinsic volumes | The reason there are exactly four things to measure and not five. A "quality" number that is not one of these is not a geometric invariant |
| **Steiner formula / Weyl's tube formula** | `Vol({x : d(x,M) ≤ r})` is a **polynomial in `r`** whose coefficients are the intrinsic volumes, valid **for `r` below the reach** | **[verified]** quoted from Taylor 2006 §1. For an SDF `f`, the offset body is `{f ≤ r}` for free — so the polynomial's breakdown radius *is* a reach measurement |
| **Crofton formula** | Area = (const) × number of intersections with a family of lines; on a lattice, count sign changes along rows in the 13 lattice directions | The mesh-free area and mean-curvature estimator. **[verified]** Legland, Kiêu & Devaux 2007: `S = 4∫χ(L) dL`, `b̄ = ∫χ(P) dP`, discretised over 3 or 13 lattice directions with Ohser–Mücklich spherical-Voronoi weights, everything reduced to `2×2×2` configuration counts |
| **Local configuration weights / `h`-function** | Assign a weight to each of the 256 `2×2×2` sign configurations; area = Σ weights. Ziegel & Kiderlen prove the configuration counts are asymptotically ∫ of an `h`-function against the surface area measure | **[verified]** *"the relative worst case error is asymptotically less than 4%"*, best possible among local weightings; only five configuration types carry weight, `λ̄ = (sr, 1.358, 2.352 − sr, 0.960, s(1 − r))`; the `h`-functions span a **50**-dimensional space, so 50 is the number of independent things a `2×2×2` census can know about the surface-area *measure*. **The same 256 configurations the case table indexes carry a provably ≤4% area estimator, and the crate has never summed it** |
| **Adjacency pair** `(F, B)` | Foreground and background adjacencies chosen *together* — `(6,26)`, `(26,6)`, `(18,6)`, `(6,18)` — so that Euler characteristic is consistent by a theorem rather than a choice | The name for what `P-145`'s `connectivity_split_blocks` (**487** on `noise_cavity`) is counting. Ohser–Nagel–Schladitz call the consistent pairs *complementary* |
| **Multigrid convergence** | An estimator on a digitised shape converges to the true value as the grid step `h → 0`, with a stated order | **The digital-geometry word for approximation order.** Cuel, Lachaud & Thibert 2014 (in the corpus) prove it for a Voronoi-covariance normal estimator with an explicit `h^{1/2}` term **[verified, Theorem 1]**; Coeurjolly, Lachaud & Levallois 2014 do the same for principal curvatures (acquired, HTML landing only — rate unread). The term retrieves forty years of results on exactly this project's data structure |
| **Directional discretisation error** | The part of a Crofton estimate that does **not** vanish with resolution, because only 3 or 13 lattice directions are sampled | **[verified]** Legland et al.: bounded by `−0.3333·S ≤ ε ≤ +0.1547·S` for 3 directions and `−0.0733·S ≤ ε ≤ +0.0227·S` for 13; an axis-aligned cube reads **1600** against a true **2400** with 3 directions and **2160.4** with 13. The *error Wulff plot* of Group A is this bound, resolved over the sphere of normals |
| **Digital plane / arithmetic plane** | `{(x,y,z) ∈ ℤ³ : μ ≤ ax+by+cz < μ+ω}`; *naive* when `ω = max(\|a\|,\|b\|,\|c\|)` | **[landing page only]** Brimkov, Coeurjolly & Klette 2007, *Digital planarity — a review*; the HTML the resolver returned carries the abstract, not the body. The staircase a flat wall leaves in the sign grid is one of these, and *recognising* one is a solved problem — which makes greedy meshing of non-axis-aligned planes a search term, not a hope |
| **Functional vs measure convergence** | The total (χ, area) converges with resolution; the *local* Euler–Poincaré and mean-breadth measures do not | **[verified]** Legland et al.'s conclusion — a lattice reconstruction has finitely many normal directions, so per-voxel curvature measures cannot converge. Any per-cell curvature instrument built on the sign grid inherits this; the per-cell *sum* does not |
| **Coarea formula** | `∫ \|∇f\| dx = ∫ Area({f = t}) dt` | Area of the zero level set from a grid integral of `\|∇f\|·δ_ε(f)` — a third mesh-free area, from the *values* rather than the signs |
| **Reach** `τ` (Federer) | The largest `r` such that every point within `r` of `M` has a unique nearest point on `M` | The regularity number that both the tube formula and the sampling theorems (Axis 6, upgraded) are stated in. For a sphere `τ = R`; for a torus `τ = r_minor`; for a polyhedron `τ = 0` |

**Field that owns it:** integral geometry (Hadwiger, Federer); stereology and image analysis (Ohser, Schladitz,
Kiderlen, Legland); digital geometry (Coeurjolly, Lachaud, Klette).
**Search terms:** `intrinsic volumes binary image lattice`, `Crofton formula discretization 3D`, `local
configuration weights surface area estimation`, `multigrid convergence digital estimator`, `digital plane
recognition`, `complementary adjacency Euler number`.

**What this buys immediately.** `✗107` could not grade `fbm_terrain` because `validate::accuracy` is
meaningless where `field.bound()` is not `Exact`. Axis 15 says area, mean curvature and χ can be estimated
from the sign grid with *proven* worst-case error, for any field whatsoever, and the mesh's own area, curvature
and χ can be graded against them. It also says the four Phase 27 fixtures that read "unmeasurable" were
unmeasurable by *one* instrument.

---

## Axis 16 — The field as a random object (random-field geometry)

**What it is:** two of the eight reference fields are noise. Everything the ledger says about them is a
per-seed reading. There is a theory that gives the *expectation* of exactly the quantities this project
measures — number of critical points, area of the level set, Euler characteristic — in closed form from the
covariance, and it has been used in brain imaging and cosmology for thirty years.

| Word | What it means | Why you want it |
|---|---|---|
| **Gaussian random field (GRF)**, **stationary**, **isotropic** | Every finite set of samples is jointly Gaussian; statistics translation-invariant; rotation-invariant | **The hypothesis every closed form below needs.** `P-176` already established `fbm_terrain` is *not* one (hash-based lattice noise). A GRF fixture is therefore a prerequisite, not an option |
| **Covariance function** `ρ(‖t−s‖²)`, **spectral moments** `ρ'(0)`, `ρ''(0)` | The first two derivatives of the covariance at zero | **[verified]** Cheng & Schwartzman: *"the obtained formulae depend on the covariance function only through its first and second derivatives at zero."* Two numbers predict the whole census |
| **Kac–Rice formula** | `E[#{t : ∇X(t) = 0, …}] = ∫ E[\|det ∇²X\| · 1_{…} \| ∇X = 0] p_{∇X}(0) dt` | The expected number of critical points — i.e. of body saddles, of ambiguous cells, of components — before any mesh exists |
| **Excursion set** `A_u = {X ≥ u}` / **nodal set** `{X = 0}` | The solid, and the surface, in this vocabulary | The `u = 0` slice of a theory built for every `u` |
| **Gaussian kinematic formula (GKF)** | `E[χ(M ∩ X⁻¹[u,∞))] = Σ_j L_j(M) · ρ_j(u)` — Lipschitz–Killing curvatures of the *domain* times **EC densities** of the *field* | **[verified]** Taylor 2006, Theorem 4.1 as stated in §1. The expected Euler characteristic of the solid, in closed form, on a cube with its boundary accounted for |
| **GOE / GOI matrices** | The Hessian of an isotropic GRF at a critical point is a Gaussian orthogonally-invariant random matrix | **[verified]** the mechanism by which Cheng & Schwartzman evaluate Kac–Rice; the `N=2` example is `E[μ₀] = 1/(√(3π) η²)` per unit area |
| **Level-set percolation**, `ℓ_c` | The threshold above which the excursion set has an infinite component | `M-489` measured `0.029340` on `noise_cavity`. For a GRF this is a *universal* number of the covariance class, and the question of whether the value-noise field sits on it is now askable |
| **Chen–Stein / Poisson approximation** | Rare, weakly dependent events are Poisson-distributed to a bounded total-variation error | Ambiguous cells are rare events at large correlation length. Their count's variance across seeds is then predicted, which is what `P-144`'s per-seed χ variance (**10.24**) lacks |

**Field that owns it:** probability (Adler, Taylor, Worsley, Azaïs–Wschebor); its consumers in neuroimaging
and cosmology (Schmalzing & Buchert's Minkowski functionals of the cosmic web are Axis 15 applied to Axis 16).
**Search terms:** `expected Euler characteristic excursion set Gaussian random field`, `Kac-Rice expected
number of critical points isotropic`, `Gaussian kinematic formula`, `level set percolation Gaussian field`.

**The precision that matters.** Every formula here is for a *Gaussian* field. `fbm_terrain` and
`noise_cavity` are not Gaussian and the theory says nothing about them until a field that *is* Gaussian is
in the roster and the value-noise fields are compared against it. That is the first registration of Group B and
it is a fixture, not a result.

---

## Axis 17 — Output representation order (the lever `✗116` named)

**What it is:** v2's Axis 11 was about the *reconstruction* filter. `✗116` then measured a filter of order 4
producing a mesh of order 2, and named the cause: the mesh is piecewise linear. The order of the *output
representation* is a separate design axis, and every extractor in this crate sits at the same point on it.

| Word | What it means | Why you want it |
|---|---|---|
| **Geometric order** vs **approximation order** | A PL surface through exact points is `O(h²)` from a smooth surface *regardless of how exact the points are* | The sentence `✗116` measured. Vertex accuracy and surface accuracy are different quantities with different orders |
| **PN triangle** (curved point-normal triangle) | A cubic Bézier triangle built from three vertices and three normals, with a quadratic normal field | Vlachos et al. 2001. Dual Contouring already carries the Hermite data it needs. Wang & Olano 2015 applied it to Marching Cubes — *paywalled, cited as prior art, unread* |
| **Bézier / Hermite patch**, **rational-quadratic patch** | Polynomial or rational surface pieces with `C⁰`/`C¹` joins | **[abstract]** Wiley et al. 2003: the isosurface of a *quadratic* element is *"a rational-quadratic patch in parameter space and a rational-quartic patch in physical space"* — an exact per-cell output for the triquadratic `P-138` priced |
| **Superconvergence** | Error at particular points (nodes, Gauss points) decays faster than the global error | If the *vertices* are `O(h⁴)` under the tricubic filter while the surface is `O(h²)`, a curved patch through those vertices can reach order 3–4. Nobody has measured vertex error separately from surface error |
| **Asymptotic error expansion / Richardson extrapolation** | `e(h) = C h² + D h⁴ + …` ⇒ `(4 v_{h/2} − v_h)/3` cancels the `h²` term | Two meshes at `h` and `h/2` share edges; extrapolating vertex positions is free if the expansion is regular. `✗42`'s crossing-position lottery is exactly a failure of regularity, so this can fail — which is why it is a registration |
| **Isogeometric** | Use the CAD representation (splines) *as* the analysis representation | Hughes, Cottrell & Bazilevs 2005. The CAD consumer's vocabulary for "give me a surface, not a triangle soup" |

**Field that owns it:** CAGD; higher-order FEM (superconvergence: Wahlbin, Zienkiewicz–Zhu); numerical
analysis (Richardson, Joyce).
**Search terms:** `curved PN triangles`, `superconvergence nodal error finite element`, `Richardson
extrapolation asymptotic expansion`, `contouring quadratic elements rational patch`, `isogeometric analysis`.

**The honest statement.** A curved output costs the consumer: a Bevy mesh is triangles, and a patch must be
tessellated *somewhere*. The claim is not "ship patches"; it is that the `h² → h⁴` lever `M-472` bounded at
7× lives here and nowhere else, and the ledger should say what it costs.

---

## Axis 18 — Time (kinetic data structures, temporal coherence)

**What it is:** v1's Axis 8 asks what happens after an edit. Nothing asks what happens *during* one — a
brush sliding through a wall changes a continuous field continuously, and the case index of a cell changes
only at discrete moments.

| Word | What it means | Why you want it |
|---|---|---|
| **Kinetic data structure (KDS)** | A structure maintained under continuous motion by a set of **certificates**, each with a **failure time**; an event queue processes certificate failures | Basch, Guibas & Hershberger 1999 (paywalled; the corpus holds Acar et al.'s robust kinetic hulls and kinetic mesh refinement). A cell's certificate is *"sign(f(corner)) is unchanged"*; its failure time is a root of `f(corner, t)` |
| **Responsiveness / locality / compactness / efficiency** | The four KDS quality measures: cost per event, certificates per object, total certificates, events vs. necessary changes | The vocabulary in which "incremental remesh" becomes a bound rather than a hope. `M-322`'s adversarial bisect is an *efficiency* failure in these terms |
| **Time-varying Reeb graph** | The Reeb graph of `f(x,t)` and the discrete events at which its topology changes | Edelsbrunner et al. 2008 (paywalled). Topology-change *times* are computable, so a sealing check need run only at events |
| **Temporal coherence** | Exploiting that frame `t+1` resembles frame `t` | The graphics word; use the KDS words to find the theorems |

**Field that owns it:** computational geometry (Guibas school); visualisation (time-varying isosurfaces:
Sutton & Hansen's T-BON, all paywalled).
**Search terms:** `kinetic data structures certificates events`, `kinetic mesh refinement`, `time-varying Reeb
graph`, `temporal coherence isosurface`.

---

# Part 2 — vocabulary upgrades on five axes

---

## Axis 1 (domain decomposition), upgraded — symmetry reduction and non-lattice sampling

| Word | What it means | Why you want it |
|---|---|---|
| **Asymmetric unit / fundamental domain** | The smallest region from which the group's action tiles the whole | `P-143` named the gyroid's space group `Ia-3d` (order 96 per cubic cell). `✗49` proved the extractor octahedrally equivariant. Together: mesh `1/48` of a symmetric field and *replicate* |
| **Orbifold seam** | A chunk boundary that lies on a mirror plane of the field's symmetry | Vertices on a mirror plane are their own images; a seam there cannot crack. `M-48`'s non-injective seams are the contrast |
| **Jittered / stratified sampling** (Cook 1986) | One sample per cell at a pseudo-random offset | The graphics remedy for aliasing that a regular lattice bakes in. `M-309`'s `0.250` mean ratio on `box_exact` is lattice aliasing with a name |
| **Rational vs irrational ellipsoid** (Landau, Jarník) | Whether the quadratic form's coefficients are commensurable | **[verified]** Ivić et al. §3: the lattice-point discrepancy has different orders for the two. The periodic gyroid sampled at an integer vs non-integer number of cells per period is this dichotomy, and `P-144`'s non-convergent χ is its symptom |

## Axis 3 (ambiguity), upgraded again — the quantum-information dictionary

`P-127` proved the discriminant is Cayley's hyperdeterminant; `P-131` reached for the W-state as a fixture. That
was the tip of a second dictionary, and it has its own continuous measures.

| Word | What it means | Why you want it |
|---|---|---|
| **3-tangle** `τ_ABC` | Coffman, Kundu & Wootters' residual three-way entanglement of three qubits; **equals `4\|Δ\|`**, the hyperdeterminant | **[abstract + folklore]** CKW define *"the three-way tangle … invariant under permutations of the qubits"*; that it is `4\|Det\|` is standard (Miyake 2003) and is what `P-134`'s normalised `\|Δ\|` already was |
| **SLOCC class** | Orbit under `GL(2)³` — *stochastic local operations and classical communication* | **The physicists' name for de Silva & Lim's orbit classification.** Dür, Vidal & Cirac 2000: three qubits have exactly two inequivalent genuinely tripartite classes over ℂ, **GHZ** (`Δ ≠ 0`) and **W** (`Δ = 0`, rank 3, border rank 2) |
| **Entanglement monotone** | A function non-increasing under SLOCC | The W class has monotones that are *not* `\|Δ\|` (which vanishes on it). `✗103` found `\|Δ\|` uninformative on 5 of 8 fields; the W-class monotones are the untested continuous measure on exactly the cells `P-131` counted |

## Axis 6 (guarantees), upgraded — the sampling theorems, and the condition of the problem

| Word | What it means | Why you want it |
|---|---|---|
| **ε-sample**, **reach condition** | Niyogi, Smale & Weinberger 2008: a sample dense enough relative to the reach `τ` recovers the homology with stated probability | **In the corpus at 0.69 and cited nowhere.** `P-140`'s adequacy ladder (`graph_cube_g5` perfect at 7³ and 11³, *nothing* at 9³) is a sampling-density phenomenon with a theorem behind it |
| **Reach estimator** | Aamari et al. 2019 **[abstract]**: *"the first investigation into the problem of how to estimate the reach"*, with minimax rates; Berenfeld et al. 2021 via the convexity defect function | For an SDF the reach is also readable from the tube polynomial (Axis 15). Two independent instruments for one number |
| **Backward error** | The smallest perturbation of the *input* that makes the computed output exact | Per vertex `v`: `\|f(v)\|`; per face: `max \|f\|`. Forward (Hausdorff) error divided by backward error is the **condition number** `≈ 1/\|∇f\|`, which is why `thin_plate` is hard and `sphere` is not |
| **Conformal prediction** | Distribution-free prediction intervals with guaranteed marginal coverage from a calibration set | **[abstract]** Angelopoulos & Bates. The statistical grade for the four fields `validate::accuracy` cannot grade: calibrate on `Exact` fields, predict on `Unbounded` ones, with a coverage number that is a theorem |

## Axis 7/14 (adaptivity and its denominator), upgraded — regularity words that predict `s`

| Word | What it means | Why you want it |
|---|---|---|
| **Besov space** `B^s_{p,q}` | Smoothness `s` measured in `L^p` with fine-scale weight `q`; characterised by wavelet coefficient decay | DeVore 1998 (paywalled): the rate of *nonlinear* approximation is a Besov regularity, the rate of *linear* approximation is a Sobolev one. `M-474`'s `s` is a Besov index measured the hard way |
| **Hurst exponent** `H` | `fbm` with Hurst `H` has Hölder regularity `H` *everywhere* | A field equally rough everywhere has nothing for adaptivity to find — the theorem behind `fbm_terrain`'s **negative** `s_difference` |
| **Box-counting dimension** | `3 − H` for the level set of `fbm` | The exponent the word *fractal* owes and has never paid: triangle count `∝ h^{−(3−H)}` instead of `M-13`'s area law |
| **Tree approximation** | `N`-term approximation constrained to a *tree* of coefficients (an octree is one) | Cohen, Dahmen, Daubechies & DeVore 2001 (paywalled): the tree constraint costs at most a constant. Measurable |
| **Entropy rate** `H(X)` | `lim H(X_n \| X_{n−1}, …)` — bits per symbol of a stationary source | The **denominator for Phase 25**: no bit-packing of the case stream beats the entropy rate. Twenty rows measured numerators |

## Axis 8 (incrementality), upgraded — see Axis 18. The kinetic vocabulary *is* the upgrade.

---

# Part 3 — how to use this

**The move that works, restated for the third time.** Find the word, then search the word *with its source
field*. Four new instances:

- "is the mesh area right" → **not** `mesh area accuracy` → `local configuration weights surface area
  estimation` (stereology)
- "how many ambiguous cells will there be" → **not** `ambiguous cells marching cubes` → `Kac-Rice expected
  number of critical points isotropic Gaussian field`
- "can I beat `h²`" → **not** `higher order marching cubes` → `superconvergence` + `curved PN triangles` +
  `contouring quadratic elements`
- "my wall has a staircase grain" → **not** `marching cubes artefacts flat surface` → `digital plane
  recognition arithmetic plane`

**Three triage signals this sweep added.**

- **Is the claim about a Gaussian field?** Half of Axis 16 is stated for Gaussian fields and none of this
  project's noise fields is one. A GRF fixture comes before any of it.
- **Which order — geometric, approximation, or multigrid?** Three fields use three words for the same
  `O(h^L)`. A paper that says "second order" without saying *of what* is not yet usable.
- **Does the DOI resolve to the paper you think?** Two of this session's twenty-eight downloads resolved to the
  wrong paper — `arXiv:1504.01720` is *Isoperimetric regions in ℝⁿ with density rᵖ*, not Letendre, and
  `10.1016/j.cviu.2014.04.005` is an action-recognition paper, not Coeurjolly. `V-51`'s trap rule is now a
  rate: **~7% of DOIs recalled from memory are wrong.** Both are in the corpus and should be removed.

**Words that now signal a paper is worth your time here:** *intrinsic volumes*, *multigrid convergence*,
*reach*, *Kac–Rice*, *EC density*, *superconvergence*, *certificate* (kinetic sense), *SLOCC*, *Besov*,
*entropy rate*, *rational ellipsoid*.

**Words that still signal folklore:** *fractal* (until Group F measures `H`), *temporal coherence* (no
certificate count), *robust* (no breakdown point), *unbiased* (no proof), *Gaussian noise* (for a hashed
lattice field).

---

# Part 4 — the six transfers I'd rank first

All fifty are registered as `P-177`…`P-226` in
`docs/research/2026-09-04-phase-29-fifty-experiments-from-twenty-six-fields.md`. Ranked by what they change if
they land.

| # | Transfer | From | The move | Falsified by |
|---|---|---|---|---|
| 1 | **A ≤4% area estimator is already in the case table** `[verified]` | Stereology | Sum Ziegel–Kiderlen weights over the 256 configurations the extractor already indexes; grade every field's mesh area against it, `Exact` bound or not | The mesh area and the estimator disagreeing by more than the proven bound on `sphere`, which would mean one of the two is not measuring area |
| 2 | **The mesh order is the output's, and curved patches move it** `[F]` | CAGD / FEM | PN triangles from DC's Hermite data; measure the Hausdorff exponent | The exponent staying at 2, which would locate the bottleneck in vertex accuracy rather than in the PL representation and retire Axis 17 |
| 3 | **A Gaussian field has a closed-form χ, and the noise fields can be compared against one** `[verified]` | Random-field geometry | Add a stationary isotropic GRF fixture with known `ρ'(0)`, `ρ''(0)`; predict `E[χ]`, `E[Area]`, `E[#saddles]` before meshing | The per-seed measurements not bracketing the prediction on the GRF itself, which would mean the fixture is not the field the theorem is about |
| 4 | **The staircase grain is an arithmetic plane and greedy meshing generalises to it** `[verified vocabulary, F transfer]` | Digital geometry | Recognise digital planes in the sign grid and emit one quad per recognised plane on rotated `box_exact` | Triangle count not falling below MC's on a 20°-rotated box, or the recognised planes not being coplanar with the true face to within a cell |
| 5 | **The adequacy ladder is the NSW sampling condition** `[F]` | Learning theory / GMT | Compute each field's reach; predict the first-correct resolution; compare with `P-140`'s ladder | `graph_cube_g5`'s 9³ hole not being explained by any reach-based bound, which would say the phenomenon is aliasing, not density |
| 6 | **Phase 25 has a denominator** `[F]` | Information theory | Entropy rate of the case-index stream vs the best measured bit-packing ratio | The ratio already sitting within 10% of the entropy floor, which would close the bit-packing line for good — the most valuable negative available |

**Why #1 first.** It costs one sum over data the extractor already has, it is provably bounded, it grades the
four fields Phase 27 called ungradeable, and if it disagrees with the mesh on `sphere` — where the mesh area
is trusted — the disagreement is a finding about one of the two instruments either way. v2 ranked the
hyperdeterminant first because it was already *true*. This is ranked first because it is already *computed*.
