# SDF corpus build-out — and five findings that change the project

**Date:** 2026-08-15
**Method:** five parallel agents over home-still. Local corpus first (`catalog_read` for presence,
never search), gaps filled by `paper_search`, acquired via `paper_download` → `scribe_convert` →
`distill_index`.
**Result:** ~60 papers acquired. The corpus was **barren** on level sets, eikonal solvers,
redistancing and distance transforms — no Osher–Sethian, no fast marching, no EDT literature at all.
That is now covered.

---

# Part 1 — The five findings

## ✗ 1. The `is_exact_distance()` bug has a name, a citation, and a sign

Marschner, Sellán, Liu & Jacobson 2023 (`10.1145/3610548.3618170`) name it exactly: min/max CSG
produces a **Pseudo-SDF** — *eikonal almost everywhere yet not a distance function*, with error
**concentrated at seams**, specifically at the union's medial axis.

That is precisely where `csg_difference`'s defects concentrate (3 non-manifold edges, 6 non-manifold
vertices, 6 inconsistently oriented). Last turn this was a hypothesis. It now has a paper.

**And the error is one-signed, which decides which operation is dangerous:**

| Operation | Effect on distance | Safe? |
|---|---|---|
| `min(f,g)` — union | **Never overestimates.** Conservative lower bound, stays 1-Lipschitz | Safe |
| `max(f,g)` / `max(f,−g)` — intersection, subtraction | **Overestimates near concave seams** | **Unsafe** |

Overestimation is the direction that lets a tracer step through a surface and mis-places an
interpolated vertex. `csg_difference` is `max(box, −sphere)` — the unsafe operation — declaring
`is_exact_distance() -> true` with the comment `// away from the seam`.

**It is wrong in the dangerous direction.** The literature says that caveat is exactly the one that
cannot be dropped.

## 2. The flag should have been a ratio, and there is a published algebra for it

Bálint, Valasek & Gergó (`10.14232/actacyb.24.1.2019.3`, acquired) prove every SDF is Lipschitz with
**smallest constant exactly 1**. Their follow-up gives a `q ∈ (0,1)` **underestimate-ratio**
formalism that composes through a CSG tree.

That is the type `is_exact_distance() -> bool` should have been. `// away from the seam` is a `q < 1`
field declared `q = 1`.

## ✅ 3. The 576-evaluations-per-cell optimisation is safe — and this is the big one

**The Lipschitz bound survives arbitrary CSG; exactness does not.**

`max`/`min` of two 1-Lipschitz functions is 1-Lipschitz. So a CSG combination of exact SDFs remains a
valid conservative bound field, forever, no matter how many brush strokes. Hart's sphere tracing
(`10.1007/s003710050084`, acquired) is built on exactly this: with `L = 1`, a **single** evaluation at
a cell centre with `|f| > half-diagonal` proves the entire cell empty.

**M-98's 70× subgrid cost is therefore attackable with a provably-correct optimisation that survives
player editing** — even though `is_exact_distance()` should become false after the first brush stroke.
Kalra & Barr give the gradient-bounded version for when the field is only `q`-bounded.

## ✗ 4. LFS-driven LOD is dead as I proposed it — and that's a theorem, not an engineering problem

Transfer #3 from the axes doc was "drive refinement by local feature size, the theoretically correct
criterion." **It cannot be computed.**

LFS(x) = dist(x, medial axis), and the medial axis is the least stable object in the theory: an
arbitrarily small Hausdorff perturbation of the boundary creates arbitrarily long spurious branches
(Attali–Boissonnat–Edelsbrunner, acquired). LFS from a sampled grid **is not even continuous in the
input**. Aamari et al. confirm it quantitatively — estimating merely the global infimum (the reach)
converges at slow minimax rates and needs C³/C⁴ regularity. Statistics-grade offline, not per-chunk.

**Three affordable replacements, in descending directness:**

- **Curvature straight from the field.** On a true SDF, `∇d` is unit and the Hessian's nonzero
  eigenvalues at a surface point are `−κᵢ/(1−κᵢd)` — principal curvatures fall out of samples you
  already take, no medial axis involved. **This is the cheap LOD driver.** Aamari & Levrard bound how
  noisy the estimate is; Cuel–Lachaud–Thibert give the grid-specific stable version.
- **λ-medial axis / scale-axis transform** — provably stable filtrations, interactive rates on grids.
- **μ-reach / weak feature size** — the theoretically correct stand-in whenever input is sampled.

## ✗ 5. Hausdorff error does not certify topological correctness — provably

Two surfaces can be arbitrarily Hausdorff-close and not homeomorphic (Schwarz lantern). Every real
theorem adds a **second hypothesis**, always either a feature-size bound or a normal-variation
condition.

The isosurface-specific results — **Plantinga & Vegter**, **Boissonnat–Cohen-Steiner–Vegter** (both
already in corpus) — certify **isotopy from a per-cell normal-variation condition**, not from a global
Hausdorff bound.

**That is local, cheap, checkable during extraction, and a natural fit for a marching pipeline.** It
is the thing to build. The alternative — Hausdorff < c·(μ-)reach — requires estimating reach and is
expensive.

Also note Morvan & Thibert: normal and area convergence require an **angle condition** on the
elements. Small `h` alone is not enough. `M-12`'s `h²` result is about position, and does not
transfer to normals for free.

---

# Part 2 — Two original measurements now available

Both fall out of the reading, both fit the existing harness, and neither exists in the literature.

**A. The degradation rate of the distance property under repeated CSG.** No paper measures error
growth as a function of operation count. The experiment: apply N random sphere subtractions to an
analytic box, sample `‖∇f‖` on a grid, plot the distribution against N. Sits alongside the manifold
and self-intersection metrics. Genuine contribution rather than a re-derivation.

**B. A head-to-head of downsampling operators for SDFs.** Mean vs min vs re-evaluate vs wavelet. **That
comparison does not exist.** And the literature predicts your result: the correct answer is *you don't
downsample, you re-sample* — every level built by evaluating the underlying distance function at that
spacing, gated by a per-cell error predicate (Frisken's ADF, Koschier's hp-adaptive). Under
re-sampling a plate thinner than a coarse cell yields all-positive corners and correctly disappears.
Under box-filter averaging the straddling ± set survives and MC keeps emitting triangles — **which is
exactly M-72's measured 4,088 → 1,016 → 248 → 56.** Your aliasing bug is the predicted failure of an
operator the literature already rejects.

---

# Part 3 — The narrower fix nobody would have guessed

Worth weighing before building redistancing at all.

**Extraction reads the sign and the zero crossing. min/max preserve the sign *exactly*** —
`{min(f,g) ≤ 0}` **is** the union. What min/max break is (a) the **interpolated crossing position**,
because the field is kinked and linear edge interpolation is wrong near a seam, and (b) gradient and
normal estimation.

That is far narrower than redistancing the whole field. The cheap fix is **stop interpolating linearly
across a known kink** — and you already know where the kinks are, because you built the CSG tree.

*(Flagged by the agent as its own reading rather than a cited claim, though consistent with Pujol &
Chica, who treat exactly this kink problem and are already in the corpus.)*

If you do want redistancing: **narrow-band PDE reinitialization** (Peng et al., acquired) scales with
edited *surface area* not chunk volume — the best structural match to a brush stroke. **Fast sweeping**
(Zhao, acquired) is O(N), a few Gauss–Seidel passes, no heap. **Jump flooding** (acquired) is the GPU
option. Watch Sussman & Fatemi's point: naive reinitialisation *moves the zero set*, which in a
destructible game means geometry creeping after every edit.

---

# Part 4 — One more correction: the standard GWN citation is obsolete

Barill et al. 2018 (`10.1145/3197517.3201337`) has been the default citation for a decade. The 2026
Antipodal Method paper's accuracy comparison states its order-0 and order-1 expansions are *"very
imprecise… can be considered not useful for applications."*

Two 2026 papers replace it, **exact and faster**: Antipodal (`10.1145/3811323`) and Xie, Hafner &
Wojtan (`10.1145/3811339`), both acquired. Shared insight: the winding number reduces to one
ray-surface intersection plus a sum over **boundary** edges — **cost scales with the number of holes,
not triangles.** A nearly-closed mesh is nearly free.

Two design notes: **two-tier is right** — pseudonormal (Bærentzen & Aanæs, acquired, and it's a proof
not a heuristic) for geometry isomesh produced itself and already guards with `V−E+F == 2`; exact GWN
for imported or carved input. And **use GWN to classify points, not to repair meshes** — Takayama et
al. 2014 (acquired) is the GWN authors' own four-page publication that the tempting orientation-repair
application is *"fundamentally flawed."*

---

# Part 5 — Corpus hygiene: a new failure mode

**AMS DOIs are a trap.** `paper_download` reports success but retrieves the *journal landing page*.
Five junk records are now catalogued and indexed with garbage text — **purge these**:

```
10.1090_s0002-9947-1983-0690039-8      (Crandall & Lions — viscosity solutions)
10.1090_s0002-9947-1984-0732102-x      (10 chunks)
10.1090_s0025-5718-07-01981-3          (7 chunks)
10.1090_s0025-5718-04-01678-3          (Zhao fast sweeping)
10.1137_060670298                       (also an UNVERIFIED DOI — guessed from memory, 1 KB stub)
10.14733_cadconfp.2022.329-333          (CAD'22 table of contents, 6 chunks)
```

All the AMS ones are free in the AMS Digital Archive (>5 years) by direct URL. This is the same
landing-page signature as the five found on 2026-08-13 — **`pdf_path` ending in `.html` with a low
chunk count** remains the detector, and it now has a second producer.

**Two tool behaviours worth writing into `CLAUDE.md`:**

- **`catalog_read` stems are case-sensitive.** `10.1109_TVCG.2006.56` misses; `10.1109_tvcg.2006.56`
  hits. An agent nearly reported two present papers as absent.
- **`paper_search` with `provider: arxiv` behaves as exact-phrase match.** Multi-word conceptual
  queries silently return `[]` while short real title fragments work. Combined with the broken CORE
  ranking already recorded, `provider` selection is now load-bearing for every search.

---

# Part 6 — Hand-acquisition list (consolidated, DOIs verified)

**Highest value first.**

| Paper | DOI | Why |
|---|---|---|
| Bálint et al., *Operations on SDF **Estimates*** | `10.14733/cadaps.2023.1154-1174` | **The closest published error calculus for composed bound fields.** Single highest-value item |
| Museth, **VDB** | `10.1145/2487228.2487235` | The storage paper everything cites |
| Balsa Rodríguez et al., *Compressed GPU-Based Direct Volume Rendering* STAR | `10.1111/cgf.12280` | The systematic rate–distortion treatment; biggest hole for bits-per-voxel |
| Federer, *Curvature measures* | `10.1090/S0002-9947-1959-0110078-1` | Origin of reach. AMS archive, free |
| Barill et al., *Fast winding numbers* | `10.1145/3197517.3201337` | Now superseded but universally cited — needed to read the 2026 critiques |
| Sethian, *Fast marching* | `10.1073/pnas.93.4.1591` | **Free on pnas.org**, Unpaywall just lacks it |
| Sethian & Vladimirsky | `10.1073/pnas.090060097` | Same |
| Kalra & Barr, *Guaranteed ray intersections* | `10.1145/74334.74364` | The `q`-bounded sphere-trace |
| Sharp & Jacobson, *Spelunking the Deep* | `10.1145/3528223.3530155` | Range analysis for cell emptiness |
| Enright et al., *Hybrid Particle Level Set* | `10.1006/jcph.2002.7166` | |
| Sussman & Fatemi, *Interface-Preserving Redistancing* | HAL: `hal.science/hal-01694576` → `redistance.pdf` | Free PDF located |
| Mullen et al., *Signing the Unsigned* | `hal.inria.fr/inria-00502473/file/signing.pdf` | Free PDF located |
| Segment Tracing | `10.1111/cgf.13951` | HAL PDF already in the catalog entry |
| Lekien & Marsden, *Tricubic interpolation* | `10.1002/nme.1296` | |
| Koschier et al., *hp-adaptive SDF generation* | `10.1109/tvcg.2017.2730202` | |
| Ricci 1973 | `10.1093/comjnl/16.2.157` | **Origin of min/max CSG** |
| Madoš et al., *CSVO* | `10.3390/sym14102114` | MDPI is gold OA — a resolver miss, retry manually |

Weakest area overall is **R-functions** — Rvachev, Shapiro, Pasko all paywalled. Only Fryazinov 2010
and Hybrid F-rep cover F-rep at all.

**Books this pipeline cannot fetch:** Cannarsa & Sinestrari *Semiconcave Functions, HJ Equations and
Optimal Control* (the canonical semiconcavity reference), Bardi & Capuzzo-Dolcetta, Lions *Generalized
Solutions of HJ Equations*, Federer *Geometric Measure Theory*, Dey *Curve and Surface Reconstruction*,
Delfour & Zolésio *Shapes and Geometries*.

---

# What I'd do

| | Action | Why |
|---|---|---|
| 1 | **Change `is_exact_distance() -> bool` to a `q ∈ (0,1)` bound**, and set `csg_difference` honestly | It is wrong in the dangerous direction *today*, and there is a published algebra for the right type |
| 2 | **Empty-cell rejection via the Lipschitz bound** | Attacks M-98's 70× with an optimisation that is provably correct under arbitrary player editing |
| 3 | **Per-cell normal-variation isotopy check** | Turns "we report Hausdorff" into "we certify topology," locally and cheaply. Plantinga–Vegter, already in corpus |
| 4 | **Measure the CSG degradation rate** | Original contribution, fits the existing harness, ~a day |
| 5 | **Purge the six polluted records** | Two producers of the same signature now; it will recur |
| 6 | Curvature-from-Hessian as the LOD driver | The surviving half of a transfer whose original form is provably impossible |
