# Unobtainable items — consolidated

**Date:** 2026-08-15
**Scope:** everything the acquisition pipeline failed to fetch across this conversation's sweeps —
incrementality (08-13), topology (08-14), SDF build-out (08-15). **Does not include** the earlier
`2026-08-10-meshing-library-target.md` §12 list, which was compiled before this chat and should be
merged separately.

**Organised by recoverability, not by topic** — that's the axis that decides what to do.

> **DOI discipline:** every identifier below was returned by `paper_search` or `catalog_read`. Items
> with no verified identifier are in §5 and are listed by title only. **Do not guess one** — an agent
> on this project once reconstructed an arXiv ID from memory and pulled an unrelated condensed-matter
> physics paper into the corpus under a meshing DOI.

---

## 1. Free — the resolver just failed. Get these first; they cost nothing.

`paper_download` resolves **arXiv + Unpaywall only**. These are all openly available by direct URL.

| Paper | DOI | Where it actually is |
|---|---|---|
| Sethian, **Fast marching level set methods** | `10.1073/pnas.93.4.1591` | Free on `pnas.org` |
| Sethian & Vladimirsky, ordered upwind methods | `10.1073/pnas.090060097` | Free on `pnas.org` |
| Federer, **Curvature measures** (origin of *reach*) | `10.1090/S0002-9947-1959-0110078-1` | AMS Digital Archive (>5 yr, free) |
| Crandall & Lions, **Viscosity solutions of HJ equations** | `10.1090/S0002-9947-1983-0690039-8` | AMS archive |
| Crandall, Evans & Lions | `10.1090/S0002-9947-1984-0732102-X` | AMS archive |
| Zhao, **Fast sweeping method** | `10.1090/S0025-5718-04-01678-3` | AMS archive |
| Cao, interpolation error estimate in ℝ² | `10.1090/S0025-5718-07-01981-3` | AMS archive |
| Sussman & Fatemi, interface-preserving redistancing | — | `hal.science/hal-01694576` → `redistance.pdf` |
| Mullen et al., **Signing the Unsigned** | — | `hal.inria.fr/inria-00502473/file/signing.pdf` |
| **Segment Tracing** | `10.1111/cgf.13951` | HAL PDF already sitting in the catalog entry |
| Madoš et al., **CSVO** | `10.3390/sym14102114` | MDPI, gold OA — pure resolver miss |
| ESA 2020 | `10.4230/lipics.esa.2020.2` | LIPIcs, fully open access |
| JoCG | `10.20382/jocg.v3i1a11` | Journal of Computational Geometry, open access |
| Thäle, *50 years sets with positive reach — a survey* | no resolvable DOI | Surveys in Math. & its Applications 3 (2008) 123–165, open access |

> ⚠️ **The AMS ones are a trap, not merely a miss.** `paper_download` reports **success** on an AMS DOI
> and retrieves the *journal landing page*. Six junk records are catalogued and indexed with garbage
> text and must be purged before re-fetching — see §6.

---

## 2. Paywalled, high value — get by author homepage or institutional access

### 2a. The single biggest hole: Acar's incremental **meshing** line

Output-sensitive incremental meshing, by the person who built the change-propagation theory. **This is
the closest existing prior art to the top-ranked transfer, and not having read it is the main risk to
the "incremental isosurface extraction doesn't exist" claim.** No arXiv version exists. Acar, Blelloch,
Hudson and Türkoğlu all post PDFs on their homepages — usually faster than institutional access.

| Paper | DOI | Venue |
|---|---|---|
| **Dynamic Well-Spaced Point Sets** | `10.1145/1810959.1811011` | SoCG 2010 |
| ↳ journal version, may resolve where ACM won't | `10.1016/j.comgeo.2012.11.007` | Comp. Geom. 2013 |
| **Kinetic Mesh Refinement in 2D** | `10.1145/1998196.1998254` | SoCG 2011 |
| Parallelism in Dynamic Well-Spaced Point Sets | `10.1145/1989493.1989498` | SPAA 2011 |
| **A Cost Semantics for Self-Adjusting Computation** — *the trace-distance paper* | `10.1145/1594834.1480907` | POPL 2009 |

### 2b. SDF / geometry — named directly by Phase 10–13 tickets

| Paper | DOI | Which ticket needs it |
|---|---|---|
| Bálint et al., **Operations on SDF *Estimates*** | `10.14733/cadaps.2023.1154-1174` | **F-001/F-003.** The closest published error calculus for composed bound fields — highest-value single item on this list |
| Ricci 1973 (origin of min/max CSG) | `10.1093/comjnl/16.2.157` | F-003 |
| Kalra & Barr, guaranteed ray intersections | `10.1145/74334.74364` | F-005 — the `q`-bounded sphere trace |
| Sharp & Jacobson, **Spelunking the Deep** | `10.1145/3528223.3530155` | F-005 — range analysis for cell emptiness |
| Koschier et al., hp-adaptive SDF generation | `10.1109/tvcg.2017.2730202` | T-014 |
| Museth, **VDB** | `10.1145/2487228.2487235` | S-004/S-005 storage |
| Enright et al., Hybrid Particle Level Set | `10.1006/jcph.2002.7166` | S-004 |
| Lekien & Marsden, tricubic interpolation | `10.1002/nme.1296` | F-007 |
| Barill et al., Fast winding numbers | `10.1145/3197517.3201337` | S-007 — **superseded, but universally cited**; needed to read the 2026 critiques |
| Balsa Rodríguez et al., Compressed GPU DVR STAR | `10.1111/cgf.12280` | Storage rate–distortion |

### 2c. Distance-function theory

| Paper | DOI |
|---|---|
| Chazal & Lieutier, **The λ-medial axis** | `10.1016/j.gmod.2005.01.002` |
| Lieutier, *Any open bounded subset of ℝⁿ has the same homotopy type as its medial axis* | `10.1016/j.cad.2004.01.011` |
| Cohen-Steiner & Morvan, restricted Delaunay and the normal cycle | `10.1145/777792.777839` |
| Chazal, Cohen-Steiner & Lieutier, normal cone approximation and offset shape isotopy | `10.1016/j.comgeo.2008.12.002` |

### 2d. Contour trees and Reeb graphs

| Paper | DOI |
|---|---|
| Carr, Snoeyink & Axen, **Computing Contour Trees in All Dimensions** — the foundational cost reference | `10.1016/s0925-7721(02)00093-7` |
| Carr et al., **Flexible Isosurfaces** | `10.1016/j.comgeo.2006.05.009` |
| Edelsbrunner et al., Time-Varying Reeb Graphs | `10.1016/j.comgeo.2007.11.001` |
| Distributed Merge Trees | `10.1145/2442516.2442526` |
| Biasotti et al., shape description survey | `10.1145/1391729.1391731` |

---

## 3. Paywalled, lower priority

**Self-adjusting computation foundations** (the theory is well covered by what we do have):

| Paper | DOI |
|---|---|
| Adaptive Functional Programming | `10.1145/1186632.1186634` |
| An Experimental Analysis of Self-Adjusting Computation | `10.1145/1596527.1596530` |
| Traceable Data Types for Self-Adjusting Computation | `10.1145/1806596.1806650` |
| Imperative Self-Adjusting Computation | `10.1145/1328438.1328476` |
| Selective Memoization | `10.1145/604131.604133` |
| Adapton | `10.1145/2594291.2594324` |
| Nominal Adapton | `10.1145/2814270.2814305` |

**Kinetic data structures and persistence classics** (pre-arXiv era):

| Paper | DOI |
|---|---|
| Basch, Guibas & Hershberger, **Data Structures for Mobile Data** | `10.1006/jagm.1998.0988` |
| Kinetic Algorithms via Self-Adjusting Computation | `10.1007/11841036_57` |
| Robust Kinetic Convex Hulls in 3D | `10.1007/978-3-540-87744-8_3` |
| Fiat & Kaplan, persistence | `10.1016/s0196-6774(03)00044-0` |
| Driscoll, Sarnak, Sleator & Tarjan, **Making Data Structures Persistent** | `10.1016/0022-0000(89)90034-2` |
| Isenburg et al., Streaming Delaunay | `10.1145/1141911.1141992` |
| T-BON, temporal isosurface indexing | `10.1109/visual.1999.809879` |
| Guibas, KDS survey | `10.1201/9781420035179.ch23` — **book chapter, no PDF located anywhere** |

**R-functions — the weakest area in the whole corpus.** Rvachev, Shapiro and Pasko are all paywalled
and no specific DOIs were captured. Only Fryazinov 2010 and *Hybrid F-rep* cover F-rep at all. If
constructive implicit modelling matters, this is the gap to fill.

---

## 4. Books — this pipeline cannot fetch any of them

| Book | Why it matters |
|---|---|
| Cannarsa & Sinestrari, *Semiconcave Functions, HJ Equations, and Optimal Control* (Birkhäuser 2004) | **The** canonical semiconcavity / singular-set reference |
| Bardi & Capuzzo-Dolcetta (Birkhäuser 1997) | Viscosity solutions |
| Lions, *Generalized Solutions of Hamilton–Jacobi Equations* (Pitman 1982) | |
| Federer, *Geometric Measure Theory* (Springer 1969) | Sets of positive reach, in full |
| Dey, *Curve and Surface Reconstruction* (CUP 2006) | Sampling guarantees |
| Delfour & Zolésio, *Shapes and Geometries* (SIAM) | Distance functions as the primary object |

---

## 5. No verified DOI — **do not download until one is confirmed**

Known by title and author only. Several of these are more valuable than items in §3.

| Paper | Why it matters | Likely route |
|---|---|---|
| **Tarasov & Vyalyi 1998**, *Construction of contour trees in 3D in O(n log n) steps* | One of **two** results the entire 3D contour-tree question hinges on | Surfaced from Agarwal's reference list |
| **Safa & Wang 2014**, *Maintaining persistence and contour trees for time-varying functions on 2- or 3-manifolds* | The other one | OSU tech report, no DOI |
| **Nielson & Hamann**, *The Asymptotic Decider* (Vis '91) | Primary source for A-002's face rule — we implemented from a secondary description | IEEE Vis proceedings |
| **Nielson**, *On Marching Cubes* (TVCG 2003) | The definitive case-table treatment | |
| **Smith, Levien & Owens**, decoupled fallback / single-pass scan | Single-pass GPU scan on **exactly this platform** (Apple M-series, WebGPU) | Raph Levien's blog or the `vello` repo |
| **Sorensen et al.**, GPU forward-progress models | Finds **Apple and ARM GPUs do not support the linear occupancy-bound model** — a hard constraint on any persistent-threads design on Metal | OOPSLA / CONCUR proceedings |
| **Sweldens et al.**, lifting scheme / integer wavelets (4 canonical papers) | Entire subtopic returned empty; all four behind Elsevier/SIAM | |

> **Why `paper_search` didn't find the last three:** its `provider: all` **CORE ranking is broken** —
> two unrelated queries both returned 19th-century Texas land-title deeds, including the same Fort
> Worth and Holland, Michigan records. Restricting to `provider: crossref` fixes it immediately. Also
> `provider: arxiv` behaves as **exact-phrase match**, so multi-word conceptual queries silently
> return `[]`.

---

## 6. "Have" but unusable — recoverable for free, from data already held

`paper_download` fetched the **HTML landing page** instead of the PDF. In several cases the catalog
entry's own `download_urls` field already contains a working direct PDF link — **a pipeline bug, not
an acquisition problem.**

| Stem | Paper | Direct PDF |
|---|---|---|
| `10.1090_s0002-9947-1983-0690039-8` | Crandall & Lions | AMS archive |
| `10.1090_s0002-9947-1984-0732102-x` | Crandall, Evans & Lions (10 junk chunks) | AMS archive |
| `10.1090_s0025-5718-07-01981-3` | Cao (7 junk chunks) | AMS archive |
| `10.1090_s0025-5718-04-01678-3` | Zhao, fast sweeping | AMS archive |
| `10.1137_060670298` | Jeong & Whitaker — **also an unverified DOI, guessed from memory; 1 KB stub** | verify first |
| `10.14733_cadconfp.2022.329-333` | CAD'22 **table of contents** (6 chunks) | — |
| `10.1090_s0025-5718-07-01959-x` | Stevenson, NVB closure bound | `ams.org/mcom/2008-77-261/S0025-5718-07-01959-X/…pdf` |
| `10.1007_978-3-642-24550-3_29` | **Shapiro et al., the canonical CRDT paper** — 3.6 KB, 1 chunk | `core.ac.uk/download/49967317.pdf` |
| `10.1111_cgf.12596` | De Floriani, *Morse Complexes for Shape Segmentation* — HTML abstract only, 3 chunks | — |
| `10.1007_s00371-007-0163-2` | 6 chunks of navigation text | check `download_urls` |
| `10.1137_060675666` | 1 chunk | check `download_urls` |
| `10.1006_acha.1997.0238` | landing page | check `download_urls` |

**Detector:** `pdf_path` ending in `.html` with a low chunk count. Two independent producers of this
signature are now confirmed (generic landing pages and AMS specifically), and `pipeline_drift` stands
at **80** against a threshold of **3** — so there are almost certainly more than these twelve.

---

## Counts

| Category | Items |
|---|---|
| Free, resolver failed | 14 |
| Paywalled, high value | 24 |
| Paywalled, lower priority | 15 |
| Books | 6 |
| No verified DOI | 7 |
| Held but unusable | 12 |
| **Total** | **~78** |

## Order I'd work it

1. **§6 first.** Costs nothing, recovers papers already paid for, and one of them is the canonical CRDT
   paper currently sitting as a 1-chunk stub. Fix the landing-page bug once and re-run the detector
   across the whole catalog.
2. **§1 next.** Fourteen papers, all free, all blocked by a resolver limitation rather than a publisher.
   Includes Sethian, Zhao, Crandall & Lions and Federer — the foundations of everything in Phase 12.
3. **§2a — the five Acar meshing papers.** Author homepages. Highest intellectual risk on the list.
4. **§2b's Bálint** — `10.14733/cadaps.2023.1154-1174` — because F-001 is the ticket that gates
   Phases 10–13 and this is the paper that specifies its type.
