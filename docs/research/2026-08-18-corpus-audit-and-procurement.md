# Corpus audit and procurement list

**Date:** 2026-08-18
**Method:** eight parallel domain sweeps over the home-still corpus (9,425 documents, 285,227 embedded chunks), then direct verification of every contradiction the sweeps produced against `catalog_read`, `paper_get` and `paper_download`.
**Tier:** the coverage tables are **M** — measured against the live corpus this session. The gap list is **V** where a DOI was resolved and **R** where it rests on a sweep agent's reading.

---

## 1. The directive's premise is false, and the true defect is worse

The directive opens: *"Audit why SDF papers are missing from home-still."*

They are not missing. Nor are the meshing papers. The corpus already holds the Jones/Bærentzen/Šrámek distance-field survey, Bærentzen's angle-weighted pseudonormal, DeepSDF, *CSG on Neural Signed Distance Fields*, *Reach For the Spheres*, *Interactive Editing of Voxel-Based SDFs*, *1-Lipschitz Neural Distance Fields*, Sethian & Vladimirsky on the eikonal equation, and *Power Diagram Enhanced Adaptive Isosurface Extraction from SDFs*.

**The contradiction is the finding, and it generalises.** Eight independent agents swept the corpus. Between them they reported roughly a dozen foundational papers ABSENT. Three of those reports were checkable in seconds, and all three were wrong:

| Paper | Reported | Actual |
|---|---|---|
| Lorensen & Cline, *Marching Cubes*, 1987 | ABSENT | On disk, converted, **15 chunks embedded** |
| van Gelder & Wilhelms, *Topological Considerations in Isosurface Generation*, 1994 | ABSENT | On disk, converted, **32 chunks embedded** |
| Kobbelt et al., *Feature Sensitive Surface Extraction*, 2001 | ABSENT | On disk, converted, **31 chunks embedded** — and indexed **twice**, also as `paper_p57-kobbelt` |

Each of the three ranks **#1** on a targeted query. The agents' searches *did* return them. Here is what a search hit for the Marching Cubes paper actually looks like:

```json
{ "doc_id": "10.1145_37401.37422", "title": null, "authors": [], "year": null,
  "doi": null, "score": 0.704, "pdf_path": "papers/10/10.1145_37401.37422.pdf" }
```

An agent cannot identify that, so it discards it and reports the paper missing.

**The corpus is not short of papers. It is short of the metadata that makes them recognisable.**

### 1.1 Quantified

Twelve diverse queries, `limit=40`, 475 raw hits:

| Metric | Count | Rate |
|---|---|---|
| Hits with `title: null` | 177 | **37.3%** |
| Hits with `doi: null` | ~166 | 34.9% |
| Hits with empty `authors` | ~177 | 37.3% |

Per-domain the rate ranges from 2% (contour trees) to **88%** (GPU compute shaders) — it tracks stem shape, not topic:

| Stem type | Hits | null title | Rate |
|---|---|---|---|
| Non-DOI stems (`s2010-advances-*`, `eth-cgl-*`, `paper_p57-kobbelt`, `scanpaper_final`) | 105 | 105 | **100%** |
| DOI-shaped stems (`10.xxxx_...`) | 370 | 72 | 19.5% |

### 1.2 Where the loss is

`catalog_read` on twelve broken stems versus two known-good controls settles it. Every broken row carries **only** `conversion` and `embedding` blocks. The whole download-provenance block — `title`, `authors`, `doi`, `downloaded_at`, `sha256`, `file_size_bytes`, `source`, `pdf_path` — is absent. The controls carry all of it.

**So the metadata never existed in the catalog; it is not stale in Qdrant.** These PDFs entered by a path that skips the provider fan-in `paper_download` performs — a manual drop, an inbox sweep, or a lookup that failed silently — and were then converted and embedded anyway. `distill_search` is reporting nulls honestly.

This also explains why `catalog_backfill_title` reports only **2** candidates: it selects rows that have a DOI and lack a title, and these rows have no DOI field at all.

### 1.3 The years are not merely missing, they are wrong

Ten hits carrying a non-null year were spot-checked. **Seven to eight are wrong**, several impossibly so:

| Stem | Reported | Actual |
|---|---|---|
| `10.1016_j.cag.2006.07.021` (Newman & Yi survey) | 1987 | 2006 |
| `10.1145_195826.195828` (van Gelder & Wilhelms) | 1980 | 1994 |
| `10.1145_1281500.1281668` | 1942 | ACM did not exist |
| `s2007-advances-course-notes-1-6-mb` | 1942 | 2007, per its own stem |
| `Bounded-Biharmonic-Weights-...` | 1968 | 2011 |
| `labelle_shewchuk_isosurface_stuffing_2007` | 1976 | 2007, per its own filename |
| `eth-cgl-sim_anim-Pfa12` | 1931 | ~2012 |

Pre-1950 dates on SIGGRAPH papers mean something other than a publication year is being written into that field. A year filter (`distill_search`'s `year: ">2020"`) applied to this corpus silently discards correct papers.

### 1.4 Pipeline health, from `catalog_repair` and `system_status`

| Condition | Count | Meaning |
|---|---|---|
| `pipeline_drift` | **79** | threshold is **3** |
| `stuck_convert` | **393** (256 pdf + 137 html) | downloaded, never converted, therefore invisible to search |
| `embedding_skipped` | 515 | converted, never embedded |
| `corrupted_pdfs` | 60 | |
| `flag_drift` | 6,672 | |
| `md_path_drift` | 188 | |
| `catalog_no_source` | 24 | phantom rows — **includes `10.1145_3197517.3201353`, TetWild**, which a sweep listed as a gap and is really a row that lost its file |
| `catalog_no_markdown` | 1 | |

Also observed: arXiv doc_ids are **case-split** — `10.48550_arXiv.2308.05371` and `10.48550_arxiv.2308.05371` are two Qdrant entries for one paper, with inconsistent metadata between them. Confirmed for FlexiCubes, the approximate-convex-decomposition paper, and the hex-mesh survey.

### 1.5 The fix, in order

1. **Backfill catalog metadata for rows that have only `conversion`+`embedding`.** For DOI-shaped stems the stem *is* the DOI — a provider lookup keyed off it fills everything. `catalog_backfill_title` cannot currently see these rows because it filters on a DOI field they lack; it needs to fall back to parsing the stem. This is the highest-value repair and it is mechanical.
2. **For non-DOI stems (100% broken), extract the title from the first page of the markdown and fuzzy-match against OpenAlex.** Roughly 60–100 unique documents.
3. **Clear the 393 stuck conversions.** These are papers you have already paid to acquire and cannot currently find.
4. **Normalise arXiv doc_id casing before re-indexing**, or the backfill runs twice per paper.
5. **Re-index after backfill** so payloads carry the repaired metadata.

Until (1) and (2) are done, **any literature review over this corpus — by a person or an agent — will produce false negatives at roughly a one-in-three rate.** That is the pipeline gap the directive asked for. It is not about SDFs.

---

## 2. Coverage by area

Verified present unless marked. Nothing in this table rests on an agent's ABSENT claim alone.

### Primal extraction — strong

Lorensen & Cline 1987 (`10.1145/37401.37422`) · van Gelder & Wilhelms 1994 (`10.1145/195826.195828`) · Nielson, *On Marching Cubes* 2003 (`10.1109/tvcg.2003.1207437`) · Montani/Scateni/Scopigno MC33 lookup table 1994 (`10.1007/bf01900830`) · Lewiner et al. 2003 (`10.1080/10867651.2003.10487582`) · Custodio et al. 2013 (`10.1016/j.cag.2013.04.004`) and 2019 (`10.1186/s13173-019-0086-6`) · Natarajan 1994 (`10.1007/bf01900699`) · Newman & Yi survey 2006 (`10.1016/j.cag.2006.07.021`) · Grosso 2016 (`10.1111/cgf.12975`) and 2017 (`10.1145/3095140.3095179`) · Kirby et al., verifiable visualization 2009 (`10.1109/tvcg.2009.194`) · Dietrich et al., edge transformations 2008 (`10.1109/tvcg.2008.60`) · Shu et al., adaptive MC 1995 (`10.1109/2945.485620`).

**Thin:** classical marching-tetrahedra decomposition-bias literature. **Genuinely absent and pre-DOI:** Nielson & Hamann's *The Asymptotic Decider* (Vis '91), Chernyaev's MC33 CERN report (1995), Doi & Koide (IEICE 1991).

### Dual extraction — strong

Ju/Losasso/Schaefer/Warren, *Dual Contouring of Hermite Data* 2002 — **present at `10.1145/566570.566586`**. Three sweeps reported it absent using the DOI `10.1145/566654.566586`, which does not resolve; that DOI is wrong, not the corpus. Also: Schaefer & Warren, *The Secret Sauce* · Schaefer/Ju/Warren, *Manifold Dual Contouring* 2007 (`10.1109/tvcg.2007.1012`, indexed twice) · Nielson, *Dual Marching Cubes* 2004 (`10.1109/VISUAL.2004.28`) · Schaefer & Warren, *Primal Contouring of Dual Grids* 2004 · Gibson, *Constrained Elastic Surface Nets* 1998 (`10.1007/bfb0056277`) · Ho et al., *Cubical Marching Squares* 2005 (`10.1111/j.1467-8659.2005.00879.x`) · Manson & Schaefer 2010 (`10.1111/j.1467-8659.2009.01607.x`) · Ju, *Robust Repair of Polygonal Models* 2004 · Kobbelt et al. 2001 (`10.1145/383259.383265`).

**Absent:** Ju's *Intersection-Free Contouring on an Octree Grid* (2006), Zhang/Hong/Kaufman (2004), Kazhdan et al. *Unconstrained Isosurface Extraction on Arbitrary Octrees* (2007), Varadhan et al. (2003).

### Adaptive, LOD, chunking, seams — deep on transitions, absent on chunk mechanics

Lengyel's Transvoxel journal paper (`10.1080/2151237X.2011.563682`) **and his full dissertation**, 13 indexed passages · Frisken et al., ADF 2000 · Hoppe progressive meshes (`10.1145/237170.237216`), view-dependent refinement, terrain LOD (`10.1109/visual.1998.745282`) · Livnat/Shen/Johnson NOISE (`10.1109/2945.489388`) · Cignoni et al. interval trees (`10.1109/2945.597798`) · Museth VDB (`10.1145/2487228.2487235`).

**Absent:** chunked/ghost-cell voxel-world literature entirely; out-of-core isosurface extraction; temporal coherence across frames; spatial-hash grids; van Kreveld et al. seed sets (**now acquired**, §3); geometry clipmaps.

**The seam answer the corpus actually supports is transition cells** (Lengyel), not dual-grid crack-avoidance. Ju 2002's "requires no crack patching" claim is now checkable — the paper is present — but has not been read against this question. Until it is, the DC-is-crack-free-by-construction claim is **R at best, not V**.

### Vertex placement and quality — strong on metrics, was thin on primary sources

Knupp et al., Verdict library (`10.2172/901967`) · Grosso & Zint (`10.1007/s00371-021-02139-w`) · Jiao & Zha (`10.1145/1364901.1364924`) · Garland & Zhou (`10.1145/1061347.1061350`) · Lindström & Turk (`10.1145/353981.353995`) · Plantinga & Vegter (`10.1145/1057432.1057465`) · Boissonnat/Cohen-Steiner/Vegter (`10.1007/s00454-007-9011-4`) · Blu/Thévenaz/Unser (`10.1109/tip.2004.826093`) · FlexiCubes, NMC, NDC, TetWeave, MeshSDF.

**Absent, now sourced (§3):** Garland & Heckbert 1997 · Cignoni et al. *Metro* 1998 · Witkin & Heckbert 1994 · Du/Faber/Gunzburger CVT 1999.

### GPU and parallel — the weakest area in the corpus

Present: Dyken et al. HistoPyramids (`10.1111/j.1467-8659.2008.01182.x`) · Schmitz et al. (`10.1111/j.1467-8659.2010.01825.x`) · Cirne & Pedrini (`10.1007/s13173-012-0097-z`) · Newman et al. SIMD MC (`10.1016/j.cag.2003.12.008`) · Crassin GigaVoxels · Aokana · Laine & Karras software rasterisation.

**Absent:** mesh/task shaders entirely · wave and subgroup intrinsics entirely · roofline and arithmetic-intensity analysis entirely · occupancy/divergence in case-table kernels entirely · cache blocking and NUMA · deterministic parallel reduction.

**Consequence for the repo:** the corpus **cannot** answer whether GPU isosurfacing is bandwidth-, latency- or compute-bound, at any grid size. This project's own banked figures (0.22 ms fixed cost, ~33³ break-even, field evaluation dominating) are currently the only measurements with that specificity anywhere in reach, and **must not be attributed to the literature**.

### Incremental and dynamic — strong theory, one genuine void

Acar's self-adjusting computation line is fully indexed (five papers) · Agarwal et al., *Maintaining Contour Trees of Dynamic Terrains* (`10.48550/arXiv.1406.4005`) · Acar/Cotter/Hudson/Turkoglu dynamic well-spaced point sets and kinetic mesh refinement · Driscoll/Sarnak/Sleator/Tarjan persistence (`10.1016/0022-0000(89)90034-2`) · Nielsen & Museth dynamic tubular grid · six modern dynamic-connectivity papers.

**The void, and it is worth naming precisely.** Neither `distill_search` over 9,425 documents nor a live provider search found **any** published scheme for giving isosurface output elements identities that survive a local edit without carrying mutable state across calls. The persistence literature solves cheap access to *old versions* and is stateful by construction. The dynamic-connectivity literature maintains *aggregates*, not output-buffer slots. **R-027a is operating in genuinely unoccupied territory** — which is the same shape of finding as V-43, and it should be recorded with the same care.

### Learned and frontier — well covered

FlexiCubes · DMTet · Neural Marching Cubes · Neural Dual Contouring · MeshSDF · nvdiffrec · VolSDF · TetWeave · fTetWild · TetGen · Labelle & Shewchuk isosurface stuffing · Ruppert · Boissonnat & Oudot · restricted power diagrams on GPU (`10.1111/cgf.142610`) · Power-diagram adaptive extraction (`10.48550/arXiv.2506.09579`) · 3DGS/2DGS/SuGaR/MILo · Subgrid Marching Tetrahedra (Baktash/Gillespie/Crane, arXiv 2606.00454 — **the paper this crate implements, present**).

---

## 3. Procurement list

### 3.1 Acquired this session

| Paper | DOI | Status |
|---|---|---|
| van Kreveld et al., *Contour Trees and Small Seed Sets for Isosurface Traversal* | `10.1145/262839.269238` | downloaded + converted |
| Lévy, *A Numerical Algorithm for L2 Semi-Discrete Optimal Transport in 3D* | `10.1051/m2an/2015055` | downloaded, conversion running |
| Aila & Laine, *Understanding the Efficiency of Ray Traversal on GPUs* | `10.1145/1572769.1572792` | downloaded, conversion running |

### 3.2 Paywalled by DOI, free copy located — needs a URL-based ingest

`paper_download` resolves DOIs only, so these six need `hs paper add <url>` or an inbox drop. Every URL below was fetched and content-verified.

| Paper | Free copy | Source |
|---|---|---|
| Garland & Heckbert, *Surface Simplification Using Quadric Error Metrics*, 1997 | `https://www.cs.cmu.edu/~garland/Papers/quadrics.pdf` | author page (CMU) |
| Treece/Prager/Gee, *Regularised Marching Tetrahedra*, 1999 | `http://mi.eng.cam.ac.uk/reports/svr-ftp/treece_tr333.pdf` | Cambridge CUED TR333 |
| Holm/de Lichtenberg/Thorup, JACM 2001 | `https://di.ku.dk/.../97-17.pdf` (Part I) and `.../97-26.pdf` (Part II) | author's institutional repo |
| Acar et al., *Parallel Batch-Dynamic Trees via Change Propagation*, ESA 2020 | `https://arxiv.org/pdf/2002.05129` | arXiv |
| Williams/Waterman/Patterson, *Roofline*, 2009 | `https://www2.eecs.berkeley.edu/Pubs/TechRpts/2008/EECS-2008-134.pdf` | Berkeley TR UCB/EECS-2008-134 |
| Hu et al., *Tetrahedral Meshing in the Wild*, 2018 | `https://www.cs.toronto.edu/~jacobson/images/tetrahedral-meshing-in-the-wild-siggraph-2018-compressed-hu-et-al.pdf` | author page (Jacobson) |

Two caveats carried honestly. The HDT JACM merged version is not posted anywhere free; the two DIKU technical reports are the same authors and jointly cover all four problems in the JACM title, so they are the legitimate equivalent, not the identical artefact. The Roofline TR is titled *"...for Floating-Point Programs and Multicore Architectures"* against the CACM's *"...for Multicore Architectures"* — same authors, same content.

This section is the working instance of **M-290's rule**: a paywalled DOI is not a missing paper. Six for six.

### 3.3 Not yet sought, ranked by leverage

1. **Ju, *Intersection-Free Contouring on an Octree Grid*, PG 2006.** A-025's open defect sits exactly in the gap between IFC (intersection-free, still non-manifold) and MDC (manifold, still self-intersecting). The corpus has one half of that pair.
2. **Kazhdan/Klein/Dalal/Hoppe, *Unconstrained Isosurface Extraction on Arbitrary Octrees*, SGP 2007.** Drops the restricted-octree assumption Ju 2002 carries.
3. **A mesh-shader primary source.** Nothing in the corpus covers the pipeline `isomesh-gpu` has already probed for.
4. **Zhang/Hong/Kaufman, DC with topology-preserving simplification, Vis 2004.** The named alternative to MDC's clustering rule.
5. **Tzeng/Patney/Owens, *Task Management for Irregular-Parallel Workloads on the GPU*, HPG 2010.** No resolvable DOI; OpenAlex W2016706026.
6. **Varadhan et al., *Feature-Sensitive Subdivision and Isosurface Reconstruction*, Vis 2003.** Thin features — the CAD consumer's case.
7. **Losasso/Fedkiw/Osher, spatially adaptive level sets, 2006** (`10.1016/j.compfluid.2005.01.006`). The nearest thing to temporal coherence for a dynamically changing field.

### 3.4 Pre-DOI, likely needs a human

Nielson & Hamann, *The Asymptotic Decider* (Vis '91) · Chernyaev, *Marching Cubes 33* (CERN CN/95-17) · Doi & Koide (IEICE Trans. E74-D, 1991). All three are primary sources for case tables this crate implements. Rule 5 says do not guess a case table; at present three of them are known only through others' paraphrase.

---

## 4. Two corrections owed to this repo's own documents

1. **`docs/research/` claims Ju 2002 is uncited/absent in places.** It is present at `10.1145/566570.566586`. Any ticket text carrying `10.1145/566654.566586` has a bad DOI.
2. **TetWild is a phantom catalog row**, not an absent paper — `catalog_no_source` lists `10.1145_3197517.3201353`. The file was lost after the row was written. A free copy is in §3.2.
