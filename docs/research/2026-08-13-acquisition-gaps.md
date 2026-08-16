# Acquisition gaps — what to get, and what we "have" that is unusable

**Date:** 2026-08-13
**Verified against home-still this session** via `catalog_read`, not from agent summaries.

---

## Part A — Broken but **recoverable from data we already hold**

This is the important finding, and it's better news than "polluted, re-source by hand."

`paper_download` fetched the **HTML landing page** instead of the PDF — but in each case the catalog
entry's own `download_urls` field already contains a working direct PDF link. **This is a pipeline
bug, not an acquisition problem.** Nothing needs hand-fetching; something needs to follow the link it
already has.

| Stem | Paper | What's there now | Direct PDF already in the catalog entry |
|---|---|---|---|
| `10.1090_s0025-5718-07-01959-x` | Stevenson, *Completion of locally refined simplicial partitions created by bisection* — **the NVB closure bound** | 11 MB of HTML → 18 KB of markdown, 7 chunks of MathML noise | `ams.org/mcom/2008-77-261/S0025-5718-07-01959-X/S0025-5718-07-01959-X.pdf` |
| `10.1007_978-3-642-24550-3_29` | Shapiro, Preguiça, Baquero, Zawirski, **the canonical CRDT paper** | 107 KB HTML → **3.6 KB, 1 chunk** — a stub | `core.ac.uk/download/49967317.pdf` and `inria.hal.science/hal-00932836` |
| `10.1007_s00371-007-0163-2` | *(per sweep — 6 chunks of nav text)* | landing page | check `download_urls` |
| `10.1137_060675666` | *(1 chunk)* | landing page | check `download_urls` |
| `10.1006_acha.1997.0238` | *(landing page)* | landing page | check `download_urls` |

**Action:** re-download from the `download_urls` PDF entry, re-convert, re-index. Then audit the whole
catalog for the same signature — **`pdf_path` ending in `.html` with a low chunk count** is the
detector. `pipeline_drift` is 80 against a threshold of 3, so there are likely more than these five.

---

## Part B — Must be obtained by hand

`paper_download` resolves **arXiv + Unpaywall only**. Every `10.1145/*` (ACM), `10.1109/*` (IEEE) and
most Elsevier/SIAM DOIs fail regardless of whether they're paywalled — LIPIcs and JoCG are fully open
access and still failed. All DOIs below were **verified via `paper_search` this session or by the
sweep agents**; none is reconstructed from memory.

Realistic route for the ACM ones: **author homepages.** Acar, Blelloch, Hudson and Türkoğlu all post
PDFs. That's usually faster than institutional access.

### Tier 1 — the closest existing prior art to the top-ranked transfer

Acar's own incremental **meshing** line. This is output-sensitive incremental meshing by the person
who built the change-propagation theory, and it is the single biggest hole in the corpus.

| Paper | DOI | Venue |
|---|---|---|
| **Dynamic Well-Spaced Point Sets** | `10.1145/1810959.1811011` | SoCG 2010 |
| ↳ *journal version, may be easier to reach* | `10.1016/j.comgeo.2012.11.007` | Comp. Geom. 2013 |
| **Kinetic Mesh Refinement in 2D** | `10.1145/1998196.1998254` | SoCG 2011 |
| Parallelism in Dynamic Well-Spaced Point Sets | `10.1145/1989493.1989498` | SPAA 2011 |
| **A Cost Semantics for Self-Adjusting Computation** — *the trace-distance paper* | `10.1145/1594834.1480907` | POPL 2009 |

### Tier 2 — self-adjusting computation foundations

| Paper | DOI |
|---|---|
| Adaptive Functional Programming | `10.1145/1186632.1186634` |
| An Experimental Analysis of Self-Adjusting Computation | `10.1145/1596527.1596530` |
| Traceable Data Types for Self-Adjusting Computation | `10.1145/1806596.1806650` |
| Imperative Self-Adjusting Computation | `10.1145/1328438.1328476` |
| Selective Memoization | `10.1145/604131.604133` |
| Adapton | `10.1145/2594291.2594324` |
| Nominal Adapton | `10.1145/2814270.2814305` |

### Tier 3 — kinetic data structures and persistence classics (pre-arXiv era)

| Paper | DOI |
|---|---|
| Basch, Guibas, Hershberger — Data Structures for Mobile Data | `10.1006/jagm.1998.0988` |
| Kinetic Algorithms via Self-Adjusting Computation | `10.1007/11841036_57` |
| Robust Kinetic Convex Hulls in 3D | `10.1007/978-3-540-87744-8_3` |
| Fiat & Kaplan — persistence | `10.1016/s0196-6774(03)00044-0` |
| Driscoll, Sarnak, Sleator, Tarjan — *Making Data Structures Persistent* | `10.1016/0022-0000(89)90034-2` |
| Isenburg et al. — Streaming Delaunay | `10.1145/1141911.1141992` |
| T-BON — temporal isosurface indexing | `10.1109/visual.1999.809879` |
| ESA 2020 (open access, still failed) | `10.4230/lipics.esa.2020.2` |
| JoCG (open access, still failed) | `10.20382/jocg.v3i1a11` |
| Guibas — KDS survey | `10.1201/9781420035179.ch23` **(book chapter, no PDF located anywhere)** |

---

## Part C — Wanted, but **do not download until the DOI is verified**

I know these by title and author. I do **not** have a verified identifier for any of them, and this
project's standing rule is never to guess one — an agent once guessed an arXiv ID and pulled an
unrelated condensed-matter physics paper into the corpus under a meshing DOI.

| Paper | Why it matters |
|---|---|
| **Nielson & Hamann — The Asymptotic Decider** (Vis '91) | The primary source for A-002's face ambiguity rule. We implemented from a secondary description |
| **Nielson — On Marching Cubes** (TVCG 2003) | The definitive case-table treatment |
| **Carr, Snoeyink, van de Panne — Flexible Isosurfaces** | Canonical "which level-set features to keep" — directly relevant to the persistence hypothesis |
| **Smith, Levien, Owens — decoupled fallback / single-pass scan** | Single-pass GPU scan on **exactly our platform** (Apple M-series, WebGPU) |
| **Sorensen et al. — GPU progress models** | Finds **Apple and ARM GPUs do not support the linear occupancy-bound model** — a hard constraint on any persistent-threads design on Metal |
| **Sweldens et al. — lifting scheme / integer wavelets** (4 canonical papers) | Entire subtopic returned empty; all four behind Elsevier/SIAM |

> ### Methodological warning — reproduced twice this session
> `paper_search` for *"The Asymptotic Decider … Nielson Hamann"* returned **three 19th-century land
> title abstracts**. A second, completely different query — *"forward progress models GPU occupancy
> bound portability Sorensen"* — returned **four more land title abstracts**, including the same
> Fort Worth and Holland, Michigan deeds. Two unrelated queries, one with no ambiguous words at all,
> both dominated by the same junk.
>
> **This is a broken relevance ranking in the CORE provider, not a keyword coincidence.** It fails
> loudly and absurdly rather than quietly, which is lucky — but it means any agent doing discovery
> through `paper_search` with the default `provider: all` is getting its result list polluted, and a
> less obviously wrong hit could be accepted.
>
> **Mitigations:** restrict `provider` to `crossref` or `openalex` for classic pre-arXiv work; use
> `search_type: "title"` with the exact title; never accept a top hit without reading its metadata.
> With `provider: crossref` the prefix-sum query returned four plausible, on-topic papers — so the
> restriction works.
>
> Neither Smith/Levien/Owens nor Sorensen et al. was located this session. **Their DOIs remain
> unverified and must not be guessed.** Likely faster routes: Raph Levien's blog/`vello` repository
> for the first, and the OOPSLA/CONCUR proceedings pages for the second.

---

## Part D — Corpus hygiene, from `system_status`

| Metric | Value | Note |
|---|---|---|
| `corrupted_pdfs` | **60** | Never enumerated. Worth listing and re-fetching |
| `embedding_skipped` | **508** | Documents converted but not searchable — invisible to `distill_search` |
| `pipeline_drift` | **80** vs threshold **3** | 27× over |
| Documents with `title: null` | many | Several sweep hits came back with no title or authors — `catalog_backfill_title` in batches of ~15 (it times out at 60) |

**The 508 skipped and the 342 previously-known invisible documents are the same class of problem**, and
they are why the standing rule exists: *presence in the corpus is decided by `catalog_read`, never by
`distill_search`.* That rule paid out again this session — **Probabilistic Quadrics**
(`10.1111/cgf.13933`), which supersedes the audit doc's `λ ≈ 0.01` regularizer, was already in the
corpus and invisible to search.

---

## Confirmed good this session

- **CoACD** — `10.48550/arXiv.2205.02961`, Wei, Liu, Ling, Su 2022. Downloaded 23 MB, converted in
  588 s, **29 chunks indexed**. This is the source for the 49% → 80% figure cited as V-15 without the
  paper. The `scribe_convert` MCP call timed out at 60 s; **the server finished asynchronously.** Do
  not retry on that error — wait and confirm with `catalog_read`.

---

## Order of operations

1. **Fix the HTML-landing-page bug and re-run Part A.** Five known, probably more. Costs nothing and
   recovers papers we already paid for — including the canonical CRDT paper, currently a 1-chunk stub.
2. **Get Tier 1 by hand** — five papers, author homepages. This is the prior art for the whole
   incremental-meshing hypothesis, and not having read it is the biggest risk to that claim.
3. **Enumerate the 60 corrupted and 508 skipped** before any further novelty search, because every
   "nobody has done this" conclusion is only as good as the corpus's searchability.
4. Tier 2 and 3 opportunistically.
