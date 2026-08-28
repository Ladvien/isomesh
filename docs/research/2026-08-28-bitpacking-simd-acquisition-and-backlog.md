# Bit-packing and SIMD: what was acquired, and twenty candidates from it

**Written:** 2026-08-28 · **Acquired:** 17 papers into home-still, verified against `M-371`'s three discriminators · **Corpus before:** 9,520 documents / 290,647 chunks.

**Why this sweep.** `✗51` closed the autovectorisation direction with a null and a rule — `%ymm` was zero in every one of eleven monomorphisations, the loop shape cost 3–7%, and the clause was unreachable anyway because the sample loop is 11.5% of extraction. That result is correct and it closed the wrong door. **It tested compiler autovectorisation of float code. It did not test bit-parallel integer work, which is where this crate's structure actually is** — `P-40`'s active-cell bitmap is already 64 cells per `u64` and already worth 5.5×, and nothing has asked what else fits in a word.

The corpus had four hard holes on exactly that question. This document reports the acquisition, and turns it into twenty candidates.

---

## Part 1 — What the corpus was missing, and what it now has

A survey across twelve topics found **four categories entirely absent**: broadword/SWAR, compressed and hierarchical bitmaps, SIMD integer codecs and bit-packing, and stream compaction without hardware `PEXT`. No Vigna, no poppy, no Gog & Petri, no Roaring, no FastPFOR, no Stream VByte. The corpus's entire succinct-data-structures holding was one external-memory graph-traversal thesis. What it did have — VDB (`10.1145/2487228.2487235`, 147k chars), sparse voxel DAGs (`10.1145/2461912.2462024`), Shewchuk (`10.1007/pl00009321`) — is the voxel and predicate side, not the word side.

### Acquired, all seventeen verified

Every row below passes all three discriminators: markdown far above the ~2,000-character landing-page floor, `chunks_indexed` never 1, and `pdf_path` ending in `.pdf`. Each was additionally confirmed with `distill_exists`.

| DOI | Paper | md chars | chunks |
|---|---|---|---|
| `10.14778/3598581.3598587` | Afroozeh & Boncz, **The FastLanes Compression Layout**, PVLDB 2023 | 69,236 | 22 |
| `10.48550/arXiv.1209.2137` | Lemire & Boytsov, **Decoding billions of integers per second through vectorization** | 122,111 | 38 |
| `10.1016/j.ipl.2017.09.011` | Lemire, Kurz & Rupp, **Stream VByte**, IPL 2018 | 26,085 | 9 |
| `10.1093/comjnl/bxx046` | Muła, Kurz & Lemire, **Faster Population Counts Using AVX2**, Comput. J. 2017 | 46,837 | 15 |
| `10.1007/s00778-019-00578-5` | Langdale & Lemire, **Parsing gigabytes of JSON per second**, VLDBJ 2019 | 99,962 | 31 |
| `10.1145/1572769.1572795` | Billeter, Olsson & Assarsson, **Efficient stream compaction on wide SIMD**, HPG 2009 | 48,658 | 16 |
| `10.48550/arXiv.1206.4300` | Vigna, **Quasi-succinct indices** (Elias–Fano), WSDM 2013 | 75,128 | 24 |
| `10.48550/arxiv.1301.5468` | Vigna, **Broadword Implementation of Parenthesis Queries** | 27,872 | 9 |
| `10.1007/978-3-031-20643-6_19` | Kurpicz, **Engineering Compact Data Structures for Rank and Select** (pasta), SPIRE 2022 | 41,829 | 14 |
| `10.1016/j.is.2021.101756` | Pibiri & Kanda, **Rank/select queries over mutable bitmaps**, Inf. Syst. 2021 | 67,403 | 22 |
| `10.1002/spe.2402` | Lemire, Ssi-Yan-Kai & Kaser, **Consistently faster and smaller compressed bitmaps with Roaring**, SPE 2016 | 107,472 | 33 |
| `10.1002/spe.2325` | Chambi, Lemire, Kaser & Godin, **Better bitmap performance with Roaring bitmaps**, SPE 2015 | 37,249 | 13 |
| `10.1145/3318464.3380588` | Lang et al., **Tree-Encoded Bitmaps**, SIGMOD 2020 | 84,626 | 26 |
| `10.1109/tpami.2021.3055337` | Bolelli et al., **One DAG to Rule Them All** (GRAPHGEN), TPAMI 2021 | 66,962 | 21 |
| `10.1117/12.596105` | Wu, Otoo & Shoshani, **Optimizing connected component labeling algorithms**, SPIE 2005 | 36,747 | 12 |
| `10.48550/arXiv.1505.05571` | Neal, **Fast exact summation using small and large superaccumulators** | 63,572 | 20 |
| `10.48550/arXiv.2401.14906` | Schroeder et al., **A High-Performance SurfaceNets Discrete Isocontouring Algorithm** | 60,007 | 19 |

**All four holes are now closed.** Broadword has Vigna's parenthesis-queries paper (the same word-level machinery, a different query) plus pasta and Pibiri; bitmaps have both Roaring papers and Tree-Encoded Bitmaps; codecs have FastLanes, Lemire & Boytsov and Stream VByte; compaction has simdjson and Billeter.

### Not acquired, and the reason is the pipeline rather than the corpus

Eight targets have **open copies that `paper_download` cannot reach**. This is worth recording as its own finding, because it is a second way for an acquisition to look impossible when it is not — a sibling of `M-371`.

| Paper | Where the open copy actually is | Why the resolver missed it |
|---|---|---|
| Vigna, **Broadword Implementation of Rank/Select Queries** (the highest-value target in the sweep) | `vigna.di.unimi.it/ftp/papers/Broadword.pdf` — fetched and confirmed to be the right document | **It has no DOI at all.** The only other identifier is ACM's `10.5555/1788888.1788900`, which no provider resolves |
| Zhou, Andersen & Kaminsky, **poppy / cs-poppy** | `ndownloader.figshare.com/files/12101855`, 302-redirecting to a live PDF; CMU KiltHub DOI `10.1184/R1/6609722.v1` is real and resolves in OpenAlex | `paper_download` traverses arXiv and Unpaywall, **not DataCite or figshare-hosted files** |
| Inoue et al., **set intersection** | `vldb.org/pvldb/vol8/p293-inoue.pdf` | The OA location is listed under **the same DOI that failed** |
| Gog & Petri | `core.ac.uk/download/342989546.pdf` | Same |
| Demmel & Nguyen | `eecs.berkeley.edu/~hdnguyen/public/papers/ARITH21_Fast_Sum.pdf` | Same |
| Flying Edges, HistoPyramids, Fujita | Kitware / MPI-Inf / publisher copies; no verified alternate identifier found | — |

**The lever for three of them is not a better DOI — it is a fetch-by-URL path**, because the metadata already lists an OA location the resolver does not follow. Four traps were identified and deliberately not taken: `arXiv:1311.1249` is not Gog & Petri's paper; Berkeley EECS-2015-229 is not Demmel & Nguyen's ARITH paper; the Hiroshima dissertation is not Fujita's; and the MPI-Inf report *HistoPyramids in Iso-Surface Extraction* is not Dyken's CGF paper. Rule 42 held.

**Two operational facts about the pipeline, both learned the expensive way and both worth a Part 5 rule.** `scribe_convert` is capped at 60 seconds by the MCP client while the job **keeps running server-side** for 82–768 seconds — so a timeout is not a failure, and the only reliable completion check is `catalog_read`. And the olmocr backend fails transiently under concurrency (`0 completed pages (failed=N)`, `workspace/markdown does not exist`) while `scribe_health` still reports `status: ok`; every one of those converted on a sequential retry. One paper timed out and then silently produced no markdown, which looked identical to a successful in-flight run — **`catalog_read` is the only thing that told them apart.**

---

## Part 2 — The strongest lead, and it directly diagnoses `✗51`

FastLanes' central result is not about SIMD. Verbatim from the paper, now in the corpus:

> *"We find it remarkable that SIMD-friendly ideas like interleaving and transposing accelerate our scalar code, rather than slow it down."*

> *"Scalar_T64 uses 64-bit scalar registers as quasi-SIMD and beats naive Scalar up to 8x."*

> *"clang++ can auto-vectorize our Scalar code, matching the performance of explicit intrinsics… when incorporating FastLanes in future systems, we recommend just using the Scalar code paths."*

Two things make this land here rather than being someone else's result. First, **`Scalar_T64` is `u64` shift, mask, AND and OR — safe, stable, `no_std`, integer, bit-exact on both targets.** It needs no `core::simd`, no intrinsics, no `unsafe`, and no target gating. It is exactly what `P-40` already did once by hand and never generalised.

Second, **the Apple M1 is in their hardware table** — ARM64, 128-bit NEON, 3.2 GHz — alongside Ice Lake, Zen3, Zen4 and two Gravitons, and they report that *"in terms of scalar performance, M1 tops Ice Lake clock-for-clock"* and that it *"clearly has more instruction level parallelism"* despite the narrow vector unit. Every other paper in this sweep is x86-only. This one measures both of this crate's machines.

**The diagnosis it offers for `✗51`.** `%ymm` was zero because the *layout* denied the compiler anything to widen, not because the compiler is weak — FastLanes' whole thesis is that the horizontal layout is the obstacle and the interleaved one is the fix, and that the fix pays in **scalar** code first. `✗51`'s null stands exactly as measured. What it did not test is a different loop over a different quantity: the sample loop is floats and 11.5% of extraction; the case-index computation is **bits**, and nothing has measured its share.

---

## Part 3 — Twenty candidates

Numbering continues from `P-102` / `R-102`. These are written as candidates rather than finished registrations — each carries the claim, the measured prior, the clause sketch and the falsifier, which is what makes turning them into `experiment.rs` entries mechanical. **Every one states the SHARE it can move, per `✗51`'s rule, and every one names the column that proves its fixture could have failed, per `P-70`'s C3.**

### Group A — the word as a vector

**`P-103` — the case index computed 64 cells at a time, from bit-sliced sign planes.** The 8 corner signs of a cell are 8 bits. Store the sign field as **eight bit-planes** rather than one array of bytes, and a whole `u64` of case indices falls out of 8 shift-and-OR pairs — `Scalar_T64` applied to the one quantity in this crate that is already boolean. `P-40` proved the packed representation is worth 5.5× for *deciding which cells*; this asks what it is worth for *deciding what each cell is*. **Share:** the case-classification stage, which has never been separately measured — C1 must measure that share before claiming a speedup, and if it is under 10% the row closes there. **Falsified by:** a share under 10%, or bit-sliced classification not beating the byte path by 2×, or any golden-hash movement (the case index is an integer and must be identical).

**`P-104` — the interleaved layout for the active-cell bitmap.** `build_inside_bits` packs 64 cells per `u64` **along X only** — FastLanes' "horizontal" layout, the one they measure as the obstacle. `M-287` found one bit of the row length was a 3.4× tax; `✗28` found the 128³ penalty is the access pattern rather than the stride. Both are layout results and neither tried a transposed layout. **Falsified by:** no improvement on the `M-287` fixture, which would mean the tax is genuinely the stride and `✗28`'s conclusion needs revisiting.

**`P-105` — Harley–Seal carry-save popcount for the per-chunk active count.** `M-349` established that the bitmap's claim was always a **count**, so the count is the quantity that matters. Harley–Seal turns *n* popcounts into ~*n*/8 using only AND, XOR and OR — safe scalar, and the reason the AVX2 paper's kernels are irrelevant while its reduction is not. **Share:** the counting pass only. **Falsified by:** under 2× on a 64³ chunk, or the count differing by one bit from the naive path.

**`P-106` — SWAR sign extraction and edge-crossing masks.** The twelve cut-edge flags of a cell are derived from the 8 sign bits by a fixed boolean circuit. Fujita's *Bitwise Parallel Bulk Computation* is the model — evolve 64 lanes at once through a hand-derived circuit — and it reports **13.4 × 10⁹ cell updates/s on an Intel Core i7 using the CPU bit-parallel technique**, not just the GPU figure. The paper did not download; the technique is documented well enough in the acquired FastLanes and popcount papers to register against. **Falsified by:** a circuit longer than 24 word operations, which is where the byte-table path wins on instruction count alone.

### Group B — rank instead of prefix sum

**`P-107` — a rank directory over the active-cell bitmap gives the output slot index in O(1).** This is the composite idea of the whole sweep. Flying Edges spends a **prefix-sum pass** to turn per-row counts into output offsets; `M-150` moved that pass to the GPU for 1.56×; `✗54` then measured the scan at **4.37% of `gpu_total_ms`**, which is why the subgroup version was not worth landing. A two-level rank directory makes the offset a lookup instead of a pass. Kurpicz's pasta gives the space figure — **3.51% overhead**, cross-confirmed in two independent papers — which on a 64³ chunk's 32 KiB bitmap is about 1.1 KiB. **Share:** the compaction/offset stage on the CPU path, which `M-135`'s stage breakdown can bound before the harness is written. **Falsified by:** the offset stage being under 5% of extraction, which makes this Amdahl-dead the way `✗51` was — *check this first*.

**`P-108` — broadword select to walk the set bits, with no `PEXT` and no table.** The dossier foreclosed `PEXT`/`PDEP` for three independent reasons and named `P-40`'s set-bit walk as the safe substitute. Vigna's broadword select is a third option: multiply, shift, mask and `count_ones`, branchless, no lookup. The honest price is on record — Pandey, Bender & Johnson measured the `PDEP`+`TZCNT` select at **2–4× the broadword one on Haswell**, which is what the constraint costs and is worth quoting rather than hiding. **Falsified by:** not beating the current walk on a bitmap that is 97% zeros, where the dossier already argued the cost is dominated by words skipped entirely.

**`P-109` — Elias–Fano for the edge→vertex map, and it is aimed at `R-027`.** `M-314` found the computation after a local edit is edit-proportional and the output buffer is not: **56–77% of vertex slots change for an edit touching 0.038% of cells**, with the ratio growing `O(n²)`. `M-318` put the ceiling at **15,706 → 346 at 129³, a 45×**, and recorded that index-is-edge-id costs **230× memory**, which is what killed the obvious encoding. Vigna's quasi-succinct indices are the encoding that does not: a monotone sequence in ~`2 + ⌈log(u/n)⌉` bits with O(1) access, no hashing, and no global build. Edge indices within a chunk are monotone once swept in order. **Falsified by:** the encoded map exceeding 4 bits per crossing, or access costing more than the current direct addressing by 1.5× — `R-027a` should run first and say where the 45× actually goes.

**`P-110` — mutable rank/select for a structure written during extraction.** `R-027` was **declined and split** on `V-45`, because its only working shape is output that depends on prior state, which converts a shipped determinism gate's failure condition into its intent. Pibiri & Kanda's contribution is rank/select that supports `flip(i)` while staying fast — the shape a cache written *during* a sweep needs rather than one built once. **Falsified by:** any dependence on insertion order in the final structure, which is `V-45`'s objection arriving again and closes the row immediately.

### Group C — compaction without hardware help

**`P-111` — table-driven scalar compaction, 8 cells per lookup, branchless.** simdjson's mechanism is a 256-entry table mapping an 8-bit mask to a permutation. Reduced to scalar it is a `[[u8; 8]; 256]` const table giving the set-bit positions of a byte — **2 KiB, `const`-evaluable, safe, deterministic, and it sidesteps `_pext_u64` entirely.** This is the most actionable single idea in the compaction literature for a crate that cannot use intrinsics. **Falsified by:** not beating `P-40`'s set-bit walk, or the table costing more L1 than it saves on a 97%-zero bitmap.

**`P-112` — count, scan, scatter, and the argument against the middle phase.** Billeter's three-phase decomposition is the standard shape and phase 2 is a prefix sum — which `P-107` replaces with an O(1) rank. Acquire the comparison, not the algorithm. **Falsified by:** the scatter phase dominating, which would mean the offsets were never the cost.

### Group D — what shape the occupancy structure should be

**`P-113` — Roaring's density thresholds as the chunk-representation decision rule, against `M-306` and `M-377`.** Roaring switches container type at **4096 set values per 2¹⁶ chunk — 6.25% density** — and says outright that *"when applications encounter integer sets with lower density (less than 0.1%), a bitmap is unlikely to be the proper data structure."* Put that beside `M-306`: `gyroid`'s rejected share is **16.8%** against `thin_plate`'s **95.1%**, so `thin_plate`'s active density is around 4.9% — **directly on Roaring's array/bitmap boundary**, and `gyroid` is far above it. And `M-377` has just moved the optimum chunk to 4³, where a chunk is *one word*. **The question is whether a near-empty chunk should carry a bitmap at all**, and there is now a published threshold to test against rather than a guess. **Falsified by:** one representation winning on all eight fields, which makes it a constant rather than a decision.

**`P-114` — a hierarchical bitmap above the active-cell bitmap.** One bit per 64 cells gives a second level at 1/64th the size; two levels skip 4,096 cells per test. VDB's node bitmask plus `popcount` child offset is the same object one level up, and it is already in the corpus. **Falsified by:** under 1.5× on `thin_plate`, the sparsest field, where it must win if it wins anywhere.

**`P-115` — Tree-Encoded Bitmaps for a subblock-empty summary.** The only RLE-family scheme that **preserves random access**, which is the property `rank` needs and which WAH and EWAH destroy. Up to a third of the space. **Falsified by:** random access costing more than 2× the flat bitmap's, which is the whole reason to prefer it over WAH.

### Group E — generating the table instead of transcribing it

**`P-116` — GRAPHGEN's decision-table pipeline applied to the 256-entry case table.** This crate's deepest recurring risk is a mistranscribed case table — `CLAUDE.md` rule 5 names it, and `✗50` is the most recent instance, where a *sampled* bound became a release-build panic. GRAPHGEN takes a decision table over boolean conditions, synthesises the **optimal decision tree** by minimum average path length, compresses it into a **DRAG**, and emits code. A marching-cubes case table is a decision table over 8 boolean corner conditions. **Generated Rust is safe, stable and deterministic by construction, and `P-64`'s Kani harness already proves the properties the generated table must keep.** The two compose: generate, then prove. **Falsified by:** the generated form being slower than the current table lookup, or `P-64`'s four properties failing on it — the second would be the interesting outcome and would mean the synthesis is unsound rather than merely unhelpful.

### Group F — determinism, and one live risk nobody has checked

**`P-117` — FMA contraction as a latent golden-hash divergence.** This is the highest-priority row in the group and it is a **risk, not an optimisation**. Shewchuk's error-free transformations — `two_sum`, `two_product` — are **broken by FMA contraction**, and aarch64 fuses far more eagerly than x86-64. `P-68`'s harness already used `two_sum` in `i128` and asserted its exactness. `M-31`'s 216 golden hashes have held across both machines, so either nothing on the hot path is contraction-sensitive or the crate has been lucky. **Nobody has checked which.** **Falsified by:** finding no contraction-sensitive expression in `crates/isomesh/src/**`, which is the outcome to hope for and is worth an afternoon to establish rather than assume.

**`P-118` — Neal's superaccumulator for cross-cell float accumulation, and it is aimed at `M-177`.** `M-177` established that **reordering cannot buy negation equivariance and the obstruction is structural**, and `M-372` found the duals stay below 48 of 48 because they accumulate crossings in an order axis relabelling permutes. A superaccumulator is not a reordering — it is a fixed-point accumulator (67 chunks of 64 bits with 32-bit overlap so carries propagate rarely) that is **exactly order-independent**, which is a stronger property than "the same on two platforms". Neal measures it at **less than twice the time of simple ordered summation**. This is a second route to `P-101`'s question and it does not depend on finding an invariant key. **Share:** the dual vertex solve only — `M-25` puts the sharp-feature solve at 3% over Surface Nets, so the ceiling is small and C1 must say so. **Falsified by:** the duals still failing 48 of 48 with an exactly order-independent accumulator, which would prove `M-177`'s obstruction covers the octahedral case and closes both `P-101` and this row together.

**`P-119` — double-buffering as the determinism mechanism, from a published isocontouring algorithm.** Parallel SurfaceNets uses double-buffered smoothing **explicitly to guarantee determinism** under parallelism, and reports one to two orders of magnitude over sequential algorithms on commodity CPU hardware. This crate's determinism today rests on the loops being sequential — `P-9` already found that a chunked consumer appending chunks in a different order gets a **different vertex count** on `noise_cavity`, and that order is deterministic only because the loops are, with nothing enforcing it. **Falsified by:** double-buffering not removing `P-9`'s spread, which would mean the order-dependence is in the weld rather than the traversal.

**`P-120` — array-based union-find for per-chunk labelling.** Wu, Otoo & Shoshani replace pointer-based equivalence trees with a **flat array** and report **5×–100× on random binary images**. Pointerless, `alloc`-only, no `unsafe`, consecutive final labels, fully deterministic. The banked union-find over the air sublevel set is a *global* structure maintained incrementally; this is the per-chunk version, for dropping isolated fragments and for `P-84`'s island detection. **Falsified by:** under 2× against the current structure, or any label depending on iteration order.

### Group G — the two that are measurements rather than mechanisms

**`P-121` — what fraction of extraction is bit work.** Before any of Group A is built, decompose extraction into float work (sampling, interpolation, the solve) and integer work (case classification, edge masks, compaction, indexing). `✗51`'s rule requires it and `✗51` is itself the example of what happens without it: the loop was 11.5% and the clause asked for a 2× that was arithmetically unreachable. **This row gates Group A and should run first.** **Falsified by:** nothing — it is a measurement, and a total under 15% closes `P-103` through `P-106` on the spot, which is a cheap and valuable outcome.

**`P-122` — Stream VByte's control/data split, applied to the case stream.** Keep the per-cell case index in one dense, branch-predictable stream and the variable-length vertex payload in another, so the classifier never branches on payload. Measured on named hardware — **over 4 billion differentially-coded integers per second from RAM to L1 on a 3.4 GHz Haswell** — though that figure is the SIMD path and only the layout transfers. **Falsified by:** no branch-misprediction reduction, measured with counters on Linux (`P-12`'s existing `perf_event_open` harness), which is the only instrument that can see it.

---

## Part 4 — Ordering, and what this does not claim

**`P-121` runs first.** It is a day, it is a measurement, and it can close six of the other rows before anyone writes a harness. `✗51` is the reason that sentence is here.

Then `P-107` (which has its own share check built into C1 and may close immediately), `P-117` (a risk, not a gain), `P-111` and `P-105` (both small), `P-103` (the largest of Group A, gated on `P-121`), `P-113` (a decision rule with a published threshold and existing fixtures), then the rest.

**What this does not claim.** Every performance figure in this sweep except FastLanes' is x86-64, and **not one paper reports bit-identical results across arm64 and x86-64** — treat all of them as upper bounds on the M5 until re-measured here. Roughly half the rank/select and connected-components literature publishes abstracts asserting superiority without ratios; where a number is missing above it is because the paper does not state one in what was retrieved, not because it was omitted. And the four foreclosures that shaped this list stand unchanged: `PEXT`/`PDEP`, Morton for the sample grid, minimal perfect hashing for the edge cache, and cache-oblivious layouts. Two more are added by this sweep — **WAH and EWAH**, because run-length encoding destroys the O(1) random access `rank` needs and the crate's 5.5× came from dense packing; and **`la_vector`'s learned rank/select**, which needs the key set up front and falls to the same objection that already killed the MPHF.
