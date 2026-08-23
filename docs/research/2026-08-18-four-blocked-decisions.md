# The four blocked decisions

**Date:** 2026-08-18
**Why this document exists:** all 14 open tickets are blocked, and four of the blocks are yours. Three are genuine design decisions the ticket text explicitly refuses to make without you (`CLAUDE.md`: *"Do NOT assume design decisions on my behalf"*). The fourth is a sequencing call.

Each section gives the question, what is already measured, the options with their real costs, a recommendation, and — importantly — **what would change my mind**. The measurements are yours; the recommendations are mine and are argued, not asserted.

---

## Decision 1 — A-026: how far does the CoACD pipeline go in-crate?

### The question

Rule 3 forbids a new dependency without written justification, so every piece of a convex decomposer lands in `crates/isomesh` or does not exist. That makes "how much of it" a real choice rather than an implementation detail.

### What is already settled, and it is a lot

- **Route (a) is decided** — a decomposer that genuinely *partitions* the interior. ✗21 killed Convex Primitive Decomposition as a substrate: it covers the **surface**, its primitives **overlap by design** (30× cost), and *"guarantees enclosure"* is exactly what disqualifies it, because a wrapper is by definition bigger than what it wraps, so cutting it yields fat pieces.
- **CoACD over V-HACD, with merging OFF (V-36).** CoACD's §6.2 gives flat boundaries between components and intersection-free hulls without voxelisation. Its §6.5 merge post-process breaks that guarantee — a merged pair's hull can reach into a third neighbour — and **merging is on by default**. The merge exists only *"to further reduce the number of components,"* so switching it off costs component count, not correctness.
- **M-297: nothing published runs at interactive rates.** VisACD 16.97 s; CoACD 36.31 s and 194–253 s per model against a 16.67 ms frame. Every "real-time" method precomputes and says so. **This is not a per-frame feature under any option.**
- **CoACD assumes 2-manifold solid input** and says so — *"pre-processing with an off-the-shelf manifold conversion algorithm."* So for **this** method the repair pass ✗20 called a quality lever is a **prerequisite**.

### What each option actually costs

**(i) The whole pipeline.** 3D convex hull · surface *and* interior point sampling · Hausdorff between point sets · solid-mesh plane cut · MCTS over candidate planes · plus manifold repair. The ticket's own words: *"several tickets and the crate's largest single addition."* The convex hull alone is a ticket.

**(ii) A narrower first cut.** Convex hull plus a greedy plane split. No MCTS, no collision-aware concavity. Lower-quality cells that still partition, at a fraction of the work.

### Recommendation: **(ii), and split it before starting.**

Three reasons, in order of weight.

**A-027 is what actually needs unblocking, and it does not care.** The decomposer's interface is *mesh in, convex cells out* under either option. A-027 — where `plane ∩ convex polyhedron = convex polygon` makes the centroid fan **provably** correct, which is why the cuboid scores 8/8 and it is not luck — is owed regardless and is currently blocked behind a ticket whose scope is undecided. Option (ii) unblocks it in a fraction of the time.

**MCTS buys quality, and quality is not what is being tested.** A-026's own hypothesis is that fragments come out closed, manifold, χ=2 and volume-conserved to 1e-3. **Those are properties of partitioning, not of partition quality.** A greedy split that partitions satisfies the hypothesis; MCTS makes the cells prettier. Building the expensive half first is the "building the wrong thing carefully" failure that A-020b already paid for.

**It keeps the falsifier reachable.** The falsifier is *any proxy fragment reporting open cut edges*, which locates the defect in plane-cell intersection. That is testable the moment cells exist. Under option (i) it is testable after several tickets.

**Suggested split:** `A-026a` 3D convex hull (self-contained, generally useful, needed by both options) · `A-026b` greedy plane split producing partitioning cells · `A-026c` manifold-repair prerequisite, **scoped by measurement first** — check how many of your actual shells are already 2-manifold before building a repair pass for a problem you may not have · `A-026d` MCTS plane selection, **deferred until A-027 has measured whether cell quality is even the binding constraint**.

**What would change my mind:** if A-027's centroid-fan correctness turns out to be sensitive to cell *shape* and not just cell *convexity*, then greedy cells are not a valid substrate for the experiment and MCTS moves onto the critical path. That is checkable early — plane-cut a deliberately ugly greedy cell and check the cap.

**Keep on the desk, as the ticket says:** Diazzi & Attene's VolumeMesher (`10.1145/3478513.3480564`) reaches classified convex cells **without** tidying the input. Unusable as a dependency (C++, rule 3), but it is the measuring stick for what CoACD's mandatory repair pass costs in fidelity — otherwise invisible.

---

## Decision 2 — A-025: the `FaceAmbiguity` default

### The question

`ManifoldDualContouring` defaults to `FaceAmbiguity::Separate`. The paper's construction is the decider-modified table, which is **20% better on `noise_cavity`** (143 → 114 non-manifold edges) and, per the module docs, **worse on `gyroid` at 25³**. Changing it re-baselines every golden hash.

### What is already settled

- **The paper's claim is falsified** (✗19, M-290). §3 says the uniform-grid dual *"is always a manifold."* Over eight fields at three resolutions MC measures **0** non-manifold edges under both face rules; MDC measures **143** / **114**. The premise holds; *"the dual preserves the topology"* does not.
- **P-17 is falsified (M-291).** The residue is *not* the interior ambiguity. `Interior::Joined` is reported by **100%** of ambiguous-face pairs on `noise_cavity` — offenders and control alike — so the any-axis test has no discriminating power. Restricted to the shared face it discriminates and points the **wrong way**.
- **M-292: no two-cell sign configuration forces the defect.** All 4,096 two-cell sign patterns × 16,384 joined-bit assignments: 512 share an ambiguous face, 18 offend under mask 0, 476 under some mask, **0 under every mask**. But a rule reading only the shared face cannot fix it — the decider leaves 25–49.
- **M-294: the minimal fixture exists.** 48 samples, a 2×2×3 column whose middle plane is the face saddle. MC: 0 defects. DC and MDC: 1. Walking the saddle across zero takes the defect 1 → 0 with the triangle count never moving — **pure connectivity**.
- **M-310: MDC earns its keep on real data.** On the bonsai CT volume it takes non-manifold edges **1,776 → 85**, a 95% cut. The synthetic fields never showed this because they never had the problem.

### Recommendation: **leave the default at `Separate`, write the reason down, and close the ticket's part (2).**

**The 20% is on one field, and that field is the outlier.** `noise_cavity` was added specifically because none of the others produces an interior ambiguity (M-208). Re-baselining 216 golden hashes to win 29 edges on the one deliberately pathological field, while losing on `gyroid`, is a bad trade — and the loss is on a *triply periodic* field, which is the shape closest to real cave systems.

**The mechanism is not understood, and a default should not move ahead of one.** P-17 is falsified and M-292 says no shared-face rule can fix it. Changing the default now would be changing a knob whose effect you cannot yet explain — and M-292's own conclusion is that the fix needs **strictly more context than the face**, which neither table provides.

**The golden-hash re-baseline is not free in a way that is easy to underrate.** Regenerating 216 hashes with `ISOMESH_BLESS=1` means reading the diff by eye. Doing that for a change you cannot explain spends the crate's single best regression instrument on a coin flip.

**What to write down instead** — this is the deliverable, and the ticket already accepts it: *the default is `Separate`; the decider table is 20% better on `noise_cavity` and worse on `gyroid`; the residue is a known defect of 85–143 edges concentrated on fields with interior ambiguities; M-292 proves no shared-face rule removes it; the candidate that could is Grosso's three-quadratic tunnel test adapted to dual cells* (novelty-table row 1).

**What would change my mind:** if the tunnel test lands and removes the residue under *both* tables, then the default choice becomes free and should follow whichever is better on `gyroid`. Also, if a consumer is feeding MDC real CT-like data rather than synthetic fields, M-310's 95% argues the decider table's real-data behaviour should be measured before deciding — **and it never has been.** That is a cheap missing measurement and arguably should precede this decision.

---

## Decision 3 — X-005: take the API break, or write the contract down?

### The question

`extract_into` takes `origin: [R; 3]`; every implementation computes `origin + cell_size · local`. A chunk at a non-zero base reaches its far sample plane as `(o + h·base) + h·local`, its neighbour as `o + h·(base + local)`. **Those are equal by algebra and not by IEEE.**

### What is already settled (M-278)

- Canonical reconstruction gives **0** unmatched seam-plane boundary edges at **every** spacing tried.
- What the crate offers today gives 0 only at a **power-of-two** spacing, and **63–348** at `0.1`, `1/12`, `1/14` — plus a hole **1.05–2.08 cells wide** in 2 of 12 rows where an ulp flipped a sign.
- **The crate's own weld hides all of it (✗18).** An unwelded consumer — M-69's collider — gets it in full.
- The working shape is *one path, not two*: replace `[R; 3]` with `(grid origin, integer base)` and compute `o + h·(base + local)`, degenerating to today's behaviour at base zero. **`TransitionCell::sample` already took this route at A-011b** — the precedent is in the tree.
- **Blast radius, counted not estimated:** 7 inherent `extract` methods behind one `forward_extractor!` macro, 39 `origin: [R; 3]` parameters under `crates/`, **294 call sites across 101 files** — 188 in 45 files under `crates/`, 106 in 56 files under `bevy_isomesh/`, a separate workspace with its own lockfile and CI.

### Recommendation: **take the break (a), now.**

**The version number will never be lower.** The crate is at 0.0.9. `X-001` exists to *stabilise* the `Extractor` trait; changing its signature after that is a different conversation entirely. Doing this at 0.0.x costs a changelog line. Doing it at 1.x costs a major version and everyone's afternoon.

**Option (b) writes down a contract the crate cannot enforce.** "Use a power-of-two cell size for a chunked world" is not checkable by any gate in the repo — nothing rejects `cell_size = 0.1`. Contracts a compiler cannot see are found by users, in the field, as a 2-cell hole in a collider. And the crate serves a **CAD** consumer, where `0.1` mm is not an exotic cell size, it is the obvious one.

**The weld argument cuts the wrong way.** ✗18 says the weld hides every hairline — so nothing a *welded* consumer sees changes either way. That is an argument that the break is **low-risk**, not an argument that the bug is unimportant. The population that gets hurt is the unwelded one, which is the collider path, which is where "the player falls through the floor" lives.

**294 call sites is a big number that is mostly mechanical.** 188 are behind one macro plus 39 parameters; 106 are Bevy examples that a compiler error walks you to one at a time. This is an afternoon of tedium, not a design problem — and the design is already validated by A-011b.

**Do not ship both paths.** The ticket says so and it is right. An added `extract_based` alongside `extract` is two behaviours for one question, and every future bug report starts with finding out which one the reporter called.

**What would change my mind:** if a real downstream consumer is already pinned to `isomesh` and would be broken by this, the calculus changes — the break is still right, but it wants coordination and a deprecation cycle rather than a single commit. **You know whether that consumer exists and I do not.** Also, if `bevy_isomesh`'s 106 sites turn out to be mostly in `examples/common/`, the real cost is far lower than the count suggests and this becomes an easy yes.

---

## Decision 4 — R-027 / R-027a: sequencing

### The question

Not really a design decision — a sequencing one. R-027 is stopped, and the ticket that would unstop it is sitting unstarted.

### What is already settled

- **M-314:** computation after an edit is edit-proportional (792 dirty cells, constant across a 64× lattice). The **output buffer** is not — **56–77%** of vertex slots change for an edit touching **0.038%** of cells, growing `O(n²)`. **The cause is a counter.**
- **M-318:** a grid-edge naming takes churn **15,706 → 346** at 129³ — **45×, flat in n**. Three shapes; only the third works. Stable *order* does not help. Index-is-edge-id costs **230× memory**. A persistent edge→slot map is what an incremental engine runs — and it is state across calls.
- **V-45 is a stop, not a preference.** `validate::determinism.rs:268` runs `check_determinism` **three** times, the third into a **reused** buffer, under a doc comment saying why in as many words: *"to catch output that depends on the buffer's prior state."* R-027's only working shape **is** that. It does not cost a migration on top of a working design; **it converts a shipped gate's failure condition into its intent.**

### Recommendation: **take R-027a immediately. It is the only unblocked ticket in the file, it is an `S`, and it may dissolve an `L`.**

R-027a decomposes the churn into three terms — vertices that genuinely moved, slots shifted by a predecessor cell's triangle count, and order alone. Its hypothesis is that the **order-only term dominates**, in which case a **canonical reorder at emission** recovers most of the 45× while needing no persistent state, leaving T-004 intact.

It needs no new field and no new extractor — just instrumenting the existing `edit_trace` bench. **Both outcomes are worth having:** if the geometric and predecessor terms dominate instead, then only the persistent map can help and R-027 is genuinely stopped rather than merely blocked, which is a real result and closes an `L` honestly.

**One thing to add to the ticket, and it is worth adding.** A sweep over 9,425 documents plus a live provider search found **no published scheme** for stable identity of isosurface output elements without carrying mutable state across calls. The persistence literature solves cheap access to *old versions* and is stateful by construction; dynamic connectivity maintains *aggregates*, not buffer slots. **This is unoccupied territory** — the same shape of finding as V-43, and it deserves the same care in the record. If the canonical-reorder shape works, it is novel, and the novelty is worth stating precisely rather than by implication.

**Do not reopen R-027 on a softened version of the persistent map.** V-45's reasoning holds: the only reopening condition is a formulation where the map is **derivable from the inputs** — a pure function of grid and field rather than of call history. Nothing in M-318's three shapes is that. A cache keyed on inputs and *rebuilt* rather than *carried* might be; that is a real question and it is downstream of R-027a's numbers, not upstream.

---

## Summary

| Decision | Recommendation | Confidence | Cost of being wrong |
|---|---|---|---|
| **A-026** scope | **(ii) narrow first cut**, split into a-d, MCTS deferred | High | Low — if greedy cells prove inadequate, A-026d is still there and A-027 has been measured meanwhile |
| **A-025** default | **Leave `Separate`**, write the reason down, close part (2) | Medium-high | Low — reversible; the hashes are only re-baselined if you change it |
| **X-005** API | **Take the break now**, at 0.0.x, one path only | High **unless** a pinned downstream consumer exists — which you know and I do not | Moderate — 294 sites is real tedium, but the design is pre-validated by A-011b |
| **R-027a** | **Start it today.** Only unblocked ticket, `S`, may dissolve an `L` | High | None — both outcomes are results |

**One cheap measurement is missing and would sharpen A-025:** MDC's decider table has never been measured on real CT data, only on synthetic fields. M-310 showed the synthetic fields systematically understate what MDC is for. That is a bench run, not a ticket.
