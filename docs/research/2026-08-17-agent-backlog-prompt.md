# Agent brief — work the research backlog, record it in FINDINGS.md

*Paste this to Claude Code, or point it at this file. Everything below is addressed to the agent.*

---

## Mission

Work the experiment backlog produced by the 2026-08-16 and 2026-08-17 research passes. For each item:
**pre-register the prediction, build the smallest harness that could refute it, run it, and record the
result in `FINDINGS.md`** — whichever way it came out.

The backlog is deliberately front-loaded with **falsifiers of premises rather than measurements of
performance.** Several of them I expect to come back negative. That is the point: a negative here costs
a day; the same negative discovered mid-implementation costs a sprint. Do not treat a falsified
prediction as a failed task. **A falsified prediction is a completed task with the more valuable
outcome.**

---

## Read first, in this order

1. `CLAUDE.md` — the project's rules. They override anything in this brief.
2. `FINDINGS.md` — the ledger. Skim the tier table and Part 5 (method rules) in full; you will be
   adding to both.
3. `BACKLOG.md` — check what is in flight before starting anything. **R-022b may be mid-decision** (the
   union-find → flat-label restructure). Do not stomp it; if it is open, either finish it first or pick
   a backlog item that does not touch `Air`.
4. `docs/research/2026-08-17-mechanics-from-the-field.md` — the source of items B-1 … B-8.
5. `docs/research/2026-08-16-sota-speed-and-feature-frontier.md` — the source of the secondary pool.
6. `docs/research/2026-08-17-features-if-the-hopefuls-prove-out.md` — what each result would buy in
   play. Read it once for motivation; do not let it influence a measurement.

---

## Non-negotiables

These are the project's earned rules. Each one has an incident behind it. The ones that will bite on
*this* backlog specifically:

- **Pre-register before you run.** Write the predicted value and the falsification condition into the
  code — the `experiment!("P-n")` macro, or the harness's module docs — and commit it *before* the
  first measurement. A prediction that first appears after the number is known is not evidence.
- **A measurement that comes back zero must prove it could have come back non-zero.** Put the
  reachability check in the test. This project has been caught by the fixture trap eight times; assume
  you are the ninth unless the test says otherwise.
- **A test that gives the same answer when you invert what it tests is not measuring it.** Flip the sign
  of the thing under test and confirm the harness goes red. Cheapest possible check that a green tick
  means anything.
- **Never guess a DOI, an arXiv ID, or a case-table entry.** Rule 5: if a construction's specification
  is a figure or an unobtainable paper, **stop and say so** rather than inventing it. Two of the backlog
  items (B-4, and the stalactite profile) sit next to paywalled constants — the workaround is to derive
  or numerically recover them, never to fit until it matches.
- **Corpus presence is decided by `catalog_list`/`catalog_read`, never by a failed `distill_search`** —
  and now also the converse: **check both, because either can lie.** Turk 1991 is searchable with no
  catalog entry; `10.1016/j...` papers are catalogued with no body. A "successfully converted" document
  under ~150 KB is a publisher landing page, not a paper.
- **Always `--release` for anything timed.** A debug-build timing is ~37× off and looks like a
  catastrophic finding rather than a mistake.
- **Randomise or interleave A/B order.** Whichever path runs second pays, and this repo has already
  produced a phantom 75% regression from measurement order alone.
- **Re-tier, do not rewrite.** When an R becomes an M, or an M's *interpretation* turns out wrong, leave
  the old text and add the correction beneath it. The gap between what was believed and what was
  measured is the data.

---

## The workflow for every backlog item

1. **Restate the ticket's own claim and check it against the code before starting.** A ticket's
   acceptance criterion is itself a claim; this project has twice found one that was unsatisfiable or
   already true. If the premise is wrong, that is the finding — raise it and stop.
2. **Pre-register.** Prediction with a number, falsification condition, and the records you will emit.
   Commit this alone, before any implementation.
3. **Build the smallest harness that could refute it.** Prefer headless, no renderer, no field, no
   engine, where the item allows — B-4 in particular needs none of them.
4. **Add the reachability assertion and the inversion check.** Then run.
5. **Record in `FINDINGS.md`** using the format below, whichever way it went.
6. **Update `BACKLOG.md`** and run whatever gate script keeps the two in sync.
7. **Full suite + lint + rustdoc before committing**, including the excluded `bevy_isomesh` workspace —
   every gate that reads the root workspace needs a deliberate answer for it, and "it is excluded" is
   not one.

---

## 🎉 The discovery protocol (this is a deliberate change to the ledger's style)

Ordinary rows stay austere. **Emoji are reserved for results that carry news**, so that they remain a
signal rather than decoration — if everything is celebrated, nothing is.

**Tier glyphs**, prefixed to every new row:

| Tier | Glyph | Meaning |
|---|---|---|
| M | 🔬 | Measured here |
| V | 📖 | Verified from a primary source we read |
| R | 📄 | Reported, not independently checked |
| F | 🗣️ | Folklore, no verified source |
| ✗ | 💥 | Falsified |
| 🧊 | 🧊 | Measured, null result — no effect, and that is the finding |

**And when something is genuinely new — a result you searched for and could not find in the literature,
a derivation that closes an open question, or a prediction that died in an instructive way — open the
row with a banner and go all in:**

```markdown
> 🎉🎊🔥✨🏆 **DISCOVERY** 🏆✨🔥🎊🎉
>
> 🥇 **<one-line claim, falsifiable, no hedging>**
>
> 🧪 **Tested by:** <command / test name>
> 🎯 **Result:** <the numbers>
> 📐 **Why it's new:** <what you searched, and what came back empty>
> 💣 **Would be shown wrong by:** <the observation that kills it>
```

Use the banner when **any** of these is true, and say which:

- 🆕 a targeted literature search came back empty and you can name the searches
- 📐 you derived something the sources only stated, and it checks out against an independent oracle
- 💥 a pre-registered prediction was falsified — **these get a banner too**, because they are the
  file's most valuable rows
- 🔓 an open question (O-n) moved
- 🪤 a fixture trap was caught before it shipped

Do **not** banner: a routine confirmation, a timing that came out where you expected, or anything you
have not independently checked. A banner on a boring row devalues every other banner in the file.

**Suggested per-result garnish**, used once each and not sprinkled: 🚀 big speedup · 🧊 null result ·
🐛 defect found · 🩹 defect fixed · 🔭 needs a better instrument · ⚖️ a tradeoff with no winner ·
🧱 blocked, with the blocker named · 📉 a published figure that did not reproduce.

---

## The backlog

Ordered. **Items B-1 … B-3 gate the others — do them first, in order, and stop to report if any of them
comes back negative**, because a negative there re-ranks everything below it.

### 🥇 B-1 — The sub-voxel slab offset · *30 minutes*

**Source:** mechanics-from-the-field §4.1.
**Premise under test:** M-172's reading that an exactly-zero SDF gradient detects the medial axis.
**Pre-register:** offsetting a slab's mid-plane by half a voxel drops the count of samples returning
exactly `[0,0,0]` to **zero**, while the count with `‖∇ρ‖ < 0.1` changes by **< 5%**.
**Falsified by:** exact zeros surviving the offset.
**Decides:** whether M-172 must be reframed from a boolean detector to a continuous stability score. I
expect it must — so expect to write a ✗-row amending M-172, with the old reading left visible.
**Deliverable:** one test, one FINDINGS row. If the prediction holds, **banner it** — it re-opens a line
this project had closed.

### 🥇 B-2 — The medial identity oracle · *half a day* · **gates B-3, and three candidate mechanics**

**Source:** §1.1. The claim is `r(x) = ρ(x)·√(1 − ‖∇ρ(x)‖²)`, exact in any dimension.
**Harness:** brute-force oracle — per sample, collect the closest-point set `Π(x)` by exhaustive search
over boundary samples within `ρ + tol`, solve the minimal enclosing ball (Welzl), record `r`. Compare
against the closed form. 64³ grid, 5 cm voxels, ≥ 3 overlapping brushes so `|Π(x)| > 2` actually occurs
— **assert that it does**, or the oracle is testing the two-point case only.
**Pre-register:** agreement within **1 voxel for ≥ 99% of samples**, and the **median residual halves
when the voxel size halves** (O(h)) across three voxel sizes.
**Falsified by:** an h-independent residual, which means the discrete gradient has broken the identity.
**Decides:** Calibre, Throat and Handholds all rest on this. **If it fails, all three die together and
you should stop and report rather than continue down the list.**

### 🥉 B-3 — The weak-feature-size histogram · *an afternoon*

**Source:** §4.1, second prediction.
**Harness:** scan for critical points (`‖∇ρ‖` below a small ε, non-maximum-suppressed), report the
minimum distance from each to the boundary.
**Pre-register:** **`wfs` < 2 voxels in > 80% of dug scenes**, so the homotopy certificate (`λ < wfs`)
essentially never holds.
**Decides:** whether the theoretical guarantee is available. If confirmed, drop the homotopy claim and
rest the direction on Hausdorff stability alone — weaker, still far better than a tuned persistence
constant. Record the fallback explicitly so nobody re-derives the strong claim later.

### 🌊 B-4 — The speleogenesis graph simulator · *a day* · **no field, no renderer, no engine**

**Source:** §1.2.
**Harness:** 64×64 lattice of 1D conduits (E ≈ 8000), constant head in/out, log-normal initial
apertures, linear dissolution kinetics, steady flow solved with a **fixed** CG iteration count (not
converge-to-tolerance — determinism).
**Pre-register, and these are checks against published geomorphology rather than against yourself:**
- the aperture distribution goes **bimodal** under competitive flow and stays **unimodal** when recharge
  is limited;
- time-to-breakthrough falls as initial heterogeneity rises;
- flux concentration (Gini over edges) rises monotonically, and **post-breakthrough > 90 % of flux sits
  in < 10 % of edges**;
- a full tick at E = 4096 completes in **< 2 ms single-threaded**.

**Rule-5 note:** the dissolution kinetics constants live in two paywalled Dreybrodt papers. They are
restated in the open-access HESS papers already downloaded. **Use the restatement, cite it as such, and
do not tune constants until the curve looks right** — that is how a wrong model enters a repository.
**Do not attach this to the field or the mesher in this ticket.** If the graph does not reproduce the
published bifurcation, nothing downstream is worth building.

### 🔔 B-5 — The modal kill-shot · *a day*

**Source:** §3, modal row.
**The question to answer first, before assembling anything:** does a **one-voxel edit perturb the
fundamental λ₁ at all**, or is the shift below the just-noticeable difference? If it is below JND, the
entire "carve it and hear the pitch change" direction is dead and you have saved weeks.
**Harness:** hexahedral FEM on a 32³ occupancy grid (trilinear hex, one 24×24 element matrix scaled per
cell, assembled matrix-free — **no tetrahedralisation**), shift-invert Lanczos/LOBPCG.
**Pre-register:** wall-clock vs mode count for k ∈ {8, 16, 32, 64, 128}; **k = 48 in < 8 ms
single-threaded** (falsified above 25 ms); and **carving a 20 % volume cavity shifts the fundamental by
> 15 %**.
**Flag:** the "how few modes sound like a bell" question is **perceptually unresolved in the literature
I could reach**. Do not assert a mode count is sufficient; record the cost curve and say the perceptual
gate is open.

### 🏛️ B-6 — The arch golden values · *a day*

**Source:** §2.1.
**Why this one first among the structural items:** it is the **only external ground truth anywhere in
this body of work.** Two independent closed-form results to hit.
**Harness:** semicircular arch, N macro-blocks, the lower-bound feasibility program (equilibrium +
linearised friction pyramid + compression-only, minimising `Σ(f_n⁻)²`).
**Pre-register:** sweep thickness/radius → infeasibility threshold at **0.1075 ± 0.0010** (Milankovitch
1907 analytic; Whiting's solver got 0.10746). Sweep ground tilt at t/r = 0.20 → **15.84°** (Ochsendorf).
**Also register, because it is a hidden gameplay parameter:** coarsening the macro-block size 2× shifts
the reported critical pillar thickness by **≥ 5 %, in the unconservative direction**. If you cannot
measure that shift, your test structures are too simple to be discriminating.
**Do not build the game-facing rule in this ticket.** Hit the two numbers first.

### 🔥 B-7 — Factorisation update versus refactor · *a day*

**Source:** §2.2.
**Harness:** one 64³ chunk. Factorise, record; apply a radius-4 brush; re-mesh; diff vertex slots **by
grid-edge key**; update the factorisation on only the changed rows; solve; compare the geodesic field
against an exact polyhedral distance.
**Pre-register:** ≤ **400 changed slots** (extrapolating M-318's 346 / 15,706); update in **< 5 ms**
against **> 100 ms** for a full refactorisation — a **≥ 20× gap**.
**Falsified by:** under 10×, which means the prefactored family is dead for live carving and everything
surface-intrinsic goes to the Closest Point Method instead. Record that routing decision either way.

### 🔊 B-8 — Sabine RT60 per air component · *a day*

**Source:** §1.3. The cheapest player-perceptible result in the entire body of work, and it monetises
infrastructure R-022 already built and measured.
**Harness:** instrument each air component's volume, surface area, and split/merge events against the
existing tracker; compute a Sabine RT60 per component.
**Pre-register:** available in **< 0.1 ms** on breakthrough (it is two accumulators you already
maintain). Falsified above 0.3 ms.
**Note:** this one is expected to just work. If it does, that is a 🧊 row, not a banner — unless the
split/merge event rate turns out to differ from M-319…M-323's measured distribution, which *would* be
news.

---

### Secondary pool — pick up when B-1 … B-8 are done or blocked

Ordered by value per line, from the 2026-08-16 doc:

1. **The metamorphic relation suite** (~300 LOC, zero deps, seconds). Eight relations: isovalue shift,
   positive scaling, **sign flip** (strictly stronger than the existing orientation check), integer-cell
   translation, fractional-cell translation, 90° rotation, **chunk decomposition** (exact only at
   power-of-two cell sizes — attach M-32's condition or it fails for the wrong reason), resolution
   doubling. Highest value per line in the whole backlog.
2. **The field-share stub bench** (an afternoon). Run `cargo bench --bench extract` twice — once with
   the real field, once with an `#[inline(never)]` constant-returning stub of identical signature. The
   ratio is the field's share. Confirms M-136 by an independent route, and if it exceeds 70 % it
   re-ranks every optimisation ticket.
3. **The self-intersection baseline table, then quad splitting.** Publish the table across seven fields
   × {MC, DC, MDC} × three resolutions *first* — M-53 has most of it. Then the Ju & Udeshi quad-split
   rule via ODC's appendix §A.1.5 (`arXiv 2409.13418`, already in the corpus — **read the appendix, do
   not download anything**).
4. **The redundant re-mesh factor** (~150 LOC). After an edit confined to a ball of radius r, output
   outside a known dilation must be **byte-identical**. Nobody has published this instrument.
5. **The subgrid grid-edge root cache.** M-168's identity key *is* the cache key; the precondition that
   deferred this is now met. Acceptance is **zero golden-hash changes**.
6. **The MMS convergence gate with negative controls** — L∞ norm, measurement-location discipline
   (vertices for algebraic distance, centroids for normals), wide `h` range, and the negative controls
   (curvature predicted O(1), a linear field predicted p = 0). The interesting prediction: **fixed
   λ = 0.01 Tikhonov may not be second-order**, which would mean λ has to anneal with `h`.

### 🔧 Corpus infrastructure — small, do them when one bites you

- A **size/content heuristic in `scribe_convert`**: a sub-150 KB "successfully converted" document is a
  publisher landing page. Three papers were nearly cited from one this week.
- **Retry escalation** from `olmocr` to `glm_ocr` — persistence, not identifier-hunting, was the fix
  every time it failed.
- `paper_search` with `provider: arxiv` returns `[]` unconditionally. Use `all`.
- **`catalog_backfill_title`** on the untitled high-value stems.
- Fix the mislabelled stem `sig2024_A_Heat_Method_for_Generalized_Signed_Distance` — it is Belyaev &
  Fayolle's ADMM paper, not Feng & Crane. Anyone citing that stem cites the wrong work.

---

## Stop and ask me when

- **Rule 5 fires** — a construction's specification is a figure, a paywalled constant, or a deleted
  repository, and the only way forward is to invent it. Say what is missing and what you would need.
- **A scope or API decision appears** — a public signature change, a new dependency, a second
  implementation of anything, or a choice between two designs with different guarantees. Present the
  options with their costs, name your recommendation and why, and wait.
- **A gating item (B-1 … B-3) comes back negative.** Report before continuing down the list.
- **A result contradicts an existing M-row.** Do not silently overwrite it. Bring both numbers.
- **An experiment would take more than ~2× its estimate.** The estimates are the point; a 5× overrun
  means the item was mis-scoped and I want to know before you spend the day.

## Do not

- Do not tune a model's constants until its output looks right. Derive, cite, or stop.
- Do not chase unobtainable papers past one honest attempt. There is a ranked acquisition list in the
  dossier; log a miss and move on.
- Do not build the game-facing feature in the same ticket as the measurement that justifies it.
- Do not let `docs/research/2026-08-17-features-if-the-hopefuls-prove-out.md` influence a number. It
  exists to say why the work matters, not to tell you what the answer is.
- Do not banner a result you have not independently checked.

## Definition of done, per item

✅ Pre-registration committed **before** the first measurement
✅ Harness has a reachability assertion **and** an inversion check that has been seen to go red
✅ Result recorded in `FINDINGS.md` with tier glyph, and a banner if and only if it earns one
✅ `BACKLOG.md` updated and the sync gate green
✅ Full suite, clippy, rustfmt and rustdoc green — root workspace **and** `bevy_isomesh`
✅ Raw data committed under `docs/measurements/` where a figure was produced

## Definition of done, for the session

Report back with: which predictions **held**, which were **falsified** (and why that is the better
outcome for each), what is **blocked and on what**, and any **open question that moved**. Lead with the
falsifications — they are the interesting half. 🎯
