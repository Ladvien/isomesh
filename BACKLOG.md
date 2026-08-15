# isomesh — BACKLOG

**Updated:** 2026-08-15
**Companions:** `CLAUDE.md` (rules), `FINDINGS.md` (what we know and how well),
`BACKLOG_ARCHIVE.md` (completed tickets + why they changed),
`docs/2026-08-11-implementation-brief.md` (the how),
`docs/2026-08-11-bevy-examples-catalog.md` (example detail), `docs/research/` (the why).

**122 tickets archived, 7 open.** Completed rows move to `BACKLOG_ARCHIVE.md` with their amendments
attached — read that before re-litigating a decision this project already made.

---

## How to work this backlog

1. Take the **topmost unblocked, unchecked ticket**. Don't cherry-pick interesting ones — the order
   encodes dependencies and the test harness exists so later work is cheap.
2. One ticket = one commit (or a short stack). Commit message starts with the ticket ID.
3. **Check the box in this file as part of that same commit.** This file is the state.
4. If a ticket can't be finished, leave it unchecked, add a `> BLOCKED:` line under it saying exactly
   what's in the way, and move to the next unblocked ticket. Do not half-finish and check the box.
5. If a ticket turns out to be wrong or to need splitting, edit it and say so in the commit.
6. **On completion, move the row to `BACKLOG_ARCHIVE.md`** with an indented annotation recording any
   amendment, deviation or falsified premise. The annotation is the point; the checkmark is not.
7. **New tickets** slot in by dependency, not by number. A ticket split after the fact takes a letter
   suffix (`T-005a`/`T-005b`, `A-002`/`A-002b`); a genuinely new one takes the next free number in its
   series even if that puts it out of numeric order (`A-015`).

### Definition of done — applies to every ticket

- Code compiles with no warnings. `cargo clippy -- -D warnings` clean.
- **`RUSTDOCFLAGS=-D warnings cargo doc --workspace --no-deps` clean.** A third of the lint job, and
  the third that clippy and fmt do not cover — it is what catches a doc link pointing at nothing, or
  at a private item. Added at A-002, which shipped one and found out from CI.
- `cargo test -p isomesh` green. `cargo fmt` applied.
- **No Bevy dependency and no Bevy in code under `crates/`** — comment-stripped manifests clean,
  non-comment `.rs` clean, resolved-graph bevy count 0. Non-negotiable — see `CLAUDE.md` rule 2. Prose
  explaining the wgpu-follows-Bevy pin is not a breach; the CI gate checks the three forms above.
- Any algorithm ticket also passes the T-001 validity suite. No exceptions, no "I'll add tests later."
- Any perf claim has a committed benchmark that produced it.
- Public items have doc comments. Anything with a sign convention, a coordinate order, or a winding
  order says so **in the doc comment**, not in a code comment.
- **`FINDINGS.md` updated in the same commit** if the ticket measured something, contradicted
  something written down, or earned a method rule. A measurement that only exists in a commit message
  is not retrievable six weeks later.

**Size key:** `S` ≈ one sitting · `M` ≈ a day · `L` ≈ multi-day, consider splitting.

---

## Phase 0 — Foundation and the test harness ✅ complete

All eleven tickets archived (I-001..I-004, T-001..T-008). The bet paid: every algorithm since has
been cheap to validate because the harness predated it.

---

## Phase 1 — The usual suspects

Each algorithm ticket is done when: T-001 reports **no unexplained violations** on all **eight** test
fields at three resolutions; T-004 determinism passes; T-005 covers it; and a benchmark exists.

> **The eighth arrived at A-002e, and it is the one that can fail things (M-208).** `noise_cavity`
> exists because none of the original seven produces a cell with an interior ambiguity — 0 of 68,385
> surface cells — so five pre-registered claims were properties of the fixtures rather than of the
> code. A-017, A-018 and A-019 own what it found. Expect a new algorithm to need a pinned census on
> this field where the other seven give zero.

> **Amended 2026-08-12.** This originally said "zero violations on all seven fields." M-4 falsified
> that as a universal gate: Surface Nets is *legitimately* non-manifold where one cell carries two
> sheets (48 edges on capped gyroid, 15 on fbm_terrain), and those counts are pinned as **non-zero
> assertions** precisely so they can't drift silently. A known defect with a pinned number and a
> ticket that owns it satisfies this gate. An unexplained one does not.
>
> **A-010 has since landed and closed that owner.** Surface Nets' and Dual Contouring's counts stay
> pinned as non-zero — they are properties of one-vertex-per-cell, not bugs — and
> `manifold_dual_contouring` is the entry that takes the zero. Its own residue is M-59's parallel-edge
> collapse, pinned at one edge on the ✗15 fixture and zero everywhere else, owned by O-16.

> **Re-scoped 2026-08-14 (A-002d).** A-002b was one `L` blocked on size. It is now a series, because
> two of the three things that made it large turned out not to be there. **(1)** Its route was Custodio
> 2019, whose non-manifold remedy is a *grid* preprocessing pass. Grosso 2016 (`10.1111/cgf.12975`)
> and Grosso 2017 (`10.1145/3095140.3095179`) — recorded as `PAYWALL` in `meshing-library-target.md`
> and in fact indexed in home-still since 2026-08-10 (V-29) — reach the same manifold result with
> interior vertices and **no grid pass**, keeping unambiguous cells on the existing table path.
> **(2)** *"Cell-interior vertices, which the grid-edge-keyed vertex cache has no slot for"* is true
> and is not a blocker: **A-015 already built that mechanism** — `table::CENTROID_BASE`, created per
> cell at `mod.rs:208-236`, uncached by design, already budgeted in the `u32` bound at `mod.rs:140`.
> What is left genuinely is the meshing, in the five tickets below.

| | ID | Ticket | Size | Blocked by |
|---|---|---|---|---|
| ☐ | **A-002g** | **Disk triangulation — Grosso §5.3.** The bulk of the win and none of the drama: on the authors' own 512²×641 CT skull, **107,012** ambiguous disk contours against 2,057 tunnels and 7 twelve-vertex contours. Contours of six or fewer vertices are a flat polygon and fan trivially. Ambiguous contours of 6, 7, 8 or 9 vertices get **one** interior vertex and fan from it — the body saddle where the two single lines cross, or the midpoint of the two saddles where there are three such lines. Selecting it needs the per-face-pair count of *single* in-range solutions, which is why A-002d solves the linear and tangent cases honestly rather than inheriting the reference's textbook formula (M-207). Extend the entry code space above `CENTROID_BASE` with an interior-vertex code whose **position comes from the classifier** rather than from averaging edge vertices, and re-derive the per-cell bound at `mod.rs:140-141` rather than reusing `MAX_CENTROIDS`. Keep `edge_position`'s split from A-015 intact — that ticket measured **~18% slower** plain Marching Cubes when the position computation moved above the cache check. | M | A-002f |
| ☐ | **A-002h** | **Tunnel and twelve-vertex contour — Grosso §5.1 and §5.2.** Where the interior ambiguity is actually meshed and where manifoldness is bought. Inner hexagon from Proposition 2; three inner vertices as the midpoints of consecutive hexagon vertices; contour vertices classified by nearest hexagon vertex; the two meshes merged by collapsing consecutive hexagon vertices. **The rule-5 hazard this was flagged for is already discharged (V-31).** The paper's own listing of the hexagon is corrupt in the copy held here — its first branch assigns `p₁` and `p₂` the same triple — and the authors' implementation, recovered from Software Heritage after GitHub deleted it, gives the order outright. A-002d ships it as `BodySaddles::inner_hexagon`, verified as a closed axis-parallel ring over 1,146 hexagons. | M | A-002g |
| ☐ | **A-002i** | **The singular case — Grosso 2017 §4.2 and its Algorithm 1.** A saddle sitting exactly *on* a face, where the standard asymptotic decider splits into two branches what is one surface. This is A-002b's own constraint (c): the 2 of M-165's 15,625 opposed configurations where a root of `F` rests on Δ's pole to within last ulps get a **defined** answer here rather than an inherited one. The reference implementation's shape is a per-face singular flag that then snaps the affected quadratic root to 0 or 1, plus an edge-coordinate comparison to choose the face pairing when the decider itself has no answer. Rare but real, and measured rather than assumed: Grosso 2017 Table 1 counts **8, 58 and 20** singular faces across three 512²×~700 CT volumes (tier V). | S | A-002h |
| ☐ | **A-002b** | **Marching Cubes 33 interior ambiguity — the trilinear body saddle.** Deliberately deferred at A-002, on evidence, not forgotten. Three reasons. (1) `catalog-v2.md:107` is explicit: *"Skip the interior test; spend the budget on chunk seams"* — a game needs topological *consistency*, which A-001 already has, over *correctness*. (2) ~~**There is no correct published table to transcribe.**~~ **Answered on 2026-08-13 by reading the sources (V-24, V-25, V-26), and two thirds of it was wrong.** Custodio et al. 2013 (`10.1016/j.cag.2013.04.004` §5.1) do prove Chernyaev's interior test tracks a quadratic where the true saddle trajectory is hyperbolic with an asymptote — and **that correction is now implemented and tested (A-002c)**, so it is no longer this ticket's. The claim that *"Lewiner's reference implementation omits disambiguation for cases 10 and 12 entirely"* is true of the code and false of the literature: §5.4 calls it *"a missing step in the implementation"* and states the rule for 10.1.1, 10.1.2 and 10.2 outright. And the missing table is not missing but unnecessary — **`10.1186/s13173-019-0086-6` (Custodio, Pesco & Silva 2019) builds the whole MC33 triangulation with no lookup table**, from groups of vertices and edges and the boundaries of their convex hulls, with case 13.3, 13.5.1 and 13.5.2 given as named combinations. Rule 5 is satisfied by following a construction rather than by inventing a table. (3) The v1 catalog prices it: the decider is *"~free"*, the guaranteed version is **730 subcases in the LUT**. Also needs cell-interior vertices for tunnels, which the grid-edge-keyed vertex cache has no slot for. **Acceptance:** a cell where the body saddle says "tunnel", meshed as a tunnel, with the sign tracked by Custodio's correction rather than Chernyaev's `F(t)` — the second half of which A-002c has already done, so what remains is the meshing. **Three things it must settle that were not visible before.** (a) **Crack-freeness is not inherited (M-166).** `ambiguity`'s two cells cannot disagree because its decider compares two *products* and IEEE multiplication is commutative; the interior test's denominator is `((A + C) − B) − D`, whose subtraction order a rotation permutes, and float addition is not associative. Two cells reading one shared face could therefore disagree about a tunnel, and that has to be closed before this is wired in. (b) **The non-manifold remedy is a grid change, not a table change** — both papers say to split the two cells at the shared ambiguous face's critical point, which is a preprocessing pass over the grid rather than anything in the extractor. (c) **Added by the 2026-08-14 review:** 2 of M-165's 15,625 opposed-family configurations rest a root of `F` on Δ's pole to within last ulps, so the decider's answer there is rounding noise rather than a decision; the meshing pass owes those configurations a defined answer, not an inherited one. **What remains under this ID is the wire-up and the acceptance**, the five tickets above having taken the rest: `set_interior_ambiguity(InteriorAmbiguity::{Ignore, Trilinear})` and the one branch in `extract`; the `marching_cubes+trilinear` golden row (192 → **216**, and the hard-coded product with it); a `check_trilinear` beside `check_mc33` with its two proptest properties; a `check_determinism` test; the `for_each_reference_field!` sweep; `bench_trilinear` beside `bench_mc33`; and a `bevy_isomesh` example, which M-40 says must use `gyroid`, `fbm_terrain` or A-002e's field or it will render two identical meshes. Leave `resolution_sweep.rs` alone — a row there rewrites the CSV ✗14, M-19–M-22 and O-11 quote. **Acceptance, three parts:** (1) a cell the classifier calls a tunnel, meshed as a tunnel; (2) on A-002e's field, `is_closed()` holds and the non-manifold census is **zero**, which the face decider alone cannot reach on a tunnel cell; (3) **Grosso's own correctness test (§7)** — mesh at `n`, `2n` and `4n` by trilinear refinement and assert the component count, boundary count and χ per component agree. This repo has no equivalent invariant yet and it is the one that grades topology rather than geometry. | M | A-002i |
> **No longer blocked (2026-08-14, A-002d).** It is now sequenced behind A-002e–A-002i above rather than stopped, and the two constraints that made it `L` are discharged: (b)'s grid-subdivision pass belongs to Custodio's route and not to Grosso's (V-29), and the vertex-cache objection was answered by A-015 before this ticket was written. Constraint (a) was already closed by M-204, and (c) is A-002i's. The note below is kept as it stood.
> BLOCKED-was: **on scope, not on a source — the rule-5 stop is lifted (V-24, V-25, V-26).** The previous note said *"there is no correct published table to transcribe"*; reading the papers showed the decider is stated in full, the case-10 rule is stated in full, and the 2019 follow-up builds the triangulation **without any table at all**. What is left is genuinely large rather than genuinely blocked: MC33's tunnel cases need vertices in the cell **interior**, and this crate's vertex cache is keyed on grid edges with no slot for one; the non-manifold case needs a grid-subdivision pass; and ~~M-166 has to be closed first, because two cells can currently disagree about a tunnel across a shared face~~ — **M-166 is closed (M-204).** The denominator was `((A + C) − B) − D`, a subtraction order a rotation permutes; grouping each diagonal as `(A + C) − (B + D)` makes every symmetry of the face exact, and two cells reading one face now agree bit-for-bit. Note the test already asserting that property **passed with the defect present** — its fixture's corners are all of similar magnitude — so the fixture that exposes it had to be searched for (`(1, 1, 1, 10⁻⁸)`), and the new test is mutation-verified in both directions. **So this ticket is now blocked on size alone**, which is where its own note already put it. `catalog-v2.md:107` still prices all of that as *"skip the interior test; spend the budget on chunk seams"* for a game, so this stays where it is in the queue — but it is now a size decision rather than a correctness one.


> **A-017, A-018 and A-019 sit below the series deliberately (2026-08-15).** They were
> discovered *by* A-002e rather than needed *for* it: nothing in the A-002 series depends on
> them, and rule 7 orders new tickets by dependency rather than by number. They are defects in
> three other components, each with a pinned census already failing in both directions, so
> nothing is lost by finishing the series that found them first.

| | ID | Ticket | Size | Blocked by |
|---|---|---|---|---|
| ☐ | **A-017** | **Manifold Dual Contouring is not manifold on a cell with an interior ambiguity.** A-010's headline property — *"`manifold_dual_contouring` is the entry that takes the zero"* — held on seven reference fields because **not one of them can produce the configuration** (M-208). On A-002e's `noise_cavity` it returns **30 non-manifold edges / 60 vertices at 17³** and **64 / 128 at 33³**, and its Euler characteristic parts company with Marching Cubes' (`−30` against `−96` at 17³), which falsifies **P-5** — pre-registered before A-010 ran. Both are pinned as whole censuses (`MDC_NON_MANIFOLD_CENSUS`, `MDC_CHI_CENSUS`) so they fail if they spread *and* if they vanish. **The exact factor of two in every row is the lead**: each offending edge carries exactly its own two endpoints, so this is one edge-level mechanism and not two. **First question, before any fix:** is this a defect in this crate's construction, or the published guarantee not covering a cell whose interior the trilinear interpolant joins? Schaefer, Ju & Warren's argument is about one vertex per cycle per cell; a tunnel is where "one cycle" stops describing the cell. If it is the latter, the honest outcome is a **documented limit** and a pinned census, not a repair — and A-010's claim needs the qualifier in its archive row. | M | A-002e |
| ☐ | **A-018** | **The positional weld can create a non-manifold edge, and one vertex pair proves it.** `subgrid/extract/tests.rs` asserted the opposite outright — *"sharing by identity is a subset of sharing by position, so every counter must be at least as good after welding as before"* — which is true of the sharing and false of the counters: merging two vertices that are geometrically coincident and topologically **distinct** fuses two sheets. On `noise_cavity` exactly **one** pair merges, 5,567 → 5,566, and it adds **2 non-manifold edges and 3 non-manifold vertices** (M-212). **One pair explains three failures**, which is why it is one ticket: the over-merge; P-7 limb (a)'s weld plateau reading `[5567, 5567, 5566]` across the four tolerances T-009 merged; and A-014h's completeness claim, where the exact identity rule removes **0** of the weld's 1 — the last being confirmation rather than a third problem, since the exact rule leaves the pair alone *correctly*. **The margin is inverted, not small:** the pair merges at exactly the `1e-4` policy and stays separate at every finer tolerance, so the usual ratio comes out below one. **Acceptance:** either the weld gains an identity-aware guard that refuses to merge crossings with distinct identities — the data is already there, `exact == raw` still holds on every field — or the claim is rewritten and the plateau restated as conditional. Do not simply loosen the tolerance. | M | A-002e |
| ☐ | **A-019** | **Orientation can raise the flipped-edge count; `after <= before` was assumed and never tested.** Propagation was expected to only ever agree edges it reaches. On `noise_cavity` it *lowers* the residue at 17³ (**1,629 → 1,015**) and *raises* it at 25³ (**1,580 → 2,422**) and 33³ (**1,477 → 3,341**) (M-213). The mechanism is the four-face edge M-187 already names as the residue's cause, now reached often enough to change what propagation *does* rather than only what it cannot fix: with 318 non-manifold edges the flood fill crosses one, commits to one side's winding, and carries a consistent orientation across a patch that is consistent with the wrong neighbour. **M-187's law itself survives** — zero non-manifold edges still implies zero residue, on every field — so what is owed here is the monotonicity claim, not the law. **Acceptance:** either propagation stops at non-manifold edges rather than crossing them, and the census shows the residue never growing; or the growth is explained and the census stands as a documented property with M-187's law restated to exclude it. | M | A-002e |

---

## Phase 2 — Game-shaped infrastructure ✅ complete

All seven tickets archived (G-001..G-007): chunk coordinates, dirty-set re-meshing, brush operations,
field-derived LOD, collider export, the frame budget, and streaming. Still zero Bevy — every one of
them lives in the core crate, so the CAD side gets them too.

The bet that paid here was **G-006's**: the ticket asked for `mesh_within_budget(ms)` and a `no_std`
crate cannot read a clock, so the budget became a predicate the caller owns. That made streaming a
composition rather than a rewrite — `ChunkStream` decides residency, `DirtySet` orders the work
nearest-first, and the caller's own clock decides when to stop.

---

## Phase 3 — `bevy_isomesh`

| | ID | Ticket | Size | Blocked by |
|---|---|---|---|---|

---

## Phase 4 — Examples

Two groups. The algorithm demos are quick and prove correctness visually. **The game-shaped ones are
the point** — they're how someone decides whether this crate is usable.

### 4a — Algorithm demos

| | ID | Example | Blocked by |
|---|---|---|---|

### 4b — Game-shaped

These use the algorithms the way a game does: chunked, edited, budgeted, collided against.

| | ID | Example | What it has to prove | Blocked by |
|---|---|---|---|---|

---

## Phase 5 — Measurement

| | ID | Ticket | Size | Blocked by |
|---|---|---|---|---|

---

## Phase 6 — GPU (do not start before Phase 5)

The speed analysis is explicit that stage placement dominates the extraction algorithm by roughly an
order of magnitude. Which means GPU work is worth doing — and worth doing *after* you know your own
numbers, or you won't be able to tell what the port bought you.

**On the mesh-shader route (E-303), which has no ticket row yet.** Kept here because it was written against GPU-013 and outlived it — GPU-013 was archived on 2026-08-14 and this was never about that row. An earlier version said the route was blocked; it is not — see M-147. **The route needs no `unsafe` in this repository at all.** `isomesh-gpu` never opens a device (its API takes `&wgpu::Device`, GPU-001's rule), and **Bevy writes the experimental token itself**: `experimental_features: unsafe { wgpu::ExperimentalFeatures::enabled() }` at `bevy_render-0.19.0/src/renderer/mod.rs:335`. `WgpuSettings`' default priority is `Functionality`, which requests every feature the adapter advertises, so **Bevy's device already reports `mesh_shader=true multiview=true points=true`** on this machine, measured. E-303 is a Bevy example and gets a mesh-shader-capable device for free; `WgpuSettings.features` is there to force it explicitly if the default ever changes.

**The probe is load-bearing, not belt-and-braces** — an earlier version of this note implied otherwise and was wrong in the opposite direction from the blocked claim it replaced. The free device is *one of three branches*: `WgpuSettings::default()` consults `settings_priority_from_env()` first, so **`WGPU_SETTINGS_PRIO` overrides it**; under any priority other than `Functionality`, `features` starts at `wgpu::Features::empty()`; and `adapter.features()` is machine-dependent. It is also **contingent upstream**: Bevy's line carries `// SAFETY: TODO, see bevyengine/bevy#22082`, an admission that a justification is owed, so if that issue lands as opt-in the default path loses mesh shaders. Track it.

**On "graceful fallback", which needs one distinction rather than a ruling.** A demo that *detects* capability and, finding none, says so plainly and shows the compute path instead is a demo reporting a capability — one path, chosen by a measurement, with the choice visible. What the one-path rule forbids is the *library* silently substituting compute for mesh shaders so a caller cannot tell which ran. The first is what this ticket should build.

**Two things still shape the work.** WGSL mesh shaders are **Vulkan-only** — wgpu's own source says *"naga is only supported on vulkan; on other platforms you will have to use passthrough shaders"* (V-23) — so on Metal a caller supplies pre-compiled MSL and the composed-WGSL pipeline does not apply, making mesh shaders a fork in the shader path rather than a flag on it. And the ticket's own wording needs revisiting: *"graceful fallback"* is a second execution path for one feature, which the one-path rule forbids; the shape that survives is a **capability check that refuses loudly**, as GPU-007's probe already does. The remaining `unsafe`-shaped gap is only `isomesh-gpu::headless::Gpu` opening its *own* mesh-shader device, which is a test convenience rather than this ticket.


---

## Phase 7 — Documentation & packaging

Added 2026-08-14 from a literature pull, per the research-first rule. The sources that shaped these
tickets: *Effective Rust* Item 27 (crates.io is for people **choosing** a crate, docs.rs for people
**using** one — two pages, two jobs), Carroll's minimalist-instruction research (users act, they
don't read: anchor the first success in a real task, and treat error recovery as content rather
than warnings), and Prana et al. 2019 (`10.1007/s10664-018-9660-3`, in home-still), whose finding
is that the "Why" is the content category most READMEs lack. This repo's READMEs have the opposite
problem: the Why is superb and the on-ramp is missing. The falsification-essay voice stays; these
tickets put a doctested front door on it and make every claim current. House rules for the phase:
no hard line breaks in new prose, no performance number without naming its machine and CSV,
absolute URLs in anything crates.io or docs.rs renders, the Vibe Coded label stays on every README,
and every README code fence must be compiled by something.

| | ID | Ticket | Size | Blocked by |
|---|---|---|---|---|

---

## Deliberately not in scope yet

Recorded so they don't get picked up early, and so it's clear they weren't forgotten.

- **`O-17` — how much does a grid-edge root cache buy?** M-98 measured subgrid Marching Tetrahedra at 70× classic MT, and the constant is field evaluations: 576 per cell at 16 samples per edge, against 8 shared corner samples for Marching Cubes. Every cell currently re-finds the roots on edges its neighbours already found, deliberately — identical endpoints through a deterministic root finder is what makes conformity hold without a cache. A cache keyed on the grid edge is the obvious optimisation and the redundancy is large, but it has a correctness precondition and the saving is **unmeasured**. Settle it by caching and re-running `cargo bench --bench extract`, with the golden hashes as the guard that the mesh did not change.
- **`O-18` — the `cycles()` recomputation in `subgrid/surface.rs`.** The 2026-08-14 review proposed restructuring extraction so the cycle set from `cycles()` (`surface.rs:252`, which allocates a `Vec<Cycle>` per call) is computed once rather than recomputed. It sits here rather than as a ticket because it is efficiency churn inside the actively-ticketed A-014 series, and because only the review's one-line summary is on record — re-derive the exact shape from the code before acting. Settle it the way O-17 says to: restructure, `cargo bench --bench extract`, golden hashes as the guard that the mesh did not change.
- Nanite-style mesh-space cluster simplification — the research concludes it can't be repaired
  edit-proportionally (no local validity certificate). Field-derived LOD is the bet instead.
- Networked/concurrent editing — **closed out, not deferred. The verdict is in and it is bounded.**
  O-4 asked whether brush operations commute. They do, conditionally: a run of same-kind *hard* edits
  reorders bit-for-bit — one result from all 40,320 orderings, all `Add` and again all `Subtract`
  (M-36) — because `min`/`max` select an argument rather than computing a value. Across an add/subtract
  boundary they do **not**: 11 distinct results, and the difference is *semantic*, so no storage format
  or arithmetic repairs it (M-37). Smooth union is worse still — 40,317 distinct results from 40,320,
  smooth-min being neither associative nor bit-commutative (M-38).

  So the coordination-free story survives inside a run and dies at every boundary, and **that is a
  protocol's problem, not this crate's.** isomesh's whole obligation was to make the truth available,
  and `BrushOp::commutes_with` already returns the honest answer rather than the optimistic one. A
  networking layer needs sockets, clocks and a session model — none of which belong in a `no_std` crate
  whose public API is `[f32; 3]`. Nothing further is owed here; reopen it as a real ticket only if a
  consumer turns up needing something the existing predicate cannot express.
- Neural / differentiable extraction (FlexiCubes, TetWeave) — different problem, different crate.
- Publishing real releases to crates.io. **`I-005 — reserve the name` is done: `isomesh 0.0.0` was
  published on 2026-08-12, and `0.0.1` on 2026-08-13 to carry the crate's README** and the name is held. `megamesh` was taken 48 hours before we checked it,
  which is what made this urgent rather than tidy. The placeholder is 82 files / 329 KiB compressed —
  source, benches, golden hashes and proptest regressions, nothing stray — and `0.0.0` is now burned
  permanently, which is the intent.

  What stays out of scope is a **real** release. That wants a `crates/isomesh/README.md` (the root one
  is outside the package directory and cannot be referenced from it, so the crates.io page currently
  shows only the one-line description), a version policy, and a decision about whether `isomesh-gpu`
  and `bevy_isomesh` publish alongside it. None of that is urgent now the name cannot be taken.

  **Amended 2026-08-14 (D-003): most of that has since happened.** `crates/isomesh/README.md` exists
  (`d9b8836`) and is the crates.io page; releases are CI-driven on version bumps; `isomesh` and
  `isomesh-gpu` are live at 0.0.4. What remains of this item is a version policy. The `bevy_isomesh`
  decision is made and ticketed: D-007 dresses it, D-010 publishes it at 0.0.4.
