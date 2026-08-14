# isomesh — BACKLOG

**Updated:** 2026-08-12
**Companions:** `CLAUDE.md` (rules), `FINDINGS.md` (what we know and how well),
`BACKLOG_ARCHIVE.md` (completed tickets + why they changed),
`docs/2026-08-11-implementation-brief.md` (the how),
`docs/2026-08-11-bevy-examples-catalog.md` (example detail), `docs/research/` (the why).

**95 tickets archived, 3 open.** Completed rows move to `BACKLOG_ARCHIVE.md` with their amendments
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
- **`grep -ri bevy crates/` returns nothing.** Non-negotiable — see `CLAUDE.md` rule 2.
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

Each algorithm ticket is done when: T-001 reports **no unexplained violations** on all seven test
fields at three resolutions; T-004 determinism passes; T-005 covers it; and a benchmark exists.

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

| | ID | Ticket | Size | Blocked by |
|---|---|---|---|---|
| ☐ | **A-014d** | **Subgrid MT — simplicial embedding (§3.2.3).** Split out of A-014b rather than dropped, because it is a distinct guarantee with its own acceptance. After §3.2.2 the mesh *"has manifold connectivity, but may be an immersed Δ-complex rather than an embedded simplicial complex"* — two adjacent tets sharing a face can each emit an oppositely-oriented copy of a polygon bound by the same non-normal loop. Three polygon types (quad, hexagon, corner pentagon); the pairs form a manifold tube or a punctured sphere, so **connectivity is already right and only *embedding* is at stake.** The fix pushes polygons into the tet interior by inserting midpoints. The paper is explicit that it is conditional — *"needed only when taking the union of the two polygons would yield a nonmanifold edge"* — which makes it a separate ticket rather than a branch inside A-014b. **Re-aimed by M-99, and the target is narrower than it looked.** Appendix A's Theorem A.1 guarantees manifold connectivity in the even-sum case, and Lemma A.2 says the polygons *before* §3.2.3 already form an edge-manifold cell complex. The precondition was checked: **0 of 98,304 tetrahedron faces violate even-sum on any of the seven reference fields**, so the guarantee is in force. Measured accordingly — **unwelded output has 0 non-manifold edges and 0 non-manifold vertices on every field**; the counts below appear only after a weld by position, which cannot tell two coincident-but-distinct polygons apart. So this ticket is **not** "make the output manifold" — it already is, as a complex. It is *separate the coincident geometry far enough that an indexed, welded mesh can still represent it*, which is precisely what the inset does. M-59 is the precedent: a manifold complex an index buffer cannot express, now seen in a second algorithm. **It owns three measured defects rather than a hypothesis (M-97).** At 17³ the subgrid extractor is clean on `sphere`, `torus`, `box_exact` and `thin_plate`, and reports `csg_difference (3 non-manifold edges, 6 non-manifold vertices, 6 inconsistently oriented)`, `fbm_terrain (4, 6, 19)` and `gyroid (0, 0, 138)`. **Re-aimed again on 2026-08-13, and it now owns one row rather than three (M-161, M-163, M-164).** `csg_difference`'s three non-manifold edges are the genuine §3.2.3 row: each is 4 faces but only **3 distinct polygons**, from two sibling tetrahedra in one cell, and removing the duplicated pair leaves exactly 2 faces. **`gyroid`'s 138 moved to A-014f** — it has 1 coincident polygon in 8,088 triangles, and of the 186 triangles on its flipped edges only 15 have the gradient in-plane while **171 have a decisive vote and disagree with their neighbour anyway**, which is A-014e's per-triangle vote, not immersion. And the premise that duplication is the defect is falsified outright: **`box_exact` carries the most coincidence of any field — 30 polygons, 348 triangles standing on their edges — and validates `(0, 0, 0)`.** Those counts are pinned in `the_validity_suite_over_every_reference_field` and in the three instruments named below. **Acceptance:** T-002 reports zero self-intersections on a two-tet fixture built from a non-normal loop that A-014b leaves immersed — but note M-83 first: the counter is blind to folds inside a Steiner fan, so the fixture has to be built from loops whose triangles do *not* all share an apex, or the zero means nothing. **The construction was retrieved in full on 2026-08-12** and is no longer blocked: *"we 'push' polygons into the tet interior: we insert the midpoints of all polygon edges contained in any edge of `f`, and move them a small distance in the inward normal direction. For pentagons, we also insert the midpoint of the edge opposite the distinguished corner. The triangulation patterns in Figure 15, right ensure that inset regions stay within the tetrahedron."* **Two implementation attempts are on record and both were reverted (M-101).** The point-insertion prose is complete, so it was applied first to every boundary-disk region — `box_exact` went from `(0,0,0)` to `(84, 84, 180)` — and then only to loops confined to a single face, the per-tet reading of "two adjacent tetrahedra sharing a face `f`": still `(42, 84, 180)`. **What the failures point at:** a polygon edge lying along an edge of `f` is shared with whatever else meets at that tet edge, so inserting a midpoint detaches the disk from every neighbour that did not also insert one. Which polygons are duplicated is a property of a **pair** of tetrahedra, and this implementation is strictly per-tet. **That question has been measured and answered (M-162): it does, and the neighbours are not all in this cell.** 27 of `csg_difference`'s 33 coincident polygons have other tetrahedra standing on their boundary edges — 312 triangles, up to 12 on one polygon — and **150 of the 312 are in a different cell**, the same 48% as `box_exact`'s 168 of 348. So neither a per-tetrahedron rule nor a per-cell rewrite can carry the condition, and what remains is a scoped architectural decision rather than an open question. Beyond that, Figure 15's patterns themselves are still unread, being a picture rather than prose — so the quad/pentagon/hexagon triangulations need deriving or the figure needs looking at, and rule 5 applies to that step alone. | M | A-014b |
> BLOCKED: **on an architectural decision, not on a missing measurement — the measurement the previous note asked for has been made.** Three instruments now live in `subgrid::extract::tests` and pin what they found: `which_polygons_coincide_across_a_shared_face`, `the_defects_traced_back_to_the_tetrahedra_that_made_them` and `the_surviving_non_manifold_edges_are_not_duplicated_polygons`. **The answer is that the inset needs information no per-tetrahedron and no per-cell rule carries** (M-162): 27 of `csg_difference`'s 33 coincident polygons have foreign tetrahedra on their boundary edges, 312 triangles in total and **150 of them in a different cell**. §3.2.3 moves those edges' midpoints, so each of those triangles must move too or the disk detaches — exactly what M-101's two reverted attempts measured, now with the count that explains them. **What is left to decide is the shape of the change**, and it is a design choice this repo does not make silently: either the extractor gains a cross-cell shared-vertex table so an inset midpoint is inserted once and seen by every triangle on that edge, or the inset becomes a post-pass over the welded mesh, which trades §3.2.3's guarantee for a repair. **Meanwhile the ticket owns exactly one row rather than three** (M-161, M-164): `csg_difference`'s 3 non-manifold edges. `gyroid`'s 138 went to A-014f, and `box_exact` — the field with the most coincidence — was already clean, which is what falsified the premise.
| ☐ | **A-014f** | **Subgrid MT — the orientation vote has no consistency guarantee between neighbours.** Split out of A-014d on measurement, not on taste: A-014d claimed `gyroid`'s `(0, 0, 138)` as *"the same cause wearing a different face"* as §3.2.3's immersion, and M-164 falsified that. `gyroid` has **1 coincident polygon in 8,088 triangles**, so there is essentially nothing there to inset. Of the **186 triangles standing on the 138 inconsistently-oriented edges: 0 are zero-area, 15 have the gradient lying in the triangle's own plane, and 171 have a decisive vote.** So for 92% of them A-014e's `dot(face_normal, gradient(centroid))` returned an unambiguous answer on *both* sides of the edge and the two answers disagreed — two triangles of one sheet, each correctly oriented against the gradient at its own centroid, wound opposite ways. **The tension this ticket has to resolve is A-014e's own rationale.** The vote is per triangle rather than per patch deliberately, because a sheet thinner than a cell puts two oppositely-facing surfaces inside one tetrahedron and a per-patch vote would flip one of them wrongly — that reasoning is sound and `thin_plate` is the case it protects. But a *local* answer that is right at every triangle is still not a *consistent* orientation across the mesh, and `gyroid` is where the difference shows: a triply-periodic surface passes close to itself, so the nearest sheet to a centroid is not always the sheet the triangle is on. The obvious remedy — orient one connected component from its most confident triangle and propagate — is a design choice with a real cost (it needs connectivity before orientation, which the per-tetrahedron soup does not have until it is welded) and it is **not** to be adopted without deciding that trade deliberately. **Acceptance:** `gyroid`'s inconsistently-oriented count driven to 0 in `the_validity_suite_over_every_reference_field` **without** regressing `thin_plate`, whose two sheets 0.4 cells apart are the reason the vote is per triangle in the first place; and the 15 in-plane-gradient triangles given a defined answer rather than an arbitrary one, since a vote of exactly zero currently falls to the `else` branch and keeps §3.2's own winding. | M | A-014e |
> BLOCKED: **the remedy is a design choice, and picking it silently is what this repo's rules forbid.** Orientation-by-propagation needs connectivity the per-tetrahedron soup does not have before welding, so adopting it means either welding earlier or carrying a second structure — and it puts at risk the exact case A-014e's per-triangle vote was written to protect (`thin_plate`'s two sheets inside one tetrahedron). The measurement that scopes it is done (M-164); what is missing is the decision about which guarantee to keep when they conflict.
| ☐ | **A-002b** | **Marching Cubes 33 interior ambiguity — the trilinear body saddle.** Deliberately deferred at A-002, on evidence, not forgotten. Three reasons. (1) `catalog-v2.md:107` is explicit: *"Skip the interior test; spend the budget on chunk seams"* — a game needs topological *consistency*, which A-001 already has, over *correctness*. (2) ~~**There is no correct published table to transcribe.**~~ **Answered on 2026-08-13 by reading the sources (V-24, V-25, V-26), and two thirds of it was wrong.** Custodio et al. 2013 (`10.1016/j.cag.2013.04.004` §5.1) do prove Chernyaev's interior test tracks a quadratic where the true saddle trajectory is hyperbolic with an asymptote — and **that correction is now implemented and tested (A-002c)**, so it is no longer this ticket's. The claim that *"Lewiner's reference implementation omits disambiguation for cases 10 and 12 entirely"* is true of the code and false of the literature: §5.4 calls it *"a missing step in the implementation"* and states the rule for 10.1.1, 10.1.2 and 10.2 outright. And the missing table is not missing but unnecessary — **`10.1186/s13173-019-0086-6` (Custodio, Pesco & Silva 2019) builds the whole MC33 triangulation with no lookup table**, from groups of vertices and edges and the boundaries of their convex hulls, with case 13.3, 13.5.1 and 13.5.2 given as named combinations. Rule 5 is satisfied by following a construction rather than by inventing a table. (3) The v1 catalog prices it: the decider is *"~free"*, the guaranteed version is **730 subcases in the LUT**. Also needs cell-interior vertices for tunnels, which the grid-edge-keyed vertex cache has no slot for. **Acceptance:** a cell where the body saddle says "tunnel", meshed as a tunnel, with the sign tracked by Custodio's correction rather than Chernyaev's `F(t)` — the second half of which A-002c has already done, so what remains is the meshing. **Two things it must settle that were not visible before.** (a) **Crack-freeness is not inherited (M-166).** `ambiguity`'s two cells cannot disagree because its decider compares two *products* and IEEE multiplication is commutative; the interior test's denominator is `((A + C) − B) − D`, whose subtraction order a rotation permutes, and float addition is not associative. Two cells reading one shared face could therefore disagree about a tunnel, and that has to be closed before this is wired in. (b) **The non-manifold remedy is a grid change, not a table change** — both papers say to split the two cells at the shared ambiguous face's critical point, which is a preprocessing pass over the grid rather than anything in the extractor. | L | A-002 |
> BLOCKED: **on scope, not on a source — the rule-5 stop is lifted (V-24, V-25, V-26).** The previous note said *"there is no correct published table to transcribe"*; reading the papers showed the decider is stated in full, the case-10 rule is stated in full, and the 2019 follow-up builds the triangulation **without any table at all**. What is left is genuinely large rather than genuinely blocked: MC33's tunnel cases need vertices in the cell **interior**, and this crate's vertex cache is keyed on grid edges with no slot for one; the non-manifold case needs a grid-subdivision pass; and M-166 has to be closed first, because two cells can currently disagree about a tunnel across a shared face. `catalog-v2.md:107` still prices all of that as *"skip the interior test; spend the budget on chunk seams"* for a game, so this stays where it is in the queue — but it is now a size decision rather than a correctness one.

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

| | ID | Ticket | Size | Blocked by |
|---|---|---|---|---|

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

| | ID | Ticket | Size | Blocked by |
|---|---|---|---|---|
> NOT BLOCKED, and an earlier note here said it was — see M-147. **The route needs no `unsafe` in this repository at all.** `isomesh-gpu` never opens a device (its API takes `&wgpu::Device`, GPU-001's rule), and **Bevy writes the experimental token itself**: `experimental_features: unsafe { wgpu::ExperimentalFeatures::enabled() }` at `bevy_render-0.19.0/src/renderer/mod.rs:335`. `WgpuSettings`' default priority is `Functionality`, which requests every feature the adapter advertises, so **Bevy's device already reports `mesh_shader=true multiview=true points=true`** on this machine, measured. E-303 is a Bevy example and gets a mesh-shader-capable device for free; `WgpuSettings.features` is there to force it explicitly if the default ever changes.
>
> **The probe is load-bearing, not belt-and-braces** — an earlier version of this note implied otherwise and was wrong in the opposite direction from the blocked claim it replaced. The free device is *one of three branches*: `WgpuSettings::default()` consults `settings_priority_from_env()` first, so **`WGPU_SETTINGS_PRIO` overrides it**; under any priority other than `Functionality`, `features` starts at `wgpu::Features::empty()`; and `adapter.features()` is machine-dependent. It is also **contingent upstream**: Bevy's line carries `// SAFETY: TODO, see bevyengine/bevy#22082`, an admission that a justification is owed, so if that issue lands as opt-in the default path loses mesh shaders. Track it.
>
> **On "graceful fallback", which needs one distinction rather than a ruling.** A demo that *detects* capability and, finding none, says so plainly and shows the compute path instead is a demo reporting a capability — one path, chosen by a measurement, with the choice visible. What the one-path rule forbids is the *library* silently substituting compute for mesh shaders so a caller cannot tell which ran. The first is what this ticket should build.
>
> **Two things still shape the work.** WGSL mesh shaders are **Vulkan-only** — wgpu's own source says *"naga is only supported on vulkan; on other platforms you will have to use passthrough shaders"* (V-23) — so on Metal a caller supplies pre-compiled MSL and the composed-WGSL pipeline does not apply, making mesh shaders a fork in the shader path rather than a flag on it. And the ticket's own wording needs revisiting: *"graceful fallback"* is a second execution path for one feature, which the one-path rule forbids; the shape that survives is a **capability check that refuses loudly**, as GPU-007's probe already does. The remaining `unsafe`-shaped gap is only `isomesh-gpu::headless::Gpu` opening its *own* mesh-shader device, which is a test convenience rather than this ticket.

---

## Deliberately not in scope yet

Recorded so they don't get picked up early, and so it's clear they weren't forgotten.

- **`O-17` — how much does a grid-edge root cache buy?** M-98 measured subgrid Marching Tetrahedra at 70× classic MT, and the constant is field evaluations: 576 per cell at 16 samples per edge, against 8 shared corner samples for Marching Cubes. Every cell currently re-finds the roots on edges its neighbours already found, deliberately — identical endpoints through a deterministic root finder is what makes conformity hold without a cache. A cache keyed on the grid edge is the obvious optimisation and the redundancy is large, but it has a correctness precondition and the saving is **unmeasured**. Settle it by caching and re-running `cargo bench --bench extract`, with the golden hashes as the guard that the mesh did not change.
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
