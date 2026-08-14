# isomesh — BACKLOG

**Updated:** 2026-08-12
**Companions:** `CLAUDE.md` (rules), `FINDINGS.md` (what we know and how well),
`BACKLOG_ARCHIVE.md` (completed tickets + why they changed),
`docs/2026-08-11-implementation-brief.md` (the how),
`docs/2026-08-11-bevy-examples-catalog.md` (example detail), `docs/research/` (the why).

**93 tickets archived, 2 open.** Completed rows move to `BACKLOG_ARCHIVE.md` with their amendments
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
| ☐ | **A-014d** | **Subgrid MT — simplicial embedding (§3.2.3).** Split out of A-014b rather than dropped, because it is a distinct guarantee with its own acceptance. After §3.2.2 the mesh *"has manifold connectivity, but may be an immersed Δ-complex rather than an embedded simplicial complex"* — two adjacent tets sharing a face can each emit an oppositely-oriented copy of a polygon bound by the same non-normal loop. Three polygon types (quad, hexagon, corner pentagon); the pairs form a manifold tube or a punctured sphere, so **connectivity is already right and only *embedding* is at stake.** The fix pushes polygons into the tet interior by inserting midpoints. The paper is explicit that it is conditional — *"needed only when taking the union of the two polygons would yield a nonmanifold edge"* — which makes it a separate ticket rather than a branch inside A-014b. **Re-aimed by M-99, and the target is narrower than it looked.** Appendix A's Theorem A.1 guarantees manifold connectivity in the even-sum case, and Lemma A.2 says the polygons *before* §3.2.3 already form an edge-manifold cell complex. The precondition was checked: **0 of 98,304 tetrahedron faces violate even-sum on any of the seven reference fields**, so the guarantee is in force. Measured accordingly — **unwelded output has 0 non-manifold edges and 0 non-manifold vertices on every field**; the counts below appear only after a weld by position, which cannot tell two coincident-but-distinct polygons apart. So this ticket is **not** "make the output manifold" — it already is, as a complex. It is *separate the coincident geometry far enough that an indexed, welded mesh can still represent it*, which is precisely what the inset does. M-59 is the precedent: a manifold complex an index buffer cannot express, now seen in a second algorithm. **It owns three measured defects rather than a hypothesis (M-97).** At 17³ the subgrid extractor is clean on `sphere`, `torus`, `box_exact` and `thin_plate`, and reports `csg_difference (3 non-manifold edges, 6 non-manifold vertices, 6 inconsistently oriented)`, `fbm_terrain (4, 6, 19)` and `gyroid (0, 0, 138)`. The first two are the pinching §3.2.3 names. **`gyroid`'s is the same cause wearing a different face**: 12 of its triangles have zero area and 24 have the gradient lying *exactly in the triangle's plane*, so A-014e's orientation vote is undecided rather than wrong — insetting gives those polygons area and a normal transverse to the gradient, which fixes both symptoms at once. Those counts are pinned in `the_validity_suite_over_every_reference_field`, so this ticket's real acceptance is **driving all three rows to zero**. **Acceptance:** T-002 reports zero self-intersections on a two-tet fixture built from a non-normal loop that A-014b leaves immersed — but note M-83 first: the counter is blind to folds inside a Steiner fan, so the fixture has to be built from loops whose triangles do *not* all share an apex, or the zero means nothing. **The construction was retrieved in full on 2026-08-12** and is no longer blocked: *"we 'push' polygons into the tet interior: we insert the midpoints of all polygon edges contained in any edge of `f`, and move them a small distance in the inward normal direction. For pentagons, we also insert the midpoint of the edge opposite the distinguished corner. The triangulation patterns in Figure 15, right ensure that inset regions stay within the tetrahedron."* **Two implementation attempts are on record and both were reverted (M-101).** The point-insertion prose is complete, so it was applied first to every boundary-disk region — `box_exact` went from `(0,0,0)` to `(84, 84, 180)` — and then only to loops confined to a single face, the per-tet reading of "two adjacent tetrahedra sharing a face `f`": still `(42, 84, 180)`. **What the failures point at:** a polygon edge lying along an edge of `f` is shared with whatever else meets at that tet edge, so inserting a midpoint detaches the disk from every neighbour that did not also insert one. Which polygons are duplicated is a property of a **pair** of tetrahedra, and this implementation is strictly per-tet. So the first open question is no longer "what are Figure 15's diagonals" but **"which polygons does the inset apply to, and does answering that need neighbour information the per-tet architecture does not carry"** — and if it does, this is an architectural change rather than a triangulation table. Beyond that, Figure 15's patterns themselves are still unread, being a picture rather than prose — so the quad/pentagon/hexagon triangulations need deriving or the figure needs looking at, and rule 5 applies to that step alone. | M | A-014b |
> BLOCKED: **which polygons the inset applies to.** The point-insertion prose is complete and was implemented twice; both attempts made the measured result *worse* and both were reverted (M-101). Applied to every boundary-disk region, `box_exact` goes from a clean `(0,0,0)` to `(84, 84, 180)`; restricted to loops confined to one face — the per-tet reading of "two adjacent tetrahedra sharing a face `f`" — it is still `(42, 84, 180)`. A polygon edge lying along an edge of `f` is shared with whatever else meets at that tet edge, so inserting a midpoint detaches the disk from every neighbour that did not also insert one, and **which polygons are duplicated is a property of a pair of tetrahedra** while this implementation is strictly per-tet. Unblocking it needs one measurement, not another attempt: instrument the extractor to report which polygons actually coincide across a shared face on `csg_difference`, and see whether that condition is expressible without neighbour information. If it is not, this is an architectural change rather than a triangulation table. **The three defects it owns are pinned and harmless meanwhile** (M-97, M-99): connectivity is provably manifold under Theorem A.1, and the counts only appear after a weld.
| ☐ | **A-002b** | **Marching Cubes 33 interior ambiguity — the trilinear body saddle.** Deliberately deferred at A-002, on evidence, not forgotten. Three reasons. (1) `catalog-v2.md:107` is explicit: *"Skip the interior test; spend the budget on chunk seams"* — a game needs topological *consistency*, which A-001 already has, over *correctness*. (2) **There is no correct published table to transcribe.** Custodio et al. 2013 (`10.1016/j.cag.2013.04.004` §5.1) prove Chernyaev's interior test tracks a quadratic where the true saddle trajectory is hyperbolic with an asymptote, so case 13.5.2 is misread as 13.5.1 — counterexample values in their Appendix A — and Lewiner's reference implementation omits disambiguation for cases 10 and 12 entirely. Rule 5 forbids inventing the missing one. (3) The v1 catalog prices it: the decider is *"~free"*, the guaranteed version is **730 subcases in the LUT**. Also needs cell-interior vertices for tunnels, which the grid-edge-keyed vertex cache has no slot for. **Acceptance:** a cell where the body saddle says "tunnel", meshed as a tunnel, with the sign tracked by Custodio's correction rather than Chernyaev's `F(t)`. | L | A-002 |
> BLOCKED: **there is no correct published table to transcribe**, which is stated on the row and is a rule-5 stop rather than a difficulty. Custodio et al. prove Chernyaev's interior test tracks a quadratic where the true saddle trajectory is hyperbolic, and Lewiner's reference implementation omits cases 10 and 12 entirely. Unblocking it means deriving the corrected decider from Custodio §5.1 and its Appendix A counterexamples — real work with a real source, not a transcription — and `catalog-v2.md:107` prices it as *"skip the interior test; spend the budget on chunk seams"* for a game.

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
