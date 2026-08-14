# isomesh — BACKLOG

**Updated:** 2026-08-14
**Companions:** `CLAUDE.md` (rules), `FINDINGS.md` (what we know and how well),
`BACKLOG_ARCHIVE.md` (completed tickets + why they changed),
`docs/2026-08-11-implementation-brief.md` (the how),
`docs/2026-08-11-bevy-examples-catalog.md` (example detail), `docs/research/` (the why).

**102 tickets archived, 8 open.** Completed rows move to `BACKLOG_ARCHIVE.md` with their amendments
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
| ☐ | **A-014d** | **Subgrid MT — §3.2.3's simplicial embedding, now a per-cell problem on three edges of one field.** **Unblocked by A-014h, which falsified the measurement that blocked it (M-186).** M-162 had this ticket stopped on the finding that *"the inset needs information no per-tetrahedron and no per-cell rule carries"* — 33 coincident polygons on `csg_difference`, **27 with foreign tetrahedra on their boundary edges, 312 such triangles and 150 of them in a different cell**. Re-measured once every crossing had a complete identity: **3 coincident polygons, 0 with foreign edge users, 0 cross-cell, all three inside one cell.** `fbm_terrain` has 4, likewise all within one cell; every other field including `box_exact` is at **0** — and `box_exact` is the field M-161 called *"the most coincidence of any field, 30 polygons and 348 triangles standing on their edges"*. So most of what M-162 measured was never §3.2.3's immersion at all: it was one grid point wearing many names, and it went when the names were fixed rather than when any geometry was inset. **The target is `csg_difference`'s 3 non-manifold edges**, each 4 faces but only **3 distinct** polygons from two sibling tetrahedra in one cell (M-163) — and the 3 *duplication-only* edges that used to sit beside them are already gone, along with `thin_plate`'s 2. **What is left to do:** the inset itself, per §3.2.3, over a pair of tetrahedra sharing a face inside one cell. **Rule 5 still applies to one step** — Figure 15's quad/pentagon/hexagon triangulation patterns are a picture rather than prose, so they must be derived or the figure read, and **not invented**; stop and say so if they cannot be established. The point-insertion prose itself is on record in full: *"we insert the midpoints of all polygon edges contained in any edge of `f`, and move them a small distance in the inward normal direction. For pentagons, we also insert the midpoint of the edge opposite the distinguished corner."* **Acceptance:** T-002 reports zero self-intersections on a two-tet fixture built from a non-normal loop — heeding M-83, whose warning is that the counter is blind to folds inside a Steiner fan, so the fixture must use loops whose triangles do **not** all share an apex or the zero means nothing; and `csg_difference`'s `(3, 6, 6)` row in `the_validity_suite_over_every_reference_field` goes to `(0, 0, 0)` with no other row moving. | M | — |
| ☐ | **A-002b** | **Marching Cubes 33 interior ambiguity — the trilinear body saddle.** Deliberately deferred at A-002, on evidence, not forgotten. Three reasons. (1) `catalog-v2.md:107` is explicit: *"Skip the interior test; spend the budget on chunk seams"* — a game needs topological *consistency*, which A-001 already has, over *correctness*. (2) ~~**There is no correct published table to transcribe.**~~ **Answered on 2026-08-13 by reading the sources (V-24, V-25, V-26), and two thirds of it was wrong.** Custodio et al. 2013 (`10.1016/j.cag.2013.04.004` §5.1) do prove Chernyaev's interior test tracks a quadratic where the true saddle trajectory is hyperbolic with an asymptote — and **that correction is now implemented and tested (A-002c)**, so it is no longer this ticket's. The claim that *"Lewiner's reference implementation omits disambiguation for cases 10 and 12 entirely"* is true of the code and false of the literature: §5.4 calls it *"a missing step in the implementation"* and states the rule for 10.1.1, 10.1.2 and 10.2 outright. And the missing table is not missing but unnecessary — **`10.1186/s13173-019-0086-6` (Custodio, Pesco & Silva 2019) builds the whole MC33 triangulation with no lookup table**, from groups of vertices and edges and the boundaries of their convex hulls, with case 13.3, 13.5.1 and 13.5.2 given as named combinations. Rule 5 is satisfied by following a construction rather than by inventing a table. (3) The v1 catalog prices it: the decider is *"~free"*, the guaranteed version is **730 subcases in the LUT**. Also needs cell-interior vertices for tunnels, which the grid-edge-keyed vertex cache has no slot for. **Acceptance:** a cell where the body saddle says "tunnel", meshed as a tunnel, with the sign tracked by Custodio's correction rather than Chernyaev's `F(t)` — the second half of which A-002c has already done, so what remains is the meshing. **Three things it must settle that were not visible before.** (a) **Crack-freeness is not inherited (M-166).** `ambiguity`'s two cells cannot disagree because its decider compares two *products* and IEEE multiplication is commutative; the interior test's denominator is `((A + C) − B) − D`, whose subtraction order a rotation permutes, and float addition is not associative. Two cells reading one shared face could therefore disagree about a tunnel, and that has to be closed before this is wired in. (b) **The non-manifold remedy is a grid change, not a table change** — both papers say to split the two cells at the shared ambiguous face's critical point, which is a preprocessing pass over the grid rather than anything in the extractor. (c) **Added by the 2026-08-14 review:** 2 of M-165's 15,625 opposed-family configurations rest a root of `F` on Δ's pole to within last ulps, so the decider's answer there is rounding noise rather than a decision; the meshing pass owes those configurations a defined answer, not an inherited one. | L | A-002 |
> BLOCKED: **on scope, not on a source — the rule-5 stop is lifted (V-24, V-25, V-26).** The previous note said *"there is no correct published table to transcribe"*; reading the papers showed the decider is stated in full, the case-10 rule is stated in full, and the 2019 follow-up builds the triangulation **without any table at all**. What is left is genuinely large rather than genuinely blocked: MC33's tunnel cases need vertices in the cell **interior**, and this crate's vertex cache is keyed on grid edges with no slot for one; the non-manifold case needs a grid-subdivision pass; and M-166 has to be closed first, because two cells can currently disagree about a tunnel across a shared face. `catalog-v2.md:107` still prices all of that as *"skip the interior test; spend the budget on chunk seams"* for a game, so this stays where it is in the queue — but it is now a size decision rather than a correctness one.
| ☐ | **A-014i** | **Subgrid MT — the subdivision path's three recorded defects.** The 2026-08-14 review found three defects in §3.2.1 case (5)'s subdivision path (`subgrid/surface.rs`, `fill_subdivision` and the `Subdivision` stencil): orphaned vertices, a child re-lerp that should use `1 − t`, and orphaned bigons. It fixed none of them, for a measured reason: **the path is unreachable on every reference field** — `subdivision == 0` is pinned on all seven — so no current output is wrong, and the bigon fix additionally hinges on A-014d's inset-geometry decision. The defect names above are the review's one-line summaries, not verified re-derivations, so rule 5 applies: **the first task is a fixture that actually forces case (5)** (`coordinates.rs`: subdivision requires ≥ 24 edge intersections) and reproduces each defect as a failing test, before any fix is written. Then the standard gates — T-001 validity on the fixture, determinism — with the pinned `subdivision == 0` counts left in place as the guard that reference-field output did not change. | M | A-014d |
| ☐ | **I-007** | **`bevy_isomesh`'s MSRV claim is asserted, not checked — the same hole M-174 just closed for the root.** The `msrv` job reads `cargo metadata --no-deps` from the repository root, and the root workspace *excludes* `bevy_isomesh` by design (feature unification, see `CLAUDE.md`). So the `rust-version = "1.95"` that landed with M-174 is derived from Bevy 0.19's own floor and nothing enforces it — exactly the class of claim M-174 found had quietly stopped being true for the root. This is the third hole the excluded workspace has produced after rustfmt and rustdoc, which is itself the finding: **every check that reads the root workspace needs a deliberate answer for `bevy_isomesh`, and "it is excluded" is not one.** Fix by giving the `msrv` job a second leg that reads `bevy_isomesh/Cargo.toml` and runs `cargo +$MSRV check --all-targets` in that directory, the same shape the existing leg uses. **Acceptance:** the claimed version is what the job installs, and lowering `rust-version` by one minor makes the job go red. | S | — |

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
| ☐ | **B-008** | **Pool the plugin's per-chunk scratch allocations.** Recorded (and deliberately not applied) by the 2026-08-14 review: `bevy_isomesh/src/mesh.rs`'s shading path allocates a scratch buffer per call to split the borrow ("shade into a scratch buffer, then swap"), which at re-meshing scale is an allocation per chunk per edit — the pattern core rule 6 exists to prevent, one workspace over. Keep the scratch in a plugin-owned resource and reuse it across chunks. Any perf claim needs a committed measurement; without one this is an allocation-count cleanup and the commit should say so. | S | — |

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
| ☐ | **E-211** | **`game_lod_flyover` — the HUD crack counter must scan the inset band, not just the seam plane.** The counter looks for boundary edges only on the seam plane itself, but Eq-4.2's taper places transition geometry at seam ± w — which is exactly where the 2026-08-14 review's mirrored-seam bug lived: 44 boundary edges on a genuinely open seam, invisible to the example's own instrument. | The counter reports 0 on the fixed build and goes non-zero if the mirrored-seam signed-width fix is reverted — i.e. it would now catch the bug it missed. | — |
| ☐ | **E-212** | **The game examples consume `fields::`, not hand-rolled SDF copies.** `game_destruction`, `game_csg_props` and `game_showcase` each re-implement field combinators the crate already ships, without the analytic gradients the real ones carry. The 2026-08-14 review skipped this as a behavior-changing demo refactor, correctly. Note `game_destruction` had concurrent active work land mid-review (E-204) — rebase on it, don't fight it. | The demos read as consumers of the crate's field API rather than private forks of it, meshes still look right by eye, and recorded framerates hold. | — |

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
| ☐ | **GPU-013** | **Two recorded cleanups in `isomesh-gpu`: the duplicated entry-point prologue, and per-pass submits.** From the 2026-08-14 review, skipped there as perf/architecture churn inside actively-ticketed code. (1) `extract` and `extract_indirect` duplicate their setup prologue; extract it once so the two entry points cannot drift — the indirect-draw budget-clamp fix touched exactly this duplicated region. (2) The extraction path submits per pass where one batched submission may do; locate the actual submit sites first, then decide. M-159/M-167 are the reminder that synchronisation, not arithmetic, is where this crate's time goes — so if batching claims a win, commit the measurement that shows it, and if it measures nothing, record that and close the second half as answered. | S | — |
> NOT BLOCKED — this note is about the mesh-shader route (E-303), not the GPU-013 row above it, and an earlier version of it said the route was blocked — see M-147. **The route needs no `unsafe` in this repository at all.** `isomesh-gpu` never opens a device (its API takes `&wgpu::Device`, GPU-001's rule), and **Bevy writes the experimental token itself**: `experimental_features: unsafe { wgpu::ExperimentalFeatures::enabled() }` at `bevy_render-0.19.0/src/renderer/mod.rs:335`. `WgpuSettings`' default priority is `Functionality`, which requests every feature the adapter advertises, so **Bevy's device already reports `mesh_shader=true multiview=true points=true`** on this machine, measured. E-303 is a Bevy example and gets a mesh-shader-capable device for free; `WgpuSettings.features` is there to force it explicitly if the default ever changes.
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
