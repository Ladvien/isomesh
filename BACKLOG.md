# isomesh — BACKLOG

**Updated:** 2026-08-14
**Companions:** `CLAUDE.md` (rules), `FINDINGS.md` (what we know and how well),
`BACKLOG_ARCHIVE.md` (completed tickets + why they changed),
`docs/2026-08-11-implementation-brief.md` (the how),
`docs/2026-08-11-bevy-examples-catalog.md` (example detail), `docs/research/` (the why).

**110 tickets archived, 10 open.** Completed rows move to `BACKLOG_ARCHIVE.md` with their amendments
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
| ☐ | **A-002b** | **Marching Cubes 33 interior ambiguity — the trilinear body saddle.** Deliberately deferred at A-002, on evidence, not forgotten. Three reasons. (1) `catalog-v2.md:107` is explicit: *"Skip the interior test; spend the budget on chunk seams"* — a game needs topological *consistency*, which A-001 already has, over *correctness*. (2) ~~**There is no correct published table to transcribe.**~~ **Answered on 2026-08-13 by reading the sources (V-24, V-25, V-26), and two thirds of it was wrong.** Custodio et al. 2013 (`10.1016/j.cag.2013.04.004` §5.1) do prove Chernyaev's interior test tracks a quadratic where the true saddle trajectory is hyperbolic with an asymptote — and **that correction is now implemented and tested (A-002c)**, so it is no longer this ticket's. The claim that *"Lewiner's reference implementation omits disambiguation for cases 10 and 12 entirely"* is true of the code and false of the literature: §5.4 calls it *"a missing step in the implementation"* and states the rule for 10.1.1, 10.1.2 and 10.2 outright. And the missing table is not missing but unnecessary — **`10.1186/s13173-019-0086-6` (Custodio, Pesco & Silva 2019) builds the whole MC33 triangulation with no lookup table**, from groups of vertices and edges and the boundaries of their convex hulls, with case 13.3, 13.5.1 and 13.5.2 given as named combinations. Rule 5 is satisfied by following a construction rather than by inventing a table. (3) The v1 catalog prices it: the decider is *"~free"*, the guaranteed version is **730 subcases in the LUT**. Also needs cell-interior vertices for tunnels, which the grid-edge-keyed vertex cache has no slot for. **Acceptance:** a cell where the body saddle says "tunnel", meshed as a tunnel, with the sign tracked by Custodio's correction rather than Chernyaev's `F(t)` — the second half of which A-002c has already done, so what remains is the meshing. **Three things it must settle that were not visible before.** (a) **Crack-freeness is not inherited (M-166).** `ambiguity`'s two cells cannot disagree because its decider compares two *products* and IEEE multiplication is commutative; the interior test's denominator is `((A + C) − B) − D`, whose subtraction order a rotation permutes, and float addition is not associative. Two cells reading one shared face could therefore disagree about a tunnel, and that has to be closed before this is wired in. (b) **The non-manifold remedy is a grid change, not a table change** — both papers say to split the two cells at the shared ambiguous face's critical point, which is a preprocessing pass over the grid rather than anything in the extractor. (c) **Added by the 2026-08-14 review:** 2 of M-165's 15,625 opposed-family configurations rest a root of `F` on Δ's pole to within last ulps, so the decider's answer there is rounding noise rather than a decision; the meshing pass owes those configurations a defined answer, not an inherited one. | L | A-002 |
> BLOCKED: **on scope, not on a source — the rule-5 stop is lifted (V-24, V-25, V-26).** The previous note said *"there is no correct published table to transcribe"*; reading the papers showed the decider is stated in full, the case-10 rule is stated in full, and the 2019 follow-up builds the triangulation **without any table at all**. What is left is genuinely large rather than genuinely blocked: MC33's tunnel cases need vertices in the cell **interior**, and this crate's vertex cache is keyed on grid edges with no slot for one; the non-manifold case needs a grid-subdivision pass; and ~~M-166 has to be closed first, because two cells can currently disagree about a tunnel across a shared face~~ — **M-166 is closed (M-204).** The denominator was `((A + C) − B) − D`, a subtraction order a rotation permutes; grouping each diagonal as `(A + C) − (B + D)` makes every symmetry of the face exact, and two cells reading one face now agree bit-for-bit. Note the test already asserting that property **passed with the defect present** — its fixture's corners are all of similar magnitude — so the fixture that exposes it had to be searched for (`(1, 1, 1, 10⁻⁸)`), and the new test is mutation-verified in both directions. **So this ticket is now blocked on size alone**, which is where its own note already put it. `catalog-v2.md:107` still prices all of that as *"skip the interior test; spend the budget on chunk seams"* for a game, so this stays where it is in the queue — but it is now a size decision rather than a correctness one.

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
| ☐ | **D-002** | **Remove the `gpu` feature that gates nothing.** `bevy_isomesh` declares `gpu = ["dep:bevy_render"]` and no `cfg(feature = "gpu")` exists in src or examples — the three GPU examples reach render types through the `bevy` umbrella dev-dependency, verified by grep. A feature that changes nothing is a second path that lies about the first. Delete the feature, the optional `bevy_render` dependency, and the stale comment claiming the plugin hasn't landed yet. | S | — |
| ☐ | **D-003** | **Truth pass, claims.** The root README says the crates.io name is a 0.0.2 placeholder; crates.io serves five isomesh releases and two isomesh-gpu releases, CI-published. It says Rust 1.85 twice across the READMEs; the workspace pins 1.89 and `bevy_isomesh` 1.95. It says 147 golden combinations; the fixture holds 168. Running totals that nothing gates (ticket counts, test counts, demo counts) become pointers at their gated sources or the commands that print them; measured-once figures keep their numbers. Add the three GPU examples to the run list. Add a **Machines** block to `FINDINGS.md`'s header naming both benchmark machines and which CSV each produced — the root README already claims that block exists, which today it does not. Refresh this file's stale deferred-publishing paragraph (the crate README it says is missing has existed since `d9b8836`). Make every FINDINGS.md name-check an actual link, absolute in the crate READMEs. | M | — |
| ☐ | **D-004** | **Truth pass, demo surfaces.** The examples catalog says "26 examples"; 31 live in `bevy_isomesh/examples/` plus 3 headless in `crates/isomesh-gpu/examples/` — correct it with a dated addendum rather than rewriting a 2026-08-11 planning document. `docs/demos/gameplay.md` gets the top back-link the other two demo pages have, and its LOD open-edge figures predate the E-211/M-195 counter fix (the taper moves cracks off the seam plane; the old counter read 0 there) — reconcile the section with the corrected instrument's numbers. Normalize the four module headers that lead with a ticket ID (`greedy_quads`, `marching_tetrahedra`, `weld`, `bevy_isomesh::plugin`) to reader-facing first lines; the ID moves to a "Ticket:" line below, because the first line is the docs.rs module summary. | M | — |
| ☐ | **D-005** | **`crates/isomesh/README.md` — the chooser's page gets a chooser's tools.** Add the minimal badge row (crates.io, docs.rs, CI, license). Add a "Choosing an extractor" table: seven rows, sharp features / chunk tiling / cost, sourced from `Extractor::chunk_seams` and `docs/measurements/shootout.csv`. Add two troubleshooting entries — the debug build (M-152's 37–62×) and zero triangles (negative-is-inside, or the domain misses the surface) — because error recovery is content, not a warning. The `cfg(doctest)` README wiring stays exactly as it is; its quickstart remains the canonical snippet other pages copy. | S | D-003 |
| ☐ | **D-006** | **`crates/isomesh-gpu/README.md` — the snippet nothing compiles.** Wire the README into doctests with the same `cfg(doctest)` `include_str!` block core uses; the snippet's fence becomes `rust,no_run` with a hidden `Ok` tail, so a GPU-less CI runner compiles it and never runs it. Add the badge row, one screenshot by absolute URL, the License section this README never had, and fix the `(../isomesh)` link that 404s on crates.io. Explicit `readme` and `documentation` manifest keys here, plus `documentation` for core in the same commit. `[package.metadata.docs.rs]` is deliberately not added anywhere: zero cargo features means there is nothing true to configure, and an empty stanza is a claim-shaped nothing. | M | D-003 |
| ☐ | **D-007** | **`bevy_isomesh/README.md` — rebuilt as the front page it is about to become.** Full rewrite to the `docs/bevy_plugins.md` rules. Bevy compat table. A `rust,no_run` quickstart doctest: `DefaultPlugins` + `IsomeshPlugin`, spawn a `VoxelVolume` and its chunks, attach `Mesh3d` yourself — that last step is the design boundary, stated as such. A runnable `MeshBuilder` snippet (CPU-only, so it actually executes). The exposed contract, written down at last: four components, two resources, `(spawn_meshing_tasks, apply_finished_meshes).chain()` in `Update`, no public SystemSet yet, no cargo features. The chunk-seam table (dual methods are `Gapped`, structurally — a chunked quickstart must use MarchingCubes or Subgrid). The why-a-separate-workspace section absorbs `src/lib.rs`'s current `//!` docs, because `src/lib.rs` becomes `#![doc = include_str!("../README.md")]` — docs.rs renders this page and `cargo test` compiles its fences, which is why every URL in it must be absolute and no intra-doc bracket shorthand may appear. A curated eight-example table with screenshots, then the pointer to all 31. Copy `LICENSE-MIT` and `LICENSE-APACHE` into the package — the README promises them, and a `license` field with no file beside it already shipped three isomesh releases before anyone noticed. Manifest: `readme`, `keywords`, `categories`, `documentation`. A public `IsomeshSystems` set is a code change and a future ticket, not this one. | M | D-002 |
| ☐ | **D-008** | **Root README — the on-ramp above the essay.** Badge row. "In 60 seconds": `cargo add isomesh`, the canonical snippet copied verbatim from `crates/isomesh/README.md`, three next-step links. An honest fit table with a when-not-to-use list: dual methods don't tile chunks, GPU readback is not a speedup, MC33 tunnels unmeshed. The claim-headline essays follow in their current order, tightened where D-003 touched them — the voice is the differentiation and it stays. The snippet copy is gated: new `scripts/readme_sync.sh` diffs the `rust` fences of the two files, and CI's lint job runs it beside the backlog gate. Chosen over a second `include_str!` because the root README is outside the package directory, and a published crate must never reference a file it does not ship. | M | D-003, D-005 |
| ☐ | **D-009** | **CHANGELOG and templates.** `CHANGELOG.md` in Keep a Changelog form, backfilled from the five version-bump commits (0.0.0 reserved the name on 2026-08-12; 0.0.3 shipped the license files three releases forgot; 0.0.4 shipped isomesh-gpu's demos). Issue templates — the bug form asks for the exact command, whether `--release`, and the machine, because machine-naming is a house value — and a PR template that asks whether the gates ran, whether any perf claim has its committed benchmark, and whether FINDINGS changed. Plain markdown, no YAML forms. CONTRIBUTING and a code of conduct were considered and declined on 2026-08-14 — recorded here so it isn't re-litigated. | S | — |
| ☐ | **D-010** | **`bevy_isomesh` 0.0.4 — the release wiring. Merging this publishes.** Version to 0.0.4, matching the workspace train. `isomesh = { path = "../crates/isomesh", version = "0.0.4" }`, because a path-only dependency cannot publish. `scripts/publish.sh` gains an explicit bevy_isomesh stanza — the script's member discovery runs `cargo metadata` from the root workspace, which excludes this crate, so the stanza is spelled out the way `ORDER` is, same 200/404 logic, after isomesh-gpu. The CI package job dry-runs it. CHANGELOG entry. `bevy_isomesh 0.0.4` burns permanently on the next push to main; that is the intent, decided 2026-08-14. | S | D-007 |

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
