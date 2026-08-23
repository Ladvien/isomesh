# Backlog and harness audit — the three axes

**Date:** 2026-08-18
**Against:** `origin/main` at `db87e56` (Phase 18 closed, M-335/P-37 HELD), 212 tickets archived, 14 open.
**Tier:** structural claims are **M** — read from the tree this session, with file paths. Judgements about what *would* help are **R**.

---

## 0. The finding that reframes the other three

**Every one of the 14 open tickets is blocked, and 4 of the blocks are decisions only you can make.**

| Ticket | Blocked on |
|---|---|
| A-026 | **Your decision** — how far the CoACD pipeline goes in-crate |
| A-027, T-022b, B-013 | A-026 |
| T-022a | Nothing, but its architecture is downstream of A-026 |
| A-025 | **Your decision** — the `FaceAmbiguity` default and its hash re-baseline |
| X-005 | **Your decision** — take the API break or write the contract down |
| R-027 | Stopped on V-45; **blocked on R-027a** |
| R-027a | Nothing. **`S`. The only genuinely unblocked ticket in the file.** |
| R-021 | R-020, and re-blocked on V-44 |
| R-026 | Nothing structural — a writeup |
| A-002i, A-020b | Architecture, on measurements already in |
| M-005 | The M5 being quiet |

Rule 1 says take the topmost unblocked unchecked ticket. Applied honestly to this file today, that is **R-027a** — an `S`, no new field, no new extractor, instrument the existing `edit_trace` bench. Everything above it in the file is waiting on you.

The companion memo (`2026-08-18-four-blocked-decisions.md`) works the four decisions with the evidence for each. **Nothing else in this audit matters as much as unblocking those**, because the queue is currently a queue of one.

---

## Axis 1 — Preparedness for novel experimentation

### What is genuinely strong, and rare

The pre-registration machinery is real, not aspirational. `crates/isomesh/src/experiment.rs` holds **30 `Preregistration` entries (P-8 … P-37)** and `experiment!("P-n")` is a **compile-time** assertion that `n` is in `PREREGISTERED` (experiment.rs:1128-1140). A harness cannot be written for an unregistered hypothesis; the build refuses. Naming is 1:1 — `P-n` → `benches/experiment_pN.rs` → `docs/experiments/p-N.csv` — and `benches/common/experiment.rs`'s `Run::record` **panics if a column the registration promised is not written**. That is a promise enforced by the machine rather than by discipline.

`scripts/backlog_gate.sh` checks that every `P-` in `experiment.rs` is mentioned in `FINDINGS.md`, that no ticket sits in both BACKLOG and ARCHIVE, that header counts match row counts, and that no "Blocked by" names a ticket that does not exist.

**The mutation check is the best thing in the repo and nobody would notice it.** `property/extraction.rs:1-30` re-implements Marching Cubes reading the real case table, proves it byte-identical to the production marcher, then runs it against a *corrupted* table to demonstrate that the manifoldness property test can actually go red. **That is a negative control on a test.** Almost no crate has one.

### The gaps, in order of how much they cost

**1. There is essentially one true oracle in the entire suite.** `dual_contouring/solve/tests.rs:43` checks the QEF solve against an independently derived closed form for a perfect orthogonal corner. Everything else is either mesh-versus-itself (manifoldness, orientation, χ, determinism, self-intersection) or mesh-versus-the-same-field-that-generated-it (`accuracy.rs` Hausdorff against `|f(p)|`, `sealing.rs` against the field's sign).

The second class has a blind spot with a name: **a bug shared between a field's `sample`/`gradient` and whatever trusts it is invisible.** That is not hypothetical — **M-289 is exactly that bug**, and it survived because the reference gradient and the thing it checked were the same code's assumptions. The repo earned the rule *"a reference implementation needs the same scrutiny as the thing it checks"* and then did not add an instrument that would enforce it.

**The cheapest genuine second oracle is the divergence theorem.** The signed volume of a closed mesh is `Σ (1/6)·(a · (b × c))` over its triangles. For `sphere`, `torus` and `box_exact` the true volume is a closed form. That is a dozen lines, no dependency, and it catches a *systematic* bias — every vertex offset slightly the same way — which Hausdorff-with-a-tolerance hides by construction. It is independent of the field's gradient, which is precisely the axis M-289 failed on.

**2. `regress.sh` never runs for real.** It exists, it does exact comparison on structural columns and tolerance comparison on timing columns, and it compares against a per-machine-slugged baseline. **CI runs only `--self-test`** — proving the checker still works — because no runner matches `amd-ryzen-9-5900x-12-core`. So regression checking depends entirely on one person running it by hand on one machine. Given **M-280** (the committed sweep was 1.45× stale and nobody knew) and **M-281** (adding an unrelated function moved a timing 152.5 → 130.8 ms), this is the single largest hole in the harness.

**3. 26 of 29 measurement CSVs have no baseline at all.** Only `ablation`, `field_quality` and `shootout` have `docs/measurements/baseline/*-amd-ryzen-9-5900x-12-core.csv` counterparts. The other 26 are regenerated every run with nothing diffing them. They can rot silently, and **M-295 already recorded a count rotting twice.**

**4. `docs/experiments.md` is stale by 16 entries.** Its scorecard stops at P-21; P-22 … P-37 are registered and resolved. `backlog_gate.sh` checks `FINDINGS.md` but not this page, and `doc_facts.sh` does not derive its count. This is precisely the class of drift `doc_facts.sh` was built for, applied everywhere except here.

**5. Adding a new algorithm touches 6–7 files by hand with no single seam.** Walked from the code: implement `Extractor` (extractor.rs:71) → add to `for_each_extractor!` → hand-edit `extract.rs`, `resolution_sweep.rs`, `stage_breakdown.rs` (bench files construct extractors by name, not through the macro — X-001's own header concedes this) → add to `golden.rs` and regenerate hashes with `ISOMESH_BLESS=1`, reading the diff by eye → mirror the six `check_<algo>` functions in `property/extraction.rs` and choose a `SurfaceGate` → decide a CSV and a baseline → add a Bevy example by convention.

Compare against `P-` (one list) and `ReferenceField` (one macro). The extractor path is the least-served registry in a repo that is otherwise good at registries — and it is the path a *novel algorithm* takes.

**6. The harness assumes triangle soup on a uniform grid, and this is already binding.** `MeshSink::vertex(position, normal) -> u32` / `triangle(a,b,c)` (mesh.rs:44-48) has no room for per-vertex UV, material or colour. `Extractor::extract_into(field, shape, origin, cell_size, out)` (extractor.rs:80) has no way to express an adaptive or multi-resolution extraction.

The evidence that this is a real constraint and not a hypothetical: **`Transvoxel` needed a non-`Extractor` path and is absent from the registry, the property suite and the golden sweep entirely.** `GreedyQuads` triangulates before handing off, so there is no native quad path anywhere in validate, golden or property — despite Grosso & Zint's parallel dual marching cubes (in corpus) being a quad-only method, and despite quad output being on the near horizon.

**7. Four registered experiments have no bench file.** P-10, P-18, P-24, P-25. P-10 is marked "never ran"; the others' status is invisible without reading `FINDINGS.md`.

**8. `subgrid_marching_tetrahedra` is registered but skipped by the property suite** — so the algorithm with the most unusual internals (it samples *inside* cells, and **M-307 measured it worst at 18/24 with interior sealing failures**) gets the least random-input coverage of anything shipped.

### Ranked additions — all zero-dependency

| # | Addition | What it catches that nothing catches | Effort |
|---|---|---|---|
| 1 | **Run `regress.sh` for real in CI**, skipping cleanly when no baseline matches the runner | Perf and metric regressions shipping silently. The script already exists — this is a CI-yaml change | XS |
| 2 | **Divergence-theorem volume oracle** against closed-form volumes for `sphere`, `torus`, `box_exact` | Systematic vertex bias; a second *independent* oracle on the axis M-289 failed | S |
| 3 | **Extend `doc_facts.sh` to derive `docs/experiments.md`'s count from `PREREGISTERED`** | The 16-entry drift that already exists | XS |
| 4 | **Baseline-coverage assertion**: every `docs/measurements/*.csv` either has a baseline or is explicitly listed as record-only | 27 CSVs that can rot unnoticed | S |
| 5 | **Mutation check for Surface Nets and Dual Contouring**, mirroring `march_with_table` | A corrupted SN/DC adjacency table currently passes everything; only MC has a demonstrated-failing negative control | M |
| 6 | **Put `subgrid_marching_tetrahedra` into the property sweep**; wire `Transvoxel` and `GreedyQuads` into `for_each_extractor!` even in a degenerate configuration | Three shipped algorithms outside the shared machinery | M |
| 7 | **A second closed-form QEF oracle case** (a planar two-crossing cell has a known least-squares answer) | The only true oracle exercises exactly one geometric configuration | S |
| 8 | **A `register_extractor!` seam** that drives benches, golden and property from one list | The 6–7-file hand-edit that every new algorithm pays | M |

---

## Axis 2 — Usability in the Bevy ecosystem

### Better than I expected, and better than the alternatives

`bevy_isomesh/src/plugin.rs` is a real plugin, not a bag of functions: `IsomeshPlugin: Plugin`, components `VoxelVolume` / `VoxelChunk` / `NeedsRemesh` / `ChunkMesh(Handle<Mesh>)`, resources `MeshBudget` / `MeshStats`, two chained `Update` systems (`spawn_meshing_tasks` → `apply_finished_meshes`) using `AsyncComputeTaskPool` correctly so extraction never runs on a system thread. It depends on **leaf crates** (`bevy_app`, `bevy_ecs`, `bevy_asset`, `bevy_mesh`, `bevy_tasks`) rather than the `bevy` umbrella, and is deliberately silent on `bevy_render` so headless and server users do not compile a renderer. `mesh.rs`'s `MeshBuilder` is a `MeshSink` writing straight into Bevy's own `Vec` layout — zero-copy into `Mesh`.

Against what a Bevy user would otherwise reach for (`fast-surface-nets`, `block-mesh`): generic over `f32` **and** `f64` where both competitors are `f32`-only; seven extractors with a *measured* decision table instead of one; sharp-feature dual contouring, which neither has; the validity and manifold harness shipped as **public API** rather than test-only; and an async, frame-budgeted plugin, where neither competitor ships a Bevy plugin at all. The core README's "In 60 seconds" is five lines and is compiled as a doctest, so it cannot rot.

### The barriers, bluntest first

**1. Rule 1's own promise is unimplemented.** `CLAUDE.md` says *"Offer `From`/`Into` behind optional features, never in the core API."* Grepping the tree for `From<Vec3>` / `Into<Vec3>` under `crates/` and `bevy_isomesh/src` returns **nothing**. Every position, origin and layout argument a Bevy user touches must be destructured and rebuilt by hand against `Vec3` and `Transform`. In the examples this shows up as ad-hoc `[t.translation.x, t.translation.y, t.translation.z]` scattered per-example.

This is the single biggest adoption barrier, and it is a **half-day fix that violates no constraint** — the feature belongs in `bevy_isomesh`, not `crates/isomesh`, so rule 1 and rule 2 both stay intact.

**2. Every consumer rewrites the same two systems.** Nothing generates the chunk grid or attaches `Mesh3d`/`MeshMaterial3d`. The documented Bevy path is: add the plugin → build a `ChunkLayout` → **hand-write a nested `for z/y/x` loop** spawning `VoxelChunk` + `NeedsRemesh` per index → **hand-write an `attach` system**. That is boilerplate a plugin normally absorbs.

**3. Streaming policy is example-only.** `ChunkLayout` gives the coordinate maths (`chunk_of`, `world_of_sample`, `sample_origin`, `at_lod`) and the plugin gives budgeted async meshing — genuinely more than the competitors. But deciding which chunks to load and unload as the player moves lives in `examples/common/mod.rs` (715 lines) and in `game_terrain_stream.rs`, which a real integrator must copy wholesale. `examples/common/` is undocumented, unpublished glue.

**4. No `SystemSet`.** Consumers cannot order their own render or physics work against meshing without guessing.

**5. `isomesh-gpu` has no `[features]` section at all**, so a consumer compiles every shader path whether or not they use it.

**6. Example balance.** 36 examples, roughly 25 demonstrating an algorithm or property against 11 showing a user how to build something. Reasonable for a research crate; thin for "dig a hole in terrain and render it." (The counts in `README.md` and `DEMOS.md` are **correct** at 36 — I expected drift here and did not find any.)

### Ranked additions

| # | Addition | Effort |
|---|---|---|
| 1 | **`glam` / `bevy_math` `From`/`Into`, feature-gated in `bevy_isomesh`** | S |
| 2 | **`spawn_chunk_grid` helper** replacing the hand-written nested loop | S |
| 3 | **A `SystemSet`** for the plugin's two systems | S |
| 4 | **Distance-based load/unload built on `ChunkLayout`**, promoting `game_terrain_stream`'s logic into the plugin | M |
| 5 | **Feature flags for `isomesh-gpu`** (mesh-shader separate from compute-MC) | S |
| 6 | **Extend `doc_facts.sh` over `bevy_isomesh/DEMOS.md`** so its example count stays derived rather than hand-maintained — it is right today | XS |

---

## Axis 3 — Harnesses that make future experimentation iterative

### The note-taking and gating discipline is the crate's real asset

Three mechanisms carry it, and each was paid for by a specific incident:

- **`FINDINGS.md` is append-only and falsified entries are never deleted.** It records which *sources* to distrust, which outlives the individual fact. 7,377 lines.
- **Golden hashes** — 216 committed FNV-1a hashes over every position, normal and index bit pattern, across 7 fields × 3 resolutions × N algorithms, gated bit-for-bit and **verified identical across macOS/arm64 and Linux/x86-64** (M-31). This is the only check that catches a regression that is geometrically and topologically invisible, and it is what makes an optimisation like A-023 provable rather than plausible: *3.09× faster, byte-identical output.*
- **`preflight.sh`** — 13 steps, about three seconds warm, including `doc_facts.sh`, which derives rottable counts mechanically rather than trusting prose. Added after **M-293**, where an example broke and 58 commits passed with every local gate green.

### What is missing for the *next* class of experiment

**1. Nothing measures the cost of a *sequence* of edits.** Every bench measures a single extraction or a single edit. The consumer's real workload is thousands of chunks re-meshed under a stream of edits with a frame budget. `edit_trace` is the closest thing and it is one trace. **M-124's rule — amortised is the wrong statistic for the frame a breakthrough lands on — has no harness that enforces it.** A "worst single frame over a scripted edit sequence" bench would be the first instrument that measures what the game actually experiences, and it is the natural home for candidate 15 in the novelty table.

**2. There is no cross-implementation check.** Nothing compares against VTK, CGAL or PyMCubes on the same input. This is the strongest available oracle class and it costs nothing at runtime — generate the reference once, commit the hash. Worth weighing against the honest objection that a second implementation may share the first's assumptions.

**3. The self-intersection detector rewrote itself once and the consequence was recorded only in a commit message** (✗25). Detectors that change meaning need the same golden-fixture treatment as the algorithms they judge: a committed fixture with a known self-intersection count, gated.

**4. `FieldBound` is good and under-used.** `Exact` / `Lipschitz{l}` / `Underestimate{q}` / `Unbounded` replaced a `bool` that had been silently wrong for months. It currently gates only whether `accuracy.rs` may compare `|sample|` to a true distance. It could equally gate empty-cell rejection, octree traversal certificates, and sphere-tracing step sizes — several novelty-table rows want exactly this information and would otherwise re-derive it.

**5. No harness expresses "the same field at two resolutions must agree at the seam."** Transvoxel exists and is outside the property machinery entirely. Given that transition cells are the answer the corpus actually supports for seams, the LOD-seam property is untested by anything except an example.

---

## Summary — the six things I would do, in order

1. **Answer the four blocked decisions.** The queue is a queue of one until you do. (Memo: `2026-08-18-four-blocked-decisions.md`.)
2. **Take R-027a** — an `S`, unblocked, and it may dissolve R-027's `L` entirely.
3. **Run `regress.sh` for real in CI.** Largest hole, smallest fix.
4. **Add the divergence-theorem volume oracle.** A dozen lines for the repo's second independent ground truth, on the exact axis M-289 failed.
5. **Ship the `glam` conversions in `bevy_isomesh`.** Half a day, removes the biggest adoption barrier, breaks no rule.
6. **Repair the corpus metadata.** 37.3% of search hits are unidentifiable, which makes rule 5 — never guess algorithm details, look it up — unreliable in practice. See the companion audit.
