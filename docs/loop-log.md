# isomesh loop log

One entry per pass of the Wiggum loop (`.omp/commands/build-loop.md`). Append; never rewrite another pass's entry.

## Iteration 0 — 2026-10-30 — loop setup

By the owner's session, not a loop pass. Added `.omp/commands/build-loop.md` (from the bevy_stigmergy loop command, rewritten for this workspace) and this log. No code, test, manifest or CI file changed.

**Published state.** `isomesh`, `isomesh-gpu` and `bevy_isomesh` are all at **0.0.11** on crates.io (published 2026-10-14). The repository has **no git tags**; the version was set at commit `7016263` (`Release 0.0.11`), which is the diff base for the first release pass. `CHANGELOG.md` exists; its `## [Unreleased]` section still holds the 0.0.11 notes, so the first release pass moves those under `## [0.0.11] - 2026-10-14` before adding its own version.

**Host.** big_mac, macOS on Apple silicon, native Metal: `isomesh-gpu`'s device tests and the Bevy examples run here. `CARGO_TARGET_DIR` is set by the supervisor. Scratch is `/tmp/iN/`, baselines `~/iN_base` and `~/iN_base-target`.

**Gate.** Section 7 of the command: `preflight.sh --full` plus what CI adds (workspace lib/tests and doc tests, the `experimental` feature, `regress.sh --self-test`, the resolved-graph bevy count, `bevy_isomesh`'s MSRV check, `publish.sh --dry-run`, `build_web.sh`). The site step needs the `~/.venvs/isomesh` venv; the command says how to make it.

**Found while setting up, for the owner (listed, not touched):** `0001-docs-research-2026-08-18-sweep.patch` (98 KB) and `_to_delete/git-index.lock` are tracked at the repository root. `.gitignore` ignores `.omp/*`, so the command file was added with `git add -f`; it is tracked from here on.

### Ranked backlog

1. **#10 — an empty chunk's mesh is extracted to the render world** (Sanctuary: 102 `contains no vertices` warnings per world build, 102 of 856 chunks). A documented-behaviour bug filed by a consumer, so it ranks first. The issue's own proposal: `MeshBuilder::into_mesh` (`bevy_isomesh/src/mesh.rs`) gives a zero-vertex mesh `RenderAssetUsages::MAIN_WORLD` only, and every finished chunk keeps its `ChunkMesh`. Measurable by its acceptance 1-4: no published `ChunkMesh` has zero vertices and `RENDER_WORLD`; the existing "chunk never finished" test passes unchanged; an empty mesh's `count_vertices()` is 0 and `try_attributes()` is `Ok`; a chunk that gains or loses its surface after an edit flips usages and drops the old asset.
2. **B-020 is unblocked: re-enable the avian demos.** `BACKLOG.md` blocks it on "an `avian3d` release for Bevy 0.20"; crates.io now has `avian3d` 0.8.0 depending on `bevy ^0.20` (checked 2026-10-30). Done when `game_destruction` and `game_capsule_walk` build under `cargo check --all-targets` in `bevy_isomesh` with no feature flag, the `avian-demos` feature and the two `[[example]]` entries are gone, and the row moves to `BACKLOG_ARCHIVE.md`.
3. **A storable edit list for a carved volume** (what bevy_emerge's F5, #101, would otherwise hand-roll: its Volume tool stores a volume "as its edit list (primitives and operations), not as voxels", with sphere, box and capsule brushes, smooth union and subtraction, and one undo step per carve). Nothing in `crates/` or `bevy_isomesh/src` derives `Serialize` today. Build the general form in `bevy_isomesh` (the core stays at `libm` only, hard rule 3) behind an optional `serde` feature. Measurable: an edit list round-trips through RON bit for bit, and the volume rebuilt from it meshes to the same hash as the original; a carve across a chunk boundary is seam-free by `ChunkSeams`.
4. **Measure Metal for the mesh-shader path.** `CLAUDE.md` says "Metal itself is still unmeasured" and asks for `cargo run -p isomesh-gpu --example mesh_shader_probe` on the Mac. This host is that Mac. Done when its output is recorded as an **M** entry in `FINDINGS.md` and `CLAUDE.md`'s sentence is corrected to the measurement.
5. **A regression baseline for this machine.** `docs/measurements/baseline/` holds only `amd-ryzen-9-5900x-12-core` files (`ablation`, `field_quality`, `shootout`), so `./scripts/regress.sh shootout` on big_mac finds no baseline and gates nothing. Record `shootout-<this machine's slug>.csv` from a clean `main` (an additive `assumed owner exception: <machine> baseline` commit), so every later performance claim here has something to be compared against.
6. **The topmost unblocked `BACKLOG.md` ticket**, R-238 (Phase 31), when items 1-5 are done or blocked.

## Iteration 2 — 2026-10-30 — an empty chunk's mesh stays out of the render world (#10)

Numbered 2 because pass 1 (the driver's) worked #10 for four hours and was stopped with no log entry: it committed the corpus and the research note to `loop/iter-1` and left two files (`plugin/tests.rs`, `CHANGELOG.md`) that the driver parked as `5c47ca7 WIP(pass 1, parked)`. Its checkpoints in `refs/wiggum/pass-1-*` show it never got past writing the tests: no implementation, nothing built. This pass adopted it — cherry-picked `464121c` and `cb77d33`, took the two parked files with `git checkout 5c47ca7 --` — and finished it on `loop/iter-2`. `loop/iter-1` (local and on origin) is left as it is.

**Why #10.** Triage: five open issues, all bugs (#10, #11, #12, #13, #14). #10 is the oldest, a consumer-reported wrong behaviour (102 allocator warnings per world build in Sanctuary, 102 of 856 chunks) and the parked tree was for it. Runner-up: #11 (Subgrid chunks never re-marked after an edit, a stale-mesh bug), which is newer and had no started work. Skipped this pass, in order, with why: #11, #13, #12, #14 — one item per pass; all four are workable and rank first for the next.

**What changed.** `MeshBuilder::into_mesh` and `to_bevy_mesh` build a mesh with no vertices as `RenderAssetUsages::MAIN_WORLD` (one private `usages_for`); every other mesh keeps the default. Rustdoc on `into_mesh` gives the why and the one way to defeat it (filling an empty mesh in place with `insert_attribute`); `to_bevy_mesh`, `ChunkMesh` and the `bevy_isomesh` README contract point to it. The README's Systems line claimed "no public `SystemSet` yet"; `IsomeshSystems` has been public since 0.0.11, so that line is corrected beside it. `Assumed decision` (pass 1's, on #10, kept): the rule lives in both constructors, not the plugin only.

**Found while verifying, not in pass 1's note:** a drawn mesh overwritten in place by an empty one (the re-mesh pattern of `tape_pruning`, `aperture_gate`, `game_edit_tape_trim`) drops `RENDER_WORLD`, and Bevy never evicts its `RenderMesh` on modification. It still stops drawing, because `MeshAllocator::free_meshes` frees every modified id (V-52, from `bevy_render-0.20.0` source and confirmed on Metal).

**Sources.** Algorithm: `transvoxel_dissertation_lengyel2010` p. 28 (quote re-checked verbatim in distill this pass). Shipping game: the same dissertation, ch. 7 p. 99 (C4 Engine). Technical: `book_real_time_rendering_4th_akenine` pp. 817-818, plus Bevy 0.20 source (`render_asset.rs`, `mesh/allocator.rs`, `bevy_pbr` `render/mesh.rs`).

**Home-still, live.** distill MCP worked (`distill_search`, `catalog_recent`); the `hs` CLI's distill search answered `401 missing bearer token`; `system_status` MCP was not connected; `hs status`: scribe on bmb up (3 converting), scribe on 192.168.1.110 backend unavailable, distill idle (CUDA), qdrant up. Corpus: five documents from `catalog_recent` added (`c0b4040`).

**The host.** Every pass process runs at nice 10 / priority 4 (background QoS) beside about a dozen other loops' rustc processes: load ~22, swap 5.9 of 6 GB, and a CPU-bound probe got 4.1 s of CPU in 39.7 s of wall (~10% of a core). This pass's first release build of two Bevy examples hit a 60-minute timeout, and the resumed build took another 85 minutes. Window screenshots (`ISOMESH_SCREENSHOT`) come back all (0,0,0) from a pass on this host, before and after alike, and `screencapture` says "could not create image from display"; the visual evidence below renders headless to an `Image` target instead.

**Evidence** (scratch under `/tmp/i2`, removed at the end of the pass; the numbers are in M-496):

| Check | Before (`41e64d3`) | After |
|---|---|---|
| New tests (2 plugin, 2 mesh) | 4 FAILED: empty mesh is `MAIN_WORLD \| RENDER_WORLD` | pass; `bevy_isomesh` lib 42/42 (release) |
| `an_edit_during_extraction_is_requeued_rather_than_swallowed` (acceptance 2) | ok | ok |
| Smoke `contains no vertices` (8 empty plugin chunks + 1 direct + 1 in-place swap) | 10 | 0 |
| `game_budget` warnings, 120 frames, 2 ms budget | 836 | 0 |
| Empty `ChunkMesh` usages | 8 of 8 `MAIN_WORLD \| RENDER_WORLD` | 8 of 8 `MAIN_WORLD` |
| Smoke frames A (both spheres) and B (in-place sphere emptied) | reference | pixel-identical, max channel diff 0; B's right half uniform 43/44/47 in both |
| Frame C: refill in place, sphere moved into the 4 empty chunks (after only) | — | both spheres drawn, no ghost (old region uniform 43/44/47), no crack at 4x zoom, 0 warnings, 4 empty chunks all `MAIN_WORLD` |
| `game_budget` `budget,` mean ms/frame at 40 / 80 frames | 4.423 / 4.974 | 2.543 / 3.504 |

The `game_budget` figure is noise on this host, not a speed-up: the code it times (extraction) is unchanged and asset creation is outside its timed region. The change adds one integer comparison per mesh and no allocation.

**Mesh checklist.** `git diff main -- crates/` is empty, so validity by field, orientation, the link walk, index bounds, seams, NaN, golden hashes and `cargo tree -p isomesh -e normal` are unchanged by construction; the `extract` and `shootout` benches were not re-run, because the code they compile is byte-identical and timings on this host are not comparable run to run (the row above). Seam and round-trip tests in `bevy_isomesh` pass.

**Critique and verdict.** GREAT. Flaws kept on purpose: filling an empty mesh in place never draws (documented; a plugin-only rule would leave direct callers warning), and Bevy keeps a stale CPU-side `RenderMesh` after an in-place swap to empty until the handle drops (GPU slabs freed, draw skipped). Not measured: Sanctuary's own 102 → 0 (its acceptance 5 is a host check). Found and filed: #16, the README and quickstart attach `Mesh3d` with `Without<Mesh3d>`, so a re-meshed chunk keeps drawing its first mesh (read from source; not yet run).

**Gate.** Local: fmt (bevy), backlog gate, findings index, doc facts, version gate, csv provenance, readme sync, toolchain drift, the site's prose and link check, `bevy_isomesh` lib tests (release). The rest ran as CI on PR #15 on `big`: lint (fmt, clippy, docs, dependency gates), `bevy_isomesh` workspace, MSRV, `test (big)`, package dry run all pass; `pages site` failed with rustc SIGKILL while linking `game_dig` for wasm, the same OOM kill (exit 137) `main`'s own CI run at `41e64d3` had in that job before this change.

**Commits.** `fc56346`, `cd310e3` (pass 1's research, cherry-picked), `c0b4040` (corpus), `cde82ef` (implementation and tests), `a95c2b6` (docs, M-496, V-52), `14261ac` (M-496, frame C). No sibling-repo commits.

### Ranked backlog

1. **#11 — Subgrid chunks are never re-marked after an edit** that adds or removes sub-cell geometry: `mark_edit` decides "surface cell" from corner signs, which Subgrid does not use. A wrong (stale) mesh, so it ranks first.
2. **#13 — `EditReport::dirty_chunks` counts the caller's whole backlog**, not this edit's chunks, and scans the queue quadratically: a wrong number and a per-edit cost.
3. **#16 — the README and quickstart attach pattern never refreshes `Mesh3d` after a remesh.** Confirm with a test or a `Without<Mesh3d>` smoke first, then attach on `Changed<ChunkMesh>` everywhere the pattern is copied.
4. **#12 — `an_edit_during_extraction_is_requeued_rather_than_swallowed` passes when the edit is swallowed**: make it count extractions.
5. **#14 — `CLAUDE.md` hard rule 3 names `glam`**; everything else says `libm`. A one-line fix.

Then iteration 0's items 2-5 (B-020 avian demos, a storable edit list for bevy_emerge F5, the Metal mesh-shader probe, this machine's regression baseline). **Release:** this pass changed `bevy_isomesh`, so 0.0.12 is due; see the release line below.

## Iteration 3 — 2026-10-31 — Subgrid volumes mark the chunks their thin features changed (#11)

**Why #11.** Triage: five open issues, all bugs (#11, #12, #13, #14, #16; #16 was filed by iteration 2). #11 is the oldest and reproduces: traced from `dirty.rs:412-432` and then shown by two failing tests. It is a stale-mesh bug, the first class in the ranking: a thin wall cut from a `Subgrid` volume stays on screen and in the collider, one added never appears, and `EditReport::output_changed_cells` is 0 with no error. Runner-up: #16 (the README's `Without<Mesh3d>` attach never refreshes a re-meshed chunk), newer and not yet confirmed by a run. Skipped, workable, for the next pass in this order: #16, #13, #12, #14. No parked work to adopt: `git branch --list 'wiggum/parked-*'` is empty, and the one `parked` commit, `5c47ca7` from pass 1, was consumed by iteration 2.

**What changed.**

- `isomesh::chunk::dirty::mark_edit_subgrid(layout, before, after, min, max, samples, dirty)`. `mark_edit` decides "surface cell" from the signs of a cell's eight corners, which is all Marching Cubes reads and less than Subgrid does: a plate thinner than a cell leaves every corner outside and still puts two roots on an edge. The new function compares the fields at the points the extractor reads: `samples` cuts along each of a cell's 19 tetrahedron edges (the 12 cube edges, 6 face diagonals, the main diagonal, read off `TETS` rather than written out), each lattice point evaluated once, one layer of cells in memory at a time. `samples == 1` is `mark_edit`'s rule counter for counter. The shared counters moved into one private `tally`.
- `bevy_isomesh::mark_changed_chunks` picks the rule from `volume.extractor` with a `match` that has no wildcard, so a new extractor must say what it reads. `Marching Cubes`, `Surface Nets` and `Dual Contouring` keep `mark_edit`.
- `bevy_isomesh::mark_changed_chunks_in(…, region)`: the same over a cell box the caller passes, because the new rule costs `1 + 7 (samples - 1)` evaluations per cell and `mark_changed_chunks` compares the whole volume's box. `mark_changed_chunks` is now this over the spawned chunks' box; one path. Assumed decision (owner may override), also in the commit body: the region is an inclusive `(min_cell, max_cell)` pair of global cell indices, as `mark_edit` takes it, over a world-space box (rounds in the crate and hides the widen-by-a-cell rule) or a region struct (a new type for two arrays).
- `subgrid/roots.rs`: `edge_point`, `edge_step`, `sample_parameter` shared by `all_roots` and the marker (same expressions as before).
- `spawning_the_work_does_not_do_the_work` rewritten from a wall-clock ratio to a structural assertion (see Gate).

**Sources.** Algorithm: `sig2026_Subgrid_Marching_Tetrahedra` §1, §1.3 (what the extractor reads, and that 1D sampling "can of course miss intersections"). Shipping game: `transvoxel_dissertation_lengyel2010` p. 28 (C4 regenerates every block in an edit's box, no filter) and p. 99. Technical: `10.1145_3730845`, *Sphere Carving* §1 (emptiness can be proven only from a conservative distance bound, which `Sdf` does not promise), plus this repo's `subgrid/roots.rs`, `extract.rs` and M-32. Note: `docs/research/subgrid-dirty-marking.md`. Corpus: six stems added from `catalog_recent` (`6772092`).

**Home-still, live.** `catalog_recent` and `markdown_read` worked. `distill_search` answered `500 gpu busy` (bge-m3 needs 5000 MB, 3.3 GB free; holders `VLLM::EngineCore` 16 GB and others) on all four tries, so the passages were read with `markdown_read` of the same stems and no new document was searched for; listed under the note's `Pending index` with one open question (no source on incremental dirty-region tracking for sparse voxel edits was found beyond C4's whole-box rule). `system_status`: two scribes healthy (one converting 4), distill healthy on CUDA, 852 conversions queued.

**The host.** Load 10-12 on 10 cores at nice 10. A debug build of `isomesh` tests took 10 min, a bench-profile build 16-17 min, `bevy_isomesh` clippy 19 min and rustdoc 22 min; the debug `isomesh` lib suite did not finish in 50 min (the Subgrid experiment tests), so the full gate ran on CI as PR #17, and local runs were the targeted ones below. A wall-clock bench round disagreed with itself by 44% (below).

**Evidence.**

| Check | Before (`a0007b4`) | After |
|---|---|---|
| Issue's plate, plugin test, both directions | would fail (mutated to the corner rule it fails: "chunk 2's Subgrid mesh changed and the chunk was not marked") | pass |
| Oracle: chunks whose Subgrid mesh changed over 8 edits, extracted before and after and compared as buffers | 10 | 10 |
| …of those marked by `mark_edit` (corner rule) | 6 (missed: plate added, plate removed, crack cut, crack closed) | n/a |
| …of those marked by `mark_edit_subgrid` | n/a | 10, at 11 marked (1 re-meshed to the same bytes, which the corner rule also marks) |
| A crack in solid: `value_changed_cells` under the corner rule | 0, mesh changes | `mark_edit_subgrid` marks chunk 0 |
| `samples == 1` vs `mark_edit` | n/a | reports and marked sets equal on all 8 edits |
| Evaluations per cell per field (counting field, samples 1, 2, 4, 5) | 1 | `1 + 7 (samples - 1)`, exact count of distinct points |
| `bevy_isomesh` lib tests, release, local | 42 | 46 of 46 |
| Mutation: `spawning…` test with a 200 ms sleep after each spawn | n/a | fails `(12, 0)` vs `(0, 12)` |

**Mesh checklist.** `git diff main -- crates/` is `chunk/dirty.rs`, the new `chunk/dirty/subgrid.rs` and its tests, and the three helpers in `subgrid/roots.rs`. Validity by field, orientation, the link walk, index bounds, seams, NaN and the 378 golden hashes: no extraction expression changed; `test (big)` (Linux) and `test (macos-latest)` pass `cargo test --workspace --lib --tests` on PR #17 in every run. Refusals (no panic on caller data): `samples == 0`, an inverted region, a 1001³-cell region at 4 samples and an `i64`-extreme region all return an error or the empty report. Footprint: `cargo tree -p isomesh -e normal` is `isomesh` and `libm`; `bevy` count in the root graph 0.

**Bench.** Wall time is not usable here. One round of `extract --bench subgrid_marching_tetrahedra` on the same code: sphere/33 1.18 s before and 0.66 s after, torus/33 1.61 vs 1.33 s, box_exact/33 2.51 vs 1.56 s; I stopped the other rounds. Instructions retired (`/usr/bin/time -l`, `extract --test <row>`, two runs per side) do not depend on load: sphere/33 26.55e9 and 26.46e9 before, 26.45e9 and 26.46e9 after; torus/33 26.74e9, 26.70e9 vs 26.67e9, 26.64e9; box_exact/33 26.83e9, 26.82e9 vs 26.77e9, 26.80e9. After/before of the two-run means 0.9982, 0.9977, 0.9986, with a 0.3% spread on one side: no regression, none resolvable. `regress.sh shootout` was not run: no baseline for this machine exists (`docs/measurements/baseline/` holds the Ryzen's), and the shootout code is unchanged.

**Critique and verdict.** Written before the verdict in `/tmp/i3/critique.md` (removed with the scratch): the oracle is eight hand-picked edits, not a randomised sweep; the cost is a count and no bench times `mark_edit_subgrid`; marking allocates per call; Dual Contouring's gradient reads and Subgrid's bisection and end-probe ulp are residuals, documented and unmeasured; no example shows it, so no still. **GREAT** for #11. The first draft was NOT YET for one reason, that a Subgrid edit's cost grew 22x with no way to bound it through the plugin; `mark_changed_chunks_in` moved it. `FINDINGS.md` M-497 and a Part 5 rule (a filter must be stated against what the consumer reads).

**Gate.** PR #17 on `big` and `macos-latest`. Lint (fmt, clippy `--workspace --all-targets -D warnings`, rustdoc), MSRV, package dry run, `test (big)`, `test (macos-latest)` and the pages site passed on every run they ran in (the pages site did not die with SIGKILL this time: it passed all four runs of the workflow). `bevy_isomesh (separate workspace)` failed three times on `spawning_the_work_does_not_do_the_work` (frame 11.6-21.9 ms against 25-31 ms for one extraction, needing a 4x margin), in a test whose systems this change does not touch; whether my two new plugin tests, which add CPU to the same process, tipped it is not established. The old assertion was a wall-clock ratio its own comment calls a hostage to scheduling; it is now structural (an update that queues twelve slow extractions returns with all twelve running and none meshed). It then failed once on `clippy::type_complexity` in my refactor, fixed with two type aliases, and passed `cargo clippy --all-targets -D warnings` and `cargo doc` with `-D warnings` locally. Final run, `f8f24c0` (run 41178282235): every job green after one rerun of `test (big)`, which a runner terminated with exit code 143 (SIGTERM) in the `experimental feature` step after every test before it had passed; the rerun passed, so that is a killed job and not a failed test.

**Release.** None: the last release, `a0007b4` (v0.0.12), was about an hour before this pass started and the rule is one every 6 hours; there are no tags. The changes ride in the next release (the `CHANGELOG.md` entries are under `[Unreleased]`).

**Commits.** `6772092` (corpus), `43b98a7` (research note), `ba8ad9c` (implementation and tests), `27304a8` (docs, M-497, Part 5 rule), `12e6668` (the timing test), `083f117` (`mark_changed_chunks_in`), `f8f24c0` (type aliases), then the log. No sibling-repo commits.

**Cleanup.** Removed: worktree `~/i3_base` (156 MB), `~/i3_base-target` (925 MB), `/tmp/i3`. `git worktree list` showed no other worktree and no `~/i*_base*` directory, so there are no leftovers from earlier passes to list. `df -h .`: 44 GB free at the start, 38 GB before cleanup (the shared target dir grew with this pass's builds and is not removed).

### Ranked backlog

1. **#16 — the README and quickstart attach `Mesh3d` with `Without<Mesh3d>`, so a re-meshed chunk keeps drawing its first mesh.** Confirm with a test (a re-meshed chunk's `ChunkMesh` change must be visible through `Changed<ChunkMesh>`), then attach on `Changed<ChunkMesh>` in the README and the five examples that copy it. Every game that copies the README into a world it edits hits it.
2. **#13 — `EditReport::dirty_chunks` counts the caller's backlog and scans it quadratically in `mark_edit`.** `mark_edit_subgrid` already counts distinct owners it inserted; `mark_edit` should do the same, in the same `tally` neighbourhood, and the issue's test goes in `dirty/tests.rs`.
3. **A randomised oracle and a bench for marking.** Random sphere, plate and crack edits against the extract-before-and-after oracle (soundness and over-mark rate on more than eight edits), and a committed bench of `mark_edit` against `mark_edit_subgrid` at samples 1, 2, 4 so the cost has a time and not only a count. Also a caller-provided scratch so a drag that edits every frame allocates nothing.
4. **Is the corner rule exact for Dual Contouring and Surface Nets?** DC reads `Sdf::gradient` at each crossing (`hermite.rs`), which is not a corner. The oracle already exists; point it at `DualContouring`.
5. **#12 and #14.** #12: make `an_edit_during_extraction_is_requeued_rather_than_swallowed` count extractions. #14: one line in `CLAUDE.md` (rule 3 says `glam`; it is `libm`) and the `doc_facts.sh` check the issue gives.

Then iteration 0's items: B-020 avian demos, a storable edit list for bevy_emerge's F5, the Metal mesh-shader probe, this machine's regression baseline. **Release:** due after 2026-10-31 03:01 -05:00 (6 h after v0.0.12).
