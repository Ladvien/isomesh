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
