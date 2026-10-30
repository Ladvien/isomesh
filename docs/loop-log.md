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
