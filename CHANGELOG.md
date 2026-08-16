# Changelog

All notable changes to the isomesh workspace. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/); versions are the shared workspace version
that `isomesh` and `isomesh-gpu` release under. Everything here is pre-1.0: **0.0.x versions make no
stability promise**, and any release may break the API. Releases are CI-driven — a manifest version
bump landing on `main` is the release (`scripts/publish.sh`, version-driven).

## [Unreleased]

The dual path got **4.26× faster with byte-identical output**, a published manifoldness claim was
falsified, and a bug in a reference implementation was found to have recorded two true hypotheses as
false. Exact geometric predicates arrived.

### Added

- **`predicates::orient2d` — a 2D orientation test whose sign is never wrong.** Shewchuk's adaptive
  method ([`10.1007/pl00009321`](https://doi.org/10.1007/pl00009321)): a floating-point estimate,
  returned only where a proven error bound shows its sign cannot be wrong, over an exact expansion
  otherwise. `no_std`, no allocation, no new dependency, and no panic path. Both branches answer the
  same question — this is one algorithm with an early exit, not a fast path with a fallback, and the
  module documents the distinction (T-024a).

- **`predicates::incircle` — whether a point is inside the circle through three others, exactly.**
  Same two-stage shape: filtered estimate over an exact determinant. Its stage-A bound is Shewchuk's
  Table 5, `(10ε + 96ε²)` — deliberately not `orient2d`'s, since the incircle determinant is 3×3 with
  squared entries. The exact path forms **no coordinate differences**, because differences round; it
  expands the lifted 4×4 by cofactors instead. The sign convention is stated for counterclockwise
  `a, b, c` and inverts with their winding, which the docs say and a test pins (T-024b).

- **`Real::UNIT_ROUNDOFF` and `Real::SPLITTER`.** The predicates' error bounds are written in terms of
  the unit roundoff `2⁻ᵖ`, which is **half** `Real::EPSILON` — the latter is the 1.0-to-next gap.
  Additive on a sealed trait, so no downstream implementation can break (T-024a).

- **The naive orientation determinant fails by reporting *collinear*, not by reporting a false
  crossing.** Exactly collinear input cannot break it: `fl(x·y)` depends only on the real product, so
  two equal real products round identically and the difference is exactly zero. The reachable defect
  is a **false zero** — a fixture with exact determinant `1` for which the naive form returns `0.0` —
  which is the worse mode, because "collinear" is the reading a triangulator trusts (M-302).

### Changed

- **`signed_distance_from_mesh_winding` is no longer a second implementation of "distance to the
  nearest triangle".** It carried an unaccelerated scan over every triangle while its own doc comment
  claimed the magnitude was computed *"exactly as `signed_distance_from_mesh` computes it"*. It now
  calls `MeshField`, the same code the pseudonormal path runs, and only substitutes the sign — which
  is the one thing that genuinely differs. **Output is bit-identical**: the old loop took
  `minₜ √(r·r)` and `MeshField` takes `√(minₜ r·r)`, and `√` is monotone and correctly rounded, so
  both select the same triangle and root the same value. Golden hashes are unchanged. The speedup is
  the shape M-260 measured at 3.9× (T-025).

- **`ColliderReadiness` gained `non_manifold_vertices`, and `supports_inside_outside()` now requires
  it to be zero.** `MeshReport` computed it with the link walk; `collider::from_report` forwarded ten
  fields and not that one, so a **bowtie** — two cones sharing an apex — reported zero on every edge
  counter and passed. Two tetrahedra point-reflected through a shared apex give `χ = 3`, and an odd
  `χ` is impossible for a closed orientable surface. **Breaking twice over:** a new public field on a
  struct with public fields, and a mesh that passed the predicate yesterday fails it today (T-021,
  M-300, M-301).

- **Surface Nets, Dual Contouring and Manifold Dual Contouring are 2.5–4.3× faster.** Two changes to
  `DualMesher`, which all three share. `emit_quads` took its loop axis as a runtime value, making
  `p[axis] = a` a dynamically indexed store that kept the coordinate array out of registers and
  created a store-to-load forwarding chain; it is now three `const`-generic monomorphisations. And
  `values` was laid out by the caller's shape, so at 128 samples per axis the plane stride was exactly
  64 KiB — a cache-set aliasing period, worth **3.37×** at the canonical voxel chunk size. The row
  length is now forced odd, unconditionally and idempotently.

  At 256³: Surface Nets **693.8 → 162.7 ms**, IPC **1.20 → 4.09**, `SN/MC` **5.43× → 1.26×**. Below
  32³ Surface Nets is now *faster* than Marching Cubes, which falsifies ✗14 — this repository's own
  earlier claim. **The golden hashes are unchanged**: an optimisation that changes the mesh is a bug in
  the optimisation. (A-023/M-285, A-024/M-287.)
- **`bevy_isomesh` re-exports `ChunkSeams`.** Its README named `ChunkSeams::Gapped` in prose while the
  type was reachable only through `bevy_isomesh::plugin`.

### Documentation

- **`docs/research/2026-08-16-decomposition-preconditions.md`** — what four published convex
  decomposition methods actually require of their input, from their own papers. Two of the four
  require nothing: V-HACD voxelizes, and CPD *"handles non-manifold, non-watertight meshes directly
  without preprocessing."* Input cleanliness is a quality axis, not a gate (R-011, M-300, ✗20).

- **`FINDINGS.md` M-297, M-298, M-299.** No published convex decomposition runs at interactive rates,
  which corroborates M-116's 241–272 ms per fragment from the literature side; a fan's
  `inconsistently_oriented_edges` count is exactly its flip-state changes, so it detects a fan that
  *reverses* and is blind to a star polygon that folds without reversing; and the on-demand/batch
  split is free for the pseudonormal backend and costs a factor of `N` for the winding one.

- **`isomesh-gpu`'s docs.rs landing page said the shaders were still to be written.** They shipped;
  the module docs now list what is exported, and state the finding the crate had been burying — with a
  readback the GPU path is *slower* than the CPU at every resolution measured, and it pays off only
  when you render from GPU memory and never read back.
- **A demo page**, [`bevy_isomesh/DEMOS.md`](bevy_isomesh/DEMOS.md): all 34 examples, an animated
  capture of each, the command, the controls, and the finding it demonstrates.
- **An experiments page**, [`docs/experiments.md`](docs/experiments.md): seventeen predictions
  registered before their measurement, ten of them compiler-enforced, with verdicts. Five of the ten
  held and four were falsified; all are listed.
- **Every GIF re-recorded** at the current commit, and fourteen added for examples that had carried
  capture-driven animation code which had never been run.

### Added

- **`construct::from_mesh::MeshField<'a, R>`** — a mesh's signed distance field, evaluated **on
  demand** rather than sampled onto a grid. Implements `Sdf`, builds the angle-weighted pseudonormals
  and blocked bounding boxes once, and borrows the mesh. `signed_distance_from_mesh` is now this type
  in a loop, so there is one implementation of the query and
  `the_grid_path_is_this_field_in_a_loop` pins the two to **identical bits** across all 35,937 samples
  of the round trip. Deliberately pseudonormal-only: `winding_numbers` casts one ray per grid row and
  shares it across the row, so an on-demand winding twin would cast `N³` rays where the batch casts
  `N²` (S-008, M-299).

- **`scripts/doc_facts.sh`** — derives the counts that keep rotting from their source and fails when a
  document states a different one. The golden-hash count had by then been wrong in the root README
  twice; 0.0.4 below records fixing the same number from 147 to 168, and it was 216. Wired into
  `scripts/preflight.sh`.
- **`scripts/record_all_gifs.sh`** — the per-example capture parameters, which existed only as shell
  history.

### Fixed

- **`scripts/readme_sync.sh` compared every Rust fence rather than the first**, so adding a second
  snippet to either README failed the gate with a diff about the new snippet instead of about drift in
  the quickstart.

## [0.0.5] — 2026-08-15

MC33's interior ambiguity, corrected where it was wrong and refused where the published
construction stops. Everything here is measured; see `FINDINGS.md` M-229…M-233.

### Added

- **`marching_cubes::trilinear::BodySaddles::same_asymptote_side`** — Grosso 2016's Proposition 1,
  derived from the paper's own normal form for the face bilinear rather than transcribed. `FINDINGS`
  V-31 had twice recorded that predicate as prose too vague to pin down; it is not. Measured against
  the contour count over 400,000 random cells it agrees exactly — false precisely when the cell has
  one contour, which is Corollary 1 (M-230).
- **`marching_cubes::trilinear::singular_face_mask`** — detects Grosso 2017 §3's singular face, where
  the bilinear saddle sits on the level set and the asymptotic decider has no answer to give.
  Detection only; the resolution is still owed, and `ambiguity` is deliberately untouched (A-002i).
- **`marching_cubes::trilinear::MAX_TUNNEL_CONTOUR`** — Corollary 6's bound of six, used
  contrapositively as the tunnel test.

### Fixed

- **A case-13 cell with contours of nine and three was meshed as a tunnel and left a hole** (M-228,
  M-229, M-230). The ring count admitted a cell Corollary 6 excludes; sent to the tunnel rule it met
  a contour edge whose endpoints land three steps apart on the inner hexagon, which the construction
  has no rule for and which the authors' implementation silently emits nothing at. Such a cell is now
  classified `Topology::SeparateDisks` and refused with `Error::UnresolvedSixSaddle` **before any
  vertex is emitted**, rather than the gap being discovered partway through the tunnel rule and
  leaving the caller's buffer half-filled.

### Changed

- **BREAKING: `marching_cubes::trilinear::Topology` gains a `SeparateDisks` variant.** The enum is
  not `#[non_exhaustive]`, so an exhaustive `match` on it needs a new arm. `Error` *is*
  `#[non_exhaustive]`, so `Error::UnresolvedSixSaddle` is additive.
- `Error::UnresolvedTunnel` is now a live guard on a case nothing classified as a tunnel has reached,
  rather than a case with a known trigger. Kept rather than deleted because that is a measurement
  over a sample, not a proof.

### Documentation

- **`ManifoldDualContouring`'s manifoldness guarantee now states its precondition** (A-017). It was
  described as the entry that takes the non-manifold count to zero; that held because none of the
  seven original reference fields can produce a cell with an interior ambiguity. On one that can it
  returns 30 non-manifold edges at 17³ and 64 at 33³. This is outside Schaefer, Ju and Warren's
  guarantee rather than a defect against it — they separate sheets *within* a cell and never claimed
  two crossed edges of one **shared** face resolving to the same cycle pair — so it is documented,
  not chased, and the counts stay pinned as whole censuses. No code change.
- `docs/research/2026-08-10-meshing-library-target.md`: Nielson 2003 and Lopes & Brodlie 2003
  corrected from `PAYWALL` to `HAVE` (V-32), and a callout added above the status legend — that code
  records what the resolver could not reach, has now been wrong four times, and correcting
  individual rows did not stop the next pair.
- Half of Grosso's Corollary 6 is recorded as falsified by measurement: it states a tunnel's contours
  are "at most 6 vertices and the other 3", and `[4,4]` tunnels are common (M-230).

## [0.0.4] — 2026-08-14

> **Folded up from `[Unreleased]` on 2026-08-15.** These items shipped in the 0.0.4 wave —
> `bevy_isomesh` 0.0.4 is live on crates.io — but nothing moves the heading, because the release
> is driven by a manifest version bump rather than by editing this file. Check the section before
> cutting a release.

### Added

- **`bevy_isomesh` prepared for its first crates.io release**, at 0.0.4 in lockstep with the
  workspace: crates.io metadata and LICENSE files (D-007), version-carrying dependencies, and the
  release wiring — `scripts/publish.sh` now carries it as the third, last leg (D-010). The upload
  executes on the next push to `main` where the version is absent from crates.io, and fails loudly
  until `CARGO_REGISTRY_TOKEN` exists in the `crates-io` environment (FINDINGS M-198).
- The documentation overhaul (Phase 7, D-001…D-009): a doctested 60-second quickstart on the root
  README with a CI gate holding it identical to the crate README's compiled example
  (`scripts/readme_sync.sh`); an "Is this for you?" fit table where the refusals are stated as
  plainly as the fits; a seven-extractor tradeoff table on the crate page; badge rows; a
  troubleshooting section; `bevy_isomesh`'s README rebuilt as its docs.rs front page
  (`#![doc = include_str!("../README.md")]`) with the plugin's full exposed contract, a compat
  table, and LICENSE files that the manifest had promised but the package never carried; the
  `isomesh-gpu` README example wired into doctests; a **Machines** block in `FINDINGS.md`'s header;
  this changelog; issue and PR templates.

### Removed

- `bevy_isomesh`'s `gpu` cargo feature and its optional `bevy_render` dependency — the feature
  gated zero code, and the examples had always reached render types through the `bevy` umbrella
  dev-dependency (D-002).

### Fixed

- Stale claims across the READMEs: Rust 1.85 → 1.89, the 0.0.2-placeholder story replaced by the
  live releases, 147 → 168 golden-hash combinations, ticket and demo counts converted from rotting
  totals into pointers at their gated sources (D-003, D-004).

### Also added in this release

- `isomesh-gpu`'s two demonstration examples shipped in the package: `extract_a_sphere` (the whole
  GPU loop, checked against the CPU) and `gpu_vs_cpu` (the resolution sweep behind the README's
  timing table). Commit `24821e6`.
- The publish pipeline's token check moved into `scripts/publish.sh`, next to the upload it guards —
  a push that publishes nothing no longer needs a registry secret to stay green, and a real release
  still fails loudly if the secret is missing (`958160d`, FINDINGS M-198).

## [0.0.3] — 2026-08-14

### Added

- First `isomesh-gpu` release on crates.io.
- E-116: the Marching Cubes 33 interior-decider analysis example (`marching_cubes_interior`).

### Fixed

- **The license files that were never shipped**: `LICENSE-MIT` and `LICENSE-APACHE` now travel in
  each published package. isomesh 0.0.0 through 0.0.2 declared `MIT OR Apache-2.0` with no license
  file beside it, unnoticed for three releases. Commit `aa46e35`.
- GPU tests share one device instead of opening 67.

## [0.0.2] — 2026-08-13

### Added

- The crates.io page gets the seam-walk GIF and the demo links. Commit `e4df213`.

## [0.0.1] — 2026-08-13

### Added

- `crates/isomesh/README.md`, written so the crates.io page carries more than a one-line
  description; the version bump exists to publish it. Commit `ae238cd`.

## [0.0.0] — 2026-08-12

### Added

- Name reservation on crates.io (I-005), explicitly a placeholder and explicitly burned forever —
  `megamesh` had been taken 48 hours before this project checked it. 82 files, 329 KiB compressed:
  source, benches, golden hashes, proptest regressions, nothing stray. Commit `7c92084`.

[Unreleased]: https://github.com/ladvien/isomesh/compare/main
[0.0.4]: https://github.com/ladvien/isomesh/commit/24821e6
[0.0.3]: https://github.com/ladvien/isomesh/commit/aa46e35
[0.0.2]: https://github.com/ladvien/isomesh/commit/e4df213
[0.0.1]: https://github.com/ladvien/isomesh/commit/ae238cd
[0.0.0]: https://github.com/ladvien/isomesh/commit/7c92084
