# Changelog

All notable changes to the isomesh workspace. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/); versions are the shared workspace version
that `isomesh` and `isomesh-gpu` release under. Everything here is pre-1.0: **0.0.x versions make no
stability promise**, and any release may break the API. Releases are CI-driven — a manifest version
bump landing on `main` is the release (`scripts/publish.sh`, version-driven).

## [Unreleased]

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
