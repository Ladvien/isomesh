# Changelog

All notable changes to the isomesh workspace. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/); versions are the shared workspace version
that `isomesh` and `isomesh-gpu` release under. Everything here is pre-1.0: **0.0.x versions make no
stability promise**, and any release may break the API. Releases are CI-driven — a manifest version
bump landing on `main` is the release (`scripts/publish.sh`, version-driven).

## [Unreleased]

### Added

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

## [0.0.4] — 2026-08-14

### Added

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
