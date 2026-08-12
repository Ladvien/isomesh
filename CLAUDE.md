# isomesh

Engine-agnostic isosurface extraction in Rust. Field in, triangles out. Must serve **both** a Bevy
voxel game and a CAD tool, which is the constraint that decides almost every design question here.

## Rules
- Use home still and lookup research on whatever you are implementing BEFORE you implement it.  We want SOTA and best practices.  And it is filled with game development research and best practices.
- Do not use hardbreaks when writing markdown.  Markdown should be easy for humans to read when rendered.
- When given git commands, it should be for all edittable repos in this project.
- When creating new features, attempt to use Bevy's plugin pattern as much as possible.  Create separate workspace crates.  Create their own Github repo with idiomatic name.  This is to ensure reusable components are generated during our work. Ensure each respective crate has a warning label of "Vibe Coded" a the top of the README.md. Please refer to [bevy_plugins.md](docs/bevy_plugins.md).  Each separate create must include examples (1-3) demoing the crate.  This is allow inspection of crate behavior without inclusion in this game.  If you discover debugging needs, make recommendations on adding it to this plugin, ever evolving.
- Do not use unwrap() or anything that'd lead to a panic.  Code safe.  Handle errors.
- Leave academic paper references in comments, if a paper was used in writing the code.
- Rember compilation cost time; try to bunch changes and use `cargo check` to spot issues
- Always run the full test suite (including determinism and headless behavioral tests) after modifying gameplay/simulation code, and verify determinism before shipping.
- Do NOT assume design decisions on my behalf. When a design or scope choice is ambiguous (colors, coverage %, approach), stop and ask before implementing. Prefer focused/concrete changes over global post-process filters or over-engineered solutions.
- When investigating whether an issue is fixed, actually inspect the underlying data/code first before offering explanations; do not assume a file is broken or blame viewport/version.
- When items are complete, move them from `BACKLOG.md` to `BACKLOG_ARCHIVE.md`.  DO NOT DELETE `BACKLOG_ARCHIVE.md`.

## Where the work comes from

**`FINDINGS.md` is the epistemic state.** Before acting on anything "known" — a performance figure,
an algorithm property, a claim about what some crate does — look for it there. If it isn't there, it
hasn't been checked, and you should say so rather than assume. Cite the tier (**M** measured here,
**V** verified from a primary source, **R** reported, **F** folklore) whenever a claim justifies a
decision.

**Adding to it is part of the work, not overhead.** Any of these earns an entry, in the same commit:

- A measurement contradicts something written down — in a doc, a ticket, a comment, or this file.
  **The contradiction is the finding.**
- You verify or fail to verify an external claim.
- You discover a property that wasn't predicted (the Euler identity in ✗1 is the model).
- A failure produces a rule that should stop it recurring — add it to Part 5, with the incident.

Falsified entries are never deleted. They record which *sources* to distrust, which is worth more
than the individual fact.

**`BACKLOG.md` is the work queue and the state.** Take the topmost unblocked, unchecked ticket; don't
cherry-pick. Check the box in the same commit that does the work. If you can't finish a ticket, leave
it unchecked, add a `> BLOCKED:` line saying exactly what's in the way, and move to the next unblocked
one. Never half-finish and check the box.

`docs/2026-08-11-implementation-brief.md` has the staged detail — algorithm specifics, exit criteria,
and the math you'd otherwise have to guess at. Read the relevant stage before starting its tickets.

---

The research that justifies this crate lives in `docs/research/`. **Read the relevant one before
implementing an algorithm** — they contain measured figures, corrected folklore, and the specific
reasons certain designs were rejected. Do not re-derive from memory.

| Doc | Read it before |
|---|---|
| `2026-08-11-meshing-crate-architecture.md` | touching `Cargo.toml`, dependencies, or crate layout |
| `2026-08-10-meshing-algorithm-catalog-v2.md` | implementing any extraction algorithm |
| `2026-08-10-adjacent-math-transfer-audit.md` | writing vertex placement / QEF code |
| `2026-08-11-meshing-speed-analysis.md` | any performance claim or optimization |
| `2026-08-11-novel-gameplay-opportunities.md` | deciding what to build next |

---

## Hard version pins — get these wrong and nothing links

Bevy 0.19 pins these exactly. A mismatch does **not** produce a resolution error; Cargo silently
compiles two copies and you get `expected TextureFormat, found a different TextureFormat` much later.

| Crate | Version | Why |
|---|---|---|
| `bevy` | **0.19** | target |
| `wgpu` / `wgpu-types` / `naga` | **29.0.3** | exactly what `bevy_render` 0.19 pins |
| `glam` | **0.32** | exactly what Bevy 0.19 pins |
| `encase` | **0.12** | same |
| edition | **2024** | |

**Never bump `wgpu` or `glam` independently of Bevy.** If Bevy moves, they move together, in one commit.

Platform here is **macOS / arm64 → Metal**. Two consequences:

- Do **not** put `"x11"` in Bevy's feature list. That's Linux. macOS needs no windowing feature flag
  beyond `bevy_winit`.
- **Mesh shaders are unverified on Metal.** wgpu's spec table lists MSL as *planned* while the tracking
  issue says the Metal HAL backend merged. These disagree. Do not build anything on mesh shaders
  without first running a capability probe and reporting what it actually says.

---

## Layout

```
isomesh/
  Cargo.toml            [workspace] members = ["crates/*"], exclude = ["bevy_isomesh"]
  crates/
    isomesh/            core. no_std + alloc, unconditionally. examples/ = headless, write OBJ.
                        Deps: libm. Only libm — see ✗16; glam cannot serve a crate
                        generic over f32 and f64.
    isomesh-gpu/        + wgpu 29.0.3. API takes &wgpu::Device / &Queue / &mut CommandEncoder.
  bevy_isomesh/         EXCLUDED from the root workspace. Own Cargo.lock. ALL Bevy examples here.
  docs/research/        the papers-derived research. Read-only unless asked.
```

`bevy_isomesh` is excluded deliberately. In a shared workspace, Cargo's feature unification leaks:
building the whole workspace gives `glam` the `std`, `serde`, `bytemuck` and `encase` features, so
`cargo test` in the root stops testing what consumers actually get. Excluding it keeps the root lock
pristine. The cost is no `cargo test --workspace` across the boundary — run both.

---

## Hard rules

1. **No math library in a public signature.** `[f32; 3]`, `[u32; 3]`. glam is an internal
   implementation detail. Reason: Bevy 0.19 wants glam 0.32, `parry3d` wants 0.33, `fast-surface-nets`
   wants 0.29 — a consumer using two of those compiles incompatible `Vec3` types. Arrays are the only
   thing that survives. Offer `From`/`Into` behind optional features, never in the core API.
2. **`grep -r "bevy" crates/` must return nothing.** Not in `src`, not in `Cargo.toml`, not in
   `dev-dependencies`. Bevy exists only inside `bevy_isomesh/`.
3. **Nothing goes in `crates/isomesh/Cargo.toml` `[dependencies]` except `glam`** without a written
   justification added to this file in the same commit. The pitch is "as light as possible"; every dep
   is a decision, not a convenience. Justifications live in the section below.
4. **Never claim a performance number without the benchmark that produced it, in the repo.** The
   research docs list several published figures that failed verification. Don't add to the pile.
5. **Never guess algorithm details.** The papers are indexed in home-still and summarized in
   `docs/research/`. If a case table, a sign convention, or an edge ordering is uncertain, say so and
   stop rather than inventing one — wrong case tables produce meshes that look fine and are subtly
   non-manifold.
6. **Output buffers are caller-provided and reusable.** `fn extract(..., out: &mut MeshBuffer)` with a
   `reset()` that resizes without reallocating. Never return a freshly-allocated mesh from a hot path —
   the real workload re-meshes thousands of chunks per edit.

---

## Dependency justifications

Required by rule 3. One entry per dependency in `crates/isomesh/Cargo.toml` `[dependencies]`, added in
the commit that adds the dependency.

**`libm` 0.2 — I-001.** `core` has no `sqrt` or `floor` on stable (they are `core_float_math`, issue
137578) and no `sin`/`cos` at all — verified in rustc 1.96.1's `library/core/src/num/f32.rs`. `Real`
needs all four: `gyroid` is `sin·cos`, `sphere`/`torus` need `sqrt`, `fbm_terrain` needs `floor`.

It is used **unconditionally** rather than behind a `std` feature switch, for two reasons.

- *One path.* A `#[cfg(feature = "std")]` fork would mean two float backends and two sets of results
  for the same input, which is exactly the class of thing that makes a bug take hours to trace.
- *Determinism.* `std`'s `sin`/`cos` are the platform's libm, are not correctly rounded, and differ
  between macOS and Linux. T-007 commits golden hashes; with `std` those hashes would be
  platform-specific and CI would disagree with the dev machine. `libm` is pure Rust and
  bit-reproducible everywhere.

  **Verified at T-007, not merely argued (M-31):** the 63 golden hashes are generated on
  macOS/arm64 and pass unchanged on Linux/x86-64 in CI — every position, normal and index
  bit-for-bit equal across both.

It costs nothing at run time: `libm::sqrtf` compiles to `fsqrt` on aarch64+neon and `sqrtss` on
x86-64+sse2 (verified in `libm-0.2.16/src/math/arch/{aarch64,x86}.rs`). `libm` itself has zero
dependencies and is maintained under `rust-lang/compiler-builtins`.

**`glam` — not declared, and no longer expected to be.** Rule 1 sanctions it as the internal math
library and ✗10 deferred it to A-007's solve, on the grounds that a 3×3 solve is the first thing that
genuinely wants matrix types. The premise held; the conclusion did not.

**glam has no scalar abstraction** (✗16, verified in glam 0.32.1's source): `Mat3` is `f32`, `DMat3` is
`f64`, they live in separate modules, and there is no generic `Mat3<T>` nor any trait spanning them.
This crate is generic over `Real`, which is both. Using glam would mean a bridge trait with two impls
forwarding every operation — more code than the 3×3 adjugate it would wrap, plus a dependency, plus two
float backends inside one solve, which is exactly what the `libm` justification rejects.

So A-007's solve is a six-entry symmetric matrix over `[R; 3]` in `dc/solve.rs`, and **the crate stays
at one dependency**. Revisit only if glam gains a generic scalar or this crate drops `f64`.

---

## Verification — required, not optional

Every extraction algorithm ships with these before it counts as done:

- **Validity gate, chosen by the field — not one blanket rule.** A blanket `V - E + F == 2` is
  unsatisfiable on two of the seven test fields and must not be written that way: `gyroid` is triply
  periodic, so any finite sampling box cuts it and the result *has boundary*, and `fbm_terrain` is a
  heightfield that leaves through the sides by construction. Neither has an Euler characteristic
  derivable a priori, and inventing one violates rule 5. Each field publishes `ReferenceField`
  metadata and the harness reads it:

  ```
  closed fields (closed_in_domain)  → report.is_closed()
  open fields                       → report.is_manifold()
  χ asserted                        → only where expected_euler() is Some
  otherwise                         → record the observed χ in the golden fixture
  ```

  `is_closed()` already folds in the parity check: for any closed orientable surface `χ = 2 - 2g`, so
  χ must be even even when the genus itself is unknown. There must be no `if field == gyroid` anywhere
  in test code.
- **Index bounds + degenerate triangle test.** No index ≥ vertex count; no triangle with two equal
  indices. **Degenerate (near-zero-area) triangles are a recorded metric, not a gate** — Marching Cubes
  genuinely emits slivers whenever a grid corner value sits near zero. That's the algorithm, not a bug.
- **Manifold test.** Report counts, don't just assert — the numbers are interesting. Note the split:
  `non_manifold_edges` counts edges with **≥3** faces and `boundary_edges` counts edges with **exactly
  1**. Lumping them together as "≠ 2" double-counts and makes zero unachievable for any open mesh,
  which includes every individual chunk.
- **Orientation test.** Every interior edge traversed in *opposite* directions by its two faces. A
  transcribed case table with one flipped triangle passes χ, edge-manifoldness *and* vertex-manifoldness
  while being inside out — this is the only check that sees it.
- **Non-manifold vertices need the link walk, not a count.** "Incident faces == incident edges" reports
  a bowtie as clean: two cones sharing an apex have `2k` of each, every edge has exactly 2 faces, and χ
  can come out right. Walk the connected components of the incident-face link.
- **Self-intersection count per 1,000 triangles.** Not a pass/fail; a recorded metric. Dual methods
  are expected to be non-zero, and the measured question is how much the cell clamp reduces it.
- **Determinism.** Same input twice → byte-identical output. Guard against `HashMap` iteration order
  leaking into vertex order.

Add a shader validation test the moment there is a first shader:

```rust
let module = naga::front::wgsl::parse_str(src)?;
naga::valid::Validator::new(ValidationFlags::all(), Capabilities::all()).validate(&module)?;
```

No GPU needed, runs in CI, catches "works on Metal, explodes on DX12."

---

## Commands

```bash
# core — must stay tiny
cargo test -p isomesh
cargo tree -p isomesh -e normal      # expect: isomesh + libm. Nothing else.
                                     # `-e normal` is load-bearing: cargo's default edges are
                                     # normal,build,dev, so the bare command shows proptest's and
                                     # criterion's whole trees once T-005 and T-006 land.
cargo metadata --format-version 1 | grep -c '"name":"bevy' # expect: 0

# benchmarks -- always release, and always via `cargo bench`
cargo bench --bench extract              # criterion, per-algorithm regression timings
cargo bench --bench resolution_sweep     # 16^3..256^3, fits t = a + b*n^3, writes
                                         # docs/measurements/resolution_sweep.csv
# The sweep is a no-op under `cargo test`. It guards on the `--bench` argument
# cargo passes, because `--all-targets` re-selects bench targets even though the
# manifest sets `test = false`, and a debug-build sweep takes minutes and would
# overwrite the committed CSV.

# gpu
cargo test -p isomesh-gpu

# bevy side — separate workspace, run from its directory
cd bevy_isomesh && cargo run --example mc_sphere --release
cd bevy_isomesh && cargo build --examples
```

**Always `--release` for examples.** Debug-build meshing is 20–50× slower and will make you think
something is wrong with the algorithm.

---

## Working style

- Small commits, one concept each. A new algorithm, its tests, and its example are one commit.
- When a research doc and your intuition disagree, the doc wins or you raise it explicitly. Several of
  its claims are corrections of things that "everyone knows" and are wrong.
- Prefer adding an example over adding prose. This project's output is demonstrable behavior.
- If something is blocked or ambiguous, stop and say so. Do not paper over it with a plausible guess.
