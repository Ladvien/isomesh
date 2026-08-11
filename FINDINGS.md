# isomesh — findings ledger

**Started:** 2026-08-11 · **Append-only.** Entries are never deleted, only re-tiered when new evidence
arrives — with the old verdict left visible.

This is the project's epistemic state: what we believe, how strongly, and **on what evidence**. It
exists because this project has already been wrong six times in ways that would have silently
propagated into code, and because the research corpus contains several published figures that failed
verification. A belief with no recorded falsification method is not a finding, it's a preference.

## How to use this file

- **Before acting on a "known" fact, look for it here.** If it's not here, it hasn't been checked.
- **When a measurement contradicts something written down, that's an entry** — the contradiction is
  the finding, not an inconvenience.
- **Falsified entries stay.** They're the most valuable rows, because they tell you which *sources* to
  distrust, not just which facts.
- Every entry names how it could be shown wrong. If you can't write that line, you have an opinion.

## Confidence tiers

| Tier | Meaning | Bar |
|---|---|---|
| **M** | **Measured here** | We ran it. Code and numbers in this repo, reproducible by checkout. |
| **V** | **Verified externally** | We read the primary source ourselves. DOI or file attached. |
| **R** | **Reported** | A credible source asserts it; we have not independently checked. |
| **F** | **Folklore** | Widely repeated, no verified primary source found. |
| **✗** | **Falsified** | Tested and found false. Includes what we believed and why. |

**Never cite an R or F claim as justification for a design decision without saying which tier it is.**

---

## Part 1 — Falsified

The most valuable section. Each row includes where the wrong belief came from, because provenance
predicts the next error.

### ✗1 — "Surface Nets produces substantially fewer triangles than Marching Cubes"

**Believed because:** stated in this repo's own implementation brief, and near-universal folklore.
**Falsified by:** A-004 measurement, then derivation, then independent numeric check.
**What's true instead:** on a closed manifold surface meshed on the same grid, both counts are pinned
by Euler's formula:

```
V_sn = V_mc + χ            F_sn = F_mc + 2χ
```

Exactly, at every resolution. `V_mc = C` (one vertex per crossed grid edge), `F_sn = 2C` (two
triangles per crossed edge), `F = 2V − 2χ` on any closed triangulated manifold. The middle step is a
combinatorial identity worth naming on its own: **surface cells = crossed edges + χ.**

| field | χ | crossed edges C | surface cells S | S−C |
|---|---|---|---|---|
| sphere 17³ / 49³ | 2 | 414 / 4206 | 416 / 4208 | 2 / 2 |
| torus 33³ / 49³ | 0 | 1296 / 3216 | 1296 / 3216 | 0 / 0 |
| box 49³ | 2 | 4374 | 4376 | 2 |
| two disjoint spheres 49³ | 4 | 1300 | 1304 | 4 |

**Consequence:** SN's case must rest on quad connectivity and inner-loop cost, not output size. In
M-001 the count columns are a **checksum with a predicted value**, not a result.
**Would be shown wrong by:** any closed-manifold pair on the same grid where the difference ≠ 2χ.
**Legitimately breaks at:** boundary-clipped meshes (**incoming at G-001 chunking** — expect the
assertion to fail there and do not "fix" it), A-013 welding, MC vs MC33 differing in χ.

### ✗2 — "You can have a manifold mesh or an intersection-free one, not both"

**Believed because:** folklore, repeated in several secondary sources.
**Falsified by:** literature review round 1. Manson & Schaefer 2010 achieved both. ODC (2024) measured
Manifold DC at **100% of models self-intersecting** against ODC at **0 of 1500**.
**Consequence:** guaranteed intersection-free extraction is on the table, which is the premise under
A-009 and the runtime-convex-decomposition opportunity.

### ✗3 — "Every interior Surface Nets vertex has four neighbours"

**Believed because:** written into isomesh's own module docs before measuring.
**Falsified by:** A-004. Measured max degree **10** — higher than MC's **9**.
**Method rule earned:** a doc comment the test suite disproves is worse than no doc comment.

### ✗4 — "Dual Contouring is absent from the home-still corpus"

**Believed because:** `distill_search` returned nothing.
**Falsified by:** `catalog_read` on `10.1145_566570.566586` — present, with zero Qdrant chunks.
**Root cause:** **342 documents are readable but invisible to `distill_search`.**
**Method rule earned:** presence in the corpus is decided by `catalog_read` / `catalog_list`, **never**
by search. Any claim of absence made from a search result is unfounded.

### ✗5 — "naga_oil is the shader composition path for Bevy"

**Believed because:** it was, and it's the name everyone knows. Carried into this project's premise.
**Falsified by:** naga_oil's last release is **v0.18, 2025-06-26** (14 months stale), and Bevy 0.19's
own release notes say they will *"port our existing internal shaders to use WESL, and endorse it as
the shader language of choice for Bevy."*
**Consequence:** GPU-002 uses a ~40-line preprocessor; WESL revisited when shader count justifies it.

### ✗6 — "Mesh shaders aren't reachable from inside Bevy"

**Believed because:** an imprecise claim made in this project's own speed analysis.
**Falsified by:** wgpu **v28** shipped `Features::EXPERIMENTAL_MESH_SHADER`; Bevy 0.19 pins wgpu
29.0.3; `WgpuSettings.features` requests it at device creation and
`RenderDevice::wgpu_device() -> &wgpu::Device` gives raw access.
**What's actually true:** wgpu's *maturity* is the blocker, not Bevy. Experimental, redesign issue
open, no browser support, **Metal status contradictory** (see O-5).

### ✗7 — "CBT 5.78 → 0.40 ms is from Dupuy 2020"

**Numbers confirmed. Attribution wrong** — it's the Unity SIGGRAPH 2021 talk.
**Method rule earned:** verify the *source* separately from the *number*. A correct figure with a
wrong citation propagates as an uncheckable claim.

### ✗8 — Velo3D "93M vertices / 31M faces"

**Unsupported.** No primary source located. Removed from the catalog rather than softened.

### ✗9 — "MC's cost inside the volumetric loop was never measured" / "navmesh rebuild cost was never measured"

Both false. Dong 2018 measured meshing at **76.5–89.6%** of the pipeline; van Toll 2012 measured
navmesh rebuild. Both were asserted as gaps in this project's v2 catalog and corrected in round 1.
**Method rule earned:** "nobody measured X" is a claim requiring the same evidence as any other.

### ✗10 — "glam should be the internal math library from day one"

**Believed because:** stated in this repo's architecture doc.
**Refined by:** the agent at I-001. The core took **libm only**; glam is deferred to A-007, where the
3×3 solve is the first thing that actually needs matrix math.
**Why the refinement is better:** the public API is arrays, so glam was never load-bearing; deferring
it means the crate carries zero exposure to glam's ~quarterly breaking releases until it buys
something. Consumers need no conversion impls either — `glam::Vec3::from([f32; 3])` already exists.

---

## Part 2 — Measured here (tier M)

| # | Finding | Evidence |
|---|---|---|
| M-1 | **surface cells = crossed edges + χ** | 4 topologies × 3 resolutions, table in ✗1 |
| M-2 | `V_sn = V_mc + χ`, `F_sn = F_mc + 2χ` | A-004 tests, all four clean fields |
| M-3 | Surface Nets max vertex degree **10**; MC **9** | A-004 |
| M-4 | SN is non-manifold where one cell carries two sheets: **48** non-manifold edges on capped gyroid, **15** on fbm_terrain at 33³ | A-004; pinned as non-zero assertions, not excluded silently |
| M-5 | On `box_exact`, SN's nearest vertex to the corner (1,1,1) is **1.15 cells** away | A-004 — this gap is what E-104 exists to show |
| M-6 | `libm::sqrtf` lowers to hardware `fsqrt` (aarch64+neon) / `sqrtss` (x86-64+sse2) | libm 0.2.16 source: `src/math/arch/aarch64.rs` raw asm, dispatched by `select_implementation!` on `target_feature` |
| M-7 | dev-dependencies do not propagate: consumer resolves **3 packages**, the crate's own lockfile has **137** | Experiment, cloud container |
| M-8 | Cargo silently co-resolves two wgpu majors — **317 packages, both 29.0.4 and 30.0.0**, no resolution error; fails later as `expected TextureFormat, found a different TextureFormat` | Experiment |
| M-9 | Workspace feature unification leaks: `-p isomesh` alone gives glam `libm`; whole-workspace gives it `std`, `serde`, `bytemuck`, `encase`, `rand` | Experiment — the reason `bevy_isomesh` is excluded |

---

## Part 3 — Verified from primary sources (tier V)

| # | Finding | Source |
|---|---|---|
| V-1 | Bevy 0.19 pins **wgpu / wgpu-types / naga 29.0.3, glam 0.32.0, encase 0.12** | `bevy_render/Cargo.toml` @ v0.19.0 |
| V-2 | Bevy 0.19 **removed `RenderGraph`**; passes are systems in ECS schedules; non-camera work targets the `RenderGraph` schedule | 0.18→0.19 migration guide |
| V-3 | MC peak: **5.42 G voxel/s, 330 M tri/s** (RTX 2080 Ti). DMC costs 1.52–3.50×; FlexiCubes 2.77–3.92× | Grosso & Zint, `10.1007/s00371-021-02139-w` |
| V-4 | **Contouring 68 ms vs halfedge construction 58 ms** — extraction is 54% of a usable mesh | same, Table 5 |
| V-5 | On unstructured grids, Delaunay/MT ratio **15.3×–81.5×** — contouring is 1–2% of the pipeline | TetWeave Table 3 |
| V-6 | **73% of FlexiCubes' 64³ MC timing is fixed launch overhead** (fitted a ≈ 1.88 ms) | fit over FlexiCubes' own resolution series |
| V-7 | Cross-paper reproducibility floor is ~1.5× **in opposite directions**: TetWeave re-measured FlexiCubes at 128³ as 9.63/15.25 ms vs FlexiCubes' own 14.06/9.53 | both papers |
| V-8 | GPU MC throughput has not tracked hardware: **10.7× more bandwidth bought ~1.7× more throughput** (GTS 450 → 2080 Ti) | speed analysis |
| V-9 | Same MC, compute shader → mesh shader: **114.2 → 2679.4 fps (23.4×)** | Elliott MSc, Waikato 2022 |
| V-10 | CBT sum-reduction, atomics → LDS staging: **5.78 → 0.40 ms** | Unity SIGGRAPH 2021 (see ✗7) |
| V-11 | Meshlet compression: **15.5 M tri in 0.59 ms** (RX 7900 XTX) | `10.2312/vmv20241204` |
| V-12 | Work graphs: 79,710 instances in 3.74 ms — **but 2.8–3.4× slower** on classification workloads | `10.1145/3675376` + independent profile |
| V-13 | nvblox: meshing is the least GPU-accelerable stage, **×3–13 vs fusion's ×174–177** | nvblox |
| V-14 | Aokana renders **10¹⁰ voxels at 6 ms**, 5% resident, RTX 3060 Ti — **explicitly not editable** | Aokana |
| V-15 | CoACD vs V-HACD: **49% → 80%** downstream manipulation success | CoACD |
| V-16 | Dimforge migrated parry (0.26.0) and rapier (0.32) **off nalgebra onto glam**, citing rust-gpu support; performance *"nothing changed, at all"* | dimforge.com, 2026-01-09 |
| V-17 | **No paper since 2020 benchmarks MC vs Surface Nets vs DC against each other.** Surface Nets has no credible published timings at all | literature review round 1 |
| V-18 | **DC's own paper quantifies the f32 QEF failure.** At 256³, `bᵀb` reaches ~10⁶; f32 carries six decimal digits, so `E[x]` evaluated on a flat region — where it should be zero — has error **on the order of 1**. The paper's own remedy is double precision | Ju, Losasso, Schaefer & Warren 2002, `10.1145/566570.566586`, §2.3, read this session |
| V-19 | **DC's topology is Surface Nets' topology.** The paper's algorithm is literally: vertex at the QEF minimizer for each sign-changing cube, quad joining the four cubes of each sign-changing edge. Only vertex *placement* differs | same, §2.2 |
| V-20 | A QEF is stored as `AᵀA` (symmetric 3×3), `Aᵀb` (3-vector) and `bᵀb` (scalar) — 10 floats — rather than as `A` and `b` | same, §2.3 |

---

## Part 4 — Open questions

Each has the test that would settle it. **An open question with no proposed test is a wish.**

| # | Question | Settled by | Why it matters |
|---|---|---|---|
| O-1 | What fraction of cells actually change per brush stroke? | G-002 instrumentation; hash cell slabs, log per stroke | **Unpublished.** Ceiling on every incremental-repair idea in the opportunities doc |
| O-2 | Does clamping the QEF vertex to (1−ε) inside its cell eliminate self-intersections? | A-009: measure per 1,000 triangles, clamp on vs off, all seven fields | Decides whether guaranteed intersection-free extraction is free → whether runtime convex decomposition can stop failing |
| O-3 | MC vs SN vs DC vs MT — actual relative speed on one machine? | M-001 | The comparison does not exist (V-17). We'd have the only apples-to-apples measurement |
| O-4 | Do brush operations commute? | G-003: 8 ops × 40,320 orderings, count distinct results. Expect 1 | If not 1, the coordination-free multiplayer story dies — cheaply, before anything is built on it |
| O-5 | Do mesh shaders work on macOS/Metal? | GPU-007 capability probe | **Sources contradict:** wgpu's spec table lists MSL as *planned*; the tracking issue says the Metal HAL backend merged. Neither is trustworthy until probed |
| O-6 | What is amortized meshing cost per frame under continuous editing? | E-206 under a deliberately overloaded queue | The only number a game cares about, and no paper reports it |
| O-7 | What fraction of *our* pipeline is contouring vs everything else? | M-003 | V-4 says 54% for someone else's code with no physics. Ours is probably worse |
| O-8 | Does DC's vertex placement need f64 in practice, or is f32 enough? | E-112, with the QEF condition number in the HUD | `M = AᵀA` squares the condition number. **Half answered by V-18**: the original paper measures f32 error ~1 on flat regions at 256³, and recommends f64. Still open for *our* fields, *our* resolutions, and the closed-form three-plane rule — which sidesteps `AᵀA` entirely and may not degrade the same way |

---

## Part 5 — Method rules, and the failure that earned each

Rules with no incident behind them get ignored. These all have one.

| Rule | Earned from |
|---|---|
| A typed error at the call site is louder than an abort — make the invalid state unrepresentable where you can, report it where you can't, and never substitute a default | The no-panic rule, reconciled with "fail loudly": `ValidateConfig` has private fields and one checked constructor, so the validator needs no runtime guard at all |
| Corpus presence is decided by `catalog_read`, never by `distill_search` | ✗4 — 342 documents readable but unsearchable |
| **Never guess a DOI or arXiv ID.** Look it up or stop | A subagent guessed an ID from memory and downloaded an unrelated condensed-matter physics paper under a meshing DOI |
| Verify the *source* separately from the *number* | ✗7 — right figure, wrong attribution |
| "Nobody measured X" needs the same evidence as any other claim | ✗9 — asserted twice, false twice |
| No performance number without the benchmark that produced it, in the repo | The corpus contains several published figures that failed verification |
| A doc comment the test suite disproves is worse than no doc comment | ✗3 |
| Assert the identity, not the inequality — a weak assertion hides a strong fact | ✗1 |
| **Verify that a property test can actually fail.** Corrupt an input and confirm red | A test that cannot fail is decoration, not evidence |
| Record an assertion's break conditions *next to it* | ✗1 — G-001 chunking will break it correctly, and it will look like a regression |
| Pin known defects as non-zero assertions rather than excluding them | M-4 — the numbers only move when someone means them to |
| Single-grid timings measure dispatch latency; sweep resolution and report the fixed cost | V-6 |
| Treat any published cross-paper ratio below ~2× as noise | V-7 |

---

## Part 6 — Adding an entry

```markdown
### <tier><n> — <one-line claim>

**Believed because:** <where it came from — and be specific, "folklore" is a real answer>
**Tested by:** <the command, test name, or source read. Must be repeatable.>
**Result:** <the numbers>
**Consequence:** <what changed as a result — a decision, a ticket, an assertion>
**Would be shown wrong by:** <the observation that would falsify this>
```

Two things that make this file worth keeping rather than a chore:

**Re-tier rather than rewrite.** When an R becomes an M, leave the R text and add the measurement
below it. The gap between what was reported and what was measured is itself data.

**Record the ones we got right for the wrong reason.** ✗10 is in the falsified section even though
the outcome was fine, because the reasoning was wrong and the reasoning is what generalizes.
