# Is this for you, and its status

Moved from the root README on 2026-10-31, with its link paths updated and stale Bevy versions corrected, when the README was brought to the [org documentation guide](https://github.com/bevy-coven/.github/blob/main/DOCUMENTATION.md).

## Is this for you?

| You want | Honest answer |
|---|---|
| a chunked voxel world — dig, build, stream | **yes**, the primary target; use `MarchingCubes` or `SubgridMarchingTetrahedra`, the two extractors measured to tile across chunk seams with zero open edges |
| CAD-grade `f64` end to end | **yes**, the second target; every algorithm is generic over the scalar, and `precision_f32_vs_f64` shows exactly where `f32` tears |
| `no_std` | **yes** — the core is `no_std + alloc` unconditionally, with one dependency |
| sharp corners *and* chunk streaming together | **not yet** — `DualContouring` holds corners but does not tile across seams (measured, structural); pick one per volume |
| GPU extraction to make CPU meshing faster | **yes, above about 33³** — and by **37×** at 129³ (0.54 ms against a single-threaded CPU's 20.14). Below that the fixed cost wins and the CPU is faster. The catch is *where the field is evaluated*: sample on the CPU and upload, and you get 2.4× at 129³ because the upload is 87% of the path; evaluate it in the shader and the upload disappears (M-155) |
| MC33 tunnels through a single cell | **yes** — opt in with `set_interior_ambiguity(InteriorAmbiguity::Trilinear)`; a tunnel is meshed against the inner hexagon and χ falls by exactly two per tunnel, checked against a tunnel count taken from the classifier rather than the mesh. **One configuration is refused rather than meshed**: a cell whose contours run past Grosso's Corollary 6 bound has no published triangulation, so `extract` returns `Error::UnresolvedSixSaddle` instead of emitting a hole (A-002b, A-020) |
| convex decomposition for physics | **not here** — export the mesh and decompose downstream; `game_destruction` shows the handoff |

## Status

Early. **Seven** extraction algorithms — including one that resolves features thinner than a voxel, which nothing else here can do — three normal-estimation strategies, a validity harness, an accuracy harness, a measured shootout between them, collider readiness, field-derived LOD, Transvoxel seams, chunk streaming, and a Bevy plugin that meshes off the main thread. [`BACKLOG.md`](../../BACKLOG.md) is the queue and the state; its archive carries every completed ticket with what implementation changed about it.

| | |
|---|---|
| **Working** | Marching Cubes · **Marching Cubes 33's asymptotic decider** · **MC33's interior ambiguity (tunnels and the twelve-vertex contour)** · Marching Tetrahedra · Surface Nets · **Dual Contouring** · **Manifold Dual Contouring** · greedy quads · Hermite data · mesh validity harness · accuracy harness · **seven-algorithm shootout** · chunk coordinates · dirty-set re-meshing · brushes · self-intersection counter · determinism harness · fourteen reference fields · property tests · vertex welding · **collider readiness** · **field-derived LOD** · **Transvoxel transition cells** · **frame-budget scheduling** · **subgrid Marching Tetrahedra** · **chunk streaming with hysteresis** · Bevy 0.19 bridge and plugin · **GPU compute Marching Cubes** · **GPU prefix scan** · **GPU field evaluation** · **mesh-shader rendering** · distance-field construction · **paint that lives in the edit log** · orientation repair · **pre-registration the compiler enforces** · **exact geometric predicates (`orient2d`, `incircle`)** · **a public validity gate you name rather than guess** · **attribute-preserving welds** · **Bevy `Mesh` → triangle soup** |
| **Not yet** | the singular face — a saddle lying *on* a cell face, which quantised input reaches and continuous `f64` reaches wherever the field's mirror symmetry is the grid's — 72 faces on two of the fourteen reference fields at 17³ (A-002i, R-187) · convex decomposition |
| **Deliberately absent** | any math library in the public API · any `bevy` mention under `crates/` · any performance number without a committed benchmark |

One optional feature: `experimental`, which adds `ProbabilisticQuadric`. It is off by default and on for docs.rs, so you can read what is behind it without checking out the source to find that it exists.

Published on crates.io: [`isomesh`](https://crates.io/crates/isomesh), [`isomesh-gpu`](https://crates.io/crates/isomesh-gpu) and [`bevy_isomesh`](https://crates.io/crates/bevy_isomesh). Releases are CI-driven: `scripts/publish.sh` is version-driven and uploads only what crates.io does not already have, so a version bump landing on `main` is the release.
