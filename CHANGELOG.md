# Changelog

All notable changes to the isomesh workspace. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/); versions are the shared workspace version
that `isomesh` and `isomesh-gpu` release under. Everything here is pre-1.0: **0.0.x versions make no
stability promise**, and any release may break the API. Releases are CI-driven — a manifest version
bump landing on `main` is the release (`scripts/publish.sh`, version-driven).

## [Unreleased]

## [0.0.10] — 2026-08-27

**The project has a site, nine of the demos are playable in it, and the front page runs `isomesh` itself.**
[ladvien.github.io/isomesh](https://ladvien.github.io/isomesh/) renders this repository's own markdown,
serves every GIF and screenshot from itself rather than hotlinking `raw.githubusercontent.com`, and
carries nine of the examples as real WebAssembly builds. The three Phase 21 demos each still print their
cross-check against their committed CSV to the browser console, so a hosted build cannot drift away from
the artefact it illustrates. The front page carries a tenth module that is not a Bevy build at all — the
core crate with a hand-written WebGL2 renderer, **133,115 bytes against the Bevy modules' 36 MB** — which
is the size claim made checkable rather than asserted.

**It serves over HTTPS, and that is load-bearing rather than hygiene.** WebGPU is a
secure-context-only API, so while the site was reachable only over `http://` — the custom domain
`isomesh.ladvien.com` never got a certificate, and the `github.io` URL redirected to it — `navigator.gpu`
was `undefined` and all nine modules refused to start in *every* browser, desktop Chrome included. The
custom domain is gone and there is deliberately no `CNAME` in the artifact (`✗46`, `M-366`).

**One of those nine could not run at all before this**, and the reason was a defect in the published
crate rather than in the demo: `IsomeshPlugin`'s frame budget called `std::time::Instant::now()`, which
compiles on `wasm32-unknown-unknown` and then panics, so every browser game built on the plugin broke on
the first frame a chunk landed (`✗44`).

**`game_dig` is a game you can play on a phone.** Its terrain is a four-layer texture array blended by
slope and depth — grass on the up-facing surface, leafy dirt on the slopes, deep dirt inside a fresh
tunnel — the sandbox is visibly walled in concrete and the walls stop the player, the body is a
human-proportioned `1.70 × 0.50` capsule instead of two spheres `0.8` wide with a pinched waist, and a
hole you dig is escapable: the ground probe was a single point that read air whenever the body was wedged
with its foot on a ledge, which refused the jump (`✗45`, `M-363`). Touch controls feed the same movement
and edit paths the keyboard and mouse drive.

### Added

- **`scripts/build_web.sh`** — the one command that builds the site: nine
  `--profile wasm-release --target wasm32-unknown-unknown` examples through `wasm-bindgen`, then
  `isomesh_web`, then the prose. The `DEMOS=(…)` array is the single place that decides what is reachable,
  and `doc_facts.sh` holds `play.html`'s allow-list, its `#notes-` blocks and three prose sites against it.
  It reads the required `wasm-bindgen` CLI version out of `bevy_isomesh/Cargo.lock` and refuses to build
  with any other, because a CLI that disagrees with the crate emits glue for a different ABI and the module
  fails to instantiate in the browser naming neither tool.
- **`scripts/build_site.py`** — renders seven markdown sources plus a hand-written front door into
  `web/dist`, copies `docs/gifs`, `docs/screenshots` and `docs/experiments` beside them, and marks every
  `<img>` `loading="lazy"` (`DEMOS.md` alone references 34 GIFs, which is tens of MB on one scroll-free
  load). **It is also the link checker**: all 91 relative targets are resolved against the repository and
  anything that resolves to nothing fails the build naming its source, which is what makes rendering only
  part of the repository safe — links into `BACKLOG.md`, `docs/research/` and `crates/` are rewritten to
  github.com blob URLs instead.
- **A `site` CI job on every push and pull request**, and a `pages` job that deploys `main`. The `site`
  job is what turns the link checker and the wasm build into gates.
- **`isomesh_web`** — a third workspace whose only product is the front page's wasm module. `extern "C"`
  over linear memory with no `wasm-bindgen`, because `MeshBuffer`'s fields are `pub` precisely so a
  consumer can read them without a copy; sixteen exports, **zero wasm imports**, eight fields, five
  extractors, and `isomesh::validate`'s report recomputed per re-mesh. `#[unsafe(no_mangle)]` is the only
  `unsafe` token it may contain and `build_web.sh` refuses to build if an `unsafe` block, `fn`, `impl` or
  `trait` appears in `isomesh_web/src` — a checkable rule where the `unsafe_code` lint cannot draw one,
  since edition 2024 has no way to emit a wasm export without that attribute.
- **`web/lite.js`** — 300 lines of hand-written WebGL2 driving that module: one program, an orbit camera,
  a wireframe pass over a line index buffer built from the triangles, and a HUD of the validity counters.
  `UNSIGNED_INT` indices are why it is WebGL2 and there is deliberately no WebGL1 path. Every typed-array
  view is re-created after every `iso_mesh` call, because memory growth detaches `memory.buffer` and a view
  captured earlier reads a dead buffer as zeros.
- **Six more playable demos** — `quickstart`, `marching_cubes_tunnel`, `dual_contouring_cube`,
  `surface_nets_vs_marching_cubes`, `game_dig` and `game_showcase`, each with its controls and its claim on
  `play.html`. `game_dig` and `game_showcase` read `ISOMESH_*` environment variables natively and a browser
  has no environment, so both run at their own documented defaults; the page says so rather than pretending
  otherwise.
- **Four gates in `scripts/doc_facts.sh`**, all four demonstrated failing before being left passing: the
  playable-demo count against three prose sites, `play.html`'s allow-list against `build_web.sh`'s array,
  a `#notes-<name>` block per allow-listed demo, and the front page's "seven extractors sit behind one
  trait" claim, which was ungated until now.
- **`scripts/preflight.sh` runs the site build and `isomesh_web`'s own three gates.** `build_site.py` is
  the repository's link checker and ran only in CI, so the gate most likely to break on a docs edit was the
  one a local run could not see. `isomesh_web` is its own workspace, so `cargo clippy --workspace` cannot
  reach it — E-111's lesson, applied to a third workspace.

### Fixed

- **`IsomeshPlugin` no longer panics in a browser (`✗44`).** `bevy_isomesh/src/plugin.rs` called
  `std::time::Instant::now()` in `apply_finished_meshes`, the system that spends the per-frame mesh budget
  — so the plugin died on the first frame a chunk landed, taking `quickstart`, `game_showcase`,
  `game_terrain_stream`, `game_walk` and `game_capsule_walk` with it. `bevy_platform = "0.19"` is now a
  direct dependency with its `web` feature enabled for `cfg(target_arch = "wasm32")` only, and `Instant`
  comes from `bevy_platform::time`. `Duration` stays on `std`; it carries no clock. Verified in the
  resolved graph rather than by inspection: `web-time v1.1.0` is in the **normal** wasm dependency tree
  under `bevy_isomesh` and absent from the native one, and `bevy_platform` appears exactly once in the
  lockfile. **This is a fix to the published crate**, not to a demo — nothing in its README or in the
  plugin's own documentation said "native only".
- **The wire-size figure on the site was wrong.** `web/index.md` and `web/play.html` both said "about
  8.4 MB on the wire"; the modules gzip to 8.79/8.76/8.76 MB and the live origin transfers **8,844,921
  bytes** for one of them. Both strings now say 8.8 MB.
- **`web/play.html` opened with no `?demo=` printed `unknown demo: null`.** `URLSearchParams.get` returns
  `null` when the key is absent and that was folded into the unknown-demo branch, so a bare visit — a URL
  shared without its query string, a bookmark, or the verification URL this repository's own docs quote —
  landed on a null nobody typed, over a black canvas, with no link to any demo. It is a landing page now:
  absent renders an index of the nine, an unrecognised name renders an accurate message *and* the index.
  Built from the allow-list with `createElement`/`textContent`, never from the query string.
- **The WebGPU gate passed on a browser with no GPU adapter and left a black canvas.** It tested only that
  the `navigator.gpu` object exists; `requestAdapter()` returning `null` — a blocklisted driver, a VM, a
  headless run — passed it, `await init()` resolved because `spawn_app` returns immediately on wasm, and
  Bevy then failed asynchronously with nothing on screen to read. Reproduced on the live site and fixed by
  awaiting the adapter (`✗46`).
- **`game_dig`'s dug holes were inescapable.** The ground probe sampled one point below the foot, which
  reads air whenever the body is wedged with its centre over a void and part of its foot on rock — so
  `grounded` went false and `Space` was refused. It now samples the foot sphere's own lower surface
  (`✗45`).
- **Digging out the bottom of `game_dig`'s sandbox dropped the player through the floor (`✗48`).** `aim`
  refuses a brush *centred* outside the box, but a brush centred **on** the floor plane reaches its own
  radius below it, so a shaft to the bottom removed the field under a box whose meshes stop at
  `y = -5.4`. The five boundary slabs are now part of the field the body resolves against rather than
  scenery with a position clamp beside it — which is also what makes the floor *standable*: a `y` clamp
  stops a fall without ever reporting ground, so it would have left the body hanging with the jump
  refused.
- **`GPU_JOBS_MAX` did not bound the GPU job queue (`✗47`).** The cap was the second clause of a budget
  predicate that `DirtySet::mesh_within_budget` consults *after* meshing a chunk — deliberately, so a
  too-small budget still progresses — so a frame beginning at the cap added one anyway, every frame in
  which no readback retired. Measured 17, then 18. The cap is now checked before the dispatch and the
  refused chunk goes back in the dirty set. Separately, `gpu_collect` now runs **before** `drain_dirty`,
  so the in-flight count is not read before the thing that shortens it (`M-368`).
- **`game_dig`'s frame rate vanished with its numbers panel.** It opens with the panel hidden, so the
  steady state showed the mesher on its banner and nothing about what it cost. The frame rate now joins
  the one line a hidden panel leaves. `ISOMESH_VIEW=nohud` still produces an empty frame, which the
  committed GIFs depend on.

### Changed

- **`game_edit_tape_trim` measures one chunk per frame instead of fanning across `std::thread::scope`.**
  Thread spawn panics on wasm and a static host cannot send the COOP/COEP headers `SharedArrayBuffer`
  needs, so the threaded sweep could not exist there; doing all 1,571 re-meshes inside `Startup` would
  freeze the frame for twenty seconds instead. `measure_chunk` is untouched and each chunk is independent,
  and the cross-check is what says the result is the same: **all 13 comparisons reproduce, including all
  64 per-chunk mesh hashes**, native and in a browser. The startup report drops the thread count and its
  17 `println!` calls became `info!`, because `println!` writes to an unsupported stdout on wasm and is
  discarded while Bevy's `LogPlugin` routes `tracing` to `console.log`.
- **`Instant` in every web-built example is `bevy::platform::time::Instant`**, which is
  `std::time::Instant` natively and `web_time::Instant` on wasm. `std::time::Instant::now()` panics there,
  and Bevy's `web` feature — now enabled for `cfg(target_arch = "wasm32")` only — is what makes the
  re-export resolve to the working one. `Duration` stays on `std` everywhere: it is arithmetic, not a clock.
- **The examples' shared `F12` screenshot names files from a counter, not `std::process::id()`**, which
  compiles on wasm and then panics with "no pids on this platform" — and
  `Window::prevent_default_event_handling` defaults to `true`, so `F12` in a browser reaches the app rather
  than devtools. Repeated presses now write `screenshot-1.png`, `screenshot-2.png`, … instead of one name
  per process.
- **`avian3d`, `parry3d`, `wgpu` and `isomesh-gpu` are native-only dev-dependencies.** None is reachable
  from the nine examples that get a web build, and `isomesh-gpu` asks wgpu for `dx12`, `metal` and
  `vulkan` and for no web backend at all. Dev-dependencies are not propagated, so consumers see no change.
- **`quickstart` binds the page's canvas**, which is the one web-specific line in the one example with no
  `WindowPlugin` — it had none because its whole point is being the shortest path, and it stays the
  shortest path with six lines added and nothing removed.
- **`game_dig`'s terrain is a four-layer texture array, and its player is a person.** `triplanar.wgsl`
  samples `texture_2d_array<f32>` and blends grass, leafy surface dirt and deep dirt by slope and world
  height, with `Triplanar.settings.z` — previously unused — selecting one layer outright for the five
  concrete slabs that now line the sandbox and stop the player at its boundary. The body went from two
  `0.4` spheres touching at a point (1.6 tall, **0.8 wide**, with a pinched waist) to four `0.25` spheres
  that overlap: **1.70 × 0.50**, and narrower than the cavity its own default brush carves.
- **`game_dig` is playable by touch.** The left third of the screen is a virtual stick, the rest drags to
  look, and three on-screen buttons carve, fill and jump — all feeding the same movement and edit paths
  the keyboard and mouse drive, with no second edit route to keep in step. The buttons stay hidden until a
  touch is seen, so a desktop run is visually unchanged.

### The four things a consumer upgrading from 0.0.9 is actually upgrading into

Everything above is the site and the demos. These are the library, and they are what a `cargo update`
lands.

- **Vertex positions move, on the default path, on every field.** A crossing is now stored as a **signed
  offset from the edge midpoint** (`cube::edge_offset`) rather than as a parameter from the lower corner.
  `d = ((a + b)/2)/(a − b)` is exactly antisymmetric under the simultaneous endpoint-and-sign swap by four
  IEEE 754 guarantees, which makes plain Marching Cubes bit-exactly equivariant under **all 48 octahedral
  elements instead of 6** — 0 mismatches on 9.2 M straddling pairs against 1,035,808 for the old form.
  **135 of 216 golden hashes were rebaselined.** Triangle counts are unchanged, Hausdorff distance and
  self-intersection counts are identical to twelve digits, and 2,285 of 28,124 cut edges moved by at most
  268 ULP. If you hash meshes, **0.0.9's hashes will not match.** There is no second path and no flag:
  `cube::edge_crossing` is gone, all six placements and the WGSL shader are in the centred frame.
- **A release-mode panic is fixed.** `MAX_PATCH_TRIANGLES` was a **sampled** maximum of 24 while the
  triangulator's own buffer was the **derived** 40, so an ordinary trilinear cell emitting 26 triangles
  indexed out of bounds — a panic in release, not an `Error` and not a hole. Found by an exhaustive sweep
  over all 2¹⁸ sign patterns of a four-cell block.
- **New public items.** `cube::edge_offset`, also re-exported from `marching_cubes::table` beside
  `is_inside`. In `isomesh-gpu`: a `wgpu` re-export so a consumer cannot end up with two `wgpu` majors;
  `StageTimestamps`, `Span`, `Spans` and `MAX_PASSES` for GPU-side pass attribution; `Gpu::with_timestamps`
  and `Gpu::with_subgroups`; and `DeferredGeometry`, a keyed queue of in-flight geometry read-backs. That
  last one is a **third** extraction contract beside the two that already ship — `extract_buffers` promises
  geometry now, `extract_indirect` promises it never leaves the device, `DeferredGeometry` promises it a
  frame or two later with the wait amortised — not a fallback for either. Measured under
  `DirtySet::mesh_within_budget` at **1.41 frames of latency, worst case 2**.
- **Breaking within `0.0.x`.** `isomesh_gpu::Error` gains `FeaturesUnsupported { missing }` and
  `DeferredQueueFull { capacity }`, so a match on it is no longer exhaustive. `AdapterReport` gains
  `subgroup_min_size` and `subgroup_max_size`. And `Gpu::open`'s feature check no longer reports
  `TimestampsUnsupported` for a missing **non-timestamp** feature — it was correct while timestamps were
  the only requestable capability and became a lie the moment `with_subgroups` existed.

## [0.0.9] — 2026-08-17

**The air region's connectivity is now maintained in both directions, and across chunk seams.** 0.0.8
shipped `connectivity::Air` as a dig-only structure with a note saying a `fill` could not exist. It can;
the note was right about the union-find and wrong about the problem.

**Breaking:** `Repair::unions` is now `Repair::relabels` and `Repair::unions_per_dirty` is
`relabels_per_dirty`. There are no union calls any more, so the old names stopped denoting anything.
`Air::connected` and `Air::components` now take `&self` rather than `&mut self`, which is a relaxation
and breaks nothing.

### Added

- **`connectivity::Air` now fills as well as digs**, and `connected` is `O(1)` and takes `&self`.

  `Air` was a union-find, and a union-find **cannot** absorb deletion — parent pointers encode union
  *history*, not spatial adjacency, so re-rooting a shed piece severs its **descendants** from a
  component they are still part of (✗26). The structure is now a flat label array, which is what makes
  a fill one write per shed member. `components()` is `O(1)` maintained rather than an `O(n)` scan, and
  `Repair::unions` became `Repair::relabels` because the old name stopped denoting anything.

  **`fill` costs the shed volume, not the component** (M-321): **3.09 → 3.38 voxels visited per seed**
  while the lattice grows 7.6×, and **436×** less work than rebuilding at 65³. The replacement search is
  lockstep — every seed grows a frontier, frontiers that meet merge, and the walk stops when all but one
  exhausts, so the surviving piece is never walked to completion.

  **The tail is real and is documented rather than hidden.** Bisecting a tunnel between two equal
  caverns costs **1.1× a full rebuild**: both frontiers are half the component, and there is no
  replacement edge to find because the component genuinely split. HDT's levels are not the remedy — they
  bound a *search*, and the remedy is decomposition, which is what `connectivity::AirWorld` below is.

- **`connectivity::AirWorld` — many `Air` grids over one `ChunkLayout`, so `connected` answers across a
  chunk seam.** Adjacent chunks share exactly one sample plane (`sample_shape` is `cells + 1`), so two
  components are joined where that shared sample is air on both sides: nothing interpolated, nothing
  matched by tolerance.

  **This is what makes `fill` usable, not an optimisation on top of it.** Without it a consumer chooses
  between one large `Air` — which has the bisect tail — and one `Air` per chunk, which cannot answer
  across a seam at all. Every extractor here is driven per chunk, so the second is the natural shape and
  it silently cannot answer the question the module exists for.

  **Measured (M-322): the bisect becomes bounded rather than cheap.** A single grid visits every air
  sample it has — 557,568 at 16 chunks wide, 0.998× its own rebuild — while the chunked world is **flat
  at 34,848** as the world grows 8×. The advantage is exactly the chunk count and grows without bound.
  But 34,848 is still 0.970× a *chunk* rebuild: chunking converts an unbounded cost into one bounded by
  a unit the mesher already budgets per edit.

  **A `sealed_cave` example demonstrates it** — two chambers in different chunks, a tunnel through the
  chunk between them, `F` to plug it. It is also what produced M-323: the bound and the cost are
  different claims. Chunking bounds the search by the chunk, but what it *costs* is the edited chunk's
  share of the severed component — **0.03× a chunk** here, where the chambers live elsewhere and the
  boundary graph resolves the global split, against **0.97×** when both halves sit inside the edited
  chunk. Same operation, 35× apart, decided by geometry.

  The global component graph is **rebuilt from scratch** on every restitch, which is the one thing a
  union-find is safe to do after ✗26 — it only ever unions — and is affordable because its nodes are
  components rather than samples. The `O(cells²)` seam scan is cached and recomputed only for seams
  touching a chunk that changed.

  **Repair is budgeted, not synchronous**, taking the same `spend: FnMut() -> bool` predicate
  `mesh_within_budget` uses — because amortised is the wrong statistic for the frame a breakthrough
  lands on. **Both directions of staleness are conservative:** an unfinished fill reads *"not sealed
  yet"* and an unfinished dig reads *"not connected yet"*. Water leaking for three frames is
  recoverable; water not leaking out of a room the engine wrongly believes is sealed is a broken game
  rule.

## [0.0.8] — 2026-08-17

Phase 17's measured results, shipped. Two of them answer questions this crate had been asking without
being able to check: **does the mesh separate what the field separates** (`validate::sealing`), and
**what does the air region's connectivity cost to maintain rather than rebuild** (`connectivity::Air`).

**One behaviour change, and it is not confined to the type it touches.** `construct::SampledField`
now supplies the exact trilinear gradient, which moves **Dual Contouring's output on sampled fields**.
Reference-field extraction is untouched and every golden hash is unchanged — those fields carry their
own analytic gradients and never used the default — so this is visible only to callers who mesh a
sampled volume. Numbers below.

### Added

- **`isomesh::connectivity::Air` — connected components of the air region, repaired as you dig.**
  *Is this cave sealed? Did I just break through?* are questions about the connected components of the
  air sublevel set, asked after every edit. Digging removes solid, so air samples only ever **appear**,
  and an insert joins at most two trees with no replacement-edge search — so a union-find is the entire
  structure.

  Measured (M-311): one brush of fixed radius into lattices of 33³, 65³ and 129³ costs **4,872 union
  operations at every one of them**, while the lattice grows 59.7×. A rebuild pays **2,146,689 samples
  scanned to discover the 925 that changed**, which is a 104× wall-clock gap at 129³ that widens with
  the lattice.

  **There is deliberately no `fill`.** Removing air is a deletion, a union-find cannot do deletions at
  any price, and an API that offered one and silently rebuilt would be a second execution path.

- **`MeshReport::mean_ratio` and `MeshReport::irregular_vertices`** — triangle-shape quality, in the
  definition the isosurfacing literature reports so the columns can be read beside it. Recorded, never
  gated, the same standing `degenerate_triangles` has.

- **`SweptFaces::margin`** — the value `Interior::test` is the sign of. `test()` now literally calls it
  and compares with zero, so the interior ambiguity decider is the `ε = 0` member of a one-parameter
  family **by construction**. Bounded by half the field's scale (M-312). Note it is a **decision
  margin**, not a persistence, and thresholding it does **not** resolve the cells where the published
  algorithms disagree — measured, and the overlap is *below chance* (M-313).

- **`marching_cubes::interior::chernyaev_numerator_test` is public** and no longer test-only, so the
  comparison between the corrected interior test and the construction it corrects is reproducible
  outside this crate's own suite.

- **`scripts/fetch_volumes.sh`** — fetches real scanned volumes from Open SciVis and verifies them
  against the **publisher's own SHA-512**. The data is gitignored; `docs/measurements/volumes/PROVENANCE.md`
  is committed. Benchmarks that read them skip cleanly when they are absent, so a clean clone with no
  network still builds.

- **`isomesh::validate::sealing` — does the mesh separate what the field separates?** Every other
  validity metric in this crate judges a mesh against itself (manifoldness, orientation, Euler
  characteristic) or against the field's geometry (`validate::accuracy`). None asked whether the mesh
  partitions *space* the way the field's sign does, and neither claim implies the other: a mesh can be
  closed, manifold, correctly wound and Hausdorff-close while sealing a passage the field leaves open.
  `SealingReport` reports holes, membranes, the two air-component counts, and how many holes touch a
  face of the sampled domain.

  The measurement it produced (M-307, `docs/experiments/p-21.csv`): **Marching Cubes, Marching Cubes +
  decider and Marching Tetrahedra seal all eight reference fields at all three resolutions.** All three
  dual methods leave `fbm_terrain`'s domain boundary open, with identical counts, and **every one of
  those holes is on a domain face** — a dual emits one quad per sign-changing grid edge and that quad
  needs all four cells around it. **For a chunked world that face is the chunk seam.**

  The test itself is Wojtan, Thürey, Gross & Turk's complex-edge test (`10.1145/1778765.1778787`) and
  is cited as theirs; running it as a correctness audit of extraction is what is new.

### Fixed

- **Subgrid Marching Tetrahedra no longer refuses a whole volume for want of one normal.** Where the
  field has a critical point *on* the isosurface, no normal exists — and if the isosurface passes
  through it, the level set is genuinely singular rather than merely awkward. That tetrahedron is now
  **skipped**, leaving a hole its size, and `SubgridMarchingTetrahedra::report()` returns a
  `NormalReport` giving the count, the cause, and **the position of every one**.

  **Positions, not just a count**: a count can stay the same while the sites move, so a count alone is
  not a regression test, and a caller repairing the mesh needs to know where. The cause separates
  `Degenerate` (gradient exactly zero — a critical point, no normal exists) from `IllConditioned`
  (non-zero but below the conditioning floor — a normal exists and is not trustworthy), because those
  have different remedies and folding them together would let a precision bug hide inside a topology
  count.

  **Nothing is substituted.** A wider stencil would return the gradient of a *smoothed* field, which is
  a different field, and at a saddle there may be no correct normal at all.

  **If your data is integer, contour at a half-offset isovalue** — `127.5` rather than `127`. Integer
  samples cannot equal a half-integer, so no sample sits on the isosurface and the degeneracy never
  arises. On `bonsai`, **3% of surface-cell corners sit exactly on the surface** against an integer
  isovalue (M-316). Standard practice in volume rendering, for this reason.

- **`construct::SampledField` now supplies the exact trilinear gradient** instead of inheriting
  `Sdf::gradient`'s central difference. **A central difference is identically zero at a local extremum,
  however steep the field is around it** — and quantised data manufactures those: on the `bonsai` CT
  volume, corners with neighbour slopes of ∓19 came out exactly symmetric because `u8` quantisation put
  both neighbours on the same integer. The subgrid extractor asks for a normal there and **refused the
  whole volume** (A-028, M-316).

  **This changes Dual Contouring's output on sampled fields**, since its QEF uses gradients: on `bonsai`,
  529,488 → 529,383 vertices and 1,776 → 1,770 non-manifold edges. Reference-field extraction is
  untouched — those fields carry their own analytic gradients and never used the default — and all
  golden hashes are unchanged.

- **A wall-clock assertion inside a unit test failed the 0.0.7 release on a macOS runner.**
  `empty_cell_rejection_is_measured_per_field` asserted a speedup above 1.0 on every reference field,
  under a doc comment that called itself "not a regression gate". `gyroid` is triply periodic, its
  surface reaches nearly every cell, and empty-cell rejection has almost nothing to reject there — 16.8%
  of cells against 80.6–95.1% on every other field. The gate is now that count, which is an integer and
  the same on every machine; the ratio is printed.

## [0.0.7] — 2026-08-16

### Documentation

- **Every README, both demo pages and the experiments page brought up to the shipped API.** None of
  them mentioned the exact predicates, the public validity gate, the attribute-preserving weld or the
  Bevy `Mesh` interop. The experiments page stopped at P-17 and now carries P-18, P-19 and P-20 —
  including the falsification, which is the point of that page existing.

- **A `weld_creases` example** (`cargo run --example weld_creases --release`): two cubes, same input,
  same tolerance, one welded on position alone and one with a normal key. The left loses the flat
  shading that made it read as a cube; the right keeps it. Both welds are correct, which is the thing
  the example exists to show.

- **`doc_facts.sh` now checks the example count it has always derived.** It computed the number,
  printed it, and never compared it to the prose — so both READMEs said "34 examples" while the
  directory held 35. That is precisely the rot the script was written to stop, sitting inside the
  script. The first version of the fix matched a bare `examples` and fired on `isomesh-gpu`'s
  perfectly true "Three examples", which is the script's own header warning — *"adding a loose one
  costs everybody who runs this"* — landing on its author.

### Added

- **`bevy_isomesh::weld_keys` — a weld key from a `Mesh`'s normals and UVs.** Feeds
  `Welder::weld_split_by`, so the crease that makes a cube look like a cube survives a weld:
  `Mesh::from(Cuboid)` keeps all 24 vertices with the key and collapses to 8 without it. **It takes a
  quantum, not a smoothing angle, and that is forced rather than preferred** — the conventional
  "within 30°" test is not transitive, so it is not an equivalence relation, and applying it to a
  `k`-way coincidence class is E×4's manufactured-bowtie failure. Quantising to a lattice is
  transitive; the cost is a missed merge at a bucket boundary, which is a seam rather than a topology
  defect. Defaults are stated as conventional, not derived. Hashing is FNV-1a spelled out, because
  `DefaultHasher` is not stable across Rust releases and a drifting key silently changes the mesh
  (B-014).

- **`Welder::weld_split_by` — a weld that refuses to merge vertices whose caller-supplied key
  differs.** One `u64` per vertex, or an empty slice meaning "one class", which is exactly what
  `weld` passes; one implementation, two entry points. **The parameter is a key and not a predicate on
  purpose**: E×4 gated the weld on the pairwise link condition and was reverted as strictly worse —
  over 56 configurations it removed at most 4 non-manifold edges and *added up to 791 non-manifold
  vertices*, because a `k`-way coincidence is manifold only if all `k` merge and a pairwise test
  leaves the odd one out a bowtie. Equality on a key is an equivalence relation, so it partitions each
  class into complete sub-classes and that failure is **unrepresentable in the signature** rather than
  merely discouraged (R-010).


- **`bevy_isomesh::from_bevy_mesh` — a Bevy `Mesh` as the `(positions, indices)` pair every `isomesh`
  entry point takes.** The inbound half of a conversion whose outbound half already shipped.
  `Indices::U16` widens to `u32`, so callers never handle both. **It does not weld**:
  `Mesh::from(Cuboid)` returns 24 positions over 8 distinct corners, and collapsing them would destroy
  the crease that makes a cube look like a cube — B-014 exposes the predicate that decides it. Returns
  a `#[non_exhaustive]` `SoupError` rather than logging and skipping, because a scene walker needs to
  know which mesh it skipped and why. Nothing is repaired: a `TriangleStrip` is refused rather than
  expanded, since expanding one silently flips every other triangle's winding (B-012).

- **`validate::SurfaceGate` and `MeshReport::satisfies(gate)` — the rule for *which* validity check
  applies, which until now was compiled out of every shipped build.** `MeshReport` offered three
  predicates and no reachable statement of which one belongs to what, so consumers re-derived it and
  got it wrong in the obvious way: calling `is_closed()` on a render mesh that was never a solid, and
  reading the failure as a mesher defect. The tag is data, the method is policy; both ship. The enum
  is `#[non_exhaustive]`, so a fourth case later is not a breaking change (T-023, ✗22).

- **`manifold_check` now says whether a mesh met the gate it was *supposed* to.** Its descriptive
  cascade is unchanged — "what is this mesh" is a correct use of the three predicates — but a verdict
  sits beside it now. Surface Nets and plain Dual Contouring earn `ClosedAllowingUnresolvedTopology`
  on a closed field rather than `Closed`, because one-vertex-per-cell is *legitimately* non-manifold
  at coarse resolutions (M-4, M-15).

### Fixed

- **The generalized-winding backend is `O(N³·B)` in boundary-edge count, and nothing said so.** A
  hole-punched sphere at `65³` takes **43.6 seconds** to batch, against 0.35 s for the same mesh
  closed. Both factors grow, so the cost on genuinely damaged input — which is the input this backend
  exists for — is far worse than any closed-mesh benchmark shows. Measured, not yet mitigated
  (M-303).

- **`MeshField`'s docs claimed Manifold Dual Contouring is a sparse consumer of it. It is not.**
  `ManifoldDualContouring::extract` pre-samples all N³ grid points into a buffer before visiting any
  cell, so it reads a grid like every other extractor here. The claim was lifted from a summary of the
  paper rather than checked against the code, and it had reached shipped documentation on a public
  type. The type is unaffected and stays — it is still the one implementation of the query, and T-025
  routed the winding path through it — but the motivating consumer named in its docs was wrong
  (D-011, ✗23).

## [0.0.6] — 2026-08-16

The dual path got **4.26× faster with byte-identical output**, a published manifoldness claim was
falsified, and a bug in a reference implementation was found to have recorded two true hypotheses as
false. Exact geometric predicates arrived.

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
