# isomesh

**Engine-agnostic isosurface extraction in Rust. Signed distance field in, triangles out.**

`isomesh` has to serve both a real-time voxel game and a CAD tool. That single constraint decides almost
everything about it: no math library appears in a public signature, output buffers are caller-provided and
reusable, the scalar type is generic over `f32` and `f64`, the core crate is `no_std` with exactly one
dependency, and seven extractors sit behind one trait. Every number on this site was produced by a command
in [the repository](https://github.com/ladvien/isomesh), on a machine named in
[`FINDINGS.md`](FINDINGS.md)'s header, rather than quoted from a paper.

> ⚠️ **Vibe Coded.** This crate was written by an AI agent working from a research corpus and a ticket
> queue. Every number here is produced by a test in the repository and every algorithm cites its source,
> but it has not been through human code review. Read the tests before trusting it with anything that
> matters.

---

## Try it right now

This is `isomesh` itself, compiled to WebAssembly, meshing in your browser at whatever resolution you drag
the slider to. **About 130 KB**, because the core crate is `no_std` with one dependency and the renderer
below it is 300 lines of hand-written WebGL2 rather than a game engine. Every counter in the panel comes
from [`isomesh::validate`](https://github.com/ladvien/isomesh/blob/main/crates/isomesh/src/validate.rs) —
the same report that gates every test in the repository — recomputed on your machine on every re-mesh.

Two of them are worth watching rather than reading. **χ** is the Euler characteristic: `2` for a sphere,
`0` for a torus, and steeply negative for the gyroid, which is what a high-genus surface means numerically.
**Non-manifold edges** should be zero, and on `gyroid` under Surface Nets or Dual Contouring it is not —
that is the measured, documented cost of placing one vertex per cell rather than a defect in this page.

<div class="lite-controls">
  <label>field <select id="lite-field"></select></label>
  <label>extractor <select id="lite-extractor"></select></label>
  <label>resolution
    <input type="range" id="lite-samples" min="9" max="49" step="2" value="33">
    <span id="lite-samples-value">33³</span>
  </label>
  <label><input type="checkbox" id="lite-wireframe"> wireframe</label>
  <span>drag to orbit · scroll to zoom</span>
</div>
<canvas id="isomesh-lite" width="1280" height="720"></canvas>
<p class="lite-hud">
  <span>vertices <b id="lite-vertices">–</b></span>
  <span>triangles <b id="lite-triangles">–</b></span>
  <span>χ <b id="lite-euler">–</b></span>
  <span>non-manifold edges <b id="lite-non-manifold">–</b></span>
  <span>boundary edges <b id="lite-boundary">–</b></span>
  <span>degenerate triangles <b id="lite-degenerate">–</b></span>
  <span>extraction <b id="lite-ms">–</b></span>
</p>
<p id="lite-error" hidden></p>
<script type="module" src="lite.js"></script>

The extraction time is measured by the browser around the one call into the module, which is why the module
itself keeps no clock and imports nothing at all — not one function, from anywhere. `fbm_terrain` reports
non-zero **boundary edges** on purpose: it is a heightfield and it exits through the sides of its own
domain, so a closed-surface gate would be the wrong gate to hold it to.

---

## Play in your browser

Nine of the demos are playable in your browser: real WebAssembly builds of the same examples the
repository ships, no video and no replay. Three of them — the Phase 21 ones — **measure their own figures
in your browser** and then print a cross-check against the CSV they were registered with; open the
developer console and you can read the comparison, line by line, on your own hardware. A demo whose numbers
stopped agreeing with its committed artefact would say so on screen.

Each module is 36 MB, about 8.8 MB on the wire once the host has gzipped it, and it is cached after the
first visit. Each also carries its own copy of the Bevy runtime, which is where almost all of that goes —
the demo at the top of this page is the same library without it, at 130 KB.

### [▶ Play `game_showcase`](site:play.html?demo=game_showcase)

![Flying through a landscape riddled with caves, arches and tunnels](docs/gifs/flying-through-the-rock.gif)

Every other gameplay demo runs on a heightfield, which is exactly the case you do *not* need a voxel mesher
for. This one is fBm terrain **intersected** with a thickened gyroid, so the camera flies through arches and
out the far side of walls. It is also the demo that streams chunks through `IsomeshPlugin`, and the plugin
called `std::time::Instant::now()` — which panics in a browser — until `✗44` was found building this page.

`M-128` · `E-210` · the chunk seam is why this streams on Marching Cubes: **0 boundary edges in the shared
plane against Surface Nets' 5 and Dual Contouring's 4**

### [▶ Play `dual_contouring_cube`](site:play.html?demo=dual_contouring_cube)

![Surface Nets rounding a box corner beside Dual Contouring holding it, across a resolution sweep](docs/gifs/dual-contouring-vs-surface-nets.gif)

The same field, the same grid, the same crossings, and one function different between the two meshes: where
a cell's vertex goes. Distance from the true corner to the nearest vertex, recomputed live on every
re-mesh — **Surface Nets 0.58 cells, Dual Contouring 0.01**. That gap is the entire case for dual
contouring, and it is the image the README leads with.

`A-009` · `E-104` · press <kbd>[</kbd> and <kbd>]</kbd> and watch the corner snap in and out as the grid
crosses an aligned resolution

### [▶ Play `surface_nets_vs_marching_cubes`](site:play.html?demo=surface_nets_vs_marching_cubes)

![The same box meshed by Marching Cubes and Surface Nets side by side](docs/gifs/surface-nets-vs-marching-cubes-box.gif)

The catalog used to bill this as "the triangle counts differ". They do not, and the reason is better than
the claim: both methods are pinned to the number of crossed grid edges, so on every closed field
`F_sn − F_mc = 2χ` **exactly**. The HUD recomputes the identity on every re-mesh; change field and
resolution and watch it hold.

`M-2` · `E-103` · what does differ is visible rather than numeric, and the wireframe is where to look

### [▶ Play `marching_cubes_tunnel`](site:play.html?demo=marching_cubes_tunnel)

![One cell meshed as two separate discs and as a single tunnel](docs/gifs/the-tunnel-meshed-as-a-tunnel.gif)

One cell, meshed twice. The face rule alone gives two separate discs; adding the interior rule gives one
cylinder passing through the cell, which is what the trilinear interpolant actually does there. A tunnel is
a handle and a handle costs exactly two of χ, so the claim is arithmetic: **χ falls by two per tunnel and by
nothing else**.

`M-222` · `E-213` · the gold ring is the six body saddles Grosso's construction is built out of

### [▶ Play `game_dig`](site:play.html?demo=game_dig)

![A first-person camera carving tunnels through chunked terrain](docs/gifs/digging-a-tunnel.gif)

The first demo where the mesh is rebuilt **while you are holding the mouse down**. A brush changes an SDF
everywhere, because an SDF is global; what it changes *visibly* is a shell, and the number the whole
incremental story rests on is how thin that shell is. **E1 — the fraction of the brush's own bounding box
that actually moved — was measured at 15–36% offline; here it is on screen, per edit.**

`M-33` · `G-002` · `E-202` · <kbd>C</kbd> outlines exactly the chunks the last edit re-meshed

### [▶ Play `quickstart`](site:play.html?demo=quickstart)

A sphere meshed as eight chunks and put on screen, in one file, with no HUD and nothing to filter out. There
is no animation of it here because there is nothing to animate: the artefact is
[the source](https://github.com/ladvien/isomesh/blob/main/bevy_isomesh/examples/quickstart.rs), which is
108 lines of which 59 are code. What it demonstrates is that the eight chunks are meshed independently and
the result is still watertight.

The one example with no ticket number and no finding, because it teaches nothing — it is the shape of a
working app, and the file to copy.

### [▶ Play `game_mirror_dedup`](site:play.html?demo=game_mirror_dedup)

![The same carved chunk rotated and mirrored, the mirrored copy speckled with markers on every moved vertex](docs/gifs/mirrored-is-not-the-same-mesh.gif)

Of the 48 ways to reorient a chunk on the cubic lattice, exactly **six** give a bit-identical mesh, and they
are precisely the six that never negate a coordinate. So mirroring an asset to reuse it silently misses on
anything keyed to a content hash — instancing, collider caches, network deltas — while rotating it hits.

`✗39 / M-356` · `E-314` · cross-checked against
[`docs/experiments/p-57.csv`](docs/experiments/p-57.csv)

### [▶ Play `game_edit_tape_trim`](site:play.html?demo=game_edit_tape_trim)

![A chunked destructible world, brush gizmos emptying from the frame while the rock does not move](docs/gifs/the-tape-you-keep-is-twenty-times-too-big.gif)

A destructible world is a tape of brush edits folded over a base field, and a Lipschitz bound already cuts a
64-brush tape to a median of 19. **Of the 1,507 brushes that survive that bound world-wide, 1,434 can be
dropped and the mesh does not move by one bit** — a further 20.6×, byte-identical on 64 of 64 chunks.

The first sixty-four frames are the ablation itself, one chunk per frame, because the 1,571 re-meshes it
costs are the point: this is headroom, not a shippable pruner.

`✗41 / M-358` · `E-315` · cross-checked against
[`docs/experiments/p-59.csv`](docs/experiments/p-59.csv)

### [▶ Play `shifted_linear_root`](site:play.html?demo=shifted_linear_root)

![A grid line with standard and shifted reconstructions, the error ratio tracing a closed-form curve as the crossing slides](docs/gifs/where-the-root-falls-decides-the-gain.gif)

Blu, Thévenaz & Unser buy *"about 8 dB asymptotically"* from linear interpolation by shifting the sampling
knots. Does that move the **root** a mesher actually uses? It does, by exactly `|σ − 2τ| / σ` — so the gain
is a lottery over where the crossing falls, and the pre-registered *"at least 30% lower"* is falsified at a
median of 1.486. The honest number is that it clears the bar on 82% of positions.

`✗42 / M-359` · `E-316` · cross-checked against
[`docs/experiments/p-60.csv`](docs/experiments/p-60.csv)

---

## Read

| | |
|---|---|
| **[Findings](FINDINGS.md)** | the ledger: every measurement, every falsified registration, in the order it was made |
| **[The demos](bevy_isomesh/DEMOS.md)** | every runnable example, with an animated capture, the exact command and the finding it demonstrates |
| **[Gameplay](docs/demos/gameplay.md)** | streaming a world past a camera · walking chunk seams with a rigid body · digging · debris that *is* the boolean · LOD ladders |
| **[Algorithms](docs/demos/algorithms.md)** | Marching Cubes · Surface Nets · Dual Contouring · Manifold Dual Contouring · Marching Tetrahedra · greedy quads · a seven-way shootout in one process |
| **[Correctness](docs/demos/correctness.md)** | where a mesh stops being a manifold · what splitting the vertex costs · which way the surface faces · ambiguous faces · the crack between two chunks |
| **[Experiments](docs/experiments.md)** | the ones that held, and what each buys you |
| **[Readme](README.md)** | the repository's own front page, rendered here |

---

## What is not on this site

The **[backlog](BACKLOG.md)**, its **[archive](BACKLOG_ARCHIVE.md)** and the
**[research memos](docs/research)** stay on github.com, and the reason is mechanical rather than editorial:
the archive navigates itself with line-number anchors, which are a feature of GitHub's web view of a file
and have no equivalent in rendered HTML. Every link into them from these pages goes to github.com.
