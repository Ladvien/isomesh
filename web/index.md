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

## Play in your browser

Three of the demos are real WebAssembly builds of the same examples the repository ships — no video, no
replay. Each one **measures its own figures in your browser** and then prints a cross-check against the CSV
it was registered with; open the developer console and you can read the comparison, line by line, on your
own hardware. A demo whose numbers stopped agreeing with its committed artefact would say so on screen.

Each module is 36 MB, about 8.4 MB on the wire once the host has gzipped it, and it is cached after the
first visit. Each also carries its own copy of the Bevy runtime, which is where almost all of that goes.

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
