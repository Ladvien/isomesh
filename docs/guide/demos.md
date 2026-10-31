# Demos

Moved from the root README on 2026-10-31, with its link paths updated and stale Bevy versions corrected, when the README was brought to the [org documentation guide](https://github.com/bevy-coven/.github/blob/main/DOCUMENTATION.md).

**[→ The demo page](../../bevy_isomesh/DEMOS.md)** — all 54 examples, each with an animated capture, the
exact command, its controls, and the finding it demonstrates.

**[→ bevy-coven.github.io/isomesh](https://bevy-coven.github.io/isomesh/)** is the same prose as a site, with
every clip served from it rather than hotlinked. Eleven of the demos are **playable in the browser** there
as WebAssembly builds, and the three Phase 21 ones each still print their cross-check against their
committed CSV to the console. The front page carries a tenth module that is not a Bevy build: `isomesh`
itself, with a hand-written WebGL2 renderer and the validity report recomputed live, at **133 KB against
the Bevy modules' 36 MB**.

It is served over **HTTPS**, and that is load-bearing rather than hygiene: WebGPU is a
secure-context-only API, so on an `http://` page `navigator.gpu` is `undefined` and every one of those
eleven modules refuses to start — in Chrome on a desktop exactly as on a phone.

The longer write-ups live on three further pages, so this one stays short. Every figure on them came
from a command you can run, and every number is measured on one of the two machines in
[`FINDINGS.md`](../../FINDINGS.md)'s header rather than quoted from a paper.

| | |
|---|---|
| **[Gameplay](../demos/gameplay.md)** | the whole extract-and-draw pipeline on the GPU · streaming a world past a camera · walking every chunk seam with a rigid body · digging tunnels · shooting a wall so the debris *is* the boolean · spraying graffiti that survives the wall being blown open · flying an LOD ladder and counting what opens up · handing a mesh to a physics engine |
| **[Algorithms](../demos/algorithms.md)** | Marching Cubes · Surface Nets · Dual Contouring · Manifold Dual Contouring · Marching Tetrahedra · greedy quads · subgrid Marching Tetrahedra, and a seven-way shootout in one process |
| **[Correctness](../demos/correctness.md)** | where a mesh stops being a manifold · what splitting the vertex costs · which way the surface faces · ambiguous faces · the crack between two chunks |

Between them they carry every measured figure this crate makes a claim about, and each one names the command that regenerates it.
