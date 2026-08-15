---
name: Bug report
about: Something extracted wrong, crashed, or measured strangely
title: ""
labels: bug
---

**The exact command**

```bash
# e.g. cd bevy_isomesh && cargo run --example game_walk --release
```

**Was it `--release`?**

A debug build meshes 37–62× slower (measured — FINDINGS M-152) and makes healthy algorithms look
broken. If the problem is speed and the build was debug, that is probably the whole story.

**What happened, and what you expected**

If a mesh looks wrong, the validity harness usually says how it is wrong:
`isomesh::validate::validate` reports Euler characteristic, non-manifold counts, orientation and
boundary edges — paste the report if you can.

**Machine**

OS, CPU, and (for anything GPU or timing-related) the GPU and backend — e.g. `RTX 3090 / Vulkan`.
Naming the machine is a house rule here; every number in the repo does it.

**Field**

If the field is yours rather than one of the seven reference fields: does `sample` return *signed*
distance (negative inside)? Sign convention is the most common zero-triangles cause — see the
Troubleshooting section of `crates/isomesh/README.md`.
