# isomesh — implementation brief

**Date:** 2026-08-11
**Audience:** the coding agent working in this repo.
**Companion docs:** `CLAUDE.md` (standing rules), `2026-08-11-bevy-examples-catalog.md` (what to demo),
`docs/research/` (why).

Stages are ordered so each one is falsifiable before the next depends on it. **Do not start stage N+1
until stage N's exit criterion passes.** If a stage can't pass, say so — don't proceed and hope.

---

## Stage 0 — Skeleton and the two traits everything hangs off

**Goal:** the workspace exists, the core crate has exactly one dependency, and the API shape is
committed before any algorithm exists.

### The four types

```rust
/// A scalar field. This is the ONLY thing an algorithm needs to see.
pub trait Sdf {
    type Scalar: Real;
    fn sample(&self, p: [Self::Scalar; 3]) -> Self::Scalar;
    /// Optional analytic gradient. Default: central differences.
    fn gradient(&self, p: [Self::Scalar; 3]) -> [Self::Scalar; 3] { /* central diff */ }
}

/// Linearization of a 3D index space. Deliberately array-based — see ndshape.
pub trait Shape3 {
    fn linearize(&self, p: [u32; 3]) -> u32;
    fn delinearize(&self, i: u32) -> [u32; 3];
    fn size(&self) -> [u32; 3];
}

/// Where triangles go. Implemented by MeshBuffer, by bevy_isomesh's Mesh writer,
/// and by anything a CAD consumer wants.
pub trait MeshSink {
    fn vertex(&mut self, position: [f32; 3], normal: [f32; 3]) -> u32;
    fn triangle(&mut self, a: u32, b: u32, c: u32);
    fn reserve(&mut self, _verts: usize, _tris: usize) {}
}

/// The default sink. Reusable across chunks — this is the point.
pub struct MeshBuffer {
    pub positions: Vec<[f32; 3]>,
    pub normals:   Vec<[f32; 3]>,
    pub indices:   Vec<u32>,
}
impl MeshBuffer {
    /// Truncate to zero WITHOUT releasing capacity.
    pub fn reset(&mut self) { /* .clear() on each, never shrink */ }
}
```

`Real` is a small sealed trait over `f32`/`f64` — enough arithmetic for the solvers, nothing more.
Default `f32`; `f64` is what makes the crate viable for CAD, and `fast-surface-nets` is disqualified
for CAD precisely because its `SignedDistance: Into<f32>` bound forecloses it. Don't repeat that.

### Exit criterion

```bash
cargo tree -p isomesh        # exactly two lines: isomesh, glam
cargo test -p isomesh        # traits compile; MeshBuffer::reset preserves capacity (assert it)
grep -ri bevy crates/        # no output
```

---

## Stage 1 — Marching Cubes, as the reference implementation

Everything else is compared against this, so correctness matters more than speed here.

### Do not type the case table from memory

The 256-entry triangle table is where transcription errors hide, and a wrong entry produces a mesh that
*looks* fine and is silently non-manifold. Two defenses, use both:

1. Take the table from a verified published source and cite the source in a comment.
2. **Write a validator that checks all 256 cases without a reference table.** For each configuration,
   the generated triangles' boundary edges must lie entirely on cube faces, and each face's boundary
   edge set must be consistent with that face's 4 corner signs. This catches transcription errors
   structurally. Run it as a unit test, not a one-off.

### Sign convention — pick one and enforce it

**Negative is inside.** Write it in a doc comment on `Sdf`, and add a debug assertion in the example
fields. Half of all "my mesh is inside out" bugs are a flipped convention crossing a module boundary.

### Vertex placement

Linear interpolation along the edge: for an edge with endpoint values `a`, `b` and positions `pa`,
`pb`, `t = a / (a - b)`, position `pa + t*(pb - pa)`. Guard `|a - b| < eps` → `t = 0.5`.

### Normals

Two paths, both worth having because the examples compare them:

- **Gradient** — central differences on the field at the vertex. Smooth, correct, costs 6 samples.
- **Area-weighted face normals** — accumulate per-triangle normals into vertices. Cheap, and it's what
  most engines actually do.

### Exit criterion

- Euler characteristic `V - E + F == 2` on a closed sphere at 3 resolutions.
- Every edge shared by exactly 2 faces; violation count reported and equal to 0.
- Deterministic: same input twice → identical buffers.
- Example `mc_sphere` runs.

---

## Stage 2 — MC33 / the asymptotic decider

Plain MC has ambiguous faces and produces holes. This is the fix, and it's small.

For a face with corner values `f00, f10, f11, f01`, the bilinear interpolant's saddle value is

```
S = (f00·f11 − f10·f01) / (f00 + f11 − f10 − f01)
```

The sign of `S` decides whether the two positive corners connect across the face or are separated.
Guard the denominator near zero. Interior ambiguity needs the trilinear body saddle as well — the
catalog doc covers which of the 6 ambiguous cases need it.

**Exit criterion:** the `mc33_ambiguity` example shows a field that holes under plain MC and is closed
under MC33, and the Euler test passes on that field *only* with MC33 enabled. A test that passes both
ways isn't testing this.

---

## Stage 3 — Surface Nets

One vertex per cell containing a sign change, placed at the centroid of the edge crossings, then
smoothed. Simpler than MC and it's what most engines actually ship — and, per the speed analysis,
**it has no credible published timings anywhere**, which is what makes stage 7 worth doing.

**Exit criterion:** same three property tests. Plus `surface_nets_sphere` renders next to MC on the
same field with a triangle-count readout — the counts should differ substantially and both meshes
should be closed.

---

## Stage 4 — Dual Contouring, sharp features, and the vertex solve

This is the CAD-differentiating stage and the one with real math in it. Read
`docs/research/2026-08-10-adjacent-math-transfer-audit.md` first.

### Hermite data

DC needs, per edge crossing: the position **and the surface normal there**. Store as
`(position: [Real;3], normal: [Real;3])` per crossing, up to 12 per cell. This is the input to the
vertex solve and the thing that makes sharp features recoverable at all.

### Default vertex rule — closed-form, three planes

Use this before reaching for a general QEF. It is exactly rotation-equivariant, needs no iterative
solve, and avoids squaring the condition number:

```
c  = (p₁+p₂+p₃)/3          dᵢ = nᵢ·(pᵢ−c)
x  = c + [ d₁(n₂×n₃) + d₂(n₃×n₁) + d₃(n₁×n₂) ] / [ n₁·(n₂×n₃) ]
```

Falls back when the triple product is near zero (near-parallel planes) — that's the degenerate case,
and it should route to the regularized form below rather than producing a huge coordinate.

### General case — regularized normal equations

For >3 planes, or when the closed form degenerates:

```
M = Σ nᵢnᵢᵀ        g = Σ dᵢnᵢ        λ ≈ 0.01
x = c + adj(M + λI)·g / det(M + λI)
```

`λ` is the regularizer that keeps under-determined cells (flat regions, where `M` is rank 1) from
flying off. Note that `M = AᵀA` squares the condition number — this is the exact reason the QR/Givens
formulation exists in the literature. In `f32` it will bite; in `f64` it mostly won't. **Measure it**:
the `precision_f32_vs_f64` example exists for this.

### The clamp — do not skip this

**Clamp the solved vertex to (1−ε) inside its own cell**, ε ≈ 1e-4. Then instrument
self-intersections per 1,000 triangles, with and without.

This is item 2 on the opportunities list, and it is the cheapest experiment in the whole project.
Measured context from the research: ODC (2024) reports Manifold DC at **100% of models
self-intersecting** vs ODC at **0 of 1500**. If the clamp gets you most of the way there, guaranteed
intersection-free extraction feeds convex decomposition that can't fail — and CoACD's measured
**49% → 80%** improvement in downstream manipulation success is what's on the other side of that.

**Exit criterion:** `dual_contouring_cube` shows a box with genuinely sharp edges where Surface Nets
rounds them. `qef_clamp` shows the self-intersection counter changing when the clamp is toggled, live.
Both numbers recorded in the commit message.

---

## Stage 5 — Chunking and seams

The real workload is not "mesh one volume." It's "re-mesh the 3% of cells a brush touched."

- Chunks with a **1-cell overlap** on the positive faces so neighbouring chunks agree on shared cells.
- Vertex welding across chunk boundaries, or prove it isn't needed.
- A dirty-set API: `mesh_dirty(&mut self, chunks: &[ChunkId], out: &mut ...)`.

**Exit criterion:** `chunked_terrain` shows a multi-chunk field with **no visible cracks** and a
chunk-boundary overlay toggle. `edit_brush` re-meshes only touched chunks and displays the count.

**Also instrument E1 here** — hash each chunk's cell slab and log what fraction actually changes per
brush stroke. That number is the ceiling on every incremental-repair idea in the opportunities doc, and
nobody has published it. Hours of work, and it's the first thing the research says to measure.

---

## Stage 6 — LOD and Transvoxel transitions

Adjacent chunks at different resolutions crack. Transvoxel's transition cells are the published fix.

**Exit criterion:** `transvoxel_seams` shows two LOD levels adjacent, with a toggle. Cracks visible
when off, gone when on. This is a visual acceptance test — take a screenshot both ways and commit them.

---

## Stage 7 — The shootout

**This is the highest-value stage in the brief and it is mostly plumbing.**

Marching Cubes vs Surface Nets vs Dual Contouring vs Marching Tetrahedra. Same field, same grid, same
hardware, same codebase, same machine, same run. Report per algorithm: extraction ms, triangle count,
vertex count, self-intersections/1k, non-manifold edges, and Hausdorff error vs the analytic surface.

Per `docs/research/2026-08-11-meshing-speed-analysis.md`: **no paper published since 2020 benchmarks
these against each other.** Every "X is faster than Y" claim in circulation traces to different
hardware, different grids, different years. You would have the only apples-to-apples measurement.

Two methodology requirements the research makes non-negotiable:

- **Sweep resolution, don't report one grid size.** Fit `t = a + b·n³`. Report `a` explicitly — in the
  published FlexiCubes numbers, **73% of the 64³ timing is fixed launch overhead**, not meshing. A
  single small-grid number measures dispatch latency.
- **Report the stage breakdown, not just extraction.** Grosso & Zint measured contouring at 68 ms
  against 58 ms of halfedge construction — the contour is **54%** of getting to a usable mesh. Time
  contour / normals / weld / collider / upload separately or the headline number is misleading.

**Exit criterion:** `bench_shootout` prints a table and `bench_resolution_sweep` prints the fitted `a`.
Commit the raw numbers as a CSV in `docs/measurements/`.

---

## Stage 8 — GPU compute path

`isomesh-gpu`, wgpu 29.0.3, public API takes `&wgpu::Device` / `&wgpu::Queue` /
`&mut wgpu::CommandEncoder`. **Build the standalone headless harness first, with no Bevy in the room.**
If it can't run against raw wgpu, the abstraction leaked and the CAD story is already broken.

Shader composition: `include_str!` plus a ~40-line preprocessor for `#include` and a boolean `#ifdef`.
Do **not** add `naga_oil` — it's Bevy-owned, 14 months stale, and Bevy's own 0.19 notes say they're
moving to WESL. Revisit `wesl-rs` when the shader count justifies it.

Then `bevy_isomesh` adds a `RenderGraph`-schedule system:

```rust
render_app.add_systems(RenderGraph, mesh_dirty_chunks.before(camera_driver));
```

**Exit criterion:** `gpu_compute_mc` produces output **bit-identical** to the CPU path on the same
field (or documents exactly where and why it differs — fast-math reassociation is a legitimate answer,
"close enough" is not). `gpu_vs_cpu` shows both live with a timing HUD.

---

## Stage 9 — Mesh shaders (gated, exploratory, do not block on this)

Behind a cargo feature, off by default, with a runtime capability probe and a graceful fallback.

```rust
app.add_plugins(DefaultPlugins.set(RenderPlugin {
    render_creation: WgpuSettings {
        features: WgpuFeatures::EXPERIMENTAL_MESH_SHADER,
        ..default()
    }.into(),
    ..default()
}));
// then: render_device.features().contains(...) before using it
// raw device: render_device.wgpu_device() -> &wgpu::Device
```

Bevy's own docs warn that setting `WgpuSettings.features` *"may cause renderer initialization to
fail"* — hence the probe. **On this machine (macOS/arm64/Metal) support is unverified.** wgpu's spec
table lists MSL as *planned* while the tracking issue says the Metal backend merged. **First task in
this stage is a probe that prints what the adapter actually reports.** Report the result before writing
any shader.

Why it's worth eventually doing: the same Marching Cubes went **114.2 → 2679.4 fps (23.4×)** moving
from a compute shader to a mesh shader, purely because the mesh-shader version never writes vertices
out to VRAM. That is a bigger effect than the entire measured spread between extraction algorithms
(1.5–3.9×). But it's experimental upstream with an open redesign issue, so: gated, probed, optional.

---

## The ordering rationale, in one line

Stages 1–4 build the algorithms. Stage 5 makes them useful. Stage 7 produces the measurement that
doesn't exist in the literature and costs almost nothing once 1–4 are done. Stages 8–9 are speed, and
the speed analysis says speed is the *last* term that matters — so it's last.
