# isomesh

> ⚠️ **Vibe Coded.** This crate was written by an AI agent working from a research corpus and a ticket
> queue. Every number below is produced by a test in this repository and every algorithm cites its
> source, but it has not been through human code review. Read the tests before trusting it with anything
> that matters.

[![crates.io](https://img.shields.io/crates/v/isomesh.svg)](https://crates.io/crates/isomesh) [![docs.rs](https://img.shields.io/docsrs/isomesh)](https://docs.rs/isomesh) [![CI](https://github.com/ladvien/isomesh/actions/workflows/ci.yml/badge.svg)](https://github.com/ladvien/isomesh/actions/workflows/ci.yml) [![license](https://img.shields.io/badge/license-MIT%20OR%20Apache--2.0-blue.svg)](https://github.com/ladvien/isomesh/blob/main/LICENSE-MIT)

**Engine-agnostic isosurface extraction in Rust. Signed distance field in, triangles out.**

`isomesh` has to serve both a real-time voxel game and a CAD tool. That single constraint decides almost
everything about it: no math library appears in a public signature, output buffers are caller-provided
and reusable, the scalar type is generic over `f32` and `f64`, and the crate has exactly one dependency
(`libm`, for `sqrt`/`floor`/`sin`/`cos`, which `core` does not provide on stable).

`no_std` + `alloc`, unconditionally.

![A ball walking across streamed terrain, chunks loading continuously around it](https://raw.githubusercontent.com/ladvien/isomesh/main/docs/gifs/walking-the-seams.gif)

*Nothing there is pre-baked. Every chunk under that ball was extracted while the camera flew toward it,
on a background thread, under a frame budget — and the ball is standing on the **triangles**, not on the
field they came from.*

That last part is the test that matters. Chunks are meshed independently, so whether two of them
actually *meet* is decided by the overlap the chunk layout chose; get it wrong and a player falls
through the world at a boundary. Measured over **495 seam crossings: 0 holes**, worst vertical step
across a seam **0.412 cells** against **0.539 cells** within a single chunk — the joins are smoother than
the terrain they join.

## And one thing nothing else here can do

A feature thinner than a voxel does not exist to a method that asks *"what sign is this grid corner"* —
one bit per edge, and a thin sheet fits between the samples. **Subgrid Marching Tetrahedra** asks
instead for *every zero along the edge*, and triangulates whatever comes back.

On letters **0.35 voxels thick**, Marching Cubes returns **0** triangles and subgrid returns **1,340** —
on the same grid, at any resolution you like.

## Extracting a surface

```rust
use isomesh::marching_cubes::MarchingCubes;
use isomesh::{MeshBuffer, RuntimeShape3, Sdf};

struct Sphere;
impl Sdf for Sphere {
    type Scalar = f32;
    fn sample(&self, p: [f32; 3]) -> f32 {
        (p[0] * p[0] + p[1] * p[1] + p[2] * p[2]).sqrt() - 1.0
    }
}

let shape = RuntimeShape3::new([33; 3])?;
let mut mesh = MeshBuffer::<f32>::new();

// The buffer is caller-provided and reusable: `reset()` clears it without
// releasing capacity, because a real workload re-meshes thousands of chunks
// per edit and an allocation per chunk is the whole budget.
MarchingCubes::<f32>::new().extract(&Sphere, &shape, [-2.0; 3], 0.125, &mut mesh)?;

assert!(mesh.triangle_count() > 0);
# Ok::<(), isomesh::Error>(())
```

Positions and indices come out as `[f32; 3]` and `u32` — no `Vec3` from any particular math library,
because a consumer using two crates that each pin a different `glam` compiles two incompatible `Vec3`
types and the failure surfaces a long way from the cause.

## What is in it

Seven extraction algorithms: **Marching Cubes** (with Marching Cubes 33's asymptotic decider),
**Marching Tetrahedra**, **Surface Nets**, **Dual Contouring**, **Manifold Dual Contouring**, **greedy
quads**, and **subgrid Marching Tetrahedra** — the last of which resolves features thinner than a voxel,
which no sign-based method can do at any grid resolution.

Around them: a mesh validity harness (Euler characteristic, manifoldness, orientation, self-intersection
counts), an accuracy harness, Hermite data and QEF vertex placement, vertex welding, chunk coordinates
with dirty-set re-meshing, brush operations, field-derived LOD, Transvoxel transition cells, collider
readiness checks, frame-budget scheduling, and chunk streaming with hysteresis.

## Choosing an extractor

The honest tradeoff table, from this repo's own measurements (`docs/measurements/shootout.csv` and the
demo pages). "Tiles" means two independently meshed chunks meet with zero boundary edges on the seam —
a structural property of where each method places vertices, not something a future release fixes.

| Extractor | Sharp corners | Tiles across chunks | The tradeoff |
|---|---|---|---|
| `MarchingCubes` | rounded | **yes** (measured 0 seam edges) | the baseline; MC33's asymptotic decider available; emits slivers near zero-valued corners |
| `MarchingTetrahedra` | rounded | unmeasured here | 2.87–3.91× the triangles, measured — not the folklore constant — and more accurate on sharp fields |
| `SurfaceNets` | rounded (0.58 cells off a true corner at 27³) | **no** — gapped, structurally | smoother output, optional smoothing passes |
| `DualContouring` | **held** (0.01 cells at 27³) | **no** — gapped, structurally | QEF with Tikhonov λ and a cell clamp; self-intersections are a measured, non-zero rate |
| `ManifoldDualContouring` | **held** | unmeasured here | one vertex per surface component: takes the non-manifold counts to zero where `DualContouring` cannot |
| `GreedyQuads` | n/a — blocky | open at boundaries by design | Minecraft surface; quad merge measured 1.70×–256× savings depending on the field |
| `SubgridMarchingTetrahedra` | rounded | **yes** (measured 0 seam edges) | resolves features thinner than a voxel, which no sign-based method can; ~70× classic MT (M-98) |

For a chunked world the seam column decides: use `MarchingCubes` or `SubgridMarchingTetrahedra`.

## Verification

Every algorithm ships with a validity gate chosen by the field rather than a blanket rule, a determinism
check, golden hashes that are bit-identical across macOS and Linux, and property tests. Nothing claims a
performance number without a committed benchmark that produced it.

The repository keeps a [`FINDINGS.md`](https://github.com/ladvien/isomesh/blob/main/FINDINGS.md) recording what is known and how well — measured here, verified from
a primary source, reported, or folklore — including the published figures that failed verification.
Falsified entries are never deleted, because which *sources* to distrust is worth more than the
individual fact.

## Troubleshooting

- **It is slow.** You are in a debug build. A debug build meshes **37–62× slower** — both ends measured
  here, not folklore (FINDINGS M-152) — and will convince you something is wrong with the algorithm
  rather than with the profile. `--release`, always.
- **Zero triangles.** Negative is inside, and a sample of exactly zero counts as outside. If your field
  returns plain distance instead of *signed* distance, or the sampled box never crosses the surface,
  every cell classifies the same way and the extractor correctly emits nothing. Check the sign at a
  point you know is inside, and that `origin` and `cell_size` actually span the surface.

## More

The demos are worth looking at before the API:

- **[Gameplay](https://github.com/ladvien/isomesh/blob/main/docs/demos/gameplay.md)** — streaming a
  world past a camera, walking every chunk seam, digging tunnels, LOD cracks and the transition cells
  that close them
- **[Algorithms](https://github.com/ladvien/isomesh/blob/main/docs/demos/algorithms.md)** — the
  extractors side by side, in one process, on the same grids
- **[Correctness](https://github.com/ladvien/isomesh/blob/main/docs/demos/correctness.md)** — where a
  mesh stops being a manifold, and what each repair costs

Bevy integration lives in
[`bevy_isomesh`](https://github.com/ladvien/isomesh/tree/main/bevy_isomesh), in its own workspace so
that Cargo's feature unification cannot leak Bevy's feature choices into this crate's lockfile.

Full README: <https://github.com/ladvien/isomesh>

## Requirements

Rust 1.89, edition 2024.

## License

MIT OR Apache-2.0.
