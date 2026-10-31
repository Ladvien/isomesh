# isomesh

> ⚠️ **Vibe Coded.** Written by AI agents working from a human's direction, covered by tests rather than a line-by-line human review. Read the code before you rely on it for something you care about.

Signed distance field in, triangles out: isosurface extraction for Rust voxel games and CAD tools, with a Bevy plugin.

[![crates.io](https://img.shields.io/crates/v/isomesh.svg)](https://crates.io/crates/isomesh) [![docs.rs](https://img.shields.io/docsrs/isomesh)](https://docs.rs/isomesh) [![CI](https://github.com/bevy-coven/isomesh/actions/workflows/ci.yml/badge.svg)](https://github.com/bevy-coven/isomesh/actions/workflows/ci.yml) [![Bevy 0.20](https://img.shields.io/badge/bevy-0.20-purple.svg)](https://bevyengine.org) [![license](https://img.shields.io/badge/license-MIT%20OR%20Apache--2.0-blue.svg)](https://github.com/bevy-coven/isomesh/blob/main/LICENSE-MIT)

![Eight of the examples running at once: a field built from primitives and smooth-unioned into a mushroom, caves and arches, letters thinner than a single voxel, one cell meshed as two discs and as a tunnel, Surface Nets against Marching Cubes, terrain being carved, a ball crossing streamed chunk seams, and a slab blown apart into debris](https://raw.githubusercontent.com/bevy-coven/isomesh/main/docs/media/kitchen-sink.gif)

*Eight examples, all recorded from this repository: a mushroom built by smooth union, caves and arches, letters thinner than a voxel, one cell meshed two ways, Surface Nets against Marching Cubes, a tunnel being carved, 495 chunk-seam crossings with zero holes, and a slab blown apart.*

**Status:** experimental · Bevy 0.20 · MSRV 1.89 (`isomesh`), 1.97.1 (`bevy_isomesh`)

**Why:** it is the one isosurface crate that is `no_std`, generic over `f32` **and** `f64`, extracts sharp features, and keeps every math library out of its public API, so a Bevy game and a CAD tool can share it without fighting over `glam` versions ([the comparison](https://github.com/bevy-coven/isomesh/blob/main/docs/guide/why.md)).

Three crates live here: [`isomesh`](https://github.com/bevy-coven/isomesh/tree/main/crates/isomesh) (the engine-free core, one dependency), [`bevy_isomesh`](https://github.com/bevy-coven/isomesh/tree/main/bevy_isomesh) (the Bevy plugin and every example) and [`isomesh-gpu`](https://github.com/bevy-coven/isomesh/tree/main/crates/isomesh-gpu) (GPU extraction). An SDF (signed distance field) is a function that returns how far a point is from a surface, negative inside.

## Quick start

```sh
cargo add isomesh
```

Implement `Sdf` for your field, or use one of the fourteen shipped reference fields, and extract:

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

The lines starting with `#` are doctest plumbing: this block is the crate README's example, compiled on every `cargo test`, and `scripts/readme_sync.sh` holds this copy identical to it in CI. For Bevy, `cargo add bevy_isomesh` and see [its README](https://github.com/bevy-coven/isomesh/blob/main/bevy_isomesh/README.md).

## Try it

```sh
cd bevy_isomesh && cargo run --example game_showcase --release
```

It flies a camera through caves and arches that are meshed while it moves. Always `--release`: a debug build meshes 37–62× slower. [Every example's command](https://github.com/bevy-coven/isomesh/blob/main/docs/guide/running.md) is in the guide.

![A sphere meshed by Marching Cubes, with the HUD reporting its triangle count and validity](https://raw.githubusercontent.com/bevy-coven/isomesh/main/docs/media/e101-marching-cubes-sphere.png)

*A still from `marching_cubes_sphere`: every picture here is regenerated from a command line ([how](https://github.com/bevy-coven/isomesh/blob/main/docs/guide/running.md)).*

## What you can do

- **Mesh a chunked voxel world that streams, digs and never cracks at a seam**: `MarchingCubes` and `SubgridMarchingTetrahedra` tile across chunk seams with zero open edges, measured over 495 seam crossings: [what it can do](https://github.com/bevy-coven/isomesh/blob/main/docs/guide/showcase.md).
- **Pick from seven extractors**: Marching Cubes (with MC33's tunnels), Surface Nets, Dual Contouring, Manifold Dual Contouring, Marching Tetrahedra, greedy quads, and subgrid Marching Tetrahedra, which sees features thinner than a voxel: [the algorithms](https://github.com/bevy-coven/isomesh/blob/main/docs/demos/algorithms.md).
- **Build fields from primitives and operators**: smooth union, difference and fourteen reference fields, re-meshed as you edit: [building a field](https://github.com/bevy-coven/isomesh/blob/main/docs/guide/showcase.md#building-the-field-not-meshing-it).
- **Check a mesh before you trust it**: Euler characteristic, manifoldness, orientation, self-intersections, determinism and Hausdorff distance are public API: [what gets checked](https://github.com/bevy-coven/isomesh/blob/main/docs/guide/checks.md).
- **Extract on the GPU**: above about 33³ the GPU path wins, by 37× at 129³ when the field is evaluated in the shader: [gameplay](https://github.com/bevy-coven/isomesh/blob/main/docs/demos/gameplay.md).
- **Run it in Bevy off the main thread**: `bevy_isomesh` meshes on the task pool under a frame budget, with LOD and Transvoxel seams: [the demos](https://github.com/bevy-coven/isomesh/blob/main/docs/guide/demos.md).
- **Know what is and is not done**: what works, what does not yet, and the preregistered experiments behind the numbers: [status](https://github.com/bevy-coven/isomesh/blob/main/docs/guide/fit.md), [experiments](https://github.com/bevy-coven/isomesh/blob/main/docs/guide/experiments.md).

## Part of the Coven

`isomesh` is one of the [Bevy Coven](https://github.com/bevy-coven) crates: the family table there says what each crate gives a game. Grimoire, the Coven's editor, takes `bevy_isomesh` for caves, overhangs and carving: see [use it in Grimoire](https://bevy-coven.github.io/crates/isomesh.html#use-it-in-grimoire).

## Bevy compatibility

| `bevy_isomesh` | Bevy |
|---|---|
| 0.0.13 | 0.20 |

`isomesh` and `isomesh-gpu` do not depend on Bevy.

## Learn more

- [docs.rs/isomesh](https://docs.rs/isomesh) and [docs.rs/bevy_isomesh](https://docs.rs/bevy_isomesh): the API.
- [The book page](https://bevy-coven.github.io/crates/isomesh.html), and [the isomesh site](https://bevy-coven.github.io/isomesh/) with eleven demos playable in the browser.
- [`docs/guide/`](https://github.com/bevy-coven/isomesh/blob/main/docs/guide/README.md): everything this README used to carry, in eight pages.
- [The demo page](https://github.com/bevy-coven/isomesh/blob/main/bevy_isomesh/DEMOS.md): all 54 examples, each with a capture, its command and its finding.
- [`FINDINGS.md`](https://github.com/bevy-coven/isomesh/blob/main/FINDINGS.md): what is believed, how strongly, and on what evidence; [`CHANGELOG.md`](https://github.com/bevy-coven/isomesh/blob/main/CHANGELOG.md).

## Who and licence

Maintained by [@Ladvien](https://github.com/Ladvien). Issues and pull requests are welcome at [github.com/bevy-coven/isomesh/issues](https://github.com/bevy-coven/isomesh/issues); read the org's [CONTRIBUTING.md](https://github.com/bevy-coven/.github/blob/main/CONTRIBUTING.md) first.

MIT OR Apache-2.0, at your option.
