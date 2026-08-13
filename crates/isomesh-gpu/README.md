# isomesh-gpu

> ⚠️ **Vibe Coded.** Written by an AI agent from a research corpus and a ticket queue. Every number
> here is produced by a test in this repository, but it has not been through human code review.

GPU isosurface extraction for [`isomesh`](../isomesh), on raw `wgpu`. No engine types anywhere in the
public API.

## The rule

**Every entry point takes `&wgpu::Device`, `&wgpu::Queue` or `&mut wgpu::CommandEncoder`.**

A Bevy consumer reaches the raw device through `RenderDevice::wgpu_device()` and hands it in. A CAD
tool passes its own. A test uses `headless::Gpu`. There is no second entry point that takes a
`RenderDevice`, because that would make the engine a dependency of the algorithm instead of a caller
of it.

## The version pin is load-bearing

`wgpu 29.0.3`, matching Bevy 0.19 — and it has to be exact. Cargo resolves two wgpu majors **side by
side with no resolution error**; the failure arrives much later as `expected TextureFormat, found a
different TextureFormat`. Verified: this workspace's lockfile and `bevy_isomesh`'s independently
resolve `wgpu 29.0.4` and `wgpu-types 29.0.4`, the same patch in both.

| `isomesh-gpu` | `wgpu` | `bevy` |
|---|---|---|
| 0.0.x | 29.0.4 | 0.19 |

## No software fallback

`headless::Gpu::new` asks for a high-performance adapter with `force_fallback_adapter: false`, and a
missing adapter is an error rather than a CPU reference driver. The entire reason this crate exists is
to be compared against the CPU path — a benchmark that silently ran on lavapipe reports numbers orders
of magnitude off and looks merely slow.

So a machine with no GPU **fails** these tests. That is the intended behaviour.

## What is here

```rust
use isomesh::fields::Sphere;
use isomesh_gpu::{FieldBuffer, GridParams, headless, read_buffer};

let gpu = headless::Gpu::new()?;
let grid = GridParams::new([33; 3], [-2.0; 3], 0.125)?;

let field = FieldBuffer::sampled(gpu.device(), gpu.queue(), grid, &Sphere::<f32>::canonical())?;
let back = read_buffer(gpu.device(), gpu.queue(), field.buffer(), grid.field_buffer_size())?;
```

- **`GridParams`** — a validated sampling grid, with the 32-byte two-`vec4` layout a shader reads it
  as. Positions are the index *multiplied*, never accumulated; `isomesh`'s M-70 and M-73 both record
  cracks caused by the other choice.
- **`FieldBuffer`** — `f32` samples in GPU memory, uploaded from a slice or sampled from any
  `isomesh::Sdf`.
- **`read_buffer`** — results back to the CPU, blocking.
- **`headless::Gpu`** — a device with no window, and an `AdapterReport` of what it turned out to be.

Shaders, their composition, and Marching Cubes itself are GPU-002 through GPU-004.

```bash
cargo test -p isomesh-gpu -- --nocapture   # prints the adapter it ran on
```
