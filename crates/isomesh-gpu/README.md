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
- **`Composer`** — WGSL `#include` and `#ifdef`, over modules compiled in with `include_str!`. Not
  `naga_oil`; see the module docs for why.
- **`FEATURES`** — every compile-time flag any shader here reads. The validation sweep covers the
  cross product of modules with **every subset** of it, so a flag missing from this list is a branch
  nothing ever compiles.

- **`MarchingCubesGpu`** — the compute kernel. Two passes (count, prefix-sum, emit) so the output is
  dense *and* in cell order, which is what makes it comparable with the CPU and with itself.

## The case table is uploaded, not transcribed

The usual way a GPU port breaks "never guess a case table" is by pasting a second copy into WGSL.
There isn't one. `case_table_bytes()` packs `isomesh::marching_cubes::table::CASES` — itself derived
by a `const fn` rather than copied from a paper — and a test unpacks the bytes and compares all 256
entries.

## Measured against the CPU

Same grid, same samples, same table, 33³ sphere:

| | |
|---|---|
| triangle count vs CPU | **equal**, on `sphere`, `torus` and `box_exact` |
| vertices bit-identical to a CPU vertex | **6,507 of 6,936** |
| within one ULP per axis | 429 |
| further than one ULP | **0** |

The 1-ULP miss is float contraction, not the algorithm: WGSL permits a multiply-add to be fused, and
this driver takes it, rounding once where the CPU rounds twice. The test asserts the bound, so two
ULPs would fail it.

Normals are the one deliberate divergence — a shader cannot call `Sdf::gradient`, so it central-
differences the sample grid. `isomesh`'s M-65 measures that at 0.460° worst against the analytic
gradient at 17³, converging at `h²`.

## Validation

Two layers, and they catch different things.

```bash
cargo test -p isomesh-gpu every_shader_permutation_validates   # naga, no GPU, belongs in CI
cargo test -p isomesh-gpu is_valid_wgsl_on_a_real_device       # the driver's own opinion
```

The `naga` sweep needs no adapter and catches "compiles on my Vulkan driver, explodes on DX12". A
companion test feeds the validator invalid WGSL and asserts it *rejects* it — a gate that has only
ever passed is indistinguishable from one that cannot fail.

```bash
cargo test -p isomesh-gpu -- --nocapture   # prints the adapter it ran on
```
