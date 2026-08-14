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

## Run it

Three examples, none of which needs a window — this crate's whole API takes `&wgpu::Device`, so a
demo that required an engine would be demonstrating the engine.

```bash
cargo run -p isomesh-gpu --example extract_a_sphere --release   # the whole loop, checked against the CPU
cargo run -p isomesh-gpu --example gpu_vs_cpu --release         # where the GPU starts winning, and by how much
cargo run -p isomesh-gpu --example mesh_shader_probe            # what this adapter says about mesh shaders
```

**`--release` matters for the second one.** A debug-build CPU extraction is 20–50× slower, which
would flatter the GPU by roughly the factor the example exists to measure.

`gpu_vs_cpu` on an RTX 3090 / Vulkan, sphere, field evaluation included on both sides:

| n | triangles | cpu ms | gpu+read ms | gpu ms | vs cpu | no-read |
|---|---|---|---|---|---|---|
| 17 | 1,064 | 0.056 | 0.290 | 0.172 | **5.2× slower** | 3.1× slower |
| 33 | 4,280 | 0.342 | 0.326 | 0.174 | 1.1× | 2.0× |
| 49 | 9,512 | 1.011 | 0.410 | 0.223 | 2.5× | 4.5× |
| 65 | 17,192 | 2.269 | 0.516 | 0.274 | 4.4× | 8.3× |
| 97 | 38,456 | 6.954 | 0.831 | 0.361 | 8.4× | 19.3× |
| 129 | 68,648 | 16.471 | 1.315 | **0.547** | 12.5× | **30.1×** |

Two things there matter more than the ratio. **The GPU loses below about 25³** — the top row is real
and the example prints it rather than starting the table where the story improves. And the last
column is **nearly flat**: 0.172 → 0.547 ms across a **420× increase in cells**, because the samples
are produced where they are read and the bus is never touched. The extractor is not remotely
saturated at 129³.

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

## Mesh shaders

```bash
cargo run -p isomesh-gpu --example mesh_shader_probe
```

Reports every adapter's `EXPERIMENTAL_MESH_SHADER` bits, then the three things the bits do not say.
Measured here (RTX 3090 / Vulkan): advertised, multiview, points.

**Reaching them needs no `unsafe` from this repository.** Enabling the feature does require a token
whose constructor is a `const unsafe fn` — but this crate never opens a device, only borrows one, and
Bevy writes that `unsafe` itself (`bevy_render` 0.19, `renderer/mod.rs:335`). Bevy's default
`Functionality` priority requests every advertised feature, so its device already reports
`mesh_shader=true` here. These crates stay 100% safe Rust.

wgpu's own source settles the Metal question: the feature reaches Vulkan, DX12 and Metal, but *"naga
is only supported on vulkan; on other platforms you will have to use passthrough shaders."* So mesh
shaders are a fork in the shader pipeline, not a flag on it. **Metal is still unmeasured** — run the
probe there.
