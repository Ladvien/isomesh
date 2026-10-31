# Layout, running it, and requirements

Moved from the root README on 2026-10-31, with its link paths updated and stale Bevy versions corrected, when the README was brought to the [org documentation guide](https://github.com/bevy-coven/.github/blob/main/DOCUMENTATION.md).

## Layout

```
crates/isomesh/      core. no_std + alloc, one dependency, arrays in every public signature
bevy_isomesh/        Bevy 0.20 bridge and ALL examples. Its own workspace and lockfile.
docs/research/       the papers-derived research this is built on
FINDINGS.md          what we believe and on what evidence
BACKLOG.md           the work queue, and the state
```

`bevy_isomesh` is excluded from the root workspace deliberately: with a shared one, Cargo's feature unification hands glam the `std`, `serde`, `bytemuck` and `encase` features Bevy asks for, and `cargo test` stops testing what a consumer actually gets.

The examples live in `bevy_isomesh` and CI compiles them on every push. That is not tidiness — it is the one thing that stops them rotting. `block-mesh`'s examples are pinned to Bevy 0.13 and `fast-surface-nets`' to Bevy 0.7, both because they lived somewhere nothing ever built.

## Running it

```bash
cargo test -p isomesh                    # the suite prints its own count; doctests included
cargo tree -p isomesh -e normal          # exactly two packages: isomesh, libm

cd bevy_isomesh
cargo run --example quickstart --release                        # start here: an SDF on screen, nothing else
cargo run --example sdf_authoring --release                     # building a field: primitives, operators, one asset
cargo run --example marching_cubes_sphere --release             # the first GIF
cargo run --example surface_nets_vs_marching_cubes --release    # the second and third
cargo run --example dual_contouring_cube --release              # the sharp-feature comparison
cargo run --example sharp_features --release                    # lambda, and what it costs at both ends
cargo run --example qef_clamp --release                         # the clamp, and the red it does not remove
cargo run --example precision_f32_vs_f64 --release               # why CAD needs f64, in two laws
cargo run --example hermite_debug --release                      # what the QEF actually sees
cargo run --example manifold_check --release                    # the red marks, and A to make them go
cargo run --example normal_estimation --release                 # three identical meshes, three shadings
cargo run --example marching_tetrahedra --release               # 3x the triangles, and what they buy
cargo run --example greedy_quads --release                      # 5014 quads down to 1089
cargo run --example transvoxel_seams --release                  # the LOD crack, and T to close it
cargo run --example subgrid_features --release                   # letters thinner than a voxel
cargo run --example game_terrain_stream --release                # a world streaming past you
cargo run --example game_walk --release                          # the acid test: walk every seam
cargo run --example game_capsule_walk --release                  # the same seams, with a body that slides
cargo run --example game_destruction --release                   # shoot it, and the debris is the boolean
cargo run --example game_lod_flyover --release                   # LOD ladder, and the crack count as you fly
cargo run --example game_showcase --release                      # caves, arches, and a roof over your head
cargo run --example game_budget --release                        # 288 chunks re-meshed without missing a frame
cargo run --example game_editor --release                        # undo as a re-fold, and the log's order audited
cargo run --example game_csg_props --release                     # a concave edge, moving, measured every frame
cargo run --example game_paint --release                         # graffiti, then a hole through it, drift 0.000000
cargo run --example resolution_plot --release                    # the fit, and where the model stops describing it
cargo run --example game_dig --release                          # carve, and watch the chunk count
cargo run --example chunk_seam_weld --release                   # the seam, and welding it
cargo run --example marching_cubes_ambiguity --release          # the decider, and how rarely it fires
cargo run --example marching_cubes_interior --release            # the saddle's hyperbola, and the 12.6% MC33 gets wrong
cargo run --example marching_cubes_tunnel --release              # two discs become one cylinder; components 2->1, chi 2->0
cargo run --example gpu_compute_mc --release                     # GPU and CPU Marching Cubes, vertex for vertex
cargo run --example gpu_vs_cpu --release                         # one GPU extraction, five timed parts
cargo run --example gpu_mesh_shader --release                    # field to pixels; four bytes come home
```

**Always `--release`.** A debug build meshes 37–62× slower — both ends measured here, not folklore (FINDINGS M-152; the resolution-sweep incident in `crates/isomesh/Cargo.toml`) — and will convince you something is wrong with the algorithm.

Keys: `W` wireframe · `N` normals · `G` grid · `[` `]` resolution · `1`–`5` field · `S` smoothing · `F12` screenshot · `Esc` quit. Drag to orbit, scroll to zoom.

**Every image in this README and on the demo pages was produced from a command line, and can be
regenerated from one.** No screenshot here was framed by hand, which is what stops the pictures drifting
from the code that made them.

A frame sequence, which is where the GIFs come from:

```bash
ISOMESH_VIEW=nohud ISOMESH_CAPTURE=/tmp/frames ISOMESH_CAPTURE_FRAMES=56 \
  cargo run --example subgrid_features --release
ffmpeg -framerate 12 -i /tmp/frames/frame_%04d.png \
  -vf "scale=780:-1,split[a][b];[a]palettegen=max_colors=128:reserve_transparent=0[p];[b][p]paletteuse" \
  docs/media/subgrid-letters-thinner-than-a-voxel.gif
```

An example driving its own parameter sweep reads `Capture::taken` rather than the clock, so a sequence
is reproducible frame for frame instead of depending on how fast the machine ran.

A single still, which is where the numbers come from — the HUD *is* the evidence, so stills keep it and
GIFs drop it with `nohud`:

```bash
ISOMESH_ALGORITHM=sn ISOMESH_FIELD=5 ISOMESH_SAMPLES=19 \
  ISOMESH_SCREENSHOT=../docs/media/e111-manifold-check-gyroid-surface-nets.png \
  cargo run --example manifold_check --release
```

`ISOMESH_SCREENSHOT` takes one shot and exits. `ISOMESH_FIELD`, `ISOMESH_SAMPLES`, `ISOMESH_VIEW`
(`wire`, `normals`, `nogrid`, `nohud`, `hud`), `ISOMESH_ALGORITHM`, `ISOMESH_CLAMP`, `ISOMESH_OFFSET` and `ISOMESH_WELD`
set what it is a shot *of*. `hud` is the positive form and exists for one demo: `game_dig` opens with the
panel hidden because it is a game first, so its still is taken with `ISOMESH_VIEW=hud`.

## Requirements

Rust **1.89** (edition 2024) for the root workspace, checked in CI against the declared MSRV. The Bevy bridge declares **1.97.1** and pins **Bevy 0.20**, which pins wgpu 29.0.3, glam 0.32 and encase 0.12; those move together or not at all, because Cargo will silently resolve two wgpu majors side by side and the failure only surfaces later as `expected TextureFormat, found a different TextureFormat`.

Developed and measured on the two machines in [`FINDINGS.md`](../../FINDINGS.md)'s header — an Apple M5 (macOS/arm64/Metal) and a Ryzen 9 5900X with an RTX 3090 (Linux/x86-64/Vulkan) — with one resolution-sweep CSV per machine in `docs/measurements/`. CI runs Linux and macOS, which is what makes the bit-reproducibility claim checkable rather than asserted.
