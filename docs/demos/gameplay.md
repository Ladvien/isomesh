# Gameplay demos

The examples that use the crate the way a game does: chunked, edited, streamed, budgeted and collided
against. These are the ones that decide whether the crate is usable, because a meshing library that
produces a beautiful sphere and falls apart at a chunk seam is not one you can ship.

Every figure here came from a command you can run, and the command is under each demo.

[← back to the README](../../README.md)

---

## A world that streams past you

![Terrain streaming past a flying camera](../gifs/terrain-streaming.gif)

*Chunks loading and unloading continuously as the camera flies. Nothing here is pre-baked — every chunk
in frame was extracted while you were watching.*

![Unbroken fBm terrain to the horizon, with a HUD reading 234 chunks resident, 0 waiting, 60 fps](../screenshots/e201-terrain-stream.png)

*`game_terrain_stream` — 234 chunks resident, 117,792 triangles, 3.2 MB, **60 fps at 16.65 ms/frame** while the camera flies and chunks load and unload continuously.*

This is the first example where none of the pieces are visible on their own. **G-007** decides which chunks exist, with a hysteresis band so a camera drifting across the boundary does not re-mesh the same chunk every frame. **B-003** extracts them on the async task pool — never in a system — and applies finished meshes under a frame budget. **G-001**'s layout is what makes each chunk's world position exact rather than merely close, which is what stops the seams (M-32).

The number to watch is not the triangle count. It is **ms/frame while chunks are landing**, because a streaming world that hitches is one doing its meshing on the main thread.

One correction is worth repeating, because it is the mistake the API invites. A radius-based residency rule loads a **ball** of chunks, and a heightfield does not need one: the first version held 952 chunks with 606 permanently waiting, and rendered as holes that never filled. Bounding the vertical extent to the two layers that can contain the surface — which is what a real game does — takes it to 234 resident and nothing waiting (M-104).

```bash
cd bevy_isomesh && cargo run --example game_terrain_stream --release
```

`Space` fly/pause · `[` `]` view distance · `W` wireframe.

---

---

## Walking every seam

![A walker crossing chunk seams on streamed terrain](../gifs/walking-the-seams.gif)

*The orange marker is the walker; every frame casts a dense transect of rays down against the meshed
triangles and counts holes and lips. It crosses a chunk boundary roughly every second.*

![Terrain with a walker on it and a HUD reading 495 seam crossings, 0 holes, seam lip 0.412 cells against terrain roughness 0.539](../screenshots/e203-game-walk.png)

*`game_walk` — **495 seam crossings tested, 0 holes.** The worst vertical discontinuity at a seam is **0.412 cells**, against **0.539 cells** within a single chunk: the joins are smoother than the terrain they join.*

This example is designed to fail. Chunks are meshed independently, and whether two of them actually *meet* is decided by the overlap G-001 chose — get it wrong and you fall through the world at a boundary. So every frame casts a dense transect of rays straight down against the **meshed triangles**, through `parry3d`, and counts holes and lips.

Two details make the answer trustworthy. The ray hits the mesh, not the field — asking the field would test the field, which was never in doubt. And a lip is compared against the terrain's own roughness rather than against zero, because real terrain has real steps and a fixed threshold would be measuring the landscape.

The first version of that test reported **439 holes** and declared the overlap broken. The bug was one operator wide in the test itself: a probe must only count as a hole once *every* chunk layer that could hold the surface has meshed, and the guard said `||` where it needed `&&` (M-105).

```bash
cd bevy_isomesh && cargo run --example game_walk --release
```

`Space` walk/pause · `[` `]` view distance · `W` wireframe.

---

---

## Digging, with the numbers a game actually cares about

![Carving a tunnel into terrain, with the re-meshed chunks outlined](../screenshots/e202-game-dig.png)

*`game_dig` — first person, left click to carve. The blue boxes are the chunks the **last edit** re-meshed: 3 of them, in `0.41 ms`. Nine chunks are resident; the other six were not touched and were not looked at.*

This is the first example where the mesh is rebuilt while someone is holding the mouse down, and it exists to put two numbers on screen that no benchmark can produce:

- **E1 — `265 of 1,728 cells in the brush's bounding box actually re-mesh, 15.3%.`** That is the number the entire incremental story rests on. If it were 100%, being clever about which cells changed would buy nothing over re-meshing the whole box.
- **The trap next to it: `756 cells moved a sample.`** Counting *value* changes rather than *output* changes reads 43% and says incremental meshing is barely worth it; counting output says 15%. The ratio here is `2.85×`, and it was measured offline at 2.8–3.7× before anyone drove it with a mouse.

Edits compose rather than mutate — the field is a stack of brushes over the terrain, which is what makes undo a re-fold of the log rather than a snapshot. So every field sample walks every brush, and the cost grows. Measured over a scripted 60-carve run, median ms per re-meshed chunk:

| edits in the log | 1–15 | 16–30 | 31–45 | 46–60 |
|---|---|---|---|---|
| ms per chunk | 0.158 | 0.354 | 0.525 | 0.589 |

**3.7× for 7× the log, and flattening** — real, and not proportional, which is weaker than "every sample walks every brush" makes it sound.

```bash
cd bevy_isomesh && cargo run --example game_dig --release
```

`LMB` carve · `RMB` fill · `WASD`/`QE` move · `[` `]` radius · `X` clear the log · `C` chunk outlines.

---

---

## Two levels of detail, and the crack between them

![A gyroid meshed at two resolutions with a jagged gap down the seam](../screenshots/e107-transvoxel-seam-cracked.png)

![The same pair with transition cells bridging the seam in orange](../screenshots/e107-transvoxel-seam-stitched.png)

*`transvoxel_seams` — one field, two blocks, the left at full resolution and the right at half. Meshed independently they do not meet: **184 unmatched boundary edges** lie in the seam plane and you can see straight through them. Transition cells take that to **0**, and the orange band is the 310 triangles doing it.*

Both counts are taken **in the seam plane only** — each block is legitimately open at its outer borders, so a global boundary count would drown the signal.

The stitch has a **width**, and that is not cosmetic. A zero-width transition patch also closes the crack — Lengyel says so explicitly — and closing it is *all* it does: every one of its vertices lies in the seam plane, so the patch stands **exactly** perpendicular to the surface it is stitching and shades as a hard crease. Measured: `|cos|` against the surface normal is `0.000` at zero width and `1.000` at `w = 2^(k−2)`. Giving it a width means the coarse block's boundary cells have to be scaled inward by the same amount, which is Lengyel's Equation 4.2 — and written in the block's own cells rather than level-0 cells, the level index cancels out of it entirely.

The property underneath all of it is bit-exactness: a crossing on a half-resolution edge lands on **precisely** the vertex the coarse neighbour's own Marching Cubes pass produced, at every spacing tried including `4/14`. An earlier version of that arithmetic was off by `1.11e-16`, which is a crack no weld can close — a weld merges vertices it can see are the same, and those two are not.

```bash
cd bevy_isomesh && cargo run --example transvoxel_seams --release
```

`T` transitions · `1`–`4` field · `[` `]` resolution.

---

---

## Handing a mesh to a physics engine

`parry3d`'s constructor is not a validity check. Its only documented failure is an empty index buffer — measured, it accepts a zero-area triangle and it accepts a two-chunk mesh with an unwelded seam. A renderer draws that seam correctly; a physics engine reads it as a hole and a character walks through the floor.

So `collider::readiness` reads the validator's report through a collider's eyes and says which of three different things you have:

| | means |
|---|---|
| `is_usable()` | parry will take it and behave — no trailing index, no out-of-range index, no non-finite position |
| `is_seam_free()` | no duplicate vertices, so nothing was assembled from chunks without welding |
| `supports_inside_outside()` | closed, manifold and consistently oriented, so parry's `ORIENTED` pseudo-normals mean something |

They are three predicates rather than one because **a single chunk of a streamed world is open by construction** and is still a perfectly good collider — it just cannot answer "is this point inside". Folding them together would make every chunk in a real world read as broken.

The seam, measured: two adjacent chunks of a torus concatenated give **36 duplicate vertices and 180 boundary edges**. After welding, **0 and 108** — and those 108 are the slab's own outer border, which is genuinely open. The weld closed exactly the 72 that were the seam.

The crate takes no dependency on parry to do this. It emits the `Vec<[u32; 3]>` parry already wants, and the conversion is one line at the call site:

```rust
let indices = isomesh::collider::triangle_indices(&mesh);
let vertices = mesh.positions.iter().map(|p| Vector::new(p[0], p[1], p[2])).collect();
let trimesh = TriMesh::new(vertices, indices)?;
```

---

---

[← back to the README](../../README.md)
