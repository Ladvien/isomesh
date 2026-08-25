# The demos

Thirty-four runnable examples. This page shows what each one looks like, what it proves, and the exact
line to run it.

**[The hosted version of this page](https://isomesh.ladvien.com/demos.html)** carries the same clips in one
link, if you would rather send someone a URL than a repository. Nine of the demos are playable in your
browser there as WebAssembly builds, and the three Phase 21 ones each still print their cross-check against
their committed CSV to the console:
[`quickstart`](https://isomesh.ladvien.com/play.html?demo=quickstart),
[`marching_cubes_tunnel`](https://isomesh.ladvien.com/play.html?demo=marching_cubes_tunnel),
[`dual_contouring_cube`](https://isomesh.ladvien.com/play.html?demo=dual_contouring_cube),
[`surface_nets_vs_marching_cubes`](https://isomesh.ladvien.com/play.html?demo=surface_nets_vs_marching_cubes),
[`game_dig`](https://isomesh.ladvien.com/play.html?demo=game_dig),
[`game_showcase`](https://isomesh.ladvien.com/play.html?demo=game_showcase),
[`game_mirror_dedup`](https://isomesh.ladvien.com/play.html?demo=game_mirror_dedup),
[`game_edit_tape_trim`](https://isomesh.ladvien.com/play.html?demo=game_edit_tape_trim),
[`shifted_linear_root`](https://isomesh.ladvien.com/play.html?demo=shifted_linear_root).

The front page carries a tenth that is not a Bevy build at all: `isomesh` itself compiled to WebAssembly
with a hand-written WebGL2 renderer, **130 KB against these modules' 36 MB**, with the validity report
recomputed live. That is the size the library costs on its own.

Everything here is `cargo run --example <name> --release` from this directory. **Use `--release`** —
debug meshing is 37–62× slower and you will think something is broken.

```bash
git clone https://github.com/ladvien/isomesh
cd isomesh/bevy_isomesh
cargo run --example game_showcase --release
```

## The keys, once

Every example except `quickstart` shares one harness, so these work everywhere:

| | |
|---|---|
| **left-drag** orbit · **scroll** zoom | |
| `W` wireframe | `N` normals | 
| `G` grid and domain box | `Space` pause |
| `1`–`7` switch reference field | `R` re-mesh |
| `F12` screenshot to the working directory | `Esc` quit |

Per-example keys are listed with each demo. The HUD always shows vertices, triangles, extraction
milliseconds, frame time and FPS.

---

## See it work

### A world with a roof over your head

![Flying through a landscape riddled with caves, arches and tunnels](https://raw.githubusercontent.com/ladvien/isomesh/main/docs/gifs/flying-through-the-rock.gif)

A camera flying **through** rock — under arches, into tunnels, out the far side. The field is nine
lines and nothing in it is authored:

```text
solid(p) = max( p.y − height(x, z),  |gyroid(p)| − thickness )
```

A `max` is an intersection, so rock exists only where a point is **below the terrain** *and* **inside a
thickened gyroid**. The gyroid is triply periodic, so it tunnels in `x`, `y` and `z` by construction.

A heightfield stores one number per column and cannot represent any of it. That is the whole reason to
reach for a voxel mesher.

```bash
cargo run --example game_showcase --release     # Space pause · [ ] fly speed · 1-3 how much cave
```

### A world that streams past you

![Endless fBm terrain streaming past a flying camera, chunks appearing at the frontier](https://raw.githubusercontent.com/ladvien/isomesh/main/docs/gifs/terrain-streaming.gif)

`fbm_terrain` sampled without bound. It is a function, so there is no edge to reach — fly long enough
and every chunk you can see was meshed while you were watching. The number to watch is `ms/frame`
*while chunks are landing*, not after.

```bash
cargo run --example game_terrain_stream --release   # Space pause · [ ] view distance
```

### Walking every seam

![A ball walking across streamed terrain, chunks loading continuously around it](https://raw.githubusercontent.com/ladvien/isomesh/main/docs/gifs/walking-the-seams.gif)

The acid test, and it is designed to fail: walk every chunk seam, ray-cast down at every step, count
what happens. **495 crossings, zero fall-throughs** (M-106).

```bash
cargo run --example game_walk --release         # Space pause the walk · [ ] view distance
```

### Building the field, not meshing it

![Four SDF primitives on a shelf and a mushroom assembled from them by smooth union and difference](https://raw.githubusercontent.com/ladvien/isomesh/main/docs/gifs/building-a-field.gif)

The only example with no meshing content at all — the extractor is the default and the resolution is
fixed. The only thing that changes is the expression.

Left, four primitives. Right,
`Union { SmoothUnion { stem, Difference { cap, flat }, k }, gills }`. **`k` is the knob a level
designer reaches for**: at zero the stem meets the cap in a crease, and by 0.25 it is a fillet. Watch
the junction.

```bash
cargo run --example sdf_authoring --release     # [ ] sweep k · 1-4 isolate a primitive · 0 assembly
```

### The shortest path from a field to a mesh

`quickstart` is one file, no HUD, no keys, nothing to filter out. A sphere meshed as eight chunks and
put on screen. Copy it.

```bash
cargo run --example quickstart --release
```

---

## Algorithms, side by side

### Proving a cell empty when the ball cannot (M-354)

![A gyroid with amber cages on cells only the affine bound rejects](https://raw.githubusercontent.com/ladvien/isomesh/main/docs/gifs/affine-rejects-where-lipschitz-cannot.gif)

Hart's Lipschitz test answers "can the surface pass through this cell?" with a **ball**, and a ball throws
away everything but one number. A **revised affine form** — five `f64`, fixed size, no heap — bounds the
whole expression over the cell's **box** and keeps the correlation between terms.

On `gyroid`, whose three trig terms cannot peak at once, that is **688 → 2,647** cells rejected, **3.85×**,
with **zero** cells lost. On `gyroid_uncapped` the Lipschitz test rejects **nothing at all** and the affine
form rejects 1,098.

Then cycle to `box_exact` and every amber cage vanishes: the ratio is exactly `1.000000` and the two
rejected **sets are identical, cell for cell**, because `min`/`max` collapses the form to an interval.
**A mechanism that wins everywhere is usually a bug; this one wins exactly where its story says.**

```bash
cargo run --example affine_rejection --release
```

### Where the root falls decides the whole gain (M-359)

![A grid line with standard and shifted reconstructions, the error ratio tracing a closed-form curve as the crossing slides](https://raw.githubusercontent.com/ladvien/isomesh/main/docs/gifs/where-the-root-falls-decides-the-gain.gif)

Blu, Thévenaz & Unser get *"about 8 dB asymptotically"* out of standard linear interpolation by shifting
the sampling knots by a fixed `τ` and enforcing the interpolation property. At `τ = 1/5` the prefilter is
`cₙ = −2⁻²·cₙ₋₁ + (1 + 2⁻²)·fₙ` — **multiplication-free**. So does it move the *root* a mesher uses?

**It does, by exactly `|σ − 2τ| / σ`.** Slide the crossing across the cell and watch the measured ratio
trace the closed form: exactly **1** when the root lands on a knot, exactly **0** at `σ = 2τ = 0.4`, never
worse than the standard rule in the quadratic regime, clearing a 30% bar on **82%** of positions. A mesher
cannot choose `σ`, so the pre-registered *"at least 30% lower"* is **falsified at a median of 1.486** and
the honest number is the 82%, not the 8 dB.

What *does* ship is the guard band: the prefilter is causal with `(1/4)^k` decay, and **10 preceding
samples** buys a root within `1e-6` cells of the whole-line answer on 8 of 8 fields, where 5 samples
fails. Two controls are on screen — a synthetic quadratic whose predicted ratio is exactly `1/5` and
measures `0.1984`, and a sabotage run that applies the shift *without* the prefilter and comes back **39×
worse**, biased by exactly `τ`.

```bash
cargo run --example shifted_linear_root --release
```

### The spheres a mesh never touches (M-355)

![Wireframe spheres around SDF samples, most never reached by the mesh](https://raw.githubusercontent.com/ladvien/isomesh/main/docs/gifs/untouched-tangency-spheres.gif)

A distance sample asserts a sphere the surface must **touch**. Nothing in this crate checks that. Between
**3.13% and 80.47%** of samples have a sphere no vertex ever reaches — and on `thin_plate` switching
Marching Cubes to Dual Contouring drops it from 717 to 45 per 1,000, a measured **15.9984×**.

The worst misses are exact constants: `box_exact` under Marching Cubes misses by **`√2` cells**, because
MC places no vertex on a box edge at all. Needs no ground truth — a reference-free detail-loss signal.

```bash
cargo run --example untouched_spheres --release
```

### The repair that welds 520 things shut (M-352)

![A bonsai CT scan covered in magenta pinch markers](https://raw.githubusercontent.com/ladvien/isomesh/main/docs/gifs/pinch-repair-welds-520-components.gif)

Integer CT data puts samples exactly on the isovalue, and **every** degenerate triangle in the mesh comes
from one — 58,097 of 58,097 on `bonsai`. The repair removes all of them and **moves no geometry at all**.

It is still only safe on one of the two volumes. On `bonsai` it welds **520 separate components**, because
516 of 17,201 collapse groups join vertices that share no triangle. **Ship the precondition, not the
repair** — one union-find decides it before anything changes.

```bash
cargo run --example pinch_repair --release
```

### An error bound that is reached, not merely respected (M-350)

![Normal ray pairs fanning at a CSG crease, with a meter tracking the predicted bound](https://raw.githubusercontent.com/ladvien/isomesh/main/docs/gifs/seam-normal-bound-is-attained.gif)

A central difference across a CSG seam averages two gradients, so its error is at most **half the crease
angle** — provable, because the difference returns `u − Λw` for a diagonal `Λ` with entries `≤ ½`.

It is **attained**: median tightness **0.9748** over 24 rows. At 175° the worst error is 2.498° against a
2.5° bound. Everywhere else normals are exact to `6.75e-10°`. The defect is a curve one cell wide, and it
does not shrink with resolution.

```bash
cargo run --example seam_normal_bound --release
```

### The corner, and the 101× (M-54)

![Surface Nets rounding a box corner beside Dual Contouring holding it, across a resolution sweep](https://raw.githubusercontent.com/ladvien/isomesh/main/docs/gifs/dual-contouring-vs-surface-nets.gif)

Same field, same grid, same edge crossings. The only difference between these two meshes is one
function: **where a cell's vertex goes.** Surface Nets takes the centroid of the crossings; Dual
Contouring solves for the point that best fits the crossing planes.

On a sharp field that is worth **101×** in Hausdorff distance — `7.217e-2` against `7.145e-4`. On a
smooth one it is worth 1.2×, which is why this is a choice and not an upgrade.

The sweep deliberately steps *around* grid-aligned resolutions: on an aligned grid the comparison
inverts, for a reason no viewer could guess from the picture.

```bash
cargo run --example dual_contouring_cube --release    # [ ] resolution · C cell clamp · 1-2 field
```

### Marching Cubes, refining

![A sphere meshed by Marching Cubes as resolution sweeps up and back down](https://raw.githubusercontent.com/ladvien/isomesh/main/docs/gifs/marching-cubes-sphere-resolution-sweep.gif)

The baseline, and the first thing in this project you can look at. The HUD recomputes the Euler
characteristic and manifoldness on every re-mesh, so you can watch χ stay at 2 while the mesh changes
underneath it.

```bash
cargo run --example marching_cubes_sphere --release   # [ ] resolution, 5 to 129 samples
```

### Surface Nets against Marching Cubes

![The same box meshed by Marching Cubes and Surface Nets side by side](https://raw.githubusercontent.com/ladvien/isomesh/main/docs/gifs/surface-nets-vs-marching-cubes-box.gif)

Marching Cubes grey on the left, Surface Nets tan on the right. The HUD computes `F_sn − F_mc` live and
compares it against `2χ`, so you watch ✗1's identity hold at every resolution:

```text
V_sn = V_mc + χ        F_sn = F_mc + 2χ
```

Exactly — not approximately. Folklore says Surface Nets produces "substantially fewer triangles". It
produces `2χ` fewer, which on a sphere is four.

![The same comparison on a capped gyroid, a high-genus field](https://raw.githubusercontent.com/ladvien/isomesh/main/docs/gifs/surface-nets-vs-marching-cubes-gyroid.gif)

On the gyroid the same identity holds with χ far from 2, which is the version worth trusting.

```bash
cargo run --example surface_nets_vs_marching_cubes --release   # 1-5 field · [ ] resolution · S smoothing
```

### Letters thinner than a voxel

![The word ISO meshed by two extractors as the letters thin; the Marching Cubes panel loses them entirely](https://raw.githubusercontent.com/ladvien/isomesh/main/docs/gifs/subgrid-letters-thinner-than-a-voxel.gif)

One field, one grid, two extractors, and a sweep driving the letters from 1.6 voxels thick down to 0.2.

The left panel is Marching Cubes. As thickness drops below one voxel it disintegrates — first a
scatter, then nothing. The right panel is subgrid Marching Tetrahedra and is unchanged.

**Marching Cubes returns 0 triangles where subgrid returns 1,340.** It is also the most expensive
extractor in the crate by a wide margin (M-98: 196× Marching Cubes), which is the trade.

```bash
cargo run --example subgrid_features --release        # - = thickness · [ ] resolution
```

### The rest

| example | what it shows |
|---|---|
| `marching_tetrahedra` | Marching Cubes against Marching Tetrahedra, wireframe **on** by default — the finding is the triangle count, 2.87–3.91× |
| `greedy_quads` | The blocky path: one quad per face against greedy merging of coplanar runs. 1.70–4.60× saving by field, against a widely-quoted 2.76× from one benchmark |
| `normal_estimation` | Three normal strategies, geometrically identical meshes, glossy material — the demo is entirely in the speculars |
| `resolution_plot` | The `t = a + b·n³` fit drawn with its residuals. The two-term model describes Marching Cubes and **does not** describe Surface Nets |

---

## Game-shaped

### The rim of your own tunnel is lit wrong (M-350)

![A shaft bored into rock, its rim glowing red with normal error](https://raw.githubusercontent.com/ladvien/isomesh/main/docs/gifs/carve-seam-lit-wrong.gif)

Dig a tunnel and the rim is a CSG crease. A central difference across it averages two gradients, so the
normal can be **44.3° wrong** at a right-angled rim — a bright arc exactly where the player just dug.

Then the punchline: **30.2× the triangles and the rim got 0.66° worse.** Refining narrows the wrong band
without making it less wrong, because the stencil step scales with position and not with the cell. The fix
is the analytic gradient on **38 of 1,374 vertices** — 2.8%, since straddling vertices lie on a curve — and
that takes the rim from 44.28° to **0.0012°** at the *cheapest* resolution.

```bash
cargo run --example game_carve_seams --release
```

### The chunk knows it is losing detail (M-355)

![Terrain of spires and boulders under a chunk heat map, cycling LOD strategies](https://raw.githubusercontent.com/ladvien/isomesh/main/docs/gifs/detail-oracle-half-the-triangles.gif)

Camera distance does not know what is *in* a chunk. At 10 px/cell it spends 33³ on plain rolling ground and
drops the chunks holding thin spires to 17³ — **24,076 triangles and 9 of 16 chunks still losing detail**,
worse on the metric than uniform 9³ at 8.8× the cost.

Counting the tangency spheres the mesh never touches answers the question from the field alone, no ground
truth. Each chunk climbs its own ladder to a 15% gate: **81,276 triangles, 0 chunks losing detail**, against
**147,504** for the range policy tuned to the same bar — **1.81× fewer**. A bake, not a frame: 875 ms on
the worst rung.

```bash
cargo run --example game_lod_oracle --release
```

### The tape you keep is twenty times too big (M-358)

![A chunked destructible world, brush gizmos emptying from the frame while the rock does not move](https://raw.githubusercontent.com/ladvien/isomesh/main/docs/gifs/the-tape-you-keep-is-twenty-times-too-big.gif)

A destructible world is a **tape** of brush edits folded over a base field. A Lipschitz interval bound
deletes the brushes that provably cannot win in a chunk, taking a 64-brush tape to a median of **19**
survivors with the mesh byte-identical (M-341). This asks how many of the 19 are *needed*.

**1,434 of 1,507 surviving brushes — 95.16% — can be dropped and the mesh does not move by one bit.** Not
one at a time, either: re-mesh each chunk from its necessary brushes *alone* and the hash is identical on
**64 of 64** chunks, so the world's tape cuts from 1,507 to **73**. A further **20.6×** on top of the
bound. Four hashes on screen, three tapes, one rock that never changes shape.

The cost is on screen too, because it has to be: deciding necessity took **1,571 re-meshes**. This is
headroom, not a shippable pruner.

```bash
cargo run --example game_edit_tape_trim --release
```

### Mirroring an asset breaks mesh-hash dedup; rotating it does not (M-356)

![The same carved chunk rotated and mirrored, the mirrored copy speckled with markers on every moved vertex](https://raw.githubusercontent.com/ladvien/isomesh/main/docs/gifs/mirrored-is-not-the-same-mesh.gif)

Of the **48** ways to reorient a chunk on the cubic lattice, exactly **six** give a bit-identical mesh —
and they are precisely the six that never negate a coordinate. Negate any axis and the mesh moves, because
every grid edge is oriented along *increasing* index: a sign flip makes the extractor compute `b/(b−a)`
where it computed `a/(a−b)`. Same point, different bits.

So mirroring an asset to reuse it — the oldest trick in level art — silently misses on anything keyed to a
content hash: instancing, collision-mesh caching, chunk mesh caches, network deltas. Rotating it hits.
**And the count is predicted from the field before any extractor runs**: cut edges whose crossing differs
when computed from the far endpoint — 72 on `sphere`, 643 on `noise_cavity` — equals the number of moved
vertices *exactly*. `box_exact` is the control at **0 of 1,350**, and it is the one field where mirroring
is bit-exact.

```bash
cargo run --example game_mirror_dedup --release
```

### Digging, the way a game does it

![A first-person camera carving tunnels through chunked terrain](https://raw.githubusercontent.com/ladvien/isomesh/main/docs/gifs/digging-a-tunnel.gif)

`WASD` to move, mouse to look, left click to carve, right click to fill. `C` outlines exactly the
chunks the last edit re-meshed.

The number this example exists for is **E1**: a brush changes **15–36%** of the cells inside its own
bounding box (M-33). Counting value changes overstates the re-mesh set by 2.8–3.7× (M-34).

```bash
cargo run --example game_dig --release          # WASD Q E move · click carve · right click fill · [ ] radius · C outlines
```

### Sixty-four edits deep, and the mesher only looks at twenty-one

![A 64-brush sculpting stroke where the pruned per-chunk cost climbs at a third of the full tape's rate](https://raw.githubusercontent.com/ladvien/isomesh/main/docs/gifs/tape-pruning-the-tape-shrinks.gif)

Every sample of every chunk evaluates every brush, so an edit history is a tax the mesher pays forever.
Bound each brush over the chunk box — one sample per brush, `f(c) ± l·r` — and drop the ones that
provably cannot win the `min`/`max` chain there.

**Median 21 of 64 brushes survive; 6.46 ms per chunk against 15.84 ms; 2.45× on the world; 25× on the
best chunk.** And the mesh does not move — IEEE `min`/`max` *select* an operand rather than computing
one, so dropping a loser costs zero ULP, and the HUD checks all 64 chunks bit-exact every sweep (M-341).

`P` turns pruning off so you can watch the cost jump while nothing else changes. The demo also refuses
to oversell itself: a *uniform* tape keeps a constant surviving fraction, so the green line climbs too —
2.5× shallower, not flat. The flat line belongs to a moving stroke.

```bash
cargo run --example tape_pruning --release   # P prune · H heat · C chunk boxes · X restart
```

### The debris is the boolean

![A hollow shell being shot, cratering, with debris falling away](https://raw.githubusercontent.com/ladvien/isomesh/main/docs/gifs/the-debris-is-the-boolean.gif)

A shot appends `Brush::subtract(sphere)` to crater the target, *and* meshes the **intersection** of the
pre-shot solid with that same sphere as the debris. Nothing is authored: the crater and the fragment
are two views of one boolean.

It fires by itself, every 0.9 seconds, up to 24 shots.

```bash
cargo run --example game_destruction --release  # Space fire · 1-3 target (wall/shell/spiral) · [ ] charge
```

### Paint that survives the wall it was sprayed on

![Graffiti sprayed on a wall, then a hole blown through it, the surviving paint exactly where it was](https://raw.githubusercontent.com/ladvien/isomesh/main/docs/gifs/paint-that-survives-the-wall.gif)

Spray, then blow a hole through the wall. The paint on what remains is exactly where you put it — not
smeared, not reset. The drift readout is **continuously zero**, because paint lives in the edit log and
was never on the surface.

```bash
cargo run --example game_paint --release        # click spray · right click carve · 1-5 colour · [ ] nozzle
```

### Undo without a snapshot

![A CSG solid morphing backwards and forwards as an edit log cursor moves](https://raw.githubusercontent.com/ladvien/isomesh/main/docs/gifs/undo-is-a-refold.gif)

The edits are a **log**, and the field is a fold of that log over a base. Undo moves the cursor back by
one and re-folds. Nothing is stored and nothing is copied.

The order is load-bearing, and the crate knows exactly how much: same-kind hard edits commute
bit-identically across all 8! = 40,320 orderings (M-36), mixed add-and-subtract gives **11** distinct
results (M-37), and smooth union gives **40,317** (M-38). Press `S` to swap the last two ops and see
whether `commutes_with` predicted it.

```bash
cargo run --example game_editor --release       # Z undo · Y redo · E edit · S swap last two · X clear
```

### A concave edge, moving, measured every frame

![A CSG solid whose cutter orbits, re-meshed every frame](https://raw.githubusercontent.com/ladvien/isomesh/main/docs/gifs/a-boolean-remeshed-every-frame.gif)

A **concave** edge is where a CAD tool lives, and it is the case where the vertex solve wants to sit
outside the cell that produced it. Here it also moves: the cutter orbits on a slow ellipse, sweeping
the edge across cell boundaries continuously.

The HUD reports the worst error **over the whole sweep**, because a single-position measurement of a
sharp feature is a measurement of that position.

```bash
cargo run --example game_csg_props --release    # Space pause · A extractor · [ ] resolution
```

### Two levels of detail, and the crack between them

![A camera flying along an LOD ladder, blocks changing level as it passes](https://raw.githubusercontent.com/ladvien/isomesh/main/docs/gifs/lod-flyover.gif)

Slabs at descending resolutions, and a camera flying out and back so coarse blocks end up on both
sides. `T` toggles transition cells; the HUD counts seam-plane boundary edges **per side** and reports
the worst vertex displacement at each level switch — the pop size, in cells.

```bash
cargo run --example game_lod_flyover --release  # Space pause · T transitions · [ ] speed · R reset
```

### The rest

| example | what it shows |
|---|---|
| `game_budget` | The same work spread across frames, and what the guarantee costs. Deliberately overloaded: every chunk dirty at once, queue refilling as it drains. `U` drains unbudgeted as the control |
| `game_capsule_walk` | An `avian3d` rigid body driven across the same terrain, comparing commanded distance against travelled — what the seams stole. A ray cannot answer "can a body move through here" |

---

## Where it breaks, and what that proves

### Where a mesh stops being a manifold

![Non-manifold edges drawn in red on a gyroid mesh, appearing and vanishing as resolution changes](https://raw.githubusercontent.com/ladvien/isomesh/main/docs/gifs/manifold-check-resolution.gif)

Non-manifold edges as thick red lines, non-manifold vertices as red spheres, boundary edges in amber —
drawn **in place on the mesh**, not counted in a corner.

Cycle `A` from Surface Nets to Dual Contouring to Manifold Dual Contouring on the gyroid and the red
marks disappear. Then read [`docs/experiments.md`](../docs/experiments.md#the-paper-is-wrong-and-here-is-the-fixture),
because on `noise_cavity` they do not, and the paper says they should.

Surface Nets is non-manifold where one cell carries two sheets: 48 non-manifold edges on the capped
gyroid, 15 on `fbm_terrain` (M-4). It is a **resolution** effect, not a topology one (M-15).

```bash
cargo run --example manifold_check --release    # A algorithm · 1-7 field · B boundary · [ ] resolution
```

### Counting the defects before the mesh exists

![noise_cavity at 65 cubed, cyan cages on critical cells with one yellow non-manifold vertex inside each](https://raw.githubusercontent.com/ladvien/isomesh/main/docs/gifs/critical-cells-count-the-defects.gif)

Cyan cages are **2D-critical** cells, magenta are **3D-critical**, yellow dots are the non-manifold
vertices Dual Contouring produced. Look at the cluster: one dot per cage, no strays, no empty cages.

That is M-338. A cell is critical when its eight corner signs host one of Latecki's two configurations,
and on `noise_cavity` at 65³ the census reads **567 + 35 = 602** against **602** non-manifold vertices
and **602** critical cells hosting one — 100% co-location, against a chance baseline near 1%. Cycle to
`sphere` or `csg_difference` and every number is zero, from both directions.

The census is a function of the sign bytes alone, so it is available **before** extraction: 13.5 ms
against 63 ms to mesh the same grid. The 256-entry classification is enumerated from the definitions at
startup — the example logs `120 2D-critical, 8 3D-critical, 0 in both` before it opens a window, so a
transcription error could not survive the first line.

```bash
cargo run --example critical_cells --release   # 1-5 field · F fly to the cluster · H surface
```

### An ambiguous face, and how rarely one turns up

![Cells with ambiguous faces boxed in amber and magenta, changing as resolution steps](https://raw.githubusercontent.com/ladvien/isomesh/main/docs/gifs/ambiguous-faces-are-rare.gif)

Plain Marching Cubes on the left, the asymptotic decider on the right, and every cell with an ambiguous
face boxed: **amber** where the decider agreed, **magenta** where it disagreed and joined the corners.

Magenta is rare, and the rarity is the finding. On five of the eight reference fields the ambiguous
face **never occurs at all** — so Marching Cubes 33 and plain Marching Cubes are bit-identical on them
at every resolution tested (M-40). Only `gyroid` (0.515% of surface cells) and `fbm_terrain` (1.532%)
reach it.

```bash
cargo run --example marching_cubes_ambiguity --release   # 1-3 field · A cell markers · [ ] resolution
```

### One cell, meshed twice

![One cell meshed as two separate discs and as a single tunnel](https://raw.githubusercontent.com/ladvien/isomesh/main/docs/gifs/the-tunnel-meshed-as-a-tunnel.gif)

Left, the face rule alone: two separate discs. Right, the same cell with the interior rule: one
cylinder passing through. Plus the closed gold ring of six body saddles — the inner hexagon — that
decides which it is.

χ falls by exactly two per tunnel and by nothing else, because a handle costs a closed surface exactly
two and everything else the rule does is topology-neutral (M-222).

**One configuration is refused rather than meshed.** A cell whose contours run past Grosso's Corollary
6 bound has no published triangulation, so `extract` returns `Error::UnresolvedSixSaddle` instead of
emitting a hole.

```bash
cargo run --example marching_cubes_tunnel --release      # 1-2 configuration · H hexagon · C contours
```

### The bilinear saddle, swept

![A plane sweeping through one cell while a dot traces a hyperbola to infinity and back](https://raw.githubusercontent.com/ladvien/isomesh/main/docs/gifs/the-interior-decider-sweep.gif)

No mesh at all — one wire cell, a plane sweeping bottom to top, and a white dot tracking the bilinear
saddle. The dot's trail draws a hyperbola: it runs off to infinity and returns from the other side at
the magenta pole plane.

Two verdicts are printed side by side: Chernyaev's numerator-only test, and the corrected one. On
configuration `1` they disagree, and **12.6%** of the time the classic test is the one that is wrong.

```bash
cargo run --example marching_cubes_interior --release    # Space pause · [ ] scrub by hand · 1-2 config
```

### The sharpness knob, at both ends

![A box with a sphere bitten out of it, its edges rounding over and snapping back as lambda sweeps](https://raw.githubusercontent.com/ladvien/isomesh/main/docs/gifs/sharp-features-lambda-sweep.gif)

What trades sharpness against stability is **λ, the Tikhonov regularizer** in the vertex solve — a
number that was a compile-time constant until this example needed to turn it.

Toward λ→0 corners come out exactly, and flat cells fling their vertices out: M-30 measured one landing
**3.18 cells** outside the cell that produced it. Toward large λ every edge rounds over into Surface
Nets. Watch both rims.

```bash
cargo run --example sharp_features --release    # - = lambda · C cell clamp · 1-4 field
```

### What the clamp removes, and the half it cannot reach

![A gyroid with self-intersecting triangles highlighted in red, appearing and vanishing as the clamp toggles](https://raw.githubusercontent.com/ladvien/isomesh/main/docs/gifs/qef-clamp-self-intersections.gif)

`I` draws self-intersecting triangles in red. Toggle the clamp and on five of the eight fields they go
to **exactly zero**.

On `gyroid` and `fbm_terrain` they do not — 3.12 and 13.84 pairs per 1,000 triangles survive. What is
left is a **connectivity** failure rather than a placement one, which is why pressing `A` for Manifold
Dual Contouring makes it *worse* (3.118 → 5.669), not better.

```bash
cargo run --example qef_clamp --release         # C clamp · A algorithm · I highlight · 1-4 field
```

### Where `f32` blurs, and where it tears

![The same field at a large coordinate offset in f32 and f64, blurring and then tearing open](https://raw.githubusercontent.com/ladvien/isomesh/main/docs/gifs/precision-f32-tears.gif)

Two failures, two laws, and neither is where you would guess. At 10⁶ `f32` does not crack — what moves
is accuracy. The crack is an order of magnitude further out, and at **2²³** it tears: χ drops from 2 to
1, the vertex count collapses, and real holes with boundary edges open.

```bash
cargo run --example precision_f32_vs_f64 --release   # - = offset exponent · 1-3 field · [ ] resolution
```

### The rest

| example | what it shows |
|---|---|
| `chunk_seam_weld` | Two chunks meshed independently — as a real game does, since an edit only dirties what it touches. Red lines are the crack; `V` welds it shut, `E` slides the chunks apart to prove they are separate meshes |
| `transvoxel_seams` | Two LOD blocks meeting on a plane, with and without transition cells. `T` toggles. Reproduce the committed stills with `ISOMESH_FIELD=3 ISOMESH_SAMPLES=6 ISOMESH_TRANSITIONS=0` |
| `hermite_debug` | What Dual Contouring actually operates on: one cell blown up, its eight corners, its edge crossings, the normal at each, where the QEF put the vertex — and in red, where it was before the clamp dragged it back |
| `weld_creases` | Two cubes, same input, same tolerance. The left is welded on position alone and its eight corners swallow all 24 vertices, so the flat shading that made it read as a cube is gone; the right is welded with a key built from the normals and keeps all 24. **Both welds are correct** — only the caller knows which was wanted. The readout also shows the right mesh gaining 24 boundary edges, which is the split seen from the edge column rather than damage |
| `sealed_cave` | **Did I just seal the cave?** Two chambers in *different chunks*, joined by a tunnel through the chunk between them; `F` plugs it, `G` re-opens it. The mesh answers none of this — a mesh can be closed, manifold and Hausdorff-close while sealing a passage the field leaves open — so the readout is `connectivity::AirWorld`: components go `1 → 2`, `A↔B` flips to **SEALED**, and `visited` shows what the repair cost. Watch that against `chunk` rather than `world`: the search runs inside one grid and cannot leave it (M-322), and here it is **0.03× a chunk** because the chambers live elsewhere and the boundary graph resolves the global split (M-323) |

---

## On the GPU

### The whole pipeline, and the bus is never touched

![A field with three brushes moving through it, extracted and drawn entirely on the GPU](https://raw.githubusercontent.com/ladvien/isomesh/main/docs/gifs/gpu-resident-mesh-shader.gif)

Field evaluation, extraction and draw, all on the GPU. Per frame the CPU sends a camera matrix and
three brushes, and waits. There is no vertex buffer, no index buffer, and **zero `Mesh3d` entities** —
nothing goes through Bevy's mesh pipeline at all.

Needs an adapter advertising `EXPERIMENTAL_MESH_SHADER`. Without one the pipeline is not built, the
HUD says so, and nothing is drawn.

```bash
cargo run --example gpu_mesh_shader --release   # [ ] resolution · Space pause brushes · T triangle count
```

### The honest GPU result

| example | what it shows |
|---|---|
| `gpu_compute_mc` | The same Marching Cubes on both paths. Two meshes look identical at any resolution including when one is wrong, so `V` colours the surface by CPU/GPU agreement instead. Triangle counts match; most vertices are bit-identical and the rest differ by one ULP, because WGSL permits FMA contraction |
| `gpu_vs_cpu` | One extraction broken into its five parts, which is the only way to see that **the extraction is not what costs anything** — `count + emit` is 0.045 ms at 129³ and barely moves across a 420× rise in cells. Everything the GPU path costs is data movement around it |

### So is the GPU faster?

**Yes, above about 33 samples per axis — and by 37× at 129³.** Below that it is not, and the reason is
a fixed cost of roughly 0.22 ms that does not care how big the grid is.

| samples per axis | CPU, single thread | GPU, field evaluated on the GPU | |
|---|---|---|---|
| 17³ | **0.06 ms** | 0.22 ms | CPU ahead 3.7× |
| 33³ | 0.34 ms | **0.23 ms** | GPU ahead 1.5× |
| 65³ | 2.44 ms | **0.27 ms** | GPU ahead 9× |
| 129³ | 20.14 ms | **0.54 ms** | GPU ahead **37×** |

The shape matters more than any single row: the GPU column is **nearly flat** — 0.22 to 0.54 ms across
a 420× rise in cell count — because the extraction itself was never the cost. `count + emit` is
0.045 ms at 129³ and does not move.

**Where the field is evaluated decides everything.** Sample on the CPU and upload, and the upload is
87% of the path — 8.37 ms at 129³, or 2.4× ahead of the CPU instead of 37×. Evaluate it in the shader
and that entire cost disappears, because the samples are produced where they are read. Three tickets
took this path from 15.01 ms to 0.54 ms at 129³ and **none of them made the extractor faster**; every
gain was data movement removed.

Numbers from `docs/measurements/gpu_vs_cpu.csv`, warmed, median of three, RTX 3090 over Vulkan
(M-145, M-150, M-155).

---

## Recording these yourself

Every clip on this page was produced by one command:

```bash
scripts/record_all_gifs.sh                    # all of them
scripts/record_all_gifs.sh dual_ qef_         # only clips whose name matches
```

It drives `ISOMESH_CAPTURE`, which writes numbered PNGs through Bevy's own screenshot path — a GPU
readback, so it works over a window the compositor never mapped and cannot catch another window
passing over the top. The harness also takes `ISOMESH_VIEW=nohud` for clips meant to be looked at
rather than read, and `ISOMESH_SPIN` to rotate a static scene.

Needs an X display. The parameters each example wants are in the script's header, along with the four
that were wrong on the first attempt and why.

---

## More

| | |
|---|---|
| The crate | [`bevy_isomesh` on crates.io](https://crates.io/crates/bevy_isomesh) · [docs.rs](https://docs.rs/bevy_isomesh) |
| The core | [`isomesh`](https://crates.io/crates/isomesh) — engine-agnostic, `no_std`, one dependency |
| Longer write-ups | [`docs/demos/`](../docs/demos/) — algorithms, correctness, gameplay |
| Every measurement | [`FINDINGS.md`](../FINDINGS.md) · [`docs/experiments.md`](../docs/experiments.md) |
