# Correctness demos

Where meshes stop being meshes. Non-manifold edges, inside-out triangles, ambiguous faces, cracks at a
seam — each shown with the count that makes it real rather than described in the abstract.

These are the least glamorous demos and the most load-bearing. A renderer forgives every defect on this
page. A physics engine forgives none of them.

[← back to the README](../../README.md)

---

## Where a mesh stops being a manifold

![Marching Cubes and Surface Nets on the capped gyroid, non-manifold features marked](https://raw.githubusercontent.com/ladvien/isomesh/main/docs/gifs/manifold-check-resolution.gif)

*`manifold_check` — the capped gyroid at 19³ under Surface Nets. Every red sphere is a non-manifold vertex and every red line a non-manifold edge, drawn where the validator found them: 39 edges and 61 vertices, clustered around the tunnel mouths where two sheets of surface share a cell. The same field and grid under Marching Cubes reports `0` on every counter and `MANIFOLD, CLOSED` ([screenshot](../screenshots/e111-manifold-check-gyroid-marching-cubes.png)).*

A count tells you a mesh is broken without telling you where, and the two most useful findings in this project were both about *where*. The marks come from `validate_features`, which returns the offending edges and vertices from the **same pass** that produced the numbers beside them — so the picture and the caption cannot drift apart.

```bash
cd bevy_isomesh && cargo run --example manifold_check --release
```

`1`–`7` field · `A` algorithm · `B` boundary overlay · `[` `]` resolution.

---

---

## Splitting the vertex, and what it costs

![The capped gyroid under Dual Contouring, covered in red non-manifold marks](../screenshots/e111-manifold-check-gyroid-dual-contouring.png)

![The same field and grid under Manifold Dual Contouring, with no marks at all](../screenshots/e111-manifold-check-gyroid-manifold-dual-contouring.png)

*The same field, the same 19³ grid, one algorithm apart. **Dual Contouring: 39 non-manifold edges, 61 non-manifold vertices, `χ = -10`, one component.** **Manifold Dual Contouring: zero, zero, `χ = -2`, seven components** — and the same 3,276 triangles. Press `A` in `manifold_check` to switch between them.*

One vertex per cell is the defect. The fix is one vertex per **surface component**: the cell's cut edges are partitioned into the cycles the Marching Cubes table already links them into, and each cycle gets its own QEF solve. Ju's own paper describes it and credits it to Nielson — the output is the *dual* of the Marching Cubes surface, so it inherits Marching Cubes' topology.

That inheritance is asserted, not assumed. On every closed field at 17³, 25³ and 33³ the dual reproduces Marching Cubes' Euler characteristic **and its component count** exactly. Look again at the two captions: the pinch was not only breaking the index buffer, it was **fusing seven pieces into one and misreporting `χ` by eight**.

Three things this measured that were not what the tickets predicted:

- **The cost is zero on five of the seven fields**, and about **5%** of the run time on the other two. Only `gyroid` and `fbm_terrain` ever need a second vertex in a cell, and their rate *falls* with resolution — 3.13% → 2.05% → 0.53% at 17³/25³/33³. Nielson's published *"about 1.3%"* counts entries in the case table, not cells in a scene. Triangle counts are unchanged: splitting moves vertices without adding quads.
- **Self-intersections get worse, not better** — `gyroid` 3.118 → 5.669 per 1,000 triangles, `fbm_terrain` 13.837 → 15.434. The prediction registered before the run said the opposite. Two vertices in one cell is exactly what breaks the within-cell partition the clamp's guarantee rests on, and a 2024 result reporting Manifold Dual Contouring as 100% self-intersecting was on the record the whole time.
- **A second, unrelated non-manifold mechanism exists.** The dual of a manifold surface is a manifold *complex*, and an index buffer cannot hold two distinct edges between one pair of vertices, so parallel dual edges collapse into one edge with four faces. The property suite found it by shrinking to the exact same three-sphere fixture that falsified unconditional manifoldness for Marching Cubes. It is identified by arithmetic rather than by eye: a collapse costs exactly one edge, so `χ_dual − χ_mc == non_manifold_edges`, measured `1 − 0 == 1` at `h = 2/3` and `0 == 0` at every finer grid.

---

---

## Which way does the surface face

Three answers, and they are not the same answer. Ask the **field** for its gradient; **difference** the field, which is all a sampled voxel buffer can offer; or take the **area-weighted average of the incident triangles**, which is all a mesh can offer. `normals::recompute` re-derives a finished `MeshBuffer` under any of them, so the choice survives welding and merging instead of being baked in by whichever extractor ran.

Differencing at the cell size — the case a game without an analytic field is stuck with — costs under half a degree, and converges the way it must:

| grid | worst | mean |
|---|---|---|
| 17³ | 0.460° | 0.299° |
| 33³ | 0.121° | 0.079° |
| 65³ | 0.031° | 0.020° |

Successive ratios 3.76 and 3.92, so `h²`, asserted as a range rather than admired in a log.

![The same CSG solid shaded three ways: crisp, slightly softened, and smeared](../screenshots/e113-normal-estimation-csg.png)

*`normal_estimation` — a sphere bitten out of a box at 41³. **All three panels are the same 2,244 vertices and 4,484 triangles**; only the normal buffer differs. Left, the field's gradient keeps the bite's rim crisp and its staircase legible. Right, area-weighted normals smear that rim and round the steps into blobs — 46.426° off at the worst. The middle is differencing the field at the cell size, which is what a sampled voxel buffer can offer: 17.974° worst, and 0.450° on average.*

```bash
cd bevy_isomesh && cargo run --example normal_estimation --release
```

`1`–`5` field · `[` `]` resolution · `W` wireframe.

The third strategy is where it gets interesting. Area-weighted normals track the field closely on smooth geometry and **cannot** on sharp geometry, because a corner vertex gets the average of three face normals where the field's gradient gives one of them. On a sphere the mean disagreement falls 3.25° → 2.16° → 1.08° across those grids and on a torus 11.65° → 6.07° → 2.45°. On `box_exact` the *worst* disagreement is **35.796° at all three resolutions, identical to six figures** — refining a grid does not soften a corner. That invariance is the assertion; the constant is just the box's corner.

---

---

## An ambiguous face, and how rarely one turns up

![Marching Cubes beside the asymptotic decider on a capped gyroid, ambiguous cells marked](https://raw.githubusercontent.com/ladvien/isomesh/main/docs/gifs/ambiguous-faces-are-rare.gif)

*`marching_cubes_ambiguity` — plain Marching Cubes on the left, the same extraction with `FaceAmbiguity::AsymptoticDecider` on the right. Every box is a cell with an ambiguous face: **amber where the decider agreed and separated the corners, magenta where it disagreed and joined them.** Magenta is the only place the two meshes can differ, and on the gyroid at 33³ there are nine such cells out of 5,240.*

The catalog originally specified this example as "holes on the left, closed on the right". **That cannot be shown, because it does not happen** — this crate's case table is derived at compile time by walking each face counter-clockwise, so two cells sharing a face cannot disagree about it and neither side ever holes. What the decider changes is *which* surface gets built on an ambiguous face, which the HUD reads off as a difference in Euler characteristic.

![The same comparison on a sphere, where the two meshes are byte-identical](../screenshots/e102-ambiguity-sphere-identical.png)

*Press `3`. On a sphere there is **not one ambiguous face in 1,160 surface cells**, and the two meshes are byte-identical — which the committed golden fixture also pins. Five of the eight reference fields behave this way; only the gyroid (0.515% of cells) and `fbm_terrain` (1.532%) reach the configuration at all. An example that only ever showed the interesting case would misrepresent how often the interesting case arrives.*

```bash
cd bevy_isomesh && cargo run --example marching_cubes_ambiguity --release
```

`1`–`3` field · `A` cell markers · `[` `]` resolution.

---

---

## A crack between two chunks, and welding it shut

![Two chunks of a torus, meshed independently, with the open seam marked in red](../screenshots/e115-chunk-seam-unwelded.png)

*`chunk_seam_weld` — one torus, two chunks, meshed **independently**, exactly as a game does when an edit dirties only the chunks it touches. Every red line is a boundary edge on the shared plane: a triangle with no neighbour. `80` of them, and `40` duplicated vertices. The surface looks continuous and is not.*

![The same two chunks after welding, with no seam](../screenshots/e115-chunk-seam-welded.png)

*The same two chunks after `V`. **`1328 → 1288` vertices, 40 merged, 0 triangles collapsed, and the seam carries no boundary at all.** χ stays `0` — it is a torus either way; what changed is that it is now one surface.*

The spacing selector is the part worth pressing. `1` is `h = 0.125` and `2` is `h = 4/35`, and only one of those is arbitrary: two chunks agree on their shared sample plane bit-for-bit **only when the cell size is a power of two**, because one computes `(o + h·cn) + h·n` and the other `o + h·(c+1)n` — equal by algebra, not by IEEE. 22% of random `(origin, h, cells, chunk)` combinations disagree by an ulp, and `4/35` came out of that search. A weld keyed on exact equality closes the seam under `1` and silently leaves it open under `2`; this one is an epsilon weld for exactly that reason.

```bash
cd bevy_isomesh && cargo run --example chunk_seam_weld --release
```

`V` weld · `E` explode the chunks apart · `1` `2` spacing · `[` `]` resolution.

---

---

[← back to the README](../../README.md)
