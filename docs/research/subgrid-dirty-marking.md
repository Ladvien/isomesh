# Marking the chunks an edit changed, for an extractor that reads more than the corners

**Issue #11** (owner's review, 2026-10-30): `mark_changed_chunks` never marks a Subgrid chunk whose sub-cell geometry an edit adds or removes. It decides "is this a surface cell" from the signs of a cell's eight corners, and `Extractor::Subgrid` does not use corner signs. The chunk keeps its old `ChunkMesh` for good: a removed thin wall stays on screen and keeps its collider, an added one never appears, and nothing reports it (`EditReport::output_changed_cells` is 0, no error is raised).

## Sources

| Kind | Stem | Where |
|---|---|---|
| Algorithm | `sig2026_Subgrid_Marching_Tetrahedra` (Baktash, Gillespie and Crane, *Subgrid Marching Tetrahedra*, ACM TOG 45(4), 2026) | §1 and §1.3: edge intersections, not corner signs, are the input; 1D sampling along an edge is the drop-in way to find them |
| Shipping game | `transvoxel_dissertation_lengyel2010` (E. Lengyel, *Voxel-Based Terrain for Real-Time Virtual Simulations*, UC Davis, 2010) | §3 p. 28: C4 regenerates **every** block in an edit's bounding box; p. 99: shipped in the C4 Engine |
| Technical | `10.1145_3730845` (Schott et al., *Sphere Carving: Bounding Volumes for Signed Distance Fields*, ACM TOG 44(4), 2025) | §1: a region can be proven empty only from a field that is a conservative distance bound; a black-box field gives no such proof |
| Technical, repository | `crates/isomesh/src/subgrid/roots.rs`, `extract.rs`, `chunk/dirty.rs`; FINDINGS M-32 | which points the extractor reads, and the rounding caveat for sample positions |

## Craft principle

Baktash, Gillespie and Crane, §1, on what the extractor consumes:

> "We replace 0-dimensional sampling (evaluate f at each grid node), with 1-dimensional root finding (find all zeros of f along each grid edge)."

and on how, for a field that is a black box (this crate's `Sdf`), the zeros are found, §1.3:

> "For instance, we can uniformly sample f along an edge, and use sign changes to identify intersection points—essentially a 1D 'marching segments' approach. While this approach can of course miss intersections, we are no worse off than classic marching (which can also miss intersections)."

So the mesh a Subgrid chunk gets is a function of the field at **points along every edge of the cell's tetrahedra**, not at the eight corners. The marking rule in `mark_edit` is a statement about which points an extractor reads: it compares `before` and `after` at the eight corners, and decides "surface cell" from their signs. For Marching Cubes that is exactly the extractor's input, so the rule is sound. For Subgrid it is a strict subset of the input, so an edit can change the mesh without changing any corner. A plate 0.6 cells thick inside a cell whose eight corners are all outside is the issue's witness: the corners stay outside, and the sample points on the tetrahedra's edges find two roots.

The cell's probe points are the **19 tetrahedron edges** of the Kuhn split (`marching_tetrahedra/table.rs`: the six monotone paths from corner 0 to corner 7): the 12 cube edges, the 6 face diagonals (one per face) and the main diagonal. Each edge is cut into `samples` intervals, so `samples + 1` points per edge, endpoints shared. Sharing across the lattice leaves `1 + 7 (samples - 1)` distinct points per cell: the corner, and `samples - 1` interior points on each of seven edge directions (x, y, z, xy, xz, yz, xyz), counting each lattice edge once however many cells share it: per cell, one x edge, one y edge, one z edge, one xy diagonal, one xz, one yz and one body diagonal.

## How shipping games do it

Lengyel, p. 28, whose method C4 deployed:

> "Whenever a region of the voxel map is edited, all blocks intersecting the region's bounding box have their triangle meshes regenerated."

C4 applies **no** filter at all: every block in the box is regenerated. That is always sound and never cheap, and the filtering this crate does (skip cells that emit nothing before and after) is an optimisation beyond shipped practice, so it has to be sound for the extractor it is paired with or it is a regression against the shipped baseline, which is the case here. The two ways to be sound are to be coarser (C4's) or to compare what the extractor reads.

The tempting middle, "skip a cell if it is provably empty", needs a bound on the field. Schott et al., §1:

> "Our method is agnostic of the function representation and only requires a conservative signed distance field as input."

`Sdf` in this crate is a black box and makes no such promise, and the extractor's own cull (`cell_is_provably_empty`) is off unless the caller supplies a Lipschitz constant, which `VoxelVolume` does not carry. Without a bound, an all-outside cell cannot be called empty from its samples, which is the whole content of #11.

## Model adopted

`mark_edit_subgrid(layout, before, after, min_cell, max_cell, samples, dirty)` in `isomesh::chunk::dirty`, beside `mark_edit`, with one rule per cell over the probes above:

- **value changed**: any probe differs between `before` and `after` (bitwise `total_cmp`, as `mark_edit`);
- **surface before / after**: some edge has a sign change between two consecutive probes (`f < 0` is inside, zero is outside, as `all_roots`), i.e. the extractor finds at least one root there;
- **output changed**: value changed and (surface before or surface after). Only these dirty the owner chunk;
- sign changed and swept: the same definitions over probes.

`samples == 1` is the corner rule exactly (the cube edges connect the eight corners, so some edge changes sign if and only if the corners disagree). That is an invariant the tests hold: the new function with `samples = 1` returns `mark_edit`'s report, field by field.

`mark_changed_chunks` picks the function from `volume.extractor`: `Subgrid { samples }` takes the new one, everything else keeps `mark_edit`. The extractor decides the rule, so a caller cannot pair the wrong one.

Interior probe positions use the extractor's own expressions, `from + (to - from) * t` with `t = i * (1 / samples)`, shared with `all_roots` through `edge_point` and `sample_parameter` so the two cannot drift. The two end points of an edge are the lattice points, `ChunkLayout::world_of_sample` of the global index, evaluated once and shared by every edge that meets there, as in `mark_edit`. The extractor evaluates an edge's end as `from + (to - from) * 1`, which can differ from `to` in the last place, so the end probes are an ulp-level approximation of what it reads (listed under Deviations).

Memory is one layer at a time (two lattice planes and the slab between them, packed one byte per probe), so a whole-volume region does not allocate a volume-sized probe grid.

## Parameters

| Parameter | Value | Why |
|---|---|---|
| `samples` | the extractor's own | the point set is the extractor's; a mismatch is the bug |
| probe evaluations | `1 + 7 (samples - 1)` per cell per field, amortised | counted by a test with a counting field, not asserted |
| size gate | total probe evaluations at most `u32::MAX`, else `IndexSpaceExhausted` | same bound and error as `mark_edit`, applied to the larger count |
| `samples == 0` | `InvalidCellSize`, as the extractor's constructor | the extractor cannot exist, so no rule can match it |

## Validation targets

1. The issue's own test: the plate marks chunk 2, and so does its mirror (plate removed).
2. **A soundness oracle, which is what makes this a measurement and not an argument**: over a sweep of edits (sphere moves, thin plates added and removed, plates placed against an existing surface where no corner changes, brushes straddling chunk faces), extract every chunk with Subgrid before and after; every chunk whose vertex or index buffer differs must be marked. Report how many marked chunks did not change (the over-mark).
3. `samples == 1` equals `mark_edit`'s report and marked set on the existing test fields.
4. The probe count per cell is `1 + 7 (samples - 1)` on a counting field.
5. An edit that changes nothing marks nothing; `samples == 0` is refused; the region gate refuses an enormous box rather than allocating.
6. A Marching Cubes volume marks exactly what it marked before: the existing plugin tests pass unchanged.

## What great looks like

A game that carves a Subgrid volume, or snows on it, sees the mesh change on the frame the edit lands: a thin wall added next to an existing surface appears, one removed disappears, and the collider built from the mesh follows. The cost is visible and documented: Subgrid pays for the sub-cell probes at marking time in proportion to the region (`1 + 7 (samples - 1)` evaluations per cell per field), while a Marching Cubes volume pays what it paid before. No mesh is stale; some chunks are remeshed to the same bytes, and the number is measured.

## Deviations

- The rule is exact with respect to the probe values, **not** to everything the extractor reads. Roots are refined by bisection, and vertex normals come from the field's gradient at each root, and neither is a probe. An edit that lies strictly between two probes, moves no probe value, and moves a root only by sub-sample-spacing, is invisible to marking. The extractor has the same aliasing (`roots.rs`: "a pair of roots closer together than the sample spacing is invisible"), and the paper prices it the same way (§1.3). Marked as a design default and documented on the function.
- Sample positions come from `world_of_sample` (one expression from the global index), where the extractor computes `origin + cell_size * local`, and an edge's two ends are the shared lattice points where the extractor computes `from + (to - from) * t` (`t = 0` and `t = 1`). The expressions are bit-equal for a power-of-two cell size and can differ in the last ulp otherwise (M-32). Only a probe within an ulp of a zero could land on the other side; the same caveat `mark_edit` already carries, and a design default here.
- Not examined here: Dual Contouring reads the field's gradient at each crossing (`hermite.rs`: one `Sdf::gradient` per crossing), which is not a corner, so the corner rule may not be exact for it either. #11 is about Subgrid and this pass leaves Dual Contouring and Surface Nets on `mark_edit`. Logged as a candidate.

## Pending index

- The corpus search (`distill_search`) answered `500 gpu busy` throughout this pass (`bge-m3` needs 5000 MB free, 3.3 GB were), so the passages above are quoted from `markdown_read` of the same stems, not found by a query; no new document was searched for or converted.
- No source on **incremental dirty-region tracking for sparse voxel edits** was found; C4's whole-box rule is the only one this pass found, from the passages it read. A search for it belongs on the list when distill is up.
