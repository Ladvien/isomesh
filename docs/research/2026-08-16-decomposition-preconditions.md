# What a convex decomposer actually requires of its input

**R-011 / P-18.** Numbers in `docs/experiments/p-18.csv`. Findings: M-300, ✗20.

The question this answers: **if isomesh's meshes are going to be decomposed into colliders downstream,
which of this crate's validity checks are entry conditions for that, and does `ColliderReadiness`
already report them?**

The premise it was written against was that all of them are:

> Every ACD method assumes closed, watertight, 2-manifold, self-intersection-free, consistently
> oriented input.

**That is false, and two of the four methods audited say so in their own words.** What is true is
narrower, more useful, and points at a different gap than the one predicted.

---

## The four methods, from their own papers

| method | what it requires of the input | source |
|---|---|---|
| **V-HACD** (Mamou 2016) | **Nothing.** It voxelizes first — *"V-HACD addresses this by voxelizing the input, guaranteeing a solid input for decomposition **at the cost of some accuracy**"* | Andrews 2024, `10.1145/3641519.3657479`, §2.1.3 |
| **CoACD** (Wei et al. 2022) | A **manifold mesh**: *"we follow [Thul et al. 2018] to utilize triangle meshes … and directly cut the **manifold meshes** with 3D planes… Each resulting part is still a manifold mesh."* Its published implementation repairs bad input rather than refusing it — *"runs an implicit surface reconstruction on problematic inputs – creating a slightly offset version of bad inputs, while still being more precisely accurate on clean inputs"* | `10.48550/arXiv.2205.02961` §5; Andrews 2024 §2.1.3 |
| **VisACD** (Fokin & Savva 2026) | **Watertight**, obtained by preprocessing: *"we preprocess the meshes using SDF remeshing to make the meshes watertight, limit the number of vertices, and make vertex density more uniform"* | `10.48550/arXiv.2604.04244` §4 |
| **CPD** (Knodt & Gao 2026) | **Nothing.** *"Unlike prior work, our approach handles **non-manifold, non-watertight meshes directly without preprocessing**."* Fig. 1 is *"a complex non-manifold mesh with boundaries"* | `10.48550/arXiv.2602.07369` §4 |

Andrews 2024 gives the taxonomy that explains the spread: methods differ along three axes, the third
being *"whether to perform **preprocessing or discretization** of the input."* Input cleanliness is a
**design axis**, not a shared requirement.

### So the requirement is not a gate, it is a quality axis

Every method either needs nothing or repairs the input itself. What clean input buys is **skipping the
repair pass and the accuracy it costs** — V-HACD's voxelization error, CoACD's *"slightly offset
version"*, VisACD's SDF remesh. That is a real and quantified benefit, and it is a different claim
from "the decomposer will refuse a mesh that fails T-001."

### The 35% figure is about output, not input

VisACD's *"merging produces intersecting convex hulls in 35% of cases"* describes the **hulls CoACD
emits**, not the mesh it is given — *"We do not evaluate methods that produce decompositions with
intersecting convex hulls."* Reading it as an input precondition is what put self-intersection-freedom
on the list. **No method in this audit requires a self-intersection-free input.** CoACD instead
*guarantees* intersection-free output by construction: *"In this way, we ensure convex hulls of the
decomposed components are intersection-free."*

---

## Coverage against `ColliderReadiness`

**7 required preconditions across the audit. 5 covered. 2 not — and they are the same field twice.**

| precondition | `ColliderReadiness` field | covered |
|---|---|---|
| closed / watertight | `boundary_edges` | ✅ |
| 2-manifold **edges** | `non_manifold_edges` | ✅ |
| consistent orientation | `inconsistently_oriented_edges` | ✅ |
| 2-manifold **vertices** | — | ❌ |

`MeshReport` computes `non_manifold_vertices` with the link walk. `collider::from_report` does not
copy it across, so a **bowtie** — two cones sharing an apex — passes `supports_inside_outside()` with
every edge count at zero. `CLAUDE.md` already states why that configuration is invisible to edge
counts: *"two cones sharing an apex have 2k of each, every edge has exactly 2 faces, and χ can come
out right."* A bowtie is not a 2-manifold, and CoACD's plane cutting is stated over manifold meshes.

**That is P-18's falsifier, and it is not the one that was registered.** The prediction named
self-intersection-freedom as the standing candidate; that turned out not to be an input precondition
of anything. The gap is one line in `from_report` — tracked as **T-021**.

---

## What this means for isomesh

1. **`collider::readiness()` per shard remains the right shape** (M-297), and this sharpens why: it
   reports what a decomposer would otherwise have to repair.
2. **Self-intersections per 1,000 triangles stays a recorded metric, not a gate.** The conditional
   that would have promoted it — *"if decomposition ever enters the pipeline"* — does not fire,
   because no audited method requires the input property.
3. **CPD remains the best fit for plane-cut shards**, and now for a second reason: it is the one
   method that needs nothing at all from the input, which suits geometry produced by cutting where
   caps may be imperfect. Its two costs are unchanged — it does not guarantee non-overlapping
   primitives, and it needs a caller-supplied target primitive count.
