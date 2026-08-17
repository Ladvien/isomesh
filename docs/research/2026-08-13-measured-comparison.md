# Seven isosurface extractors, measured

**2026-08-13, re-derived 2026-08-16.** Marching Cubes, Marching Cubes + asymptotic decider, Marching
Tetrahedra, Surface Nets, Dual Contouring, Manifold Dual Contouring and subgrid Marching Tetrahedra, on
**eight** fields plus **two real scanned volumes**, from one codebase. Speed, accuracy, topology,
triangle quality, whether the mesh separates what the field separates, what it costs to keep the air
region's connectivity up to date as you dig, and — the part nobody publishes — what fraction of a
*usable* mesh the extraction actually is.

This exists because the comparison does not. **No paper since 2020 benchmarks Marching Cubes against
Surface Nets against Dual Contouring, and Surface Nets has no credible published timings at all**
(V-17, literature review). What circulates instead is folklore, and this project has now falsified
enough of it to be worth writing down — **including two of its own earlier conclusions**, which is the
most interesting thing in here and has a section of its own.

Every number below has a committed benchmark that produced it and a `FINDINGS.md` entry that owns it.

> **What the 2026-08-16 re-derivation changed, in one place.**
>
> - **§1's headline is retired (✗25).** *"Marching Cubes is alone in the good corner"* was true of a
>   self-intersection detector with a defect. Marching Tetrahedra's `3.405 per 1k` was a **false
>   positive**; the good corner has **three** occupants.
> - **§4's conclusion is retired (M-308), and not by a better measurement.** *"On Zen 3 there is no
>   crossover at any resolution"* is false on that same machine today, because **we optimised Surface
>   Nets** (A-023/A-024, 4.26×, byte-identical output). Nothing about the machine changed.
> - **§2 loses its `csg_difference` row.** Phase 11 established that field is a Pseudo-SDF, so the
>   accuracy harness now correctly declines to publish a distance for it.
> - **Eight fields, not seven.** `noise_cavity` landed at A-002e and is the worst field in every
>   topology column for every method. Every figure here is re-derived over eight.
> - **One machine, not two.** The Apple M5 half is stale and is not quoted (M-005).
> - **§6 through §10 are new** — R-024's sealing audit, T-026's triangle quality against a published
>   baseline, R-022a's incremental-connectivity cost, R-020's edit-proportionality split and R-025's
>   accuracy floor. No published comparison contains any of them.
> - **Real volumes, for the first time.** `fuel` and `bonsai` from Open SciVis, fetched rather than
>   vendored (M-006). §7 is where they change a conclusion.

---

## Method

| | |
|---|---|
| **Fields** | the **eight** reference fields — `sphere`, `torus`, `box_exact`, `csg_difference`, `thin_plate`, `gyroid`, `fbm_terrain`, `noise_cavity` |
| **Grids** | 33³ and 65³ for topology and accuracy; 16³–256³ for the speed sweep; 17³/25³/33³ for sealing |
| **Timing** | median of repeated runs in one process, single-threaded, after warmup |
| **Machine** | AMD Ryzen 9 5900X (Zen 3, Linux, x86-64). **The Apple M5 half of every timing here is stale** and is not quoted — M-005 |
| **Benches** | `benches/shootout.rs`, `benches/family.rs`, `benches/stage_breakdown.rs`, `benches/experiment_p21.rs` |
| **Data** | `docs/measurements/*.csv` and `docs/experiments/p-21.csv`, committed |

Three deliberate restrictions, each of which removes a way to publish a number that is not there.

**Hausdorff is reported only where the field publishes an exact distance.** That excludes `gyroid`,
`fbm_terrain`, `noise_cavity` — and, since Phase 11, **`csg_difference`**. It used to answer
`is_exact_distance() -> true // away from the seam`, a comment admitting the invariant is false on a
function returning true. `FieldBound` replaced that predicate, `csg_difference` is a **Pseudo-SDF**
(Marschner et al., `10.1145/3610548.3618170`), and its accuracy column is now blank rather than wrong.
**Every `csg_difference` accuracy figure in the pre-2026-08-16 version of this document was measured
against a field that is not a distance function.**

**The shootout's timing column is not used for any speed claim here.** It runs three timed passes per
point and is a few percent noisy. The shootout is the *topology and accuracy* authority; `family` is
the timing authority.

**Timings are one machine and one field.** `family` sweeps `sphere` in `f32`. §4 says so where it
matters, because the claim it retires was a two-machine, eight-field claim and the replacement is
narrower.

---

## 1. The good corner has three occupants, and one of them was wrongly excluded

Two properties matter to anything downstream — is the mesh manifold, and does it intersect itself.
They are independent, so there are four corners.

Eight fields, two grids, one process (`docs/measurements/shootout.csv`):

| | non-manifold edges | worst self-intersections / 1k |
|---|---:|---:|
| **marching cubes** | **0** | **0** |
| **marching cubes + decider** | **0** | **0** |
| **marching tetrahedra** | **0** | **0** |
| surface nets | 747 | 0 |
| dual contouring | 747 | 29.745 (`noise_cavity`) |
| manifold dual contouring | 90 | 45.232 (`noise_cavity`) |
| subgrid marching tetrahedra | 446 | 0 |

**The earlier version of this table said Marching Tetrahedra scored `3.405 per 1k` on `csg_difference`
and that Marching Cubes was alone in the good corner. Both are wrong (✗25).** `99415af`, landed the day
after this document was first written, rewrote the straddle test in `validate::self_intersection`: a
triangle merely **touching** another's plane had been counted as a transverse crossing. Marching
Tetrahedra splits each cell into six tetrahedra that **share faces**, so tangential contact is not an
edge case for it — it is the common case, and it was the whole of the case against it.

Read the rest carefully, because it is still not "Marching Cubes is best". The non-manifold edges under
the dual methods are **not a bug** — they are a property of one-vertex-per-cell, which cannot represent
a cell carrying two sheets of surface (M-4, M-15). What the table says is narrower and more useful: *if
you need both properties and will not pay for a manifold variant, the primal methods have them.*

**Two things the eighth field changed.** `noise_cavity` is the worst field in **every** column for
**every** method — it exists because none of the other seven produces a cell with an interior ambiguity,
0 of 68,385 surface cells (M-208). And **Manifold Dual Contouring no longer takes the zero**: 90
non-manifold edges, all of them `noise_cavity`. That is ✗19 — its own paper says the uniform-grid dual
*"is always a manifold because … the dual preserves the topology of the surface"*, and the premise holds
where the conclusion does not.

**Self-intersections appear only at 33³ and never at 65³.** All six non-zero rows are the coarse grid,
for both duals. Refining removes them, which is M-15's mechanism — *any feature thinner than one cell
forces two sheets through it* — seen in a third place, after M-60's falling manifold-split rate.

**The manifold guarantee is cheaper than its reputation.** Extra vertices over plain Dual Contouring:
**exactly 0 on five of seven fields at every resolution**, and on the other two `gyroid` 3.13% → 2.05%
→ 0.53% and `fbm_terrain` 1.70% → 0.84% → 0.77% across 17³/25³/33³ — *falling* with resolution. Run-time
cost is a median `1.046×` (M-60). Nielson's *"typically comprise about 1.3% of all configurations"* is a
statement about the 256-case table, not about a scene.

---

## 2. Dual Contouring is 101× more accurate where features are sharp, and nowhere else

Symmetric Hausdorff at 65³ (M-54, re-derived):

| field | marching cubes | marching tetrahedra | surface nets | dual contouring | subgrid MT |
|---|---:|---:|---:|---:|---:|
| `box_exact` | 7.217e-2 | 5.103e-2 | 7.217e-2 | **7.145e-4** | 5.103e-2 |
| `thin_plate` | 4.593e-2 | 4.593e-2 | 5.951e-2 | **5.892e-4** | 4.443e-2 |
| `sphere` | 1.329e-3 | 1.413e-3 | 2.251e-3 | 1.093e-3 | **1.083e-3** |
| `torus` | 3.976e-3 | 4.692e-3 | 5.420e-3 | **2.478e-3** | 3.412e-3 |

`box_exact` is **101×** and `thin_plate` **77.9×**; `sphere` is 1.2× and `torus` 1.6×. Two orders of
magnitude exactly where the features are sharp, and nothing at all where they are not. That is the whole
case for the QEF, stated as a number: **it buys a corner and it buys nothing else.** On a smooth field
you are paying a 3×3 solve per surface cell for a 20% improvement.

Manifold Dual Contouring's column is omitted because it is **identical to Dual Contouring's on every
field** — splitting a cell by cycle moves vertices without moving the ones that were already right.

**`csg_difference` is not in this table and used to be.** It scored `3.7×` in the earlier version. That
number was computed against a Pseudo-SDF whose error is concentrated at exactly the seam the figure was
about, and the harness now declines to produce it. See *Method*.

**The corresponding cost used to be surprisingly small and is no longer small.** M-25 measured Dual
Contouring at **3% over Surface Nets** at 256³, and concluded that both are dominated by the *shared*
dual topology rather than by vertex placement. That was true and it is now **33%** — 13.06 against 9.80
ns per sample at 256³.

**Nothing about the QEF changed; the thing it was a small fraction of got 4.26× smaller.** A-023 and
A-024 optimised exactly the shared topology M-25 named as dominant, so the solve is now a much larger
share of a much smaller total. This is §4's lesson arriving in a second section: *"the feature-resolving
solve barely registers"* was a statement about one implementation's balance, not about the algorithm's,
and it inverted without anyone touching the solve.

**And the clamp that makes Dual Contouring safe is free.** Constraining the vertex to its own cell takes
self-intersections to **exactly 0 on five of seven fields** and cuts `gyroid` 23× and `fbm_terrain`
13.7× — while the corner gap on `box_exact` is **0.0057 cells either way, identical**, because a convex
corner's solution is interior to its cell and the constraint never binds where the feature is (M-28).

---

## 3. A pre-registered prediction, half confirmed and half falsified

Two predictions about Marching Tetrahedra were written down *before* the measurement, which is the only
arrangement under which a prediction is worth anything.

**O-13 — "3.0× the vertices of Marching Cubes, converging from above" — confirmed exactly.** Measured on
`sphere`: **3.036 / 3.026 / 3.003** at 33³/49³/65³.

**O-14 — "≈2.6e-3 Hausdorff, 1.86× worse than Marching Cubes" — falsified.** Measured **1.413e-3** and
**1.06×** — M-55 recorded 1.4386e-3 and 1.043× over the seven-field grid, and this is the same result
re-derived over eight. And the clause flagged in advance as counterintuitive, *"more vertices **and** worse
accuracy"*, is the half that does not hold: Marching Tetrahedra beats Surface Nets (2.251e-3) rather
than losing to it, and on the sharp fields it beats **Marching Cubes** — `box_exact` 5.103e-2 against
7.217e-2 — because its extra edge families sample a corner from more directions.

The vertex ratio has a mechanism the prediction did not need and turns out to have: it is **exactly 4.0**
for any surface orientation inside one octant and **2.0** across a sign change, so `2.992` is an average
hiding a factor-of-two spread (M-52).

**✗25 adds a third result to this section**: Marching Tetrahedra is also intersection-free on every
field, which the earlier measurement denied it. It costs ~3× the triangles for accuracy within 6% of
Marching Cubes and better on sharp fields — a worse trade than it sounds, and a better method than this
document used to say.

---

## 4. Speed: the section our own optimisation falsified

This is the most important section here and it is not about speed.

The earlier version of this document said:

> *"Surface Nets is slower than Marching Cubes on every field measured, on both machines. … On Zen 3
> there is no crossover at any resolution — Surface Nets is 2.46× behind even at 16³. **So the crossover
> is a property of one machine's cache behaviour**, and the degradation is a property of the algorithm."*

That was a correct measurement and a falsification of real folklore (✗14, *"Surface Nets is the cheapest
thing in the family and the natural default"*). **On the same machine today it is false** (M-308).
`docs/measurements/family.csv`, `f32`, `sphere`, one binary, one run:

| samples per axis | marching cubes | surface nets | SN / MC |
|---:|---:|---:|---:|
| 16³ | 9.24 ns | **8.63** | **0.93×** |
| 24³ | 8.68 | **8.32** | **0.96×** |
| 32³ | 8.11 | **7.96** | **0.98×** |
| 48³ | 7.86 | 8.10 | 1.03× |
| 64³ | 7.64 | 8.41 | 1.10× |
| 128³ | 7.56 | 8.78 | 1.16× |
| 256³ | 7.60 | 9.80 | **1.29×** |

**The crossover exists on Zen 3, at about 40³.** `SN/MC` at 256³ has gone **5.43× → 1.29×**.

**Nothing about the machine changed.** A-023 made the dual mesher's loop axis a const-generic parameter;
A-024 forced the row length odd to break a 64 KiB cache-set aliasing period. Together, **4.26× on the
dual path with byte-identical output** (M-285, M-287) — this document's own §1 and §2 numbers are
unchanged to the last digit, which is how we know.

So the published claim was never about the algorithms. ***"Surface Nets is slower than Marching Cubes"*
was a measurement of one implementation of Surface Nets at one moment, and it read as a statement about
the method — including to us, who wrote it.** ✗14 was folklore that measurement falsified; this section
is measurement that **our own optimisation** falsified. Same failure, one level up.

**What survives is the shape, and the shape was always the finding.** Surface Nets' per-sample cost still
*rises* with resolution — 7.96 ns at 32³ to 9.80 at 256³, **+23%** — while Marching Cubes' falls and then
flattens, 9.24 → 7.60. That is the cache term M-21 named and M-286 can now see, and the 4.26× did not
touch it. **The ratio was never the result; the derivative was.**

**A number nobody has quoted.** Subgrid Marching Tetrahedra runs at **1,553–1,704 ns/sample** across the
same sweep — **200× Marching Cubes**, and about 100× classic Marching Tetrahedra. It samples 16 points
per tet edge, 576 field evaluations per cell against Marching Cubes' 8 shared corner samples (M-98).
That constant is the method, not the implementation.

**`f64` costs 8–10%, not 2×**, on extraction paths with no matrix solve in them (M-23). **The asymptotic
decider is free**: −0.3% to +3.3% across four (field, grid) points (M-42).

### The fixed cost that is not there

`t = a + b·n³` was fitted expecting a large intercept at small grids. There isn't one. Marching Cubes
fits **a = 0.5118 ms — 0.64% of the largest run — r² = 0.99976** (M-62), and the dual methods fit `a < 0`,
which is physically impossible and is the model reporting that it does not describe them.

That diagnostic then failed to reproduce, and **the failure is the useful part**. Fitting live over
17³–89³ across eight independent runs: Surface Nets' `a` is negative in **7 of 8**, Marching Cubes
positive in 7 of 8, Dual Contouring positive in **8 of 8** — every value under **0.18 ms** in absolute
value, against the originally committed **−3.13 ms** (M-134). Forty times smaller.

So the amended rule is: **report the per-sample curve, not the intercept.** The sign carried the meaning;
the number never did. It is the same lesson as the table above, arrived at independently.

---

## 5. The contour is 36% of the job, and the biggest stage is the collider check

This is the part with no published counterpart, and it is the one that should change what anyone
optimises.

A paper reports the *contour*. A consumer pays for a **usable mesh**: the contour plus everything that
has to happen before a renderer or a physics engine will accept it. The one published comparison this
project has of the two puts contouring at 68 ms against halfedge construction at 58 ms — **54%**.

Measured here, **eight** fields at 33³ and 65³, `f64` (M-135, re-derived):

| stage | mean share |
|---|---:|
| contour | **35.9%** |
| normals | 0.3% |
| weld | 22.8% |
| collider check | **41.0%** |

**Optimising the extractor alone therefore buys at most 1.56× on the whole job.** That is the ceiling,
and it is worth knowing before starting rather than after.

Three qualifications, all of which survive.

**Normals are effectively free (0.3%)** because every reference field here overrides `Sdf::gradient`. The
measured pass is the cost a consumer pays whose field has *no* analytic gradient — a sampled volume, an
imported scan — and even then it is noise.

**The collider figure times `collider::readiness`, which validates.** A physics engine building its own
structure is a different cost, and a shipping game may validate once rather than per chunk. But
welded-manifold-correctly-wound is the contract, and checking that contract is the single largest line
in the budget.

**There is no bar for upload, deliberately.** Its CPU half is a *move* — the Bevy sink writes straight
into the arrays a `Mesh` owns — so there is no pass to time, and its GPU half needs a device the core
crate must not have.

### And the single number is a fiction

The contour's share runs **14.2% to 77.5%** across the eight fields — a **5.5× spread on identical code**
(M-136):

| field | 33³ | 65³ |
|---|---:|---:|
| `gyroid` | 14.2% | 24.1% |
| `sphere` | 15.0% | 21.4% |
| `box_exact` | 15.0% | 21.5% |
| `noise_cavity` | 16.9% | 19.2% |
| `torus` | 17.6% | 26.2% |
| `csg_difference` | 24.5% | 23.6% |
| `thin_plate` | 33.4% | 38.8% |
| `fbm_terrain` | **65.6%** | **77.5%** |

The variable is **how expensive the field is to sample**, not anything about the mesher. fBm noise swamps
the fixed post-processing; a cheap analytic field leaves the post-passes dominant.

So *"is it worth optimising the contour"* cannot be answered in general. **A game on procedural terrain is
in `fbm_terrain`'s regime and should optimise the extractor. One on authored CSG is in `sphere`'s and
should not.**

---

## 6. Does the mesh seal what the field seals?

**New in the 2026-08-16 re-derivation, and no published comparison contains it** (R-024, M-307).

Every metric above judges a mesh against **itself** — manifoldness, orientation, Euler characteristic —
or against the field's **geometry**. None asks whether the mesh partitions *space* the way the field's
sign does, and neither claim implies the other: a mesh can be closed, manifold, correctly wound and
Hausdorff-close while sealing a passage the field leaves open, or opening one it seals. For a game that
is the whole question, because *is this cave sealed* is asked of the collider and answered by the field.

The probe is a grid edge. The field says two adjacent samples straddle the surface iff their signs
differ; the mesh says they are separated iff the segment between them is crossed an **odd** number of
times. **The test is Wojtan, Thürey, Gross & Turk's complex-edge test** (`10.1145/1778765.1778787`) and
is theirs, not ours (V-37); running it as a correctness audit of *extraction* is what is new.

Eight fields × three resolutions, 9,151,296 probes (`docs/experiments/p-21.csv`):

| | sealed | worst hole count | on a domain face |
|---|---:|---:|---:|
| **marching cubes** | **24 / 24** | — | — |
| **marching cubes + decider** | **24 / 24** | — | — |
| **marching tetrahedra** | **24 / 24** | — | — |
| surface nets | 21 / 24 | 190 (`fbm_terrain` 33³) | **190 / 190** |
| dual contouring | 21 / 24 | 190 (`fbm_terrain` 33³) | **190 / 190** |
| manifold dual contouring | 21 / 24 | 190 (`fbm_terrain` 33³) | **190 / 190** |
| subgrid marching tetrahedra | 18 / 24 | 4 (`fbm_terrain` 17³) | 1 / 4 |

**All three duals report the identical count** — 92 / 138 / 190 at 17³/25³/33³ — and every hole is on a
face of the sampled domain. Three independent implementations agreeing to the unit is what says the
mechanism is *one vertex per cell* rather than any particular solve: a dual emits one quad per
sign-changing grid edge, and that quad needs all **four** cells around the edge. At the domain face only
one or two exist, so no quad is emitted. A primal method emits per *cell*, meshes every cell it has, and
seals the same edge.

**For a chunked world that face is the chunk seam**, so this is the measured statement of why a dual
chunk's collider is not watertight on its own.

**The property that splits the family is not primal versus dual.** It is whether the method puts its
crossing **on the probed grid edge**. Subgrid Marching Tetrahedra is primal by family and scores worst,
because it samples along tet edges *inside* the cell. Its residue is recorded and deliberately **not**
called a defect: this audit compares against the sign pattern at the grid corners, which is all the other
six methods ever see, and subgrid MT looks between them — so where it disagrees it may be the more
faithful of the two.

`fbm_terrain` is the only reference field whose surface leaves through the sides, which is why it is the
only field where any of this is reachable.

---

## 7. Triangle quality, and a band that only transfers on the right kind of field

**New in the 2026-08-16 re-derivation** (T-026, M-309, M-310). The metric is Grosso & Zint's **mean
ratio** (`10.1007/s00371-021-02139-w` §5):

```text
q = 4√3 · A / Σᵢ lᵢ²        1 for equilateral, 0 for degenerate
```

Chosen over the `AR > 4` figures the differentiable-isosurfacing line reports, because those are
measured on meshes extracted from **learned** fields inside an optimisation loop and are not comparable
with meshing on a uniform grid (V-38). Grosso & Zint mesh uniform grids and report Marching Cubes,
topologically correct Marching Cubes and Dual Contouring by name.

Averaged over 17³/25³/33³:

| extractor | `sphere` | `torus` | `box_exact` | `thin_plate` | `gyroid` | `fbm_terrain` | `noise_cavity` |
|---|---:|---:|---:|---:|---:|---:|---:|
| marching cubes | 0.728 | 0.714 | 0.858 | 0.866 | 0.687 | 0.691 | 0.688 |
| surface nets | **0.814** | **0.797** | **0.863** | 0.743 | **0.795** | **0.770** | **0.758** |
| dual contouring | 0.815 | 0.736 | 0.845 | 0.743 | 0.727 | 0.600 | 0.621 |
| marching tetrahedra | 0.674 | 0.675 | **0.250** | 0.499 | 0.652 | 0.652 | 0.637 |

**Surface Nets beats Marching Cubes on all eight fields. Dual Contouring beats it on the smooth ones and
*loses* on the rough ones** — 0.600 against 0.691 on `fbm_terrain`. The QEF places a vertex to fit
planes, and where there is no feature to fit it places it badly; a centroid has no such failure mode.
**The published claim that the dual is better is a claim about smooth data.**

**Marching Tetrahedra scores 0.250 on `box_exact`** against Marching Cubes' 0.858 — while being the
*most accurate* primal method on that same field, §2's 5.103e-2 against 7.217e-2. A 6-tetrahedron
decomposition cuts an axis-aligned face along tet diagonals and the triangles that result are needles.
**Accuracy and triangle quality point in opposite directions there by a factor of three**, which is not
a trade-off anyone lists.

### The band transfers on real volumes and not on analytic ones

Their Marching Cubes sits at **0.65–0.71**. Ours measures **outside** that on analytic fields — 0.7785,
0.7131, 0.7510 — which read as a falsification until the same metric ran on real scanned data (M-006,
M-310):

| volume | marching cubes | surface nets | dual contouring | manifold DC | subgrid MT |
|---|---:|---:|---:|---:|---:|
| `fuel` 64³ | **0.7006** | 0.7976 | 0.6876 | 0.6891 | 0.6622 |
| `bonsai` 256³ | **0.6888** | 0.7957 | 0.6443 | 0.6444 | 0.6770 |

**Both inside 0.65–0.71.** Their band was measured on CT and simulation data and it reproduces on CT and
simulation data — **the band is a property of the input class, not of the implementation**, and the
earlier reading was measuring the wrong kind of field. Their *Dual Contouring* figure does not reproduce
at all (0.6649 here against their 0.82–0.86, below our own Marching Cubes), and that asymmetry is
recorded unexplained.

**On real data, Marching Cubes emits 0 non-manifold edges on a million-triangle CT surface**, and
Manifold Dual Contouring takes `bonsai` from **1,776 non-manifold edges to 85** — the manifold
construction earning its keep on the input class it was designed for, where on the analytic fields it
was 0 extra vertices on five of seven (M-60).

### Subgrid Marching Tetrahedra costs 471× on real data, and that is the number to decide on

**It refused `bonsai` outright until A-028**, on a zero-gradient error, and the diagnosis is worth more
than the fix (M-316, M-317). The cause was not a plateau: the failing corners are **local extrema in
every axis**, one with neighbour slopes of **∓19** — steep, and *exactly* symmetric because `u8`
quantisation put both neighbours on the same integer. **A central difference is identically zero at a
local extremum however steep the field is**, and `SampledField` was inheriting one instead of supplying
the trilinear gradient it had in closed form.

It now meshes: **1,572,901 vertices, 3,138,925 triangles, 879 non-manifold edges** — and **232.3 s**
against Marching Cubes' 0.49 s on the same volume, which is **471×**. §4's sphere figure was ~200×
(M-308); **the constant is more than twice as bad on real data**, and that is the number anyone
considering this extractor for a scan should be deciding on rather than the sphere one.

**483 tetrahedra are declined, clustered around 33 points.** Where the field has a critical point *on*
the isosurface there is no normal — not a missing one, an absent one — so those tetrahedra are skipped
and the report says where. 475 are `Degenerate` and **8 are `IllConditioned`**, a precision problem
rather than a topological one, which a single count would have hidden. 483 / 33 = **14.6** tetrahedra
per singular corner, so the holes are 33 small clusters rather than 483 scattered ones.

**If your data is integer, contour at a half-offset isovalue.** `127.5` rather than `127`: integer
samples cannot equal a half-integer, so no sample sits on the isosurface and the degeneracy never
arises. On `bonsai` against an integer isovalue, **16,284 of 529,508 surface-cell corners — 3% — sit
exactly on the surface**, which is where this extractor asks for a normal at a grid point. One line,
and it removes the case rather than reporting it.

---

## 8. Repairing connectivity costs the edit, not the lattice

**New** (R-022a, M-311). *Is this cave sealed? Did I just break through?* are questions about the
connected components of the air region, asked after every edit. Digging removes solid, so air samples
only ever **appear** — and Durfee et al. note that *"an insert can cause at most two trees in `F` to be
joined"*, with no replacement-edge search. A union-find is the entire structure.

One spherical brush of radius 6, **identical at every resolution**, into a solid lattice:

| n | samples scanned | dirty | incremental unions | rebuild unions | incr ms | rebuild ms |
|---:|---:|---:|---:|---:|---:|---:|
| 33 | 35,937 | 925 | **4,872** | 2,436 | 0.028 | 0.2 |
| 65 | 274,625 | 925 | **4,872** | 2,436 | 0.032 | 1.3 |
| 129 | 2,146,689 | 925 | **4,872** | 2,436 | 0.051 | 5.3 |

**The lattice grows 59.7× and the incremental union count does not move by one.**

**The `n³` is in the scan, not in the unions**, and that was worth getting wrong to learn: the rebuild's
union count is *also* flat, because a union-find build unions only air-air edges and the air volume is
the brush. What a rebuild actually pays is **visiting 2,146,689 samples to discover that 925 changed** —
2,321 touched per sample that mattered, and a **104×** wall-clock gap at 129³ that widens with the
lattice.

**Filling is not this problem.** Removing air is a deletion, which needs a replacement-edge search and
which a union-find cannot do at any price. `connectivity::Air` therefore ships `dig` and no `fill`,
and the other half is a different data structure rather than a longer version of this one.

---

## 9. What a local edit costs, and why the answer is two answers

**New** (R-020, M-314). Every published extractor's incrementality axis reads *"full re-mesh"*, and a
voxel game re-meshes after every brush stroke. So: after a local field edit, does the work track the
**edit** or the **grid**?

One spherical brush of radius 5 carved at a sphere's equator, **identical at every resolution**:

| n | cells | dirty cells | vertices | **buffer moved** | **geometric moved** | first moved |
|---:|---:|---:|---:|---:|---:|---:|
| 33 | 32,768 | **792** | 1,758 | 1,348 (77%) | **330** | 414 |
| 65 | 262,144 | **792** | 6,918 | 4,257 (62%) | **322** | 2,609 |
| 129 | 2,097,152 | **792** | 27,822 | 15,706 (56%) | **346** | 12,257 |

**The computation is edit-proportional and exactly so.** The lattice grows **64×** and the dirty set
does not move by one: **792 cells**, from 515 changed samples, comfortably inside the `8k = 4,120` bound
arithmetic gives since each sample is a corner of at most eight cells.

**The output encoding is not.** Vertices that genuinely appeared or vanished — compared as a *set of
positions*, not by index — are **330, 322, 346**: flat, and about the size of the dirty set. Vertices
whose **buffer slot** changed are **1,348, 4,257, 15,706**: growing with the `O(n²)` vertex count, and
**56–77% of the whole buffer**. So more than half the output is rewritten for an edit touching **0.038%**
of the cells.

**The cause is a counter, not the algorithm.** `first moved` is 12,257 of 27,822 at 129³ and essentially
everything after it differs: vertices are appended in scan order and indices name buffer positions, so a
cell emitting a different *number* of triangles shifts every index after it.

**The control is what makes either number mean anything.** An index-wise diff reports the same large
number whether the surface moved or the buffer was merely reshuffled. Separating them turns *"re-meshing
is expensive"* into ***"re-meshing is cheap and the encoding is expensive"*** — and the crate already
names vertices stably *internally*, on `(lower sample, axis)`, then discards that naming when packing.
**R-027** is that fix.

---

## 10. Accuracy has a floor, and it decides which method a better vertex rule could help

**New** (R-025, M-315). The reported Hausdorff is the worst mesh sample, and the harness samples
**vertices and triangle centroids** — so it is `max(vertex, centroid)`. Project every vertex exactly onto
the true surface and the vertex term vanishes, leaving the centroid term on perfect vertices. **That
residue is a floor no vertex placement can go below**, because a flat triangle inscribed in a curved
surface deviates at its interior whatever its corners do.

Over the two smooth fields that publish an exact distance:

| | sphere 17³/33³/65³/129³ | torus 17³/33³/65³/129³ |
|---|---|---|
| dual contouring | 9.2, 17.4, 10.7, 12.3 % | 12.2, 16.8, **1.5**, 8.0 % |
| marching cubes | **21.5**, 14.5, 13.8, 9.0 % | 16.0, 12.1, 12.3, 5.7 % |

**A perfect placement rule is worth 1.5–21.5%, median 12.3%, and only one of sixteen rows reaches 20%.**

**Which term *is* the Hausdorff is the crisper result**, and it is clean:

- **Dual Contouring — the vertex term, 8 of 8 rows. Placement-limited.**
- **Marching Cubes — the centroid term, 8 of 8 rows. Tessellation-limited.**

So a better vertex rule is aimed at Dual Contouring and **cannot help Marching Cubes at all**; refining
Marching Cubes' accuracy means more triangles, not better-placed ones.

**And the QEF is already trading in the direction such a rule would want.** Dual Contouring's *centroid*
error is **better than the floor** at 7 of 8 rows — by **2.9–3.6×** on `sphere` — because the QEF
minimises distance to the tangent planes rather than putting a vertex on the surface. It buys
better-centred facets at the cost of worse-placed vertices, and **pushing its vertices onto the surface
would make the triangles fit worse**. §2's *"the whole case for the QEF is that it buys a corner"* needs
this beside it: on a smooth field it is buying facet fit, quietly.

---

## What surprised us

1. **A conclusion of ours was falsified by an optimisation of ours** (§4). Not by a better instrument,
   not by another machine. The family comparison was measuring an implementation and reading as a
   statement about an algorithm.
2. **A defect in the instrument put a method in the wrong corner for two days** (§1), and the fix's own
   commit message said it had *"inflated a metric this repo quotes"* — while the metric stayed quoted.
3. **The collider check costs more than the contour.** Nothing in the literature suggests looking, and it
   is 41% of the job.
4. **The field decides the ratio, by 5.5×.** "What fraction is the extraction" is a property of the
   workload, not of the code.
5. **Every dual method leaves the domain boundary unsealed, identically** (§6), and no primal method
   does.
6. **A borrowed band failed on our fields and held on real ones** (§7). Grosso & Zint's 0.65–0.71 for
   Marching Cubes reads as falsified on analytic fields and reproduces on CT data, because it was a
   property of the **input class** all along. The eighth field-dependent figure this project has caught,
   and the first that was somebody else's.
7. **Marching Tetrahedra is the most accurate primal method on `box_exact` and by far the worst-shaped**
   — 5.103e-2 against Marching Cubes' 7.217e-2 on Hausdorff, and 0.250 against 0.858 on triangle
   quality. Two metrics, opposite verdicts, factor of three.
8. **The `n³` in a connectivity rebuild is the scan, not the work** (§8). Both union counts are flat;
   what a rebuild pays is visiting 2.1 M samples to find the 925 that changed.
9. **Re-meshing after an edit is cheap and the *encoding* is expensive** (§9). The dirty set is 792 cells
   whether the grid is 32 thousand cells or 2 million; the buffer rewrites 56–77% of itself either way,
   because indices name positions in a sequentially packed array.
10. **A perfect vertex placement is worth about 12%** (§10), and only to Dual Contouring. Marching Cubes'
    error is in its triangles, not its vertices, at every resolution on both smooth fields.
11. **Dual Contouring's facets fit better than perfectly-placed vertices would** (§10), by 2.9–3.6× on a
    sphere. The QEF is not trying to put a vertex on the surface, and the accident is a good one.
12. **Marching Tetrahedra is more accurate than Marching Cubes on sharp fields** — the opposite of the
   registered prediction, and for a reason (more edge families) the prediction did not consider.
13. **The sharp-feature solve is nearly free**, and so is the clamp that makes it safe, and so is the
   asymptotic decider. Three "expensive" features that are not.

## What we would now distrust

Not individual facts — *sources*. Each of these published a single number for a quantity that is
field-dependent, implementation-dependent, or instrument-dependent:

- The `2.76×` greedy-quads merge ratio. Measured range: **1.70× to 256×** (M-56).
- Nielson's `1.3%` manifold-split rate. Measured: **0 on five of seven fields** (M-60).
- The `54%` contour share. Measured: **14.2%–77.5%** (M-136).
- ✗14's "cheapest in the family" — **and its refutation**, which lasted three days (M-308).
- **Our own self-intersection figures before 2026-08-14** (✗25). A metric is only as good as the
  predicate under it, and a tangential-contact bug inflates exactly the methods that generate tangential
  contact.

The pattern is consistent enough to be a rule, and the re-derivation added a clause to it: **a single
figure quoted for a mesher's behaviour is usually one scene's number, on one build, from one
instrument.** Ask which field, which commit, and which predicate — including of the numbers in this
document. Every one of them names its field, and this version names its date.

---

## Reproducing this

```bash
cargo bench --bench shootout          # topology + accuracy -> docs/measurements/shootout.csv
cargo bench --bench family            # per-sample timing   -> docs/measurements/family.csv
cargo bench --bench stage_breakdown   # stage shares        -> docs/measurements/stage_breakdown.csv
cargo bench --bench experiment_p21    # sealing audit       -> docs/experiments/p-21.csv
cargo bench --bench experiment_p22    # triangle quality    -> docs/experiments/p-22.csv
cargo bench --bench experiment_p23    # connectivity repair -> docs/experiments/p-23.csv
cargo bench --bench edit_trace        # edit-proportionality -> docs/measurements/edit_trace.csv
cargo bench --bench placement_ceiling # the accuracy floor   -> docs/measurements/placement_ceiling.csv
cargo bench --bench a028_diagnose     # the zero-gradient diagnosis, on a real volume
cargo bench --bench interior_margin   # the interior decider's margin

./scripts/fetch_volumes.sh            # real volumes; not committed, verified by published SHA-512
cargo bench --bench volumes           # every extractor on them -> docs/measurements/volumes.csv
```

Always `--release`; `cargo bench` handles that.

---

*Ledger entries behind this document: M-4, M-15, M-20, M-21, M-23, M-25, M-28, M-42, M-52, M-54, M-55,
M-56, M-60, M-62, M-98, M-134, M-135, M-136, M-208, M-285, M-286, M-287, M-290, M-307, M-308, M-309,
M-310, M-311, M-314, M-315, M-316, M-317, V-37, V-38, V-39, V-40, V-41, V-43, ✗14, ✗19 and ✗25. Tickets: M-001a,
M-001b, M-002, M-003, M-004, M-006, A-023, A-024, A-028, R-020, R-022a, R-024, R-025, R-026, T-026.*
