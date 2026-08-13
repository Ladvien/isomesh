# Five isosurface extractors, measured

**2026-08-13.** Marching Cubes, Marching Cubes + asymptotic decider, Marching Tetrahedra, Surface
Nets, Dual Contouring and Manifold Dual Contouring, on seven fields, on two machines, from one
codebase. Speed, accuracy, topology, and — the part nobody publishes — what fraction of a *usable*
mesh the extraction actually is.

This exists because the comparison does not. **No paper since 2020 benchmarks Marching Cubes against
Surface Nets against Dual Contouring, and Surface Nets has no credible published timings at all**
(V-17, literature review). What circulates instead is folklore, and this project has now falsified
enough of it to be worth writing down.

Every number below has a committed benchmark that produced it and a `FINDINGS.md` entry that owns it.
Where a figure is machine-specific, it says so, because three of the most interesting results here are
about exactly that.

---

## Method

| | |
|---|---|
| **Fields** | seven reference fields — `sphere`, `torus`, `box_exact`, `csg_difference`, `thin_plate`, `gyroid`, `fbm_terrain` |
| **Grids** | 33³ and 65³ for topology and accuracy; 16³–256³ for the speed sweep |
| **Timing** | median of repeated runs in one process, single-threaded, after warmup |
| **Machines** | Apple M5 (macOS, arm64) and AMD Ryzen 9 5900X (Zen 3, Linux, x86-64) |
| **Benches** | `benches/shootout.rs`, `benches/resolution_sweep.rs`, `benches/stage_breakdown.rs` |
| **Data** | `docs/measurements/*.csv`, committed |

Two deliberate restrictions, both of which remove a way to publish a number that is not there.

**Hausdorff is reported only where the field publishes an exact distance.** What the accuracy harness
computes against `gyroid` is not a distance, and printing it would invent a quantity. So `gyroid` and
`fbm_terrain` have no accuracy column, rather than a wrong one.

**The shootout's timing column is not used for any speed claim here.** It runs three timed passes per
point and is a few percent noisy — noisy enough that the decider reads `0.89×` on `sphere` where it
must be `1.00×`. The shootout is the *topology and accuracy* authority; `resolution_sweep` is the
timing authority. Its topology and accuracy columns are hardware-independent, which is why they carry
no machine attribution.

---

## 1. Marching Cubes wins the property nobody credits it with

Two properties matter to anything downstream — is the mesh manifold, and does it intersect itself.
They are independent, so there are four corners. **Three are occupied, and the algorithm the folklore
treats as the crude baseline is alone in the good one.**

| | non-manifold edges | self-intersections / 1k |
|---|---:|---:|
| **marching cubes** | **0** | **0** |
| **marching cubes + decider** | **0** | **0** |
| marching tetrahedra | 0 | 3.405 (`csg_difference`) |
| surface nets | 128 | 0 |
| dual contouring | 128 | 13.837 (`fbm_terrain`) |
| manifold dual contouring | 0 | 15.434 (`fbm_terrain`) |

Seven fields, two grids, one process (M-53). Each non-Marching-Cubes method fails exactly one property
or both.

Read this carefully, because it is not "Marching Cubes is best". The 128 non-manifold edges under the
dual methods are **not a bug** — they are a property of one-vertex-per-cell, which cannot represent a
cell carrying two sheets of surface, and they appear on exactly the two fields that have such cells
(M-4, M-15). Manifold Dual Contouring is the entry that takes the zero, and it buys that by splitting
cells. What the table says is narrower and more useful: *if you need both properties and will not pay
for the manifold variant, the crude baseline is the one that has them.*

**The manifold guarantee is cheaper than its reputation.** Extra vertices over plain Dual Contouring:
**exactly 0 on five of seven fields at every resolution**, and on the other two `gyroid` 3.13% → 2.05%
→ 0.53% and `fbm_terrain` 1.70% → 0.84% → 0.77% across 17³/25³/33³ — *falling* with resolution, because
refining removes multi-sheet cells rather than creating them. Run-time cost is a median `1.046×`
(M-60). Nielson's *"typically comprise about 1.3% of all configurations"* is a statement about the
256-case table, not about a scene; on a real field the answer is field-dependent and usually zero.

---

## 2. Dual Contouring is 101× more accurate where features are sharp, and nowhere else

Symmetric Hausdorff at 65³, Marching Cubes → Dual Contouring (M-54):

| field | marching cubes | dual contouring | ratio |
|---|---:|---:|---:|
| `box_exact` | 7.217e-2 | **7.145e-4** | **101×** |
| `thin_plate` | 4.593e-2 | **5.892e-4** | **77.9×** |
| `csg_difference` | 7.655e-2 | 2.057e-2 | 3.7× |
| `sphere` | — | — | 1.2× |
| `torus` | — | — | 1.6× |

Two orders of magnitude exactly where the features are sharp, and nothing at all where they are not.
That is the whole case for the QEF, stated as a number: it buys a corner and it buys nothing else. On a
smooth field you are paying a 3×3 solve per surface cell for a 20% improvement.

The corresponding cost is small enough to be surprising. **Dual Contouring costs 3% more than Surface
Nets at 256³ on the M5** (218.9 ms against 212.5) and **6.5% on Zen 3** (877.4 against 823.4) — M-25.
Both methods are dominated by the *shared* dual topology (sampling, the quad walk), not by vertex
placement, so the feature-resolving solve barely registers.

**And the clamp that makes Dual Contouring safe is free.** Constraining the vertex to its own cell
takes self-intersections to **exactly 0 on five of seven fields** and cuts `gyroid` 23× and
`fbm_terrain` 13.7× — while the corner gap on `box_exact` is **0.0057 cells either way, identical**,
because a convex corner's solution is interior to its cell and the constraint never binds where the
feature is (M-28).

---

## 3. A pre-registered prediction, half confirmed and half falsified

Two predictions about Marching Tetrahedra were written down *before* the measurement, which is the only
arrangement under which a prediction is worth anything.

**O-13 — "3.0× the vertices of Marching Cubes, converging from above" — confirmed exactly.** Measured
on `sphere`: **3.036 / 3.026 / 3.003** at 33³/49³/65³.

**O-14 — "≈2.6e-3 Hausdorff, 1.86× worse than Marching Cubes" — falsified.** Measured **1.4386e-3** and
**1.043×** (M-55). And the clause flagged in advance as counterintuitive, *"more vertices **and** worse
accuracy"*, is the half that does not hold: Marching Tetrahedra beats Surface Nets (2.251e-3) rather
than losing to it, and on the sharp fields it beats **Marching Cubes** — `box_exact` 5.103e-2 against
7.217e-2 — because its extra edge families sample a corner from more directions.

The vertex ratio also has a mechanism the prediction did not need and turns out to have: it is
**exactly 4.0** for any surface orientation inside one octant and **2.0** across a sign change, so
`2.992` is an average hiding a factor-of-two spread (M-52). That is why the CSV carries every field
rather than a headline.

---

## 4. Speed: the folklore is wrong in the same direction on both machines

> *"Surface Nets is the cheapest thing in the family and the natural default."* — ✗14

**Surface Nets is slower than Marching Cubes on every field measured, on both machines.** Per-sample
cost from 16³ to 256³ (M-45):

| | 16³ | 256³ |
|---|---:|---:|
| marching cubes, M5 | 24.99 ns | **4.78 ns** |
| marching cubes, Zen 3 | 15.18 ns | 13.19 ns |
| surface nets, M5 | 8.40 ns | 12.66 ns |
| surface nets, Zen 3 | 37.38 ns | 49.08 ns |

**Surface Nets degrades on both machines** — its per-sample cost *rises* while Marching Cubes' falls —
so the effect is the algorithm's memory pattern rather than one cache hierarchy. At 256³ the ratio is
**3.72× on Zen 3 against 2.65× on the M5**.

The one place the two machines disagree is instructive: on the M5, Surface Nets *does* win below a
crossover, because M5 Marching Cubes starts expensive and converges. On Zen 3 there is no crossover at
any resolution — Surface Nets is 2.46× behind even at 16³. **So the crossover is a property of one
machine's cache behaviour, and the degradation is a property of the algorithm.** Quoting the first as
if it were the second is how the folklore got made.

Marching Cubes' marginal cost is **4.75 ns/sample — 211 M samples/s, single-threaded, f32, M5** (M-20).
Against a published `5.42 G voxel/s` on an RTX 2080 Ti that is a **~26× gap**, which is the number a
GPU decision should be argued from rather than from enthusiasm.

**`f64` costs 8–10%, not 2×**, on extraction paths with no matrix solve in them (M-23) — the work is
dominated by field evaluation and branchy table lookup, not memory bandwidth.

**The asymptotic decider is free**: −0.3% to +3.3% across four (field, grid) points (M-42). `sphere`
has no ambiguous face at all, so its difference is the price of *asking*; `gyroid`'s extra ~1.7 points
is the price of *answering*.

### The fixed cost that is not there

`t = a + b·n³` was fitted expecting a large intercept at small grids. There isn't one. Marching Cubes
fits **a = 0.5118 ms — 0.64% of the largest run — b = 4.7389 ns/sample, r² = 0.99976** (M-62), and the
dual methods fit `a < 0`, which is physically impossible and is the model reporting that it does not
describe them.

That last diagnostic then failed to reproduce, and **the failure is the most useful thing in this
section**. On Zen 3, fitting live over 17³–89³ across eight independent runs: Surface Nets' `a` is
negative in **7 of 8**, Marching Cubes positive in 7 of 8, Dual Contouring positive in **8 of 8** — and
every value is under **0.18 ms in absolute value**, against the originally committed **−3.13 ms**
(M-134). Forty times smaller, and on Zen 3 *both* intercepts come back negative and numerically
negligible.

So the amended rule is: **report the per-sample curve, not the intercept.** The sign carried the
meaning; the number never did. `r²` even inverts over the short range — Marching Cubes fits *worst*
there (0.909–0.996) because its timings are smaller and therefore noisier relative to themselves.

---

## 5. The contour is 29% of the job, not 54% — and the biggest stage is the collider check

This is the part with no published counterpart, and it is the one that should change what anyone
optimises.

A paper reports the *contour*. A consumer pays for a **usable mesh**: the contour plus everything that
has to happen before a renderer or a physics engine will accept it. The one published comparison this
project has of the two puts contouring at 68 ms against halfedge construction at 58 ms — **54%**.

Measured here, seven fields at 33³ and 65³, `f64`, median of seven runs after two warmups (M-135):

| stage | mean share |
|---|---:|
| contour | **29.0%** |
| normals | 0.4% |
| weld | 25.5% |
| collider check | **45.0%** |

**Optimising the extractor alone therefore buys at most 1.41× on the whole job.** That is the ceiling,
and it is worth knowing before starting rather than after.

Three qualifications, all of which survive.

**Normals are effectively free (0.4%)** because every reference field here overrides `Sdf::gradient`.
The measured `AreaWeightedFaces` pass is the cost a consumer pays whose field has *no* analytic
gradient — a sampled volume, an imported scan — and even then it is noise.

**The collider figure times `collider::readiness`, which validates.** A physics engine building its own
structure is a different cost, and a shipping game may validate once rather than per chunk. But
welded-manifold-correctly-wound is the contract, and checking that contract is the single largest line
in the budget.

**There is no bar for upload, deliberately.** Its CPU half is a *move* — the Bevy sink writes straight
into the arrays a `Mesh` owns — so there is no pass to time, and its GPU half needs a device the core
crate must not have. A bar labelled "upload" measuring a `Vec` move would imply a cost that is not
there.

### And the single number is a fiction

The contour's share runs **13.1% to 74.3%** across the seven fields — a **5.7× spread on identical
code** (M-136):

| field | 33³ | 65³ |
|---|---:|---:|
| `sphere` | 13.1% | 19.2% |
| `box_exact` | 15.1% | 21.1% |
| `csg_difference` | 16.1% | 22.5% |
| `torus` | 16.4% | 25.5% |
| `gyroid` | 17.5% | 22.9% |
| `thin_plate` | 33.0% | 44.0% |
| `fbm_terrain` | **65.2%** | **74.3%** |

The variable is **how expensive the field is to sample**, not anything about the mesher. fBm noise
swamps the fixed post-processing; a cheap analytic field leaves the post-passes dominant.

So *"is it worth optimising the contour"* cannot be answered in general. **A game on procedural terrain
is in `fbm_terrain`'s regime and should optimise the extractor. One on authored CSG is in `sphere`'s
and should not.** The share also rises with resolution on every field, because the contour is `O(n³)`
while the post-passes scale with triangle count.

---

## What surprised us

1. **The crude baseline holds the good corner.** Marching Cubes being the only entry with both zero
   non-manifold edges and zero self-intersections was not the expected result.
2. **The collider check costs more than the contour.** Nothing in the literature suggests looking, and
   it is 45% of the job.
3. **The field decides the ratio, by 5.7×.** "What fraction is the extraction" is a property of the
   workload, not of the code.
4. **Marching Tetrahedra is more accurate than Marching Cubes on sharp fields** — the opposite of the
   registered prediction, and for a reason (more edge families) that the prediction did not consider.
5. **Surface Nets gets *worse* with resolution.** The folklore has it as the cheap default; its
   per-sample cost rises on both machines while Marching Cubes' falls.
6. **The sharp-feature solve is nearly free** (3–6.5%), and so is the clamp that makes it safe, and so
   is the asymptotic decider. Three "expensive" features that are not.

## What we would now distrust

Not individual facts — *sources*. Each of these published a single number for a quantity that is
field-dependent:

- The `2.76×` greedy-quads merge ratio. Measured range: **1.70× to 256×** (M-56).
- Nielson's `1.3%` manifold-split rate. Measured: **0 on five of seven fields** (M-60).
- The `54%` contour share. Measured: **13.1%–74.3%** (M-136).
- ✗14's "cheapest in the family". Measured: **slower on every field, on two machines** (M-45).

The pattern is consistent enough to be a rule: **a single figure quoted for a mesher's behaviour is
usually one scene's number.** Ask which field it was measured on before believing it, including of the
numbers in this document — every one of them names its field.

---

## Reproducing this

```bash
cargo bench --bench shootout          # topology + accuracy -> docs/measurements/shootout.csv
cargo bench --bench resolution_sweep  # 16^3..256^3 -> docs/measurements/resolution_sweep.csv
cargo bench --bench stage_breakdown   # stage shares -> docs/measurements/stage_breakdown.csv
cargo bench --bench extract           # per-algorithm regression timings
```

Always `--release`; `cargo bench` handles that. The sweep is a no-op under `cargo test` so a debug
build cannot overwrite the committed CSV.

The live version of §4's fit is a Bevy example, which is the only way to see a two-term model failing
rather than being told it fails:

```bash
cd bevy_isomesh && cargo run --example resolution_plot --release
```

---

*Ledger entries behind this document: M-20, M-21, M-23, M-25, M-28, M-42, M-45, M-52, M-53, M-54,
M-55, M-56, M-60, M-62, M-134, M-135, M-136, and ✗14. Tickets: M-001a, M-001b, M-002, M-003, M-004.*
