# Experiments that held, and what they buy you

Seventeen times, this project wrote down what it expected to measure — and what result would prove it
wrong — **before running the measurement**. Ten of those are enforced by the compiler.

This page is the scorecard. Five of the ten compiler-enforced predictions held, four were falsified,
one never ran, and the most useful thing in the whole record is a bug that made two true hypotheses
look false for a day.

If you only read one section, read [the instrument story](#when-a-measurement-is-impossible-suspect-the-instrument).

---

## Why pre-register at all

A prediction written after the number is known is not evidence, and it reads exactly like one written
before. This project has already caught itself writing expectations into docs that the measurement
then disproved — twice, at ✗1 and ✗3.

So the prediction moved somewhere the compiler can see it:

```rust
use isomesh::experiment;

// This line does not compile if "P-8" is not in `PREREGISTERED`.
let p = experiment!("P-8");
assert!(p.hypothesis.starts_with("A weld gated"));
```

Three properties make it a gate rather than a habit:

- **An unregistered id is a compile error**, not a test failure you can ignore.
- **`Preregistration` is `#[non_exhaustive]`**, so nothing outside the crate can construct one. The
  macro is the only door.
- **`falsified_by` is not optional.** There is no way to express a registration without it, which is
  the point of the field existing: *a hypothesis with no falsifier is not a hypothesis.*

Registering is a commit, so git carries the ordering that prose cannot. `scripts/backlog_gate.sh`
checks that every registered id reaches `FINDINGS.md`, so the prose can elaborate and cannot
contradict.

| era | how it was registered | held | falsified | never ran |
|---|---|---|---|---|
| P-1 … P-7 | prose in `FINDINGS.md` | 5 | 2 | — |
| P-8 … P-17 | `crates/isomesh/src/experiment.rs` | 5 | 4 | 1 |
| P-18 … P-21 | `crates/isomesh/src/experiment.rs` | 3 | 1 | — |

A page called "experiments that held" which quietly dropped the four falsifications would be exactly
the failure this machinery exists to prevent. They are all here.

---

## The scorecard

| prediction | verdict | what it means if you depend on this crate |
|---|---|---|
| **P-11** — one canonical `world_of_sample` closes seam cracks at every cell size, not just powers of two | **held**, both clauses (M-278) | Chunk seams are exact **only at power-of-two cell sizes** today. Canonical reconstruction gives 0 unmatched seam edges at every spacing tried; what the crate can offer now gives 0 only at powers of two and 63–348 at `0.1`, `1/12`, `1/14`. The weld hides it; an unwelded consumer gets it in full. The fix is priced and open as X-005 |
| **P-15** — more than half the dual mesher's cycles are in `emit_quads` | **held**, and by more than it claimed (M-284) | 45% of the instructions and **82% of the time**, at `r² = 0.995`, predicted to +0.12% on a held-out shape. This is the measurement that made the 4.26× possible |
| **P-9** — k-way weld order changes the output | **held** | The weld is **not confluent**, so vertex order is pinned deliberately rather than left to `HashMap` iteration. Your meshes reproduce byte-for-byte |
| **P-13** — the non-convergent angle at a sharp feature is bounded below by the feature's dihedral | **held** — after M-289 | It is a **property of sharp edges**, not a defect with a location. Nothing to fix; something to know |
| **P-16** — every past-90° vertex sits on a cell straddling the crease | **held, at 0%** — after M-289 | Same conclusion from the other side: two faces meeting inside one cell, not a winding bug |
| **P-8** — a weld gated on the link condition reaches zero non-manifold output | **falsified**, both clauses (E×4) | **The gate makes things worse.** The crate does not ship one, and now you know why rather than wondering |
| **P-12** — the dual's superlinearity is a cache-miss effect | **falsified** (M-279) … and later true | The best methodological lesson here. See below |
| **P-14** — the residual non-manifold vertices are one-vertex-per-cell meeting two surface components | **falsified** on clause one (M-275) | The vertex-to-edge ratio came out **exactly 2:1 on every row** — 597/314, 99/48, 38/19, 106/53. There are not 597 defects plus 314 defects; there are 314, each counted twice more by the vertex walk |
| **P-17** — Manifold Dual Contouring's residue is an interior ambiguity | **falsified** (M-291) | `Interior::Joined` fires on **100%** of ambiguous-face pairs, offenders and control alike — no discriminating power at all. Bounded exhaustively instead (M-292), then reduced to a 48-sample fixture (M-294) |
| **P-10** — vertex inflation from gate-plus-split is under 1% | **never ran** | Superseded. P-8 killed the gate P-10 presupposed before P-10 could measure its cost |
| **P-19** — the winding number's on-demand crossover is set by row sharing, so `Q*` is of order `N²` not `N³` | **held**, both arms (M-303) | `Q*/N²` came out **1.13–1.43** at 17³/33³/65³ against a registered 0.5×–4× band. The naive reading overstates the win by a factor of `N`: at 65³ the real crossover is **4,791 queries against 274,625 samples**, so an on-demand field pays below about **1.7% of the grid**. The hole-punched control moved `Q*/N²` to 13–42, and per-point cost tracked boundary-edge count to **×15.7 against ×16.2** — the mechanism confirmed by number, not by direction |
| **P-20** — a weld key moves no topology metric beyond the splits it names | **held on both falsifiers, wording amended** (M-305) | A constant key moves nothing on **0 of 51** configurations, and non-manifold vertices never rise by more than the split count — so E×4's manufactured bowties do not come back. But boundary edges go **0 → 24** on a creased cube, which H said would not happen. That is not a defect: it **is** the split, seen from the edge column. The exact mirror of E×4, where a *pairwise* refusal did its damage in the vertex column |
| **P-21** — a freshly extracted mesh separates exactly the sample pairs the field's own sign separates | **held**, both clauses (M-307) | **Marching Cubes, MC+decider and Marching Tetrahedra seal 24 of 24** — eight fields × three resolutions, 9.15 M probes, zero holes and zero membranes. **All three duals fail `fbm_terrain` with the identical count** (92 / 138 / 190) and **every one of those holes is on a face of the sampled domain**: a dual's quad needs all four cells around a sign-changing grid edge and at the domain face only one or two exist. **For a chunked world that face is the chunk seam**, so a dual chunk's collider is not watertight on its own. The property that splits the family is not primal-versus-dual but *does the method put its crossing on the probed grid edge* — subgrid Marching Tetrahedra is primal and scores worst, 18 of 24 |
| **P-18** — every convex-decomposition precondition is already reported by `ColliderReadiness` | **falsified**, and not where predicted (M-300) | Registered *expecting* to die on self-intersection-freedom. It died on **non-manifold vertices** instead, and self-intersection-freedom turned out not to be an input precondition of anything (✗20). The audit was the product either way |

---

## The 4.26×, and not one triangle changed

Two changes, three days apart, on the same function.

**A-023 made the loop axis a constant.** `DualMesher::emit_quads` took `axis`, `u` and `v` as runtime
values, so `p[axis] = a` was a dynamically indexed store. The coordinate array could not live in
registers, and every iteration wrote three coordinates to the stack that `linearize` read straight
back — a store-to-load forwarding chain. Splitting it into three `const`-generic monomorphisations
left the same three passes in the same order over the same bounds:

| `emit_quads`, per iteration | before | after |
|---|---|---|
| cycles | 43.26 | **3.33** |
| instructions | 31.24 | **12.02** |
| share of the mesher's cycles at 193³ | 82% | **27–32%** |

**A-024 forced the row length odd.** `values` was laid out by the caller's shape, so at 128 samples
per axis the row stride was **512 bytes** and the plane stride **exactly 64 KiB** — a cache-set
aliasing period on this machine, twice over. The remedy is `size[0] | 1`, unconditional and
idempotent. Unconditional because a pad applied only when the stride looks bad is a second layout
reachable from one call; idempotent because a *fixed* pad of one would map every `size[0] = 2ᵏ − 1`
onto the very stride it is avoiding.

Together, on a sphere at 256³, `f32`, one thread, one binary, one run:

| | before | after |
|---|---|---|
| Surface Nets | 693.8 ms | **162.7 ms** — 4.26× |
| Surface Nets IPC | 1.20 | **4.09** |
| `SN / MC` | 5.43× | **1.26×** |
| per-sample cost, 16³ → 256³ | 29.63 → 41.35 ns (+40%) | **8.71 → 9.70 ns (+11%)** |
| cache misses per sample | 3.72 | **1.56** |

**Not one triangle changed.** The golden hashes pass untouched — an optimisation that changes the mesh
is a bug in the optimisation.

And it cost this repository one of its own beliefs. ✗14 read *"Surface Nets never wins on Zen 3, at
any resolution: 2.46× behind even at 16³."* Surface Nets is now **faster** than Marching Cubes at 16³,
24³ and 32³:

| n | 16 | 24 | 32 | 48 | 64 | 96 | 128 | 192 | 256 |
|---|---|---|---|---|---|---|---|---|---|
| Marching Cubes ns/sample | 9.37 | 8.87 | 8.29 | 8.03 | 7.74 | 7.41 | 7.48 | 7.62 | 7.67 |
| Surface Nets ns/sample | **8.71** | **8.59** | **8.05** | 8.24 | 8.42 | 8.53 | 8.75 | 9.53 | 9.70 |
| ratio | **0.93** | **0.97** | **0.97** | 1.03 | 1.09 | 1.15 | 1.17 | 1.25 | 1.26 |

The reason it had been true was one missing `const`.

---

## One bit of a row length, at the chunk size everybody uses

128³ is the canonical chunk size in voxel engines, and it was the worst possible one.

The diagnosis came before any code changed, by letting the *caller* arrange the pad — adding one
sample moves the stride while changing the work by under 1%:

| shape | plane bytes | cycles/sample | vs 127³ and 129³ |
|---|---|---|---|
| 127³ | 64,516 | 33.10 | — |
| **128³** | **65,536 = 2¹⁶** | **108.51** | **3.37×** |
| 129³ | 66,564 | 31.39 | — |
| 129×128×128 — pad `x` | 66,048 | 31.48 | 0.98× |
| 128×129×128 — pad `y` | 66,048 | 36.45 | 1.13× |
| **128×128×129 — pad `z`, the control** | **65,536 = 2¹⁶** | **107.89** | **3.35×** |

**The `z` row is what makes this a measurement rather than a story.** It adds the same 0.8% of work and
touches neither stride, and it keeps the entire penalty. Without it, "adding a sample helped" would be
indistinguishable from "the fixture got slightly different."

Measured on a field with no surface in it at all — a sphere entirely inside the domain — so this is the
scaffolding, not the geometry.

---

## When a measurement is impossible, suspect the instrument

This is the most useful entry in the file, and it is an entry about being wrong.

R-006 and R-008 both compared vertex normals against a *reference* gradient computed from an exact
wedge SDF written for the experiment. Both came back falsified. P-16 missed its threshold by four
times over. Two independent hypotheses about sharp features, both dead, on the same afternoon.

The reference gradient was returning noise.

Its exterior branch normalised `away = q − dir·t` guarded only by `e > 0.0`. For a point lying *on* the
ray, that difference is a catastrophic-cancellation residue of order `ε·|q|` — a random direction with
a plausible magnitude. About half of all Marching Cubes vertices land epsilon-outside the wedge, so
about half were being compared against a random unit vector.

| | as reported | corrected |
|---|---|---|
| past-90° vertices, Marching Cubes | 6,959 | **472** |
| past-90° vertices, Dual Contouring | 4,868 | **232** |
| rows containing one, of 168 | 75 | **8** |
| **P-13** | falsified | **held** |
| **P-16** | falsified at 4× its threshold | **held, at 0%** |

Both falsified entries are still in `FINDINGS.md`, with their original wrong text, under a banner
pointing at the correction. Deleting them would delete the evidence that the instrument can be the
thing that is broken.

The rule it earned, now in the method section: **a reference implementation used as ground truth needs
the same scrutiny as the thing it checks.** And its shorter form — *when a measurement is impossible,
suspect the instrument before the world.*

P-12 is the same lesson wearing different clothes. It predicted the dual's superlinearity was a cache
effect; the miss rates came back flat and it was recorded falsified. That was **a correct measurement
of a machine whose bottleneck was somewhere else.** Once A-023 removed the store-to-load stall that was
covering it, the 128³ spike went from 1.24× — indistinguishable from noise — to 2.6×, and A-024 then
took it to 1.01×. A null measured under a dominant confound is a statement about the confound.

---

## The paper is wrong, and here is the fixture

Schaefer, Ju & Warren, *Manifold Dual Contouring* (IEEE TVCG 13(3), 2007), §3, on the uniform grid:

> *"this surface is always a manifold because the original MC algorithm always constructs a manifold
> and the dual preserves the topology of the surface."*

Over eight reference fields at three resolutions, no chunking and no weld, so there is only one
mechanism in the count:

| | non-manifold edges | non-manifold vertices |
|---|---|---|
| Marching Cubes | **0** | **0** |
| Marching Cubes + asymptotic decider | **0** | **0** |
| Manifold Dual Contouring, this crate's default table | 143 | — |
| Manifold Dual Contouring, the paper's decider-modified table | **114** | 222 |

The premise holds here. The conclusion does not.

Every one of those defects is on `noise_cavity`, and the mechanism now fits in 48 samples. A `2×2×3`
column of `±1` inside a `4×4×3` lattice, trilinear between the samples, everything else outside — its
middle plane reads `out, in, in, out`, the face saddle. On those twelve hand-written values Marching
Cubes emits 40 triangles and **0** non-manifold edges; both duals emit 20 triangles and **1**, carrying
four distinct faces.

Four quads meet on one dual edge because all four cut edges of the shared face land in one cycle on
each side, so all four connect the same pair of vertices. And the manifold construction buys nothing
here: splitting a cell by cycle cannot split a cell that has one cycle.

The sharper half is that **the offending set is not a set of sign configurations.** Scale the shared
face's two inside corners — no sign moves — and the asymptotic decider's saddle crosses zero with the
defect:

| inside corners | saddle `s` | the face | non-manifold edges | triangles |
|---|---|---|---|---|
| `−0.25` | `+0.375` | separated | **1** | 20 |
| `−1` | `0`, an exact tie | separated | **1** | 20 |
| `−4` | `−1.5` | **joined** | **0** | 20 |

One configuration, two answers, chosen by the field rather than by the rule. Which is why a face rule
still leaves 25–49 offending faces per resolution while being combinatorially capable of avoiding all
of them — over all 4,096 two-cell sign patterns, **not one** offends under every consistent mask.

---

## Reading the tiers

Every number in this repository's documentation carries one of these, and the legend matters because
the tiers are not decoration:

| tier | means | the bar |
|---|---|---|
| **M** | measured here | We ran it. Code and numbers in this repo, reproducible by checkout |
| **V** | verified externally | We read the primary source ourselves. DOI or file attached |
| **R** | reported | A credible source asserts it; we have not checked |
| **F** | folklore | Widely repeated, no verified primary source found |
| **✗** | falsified | Tested and found false — including what we believed, and why |

**Never cite an R or F claim as justification for a design decision without saying which tier it is.**

Two standing rules keep the ledger honest. **Re-tier rather than rewrite** — when an R becomes an M,
the R text stays and the measurement goes below it, because the gap between what was reported and what
was measured is itself data. And **record the ones we got right for the wrong reason**: ✗10 sits in the
falsified section even though its outcome was fine, because the reasoning was wrong and the reasoning
is what generalises.

---

## Where the raw results live

| | |
|---|---|
| the registrations | [`crates/isomesh/src/experiment.rs`](../crates/isomesh/src/experiment.rs) |
| the full ledger, 383 entries | [`FINDINGS.md`](../FINDINGS.md) |
| per-experiment result data | [`docs/experiments/`](experiments/) — each CSV carries its `# hypothesis:` and `# falsified by:` header |
| timing and quality measurements | [`docs/measurements/`](measurements/) |
| what falsified beliefs cost | `FINDINGS.md` Part 1 — 24 entries, each recording where the wrong belief came from, because provenance predicts the next error |
