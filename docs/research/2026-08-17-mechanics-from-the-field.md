# Mechanics from the field — a novelty dossier

**Date:** 2026-08-17
**Question asked:** what will make games more interesting, sourced from academia, that nobody has ever
shipped.
**Method:** five parallel corpus hunts (surface-intrinsic computation; modal analysis; structural
mechanics; shape semantics; volumetric processes), each authorised to download, convert and index
missing papers. 20 papers acquired this session. Every candidate had to answer **"name the shipped
game that does this"** and say so when the answer was a game.

**Evidence grading, preserved from the hunts.** `READ` = the passage was read in the converted
document. `ABSTRACT` = metadata only. `SNIPPET` = quoted inside another paper. Every DOI and arXiv ID
below was returned by a tool; none was reconstructed. Where a hunt's number is its own arithmetic on a
paper's table rather than a figure the paper states, it says so.

**One warning that applies throughout.** This document is tier R. It is five agents' reading, not this
project's measurement. Where a claim touches something `FINDINGS.md` measured, the M-row wins — and
§4.1 is a case where an M-row's *interpretation* turns out to be wrong in a way that opens a line
rather than closing one.

---

# Part 0 — Two laws that decide which of these ship

The five hunts ran independently and converged on the same two structural facts. They are worth more
than any individual candidate, because they let you triage the next hundred ideas without a literature
search.

## Law 1 — A process is expressible in the edit log exactly when its state is sparse

The volumetric hunt states it outright: *a process is log-expressible exactly when its natural state
variable is sparse.* Dissolution's state is one aperture per conduit edge. A speleothem's state is one
length scalar. Both are O(sparse), both emit same-kind hard brushes, and both therefore inherit
**M-36's measured commutativity for free** — one result from all 40,320 orderings of eight brushes.

Generalised across all five hunts, it sorts every candidate in this document:

| Candidate | Natural state | Sparse? | Verdict |
|---|---|---|---|
| Speleogenesis | aperture per conduit edge | yes | log-expressible, exactly |
| Speleothems | length per speleothem | yes | log-expressible, exactly |
| λ-medial axis / Calibre | **none** — pointwise function of `ρ`, `∇ρ` | free | no state at all |
| Structural feasibility | force per interface, O(blocks) | yes | sparse LP, 12 nonzeros/column |
| Surface reaction–diffusion | 2 floats per surface vertex | O(n²), not O(n³) | survives, keyed by grid edge |
| Modal sound / fracture | k modes per chunk | yes | sparse, but the *solve* is the problem |
| Air acoustics | one reverb per air component | yes | already built |
| Hydraulic erosion, repose relaxation | dense field | **no** | forces a grid and a bake |

**And the cost of the "no" row is precisely stated.** A dense-state process must materialise and mutate
a sample grid, then *bake* — collapse base + log into a stored field and resume logging. That costs
four things at once: the analytic base is gone (memory goes O(edits) → O(volume)); commutativity is
lost across the bake boundary, because a bake is not a min/max brush; multiplayer can no longer ship
the log, so late joiners download a baked field; and re-mesh degrades from edit-proportional to
region-proportional, which is M-314/M-318's whole result thrown away.

**So: choose processes by the dimensionality of their state, not by how physical they feel.** The
intuition inverts — the most cinematic candidate (a cave dissolving itself into existence) is exactly
log-expressible, and the one that feels most "voxel game" (rubble slumping to its repose angle) is the
one that breaks the representation.

## Law 2 — Classify by whether a local edit has a *local answer*, not by whether the operator is cheap

Three hunts found the same wall from three directions and none of them was looking for it.

- **Surface-intrinsic:** the heat method is 21–68× faster than fast marching *per query* — and its
  precompute is 0.21 s at 28k triangles rising to **63.4 s at 1.6M** (Crane, Weischedel, Wardetzky,
  Table I, `READ`). Read the precompute column as **the carve cost**. A prefactored operator is
  disqualified for a live carve unless the factorisation can be updated rather than rebuilt.
- **Structural:** feasibility is a global property of a connected structure. There is no Saint-Venant
  decay for a feasibility question. Removing a keystone changes the admissible force state of the
  entire arch, so the re-solve scope is *the structural component containing the edit* — which, if the
  player has built one continuous cathedral, is the cathedral.
- **Volumetric:** the explicit-surface CFL — *"the time step is chosen such that the front moves less
  than one grid point"* (Peng et al., `READ`) — and, independently, a 2024 heightfield-erosion paper
  reporting *"the convergence speed is limited by the fact that basin boundaries can only move by one
  cell per iteration"* (`READ`). **Two literatures, twenty-five years apart, discovering the same
  constant.** One cell per step is the real cost law of any explicit surface process.

The one place this wall has a door, and it is the most interesting technical finding in the sweep:

> **Combined Approximations reanalysis converges iff `‖K̄⁻¹ΔK‖₂ < 1`** — Dahlberg, Dalklint, Spicer,
> Amir, Wallin 2023, DOI `10.1007/s00158-023-03616-7`, **downloaded, converted, indexed, `READ`**.

A brush edit that changes 15–36% of cells in a small bounding box (**M-33**) gives a small `ΔK` and the
series converges in a handful of terms. A brush edit that removes a load-bearing pillar gives a large
`ΔK` and **the series diverges — and the divergence is the collapse signal.** The convergence test of
the numerical method *is* the gameplay predicate, at roughly 1% of the cost of the full solve. Nothing
in the games literature has noticed this.

---

# Part 1 — Tier 1: build these

Four candidates that are novel, evidenced, and whose incremental answer is local.

---

## 1.1 CALIBRE — width is a resource, and the world adjudicates it

**The mechanic.** Every creature and projectile carries one number: `λ`, its half-width. The engine
maintains one scalar field over the void. For any `λ`, the space reachable by that agent is the
connected component of `{r ≥ λ}` containing it — one structure, all agent sizes, no rebuild. In play:
you cut a 0.4 m rathole to escape a 0.9 m brute. It cannot follow — not because a designer tagged the
tunnel, but because `{r ≥ 0.45}` is connected across it and `{r ≥ 0.9}` is not. The brute widens the
hole, and **the instant the widest inscribed ball at the throat passes 0.9 m the components merge and
it comes through.** Dig narrow to survive, dig wide to pursue, and the door is a measured length.

**Why nobody has it.** Recast/Detour — the industry standard — voxelises, classifies walkable spans,
**erodes by agent radius at build time**, and contours. Agent radius is *baked*; N agent sizes = N
navmeshes = N full builds. It is 2.5D and derives walkability and nothing else. The genuinely close
academic prior art is the **Explicit Corridor Map** (van Toll, Cook, van Kreveld, Geraerts, DOI
`10.1145/3204456`, `ABSTRACT`), a medial-axis navmesh that explicitly *"enables path planning for
disk-shaped characters of any radius"* from one structure, with local dynamic updates (DOI
`10.1002/cav.1468`, `ABSTRACT`). But the ECM's medial axis is defined *"based on projected distances on
the ground plane"*, for characters *"constrained to walkable surfaces"*, from a **polygonal** boundary,
with obstacles as polygons that appear and disappear. It is a 2.5D floor-plan skeleton. **The
volumetric λ-medial axis of the void, maintained incrementally from a live SDF, is unoccupied.**

**The math — and it is a derivation, not a citation.** From Attali, Boissonnat & Edelsbrunner,
*Stability and Computation of Medial Axes: a State-of-the-Art Report*, DOI `10.1007/b106657_6` (in
corpus, **`READ` in full**): Eq. (1) gives `∇ρ(x) = (x − c(x))/ρ(x)` where `c(x)` is the centre of the
smallest ball enclosing the closest-point set `Π(x)`; and the λ-medial axis is
`M_λ[X] = {x : r(x) ≥ λ}` where `r(x)` is that ball's radius (Chazal & Lieutier). Every point of `Π(x)`
lies at distance exactly `ρ(x)` from `x`, so `Π(x) ⊂ ∂B(x,ρ)`, `x` sits on the normal axis through `c`,
and Pythagoras gives:

```
r(x) = ρ(x) · sqrt( 1 − ‖∇ρ(x)‖² )        exact, any dimension, any |Π(x)|
```

**λ-medial-axis membership is a pointwise O(1) test on the SDF value and its gradient.** No Voronoi
diagram, no Delaunay, no poles, no global structure. Both inputs are already computed per sample —
`ρ` is the field, `∇ρ` is the shading normal. Cross-checked against the two-point case printed in the
survey: `Π = {a,b}`, separation angle `θ`, gives `r = h`, `‖∇ρ‖ = cos(θ/2)`, `ρ = h/sin(θ/2)` — exactly
the relation the survey states.

**Prior art, stated honestly.** Separation-angle pruning of a distance field with its gradient is
23-year-old work: Foskey, Lin & Manocha's θ-SMA, DOI `10.1115/1.1631582` (`ABSTRACT`), *"relies on
computation of the distance field and its gradient using graphics hardware"*, seconds-to-minutes on a
2003 machine, **from a static field**. The novelty is not the estimator. It is that (a) it is
incremental under destruction and (b) **its parameter is a gameplay quantity**, so the estimator
becomes a mechanic rather than a shape-analysis tool.

**The reframe that makes it work, and it is the best idea in this document.** The 3D Skeletons STAR
(Tagliasacchi et al., DOI `10.1111/cgf.12865`, in corpus, §2.4.7 `READ`) lists two complaints about the
λ-medial axis: thresholding local measures can **disconnect skeletons** ("reconnection needs extra
work"), and *"this metric does not allow capturing details at different scales, as salient shape
features are removed before noise."* In a destructible game where `λ` is an agent's half-width,
**disconnection is the correct answer** ("that agent cannot get through") and deleting sub-λ features
is **the correct answer** ("that gap is not a passage"). *The academic defect is the game's semantics.*

**And it answers the objection that closed the Morse–Smale line.** That line was closed because *"the
persistence threshold is a hand-tuned magic constant, and for a game that is fatal, since the threshold
decides what counts as a room."* Here the threshold is **λ, a length in metres, equal to a creature's
shoulder width.** It is not tuned; it is measured off the character.

**Cost, incremental question first.** `r(x)` is free per sample. What is not local is *connectivity* —
but digging is monotone for the void, so `{r ≥ λ}` only ever grows and components only ever **merge**
→ union-find, near-O(1) amortised, no deletions. **The cheap direction is the common verb.** Filling
splits components and needs the decremental structure R-022 already owns. The genuinely hard
sub-problem is maintaining `ρ` itself under material removal; the nearest off-the-shelf answer is the
limited incremental distance transform (Scherer, Ferguson & Singh, DOI `10.1109/robot.2009.5152790`,
`ABSTRACT` — *"computation time reduced by an order of magnitude"* versus a full EDT).

**Pre-registered predictions.**
- On a 64³ grid at 5 cm voxels with ≥3 overlapping brushes, `ρ√(1−‖∇ρ‖²)` matches a brute-force
  minimal-enclosing-ball radius to within **1 voxel for ≥99% of samples**, and the median residual
  **halves when the voxel size halves** (O(h)). If the residual is h-independent, the discrete gradient
  has broken the identity and every candidate in §1.1, §2.4 and §3.3 dies with it.
- For a 0.5 m dig brush at 5 cm voxels: λ-membership flips ≤ **4×** the samples whose `ρ` changed;
  merge events ≤ **8** per edit on average; total maintenance **< 0.3 ms** on one core.
- Across 20 hand-built scenes, `{r ≥ λ}` component count changes at a doorway **only** when measured
  clearance crosses λ, hysteresis **< 1 voxel**, with **no parameter but λ**.

**Harness.** A brute-force oracle: per sample collect `Π(x)` by exhaustive search within `ρ+tol`, solve
the minimal enclosing ball (Welzl), compare against the closed form at three voxel sizes. Then an
edit-replay harness on the brush log instrumented for dirty samples, membership flips, merge events and
wall-clock.

---

## 1.2 SPELEOGENESIS — the world digs its own caves, along the route you chose

**The mechanic.** The world carries a hidden **proto-conduit graph**: a lattice of hairline fractures
with sub-millimetre apertures, generated with the base field, invisible and inert. Water enters at
sinks and leaves at springs. Nothing happens.

Then the player breaches an aquifer, or cuts a channel that redirects surface water into a fracture, or
dams a stream. **Flow reorganises. The conduits carrying the most flow dissolve fastest, and because
widening lowers resistance, a widening conduit steals flow from its neighbours.** Positive feedback.
One pathway wins, runs away, and breaks through: hairline crack → crawlspace → passage → chamber. The
losers stop growing and stay hairline forever. Dam it again and the system re-competes, the old passage
goes dry, a new one starts elsewhere. **You are doing landscape engineering, not level design.**

**Why nobody has it.** Three reasons compound. A heightfield cannot represent a conduit, and
heightfields are ~100% of shipped terrain-erosion tech — the entire industry toolchain is *structurally
incapable* of this. The relevant literature is in hydrogeology journals (HESS, Water Resources
Research), so a graphics programmer searching "erosion" never finds it. And everyone assumes a
self-modifying voxel world means mutating voxels, which is expensive and non-deterministic, so it is
rejected at the design stage. **The graph formulation dissolves all three at once.** *From Dust* runs
genuine real-time erosion on 2.5D material layers and cannot make a cave, an arch or an overhang.
Nothing has shipped flow-coupled dissolution.

**The math.** Dreybrodt & Gabrovšek, *Dynamics of wormhole formation in fractured limestones*, HESS
2019, DOI `10.5194/hess-23-1995-2019` — **downloaded**, conversion failed, `ABSTRACT` (unusually
detailed, carries the mechanism). A 2D net of 1D fractures, constant-head boundary conditions, linear
dissolution kinetics. Two phases: uniform widening, then **an instability** — *"due to small
perturbations, some of the foremost fractures gain length compared to the neighboring ones… they
attract more fresh aggressive water and their propagation is enhanced."* And the gameplay-critical
sentence: *"Several wormholes (caves) are penetrating into the aquifer but only one reaches the output,
whereas the others stop growing due to the redistribution of hydraulic heads caused by the leading
wormhole… there is a critical distance between the wormholes. Within this distance only one wormhole
survives."* Plus the tuning knob, stated by the authors: *"if one uses a heterogeneous net, the first
step of evolution is suppressed because of the large perturbations, and wormholes start to grow
immediately."* Seed heterogeneous apertures and you skip the boring uniform phase.

Supporting, also downloaded: Perne, Covington & Gabrovšek, DOI `10.5194/hess-18-4617-2014`; Covington &
Perne, DOI `10.3986/ac.v44i3.1925`.

**And the representation is already published — by a graphics group.** Paris, Guérin, Peytavie, Collon
& Galin, *Synthesizing Geologically Coherent Cave Networks*, CGF 2021, DOI `10.1111/cgf.14420`
(`ABSTRACT`; the download resolved to a Wiley landing page): karst skeletons from a gridless anisotropic
shortest-path, and then — *"From this skeleton, we define the geometry of the conduits as a **signed
distance function construction tree** combining primitives with blending and warping operators."*
**They already have your representation. They generate statically. Driving it dynamically from flow is
the entire gap, and it is one sentence wide.**

**Cost — and this is the candidate that costs almost nothing.** State is `E` scalars, one aperture per
edge. Per tick: solve steady flow on the graph (sparse SPD, E = 10³–10⁴, **sub-ms to ~1 ms**); update
apertures, O(E), trivially parallel; rewrite the capsule brush for each edge that crossed a voxel
threshold. **Coaxial nested subtractive capsules collapse under `max` — only the largest matters — so a
conduit widening for 10,000 ticks rewrites one brush rather than appending 10,000. Log size stays
O(E), not O(E × ticks).** Determinism: fix the CG iteration count rather than converging to a
tolerance, fix the reduction order, and it is bit-reproducible. Multiplayer replicates a few KB of
changed scalars. **Nothing in the representation breaks** — a dissolution capsule and a player dig are
the same kind of subtractive brush, so M-36's commutativity covers both.

**Pre-registered predictions.** At E = 4096: a full tick (32 fixed CG iterations + aperture update +
brush emission) in **< 2 ms single-threaded**; fewer than **5% of edges** cross a voxel threshold per
tick at a rate tuned to visible change every 20 s; **post-breakthrough, >90% of dissolution flux in
<10% of edges**, so per-tick re-mesh cost *decreases* monotonically after breakthrough — if it
increases, the kinetics is wrong; log size **O(E)**, and if it grows linearly in ticks the coaxial
collapse is not firing.

**Harness — and it needs no renderer, no field and no engine.** A headless graph simulator: 64×64
lattice (E ≈ 8000), constant head in/out, log-normal initial apertures, linear kinetics. Plot the
aperture distribution over time — the literature predicts it goes **bimodal** under competitive flow
and stays **unimodal** when recharge is limited, so reproducing that bifurcation validates your
kinetics against published geomorphology before a single voxel is touched. A few hundred lines.

---

## 1.3 THE ROOM YOU DUG — acoustics of the air volume, on the frame you break through

**The mechanic.** Dig a tunnel; it *sounds* like a tunnel, immediately. Breach into a cavern and the
reverb tail opens up on the frame of the breach. Seal a doorway and the room goes dead. Your footsteps
tell you the size of a space you cannot see. Wall a chamber in two and each half gets its own reverb.

**Why nobody has it.** The dominant shipped technology is a **static bake**. Microsoft Project
Acoustics, shipped in *Gears of War 4*, describes its own philosophy in its documentation (`SNIPPET`):
*"similar to static lighting: bake detailed physics offline… and use a lightweight runtime."* Static
lighting is exactly the right analogy and exactly the limitation — a destructible world invalidates the
bake on every dig. The underlying Raghuvanshi & Snyder line is explicitly for *static* scenes.

**The escape hatch, and it has public source.** Rosen, Godin & Raghuvanshi, *Interactive sound
propagation for dynamic scenes using 2D wave simulation*, DOI `10.1111/cgf.14099` (`ABSTRACT`,
paywalled): 2D wave simulation capturing geometry-based diffraction, *"real-time performance on a
single CPU core"*, on *"static scene snapshots"*, results *"robust to motion and geometry changes"* —
and *"We share the complete C++ code of our 'Planeverb' system."*

**Why this project specifically should build it first.** R-022's incremental air connectivity is already
measured, and it hands you the three things a dynamic-acoustics system needs and no shipped engine has.
**The reverb domain is the connected air component** — not a designer-placed reverb volume; the
component *is* the room. **The invalidation trigger is exact and cheap** — M-319..M-323 measured that
one fill in six splits the air region, so a split is a re-bake of two rooms and everything else is a
parameter nudge; you are not re-baking per edit, only on topology change. **The median split sheds one
voxel**, so the overwhelming majority of splits are trivial slivers you discard without simulating.
And component volume and surface area give a Sabine RT60 instantly, as a stopgap while a snapshot
re-bake runs asynchronously.

**Pre-registered predictions.** Breaking through changes the player's air component ID on the frame of
the breach, and a Sabine RT60 from that component is available in **< 0.1 ms** (two accumulators you
already maintain). A Planeverb-style 2D re-bake of a 64×64 slice completes in **< 30 ms** on one core,
amortisable over two frames. Falsified if either exceeds 3×.

**Harness.** Instrument air-component volume, surface area and split/merge events against the existing
tracker; compute Sabine RT60 per component; A/B against a 2D FDTD re-bake. **This is the cheapest ship
in the document: no eigensolver, a public reference implementation, and it monetises infrastructure
that is already built and measured.**

---

## 1.4 THE SAFE-TO-DIG FIELD — a second scalar field, on the same grid, that says what will fall

**The mechanic.** Before the player cuts, the game paints on the rock face **how much closer to
collapse removing this cell moves the structure**. Not a binary safe/unsafe — a scalar field, in the
same grid, same data type, same trilinear interpolation as the SDF.

**The math.** This is the **topological derivative**: the sensitivity of a shape functional to
nucleating an infinitesimal hole. Lopes, dos Santos & Novotny, DOI `10.1590/1679-78251252`
(**downloaded**, conversion failing, `ABSTRACT`) — *"measures the sensitivity of a shape functional with
respect to an infinitesimal singular domain perturbation, such as the insertion of holes"*, combined
*"with a level-set domain representation method."* Bonnet & Delgado, DOI `10.1093/qjmam/hbt018`
(`ABSTRACT`): closed-form topological derivative of compliance via **the adjoint solution approach**,
*"simple to implement and computationally efficient"*. Wang, Yang & Kang, DOI `10.1115/1.4053989`
(`ABSTRACT`): the field is evaluated everywhere on the level-set grid every iteration — so per-cell
evaluation is known-cheap.

**The direct structural analogue already exists and stopped one step short.** Whiting et al. 2012,
*Structural optimization of 3D masonry buildings*, DOI `10.1145/2366145.2366178` (`ABSTRACT` only,
download failed): *"We define a new measure of structural soundness… and derive its closed-form
derivative with respect to the displacement of all the vertices describing the geometry… The gradients
are visualized as interaction tools, giving user-guidance for effectively modifying a structural
design."* **They built the tool. For architects, on vertices, offline.** The move is vertices → SDF
cells, offline → after every brush stroke, architect → player.

**Why nobody has it.** No game exposes a structural sensitivity field. The blocker is **conceptual, not
computational** — nobody in games has connected "topological derivative" to "which voxel is safe to
mine", because the two communities do not overlap.

**Cost.** One adjoint solve per structural component per edit — a pair of triangular solves against a
factorisation you already have, not a refactorisation. Dahlberg et al. (`READ`) measured triangular
solves as cheap relative to factorisations (53k–63k triangular solves against 121–400 factorisations,
with factorisation dominating total effort). Under Law 2's convergence test, if `‖K̄⁻¹ΔK‖₂ < 1` you
reuse the old factor and the field update costs triangular solves only.

**What the SDF substrate buys — and I looked for an impedance mismatch and did not find one.** The
topological derivative is a scalar field on a fixed Eulerian grid. So is the SDF. Same grid, same
chunking, same LOD, same incremental invalidation, same GPU upload path, same trilinear sample in the
shader. **You can pack it into a second channel and shade the mined surface by it with a one-line
change.**

**Pre-registered prediction.** Given a warm factorisation, computing the field over a 32³ chunk costs
**< 1.5× the primal solve** and is dominated by triangular solves. And the field's **top-1% cells
contain ≥80% of the cells whose individual removal actually flips feasibility**, verified by brute-force
ablation on a 1000-block test structure. **If the top-1% capture rate is below 50%, the linearisation
is too weak to be a game signal and this candidate dies** while §2.1 survives.

---

# Part 2 — Tier 2: strong, with a named obstacle

## 2.1 THE ADMISSIBILITY GATE — the ceiling falls because it was load-bearing

**The finding that reframes this whole direction: the cheap check is not FEM.** It is a **lower-bound
limit-analysis feasibility program** — an LP/QP with **no material stiffness in it at all**, whose only
physical parameters are density and friction coefficient. Pepe et al. 2020, DOI
`10.3389/fbuil.2020.00043` (**downloaded, converted, indexed, `READ`**): *"A specific advantage of the
Limit Analysis approach is the low number of parameters required… friction and self-weight per unit
volume are the only pieces of mechanical information needed."*

**A rule with two parameters is a rule a player can build a mental model of.** That is a rare property
and it is the answer to the legibility objection. Compare an XPBD or mass-spring substitute, whose
effective stiffness depends on iteration count and timestep — unlearnable.

**And elasticity is not merely expensive here, it is the wrong question.** Whiting, Ochsendorf & Durand
2009, DOI `10.1145/1618452.1618458` (**downloaded, converted, indexed, `READ` in full**), §1: *"current
engineering tools based on finite element methods and elasticity theory are not appropriate in this
context because they focus on material failure and stress, and because the high stiffness of stone
results in poorly conditioned numerical systems… Block et al. [2006] demonstrated that linear elastic
theory was unable to differentiate between a feasible masonry arch and an infeasible arch."*

The program (§3.1, `READ`): find internal forces `f` satisfying equilibrium `A_eq·f = −w` (6 rows per
block), a linearised friction pyramid `A_fr·f ≤ 0`, and compression-only `f_n ≥ 0`. Their contribution
is a *continuous* infeasibility measure: split `f_n = f_n⁺ − f_n⁻` and minimise `Σ(f_n⁻)²`. **`y = 0`
means it stands (Heyman's safe theorem); `y > 0` is how much imaginary glue it needs, and the solver
tells you which interfaces need tension — those are the hinge lines.** The ceiling comes down along
hinges the mathematics produced. Then re-solve on the debris: that is literally the alternate-load-path
/ column-removal method of progressive-collapse engineering.

The paper contains the whole thing already — including Figure 14(b), *"the Sainte Chapelle model
collapses after a central column is broken."* **What it does not do is the gap: the analysis runs
offline at authoring time to pick grammar parameters, and the runtime collapse is handed to Bullet,
which has no idea whether the structure was admissible.** Nobody has moved the analysis to after the
player's edit and used its verdict as the rule.

**Calibration you can hold it to.** §3.3 (`READ`): minimum thickness/radius of a semicircular arch —
Milankovitch's 1907 analytic value **0.1075**, their solver **0.10746**. Critical ground-tilt angle at
t/r = 0.20 — Ochsendorf's **15.84°**, matched exactly. **Two independent external ground truths. This
is not a heuristic.**

**What every shipped game does instead, precisely.** PhysX Destruction/Blast does **island detection**
on a convex compound — *"island detection is simple and fast because the visual mesh provides
connectivity information"* (`READ`) — i.e. pure graph connectivity, so a one-voxel thread of rock holds
up a cathedral. Teardown, Space Engineers, Dwarf Fortress: connectivity plus radius. *Red Faction:
Guerrilla* is the closest shipped thing and genuinely propagates load — but over a **pre-authored**
element graph with scalar joint strengths, so it has no notion of "an admissible force state exists"
and cannot reproduce a hinging arch. Poly Bridge and Besiege solve real member statics, but on 1-D
members assembled from a parts palette. **No shipped game solves a feasibility program on
runtime-generated geometry.**

**Cost, and the honest incremental answer.** From Whiting's Table 1 and their stated cost model — this
arithmetic is the hunt's, not the paper's — one QP solve is ~3.1–4.4 s at 486 blocks and ~8.8–10.7 s at
986 blocks, on 2009 MATLAB + BPMPD, scaling ≈ O(n^1.7). **One voxel = one block is dead on arrival**
(a 32³ chunk is 33× Cluny). The viable shape is **macro-blocks**: O(200–2000) per structural island,
which puts a modern sparse solver at an estimated **30–200 ms** — a background-thread budget, not a
frame budget, which is *also better game feel*, because the ceiling should groan before it falls.

Three warnings the hunt surfaced and that must not be skipped. **The lower-bound theorem is one-sided**
— it gives existence of an admissible force state, not dynamic stability; in a game that is acceptable
(false negatives are things that stand which maybe shouldn't) but you must know it. **The block
decomposition is a hidden gameplay parameter** — Whiting measured that coarser blocks *over-estimate*
stability, and that rotating the friction pyramid 45° moved optimal column thickness by 1%–10.5%; for
multiplayer the decomposition must be a pure function of the edit log. And **warm-starting is a solver-
class decision**: Whiting used an interior-point solver, which is the worst class for warm starts. If
you want incrementality you must use dual simplex with a hot basis, and even then there is no
theoretical bound.

**Pre-registered predictions.** A Rust reimplementation with a modern sparse solver on 1000 macro-blocks
solves in **< 250 ms single-threaded**; a 100-block semicircular arch reports the infeasibility
threshold at **0.1075 ± 0.0010**; coarsening blocks 2× shifts the reported critical pillar thickness by
**≥ 5%, in the unconservative direction**; and warm-started re-solve is **bimodal** — under 15% of cold
time for edits that do not sever a load path, under 2× speedup for edits that do. *The bimodality is
the prediction; a unimodal result falsifies the cheap-incremental story.*

**The legibility harness, which is not optional.** Show 20 subjects a pillar and ask "if I cut here,
does it fall?" **If they are below ~70% correct, the mechanic is correct and unshippable** — and the
fix is §1.4's field plus rendering the thrust network diegetically, as veins of strain following the
compression path, so the player *sees* which rock not to cut. Panozzo, Block & Sorkine-Hornung (in
corpus, `READ`) give the form: `Cᵀ D_w C` — **the equilibrium operator of a compression-only structure
is a weighted graph Laplacian**, and `w ≥ 0` is exactly what makes it a proper one.

## 2.2 THE FIREBREAK — spread whose metric is geodesic distance on the carved surface

**The mechanic.** A blight crawls across the cavern skin at a fixed speed *measured along the surface*,
so it climbs walls, crosses ceilings and pours into your tunnel. You stop it by cutting a trench — and
you can see whether the trench is wide enough, because the front slows and detours; sever the last land
bridge and it stops dead though the far side is two metres away in straight-line space.

**Why nobody has it.** Games measure spread by Euclidean radius, by a navmesh (floor-only, 2.5D, does
not exist on walls or ceilings), or by voxel-adjacency BFS (Minecraft fire) — a *volumetric* graph
distance that cannot tell "around a 2 m ledge" from "through a 2 cm wall". True geodesic distance has
been offline: fast marching needs non-obtuse triangulations *"notoriously difficult to obtain"* and has
an inherently serial priority queue; exact is O(n² log n) (`READ`). **Marching-cubes output is maximally
obtuse**, so the classical wavefront methods are the worst possible fit for the one substrate that
could afford them — while the heat method is demonstrated on *"an extremely poor triangulation with
significant noise"*.

**The math.** Crane, Weischedel & Wardetzky, *Geodesics in Heat*, DOI `10.1145/2516971.2516977`
(**downloaded, converted, indexed, `READ` in full**): solve `(id − tΔ)u = δ`, normalise `X = −∇u/|∇u|`,
solve `Δφ = ∇·X`. **The one free parameter is `t = h²`** — and on an authored mesh `h` varies by orders
of magnitude, while on a uniform lattice `h` *is the voxel size*. **The method has zero tuning
parameters on this substrate.** Extended to transport and the logarithmic map by the vector heat method,
DOI `10.1145/3243651` (**downloaded, converted, indexed, `READ`**), which also gives geodesic polar
coordinates over the whole surface — a decal system with no UV atlas.

**The obstacle, named.** The factorisation is the carve cost (Law 2). Three escapes, in order of
preference: **chunk-scoped factorisations** (only re-chunked regions refactor); **CHOLMOD
update/downdate** at cost proportional to the changed entries of `L` — and **M-318 is exactly the
statement that those are O(brush support) rather than O(N)**, which converts this from disqualified to
incremental; or drop the prefactored family entirely for the **Closest Point Method**, where Dziuk &
Elliott (DOI `10.4171/ifb/182`, **downloaded, converted, indexed, `READ`**) show the matrices *"depend
only on the evaluation of the gradient of the level set function Φ"* — a brush touching k voxels
changes exactly k stencils. A grid-free fallback exists too: projected walk-on-spheres, DOI
`10.1145/3680528.3687599` (**downloaded, `READ`**) — no matrix, nothing to invalidate, but 17 s for
28k evaluation points on 40 cores, so not yet real-time.

**Pre-registered prediction.** On a 64³ chunk (~15k Surface Nets vertices), a radius-4 brush changes
≤ 400 vertex slots (extrapolating M-318's 346/15,706), and a CHOLMOD update restricted to those rows
re-establishes a valid factorisation in **< 5 ms** against **> 100 ms** for a full refactorisation — a
≥ 20× gap. **If the update is not ≥ 10× cheaper, the prefactored family is dead for live carving and
you go to CPM.**

## 2.3 SKIN THAT GROWS ON A FRESH CUT — reaction–diffusion as surface state, not texture

**The mechanic.** Every surface carries two chemical concentrations running a Gray–Scott system *on the
surface itself*, so exposed rock grows spots, stripes or lichen over play-minutes with **no UVs and no
seams** — the pattern wraps a corner and continues into the tunnel you just dug because there is no
texture to wrap. Cutting does not erase it: it creates a fresh boundary of unpatterned substrate, and
the old pattern grows into the cut and **collides with itself at a visible merge seam**, so a wall
carved twice a minute apart is legibly two different ages.

**Why nobody has it.** Turk 1991 (DOI `10.1145/122718.122749`, in corpus, `READ` via chunks) could not
run RD on the render mesh at all — he had to synthesise a separate computation mesh, requiring *"only
that the model be divided up into relatively evenly-spaced regions."* Authored meshes are neither evenly
spaced nor seam-free, so shipped games bake RD offline into a 2D texture and lose all coupling to
geometry. And the concentrations live on vertices, which a re-mesh renumbers — annihilating the
chemistry.

**What the substrate buys, and it is the cleanest fit in the document.** Surface Nets on a uniform
lattice *is* Turk's computation mesh — one vertex per crossed cell, spacing bounded by the voxel size,
adjacency from cell adjacency — generated as a side effect of rendering. And **the grid-edge vertex key
(M-318) supplies the second piece Turk never needed and authored pipelines cannot give: a concentration
attached to a grid edge survives arbitrary re-meshing, because the edge is a property of the lattice,
not of the triangulation.** Finally, the ordered brush log gives per-vertex time-since-exposure for
free — which becomes the spatially varying diffusion coefficient Turk used to widen stripes, so
**carving history determines pattern morphology at zero extra storage**.

**Cost — this is the mechanic immune to Law 2.** RD is an explicit local stencil, O(V) per step, no
solve, no factorisation. A local edit changes the adjacency of the changed slots and nothing else; the
timestep after a carve costs the same as the timestep before. Storage is 2 floats per *surface* vertex,
O(n²), against O(n³) for the 3D-texture approach games use now. The one real risk is stability: the
cotan Laplacian on marching-cubes slivers has **negative weights** and explicit diffusion will diverge —
the fix is the tufted intrinsic-Delaunay Laplacian (Sharp & Crane, DOI `10.1111/cgf.14069`, `ABSTRACT`,
validated on all of Thingi10k), and **whether it is mandatory or merely prudent is the go/no-go metric
of the harness.**

**Pre-registered predictions.** After 1,000 random brush edits, the fraction of surface vertices whose
concentrations were spuriously reset is **< 3%** with grid-edge keys against **> 90%** with append-order
keys — the same regime as 15,706 → 346. One explicit RD step over ~15k vertices costs **< 0.15 ms**
single-threaded.

## 2.4 THE THROAT — chokepoint width as a live number in metres

On the λ-skeleton, `ρ` has local minima at saddles: those are the throats, and `2ρ` is the passage
diameter **in metres**, so "how many abreast" is `floor(2ρ / shoulder_width)`. Blow a hole in the flank
and a new skeleton branch appears, the old throat's betweenness collapses, and **the holding position
becomes automatically invalid** — the AI relocates because the derived graph changed, not because a
designer volume was invalidated. No shipped game derives chokepoints; every tactical AI I know of uses
designer-placed cover volumes or influence maps over a static navmesh, all invalidated the instant the
player digs. **Prediction:** derived throat diameter matches a tape measure within ±1 voxel on ≥95% of
throats across 20 scenes; the tactical graph updates within ≤2 frames at ≤0.2 ms amortised. **Risk:**
spurious throats from destruction noise — filtering by `‖∇ρ‖ < 0.5` (θ > 120°) should remove ≥90% of
throats a human would not call a doorway, and if it does not, throat detection needs a global collapse
metric and falls to the background budget.

## 2.5 SPELEOTHEM MORPHOGENESIS — stalactites grow where you made it drip

Water seeping through ceilings **the player created** grows stalactites with the *correct shape* — a
specific tapering profile real stalactites converge to — with stalagmites rising to meet them. Divert
the water and growth stops, leaving a dead stub. Minecraft ships pointed dripstone, so this is the one
candidate with a real shipped precedent: it grows one *block* per random tick on a lattice, with no
shape law and no dependence on drip rate. **The novelty is morphogenesis** — the shape is derived, not
authored.

The physics is Short, Baygents, Beck, Stone, Toomey & Goldstein, PRL 94, 018501, DOI
`10.1103/physrevlett.94.018501` (`ABSTRACT`; **unobtainable — no OA PDF and no arXiv preprint**): thin-film
fluid dynamics plus calcite chemistry plus CO₂ outgassing give *"a local geometric growth law… the local
growth rate is proportional to the local thickness of the fluid layer"*, with *"extreme amplification at
the tip"*. And the result that makes it nearly free: *"a broad class of initial conditions is attracted
to an ideal universal shape, whose mathematical form is found analytically"*, validated against
stalactites measured at Kartchner Caverns.

**Because there is a universal attractor with a closed analytic form, you do not simulate a
free-boundary problem at runtime — you evaluate one profile function with one parameter.** State is a
position, an axis, a length and a drip rate; the cost is a scalar increment and an occasional brush
rewrite. **The paywall blocks the closed form, and the workaround is the harness itself**: a ~100-line
1D free-boundary integrator started from a cone, a cylinder, a bulge and a kink, confirming they *all
converge to the same shape* — which both tests the paper's central claim and recovers the profile
numerically.

---

# Part 3 — Tier 3: the honest losers, and why

Recording these matters as much as the winners, because each has an attractive surface and a specific
number underneath it.

| Candidate | The number that kills it |
|---|---|
| **Fracture modes on carved geometry** | Breaking Good (`10.1145/3549540`, in corpus, `READ`): *"Each mode takes between 0.5 and 12 seconds to compute"* on 3k–15k tet meshes, and *"Computing a new set of fracture modes for each piece would exceed realtime constraints."* Ten modes = 5–120 s per shape. The bottleneck is a **conic program**, not an eigensolve. Their own optimisation claims "several orders of magnitude faster", which still lands at ~1–100 ms *per mode*. **A background job, not a frame budget — do not promise it at 60 Hz.** |
| **Modal sound from carved shapes** | Picard et al. (`10.1155/2010/392782`, downloaded, OCR failed, table read via publisher fetch — *single-sourced, verify before betting*): 81 modes in 0.05 s, 1191 modes in 21.9 s, scaling ≈ O(m^2.8). Extrapolating to a 5 ms budget gives **≈40 modes**. Full-bandwidth modal audio for arbitrary carved shapes is dead; a 40-mode bell is credible. The good news is the discretisation: **hexahedral FEM on a voxel grid**, no tetrahedralisation, matching the substrate natively — where Breaking Good needs TetGen plus an 85-second caging step. |
| **Handholds from solid thickness** | The λ-medial axis *"removes salient shape features before noise"* (`READ`) — and a thin rib on a thick slab is exactly that case. Predicted rib retention **< 50%**. If confirmed, this needs the Scale Axis Transform, which requires **two full medial-axis computations per evaluation** and is not a per-edit operation under any budget. |
| **Ligament / severance as collapse criterion** | Cutting *splits* solid components, and union-find cannot split. Costs a decremental structure or a localised rebuild; predicted **> 1 ms for structures above ~10⁵ solid voxels**, which forces it permanently off the hot path. |
| **3D repose-angle relaxation (talus, slumping)** | Class B: dense state, forces a grid and a bake, with all four costs in Law 1. The mathematics is the best-grounded in the sweep — Musgrave 1989 (`10.1145/74334.74337`, in corpus, **`READ` in full**) gives the operator and its asymptotic talus-angle property; Peng et al. (`10.1006/jcph.1999.6345`, `READ`) give the O(N) narrow-band cost law and the CFL. **Ship it as a bounded, player-triggered, locally-baked region** — you get roof collapse and talus cones without globally sacrificing the representation. |
| **Buckling as a rule** | Dahlberg (`READ`): *"Eigenvalue problems are notoriously expensive… This cost becomes prohibitive for medium to large problems."* **Never solve a generalized eigenproblem at runtime.** Use a slenderness proxy: effective length over minimum radius of gyration from SDF moments, which are cheap voxel sums, with an Euler-style threshold. O(cells in path), no solve. No shipped game distinguishes buckling from crushing. |
| **Topology whorls (Poincaré–Hopf as a law)** | Beautiful and mostly unreachable: globally optimal direction fields need a **sparse eigenproblem per stroke**. But the cheap half is real and free — χ = V − E + F is computed during meshing, Σ(singularity indices) must equal 2 − 2g, and **drilling one through-hole changes Σindex from +2 to 0**. Use χ to *budget* how many whorls a chunk's field may have, and let the field itself be local. It is also the cheapest correctness probe of the whole intrinsic-geometry substrate: if the discrete index sum does not track χ, your Laplacian is wrong and §2.2 and §2.3 are both suspect. |
| **Freeze–thaw** | Maximum novelty, no math to stand on. Searches returned laboratory rock mechanics and Alpine rockfall monitoring, not a simulation model. You would be inventing it. One usable fact survives: rockfall peaks during snowmelt with a measured **2 h lag between peak air temperature and peak rockfall** (`10.5194/nhess-22-445-2022`). |

---

# Part 4 — Corrections

## 4.1 The medial-axis foothold is half wrong, and the wrong half opens the line

The hunt was given M-172 as a foothold: `BrushStack::gradient` returns exactly `[0.0, 0.0, 0.0]` on a
slab's mid-plane, therefore the medial axis is free.

**The exact zero is not the medial axis. It is one measure-zero point of it.** From the identity in
§1.1, on the medial axis `‖∇ρ‖ = cos(θ/2)` where θ is the separation angle. Exact zero requires
θ = π exactly — the two nearest boundary points antipodal through `x`. That is the mid-plane of a
slab, which is precisely what M-172 measured. **In a wedge, a corner, a tunnel junction, or anything a
player actually digs, `‖∇ρ‖` is strictly between 0 and 1 and exact-zero detection finds nothing.**
Worse, a sub-voxel translation of the slab destroys even the slab case.

**But the finding survives in a stronger form, which is why this is a re-tier rather than a
falsification.** `‖∇ρ‖` is a *continuous* readout of exactly the quantity the instability literature
says you need — small `‖∇ρ‖` and large `ρ` is medial and stable; `‖∇ρ‖ → 1` is noise. And Attali &
Montanvert's experimentally-derived noise criterion (`READ` in the survey) — keep the "upper right
quadrant" of the (θ, ρ) parameter graph — *is* "small `‖∇ρ‖`, large `ρ`", with boundary noise appearing
as *"a hyperbola-like point cloud located near the coordinate axes."* The literature validated the
filter direction; the identity makes it one multiply and one square root on data already resident.

**The cheapest experiment in this entire document tests it, and it takes thirty minutes.** Offset a
slab's mid-plane by half a voxel. Predicted: the count of samples returning exactly `[0,0,0]` drops to
**zero**, while the count with `‖∇ρ‖ < 0.1` changes by **< 5%**. If confirmed, retire "exact zero
detects the medial axis" and adopt "`‖∇ρ‖` magnitude scores medial stability" — and the M-172 row gets
its correction the way ✗-rows always do here, with the old reading left visible.

**A second prediction worth registering before anything is built**, because it decides whether the
theoretical guarantee is available at all: the homotopy certificate needs `λ < wfs`, where the weak
feature size is the smallest distance from a critical point of `∇ρ` to the boundary. **Destruction
creates near-tangential features, so critical points sit close to the boundary — predicted `wfs < 2
voxels in > 80% of dug scenes**, meaning `λ < wfs` essentially never holds. If confirmed, drop the
homotopy guarantee and rest on Hausdorff stability alone (Chazal–Lieutier; Lieutier & Wintraecken's
`d_H ≲ d_H^{1/2}`, `d_GH ≲ d_H^{1/4}`, DOI `10.1145/3564246.3585113`, `ABSTRACT`) — weaker, and still
far better than a tuned persistence constant.

## 4.2 A definition discrepancy that must be resolved before implementing

The Attali survey defines the λ-medial axis as `M_λ = {x : r(x) ≥ λ}` — **keep** large `r`. The 3D
Skeletons STAR paraphrases it as discarding medial samples whose circumradius is **larger** than λ. One
of these is inverted, or they describe dual constructions. **The primary source that settles it,
Chazal & Lieutier's `10.1016/j.gmod.2005.01.002`, is currently unobtainable** and is the highest-value
acquisition in this document.

Related, and worth the same caution: the STAR attributes Foskey–Lin–Manocha to "distance-to-boundary"
measures, while FLM's own abstract says the θ-SMA is parameterised by a **separation angle**. Trust the
primary abstract.

## 4.3 One mislabelled document in the corpus

The stem `sig2024_A_Heat_Method_for_Generalized_Signed_Distance` is **not** Feng & Crane's signed heat
method. Page 1 reads *"An ADMM-based scheme for distance function approximation, Alexander Belyaev,
Pierre-Alain Fayolle"* (`READ`). Anyone citing that stem for the signed heat method cites the wrong
paper; the real one is DOI `10.1145/3658220`.

---

# Part 5 — What to run first

Ranked by what each decides per hour spent. The first four are all under a day and three of them can
kill a whole branch.

| | Experiment | Cost | What it decides |
|---|---|---|---|
| 1 | **Sub-voxel slab offset** (§4.1) | 30 min | Whether M-172's reading has to be reframed from a boolean to a magnitude. I expect it does |
| 2 | **The `r = ρ√(1−‖∇ρ‖²)` oracle check** at three voxel sizes (§1.1) | half a day | Whether the identity survives discrete gradients with O(h) error. **If it does not, §1.1, §2.4 and §3.3 all die together** — check it before any of them |
| 3 | **Speleogenesis graph simulator**, headless, no field, no renderer (§1.2) | a day | Whether the wormhole competition and the bimodal aperture distribution reproduce. Validates the kinetics against published geomorphology before a voxel is touched |
| 4 | **`wfs` histogram on dug scenes** (§4.1) | an afternoon | Whether the homotopy certificate is available at all. I expect not, and the direction survives on Hausdorff stability |
| 5 | **Mode-shape drift under a one-voxel edit** (§3, modal) | a day | Whether a single dig perturbs λ₁ above the just-noticeable difference at all. **The cheapest possible kill-shot on the modal direction — run it before assembling anything** |
| 6 | **Arch golden values: 0.1075 and 15.84°** (§2.1) | a day | Whether a feasibility solver is implemented correctly, against the only external ground truth in this whole document |
| 7 | **CHOLMOD update vs refactor on a 64³ chunk** (§2.2) | a day | Whether the prefactored family is reachable for live carving, or whether everything intrinsic goes to the Closest Point Method |
| 8 | **Sabine RT60 per air component** (§1.3) | a day | Nothing — it just works, and it is the cheapest thing in the document that a player can hear |

**One meta-note on ordering.** Experiments 1, 2, 4 and 5 are all *falsifiers of premises*, not
measurements of performance. Three of them I expect to come back negative, and that is why they are at
the top: a negative here costs a day, and the same negative discovered during implementation costs a
sprint. This is the same discipline that caught ✗26 before the harness ran.

---

# Part 6 — Acquisitions, and four corpus problems worth fixing

## Downloaded and indexed this session

| Identifier | Work | State |
|---|---|---|
| `10.1145/2516971.2516977` | Geodesics in Heat | converted, indexed (15 chunks), **`READ` in full** |
| `10.1145/3243651` | The Vector Heat Method | converted, indexed (49 chunks), `READ` |
| `10.1145/3680528.3687599` | Projected Walk on Spheres | converted, indexed (20 chunks), **`READ` in full** |
| `10.4171/ifb/182` | Dziuk & Elliott, Eulerian FEM on implicit surfaces | converted, indexed (17 chunks), `READ` |
| `10.1145/1618452.1618458` | Whiting et al., structurally-sound masonry | converted, indexed (15 chunks), **`READ` in full** |
| `10.1007/s00158-023-03616-7` | Dahlberg et al., buckling TO with reanalysis ROM | converted, indexed (15 chunks), `READ` |
| `10.3389/fbuil.2020.00043` | Pepe et al., limit analysis vs FEM/DEM | converted, indexed (20 chunks), `READ` |
| `10.1145/1964921.1964933` | Toward high-quality modal contact sound | converted, indexed (20 chunks), `READ` |

Also downloaded, conversion failed or pending: `10.5194/hess-23-1995-2019`, `10.5194/hess-18-4617-2014`,
`10.3986/ac.v44i3.1925`, `10.21203/rs.3.rs-4256563/v1`, `10.1016/j.patrec.2010.09.002`,
`10.1590/1679-78251252`, `10.1109/tac.2020.3007688`, `10.1007/s11831-019-09351-x`,
`10.1155/2010/392782`, `10.1145/1833349.1778806`, `10.1145/3641519.3657493`, `10.1145/3592404`,
`10.1371/journal.pone.0190666`.

## The acquisitions that would most change the analysis

1. **`10.1016/j.gmod.2005.01.002`** — Chazal & Lieutier, *The "λ-medial axis"*. Resolves §4.2's
   definition discrepancy, which blocks implementing §1.1 correctly.
2. **`10.1145/2856317`** — Yeung, Crouch & Pothen, *Interactively Cutting and Constraining Vertices in
   Meshes Using Augmented Matrices*. Real-time FEM under **cutting**, augmented-matrix Schur complement
   on the original Cholesky factor, demonstrated to 167,000 elements. **The single most on-point paper
   for the incremental question anywhere in this sweep.**
3. **`10.1145/3673652`** — King et al., Closest Point Method with interior boundary conditions. Decides
   whether the CPM escape route from Law 2 actually hits interactive rates.
4. **`10.1145/3414685.3417828`** and **`10.1145/3130800.3130849`** — Herholz et al., sparse Cholesky
   updates for interactive mesh work. The numbers behind §2.2's prediction.
5. **`10.1111/cgf.14099`** — Planeverb. §1.3's reference implementation.
6. **`10.1145/2366145.2366178`** — Whiting et al. 2012, closed-form structural gradients. §1.4's source.
7. **`10.1111/cgf.14069`** — Sharp & Crane, nonmanifold Laplacian. Decides whether the tufted operator
   is mandatory on marching-cubes output or merely prudent, which is §2.3's go/no-go.
8. **`10.1103/physrevlett.94.018501`** — the stalactite profile. No OA PDF and no arXiv preprint.
9. **`10.1145/2835173`** — Shin et al., *Reconciling Elastic and Equilibrium Methods*. Settles whether
   FEM is *wrong* here or merely expensive.

## Four infrastructure problems, each reproducible

- **`scribe`'s olmocr backend failed on the majority of PDFs this session** — three distinct modes
  (missing workspace output, 0/N pages completed, 60 s MCP timeout). Retrying sometimes routed to
  `glm_ocr` and succeeded, so persistence is the workaround. `pipeline_drift` reported **79–84** against
  a threshold of 3.
- **Wiley, ACM and HAL DOIs silently resolve to landing pages that convert "successfully" at ~100 KB.**
  Three papers in this sweep were nearly cited from a landing page. **A size/content heuristic — a
  sub-150 KB "paper" is always a landing page — would have caught all three.** `hal.science` is
  additionally blocked by the agent proxy at CONNECT, so every HAL-hosted LIRIS paper fails this way.
- **`paper_search` with `provider: arxiv` returned `[]` for every query tried**, including for papers
  whose arXiv PDFs the `all` provider had just surfaced. Treat the arXiv provider as non-functional and
  use `all`.
- **Presence and searchability have come apart in both directions.** Turk 1991 returns chunks from
  `distill_search` but `catalog_read` says "No catalog entry" and `markdown_read` says "Markdown not
  found" — the *inverse* of ✗4's failure mode. So the rule "presence is decided by `catalog_list`" needs
  its converse added: **check both, because either can be the one that is lying.**

---

# Part 7 — The one-paragraph version

**The strongest single finding is a derivation, not a paper**: `r(x) = ρ(x)·√(1 − ‖∇ρ(x)‖²)` makes
λ-medial-axis membership an O(1) test on the SDF value and its gradient — both already computed at
every sample — and the threshold λ is *a length in metres equal to a creature's shoulder width*, which
is the principled threshold the closed Morse–Smale line lacked. The academic complaints about that
skeleton (it disconnects, it deletes small features) are the game's correct answers. **The second
strongest is that a graphics group has already published this project's exact representation for cave
networks** — a signed-distance construction tree over a karst skeleton — and generates it statically,
so driving it dynamically from flow is a one-sentence gap that yields a world which digs its own caves
along the route the player's excavation chose, entirely inside the edit log, inheriting M-36's
commutativity for free. **The third is that the structural question every destructible game gets wrong
is not a stress question but a feasibility question** — an LP with two physical parameters and no
stiffness matrix, calibrated against two closed-form masonry results, which Whiting solved in 2009 for
authoring and nobody has moved to after the player's edit. Underneath all of it sit two laws worth more
than any candidate: **a process belongs in the edit log exactly when its state is sparse**, and
**classify by whether a local edit has a local answer, not by whether the operator is cheap** — with
the one door through that wall being a reanalysis convergence test, `‖K̄⁻¹ΔK‖₂ < 1`, whose *divergence
is the collapse signal*.
