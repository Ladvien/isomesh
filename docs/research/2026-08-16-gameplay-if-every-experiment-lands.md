# What the game looks like if every experiment comes back positive

**Date:** 2026-08-16
**Companion to:** `2026-08-16-sota-speed-and-feature-frontier.md`
**What this is:** each candidate and harness from that document, translated into what a player would
actually see and do. Grounded in the measured numbers where they exist, and explicit about the four
that would produce **nothing a player can perceive** — which is as useful to know as the rest.

**How to read the estimates.** Where a speedup converts into view distance I use the cube-root rule: a
`k×` throughput win at fixed frame budget buys `k^(1/3)` more radius, because the resident set is a
volume. A 1.8× win is 22% more distance, not 80%. That constant is why several of these read smaller
than the headline suggests.

---

# Part 0 — The composite: one game, everything positive

Before the itemised list, the thing they add up to. Not a wishlist — every sentence below traces to
a specific experiment, named in brackets.

You are standing in a cave you dug. The walls carry **tool marks** — not a texture, actual geometry,
kerfs thinner than a voxel left by a chisel [1.9, 2.6]. Where you cut into a seam of ore, the rock has
**strata**: the material changes at a real boundary because the world is a multiphase field, not one
scalar with paint on it [1.2/PhaseTree]. The cut edges are **sharp** — a chisel leaves a corner, not a
rounded blob — and they are sharp *everywhere*, on the machined surfaces you carved and on the natural
rock around them, because dual contouring finally ships without self-intersecting on rough terrain
[2.1, 2.8a].

You have been playing on this server for four months. The region around the settlement has tens of
thousands of accumulated edits, and it re-meshes at the same cost as pristine wilderness [1.2]. The
game does not punish you for having played.

You seal the cave mouth. The game **knows** it is sealed — not because a designer placed a trigger
volume, but because the mesh's topology provably matches the field's, so "this volume is closed" is a
query rather than a guess [2.3, 2.4, 3.1]. Water pools instead of leaking. The gas you released stays
in. When you later breach the far wall, the engine fires a **breakthrough event** the instant two
regions merge, and the water finds the new opening on the same frame [2.3 + round-2's dynamic
connectivity].

You set charges. The building comes down as a *building* — glass spiders, timber splinters along its
grain, stone cleaves on its bedding planes — because fracture follows precomputed modes with an
anisotropy field rather than Voronoi shards [2.8c + 2.8j]. Every fragment has correct collision
[1.3, 2.1]. Nothing hitches: the frame that the charge lands on is the same length as every other
frame, because the re-mesh queue is budgeted and its p99.9 is a number somebody measured [3.3].

You fly out over the world at speed. Terrain does not pop in ahead of you and does not pop *at all* as
LOD changes — the switch distances were chosen so that the worst vertex displacement subtends less than
a pixel [3.8, 2.8l]. Look down: the cave system you are flying over is a real cave system, with the
topology the generator intended, because no cell was resolved by fiat [2.4]. Look at the horizon: it is
22–40% further away than it would otherwise be [1.1, 1.3, 1.5+parallelism].

You point a cutting beam at a wall and hold the trigger. The wall opens **continuously**, at the rate
your hand moves, not in discrete chunk-sized pulses — because a 3 cm cut costs 3 cm of work [3.7].
Vehicle tracks leave ruts. Water carves a channel over minutes of play.

None of that is a new rendering technique. It is six speed results and four correctness results.

---

# Part 1 — Speed experiments

## 1.1 Conservative-SDF empty-space skipping

**Result if positive:** >60% fewer field evaluations on conservative fields with bit-identical output;
`fbm_terrain` gets ≥1.8× from the brick pass.

**What the player sees.** **Distance, and the absence of holes.** M-124 measured 15.15 chunks/frame at a
2 ms budget; 1.8× makes it ~27. In view-distance terms that is +22% radius — noticeable but not
transformative on its own. The transformative part is elsewhere: **M-104 measured 606 chunks
permanently waiting** in the naive residency case, and the visible symptom was *holes in the terrain
that never fill*. A 1.8× win moves the boundary of that regime outward, which means **the game can
support a faster player**.

**The mechanic it unlocks.** Vehicles and flight. Not "a flying camera" — a flight speed at which the
streaming queue still drains. Right now the design constraint is that a player who moves faster than
the mesher outruns the world, and the failure is ugly and immediate. Positive here converts "walking
game with a flying debug camera" into "you can own an aircraft."

**Honest ceiling.** Nothing on `gyroid`-class dense fields — the prediction is explicitly ≥5% cost and
<10% gain there. A world made of dense triply-periodic structure gets nothing.

---

## 1.2 In-tree SDF acceleration / proxy nodes

**Result if positive:** field evaluation is >70% of extraction (confirming M-136), and collapsing a
60-brush edit log into a conservative proxy cuts per-chunk cost ≥2× for long logs.

**What the player sees.** **Nothing — and that is exactly the point.** This is the experiment whose
positive result is defined by an *absence*.

M-50 measured the current behaviour: median cost per re-meshed chunk of **0.158 / 0.354 / 0.525 /
0.589 ms** for edit logs of 1–15 / 16–30 / 31–45 / 46–60. The world gets slower the more you play in
it. On a persistent server that is a slow, unavoidable decay: the town square, the mine everyone uses,
the base you have been extending for a month — all of them the most expensive places in the world,
precisely because they are the places people care about.

**The mechanic it unlocks.** **Persistence.** Months-long servers. Shared worlds where the accumulated
history of every player is permanent and free. "Your changes are saved forever and cost nothing" is a
promise the current numbers cannot keep, and it is the promise that separates a sandbox from a session.

**The second half — PhaseTree, at <25% overhead.** A multiphase field *is* a multi-material world.
Positive here means stone/dirt/ore/ice are **phases of the field**, not colours painted on a surface.
Player-visible: you dig through a hillside and hit an actual layer boundary — topsoil, clay, bedrock —
with different hardness, different tool requirements, different sound. Ore veins have shape and follow
the geology instead of being a texture swap. This is the difference between mining and recolouring.

**Caveat that belongs to you, not to the research.** Multi-material is a public-API change. The doc
flags it as a scope question rather than a ticket.

---

## 1.3 Collider readiness (the unexamined 45%)

**Result if positive:** most chunks skip the readiness census because their extractor carries the
guarantee; the 45% collapses for the default path.

**What the player sees.** **How much of the world is physical.** Nearly 2× the chunks per frame is
+26% radius, but the sharper consequence is qualitative: right now the cost of making a chunk
*collidable* is larger than the cost of making it *visible*, so there is structural pressure toward a
world where distant geometry is scenery. Positive here removes that split.

**The mechanic it unlocks.** **Destruction at a distance.** Artillery. A demolition you watch from
across the valley where the rubble actually lands on the rubble. Snipers shooting through walls a
kilometre out. NPCs that path across terrain you carved this frame rather than next frame.

**And the adjacent result changes the category.** SDF↔SDF collision (`10.1016/j.cagd.2024.102305`)
skips collider generation entirely for field-defined geometry. If that lands, the number is not 45%
smaller — it is **zero**, and collision becomes valid *mid-carve*. Player-visible: you can stand on a
surface while it is being cut out from under you, continuously, with no frame where the collider is
stale. Drills you ride. Floors that fail progressively. Digging while falling.

---

## 1.4 Decoupled Fallback prefix scan

**Result if positive:** the scan is 1.4× faster than reduce-then-scan.

**What the player sees.** **Nothing.** M-155 already took the whole 129³ GPU path to 0.54 ms; 1.4× on
a fraction of that is microseconds.

**What the *negative* case looks like, which is why this is in the doc at all.** If somebody ports a
CUDA decoupled-lookback scan, the M5 hangs. The player-visible outcome is a **two-second freeze and a
driver reset**, on Apple hardware only, passing CI on the Ryzen. The gameplay consequence of this
experiment is "the game runs on a Mac," which is a platform, not a feature.

---

## 1.5 Stage-streamed extraction (Flying Edges lineage)

**Result if positive:** ≥1.6× on the extraction stage, classify pass under 4 cycles/grid-point.

**What the player sees, taken alone.** Very little. §0.1 of the research doc caps the whole extractor at
**1.41× on a usable mesh**, which is +12% view distance. Marginal.

**What the player sees, taken as what it actually is.** The restructuring's side effect is that the
edge-indexed intersection array **is** the vertex dedup — no hash map. That removes the
iteration-order determinism hazard, which is the thing standing between this crate and **safe
multi-core meshing**. On a 12-core desktop that is not 1.41×, it is closer to 8×, which is +100% view
distance and roughly 8× the chunks per frame.

**The mechanic it unlocks.** Scale. A world where the streaming radius is bounded by memory rather than
by a single core. And because determinism survives (M-31, M-36), **destruction stays identical across
clients while only the edit log is synchronised** — no mesh replication, no desync, coordination-free
multiplayer editing inside a run of same-kind brushes.

**Honest framing.** The 1.6× is the advertised result and it is the boring half. The parallelism
unlock is the reason to do it, and it is a consequence rather than a measurement — which means it
needs its own experiment.

---

## 1.6 Workgroup-shared GPU reduction staging

**Result if positive:** 14.5× on a reduction pass.

**What the player sees.** **Nothing.** GPU-010a already banked this shape of win (5.24 → 0.37 ms) and
M-155 took the path to 0.54 ms total. Kept as a checklist for the next GPU reduction, not as a ticket.

---

## 1.7 Output-driven vs scatter compaction (the crossover)

**Result if positive:** ≥1.5× on chunks below 3% active cells; scan still wins above 25%.

**What the player sees.** **Vertical worlds.** M-104's 606 permanently-waiting chunks were "empty air
or solid rock" — exactly the sparse regime this wins in. Positive here makes the sky and the deep
underground *cheap*, which is the difference between a world that is a shell around a heightfield and
a world with real volume above and below it.

**The mechanic it unlocks.** Deep cave systems that go down for kilometres, and altitude that is not a
skybox. Both currently pay full price for chunks that produce almost no triangles.

---

## 1.8 Morton across chunks, row-major within

**Result if positive:** ≥10% on a 16-chunk re-mesh, and **±5% (i.e. nothing) at 32³** — the null result
is the finding.

**What the player sees.** Marginally larger explosions inside the same frame budget. ~3% more view
distance. This is a tidy-up, and the doc says so.

---

## 1.9 Subgrid MT: fewer triangles at equal error, plus the root cache

**Result if positive:** ≥20% fewer triangles at matched Hausdorff, and the grid-edge root cache cuts
field evaluations ≥5× with zero golden-hash changes.

**What the player sees.** **This is the biggest feature result in the entire document, and it is
disguised as a speed result.**

M-98 measured subgrid at **70× classic MT and 196× Marching Cubes** — 22.5 ms per chunk at 33³, which
is a demo, not a shipping path. A 5× cut takes it to ~4.5 ms, which is a budget you can actually spend
on the chunks that need it. And what it buys is **M-95**: subgrid resolves a slab **1/20 of a cell
thick**, on a grid where a sign-based method returns *nothing at all*.

**The mechanic it unlocks: sub-voxel geometry.**

- **Engraving.** Carve your name into a wall with a knife. The letters are geometry — collidable,
  shadow-casting, visible in silhouette — at a scale far below the voxel grid. Currently impossible:
  A-005 measured greedy quads returning **zero triangles** on `thin_plate`, and M-72 measured Marching
  Cubes turning the same feature into a resolution-dependent scatter.
- **Thin structures as first-class objects.** Wires, chains, grates, railings, sheet metal, glass
  panes, ice. Things that are currently faked with props because the field cannot represent them.
- **Kerf.** A sword slash leaves a cut narrower than a voxel, and the cut is real. So is the crack
  that propagates from it.
- **Detail that survives distance on a curve you choose.** M-72's framing, from the round-2 doc: a
  feature that vanishes at a *known* distance can be faded; one that disintegrates into scatter pops.
  Positive here means you author the distance at which the engraving stops existing.

**Why the cache is now reachable.** The round-2 doc deferred it because the correctness precondition
was unstated. A-014g stated it — M-168 gave every crossing a global identity keyed on `(the tet edge's
two grid points, the root's ordinal)` — and M-184 completed it. **That key is the cache key.** The
acceptance test is zero golden-hash changes.

---

# Part 2 — Feature and correctness experiments

## 2.1 Quad splitting → intersection-free dual contouring

**Result if positive:** self-intersections go to **exactly 0** on `gyroid` and `fbm_terrain`, triangle
count up <15%, Hausdorff unchanged.

**What the player sees.** **Sharp edges that ship.** M-54 measured Dual Contouring at **101× more
accurate than Marching Cubes on `box_exact`** and 77.9× on `thin_plate` — corners reproduced to 0.01
cells against 0.58. That is already true and already unusable as a default, because M-53 measured DC
self-intersecting at 13.837 per 1k triangles on `fbm_terrain`, and the SIGGRAPH-deck testimony is that
the shipping symptom is **broken lighting** rather than broken geometry. Glitched shading on natural
terrain is not something you can ask a player to ignore.

Positive here means one extractor gives you both: crisp machined geometry *and* clean rough terrain.

**The mechanic it unlocks.**

- **Tools that leave the mark of the tool.** A chisel leaves a chisel edge. A plasma cutter leaves a
  straight kerf with corners. A blast leaves a rough face. The *shape* of the cut communicates what
  made it — which is a legibility mechanic, not just a visual one.
- **Building.** Player-placed structures that read as constructed rather than extruded blobs. Right
  angles that stay right.
- **Reliable convex decomposition** for debris. M-116 measured 241–272 ms per fragment — 14–22 whole
  frames — and the decomposer's failure modes get worse with self-intersecting input. Clean input is a
  precondition for making that offline path trustworthy at all.

---

## 2.2 MMS convergence, and whether λ = 0.01 may be constant

**Result if positive (in the interesting direction):** fixed-λ Tikhonov is **not** second-order, and
annealing λ with `h` restores it.

**What the player sees.** **Detail that improves as you approach.** Right now LOD refinement makes a
smooth surface smoother, but a carved corner stays exactly as rounded as λ makes it — the bias does not
shrink with the grid. Positive here means walking toward a cut edge makes it visibly, measurably
crisper.

**The mechanic it unlocks.** Anything built on inspection. Forensics, archaeology, appraisal, "lean in
and look at the tool marks to identify who made this cut." A world where getting closer is rewarded
with information rather than with the same blob at higher resolution.

**And the negative result is a shipping constraint worth having.** If λ must anneal, then λ is a
function of chunk LOD, and every LOD level has a different vertex rule — which interacts with the
transition cells A-011 already owns. Better to know before it is load-bearing.

---

## 2.3 The stratified-Morse-theory χ oracle

**Result if positive:** the extractor's topology is provably correct on random trilinear fields, not
merely self-consistent.

**What the player sees.** **This is the one that unlocks an entire mechanic family, and none of it is
visual.**

The round-2 doc named the gap precisely: *"Nobody has established that a region sealed in the scalar
field is sealed in the extracted mesh."* Every mechanic below is currently unshippable for that reason
— you cannot build a game rule on a guarantee you do not have, because the failure mode is a player
drowning in a room the engine thought was open.

**The mechanic family it unlocks — sealed as a queryable predicate:**

- **Airtightness.** Build a base, seal it, and the atmosphere stays. Breach it and the atmosphere
  leaves, through the hole you made, at a rate the hole determines. Airlocks that are airlocks because
  of their geometry, not because a designer tagged them.
- **Flooding.** Water that respects the volume it is in. Dig into an aquifer and the water goes where
  the space is. Wall it off and it stops.
- **Containment.** Gas, fire, spores, a creature. "It cannot get out" as a physical fact.
- **Breakthrough as an engine event.** The instant your tunnel connects to an unexplored chamber, the
  game knows two regions merged — one union-find merge in the air sublevel set. That is a *moment*:
  the sound changes, the air moves, something on the other side notices. No engine currently has this
  signal.
- **AI that understands rooms.** "This chamber has exactly one exit" is derivable, so an ambush, a
  siege, a chokepoint defence, or a "you are trapped" reveal all become structural rather than authored
  — and they survive the player rearranging the terrain, which a baked navmesh does not.

**Why it belongs to this experiment.** All of it rests on field topology and mesh topology agreeing.
The χ oracle is what makes that an assertion rather than an assumption.

---

## 2.4 Isotopy / interval arithmetic, and the genus-stability sweep

**Result if positive:** the fraction of cells "resolved by fiat" is known and small; a sub-voxel origin
sweep on `gyroid` produces exactly **one** distinct χ instead of ≥3.

**What the player sees.** **Generated worlds that are the worlds that were generated.**

Right now a procedural cave system's *topology* is partly an accident of where the grid landed. Two
chambers the generator meant to connect can come out sealed; two it meant to separate can come out
joined. At 64³ the sweep predicts ≥3 distinct Euler characteristics from nothing but sub-voxel origin
offsets — the same cave, meshed at three different genera.

**The mechanic it unlocks.** **Quest and level design on procedural terrain.** "There is exactly one
way into the vault." "This cave is a loop, not a dead end." "The river reaches the sea." Those are
statements a designer can only make if the mesher preserves what the generator built. Positive here
turns procedural generation from decoration into structure — and it composes directly with 2.3, since
a sealed-volume mechanic is worthless if sealing is a coin flip at the grid scale.

---

## 2.5 Off-surface deviation (`mean |f(v)|`)

**Result if positive:** the published 39.09-vs-0.00 gap does **not** reproduce; DC and MC differ by <2×
on smooth fields and >5× on sharp ones.

**What the player sees.** Almost nothing directly. Its value is as a **shading canary**: a vertex that
is off the isosurface has a normal that disagrees with the field's gradient, and disagreeing normals
are shimmer. Two lines of code, and it is the only direct measurement of the Tikhonov bias.

---

## 2.6 `thin_slab(t)` — the thickness-parameterised field

**Result if positive:** three distinct signatures — greedy quads collapse discontinuously at
`t/h ≈ 1`; Marching Cubes degrades into scatter with no clean cutoff; subgrid holds to `t/h ≈ 0.05`.

**What the player sees.** **Thin things that fade instead of shattering.** The measurement itself is
invisible; what it produces is the *number* you need to author a fade. A fence, a railing, a wire, a
grate currently disintegrates into a resolution-dependent scatter as it recedes (M-72's mechanism).
Knowing the exact `t/h` at which each extractor gives up lets you dissolve the feature deliberately at
a chosen distance.

**The mechanic it unlocks.** Nothing new — it makes 1.9's sub-voxel detail *shippable* by giving it an
LOD story. Without this, engraving looks spectacular up close and boils at range.

---

## 2.7 Comparable quality metrics (SA<10°, IN>5°, NV/NE/SI%)

**What the player sees.** **Nothing.** This buys the ability to place this project's numbers beside
FlexiCubes' and TetWeave's tables. Worth doing for the project; invisible in the game.

One thing that could become visible: reconciling the published MC NV = 47–52% against M-53's 0 might
surface a real defect class. If it does, that is a bug fix, not a feature.

---

## 2.8 The smaller items, by what they'd change in play

| Experiment | Positive result | What the player sees |
|---|---|---|
| **a. Gradient-free DC** (arXiv `2604.00157`) | Sharp features from sampled values alone, no Hermite data | **Sharp edges survive editing.** After a carve, the analytic field is gone and only samples remain — so today's sharp-feature path is strongest exactly where the player has not been. Positive means **player-built** structures get crisp edges, and so does voxel data streamed from a server or imported from a scan. This is the UGC unlock |
| **b. DCx over expanded cubes** (DOI `10.1145/3811388`) | Manifold *and* non-manifold surfaces from unsigned fields | **Walls that tear rather than only being carved.** Cloth, flags, leaves, paper, torn sheet metal, membranes. A crack that is a *surface* instead of a thin volume. This is the only published route in the sweep to open/non-orientable output from a grid |
| **c + j. Isosurface stuffing → fracture modes** (DOI `10.1145/1275808.1276448`, `10.1145/3549540`) | A tet mesh with bounded dihedrals; `k` natural break modes, impact response by linear projection at *"no runtime cost"* | **Breakage with a material identity.** Glass spiders. Timber splinters along the grain — literally, via the anisotropy field `η = (10,10,1)` "to favour vertical faults over horizontal ones." Stone cleaves on bedding planes. Instead of every object shattering into the same recognisably-convex Voronoi shards, each material breaks like itself. Precomputation is per-shape and offline, so this applies to **authored props**, not to freshly carved chunks — a real limit worth stating |
| **f. Winding-number solidity probe** | A validity check stronger than manifoldness | **You do not fall through the world.** Catches the "looks solid, physics disagrees" class that manifoldness passes |
| **g. Angle-weighted pseudonormal** | <2° difference on good triangles, >15° on radius-ratio <0.15 | Fewer shading artifacts **exactly where the slivers are**, which is exactly where a sharp carve is. Subtle, cheap (~40 lines), and it makes cut faces read cleanly |
| **h. LOD by re-tagging, not re-meshing** (DOI `10.1109/tvcg.2007.1012`) | Changing the LOD threshold costs <5% of a re-extract | **Detail that follows your attention, live.** The game can re-budget quality every frame: push resolution into what you are looking at, pull it from what you are not. Scope zoom that *adds* geometry. Photo mode. Slow-motion that gets sharper. Dynamic quality that holds frame rate on weaker hardware without a menu setting |
| **k. SDF↔SDF collision** (DOI `10.1016/j.cagd.2024.102305`) | Penetration depth, contact points and normals directly from fields | Covered in 1.3 — **collision valid mid-carve**, and no collider generation at all for field-defined geometry |
| **l. Lengyel's one-scalar geomorph** | The failure rate is 0 on the target fields | **Pop-free LOD at half the storage** (6 bytes/vertex instead of 12). Covered in 3.6 |

---

# Part 3 — Harness results, and what they let you ship

Harnesses do not add features. They **convert a thing that works into a thing you can promise**, which
is what actually decides whether a mechanic ships.

## 3.1 The metamorphic relation suite

**Result if positive:** all eight relations hold, including chunk decomposition at power-of-two cell
sizes.

**What the player sees.** **No cracks. Anywhere. Ever.** Not "we tested seven fields and found none" —
the chunk-decomposition relation asserts that meshing a world in pieces gives the same surface as
meshing it whole, on *any* field. B-006/B-007 and M-128 already found the dual methods leaving gaps at
seams, and M-195 found a crack that the seam counter was **structurally blind to** — 0 edges on the
plane, 28 on each side. That is the class of defect a player finds by walking somewhere you did not
test.

**What it lets you ship.** Any extractor as the default, rather than only the one that happens to be
seam-safe. Right now `Extractor::DualContouring` on a chunked volume silently produces a cracked world
(M-128), which means the sharp-feature path — the crate's headline capability — cannot be the default.

## 3.3 The per-chunk latency tail

**Result if positive:** p99.9 is known and inside budget.

**What the player sees.** **No hitch when it matters.** M-124 already measured the amortised story
beautifully — and the control is the important number: the same queue *unbudgeted* costs **20.62 ms in
one frame** against a budgeted peak of **2.10 ms**. That is a dropped frame at the exact moment the
player set off the charge, which is the most-watched frame in the session.

**What it lets you ship.** Big, dramatic, multi-charge destruction as a *designed* moment rather than
something the tech team asks designers to keep small. A mean with a confidence interval structurally
cannot tell you whether a ten-charge demolition stutters; a p99.9 can.

## 3.6 Lengyel's open question, settled

**Result if positive:** the ray `p + t·n` always finds its coarse-mesh counterpart, failure rate 0.

**What the player sees.** **The world stops changing while you look at it.** LOD transitions morph
continuously instead of snapping — and at 6 bytes/vertex instead of 12, so it is affordable on a large
resident set.

**Honest note.** The prediction in the research doc is that it **fails on `gyroid`** (>1%) and succeeds
on `fbm_terrain` (≈0%), because a heightfield is star-shaped along the up axis and a triply-periodic
minimal surface is not. So the likely shipping answer is "geomorph on terrain, not on caves" — which is
still a result, and it is one Lengyel himself flagged as unknown in 2010.

## 3.7 The redundant re-mesh factor

**Result if positive:** the factor is near **1.0** — an edit costs only what it must, and output outside
the dilated support is byte-identical.

**What the player sees.** **Continuous cutting.** This is the difference between a tool that fires in
discrete pulses and one that *cuts*. If a 3 cm carve costs a chunk's worth of re-meshing, the game has
to quantise editing into steps large enough to be worth the work. If it costs 3 cm of work, it does not.

**The mechanic it unlocks:**

- **Beam and drill tools** that open a surface at the rate your hand moves.
- **Sculpting** at a granularity below the brush — fine shaping, smoothing, detail passes.
- **Ambient deformation.** Vehicle tracks that rut. Footprints in snow and mud. Water that carves a
  channel over minutes of play. Erosion. All of these are thousands of tiny edits per second, which is
  affordable only if each one costs its own size.

**And the correctness half is the real payload.** Any change outside the dilated support is a **bug** —
a piece of the world silently rewriting itself somewhere the player was not looking. Nobody has
published this instrument, and it is the one that decides whether the incremental path is trustworthy
at all rather than merely fast.

## 3.8 Pixels-of-pop

**Result if positive:** LOD switch distances derived so that p99 displacement is under one pixel.

**What the player sees.** **Nothing — which is the entire product.** M-121 measured the current pop at
up to **3.136 cells**, and recorded that its size is what decides whether it can be hidden by a fade, a
morph, or nothing at all. Converting that into projected pixels turns a number into a rule: *switch at
the distance where p99 pixels-of-pop < 1*.

**What it lets you ship.** LOD that is invisible without hand-tuning per biome, per field, per art pass
— derived rather than eyeballed. And it composes with 1.9: sub-voxel engraving needs a principled
disappearance distance or it boils at range.

## 3.9 / 3.10 The mutation-coverage map and the benchmark corpus

**What the player sees.** **Nothing, in either case.** The mutation map tells you which instrument
catches which defect class — a project artifact. The corpus is a paper. Both are worth doing and
neither should be justified on gameplay grounds.

---

# Part 4 — The honest ledger: what returns nothing visible

Four experiments produce no player-perceivable change even when fully positive. Saying so protects the
ranking of the ones that do.

| Experiment | Positive result | Player-visible effect |
|---|---|---|
| 1.4 Decoupled Fallback scan | 1.4× on a stage already at ~0.37 ms | **None.** Its value is avoiding a device-class crash on Apple hardware |
| 1.6 Workgroup-shared reduction | 14.5× on a reduction | **None.** Already banked by GPU-010a |
| 2.7 Comparable quality metrics | Numbers placed beside published tables | **None**, unless the NV(%) reconciliation surfaces a real defect |
| 3.9 / 3.10 Mutation map, benchmark corpus | Coverage matrix; a distributable corpus | **None.** Project and publication value only |

And two more that are near-invisible on their own but load-bearing underneath something else: **1.8
Morton** (~3% radius; its real content is a null result) and **2.5 off-surface deviation** (a canary for
2.2's convergence answer).

---

# Part 5 — The combinations, because none of these ships alone

The interesting mechanics need two or three results at once. This is the part a roadmap should be built
from.

| Mechanic | Needs | Why the combination |
|---|---|---|
| **Engraving / sub-voxel tool marks** | 1.9 (root cache) + 2.6 (`thin_slab`) + 3.8 (pixels-of-pop) | 1.9 makes the geometry possible, 2.6 tells you when it stops being representable, 3.8 tells you at what distance to dissolve it. Without the last two it looks superb at 1 m and boils at 20 m |
| **Sealed volumes, flooding, breakthrough** | 2.3 (χ oracle) + 2.4 (isotopy) + 3.1 (metamorphic) + round-2's dynamic connectivity | 2.3 proves mesh topology matches field topology, 2.4 proves the field's topology is the generator's, 3.1 proves chunking did not change it. Drop any one and "sealed" becomes a coin flip — and the failure mode is a drowned player |
| **Continuous cutting / erosion / tracks** | 3.7 (re-mesh factor) + 1.2 (edit-log proxy) + 3.3 (latency tail) | 3.7 makes a small edit cheap, 1.2 stops the thousandth edit costing more than the first, 3.3 proves no individual edit hitches |
| **Sharp tools on natural terrain** | 2.1 (quad split) + 2.8a (gradient-free DC) + 3.1 (seams) | 2.1 makes DC shippable on rough fields, 2.8a makes it work on *edited* fields, 3.1 makes it safe as a chunked default. Today all three are open and DC is not the default for exactly that reason |
| **Flight, vehicles, long sightlines** | 1.1 + 1.3 + 1.5 (parallelism) + 1.7 (sparse chunks) | Each is +12–26% radius alone; together they are roughly 2× radius and a working sky and underground |
| **Material-authentic destruction** | 2.8c+j (tet mesh, fracture modes) + 2.1 (clean input) + 1.3 (collider budget) | Fracture modes need a tet mesh, the decomposer needs intersection-free input, and the fragments need colliders inside a frame budget |
| **Persistent shared worlds** | 1.2 (edit-log proxy) + 1.5 (determinism under parallelism) + M-36 (already measured) | 1.2 stops history costing anything, 1.5 keeps clients bit-identical while syncing only the log. M-36's ordering guarantee is already in hand for same-kind hard edits |

---

# Part 6 — The one-line version

If **everything** lands, the game is a persistent, fully destructible world you can fly across, in which
your tools leave marks finer than the voxel grid, cuts have real corners, the engine knows which volumes
are sealed and tells you the instant you break through, materials break like themselves, and nothing
pops, cracks, hitches, or gets slower the longer you play in it.

If only the **top four** in the research doc's ordering land — the field-share bench, the metamorphic
suite, quad splitting, and the re-mesh factor — the game is a destructible world with sharp tool marks,
no seams, and continuous cutting. That is already a game nobody has shipped.

If only **1.9** lands, you have sub-voxel engraving, and that alone is a visual signature no voxel game
has.
