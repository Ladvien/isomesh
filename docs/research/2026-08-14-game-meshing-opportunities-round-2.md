# Game meshing opportunities — round 2

**Date:** 2026-08-14
**What's different this time:** the previous opportunity passes (2026-08-11) mined the literature
before the crate existed. This one cross-references the corpus against **107+ of our own measurements**,
which changes which opportunities are real. Two halves: a new unexplored direction, and the openings
visible in numbers we already own.

---

# Part 1 — The new direction: topology as gameplay data

## The inversion nobody has done

Every PCG paper in the corpus — the dungeon-generation survey, PCG Book ch3, Marahel, ASP-PCG,
cyclic-graph dungeons — is **generative and graph-first**: author a graph, then realise geometry from
it. Not one is **analytic**: extract the graph from a field that already exists.

A destructible voxel game already computes a scalar field over the whole world. Its topology — every
chamber, tunnel, handle and chokepoint — is derivable from data the game already has, and no engine
does it. Games fake all of it with hand-placed trigger volumes.

## The reframe that makes it buildable

**The full contour tree is the wrong structure for the questions games actually ask.** "Is this
sealed?", "did I just break through?", "is this a chokepoint?" are not all-thresholds queries. They
are single-threshold questions about the connected components and bridges of the **air sublevel set**.

That is **dynamic connectivity**, and unlike contour trees it is measured and it is cheap:

| | Measured |
|---|---|
| CUPCaKE (`10.48550/arXiv.2509.14433`) | **microsecond query latency, millions of updates/sec** |
| Query complexity | `O(log V / log log V)` |
| Update, spanning forest unchanged (the common case — most digging doesn't alter connectivity) | `O(log V)` depth / `O(log²V)` work |
| Update, worst case | `O(log²V)` depth / `O(log⁴V)` work |
| Bridges = chokepoints (2-edge-connectivity) | polylog **amortized**; polynomial worst case (`2001.00336`) |
| k-edge connectivity, k ≥ 3 | jumps to `Õ(√n)` — don't go here |

**So: dynamic connectivity as the per-frame queryable layer; contour tree / Morse–Smale as a
background, amortized world-understanding pass.** Those are different budgets and different tickets.

## Three mechanics current engines cannot do

**1. Breakthrough as a first-class engine event.** *"You just connected two previously separate
regions."* No engine knows this happened. The signal is exactly a union-find merge in the air sublevel
set. Agarwal et al. (`1406.4005`) supply the ready-made combinatorial vocabulary for how components
"appear, merge, split and disappear" as the field varies, with `O(log n)` per certificate failure.

**2. Sealed-volume as a queryable predicate.** Airtightness, flooding, containment, "am I trapped."
One connectivity query instead of a designer-authored volume. Water and gas that respect the space
they're in; base-building with real airlocks; AI that knows a room has exactly one exit. The merge
tree additionally answers *"sealed at what density?"* rather than only yes/no.

**3. Navigation and encounter graph from field topology, with no baked navmesh.** A destructible world
invalidates a baked navmesh the moment the player digs. Morse–Smale segmentation partitions the domain
by critical point; saddles and bridges **are** the chokepoints. Cover selection, ambush points, loot
placement and difficulty gradients fall out of structure the game already computes. The demand side is
already in our corpus — GameAIPro ch27 (tactical pathfinding), ch26 (tactical position selection),
ch31 (spatial reasoning), and Hale & Youngblood (`10.1609/aiide.v5i1.12376`) all derive tactical
structure from **geometry**, never from field topology.

## What is measured, and what nobody has measured

**The only hard 3D numbers in the corpus** — PLMSS, `10.1109/tvcg.2023.3261981`, Table 2 cross-checked
against the prose on p.10 to resolve OCR-scrambled headers:

| Volume | Task | 1 thread | 24 threads |
|---|---|---|---|
| Foot CT, 256³ (16.8M verts) | asc+desc segmentation | 4.40 s | **0.36 s** |
| Foot CT, 256³ | full pipeline | 9.21 s | 0.65 s |
| Miranda, 512³ (134M verts) | asc+desc segmentation | 39.61 s | 2.31 s |
| Miranda, 512³ | full pipeline | 181.72 s | 4.20 s |

For comparison TTK takes 214.50 s / 34.78 s on the Foot volume — PLMSS is ~50× faster.

**0.36 s for a 256³ recompute on 24 cores is 20–40× over a 16 ms frame.** Background pass only.

Three things nobody has measured, which is where the opportunity is:

- **Incremental Morse–Smale in 3D does not exist.** The 2023 paper recomputes from scratch every time.
  The corpus's one maintenance result (`1406.4005`) is **2-manifolds only** — `h: ℝ² → ℝ`, a
  triangulated terrain, not a voxel volume.
- **Contour trees have no absolute wall-clock in anything convertible.** Only relative figures: ≤10×
  OpenMP, up to 70× vs serial, "6× faster than TTK." Nobody publishes seconds.
- **Dynamic connectivity has never been run on a voxel lattice at game scale.** Every measured system
  was benchmarked on social/web graphs — Twitter 81K vertices, Stanford 280K. A voxel air-graph is a
  bounded-degree 6-connected lattice with 10⁶–10⁹ vertices. Bounded degree should help; sheer V may
  not. And **batching is untouched** — games edit thousands of voxels per explosion, not one.
  `2002.05129` (batch-dynamic trees) is the right tool, is in the corpus, and has never been pointed
  at this.

## The correctness gap that lands on this crate specifically

> **Nobody has established that a region sealed in the scalar field is sealed in the extracted mesh.**

A cell can be topologically connected in the field and still produce a watertight surface, or the
reverse, depending on the case table and the interpolant. Every paper found treats field topology and
mesh topology as interchangeable. They are not, and the gap sits exactly on the seam this crate
occupies.

**This should be a property test before it is a mechanic** — and it is a good one: extract a mesh,
compute connected components of the air sublevel set, compute connected components of the mesh's
complement, assert they agree. That test does not exist anywhere and would be publishable on its own.

## Acquisition status

Downloaded and indexed: **Topology ToolKit** (`10.1109/tvcg.2017.2743938`, 30 chunks).
Downloaded, catalogued, **conversion failed** (olmocr 0 pages, retryable): Data Parallel Contour Trees
(`10.1109/tvcg.2021.3064385`), Parallel Peak Pruning (`10.1109/ldav.2016.7874312`), Dey & Wang Reeb
Graphs (`10.1007/s00454-012-9463-z`).

**The two papers the whole question hinges on are unobtainable.** Tarasov & Vyalyi 1998, *Construction
of contour trees in 3D in O(n log n) steps*, and Safa & Wang 2014, *Maintaining persistence and contour
trees for time-varying functions on 2 or 3-manifolds* (OSU tech report, no DOI). Both surfaced from
Agarwal's reference list. Add to the hand-acquisition list.

Also absent with no OA route: Carr, Snoeyink & Axen `10.1016/s0925-7721(02)00093-7` (*Computing Contour
Trees in All Dimensions* — the foundational cost reference), Flexible Isosurfaces
`10.1016/j.comgeo.2006.05.009`, Edelsbrunner et al. Time-Varying Reeb Graphs
`10.1016/j.comgeo.2007.11.001`.

One corpus repair: `10.1111/cgf.12596` (De Floriani, *Morse Complexes for Shape Segmentation*) is
indexed at **3 chunks / 5,162 chars** — that's the HTML abstract, not the paper. Same landing-page
signature as the five already identified.

---

# Part 2 — Openings visible in our own numbers

These need no new literature. They are consequences of measurements already in `FINDINGS.md`.

## A. Subgrid MT's 70× is one cache away — and the ticket says so

**M-98** measured subgrid at **70× classic MT and 196× Marching Cubes**, and identified the constant
exactly: `6 tets × 6 edges × 16 samples = 576 field evaluations per cell` against Marching Cubes' 8
shared corner samples. ~72× predicted, 70× measured.

The row already names the fix and why it wasn't taken: *"every cell re-finds the roots on edges its
neighbours already found, deliberately, because identical endpoints through a deterministic root
finder is what makes conformity hold without a cache — a grid-edge cache is the obvious optimisation
and the redundancy is large, but it has a **correctness precondition** and is unmeasured."*

**The precondition is stateable and testable:** a grid-edge cache is safe iff the root finder is a pure
function of `(edge endpoints, samples)`. M-95 already establishes the boundary — *"same arguments, same
output, not same field, same output"* — so the cache key must include `samples`. That is a bounded,
falsifiable change to the single largest measured cost in the crate, and it is the difference between
subgrid being a demo and being shippable.

## B. The GPU result says the boundary, not the algorithm

**GPU-010b**: the `poll(Wait)` cost **0.375 ms of a 0.454 ms extraction** — 83% — while the four bytes
it was waiting on cost 0.033 ms. Removing the sync point made CPU time **flat at ~0.17 ms from 33³ to
129³**, against 0.170 / 0.245 / 0.722.

That is `V-9` and `V-10`'s ordering — stage placement and data movement dominating arithmetic —
reproducing on your own hardware, in your own code. It deserves an `M-` entry framed as *confirmation
of the speed thesis*, not as a GPU optimisation note, because it is the strongest in-repo evidence
that the extraction algorithm is not the term to optimise.

## C. Aliasing is a shippable feature, not just a defect

**M-72**: sub-cell features **alias** under coarsening rather than vanish — `thin_plate` goes
4,088 → 1,016 → 248 → **56** triangles across LOD 0–3, still 56 where the plate is a fraction of a cell
thick. The row's own framing: *"a feature that vanishes at a known distance can be faded, one that
disintegrates into a resolution-dependent scatter pops."*

That is a **gameplay-visible** finding, and it points at a mechanic nobody ships: because subgrid MT
resolves features the grid cannot, you can **choose** the distance at which a feature stops existing,
rather than discovering it. Carved detail that fades on a designed curve instead of shattering is a
visual signature no voxel game has, and it is a consequence of already-landed work.

---

# What I'd do with this

| | Action | Why now | Size |
|---|---|---|---|
| 1 | **The field-vs-mesh sealing property test** | Costs a day, tests a claim nobody has stated, and is a precondition for every mechanic in Part 1. Publishable alone | S |
| 2 | **Grid-edge root cache for subgrid MT** | Largest measured cost in the crate; the precondition is now stateable; M-95 already fixes the cache key | M |
| 3 | **Dynamic connectivity over the air sublevel set** | The cheap half of Part 1, on measured foundations. Also: nobody has run DC on a voxel lattice, so the benchmark itself is a contribution | M |
| 4 | Re-frame GPU-010b as speed-thesis confirmation in FINDINGS | Five minutes; it's the best in-repo evidence for the project's central claim | S |
| 5 | Morse–Smale as a background pass | Only after 3. 0.36 s/256³ means amortized, and the persistence threshold is a hand-tuned magic constant in every paper — for a game that is fatal, since **the threshold decides what counts as a room** | L |

The honest ranking argument: items 1–3 are all bounded, all falsifiable, and all rest on things already
measured. Item 5 is the exciting one and is the one with an unsolved parameter-selection problem
sitting in the middle of it.
