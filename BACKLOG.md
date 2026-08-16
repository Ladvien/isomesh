# isomesh — BACKLOG

**Updated:** 2026-08-16
**Companions:** `CLAUDE.md` (rules), `FINDINGS.md` (what we know and how well),
`BACKLOG_ARCHIVE.md` (completed tickets + why they changed),
`docs/2026-08-11-implementation-brief.md` (the how),
`docs/2026-08-11-bevy-examples-catalog.md` (example detail), `docs/research/` (the why).

**187 tickets archived, 14 open.** Completed rows move to `BACKLOG_ARCHIVE.md` with their amendments
attached — read that before re-litigating a decision this project already made.

---

## How to work this backlog

1. Take the **topmost unblocked, unchecked ticket**. Don't cherry-pick interesting ones — the order
   encodes dependencies and the test harness exists so later work is cheap.
2. One ticket = one commit (or a short stack). Commit message starts with the ticket ID.
3. **Check the box in this file as part of that same commit.** This file is the state.
4. If a ticket can't be finished, leave it unchecked, add a `> BLOCKED:` line under it saying exactly
   what's in the way, and move to the next unblocked ticket. Do not half-finish and check the box.
5. If a ticket turns out to be wrong or to need splitting, edit it and say so in the commit.
6. **On completion, move the row to `BACKLOG_ARCHIVE.md`** with an indented annotation recording any
   amendment, deviation or falsified premise. The annotation is the point; the checkmark is not.
7. **New tickets** slot in by dependency, not by number. A ticket split after the fact takes a letter
   suffix (`T-005a`/`T-005b`, `A-002`/`A-002b`); a genuinely new one takes the next free number in its
   series even if that puts it out of numeric order (`A-015`).

### Definition of done — applies to every ticket

- Code compiles with no warnings. `cargo clippy -- -D warnings` clean.
- **`RUSTDOCFLAGS=-D warnings cargo doc --workspace --no-deps` clean.** A third of the lint job, and
  the third that clippy and fmt do not cover — it is what catches a doc link pointing at nothing, or
  at a private item. Added at A-002, which shipped one and found out from CI.
- `cargo test -p isomesh` green. `cargo fmt` applied.
- **No Bevy dependency and no Bevy in code under `crates/`** — comment-stripped manifests clean,
  non-comment `.rs` clean, resolved-graph bevy count 0. Non-negotiable — see `CLAUDE.md` rule 2. Prose
  explaining the wgpu-follows-Bevy pin is not a breach; the CI gate checks the three forms above.
- Any algorithm ticket also passes the T-001 validity suite. No exceptions, no "I'll add tests later."
- Any perf claim has a committed benchmark that produced it.
- Public items have doc comments. Anything with a sign convention, a coordinate order, or a winding
  order says so **in the doc comment**, not in a code comment.
- **`FINDINGS.md` updated in the same commit** if the ticket measured something, contradicted
  something written down, or earned a method rule. A measurement that only exists in a commit message
  is not retrievable six weeks later.

**Size key:** `S` ≈ one sitting · `M` ≈ a day · `L` ≈ multi-day, consider splitting.

---

## Phase 16 — The fracture substrate

**Added 2026-08-16, and placed above Phase 0 deliberately.** Rule 1 of this file reads top-down, so the
topmost unblocked row is what gets taken next; this is the current work front and the phases below it
are history. Nothing here replaces an existing ticket except **D-011**, which retires a premise.

**Phase 15's experimental protocol applies to every ticket in this phase** — **H** pre-registered as a
`P-` entry in `FINDINGS.md` *in the commit before* the measuring commit, a committed **Harness** behind
one documented command, named **Records** to `docs/measurements/*.csv`, an explicit **Falsified by**,
and the **FINDINGS obligation** discharged in the same commit. See the table under Phase 15 for the
exact requirements; they are not restated here.

### Why this phase exists — the reframe, in one paragraph

A downstream consumer (`bevy_autogib`, a plane-cut fracture crate) measured its own output for the
first time and found two independent causes, only one of them about manifoldness: a **non-convex
cross-section** breaks its centroid-fan capper (a closed manifold U-prism fails; a cuboid passes only
because it is *convex*), and **non-manifold multi-shell input** breaks loop recovery (22 open cut
edges across 12 shards). The literature's answer is not to repair either: **production fracture does
not cut the triangle soup at all.** Müller, Chentanez & Kim (`10.1145/2461912.2461934` — the NVIDIA
lineage behind PhysX Blast) cut a *volumetric convex decomposition* and carry the visual triangles as
a payload assigned to a cell. Because plane ∩ convex polyhedron = convex polygon, the centroid fan is
**provably correct for every cap** — which is why the cuboid scores 8/8 and is not luck. Sellán et al.
(`10.1145/3549540`) reach the same architecture independently.

Three consequences for *this* crate: convex decomposition stops being "the collider answer that blocks
nothing" and becomes **the cutting substrate**; the SDF/GWN backend drops off the fracture critical
path (S-009 is re-parented, D-011 retires its premise); and the boundary-loop triangulation problem
dissolves into a solved one — Shewchuk's PSLG, not a polygon.

**Union-first is ruled out by measurement, not only by Takayama et al.** Sacht et al. ran exactly this
experiment on interpenetrating character limbs and report the legs sticking together and the arms
sticking to the belly and head. For fracture that is a *correctness* loss — you lose the ability to
separate head from torso, which is the whole point.

### 16a — isomesh core

| | ID | Ticket | Size | Blocked by |
|---|---|---|---|---|
| ☐ | **A-026** | **Convex decomposition as a cutting substrate. RE-PRIORITISED — this was ranked last twice and that was wrong.** It was scoped as "the collider *answer*, not a dependency." Under the production architecture it is the thing you cut, and every downstream defect in a plane-cut fracture pipeline is downstream of not having it. ~~Use Convex Primitive Decomposition over the V-HACD/CoACD line.~~ **Reversed — see ✗21 and the DECIDED note below.** Use **V-HACD or CoACD**, which *partition* the interior: per CPD's own §2, prior ACD work *"remesh or voxelize the input to make it manifold, then **partition** the manifold mesh top-down along cutting planes"*, and partitioning is the property a cutting substrate needs. CPD's three cited virtues are real and are collider virtues; **enclosure is disqualifying here**, since a wrapper strictly larger than the solid cannot conserve its volume. Correction to an earlier claim in this repo's research: input cleanliness is a **quality axis, not an entry condition** — V-HACD and CPD require nothing; VisACD's 35% intersecting-hull figure describes the hulls CoACD *emits*, not the mesh it is *given*. **H:** decomposing a closed shell and plane-cutting the *cells* yields fragments that are closed, manifold, χ=2 and volume-conserved to 1e-3 — matching the convex-cuboid baseline — on input where cutting the soup does not. **Harness:** a torso+head fixture, per-shell decomposition, same 12-plane sequence. **Records:** per-fragment closed/manifold/χ/volume, cell count, decomposition wall-clock per shell. **Falsified by:** any proxy fragment reporting open cut edges — which would locate the defect in plane-cell intersection rather than in the shells. **FINDINGS:** `M-`. **Renumbered from A-030 as written in the 2026-08-16 brief**, per rule 7 — A-026 is the next free number in the series. A consumer backlog citing "isomesh A-030" means this ticket. | L | — |
> **DECIDED 2026-08-16 (user).** Route **(a)**: a decomposer that genuinely partitions the interior — **V-HACD or CoACD**. ✗21 stands and CPD is out as a substrate, for the reason the ticket itself named as a virtue: *"guarantees enclosure"* is **exactly what disqualifies it**, because an enclosing wrapper is by definition bigger than the shape it wraps, so cutting it yields fat pieces. The three CPD claims are real and remain real *for colliders*; they were read as virtues for cutting without checking what the method produces.
>
> **Müller 2013 is parallel reading, not a gate.** The decomposer interface is the same either way — mesh in, convex cells out — so integration does not change based on what that paper says. What it decides is the **cut-and-assign** step, which is **A-027** and is owed regardless. Running the two serially would cost a day for information that arrives anyway.
>
> **One reference point to keep on the desk while doing this: Diazzi & Attene's VolumeMesher** (`10.1145/3478513.3480564`, in corpus). It reaches convex cells classified internal/external **without** tidying the input first — tolerant of self-intersecting non-manifold soup with holes. **Almost certainly unusable as a dependency** (C++, and rule 3 keeps this crate at one dep), so it is not a candidate. Its value is as a *measuring stick*: it tells you what V-HACD's or CoACD's mandatory tidy-up pass is actually costing in fidelity, which is otherwise invisible.
>
> **Scoped 2026-08-16 from both papers, and the pick between the two is now made on a measured property rather than taste: CoACD, with merging OFF (V-36).** A substrate needs cells that *partition* the interior — overlapping hulls double-count volume, which is exactly what disqualified CPD. CoACD guarantees that and V-HACD does not: cutting solid meshes with planes *"results in flat boundaries between components. It ensures intersection-free convex hulls and avoids the defects caused by voxelization"* (§6.2). **But its own §6.5 merge post-process breaks the guarantee** — a merged pair's hull is the hull of their union and can reach into a third neighbour — which is where VisACD's *"merging produces intersecting convex hulls in 35% of cases"* comes from, and **merging is on by default**. The merge exists only to *"further reduce the number of components"*, so switching it off costs component count and not correctness — the same economy-versus-correctness split as Delaunay-ness in T-024b.
>
> **What this ticket therefore owes, and it is more than one sitting.** CoACD's pipeline needs: a 3D convex hull; surface *and interior* point sampling; Hausdorff distance between point sets; a plane-cut of a solid mesh; and MCTS over candidate planes. **It also assumes 2-manifold solid input** — *"we can convert imperfect input … by pre-processing with an off-the-shelf manifold conversion algorithm [Huang et al. 2018]"* — so the repair pass ✗20 called a quality lever is, for **this** method, a prerequisite. That is the concrete answer to what the tidy-up costs, and it is what VolumeMesher would have avoided. **Split before starting** — the convex hull alone is a ticket, and rule 3 means every one of these lands with no new dependency.
| ☐ | **A-027** | **Cut-and-assign: plane-cut the cells, carry the triangles as payload.** **Split out of A-026 2026-08-16** — the decomposer's interface is *mesh in, convex cells out*, and that is A-026 whichever method wins. This is the half that Müller, Chentanez & Kim (`10.1145/2461912.2461934`) actually decides, and it is owed regardless of A-026's outcome, so the reading runs **alongside** A-026 rather than in front of it. Two halves: **(1)** recursively plane-cut the *cells*, where `plane ∩ convex polyhedron = convex polygon` makes the cap provably a convex polygon and a centroid fan provably correct — this is why a cuboid scores 8/8 and it is not luck; **(2)** assign each input triangle to the fragment whose cell contains its centroid, splitting only the *straddling* ones against the plane. A triangle-plane split is exact and **needs no loop recovery at all**, which is what dissolves the capper problem rather than solving it. **H:** proxy fragments are closed, manifold, χ=2 and volume-conserved to 1e-3; the *render* fragments still carry nonzero open edges, and **that is correct, not a failure** — see T-023. **Falsified by:** any proxy fragment reporting open cut edges, which would locate the defect in plane-cell intersection rather than in the cells. **FINDINGS:** `M-`. | L | A-026 |
| ☐ | **T-022** | **A 2D constrained-Delaunay capper with exact predicates, replacing loop-recovery-plus-fan.** The decisive result is that **the loop is the wrong data structure.** Shewchuk's `Triangle` (`10.1007/bfb0014497`, in corpus) takes a **PSLG** — vertices and segments — not a polygon, and its own parenthetical answers the nesting question: holes are handled by a flood fill *halted at constrained edges*, which "saves both the user and the implementation from a common outlook wherein one must define oriented curves whose insides are clearly distinguishable from their outsides." That kills four failure modes at once: a **figure-eight cannot be constructed** (a self-touching vertex is just degree-4); **crossing segments** resolve by inserting the intersection vertex; **non-convex sections** need no star-shapedness anywhere; and **nested loops stay holes with no containment query**, which is what makes it robust to welded shared vertices. Pipeline: collect segments in the plane's 2D frame → weld with tolerance **relative to the model bounding box, not an absolute epsilon** → resolve crossings → CDT → flood-fill from outside the bounding box, halting at constrained edges → emit inside-labelled triangles only. **One asymmetry in our favour:** the CDT-existence pathology Diazzi & Attene name (*"the CDT is not guaranteed to exist for arbitrary input triangles"*) is **3D-only** — in a cut plane it never arises. **H:** on a U-prism, CDT+flood-fill matches analytic cap area to 1e-6 with zero inconsistently-oriented edges, where a centroid fan overshoots **by exactly the notch area**. **Records:** cap area vs analytic, oriented-edge violations, both methods. **Falsified by:** the fan *not* overshooting by the notch area — meaning the star-shaped diagnosis is incomplete. **FINDINGS:** `M-`, plus `✗` against "the capper is correct on manifold input" (it is correct on *convex* input). **Unblocked by T-024a alone (2026-08-16): Delaunay-ness is triangle shape, not cap area or topology, so `incircle` (T-024b) is a quality lever here rather than a gate.** **Predicates are now done — T-024a and T-024b both landed 2026-08-16, so nothing blocks this.** | L | — |
> NOTE 2026-08-16: **two scoping problems, neither blocking, and this needs splitting before it is started.** (1) **The falsifier cannot be evaluated in this repo.** *"the fan not overshooting by the notch area"* compares against a **plane-cut capper**, and there is none here — `subgrid/surface.rs:731`'s `fill_centroid_fan` is a different thing, filling a contour cycle inside a tet rather than capping a cut, and it is **measured intersection-free on every reference field at every resolution** (M-199), which is the paper's own guarantee holding. So the primary acceptance — cap area against analytic to 1e-6 — is testable here and the comparative half is not, until a capper exists to compare with. (2) **It is an `L` whose architecture is downstream of A-026**, which is blocked on the decomposer choice; a capper for an undecided substrate is speculative. **Suggested split, by the ticket's own pipeline:** `T-022a` the constrained Delaunay triangulator proper (PSLG in, triangles out, flood-fill labelling), which is self-contained, generally useful, and now has its predicates; `T-022b` the cap pipeline around it (plane frame projection, bbox-relative weld, crossing resolution, emit). Only `T-022b` depends on A-026's outcome.
| ☐ | **S-009** | **RE-PARENTED — on-demand GWN field.** The original justification was that Manifold Dual Contouring "queries where it needs to." **That is false for this codebase** — `extract` calls `self.sample(sdf, shape, origin, cell_size)` at `crates/isomesh/src/dual.rs:257`, which loops all N³ points into a buffer before anything else runs — so the dense N³ ray cost is the right price *for extractors*. The question survives for genuinely sparse consumers, and the real one is **internal**: F-005's empty-cell rejection by sphere tracing, and any future collision query against the field. **H:** on-demand GWN beats the batch path below a query-count crossover Q\* < N³; pre-register Q\*. **Falsified by:** no crossover existing above the point where batching is trivially cheaper. **Not blocked by any downstream consumer.** | M | F-005 |
| ☐ | **D-011** | **Retire the premise that Manifold Dual Contouring "queries where it needs to."** It is false for this codebase and now verified false at `dual.rs:257` (see S-009). Delete it from the `S-007`/`S-008` archive notes in `BACKLOG_ARCHIVE.md` and from the Ask-2 commentary wherever it appears, and record the retirement as a `✗` in `FINDINGS.md` — a falsified premise is never deleted from `FINDINGS.md`, only from the places that still act on it. **This is the D-003/D-004 truth-pass pattern**, applied to a claim rather than a figure. **Acceptance:** `grep -rn "queries where it needs to"` returns only the `✗` entry. | S | — |

### 16b — bevy_isomesh

Thin by design: the downstream fracture crate depends on `bevy` directly and does its own `Mesh`
handling, so the Bevy layer is not on its critical path. These exist so the 16a work is reachable from
an app without re-implementing the glue.

**The three rows below were written as B-008/B-009/B-010 in the 2026-08-16 brief. All three numbers are
already taken** by archived tickets that mean something else entirely (scratch pooling, the quickstart
example, publishing metadata), so they are renumbered to the next free block per rule 7. A consumer
backlog citing "isomesh B-010" means **B-014**.

| | ID | Ticket | Size | Blocked by |
|---|---|---|---|---|
| ☐ | **B-012** | **`Mesh` → triangle soup, and back.** Downstream consumers merge scene meshes into one soup by hand today. Expose the conversion once, handling `Indices::U16`/`U32`, non-`TriangleList` topology refusal, and missing `Float32x3` positions — the same `warn!`-and-skip discipline the consumers already use. **Acceptance:** round-trip a `Mesh::from(Cuboid)` and get 24 vertices back, not 8 — i.e. the conversion does **not** weld, because welding is the caller's decision (see B-014). | M | — |
| ☐ | **B-013** | **`proxy_cells` example.** Render A-026's convex decomposition as wireframe cells over the source mesh, with a slider for cell count and a readout of per-cell volume vs source volume. **This is the example that makes the Tier A/Tier B architecture legible** — nobody believes "cut the proxy, not the mesh" until they see the cells. | M | A-026, B-012 |
| ☐ | **B-014** | **Expose the merge predicate over Bevy attributes.** R-010's `MergePredicate` gates a weld on the link condition (`Lk u ∩ Lk v = ∅`); a Bevy consumer additionally needs to refuse a merge when normals or UVs differ, or a cube corner's three normals collapse to one arbitrary one. **Same mechanism, two predicates** — topological safety and attribute preservation compose. **Acceptance:** welding `Mesh::from(Cuboid)` preserves all 24 vertices under the composite predicate and collapses to 8 without it. | M | R-010 |

### Reading order, for whoever picks this up

1. **Müller, Chentanez & Kim 2013** — `10.1145/2461912.2461934`. §1–2 and the VACD section. The
   production answer; dissolves the capper problem as a side effect.
2. **Shewchuk 1996, *Triangle*** — `10.1007/bfb0014497`, in the corpus. The PSLG definition and the
   hole/concavity flood fill. Answers figure-eights and nesting together.
3. **Diazzi & Attene 2021** — `10.1145/3478513.3480564`. The only method whose *stated* input tolerance
   matches a glTF character: self-intersecting, non-manifold, disconnected, holes and gaps. Reference
   implementation exists.

Ten-minute runner-up: **Sacht et al.**, *Consistent Volumetric Discretizations Inside Self-Intersecting
Surfaces*, Figs. 10–11 — the picture of a GWN union welding a character's limbs to its torso.

---

## Phase 0 — Foundation and the test harness ✅ complete

All eleven tickets archived (I-001..I-004, T-001..T-008). The bet paid: every algorithm since has
been cheap to validate because the harness predated it.

---

## Phase 1 — The usual suspects

Each algorithm ticket is done when: T-001 reports **no unexplained violations** on all **eight** test
fields at three resolutions; T-004 determinism passes; T-005 covers it; and a benchmark exists.

> **The eighth arrived at A-002e, and it is the one that can fail things (M-208).** `noise_cavity`
> exists because none of the original seven produces a cell with an interior ambiguity — 0 of 68,385
> surface cells — so five pre-registered claims were properties of the fixtures rather than of the
> code. A-017, A-018 and A-019 own what it found. Expect a new algorithm to need a pinned census on
> this field where the other seven give zero.

> **Amended 2026-08-12.** This originally said "zero violations on all seven fields." M-4 falsified
> that as a universal gate: Surface Nets is *legitimately* non-manifold where one cell carries two
> sheets (48 edges on capped gyroid, 15 on fbm_terrain), and those counts are pinned as **non-zero
> assertions** precisely so they can't drift silently. A known defect with a pinned number and a
> ticket that owns it satisfies this gate. An unexplained one does not.
>
> **A-010 has since landed and closed that owner.** Surface Nets' and Dual Contouring's counts stay
> pinned as non-zero — they are properties of one-vertex-per-cell, not bugs — and
> `manifold_dual_contouring` is the entry that takes the zero. Its own residue is M-59's parallel-edge
> collapse, pinned at one edge on the ✗15 fixture and zero everywhere else, owned by O-16.

> **Re-scoped 2026-08-14 (A-002d).** A-002b was one `L` blocked on size. It is now a series, because
> two of the three things that made it large turned out not to be there. **(1)** Its route was Custodio
> 2019, whose non-manifold remedy is a *grid* preprocessing pass. Grosso 2016 (`10.1111/cgf.12975`)
> and Grosso 2017 (`10.1145/3095140.3095179`) — recorded as `PAYWALL` in `meshing-library-target.md`
> and in fact indexed in home-still since 2026-08-10 (V-29) — reach the same manifold result with
> interior vertices and **no grid pass**, keeping unambiguous cells on the existing table path.
> **(2)** *"Cell-interior vertices, which the grid-edge-keyed vertex cache has no slot for"* is true
> and is not a blocker: **A-015 already built that mechanism** — `table::CENTROID_BASE`, created per
> cell at `mod.rs:208-236`, uncached by design, already budgeted in the `u32` bound at `mod.rs:140`.
> What is left genuinely is the meshing, in the five tickets below.

| | ID | Ticket | Size | Blocked by |
|---|---|---|---|---|
| ☐ | **A-002i** | **The singular case — Grosso 2017 §4.2 and its Algorithm 1.** **Re-sequenced 2026-08-15 on a measurement, and no longer blocks A-002b (M-220): it is 0 of 1,838 ambiguous faces on all eight reference fields and 0 of 299,215 over 400,000 random cells.** A singular face needs `v₀·v₂` and `v₁·v₃` bit-identical, which quantised CT voxels collide into readily — Grosso 2017 counts 8, 58 and 20 per volume — and a continuous `f64` field essentially never does. So it cannot change any mesh A-002b's acceptance measures. **It is still worth doing**, because a consumer feeding `u8` density reaches it immediately and that is this crate's audience; the fixture must be *constructed* rather than sampled, as ✗22's was. Note also that this crate already gives ties a defined answer — `ambiguity::face_is_joined` resolves them to *separated* — so what is owed is topological correctness, not a missing branch, and `ambiguity` should not be touched: handle it inside the `trilinear` path, which is opt-in. A saddle sitting exactly *on* a face, where the standard asymptotic decider splits into two branches what is one surface. This is A-002b's own constraint (c): the 2 of M-165's 15,625 opposed configurations where a root of `F` rests on Δ's pole to within last ulps get a **defined** answer here rather than an inherited one. The reference implementation's shape is a per-face singular flag that then snaps the affected quadratic root to 0 or 1, plus an edge-coordinate comparison to choose the face pairing when the decider itself has no answer. Rare but real, and measured rather than assumed: Grosso 2017 Table 1 counts **8, 58 and 20** singular faces across three 512²×~700 CT volumes (tier V). **A second route to it, found at A-020b (M-231).** The `[9,3]` case-13 cells A-020 refuses are singular faces seen from the other side: 261 of 261 have a body saddle within `1e-12` of a cell face, and continuous corner values produce none at all. **That also sharpens this ticket's own reachability claim** — the 0-of-1,838 figure comes from a *bit-exact* product comparison, and 86–100% of those cells have a bit-exact singular face while the rest are the same configuration one rounding away, so the exact test undercounts the phenomenon. **A-020b is now blocked on this ticket and will most likely be closed by it.** | M | A-002h |
> **PROGRESS 2026-08-15 — detection landed, and the blocker is characterised rather than removed.** `trilinear::singular_face_mask` now says which of a cell's six faces are singular, `how_often_a_face_is_singular` reuses it so the census cannot drift from the extractor, and `a_singular_face_needs_quantised_data` pins the reachability at both ends (M-232): **0** singular ambiguous faces from continuous `f64` over 400,000 cells, **6,658** at quantum 0.1, **20** at 1/255 — the same order as Grosso's 8/58/20 per CT volume. `ambiguity` is untouched, as the ticket requires. **Two things remain, and the second is a rule-5 stop.** (1) The face-keyed vertex cache, which is now a determined piece of work rather than an open design question: Grosso 2017 §4.2 says *"three saddle points will lie on a singular face, but only **one** will be shared with the neighbor cell"*, so one slot per grid face is enough, and a grid face is named by its min-corner sample plus its normal axis exactly as an edge is named by its lower sample plus its direction axis — a structural mirror of `edge_vertices`, same size and shape. (2) **A third face state carried through `Contours`, which is the blocker and is not the cache (M-233).** Definition 3.2 requires a singular face *not* to divide the surface into two branches, so its four cut edges must meet at the hyperbola **centre** — a four-valent junction. `segment_links` takes `joined` as one bit per face: exactly two routings exist and both are permutations of the cut edges, asserted over all 384 (case, ambiguous face, bit) combinations. So the change is to the contour representation the whole trilinear path and A-002's 16,384-pair decider validation rest on, which is larger and more delicate than a second cache and means **this ticket needs splitting, not just re-sizing**. (3) The triangulation. §4.2's fewer-than-six-saddle arm is fully specified — singular saddles become inner vertices, then *"triangles containing edges of the contour which are on singular faces are eliminated"*. **Its six-saddle arm is not**: *"the other two points are **slightly moved** towards the interior of the cell"*, with no distance given, and the recovered reference is the 2016 code whose singular handling is the face-pairing choice rather than §4.2's construction — so no artefact supplies the constant. Deriving or bounding that displacement is what this ticket now turns on.
> **The artefact search is now closed too (V-35, 2026-08-16).** V-31 recovered the code both Grosso
> papers cite and found it is the **2016** one. Three further routes are empty: `github.com/rogrosso`
> has one public repository and it is a lecture course, a GitHub code search returns nothing, and
> `github.com/reproducibilitystamp/tmc` — a live mirror, so not a deletion — is `pushed_at`
> 2016-06-06, the same artefact. **The displacement constant is in no published artefact**, so the
> rule-5 stop is confirmed rather than unresolved, and deriving or bounding it really is what this
> ticket turns on.
> BLOCKED: **on architecture, and the size was wrong — it is `M`, not `S` (2026-08-15).** The ticket assumed the fix is a pairing choice inside the cell, which is what the reference implementation does. Grosso 2017's actual rule is not that. Definition 3.2: *"A topologically correct triangulation across singular cell faces will not divide the surface into two branches. **The asymptotes of the hyperbolas at the singular face including the hyperbola center are part of the isosurface.**"* — and §4.2 makes the singular saddle an inner vertex, then *eliminates* the triangles whose edges lie on the singular face. Both cells sharing that face do the same, and the two patches join **through the saddle point on the face**. **That point is shared between two cells, and this crate has nowhere to put it.** `edge_vertices` is keyed on `(lower sample, axis)` — a grid *edge*. A face-interior vertex needs a `(lower sample, face)` slot, or the two cells emit coincident vertices with different indices and the index buffer carries a seam that only `weld` closes; Marching Cubes here does not rely on welding, and A-015's interior vertices are cell-local *precisely* because nothing else can name them, which is the opposite case. So the work is a new cache keyed on faces, not a branch in a decider. **Not urgent, and the measurement says why (M-220):** 0 of 1,838 ambiguous faces on all eight reference fields and 0 of 299,215 over 400,000 random cells — a singular face needs `v₀·v₂` and `v₁·v₃` bit-identical, which continuous `f64` fields do not produce. It stays open because a consumer feeding **quantised** density reaches it immediately, which is where Grosso's 8, 58 and 20 per CT volume come from, and that consumer is this crate's audience.
| ☐ | **A-020b** | **The disk triangulation for a six-saddle cell that is not a tunnel.** ~~Grosso does not give one; derive it.~~ **Re-scoped on the day it was written, and the premise is gone (M-231).** A-020 classified these cells — an inner hexagon with a contour past Corollary 6's bound of six — as `Topology::SeparateDisks`, and `extract` refuses them with `Error::UnresolvedSixSaddle`. The refusal is right and stays. What is wrong is the assumption that a **new triangulation rule** is what is owed. Two measurements: **continuous corner values produce zero such cells** in 11,354 six-saddle cells drawn from 2,000,000 random ones, and **every one that quantised values produce has a body saddle within `1e-12` of a cell face** — 261 of 261 across four quanta, no exceptions, against a background degeneracy rate among other six-saddle cells that swings between 8% and 79% with the quantum. A saddle *on* a face is Grosso 2017 §4.2's **singular case**, which is **A-002i**; these cells are singular faces that `has_inner_hexagon`'s strict `0 < x < 1` test admits because floating point puts the root a few ulps inside. So this ticket is **blocked on A-002i** and will most likely be closed by it rather than needing work of its own. `every_separate_disks_cell_has_a_saddle_on_a_face` pins both halves, and its continuous arm fails loudly if a non-degenerate one ever appears — which is the only event that puts a triangulation back in scope. **One option deliberately not taken:** widening `has_inner_hexagon` to a tolerance would reclassify these cells directly, and it is not done here because it changes the classification of every six-saddle cell in the crate and that is a design decision, not a bug fix. | S | A-002i |
> BLOCKED: **on A-002i, and on a measurement rather than on effort (M-231).** The cells this ticket exists to triangulate are singular faces, not a topological subcase — 261 of 261 have a body saddle on a cell face and continuous values produce none at all. A-002i owns the singular case and is itself blocked on architecture: the saddle sits on a face *shared between two cells*, and `edge_vertices` is keyed on `(lower sample, axis)` with no slot for a face-interior vertex. Deriving a disk rule here before that lands would be building the wrong thing carefully.
| ☐ | **A-025** | **Manifold Dual Contouring is not manifold on `noise_cavity`, and the paper says it should be.** A-022 (✗19, M-290) obtained the source and falsified its claim: §3 says the uniform-grid dual *"is always a manifold because the original MC algorithm always constructs a manifold and the dual preserves the topology of the surface"*, and over eight fields at three resolutions **Marching Cubes measures 0 non-manifold edges under both face rules** while `manifold_dual_contouring` measures **143** with the crate's default table and **114** with the decider-modified one the paper specifies. The premise holds; *"the dual preserves the topology"* does not. **Every one of them is `noise_cavity`** — MDC is manifold on the other seven fields — and that is the field A-002e added because none of the others produces a cell with an **interior** ambiguity (M-208), which a *face* decider cannot see by construction. **H, to pre-register as P-17:** the residue is cells whose two sides resolve one shared ambiguous face to the same cycle pair **because of an interior ambiguity**, so it falls to zero on cells where `InteriorAmbiguity` changes the cycle set and nowhere else. **Two things to settle, and the second is a decision not a measurement.** (1) The mechanism, against `the_defect_count_is_predicted_from_the_grid_alone`, which already predicts the count from the grid — extend it to say *which* cells and check the interior test on each. (2) **The default.** `ManifoldDualContouring` defaults to `FaceAmbiguity::Separate`; the paper's construction is the decider-modified table, which is 20% better on `noise_cavity` and, per the module docs, *worse* on `gyroid` at 25³. Changing it re-baselines every golden hash. **Do not change it as a side effect of (1).** **Acceptance:** the mechanism named with a constructed minimal fixture (A-021's method, not a wider census), and the default either changed with the hash diff explained or left with the reason written down. | M | A-022 |
> PROGRESS 2026-08-16 — **P-17 falsified, and one candidate is off the list (M-291).** The residue is
> **not** the interior ambiguity. `Interior::Joined` is reported by **100% of ambiguous-face pairs** on
> `noise_cavity` — offenders and control alike, all four resolutions, both face rules — so the any-axis
> test has no discriminating power at all. Restricted to the sweep across the **shared** face it does
> discriminate and points the *wrong way*: under `Separate` the offenders carry the join **less** often
> than the control, 0.58–0.73 against 0.95–0.99. The harness reproduces the crate's pinned counts
> (30/64 and 8/40 at 17³/33³) before reporting anything new, and extends them to 53/26 and 49/25 at
> 49³/65³ — matching M-290's mesh-derived numbers from the other direction. **What is left is naming
> the mechanism**, and the next step is A-021's method rather than another census: the offending set is
> **26 pairs at 65³** and their sign configurations can be printed.
> PROGRESS 2026-08-16 — **the mechanism is bounded, and the bound is exhaustive (M-292).** All 4,096
> two-cell sign patterns, with every *consistent* joined-mask assignment on top — the two cells
> required to agree about the shared face. **512** share an ambiguous face; **18** offend under mask 0,
> which is exactly what `Separate` does; **476** offend under some consistent mask; and **0** offend
> under every one. So the defect is **never forced by the sign configuration**. That does not license
> "a face rule can fix it" — a rule reads the face's values and has none of this enumeration's freedom,
> and the decider still leaves 25–49 pairs per resolution. **Combinatorially always avoidable; with a
> rule that is a function of the shared face alone, not.** Anything that fixes it needs strictly more
> context than the face, which is A-017's two rejected alternatives. Two exact structures fell out: the
> default's 18 are precisely the `(1, 1)` cycle-count bucket, and `(1, 2)`/`(2, 1)` are the only
> buckets the mask does not always control (0.700 against 1.0000 everywhere else).
> PROGRESS 2026-08-16 — **the mechanism, constructed rather than sampled, and half the acceptance is
> met (M-294).** Two tests on a hand-built `4×4×3` lattice — **48 samples, no field**. On the same
> samples Marching Cubes measures **0** non-manifold edges and both duals measure **1**, carrying four
> distinct faces: ✗19 in a single fixture, and the manifold construction priced at nothing on the
> `(1, 1)` bucket, since splitting a cell by cycle cannot split a cell that has one. The sharper half
> is that scaling the shared face's two inside corners, **with every sign held fixed**, walks the
> asymptotic decider's saddle across zero and takes the defect with it — `−0.25` and `−1` separate and
> offend, `−4` joins and does not, at 20 triangles throughout. So the offending set is not a set of
> sign configurations at all, which is M-292 seen from the other side. Mutation-tested four ways.
> BLOCKED: **the second half is a decision, not a measurement.** `ManifoldDualContouring` defaults to
> `FaceAmbiguity::Separate`; the paper's construction is the decider-modified table (V-34), which is
> 20% better on `noise_cavity` and, per the module docs, *worse* on `gyroid` at 25³. Changing it
> re-baselines every golden hash. That is the crate owner's call and the ticket says explicitly not to
> make it as a side effect of the mechanism work.

### 4b — Game-shaped

These use the algorithms the way a game does: chunked, edited, budgeted, collided against.

| | ID | Example | What it has to prove | Blocked by |
|---|---|---|---|---|

---

## Phase 5 — Measurement

| | ID | Ticket | Size | Blocked by |
|---|---|---|---|---|
| ☐ | **M-005** | **The Apple half of the family measurement.** M-001 landed `benches/family` and ran it on the Ryzen (M-282); the same run is owed on the M5, because six findings quote M5 figures that nothing has re-measured. **What it settles, and it is not cosmetic:** M-19's fitted intercept, M-20's *"4.75 ns/sample, 211 M samples/s"* marginal cost, M-22, and M-45's *"the M5 is 2.76× faster than the Ryzen on Marching Cubes at 256³"* — the last of which is currently **unquotable**, since its Ryzen half moved 1.74× and its Apple half did not. **Acceptance:** `cargo bench --bench family` on the M5 at the current commit, the CSV committed as `docs/measurements/family-<slug>.csv`, and those four findings amended against it. Note the counter columns will read `unavailable`: `perf_event_open` is Linux, so the Apple rows carry milliseconds only, and M-281 says a millisecond is comparable **only within one binary and one build** — so the cross-machine comparison must be made on `family` against `family`, never `family` against `resolution_sweep`. | S | M-001 |
> BLOCKED: **on the machine being quiet, and the earlier reason for this was wrong (2026-08-16).**
> What was written here first was *"it needs someone else's working tree"* — `mac_air`'s `isomesh`
> checkout sits at `4369e3c`, over a hundred commits behind, with `BACKLOG.md` modified and two
> untracked docs in it. **That was never the obstacle**, and it was asserted rather than checked: a
> clone into a scratch directory does not touch their checkout at all, and the host is reachable, on
> AC, and carries the same toolchain this branch is built with (`cargo 1.96.1`, `rustc 1.96.1`).
> **What actually blocks it is contention.** `mac_air` has been running another job at a steady
> 42–48% of a core for over four hours — sampled six times a minute apart, load average `1.4–1.7`, no
> sign of ending. `family` is a **single-threaded, memory-bound timing**, and this ticket exists to
> make four findings quotable again; a figure taken beside a persistent competitor for last-level
> cache and memory bandwidth is not one. Worse, the release build needed first would take every core
> the other job is using. So the run is owed a quiet machine, not a working tree, and it remains a
> ten-minute ticket the moment there is one. Until then the Apple numbers stay marked stale rather
> than quietly re-used.

---

## Phase 6 — GPU (do not start before Phase 5)

The speed analysis is explicit that stage placement dominates the extraction algorithm by roughly an
order of magnitude. Which means GPU work is worth doing — and worth doing *after* you know your own
numbers, or you won't be able to tell what the port bought you.

**On the mesh-shader route (E-303), which has no ticket row yet.** Kept here because it was written against GPU-013 and outlived it — GPU-013 was archived on 2026-08-14 and this was never about that row. An earlier version said the route was blocked; it is not — see M-147. **The route needs no `unsafe` in this repository at all.** `isomesh-gpu` never opens a device (its API takes `&wgpu::Device`, GPU-001's rule), and **Bevy writes the experimental token itself**: `experimental_features: unsafe { wgpu::ExperimentalFeatures::enabled() }` at `bevy_render-0.19.0/src/renderer/mod.rs:335`. `WgpuSettings`' default priority is `Functionality`, which requests every feature the adapter advertises, so **Bevy's device already reports `mesh_shader=true multiview=true points=true`** on this machine, measured. E-303 is a Bevy example and gets a mesh-shader-capable device for free; `WgpuSettings.features` is there to force it explicitly if the default ever changes.

**The probe is load-bearing, not belt-and-braces** — an earlier version of this note implied otherwise and was wrong in the opposite direction from the blocked claim it replaced. The free device is *one of three branches*: `WgpuSettings::default()` consults `settings_priority_from_env()` first, so **`WGPU_SETTINGS_PRIO` overrides it**; under any priority other than `Functionality`, `features` starts at `wgpu::Features::empty()`; and `adapter.features()` is machine-dependent. It is also **contingent upstream**: Bevy's line carries `// SAFETY: TODO, see bevyengine/bevy#22082`, an admission that a justification is owed, so if that issue lands as opt-in the default path loses mesh shaders. Track it.

**On "graceful fallback", which needs one distinction rather than a ruling.** A demo that *detects* capability and, finding none, says so plainly and shows the compute path instead is a demo reporting a capability — one path, chosen by a measurement, with the choice visible. What the one-path rule forbids is the *library* silently substituting compute for mesh shaders so a caller cannot tell which ran. The first is what this ticket should build.

**Two things still shape the work.** WGSL mesh shaders are **Vulkan-only** — wgpu's own source says *"naga is only supported on vulkan; on other platforms you will have to use passthrough shaders"* (V-23) — so on Metal a caller supplies pre-compiled MSL and the composed-WGSL pipeline does not apply, making mesh shaders a fork in the shader path rather than a flag on it. And the ticket's own wording needs revisiting: *"graceful fallback"* is a second execution path for one feature, which the one-path rule forbids; the shape that survives is a **capability check that refuses loudly**, as GPU-007's probe already does. The remaining `unsafe`-shaped gap is only `isomesh-gpu::headless::Gpu` opening its *own* mesh-shader device, which is a test convenience rather than this ticket.


---

## Phase 7 — Documentation & packaging

Added 2026-08-14 from a literature pull, per the research-first rule. The sources that shaped these
tickets: *Effective Rust* Item 27 (crates.io is for people **choosing** a crate, docs.rs for people
**using** one — two pages, two jobs), Carroll's minimalist-instruction research (users act, they
don't read: anchor the first success in a real task, and treat error recovery as content rather
than warnings), and Prana et al. 2019 (`10.1007/s10664-018-9660-3`, in home-still), whose finding
is that the "Why" is the content category most READMEs lack. This repo's READMEs have the opposite
problem: the Why is superb and the on-ramp is missing. The falsification-essay voice stays; these
tickets put a doctested front door on it and make every claim current. House rules for the phase:
no hard line breaks in new prose, no performance number without naming its machine and CSV,
absolute URLs in anything crates.io or docs.rs renders, the Vibe Coded label stays on every README,
and every README code fence must be compiled by something.

| | ID | Ticket | Size | Blocked by |
|---|---|---|---|---|

---

## Phase 8 — Experiment infrastructure

**Added 2026-08-13 after a re-evaluation against three questions: is the crate ready for novel
experimentation, is it usable from the Bevy ecosystem, and will the harness make experiments
iterative?** Phases 0–6 built algorithms and proved them correct. Nothing in them was built for
*swapping a rule and measuring the difference*, which is the entire shape of the work the research
docs now point at. This phase is the cost of that gap, paid deliberately.

The evidence it is real, re-verified on 2026-08-15 against the current tree: `benches/shootout.rs` and
`src/property/extraction.rs` both hand-enumerate every algorithm by name — 26 and 12 references
respectively, with 9 more in `resolution_sweep.rs`, 10 in `extract.rs` and 2 in `stage_breakdown.rs`.
There is no library-level `Extractor` trait; the public traits are `MeshSink`, `Real`, `Sdf`,
`Shape3`, `ReferenceField`, and each extractor is an unrelated struct. **Adding algorithm #9 costs an
O(N) edit across benches, property tests and examples instead of O(1).**

| | ID | Ticket | Size | Blocked by |
|---|---|---|---|---|
| ☐ | **X-005** | **Give `Extractor` the global sample base its callers cannot supply, and decide whether that is worth the API break.** `extract_into` takes `origin: [R; 3]` and every implementation computes `origin + cell_size · local`. `ChunkLayout::world_of_sample`'s doc calls itself *"the single place a sample's world position is defined — everything else routes through it"*, and **no extractor does**: a chunk at a non-zero base reaches its far sample plane as `(o + h·base) + h·local` where its neighbour reaches the same plane as `o + h·(base + local)`, and those are equal by algebra and not by IEEE. **R-004 priced it (M-278).** Canonical reconstruction gives **0** unmatched seam-plane boundary edges at every spacing tried; what the crate can offer today gives 0 only at a power-of-two spacing and **63–348** at `0.1`, `1/12` and `1/14`, plus a hole 1.05–2.08 cells wide in 2 of 12 rows where an ulp flipped a sign. The crate's weld hides all of it (✗18) and an unwelded consumer — M-69's collider — gets it in full. **The shape that works is one path, not two:** replace the `[R; 3]` origin with a pair `(grid origin, integer base)` and compute `o + h·(base + local)`, which degenerates to today's behaviour at base zero. `TransitionCell::sample` already took exactly this route at A-011b, so the precedent is in the tree. **This is a decision, not a fix, and it is the reason this ticket is unstarted:** it changes the signature of the crate's central trait and every one of its call sites — eight extractors, five benches, the property suite and 32 Bevy examples. **Acceptance:** either the change lands with R-004's harness re-run and the offset arm gone from the crate entirely, or the ticket is closed with a written decision to keep the API and treat power-of-two cell sizes as a documented input contract. Do not ship both paths. **Blast radius, counted rather than estimated:** 7 inherent `extract` methods behind one `forward_extractor!` macro, 39 `origin: [R; 3]` parameters under `crates/`, and **294 call sites across 101 files** — 188 in 45 files under `crates/`, 106 in 56 files under `bevy_isomesh/`, which is a separate workspace with its own lockfile and CI. | L | R-004 |
> BLOCKED: **on a decision that is the crate owner's, not the implementer's — and the measurement that would settle it is already in (M-278).** Both answers are defensible and they are not close together. **(a) Take the break.** `isomesh` is at 0.0.5, pre-1.0, and the fix makes vertex sharing structural at every cell size instead of at half of them; the cost is 294 call sites and a signature change on the trait X-001 exists to stabilise. **(b) Keep the API and write the contract down.** The crate's own weld closes every hairline (✗18), so nothing a welded consumer sees changes; what is owed then is a documented input contract — *use a power-of-two cell size for a chunked world* — plus the 1.05–2.08-cell holes in 2 of 12 rows stated as a known limit rather than left in a findings file. **What is not acceptable is both**, which is what an added `extract_based` alongside `extract` would be. Ask before starting.

---

## Phase 9 — Usable by someone who is not us

Every example in `bevy_isomesh/examples/` demonstrates an *algorithm*. **None of the 32 shows a person
how to put a meshed SDF into their own Bevy app.** That is the first thing a prospective user looks
for, and it does not exist.

> **Re-scoped 2026-08-15 on measurement, and one ticket was deleted.** This phase was written against
> a base 107 commits behind, when neither crate was on crates.io and `bevy_isomesh` had no README.
> Phase 7's D-001…D-011 have since landed. **`I-005` (reserve the names on crates.io) is dropped as
> done** — `isomesh`, `isomesh-gpu` and `bevy_isomesh` are all published. The README, its `readme`
> key, the compatibility matrix and `CHANGELOG.md` all exist. What survives below is the residue,
> and each row says what was verified present rather than assuming the original scope.

| | ID | Ticket | Size | Blocked by |
|---|---|---|---|---|

---

## Phase 10 — Keeping the harness honest as it grows

`FINDINGS.md` is the most valuable artefact in the repo and it is now **387 KB / 945 lines / 231
measurements**, with no index. **The figures this phase was written against — 166 KB, 730 lines, 107
measurements — were already stale when it was written and the file has since more than doubled its
measurement count.** It is past the size at which anyone reads it end to end.

| | ID | Ticket | Size | Blocked by |
|---|---|---|---|---|

---

## Phase 11 — The field contract

**Added 2026-08-15 from the SDF corpus build-out.** The crate has an *input contract it never wrote
down*, and one reference field already violates it in the exact region where its defects live.

`csg_difference` declares `is_exact_distance() -> bool { true } // away from the seam` — a comment
admitting the invariant is false, on a function returning true. Marschner, Sellán, Liu & Jacobson 2023
(`10.1145/3610548.3618170`) name this object: a **Pseudo-SDF**, *eikonal almost everywhere yet not a
distance function*, with error **concentrated at seams** — the union's medial axis. That is exactly
where A-014d located `csg_difference`'s coincident polygons.

**And the error is one-signed.** `min` (union) never overestimates — a conservative lower bound, safe.
`max` (intersection, subtraction) **overestimates near concave seams** — the direction that lets a
tracer step through a surface and mis-places an interpolated vertex. `csg_difference` is
`max(box, −sphere)`. **It is wrong in the dangerous direction.**

**The load-bearing distinction for everything below:** `min`/`max` of 1-Lipschitz functions is
1-Lipschitz. **The Lipschitz bound survives arbitrary CSG; exactness does not.** So a field stays a
valid conservative bound forever, no matter how many brush strokes — which is what makes Phase 12
provably correct under unlimited player editing.

| | ID | Ticket | Size | Blocked by |
|---|---|---|---|---|

---

## Phase 12 — Exploiting the bound

Everything here rests on Phase 11's finding that **the Lipschitz bound survives editing.** These are
correct under unlimited player carving; nothing here assumes exactness.

| | ID | Ticket | Size | Blocked by |
|---|---|---|---|---|

---

## Phase 13 — SDF construction

The crate consumes fields and has never built one. Every ticket here also gives the harness a *second
source of truth* to check the first against.

| | ID | Ticket | Size | Blocked by |
|---|---|---|---|---|

---

## Phase 14 — Certificates and field harness

| | ID | Ticket | Size | Blocked by |
|---|---|---|---|---|

---

## Phase 15 — Research tickets

**Added 2026-08-15 after a full pass over `FINDINGS.md` followed by a literature check on the
patterns.** Three measurements recur as **mechanisms** rather than incidents. Those are the research
directions; the rest is history.

### The experimental protocol — mandatory for every ticket in this phase

A ticket here is an **experiment**, not a feature. It is done when the question is answered, including
when the answer is "no." Each one must carry all five fields, and **the hypothesis must be written
into `FINDINGS.md` as a `P-` entry before the measurement runs** — a prediction that first appears
after the number is known is not evidence, and this project has already caught itself writing
expectations into docs that measurement then disproved (✗1, ✗3, ✗14, O-14).

| Field | Requirement |
|---|---|
| **H** — Hypothesis | Falsifiable, numeric where possible, pre-registered as `P-n` in `FINDINGS.md` **in the commit before** the measuring commit |
| **Harness** | Committed code. Runs in CI or by one documented command. No throwaway probes — M-89's census had to be re-run because the first one wasn't committed |
| **Records** | Named metrics, to `docs/measurements/*.csv`, wired into T-011's regression baseline |
| **Falsified by** | The specific observation that kills H. **A ticket with no falsifier is not an experiment** |
| **FINDINGS obligation** | `M-` if measured, `✗` if a written claim died, `E×-` if the change was reverted (T-013's format). **Same commit.** A result only in a commit message is not retrievable in six weeks |

| | ID | Ticket | Size | Blocked by |
|---|---|---|---|---|


---

### 15a — Welding is a topology-destroying operation, and the predicate exists

**The strongest pattern in `FINDINGS.md`.** Five measurements, two independent algorithms, one
mechanism:

- **M-59** — *"The dual of a manifold surface is a manifold complex; the index buffer is where it stops being a manifold mesh."*
- **M-99** — *"provably manifold and my weld is what breaks it — the same mechanism as M-59, in a second algorithm. Unwelded: 0 non-manifold."*
- **M-96** — unwelded output has no topology to check at all (2,240 boundary edges / 896 triangles). **The weld is a precondition, not a tidy-up.**
- **M-93** — 30 reported self-intersections were *all* vertex-duplication artefacts.
- **M-48** — the edge cache "does not share everything."

So the weld is simultaneously **required** and **destructive**, and the crate has no theory of when a
merge is safe. **The literature has one, it is already in the corpus, and it is cheap.**

Dey, Fan & Wang (`10.48550/arXiv.1208.5018`, in corpus) give the **link condition**:
`Lk u ∩ Lk v = Lk{u,v}`. For two *non-adjacent* coincident vertices `Lk{u,v} = ∅`, so on a triangle
surface it reduces to a one-ring test:

> **merge (u,v) is safe ⟺ `Lk u ∩ Lk v = ∅`** — their one-rings share no vertex. **O(deg u + deg v).**

They also prove a k-way merge decomposes into **k−1 pairwise merges evaluated in the intermediate
complex** — so a bucket of ≥3 coincident vertices is **not atomic**, which is what R-002 is about.
Guéziec et al. 1998 (`10.1145/280953.281628`, acquired) state M-59's framing verbatim 28 years early:
*"Several manifolds can be mapped to the original non-manifold by identifying vertices."*

**What is unclaimed:** nobody states this predicate for **index-buffer welding of coincident vertices
emitted by an isosurface extractor**, and nobody publishes the measured rejection rate. That is the
contribution — modest, real, and a paragraph rather than a paper.

> **The predicate died, and this section's contribution claim died with it (2026-08-16).** R-001 ran
> it. P-8 is falsified in both clauses and the gated weld is recorded as **strictly worse than no
> gate** (E×4): across 56 configurations it removed **at most 4 non-manifold edges and added up to
> 791 non-manifold vertices**, taking `noise_cavity` + subgrid from 301 to **1,092** and `sphere` +
> Marching Cubes from 0 to **96**. The mechanism is the k-way sentence above, read the other way
> round: a bucket of ≥3 coincident vertices is not atomic, so refusing one pair of `k` leaves the
> representative a **bowtie** — which is why the damage lands in the vertex column while the edge
> column barely moves. The rejection rate did get measured; it simply does not buy what this section
> predicted it would. **R-010 is what survives** — the same hook, an equivalence-relation key instead
> of a pairwise test, and no topological claim attached to it.

**Note the interaction with A-018.** That ticket already established, on `noise_cavity`, that the
positional weld can *create* a non-manifold edge and that the subgrid validity suite therefore stopped
welding before judging (M-226). R-001 is the general form of the same mechanism; read A-018's archive
row before starting, because half the evidence is already there.

| | ID | Ticket | Size | Blocked by |
|---|---|---|---|---|
| ☐ | **R-010** | **A merge predicate on `Welder`, for attribute preservation — and only that.** **Read E×4 before starting**, and the note above it. The link-condition gate this section was written around is dead: measured at R-001, strictly worse than no gate, reverted. The reason is what makes *this* predicate a different object rather than the same one with a new key. The link condition is a **pairwise** test applied to a **k-way** coincidence class, so it refuses one member of a set that would otherwise merge whole, and the leftover representative is a bowtie. A composite key such as `(class, normal, uv)` is an **equivalence relation**: it partitions the class into complete sub-classes, every member of each merges, and no proper subset is ever refused — so it cannot reproduce E×4's failure, and it is the only instantiation this ticket ships. `MeshBuffer` carries no UVs and no vertex class, so the key is **caller-supplied**, which also keeps hard rule 1; the hook is the general form of `Welder::remap()`, already documented as how consumers move their own parallel data. **H:** splitting on a caller-supplied key moves no topology metric relative to the unconditional weld beyond the splits the key itself names, on all eight fields × all extractors. **Harness:** the P-8 bench shape in `benches/experiment_p8.rs`, both welds in one pass. **Records:** non-manifold edges and vertices, boundary edges, Δ vertex count and split count, both ways, to `docs/measurements/`. **Falsified by:** any topology metric moving where the key is constant — which would mean the hook itself, not the key, is doing something. **FINDINGS:** `M-` either way. | M | — |

---

### 15b — Coordinate reconstruction is the crack source, not the algorithm

**Second recurring mechanism. Three measurements, three different subsystems, one cause:**

- **M-32** — *"Chunk seams are bit-exact only when the cell size is a power of two."*
- **M-49** — *"`ChunkLayout::cell_of` inverts `world_of_sample` inside a cell and not reliably on its corner — M-32 in a second place."*
- **M-73** — *"a transition cell that computes its sample positions by offsetting from a face origin puts a hairline crack in the seam"* — its *"and no weld can close it"* is ✗18, falsified at R-004.

Every one is floating-point coordinate reconstruction, not extraction. **Nobody has published what
fraction of "seam cracking" in shipped voxel engines is this rather than algorithmic.**

> **R-004 answered that for this crate, and the split is clean (M-278).** The **algorithm** owns the
> whole visible budget — remove the transition cells and the seam opens to 32–184 boundary edges,
> 1.03–3.01 cells wide, identically under both arithmetics. The **arithmetic** owns the invisible one:
> `1.44e-15` world units against a weld epsilon of `h · 1e-4`, so it is 0 cracks welded and 63–348
> under bit-identity, with a 1.05–2.08-cell hole in 2 of 12 rows where an ulp flipped a sign.
> Canonical reconstruction takes every column to zero at every spacing; **X-005** is what it would
> cost to have it.

| | ID | Ticket | Size | Blocked by |
|---|---|---|---|---|

---

### 15c — Two mechanisms nobody has explained

> **One of the two is now explained, and it produced a third (M-279).** R-005 asked why the dual goes
> superlinear and the answer is **IPC**: Surface Nets runs 1.57× Marching Cubes' instructions and
> 5.24× its cycles, and the growth is a 16% IPC decline on an instruction stream that is flat per
> sample. The gather everyone suspected is `O(n²)` and the cost is `O(n³)` — a field with **no
> surface at all** costs the same to within 0.9%. What is left is *where* the IPC goes, which is
> **R-007**.

| | ID | Ticket | Size | Blocked by |
|---|---|---|---|---|

---

## Deliberately not in scope yet

Recorded so they don't get picked up early, and so it's clear they weren't forgotten.

- **`O-17` — how much does a grid-edge root cache buy?** M-98 measured subgrid Marching Tetrahedra at 70× classic MT, and the constant is field evaluations: 576 per cell at 16 samples per edge, against 8 shared corner samples for Marching Cubes. Every cell currently re-finds the roots on edges its neighbours already found, deliberately — identical endpoints through a deterministic root finder is what makes conformity hold without a cache. A cache keyed on the grid edge is the obvious optimisation and the redundancy is large, but it has a correctness precondition and the saving is **unmeasured**. Settle it by caching and re-running `cargo bench --bench extract`, with the golden hashes as the guard that the mesh did not change.
- **`O-18` — the `cycles()` recomputation in `subgrid/surface.rs`.** The 2026-08-14 review proposed restructuring extraction so the cycle set from `cycles()` (`surface.rs:252`, which allocates a `Vec<Cycle>` per call) is computed once rather than recomputed. It sits here rather than as a ticket because it is efficiency churn inside the actively-ticketed A-014 series, and because only the review's one-line summary is on record — re-derive the exact shape from the code before acting. Settle it the way O-17 says to: restructure, `cargo bench --bench extract`, golden hashes as the guard that the mesh did not change.
- Nanite-style mesh-space cluster simplification — the research concludes it can't be repaired
  edit-proportionally (no local validity certificate). Field-derived LOD is the bet instead.
- Networked/concurrent editing — **closed out, not deferred. The verdict is in and it is bounded.**
  O-4 asked whether brush operations commute. They do, conditionally: a run of same-kind *hard* edits
  reorders bit-for-bit — one result from all 40,320 orderings, all `Add` and again all `Subtract`
  (M-36) — because `min`/`max` select an argument rather than computing a value. Across an add/subtract
  boundary they do **not**: 11 distinct results, and the difference is *semantic*, so no storage format
  or arithmetic repairs it (M-37). Smooth union is worse still — 40,317 distinct results from 40,320,
  smooth-min being neither associative nor bit-commutative (M-38).

  So the coordination-free story survives inside a run and dies at every boundary, and **that is a
  protocol's problem, not this crate's.** isomesh's whole obligation was to make the truth available,
  and `BrushOp::commutes_with` already returns the honest answer rather than the optimistic one. A
  networking layer needs sockets, clocks and a session model — none of which belong in a `no_std` crate
  whose public API is `[f32; 3]`. Nothing further is owed here; reopen it as a real ticket only if a
  consumer turns up needing something the existing predicate cannot express.
- Neural / differentiable extraction (FlexiCubes, TetWeave) — different problem, different crate.
- Publishing real releases to crates.io. **`I-005 — reserve the name` is done: `isomesh 0.0.0` was
  published on 2026-08-12, and `0.0.1` on 2026-08-13 to carry the crate's README** and the name is held. `megamesh` was taken 48 hours before we checked it,
  which is what made this urgent rather than tidy. The placeholder is 82 files / 329 KiB compressed —
  source, benches, golden hashes and proptest regressions, nothing stray — and `0.0.0` is now burned
  permanently, which is the intent.

  What stays out of scope is a **real** release. That wants a `crates/isomesh/README.md` (the root one
  is outside the package directory and cannot be referenced from it, so the crates.io page currently
  shows only the one-line description), a version policy, and a decision about whether `isomesh-gpu`
  and `bevy_isomesh` publish alongside it. None of that is urgent now the name cannot be taken.

  **Amended 2026-08-14 (D-003): most of that has since happened.** `crates/isomesh/README.md` exists
  (`d9b8836`) and is the crates.io page; releases are CI-driven on version bumps; `isomesh` and
  `isomesh-gpu` are live at 0.0.4. What remains of this item is a version policy. The `bevy_isomesh`
  decision is made and ticketed: D-007 dresses it, D-010 publishes it at 0.0.4.
