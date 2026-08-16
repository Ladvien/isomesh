# isomesh — BACKLOG

**Updated:** 2026-08-15
**Companions:** `CLAUDE.md` (rules), `FINDINGS.md` (what we know and how well),
`BACKLOG_ARCHIVE.md` (completed tickets + why they changed),
`docs/2026-08-11-implementation-brief.md` (the how),
`docs/2026-08-11-bevy-examples-catalog.md` (example detail), `docs/research/` (the why).

**176 tickets archived, 6 open.** Completed rows move to `BACKLOG_ARCHIVE.md` with their amendments
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
> BLOCKED: **on architecture, and the size was wrong — it is `M`, not `S` (2026-08-15).** The ticket assumed the fix is a pairing choice inside the cell, which is what the reference implementation does. Grosso 2017's actual rule is not that. Definition 3.2: *"A topologically correct triangulation across singular cell faces will not divide the surface into two branches. **The asymptotes of the hyperbolas at the singular face including the hyperbola center are part of the isosurface.**"* — and §4.2 makes the singular saddle an inner vertex, then *eliminates* the triangles whose edges lie on the singular face. Both cells sharing that face do the same, and the two patches join **through the saddle point on the face**. **That point is shared between two cells, and this crate has nowhere to put it.** `edge_vertices` is keyed on `(lower sample, axis)` — a grid *edge*. A face-interior vertex needs a `(lower sample, face)` slot, or the two cells emit coincident vertices with different indices and the index buffer carries a seam that only `weld` closes; Marching Cubes here does not rely on welding, and A-015's interior vertices are cell-local *precisely* because nothing else can name them, which is the opposite case. So the work is a new cache keyed on faces, not a branch in a decider. **Not urgent, and the measurement says why (M-220):** 0 of 1,838 ambiguous faces on all eight reference fields and 0 of 299,215 over 400,000 random cells — a singular face needs `v₀·v₂` and `v₁·v₃` bit-identical, which continuous `f64` fields do not produce. It stays open because a consumer feeding **quantised** density reaches it immediately, which is where Grosso's 8, 58 and 20 per CT volume come from, and that consumer is this crate's audience.
| ☐ | **A-020b** | **The disk triangulation for a six-saddle cell that is not a tunnel.** ~~Grosso does not give one; derive it.~~ **Re-scoped on the day it was written, and the premise is gone (M-231).** A-020 classified these cells — an inner hexagon with a contour past Corollary 6's bound of six — as `Topology::SeparateDisks`, and `extract` refuses them with `Error::UnresolvedSixSaddle`. The refusal is right and stays. What is wrong is the assumption that a **new triangulation rule** is what is owed. Two measurements: **continuous corner values produce zero such cells** in 11,354 six-saddle cells drawn from 2,000,000 random ones, and **every one that quantised values produce has a body saddle within `1e-12` of a cell face** — 261 of 261 across four quanta, no exceptions, against a background degeneracy rate among other six-saddle cells that swings between 8% and 79% with the quantum. A saddle *on* a face is Grosso 2017 §4.2's **singular case**, which is **A-002i**; these cells are singular faces that `has_inner_hexagon`'s strict `0 < x < 1` test admits because floating point puts the root a few ulps inside. So this ticket is **blocked on A-002i** and will most likely be closed by it rather than needing work of its own. `every_separate_disks_cell_has_a_saddle_on_a_face` pins both halves, and its continuous arm fails loudly if a non-degenerate one ever appears — which is the only event that puts a triangulation back in scope. **One option deliberately not taken:** widening `has_inner_hexagon` to a tolerance would reclassify these cells directly, and it is not done here because it changes the classification of every six-saddle cell in the crate and that is a design decision, not a bug fix. | S | A-002i |
> BLOCKED: **on A-002i, and on a measurement rather than on effort (M-231).** The cells this ticket exists to triangulate are singular faces, not a topological subcase — 261 of 261 have a body saddle on a cell face and continuous values produce none at all. A-002i owns the singular case and is itself blocked on architecture: the saddle sits on a face *shared between two cells*, and `edge_vertices` is keyed on `(lower sample, axis)` with no slot for a face-interior vertex. Deriving a disk rule here before that lands would be building the wrong thing carefully.

### 4b — Game-shaped

These use the algorithms the way a game does: chunked, edited, budgeted, collided against.

| | ID | Example | What it has to prove | Blocked by |
|---|---|---|---|---|

---

## Phase 5 — Measurement

| | ID | Ticket | Size | Blocked by |
|---|---|---|---|---|
| ☐ | **M-005** | **The Apple half of the family measurement.** M-001 landed `benches/family` and ran it on the Ryzen (M-282); the same run is owed on the M5, because six findings quote M5 figures that nothing has re-measured. **What it settles, and it is not cosmetic:** M-19's fitted intercept, M-20's *"4.75 ns/sample, 211 M samples/s"* marginal cost, M-22, and M-45's *"the M5 is 2.76× faster than the Ryzen on Marching Cubes at 256³"* — the last of which is currently **unquotable**, since its Ryzen half moved 1.74× and its Apple half did not. **Acceptance:** `cargo bench --bench family` on the M5 at the current commit, the CSV committed as `docs/measurements/family-<slug>.csv`, and those four findings amended against it. Note the counter columns will read `unavailable`: `perf_event_open` is Linux, so the Apple rows carry milliseconds only, and M-281 says a millisecond is comparable **only within one binary and one build** — so the cross-machine comparison must be made on `family` against `family`, never `family` against `resolution_sweep`. | S | M-001 |
> BLOCKED: **it needs someone else's working tree, which is their call.** `mac_air` has an `isomesh`
> checkout at `4369e3c` — over a hundred commits behind — with `BACKLOG.md` modified and two untracked
> docs in it. Fetching and checking out the current commit there would touch uncommitted work that is
> not this repository's to move. Say the word and it is a ten-minute ticket; until then the Apple
> numbers stay as they are, marked stale rather than quietly re-used.

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

**Note the interaction with A-018.** That ticket already established, on `noise_cavity`, that the
positional weld can *create* a non-manifold edge and that the subgrid validity suite therefore stopped
welding before judging (M-226). R-001 is the general form of the same mechanism; read A-018's archive
row before starting, because half the evidence is already there.

| | ID | Ticket | Size | Blocked by |
|---|---|---|---|---|

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
| ☐ | **A-022** | **Disambiguate the face in the dual path, as the primal path already does.** A-021 (M-276) measured that **all 314** of `noise_cavity`'s non-manifold edges under Surface Nets have exactly **4** crossed boundary edges on their shared grid face, against **2** on all 30,891 manifold ones — the *ambiguous face*, in the literature's own words. The dual quad walk emits a quad per crossed edge and so makes **both** choices at once, where `FaceAmbiguity::AsymptoticDecider` makes one for Marching Cubes. **Read `dualsimp_tvcg.pdf` (Schaefer, Ju & Warren) before writing anything** — whether MDC's own criterion is face-based or component-based decides whether this is a new rule or a bug in the existing one, and rule 5 forbids guessing which. **Acceptance:** the 314 goes to 0 on `noise_cavity`, `gyroid` and `fbm_terrain` under Surface Nets and Dual Contouring, with T-007's golden hashes regenerated in the same commit and the diff explained rather than accepted. | L | A-021 |
> BLOCKED: **the source is paywalled and not in the corpus.** Schaefer, Ju & Warren, *Manifold Dual
> Contouring*, IEEE TVCG 13(3) (2007), `10.1109/TVCG.2007.1012` — `paper_download` reports *"No
> open-access PDF found"*, and `distill_search` finds only a slide deck that cites it. The nearest
> thing present is Ju, Losasso, Schaefer & Warren, *Dual Contouring of Hermite Data*
> (`10.1145/566570.566586`), whose manifold criterion is about **octree simplification** — *"test
> whether the dual contour for each individual fine cube is a manifold"* — and says nothing about the
> base-level ambiguous face this ticket is about. **Rule 5 applies exactly:** a wrong disambiguation
> rule produces meshes that look fine and are subtly non-manifold, which is the failure this crate
> exists to avoid, so the rule is not being invented. A-021's diagnosis (M-276) stands and needs
> nothing from the paper; only the remedy does.

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
| ☐ | **R-009** | **Find the one-to-three cells per cross-section that face backwards.** R-008 (M-288) failed to locate 80% of the past-90° normals but bounded them tightly. Counts are **31 and 62 at n=33, 63 and 126 at n=65, 127, 254 and 381 at n=129** — every one an exact multiple of `n − 2`. The wedge is extruded along `z`, so the sampling repeats identically in every layer and any cross-section feature repeats once per layer: **there are one, two or three offending locations in the whole two-dimensional cross-section**, and they sit **11 to 94 cells from the crease** on a surface that is otherwise two exact planes, which Marching Cubes reproduces exactly (0.0000° at a 180° wedge). **So a plane is being meshed wrongly somewhere, at one specific place, reproducibly.** **Approach, and it is not another sweep:** dump the offending vertices' positions and their incident triangles for one configuration — `θ = 150°`, rotation 17°, n = 65, which has 126 of them — and read the cell. A-021 is the model: it found its answer by printing a face-count histogram for a **plain half-space**, not by widening a census. **Pre-register as P-17 once there is a hypothesis worth registering; there is not one yet, and inventing one to satisfy the protocol is the failure the protocol exists to prevent.** **Falsified by:** nothing, until it has a hypothesis — this is the locating step, and it is done when the configuration has a name and a constructed minimal fixture. | M | R-008 |

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
