# isomesh — BACKLOG

**Updated:** 2026-08-15
**Companions:** `CLAUDE.md` (rules), `FINDINGS.md` (what we know and how well),
`BACKLOG_ARCHIVE.md` (completed tickets + why they changed),
`docs/2026-08-11-implementation-brief.md` (the how),
`docs/2026-08-11-bevy-examples-catalog.md` (example detail), `docs/research/` (the why).

**128 tickets archived, 3 open.** Completed rows move to `BACKLOG_ARCHIVE.md` with their amendments
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
| ☐ | **A-002i** | **The singular case — Grosso 2017 §4.2 and its Algorithm 1.** **Re-sequenced 2026-08-15 on a measurement, and no longer blocks A-002b (M-220): it is 0 of 1,838 ambiguous faces on all eight reference fields and 0 of 299,215 over 400,000 random cells.** A singular face needs `v₀·v₂` and `v₁·v₃` bit-identical, which quantised CT voxels collide into readily — Grosso 2017 counts 8, 58 and 20 per volume — and a continuous `f64` field essentially never does. So it cannot change any mesh A-002b's acceptance measures. **It is still worth doing**, because a consumer feeding `u8` density reaches it immediately and that is this crate's audience; the fixture must be *constructed* rather than sampled, as ✗22's was. Note also that this crate already gives ties a defined answer — `ambiguity::face_is_joined` resolves them to *separated* — so what is owed is topological correctness, not a missing branch, and `ambiguity` should not be touched: handle it inside the `trilinear` path, which is opt-in. A saddle sitting exactly *on* a face, where the standard asymptotic decider splits into two branches what is one surface. This is A-002b's own constraint (c): the 2 of M-165's 15,625 opposed configurations where a root of `F` rests on Δ's pole to within last ulps get a **defined** answer here rather than an inherited one. The reference implementation's shape is a per-face singular flag that then snaps the affected quadratic root to 0 or 1, plus an edge-coordinate comparison to choose the face pairing when the decider itself has no answer. Rare but real, and measured rather than assumed: Grosso 2017 Table 1 counts **8, 58 and 20** singular faces across three 512²×~700 CT volumes (tier V). | M | A-002h |
> BLOCKED: **on architecture, and the size was wrong — it is `M`, not `S` (2026-08-15).** The ticket assumed the fix is a pairing choice inside the cell, which is what the reference implementation does. Grosso 2017's actual rule is not that. Definition 3.2: *"A topologically correct triangulation across singular cell faces will not divide the surface into two branches. **The asymptotes of the hyperbolas at the singular face including the hyperbola center are part of the isosurface.**"* — and §4.2 makes the singular saddle an inner vertex, then *eliminates* the triangles whose edges lie on the singular face. Both cells sharing that face do the same, and the two patches join **through the saddle point on the face**. **That point is shared between two cells, and this crate has nowhere to put it.** `edge_vertices` is keyed on `(lower sample, axis)` — a grid *edge*. A face-interior vertex needs a `(lower sample, face)` slot, or the two cells emit coincident vertices with different indices and the index buffer carries a seam that only `weld` closes; Marching Cubes here does not rely on welding, and A-015's interior vertices are cell-local *precisely* because nothing else can name them, which is the opposite case. So the work is a new cache keyed on faces, not a branch in a decider. **Not urgent, and the measurement says why (M-220):** 0 of 1,838 ambiguous faces on all eight reference fields and 0 of 299,215 over 400,000 random cells — a singular face needs `v₀·v₂` and `v₁·v₃` bit-identical, which continuous `f64` fields do not produce. It stays open because a consumer feeding **quantised** density reaches it immediately, which is where Grosso's 8, 58 and 20 per CT volume come from, and that consumer is this crate's audience.
| ☐ | **A-020** | **The tunnel triangulation Grosso does not define, and it is reachable.** The rule closes each contour edge by how many steps its endpoints are apart around the inner hexagon — one triangle for zero, two for one, three for two. **Three steps has no rule**: the paper gives none, and the authors' implementation runs `case 0/1/2` with no `default`, so it silently emits nothing and leaves a hole (M-228). `extract` now returns `Error::UnresolvedTunnel` there rather than the hole, because inventing the missing case is what rule 5 forbids. **Where it lives:** Marching Cubes' **case 13** — the four alternating corners, the only case with all six faces ambiguous — at particular face resolutions, giving a tunnel whose contours are **nine and three** vertices. That is also outside **Corollary 6**, which says a tunnel's contours are at most six and three, so either the corollary is narrower than it reads or the classification of these cells is. ~~**Settle that first**~~ **— settled on 2026-08-15, and it is the classification (M-229).** Flood-filling the cell's inside region on a 96³ grid and counting how many components its inside *corners* fall into: a genuine tunnel joins its same-signed corners through the interior and lands in **one** component, and both shipped tunnel fixtures do. The `[9,3]` cells land in **two** — two separate blobs, not a cylinder. **So Corollary 6 was right and was being read as a description rather than a test**, and the ring count admits a cell the corollary excludes. **The fix is therefore in `Contours::topology`, and the remaining work is to make the test the corollary implies.** Proposition 1's asymptote-side predicate is what the paper offers and its prose does not pin it down (V-31, now amended) — so derive it, or derive an equivalent from the corollary's own statement, and do **not** guess. A triangulation for a nine-vertex contour against a six-ring is needed only if some genuinely-tunnel cell also reaches three steps, which nothing has yet shown; check that before deriving one. `a_tunnel_can_span_three_hexagon_steps_and_is_refused` carries the fixture, found by search over rounded corner values. | M | A-002b |
| ☐ | **A-017** | **Manifold Dual Contouring is not manifold on a cell with an interior ambiguity.** A-010's headline property — *"`manifold_dual_contouring` is the entry that takes the zero"* — held on seven reference fields because **not one of them can produce the configuration** (M-208). On A-002e's `noise_cavity` it returns **30 non-manifold edges / 60 vertices at 17³** and **64 / 128 at 33³**, and its Euler characteristic parts company with Marching Cubes' (`−30` against `−96` at 17³), which falsifies **P-5** — pre-registered before A-010 ran. Both are pinned as whole censuses (`MDC_NON_MANIFOLD_CENSUS`, `MDC_CHI_CENSUS`) so they fail if they spread *and* if they vanish. **The exact factor of two in every row is the lead**: each offending edge carries exactly its own two endpoints, so this is one edge-level mechanism and not two. ~~**First question, before any fix:** is this a defect in this crate's construction, or the published guarantee not covering a cell whose interior the trilinear interpolant joins?~~ **Asked and answered on 2026-08-15, and it is neither (M-224).** Three measurements, each killing an explanation. **(a) Not tunnels.** Exactly *one* of the 30 offending edges at 17³ and one of the 64 at 33³ lies within `1.5h` of a tunnel cell; **all** of them lie within `1.5h` of an ambiguous cell, and the field has 193 and 502 of those against 3 and 2 tunnels — a rate of 13–15% of ambiguous cells. So interior topology is not involved and Schaefer, Ju & Warren's one-vertex-per-cycle argument is not what is failing. **(b) Not duplication.** Every offending edge carries exactly **four faces and four distinct triangles**, uniformly — the opposite of ✗17, where Marching Cubes' four faces were two cells each emitting the same triangle twice. Four distinct triangles on one edge is two sheets genuinely meeting. **(c) Not the face pairing.** `AsymptoticDecider` drops it from 30 to 8 at 17³ without reaching zero, and *introduces* 3 offending edges on `gyroid` at 25³ where `Separate` gives none — so no face rule is uniformly better and this cannot be fixed by choosing one. **The mechanism is now exact, and it is a limit of the algorithm rather than of this transcription (M-225).** An ambiguous face has **all four** of its edges cut; one vertex per cycle per cell means that when all four lie in *one* cycle on each side, all four dual quads connect the **same pair** of cell vertices — four quads on one dual edge, one triangle each, which is exactly the four distinct faces. The identity `non_manifold_edges == shared ambiguous faces whose four cut edges lie in one cycle on both sides` is asserted from the grid against the mesh and holds with **zero error** under both face rules: `Separate` 30 and 64, `AsymptoticDecider` 8 and 40. **What remains is the decision, not the diagnosis.** Three options, and none is obviously right. (1) **Document the limit** — Schaefer, Ju & Warren separate sheets *within* a cell and never claimed to handle two crossed edges of one shared face resolving to the same cycle pair; amend A-010's archive row, keep the census, and stop calling the entry the one that takes the zero. (2) **Split the cycle** on such a face, which needs a second vertex where the paper gives one and so is no longer their algorithm. (3) **Take the dual of the trilinear surface instead** — A-002b's rule already meshes those cells correctly, and the ambiguous face is exactly where its contours differ. Option 3 is the one this repo is now unusually well placed to try. Do **not** pick by face rule: the decider lowers the count on `noise_cavity` and raises it on `gyroid`. | M | A-002e |

### 4b — Game-shaped

These use the algorithms the way a game does: chunked, edited, budgeted, collided against.

| | ID | Example | What it has to prove | Blocked by |
|---|---|---|---|---|

---

## Phase 5 — Measurement

| | ID | Ticket | Size | Blocked by |
|---|---|---|---|---|

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
