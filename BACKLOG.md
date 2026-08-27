# isomesh — BACKLOG

**Updated:** 2026-08-26
**Companions:** `CLAUDE.md` (rules), `FINDINGS.md` (what we know and how well),
`BACKLOG_ARCHIVE.md` (completed tickets + why they changed),
`docs/2026-08-11-implementation-brief.md` (the how),
`docs/2026-08-11-bevy-examples-catalog.md` (example detail), `docs/research/` (the why).

**257 tickets archived, 21 open.** Completed rows move to `BACKLOG_ARCHIVE.md` with their amendments
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

## Phase 23 — twelve registrations from the 2026-08-26 audit, and the corrections it earned

**Added 2026-08-26, above Phase 22 for the reason every phase goes on top: rule 1 reads top-down.**
Phase 22 is closed. Nothing here supersedes Phase 17's or Phase 18's open rows.

**Source: `docs/research/2026-08-26-audit-and-phase-23-registrations.md`** — a skeptical audit of the
four Phase 21 experiments and the ✗43 entry, followed by twelve new pre-registrations drawn from
mathematics, formal logic and systems results the 2026-08-23 sweeps did not reach. The audit's own
framing is the reason the chores come first: five of the defects it found are the `✗35` failure mode
recurring, and one of them means several Phase 21 datasets correspond to no commit in this branch.

**Phase 15's protocol applies in full**, and all twelve `P-` entries are registered in
`crates/isomesh/src/experiment.rs` **before** any harness commit. Experiments are **bench-local**:
`crates/isomesh/src/**` is read-only apart from the registrations themselves, except where an
experiment *is* a source change and says so at registration (P-61, P-68's feature, P-69).

**The audit's central finding is a rule rather than a defect, and it changes how these are written.**
Four Phase-21 clauses could not have discriminated anything — P-58's C1 and C2, P-59's C2, P-60's C2 —
and every one was catchable before its harness existed. So each clause below names the rows on which
it *can* fire, and that set is shown non-empty from the CSV rather than assumed.

**R-059 — Is the crossing bit-exactly antisymmetric if it is stored as an offset from the edge midpoint? (P-61)**
`✗39` found bit-exact octahedral equivariance available for the **six pure axis permutations on plain
`marching_cubes`** and attributed the reflection failure to `a/(a−b)` and `b/(b−a)` being two divisions of
the same two values. That attribution is right and incomplete: `fl(b−a) = −fl(a−b)` exactly, so the
subtraction is innocent, and what breaks is the **anchor** — `cube.rs`'s parameter is measured from the
*lower* corner and a reflection swaps which corner is lower, so the correct reflected parameter is `1 − t`,
which is not `b/(b−a)`. Storing `d = ((a + b)/2)/(a − b)` as a signed offset from the edge *midpoint* makes
reflection a **sign flip** rather than the affine map `0 ↔ 1`, and floating point respects sign flips
exactly: `fl(a+b) = fl(b+a)`, halving is exact, `fl(b−a) = −fl(a−b)`, `fl(S/−D) = −fl(S/D)`.
**This one is a `src/` change and says so at registration** — `cube::edge_crossing` and its five
placements — and its C2 is a **cost** clause: `T-007`'s 216 golden hashes are rebaselined in the same
commit, and a hash that survives falsifies it.

**R-061 — Is `O-12` finite at 2¹⁸, and does one sweep settle it for Marching Cubes? (P-63)**
`O-12` — *"is Marching Cubes unconditionally manifold now?"* — is the oldest open question in the ledger,
and its own text says what would settle it: an exhaustive search over configurations spanning more than
two cells, or a proof that a cell-local cycle triangulation plus shared face segments cannot produce a
non-manifold **vertex**. `✗43` found the first counterexample and it was inside one cell; the question
stands because *"a vertex whose two face groups sit in different cells would be"* a third mechanism.
**The search space is much smaller than it looks.** Every Marching Cubes vertex sits on a grid **edge**,
so every face incident to it comes from one of the **four cells sharing that edge** — a 3 × 3 × 2 block of
grid corners, **18 corners, 2¹⁸ = 262,144 sign patterns**. That is not a sample, it is the whole space,
and it runs in seconds. **C2 is the fixture-can-fail control and it already exists**: inject `✗43`'s
pre-fix single-apex fan and require a non-zero count. **C3 is a necessary-condition sweep only** — a dual
vertex lives at a cell centre and its link involves 4³ = 64 corners, so the same 18-corner block cannot
decide the dual family and the entry must not claim it does. The full 2²⁷ dual sweep is a nightly gate and
a separate ticket, deliberately not registered.

**R-067 — Does restructuring the sample loop autovectorise, with bit-identity as the gate? (P-69)**
`core::simd` is nightly and staying nightly, so the lever is **autovectorisation**, and the measured prior
says that is enough: Wilcox's AArch64 study on 100k `f32` samples measured scalar 77.67 µs, hand-written
intrinsics 25.78 µs and **autovectorised safe Rust 25.54 µs** — safe code matched intrinsics. The patterns
that decide it are shape rather than machinery: struct-of-fields rather than index arithmetic, pre-slicing
once outside the loop so LLVM can prove the bound, and `chunks_exact`/`zip` iterators. `dual.rs`'s
`sample` currently `push`es into a `Vec` inside a triple loop, so the bound is re-proved on every element.
**C2 is the gate and it is the opposite of P-61's:** all 216 golden hashes **unchanged**, because IEEE
elementwise operations are exact per lane and a hash movement means LLVM reassociated something — so a
movement is a **defect** and the change is rejected, not rebaselined. Two facts established before the
harness: `libm::sqrtf` selects a hardware instruction under `target_feature = "sse2"`, so `sphere` can
vectorise, while `libm::sinf`/`cosf` have **no** arch selection at all — pure software with
argument-reduction branches — so `gyroid` cannot, at any loop shape, while `libm` is the transcendental
path. And the M5 that C1's threshold comes from is **contended**, which is `M-005`'s block, so the ratio
is measured here within one binary and one run (`M-281`) against this repo's own committed Zen 3 baseline.

**R-070 — Is the granularity of the active-cell structure a parameter, and where is its optimum? (P-72)**
`P-40` chose **64 cells per word** and never asked whether 64 was right; `M-337` measured the stage at
5.5× and 12/12 bit-identical, which settles that the bitmap works and says nothing about granularity.
GVDB (Hoetzlein, HPG 2016) is a **256× spread from one knob**: at 2048³, ⟨3,3,3,3⟩ builds 616,444 bricks
in 461 ms and ⟨3,3,3,6⟩ builds 2,036 bricks in **1.8 ms**, same data, and its own conclusion is *"larger
brick sizes produce a fewer number of bricks resulting in faster tree changes."*
**What the cubic 8³–64³ knob is in this crate, read from the source before registering:** it is the
**chunk**, not the word. `build_inside_bits` packs 64 cells per `u64` along **x only** — a flat per-row
word array with no block or brick layer — so there is no cubic granularity inside the bitmap to sweep, and
`u8`/`u16`/`u32`/`u64` would be a 1-D sweep of a packing width, not GVDB's knob. The unit that is
*rebuilt on an edit* is the chunk: `mark_edit` marks chunks, `DirtySet` holds chunks, `mesh_dirty`
re-meshes chunks. That is GVDB's leaf brick with the same semantics, and C1 is denominated in
**edit-plus-remesh** time, which only the chunk path has.
**The tradeoff is arithmetic and is why an optimum should exist at all.** Total world cells are held
fixed, so a chunk of `c` cells re-samples its shared corner planes: `((c+1)/c)³` of the field
evaluations, i.e. **1.42× at c = 8** against **1.05× at c = 64**. Against that, a finer dirty set
re-meshes fewer cells per edit. Small chunks pay in duplicated samples and save in wasted remesh.

| | Ticket | Size | Blocked by |
|---|---|---|---|

**R-060 — Is the Plantinga-Vegter certificate sound against this crate's own tunnel classifier? (P-62)**
`P-48` certifies a cell **empty** (`M-347`: zero unsound over 1.07e9 evaluations); `P-54` tightens it
(`M-354`: 3.85x more rejections on `gyroid`). Nothing certifies the other direction - *"this cell's
surface patch has no hidden topology"* - which is the difference between a mesher that is correct and one
that can **state where** it is correct.
**The predicate is already in the tree, so this is a measurement and not a build.**
`validate::isotopy::cell_is_certified` shipped under `T-015` and is exactly the registered form:
`0 in-not box-F(C) OR <box-grad-F(C), box-grad-F(C)> > 0`, with both bounds **exact** rather than interval
approximations, because the surface Marching Cubes approximates is the trilinear interpolant - `F` is a
convex combination of the eight corners, so clause one is *"all eight corners share a sign"*, and each
partial is bilinear in the other two coordinates, so its exact range is the min and max of four corner
differences. `h` cancels: the predicate tests the sign of a sum of three squares.
**What has never been done is the soundness check, and this crate owns a ground truth the PV literature
does not.** `A-020`'s classifier counts tunnels and twelve-vertex contours from the trilinear itself -
`M-214` recorded 2,053 and 173 in 396,000 cells, and `M-222` established that chi falls by exactly two
per tunnel. A cell containing a tunnel is a cell whose patch is **not** a graph, so a certificate on such
a cell is unsound. `M-214`'s counts are what make C1 a kill-shot rather than an `M-44` pass over an
unreached case.
**Registered caveat, not a discovery:** `C1` guarantees the patch is a *graph*, not that its planar domain
is connected, so the honest claim is *"no hidden topology in this cell"* and not *"exactly one
component"*. PV close that gap with a balanced octree this crate does not have; Lin & Yap document the
same gap (`10.1007/s00454-011-9345-9`).

| | Ticket | Size | Blocked by |
|---|---|---|---|

**R-062 — Can bounded model checking prove the combinatorics the property tests only sample? (P-64)**
`CLAUDE.md` rule 5 names this crate's correctness risk exactly: *"wrong case tables produce meshes that
look fine and are subtly non-manifold."* That is **combinatorics over eight sign bits** - 256 states -
and trivial for BMC, where bit-blasting IEEE 754 to SAT is the adversarial case a model checker is worst
at. So the split is: verify the combinatorics, keep testing the arithmetic.
**Kani** (`arXiv:2607.01504`) ranges over the eight sign bits and proves, for all 256 patterns, that no
case-table index goes out of range, no emitted index is at or past the vertex count, no triangle carries
two equal indices, and nothing panics.
**Both are dev tools with no runtime footprint, so hard rule 3 is not engaged**, and no published use of
either on geometry or graphics code was found - which makes this novel as well as useful.
**Scope note, registered:** neither tool touches vertex placement. Placement stays under proptest and
golden hashes. The honest scope is *"the table cannot be indexed wrongly"*, not *"the mesh is correct"*.

| | Ticket | Size | Blocked by |
|---|---|---|---|

**R-064 — Is a derivative sign test a usable under-resolution witness, scored against the root finder? (P-66)**
**This line died twice and the third attempt has to be a different mechanism.** `P-43`/x29 tried one
evaluation at the cell centre; `P-44`/x31 tried the mean residual. Both were **value** witnesses - they
asked whether the trilinear's value at a point disagrees with the field's - and the failure mode they were
chasing is not a value disagreement, it is a **missed root**: an edge whose endpoints share a sign while
the field crosses zero twice between them.
**The witness is a derivative sign test** from Finken, Li, Wang, Guo & Levine, *Topology-Preserving
Meshing of Implicit Scalar Fields via Monotonicity Constraints*, `arXiv:2608.12142`, IEEE Vis 2026 short
paper: if every edge of a PL mesh is monotonic with respect to `f`, the PL approximation is topologically
consistent with `f`'s critical points. The test is one line - sample the directional derivative at `k`
points along the edge, non-monotonic when any two disagree in sign.
**Do not port the paper.** It is explicitly 2D and the authors say so; its sampling-density argument, its
Theorem 1 case analysis (3D Morse theory has four critical-point types and a spherical link, not three
and a circle) and its separatrix refinement all fail to generalise, and it wants a Hessian this crate's
field trait does not expose. **Take the edge test alone, as a diagnostic rather than an extraction rule.**
**What makes this registrable here rather than anywhere else is the oracle.** `subgrid::roots::all_roots`
finds *all* roots along an edge - `M-94` resolved a slab at 1/1000 of edge length, `M-168` gave each
crossing an identity - so the root count per edge is a **known quantity in this repository** and the test
can be scored against it rather than against a hunch.
**Consequence if it holds:** a chunk can report *"this grid under-resolves this field, here"* as a number,
cheaply, per chunk - the missing input to the LOD decision `M-121`'s 3.14-cell surface pop and `M-72`'s
aliasing both want and neither has.

| | Ticket | Size | Blocked by |
|---|---|---|---|

**R-068 — Does subgroup ballot compaction beat the Hillis-Steele scan, and can it reach C1's number? (P-70)**
**The measured prior is on this hardware class and this API.** Smith, Levien & Owens, *Decoupled Fallback:
A Portable Single-Pass GPU Scan*, SPAA '25, `10.1145/3694906.3743326`: inclusive prefix sum over 2^25
elements on WebGPU/Dawn - M1 Max **1.43x**, M3 1.46x, RTX 2080 Super 1.49x, RX 7900 XT 1.33x, Mali-G78
1.35x, Intel HD 620 1.43x. Two of their findings foreclose the obvious plan whatever this returns: ARM and
Apple GPUs give **no forward-progress guarantee**, so plain decoupled look-back **times out on M1 Max and
M3**; and the ceiling is structural - reduce-then-scan moves O(3n), chained scan O(2n), so **50% is the
theoretical maximum**.
**C1's number was computed from the committed CSV before any shader was written, and it is arithmetically
unreachable.** `docs/measurements/gpu_vs_cpu.csv` at 129³ reads `gpu_total_ms` **8.3694**, of which
`scan_ms` is **0.3657** - **4.37%**. Reaching C1's 7.0 ms needs **1.3694 ms** removed. A **free** scan
leaves 8.0037; the literature's 1.5x leaves 8.2475. The residue is `upload_ms` at **7.324 ms**, i.e.
**87.5%**. This is Part 5's rule - *a clause stated as a ratio of a total must name the share of that total
it can move* - applied **prospectively** for the first time, which is what `M-375` earned it for.
**So the experiment measures the mechanism and does not land it.** A second WGSL path in the shipped crate
is what `CLAUDE.md`'s one-path rule forbids, and 4.37% capped at 1.5x does not buy an exception. The
subgroup scan is compiled **in the bench**, from inline WGSL, and its output is required to match both the
shipped `PrefixScan` and `cpu_prefix_sum` - three-way, so a transcription that drifted is caught by the
crate's own oracle rather than by reading two shaders side by side.
**Blocker resolved before writing anything:** naga validates `subgroup_invocation_id` only in **1-D
workgroup** compute shaders. `scan.wgsl` is already `@workgroup_size(256)` dispatched `(groups, 1, 1)`, so
no flattening is needed and its hash risk does not arise.

| | Ticket | Size | Blocked by |
|---|---|---|---|

**R-066 — Can a running error bound turn a vertex position into a position *and a certified interval*? (P-68)**
**The gap is stated in the crate's own mandate.** A CAD tool wants to know how much to trust a coordinate.
`M-142` found GPU and CPU agree on every triangle and disagree on 6% of vertices by exactly one ULP;
`M-144` found bit-identity is a property of the cell size, not the port; `M-30` found an unclamped solve
can fling a vertex 3.18 cells out of its own cell. **Every one of those had to be measured after the
fact.** Nothing reports, per vertex, at run time, how wide the interval containing the true crossing is.
**The construction is two extra flops per crossing.** A running error bound propagates a first-order term
alongside the value - Shewchuk's machinery (`10.1007/PL00009321`, already mined in five files) and the
filter hierarchy of Bartels, Fisikopoulos & Weiser (`10.1007/s10543-023-00975-x`, in corpus, uncited) -
but used to **report** rather than to **branch**. `P-61`'s centred form makes this cheap in a way the
parameter form does not: `d` is a single quotient with no cancellation, so its bound is
`|d| * (2 + |a-b|_err / |a-b|) * u` rather than the compounded bound a subtract-then-lerp accumulates.
**The ground truth is exact.** `a` and `b` are `f64` samples, hence dyadic rationals, so `(a+b)/2` and
`a-b` are exact integers over a common power of two and the true `d` is an exact rational. The error of the
`f64` result is therefore computable in `i128` with no floating point in the loop - which is what makes C1
a soundness statement rather than a comparison of two approximations.

| | Ticket | Size | Blocked by |
|---|---|---|---|

**R-063 — Does MCPro's procedural construction resolve `UnresolvedSixSaddle`? (P-65) — BLOCKED**
`Error::UnresolvedSixSaddle` is the one configuration this crate refuses to mesh; `M-231` found the `[9,3]`
cell is a **singular face** the strict interior test lets through, and `M-233` recorded the gap - a
singular face needs a **third routing** and the resolution mask has two. MCPro (Stahl & Grosso, GRAPP 2025,
`10.5220/0013309800003912`) says it built that routing, with **no lookup table at all**, and passes all
20,000 Etiene et al. cases on Betti numbers and Euler characteristic.
**Blocked at the door.** Six acquisition routes through home-still; the catalog entry is a **383-character
SciTePress landing page** (`conversion.server = html-parser`, `total_pages = 1`), no open-access PDF
(`M-371`). Running it would require inventing the quadrant subdivision, halfedge assembly and third routing
from an abstract, which rule 5 forbids. **Registered in `experiment.rs`; the harness waits on the paper.**
**Carry this if it ever lands:** the paper's own disclosure is that a trilinear isosurface **can genuinely
be non-manifold**, so `is_manifold()` failing on a singular face is *correct behaviour* - which changes
what the validity gate asserts. And `M-43`'s division-free, epsilon-free decider must be shown to survive,
or the golden hashes become scalar-dependent.

**R-065 — Does reduced affine arithmetic keep P-54's rejection rate in a fixed-size struct? (P-67) — BLOCKED**
`P-54` held (`M-354`: **3.85x** more cells rejected on `gyroid`) and left a structural problem: the form
grows a noise symbol per non-affine operation, so its size depends on the tape - an allocation whose size
depends on the scene, in a `no_std` crate. RAA folds all condensed error into one term, giving `[R; N+1]`
with no allocation and a fixed operation sequence.
**Blocked on C3, and splitting is refused.** C1 and C2 are runnable against this crate's own `P-54`
baseline. C3 exists to reproduce Knoll's measured **1.5-2x / 3-4x** band *and* his superquadric inversion
where intervals **win**, and `10.1111/j.1467-8659.2008.01189.x` has no open-access PDF (`M-371`).
Reproducing a band from a summary is `x21`'s failure. C3 is also the most informative clause, because the
exception is a **mechanism** rather than a measurement - so a P-67 without it would measure that RAA is
cheaper without ever testing where it is not. **Registered in `experiment.rs`; the harness waits.**

| | Ticket | Size | Blocked by |
|---|---|---|---|
| ☐ | **R-063** | L | — |
| ☐ | **R-065** | M | — |

**Both are blocked on acquisition, not on a ticket.** The `Blocked by` column names tickets and no ticket
can unblock these: `M-371` records that six routes were attempted and neither paper has an open-access PDF.
They are registered in `crates/isomesh/src/experiment.rs`, their questions are on record in `FINDINGS.md`,
and the unblocking event is external - a PDF appearing, or the authors publishing one.


**R-069 — Is the 83% a blocking round-trip, and can both targets avoid it? (P-71)**
`M-167` is the largest single number this project owns about its own GPU path: synchronisation was **83%**
of an extraction. `M-159` localised it — the last four bytes cost 0.033 ms to move and **0.375 ms to wait
for**, because `poll(Wait)` drains every dispatch queued before it — and `M-160` showed what removing it
buys: CPU time flat at ~0.17 ms from 33³ to 129³. **Two mechanisms already exist in the tree and this is
partly a measurement of them rather than a build:** `extract_buffers` waits once for the four-byte count,
and `extract_indirect` waits **not at all**, sizing from a budget and turning the total into indirect draw
arguments on the device. C2 is therefore the difference between two shipped entry points.
**Probed before registering, on this host's RTX 3090 / Vulkan:** `TIMESTAMP_QUERY`,
`TIMESTAMP_QUERY_INSIDE_PASSES` and `TIMESTAMP_QUERY_INSIDE_ENCODERS` are all **available**, so C1's
instrument exists — the crate currently requests `Features::empty()` and that is why `ExtractTimings`'
own doc says timestamp attribution "needs a device feature this crate does not request". C3's staging ring
is the only genuinely new capability, and its latency question is **the owner's**, not the harness's.

| | Ticket | Size | Blocked by |
|---|---|---|---|

---

## Phase 22 — the project has a site, and it is played in

**Added 2026-08-24, above Phase 21 for the reason every phase goes on top: rule 1 reads top-down.**
Phase 21 is closed. Nothing here supersedes Phase 17's or Phase 18's open rows.

**Source: a question — "is it possible to create actual WASM demos for these, and host them through
GitHub?"** It was, and half of it was already built and unmerged. What the other half cost is the part
worth recording: the roster went from three demos to nine, and the flagship could not run at all until a
defect in the *published* `bevy_isomesh` was found and fixed. `IsomeshPlugin` — the crate's entire public
reason to exist from Bevy — called `std::time::Instant::now()` in its frame-budget system, which compiles
on `wasm32-unknown-unknown` and then **panics on the first frame a chunk lands**. Nothing in the crate's
README or in the plugin's own documentation said "native only", so every browser game built on it broke
the moment it worked. That is `✗44`, and it was found by trying to put `game_showcase` on a web page
rather than by any test in the repository.

**The tenth demo is not a Bevy build**, and it exists because the nine make an argument the nine cannot
settle. Each Bevy module is 36 MB — about 8.8 MB gzipped on the wire — and essentially all of it is the
engine. `isomesh_web` is the same library with a 300-line hand-written WebGL2 renderer instead: **133,115
bytes**, eight fields, five extractors, and `isomesh::validate`'s report recomputed in the reader's
browser on every re-mesh. It is the front page, and it is the size claim made checkable.

**The gates are the deliverable as much as the demos are.** A demo built but not allow-listed is a 36 MB
module nothing can reach; one allow-listed but not built is a link to a 404; and the `site` CI job is
green either way. `doc_facts.sh` now derives the playable count from `build_web.sh`'s array and holds
three prose sites, `play.html`'s allow-list and every `#notes-` block against it — and each of those four
clauses was demonstrated failing before being left passing, because a gate never seen failing is
decoration.

**D-012 is closed.** The row is in `BACKLOG_ARCHIVE.md` with what it cost; the prose above stays here
because it is what the work was for.

| | Ticket | Size | Blocked by |
|---|---|---|---|

---

## Phase 21 — four registrations from the H1–H5 review and the corpus sweep

**Added 2026-08-24, above Phase 20 for the reason every phase goes on top: rule 1 reads top-down.**
Phase 20 is closed. Nothing here supersedes Phase 17's or Phase 18's open rows.

**Source: a critique of five hypotheses (H1–H5) plus a corpus sweep** that found two discrete-Morse and
shifted-linear papers sitting in the library uncited. Verification against the ledger and the source
**killed three of the five outright**: H1's tropical/min-plus framing is what `BrushStack`'s fold already
is and P-39 already exploits; H2's Clarke hull describes a computation `box_gradient` does not perform —
it returns *"the tied axis of lowest index"* (`fields/mod.rs:494-501`) — and M-66's constant already
carries an R-006 amendment saying it is not a property of the corner angle; H3's Krawczyk/interval-Newton
certificate is **foreclosed by name** at `docs/research/2026-08-23-discovery-dossier.md:294-298`, and
reopening named foreclosed ground is an owner's decision, not something to smuggle into a plan.

**Four of the review's own load-bearing claims were also wrong**, and the corrections cost nothing now
and a harness each later: M-24's 72/72 is 24 rotations × **3 hand-written Hermite crossing sets**
(`dual_contouring/solve/tests.rs:176-195`), not fields, and its own doc records that it is **structurally
blind to accumulation order**; the `20.1°–128.0°` figure is the span over 17³ *and* 37³ and M-283's
verdict is formally **reversed** by M-289; M-311's "792 dirty cells" is **925** by the entry's own table.
A certified cell-*emptiness* experiment also already exists as M-354, so nothing here re-derives it.

**Phase 15's protocol applies in full**, and all four `P-` entries are registered in
`crates/isomesh/src/experiment.rs` **before** any harness commit. Experiments are **bench-local**:
`crates/isomesh/src/**` is read-only apart from registrations, exactly as P-48 reimplemented interval
arithmetic and P-54 the affine form. A held result is evidence that a feature is worth landing; landing
it is a later ticket.

**All four are closed.** The rows are in `BACKLOG_ARCHIVE.md` with their verdicts; the prose below stays
here because it is what was predicted, and the archive annotations are what happened.

| | Ticket | Size | Blocked by |
|---|---|---|---|

**R-055 — Is a mesh bit-exactly equivariant under the octahedral group? (P-57)**
All 48 elements, not the 24 the existing generator at `dual_contouring/solve/tests.rs:130-165` filters to:
a signed coordinate permutation is exact in `f64`, so `mesh(g·f)` and `g·mesh(f)` are comparable
bit-for-bit. Compared as **sorted vertex-position multisets**, because the dossier already considered this
relation at `docs/research/2026-08-23-discovery-dossier.md:267` and held it — `table.rs` picks `safe_apex`
by **lowest edge index**, which is not invariant under axis relabelling, so a triangle-level statement
*"manufactures 2,688 false positives"*. C3 quantifies that warning instead of repeating it. M-178's
vacuous-fixture trap is the reason `fixture_can_fail` is a column and not an assumption.

**R-056 — Robins' `ProcessLowerStars`, with a chunk-local tie-break this crate has to invent. (P-58)**
Robins, Wood & Sheppard `10.1109/tpami.2011.95` is in the corpus, converted, embedded and **cited in zero
repo files**. Stage 1 is per-voxel local over a ≤27-cell lower star; stage 2 (`ExtractMorseComplex`) is
not local and is out of scope. **A fourth review claim was wrong and was checked before the harness
existed:** the review said the corpus markdown *"terminates mid-§4 before Theorem 11"*. It does not —
§4 is complete, Theorem 11 is stated in full, and so are Theorem 3, Theorem 6, Propositions 4–5, Lemma 12
and Algorithms 1–3. The clause C1 tests is the paper's own sentence, that the critical-cell census is
*"independent of this ordering"*. What the paper does not have is any test of that against **exact
ties**: Eq. (8) perturbs them away with a **global ramp depending on `I, J, K`**, which is
chunk-dependent and hash-breaking. This crate's fields tie *exactly*, so a chunk-local exact order is
**registered, not discovered**.

**R-057 — How many of P-39's 19 survivors are necessary? (P-59)**
M-341 measured a median **19 of 64** survivors with the mesh byte-identical on 64 of 64 chunks, and
**nothing measures the minimal correct survivor set** — no ablation, no per-brush necessity test.
Leave-one-out decides it. C1 is a soundness control reported first: if removing all non-survivors moves a
hash, the bound is unsound and every other number is void.

**R-058 — Does a shifted-linear reconstruction move a root, on one grid line? (P-60)**
Blu, Thévenaz & Unser `10.1109/tip.2004.826093` is in the corpus, named once in an inventory line
(`docs/research/2026-08-18-corpus-audit-and-procurement.md:128`) and never used. At `τ = 1/5` the causal
one-pole prefilter `c_n = −2⁻²·c_{n−1} + (1 + 2⁻²)·f_n` is **multiplication-free**. Scoped to a single
grid line and no extractor, because the prefilter is global with `(1/4)^k` decay — slab-local, not
chunk-local. C2 is a **pre-registered failure**: the paper states a Gibbs phenomenon on a step, and a
sharp CSG boundary is one.

## Phase 20 — the audit's registrations, with four citations corrected before a harness existed

**Added 2026-08-23, above Phase 19 for the reason every phase goes on top: rule 1 reads top-down.**
Phase 19 is closed. Nothing here supersedes Phase 17's or Phase 18's open rows.

**Source: `docs/research/2026-08-23-findings-audit-and-phase-20-registrations.md`** — an external audit
of the ledger against its own committed artefacts, plus
`docs/research/2026-08-23-phase-20-source-corrections.md`, which records what four parallel reads of the
primary literature found before any of this ran. The audit proposed eight registrations; one shipped
already as P-50/M-349, one is deferred, and **three of the remaining six rested on a claim the source
does not make**. Those corrections are the cheapest findings in the phase because they cost nothing to
make now and would have cost a harness each to discover later.

**Phase 15's protocol applies in full**, and all six `P-` entries are registered in
`crates/isomesh/src/experiment.rs` **before** any harness commit. Experiments are **bench-local**:
`crates/isomesh/src/**` is read-only apart from registrations, exactly as P-48 reimplemented interval
arithmetic and P-49 its union-find. A held result is evidence that a feature is worth landing; landing
it is a later ticket.

| | Ticket | Size | Blocked by |
|---|---|---|---|
| ☐ | **R-054** | S | |
| ☐ | **R-052** | M | |
| ☐ | **R-053** | S | |

**R-046 — Does either half of the tangent-sphere constraint survive this crate's extractors? (P-51)**
Sellán, Batty & Stein state both halves: the surface excludes every positive sphere and is *tangent to
every sphere at least once*. Count piercing vertices and untouched spheres as integers over the five
`FieldBound::Exact` fields. C2 is a **ratio** (DC ≥ 20× MC), not a per-1k bar, because M-27 already
measured ~150 per 1,000 and a bar of 20 could not fail. The untouched count is the half nobody has
measured and is a reference-free detail-loss signal if it is non-zero.

**R-047 — A tangency vertex rule, on this crate's own baseline. (P-52)**
Eq. (8) alone — one normalize and one fma per sample — as a third `VertexRule`, two iterations, clamped.
**Not their algorithm**: theirs is a global sparse solve with per-iteration remeshing over hundreds of
iterations, and its own Fig. 17 ablation measures clamping as detail loss. C4 is the clause that pays
either way: M-315 measured the QEF's centroid error *better than the perfect-placement floor* by 2.9–3.6×,
so a rule that pulls vertices onto spheres should spend that trade in reverse — and if it does not, the
QEF is not doing what M-315 says.

**R-048 — Custódio's third corner label, on the CT data that actually hits it. (P-53)**
The label is a **pure pre-pass** over the eight-corner classification; the paper's convex-hull
triangulator is not reproduced and not claimed. M-316 measured 3% of `bonsai` surface-cell corners
exactly on the isovalue. C3's half-offset isovalue is the control that says the label touches only the
equal-corner case.

**R-049 — Affine arithmetic where correlation actually lives. (P-54)**
Five stored reals, fixed size, no heap. The prediction is deliberately **non-uniform**: ≥1.5× more cells
rejected on `gyroid`, whose three trig terms cannot be extremal at once, and **<5%** on `box_exact` and
`csg_difference`, built from min/max for which the source gives no affine rule at all. Counted, not
timed.

**R-050 — A monotone-edge certificate for the two fields with no topological gate. (P-55)**
Finken et al.'s Theorem 1 is **proved in 2D** and its pigeonhole step has no hexahedral analogue, so this
is a labelled 3D port, not a transported proof. The tolerance is this crate's invention, fixed in the
registration and reported at two neighbouring values so its sensitivity is visible.

**R-053 — Ship the pinch test, not the label. (successor to R-048)**
M-352 measured the `=`-corner repair removing every degenerate triangle on both CT volumes and changing
the topology of one: 516 of 17,201 collapse groups on `bonsai` are pinches, welding 520 components,
against 0 of 50 on `fuel`. The decision is a graph property computable **before** the repair — a
union-find over the baseline triangles, asking whether two vertices snapping to the same corner already
share a triangle. Expose that as a `validate` report so a caller can ask whether the repair is safe on
its data, rather than shipping a repair that is safe on some data. The repair itself moves no geometry
(`max_snap_distance` is exactly 0), so this is a pure connectivity decision and belongs beside
`validate::sealing`.

**R-054 — `csg_difference` declares an underestimate ratio it does not meet. (from M-355)**
`fields/mod.rs:965` declares `Underestimate { q: 0.5 }`. P-51 measured ≈**0.11** at sample `(1.0625)³`,
where `|f| = h√3` points at a box corner the subtracted sphere has carved away, and three independent
extractors agree on the true distance to within 7%. Conservative rejection trusts declared bounds, so a
`q` that is optimistic by 5× is a correctness hazard and not a documentation nit. Either derive the
honest `q` for `max(box, −sphere)` or stop declaring one. Sweep the other composite fields for the same
defect while the harness exists.

**R-052 — The monotone-edge condition, on the complex it is actually about. (successor to R-050)**
✗36 proved P-55's zero unreachable: a mesh edge is a chord of the zero set, so the predicate is saturated
by geometry before it measures anything. Finken et al.'s PL function lives on the **ambient** complex.
Port it there — a fixed simplicial subdivision of the grid, monotonicity over tet edges including
diagonals, with the tolerance scaled by a quantity that is not identically zero at the endpoints. The
falsifier `box_exact` already handed over: its corner population is `O(n)` and halves per refinement,
which is what a resolution witness should look like. Gate on the contrapositive of Theorem 1 part 2b —
all-edges-monotone implies zero or ≥2 interior critical points — as an integer, where one failing cell
suffices.

**R-051 — The one vertex in 57,470, and whether it is a seam. (P-56)**
P-47's dead accuracy clause, re-asked as a bound: `(180° − θ)/2` where the six-sample stencil straddles a
CSG seam, and nowhere else. Tightness is recorded per fixture, because a bound that holds only by being
loose everywhere is not evidence.

---


## Phase 19 — six lenses nobody had pointed at a mesher

**Added 2026-08-23, above Phase 18 for the reason every phase goes on top: rule 1 reads top-down and this
is the current work front.** Phase 18 is closed. Nothing here supersedes Phase 17's open rows.

**Source: `docs/research/2026-08-23-discovery-dossier.md`** — eight parallel corpus sweeps aimed at
lenses the 2026-08-18 novelty table never used: non-cubic sampling lattices, certified root isolation,
bit-level and data-layout kernels, machine-checked combinatorial verification, discrete differential
geometry and digital topology, staged computation, kinetic geometry, and player-visible by-products.
**The whole dossier is tier R** except where it cites an M row, in which case the M row wins. Its
foreclosed list is as load-bearing as its candidates — eleven attractive directions are killed there with
the specific reason, and re-proposing one costs a sprint.

**Phase 15's protocol applies in full**, and all six `P-` entries are already registered — P-38 … P-43 in
`crates/isomesh/src/experiment.rs`, in the commit before any of them is measured. Committed harness behind
one documented command, named records to `docs/experiments/p-NN.csv`, explicit falsifier, FINDINGS
obligation in the same commit as the result.

**Ordering within this phase — this table governs rule 1.**

| Order | Ticket | Why here |
|---|---|---|

| | ID | Ticket | Size | Blocked by |
|---|---|---|---|---|
| ☐ | **F-009** | **Nothing gates a `FINDINGS.md` table against the CSV it names.** `doc_facts.sh` gates counts, `findings_index.sh` gates the index, `backlog_gate.sh` gates rows and `P-` reachability. The one systematic error Phase 19 produced — twice — is an entry quoting numbers from a run that is not the committed artefact (M-348), and no gate sees it. The check is mechanical for the subset that matters: an entry naming `docs/experiments/p-NN.csv` and printing a markdown table can have those numbers matched against that file's columns. **It would have failed on Phase 19 twice.** Scope it to tables whose header names registered record columns, so prose and derived quantities do not raise false positives. | M | — |
| ☐ | **F-008** | **The noise constants are private, and a certificate outside the crate has to transcribe them.** M-347's inclusion function for `NoiseVolume` and `FbmTerrain` re-implements `hash3`, `GRAD12` and `OCTAVE_OFFSET` because `fields/noise.rs` keeps them private, and it is guarded rather than trusted — a bit-exact comparison against `sample` at 137,842 points, max `|Δ| = 0`, runs before any certificate is issued. **The guard catches drift today and nothing prevents it tomorrow.** Two shapes: expose the three items (smallest change, but they are implementation detail and rule 3's spirit is against widening the surface for one consumer), or give `Sdf` an optional `enclose(lo, hi) -> (R, R)` with a default that returns `(-inf, +inf)` and let each field answer for itself — which is what P-48 wanted all along and is the shape `isotopy.rs`'s header implies. **Decide before anything else builds on M-347**, because both consumers of a transcribed constant are silent when it moves. | S | — |

---

## Phase 18 — Mechanics from the field

**Added 2026-08-17, and placed above Phase 17 for the same reason Phase 17 sits above Phase 16:** rule
1 reads top-down and this is the current work front. Nothing here supersedes Phase 17's open rows —
R-027a in particular stays live and stays cheap; if a sitting is short, it is still the best S in the
file.

**Source: `docs/research/2026-08-17-mechanics-from-the-field.md`** — five parallel corpus hunts
(surface-intrinsic computation, modal analysis, structural mechanics, shape semantics, volumetric
processes), 20 papers acquired, every candidate required to name the shipped game that already does it.
**The whole document is tier R** — five agents' reading, not this project's measurement; where a claim
touches an M-row, the M-row wins. Its two triage laws are worth reading before any ticket here: a
process is log-expressible exactly when its state is sparse, and candidates are classified by whether a
local edit has a *local answer*, not by whether the operator is cheap.

**Phase 15's protocol applies in full** — `P-` entry in the commit *before* the measuring commit,
committed harness behind one documented command, named records to `docs/measurements/*.csv`, explicit
falsifier, FINDINGS obligation in the same commit.

**These eight tickets are the dossier's Part 5, in its order: premise falsifiers first.** The dossier
itself expects three of the first five to come back negative, and that is why they are on top — a
negative here costs a day, and the same negative discovered during implementation costs a sprint (the
✗26 discipline). **Do not treat a null as a wasted ticket.** Build tickets for the Tier 1 candidates
(Calibre §1.1, speleogenesis §1.2, dynamic acoustics §1.3, the safe-to-dig field §1.4) are
**deliberately not written yet** — each is gated on its premise experiment below, and writing them
first would be building the wrong thing carefully (A-020b's lesson). The Tier 3 losers are recorded in
the dossier with the specific number that kills each one; they get no rows here on purpose.

**Ordering within this phase — this table governs rule 1.**

| Order | Ticket | Why here |
|---|---|---|
| 1 | **R-029** | Thirty minutes, and it re-tiers M-172's reading either way |
| 2 | **R-030** | The identity gate — three candidates die together if it fails |
| 3 | **R-031** | Validates the kinetics against published geomorphology before a voxel is touched |
| 4 | **R-032** | Decides whether the homotopy certificate is available at all |
| 5 | **R-033** | The cheapest possible kill-shot on the modal direction — run before assembling anything |
| 6 | **R-034** | The only external ground truth in the dossier |
| 7 | **R-035** | Decides prefactored-family vs Closest Point Method for everything intrinsic |
| 8 | **R-036** | Nothing left to falsify — cheapest audible ship, monetises R-022a/R-028 |

| | ID | Ticket | Size | Blocked by |
|---|---|---|---|---|

---

## Phase 17 — SOTA

**Added 2026-08-16, and placed above Phase 16 for the same reason Phase 16 was placed above Phase 0.**
Rule 1 reads top-down, so the newest work front goes on top. Phase 16 remains live and its topmost row
is blocked on a scope decision that is the crate owner's; nothing here supersedes it.

Every prior `R-` ticket is archived. These are the directions the research turned up that **never
became work** — verified before this phase was written: *self-adjusting*, *change propagation*,
*contour tree*, *persistent homology*, *dynamic connectivity* and *second fundamental form* returned
**zero** occurrences across `BACKLOG.md` and `BACKLOG_ARCHIVE.md`.

Each ticket states **the gap and why it is unoccupied**, a pre-registered hypothesis, the harness, the
falsifier, and **what it is worth if it holds**. Phase 15's protocol applies in full: `P-` entry in the
commit *before* the measuring commit; `M-`/`✗`/`E×-` in the same commit as the result.

**A note on tier.** These are tier **F** — hypotheses. Several are one negative result away from being
closed, and closing one is a finding. **Do not treat a null as a wasted ticket.**

**Ordering within this phase — this table governs rule 1, because the sections below are thematic
rather than sequential.**

| Order | Ticket | Why here |
|---|---|---|
| 1 | **R-024** | One day, publishable alone, and it gates R-022 |
| 2 | **R-026** | The result exists; only the writing is missing |
| 2a | **T-026**, **M-006** | R-026's two companions, re-created from prose references to tickets that never existed. Neither blocks it |
| 3 | **R-022a** / **R-022b** | Buildable now on measured foundations; the benchmark itself is unpublished. Split on V-41 — digging is insertion-only and cheap, filling needs a replacement search |
| 4 | **R-023** | Retires three blocked tickets at once, and nobody has even a negative result |
| 5 | **R-020** | The biggest space, and the one most at risk from unread prior art — get Acar's two papers first |
| 6 | **R-025** | Cleanest hypothesis, most likely to null out honestly |
| 7 | **R-021** | Highest ceiling, worst evidence base — two load-bearing papers unobtainable |

### 17c — Two things sitting on this crate's own seam

| | ID | Ticket | Size | Blocked by |
|---|---|---|---|---|

### 17d — The result you already have and have not claimed

| | ID | Ticket | Size | Blocked by |
|---|---|---|---|---|
| ☐ | **R-026** | **Write up the head-to-head.** **M-001 produced the comparison that does not exist in the literature**, and M-004's writeup ticket is archived while the paper is not written. Verified: **no paper since 2020 benchmarks Marching Cubes vs Surface Nets vs Dual Contouring against each other**, and Surface Nets — the thing engines actually ship — **has no credible published timings at all.** You additionally hold results that **contradict** published figures: M-51 and M-55 falsify the literature's `2–3×` Marching Tetrahedra ratio (measured `~3×` triangles for `4.3%` worse geometry, not 86%), M-1's `V_sn = V_mc + χ` identity, M-53's four-corner table of manifold × intersection-free, and M-54's `101×` Dual Contouring accuracy advantage on sharp fields. **This is the least speculative item in the phase and the only one whose result is already in hand.** The remaining work is Open SciVis volumes for comparability (H-005), mesh-quality metrics for the table reviewers expect (H-003), and prose. | L | — |
> **PROGRESS 2026-08-16 — the prose is delivered and the re-derivation found two of the document's own headlines were wrong.** `docs/research/2026-08-13-measured-comparison.md` is current: eight fields, seven extractors, one machine, every figure re-derived from a fresh `shootout` and `stage_breakdown` run, plus a new §6 carrying R-024's sealing audit. **Two conclusions did not survive.** ✗25 — *"Marching Cubes is alone in the good corner"* was true of a self-intersection detector whose straddle test counted tangential contact as a transverse crossing; the fix landed the day after the document was written, its commit message said it had *"inflated a metric this repo quotes"*, and nobody propagated it. M-308 — §4's *"the crossover is a property of one machine's cache behaviour"* is false on that same machine today, because **we optimised Surface Nets**; the comparison was measuring an implementation and reading as a statement about an algorithm. A third figure moved the same way: M-25's *"the feature-resolving solve barely registers"* at 3% is now 33%, because A-023/A-024 shrank the shared topology it was a fraction of.
> 
> **What is left is the two companions, and that is the whole of it.** T-026's AR>4 and sliver columns and M-006's real volumes are the two things a reviewer would ask for that the table still cannot show. The prose does not need rewriting again when they land — both are columns, not sections.
> **NOTE 2026-08-16: `H-005` and `H-003` do not exist and never have.** This ticket names them as its remaining work; there is **no `H-` ticket series anywhere in this repo** — not in `BACKLOG.md`, not in `BACKLOG_ARCHIVE.md`, not in `docs/`. The backlog gate cannot see it because the reference is in prose rather than in the *Blocked by* column. They are re-created below under the next free numbers in their own series, per rule 7: **`T-026`** is the mesh-quality metrics and **`M-006`** is the Open SciVis volumes. A consumer backlog citing "isomesh H-003" means T-026, and "H-005" means M-006.
> 
> **Neither blocks this ticket, and the scoping is deliberate.** The ticket's own claim is that the result is already in hand, and it is: `docs/research/2026-08-13-measured-comparison.md` is the draft, 325 lines of it, with every figure owned by a `FINDINGS.md` entry. What that document needs first is **re-derivation, not new data** — its own banner says the dual timings are high by ~4.26× since A-023/A-024, that it covers **five** extractors where seven ship, and **seven** fields where eight do. **R-024 has since added a result it does not carry at all** (M-307): the primal family seals all eight fields and all three duals leave the domain boundary open, which is a head-to-head axis no paper reports. T-026 and M-006 make the table more comparable; they do not gate the prose.

### 17a — The empty column

| | ID | Ticket | Size | Blocked by |
|---|---|---|---|---|
| ☐ | **R-027** | **A stable vertex naming, so the output is edit-proportional too.** **Split out of R-020 2026-08-16 on M-314**, which measured the gap: after a brush edit the *computation* is edit-proportional — 792 dirty cells at 33³, 65³ and 129³ alike — while the *buffer* is not, with **56–77%** of vertex slots changing for an edit touching **0.038%** of cells, and the ratio growing `O(n²)`. **The cause is a counter.** Vertices are appended in scan order and indices name buffer positions, so a cell emitting a different triangle count shifts every index after it. **The crate already has the stable name and throws it away**: the edge cache is keyed on `(lower sample, axis)`, which is a grid-global identity independent of emission order, and the packing step replaces it with a sequential slot. **H:** naming vertices by grid edge rather than by emission order takes `buffer_moved` to within a constant factor of `geometric_moved` — 330-ish rather than 15,706 at 129³ — with the *same triangles*, byte-identical after a canonical reorder. **Falsified by:** buffer churn that stays `O(n²)` under a stable naming, which would locate the instability somewhere other than the counter and would be the more interesting result. **This is what would make R-020's claim about a shipped extractor rather than about its dependency structure**, and it is the one thing standing between the measurement and *"the first isosurface extractor with an edit-proportional output"* — a claim V-43 says is available for scalar-field isosurfacing and for nothing wider. **Cost to weigh before starting:** a grid-keyed buffer is sparser than a packed one and every consumer's index buffer changes meaning, so this is an API question as much as an implementation one — read X-005's blast-radius note first, since it counted the same kind of cost for a different change. **FINDINGS:** `M-`. | L | — |
> **DECLINED AND SPLIT 2026-08-17 on V-45. The reason is harder than the API cost this row weighs.** M-318 found only one of three shapes delivers the 45×: a **persistent edge→slot map**, i.e. state carried across extractions. This row prices that as `extract_into` losing purity plus X-005's 294 call sites. **That accounting misses the binding constraint.** `validate::determinism.rs:268` runs `check_determinism` **three** times, the third into a **reused** buffer, under a doc comment saying why in as many words — *"to catch output that depends on the buffer's prior state… nothing else checks that it survives being driven that way."* R-027's only working shape **is** output that depends on prior state. So this does not cost a migration on top of a working design; **it converts a shipped gate's failure condition into its intent**, and T-004 is committed rather than preferred. Not a preference — a stop.
> 
> Reopen only on a formulation where the map is **derivable from the inputs** — a pure function of grid and field rather than of call history — which would keep the third run meaningful. Nothing in M-318's three shapes is that. **Blocked on R-027a**, which is the measurement this row should have had first: M-318 already says the encoding is not the cost, so locating where the 45× actually goes may dissolve the L entirely.
| ☐ | **R-027a** | **Where does M-318's 45× actually go?** **Split out of R-027 2026-08-17 on V-45**, which stopped the parent: its only working shape breaks T-004's reused-buffer run. This is the measurement that should have preceded the design. M-318 established the *ceiling* — a grid-edge naming takes buffer churn 15,706 → 346 at 129³, flat in `n` — and simultaneously established that **the encoding is not the obstacle**. Those two facts together mean nobody has yet asked where the 45× is spent: how much is vertices that genuinely move geometrically, how much is slots shifted by a predecessor cell's triangle count, and how much is order alone. **H:** the churn decomposes with the order-only term dominant, so a *canonical reorder at emission* — which needs no persistent state and so leaves T-004 intact — recovers most of the 45× without touching `extract_into`'s contract. **Falsified by:** a decomposition where the geometric and predecessor-shift terms dominate, which would mean only the persistent map can help and R-027 stays stopped rather than merely blocked — a result worth having either way. **Harness:** instrument the existing `edit_trace` bench; no new field, no new extractor. **FINDINGS:** `M-`. | S | — |
> **CEILING MEASURED 2026-08-17, before any API was touched (M-318).** A grid-edge naming would take buffer churn from **15,706 to 346** at 129³ — **45×**, flat in `n`, and equal to the true geometric change. **The prize is the whole gap**; M-314's split has no residue once the naming changes. Counted from the value arrays alone, keyed on the *edge* rather than on position, since position-keying would make the answer equal `geometric_moved` by construction and measure nothing.
> 
> **But the encoding is not the hard part, and that changes the scoping.** Three shapes and only the third works: a stable *order* does not help, because a crossing appearing still shifts every index after it; **index-is-edge-id** is stable and costs **230× the memory** (6.4 M slots for 27,822 vertices at 129³); and a **persistent edge → slot map**, allocating on first use and compacted occasionally, is what an incremental engine actually runs — and it is **state carried across extractions**, so `extract_into` stops being a pure function of its inputs. **X-005's 294 call sites are the smaller half of the cost.** The measurement says the prize is real; it does not say the shape is cheap.
| ☐ | **R-021** | **Maintain the contour tree, not the triangles.** **The reframe from `10.48550/arXiv.1406.4005`:** the maintainable object is scalar-field level-set topology, `O(log n)` per certificate failure, with certificates failing only on adjacent-vertex value swaps or saddle collisions — and it handles general update operations, not just continuous motion. **Caveat that must be carried:** that paper is **2-manifolds only** (`h: ℝ² → ℝ`, a triangulated terrain). The two 3D results the question hinges on — Tarasov & Vyalyi 1998 and Safa & Wang 2014 — are **both unobtainable**, and Edelsbrunner's 3-manifold Reeb maintenance is `O(n)` per certificate failure, asymptotically no better than rebuilding. **H:** the contour tree of a chunk can be maintained under a brush edit in time proportional to the dirty set, where full recomputation is not. **Falsified by:** maintenance cost tracking chunk volume. **Worth if it holds:** it changes axis 8 from "re-mesh fewer triangles" to "maintain topology, re-derive geometry" — a different algorithm class. **FINDINGS:** `M-`. | L | R-020 |
> **CORRECTED AND RE-BLOCKED 2026-08-17 on V-44, and the new reason is better than the old one.** **This row's caveat is contradicted by the paper it is built on.** `10.48550/arXiv.1406.4005`'s own related-work paragraph says a prior algorithm *"handles certificate failures in **`O(log(n))`** time"* and that *"their algorithm also works for **simple 3-manifolds** where the Reeb Graph is a contour tree"* — so *"`O(n)` per certificate failure, asymptotically no better than rebuilding"* cannot stand as written. The real weakness that paragraph names is the **event count**, not the per-event cost: certificates fail whenever *any* two vertices share a contour rather than only on adjacent-vertex swaps. Different objection, different remedy. **Tarasov & Vyalyi is findable** (`10.1145/276884.276892`, 60 citations, no open link) so *paywalled* is the accurate word, not *unobtainable* — and its title is **construction**, not maintenance, which is weaker support than this row implies. The 2-manifold restriction on the 2014 paper *is* confirmed.
>
> BLOCKED: **on a consumer, not on evidence — and that is the honest blocker.** R-022 called itself *"the cheap half of R-021"* and **that half is now delivered and measured**: M-311 has connectivity repair on a lattice at `O(|edit|)` with no logarithm, 792 dirty cells constant across a 64× lattice. R-022's own framing is that the questions a game asks — *is this sealed, did I break through, is this a chokepoint* — are **single-threshold**, and a contour tree answers **all** thresholds. Nothing in this backlog names a mechanic that needs more than one threshold at once. Add one, or leave this closed: **the ceiling is high and the demand is zero**, and M-314 additionally shows the extractor's trace is already edit-proportional with the bottleneck in the output encoding (R-027) rather than in topology maintenance.

### 17b — One disputed table, replaced by one scalar

| | ID | Ticket | Size | Blocked by |
|---|---|---|---|---|

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
> BLOCKED: **on how far the CoACD pipeline goes in-crate — a scope decision, asked 2026-08-16 and not yet answered.** Rule 3 forbids a new dependency, so every piece lands here or not at all, which makes "how much of it" a real choice rather than a detail. Two shapes: **(i) the whole pipeline** — 3D convex hull, surface *and* interior sampling, Hausdorff between point sets, solid-mesh plane cut, MCTS over candidate planes, plus a manifold-repair prerequisite — which is several tickets and the crate's largest single addition; or **(ii) a narrower first cut** — convex hull plus a greedy plane split, no MCTS, no collision-aware concavity — which yields partitioning cells of lower quality but is a fraction of the work and enough for A-027 to be built and measured against. **The measured facts are already banked** (V-36 and the scoping note above), so the answer changes only how this splits, not what it must achieve.
| ☐ | **A-027** | **Cut-and-assign: plane-cut the cells, carry the triangles as payload.** **Split out of A-026 2026-08-16** — the decomposer's interface is *mesh in, convex cells out*, and that is A-026 whichever method wins. This is the half that Müller, Chentanez & Kim (`10.1145/2461912.2461934`) actually decides, and it is owed regardless of A-026's outcome, so the reading runs **alongside** A-026 rather than in front of it. Two halves: **(1)** recursively plane-cut the *cells*, where `plane ∩ convex polyhedron = convex polygon` makes the cap provably a convex polygon and a centroid fan provably correct — this is why a cuboid scores 8/8 and it is not luck; **(2)** assign each input triangle to the fragment whose cell contains its centroid, splitting only the *straddling* ones against the plane. A triangle-plane split is exact and **needs no loop recovery at all**, which is what dissolves the capper problem rather than solving it. **H:** proxy fragments are closed, manifold, χ=2 and volume-conserved to 1e-3; the *render* fragments still carry nonzero open edges, and **that is correct, not a failure** — see T-023. **Falsified by:** any proxy fragment reporting open cut edges, which would locate the defect in plane-cell intersection rather than in the cells. **FINDINGS:** `M-`. | L | A-026 |
| ☐ | **T-022a** | **A 2D constrained Delaunay triangulator: PSLG in, inside-labelled triangles out.** **Split from T-022 2026-08-16** — the triangulator is self-contained and generally useful; the cap pipeline around it is not, and depends on A-026's outcome (that half is T-022b). **The decisive result is that the loop is the wrong data structure.** Shewchuk's `Triangle` (`10.1007/bfb0014497`, in corpus) takes a **PSLG** — vertices and segments — not a polygon, and its own parenthetical answers the nesting question: holes are handled by a flood fill *halted at constrained edges*, which "saves both the user and the implementation from a common outlook wherein one must define oriented curves whose insides are clearly distinguishable from their outsides." That kills four failure modes at once: a **figure-eight cannot be constructed** (a self-touching vertex is just degree-4); **crossing segments** resolve by inserting the intersection vertex; **non-convex sections** need no star-shapedness anywhere; and **nested loops stay holes with no containment query**. Scope here: segment input → resolve crossings → CDT → flood-fill from outside the bounding box, halting at constrained edges → emit inside-labelled triangles. **Predicates are done** — T-024a's `orient2d` is the correctness path and T-024b's `incircle` is the Delaunay quality lever. **One asymmetry in our favour:** the CDT-existence pathology Diazzi & Attene name (*"the CDT is not guaranteed to exist for arbitrary input triangles"*) is **3D-only**; in a plane it never arises. **Acceptance:** on a U-prism cross-section and a nested-loop cross-section, the triangulation covers the region with area matching analytic to 1e-6, zero inconsistently-oriented edges, and holes left empty. | L | — |
> NOTE 2026-08-16: **the split is done; what survives from the pre-split note belongs to T-022b.** That ticket's falsifier — *"the fan not overshooting by the notch area"* — compares against a **plane-cut capper**, and there is none in this repo. `subgrid/surface.rs:731`'s `fill_centroid_fan` is a different thing, filling a contour cycle inside a tet rather than capping a cut, and it is **measured intersection-free on every reference field at every resolution** (M-199), which is the paper's own guarantee holding rather than a comparison target. So T-022b's primary acceptance is testable here and its comparative half is not, until a capper exists to compare with.
| ☐ | **T-022b** | **The cap pipeline around T-022a.** **Split from T-022 2026-08-16.** Plane-frame projection of the cut segments → weld with tolerance **relative to the model bounding box, not an absolute epsilon** → hand to T-022a → lift the result back to 3D with consistent winding. **H:** on a U-prism, CDT+flood-fill matches analytic cap area to 1e-6 with zero inconsistently-oriented edges, where a centroid fan overshoots **by exactly the notch area**. **Falsified by:** the fan *not* overshooting by the notch area — meaning the star-shaped diagnosis is incomplete. **FINDINGS:** `M-`, plus `✗` against "the capper is correct on manifold input" (it is correct on *convex* input). **Blocked on A-026 in substance rather than in code:** under the Tier A architecture a cap is a plane ∩ convex cell, which is *provably* a convex polygon and needs no CDT at all — so whether this ticket is needed depends on which substrate A-026 lands. Do not start it before that is decided. | M | A-026, T-022a |

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
| ☐ | **B-013** | **`proxy_cells` example.** Render A-026's convex decomposition as wireframe cells over the source mesh, with a slider for cell count and a readout of per-cell volume vs source volume. **This is the example that makes the Tier A/Tier B architecture legible** — nobody believes "cut the proxy, not the mesh" until they see the cells. | M | A-026, B-012 |

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
>
> BLOCKED: **on a second mechanism, now characterised precisely, which is a change to the extractor rather than to the field.** The refusal moved `vertex 1` → `vertex 3`. **At a cell corner `frac = (0, 0, 0)`, so the analytic gradient reduces to the three *forward* differences alone** — `c₁₀₀ − c₀₀₀`, `c₀₁₀ − c₀₀₀`, `c₀₀₁ − c₀₀₀` — and is zero whenever those three neighbours match the corner. That is common in quantised data and is **not** a plateau: the `+`-side cell is not uniform at any of the three residual failures. **The extractor asks for a gradient at a cell corner**, where the interpolant's gradient is one-sided by construction. It snaps the *position* to the corner deliberately, to keep it bit-identical across cells (M-32, M-180); **the normal does not need that snap and is what breaks.** Evaluating the normal at the unsnapped crossing while keeping the snapped position is the obvious shape, and it is a change to `subgrid/extract.rs`'s position/normal split rather than a one-liner — and it needs a decision about whether two cells sharing a vertex may then disagree about its normal. **Counted against the denominator that matters:** the extractor evaluates a gradient at a corner only where the surface passes *exactly through* it, and on `bonsai` **16,284 of 529,508** surface-cell corners have value exactly zero — **3%**, because `u8` data with an integer isovalue lands *on* the isosurface constantly. **33 of those 16,284 also have a zero gradient**, so the surface passes through a **critical point of the field** and is singular there. **At such a point there is no normal** — not a missing one, an absent one — so refusing is correct and what is actually in question is only whether **33 singular points should refuse a 16 MB volume**. That is a granularity decision, not a correctness one, and it is the crate owner's.
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

## Phase 17 — SOTA

**Written 2026-08-16.** Every prior `R-` ticket is archived. These are the directions the research
turned up that **never became work** — verified: *self-adjusting*, *change propagation*, *contour
tree*, *persistent homology*, *dynamic connectivity* and *second fundamental form* return **zero**
occurrences across `BACKLOG.md` and `BACKLOG_ARCHIVE.md`.

Each ticket states **the gap and why it is unoccupied**, a pre-registered hypothesis, the harness,
the falsifier, and **what it is worth if it holds**. Phase 15's protocol applies: `P-` entry in the
commit *before* the measuring commit; `M-`/`✗`/`E×-` in the same commit as the result.

**A note on tier.** These are tier **F** — hypotheses. Several are one negative result away from being
closed, and closing one is a finding. Do not treat a null as a wasted ticket.

---

## 17a — The empty column

### R-020 · Isosurface extraction as a self-adjusting computation

**The gap, and it is the sharpest one found.** Across the ten-axis decomposition, axis 8 —
incrementality — reads **"full re-mesh" for every published algorithm.** An independent sweep of arXiv,
OpenAlex, CrossRef, Semantic Scholar and CORE found **incremental isosurface extraction does not
exist**. What exists is three adjacent things, none of which is it: time-varying *indexing* of a
precomputed series (T-BON), incremental construction toward convergence on a *fixed* field, and
dirty-region invalidation with no output-sensitivity bound — which is what every voxel engine and TSDF
system actually ships, at chunk granularity, with nothing proved about it.

**Axes where everyone differs are crowded. The axis where everyone is identical is where the space is.**

**The theory is in the corpus.** `10.48550/arXiv.2105.06712` (Acar, Anderson, Blelloch, Baweja) bounds
update cost by a **computation distance** between two executions and gives **work *and* span** — the
proof shape for edit-proportional re-meshing, already parallel, which matters because chunk meshing
already is.

**H:** re-meshing after a local field edit has computation distance `O(|edit|)` — i.e. the recorded
trace changes proportionally to the number of cells touched, not to grid size. **E1 is already
measured and supports this**: M-33/M-50 put a brush at **15–36%** of cells in its own bounding box,
reproduced live under a mouse.

**Harness:** instrument the extractor to record a dependency trace; apply a single-cell edit; measure
trace delta against edit size across grid sizes. **Records:** trace delta vs `|edit|` vs `n`.

**Falsified by:** trace delta scaling with `n` rather than `|edit|` — which would mean the extractor is
*unstable* in Acar's sense and edit-proportional repair is provably unavailable, closing the direction
with a reason. That is a real finding, not a failure.

**Worth if it holds:** the first isosurface extractor with a proved incremental bound. **Risk to
flag:** the closest prior art — Acar's own *Dynamic Well-Spaced Point Sets* (`10.1145/1810959.1811011`)
and *Kinetic Mesh Refinement in 2D* (`10.1145/1998196.1998254`) — is **unobtainable through the
pipeline** and does output-sensitive incremental *meshing*. Read before claiming novelty.

### R-021 · Maintain the contour tree, not the triangles

**The reframe from `10.48550/arXiv.1406.4005`:** the maintainable object is scalar-field level-set
topology, `O(log n)` per certificate failure, with certificates failing only on adjacent-vertex value
swaps or saddle collisions — and it handles general update operations, not just continuous motion.

**Caveat that must be carried:** that paper is **2-manifolds only** (`h: ℝ² → ℝ`, a triangulated
terrain). The two 3D results the question hinges on — Tarasov & Vyalyi 1998 and Safa & Wang 2014 — are
**both unobtainable**, and Edelsbrunner's 3-manifold Reeb maintenance is `O(n)` per certificate
failure, asymptotically no better than rebuilding.

**H:** the contour tree of a chunk can be maintained under a brush edit in time proportional to the
dirty set, where full recomputation is not. **Falsified by:** maintenance cost tracking chunk volume.

**Worth if it holds:** it changes axis 8 from "re-mesh fewer triangles" to "maintain topology,
re-derive geometry" — a different algorithm class.

### R-022 · Dynamic connectivity on the air sublevel set

**The cheap half of R-021, and it is buildable now.** The questions a game asks — *is this sealed? did
I break through? is this a chokepoint?* — are **not** all-thresholds queries. They are single-threshold
questions about connected components and bridges of the air sublevel set. That is dynamic
connectivity, which **is** measured and **is** cheap: microsecond queries, millions of updates/sec,
`O(log V)` depth when the spanning forest is unchanged — the common case, since most digging does not
alter connectivity. Bridges (chokepoints) are polylog amortized.

**The unoccupied part: dynamic connectivity has never been run on a voxel lattice.** Every measured
system was benchmarked on social/web graphs — Twitter 81K vertices, Stanford 280K. A voxel air-graph is
a **bounded-degree 6-connected lattice with 10⁶–10⁹ vertices.** Bounded degree should help; sheer V may
not. **And batching is untouched** — games edit thousands of voxels per explosion, not one;
`10.48550/arXiv.2002.05129` (batch-dynamic trees) is the right tool, is in the corpus, and has never
been pointed at this.

**H:** batched dynamic connectivity sustains a brush-sized edit under 1 ms at 128³.
**Falsified by:** per-edit cost scaling with lattice size rather than with the dirty set.

**Worth if it holds:** breakthrough-as-an-engine-event and sealed-volume-as-a-predicate, neither of
which any engine has. The benchmark alone is a contribution, since nobody has published it.

---

## 17b — One disputed table, replaced by one scalar

### R-023 · Persistence-thresholded ambiguity resolution

**arXiv returned literally zero** for persistence applied to Marching Cubes ambiguity. The nearest work
(Kissi & Tierny 2024; Brüel-Gabrielsson 2018) simplifies **the field globally**, then contours. This
hypothesis is a **local decision inside the cell**, and nobody has tried it.

**It is live right now.** A-002i and A-020b are blocked on architecture for exactly the singular and
tunnel cases, and A-002b's own reasoning is that *there is no correct published table to transcribe* —
Custodio et al. proved Chernyaev's interior test tracks a quadratic where the true saddle trajectory is
hyperbolic, and Lewiner's reference implementation omits cases 10 and 12 entirely. **The guaranteed
version is 730 subcases.**

**The move:** stop asking *"is there a tunnel"* and ask *"does this tunnel have persistence above ε."*
Below threshold, mesh closed; above, mesh open. A disputed 730-entry table becomes a **computable
scalar with a stability theorem behind it** — and a knob a game wants.

**Scope it honestly.** Persistence is defined for a filtration of a *space*; the ambiguous cases are
ambiguous precisely because eight corner samples underdetermine the trilinear interpolant, so a
tunnel's persistence depends on which interpolant you assume. **This does not remove the modelling
choice — it relocates it**, from a hardcoded table into one tunable, stability-backed threshold. Claim
exactly that.

**H:** a persistence threshold reproduces MC33's topology on the fields where MC33 is agreed correct,
and differs only on cells where the published algorithms disagree with each other.
**Falsified by:** disagreeing on cells where Chernyaev and Lewiner agree.

**Worth if it holds:** it retires A-002b, A-002i and A-020b together, and it is the one idea here
nobody has published even a negative result on.

---

## 17c — Two things sitting on this crate's own seam

### R-024 · Does field-sealed imply mesh-sealed?

**Nobody has established this, and every paper treats the two as interchangeable.** A cell can be
topologically connected in the field and still produce a watertight surface, or the reverse, depending
on the case table and the interpolant.

**This is a day of work and it is publishable alone:** extract a mesh; compute connected components of
the air sublevel set; compute connected components of the mesh complement; **assert they agree.**

**H:** they agree for Marching Cubes on all eight reference fields, and **disagree for at least one
dual method** — the duals place vertices by solve rather than on the interpolant, which is exactly
where the correspondence could break.

**Falsified by:** universal agreement — a stronger correctness statement than the crate currently
makes, and worth saying so.

**Worth either way:** it is the precondition for every mechanic in R-022, and the gap sits exactly on
the seam this crate occupies.

### R-025 · Second-order vertex placement

**Both ingredients are published separately and have never been composed.** Jet/Hessian fitting
(Cazals & Pouget; Jiao & Zha, in corpus) and QEF placement (Ju et al.) — nobody fits the **second
fundamental form** per cell instead of tangent planes.

**It attacks a measured term.** P-2's error model is `O(|e|²κ)`, and on a true SDF the Hessian's
nonzero eigenvalues at a surface point are `−κᵢ/(1−κᵢd)` — **principal curvatures fall out of samples
already taken**, no medial axis involved.

**H:** curvature-aware placement improves Hausdorff on smooth fields (`sphere`, `torus`, `gyroid`) by
>20% over planar QEF at fixed resolution, and **does not** improve it on `box_exact` — where the
surface is flat and the second-order term is zero.

**Falsified by:** no improvement on smooth fields, meaning curvature estimation is too noisy at game
resolutions — which is Aamari & Levrard's minimax bound biting, and worth recording as such.

---

## 17d — The result you already have and have not claimed

### R-026 · Write up the head-to-head

**M-001 produced the comparison that does not exist in the literature**, and M-004's writeup ticket is
archived while the paper is not written. Verified: **no paper since 2020 benchmarks Marching Cubes vs
Surface Nets vs Dual Contouring against each other**, and Surface Nets — the thing engines actually
ship — **has no credible published timings at all.**

You additionally hold results that **contradict** published figures: M-51 and M-55 falsify the
literature's `2–3×` Marching Tetrahedra ratio (measured `~3×` triangles for `4.3%` worse geometry, not
86%), M-1's `V_sn = V_mc + χ` identity, M-53's four-corner table of manifold × intersection-free, and
M-54's `101×` Dual Contouring accuracy advantage on sharp fields.

**This is the least speculative item in the phase and the only one whose result is already in hand.**
The remaining work is Open SciVis volumes for comparability (H-005), mesh-quality metrics for the
table reviewers expect (H-003), and prose.

---

## Ordering

| | Why |
|---|---|
| **R-024** | One day, publishable alone, and it gates R-022 |
| **R-026** | The result exists; only the writing is missing |
| **R-022** | Buildable now on measured foundations; the benchmark itself is unpublished |
| **R-023** | Retires three blocked tickets at once, and nobody has even a negative result |
| **R-020** | The biggest space, and the one most at risk from unread prior art — get Acar's two papers first |
| **R-025** | Cleanest hypothesis, most likely to null out honestly |
| **R-021** | Highest ceiling, worst evidence base — two load-bearing papers unobtainable |

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
