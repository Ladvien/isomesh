---
description: One iteration of the isomesh improvement loop (SDF to mesh, caves, overhangs, carving)
argument-hint: "[optional focus, e.g. 'issue #10 only']"
---

One run of this command is one iteration of a long-running loop on isomesh, an engine-agnostic isosurface workspace: signed distance field in, triangles out. `crates/isomesh` is the `no_std` core (Marching Cubes, Dual Contouring, Surface Nets, Transvoxel seams, SDF construction and CSG, chunked edits; `libm` its only dependency), `crates/isomesh-gpu` extracts on wgpu, and `bevy_isomesh` (its own workspace, with every example) turns a chunked `VoxelVolume` into Bevy meshes. Games use it for what a heightfield cannot hold: caves, arches, overhangs and terrain you carve. Sanctuary's destructible wall is 20 `bevy_isomesh` volumes of 856 chunks; bevy_emerge's family issue F5 (#101, "Caves, overhangs and carving") builds its Volume tool on it; bevy_elements' voxel snow reads its volumes. Do ONE improvement per iteration, the single most valuable one from a game developer's point of view, finish it completely, then stop.

This crate must be correct AND deterministic AND cheap on the caller's frame. A crack at a chunk seam, a flipped triangle, a non-manifold edge on a closed field, a mesh that changes between runs or operating systems, and a carve that stalls the frame it lands on are all failures, however good the render looks. The judgement here is a geometry-processing engineer's who also ships games: a mesh that passes the validity gate its field calls for, golden hashes that hold bit for bit on macOS and Linux, a remesh that costs only the chunks an edit touched and allocates nothing in the hot path, and a core API that stays free of any engine or math library (`CLAUDE.md`, "Hard rules").

You do all the work yourself, including looking at any images and judging them. Do not delegate
visual review.

Extra instructions for this run (a focus such as "issue #10 only"): $ARGUMENTS
They override the choice of item in section 4, not the gates.

State lives in `docs/loop-log.md`; rules and traps live in `CLAUDE.md`. Follow them; do not
re-derive them.

## Your work is saved even if you are killed — read this before you start

**The driver saves whatever you leave uncommitted, without touching your branch.** Every ten minutes
it records the whole tree to a snapshot ref (`refs/wiggum/pass-N-...`), and again when the pass ends.
Your `git status`, your index and your branch are unchanged — so this is not a substitute for your own
commit step, which is still the only thing that puts work on `main`.
Recover a killed pass's tree with `git -C <repo> log --all --oneline --grep='WIP(pass'` then
`git show <ref>`.

1. **Commit as you go, not once at the end.** Your commit step is last, behind implement, gate and
   benchmark, so a pass killed at the ceiling never reaches it. Commit the research note alone before
   writing code, as your sections already say; then commit the implementation, then the docs and the
   log. Each commit that lands is work a killed pass cannot lose.
2. **A snapshot ref is not a commit on your branch.** `git log` will not show it. Do not look for
   one; if an earlier pass was killed, its tree is in `refs/wiggum/` and nobody's gate ran on it.
3. **Parked work is a branch, and you look for it before you pick an item.** When a pass is stopped or
   killed with its tree uncommitted, the driver's next start commits that tree to a branch
   (`wiggum/parked-pass-N-<time>`, or on the pass's own `loop/iter-N`) in a commit titled
   `WIP(pass N, parked)` and puts the clone back on `main`: you start clean and nobody tells you. So run
   `git branch --list 'wiggum/parked-*'` and `git log --all --oneline -n 5 --grep='parked'` once, after
   you have triaged the filed issues and before you commit to an item. If the newest parked tree is for
   the item you would pick, or for an issue or bug that outranks your pick, **adopt it**: it is a first
   draft that was never built, tested or gated and was cut off mid-work, so read it, bring its files into
   your own `loop/iter-N` branch (`git checkout <parked-branch> -- <paths>` or `git cherry-pick`),
   finish it, and run the real gate on the result before it ships under your name. Never merge the
   parked branch itself and never delete it. If it is stale (its files conflict with `main`, or the item
   was fixed meanwhile), say so in one line in the log and ignore it. A parked branch that exists only on
   this clone holds a lot of data: the driver keeps anything over 10 MB of new data off origin.


## 1. Orient

1. `CLAUDE.md`: the rules, the six hard rules, the version pins, the verification the validity gate
   needs, the commands. Then `FINDINGS.md`'s index (the epistemic state: look a claim up there before
   acting on it, and cite its tier) and the top of `BACKLOG.md` (the ticket queue and how it is worked).
2. `docs/loop-log.md` (create it if missing): what earlier iterations did, found and deferred, and the
   newest iteration's ranked backlog. The next iteration number is the last one plus one.
3. `git log` since the last logged iteration.
4. The `docs/research/<system>.md` note of any system you may touch, and `docs/research/corpus.md` for
   the reading list per domain (create both the first time a pass researches).
5. `git status --short`. The driver (`wiggum.sh`) starts a pass only on `main` with no
   tracked changes. If you find tracked changes you did not make, stop and log them; do not commit or
   revert them.
6. **Re-scan the corpus for new sources.** Run `catalog_recent` (limit 200) and add any isosurface
   extraction (Marching Cubes, Dual Contouring, Surface Nets, Transvoxel, adaptive or octree meshing),
   signed-distance-field and CSG, voxel or destructible terrain, mesh validity (manifoldness, orientation,
   self-intersection), level-of-detail seam, GPU meshing or chunked-streaming document that
   `docs/research/corpus.md` does not list, as its own commit (`Research: corpus`).

Preflight home-still now, so research is not blocked later:

- `system_status` (MCP) or `hs status`. HomeStill's servers run on other machines; never restart them
  from here. If distill is down, skip `distill_search` and `distill_index`, search with the CLI
  (`hs distill search "<query>" --output json`) or `paper_search`, and list each unindexed note under
  `Pending index`. Record the live state in the log; do not copy an earlier iteration's outage note.
- Scribe: check `scribe_health`. `scribe_convert` may queue or fail; then rely on the corpus.
- **An outage is not a reason to idle.** `paper_search`, the CLI and the papers already quoted in this
  repository need no GPU; research from them and build.

## 2. Capture the current state

`CARGO_TARGET_DIR` is set by the supervisor and shared by all three workspaces (the root,
`bevy_isomesh/` and `isomesh_web/`), so every binary is under `$CARGO_TARGET_DIR/release/`, never
`./target`. Build first: `cargo build --release --workspace --all-targets` from the root, then from
`bevy_isomesh/` only the examples this pass will run (`cargo build --release --example <name>`; one
optimised Bevy example is about 2 GB of output, so never build all 49 at once). Keep everything under
`/tmp/iN/` (`N` the iteration), the baseline under `/tmp/iN/before/`. Keep a worktree of the starting
commit for later comparisons, **detached, at a commit, and off `/tmp`**: `git worktree add --detach
~/iN_base <sha>`, built with `CARGO_TARGET_DIR=~/iN_base-target`.

Section 9 says how to name, build and remove that baseline. A baseline left behind costs the next pass
a build cache's worth of disk, and one with `main` checked out stops anything else checking `main` out.

**The numbers are the state.** This crate's cost lands on its caller's frame and its correctness on the
caller's collider and renderer, so measure both the way a caller pays them:

- **Extraction cost:** `cargo bench --bench extract` (criterion, the re-mesh path per algorithm with a
  reused buffer). Record each case's median from criterion's `time: [low median high]` line.
- **Structure and quality:** `cargo bench --bench shootout` writes `docs/measurements/shootout.csv`
  (five extractors on identical fields: vertices, triangles, non-manifold edges, Hausdorff,
  self-intersections per 1,000 triangles, `median_ms`). Then `./scripts/regress.sh shootout`: the
  structural columns compare exactly, the measured ones within the script's stated tolerances. Its
  baselines are per machine (`docs/measurements/baseline/<name>-<machine>.csv`, `scripts/machine.sh
  --slug`); when this host has none, run `./scripts/regress.sh --accept shootout` on the baseline
  worktree and copy the file into `/tmp/iN/before/`, never into the repository unless a pass commits it on
  purpose as its own `assumed owner exception: <machine> baseline` commit.
- **Edit cost:** for a change to edits, chunks or the plugin, `cargo bench --bench edit_trace` (what a
  local field edit costs, chunk by chunk) and `cargo bench --bench stage_breakdown` (where a usable mesh's
  time goes). On the Bevy side, `game_budget` (`ISOMESH_BUDGET_US`) is the frame-cost view of the remesh
  queue: run it in release and quote its logged per-frame figure.
- **Determinism:** `cargo test -p isomesh golden` (T-007's committed hashes, generated on one OS and
  verified on the other). A change that moves a hash on purpose says which and why.
- **GPU:** this host is macOS on Apple silicon, so Metal is real here: `cargo test -p isomesh-gpu`, and
  `cargo run --release -p isomesh-gpu --example gpu_vs_cpu` for a change on that path.

Every bench rewrites its committed CSV under `docs/measurements/`. Copy the output to `/tmp/iN/` and
restore the committed file (`git restore docs/measurements/<name>.csv`, your own write) unless the pass
updates that measurement on purpose, in which case `./scripts/csv_provenance.sh` must pass on it. Record
every figure in `/tmp/iN/before/bench.txt`.

The loop's host shares its cores with other loops, so a single run is noise: run the bench three times
before and three times after, back to back on the same machine, and compare medians of medians.

The visible behaviour is the `bevy_isomesh` examples (`quickstart`, `marching_cubes_sphere`,
`game_dig`, `game_carve_seams`, `sealed_cave`, `transvoxel_seams`, `game_terrain_stream`,
`game_showcase` and the rest; `DEMOS.md` lists them). When a change can be seen there, capture a still
with `ISOMESH_SCREENSHOT` or frames with `ISOMESH_CAPTURE` / `ISOMESH_CAPTURE_FRAMES` (`ISOMESH_WINDOW`
sets the size), or a GIF with `./scripts/record_gif.sh <example> /tmp/iN/<name>.gif`. Run examples with
`--release` only. The host renders natively on Metal; a window opens on its display.

## 3. Inspect closely (you judge the results yourself)

Read every number and every image yourself. Write the **mesh checklist** to `/tmp/iN/findings.md`,
each line with a measurement, not an impression:

- **validity, by field:** for each reference field the change touches, the gate its `ReferenceField`
  metadata calls for (`is_closed()` for closed fields, `is_manifold()` for open ones, χ only where
  `expected_euler()` is `Some`), with `boundary_edges` and `non_manifold_edges` quoted as counts;
- **orientation:** interior edges traversed in the same direction by both faces: 0 (the only check that
  sees an inside-out triangle);
- **non-manifold vertices** by the link walk, and **index bounds**: no index at or past the vertex
  count, no triangle with two equal indices; degenerate (near-zero-area) triangles quoted as a count;
- **seams:** a carve or field edit that crosses chunk boundaries meshes with zero cracks
  (`ChunkSeams` / `is_seam_free`), and Transvoxel LOD seams likewise when LOD is in play;
- **no NaN or infinity** in any position or normal, and no panic on a caller's field, grid, edit or
  chunk index (refused with an error, not panicked on), including empty chunks and edits outside the
  volume;
- **determinism:** the golden hashes unchanged (or each moved one named), and the same input twice
  gives byte-identical buffers;
- **cost:** extract medians, shootout `median_ms`, and for an edit change the edit-trace and
  game_budget figures; no allocation after warm-up on the re-mesh path (hard rule 6);
- **footprint:** `cargo tree -p isomesh -e normal` resolves to exactly `isomesh` and `libm`;
- **the examples:** a still of the changed example, read at full frame and zoomed into seams, silhouettes
  and cave mouths (cracks, slivers, flipped shading, holes).

## Owner input: recommend, record, proceed

**This loop runs fully automated (owner, 2026-10-29): never stop to wait for the owner.** When a pass needs the owner's input (a design choice the issue or item does not settle, such as a colour, a coverage, an API shape, a schema change or which of two behaviours; or a protected file that must change), it answers the question itself:

1. **Recommend.** Name the two or three real options and pick one, with one line of why. Prefer the conservative choice: reversible, conventional, consistent with the existing code and the issue's stated intent, and the smallest that meets the acceptance.
2. **Assume the owner chose it, and build it this pass.** A decision is never a reason to skip an item or end a pass.
3. **Record it where the owner reviews.** On the issue, a comment beginning `Assumed decision (owner may override):` with the options, the pick and why, and the label `assumed-decision` (`gh label create assumed-decision --color FBCA04 2>/dev/null; gh issue edit <N> --add-label assumed-decision`). With no issue, in the commit body and the log line, beginning `Assumed decision:`.
4. **The owner's reply wins.** An owner comment that contradicts an assumed decision overrides it: the next pass that reads it reworks the item to the owner's choice before anything else.
5. **Earlier blocks are lifted.** An issue or backlog item logged as blocked on, or waiting for, an owner decision is workable now: decide it this way.
6. **Protected checks change only additively.** An assumed decision may add to a protected test or check (a new row or entry, a new field's empty default in a literal, a count or version raised by exactly what the change adds), with the commit message `assumed owner exception: <what>`. It never removes, loosens, skips or re-pins an existing assertion, golden, eval answer or threshold; that is the one block left, so say so on the issue and go on.
7. **Never, whatever you assume.** These need the owner's word in a conversation, and no prompt can give it: yanking a release or changing crate ownership; force-pushing, rewriting history or deleting branches; deleting files or work this pass did not create; editing this command file; spending money or sharing credentials.

## Filed issues come before new features

**A filed issue is a consumer telling you what broke, which outranks anything this pass invented.**
Before you pick an item, list the repository's open issues and triage them:

```sh
gh issue list --repo Ladvien/isomesh --state open --limit 100 \
  --json number,title,labels,createdAt --jq '.[] | "#\(.number) \(.title)"'
```

(If `gh` is absent, the same list is
`https://api.github.com/repos/Ladvien/isomesh/issues?state=open&per_page=100`.)

Then, in this order:

1. **A bug issue is the pass's item, ahead of any enhancement.** Anything reporting wrong behaviour —
   a crash, a panic, a visual artifact, a value that is implausible, a documented promise that is not
   kept — is a bug. Take the oldest such issue that is still reproducible, and say why it beats the
   others.
2. **An issue labelled `consumer` is the pass's item next, ahead of the series and your own candidates.**
   A project built on this crate needs it and is waiting. The `consumer:<name>` label says which project;
   `consumer:bevy_emerge` is the owner's editor. Its own work resumes when the issue is closed as completed
   and it moves its git pin to a `main` commit containing the fix. Take the oldest one you can work this pass.
   - **Build the general form.** The requester is the example, never a name in this crate's API, a path or a
     default: this crate serves other projects too. Where an issue asks for something shaped to one project,
     build the general form and say so on the issue.
   - **Prove it here**, with a test in this repository that fails before and passes after against the issue's
     acceptance. Never build or run the requesting project to prove it.
   - **Pushing to `main` is the delivery.** Requesters pin a git revision, so a consumer fix never needs a
     publish or a version bump.
   - **One that needs a decision the issue does not settle is decided, not skipped:** recommend, record and
     build it (see "Owner input: recommend, record, proceed"). Only a change that would weaken a protected
     check cannot be worked: say exactly what is needed on the issue (`gh issue comment`), write one log
     line, and go on. A blocked consumer issue never stops the pass.
   - **Read `git log` first.** While this loop is stopped, the requester's loop may fix the issue here at the
     source, in its own commits with tests; check whether that has happened so the work is not done twice.
3. **Every other open enhancement issue is the pass's item next, ahead of your own candidates** (the
   owner's rule, 2026-10-29): a filed request was asked for; a candidate this pass invents was not. The
   order of work is: **bugs first** (a bug issue, then a bug already recorded in the backlog or the log),
   then `consumer` issues, then every other open enhancement issue, oldest first. Only when each of them
   is done, blocked or a release issue does "Pick the single most valuable item" choose a feature of your
   own. A blocked issue gets a `gh issue comment` saying exactly what it needs and is skipped; it never
   stops the pass. Name the issues you skipped, and why, in the log line, so the next pass sees the list
   was worked through.
4. **Close the loop.** An issue your pass fixes is referenced in the commit body (`Fixes #N`) **and**
   in the report, so the next pass can see it was addressed. Do not close an issue you did not fix,
   and do not close one whose fix this pass did not verify. For a `consumer` issue the close is the signal
   the requester waits for: close it only through a verified fix pushed to `main`, never as not planned
   without the owner's word.
5. **A release issue is done by the next release.** An issue asking to publish or bump a version is closed by the pass whose release (the "Release" block in the ship section) publishes it, with the crates.io URL.

**When the API fallback answers `403` or `429`** (rate limit — 60 requests/hour per IP unauthenticated,
shared by everything on this network), do not treat it as "no issues": write one line saying the list
could not be read and continue. Never spend a second call to retry it in the same pass, and never infer
from silence that nothing is filed.

**An empty issue list costs nothing:** if it is empty or `gh` is unreachable, write one line saying so
and continue — never spend a pass on the lookup.

## 4. Pick the single most valuable item

**You are here only when no bug and no open issue can be worked this pass** (see "Filed issues come
before new features"); otherwise that issue is the item. Rank candidates from your findings, the loop log's backlog, the topmost unblocked tickets of `BACKLOG.md` and the extra instructions, always asking
what matters most to a game developer who will ship this:

1. A wrong mesh: a crack at a chunk or LOD seam, an inside-out or flipped triangle, a non-manifold edge
   or vertex on a field that should have none, a NaN, a panic on caller data, or output that changes
   between runs or operating systems.
2. A cost the caller cannot afford: a carve whose remesh blows the frame budget, a remesh of chunks the
   edit never touched, an allocation per extract, a stall on the main thread while finished meshes are
   applied, a render-world cost for a chunk that draws nothing.
3. API friction for a game carving caves: setup that takes more lines than it should, an edit history a
   game cannot store, replay or undo, a query it needs and cannot make (inside test, surface height, ray
   hit, seam check, changed chunks), a value it must recompute because the crate hides it.
4. Missing features a game with caves and overhangs expects: smooth CSG brushes (sphere, box, capsule)
   with bounded edit regions, crack-free LOD at distance, per-material or painted volumes, collider
   meshes, streaming without hitches.
5. The examples, `DEMOS.md` and the site's demos looking unconvincing next to what the extractors can show.

Ownership is a gate: `crates/` never depends on Bevy or a math library and `crates/isomesh` adds no
dependency without a written justification (`CLAUDE.md`, hard rules 1-3); everything that touches Bevy
lives in `bevy_isomesh`. Physics and colliders' simulation are avian's, snow and wetness on volumes are
bevy_elements', the editor's tools, undo stack and panels are bevy_emerge's, and a game's world layout
and when it remeshes are the game's. Do not build those here; add a hook instead.

State in one sentence why your pick beats the runner-up. If it is too big, carve off the most valuable
self-contained slice. Record the ranked candidate list in the log, not only the winner.

**A pass never ends idle** (owner, 2026-10-29: "If they're not addressing bugs or enhancing features, they need to evolve their own code base to support Bevy Emerge better or add new features. Research these things in home still, implement them. They shouldn't stay idle."). An empty issue list or a stale backlog is this step's input, never a reason to ship nothing, and a candidate is not "speculation" once it is researched. When the ranked list above gives no workable item, take the first of these that applies:
1. **What the crate's consumers hand-roll.** Read how Sanctuary uses this crate (`gh repo clone Ladvien/sanctuary /tmp/iN/sanctuary -- --depth 1`, then `rg -l isomesh`; its `isowall` module drives the destructible wall; removed afterwards) and bevy_emerge's family issues (`gh issue list -R Ladvien/bevy_emerge --label family`, above all F5, `gh issue view 101 -R Ladvien/bevy_emerge`: a Volume tool storing its volume as an edit list of SDF primitives and operations, carving without cracks at chunk seams, one undo step per carve). Find a formula, a wrapper or a query either writes by hand around this crate, and build its general form here, never named for either. Prove it with a test here.
2. **A feature a shipping game in this domain expects and this crate lacks**, found by research: home-still first (`distill_search`; for a paper the corpus lacks, `paper_search`, `paper_download`, `scribe_convert`, then `distill_index`), then shipped-game practice, then section 5 before any code.
3. **Depth on the weakest existing feature**: measure it against its research note's validation targets and improve the worst gap.

Record the ranked candidates and the pick in the log. A pass that ships nothing names the concrete blocker that stopped all three; "nothing new to do" is not a blocker. **"Not researched yet" is not a blocker either: it is this pass's work.** Take the top item, do section 5 for it now (home-still preflight included), commit the note, and build its first slice. "Building it would be speculation" means the research is missing, so do the research. A log-only iteration is a failed pass, not an outcome, and a pass never copies an earlier iteration's reasons for shipping nothing (measured 2026-10-29: iterations 194 and 195 each shipped a log line saying every item was "not researched", and neither ran the research).

## 5. Research thoroughly (before any code)

Every feature gets researched, new or improved, with no exceptions. A mesher is published geometry
put to work under a game's frame budget, so every iteration cites **three** kinds of source:

1. an **algorithm** source: the extraction or field method itself (Lorensen and Cline's Marching Cubes
   and its ambiguity resolutions, Ju et al.'s Dual Contouring, Surface Nets, Lengyel's Transvoxel, SDF
   construction and CSG); `docs/research/2026-08-10-meshing-algorithm-catalog-v2.md` and the adjacent-math
   audit are the starting points `CLAUDE.md` names;
2. a **shipping-game** source: voxel, cave and destructible terrain in shipped games and engines (GDC
   talks, Game AI Pro / GPU Gems chapters, studio and engine write-ups such as C4's Transvoxel, Astroneer,
   Deep Rock Galactic, Teardown);
3. a **technical** source: cache-friendly chunk layout and edit regions, SIMD, GPU compute meshing,
   frame-budgeted task scheduling, allocation-free buffers, robust predicates and floating-point
   determinism.

When the home-still MCP bridge answers `401 Refresh token expired`, search with the CLI,
`hs distill search "<query>" --output json`, which refreshes its own token.

If the corpus lacks one of the three, `paper_search` for it, and record the gap under the note's
`## Pending index`. Write the note with the headings this loop requires (`## Sources`,
`## Craft principle`, `## How shipping games do it`, `## Model adopted`, `## Parameters`,
`## Validation targets`, `## Deviations`, `## Pending index`) plus `## What great looks like`: the
visual qualities and reference images you are aiming for. `## Craft principle` (the model source) and
`## How shipping games do it` must each **quote** the source they rest on, not paraphrase it from
memory.

Any changed default goes under `Deviations` with a citation, or is marked a design default. Commit the
note alone before writing code: `git add docs/research/<system>.md && git commit -m "Research: <system>"`.

## 6. Implement within the crate's rules

- Keep every rule in `CLAUDE.md`. In particular the six hard rules: no math library in a public
  signature (`[f32; 3]`, `[u32; 3]`); no Bevy in `crates/` (manifests, code, resolved graph); nothing in
  `crates/isomesh` `[dependencies]` but `libm` without a justification in `CLAUDE.md` in the same commit;
  no performance number without the committed bench that produced it; no guessed algorithm detail (a
  case table, a sign convention or an edge order comes from a cited source); output buffers are
  caller-provided and reusable.
- No `unwrap()`, no `expect` on caller data, no panicking index: fields, grids, edits and chunk indices
  come from callers. Leave the paper a piece of code rests on in a comment beside it.
- **Determinism is checked, not argued:** transcendentals go through `libm`, never `std`; no `HashMap`
  order leaks into vertex order. A change that could move a single output bit proves it did not with the
  golden hashes, or names each hash it moved on purpose.
- An extraction change passes the validity gate its field calls for (never `if field == gyroid` in a
  test), plus orientation, the link-walk vertex check, index bounds, and a recorded self-intersection
  count. The first shader gets a naga validation test.
- Version pins: `bevy` 0.20 with exactly the `wgpu` and `glam` it pins, moved together in one commit,
  never alone; no `"x11"` feature (this is macOS / Metal).
- `bevy_isomesh` and `isomesh_web` are their own workspaces: a change in `crates/` is gated in all
  three. Every example is a Bevy one in `bevy_isomesh/examples/`, run with `--release`.
- `FINDINGS.md` gets an entry, with its tier, in the same commit when a pass measures something,
  contradicts something written down or earns a rule; `./scripts/findings_index.sh` regenerates its
  index. A `BACKLOG.md` ticket a pass finishes is checked and moved to `BACKLOG_ARCHIVE.md` in the same
  commit, with its annotation; `./scripts/backlog_gate.sh` holds the counts.
- `CLAUDE.md` asks to stop and ask on an ambiguous design choice; in a loop pass that question is
  answered by "Owner input: recommend, record, proceed" above (owner, 2026-10-29). Its force-push and
  deletion rules hold unchanged.
- Plain functions with their own tests for any new maths, testable without an app; tests are inline
  `#[cfg(test)] mod tests;` siblings, as the crate has them.
- Re-export new public items from the crate root; `missing_docs` holds, and a sign convention, a
  coordinate order or a winding order is stated in the doc comment.
- One path per feature, no fallbacks or shims; delete the code the change replaces.
- Tests: for a bug, a test that fails before and passes after. Permanent tests check what a caller sees
  (vertices, indices, validity counts, hashes, refusals); throwaway checks go in `/tmp`.

## 7. Verify and sign off: it must be RIGHT and FAST

You are grading your own work, so grade it as a hostile reviewer would — a geometry-processing engineer
who will run the validity gate, the seam check and the golden hashes before reading anything you are
proud of.

- Rebuild, re-run the bench three times and the mesh checklist under `/tmp/iN/after/`, and compare
  every line with the "before".
- Write the critique BEFORE the verdict (`/tmp/iN/critique.md`): every flaw you can find, against the
  research note's `## What great looks like` and `## Validation targets`. Look for what is wrong, not for
  confirmation that your change worked.
- Then give each feature exactly one verdict: GREAT, NOT YET or FAIL. GREAT means a game could ship on it
  today. "Fine", "acceptable" and "good" are NOT YET.
- On NOT YET or FAIL, write the specific changes that would move it to GREAT, apply them, rebuild,
  re-measure and re-judge. Do not lower the bar to pass.
- After 5 failed rounds on one feature, stop: log the blocker, the remaining flaws and your best
  hypothesis, and revert to the last GREAT state.
- **Performance gate**, equal in weight: median of the three `extract` bench medians before and after,
  same machine; no regression over 5% on any case without a logged, justified trade. `regress.sh
  shootout` passes against this machine's baseline (structural columns exactly). A change aimed at edit
  or frame cost shows `edit_trace` and the `game_budget` frame figure before and after.
- Run the full gate (it is CI's, across all three workspaces; `preflight.sh --full` is the repository's
  own copy of most of it, and the rest is what CI adds). One-time setup on a host without it:
  `python3 -m venv ~/.venvs/isomesh && ~/.venvs/isomesh/bin/pip install 'markdown>=3.7,<4'
  'pymdown-extensions>=10,<11'`.
  1. `PYTHON=~/.venvs/isomesh/bin/python ./scripts/preflight.sh --full` (fmt in all three workspaces,
     backlog gate, findings index, doc facts, version gate, csv provenance, readme sync, toolchain drift,
     `isomesh_web` clippy and test, the site's link check, `cargo clippy --workspace --all-targets -- -D
     warnings`, `bevy_isomesh` check, clippy and rustdoc, root rustdoc, rule 2, the two-package
     dependency count, `cargo test -p isomesh`, `cargo test -p isomesh-gpu`, `bevy_isomesh` lib and doc
     tests, MSRV 1.89)
  2. `cargo test --workspace --lib --tests && cargo test --workspace --doc && cargo test -p isomesh --features experimental --lib --tests`
  3. `./scripts/regress.sh --self-test && [ "$(cargo metadata --format-version 1 | grep -c '"name":"bevy')" = 0 ]`
  4. `cd bevy_isomesh && cargo +<its rust-version> check --all-targets` (the `rust-version` in
     `bevy_isomesh/Cargo.toml`, 1.97.1 at setup)
  5. `./scripts/publish.sh --dry-run`
  6. `PYTHON=~/.venvs/isomesh/bin/python ./scripts/build_web.sh` (the wasm32 demos and the site, as CI's
     `site` job; needs the `wasm32-unknown-unknown` target and the `wasm-bindgen-cli` version pinned in
     `bevy_isomesh/Cargo.lock`)

## 8. Ship

Start only after section 7 passes. Then release (the "Release" block at the end of this section).

1. Document: rustdoc for every changed public item, the README for user-visible changes (the root
   `README.md`, `crates/*/README.md`, `bevy_isomesh/README.md` and `DEMOS.md`, as each applies), the
   research note's `Deviations` and `Pending index`, and an entry under `## [Unreleased]` in
   `CHANGELOG.md` (`Added`, `Changed`, `Fixed`, `Removed`).
2. Regenerate the GIF or still under `docs/gifs/` or `docs/screenshots/` that shows the changed
   behaviour when it has moved, with `./scripts/record_gif.sh` or the example's `ISOMESH_CAPTURE` knob.
3. Log: append `## Iteration N — <date> — <title>` to `docs/loop-log.md`:
   - what changed and why it ranked first (and the runner-up)
   - the sources used: the algorithm, shipping-game and technical stems, by stem
   - your critique and verdict per capture, with evidence filenames
   - the mesh checklist's numbers, before and after
   - the bench's medians, before and after
   - the commit hash, and any sibling-repo commits
   - the top 5 remaining items, ranked, for the next iteration
4. Branch and commit: `git checkout -b loop/iter-N` from `main`. Commit by path with messages that
   state the user-visible change: the research note (already committed in section 5), the
   implementation, the docs and changelog, then the log (`Loop log: iteration N, <title>`).
   `git status --short` must be empty afterwards.
5. Push and merge: `git push -u origin loop/iter-N`, then `git checkout main && git pull --ff-only &&
   git merge --no-ff loop/iter-N -m "Merge loop/iter-N: <title>"`. On a conflict, `git merge --abort`,
   log it, and stop.
6. Re-run the gate on the merged tree; a failure stops the pass before the push. Then
   `git push origin main`.

**Release, every pass that changed the package** (owner, 2026-10-29: loops publish, and never wait on CI to do it). After the push, from that pushed commit, on a clean tree with the gate green: if a file `cargo package --list` ships changed since the last published version (`git diff --quiet v<last> HEAD -- <those files>` fails; with no tag `v<last>`, diff from the commit that set that version in `Cargo.toml`, found with `git log -S'version = "<last>"' -- Cargo.toml`), release:
1. Decide the version with `cargo semver-checks check-release` (`cargo install cargo-semver-checks --locked` if missing). Pre-1.0: a breaking change bumps MINOR and resets PATCH, anything else bumps PATCH; never 1.0.0. Move `## [Unreleased]` in `CHANGELOG.md` to `## [<version>] - <date>`, bump `Cargo.toml` and `Cargo.lock`, and commit `release: v<version>`.
2. `cargo publish --dry-run`, then `cargo publish`, locally, with the credentials `cargo login` stored; never print them. Never wait for, watch or re-run a GitHub workflow to release.
3. `git tag -a v<version> -m v<version>`, then push the commit and the tag (merge, never force, on a refused push).
4. `cargo info <crate>@<version>` resolves; put the crates.io URL in the log entry.

**At most one release every 6 hours** (owner, 2026-10-30: bevy_lights published v0.2.1 through v0.2.20 in one day, and then crates.io answered 429 for three passes in a row, each of which did nothing but retry). Release only when the last tag is more than 6 hours old (`git log -1 --format=%ct <last tag>`); until then the changes ride in the next release. A publish refused with 429 is not retried in this pass or the next: log it once and build, and the next release carries it. A pass that only retries a release is an idle pass. Any other failed publish is logged with its error and retried at the next release. Never yank, never reuse a version; fix a bad release forward.

**This workspace releases three crates under one version** (`isomesh`, `isomesh-gpu`, `bevy_isomesh`). The version lives at the six manifest sites across the three workspaces that `./scripts/version_gate.sh` checks; bump them all and the three lockfiles, run the gate, and publish with `./scripts/publish.sh --dry-run` then `./scripts/publish.sh`, which uploads in dependency order only what crates.io lacks. "Changed the package" means any of the three. CI's `publish` job runs the same script on a version-bump push and finds nothing left to upload; never wait for it. With no tag yet, the last published version is 0.0.11, set at commit `7016263`.

### README, demos and docs stay true

A README is read before anything else, so a false line there costs more than a false line anywhere else.
An audit on 2026-10-27 found one in each of the four loop repos: a `tour.gif` captioned "at noon" over a HUD
reading 10.00 h with a running clock; a `ties.gif` captioned for rivalry colours that appear in none of its 240
frames; a hero caption promising an effect its file never shows; and a demo video captioned as the full
recording of the GIF above it when it shows different objects. Before you push a change a user can see:

1. **Put it in the section it belongs to.** Do not append to the introduction or to the end, and do not
   restructure the README in a feature pass. It reads in this order — what it is in one line, the demo, what
   you get, install, quick start, then the detail sections — and an item added in the wrong place undoes it.
2. **A caption describes the file as it exists.** Open the GIF or still (frame count, size, any HUD burned into
   it) and write only what you can read off it; a figure quoted from program output comes from running that
   program now, not from the log. Recipe: the example's capture knob (section 2).
3. **Media you cannot regenerate here is logged, not left lying.** With no GPU or display on this machine, name
   the stale file and its path in the report and change its caption to what it still shows.
4. **The README's code compiles.** Edit a README snippet only together with the code it shows, and check it against the crate-root doctest it mirrors.

## 9. Clean up when you are done

Worktrees and build directories outside the clone are not free: a `target/` for this crate is tens of GB, a
loop host's `/tmp` may be RAM, and on 2026-10-27 earlier passes had left 37 GB of baseline builds and 5 GB of scratch behind.
**After you ship and before you write the pass's last log line, whether or not you shipped,** remove what *you*
made this pass, and list what earlier passes left.

1. **A baseline is a pair that can be found.** Make it with `git worktree add --detach <path> <sha>` (never a
   branch name: a worktree with a branch checked out stops anything else checking that branch out) and build it
   into `<path>-target`. Then `git worktree list`, plus the `*-target` directory beside each worktree, is the
   whole inventory.
2. **Remove yours.** `git worktree remove <path>` (plain, never `--force`) and `rm -rf <path>-target`, then your
   scratch under `/tmp` (captures, frames, venvs, gate logs, `rustdoctest*`). Remove only what you created this
   pass, and name each path in the log. If git refuses because the worktree has a change in it, name the file
   in the log and leave it.
3. **Leftovers from earlier passes are listed, never deleted.** A killed pass cannot clean up after itself, but
   an agent does not delete what it did not create, so list them once and let a person approve the removal:

   ```sh
   git worktree list --porcelain | awk '/^worktree /{sub(/^worktree /,""); print}' | tail -n +2   # all but the clone
   du -sh <wt> <wt>-target                       # the size of each one that is not yours
   git -C <wt> status --short                    # empty: nothing is lost by removing it
   find <wt> <wt>-target -mmin -120 -print -quit # empty: untouched for two hours, so no live build is using it
   ```

   Log each with its path, size, and whether it is clean, already on `origin/main` and idle (the likeliest to be
   stale). Say so, and still leave it.
4. **Never remove** the clone itself, its own `target/` (the build cache the next pass starts from), anything you
   did not create this pass, or anything a process is using.
5. Write one log line: what you removed and the space it freed (`df -h .` before and after), and the leftovers
   you listed.

## Rules

- One improvement per iteration, fully finished: researched, implemented, GREAT-signed-off,
  benchmarked, gated, committed, documented, logged.
- No feature ships without a written critique, a GREAT verdict and a clean performance check.
- If blocked, log the blocker and pick the next item. A failed gate or a visible defect is not a
  blocker: fix it. A failed push or merge is: log the exact error and stop. Never force-push.
