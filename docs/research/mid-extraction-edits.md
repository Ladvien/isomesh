# An edit that arrives mid-extraction is re-extracted, and the test that says so must be able to fail

**Issue #12** (owner's review, 2026-10-30): `an_edit_during_extraction_is_requeued_rather_than_swallowed` in `bevy_isomesh/src/plugin/tests.rs` passes when the edit it names is swallowed. Its last three assertions are a `ChunkMesh`, no `NeedsRemesh` and no `MeshingTask`, and the swallowing implementation (one `.remove::<NeedsRemesh>()` in the `Ok` arm of `apply_finished_meshes`) leaves exactly that state. Nothing in the test counts a second extraction. The claim is in the README (`bevy_isomesh/README.md`, "Extraction runs on the `AsyncComputeTaskPool`…") and the plugin's own docs, so a regression of it would ship green: a carved wall that stays as it was.

## Sources

| Kind | Stem | Where |
|---|---|---|
| Algorithm | `transvoxel_dissertation_lengyel2010` (E. Lengyel, *Voxel-Based Terrain for Real-Time Virtual Simulations*, UC Davis, 2010) | Chapter 1, before §1.2 (Contributions): the mesh is regenerated from the voxel data after an edit, and the data "can often be modified automatically in unpredictable ways by an application" |
| Shipping game | `transvoxel_dissertation_lengyel2010` | Chapter 1, before §1.2: the C4 Engine's requirement that terrain be "built or damaged in real time in response to events taking place during a simulation"; p. 99, deployed |
| Technical | `10.1109_visual.2001.964533` (Lindstrom and Pascucci, *Visualization of Large Terrains Made Easy*, IEEE Visualization 2001) | §3, before §3.2 (View Culling): a refinement thread supplies the render thread a display list that is reused until a newly refined mesh is obtained |
| Technical | `10.1109_aitest.2019.00019` (Saha and Kanewala, *Fault Detection Effectiveness of Metamorphic Relations Developed for Testing Supervised Classifiers*, 2019) | the mutation-analysis section (PDF p. 4): mutation analysis as the measure of whether a check can fail |
| Technical, repository | `bevy_isomesh/src/plugin.rs` (`spawn_meshing_tasks`, `apply_finished_meshes`, `WaitingChunks`), `plugin/tests.rs`; FINDINGS M-171, M-497 and Part 5 | the contract under test and the rule about checks that cannot fail |

## Craft principle

Lengyel, chapter 1, on what an edit must produce:

> "Second, it would allow terrain to be built or damaged in real time in response to events taking place during a simulation. For example, an explosion could create a spherical hole in the voxel map that produces a crater in the regenerated terrain."

and on what the voxel data is allowed to be when the mesher reads it:

> "Our methods must function correctly for arbitrary input data. We can place no constraints on the voxel data used to generate a terrain surface because it can often be modified automatically in unpredictable ways by an application"

So the mesh is a function of the voxel data **at the moment it is extracted**, and a mesher cannot know which edits happened while it was running. The caller can: it edits, and it marks the chunk. What the crate owes in return is the one thing a function of the data cannot supply, that **an extraction which starts after the mark reads the data the mark was about.** `bevy_isomesh` pays this in three lines of `plugin.rs`: `NeedsRemesh` comes off when the task is spawned (not when the mesh lands), each task holds its own `Arc` to the field it was spawned with, and a re-marked chunk with no task in flight spawns a new one on the next frame with whatever field the volume holds then.

## How shipping games do it

Lindstrom and Pascucci, §3, on the pattern every asynchronous mesher uses:

> "it is common to parallelize the otherwise sequential, interleaved tasks of refinement and rendering as two asynchronous processes or threads [16, 22]. In this model, the render thread is periodically and asynchronously supplied with a list of geometry to render by the refinement thread. This “display list” is then used, and potentially reused over several frames, until a newly refined mesh is obtained."

That is stale-while-revalidate: the old mesh draws until a new one lands, and the guarantee that matters is the *until*. An edit that arrives while a mesh is being built produces a mesh that is stale on arrival; if the re-mark is cleared by that mesh's landing, the *until* never comes and the chunk is stale for good, with no error and a green mesh count. That is the failure issue #12 names. It is silent in the way M-497's was: nothing is wrong with any single mesh, only with which field it was built from.

C4 (Lengyel, p. 99) shipped destructible voxel terrain on this contract, and the dissertation does not need the case spelled out because its only mechanism is "regenerate the block from the current voxels".

## Model adopted

No change to the plugin. The contract is already as above, and the issue's reading of the code is right: the swallow is one line, and the existing assertions do not see it. The change is to the test, so that it observes the extraction and the data, not only the end state:

1. **Count the extractions after the re-mark.** `MeshStats::spawned` resets every frame, so the test sums it over the frames it drives. Exactly one spawn follows the re-mark: zero is the swallow, two would be a duplicate queued behind `WaitingChunks`' `Without<MeshingTask>`.
2. **Check the data, not only the count.** The field is replaced (`VoxelVolume::replace_field`) between the two extractions, as a game does, and the test compares each mesh the app shows with a direct `extract_chunk` of the field it should have been built from: the first mesh equals the **old** field's extraction (a running task finishes against the field it was spawned with, `replace_field`'s documented promise), the final mesh equals the **new** field's.
3. **Show the fixture can tell them apart.** The old and new reference meshes are both non-empty and differ, asserted in the test, so a field swap that changes nothing in this chunk cannot make assertion 2 vacuous. (M-107, M-112, M-113, M-115, M-117, M-119, M-120 and M-127 are the eight earlier instances of a fixture that could not exhibit the property it was chosen for.)
4. **Prove it can fail** by running the test against the mutation the issue gives, and against a second one that stales the data and not the count (the second extraction is spawned and reads the first extraction's field). Saha and Kanewala's measure: "mutants are a proxy to the real faults for comparing testing techniques", which is the standard under which a test that kills no mutant has not been shown to test anything. Both outputs go in FINDINGS.

The sequence the test relies on is fixed by the schedule and not by timing: the re-mark lands while the task component is on the entity, `apply_finished_meshes` only runs inside `App::update`, so the first task cannot land before the test's next `update` however fast the worker is; `spawn_meshing_tasks` then sees the chunk without a task and spawns the second in the following frame. The only clock is `drain_until`'s 10 s deadline.

## Parameters

- Fixture: `app(spin, 1)`: one chunk `[0, 0, 0]`, an 8-cell chunk of cell size 0.25 (spans 0 to 2 per axis), `SlowSphere` of radius 1.0 at `spin` 2000 for the old field, a plain `Sphere` of radius 1.5 for the new one. Both surfaces cross the chunk, at different places. Design default, not from a source.
- Extractor: the volume's default (`Extractor::MarchingCubes`); the reference meshes call `extract_chunk` with the same extractor, so the comparison is bit-for-bit on `f32` positions, which is sound because `isomesh` extraction is deterministic (golden hashes, M-31).
- No change to `MeshBudget` defaults.

## Validation targets

- V1: the strengthened test passes on today's plugin.
- V2: it fails on the issue's mutation (`.remove::<NeedsRemesh>()` in the `Ok` arm), with the assertion that names the missing second extraction.
- V3: it fails on a mutation that spawns the second extraction against the first extraction's field, with the assertion that names the data.
- V4: it does not depend on wall time: three consecutive full runs of the plugin test module pass, and the mutations fail every time.
- V5: the other plugin tests are unchanged and still pass.

## What great looks like

A regression anywhere in the re-queue path (a cleared mark, a cached field, a double queue) turns this test red with a message that says which of the three it was, and the README's sentence "Editing a chunk mid-extraction re-queues it rather than swallowing the edit" is backed by a test that has been seen to fail. On screen the property is the one every carving game depends on: dig while the wall is still meshing and the hole appears a frame or two later, never not at all.

## Deviations

None to a default. The issue offers the count (`respawned == 1`) and a "stronger form" that also checks content; both are built. The count uses `MeshStats::spawned` summed across frames, as the issue proposes, and no new counter is added to the plugin. Design default, owner may override: the new field is a different-radius sphere and not an empty field; an empty second mesh would also distinguish the two but gives a weaker comparison (zero vertices equals zero vertices for any field that misses the chunk).

## Pending index

- No document in the corpus on how shipped engines other than C4 queue, version or discard chunk-mesh jobs (Minecraft's chunk compile queue, Astroneer, Teardown). A search for "Minecraft chunk section mesh rebuilt when a block is modified" returned only unrelated hits; Aokana (`10.48550_arXiv.2505.02017`) describes the chunking and not the rebuild queue. Open question for a source: whether any engine tags a mesh job with the edit's version and drops a landing whose version is older than the chunk's mark, instead of applying it and re-queueing as this plugin does (this plugin shows the stale mesh for a frame or two; a dropped landing would show the older one for longer).
- No document on mutation testing for Rust or for game code; Saha and Kanewala is on classifiers and is cited for the measure, not for a method.
