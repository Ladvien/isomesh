# What an edit reports: its own chunks, not the caller's queue

**Issue #13** (owner's review, 2026-10-30): `EditReport::dirty_chunks` is documented as the chunks containing at least one output-changed cell, "what actually goes on the dirty queue". `mark_edit` computes it as the region's chunks that are **in the caller's whole `DirtySet`** (`dirty.rs:457-460`). A game that drains its queue under a frame budget (`DirtySet::mesh_within_budget` keeps whatever it did not reach) carries a backlog into the next edit, and every backlog chunk inside the new edit's box is then counted as dirtied by an edit that changed nothing in it. `dirty_chunk_fraction()` is inflated the same way. The scan is also `region_chunks × queue length` per edit (`dirty.iter().any(..)` over a set that is kept sorted for binary search).

`mark_edit_subgrid` (iteration 3) already counts the distinct owners it dirtied in a local set. The two markers define the same field of the same struct two ways.

## Sources

| Kind | Stem | Where |
|---|---|---|
| Algorithm | `transvoxel_dissertation_lengyel2010` (E. Lengyel, *Voxel-Based Terrain for Real-Time Virtual Simulations*, UC Davis, 2010) | §3 p. 28: the set of blocks to regenerate is a function of the edited region |
| Shipping game | `s2015-advances-pdf-230-mb` p. 25 and `s2015-advances-compressed-34mb-version-also-availablehere` p. 36 (A. Evans, Media Molecule, *Learning from Failure*, the *Dreams* talk, SIGGRAPH 2015 Advances course notes; one talk, two indexed files) | the shipped renderer keeps, per voxel, the list of **edits** that touch it, and counts edits per voxel |
| Technical | `10.1017_cbo9780511809071` (Manning, Raghavan and Schütze, *Introduction to Information Retrieval*, 2008) | end of §1.2: a sorted list kept as a contiguous variable-length array, the structure `DirtySet` is, and what that buys |
| Technical, repository | `crates/isomesh/src/chunk/dirty.rs` (`mark_edit`, `tally`, `DirtySet`), `dirty/subgrid.rs` (`mark_edit_subgrid`), `dirty/tests.rs`; FINDINGS M-497, M-32 | what the two markers compute today, and the test that already pins the Subgrid one |

## Craft principle

Lengyel, p. 28, on what an edit obliges the mesher to do (re-read from the indexed text in this pass):

> "Whenever a region of the voxel map is edited, all blocks intersecting the region's bounding box have their triangle meshes regenerated."

The obligation is a function of the edit's region and the data: a block is regenerated because the edit reaches it, and nothing in the sentence depends on which blocks were already awaiting regeneration. A queue is the union of obligations (the same block asked for twice is regenerated once, which is why `DirtySet` is a set). So a number that says "how much did **this** edit add to the obligation" is a function of the edit and the two fields; a number computed against the queue is a function of the edit **and of everything the caller has not yet got to**, which is a different quantity under the same name.

Manning, Raghavan and Schütze, end of §1.2, on the structure a sorted list should have in memory (`DirtySet` is a sorted, deduplicated `Vec<ChunkId>`):

> "Variable length arrays win in space requirements by avoiding the overhead for pointers and in time requirements because their use of contiguous memory increases speed on modern processors with memory caches."

and

> "If updates are relatively infrequent, variable length arrays will be more compact and faster to traverse."

The representation is the right one for a queue that is walked each frame and edited once per carve. It also fixes what a membership test costs: one binary search, `O(log Q)`, because the array is sorted. `dirty.iter().any(|d| d == *id)` ignores that order and pays `O(Q)` per test, `O(R × Q)` per edit for `R` region chunks and a queue of `Q`. The fix removes the test rather than speeding it up, since what it tested was the wrong question.

## How shipping games do it

Evans, `s2015-advances-pdf-230-mb` p. 25, on how *Dreams* keeps a sculpted world as a list of edits and finds what each part of space has to evaluate:

> "building a list of edits that could possibly overlap each voxel, and then iteratively refining the voxels by splitting them and shortening the lists for each."

and `s2015-advances-compressed-34mb-version-also-availablehere` p. 36, the figure caption for the same structure:

> "...and this is a visualisation of the number of edits that touch the leaf voxels"

The shipped system counts edits **per voxel, against the edit list**: a leaf's number is how many edits reach it, and it is the same whether or not other work is queued. That is the shape `dirty_chunks` has to have (per edit, against the fields). This is the nearest shipped-practice evidence the corpus holds; neither this talk nor the C4 dissertation says how a re-mesh queue's own statistics are defined, so this note takes the *shape* from them and the definition from the field's documentation and from `mark_edit_subgrid` (see Pending index).

## Model adopted

`EditReport::dirty_chunks` is **the number of distinct chunks that own at least one output-changed cell of this edit**, whatever the caller's `DirtySet` already holds. Both markers do it the same way:

- a local `DirtySet` (`dirtied`) collects the owners of output-changed cells while the cells are walked, exactly as `region_chunks` collects the owners of every cell;
- `dirty_chunks = dirtied.len()`, `region_chunks = region_chunks.len()`;
- the caller's queue then receives `dirtied`, each id through `DirtySet::insert`, so the queue ends as the union of what it held and what the edit dirtied, sorted and deduplicated as before.

The tail that does this is one private function used by both markers, so the definition cannot drift apart a second time. The per-cell counters stay in the shared `tally`.

Options considered. The issue's own fix direction is the first.

| Option | `dirty_chunks` means | Verdict |
|---|---|---|
| A. This edit's dirtied owners (adopted) | chunks the edit made stale | the documented meaning, `dirty_chunk_fraction` stays "the fraction of the region's chunks that need re-meshing", the report is a function of `(before, after, region)` as every other counter in it is, and it is what `mark_edit_subgrid` already returns |
| B. Chunks newly added to the queue | queue growth | answers "how much work did I add", which a caller can read as `dirty.len()` before and after; makes the same edit report differently depending on the backlog, so the fraction stops being a property of the edit, and the two markers would then need to agree on a number that is not part of this struct's contract |
| C. The old rule | region chunks that are queued | counts chunks the edit did not touch, and costs `region × queue` |

No owner decision is open: the issue names A, the field's documentation says A and the Subgrid marker does A. B stays available to the caller without a new API.

## Parameters

None. The change adds one `DirtySet` allocation per call (the Subgrid marker already makes it) and removes the `region_chunks × queue length` scan. Nothing about which cells are marked, which chunks are queued or their order changes, so no extraction, golden hash or queue content moves.

## Validation targets

1. The issue's own test: a no-op edit over a region whose first chunk is already queued reports `dirty_chunks == 0` and `dirty_chunk_fraction() == 0.0` and keeps the queue (fails before the change: `dirty_chunks == 1`).
2. The worst case of the old rule: with **every** chunk the region overlaps queued, a no-op edit still reports 0 (the old rule reported `region_chunks`).
3. The report is a function of the edit: one real edit into an empty queue and into a queue holding half of the chunks it dirties plus a chunk outside the region gives equal `EditReport`s, and the queue ends as the union of the two, each id once (so a marker that counts but does not queue, or that queues but drops what was waiting, fails).
4. Every existing marking test passes unchanged: the queue contents, the Subgrid marker's own `dirty_chunks_counts_this_edit_and_not_the_queue`, `samples == 1` equal to `mark_edit`'s report field by field, and the plugin tests that count `report.dirty_chunks`.
5. Before and after, the 378 golden hashes, `cargo tree -p isomesh -e normal` and the benches' structural columns are unchanged, because no extractor expression is touched.
6. The doc comment of `dirty_chunks` states the definition and says what it is not (the length of the caller's queue).

## What great looks like

A game that drains a 300-chunk backlog over several frames and digs in the meantime reads `"3 of 8 chunks in the box were dirty"` from the HUD of its own dig, the same number it would read with an empty queue, and the cost of reading it is the cost of the cells it walked, not of the queue it holds. `game_dig` is the example that prints it.

## Deviations

- `dirty_chunks` gets **smaller** for a caller that passes a non-empty queue. That is the fix, and the only behavioural change: a caller that read the old value as "region chunks that are queued" read a number the field's documentation never promised. Marked as a bug fix, with the issue's evidence.
- `mark_edit` used to insert into the caller's queue while walking the cells; it now inserts at the end, as `mark_edit_subgrid` does. The queue after the call is identical (a sorted set), and the function returns no error after the size gate, so no caller can observe a half-marked queue either way.
- Not in scope: `game_dig` also reads `dirty.iter()` for its outline gizmo and its halo expansion, which includes the backlog by design (it shows what the queue will re-mesh).

## Pending index

- No source in the corpus defines the statistics of a re-mesh queue (what a "dirty count" per edit should mean when work is carried over a frame boundary). C4 and *Dreams* give the per-edit, per-region shape; the definition here is this crate's own and is a design default backed by its documentation and by the Subgrid marker.
- No corpus source on the merge-intersection of two sorted sets (Manning §1.3's postings merge is the standard reference and was not read in this pass). It is not needed here: the fix removes the intersection rather than speeding it up.
- New SIGGRAPH stems from `catalog_recent` that are about implicit shape editing (`sig2024_CNS-Edit_3D_Shape_Editing_via_Coupled_Neural_Shape_Optimization`, `sig2026_Implicit_Minimal_Surfaces_for_Bijective_Correspondences`) are listed in `corpus.md` but were not read for this note.
