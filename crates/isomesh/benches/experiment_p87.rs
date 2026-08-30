//! **P-87 — an octree navigation graph that repairs itself locally, built on the
//! sign field rather than on triangles.**
//!
//! Ticket: R-087. Pre-registered before this harness existed.
//!
//! ```bash
//! cargo bench --bench experiment_p87
//! ```
//!
//! Writes `docs/experiments/p-87.csv`.
//!
//! # The source, and what it actually claims
//!
//! Massonnat & Verbrugge, *Efficient Octree-based 3D Pathfinding*, IEEE CoG 2024
//! (`10.1109/CoG60054.2024.10645669`): an octree splits cells containing
//! obstacles; adjacent cells merge by a Hertel–Mehlhorn-inspired greedy method
//! **while preserving convexity**; A\* runs on the resulting coarse graph; and
//! dynamic environments get **local** octree and graph updates plus a
//! cell-repairing strategy. Measured on an Intel Core i7-12700H: octree update
//! 0.22–1.36 ms, local graph update 0.03 ms, about 1 ms total, with cell
//! reduction up to an order of magnitude (28,190 → 303).
//!
//! The two alternatives this repository already costed: PLMSS Morse–Smale at
//! 256³ is 4.40 s single-threaded and 0.36 s on 24 threads (the docs' own
//! "20–40× over a 16 ms frame"), and Recast rebuilds a tile by **voxelising
//! collision geometry** — so a voxel game would voxelise a mesh it generated
//! from voxels. This harness therefore never touches a triangle: every
//! classification below comes from the packed sign bitmap, which is the same
//! representation `dual.rs`'s `inside` prepass uses (R-039), and the same
//! `is_inside(value)` comparison rather than the IEEE sign bit.
//!
//! # SHARE, recomputed before a line of this was written
//!
//! The registration's SHARE line is *"C1 moves the whole navigation-rebuild
//! stage, currently unbuilt"* — share **1.0**. `✗51`'s arithmetic
//! (`1/(1 − s + s/factor)`) is what forecloses a ratio clause denominated in a
//! small share, and it does not bind here for a plainer reason than share 1.0:
//! **none of the three clauses is a speedup against a baseline.** C1 is an
//! absolute latency bar (2 ms) on a stage that does not exist yet, C2 is the
//! ratio of two directly counted sets, and C3 is the ratio of two directly
//! counted sets. All three were arithmetically reachable before the run, and
//! that is recorded here rather than discovered afterwards.
//!
//! What is **not** reachable is a comparison with the paper's ~1 ms: that is an
//! i7-12700H and this is a Zen 3, so `M-281` forbids it as a gate. The gate is
//! the registered 2 ms on this machine, `clock_mhz` is on every row (`M-280`),
//! and C2 and C3 are integer ratios and therefore machine-independent.
//!
//! # The structure, in one screen
//!
//! | stage | what it is | unit |
//! |---|---|---|
//! | sign bitmap | one bit per **sample**, packed 64-to-a-word along `x` | sample |
//! | cell state | 8 corner bits → `AIR` / `SOLID` / `MIXED` | cell |
//! | pyramid | 8 children → parent; uniform or `MIXED` | node |
//! | leaves | descend from the root, stop at the first uniform node | nav cell |
//! | regions | greedy axis merge of **free** leaves, union must be a box | nav cell |
//! | graph | face-adjacent leaves in different regions | edge |
//!
//! **A `MIXED` unit cell is an obstacle.** That is the paper's "cells containing
//! obstacles are split", taken to the leaf: at maximum depth there is nothing
//! left to split, so a cell the surface passes through cannot be navigable.
//!
//! **Convexity is preserved by construction and it is asserted, not argued.**
//! Two boxes merge only when they are face-contiguous on one axis and their
//! extents agree exactly on the other two, so the union is a box. The control is
//! that every region's box volume equals the summed volume of its member leaves
//! — a merge that produced a non-convex or non-tiling region fails that
//! immediately.
//!
//! # The trace is P-72's, inherited literally
//!
//! Eleven spherical brushes subtracted along a straight path across the middle of
//! the world, radius **6 cells**, each a separate edit, with the height probed
//! **per edit at that edit's own `x`** — P-72 recorded two void runs before
//! arriving at that shape, and the reason is in its header. The brush is in
//! cells and the path is in world coordinates, exactly as there.
//!
//! The world-size sweep holds the brush at 6 cells across 32³ … 256³, which is
//! **M-311's own protocol**: radius 6, identical at 33³, 65³ and 129³, dirty
//! 925 three times while the lattice grew 59.7×. C2's *"at fixed edit size"* is
//! therefore fixed in the octree's own unit, and the world grows 512× in cells.
//!
//! # Controls, each an assertion rather than a printed number
//!
//! - **VACUITY (registered).** `topology_changes` — the dig must change the
//!   free-space topology at least once, read from `isomesh::connectivity::Air`
//!   (crate source, maintained independently of everything above) and never from
//!   the brush. Asserted `>= 1`. A trace that dug in air would report a fast
//!   time for a no-op, which is `M-44` exactly, and which is how P-72's first
//!   two fixtures were caught.
//! - **The dirty set is the crate's, not the harness's.** `dirty_cells` is the
//!   sum of `EditReport::sign_changed_cells` from `isomesh::chunk::dirty::mark_edit`,
//!   and the harness's own count of cells incident to a flipped sign bit is
//!   asserted **equal** to it per edit. Two independent walks over the same
//!   edit; a harness that mis-locates the brush box disagrees at once (which is
//!   the defect P-72 found in itself).
//! - **The leaves partition the world.** Σ `(2^level)³` over live leaves equals
//!   `world_cells³` after the build and after every edit.
//! - **The regions tile the free leaves.** Σ region box volumes equals Σ free
//!   leaf volumes, and each region's box volume equals its own leaves' volumes.
//! - **The repair is deterministic.** Every integer count from rep 0 is asserted
//!   equal to reps 1 and 2, which re-run the identical trace from the snapshot
//!   (`M-36`).
//! - **Local repair is scored against a global re-merge.** `post_regions_incremental`
//!   against `post_regions_global` on the same final field: a local strategy that
//!   silently shreds the merge would still be fast, and that would be a worse
//!   result reported as a better one.

#![allow(
    clippy::cast_possible_truncation,
    clippy::cast_precision_loss,
    clippy::cast_sign_loss,
    clippy::needless_range_loop,
    clippy::too_many_lines
)]

mod common;

use std::collections::{BTreeMap, BTreeSet};
use std::time::Instant;

use isomesh::Sdf;
use isomesh::chunk::ChunkLayout;
use isomesh::chunk::dirty::{DirtySet, mark_edit};
use isomesh::connectivity::Air;
use isomesh::fields::{FbmTerrain, Gyroid};

/// World sizes in cells per axis, powers of two and multiples of
/// [`CHUNK_CELLS`].
///
/// 32 to 256 is **512× in cells**, which is what C2's *"does not grow with
/// world size"* has to be tested over. M-311 tested its own version of this
/// claim over 59.7×.
const WORLD_SIZES: [u32; 4] = [32, 64, 128, 256];

/// Chunk granularity, in cells per axis. C1 names 4³ and this is it.
///
/// It is load-bearing rather than decorative: the repair re-derives cell states
/// over **whole chunks** containing a flipped sign bit, because that is the
/// granularity a chunked world hands its consumers, and `reclassified_cells` is
/// therefore `nav_dirty_chunks × 64`.
const CHUNK_CELLS: u32 = 4;

/// World extent per axis, fixed across the sweep so a coarser world is a coarser
/// sampling of the *same* field.
const EXTENT: f64 = 4.0;

/// World origin. Centred on the reference fields' own domain centre, for the
/// reason P-72 records: the positive octant alone does not contain
/// `fbm_terrain`'s sheet.
const ORIGIN: f64 = -EXTENT * 0.5;

/// Brush radius in **cells**, so the brush is 6 cells at every world size —
/// M-311's protocol.
const BRUSH_CELLS: f64 = 6.0;

/// Edits in the trace. P-72's eleven.
const EDITS: usize = 11;

/// Traces per arm. The counts are integers and identical across reps (asserted);
/// the reps exist only so the millisecond figure is a median rather than one
/// sample on a governed CPU (`M-337`).
const REPS: usize = 3;

/// M-311's dirty-cell count for one radius-6 brush in a solid lattice, quoted
/// from `docs/experiments/p-23.csv` via FINDINGS `M-311`. Used only to report
/// `repair_per_m311`; the clause's own denominator is this run's `dirty_cells`.
const M311_DIRTY: f64 = 925.0;

/// No leaf / no region / no label.
const NONE: u32 = u32::MAX;

/// Every corner sample of the cell is air.
const AIR: u8 = 0;
/// Every corner sample of the cell is solid.
const SOLID: u8 = 1;
/// The surface passes through: navigable only if it can be split further, and at
/// the leaf it cannot, so it is an obstacle.
const MIXED: u8 = 2;

// ── the sign field ──────────────────────────────────────────────────────────

/// One bit per sample, packed 64-to-a-word along `x`: is that sample inside?
///
/// The same layout and the same `is_inside(value)` predicate as `dual.rs`'s
/// active-cell prepass (R-039), and for the same reason it is not the IEEE sign
/// bit: `-0.0` has the sign bit set while `-0.0 < 0.0` is false, and this
/// crate's convention is that exactly zero is **outside**.
struct Signs {
    words: Vec<u64>,
    /// Words per bitmap row, `dims[0].div_ceil(64)`.
    bit_row: usize,
    dims: [u32; 3],
}

impl Signs {
    fn new(dims: [u32; 3]) -> Self {
        let bit_row = (dims[0] as usize).div_ceil(64);
        let words = vec![0u64; bit_row * dims[1] as usize * dims[2] as usize];
        Self {
            words,
            bit_row,
            dims,
        }
    }

    #[inline]
    fn slot(&self, s: [u32; 3]) -> (usize, u64) {
        let row = (s[2] as usize * self.dims[1] as usize + s[1] as usize) * self.bit_row;
        (row + (s[0] as usize >> 6), 1u64 << (s[0] & 63))
    }

    #[inline]
    fn solid(&self, s: [u32; 3]) -> bool {
        let (w, m) = self.slot(s);
        self.words[w] & m != 0
    }

    #[inline]
    fn set(&mut self, s: [u32; 3], solid: bool) {
        let (w, m) = self.slot(s);
        if solid {
            self.words[w] |= m;
        } else {
            self.words[w] &= !m;
        }
    }

    /// `AIR`, `SOLID` or `MIXED` for the cell whose minimum corner is `c`.
    #[inline]
    fn cell_state(&self, c: [u32; 3]) -> u8 {
        let mut solid = 0u32;
        for dz in 0..2 {
            for dy in 0..2 {
                for dx in 0..2 {
                    if self.solid([c[0] + dx, c[1] + dy, c[2] + dz]) {
                        solid += 1;
                    }
                }
            }
        }
        match solid {
            0 => AIR,
            8 => SOLID,
            _ => MIXED,
        }
    }
}

// ── the octree, as a state pyramid over cells ───────────────────────────────

/// `levels[0]` is one state per cell; `levels[k]` is one state per `2^k` node.
///
/// A pyramid rather than pointers because the whole point is that a node's state
/// is a pure function of its eight children, so the update after an edit is a
/// bottom-up walk over the dirty cells' ancestors and nothing else — and because
/// a canonical, coordinate-derived subdivision is *"the same as the one built
/// from scratch"* after any sequence of updates, which is the property
/// `ChunkLayout::at_lod` is documented against.
struct Tree {
    levels: Vec<Vec<u8>>,
    n: u32,
    depth: usize,
}

impl Tree {
    fn new(n: u32) -> Self {
        assert!(n.is_power_of_two(), "world must be a power of two in cells");
        let depth = n.trailing_zeros() as usize;
        let mut levels = Vec::with_capacity(depth + 1);
        for k in 0..=depth {
            let m = (n >> k) as usize;
            levels.push(vec![MIXED; m * m * m]);
        }
        Self { levels, n, depth }
    }

    #[inline]
    fn side(&self, level: usize) -> u32 {
        self.n >> level
    }

    #[inline]
    fn idx(&self, level: usize, c: [u32; 3]) -> usize {
        let m = (self.n >> level) as usize;
        (c[2] as usize * m + c[1] as usize) * m + c[0] as usize
    }

    #[inline]
    fn state(&self, level: usize, c: [u32; 3]) -> u8 {
        self.levels[level][self.idx(level, c)]
    }

    /// The whole pyramid from the bitmap. `O(cells)`, once per rep.
    fn rebuild(&mut self, signs: &Signs) {
        let n = self.n as usize;
        for z in 0..n {
            for y in 0..n {
                for x in 0..n {
                    self.levels[0][(z * n + y) * n + x] =
                        signs.cell_state([x as u32, y as u32, z as u32]);
                }
            }
        }
        for l in 1..=self.depth {
            let m = (self.n >> l) as usize;
            let mc = m * 2;
            let (lower, upper) = self.levels.split_at_mut(l);
            let child = &lower[l - 1];
            let parent = &mut upper[0];
            for z in 0..m {
                for y in 0..m {
                    for x in 0..m {
                        let mut kids = [MIXED; 8];
                        for dz in 0..2 {
                            for dy in 0..2 {
                                for dx in 0..2 {
                                    kids[(dz * 2 + dy) * 2 + dx] = child
                                        [((2 * z + dz) * mc + (2 * y + dy)) * mc + (2 * x + dx)];
                                }
                            }
                        }
                        parent[(z * m + y) * m + x] = combine(kids);
                    }
                }
            }
        }
    }
}

/// A node is uniform only when all eight children agree and are not `MIXED`.
#[inline]
fn combine(kids: [u8; 8]) -> u8 {
    let f = kids[0];
    if f != MIXED && kids.iter().all(|&s| s == f) {
        f
    } else {
        MIXED
    }
}

// ── boxes, leaves, regions ──────────────────────────────────────────────────

/// A half-open box in cell coordinates. `hi` is exclusive.
#[derive(Clone, Copy, PartialEq, Eq)]
struct Box3 {
    lo: [u32; 3],
    hi: [u32; 3],
}

impl Box3 {
    #[inline]
    fn volume(&self) -> u64 {
        let mut v = 1u64;
        for a in 0..3 {
            v *= u64::from(self.hi[a] - self.lo[a]);
        }
        v
    }
}

/// One octree leaf: a uniform node, or a `MIXED` unit cell.
struct Leaf {
    level: u8,
    state: u8,
    coords: [u32; 3],
    /// Owning merged region, or [`NONE`] for a leaf that is not free.
    region: u32,
    alive: bool,
}

/// One merged convex navigation cell.
struct Region {
    bx: Box3,
    leaves: Vec<u32>,
    alive: bool,
}

/// What one edit's local repair touched.
#[derive(Clone, Copy, Default, PartialEq, Eq, Debug)]
struct Cost {
    /// Sample bits that flipped.
    flipped: u64,
    /// Cells incident to a flipped bit — the harness's own dirty set, asserted
    /// equal to `mark_edit`'s `sign_changed_cells`.
    pattern_cells: u64,
    /// 4³ chunks containing a flipped bit.
    nav_chunks: u64,
    /// Cells whose state was recomputed: `nav_chunks × 64`.
    reclassified: u64,
    /// Cells whose `AIR`/`SOLID`/`MIXED` state actually moved.
    state_cells: u64,
    /// Pyramid nodes above level 0 whose state was recomputed.
    nodes: u64,
    leaves_removed: u64,
    leaves_added: u64,
    regions_removed: u64,
    regions_added: u64,
    /// Leaves handed to the local greedy re-merge — the cost driver R-088 needs
    /// to know about.
    pool: u64,
}

impl Cost {
    /// The registered repair set, in navigation cells.
    fn repair_cells(&self) -> u64 {
        self.leaves_removed + self.leaves_added + self.regions_removed + self.regions_added
    }
}

// ── the navigation structure ────────────────────────────────────────────────

/// Octree leaves, merged convex regions, and the adjacency graph over them.
struct Nav {
    tree: Tree,
    /// Per level, node index → leaf id, or [`NONE`] when the node is not a leaf.
    ///
    /// This is what makes point location and neighbour-finding `O(depth)`
    /// **without** a dense cell → region map. A dense per-cell map would make
    /// splitting one big region a world-proportional rewrite, which is precisely
    /// the growth C2 forbids, so it is not used.
    leaf_at: Vec<Vec<u32>>,
    leaves: Vec<Leaf>,
    leaf_free: Vec<u32>,
    /// Leaf ids killed by the repair in progress, held back from
    /// [`Nav::leaf_free`] until it finishes.
    ///
    /// **This is a bug the region-tiling audit caught on the first run.**
    /// Recycling an id inside one repair makes a dissolved region's member list
    /// point at a *different* node — the fresh leaf that took the id — so the
    /// surviving-leaf collection picks up a leaf that was never in that region
    /// and, when the fresh leaf happens to be `SOLID`, hands a non-free leaf to
    /// the free-space merge. It reported "regions cover 10854 of 10831 free
    /// cells" on `fbm_terrain` at 32³, which is the tiling control doing exactly
    /// what it exists for.
    leaf_retire: Vec<u32>,
    regions: Vec<Region>,
    region_free: Vec<u32>,
    region_edges: Vec<BTreeSet<u32>>,
    live_leaves: u64,
    live_regions: u64,
    edges: u64,
}

impl Nav {
    fn new(n: u32) -> Self {
        let tree = Tree::new(n);
        let mut leaf_at = Vec::with_capacity(tree.depth + 1);
        for k in 0..=tree.depth {
            let m = (n >> k) as usize;
            leaf_at.push(vec![NONE; m * m * m]);
        }
        Self {
            tree,
            leaf_at,
            leaves: Vec::new(),
            leaf_free: Vec::new(),
            leaf_retire: Vec::new(),
            regions: Vec::new(),
            region_free: Vec::new(),
            region_edges: Vec::new(),
            live_leaves: 0,
            live_regions: 0,
            edges: 0,
        }
    }

    fn clear(&mut self) {
        for l in &mut self.leaf_at {
            l.fill(NONE);
        }
        self.leaves.clear();
        self.leaf_free.clear();
        self.leaf_retire.clear();
        self.regions.clear();
        self.region_free.clear();
        self.region_edges.clear();
        self.live_leaves = 0;
        self.live_regions = 0;
        self.edges = 0;
    }

    #[inline]
    fn leaf_box(&self, id: u32) -> Box3 {
        let leaf = &self.leaves[id as usize];
        let s = 1u32 << leaf.level;
        Box3 {
            lo: [leaf.coords[0] * s, leaf.coords[1] * s, leaf.coords[2] * s],
            hi: [
                (leaf.coords[0] + 1) * s,
                (leaf.coords[1] + 1) * s,
                (leaf.coords[2] + 1) * s,
            ],
        }
    }

    fn new_leaf(&mut self, level: usize, coords: [u32; 3], state: u8) -> u32 {
        let id = match self.leaf_free.pop() {
            Some(id) => {
                self.leaves[id as usize] = Leaf {
                    level: level as u8,
                    state,
                    coords,
                    region: NONE,
                    alive: true,
                };
                id
            }
            None => {
                self.leaves.push(Leaf {
                    level: level as u8,
                    state,
                    coords,
                    region: NONE,
                    alive: true,
                });
                (self.leaves.len() - 1) as u32
            }
        };
        let i = self.tree.idx(level, coords);
        self.leaf_at[level][i] = id;
        self.live_leaves += 1;
        id
    }

    fn kill_leaf(&mut self, id: u32) {
        let (level, coords) = {
            let leaf = &self.leaves[id as usize];
            (leaf.level as usize, leaf.coords)
        };
        let i = self.tree.idx(level, coords);
        self.leaf_at[level][i] = NONE;
        self.leaves[id as usize].alive = false;
        self.leaves[id as usize].region = NONE;
        self.leaf_retire.push(id);
        self.live_leaves -= 1;
    }

    /// Create leaves for the subtree at `(level, coords)`, stopping at the first
    /// uniform node. Level 0 is always a leaf.
    fn derive_under(&mut self, level: usize, coords: [u32; 3], added: &mut Vec<u32>) {
        let mut stack = vec![(level, coords)];
        while let Some((l, c)) = stack.pop() {
            let st = self.tree.state(l, c);
            if st != MIXED || l == 0 {
                let id = self.new_leaf(l, c, st);
                added.push(id);
                continue;
            }
            for dz in 0..2 {
                for dy in 0..2 {
                    for dx in 0..2 {
                        stack.push((l - 1, [2 * c[0] + dx, 2 * c[1] + dy, 2 * c[2] + dz]));
                    }
                }
            }
        }
    }

    /// The existing leaves covering the subtree at `(level, coords)`.
    fn collect_under(&self, level: usize, coords: [u32; 3], out: &mut Vec<u32>) {
        let mut stack = vec![(level, coords)];
        while let Some((l, c)) = stack.pop() {
            let id = self.leaf_at[l][self.tree.idx(l, c)];
            if id != NONE {
                out.push(id);
                continue;
            }
            assert!(
                l > 0,
                "a level-0 node with no leaf and no uniform ancestor: the leaf \
                 structure is not a partition, which means the repair is wrong"
            );
            for dz in 0..2 {
                for dy in 0..2 {
                    for dx in 0..2 {
                        stack.push((l - 1, [2 * c[0] + dx, 2 * c[1] + dy, 2 * c[2] + dz]));
                    }
                }
            }
        }
    }

    /// Leaves sharing a face with `id` in the `+axis` or `-axis` direction.
    ///
    /// Three cases and they are exhaustive: the neighbour node is itself a leaf;
    /// it is uniform but covered by a coarser leaf (walk up); or it is `MIXED`
    /// above level 0, which — because a uniform node's children are all uniform
    /// and therefore no `MIXED` node can have a uniform ancestor — means it is
    /// subdivided, so descend into the four children touching the shared face.
    fn neighbours(&self, id: u32, axis: usize, positive: bool, out: &mut Vec<u32>) {
        let (level, coords) = {
            let leaf = &self.leaves[id as usize];
            (leaf.level as usize, leaf.coords)
        };
        let mut c = coords;
        if positive {
            if c[axis] + 1 >= self.tree.side(level) {
                return;
            }
            c[axis] += 1;
        } else {
            if c[axis] == 0 {
                return;
            }
            c[axis] -= 1;
        }

        let here = self.leaf_at[level][self.tree.idx(level, c)];
        if here != NONE {
            out.push(here);
            return;
        }
        if self.tree.state(level, c) != MIXED {
            // Uniform but not a leaf: an ancestor is the leaf.
            let mut l = level;
            let mut cc = c;
            while l < self.tree.depth {
                l += 1;
                cc = [cc[0] >> 1, cc[1] >> 1, cc[2] >> 1];
                let up = self.leaf_at[l][self.tree.idx(l, cc)];
                if up != NONE {
                    out.push(up);
                    return;
                }
            }
            panic!("a uniform node with neither a leaf nor a leaf ancestor");
        }
        // Subdivided: descend into the children on the shared face.
        let mut stack = vec![(level, c)];
        while let Some((l, cc)) = stack.pop() {
            let leaf = self.leaf_at[l][self.tree.idx(l, cc)];
            if leaf != NONE {
                out.push(leaf);
                continue;
            }
            assert!(l > 0, "descended past level 0 looking for a face neighbour");
            // The shared face is the neighbour's low face when we stepped in the
            // positive direction, and its high face otherwise.
            let want = u32::from(!positive);
            for d0 in 0..2 {
                for d1 in 0..2 {
                    let mut off = [0u32; 3];
                    off[axis] = want;
                    off[(axis + 1) % 3] = d0;
                    off[(axis + 2) % 3] = d1;
                    stack.push((
                        l - 1,
                        [2 * cc[0] + off[0], 2 * cc[1] + off[1], 2 * cc[2] + off[2]],
                    ));
                }
            }
        }
    }

    fn new_region(&mut self, bx: Box3, leaves: Vec<u32>) -> u32 {
        let id = match self.region_free.pop() {
            Some(id) => {
                self.regions[id as usize] = Region {
                    bx,
                    leaves,
                    alive: true,
                };
                self.region_edges[id as usize].clear();
                id
            }
            None => {
                self.regions.push(Region {
                    bx,
                    leaves,
                    alive: true,
                });
                self.region_edges.push(BTreeSet::new());
                (self.regions.len() - 1) as u32
            }
        };
        let members = core::mem::take(&mut self.regions[id as usize].leaves);
        for &l in &members {
            self.leaves[l as usize].region = id;
        }
        self.regions[id as usize].leaves = members;
        self.live_regions += 1;
        id
    }

    fn kill_region(&mut self, id: u32) {
        let peers: Vec<u32> = self.region_edges[id as usize].iter().copied().collect();
        for p in peers {
            self.region_edges[p as usize].remove(&id);
            self.region_edges[id as usize].remove(&p);
            self.edges -= 1;
        }
        let members = core::mem::take(&mut self.regions[id as usize].leaves);
        for l in members {
            if self.leaves[l as usize].alive && self.leaves[l as usize].region == id {
                self.leaves[l as usize].region = NONE;
            }
        }
        self.regions[id as usize].alive = false;
        self.region_free.push(id);
        self.live_regions -= 1;
    }

    fn add_edge(&mut self, a: u32, b: u32) {
        if a == b {
            return;
        }
        if self.region_edges[a as usize].insert(b) {
            self.region_edges[b as usize].insert(a);
            self.edges += 1;
        }
    }

    /// Merge a pool of free leaves into convex regions and link them.
    fn merge_pool(&mut self, pool: &[u32], link_both_ways: bool) -> u64 {
        let mut items: Vec<(Box3, Vec<u32>)> =
            pool.iter().map(|&l| (self.leaf_box(l), vec![l])).collect();
        greedy_merge(&mut items);
        let mut made = Vec::with_capacity(items.len());
        for (bx, leaves) in items {
            made.push(self.new_region(bx, leaves));
        }
        let mut scratch = Vec::new();
        for &r in &made {
            let members: Vec<u32> = self.regions[r as usize].leaves.clone();
            for l in members {
                for axis in 0..3 {
                    for &positive in if link_both_ways {
                        &[true, false][..]
                    } else {
                        &[true][..]
                    } {
                        scratch.clear();
                        self.neighbours(l, axis, positive, &mut scratch);
                        for &m in &scratch {
                            let other = &self.leaves[m as usize];
                            if other.state != AIR || other.region == NONE {
                                continue;
                            }
                            let o = other.region;
                            self.add_edge(r, o);
                        }
                    }
                }
            }
        }
        made.len() as u64
    }

    /// Full build from the pyramid: leaves, global merge, graph.
    fn build(&mut self) {
        self.clear();
        let mut added = Vec::new();
        let (depth, root) = (self.tree.depth, [0u32; 3]);
        self.derive_under(depth, root, &mut added);
        let pool: Vec<u32> = added
            .into_iter()
            .filter(|&l| self.leaves[l as usize].state == AIR)
            .collect();
        // `link_both_ways = false`: every leaf is scanned, so each face contact
        // is seen exactly once, from its lower side.
        self.merge_pool(&pool, false);
    }

    /// Connected components of the region graph, by union-find over live regions.
    fn components(&self) -> u64 {
        let mut parent: Vec<u32> = (0..self.regions.len() as u32).collect();
        fn find(parent: &mut [u32], mut x: u32) -> u32 {
            while parent[x as usize] != x {
                let g = parent[parent[x as usize] as usize];
                parent[x as usize] = g;
                x = g;
            }
            x
        }
        for (i, edges) in self.region_edges.iter().enumerate() {
            if !self.regions[i].alive {
                continue;
            }
            for &j in edges {
                let (a, b) = (find(&mut parent, i as u32), find(&mut parent, j));
                if a != b {
                    parent[a as usize] = b;
                }
            }
        }
        let mut roots = BTreeSet::new();
        for i in 0..self.regions.len() {
            if self.regions[i].alive {
                roots.insert(find(&mut parent, i as u32));
            }
        }
        roots.len() as u64
    }

    /// The static build's counts, for C3 and for the mechanism.
    fn snapshot(&self) -> Static {
        let mut free_leaves = 0u64;
        let mut free_cells = 0u64;
        for (i, leaf) in self.leaves.iter().enumerate() {
            if leaf.alive && leaf.state == AIR {
                free_leaves += 1;
                free_cells += self.leaf_box(i as u32).volume();
            }
        }
        Static {
            all_leaves: self.live_leaves,
            free_leaves,
            regions: self.live_regions,
            edges: self.edges,
            components: self.components(),
            free_cells,
        }
    }

    /// Every invariant the structure has, checked at once. Untimed.
    fn audit(&self, where_: &str) {
        let mut leaf_volume = 0u64;
        let mut free_volume = 0u64;
        for (i, leaf) in self.leaves.iter().enumerate() {
            if !leaf.alive {
                continue;
            }
            let v = self.leaf_box(i as u32).volume();
            leaf_volume += v;
            if leaf.state == AIR {
                free_volume += v;
                assert!(
                    leaf.region != NONE && self.regions[leaf.region as usize].alive,
                    "{where_}: free leaf {i} has no live region"
                );
            }
        }
        let world = u64::from(self.tree.n).pow(3);
        assert_eq!(
            leaf_volume, world,
            "{where_}: leaves cover {leaf_volume} of {world} cells, so they are not a partition"
        );
        let mut region_volume = 0u64;
        for (i, region) in self.regions.iter().enumerate() {
            if !region.alive {
                continue;
            }
            let members: u64 = region
                .leaves
                .iter()
                .map(|&l| self.leaf_box(l).volume())
                .sum();
            // The convexity instrument: a region is a box, and its box volume
            // equals its members' volumes exactly, so the merge produced a
            // convex region that its leaves tile without gap or overlap.
            assert_eq!(
                region.bx.volume(),
                members,
                "{where_}: region {i} has box volume {} against member volume {members}: the \
                 greedy merge produced a region that is not an exactly tiled box, so it is not \
                 convexity-preserving",
                region.bx.volume()
            );
            region_volume += members;
        }
        assert_eq!(
            region_volume, free_volume,
            "{where_}: regions cover {region_volume} of {free_volume} free cells"
        );
    }
}

/// Hertel–Mehlhorn-inspired greedy merge, restricted to axis-aligned boxes.
///
/// Hertel & Mehlhorn's rule for polygons is *"remove a diagonal if the result is
/// still convex"*. For boxes the convexity test is exact and cheap: the union of
/// two boxes is convex **iff** they are contiguous on one axis and their extents
/// agree on the other two, in which case the union is itself a box. So the pass
/// groups boxes by their cross-section on the two other axes and merges maximal
/// contiguous runs, repeating over the three axes until a whole sweep merges
/// nothing.
///
/// Deterministic: the groups come out of a `BTreeMap` in key order and each
/// group is sorted by its position on the merge axis. Determinism is load-bearing
/// here for the same reason it is in `connectivity` (`M-36`) — the region ids
/// this hands `R-088` have to be reproducible.
fn greedy_merge(items: &mut Vec<(Box3, Vec<u32>)>) {
    loop {
        let mut merged_any = false;
        for axis in 0..3 {
            let b = (axis + 1) % 3;
            let c = (axis + 2) % 3;
            let mut groups: BTreeMap<[u32; 4], Vec<usize>> = BTreeMap::new();
            for (i, it) in items.iter().enumerate() {
                groups
                    .entry([it.0.lo[b], it.0.hi[b], it.0.lo[c], it.0.hi[c]])
                    .or_default()
                    .push(i);
            }
            let mut dead = vec![false; items.len()];
            let mut any_here = false;
            for (_, mut group) in groups {
                if group.len() < 2 {
                    continue;
                }
                group.sort_unstable_by_key(|&i| items[i].0.lo[axis]);
                let mut head = group[0];
                for &i in &group[1..] {
                    if items[head].0.hi[axis] == items[i].0.lo[axis] {
                        items[head].0.hi[axis] = items[i].0.hi[axis];
                        let taken = core::mem::take(&mut items[i].1);
                        items[head].1.extend(taken);
                        dead[i] = true;
                        any_here = true;
                    } else {
                        head = i;
                    }
                }
            }
            if any_here {
                merged_any = true;
                let old = core::mem::take(items);
                let mut out = Vec::with_capacity(old.len());
                for (i, it) in old.into_iter().enumerate() {
                    if !dead[i] {
                        out.push(it);
                    }
                }
                *items = out;
            }
        }
        if !merged_any {
            return;
        }
    }
}

// ── the local repair ────────────────────────────────────────────────────────

/// Apply one edit's new sign bits and repair the octree and the graph locally.
///
/// Returns `(octree_ns, graph_ns)`. The field evaluations that produced
/// `new_solid` are **not** charged here: a chunked world re-samples an edited
/// chunk to re-mesh it whatever the navigation layer does, so charging the
/// navigation repair for them would overstate a cost it does not add. Everything
/// else — the bit diff, the chunk marking, the cell reclassification, the
/// pyramid walk, the leaf and region surgery, the relinking — is inside the
/// timers.
#[allow(clippy::too_many_arguments)]
fn repair(
    nav: &mut Nav,
    signs: &mut Signs,
    lo: [u32; 3],
    ext: [u32; 3],
    new_solid: &[bool],
    flipped: &mut Vec<[u32; 3]>,
    cost: &mut Cost,
) -> (u128, u128) {
    let n = nav.tree.n;
    let chunks_per_axis = n / CHUNK_CELLS;

    let t = Instant::now();

    // ── bits ────────────────────────────────────────────────────────────────
    flipped.clear();
    for z in 0..ext[2] {
        for y in 0..ext[1] {
            for x in 0..ext[0] {
                let s = [lo[0] + x, lo[1] + y, lo[2] + z];
                let want = new_solid[((z * ext[1] + y) * ext[0] + x) as usize];
                if signs.solid(s) != want {
                    signs.set(s, want);
                    flipped.push(s);
                }
            }
        }
    }
    cost.flipped = flipped.len() as u64;

    // ── the chunks a flipped bit touches, and the cells incident to one ─────
    //
    // A sample belongs to up to eight cells, indices `s - 1` and `s` per axis.
    let mut chunk_ids: Vec<u32> = Vec::new();
    let mut pattern: Vec<u32> = Vec::new();
    for &s in flipped.iter() {
        for dz in 0..2u32 {
            for dy in 0..2u32 {
                for dx in 0..2u32 {
                    let mut cell = [0u32; 3];
                    let mut ok = true;
                    for (a, d) in [dx, dy, dz].into_iter().enumerate() {
                        let v = s[a] + d;
                        if v == 0 || v > n {
                            ok = false;
                            break;
                        }
                        cell[a] = v - 1;
                    }
                    if !ok {
                        continue;
                    }
                    pattern.push(nav.tree.idx(0, cell) as u32);
                    let ch = [
                        cell[0] / CHUNK_CELLS,
                        cell[1] / CHUNK_CELLS,
                        cell[2] / CHUNK_CELLS,
                    ];
                    chunk_ids.push((ch[2] * chunks_per_axis + ch[1]) * chunks_per_axis + ch[0]);
                }
            }
        }
    }
    pattern.sort_unstable();
    pattern.dedup();
    cost.pattern_cells = pattern.len() as u64;
    chunk_ids.sort_unstable();
    chunk_ids.dedup();
    cost.nav_chunks = chunk_ids.len() as u64;
    cost.reclassified = cost.nav_chunks * u64::from(CHUNK_CELLS).pow(3);

    // ── reclassify whole chunks ─────────────────────────────────────────────
    let mut changed: Vec<u32> = Vec::new();
    for &ch in &chunk_ids {
        let cz = ch / (chunks_per_axis * chunks_per_axis);
        let cy = (ch / chunks_per_axis) % chunks_per_axis;
        let cx = ch % chunks_per_axis;
        let base = [cx * CHUNK_CELLS, cy * CHUNK_CELLS, cz * CHUNK_CELLS];
        for dz in 0..CHUNK_CELLS {
            for dy in 0..CHUNK_CELLS {
                for dx in 0..CHUNK_CELLS {
                    let cell = [base[0] + dx, base[1] + dy, base[2] + dz];
                    let st = signs.cell_state(cell);
                    let i = nav.tree.idx(0, cell);
                    if nav.tree.levels[0][i] != st {
                        nav.tree.levels[0][i] = st;
                        changed.push(i as u32);
                    }
                }
            }
        }
    }
    cost.state_cells = changed.len() as u64;

    // ── bottom-up pyramid update ────────────────────────────────────────────
    //
    // The candidates at level `l` are the parents of the level `l-1` nodes whose
    // state actually moved, so an edit that does not reach a node never touches
    // it. This is the whole of C2's mechanism: the walk is over the edit's
    // ancestors, and the world contributes only `depth`.
    let depth = nav.tree.depth;
    let mut changed_levels: Vec<Vec<u32>> = vec![Vec::new(); depth + 1];
    changed_levels[0] = changed;
    for l in 1..=depth {
        let m = nav.tree.side(l) as usize;
        let mc = m * 2;
        let mut cand: Vec<u32> = Vec::with_capacity(changed_levels[l - 1].len());
        for &i in &changed_levels[l - 1] {
            let i = i as usize;
            let x = i % mc;
            let y = (i / mc) % mc;
            let z = i / (mc * mc);
            cand.push((((z / 2) * m + y / 2) * m + x / 2) as u32);
        }
        cand.sort_unstable();
        cand.dedup();
        cost.nodes += cand.len() as u64;
        let mut moved = Vec::with_capacity(cand.len());
        for &i in &cand {
            let i = i as usize;
            let x = i % m;
            let y = (i / m) % m;
            let z = i / (m * m);
            let mut kids = [MIXED; 8];
            for dz in 0..2 {
                for dy in 0..2 {
                    for dx in 0..2 {
                        kids[(dz * 2 + dy) * 2 + dx] = nav.tree.levels[l - 1]
                            [((2 * z + dz) * mc + (2 * y + dy)) * mc + (2 * x + dx)];
                    }
                }
            }
            let st = combine(kids);
            if nav.tree.levels[l][i] != st {
                nav.tree.levels[l][i] = st;
                moved.push(i as u32);
            }
        }
        changed_levels[l] = moved;
    }

    // ── the top-most nodes whose state moved ────────────────────────────────
    //
    // A node changes state only when its parent either changed with it (uniform
    // parent absorbing a uniform child) or was `MIXED` and stayed `MIXED`. So
    // the nodes whose *leaf* structure has to be rebuilt are exactly those with
    // no changed parent, and everything below them follows.
    let mut tops: Vec<(usize, [u32; 3])> = Vec::new();
    for l in 0..=depth {
        let m = nav.tree.side(l) as usize;
        for &i in &changed_levels[l] {
            let i = i as usize;
            let x = i % m;
            let y = (i / m) % m;
            let z = i / (m * m);
            if l < depth {
                let mp = m / 2;
                let pi = (((z / 2) * mp + y / 2) * mp + x / 2) as u32;
                if changed_levels[l + 1].binary_search(&pi).is_ok() {
                    continue;
                }
            }
            tops.push((l, [x as u32, y as u32, z as u32]));
        }
    }
    // `changed_levels[l+1]` is only sorted for `l + 1 >= 1`, which is every
    // level the lookup above reaches, because level 0 is never a parent.
    debug_assert!(changed_levels.iter().skip(1).all(|v| v.is_sorted()));

    // ── leaf surgery ────────────────────────────────────────────────────────
    let mut old_leaves: Vec<u32> = Vec::new();
    for &(l, c) in &tops {
        nav.collect_under(l, c, &mut old_leaves);
    }
    old_leaves.sort_unstable();
    old_leaves.dedup();
    cost.leaves_removed = old_leaves.len() as u64;

    // Regions that lose a leaf must be dissolved: a region is a set of leaves
    // and it cannot survive one of them ceasing to exist.
    let mut doomed: Vec<u32> = Vec::new();
    for &l in &old_leaves {
        let r = nav.leaves[l as usize].region;
        if r != NONE {
            doomed.push(r);
        }
    }
    doomed.sort_unstable();
    doomed.dedup();

    for &l in &old_leaves {
        nav.kill_leaf(l);
    }
    let mut new_leaves: Vec<u32> = Vec::new();
    for &(l, c) in &tops {
        nav.derive_under(l, c, &mut new_leaves);
    }
    cost.leaves_added = new_leaves.len() as u64;

    let octree_ns = t.elapsed().as_nanos();

    // ── graph ───────────────────────────────────────────────────────────────
    let t = Instant::now();
    let mut pool: Vec<u32> = new_leaves
        .iter()
        .copied()
        .filter(|&l| nav.leaves[l as usize].state == AIR)
        .collect();
    for &r in &doomed {
        if !nav.regions[r as usize].alive {
            continue;
        }
        for l in nav.regions[r as usize].leaves.clone() {
            if nav.leaves[l as usize].alive {
                pool.push(l);
            }
        }
        nav.kill_region(r);
        cost.regions_removed += 1;
    }
    pool.sort_unstable();
    pool.dedup();
    cost.pool = pool.len() as u64;
    cost.regions_added = nav.merge_pool(&pool, true);
    let graph_ns = t.elapsed().as_nanos();

    // Only now may killed ids be handed out again: within one repair a
    // dissolved region's member list still names them.
    let mut retired = core::mem::take(&mut nav.leaf_retire);
    nav.leaf_free.append(&mut retired);
    nav.leaf_retire = retired;

    (octree_ns, graph_ns)
}

// ── the trace ───────────────────────────────────────────────────────────────

/// A sphere subtracted from a field: `max(field, -(|p - c| - r))`.
///
/// P-72's, verbatim, because the trace has to be the same one: the edit needs to
/// exist as two separate `Sdf`s, before and after, for `mark_edit`.
struct Dug<'a, F> {
    field: &'a F,
    centres: &'a [[f64; 3]],
    radius: f64,
}

impl<F: Sdf<Scalar = f64>> Sdf for Dug<'_, F> {
    type Scalar = f64;

    #[inline]
    fn sample(&self, p: [f64; 3]) -> f64 {
        let mut v = self.field.sample(p);
        for c in self.centres {
            let d = [p[0] - c[0], p[1] - c[1], p[2] - c[2]];
            let sphere = (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt() - self.radius;
            v = v.max(-sphere);
        }
        v
    }
}

/// The structure one arm builds before any edit.
///
/// Identical for both traces, because both start from the same pristine bitmap,
/// and asserted equal across them — a free control on the rebuild path.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
struct Static {
    all_leaves: u64,
    free_leaves: u64,
    regions: u64,
    edges: u64,
    components: u64,
    /// Cells covered by free leaves. With `free_leaves` this is the merge's real
    /// cost driver: leaves per free cell is surface density.
    free_cells: u64,
}

/// One trace on one field at one world size.
struct Arm {
    field: &'static str,
    trace: &'static str,
    world_cells: u32,
    depth: usize,
    brush_cells: f64,
    stat: Static,
    dirty_cells: u64,
    state_cells: u64,
    pattern_cells: u64,
    nav_chunks: u64,
    mesh_chunks: u64,
    reclassified: u64,
    nodes: u64,
    leaves_removed: u64,
    leaves_added: u64,
    regions_removed: u64,
    regions_added: u64,
    max_pool: u64,
    repair_cells: u64,
    /// Worst single edit, and its two components, from one rep.
    worst_edit: usize,
    octree_ms: f64,
    graph_ms: f64,
    total_ms: f64,
    median_edit_ms: f64,
    trace_total_ms: f64,
    graph_components_after: u64,
    post_regions_incremental: u64,
    post_regions_global: u64,
    topology_changes: u64,
    edits_joining: u64,
    edits_pocketing: u64,
    air_merges: u64,
    air_newly_air: u64,
    air_components_before: u64,
    air_components_after: u64,
    air_components_max: u64,
    /// From-scratch cost of the whole structure: pyramid, leaves, global merge,
    /// graph. The denominator of the incremental-against-rebuild ratio.
    build_ms: f64,
}

/// `Air` over exactly this bitmap.
///
/// `Air::build` reads nothing from a value but `is_inside`, so ±1 reproduces the
/// bitmap's notion of air **exactly** and there is no second sampling of the
/// field that could disagree with the first. It also means the vacuity control
/// and the octree cannot drift apart, which matters more than the memory it
/// saves at 257³.
fn air_of(signs: &Signs, sdim: u32) -> Air {
    let count = (sdim as usize).pow(3);
    let mut values = vec![0f32; count];
    for z in 0..sdim {
        for y in 0..sdim {
            for x in 0..sdim {
                let i = ((z as usize * sdim as usize) + y as usize) * sdim as usize + x as usize;
                values[i] = if signs.solid([x, y, z]) { -1.0 } else { 1.0 };
            }
        }
    }
    let shape = isomesh::RuntimeShape3::new([sdim; 3]).expect("shape");
    let (air, _) = Air::build(&values, &shape).expect("air");
    air
}

/// P-72's dig path: straight across `x` through the middle of the world, with
/// the height probed **per edit at that edit's own `x`**.
///
/// P-72 recorded two void runs before arriving at this shape — a path through
/// the world centre missed `fbm_terrain`'s sheet, and a path at one probed
/// height missed `gyroid`'s surface everywhere but the probe's own `x` — and
/// both were caught by the `M-44` control refusing to time a trace that marked
/// nothing. Inherited rather than reinvented.
fn p72_path<F: Sdf<Scalar = f64>>(field: &F) -> Vec<[f64; 3]> {
    let mid = ORIGIN + EXTENT * 0.5;
    let surface_y = |x: f64| -> f64 {
        let steps = 1024;
        let mut prev = field.sample([x, ORIGIN, mid]);
        for i in 1..=steps {
            let y = ORIGIN + EXTENT * (f64::from(i) / f64::from(steps));
            let v = field.sample([x, y, mid]);
            if (prev < 0.0) != (v < 0.0) {
                return y;
            }
            prev = v;
        }
        panic!("no surface crossing along y at x = {x}: the trace would dig in empty space");
    };
    (0..EDITS)
        .map(|i| {
            let t = (i as f64 + 0.5) / EDITS as f64;
            let x = ORIGIN + EXTENT * t;
            [x, surface_y(x), mid]
        })
        .collect()
}

/// A sealed chamber, then a shaft that breaks into it.
///
/// # Why this arm exists, and it is the most important thing in this harness
///
/// The registered vacuity control is *"the dig trace must change the free-space
/// topology at least once, asserted from the union-find"*. Read strictly — as a
/// change in the number of air components, which is the only topological
/// quantity `Air` maintains — **P-72's trace cannot satisfy it**, and the reason
/// is structural rather than accidental: every brush in that trace straddles the
/// surface, so every newly-air sample is adjacent to the open air that is
/// already there. The dig *widens* one component eleven times. `Air::dig` does
/// report `merges > 0` on those edits, but that is the blob's own fresh label
/// being absorbed within the same call, not two pre-existing components becoming
/// one, and scoring it as a topology change would be answering an easier
/// question (`P-70`'s C3).
///
/// So this arm is added, on `P-62`'s precedent — that row added a
/// 400,000-cell random arm because eight reference fields gave only seven tunnel
/// cells, *"a hair from `M-44`'s vacuous zero"*. One brush entirely inside solid
/// creates a **sealed pocket** (components +1); ten shaft brushes dug top-down
/// from the open air break into it (components −1). Both events are the game
/// event that matters — *"did I just break through"* — and both are read from
/// `Air::components()`.
///
/// # The fixture is searched, not assumed
///
/// A sealed chamber needs a sphere of radius `r` with a solid one-cell shell
/// around it, and `gyroid`'s solid labyrinth branch is about 0.65 world units
/// half-thick, so a 6-cell brush is 0.75 at 32³ and **does not fit**. The brush
/// is therefore the largest of 6, 4, 3 and 2 cells for which a chamber exists,
/// found by scanning the bitmap; among candidates the one closest to air wins,
/// so the shaft is short. Nothing is defaulted: if no chamber fits at any of the
/// four radii the harness panics rather than run a fixture that cannot show the
/// effect.
fn breakthrough_path(
    signs: &Signs,
    layout: &ChunkLayout<f64>,
    n: u32,
) -> (Vec<[f64; 3]>, f64, u32) {
    let cell = EXTENT / f64::from(n);
    for &brush in &[6u32, 4, 3, 2] {
        // The sphere plus a one-cell solid shell: what has to be solid for the
        // chamber to be sealed after the brush has run.
        let need = brush + 1;
        let reach = i64::from(need);
        let mut best: Option<([u32; 3], u32)> = None;
        let mut cz = need;
        while cz + need <= n {
            let mut cy = need;
            while cy + need <= n {
                let mut cx = need;
                while cx + need <= n {
                    let c = [cx, cy, cz];
                    if signs.solid(c) && sphere_is_solid(signs, c, reach) {
                        // Distance up to the nearest air sample.
                        let mut gap = 0u32;
                        let mut y = cy;
                        while y < n && signs.solid([cx, y, cz]) {
                            y += 1;
                            gap += 1;
                        }
                        if y < n && best.is_none_or(|(_, g)| gap < g) {
                            best = Some((c, gap));
                        }
                    }
                    cx += 2;
                }
                cy += 2;
            }
            cz += 2;
        }
        let Some((c, gap)) = best else { continue };
        let radius = f64::from(brush) * cell;
        let gap_world = f64::from(gap) * cell;
        // The shaft has to be a connected chain: consecutive brushes are `step`
        // apart and overlap only while `step < 2 * radius`.
        let total = gap_world + radius;
        let step = total / (EDITS - 1) as f64;
        assert!(
            step < 2.0 * radius,
            "shaft step {step} is not under one brush diameter {}: the {} shaft brushes would \
             not form a connected chain and the breakthrough would never happen",
            2.0 * radius,
            EDITS - 1
        );
        let centre = layout.world_of_sample([i64::from(c[0]), i64::from(c[1]), i64::from(c[2])]);
        let mut path = Vec::with_capacity(EDITS);
        path.push(centre);
        for k in (1..EDITS).rev() {
            path.push([centre[0], centre[1] + k as f64 * step, centre[2]]);
        }
        assert_eq!(
            path.len(),
            EDITS,
            "the breakthrough trace must have {EDITS} edits"
        );
        return (path, radius, brush);
    }
    panic!(
        "no sealed chamber of radius 6, 4, 3 or 2 cells exists at {n}³: the vacuity control \
         cannot be made to fire on this field and world size, and a fixture that cannot show \
         the effect must not be run"
    );
}

/// Every sample within `reach` cells of `c` is solid.
fn sphere_is_solid(signs: &Signs, c: [u32; 3], reach: i64) -> bool {
    let r2 = reach * reach;
    for dz in -reach..=reach {
        for dy in -reach..=reach {
            for dx in -reach..=reach {
                if dx * dx + dy * dy + dz * dz > r2 {
                    continue;
                }
                let s = [
                    i64::from(c[0]) + dx,
                    i64::from(c[1]) + dy,
                    i64::from(c[2]) + dz,
                ];
                for a in 0..3 {
                    if s[a] < 0 || s[a] >= i64::from(signs.dims[a]) {
                        return false;
                    }
                }
                if !signs.solid([s[0] as u32, s[1] as u32, s[2] as u32]) {
                    return false;
                }
            }
        }
    }
    true
}

/// Run one eleven-edit trace, `REPS` times, on a world already sampled.
#[allow(clippy::too_many_arguments)]
fn run_trace<F: Sdf<Scalar = f64>>(
    field: &F,
    name: &'static str,
    trace: &'static str,
    n: u32,
    layout: &ChunkLayout<f64>,
    signs: &mut Signs,
    pristine: &[u64],
    nav: &mut Nav,
    centres: &[[f64; 3]],
    radius: f64,
    brush_cells: f64,
) -> Arm {
    let sdim = n + 1;
    let mut stat = Static {
        all_leaves: 0,
        free_leaves: 0,
        regions: 0,
        edges: 0,
        components: 0,
        free_cells: 0,
    };
    let mut graph_components_after = 0u64;
    let mut post_regions_incremental = 0u64;
    let mut topology_changes = 0u64;
    let mut edits_joining = 0u64;
    let mut edits_pocketing = 0u64;
    let mut air_merges = 0u64;
    let mut air_newly_air = 0u64;
    let mut air_components_before = 0u64;
    let mut air_components_after = 0u64;
    let mut air_components_max = 0u64;
    let mut mesh_chunks = 0u64;
    let mut dirty_cells = 0u64;

    let mut rep_costs: Vec<Vec<Cost>> = Vec::with_capacity(REPS);
    let mut rep_times: Vec<Vec<(u128, u128)>> = Vec::with_capacity(REPS);
    let mut flipped: Vec<[u32; 3]> = Vec::new();
    let mut newly_air: Vec<[u32; 3]> = Vec::new();
    let mut seeds: Vec<u32> = Vec::new();
    // The from-scratch cost of the whole structure, so "incremental against
    // rebuild" is a ratio taken inside one build and one run (`M-281`) rather
    // than against a number from another machine's paper.
    let mut build_reps: Vec<f64> = Vec::with_capacity(REPS);

    for rep in 0..REPS {
        signs.words.copy_from_slice(pristine);
        let t = Instant::now();
        nav.tree.rebuild(signs);
        nav.build();
        build_reps.push(t.elapsed().as_nanos() as f64 / 1e6);
        let mut air = None;
        if rep == 0 {
            nav.audit("after build");
            stat = nav.snapshot();
            let a = air_of(signs, sdim);
            air_components_before = a.components();
            air_components_max = air_components_before;
            air = Some(a);
        }

        let mut costs = Vec::with_capacity(EDITS);
        let mut times = Vec::with_capacity(EDITS);
        let mut dirty = DirtySet::new();

        for step in 0..EDITS {
            let before = Dug {
                field,
                centres: &centres[..step],
                radius,
            };
            let after = Dug {
                field,
                centres: &centres[..=step],
                radius,
            };
            let c = centres[step];

            // The cell box the brush can touch, through `layout.cell_of` and not
            // by hand — P-72's second self-caught defect was exactly the
            // hand-rolled form, which assumes a zero origin. One cell of margin,
            // then clamped to the world.
            let lo_world = [0, 1, 2].map(|a| c[a] - radius);
            let hi_world = [0, 1, 2].map(|a| c[a] + radius);
            let lo_i = layout
                .cell_of(lo_world)
                .map(|v| (v - 1).clamp(0, i64::from(n) - 1));
            let hi_i = layout
                .cell_of(hi_world)
                .map(|v| (v + 1).clamp(0, i64::from(n) - 1));

            // Untimed: the crate's own instrument on the same box, as the
            // control on the harness's dirty set.
            let report = mark_edit(layout, &before, &after, lo_i, hi_i, &mut dirty).expect("mark");
            dirty.clear();

            // Untimed: the field evaluations. A chunked world pays these to
            // re-mesh; they are not a cost this stage adds.
            let lo = [0, 1, 2].map(|a| lo_i[a] as u32);
            let ext = [0, 1, 2].map(|a| (hi_i[a] - lo_i[a] + 2) as u32);
            let mut new_solid = vec![false; (ext[0] * ext[1] * ext[2]) as usize];
            for z in 0..ext[2] {
                for y in 0..ext[1] {
                    for x in 0..ext[0] {
                        let s = [
                            i64::from(lo[0] + x),
                            i64::from(lo[1] + y),
                            i64::from(lo[2] + z),
                        ];
                        let v = after.sample(layout.world_of_sample(s));
                        new_solid[((z * ext[1] + y) * ext[0] + x) as usize] = v < 0.0;
                    }
                }
            }

            let mut cost = Cost::default();
            let t = repair(nav, signs, lo, ext, &new_solid, &mut flipped, &mut cost);

            // ── control: two independent walks over the same edit ───────────
            assert_eq!(
                cost.pattern_cells, report.sign_changed_cells,
                "{name}/{trace} at {n}³, edit {step}: the harness found {} cells incident to a \
                 flipped sign bit and `mark_edit` found {} sign-changed cells over the same box, \
                 so one of the two is looking at the wrong region",
                cost.pattern_cells, report.sign_changed_cells
            );

            if let Some(a) = air.as_mut() {
                mesh_chunks += report.dirty_chunks;
                dirty_cells += report.sign_changed_cells;
                nav.audit("after edit");

                // ── the registered vacuity control ─────────────────────────
                //
                // `Air::components()` before and after, and nothing else: the
                // exact number of air components, maintained by crate source
                // that knows nothing about the octree, the merge or the brush.
                // `edits_joining` and `edits_pocketing` are reported beside it
                // as the mechanism, read from the pre-existing labels around the
                // new air, but the verdict column is the component delta because
                // that is the one quantity with no false positives.
                newly_air.clear();
                for &s in &flipped {
                    if !signs.solid(s) {
                        newly_air.push(s);
                    }
                }
                air_newly_air += newly_air.len() as u64;
                seeds.clear();
                for &s in &newly_air {
                    for axis in 0..3 {
                        for d in [-1i64, 1] {
                            let mut nb = [i64::from(s[0]), i64::from(s[1]), i64::from(s[2])];
                            nb[axis] += d;
                            if nb.iter().any(|&v| v < 0 || v >= i64::from(sdim)) {
                                continue;
                            }
                            let p = [nb[0] as u32, nb[1] as u32, nb[2] as u32];
                            if let Some(l) = a.label_of(p) {
                                seeds.push(l);
                            }
                        }
                    }
                }
                seeds.sort_unstable();
                seeds.dedup();
                if seeds.len() >= 2 {
                    edits_joining += 1;
                } else if seeds.is_empty() && !newly_air.is_empty() {
                    edits_pocketing += 1;
                }

                let comps_before = a.components();
                let counts = a.dig(&newly_air, || true);
                let comps_after = a.components();
                air_merges += counts.merges;
                air_components_max = air_components_max.max(comps_after);
                if comps_after != comps_before {
                    topology_changes += 1;
                }
            }

            costs.push(cost);
            times.push(t);
        }

        if let Some(a) = air {
            air_components_after = a.components();
            graph_components_after = nav.components();
            post_regions_incremental = nav.live_regions;
        }
        rep_costs.push(costs);
        rep_times.push(times);
    }

    // ── control: the integers are the same in every rep ─────────────────────
    for rep in 1..REPS {
        assert_eq!(
            rep_costs[rep], rep_costs[0],
            "{name}/{trace} at {n}³: rep {rep} produced different counts from rep 0 on an \
             identical trace, so the repair is not deterministic"
        );
    }

    // ── the local merge, scored against a global one on the same field ──────
    let mut fresh = Nav::new(n);
    fresh.tree.rebuild(signs);
    fresh.build();
    fresh.audit("global re-merge");
    let post_regions_global = fresh.live_regions;
    drop(fresh);

    // ── the timing row: the worst edit of the median rep ────────────────────
    //
    // The worst edit rather than the mean, because `M-124` is the whole reason
    // repair in this crate is budgeted: amortised is not the statistic for the
    // frame the breakthrough lands on. Both components come from that one edit,
    // so `octree_update_ms + graph_update_ms = total_ms` exactly.
    let mut per_rep: Vec<(f64, usize, f64, f64, f64, f64)> = Vec::with_capacity(REPS);
    for times in &rep_times {
        let mut totals: Vec<f64> = times.iter().map(|&(o, g)| (o + g) as f64 / 1e6).collect();
        let trace_total: f64 = totals.iter().sum();
        let (worst, &worst_ms) = totals
            .iter()
            .enumerate()
            .max_by(|a, b| a.1.partial_cmp(b.1).expect("finite"))
            .expect("edits");
        let (o, g) = times[worst];
        totals.sort_unstable_by(|a, b| a.partial_cmp(b).expect("finite"));
        per_rep.push((
            worst_ms,
            worst,
            o as f64 / 1e6,
            g as f64 / 1e6,
            totals[EDITS / 2],
            trace_total,
        ));
    }
    per_rep.sort_unstable_by(|a, b| a.0.partial_cmp(&b.0).expect("finite"));
    let (total_ms, worst_edit, octree_ms, graph_ms, median_edit_ms, trace_total_ms) =
        per_rep[REPS / 2];

    build_reps.sort_unstable_by(|a, b| a.partial_cmp(b).expect("finite"));
    let build_ms = build_reps[REPS / 2];

    let sum = |f: fn(&Cost) -> u64| rep_costs[0].iter().map(f).sum::<u64>();
    Arm {
        field: name,
        trace,
        world_cells: n,
        depth: nav.tree.depth,
        brush_cells,
        stat,
        dirty_cells,
        state_cells: sum(|c| c.state_cells),
        pattern_cells: sum(|c| c.pattern_cells),
        nav_chunks: sum(|c| c.nav_chunks),
        mesh_chunks,
        reclassified: sum(|c| c.reclassified),
        nodes: sum(|c| c.nodes),
        leaves_removed: sum(|c| c.leaves_removed),
        leaves_added: sum(|c| c.leaves_added),
        regions_removed: sum(|c| c.regions_removed),
        regions_added: sum(|c| c.regions_added),
        max_pool: rep_costs[0].iter().map(|c| c.pool).max().unwrap_or(0),
        repair_cells: rep_costs[0].iter().map(Cost::repair_cells).sum(),
        worst_edit,
        octree_ms,
        graph_ms,
        total_ms,
        median_edit_ms,
        trace_total_ms,
        graph_components_after,
        post_regions_incremental,
        post_regions_global,
        topology_changes,
        edits_joining,
        edits_pocketing,
        air_merges,
        air_newly_air,
        air_components_before,
        air_components_after,
        air_components_max,
        build_ms,
    }
}

/// Sample one field at one world size, then run both traces over it.
fn run_arm<F: Sdf<Scalar = f64>>(field: &F, name: &'static str, n: u32) -> [Arm; 2] {
    let cell_size = EXTENT / f64::from(n);
    let layout = ChunkLayout::<f64>::new(CHUNK_CELLS, cell_size, [ORIGIN; 3]).expect("layout");

    let sdim = n + 1;
    let mut signs = Signs::new([sdim; 3]);
    for z in 0..sdim {
        for y in 0..sdim {
            for x in 0..sdim {
                let p = layout.world_of_sample([i64::from(x), i64::from(y), i64::from(z)]);
                signs.set([x, y, z], field.sample(p) < 0.0);
            }
        }
    }
    let pristine = signs.words.clone();

    let dig = p72_path(field);
    let (shaft, shaft_radius, shaft_cells) = breakthrough_path(&signs, &layout, n);

    let mut nav = Nav::new(n);
    let a = run_trace(
        field,
        name,
        "p72_dig",
        n,
        &layout,
        &mut signs,
        &pristine,
        &mut nav,
        &dig,
        BRUSH_CELLS * cell_size,
        BRUSH_CELLS,
    );
    let b = run_trace(
        field,
        name,
        "breakthrough",
        n,
        &layout,
        &mut signs,
        &pristine,
        &mut nav,
        &shaft,
        shaft_radius,
        f64::from(shaft_cells),
    );
    assert_eq!(
        a.stat, b.stat,
        "{name} at {n}³: the two traces disagree about the world they started from, so the \
         rebuild is not a pure function of the bitmap"
    );
    [a, b]
}

/// The clock, on the row. `M-280`: on a governed CPU a nanosecond is not a unit.
fn clock_mhz() -> f64 {
    let text = std::fs::read_to_string("/proc/cpuinfo").unwrap_or_default();
    let mut best = 0f64;
    for line in text.lines() {
        let Some(rest) = line.strip_prefix("cpu MHz") else {
            continue;
        };
        let Some(value) = rest.split(':').nth(1) else {
            continue;
        };
        let Ok(mhz) = value.trim().parse::<f64>() else {
            continue;
        };
        if mhz > best {
            best = mhz;
        }
    }
    best
}

type Row = Vec<(&'static str, String)>;

fn main() {
    if !std::env::args().any(|arg| arg == "--bench") {
        return;
    }

    let prereg = isomesh::experiment!("P-87");
    let clock = clock_mhz();

    println!(
        "SHARE: the registration's share is 1.0 — the navigation-rebuild stage does not exist \
         yet. None of C1, C2, C3 is a speedup against a baseline, so ✗51's \
         1/(1 - s + s/factor) ceiling does not bind: C1 is an absolute 2 ms bar, C2 and C3 are \
         ratios of directly counted sets. All three were reachable before this run.\n"
    );
    println!(
        "{:>12} {:>13} {:>6} {:>6} {:>7} {:>8} {:>7} {:>8} {:>8} {:>8} {:>9} {:>8} {:>7} {:>5}",
        "field",
        "trace",
        "cells",
        "brush",
        "dirty",
        "repair",
        "factor",
        "octree",
        "graph",
        "total",
        "before",
        "after",
        "reduce",
        "topo"
    );

    let mut arms: Vec<Arm> = Vec::new();
    for name in ["fbm_terrain", "gyroid"] {
        for n in WORLD_SIZES {
            let pair = match name {
                "gyroid" => run_arm(&Gyroid::<f64>::canonical(), name, n),
                _ => run_arm(&FbmTerrain::<f64>::canonical(), name, n),
            };
            for arm in pair {
                println!(
                    "{:>12} {:>13} {:>6} {:>6.0} {:>7} {:>8} {:>7.3} {:>8.4} {:>8.4} {:>8.4} \
                     {:>9} {:>8} {:>7.3} {:>5}",
                    arm.field,
                    arm.trace,
                    arm.world_cells,
                    arm.brush_cells,
                    arm.dirty_cells,
                    arm.repair_cells,
                    arm.repair_cells as f64 / arm.dirty_cells as f64,
                    arm.octree_ms,
                    arm.graph_ms,
                    arm.total_ms,
                    arm.stat.free_leaves,
                    arm.stat.regions,
                    arm.stat.free_leaves as f64 / arm.stat.regions as f64,
                    arm.topology_changes
                );
                arms.push(arm);
            }
        }
    }

    // ── controls, asserted ──────────────────────────────────────────────────
    for a in &arms {
        assert!(
            a.dirty_cells > 0,
            "VOID: {}/{} at {}³ marked no dirty cell in {EDITS} edits",
            a.field,
            a.trace,
            a.world_cells
        );
        assert!(
            a.stat.free_leaves > 0 && a.stat.regions > 0,
            "VOID: {}/{} at {}³ has no free navigation cells at all",
            a.field,
            a.trace,
            a.world_cells
        );
    }

    // ── the registered vacuity control ──────────────────────────────────────
    //
    // Asserted on the breakthrough arm, which is the arm built so that it can
    // fire: one brush inside solid makes a sealed pocket and the shaft breaks
    // into it, so `Air::components()` must move at least twice. The P-72 arm's
    // value is reported rather than asserted, because a surface dig widens one
    // component and never joins two — see `breakthrough_path`'s docs. Scoring
    // that arm's `merges > 0` as a topology change would be `P-70`'s C3.
    for a in arms.iter().filter(|a| a.trace == "breakthrough") {
        assert!(
            a.topology_changes >= 2,
            "VACUOUS: the breakthrough trace on {} at {}³ changed the air-component count {} \
             times, and it is built to do it twice — a sealed pocket then a break-in. `Air` \
             reports {} components before, {} after, {} at the peak, {} joining edits, {} \
             pocketing edits, {} merges over {} newly-air samples. The repair is timing a no-op \
             (M-44).",
            a.field,
            a.world_cells,
            a.topology_changes,
            a.air_components_before,
            a.air_components_after,
            a.air_components_max,
            a.edits_joining,
            a.edits_pocketing,
            a.air_merges,
            a.air_newly_air
        );
    }

    // ── rows ────────────────────────────────────────────────────────────────
    let mut rows: Vec<Row> = Vec::new();
    for a in &arms {
        let reduction = a.stat.free_leaves as f64 / a.stat.regions as f64;
        let factor = a.repair_cells as f64 / a.dirty_cells as f64;
        // The other field at the same world size and trace, for C3's ordering
        // clause. The merge is a property of the static build, so both traces of
        // one arm carry the same reduction.
        let peer = arms
            .iter()
            .find(|b| b.world_cells == a.world_cells && b.trace == a.trace && b.field != a.field)
            .expect("both fields at every size");
        let peer_reduction = peer.stat.free_leaves as f64 / peer.stat.regions as f64;
        let gyroid_worst = if a.field == "gyroid" {
            reduction < peer_reduction
        } else {
            peer_reduction < reduction
        };
        let c3 = reduction >= 5.0 && peer_reduction >= 5.0 && gyroid_worst;
        // Growth against the smallest world at fixed edit size: world-proportional
        // is the registered falsifier, and 32³ → 256³ is 512× in cells.
        let base = arms
            .iter()
            .find(|b| b.field == a.field && b.trace == a.trace && b.world_cells == WORLD_SIZES[0])
            .expect("smallest world");
        let growth = a.repair_cells as f64 / base.repair_cells as f64;
        let c2 = factor < 3.0 && growth < 2.0;
        let c1 = a.total_ms < 2.0;

        rows.push(vec![
            ("field", a.field.to_string()),
            ("world_cells", a.world_cells.to_string()),
            ("chunk_cells", CHUNK_CELLS.to_string()),
            ("edits", EDITS.to_string()),
            ("dirty_cells", a.dirty_cells.to_string()),
            ("repair_cells", a.repair_cells.to_string()),
            ("repair_factor", format!("{factor:.4}")),
            ("octree_update_ms", format!("{:.6}", a.octree_ms)),
            ("graph_update_ms", format!("{:.6}", a.graph_ms)),
            ("total_ms", format!("{:.6}", a.total_ms)),
            ("cells_before_merge", a.stat.free_leaves.to_string()),
            ("cells_after_merge", a.stat.regions.to_string()),
            ("reduction", format!("{reduction:.4}")),
            ("topology_changes", a.topology_changes.to_string()),
            ("c1_holds", c1.to_string()),
            ("c2_holds", c2.to_string()),
            ("c3_holds", c3.to_string()),
            // ── extras ──────────────────────────────────────────────────────
            ("trace", a.trace.to_string()),
            ("brush_cells", format!("{:.0}", a.brush_cells)),
            ("clock_mhz", format!("{clock:.0}")),
            ("octree_depth", a.depth.to_string()),
            (
                "world_total_cells",
                u64::from(a.world_cells).pow(3).to_string(),
            ),
            ("repair_growth", format!("{growth:.4}")),
            (
                "repair_per_m311",
                format!("{:.4}", a.repair_cells as f64 / (EDITS as f64 * M311_DIRTY)),
            ),
            (
                "dirty_per_edit",
                format!("{:.1}", a.dirty_cells as f64 / EDITS as f64),
            ),
            ("state_changed_cells", a.state_cells.to_string()),
            ("pattern_changed_cells", a.pattern_cells.to_string()),
            ("nav_dirty_chunks", a.nav_chunks.to_string()),
            ("mesh_dirty_chunks", a.mesh_chunks.to_string()),
            ("reclassified_cells", a.reclassified.to_string()),
            ("nodes_reevaluated", a.nodes.to_string()),
            ("leaves_removed", a.leaves_removed.to_string()),
            ("leaves_added", a.leaves_added.to_string()),
            ("regions_removed", a.regions_removed.to_string()),
            ("regions_added", a.regions_added.to_string()),
            ("max_pool_leaves", a.max_pool.to_string()),
            ("all_leaves", a.stat.all_leaves.to_string()),
            ("free_cells", a.stat.free_cells.to_string()),
            (
                "leaves_per_kilo_free_cell",
                format!(
                    "{:.3}",
                    1000.0 * a.stat.free_leaves as f64 / a.stat.free_cells as f64
                ),
            ),
            ("graph_edges", a.stat.edges.to_string()),
            ("graph_components_before", a.stat.components.to_string()),
            (
                "graph_components_after",
                a.graph_components_after.to_string(),
            ),
            (
                "post_regions_incremental",
                a.post_regions_incremental.to_string(),
            ),
            ("post_regions_global", a.post_regions_global.to_string()),
            (
                "merge_drift",
                format!(
                    "{:.4}",
                    a.post_regions_incremental as f64 / a.post_regions_global as f64
                ),
            ),
            ("worst_edit", a.worst_edit.to_string()),
            ("median_edit_ms", format!("{:.6}", a.median_edit_ms)),
            ("trace_total_ms", format!("{:.6}", a.trace_total_ms)),
            ("build_ms", format!("{:.6}", a.build_ms)),
            (
                "rebuild_over_repair",
                format!("{:.2}", a.build_ms / a.total_ms),
            ),
            (
                "rebuild_over_trace",
                format!("{:.2}", a.build_ms / a.trace_total_ms),
            ),
            ("air_components_before", a.air_components_before.to_string()),
            ("air_components_after", a.air_components_after.to_string()),
            ("air_components_max", a.air_components_max.to_string()),
            ("edits_joining", a.edits_joining.to_string()),
            ("edits_pocketing", a.edits_pocketing.to_string()),
            ("air_merges", a.air_merges.to_string()),
            ("air_newly_air", a.air_newly_air.to_string()),
        ]);
    }

    common::experiment::run(prereg, |run| {
        for row in &rows {
            run.record(row);
        }
    });
}
