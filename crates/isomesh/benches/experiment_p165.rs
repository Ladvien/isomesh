//! **P-165 — how far greedy meshing is from the optimum, which is computable in
//! polynomial time.**
//!
//! Ticket: R-165. Pre-registered before this harness existed.
//!
//! ```bash
//! cargo bench --bench experiment_p165
//! ```
//!
//! Writes `docs/experiments/p-165.csv`.
//!
//! # What was missing
//!
//! `M-56` (`FINDINGS.md:1172`) measured greedy meshing's saving **over face
//! culling** — `Merge::Greedy` against `Merge::Off` on the same occupancy — and
//! found it is a property of the scene rather than of the algorithm: `gyroid`
//! **1.70×**, `sphere` **1.94×**, `torus` **2.69×**, `fbm_terrain` **4.60×**,
//! `csg_difference` **10.64×**, `box_exact` **256×**, all at 33³. Every one of
//! those is a ratio against a **worse** algorithm. Nobody has ever asked the
//! other direction: how far is greedy from the **best possible** rectangle
//! partition of the very same masks?
//!
//! That question has had an exact polynomial-time answer since 1979 and this
//! crate has never used it. For a rectilinear region the minimum number of
//! rectangles in a partition is
//!
//! ```text
//! K_min = n/2 + h − g − 1
//! ```
//!
//! where `n` is the **total** number of boundary vertices, `h` the number of
//! holes, and `g` the maximum number of pairwise-disjoint **good diagonals** —
//! axis-parallel segments interior to the region joining two concave vertices.
//! Eppstein, *Graph-Theoretic Solutions to Computational Geometry Problems*
//! (`arXiv:0908.3916`) states it in this form and supplies the reduction that
//! makes `g` cheap: the intersection graph of horizontal against vertical chords
//! is **bipartite**, so a maximum independent set of chords is
//! `(#chords) − (maximum matching)` by König–Gallai. Independently discovered by
//! Lipski, Lodi, Luccio, Mugnai & Pagli (1979), Ohtsuki (1982), and Ferrari,
//! Sankar & Sklansky (1984), which is why no single citation owns it.
//!
//! Two equivalent readings of the same formula are used below, because their
//! *disagreement* is this harness's main correctness gate. Writing `r` for the
//! reflex (concave) vertex count and using the discrete turning-number theorem
//! for a rectilinear region — `convex − reflex = 4 − 4h`, hence
//! `n = 2r + 4 − 4h` — the formula is also `K_min = r − h − g + 1`. The two are
//! computed from **separately measured** quantities (`n` and `h` from a boundary
//! walk, `r` from a right-turn count) and asserted equal per component.
//!
//! # The region, and why it is one binary slice rather than the volume
//!
//! `GreedyQuads` is not an isosurface extractor (`extractor.rs:243-256`). It
//! classifies whole cells by the sign at the **cell centre** and emits the
//! axis-aligned faces between a solid cell and an empty one
//! (`greedy_quads.rs:178-193`). Its merge then runs on a strictly 2D object: for
//! each `axis` in `0..3`, each `step` in `[1, −1]` and each `slice` in
//! `0..cells[axis]`, it builds the mask `at(c) && !at(n)` over a `du × dv`
//! grid and walks it — run along `u`, extend along `v`, emit one quad
//! (`greedy_quads.rs:214-272`).
//!
//! **So the rectilinear region is one such mask, and there are `6·(n−1)` of
//! them.** That is not an approximation of greedy's problem, it *is* greedy's
//! problem: the shipped mesher partitions each mask into rectangles
//! independently and can never merge across two masks, because two masks are two
//! different planes with two different normals — and a shared vertex there would
//! have to average two normals, which the module doc rejects at
//! `greedy_quads.rs:36-49`. The optimum is therefore computed per mask and
//! summed, and the comparison is exact rather than analogical.
//!
//! A **whole-volume** claim from this number would be a category error and is not
//! made anywhere here. "The minimum number of axis-aligned quads needed to
//! represent this occupancy grid's surface" is a different and harder problem
//! (it would allow a quad to span what greedy calls two sweeps); what is measured
//! is "the minimum number of rectangles in a partition of each of the `6·(n−1)`
//! masks greedy actually partitions".
//!
//! Two exactnesses make the per-mask sum legitimate:
//!
//! - **Rectangle partition decomposes over 4-connected components**, with no
//!   theorem required: every rectangle of set cells is 4-connected, so it lies
//!   inside one component. `K_min(mask) = Σ K_min(component)` identically.
//! - **`g` decomposes too.** A chord's interior is interior to the region, so all
//!   its cells are set and belong to one component; a horizontal chord in one
//!   component cannot meet a vertical chord in another. The matching is therefore
//!   run once per mask and `Σ_c g_c = g_mask` exactly, which is why the CSV
//!   carries a single `good_diagonals_g` per row and a `components` count to
//!   recover the arithmetic: `optimum_rectangles = vertices_n/2 + holes_h −
//!   good_diagonals_g − components`.
//!
//! # Degeneracy, which is what C1's falsifier is about
//!
//! Eppstein's theorem is about a polygon. A binary grid mask is not always one:
//! at a **checkerboard vertex** — two diagonally opposite cells set, the other
//! two clear — the boundary passes through one point twice and the region is not
//! locally a disc. `pinch_vertices` counts them and they are common (1,176 on
//! `gyroid` at 33³ in the prototype that sized this harness). The falsifier
//! *"C1 by the reduction not applying — for instance if the merged regions have
//! holes the formula's `h` does not model"* is exactly this hazard, so it is
//! resolved by construction and then **checked against ground truth** rather than
//! assumed:
//!
//! 1. **4-connected components first.** A checkerboard vertex whose two solid
//!    quadrants are in different components is fully handled by this step alone:
//!    two unit cells touching only at a corner are two regions, `1 + 1 = 2`
//!    rectangles, and no formula is asked to see them as one.
//! 2. **The boundary is walked with the solid on the left and the leftmost turn
//!    taken at every vertex.** At an ordinary vertex there is exactly one
//!    outgoing edge and the rule is vacuous. At a checkerboard vertex there are
//!    two in and two out, and the leftmost rule pairs each solid quadrant's two
//!    rays with each other — the *cell's own corner* — so the pinch contributes
//!    **two convex vertices** and the walk never crosses into the other quadrant.
//! 3. **`h` is `cycles − 1` per component**, from that walk. This is the step
//!    that a naive enclosure count gets wrong. A 4-connected component with a
//!    pinch always encloses at least one clear cell (the path joining the two
//!    solid quadrants plus the pinch is a closed curve, and the two clear
//!    quadrants lie on opposite sides of it), yet the walk merges the outer curve
//!    and that enclosure into **one** cycle, giving `h = 0`. The prototype's
//!    seven-cell pinched ring is the smallest witness: `n = 10`, `r = 3`,
//!    `h = 0`, `g = 0`, so `10/2 + 0 − 0 − 1 = 4`, and 4 is the true minimum.
//!    Counting the enclosure as a hole gives 3 and is wrong.
//! 4. **The turning identity is asserted per component**, in the form that bites:
//!    `n/2 + h − 1` must equal `r − h + 1`. `identity_failures` counts violations
//!    and a non-zero count falsifies C1 for that row. It is not a restatement —
//!    `n` and `h` come from the walk, `r` from the right-turn tally.
//! 5. **The formula is compared against exhaustive search** on every component of
//!    at most `BRUTE_LIMIT` cells, by the memoised branch-and-bound in
//!    [`brute_min_rectangles`], which knows nothing about Eppstein.
//!    `brute_disagreements` is the column that would end this row's C1.
//!
//! # Why the chord graph is bipartite, which the harness relies on and does not
//! merely cite
//!
//! A reflex grid vertex has exactly **three** of its four quadrant cells solid.
//! The single clear quadrant blocks one horizontal direction and one vertical
//! direction: with `NE` clear, the segment leaving west has solid above and below
//! and is interior, while the one leaving east has the clear `NE` cell on one
//! side and is boundary. Three consequences, all used:
//!
//! - **No axis-parallel interior segment passes *through* a reflex vertex**, so a
//!   good diagonal has a reflex vertex only at its two ends and no chord contains
//!   a third concave vertex.
//! - **No two horizontal chords share an endpoint or overlap**, because a reflex
//!   vertex admits an interior horizontal segment on exactly one side. Hence
//!   there are no H–H edges, and by symmetry no V–V edges: the intersection graph
//!   really is bipartite, and König's theorem really does apply.
//! - **A pinch vertex is never a chord endpoint** — both its horizontal sides and
//!   both its vertical sides are boundary — so the degeneracy of §"Degeneracy"
//!   does not leak into `g`.
//!
//! Chords are therefore enumerated in one linear scan per grid line: on line `y`
//! the predicate `solid(k, y−1) && solid(k, y)` marks the cells an interior
//! horizontal segment can cross, and each maximal run of it yields **at most one**
//! chord, namely when both of its endpoints are reflex with the matching
//! direction. Two axis-parallel segments on grid lines can only meet at a grid
//! point, so adjacency is built by marking each horizontal chord's grid points
//! and walking each vertical chord's — `O(cells)`, no `O(H·V)` product.
//!
//! The matching is **Hopcroft–Karp**, written in this file: BFS layering from the
//! free horizontal chords, then layered DFS augmentation, iterated until a phase
//! augments nothing (Berge's lemma is the termination condition, the layering is
//! the `O(E√V)`).
//!
//! # How `greedy_rectangles` was verified, before being relied on
//!
//! `emit_quad` (`greedy_quads.rs:293-341`) pushes exactly **four vertices and two
//! triangles** per merged quad, so `triangle_count() / 2` is the quad count and
//! `vertex_count() / 4` is the same number. That is read off the source; it is
//! *checked* three independent ways, and the row carries the evidence:
//!
//! 1. **The shipped mesh's own arithmetic.** `vertex_count() == 2 ·
//!    triangle_count()` and `triangle_count()` even, asserted on every extraction.
//!    Recorded as `greedy_vertices` and `greedy_triangles`.
//! 2. **This file's independent replication.** [`greedy_quad_count`] re-walks the
//!    same masks with the same merge and is asserted equal to
//!    `triangle_count() / 2`. This is the load-bearing one: it is what proves the
//!    masks analysed here **are** the masks the shipped extractor partitioned,
//!    down to the occupancy sample points. `mask_quads_agree` records it.
//! 3. **The committed golden fixture.** `crates/isomesh/golden_hashes.json`
//!    carries `greedy_quads` rows at 17, 25 and 33 samples for all eight fields;
//!    `golden_check` reads the `triangles` field at 17 and 33 and asserts equality,
//!    and reads `absent` at 65 where the fixture has no row. This is the only
//!    check that is independent of *this run*.
//!
//! The replication was additionally validated against `M-56`'s published ratios
//! before this file was written, in a prototype using the same occupancy rule:
//! `box_exact` 1536/6 = **256×**, `csg_difference` 1596/150 = **10.64×**,
//! `torus` 984/366 = **2.69×**, `sphere` 1248/642 = **1.94×**, `gyroid`
//! 5250/3093 = **1.70×** at 33³. All five reproduce `M-56` exactly, so
//! `culled_faces / greedy_rectangles` is recorded as `merge_saving` and this row
//! re-derives `M-56` in the same line that measures the distance to the optimum.
//!
//! # Seven fields or eight
//!
//! The registration says *"the optimum is computed on all **seven** fields"*.
//! `for_each_reference_field!` has **eight** (`fields/mod.rs:211-256`). The
//! registration cannot be amended, so the eighth field is surplus rather than a
//! violation and **all eight are run**.
//!
//! The discrepancy is not a slip: it is inherited from `M-56`, whose sweep is
//! *"seven reference fields at 33³"* — and the missing one is `thin_plate`,
//! because **`thin_plate` has no blocky mesh at all**. It is a slab of
//! half-thickness `0.0125` centred on `y = 0` (`fields/mod.rs:617-637`), the grid
//! is `n` samples over `[−2, 2]` with `n` odd, and a cell centre sits at
//! `−2 + h·j + h/2`, which equals `0` only for half-integer `j`. No cell centre
//! is ever inside the plate, at any odd resolution. The committed fixture agrees:
//! `greedy_quads / thin_plate` reads `"vertices":0,"triangles":0` at 17, 25 **and**
//! 33. So `thin_plate`'s three rows carry `greedy_rectangles = 0`,
//! `optimum_rectangles = 0`, `components = 0` and `ratio_defined = false` — the
//! eighth field is the one that shows why the registration says seven, which is a
//! better outcome than dropping it. It is M-266's sub-voxel aliasing seen through
//! a centre-sampled mesher.
//!
//! # Arms
//!
//! | arm | what it varies | is_control |
//! |---|---|---|
//! | `GreedyQuads` with `Merge::Greedy`, per field per resolution | the shipped algorithm — this is `greedy_rectangles` | no |
//! | Eppstein's `n/2 + h − g − 1` over the same masks | the optimum — `optimum_rectangles` | no |
//! | `GreedyQuads` with `Merge::Off` | face culling, `M-56`'s baseline — `culled_faces`, `merge_saving` | **yes** (re-derives `M-56`) |
//! | this file's mask replication | `mask_quads_agree` — proves the masks are greedy's masks | **yes** |
//! | `golden_hashes.json` at 17 and 33 | `golden_check` — proves the count is the blessed shipped one | **yes** |
//! | exhaustive search on components ≤ `BRUTE_LIMIT` cells | `brute_disagreements` — proves the formula, not just its citation | **yes** |
//! | the turning identity `n/2 + h − 1 == r − h + 1` | `identity_failures` — proves `n`, `h` and `r` are mutually consistent | **yes** |
//!
//! Three resolutions: **17, 33, 65** samples. 17 and 33 are the golden fixture's
//! own (`golden.rs:72` is `[17, 25, 33]`), which is what makes `golden_check`
//! possible on two of the three; 65 straddles the `u64` word boundary the repo
//! deliberately tests at and supplies a workload where the `O(cells)` claims
//! about the analysis can be seen in the clock. `resolution` counts **samples**,
//! so `n` samples span `n − 1` cells (`benches/common/mod.rs:40-43`).
//!
//! # SHARE, recomputed before the numbers
//!
//! The registration's SHARE is *"C2 moves the greedy-meshing stage, whose saving
//! `M-56` bounds."* Discharged, and narrowed: what C2 can move is the
//! **expectation**, not the code. `M-56` left greedy's quality bracketed only from
//! below (`1.70×`–`256×` better than culling) and unbounded from above, so any
//! future work on the merge could claim an arbitrary headroom. `ratio` closes
//! that: it is the exact multiplicative distance from optimal on the very masks
//! the shipped mesher partitions. If it is at or near `1`, the greedy-meshing
//! stage is **finished** and the correct engineering action is to stop working on
//! it — which is a change to the plan for that stage, and the only kind of change
//! a measurement is entitled to make. This row proposes no source change; it
//! prices one.
//!
//! `optimum_ms` is the second half of the same sentence. If the optimum is cheap
//! it is a candidate; if it is not, C3 says so and the formula is labelled a
//! measurement instrument. Note the deliberate non-claim: this harness computes
//! the optimal **count**, not an optimal **partition**. A shippable optimal
//! mesher would additionally have to construct one, which is not measured here
//! and must not be read out of `optimum_ms`.
//!
//! # Vacuity controls
//!
//! - **At least one field must report `good_diagonals_g > 0`.** With `g == 0`
//!   everywhere the formula degenerates to `n/2 + h − 1`, the whole bipartite
//!   reduction — the part of Eppstein that is actually a contribution — is never
//!   exercised, and the Hopcroft–Karp code could be deleted without changing a
//!   number. Column: `good_diagonals_g`.
//! - **At least one field must report `holes_h > 0`.** `h` is the term C1's
//!   falsifier names by name. Asserted separately from `g`, because a fixture can
//!   easily have chords and no holes. Column: `holes_h`.
//! - **At least one field must report `pinch_vertices > 0`.** Otherwise every
//!   word of §"Degeneracy" is untested and the leftmost-turn rule is being
//!   asserted over a population where it cannot differ from the naive one.
//! - **At least one non-rectangular component must reach the exhaustive
//!   oracle.** A brute-force check that only ever saw rectangles would agree with
//!   any formula that returns 1 for a rectangle. Column:
//!   `brute_nonrectangular`.
//! - **`culled_faces` must exceed `greedy_rectangles` on at least one field**, or
//!   `merge_saving` is `1` everywhere and the control is not reproducing `M-56`.
//!
//! Two further asserts are not vacuity controls but hard stops, because a
//! violation means the numbers describe something other than the shipped
//! algorithm: `mask_quads_agree` (the replication) and `golden_check` (the
//! fixture). Both abort rather than record.
//!
//! # Verdicts
//!
//! - `c1_holds` is **per row**: `identity_failures == 0`,
//!   `brute_disagreements == 0`, and `optimum_rectangles <= greedy_rectangles`
//!   (greedy exhibits a partition, so an optimum above it is arithmetically
//!   impossible and would mean the reduction is misapplied). On `thin_plate` all
//!   three hold over an empty population; `components = 0` on the row is what
//!   makes that visible rather than hidden.
//! - `c2_holds` is **global** and carries the same value on every row, because
//!   the registration's falsifier is global: *"C2 by greedy already achieving the
//!   optimum on every field"*. So C2 holds iff `ratio > 1` on at least one row
//!   with a defined ratio. The per-row reading is not lost — it is `ratio`
//!   itself, and `ratio_defined` marks the rows that cannot speak.
//! - `c3_holds` is **per row**: `cost_ratio = optimum_ms / greedy_ms < 10`.
//!   `greedy_ms` times the shipped `extract`, which includes the `(n−1)³` SDF
//!   samples it cannot avoid; `optimum_ms` times occupancy → optimum, which
//!   re-uses those samples and pays for none of them. The ratio is therefore
//!   measured in the direction that makes C3 **easiest** to pass, and
//!   `total_cost_ratio = (greedy_ms + optimum_ms) / greedy_ms` is recorded so the
//!   other reading is on the artefact too. `REPEATS` repeats, median as the
//!   headline, min and max recorded; `cost_ratio_worst =
//!   optimum_max_ms / greedy_min_ms` is the least favourable pairing, so a row
//!   whose scatter crosses the threshold while its median does not is reported
//!   with that scatter rather than averaged into a pass (M-280: this host's
//!   `amd-pstate-epp` governor swings the same binary 1.45×).

#![allow(
    clippy::cast_possible_truncation,
    clippy::cast_precision_loss,
    clippy::cast_sign_loss,
    clippy::too_many_lines
)]

mod common;

use std::time::Instant;

use isomesh::greedy_quads::{GreedyQuads, Merge};
use isomesh::marching_cubes::table::is_inside;
use isomesh::{MeshBuffer, Sdf, Shape3, for_each_reference_field};

/// Samples per axis. 17 and 33 are the golden fixture's own resolutions, which
/// is what lets `golden_check` anchor two of the three rows per field against a
/// committed number; 65 straddles the `u64` word boundary and supplies a real
/// workload for the timing clause.
const RESOLUTIONS: [u32; 3] = [17, 33, 65];

/// Timed repeats per row. The contract's floor for a registration that names a
/// ratio threshold, which C3 does.
const REPEATS: usize = 5;

/// Largest component handed to the exhaustive oracle, in cells.
///
/// The memo is a dense table over the `2^BRUTE_LIMIT` subsets of one component,
/// so this constant is also the size of that table. Twelve keeps it at 4096
/// entries and still covers every topologically interesting small case: the
/// 3×3 square annulus is 8 cells, the smallest pinched ring is 7, and every
/// L, T, plus and staircase is smaller than both.
const BRUTE_LIMIT: usize = 12;

/// C3's threshold, from the registration: *"C3 by above 10x"*.
const COST_LIMIT: f64 = 10.0;

// ─── the occupancy grid, sampled exactly the way `GreedyQuads` samples it ────

/// Cell-centre occupancy, byte-for-byte the grid `greedy_quads.rs:178-193`
/// builds.
///
/// One sample per cell, at the centre, `is_inside` deciding — not eight corners.
/// The arithmetic is transcribed in the same order as the source
/// (`origin + cell_size · index + half`) because a reassociated sum is a
/// different `f64` and a single flipped cell would put this file's masks out of
/// step with the extractor's.
#[derive(Debug)]
struct Occupancy {
    /// Cells per axis, which is `samples − 1`.
    cells: usize,
    /// Row-major with `x` fastest, matching `Shape3`'s convention.
    solid: Vec<bool>,
}

impl Occupancy {
    /// Sample `sdf` at every cell centre of the grid `common::grid` describes.
    fn of<S: Sdf<Scalar = f64>>(sdf: &S, cells: usize, origin: [f64; 3], cell_size: f64) -> Self {
        let half = cell_size * 0.5;
        let mut solid = Vec::with_capacity(cells * cells * cells);
        for z in 0..cells {
            for y in 0..cells {
                for x in 0..cells {
                    let p = [
                        origin[0] + cell_size * (x as f64) + half,
                        origin[1] + cell_size * (y as f64) + half,
                        origin[2] + cell_size * (z as f64) + half,
                    ];
                    solid.push(is_inside(sdf.sample(p)));
                }
            }
        }
        Self { cells, solid }
    }

    /// Whether cell `c` is solid, with **outside the grid counting as empty** —
    /// the capping rule at `greedy_quads.rs:195-205`, which is why a solid cell
    /// at the domain edge emits a face there.
    fn at(&self, c: [isize; 3]) -> bool {
        let n = self.cells as isize;
        if c[0] < 0 || c[1] < 0 || c[2] < 0 || c[0] >= n || c[1] >= n || c[2] >= n {
            return false;
        }
        self.solid[c[0] as usize + self.cells * (c[1] as usize + self.cells * c[2] as usize)]
    }
}

// ─── one binary slice: the rectilinear region ────────────────────────────────

/// One sweep's visible-face mask — the 2D binary region greedy partitions into
/// quads and the region Eppstein's formula is applied to.
#[derive(Debug, Default)]
struct Mask {
    /// Extent along the sweep's `u` axis.
    du: usize,
    /// Extent along the sweep's `v` axis.
    dv: usize,
    /// `cell[a + du·b]`, set where a face exists.
    cell: Vec<bool>,
}

impl Mask {
    /// Whether mask cell `(i, j)` is set, with outside the mask counting as
    /// clear. Signed, because every boundary and chord predicate below reads one
    /// cell off each edge.
    fn at(&self, i: isize, j: isize) -> bool {
        if i < 0 || j < 0 || i >= self.du as isize || j >= self.dv as isize {
            return false;
        }
        self.cell[i as usize + self.du * j as usize]
    }

    /// The four cells around grid point `(px, py)`, as `[SW, SE, NW, NE]`.
    ///
    /// Every vertex predicate in this file is a function of these four bits and
    /// nothing else, which is what makes the boundary walk and the chord scan
    /// agree by construction rather than by care.
    fn quadrants(&self, px: isize, py: isize) -> [bool; 4] {
        [
            self.at(px - 1, py - 1),
            self.at(px, py - 1),
            self.at(px - 1, py),
            self.at(px, py),
        ]
    }
}

/// Build every mask in `GreedyQuads`' own sweep order and hand each to `body`.
///
/// `axis`, `u = (axis+1)%3`, `v = (axis+2)%3`, `step ∈ [1, −1]`,
/// `slice ∈ 0..cells[axis]` — `greedy_quads.rs:208-230`, in that order, so the
/// masks arrive in the same sequence the extractor consumed them.
fn for_each_mask(occ: &Occupancy, mask: &mut Mask, mut body: impl FnMut(&Mask)) {
    for axis in 0..3usize {
        let u = (axis + 1) % 3;
        let v = (axis + 2) % 3;
        // `common::grid` builds `[samples; 3]`, so every axis has the same cell
        // count and `cells[u] == cells[v] == cells[axis]`.
        mask.du = occ.cells;
        mask.dv = occ.cells;
        for step in [1isize, -1] {
            for slice in 0..occ.cells as isize {
                mask.cell.clear();
                mask.cell.resize(mask.du * mask.dv, false);
                for b in 0..mask.dv {
                    for a in 0..mask.du {
                        let mut c = [0isize; 3];
                        c[axis] = slice;
                        c[u] = a as isize;
                        c[v] = b as isize;
                        let mut n = c;
                        n[axis] = slice + step;
                        mask.cell[a + mask.du * b] = occ.at(c) && !occ.at(n);
                    }
                }
                body(mask);
            }
        }
    }
}

/// The shipped greedy merge, counting quads instead of emitting them.
///
/// `greedy_quads.rs:238-272` in structure and in loop order: walk the mask, take
/// the widest run along `u`, extend along `v` for as many rows as keep that full
/// width, clear what was consumed. The only difference is that `emit_quad` is
/// replaced by `quads += 1`.
///
/// This is a **control**, not a second implementation of anything shipped: its
/// output is asserted equal to `MeshBuffer::triangle_count() / 2`, and that
/// equality is what establishes that the masks analysed in this file are the
/// masks the extractor partitioned — same occupancy samples, same sweep order,
/// same visibility test. Without it, every number below would rest on the
/// assumption that this file reproduced `greedy_quads.rs` correctly.
fn greedy_quad_count(mask: &Mask, scratch: &mut Vec<bool>) -> u64 {
    scratch.clear();
    scratch.extend_from_slice(&mask.cell);
    let (du, dv) = (mask.du, mask.dv);
    let mut quads = 0u64;
    for b in 0..dv {
        let mut a = 0usize;
        while a < du {
            if !scratch[a + du * b] {
                a += 1;
                continue;
            }
            let mut width = 1usize;
            let mut height = 1usize;
            while a + width < du && scratch[a + width + du * b] {
                width += 1;
            }
            'grow: while b + height < dv {
                for k in 0..width {
                    if !scratch[a + k + du * (b + height)] {
                        break 'grow;
                    }
                }
                height += 1;
            }
            for row in 0..height {
                for k in 0..width {
                    scratch[a + k + du * (b + row)] = false;
                }
            }
            quads += 1;
            a += width;
        }
    }
    quads
}

// ─── Eppstein's three numbers ───────────────────────────────────────────────

/// Direction codes `0 = +u (E)`, `1 = +v (N)`, `2 = −u (W)`, `3 = −v (S)`, so
/// `(d + 1) % 4` is a left turn and `(d + 3) % 4` a right one.
const STEP: [[isize; 2]; 4] = [[1, 0], [0, 1], [-1, 0], [0, -1]];

/// The four 4-connected neighbour offsets.
const NEIGHBOURS: [[isize; 2]; 4] = [[-1, 0], [1, 0], [0, -1], [0, 1]];

/// Is there a boundary edge leaving grid point `q`'s vertex in direction `d`,
/// oriented with the **solid on its left**?
///
/// Each unit boundary segment separates one solid cell from one clear one, so it
/// admits exactly one such direction; these four predicates therefore *are* the
/// directed boundary graph, with in-degree equal to out-degree at every vertex.
fn leaves(q: [bool; 4], d: usize) -> bool {
    let [sw, se, nw, ne] = q;
    match d {
        0 => ne && !se,
        1 => nw && !ne,
        2 => sw && !nw,
        _ => se && !sw,
    }
}

/// The cell immediately left of the edge leaving `(px, py)` in direction `d`.
///
/// `leaves` guarantees that cell is solid, hence in range, hence labelled — so
/// this is how a cycle is attributed to a component.
fn left_cell(px: isize, py: isize, d: usize) -> [isize; 2] {
    match d {
        0 => [px, py],
        1 => [px - 1, py],
        2 => [px - 1, py - 1],
        _ => [px, py - 1],
    }
}

/// Which quadrant of a reflex grid vertex is clear — `0 = SW`, `1 = SE`,
/// `2 = NW`, `3 = NE` — or `None` if the vertex is not reflex.
///
/// A reflex vertex is exactly a grid point with three of four quadrants solid.
/// Two solid quadrants means flat (adjacent) or a pinch (diagonal), one means
/// convex, four means interior; none of those can end a good diagonal.
fn reflex_clear(q: [bool; 4]) -> Option<usize> {
    if q.iter().filter(|&&b| b).count() != 3 {
        return None;
    }
    q.iter().position(|&b| !b)
}

/// The one horizontal and one vertical direction along which an interior
/// segment can leave a reflex vertex whose clear quadrant is `clear`: `.0` is
/// true when the horizontal one runs east, `.1` when the vertical one runs
/// north.
///
/// With `NE` clear, for instance, the westward segment has solid above and below
/// and is interior while the eastward one has `NE` on one side and is boundary;
/// likewise southward is interior and northward is boundary.
fn chord_dirs(clear: usize) -> (bool, bool) {
    (clear == 0 || clear == 2, clear == 0 || clear == 1)
}

/// Eppstein's numbers for one mask, summed over its 4-connected components.
///
/// `rectangles` is `Σ_c (n_c/2 + h_c − 1) − g`, which is the formula's
/// `n/2 + h − g − 1` per component with the single mask-wide `g` subtracted once
/// — legitimate because the chord graph decomposes over components.
#[derive(Clone, Copy, Debug, Default)]
struct Optimum {
    /// Total boundary vertices, a pinch counted once per solid quadrant.
    vertices_n: i64,
    /// Reflex (concave) vertices, counted as right turns on the walk.
    reflex_r: i64,
    /// Holes, `cycles − 1` per component.
    holes_h: i64,
    /// Maximum pairwise-disjoint good diagonals.
    good_diagonals_g: i64,
    /// 4-connected components of set cells.
    components: i64,
    /// The optimum itself.
    rectangles: i64,
    /// Horizontal chords found.
    chords_h: i64,
    /// Vertical chords found.
    chords_v: i64,
    /// Maximum matching in the chord intersection graph.
    matched: i64,
    /// Checkerboard grid vertices.
    pinch_vertices: i64,
    /// Components whose `n/2 + h − 1` disagreed with `r − h + 1`.
    identity_failures: i64,
    /// Masks with at least one set cell.
    nonempty_masks: i64,
    /// Set mask cells, which is the face-culling quad count.
    faces: i64,
}

impl Optimum {
    /// Accumulate one mask's numbers into a running total.
    fn add(&mut self, other: &Self) {
        self.vertices_n += other.vertices_n;
        self.reflex_r += other.reflex_r;
        self.holes_h += other.holes_h;
        self.good_diagonals_g += other.good_diagonals_g;
        self.components += other.components;
        self.rectangles += other.rectangles;
        self.chords_h += other.chords_h;
        self.chords_v += other.chords_v;
        self.matched += other.matched;
        self.pinch_vertices += other.pinch_vertices;
        self.identity_failures += other.identity_failures;
        self.nonempty_masks += other.nonempty_masks;
        self.faces += other.faces;
    }
}

/// The mutable half of Hopcroft–Karp, held apart from the read-only CSR so the
/// two can be borrowed from the same [`Scratch`] at once.
#[derive(Debug, Default)]
struct MatchState {
    /// Partner of each horizontal chord, `−1` for free.
    match_h: Vec<i32>,
    /// Partner of each vertical chord, `−1` for free.
    match_v: Vec<i32>,
    /// BFS layer of each horizontal chord, `−1` for unreached or pruned.
    dist: Vec<i32>,
    /// BFS frontier, used as a flat queue with a read cursor.
    queue: Vec<u32>,
}

/// Every buffer the per-mask analysis needs, allocated once for the whole sweep.
///
/// A `6·(n−1)`-mask sweep at 65³ visits 384 masks five times over; allocating
/// per mask would put the allocator in the middle of the timing clause.
#[derive(Debug, Default)]
struct Scratch {
    /// Component id per mask cell, `−1` where clear.
    label: Vec<i32>,
    /// Flood-fill stack.
    stack: Vec<u32>,
    /// Boundary cycles per component.
    cycles: Vec<u32>,
    /// Left turns per component.
    convex: Vec<u32>,
    /// Right turns per component.
    reflex: Vec<u32>,
    /// `(grid point · 4 + direction)`, set once that directed edge is walked.
    walked: Vec<bool>,
    /// Horizontal chord index covering each grid point, `−1` for none.
    hcov: Vec<i32>,
    /// Horizontal chords as `[y, x_start, x_end]`.
    hchord: Vec<[i32; 3]>,
    /// Vertical chords as `[x, y_start, y_end]`.
    vchord: Vec<[i32; 3]>,
    /// Intersecting pairs as `[horizontal, vertical]`.
    edges: Vec<[u32; 2]>,
    /// CSR row starts over the horizontal chords.
    head: Vec<u32>,
    /// CSR column entries: vertical chord indices.
    adj: Vec<u32>,
    /// CSR fill cursor.
    cursor: Vec<u32>,
    /// Hopcroft–Karp state.
    matching: MatchState,
    /// Working copy of a mask for [`greedy_quad_count`].
    merge: Vec<bool>,
}

/// Label the 4-connected components of the set cells, returning how many.
///
/// 4-connected and not 8-connected, and that is the whole reason the pinch
/// hazard is survivable: a rectangle of set cells is 4-connected, so
/// `K_min(mask) = Σ K_min(component)` needs no theorem, and two cells touching
/// only at a corner are correctly two problems.
fn label_components(mask: &Mask, s: &mut Scratch) -> usize {
    let (du, dv) = (mask.du, mask.dv);
    s.label.clear();
    s.label.resize(du * dv, -1);
    let mut next = 0i32;
    for start in 0..du * dv {
        if !mask.cell[start] || s.label[start] >= 0 {
            continue;
        }
        s.label[start] = next;
        s.stack.clear();
        s.stack.push(start as u32);
        while let Some(k) = s.stack.pop() {
            let k = k as usize;
            let x = (k % du) as isize;
            let y = (k / du) as isize;
            for [dx, dy] in NEIGHBOURS {
                let nx = x + dx;
                let ny = y + dy;
                if nx < 0 || ny < 0 || nx >= du as isize || ny >= dv as isize {
                    continue;
                }
                let m = nx as usize + du * ny as usize;
                if mask.cell[m] && s.label[m] < 0 {
                    s.label[m] = next;
                    s.stack.push(m as u32);
                }
            }
        }
        next += 1;
    }
    next as usize
}

/// Walk every boundary cycle, tallying turns per component; return the pinch
/// count.
///
/// The rule is **solid on the left, leftmost turn available at each vertex**. At
/// an ordinary vertex out-degree is one and the rule is vacuous. At a
/// checkerboard vertex there are two incoming and two outgoing edges and the
/// leftmost rule pairs each solid quadrant's rays with each other — its own
/// corner — so the pinch yields two convex vertices, the walk stays inside one
/// component, and a component with a pinch reports `cycles = 1` rather than
/// counting the enclosure it necessarily has as a hole.
fn trace_boundary(mask: &Mask, s: &mut Scratch, ncomp: usize) -> i64 {
    let (du, dv) = (mask.du, mask.dv);
    let gw = du + 1;
    s.walked.clear();
    s.walked.resize(gw * (dv + 1) * 4, false);
    s.cycles.clear();
    s.cycles.resize(ncomp, 0);
    s.convex.clear();
    s.convex.resize(ncomp, 0);
    s.reflex.clear();
    s.reflex.resize(ncomp, 0);
    let mut pinches = 0i64;

    for py in 0..=dv as isize {
        for px in 0..=du as isize {
            let q = mask.quadrants(px, py);
            let [sw, se, nw, ne] = q;
            if (sw && ne && !se && !nw) || (se && nw && !sw && !ne) {
                pinches += 1;
            }
            for d in 0..4usize {
                if !leaves(q, d) || s.walked[(px as usize + gw * py as usize) * 4 + d] {
                    continue;
                }
                let lc = left_cell(px, py, d);
                let comp = s.label[lc[0] as usize + du * lc[1] as usize] as usize;
                s.cycles[comp] += 1;
                let (mut cx, mut cy, mut cd) = (px, py, d);
                loop {
                    let slot = (cx as usize + gw * cy as usize) * 4 + cd;
                    if s.walked[slot] {
                        break;
                    }
                    s.walked[slot] = true;
                    let nx = cx + STEP[cd][0];
                    let ny = cy + STEP[cd][1];
                    let nq = mask.quadrants(nx, ny);
                    let left = (cd + 1) % 4;
                    let right = (cd + 3) % 4;
                    let nd = if leaves(nq, left) {
                        s.convex[comp] += 1;
                        left
                    } else if leaves(nq, cd) {
                        cd
                    } else {
                        assert!(
                            leaves(nq, right),
                            "boundary walk left ({nx}, {ny}) with no outgoing edge, so the \
                             directed boundary is not in-degree = out-degree and every number \
                             below is meaningless"
                        );
                        s.reflex[comp] += 1;
                        right
                    };
                    cx = nx;
                    cy = ny;
                    cd = nd;
                }
            }
        }
    }
    pinches
}

/// Enumerate the good diagonals of a mask.
///
/// On grid line `y`, the predicate `solid(k, y−1) && solid(k, y)` marks the mask
/// cells an interior horizontal segment can cross. Each **maximal run** of it
/// yields at most one chord, because a grid point strictly inside a run has all
/// four quadrants solid and so cannot be reflex — which is also the proof that no
/// chord contains a third concave vertex and that no two horizontal chords
/// overlap or share an endpoint. The run's two ends are chord endpoints exactly
/// when both are reflex with the direction that points into the run.
fn find_chords(mask: &Mask, s: &mut Scratch) {
    let (du, dv) = (mask.du, mask.dv);
    s.hchord.clear();
    s.vchord.clear();

    for y in 0..=dv as isize {
        let mut k = 0isize;
        while k < du as isize {
            if !(mask.at(k, y - 1) && mask.at(k, y)) {
                k += 1;
                continue;
            }
            let a = k;
            while k < du as isize && mask.at(k, y - 1) && mask.at(k, y) {
                k += 1;
            }
            let starts_east = reflex_clear(mask.quadrants(a, y)).is_some_and(|c| chord_dirs(c).0);
            let ends_west = reflex_clear(mask.quadrants(k, y)).is_some_and(|c| !chord_dirs(c).0);
            if starts_east && ends_west {
                s.hchord.push([y as i32, a as i32, k as i32]);
            }
        }
    }

    for x in 0..=du as isize {
        let mut k = 0isize;
        while k < dv as isize {
            if !(mask.at(x - 1, k) && mask.at(x, k)) {
                k += 1;
                continue;
            }
            let a = k;
            while k < dv as isize && mask.at(x - 1, k) && mask.at(x, k) {
                k += 1;
            }
            let starts_north = reflex_clear(mask.quadrants(x, a)).is_some_and(|c| chord_dirs(c).1);
            let ends_south = reflex_clear(mask.quadrants(x, k)).is_some_and(|c| !chord_dirs(c).1);
            if starts_north && ends_south {
                s.vchord.push([x as i32, a as i32, k as i32]);
            }
        }
    }
}

/// Build the chord intersection graph as a CSR over the horizontal chords.
///
/// Two axis-parallel segments lying on grid lines can only meet at a grid point,
/// so the whole graph is found by marking each horizontal chord's grid points and
/// then walking each vertical chord's — `O(cells)` rather than the `O(H · V)` of
/// testing every pair. Horizontal chords are pairwise disjoint, so at most one
/// covers any given grid point and no edge is discovered twice.
fn build_chord_graph(mask: &Mask, s: &mut Scratch) {
    let gw = mask.du + 1;
    s.hcov.clear();
    s.hcov.resize(gw * (mask.dv + 1), -1);
    for (i, &[y, xa, xb]) in s.hchord.iter().enumerate() {
        for x in xa..=xb {
            s.hcov[x as usize + gw * y as usize] = i as i32;
        }
    }

    s.edges.clear();
    for (j, &[x, ya, yb]) in s.vchord.iter().enumerate() {
        for y in ya..=yb {
            let i = s.hcov[x as usize + gw * y as usize];
            if i >= 0 {
                s.edges.push([i as u32, j as u32]);
            }
        }
    }

    let nh = s.hchord.len();
    s.head.clear();
    s.head.resize(nh + 1, 0);
    for e in &s.edges {
        s.head[e[0] as usize + 1] += 1;
    }
    for i in 0..nh {
        s.head[i + 1] += s.head[i];
    }
    s.cursor.clear();
    s.cursor.extend_from_slice(&s.head);
    s.adj.clear();
    s.adj.resize(s.edges.len(), 0);
    for e in &s.edges {
        let slot = &mut s.cursor[e[0] as usize];
        s.adj[*slot as usize] = e[1];
        *slot += 1;
    }
}

/// One layered DFS augmentation from horizontal chord `i`.
///
/// Follows only edges that advance a BFS layer, and sets `dist[i] = −1` on
/// failure so no later augmentation in the same phase revisits it. Depth is
/// bounded by the phase's shortest augmenting path length, which Hopcroft–Karp
/// keeps at `O(√V)`.
fn augment(i: usize, head: &[u32], adj: &[u32], m: &mut MatchState) -> bool {
    for k in head[i]..head[i + 1] {
        let j = adj[k as usize] as usize;
        let partner = m.match_v[j];
        let advances = partner < 0
            || (m.dist[partner as usize] == m.dist[i] + 1
                && augment(partner as usize, head, adj, m));
        if advances {
            m.match_v[j] = i as i32;
            m.match_h[i] = j as i32;
            return true;
        }
    }
    m.dist[i] = -1;
    false
}

/// Maximum matching in the bipartite chord intersection graph, Hopcroft–Karp.
///
/// This is the number Eppstein's formula actually needs. `g` is a maximum
/// independent set of chords; in a bipartite graph a maximum independent set has
/// size `|V| − |maximum matching|` (König–Gallai), so
/// `g = (horizontal + vertical chords) − matching`. Each phase layers the graph
/// by BFS from the free horizontal chords and then augments greedily along those
/// layers; the loop ends when a phase augments nothing, which is Berge's lemma
/// and therefore an exact termination condition rather than a heuristic one.
fn hopcroft_karp(head: &[u32], adj: &[u32], nv: usize, m: &mut MatchState) -> i64 {
    let nh = head.len().saturating_sub(1);
    m.match_h.clear();
    m.match_h.resize(nh, -1);
    m.match_v.clear();
    m.match_v.resize(nv, -1);
    let mut matched = 0i64;
    if nh == 0 || nv == 0 {
        return 0;
    }
    loop {
        m.dist.clear();
        m.dist.resize(nh, -1);
        m.queue.clear();
        for i in 0..nh {
            if m.match_h[i] < 0 {
                m.dist[i] = 0;
                m.queue.push(i as u32);
            }
        }
        let mut read = 0usize;
        while read < m.queue.len() {
            let i = m.queue[read] as usize;
            read += 1;
            for k in head[i]..head[i + 1] {
                let j = adj[k as usize] as usize;
                let partner = m.match_v[j];
                if partner >= 0 && m.dist[partner as usize] < 0 {
                    m.dist[partner as usize] = m.dist[i] + 1;
                    m.queue.push(partner as u32);
                }
            }
        }
        let mut grew = 0i64;
        for i in 0..nh {
            if m.match_h[i] < 0 && augment(i, head, adj, m) {
                grew += 1;
            }
        }
        if grew == 0 {
            break;
        }
        matched += grew;
    }
    matched
}

/// Eppstein's `n/2 + h − g − 1`, summed over one mask's components.
///
/// `matching_ns` accumulates only the [`hopcroft_karp`] call, which is the
/// registered `matching_ms`; everything else in this function is charged to
/// `optimum_ms`.
fn optimum_of(mask: &Mask, s: &mut Scratch, matching_ns: &mut u128) -> Optimum {
    let faces = mask.cell.iter().filter(|&&b| b).count() as i64;
    if faces == 0 {
        return Optimum::default();
    }

    let ncomp = label_components(mask, s);
    let mut out = Optimum {
        faces,
        nonempty_masks: 1,
        components: ncomp as i64,
        ..Optimum::default()
    };
    out.pinch_vertices = trace_boundary(mask, s, ncomp);
    find_chords(mask, s);
    build_chord_graph(mask, s);

    let nv = s.vchord.len();
    let started = Instant::now();
    out.matched = hopcroft_karp(&s.head, &s.adj, nv, &mut s.matching);
    *matching_ns += started.elapsed().as_nanos();

    out.chords_h = s.hchord.len() as i64;
    out.chords_v = nv as i64;
    out.good_diagonals_g = out.chords_h + out.chords_v - out.matched;

    let mut per_component = 0i64;
    for c in 0..ncomp {
        let convex = i64::from(s.convex[c]);
        let reflex = i64::from(s.reflex[c]);
        let holes = i64::from(s.cycles[c]) - 1;
        // The two readings of the same formula, from separately measured
        // quantities. Equal iff `convex − reflex == 4 − 4h`, the discrete
        // turning-number theorem for a rectilinear region.
        let from_vertices = (convex + reflex) / 2 + holes - 1;
        let from_reflex = reflex - holes + 1;
        if from_vertices != from_reflex {
            out.identity_failures += 1;
        }
        out.vertices_n += convex + reflex;
        out.reflex_r += reflex;
        out.holes_h += holes;
        per_component += from_vertices;
    }
    out.rectangles = per_component - out.good_diagonals_g;
    out
}

// ─── the independent oracle ─────────────────────────────────────────────────

/// Everything the exhaustive cross-check needs, allocated once.
#[derive(Debug, Default)]
struct Verify {
    /// CSR row starts over components.
    comp_start: Vec<u32>,
    /// CSR fill cursor.
    cursor: Vec<u32>,
    /// Component cells as `[x, y]`, grouped by component.
    comp_cells: Vec<[i32; 2]>,
    /// One component in isolation, so the real analysis can be run on it.
    sub: Mask,
    /// The component under test, in row-major `(y, x)` order.
    cells: Vec<[i32; 2]>,
    /// Bit index of each cell of the component's bounding box, `−1` for absent.
    box_idx: Vec<i8>,
    /// Bounding-box extent along `x`.
    box_w: i32,
    /// Bounding-box extent along `y`.
    box_h: i32,
    /// Bounding-box origin along `x`.
    box_x0: i32,
    /// Bounding-box origin along `y`.
    box_y0: i32,
    /// Memoised partition size per remaining-cell subset.
    memo_val: Vec<u8>,
    /// Generation stamp per subset, so the memo needs no clearing.
    memo_gen: Vec<u64>,
    /// Current generation.
    generation: u64,
    /// Components handed to the oracle.
    checked: i64,
    /// Of those, the ones that were not already a full rectangle.
    nonrectangular: i64,
    /// Of those, the ones where the formula and the oracle disagreed.
    disagreements: i64,
}

impl Verify {
    /// A verifier with its memo sized for `BRUTE_LIMIT` cells.
    fn new() -> Self {
        Self {
            memo_val: vec![0; 1 << BRUTE_LIMIT],
            memo_gen: vec![0; 1 << BRUTE_LIMIT],
            generation: 0,
            ..Self::default()
        }
    }

    /// The subset bit of cell `(x, y)`, or `None` if the component has no such
    /// cell.
    fn bit_at(&self, x: i32, y: i32) -> Option<u32> {
        let bx = x - self.box_x0;
        let by = y - self.box_y0;
        if bx < 0 || by < 0 || bx >= self.box_w || by >= self.box_h {
            return None;
        }
        let slot = self.box_idx[(bx + self.box_w * by) as usize];
        if slot < 0 { None } else { Some(1u32 << slot) }
    }
}

/// The exact minimum rectangle partition of the loaded component, by exhaustive
/// search with memoisation.
///
/// This function knows nothing about Eppstein, holes, chords or matchings, which
/// is the point: it is the ground truth `brute_disagreements` is measured
/// against. The search is complete because the cells are in row-major order, so
/// the lowest set bit of `rem` is the minimal uncovered cell and **must** be the
/// bottom-left corner of whichever rectangle covers it; enumerating every width
/// and height from there enumerates every partition.
fn brute_min_rectangles(rem: u32, v: &mut Verify) -> u8 {
    if rem == 0 {
        return 0;
    }
    let slot = rem as usize;
    if v.memo_gen[slot] == v.generation {
        return v.memo_val[slot];
    }
    let [ox, oy] = v.cells[rem.trailing_zeros() as usize];
    let mut best = u8::MAX;
    let mut width = 0i32;
    'width: loop {
        width += 1;
        let mut bits = 0u32;
        for dx in 0..width {
            match v.bit_at(ox + dx, oy) {
                Some(b) if b & rem != 0 => bits |= b,
                _ => break 'width,
            }
        }
        let mut height = 1i32;
        loop {
            let deeper = brute_min_rectangles(rem & !bits, v);
            best = best.min(deeper + 1);
            let mut next = 0u32;
            let mut grew = true;
            for dx in 0..width {
                match v.bit_at(ox + dx, oy + height) {
                    Some(b) if b & rem != 0 => next |= b,
                    _ => {
                        grew = false;
                        break;
                    }
                }
            }
            if !grew {
                break;
            }
            bits |= next;
            height += 1;
        }
    }
    v.memo_gen[slot] = v.generation;
    v.memo_val[slot] = best;
    best
}

/// Cross-check Eppstein's formula against exhaustive search on every small
/// component of one mask.
///
/// The component's cells are collected **before** anything else runs, because
/// isolating one of them re-labels `s.label` and would otherwise destroy the
/// list being iterated. Each small component is then rebuilt as a mask of its
/// own and put through [`optimum_of`] — the real code path, not a copy of it —
/// and its answer compared with [`brute_min_rectangles`].
fn verify_mask(mask: &Mask, s: &mut Scratch, v: &mut Verify) {
    let (du, dv) = (mask.du, mask.dv);
    let ncomp = label_components(mask, s);
    if ncomp == 0 {
        return;
    }

    v.comp_start.clear();
    v.comp_start.resize(ncomp + 1, 0);
    for k in 0..du * dv {
        if s.label[k] >= 0 {
            v.comp_start[s.label[k] as usize + 1] += 1;
        }
    }
    for i in 0..ncomp {
        v.comp_start[i + 1] += v.comp_start[i];
    }
    v.cursor.clear();
    v.cursor.extend_from_slice(&v.comp_start);
    v.comp_cells.clear();
    v.comp_cells.resize(v.comp_start[ncomp] as usize, [0i32; 2]);
    for k in 0..du * dv {
        if s.label[k] >= 0 {
            let slot = &mut v.cursor[s.label[k] as usize];
            v.comp_cells[*slot as usize] = [(k % du) as i32, (k / du) as i32];
            *slot += 1;
        }
    }

    for c in 0..ncomp {
        let lo = v.comp_start[c] as usize;
        let hi = v.comp_start[c + 1] as usize;
        if hi - lo > BRUTE_LIMIT {
            continue;
        }
        v.cells.clear();
        v.cells.extend_from_slice(&v.comp_cells[lo..hi]);
        // Row-major, so `rem`'s lowest set bit is the minimal uncovered cell.
        v.cells.sort_unstable_by_key(|&[x, y]| (y, x));

        let x0 = v.cells.iter().map(|c| c[0]).min().unwrap_or(0);
        let y0 = v.cells.iter().map(|c| c[1]).min().unwrap_or(0);
        let x1 = v.cells.iter().map(|c| c[0]).max().unwrap_or(0);
        let y1 = v.cells.iter().map(|c| c[1]).max().unwrap_or(0);
        v.box_x0 = x0;
        v.box_y0 = y0;
        v.box_w = x1 - x0 + 1;
        v.box_h = y1 - y0 + 1;
        v.box_idx.clear();
        v.box_idx.resize((v.box_w * v.box_h) as usize, -1);
        for (i, &[x, y]) in v.cells.iter().enumerate() {
            v.box_idx[((x - x0) + v.box_w * (y - y0)) as usize] = i as i8;
        }

        // The same component in isolation, so `optimum_of` sees exactly one
        // component and its answer is comparable with the oracle's.
        v.sub.du = v.box_w as usize;
        v.sub.dv = v.box_h as usize;
        v.sub.cell.clear();
        v.sub.cell.resize(v.sub.du * v.sub.dv, false);
        for &[x, y] in &v.cells {
            let i = (x - x0) as usize + v.sub.du * (y - y0) as usize;
            v.sub.cell[i] = true;
        }
        let mut ignored = 0u128;
        let formula = optimum_of(&v.sub, s, &mut ignored).rectangles;

        v.generation += 1;
        let full = if v.cells.len() == 32 {
            u32::MAX
        } else {
            (1u32 << v.cells.len()) - 1
        };
        let exact = i64::from(brute_min_rectangles(full, v));

        v.checked += 1;
        if v.cells.len() != (v.box_w * v.box_h) as usize {
            v.nonrectangular += 1;
        }
        if formula != exact {
            v.disagreements += 1;
        }
    }
}

// ─── the golden fixture ─────────────────────────────────────────────────────

/// The committed `greedy_quads` triangle count for one fixture row, or `None`
/// when the fixture carries no row at that resolution.
///
/// One-line scanner in the shape of `golden.rs:245`'s `field_of`, which is
/// `#[cfg(test)]` and so not reachable from a bench: find the key, then cut at
/// the first non-digit. The fixture is one object per line with fixed key order
/// and no nesting or escapes (`golden_hashes.json`), so this is exact rather than
/// lenient.
fn golden_triangles(fixture: &str, field: &str, samples: u32) -> Option<u64> {
    let needle =
        format!("\"algorithm\":\"greedy_quads\",\"field\":\"{field}\",\"samples\":{samples},");
    let line = fixture.lines().find(|l| l.contains(&needle))?;
    let key = "\"triangles\":";
    let rest = &line[line.find(key)? + key.len()..];
    let end = rest.find(|c: char| !c.is_ascii_digit())?;
    rest[..end].parse().ok()
}

// ─── one row ────────────────────────────────────────────────────────────────

/// Everything measured for one `(field, resolution)` pair.
#[derive(Clone, Debug)]
struct Row {
    /// `ReferenceField::NAME`.
    field: &'static str,
    /// Samples per axis.
    resolution: u32,
    /// Eppstein's totals over all `6·(n−1)` masks.
    totals: Optimum,
    /// Masks swept, which is `6·(n−1)`.
    masks: i64,
    /// Merged quads from the shipped extractor: `triangle_count() / 2`.
    greedy_rectangles: i64,
    /// `MeshBuffer::triangle_count()` of the merged mesh.
    greedy_triangles: i64,
    /// `MeshBuffer::vertex_count()` of the merged mesh.
    greedy_vertices: i64,
    /// Unmerged quads from `Merge::Off`: `M-56`'s face-culling baseline.
    culled_faces: i64,
    /// Whether this file's replication matched the shipped quad count.
    mask_quads_agree: bool,
    /// `match` where the golden fixture confirmed the count, `absent` where it
    /// has no row.
    golden_check: &'static str,
    /// Median, min and max of `GreedyQuads::extract`, in milliseconds.
    greedy_ms: [f64; 3],
    /// Median, min and max of occupancy → optimum, in milliseconds.
    optimum_ms: [f64; 3],
    /// Median, min and max of the Hopcroft–Karp calls alone, in milliseconds.
    matching_ms: [f64; 3],
    /// Components handed to the exhaustive oracle.
    brute_checked: i64,
    /// Of those, the ones that were not already a rectangle.
    brute_nonrectangular: i64,
    /// Of those, the ones where the formula and the oracle disagreed.
    brute_disagreements: i64,
}

impl Row {
    /// `greedy_rectangles / optimum_rectangles`, and whether it means anything.
    ///
    /// Both are zero exactly when the field has no visible face at all, and then
    /// the ratio is `0/0`. It is recorded as `0` with `ratio_defined = false`
    /// rather than as `1`, because `1` is the value that would say "greedy
    /// achieved the optimum" — the sentence that falsifies C2.
    fn ratio(&self) -> (f64, bool) {
        if self.totals.rectangles == 0 {
            (0.0, false)
        } else {
            (
                self.greedy_rectangles as f64 / self.totals.rectangles as f64,
                true,
            )
        }
    }

    /// `optimum_ms / greedy_ms` at the median, and the least favourable pairing
    /// of the recorded extremes.
    fn cost_ratio(&self) -> (f64, f64) {
        (
            self.optimum_ms[0] / self.greedy_ms[0],
            self.optimum_ms[2] / self.greedy_ms[1],
        )
    }

    /// Whether the reduction applied on this row.
    fn c1(&self) -> bool {
        self.totals.identity_failures == 0
            && self.brute_disagreements == 0
            && self.totals.rectangles <= self.greedy_rectangles
    }
}

/// Median, min and max of a small sample, in that order.
///
/// Sorted with `f64::total_cmp` and never `partial_cmp().unwrap()`: a NaN
/// comparison is a determinism leak (T-004), and `REPEATS` is odd so the median
/// is an element rather than an average.
fn spread(mut xs: Vec<f64>) -> [f64; 3] {
    xs.sort_by(f64::total_cmp);
    [xs[xs.len() / 2], xs[0], xs[xs.len() - 1]]
}

/// One full `(field, resolution)` measurement.
fn measure<F>(
    field: &F,
    name: &'static str,
    samples: u32,
    fixture: &str,
    s: &mut Scratch,
    v: &mut Verify,
) -> Row
where
    F: isomesh::fields::ReferenceField + Sdf<Scalar = f64>,
{
    let (shape, origin, cell_size) = common::grid::<f64, _>(field, samples);
    let cells = (shape.size()[0] - 1) as usize;

    // ── the shipped algorithm, and the face-culling control ─────────────────
    let mut greedy = GreedyQuads::<f64>::new();
    let mut mesh = MeshBuffer::<f64>::new();
    greedy
        .extract(field, &shape, origin, cell_size, &mut mesh)
        .expect("greedy quads accepts the reference grid");
    assert_eq!(
        mesh.vertex_count(),
        2 * mesh.triangle_count(),
        "{name} at {samples}: `emit_quad` is documented to push four vertices and two \
         triangles per quad, so `triangle_count() / 2` is only the quad count while this holds"
    );
    assert_eq!(
        mesh.triangle_count() % 2,
        0,
        "{name} at {samples}: an odd triangle count cannot be a quad mesh"
    );
    let greedy_rectangles = (mesh.triangle_count() / 2) as i64;

    let mut culling = GreedyQuads::<f64>::new();
    culling.set_merge(Merge::Off);
    let mut culled = MeshBuffer::<f64>::new();
    culling
        .extract(field, &shape, origin, cell_size, &mut culled)
        .expect("greedy quads accepts the reference grid with merging off");
    let culled_faces = (culled.triangle_count() / 2) as i64;

    let golden_check = match golden_triangles(fixture, name, samples) {
        Some(triangles) => {
            assert_eq!(
                triangles as usize,
                mesh.triangle_count(),
                "{name} at {samples}: the shipped extractor disagrees with the committed \
                 golden fixture, so `greedy_rectangles` is not the blessed number"
            );
            "match"
        }
        None => "absent",
    };

    // ── the occupancy both arms share ───────────────────────────────────────
    let occupancy = Occupancy::of(field, cells, origin, cell_size);
    let mut mask = Mask::default();

    // ── the replication control: are these greedy's masks? ──────────────────
    let mut replicated = 0i64;
    let mut faces = 0i64;
    for_each_mask(&occupancy, &mut mask, |m| {
        replicated += greedy_quad_count(m, &mut s.merge) as i64;
        faces += m.cell.iter().filter(|&&b| b).count() as i64;
    });
    assert_eq!(
        replicated, greedy_rectangles,
        "{name} at {samples}: this file's mask replication disagrees with the shipped quad \
         count, so the regions Eppstein's formula is being applied to are not the regions the \
         extractor partitioned"
    );
    assert_eq!(
        faces, culled_faces,
        "{name} at {samples}: this file's visible-face count disagrees with `Merge::Off`, so \
         the masks differ from the extractor's before any merging happens"
    );

    // ── the optimum, and the oracle that checks it ──────────────────────────
    let mut totals = Optimum::default();
    let mut ignored = 0u128;
    for_each_mask(&occupancy, &mut mask, |m| {
        totals.add(&optimum_of(m, s, &mut ignored));
    });

    let before = (v.checked, v.nonrectangular, v.disagreements);
    for_each_mask(&occupancy, &mut mask, |m| verify_mask(m, s, v));
    let brute_checked = v.checked - before.0;
    let brute_nonrectangular = v.nonrectangular - before.1;
    let brute_disagreements = v.disagreements - before.2;

    // ── the clocks ─────────────────────────────────────────────────────────
    let mut warm = MeshBuffer::<f64>::new();
    greedy
        .extract(field, &shape, origin, cell_size, &mut warm)
        .expect("greedy quads accepts the reference grid");
    let mut greedy_samples = Vec::with_capacity(REPEATS);
    for _ in 0..REPEATS {
        warm.reset();
        let started = Instant::now();
        greedy
            .extract(field, &shape, origin, cell_size, &mut warm)
            .expect("greedy quads accepts the reference grid");
        greedy_samples.push(started.elapsed().as_secs_f64() * 1e3);
    }

    let mut optimum_samples = Vec::with_capacity(REPEATS);
    let mut matching_samples = Vec::with_capacity(REPEATS);
    for _ in 0..REPEATS {
        let mut matching_ns = 0u128;
        let started = Instant::now();
        let mut repeat = Optimum::default();
        for_each_mask(&occupancy, &mut mask, |m| {
            repeat.add(&optimum_of(m, s, &mut matching_ns));
        });
        optimum_samples.push(started.elapsed().as_secs_f64() * 1e3);
        matching_samples.push(matching_ns as f64 * 1e-6);
        assert_eq!(
            repeat.rectangles, totals.rectangles,
            "{name} at {samples}: the optimum is not deterministic across repeats"
        );
    }

    Row {
        field: name,
        resolution: samples,
        totals,
        masks: 6 * cells as i64,
        greedy_rectangles,
        greedy_triangles: mesh.triangle_count() as i64,
        greedy_vertices: mesh.vertex_count() as i64,
        culled_faces,
        mask_quads_agree: true,
        golden_check,
        greedy_ms: spread(greedy_samples),
        optimum_ms: spread(optimum_samples),
        matching_ms: spread(matching_samples),
        brute_checked,
        brute_nonrectangular,
        brute_disagreements,
    }
}

fn main() {
    if !std::env::args().any(|arg| arg == "--bench") {
        return;
    }
    let prereg = isomesh::experiment!("P-165");

    common::experiment::run(prereg, |run| {
        let fixture_path =
            std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("golden_hashes.json");
        let fixture = std::fs::read_to_string(&fixture_path)
            .expect("golden_hashes.json is committed beside the crate's Cargo.toml");

        let mut scratch = Scratch::default();
        let mut verify = Verify::new();
        let mut rows: Vec<Row> = Vec::new();

        for_each_reference_field!(f64, |name, field| {
            for &samples in &RESOLUTIONS {
                let row = measure(&field, name, samples, &fixture, &mut scratch, &mut verify);
                println!(
                    "  {:15} {:3}³  n={:7} h={:5} g={:5} comps={:6} opt={:7} greedy={:7} \
                     ratio={:.6} pinch={:6}",
                    row.field,
                    row.resolution,
                    row.totals.vertices_n,
                    row.totals.holes_h,
                    row.totals.good_diagonals_g,
                    row.totals.components,
                    row.totals.rectangles,
                    row.greedy_rectangles,
                    row.ratio().0,
                    row.totals.pinch_vertices,
                );
                rows.push(row);
            }
        });

        // ── vacuity controls, every one before the first `record` ───────────
        let fields_with_g = rows
            .iter()
            .filter(|r| r.totals.good_diagonals_g > 0)
            .count();
        assert!(
            fields_with_g > 0,
            "VOID: not one row has a good diagonal, so `g` is identically zero, Eppstein's \
             bipartite reduction -- the part of the result that is a contribution -- is never \
             exercised, and `hopcroft_karp` could be deleted without moving a number"
        );
        let fields_with_h = rows.iter().filter(|r| r.totals.holes_h > 0).count();
        assert!(
            fields_with_h > 0,
            "VOID: not one row has a hole, so `h` is identically zero and the formula is being \
             checked in exactly the degenerate case C1's falsifier names"
        );
        let rows_with_pinch = rows.iter().filter(|r| r.totals.pinch_vertices > 0).count();
        assert!(
            rows_with_pinch > 0,
            "VOID: not one mask has a checkerboard vertex, so the leftmost-turn rule is asserted \
             over a population where it cannot differ from the naive pairing and the whole \
             degeneracy argument is unmeasured"
        );
        let nonrectangular: i64 = rows.iter().map(|r| r.brute_nonrectangular).sum();
        assert!(
            nonrectangular > 0,
            "VOID: the exhaustive oracle only ever saw rectangles, and it agrees with any \
             formula that returns 1 for a rectangle, so `brute_disagreements == 0` is a zero \
             that could not have been non-zero (M-44)"
        );
        let merging_pays = rows
            .iter()
            .filter(|r| r.culled_faces > r.greedy_rectangles)
            .count();
        assert!(
            merging_pays > 0,
            "VOID: merging saved nothing on any field, so `merge_saving` is 1 everywhere and the \
             face-culling control is not reproducing M-56"
        );

        // ── C2 is global: its falsifier is `every field` ────────────────────
        let c2 = rows.iter().any(|r| {
            let (ratio, defined) = r.ratio();
            defined && ratio > 1.0
        });

        for row in &rows {
            let (ratio, ratio_defined) = row.ratio();
            let (cost_ratio, cost_ratio_worst) = row.cost_ratio();
            let t = &row.totals;
            run.record(&[
                ("field", row.field.to_string()),
                ("resolution", row.resolution.to_string()),
                ("vertices_n", t.vertices_n.to_string()),
                ("holes_h", t.holes_h.to_string()),
                ("good_diagonals_g", t.good_diagonals_g.to_string()),
                ("optimum_rectangles", t.rectangles.to_string()),
                ("greedy_rectangles", row.greedy_rectangles.to_string()),
                ("ratio", format!("{ratio:.6}")),
                ("matching_ms", format!("{:.6}", row.matching_ms[0])),
                ("optimum_ms", format!("{:.6}", row.optimum_ms[0])),
                ("c1_holds", row.c1().to_string()),
                ("c2_holds", c2.to_string()),
                ("c3_holds", (cost_ratio < COST_LIMIT).to_string()),
                // ── extras (M-273) ──
                ("brute_checked", row.brute_checked.to_string()),
                ("brute_disagreements", row.brute_disagreements.to_string()),
                ("brute_nonrectangular", row.brute_nonrectangular.to_string()),
                ("chords_horizontal", t.chords_h.to_string()),
                ("chords_vertical", t.chords_v.to_string()),
                ("components", t.components.to_string()),
                ("cost_ratio", format!("{cost_ratio:.6}")),
                ("cost_ratio_worst", format!("{cost_ratio_worst:.6}")),
                ("culled_faces", row.culled_faces.to_string()),
                ("golden_check", row.golden_check.to_string()),
                ("greedy_max_ms", format!("{:.6}", row.greedy_ms[2])),
                ("greedy_min_ms", format!("{:.6}", row.greedy_ms[1])),
                ("greedy_ms", format!("{:.6}", row.greedy_ms[0])),
                ("greedy_triangles", row.greedy_triangles.to_string()),
                ("greedy_vertices", row.greedy_vertices.to_string()),
                ("identity_failures", t.identity_failures.to_string()),
                ("mask_quads_agree", row.mask_quads_agree.to_string()),
                ("masks", row.masks.to_string()),
                ("matching_max_ms", format!("{:.6}", row.matching_ms[2])),
                ("matching_min_ms", format!("{:.6}", row.matching_ms[1])),
                ("matching_pairs", t.matched.to_string()),
                (
                    "merge_saving",
                    format!(
                        "{:.6}",
                        if row.greedy_rectangles == 0 {
                            0.0
                        } else {
                            row.culled_faces as f64 / row.greedy_rectangles as f64
                        }
                    ),
                ),
                ("nonempty_masks", t.nonempty_masks.to_string()),
                ("optimum_max_ms", format!("{:.6}", row.optimum_ms[2])),
                ("optimum_min_ms", format!("{:.6}", row.optimum_ms[1])),
                (
                    "optimum_le_greedy",
                    (t.rectangles <= row.greedy_rectangles).to_string(),
                ),
                ("pinch_vertices", t.pinch_vertices.to_string()),
                ("ratio_defined", ratio_defined.to_string()),
                ("reflex_r", t.reflex_r.to_string()),
                ("repeats", REPEATS.to_string()),
                ("total_cost_ratio", format!("{:.6}", 1.0 + cost_ratio)),
            ]);
        }

        println!(
            "\n  vacuity: {fields_with_g} rows with g>0, {fields_with_h} with h>0, \
             {rows_with_pinch} with a pinch, {nonrectangular} non-rectangular components \
             checked exhaustively, {merging_pays} rows where merging paid"
        );
        println!(
            "  C2 (global): greedy {} the optimum on every field",
            if c2 { "does NOT achieve" } else { "achieves" }
        );
    });
}
