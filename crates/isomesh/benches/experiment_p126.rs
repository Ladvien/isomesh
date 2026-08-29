//! **P-126 — `O-12`'s remaining half — the dual vertex link at 2^27, as a nightly gate.**
//!
//! Ticket: R-072. Pre-registered before this harness existed.
//!
//! ```bash
//! cargo bench --bench experiment_p126
//! ```
//!
//! Writes `docs/experiments/p-126.csv`. **This bench is the gate and there is no
//! second path.** After the CSV is written it `exit(1)`s, naming the arm and the
//! pattern, if any non-control arm read `worst_link_components` above 1.
//! `.github/workflows/nightly.yml` runs it nightly and reports its exit code; it
//! contains no verdict logic, because one binary answers the question once and
//! then gates it forever.
//!
//! # What was missing
//!
//! **`P-63` closed the primal half of `O-12` by exhaustion and reported the dual
//! half `VACUOUS`, and the reason was the block rather than the algorithm.**
//! `M-374` ran eleven arms over 2^18 = 262,144 sign patterns of a **3 × 3 × 2**
//! block — the four cells around one grid edge, which is every cell that can
//! contribute a face to a Marching Cubes *edge* vertex — and read
//! `worst_link_components` 1 on every non-control arm against a control that
//! produced 5,302 link-defective vertices. Then it ran `surface_nets`,
//! `dual_contouring` and `manifold_dual_contouring` on that same block and came
//! back with `max_incident_faces` **2** on all three and 0 defects against
//! 524,288 critical cells, and `FINDINGS.md:10913` records the verdict as
//! *"FALSIFIED, and vacuously … because a dual vertex in this block carries at
//! most 2 faces"*.
//!
//! A dual vertex lives at a **cell centre**, not on a grid edge. Its incident
//! faces are the quads of the **twelve grid edges of its own cell**
//! (`dual.rs:648-651`: one quad per crossed grid edge, joining the four cells
//! around it), and the link's *connectivity* is decided by the cells across its
//! **six faces**, because two quads meet along a dual edge only if both cells
//! agree about which sheet owns the two grid edges. So the block has to be
//! **3 × 3 × 3 cells = 4^3 = 64 corners**: the smallest one in which a cell has
//! all 26 neighbours. An 18-corner block cannot host the configuration at all,
//! and a zero over a population that cannot exist is `M-44`'s rule, not a
//! measurement.
//!
//! **The named error this harness exists not to repeat.** `P-63`'s C3 derived
//! its expected population from the **critical** census — 524,288 cells that are
//! 2D- or 3D-ambiguous — when the population it needed was dual vertices
//! carrying a **complete** link, which in an 18-corner block is *empty*. So it
//! divided a zero by the wrong denominator. Every population here is derived
//! **from the block this sweep sweeps**, in closed form, before any number is
//! read: `expected_active_cells`, `expected_dual_vertices` and, on the
//! one-vertex-per-cell control, `expected_link_defective_vertices` are all exact
//! integers computed from the 256 sign bytes at start-up and **asserted** against
//! what the sweep counted.
//!
//! # The block, and the one place the registration's arithmetic does not close
//!
//! The registration fixes two numbers that cannot both describe a free sweep:
//! **64 corners** and **2^27 = 134,217,728 patterns**. Sixty-four free corners is
//! 2^64. The harness does not amend the registration
//! (`crates/isomesh/src/experiment.rs:27-31` forbids it); it reconciles the two
//! by naming the identification that makes them simultaneously true and then
//! reporting it as a column:
//!
//! **The sign field is 3-periodic on each axis.** The block as drawn is
//! `4 × 4 × 4 = 64` corners at integer coordinates `0..=3`; corner 3 *is* corner
//! 0 on every axis, so the 64 corners take **27 distinct values** and the sweep
//! is exactly `2^27 = 134,217,728` patterns. `patterns_swept` is asserted equal
//! to 134,217,728 — a sweep that silently narrowed its own domain is a failure
//! rather than a fast pass — and `corner_identification` records
//! `period_3_per_axis` on every row so no reader has to infer it.
//!
//! **Why this identification and not a 27-corner sub-lattice.** A 4^3 corner
//! lattice has orbits of size 8, 24, 24 and 8 under the cube group, and no union
//! of them is 27 — so *every* 27-bit sub-domain of a 64-corner block breaks the
//! cube symmetry and leaves some cell of the block with corners frozen.
//! 3-periodicity freezes nothing: the field extends to all of Z^3, so **every one
//! of the 27 cells has all 26 neighbours** and every dual vertex in the block
//! carries a complete link. That is the property `P-63`'s block did not have, and
//! it is the property C3 is about. The cost is stated rather than hidden: the
//! sweep is exhaustive over **3-periodic** sign fields, which are 2^27 of the
//! 2^64 configurations a free 4^3 block admits. Both consequences fall out of one
//! number the printed summary reports: **the link of one cell's vertex depends on
//! 20 of the 27 sites**, not 32, because the `+axis` and `-axis` neighbours of a
//! cell share their outer slab under the wrap and their sheet structures are
//! therefore correlated — and the 7 free sites mean each distinct link
//! configuration is visited 2^7 = **128** times.
//!
//! **Translation invariance is why all 27 cells are censused rather than one.**
//! Shifting a 3-periodic pattern is a bijection of the 2^27 domain and carries
//! cell `(i, j, k)` onto cell `(0, 0, 0)`, so the census over
//! (all patterns × 27 cells) is exactly 27 times the census over
//! (all patterns × 1 cell) and sees not one configuration more. Twenty-seven
//! cells is chosen for the population, not the coverage: `M-44` wants the
//! denominator to be the block, and `dual_vertices` then runs to
//! 5,067,767,808 instead of 187,695,104.
//!
//! # The arms
//!
//! All six run over the same 2^27 patterns in one binary (`M-281`), and the
//! difference between any two of them is one rule and not one compilation.
//!
//! | arm | stitching | vertex rule | ambiguous-face pairing | control | what it is for |
//! |---|---|---|---|---|---|
//! | `mdc/separate` | periodic | one per Marching Cubes cycle | separated everywhere — the shipped default (`manifold_dual_contouring.rs:270`) | no | **C1 on the shipped path** |
//! | `mdc/every_ambiguous_face_joined` | periodic | one per cycle | joined on every ambiguous face | no | the registration's `FaceAmbiguity::Connected`, which the crate does not ship |
//! | `mdc/mixed_consistent` | periodic | one per cycle | one bit per **grid face**, drawn from the pattern, so the two cells agree | no | stands in for `FaceAmbiguity::AsymptoticDecider` |
//! | `control/one_vertex_per_cell` | periodic | one vertex owning the whole cell (`dual.rs:90-99`) | separated | **yes** | `DualContouring` and `SurfaceNets`, whose pinch `M-53` measured at 128 non-manifold edges |
//! | `control/inconsistent_join_x` | periodic | one per cycle | each cell joins its **high x** face and separates its low one, so the two cells sharing an x-face always disagree | **yes** | breaks the one property C1's proof uses, and nothing else |
//! | `control/open_block` | **open** | one per cycle, separated | separated | **yes** | the shipped rule on the **unstitched** block: `P-63`'s truncation, at the right block size |
//!
//! **`FaceAmbiguity::Connected` does not exist.** `marching_cubes/ambiguity.rs:75-86`
//! ships exactly `Separate` and `AsymptoticDecider`. The registration names a
//! third variant; the harness reports that rather than inventing one, and
//! `mdc/every_ambiguous_face_joined` is what such a variant would do. It is a
//! **non-control** arm, because the registration's C2 as literally written cannot
//! fire: an always-joined rule is still a function of the shared face's four
//! corner signs, so both cells reach the same answer, and the argument below
//! shows a link cannot split when they do. That is why the control that carries
//! C2 is `control/one_vertex_per_cell` — a shipped extractor, not a synthetic
//! mangling — with `control/inconsistent_join_x` beside it as the arm that
//! isolates *which* property is load-bearing.
//!
//! **`mdc/mixed_consistent` is how the decider is covered without a float.** The
//! asymptotic decider's answer on a face is a function of the four corner
//! *magnitudes* (`ambiguity.rs:100-117`), which a sign sweep does not carry. What
//! the connectivity argument uses is not *which* answer but that **both cells
//! reach the same one** — `ambiguity.rs:50-57` is that guarantee. So this arm
//! draws one bit per grid face from the pattern itself and hands the same bit to
//! both cells, sweeping the whole family of crack-free pairings of which the
//! decider is one member. It is the only arm that reaches the length-12 cycles:
//! `max_incident_faces` is **7** under all-separate and all-joined and **12**
//! under a mixed mask, which is the cycle `table.rs:93-99`'s `MAX_TRIANGLES` and
//! `table.rs:104-140`'s `CENTROID_BASE` are about, and it is unreachable from
//! either extreme.
//!
//! # C1 has a proof, and the sweep is the check on the transcription as much as
//! # on the claim
//!
//! Two lemmas, both verified over their whole domains at start-up rather than
//! argued:
//!
//! 1. **Consecutive edges of a cycle share a face of the cell.** `next[e]` is the
//!    segment leaving `e` on one of `e`'s two faces (`table.rs:233-241`), so the
//!    pair lies on that face. Checked over all 256 cases × all 64 masks:
//!    `cycle_adjacency_violations` must be 0.
//! 2. **With the joined bit agreed, the two cells sharing a face induce the
//!    identical pairing on that face's cut edges.** Checked over
//!    3 axes × 2 sides × 256 near cases × 16 far patterns × both bits = 49,152
//!    combinations: `pairing_disagreements` must be 0.
//!
//! Together: every consecutive pair of a cycle is linked in the neighbour too, so
//! the cycle's link is connected and `worst_link_components` is 1. C1 is
//! therefore **decided before the harness runs**, and saying so is the point —
//! what the 2^27 sweep buys is that the bench-local integer model of the dual
//! topology is the shipped one, which the replica check below measures directly,
//! and that no reachable pattern escapes the argument.
//!
//! **The link is the quad mesh's, and that is a choice with a reason.**
//! `dual.rs:752-758` triangulates each quad on the `slot 0 – slot 2` diagonal.
//! That diagonal is an output-format artefact of the triangulator, and a link
//! walked through it answers a question about `emit_quad_axis` rather than about
//! the dual. Adding it can only *merge* link components, so a quad-level 1
//! implies a triangulated 1 — the direction C1 needs — and the replica check
//! measures both and asserts the inequality rather than asserting the argument.
//!
//! # The replica check: the model is asserted to be the shipped path
//!
//! Before the sweep, the harness runs the **shipped** `DualContouring` and
//! `ManifoldDualContouring` (under both shipped face rules) on an open
//! 5 × 5 × 5-sample grid over 1,024 sign patterns with SplitMix64 magnitudes —
//! `experiment_p63.rs:112-270`'s fixture at a block that can host a complete
//! link — and asserts, per pattern and per arm:
//!
//! * `mesh.positions.len()` equals the model's dual-vertex count, exactly;
//! * `mesh.indices.len()` equals `6 ×` the model's quad count, exactly, and the
//!   per-cell owned-and-available edge total equals `4 ×` that quad count;
//! * the **per-vertex incident-triangle count**, element for element in emission
//!   order, equals the model's prediction — 2 triangles where the owning cell is
//!   the quad's slot 0 or slot 2 and 1 where it is slot 1 or slot 3;
//! * on the decider arm, the two cells' joined bits agree on every interior grid
//!   face: `replica_decider_face_disagreements` is 0.
//!
//! `replica_bit_identical` carries the third of those and is on every row.
//! `experiment_p101.rs`'s `edge_slot` arm is the model for this discipline: a
//! transcription of shipped arithmetic is not attributable until it has been
//! asserted identical to the shipped path.
//!
//! # What is imported and what is transcribed
//!
//! Nothing in the cycle decomposition is transcribed, because none of it is
//! private: `marching_cubes::table` re-exports the cube primitives at
//! `table.rs:88-91` and exposes `segment_links` (`:244`), `AMBIGUOUS_FACES`
//! (`:202`) and `face_bit` (`:173`), and `marching_cubes::ambiguity::joined_mask`
//! (`:135`) is public too. The **cycle walk** is copied from
//! `manifold_dual_contouring.rs:225-240` and the **whole-cell rule** from
//! `dual.rs:90-99`; the **quad walk's cell-to-edge correspondence** is copied
//! from `dual.rs:158-177` and `dual.rs:715-758`, which are private, with the
//! source line on the row that uses it — `experiment_p117.rs:53-56` states that
//! discipline verbatim and `crates/isomesh/src/` is read-only for this row.
//!
//! # SHARE, recomputed before the numbers
//!
//! **None, and it is registered that way.** This row moves nothing, proposes no
//! source change, and no clause here is or may become a speedup claim, so `✗51`'s
//! share bar does not apply and there is no Amdahl ceiling to compute. What
//! stands in a share's place is the column each clause is read from:
//!
//! * **C1** — `worst_link_components`, over `dual_vertices` on the three
//!   non-control arms. The denominator is 5,067,767,808 dual vertices on
//!   `mdc/separate` and the same on `mdc/every_ambiguous_face_joined`; the
//!   quantity is an integer count of connected components and there is no
//!   tolerance in it.
//! * **C2** — `worst_link_components` and `link_defective_vertices` on the three
//!   control arms. `control/one_vertex_per_cell`'s figure is predicted in closed
//!   form — **113,246,208** = 27 cells × 2^19 × the 8 of 256 sign bytes whose cut
//!   set is face-disconnected — and asserted, so the control is not merely
//!   non-zero but the *right* non-zero.
//! * **C3** — `dual_vertices_with_complete_link`, per arm, never pooled, against
//!   `expected_dual_vertices` computed from the 256 sign bytes at start-up.
//!
//! # Which unit carries each verdict, and it is never a nanosecond
//!
//! Every clause is an **integer count over an enumerated population whose
//! denominator is exact by construction**. `wall_seconds` is recorded because the
//! registration names it and because a nightly gate's schedule is a fact worth
//! keeping, and it is **read by nothing**: no clause, no assertion and no exit
//! code depends on it. `M-280` and `✗24` are why that is said out loud — this
//! machine spans 1.96–5.62 GHz under `powersave`/`balance_performance`, and a
//! quantity that moves with the governor cannot be a gate. There is no
//! `common::counters::Probe` here either, and that is not an omission: a
//! retired-instruction count is the right unit for a *cost* clause and there is
//! no cost clause on this row.
//!
//! # The vacuity controls, all four of them
//!
//! 1. **`patterns_swept == 134_217_728`**, asserted, on every arm.
//! 2. **`dual_vertices_with_complete_link > 0` per arm and never pooled**,
//!    asserted — and the counter is shown able to read *less*:
//!    `control/open_block` runs the identical rule on the unstitched block, where
//!    only cell `(1, 1, 1)` has all twelve of its edges' quads, and
//!    `incomplete_link_vertices` there is far larger than
//!    `dual_vertices_with_complete_link`. On the five stitched arms
//!    `incomplete_link_vertices` is asserted **0**. A counter that could only ever
//!    report "complete" would be `P-63`'s error inverted.
//! 3. **The control arms' `link_defective_vertices > 0`**, asserted, and
//!    `worst_link_components > 1` — a control reading 1 **aborts the run** rather
//!    than recording a pass, which is the registration's own instruction.
//! 4. **The closed forms**: `active_cells` must equal
//!    27 × 2^19 × 254 = **3,595,567,104** on every arm, and `dual_vertices` must
//!    equal its per-rule closed form on the five arms that have one. These are
//!    exact integer identities over the block, so a fixture that drifted by one
//!    cell cannot pass them.
//!
//! # What one iteration costs, and how that number was arrived at
//!
//! Per pattern per arm the loop does 27 cell evaluations and nothing else — no
//! allocation, no mesh, no floating point, no field evaluation. One cell is: two
//! L1 loads out of a 4.6 KB nibble table for the sign byte (the three 9-bit
//! z-planes of the pattern are extracted once per pattern, and the byte is
//! `nibble[plane[k]][ij] | nibble[plane[k+1]][ij] << 4`), one 24-byte row copy out
//! of a component table of at most 393 KB, six 12-bit face masks of which about
//! 0.75 take the neighbour branch, and about 1.4 bitmask closures over at most
//! twelve clique masks. `u64::count_ones` is **not** used: `target_feature_popcnt`
//! is false in this build (`P-105` measured it: no `target-cpu`, no
//! `RUSTFLAGS`), so it would lower to a ~12-instruction SWAR sequence, and a
//! 4 KB 12-bit population table replaces it. A shipped extractor is invoked
//! **3,072 times in the replica check and zero times in the sweep** — the sweep
//! builds no mesh, because 134,217,728 extractions is not the experiment.
//!
//! That is roughly **15 cycles per cell, 400 per pattern per arm** and
//! `6 × 27 × 2^27 ≈ 2.2 × 10^10` cell evaluations for the whole sweep, so
//! `3.2 × 10^11` cycles — **57 to 164 seconds** of issue-limited work across this
//! machine's clock range. **That is a static accounting and not a measurement**:
//! the author of this file does not run it, per the phase's rule that
//! twenty-five concurrently building harnesses are what moved five of Phase 24's
//! results (`FINDINGS.md:12639-12648`). The run itself supplies the number, per
//! arm, in `wall_seconds`, and the registration's own estimate — 4 to 45 minutes
//! single-threaded — is the band the accounting sits inside.

mod common;

use std::time::Instant;

use isomesh::dual_contouring::DualContouring;
use isomesh::extractor::Extractor;
use isomesh::manifold_dual_contouring::ManifoldDualContouring;
use isomesh::marching_cubes::FaceAmbiguity;
use isomesh::marching_cubes::ambiguity::joined_mask;
use isomesh::marching_cubes::table::{
    AMBIGUOUS_FACES, EDGE_AXIS, EDGE_CORNERS, EDGE_COUNT, NO_EDGE, corner_inside, edge_index,
    edge_on_face, face_bit, is_inside, segment_links,
};
use isomesh::{MeshBuffer, RuntimeShape3, Sdf};

// ─── the block ──────────────────────────────────────────────────────────────

/// Cells per axis of the swept block: 3, so a cell has all 26 neighbours.
const CELLS_PER_AXIS: usize = 3;
/// 27 cells.
const BLOCK_CELLS: usize = CELLS_PER_AXIS * CELLS_PER_AXIS * CELLS_PER_AXIS;
/// The block as drawn: `4^3` corners at integer coordinates `0..=3`.
const BLOCK_CORNERS: usize = 64;
/// Distinct corner values under the period-3 identification.
const CORNER_SITES: usize = 27;
/// `2^27 = 134,217,728`, the registered domain size.
const PATTERNS: u32 = 1 << CORNER_SITES;
/// Cells in the largest block any code here builds — the replica grid's `4^3`.
const MAX_CELLS: usize = 64;
/// Sign patterns the replica check draws.
const REPLICA_PATTERNS: u32 = 1024;
/// No neighbour across this face.
const NO_CELL: u16 = u16::MAX;
/// No component owns this edge.
const NO_OWNER: u8 = u8::MAX;
/// Progress is printed every this many patterns.
const PROGRESS_STRIDE: u32 = 1 << 24;

/// One cell's surface components, as the quad walk needs them.
///
/// `edges[c]` is component `c`'s twelve-bit set of owned cube edges and
/// `owner[e]` is the component that owns edge `e` — the two halves of
/// `dual::CellVertices` (`dual.rs:59-71`), which is what the quad walk reads.
#[derive(Clone, Copy)]
struct Components {
    /// The cell's sign byte, carried so a defect can name it.
    case: u8,
    /// Twelve-bit set of cut edges.
    cut: u16,
    /// How many components the rule placed.
    ncomp: u8,
    /// Edges owned by each component, in the order the rule pushed them.
    edges: [u16; 4],
    /// Owning component per edge, or [`NO_OWNER`].
    owner: [u8; EDGE_COUNT],
}

impl Components {
    /// The cut set of a sign byte: an edge whose two corners disagree.
    fn cut_of(case: u8) -> u16 {
        let mut cut = 0u16;
        for (e, [lo, hi]) in EDGE_CORNERS.into_iter().enumerate() {
            if corner_inside(case, lo) != corner_inside(case, hi) {
                cut |= 1 << e;
            }
        }
        cut
    }

    /// One vertex per Marching Cubes cycle — `manifold_dual_contouring.rs:225-240`
    /// verbatim, over the shipped `segment_links`.
    fn per_cycle(case: u8, joined: u8) -> Self {
        let next = segment_links(case, joined);
        let mut out = Self {
            case,
            cut: Self::cut_of(case),
            ncomp: 0,
            edges: [0; 4],
            owner: [NO_OWNER; EDGE_COUNT],
        };
        let mut visited = 0u16;
        for start in 0..EDGE_COUNT as u8 {
            if next[start as usize] == NO_EDGE || visited & (1 << start) != 0 {
                continue;
            }
            let mut edges = 0u16;
            let mut current = start;
            while visited & (1 << current) == 0 {
                visited |= 1 << current;
                edges |= 1 << current;
                current = next[current as usize];
            }
            let slot = out.ncomp;
            out.edges[slot as usize] = edges;
            for (e, owner) in out.owner.iter_mut().enumerate() {
                if edges & (1 << e) != 0 {
                    *owner = slot;
                }
            }
            out.ncomp += 1;
        }
        assert_eq!(
            visited, out.cut,
            "case {case:#010b} mask {joined:#08b}: the cycles must cover exactly the cut edges"
        );
        out
    }

    /// One vertex owning the whole cell — `dual.rs:90-99`'s `push_whole_cell`,
    /// which is `DualContouring`'s and `SurfaceNets`' rule.
    fn whole_cell(case: u8) -> Self {
        let cut = Self::cut_of(case);
        let mut out = Self {
            case,
            cut,
            ncomp: u8::from(cut != 0),
            edges: [cut, 0, 0, 0],
            owner: [NO_OWNER; EDGE_COUNT],
        };
        for (e, owner) in out.owner.iter_mut().enumerate() {
            if cut & (1 << e) != 0 {
                *owner = 0;
            }
        }
        out
    }
}

/// The cube tables the quad walk needs, derived once.
struct CubeTables {
    /// `mirror[axis][edge]` is the same grid edge seen from the cell across a
    /// face of that axis. Two cells sharing a face of axis `a` are offset by
    /// `±e_a`, so their local corner indices differ by `1 << a`
    /// (`dual.rs:721-723` subtracts the offset that makes this true).
    mirror: [[u8; EDGE_COUNT]; 3],
    /// Per `(axis, side)`: the face's twelve-bit edge mask and its four edges.
    face: [[(u16, [u8; 4]); 2]; 3],
    /// Which corner of its quad the owning cell occupies, per local edge.
    ///
    /// `dual.rs:718-723` walks the quad's corners as
    /// `(du, dv) ∈ {(0,0), (1,0), (1,1), (0,1)}`, reaching each cell by
    /// subtracting `du` on `u = (axis+1) % 3` and `dv` on `v = (axis+2) % 3` from
    /// the edge's low grid point, and `dual.rs:168` names the cell's local edge
    /// from `lo = (du << u) | (dv << v)`. So the owning cell sits at the slot
    /// whose `(du, dv)` are the edge's low corner's `u` and `v` bits, and
    /// `dual.rs:753-757` puts slots 0 and 2 in two triangles and slots 1 and 3 in
    /// one.
    quad_slot: [u8; EDGE_COUNT],
    /// Population of a twelve-bit mask, because `target_feature_popcnt` is false.
    pop12: Vec<u8>,
}

impl CubeTables {
    fn new() -> Self {
        let mut mirror = [[0u8; EDGE_COUNT]; 3];
        for (axis, row) in mirror.iter_mut().enumerate() {
            for (e, slot) in row.iter_mut().enumerate() {
                let [lo, hi] = EDGE_CORNERS[e];
                *slot = edge_index(lo ^ (1 << axis), hi ^ (1 << axis));
            }
        }

        let mut face = [[(0u16, [0u8; 4]); 2]; 3];
        for (axis, sides) in face.iter_mut().enumerate() {
            for (side, entry) in sides.iter_mut().enumerate() {
                let mut mask = 0u16;
                let mut edges = [0u8; 4];
                let mut n = 0usize;
                for e in 0..EDGE_COUNT as u8 {
                    if edge_on_face(e, axis, side as u8) {
                        mask |= 1 << e;
                        edges[n] = e;
                        n += 1;
                    }
                }
                assert_eq!(n, 4, "a cube face has four edges");
                *entry = (mask, edges);
            }
        }

        let mut quad_slot = [0u8; EDGE_COUNT];
        for (e, slot) in quad_slot.iter_mut().enumerate() {
            let axis = EDGE_AXIS[e] as usize;
            let u = (axis + 1) % 3;
            let v = (axis + 2) % 3;
            let lo = EDGE_CORNERS[e][0];
            *slot = match ((lo >> u) & 1, (lo >> v) & 1) {
                (0, 0) => 0,
                (1, 0) => 1,
                (1, 1) => 2,
                _ => 3,
            };
        }

        let pop12 = (0u32..4096).map(|m| m.count_ones() as u8).collect();

        Self {
            mirror,
            face,
            quad_slot,
            pop12,
        }
    }
}

/// The cells of a block, their neighbours, and which of their quads exist.
///
/// Two blocks are built: the **stitched** 3 × 3 × 3 torus the sweep runs on,
/// where nothing is missing, and **open** blocks, where a quad exists only when
/// all four of its cells do — which is `dual.rs:697-699`'s `1..cells[u]` bound
/// expressed as a mask.
struct Block {
    cells: usize,
    /// `neighbour[cell][axis][side]`, or [`NO_CELL`].
    neighbour: Vec<[[u16; 2]; 3]>,
    /// Twelve-bit mask of local edges whose quad exists.
    edge_available: Vec<u16>,
    /// Six-bit mask of faces whose neighbour exists, in [`face_bit`] order.
    face_available: Vec<u8>,
}

impl Block {
    /// `dims` counts **cells**. `wrap` stitches the block into a torus.
    fn new(dims: [usize; 3], wrap: bool) -> Self {
        let cells = dims[0] * dims[1] * dims[2];
        assert!(
            cells <= MAX_CELLS,
            "the scratch arrays are sized for {MAX_CELLS}"
        );
        let index = |c: [usize; 3]| c[0] + dims[0] * (c[1] + dims[1] * c[2]);
        let mut neighbour = vec![[[NO_CELL; 2]; 3]; cells];
        let mut edge_available = vec![0u16; cells];
        let mut face_available = vec![0u8; cells];

        for z in 0..dims[2] {
            for y in 0..dims[1] {
                for x in 0..dims[0] {
                    let c = [x, y, z];
                    let n = index(c);
                    for axis in 0..3 {
                        for (side, step) in [-1i64, 1].into_iter().enumerate() {
                            let t = c[axis] as i64 + step;
                            let inside = if wrap {
                                true
                            } else {
                                t >= 0 && t < dims[axis] as i64
                            };
                            if inside {
                                let mut d = c;
                                d[axis] = t.rem_euclid(dims[axis] as i64) as usize;
                                neighbour[n][axis][side] = index(d) as u16;
                                face_available[n] |= face_bit(axis, side as u8);
                            }
                        }
                    }
                    // A quad exists when all four cells around its grid edge do.
                    // The edge's low grid point sits at `cell + du * e_u + dv * e_v`
                    // and the four cells are that point minus 0 or 1 on each of
                    // `u` and `v` (`dual.rs:715-723`), so the point's `u` and `v`
                    // coordinates must both lie in `1..dims`.
                    for e in 0..EDGE_COUNT {
                        let axis = EDGE_AXIS[e] as usize;
                        let u = (axis + 1) % 3;
                        let v = (axis + 2) % 3;
                        let lo = EDGE_CORNERS[e][0];
                        let pu = c[u] + usize::from((lo >> u) & 1 == 1);
                        let pv = c[v] + usize::from((lo >> v) & 1 == 1);
                        let ok = wrap || (pu >= 1 && pu < dims[u] && pv >= 1 && pv < dims[v]);
                        if ok {
                            edge_available[n] |= 1 << e;
                        }
                    }
                }
            }
        }

        Self {
            cells,
            neighbour,
            edge_available,
            face_available,
        }
    }
}

// ─── the link census ────────────────────────────────────────────────────────

/// What one arm's sweep counted.
#[derive(Clone, Copy, Default)]
struct Census {
    active_cells: u64,
    dual_vertices: u64,
    complete: u64,
    incomplete: u64,
    isolated: u64,
    defective: u64,
    worst: u32,
    worst_pattern: u32,
    worst_cell: u32,
    worst_case: u8,
    max_incident: u32,
}

/// Per-cell working state, reused across patterns so the loop allocates nothing.
struct Scratch {
    case: [u8; MAX_CELLS],
    comps: [Components; MAX_CELLS],
}

impl Scratch {
    fn new() -> Self {
        Self {
            case: [0; MAX_CELLS],
            comps: [Components::whole_cell(0); MAX_CELLS],
        }
    }
}

/// Connected components of every dual vertex's incident-quad link, folded into
/// `out`.
///
/// Two quads incident to a dual vertex meet along a **dual edge** exactly when
/// their two grid edges lie on a common face of the cell *and* the cell across
/// that face gives both of them to the same one of its own components — the quad
/// walk asks each of the four cells *which of its vertices owns this edge*
/// (`dual.rs:648-651`, `:725-749`), so two quads share a mesh edge only when both
/// endpoints agree. Counting components of that relation is precisely what an
/// edge census cannot see: two cones glued at a point give every edge two faces.
///
/// When `RECORD` is set, one component count per dual vertex is pushed in
/// emission order — cells in `z`, `y`, `x` order and components in the order the
/// rule pushed them, which is `dual.rs:487-497` and `:540-543`.
fn census<const RECORD: bool>(
    block: &Block,
    tables: &CubeTables,
    scratch: &Scratch,
    pattern: u32,
    out: &mut Census,
    per_vertex: &mut Vec<u32>,
) {
    for n in 0..block.cells {
        let cell = &scratch.comps[n];
        if cell.cut == 0 {
            continue;
        }
        out.active_cells += 1;

        let avail = block.edge_available[n];
        let faces = block.face_available[n];

        // The relation, as at most twelve clique masks. A face with exactly two
        // cut edges needs no lookup: `segment_links` pairs a face's single entry
        // with its single exit in *both* cells, so those two edges always land in
        // one component of the neighbour — checked over all 256 cases and all 64
        // masks by `two_cut_edge_face_splits` at start-up. Only an
        // **ambiguous** face, whose four edges are all cut, has a pairing to ask
        // about, and that is 8 of 64 corner-sign patterns per face.
        let mut clique = [0u16; 16];
        let mut cliques = 0usize;
        for axis in 0..3 {
            for side in 0..2usize {
                let (fmask, fedges) = tables.face[axis][side];
                let fcut = cell.cut & fmask;
                if fcut & fcut.wrapping_sub(1) == 0 {
                    continue;
                }
                if fcut != fmask {
                    let k = fcut & avail;
                    if k & k.wrapping_sub(1) != 0 {
                        clique[cliques] = k;
                        cliques += 1;
                    }
                    continue;
                }
                if faces & face_bit(axis, side as u8) == 0 {
                    continue;
                }
                let d = block.neighbour[n][axis][side] as usize;
                let across = &scratch.comps[d].owner;
                let mut blocks = [0u16; 4];
                for &e in &fedges {
                    let t = across[tables.mirror[axis][e as usize] as usize];
                    if t != NO_OWNER {
                        blocks[t as usize] |= 1 << e;
                    }
                }
                for b in blocks {
                    let k = b & avail;
                    if k & k.wrapping_sub(1) != 0 {
                        clique[cliques] = k;
                        cliques += 1;
                    }
                }
            }
        }
        let clique = &clique[..cliques];

        for ci in 0..cell.ncomp as usize {
            let owned = cell.edges[ci];
            out.dual_vertices += 1;
            if owned & !avail == 0 {
                out.complete += 1;
            } else {
                out.incomplete += 1;
            }
            let live = owned & avail;
            if live == 0 {
                out.isolated += 1;
                if RECORD {
                    per_vertex.push(0);
                }
                continue;
            }
            let incident = u32::from(tables.pop12[live as usize]);
            if incident > out.max_incident {
                out.max_incident = incident;
            }

            let mut remaining = live;
            let mut parts = 0u32;
            while remaining != 0 {
                let mut reached = remaining & remaining.wrapping_neg();
                loop {
                    let before = reached;
                    for &k in clique {
                        if k & reached != 0 {
                            reached |= k & live;
                        }
                    }
                    if reached == before {
                        break;
                    }
                }
                remaining &= !reached;
                parts += 1;
            }

            if parts > 1 {
                out.defective += 1;
            }
            if parts > out.worst {
                out.worst = parts;
                out.worst_pattern = pattern;
                out.worst_cell = n as u32;
                out.worst_case = cell.case;
            }
            if RECORD {
                per_vertex.push(parts);
            }
        }
    }
}

// ─── the arms ───────────────────────────────────────────────────────────────

/// How an ambiguous face's four cut edges are paired.
#[derive(Clone, Copy, PartialEq, Eq)]
enum Pairing {
    /// Separated everywhere — `FaceAmbiguity::Separate`, the shipped default.
    Separate,
    /// Joined on every ambiguous face. The registration's `Connected`, which
    /// `marching_cubes/ambiguity.rs:75-86` does not ship.
    EveryAmbiguousFace,
    /// One bit per **grid** face, drawn from the pattern and handed to both cells.
    ConsistentPerGridFace,
    /// Each cell joins its high `x` face and separates its low one, so the two
    /// cells sharing an `x`-face always disagree and the other four faces agree.
    HighXFaceOfEachCell,
}

/// What a cell's vertices are.
#[derive(Clone, Copy, PartialEq, Eq)]
enum Placement {
    /// One per Marching Cubes cycle — `CycleQef`.
    PerCycle,
    /// One owning the whole cell — `Qef`, and `SurfaceNets`.
    WholeCell,
}

struct Arm {
    name: &'static str,
    placement: Placement,
    pairing: Pairing,
    /// Whether the block is stitched into a torus.
    stitched: bool,
    is_control: bool,
}

const ARMS: [Arm; 6] = [
    Arm {
        name: "mdc/separate",
        placement: Placement::PerCycle,
        pairing: Pairing::Separate,
        stitched: true,
        is_control: false,
    },
    Arm {
        name: "mdc/every_ambiguous_face_joined",
        placement: Placement::PerCycle,
        pairing: Pairing::EveryAmbiguousFace,
        stitched: true,
        is_control: false,
    },
    Arm {
        name: "mdc/mixed_consistent",
        placement: Placement::PerCycle,
        pairing: Pairing::ConsistentPerGridFace,
        stitched: true,
        is_control: false,
    },
    Arm {
        name: "control/one_vertex_per_cell",
        placement: Placement::WholeCell,
        pairing: Pairing::Separate,
        stitched: true,
        is_control: true,
    },
    Arm {
        name: "control/inconsistent_join_x",
        placement: Placement::PerCycle,
        pairing: Pairing::HighXFaceOfEachCell,
        stitched: true,
        is_control: true,
    },
    Arm {
        name: "control/open_block",
        placement: Placement::PerCycle,
        pairing: Pairing::Separate,
        stitched: false,
        is_control: true,
    },
];

impl Arm {
    fn pairing_name(&self) -> &'static str {
        match self.pairing {
            Pairing::Separate => "separate",
            Pairing::EveryAmbiguousFace => "every_ambiguous_face_joined",
            Pairing::ConsistentPerGridFace => "consistent_per_grid_face",
            Pairing::HighXFaceOfEachCell => "high_x_face_of_each_cell",
        }
    }

    fn placement_name(&self) -> &'static str {
        match self.placement {
            Placement::PerCycle => "one_per_marching_cubes_cycle",
            Placement::WholeCell => "one_per_cell",
        }
    }

    /// Do the two cells sharing a face always agree about the pairing?
    fn pairing_is_consistent(&self) -> bool {
        self.pairing != Pairing::HighXFaceOfEachCell
    }

    /// The joined mask, where it is a function of the cell's own sign byte.
    fn joined_of_case(&self, case: u8) -> Option<u8> {
        let ambiguous = AMBIGUOUS_FACES[case as usize];
        match self.pairing {
            Pairing::Separate => Some(0),
            Pairing::EveryAmbiguousFace => Some(ambiguous),
            Pairing::HighXFaceOfEachCell => Some(face_bit(0, 1) & ambiguous),
            Pairing::ConsistentPerGridFace => None,
        }
    }
}

/// One arm's component lookup, prebuilt so the sweep does no arithmetic per cell.
enum ArmTable {
    /// 256 rows: the joined mask is a function of the cell's own sign byte.
    PerCase(Vec<Components>),
    /// `256 × 64` rows, indexed `case << 6 | joined`.
    PerCaseAndMask(Vec<Components>),
}

impl ArmTable {
    fn build(arm: &Arm) -> Self {
        match (arm.placement, arm.joined_of_case(0)) {
            (Placement::WholeCell, _) => {
                Self::PerCase((0..256).map(|c| Components::whole_cell(c as u8)).collect())
            }
            (Placement::PerCycle, Some(_)) => Self::PerCase(
                (0..256u32)
                    .map(|c| {
                        let case = c as u8;
                        let joined = arm
                            .joined_of_case(case)
                            .expect("this arm's mask is a function of the case");
                        Components::per_cycle(case, joined)
                    })
                    .collect(),
            ),
            (Placement::PerCycle, None) => Self::PerCaseAndMask(
                (0..256u32 * 64)
                    .map(|i| {
                        let case = (i >> 6) as u8;
                        let joined = (i & 63) as u8;
                        Components::per_cycle(case, joined & AMBIGUOUS_FACES[case as usize])
                    })
                    .collect(),
            ),
        }
    }

    /// Total components over the 256 sign bytes, where the mask is a function of
    /// the byte — the closed form `expected_dual_vertices` is built from.
    fn cycles_per_256_cases(&self) -> Option<u64> {
        match self {
            Self::PerCase(rows) => Some(rows.iter().map(|r| u64::from(r.ncomp)).sum()),
            Self::PerCaseAndMask(_) => None,
        }
    }
}

// ─── the sweep ──────────────────────────────────────────────────────────────

/// The three 9-bit `z`-planes of the pattern, and which nibble each cell reads.
///
/// Site `(a, b, c)` of the period-3 lattice is bit `a + 3b + 9c`, so plane `c` is
/// bits `9c..9c+9`. Cell `(i, j, k)` reads plane `k` for its `dz = 0` nibble and
/// plane `(k + 1) % 3` for its `dz = 1` nibble, both at row index `i + 3j`.
struct PlaneTables {
    /// `nibble[plane_value][i + 3j]`, four bits in cube-corner order
    /// `dx + 2 dy`.
    nibble: Vec<[u8; 9]>,
    /// Per cell: `(i + 3j, k, (k + 1) % 3)`.
    cell: [[u8; 3]; BLOCK_CELLS],
}

impl PlaneTables {
    fn new() -> Self {
        let mut nibble = vec![[0u8; 9]; 512];
        for (value, row) in nibble.iter_mut().enumerate() {
            for j in 0..3usize {
                for i in 0..3usize {
                    let bit = |a: usize, b: usize| ((value >> ((a % 3) + 3 * (b % 3))) & 1) as u8;
                    row[i + 3 * j] = bit(i, j)
                        | (bit(i + 1, j) << 1)
                        | (bit(i, j + 1) << 2)
                        | (bit(i + 1, j + 1) << 3);
                }
            }
        }
        let mut cell = [[0u8; 3]; BLOCK_CELLS];
        for k in 0..3usize {
            for j in 0..3usize {
                for i in 0..3usize {
                    cell[i + 3 * j + 9 * k] = [(i + 3 * j) as u8, k as u8, ((k + 1) % 3) as u8];
                }
            }
        }
        Self { nibble, cell }
    }
}

/// The grid face each of a cell's six faces belongs to, for the consistent
/// pairing rule.
///
/// A grid face is named by its axis and its **upper** cell, so the two cells
/// sharing it compute the same key and therefore the same bit.
fn grid_face_table(block: &Block) -> [[u8; 6]; BLOCK_CELLS] {
    let mut out = [[0u8; 6]; BLOCK_CELLS];
    for (n, row) in out.iter_mut().enumerate() {
        for axis in 0..3usize {
            for side in 0..2usize {
                let upper = if side == 1 {
                    block.neighbour[n][axis][1] as usize
                } else {
                    n
                };
                row[axis * 2 + side] = (axis * BLOCK_CELLS + upper) as u8;
            }
        }
    }
    out
}

/// SplitMix64, so the consistent pairing is a pure function of the pattern and
/// the sweep is byte-identical on every machine — `experiment_p63.rs:168-176`.
fn splitmix64(mut z: u64) -> u64 {
    z = z.wrapping_add(0x9E37_79B9_7F4A_7C15);
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    z ^ (z >> 31)
}

/// 81 grid-face bits for one pattern: three axes × 27 upper cells.
fn grid_face_bits(pattern: u32) -> u128 {
    let lo = splitmix64(u64::from(pattern));
    let hi = splitmix64(u64::from(pattern) ^ 0xA076_1D64_78BD_642F);
    (u128::from(hi) << 64) | u128::from(lo)
}

fn sweep(
    arm: &Arm,
    block: &Block,
    tables: &CubeTables,
    planes: &PlaneTables,
    grid_face: &[[u8; 6]; BLOCK_CELLS],
    table: &ArmTable,
) -> (Census, f64) {
    let mut out = Census::default();
    let mut scratch = Scratch::new();
    let mut sink = Vec::new();
    let started = Instant::now();

    for pattern in 0..PATTERNS {
        if pattern % PROGRESS_STRIDE == 0 {
            println!(
                "  {}: {:>4}% at {:>7.1} s",
                arm.name,
                pattern / (PATTERNS / 100),
                started.elapsed().as_secs_f64()
            );
        }

        let plane = [
            (pattern & 0x1FF) as usize,
            ((pattern >> 9) & 0x1FF) as usize,
            ((pattern >> 18) & 0x1FF) as usize,
        ];
        for (n, key) in planes.cell.iter().enumerate() {
            let [ij, k, k1] = *key;
            scratch.case[n] = planes.nibble[plane[k as usize]][ij as usize]
                | (planes.nibble[plane[k1 as usize]][ij as usize] << 4);
        }

        match table {
            ArmTable::PerCase(rows) => {
                for n in 0..BLOCK_CELLS {
                    scratch.comps[n] = rows[scratch.case[n] as usize];
                }
            }
            ArmTable::PerCaseAndMask(rows) => {
                let bits = grid_face_bits(pattern);
                for (n, faces) in grid_face.iter().enumerate() {
                    let mut joined = 0u8;
                    for (q, &site) in faces.iter().enumerate() {
                        joined |= (((bits >> site) & 1) as u8) << q;
                    }
                    let case = scratch.case[n];
                    scratch.comps[n] = rows[(usize::from(case) << 6) | usize::from(joined)];
                }
            }
        }

        census::<false>(block, tables, &scratch, pattern, &mut out, &mut sink);
    }

    (out, started.elapsed().as_secs_f64())
}

// ─── the fixture's own lemmas, checked over their whole domains ──────────────

/// Every consecutive pair of a cycle lies on a common face of the cell.
///
/// `next[e]` is the segment leaving `e` on one of `e`'s two faces
/// (`table.rs:233-241`), so this must hold identically. It is the first half of
/// C1's proof.
fn cycle_adjacency_violations() -> u64 {
    let mut bad = 0u64;
    for case in 0..256u32 {
        let case = case as u8;
        for joined in 0..64u8 {
            if joined & !AMBIGUOUS_FACES[case as usize] != 0 {
                continue;
            }
            let next = segment_links(case, joined);
            for (e, &f) in next.iter().enumerate() {
                if f == NO_EDGE {
                    continue;
                }
                let shared = (0..3).any(|axis| {
                    (0..2).any(|side| {
                        edge_on_face(e as u8, axis, side) && edge_on_face(f, axis, side)
                    })
                });
                if !shared {
                    bad += 1;
                }
            }
        }
    }
    bad
}

/// A face with exactly two cut edges puts both of them in one component.
///
/// The face's single entry is paired with its single exit in either cell, so the
/// two are consecutive in a cycle. This is what lets the census skip the
/// neighbour lookup on 7 of every 8 cut faces.
fn two_cut_edge_face_splits(tables: &CubeTables) -> u64 {
    let mut bad = 0u64;
    for case in 0..256u32 {
        let case = case as u8;
        for joined in 0..64u8 {
            if joined & !AMBIGUOUS_FACES[case as usize] != 0 {
                continue;
            }
            let c = Components::per_cycle(case, joined);
            for axis in 0..3 {
                for side in 0..2usize {
                    let (fmask, fedges) = tables.face[axis][side];
                    let fcut = c.cut & fmask;
                    if fcut == fmask || fcut & fcut.wrapping_sub(1) == 0 {
                        continue;
                    }
                    let owners: Vec<u8> = fedges
                        .iter()
                        .filter(|&&e| fcut & (1 << e) != 0)
                        .map(|&e| c.owner[e as usize])
                        .collect();
                    if owners.windows(2).any(|w| w[0] != w[1]) {
                        bad += 1;
                    }
                }
            }
        }
    }
    bad
}

/// With the joined bit agreed, the two cells sharing a face induce the identical
/// **pairing** on that face's cut edges.
///
/// `ambiguity.rs:50-57` is the argument — the two cells read the same four sample
/// values in orders that differ by a rotation and a reflection, and neither
/// transformation changes which diagonal is which. This is that claim as an
/// exhaustive check over 3 axes × 2 sides × 256 near bytes × 16 far patterns ×
/// both bits, and it is the second half of C1's proof.
///
/// What is compared is the **pairing**, not cycle membership. `segment_links`
/// links a face's entry to its exit, and two edges of one face can perfectly well
/// sit in the same *cycle* of one cell and in different cycles of the other — the
/// rest of each cell's faces differ. The pairing is what has to agree, because it
/// is what makes two consecutive edges of a cycle land in one component across
/// the face, and that is what `control/inconsistent_join_x` breaks on purpose.
///
/// Two distinct faces of a cube share exactly one edge, so a link whose two edges
/// both lie on this face was made *on* this face and nowhere else.
fn pairing_disagreements(tables: &CubeTables) -> u64 {
    let mut bad = 0u64;
    for axis in 0..3usize {
        for side in 0..2usize {
            let (near_mask, near_edges) = tables.face[axis][side];
            let (far_mask, far_edges) = tables.face[axis][1 - side];
            let far_corners: Vec<u8> = (0..8u8)
                .filter(|k| usize::from((k >> axis) & 1) == side)
                .collect();
            for near in 0..256u32 {
                let near = near as u8;
                for far in 0..16u32 {
                    // The cell across the face: its own shared-face corners are
                    // the near cell's, mirrored, and its four far corners are
                    // free.
                    let mut across = 0u8;
                    for k in 0..8u8 {
                        if usize::from((k >> axis) & 1) == 1 - side
                            && corner_inside(near, k ^ (1 << axis))
                        {
                            across |= 1 << k;
                        }
                    }
                    for (t, &k) in far_corners.iter().enumerate() {
                        if (far >> t) & 1 == 1 {
                            across |= 1 << k;
                        }
                    }
                    for bit in 0..2u8 {
                        let jn = if bit == 1 {
                            face_bit(axis, side as u8)
                        } else {
                            0
                        };
                        let jf = if bit == 1 {
                            face_bit(axis, 1 - side as u8)
                        } else {
                            0
                        };
                        let nn = segment_links(near, jn);
                        let ff = segment_links(across, jf);
                        // The pairing as an involution over the near cell's edge
                        // labels, from each side.
                        let mut pa = [NO_EDGE; EDGE_COUNT];
                        let mut pb = [NO_EDGE; EDGE_COUNT];
                        for &e in &near_edges {
                            let p = nn[e as usize];
                            if p != NO_EDGE && near_mask & (1 << p) != 0 {
                                pa[e as usize] = p;
                                pa[p as usize] = e;
                            }
                        }
                        for &e in &far_edges {
                            let p = ff[e as usize];
                            if p != NO_EDGE && far_mask & (1 << p) != 0 {
                                let me = tables.mirror[axis][e as usize];
                                let mp = tables.mirror[axis][p as usize];
                                pb[me as usize] = mp;
                                pb[mp as usize] = me;
                            }
                        }
                        if pa != pb {
                            bad += 1;
                        }
                    }
                }
            }
        }
    }
    bad
}

/// Sign bytes whose cut set is disconnected under "share a face of the cell".
///
/// One vertex per cell gives every cut edge to the same vertex, and the
/// neighbour then always agrees, so the link's components are exactly the
/// face-connected components of the cut set — a function of the byte alone. This
/// count is `control/one_vertex_per_cell`'s closed form.
fn face_disconnected_cases(tables: &CubeTables) -> u64 {
    let mut n = 0u64;
    for case in 0..256u32 {
        let cut = Components::cut_of(case as u8);
        if cut == 0 {
            continue;
        }
        let mut clique = Vec::new();
        for axis in 0..3 {
            for side in 0..2usize {
                let k = cut & tables.face[axis][side].0;
                if k & k.wrapping_sub(1) != 0 {
                    clique.push(k);
                }
            }
        }
        let mut remaining = cut;
        let mut parts = 0u32;
        while remaining != 0 {
            let mut reached = remaining & remaining.wrapping_neg();
            loop {
                let before = reached;
                for &k in &clique {
                    if k & reached != 0 {
                        reached |= k & cut;
                    }
                }
                if reached == before {
                    break;
                }
            }
            remaining &= !reached;
            parts += 1;
        }
        if parts > 1 {
            n += 1;
        }
    }
    n
}

// ─── the replica: the model asserted to be the shipped path ─────────────────

/// The trilinear interpolant of one sign pattern over an open sample grid.
///
/// `experiment_p63.rs:112-270`'s fixture at a block that can host a complete
/// dual link. Magnitudes come from SplitMix64 in `[1/4, 5/4)` so no corner sits
/// on the surface and the asymptotic decider has something to decide — `P-63`'s
/// own `M-44` control is why they are not all `±1`.
struct Lattice {
    size: [usize; 3],
    value: Vec<f64>,
}

impl Lattice {
    fn new(size: [usize; 3], pattern: u128, seed: u64) -> Self {
        let samples = size[0] * size[1] * size[2];
        let value = (0..samples)
            .map(|i| {
                let sign = if (pattern >> i) & 1 == 1 { -1.0 } else { 1.0 };
                let z = splitmix64(seed ^ (i as u64).wrapping_mul(0x9E37_79B9_7F4A_7C15));
                sign * (0.25 + (z >> 11) as f64 * (1.0 / 9_007_199_254_740_992.0))
            })
            .collect();
        Self { size, value }
    }

    fn at(&self, x: usize, y: usize, z: usize) -> f64 {
        self.value[x + self.size[0] * (y + self.size[1] * z)]
    }

    /// The containing cell and the fraction inside it, clamped to the block.
    /// Clamping is a boundary condition and not a fallback: the extractor only
    /// asks inside, and the decider's saddle solve can ask on a face.
    fn locate(&self, p: [f64; 3]) -> ([usize; 3], [f64; 3]) {
        let mut base = [0usize; 3];
        let mut frac = [0.0f64; 3];
        for k in 0..3 {
            let limit = self.size[k] - 2;
            let f = p[k].floor();
            let i = if f < 0.0 { 0 } else { (f as usize).min(limit) };
            base[k] = i;
            frac[k] = p[k] - i as f64;
        }
        (base, frac)
    }
}

impl Sdf for Lattice {
    type Scalar = f64;

    fn sample(&self, p: [f64; 3]) -> f64 {
        let (b, f) = self.locate(p);
        let mut acc = 0.0;
        for k in 0..8u32 {
            let d = [
                (k & 1) as usize,
                ((k >> 1) & 1) as usize,
                ((k >> 2) & 1) as usize,
            ];
            let w = [f[0], f[1], f[2]];
            let mut weight = 1.0;
            for axis in 0..3 {
                weight *= if d[axis] == 1 { w[axis] } else { 1.0 - w[axis] };
            }
            acc += weight * self.at(b[0] + d[0], b[1] + d[1], b[2] + d[2]);
        }
        acc
    }

    /// The trilinear's own analytic gradient. Not a central difference: the
    /// stencil would reach outside the block and measure the clamp.
    fn gradient(&self, p: [f64; 3]) -> [f64; 3] {
        let (b, f) = self.locate(p);
        let mut g = [0.0f64; 3];
        for k in 0..8u32 {
            let d = [
                (k & 1) as usize,
                ((k >> 1) & 1) as usize,
                ((k >> 2) & 1) as usize,
            ];
            let v = self.at(b[0] + d[0], b[1] + d[1], b[2] + d[2]);
            for axis in 0..3 {
                let mut w = if d[axis] == 1 { 1.0 } else { -1.0 };
                for other in 0..3 {
                    if other != axis {
                        w *= if d[other] == 1 {
                            f[other]
                        } else {
                            1.0 - f[other]
                        };
                    }
                }
                g[axis] += w * v;
            }
        }
        g
    }
}

/// What the replica check established.
struct Replica {
    patterns: u32,
    arms: usize,
    bit_identical: bool,
    vertices: u64,
    triangles: u64,
    decider_disagreements: u64,
    /// Vertices where the triangulated link has fewer components than the quad
    /// link, i.e. where `dual.rs:753-757`'s diagonal merged two of them.
    diagonal_merges: u64,
    /// Vertices where the two agree.
    link_agreements: u64,
}

/// Run the shipped extractors beside the model and assert they are the same mesh.
fn replica(tables: &CubeTables) -> Replica {
    let cells = [4usize, 4, 4];
    let size = [cells[0] + 1, cells[1] + 1, cells[2] + 1];
    let samples = size[0] * size[1] * size[2];
    assert!(samples <= 128, "the replica pattern is one u128");
    let block = Block::new(cells, false);
    let shape = RuntimeShape3::new([size[0] as u32, size[1] as u32, size[2] as u32])
        .expect("5x5x5 samples is a legal grid");

    let mut dc = DualContouring::<f64>::new();
    let mut mdc_separate = ManifoldDualContouring::<f64>::new();
    let mut mdc_decider = ManifoldDualContouring::<f64>::new();
    mdc_decider.set_face_ambiguity(FaceAmbiguity::AsymptoticDecider);

    let mut mesh = MeshBuffer::<f64>::new();
    let mut scratch = Scratch::new();
    let mut model_link = Vec::new();
    let mut incident = Vec::new();

    let mut out = Replica {
        patterns: REPLICA_PATTERNS,
        arms: 3,
        bit_identical: true,
        vertices: 0,
        triangles: 0,
        decider_disagreements: 0,
        diagonal_merges: 0,
        link_agreements: 0,
    };

    for trial in 0..REPLICA_PATTERNS {
        let bits = (u128::from(splitmix64(u64::from(trial))) << 64)
            | u128::from(splitmix64(u64::from(trial) ^ 0xDEAD_BEEF_CAFE_F00D));
        let field = Lattice::new(size, bits, 0x0000_2026 ^ u64::from(trial));

        // The eight corner values of every cell, and its sign byte, exactly as
        // `dual.rs:500-513` gathers them.
        let mut corners = vec![[0.0f64; 8]; block.cells];
        for z in 0..cells[2] {
            for y in 0..cells[1] {
                for x in 0..cells[0] {
                    let n = x + cells[0] * (y + cells[1] * z);
                    let mut case = 0u8;
                    for (c, slot) in corners[n].iter_mut().enumerate() {
                        let v = field.at(x + (c & 1), y + ((c >> 1) & 1), z + ((c >> 2) & 1));
                        *slot = v;
                        if is_inside(v) {
                            case |= 1 << c;
                        }
                    }
                    scratch.case[n] = case;
                }
            }
        }

        for arm in 0..3usize {
            for (n, corner) in corners.iter().enumerate() {
                let case = scratch.case[n];
                scratch.comps[n] = match arm {
                    0 => Components::whole_cell(case),
                    1 => Components::per_cycle(case, 0),
                    // `CycleQef::place` at `manifold_dual_contouring.rs:214-218`.
                    _ => Components::per_cycle(
                        case,
                        joined_mask(corner, AMBIGUOUS_FACES[case as usize]),
                    ),
                };
            }

            // `Extractor::extract_into` is generic over the field and the sink
            // (`extractor.rs:80-90`), so the trait is not object-safe and the
            // three arms are three calls rather than one through a `dyn`.
            mesh.reset();
            match arm {
                0 => dc.extract_into(&field, &shape, [0.0; 3], 1.0, &mut mesh),
                1 => mdc_separate.extract_into(&field, &shape, [0.0; 3], 1.0, &mut mesh),
                _ => mdc_decider.extract_into(&field, &shape, [0.0; 3], 1.0, &mut mesh),
            }
            .expect("extraction on a 5x5x5 grid");

            // The model's vertices, in emission order, with the number of
            // triangles each one must appear in.
            let mut predicted = Vec::new();
            let mut owned_and_available = 0u64;
            for n in 0..block.cells {
                let cell = &scratch.comps[n];
                let avail = block.edge_available[n];
                for ci in 0..cell.ncomp as usize {
                    let live = cell.edges[ci] & avail;
                    owned_and_available += u64::from(tables.pop12[live as usize]);
                    let mut triangles = 0u32;
                    for (e, &slot) in tables.quad_slot.iter().enumerate() {
                        if live & (1 << e) != 0 {
                            triangles += if slot % 2 == 0 { 2 } else { 1 };
                        }
                    }
                    predicted.push(triangles);
                }
            }
            let quads = owned_and_available / 4;
            assert_eq!(
                owned_and_available % 4,
                0,
                "every quad is a local edge of exactly four cells"
            );
            assert_eq!(
                mesh.positions.len(),
                predicted.len(),
                "arm {arm}: the shipped mesher placed a different number of vertices than the \
                 model — one vertex per surface component of every active cell"
            );
            assert_eq!(
                mesh.indices.len() as u64,
                6 * quads,
                "arm {arm}: two triangles per quad, three indices each"
            );

            incident.clear();
            incident.resize(mesh.positions.len(), 0u32);
            for &v in &mesh.indices {
                incident[v as usize] += 1;
            }
            out.bit_identical &= incident == predicted;
            assert!(
                out.bit_identical,
                "arm {arm}: the per-vertex incident-triangle counts differ from the model, so the \
                 bench-local quad walk is not the shipped one and no other arm's number is \
                 attributable"
            );

            out.vertices += mesh.positions.len() as u64;
            out.triangles += mesh.indices.len() as u64 / 3;

            // The quad link, from the model, beside the triangulated link, from
            // the mesh the consumer actually receives.
            model_link.clear();
            let mut ignored = Census::default();
            census::<true>(&block, tables, &scratch, 0, &mut ignored, &mut model_link);
            let walked = triangulated_link_components(&mesh);
            assert_eq!(walked.len(), model_link.len(), "one link count per vertex");
            for (mesh_parts, quad_parts) in walked.iter().zip(&model_link) {
                assert!(
                    mesh_parts <= quad_parts,
                    "arm {arm}: the quad diagonal at dual.rs:753-757 can only merge link \
                     components, never split them"
                );
                if mesh_parts == quad_parts {
                    out.link_agreements += 1;
                } else {
                    out.diagonal_merges += 1;
                }
            }

            if arm == 2 {
                out.decider_disagreements += decider_face_disagreements(&block, &corners);
            }
        }
    }

    out
}

/// Connected components of every vertex's incident-face link, walked on the
/// emitted triangles — `experiment_p63.rs:299-344`'s walk.
///
/// Two incident faces are adjacent when they share an edge **through** the
/// vertex, which is the vertex plus one other corner.
fn triangulated_link_components(mesh: &MeshBuffer<f64>) -> Vec<u32> {
    let mut incident: Vec<Vec<[u32; 2]>> = vec![Vec::new(); mesh.positions.len()];
    for tri in mesh.indices.as_chunks::<3>().0 {
        for (k, &v) in tri.iter().enumerate() {
            incident[v as usize].push([tri[(k + 1) % 3], tri[(k + 2) % 3]]);
        }
    }
    incident
        .iter()
        .map(|faces| {
            let n = faces.len();
            if n == 0 {
                return 0;
            }
            let mut seen = vec![false; n];
            let mut parts = 0u32;
            let mut stack = Vec::new();
            for start in 0..n {
                if seen[start] {
                    continue;
                }
                parts += 1;
                seen[start] = true;
                stack.push(start);
                while let Some(i) = stack.pop() {
                    for j in 0..n {
                        if seen[j] {
                            continue;
                        }
                        if faces[i].iter().any(|a| faces[j].contains(a)) {
                            seen[j] = true;
                            stack.push(j);
                        }
                    }
                }
            }
            parts
        })
        .collect()
}

/// Interior grid faces on which the two cells' decider bits disagree.
///
/// `ambiguity.rs:50-57` says they cannot. This is that claim read off the
/// shipped `joined_mask` on real corner values, and it is what licenses
/// `mdc/mixed_consistent` to sweep consistent bits rather than run the decider.
fn decider_face_disagreements(block: &Block, corners: &[[f64; 8]]) -> u64 {
    let mut bad = 0u64;
    for n in 0..block.cells {
        let near = &corners[n];
        let mut near_case = 0u8;
        for (c, v) in near.iter().enumerate() {
            if is_inside(*v) {
                near_case |= 1 << c;
            }
        }
        let near_mask = joined_mask(near, AMBIGUOUS_FACES[near_case as usize]);
        for axis in 0..3usize {
            // Each interior grid face once, from its lower cell's high side.
            let d = block.neighbour[n][axis][1];
            if d == NO_CELL {
                continue;
            }
            let far = &corners[d as usize];
            let mut far_case = 0u8;
            for (c, v) in far.iter().enumerate() {
                if is_inside(*v) {
                    far_case |= 1 << c;
                }
            }
            let far_mask = joined_mask(far, AMBIGUOUS_FACES[far_case as usize]);
            let near_bit = near_mask & face_bit(axis, 1) != 0;
            let far_bit = far_mask & face_bit(axis, 0) != 0;
            if near_bit != far_bit {
                bad += 1;
            }
        }
    }
    bad
}

/// How many of the block's 27 corner sites decide one cell's link.
///
/// A dual vertex's link needs its own cell's cycles and, across each of its six
/// faces, the neighbour's — so seven cells' sign bytes. Under the period-3 wrap
/// the `+axis` and `-axis` neighbours of a cell **share their outer slab**, so
/// the seven cells span 8 + 3 × 4 = 20 sites rather than 8 + 6 × 4 = 32. The
/// remaining 7 sites are free, which is the sweep's 2^7 = 128-fold multiplicity.
fn link_determining_sites() -> usize {
    let site = |a: usize, b: usize, c: usize| a % 3 + 3 * (b % 3) + 9 * (c % 3);
    let mut cells = vec![[1usize; 3]];
    for axis in 0..3 {
        for step in [1usize, 2] {
            let mut c = [1usize; 3];
            c[axis] = (c[axis] + step) % 3;
            cells.push(c);
        }
    }
    let mut seen = [false; CORNER_SITES];
    for c in cells {
        for corner in 0..8usize {
            seen[site(
                c[0] + (corner & 1),
                c[1] + ((corner >> 1) & 1),
                c[2] + ((corner >> 2) & 1),
            )] = true;
        }
    }
    seen.iter().filter(|s| **s).count()
}

// ─── main ───────────────────────────────────────────────────────────────────

fn main() {
    if !std::env::args().any(|arg| arg == "--bench") {
        return;
    }
    let prereg = isomesh::experiment!("P-126");
    let mut breaches: Vec<(&'static str, u32, String)> = Vec::new();

    common::experiment::run(prereg, |run| {
        // The block's arithmetic, asserted rather than trusted.
        assert_eq!(BLOCK_CELLS, 27, "the block is 3 x 3 x 3 cells");
        assert_eq!(BLOCK_CORNERS, 64, "which is 4^3 = 64 corners as drawn");
        assert_eq!(
            PATTERNS, 134_217_728,
            "and 2^27 patterns, because corner 3 is corner 0 on every axis"
        );

        let tables = CubeTables::new();
        let planes = PlaneTables::new();
        let stitched = Block::new([CELLS_PER_AXIS; 3], true);
        let open = Block::new([CELLS_PER_AXIS; 3], false);
        let grid_face = grid_face_table(&stitched);

        // Every cell of the stitched block has all 26 neighbours, which is the
        // whole reason for the period-3 identification; the open block has one.
        let complete_stitched = (0..BLOCK_CELLS)
            .filter(|&n| stitched.edge_available[n] == 0xFFF && stitched.face_available[n] == 0x3F)
            .count();
        let complete_open = (0..BLOCK_CELLS)
            .filter(|&n| open.edge_available[n] == 0xFFF && open.face_available[n] == 0x3F)
            .count();
        assert_eq!(
            complete_stitched, BLOCK_CELLS,
            "stitched, every one of the 27 cells must carry all twelve of its quads"
        );
        assert_eq!(
            complete_open, 1,
            "open, exactly the centre cell may — this is P-63's truncation at the right block size"
        );

        // C1's two lemmas, over their whole domains.
        let adjacency = cycle_adjacency_violations();
        let two_cut = two_cut_edge_face_splits(&tables);
        let disagreements = pairing_disagreements(&tables);
        let disconnected = face_disconnected_cases(&tables);
        assert_eq!(
            adjacency, 0,
            "a cycle's consecutive edges must share a face of the cell"
        );
        assert_eq!(
            two_cut, 0,
            "a face with two cut edges must put both in one component"
        );
        assert_eq!(
            disagreements, 0,
            "with the joined bit agreed, two cells must induce the same pairing on their shared \
             face — this is the property C1's proof rests on"
        );
        assert!(
            disconnected > 0,
            "VOID: no sign byte has a face-disconnected cut set, so the one-vertex-per-cell \
             control cannot fire and C1's ones would prove nothing"
        );

        let expected_active = 27 * (1u64 << 19) * 254;
        println!(
            "\nblock: {BLOCK_CELLS} cells, {BLOCK_CORNERS} corners as drawn, {CORNER_SITES} \
             distinct under period 3, {PATTERNS} patterns"
        );
        println!(
            "lemmas: cycle adjacency violations {adjacency}, two-cut-edge splits {two_cut}, \
             pairing disagreements {disagreements}, face-disconnected sign bytes {disconnected} \
             of 256"
        );
        println!("expected active cells: {expected_active}");

        // The redundancy the period-3 identification buys, reported rather than
        // hidden: one cell's link depends on its own eight corner sites and its
        // six face-neighbours', which under the wrap is 20 of the 27 — so each
        // distinct link configuration is visited 2^(27-20) = 128 times, and the
        // sweep is exhaustive with a multiplicity rather than exhaustive over
        // more.
        let determining = link_determining_sites();
        assert_eq!(
            determining, 20,
            "a cell and its six face-neighbours span 20 of the block's 27 corner sites"
        );
        println!(
            "sweep: a cell's link is decided by {determining} of {CORNER_SITES} sites, so each \
             distinct configuration is visited {} times; all {BLOCK_CELLS} cells are censused \
             because shifting a 3-periodic pattern is a bijection of the domain",
            1u32 << (CORNER_SITES - determining)
        );

        let replica = replica(&tables);
        assert_eq!(
            replica.decider_disagreements, 0,
            "the shipped asymptotic decider must agree across every interior grid face"
        );
        assert!(
            replica.link_agreements > 0,
            "VOID: the quad link and the triangulated link never agreed, so one of the two walks \
             is not measuring a link"
        );
        println!(
            "replica: {} patterns x {} shipped arms, {} vertices, {} triangles, bit identical {}, \
             decider face disagreements {}, quad-link agreements {} against {} diagonal merges",
            replica.patterns,
            replica.arms,
            replica.vertices,
            replica.triangles,
            replica.bit_identical,
            replica.decider_disagreements,
            replica.link_agreements,
            replica.diagonal_merges
        );

        println!(
            "\n{:<32} {:>4} {:>13} {:>13} {:>6} {:>13} {:>4} {:>9}",
            "arm", "ctrl", "vertices", "complete", "maxInc", "defective", "cmp", "s"
        );

        for arm in &ARMS {
            let block = if arm.stitched { &stitched } else { &open };
            let table = ArmTable::build(arm);
            let cycles = table.cycles_per_256_cases();
            let (c, wall) = sweep(arm, block, &tables, &planes, &grid_face, &table);

            // Populations derived from the block, before the numbers are read.
            assert_eq!(
                c.active_cells, expected_active,
                "{}: 27 cells x 2^19 x 254 active sign bytes is exact by construction",
                arm.name
            );
            let expected_vertices = cycles.map(|n| 27 * (1u64 << 19) * n);
            if let Some(expected) = expected_vertices {
                assert_eq!(
                    c.dual_vertices, expected,
                    "{}: this arm's joined mask is a function of the cell's own sign byte, so its \
                     dual-vertex count is an exact integer over the block",
                    arm.name
                );
            } else {
                assert!(
                    c.dual_vertices >= c.active_cells,
                    "{}: every active cell carries at least one component",
                    arm.name
                );
            }

            // The registered vacuity control: per arm, never pooled.
            assert!(
                c.complete > 0,
                "VOID: {} reported no dual vertex with a complete link, which is P-63's empty \
                 population arriving again",
                arm.name
            );
            if arm.stitched {
                assert_eq!(
                    c.incomplete, 0,
                    "{}: the stitched block truncates nothing",
                    arm.name
                );
            } else {
                assert!(
                    c.incomplete > c.complete,
                    "VOID: the open block must produce more truncated links than complete ones, \
                     or the completeness counter cannot report bad news"
                );
            }

            let expected_defective = if arm.placement == Placement::WholeCell && arm.stitched {
                Some(27 * (1u64 << 19) * disconnected)
            } else {
                None
            };
            if let Some(expected) = expected_defective {
                assert_eq!(
                    c.defective, expected,
                    "{}: one vertex per cell makes the link's components the face-connected \
                     components of the cut set, which is a function of the sign byte alone",
                    arm.name
                );
            }
            if arm.is_control {
                assert!(
                    c.worst > 1,
                    "VOID: control arm {} read worst_link_components {}, which is the P-63-C3 \
                     failure mode arriving again — a zero over an unreachable population — and \
                     aborts the run rather than recording a pass",
                    arm.name,
                    c.worst
                );
                assert!(
                    c.defective > 0,
                    "VOID: control arm {} produced no link-defective vertex",
                    arm.name
                );
            } else {
                breaches.push((
                    arm.name,
                    c.worst,
                    if c.worst > 1 {
                        format!("{:#09x}", c.worst_pattern)
                    } else {
                        String::from("none")
                    },
                ));
            }

            let worst_pattern = if c.worst > 1 {
                format!("{:#09x}", c.worst_pattern)
            } else {
                String::from("none")
            };
            let c1 = !arm.is_control && c.worst == 1;
            let c2 = arm.is_control && c.worst > 1 && c.defective > 0;
            let c3 = c.complete > 0
                && c.active_cells == expected_active
                && expected_vertices.is_none_or(|e| c.dual_vertices == e);

            println!(
                "{:<32} {:>4} {:>13} {:>13} {:>6} {:>13} {:>4} {:>9.1}",
                arm.name,
                arm.is_control,
                c.dual_vertices,
                c.complete,
                c.max_incident,
                c.defective,
                c.worst,
                wall
            );

            run.record(&[
                ("arm", arm.name.to_string()),
                ("block_cells", BLOCK_CELLS.to_string()),
                ("block_corners", BLOCK_CORNERS.to_string()),
                ("patterns_swept", PATTERNS.to_string()),
                ("dual_vertices", c.dual_vertices.to_string()),
                ("dual_vertices_with_complete_link", c.complete.to_string()),
                ("worst_link_components", c.worst.to_string()),
                ("worst_link_pattern", worst_pattern),
                ("max_incident_faces", c.max_incident.to_string()),
                ("link_defective_vertices", c.defective.to_string()),
                ("is_control", arm.is_control.to_string()),
                ("wall_seconds", format!("{wall:.3}")),
                ("c1_holds", c1.to_string()),
                ("c2_holds", c2.to_string()),
                ("c3_holds", c3.to_string()),
                // Extras (M-273): the rule the row is about, the closed forms it
                // was checked against, and the provenance of the instrument.
                (
                    "boundary",
                    String::from(if arm.stitched { "periodic" } else { "open" }),
                ),
                ("vertex_rule", arm.placement_name().to_string()),
                ("ambiguous_face_pairing", arm.pairing_name().to_string()),
                (
                    "pairing_is_consistent",
                    arm.pairing_is_consistent().to_string(),
                ),
                ("corner_identification", String::from("period_3_per_axis")),
                ("distinct_corner_sites", CORNER_SITES.to_string()),
                ("active_cells", c.active_cells.to_string()),
                ("expected_active_cells", expected_active.to_string()),
                (
                    "cycles_per_256_sign_bytes",
                    cycles.map_or_else(|| String::from("unavailable"), |n| n.to_string()),
                ),
                (
                    "expected_dual_vertices",
                    expected_vertices
                        .map_or_else(|| String::from("unavailable"), |n| n.to_string()),
                ),
                (
                    "expected_link_defective_vertices",
                    expected_defective
                        .map_or_else(|| String::from("unavailable"), |n| n.to_string()),
                ),
                ("incomplete_link_vertices", c.incomplete.to_string()),
                ("isolated_vertices", c.isolated.to_string()),
                ("worst_link_cell", c.worst_cell.to_string()),
                ("worst_link_case", format!("{:#04x}", c.worst_case)),
                ("cycle_adjacency_violations", adjacency.to_string()),
                ("two_cut_edge_face_splits", two_cut.to_string()),
                ("pairing_disagreements", disagreements.to_string()),
                ("face_disconnected_sign_bytes", disconnected.to_string()),
                ("replica_patterns", replica.patterns.to_string()),
                ("replica_shipped_arms", replica.arms.to_string()),
                ("replica_bit_identical", replica.bit_identical.to_string()),
                ("replica_vertices", replica.vertices.to_string()),
                ("replica_triangles", replica.triangles.to_string()),
                (
                    "replica_decider_face_disagreements",
                    replica.decider_disagreements.to_string(),
                ),
                (
                    "replica_quad_link_agreements",
                    replica.link_agreements.to_string(),
                ),
                (
                    "replica_diagonal_merges",
                    replica.diagonal_merges.to_string(),
                ),
                (
                    "target_feature_popcnt",
                    cfg!(target_feature = "popcnt").to_string(),
                ),
            ]);
        }
    });

    // **The gate, after the CSV.** The bench is the single implementation and
    // there is no second path: `nightly.yml` runs it and reports this exit code.
    let failed: Vec<&(&str, u32, String)> = breaches.iter().filter(|(_, w, _)| *w > 1).collect();
    if failed.is_empty() {
        println!(
            "\nP-126 GATE: PASS — every non-control arm read worst_link_components 1 over all \
             {PATTERNS} patterns of the 3 x 3 x 3 block."
        );
        return;
    }
    for (arm, worst, pattern) in &failed {
        eprintln!(
            "P-126 GATE: FAIL — non-control arm `{arm}` read worst_link_components {worst} at \
             pattern {pattern}. The dual output is non-manifold at a vertex on a reachable sign \
             pattern; A-017 documents MDC's non-manifoldness as a limit of the algorithm and this \
             names the configuration."
        );
    }
    std::process::exit(1);
}
