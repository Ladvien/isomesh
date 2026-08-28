//! **P-119 — double-buffering as the determinism mechanism.**
//!
//! Ticket: R-119. Pre-registered before this harness existed.
//!
//! ```bash
//! cargo bench --bench experiment_p119
//! ```
//!
//! Writes `docs/experiments/p-119.csv`. **Linux only**, for `experiment_p12`'s
//! reason: C2 is a cost bound and `M-280`/`M-281` put the verdict on the
//! instruction form, whose only instrument is `perf_event_open`. Off Linux
//! there is nothing to degrade to and a recorded zero would be a fabrication.
//!
//! # What was missing
//!
//! Schroeder, Maynard, Geveci et al., *A High-Performance SurfaceNets Discrete
//! Isocontouring Algorithm* (`10.48550/arXiv.2401.14906`, in the corpus),
//! double-buffers the smoothing pass on purpose: every iteration reads the
//! *previous* position array and writes a *different* one, so a vertex's new
//! position cannot depend on whether its neighbours have already been updated
//! in this pass. That is a determinism mechanism, not an optimisation, and it
//! is the only structural device in the sweep aimed squarely at the property
//! this crate currently gets by accident.
//!
//! **Their one-to-two-orders-of-magnitude figure is a *parallel* result.** It
//! is speedup from threading a smoothing pass that double-buffering made safe
//! to thread. Every loop in `crates/isomesh/src/**` is sequential. **That figure
//! is therefore not a comparand and nothing here is scored against it.** Only
//! the *structure* transfers: snapshot reads, no read-your-own-writes.
//!
//! ## What P-9 found, and where the defect lives
//!
//! `P-9` (`R-002`, `docs/experiments/p-9.csv`) permuted within-bucket merge
//! order on 56 configurations and found nine that produce more than one
//! distinct output. Two of them move the vertex **count**: `noise_cavity` under
//! `dual_contouring` spans **4** vertices across eight permutations and under
//! `manifold_dual_contouring` **2** (`FINDINGS.md:1641-1642`, and the rows are
//! in the committed CSV, which is what [`P9_PUBLISHED`] holds).
//!
//! **The implicated path is the weld, and the registration says so.**
//! `weld.rs:232` `Welder`, `:296` `weld`, `:338` `weld_split_by`, `:220`
//! `epsilon_for`. Epsilon-closeness is **not transitive**, so "weld everything
//! within ε" does not define equivalence classes; what `weld.rs:39-56` defines
//! instead is **first fit** — walk vertices in input index order, join the
//! lowest-indexed *representative* within ε, or become one. The inner test at
//! `weld.rs:424` reads `self.kept[u]`, the **partially built output of this
//! very pass**. That is read-your-own-writes, and it is exactly the dependence
//! Parallel SurfaceNets' double buffer removes. A chain `a ~ b ~ c` with
//! `a`–`c` further than ε apart elects one representative or two depending on
//! which end the walk reaches first, so permuting the input changes **how
//! many** representatives are elected, not merely which.
//!
//! That is why C1's falsifier points at the weld rather than at the traversal:
//! if double-buffering the weld does not remove the spread, the order
//! dependence is somewhere else and the work goes back to `weld.rs`.
//!
//! # The three arms, and why there are three rather than two
//!
//! Every arm is **bench-local** — `crates/isomesh/src/**` is read-only for
//! Phase 25 — and all three share one broadphase ([`Lattice`], [`sorted_cells`])
//! so the only thing that differs between them is the election.
//!
//! - **`single`** — [`SingleBuffer`], a line-by-line mirror of
//!   `Welder::weld_split_by` with empty keys: same lattice, same 27-cell probe,
//!   same "lowest-indexed representative", same `kept` read-your-own-writes,
//!   same in-place compaction and the same degenerate-triangle drop. It exists
//!   because the shipped `Welder`'s scratch capacities are private and
//!   `peak_bytes_single` is a registered column, and because two arms compared
//!   for cost must be structurally identical except for the mechanism.
//!
//!   **The mirror is proved against the original in the same run.** Column
//!   `mirror_matches_shipped` compares positions, normals and indices as bit
//!   patterns *and* the full remap against `isomesh::weld::Welder` on every
//!   row, and it is asserted, not merely recorded. `R-120` and `R-121` both
//!   caught real defects that way, and it is what licenses reading a cost or a
//!   spread off a mirror at all (`M-279`).
//!
//! - **`double_indexed`** — [`DoubleBuffer`] with `Order::InputIndex`. This is
//!   double-buffering *and nothing else*. Pass A reads the snapshot and writes
//!   a CSR ε-adjacency; pass B propagates labels `next[v] = min(cur[v],
//!   min_{u~v} cur[u])` reading `cur` and writing `next`, swapping between
//!   iterations, to a fixed point. Nothing ever reads what this pass wrote.
//!   The label is the component's minimum **input index**.
//!
//! - **`double`** — [`DoubleBuffer`] with `Order::Canonical`. Identical, except
//!   that "minimum" is over the vertex's own `(position, normal)` bit pattern
//!   rather than its buffer slot, and the output vertices are emitted in that
//!   same order.
//!
//! **The third arm is the finding, not padding.** Double-buffering alone makes
//! the *partition* a function of the ε-graph, which is a set of unordered
//! pairs — so the component count, and therefore the vertex count, is
//! order-free by construction and C1 is settled by `double_indexed` already.
//! It does **not** make the *mesh* order-free: "minimum input index" is a
//! statement about the buffer, so which member of a component survives, and in
//! what slot, still moves. C3 needs the canonical tie-break, which is the
//! second half of the mechanism and is not double-buffering. Columns
//! `spread_double_indexed` and `mesh_identical_across_permutations_indexed`
//! separate the two halves, and a reader can see which clause each bought.
//!
//! ## What the double buffer costs, stated rather than hidden
//!
//! Label propagation to a fixed point over the ε-graph converges to the
//! **connected components** of that graph — that is, to the transitive closure
//! of ε-closeness. `weld.rs:41-52` refuses the transitive closure on purpose,
//! for a stated reason: under first fit "no vertex is ever moved further than ε
//! from where the extractor put it", and a chain can drag one arbitrarily far.
//!
//! **The order dependence and the non-transitivity are the same property.**
//! First fit is order-dependent *because* it breaks chains, and any rule that
//! is a function of the ε-graph alone cannot break them. So this is not a
//! deterministic version of the shipped weld; it is a different weld that is
//! deterministic. Three columns price that rather than leaving it as prose:
//! `vertices_single`, `vertices_double` and `vertex_delta` say how much the
//! answer moves, and `max_move_over_epsilon` — the furthest any vertex is
//! carried from where the extractor put it, in units of ε — says whether the
//! chain-dragging `weld.rs` warns about is real on this fixture or theoretical.
//!
//! # SHARE
//!
//! Each clause's reachable share, as a column.
//!
//! - **C1's share is `baseline_spread`**, the single-buffer arm's vertex-count
//!   spread, which is the entire quantity available to remove. It is the
//!   registered **VACUITY CONTROL** and it is checked rather than hoped: a
//!   harness that cannot see `P-9`'s defect cannot measure its removal. Two
//!   perturbation families are run, because the registration names one and the
//!   number 4/2 came from the other:
//!   - `baseline_spread` is over **chunk append order**, the registered fixture
//!     and the live consumer scenario `P-9`'s own conclusion names.
//!   - `baseline_spread_within_bucket` is `P-9`'s **exact** protocol —
//!     within-bucket permutation, the shipped `Welder`, the same eight seeds —
//!     scored against [`P9_PUBLISHED`], every one of the sixteen rows of
//!     `docs/experiments/p-9.csv` this fixture covers, not merely the two that
//!     move. `p9_reproduced` is that comparison and it is asserted.
//!
//!   `c1_holds` reads `vacuous` on a row whose control did not fire, never
//!   `true`. A zero spread removed from a zero spread is not a measurement.
//!
//! - **C2's share is `instructions_per_cell_single`**, the shipped weld's
//!   instruction rate over the 46,656 cells of the eight-chunk block, and the
//!   bar is a ratio of 1.25. **The verdict reads the instruction form**
//!   (`instruction_cost_ratio`), not the clock: `M-280`, `M-281`, and `R-105`
//!   watching one binary's cycle band drift 0.984 → 1.035 across three runs
//!   while its instruction counts held to four figures. `ns_per_cell_single`,
//!   `ns_per_cell_double`, `cost_ratio` and `ghz` are on every row as the
//!   registration requires, `cycles_per_cell_*` beside them, and `c2_holds_ns`
//!   records what the clock would have said so the two can be compared.
//!
//!   **The window is the weld, not the weld plus the append.** Both arms are
//!   preceded by the identical `restore` of a pristine buffer, so a third
//!   sibling window measures `restore` alone and the reported figures are the
//!   **prefix difference**. Diluting the ratio with work common to both arms
//!   can only push it toward 1.0, which is the direction that would make a
//!   1.25 bar hold by arithmetic rather than by measurement;
//!   `cost_ratio_undiluted` is the honest reading and `cost_ratio_with_restore`
//!   is beside it so the size of the dilution is visible.
//!
//!   Windows are **siblings, never nested**. Zen 3 has six general-purpose
//!   counters and `Probe` opens six, so a nested pair multiplexes and
//!   `Probe::worst_ratio` refuses — `R-121` paid for that discovery. The three
//!   windows rotate their order by repetition index so no arm is permanently
//!   the one that inherits another's cache state.
//!
//! - **C3's share is the mesh**, and it is an equality rather than a fraction.
//!   Comparing raw buffers across a permutation of chunk append order would be
//!   trivially false for every arm — the vertex array *is* permuted — and would
//!   measure nothing. So the comparison is the **canonical** mesh
//!   ([`canonical_digest`]): the multiset of `(position, normal)` bit patterns,
//!   and the multiset of triangles written as canonically-ranked corner
//!   triples, rotated to put the smallest corner first so winding survives.
//!   The identical comparison is applied to all three arms, and
//!   `mesh_identical_across_permutations_single` is the control that gives the
//!   clause meaning.
//!
//! - **No popcount is involved on either side.** This repository emits zero
//!   `popcnt` — no `.cargo/config.toml`, no `target-cpu`, so
//!   `cfg!(target_feature = "popcnt")` is false and `u64::count_ones()` lowers
//!   to a ~12-instruction SWAR sequence. Neither arm is a bitmap and neither
//!   calls `count_ones`, so no verdict here is contingent on that lowering.
//!   `target_feature_popcnt` is a column so a reader comparing this row against
//!   a rank/select row can see in one cell which of them pays the SWAR tax.

#![allow(
    clippy::cast_possible_truncation,
    clippy::cast_precision_loss,
    clippy::too_many_lines
)]

mod common;

#[cfg(target_os = "linux")]
mod experiment {
    use std::collections::{BTreeMap, BTreeSet};
    use std::hint::black_box;
    use std::time::Instant;

    use isomesh::chunk::{ChunkId, ChunkLayout};
    use isomesh::extractor::Extractor;
    use isomesh::weld::{Welder, epsilon_for};
    use isomesh::{MeshBuffer, Sdf};

    use crate::common::counters::{MIN_TIME_RATIO, Probe};

    // ─── the fixture, which is P-9's own ───────────────────────────────────

    /// Cells per chunk. `P-9`'s value, and `P-8`'s before it (`M-274`).
    const CELLS: u32 = 18;
    /// `P-9`'s cell size. `4/35` because a seam is bit-exact only at a power of
    /// two (`M-32`), so this is the size at which a seam bucket exists at all.
    const CELL_SIZE: f64 = 4.0 / 35.0;
    /// Block origin, placed so the 2×2×2 block is centred on the field.
    const ORIGIN: f64 = -(2.0 * CELLS as f64) * CELL_SIZE / 2.0;
    /// Chunks in the block.
    const CHUNKS: usize = 8;
    /// Permutations per row, and **`P-9`'s eight**: the vacuity control is a
    /// comparison against `P-9`'s numbers, and a spread is monotone
    /// non-decreasing in the permutation count, so a different count would not
    /// be the same measurement.
    const PERMUTATIONS: u32 = 8;
    /// A second, larger chunk-order sweep, recorded beside the registered one.
    ///
    /// Not a substitute for it. A spread that is already saturated at eight
    /// permutations and a spread that merely has not been sampled enough look
    /// identical from one number, and `baseline_spread_extended` is what tells
    /// them apart.
    const EXTENDED: u32 = 32;
    /// Cells in the block: eight chunks of `CELLS³`. The denominator of every
    /// per-cell column.
    const CELLS_TOTAL: usize = CHUNKS * (CELLS as usize).pow(3);

    /// Counted repetitions per row. A multiple of three, because the three
    /// sibling windows rotate their order by repetition index and an
    /// off-multiple count would give one ordering an extra vote.
    const REPS: usize = 6;
    /// Target wall time for one counter window, in nanoseconds.
    const TARGET_BATCH_NS: f64 = 30_000_000.0;
    /// Ceiling on the inner repetition count.
    const MAX_INNER: usize = 4096;

    /// `docs/experiments/p-9.csv`, restricted to the sixteen rows this fixture
    /// covers: `(field, extractor, distinct_outputs, vertex_count_spread,
    /// buckets_of_three_or_more)`.
    ///
    /// Transcribed from the committed dataset rather than from the prose,
    /// because the prose quotes a seven-row selection and the control is
    /// stronger over all sixteen — and over all **three** of `P-9`'s recorded
    /// columns rather than the one the registration names.
    /// `FINDINGS.md:1642-1643` is the `noise_cavity` pair — **4** and **2** —
    /// which is the two numbers the registration's VACUITY CONTROL names.
    const P9_PUBLISHED: [(&str, &str, usize, usize, usize); 16] = [
        ("sphere", "dual_contouring", 1, 0, 0),
        ("sphere", "manifold_dual_contouring", 1, 0, 0),
        ("torus", "dual_contouring", 5, 0, 0),
        ("torus", "manifold_dual_contouring", 5, 0, 0),
        ("box_exact", "dual_contouring", 1, 0, 0),
        ("box_exact", "manifold_dual_contouring", 1, 0, 0),
        ("csg_difference", "dual_contouring", 1, 0, 0),
        ("csg_difference", "manifold_dual_contouring", 1, 0, 0),
        ("thin_plate", "dual_contouring", 1, 0, 0),
        ("thin_plate", "manifold_dual_contouring", 1, 0, 0),
        ("gyroid", "dual_contouring", 1, 0, 0),
        ("gyroid", "manifold_dual_contouring", 1, 0, 0),
        ("fbm_terrain", "dual_contouring", 1, 0, 0),
        ("fbm_terrain", "manifold_dual_contouring", 1, 0, 0),
        ("noise_cavity", "dual_contouring", 8, 4, 1),
        ("noise_cavity", "manifold_dual_contouring", 8, 2, 2),
    ];

    /// What `P-9` published for one configuration.
    fn published(field: &str, extractor: &str) -> (usize, usize, usize) {
        P9_PUBLISHED
            .iter()
            .find(|(f, e, ..)| *f == field && *e == extractor)
            .map(|&(_, _, distinct, spread, big)| (distinct, spread, big))
            .expect("every row of this fixture is a row of docs/experiments/p-9.csv")
    }

    // ─── the reproducible shuffle, copied from `experiment_p9` ─────────────

    /// A seeded xorshift64\*, byte-for-byte `experiment_p9`'s.
    ///
    /// Copied rather than shared because the two harnesses must draw the *same*
    /// sequence from the same seed or the within-bucket arm is not `P-9`'s
    /// experiment, and a shared bench helper that later drifts would break that
    /// silently.
    struct Rng(u64);

    impl Rng {
        const fn new(seed: u64) -> Self {
            Self(seed | 1)
        }

        fn next(&mut self) -> u64 {
            let mut x = self.0;
            x ^= x >> 12;
            x ^= x << 25;
            x ^= x >> 27;
            self.0 = x;
            x.wrapping_mul(0x2545_F491_4F6C_DD1D)
        }

        /// Fisher–Yates, so every ordering is reachable.
        fn shuffle<T>(&mut self, slice: &mut [T]) {
            for i in (1..slice.len()).rev() {
                let j = (self.next() % (i as u64 + 1)) as usize;
                slice.swap(i, j);
            }
        }
    }

    /// `P-9`'s seed schedule, so seed `k` here is seed `k` there.
    fn seeded(seed: u32) -> Rng {
        Rng::new(u64::from(seed).wrapping_mul(0x9E37_79B9_7F4A_7C15) + 1)
    }

    // ─── the broadphase, mirrored from `weld.rs` ───────────────────────────

    /// `weld::quantise` for `f64`, including the narrowing through `f32`.
    ///
    /// The narrowing is not incidental — `weld.rs:88-101` explains that it is
    /// what forces the lattice anchor to be the mesh's own minimum — so a
    /// mirror that widened it would bucket differently at the margins and the
    /// `mirror_matches_shipped` assertion would fire.
    fn quantise(scaled: f64) -> i64 {
        let f = scaled.floor();
        if f.is_finite() { f as f32 as i64 } else { 0 }
    }

    /// `weld::Lattice` for `f64`.
    struct Lattice {
        anchor: [f64; 3],
        inv_epsilon: f64,
    }

    impl Lattice {
        fn new(positions: &[[f64; 3]], epsilon: f64) -> Self {
            let mut anchor = [f64::INFINITY; 3];
            for p in positions {
                for (a, slot) in anchor.iter_mut().enumerate() {
                    if p[a].is_finite() && p[a] < *slot {
                        *slot = p[a];
                    }
                }
            }
            for slot in &mut anchor {
                if !slot.is_finite() {
                    *slot = 0.0;
                }
            }
            Self {
                anchor,
                inv_epsilon: epsilon.recip(),
            }
        }

        fn key_of(&self, p: [f64; 3]) -> [i64; 3] {
            [
                quantise((p[0] - self.anchor[0]) * self.inv_epsilon),
                quantise((p[1] - self.anchor[1]) * self.inv_epsilon),
                quantise((p[2] - self.anchor[2]) * self.inv_epsilon),
            ]
        }
    }

    /// The sorted `(lattice cell, vertex)` broadphase, shared by all three arms.
    ///
    /// The vertex index is part of the key and is unique, so no two entries
    /// compare equal and an unstable sort is a deterministic one —
    /// `weld.rs:394-396`'s reasoning, kept because the mirror must reproduce
    /// its result exactly.
    fn sorted_cells(cells: &mut Vec<([i64; 3], u32)>, lattice: &Lattice, positions: &[[f64; 3]]) {
        cells.clear();
        cells.reserve(positions.len());
        for (i, p) in positions.iter().enumerate() {
            cells.push((lattice.key_of(*p), i as u32));
        }
        cells.sort_unstable();
    }

    /// Squared distance, in the order `weld.rs:437-438` writes it.
    fn dist_sq(p: [f64; 3], q: [f64; 3]) -> f64 {
        let d = [p[0] - q[0], p[1] - q[1], p[2] - q[2]];
        d[0] * d[0] + d[1] * d[1] + d[2] * d[2]
    }

    // ─── arm one: the single-buffer mirror ─────────────────────────────────

    /// A line-by-line mirror of `Welder::weld_split_by` with empty keys.
    ///
    /// Owns its scratch for the same reason the original does, and exposes
    /// [`scratch_bytes`](Self::scratch_bytes) because the original's capacities
    /// are private and `peak_bytes_single` is a registered column.
    #[derive(Default)]
    struct SingleBuffer {
        cells: Vec<([i64; 3], u32)>,
        remap: Vec<u32>,
        kept: Vec<bool>,
    }

    impl SingleBuffer {
        /// `weld.rs:338-498`, for `f64` and with no split keys.
        fn weld(&mut self, mesh: &mut MeshBuffer<f64>, epsilon: f64) {
            let n = mesh.positions.len();
            // `weld.rs:367-375`. Kept because it is real work on the shipped
            // path and dropping it would make the mirror cheaper than the
            // thing it stands in for.
            for &i in &mesh.indices {
                assert!((i as usize) < n, "index out of range");
            }
            let triangles_before = mesh.indices.len() / 3;
            if n == 0 {
                return;
            }

            let lattice = Lattice::new(&mesh.positions, epsilon);
            sorted_cells(&mut self.cells, &lattice, &mesh.positions);

            self.remap.clear();
            self.remap.resize(n, u32::MAX);
            self.kept.clear();
            self.kept.resize(n, false);

            let eps_sq = epsilon * epsilon;
            let mut next_output = 0u32;

            for v in 0..n {
                let p = mesh.positions[v];
                let base = lattice.key_of(p);

                // The *lowest-indexed* representative within epsilon, so the
                // answer does not depend on the order the 27 cells are visited
                // in — and `self.kept` is the output of this same pass, which
                // is the read-your-own-writes the double buffer removes.
                let mut best: Option<u32> = None;
                for dz in -1..=1i64 {
                    for dy in -1..=1i64 {
                        for dx in -1..=1i64 {
                            let key = [base[0] + dx, base[1] + dy, base[2] + dz];
                            let lo = self.cells.partition_point(|(k, _)| *k < key);
                            for &(k, u) in &self.cells[lo..] {
                                if k != key {
                                    break;
                                }
                                if u as usize >= v || !self.kept[u as usize] {
                                    continue;
                                }
                                if best.is_some_and(|b| u >= b) {
                                    continue;
                                }
                                if dist_sq(p, mesh.positions[u as usize]) <= eps_sq {
                                    best = Some(u);
                                }
                            }
                        }
                    }
                }

                match best {
                    Some(u) => self.remap[v] = self.remap[u as usize],
                    None => {
                        self.kept[v] = true;
                        self.remap[v] = next_output;
                        next_output += 1;
                    }
                }
            }

            if next_output as usize == n {
                // `weld.rs:457-462`: nothing coincided, so no index moved.
                return;
            }

            let mut w = 0usize;
            for v in 0..n {
                if self.kept[v] {
                    mesh.positions[w] = mesh.positions[v];
                    mesh.normals[w] = mesh.normals[v];
                    w += 1;
                }
            }
            mesh.positions.truncate(w);
            mesh.normals.truncate(w);

            let mut t = 0usize;
            for tri in 0..triangles_before {
                let a = self.remap[mesh.indices[tri * 3] as usize];
                let b = self.remap[mesh.indices[tri * 3 + 1] as usize];
                let c = self.remap[mesh.indices[tri * 3 + 2] as usize];
                if a == b || b == c || a == c {
                    continue;
                }
                mesh.indices[t * 3] = a;
                mesh.indices[t * 3 + 1] = b;
                mesh.indices[t * 3 + 2] = c;
                t += 1;
            }
            mesh.indices.truncate(t * 3);
        }

        /// Peak scratch, from the real capacities rather than from arithmetic.
        fn scratch_bytes(&self) -> usize {
            self.cells.capacity() * size_of::<([i64; 3], u32)>()
                + self.remap.capacity() * size_of::<u32>()
                + self.kept.capacity()
        }
    }

    // ─── arm two and three: the double buffer ──────────────────────────────

    /// Which total order the component's surviving member is the minimum of.
    #[derive(Clone, Copy, PartialEq, Eq)]
    enum Order {
        /// The buffer slot. Double-buffering and nothing else.
        InputIndex,
        /// The vertex's own `(position, normal)` bit pattern.
        Canonical,
    }

    /// The double-buffered weld.
    ///
    /// Pass A reads the snapshot and writes the ε-adjacency; pass B iterates
    /// `next[v] = min(cur[v], min over neighbours of cur[u])`, reading `cur`
    /// and writing `next` and swapping between iterations, to a fixed point.
    /// No pass ever reads a value it wrote.
    #[derive(Default)]
    struct DoubleBuffer {
        cells: Vec<([i64; 3], u32)>,
        /// CSR row starts, `n + 1` long.
        offsets: Vec<u32>,
        /// CSR neighbour lists.
        neighbours: Vec<u32>,
        /// The two label buffers. `cur` is read, `next` is written.
        cur: Vec<u32>,
        next: Vec<u32>,
        /// Output slot of each representative, `u32::MAX` for a non-root.
        out_index: Vec<u32>,
        /// Output slot of every input vertex.
        remap: Vec<u32>,
        /// Representatives, in output order.
        roots: Vec<u32>,
        out_positions: Vec<[f64; 3]>,
        out_normals: Vec<[f64; 3]>,
    }

    /// What one double-buffered weld did, beyond the mesh it rewrote.
    #[derive(Clone, Copy, Default)]
    struct DoubleReport {
        /// Label-propagation passes to the fixed point, including the pass
        /// that changed nothing. Two on a graph of isolated pairs.
        iterations: u32,
        /// The furthest any vertex was carried from its own position, in units
        /// of ε. Above 1.0 is the chain-dragging `weld.rs:41-52` warns about.
        max_move_over_epsilon: f64,
        /// Vertices in components of three or more, which is where first fit
        /// and the transitive closure can disagree.
        vertices_in_large_components: usize,
    }

    impl DoubleBuffer {
        /// `true` when `a` should survive in preference to `b`.
        fn beats(order: Order, mesh: &MeshBuffer<f64>, a: u32, b: u32) -> bool {
            match order {
                Order::InputIndex => a < b,
                Order::Canonical => vertex_key(mesh, a) < vertex_key(mesh, b),
            }
        }

        fn weld(&mut self, mesh: &mut MeshBuffer<f64>, epsilon: f64, order: Order) -> DoubleReport {
            let n = mesh.positions.len();
            for &i in &mesh.indices {
                assert!((i as usize) < n, "index out of range");
            }
            let triangles_before = mesh.indices.len() / 3;
            if n == 0 {
                return DoubleReport::default();
            }

            let lattice = Lattice::new(&mesh.positions, epsilon);
            sorted_cells(&mut self.cells, &lattice, &mesh.positions);
            let eps_sq = epsilon * epsilon;

            // ── pass A: the snapshot, read-only, written to the CSR ─────────
            self.offsets.clear();
            self.offsets.push(0);
            self.neighbours.clear();
            for v in 0..n {
                let p = mesh.positions[v];
                let base = lattice.key_of(p);
                for dz in -1..=1i64 {
                    for dy in -1..=1i64 {
                        for dx in -1..=1i64 {
                            let key = [base[0] + dx, base[1] + dy, base[2] + dz];
                            let lo = self.cells.partition_point(|(k, _)| *k < key);
                            for &(k, u) in &self.cells[lo..] {
                                if k != key {
                                    break;
                                }
                                if u as usize == v {
                                    continue;
                                }
                                if dist_sq(p, mesh.positions[u as usize]) <= eps_sq {
                                    self.neighbours.push(u);
                                }
                            }
                        }
                    }
                }
                self.offsets.push(self.neighbours.len() as u32);
            }

            // ── pass B: double-buffered label propagation ───────────────────
            self.cur.clear();
            self.cur.extend(0..n as u32);
            self.next.clear();
            self.next.resize(n, 0);

            let mut iterations = 0u32;
            loop {
                iterations += 1;
                let mut changed = false;
                for v in 0..n {
                    // Every read is from `cur`; every write is to `next`. That
                    // is the whole mechanism.
                    let mut best = self.cur[v];
                    let lo = self.offsets[v] as usize;
                    let hi = self.offsets[v + 1] as usize;
                    for &u in &self.neighbours[lo..hi] {
                        let candidate = self.cur[u as usize];
                        if Self::beats(order, mesh, candidate, best) {
                            best = candidate;
                        }
                    }
                    if best != self.cur[v] {
                        changed = true;
                    }
                    self.next[v] = best;
                }
                core::mem::swap(&mut self.cur, &mut self.next);
                if !changed {
                    break;
                }
            }

            // ── elect, in the order the rule names ──────────────────────────
            self.roots.clear();
            for v in 0..n {
                if self.cur[v] as usize == v {
                    self.roots.push(v as u32);
                }
            }
            if order == Order::Canonical {
                self.roots.sort_unstable_by_key(|&r| vertex_key(mesh, r));
            }

            self.out_index.clear();
            self.out_index.resize(n, u32::MAX);
            for (slot, &root) in self.roots.iter().enumerate() {
                self.out_index[root as usize] = slot as u32;
            }
            self.remap.clear();
            self.remap.resize(n, u32::MAX);
            let mut max_move_sq = 0.0f64;
            let mut in_large = 0usize;
            let mut members = vec![0u32; self.roots.len()];
            for v in 0..n {
                let root = self.cur[v] as usize;
                let slot = self.out_index[root];
                debug_assert_ne!(slot, u32::MAX, "every label is a root");
                self.remap[v] = slot;
                members[slot as usize] += 1;
                let d = dist_sq(mesh.positions[v], mesh.positions[root]);
                if d > max_move_sq {
                    max_move_sq = d;
                }
            }
            for v in 0..n {
                if members[self.remap[v] as usize] >= 3 {
                    in_large += 1;
                }
            }

            // ── write the output, in canonical (or index) order ─────────────
            self.out_positions.clear();
            self.out_normals.clear();
            for &root in &self.roots {
                self.out_positions.push(mesh.positions[root as usize]);
                self.out_normals.push(mesh.normals[root as usize]);
            }
            mesh.positions.clear();
            mesh.positions.extend_from_slice(&self.out_positions);
            mesh.normals.clear();
            mesh.normals.extend_from_slice(&self.out_normals);

            let mut t = 0usize;
            for tri in 0..triangles_before {
                let a = self.remap[mesh.indices[tri * 3] as usize];
                let b = self.remap[mesh.indices[tri * 3 + 1] as usize];
                let c = self.remap[mesh.indices[tri * 3 + 2] as usize];
                if a == b || b == c || a == c {
                    continue;
                }
                mesh.indices[t * 3] = a;
                mesh.indices[t * 3 + 1] = b;
                mesh.indices[t * 3 + 2] = c;
                t += 1;
            }
            mesh.indices.truncate(t * 3);

            DoubleReport {
                iterations,
                max_move_over_epsilon: max_move_sq.sqrt() / epsilon,
                vertices_in_large_components: in_large,
            }
        }

        /// Peak scratch, from the real capacities.
        fn scratch_bytes(&self) -> usize {
            self.cells.capacity() * size_of::<([i64; 3], u32)>()
                + (self.offsets.capacity()
                    + self.neighbours.capacity()
                    + self.cur.capacity()
                    + self.next.capacity()
                    + self.out_index.capacity()
                    + self.remap.capacity()
                    + self.roots.capacity())
                    * size_of::<u32>()
                + (self.out_positions.capacity() + self.out_normals.capacity())
                    * size_of::<[f64; 3]>()
        }
    }

    // ─── identity, raw and canonical ───────────────────────────────────────

    /// A vertex's `(position, normal)` bit pattern, as a total order.
    ///
    /// Bit patterns rather than floats: this is an identity, and `float_cmp` is
    /// the right lint to obey here rather than silence. Any fixed total order
    /// works so long as it is a function of the value alone, which `to_bits`
    /// is — including for `-0.0`, which it distinguishes, as byte-identity
    /// must.
    fn vertex_key(mesh: &MeshBuffer<f64>, v: u32) -> [u64; 6] {
        let p = mesh.positions[v as usize];
        let q = mesh.normals[v as usize];
        [
            p[0].to_bits(),
            p[1].to_bits(),
            p[2].to_bits(),
            q[0].to_bits(),
            q[1].to_bits(),
            q[2].to_bits(),
        ]
    }

    /// FNV-1a, `experiment_p9`'s.
    struct Fnv(u64);

    impl Fnv {
        const fn new() -> Self {
            Self(0xcbf2_9ce4_8422_2325)
        }
        fn eat(&mut self, word: u64) {
            for b in word.to_le_bytes() {
                self.0 ^= u64::from(b);
                self.0 = self.0.wrapping_mul(0x1000_0000_01b3);
            }
        }
    }

    /// `experiment_p9`'s byte digest: the buffers exactly as they lie.
    ///
    /// Kept because the within-bucket arm has to reproduce `P-9`'s
    /// `distinct_outputs`, and that number is defined by *this* comparison.
    fn raw_digest(mesh: &MeshBuffer<f64>) -> u64 {
        let mut h = Fnv::new();
        for p in &mesh.positions {
            for c in p {
                h.eat(c.to_bits());
            }
        }
        for nrm in &mesh.normals {
            for c in nrm {
                h.eat(c.to_bits());
            }
        }
        for i in &mesh.indices {
            h.eat(u64::from(*i));
        }
        h.0
    }

    /// The mesh as a geometric object, independent of buffer labelling.
    ///
    /// Under a permutation of chunk append order the vertex array *is*
    /// permuted, so [`raw_digest`] would differ for every arm and would measure
    /// nothing. This hashes the sorted multiset of vertex keys and the sorted
    /// multiset of triangles written as canonically-ranked corner triples,
    /// rotated so the smallest corner is first — which preserves winding, so
    /// two triangles that differ only in orientation still differ here.
    ///
    /// # Panics
    ///
    /// If two output vertices share a key. Both welds elect representatives
    /// that are pairwise more than ε apart, so a collision means one of them
    /// is broken and the ranking below would be ambiguous — which is exactly
    /// the kind of thing that must fail loudly rather than hash to something
    /// plausible.
    fn canonical_digest(mesh: &MeshBuffer<f64>) -> u64 {
        let n = mesh.positions.len();
        let mut rank_of: BTreeMap<[u64; 6], u32> = BTreeMap::new();
        for v in 0..n as u32 {
            assert!(
                rank_of.insert(vertex_key(mesh, v), 0).is_none(),
                "two output vertices share a position and normal, so the canonical rank is \
                 ambiguous and the weld that produced them is wrong"
            );
        }
        for (slot, value) in rank_of.values_mut().enumerate() {
            *value = slot as u32;
        }
        let mut rank = vec![0u32; n];
        for v in 0..n as u32 {
            rank[v as usize] = rank_of[&vertex_key(mesh, v)];
        }

        let mut triangles: Vec<[u32; 3]> = Vec::with_capacity(mesh.indices.len() / 3);
        for tri in mesh.indices.as_chunks::<3>().0 {
            let mut t = [
                rank[tri[0] as usize],
                rank[tri[1] as usize],
                rank[tri[2] as usize],
            ];
            // Rotate, never sort: a sort would identify a triangle with its
            // mirror image and C3 would stop being able to see a flipped
            // winding.
            let smallest = (0..3).min_by_key(|&i| t[i]).expect("three corners");
            t.rotate_left(smallest);
            triangles.push(t);
        }
        triangles.sort_unstable();

        let mut h = Fnv::new();
        for key in rank_of.keys() {
            for word in key {
                h.eat(*word);
            }
        }
        for t in &triangles {
            for c in t {
                h.eat(u64::from(*c));
            }
        }
        h.0
    }

    // ─── the block, and the two perturbation families ──────────────────────

    /// The eight chunk meshes, extracted once and appended many times.
    fn eight_pieces<E: Extractor<f64>>(
        field: &impl Sdf<Scalar = f64>,
        layout: &ChunkLayout<f64>,
        extractor: &mut E,
    ) -> Vec<MeshBuffer<f64>> {
        let shape = layout.sample_shape().expect("valid shape");
        let mut pieces = Vec::with_capacity(CHUNKS);
        for z in 0..2 {
            for y in 0..2 {
                for x in 0..2 {
                    let mut piece = MeshBuffer::<f64>::new();
                    extractor
                        .extract_into(
                            field,
                            &shape,
                            layout.sample_origin(ChunkId::new([x, y, z])),
                            layout.cell_size(),
                            &mut piece,
                        )
                        .expect("extraction");
                    pieces.push(piece);
                }
            }
        }
        pieces
    }

    /// Append the eight pieces in `order`.
    fn joined(pieces: &[MeshBuffer<f64>], order: &[usize]) -> MeshBuffer<f64> {
        let mut out = MeshBuffer::<f64>::new();
        for &c in order {
            out.append(&pieces[c]).expect("the meshes fit u32");
        }
        out
    }

    /// `experiment_p9`'s within-bucket permutation, unchanged.
    fn permute_within_buckets(
        mesh: &MeshBuffer<f64>,
        buckets: &[Vec<u32>],
        rng: &mut Rng,
    ) -> MeshBuffer<f64> {
        let n = mesh.positions.len();
        let mut to: Vec<u32> = (0..n as u32).collect();
        for members in buckets {
            if members.len() < 2 {
                continue;
            }
            let mut shuffled = members.clone();
            rng.shuffle(&mut shuffled);
            for (slot, &member) in members.iter().zip(&shuffled) {
                to[member as usize] = *slot;
            }
        }

        let mut out = MeshBuffer::<f64>::new();
        out.positions = vec![[0.0; 3]; n];
        out.normals = vec![[0.0; 3]; n];
        for (old, &slot) in to.iter().enumerate() {
            out.positions[slot as usize] = mesh.positions[old];
            out.normals[slot as usize] = mesh.normals[old];
        }
        out.indices = mesh.indices.iter().map(|&i| to[i as usize]).collect();
        out
    }

    /// The buckets one reference weld of `mesh` finds, `experiment_p9`'s way.
    fn buckets_of(mesh: &MeshBuffer<f64>, epsilon: f64) -> Vec<Vec<u32>> {
        let mut probe = mesh.clone();
        let mut welder = Welder::<f64>::new();
        welder.weld(&mut probe, epsilon).expect("valid epsilon");
        let mut groups: BTreeMap<u32, Vec<u32>> = BTreeMap::new();
        for (v, &to) in welder.remap().iter().enumerate() {
            groups.entry(to).or_default().push(v as u32);
        }
        groups.into_values().filter(|g| g.len() > 1).collect()
    }

    /// What one arm did across one family of permutations.
    #[derive(Clone, Copy)]
    struct Spread {
        lo: usize,
        hi: usize,
        distinct_raw: usize,
        distinct_canonical: usize,
    }

    impl Spread {
        fn spread(self) -> usize {
            self.hi - self.lo
        }
    }

    /// Accumulates one arm's outputs over a family of permutations.
    #[derive(Default)]
    struct Watch {
        lo: usize,
        hi: usize,
        raw: BTreeSet<u64>,
        canonical: BTreeSet<u64>,
    }

    impl Watch {
        fn new() -> Self {
            Self {
                lo: usize::MAX,
                hi: 0,
                raw: BTreeSet::new(),
                canonical: BTreeSet::new(),
            }
        }

        fn see(&mut self, mesh: &MeshBuffer<f64>) {
            self.lo = self.lo.min(mesh.positions.len());
            self.hi = self.hi.max(mesh.positions.len());
            self.raw.insert(raw_digest(mesh));
            self.canonical.insert(canonical_digest(mesh));
        }

        fn finish(&self) -> Spread {
            Spread {
                lo: self.lo,
                hi: self.hi,
                distinct_raw: self.raw.len(),
                distinct_canonical: self.canonical.len(),
            }
        }
    }

    // ─── the counter windows ───────────────────────────────────────────────

    /// One counter window's totals, per repetition of the body.
    #[derive(Clone, Copy, Default)]
    struct Counted {
        cycles: f64,
        instructions: f64,
        nanos: f64,
    }

    impl Counted {
        fn minus(self, base: Self) -> Self {
            Self {
                cycles: self.cycles - base.cycles,
                instructions: self.instructions - base.instructions,
                nanos: self.nanos - base.nanos,
            }
        }
    }

    /// One counter window over `inner` repetitions of `body`, divided by
    /// `inner`.
    ///
    /// The `perf_event` calls are outside the counted region, and windows are
    /// never nested: `Probe` opens six hardware counters and Zen 3 has six, so
    /// a nested pair would multiplex and `worst_ratio` refuses (`R-121`).
    fn window(probe: &mut Probe, inner: usize, mut body: impl FnMut()) -> Counted {
        probe.reset_and_enable();
        let started = Instant::now();
        for _ in 0..inner {
            body();
        }
        let nanos = started.elapsed().as_nanos() as f64;
        probe.disable();
        let counts = probe.read();
        assert!(
            counts.worst_ratio() >= MIN_TIME_RATIO,
            "a counter ran only {:.1}% of the time it was enabled, so its value is an \
             extrapolation rather than a measurement",
            counts.worst_ratio() * 100.0
        );
        let scale = 1.0 / inner as f64;
        Counted {
            cycles: counts.cycles.count as f64 * scale,
            instructions: counts.instructions.count as f64 * scale,
            nanos: nanos * scale,
        }
    }

    /// The median of an odd-or-even sample, taken as the lower middle so the
    /// reported value is one that was actually observed.
    fn median(mut xs: Vec<f64>) -> f64 {
        xs.sort_by(f64::total_cmp);
        xs[(xs.len() - 1) / 2]
    }

    /// Which of the three sibling windows a repetition is measuring.
    ///
    /// `Restore` is the prefix: the buffer copy both welds are preceded by,
    /// measured on its own so the reported figures can be the difference.
    /// Diluting the ratio with work common to both arms can only push it
    /// toward 1.0, which is the direction that would make C2's bar hold by
    /// arithmetic rather than by measurement.
    #[derive(Clone, Copy)]
    enum Arm {
        Restore,
        Single,
        Double,
    }

    /// A pristine copy of the mesh, restored before every weld.
    struct Template {
        positions: Vec<[f64; 3]>,
        normals: Vec<[f64; 3]>,
        indices: Vec<u32>,
    }

    impl Template {
        fn of(mesh: &MeshBuffer<f64>) -> Self {
            Self {
                positions: mesh.positions.clone(),
                normals: mesh.normals.clone(),
                indices: mesh.indices.clone(),
            }
        }

        fn restore(&self, into: &mut MeshBuffer<f64>) {
            into.positions.clear();
            into.positions.extend_from_slice(&self.positions);
            into.normals.clear();
            into.normals.extend_from_slice(&self.normals);
            into.indices.clear();
            into.indices.extend_from_slice(&self.indices);
        }
    }

    /// What the coincidence buckets of one block look like, and — the question
    /// the registered fixture turns on — whether any of them spans two chunks.
    ///
    /// `MeshBuffer::append` is a concatenation, so permuting chunk append
    /// order moves whole blocks of slots and **preserves the relative order of
    /// every pair of vertices inside one piece**. A bucket whose members all
    /// come from one chunk therefore cannot have its first-fit outcome changed
    /// by that perturbation, however the chunks are shuffled. Only a
    /// **cross-chunk** bucket can move at all, and only a cross-chunk bucket of
    /// three or more can move the vertex *count*, because a two-member bucket
    /// elects one representative whichever of the two is first.
    ///
    /// This is the census that says whether the registered fixture's
    /// perturbation can reach `P-9`'s defect, rather than leaving a zero
    /// spread to be read as either "no defect" or "no perturbation".
    #[derive(Clone, Copy, Default)]
    struct BucketCensus {
        /// Buckets of two or more, `experiment_p9`'s definition.
        buckets: usize,
        /// Of those, how many span more than one chunk.
        cross_chunk: usize,
        /// `P-9`'s `buckets_of_three_or_more`.
        big: usize,
        /// Buckets of three or more that span more than one chunk. **The
        /// decisive number**: it is the only kind of bucket whose member count
        /// chunk append order can change.
        cross_chunk_big: usize,
        /// The largest bucket, so a reader can see there was something to
        /// permute.
        largest: usize,
    }

    /// Take the census.
    fn census(buckets: &[Vec<u32>], chunk_of: &[u32]) -> BucketCensus {
        let mut out = BucketCensus::default();
        for members in buckets {
            out.buckets += 1;
            out.largest = out.largest.max(members.len());
            let first = chunk_of[members[0] as usize];
            let spans = members.iter().any(|&m| chunk_of[m as usize] != first);
            if spans {
                out.cross_chunk += 1;
            }
            if members.len() >= 3 {
                out.big += 1;
                if spans {
                    out.cross_chunk_big += 1;
                }
            }
        }
        out
    }

    /// Which chunk each vertex of the identity-order block came from.
    fn chunk_of_vertex(pieces: &[MeshBuffer<f64>]) -> Vec<u32> {
        let mut out = Vec::new();
        for (c, piece) in pieces.iter().enumerate() {
            out.extend(core::iter::repeat_n(c as u32, piece.positions.len()));
        }
        out
    }

    // ─── one row ───────────────────────────────────────────────────────────

    /// Everything one `(field, extractor)` produced.
    struct Measured {
        field: &'static str,
        extractor: &'static str,
        vertices_before_weld: usize,
        vertices_single: usize,
        vertices_double: usize,
        mirror_matches_shipped: bool,
        chunk_single: Spread,
        chunk_double: Spread,
        chunk_indexed: Spread,
        chunk_single_extended: Spread,
        chunk_double_extended: Spread,
        bucket_single: Spread,
        bucket_double: Spread,
        published_distinct: usize,
        published_spread: usize,
        published_big_buckets: usize,
        buckets: BucketCensus,
        chunk_orders_distinct: usize,
        report: DoubleReport,
        peak_bytes_single: usize,
        peak_bytes_double: usize,
        single: Counted,
        double: Counted,
        restore: Counted,
        gross_single: Counted,
        gross_double: Counted,
    }

    /// Measure one `(field, extractor)`.
    fn measure<E: Extractor<f64>>(
        field: &'static str,
        extractor_name: &'static str,
        sdf: &impl Sdf<Scalar = f64>,
        extractor: &mut E,
        probe: &mut Probe,
    ) -> Measured {
        let layout = ChunkLayout::<f64>::new(CELLS, CELL_SIZE, [ORIGIN; 3]).expect("valid layout");
        let epsilon = epsilon_for(CELL_SIZE);
        let pieces = eight_pieces(sdf, &layout, extractor);
        let identity: Vec<usize> = (0..CHUNKS).collect();
        let base = joined(&pieces, &identity);
        assert!(
            !base.indices.is_empty(),
            "{field}/{extractor_name} meshed the block to nothing, so there is no weld to permute"
        );

        let mut single = SingleBuffer::default();
        let mut double = DoubleBuffer::default();

        // ── the mirror, proved against the shipped Welder ────────────────────
        let mut mine = base.clone();
        single.weld(&mut mine, epsilon);
        let mut theirs = base.clone();
        let mut welder = Welder::<f64>::new();
        welder.weld(&mut theirs, epsilon).expect("valid epsilon");
        let mirror_matches_shipped = mine.indices == theirs.indices
            && mine.positions.len() == theirs.positions.len()
            && mine
                .positions
                .iter()
                .zip(&theirs.positions)
                .all(|(a, b)| a.map(f64::to_bits) == b.map(f64::to_bits))
            && mine
                .normals
                .iter()
                .zip(&theirs.normals)
                .all(|(a, b)| a.map(f64::to_bits) == b.map(f64::to_bits))
            && single.remap == welder.remap();
        assert!(
            mirror_matches_shipped,
            "{field}/{extractor_name}: the bench-local single-buffer mirror does not reproduce \
             isomesh::weld::Welder, so nothing measured against it means anything"
        );
        let vertices_single = mine.positions.len();

        // ── family one: chunk append order, the registered fixture ───────────
        let sweep = |double: &mut DoubleBuffer,
                     single: &mut SingleBuffer,
                     count: u32|
         -> (Spread, Spread, Spread, DoubleReport) {
            let mut w_single = Watch::new();
            let mut w_double = Watch::new();
            let mut w_indexed = Watch::new();
            let mut report = DoubleReport::default();
            for seed in 0..count {
                let mut order = identity.clone();
                seeded(seed).shuffle(&mut order);
                let start = joined(&pieces, &order);

                let mut a = start.clone();
                single.weld(&mut a, epsilon);
                w_single.see(&a);

                let mut b = start.clone();
                report = double.weld(&mut b, epsilon, Order::Canonical);
                w_double.see(&b);

                let mut c = start;
                double.weld(&mut c, epsilon, Order::InputIndex);
                w_indexed.see(&c);
            }
            (
                w_single.finish(),
                w_double.finish(),
                w_indexed.finish(),
                report,
            )
        };

        let (chunk_single, chunk_double, chunk_indexed, report) =
            sweep(&mut double, &mut single, PERMUTATIONS);
        let (chunk_single_extended, chunk_double_extended, _, _) =
            sweep(&mut double, &mut single, EXTENDED);

        // ── family two: P-9's own protocol, reproduced exactly ───────────────
        let buckets = buckets_of(&base, epsilon);
        let mut w_bucket_single = Watch::new();
        let mut w_bucket_double = Watch::new();
        for seed in 0..PERMUTATIONS {
            let mut rng = seeded(seed);
            let shuffled = permute_within_buckets(&base, &buckets, &mut rng);

            // The shipped `Welder`, not the mirror: this arm's job is to
            // reproduce `P-9`'s published numbers, and `P-9` used that type.
            let mut a = shuffled.clone();
            let mut w = Welder::<f64>::new();
            w.weld(&mut a, epsilon).expect("valid epsilon");
            w_bucket_single.see(&a);

            let mut b = shuffled;
            double.weld(&mut b, epsilon, Order::Canonical);
            w_bucket_double.see(&b);
        }
        let (published_distinct, published_spread, published_big_buckets) =
            published(field, extractor_name);

        // The census, and the count of distinct chunk orders the sweep actually
        // drew. Neither is a clause; both exist so a zero `double_spread`
        // cannot be read as "no defect" when the truth is "no perturbation"
        // (`BucketCensus`'s doc, above).
        let bucket_census = census(&buckets, &chunk_of_vertex(&pieces));
        let chunk_orders_distinct = {
            let mut seen: Vec<Vec<usize>> = Vec::new();
            for seed in 0..PERMUTATIONS {
                let mut order = identity.clone();
                seeded(seed).shuffle(&mut order);
                if !seen.contains(&order) {
                    seen.push(order);
                }
            }
            seen.len()
        };

        let mut sample = base.clone();
        let vertices_double = {
            double.weld(&mut sample, epsilon, Order::Canonical);
            sample.positions.len()
        };

        // ── cost: three sibling windows, rotated ─────────────────────────────
        let template = Template::of(&base);
        let mut scratch = base.clone();
        let inner = {
            let started = Instant::now();
            template.restore(&mut scratch);
            single.weld(&mut scratch, epsilon);
            let pass = started.elapsed().as_nanos() as f64;
            ((TARGET_BATCH_NS / pass.max(1.0)).ceil() as usize).clamp(1, MAX_INNER)
        };

        let mut singles = Vec::with_capacity(REPS);
        let mut doubles = Vec::with_capacity(REPS);
        let mut restores = Vec::with_capacity(REPS);
        for rep in 0..REPS {
            // The three windows are siblings, never nested, and their order
            // rotates with the repetition index so no arm is permanently the
            // one that inherits another's cache state.
            for k in 0..3 {
                let arm = [Arm::Restore, Arm::Single, Arm::Double][(rep + k) % 3];
                let counted = window(probe, inner, || {
                    template.restore(&mut scratch);
                    match arm {
                        Arm::Restore => {}
                        Arm::Single => single.weld(&mut scratch, epsilon),
                        Arm::Double => {
                            double.weld(&mut scratch, epsilon, Order::Canonical);
                        }
                    }
                    black_box(&scratch);
                });
                match arm {
                    Arm::Restore => restores.push(counted),
                    Arm::Single => singles.push(counted),
                    Arm::Double => doubles.push(counted),
                }
            }
        }

        let pick = |xs: &[Counted]| Counted {
            cycles: median(xs.iter().map(|c| c.cycles).collect()),
            instructions: median(xs.iter().map(|c| c.instructions).collect()),
            nanos: median(xs.iter().map(|c| c.nanos).collect()),
        };
        let gross_single = pick(&singles);
        let gross_double = pick(&doubles);
        let restore = pick(&restores);

        Measured {
            field,
            extractor: extractor_name,
            vertices_before_weld: base.positions.len(),
            vertices_single,
            vertices_double,
            mirror_matches_shipped,
            chunk_single,
            chunk_double,
            chunk_indexed,
            chunk_single_extended,
            chunk_double_extended,
            bucket_single: w_bucket_single.finish(),
            bucket_double: w_bucket_double.finish(),
            published_distinct,
            published_spread,
            published_big_buckets,
            buckets: bucket_census,
            chunk_orders_distinct,
            report,
            peak_bytes_single: single.scratch_bytes(),
            peak_bytes_double: double.scratch_bytes(),
            single: gross_single.minus(restore),
            double: gross_double.minus(restore),
            restore,
            gross_single,
            gross_double,
        }
    }

    // ─── the sweep ─────────────────────────────────────────────────────────

    fn sweep() -> Vec<Measured> {
        let mut probe = Probe::open();
        let mut rows = Vec::new();
        isomesh::for_each_reference_field!(f64, |name, field| {
            {
                let mut e = isomesh::dual_contouring::DualContouring::<f64>::new();
                rows.push(measure(name, "dual_contouring", &field, &mut e, &mut probe));
            }
            {
                let mut e = isomesh::manifold_dual_contouring::ManifoldDualContouring::<f64>::new();
                rows.push(measure(
                    name,
                    "manifold_dual_contouring",
                    &field,
                    &mut e,
                    &mut probe,
                ));
            }
        });
        rows
    }

    pub(crate) fn run(run: &mut crate::common::experiment::Run) {
        let rows = sweep();

        println!(
            "{:<16} {:<26} {:>9} {:>9} {:>9} {:>8} {:>8} {:>7}",
            "field", "extractor", "base", "dbl", "p9(pub)", "instr", "ns", "verts"
        );

        let mut c1_all = true;
        let mut c2_all = true;
        let mut c3_all = true;
        let mut control_all = true;

        for m in &rows {
            let baseline_spread = m.chunk_single.spread();
            let double_spread = m.chunk_double.spread();
            let indexed_spread = m.chunk_indexed.spread();
            let bucket_spread = m.bucket_single.spread();

            // The registered VACUITY CONTROL, and the wider form of it: the
            // within-bucket arm has to reproduce `P-9` exactly — spread and
            // distinct-output count — on all sixteen rows.
            let p9_reproduced = bucket_spread == m.published_spread
                && m.bucket_single.distinct_raw == m.published_distinct;
            // A row can only *demonstrate* removal if there was something to
            // remove. `baseline_spread` is the registered column and the
            // registered control names it.
            let control_fires = baseline_spread > 0;
            control_all &= p9_reproduced;

            let c1 = if control_fires {
                if double_spread == 0 { "true" } else { "false" }
            } else {
                "vacuous"
            };
            if c1 == "false" {
                c1_all = false;
            }

            let instruction_cost_ratio = m.double.instructions / m.single.instructions;
            let cost_ratio = m.double.nanos / m.single.nanos;
            let c2 = instruction_cost_ratio <= 1.25;
            c2_all &= c2;

            let c3 = m.chunk_double.distinct_canonical == 1
                && m.bucket_double.distinct_canonical == 1
                && m.chunk_double_extended.distinct_canonical == 1;
            c3_all &= c3;

            let cells = CELLS_TOTAL as f64;
            let ghz = m.single.cycles / m.single.nanos;

            println!(
                "{:<16} {:<26} {baseline_spread:>9} {double_spread:>9} {:>9} \
                 {instruction_cost_ratio:>8.3} {cost_ratio:>8.3} {:>7}",
                m.field,
                m.extractor,
                format!("{bucket_spread}({})", m.published_spread),
                m.vertices_single,
            );

            run.record(&[
                ("field", m.field.to_string()),
                ("extractor", m.extractor.to_string()),
                ("permutations", PERMUTATIONS.to_string()),
                ("baseline_spread", baseline_spread.to_string()),
                ("double_buffered_spread", double_spread.to_string()),
                (
                    "mesh_identical_across_permutations",
                    (m.chunk_double.distinct_canonical == 1).to_string(),
                ),
                (
                    "ns_per_cell_single",
                    format!("{:.6}", m.single.nanos / cells),
                ),
                (
                    "ns_per_cell_double",
                    format!("{:.6}", m.double.nanos / cells),
                ),
                ("cost_ratio", format!("{cost_ratio:.4}")),
                ("peak_bytes_single", m.peak_bytes_single.to_string()),
                ("peak_bytes_double", m.peak_bytes_double.to_string()),
                ("c1_holds", c1.to_string()),
                ("c2_holds", c2.to_string()),
                ("c3_holds", c3.to_string()),
                // ── the vacuity control, in full ─────────────────────────────
                ("baseline_spread_within_bucket", bucket_spread.to_string()),
                ("p9_published_spread", m.published_spread.to_string()),
                (
                    "p9_published_distinct_outputs",
                    m.published_distinct.to_string(),
                ),
                (
                    "distinct_outputs_within_bucket",
                    m.bucket_single.distinct_raw.to_string(),
                ),
                ("p9_reproduced", p9_reproduced.to_string()),
                ("control_fires", control_fires.to_string()),
                (
                    "double_buffered_spread_within_bucket",
                    m.bucket_double.spread().to_string(),
                ),
                // ── the census: was there anything to permute? ──────────────
                ("buckets_of_two_or_more", m.buckets.buckets.to_string()),
                ("buckets_of_three_or_more", m.buckets.big.to_string()),
                (
                    "p9_published_buckets_of_three_or_more",
                    m.published_big_buckets.to_string(),
                ),
                ("buckets_spanning_chunks", m.buckets.cross_chunk.to_string()),
                (
                    "buckets_of_three_or_more_spanning_chunks",
                    m.buckets.cross_chunk_big.to_string(),
                ),
                ("largest_bucket", m.buckets.largest.to_string()),
                ("chunk_orders_distinct", m.chunk_orders_distinct.to_string()),
                // ── double-buffering alone, without the canonical tie-break ──
                ("spread_double_indexed", indexed_spread.to_string()),
                (
                    "mesh_identical_across_permutations_indexed",
                    (m.chunk_indexed.distinct_canonical == 1).to_string(),
                ),
                (
                    "mesh_identical_across_permutations_single",
                    (m.chunk_single.distinct_canonical == 1).to_string(),
                ),
                (
                    "mesh_identical_across_permutations_within_bucket",
                    (m.bucket_double.distinct_canonical == 1).to_string(),
                ),
                (
                    "distinct_canonical_meshes_single",
                    m.chunk_single.distinct_canonical.to_string(),
                ),
                (
                    "distinct_canonical_meshes_double",
                    m.chunk_double.distinct_canonical.to_string(),
                ),
                // ── the extended chunk-order sweep ───────────────────────────
                ("permutations_extended", EXTENDED.to_string()),
                (
                    "baseline_spread_extended",
                    m.chunk_single_extended.spread().to_string(),
                ),
                (
                    "double_buffered_spread_extended",
                    m.chunk_double_extended.spread().to_string(),
                ),
                // ── what the mechanism changes about the answer ──────────────
                ("vertices_before_weld", m.vertices_before_weld.to_string()),
                ("vertices_single", m.vertices_single.to_string()),
                ("vertices_double", m.vertices_double.to_string()),
                (
                    "vertex_delta",
                    (m.vertices_double as i64 - m.vertices_single as i64).to_string(),
                ),
                (
                    "max_move_over_epsilon",
                    format!("{:.6}", m.report.max_move_over_epsilon),
                ),
                (
                    "vertices_in_large_components",
                    m.report.vertices_in_large_components.to_string(),
                ),
                (
                    "label_propagation_iterations",
                    m.report.iterations.to_string(),
                ),
                (
                    "mirror_matches_shipped",
                    m.mirror_matches_shipped.to_string(),
                ),
                // ── cost, in every form ──────────────────────────────────────
                ("cells", CELLS_TOTAL.to_string()),
                (
                    "instructions_per_cell_single",
                    format!("{:.4}", m.single.instructions / cells),
                ),
                (
                    "instructions_per_cell_double",
                    format!("{:.4}", m.double.instructions / cells),
                ),
                (
                    "instruction_cost_ratio",
                    format!("{instruction_cost_ratio:.4}"),
                ),
                (
                    "cycles_per_cell_single",
                    format!("{:.4}", m.single.cycles / cells),
                ),
                (
                    "cycles_per_cell_double",
                    format!("{:.4}", m.double.cycles / cells),
                ),
                (
                    "cycle_cost_ratio",
                    format!("{:.4}", m.double.cycles / m.single.cycles),
                ),
                ("cost_ratio_undiluted", format!("{cost_ratio:.4}")),
                (
                    "cost_ratio_with_restore",
                    format!("{:.4}", m.gross_double.nanos / m.gross_single.nanos),
                ),
                (
                    "instruction_cost_ratio_with_restore",
                    format!(
                        "{:.4}",
                        m.gross_double.instructions / m.gross_single.instructions
                    ),
                ),
                (
                    "ns_per_cell_restore",
                    format!("{:.6}", m.restore.nanos / cells),
                ),
                (
                    "instructions_per_cell_restore",
                    format!("{:.4}", m.restore.instructions / cells),
                ),
                (
                    "ns_per_vertex_single",
                    format!("{:.4}", m.single.nanos / m.vertices_before_weld as f64),
                ),
                (
                    "ns_per_vertex_double",
                    format!("{:.4}", m.double.nanos / m.vertices_before_weld as f64),
                ),
                ("ghz", format!("{ghz:.4}")),
                ("c2_holds_ns", (cost_ratio <= 1.25).to_string()),
                (
                    "peak_bytes_ratio",
                    format!(
                        "{:.4}",
                        m.peak_bytes_double as f64 / m.peak_bytes_single as f64
                    ),
                ),
                (
                    "target_feature_popcnt",
                    cfg!(target_feature = "popcnt").to_string(),
                ),
                ("count_ones_calls_per_vertex", "0".to_string()),
            ]);
        }

        // ── the aggregate the registration is scored on ──────────────────────
        let p9_rows: Vec<&Measured> = rows.iter().filter(|m| m.field == "noise_cavity").collect();
        println!();
        for m in &p9_rows {
            println!(
                "P-9 row {}/{}: baseline_spread {} (chunk order), {} (within bucket, P-9 \
                 published {}), double-buffered {}",
                m.field,
                m.extractor,
                m.chunk_single.spread(),
                m.bucket_single.spread(),
                m.published_spread,
                m.chunk_double.spread(),
            );
        }
        println!(
            "\nVACUITY CONTROL: P-9 reproduced on all 16 rows: {control_all}\n\
             C1 (spread removed where the control fires): {c1_all}\n\
             C2 (instruction cost ratio <= 1.25 on every row): {c2_all}\n\
             C3 (canonical mesh invariant under every permutation): {c3_all}"
        );
    }
}

fn main() {
    if !std::env::args().any(|arg| arg == "--bench") {
        return;
    }

    let prereg = isomesh::experiment!("P-119");

    #[cfg(target_os = "linux")]
    common::experiment::run(prereg, experiment::run);

    // `experiment_p12`'s precedent. C2's verdict is an instruction ratio and
    // the only instrument that can read one is `perf_event_open`; `M-281`
    // forbids a clock carrying it, and a recorded zero would be a fabrication.
    #[cfg(not(target_os = "linux"))]
    {
        eprintln!(
            "{} scores its cost clause on instruction counts, and this platform has no \
             `perf_event_open`. There is no clock substitute.",
            prereg.id
        );
        std::process::exit(1);
    }
}
