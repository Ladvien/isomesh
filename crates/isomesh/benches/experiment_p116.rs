//! **P-116 — GRAPHGEN's decision-table pipeline against a case table that is already `const`-derived.**
//!
//! Ticket: R-116. Pre-registered before this harness existed.
//!
//! ```bash
//! cargo bench --bench experiment_p116
//! ```
//!
//! Writes `docs/experiments/p-116.csv`. **Linux only**: C2 is *instructions per
//! cell*, which comes from `perf_event_open`, and off Linux there is nothing to
//! degrade to. `experiment_p12` is the precedent — the harness refuses and exits
//! 1 rather than record a fabricated zero.
//!
//! # What was missing — and the premise, which was wrong
//!
//! **The motivation the research doc gives for this row does not exist.** The doc
//! frames GRAPHGEN as a defence against *mistranscribing* the 256-entry Marching
//! Cubes case table. `CASES: [McCase; 256]` at
//! `crates/isomesh/src/marching_cubes/table.rs:180` is **already derived at
//! compile time** by `const fn build_cases()` (`:182-194`), which calls
//! `triangulate(segment_links(case, 0))` on all 256 sign patterns. The two
//! transcribed Bourke tables — `reference.rs:29` `BOURKE_EDGE_TABLE` and `:53`
//! `BOURKE_TRI_TABLE` — exist **only** as a test cross-check, imported at
//! `marching_cubes/tests.rs:15` and nowhere on the extraction path. So the
//! mistranscription risk **is already mitigated in the shipped path**, and this
//! row does not compare a synthesised decision tree against a transcription. It
//! compares one against a `const fn` derivation.
//!
//! `✗50` remains the real incident in this area and **it was a bound, not a
//! table**: a *sampled* maximum that became a release-build panic. A code
//! generator would not have caught it, because the table was never the thing that
//! was wrong.
//!
//! What was genuinely missing is the number this row exists to produce: **how
//! much GRAPHGEN's pipeline actually compresses a decision table whose entries
//! are all distinct.** Bolelli et al.'s DRAG compression (`10.1109/tpami.2021.
//! 3055337`, in the corpus) is driven by *don't-care* entries and by conditions
//! whose combinations share an action. A Marching Cubes case table is a decision
//! table over eight boolean corner conditions with **no don't-cares at all** —
//! all 256 assignments are specified — so the question is how many actions are
//! shared. That is measured here rather than assumed, and it decides the whole
//! row: [`experiment::run`] prints `distinct_actions` before it prints anything
//! else.
//!
//! # The pipeline, stated before it runs
//!
//! Four stages, each producing a column:
//!
//! 1. **The decision table.** 256 rows. Conditions `c0…c7` are
//!    `is_inside(corner_value[k])` — `cube.rs:171`, `value < 0`, so an exact zero
//!    is *outside*. The action is the triangulation, and two rows share an action
//!    when they agree on `count`, `centroids` and `triangles[..count]`.
//!    `dont_care_entries` is 0 by construction and is a column so a reader does
//!    not have to take that on trust.
//! 2. **The optimal decision tree, by minimum average path length.** A dynamic
//!    program over subcubes: state `(mask, values)` — which conditions have been
//!    tested and to what — memoised over all 6,561 reachable states. A state
//!    whose members all carry one action is a leaf; otherwise the cost is the
//!    number of members plus the best split over the untested conditions. **The
//!    same program is run for the maximum**, and `min_average_path_len` and
//!    `max_average_path_len` are both columns. When they are equal, *every*
//!    ordering of the conditions is optimal and the synthesis has no freedom —
//!    which is the interesting outcome and is asserted rather than hoped, because
//!    the emitted code below uses a **fixed** condition order and would not be
//!    optimal if the two differed.
//! 3. **The DRAG.** The tree is hash-consed bottom-up: a node's identity is
//!    `Leaf(action)` or `Branch(condition, false_child, true_child)`, and two
//!    nodes merge exactly when their whole subtrees are identical. `tree_nodes`,
//!    `drag_nodes` and the split into `*_internal_nodes` / `*_leaf_nodes` are all
//!    columns, because *where* the compression lands is the finding.
//! 4. **Rust.** Emitted as text — a self-contained module with the leaf action
//!    array and the decision tree as nested `if`s — and **compiled by `rustc`**
//!    into `CARGO_TARGET_TMPDIR`, `emitted_rust_compiles`. A generator whose
//!    output has never been through a compiler has not emitted Rust.
//!
//! # How the *measured* generated form gets into this binary
//!
//! `M-281` forbids comparing figures across binaries, so the form C2 counts
//! instructions on has to be compiled into **this** bench, with these flags.
//! `rustc`-ing the emitted text into a separate artefact and timing that would
//! break exactly that rule, and `dlopen` is not available — the workspace
//! `forbid`s `unsafe_code`.
//!
//! So the emitted tree is *also* expanded in place by [`decision_tree`], a
//! recursive `macro_rules!` that takes the condition order as literal arguments
//! and emits the same complete nested-`if` tree, each leaf a
//! `const { LEAF_OF_CASE[…] }` block — a compile-time constant, guaranteed folded
//! rather than hoped folded. Three things are checked so that the compiled
//! expansion is the pipeline's output and not a second, drifting copy:
//!
//! - **`emitted_matches_compiled`** — the leaf ids are parsed out of the emitted
//!   *text* in source order and compared, one by one, against the compiled
//!   function evaluated on the 256 patterns in the tree's own depth-first order.
//!   All 256 must agree. Text and binary are then the same decision function
//!   leaf for leaf, not merely the same on a sample.
//! - **`triangulations_identical`** — 256, the whole input space, against `CASES`.
//!   The tree has exactly 256 root-to-leaf paths and every one of them is walked,
//!   so behavioural equality here *is* structural equality; there is no sampling
//!   gap to argue about.
//! - **`drag_internal_nodes == tree_internal_nodes`**, asserted. The macro emits a
//!   pure tree, so it can only *be* the DRAG if the DRAG shares nothing above the
//!   leaves. Leaf sharing it does realise, through the `LEAF_OF_CASE` indirection
//!   into a `GENERATED_LEAVES` array sized by the pipeline at compile time.
//!
//! # The two arms, held identical except for the mechanism
//!
//! `P-122`'s discipline. Both arms are the same Marching Cubes march over the
//! same buffers, sharing [`experiment::gather`] (the eight-corner load) and
//! `emit_cell` (the whole payload), so the *only* difference between them is how
//! a `&'static McCase` is reached from eight corner values:
//!
//! - **table** — `case_of`, eight `case |= 1 << c` under `is_inside`
//!   (`marching_cubes/mod.rs:259-268`), then `&CASES[case]`. Branchless, one
//!   indexed load.
//! - **generated** — the eight-deep emitted branch chain, then
//!   `&GENERATED_LEAVES[leaf]`. Eight compares, eight branches, one indexed load.
//!
//! Both meshes are checked bit-identical to each other **and to
//! `MarchingCubes::extract`** on every row (`mesh_identical`,
//! `mesh_identical_to_shipped`). `R-120` and `R-121` both caught real defects that
//! way, and it is what licenses reading a ratio off a mirror at all (`M-279`).
//!
//! # Which instrument produced the C3 verdict
//!
//! `crates/isomesh/src/**` is read-only for the whole of Phase 25, and
//! `marching_cubes::proofs` is `cfg(kani)` inside it — so **a `#[kani::proof]`
//! cannot be placed over the generated form.** That is a fact about the phase,
//! not about Kani, and it is stated rather than worked around.
//!
//! C3's verdict therefore reads an **exhaustive 256-pattern Rust check** of
//! `proofs.rs:65`'s four properties, copied here with the assertions turned into
//! recorded violations so the count of checks made is itself a column
//! (`property_checks`) — `experiment_p64.rs:169-183`'s discipline, because
//! *"SUCCESSFUL over zero checks"* prints the same word as a proof. 256 patterns
//! is the entire input space, so the Rust check and a bounded model check settle
//! the same set of states; what Kani adds is that it settles them symbolically.
//!
//! Kani **is** installed here, so it is run and reported: `kani_checks` is the
//! check count from `shipped_case_table_is_indexable_for_every_sign_pattern`, a
//! proof about **`CASES`**, which transfers to the generated form through C1's
//! pointwise identity over all 256 patterns and through nothing else. The
//! registration's named sabotage precedent, `the_properties_can_fail`, is run
//! beside it. Neither is the verdict; both are corroboration, and the columns say
//! which is which.
//!
//! # SHARE
//!
//! Every clause's reachable share, as a column.
//!
//! - **C1 has no share.** It is an equality over an enumerated population:
//!   `triangulations_identical / patterns_tested`, and the bar is 1. Both are
//!   columns; the denominator is exact by construction, 256.
//! - **C2 is a *not-slower* bound, not a speedup, so no share is claimed** — the
//!   registration says so and this harness does not smuggle one back in. The
//!   quantity is `ratio = instructions_per_cell_generated /
//!   instructions_per_cell_table` against a bar of **1.0**. Instructions are
//!   deterministic on this machine, so that bar is exact and needs no tolerance;
//!   `instruction_ratio_spread` over the repetitions is a column and is what
//!   says so. `cycle_ratio`, `branch_miss_ratio` and `ns_per_cell_*` are beside
//!   it and **carry no verdict** (`M-280`, `M-281`): `R-105` watched an identical
//!   binary's cycle ratio move from 0.984 to 1.035 across three runs while its
//!   instruction counts held to four figures. `ghz` is on every row for that
//!   reason.
//! - **C3 has no share.** `properties_held` out of 4, over
//!   `property_checks` assertions on all 256 patterns, and the bar is 4 with
//!   `property_violations` at 0.
//!
//! Neither arm popcounts and neither reads a `u64` bitmap, so this build's
//! missing `popcnt` — no `.cargo/config.toml`, no `target-cpu`, so
//! `u64::count_ones()` lowers to a ~12-instruction SWAR sequence — cannot reach
//! any verdict here. There are **zero** `count_ones` calls per cell in either arm.
//!
//! # VACUITY CONTROL, asserted rather than recorded
//!
//! - `patterns_tested` must be 256 and `triangulations_identical` must equal it.
//! - **`sabotage_failed`** — four sabotages, **one per property**, each mutating a
//!   real entry the way `proofs.rs:221`'s `the_properties_can_fail` does, and each
//!   required to trip *its own* property and not merely some property. Four of
//!   four, `sabotages_caught`. Without it, C3 is four assertions nobody has ever
//!   seen fire.
//! - **`property_checks` must be non-zero.** `experiment_p64`'s rule in its
//!   original form.
//! - **`mutant_tree_mismatches` must be non-zero** — two leaf ids swapped in a
//!   copy of `LEAF_OF_CASE`, and the C1 comparator required to see it. C1
//!   compares the generated form against `CASES`, and both are derived from the
//!   same `triangulate(segment_links(…))`; a comparator that cannot fail would
//!   make C1 an equality between two names for one computation.
//! - **`emitted_rust_compiles` and `emitted_matches_compiled`** must both be true,
//!   or "emitted as Rust" is a string nobody compiled and a shape nobody checked.
//! - **`min_average_path_len` must equal `max_average_path_len`**, or the fixed
//!   condition order the emitted code uses is not known to be optimal.

#![allow(
    clippy::cast_possible_truncation,
    clippy::cast_precision_loss,
    clippy::too_many_lines
)]

mod common;

#[cfg(target_os = "linux")]
mod experiment {
    use std::collections::HashMap;
    use std::fmt::Write as _;
    use std::hint::black_box;
    use std::path::PathBuf;
    use std::process::Command;
    use std::time::Instant;

    use isomesh::marching_cubes::MarchingCubes;
    use isomesh::marching_cubes::table::{
        CASES, CENTROID_BASE, EDGE_AXIS, EDGE_CORNERS, EDGE_COUNT, MAX_CENTROIDS, MAX_TRIANGLES,
        McCase, NO_EDGE, edge_offset, is_centroid, is_inside, segment_links, triangulate,
    };
    use isomesh::{MeshBuffer, Real, RuntimeShape3, Sdf};

    use crate::common::counters::{MIN_TIME_RATIO, Probe};

    // ─── the fixture ───────────────────────────────────────────────────────

    /// Sign patterns over eight boolean corner conditions. The whole input
    /// space, and the denominator of every C1 and C3 figure.
    const PATTERNS: usize = 256;

    /// The registered C2 resolution, in samples per axis. Eight reference
    /// fields at this one grid.
    const RESOLUTION: u32 = 65;

    /// Measured repetitions per row, **even on purpose**: the arms alternate
    /// which of them runs second, so an even count gives each ordering exactly
    /// half the repetitions.
    ///
    /// Sixteen is generous for the clause that carries the verdict —
    /// instructions are deterministic and `instruction_ratio_spread` is the
    /// column that proves it — and is chosen for the two figures beside it.
    /// Branch mispredictions are a few hundredths of an event per cell and a
    /// cycle count on a governed CPU is not a unit at all (`M-281`).
    const REPS: usize = 16;

    /// Untimed passes of both arms before anything is counted, so the buffers
    /// are at final capacity and the pages are faulted in.
    const WARMUP: usize = 2;

    /// About this long per counter window, so the ~28 `perf_event` system calls
    /// a window costs land outside it and cannot inflate anything.
    const TARGET_BATCH_NS: f64 = 40_000_000.0;

    /// Ceiling on the batch, so a cheap row cannot run for a minute.
    const MAX_INNER: usize = 512;

    /// C2's bar: *not slower*, so instructions per cell may not rise at all.
    /// Exact rather than tolerant, because instructions are deterministic here.
    const INSTRUCTION_BAR: f64 = 1.0;

    /// The condition order the emitted tree tests in, and the literal arguments
    /// [`decision_tree`] is invoked with below.
    ///
    /// Identity, and [`Synthesis`] asserts that identity is optimal rather than
    /// assuming it: when the minimum and maximum average path lengths coincide,
    /// every ordering is optimal and this one is free.
    const CONDITION_ORDER: [usize; 8] = [0, 1, 2, 3, 4, 5, 6, 7];

    // ─── stage 1: the decision table, derived the way the crate derives it ──

    /// Every action, `const`-evaluated by the crate's own public constructors.
    ///
    /// This is `table.rs:182-194` `build_cases()` spelled again in the bench,
    /// which is legitimate because `segment_links` and `triangulate` are
    /// `pub const fn` — so the derivation is *shared*, not copied, and only the
    /// loop around it is duplicated. C1 then compares the pipeline's leaves
    /// against the shipped `static CASES`, which is the object `extract` reads.
    const ALL: [McCase; PATTERNS] = build_all();

    const fn empty_case() -> McCase {
        McCase {
            count: 0,
            centroids: 0,
            triangles: [[0u8; 3]; MAX_TRIANGLES],
        }
    }

    const fn build_all() -> [McCase; PATTERNS] {
        let mut out = [empty_case(); PATTERNS];
        let mut case = 0usize;
        while case < PATTERNS {
            out[case] = triangulate(segment_links(case as u8, 0));
            case += 1;
        }
        out
    }

    /// Do two rows of the decision table carry the same action?
    ///
    /// Only the live prefix counts. `triangulate` zeroes the tail of
    /// `triangles`, so the comparison would give the same answer either way —
    /// but an action is `triangles[..count]` and saying so here is what makes
    /// the dedup below a statement about actions rather than about padding.
    const fn same_action(a: &McCase, b: &McCase) -> bool {
        if a.count != b.count || a.centroids != b.centroids {
            return false;
        }
        let mut t = 0usize;
        while t < a.count as usize {
            let mut k = 0usize;
            while k < 3 {
                if a.triangles[t][k] != b.triangles[t][k] {
                    return false;
                }
                k += 1;
            }
            t += 1;
        }
        true
    }

    /// Sign pattern → DRAG leaf id, ids assigned in first-appearance order.
    ///
    /// This *is* the DRAG's leaf sharing, materialised at compile time: two
    /// patterns carrying one action get one id and therefore one entry in
    /// [`GENERATED_LEAVES`].
    const LEAF_OF_CASE: [u8; PATTERNS] = build_leaf_of_case();

    const fn build_leaf_of_case() -> [u8; PATTERNS] {
        let mut out = [0u8; PATTERNS];
        let mut next = 0u8;
        let mut case = 0usize;
        while case < PATTERNS {
            let mut found = false;
            let mut earlier = 0usize;
            while earlier < case {
                if same_action(&ALL[case], &ALL[earlier]) {
                    out[case] = out[earlier];
                    found = true;
                    break;
                }
                earlier += 1;
            }
            if !found {
                out[case] = next;
                next += 1;
            }
            case += 1;
        }
        out
    }

    /// Distinct actions in the decision table, **sized by the pipeline**.
    ///
    /// The whole DRAG compression this row can find is `PATTERNS - DISTINCT`
    /// leaves plus whatever the internal hash-consing merges, and this constant
    /// is the first half of that number arriving at compile time.
    const DISTINCT: usize = count_distinct();

    const fn count_distinct() -> usize {
        let mut max = 0u8;
        let mut case = 0usize;
        while case < PATTERNS {
            if LEAF_OF_CASE[case] > max {
                max = LEAF_OF_CASE[case];
            }
            case += 1;
        }
        max as usize + 1
    }

    /// The DRAG's leaf actions, one per distinct action.
    ///
    /// A `static` and not a `const`: `&GENERATED_LEAVES[i]` has to be an indexed
    /// address into rodata, exactly as `&CASES[i]` is, or the generated arm
    /// would pay a stack copy of the whole array per cell and C2 would be
    /// measuring the difference between a `static` and a `const`.
    static GENERATED_LEAVES: [McCase; DISTINCT] = build_leaves();

    const fn build_leaves() -> [McCase; DISTINCT] {
        let mut out = [empty_case(); DISTINCT];
        let mut case = 0usize;
        while case < PATTERNS {
            out[LEAF_OF_CASE[case] as usize] = ALL[case];
            case += 1;
        }
        out
    }

    // ─── stage 4b: the emitted tree, expanded into this binary ──────────────

    /// The emitted decision tree, as a recursive expansion.
    ///
    /// `$acc` accumulates the packed sign pattern as a constant expression, so
    /// the leaf is `const { LEAF_OF_CASE[<literal>] }` — evaluated at compile
    /// time by construction rather than by hoping the optimiser folds a lookup.
    /// The condition order is the literal argument list, which is the one degree
    /// of freedom the emission has and the one [`Synthesis`] proves optimal.
    macro_rules! decision_tree {
        ($v:ident, $acc:expr,) => { const { LEAF_OF_CASE[$acc] } };
        ($v:ident, $acc:expr, $c:literal $($rest:literal)*) => {
            if is_inside($v[$c]) {
                decision_tree!($v, $acc | (1usize << $c), $($rest)*)
            } else {
                decision_tree!($v, $acc, $($rest)*)
            }
        };
    }

    /// Eight corner values → a DRAG leaf id, by the emitted branch chain.
    ///
    /// Eight compares and eight branches, against `case_of`'s eight branchless
    /// `or`s. The literal order here must match [`CONDITION_ORDER`], which
    /// [`Synthesis::assert_sound`] checks.
    #[inline]
    fn generated_leaf<R: Real>(v: &[R; 8]) -> u8 {
        decision_tree!(v, 0usize, 0 1 2 3 4 5 6 7)
    }

    /// `marching_cubes/mod.rs:259-268`'s index build: the comparand.
    #[inline]
    fn case_of<R: Real>(v: &[R; 8]) -> u8 {
        let mut case = 0u8;
        for (c, &value) in v.iter().enumerate() {
            if is_inside(value) {
                case |= 1 << c;
            }
        }
        case
    }

    /// The eight corner values for one sign pattern, for the 256-pattern sweep.
    ///
    /// `-1.0` for inside and `+1.0` for outside, which is what `is_inside`
    /// (`value < 0`) reads. Never `0.0` on either side: an exact zero is
    /// *outside* by the crate's convention and a fixture that used it would be
    /// testing the convention instead of the tree.
    fn values_for(pattern: u8) -> [f32; 8] {
        let mut v = [1.0f32; 8];
        for (c, slot) in v.iter_mut().enumerate() {
            if (pattern >> c) & 1 == 1 {
                *slot = -1.0;
            }
        }
        v
    }

    // ─── stages 2 and 3: synthesis and DRAG compression ────────────────────

    /// One node of the synthesised tree.
    enum Tree {
        Leaf(u8),
        Branch {
            condition: u8,
            /// The condition false.
            lo: Box<Tree>,
            /// The condition true.
            hi: Box<Tree>,
        },
    }

    /// A hash-consed DRAG node. Two nodes are the same node exactly when their
    /// whole subtrees are identical, which for a `Branch` means the same
    /// condition and the same two already-consed children.
    #[derive(PartialEq, Eq, Hash)]
    enum Node {
        Leaf(u8),
        Branch(u8, usize, usize),
    }

    /// Members of the subcube `(mask, values)`: `2^(8 - |mask|)` patterns.
    fn subcube_size(mask: u8) -> u64 {
        1u64 << (8 - mask.count_ones())
    }

    /// Do all patterns in this subcube carry one action?
    fn subcube_uniform(leaf_of: &[u8; PATTERNS], mask: u8, values: u8) -> bool {
        let mut first: Option<u8> = None;
        for (pattern, &leaf) in leaf_of.iter().enumerate() {
            if (pattern as u8) & mask != values {
                continue;
            }
            match first {
                None => first = Some(leaf),
                Some(a) if a != leaf => return false,
                Some(_) => {}
            }
        }
        true
    }

    /// Total condition tests summed over the patterns in this subcube, over
    /// **every** tree shape — the minimum when `minimise`, the maximum when not.
    ///
    /// Divided by 256 at the root this is the average path length GRAPHGEN
    /// optimises. Running it in both directions is what turns "the synthesis
    /// found an optimum" into "the objective has no preference", which is a
    /// stronger and much more useful statement about this particular table.
    fn dp(
        leaf_of: &[u8; PATTERNS],
        memo: &mut [Option<u64>],
        mask: u8,
        values: u8,
        minimise: bool,
    ) -> u64 {
        let key = (mask as usize) << 8 | values as usize;
        if let Some(v) = memo[key] {
            return v;
        }
        let total = if subcube_uniform(leaf_of, mask, values) {
            0
        } else {
            let mut best = if minimise { u64::MAX } else { 0 };
            for c in 0..8u8 {
                let bit = 1u8 << c;
                if mask & bit != 0 {
                    continue;
                }
                let split = dp(leaf_of, memo, mask | bit, values, minimise)
                    + dp(leaf_of, memo, mask | bit, values | bit, minimise);
                best = if minimise {
                    best.min(split)
                } else {
                    best.max(split)
                };
            }
            subcube_size(mask) + best
        };
        memo[key] = Some(total);
        total
    }

    /// The tree for a **fixed** condition order — the shape the emitted code has.
    fn build_fixed(
        leaf_of: &[u8; PATTERNS],
        order: [usize; 8],
        depth: usize,
        mask: u8,
        values: u8,
    ) -> Tree {
        if subcube_uniform(leaf_of, mask, values) {
            let leaf = leaf_of
                .iter()
                .enumerate()
                .find(|(p, _)| (*p as u8) & mask == values)
                .map(|(_, &l)| l)
                .expect("a subcube has at least one member");
            return Tree::Leaf(leaf);
        }
        assert!(
            depth < 8,
            "a fully tested subcube has one member and must be uniform"
        );
        let condition = order[depth] as u8;
        let bit = 1u8 << condition;
        Tree::Branch {
            condition,
            lo: Box::new(build_fixed(leaf_of, order, depth + 1, mask | bit, values)),
            hi: Box::new(build_fixed(
                leaf_of,
                order,
                depth + 1,
                mask | bit,
                values | bit,
            )),
        }
    }

    fn count_nodes(t: &Tree) -> (usize, usize) {
        match t {
            Tree::Leaf(_) => (0, 1),
            Tree::Branch { lo, hi, .. } => {
                let (li, ll) = count_nodes(lo);
                let (ri, rl) = count_nodes(hi);
                (1 + li + ri, ll + rl)
            }
        }
    }

    /// Total condition tests summed over the 256 patterns of this tree.
    fn total_path_len(t: &Tree, depth: u64, weight: u64) -> u64 {
        match t {
            Tree::Leaf(_) => depth * weight,
            Tree::Branch { lo, hi, .. } => {
                total_path_len(lo, depth + 1, weight / 2)
                    + total_path_len(hi, depth + 1, weight / 2)
            }
        }
    }

    /// Hash-cons one subtree, returning its DRAG node id.
    fn hash_cons(t: &Tree, seen: &mut HashMap<Node, usize>) -> usize {
        let node = match t {
            Tree::Leaf(a) => Node::Leaf(*a),
            Tree::Branch { condition, lo, hi } => {
                let l = hash_cons(lo, seen);
                let r = hash_cons(hi, seen);
                Node::Branch(*condition, l, r)
            }
        };
        if let Some(&id) = seen.get(&node) {
            return id;
        }
        let id = seen.len();
        seen.insert(node, id);
        id
    }

    /// The leaf ids of a tree in depth-first, **true-branch-first** order —
    /// which is the order the emitted text lists them in.
    fn leaf_sequence(t: &Tree, out: &mut Vec<u8>) {
        match t {
            Tree::Leaf(a) => out.push(*a),
            Tree::Branch { lo, hi, .. } => {
                leaf_sequence(hi, out);
                leaf_sequence(lo, out);
            }
        }
    }

    /// Everything stages 1 to 3 produced.
    struct Synthesis {
        distinct_actions: usize,
        shared_action_patterns: usize,
        tree_internal: usize,
        tree_leaves: usize,
        drag_internal: usize,
        drag_leaves: usize,
        emitted_total_path_len: u64,
        min_total_path_len: u64,
        max_total_path_len: u64,
        dp_states: usize,
        emitted_leaf_order: Vec<u8>,
    }

    impl Synthesis {
        fn new() -> Self {
            let leaf_of = LEAF_OF_CASE;
            let distinct_actions = DISTINCT;
            let shared_action_patterns = PATTERNS - distinct_actions;

            let mut memo_min = vec![None; 1 << 16];
            let mut memo_max = vec![None; 1 << 16];
            let min_total_path_len = dp(&leaf_of, &mut memo_min, 0, 0, true);
            let max_total_path_len = dp(&leaf_of, &mut memo_max, 0, 0, false);
            let dp_states = memo_min.iter().filter(|s| s.is_some()).count();

            let tree = build_fixed(&leaf_of, CONDITION_ORDER, 0, 0, 0);
            let (tree_internal, tree_leaves) = count_nodes(&tree);
            let emitted_total_path_len = total_path_len(&tree, 0, PATTERNS as u64);

            let mut seen: HashMap<Node, usize> = HashMap::new();
            hash_cons(&tree, &mut seen);
            let drag_leaves = seen.keys().filter(|n| matches!(n, Node::Leaf(_))).count();
            let drag_internal = seen.len() - drag_leaves;

            let mut emitted_leaf_order = Vec::with_capacity(PATTERNS);
            leaf_sequence(&tree, &mut emitted_leaf_order);

            Self {
                distinct_actions,
                shared_action_patterns,
                tree_internal,
                tree_leaves,
                drag_internal,
                drag_leaves,
                emitted_total_path_len,
                min_total_path_len,
                max_total_path_len,
                dp_states,
                emitted_leaf_order,
            }
        }

        fn tree_nodes(&self) -> usize {
            self.tree_internal + self.tree_leaves
        }

        fn drag_nodes(&self) -> usize {
            self.drag_internal + self.drag_leaves
        }

        fn average_path_len(&self) -> f64 {
            self.emitted_total_path_len as f64 / PATTERNS as f64
        }

        /// The two things that license the emission below, checked rather than
        /// assumed. Both were live possibilities before the run.
        fn assert_sound(&self) {
            assert_eq!(
                self.min_total_path_len, self.max_total_path_len,
                "the minimum and maximum average path lengths differ ({} vs {} tests over 256 \
                 patterns), so the objective PREFERS an ordering and the emitted code's fixed \
                 CONDITION_ORDER is not known to be optimal -- the emission would have to follow \
                 the dynamic program's per-node choice instead",
                self.min_total_path_len, self.max_total_path_len
            );
            assert_eq!(
                self.emitted_total_path_len, self.min_total_path_len,
                "the emitted fixed-order tree costs {} tests over 256 patterns against the \
                 optimum's {}, so it is not the optimal tree the registration names",
                self.emitted_total_path_len, self.min_total_path_len
            );
            assert_eq!(
                self.drag_internal,
                self.tree_internal,
                "the DRAG merges {} internal nodes, so it is not a pure tree and the \
                 macro_rules expansion -- which emits one -- is not the DRAG",
                self.tree_internal - self.drag_internal
            );
            assert_eq!(
                self.drag_leaves, self.distinct_actions,
                "the DRAG has {} leaf nodes against {} distinct actions; leaf sharing and action \
                 sharing are the same thing and disagreeing means the hash-consing and the \
                 compile-time dedup do not",
                self.drag_leaves, self.distinct_actions
            );
        }
    }

    // ─── stage 4: emitting Rust, and compiling it ──────────────────────────

    /// Delimiters around the emitted tree, so its leaf ids can be read back out
    /// of the text without parsing Rust.
    const TREE_BEGIN: &str = "// ---- decision tree begins ----";
    const TREE_END: &str = "// ---- decision tree ends ----";

    fn emit_node(out: &mut String, order: [usize; 8], t: &Tree, indent: usize) {
        let pad = " ".repeat(indent);
        match t {
            Tree::Leaf(a) => {
                let _ = writeln!(out, "{pad}{a}u8");
            }
            Tree::Branch { condition, lo, hi } => {
                debug_assert!(order.contains(&(*condition as usize)));
                let _ = writeln!(out, "{pad}if is_inside(v[{condition}]) {{");
                emit_node(out, order, hi, indent + 4);
                let _ = writeln!(out, "{pad}}} else {{");
                emit_node(out, order, lo, indent + 4);
                let _ = writeln!(out, "{pad}}}");
            }
        }
    }

    /// The whole generated module, as text: the leaf action array, the decision
    /// tree, and a `classify` that puts them together.
    ///
    /// Self-contained on purpose. A generator that emits only the tree and
    /// leaves the actions behind has emitted half a program, and `rustc` could
    /// not be asked whether the half compiles.
    fn emit_rust(order: [usize; 8], tree: &Tree, s: &Synthesis) -> String {
        let mut out = String::with_capacity(128 * 1024);
        let _ = writeln!(
            out,
            "// GENERATED by P-116's GRAPHGEN pipeline in \
             crates/isomesh/benches/experiment_p116.rs. Do not edit."
        );
        let _ = writeln!(
            out,
            "// decision table: {PATTERNS} rows, 8 boolean corner conditions, 0 don't-cares, \
             {} distinct actions.",
            s.distinct_actions
        );
        let _ = writeln!(
            out,
            "// tree_nodes = {}, drag_nodes = {}, average_path_len = {:.6}, condition order = {}.",
            s.tree_nodes(),
            s.drag_nodes(),
            s.average_path_len(),
            order.iter().map(|c| c.to_string()).collect::<String>()
        );
        let _ = writeln!(out, "#![allow(dead_code, clippy::all)]\n");
        let _ = writeln!(
            out,
            "/// `isomesh::marching_cubes::table::McCase`, restated so this module stands alone."
        );
        let _ = writeln!(out, "#[derive(Clone, Copy, Debug)]");
        let _ = writeln!(out, "pub struct McCase {{");
        let _ = writeln!(out, "    pub count: u8,");
        let _ = writeln!(out, "    pub centroids: u8,");
        let _ = writeln!(out, "    pub triangles: [[u8; 3]; {MAX_TRIANGLES}],\n}}\n");
        let _ = writeln!(
            out,
            "/// `cube::is_inside`: a negative value is inside, so an exact zero is outside."
        );
        let _ = writeln!(out, "#[inline]");
        let _ = writeln!(
            out,
            "pub fn is_inside(value: f32) -> bool {{\n    value < 0.0\n}}\n"
        );
        let _ = writeln!(
            out,
            "/// The DRAG's leaf actions, one per distinct triangulation."
        );
        let _ = writeln!(
            out,
            "pub static LEAVES: [McCase; {}] = [",
            s.distinct_actions
        );
        for leaf in &GENERATED_LEAVES {
            let tris = leaf
                .triangles
                .iter()
                .map(|t| format!("[{}, {}, {}]", t[0], t[1], t[2]))
                .collect::<Vec<_>>()
                .join(", ");
            let _ = writeln!(
                out,
                "    McCase {{ count: {}, centroids: {}, triangles: [{tris}] }},",
                leaf.count, leaf.centroids
            );
        }
        let _ = writeln!(out, "];\n");
        let _ = writeln!(
            out,
            "/// Eight corner values -> a leaf id, by the synthesised tree."
        );
        let _ = writeln!(out, "{TREE_BEGIN}");
        let _ = writeln!(out, "#[inline]");
        let _ = writeln!(out, "pub fn leaf_of(v: &[f32; 8]) -> u8 {{");
        emit_node(&mut out, order, tree, 4);
        let _ = writeln!(out, "}}");
        let _ = writeln!(out, "{TREE_END}\n");
        let _ = writeln!(
            out,
            "/// The action for one cell, which is what the table lookup returns."
        );
        let _ = writeln!(out, "#[inline]");
        let _ = writeln!(
            out,
            "pub fn classify(v: &[f32; 8]) -> &'static McCase {{\n    \
             &LEAVES[leaf_of(v) as usize]\n}}"
        );
        out
    }

    /// The leaf ids the emitted **text** lists, in source order.
    fn emitted_leaf_sequence(source: &str) -> Vec<u8> {
        let body = source
            .split_once(TREE_BEGIN)
            .expect("the emitted tree carries its opening marker")
            .1
            .split_once(TREE_END)
            .expect("the emitted tree carries its closing marker")
            .0;
        body.lines()
            .filter_map(|line| line.trim().strip_suffix("u8"))
            .map(|digits| {
                digits
                    .parse::<u8>()
                    .expect("a leaf line is an integer literal and nothing else")
            })
            .collect()
    }

    /// Write the emitted module out and put `rustc` over it.
    ///
    /// `--emit metadata` type-checks and lowers without codegen, which is all
    /// that is being claimed: the pipeline emitted **Rust**, and here is the
    /// compiler agreeing. The artefact goes to `CARGO_TARGET_TMPDIR` and is not
    /// committed; `M-281` is why the *measured* form is the in-binary expansion
    /// instead of this one.
    fn compile_emitted(source: &str) -> (bool, String, PathBuf) {
        let dir = PathBuf::from(env!("CARGO_TARGET_TMPDIR"));
        let path = dir.join("p116_generated.rs");
        std::fs::create_dir_all(&dir).expect("the bench temp dir");
        std::fs::write(&path, source).expect("writing the emitted module");
        let out = Command::new("rustc")
            .args(["--edition", "2024", "--crate-type", "lib", "--emit"])
            .arg("metadata")
            .arg("-o")
            .arg(dir.join("p116_generated.rmeta"))
            .arg(&path)
            .output();
        match out {
            Ok(o) if o.status.success() => (true, String::from("ok"), path),
            Ok(o) => (false, csv_safe(&String::from_utf8_lossy(&o.stderr)), path),
            Err(e) => (false, csv_safe(&e.to_string()), path),
        }
    }

    /// One line, no CSV separators. `Run::record` panics on a comma, a quote or
    /// a newline, and a compiler diagnostic is full of all three.
    fn csv_safe(s: &str) -> String {
        let one_line: String = s
            .chars()
            .map(|c| match c {
                ',' => ';',
                '"' | '\n' | '\r' | '\t' => ' ',
                other => other,
            })
            .collect();
        let trimmed = one_line.split_whitespace().collect::<Vec<_>>().join(" ");
        if trimmed.len() > 160 {
            trimmed.chars().take(160).collect()
        } else {
            trimmed
        }
    }

    // ─── C3: P-64's four properties, on the generated form ─────────────────

    /// What one sweep of [`check_all_properties`] saw.
    #[derive(Default, Clone, Copy)]
    struct Properties {
        /// Assertions made. Non-zero is a control: a verdict over an empty check
        /// set prints the same word as a proof (`experiment_p64.rs:169-183`).
        checks: u64,
        /// Violations per property, in `proofs.rs`' own numbering: shape,
        /// nameable, edge-is-cut, non-degenerate.
        violations: [u64; 4],
    }

    impl Properties {
        fn total_violations(&self) -> u64 {
            self.violations.iter().sum()
        }

        fn held(&self) -> usize {
            self.violations.iter().filter(|v| **v == 0).count()
        }
    }

    /// `proofs.rs:55` `edge_is_cut`. Over sign bits, "the two corner values
    /// straddle zero" is "the bits differ".
    fn edge_is_cut(case: u8, e: usize) -> bool {
        let [lo, hi] = EDGE_CORNERS[e];
        ((case >> lo) & 1) != ((case >> hi) & 1)
    }

    /// `proofs.rs:65` `check_all_properties`, with the assertions turned into
    /// **recorded** violations.
    ///
    /// One statement of the four properties, as in `proofs.rs` — and the
    /// assertion moved to the call site, which is what lets the sabotage arm be
    /// a check rather than a caught panic, and what makes `property_checks` a
    /// column instead of an unobservable count inside a proof.
    fn check_all_properties(case: u8, c: &McCase, out: &mut Properties) {
        // Property 1: shape. The fixed-size arrays a caller stack-allocates
        // cannot overflow.
        out.checks += 1;
        if c.count as usize > MAX_TRIANGLES {
            out.violations[0] += 1;
        }
        out.checks += 1;
        if c.centroids as usize > MAX_CENTROIDS {
            out.violations[0] += 1;
        }

        for t in 0..MAX_TRIANGLES {
            if t >= c.count as usize {
                break;
            }
            let tri = c.triangles[t];

            for &code in &tri {
                // Property 2: every code is nameable.
                out.checks += 1;
                if code == NO_EDGE {
                    out.violations[1] += 1;
                    continue;
                }
                if is_centroid(code) {
                    out.checks += 1;
                    if code - CENTROID_BASE >= c.centroids {
                        out.violations[1] += 1;
                    }
                } else {
                    out.checks += 1;
                    if code as usize >= EDGE_COUNT {
                        out.violations[1] += 1;
                        continue;
                    }
                    // Property 3: the named edge is cut. The load-bearing one,
                    // and a property of the PAIR (table, sign pattern).
                    out.checks += 1;
                    if !edge_is_cut(case, code as usize) {
                        out.violations[2] += 1;
                    }
                }
            }

            // Property 4: no degenerate triangle.
            out.checks += 1;
            if tri[0] == tri[1] || tri[1] == tri[2] || tri[0] == tri[2] {
                out.violations[3] += 1;
            }
        }
    }

    /// One sabotage: a real entry, mutated so that exactly one property must
    /// fail.
    struct Sabotage {
        what: &'static str,
        case: u8,
        entry: McCase,
        /// Which of the four properties this mutation must trip.
        property: usize,
    }

    /// Four sabotages, one per property, in `proofs.rs:221`'s shape.
    ///
    /// `the_properties_can_fail` corrupts one entry and requires the check to
    /// fire; this is that, once per property, so a single tautology cannot hide
    /// behind three live assertions.
    fn sabotages() -> Vec<Sabotage> {
        // A pattern with exactly one triangle, all three codes cube edges.
        let case = 1u8;
        let base = GENERATED_LEAVES[LEAF_OF_CASE[case as usize] as usize];
        assert!(
            base.count >= 1 && base.centroids == 0,
            "the sabotage fixture wants a case with at least one edge-only triangle; case {case} \
             has count {} and {} centroids",
            base.count,
            base.centroids
        );

        let uncut = (0..EDGE_COUNT)
            .find(|&e| !edge_is_cut(case, e) && !base.triangles[0].contains(&(e as u8)))
            .expect("a one-corner pattern cuts three of twelve edges, so nine are uncut")
            as u8;

        let mut over_count = base;
        over_count.count = MAX_TRIANGLES as u8 + 1;

        let mut unnameable = base;
        unnameable.triangles[0][0] = NO_EDGE;

        let mut uncut_edge = base;
        uncut_edge.triangles[0][0] = uncut;

        let mut degenerate = base;
        degenerate.triangles[0][1] = degenerate.triangles[0][0];

        vec![
            Sabotage {
                what: "count past MAX_TRIANGLES",
                case,
                entry: over_count,
                property: 0,
            },
            Sabotage {
                what: "a live triangle carries NO_EDGE",
                case,
                entry: unnameable,
                property: 1,
            },
            Sabotage {
                what: "a triangle names an uncut edge",
                case,
                entry: uncut_edge,
                property: 2,
            },
            Sabotage {
                what: "a triangle carries two equal codes",
                case,
                entry: degenerate,
                property: 3,
            },
        ]
    }

    // ─── the Kani arm: corroboration, and it says so ───────────────────────

    /// What one `cargo kani --harness` run reported.
    struct Verdict {
        checks: u64,
        failed: u64,
        status: String,
        solver_seconds: f64,
    }

    /// `experiment_p64.rs:79-133`'s parser: the check count, the failure count,
    /// the verdict and the solver time out of Kani's own output.
    fn run_kani_harness(name: &str) -> Verdict {
        let out = Command::new("cargo")
            .args(["kani", "--harness", name])
            .current_dir(env!("CARGO_MANIFEST_DIR"))
            .output();
        let Ok(out) = out else {
            return Verdict {
                checks: 0,
                failed: 0,
                status: String::from("cargo kani could not be launched"),
                solver_seconds: f64::NAN,
            };
        };
        let text = format!(
            "{}{}",
            String::from_utf8_lossy(&out.stdout),
            String::from_utf8_lossy(&out.stderr)
        );

        let mut checks = 0u64;
        let mut failed = 0u64;
        for line in text.lines() {
            if let Some(rest) = line.trim().strip_prefix("** ")
                && let Some((n, rest)) = rest.split_once(" of ")
                && let Some((m, _)) = rest.split_once(" failed")
                && let (Ok(n), Ok(m)) = (n.parse::<u64>(), m.parse::<u64>())
            {
                failed = n;
                checks = m;
            }
        }
        let status = text
            .lines()
            .find_map(|l| l.trim().strip_prefix("VERIFICATION:- "))
            .unwrap_or("NO VERDICT")
            .trim()
            .to_string();
        let solver_seconds = text
            .lines()
            .find_map(|l| {
                l.trim()
                    .strip_prefix("Verification Time: ")
                    .and_then(|v| v.trim_end_matches('s').parse::<f64>().ok())
            })
            .unwrap_or(f64::NAN);

        Verdict {
            checks,
            failed,
            status: csv_safe(&status),
            solver_seconds,
        }
    }

    /// The Kani corroboration, or an honest account of why there is none.
    struct Kani {
        version: String,
        proof: Verdict,
        sabotage: Verdict,
    }

    fn run_kani() -> Kani {
        let version = Command::new("cargo")
            .args(["kani", "--version"])
            .output()
            .ok()
            .filter(|o| o.status.success())
            .and_then(|o| String::from_utf8(o.stdout).ok())
            .map(|s| csv_safe(s.trim()));
        let Some(version) = version else {
            // The registration's own instruction: do not fake a Kani run.
            let absent = || Verdict {
                checks: 0,
                failed: 0,
                status: String::from("cargo kani is not installed on this machine"),
                solver_seconds: f64::NAN,
            };
            return Kani {
                version: String::from("unavailable"),
                proof: absent(),
                sabotage: absent(),
            };
        };
        Kani {
            version,
            proof: run_kani_harness("shipped_case_table_is_indexable_for_every_sign_pattern"),
            sabotage: run_kani_harness("the_properties_can_fail"),
        }
    }

    // ─── the pattern sweep: C1, C3 and their controls ──────────────────────

    /// Everything the 256-pattern sweep settled. Field-independent, so it is
    /// computed once and repeated on every row — the table is a property of the
    /// crate, not of the fixture.
    struct Sweep {
        synthesis: Synthesis,
        triangulations_identical: usize,
        mutant_tree_mismatches: usize,
        properties: Properties,
        sabotages_caught: usize,
        emitted_bytes: usize,
        emitted_lines: usize,
        emitted_rust_compiles: bool,
        emitted_rust_message: String,
        emitted_matches_compiled: bool,
        kani: Kani,
    }

    impl Sweep {
        fn c1_holds(&self) -> bool {
            self.triangulations_identical == PATTERNS
                && self.emitted_rust_compiles
                && self.emitted_matches_compiled
        }

        fn c3_holds(&self) -> bool {
            self.properties.held() == 4
                && self.properties.total_violations() == 0
                && self.properties.checks > 0
                && self.sabotages_caught == 4
        }
    }

    fn sweep_patterns() -> Sweep {
        let synthesis = Synthesis::new();
        synthesis.assert_sound();

        let tree = build_fixed(&LEAF_OF_CASE, CONDITION_ORDER, 0, 0, 0);
        let source = emit_rust(CONDITION_ORDER, &tree, &synthesis);
        let (emitted_rust_compiles, emitted_rust_message, path) = compile_emitted(&source);
        assert!(
            emitted_rust_compiles,
            "the emitted module does not compile, so the pipeline did not emit Rust: {} \
             (kept at {})",
            emitted_rust_message,
            path.display()
        );

        // ── C1: the compiled expansion against the shipped table ─────────────
        let triangulations_identical = CASES
            .iter()
            .take(PATTERNS)
            .enumerate()
            .filter(|&(pattern, case)| {
                let leaf = generated_leaf(&values_for(pattern as u8));
                same_action(&GENERATED_LEAVES[leaf as usize], case)
            })
            .count();

        // ── C1's comparator must be able to fail ─────────────────────────────
        let mut mutant = LEAF_OF_CASE;
        mutant.swap(1, 2);
        let mutant_tree_mismatches = (0..PATTERNS)
            .filter(|&p| !same_action(&GENERATED_LEAVES[mutant[p] as usize], &CASES[p]))
            .count();

        // ── the emitted text against the compiled expansion, leaf for leaf ───
        let from_text = emitted_leaf_sequence(&source);
        let from_binary: Vec<u8> = (0..PATTERNS)
            .map(|t| {
                let mut pattern = 0u8;
                for (d, &condition) in CONDITION_ORDER.iter().enumerate() {
                    if (t >> (7 - d)) & 1 == 1 {
                        pattern |= 1 << condition;
                    }
                }
                generated_leaf(&values_for(pattern))
            })
            .collect();
        let emitted_matches_compiled = from_text.len() == PATTERNS
            && from_text == from_binary
            && from_text == synthesis.emitted_leaf_order;

        // ── C3: the four properties on the generated form ────────────────────
        let mut properties = Properties::default();
        for pattern in 0..PATTERNS {
            let leaf = generated_leaf(&values_for(pattern as u8));
            check_all_properties(
                pattern as u8,
                &GENERATED_LEAVES[leaf as usize],
                &mut properties,
            );
        }

        // ── C3's control: each property, individually, must be able to fail ──
        let mut sabotages_caught = 0usize;
        for s in sabotages() {
            let mut got = Properties::default();
            check_all_properties(s.case, &s.entry, &mut got);
            let caught = got.violations[s.property] > 0;
            println!(
                "  sabotage: {:<38} property {} violations {} -> {}",
                s.what,
                s.property + 1,
                got.violations[s.property],
                if caught { "CAUGHT" } else { "MISSED" }
            );
            if caught {
                sabotages_caught += 1;
            }
        }

        Sweep {
            synthesis,
            triangulations_identical,
            mutant_tree_mismatches,
            properties,
            sabotages_caught,
            emitted_bytes: source.len(),
            emitted_lines: source.lines().count(),
            emitted_rust_compiles,
            emitted_rust_message,
            emitted_matches_compiled,
            kani: run_kani(),
        }
    }

    // ─── private crate mechanisms, copied rather than made `pub` ────────────

    /// `cube::corner_offset`. Private, and `src/**` is read-only this phase.
    #[inline]
    const fn corner_offset(corner: u8) -> [u32; 3] {
        [
            (corner & 1) as u32,
            ((corner >> 1) & 1) as u32,
            ((corner >> 2) & 1) as u32,
        ]
    }

    /// `cube::place`: the centred frame.
    #[inline]
    fn place<R: Real>(lo: R, hi: R, d: R) -> R {
        (lo + hi) * R::HALF + (hi - lo) * d
    }

    /// `vec3::length`. Left-to-right summation, which is what makes the mirror's
    /// normals bit-identical to the crate's rather than merely close.
    #[inline]
    fn length<R: Real>(a: [R; 3]) -> R {
        (a[0] * a[0] + a[1] * a[1] + a[2] * a[2]).sqrt()
    }

    /// `marching_cubes::unit_gradient`, over `vec3::scale`.
    #[inline]
    fn unit_gradient<R: Real, S: Sdf<Scalar = R>>(sdf: &S, position: [R; 3]) -> [R; 3] {
        let g = sdf.gradient(position);
        let inv = length(g).recip();
        [g[0] * inv, g[1] * inv, g[2] * inv]
    }

    /// `marching_cubes::corner_position`.
    #[inline]
    fn corner_position<R: Real>(
        base: [u32; 3],
        corner: u8,
        origin: [R; 3],
        cell_size: R,
    ) -> [R; 3] {
        let o = corner_offset(corner);
        [
            origin[0] + cell_size * R::from_f64(f64::from(base[0] + o[0])),
            origin[1] + cell_size * R::from_f64(f64::from(base[1] + o[1])),
            origin[2] + cell_size * R::from_f64(f64::from(base[2] + o[2])),
        ]
    }

    /// `marching_cubes::edge_position` with `crossing_refinement == 0`, which is
    /// `MarchingCubes::new`'s default.
    #[inline]
    fn edge_position<R: Real>(
        base: [u32; 3],
        edge: u8,
        corner_value: &[R; 8],
        origin: [R; 3],
        cell_size: R,
    ) -> [R; 3] {
        let [lo_corner, hi_corner] = EDGE_CORNERS[edge as usize];
        let d = edge_offset(
            corner_value[lo_corner as usize],
            corner_value[hi_corner as usize],
        );
        let lo_pos = corner_position(base, lo_corner, origin, cell_size);
        let hi_pos = corner_position(base, hi_corner, origin, cell_size);
        [
            place(lo_pos[0], hi_pos[0], d),
            place(lo_pos[1], hi_pos[1], d),
            place(lo_pos[2], hi_pos[2], d),
        ]
    }

    // ─── the grid ──────────────────────────────────────────────────────────

    /// One cubic grid: `n` samples per axis, spanning `n − 1` cells.
    #[derive(Clone, Copy)]
    struct Grid<R: Real> {
        n: u32,
        origin: [R; 3],
        cell_size: R,
    }

    impl<R: Real> Grid<R> {
        #[inline]
        fn samples(self) -> usize {
            let n = self.n as usize;
            n * n * n
        }

        #[inline]
        fn cells(self) -> u32 {
            self.n - 1
        }

        #[inline]
        fn cell_count(self) -> usize {
            let c = self.cells() as usize;
            c * c * c
        }

        /// `RuntimeShape3::linearize`, which is the layout Marching Cubes
        /// samples into: no row padding, so the stride is the row itself.
        #[inline]
        fn sample_index(self, p: [u32; 3]) -> usize {
            let n = self.n as usize;
            p[0] as usize + n * (p[1] as usize + n * p[2] as usize)
        }
    }

    /// The eight corner values for one cell — `mod.rs:259-268`'s gather, in
    /// **one** place, called by both arms, so neither has a cheaper way to load
    /// a cell than the other.
    #[inline]
    fn gather<R: Real>(values: &[R], g: Grid<R>, base: [u32; 3]) -> [R; 8] {
        let mut corner_value = [R::ZERO; 8];
        for (c, slot) in corner_value.iter_mut().enumerate() {
            let o = corner_offset(c as u8);
            *slot = values[g.sample_index([base[0] + o[0], base[1] + o[1], base[2] + o[2]])];
        }
        corner_value
    }

    // ─── the march, in two dispatch mechanisms ─────────────────────────────

    /// Both arms, sharing every buffer and sharing [`Self::emit_cell`], so the
    /// payload code is literally the same instructions in both.
    struct March<R: Real> {
        /// Filled once per row, **outside** every counter window.
        values: Vec<R>,
        /// `MarchingCubes::edge_vertices`: one `u32` slot per (sample, axis).
        edge_vertices: Vec<u32>,
        positions: Vec<[R; 3]>,
        normals: Vec<[R; 3]>,
        indices: Vec<u32>,
    }

    impl<R: Real> March<R> {
        fn new() -> Self {
            Self {
                values: Vec::new(),
                edge_vertices: Vec::new(),
                positions: Vec::new(),
                normals: Vec::new(),
                indices: Vec::new(),
            }
        }

        /// `sdf::sample_grid` with `row_stride == size[0]`. Outside every window.
        fn sample<S: Sdf<Scalar = R>>(&mut self, sdf: &S, g: Grid<R>) {
            self.values.clear();
            self.values.reserve(g.samples());
            for z in 0..g.n {
                for y in 0..g.n {
                    for x in 0..g.n {
                        self.values.push(sdf.sample([
                            g.origin[0] + g.cell_size * R::from_f64(f64::from(x)),
                            g.origin[1] + g.cell_size * R::from_f64(f64::from(y)),
                            g.origin[2] + g.cell_size * R::from_f64(f64::from(z)),
                        ]));
                    }
                }
            }
        }

        /// `mod.rs:250-251`, and it runs inside both arms' windows because it
        /// runs on every shipped extract.
        #[inline]
        fn prepare(&mut self, g: Grid<R>) {
            self.edge_vertices.clear();
            self.edge_vertices.resize(g.samples() * 3, u32::MAX);
            self.positions.clear();
            self.normals.clear();
            self.indices.clear();
        }

        /// `MeshSink::vertex`.
        #[inline]
        fn vertex<S: Sdf<Scalar = R>>(&mut self, sdf: &S, position: [R; 3]) -> u32 {
            let index = self.positions.len() as u32;
            self.positions.push(position);
            self.normals.push(unit_gradient(sdf, position));
            index
        }

        /// `MarchingCubes::vertex_on_edge` — the edge-cache probe.
        #[inline]
        fn vertex_on_edge<S: Sdf<Scalar = R>>(
            &mut self,
            sdf: &S,
            g: Grid<R>,
            base: [u32; 3],
            edge: u8,
            corner_value: &[R; 8],
        ) -> u32 {
            let axis = EDGE_AXIS[edge as usize] as usize;
            let o = corner_offset(EDGE_CORNERS[edge as usize][0]);
            let lo_sample = g.sample_index([base[0] + o[0], base[1] + o[1], base[2] + o[2]]);
            let key = lo_sample * 3 + axis;
            let cached = self.edge_vertices[key];
            if cached != u32::MAX {
                return cached;
            }
            let position = edge_position(base, edge, corner_value, g.origin, g.cell_size);
            let index = self.vertex(sdf, position);
            self.edge_vertices[key] = index;
            index
        }

        /// **The payload for one cell**, `mod.rs:307-377` under
        /// `MarchingCubes::new`'s defaults — `FaceAmbiguity::Separate`,
        /// `InteriorAmbiguity::Ignore`, `crossing_refinement == 0`.
        ///
        /// Takes the **entry**, not the case index, which is the whole point:
        /// both arms reach an `&McCase` and then run identical code, so the only
        /// thing under comparison is how the entry was reached.
        #[inline]
        fn emit_cell<S: Sdf<Scalar = R>>(
            &mut self,
            sdf: &S,
            g: Grid<R>,
            base: [u32; 3],
            entry: &McCase,
            corner_value: &[R; 8],
        ) {
            if entry.count == 0 {
                return;
            }
            // Cycle centroids first (A-015). Cell-local, so never cached.
            let mut centroid = [0u32; MAX_CENTROIDS];
            for (c, slot) in centroid
                .iter_mut()
                .enumerate()
                .take(entry.centroids as usize)
            {
                let code = CENTROID_BASE + c as u8;
                let mut sum = [R::ZERO; 3];
                let mut n = 0u32;
                for tri in &entry.triangles[..entry.count as usize] {
                    if tri[0] != code {
                        continue;
                    }
                    let p = edge_position(base, tri[1], corner_value, g.origin, g.cell_size);
                    sum = [sum[0] + p[0], sum[1] + p[1], sum[2] + p[2]];
                    n += 1;
                }
                let scale = R::from_f64(f64::from(n)).recip();
                let position = [sum[0] * scale, sum[1] * scale, sum[2] * scale];
                *slot = self.vertex(sdf, position);
            }
            for tri in &entry.triangles[..entry.count as usize] {
                let mut idx = [0u32; 3];
                for (k, &code) in tri.iter().enumerate() {
                    idx[k] = if is_centroid(code) {
                        centroid[(code - CENTROID_BASE) as usize]
                    } else {
                        self.vertex_on_edge(sdf, g, base, code, corner_value)
                    };
                }
                self.indices.extend_from_slice(&idx);
            }
        }

        /// **The table arm.** Eight branchless `case |= 1 << c`, one indexed
        /// load out of the shipped `static CASES`.
        fn march_table<S: Sdf<Scalar = R>>(&mut self, sdf: &S, g: Grid<R>) {
            self.prepare(g);
            let c = g.cells();
            for z in 0..c {
                for y in 0..c {
                    for x in 0..c {
                        let base = [x, y, z];
                        let corner_value = gather(&self.values, g, base);
                        let entry = &CASES[case_of(&corner_value) as usize];
                        self.emit_cell(sdf, g, base, entry, &corner_value);
                    }
                }
            }
            black_box(&self.indices);
        }

        /// **The generated arm.** The emitted eight-deep branch chain, one
        /// indexed load out of the pipeline's own leaf array.
        fn march_generated<S: Sdf<Scalar = R>>(&mut self, sdf: &S, g: Grid<R>) {
            self.prepare(g);
            let c = g.cells();
            for z in 0..c {
                for y in 0..c {
                    for x in 0..c {
                        let base = [x, y, z];
                        let corner_value = gather(&self.values, g, base);
                        let entry = &GENERATED_LEAVES[generated_leaf(&corner_value) as usize];
                        self.emit_cell(sdf, g, base, entry, &corner_value);
                    }
                }
            }
            black_box(&self.indices);
        }

        /// Cells the case table gives triangles for. Outside every window.
        fn active_cells(&self, g: Grid<R>) -> usize {
            let c = g.cells();
            let mut n = 0usize;
            for z in 0..c {
                for y in 0..c {
                    for x in 0..c {
                        let corner_value = gather(&self.values, g, [x, y, z]);
                        if CASES[case_of(&corner_value) as usize].count != 0 {
                            n += 1;
                        }
                    }
                }
            }
            n
        }
    }

    // ─── counting ──────────────────────────────────────────────────────────

    /// What one or more counter windows read.
    #[derive(Clone, Copy, Default)]
    struct Counted {
        cycles: f64,
        instructions: f64,
        branch_misses: f64,
    }

    /// One counter window over `inner` repetitions of `body`, divided by
    /// `inner`.
    ///
    /// The `perf_event` system calls are all **outside** the counted region, and
    /// windows are never nested: `Probe` opens six hardware counters and Zen 3
    /// has six, so a nested pair would multiplex and `worst_ratio` refuses —
    /// `R-121` paid for that discovery.
    fn window(probe: &mut Probe, inner: usize, mut body: impl FnMut()) -> (Counted, f64) {
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
        (
            Counted {
                cycles: counts.cycles.count as f64 * scale,
                instructions: counts.instructions.count as f64 * scale,
                branch_misses: counts.branch_misses.count as f64 * scale,
            },
            nanos * scale,
        )
    }

    /// One repetition: one sibling window per arm, in one of the two orderings.
    #[derive(Clone, Copy, Default)]
    struct Rep {
        table: Counted,
        generated: Counted,
        table_ns: f64,
        generated_ns: f64,
        /// Which arm ran first. Recorded rather than assumed away: the second
        /// arm inherits the first's cache and predictor state.
        table_first: bool,
    }

    /// The median of a set of readings, taken **per quantity** rather than per
    /// repetition: one repetition disturbed by another process on the machine
    /// should move one number, not a whole row.
    fn median(pick: &dyn Fn(&Rep) -> f64, reps: &[Rep]) -> f64 {
        let mut values: Vec<f64> = reps.iter().map(pick).collect();
        values.sort_by(|a, b| a.total_cmp(b));
        values[values.len() / 2]
    }

    fn median_counted(pick: &dyn Fn(&Rep) -> Counted, reps: &[Rep]) -> Counted {
        Counted {
            cycles: median(&|r| pick(r).cycles, reps),
            instructions: median(&|r| pick(r).instructions, reps),
            branch_misses: median(&|r| pick(r).branch_misses, reps),
        }
    }

    // ─── one row ───────────────────────────────────────────────────────────

    /// One `(field, 65³)` measurement.
    struct Measured {
        field: &'static str,
        inner: usize,
        cells: usize,
        active_cells: usize,
        vertices: usize,
        triangles: usize,
        mesh_identical: bool,
        mesh_identical_to_shipped: bool,
        table: Counted,
        generated: Counted,
        table_ns: f64,
        generated_ns: f64,
        instruction_ratio_spread: f64,
    }

    impl Measured {
        fn per_cell(&self, pick: impl Fn(Counted) -> f64, generated: bool) -> f64 {
            let c = if generated {
                self.generated
            } else {
                self.table
            };
            pick(c) / self.cells as f64
        }

        fn instruction_ratio(&self) -> f64 {
            self.generated.instructions / self.table.instructions
        }

        fn cycle_ratio(&self) -> f64 {
            self.generated.cycles / self.table.cycles
        }

        fn branch_miss_ratio(&self) -> f64 {
            self.generated.branch_misses / self.table.branch_misses
        }

        fn ghz(&self) -> f64 {
            self.table.cycles / self.table_ns
        }

        fn c2_holds(&self) -> bool {
            self.instruction_ratio() <= INSTRUCTION_BAR
                && self.mesh_identical
                && self.mesh_identical_to_shipped
        }
    }

    /// Bit for bit, as bit patterns rather than as values.
    fn same<R: Real>(a: &[[R; 3]], b: &[[R; 3]]) -> bool {
        a.len() == b.len()
            && a.iter()
                .zip(b)
                .all(|(p, q)| (0..3).all(|k| p[k].as_f64().to_bits() == q[k].as_f64().to_bits()))
    }

    /// Measure one `(field, 65³)`.
    fn measure<R, S>(field: &'static str, n: u32, sdf: &S, origin: [R; 3], cell_size: R) -> Measured
    where
        R: Real,
        S: Sdf<Scalar = R>,
    {
        let g = Grid {
            n,
            origin,
            cell_size,
        };
        let shape = RuntimeShape3::new([n; 3]).expect("the fixture fits u32");

        let mut m = March::<R>::new();
        m.sample(sdf, g);
        let active = m.active_cells(g);

        for _ in 0..WARMUP {
            m.march_table(sdf, g);
            m.march_generated(sdf, g);
        }

        let started = Instant::now();
        m.march_table(sdf, g);
        let pass_ns = started.elapsed().as_nanos() as f64;
        let inner = ((TARGET_BATCH_NS / pass_ns.max(1.0)).ceil() as usize).clamp(1, MAX_INNER);

        // ── the mirror agreement that licenses the ratio (M-279) ─────────────
        m.march_table(sdf, g);
        let table_positions = m.positions.clone();
        let table_normals = m.normals.clone();
        let table_indices = m.indices.clone();
        m.march_generated(sdf, g);
        let mesh_identical = same(&table_positions, &m.positions)
            && same(&table_normals, &m.normals)
            && table_indices.as_slice() == m.indices.as_slice();
        assert!(
            mesh_identical,
            "{field} {n}^3: the generated arm's mesh differs from the table arm's ({} vs {} \
             vertices, {} vs {} indices) -- C1 is an equality over 256 patterns and a mesh that \
             disagrees means the sweep missed a pattern the fixture reaches",
            table_positions.len(),
            m.positions.len(),
            table_indices.len(),
            m.indices.len()
        );

        let mut shipped = MarchingCubes::<R>::new();
        let mut out = MeshBuffer::<R>::new();
        shipped
            .extract(sdf, &shape, origin, cell_size, &mut out)
            .expect("extraction");
        let mesh_identical_to_shipped = same(&table_positions, &out.positions)
            && same(&table_normals, &out.normals)
            && table_indices.as_slice() == out.indices.as_slice();
        assert!(
            mesh_identical_to_shipped,
            "{field} {n}^3: the table arm's mesh differs from MarchingCubes::extract's ({} vs {} \
             vertices, {} vs {} indices) -- then it is not the shipped shape and C2 compares two \
             things neither of which ships",
            table_positions.len(),
            out.positions.len(),
            table_indices.len(),
            out.indices.len()
        );

        // ── REPS repetitions, one sibling window per arm ─────────────────────
        let mut probe = Probe::open();
        let mut reps: Vec<Rep> = Vec::with_capacity(REPS);
        for rep in 0..REPS {
            let mut r = Rep {
                table_first: rep % 2 == 0,
                ..Rep::default()
            };
            if r.table_first {
                let (c, ns) = window(&mut probe, inner, || m.march_table(sdf, g));
                r.table = c;
                r.table_ns = ns;
                let (c, ns) = window(&mut probe, inner, || m.march_generated(sdf, g));
                r.generated = c;
                r.generated_ns = ns;
            } else {
                let (c, ns) = window(&mut probe, inner, || m.march_generated(sdf, g));
                r.generated = c;
                r.generated_ns = ns;
                let (c, ns) = window(&mut probe, inner, || m.march_table(sdf, g));
                r.table = c;
                r.table_ns = ns;
            }
            reps.push(r);
        }

        let ratios: Vec<f64> = reps
            .iter()
            .map(|r| r.generated.instructions / r.table.instructions)
            .collect();
        let instruction_ratio_spread = ratios.iter().copied().fold(0.0f64, f64::max)
            - ratios.iter().copied().fold(f64::MAX, f64::min);

        let table = median_counted(&|r| r.table, &reps);
        let generated = median_counted(&|r| r.generated, &reps);
        assert!(
            table.instructions > 0.0 && generated.instructions > 0.0,
            "{field} {n}^3: an arm read zero instructions, so the ratio is a division by a floor"
        );

        Measured {
            field,
            inner,
            cells: g.cell_count(),
            active_cells: active,
            vertices: table_positions.len(),
            triangles: table_indices.len() / 3,
            mesh_identical,
            mesh_identical_to_shipped,
            table,
            generated,
            table_ns: median(&|r| r.table_ns, &reps),
            generated_ns: median(&|r| r.generated_ns, &reps),
            instruction_ratio_spread,
        }
    }

    /// The registered C2 fixture: eight reference fields × 65³.
    fn sweep_fields() -> Vec<Measured> {
        let mut rows = Vec::new();
        isomesh::for_each_reference_field!(f32, |name, field| {
            let (_, origin, cell_size) = crate::common::grid(&field, RESOLUTION);
            rows.push(measure(name, RESOLUTION, &field, origin, cell_size));
        });
        rows
    }

    pub(crate) fn run(run: &mut crate::common::experiment::Run) {
        // ── stages 1 to 4, and C1 and C3 ─────────────────────────────────────
        println!("GRAPHGEN's pipeline over the 256-entry case table\n");
        let s = sweep_patterns();
        let syn = &s.synthesis;

        println!(
            "\n  decision table   {PATTERNS} rows, 8 boolean conditions, 0 don't-cares, {} \
             distinct actions ({} patterns share one)",
            syn.distinct_actions, syn.shared_action_patterns
        );
        println!(
            "  optimal tree     {} nodes ({} internal + {} leaves), average path length {:.6}",
            syn.tree_nodes(),
            syn.tree_internal,
            syn.tree_leaves,
            syn.average_path_len()
        );
        println!(
            "  objective        min {:.6} / max {:.6} average path length over {} memoised \
             subcube states",
            syn.min_total_path_len as f64 / PATTERNS as f64,
            syn.max_total_path_len as f64 / PATTERNS as f64,
            syn.dp_states
        );
        println!(
            "  DRAG             {} nodes ({} internal + {} leaves), {} node(s) merged, {:.4}% of \
             the tree",
            syn.drag_nodes(),
            syn.drag_internal,
            syn.drag_leaves,
            syn.tree_nodes() - syn.drag_nodes(),
            100.0 * syn.drag_nodes() as f64 / syn.tree_nodes() as f64
        );
        println!(
            "  emitted Rust     {} bytes, {} lines, rustc --emit metadata -> {}, text matches \
             the compiled expansion -> {}",
            s.emitted_bytes, s.emitted_lines, s.emitted_rust_compiles, s.emitted_matches_compiled
        );
        println!(
            "  C1               {}/{} triangulations identical to CASES; mutant comparator caught \
             {} mismatch(es)",
            s.triangulations_identical, PATTERNS, s.mutant_tree_mismatches
        );
        println!(
            "  C3               {}/4 properties held over {} checks on {PATTERNS} patterns, {} \
             violation(s); {}/4 sabotages caught",
            s.properties.held(),
            s.properties.checks,
            s.properties.total_violations(),
            s.sabotages_caught
        );
        println!(
            "  Kani             {} -- {} over {} checks ({} failed, {:.4}s); sabotage {} over {} \
             checks ({} failed)",
            s.kani.version,
            s.kani.proof.status,
            s.kani.proof.checks,
            s.kani.proof.failed,
            s.kani.proof.solver_seconds,
            s.kani.sabotage.status,
            s.kani.sabotage.checks,
            s.kani.sabotage.failed
        );

        // ── the registered vacuity controls, asserted ────────────────────────
        assert_eq!(
            s.triangulations_identical,
            PATTERNS,
            "C1 is an equality over the whole input space and {} of {PATTERNS} patterns disagree \
             with CASES",
            PATTERNS - s.triangulations_identical
        );
        assert!(
            s.emitted_matches_compiled,
            "the emitted text's leaf sequence is not the compiled expansion's, so the measured \
             form is not the pipeline's output"
        );
        assert!(
            s.mutant_tree_mismatches > 0,
            "VOID: two leaf ids were swapped and the C1 comparator saw nothing, so C1 is an \
             equality between two names for one computation"
        );
        assert!(
            s.properties.checks > 0,
            "VOID: the property sweep made ZERO checks. Four properties holding over an empty \
             check set is M-44's vacuous zero (experiment_p64.rs:169-183)"
        );
        assert_eq!(
            s.sabotages_caught,
            4,
            "VOID: {} of 4 sabotages went uncaught, so at least one of the four properties \
             cannot fail and the verdict beside it is a tautology",
            4 - s.sabotages_caught
        );

        // ── C2 ───────────────────────────────────────────────────────────────
        let rows = sweep_fields();
        println!(
            "\n{:<16} {:>5} {:>12} {:>12} {:>9} {:>8} {:>12} {:>12} {:>9} {:>10} {:>10} {:>7} \
             {:>6} {:>4}",
            "field",
            "inner",
            "instr/cell",
            "instr/cell",
            "instr",
            "spread",
            "bmiss/cell",
            "bmiss/cell",
            "bmiss",
            "cyc/cell",
            "cyc/cell",
            "cycle",
            "ghz",
            "C2"
        );
        println!(
            "{:<16} {:>5} {:>12} {:>12} {:>9} {:>8} {:>12} {:>12} {:>9} {:>10} {:>10} {:>7} \
             {:>6} {:>4}",
            "",
            "",
            "table",
            "generated",
            "ratio",
            "",
            "table",
            "generated",
            "ratio",
            "table",
            "generated",
            "ratio",
            "",
            ""
        );
        for r in &rows {
            println!(
                "{:<16} {:>5} {:>12.4} {:>12.4} {:>9.6} {:>8.6} {:>12.6} {:>12.6} {:>9.4} \
                 {:>10.4} {:>10.4} {:>7.4} {:>6.3} {:>4}",
                r.field,
                r.inner,
                r.per_cell(|c| c.instructions, false),
                r.per_cell(|c| c.instructions, true),
                r.instruction_ratio(),
                r.instruction_ratio_spread,
                r.per_cell(|c| c.branch_misses, false),
                r.per_cell(|c| c.branch_misses, true),
                r.branch_miss_ratio(),
                r.per_cell(|c| c.cycles, false),
                r.per_cell(|c| c.cycles, true),
                r.cycle_ratio(),
                r.ghz(),
                if r.c2_holds() { "HELD" } else { "FALS" }
            );
        }

        let c1 = s.c1_holds();
        let c3 = s.c3_holds();
        let c2_all = rows.iter().all(Measured::c2_holds);
        let c2_rows = rows.iter().filter(|r| r.c2_holds()).count();
        let worst = rows
            .iter()
            .map(Measured::instruction_ratio)
            .fold(0.0f64, f64::max);
        let best = rows
            .iter()
            .map(Measured::instruction_ratio)
            .fold(f64::MAX, f64::min);

        println!(
            "\nC1 the emitted DRAG's triangulations identical to CASES on all 256 patterns -> {}",
            if c1 { "HELD" } else { "FALSIFIED" }
        );
        println!(
            "C2 instructions per cell not above the table lookup's (bar {INSTRUCTION_BAR:.2}) -> \
             {} on {c2_rows}/{} rows, ratio {best:.6}-{worst:.6}",
            if c2_all { "HELD" } else { "FALSIFIED" },
            rows.len()
        );
        println!(
            "C3 P-64's four properties hold on the generated form -> {}",
            if c3 { "HELD" } else { "FALSIFIED" }
        );
        println!(
            "\nINSTRUMENT, per the registration: C3's verdict is the EXHAUSTIVE 256-pattern Rust \
             check above, because crates/isomesh/src/** is read-only this phase and \
             marching_cubes::proofs is cfg(kani) inside it -- a #[kani::proof] cannot be placed \
             over a bench-local form. Kani's {} checks are a proof about the SHIPPED CASES and \
             reach the generated form only through C1's pointwise identity.",
            s.kani.proof.checks
        );
        println!(
            "PREMISE, per the registration: CASES is already const-derived by build_cases() at \
             table.rs:182-194, so the mistranscription risk this row was motivated by is already \
             mitigated in the shipped path. The Bourke tables are a test cross-check only."
        );

        for r in &rows {
            run.record(&[
                // ── the registration's columns ───────────────────────────────
                ("patterns_tested", PATTERNS.to_string()),
                (
                    "triangulations_identical",
                    s.triangulations_identical.to_string(),
                ),
                ("tree_nodes", syn.tree_nodes().to_string()),
                ("drag_nodes", syn.drag_nodes().to_string()),
                ("average_path_len", format!("{:.6}", syn.average_path_len())),
                (
                    "instructions_per_cell_table",
                    format!("{:.4}", r.per_cell(|c| c.instructions, false)),
                ),
                (
                    "instructions_per_cell_generated",
                    format!("{:.4}", r.per_cell(|c| c.instructions, true)),
                ),
                ("ratio", format!("{:.6}", r.instruction_ratio())),
                ("properties_held", s.properties.held().to_string()),
                ("sabotage_failed", (s.sabotages_caught == 4).to_string()),
                ("kani_checks", s.kani.proof.checks.to_string()),
                ("c1_holds", c1.to_string()),
                ("c2_holds", r.c2_holds().to_string()),
                ("c3_holds", c3.to_string()),
                // ── the fixture ──────────────────────────────────────────────
                ("field", r.field.to_string()),
                ("resolution", RESOLUTION.to_string()),
                ("cells", r.cells.to_string()),
                ("active_cells", r.active_cells.to_string()),
                ("vertices", r.vertices.to_string()),
                ("triangles", r.triangles.to_string()),
                // ── stages 1 to 3 ────────────────────────────────────────────
                ("distinct_actions", syn.distinct_actions.to_string()),
                (
                    "shared_action_patterns",
                    syn.shared_action_patterns.to_string(),
                ),
                ("dont_care_entries", 0.to_string()),
                ("tree_internal_nodes", syn.tree_internal.to_string()),
                ("tree_leaf_nodes", syn.tree_leaves.to_string()),
                ("drag_internal_nodes", syn.drag_internal.to_string()),
                ("drag_leaf_nodes", syn.drag_leaves.to_string()),
                (
                    "drag_nodes_merged",
                    (syn.tree_nodes() - syn.drag_nodes()).to_string(),
                ),
                (
                    "drag_compression",
                    format!("{:.6}", syn.drag_nodes() as f64 / syn.tree_nodes() as f64),
                ),
                (
                    "min_average_path_len",
                    format!("{:.6}", syn.min_total_path_len as f64 / PATTERNS as f64),
                ),
                (
                    "max_average_path_len",
                    format!("{:.6}", syn.max_total_path_len as f64 / PATTERNS as f64),
                ),
                ("dp_subcube_states", syn.dp_states.to_string()),
                (
                    "condition_order",
                    CONDITION_ORDER
                        .iter()
                        .map(usize::to_string)
                        .collect::<String>(),
                ),
                // ── stage 4 ──────────────────────────────────────────────────
                ("emitted_bytes", s.emitted_bytes.to_string()),
                ("emitted_lines", s.emitted_lines.to_string()),
                ("emitted_rust_compiles", s.emitted_rust_compiles.to_string()),
                ("emitted_rustc_message", s.emitted_rust_message.clone()),
                (
                    "emitted_matches_compiled",
                    s.emitted_matches_compiled.to_string(),
                ),
                // ── controls ─────────────────────────────────────────────────
                (
                    "mutant_tree_mismatches",
                    s.mutant_tree_mismatches.to_string(),
                ),
                ("property_checks", s.properties.checks.to_string()),
                (
                    "property_violations",
                    s.properties.total_violations().to_string(),
                ),
                ("sabotages_caught", s.sabotages_caught.to_string()),
                ("mesh_identical", r.mesh_identical.to_string()),
                (
                    "mesh_identical_to_shipped",
                    r.mesh_identical_to_shipped.to_string(),
                ),
                // ── Kani, corroboration only ─────────────────────────────────
                ("kani_version", s.kani.version.clone()),
                ("kani_status", s.kani.proof.status.clone()),
                ("kani_failed_checks", s.kani.proof.failed.to_string()),
                (
                    "kani_solver_seconds",
                    format!("{:.4}", s.kani.proof.solver_seconds),
                ),
                ("kani_sabotage_status", s.kani.sabotage.status.clone()),
                ("kani_sabotage_checks", s.kani.sabotage.checks.to_string()),
                (
                    "kani_sabotage_failed_checks",
                    s.kani.sabotage.failed.to_string(),
                ),
                // ── C2's provenance, none of which carries a verdict ─────────
                (
                    "cycles_per_cell_table",
                    format!("{:.4}", r.per_cell(|c| c.cycles, false)),
                ),
                (
                    "cycles_per_cell_generated",
                    format!("{:.4}", r.per_cell(|c| c.cycles, true)),
                ),
                ("cycle_ratio", format!("{:.6}", r.cycle_ratio())),
                (
                    "branch_misses_per_cell_table",
                    format!("{:.6}", r.per_cell(|c| c.branch_misses, false)),
                ),
                (
                    "branch_misses_per_cell_generated",
                    format!("{:.6}", r.per_cell(|c| c.branch_misses, true)),
                ),
                ("branch_miss_ratio", format!("{:.6}", r.branch_miss_ratio())),
                (
                    "instruction_ratio_spread",
                    format!("{:.6}", r.instruction_ratio_spread),
                ),
                (
                    "ns_per_cell_table",
                    format!("{:.6}", r.table_ns / r.cells as f64),
                ),
                (
                    "ns_per_cell_generated",
                    format!("{:.6}", r.generated_ns / r.cells as f64),
                ),
                ("ghz", format!("{:.4}", r.ghz())),
                ("inner", r.inner.to_string()),
                ("reps", REPS.to_string()),
            ]);
        }
    }
}

fn main() {
    if !std::env::args().any(|arg| arg == "--bench") {
        return;
    }

    let prereg = isomesh::experiment!("P-116");

    #[cfg(target_os = "linux")]
    common::experiment::run(prereg, experiment::run);

    // `experiment_p12`'s precedent. C2 is instructions per cell, which comes
    // from `perf_event_open`; a recorded zero would be a fabricated measurement
    // and a clock cannot substitute for an instruction count (`M-281`).
    #[cfg(not(target_os = "linux"))]
    {
        eprintln!(
            "{} scores C2 on instructions per cell from hardware performance counters, and this \
             platform has no `perf_event_open`. There is no substitute.",
            prereg.id
        );
        std::process::exit(1);
    }
}
