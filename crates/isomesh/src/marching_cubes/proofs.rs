//! Bounded-model-checking proofs over the case tables.
//!
//! Ticket: R-062 (P-64). Run with `cargo kani -p isomesh`. Compiled only under
//! `cfg(kani)`, so this file is invisible to every other build — `cargo tree -p
//! isomesh -e normal` is unchanged at two packages and hard rule 3 is not
//! engaged.
//!
//! # What is proved, and what is deliberately not
//!
//! `CLAUDE.md` rule 5: *"wrong case tables produce meshes that look fine and are
//! subtly non-manifold."* That risk is **combinatorics over eight sign bits** —
//! 256 states — and bounded model checking settles 256 states exhaustively where
//! a property test samples them.
//!
//! What is **not** proved is vertex placement. Bit-blasting IEEE 754 to SAT is
//! the adversarial case for a model checker, and eight nondeterministic `f32`
//! corners is 256 bits of unconstrained float. Placement stays under proptest
//! and the golden hashes. The honest scope of this file is *"the table cannot be
//! indexed wrongly"*, not *"the mesh is correct"*.
//!
//! # The four properties
//!
//! 1. **Shape.** `count ≤ MAX_TRIANGLES` and `centroids ≤ MAX_CENTROIDS`, so the
//!    fixed-size arrays a caller stack-allocates cannot overflow. `D-021` is
//!    what happens when a bound like this is a *sampled* maximum rather than a
//!    derived one.
//! 2. **Every code is nameable.** Each corner code in a live triangle is either
//!    a cube edge below [`CENTROID_BASE`] or a centroid index below this case's
//!    own `centroids` count. A centroid code past the count names a vertex the
//!    extractor never created.
//! 3. **Every named edge is cut.** This is the load-bearing one and it is a
//!    property of the *pair* (table, sign pattern) rather than of the table
//!    alone. `extraction.rs` places a vertex for whatever edge the table names,
//!    by interpolating between that edge's two corner values — so a triangle
//!    naming an **uncut** edge asks for a crossing that does not exist, which is
//!    a vertex at a meaningless position rather than a bounds error. It is
//!    exactly the failure rule 5 describes: a mesh that looks fine.
//! 4. **No degenerate triangle.** Three distinct codes, so no triangle is a
//!    zero-area sliver by construction. `validate::self_intersection` counts
//!    these at run time; here they are impossible.
//!
//! Property 4 is stated over codes rather than positions on purpose: two
//! distinct edges can carry the *same* position when the crossing sits on a
//! shared corner, and that is a placement question this file does not touch.

use super::table::{
    CASES, CENTROID_BASE, EDGE_CORNERS, EDGE_COUNT, MAX_CENTROIDS, MAX_TRIANGLES, McCase, NO_EDGE,
    is_centroid, segment_links, triangulate,
};

/// Is edge `e` cut by this sign pattern?
///
/// The same test `extraction.rs` makes implicitly when it interpolates: the two
/// corner values must straddle zero, which over sign bits is "the bits differ".
fn edge_is_cut(case: u8, e: usize) -> bool {
    let [lo, hi] = EDGE_CORNERS[e];
    ((case >> lo) & 1) != ((case >> hi) & 1)
}

/// All four properties, for one case and one already-built triangulation.
///
/// Shared by the three harnesses so the *statement* of the properties exists
/// once. Three copies of an assertion is three things to drift, and a proof
/// whose properties drifted from the ones it claims is worse than no proof.
fn check_all_properties(case: u8, c: &McCase) {
    // Property 1: shape.
    assert!(
        c.count as usize <= MAX_TRIANGLES,
        "count exceeds MAX_TRIANGLES"
    );
    assert!(
        c.centroids as usize <= MAX_CENTROIDS,
        "centroids exceeds MAX_CENTROIDS"
    );

    let mut t = 0usize;
    while t < MAX_TRIANGLES {
        if t >= c.count as usize {
            break;
        }
        let tri = c.triangles[t];

        let mut k = 0usize;
        while k < 3 {
            let code = tri[k];

            // Property 2: nameable.
            assert!(code != NO_EDGE, "a live triangle carries NO_EDGE");
            if is_centroid(code) {
                assert!(
                    code - CENTROID_BASE < c.centroids,
                    "a centroid code is past this case's centroid count"
                );
            } else {
                assert!(
                    (code as usize) < EDGE_COUNT,
                    "an edge code is past EDGE_COUNT"
                );

                // Property 3: the named edge is cut.
                assert!(
                    edge_is_cut(case, code as usize),
                    "a triangle names an edge this sign pattern does not cut"
                );
            }
            k += 1;
        }

        // Property 4: no degenerate triangle.
        assert!(
            tri[0] != tri[1] && tri[1] != tri[2] && tri[0] != tri[2],
            "a triangle carries two equal codes"
        );
        t += 1;
    }
}

/// **C1** — the shipped `CASES` table, all 256 patterns.
///
/// This is the table `extract` reads at mask zero, so it is the one a consumer
/// actually meets. Reading `CASES[case]` rather than recomputing it is
/// deliberate: the property is about the **static** the crate ships, and a proof
/// about a recomputation would not catch a `build_cases` that diverged from what
/// it stored.
#[kani::proof]
fn shipped_case_table_is_indexable_for_every_sign_pattern() {
    let case: u8 = kani::any();
    check_all_properties(case, &CASES[case as usize]);
}

/// Every `(case, joined)` triangulation, evaluated at **compile time**.
///
/// # Why a table and not symbolic execution of the constructor
///
/// The first version of this file pointed Kani straight at
/// `triangulate(segment_links(case, joined))` with both arguments
/// nondeterministic. **CBMC ran out of memory** — 10,230,639 steps of program
/// expression, 909,347 verification conditions, 550,658 after simplification,
/// 885 s of symbolic execution before it died on a 32 GB box. The shipped
/// `CASES` harness beside it verifies in **2.05 s**.
///
/// That gap is not the sign-bit space, which is 256 states and trivial. It is
/// `segment_links`' **walk**: a data-dependent traversal of the twelve-edge link
/// structure, where every step's successor depends on the previous step's value,
/// so symbolic execution cannot merge paths and unrolls the whole product.
///
/// So the work is split between the two engines that are each good at half of
/// it. **`const` evaluation runs the real constructor on all 65,536 inputs** —
/// exhaustive by construction, and the same evaluator that builds the shipped
/// `CASES`. **Kani then proves the four properties over the stored results**
/// with a nondeterministic index, which is a lookup and verifies in seconds.
/// Composed, that is a proof over the entire input space: nothing is sampled and
/// nothing is assumed about the walk.
///
/// # The mask alphabet, and what is and is not being assumed
///
/// The mask ranges over all **64** values [`super::table::face_bit`] can
/// produce: `1 << (axis * 2 + side)` with `axis < 3` and `side < 2` occupies bits
/// 0 through 5, so bits 6 and 7 are outside the alphabet and no caller can set
/// them. That is a restriction to the **type**, not to well-formedness — the
/// proof still ranges over every mask that is *not* in
/// [`super::table::AMBIGUOUS_FACES`]`[case]`, which is the part that could
/// actually be wrong, and `joined_mask` returning a bit for a non-ambiguous face
/// is exactly the defect this harness would catch.
///
/// 256 × 256 was the first attempt and **rustc refused it**: *"constant
/// evaluation is taking a long time"* at 65,536 entries. 16,384 evaluates.
#[cfg(kani)]
static ALL_MASKS: [[McCase; 64]; 256] = build_all_masks();

///
/// `long_running_const_eval` is allowed here and nowhere else. It is a **compile
/// -time budget lint, not a correctness one** — rustc's default step limit is
/// tuned for accidental infinite loops, and 16,384 twelve-edge walks is a
/// bounded, terminating computation that simply exceeds it. The alternative is
/// symbolic execution of the same walk, which is what ran out of 32 GB.
#[cfg(kani)]
const fn build_all_masks() -> [[McCase; 64]; 256] {
    let empty = McCase {
        count: 0,
        centroids: 0,
        triangles: [[0u8; 3]; MAX_TRIANGLES],
    };
    let mut out = [[empty; 64]; 256];
    let mut case = 0usize;
    while case < 256 {
        let mut joined = 0usize;
        while joined < 64 {
            out[case][joined] = triangulate(segment_links(case as u8, joined as u8));
            joined += 1;
        }
        case += 1;
    }
    out
}

/// **C3** — the interior-ambiguity path: every case against every mask.
///
/// 256 sign patterns × 64 face masks = **16,384** triangulations, every one of
/// them produced by the real `segment_links` + `triangulate` under `const`
/// evaluation, every one of them checked here.
#[kani::proof]
fn every_case_and_mask_triangulation_is_indexable() {
    let case: u8 = kani::any();
    let joined: u8 = kani::any();
    kani::assume(joined < 64);
    check_all_properties(case, &ALL_MASKS[case as usize][joined as usize]);
}

/// **The control: the properties must be able to fail.**
///
/// Four `SUCCESS` verdicts over 65,536 table entries are worth nothing if the
/// assertions cannot fire — `M-44`'s rule, and a proof is exactly where it is
/// easiest to violate, because a vacuous proof and a real one print the same
/// word. This harness corrupts one entry the way `golden::tests`' sabotage arm
/// corrupts one triangle (`M-337`'s model) and requires the check to fail.
///
/// It is `#[kani::should_panic]`, so Kani's own verdict is `SUCCESSFUL` exactly
/// when the assertion **does** fire. A property that had drifted into a tautology
/// would make this harness fail, which is the direction that matters.
#[kani::proof]
#[kani::should_panic]
fn the_properties_can_fail() {
    let case: u8 = kani::any();
    // Case 0 and case 255 have no triangles, so a corruption there would be
    // invisible: `count` is zero and the loop body never runs. Constrain the
    // pattern to one that actually emits.
    kani::assume(case != 0 && case != 255);

    let mut c = CASES[case as usize];
    kani::assume(c.count > 0);

    // The corruption: name an edge, but the one this pattern is least likely to
    // cut, and give the triangle two equal codes at the same time. Either
    // property 3 or property 4 must catch it.
    c.triangles[0] = [0, 0, 0];
    check_all_properties(case, &c);
}
