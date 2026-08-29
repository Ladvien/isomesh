//! **P-125 — the pinch predicate, shipped as a `validate` report.**
//!
//! Ticket: R-053. Pre-registered before this harness existed.
//!
//! ```bash
//! cargo bench --bench experiment_p125
//! ```
//!
//! Writes `docs/experiments/p-125.csv`.
//!
//! # What was missing
//!
//! `M-352` measured the `=`-corner repair on two real CT volumes and found it
//! **safe on one and topology-changing on the other**. Both lost every
//! degenerate triangle — 164 of 164 on `fuel` 64³, 58,097 of 58,097 on `bonsai`
//! 256³ — and `max_snap_distance` was **exactly 0** on both, so the repair moved
//! no geometry whatsoever and is purely combinatorial. On `fuel`, χ 19 → 19,
//! non-manifold edges 0 → 0, boundary edges 24 → 24. On `bonsai`, the identical
//! repair took χ 517 → 585, non-manifold edges **0 → 561** and boundary edges
//! 4,366 → 3,716, because **516 of its 17,201** collapse groups are *pinches*
//! and identifying them welded **520** separate pieces of a plant scan together.
//!
//! `M-352`'s own closing sentence is the thing that had never been built: *"The
//! precondition is cheap to test — one union-find over the baseline triangles —
//! and that is the shippable result here, not the label."* Until this row the
//! crate had no way for a caller to ask it. A caller could only be handed a
//! repair that is safe on somebody else's scan, and would find out otherwise as
//! a moved Euler characteristic and 561 non-manifold edges — passing every gate
//! except `euler_characteristic`.
//!
//! **This row therefore lands a source change, and the change was registered in
//! advance** — `crates/isomesh/src/validate/pinch.rs`, beside `validate::sealing`
//! — the `P-61` and `P-69` way, because a landing that was not registered in
//! advance is exactly what `V-45` cost.
//!
//! # The predicate, and what a pinch is
//!
//! Two vertices a collapse identifies either **already share a triangle**, in
//! which case that triangle is one of the degenerates and the merge flattens a
//! fold — no edge, boundary or component can move — or they **share no
//! triangle**, in which case they lie on pieces of the surface that meet only at
//! that point and identifying them welds those pieces. The transitive closure of
//! "shares a triangle" *within* a group gives its **sharing clusters**; one
//! cluster is a fold, two or more is a pinch joining `clusters − 1` pieces.
//!
//! # Six fixtures, one build, one run (`M-281`)
//!
//! | `fixture` | shape | `is_control` | what it is for |
//! |---|---|---|---|
//! | `fuel_iso32` | 64³ `uint8` CT, isovalue 32 | no | C1's zero: **0 of 50** |
//! | `fuel_half_offset` | the same volume at 32.5 | no | C1's *other* zero — a `u8` cannot equal 32.5, so there is nothing to collapse and the instrument must say 0 groups rather than invent them |
//! | `bonsai_iso32` | 256³ `uint8` CT, isovalue 32 | no | C1's non-zero: **516 of 17,201**, and the vacuity control for every zero above |
//! | `bonsai_half_offset` | the same volume at 32.5 | no | as `fuel_half_offset`, at 21× the vertex count |
//! | `constructed_pinch` | hand-built, 38 vertices, 14 faces | **yes** | the deliberately pinched fixture the registration demands, constructed rather than sampled the way `✗22`'s was: 3 two-cluster pinches, 1 three-cluster pinch, 2 folds and 1 untouched triangle, so `collapse_groups` **6**, `pinch_groups` **4** and `pinch_excess_components` **5** are known before the run. Also **C3's positive control** |
//! | `constructed_long_range` | hand-built, 10 vertices, 8 faces | **yes** | **C3's negative control.** Two vertices at one position that share no triangle — a pinch by the predicate — but joined by a six-triangle strip, so identifying them welds *nothing globally*. `components_welded` **0** against `pinch_excess_components` **1**, so `c3_holds` is **false by construction** |
//!
//! The two `constructed_` rows carry `is_control = true` and are excluded from
//! the run-level verdict, which is taken over the four CT rows. Every registered
//! column is present on every row; the two `constructed_pinch_*` columns are
//! properties of the *run* and are therefore repeated on all six, so the control
//! is visible beside every number it licenses.
//!
//! # SHARE, recomputed before the numbers
//!
//! **None, and it is registered as none.** The report is a validate-time query
//! and is on no extraction path: nothing in `crates/isomesh/src/extractor.rs`,
//! `marching_cubes/`, `dual_contouring/` or `chunk/` calls it, and it did not
//! exist until this ticket. So there is no extraction fraction for a ceiling to
//! be computed against, `✗51`'s share bar does not apply, and **no clause here
//! is or may become a speedup claim**.
//!
//! Which is also why `M-280`/`✗24` has nothing to bite on. Every clause is an
//! **integer equality over an enumerated population**:
//!
//! * **C1** — `collapse_groups` and `pinch_groups` against `M-352`'s own
//!   figures. Integers over a population the fixture fixes exactly; there is no
//!   tolerance to hide in, and the registration says so.
//! * **C2** — `distinct_censuses` over `order_permutations` = **128** face
//!   permutations, and whether the returned buffers were reserved exactly once.
//!   Both integers.
//! * **C3** — `components_welded` against `pinch_excess_components`. Integers.
//!
//! No clause is a cost or a ratio, so **nothing here is timed and nothing here is
//! counted by a hardware counter**, and this is the one Phase 25/26 harness that
//! needs no `perf_event_open` and therefore runs on any platform. `M-280` says a
//! cost clause must read retired instructions; the corollary, which `P-112`
//! learned by deleting the column, is that a figure no clause reads should not be
//! recorded at all.
//!
//! # The registration's C2 overstates two things, and the harness says which
//!
//! Neither is amended — `crates/isomesh/src/experiment.rs:27-31` forbids that —
//! and neither is in C2's falsifier, which gates only on *face-order dependence*
//! and on *per-group allocation*.
//!
//! 1. **`O(V + F)` is `O(V log V + V·k + F)`.** The `log V` is a sort and the
//!    `k` is the 27-cell probe, and both are `validate`'s own conventions rather
//!    than accidents: `validate.rs:13-24` states why every structural pass in
//!    that module is *sort a flat `Vec`, then scan runs of equal keys*, and
//!    `validate.rs:954-961` states why the coincidence neighbourhood has to be
//!    the one `weld` probes. A sort-free grouping would be a second convention
//!    and a narrower neighbourhood would describe a mesh no weld produced.
//!
//!    This is settled by **reading the source**, not by a column, and
//!    deliberately so. Two fixtures whose vertex counts are 3,200 and 533,221
//!    differ by 1.63× in `log₂ V` against a mixed expression with a large
//!    per-fixture constant, so an instruction count over them cannot separate
//!    `V log V` from `V` at all. A number that looks like evidence and is not is
//!    worse than its absence.
//! 2. **"allocates once" is seven buffers, allocated once each and none per
//!    group.** `buffers_reserved_exactly` is the measurable half: every one is
//!    `Vec::with_capacity`-reserved to a size known before its first push, so
//!    `len == capacity` on all three returned buffers is only true of a `Vec`
//!    that never reallocated — checked on `bonsai`'s 17,201 groups, where any
//!    per-group growth would have over-reserved by the usual doubling. A
//!    counting global allocator would settle it outright and **cannot exist
//!    here** — `Cargo.toml:42` sets `unsafe_code = "forbid"` workspace-wide,
//!    which `#[allow]` cannot override, and `impl GlobalAlloc` is an
//!    `unsafe impl`.
//!
//! # C3 has two readings and only one of them can hold
//!
//! The registration's C3 reads *"applying the repair and counting the components
//! it welds reproduces the report's pinch count exactly"*. The report has two
//! pinch figures and they differ: **516** pinch groups and **520** pieces joined
//! on `bonsai`, so four of its groups span three or more sharing clusters.
//!
//! A group of three clusters welds **two** components, so an identity against the
//! group *count* is arithmetically impossible on any fixture where a group spans
//! three — and `M-352`'s own 520 ≠ 516 says four `bonsai` groups do. `c3_holds` is
//! therefore scored on the reading that can hold, `components_welded ==
//! pinch_excess_components`, and the other reading is recorded beside it as
//! `components_welded_equals_pinch_groups` so the entry can quote either. That is
//! a clause decided by arithmetic that was available before this harness ran, and
//! it is reported rather than resolved silently.
//!
//! **What `components_welded` measures.** A repair does two things — it
//! identifies vertices and it drops the faces that identification degenerates —
//! and only the first can weld. So `components_welded` is the number of
//! *baseline* connected components the identification joins, counted as the
//! successful unions when every group's members are unioned into a disjoint-set
//! built from the baseline faces. `components_before` and `components_after` are
//! the measured endpoints over the actual repaired index buffer, and
//! `components_drift_from_dropped_faces` closes the accounting between them:
//! dropping a fold's only face **removes** a component, and dropping a face that
//! was a component's only bridge **splits** one, neither of which is a weld.
//!
//! Components here are **vertex**-connected components over triangle sides,
//! which is the transitive closure of the very relation the predicate tests
//! locally. Deliberately not `MeshReport::components`, which counts
//! **edge**-connected components of the *face* set (`validate.rs:184`): two
//! triangles meeting at a pinch point share a vertex and no edge, so that count
//! is structurally unable to see the thing C3 is about.
//!
//! # The vacuity controls, all four asserted rather than recorded
//!
//! `M-44`: a zero that could not have been non-zero is not a measurement, and
//! `fuel`'s **0 of 50** is a zero. A run that cannot fire aborts instead of
//! recording a pass.
//!
//! 1. **The same instrument returns non-zero in the same run.**
//!    `bonsai_iso32`'s `pinch_groups` is asserted `> 0`. C1 already requires it;
//!    it is asserted anyway, because a zero everywhere would otherwise read as
//!    four passes.
//! 2. **The deliberately pinched fixture reports its construction.**
//!    `constructed_pinch_reported == constructed_pinch_expected`, both **4**,
//!    and `pinch_excess_components` asserted **5** on that row — the fixture is
//!    built so the two differ, so a harness that conflated them could not pass.
//! 3. **The predicate can answer *yes, they share*.** For every face of every
//!    fixture, every edge whose two ends are distinct members of one collapse
//!    group must carry the **same** cluster label. `share_control_all_true` is
//!    asserted on all six rows — a union-find that could not answer *shared*
//!    would call every group a pinch and `M-352`'s 516 would be an artefact.
//!
//!    The population is asserted non-empty in two places, and the shape of that
//!    assertion is itself a finding this harness had to get right. A group whose
//!    members form one cluster got there by a union, and a union comes from a
//!    face edge — so **a fold implies a shared pair exactly**, and
//!    `share_control_pairs` is asserted non-zero on every row where
//!    `collapse_groups > pinch_groups`. Where every group is a pinch there is
//!    legitimately nothing to ask, and `constructed_long_range` is built to be
//!    precisely that row; asserting non-emptiness there would void a correct
//!    run. The run-level tally is asserted non-zero instead, so the control
//!    cannot be vacuous everywhere at once.
//! 4. **C3 is proven able to read both ways, in the same run.**
//!    `constructed_pinch` is built so C3 holds and its `c3_holds` is asserted
//!    **true**: a C3 that cannot hold on a fixture constructed to satisfy it is
//!    an instrument fault, not a finding. `constructed_long_range` is built so
//!    C3 fails and its `c3_holds` is asserted **false**: a control that reads
//!    true is `P-63`-C3's failure mode arriving again, a zero over a population
//!    that cannot host the phenomenon.
//!
//! And one fixture-integrity cross-check, asserted on every row:
//! `repair_matches_report` — the faces this harness actually drops when it
//! applies the repair must equal the shipped report's `folding_faces`. Without
//! it, C3 could be measuring a repair the census never described, which is
//! `experiment_p101.rs`'s `replica_bit_identical` rule in a different costume.
//!
//! # What is bench-local, and what is the crate's
//!
//! The census, the groups and the cluster labels are the crate's own
//! `isomesh::validate::pinch_features` — this harness does not transcribe the
//! predicate, which is the whole point of a row that lands a source change. The
//! repair, the component counts and the permutation sweep are bench-local, and
//! **every mechanism copied out of a private module is copied here with the
//! source line it came from on the row that uses it**, per
//! `experiment_p117.rs:53-56`: `crates/isomesh/src/` is read-only for this row
//! apart from the registered `validate` addition, and a copy whose line number is
//! a comment is auditable in a way a `pub` would not be.
//!
//! The repair applies to **the shipped groups**, never to a re-derived
//! partition. Two groupings that disagreed by one member would make C3 compare a
//! census of one collapse against the components of another.

// Exact comparisons on purpose: `max_snap_distance` is a claim about bits, and
// `M-352` measured it as exactly zero.
#![allow(
    clippy::float_cmp,
    reason = "the snap distance is asserted against exactly zero, which is what M-352 measured"
)]

mod common;

use std::path::{Path, PathBuf};

use isomesh::construct::SampledField;
use isomesh::marching_cubes::MarchingCubes;
use isomesh::validate::{PinchGroups, PinchReport, ValidateConfig, pinch_features};
use isomesh::{MeshBuffer, RuntimeShape3};

/// One world unit per voxel: both datasets declare `spacing 1x1x1`, and it is
/// what makes `weld_epsilon` `1e-4` absolute — `ValidateConfig::from_cell_size`
/// derives it as `cell_size · WELD_EPSILON_REL` (`validate.rs:106`, `:77`).
const CELL: f64 = 1.0;

/// Sample `[0, 0, 0]` sits at the world origin, so a corner's world position is
/// its integer grid coordinate exactly.
const ORIGIN: [f64; 3] = [0.0; 3];

/// Face permutations per fixture, including the identity.
///
/// The number the registration names. `✗26` was a face-iteration-order leak
/// found after a landing; this asks the question before one.
const ORDER_PERMUTATIONS: u64 = 128;

/// A `uint8` CT volume and the isovalue to contour it at.
struct Volume {
    file: &'static str,
    short: &'static str,
    /// Samples per axis; the files are cubes.
    n: u32,
    iso: f64,
    /// The label the row carries, so `32.5` never has to be re-formatted.
    iso_label: &'static str,
    /// `M-352`'s own figures for this row: collapse groups, then pinches.
    expect: (u64, u64),
}

/// The four CT rows.
///
/// The half-offset rows are not decoration. `M-317`'s guidance is to contour at
/// a half-integer precisely because integer data cannot reach it, so
/// `equal_corners` is 0 there (`p-53.csv`) and the census must report **0
/// groups** rather than invent them out of the lattice — which is the control
/// that says `17,201` is a property of the data and not of the tolerance.
const VOLUMES: [Volume; 4] = [
    Volume {
        file: "fuel_64x64x64_uint8.raw",
        short: "fuel_iso32",
        n: 64,
        iso: 32.0,
        iso_label: "32",
        expect: (50, 0),
    },
    Volume {
        file: "fuel_64x64x64_uint8.raw",
        short: "fuel_half_offset",
        n: 64,
        iso: 32.5,
        iso_label: "32.5",
        expect: (0, 0),
    },
    Volume {
        file: "bonsai_256x256x256_uint8.raw",
        short: "bonsai_iso32",
        n: 256,
        iso: 32.0,
        iso_label: "32",
        expect: (17_201, 516),
    },
    Volume {
        file: "bonsai_256x256x256_uint8.raw",
        short: "bonsai_half_offset",
        n: 256,
        iso: 32.5,
        iso_label: "32.5",
        expect: (0, 0),
    },
];

fn dir() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../docs/measurements/volumes")
}

/// Read a `uint8` volume, raw. The length is checked against the dimensions
/// rather than trusted from the filename, as `benches/volumes.rs:79-92` does.
fn read_u8(path: &Path, n: u32) -> Result<Vec<u8>, String> {
    let bytes = std::fs::read(path).map_err(|e| format!("{}: {e}", path.display()))?;
    let want = (n as usize).pow(3);
    if bytes.len() != want {
        return Err(format!(
            "{}: {} bytes, expected {want} for {n}³",
            path.display(),
            bytes.len()
        ));
    }
    Ok(bytes)
}

/// A mesh to census, and what its census must say.
struct Fixture {
    name: String,
    isovalue: &'static str,
    positions: Vec<[f64; 3]>,
    indices: Vec<u32>,
    /// Collapse groups, pinch groups, pieces joined — from `M-352` for the CT
    /// rows and from the construction for the built ones.
    expect: (u64, u64, u64),
    /// Excluded from the run-level verdict: a fixture built to make a clause
    /// read a particular way is evidence about the instrument, not about the
    /// data.
    is_control: bool,
    /// `Some(false)` on `constructed_long_range`, whose C3 is false by
    /// construction; `Some(true)` on `constructed_pinch`. `None` where the
    /// answer is the measurement.
    c3_must_be: Option<bool>,
}

/// The deliberately pinched fixture: **constructed, not sampled**.
///
/// `✗22` is why. Six anchors ten units apart, so no two anchors are within a
/// thousand `weld_epsilon`s of each other and every group is exactly the
/// coincidence it was built to be:
///
/// * three anchors carry **two** disjoint triangles sharing one position — a
///   two-cluster pinch each, joining one piece each;
/// * one anchor carries **three** disjoint triangles sharing one position — a
///   three-cluster pinch, joining **two** pieces, which is what makes
///   `pinch_groups` (4) and `pinch_excess_components` (5) differ on this row
///   the way they differ on `bonsai` (516 and 520);
/// * two anchors carry a triangle with **two of its own corners at one point**
///   plus a neighbour sharing that triangle's surviving edge — a fold, so not a
///   pinch, and the neighbour is there so that dropping the degenerate face does
///   not also delete the component and contaminate C3;
/// * one anchor carries a plain triangle with no coincidence at all, so the
///   fixture is not made entirely of degeneracies.
///
/// Known by construction: 38 vertices, 14 faces, 6 collapse groups, 4 pinches,
/// 5 pieces joined, 2 folding faces, 12 components before and 7 after.
fn constructed_pinch() -> Fixture {
    let mut positions: Vec<[f64; 3]> = Vec::new();
    let mut indices: Vec<u32> = Vec::new();
    let push = |tri: [[f64; 3]; 3], positions: &mut Vec<[f64; 3]>, indices: &mut Vec<u32>| {
        let base = u32::try_from(positions.len()).expect("a small fixture");
        positions.extend(tri);
        indices.extend([base, base + 1, base + 2]);
    };

    // Three two-cluster pinches.
    for k in 0..3 {
        let x = f64::from(k) * 10.0;
        push(
            [[x, 0.0, 0.0], [x + 1.0, 0.0, 0.0], [x, 1.0, 0.0]],
            &mut positions,
            &mut indices,
        );
        push(
            [[x, 0.0, 0.0], [x - 1.0, 0.0, 0.0], [x, -1.0, 0.0]],
            &mut positions,
            &mut indices,
        );
    }
    // One three-cluster pinch.
    {
        let x = 100.0;
        push(
            [[x, 0.0, 0.0], [x + 1.0, 0.0, 0.0], [x, 1.0, 0.0]],
            &mut positions,
            &mut indices,
        );
        push(
            [[x, 0.0, 0.0], [x - 1.0, 0.0, 0.0], [x, -1.0, 0.0]],
            &mut positions,
            &mut indices,
        );
        push(
            [[x, 0.0, 0.0], [x, 0.0, 1.0], [x + 1.0, 0.0, 1.0]],
            &mut positions,
            &mut indices,
        );
    }
    // Two folds, each with a neighbour across its surviving edge.
    for k in 0..2 {
        let x = 200.0 + f64::from(k) * 10.0;
        let base = u32::try_from(positions.len()).expect("a small fixture");
        positions.extend([
            [x, 0.0, 0.0],
            [x, 0.0, 0.0],
            [x, 1.0, 0.0],
            [x + 1.0, 1.0, 0.0],
        ]);
        indices.extend([base, base + 1, base + 2, base, base + 2, base + 3]);
    }
    // One triangle that coincides with nothing.
    push(
        [[300.0, 0.0, 0.0], [301.0, 0.0, 0.0], [300.0, 1.0, 0.0]],
        &mut positions,
        &mut indices,
    );

    Fixture {
        name: String::from("constructed_pinch"),
        isovalue: "none",
        positions,
        indices,
        expect: (6, 4, 5),
        is_control: true,
        c3_must_be: Some(true),
    }
}

/// **C3's negative control: a pinch that welds nothing.**
///
/// A six-triangle strip, plus two vertices at one position hung off its two
/// ends. They share no triangle, so the predicate calls the group a pinch and
/// `pinch_excess_components` reads **1** — but the strip already connects them,
/// so identifying them merges no component and `components_welded` reads **0**.
///
/// This is the failure mode C3 exists to detect, built so that it must occur:
/// the predicate is local and the question is global. A run in which this row's
/// `c3_holds` came out **true** would mean the harness cannot tell the two
/// apart, and it aborts rather than reporting four passes.
fn constructed_long_range() -> Fixture {
    let positions = vec![
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [2.0, 0.0, 0.0],
        [3.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [1.0, 1.0, 0.0],
        [2.0, 1.0, 0.0],
        [3.0, 1.0, 0.0],
        [10.0, 10.0, 0.0],
        [10.0, 10.0, 0.0],
    ];
    let indices = vec![
        0, 1, 4, 1, 5, 4, 1, 2, 5, 2, 6, 5, 2, 3, 6, 3, 7, 6, // the strip
        0, 4, 8, // vertex 8 on its left end
        3, 7, 9, // vertex 9 on its right end
    ];
    Fixture {
        name: String::from("constructed_long_range"),
        isovalue: "none",
        positions,
        indices,
        expect: (1, 1, 1),
        is_control: true,
        c3_must_be: Some(false),
    }
}

/// Marching Cubes over one `uint8` volume at one isovalue.
///
/// `iso − value`, so a dense voxel is negative and the crate's sign convention
/// holds unchanged — `benches/volumes.rs:31-36` and `experiment_p53.rs:736-738`.
fn ct_fixture(v: &Volume) -> Fixture {
    let raw = match read_u8(&dir().join(v.file), v.n) {
        Ok(raw) => raw,
        Err(e) => {
            // Not a skip. C1 is denominated in these two volumes, so a missing
            // file makes the row unmeasurable rather than smaller.
            println!("::error:: {e}");
            std::process::exit(1);
        }
    };
    let shape = match RuntimeShape3::new([v.n; 3]) {
        Ok(shape) => shape,
        Err(e) => {
            println!("::error:: {}: {e}", v.file);
            std::process::exit(1);
        }
    };
    let values: Vec<f64> = raw.iter().map(|b| v.iso - f64::from(*b)).collect();
    let field = match SampledField::new(&values, &shape, ORIGIN, CELL) {
        Ok(field) => field,
        Err(e) => {
            println!("::error:: {}: {e}", v.file);
            std::process::exit(1);
        }
    };
    let mut mesh = MeshBuffer::<f64>::new();
    if let Err(e) = MarchingCubes::<f64>::new().extract(&field, &shape, ORIGIN, CELL, &mut mesh) {
        println!("::error:: {} at iso {}: {e}", v.short, v.iso_label);
        std::process::exit(1);
    }
    Fixture {
        name: String::from(v.short),
        isovalue: v.iso_label,
        positions: mesh.positions,
        indices: mesh.indices,
        expect: (v.expect.0, v.expect.1, 0),
        is_control: false,
        c3_must_be: None,
    }
}

// ── bench-local mechanisms, each with the line it was copied from ────────────

/// Root of `x`, halving the path on the way.
///
/// `validate::Dsu` (`validate.rs:541`) is private and so is
/// `validate::pinch`'s own copy. Unioned to the **lower** root, as
/// `experiment_p53.rs:454-484` did, so the representative is the set's minimum
/// and the component count below is a pure function of the mesh rather than of
/// the order the unions arrived in.
fn find(parent: &mut [u32], mut x: u32) -> u32 {
    while parent[x as usize] != x {
        let grand = parent[parent[x as usize] as usize];
        parent[x as usize] = grand;
        x = grand;
    }
    x
}

/// Join two sets, keeping the lower index as the root. `true` if they were
/// separate — which is what makes the merge *count*.
fn union_to_lower(parent: &mut [u32], a: u32, b: u32) -> bool {
    let (ra, rb) = (find(parent, a), find(parent, b));
    if ra == rb {
        return false;
    }
    let (lo, hi) = if ra < rb { (ra, rb) } else { (rb, ra) };
    parent[hi as usize] = lo;
    true
}

/// One disjoint-set over vertex indices, with every triangle side unioned.
///
/// **Vertex** connectivity, not the face-edge connectivity
/// `MeshReport::components` reports (`validate.rs:184`): two triangles meeting
/// at a pinch share a vertex and no edge, so a face-edge count cannot see a
/// weld at a point.
fn vertex_dsu(n: usize, indices: &[u32]) -> Vec<u32> {
    let mut parent: Vec<u32> = (0..n as u32).collect();
    for tri in indices.as_chunks::<3>().0 {
        union_to_lower(&mut parent, tri[0], tri[1]);
        union_to_lower(&mut parent, tri[1], tri[2]);
    }
    parent
}

/// Connected components over the vertices at least one face names.
fn components(n: usize, indices: &[u32]) -> u64 {
    let mut parent = vertex_dsu(n, indices);
    let mut referenced = vec![false; n];
    for &v in indices {
        referenced[v as usize] = true;
    }
    let mut count = 0u64;
    for v in 0..n as u32 {
        if referenced[v as usize] && find(&mut parent, v) == v {
            count += 1;
        }
    }
    count
}

/// How many baseline components the identification joins.
///
/// The repair does two separable things and only this one can weld: it
/// identifies each group's members, and it drops the faces that identification
/// degenerates. Every successful union reduces the component count by exactly
/// one, so the tally is the weld count exactly — and it is order-independent as
/// a total even though its attribution to a particular group would not be.
fn welded_components(n: usize, indices: &[u32], groups: &PinchGroups) -> u64 {
    let mut parent = vertex_dsu(n, indices);
    let mut merges = 0u64;
    for g in 0..groups.group_count() {
        let members = groups.members(g);
        for &v in &members[1..] {
            if union_to_lower(&mut parent, members[0], v) {
                merges += 1;
            }
        }
    }
    merges
}

/// Apply the collapse to **the shipped groups**, returning the repaired index
/// buffer and the faces it dropped.
///
/// The representative is the group's lowest-indexed member, which is the
/// tie-break `weld` documents (`weld.rs:45-46`) and the one
/// `experiment_p53.rs:574-576` used.
fn repair(n: usize, indices: &[u32], groups: &PinchGroups) -> (Vec<u32>, u64) {
    let mut remap: Vec<u32> = (0..n as u32).collect();
    for g in 0..groups.group_count() {
        let members = groups.members(g);
        for &v in members {
            remap[v as usize] = members[0];
        }
    }
    let mut out: Vec<u32> = Vec::with_capacity(indices.len());
    let mut dropped = 0u64;
    for tri in indices.as_chunks::<3>().0 {
        let a = remap[tri[0] as usize];
        let b = remap[tri[1] as usize];
        let c = remap[tri[2] as usize];
        if a == b || b == c || c == a {
            dropped += 1;
            continue;
        }
        out.extend_from_slice(&[a, b, c]);
    }
    (out, dropped)
}

/// The largest distance any vertex moves, as a Chebyshev norm over the axes.
///
/// `M-352` measured exactly `0` on both volumes because the vertices are already
/// at the corner. Measured rather than assumed: if it is not zero, the collapse
/// is moving geometry and the registration's "pure connectivity decision" does
/// not hold for that row.
fn max_snap(positions: &[[f64; 3]], groups: &PinchGroups) -> f64 {
    let mut worst = 0.0f64;
    for g in 0..groups.group_count() {
        let members = groups.members(g);
        let to = positions[members[0] as usize];
        for &v in &members[1..] {
            let from = positions[v as usize];
            for a in 0..3 {
                worst = worst.max((to[a] - from[a]).abs());
            }
        }
    }
    worst
}

/// Ask the shipped predicate of every pair that demonstrably shares a triangle.
///
/// Returns `(pairs asked, every answer was "shared")`. A pair is asked only when
/// both ends are distinct members of the **same** collapse group, which is
/// exactly the population the pinch decision is made over.
fn share_control(n: usize, indices: &[u32], groups: &PinchGroups) -> (u64, bool) {
    let mut slot = vec![u32::MAX; n];
    for (s, &v) in groups.vertices.iter().enumerate() {
        slot[v as usize] = u32::try_from(s).expect("one slot per vertex");
    }
    let mut group_of_slot = vec![0u32; groups.vertices.len()];
    for g in 0..groups.group_count() {
        let (lo, hi) = (groups.starts[g] as usize, groups.starts[g + 1] as usize);
        for entry in &mut group_of_slot[lo..hi] {
            *entry = u32::try_from(g).expect("a group index");
        }
    }

    let mut pairs = 0u64;
    let mut all_shared = true;
    for tri in indices.as_chunks::<3>().0 {
        for (a, b) in [(0usize, 1usize), (1, 2), (2, 0)] {
            let (u, w) = (tri[a], tri[b]);
            if u == w {
                continue;
            }
            let (su, sw) = (slot[u as usize], slot[w as usize]);
            if su == u32::MAX || sw == u32::MAX {
                continue;
            }
            if group_of_slot[su as usize] != group_of_slot[sw as usize] {
                continue;
            }
            pairs += 1;
            if groups.clusters[su as usize] != groups.clusters[sw as usize] {
                all_shared = false;
            }
        }
    }
    (pairs, all_shared)
}

/// SplitMix64, so the permutations are reproducible from the seed alone.
///
/// Steele, Lea & Flood, *Fast splittable pseudorandom number generators*, OOPSLA
/// 2014 (`10.1145/2660193.2660195`). A shuffled face list has to be the same
/// shuffled face list on the next run, or `distinct_censuses` would be measuring
/// the generator.
fn splitmix64(state: &mut u64) -> u64 {
    *state = state.wrapping_add(0x9E37_79B9_7F4A_7C15);
    let mut z = *state;
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    z ^ (z >> 31)
}

/// How many *distinct* censuses 128 permutations of the face list produce.
///
/// The whole `(PinchReport, PinchGroups)` pair is compared, not only the counts:
/// the cluster labels are a public artefact, so a face-order leak in them is as
/// much a defect as one in the pinch count. Permutation 0 is the identity, so a
/// reading of 1 means every one of the other 127 reproduced it exactly.
fn distinct_censuses(
    positions: &[[f64; 3]],
    indices: &[u32],
    cfg: &ValidateConfig,
    baseline: &(PinchReport, PinchGroups),
) -> u64 {
    let faces = indices.len() / 3;
    let mut distinct: Vec<(PinchReport, PinchGroups)> = vec![baseline.clone()];
    let mut order: Vec<u32> = Vec::with_capacity(faces);
    let mut permuted: Vec<u32> = Vec::with_capacity(indices.len());
    for p in 1..ORDER_PERMUTATIONS {
        let mut state = 0x0000_0000_0000_007D ^ p;
        order.clear();
        order.extend(0..u32::try_from(faces).expect("a face count"));
        for i in (1..faces).rev() {
            let j = (splitmix64(&mut state) % (i as u64 + 1)) as usize;
            order.swap(i, j);
        }
        permuted.clear();
        for &f in &order {
            let at = f as usize * 3;
            permuted.extend_from_slice(&indices[at..at + 3]);
        }
        let census = pinch_features(positions, &permuted, cfg);
        if !distinct.contains(&census) {
            distinct.push(census);
        }
    }
    u64::try_from(distinct.len()).expect("a small count")
}

/// One fixture, measured.
struct Measured {
    report: PinchReport,
    max_snap: f64,
    distinct: u64,
    components_before: u64,
    components_after: u64,
    components_welded: u64,
    drift: i64,
    triangles_after: u64,
    share_pairs: u64,
    share_all_shared: bool,
    buffers_exact: bool,
    repair_matches: bool,
    c1: bool,
    c2: bool,
    c3: bool,
    c3_literal: bool,
}

fn measure(fixture: &Fixture, cfg: &ValidateConfig) -> Measured {
    let n = fixture.positions.len();
    let positions = &fixture.positions;
    let indices = &fixture.indices;

    let baseline = pinch_features(positions, indices, cfg);
    let (report, groups) = (baseline.0, baseline.1.clone());
    // `Vec::with_capacity(k)` followed by exactly `k` pushes is the only way a
    // `Vec` ends with `len == capacity` on a non-empty buffer; anything that
    // grew would have over-reserved.
    let buffers_exact = groups.vertices.len() == groups.vertices.capacity()
        && groups.clusters.len() == groups.clusters.capacity()
        && groups.starts.len() == groups.starts.capacity();

    let (repaired, dropped) = repair(n, indices, &groups);
    let components_before = components(n, indices);
    let components_after = components(n, &repaired);
    let components_welded = welded_components(n, indices, &groups);
    let (share_pairs, share_all_shared) = share_control(n, indices, &groups);
    let distinct = distinct_censuses(positions, indices, cfg, &baseline);

    let c1 = report.collapse_groups == fixture.expect.0 && report.pinch_groups == fixture.expect.1;
    let c2 = distinct == 1 && buffers_exact;
    let c3 = components_welded == report.pieces_joined;

    Measured {
        report,
        max_snap: max_snap(positions, &groups),
        distinct,
        components_before,
        components_after,
        components_welded,
        drift: components_after as i64 - (components_before as i64 - components_welded as i64),
        triangles_after: (repaired.len() / 3) as u64,
        share_pairs,
        share_all_shared,
        buffers_exact,
        repair_matches: dropped == report.folding_faces,
        c1,
        c2,
        c3,
        c3_literal: components_welded == report.pinch_groups,
    }
}

fn main() {
    if !std::env::args().any(|arg| arg == "--bench") {
        return;
    }
    let prereg = isomesh::experiment!("P-125");

    common::experiment::run(prereg, |run| {
        let cfg = match ValidateConfig::from_cell_size(CELL) {
            Ok(cfg) => cfg,
            Err(e) => {
                println!("::error:: {e}");
                std::process::exit(1);
            }
        };

        let mut fixtures: Vec<Fixture> = VOLUMES.iter().map(ct_fixture).collect();
        fixtures.push(constructed_pinch());
        fixtures.push(constructed_long_range());

        let measured: Vec<Measured> = fixtures.iter().map(|f| measure(f, &cfg)).collect();

        // ── the vacuity controls, before any verdict is reported ────────────
        let at = |name: &str| -> &Measured {
            let i = fixtures
                .iter()
                .position(|f| f.name == name)
                .expect("a fixture this harness built");
            &measured[i]
        };

        // 1. The same instrument, non-zero, in the same run.
        let bonsai = at("bonsai_iso32");
        assert!(
            bonsai.report.pinch_groups > 0,
            "VOID: the predicate reports no pinch on bonsai, so fuel's zero is a \
             zero that could not have been non-zero (M-44) and every clause here \
             is unmeasured:\n{}",
            bonsai.report
        );

        // 2. The deliberately pinched fixture reports its construction.
        let built = at("constructed_pinch");
        assert_eq!(
            built.report.pinch_groups, 4,
            "VOID: the constructed fixture has four pinches by construction and \
             the instrument found {}:\n{}",
            built.report.pinch_groups, built.report
        );
        assert_eq!(
            built.report.pieces_joined, 5,
            "VOID: the constructed fixture joins five pieces by construction — \
             three two-cluster pinches and one three-cluster — and the \
             instrument found {}",
            built.report.pieces_joined
        );

        // 3. The predicate can answer "yes, they share", and was asked.
        for (f, m) in fixtures.iter().zip(&measured) {
            assert!(
                m.share_all_shared,
                "VOID: {}: two vertices that demonstrably share a triangle came \
                 out in different sharing clusters, so every pinch count in this \
                 run is an artefact of a union-find that cannot answer",
                f.name
            );
            // A group whose members are one cluster got there by a union, and a
            // union comes from a face edge — so a fold implies a shared pair
            // exactly. Where every group is a pinch there is legitimately
            // nothing to ask, `constructed_long_range` being the case built to
            // be exactly that, and the run-level tally below is what keeps the
            // control from being vacuous everywhere at once.
            assert!(
                m.share_pairs > 0 || m.report.collapse_groups == m.report.pinch_groups,
                "VOID: {}: {} collapse groups of which {} are pinches, so at \
                 least one is a fold, and yet not one pair of members shares a \
                 triangle — the folds came from nowhere",
                f.name,
                m.report.collapse_groups,
                m.report.pinch_groups
            );
            assert!(
                m.repair_matches,
                "VOID: {}: applying the collapse dropped a different number of \
                 faces than the census predicted, so C3 would be comparing a \
                 repair the census never described",
                f.name
            );
        }
        let asked: u64 = measured.iter().map(|m| m.share_pairs).sum();
        assert!(
            asked > 0,
            "VOID: not one pair of vertices that demonstrably share a triangle \
             was asked of the predicate anywhere in this run, so nothing shows it \
             able to answer `shared` and every pinch count is a default"
        );

        // 4. C3 is shown able to read both ways, in this run.
        for (f, m) in fixtures.iter().zip(&measured) {
            if let Some(want) = f.c3_must_be {
                assert_eq!(
                    m.c3, want,
                    "VOID: {}: C3 must read {want} on this fixture by \
                     construction — components_welded {} against \
                     pinch_excess_components {}",
                    f.name, m.components_welded, m.report.pieces_joined
                );
            }
        }

        // ── the rows ───────────────────────────────────────────────────────
        for (f, m) in fixtures.iter().zip(&measured) {
            let r = &m.report;
            println!(
                "{:>22} iso {:<5} v {:>7} t {:>8}  groups {:>6} (expected {:>6})  \
                 pinches {:>5} (expected {:>5})  pieces {:>5}",
                f.name,
                f.isovalue,
                r.vertices,
                r.triangles,
                r.collapse_groups,
                f.expect.0,
                r.pinch_groups,
                f.expect.1,
                r.pieces_joined
            );
            println!(
                "                       components {} -> {} welded {} drift {}  \
                 folds {} drop {} faces  snap {:e}  share pairs {}  censuses {} of {}",
                m.components_before,
                m.components_after,
                m.components_welded,
                m.drift,
                r.collapse_groups - r.pinch_groups,
                r.folding_faces,
                m.max_snap,
                m.share_pairs,
                m.distinct,
                ORDER_PERMUTATIONS
            );
            println!(
                "                       C1 {} C2 {} C3 {} (literal reading {})  control {}",
                m.c1, m.c2, m.c3, m.c3_literal, f.is_control
            );

            run.record(&[
                ("fixture", f.name.clone()),
                ("vertices", r.vertices.to_string()),
                ("triangles", r.triangles.to_string()),
                ("collapse_groups", r.collapse_groups.to_string()),
                ("pinch_groups", r.pinch_groups.to_string()),
                ("components_welded", m.components_welded.to_string()),
                ("components_before", m.components_before.to_string()),
                ("components_after", m.components_after.to_string()),
                ("max_snap_distance", format!("{:e}", m.max_snap)),
                ("order_permutations", ORDER_PERMUTATIONS.to_string()),
                ("distinct_censuses", m.distinct.to_string()),
                ("constructed_pinch_expected", String::from("4")),
                (
                    "constructed_pinch_reported",
                    built.report.pinch_groups.to_string(),
                ),
                ("share_control_all_true", m.share_all_shared.to_string()),
                ("c1_holds", m.c1.to_string()),
                ("c2_holds", m.c2.to_string()),
                ("c3_holds", m.c3.to_string()),
                // ── extras (M-273) ──────────────────────────────────────────
                ("isovalue", String::from(f.isovalue)),
                ("is_control", f.is_control.to_string()),
                ("expected_collapse_groups", f.expect.0.to_string()),
                ("expected_pinch_groups", f.expect.1.to_string()),
                ("pinch_excess_components", r.pieces_joined.to_string()),
                ("collapsing_vertices", r.collapsing_vertices.to_string()),
                (
                    "vertices_removed",
                    (r.collapsing_vertices - r.collapse_groups).to_string(),
                ),
                ("folding_faces", r.folding_faces.to_string()),
                ("triangles_after", m.triangles_after.to_string()),
                ("sharing_edges", r.sharing_edges.to_string()),
                (
                    "groups_moving_geometry",
                    r.groups_moving_geometry.to_string(),
                ),
                ("faces_skipped", r.faces_skipped.to_string()),
                ("share_control_pairs", m.share_pairs.to_string()),
                ("components_drift_from_dropped_faces", m.drift.to_string()),
                (
                    "components_welded_equals_pinch_groups",
                    m.c3_literal.to_string(),
                ),
                ("repair_matches_report", m.repair_matches.to_string()),
                ("buffers_reserved_exactly", m.buffers_exact.to_string()),
            ]);
        }

        // ── the run-level verdict, over the four measurement rows ───────────
        let verdict = |clause: fn(&Measured) -> bool| -> (usize, usize) {
            let rows: Vec<&Measured> = fixtures
                .iter()
                .zip(&measured)
                .filter(|(f, _)| !f.is_control)
                .map(|(_, m)| m)
                .collect();
            (rows.iter().filter(|m| clause(m)).count(), rows.len())
        };
        let (c1_held, rows) = verdict(|m| m.c1);
        let (c2_held, _) = verdict(|m| m.c2);
        let (c3_held, _) = verdict(|m| m.c3);
        println!();
        println!(
            "C1 held on {c1_held} of {rows} measurement rows, C2 on {c2_held}, C3 on {c3_held}; \
             the two constructed rows are controls and are excluded"
        );
        println!(
            "C3's literal reading (components_welded == pinch_groups) held on {} of {rows} — \
             it is arithmetically unreachable wherever a group spans three sharing clusters",
            verdict(|m| m.c3_literal).0
        );
    });
}
