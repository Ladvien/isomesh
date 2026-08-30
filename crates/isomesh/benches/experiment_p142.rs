//! **P-142 — the gyroid's Euler characteristic is `-8` per cubic cell, and the hole `CLAUDE.md` records is closable.**
//!
//! Ticket: R-142. Pre-registered before this harness existed.
//!
//! ```bash
//! cargo bench --bench experiment_p142
//! ```
//!
//! Writes `docs/experiments/p-142.csv`.
//!
//! # What was missing
//!
//! The crate has written down twice, as a standing rule, that the gyroid has no
//! assertable Euler characteristic:
//!
//! * `CLAUDE.md:199-203` — *"a blanket `V - E + F == 2` is unsatisfiable on two of
//!   the seven test fields … `gyroid` is triply periodic, so any finite sampling
//!   box cuts it and the result **has boundary** … Neither has an Euler
//!   characteristic derivable a priori, and inventing one violates rule 5."*
//! * `crates/isomesh/src/fields/mod.rs:23-28` — the same sentence in the module
//!   docs, and `capped_gyroid`'s own doc at `fields/mod.rs:1044-1052` gives the
//!   consequence: the shipped `gyroid` entry is `Gyroid ∩ Sphere{r:6}` over
//!   `[-7,7]³` **because** an uncapped periodic surface has boundary, and even
//!   capped *"its genus still is not known in closed form, which is why
//!   `expected_euler` returns `None`"* (`fields/mod.rs:1078`).
//!
//! Both statements are true of the box the crate samples and false of the
//! surface. The missing step is not a genus formula, it is a **periodic-conforming
//! box plus a seam identification**: a triply periodic minimal surface has no
//! boundary on the 3-torus, and the boundary the rule describes is an artefact of
//! cutting the torus open. `M-4` read `gyroid`'s non-manifold edges as a
//! "high-genus/open-field effect" and `M-15` corrected the *manifoldness* half of
//! that (`validate.rs:1111-1113`); the χ half has stood uncorrected since.
//!
//! So this row does not measure a new algorithm. It measures whether the number
//! the project declared underivable is `-8N³`, and it is the first row in the
//! repo to extract a reference-class field over a box that is an integer number
//! of its own periods.
//!
//! # The prediction is arithmetic, and the factor of two is the whole content
//!
//! Every triply periodic minimal surface has genus 3 in its **own primitive
//! translational cell**, so `χ = 2 - 2·3 = -4` there. What differs between the
//! three classical surfaces is how many primitive cells fit in the *conventional
//! cubic cell*, and for the gyroid that is **two**, because the body-centring
//! shift `(π,π,π)` leaves `F_G` **invariant** rather than negating it — two sign
//! flips per term. Hence `-8` per cubic cell and `-8N³` over `N³` of them, and
//! `genus = (2 - χ)/2 = 1 + 4N³`.
//!
//! That derivation is `common::tpms`'s, not this bench's, and its first vacuity
//! control below checks the invariance numerically rather than believing the
//! table — because if `(π,π,π)` merely negated `F_G` the prediction would be `-4`
//! and a harness that assumed `-8` would be grading its own arithmetic.
//!
//! # Why `voxels_per_period` is **odd**, and why 33/65/97 rather than 32/64/96
//!
//! The registration names *"32, 64, 96 and 128 voxels per period"*. Every one of
//! those is a multiple of 8, and `common::tpms`'s author measured 168
//! configurations and found that a multiple of 8 puts samples exactly on the
//! `π/4` lattice, where `F_D`'s four terms cancel to **exactly `0.0`** — `M-48`'s
//! degenerate crossing — so the crossing parameter is 0 or 1, one cell places
//! coincident vertices, and the weld turns them into a pinch. That family is
//! `-12` instead of `-16` at 32 and 56, `-9` at 64, `-7` at 96 and `+1` at 128 on
//! Schwarz D, and it bit Schwarz P once at `N = 3`.
//!
//! The gyroid was not observed to fail there — but the mechanism is a property of
//! the *grid*, not of the surface, and running the one field this row owns on the
//! one grid family known to manufacture exact zeros would be measuring the
//! extractor's degenerate-crossing handling under the name of a topological
//! oracle. So: **33, 65 and 97**, the odd neighbours of the registration's own
//! numbers, on which the module measured every surface exact. `nonmanifold_edges`
//! is recorded on every row because in all twelve of the module's pinching runs
//! `chi_measured - chi_predicted` equalled `non_manifold_edges` **exactly**; a
//! zero there is the positive statement that no sample landed on the isosurface,
//! and a non-zero one names how much of any χ gap is pinching rather than wrap.
//!
//! # Arms
//!
//! Five grids × seven extractors × two wrap modes = **70 rows**. Both wrap modes
//! come from **one** extraction: the open reading is taken before `wrap_seams`
//! and the periodic reading after it, so the two arms differ in exactly one
//! operation and in nothing else — no second extraction, no second grid, no
//! chance of comparing two different meshes.
//!
//! | arm | `periods_per_axis` | `voxels_per_period` | `resolution` | `cells` | `chi_predicted` | `is_control` |
//! |---|---|---|---|---|---|---|
//! | `period_sweep` | 1 | 33 | 34 | 35,937 | `-8` | no |
//! | `period_sweep` | 2 | 33 | 67 | 287,496 | `-64` | no |
//! | `period_sweep` | 3 | 33 | 100 | 970,299 | `-216` | no |
//! | `resolution_sweep` | 1 | 65 | 66 | 274,625 | `-8` | no |
//! | `resolution_sweep` | 1 | 97 | 98 | 912,673 | `-8` | no |
//!
//! `period_sweep` is C1's registered range `N ∈ {1,2,3}`. `resolution_sweep` is
//! the hypothesis's *"stable across resolution"* clause, taken at `N = 1` because
//! that is the claim — the same surface read on three grids — and because
//! `subgrid_marching_tetrahedra` costs 576 field evaluations per cell (`M-98`,
//! `M-248`) and a full cross product of periods and resolutions on it would be
//! twenty minutes of the phase's budget for no additional question.
//!
//! Every row carries `wrap_mode = open` or `periodic`; the `open` half **is** the
//! registered vacuity control and is marked `is_control = true`.
//!
//! # C1's scope is a measurement, not a courtesy
//!
//! C1 reads *"on every extractor in the crate that produces a closed surface"*,
//! so which extractors are in scope is itself a result, and this harness records
//! it per row in `is_closed` and `c1_scope` rather than dropping rows quietly.
//! The wrap is a **primal** operation: `wrap_seams` folds a coordinate within
//! `tol` of the far face onto the near one and welds, so it can only identify
//! vertices that actually lie **on** the domain boundary plane.
//!
//! * `marching_cubes`, `marching_cubes+decider`, `marching_tetrahedra` and
//!   `subgrid_marching_tetrahedra` interpolate on grid edges and tetrahedron
//!   edges, and their boundary vertices sit exactly on the plane.
//! * `surface_nets`, `dual_contouring` and `manifold_dual_contouring` place one
//!   vertex per **cell** and emit a quad per interior grid edge, so no vertex
//!   lands on the boundary plane and no quad spans it. `seam_pairs_identified`
//!   is expected to be **0** for all three, the wrap is a no-op, and their
//!   periodic row is numerically their open row — which is exactly why the
//!   column is recorded: a wrap that identified nothing must be visible as a
//!   zero rather than inferred from a disagreeing χ.
//!
//! An extractor out of scope is recorded with its `boundary_edges`, not skipped.
//! A skip that is recorded is a result; a skip that is silent is not.
//!
//! # SHARE, recomputed before the numbers
//!
//! **None, and the registration says so: *"SHARE: none — this is a correctness
//! gate where there was none."*** Nothing here is on an extraction path and
//! nothing here is a speedup claim, so `✗51`'s share bar does not apply and
//! `M-280`/`✗24`'s governor scatter has nothing to bite on: **every clause is an
//! integer equality over an enumerated population** — χ against `-8N³`, genus
//! against `1 + 4N³`, and `2 - 2g` against χ. Following `P-112`'s lesson that a
//! figure no clause reads should not be recorded at all, **no wall clock is
//! recorded**; the harness prints its own elapsed time for the operator's budget
//! and puts no timing column in the CSV.
//!
//! # C3 is answered as a reported result, and the falsifier's reason is too strong
//!
//! C3 asks that the oracle *"is added to the validity suite as a gate for
//! `gyroid`"*. It is **not** added, so `c3_holds` is `false` on every row, and
//! `crates/isomesh/src/**` is unchanged by this row. Its falsifier anticipates
//! this — *"C3 by the gate being unimplementable within the existing suite"* — but
//! that reason is **too strong**, and naming the change is the useful half of the
//! answer. Two things are needed and neither is what the falsifier implies:
//!
//! 1. **A periodic domain on `ValidateConfig`.** `validate_indexed(positions,
//!    indices, cfg)` counts `V - E + F` over the buffer it is handed, and
//!    `ValidateConfig::from_cell_size` (`validate.rs:102`) is its only
//!    constructor: it carries a cell size and nothing about the box. The gate
//!    needs `(lo, hi, wrap: [bool; 3])` and a fold in front of the weld it
//!    already performs — a second constructor plus one pass.
//!    **It is not a change to `Extractor::extract_into`**, and this row measures
//!    that: `wrap_seams` is a pure post-pass over `MeshBuffer` that needs no
//!    cooperation from any extractor and closes the primal extractors' output
//!    unmodified. `extractor.rs:63-70` states that the extraction signature *"is
//!    not negotiable"*; this row's finding is that it does not have to be
//!    negotiated.
//! 2. **An uncapped, period-conforming reference field.** The shipped `gyroid`
//!    entry is `capped_gyroid()` (`fields/mod.rs:1044-1060`), whose χ is not
//!    `-8N³` at all — the cap exists precisely because the uncapped surface has
//!    boundary in a finite box. A gate needs a new field on `[0, 2πN]³` with
//!    `expected_euler() = Some(-8N³)`, and one new reference field adds 27 rows to
//!    `golden_hashes.json` (216 = 8 fields × 9 algorithms × 3 resolutions) and
//!    moves `doc_facts.sh`'s gated `FIELDS` and `HASHES` counts in twelve
//!    documents.
//!
//! So the gate is implementable in `validate`; what it is not is implementable
//! inside a measurement commit. That is a Phase 28 landing ticket with the
//! golden-hash ripple priced into it, not a shortfall of this row.
//!
//! # What the wrap costs, and which readings survive it
//!
//! `common::tpms`'s header states it and this harness obeys it: a closed surface
//! on the 3-torus cannot be embedded in a box, so after `wrap_seams` some
//! triangles have one corner at each side of the domain. The buffer is then a
//! valid *simplicial complex* and an invalid *geometric* mesh — **connectivity
//! readings (χ, components, genus, boundary and non-manifold edges) are exact and
//! every metric reading (area, mean ratio, Hausdorff, self-intersections) is
//! nonsense.** This bench therefore reads only connectivity, and calls no
//! `accuracy`, no `field_bound_report` and no `self_intersections` — which would
//! be meaningless here twice over, since the nodal function is a level set and
//! not a distance (`|∇F|` vanishes on the whole singular skeleton).
//!
//! Two independent readers are run on every row and cross-checked:
//! `common::tpms::euler` and the crate's own `validate::validate_indexed`, both
//! welding at `weld::epsilon_for(cell_size)`. Their χ must agree exactly or the
//! instrument, not the hypothesis, is what failed — asserted, not recorded.
//!
//! # Vacuity controls
//!
//! * **The `-8` is derived, not assumed.** `body_centring_check(Gyroid)` must
//!   report the `(π,π,π)` invariance holding *and* the negation failing, and
//!   `primitive_cells_per_cubic_cell` must be 2. Proves the column
//!   `chi_per_cubic_cell = -8` is earned rather than transcribed; if the shift
//!   merely negated `F_G` the prediction would be `-4`.
//! * **The surface exists.** Every row must have `faces > 0`, or every count is
//!   trivially zero and χ is 0 by vacuity rather than by measurement.
//! * **C1's scope is non-empty, and the wrap is what put it there.** At least one
//!   closed periodic row must exist, or C1 is a claim about the empty set
//!   (`M-44`); and every closed periodic row must have `seam_pairs_identified >
//!   0`, or `periodic` and `open` are one computation and the contrast is vacuous.
//! * **The registered control: the non-wrapped arm must be run and must fail.**
//!   One `open` row per `periodic` row; **every** open row must have
//!   `boundary_edges > 0` — the module measured that the open arm is recognised
//!   by its boundary and *not* by its χ, because non-wrapped Schwarz P hits its
//!   own prediction by coincidence — and **no** open row may agree with `-8N³`.

#![allow(
    clippy::cast_possible_truncation,
    clippy::cast_precision_loss,
    clippy::too_many_lines
)]

mod common;

use std::time::Instant;

use isomesh::MeshBuffer;
use isomesh::extractor::Extractor;
use isomesh::validate::{ValidateConfig, validate_indexed};

use common::tpms::{self, EulerCount, NodalTpms, Tpms};

/// The `wrap_mode` value for a seam-identified extraction.
const WRAP_PERIODIC: &str = "periodic";
/// The `wrap_mode` value for the registered non-wrapped control arm.
const WRAP_OPEN: &str = "open";

/// Entries `isomesh::for_each_extractor!` visits, used only to size the row
/// buffer; the macro itself is the roster and this number cannot select from it.
const EXTRACTOR_COUNT: usize = isomesh::extractor::ALL_EXTRACTORS.len();

/// One `(periods, voxels_per_period)` grid and the arm it answers for.
#[derive(Clone, Copy, Debug)]
struct Config {
    /// `N` — periods per axis, so the box is `N³` conventional cubic cells.
    periods: u32,
    /// Cells spanning one period. Odd, for the reason in the header.
    voxels: u32,
    /// `period_sweep` for C1's `N ∈ {1,2,3}`, `resolution_sweep` for stability.
    arm: &'static str,
}

/// The five grids. `period_sweep` is C1's registered range; `resolution_sweep`
/// is the hypothesis's resolution-stability clause at fixed `N`.
const CONFIGS: [Config; 5] = [
    Config {
        periods: 1,
        voxels: 33,
        arm: "period_sweep",
    },
    Config {
        periods: 2,
        voxels: 33,
        arm: "period_sweep",
    },
    Config {
        periods: 3,
        voxels: 33,
        arm: "period_sweep",
    },
    Config {
        periods: 1,
        voxels: 65,
        arm: "resolution_sweep",
    },
    Config {
        periods: 1,
        voxels: 97,
        arm: "resolution_sweep",
    },
];

/// One CSV row: one extractor's reading of one grid under one wrap mode.
#[derive(Clone, Debug)]
struct Row {
    /// Which sweep this grid belongs to.
    arm: &'static str,
    /// `for_each_extractor!`'s entry name.
    extractor: &'static str,
    /// `N`.
    periods: u32,
    /// Cells per period.
    voxels: u32,
    /// Samples per axis, `voxels * periods + 1`.
    samples: u32,
    /// `(voxels * periods)³`.
    cells: u64,
    /// `2π / voxels`.
    cell_size: f64,
    /// `weld::epsilon_for(cell_size)`, the seam and weld tolerance.
    weld_tolerance: f64,
    /// `periodic` or `open`.
    wrap_mode: &'static str,
    /// Vertex pairs `wrap_seams` identified; `0` on every open row by definition.
    seam_pairs: u64,
    /// `-8N³`.
    chi_predicted: i64,
    /// `common::tpms::euler`'s reading — this row's χ and edge counts.
    counted: EulerCount,
    /// `MeshReport::euler_characteristic`, the independent reader.
    validate_chi: i64,
    /// `MeshReport::genus`; `None` where the crate declines to name one.
    genus: Option<i64>,
    /// `MeshReport::components`.
    components: u64,
    /// `MeshReport::boundary_edges`.
    validate_boundary_edges: u64,
    /// `MeshReport::non_manifold_edges`.
    validate_nonmanifold_edges: u64,
    /// `MeshReport::is_closed()` — C1's own scope predicate.
    closed: bool,
    /// `MeshReport::is_manifold()`.
    manifold: bool,
    /// Buffer vertices at the time of this reading.
    mesh_vertices: usize,
    /// Buffer triangles at the time of this reading.
    mesh_triangles: usize,
}

impl Row {
    /// Is this the periodic arm?
    fn periodic(&self) -> bool {
        self.wrap_mode == WRAP_PERIODIC
    }

    /// `chi_measured == chi_predicted`, the registered per-row column.
    fn chi_agreement(&self) -> bool {
        self.counted.chi == self.chi_predicted
    }

    /// `chi_measured - chi_predicted`, whose comparison with `nonmanifold_edges`
    /// is `common::tpms`'s pinching diagnostic.
    fn chi_gap(&self) -> i64 {
        self.counted.chi - self.chi_predicted
    }

    /// Is this row inside C1's registered scope — a periodic reading of an
    /// extractor that produced a closed surface?
    fn in_c1_scope(&self) -> bool {
        self.periodic() && self.closed
    }

    /// Why this row is or is not in C1's scope.
    fn scope(&self) -> &'static str {
        if !self.periodic() {
            "control_open_arm"
        } else if self.closed {
            "in_scope"
        } else {
            "excluded_not_closed"
        }
    }

    /// `1 + 4N³`, C2's prediction.
    fn genus_predicted(&self) -> i64 {
        let n = i64::from(self.periods);
        1 + 4 * n * n * n
    }

    /// `2 - 2g` from the measured genus — C2's consistency relation with C1.
    fn chi_from_genus(&self) -> Option<i64> {
        self.genus.map(|g| 2 - 2 * g)
    }

    /// C2 for this row: the genus is `1 + 4N³` **and** `χ = 2 - 2g` closes back
    /// onto the measured χ. Both halves, because the registration asks for the
    /// genus *and* for its consistency with C1.
    fn c2_consistent(&self) -> bool {
        self.genus == Some(self.genus_predicted())
            && self.chi_from_genus() == Some(self.counted.chi)
    }
}

/// Extract every grid with every extractor, reading χ before and after the wrap.
fn sweep() -> Vec<Row> {
    let mut rows: Vec<Row> = Vec::with_capacity(CONFIGS.len() * EXTRACTOR_COUNT * 2);

    for config in CONFIGS {
        let field = NodalTpms::new(Tpms::Gyroid, config.periods);
        let (lo, hi) = field.domain();
        let (shape, origin, cell_size) = field.periodic_grid(config.voxels);
        let samples = config.voxels * config.periods + 1;
        let cells = u64::from(config.voxels * config.periods).pow(3);
        let tol = isomesh::weld::epsilon_for(cell_size);
        let cfg = ValidateConfig::from_cell_size(cell_size)
            .expect("a periodic cell size is finite and positive");
        let chi_predicted = field.chi_predicted();

        // Inline blocks, so no `return` in here (M-253).
        isomesh::for_each_extractor!(f64, |name, extractor| {
            let mut mesh = MeshBuffer::<f64>::new();
            extractor
                .extract_into(&field, &shape, origin, cell_size, &mut mesh)
                .expect("the periodic grid holds at least two samples on every axis");

            // Both arms from one extraction: the open reading is taken with the
            // wrap withheld, so the arms differ in one operation and no more.
            let push = |rows: &mut Vec<Row>,
                        wrap_mode: &'static str,
                        seam_pairs: u64,
                        mesh: &MeshBuffer<f64>| {
                let counted = tpms::euler(&mesh.positions, &mesh.indices, tol);
                let report = validate_indexed(&mesh.positions, &mesh.indices, &cfg);
                assert_eq!(
                    counted.chi, report.euler_characteristic,
                    "{name} {wrap_mode} N={} v={}: common::tpms::euler and \
                     validate_indexed disagree about V-E+F at one weld tolerance, \
                     so the instrument is what failed, not the hypothesis",
                    config.periods, config.voxels
                );
                rows.push(Row {
                    arm: config.arm,
                    extractor: name,
                    periods: config.periods,
                    voxels: config.voxels,
                    samples,
                    cells,
                    cell_size,
                    weld_tolerance: tol,
                    wrap_mode,
                    seam_pairs,
                    chi_predicted,
                    counted,
                    validate_chi: report.euler_characteristic,
                    genus: report.genus,
                    components: report.components,
                    validate_boundary_edges: report.boundary_edges,
                    validate_nonmanifold_edges: report.non_manifold_edges,
                    closed: report.is_closed(),
                    manifold: report.is_manifold(),
                    mesh_vertices: mesh.vertex_count(),
                    mesh_triangles: mesh.triangle_count(),
                });
            };

            push(&mut rows, WRAP_OPEN, 0, &mesh);
            let seam_pairs = tpms::wrap_seams(&mut mesh, lo, hi, tol);
            push(&mut rows, WRAP_PERIODIC, seam_pairs, &mesh);
        });
    }

    rows
}

fn main() {
    if !std::env::args().any(|arg| arg == "--bench") {
        return;
    }
    let prereg = isomesh::experiment!("P-142");

    common::experiment::run(prereg, |run| {
        let started = Instant::now();
        let rows = sweep();
        println!(
            "-- swept {} rows over {} grids x {EXTRACTOR_COUNT} extractors x 2 wrap modes in {:.1} s",
            rows.len(),
            CONFIGS.len(),
            started.elapsed().as_secs_f64()
        );

        // ── vacuity controls, before any record ─────────────────────────────
        //
        // 1. The `-8` is derived from the body-centring invariance, not assumed.
        let kind = Tpms::Gyroid;
        let (residual, centring_ok) = tpms::body_centring_check(kind);
        assert!(
            centring_ok && residual <= tpms::SHIFT_RESIDUAL_TOLERANCE,
            "VOID: the gyroid's body-centring shift (pi,pi,pi) does not check out \
             (residual {residual:.3e}, claim held {centring_ok}), so the two \
             primitive cells per conventional cubic cell are unproven and every \
             `chi_predicted = -8N^3` in this CSV is a transcription rather than a \
             derivation -- if the shift merely negated F_G the prediction would be \
             -4 and C1 would be graded against the wrong number"
        );
        assert_eq!(
            kind.primitive_cells_per_cubic_cell(),
            2,
            "VOID: the gyroid's conventional cubic cell must hold two primitive \
             cells for chi = -4 * 2 = -8; one would make the prediction -4"
        );
        assert_eq!(
            kind.chi_per_cubic_cell(),
            -8,
            "VOID: chi per conventional cubic cell must be -8, or this harness is \
             not measuring P-142's prediction"
        );

        // 2. Every row read a surface that exists.
        for row in &rows {
            assert!(
                row.counted.faces > 0,
                "VOID: {} {} at N={} v={} meshed the gyroid to zero faces, so its \
                 V, E and F are all zero and its chi of {} is vacuous rather than \
                 measured (M-44)",
                row.extractor,
                row.wrap_mode,
                row.periods,
                row.voxels,
                row.counted.chi
            );
        }

        // 3. C1's scope is non-empty, and the wrap is what put it there.
        let closed_periodic: Vec<&Row> = rows.iter().filter(|r| r.in_c1_scope()).collect();
        assert!(
            !closed_periodic.is_empty(),
            "VOID: not one extractor produced a closed surface under periodic \
             wrap, so C1's population -- 'every extractor in the crate that \
             produces a closed surface' -- is empty and C1 would hold for the \
             wrong reason (M-44)"
        );
        for row in &closed_periodic {
            assert!(
                row.seam_pairs > 0,
                "VOID: {} closed at N={} v={} while wrap_seams identified zero \
                 vertex pairs, so the periodic arm is byte-for-byte the open arm \
                 and this row's agreement with -8N^3 says nothing about \
                 periodicity",
                row.extractor,
                row.periods,
                row.voxels
            );
        }

        // 4. The registered vacuity control: the non-wrapped arm ran, is
        //    recognisable as non-wrapped, and fails the -8N^3 prediction.
        let open: Vec<&Row> = rows.iter().filter(|r| !r.periodic()).collect();
        let periodic: Vec<&Row> = rows.iter().filter(|r| r.periodic()).collect();
        assert!(
            !open.is_empty() && open.len() == periodic.len(),
            "VOID: the non-wrapped control arm must be run once per periodic row \
             ({} open against {} periodic), or the experiment has not shown that \
             periodicity is what matters",
            open.len(),
            periodic.len()
        );
        for row in &open {
            assert!(
                row.counted.boundary_edges > 0,
                "VOID: the non-wrapped {} at N={} v={} has zero boundary edges, so \
                 it is not recognisable as non-wrapped and is not a control; \
                 common::tpms measured that the open arm must be identified by its \
                 boundary and never by its chi, because non-wrapped Schwarz P hits \
                 its own prediction by coincidence",
                row.extractor,
                row.periods,
                row.voxels
            );
            assert!(
                !row.chi_agreement(),
                "VOID: the non-wrapped {} at N={} v={} reproduced chi = {} = -8N^3 \
                 without any seam identification, so this fixture does not \
                 discriminate the wrap and the registered vacuity control has \
                 failed",
                row.extractor,
                row.periods,
                row.voxels,
                row.counted.chi
            );
        }

        // ── verdicts ────────────────────────────────────────────────────────
        //
        // All three clauses are global claims over an enumerated population, so
        // the same verdict is stamped on every row and the per-row facts live in
        // `chi_agreement`, `genus_agreement`, `is_closed` and `c1_scope`.
        let c1 = closed_periodic.iter().all(|r| r.chi_agreement());
        let c2 = closed_periodic.iter().all(|r| r.c2_consistent());
        // C3: the oracle is NOT added to the validity suite. Reported, with the
        // two-part change named in this file's header. `crates/isomesh/src/**` is
        // unchanged by this row, so this is false by construction and the
        // constant says so in exactly one place.
        let c3 = false;

        let non_closed: Vec<&Row> = periodic.iter().filter(|r| !r.closed).copied().collect();
        println!(
            "\n-- C1 population: {} closed periodic rows",
            closed_periodic.len()
        );
        for row in &closed_periodic {
            println!(
                "   {:<28} N={} v={:<3} chi {:>5} vs {:>5}  genus {:>4} vs {:>4}  pairs {:>5}",
                row.extractor,
                row.periods,
                row.voxels,
                row.counted.chi,
                row.chi_predicted,
                row.genus
                    .map_or_else(|| "none".to_string(), |g| g.to_string()),
                row.genus_predicted(),
                row.seam_pairs
            );
        }
        println!(
            "\n-- out of C1's scope: {} periodic rows did not close",
            non_closed.len()
        );
        for row in &non_closed {
            println!(
                "   {:<28} N={} v={:<3} chi {:>5} vs {:>5}  boundary {:>6}  nm {:>4}  pairs {:>5}",
                row.extractor,
                row.periods,
                row.voxels,
                row.counted.chi,
                row.chi_predicted,
                row.counted.boundary_edges,
                row.counted.non_manifold_edges,
                row.seam_pairs
            );
        }
        println!(
            "\n-- verdicts: c1_holds={c1} c2_holds={c2} c3_holds={c3} (C3 is not \
             added to the validity suite; see this bench's header for the \
             ValidateConfig constructor and the uncapped reference field it would \
             need)"
        );

        // ── rows ────────────────────────────────────────────────────────────
        for row in &rows {
            run.record(&[
                ("field", kind.name().to_string()),
                ("periods_per_axis", row.periods.to_string()),
                ("resolution", row.samples.to_string()),
                ("wrap_mode", row.wrap_mode.to_string()),
                ("chi_predicted", row.chi_predicted.to_string()),
                ("chi_measured", row.counted.chi.to_string()),
                ("chi_agreement", row.chi_agreement().to_string()),
                (
                    "genus_measured",
                    row.genus
                        .map_or_else(|| "none".to_string(), |g| g.to_string()),
                ),
                ("boundary_edges", row.counted.boundary_edges.to_string()),
                (
                    "nonmanifold_edges",
                    row.counted.non_manifold_edges.to_string(),
                ),
                ("cells", row.cells.to_string()),
                ("c1_holds", c1.to_string()),
                ("c2_holds", c2.to_string()),
                ("c3_holds", c3.to_string()),
                // ── extras (M-273) ──
                ("arm", row.arm.to_string()),
                ("body_centring_ok", centring_ok.to_string()),
                ("body_centring_residual", format!("{residual:.3e}")),
                ("c1_scope", row.scope().to_string()),
                ("cell_size", format!("{:.9}", row.cell_size)),
                (
                    "chi_from_genus",
                    row.chi_from_genus()
                        .map_or_else(|| "none".to_string(), |c| c.to_string()),
                ),
                ("chi_gap", row.chi_gap().to_string()),
                ("chi_per_cubic_cell", kind.chi_per_cubic_cell().to_string()),
                ("components", row.components.to_string()),
                ("edges", row.counted.edges.to_string()),
                ("extractor", row.extractor.to_string()),
                ("faces", row.counted.faces.to_string()),
                ("genus_agreement", row.c2_consistent().to_string()),
                ("genus_predicted", row.genus_predicted().to_string()),
                ("is_closed", row.closed.to_string()),
                ("is_control", (!row.periodic()).to_string()),
                ("is_manifold", row.manifold.to_string()),
                ("mesh_triangles", row.mesh_triangles.to_string()),
                ("mesh_vertices", row.mesh_vertices.to_string()),
                (
                    "primitive_cells_per_cubic_cell",
                    kind.primitive_cells_per_cubic_cell().to_string(),
                ),
                ("primitive_lattice", kind.primitive_lattice().to_string()),
                ("seam_pairs_identified", row.seam_pairs.to_string()),
                ("space_group", kind.space_group().to_string()),
                (
                    "validate_boundary_edges",
                    row.validate_boundary_edges.to_string(),
                ),
                ("validate_chi", row.validate_chi.to_string()),
                (
                    "validate_nonmanifold_edges",
                    row.validate_nonmanifold_edges.to_string(),
                ),
                ("vertices", row.counted.vertices.to_string()),
                ("voxels_per_period", row.voxels.to_string()),
                ("weld_tolerance", format!("{:.3e}", row.weld_tolerance)),
            ]);
        }

        println!(
            "\n-- total wall {:.1} s (no timing column is recorded: no clause \
             reads one, and P-112 retired the figure that none did)",
            started.elapsed().as_secs_f64()
        );
    });
}
