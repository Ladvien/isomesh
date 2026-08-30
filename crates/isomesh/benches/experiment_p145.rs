//! **P-145 — a ground-truth Euler characteristic read off the field's own signs, on the three fields the crate declares unknowable.**
//!
//! Ticket: R-145. Pre-registered before this harness existed.
//!
//! ```bash
//! cargo bench --bench experiment_p145
//! ```
//!
//! Writes `docs/experiments/p-145.csv`.
//!
//! # What was missing
//!
//! **Every χ in this crate is read off a mesh.** `validate.rs:1020-1021` computes
//! `euler_characteristic = referenced_vertices − edges + faces`, and that is the
//! only χ the crate can produce; the only other source is a literal a field
//! *declares*, `ReferenceField::expected_euler` (`fields/mod.rs:169`). Where the
//! literal is `None` there is nothing to compare a mesh against, and the crate
//! says so out loud:
//!
//! * `fields/mod.rs:1078` — `capped_gyroid`: *"genus depends on how many tunnels
//!   the cap encloses"*.
//! * `fields/mod.rs:1256` — `noise_cavity`: *"genus depends on how many noise
//!   blobs the cap encloses"*.
//! * `fields/mod.rs:1387` — `fbm_terrain`: *"not closed, so there is nothing to
//!   assert"*.
//! * `marching_cubes/tests.rs:347-353` states the consequence as policy: the
//!   capped gyroid's χ *"is not known analytically, so it is recorded rather than
//!   asserted — exactly what `expected_euler() == None` is telling the harness to
//!   do"*, and `:362-366` asserts only that it is **even**.
//!
//! So on three of eight reference fields the topological gate is "record whatever
//! the mesh said". A defect that changes the genus on those fields is invisible to
//! everything except `golden_hashes.json`, which cannot say *which* number moved
//! or whether the old one was right.
//!
//! `P-142`/`R-142` closed exactly one of those holes and did it from the other
//! side: it derived the gyroid's χ **analytically** (`-8N³` per `N³` conventional
//! cubic cells, from genus 3 per primitive cell and two primitive cells per cubic
//! cell) and measured it, `docs/experiments/p-142.csv:7,21,35` giving
//! `chi_measured = -8, -64, -216` against `chi_predicted = -8, -64, -216` on the
//! wrapped `marching_cubes` rows at `N = 1, 2, 3`. That is a closed-form
//! derivation for one surface family. It does not generalise: there is no closed
//! form for `noise_cavity`, and there never will be.
//!
//! This row supplies the missing instrument — a χ computed from **the field's
//! signs on the grid**, with no mesh involved — and validates it against the one
//! number `P-142` already established.
//!
//! ## The citation, and the flag on it
//!
//! `docs/research/2026-08-26-audit-and-phase-23-registrations.md:365` lists
//! *"Etiene et al., Topology Verification for Isosurface Extraction, IEEE TVCG
//! 2012 (**DOI unverified**; PubMed 21690649)"*. The registration asserts the DOI
//! is `10.1109/tvcg.2011.109`. **It is, and the resolution was performed rather
//! than quoted** — three independent routes, all agreeing:
//!
//! | route | returned |
//! |---|---|
//! | Crossref + Europe PMC + OpenAlex + CORE, by DOI | *Topology verification for isosurface extraction*, Etiene T, Nonato LG, Scheidegger C, Tierny J, Peters TJ, Pascucci V, Kirby RM, Silva CT |
//! | OpenAlex `W2086220855`, by DOI | *Topology Verification for Isosurface Extraction*, 2011 online, `is_oa`, primary source `S84775595` |
//! | IEEE / NYU Scholars / PubMed `21690649` | IEEE TVCG **18(6), 952–965, 2012** |
//!
//! The abstract names the method this harness implements: *"we use stratified
//! Morse theory and digital topology to design algorithms which verify topological
//! invariants."* `doi_verified` is `true` and `doi_resolved_title` carries the
//! title, so the CSV is the record rather than this comment.
//!
//! # The method, and which half of it this is
//!
//! Etiene et al. give **two** oracles (their §4.1 and §4.2, and their Figure 1
//! pipeline runs one or the other):
//!
//! 1. **Digital topology** — build a digital object from the field's signs and
//!    read its invariants off the voxel structure. Their Theorem 4.1: if no cubic
//!    cell of the grid is *ambiguous*, the digital surface is **homeomorphic** to
//!    the level set of the trilinear interpolant. Restricted, in their words, to
//!    *"isosurfaces without boundaries"*.
//! 2. **Stratified Morse theory** — sum `Δχ` over the critical points of the
//!    stratified trilinear function. Handles surfaces **with** boundary, and
//!    yields χ only, not the Betti numbers.
//!
//! **This harness implements (1).** The registered column is named
//! `chi_stratified_morse` and that name is the registration's, not this file's; it
//! is honoured as *the ground-truth χ of Etiene et al.'s method*, and the row says
//! in `oracle_method` which of the paper's two halves produced it. Where (1) does
//! not apply — a level set with boundary — the column is recorded `none` rather
//! than filled with a number from a different space. That exclusion is C2's
//! subject, not a shortfall: C2 asks for *"the set of fields where it applies"*.
//!
//! Two deliberate departures from the paper's §4.1, both in the direction of
//! fewer moving parts:
//!
//! * **No Majority Interpolation, and no refinement.** The paper builds its
//!   digital object on a **doubled** grid via MI (their Figure 5) and refines
//!   until no cell is ambiguous, because it needs the digital *surface* to be a
//!   2-manifold so that Chen & Rong's neighbour-count genus formula (their
//!   Figure 4, `g = 1 + (|N₅| + 2|N₆| − |N₃|)/8`) and `χ = 2 − 2g` per component
//!   apply. This harness never computes a genus, so it never needs a manifold: it
//!   sums χ **additively over local configurations**, which is defined for any
//!   voxel set whatever. What the refinement bought is instead *reported*, as
//!   `disagreement_cells` — see below.
//! * **χ of the solid, doubled.** The oracle reads `χ` of the digital **solid**
//!   `{f < 0}`. For a compact 3-manifold `M` with boundary, `χ(∂M) = 2·χ(M)`, so
//!   the surface χ this row compares against a mesh is `2 · chi_solid_26`. The
//!   gyroid checks this arithmetic independently: its labyrinth per conventional
//!   cubic cell is the `srs` (Laves) net, 8 nodes and 12 edges, `χ = 8 − 12 = -4`,
//!   and `2 · (-4) = -8` is exactly `P-142`'s number.
//!
//! # The weights are derived here, by enumeration, and never transcribed
//!
//! The classical result (Ohser, Nagel & Schladitz) is that χ of a 3D digital
//! object is a **fixed integer combination of the counts of the 256 `2×2×2`
//! binary patterns**. The combination is not looked up. For each of the 256
//! patterns this file computes
//!
//! ```text
//! w(pattern) = Σ over the cells of the cubical complex that meet the block
//!                (-1)^dim(cell) / (number of blocks sharing that cell)
//! ```
//!
//! and `χ = Σ_blocks w(pattern of that block)` is then exactly `V − E + F − C` of
//! the whole complex, because every cell is counted once: `1/m` from each of its
//! `m` blocks. Blocks are indexed by lattice vertex, one per vertex, and the
//! sharing multiplicities are geometry, not convention —
//!
//! | connectivity | 0-cells | 1-cells | 2-cells | 3-cells |
//! |---|---|---|---|---|
//! | **foreground 26** (union of closed voxels) | the block's central lattice vertex, 8 voxels, `m = 1` | 6 lattice edges through it, 4 voxels each, `m = 2` | 12 lattice faces through it, 2 voxels each, `m = 4` | the 8 voxels, `m = 8` |
//! | **foreground 6** (complex on voxel centres) | the 8 voxels, `m = 8` | 12 axis-adjacent voxel pairs, `m = 4` | 6 planar `2×2` voxel squares, `m = 2` | the full block, `m = 1` |
//!
//! A 26-cell is present when **any** voxel sharing it is occupied (a union of
//! closed cubes); a 6-cell is present when **every** voxel it spans is occupied (a
//! clique complex on centres). Every multiplicity divides 8, so every weight is an
//! exact multiple of `1/8` and the whole computation is integer arithmetic over
//! `WEIGHT_DENOMINATOR`; the divisibility of the total is asserted, which is a
//! global check that the fractional parts cancelled.
//!
//! # The connectivity is stated, and the two answers differ
//!
//! **The headline is foreground 26 / background 6**, `connectivity = f26_b6`,
//! because that is the model the paper's digital object lives in: `O_a` is a *set
//! of voxels* and `∂O_a` is the boundary of their **union**, each voxel being the
//! Voronoi cell of its sample (their Figure 2). Under that model two voxels
//! touching at a corner are one object.
//!
//! **Foreground 6 is a different number on the same occupancy, and the row carries
//! both.** The smallest witness is two voxels sharing only a corner: as a union of
//! closed cubes that is a wedge of two cubes at a point, `χ = 1 + 1 − 1 = 1`; as a
//! complex on centres it is two components, `χ = 2`. The derived tables disagree
//! on **96 of the 256** configurations, and `connectivity_split_blocks` counts how
//! many blocks of each actual grid landed on one of those 96 — so a reader can see
//! whether the choice mattered *here* rather than in principle. Both are recorded:
//! `chi_solid_26` / `chi_surface_26` and `chi_solid_6` / `chi_surface_6`.
//!
//! # `disagreement_cells` is Etiene's own hypothesis, counted
//!
//! Theorem 4.1 holds *"if no cubic cell of `G` is ambiguous"*, and their §4.1
//! defines ambiguous exactly: a cell is **un**ambiguous when the corners with
//! `e < a` are joined by a path of negative edges inside the cell **and** the
//! corners with `e > a` are joined by a path of positive edges. `disagreement_cells`
//! is the number of cells of the sampled grid that fail that test — the cells
//! where the digital reconstruction carries **no** homeomorphism guarantee, and
//! equally the cells where an extractor's face and interior rules are free to
//! choose a different topology from the trilinear interpolant's. It is therefore
//! the column that *explains* a `chi_stratified_morse ≠ chi_extracted`, and its
//! being non-zero somewhere is what makes `agreement` a measurement rather than a
//! tautology.
//!
//! The 256-entry table is derived from the crate's **own** cube graph,
//! `marching_cubes::table::EDGE_CORNERS`, so it cannot disagree with the extractor
//! about which corners are adjacent; `corner_offset` is `pub(crate)` so the three
//! lines are rewritten here and then *checked* against the public `EDGE_CORNERS`
//! and `EDGE_AXIS` rather than trusted. The derived table is additionally
//! cross-checked against the shipped `AMBIGUOUS_FACES`: every case the crate marks
//! with an ambiguous face must be ambiguous by Etiene's definition too.
//!
//! # Arms
//!
//! | arm | fields | grids | `is_control` | what it answers |
//! |---|---|---|---|---|
//! | `calibration` | the five reference fields with `expected_euler() == Some(χ)`: `sphere`, `torus`, `box_exact`, `csg_difference`, `thin_plate` | 33³, 65³ | **yes** | the registered vacuity control (`sphere → 2`) and the assigned second one (`torus → 0`), plus three further known χ recorded |
//! | `new_ground_truth` | `gyroid` (`capped_gyroid`) and `noise_cavity` — `expected_euler() == None`, closed in domain | 33³, 65³ | no | C1's *"at least two fields where the crate currently has none"* |
//! | `out_of_scope` | `fbm_terrain` — `expected_euler() == None`, **not** closed | 33³, 65³ | no | C2's *"the set of fields where it applies is stated"*, with the exclusion **measured** |
//! | `p142_cross_check` | `gyroid_nodal`, the nodal gyroid on the 3-torus over `N` periods | `N ∈ {1,2,3}` at 33 voxels/period | no | C1's *"agrees with `P-142`'s analytic prediction on `gyroid`"* |
//!
//! Every row is meshed twice, `extractor ∈ {marching_cubes, marching_cubes+trilinear}`
//! — the crate's default (`FaceAmbiguity::Separate` + `InteriorAmbiguity::Ignore`)
//! and its ambiguity-resolving configuration (`AsymptoticDecider` + `Trilinear`),
//! the two names golden.rs uses. That is the comparison Etiene's paper exists to
//! make: an extractor that claims to reproduce the trilinear topology should agree
//! with the trilinear ground truth where the default is free not to.
//!
//! 33 and 65 are the contract's default pair: word-boundary sample counts, one
//! below and one above `thin_plate`'s construction cell size of `4/64`. The nodal
//! arm uses **33 voxels per period, odd on purpose** — `common::tpms`'s author
//! measured that a `voxels_per_period` divisible by 8 puts samples on the `π/4`
//! lattice where a nodal function cancels to exactly `0.0` (`M-48`'s degenerate
//! crossing) and the weld turns coincident vertices into a pinch. `33 · N + 1`
//! samples per axis reproduces `p-142.csv`'s `resolution` column exactly — 34, 67,
//! 100 — so the two CSVs' rows line up.
//!
//! # SHARE, recomputed before the numbers
//!
//! **None, and the registration says so: *"SHARE: none — verification cost, not
//! runtime cost."*** Nothing here runs on an extraction path. `✗51`'s share bar
//! therefore does not apply and no `*_share` column exists.
//!
//! `method_cost_ms` is nevertheless a wall clock, because C2 asks for a *cost*
//! and there is no integer proxy for "would this oracle be affordable in a test".
//! `M-280`/`✗24` measured this host's `amd-pstate-epp` governor swinging the same
//! binary 1.45×, so the clock is taken as the **median of 5 repeats** after one
//! warm-up, with `method_cost_ms_min`/`method_cost_ms_max` beside it so the
//! scatter is visible rather than averaged away. **No clause reads a threshold on
//! it** — C2 is falsified by inapplicability, not by slowness — so the governor
//! has nothing to bite on. `extract_cost_ms` and `cost_ratio` are recorded next to
//! it, because C2's sentence is *"an oracle that costs more than the mesh is a
//! test-only instrument and should be labelled one"* and that comparison needs
//! both sides.
//!
//! # Vacuity controls
//!
//! * **The registered one: `sphere` must return `χ = 2`.** Asserted at every
//!   resolution, on `chi_surface_26` and `chi_surface_6` both. Column:
//!   `chi_stratified_morse` on the `sphere` rows.
//! * **The assigned second: `torus` must return `χ = 0`.** Same, and it is the
//!   stronger of the two — a method that reported `2` for everything would pass
//!   the sphere. Column: `chi_stratified_morse` on the `torus` rows.
//! * **The oracle is not a constant.** At least three distinct
//!   `chi_stratified_morse` values across the in-scope rows, or `2` and `0` were
//!   reached by an instrument that cannot tell fields apart (`M-44`).
//! * **The weights are derived and the derivation is checkable.** `w(0x00) == 0`
//!   (empty block), `w(0xFF) == 0` (interior block — χ is a boundary phenomenon),
//!   a single occupied voxel weighs exactly `1/8` and eight such blocks make the
//!   `χ = 1` of a cube. Proves `WEIGHTS_26`/`WEIGHTS_6` are not a transcribed
//!   table with a typo in it.
//! * **The two connectivities genuinely differ.** The derived tables must disagree
//!   on at least one configuration, **and** the corner-touching pair fixture must
//!   read `χ_26 = 1` and `χ_6 = 2`. A `connectivity` column naming a choice that
//!   changes no answer is not a statement about anything.
//! * **The digital object is neither empty nor full.** `inside_samples` strictly
//!   between `0` and the sample count on every row, or every χ is `0` by vacuity.
//! * **The mesh exists.** `mesh_triangles > 0` on every row.
//! * **`disagreement_cells` could have been non-zero.** At least one row must have
//!   an ambiguous cell, or the column is measuring nothing and Theorem 4.1's
//!   hypothesis was never tested.
//! * **The exclusion is measured, not declared.** `fbm_terrain` must have
//!   `boundary_inside > 0` (its digital object reaches the sampling box) **and**
//!   `boundary_edges > 0` (its mesh is open); every in-scope field must have both
//!   at zero. Two independent readings of "closed", from the field side and the
//!   mesh side, and they must agree — otherwise `applicable_fields` is an opinion.
//! * **The periodic identification is exact.** `wrap_sign_mismatches == 0`: the
//!   sign at sample `n` must equal the sign at sample `0` on every axis, or the
//!   torus occupancy is glued to the wrong thing. `sin(2πN)` is `-2.4e-16` and not
//!   `0.0` in `f64`, so this is a real hazard and the occupancy is built by
//!   wrapping indices rather than by trusting it.
//! * **The DOI column is the registration's DOI.** `prereg.hypothesis` must
//!   contain the string recorded in `doi`, so the resolution reported here is a
//!   resolution *of the registered claim* and not of a retyped one.
//! * **The cube-corner convention is checked, not assumed.** The rewritten
//!   `corner_offset` must reproduce `EDGE_CORNERS`/`EDGE_AXIS`, and every case the
//!   crate's `AMBIGUOUS_FACES` marks must be ambiguous by the derived table.

#![allow(
    clippy::cast_possible_truncation,
    clippy::cast_precision_loss,
    clippy::too_many_lines
)]

mod common;

use std::time::Instant;

use isomesh::fields::ReferenceField;
use isomesh::marching_cubes::table::{
    AMBIGUOUS_FACES, CORNER_COUNT, EDGE_AXIS, EDGE_CORNERS, EDGE_COUNT, is_inside,
};
use isomesh::marching_cubes::{FaceAmbiguity, InteriorAmbiguity, MarchingCubes};
use isomesh::validate::{ValidateConfig, validate_indexed};
use isomesh::{MeshBuffer, Sdf, Shape3};

use common::tpms::{self, NodalTpms, Tpms};

// ─── the citation, resolved ──────────────────────────────────────────────────

/// Etiene et al. 2012's DOI, as the registration asserts it.
///
/// Checked against `prereg.hypothesis` before any row is written, so that the
/// column reports a resolution of the *registered* claim.
const DOI: &str = "10.1109/tvcg.2011.109";

/// The title the DOI resolves to. Recorded rather than described, because
/// `2026-08-26-audit-and-phase-23-registrations.md:365` flagged the DOI unverified
/// and a finding that says "verified" without the title has verified nothing.
const DOI_RESOLVED_TITLE: &str = "Topology Verification for Isosurface Extraction";

/// Venue, formatted without a comma because the CSV writer refuses one.
const DOI_VENUE: &str = "IEEE-TVCG-18(6)-952-965-2012";

/// OpenAlex work id and PubMed id, the two secondary handles that agreed.
const DOI_HANDLES: &str = "openalex:W2086220855|pubmed:21690649";

// ─── the enumerated tables ───────────────────────────────────────────────────

/// Common denominator of every Ohser–Nagel–Schladitz weight.
///
/// The sharing multiplicities of the cells of a cubical complex are 1, 2, 4 and
/// 8, so every weight is an exact multiple of an eighth and the whole sum can be
/// integer arithmetic. Nothing else about the number is chosen.
const WEIGHT_DENOMINATOR: i64 = 8;

/// The voxel offsets of one `2×2×2` block, in bit order.
///
/// Bit `b` of a configuration index is the voxel at `BLOCK_OFFSETS[b]` from the
/// block's base voxel. Both the weight derivation and the block census read this
/// one constant, so they cannot disagree about which bit is which voxel.
const BLOCK_OFFSETS: [[usize; 3]; 8] = [
    [0, 0, 0],
    [0, 0, 1],
    [0, 1, 0],
    [0, 1, 1],
    [1, 0, 0],
    [1, 0, 1],
    [1, 1, 0],
    [1, 1, 1],
];

/// Which digital connectivity the foreground `{f < 0}` is given.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Connectivity {
    /// Foreground 26-connected, background 6-connected: the object is the union
    /// of the closed Voronoi voxels of the inside samples. Etiene et al.'s model.
    Foreground26,
    /// Foreground 6-connected, background 26-connected: the object is the cubical
    /// complex on the inside samples' centres.
    Foreground6,
}

impl Connectivity {
    /// The CSV token.
    fn name(self) -> &'static str {
        match self {
            Connectivity::Foreground26 => "f26_b6",
            Connectivity::Foreground6 => "f6_b26",
        }
    }
}

/// One cell of the cubical complex, as it appears in a `2×2×2` block.
#[derive(Clone, Copy, Debug)]
struct LocalCell {
    /// Which of the block's eight voxels the cell is incident to.
    voxels: u8,
    /// The cell's dimension; the Euler sum takes `(-1)^dimension`.
    dimension: u32,
    /// How many blocks of the whole lattice share this cell.
    shared: i64,
}

/// The bits of `BLOCK_OFFSETS` whose `axis` offset is `side`.
fn mask_where(axis: usize, side: usize) -> u8 {
    let mut mask = 0u8;
    for (bit, offset) in BLOCK_OFFSETS.iter().enumerate() {
        if offset[axis] == side {
            mask |= 1 << bit;
        }
    }
    mask
}

/// The bit of the voxel matching `bit`'s offsets except along `axis`, where it is
/// one rather than zero.
///
/// Found by searching `BLOCK_OFFSETS` so that nothing here depends on how the
/// offsets happen to be packed into a bit index.
fn axis_partner(bit: usize, axis: usize) -> usize {
    let mut want = BLOCK_OFFSETS[bit];
    want[axis] = 1;
    BLOCK_OFFSETS
        .iter()
        .position(|offset| *offset == want)
        .expect("every offset triple in {0,1}^3 is one of the block's voxels")
}

/// The cells of the union-of-closed-voxels complex that a block owns exactly one
/// `1/shared` of.
///
/// Blocks are indexed by lattice vertex: the block's own central lattice vertex
/// is the common corner of all eight of its voxels, and every lattice vertex is
/// the centre of exactly one block. So the vertex is owned outright, each of the
/// six lattice edges through it by two blocks, each of the twelve lattice faces
/// through it by four, and each voxel by the eight blocks at its corners.
fn cells_26() -> Vec<LocalCell> {
    let mut out = Vec::with_capacity(1 + 6 + 12 + 8);
    out.push(LocalCell {
        voxels: u8::MAX,
        dimension: 0,
        shared: 1,
    });
    for axis in 0..3 {
        for side in 0..2 {
            out.push(LocalCell {
                voxels: mask_where(axis, side),
                dimension: 1,
                shared: 2,
            });
        }
    }
    for normal in 0..3 {
        let u = (normal + 1) % 3;
        let v = (normal + 2) % 3;
        for su in 0..2 {
            for sv in 0..2 {
                out.push(LocalCell {
                    voxels: mask_where(u, su) & mask_where(v, sv),
                    dimension: 2,
                    shared: 4,
                });
            }
        }
    }
    for bit in 0..8 {
        out.push(LocalCell {
            voxels: 1 << bit,
            dimension: 3,
            shared: 8,
        });
    }
    out
}

/// The cells of the complex on voxel centres, with the same ownership accounting.
///
/// A voxel centre is a corner of the eight blocks around it, an axis-adjacent
/// pair is an edge of four, a planar `2×2` square is a face of two, and the full
/// `2×2×2` is a cube of exactly one.
fn cells_6() -> Vec<LocalCell> {
    let mut out = Vec::with_capacity(8 + 12 + 6 + 1);
    for bit in 0..8 {
        out.push(LocalCell {
            voxels: 1 << bit,
            dimension: 0,
            shared: 8,
        });
    }
    for axis in 0..3 {
        for bit in 0..8 {
            if BLOCK_OFFSETS[bit][axis] == 0 {
                out.push(LocalCell {
                    voxels: (1 << bit) | (1 << axis_partner(bit, axis)),
                    dimension: 1,
                    shared: 4,
                });
            }
        }
    }
    for axis in 0..3 {
        for side in 0..2 {
            out.push(LocalCell {
                voxels: mask_where(axis, side),
                dimension: 2,
                shared: 2,
            });
        }
    }
    out.push(LocalCell {
        voxels: u8::MAX,
        dimension: 3,
        shared: 1,
    });
    out
}

/// The 256 weights of `connectivity`, as numerators over [`WEIGHT_DENOMINATOR`].
///
/// This is the enumeration the header describes: `(-1)^dim / shared`, summed over
/// the cells the configuration activates. Nothing is transcribed.
fn weight_numerators(connectivity: Connectivity) -> [i64; 256] {
    let cells = match connectivity {
        Connectivity::Foreground26 => cells_26(),
        Connectivity::Foreground6 => cells_6(),
    };
    let mut out = [0i64; 256];
    for (configuration, weight) in out.iter_mut().enumerate() {
        let occupied = configuration as u8;
        for cell in &cells {
            let present = match connectivity {
                // A union of closed cubes holds a cell as soon as one of the
                // voxels sharing it is occupied.
                Connectivity::Foreground26 => occupied & cell.voxels != 0,
                // A complex on centres holds a cell only when every voxel it
                // spans is occupied.
                Connectivity::Foreground6 => occupied & cell.voxels == cell.voxels,
            };
            if present {
                let sign = if cell.dimension % 2 == 0 { 1 } else { -1 };
                *weight += sign * (WEIGHT_DENOMINATOR / cell.shared);
            }
        }
    }
    out
}

/// The `(x, y, z)` grid offset of cube corner `i`.
///
/// `isomesh`'s own `cube::corner_offset` is `pub(crate)`, so these are the three
/// lines the cheat sheet names. They are **checked** against the public
/// `EDGE_CORNERS` and `EDGE_AXIS` by [`check_corner_convention`] rather than
/// trusted, because a wrong convention here would silently build a different cube
/// graph from the extractor's.
fn corner_offset(corner: usize) -> [usize; 3] {
    [corner & 1, (corner >> 1) & 1, (corner >> 2) & 1]
}

/// Assert that [`corner_offset`] is the convention the crate's edge list encodes.
///
/// Each of the twelve edges joins two corners differing by one step along
/// `EDGE_AXIS[e]` and in nothing else. If that holds for all twelve, the
/// numbering is pinned.
fn check_corner_convention() {
    for edge in 0..EDGE_COUNT {
        let [lo, hi] = EDGE_CORNERS[edge];
        let (a, b) = (corner_offset(lo as usize), corner_offset(hi as usize));
        let axis = EDGE_AXIS[edge] as usize;
        for k in 0..3 {
            let expected = a[k] + usize::from(k == axis);
            assert_eq!(
                b[k], expected,
                "VOID: the rewritten corner_offset disagrees with the crate's \
                 EDGE_CORNERS/EDGE_AXIS on edge {edge}, so the derived cube graph \
                 is not the extractor's and disagreement_cells counts the wrong \
                 cells"
            );
        }
    }
}

/// The 256 cases in which the inside corners form one connected set on the cube
/// graph **and** the outside corners do — Etiene et al. §4.1's *unambiguous*
/// cells, the hypothesis of their Theorem 4.1.
///
/// The graph is built from `EDGE_CORNERS`, the crate's own edge list, so this
/// table cannot disagree with the extractor about adjacency.
fn unambiguous_cases() -> [bool; 256] {
    let mut adjacency = [0u8; CORNER_COUNT];
    for [a, b] in EDGE_CORNERS {
        adjacency[a as usize] |= 1 << b;
        adjacency[b as usize] |= 1 << a;
    }
    let mut out = [false; 256];
    for (case, unambiguous) in out.iter_mut().enumerate() {
        let inside = case as u8;
        *unambiguous =
            is_one_component(inside, &adjacency) && is_one_component(!inside, &adjacency);
    }
    out
}

/// Is `members` a single connected set on the cube graph? The empty set counts.
///
/// `x & x.wrapping_neg()` rather than `isolate_lowest_one`, which is 1.97 and the
/// MSRV is 1.89.
fn is_one_component(members: u8, adjacency: &[u8; CORNER_COUNT]) -> bool {
    if members == 0 {
        return true;
    }
    let mut reached = members & members.wrapping_neg();
    loop {
        let mut grown = reached;
        for corner in 0..CORNER_COUNT {
            if reached & (1 << corner) != 0 {
                grown |= adjacency[corner] & members;
            }
        }
        if grown == reached {
            return reached == members;
        }
        reached = grown;
    }
}

// ─── the digital object ──────────────────────────────────────────────────────

/// What the sampled box is topologically.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Domain {
    /// A finite box. Voxels outside it are empty, which is the paper's own device
    /// for keeping the digital object clear of the wall — and applies only when
    /// the object does not reach the wall, which `boundary_inside` measures.
    Bounded,
    /// The 3-torus. Opposite faces are the same period point, so the last sample
    /// on each axis is the first one and voxel indices wrap.
    Periodic,
}

/// The digitisation of `{f < 0}`, extended by one voxel per axis so that every
/// lattice vertex is the base of one in-range `2×2×2` block.
#[derive(Clone, Debug)]
struct Digital {
    /// Extended voxel counts per axis: `n + 2` bounded, `n + 1` periodic.
    extent: [usize; 3],
    /// Blocks per axis: `n + 1` bounded, `n` periodic.
    blocks: [usize; 3],
    /// Occupancy of the extended array, `x` fastest — the crate's own ordering.
    occupied: Vec<bool>,
    /// Samples of the original grid that are inside.
    inside: u64,
    /// Samples on the outermost shell of the original grid that are inside. A
    /// non-zero count is the field-side proof that the object reaches the box.
    boundary_inside: u64,
    /// Samples whose sign at index `n` differs from index `0` on some axis.
    /// Periodic domains only, and it must be zero.
    wrap_sign_mismatches: u64,
}

/// Digitise `{f < 0}` from a sample grid.
///
/// `samples` are the dimensions of `values`. A bounded domain has one voxel per
/// sample; a periodic one has `samples - 1`, because the last sample on each axis
/// is the periodic image of the first and counting it twice would double one
/// layer of the torus.
fn digitise(values: &[f64], samples: [usize; 3], domain: Domain) -> Digital {
    let voxels = match domain {
        Domain::Bounded => samples,
        Domain::Periodic => [samples[0] - 1, samples[1] - 1, samples[2] - 1],
    };
    let pad = match domain {
        Domain::Bounded => 2,
        Domain::Periodic => 1,
    };
    let extent = [voxels[0] + pad, voxels[1] + pad, voxels[2] + pad];
    let blocks = [
        voxels[0] + pad - 1,
        voxels[1] + pad - 1,
        voxels[2] + pad - 1,
    ];
    let mut occupied = vec![false; extent[0] * extent[1] * extent[2]];

    let sample =
        |x: usize, y: usize, z: usize| is_inside(values[x + samples[0] * (y + samples[1] * z)]);
    let mut inside = 0u64;
    let mut boundary_inside = 0u64;
    for z in 0..voxels[2] {
        for y in 0..voxels[1] {
            for x in 0..voxels[0] {
                if !sample(x, y, z) {
                    continue;
                }
                inside += 1;
                let shell = x == 0
                    || y == 0
                    || z == 0
                    || x + 1 == voxels[0]
                    || y + 1 == voxels[1]
                    || z + 1 == voxels[2];
                if shell && domain == Domain::Bounded {
                    boundary_inside += 1;
                }
                // Bounded: voxel `v` sits at extended index `v + 1`, leaving an
                // empty shell at each end. Periodic: at `v`, and the wrap below
                // fills index `n`.
                let base = usize::from(domain == Domain::Bounded);
                let index = (x + base) + extent[0] * ((y + base) + extent[1] * (z + base));
                occupied[index] = true;
            }
        }
    }

    let mut wrap_sign_mismatches = 0u64;
    if domain == Domain::Periodic {
        // The last layer of the extended array is the first layer of the torus.
        // Written from the occupancy rather than from `values[n]`, because
        // `sin(2*pi*N)` is `-2.4e-16` and not `0.0` and the sign of a sample near
        // the surface can flip between the two ends of the same period.
        for z in 0..extent[2] {
            for y in 0..extent[1] {
                for x in 0..extent[0] {
                    if x < voxels[0] && y < voxels[1] && z < voxels[2] {
                        continue;
                    }
                    let src = (x % voxels[0], y % voxels[1], z % voxels[2]);
                    let from = src.0 + extent[0] * (src.1 + extent[1] * src.2);
                    let index = x + extent[0] * (y + extent[1] * z);
                    occupied[index] = occupied[from];
                }
            }
        }
        for axis in 0..3 {
            let mut far = [0usize; 3];
            far[axis] = voxels[axis];
            for b in 0..samples[(axis + 1) % 3] {
                for c in 0..samples[(axis + 2) % 3] {
                    let mut near = [0usize; 3];
                    near[(axis + 1) % 3] = b;
                    near[(axis + 2) % 3] = c;
                    let mut away = near;
                    away[axis] = voxels[axis];
                    if sample(near[0], near[1], near[2]) != sample(away[0], away[1], away[2]) {
                        wrap_sign_mismatches += 1;
                    }
                }
            }
            let _ = far;
        }
    }

    Digital {
        extent,
        blocks,
        occupied,
        inside,
        boundary_inside,
        wrap_sign_mismatches,
    }
}

impl Digital {
    /// How many blocks carry each of the 256 configurations.
    ///
    /// One block per lattice vertex, and the extension above is sized so that
    /// every read is in range without a bounds test in the inner loop.
    fn configuration_counts(&self) -> [u64; 256] {
        let stride = [1usize, self.extent[0], self.extent[0] * self.extent[1]];
        let offsets: [usize; 8] = std::array::from_fn(|bit| {
            let o = BLOCK_OFFSETS[bit];
            o[0] * stride[0] + o[1] * stride[1] + o[2] * stride[2]
        });
        let mut counts = [0u64; 256];
        for z in 0..self.blocks[2] {
            for y in 0..self.blocks[1] {
                let row = y * stride[1] + z * stride[2];
                for x in 0..self.blocks[0] {
                    let base = row + x;
                    let mut configuration = 0u8;
                    for (bit, offset) in offsets.iter().enumerate() {
                        configuration |= u8::from(self.occupied[base + offset]) << bit;
                    }
                    counts[configuration as usize] += 1;
                }
            }
        }
        counts
    }

    /// Total blocks censused.
    fn block_count(&self) -> u64 {
        (self.blocks[0] * self.blocks[1] * self.blocks[2]) as u64
    }
}

/// `V − E + F − C` of the digital solid, from the block census.
///
/// The divisibility assert is the global check that the `1/8`s cancelled: a
/// weight table that mis-accounted a sharing multiplicity leaves a remainder.
fn chi_from_counts(counts: &[u64; 256], numerators: &[i64; 256]) -> i64 {
    let total: i64 = counts
        .iter()
        .zip(numerators.iter())
        .map(|(count, weight)| *count as i64 * weight)
        .sum();
    assert_eq!(
        total % WEIGHT_DENOMINATOR,
        0,
        "VOID: the weighted block census came to {total}/{WEIGHT_DENOMINATOR}, \
         which is not an integer chi -- the cell-sharing multiplicities in the \
         derived weight table do not account for every cell exactly once"
    );
    total / WEIGHT_DENOMINATOR
}

// ─── the ambiguity census ────────────────────────────────────────────────────

/// `(cells, ambiguous cells)` of a sample grid, by Etiene et al. §4.1's test.
fn ambiguity_census(values: &[f64], samples: [usize; 3], unambiguous: &[bool; 256]) -> (u64, u64) {
    let cells = [samples[0] - 1, samples[1] - 1, samples[2] - 1];
    let corner: [usize; CORNER_COUNT] = std::array::from_fn(|c| {
        let o = corner_offset(c);
        o[0] + samples[0] * (o[1] + samples[1] * o[2])
    });
    let mut ambiguous = 0u64;
    for z in 0..cells[2] {
        for y in 0..cells[1] {
            for x in 0..cells[0] {
                let base = x + samples[0] * (y + samples[1] * z);
                let mut case = 0u8;
                for (c, step) in corner.iter().enumerate() {
                    if is_inside(values[base + step]) {
                        case |= 1 << c;
                    }
                }
                if !unambiguous[case as usize] {
                    ambiguous += 1;
                }
            }
        }
    }
    ((cells[0] * cells[1] * cells[2]) as u64, ambiguous)
}

// ─── the oracle ──────────────────────────────────────────────────────────────

/// The three enumerated tables, derived once and read by every row.
#[derive(Clone, Debug)]
struct Tables {
    /// Foreground-26 weights, numerators over [`WEIGHT_DENOMINATOR`].
    weights_26: [i64; 256],
    /// Foreground-6 weights, same denominator.
    weights_6: [i64; 256],
    /// Etiene et al. §4.1's unambiguous cases.
    unambiguous: [bool; 256],
}

impl Tables {
    /// Derive all three.
    fn derive() -> Self {
        Self {
            weights_26: weight_numerators(Connectivity::Foreground26),
            weights_6: weight_numerators(Connectivity::Foreground6),
            unambiguous: unambiguous_cases(),
        }
    }

    /// How many of the 256 configurations the two connectivities weigh
    /// differently.
    fn split_configurations(&self) -> u64 {
        self.weights_26
            .iter()
            .zip(self.weights_6.iter())
            .filter(|(a, b)| a != b)
            .count() as u64
    }

    /// How many blocks of one grid landed on a configuration where the two
    /// connectivities disagree.
    fn split_blocks(&self, counts: &[u64; 256]) -> u64 {
        counts
            .iter()
            .enumerate()
            .filter(|(configuration, _)| {
                self.weights_26[*configuration] != self.weights_6[*configuration]
            })
            .map(|(_, count)| *count)
            .sum()
    }
}

/// One complete run of the digital-topology oracle.
#[derive(Clone, Copy, Debug)]
struct Oracle {
    /// `χ` of the digital solid under foreground 26.
    chi_solid_26: i64,
    /// `χ` of the same occupancy under foreground 6.
    chi_solid_6: i64,
    /// Inside samples.
    inside: u64,
    /// Inside samples on the grid's outermost shell.
    boundary_inside: u64,
    /// Blocks censused.
    blocks: u64,
    /// Blocks on a configuration the two connectivities weigh differently.
    split_blocks: u64,
    /// Cells of the sample grid.
    cells: u64,
    /// Cells failing Etiene et al. §4.1's unambiguity test.
    ambiguous_cells: u64,
    /// Periodic sign mismatches across the period; must be zero.
    wrap_sign_mismatches: u64,
}

/// Sample the field on the grid, digitise `{f < 0}`, and read χ off the block
/// census — the whole cost a test would pay, timed as one thing.
fn run_oracle<S: Sdf<Scalar = f64>>(
    field: &S,
    samples: [usize; 3],
    origin: [f64; 3],
    cell_size: f64,
    domain: Domain,
    tables: &Tables,
) -> Oracle {
    let mut values = Vec::with_capacity(samples[0] * samples[1] * samples[2]);
    for z in 0..samples[2] {
        for y in 0..samples[1] {
            for x in 0..samples[0] {
                values.push(field.sample([
                    origin[0] + cell_size * x as f64,
                    origin[1] + cell_size * y as f64,
                    origin[2] + cell_size * z as f64,
                ]));
            }
        }
    }
    let digital = digitise(&values, samples, domain);
    let counts = digital.configuration_counts();
    let (cells, ambiguous_cells) = ambiguity_census(&values, samples, &tables.unambiguous);
    Oracle {
        chi_solid_26: chi_from_counts(&counts, &tables.weights_26),
        chi_solid_6: chi_from_counts(&counts, &tables.weights_6),
        inside: digital.inside,
        boundary_inside: digital.boundary_inside,
        blocks: digital.block_count(),
        split_blocks: tables.split_blocks(&counts),
        cells,
        ambiguous_cells,
        wrap_sign_mismatches: digital.wrap_sign_mismatches,
    }
}

/// Repeats of the timed oracle. Five is `P-163`'s floor and the least that has a
/// median distinct from its extremes.
const REPEATS: usize = 5;

/// `(oracle, median ms, min ms, max ms)` over [`REPEATS`] runs after one warm-up.
fn timed_oracle<S: Sdf<Scalar = f64>>(
    field: &S,
    samples: [usize; 3],
    origin: [f64; 3],
    cell_size: f64,
    domain: Domain,
    tables: &Tables,
) -> (Oracle, f64, f64, f64) {
    let mut oracle = run_oracle(field, samples, origin, cell_size, domain, tables);
    let mut timings = Vec::with_capacity(REPEATS);
    for _ in 0..REPEATS {
        let started = Instant::now();
        oracle = run_oracle(field, samples, origin, cell_size, domain, tables);
        timings.push(started.elapsed().as_secs_f64() * 1e3);
    }
    timings.sort_unstable_by(f64::total_cmp);
    (
        oracle,
        timings[REPEATS / 2],
        timings[0],
        timings[REPEATS - 1],
    )
}

// ─── the extractors under test ───────────────────────────────────────────────

/// One Marching Cubes configuration, named as `golden.rs`'s roster names it.
#[derive(Clone, Copy, Debug)]
struct Arm {
    /// The `extractor` column.
    name: &'static str,
    /// Face-ambiguity rule.
    face: FaceAmbiguity,
    /// Interior-ambiguity rule.
    interior: InteriorAmbiguity,
}

/// The crate's default, and its ambiguity-resolving configuration.
///
/// The second is the one that claims to reproduce the trilinear interpolant's
/// topology, which is precisely the claim Etiene et al.'s oracle exists to test.
const ARMS: [Arm; 2] = [
    Arm {
        name: "marching_cubes",
        face: FaceAmbiguity::Separate,
        interior: InteriorAmbiguity::Ignore,
    },
    Arm {
        name: "marching_cubes+trilinear",
        face: FaceAmbiguity::AsymptoticDecider,
        interior: InteriorAmbiguity::Trilinear,
    },
];

/// The two grids every reference field is read on.
///
/// Word-boundary sample counts, one below and one above `thin_plate`'s
/// construction cell size of `4/64`.
const RESOLUTIONS: [u32; 2] = [33, 65];

/// Voxels per period for the nodal arm. **Odd on purpose** — see the header.
const VOXELS_PER_PERIOD: u32 = 33;

/// Periods per axis for the nodal arm: `P-142` C1's registered range.
const NODAL_PERIODS: [u32; 3] = [1, 2, 3];

// ─── rows ────────────────────────────────────────────────────────────────────

/// One CSV row: one field at one resolution read by one extractor.
#[derive(Clone, Debug)]
struct Row {
    /// `calibration`, `new_ground_truth`, `out_of_scope` or `p142_cross_check`.
    arm: &'static str,
    /// The `field` column.
    field: &'static str,
    /// Samples per axis.
    resolution: u32,
    /// Sample spacing.
    cell_size: f64,
    /// Is this a registered or assigned calibration point?
    is_control: bool,
    /// Does the digital-topology half of the method apply at all?
    in_scope: bool,
    /// `expected_euler()`, or `None` where the crate declares it unknown.
    chi_expected: Option<i64>,
    /// The analytic prediction `P-142` established, on the nodal arm only.
    chi_predicted_p142: Option<i64>,
    /// Periods per axis on the nodal arm.
    periods: Option<u32>,
    /// Seam vertex pairs the periodic wrap identified.
    seam_pairs: Option<u64>,
    /// What the oracle produced.
    oracle: Oracle,
    /// Median, min and max oracle wall clock in milliseconds.
    method_ms: (f64, f64, f64),
    /// The extractor.
    extractor: &'static str,
    /// Extraction wall clock in milliseconds.
    extract_ms: f64,
    /// `V − E + F` of the extracted mesh, from `common::tpms::euler`.
    chi_extracted: i64,
    /// The same from the crate's `validate_indexed`, asserted equal.
    chi_extracted_validate: i64,
    /// Mesh readings the scope decision is made from.
    boundary_edges: u64,
    /// Edges shared by three or more triangles.
    nonmanifold_edges: u64,
    /// Mesh size.
    mesh_vertices: usize,
    /// Mesh size.
    mesh_triangles: usize,
}

impl Row {
    /// The surface χ the oracle asserts, `2 · χ(solid)`, where it applies.
    fn chi_surface_26(&self) -> Option<i64> {
        self.in_scope.then_some(2 * self.oracle.chi_solid_26)
    }

    /// The same under foreground 6.
    fn chi_surface_6(&self) -> Option<i64> {
        self.in_scope.then_some(2 * self.oracle.chi_solid_6)
    }

    /// Does the oracle agree with the mesh? `None` where the oracle does not
    /// apply, which is recorded as `not_applicable` rather than as a false.
    fn agreement(&self) -> Option<bool> {
        self.chi_surface_26().map(|chi| chi == self.chi_extracted)
    }
}

/// Format an optional integer for the CSV: `none` where the value does not exist.
fn opt_i64(value: Option<i64>) -> String {
    value.map_or_else(|| String::from("none"), |v| v.to_string())
}

/// Format an optional unsigned integer the same way.
fn opt_u64(value: Option<u64>) -> String {
    value.map_or_else(|| String::from("none"), |v| v.to_string())
}

/// Extract one mesh and read its topology two ways.
///
/// `wrap` is the periodic identification, applied before either reading; on a
/// bounded domain it is `None` and no seam is touched. The two readers must agree
/// about `V − E + F` at one weld tolerance, or the instrument failed rather than
/// the hypothesis.
fn extract_and_read<S: Sdf<Scalar = f64>>(
    field: &S,
    arm: Arm,
    shape: &impl Shape3,
    origin: [f64; 3],
    cell_size: f64,
    wrap: Option<([f64; 3], [f64; 3])>,
) -> (tpms::EulerCount, i64, Option<u64>, f64, MeshBuffer<f64>) {
    let mut extractor = MarchingCubes::<f64>::new();
    extractor.set_face_ambiguity(arm.face);
    extractor.set_interior_ambiguity(arm.interior);
    let mut mesh = MeshBuffer::<f64>::new();
    let started = Instant::now();
    extractor
        .extract(field, shape, origin, cell_size, &mut mesh)
        .expect("every grid here holds at least two samples on every axis");
    let extract_ms = started.elapsed().as_secs_f64() * 1e3;

    let tol = isomesh::weld::epsilon_for(cell_size);
    let seam_pairs = wrap.map(|(lo, hi)| tpms::wrap_seams(&mut mesh, lo, hi, tol));
    let counted = tpms::euler(&mesh.positions, &mesh.indices, tol);
    let cfg = ValidateConfig::from_cell_size(cell_size)
        .expect("every cell size here is finite and positive");
    let report = validate_indexed(&mesh.positions, &mesh.indices, &cfg);
    assert_eq!(
        counted.chi, report.euler_characteristic,
        "VOID: common::tpms::euler and validate_indexed disagree about V-E+F at \
         one weld tolerance on {}, so the instrument is what failed",
        arm.name
    );
    (
        counted,
        report.euler_characteristic,
        seam_pairs,
        extract_ms,
        mesh,
    )
}

/// Every row for one reference field.
fn reference_rows<F>(field: &F, tables: &Tables) -> Vec<Row>
where
    F: ReferenceField + Sdf<Scalar = f64>,
{
    let mut rows = Vec::with_capacity(RESOLUTIONS.len() * ARMS.len());
    let chi_expected = field.expected_euler();
    // A field the crate declares closed is one whose level set has no boundary,
    // which is Etiene et al. §4.1's own precondition. `fbm_terrain` is the only
    // reference field that fails it, and both `boundary_inside` and
    // `boundary_edges` are recorded so the exclusion is measured twice.
    let in_scope = field.closed_in_domain();
    let arm = match (in_scope, chi_expected) {
        (true, Some(_)) => "calibration",
        (true, None) => "new_ground_truth",
        (false, _) => "out_of_scope",
    };
    for samples in RESOLUTIONS {
        let (shape, origin, cell_size) = common::grid::<f64, _>(field, samples);
        let dims = [samples as usize; 3];
        let (oracle, median, low, high) =
            timed_oracle(field, dims, origin, cell_size, Domain::Bounded, tables);
        for extract_arm in ARMS {
            let (counted, validate_chi, _, extract_ms, mesh) =
                extract_and_read(field, extract_arm, &shape, origin, cell_size, None);
            rows.push(Row {
                arm,
                field: F::NAME,
                resolution: samples,
                cell_size,
                is_control: chi_expected.is_some(),
                in_scope,
                chi_expected,
                chi_predicted_p142: None,
                periods: None,
                seam_pairs: None,
                oracle,
                method_ms: (median, low, high),
                extractor: extract_arm.name,
                extract_ms,
                chi_extracted: counted.chi,
                chi_extracted_validate: validate_chi,
                boundary_edges: counted.boundary_edges,
                nonmanifold_edges: counted.non_manifold_edges,
                mesh_vertices: mesh.vertex_count(),
                mesh_triangles: mesh.triangle_count(),
            });
        }
    }
    rows
}

/// Every row for the nodal gyroid on the 3-torus — C1's cross-check.
fn nodal_rows(tables: &Tables) -> Vec<Row> {
    let mut rows = Vec::with_capacity(NODAL_PERIODS.len() * ARMS.len());
    for periods in NODAL_PERIODS {
        let field = NodalTpms::new(Tpms::Gyroid, periods);
        let (lo, hi) = field.domain();
        let (shape, origin, cell_size) = field.periodic_grid(VOXELS_PER_PERIOD);
        let samples = VOXELS_PER_PERIOD * periods + 1;
        let dims = [samples as usize; 3];
        let (oracle, median, low, high) =
            timed_oracle(&field, dims, origin, cell_size, Domain::Periodic, tables);
        for extract_arm in ARMS {
            let (counted, validate_chi, seam_pairs, extract_ms, mesh) = extract_and_read(
                &field,
                extract_arm,
                &shape,
                origin,
                cell_size,
                Some((lo, hi)),
            );
            rows.push(Row {
                arm: "p142_cross_check",
                field: "gyroid_nodal",
                resolution: samples,
                cell_size,
                is_control: false,
                in_scope: true,
                chi_expected: None,
                chi_predicted_p142: Some(field.chi_predicted()),
                periods: Some(periods),
                seam_pairs,
                oracle,
                method_ms: (median, low, high),
                extractor: extract_arm.name,
                extract_ms,
                chi_extracted: counted.chi,
                chi_extracted_validate: validate_chi,
                boundary_edges: counted.boundary_edges,
                nonmanifold_edges: counted.non_manifold_edges,
                mesh_vertices: mesh.vertex_count(),
                mesh_triangles: mesh.triangle_count(),
            });
        }
    }
    rows
}

/// `χ` of a hand-built occupancy, for the micro-fixtures the weight tables are
/// calibrated on.
fn chi_of_fixture(occupied: &[[usize; 3]], extent: [usize; 3], tables: &Tables) -> (i64, i64) {
    let mut values = vec![1.0f64; extent[0] * extent[1] * extent[2]];
    for v in occupied {
        values[v[0] + extent[0] * (v[1] + extent[1] * v[2])] = -1.0;
    }
    let digital = digitise(&values, extent, Domain::Bounded);
    let counts = digital.configuration_counts();
    (
        chi_from_counts(&counts, &tables.weights_26),
        chi_from_counts(&counts, &tables.weights_6),
    )
}

fn main() {
    if !std::env::args().any(|arg| arg == "--bench") {
        return;
    }
    let prereg = isomesh::experiment!("P-145");

    common::experiment::run(prereg, |run| {
        let started = Instant::now();
        check_corner_convention();
        let tables = Tables::derive();

        let mut rows: Vec<Row> = Vec::new();
        isomesh::for_each_reference_field!(f64, |_name, field| {
            rows.extend(reference_rows(&field, &tables));
        });
        rows.extend(nodal_rows(&tables));
        println!(
            "-- {} rows over {} fields x {} extractors in {:.1} s",
            rows.len(),
            rows.len() / ARMS.len(),
            ARMS.len(),
            started.elapsed().as_secs_f64()
        );

        // ── vacuity controls, before any record ─────────────────────────────

        // The weight tables are derived, and the derivation is checkable.
        for (connectivity, weights) in [
            (Connectivity::Foreground26, &tables.weights_26),
            (Connectivity::Foreground6, &tables.weights_6),
        ] {
            let single = weights[1 << 7];
            assert!(
                weights[0] == 0 && weights[255] == 0 && single == 1,
                "VOID: {}'s derived weights read empty={}, full={}, single voxel={}/{} \
                 against the required 0, 0 and 1/{}: an empty block and an interior \
                 block must both contribute nothing (chi is a boundary phenomenon) \
                 and one isolated voxel must be an eighth in each of its eight \
                 corner blocks. A table failing these is transcribed, not derived.",
                connectivity.name(),
                weights[0],
                weights[255],
                single,
                WEIGHT_DENOMINATOR,
                WEIGHT_DENOMINATOR
            );
        }

        // The two connectivities are a real distinction, in the table and on an
        // actual occupancy.
        let split = tables.split_configurations();
        let (corner_26, corner_6) = chi_of_fixture(&[[1, 1, 1], [2, 2, 2]], [4, 4, 4], &tables);
        let (single_26, single_6) = chi_of_fixture(&[[1, 1, 1]], [3, 3, 3], &tables);
        assert!(
            split > 0 && corner_26 == 1 && corner_6 == 2 && single_26 == 1 && single_6 == 1,
            "VOID: the connectivity choice changes nothing measurable -- {split} of 256 \
             configurations differ, the corner-touching pair reads chi_26={corner_26} \
             chi_6={corner_6} against the required 1 and 2, and a single voxel reads \
             {single_26}/{single_6} against 1 and 1. Without all four the `connectivity` \
             column names a choice with no consequence."
        );

        // The registered control, and the assigned second one.
        for (name, want) in [("sphere", 2i64), ("torus", 0i64)] {
            let mut seen = 0;
            for row in rows.iter().filter(|r| r.field == name) {
                seen += 1;
                assert!(
                    row.chi_surface_26() == Some(want) && row.chi_surface_6() == Some(want),
                    "VOID: the oracle reads chi {} / {} on {name} at {}^3 instead of \
                     {want}, so it is not calibrated and no number it reports on a \
                     field with no known chi means anything",
                    opt_i64(row.chi_surface_26()),
                    opt_i64(row.chi_surface_6()),
                    row.resolution
                );
            }
            assert!(
                seen == RESOLUTIONS.len() * ARMS.len(),
                "VOID: {seen} {name} rows, so the registered calibration ran on a \
                 different sweep from the one being reported"
            );
        }

        // A calibrated constant is still a constant.
        let mut distinct: Vec<i64> = rows.iter().filter_map(Row::chi_surface_26).collect();
        distinct.sort_unstable();
        distinct.dedup();
        assert!(
            distinct.len() >= 3,
            "VOID: the oracle produced only {} distinct chi over the whole sweep \
             ({distinct:?}), so agreeing with 2 on the sphere and 0 on the torus is \
             consistent with an instrument that cannot tell fields apart (M-44)",
            distinct.len()
        );

        // Neither empty nor full, and the mesh exists.
        for row in &rows {
            let samples = u64::from(row.resolution).pow(3);
            assert!(
                row.oracle.inside > 0 && row.oracle.inside < samples,
                "VOID: {} at {}^3 has {} of {samples} samples inside, so its digital \
                 object is empty or fills the grid and every chi it reports is zero \
                 by vacuity",
                row.field,
                row.resolution,
                row.oracle.inside
            );
            assert!(
                row.mesh_triangles > 0,
                "VOID: {} at {}^3 with {} produced no triangle, so chi_extracted is \
                 0 by vacuity rather than by measurement",
                row.field,
                row.resolution,
                row.extractor
            );
        }

        // Theorem 4.1's hypothesis was actually tested.
        let ambiguous_total: u64 = rows.iter().map(|r| r.oracle.ambiguous_cells).sum();
        assert!(
            ambiguous_total > 0,
            "VOID: not one cell in the whole sweep is ambiguous, so disagreement_cells \
             is a column of zeros that could not have been non-zero and Etiene's \
             Theorem 4.1 hypothesis was never exercised (M-44)"
        );

        // The derived ambiguity table must cover the crate's own face ambiguity.
        for case in 0..256usize {
            if AMBIGUOUS_FACES[case] != 0 {
                assert!(
                    !tables.unambiguous[case],
                    "VOID: case {case:#04x} carries an ambiguous face by the crate's \
                     AMBIGUOUS_FACES yet the derived table calls it unambiguous, so \
                     disagreement_cells undercounts exactly the cells the extractor \
                     has a free choice in"
                );
            }
        }

        // The exclusion is measured from both sides, and they agree.
        for row in &rows {
            let field_side = row.oracle.boundary_inside > 0;
            let mesh_side = row.boundary_edges > 0;
            assert!(
                field_side == !row.in_scope && mesh_side == !row.in_scope,
                "VOID: {} at {}^3 with {} says boundary_inside={} boundary_edges={} \
                 against in_scope={}: the field-side and mesh-side readings of \
                 'closed' must agree with each other and with the scope decision, \
                 or applicable_fields is an opinion rather than a measurement",
                row.field,
                row.resolution,
                row.extractor,
                row.oracle.boundary_inside,
                row.boundary_edges,
                row.in_scope
            );
        }

        // The periodic identification is exact.
        for row in rows.iter().filter(|r| r.periods.is_some()) {
            assert!(
                row.oracle.wrap_sign_mismatches == 0,
                "VOID: {} signs disagree between sample 0 and sample n on \
                 gyroid_nodal N={}, so the torus occupancy is glued to the wrong \
                 layer and its chi is not the chi of any closed surface",
                row.oracle.wrap_sign_mismatches,
                opt_u64(row.periods.map(u64::from))
            );
            assert!(
                row.seam_pairs.unwrap_or(0) > 0,
                "VOID: the periodic wrap identified no seam pair on gyroid_nodal \
                 N={}, so the mesh reading is the non-wrapped one and the \
                 cross-check compares the oracle against an open surface",
                opt_u64(row.periods.map(u64::from))
            );
        }

        // The DOI reported is the DOI registered.
        assert!(
            prereg.hypothesis.contains(DOI),
            "VOID: the registration's hypothesis does not contain {DOI}, so the \
             doi_verified column reports the resolution of a retyped DOI rather \
             than of the registered claim"
        );

        // ── the clause verdicts, global ─────────────────────────────────────

        // C1(i): at least two fields the crate declares unknown now have one.
        let mut new_fields: Vec<&'static str> = rows
            .iter()
            .filter(|r| r.arm == "new_ground_truth" && r.chi_surface_26().is_some())
            .map(|r| r.field)
            .collect();
        new_fields.sort_unstable();
        new_fields.dedup();

        // C1(ii): the nodal arm must reproduce P-142's analytic -8N^3.
        let cross_check_rows: Vec<&Row> = rows.iter().filter(|r| r.periods.is_some()).collect();
        let cross_check = !cross_check_rows.is_empty()
            && cross_check_rows
                .iter()
                .all(|r| r.chi_surface_26() == r.chi_predicted_p142);
        let c1 = new_fields.len() >= 2 && cross_check;

        // C2: the cost is reported and the applicable set is non-empty. Falsified
        // by inapplicability to every field, and by nothing else.
        let mut applicable: Vec<&'static str> = rows
            .iter()
            .filter(|r| r.in_scope)
            .map(|r| r.field)
            .collect();
        applicable.sort_unstable();
        applicable.dedup();
        let costed = rows
            .iter()
            .all(|r| r.method_ms.0.is_finite() && r.method_ms.0 > 0.0);
        let c2 = !applicable.is_empty() && costed;
        let applicable_fields = applicable.join("|");

        println!(
            "-- C1 {}: {} new ground truths ({}), P-142 cross-check {}",
            if c1 { "HELD" } else { "FALSIFIED" },
            new_fields.len(),
            new_fields.join("|"),
            if cross_check { "agrees" } else { "DISAGREES" }
        );
        println!(
            "-- C2 {}: applicable on {applicable_fields}",
            if c2 { "HELD" } else { "FALSIFIED" }
        );

        for row in &rows {
            let (median, low, high) = row.method_ms;
            run.record(&[
                ("field", row.field.to_string()),
                ("resolution", row.resolution.to_string()),
                ("chi_stratified_morse", opt_i64(row.chi_surface_26())),
                ("chi_extracted", row.chi_extracted.to_string()),
                (
                    "agreement",
                    row.agreement()
                        .map_or_else(|| String::from("not_applicable"), |a| a.to_string()),
                ),
                ("disagreement_cells", row.oracle.ambiguous_cells.to_string()),
                ("method_cost_ms", format!("{median:.4}")),
                ("applicable_fields", applicable_fields.clone()),
                ("doi_verified", true.to_string()),
                ("c1_holds", c1.to_string()),
                ("c2_holds", c2.to_string()),
                // ── extras (M-273) ──
                ("arm", row.arm.to_string()),
                ("blocks", row.oracle.blocks.to_string()),
                ("boundary_edges", row.boundary_edges.to_string()),
                ("boundary_inside", row.oracle.boundary_inside.to_string()),
                ("cell_size", format!("{:.9}", row.cell_size)),
                ("cells", row.oracle.cells.to_string()),
                ("chi_expected", opt_i64(row.chi_expected)),
                (
                    "chi_extracted_validate",
                    row.chi_extracted_validate.to_string(),
                ),
                ("chi_predicted_p142", opt_i64(row.chi_predicted_p142)),
                ("chi_solid_26", row.oracle.chi_solid_26.to_string()),
                ("chi_solid_6", row.oracle.chi_solid_6.to_string()),
                ("chi_surface_6", opt_i64(row.chi_surface_6())),
                (
                    "connectivity",
                    Connectivity::Foreground26.name().to_string(),
                ),
                (
                    "connectivity_split_blocks",
                    row.oracle.split_blocks.to_string(),
                ),
                ("connectivity_split_configurations", split.to_string()),
                (
                    "cost_ratio",
                    format!("{:.4}", median / row.extract_ms.max(f64::MIN_POSITIVE)),
                ),
                ("doi", DOI.to_string()),
                ("doi_handles", DOI_HANDLES.to_string()),
                ("doi_resolved_title", DOI_RESOLVED_TITLE.to_string()),
                ("doi_venue", DOI_VENUE.to_string()),
                ("extract_cost_ms", format!("{:.4}", row.extract_ms)),
                ("extractor", row.extractor.to_string()),
                ("inside_samples", row.oracle.inside.to_string()),
                ("is_control", row.is_control.to_string()),
                (
                    "matches_chi_expected",
                    match (row.chi_expected, row.chi_surface_26()) {
                        (Some(want), Some(got)) => (want == got).to_string(),
                        _ => String::from("none"),
                    },
                ),
                ("mesh_triangles", row.mesh_triangles.to_string()),
                ("mesh_vertices", row.mesh_vertices.to_string()),
                ("method_cost_ms_max", format!("{high:.4}")),
                ("method_cost_ms_min", format!("{low:.4}")),
                ("method_repeats", REPEATS.to_string()),
                ("nonmanifold_edges", row.nonmanifold_edges.to_string()),
                (
                    "oracle_method",
                    String::from("digital_topology_etiene_2012_sec_4_1"),
                ),
                (
                    "oracle_scope",
                    String::from(if row.in_scope {
                        "in_scope"
                    } else {
                        "out_of_scope_level_set_has_boundary"
                    }),
                ),
                ("periods", opt_u64(row.periods.map(u64::from))),
                ("seam_pairs", opt_u64(row.seam_pairs)),
                (
                    "wrap_sign_mismatches",
                    row.oracle.wrap_sign_mismatches.to_string(),
                ),
            ]);
        }
    });
}
