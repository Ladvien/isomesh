//! **P-102 — ✗43's withdrawn prevalence sweep, rebuilt as a bench.**
//!
//! Ticket: R-102. Pre-registered before this harness existed.
//!
//! ```bash
//! cargo bench --bench experiment_p102
//! ```
//!
//! Writes `docs/experiments/p-102.csv`.
//!
//! # What was missing
//!
//! `✗43` reported *"2 of 8,064 not closed before the fix, 0 of 8,064 after"* from
//! a deterministic sweep of **1,152 single-plane caps (16 θ × 8 cos φ × 9
//! offsets over the generator's own `-0.8..=0.8`) at each of 6³–12³**, and that
//! sweep *"exists in no bench, no test and no CSV"* — the entry says so itself,
//! withdraws the rate, and declines to invent a `P-` id for re-adding it. This
//! is that id. Nothing here changes the crate: the fix is already in, and the
//! only question is whether its generalisation has evidence behind it.
//!
//! # SHARE
//!
//! **No clause here is a ratio of times, so `✗51`'s `1/(1 − share/factor)`
//! bar does not apply.** Both rate clauses are counts over an *enumerated*
//! population, so their denominators are exact by construction: 1,152 per size
//! and 8,064 over the seven sizes. What is not exact by construction is
//! C2's reachable share, and that is the number this harness has to name before
//! it runs:
//!
//! - **C1's share is 8,064 of 8,064.** Every configuration is meshed and
//!   validated, so a single unclosed mesh anywhere moves the clause. The zero is
//!   only worth having if the meshes are real, so `triangles` and `empty_meshes`
//!   are columns and `empty_meshes` is asserted zero.
//! - **C2's share is `fan_configurations`, not 1,152.** The pre-fix defect can
//!   only exist in a configuration that has a cell fanning **two or more** rings
//!   from one apex; every other configuration is bit-identical before and after
//!   the fix, because a single fanned ring took the body saddle either way. So
//!   `unclosed_pre_fix ≤ fan_configurations`, and *"exactly 2 at 6³"* is
//!   arithmetically unreachable unless `fan_configurations ≥ 2` there. That
//!   share is a column, and it is the one number that decides whether C2 could
//!   have come out at all.
//! - **C3 moves nothing and is a statement about a fixture**, so it has no share.
//!
//! # The four factors, and what the entry under-determines
//!
//! `✗43` pins the factors themselves — `16 θ × 8 cos φ × 9 offsets` over
//! `property::unit_vector`'s `(0..τ, −1..=1)` and `convex_body`'s `−0.8..=0.8`,
//! at 7 sizes — so they are recovered rather than invented, and they are an
//! **orientation and offset** enumeration rather than a case/ring one. What the
//! entry does **not** pin is the discretisation of each interval: 16 values of θ
//! over a half-open interval could be `i·τ/16` or `(i+½)·τ/16`, and 8 values of
//! `cos φ` over a **closed** `[−1, 1]` could be the endpoints (`−1 + 2j/7`), the
//! cell centres (`−1 + (2j+1)/8`) or the half-open ladder (`−1 + 2j/8`). Only
//! the offsets are forced: 9 values over a closed `−0.8..=0.8` is the step-0.2
//! ladder, endpoints included.
//!
//! The harness therefore **prints the pre-fix count at 6³ under all six
//! readings** before it sweeps, and the sweep runs under [`READING`], chosen —
//! per the ticket — as the reading that reproduces C2's registered count of 2.
//! The calibration table is the evidence for that choice and is printed on every
//! run, so the choice can be re-checked rather than believed.
//!
//! `✗43`'s own reduced counterexample is **not** on any of these grids
//! (its `θ ≈ 1.3301`, `cos φ ≈ 0.0624`), which is consistent with the entry:
//! that plane came from proptest and the sweep was a separate prevalence scan.
//!
//! # The pre-fix fan, bench-local
//!
//! `✗43`'s defect was `Contours::fan` naming one shared `INTERIOR` apex for
//! every ring of a cell. The fix names ring `r`'s apex `INTERIOR + r` and places
//! it at the body saddle when exactly one ring is fanned and at that ring's own
//! centroid when several are. So the pre-fix mesh is recovered from the post-fix
//! one by **identifying a cell's apexes into a single vertex and putting it back
//! on the saddle** — the exact operation the fix reversed, in both topology and
//! geometry. This is `P-63`'s control, and here it is calibrated against
//! published ground truth rather than argued: `✗43` states the pre-fix report for
//! its own cell at 6³ (`V 27, E 72, F 48, χ 3, components 2, non-manifold
//! vertices 1, degenerate triangles 1`, `is_closed()` false) and this harness
//! reproduces that report, number for number, before it sweeps.
//!
//! Which cells those are is not guessed: the crate's own public trilinear
//! machinery is re-run per cell — `table::AMBIGUOUS_FACES`,
//! `ambiguity::joined_mask`, `Contours::of` and `BodySaddles::of` — and the
//! number of apexes it predicts is checked against the number of cell-interior
//! vertices found in the mesh. A disagreement is `classification_mismatches`,
//! and it is asserted zero on every live arm rather than assumed.
//!
//! One documented difference from the pre-fix code: pre-fix also created an apex
//! in a cell whose interior vertex existed but whose every ring was a triangle —
//! *"an unreferenced vertex in the output"*. That is reproduced too
//! (`unreferenced_pre_fix`), and it cannot affect `is_closed()`, because χ is
//! computed from **referenced** vertices.
//!
//! # The magnitude arms, and why `±1` is void here
//!
//! `P-63`'s protocol, inherited: the sign pattern is the geometry's, and the
//! **magnitudes** are a sampled parameter reported per seed and never pooled.
//! `exact` is `✗43`'s fixture verbatim — the analytic cap. `unit` and
//! `generic/s0..s3` replace each sample's magnitude while keeping its sign,
//! through the trilinear interpolant of the re-magnituded lattice, which is the
//! same object Marching Cubes reads: with `crossing_refinement` at its default
//! of zero the extractor touches the field **only** at lattice points, so the
//! interpolant of the cap's own samples must reproduce the cap's mesh bit for
//! bit. That is asserted over all 1,152 configurations at both ends of the size
//! range before anything else runs, which is what makes the generic arms the
//! same experiment rather than a different one.
//!
//! `M-374` found what `±1` does: every corner at unit magnitude makes the
//! trilinear's saddles symmetric, the strict `0 < x < 1` test behind
//! `has_inner_hexagon` and `interior_vertex` rejects, and the interior rule
//! never fires. Such an arm's zero is not a pass — it is the interior-rule-off
//! arm wearing a different name — so it is reported **VOID**, and C3 is the
//! clause that says so.

#![allow(
    clippy::cast_possible_truncation,
    clippy::cast_precision_loss,
    clippy::float_cmp,
    clippy::too_many_lines
)]

mod common;

use std::time::Instant;

use isomesh::fields::Sphere;
use isomesh::marching_cubes::ambiguity::joined_mask;
use isomesh::marching_cubes::table::AMBIGUOUS_FACES;
use isomesh::marching_cubes::trilinear::{BodySaddles, Contours};
use isomesh::marching_cubes::{FaceAmbiguity, InteriorAmbiguity, MarchingCubes};
use isomesh::validate::{MeshReport, ValidateConfig, validate_indexed};
use isomesh::{MeshBuffer, Real, RuntimeShape3, Sdf};

// ─── the population ────────────────────────────────────────────────────────

/// `property::DOMAIN`: the half-extent every generated field lives inside.
const DOMAIN: f64 = 2.0;
/// `convex_body`'s own bounding sphere radius.
const BOUND: f64 = 1.5;
/// `convex_body`'s own offset limit: the range is `-0.8..=0.8`.
const OFFSET_LIMIT: f64 = 0.8;

/// θ values. `✗43`'s first factor.
const THETAS: u32 = 16;
/// `cos φ` values. `✗43`'s second factor.
const COS_PHIS: u32 = 8;
/// Plane offsets. `✗43`'s third factor.
const OFFSETS: u32 = 9;
/// 1,152 caps.
const CONFIGURATIONS: u32 = THETAS * COS_PHIS * OFFSETS;
/// `✗43`'s fourth factor: 6³ to 12³, which is where its two failures were.
const SIZES: [u32; 7] = [6, 7, 8, 9, 10, 11, 12];

/// The seeds the generic arms are drawn from. `P-63`'s, unchanged, so the two
/// experiments' magnitude arms are the same four draws.
const SEEDS: [u64; 4] = [0x0000_2026, 0x0005_EED1, 0x00C0_FFEE, 0xDEAD_BEEF];

/// How θ is laid out over `0..τ`.
#[derive(Clone, Copy, PartialEq, Eq)]
enum Theta {
    /// `i·τ/16`. The half-open interval's own ladder, starting at zero.
    Ladder,
    /// `(i+½)·τ/16`. Cell centres.
    Centres,
}

/// How `cos φ` is laid out over the closed `-1..=1`.
#[derive(Clone, Copy, PartialEq, Eq)]
enum CosPhi {
    /// `-1 + 2j/7`. Both endpoints included, so two normals are axis-aligned.
    Endpoints,
    /// `-1 + (2j+1)/8`. Cell centres, so no normal is axis-aligned.
    Centres,
    /// `-1 + 2j/8`. The lower endpoint only.
    Ladder,
}

/// One reading of `✗43`'s *"16 θ × 8 cos φ × 9 offsets"*.
#[derive(Clone, Copy, PartialEq, Eq)]
struct Reading {
    theta: Theta,
    cos_phi: CosPhi,
}

impl Reading {
    /// A short label with no comma in it, since the CSV writer does not quote.
    fn label(self) -> &'static str {
        match (self.theta, self.cos_phi) {
            (Theta::Ladder, CosPhi::Endpoints) => "theta_ladder+cosphi_endpoints",
            (Theta::Ladder, CosPhi::Centres) => "theta_ladder+cosphi_centres",
            (Theta::Ladder, CosPhi::Ladder) => "theta_ladder+cosphi_ladder",
            (Theta::Centres, CosPhi::Endpoints) => "theta_centres+cosphi_endpoints",
            (Theta::Centres, CosPhi::Centres) => "theta_centres+cosphi_centres",
            (Theta::Centres, CosPhi::Ladder) => "theta_centres+cosphi_ladder",
        }
    }
}

/// Every reading the entry's wording admits. Printed at 6³ before the sweep, so
/// the choice below is visible rather than asserted.
const READINGS: [Reading; 6] = [
    Reading {
        theta: Theta::Ladder,
        cos_phi: CosPhi::Endpoints,
    },
    Reading {
        theta: Theta::Ladder,
        cos_phi: CosPhi::Centres,
    },
    Reading {
        theta: Theta::Ladder,
        cos_phi: CosPhi::Ladder,
    },
    Reading {
        theta: Theta::Centres,
        cos_phi: CosPhi::Endpoints,
    },
    Reading {
        theta: Theta::Centres,
        cos_phi: CosPhi::Centres,
    },
    Reading {
        theta: Theta::Centres,
        cos_phi: CosPhi::Ladder,
    },
];

/// The reading the sweep runs under.
///
/// Chosen by the calibration table this harness prints: the ticket's instruction
/// is to take the reading that reproduces C2's registered count of 2 at 6³, and
/// to say that this is what was done. It is said here and in the report.
const READING: Reading = Reading {
    theta: Theta::Ladder,
    cos_phi: CosPhi::Endpoints,
};

/// One cap: an index triple into the three ladders.
#[derive(Clone, Copy)]
struct Config {
    theta: u32,
    cos_phi: u32,
    offset: u32,
}

impl Config {
    /// Every configuration, in a fixed order, so `unclosed_*_where` names the
    /// same cap on every machine.
    fn all() -> Vec<Self> {
        let mut out = Vec::with_capacity(CONFIGURATIONS as usize);
        for theta in 0..THETAS {
            for cos_phi in 0..COS_PHIS {
                for offset in 0..OFFSETS {
                    out.push(Self {
                        theta,
                        cos_phi,
                        offset,
                    });
                }
            }
        }
        out
    }

    /// A comma-free identifier for one cap.
    fn label(self) -> String {
        format!("t{}c{}o{}", self.theta, self.cos_phi, self.offset)
    }

    /// The plane this configuration names, under one reading.
    ///
    /// `unit_vector`'s construction exactly — `cos φ` sampled rather than `φ`,
    /// and `Real::cos`/`Real::sin` rather than `std`'s, because those are the
    /// functions the generator calls and they need not agree in the last bit.
    fn plane(self, reading: Reading) -> ([f64; 3], f64) {
        let theta = match reading.theta {
            Theta::Ladder => core::f64::consts::TAU * f64::from(self.theta) / f64::from(THETAS),
            Theta::Centres => {
                core::f64::consts::TAU * (f64::from(self.theta) + 0.5) / f64::from(THETAS)
            }
        };
        let cos_phi = match reading.cos_phi {
            CosPhi::Endpoints => -1.0 + 2.0 * f64::from(self.cos_phi) / f64::from(COS_PHIS - 1),
            CosPhi::Centres => -1.0 + (2.0 * f64::from(self.cos_phi) + 1.0) / f64::from(COS_PHIS),
            CosPhi::Ladder => -1.0 + 2.0 * f64::from(self.cos_phi) / f64::from(COS_PHIS),
        };
        let sin_phi = (1.0 - cos_phi * cos_phi).sqrt();
        let normal = [
            sin_phi * Real::cos(theta),
            sin_phi * Real::sin(theta),
            cos_phi,
        ];
        // 9 values over the closed `-0.8..=0.8`: the step-0.2 ladder.
        let offset =
            2.0 * OFFSET_LIMIT * f64::from(self.offset) / f64::from(OFFSETS - 1) - OFFSET_LIMIT;
        (normal, offset)
    }
}

// ─── the field ─────────────────────────────────────────────────────────────

/// `✗43`'s fixture: one half-space against `convex_body`'s radius-1.5 bound.
///
/// `ConvexBody::sample`'s own arithmetic, with the bound taken from the crate's
/// own [`Sphere`] rather than re-derived, and deliberately **no** `gradient`
/// override — `ConvexBody` has none either, so both take the trait's central
/// difference and the normals are the same normals.
struct Cap {
    normal: [f64; 3],
    offset: f64,
}

impl Cap {
    fn of(config: Config, reading: Reading) -> Self {
        let (normal, offset) = config.plane(reading);
        Self { normal, offset }
    }
}

impl Sdf for Cap {
    type Scalar = f64;

    fn sample(&self, p: [f64; 3]) -> f64 {
        let bound = Sphere {
            center: [0.0; 3],
            radius: BOUND,
        }
        .sample(p);
        // `vec3::dot`'s order, which is the order `ConvexBody` evaluates in.
        let plane =
            self.normal[0] * p[0] + self.normal[1] * p[1] + self.normal[2] * p[2] - self.offset;
        bound.max(plane)
    }
}

/// What a sample's magnitude becomes once its sign has been read off the cap.
#[derive(Clone, Copy, PartialEq, Eq)]
enum Magnitudes {
    /// The cap's own value. Used only by the faithfulness control, which shows
    /// that the interpolant of these values *is* the cap as far as the
    /// extractor is concerned.
    Cap,
    /// `±1`. `M-374`'s bad fixture, kept because its result is C3.
    Unit,
    /// `sign · SplitMix64(seed, config, size, sample)` in `[1/4, 5/4)`.
    Generic(u64),
}

/// SplitMix64, so a magnitude is a pure function of its key and the sweep is
/// byte-identical on every machine and every run.
fn magnitude(seed: u64, key: u64, sample: usize) -> f64 {
    let mut z = seed
        .wrapping_add(key)
        .wrapping_mul(0x9E37_79B9_7F4A_7C15)
        .wrapping_add(sample as u64 + 1)
        .wrapping_mul(0xBF58_476D_1CE4_E5B9);
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    z ^= z >> 31;
    // Bounded away from zero, so no corner sits on the surface and the sign
    // pattern the cap defines is the sign pattern meshed.
    0.25 + (z >> 11) as f64 * (1.0 / 9_007_199_254_740_992.0)
}

/// The exact sample coordinates of one grid, per axis.
///
/// `origin + cell_size · k`, which is the expression `sdf::sample_grid` and
/// `marching_cubes::corner_position` both evaluate (`M-143`: `mul_add` is a
/// different number). Holding them means two things are exact rather than
/// tolerant: the interpolant reproduces a corner value at a corner, and "does
/// this coordinate lie on a sample plane" is an equality rather than an epsilon.
struct Lattice {
    coord: [Vec<f64>; 3],
}

impl Lattice {
    fn new(size: [u32; 3], origin: [f64; 3], cell_size: f64) -> Self {
        Self {
            coord: core::array::from_fn(|axis| {
                (0..size[axis])
                    .map(|k| origin[axis] + cell_size * f64::from(k))
                    .collect()
            }),
        }
    }

    /// Does `x` sit exactly on a sample plane of `axis`?
    ///
    /// An edge vertex's two constant coordinates are exactly a corner's, because
    /// `cube::place` reduces to `(lo + hi)·½` when `lo == hi` and that is exact.
    /// So this is an equality test and needs no tolerance.
    fn on_plane(&self, axis: usize, x: f64) -> bool {
        self.coord[axis].contains(&x)
    }

    /// The cell index containing an `x` that is strictly between two planes.
    fn cell(&self, axis: usize, x: f64) -> u32 {
        let planes = &self.coord[axis];
        let mut k = 0usize;
        while k + 2 < planes.len() && x > planes[k + 1] {
            k += 1;
        }
        k as u32
    }
}

/// The trilinear interpolant of one cap's **sign** pattern, with magnitudes
/// substituted.
///
/// The extractor reads a field only at lattice points (`crossing_refinement` is
/// zero by default, and normals do not enter any clause here), so this is the
/// same object Marching Cubes sees when handed the cap — which
/// [`faithfulness_control`] checks bit for bit rather than asserting.
struct Interpolant {
    lattice: Lattice,
    size: [usize; 3],
    cell_size: f64,
    value: Vec<f64>,
    /// Lattice points where the cap sampled exactly zero.
    ///
    /// A degeneracy detector, not a nuisance count: at such a point a cut edge's
    /// crossing lands exactly on the corner, and a mesh whose vertices coincide
    /// is a fixture artefact rather than an extractor defect. Reported so a C1
    /// failure can be told apart from one.
    zero_samples: u64,
}

impl Interpolant {
    fn of(cap: &Cap, grid: &Grid, magnitudes: Magnitudes, key: u64) -> Self {
        let n = grid.size as usize;
        let lattice = Lattice::new([grid.size; 3], grid.origin, grid.cell_size);
        let mut value = Vec::with_capacity(n * n * n);
        let mut zero_samples = 0u64;
        for z in 0..n {
            for y in 0..n {
                for x in 0..n {
                    let truth = cap.sample([
                        lattice.coord[0][x],
                        lattice.coord[1][y],
                        lattice.coord[2][z],
                    ]);
                    if truth == 0.0 {
                        zero_samples += 1;
                    }
                    // `cube::is_inside`: a sample of exactly zero is outside.
                    let sign = if truth < 0.0 { -1.0 } else { 1.0 };
                    let index = x + n * (y + n * z);
                    value.push(match magnitudes {
                        Magnitudes::Cap => truth,
                        Magnitudes::Unit => sign,
                        Magnitudes::Generic(seed) => sign * magnitude(seed, key, index),
                    });
                }
            }
        }
        Self {
            lattice,
            size: [n; 3],
            cell_size: grid.cell_size,
            value,
            zero_samples,
        }
    }

    fn at(&self, x: usize, y: usize, z: usize) -> f64 {
        self.value[x + self.size[0] * (y + self.size[1] * z)]
    }

    /// The containing cell's base index and the fraction inside it.
    ///
    /// **A sample plane has to weigh exactly 0 or 1, and getting that wrong was
    /// worth a run.** The first version of this deferred to
    /// [`Lattice::cell`], which answers "which cell is `x` strictly inside" and
    /// therefore returns the *lower* cell for an `x` sitting exactly on a plane.
    /// The fraction was then `(planes[j] − planes[j−1])/h`, which is `1.0` to
    /// within an ulp and not `1.0` — so the interpolant blended two corners
    /// where it should have returned one, and [`faithfulness_control`] caught it
    /// as a one-ulp vertex displacement on the first configuration it tried.
    /// That is the control doing its job: the defect is invisible in every
    /// aggregate this experiment reports.
    ///
    /// Outside the grid the base clamps — a boundary condition rather than a
    /// fallback: only the trait's central-difference gradient ever asks, and it
    /// asks for a normal.
    fn locate(&self, axis: usize, x: f64) -> (usize, f64) {
        let planes = &self.lattice.coord[axis];
        let last = planes.len() - 2;
        if let Some(j) = planes.iter().position(|&c| c == x) {
            return if j <= last { (j, 0.0) } else { (last, 1.0) };
        }
        let mut k = 0usize;
        while k < last && x > planes[k + 1] {
            k += 1;
        }
        (k, (x - planes[k]) / self.cell_size)
    }
}

impl Sdf for Interpolant {
    type Scalar = f64;

    fn sample(&self, p: [f64; 3]) -> f64 {
        let (bx, fx) = self.locate(0, p[0]);
        let (by, fy) = self.locate(1, p[1]);
        let (bz, fz) = self.locate(2, p[2]);
        let mut acc = 0.0;
        for k in 0..8u32 {
            let (dx, dy, dz) = (
                (k & 1) as usize,
                ((k >> 1) & 1) as usize,
                ((k >> 2) & 1) as usize,
            );
            let wx = if dx == 1 { fx } else { 1.0 - fx };
            let wy = if dy == 1 { fy } else { 1.0 - fy };
            let wz = if dz == 1 { fz } else { 1.0 - fz };
            acc += wx * wy * wz * self.at(bx + dx, by + dy, bz + dz);
        }
        acc
    }
}

// ─── the grid ──────────────────────────────────────────────────────────────

/// One grid of `✗43`'s sweep: `property::extraction::grid_for`'s own arithmetic.
struct Grid {
    size: u32,
    shape: RuntimeShape3,
    origin: [f64; 3],
    cell_size: f64,
    validate: ValidateConfig,
    lattice: Lattice,
}

impl Grid {
    fn of(size: u32) -> Self {
        let cell_size = 2.0 * DOMAIN / f64::from(size - 1);
        Self {
            size,
            shape: RuntimeShape3::new([size; 3]).expect("a cubic grid of 6 to 12 samples"),
            origin: [-DOMAIN; 3],
            cell_size,
            validate: ValidateConfig::from_cell_size(cell_size)
                .expect("a positive finite cell size"),
            lattice: Lattice::new([size; 3], [-DOMAIN; 3], cell_size),
        }
    }

    /// `emit_trilinear`'s own `to_world`, so a saddle lands where the extractor
    /// would have put it.
    fn to_world(&self, base: [u32; 3], local: [f64; 3]) -> [f64; 3] {
        core::array::from_fn(|k| self.origin[k] + self.cell_size * (f64::from(base[k]) + local[k]))
    }

    /// The cells of this grid, in the extractor's own order.
    fn cells(&self) -> impl Iterator<Item = [u32; 3]> + '_ {
        let n = self.size - 1;
        (0..n).flat_map(move |z| (0..n).flat_map(move |y| (0..n).map(move |x| [x, y, z])))
    }
}

// ─── what the crate's own trilinear machinery says, cell by cell ───────────

/// A cell that fans two or more rings, and therefore differs before and after
/// the fix.
struct Fanned {
    base: [u32; 3],
    /// Rings longer than three: the apexes `Contours::fan` names.
    rings: usize,
    /// `BodySaddles::interior_vertex` in world coordinates — the one apex the
    /// pre-fix code created for all of them.
    saddle: [f64; 3],
}

/// Per-cell analysis of one grid, from the crate's own public tables.
struct Analysis {
    ambiguous_cells: u64,
    trilinear_cells: u64,
    hexagon_cells: u64,
    /// Trilinear cells where the strict `0 < x < 1` test rejected the hexagon.
    hexagon_rejects: u64,
    /// Trilinear cells with no inner hexagon and no interior vertex either:
    /// `M-374`'s mechanism, counted.
    interior_vertex_rejects: u64,
    /// Apexes the disk path creates: one per fanned ring.
    apexes: u64,
    fanned: Vec<Fanned>,
    /// Cells whose interior vertex exists but whose every ring is a triangle.
    /// Pre-fix these were created anyway and referenced by nothing.
    unfanned: Vec<[f64; 3]>,
}

/// Re-run the extractor's own cell classification, using only public API.
///
/// Not a second implementation of the extractor: it names the same tables and
/// the same two constructors the extractor names, and the *only* thing it
/// derives is how many apexes each cell asks for. That number is then checked
/// against the mesh.
fn analyse(field: &Interpolant, grid: &Grid) -> Analysis {
    let mut out = Analysis {
        ambiguous_cells: 0,
        trilinear_cells: 0,
        hexagon_cells: 0,
        hexagon_rejects: 0,
        interior_vertex_rejects: 0,
        apexes: 0,
        fanned: Vec::new(),
        unfanned: Vec::new(),
    };
    for base in grid.cells() {
        let mut case = 0u8;
        let mut corner = [0.0f64; 8];
        for (c, slot) in corner.iter_mut().enumerate() {
            // `cube::corner_offset`'s bit layout.
            let o = [(c & 1) as u32, ((c >> 1) & 1) as u32, ((c >> 2) & 1) as u32];
            let v = field.at(
                (base[0] + o[0]) as usize,
                (base[1] + o[1]) as usize,
                (base[2] + o[2]) as usize,
            );
            *slot = v;
            // `cube::is_inside`: exactly zero is outside.
            if v < 0.0 {
                case |= 1 << c;
            }
        }

        let ambiguous = AMBIGUOUS_FACES[case as usize];
        if ambiguous == 0 {
            continue;
        }
        out.ambiguous_cells += 1;
        let mask = joined_mask(&corner, ambiguous);
        let contours = Contours::of(case, mask);
        if contours.count() == 0 {
            continue;
        }
        out.trilinear_cells += 1;

        let saddles = BodySaddles::of(&corner);
        if saddles.has_inner_hexagon() {
            // The tunnel path: six interior vertices, and the fix did not touch
            // how they are named.
            out.hexagon_cells += 1;
            continue;
        }
        out.hexagon_rejects += 1;
        let Some(local) = saddles.interior_vertex() else {
            out.interior_vertex_rejects += 1;
            continue;
        };
        let rings = (0..contours.count())
            .filter(|&r| contours.ring(r).len() > 3)
            .count();
        out.apexes += rings as u64;
        let saddle = grid.to_world(base, local);
        if rings >= 2 {
            out.fanned.push(Fanned {
                base,
                rings,
                saddle,
            });
        } else if rings == 0 {
            out.unfanned.push(saddle);
        }
    }
    out
}

// ─── the pre-fix fan ───────────────────────────────────────────────────────

/// The mesh `✗43`'s code would have produced, and how much of it moved.
struct PreFix {
    positions: Vec<[f64; 3]>,
    indices: Vec<u32>,
    /// Apexes identified away: `k` rings sharing one apex costs `k − 1`.
    merges: u64,
    /// Cells where the analysis predicted `rings` apexes and the mesh did not
    /// hold that many. Zero, or the instrument is wrong.
    mismatches: u64,
}

/// Which cell each cell-interior vertex of a mesh belongs to.
///
/// A vertex is cell-interior when **no** coordinate lies on a sample plane. That
/// is `✗43`'s own correction applied to the classification: an edge vertex has
/// two coordinates on the lattice and a face vertex one, so counting them is the
/// test and the count alone is enough here because the only thing needed is
/// "strictly inside some cell".
fn interior_cells(positions: &[[f64; 3]], lattice: &Lattice) -> Vec<Option<[u32; 3]>> {
    positions
        .iter()
        .map(|p| {
            if (0..3).any(|axis| lattice.on_plane(axis, p[axis])) {
                None
            } else {
                Some(core::array::from_fn(|axis| lattice.cell(axis, p[axis])))
            }
        })
        .collect()
}

/// Undo the per-ring apex: one apex per cell, back on the body saddle.
fn pre_fix(mesh: &MeshBuffer<f64>, cells: &[Option<[u32; 3]>], analysis: &Analysis) -> PreFix {
    let mut positions = mesh.positions.clone();
    let mut remap: Vec<u32> = (0..positions.len() as u32).collect();
    let mut merges = 0u64;
    let mut mismatches = 0u64;

    for cell in &analysis.fanned {
        let apexes: Vec<u32> = cells
            .iter()
            .enumerate()
            .filter(|(_, c)| **c == Some(cell.base))
            .map(|(v, _)| v as u32)
            .collect();
        if apexes.len() != cell.rings {
            mismatches += 1;
        }
        let Some((&keep, rest)) = apexes.split_first() else {
            continue;
        };
        // The pre-fix apex is the body saddle for every ring, not a centroid.
        positions[keep as usize] = cell.saddle;
        for &v in rest {
            remap[v as usize] = keep;
            merges += 1;
        }
    }

    let indices = mesh.indices.iter().map(|&i| remap[i as usize]).collect();
    // Pre-fix, an apex was created whenever the interior vertex existed, even
    // where no ring was long enough to fan it. Reproduced because it is part of
    // what that code did; it cannot move `is_closed()`, since χ counts
    // referenced vertices.
    positions.extend_from_slice(&analysis.unfanned);
    PreFix {
        positions,
        indices,
        merges,
        mismatches,
    }
}

// ─── one arm ───────────────────────────────────────────────────────────────

/// One magnitude arm.
#[derive(Clone, Copy)]
struct Arm {
    name: &'static str,
    seed: &'static str,
    magnitudes: Magnitudes,
}

/// Everything one `(arm, size)` row reports.
struct Sweep {
    arm: Arm,
    size: u32,
    cell_size: f64,
    configurations: u64,
    cells: u64,
    ambiguous_cells: u64,
    trilinear_cells: u64,
    hexagon_cells: u64,
    hexagon_rejects: u64,
    interior_vertex_rejects: u64,
    apexes: u64,
    interior_vertices: u64,
    fan_configurations: u64,
    apex_merges: u64,
    unfanned_apex_cells: u64,
    unclosed_post_fix: u64,
    unclosed_pre_fix: u64,
    nmv_post_fix: u64,
    nmv_pre_fix: u64,
    unreferenced_post_fix: u64,
    unreferenced_pre_fix: u64,
    triangles: u64,
    empty_meshes: u64,
    zero_sample_configurations: u64,
    classification_mismatches: u64,
    /// Up to four configuration labels, `+`-joined, or `none`.
    ///
    /// A count says how many and a label says **which**, and which is what
    /// another run needs to replay one cap. Four because the interesting answers
    /// here are 0 and 2.
    unclosed_post_fix_where: String,
    unclosed_pre_fix_where: String,
    wall_ms: u128,
}

impl Sweep {
    /// The arm never fired the interior rule, so its zeros are the
    /// interior-rule-off arm's zeros wearing a different name (`M-374`).
    fn void(&self) -> bool {
        self.apexes == 0
    }

    /// Zero unclosed meshes, on an arm that could have reported one.
    fn c1(&self) -> bool {
        !self.void() && self.empty_meshes == 0 && self.unclosed_post_fix == 0
    }

    /// `✗43`'s figure: two at 6³, none elsewhere.
    ///
    /// **Gated on `fan_configurations`, which is C2's share, and that gate is
    /// the whole reason this column is worth reading.** A row where no cell
    /// fanned two rings is a row where the pre-fix and post-fix meshes are
    /// identical by construction, so its `unclosed_pre_fix` of zero is not a
    /// measurement of the pre-fix fan — it is a measurement of the fan's
    /// absence. Reporting that as a held clause is `M-44` in the direction
    /// nobody checks, so it reads `false`, and the share it could have moved is
    /// on the row beside it.
    fn c2(&self) -> bool {
        if self.void() || self.fan_configurations == 0 {
            return false;
        }
        if self.size == 6 {
            self.unclosed_pre_fix == 2
        } else {
            self.unclosed_pre_fix == 0
        }
    }

    /// The share bound C2's arithmetic rests on, as a check rather than a claim.
    ///
    /// Undoing the fix can only unclose a mesh in a configuration that has a
    /// cell fanning two or more rings; every other configuration comes out of
    /// `pre_fix` unchanged. So `unclosed_pre_fix` cannot exceed
    /// `unclosed_post_fix + fan_configurations`, and if it does the
    /// reconstruction is touching vertices it has no business touching.
    fn share_bound_holds(&self) -> bool {
        self.unclosed_pre_fix <= self.unclosed_post_fix + self.fan_configurations
    }

    /// The unit arm must report VOID; every other arm must be live.
    fn c3(&self) -> bool {
        if self.arm.magnitudes == Magnitudes::Unit {
            self.void()
        } else {
            !self.void()
        }
    }

    fn row(&self, reading: Reading) -> Vec<(&'static str, String)> {
        vec![
            ("arm", self.arm.name.to_string()),
            ("seed", self.arm.seed.to_string()),
            ("size", self.size.to_string()),
            ("reading", reading.label().to_string()),
            ("cell_size", format!("{:.6}", self.cell_size)),
            ("configurations", self.configurations.to_string()),
            ("cells", self.cells.to_string()),
            ("ambiguous_cells", self.ambiguous_cells.to_string()),
            ("trilinear_cells", self.trilinear_cells.to_string()),
            ("hexagon_cells", self.hexagon_cells.to_string()),
            (
                "has_inner_hexagon_rejects",
                self.hexagon_rejects.to_string(),
            ),
            (
                "interior_vertex_rejects",
                self.interior_vertex_rejects.to_string(),
            ),
            ("interior_apexes", self.apexes.to_string()),
            ("interior_vertices", self.interior_vertices.to_string()),
            ("fan_configurations", self.fan_configurations.to_string()),
            ("apex_merges", self.apex_merges.to_string()),
            ("unfanned_apex_cells", self.unfanned_apex_cells.to_string()),
            ("unclosed_post_fix", self.unclosed_post_fix.to_string()),
            ("unclosed_pre_fix", self.unclosed_pre_fix.to_string()),
            (
                "non_manifold_vertices_post_fix",
                self.nmv_post_fix.to_string(),
            ),
            (
                "non_manifold_vertices_pre_fix",
                self.nmv_pre_fix.to_string(),
            ),
            (
                "unreferenced_post_fix",
                self.unreferenced_post_fix.to_string(),
            ),
            // **What this counts, exactly.** `pre_fix` leaves a merged-away
            // apex in `positions` and stops referencing it, and it appends the
            // apexes pre-fix created for cells no ring fanned. So this is
            // `apex_merges + unfanned_apex_cells`, and only the second term is
            // the artefact `✗43` describes. Neither can move `is_closed()`,
            // because χ counts referenced vertices — which is why the pre-fix
            // report can be trusted to differ from the post-fix one only where
            // the fan itself differs.
            (
                "unreferenced_pre_fix",
                self.unreferenced_pre_fix.to_string(),
            ),
            ("triangles", self.triangles.to_string()),
            ("empty_meshes", self.empty_meshes.to_string()),
            (
                "zero_sample_configurations",
                self.zero_sample_configurations.to_string(),
            ),
            (
                "classification_mismatches",
                self.classification_mismatches.to_string(),
            ),
            (
                "unclosed_post_fix_where",
                self.unclosed_post_fix_where.clone(),
            ),
            (
                "unclosed_pre_fix_where",
                self.unclosed_pre_fix_where.clone(),
            ),
            ("void", self.void().to_string()),
            ("c1_holds", self.c1().to_string()),
            ("c2_holds", self.c2().to_string()),
            ("c3_holds", self.c3().to_string()),
            ("wall_ms", self.wall_ms.to_string()),
        ]
    }
}

/// A configured extractor: `check_trilinear`'s settings, which are the ones
/// `✗43` failed and fixed under.
fn extractor() -> MarchingCubes<f64> {
    let mut mc = MarchingCubes::<f64>::new();
    mc.set_face_ambiguity(FaceAmbiguity::AsymptoticDecider);
    mc.set_interior_ambiguity(InteriorAmbiguity::Trilinear);
    mc
}

/// Mesh one cap at one size, both ways, and report both.
struct One {
    post: MeshReport,
    pre: MeshReport,
    triangles: u64,
    interior_vertices: u64,
    apexes: u64,
    ambiguous_cells: u64,
    trilinear_cells: u64,
    hexagon_cells: u64,
    hexagon_rejects: u64,
    interior_vertex_rejects: u64,
    fanned_cells: u64,
    merges: u64,
    unfanned: u64,
    mismatches: u64,
    zero_samples: u64,
}

fn one(
    cap: &Cap,
    grid: &Grid,
    magnitudes: Magnitudes,
    key: u64,
    mc: &mut MarchingCubes<f64>,
    mesh: &mut MeshBuffer<f64>,
) -> One {
    let field = Interpolant::of(cap, grid, magnitudes, key);
    let analysis = analyse(&field, grid);

    mesh.reset();
    mc.extract(&field, &grid.shape, grid.origin, grid.cell_size, mesh)
        .expect("marching cubes on a single-plane cap inside its own bound");

    let cells = interior_cells(&mesh.positions, &grid.lattice);
    let interior_vertices = cells.iter().filter(|c| c.is_some()).count() as u64;
    let post = validate_indexed(&mesh.positions, &mesh.indices, &grid.validate);
    let fixed = pre_fix(mesh, &cells, &analysis);
    let pre = validate_indexed(&fixed.positions, &fixed.indices, &grid.validate);

    One {
        post,
        pre,
        triangles: (mesh.indices.len() / 3) as u64,
        interior_vertices,
        apexes: analysis.apexes,
        ambiguous_cells: analysis.ambiguous_cells,
        trilinear_cells: analysis.trilinear_cells,
        hexagon_cells: analysis.hexagon_cells,
        hexagon_rejects: analysis.hexagon_rejects,
        interior_vertex_rejects: analysis.interior_vertex_rejects,
        fanned_cells: analysis.fanned.len() as u64,
        merges: fixed.merges,
        unfanned: analysis.unfanned.len() as u64,
        mismatches: fixed.mismatches,
        zero_samples: field.zero_samples,
    }
}

/// Name a configuration in a `+`-joined list, up to four of them.
///
/// The separator is `+` and not a comma because the CSV writer does not quote,
/// and it refuses a value that would shift every column after it.
fn note(slot: &mut String, config: &Config) {
    if slot == "none" {
        *slot = config.label();
    } else if slot.matches('+').count() < 3 {
        slot.push('+');
        slot.push_str(&config.label());
    }
}

/// Sweep every configuration at one size, on one magnitude arm.
fn sweep(arm: Arm, size: u32, reading: Reading, configs: &[Config]) -> Sweep {
    let grid = Grid::of(size);
    let mut mc = extractor();
    let mut mesh = MeshBuffer::<f64>::new();
    let started = Instant::now();

    let mut out = Sweep {
        arm,
        size,
        cell_size: grid.cell_size,
        configurations: configs.len() as u64,
        cells: u64::from(size - 1).pow(3) * configs.len() as u64,
        ambiguous_cells: 0,
        trilinear_cells: 0,
        hexagon_cells: 0,
        hexagon_rejects: 0,
        interior_vertex_rejects: 0,
        apexes: 0,
        interior_vertices: 0,
        fan_configurations: 0,
        apex_merges: 0,
        unfanned_apex_cells: 0,
        unclosed_post_fix: 0,
        unclosed_pre_fix: 0,
        nmv_post_fix: 0,
        nmv_pre_fix: 0,
        unreferenced_post_fix: 0,
        unreferenced_pre_fix: 0,
        triangles: 0,
        empty_meshes: 0,
        zero_sample_configurations: 0,
        classification_mismatches: 0,
        unclosed_post_fix_where: String::from("none"),
        unclosed_pre_fix_where: String::from("none"),
        wall_ms: 0,
    };

    for (id, config) in configs.iter().enumerate() {
        let cap = Cap::of(*config, reading);
        let key = id as u64 | (u64::from(size) << 32);
        let r = one(&cap, &grid, arm.magnitudes, key, &mut mc, &mut mesh);

        out.ambiguous_cells += r.ambiguous_cells;
        out.trilinear_cells += r.trilinear_cells;
        out.hexagon_cells += r.hexagon_cells;
        out.hexagon_rejects += r.hexagon_rejects;
        out.interior_vertex_rejects += r.interior_vertex_rejects;
        out.apexes += r.apexes;
        out.interior_vertices += r.interior_vertices;
        out.apex_merges += r.merges;
        out.unfanned_apex_cells += r.unfanned;
        out.classification_mismatches += r.mismatches;
        out.triangles += r.triangles;
        out.nmv_post_fix += r.post.non_manifold_vertices;
        out.nmv_pre_fix += r.pre.non_manifold_vertices;
        out.unreferenced_post_fix += r.post.unreferenced_vertices;
        out.unreferenced_pre_fix += r.pre.unreferenced_vertices;
        out.fan_configurations += u64::from(r.fanned_cells > 0);
        out.zero_sample_configurations += u64::from(r.zero_samples > 0);
        out.empty_meshes += u64::from(r.triangles == 0);
        if !r.post.is_closed() {
            out.unclosed_post_fix += 1;
            note(&mut out.unclosed_post_fix_where, config);
        }
        if !r.pre.is_closed() {
            out.unclosed_pre_fix += 1;
            note(&mut out.unclosed_pre_fix_where, config);
        }
    }

    out.wall_ms = started.elapsed().as_millis();
    out
}

// ─── the controls that run before the sweep ────────────────────────────────

/// The interpolant of the cap's own samples **is** the cap, to the extractor.
///
/// Positions and indices compared bit for bit over every configuration at both
/// ends of the size range. This is what licenses the magnitude arms: if the
/// requantisation changed the mesh, a generic arm would be a different
/// experiment rather than the same one at another magnitude draw.
fn faithfulness_control(reading: Reading, configs: &[Config]) -> u64 {
    let mut identical = 0u64;
    let mut mc = extractor();
    let mut cap_mesh = MeshBuffer::<f64>::new();
    let mut interp_mesh = MeshBuffer::<f64>::new();
    for size in [SIZES[0], SIZES[SIZES.len() - 1]] {
        let grid = Grid::of(size);
        for config in configs {
            let cap = Cap::of(*config, reading);
            let field = Interpolant::of(&cap, &grid, Magnitudes::Cap, 0);
            cap_mesh.reset();
            mc.extract(
                &cap,
                &grid.shape,
                grid.origin,
                grid.cell_size,
                &mut cap_mesh,
            )
            .expect("marching cubes on the cap");
            interp_mesh.reset();
            mc.extract(
                &field,
                &grid.shape,
                grid.origin,
                grid.cell_size,
                &mut interp_mesh,
            )
            .expect("marching cubes on the interpolant of the cap");
            assert_eq!(
                cap_mesh.positions,
                interp_mesh.positions,
                "the interpolant of the cap's samples moved a vertex at {} on {size}^3, so the \
                 magnitude arms would not be the same experiment",
                config.label()
            );
            assert_eq!(
                cap_mesh.indices,
                interp_mesh.indices,
                "the interpolant of the cap's samples changed the triangulation at {} on {size}^3",
                config.label()
            );
            identical += 1;
        }
    }
    identical
}

/// `✗43`'s own reduced counterexample, before and after.
///
/// **This is the control the entry lacked and the reason the pre-fix arm can be
/// believed.** The entry publishes both reports for this cell at 6³, so the
/// bench-local pre-fix fan is calibrated against a number it did not choose:
/// pre-fix `V 27, E 72, F 48, χ 3, components 2, non-manifold vertices 1,
/// degenerate triangles 1`, `is_closed()` false; post-fix the same with
/// `V 28, χ 4`, `non_manifold_vertices 0`, `degenerate_triangles 0` and
/// `is_closed()` true. If the reconstruction cannot reproduce that, the sweep's
/// pre-fix column means nothing — and `M-44` is satisfied by demonstration
/// rather than by argument, because the instrument is shown returning non-zero.
fn x43_control() -> Sweep {
    let grid = Grid::of(6);
    let mut mc = extractor();
    let mut mesh = MeshBuffer::<f64>::new();
    let started = Instant::now();
    let cap = Cap {
        normal: [
            0.237_655_829_452_438_93,
            0.969_345_284_197_866_1,
            0.062_365_268_624_701_24,
        ],
        offset: -0.603_027_951_078_954_2,
    };
    let r = one(&cap, &grid, Magnitudes::Cap, 0, &mut mc, &mut mesh);

    println!(
        "x43 fixture at 6^3, post-fix: V {} E {} F {} chi {} components {} nmv {} degenerate {} \
         closed {}",
        r.post.referenced_vertices,
        r.post.edges,
        r.post.faces,
        r.post.euler_characteristic,
        r.post.components,
        r.post.non_manifold_vertices,
        r.post.degenerate_triangles,
        r.post.is_closed()
    );
    println!(
        "x43 fixture at 6^3, pre-fix:  V {} E {} F {} chi {} components {} nmv {} degenerate {} \
         closed {}",
        r.pre.referenced_vertices,
        r.pre.edges,
        r.pre.faces,
        r.pre.euler_characteristic,
        r.pre.components,
        r.pre.non_manifold_vertices,
        r.pre.degenerate_triangles,
        r.pre.is_closed()
    );

    // The entry's published post-fix numbers, which the committed test asserts.
    assert_eq!(r.post.referenced_vertices, 28, "x43 post-fix vertices");
    assert_eq!(r.post.edges, 72, "x43 post-fix edges");
    assert_eq!(r.post.faces, 48, "x43 post-fix faces");
    assert_eq!(r.post.euler_characteristic, 4, "x43 post-fix chi");
    assert_eq!(r.post.components, 2, "x43 post-fix components");
    assert_eq!(r.post.non_manifold_vertices, 0, "x43 post-fix nmv");
    assert!(r.post.is_closed(), "x43 post-fix closed");
    // The entry's published *pre-fix* numbers, which nothing has asserted until
    // now. This is the fixture-can-fail control.
    assert_eq!(r.pre.referenced_vertices, 27, "x43 pre-fix vertices");
    assert_eq!(r.pre.edges, 72, "x43 pre-fix edges");
    assert_eq!(r.pre.faces, 48, "x43 pre-fix faces");
    assert_eq!(r.pre.euler_characteristic, 3, "x43 pre-fix chi");
    assert_eq!(r.pre.components, 2, "x43 pre-fix components");
    assert_eq!(r.pre.non_manifold_vertices, 1, "x43 pre-fix nmv");
    assert_eq!(r.pre.degenerate_triangles, 1, "x43 pre-fix degenerate");
    assert!(
        !r.pre.is_closed(),
        "VOID: the pre-fix fan produced a closed mesh on the one cell known to \
         have failed, so the pre-fix column cannot report bad news"
    );

    Sweep {
        arm: Arm {
            name: "x43_fixture",
            seed: "none",
            magnitudes: Magnitudes::Cap,
        },
        size: 6,
        cell_size: grid.cell_size,
        configurations: 1,
        cells: 125,
        ambiguous_cells: r.ambiguous_cells,
        trilinear_cells: r.trilinear_cells,
        hexagon_cells: r.hexagon_cells,
        hexagon_rejects: r.hexagon_rejects,
        interior_vertex_rejects: r.interior_vertex_rejects,
        apexes: r.apexes,
        interior_vertices: r.interior_vertices,
        fan_configurations: u64::from(r.fanned_cells > 0),
        apex_merges: r.merges,
        unfanned_apex_cells: r.unfanned,
        unclosed_post_fix: u64::from(!r.post.is_closed()),
        unclosed_pre_fix: u64::from(!r.pre.is_closed()),
        nmv_post_fix: r.post.non_manifold_vertices,
        nmv_pre_fix: r.pre.non_manifold_vertices,
        unreferenced_post_fix: r.post.unreferenced_vertices,
        unreferenced_pre_fix: r.pre.unreferenced_vertices,
        triangles: r.triangles,
        empty_meshes: u64::from(r.triangles == 0),
        zero_sample_configurations: u64::from(r.zero_samples > 0),
        classification_mismatches: r.mismatches,
        unclosed_post_fix_where: String::from("none"),
        unclosed_pre_fix_where: String::from("x43_plane"),
        wall_ms: started.elapsed().as_millis(),
    }
}

/// The pre-fix count at 6³ under every reading the entry's wording admits.
///
/// Printed rather than recorded: the CSV's population is the chosen reading's,
/// and a table of six readings is evidence for the choice rather than a result.
fn calibrate(configs: &[Config]) {
    println!(
        "\n{:<32} {:>8} {:>10} {:>10} {:>10} {:>10}",
        "reading", "fanCfg", "preUncl", "postUncl", "apexes", "zeroSmpl"
    );
    for reading in READINGS {
        let s = sweep(
            Arm {
                name: "exact",
                seed: "none",
                magnitudes: Magnitudes::Cap,
            },
            6,
            reading,
            configs,
        );
        println!(
            "{:<32} {:>8} {:>10} {:>10} {:>10} {:>10}{}",
            reading.label(),
            s.fan_configurations,
            s.unclosed_pre_fix,
            s.unclosed_post_fix,
            s.apexes,
            s.zero_sample_configurations,
            if reading == READING {
                "  <- chosen"
            } else {
                ""
            }
        );
    }
}

// ─── the run ───────────────────────────────────────────────────────────────

fn main() {
    if !std::env::args().any(|arg| arg == "--bench") {
        return;
    }
    let prereg = isomesh::experiment!("P-102");
    common::experiment::run(prereg, |run| {
        // The registered population, asserted rather than trusted.
        assert_eq!(CONFIGURATIONS, 1_152, "16 x 8 x 9 caps");
        assert_eq!(
            CONFIGURATIONS as usize * SIZES.len(),
            8_064,
            "1,152 caps at each of seven sizes"
        );
        assert_eq!(SIZES[0], 6, "the sweep starts at 6^3");
        assert_eq!(SIZES[SIZES.len() - 1], 12, "the sweep ends at 12^3");
        let configs = Config::all();
        assert_eq!(configs.len(), CONFIGURATIONS as usize);

        println!(
            "P-102: {} caps ({THETAS} theta x {COS_PHIS} cos phi x {OFFSETS} offsets) at each of \
             {:?}, reading {}",
            configs.len(),
            SIZES,
            READING.label()
        );

        let identical = faithfulness_control(READING, &configs);
        println!(
            "requantisation control: {identical}/{} meshes bit-identical to the analytic cap",
            2 * CONFIGURATIONS
        );
        assert_eq!(
            identical,
            u64::from(2 * CONFIGURATIONS),
            "the interpolant must reproduce the cap on every configuration"
        );

        let control = x43_control();
        calibrate(&configs);

        let arms = [
            Arm {
                name: "exact",
                seed: "none",
                magnitudes: Magnitudes::Cap,
            },
            Arm {
                name: "unit",
                seed: "none",
                magnitudes: Magnitudes::Unit,
            },
            Arm {
                name: "generic",
                seed: "0x00002026",
                magnitudes: Magnitudes::Generic(SEEDS[0]),
            },
            Arm {
                name: "generic",
                seed: "0x0005eed1",
                magnitudes: Magnitudes::Generic(SEEDS[1]),
            },
            Arm {
                name: "generic",
                seed: "0x00c0ffee",
                magnitudes: Magnitudes::Generic(SEEDS[2]),
            },
            Arm {
                name: "generic",
                seed: "0xdeadbeef",
                magnitudes: Magnitudes::Generic(SEEDS[3]),
            },
        ];

        println!(
            "\n{:<8} {:<11} {:>5} {:>7} {:>7} {:>8} {:>8} {:>7} {:>7} {:>6} {:>7}",
            "arm",
            "seed",
            "size",
            "postUnc",
            "preUnc",
            "fanCfg",
            "apexes",
            "intV",
            "hexRej",
            "void",
            "ms"
        );
        let mut sweeps = Vec::new();
        for arm in arms {
            for size in SIZES {
                let s = sweep(arm, size, READING, &configs);
                println!(
                    "{:<8} {:<11} {:>5} {:>7} {:>7} {:>8} {:>8} {:>7} {:>7} {:>6} {:>7}",
                    s.arm.name,
                    s.arm.seed,
                    s.size,
                    s.unclosed_post_fix,
                    s.unclosed_pre_fix,
                    s.fan_configurations,
                    s.apexes,
                    s.interior_vertices,
                    s.hexagon_rejects,
                    s.void(),
                    s.wall_ms
                );
                sweeps.push(s);
            }
        }

        // ── per-arm totals, reported per seed and never pooled ──────────────
        println!(
            "\n{:<8} {:<11} {:>9} {:>9} {:>9} {:>9} {:>9} {:>9}",
            "arm", "seed", "configs", "postUnc", "preUnc", "apexes", "fanCfg", "tris"
        );
        for arm in arms {
            let rows: Vec<&Sweep> = sweeps
                .iter()
                .filter(|s| s.arm.name == arm.name && s.arm.seed == arm.seed)
                .collect();
            let total = |f: fn(&Sweep) -> u64| rows.iter().map(|s| f(s)).sum::<u64>();
            println!(
                "{:<8} {:<11} {:>9} {:>9} {:>9} {:>9} {:>9} {:>9}",
                arm.name,
                arm.seed,
                total(|s| s.configurations),
                total(|s| s.unclosed_post_fix),
                total(|s| s.unclosed_pre_fix),
                total(|s| s.apexes),
                total(|s| s.fan_configurations),
                total(|s| s.triangles)
            );
        }

        // ── the vacuity control, per seed, exactly as registered ────────────
        //
        // **A measurement that comes back zero must prove it could have come
        // back non-zero.** Three separate ways here, and none of them is implied
        // by another assertion in this file:
        //
        // - the interior rule must have fired on every live arm, or C1's zero is
        //   the interior-rule-off arm's zero (M-374, and the registered control);
        // - the sweep must have meshed something, or "no unclosed mesh" is a
        //   statement about an empty output;
        // - the per-cell analysis must agree with the mesh about how many apexes
        //   a fanned cell holds, or the pre-fix column is undoing the wrong
        //   vertices.
        for arm in arms {
            let rows: Vec<&Sweep> = sweeps
                .iter()
                .filter(|s| s.arm.name == arm.name && s.arm.seed == arm.seed)
                .collect();
            let apexes: u64 = rows.iter().map(|s| s.apexes).sum();
            let empty: u64 = rows.iter().map(|s| s.empty_meshes).sum();
            let mismatches: u64 = rows.iter().map(|s| s.classification_mismatches).sum();
            let triangles: u64 = rows.iter().map(|s| s.triangles).sum();
            assert_eq!(
                empty, 0,
                "{}/{}: {empty} of 8,064 configurations meshed nothing, so 'no unclosed \
                 mesh' would be a statement about an empty output",
                arm.name, arm.seed
            );
            assert!(
                triangles > 0,
                "{}/{}: the whole arm produced no triangles",
                arm.name,
                arm.seed
            );
            if arm.magnitudes == Magnitudes::Unit {
                // **Nothing about the unit arm's apex count is asserted here,
                // and that omission is deliberate.** C3 *is* the statement that
                // this arm comes back void; asserting `apexes == 0` would make
                // C3 a clause whose predicate an assertion in its own harness
                // already forces, which is `P-70`'s C3 — a HELD with no
                // instrument, the worst outcome available. So the unit arm's
                // `interior_apexes`, `void`, `c1_holds` and `c2_holds` are
                // measured, written to the CSV and printed, and if the interior
                // rule ever does fire at `±1` this run *completes* and the row
                // says so. What is asserted instead is the half that is a
                // control rather than a clause: the live arms above must return
                // non-zero, so the detector is shown able to answer both ways
                // inside one walk.
                continue;
            }
            assert!(
                apexes > 0,
                "VOID: {}/{} produced no interior apex over 8,064 configurations, so its \
                 zero is the interior-rule-off arm's zero wearing a different name",
                arm.name,
                arm.seed
            );
            assert_eq!(
                mismatches, 0,
                "{}/{}: the per-cell analysis and the mesh disagreed about the apex count \
                 in {mismatches} fanned cells, so the pre-fix arm is undoing the wrong \
                 vertices",
                arm.name, arm.seed
            );
        }

        // ── what the unit arm reports, measured and not asserted ────────────
        println!("\nC3, the unit arm, per size:");
        for s in sweeps.iter().filter(|s| s.arm.name == "unit") {
            println!(
                "   {}^3: trilinear cells {}  hexagon rejects {}  interior-vertex rejects {}  \
                 apexes {}  interior vertices {}  void {}  c1 {}  c2 {}",
                s.size,
                s.trilinear_cells,
                s.hexagon_rejects,
                s.interior_vertex_rejects,
                s.apexes,
                s.interior_vertices,
                s.void(),
                s.c1(),
                s.c2()
            );
        }

        // ── the share bound, on every row including the control ─────────────
        for s in std::iter::once(&control).chain(sweeps.iter()) {
            assert!(
                s.share_bound_holds(),
                "{}/{} at {}^3 unclosed {} meshes pre-fix with only {} fan configurations \
                 and {} unclosed post-fix, so the reconstruction moved a vertex the fix \
                 never touched",
                s.arm.name,
                s.arm.seed,
                s.size,
                s.unclosed_pre_fix,
                s.fan_configurations,
                s.unclosed_post_fix
            );
        }

        // ── two instruments for one quantity ────────────────────────────────
        //
        // `interior_apexes` is predicted from the crate's own public tables and
        // `interior_vertices` is measured off the mesh by coordinate, so on any
        // grid where the classification is unambiguous they must satisfy
        // `interior_vertices = interior_apexes + 6·hexagon_cells` — the tunnel
        // path names six. Where they disagree the *coordinate* instrument is the
        // one at fault, for `M-374`'s reason: an apex that lands exactly on a
        // sample plane is mis-filed as a face vertex. Printed rather than
        // asserted, because the disagreement is a fact about the lattice and the
        // verdicts are taken from the predicted column.
        let mut divergent = Vec::new();
        for s in std::iter::once(&control).chain(sweeps.iter()) {
            let predicted = s.apexes + 6 * s.hexagon_cells;
            if predicted != s.interior_vertices {
                divergent.push(format!(
                    "{}/{}@{}^3 {predicted} vs {}",
                    s.arm.name, s.arm.seed, s.size, s.interior_vertices
                ));
            }
        }
        println!(
            "\ninterior_apexes + 6*hexagon_cells vs interior_vertices: {} of {} rows agree{}",
            sweeps.len() + 1 - divergent.len(),
            sweeps.len() + 1,
            if divergent.is_empty() {
                String::new()
            } else {
                format!("; on-plane apexes at {}", divergent.join(" "))
            }
        );

        println!(
            "\nC1: post-fix unclosed meshes, per seed: {}",
            sweeps
                .iter()
                .filter(|s| s.arm.name == "generic")
                .map(|s| s.unclosed_post_fix.to_string())
                .collect::<Vec<_>>()
                .join(" ")
        );
        println!(
            "C2: pre-fix unclosed at 6^3 / elsewhere, exact arm: {} / {}",
            sweeps
                .iter()
                .filter(|s| s.arm.name == "exact" && s.size == 6)
                .map(|s| s.unclosed_pre_fix)
                .sum::<u64>(),
            sweeps
                .iter()
                .filter(|s| s.arm.name == "exact" && s.size != 6)
                .map(|s| s.unclosed_pre_fix)
                .sum::<u64>()
        );
        println!(
            "C3: unit arm interior apexes {} against exact {}",
            sweeps
                .iter()
                .filter(|s| s.arm.name == "unit")
                .map(|s| s.apexes)
                .sum::<u64>(),
            sweeps
                .iter()
                .filter(|s| s.arm.name == "exact")
                .map(|s| s.apexes)
                .sum::<u64>()
        );

        run.record(&control.row(READING));
        for s in &sweeps {
            run.record(&s.row(READING));
        }
    });
}
