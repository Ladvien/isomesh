//! E-318 — a triply periodic minimal surface has an exact Euler characteristic,
//! and whether you measure it is decided by the grid and by the seam.
//!
//! ```bash
//! cd bevy_isomesh && cargo run --example tpms_euler --release
//! ```
//!
//! Keys: `1` `2` `3` gyroid / Schwarz P / Schwarz D, `[` `]` periods per axis,
//! `,` `.` voxels per period, `V` swaps the periodic wrap for the open control
//! arm, `W` wireframe, `N` normals, `G` the extraction box, `H` HUD.
//!
//! Under `ISOMESH_CAPTURE` it needs no keyboard: three surfaces times `N = 1, 2,
//! 3` is a nine-state cycle, one state per `ISOMESH_CAPTURE_FRAMES / 9` captured
//! frames, so the clip is the prediction tracking `N^3` across three surfaces
//! rather than a still. `ISOMESH_FIELD=0|1|2` pins the surface for a screenshot,
//! `ISOMESH_PERIODS` and `ISOMESH_VOXELS` pin the grid, `ISOMESH_OPEN` starts on
//! the control arm.
//!
//! ```bash
//! # Nine states, ten captured frames each: exactly one cycle, so the GIF loops.
//! ISOMESH_CAPTURE_FRAMES=90 ISOMESH_CAPTURE_EVERY=2 \
//!   scripts/record_gif.sh tpms_euler docs/gifs/e318.gif
//! # The same cycle on the degenerate grid, where Schwarz D misses by its pinches.
//! ISOMESH_VOXELS=32 ISOMESH_CAPTURE_FRAMES=90 ISOMESH_CAPTURE_EVERY=2 \
//!   scripts/record_gif.sh tpms_euler docs/gifs/e318-degenerate.gif
//! ```
//!
//! # What is on screen
//!
//! One nodal triply periodic minimal surface over `N` periods per axis, on the
//! domain `[0, 2*pi*N]^3`. Beside it the HUD carries the **predicted** Euler
//! characteristic — `-8*N^3` for the gyroid, `-4*N^3` for Schwarz P, `-16*N^3`
//! for Schwarz D — and the **measured** one, from
//! [`isomesh::validate::validate_features`] over the extracted mesh.
//!
//! Red lines are boundary edges, yellow markers are non-manifold edges. On the
//! periodic arm there should be none of either; every one you can see is
//! accounted for in the HUD.
//!
//! # What the number means
//!
//! All three surfaces have **genus 3 in their own primitive translational
//! cell**, so `chi = 2 - 2*3 = -4` there. The `-8 / -4 / -16` spread is entirely
//! about how many primitive cells tile the conventional cubic cell of side one
//! period: two for the gyroid, one for Schwarz P, four for Schwarz D. `chi` is
//! additive over a disjoint decomposition, so `N^3` cubic cells give
//! `N^3 * chi_cell`.
//!
//! **The demo does not transcribe that count, it measures it.** Each rebuild
//! evaluates all seven non-zero half-period shifts on a 17^3 offset grid and
//! classifies each one as leaving the nodal function invariant, negating it, or
//! neither. The invariant ones are the extra centring translations, so
//! `primitive cells = 1 + invariant shifts` and `chi_cell = -4 * that`. The
//! prediction on screen is the end of that chain, which is why it cannot drift
//! away from its own justification. That census is P-143 C2, and C2 **held**:
//! residuals of `1.9e-15`, `1.4e-15` and `1.7e-15` against a `1e-12` threshold.
//!
//! # Schwarz P is not body-centred, and the demo says so
//!
//! `Im-3m` is a body-centred space-group symbol, and reading the translation
//! lattice off it is wrong here. Substituting `(x,y,z) -> (x+pi, y+pi, z+pi)`
//! into `F_P = cos x + cos y + cos z` flips **one** sign per term and so
//! **negates** `F_P` rather than leaving it alone. A negating shift maps the
//! zero set to itself but exchanges the two labyrinths, so it is a symmetry of
//! the surface and not a translation of the labelled structure, and it does not
//! shrink the translational cell. Schwarz P's conventional cubic cell therefore
//! holds **one** primitive cell, its translation lattice is `simple_cubic`, and
//! its `chi` per cubic cell is `-4` rather than `-8`. The HUD prints the shift's
//! measured verdict beside the lattice name for exactly this reason.
//!
//! The gyroid is the opposite case: every term of `F_G` is a product of two trig
//! factors, both flip, and the function is **invariant** — one extra centring
//! translation, `bcc`, two primitive cells, `-8`. Schwarz D's four terms are
//! products of three factors each, so `(pi,pi,pi)` negates it too, while the
//! three half **face** diagonals leave it invariant: `fcc`, four primitive
//! cells, `-16`.
//!
//! # Use an odd voxels-per-period
//!
//! The ladder is `32`, `33`, `65`, `97`, and `33` is the default. **`32` is in
//! the ladder as a control**, because a demo that only offers grids which work
//! cannot show the failure it is warning about.
//!
//! Schwarz D gives the wrong `chi` at every multiple of 8. R-142's harness
//! (`crates/isomesh/benches/common/tpms.rs`) swept 168 configurations and read
//! `-12` instead of `-16` at 32 and 56, `-9` at 64, `-7` at 96 and `+1` at 128,
//! and only there. Those grids put samples on the `pi/4` lattice, where `F_D`'s
//! four terms are equal in magnitude and cancel to **exactly** `0.0` (at
//! `(pi/4, pi/4, 3*pi/4)`, for instance). A sample on the isosurface puts the
//! crossing parameter at 0 or 1, two cut edges of one cell place coincident
//! vertices, and the weld turns them into a pinch. That is M-48's mechanism, not
//! the wrap's, and the fix is the grid: pick an odd `voxels_per_period`.
//!
//! **In all twelve pinching runs of that sweep, `chi_measured - chi_predicted`
//! equalled `non_manifold_edges` exactly** — never off by one, because each
//! pinch merges two sheets and costs precisely one from `chi`. The HUD prints
//! that subtraction beside that count, and it is the line that turns a wrong
//! number into a named mechanism. Two of those runs are in the committed CSVs
//! and this demo reaches both: Schwarz D at 32 voxels (gap 4, four pinches,
//! p-143.csv's `degenerate_grid_control` row) and Schwarz P at `N = 3`, 33
//! voxels (gap 6, six pinches, from ordinary floating-point cancellation rather
//! than an exact lattice).
//!
//! # The control arm is judged by `boundary_edges`, never by `chi`
//!
//! `V` drops the seam identification. The picture does not change — identifying
//! opposite faces is a **combinatorial** operation and moves no vertex — but the
//! box now cuts the surface and leaves hundreds to thousands of boundary edges,
//! and `chi` is no longer the surface's.
//!
//! It is tempting to recognise that arm by its `chi` disagreeing. **That does
//! not work, and P-142 measured why**: the non-wrapped gyroid reads `-3` against
//! `-8` and non-wrapped Schwarz D `-11` against `-16`, but non-wrapped Schwarz P
//! reads `-4` at `N = 1` and `-32` at `N = 2` — its own prediction, by
//! coincidence of the caps the box cuts. A vacuity control comparing only `chi`
//! passes the wrong arm on one field in three. So the verdict here refuses to
//! say `AGREES` while `boundary_edges` is non-zero, whatever `chi` reads.
//!
//! # The fields are example-local
//!
//! [`Tpms`] and [`NodalTpms`] live in this file. They are **not**
//! `isomesh::fields` reference fields, and adding them would be a much larger
//! change than it looks: 27 rows of `golden_hashes.json` and gated prose counts
//! in twelve documents. That is a Phase 28 ticket, and it is also exactly why
//! **P-143 C1 was FALSIFIED** — C1 asked for both surfaces to be *added as
//! reference fields* and to reproduce `N^3 * chi_cell`. They reproduce it — nine
//! of p-143.csv's ten in-scope readings are exact and the one that misses is the
//! Schwarz P pinch above — but they were not added.
//!
//! The nodal value is returned **directly**, so this is a level-set function and
//! not a signed distance field: `|grad F|` vanishes on the whole singular
//! skeleton. `isomesh::validate::accuracy` and anything else that reads a field
//! value as a distance is meaningless here. `chi` needs only the sign, which is
//! why it is the invariant this demo gates on. The gradient is analytic rather
//! than differenced, because a private field that implements only `sample` pays
//! seven evaluations per normal and loses the exact one (M-196).
//!
//! # What this shows, and what it does not
//!
//! - **P-142.** `chi = -8*N^3` for the gyroid under periodic wrap, at `N = 1, 2,
//!   3` and at 33, 65 and 97 voxels per period. C1 and C2 held. **C3 was
//!   FALSIFIED**: the gate is not implementable inside the existing validity
//!   suite, because the crate has no periodic-wrap extraction and the oracle
//!   cannot join without a trait-signature change — which is exactly what C3's
//!   own falsifier anticipated. The seam identification in this file is
//!   example-local for the same reason.
//! - **P-143.** Three surfaces, three different predictions, all reproduced on a
//!   conforming odd grid. **C1 FALSIFIED** (not added as reference fields, see
//!   above); C2 held, and the shift census on screen is C2 recomputed live.
//! - **P-144.** The neighbouring question — can a periodic *noise* field carry a
//!   `chi` oracle at all — is cited rather than drawn, because this demo has no
//!   noise arm. Its result is the interesting one: **at one octave every seed
//!   converges**, which the registration did not expect (its falsifier called
//!   continued movement at the top rung "the likely outcome"). All 35 wrapped
//!   one-octave rows carry `boundary_edges = 0` and `oracle_exists = true`, and
//!   the rung it settles on is 17, 33 or 49 depending on the seed. At three
//!   octaves there is no oracle on any of the 35 rows and four of the five seeds
//!   never converge at all. C2 held — the converged values are `-8`, `-8`, `-8`,
//!   `-12` and `-16`, a span of 8 to 10, so any such gate is per-seed — and C1
//!   reads `false`.
//!
//! And what it does not show: nothing here says the *shipped* Marching Cubes is
//! defective. Phase 27 measured zero non-manifold edges, zero non-manifold
//! vertices and zero self-intersecting pairs on all eight reference fields. The
//! pinches this demo can produce need a sample that lands exactly on the
//! isosurface, and that needs a field and a grid built to arrange it.

mod common;

use std::f64::consts::PI;

use bevy::platform::time::Instant;
use bevy::prelude::*;
use bevy_isomesh::MeshBuilder;
use common::{Capture, CommonPlugin, DemoDomain, DemoMesh, DemoStats, OrbitCamera, ViewFlags};
use isomesh::marching_cubes::MarchingCubes;
use isomesh::validate::{MeshReport, ValidateConfig, validate_features};
use isomesh::weld::{Welder, epsilon_for};
use isomesh::{MeshBuffer, MeshSink, RuntimeShape3, Sdf};

// ── the fields, example-local (see the module docs) ─────────────────────────

/// The three nodal triply periodic minimal surfaces this demo measures.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Tpms {
    /// Schoen's gyroid, `Ia-3d`.
    Gyroid,
    /// Schwarz' primitive surface, `Im-3m`.
    SchwarzP,
    /// Schwarz' diamond surface, `Pn-3m`.
    SchwarzD,
}

impl Tpms {
    /// All three, in the order P-142 and P-143 report them. The digit keys index
    /// this list.
    const ALL: [Tpms; 3] = [Tpms::Gyroid, Tpms::SchwarzP, Tpms::SchwarzD];

    /// The CSV `field` column: `gyroid`, `schwarz_p`, `schwarz_d`.
    fn name(self) -> &'static str {
        match self {
            Tpms::Gyroid => "gyroid",
            Tpms::SchwarzP => "schwarz_p",
            Tpms::SchwarzD => "schwarz_d",
        }
    }

    /// The crystallographic space group of the surface.
    ///
    /// Printed beside the *measured* lattice on purpose: `Im-3m` is a
    /// body-centred symbol and Schwarz P's translation lattice is not
    /// body-centred, which is the whole of P-143 C2.
    fn space_group(self) -> &'static str {
        match self {
            Tpms::Gyroid => "Ia-3d",
            Tpms::SchwarzP => "Im-3m",
            Tpms::SchwarzD => "Pn-3m",
        }
    }

    /// The nodal expression, for the HUD. ASCII, because the HUD font renders
    /// nothing else.
    fn expression(self) -> &'static str {
        match self {
            Tpms::Gyroid => "sin x cos y + sin y cos z + sin z cos x",
            Tpms::SchwarzP => "cos x + cos y + cos z",
            Tpms::SchwarzD => "sx sy sz + sx cy cz + cx sy cz + cx cy sz",
        }
    }
}

/// The nodal function of `kind` at `p`, with one period equal to `2*pi` per
/// axis.
///
/// Negative inside one labyrinth, positive in the other. Not a distance — see
/// the module docs.
fn nodal(kind: Tpms, p: [f64; 3]) -> f64 {
    let (sx, cx) = (p[0].sin(), p[0].cos());
    let (sy, cy) = (p[1].sin(), p[1].cos());
    let (sz, cz) = (p[2].sin(), p[2].cos());
    match kind {
        Tpms::Gyroid => sx * cy + sy * cz + sz * cx,
        Tpms::SchwarzP => cx + cy + cz,
        Tpms::SchwarzD => sx * sy * sz + sx * cy * cz + cx * sy * cz + cx * cy * sz,
    }
}

/// The analytic gradient of [`nodal`].
///
/// Overriding [`Sdf::gradient`] rather than letting it central-difference: the
/// default costs six extra [`nodal`] calls per normal and loses the exact
/// direction, which is M-196's lesson about private field copies. It changes no
/// count on screen — Marching Cubes places vertices from sample *values* — so
/// the readings still match the committed CSVs to the digit.
fn nodal_gradient(kind: Tpms, p: [f64; 3]) -> [f64; 3] {
    let (sx, cx) = (p[0].sin(), p[0].cos());
    let (sy, cy) = (p[1].sin(), p[1].cos());
    let (sz, cz) = (p[2].sin(), p[2].cos());
    match kind {
        Tpms::Gyroid => [cx * cy - sz * sx, cy * cz - sx * sy, cz * cx - sy * sz],
        Tpms::SchwarzP => [-sx, -sy, -sz],
        Tpms::SchwarzD => [
            cx * sy * sz + cx * cy * cz - sx * sy * cz - sx * cy * sz,
            sx * cy * sz - sx * sy * cz + cx * cy * cz - cx * sy * sz,
            sx * sy * cz - sx * cy * sz - cx * sy * sz + cx * cy * cz,
        ],
    }
}

/// Samples per axis in the deterministic grid the shift identities are checked
/// on. `17^3 = 4913` points, the same grid P-143 used.
const SHIFT_GRID: u32 = 17;

/// Residual below which a shift identity is called exact.
///
/// The nodal functions are sums of at most four products of unit-magnitude
/// terms, so an exact identity lands within a few ulp of `4.0` — under `1e-15` —
/// and a broken one is `O(1)`. `1e-12` is three orders clear of both.
const SHIFT_RESIDUAL_TOLERANCE: f64 = 1e-12;

/// `(max |F(p + shift) - F(p)|, max |F(p + shift) + F(p)|)` over a fixed grid.
///
/// The first component is the residual of the *invariance* claim, the second of
/// the *negation* claim. Reporting both is what makes either non-vacuous: a grid
/// on which `F` happened to vanish returns two zeros, and a reader can see that
/// instead of being told "invariant" twice.
///
/// The grid is offset by `1/3`, `1/5`, `1/7` per axis so that no sample sits on
/// a symmetry plane: a coordinate is a multiple of `pi/2` only when
/// `4*(i + o_k)/17` is an integer, and `4/3`, `4/5` and `4/7` are not. On a
/// naive grid a large fraction of samples land where both residuals collapse and
/// the check stops discriminating.
fn shift_residuals(kind: Tpms, shift: [f64; 3]) -> (f64, f64) {
    let offsets = [1.0 / 3.0, 1.0 / 5.0, 1.0 / 7.0];
    let coord =
        |i: u32, axis: usize| 2.0 * PI * (f64::from(i) + offsets[axis]) / f64::from(SHIFT_GRID);
    let mut symmetric = 0.0f64;
    let mut antisymmetric = 0.0f64;
    for i in 0..SHIFT_GRID {
        let x = coord(i, 0);
        for j in 0..SHIFT_GRID {
            let y = coord(j, 1);
            for k in 0..SHIFT_GRID {
                let p = [x, y, coord(k, 2)];
                let here = nodal(kind, p);
                let there = nodal(kind, [p[0] + shift[0], p[1] + shift[1], p[2] + shift[2]]);
                symmetric = symmetric.max((there - here).abs());
                antisymmetric = antisymmetric.max((there + here).abs());
            }
        }
    }
    (symmetric, antisymmetric)
}

/// What a half-period shift does to the nodal function.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum ShiftClass {
    /// `F(p + t) = F(p)`. A genuine extra centring translation, so it shrinks the
    /// primitive cell.
    Invariant,
    /// `F(p + t) = -F(p)`. Maps the zero set to itself but exchanges the two
    /// labyrinths, so it is a symmetry of the surface and **not** a translation
    /// of the labelled structure.
    Negated,
    /// Neither, within [`SHIFT_RESIDUAL_TOLERANCE`].
    Neither,
}

/// Classify one shift, requiring the claimed relation to hold **and** the
/// opposite one to fail.
///
/// The second half is the vacuity guard: `F = 0` satisfies invariance and
/// negation at once, so a residual of zero on its own establishes nothing.
fn classify(kind: Tpms, shift: [f64; 3]) -> ShiftClass {
    let (symmetric, antisymmetric) = shift_residuals(kind, shift);
    let exact = |v: f64| v <= SHIFT_RESIDUAL_TOLERANCE;
    match (exact(symmetric), exact(antisymmetric)) {
        (true, false) => ShiftClass::Invariant,
        (false, true) => ShiftClass::Negated,
        _ => ShiftClass::Neither,
    }
}

/// The census of all seven non-zero half-period shifts, and what it implies.
///
/// This is the demo's prediction, derived rather than transcribed: the invariant
/// shifts are the cosets of the simple cubic lattice inside the surface's own
/// translation lattice, so their count plus one is the number of primitive cells
/// in the conventional cubic cell, and `chi_cell` is `-4` times that.
#[derive(Clone, Copy, Debug)]
struct Symmetry {
    /// Shifts leaving `F` invariant.
    invariant: u32,
    /// Shifts negating `F`.
    negated: u32,
    /// Shifts doing neither.
    neither: u32,
    /// What `(pi, pi, pi)` does — the body-centring operation itself.
    body_centring: ShiftClass,
    /// `max |F(p + (pi,pi,pi)) - F(p)|` over the census grid.
    body_symmetric: f64,
    /// `max |F(p + (pi,pi,pi)) + F(p)|` over the same grid.
    body_antisymmetric: f64,
}

impl Symmetry {
    /// Measure all seven shifts of `kind`.
    ///
    /// The shifts are the non-zero elements of `{0, pi}^3`, enumerated by bit
    /// pattern with bit `i` meaning axis `i` — the same corner-bit convention
    /// the extractor and `common::draw_domain` use.
    fn measure(kind: Tpms) -> Self {
        let half = |bit: u32| if bit == 0 { 0.0 } else { PI };
        let mut census = Self {
            invariant: 0,
            negated: 0,
            neither: 0,
            body_centring: ShiftClass::Neither,
            body_symmetric: 0.0,
            body_antisymmetric: 0.0,
        };
        for mask in 1..8u32 {
            let shift = [half(mask & 1), half(mask & 2), half(mask & 4)];
            match classify(kind, shift) {
                ShiftClass::Invariant => census.invariant += 1,
                ShiftClass::Negated => census.negated += 1,
                ShiftClass::Neither => census.neither += 1,
            }
        }
        let body = [PI, PI, PI];
        let (symmetric, antisymmetric) = shift_residuals(kind, body);
        census.body_centring = classify(kind, body);
        census.body_symmetric = symmetric;
        census.body_antisymmetric = antisymmetric;
        census
    }

    /// Primitive translational cells per conventional cubic cell: the identity
    /// plus every invariant centring shift.
    fn primitive_cells(self) -> i64 {
        1 + i64::from(self.invariant)
    }

    /// `chi` inside one conventional cubic cell: genus 3 per primitive cell
    /// gives `-4` there, and `chi` is additive.
    fn chi_per_cubic_cell(self) -> i64 {
        -4 * self.primitive_cells()
    }

    /// The Bravais lattice of the **translation** group, named from the census
    /// rather than read off the space-group symbol.
    fn lattice(self) -> &'static str {
        match self.invariant {
            0 => "simple_cubic",
            1 => "bcc",
            3 => "fcc",
            _ => "unclassified",
        }
    }
}

/// A nodal TPMS over `periods` periods per axis on `[0, 2*pi*periods]^3`.
///
/// The domain is exactly `periods^3` conventional cubic cells. `periods` comes
/// from the [`PERIODS`] ladder and is never zero, which would give an empty
/// domain and a predicted `chi` of zero — a number that would look like a
/// measurement.
#[derive(Clone, Copy, Debug)]
struct NodalTpms {
    /// Which surface.
    kind: Tpms,
    /// Periods per axis. The domain grows with it; the function does not scale.
    periods: u32,
}

impl NodalTpms {
    /// `(lo, hi)` of the extraction box: exactly `periods` full periods per axis,
    /// from the origin.
    fn domain(self) -> ([f64; 3], [f64; 3]) {
        let span = 2.0 * PI * f64::from(self.periods);
        ([0.0; 3], [span; 3])
    }

    /// The **periodic-conforming** sample grid: `(shape, origin, cell_size)`.
    ///
    /// `voxels_per_period` cells span one period, so `cell_size` divides the
    /// period exactly and the sample at `hi` is the same period point as the
    /// sample at `lo`. That identity is the only reason the two opposite faces
    /// carry the same sign configuration and hence the same cut, and so the only
    /// reason [`identify`] can close the mesh. A grid whose spacing does not
    /// divide the period cannot be wrapped however the seam is welded, and P-142
    /// names exactly that as the defect the project mistook for "the gyroid has
    /// no chi".
    ///
    /// Sample count is `voxels_per_period * periods + 1` per axis, because
    /// `Shape3::size` counts samples and `n` samples span `n - 1` cells.
    fn periodic_grid(self, voxels_per_period: u32) -> Option<(RuntimeShape3, [f64; 3], f64)> {
        let samples = voxels_per_period * self.periods + 1;
        let shape = RuntimeShape3::new([samples; 3]).ok()?;
        let cell_size = 2.0 * PI / f64::from(voxels_per_period);
        Some((shape, self.domain().0, cell_size))
    }
}

impl Sdf for NodalTpms {
    type Scalar = f64;

    fn sample(&self, p: [f64; 3]) -> f64 {
        nodal(self.kind, p)
    }

    fn gradient(&self, p: [f64; 3]) -> [f64; 3] {
        nodal_gradient(self.kind, p)
    }
}

// ── the seam identification, example-local ──────────────────────────────────

/// Whether opposite boundary faces of the extraction are identified.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Wrap {
    /// Identified by the period translation: the mesh closes on the 3-torus and
    /// `chi` is the surface's own.
    Periodic,
    /// The control arm. The box cuts the surface and leaves a boundary.
    Open,
}

impl Wrap {
    /// The other arm — what `V` switches to.
    fn other(self) -> Self {
        match self {
            Wrap::Periodic => Wrap::Open,
            Wrap::Open => Wrap::Periodic,
        }
    }

    /// The CSV `wrap_mode` word.
    fn name(self) -> &'static str {
        match self {
            Wrap::Periodic => "periodic",
            Wrap::Open => "open",
        }
    }
}

/// One position with every coordinate within `tol` of the far face folded onto
/// the near one.
///
/// Written per axis, so an edge of the box folds twice and a corner three times.
/// A single "is this the `x = hi` face" test would miss those.
fn fold_to_lo(p: [f64; 3], lo: [f64; 3], hi: [f64; 3], tol: f64) -> [f64; 3] {
    let mut folded = p;
    for ((slot, &l), &h) in folded.iter_mut().zip(lo.iter()).zip(hi.iter()) {
        if (h - *slot).abs() <= tol {
            *slot = l;
        }
    }
    folded
}

/// Weld `mesh` in place, folding the far faces onto the near ones first when
/// `wrap` is [`Wrap::Periodic`]. Returns the vertices the weld removed.
///
/// **One weld, both arms.** The open arm needs one too: Marching Cubes shares a
/// vertex only between cells meeting on a grid *edge*, so a sample landing on
/// the isosurface leaves genuinely coincident duplicates (M-48), and a `chi`
/// taken over an unwelded buffer counts one surface as several. Folding only
/// changes which vertices the *same* weld considers coincident, so the two arms
/// are read by one instrument and cannot disagree about what "coincident" means.
///
/// The fold lives in a key buffer rather than in the mesh: folding the mesh
/// would move a boundary vertex a whole period away from the triangle that owns
/// it, which is a different mesh, not a wrapped one. The surviving buffer keeps
/// each representative's **original** position and normal.
///
/// The welder is the crate's shipped [`Welder`] rather than a local copy, so its
/// determinism argument — a sorted broadphase, lowest-indexed representative, no
/// hash iteration — is the one already stated in `isomesh::weld`.
///
/// After a periodic identification the buffer is a valid *simplicial complex*
/// and an invalid *geometric mesh*: some triangles now have one corner at one
/// side of the box and another at the far side. That is not a defect, it is what
/// "closes on the torus" means — and it is why the picture is built from the
/// un-identified extraction while the counting is done here.
fn identify(
    mesh: &mut MeshBuffer<f64>,
    lo: [f64; 3],
    hi: [f64; 3],
    tol: f64,
    wrap: Wrap,
) -> Option<u64> {
    let count = mesh.positions.len();
    let mut keys = MeshBuffer::<f64>::new();
    keys.positions = mesh
        .positions
        .iter()
        .map(|p| match wrap {
            Wrap::Periodic => fold_to_lo(*p, lo, hi, tol),
            Wrap::Open => *p,
        })
        .collect();
    // The welder reads normals only to compact them, and this buffer's are
    // discarded; the real ones travel with the representatives below.
    keys.normals = vec![[0.0; 3]; count];
    keys.indices = mesh.indices.clone();

    let mut welder = Welder::<f64>::new();
    let report = welder.weld(&mut keys, tol).ok()?;
    let remap = welder.remap();

    let survivors = keys.positions.len();
    let mut positions = vec![[0.0f64; 3]; survivors];
    let mut normals = vec![[0.0f64; 3]; survivors];
    let mut written = vec![false; survivors];
    // Ascending input order, and the welder gives a representative an output
    // index no greater than its input index, so the first vertex to reach a
    // given output slot is that slot's representative.
    for (input, &output) in remap.iter().enumerate() {
        let output = output as usize;
        if !written[output] {
            written[output] = true;
            positions[output] = mesh.positions[input];
            normals[output] = mesh.normals[input];
        }
    }

    mesh.positions = positions;
    mesh.normals = normals;
    mesh.indices = keys.indices;
    Some(report.vertices_removed() as u64)
}

// ── the demo ────────────────────────────────────────────────────────────────

/// Voxels per period the demo offers.
///
/// `32` is the degenerate control and is deliberately first: see the module
/// docs. The rest are odd on purpose.
const VOXELS: [u32; 4] = [32, 33, 65, 97];

/// Index of the default rung, `33` — the smallest odd grid P-142 measured, and
/// one that fits every period count in [`PERIODS`].
const DEFAULT_VOXELS: usize = 1;

/// Periods per axis the demo offers.
const PERIODS: [u32; 3] = [1, 2, 3];

/// Largest sample count per axis a rung may ask for.
///
/// `100^3` is a million samples and about 285k triangles on the gyroid, which is
/// a visible hitch on a rebuild and no worse. It is also exactly the ceiling
/// P-142 and P-143 ran under, so every configuration this demo can reach has a
/// committed row to be checked against.
const MAX_SAMPLES_PER_AXIS: u32 = 100;

/// The finest rung of [`VOXELS`] at or below `wanted` that fits
/// [`MAX_SAMPLES_PER_AXIS`] at `periods` periods.
///
/// Rung 0 always fits — 32 voxels over three periods is 97 samples per axis — so
/// this always returns a usable index. Raising `N` therefore drops the grid
/// rather than asking for an extraction that would not finish, and the HUD
/// prints the rung that was actually used.
fn affordable(periods: u32, wanted: usize) -> usize {
    let mut index = wanted.min(VOXELS.len() - 1);
    while index > 0 && VOXELS[index] * periods + 1 > MAX_SAMPLES_PER_AXIS {
        index -= 1;
    }
    index
}

/// Captured frames spent on each of the nine `(surface, periods)` states.
///
/// Read from the environment rather than from [`Capture`], which keeps its
/// length private. `ISOMESH_CAPTURE_FRAMES=90` gives ten frames a state and one
/// exact cycle, so the clip loops.
fn capture_frames_per_state() -> u32 {
    let frames: u32 = std::env::var("ISOMESH_CAPTURE_FRAMES")
        .ok()
        .and_then(|v| v.parse().ok())
        .unwrap_or(60);
    let states = (Tpms::ALL.len() * PERIODS.len()) as u32;
    (frames / states).max(1)
}

/// What the reader has asked for.
#[derive(Resource)]
struct Demo {
    /// Index into [`Tpms::ALL`]. Follows `ViewFlags::field`, so the digit keys
    /// and `ISOMESH_FIELD` both drive it.
    surface: usize,
    /// Index into [`PERIODS`].
    periods: usize,
    /// Index into [`VOXELS`] — the rung *asked for*. The rung used is
    /// [`affordable`] of this and the period count.
    voxels: usize,
    /// Whether the seam is identified.
    wrap: Wrap,
}

impl Default for Demo {
    fn default() -> Self {
        Self {
            surface: 0,
            periods: 0,
            voxels: DEFAULT_VOXELS,
            wrap: Wrap::Periodic,
        }
    }
}

impl Demo {
    /// Periods per axis.
    fn period_count(&self) -> u32 {
        PERIODS[self.periods.min(PERIODS.len() - 1)]
    }

    /// Voxels per period, after [`affordable`].
    fn voxel_count(&self) -> u32 {
        VOXELS[affordable(self.period_count(), self.voxels)]
    }

    /// The surface.
    fn kind(&self) -> Tpms {
        Tpms::ALL[self.surface.min(Tpms::ALL.len() - 1)]
    }

    /// Everything a rebuild depends on. A rebuild that is not keyed on exactly
    /// this either misses a change or repeats a second of work every frame.
    fn rung(&self) -> (Tpms, u32, u32, Wrap) {
        (
            self.kind(),
            self.period_count(),
            self.voxel_count(),
            self.wrap,
        )
    }
}

/// What one rebuild measured. Written by [`rebuild`], read by [`report`] and
/// [`draw_defects`].
#[derive(Resource, Default)]
struct Measurement {
    /// `None` until the first rebuild has run.
    reading: Option<Reading>,
}

/// One complete reading, HUD-ready.
struct Reading {
    /// Which surface.
    surface: Tpms,
    /// Periods per axis.
    periods: u32,
    /// Voxels per period.
    voxels: u32,
    /// Samples per axis.
    samples: u32,
    /// Grid spacing.
    cell_size: f64,
    /// Which arm.
    wrap: Wrap,
    /// The measured symmetry census the prediction is derived from.
    symmetry: Symmetry,
    /// `periods^3 * symmetry.chi_per_cubic_cell()`.
    chi_predicted: i64,
    /// What `validate_features` read off the identified complex.
    report: MeshReport,
    /// Vertices the identification removed.
    merged: u64,
    /// Vertices Marching Cubes emitted, before identification.
    mesh_vertices: usize,
    /// Triangles Marching Cubes emitted.
    mesh_triangles: usize,
    /// Extraction time.
    extract_ms: f64,
    /// Identification plus validation time.
    measure_ms: f64,
    /// Boundary-edge segments, centred on the origin like the mesh.
    boundary: Vec<[Vec3; 2]>,
    /// Midpoints of the non-manifold edges — the pinches.
    pinches: Vec<Vec3>,
    /// Marker radius for a pinch, in world units.
    pinch_radius: f32,
}

impl Reading {
    /// Half the domain span. The mesh is centred, so this is the box.
    fn half_extent(&self) -> f32 {
        (PI * f64::from(self.periods)) as f32
    }

    /// `chi_measured - chi_predicted`.
    fn gap(&self) -> i64 {
        self.report.euler_characteristic - self.chi_predicted
    }

    /// Whether the identified complex is closed. The one question that has to be
    /// answered before `chi` is worth reading at all.
    fn closed(&self) -> bool {
        self.report.boundary_edges == 0
    }
}

/// The result of one measurement, plus the mesh the reader looks at.
struct Measured {
    /// The numbers.
    reading: Reading,
    /// The **un-identified** extraction, centred on the origin. Identifying the
    /// seam is combinatorial and would stretch triangles across the box, so the
    /// picture is built before it — see [`identify`].
    mesh: Mesh,
}

/// A `[f64; 3]` sample position as a centred Bevy vector.
fn centred(p: [f64; 3], centre: f64) -> Vec3 {
    Vec3::new(
        (p[0] - centre) as f32,
        (p[1] - centre) as f32,
        (p[2] - centre) as f32,
    )
}

/// Extract, identify the seam, validate.
///
/// `None` when any step refuses the configuration; the caller logs it. `.ok()?`
/// rather than `unwrap`, which is this crate's rule — an example that unwraps
/// teaches unwrapping.
fn measure(surface: Tpms, periods: u32, voxels: u32, wrap: Wrap) -> Option<Measured> {
    let field = NodalTpms {
        kind: surface,
        periods,
    };
    let (lo, hi) = field.domain();
    let (shape, origin, cell_size) = field.periodic_grid(voxels)?;

    let mut buffer = MeshBuffer::<f64>::new();
    let started = Instant::now();
    MarchingCubes::<f64>::new()
        .extract(&field, &shape, origin, cell_size, &mut buffer)
        .ok()?;
    let extract_ms = started.elapsed().as_secs_f64() * 1000.0;
    // An empty extraction would leave `Assets<Mesh>` logging
    // "attempted to copy element data for an unallocated key" twice a frame
    // forever (E-307). Unreachable on these fields at these grids, and cheaper
    // to refuse than to diagnose.
    if buffer.indices.is_empty() {
        return None;
    }

    let mesh_vertices = buffer.positions.len();
    let mesh_triangles = buffer.indices.len() / 3;
    let centre = PI * f64::from(periods);

    let mut builder = MeshBuilder::new();
    for (p, n) in buffer.positions.iter().zip(buffer.normals.iter()) {
        let c = centred(*p, centre);
        builder.vertex([c.x, c.y, c.z], [n[0] as f32, n[1] as f32, n[2] as f32]);
    }
    for tri in buffer.indices.as_chunks::<3>().0 {
        builder.triangle(tri[0], tri[1], tri[2]);
    }
    let mesh = builder.into_mesh();

    let measuring = Instant::now();
    let tol = epsilon_for(cell_size);
    let mut counted = buffer;
    let merged = identify(&mut counted, lo, hi, tol, wrap)?;
    let cfg = ValidateConfig::from_cell_size(cell_size).ok()?;
    let (report, features) = validate_features(&counted.positions, &counted.indices, &cfg);
    let measure_ms = measuring.elapsed().as_secs_f64() * 1000.0;

    let segment = |edge: &[u32; 2]| {
        [
            centred(counted.positions[edge[0] as usize], centre),
            centred(counted.positions[edge[1] as usize], centre),
        ]
    };
    let boundary: Vec<[Vec3; 2]> = features.boundary_edges.iter().map(segment).collect();
    let pinches: Vec<Vec3> = features
        .edges
        .iter()
        .map(|edge| {
            let [a, b] = segment(edge);
            a.midpoint(b)
        })
        .collect();

    let symmetry = Symmetry::measure(surface);
    let n = i64::from(periods);
    Some(Measured {
        reading: Reading {
            surface,
            periods,
            voxels,
            samples: voxels * periods + 1,
            cell_size,
            wrap,
            symmetry,
            chi_predicted: symmetry.chi_per_cubic_cell() * n * n * n,
            report,
            merged,
            mesh_vertices,
            mesh_triangles,
            extract_ms,
            measure_ms,
            boundary,
            pinches,
            // Three cells across. A pinch edge is one cell long and invisible at
            // the framing this demo uses; the marker is what makes the count on
            // the HUD findable in the picture.
            pinch_radius: (cell_size * 1.5) as f32,
        },
        mesh,
    })
}

/// Defect overlay: boundary edges.
#[derive(Default, Reflect, GizmoConfigGroup)]
struct BoundaryGizmos;

/// Pinch markers. Their own group, biased harder, so a marker is never lost
/// behind the surface it sits on — `manifold_check` and `critical_cells` earned
/// this the hard way.
#[derive(Default, Reflect, GizmoConfigGroup)]
struct PinchGizmos;

fn main() {
    App::new()
        .add_plugins(DefaultPlugins.set(WindowPlugin {
            primary_window: Some(Window {
                title: "isomesh - E-318 tpms euler".into(),
                // Web only, inert on native: bind to the 1280x720 canvas the
                // page supplies rather than letting Bevy append its own. The HUD
                // panel is laid out in pixels for that size.
                canvas: Some("#isomesh-canvas".into()),
                ..default()
            }),
            ..default()
        }))
        .add_plugins(CommonPlugin)
        .init_gizmo_group::<BoundaryGizmos>()
        .init_gizmo_group::<PinchGizmos>()
        .init_resource::<Demo>()
        .init_resource::<Measurement>()
        .add_systems(Startup, setup)
        // `PreUpdate` for E-306's reason: the harness's `update_hud` renders
        // `DemoStats` and its `capture_sequence` advances `Capture::taken`, both
        // in `Update` with no ordering against an example's own systems. In
        // `Update` the HUD would render a frame-old readout beside a current
        // picture, which on this demo means a `chi` from the previous surface.
        .add_systems(PreUpdate, (controls, rebuild, report).chain())
        .add_systems(Update, draw_defects)
        .run();
}

fn setup(
    mut commands: Commands,
    mut materials: ResMut<Assets<StandardMaterial>>,
    mut gizmo_config: ResMut<GizmoConfigStore>,
    mut flags: ResMut<ViewFlags>,
    mut demo: ResMut<Demo>,
) {
    if flags.field >= Tpms::ALL.len() {
        error!(
            "ISOMESH_FIELD={} is out of range: 0=gyroid 1=schwarz_p 2=schwarz_d",
            flags.field
        );
        flags.field = 0;
    }
    demo.surface = flags.field;
    demo.periods = rung("ISOMESH_PERIODS", &PERIODS, 0);
    demo.voxels = rung("ISOMESH_VOXELS", &VOXELS, DEFAULT_VOXELS);
    // `ISOMESH_OPEN` starts on the control arm, so the before/after pair can be
    // captured without a human pressing `V` -- the same reason the harness has
    // `ISOMESH_VIEW`.
    if std::env::var("ISOMESH_OPEN").is_ok() {
        demo.wrap = Wrap::Open;
    }

    let (boundary, _) = gizmo_config.config_mut::<BoundaryGizmos>();
    boundary.line.width = 1.8;
    boundary.depth_bias = -0.4;
    let (pinch, _) = gizmo_config.config_mut::<PinchGizmos>();
    pinch.line.width = 3.2;
    pinch.depth_bias = -0.8;

    commands.spawn((
        // Named by no asset until the first rebuild, which is what an empty
        // result actually wants (E-307).
        Mesh3d::default(),
        MeshMaterial3d(materials.add(StandardMaterial {
            base_color: Color::srgb(0.72, 0.76, 0.82),
            perceptual_roughness: 0.45,
            metallic: 0.05,
            // A TPMS is two interpenetrating labyrinths and the box cuts both
            // open. Culling back faces would show holes where the far wall is.
            double_sided: true,
            cull_mode: None,
            ..default()
        })),
        Transform::IDENTITY,
        DemoMesh,
        DemoDomain {
            min: Vec3::ZERO,
            max: Vec3::ZERO,
        },
    ));
}

/// The ladder index `name` asks for, or `fallback` with one `error!` when the
/// value is not a rung.
///
/// Called from `Startup` rather than from [`Demo::default`] so the log is up by
/// the time it can complain.
fn rung(name: &str, ladder: &[u32], fallback: usize) -> usize {
    let Ok(text) = std::env::var(name) else {
        return fallback;
    };
    match text
        .parse::<u32>()
        .ok()
        .and_then(|value| ladder.iter().position(|rung| *rung == value))
    {
        Some(index) => index,
        None => {
            error!("{name}={text} is not one of {ladder:?}");
            fallback
        }
    }
}

/// Keys, and the capture's nine-state cycle.
///
/// A pinned `ISOMESH_VOXELS` survives the capture — the cycle drives the surface
/// and the period count only, so recording the degenerate ladder is one variable
/// rather than a second code path.
fn controls(
    keys: Res<ButtonInput<KeyCode>>,
    flags: Res<ViewFlags>,
    capture: Res<Capture>,
    mut demo: ResMut<Demo>,
) {
    if capture.is_active() {
        let states = Tpms::ALL.len() * PERIODS.len();
        let state = (capture.taken / capture_frames_per_state()) as usize % states;
        let (surface, periods) = (state / PERIODS.len(), state % PERIODS.len());
        if (demo.surface, demo.periods) != (surface, periods) {
            demo.surface = surface;
            demo.periods = periods;
        }
        return;
    }

    // Written only when it differs: an unconditional write marks the resource
    // changed every frame, and a reader that trusts change detection would then
    // see a rebuild request on every one of them.
    let wanted = flags.field.min(Tpms::ALL.len() - 1);
    if demo.surface != wanted {
        demo.surface = wanted;
    }
    if keys.just_pressed(KeyCode::BracketRight) {
        demo.periods = (demo.periods + 1).min(PERIODS.len() - 1);
    }
    if keys.just_pressed(KeyCode::BracketLeft) {
        demo.periods = demo.periods.saturating_sub(1);
    }
    if keys.just_pressed(KeyCode::Period) {
        demo.voxels = (demo.voxels + 1).min(VOXELS.len() - 1);
    }
    if keys.just_pressed(KeyCode::Comma) {
        demo.voxels = demo.voxels.saturating_sub(1);
    }
    if keys.just_pressed(KeyCode::KeyV) {
        demo.wrap = demo.wrap.other();
    }
}

/// Where the subject sits relative to the middle of the frame.
///
/// **The HUD must not sit on the evidence.** Centring photographs the argument
/// with its numbers over the top of it, which is E-112's lesson. Applied in the
/// camera's own basis rather than as a world offset, so it survives a change of
/// yaw.
const SUBJECT_OFFSET: Vec2 = Vec2::new(0.20, 0.10);

/// Camera distance as a multiple of the domain's half extent.
///
/// Derived from the field's own domain, never hardcoded: a fixed radius put the
/// camera *inside* the gyroid once already, and the committed screenshot was a
/// picture of an inner wall (E-304).
const FRAMING: f32 = 3.6;

fn rebuild(
    demo: Res<Demo>,
    mut measurement: ResMut<Measurement>,
    mut meshes: ResMut<Assets<Mesh>>,
    mut query: Query<(&mut Mesh3d, &mut DemoDomain)>,
    mut camera: Query<&mut OrbitCamera>,
    mut last: Local<Option<(Tpms, u32, u32, Wrap)>>,
) {
    let rung = demo.rung();
    if *last == Some(rung) {
        return;
    }
    *last = Some(rung);
    let (surface, periods, voxels, wrap) = rung;

    let Some(measured) = measure(surface, periods, voxels, wrap) else {
        error!(
            "E-318: {} at N={periods} on {} voxels/period did not extract",
            surface.name(),
            voxels
        );
        return;
    };
    let reading = measured.reading;

    let handle = meshes.add(measured.mesh);
    let half = reading.half_extent();
    for (mut mesh, mut domain) in &mut query {
        mesh.0 = handle.clone();
        domain.min = Vec3::splat(-half);
        domain.max = Vec3::splat(half);
    }
    let radius = half * FRAMING;
    for mut orbit in &mut camera {
        let direction = Vec3::new(
            orbit.yaw.cos() * orbit.pitch.cos(),
            orbit.pitch.sin(),
            orbit.yaw.sin() * orbit.pitch.cos(),
        );
        let forward = -direction;
        let right = forward.cross(Vec3::Y).normalize_or_zero();
        let up = right.cross(forward).normalize_or_zero();
        orbit.focus = -right * (SUBJECT_OFFSET.x * radius) + up * (SUBJECT_OFFSET.y * radius);
        orbit.radius = radius;
    }

    // The HUD is the evidence and a headless capture has no HUD to read. One
    // line per rebuild, so `ISOMESH_CAPTURE` leaves the reading in the log where
    // a script can hold it against p-142.csv and p-143.csv.
    info!(
        "E-318 {} N={} {} voxels/period {}^3 {}: chi predicted {} measured {} \
         (gap {}), boundary {} non-manifold {}, V {} E {} F {}, {} merged, \
         {:.1} ms extract + {:.1} ms measure",
        reading.surface.name(),
        reading.periods,
        reading.voxels,
        reading.samples,
        reading.wrap.name(),
        reading.chi_predicted,
        reading.report.euler_characteristic,
        reading.gap(),
        reading.report.boundary_edges,
        reading.report.non_manifold_edges,
        reading.report.referenced_vertices,
        reading.report.edges,
        reading.report.faces,
        reading.merged,
        reading.extract_ms,
        reading.measure_ms,
    );
    measurement.reading = Some(reading);
}

/// The agreement verdict.
///
/// It refuses to say `AGREES` while there is a boundary, whatever `chi` reads —
/// P-142 measured non-wrapped Schwarz P at `-4` (`N = 1`) and `-32` (`N = 2`),
/// its own prediction, by coincidence of the caps the box cuts.
fn verdict(reading: &Reading) -> (&'static str, String) {
    if !reading.closed() {
        return (
            "NOT CLOSED",
            format!(
                "{} boundary edges: this arm is judged by them, never by chi",
                reading.report.boundary_edges
            ),
        );
    }
    if reading.gap() == 0 {
        ("AGREES", "prediction reproduced exactly".into())
    } else {
        (
            "MISSES",
            format!("measured chi is {:+} from the prediction", reading.gap()),
        )
    }
}

/// The pinch identity, which is the diagnostic that turns a wrong number into a
/// named mechanism.
fn identity_line(reading: &Reading) -> String {
    let gap = reading.gap();
    let pinches = reading.report.non_manifold_edges;
    let note = if !reading.closed() {
        "not applicable on the open arm: the gap is the caps the box cuts"
    } else if pinches == 0 && gap == 0 {
        "no pinch and no gap"
    } else if gap == pinches as i64 {
        "EQUAL <- each pinch merges two sheets and costs exactly one (P-143)"
    } else {
        "!! differs - pinches do not account for the gap"
    };
    format!("chi_measured - chi_predicted = {gap:<5} non-manifold edges = {pinches:<5} {note}")
}

fn report(measurement: Res<Measurement>, mut stats: ResMut<DemoStats>) {
    let Some(reading) = &measurement.reading else {
        return;
    };
    let (word, detail) = verdict(reading);
    stats.title = format!(
        "E-318  tpms euler - {}  N={}  {} voxels/period  {}",
        reading.surface.name(),
        reading.periods,
        reading.voxels,
        reading.wrap.name()
    );
    stats.vertices = reading.mesh_vertices;
    stats.triangles = reading.mesh_triangles;
    stats.extract_ms = reading.extract_ms;
    stats.banner = Some((
        format!(
            "chi  predicted {}   measured {}   {word}",
            reading.chi_predicted, reading.report.euler_characteristic
        ),
        match word {
            "AGREES" => Color::srgb(0.35, 0.92, 0.55),
            "MISSES" => Color::srgb(1.0, 0.42, 0.36),
            _ => Color::srgb(1.0, 0.80, 0.30),
        },
    ));
    stats.extra = lines(reading, &detail);
    stats.keys = Some(
        "[1-3] surface   [ ] periods   , . voxels/period   V wrap   \
         W wire  N normals  G box  H hud"
            .into(),
    );
}

/// The panel, one entry per line.
fn lines(reading: &Reading, detail: &str) -> Vec<String> {
    let census = reading.symmetry;
    let expected = |value: &str| {
        if reading.wrap == Wrap::Periodic {
            value.to_string()
        } else {
            "-".to_string()
        }
    };
    let parity = if reading.voxels.is_multiple_of(2) {
        "EVEN - degenerate control"
    } else {
        "odd"
    };
    let body = match census.body_centring {
        ShiftClass::Invariant => "INVARIANT",
        ShiftClass::Negated => "NEGATED",
        ShiftClass::Neither => "NEITHER",
    };

    vec![
        format!(
            "surface    {:<10} {:<7} lattice {:<13} {} primitive cell(s) per cubic cell",
            reading.surface.name(),
            reading.surface.space_group(),
            census.lattice(),
            census.primitive_cells()
        ),
        format!("           F = {}", reading.surface.expression()),
        format!(
            "shifts     7 half-period shifts measured: {} invariant, {} negated, {} neither",
            census.invariant, census.negated, census.neither
        ),
        format!(
            "           (pi,pi,pi) {body}  residual {:.1e}  opposite relation {:.1e}  (P-143 C2)",
            if census.body_centring == ShiftClass::Invariant {
                census.body_symmetric
            } else {
                census.body_antisymmetric
            },
            if census.body_centring == ShiftClass::Invariant {
                census.body_antisymmetric
            } else {
                census.body_symmetric
            }
        ),
        format!(
            "           chi per cubic cell = -4 x {} = {},  so predicted = N^3 x {} = {}",
            census.primitive_cells(),
            census.chi_per_cubic_cell(),
            census.chi_per_cubic_cell(),
            reading.chi_predicted
        ),
        String::new(),
        format!(
            "grid       N = {} periods   {} voxels/period ({})   {}^3 samples   h {:.6}",
            reading.periods, reading.voxels, parity, reading.samples, reading.cell_size
        ),
        format!(
            "wrap       {}   {} vertices merged   {:.1} ms extract + {:.1} ms measure",
            reading.wrap.name(),
            reading.merged,
            reading.extract_ms,
            reading.measure_ms
        ),
        String::new(),
        format!("{:>32}  {:>10}", "predicted", "measured"),
        format!(
            "chi{:>29}  {:>10}   {}",
            reading.chi_predicted, reading.report.euler_characteristic, detail
        ),
        format!(
            "boundary edges{:>18}  {:>10}",
            expected("0"),
            reading.report.boundary_edges
        ),
        format!(
            "non-manifold edges{:>14}  {:>10}",
            expected("0"),
            reading.report.non_manifold_edges
        ),
        String::new(),
        identity_line(reading),
        format!(
            "V {}  E {}  F {}   components {}   genus {}",
            reading.report.referenced_vertices,
            reading.report.edges,
            reading.report.faces,
            reading.report.components,
            reading
                .report
                .genus
                .map_or_else(|| "none".to_string(), |g| g.to_string())
        ),
        String::new(),
        "P-142  gyroid = -8N^3 held under wrap at N=1,2,3 and 33/65/97 voxels; C3".into(),
        "       FALSIFIED - the crate has no periodic-wrap extraction, so the oracle".into(),
        "       cannot join the validity suite without a trait-signature change".into(),
        "P-143  -8/-4/-16 all reproduced (9 of 10 in-scope rows exact, the 10th the".into(),
        "       Schwarz P pinch) and C2 held, but C1 FALSIFIED: neither new surface".into(),
        "       was added as a reference field - that is a Phase 28 ticket".into(),
        "P-144  periodic value noise, wrapped: at 1 octave all 35 rows read".into(),
        "       boundary_edges 0 and oracle_exists true - every seed converges,".into(),
        "       which the registration did not expect. At 3 octaves no oracle on".into(),
        "       any row. Converged values -8/-8/-8/-12/-16, so a noise gate is".into(),
        "       per-seed (C2 held, C1 false).".into(),
    ]
}

fn draw_defects(
    measurement: Res<Measurement>,
    mut boundary: Gizmos<BoundaryGizmos>,
    mut pinches: Gizmos<PinchGizmos>,
) {
    let Some(reading) = &measurement.reading else {
        return;
    };
    const RED: Color = Color::srgb(0.98, 0.28, 0.24);
    const YELLOW: Color = Color::srgb(1.0, 0.95, 0.20);
    for [a, b] in &reading.boundary {
        boundary.line(*a, *b, RED);
    }
    for pinch in &reading.pinches {
        pinches
            .sphere(
                Isometry3d::from_translation(*pinch),
                reading.pinch_radius,
                YELLOW,
            )
            .resolution(8);
    }
}

#[cfg(test)]
mod tests {
    use bevy::ecs::system::RunSystemOnce;

    use super::*;

    /// The demo's own rebuild in an `App` with no window and no renderer.
    ///
    /// This is the closest thing to running the demo that a machine with no
    /// display can do. `controls` is left out on purpose — it wants a keyboard
    /// and the capture rig — and the state it would have produced is inserted
    /// directly, so each test names one configuration. `report` is left out for
    /// the same reason `game_dig` leaves it out: it wants a [`DemoStats`], and
    /// the tests run it as a one-shot with one inserted.
    ///
    /// No `TimeUpdateStrategy`: nothing here reads `Time`, so there is no delta
    /// to pin.
    fn harness(surface: Tpms, periods: u32, voxels: u32, wrap: Wrap) -> App {
        let surface = Tpms::ALL
            .iter()
            .position(|k| *k == surface)
            .expect("every Tpms is in ALL");
        let periods = PERIODS
            .iter()
            .position(|p| *p == periods)
            .expect("the test asked for a rung of PERIODS");
        let voxels = VOXELS
            .iter()
            .position(|v| *v == voxels)
            .expect("the test asked for a rung of VOXELS");

        let mut app = App::new();
        app.add_plugins(MinimalPlugins)
            .add_plugins(bevy::asset::AssetPlugin::default())
            .init_asset::<Mesh>()
            .insert_resource(Demo {
                surface,
                periods,
                voxels,
                wrap,
            })
            .init_resource::<Measurement>()
            .init_resource::<DemoStats>()
            .add_systems(Update, rebuild);
        app.world_mut().spawn((
            Mesh3d::default(),
            DemoDomain {
                min: Vec3::ZERO,
                max: Vec3::ZERO,
            },
        ));
        app.world_mut().spawn(OrbitCamera::default());
        app
    }

    /// One frame, then the HUD system as a one-shot, then the panel.
    fn panel(app: &mut App) -> Vec<String> {
        app.update();
        app.world_mut()
            .run_system_once(report)
            .expect("the HUD system");
        let lines = app.world().resource::<DemoStats>().extra.clone();
        for line in &lines {
            println!("{line}");
        }
        lines
    }

    /// The line containing `needle`, or a failure naming the whole panel.
    fn find<'a>(lines: &'a [String], needle: &str) -> &'a str {
        lines
            .iter()
            .find(|line| line.contains(needle))
            .unwrap_or_else(|| panic!("no HUD line contains {needle:?}: {lines:#?}"))
    }

    /// The panel row *beginning* with `label`.
    ///
    /// Not a substring search: the panel says `chi` twice on purpose -- once
    /// deriving the prediction from the shift census and once comparing it with
    /// the measurement -- and a test that grabbed whichever came first would
    /// assert on the arithmetic rather than on the reading.
    fn row<'a>(lines: &'a [String], label: &str) -> &'a str {
        lines
            .iter()
            .find(|line| line.starts_with(label))
            .unwrap_or_else(|| panic!("no HUD row starts with {label:?}: {lines:#?}"))
    }

    /// The rightmost column of a two-column panel row: the measured value.
    fn measured(line: &str) -> &str {
        line.split_whitespace()
            .next_back()
            .expect("a panel row is never empty")
    }

    /// The conforming grid reproduces P-142's committed gyroid row, digit for
    /// digit, and says so on the two lines a reader is meant to read.
    ///
    /// **This is the only way to see this demo's screen on a machine with no
    /// display.** `V 5284  E 15876  F 10584` is `docs/experiments/p-142.csv`'s
    /// `gyroid, periods 1, resolution 34, periodic, marching_cubes` row —
    /// `vertices`, `edges`, `faces` — so this asserts that the demo's public-API
    /// pipeline lands on the committed measurement rather than merely on a
    /// plausible one.
    #[test]
    fn the_conforming_grid_reproduces_the_committed_gyroid_row() {
        let mut app = harness(Tpms::Gyroid, 1, 33, Wrap::Periodic);
        let lines = panel(&mut app);

        let chi = row(&lines, "chi ");
        assert!(
            chi.contains("-8") && chi.contains("reproduced exactly"),
            "the gyroid did not reproduce -8N^3 under wrap: {chi}"
        );
        let census = row(&lines, "surface ");
        assert!(
            census.contains("lattice bcc") && census.contains("2 primitive cell(s)"),
            "the measured shift census stopped deriving the gyroid's bcc lattice: {census}"
        );
        let counts = row(&lines, "V ");
        for column in ["V 5284", "E 15876", "F 10584", "genus 5"] {
            assert!(
                counts.contains(column),
                "the reading drifted off p-142.csv's gyroid row ({column}): {counts}"
            );
        }
        assert_eq!(
            measured(row(&lines, "boundary edges")),
            "0",
            "the wrap left the gyroid open"
        );
        assert_eq!(
            measured(row(&lines, "non-manifold edges")),
            "0",
            "the conforming grid produced a pinch"
        );
        assert!(
            row(&lines, "chi_measured").contains("no pinch and no gap"),
            "the identity line lost its clean case"
        );
    }

    /// Schwarz D on the degenerate grid misses, and the pinches name the gap
    /// exactly.
    ///
    /// This is p-143.csv's `degenerate_grid_control` row — `voxels_per_period
    /// 32`, `chi_measured -12` against `chi_predicted -16`,
    /// `validate_non_manifold_edges 4`. It is the reason the ladder offers an
    /// even rung at all, and the identity is the demo's whole diagnostic: a
    /// `chi` that is wrong by exactly the pinch count is a statement about the
    /// grid, not about the `-16` prediction.
    #[test]
    fn the_degenerate_grid_misses_by_exactly_its_pinches() {
        let mut app = harness(Tpms::SchwarzD, 1, 32, Wrap::Periodic);
        let lines = panel(&mut app);

        let chi = row(&lines, "chi ");
        assert!(
            chi.contains("-16") && chi.contains("-12") && chi.contains("+4"),
            "the degenerate grid stopped missing by 4: {chi}"
        );
        let identity = row(&lines, "chi_measured");
        assert!(
            identity.contains("= 4") && identity.contains("EQUAL"),
            "the gap is no longer accounted for by the pinch count: {identity}"
        );
        assert_eq!(
            measured(row(&lines, "non-manifold edges")),
            "4",
            "the degenerate grid stopped producing four pinches"
        );
        assert_eq!(
            measured(row(&lines, "boundary edges")),
            "0",
            "the wrap left the degenerate arm open, so the gap is not a pinch story"
        );
        let counts = row(&lines, "V ");
        for column in ["V 4236", "E 12752", "F 8504"] {
            assert!(
                counts.contains(column),
                "the reading drifted off p-143.csv's degenerate_grid_control row \
                 ({column}): {counts}"
            );
        }
        let grid = row(&lines, "grid ");
        assert!(
            grid.contains("EVEN - degenerate control"),
            "the panel stopped warning that 32 is the even rung: {grid}"
        );
    }

    /// The open arm is refused even when its `chi` agrees.
    ///
    /// Non-wrapped Schwarz P at `N = 1` reads `-4` — its own prediction, by
    /// coincidence of the caps the box cuts (P-142). A control that compared
    /// only `chi` would pass the wrong arm on one field in three, so the verdict
    /// has to be gated on `boundary_edges`. Asserting that `chi` *does* agree
    /// here is what keeps the test from passing for the wrong reason: without
    /// it, a verdict that simply always said `NOT CLOSED` would look correct.
    #[test]
    fn the_open_arm_is_refused_even_when_chi_agrees() {
        let mut app = harness(Tpms::SchwarzP, 1, 33, Wrap::Open);
        let lines = panel(&mut app);

        let chi = row(&lines, "chi ");
        assert!(
            chi.contains("judged by them, never by chi"),
            "the open arm was not refused: {chi}"
        );
        let reading = app
            .world()
            .resource::<Measurement>()
            .reading
            .as_ref()
            .map(|r| {
                (
                    r.chi_predicted,
                    r.report.euler_characteristic,
                    r.report.boundary_edges,
                )
            })
            .expect("a reading");
        assert_eq!(
            (reading.0, reading.1),
            (-4, -4),
            "non-wrapped schwarz_p no longer reproduces its own prediction by coincidence, \
             so this test no longer exercises the trap P-142 named"
        );
        assert!(
            reading.2 > 0,
            "the open arm reported no boundary edges, so it is not open"
        );
        let census = row(&lines, "surface ");
        assert!(
            census.contains("lattice simple_cubic") && census.contains("1 primitive cell(s)"),
            "Schwarz P was labelled body-centred, which P-143 C2 measured as wrong: {census}"
        );
        assert!(
            find(&lines, "(pi,pi,pi)").contains("NEGATED"),
            "the body-centring shift stopped being measured as a negation"
        );
    }
}
