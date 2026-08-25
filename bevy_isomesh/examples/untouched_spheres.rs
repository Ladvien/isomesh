//! E-309 — the spheres nobody touches.
//!
//! ```bash
//! cd bevy_isomesh && cargo run --example untouched_spheres --release
//! ```
//!
//! **Always `--release`.** The startup self-check is four exhaustive scans —
//! 274,625 samples against every mesh vertex, no cutoff — and a debug build
//! turns seconds into minutes.
//!
//! `1` `thin_plate`, `2` `box_exact`, `M` flips the extractor, `[`/`]` change the
//! drawing stride, `,`/`.` move the slice, `T` hides the touched spheres, `H`
//! hides the surface. The rest are the shared keys — `W` wireframe, `G` domain
//! box, `Space` freeze, `F12` screenshot.
//!
//! Under `ISOMESH_CAPTURE` it needs no keyboard: the extractor alternates every
//! ten captured frames and the field steps halfway through, so
//! `record_gif.sh`'s default 80 frames is four MC↔DC flips on `thin_plate` and
//! four on `box_exact`. `ISOMESH_FIELD=0|1` pins one field and leaves the flips,
//! `ISOMESH_STRIDE=4|8|16|32` and `ISOMESH_SLICE=0..64` pick the subsample, and
//! `ISOMESH_SPIN=0.003` adds yaw. Every one of those exists so a committed still
//! can be regenerated from a command line rather than by holding a key down.
//!
//! ```bash
//! ISOMESH_WINDOW=1280x720 ./scripts/record_gif.sh untouched_spheres docs/gifs/e309.gif
//! ```
//!
//! Demonstrates **M-355 / P-51 clause C3**.
//!
//! # The constraint has two halves and only one of them is ever checked
//!
//! Sellan, Batty & Stein (`10.1145/3610548.3618196`) state what a signed
//! distance sample `p` with value `d(p)` actually asserts, and it is not a
//! number to interpolate on an edge. It is a **sphere**: the surface must
//! *exclude* the open ball of radius `|d(p)|` about `p`, **and it must touch
//! that ball's boundary at least once**.
//!
//! P-51 measured both halves as integer counts over the extractors' own output.
//! The exclusion half came back already satisfied — Dual Contouring pierces 0 of
//! 22,798 vertices, which is what falsified the pre-registered C2. The tangency
//! half is the one nobody had measured, and it is wide open: over the fifteen
//! ledger rows, between **3.1%** and **80.5%** of all samples have a sphere that
//! no vertex of the mesh ever reaches.
//!
//! That range is read out of `docs/experiments/p-51.csv` at startup rather than
//! quoted. **FINDINGS M-355's own title says "between 2.9% and 80.5%", and the
//! 2.9% is not in the artefact** — the smallest `untouched_per_1k` on any of the
//! fifteen rows is `31.2863`, i.e. 3.1%, on `sphere`/`dual_contouring`. The
//! ceiling reproduces exactly. The self-check logs the range it computed and
//! says so when the floor disagrees; nothing here is adjusted to match the
//! prose.
//!
//! # What is on screen
//!
//! - **Translucent steel surface** — the extracted mesh. Context, not subject;
//!   `H` removes it.
//! - **Orange wireframe spheres** — samples whose tangency sphere **no vertex
//!   ever touches**, at exactly the radius `|d(p)|` the SDF asserts. These are
//!   the finding, so they are the loud colour.
//! - **Faint cyan wireframe spheres** — samples the mesh *does* touch, to within
//!   the registered `0.05` cells.
//! - **Yellow crosshair, with a green segment to a green dot** — the single
//!   worst miss on this row: the sample whose sphere is missed by the widest
//!   margin, and the vertex that comes closest to its shell. The crosshair is a
//!   **locator** at a fixed eight cells, because on every one of these four rows
//!   the worst witness has `|d| = 0` and its sphere is a point. The green segment
//!   is the measured quantity — its length minus the sphere's radius *is* the
//!   reported number, one subtraction apart.
//!
//! # Switching the extractor is the measurement
//!
//! On `thin_plate` Marching Cubes leaves **717.1816 per 1k** untouched against
//! Dual Contouring's **44.8284** — **15.9984×**, which is the `16.00×` in
//! M-355's table. `M` flips between them, and at the default stride the drawn
//! set goes from **8 orange of 9** to **0 of 9**: every shell in the slice turns
//! from missed to touched. That flip is what a 16× ratio looks like in one
//! frame.
//!
//! `box_exact` is here because both extremes there are **exact geometric
//! constants**, and this example measures the deltas rather than quoting them:
//!
//! - **Marching Cubes misses by `√2` cells, bit for bit.** The box corner
//!   sample `(−1, −1, −1)` has `|d| = 0`, so its "sphere" is a point lying *on*
//!   the surface — and Marching Cubes puts no vertex on a box edge or corner at
//!   all, because both endpoints of every spanning grid edge there read `f = 0`
//!   and there is no sign change to interpolate. The nearest vertex it does
//!   place is `√2 h` away, and the live `f64` equals `√2` with `delta` exactly
//!   `0`.
//! - **Dual Contouring misses by `1/√2` cells — to nine figures, and the tenth
//!   has a name.** Its vertices live on the dual lattice, half a cell off the
//!   sample lattice, so a sample sitting exactly on the surface can never be
//!   touched however good the vertex placement is. The measured worst miss is
//!   `0.707106782954315` against `1/√2 = 0.707106781186547` — `1.767767e-9`
//!   cells over. M-355 calls it "1/√2 exactly", and the six decimals `p-51.csv`
//!   prints cannot separate the two, so this example measures the difference and
//!   then **predicts it**: `Clamp::ToCell` insets the QEF vertex by
//!   `o = CLAMP_EPSILON × h/2`, that inset is *perpendicular* to the `h/√2`
//!   displacement the miss is made of, and a perpendicular offset enters a
//!   distance only at second order — `o²/2b`, which in cells is `ε²√2/8`
//!   independent of `h`. That is `1.767767e-9`, and the live delta is that to
//!   within `3e-7` relative.
//!
//! `thin_plate`/MC's `1.000000` comes back bit-exact too, so the split runs
//! between the two *extractors* rather than between the two fields: a Marching
//! Cubes vertex is an interpolation between grid corners and at these witnesses
//! the interpolant is a grid point, while a Dual Contouring vertex is a clamped
//! least-squares solve and carries the clamp's fingerprint into the ninth digit.
//!
//! # The counts are the exhaustive pass, and the picture is a stated subsample
//!
//! P-51 applies **no cutoff** to the touching search: the minimum is taken over
//! every vertex of the mesh for every one of the 274,625 samples, and `p-51.csv`
//! records `touch_search = exhaustive` on all fifteen rows. This example does
//! the same scan, in full, once per `(field, extractor)` pair before the window
//! opens — four scans, on four threads, and the numbers on the HUD are that
//! scan's, never a re-derivation.
//!
//! **Only the drawing is subsampled.** 274,625 wireframe spheres is fog, so one
//! slice of the lattice is drawn at a stride, and both the stride and the drawn
//! count are on the HUD next to the census count so the two cannot be confused.
//! A sphere's colour is the verdict the *census* recorded for that sample, read
//! out of the same `bool` the count was made of — the picture cannot disagree
//! with the number.
//!
//! # 65 samples per axis, and `ISOMESH_SAMPLES` is refused
//!
//! Every number here exists only at the registered 65³. There is no ledger row
//! at any other resolution to hold the live run against, and the scan is
//! `O(samples³ × vertices)` — at 129 it is two orders of magnitude more work.
//! So `ISOMESH_SAMPLES` is **rejected with an error** rather than honoured into
//! a HUD full of unfalsifiable numbers.
//!
//! # The self-check runs on a machine with no display
//!
//! It is `println!` before `App::new()`, not `info!` inside a system, and that
//! is measured rather than stylistic: `add_plugins(DefaultPlugins)` builds both
//! the log subscriber *and* `WinitPlugin`'s event loop, and the latter panics
//! outright where there is no X display. Anything logged after it would be
//! unreachable exactly where a terminal is the only output. So
//! `cargo run --example untouched_spheres --release` prints the whole audit —
//! every column, every witness, every delta — before it tries to open a window.
//!
//! # `f64`
//!
//! P-51 was measured in `f64`, so the counts reproduce only in `f64`; the
//! surface is cast to `f32` on its way into the [`Mesh`] asset and nothing but
//! the picture depends on that.

mod common;

use std::f64::consts::SQRT_2;
use std::time::Instant;

use bevy::prelude::*;
use bevy_isomesh::MeshBuilder;
use common::{Capture, CommonPlugin, DemoDomain, DemoMesh, DemoStats, OrbitCamera, ViewFlags};
use isomesh::dual_contouring::{CLAMP_EPSILON, DualContouring};
use isomesh::fields::{BoxExact, ReferenceField, ThinPlate};
use isomesh::marching_cubes::MarchingCubes;
use isomesh::{MeshBuffer, MeshSink, RuntimeShape3, Sdf};

// ─── the registered measurement ─────────────────────────────────────────────

/// Samples per axis. Registered in P-51, and not a knob — see the module docs.
const SAMPLES_PER_AXIS: u32 = 65;

/// The gate, in cells. Registered: a sample counts as untouched when the nearest
/// vertex misses its shell by more than this.
const THRESHOLD_CELLS: f64 = 0.05;

/// The ledger, compiled in.
///
/// `include_str!` rather than a runtime path: the artefact this example claims
/// to reproduce is then *in the binary*, so the self-check cannot silently pass
/// because a file moved, and a stale example cannot be linked against a fresh
/// CSV.
const LEDGER: &str = include_str!("../../docs/experiments/p-51.csv");

/// The two fields, in the order `1` and `2` select them.
///
/// `thin_plate` leads: it is the field with the 16× MC/DC spread, which is the
/// thing switching extractor is supposed to show.
const FIELDS: [&str; 2] = ["thin_plate", "box_exact"];

/// The two extractors, in the order `M` alternates them.
const EXTRACTORS: [&str; 2] = ["marching_cubes", "dual_contouring"];

/// `(field, extractor)` pairs, all measured at startup.
const COMBOS: usize = FIELDS.len() * EXTRACTORS.len();

/// Where `(field, extractor)` lives in the flat combo list.
fn combo_index(field: usize, extractor: usize) -> usize {
    field * EXTRACTORS.len() + extractor
}

// ─── the lattice ────────────────────────────────────────────────────────────

/// `|d(p)|` at every sample: the radius of the sphere each sample asserts.
///
/// Stored in the `k`-major order P-51 scans in, so the reduction order — and
/// therefore which sample wins a tie for "worst" — is the ledger's.
struct Lattice {
    abs_d: Vec<f64>,
    lo: [f64; 3],
    h: f64,
    n: usize,
}

impl Lattice {
    /// Sample `field` at [`SAMPLES_PER_AXIS`] per axis over its own domain.
    ///
    /// The cell size is `(hi.x − lo.x) / (samples − 1)`, which is
    /// `benches/common::grid`'s definition verbatim: `n` samples span `n − 1`
    /// cells, and a demo that spanned `n` would be measuring a different grid
    /// from the row it quotes.
    fn of<F>(field: &F) -> Self
    where
        F: ReferenceField + Sdf<Scalar = f64>,
    {
        let (lo, hi) = field.domain();
        let h = (hi[0] - lo[0]) / f64::from(SAMPLES_PER_AXIS - 1);
        let n = SAMPLES_PER_AXIS as usize;
        let mut abs_d = vec![0.0_f64; n * n * n];
        for k in 0..n {
            let z = lo[2] + k as f64 * h;
            for j in 0..n {
                let y = lo[1] + j as f64 * h;
                for i in 0..n {
                    let x = lo[0] + i as f64 * h;
                    abs_d[(k * n + j) * n + i] = field.sample([x, y, z]).abs();
                }
            }
        }
        Self { abs_d, lo, h, n }
    }

    /// Position of the sample at index `(i, j, k)`.
    fn at(&self, i: usize, j: usize, k: usize) -> [f64; 3] {
        [
            self.lo[0] + i as f64 * self.h,
            self.lo[1] + j as f64 * self.h,
            self.lo[2] + k as f64 * self.h,
        ]
    }

    /// The flat index the scan and the drawing both address samples by.
    fn index(&self, i: usize, j: usize, k: usize) -> usize {
        (k * self.n + j) * self.n + i
    }

    fn count(&self) -> usize {
        self.abs_d.len()
    }
}

// ─── the touching half ──────────────────────────────────────────────────────

/// What one exhaustive touching scan found.
struct Scan {
    untouched: u64,
    worst_cells: f64,
    /// Per sample: is its sphere untouched? **The same `bool` the count is made
    /// of**, kept so the picture colours a sphere the way the census counted it
    /// rather than re-deciding from a rounded copy.
    untouched_mask: Vec<bool>,
    /// The worst miss, with the geometry behind it — the sample, its sphere's
    /// radius, the vertex nearest that shell, and how far that vertex actually
    /// is. `distance − radius` is the reported miss, one subtraction away.
    sample: [f64; 3],
    radius: f64,
    nearest: [f64; 3],
    distance: f64,
}

/// `touch(p) = min over mesh vertices v of |‖v − p‖ − |d(p)||`, in cells, for
/// every sample, over **every** vertex.
///
/// Transcribed from `benches/experiment_p51.rs::Samples::scan`, touching half
/// only, including the pruning: a vertex can beat the running best `t` only if
/// `‖v − p‖ ∈ (r − t, r + t)`, so a squared distance outside `((r−t)², (r+t)²)`
/// is rejected with three subtractions, three multiplies, two adds and two
/// compares, and the square root is paid only on a candidate that actually
/// narrows the window. That is what makes an exhaustive search affordable, and
/// it is why there is no cutoff to blame the count on.
///
/// The `>` in the worst-case update is strict, as in the bench, so the first
/// sample in `k`-major order wins a tie and the witness on the HUD is the
/// witness in the CSV's log.
fn scan(lattice: &Lattice, verts: &[[f64; 3]]) -> Scan {
    let inv_h = lattice.h.recip();
    let n = lattice.n;
    let mut out = Scan {
        untouched: 0,
        worst_cells: 0.0,
        untouched_mask: vec![false; lattice.count()],
        sample: [0.0; 3],
        radius: 0.0,
        nearest: [0.0; 3],
        distance: 0.0,
    };
    for k in 0..n {
        for j in 0..n {
            for i in 0..n {
                let index = lattice.index(i, j, k);
                let p = lattice.at(i, j, k);
                let r = lattice.abs_d[index];
                let mut best = f64::INFINITY;
                let mut best_at = [0.0_f64; 3];
                let mut best_distance = 0.0;
                let mut window_lo = f64::NEG_INFINITY;
                let mut window_hi = f64::INFINITY;
                for v in verts {
                    let dx = v[0] - p[0];
                    let dy = v[1] - p[1];
                    let dz = v[2] - p[2];
                    let d2 = dx * dx + dy * dy + dz * dz;
                    if d2 > window_lo && d2 < window_hi {
                        let d = d2.sqrt();
                        let t = (d - r).abs();
                        if t < best {
                            best = t;
                            best_at = *v;
                            best_distance = d;
                            let inner = r - t;
                            window_lo = if inner > 0.0 {
                                inner * inner
                            } else {
                                f64::NEG_INFINITY
                            };
                            window_hi = (r + t) * (r + t);
                        }
                    }
                }
                let cells = best * inv_h;
                if cells > THRESHOLD_CELLS {
                    out.untouched += 1;
                    out.untouched_mask[index] = true;
                }
                if cells > out.worst_cells {
                    out.worst_cells = cells;
                    out.sample = p;
                    out.radius = r;
                    out.nearest = best_at;
                    out.distance = best_distance;
                }
            }
        }
    }
    out
}

// ─── one row ────────────────────────────────────────────────────────────────

/// One `(field, extractor)` row: the extraction, and the scan over it.
struct Combo {
    field: usize,
    extractor: usize,
    samples: usize,
    /// The extraction the scan minimised over. Kept rather than discarded: the
    /// mesh on screen has to be the one the numbers were computed from, and the
    /// nearest-vertex marker is one of its positions.
    buffer: MeshBuffer<f64>,
    scan: Scan,
    extract_ms: f64,
    scan_ms: f64,
}

impl Combo {
    /// Untouched samples per 1,000 samples, exactly as `Row::untouched_per_1k`
    /// computes it in the bench.
    fn per_1k(&self) -> f64 {
        if self.samples == 0 {
            return 0.0;
        }
        1000.0 * self.scan.untouched as f64 / self.samples as f64
    }
}

/// Extract with one extractor, then scan the whole lattice against it.
fn extract_and_scan<F>(field: usize, extractor: usize, sdf: &F, lattice: &Lattice) -> Combo
where
    F: ReferenceField + Sdf<Scalar = f64>,
{
    let shape = RuntimeShape3::new([SAMPLES_PER_AXIS; 3]).expect("65^3 fits u32");
    let mut buffer = MeshBuffer::<f64>::new();
    let started = Instant::now();
    // Matched rather than dispatched through `Extractor`: that trait's method is
    // generic over the field and the sink, so it is not object-safe, and both
    // inherent `extract`s have this exact signature anyway.
    let outcome = match extractor {
        0 => MarchingCubes::<f64>::new().extract(sdf, &shape, lattice.lo, lattice.h, &mut buffer),
        _ => DualContouring::<f64>::new().extract(sdf, &shape, lattice.lo, lattice.h, &mut buffer),
    };
    if let Err(error) = outcome {
        // Loud, and before the window exists. There is no degraded picture to
        // fall back to: with no mesh there is nothing for a sphere to touch and
        // every number on the HUD would be a count of nothing.
        panic!(
            "E-309: {} could not extract {} at {SAMPLES_PER_AXIS}^3: {error}",
            EXTRACTORS[extractor],
            F::NAME
        );
    }
    let extract_ms = started.elapsed().as_secs_f64() * 1000.0;

    let started = Instant::now();
    let scan = scan(lattice, &buffer.positions);
    let scan_ms = started.elapsed().as_secs_f64() * 1000.0;

    Combo {
        field,
        extractor,
        samples: lattice.count(),
        buffer,
        scan,
        extract_ms,
        scan_ms,
    }
}

/// The lattice of one field, chosen by index.
///
/// The reference fields are distinct types and both the sampler and the
/// extractors are generic over a sized field, so a runtime choice is a `match`
/// rather than a `&dyn Sdf` — the same shape `critical_cells` uses.
fn lattice_of(field: usize) -> Lattice {
    match field {
        0 => Lattice::of(&ThinPlate::<f64>::canonical()),
        _ => Lattice::of(&BoxExact::<f64>::canonical()),
    }
}

/// One row, chosen by index.
fn measure(field: usize, extractor: usize, lattice: &Lattice) -> Combo {
    match field {
        0 => extract_and_scan(field, extractor, &ThinPlate::<f64>::canonical(), lattice),
        _ => extract_and_scan(field, extractor, &BoxExact::<f64>::canonical(), lattice),
    }
}

// ─── the ledger ─────────────────────────────────────────────────────────────

/// `docs/experiments/p-51.csv`, split into a header and its rows.
///
/// Addressed **by column name**, never by position. A CSV whose columns were
/// reordered would still be the same artefact, and an example that read column
/// 9 by number would silently start comparing the wrong quantity.
struct Ledger {
    header: Vec<&'static str>,
    rows: Vec<Vec<&'static str>>,
}

impl Ledger {
    fn parse() -> Self {
        let mut lines = LEDGER
            .lines()
            .filter(|line| !line.starts_with('#') && !line.trim().is_empty());
        let header = lines
            .next()
            .map(|line| line.split(',').collect())
            .unwrap_or_default();
        let rows = lines.map(|line| line.split(',').collect()).collect();
        Self { header, rows }
    }

    fn column(&self, name: &str) -> Option<usize> {
        self.header.iter().position(|head| *head == name)
    }

    /// One cell, by field, extractor and column name.
    fn cell(&self, field: &str, extractor: &str, column: &str) -> Option<&'static str> {
        let f = self.column("field")?;
        let e = self.column("extractor")?;
        let c = self.column(column)?;
        self.rows
            .iter()
            .find(|row| {
                row.get(f).copied() == Some(field) && row.get(e).copied() == Some(extractor)
            })?
            .get(c)
            .copied()
    }

    /// The smallest and largest `untouched_per_1k` over **every** row, as
    /// percentages — the range M-355's title states.
    fn untouched_range_pct(&self) -> Option<(f64, f64)> {
        let c = self.column("untouched_per_1k")?;
        let mut lo = f64::INFINITY;
        let mut hi = f64::NEG_INFINITY;
        for row in &self.rows {
            let value: f64 = row.get(c)?.parse().ok()?;
            lo = lo.min(value);
            hi = hi.max(value);
        }
        if lo.is_finite() && hi.is_finite() {
            Some((lo / 10.0, hi / 10.0))
        } else {
            None
        }
    }
}

// ─── the self-check ─────────────────────────────────────────────────────────

/// The floor and ceiling M-355's title states, in percent of samples.
const FINDING_RANGE_PCT: (f64, f64) = (2.9, 80.5);

/// Runs the ledger comparisons and counts them.
struct Auditor<'a> {
    ledger: &'a Ledger,
    held: usize,
    total: usize,
}

impl Auditor<'_> {
    /// One live value against one committed cell, string against string.
    ///
    /// String equality rather than a float tolerance: the CSV *is* text, the
    /// bench wrote it with a stated format, and formatting the live value the
    /// same way makes "matches to the digit" the literal thing being tested.
    fn hold(&mut self, field: &str, extractor: &str, column: &str, live: &str) {
        self.total += 1;
        match self.ledger.cell(field, extractor, column) {
            Some(recorded) if recorded == live => {
                self.held += 1;
                println!("  held   {field:<10} {extractor:<15} {column} = {live}");
            }
            Some(recorded) => eprintln!(
                "  BROKE  {field:<10} {extractor:<15} {column}: live {live}, p-51.csv {recorded}"
            ),
            None => {
                eprintln!("  BROKE  {field:<10} {extractor:<15} {column} absent from p-51.csv");
            }
        }
    }
}

/// One closed form M-355 states for a row's worst miss, and where the live
/// `f64` actually lands relative to it.
struct Identity {
    /// The combo it belongs to, as [`combo_index`] numbers them.
    combo: usize,
    /// How M-355 writes the constant.
    label: &'static str,
    delta: f64,
    /// The residual predicted from the extractor's own published constants,
    /// when the delta has a closed form of its own. See [`clamp_residual`].
    predicted: Option<f64>,
}

/// The residual `box_exact`/`dual_contouring` carries on top of `1/√2`,
/// predicted from Dual Contouring's own [`CLAMP_EPSILON`].
///
/// `Clamp::ToCell` insets a QEF vertex to `half × (1 − ε)` about its cell
/// centre, so a vertex that wants to sit exactly on a cell face lands
/// `o = ε h/2` inside it. At this witness that inset is **perpendicular** to the
/// `b = h/√2` displacement the miss is made of, and a perpendicular offset
/// enters a distance only at second order:
///
/// ```text
/// sqrt(b² + o²) − b  ≈  o²/2b  =  (ε h/2)² / (2 h/√2)  =  ε² √2 h / 8
/// ```
///
/// In cells the `h` cancels, so the prediction is `ε²√2/8` — the same residual
/// at every resolution, and `1.767767e−9` at `ε = 1e−4`. Measured: the live
/// delta is that to within `3e−7` relative, which is the QEF solve's own noise
/// riding on top of a clamp whose fingerprint is exact.
fn clamp_residual() -> f64 {
    CLAMP_EPSILON * CLAMP_EPSILON * SQRT_2 / 8.0
}

/// Measure one closed form, and say how close it came.
///
/// **Reported, not gated, and deliberately kept out of the ledger audit** —
/// they are two different claims and conflating them would have hidden the more
/// interesting one. `worst_untouched_cells` matching `p-51.csv` to the digit is
/// this example reproducing the artefact, and it does, on all four rows.
/// Whether that six-decimal figure is *literally* the closed form in `f64` is a
/// statement about the extractor, and the two answer differently:
///
/// - `thin_plate`/MC hits `1` and `box_exact`/MC hits `√2` **bit for bit**,
///   `delta` exactly `0`. A Marching Cubes vertex is a linear interpolation
///   between two grid corners, and at these witnesses the interpolant *is* a
///   grid point, so the arithmetic never leaves the lattice.
/// - `box_exact`/DC lands `1.768e-9` cells off `1/√2` — and that residual is not
///   noise. It is `CLAMP_EPSILON` seen through a Pythagorean square root, to
///   within `3e-7` relative of [`clamp_residual`]'s closed form.
///
/// So M-355's "1/√2 exactly" is exact to nine significant figures rather than to
/// the last bit, and the missing part is *nameable*. Six decimals in a CSV
/// cannot tell those apart; this prints the delta and its prediction so a reader
/// does not have to guess.
fn identity(
    combo: usize,
    label: &'static str,
    live: f64,
    exact: f64,
    predicted: Option<f64>,
) -> Identity {
    let delta = (live - exact).abs();
    println!(
        "  closed form  {:<10} {:<15} worst miss vs {label}: live {live:.15}, \
         exact {exact:.15}, delta {delta:.6e} cells",
        FIELDS[combo / EXTRACTORS.len()],
        EXTRACTORS[combo % EXTRACTORS.len()],
    );
    if let Some(prediction) = predicted {
        println!(
            "               that delta predicted from CLAMP_EPSILON = {CLAMP_EPSILON:e} as \
             eps^2*sqrt(2)/8 = {prediction:.6e} cells; measured/predicted = {:.7}",
            delta / prediction
        );
    }
    Identity {
        combo,
        label,
        delta,
        predicted,
    }
}

/// What the self-check concluded, for the HUD.
struct Check {
    /// Ledger columns that matched, and how many were compared.
    held: usize,
    total: usize,
    /// The closed forms M-355 names, with their measured deltas.
    identities: Vec<Identity>,
    /// The `untouched_per_1k` range over all fifteen ledger rows, in percent.
    range_pct: (f64, f64),
    /// Does that range agree with M-355's title to one decimal?
    range_agrees: bool,
}

impl Check {
    /// The closed form for one combo, if M-355 names one.
    fn identity(&self, combo: usize) -> Option<&Identity> {
        self.identities.iter().find(|i| i.combo == combo)
    }
}

/// Reproduce every headline number, then hold the whole row against the CSV.
///
/// **Nothing here panics on a mismatch.** This is a demo a stranger runs: a
/// broken claim must be loud on the terminal and on the HUD, and must not take
/// the window down with it — the picture of the disagreement is the useful
/// artefact.
///
/// # `println!` rather than `info!`, and the reason is measured
///
/// `tracing` events go nowhere until Bevy's `LogPlugin` installs a subscriber,
/// which happens inside `add_plugins(DefaultPlugins)` — and so does
/// `WinitPlugin`'s event loop, which **panics on a machine with no X display**:
/// `Failed to build event loop: XNotSupported(XOpenDisplayFailed)`, verified on
/// this repo's own CI-shaped box. A self-check written with `info!` is therefore
/// unreachable in exactly the situation where a terminal is all you have.
///
/// So this runs *before* the plugins and writes to stdout, held lines on stdout
/// and broken ones on stderr. `cargo run --example untouched_spheres --release`
/// prints the whole audit and reproduces every number in P-51's C3 rows even
/// where no window can open. E-203's lesson, applied: a measurement that only
/// exists on screen cannot be verified from a terminal.
fn self_check(combos: &[Combo]) -> Check {
    let ledger = Ledger::parse();
    let mut audit = Auditor {
        ledger: &ledger,
        held: 0,
        total: 0,
    };

    println!(
        "E-309 self-check — reproducing P-51 clause C3 against docs/experiments/p-51.csv, \
         {SAMPLES_PER_AXIS}^3, f64, exhaustive touching search"
    );

    for combo in combos {
        let field = FIELDS[combo.field];
        let extractor = EXTRACTORS[combo.extractor];
        audit.hold(
            field,
            extractor,
            "samples_per_axis",
            &SAMPLES_PER_AXIS.to_string(),
        );
        audit.hold(field, extractor, "samples", &combo.samples.to_string());
        audit.hold(
            field,
            extractor,
            "vertices",
            &combo.buffer.vertex_count().to_string(),
        );
        audit.hold(
            field,
            extractor,
            "triangles",
            &combo.buffer.triangle_count().to_string(),
        );
        audit.hold(
            field,
            extractor,
            "threshold_cells",
            &format!("{THRESHOLD_CELLS}"),
        );
        audit.hold(field, extractor, "touch_search", "exhaustive");
        audit.hold(
            field,
            extractor,
            "vertices_probed_per_sample",
            &combo.buffer.vertex_count().to_string(),
        );
        audit.hold(
            field,
            extractor,
            "untouched",
            &combo.scan.untouched.to_string(),
        );
        audit.hold(
            field,
            extractor,
            "untouched_per_1k",
            &format!("{:.4}", combo.per_1k()),
        );
        audit.hold(
            field,
            extractor,
            "worst_untouched_cells",
            &format!("{:.6}", combo.scan.worst_cells),
        );
    }

    // The ratio is a quantity *between* rows, so it is computed after every row
    // exists — the same order the bench writes it in.
    for field in 0..FIELDS.len() {
        let mc = &combos[combo_index(field, 0)];
        let dc = &combos[combo_index(field, 1)];
        let ratio = format!("{:.4}", mc.per_1k() / dc.per_1k());
        for extractor in EXTRACTORS {
            audit.hold(FIELDS[field], extractor, "untouched_mc_over_dc", &ratio);
        }
    }

    // The three worst cases M-355 names as closed forms. Reported beside the
    // ledger audit, never folded into it — see `identity`.
    let identities = vec![
        identity(
            combo_index(0, 0),
            "1 cell",
            combos[combo_index(0, 0)].scan.worst_cells,
            1.0,
            None,
        ),
        identity(
            combo_index(1, 0),
            "sqrt(2)",
            combos[combo_index(1, 0)].scan.worst_cells,
            SQRT_2,
            None,
        ),
        identity(
            combo_index(1, 1),
            "1/sqrt(2)",
            combos[combo_index(1, 1)].scan.worst_cells,
            SQRT_2.recip(),
            Some(clamp_residual()),
        ),
    ];

    // The four numbers this example exists to put on screen, each with the
    // geometry behind it, so a headless run leaves the witness on the terminal
    // and `distance - radius` can be checked by subtracting two printed floats
    // rather than by trusting a count. The bench prints its witnesses for the
    // same reason.
    for combo in combos {
        let scan = &combo.scan;
        println!(
            "  {:<10} {:<15} untouched {:>6} ({:>8.4}/1k, {:>5.2}% of samples), \
             worst miss {:.6} cells; extract {:.1} ms, exhaustive scan {:.0} ms over {} vertices",
            FIELDS[combo.field],
            EXTRACTORS[combo.extractor],
            scan.untouched,
            combo.per_1k(),
            combo.per_1k() / 10.0,
            scan.worst_cells,
            combo.extract_ms,
            combo.scan_ms,
            combo.buffer.vertex_count(),
        );
        println!(
            "      witness: sample [{:.6}, {:.6}, {:.6}] asserts |d| = {:.6}; its \
             closest-to-shell vertex [{:.6}, {:.6}, {:.6}] is {:.6} away — miss {:.6} = \
             {:.6} cells",
            scan.sample[0],
            scan.sample[1],
            scan.sample[2],
            scan.radius,
            scan.nearest[0],
            scan.nearest[1],
            scan.nearest[2],
            scan.distance,
            scan.distance - scan.radius,
            scan.worst_cells,
        );
    }

    // M-355's title states a range over all fifteen rows, not just these four,
    // so it is checked against the whole artefact.
    let range_pct = ledger.untouched_range_pct().unwrap_or((f64::NAN, f64::NAN));
    let agrees = (range_pct.0 - FINDING_RANGE_PCT.0).abs() < 0.05
        && (range_pct.1 - FINDING_RANGE_PCT.1).abs() < 0.05;
    if agrees {
        println!(
            "  held   p-51.csv untouched_per_1k spans {:.1}% to {:.1}% of samples, \
             which is M-355's stated range",
            range_pct.0, range_pct.1
        );
    } else {
        eprintln!(
            "  DISAGREES  p-51.csv untouched_per_1k spans {:.1}% to {:.1}% of samples over its \
             {} rows; FINDINGS M-355's title says {:.1}% to {:.1}%. The artefact is what this \
             example reports.",
            range_pct.0,
            range_pct.1,
            ledger.rows.len(),
            FINDING_RANGE_PCT.0,
            FINDING_RANGE_PCT.1,
        );
    }

    if audit.held == audit.total {
        println!(
            "E-309 self-check: {}/{} ledger columns held to the digit",
            audit.held, audit.total
        );
    } else {
        eprintln!(
            "E-309 self-check: only {}/{} ledger columns held — the live run and p-51.csv \
             disagree, and the HUD says so",
            audit.held, audit.total
        );
    }

    Check {
        held: audit.held,
        total: audit.total,
        identities,
        range_pct,
        range_agrees: agrees,
    }
}

// ─── everything measured, once ──────────────────────────────────────────────

/// The four rows, their lattices, and the verdict — computed before the window
/// opens and never recomputed.
///
/// Caching all four is not an optimisation on top of a live path; it is the only
/// path. The self-check needs every row anyway, so `M` and `1`/`2` are a lookup
/// and the picture always shows a number that has been held against the ledger.
#[derive(Resource)]
struct Measured {
    lattices: Vec<Lattice>,
    combos: Vec<Combo>,
    check: Check,
}

impl Measured {
    fn combo(&self, field: usize, extractor: usize) -> &Combo {
        &self.combos[combo_index(field, extractor)]
    }

    fn lattice(&self, field: usize) -> &Lattice {
        &self.lattices[field]
    }

    /// Marching Cubes' untouched rate over Dual Contouring's, on one field.
    fn mc_over_dc(&self, field: usize) -> f64 {
        self.combo(field, 0).per_1k() / self.combo(field, 1).per_1k()
    }
}

/// Sample both fields, run all four scans, and check the lot.
///
/// The four scans go on four threads. Each owns its own field instance, its own
/// mesh buffer and its own output, reads a shared immutable lattice, and is
/// joined in a fixed order — so this is exactly as deterministic as the
/// sequential version, and startup is one scan long rather than four.
fn measure_everything() -> Measured {
    if let Some(requested) = common::samples_override()
        && requested != SAMPLES_PER_AXIS
    {
        eprintln!(
            "ISOMESH_SAMPLES={requested} refused: every number in E-309 is a P-51 row at \
             {SAMPLES_PER_AXIS}^3 and there is no committed row at any other resolution to \
             hold a live run against. Continuing at {SAMPLES_PER_AXIS}^3."
        );
    }

    let started = Instant::now();
    let lattices: Vec<Lattice> = (0..FIELDS.len()).map(lattice_of).collect();
    let combos: Vec<Combo> = std::thread::scope(|scope| {
        let handles: Vec<_> = (0..COMBOS)
            .map(|index| {
                let lattice = &lattices[index / EXTRACTORS.len()];
                scope.spawn(move || {
                    measure(index / EXTRACTORS.len(), index % EXTRACTORS.len(), lattice)
                })
            })
            .collect();
        handles
            .into_iter()
            .map(|handle| handle.join().expect("E-309: a scan thread panicked"))
            .collect()
    });
    let check = self_check(&combos);
    let total_ms = started.elapsed().as_secs_f64() * 1000.0;
    println!("E-309: four exhaustive scans in {total_ms:.0} ms wall clock");

    Measured {
        lattices,
        combos,
        check,
    }
}

// ─── the picture ────────────────────────────────────────────────────────────

/// The strides `[` and `]` cycle, in samples.
///
/// Powers of two, and all of them divide 64 — so a stride always draws the
/// slice's two edge samples as well as its centre, and the drawn set is a
/// symmetric subsample rather than a lopsided one.
const STRIDES: [usize; 4] = [4, 8, 16, 32];

/// The default, as an index into [`STRIDES`].
///
/// `32`, giving 3×3 = 9 spheres in the slice, and the number was read off the
/// frames rather than chosen. **The overlap is intrinsic rather than a rendering
/// problem**: a sample `m` cells from the surface asserts a sphere of radius
/// about `m` cells, so at stride `s` every sample further than `s/2` from the
/// surface has a shell wider than the gap to its neighbour and the shells nest.
/// Measured at 1280x720 scaled to 900: 81 shells and 25 shells are both a bird's
/// nest that swallows the surface, 9 reads as a diagram — one cyan shell resting
/// tangent on the plate ringed by eight orange ones — and the MC→DC flip takes
/// it from 8 untouched to 0, which is legible in a single frame.
///
/// `[` and `]` densify it to 16, 8 and 4 for anyone who wants the population
/// rather than the picture, and the HUD always says which stride is showing and
/// how many of the drawn spheres are untouched.
const DEFAULT_STRIDE: usize = 3;

/// How far `,` and `.` move the slice, in samples.
const SLICE_STEP: u32 = 4;

/// The slice drawn by default, as a sample index along `y`.
///
/// `40`, which on the `[−2, 2]` domain is `y = 0.5` — eight cells above the
/// plate rather than through it. See [`resolve`] for why the obvious choice,
/// the middle plane at `y = 0`, is the one plane of this family that hides the
/// finding: it is the only sample plane inside the plate, and the Marching Cubes
/// rim vertices it creates lie in the slice and touch the outer samples exactly.
const DEFAULT_SLICE: u32 = 40;

/// Radius of the orbit, in world units.
///
/// Both fields span `[−2, 2]³`, so this is one number rather than a per-field
/// framing. At Bevy's default 45° vertical field of view the visible half-height
/// is `0.414 × radius`. The subject is not the 4-unit domain but the shells,
/// and the outermost of them are centred on the domain corner with a radius of
/// about 1.5 — so `6.6` is what keeps the ring inside the frame instead of
/// letting its left and right sides leave it.
const VIEW_RADIUS: f32 = 6.6;

/// Where the subject sits in frame, right and down from centre, as a fraction of
/// the orbit radius.
///
/// The HUD and its backdrop occupy 500x524 pixels of the upper left and the
/// lattice is the subject; centring the lattice photographs the argument behind
/// its own evidence. Applied in
/// the camera's own basis, so it holds while `ISOMESH_SPIN` yaws.
const SUBJECT_OFFSET: Vec2 = Vec2::new(0.16, 0.06);

/// Yaw and pitch the slice reads best from.
///
/// The slice is a plane of constant `y`, so the camera looks **down** on it —
/// far enough down that the annulus of untouched shells reads as a ring rather
/// than as a row, and not so far that the plate and the box stop reading as
/// solids. Measured on 1280x720 frames scaled to 900: at a pitch of 0.28 the
/// ring collapsed into a band and the outer shells left the frame sideways.
const VIEW_YAW: f32 = 0.86;
const VIEW_PITCH: f32 = 0.95;

/// Captured frames per `(field, extractor)` stage.
///
/// Eight stages at ten frames is 80, which is `record_gif.sh`'s default
/// `ISOMESH_CAPTURE_FRAMES` — so the default capture is four MC↔DC flips on
/// `thin_plate` and four on `box_exact`. Flipping four times rather than once
/// matters: the finding is a *difference* between two populations, and a viewer
/// reads a difference from a repeat, not from a single cut.
const CAPTURE_FRAMES_PER_STAGE: u32 = 10;

/// Flips per field, when the capture is not pinned to one.
const FLIPS_PER_FIELD: u32 = 4;

/// Seconds per extractor when nobody is capturing.
const STAGE_SECONDS: f32 = 2.6;

/// The untouched spheres: the finding, and the loud colour.
const UNTOUCHED: Color = Color::srgb(1.0, 0.36, 0.13);

/// The touched ones: present, so the ratio has a denominator on screen, and dim,
/// so it is not competing for the eye.
const TOUCHED: Color = Color::srgba(0.30, 0.80, 0.95, 0.22);

/// The worst miss, and the vertex that came closest to its shell.
const WORST: Color = Color::srgb(1.0, 0.93, 0.28);
const NEAREST: Color = Color::srgb(0.45, 1.0, 0.45);

/// Which row is on screen, and how long it has been.
#[derive(Resource, Default)]
struct Stage {
    field: usize,
    extractor: usize,
    /// Seconds on this extractor, when nobody is capturing.
    phase: f32,
}

/// How the lattice is subsampled for drawing, and what is drawn.
#[derive(Resource)]
struct Draw {
    stride: usize,
    slice: u32,
    touched: bool,
    surface: bool,
}

impl Default for Draw {
    /// `ISOMESH_STRIDE` and `ISOMESH_SLICE` pick the subsample without a
    /// keyboard.
    ///
    /// The harness's rule is that anything a capture depends on is reachable
    /// from the environment — that is what makes a committed still
    /// regenerable from a command line rather than by holding `]` down. A
    /// stride that is not one of [`STRIDES`] is refused rather than rounded:
    /// silently drawing a different subsample than the one asked for is exactly
    /// the failure `ISOMESH_SAMPLES` was added to fix.
    fn default() -> Self {
        let stride = match std::env::var("ISOMESH_STRIDE")
            .ok()
            .and_then(|v| v.parse::<usize>().ok())
        {
            Some(requested) => match STRIDES.iter().position(|s| *s == requested) {
                Some(index) => index,
                None => {
                    eprintln!(
                        "ISOMESH_STRIDE={requested} refused: the strides are {STRIDES:?}. \
                         Using {}.",
                        STRIDES[DEFAULT_STRIDE]
                    );
                    DEFAULT_STRIDE
                }
            },
            None => DEFAULT_STRIDE,
        };
        Self {
            stride,
            slice: std::env::var("ISOMESH_SLICE")
                .ok()
                .and_then(|v| v.parse::<u32>().ok())
                .unwrap_or(DEFAULT_SLICE)
                .min(SAMPLES_PER_AXIS - 1),
            touched: true,
            surface: true,
        }
    }
}

/// The spheres to draw, resolved once per change rather than per frame.
#[derive(Resource, Default)]
struct Overlay {
    untouched: Vec<(Vec3, f32)>,
    touched: Vec<(Vec3, f32)>,
    worst_sample: Vec3,
    worst_radius: f32,
    worst_vertex: Vec3,
    cell: f32,
}

/// The shells get their own group so they can be drawn through the translucent
/// surface without dragging the shared wireframe's bias along with them.
#[derive(Default, Reflect, GizmoConfigGroup)]
struct ShellGizmos;

/// The worst-miss marker gets a harder bias still, so the one sphere that names
/// the reported number is never lost behind the nine that do not.
#[derive(Default, Reflect, GizmoConfigGroup)]
struct MarkerGizmos;

/// The four surfaces, uploaded once.
#[derive(Resource)]
struct Surfaces(Vec<Handle<Mesh>>);

/// The plate behind the HUD text.
#[derive(Component)]
struct HudBackdrop;

/// Size of that plate, in logical pixels.
///
/// **A backdrop rather than brighter text, because the subject is line art over
/// the whole frame.** Every other example in this catalog draws its overlay in
/// one region and can put the HUD beside it; here the shells reach past the
/// domain, so orange strokes cross the upper left wherever the camera sits, and
/// white 13px text over moving orange line art is unreadable at the 900px a GIF
/// is scaled to. Measured on the frames: without this the first still was
/// illegible from the second line down.
///
/// Sized against what the HUD actually emits rather than guessed — 32 lines at
/// 13px, nothing wider than 60 characters, and 13px of this font measures about
/// 8px per character. Both numbers are held to: `report` says so, and a line
/// that outgrows the plate lands back on the line art.
const BACKDROP: Vec2 = Vec2::new(500.0, 524.0);

fn main() {
    // Measured **before** any plugin, and that ordering is load-bearing.
    // `add_plugins(DefaultPlugins)` builds `WinitPlugin`'s event loop, which
    // panics outright where there is no display — so anything that ran after it
    // would be unreachable on exactly the machines where the terminal is the
    // only output. Doing it first also means every system can take a plain
    // `Res<Measured>` rather than an `Option<Res<_>>` guarding a state that
    // cannot happen.
    let measured = measure_everything();

    App::new()
        .add_plugins(DefaultPlugins.set(WindowPlugin {
            primary_window: Some(Window {
                title: "isomesh — E-309 untouched tangency spheres".into(),
                ..default()
            }),
            ..default()
        }))
        .add_plugins(CommonPlugin)
        .init_gizmo_group::<ShellGizmos>()
        .init_gizmo_group::<MarkerGizmos>()
        .insert_resource(measured)
        .init_resource::<Stage>()
        .init_resource::<Draw>()
        .init_resource::<Overlay>()
        .add_systems(Startup, setup)
        .add_systems(
            Update,
            (controls, advance, rebuild, frame_camera, draw, report).chain(),
        )
        .run();
}

fn setup(
    mut commands: Commands,
    mut meshes: ResMut<Assets<Mesh>>,
    mut materials: ResMut<Assets<StandardMaterial>>,
    mut gizmo_config: ResMut<GizmoConfigStore>,
    mut camera: Query<&mut OrbitCamera>,
    measured: Res<Measured>,
) {
    for mut orbit in &mut camera {
        orbit.yaw = VIEW_YAW;
        orbit.pitch = VIEW_PITCH;
        orbit.radius = VIEW_RADIUS;
    }

    let (shells, _) = gizmo_config.config_mut::<ShellGizmos>();
    shells.line.width = 1.4;
    shells.depth_bias = -0.3;

    let (marks, _) = gizmo_config.config_mut::<MarkerGizmos>();
    marks.line.width = 3.0;
    marks.depth_bias = -0.9;

    // Translucent and double-sided, because the spheres are centred on samples
    // that are mostly *outside* the surface but their shells pass straight
    // through it — an opaque surface would hide the tangency it is being judged
    // on. The same reason `critical_cells` is translucent.
    let material = materials.add(StandardMaterial {
        base_color: Color::srgba(0.74, 0.78, 0.84, 0.16),
        perceptual_roughness: 0.5,
        alpha_mode: AlphaMode::Blend,
        double_sided: true,
        cull_mode: None,
        ..default()
    });

    let handles: Vec<Handle<Mesh>> = measured
        .combos
        .iter()
        .map(|combo| meshes.add(to_mesh(&combo.buffer)))
        .collect();
    commands.spawn((
        Mesh3d(handles[0].clone()),
        MeshMaterial3d(material),
        DemoMesh,
    ));
    commands.insert_resource(Surfaces(handles));

    // Both fields share the `[-2, 2]^3` compact domain, so this is set once.
    commands.spawn(DemoDomain {
        min: Vec3::splat(-2.0),
        max: Vec3::splat(2.0),
    });

    // Behind the harness's HUD text. `GlobalZIndex(-1)` rather than spawn order,
    // because the text entity belongs to `CommonPlugin` and was spawned first —
    // relying on order would put this plate on top of the numbers it exists to
    // make readable. Absolutely positioned to match `spawn_hud`'s 10/12 offset,
    // inset a little so the plate has a margin around the glyphs.
    commands.spawn((
        Node {
            position_type: PositionType::Absolute,
            top: Val::Px(4.0),
            left: Val::Px(4.0),
            width: Val::Px(BACKDROP.x),
            height: Val::Px(BACKDROP.y),
            ..default()
        },
        BackgroundColor(Color::srgba(0.04, 0.05, 0.07, 0.78)),
        GlobalZIndex(-1),
        HudBackdrop,
    ));
}

/// The `f64` extraction as a Bevy mesh.
///
/// Cast rather than re-extracted in `f32`: the counts on the HUD are `f64`
/// counts over `f64` vertices, so the mesh the picture is drawn from has to be
/// the one they were computed on.
fn to_mesh(buffer: &MeshBuffer<f64>) -> Mesh {
    let mut builder = MeshBuilder::new();
    for i in 0..buffer.positions.len() {
        let (Some(p), Some(n)) = (buffer.positions.get(i), buffer.normals.get(i)) else {
            continue;
        };
        builder.vertex(
            [p[0] as f32, p[1] as f32, p[2] as f32],
            [n[0] as f32, n[1] as f32, n[2] as f32],
        );
    }
    for t in buffer.indices.as_chunks::<3>().0 {
        let (Some(a), Some(b), Some(c)) = (t.first(), t.get(1), t.get(2)) else {
            continue;
        };
        builder.triangle(*a, *b, *c);
    }
    builder.into_mesh()
}

/// The keyboard. Ignored under capture, which drives itself.
fn controls(
    keys: Res<ButtonInput<KeyCode>>,
    capture: Res<Capture>,
    mut stage: ResMut<Stage>,
    mut draw: ResMut<Draw>,
) {
    if capture.is_active() {
        return;
    }
    for (key, field) in [(KeyCode::Digit1, 0), (KeyCode::Digit2, 1)] {
        if keys.just_pressed(key) {
            stage.field = field;
        }
    }
    if keys.just_pressed(KeyCode::KeyM) {
        stage.extractor = 1 - stage.extractor;
        stage.phase = 0.0;
    }
    if keys.just_pressed(KeyCode::BracketLeft) {
        draw.stride = (draw.stride + STRIDES.len() - 1) % STRIDES.len();
    }
    if keys.just_pressed(KeyCode::BracketRight) {
        draw.stride = (draw.stride + 1) % STRIDES.len();
    }
    if keys.just_pressed(KeyCode::Comma) {
        draw.slice = draw.slice.saturating_sub(SLICE_STEP);
    }
    if keys.just_pressed(KeyCode::Period) {
        draw.slice = (draw.slice + SLICE_STEP).min(SAMPLES_PER_AXIS - 1);
    }
    if keys.just_pressed(KeyCode::KeyT) {
        draw.touched = !draw.touched;
    }
    if keys.just_pressed(KeyCode::KeyH) {
        draw.surface = !draw.surface;
    }
}

/// Decide which row is on screen this frame.
///
/// Under capture both the field and the extractor come off the captured-frame
/// counter, so a clip of any length shows the flip that is the finding.
/// Interactively the extractor alternates on a timer for the same reason — an
/// example whose subject only changes on a keypress photographs as a still.
fn advance(
    time: Res<Time>,
    capture: Res<Capture>,
    flags: Res<ViewFlags>,
    mut stage: ResMut<Stage>,
) {
    if capture.is_active() {
        let step = capture.taken / CAPTURE_FRAMES_PER_STAGE;
        // `ISOMESH_FIELD` pins the field and leaves the flips, which is what
        // makes a one-field clip reachable from a command line.
        match std::env::var("ISOMESH_FIELD")
            .ok()
            .and_then(|v| v.parse::<usize>().ok())
        {
            Some(field) => {
                stage.field = field.min(FIELDS.len() - 1);
                stage.extractor = (step % 2) as usize;
            }
            None => {
                let cycle = step % (FLIPS_PER_FIELD * FIELDS.len() as u32);
                stage.field = (cycle / FLIPS_PER_FIELD) as usize;
                stage.extractor = (cycle % 2) as usize;
            }
        }
        return;
    }
    if flags.paused {
        return;
    }
    stage.phase += time.delta_secs();
    if stage.phase >= STAGE_SECONDS {
        stage.phase = 0.0;
        stage.extractor = 1 - stage.extractor;
    }
}

/// Swap the surface and resolve the spheres — only when the picture changed.
fn rebuild(
    stage: Res<Stage>,
    draw: Res<Draw>,
    measured: Res<Measured>,
    surfaces: Res<Surfaces>,
    mut overlay: ResMut<Overlay>,
    mut query: Query<&mut Mesh3d, With<DemoMesh>>,
    mut last: Local<Option<(usize, usize, usize, u32)>>,
) {
    let key = (stage.field, stage.extractor, draw.stride, draw.slice);
    if *last == Some(key) {
        return;
    }
    *last = Some(key);

    let handle = surfaces.0[combo_index(stage.field, stage.extractor)].clone();
    for mut mesh in &mut query {
        mesh.0 = handle.clone();
    }
    *overlay = resolve(
        measured.combo(stage.field, stage.extractor),
        measured.lattice(stage.field),
        &draw,
    );
}

/// The drawn subset of one row's verdict.
///
/// A sphere's colour is `untouched_mask[index]` — the *same* `bool` the census
/// counted — so the picture and the number cannot come apart. The subsample
/// chooses which spheres appear, never what they mean.
///
/// # The slice is a plane of constant `y`, and it is deliberately not `y = 0`
///
/// On `thin_plate` the untouched set has a shape, and it is a shape in `x` and
/// `z`: the plate spans `[−1, 1]` in both, its top and bottom faces only reach
/// `±0.9375` because that is the last *inside* sample, and Marching Cubes puts
/// no vertex on the rim's corners at all. A sample outside the footprint has to
/// reach a rim or a corner its mesh never covered. Three quarters of a
/// `[−2, 2]²` square lies outside its middle quarter — `1 − (33/65)² = 74.2%`,
/// against the measured 71.7%, so that is what the census is counting.
///
/// So the slice has to be a plane of constant `y`. A plane of constant `z` cuts
/// *across* the annulus: measured on the frames, `z = 0` at stride 16 draws 8
/// untouched of 25 — 32% against a census of 71.7% — and puts them at the edge
/// of the frame where they cannot be counted.
///
/// **And `y = 0` is the worst plane of its own family**, which is the part that
/// had to be measured rather than reasoned. The plate is centred there and
/// `THICKNESS_IN_CELLS = 0.4`, so `y = 0` is the only sample plane *inside* the
/// plate — and every crossing edge leaving it in `x` or `z` roots at `t = 1`,
/// which puts a Marching Cubes vertex exactly on the rim at `y = 0`. Those
/// vertices are coplanar with the samples and touch the outer ones *exactly*, so
/// `y = 0` at stride 16 also draws only 8 untouched of 25, and the eight are the
/// four plate corners and the four domain corners. One plane up, at
/// [`DEFAULT_SLICE`], they are no longer coplanar and the annulus appears: **24
/// of 25** at stride 16, **8 of 9** at the default stride 32.
///
/// So the two candidate planes bracket the census rather than land on it — 32%
/// one side, 89% the other, against 71.7%. **A coarse lattice over a compact
/// object over-weights the periphery and there is no stride that fixes it**: at
/// stride 32 the drawn `x` values are `{−2, 0, 2}` and only one of the three is
/// over the plate at all. That is why the HUD prints the drawn count and the
/// census count as two separate numbers and says which is which, and why the
/// colour is read out of the census's own `bool` rather than recomputed on the
/// subsample.
///
/// The radii sort themselves out on this plane too. At `y = 0.5` a sample over
/// the footprint asserts a sphere of radius `0.4875` — the one cyan shell in the
/// default frame, resting tangent on the plate — while a sample past the rim
/// asserts one of `1.1` to `1.5`. The samples with something to show are exactly
/// the ones that get a visible shell.
fn resolve(combo: &Combo, lattice: &Lattice, draw: &Draw) -> Overlay {
    let stride = STRIDES[draw.stride];
    let j = (draw.slice as usize).min(lattice.n - 1);
    let mut overlay = Overlay {
        cell: lattice.h as f32,
        ..default()
    };
    for k in (0..lattice.n).step_by(stride) {
        for i in (0..lattice.n).step_by(stride) {
            let index = lattice.index(i, j, k);
            let p = lattice.at(i, j, k);
            let centre = Vec3::new(p[0] as f32, p[1] as f32, p[2] as f32);
            let radius = lattice.abs_d[index] as f32;
            if combo.scan.untouched_mask[index] {
                overlay.untouched.push((centre, radius));
            } else {
                overlay.touched.push((centre, radius));
            }
        }
    }
    let worst = &combo.scan;
    overlay.worst_sample = Vec3::new(
        worst.sample[0] as f32,
        worst.sample[1] as f32,
        worst.sample[2] as f32,
    );
    overlay.worst_radius = worst.radius as f32;
    overlay.worst_vertex = Vec3::new(
        worst.nearest[0] as f32,
        worst.nearest[1] as f32,
        worst.nearest[2] as f32,
    );
    overlay
}

/// Keep the subject off-centre and clear of the HUD.
fn frame_camera(flags: Res<ViewFlags>, mut camera: Query<&mut OrbitCamera>) {
    for mut orbit in &mut camera {
        if flags.paused {
            continue;
        }
        // The camera's own basis, from the same yaw/pitch `orbit_camera` builds
        // its transform from, so the offset is one screen-space nudge however
        // far `ISOMESH_SPIN` has turned. The eye sits at `focus + dir * radius`,
        // so a focus moved along `−right` puts the subject right of centre.
        let dir = Vec3::new(
            orbit.yaw.cos() * orbit.pitch.cos(),
            orbit.pitch.sin(),
            orbit.yaw.sin() * orbit.pitch.cos(),
        );
        let forward = -dir;
        let right = forward.cross(Vec3::Y).normalize_or_zero();
        let up = right.cross(forward).normalize_or_zero();
        orbit.focus =
            -right * (SUBJECT_OFFSET.x * VIEW_RADIUS) + up * (SUBJECT_OFFSET.y * VIEW_RADIUS);
        orbit.radius = VIEW_RADIUS;
    }
}

/// Draw the shells and the worst-miss marker.
fn draw(
    overlay: Res<Overlay>,
    settings: Res<Draw>,
    flags: Res<ViewFlags>,
    mut visibility: Query<&mut Visibility, With<DemoMesh>>,
    mut backdrop: Query<&mut Visibility, (With<HudBackdrop>, Without<DemoMesh>)>,
    mut shells: Gizmos<ShellGizmos>,
    mut marks: Gizmos<MarkerGizmos>,
) {
    // Written only when it differs. Bevy's visibility propagation is
    // change-driven, so an unconditional write turns a toggle nobody pressed
    // into per-frame work on every descendant.
    let wanted = if settings.surface {
        Visibility::Inherited
    } else {
        Visibility::Hidden
    };
    for mut visible in &mut visibility {
        if *visible != wanted {
            *visible = wanted;
        }
    }

    // `nohud` has to take the plate with it, or the harness's own toggle leaves
    // a dark rectangle over a quarter of the frame with nothing written on it.
    let wanted = if flags.hud {
        Visibility::Inherited
    } else {
        Visibility::Hidden
    };
    for mut visible in &mut backdrop {
        if *visible != wanted {
            *visible = wanted;
        }
    }

    if settings.touched {
        for (centre, radius) in &overlay.touched {
            shells
                .sphere(Isometry3d::from_translation(*centre), *radius, TOUCHED)
                .resolution(20);
        }
    }
    for (centre, radius) in &overlay.untouched {
        shells
            .sphere(Isometry3d::from_translation(*centre), *radius, UNTOUCHED)
            .resolution(20);
        // The centre, so a sample whose sphere is too small to see is still
        // visibly a sample rather than a gap in the lattice.
        cross(&mut shells, *centre, overlay.cell * 0.30, UNTOUCHED);
    }

    // The worst miss, always drawn, wherever in the lattice it is.
    //
    // **The crosshair is a locator at a fixed eight cells, and that is stated
    // rather than implied.** The worst miss on all four rows is a witness whose
    // sphere has `|d| = 0` — a sample sitting exactly on the surface — so its
    // tangency sphere is a point and there is nothing to draw at the reported
    // radius. The measured quantity is the green segment: `0.088` world units on
    // `box_exact`/MC, which is `1.4` cells in a domain 64 cells across and would
    // be three pixels at 900px if the eye were not led to it first. So the
    // crosshair says *where*, the segment is the *what*, and the HUD carries the
    // number. The sphere is still drawn at the true `|d|` for the rows where
    // that is nonzero.
    marks
        .sphere(
            Isometry3d::from_translation(overlay.worst_sample),
            overlay.worst_radius,
            WORST,
        )
        .resolution(32);
    cross(&mut marks, overlay.worst_sample, overlay.cell * 8.0, WORST);
    marks.line(overlay.worst_sample, overlay.worst_vertex, NEAREST);
    marks
        .sphere(
            Isometry3d::from_translation(overlay.worst_vertex),
            overlay.cell * 0.35,
            NEAREST,
        )
        .resolution(12);
}

/// Three axis-aligned segments through a point.
fn cross<T: GizmoConfigGroup>(gizmos: &mut Gizmos<T>, at: Vec3, arm: f32, colour: Color) {
    for axis in [Vec3::X, Vec3::Y, Vec3::Z] {
        gizmos.line(at - axis * arm, at + axis * arm, colour);
    }
}

/// The HUD. The numbers are the demo.
fn report(
    stage: Res<Stage>,
    draw: Res<Draw>,
    measured: Res<Measured>,
    overlay: Res<Overlay>,
    mut stats: ResMut<DemoStats>,
) {
    let combo = measured.combo(stage.field, stage.extractor);
    let field = FIELDS[stage.field];
    let extractor = EXTRACTORS[stage.extractor];
    let scan = &combo.scan;

    // ASCII only in every string that reaches the HUD. Measured: Bevy's default
    // font has no glyph for `U+2014 EM DASH` and draws a tofu box, which put a
    // small black rectangle in the middle of the title and of the `drawn` line
    // in the first stills. The doc comments above are free to use it; the screen
    // is not.
    stats.title =
        format!("E-309  untouched tangency spheres - P-51 C3 / M-355  {SAMPLES_PER_AXIS}^3");
    stats.vertices = combo.buffer.vertex_count();
    stats.triangles = combo.buffer.triangle_count();
    stats.extract_ms = combo.extract_ms;

    let mc = measured.combo(stage.field, 0);
    let dc = measured.combo(stage.field, 1);
    let drawn = overlay.untouched.len() + overlay.touched.len();
    let stride = STRIDES[draw.stride];

    // The closed form M-355 names for this row's worst miss, if it names one,
    // with the delta this run measured against it. `bit for bit` and
    // `+ eps^2*sqrt(2)/8` are different statements and a six-decimal CSV cannot
    // separate them, so the HUD says which one it is.
    let (headline, closed) = match measured
        .check
        .identity(combo_index(stage.field, stage.extractor))
    {
        Some(id) if id.delta == 0.0 => (
            format!("= {} exactly", id.label),
            format!("{}: exact in f64, bit for bit", id.label),
        ),
        Some(id) => match id.predicted {
            Some(prediction) => (
                format!("= {} + {:.1e}", id.label, id.delta),
                format!(
                    "{} + eps^2*sqrt2/8 (ToCell clamp, {:.4}x)",
                    id.label,
                    id.delta / prediction
                ),
            ),
            None => (
                format!("= {} to {:.1e}", id.label, id.delta),
                format!("{}: {:.1e} cells off, unexplained", id.label, id.delta),
            ),
        },
        None => (
            String::from("no closed form named"),
            String::from("M-355 names no closed form for this row"),
        ),
    };

    // Every line kept inside 60 characters. That is not a style rule: the
    // backdrop behind this text is sized in pixels, and a line that runs past
    // it lands back on the orange line art and stops being readable at 900px.
    stats.extra = vec![
        format!("    field  {field:<12}   [M]  {extractor}"),
        String::new(),
        format!(
            "{:>9} untouched spheres  {:>5.2}% of {}",
            scan.untouched,
            combo.per_1k() / 10.0,
            combo.samples
        ),
        format!(
            "{:>9.4} untouched per 1k   gate {THRESHOLD_CELLS} cells",
            combo.per_1k()
        ),
        format!("{:>9.6} worst miss, cells  {headline}", scan.worst_cells),
        format!(
            "{:>9.4} MC/DC on {field:<11} {:.2} vs {:.2} /1k",
            measured.mc_over_dc(stage.field),
            mc.per_1k(),
            dc.per_1k()
        ),
        String::new(),
        format!(
            // "every 32th" is what an ordinal suffix does when it is glued to a
            // number, so the stride is spelled as a stride.
            "    drawn  stride {stride} over the y={} slice",
            draw.slice
        ),
        format!(
            "           {drawn} spheres, {} untouched - a subsample,",
            overlay.untouched.len()
        ),
        String::from("           not the count; colour is the census verdict"),
        format!(
            "    touch  exhaustive: all {} vertices, {:.0} ms",
            combo.buffer.vertex_count(),
            combo.scan_ms
        ),
        String::new(),
        format!(
            "    worst  sample ({:.4}, {:.4}, {:.4})",
            scan.sample[0], scan.sample[1], scan.sample[2]
        ),
        format!(
            "           |d| = {:.6}, vertex at {:.6}",
            scan.radius, scan.distance
        ),
        format!(
            "           miss {:.6} = {:.6} cells",
            scan.distance - scan.radius,
            scan.worst_cells
        ),
        format!("           {closed}"),
        String::new(),
        format!(
            "   ledger  p-51.csv: {}/{} columns held to the digit",
            measured.check.held, measured.check.total
        ),
        format!(
            "           {:.1}%-{:.1}% over its 15 rows{}",
            measured.check.range_pct.0,
            measured.check.range_pct.1,
            if measured.check.range_agrees {
                String::new()
            } else {
                format!(" (M-355 says {:.1}%)", FINDING_RANGE_PCT.0)
            }
        ),
        String::new(),
        format!(
            "  [1|2] field  [ [ | ] ] stride {stride}  [,|.] slice {}",
            draw.slice
        ),
        format!(
            "  [T] touched shells {}  [H] surface {}",
            on_off(draw.touched),
            on_off(draw.surface)
        ),
    ];
}

fn on_off(value: bool) -> &'static str {
    if value { "on" } else { "off" }
}
