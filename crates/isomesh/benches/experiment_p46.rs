//! **P-46 — repair the sign lattice, not the mesh, by minimal value
//! perturbation.**
//!
//! Ticket: R-040a. Pre-registered before this harness existed.
//!
//! ```bash
//! cargo bench --bench experiment_p46
//! ```
//!
//! Writes `docs/experiments/p-46.csv`.
//!
//! # What P-41 left behind
//!
//! P-41 (M-338) measured a bijection at 65 samples per axis: on every field with
//! a non-zero census, the number of critical cells, the number of non-manifold
//! vertices, and the number of critical cells hosting one were the *same*
//! number — 602 on `noise_cavity`, 141 on `gyroid`, 58 on `fbm_terrain` — and
//! co-location was 2442/2442. That made well-composedness **necessary**. It said
//! nothing about sufficiency, which is what this measures.
//!
//! # The repair, and why it is a perturbation rather than a flip
//!
//! Latecki's characterisation is about the *sign* lattice, so the naive repair is
//! to flip a sign. That is unbounded: the corner whose sign is flipped may be
//! the one furthest from the surface, and the zero crossing then moves by up to
//! a whole cell.
//!
//! The registered repair instead moves the **grey-level function** by the least
//! it can:
//!
//! 1. Find the critical cell's corner of smallest `|value|` — the corner the
//!    surface already passes closest to.
//! 2. Move that value across zero by the smallest representable step.
//!
//! Step 2 is asymmetric because the sign test is, and that asymmetry is exactly
//! `cube.rs::is_inside`'s: inside is `value < 0.0`, so
//!
//! - inside → outside costs a landing at `0.0`, which is the smallest step there
//!   is, because `is_inside(0.0)` is already false;
//! - outside → inside costs a landing at `-f64::from_bits(1)`, the negative
//!   subnormal nearest zero, which is the largest strictly-negative `f64`.
//!
//! Both landings put the corner value at (or one subnormal past) zero, so
//! `edge_crossing(a, b) = a / (a - b)` places the crossing *at that grid corner*.
//! The surface does not move by a cell; it moves to the nearest corner it was
//! already near. That is the whole reason the cheapest corner is chosen, and it
//! is what clause three prices.
//!
//! This is Boutry, Géraud & Najman's self-dual repair (*A Tutorial on
//! Well-Composedness*, JMIV 2018, `10.1007/s10851-017-0769-6`) applied to the
//! grey-level function rather than to the binary set. Self-dual matters here:
//! the rule treats a critical configuration and its complement identically —
//! `classify` marks both, and `across_zero` moves whichever corner is cheapest
//! in whichever direction it happens to need — so the repair commits to no
//! foreground/background asymmetry that the extractors do not already have.
//!
//! # Termination is a property of the rule, not an addition to it
//!
//! A perturbed corner has `|value|` of zero or one subnormal, which makes it the
//! *cheapest* corner of every cell it touches. Re-selecting it would flip it
//! straight back, and the sweep would oscillate forever without ever being
//! wrong about any single step. So a corner already moved this run is excluded
//! from selection. That is not a second path or a fallback: it is the
//! restatement of "move it across zero" as an idempotent operation, which it has
//! to be for "sweep until fixpoint" to be a definition at all.
//!
//! When every corner of a still-critical cell has already been moved there is
//! nothing the rule permits, and that cell is counted in `stuck_cells` and
//! reported rather than repaired by some other means. **It happens.** On
//! `noise_cavity` all 118 surviving critical cells are corner-exhausted, which
//! makes the registered rule a dead end there rather than a slow one.
//!
//! # What the primary source says about exactly this, read after measuring
//!
//! The tutorial's §7.2 is specific, and it names the property this rule does not
//! have. Of the in-place *n*-D gray-level repair — Boutry et al. [28], *How to
//! make images well-composed in nD without interpolation* — it says: *"this
//! method is based on an **'increasing' procedure which avoids oscillations and
//! ensures convergence in linear time** w.r.t. the size of the domain of the
//! image."*
//!
//! The registered rule is **not** increasing. It moves the cheapest corner in
//! whichever direction that corner happens to need, and both directions occur —
//! that is what makes it minimal, and it is also exactly the non-monotonicity
//! whose absence the source credits for convergence. The never-re-move exclusion
//! above is this harness paying for that missing monotonicity by hand, and
//! `residual_exhausted_cells` is the bill.
//!
//! Two further numbers from the same section are worth having beside the result.
//! Siqueira et al.'s randomized binary repair carries *"a theoretical bound
//! [that] ensures that the maximal number of new critical configurations which
//! will appear during the elimination of the m initial configurations is lower
//! than or equal to m/2"* — a convergent geometric series, at most `2m`
//! modifications. This rule does not inherit it: `noise_cavity` starts at
//! `m = 602` and spends 5190 corner moves without converging, which is `8.6·m`.
//! And the section opens by stating that *"no method able to ... make gray-level
//! images CWC in 3D exist nowadays to the best of our knowledge"*, so the target
//! is legitimate — DWC and CWC are equivalent in 3D (the tutorial's Table 1), so
//! removing every critical configuration *does* buy a 2-manifold boundary — but
//! the mechanism that reaches it is [28]'s monotone one, not this one.
//!
//! None of this amends the prediction. The rule measured here is the rule that
//! was registered, implemented as written; this section records what the source
//! says about why it behaves as it does.
//!
//! # How the perturbed lattice is meshed without touching the crate
//!
//! [`Lattice`] is the sample grid, and [`Perturbed`] is an [`Sdf`] over it.
//! `dual.rs` calls `sample` in exactly one place — a single pass over the whole
//! grid at `origin + cell_size · (x, y, z)` — so a field that answers at grid
//! points is all either dual extractor ever asks for. [`Perturbed::sample`]
//! inverts that arithmetic, and **panics** if it is ever handed a point that is
//! not a grid point, rather than rounding to the nearest one and returning a
//! plausible number. A silent round would make every result here unfalsifiable.
//!
//! `gradient` forwards to the true field. It must: gradients are evaluated at
//! final *vertex* positions, which are off-grid by construction, and there is no
//! perturbed value there to differentiate. This is also the right answer rather
//! than merely the available one — the perturbation is one subnormal of one
//! sample, so the field's own derivative is still the limit it is approximating.
//! `wrapper_faithful` is the instrument check: the "before" mesh is extracted
//! through [`Perturbed`] over the *unperturbed* lattice and its
//! [`mesh_hash`](isomesh::validate::mesh_hash) compared against the same mesh
//! extracted from the field directly. Equal hashes mean the wrapper is not
//! what moved anything.
//!
//! # Clause three measures against the field the caller asked for
//!
//! `hausdorff_before` and `hausdorff_after` both call
//! [`accuracy`](isomesh::validate::accuracy) against the **original, unrepaired
//! field**. Measuring the repaired mesh against the perturbed samples would
//! compare the repair with itself and could only ever look good.
//!
//! # This is a capability claim, and only that
//!
//! The repair changes the mesh. It is not a speedup, it is not offered as one,
//! and golden hashes are beside the point here: nothing in `crates/isomesh/src`
//! is touched, and the perturbed grid lives in this bench. Whether the crate
//! should adopt it is a separate decision that clause three prices and this
//! experiment does not make.

mod common;

use isomesh::dual_contouring::DualContouring;
use isomesh::fields::ReferenceField;
use isomesh::surface_nets::SurfaceNets;
use isomesh::validate::{AccuracyConfig, ValidateConfig, accuracy, mesh_hash, validate_features};
use isomesh::{MeshBuffer, RuntimeShape3, Sdf};

/// The registered resolution. `65` samples span `64` cells per axis.
const SAMPLES: u32 = 65;

/// Sweep cap. Clause two claims a fixpoint in two; this is high enough that
/// hitting it is a cascade rather than a tight budget, and finite so that a
/// cascade is a number instead of a hung bench.
const MAX_SWEEPS: u32 = 64;

/// Inside, the way the extractors decide it.
///
/// `cube.rs::is_inside`: strictly negative. `-0.0` is **not** less than `0.0`,
/// so a negative zero is outside here exactly as it is there.
fn is_inside(value: f64) -> bool {
    value < 0.0
}

/// Move a value across zero by the smallest representable step.
///
/// Asymmetric because [`is_inside`] is. Landing on `0.0` is already outside, so
/// that is the cheapest possible exit from inside; coming the other way needs a
/// value *strictly* below zero, and the largest of those is the negative
/// subnormal nearest zero.
fn across_zero(value: f64) -> f64 {
    if is_inside(value) {
        0.0
    } else {
        -f64::from_bits(1)
    }
}

// ─── the 256-byte classification, enumerated (P-41's, unchanged) ────────────

/// Which of the 256 possible cell sign bytes host each critical configuration.
///
/// Corner bit layout: bit `i` is corner `(x, y, z)` with `i = x + 2y + 4z`, so
/// two corners are cell-diagonal exactly when `i ^ j == 0b111`, and two corners
/// of the face fixing axis `a` are face-diagonal exactly when
/// `i ^ j == 0b111 ^ (1 << a)`.
///
/// Enumerated from the definitions, never transcribed — `CLAUDE.md` rule 5.
/// P-41 cross-checked the result against an independently written classifier:
/// 120 of 256 bytes are 2D-critical, 8 are 3D-critical, and the classes are
/// disjoint.
struct Critical {
    /// A checkerboard `2×2` face: the inside pair shares only a cell edge.
    two_d: [bool; 256],
    /// A main-diagonal inside pair, or its complement: the pair shares only a
    /// cell vertex.
    three_d: [bool; 256],
}

impl Critical {
    fn hosts(&self, byte: u32) -> bool {
        self.two_d[byte as usize] || self.three_d[byte as usize]
    }
}

/// The inside corners of `byte`, and how many there are.
fn inside_corners(byte: u32) -> ([u32; 8], usize) {
    let mut out = [0u32; 8];
    let mut n = 0;
    for corner in 0..8u32 {
        if (byte >> corner) & 1 == 1 {
            out[n] = corner;
            n += 1;
        }
    }
    (out, n)
}

/// Exactly two inside corners, differing in all three coordinates.
fn is_vertex_diagonal_pair(byte: u32) -> bool {
    let (corners, n) = inside_corners(byte);
    n == 2 && (corners[0] ^ corners[1]) == 0b111
}

/// Some `2×2` face of the cell is a checkerboard.
fn has_checkerboard_face(byte: u32) -> bool {
    for axis in 0..3u32 {
        let diagonal = 0b111 ^ (1 << axis);
        for side in 0..2u32 {
            let mut inside = [0u32; 4];
            let mut n = 0;
            for corner in 0..8u32 {
                if (corner >> axis) & 1 == side && (byte >> corner) & 1 == 1 {
                    inside[n] = corner;
                    n += 1;
                }
            }
            if n == 2 && (inside[0] ^ inside[1]) == diagonal {
                return true;
            }
        }
    }
    false
}

/// Decide all 256 sign bytes from the definitions, once.
fn classify() -> Critical {
    let mut two_d = [false; 256];
    let mut three_d = [false; 256];
    for byte in 0..256u32 {
        two_d[byte as usize] = has_checkerboard_face(byte);
        three_d[byte as usize] =
            is_vertex_diagonal_pair(byte) || is_vertex_diagonal_pair(!byte & 0xFF);
    }
    Critical { two_d, three_d }
}

// ─── the grid (P-41's, unchanged) ───────────────────────────────────────────

/// The grid a field is sampled and meshed on.
struct Grid {
    /// World position of sample `[0, 0, 0]`.
    origin: [f64; 3],
    /// Spacing.
    cell_size: f64,
    /// Samples per axis.
    samples: u32,
    /// Cells per axis: `samples - 1`.
    cells: u32,
}

impl Grid {
    /// `i = x + y·sx + z·sx·sy`, the crate's order.
    fn sample_index(&self, x: usize, y: usize, z: usize) -> usize {
        let n = self.samples as usize;
        x + y * n + z * n * n
    }

    /// `i` over cells rather than samples, same order.
    fn cell_index(&self, cell: [usize; 3]) -> usize {
        let c = self.cells as usize;
        cell[0] + cell[1] * c + cell[2] * c * c
    }

    /// The cell a dual vertex belongs to, by the floor rule P-41 established
    /// and measured: `Clamp::ToCell` insets by `(1 − ε)` so a Dual Contouring
    /// vertex is strictly interior, and P-41 confirmed the mapping is a
    /// bijection onto the sign-active cells with zero escapes.
    fn cell_of(&self, p: [f64; 3]) -> usize {
        let mut cell = [0usize; 3];
        let last = (self.cells - 1) as f64;
        for (axis, slot) in cell.iter_mut().enumerate() {
            let t = ((p[axis] - self.origin[axis]) / self.cell_size).floor();
            *slot = t.clamp(0.0, last) as usize;
        }
        self.cell_index(cell)
    }

    /// World position of sample `(x, y, z)`, by the same arithmetic `dual.rs`
    /// uses, so the two agree bit for bit.
    fn point(&self, x: usize, y: usize, z: usize) -> [f64; 3] {
        [
            self.origin[0] + self.cell_size * x as f64,
            self.origin[1] + self.cell_size * y as f64,
            self.origin[2] + self.cell_size * z as f64,
        ]
    }

    /// The eight sample indices of cell `(cx, cy, cz)`, in corner-bit order.
    fn corners(&self, cx: usize, cy: usize, cz: usize) -> [usize; 8] {
        let mut out = [0usize; 8];
        for (corner, slot) in out.iter_mut().enumerate() {
            *slot = self.sample_index(
                cx + (corner & 1),
                cy + ((corner >> 1) & 1),
                cz + ((corner >> 2) & 1),
            );
        }
        out
    }
}

// ─── the lattice, and an Sdf over it ────────────────────────────────────────

/// The sampled grid, plus which samples the repair has already moved.
struct Lattice {
    /// One value per sample, `x + y·n + z·n²`.
    values: Vec<f64>,
    /// One flag per sample, set when the repair has moved that value. A moved
    /// value has `|v|` of zero or one subnormal and would otherwise be
    /// re-selected forever.
    moved: Vec<bool>,
}

impl Lattice {
    /// Sample the field on the whole grid.
    fn sampled<F: Sdf<Scalar = f64>>(field: &F, grid: &Grid) -> Self {
        let n = grid.samples as usize;
        let mut values = vec![0.0f64; n * n * n];
        for z in 0..n {
            for y in 0..n {
                for x in 0..n {
                    values[grid.sample_index(x, y, z)] = field.sample(grid.point(x, y, z));
                }
            }
        }
        let moved = vec![false; values.len()];
        Self { values, moved }
    }

    /// The sign byte of cell `(cx, cy, cz)` under the current values.
    fn byte(&self, grid: &Grid, cx: usize, cy: usize, cz: usize) -> u32 {
        let mut byte = 0u32;
        for (corner, index) in grid.corners(cx, cy, cz).into_iter().enumerate() {
            if is_inside(self.values[index]) {
                byte |= 1 << corner;
            }
        }
        byte
    }
}

/// An [`Sdf`] that answers from a [`Lattice`] at grid points and forwards
/// gradients to the true field.
struct Perturbed<'a, F> {
    lattice: &'a Lattice,
    grid: &'a Grid,
    field: &'a F,
}

impl<F: Sdf<Scalar = f64>> Sdf for Perturbed<'_, F> {
    type Scalar = f64;

    /// # Panics
    ///
    /// If `p` is not a grid point. `dual.rs` samples only at grid points, so
    /// this cannot fire on the intended path — and if it ever did, rounding to
    /// the nearest sample and returning a plausible value would silently turn
    /// this whole experiment into a measurement of something else.
    fn sample(&self, p: [f64; 3]) -> f64 {
        let mut index = [0usize; 3];
        let last = (self.grid.samples - 1) as f64;
        // A grid point is an exact multiple of the spacing from the origin. The
        // tolerance is a millionth of a cell: far below any real off-grid probe
        // (the central-difference default steps by `DIFF_STEP·max(|p|, 1)`,
        // which is many orders larger) and far above `f64` rounding on this
        // arithmetic.
        let tolerance = 1e-6;
        for (axis, slot) in index.iter_mut().enumerate() {
            let t = (p[axis] - self.grid.origin[axis]) / self.grid.cell_size;
            let rounded = t.round();
            assert!(
                (t - rounded).abs() < tolerance && rounded >= 0.0 && rounded <= last,
                "Perturbed::sample was asked for an off-grid point: axis {axis} of {p:?} \
                 is {t} cells from the origin, which is not a sample. The perturbed \
                 lattice has no value there and guessing one would invalidate P-46."
            );
            *slot = rounded as usize;
        }
        self.lattice.values[self.grid.sample_index(index[0], index[1], index[2])]
    }

    #[inline]
    fn gradient(&self, p: [f64; 3]) -> [f64; 3] {
        self.field.gradient(p)
    }
}

// ─── the census and the repair ──────────────────────────────────────────────

/// Critical cells under the lattice's current signs.
struct Count {
    two_d: u64,
    three_d: u64,
    either: u64,
}

/// Count the critical cells, classifying each of the two configurations.
fn count_critical(lattice: &Lattice, grid: &Grid, table: &Critical) -> Count {
    let c = grid.cells as usize;
    let mut out = Count {
        two_d: 0,
        three_d: 0,
        either: 0,
    };
    for cz in 0..c {
        for cy in 0..c {
            for cx in 0..c {
                let byte = lattice.byte(grid, cx, cy, cz);
                let flat = table.two_d[byte as usize];
                let solid = table.three_d[byte as usize];
                out.two_d += u64::from(flat);
                out.three_d += u64::from(solid);
                out.either += u64::from(flat || solid);
            }
        }
    }
    out
}

/// The critical cells that survived the repair, described rather than counted.
///
/// The registration's falsifier reads *"any non-manifold output surviving the
/// repair ... would mean well-composedness is not sufficient"*. That inference
/// is only available if the repair actually reached a well-composed lattice, so
/// when it does not the survivors have to be enumerated rather than summarised
/// into a single failure.
struct Residual {
    /// One flag per cell, `true` when the cell is still critical.
    critical: Vec<bool>,
    /// `(sign byte, count)` for every surviving byte, descending by count then
    /// ascending by byte so the order is total.
    bytes: Vec<(u32, u64)>,
    /// Survivors whose every corner had already been moved. The registered rule
    /// has nothing left to spend on these: it is a dead end, not a slow
    /// convergence.
    exhausted: u64,
    two_d: u64,
    three_d: u64,
}

/// Enumerate the surviving critical cells.
fn residual(lattice: &Lattice, grid: &Grid, table: &Critical) -> Residual {
    let c = grid.cells as usize;
    let mut critical = vec![false; c * c * c];
    let mut histogram = [0u64; 256];
    let mut out = Residual {
        critical: Vec::new(),
        bytes: Vec::new(),
        exhausted: 0,
        two_d: 0,
        three_d: 0,
    };
    for cz in 0..c {
        for cy in 0..c {
            for cx in 0..c {
                let byte = lattice.byte(grid, cx, cy, cz);
                if !table.hosts(byte) {
                    continue;
                }
                critical[grid.cell_index([cx, cy, cz])] = true;
                histogram[byte as usize] += 1;
                out.two_d += u64::from(table.two_d[byte as usize]);
                out.three_d += u64::from(table.three_d[byte as usize]);
                if grid
                    .corners(cx, cy, cz)
                    .into_iter()
                    .all(|index| lattice.moved[index])
                {
                    out.exhausted += 1;
                }
            }
        }
    }
    out.critical = critical;
    for (byte, count) in histogram.into_iter().enumerate() {
        if count > 0 {
            out.bytes.push((byte as u32, count));
        }
    }
    out.bytes.sort_by(|a, b| b.1.cmp(&a.1).then(a.0.cmp(&b.0)));
    out
}

/// What one sweep did.
struct Swept {
    /// Corners moved across zero.
    flips: u64,
    /// Critical cells the rule could not touch because all eight corners had
    /// already been moved. Never repaired by other means; reported instead.
    stuck: u64,
}

/// One sweep in ascending cell index order, applying each move immediately.
///
/// Immediate application (Gauss–Seidel rather than Jacobi) is deliberate: a move
/// often un-criticalises the neighbouring cells the same sweep is about to
/// visit, which is what keeps the cascade short. It is order-dependent, and the
/// order is `cx` fastest then `cy` then `cz` — the crate's own index order, so
/// the result is a pure function of the field and the grid.
fn sweep(lattice: &mut Lattice, grid: &Grid, table: &Critical) -> Swept {
    let c = grid.cells as usize;
    let mut out = Swept { flips: 0, stuck: 0 };
    for cz in 0..c {
        for cy in 0..c {
            for cx in 0..c {
                if !table.hosts(lattice.byte(grid, cx, cy, cz)) {
                    continue;
                }
                // The corner the surface already passes closest to, among those
                // the rule has not already spent. Ties break to the lowest
                // corner index, so the choice is deterministic.
                let mut best: Option<(f64, usize)> = None;
                for index in grid.corners(cx, cy, cz) {
                    if lattice.moved[index] {
                        continue;
                    }
                    let magnitude = lattice.values[index].abs();
                    if best.is_none_or(|(m, _)| magnitude < m) {
                        best = Some((magnitude, index));
                    }
                }
                match best {
                    Some((_, index)) => {
                        lattice.values[index] = across_zero(lattice.values[index]);
                        lattice.moved[index] = true;
                        out.flips += 1;
                    }
                    None => out.stuck += 1,
                }
            }
        }
    }
    out
}

/// Sweep until the census is zero, or until [`MAX_SWEEPS`].
struct Repair {
    /// Sweeps actually run. Clause two claims at most two.
    sweeps: u32,
    /// Critical cells left when sweeping stopped.
    residual: u64,
    /// Critical cells after exactly two sweeps, which is the number clause two
    /// is decided on whatever the fixpoint turns out to cost.
    residual_after_two: u64,
    /// Corners moved, total.
    flips: u64,
    /// Cells the rule could not touch.
    stuck: u64,
}

fn repair(lattice: &mut Lattice, grid: &Grid, table: &Critical) -> Repair {
    let mut sweeps = 0u32;
    let mut flips = 0u64;
    let mut stuck = 0u64;
    let mut residual = count_critical(lattice, grid, table).either;
    let mut residual_after_two = residual;
    while residual > 0 && sweeps < MAX_SWEEPS {
        let done = sweep(lattice, grid, table);
        flips += done.flips;
        stuck += done.stuck;
        sweeps += 1;
        residual = count_critical(lattice, grid, table).either;
        if sweeps <= 2 {
            // A sweep over a clean lattice moves nothing, so once the census is
            // zero this is also the count "after two sweeps".
            residual_after_two = residual;
        }
    }
    Repair {
        sweeps,
        residual,
        residual_after_two,
        flips,
        stuck,
    }
}

// ─── measuring one arm ──────────────────────────────────────────────────────

/// Which dual extractor.
#[derive(Clone, Copy)]
enum Which {
    DualContouring,
    SurfaceNets,
}

impl Which {
    const fn name(self) -> &'static str {
        match self {
            Self::DualContouring => "dual_contouring",
            Self::SurfaceNets => "surface_nets",
        }
    }
}

/// The two dual extractors, in the order they are reported.
const ARMS: [Which; 2] = [Which::DualContouring, Which::SurfaceNets];

/// Extract with one of the two dual extractors into `out`.
fn extract<S: Sdf<Scalar = f64>>(
    which: Which,
    field: &S,
    shape: &RuntimeShape3,
    grid: &Grid,
    out: &mut MeshBuffer<f64>,
) {
    out.reset();
    match which {
        Which::DualContouring => DualContouring::<f64>::new()
            .extract(field, shape, grid.origin, grid.cell_size, out)
            .expect("dual contouring extraction"),
        Which::SurfaceNets => SurfaceNets::<f64>::new()
            .extract(field, shape, grid.origin, grid.cell_size, out)
            .expect("surface nets extraction"),
    }
}

/// Non-manifold counts, and the symmetric Hausdorff against the true field.
struct Measured {
    edges: u64,
    vertices: u64,
    hausdorff: f64,
    mesh_vertices: u64,
    mesh_triangles: u64,
    hash: u64,
    /// `accuracy` sampled both directions. A `false` here voids the Hausdorff.
    has_coverage: bool,
    /// Which vertices are non-manifold, so a survivor can be mapped back to the
    /// cell it sits in. From the crate's own `validate_features`, not a
    /// recomputation — P-41 established the two agree exactly.
    nm_vertex_ids: Vec<u32>,
}

/// Validate a mesh and measure it against `truth` — always the original field,
/// never the perturbed lattice.
fn assess<F: Sdf<Scalar = f64>>(
    mesh: &MeshBuffer<f64>,
    truth: &F,
    shape: &RuntimeShape3,
    grid: &Grid,
) -> Measured {
    let vcfg = ValidateConfig::from_cell_size(grid.cell_size).expect("positive cell size");
    let (report, features) = validate_features(&mesh.positions, &mesh.indices, &vcfg);
    let acfg = AccuracyConfig::from_cell_size(grid.cell_size).expect("positive cell size");
    let acc = accuracy(
        &mesh.positions,
        &mesh.indices,
        truth,
        shape,
        grid.origin,
        &acfg,
    )
    .expect("accuracy is measurable on this grid");
    Measured {
        edges: report.non_manifold_edges,
        vertices: report.non_manifold_vertices,
        hausdorff: acc.symmetric_hausdorff(),
        mesh_vertices: mesh.vertex_count() as u64,
        mesh_triangles: mesh.triangle_count() as u64,
        hash: mesh_hash(mesh),
        has_coverage: acc.has_coverage(),
        nm_vertex_ids: features.vertices,
    }
}

/// Everything one `(field, extractor)` pair produced.
struct Arm {
    field: &'static str,
    extractor: &'static str,
    critical_before: u64,
    critical_after: u64,
    critical_2d_before: u64,
    critical_3d_before: u64,
    sweeps: u32,
    residual_after_two: u64,
    flips: u64,
    stuck: u64,
    before: Measured,
    after: Measured,
    /// The "before" mesh through [`Perturbed`] over the unperturbed lattice has
    /// the same hash as the same mesh straight from the field.
    wrapper_faithful: bool,
    /// Samples the repair moved, as a fraction of all samples.
    moved_samples: u64,
    total_samples: u64,
    /// The surviving critical configurations, enumerated. Empty when the repair
    /// reached a well-composed lattice.
    survivors: Vec<(u32, u64)>,
    /// Survivors with every corner already moved: a dead end for the registered
    /// rule rather than slow progress.
    exhausted: u64,
    residual_2d: u64,
    residual_3d: u64,
    /// Surviving non-manifold vertices that sit in a surviving critical cell.
    /// Equal to `after.vertices` means P-41's bijection still holds and the
    /// residual non-manifoldness is fully accounted for by the residual census —
    /// which is a very different claim from "well-composedness is insufficient".
    residual_colocated: u64,
}

impl Arm {
    fn hausdorff_ratio(&self) -> f64 {
        self.after.hausdorff / self.before.hausdorff
    }
}

/// The three fields P-41 found a non-zero census on, which are the three the
/// registration names.
fn measure<F>(field_name: &'static str, field: &F, table: &Critical, which: Which) -> Arm
where
    F: ReferenceField + Sdf<Scalar = f64>,
{
    let (shape, origin, cell_size) = common::grid::<f64, _>(field, SAMPLES);
    let grid = Grid {
        origin,
        cell_size,
        samples: SAMPLES,
        cells: SAMPLES - 1,
    };

    // ── before ──────────────────────────────────────────────────────────────
    let mut lattice = Lattice::sampled(field, &grid);
    let census_before = count_critical(&lattice, &grid, table);

    let mut mesh = MeshBuffer::<f64>::new();
    extract(which, field, &shape, &grid, &mut mesh);
    let direct = mesh_hash(&mesh);
    let before = assess(&mesh, field, &shape, &grid);

    // The same mesh, but read through the wrapper over the untouched lattice.
    // Equal hashes are what license every "after" number below.
    extract(
        which,
        &Perturbed {
            lattice: &lattice,
            grid: &grid,
            field,
        },
        &shape,
        &grid,
        &mut mesh,
    );
    let wrapper_faithful = mesh_hash(&mesh) == direct;

    // ── repair, then after ──────────────────────────────────────────────────
    let repaired = repair(&mut lattice, &grid, table);
    extract(
        which,
        &Perturbed {
            lattice: &lattice,
            grid: &grid,
            field,
        },
        &shape,
        &grid,
        &mut mesh,
    );
    let after = assess(&mesh, field, &shape, &grid);
    let survivors = residual(&lattice, &grid, table);
    let residual_colocated = after
        .nm_vertex_ids
        .iter()
        .filter(|v| survivors.critical[grid.cell_of(mesh.positions[**v as usize])])
        .count() as u64;

    let moved_samples = lattice.moved.iter().filter(|m| **m).count() as u64;
    Arm {
        field: field_name,
        extractor: which.name(),
        critical_before: census_before.either,
        critical_after: repaired.residual,
        critical_2d_before: census_before.two_d,
        critical_3d_before: census_before.three_d,
        sweeps: repaired.sweeps,
        residual_after_two: repaired.residual_after_two,
        flips: repaired.flips,
        stuck: repaired.stuck,
        before,
        after,
        wrapper_faithful,
        moved_samples,
        total_samples: lattice.values.len() as u64,
        survivors: survivors.bytes,
        exhausted: survivors.exhausted,
        residual_2d: survivors.two_d,
        residual_3d: survivors.three_d,
        residual_colocated,
    }
}

// ─── main ───────────────────────────────────────────────────────────────────

fn main() {
    let prereg = isomesh::experiment!("P-46");

    common::experiment::run(prereg, |run| {
        let table = classify();

        let mut arms: Vec<Arm> = Vec::new();
        isomesh::for_each_reference_field!(f64, |name, field| {
            // No early return in here: the macro inlines its body once per
            // field, so a `return` would leave the sweep silently short.
            // `noise_cavity`, `gyroid` and `fbm_terrain` are the three fields
            // P-41 measured a non-zero census on, and the three the
            // registration names. The other five are the control below.
            if name == "noise_cavity" || name == "gyroid" || name == "fbm_terrain" {
                for which in ARMS {
                    arms.push(measure(name, &field, &table, which));
                }
            }
        });

        for arm in &arms {
            println!(
                "{:>13} {:>16}  critical {:>4} → {:<4} in {} sweep(s), {} flips  \
                 nm {:>3}e/{:>3}v → {:>3}e/{:>3}v  hausdorff {:.6e} → {:.6e} ({:.4}×)",
                arm.field,
                arm.extractor,
                arm.critical_before,
                arm.critical_after,
                arm.sweeps,
                arm.flips,
                arm.before.edges,
                arm.before.vertices,
                arm.after.edges,
                arm.after.vertices,
                arm.before.hausdorff,
                arm.after.hausdorff,
                arm.hausdorff_ratio(),
            );
            if !arm.survivors.is_empty() {
                // Every surviving byte, with the class it belongs to and how
                // many cells carry it. This is the enumeration the falsifier
                // asks for when the repair does not reach a fixpoint.
                let listed: Vec<String> = arm
                    .survivors
                    .iter()
                    .map(|(byte, count)| {
                        let class =
                            match (table.two_d[*byte as usize], table.three_d[*byte as usize]) {
                                (true, true) => "2d+3d",
                                (true, false) => "2d",
                                (false, true) => "3d",
                                (false, false) => "??",
                            };
                        format!("{byte:#04x}/{class}×{count}")
                    })
                    .collect();
                println!(
                    "{:>13} {:>16}  SURVIVING: {} distinct bytes, {} exhausted \
                     (all 8 corners spent), residual 2d {} / 3d {}, {} of {} \
                     residual nm vertices in residual critical cells\n{:>31}  {}",
                    "",
                    "",
                    arm.survivors.len(),
                    arm.exhausted,
                    arm.residual_2d,
                    arm.residual_3d,
                    arm.residual_colocated,
                    arm.after.vertices,
                    "",
                    listed.join(" "),
                );
            }
        }

        // The control the registration does not require and the mechanism does:
        // a field whose census is already zero must be left alone by a repair
        // that only ever touches critical cells. Printed rather than filed,
        // because the registered rows are the three fields named.
        println!();
        isomesh::for_each_reference_field!(f64, |name, field| {
            if name != "noise_cavity" && name != "gyroid" && name != "fbm_terrain" {
                let (_, origin, cell_size) = common::grid::<f64, _>(&field, SAMPLES);
                let grid = Grid {
                    origin,
                    cell_size,
                    samples: SAMPLES,
                    cells: SAMPLES - 1,
                };
                let mut lattice = Lattice::sampled(&field, &grid);
                let done = repair(&mut lattice, &grid, &table);
                println!(
                    "control  {name:>14}  critical 0  sweeps {}  flips {}  \
                     samples moved {}",
                    done.sweeps,
                    done.flips,
                    lattice.moved.iter().filter(|m| **m).count(),
                );
            }
        });

        for arm in &arms {
            run.record(&[
                ("field", arm.field.to_string()),
                ("samples_per_axis", SAMPLES.to_string()),
                ("extractor", arm.extractor.to_string()),
                ("critical_before", arm.critical_before.to_string()),
                ("critical_after", arm.critical_after.to_string()),
                ("sweeps", arm.sweeps.to_string()),
                ("non_manifold_edges_before", arm.before.edges.to_string()),
                ("non_manifold_edges_after", arm.after.edges.to_string()),
                (
                    "non_manifold_vertices_before",
                    arm.before.vertices.to_string(),
                ),
                (
                    "non_manifold_vertices_after",
                    arm.after.vertices.to_string(),
                ),
                ("hausdorff_before", format!("{:.9e}", arm.before.hausdorff)),
                ("hausdorff_after", format!("{:.9e}", arm.after.hausdorff)),
                ("hausdorff_ratio", format!("{:.6}", arm.hausdorff_ratio())),
                ("critical_2d_before", arm.critical_2d_before.to_string()),
                ("critical_3d_before", arm.critical_3d_before.to_string()),
                (
                    "critical_after_two_sweeps",
                    arm.residual_after_two.to_string(),
                ),
                ("corners_moved", arm.flips.to_string()),
                ("stuck_cells", arm.stuck.to_string()),
                ("residual_distinct_bytes", arm.survivors.len().to_string()),
                (
                    "residual_bytes",
                    arm.survivors
                        .iter()
                        .map(|(byte, count)| format!("{byte:#04x}x{count}"))
                        .collect::<Vec<_>>()
                        .join(" "),
                ),
                ("residual_exhausted_cells", arm.exhausted.to_string()),
                ("residual_critical_2d", arm.residual_2d.to_string()),
                ("residual_critical_3d", arm.residual_3d.to_string()),
                (
                    "residual_nm_vertices_in_critical",
                    arm.residual_colocated.to_string(),
                ),
                ("samples_moved", arm.moved_samples.to_string()),
                ("samples_total", arm.total_samples.to_string()),
                (
                    "samples_moved_fraction",
                    format!("{:.9}", arm.moved_samples as f64 / arm.total_samples as f64),
                ),
                ("wrapper_faithful", arm.wrapper_faithful.to_string()),
                ("mesh_vertices_before", arm.before.mesh_vertices.to_string()),
                ("mesh_vertices_after", arm.after.mesh_vertices.to_string()),
                (
                    "mesh_triangles_before",
                    arm.before.mesh_triangles.to_string(),
                ),
                ("mesh_triangles_after", arm.after.mesh_triangles.to_string()),
                ("mesh_hash_before", arm.before.hash.to_string()),
                ("mesh_hash_after", arm.after.hash.to_string()),
                (
                    "accuracy_coverage_before",
                    arm.before.has_coverage.to_string(),
                ),
                (
                    "accuracy_coverage_after",
                    arm.after.has_coverage.to_string(),
                ),
            ]);
        }
    });
}
