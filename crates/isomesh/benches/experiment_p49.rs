//! **P-49 — the aperture: how big a thing fits between two chunk faces.**
//!
//! Ticket: R-045. Pre-registered before this harness existed.
//!
//! ```bash
//! cargo bench --bench experiment_p49
//! ```
//!
//! Writes `docs/experiments/p-49.csv`.
//!
//! # What is being computed
//!
//! `validate::sealing` answers *"is this sealed"*. A game asks *"can the player
//! get through"*, and that is a **bottleneck value** rather than a boolean: for a
//! pair of chunk faces, the maximum over air paths of the minimum
//! distance-to-solid along the path. A max-min path, i.e. the widest bottleneck.
//!
//! The algorithm is the standard descending-threshold union-find, which computes
//! exactly that:
//!
//! 1. Take every **air** sample — `value > 0.0`, the strict complement of
//!    `cube.rs::is_inside`'s `value < 0.0` with the surface itself excluded, so
//!    an aperture is a strictly positive clearance.
//! 2. Sort them **descending by `(field value, grid index)`**. That is a total
//!    order on distinct samples, so there is no PRNG, no atomics, no `HashMap`
//!    and no tie broken by address. Implemented as an ascending sort of
//!    `(!value.to_bits(), index)`: for a positive `f64` the IEEE bit pattern is
//!    monotone in the value, so the complement sorts descending, exactly, with
//!    integer comparisons rather than float ones.
//! 3. Activate them in that order, unioning each with its already-active
//!    6-neighbours. Each component carries a 6-bit mask of which grid faces it
//!    touches.
//! 4. The **first** moment a component's mask contains both faces of a pair, the
//!    current sample's value is that pair's aperture. First moment in descending
//!    order is the largest bottleneck, which is the definition.
//!
//! 6-connectivity, not 18 or 26: a diagonal step between two face-adjacent
//! solids passes through material, and a clearance a game gates movement on must
//! not claim it.
//!
//! The output is a 6×6 symmetric matrix plus a reachability mask — the
//! composable boundary summary, since neighbouring chunks combine matrices with
//! no global solve.
//!
//! # Why the slab is the clause that matters
//!
//! A `BoxExact` big enough to swallow the whole domain, with a
//! [`Capsule`](isomesh::brush::Capsule) of radius `r` subtracted along `x`.
//! Subtraction is `max(field, −shape)`, so inside the channel the box term is
//! deeply negative and the value is `−capsule` — which, `Capsule` being an exact
//! distance field, is *exactly* `r − ρ` at radius `ρ` from the axis. So:
//!
//! - the deepest air is the axis itself, at exactly `r`;
//! - the axis is a sample line, because `65` samples over `[−2, 2]` put a sample
//!   at `0.0` on every axis, so the exact value `r` is actually attained;
//! - the axis runs out of both `x` faces and reaches no other face, because the
//!   largest `r` tested is 8 cells `= 0.5` and the `±y`/`±z` walls are 2 away.
//!
//! Ground truth is therefore `r` with no discretisation slack, and exactly one
//! reachable pair. The instrument has nowhere to hide.
//!
//! # Two places this harness reports more than the registration asked for
//!
//! **The capped gyroid's six faces are connected by the cap's exterior, not by
//! the gyroid's channels.** `capped_gyroid` is `max(gyroid, sphere(6))` over
//! `[−7, 7]³`, so every point outside radius 6 has a positive value and is air —
//! including all six domain faces, which sit at `|x| = 7`. Worse, each pair of
//! faces *shares a corner sample*: `(7, 7, 7)` is on the `+x`, `+y` and `+z`
//! faces at once and has clearance `7√3 − 6 ≈ 6.12`. So clause two is satisfied
//! by the shell around the cap before the gyroid is consulted at all. It is
//! reported as registered, and the `gyroid_uncapped` row is added beside it:
//! the same grid, the same resolution, `Gyroid::canonical()` with no cap, where
//! the outer shell does not exist and the only way between two faces is through
//! the bicontinuous channel network the clause's reasoning actually invokes.
//!
//! **The cost clause is reported the expensive way.** *"The whole 6×6
//! computation"* is read as including sampling the grid, because that is what
//! *whole* means and because it is the reading that makes the clause harder to
//! pass — the right way to resolve an ambiguity in a cost claim. Marching Cubes'
//! own time includes its sampling too, so the comparison is like for like. The
//! marginal cost, given a chunk whose samples the mesher already has, is recorded
//! beside it as `aperture_ns_marginal`; that is the number a chunk pipeline would
//! actually budget, and it is much smaller.
//!
//! The early exit is part of the algorithm rather than a shortcut around it: once
//! all 15 pairs have an aperture, no later sample can change one, because each is
//! fixed at the first — and therefore highest — value that connected it.
//! `aperture_ns_full` reports the same computation with the exit disabled, which
//! is the worst case and the one to budget for.
//!
//! # Timing discipline
//!
//! Buffers allocated once and reused, one untimed warm-up, then the **median** of
//! [`REPS`] timed runs. Never a mean: the first run's page faults would drag it.
//! `loadavg` is recorded because a sibling agent is benching on this machine
//! concurrently and a cost ratio measured under load should say so.

mod common;

use std::time::Instant;

use isomesh::brush::{Brush, BrushStack, Capsule};
use isomesh::fields::{BoxExact, Gyroid, ReferenceField, capped_gyroid};
use isomesh::marching_cubes::MarchingCubes;
use isomesh::{MeshBuffer, RuntimeShape3, Sdf};

/// The registered resolution.
const SAMPLES: u32 = 65;

/// Channel radii, in cells.
const RADII: [u32; 3] = [2, 4, 8];

/// Timed repetitions. The median is reported.
const REPS: usize = 7;

/// Face order: `-x, +x, -y, +y, -z, +z`.
const FACE_NAMES: [&str; 6] = ["-x", "+x", "-y", "+y", "-z", "+z"];

/// All 15 pairs recorded.
const ALL_PAIRS: u16 = 0x7FFF;

// ─── the grid ───────────────────────────────────────────────────────────────

/// The grid a field is sampled on. Same convention as `common::grid`: `shape`
/// counts samples, so `n` samples span `n − 1` cells.
struct Grid {
    origin: [f64; 3],
    cell_size: f64,
    samples: u32,
}

impl Grid {
    fn shape(&self) -> RuntimeShape3 {
        RuntimeShape3::new([self.samples; 3]).expect("aperture grid fits u32")
    }

    fn point(&self, i: usize, j: usize, k: usize) -> [f64; 3] {
        [
            self.origin[0] + self.cell_size * i as f64,
            self.origin[1] + self.cell_size * j as f64,
            self.origin[2] + self.cell_size * k as f64,
        ]
    }
}

// ─── the pair bookkeeping ───────────────────────────────────────────────────

/// The 15 unordered face pairs, in a fixed order.
fn pair_list() -> [(usize, usize); 15] {
    let mut out = [(0usize, 0usize); 15];
    let mut k = 0;
    for i in 0..6usize {
        for j in (i + 1)..6usize {
            out[k] = (i, j);
            k += 1;
        }
    }
    out
}

/// For each 6-bit face mask, which of the 15 pairs it already contains.
///
/// Turns "has this union just connected a pair nobody has recorded yet" into one
/// table lookup and one mask test, rather than 15 comparisons per sample.
fn pair_table(pairs: &[(usize, usize); 15]) -> [u16; 64] {
    let mut table = [0u16; 64];
    for (mask, slot) in table.iter_mut().enumerate() {
        for (bit, (i, j)) in pairs.iter().enumerate() {
            if (mask >> i) & 1 == 1 && (mask >> j) & 1 == 1 {
                *slot |= 1 << bit;
            }
        }
    }
    table
}

// ─── the bottleneck engine ──────────────────────────────────────────────────

/// Monotone union-find over air samples, with every buffer allocated once.
struct Bottleneck {
    samples: u32,
    /// One field value per sample, `i + j·n + k·n²`.
    values: Vec<f64>,
    /// Air samples as `(!value.to_bits(), index)`, sorted ascending — which is
    /// descending by value, then ascending by index.
    order: Vec<(u64, u32)>,
    parent: Vec<u32>,
    size: Vec<u32>,
    /// Face mask, meaningful on a component root.
    mask: Vec<u8>,
    active: Vec<bool>,
    /// Aperture per pair, in world units. `None` is unreachable.
    aperture: [Option<f64>; 15],
    pairs: [(usize, usize); 15],
    table: [u16; 64],
    /// Air samples in the last solve, and how many were actually visited before
    /// the early exit fired.
    air: u64,
    visited: u64,
}

impl Bottleneck {
    fn new(samples: u32) -> Self {
        let n = samples as usize;
        let total = n * n * n;
        let pairs = pair_list();
        Self {
            samples,
            values: vec![0.0; total],
            order: Vec::with_capacity(total),
            parent: vec![0; total],
            size: vec![0; total],
            mask: vec![0; total],
            active: vec![false; total],
            aperture: [None; 15],
            table: pair_table(&pairs),
            pairs,
            air: 0,
            visited: 0,
        }
    }

    /// Fill `values` from the field. Separate from [`Bottleneck::solve`] so the
    /// marginal cost of the summary can be reported apart from the sampling a
    /// mesher pays anyway.
    fn sample_grid<F: Sdf<Scalar = f64>>(&mut self, field: &F, grid: &Grid) {
        let n = self.samples as usize;
        for k in 0..n {
            for j in 0..n {
                for i in 0..n {
                    self.values[i + j * n + k * n * n] = field.sample(grid.point(i, j, k));
                }
            }
        }
    }

    fn face_mask(&self, i: usize, j: usize, k: usize) -> u8 {
        let last = self.samples as usize - 1;
        u8::from(i == 0)
            | (u8::from(i == last) << 1)
            | (u8::from(j == 0) << 2)
            | (u8::from(j == last) << 3)
            | (u8::from(k == 0) << 4)
            | (u8::from(k == last) << 5)
    }

    fn find(&mut self, mut x: u32) -> u32 {
        while self.parent[x as usize] != x {
            let grand = self.parent[self.parent[x as usize] as usize];
            self.parent[x as usize] = grand;
            x = grand;
        }
        x
    }

    /// Union by size, so the tree stays shallow. Deterministic regardless: the
    /// processing order is a total order, so the same input gives the same tree.
    fn union(&mut self, a: u32, b: u32) {
        let (mut ra, mut rb) = (self.find(a), self.find(b));
        if ra == rb {
            return;
        }
        if self.size[ra as usize] < self.size[rb as usize] {
            core::mem::swap(&mut ra, &mut rb);
        }
        let merged = self.mask[ra as usize] | self.mask[rb as usize];
        self.parent[rb as usize] = ra;
        self.size[ra as usize] += self.size[rb as usize];
        self.mask[ra as usize] = merged;
    }

    /// The whole 6×6, from the values already in `values`.
    ///
    /// `early_exit` stops once all 15 pairs are known. Sound, not a shortcut:
    /// each aperture is fixed at the first and therefore highest value that
    /// connected its pair, so no later sample can revise one.
    fn solve(&mut self, early_exit: bool) {
        let n = self.samples as usize;
        let plane = n * n;

        self.aperture = [None; 15];
        self.order.clear();
        self.active.fill(false);

        for (index, value) in self.values.iter().enumerate() {
            if *value > 0.0 {
                self.order.push((!value.to_bits(), index as u32));
            }
        }
        self.order.sort_unstable();
        self.air = self.order.len() as u64;
        self.visited = 0;

        let mut recorded = 0u16;
        for step in 0..self.order.len() {
            let (key, index) = self.order[step];
            let value = f64::from_bits(!key);
            let flat = index as usize;
            let i = flat % n;
            let j = (flat / n) % n;
            let k = flat / plane;

            self.active[flat] = true;
            self.parent[flat] = index;
            self.size[flat] = 1;
            self.mask[flat] = self.face_mask(i, j, k);
            self.visited += 1;

            // The six face-adjacent neighbours, only where they exist.
            if i > 0 && self.active[flat - 1] {
                self.union(index, index - 1);
            }
            if i + 1 < n && self.active[flat + 1] {
                self.union(index, index + 1);
            }
            if j > 0 && self.active[flat - n] {
                self.union(index, index - n as u32);
            }
            if j + 1 < n && self.active[flat + n] {
                self.union(index, index + n as u32);
            }
            if k > 0 && self.active[flat - plane] {
                self.union(index, index - plane as u32);
            }
            if k + 1 < n && self.active[flat + plane] {
                self.union(index, index + plane as u32);
            }

            let root = self.find(index);
            let fresh = self.table[self.mask[root as usize] as usize] & !recorded;
            if fresh != 0 {
                for (bit, slot) in self.aperture.iter_mut().enumerate() {
                    if fresh & (1 << bit) != 0 {
                        *slot = Some(value);
                    }
                }
                recorded |= fresh;
                if early_exit && recorded == ALL_PAIRS {
                    break;
                }
            }
        }
    }

    /// Reachable pairs.
    fn reachable(&self) -> u64 {
        self.aperture.iter().filter(|a| a.is_some()).count() as u64
    }

    /// Solve twice from the same values and compare all 15 apertures bit for
    /// bit, leaving the second answer in place.
    ///
    /// The hypothesis claims determinism outright — "a total order and therefore
    /// deterministic with no PRNG, no atomics and no HashMap" — so it is checked
    /// rather than asserted from the absence of those three things. `Option<f64>`
    /// compares by value here, which is what "bit for bit" needs: the apertures
    /// are values copied straight out of the sample array, never arithmetic
    /// results, so an exact comparison is the right one and cannot be flaky.
    fn is_deterministic(&mut self, early_exit: bool) -> bool {
        self.solve(early_exit);
        let first = self.aperture;
        let (air, visited) = (self.air, self.visited);
        self.solve(early_exit);
        first == self.aperture && air == self.air && visited == self.visited
    }

    /// The aperture of one pair, in cells.
    fn cells(&self, bit: usize, cell_size: f64) -> Option<f64> {
        self.aperture[bit].map(|v| v / cell_size)
    }

    /// The pair bit for an unordered face pair.
    fn bit_of(&self, a: usize, b: usize) -> usize {
        let (lo, hi) = if a < b { (a, b) } else { (b, a) };
        self.pairs
            .iter()
            .position(|p| *p == (lo, hi))
            .expect("every unordered pair of six faces is in the table")
    }

    /// The 6×6, as `-x/+x=2.0000 -x/-y=unreachable ...`. Symmetric, so only the
    /// 15 upper-triangle entries are printed; the diagonal is not a pair.
    fn matrix(&self, cell_size: f64) -> String {
        self.pairs
            .iter()
            .enumerate()
            .map(|(bit, (i, j))| match self.cells(bit, cell_size) {
                Some(v) => format!("{}/{}={v:.4}", FACE_NAMES[*i], FACE_NAMES[*j]),
                None => format!("{}/{}=unreachable", FACE_NAMES[*i], FACE_NAMES[*j]),
            })
            .collect::<Vec<_>>()
            .join(" ")
    }
}

// ─── timing ─────────────────────────────────────────────────────────────────

/// Median of a list of timings. Never a mean.
fn median(mut runs: Vec<f64>) -> f64 {
    runs.sort_by(f64::total_cmp);
    runs[runs.len() / 2]
}

/// One-minute load average, or `unknown`.
///
/// Recorded because a sibling is benching concurrently, and a cost ratio
/// measured under load should carry that on the artefact rather than in
/// somebody's memory.
fn loadavg() -> String {
    std::fs::read_to_string("/proc/loadavg")
        .ok()
        .and_then(|s| s.split_whitespace().next().map(str::to_string))
        .unwrap_or_else(|| String::from("unknown"))
}

/// Median nanoseconds for one Marching Cubes extraction of this grid, and the
/// triangle count it produced.
fn extract_ns<F: Sdf<Scalar = f64>>(field: &F, grid: &Grid) -> (f64, usize) {
    let shape = grid.shape();
    let mut mc = MarchingCubes::<f64>::new();
    let mut out = MeshBuffer::<f64>::new();

    // One untimed run so every buffer is allocated and every page is resident.
    out.reset();
    mc.extract(field, &shape, grid.origin, grid.cell_size, &mut out)
        .expect("marching cubes extraction");
    let triangles = out.triangle_count();

    let mut runs = Vec::with_capacity(REPS);
    for _ in 0..REPS {
        let t = Instant::now();
        out.reset();
        mc.extract(field, &shape, grid.origin, grid.cell_size, &mut out)
            .expect("marching cubes extraction");
        runs.push(t.elapsed().as_secs_f64() * 1e9);
    }
    (median(runs), triangles)
}

/// The three timings the cost clause needs, all on buffers allocated once.
struct Cost {
    /// Sampling plus solve, with the early exit. The registered number.
    whole: f64,
    /// Solve only, given values the mesher already has.
    marginal: f64,
    /// Sampling plus solve with the early exit disabled: the worst case.
    full: f64,
}

fn aperture_cost<F: Sdf<Scalar = f64>>(engine: &mut Bottleneck, field: &F, grid: &Grid) -> Cost {
    // One untimed warm-up so every buffer is resident.
    engine.sample_grid(field, grid);
    engine.solve(true);

    let mut whole = Vec::with_capacity(REPS);
    for _ in 0..REPS {
        let t = Instant::now();
        engine.sample_grid(field, grid);
        engine.solve(true);
        whole.push(t.elapsed().as_secs_f64() * 1e9);
    }

    let mut marginal = Vec::with_capacity(REPS);
    for _ in 0..REPS {
        let t = Instant::now();
        engine.solve(true);
        marginal.push(t.elapsed().as_secs_f64() * 1e9);
    }

    let mut full = Vec::with_capacity(REPS);
    for _ in 0..REPS {
        let t = Instant::now();
        engine.sample_grid(field, grid);
        engine.solve(false);
        full.push(t.elapsed().as_secs_f64() * 1e9);
    }

    // Leave the engine holding the early-exit answer, which is the one reported.
    engine.solve(true);

    Cost {
        whole: median(whole),
        marginal: median(marginal),
        full: median(full),
    }
}

// ─── one fixture ────────────────────────────────────────────────────────────

/// Which number `aperture_reported_cells` means for this fixture.
#[derive(Clone, Copy)]
enum Scalar {
    /// The aperture of one named pair, against a known radius.
    Pair(usize, usize, u32),
    /// The worst of the 15, where there is no single channel to name.
    Worst,
}

/// The set of pairs a fixture is known to connect, and what its scalar means.
struct Truth {
    expected: u16,
    scalar: Scalar,
}

/// Air samples a straight cylinder of radius `radius` cells must produce on an
/// `n³` grid whose axis is a sample line.
///
/// A sample is air exactly when `−capsule > 0`, i.e. strictly inside the
/// cylinder, i.e. `j² + k² < radius²` in cell units — and the cylinder spans
/// every one of the `n` planes along `x`. Integer arithmetic throughout, so this
/// is a closed form and not a second approximation of the same thing.
fn cylinder_air(samples: u32, radius: u32) -> u64 {
    let r2 = i64::from(radius) * i64::from(radius);
    let span = i64::from(radius);
    let mut per_plane = 0u64;
    for j in -span..=span {
        for k in -span..=span {
            if j * j + k * k < r2 {
                per_plane += 1;
            }
        }
    }
    per_plane * u64::from(samples)
}

/// Everything one fixture produced.
struct Row {
    fixture: String,
    channel_radius: Option<u32>,
    reported_cells: Option<f64>,
    error_cells: Option<f64>,
    reachable: u64,
    expected_reachable: u64,
    false_reachable: u64,
    missed_reachable: u64,
    cost: Cost,
    extract: f64,
    triangles: usize,
    air: u64,
    visited: u64,
    matrix: String,
    min_cells: Option<f64>,
    max_cells: Option<f64>,
    analytic_gradient: bool,
    /// Air samples the geometry says there must be, where a closed form exists.
    ///
    /// On the drilled slab it does: a sample is air exactly when it is strictly
    /// inside the cylinder, so the count is `samples_per_axis` times the number
    /// of lattice points with `j² + k² < r²`. This is the instrument checking its
    /// own input rather than only its output — an aperture of exactly `r` read off
    /// the wrong air set would still be wrong.
    air_expected: Option<u64>,
    /// Two solves from the same values agreed on all 15 apertures.
    deterministic: bool,
}

impl Row {
    fn cost_ratio(&self) -> f64 {
        self.cost.whole / self.extract
    }
}

/// Measure one fixture end to end.
fn measure<F: Sdf<Scalar = f64>>(
    fixture: String,
    field: &F,
    grid: &Grid,
    truth: &Truth,
    engine: &mut Bottleneck,
    analytic_gradient: bool,
) -> Row {
    let cost = aperture_cost(engine, field, grid);
    let (extract, triangles) = extract_ns(field, grid);
    // Checked before the reported answer is read off, and it leaves the same
    // answer behind, so this costs two extra solves and nothing else.
    let deterministic = engine.is_deterministic(true);

    let mut reported = 0u16;
    for (bit, slot) in engine.aperture.iter().enumerate() {
        if slot.is_some() {
            reported |= 1 << bit;
        }
    }
    // The two directions of error, kept apart. A falsely reachable pair is the
    // unsound one: a game would gate movement on it.
    let false_reachable = u64::from((reported & !truth.expected).count_ones());
    let missed_reachable = u64::from((truth.expected & !reported).count_ones());

    let (reported_cells, error_cells) = match truth.scalar {
        Scalar::Pair(a, b, radius) => {
            let bit = engine.bit_of(a, b);
            let v = engine.cells(bit, grid.cell_size);
            (v, v.map(|v| (v - f64::from(radius)).abs()))
        }
        Scalar::Worst => {
            let worst = engine
                .aperture
                .iter()
                .flatten()
                .copied()
                .min_by(f64::total_cmp)
                .map(|v| v / grid.cell_size);
            (worst, None)
        }
    };

    let present: Vec<f64> = engine
        .aperture
        .iter()
        .flatten()
        .map(|v| v / grid.cell_size)
        .collect();

    Row {
        fixture,
        channel_radius: match truth.scalar {
            Scalar::Pair(_, _, r) => Some(r),
            Scalar::Worst => None,
        },
        reported_cells,
        error_cells,
        reachable: engine.reachable(),
        expected_reachable: u64::from(truth.expected.count_ones()),
        false_reachable,
        missed_reachable,
        cost,
        extract,
        triangles,
        air: engine.air,
        visited: engine.visited,
        matrix: engine.matrix(grid.cell_size),
        min_cells: present.iter().copied().min_by(f64::total_cmp),
        max_cells: present.iter().copied().max_by(f64::total_cmp),
        analytic_gradient,
        air_expected: match truth.scalar {
            Scalar::Pair(_, _, r) => Some(cylinder_air(grid.samples, r)),
            Scalar::Worst => None,
        },
        deterministic,
    }
}

/// A number, or `n/a` where the fixture does not define one.
fn cell_or_na(value: Option<f64>) -> String {
    value.map_or_else(|| String::from("n/a"), |v| format!("{v:.6}"))
}

/// A radius, or `n/a` on a fixture with no drilled channel.
fn radius_or_na(value: Option<u32>) -> String {
    value.map_or_else(|| String::from("n/a"), |r| r.to_string())
}

// ─── main ───────────────────────────────────────────────────────────────────

fn main() {
    let prereg = isomesh::experiment!("P-49");

    common::experiment::run(prereg, |run| {
        let load = loadavg();
        println!("loadavg (1 min) at start: {load}\n");

        let mut engine = Bottleneck::new(SAMPLES);
        let mut rows: Vec<Row> = Vec::new();

        // ── clause one: exact ground truth on a drilled slab ────────────────
        //
        // The box swallows the whole domain, so every sample is solid unless the
        // capsule carves it. The capsule runs well past both x walls, so the
        // channel is a straight cylinder rather than a capped one anywhere
        // inside the grid.
        let (lo, hi) = BoxExact::<f64>::canonical().domain();
        let cell_size = (hi[0] - lo[0]) / f64::from(SAMPLES - 1);
        let slab_grid = Grid {
            origin: lo,
            cell_size,
            samples: SAMPLES,
        };
        for radius in RADII {
            let brushes = [Brush::subtract(Capsule {
                a: [-4.0, 0.0, 0.0],
                b: [4.0, 0.0, 0.0],
                radius: cell_size * f64::from(radius),
            })];
            let field = BrushStack {
                base: BoxExact::<f64> {
                    center: [0.0; 3],
                    half_extents: [4.0; 3],
                },
                brushes: &brushes,
            };
            // Exactly one pair is truly reachable: the channel's two ends.
            let truth = Truth {
                expected: 1 << 0,
                scalar: Scalar::Pair(0, 1, radius),
            };
            rows.push(measure(
                format!("slab_r{radius}"),
                &field,
                &slab_grid,
                &truth,
                &mut engine,
                false,
            ));
        }

        // ── clause two: a real field ────────────────────────────────────────
        let capped = capped_gyroid::<f64>();
        let (glo, ghi) = capped.domain();
        let gyroid_grid = Grid {
            origin: glo,
            cell_size: (ghi[0] - glo[0]) / f64::from(SAMPLES - 1),
            samples: SAMPLES,
        };
        let all = Truth {
            expected: ALL_PAIRS,
            scalar: Scalar::Worst,
        };
        rows.push(measure(
            String::from("gyroid_capped"),
            &capped,
            &gyroid_grid,
            &all,
            &mut engine,
            true,
        ));

        // Not registered, and the measurement clause two's *reasoning* describes:
        // no cap, so no exterior shell, so the only route between two faces is
        // the bicontinuous channel network itself.
        rows.push(measure(
            String::from("gyroid_uncapped"),
            &Gyroid::<f64>::canonical(),
            &gyroid_grid,
            &all,
            &mut engine,
            true,
        ));

        for row in &rows {
            println!(
                "{:>16}  aperture {:>10} cells (truth {:>4}, err {:>10})  \
                 reachable {:>2}/{:<2} false {} missed {}  \
                 aperture {:>10.0} ns / extract {:>10.0} ns = {:.4}×",
                row.fixture,
                cell_or_na(row.reported_cells),
                radius_or_na(row.channel_radius),
                cell_or_na(row.error_cells),
                row.reachable,
                row.expected_reachable,
                row.false_reachable,
                row.missed_reachable,
                row.cost.whole,
                row.extract,
                row.cost_ratio(),
            );
            println!(
                "{:>16}  air {} samples, {} visited before exit; marginal {:.0} ns \
                 ({:.4}×), full-sweep {:.0} ns ({:.4}×); {} triangles; analytic \
                 gradient {}",
                "",
                row.air,
                row.visited,
                row.cost.marginal,
                row.cost.marginal / row.extract,
                row.cost.full,
                row.cost.full / row.extract,
                row.triangles,
                row.analytic_gradient,
            );
            println!(
                "{:>16}  air census {}closed form {}; deterministic over two solves: {}",
                "",
                match row.air_expected {
                    Some(e) if e == row.air => String::from("MATCHES "),
                    Some(_) => String::from("DISAGREES with "),
                    None => String::from("no "),
                },
                row.air_expected
                    .map_or_else(|| String::from("(none for this field)"), |e| e.to_string()),
                row.deterministic,
            );
            println!("{:>16}  {}", "", row.matrix);
        }

        for row in &rows {
            run.record(&[
                ("fixture", row.fixture.clone()),
                ("samples_per_axis", SAMPLES.to_string()),
                ("channel_radius_cells", radius_or_na(row.channel_radius)),
                ("aperture_reported_cells", cell_or_na(row.reported_cells)),
                ("aperture_error_cells", cell_or_na(row.error_cells)),
                ("reachable_pairs", row.reachable.to_string()),
                (
                    "expected_reachable_pairs",
                    row.expected_reachable.to_string(),
                ),
                ("false_reachable_pairs", row.false_reachable.to_string()),
                ("aperture_ns", format!("{:.0}", row.cost.whole)),
                ("extract_ns", format!("{:.0}", row.extract)),
                ("cost_ratio", format!("{:.6}", row.cost_ratio())),
                ("missed_reachable_pairs", row.missed_reachable.to_string()),
                ("aperture_ns_marginal", format!("{:.0}", row.cost.marginal)),
                (
                    "cost_ratio_marginal",
                    format!("{:.6}", row.cost.marginal / row.extract),
                ),
                ("aperture_ns_full", format!("{:.0}", row.cost.full)),
                (
                    "cost_ratio_full",
                    format!("{:.6}", row.cost.full / row.extract),
                ),
                ("air_samples", row.air.to_string()),
                ("air_samples_visited", row.visited.to_string()),
                (
                    "air_samples_expected",
                    row.air_expected
                        .map_or_else(|| String::from("n/a"), |e| e.to_string()),
                ),
                (
                    "air_census_agrees",
                    row.air_expected
                        .map_or_else(|| String::from("n/a"), |e| (e == row.air).to_string()),
                ),
                ("deterministic", row.deterministic.to_string()),
                ("aperture_min_cells", cell_or_na(row.min_cells)),
                ("aperture_max_cells", cell_or_na(row.max_cells)),
                ("mc_triangles", row.triangles.to_string()),
                ("analytic_gradient", row.analytic_gradient.to_string()),
                ("connectivity", String::from("6")),
                ("reps", REPS.to_string()),
                ("loadavg_1min", load.clone()),
                ("aperture_matrix_cells", row.matrix.clone()),
            ]);
        }
    });
}
