//! E-313 — the chunk that knows it is losing detail.
//!
//! ```bash
//! cd bevy_isomesh && cargo run --example game_lod_oracle --release
//! ```
//!
//! **Always `--release`.** Startup meshes 16 chunks at 4 resolutions each,
//! censuses every one of those 64 with an exhaustive tangency scan, and runs one
//! more exhaustive 65³ scan to hold the predicate against
//! `docs/experiments/p-51.csv`. One second on 17 threads here; minutes in debug.
//!
//! Game framing of **M-355 / P-51 clause C3**.
//!
//! `1` uniform, `2` range LOD, `3` range matched to the gate, `4` oracle, `5`
//! back to the scripted beats. Under `ISOMESH_CAPTURE` it needs no keyboard: the
//! four beats come off the captured-frame counter.
//!
//! ```bash
//! ISOMESH_WINDOW=1280x720 ./scripts/record_gif.sh game_lod_oracle docs/gifs/e313.gif
//! ```
//!
//! # The idea
//!
//! A game picks chunk resolution from camera range and hopes. But every SDF
//! sample asserts a **sphere** of radius `|d(p)|` that the surface has to touch,
//! and counting the spheres a mesh never reaches needs only the field and the
//! mesh — no reference model, no ground truth, nothing baked from a finer
//! version of the same chunk. So a chunk can score *itself*, and this compares
//! what it says against what range alone decides.
//!
//! # What is on screen
//!
//! A 4×4 grid of one-unit chunks, each meshed independently at its own
//! resolution, over a scene with three feature scales: **plain rolling ground**
//! (four spheres of radius 9–11), a **smooth rampart** along the far edge
//! (radius 1.05), two **boulders** (radius 0.32) and nine **thin spires**
//! (radius 0.045, which is 0.72 cells across at the coarsest level and 5.76 at
//! the finest).
//!
//! Chunks are coloured by their own untouched fraction, with the **hue breaking
//! on the gate** — cool below, warm above — so "this chunk is losing detail" is
//! literally the colour of the ground under your feet.
//!
//! # The four beats
//!
//! | beat | allocation | triangles | chunks over gate |
//! |---|---|---:|---:|
//! | 1 | uniform 9³ | 2,738 | 8 of 16 |
//! | 2 | range LOD, 10 px/cell | 24,076 | 9 of 16 |
//! | 3 | oracle | **81,276** | **0** |
//! | 4 | A/B: range LOD tightened to 4.1 px/cell | 147,504 | 0 |
//!
//! **Beat 4 is the payoff and it is a matched comparison.** Beat 2's range
//! policy is the one a game ships and it leaves nine chunks losing detail, so
//! comparing its triangle count against the oracle's would be comparing two
//! different pictures. So the policy's own pixel budget is *searched downward*
//! until it clears the gate on every chunk — 4.1 px/cell — and that costs
//! **147,504 triangles against the oracle's 81,276, 1.81× for the same detail
//! bar.** The two sides of the A/B look the same; only the counter moves.
//!
//! And the reason is structural rather than a badly chosen constant. This scene
//! spans 4.25 to 7.27 units of range, a factor of **1.7**, while the ladder spans
//! a factor of **8** in cell size. A monotone function of range cannot express
//! "this chunk needs four times its neighbour's resolution at the same
//! distance", so tightening the budget until the worst chunk is served drags
//! almost every other chunk to the top of the ladder with it.
//!
//! # The predicate, and what licenses it
//!
//! P-51's, verbatim: a sample counts as **untouched** when the nearest mesh
//! vertex misses tangency to its sphere by more than [`THRESHOLD_CELLS`] cells.
//! The search is exhaustive over every vertex of the chunk with P-51's own
//! window pruning — there is no cutoff to blame a count on.
//!
//! Before the window opens, the same predicate is run **whole-domain, no scope,
//! 65³** on `thin_plate` with Marching Cubes and held against the committed
//! artefact. It reproduces exactly: 2,046 vertices, 196,956 untouched,
//! `717.1816` per 1,000 and a worst miss of `1.000000` cells, all four matching
//! `p-51.csv` to the digit it prints. That is what makes the per-chunk numbers
//! more than a plausible-looking ratio.
//!
//! # Two things had to be true of the field, and both were learned the hard way
//!
//! **It has to be an exact distance.** A heightfield `y − H(x, z)` is not one,
//! and its error is a fixed number of *world units* — so dividing by the cell
//! size turns it into a *growing* number of cells. The first draft's smooth-only
//! chunk reported a worst miss of `0.0173, 0.0259, 0.0435, 0.0540` cells across
//! four levels: exactly `1/h`, the field's own error rising as the mesh got
//! better. Every primitive here is exact, and the census counts **exterior
//! samples only**, because distance to a union of sets is the min of the
//! distances on the outside and an understatement within (M-246).
//!
//! **It has to be smooth.** P-51 measured Marching Cubes missing `box_exact`'s
//! corner samples by `√2` cells at every resolution, because it puts no vertex
//! on a box edge at all. The first draft's spires were thin boxes and their
//! chunks plateaued near 30% untouched from the coarsest level to the finest —
//! an oracle reading that says "give me more" forever. Spheres and capsules
//! resolve; convex sharp edges never do.
//!
//! With both fixed the meter behaves. Measured, untouched fraction at
//! 9³/17³/33³/65³:
//!
//! | chunk | 9³ | 17³ | 33³ | 65³ |
//! |---|---:|---:|---:|---:|
//! | plain ground `(0,0)` | 0.00% | 0.00% | 0.00% | 0.00% |
//! | rampart `(0,3)` | 20.34% | 25.75% | 10.97% | 5.27% |
//! | boulder `(0,2)` | 61.84% | 41.83% | 26.09% | 12.31% |
//! | spires `(0,1)` | 46.36% | 45.91% | 29.56% | **14.40%** |
//!
//! Plain ground is satisfied at the cheapest level there is and stays satisfied;
//! everything with curvature falls roughly linearly in the cell size. That is
//! the whole signal, and it is three-valued: the oracle answers **9³ for eight
//! chunks, 33³ for three and 65³ for five**.
//!
//! # It is not monotone, and the rule is first-crossing rather than minimum
//!
//! Chunk `(3,3)` reads **11.86% at 9³ and 20.15% at 17³** before settling to
//! 3.98% at 65³, so a chunk can clear the gate, fail it one level up, and clear
//! it again. Six of the sixteen do something like that, and **every one of them
//! does it on the 9³ → 17³ step**, which names the cause: the denominator. At 9³
//! a chunk owns only **59 to 110** in-scope spheres against 33,000 to 53,000 at
//! 65³, so a coarse-level fraction is quantised at about 1.7% and turns on which
//! handful of samples happen to land in scope.
//!
//! [`oracle_level`] therefore climbs from the cheapest rung and stops at the
//! **first** level that clears — which is what a runtime would do — and
//! [`audit`] prints every chunk whose ladder rises by a percentage point or more
//! rather than smoothing it away. On this scene the early stops are not wrong:
//! `(3,3)` reads 3.98% at the top of its own ladder.
//!
//! # The gate is 15%, and it is pinned to the measurement
//!
//! See [`GATE`]. There is an irreducible residue — P-51 scores tangency against
//! *vertices*, and a triangle's interior is not a vertex — so the bar is set at
//! the smallest whole percent the top of the ladder clears on **every** chunk,
//! which is 15% against a measured worst of 14.47%. [`audit`] asserts that both
//! ways: 15% clears, and 14% does not.
//!
//! # What this costs, honestly
//!
//! The census is not free and the HUD says so. The heaviest single
//! `(chunk, level)` census by work is chunk `(1,2)` at 65³: **61,551 in-scope
//! spheres against 6,875 vertices**, and the slowest rung on any run lands
//! around **0.8 s**. So this is a **bake or a background budget, not a per-frame
//! query**. The counts are the figures worth quoting — they are exact and
//! reproduce, while which rung wins the stopwatch moves between runs. What the
//! finding buys is that the bake needs no reference model at all: the field the
//! chunk already has is the ground truth.

mod common;

use std::thread;
use std::time::Instant;

use bevy::prelude::*;
use bevy_isomesh::MeshBuilder;
use common::{Capture, CommonPlugin, DemoMesh, DemoStats, OrbitCamera, ViewFlags};
use isomesh::fields::{ReferenceField, ThinPlate};
use isomesh::marching_cubes::MarchingCubes;
use isomesh::{MeshBuffer, MeshSink, RuntimeShape3, Sdf};

// ─── the registered predicate ───────────────────────────────────────────────

/// The per-sample gate, in cells. Registered in P-51: a sample's sphere counts
/// as **untouched** when the nearest mesh vertex misses tangency by more than
/// this.
const THRESHOLD_CELLS: f64 = 0.05;

/// The per-chunk gate: the fraction of a chunk's in-scope spheres allowed to go
/// untouched before the chunk asks for more resolution.
///
/// # Why it is 15% and not 1%
///
/// There is an irreducible residue, and it is not error. P-51 measures tangency
/// against **vertices**; a triangle's interior is not a vertex, so a sample
/// sitting a fraction of a cell off the surface has a sphere that the mesh
/// passes through and no vertex lands on. That residue scales with how much
/// curvature the cell has to follow, and the ladder here cannot drive it to
/// zero: measured, a chunk holding nothing but the `R ≈ 9` ground reads
/// **0.00% at every level**, while the `R = 1.05` rampart still reads **6.93%
/// at 65³** and the tightest chunk in the scene reads **14.47%** there.
///
/// So the gate is set at **the smallest whole percent the top of the ladder
/// clears on every chunk**. Tighter and the oracle's answer on at least one
/// chunk would be *the ladder is not enough*, which is a different conversation
/// from *this chunk wants the next level*. [`audit`] asserts both halves of
/// that: `GATE` clears everywhere and `GATE − 1%` does not, so the constant is
/// pinned to the measurement rather than dialled to taste.
const GATE: f64 = 0.15;

/// How much a chunk's fraction may rise between levels before [`audit`] calls
/// its ladder non-monotone: one percentage point.
const MONOTONE_SLACK: f64 = 0.01;

/// `docs/experiments/p-51.csv`, compiled in.
const LEDGER: &str = include_str!("../../docs/experiments/p-51.csv");

/// The ledger row this example reproduces, and its resolution.
const LEDGER_FIELD: &str = "thin_plate";
const LEDGER_EXTRACTOR: &str = "marching_cubes";
const LEDGER_SAMPLES: u32 = 65;

// ─── the scene ──────────────────────────────────────────────────────────────

/// Chunks along `x` and `z`. One layer in `y`, so the grid the viewer sees is
/// the grid the oracle scores.
const CHUNKS_X: usize = 4;
const CHUNKS_Z: usize = 4;
const CHUNKS: usize = CHUNKS_X * CHUNKS_Z;

/// Chunk edge, world units. Cubic, so `y` needs no second index.
const CHUNK: f64 = 1.0;

/// Samples per axis, coarsest first. `n` samples span `n − 1` cells, so the
/// cell size is `CHUNK / (n − 1)`: `0.125`, `0.0625`, `0.03125`, `0.015625`.
/// Powers of two, so a level change halves the cell exactly.
const LADDER: [u32; 4] = [9, 17, 33, 65];

/// The rolling ground: four **very large spheres**, unioned.
///
/// A heightfield is the obvious way to write terrain and it is the wrong field
/// for this instrument. `y − H(x, z)` is not a distance — it overstates by up to
/// `sqrt(1 + |∇H|²)` — and dividing by that bound only turns the overstatement
/// into an understatement of the same order. Either way the error is a fixed
/// number of *world units*, so dividing it by the cell size makes it a **growing
/// number of cells**: measured, the first draft of this example reported a worst
/// miss of `0.0173, 0.0259, 0.0435, 0.0540` cells across four levels on a chunk
/// holding nothing but smooth ground — exactly `1/h`, the field's own error
/// rising as the mesh got better. A detail meter that complains harder the more
/// resolution you give it is not a detail meter.
///
/// So every primitive here is an **exact** distance, and every one of them is
/// **smooth**. Sharp convex edges are the other trap: P-51 measured Marching
/// Cubes leaving `box_exact`'s corner samples missed by `√2` cells because it
/// puts no vertex on a box edge at all, at *any* resolution — the first draft's
/// blades were thin boxes and their chunks plateaued at 30% untouched from the
/// coarsest level to the finest. Spheres and capsules resolve; box edges never
/// do.
///
/// Each entry is `[x, y, z, radius]`, and the crest of the hill it makes is at
/// `y + radius`.
const HILLS: [[f64; 4]; 4] = [
    [0.40, -8.48, 0.70, 9.00],
    [3.60, -8.50, 0.90, 9.00],
    [1.20, -10.56, 3.30, 11.00],
    [3.40, -9.52, 3.00, 10.00],
];

/// The rampart: a broad smooth ridge along the **far** edge of the grid, a
/// capsule along `x` running out past both ends of the domain.
///
/// **This is the thing that does not need resolution.** Radius `1.05`, so its
/// curvature is `0.95` per unit and the middle of the ladder already tracks it —
/// it is what the oracle answers *between* the two extremes with, so the
/// allocation is three-valued rather than cheap-or-expensive.
///
/// **Far, not near, and that was measured.** Its first home was the near chunk
/// row, where at that range it faced the camera nearly square-on and filled the
/// bottom third of the frame with one orange dome — the near row's own heat
/// tiles were invisible behind the object they were scoring. On the far edge it
/// reads as a ridge on the horizon, and it leaves the near row as **plain
/// rolling ground**, which is the sharpest possible statement of the waste: the
/// flattest, emptiest chunks in the scene are the closest ones, so a range
/// policy spends its most on exactly the chunks whose field asks for the least.
///
/// `z = 3.85` and `y = −0.34` put the crest at `0.71` and the line where it
/// leaves the ground at `z = 3.10`, so it belongs to the far chunk row and
/// nothing else.
const RAMPART: ([f64; 3], [f64; 3], f64) = ([-0.60, -0.34, 3.85], [4.60, -0.34, 3.85], 1.05);

/// Two mid-scale boulders, one at each end of the middle row.
const BOULDERS: [[f64; 4]; 2] = [[0.62, 0.50, 2.38, 0.32], [2.55, 0.52, 2.42, 0.32]];

/// The nine thin spires, in three clusters: two at the left and right edges of
/// the frame, where a range policy starves them, and one mid-field. Each is a
/// vertical capsule from `SPIRE_BASE` to its own top, so `[x, z, top]`.
const SPIRES: [[f64; 3]; 9] = [
    [0.28, 1.42, 0.86],
    [0.52, 1.52, 0.93],
    [0.78, 1.46, 0.80],
    [3.24, 1.47, 0.88],
    [3.50, 1.41, 0.78],
    [3.76, 1.51, 0.91],
    [1.26, 2.48, 0.84],
    [1.52, 2.57, 0.94],
    [1.78, 2.45, 0.82],
];

/// Spire radius and the depth its capsule starts at.
///
/// `0.045` is `0.72` cells across at the coarsest level — under one cell, so
/// Marching Cubes finds no sign change and the spire is simply **not there** —
/// and `5.76` cells across at the finest, where it is a clean pole.
const SPIRE_R: f64 = 0.045;
const SPIRE_BASE: f64 = 0.15;

/// The whole world as one signed distance field.
///
/// A union by `min` of exact primitives. **Distance to a union of sets is the
/// min of the distances**, so on the exterior — where every sphere this example
/// counts lives — `min` is exact, whatever the primitives do to each other. It
/// is only *inside* the solid that `min` understates (M-246), which is why the
/// census counts exterior samples and nothing else.
#[derive(Clone, Copy)]
struct Scene;

impl Sdf for Scene {
    type Scalar = f64;

    fn sample(&self, p: [f64; 3]) -> f64 {
        let mut d = f64::INFINITY;
        for hill in &HILLS {
            d = d.min(ball(p, [hill[0], hill[1], hill[2]], hill[3]));
        }
        for rock in &BOULDERS {
            d = d.min(ball(p, [rock[0], rock[1], rock[2]], rock[3]));
        }
        d = d.min(capsule(p, RAMPART.0, RAMPART.1, RAMPART.2));
        for spire in &SPIRES {
            d = d.min(capsule(
                p,
                [spire[0], SPIRE_BASE, spire[1]],
                [spire[0], spire[2], spire[1]],
                SPIRE_R,
            ));
        }
        d
    }
}

/// Exact distance to a sphere.
fn ball(p: [f64; 3], at: [f64; 3], radius: f64) -> f64 {
    let (dx, dy, dz) = (p[0] - at[0], p[1] - at[1], p[2] - at[2]);
    (dx * dx + dy * dy + dz * dz).sqrt() - radius
}

/// Exact distance to a capsule: the segment `a → b` swept by `radius`.
fn capsule(p: [f64; 3], a: [f64; 3], b: [f64; 3], radius: f64) -> f64 {
    let ba = [b[0] - a[0], b[1] - a[1], b[2] - a[2]];
    let pa = [p[0] - a[0], p[1] - a[1], p[2] - a[2]];
    let along = ba[0] * ba[0] + ba[1] * ba[1] + ba[2] * ba[2];
    let t = ((pa[0] * ba[0] + pa[1] * ba[1] + pa[2] * ba[2]) / along).clamp(0.0, 1.0);
    let (dx, dy, dz) = (pa[0] - ba[0] * t, pa[1] - ba[1] * t, pa[2] - ba[2] * t);
    (dx * dx + dy * dy + dz * dz).sqrt() - radius
}

// ─── one rung of one chunk's ladder ─────────────────────────────────────────

/// One `(chunk, level)` measurement: the extraction, and the oracle's verdict.
struct Rung {
    n: u32,
    h: f64,
    samples: usize,
    in_scope: usize,
    untouched: usize,
    worst_cells: f64,
    buffer: MeshBuffer<f64>,
    extract_ms: f64,
    scan_ms: f64,
}

impl Rung {
    /// Untouched spheres as a fraction of the ones this chunk owns.
    fn frac(&self) -> f64 {
        self.untouched as f64 / self.in_scope as f64
    }

    fn triangles(&self) -> usize {
        self.buffer.indices.len() / 3
    }

    fn vertices(&self) -> usize {
        self.buffer.positions.len()
    }
}

/// Mesh one chunk at one level and score it.
///
/// # Scope, in two clauses, and neither of them is a tolerance
///
/// **Outside the solid.** `min` of exact primitives is the exact distance on the
/// exterior and an understatement inside (M-246), so an interior sample asserts
/// a sphere smaller than the one the surface actually has to touch — and that
/// error is a fixed number of world units, which is a *growing* number of cells
/// as the mesh refines. Counting interior samples would make the meter complain
/// harder the better the mesh got. So the census is exterior only, `f(p) ≥ 0`,
/// which is where the field this scene is built from is exact.
///
/// **The sphere has to fit inside the chunk.** A sample whose ball reaches
/// outside its own chunk has its nearest surface point in a *neighbour*, and no
/// amount of resolution here can touch that sphere. So a sample is in scope when
/// its closed ball lies inside the chunk box: the nearest surface point is then
/// inside the box too, and the chunk's own extraction is exactly what is
/// responsible for reaching it.
///
/// Within that scope the search is **exhaustive over every vertex of the
/// chunk**, with P-51's own window pruning and no cutoff.
fn measure_rung(cx: usize, cz: usize, level: usize) -> Rung {
    let n = LADDER[level];
    let h = CHUNK / f64::from(n - 1);
    let lo = [cx as f64 * CHUNK, 0.0, cz as f64 * CHUNK];
    let hi = [lo[0] + CHUNK, lo[1] + CHUNK, lo[2] + CHUNK];
    let count = n as usize;

    let shape = RuntimeShape3::new([n; 3]).expect("ladder resolutions fit u32");
    let mut buffer = MeshBuffer::<f64>::new();
    let started = Instant::now();
    if let Err(error) = MarchingCubes::<f64>::new().extract(&Scene, &shape, lo, h, &mut buffer) {
        panic!("E-313: chunk ({cx},{cz}) at {n}^3 did not extract: {error}");
    }
    let extract_ms = started.elapsed().as_secs_f64() * 1000.0;
    assert!(
        !buffer.positions.is_empty(),
        "E-313: chunk ({cx},{cz}) at {n}^3 produced no vertices; every count over it \
         would be a count of nothing"
    );

    let started = Instant::now();
    let inv_h = h.recip();
    let mut in_scope = 0_usize;
    let mut untouched = 0_usize;
    let mut worst_cells = 0.0_f64;
    for k in 0..count {
        let pz = lo[2] + k as f64 * h;
        for j in 0..count {
            let py = lo[1] + j as f64 * h;
            for i in 0..count {
                let px = lo[0] + i as f64 * h;
                let signed = Scene.sample([px, py, pz]);
                let clearance = (px - lo[0])
                    .min(hi[0] - px)
                    .min(py - lo[1])
                    .min(hi[1] - py)
                    .min(pz - lo[2])
                    .min(hi[2] - pz);
                if signed < 0.0 || signed > clearance {
                    continue;
                }
                let r = signed;
                in_scope += 1;

                // `min over v of | ‖v − p‖ − r |`, over every vertex. A vertex
                // can only beat the running best `t` when its distance lands in
                // `(r − t, r + t)`, so squared distances outside `((r−t)²,
                // (r+t)²)` are rejected without a square root — P-51's own
                // pruning, which is what makes an exhaustive search affordable.
                let mut best = f64::INFINITY;
                let mut window_lo = f64::NEG_INFINITY;
                let mut window_hi = f64::INFINITY;
                for v in &buffer.positions {
                    let dx = v[0] - px;
                    let dy = v[1] - py;
                    let dz = v[2] - pz;
                    let d2 = dx * dx + dy * dy + dz * dz;
                    if d2 > window_lo && d2 < window_hi {
                        let t = (d2.sqrt() - r).abs();
                        if t < best {
                            best = t;
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
                    untouched += 1;
                }
                if cells > worst_cells {
                    worst_cells = cells;
                }
            }
        }
    }
    let scan_ms = started.elapsed().as_secs_f64() * 1000.0;
    assert!(
        in_scope > 0,
        "E-313: chunk ({cx},{cz}) at {n}^3 owns no sphere; its fraction would have no \
         denominator"
    );

    Rung {
        n,
        h,
        samples: count * count * count,
        in_scope,
        untouched,
        worst_cells,
        buffer,
        extract_ms,
        scan_ms,
    }
}

// ─── the ledger cross-check ─────────────────────────────────────────────────

/// `p-51.csv`, addressed by column name so a reordered artefact is still the
/// same artefact.
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

    fn cell(&self, column: &str) -> Option<&'static str> {
        let f = self.column("field")?;
        let e = self.column("extractor")?;
        let c = self.column(column)?;
        self.rows
            .iter()
            .find(|row| {
                row.get(f).copied() == Some(LEDGER_FIELD)
                    && row.get(e).copied() == Some(LEDGER_EXTRACTOR)
            })?
            .get(c)
            .copied()
    }

    fn number(&self, column: &str) -> f64 {
        self.cell(column)
            .and_then(|text| text.parse().ok())
            .unwrap_or(f64::NAN)
    }
}

/// What holding the live predicate against `p-51.csv` produced.
struct LedgerCheck {
    csv_untouched: f64,
    live_untouched: f64,
    csv_per_1k: f64,
    live_per_1k: f64,
    csv_worst: f64,
    live_worst: f64,
    csv_vertices: f64,
    live_vertices: f64,
    reproduces: bool,
    ms: f64,
}

/// Run **P-51's predicate verbatim** — whole domain, no scope, every vertex —
/// on one reference field at the registered resolution, and hold the result
/// against the committed artefact.
///
/// This is what licenses the per-chunk numbers. The chunk scan adds a scope and
/// nothing else; if the predicate underneath it did not reproduce a committed
/// row, every fraction on the HUD would be a number with no provenance.
fn check_ledger() -> LedgerCheck {
    let started = Instant::now();
    let field = ThinPlate::<f64>::canonical();
    let (lo, hi) = field.domain();
    let h = (hi[0] - lo[0]) / f64::from(LEDGER_SAMPLES - 1);
    let n = LEDGER_SAMPLES as usize;

    let shape = RuntimeShape3::new([LEDGER_SAMPLES; 3]).expect("65^3 fits u32");
    let mut buffer = MeshBuffer::<f64>::new();
    if let Err(error) = MarchingCubes::<f64>::new().extract(&field, &shape, lo, h, &mut buffer) {
        panic!("E-313: the ledger row did not extract: {error}");
    }

    let inv_h = h.recip();
    let mut untouched = 0_u64;
    let mut worst = 0.0_f64;
    for k in 0..n {
        let pz = lo[2] + k as f64 * h;
        for j in 0..n {
            let py = lo[1] + j as f64 * h;
            for i in 0..n {
                let px = lo[0] + i as f64 * h;
                let r = field.sample([px, py, pz]).abs();
                let mut best = f64::INFINITY;
                let mut window_lo = f64::NEG_INFINITY;
                let mut window_hi = f64::INFINITY;
                for v in &buffer.positions {
                    let dx = v[0] - px;
                    let dy = v[1] - py;
                    let dz = v[2] - pz;
                    let d2 = dx * dx + dy * dy + dz * dz;
                    if d2 > window_lo && d2 < window_hi {
                        let t = (d2.sqrt() - r).abs();
                        if t < best {
                            best = t;
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
                    untouched += 1;
                }
                if cells > worst {
                    worst = cells;
                }
            }
        }
    }

    let samples = (n * n * n) as f64;
    let ledger = Ledger::parse();
    let live_per_1k = 1000.0 * untouched as f64 / samples;
    let csv_untouched = ledger.number("untouched");
    let csv_per_1k = ledger.number("untouched_per_1k");
    let csv_worst = ledger.number("worst_untouched_cells");
    let csv_vertices = ledger.number("vertices");
    let live_vertices = buffer.positions.len() as f64;
    // Four decimals on the rate and six on the worst case, because that is what
    // the CSV prints; the count is exact or it is not the same measurement.
    let reproduces = untouched as f64 == csv_untouched
        && (live_per_1k - csv_per_1k).abs() < 5e-5
        && (worst - csv_worst).abs() < 5e-7
        && live_vertices == csv_vertices;

    LedgerCheck {
        csv_untouched,
        live_untouched: untouched as f64,
        csv_per_1k,
        live_per_1k,
        csv_worst,
        live_worst: worst,
        csv_vertices,
        live_vertices,
        reproduces,
        ms: started.elapsed().as_secs_f64() * 1000.0,
    }
}

// ─── everything measured, once ──────────────────────────────────────────────

/// The four allocations this example compares.
///
/// [`Strategy::Matched`] is the one that makes the triangle comparison mean
/// something: the *same* range policy as [`Strategy::Range`], with its pixel
/// budget tightened until it clears the gate on every chunk. Comparing the
/// oracle against a range policy that is still losing detail would be comparing
/// two different pictures.
#[derive(Clone, Copy, PartialEq, Eq, Default)]
enum Strategy {
    #[default]
    Uniform,
    Range,
    Matched,
    Oracle,
}

/// Every strategy, in scoreboard order.
const STRATEGIES: [Strategy; 4] = [
    Strategy::Uniform,
    Strategy::Range,
    Strategy::Matched,
    Strategy::Oracle,
];

impl Strategy {
    fn index(self) -> usize {
        match self {
            Self::Uniform => 0,
            Self::Range => 1,
            Self::Matched => 2,
            Self::Oracle => 3,
        }
    }

    fn label(self) -> &'static str {
        match self {
            Self::Uniform => "UNIFORM 9^3",
            Self::Range => "RANGE LOD",
            Self::Matched => "RANGE, MATCHED",
            Self::Oracle => "ORACLE LOD",
        }
    }
}

/// The order the scripted beats visit the strategies.
const BEAT_ORDER: [Strategy; 3] = [Strategy::Uniform, Strategy::Range, Strategy::Oracle];

/// The whole ladder, the oracle's answer, the tightened range budget, and the
/// ledger verdict.
#[derive(Resource)]
struct Measured {
    /// `[chunk][level]`, chunk index `cz * CHUNKS_X + cx`.
    rungs: Vec<Vec<Rung>>,
    oracle: [usize; CHUNKS],
    /// Pixels per cell the range policy needs to clear the gate everywhere.
    matched_pixels: f32,
    /// Slowest single `(chunk, level)` census.
    ///
    /// On screen and in the audit because it is the honest limit on the whole
    /// idea: the meter needs no reference model, but an exhaustive tangency scan
    /// over a 65³ chunk against its own 6,000 vertices is a bake-time or
    /// background cost, not something to run inside a frame.
    worst_scan_ms: f64,
    ledger: LedgerCheck,
    ms: f64,
}

impl Measured {
    fn rung(&self, chunk: usize, level: usize) -> &Rung {
        &self.rungs[chunk][level]
    }
}

/// The first level whose untouched fraction is at or under the gate.
///
/// Exactly the runtime algorithm: climb from the cheapest rung, stop when the
/// chunk stops complaining. Falls back to the top of the ladder, which is the
/// most resolution there is to give.
fn oracle_level(ladder: &[Rung]) -> usize {
    ladder
        .iter()
        .position(|rung| rung.frac() <= GATE)
        .unwrap_or(LADDER.len() - 1)
}

/// Measure every chunk at every level, and cross-check the predicate, on
/// threads.
///
/// One thread per chunk plus one for the ledger row. Each owns its own buffers
/// and touches no shared mutable state, and they are joined in a fixed order —
/// as deterministic as the sequential version, and startup is one chunk long.
fn measure_everything() -> Measured {
    let started = Instant::now();
    let (rungs, ledger) = thread::scope(|scope| {
        let ledger = scope.spawn(check_ledger);
        let chunks: Vec<_> = (0..CHUNKS)
            .map(|chunk| {
                let cx = chunk % CHUNKS_X;
                let cz = chunk / CHUNKS_X;
                scope.spawn(move || {
                    (0..LADDER.len())
                        .map(|level| measure_rung(cx, cz, level))
                        .collect::<Vec<_>>()
                })
            })
            .collect();
        let rungs: Vec<Vec<Rung>> = chunks
            .into_iter()
            .map(|handle| handle.join().expect("chunk measurement panicked"))
            .collect();
        (rungs, ledger.join().expect("ledger check panicked"))
    });

    let mut oracle = [0_usize; CHUNKS];
    for (chunk, ladder) in rungs.iter().enumerate() {
        oracle[chunk] = oracle_level(ladder);
    }
    let worst_scan_ms = rungs
        .iter()
        .flatten()
        .map(|rung| rung.scan_ms)
        .fold(0.0_f64, f64::max);

    let mut measured = Measured {
        rungs,
        oracle,
        matched_pixels: DETAIL_PIXELS,
        worst_scan_ms,
        ledger,
        ms: 0.0,
    };
    measured.matched_pixels = matched_pixels(&measured);
    measured.ms = started.elapsed().as_secs_f64() * 1000.0;
    measured
}

/// The coarsest tenth of a pixel at which the range policy clears the gate on
/// every chunk.
///
/// Searched rather than solved, and downward from the policy a game would
/// actually ship, so the answer is *verified* by the same [`Tally`] the
/// scoreboard prints rather than derived from an inequality and hoped over.
/// Computed once, from [`reference_eye`], so the number the HUD quotes is a
/// fixed property of the scene and not something that slides while the camera
/// moves.
fn matched_pixels(measured: &Measured) -> f32 {
    let eye = reference_eye();
    let mut pixels = DETAIL_PIXELS;
    while pixels > 0.1 {
        let allocation = range_allocation(eye, pixels);
        if Tally::of(measured, &allocation).over == 0 {
            return pixels;
        }
        pixels -= 0.1;
    }
    0.1
}

// ─── the audit, on stdout, before any window ────────────────────────────────

/// Print the whole ladder and both allocations.
///
/// `println!` before `App::new()`, not `info!` inside a system: `DefaultPlugins`
/// builds the log subscriber *and* `WinitPlugin`'s event loop, and the latter
/// panics where there is no X display. An audit written with `info!` would be
/// unreachable exactly where a terminal is the only output.
fn audit(measured: &Measured) {
    println!("E-313 game_lod_oracle — the chunk that knows it is losing detail");
    println!(
        "  {CHUNKS} chunks, ladder {LADDER:?} samples/axis, gate {:.1}% of in-scope spheres, \
         miss gate {THRESHOLD_CELLS} cells",
        GATE * 100.0
    );
    println!("  field: union by min of exact primitives; census is exterior samples only");
    println!("  measured in {:.0} ms", measured.ms);
    println!();

    let ledger = &measured.ledger;
    println!(
        "p-51.csv {LEDGER_FIELD}/{LEDGER_EXTRACTOR} at {LEDGER_SAMPLES}^3, whole domain, no scope, \
         every vertex — {:.0} ms",
        ledger.ms
    );
    println!(
        "  vertices        csv {:>10.0}   live {:>10.0}",
        ledger.csv_vertices, ledger.live_vertices
    );
    println!(
        "  untouched       csv {:>10.0}   live {:>10.0}",
        ledger.csv_untouched, ledger.live_untouched
    );
    println!(
        "  untouched/1k    csv {:>10.4}   live {:>10.4}",
        ledger.csv_per_1k, ledger.live_per_1k
    );
    println!(
        "  worst cells     csv {:>10.6}   live {:>10.6}",
        ledger.csv_worst, ledger.live_worst
    );
    println!(
        "  verdict         {}",
        if ledger.reproduces {
            "REPRODUCES"
        } else {
            "DISAGREES"
        }
    );
    println!();

    println!(
        "chunk  lvl    n        h   samples  in_scope  untouched     frac   verts    tris  \
         worst_cells   mesh_ms   scan_ms"
    );
    for chunk in 0..CHUNKS {
        for level in 0..LADDER.len() {
            let rung = measured.rung(chunk, level);
            println!(
                "{},{}    {level}  {:>3}  {:.5}  {:>8}  {:>8}  {:>9}  {:>6.2}%  {:>6}  {:>6}  \
                 {:>11.4}  {:>8.2}  {:>8.2}",
                chunk % CHUNKS_X,
                chunk / CHUNKS_X,
                rung.n,
                rung.h,
                rung.samples,
                rung.in_scope,
                rung.untouched,
                rung.frac() * 100.0,
                rung.vertices(),
                rung.triangles(),
                rung.worst_cells,
                rung.extract_ms,
                rung.scan_ms,
            );
        }
    }
    println!();

    let eye = reference_eye();
    let allocations = every_allocation(measured, eye);
    println!("allocation per chunk (level index into the ladder)");
    for cz in 0..CHUNKS_Z {
        for cx in 0..CHUNKS_X {
            let chunk = cz * CHUNKS_X + cx;
            println!(
                "  ({cx},{cz})  range {:>5.2}   uniform {}   range-lod {}   matched {}   oracle {}",
                chunk_distance(chunk, eye),
                allocations[0][chunk],
                allocations[1][chunk],
                allocations[2][chunk],
                allocations[3][chunk],
            );
        }
    }
    println!();

    for strategy in STRATEGIES {
        let tally = Tally::of(measured, &allocations[strategy.index()]);
        println!(
            "{:<15}  {:>8} tris  {:>7} verts   worst chunk {:>6.2}%   over gate {:>2} of {CHUNKS}",
            strategy.label(),
            tally.triangles,
            tally.vertices,
            tally.worst * 100.0,
            tally.over,
        );
    }
    println!(
        "range policy: {DETAIL_PIXELS} px/cell as shipped, {:.1} px/cell to clear the gate \
         everywhere",
        measured.matched_pixels
    );
    println!();
    // Named rather than smoothed away: a chunk whose fraction rises with
    // resolution makes `oracle_level`'s first-crossing rule a genuine choice
    // rather than a way of writing "the minimum passing level". Only rises worth
    // a whole percentage point are reported — below that the two readings are the
    // same statement about a denominator that changed by two orders of magnitude.
    for chunk in 0..CHUNKS {
        let ladder = &measured.rungs[chunk];
        let step = (1..ladder.len())
            .find(|&level| ladder[level].frac() > ladder[level - 1].frac() + MONOTONE_SLACK);
        if let Some(level) = step {
            println!(
                "non-monotone: ({},{}) reads {:.2}% at {}^3 and {:.2}% at {}^3 \
                 (in scope {} then {})",
                chunk % CHUNKS_X,
                chunk / CHUNKS_X,
                ladder[level - 1].frac() * 100.0,
                ladder[level - 1].n,
                ladder[level].frac() * 100.0,
                ladder[level].n,
                ladder[level - 1].in_scope,
                ladder[level].in_scope,
            );
        }
    }
    println!(
        "worst single census: {:.0} ms — a bake or a background budget, not a frame",
        measured.worst_scan_ms
    );
    println!();

    // The gate is claimed to be the smallest whole percent the top of the ladder
    // clears on every chunk. Both halves of that are asserted here, so a change
    // to the scene cannot quietly leave the constant behind.
    let top = LADDER.len() - 1;
    let ceiling = (0..CHUNKS)
        .map(|chunk| measured.rung(chunk, top).frac())
        .fold(0.0_f64, f64::max);
    assert!(
        ceiling <= GATE,
        "E-313: GATE {GATE} does not clear the top of the ladder — worst chunk reads \
         {ceiling:.4} at {}^3, so the oracle's answer there is 'the ladder is not enough'",
        LADDER[top]
    );
    assert!(
        ceiling > GATE - 0.01,
        "E-313: GATE {GATE} is looser than the measurement needs — the top of the ladder \
         clears {ceiling:.4} everywhere, so {:.2} would do",
        (ceiling * 100.0).ceil() / 100.0
    );
    println!(
        "gate check: top of the ladder clears {:.2}% on its worst chunk, gate {:.0}%",
        ceiling * 100.0,
        GATE * 100.0
    );
    println!();
}

// ─── allocations ────────────────────────────────────────────────────────────

/// Where the camera sits when a capture starts, in world units.
///
/// The range allocation is recomputed live from the camera, so orbiting really
/// does move the LOD; this is the eye the audit, the matched pixel budget and
/// the captured frames are all taken from.
fn reference_eye() -> Vec3 {
    let direction = Vec3::new(
        VIEW_YAW.cos() * VIEW_PITCH.cos(),
        VIEW_PITCH.sin(),
        VIEW_YAW.sin() * VIEW_PITCH.cos(),
    );
    VIEW_FOCUS + direction * VIEW_RADIUS
}

/// Centre of a chunk, in world units.
fn chunk_centre(chunk: usize) -> Vec3 {
    let cx = chunk % CHUNKS_X;
    let cz = chunk / CHUNKS_X;
    Vec3::new(
        (cx as f64 * CHUNK + CHUNK * 0.5) as f32,
        (CHUNK * 0.5) as f32,
        (cz as f64 * CHUNK + CHUNK * 0.5) as f32,
    )
}

fn chunk_distance(chunk: usize, eye: Vec3) -> f32 {
    eye.distance(chunk_centre(chunk))
}

/// The screen-space budget the range policy is shipped with, in pixels per cell.
///
/// This is the LOD a game actually ships: pick the coarsest level whose
/// projected cell stays under a fixed number of pixels. One knob, monotone in
/// range — which is the whole point, because *need* is not monotone in range.
///
/// **Derived, not dialled.** At Bevy's default 45° vertical field of view a
/// 720-pixel-high frame subtends `2·tan(22.5°) = 0.8284` world units per unit
/// of range, so one pixel is [`UNITS_PER_PIXEL`] units per unit of range. Ten
/// pixels per cell is a mainstream terrain target and it is the number quoted on
/// the HUD, so a reader can disagree with the policy rather than with an unnamed
/// constant.
const DETAIL_PIXELS: f32 = 10.0;

/// World units of vertical frame per pixel, per unit of range, at 720p and 45°.
const UNITS_PER_PIXEL: f32 = 0.828_427_1 / 720.0;

/// Coarsest level whose cell stays inside a budget of `pixels` at this range.
fn range_level(range: f32, pixels: f32) -> usize {
    let target = pixels * UNITS_PER_PIXEL * range;
    for (level, n) in LADDER.iter().enumerate() {
        if (CHUNK / f64::from(n - 1)) as f32 <= target {
            return level;
        }
    }
    LADDER.len() - 1
}

/// The range policy's allocation at a given pixel budget.
fn range_allocation(eye: Vec3, pixels: f32) -> [usize; CHUNKS] {
    let mut levels = [0_usize; CHUNKS];
    for (chunk, level) in levels.iter_mut().enumerate() {
        *level = range_level(chunk_distance(chunk, eye), pixels);
    }
    levels
}

/// All four allocations, in [`Strategy::index`] order.
fn every_allocation(measured: &Measured, eye: Vec3) -> [[usize; CHUNKS]; 4] {
    [
        [0_usize; CHUNKS],
        range_allocation(eye, DETAIL_PIXELS),
        range_allocation(eye, measured.matched_pixels),
        measured.oracle,
    ]
}

/// What one allocation costs and what it leaves broken.
#[derive(Clone, Copy, Default)]
struct Tally {
    triangles: usize,
    vertices: usize,
    worst: f64,
    worst_chunk: usize,
    over: usize,
}

impl Tally {
    fn of(measured: &Measured, allocation: &[usize; CHUNKS]) -> Self {
        let mut tally = Self::default();
        for (chunk, &level) in allocation.iter().enumerate() {
            let rung = measured.rung(chunk, level);
            tally.triangles += rung.triangles();
            tally.vertices += rung.vertices();
            let frac = rung.frac();
            if frac > GATE {
                tally.over += 1;
            }
            if frac > tally.worst {
                tally.worst = frac;
                tally.worst_chunk = chunk;
            }
        }
        tally
    }
}

// ─── the beats ──────────────────────────────────────────────────────────────

/// Captured frames each of the first three beats holds.
///
/// 16 + 22 + 22 = 60, so `record_gif.sh`'s default 80 leaves 20 frames for the
/// closing A/B, which at [`AB_FRAMES`] is three cuts. **The A/B is between the
/// range policy tightened to the same detail bar and the oracle**, so the two
/// sides look the same and the only thing that moves is the triangle count. A
/// viewer reads a difference from a repeat, not from a single cut.
const BEAT_FRAMES: [u32; 3] = [16, 22, 22];

/// Captured frames per side of the closing A/B flip.
const AB_FRAMES: u32 = 6;

/// Seconds per beat, and per A/B side, when nobody is capturing.
const BEAT_SECONDS: [f32; 3] = [3.2, 4.4, 4.4];
const AB_SECONDS: f32 = 1.5;

/// Which beat and which strategy, from a captured-frame index.
fn beat_at_frame(frame: u32) -> (usize, Strategy) {
    let mut acc = 0;
    for (index, length) in BEAT_FRAMES.iter().enumerate() {
        if frame < acc + length {
            return (index, BEAT_ORDER[index]);
        }
        acc += length;
    }
    let side = ((frame - acc) / AB_FRAMES) % 2;
    (
        BEAT_FRAMES.len(),
        if side == 0 {
            Strategy::Matched
        } else {
            Strategy::Oracle
        },
    )
}

/// The same schedule on a clock, for the interactive run.
fn beat_at_seconds(seconds: f32) -> (usize, Strategy) {
    let mut acc = 0.0;
    for (index, length) in BEAT_SECONDS.iter().enumerate() {
        if seconds < acc + length {
            return (index, BEAT_ORDER[index]);
        }
        acc += length;
    }
    let side = (((seconds - acc) / AB_SECONDS) as u32) % 2;
    (
        BEAT_SECONDS.len(),
        if side == 0 {
            Strategy::Matched
        } else {
            Strategy::Oracle
        },
    )
}

/// Which beat is on screen, and whether a key pinned it.
#[derive(Resource, Default)]
struct Beat {
    index: usize,
    strategy: Strategy,
    elapsed: f32,
    pinned: Option<Strategy>,
}

/// The current allocation and the cost of all four, recomputed each frame.
#[derive(Resource, Default)]
struct Live {
    level: [usize; CHUNKS],
    tally: [Tally; 4],
    eye: Vec3,
}

// ─── the picture ────────────────────────────────────────────────────────────

/// The view the captured frames are taken from: eye `(3.40, 3.60, −2.40)`,
/// pitched 41.3° down the `z` axis.
///
/// **Every number here was measured off the frames rather than chosen.**
///
/// *Pitch.* At 17° the terrain projected to a thin horizontal band with the near
/// chunk row filling the bottom edge-on, and neither the 4×4 chunk grid nor the
/// heat map over it could be read at all. 41° is the three-quarter view that
/// makes the grid a receding quad while still leaving the spires standing
/// against the sky rather than seen from above.
///
/// *Focus 1.15 below the terrain.* The scoreboard owns the bottom 236 pixels and
/// an orbit camera always looks at its focus, so the only way to sit the subject
/// high in frame is to aim under it. This puts the terrain's centre 26% of a
/// half-height above the frame centre, which clears the scoreboard with the far
/// row still inside the top edge.
///
/// *Focus at `x = 3.40` rather than the grid's own `2.00`.* The harness HUD's
/// plate is 516 pixels wide, and centring the grid put a whole spire cluster
/// behind it — the frame lost one of the two subjects the beat is about. Note
/// the direction: this camera looks along **+z**, so screen-right is world
/// **−x** and moving the focus to larger `x` pushes the terrain right. The first
/// attempt moved it to `0.90` and shoved the grid further under the plate.
/// Measured at `3.40`: the terrain spans screen `x = 536` to `1228`, so the
/// plate's 520 is clear. It also breaks the range policy's symmetry in `x`,
/// which is a bonus rather than a cost: chunks with identical geometry now get
/// different levels for no reason but range.
const VIEW_FOCUS: Vec3 = Vec3::new(3.40, -0.70, 2.50);
const VIEW_RADIUS: f32 = 6.52;
const VIEW_YAW: f32 = -std::f32::consts::FRAC_PI_2;
const VIEW_PITCH: f32 = 0.7203;

/// The heat ramp, as terrain rather than as a debug palette.
///
/// **The hue breaks exactly on the gate**, and that is the whole design. A
/// continuous ramp read as "everything is a bit warm" — the oracle's own answer
/// left every chunk sitting *just under* the bar, because that is what climbing
/// until you clear it means, so on a smooth ramp a passing scene photographed
/// amber and contradicted a scoreboard that said nothing was losing detail. Two
/// families instead: cool below the gate, warm above it. Depth inside each band
/// still carries the magnitude — deep moss at zero, pale sage at the gate, amber
/// just over it, bare rock at three times it.
const MOSS_DEEP: [f32; 3] = [0.13, 0.34, 0.21];
const MOSS_PALE: [f32; 3] = [0.52, 0.70, 0.43];
const AMBER: [f32; 3] = [0.88, 0.62, 0.14];
const ROCK: [f32; 3] = [0.78, 0.11, 0.06];

/// Chunk colour from untouched fraction.
fn heat(frac: f64) -> Color {
    let t = (frac / GATE) as f32;
    let (a, b, u) = if t <= 1.0 {
        (MOSS_DEEP, MOSS_PALE, t.clamp(0.0, 1.0))
    } else {
        (AMBER, ROCK, ((t - 1.0) * 0.5).clamp(0.0, 1.0))
    };
    Color::srgb(
        a[0] + (b[0] - a[0]) * u,
        a[1] + (b[1] - a[1]) * u,
        a[2] + (b[2] - a[2]) * u,
    )
}

/// Every `(chunk, level)` mesh and material, uploaded once.
#[derive(Resource)]
struct Prebuilt {
    mesh: Vec<Vec<Handle<Mesh>>>,
    material: Vec<Vec<Handle<StandardMaterial>>>,
}

#[derive(Component)]
struct Chunk(usize);

#[derive(Component)]
struct Banner;

#[derive(Component)]
struct Score;

#[derive(Component)]
struct Backdrop;

/// Size of the plate behind the harness HUD, in logical pixels.
///
/// Sized against what `report` emits: 22 lines at 13px, nothing wider than 64
/// characters, and 13px of Bevy's default mono measures about 7.8px wide and
/// 15.6px tall.
const HUD_PLATE: Vec2 = Vec2::new(516.0, 364.0);

/// Size of the plate behind the scoreboard.
///
/// Eight lines at 22px — a caption, a blank, four strategies, a blank, the
/// verdict — and nothing wider than 72 characters at 13.2px each. **A plate
/// rather than brighter text**: the scoreboard sits over lit terrain, and white
/// glyphs over a sunlit amber rampart were unreadable on the first frames.
const SCORE_PLATE: Vec2 = Vec2::new(996.0, 236.0);

fn main() {
    let measured = measure_everything();
    audit(&measured);

    App::new()
        .add_plugins(DefaultPlugins.set(WindowPlugin {
            primary_window: Some(Window {
                title: "isomesh — detail oracle LOD".into(),
                ..default()
            }),
            ..default()
        }))
        .add_plugins(CommonPlugin)
        .insert_resource(measured)
        .init_resource::<Beat>()
        .init_resource::<Live>()
        .add_systems(Startup, setup)
        .add_systems(Update, (controls, advance, allocate, apply, report).chain())
        .run();
}

#[allow(clippy::too_many_arguments)]
fn setup(
    mut commands: Commands,
    mut meshes: ResMut<Assets<Mesh>>,
    mut materials: ResMut<Assets<StandardMaterial>>,
    mut camera: Query<&mut OrbitCamera>,
    mut light: Query<(&mut DirectionalLight, &mut Transform)>,
    measured: Res<Measured>,
) {
    for mut orbit in &mut camera {
        orbit.focus = VIEW_FOCUS;
        orbit.radius = VIEW_RADIUS;
        orbit.yaw = VIEW_YAW;
        orbit.pitch = VIEW_PITCH;
    }

    // One strong key light, low and raking across the grid, and modest ambient.
    // The subject is a spire `0.09` units across: flat lighting turns it into a
    // coloured smear, and the whole point of beat 3 is that you can watch it
    // become a pole. The long shadows are load-bearing too — a spire that has
    // vanished at 9^3 takes its shadow with it.
    for (mut key, mut transform) in &mut light {
        key.illuminance = 24_000.0;
        key.shadow_maps_enabled = true;
        *transform =
            Transform::from_xyz(-2.6, 3.4, -0.4).looking_at(Vec3::new(2.2, 0.30, 2.4), Vec3::Y);
    }
    commands.insert_resource(GlobalAmbientLight {
        brightness: 110.0,
        ..default()
    });
    // A sky rather than the default grey: the heat ramp is warm at both ends, so
    // a cool background is what keeps a moss-green chunk reading as cool.
    commands.insert_resource(ClearColor(Color::srgb(0.09, 0.11, 0.17)));

    let mut mesh = Vec::with_capacity(CHUNKS);
    let mut material = Vec::with_capacity(CHUNKS);
    for chunk in 0..CHUNKS {
        let mut per_level_mesh = Vec::with_capacity(LADDER.len());
        let mut per_level_material = Vec::with_capacity(LADDER.len());
        for level in 0..LADDER.len() {
            let rung = measured.rung(chunk, level);
            per_level_mesh.push(meshes.add(to_mesh(&rung.buffer)));
            per_level_material.push(materials.add(StandardMaterial {
                base_color: heat(rung.frac()),
                perceptual_roughness: 0.82,
                ..default()
            }));
        }
        mesh.push(per_level_mesh);
        material.push(per_level_material);
        commands.spawn((
            Mesh3d(mesh[chunk][0].clone()),
            MeshMaterial3d(material[chunk][0].clone()),
            Chunk(chunk),
            DemoMesh,
        ));
    }
    commands.insert_resource(Prebuilt { mesh, material });

    // Behind the harness's HUD text and behind the scoreboard. `GlobalZIndex(-1)`
    // rather than spawn order, because the HUD text entity belongs to
    // `CommonPlugin` and was spawned first.
    for (offset, size) in [
        ((Val::Px(4.0), Val::Auto), HUD_PLATE),
        ((Val::Auto, Val::Px(4.0)), SCORE_PLATE),
    ] {
        commands.spawn((
            Node {
                position_type: PositionType::Absolute,
                top: offset.0,
                bottom: offset.1,
                left: Val::Px(4.0),
                width: Val::Px(size.x),
                height: Val::Px(size.y),
                ..default()
            },
            BackgroundColor(Color::srgba(0.03, 0.04, 0.06, 0.78)),
            GlobalZIndex(-1),
            Backdrop,
        ));
    }

    commands.spawn((
        Text::new(""),
        TextFont {
            font_size: FontSize::Px(26.0),
            ..default()
        },
        TextColor(Color::srgb(0.99, 0.97, 0.90)),
        Node {
            position_type: PositionType::Absolute,
            top: Val::Px(14.0),
            right: Val::Px(22.0),
            ..default()
        },
        Banner,
    ));

    commands.spawn((
        Text::new(""),
        TextFont {
            font_size: FontSize::Px(22.0),
            ..default()
        },
        TextColor(Color::srgb(0.96, 0.95, 0.92)),
        Node {
            position_type: PositionType::Absolute,
            bottom: Val::Px(18.0),
            left: Val::Px(18.0),
            ..default()
        },
        Score,
    ));
}

/// The `f64` extraction as a Bevy mesh.
///
/// Cast rather than re-extracted in `f32`: the fractions on screen are `f64`
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
    for t in buffer.indices.chunks_exact(3) {
        let (Some(a), Some(b), Some(c)) = (t.first(), t.get(1), t.get(2)) else {
            continue;
        };
        builder.triangle(*a, *b, *c);
    }
    builder.into_mesh()
}

/// The keyboard. Ignored under capture, which drives itself.
fn controls(keys: Res<ButtonInput<KeyCode>>, capture: Res<Capture>, mut beat: ResMut<Beat>) {
    if capture.is_active() {
        return;
    }
    for (key, pin) in [
        (KeyCode::Digit1, Some(Strategy::Uniform)),
        (KeyCode::Digit2, Some(Strategy::Range)),
        (KeyCode::Digit3, Some(Strategy::Matched)),
        (KeyCode::Digit4, Some(Strategy::Oracle)),
        (KeyCode::Digit5, None),
    ] {
        if keys.just_pressed(key) {
            beat.pinned = pin;
        }
    }
}

/// Decide which beat is on screen.
///
/// Under capture both the beat and the A/B side come off the captured-frame
/// counter, so a clip of any length shows the whole argument in order and
/// needs no keyboard.
fn advance(time: Res<Time>, capture: Res<Capture>, flags: Res<ViewFlags>, mut beat: ResMut<Beat>) {
    if let Some(pinned) = beat.pinned {
        beat.strategy = pinned;
        beat.index = BEAT_ORDER
            .iter()
            .position(|entry| *entry == pinned)
            .unwrap_or(BEAT_ORDER.len());
        return;
    }
    if capture.is_active() {
        let (index, strategy) = beat_at_frame(capture.taken);
        beat.index = index;
        beat.strategy = strategy;
        return;
    }
    if !flags.paused {
        beat.elapsed += time.delta_secs();
    }
    let (index, strategy) = beat_at_seconds(beat.elapsed);
    beat.index = index;
    beat.strategy = strategy;
}

/// Score all four allocations against the live camera.
fn allocate(
    measured: Res<Measured>,
    beat: Res<Beat>,
    camera: Query<&Transform, With<Camera3d>>,
    mut live: ResMut<Live>,
) {
    let eye = camera
        .iter()
        .next()
        .map_or_else(reference_eye, |transform| transform.translation);
    live.eye = eye;
    let allocations = every_allocation(&measured, eye);
    for (index, allocation) in allocations.iter().enumerate() {
        live.tally[index] = Tally::of(&measured, allocation);
    }
    live.level = allocations[beat.strategy.index()];
}

/// Swap in the meshes and colours the current allocation asks for.
fn apply(
    live: Res<Live>,
    prebuilt: Res<Prebuilt>,
    mut chunks: Query<(&Chunk, &mut Mesh3d, &mut MeshMaterial3d<StandardMaterial>)>,
) {
    for (chunk, mut mesh, mut material) in &mut chunks {
        let level = live.level[chunk.0];
        let wanted = &prebuilt.mesh[chunk.0][level];
        if mesh.0.id() != wanted.id() {
            mesh.0 = wanted.clone();
            material.0 = prebuilt.material[chunk.0][level].clone();
        }
    }
}

/// `12345` as `12,345`.
fn thousands(value: usize) -> String {
    let digits = value.to_string();
    let mut out = String::with_capacity(digits.len() + digits.len() / 3);
    for (index, ch) in digits.chars().enumerate() {
        if index > 0 && (digits.len() - index).is_multiple_of(3) {
            out.push(',');
        }
        out.push(ch);
    }
    out
}

/// The banner, the scoreboard, and the harness HUD.
#[allow(clippy::too_many_arguments)]
fn report(
    beat: Res<Beat>,
    live: Res<Live>,
    measured: Res<Measured>,
    flags: Res<ViewFlags>,
    mut stats: ResMut<DemoStats>,
    mut banner: Query<&mut Text, (With<Banner>, Without<Score>)>,
    mut score: Query<&mut Text, (With<Score>, Without<Banner>)>,
    mut backdrop: Query<&mut Visibility, With<Backdrop>>,
) {
    let current = beat.strategy.index();
    let tally = live.tally[current];
    let matched = live.tally[Strategy::Matched.index()];
    let oracle = live.tally[Strategy::Oracle.index()];

    // ASCII only, in every string that reaches a `Text`. Bevy's default font is
    // a FiraMono subset with no U+00B3 and no U+00B7: the first cut of this HUD
    // rendered "UNIFORM 9³" and "·" as tofu boxes on the captured frames.
    let shipped = format!("BEAT 2/4   RANGE LOD @ {DETAIL_PIXELS:.0} PX/CELL");
    let ab = format!(
        "BEAT 4/4   A/B: {}",
        match beat.strategy {
            Strategy::Oracle => "ORACLE",
            _ => "RANGE, MATCHED",
        }
    );
    let headline = match beat.index {
        0 => "BEAT 1/4   UNIFORM 9^3",
        1 => shipped.as_str(),
        2 => "BEAT 3/4   ORACLE LOD",
        _ => ab.as_str(),
    };
    for mut text in &mut banner {
        text.0 = if flags.hud {
            headline.to_string()
        } else {
            String::new()
        };
    }

    let caption = match beat.index {
        0 => "every chunk at 9^3: the spires are stubs, and their own field says so",
        1 => "resolution from camera range alone: it cannot see what a chunk needs",
        2 => "each chunk climbed the ladder until its own field stopped complaining",
        _ => "same detail bar, half the triangles: watch the count, not the terrain",
    };
    let mut board = format!("{caption}\n\n");
    for strategy in STRATEGIES {
        let entry = live.tally[strategy.index()];
        let budget = match strategy {
            Strategy::Range => format!("{DETAIL_PIXELS:>4.1} px/cell"),
            Strategy::Matched => format!("{:>4.1} px/cell", measured.matched_pixels),
            _ => String::new(),
        };
        board.push_str(&format!(
            "{} {:<15} {:>12} {:>8} tris  {:>2} chunks losing detail\n",
            if strategy == beat.strategy { ">" } else { " " },
            strategy.label(),
            budget,
            thousands(entry.triangles),
            entry.over,
        ));
    }
    let ratio = if oracle.triangles == 0 {
        0.0
    } else {
        matched.triangles as f64 / oracle.triangles as f64
    };
    board.push_str(&format!(
        "\nORACLE WINS: the same detail bar as range-at-{:.1}px, {ratio:.2}x fewer triangles",
        measured.matched_pixels
    ));
    for mut text in &mut score {
        text.0 = if flags.hud {
            board.clone()
        } else {
            String::new()
        };
    }

    for mut visibility in &mut backdrop {
        *visibility = if flags.hud {
            Visibility::Inherited
        } else {
            Visibility::Hidden
        };
    }

    let ledger = &measured.ledger;
    let worst = measured.rung(tally.worst_chunk, live.level[tally.worst_chunk]);
    stats.title = format!(
        "E-313  detail oracle  --  {CHUNKS} chunks, {} levels",
        LADDER.len()
    );
    stats.triangles = tally.triangles;
    stats.vertices = tally.vertices;
    stats.extract_ms = measured.ms;
    stats.extra = vec![
        format!(
            "{:>13}  {:.0}% of a chunk's in-scope spheres untouched",
            "gate",
            GATE * 100.0
        ),
        format!(
            "{:>13}  {THRESHOLD_CELLS} cells (P-51, registered)",
            "miss gate"
        ),
        format!(
            "{:>13}  exterior, and the ball fits inside its own chunk",
            "in scope"
        ),
        format!(
            "{:>13}  exhaustive over every chunk vertex, no cutoff",
            "search"
        ),
        String::new(),
        format!(
            "{:>13}  {:.2}% at ({},{}) n={}",
            "worst chunk",
            tally.worst * 100.0,
            tally.worst_chunk % CHUNKS_X,
            tally.worst_chunk / CHUNKS_X,
            worst.n,
        ),
        format!("{:>13}  {} of {CHUNKS} chunks", "over gate", tally.over),
        String::new(),
        format!(
            "{:>13}  {:.0} ms on the worst rung -- a bake, not a frame",
            "oracle cost", measured.worst_scan_ms
        ),
        String::new(),
        format!("p-51.csv {LEDGER_FIELD}/{LEDGER_EXTRACTOR} at {LEDGER_SAMPLES}^3"),
        format!(
            "{:>13}  csv {:.4}   live {:.4}   {}",
            "untouched/1k",
            ledger.csv_per_1k,
            ledger.live_per_1k,
            if ledger.reproduces {
                "REPRODUCES"
            } else {
                "DISAGREES"
            },
        ),
        String::new(),
        "[1] uniform  [2] range  [3] matched  [4] oracle  [5] beats".to_string(),
    ];
}
