//! E-306 — the aperture: not *is it connected* but *how big a thing fits*.
//!
//! ```bash
//! cd bevy_isomesh && cargo run --example aperture_gate --release
//! ```
//!
//! **Always `--release`.** A debug build samples 65³ points three times per
//! rebuild and the sweep crawls.
//!
//! `1` the drilled slab, `2` the uncapped gyroid, `3` the capped gyroid, `X`
//! restarts the sweep, `Space` freezes it. The rest are the shared keys — `W`
//! wireframe, `G` domain box, `F12` screenshot.
//!
//! Under `ISOMESH_CAPTURE` it needs no keyboard: the channel widens across the
//! first 62% of the captured frames, and the field then steps to the capped
//! gyroid and finally to the uncapped one, so a capture of *any* length shows the
//! gate open and both gyroids. `ISOMESH_FIELD=1` pins one field even under
//! capture.
//!
//! ```bash
//! # Measured 1.13 MB, inside the 0.7-4.8 MB the committed clips sit within.
//! # The camera never moves and only the ball and the hole do, which is exactly
//! # what a GIF's inter-frame compression wants, so this affords the script's
//! # default palette and error-diffusion dither rather than needing the cheap
//! # ordered one. `COLORS=64` takes it to 0.93 MB if a smaller file is wanted.
//! ISOMESH_CAPTURE_FRAMES=90 ISOMESH_CAPTURE_EVERY=3 ISOMESH_WINDOW=1280x720 \
//!   FPS=15 ./scripts/record_gif.sh aperture_gate docs/gifs/e306.gif
//! ```
//!
//! **1280x720 rather than smaller.** Every HUD line is kept inside 76 characters
//! so nothing wraps at 640 wide either, but at 360 tall the HUD reaches the
//! bottom edge and crosses the readout. A 6-frame 640x360 smoke capture still
//! shows the red-to-green transition — the gate opens on captured frame 2 of 6 —
//! and comes back at 0.42 MB, legitimately under the band at five frames.
//!
//! Demonstrates **M-346 / P-49**.
//!
//! # The question the crate could not answer
//!
//! `isomesh::validate::sealing` and `chunk::AirWorld` answer *are these two
//! places connected*, which is a boolean. A game asks *can the player get
//! through*, which is a **number**: for a pair of chunk faces, the **maximum over
//! air paths of the minimum distance-to-solid along the path** — the widest
//! bottleneck. Connectivity is the aperture being positive; the aperture is
//! everything connectivity throws away.
//!
//! That distinction is on screen twice. The ball is a fixed body radius and the
//! hole widens under it: the pair is *connected* the whole time and the ball
//! still does not fit until the number crosses. And in the 6×6 grid a cell is
//! **amber** when the pair is connected but too tight, **green** when the body
//! fits, dark when there is no air path at all. A boolean has only the two
//! colours it cannot tell apart.
//!
//! # The ordering is what makes it deterministic
//!
//! [`Aperture::solve`] is a monotone union-find over **air** samples — `value >
//! 0.0`, the strict complement of `cube.rs::is_inside`'s `value < 0.0`, so an
//! aperture is a strictly positive clearance. The samples are processed in
//! **descending `(field value, grid index)`**, compared with
//! [`f64::total_cmp`]. That is a **total order** on the air set: no two entries
//! compare equal, because the index breaks every value tie. So there is no PRNG,
//! no atomics, no `HashMap` and no tie broken by allocation address — the answer
//! is a function of the sample array alone, and the same array gives the same
//! matrix on every machine and every run. The `determinism` line on the HUD is
//! that claim checked live, two solves compared entry by entry.
//!
//! Each sample is activated in that order and unioned with its already-active
//! **6-neighbours**; every component carries a 6-bit mask of the grid faces it
//! touches. The first moment a component's mask contains both faces of a pair,
//! the current sample's value *is* that pair's aperture — first in descending
//! order is the largest bottleneck, which is the definition.
//!
//! 6-connectivity rather than 18 or 26: a diagonal step between two
//! face-adjacent solids passes through material, and a clearance a game gates
//! movement on must not claim it.
//!
//! P-49's bench spells the same total order as an ascending integer sort of
//! `(!value.to_bits(), index)` — for a positive `f64` the bit pattern is monotone
//! in the value, so the complement sorts descending exactly, with integer
//! comparisons. This file uses `total_cmp` because the order is the thing being
//! explained. It is the identical order and a hair slower, which is why the live
//! cost ratio sits a little above the bench's.
//!
//! # The early exit is sound, and that is checked rather than asserted
//!
//! Once all 15 pairs have an aperture no later sample can revise one, because
//! each was fixed at the first — and therefore highest — value that connected its
//! pair. So the loop may stop. Every rebuild here solves **both** ways and
//! compares all 15 entries; `early exit sound` on the HUD is that comparison, and
//! the no-early-exit time is reported beside the early-exit one as the worst case
//! to budget for. On the gyroids the exit fires after about **5%** of the air
//! samples, which is the whole of the gap between the two numbers.
//!
//! # Why the slab is the fixture with nowhere to hide
//!
//! A [`BoxExact`] big enough to swallow the entire domain, with a
//! [`Capsule`](isomesh::brush::Capsule) of radius `r` subtracted along `x`.
//! Subtraction is `max(field, −shape)`, the box term is ≤ −2 everywhere in the
//! domain, and `Capsule` is an exact distance field — so inside the channel the
//! value is *exactly* `r − ρ` at radius `ρ` from the axis. Therefore:
//!
//! - the deepest air is the axis, at exactly `r`;
//! - the axis is a **sample line**, because 65 samples over `[−2, 2]` put a
//!   sample at `0.0` on every axis, so the exact value `r` is attained rather
//!   than approached;
//! - `r` is a whole number of cells times `0.0625`, a power of two, so `r` in
//!   world units and `r` in cells are both exact and so is the quotient;
//! - the channel leaves through both `x` faces and reaches no other, because the
//!   `±y`/`±z` walls are 32 cells from the axis and the sweep stops at 12.
//!
//! Ground truth is `r` with **no discretisation slack** and exactly one reachable
//! pair. M-346 measured error `0.000000` at r = 2, 4 and 8 cells with zero
//! falsely-reachable pairs; this file replays those three at startup, logs the
//! errors and puts them on the HUD, and the live `exactness` line re-checks the
//! same identity at whatever radius the sweep is currently at.
//!
//! The **air census** line is a second, independent check on the same fixture: a
//! sample is air exactly when `j² + k² < r²` in integer cell units, over all 65
//! planes along `x`, so the air count has a closed form in integer arithmetic and
//! the union-find's count must equal it.
//!
//! # The gyroid is the *this is not a toy fixture* half
//!
//! At 65³ all six faces are mutually reachable on both gyroids and all fifteen
//! apertures are positive — **6.39 to 6.83 cells** on the uncapped one, 17.8 to
//! 28.0 on the capped one. The capped gyroid's faces are joined by the *exterior
//! of the cap* rather than by the gyroid's channels (`capped_gyroid` is
//! `max(gyroid, sphere(6))` over `[−7, 7]³`, so everything outside radius 6 is
//! air, including all six domain faces), which is why the uncapped one is here
//! too: no shell, so the only route between two faces is the bicontinuous channel
//! network itself. Both are sampled on the capped field's `[−7, 7]³` domain, as
//! P-49 measured them.
//!
//! # The cost, which is what makes it shippable
//!
//! Reported the expensive way, because that is the reading that makes the claim
//! harder to pass: `whole` includes **sampling the grid**, since Marching Cubes'
//! own time includes its sampling too and the comparison should be like for like.
//! `marginal` is the solve alone — the number a chunk pipeline would actually
//! budget, because a mesher already has the samples. M-346 measured the worst
//! case at **1.108× a Marching Cubes extraction whole** and **0.391× marginal**,
//! and the HUD carries both the live numbers and those.
//!
//! ## The live ratio is not the bench's, and the reason is the build
//!
//! Measured here, one shot per rebuild, against P-49's median of seven timed
//! reps: the **gyroid** whole ratios reproduce — `0.850×` against `0.829×`
//! uncapped and `1.164×` against `1.108×` capped, with the no-early-exit worst
//! case at `1.034×`/`1.591×` against `1.101×`/`1.486×`. The **slab** rows do
//! not: `1.33×` to `1.42×` here against `0.44×` to `0.59×` in the bench.
//!
//! The slab gap is sampling, and it is a **build** difference rather than a
//! disagreement about the algorithm. `crates/isomesh`'s own workspace sets
//! `lto = "thin"` and `codegen-units = 1`; `bevy_isomesh` is a separate
//! workspace and takes cargo's release defaults, so the
//! `BrushStack → BoxExact`/`Capsule` call chain does not inline across the crate
//! boundary. Measured: 6.7 ms to sample 65³ points of the slab here against
//! 3.0 ms in the bench, on the same machine. On the gyroids the field is one
//! trigonometric expression rather than a brush fold, the ratio survives, and
//! the whole-vs-marginal spread — the thing the claim is about — is unaffected
//! either way, because both arms pay the same profile.
//!
//! Neither number is adjusted to agree. The HUD reports what this build does and
//! labels M-346's figures as the thin-LTO ones, which is the only honest way to
//! show two measurements of one quantity taken under different conditions.
//!
//! # `f64`, and Marching Cubes
//!
//! M-346 was measured in `f64`, so the numbers on the HUD are reproducible only
//! in `f64`; the surface is cast to `f32` on its way into the [`Mesh`] asset and
//! nothing but the picture depends on that. Marching Cubes rather than a dual
//! method because that is what P-49 timed against, and because it is the
//! extractor a chunked pipeline would be paying for.
//!
//! # What is on screen
//!
//! - **Teal translucent tube** — the extracted isosurface of the drilled slab.
//!   It is the boundary between rock and air, drawn see-through so the ball
//!   inside it is visible. The rock itself has no surface to draw: every sample
//!   in the domain is solid except the channel, so `G`'s domain box is the block.
//! - **Grey surface** — the same isosurface on a gyroid, opaque, because there
//!   you are looking at rock rather than into a hole.
//! - **The ball** — a fixed body radius of 4 cells. On the slab it is red and
//!   parked against the `−x` face while the aperture is below it, and green and
//!   shuttling through once the aperture crosses; its width against the tube's is
//!   the same comparison the two numbers make. On a gyroid it sits at the
//!   `−x`/`+x` bottleneck witness, which is *inside* the structure, so it is
//!   mostly occluded — correctly, since that is where it would actually be, and
//!   the two circles are what you read instead.
//! - **Cyan circle** — the aperture, drawn at the sample that set it, in the
//!   plane across the `−x`/`+x` passage. The red-or-green circle inside or around
//!   it is the body radius, and the cyan line through both is the passage. These
//!   have a negative depth bias, so on a gyroid they show through the rock and
//!   carry the comparison the ball cannot.

mod common;

use std::time::Instant;

use bevy::prelude::*;
use bevy_isomesh::MeshBuilder;
use common::{Capture, CommonPlugin, DemoDomain, DemoMesh, DemoStats, OrbitCamera, ViewFlags};
use isomesh::brush::{Brush, BrushStack, Capsule};
use isomesh::fields::{BoxExact, Gyroid, ReferenceField, capped_gyroid};
use isomesh::marching_cubes::MarchingCubes;
use isomesh::{MeshBuffer, MeshSink, RuntimeShape3, Sdf};

// ─── the registered fixture ─────────────────────────────────────────────────

/// P-49's resolution. 65 samples span 64 cells, and it is the grid every number
/// in this file's prose was measured on.
const DEFAULT_SAMPLES: u32 = 65;

/// Below this there is no interior sample to be air.
const MIN_SAMPLES: u32 = 9;

/// Above this a solve allocates more than the demo is worth.
const MAX_SAMPLES: u32 = 129;

/// The three channel radii M-346 measured exactly, in cells.
const REPLAY_RADII: [u32; 3] = [2, 4, 8];

/// The body radius, in **cells**.
///
/// In cells rather than world units on purpose: the two scenes are 4 and 14
/// units across, so a fixed world radius would be a different fraction of each,
/// while cells are the unit P-49 reports its apertures in and the unit a chunk
/// pipeline budgets in. Four cells is also the middle of M-346's three replayed
/// radii, so the gate closes below one measured fixture and opens above another.
const BODY_CELLS: f64 = 4.0;

/// Widest the channel gets, in cells. The `±y`/`±z` walls are 32 cells from the
/// axis, so the channel never reaches a face it should not.
const CHANNEL_MAX_CELLS: f64 = 12.0;

/// The sweep is quantised to this many cells.
///
/// A quarter cell is smooth enough to read as an animation and coarse enough
/// that a rebuild — three solves and an extraction — happens about a dozen times
/// a second rather than sixty. `0.25` and `0.0625` are both powers of two, so
/// every radius on the sweep is exact in `f64` and so is its value in cells.
const CHANNEL_STEP_CELLS: f64 = 0.25;

/// Seconds for one full open-and-close of the channel, when nobody is capturing.
const SWEEP_SECONDS: f32 = 9.0;

/// Share of a capture spent on the slab, then on the **capped** gyroid. The
/// remainder is the uncapped one.
///
/// The capped gyroid comes second and the uncapped one last, deliberately. The
/// uncapped field is the one whose six faces are joined by the gyroid's own
/// channels rather than by the shell around a cap, so it is the frame to end on;
/// and a field switch is the single most expensive thing this example does — one
/// sample, three solves and an extraction, about 60 ms — so putting the last
/// switch at 0.80 rather than at the end leaves the closing fifth of the clip
/// rebuilding nothing.
const SLAB_SHARE: f32 = 0.62;
const CAPPED_SHARE: f32 = 0.80;

/// Passes the ball makes through the channel over the open part of the sweep.
const BALL_PASSES: f32 = 1.5;

/// M-346's worst-case cost ratios, for the HUD to be held against.
const MEASURED_WHOLE_RATIO: f64 = 1.108;
const MEASURED_MARGINAL_RATIO: f64 = 0.391;

/// The fields, in the order the digit keys select them.
const FIELD_COUNT: usize = 3;
const FIELD_NAMES: [&str; FIELD_COUNT] = ["drilled slab", "gyroid (uncapped)", "gyroid (capped)"];

// ─── faces and pairs ────────────────────────────────────────────────────────

/// Face order: `-x, +x, -y, +y, -z, +z`. Bit `f` of a face mask is face `f`.
const FACE_NAMES: [&str; 6] = ["-x", "+x", "-y", "+y", "-z", "+z"];

/// Unordered pairs of six faces.
const PAIRS: usize = 15;

/// All 15 pair bits set.
const ALL_PAIRS: u16 = 0x7FFF;

/// `(-x, +x)` is pair zero in the `i < j` enumeration, and it is the headline:
/// the pair the drilled channel connects and the one the ball travels along.
const HEADLINE: usize = 0;

/// The 15 unordered face pairs, in the `i < j` order the matrix uses.
fn pair_list() -> [(usize, usize); PAIRS] {
    let mut out = [(0usize, 0usize); PAIRS];
    let mut k = 0;
    for i in 0..6usize {
        for j in (i + 1)..6usize {
            out[k] = (i, j);
            k += 1;
        }
    }
    out
}

/// Where `(a, b)` sits in [`pair_list`], or `None` on the diagonal.
///
/// A closed form rather than a search, because the matrix paint asks 36 times a
/// frame. It is checked against [`pair_list`] at startup rather than trusted —
/// an off-by-one here would light the wrong cell and nothing on screen would say
/// so.
fn pair_index(a: usize, b: usize) -> Option<usize> {
    let (i, j) = if a < b { (a, b) } else { (b, a) };
    if i == j || j >= 6 {
        return None;
    }
    Some(i * (11 - i) / 2 + (j - i - 1))
}

/// For each 6-bit face mask, which of the 15 pairs it already contains.
///
/// Turns "has this union just connected a pair nobody has recorded yet" into one
/// table lookup and one mask test, rather than 15 comparisons per sample.
fn pair_table(pairs: &[(usize, usize); PAIRS]) -> [u16; 64] {
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

// ─── the grid ───────────────────────────────────────────────────────────────

/// The grid a field is sampled, solved and meshed on.
///
/// One definition, so the aperture and the surface cannot disagree about which
/// lattice they looked at. `samples` counts samples, so `n` samples span `n − 1`
/// cells.
#[derive(Clone, Copy)]
struct Grid {
    origin: [f64; 3],
    cell_size: f64,
    samples: u32,
}

impl Grid {
    /// The grid a reference field's own domain asks for.
    fn over<F: ReferenceField + Sdf<Scalar = f64>>(field: &F, samples: u32) -> Self {
        let (lo, hi) = field.domain();
        Self {
            origin: lo,
            cell_size: (hi[0] - lo[0]) / f64::from(samples.saturating_sub(1).max(1)),
            samples,
        }
    }

    fn point(&self, i: usize, j: usize, k: usize) -> [f64; 3] {
        [
            self.origin[0] + self.cell_size * i as f64,
            self.origin[1] + self.cell_size * j as f64,
            self.origin[2] + self.cell_size * k as f64,
        ]
    }

    /// World-space extent, as Bevy vectors.
    fn bounds(&self) -> (Vec3, Vec3) {
        let span = self.cell_size * f64::from(self.samples.saturating_sub(1));
        let lo = Vec3::new(
            self.origin[0] as f32,
            self.origin[1] as f32,
            self.origin[2] as f32,
        );
        (lo, lo + Vec3::splat(span as f32))
    }

    /// The world position of a flat sample index.
    fn position_of(&self, index: u32) -> Vec3 {
        let n = self.samples as usize;
        let flat = index as usize;
        let p = self.point(flat % n, (flat / n) % n, flat / (n * n));
        Vec3::new(p[0] as f32, p[1] as f32, p[2] as f32)
    }
}

// ─── the aperture engine ────────────────────────────────────────────────────

/// Monotone union-find over air samples, with every buffer allocated once.
struct Aperture {
    samples: u32,
    /// One field value per sample, `i + j·n + k·n²`.
    values: Vec<f64>,
    /// Air samples as `(value, index)`, sorted descending by value then
    /// ascending by index — a total order, which is the whole determinism
    /// argument.
    order: Vec<(f64, u32)>,
    parent: Vec<u32>,
    size: Vec<u32>,
    /// Face mask, meaningful on a component root.
    mask: Vec<u8>,
    active: Vec<bool>,
    /// Aperture per pair, in world units. `None` is unreachable.
    aperture: [Option<f64>; PAIRS],
    /// The sample whose value set each aperture — the bottleneck witness.
    witness: [Option<u32>; PAIRS],
    pairs: [(usize, usize); PAIRS],
    table: [u16; 64],
    /// Air samples in the last solve, and how many were visited before the
    /// early exit fired.
    air: u64,
    visited: u64,
}

impl Aperture {
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
            aperture: [None; PAIRS],
            witness: [None; PAIRS],
            table: pair_table(&pairs),
            pairs,
            air: 0,
            visited: 0,
        }
    }

    /// Fill `values` from the field.
    ///
    /// Separate from [`Aperture::solve`] so the marginal cost of the summary can
    /// be reported apart from the sampling a mesher pays anyway.
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
    /// `early_exit` stops once all 15 pairs are known. Sound rather than a
    /// shortcut — see the module docs — and checked on every rebuild by solving
    /// both ways and comparing.
    fn solve(&mut self, early_exit: bool) {
        let n = self.samples as usize;
        let plane = n * n;

        self.aperture = [None; PAIRS];
        self.witness = [None; PAIRS];
        self.order.clear();
        self.active.fill(false);

        // Air is `value > 0.0`, the strict complement of the extractors'
        // `value < 0.0`. Strict on both sides leaves the surface itself out, so
        // an aperture is a strictly positive clearance rather than a tangency.
        for (index, value) in self.values.iter().enumerate() {
            if *value > 0.0 {
                self.order.push((*value, index as u32));
            }
        }
        // Descending by value, then ascending by index. `total_cmp` rather than
        // `partial_cmp`, so this is a total order without an `Option` to
        // unwrap; the index breaks every value tie, so no two entries compare
        // equal and the sort has nothing left to decide.
        self.order
            .sort_unstable_by(|a, b| b.0.total_cmp(&a.0).then(a.1.cmp(&b.1)));
        self.air = self.order.len() as u64;
        self.visited = 0;

        let mut recorded = 0u16;
        for step in 0..self.order.len() {
            let (value, index) = self.order[step];
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
                for bit in 0..PAIRS {
                    if fresh & (1 << bit) != 0 {
                        self.aperture[bit] = Some(value);
                        self.witness[bit] = Some(index);
                    }
                }
                recorded |= fresh;
                if early_exit && recorded == ALL_PAIRS {
                    break;
                }
            }
        }
    }

    /// The 15 entries, as `-x/+x=6.5000 -x/-y=unreachable ...`, in cells.
    ///
    /// Symmetric, so only the upper triangle is named; the diagonal is not a
    /// pair.
    fn matrix(&self, cell_size: f64) -> String {
        self.pairs
            .iter()
            .enumerate()
            .map(|(bit, (i, j))| match self.aperture[bit] {
                Some(v) => format!("{}/{}={:.4}", FACE_NAMES[*i], FACE_NAMES[*j], v / cell_size),
                None => format!("{}/{}=unreachable", FACE_NAMES[*i], FACE_NAMES[*j]),
            })
            .collect::<Vec<_>>()
            .join(" ")
    }
}

// ─── the drilled slab ───────────────────────────────────────────────────────

/// The drilled slab, named so it can be passed to a `Sized` generic.
///
/// Marching Cubes' `extract` and the sampler are both generic over a sized
/// field, so this cannot be handed round as `&dyn Sdf`.
type Slab<'a> = BrushStack<'a, BoxExact<f64>, Capsule<f64>>;

/// The slab fixture at a channel radius of `cells` cells.
///
/// The brush slice has to outlive the [`BrushStack`] that borrows it, so the
/// caller owns it and this hands the stack to a closure rather than returning
/// it.
fn with_slab<T>(grid: &Grid, cells: f64, body: impl FnOnce(&Slab<'_>) -> T) -> T {
    // The capsule runs well past both x walls, so the channel is a straight
    // cylinder rather than a capped one anywhere inside the grid.
    let brushes = [Brush::subtract(Capsule {
        a: [-4.0, 0.0, 0.0],
        b: [4.0, 0.0, 0.0],
        radius: grid.cell_size * cells,
    })];
    // Half-extents of 4 against a domain half-extent of 2: the box swallows the
    // whole grid, so every sample is solid unless the capsule carves it.
    let field = BrushStack {
        base: BoxExact::<f64> {
            center: [0.0; 3],
            half_extents: [4.0; 3],
        },
        brushes: &brushes,
    };
    body(&field)
}

/// Air samples a straight cylinder of radius `cells` must produce on an `n³`
/// grid whose axis is a sample line.
///
/// A sample is air exactly when `−capsule > 0`, i.e. strictly inside the
/// cylinder, i.e. `j² + k² < cells²` in cell units — and the cylinder spans
/// every one of the `n` planes along `x`. Integer arithmetic for the lattice
/// count, so this is a closed form and not a second approximation of the same
/// thing.
///
/// Returns `None` when the axis is not a sample line, because then the closed
/// form is about a different lattice than the one being solved.
fn cylinder_air(samples: u32, cells: f64) -> Option<u64> {
    if !samples.is_multiple_of(2) {
        // An odd sample count puts a sample at the centre of every axis.
        let mid = i64::from(samples / 2);
        let limit = cells * cells;
        let reach = cells.ceil() as i64;
        let mut count = 0u64;
        for j in -reach..=reach {
            for k in -reach..=reach {
                if j.abs() > mid || k.abs() > mid {
                    continue;
                }
                if ((j * j + k * k) as f64) < limit {
                    count += 1;
                }
            }
        }
        return Some(count * u64::from(samples));
    }
    None
}

// ─── one measurement ────────────────────────────────────────────────────────

/// Everything one rebuild produced.
struct Measured {
    /// Aperture per pair, world units.
    aperture: [Option<f64>; PAIRS],
    /// Bottleneck witness per pair, world space.
    witness: [Option<Vec3>; PAIRS],
    air: u64,
    visited: u64,
    /// Two solves of the same values agreed on all 15 entries and both counts.
    deterministic: bool,
    /// The early-exit solve and the full sweep agreed on all 15 entries.
    early_exit_sound: bool,
    sample_ms: f64,
    solve_ms: f64,
    /// The same solve with the early exit disabled — the worst case.
    solve_full_ms: f64,
    extract_ms: f64,
    vertices: usize,
    triangles: usize,
    /// The 15 entries in cells, for the log.
    matrix: String,
}

/// Sample, solve three ways, and mesh.
///
/// The reported answer is the **early-exit** solve, which is the one a pipeline
/// would run. The second solve is the determinism check and the third is the
/// no-early-exit worst case; neither is allowed to change what is reported, and
/// both comparisons are on the HUD.
fn measure<F: Sdf<Scalar = f64>>(
    engine: &mut Aperture,
    mc: &mut MarchingCubes<f64>,
    buffer: &mut MeshBuffer<f64>,
    field: &F,
    grid: &Grid,
) -> Option<Measured> {
    let shape = match RuntimeShape3::new([grid.samples; 3]) {
        Ok(shape) => shape,
        Err(error) => {
            error!("grid {}^3 rejected: {error}", grid.samples);
            return None;
        }
    };

    let started = Instant::now();
    engine.sample_grid(field, grid);
    let sample_ms = started.elapsed().as_secs_f64() * 1e3;

    let started = Instant::now();
    engine.solve(true);
    let solve_ms = started.elapsed().as_secs_f64() * 1e3;
    let reported = engine.aperture;
    let witness = engine.witness;
    let (air, visited) = (engine.air, engine.visited);

    // Determinism, checked rather than inferred from the absence of a PRNG. The
    // apertures are values copied straight out of the sample array and never
    // arithmetic results, so an exact comparison is the right one here and
    // cannot be flaky.
    engine.solve(true);
    let deterministic = engine.aperture == reported
        && engine.witness == witness
        && engine.air == air
        && engine.visited == visited;

    let started = Instant::now();
    engine.solve(false);
    let solve_full_ms = started.elapsed().as_secs_f64() * 1e3;
    let early_exit_sound = engine.aperture == reported;
    if !early_exit_sound {
        error!(
            "the early exit changed the answer: {} against {}",
            engine.matrix(grid.cell_size),
            {
                engine.aperture = reported;
                engine.matrix(grid.cell_size)
            }
        );
    }

    engine.aperture = reported;
    engine.witness = witness;
    engine.air = air;
    engine.visited = visited;

    buffer.reset();
    let started = Instant::now();
    if let Err(error) = mc.extract(field, &shape, grid.origin, grid.cell_size, buffer) {
        error!("marching cubes failed at {}^3: {error}", grid.samples);
        return None;
    }
    let extract_ms = started.elapsed().as_secs_f64() * 1e3;

    let mut positions = [None; PAIRS];
    for (slot, index) in positions.iter_mut().zip(witness.iter()) {
        *slot = index.map(|i| grid.position_of(i));
    }

    Some(Measured {
        aperture: reported,
        witness: positions,
        air,
        visited,
        deterministic,
        early_exit_sound,
        sample_ms,
        solve_ms,
        solve_full_ms,
        extract_ms,
        vertices: buffer.positions.len(),
        triangles: buffer.indices.len() / 3,
        matrix: engine.matrix(grid.cell_size),
    })
}

// ─── the startup replay of M-346's exact rows ───────────────────────────────

/// M-346's three exactly-known slab rows, replayed once at startup.
///
/// This is the acceptance evidence, and it is checked rather than quoted. A
/// wrong answer must be **loud and must not take the window down with it**: it
/// logs at `error!` and the HUD says so, because a demo a stranger runs is not
/// the place for an assertion.
#[derive(Resource)]
struct Replay {
    /// One line per replayed radius, for the HUD.
    line: String,
    /// Whether every radius came back exact with exactly one reachable pair.
    exact: bool,
}

impl Replay {
    fn run(samples: u32) -> Self {
        // `pair_index`'s closed form against the enumeration it claims to
        // index. An off-by-one here would light the wrong matrix cell and
        // nothing on screen would say so.
        let pairs = pair_list();
        let mismatched = (0..6usize)
            .flat_map(|i| (0..6usize).map(move |j| (i, j)))
            .filter(|(i, j)| i != j)
            .filter(|(i, j)| {
                pair_index(*i, *j).is_none_or(|bit| pairs[bit] != (i.min(j).to_owned(), *i.max(j)))
            })
            .count();
        if mismatched > 0 {
            error!(
                "pair_index disagrees with pair_list on {mismatched} of 30 ordered pairs, so \
                 the reachability grid is indexing the wrong entries"
            );
        }

        let mut engine = Aperture::new(samples);
        let mut mc = MarchingCubes::<f64>::new();
        let mut buffer = MeshBuffer::<f64>::new();
        let grid = Grid::over(&BoxExact::<f64>::canonical(), samples);

        let mut parts = Vec::with_capacity(REPLAY_RADII.len());
        let mut exact = true;
        for radius in REPLAY_RADII {
            let cells = f64::from(radius);
            let Some(m) = with_slab(&grid, cells, |field| {
                measure(&mut engine, &mut mc, &mut buffer, field, &grid)
            }) else {
                exact = false;
                parts.push(format!("r={radius} FAILED"));
                continue;
            };
            let Some(reported) = m.aperture[HEADLINE] else {
                exact = false;
                parts.push(format!("r={radius} UNREACHABLE"));
                error!(
                    "the drilled channel of radius {radius} cells reported -x/+x unreachable, \
                     which contradicts M-346"
                );
                continue;
            };
            let error_cells = (reported / grid.cell_size - cells).abs();
            let others = m.aperture.len() - 1;
            let false_reachable = m
                .aperture
                .iter()
                .enumerate()
                .filter(|(bit, a)| *bit != HEADLINE && a.is_some())
                .count();
            let row_exact = error_cells == 0.0 && false_reachable == 0;
            exact &= row_exact;
            parts.push(format!("r{radius} {error_cells:.6}"));
            if row_exact {
                info!(
                    "M-346 replay, slab r={radius} cells at {samples}^3: aperture {:.6} cells, \
                     error {error_cells:.6}, 1 of 15 pairs reachable, {others} unreachable, \
                     air {} samples, deterministic {}, early exit sound {}",
                    reported / grid.cell_size,
                    m.air,
                    m.deterministic,
                    m.early_exit_sound,
                );
            } else {
                error!(
                    "M-346 replay, slab r={radius} cells at {samples}^3: aperture {:.6} cells is \
                     {error_cells:.6} cells off exact and {false_reachable} other pairs came back \
                     reachable. M-346 measured error 0.000000 with zero falsely-reachable pairs, \
                     so either this fixture or that finding is wrong.",
                    reported / grid.cell_size,
                );
            }
        }

        Self {
            line: format!(
                "replay    {}  cells of error at {samples}^3{}",
                parts.join("  "),
                if exact { "" } else { "  -- SEE THE LOG" }
            ),
            exact,
        }
    }
}

// ─── state ──────────────────────────────────────────────────────────────────

/// Samples per axis, fixed at startup from `ISOMESH_SAMPLES`.
#[derive(Resource)]
struct Resolution(u32);

/// A field pinned by `ISOMESH_FIELD`, which overrides the capture's own stepping.
///
/// The harness's contract is that anything a capture depends on is reachable
/// from the environment; without this a clip of one field could only be produced
/// by holding a key down.
#[derive(Resource)]
struct Pinned(Option<usize>);

/// Where the sweep is, in seconds, when nobody is capturing.
#[derive(Resource, Default)]
struct Sweep(f32);

/// What the demo is showing this frame.
#[derive(Resource, Default)]
struct Stage {
    field: usize,
    /// Channel radius in cells, on the slab only.
    channel_cells: Option<f64>,
    /// Progress through the slab sweep, `0..1`.
    t: f32,
}

/// The answer, and what it cost.
#[derive(Resource, Default)]
struct Solved {
    field: usize,
    grid_cell_size: f64,
    domain_min: Vec3,
    domain_max: Vec3,
    channel_cells: Option<f64>,
    aperture: [Option<f64>; PAIRS],
    witness: [Option<Vec3>; PAIRS],
    air: u64,
    air_expected: Option<u64>,
    visited: u64,
    deterministic: bool,
    early_exit_sound: bool,
    sample_ms: f64,
    solve_ms: f64,
    solve_full_ms: f64,
    extract_ms: f64,
    vertices: usize,
    triangles: usize,
}

impl Solved {
    /// Body radius in world units on this grid.
    fn body_world(&self) -> f64 {
        BODY_CELLS * self.grid_cell_size
    }

    /// Whether the body fits through the headline pair.
    fn fits(&self) -> bool {
        self.aperture[HEADLINE].is_some_and(|a| a >= self.body_world())
    }

    /// Aperture of a pair in cells.
    fn cells(&self, bit: usize) -> Option<f64> {
        self.aperture[bit].map(|a| a / self.grid_cell_size)
    }
}

/// The rig, the entities and the materials.
#[derive(Resource)]
struct Demo {
    engine: Aperture,
    mc: MarchingCubes<f64>,
    buffer: MeshBuffer<f64>,
    /// The surface asset, overwritten in place rather than replaced.
    mesh: Option<Handle<Mesh>>,
    surface: Entity,
    body: Entity,
    /// Translucent, for the slab: the tube has to be see-through or it hides
    /// the ball inside it.
    channel_material: Handle<StandardMaterial>,
    /// Opaque, for the gyroids: there you are looking at rock.
    rock_material: Handle<StandardMaterial>,
    body_material: Handle<StandardMaterial>,
    /// The rebuild this state was produced for.
    last_key: Option<(usize, i64)>,
    /// The pass/fail the ball is currently painted for.
    painted_fits: Option<bool>,
}

/// The gate overlay draws in front of the translucent tube, so it needs its own
/// depth bias.
#[derive(Default, Reflect, GizmoConfigGroup)]
struct GateGizmos;

// ─── app ────────────────────────────────────────────────────────────────────

fn main() {
    let samples = common::samples_override()
        .unwrap_or(DEFAULT_SAMPLES)
        .clamp(MIN_SAMPLES, MAX_SAMPLES);
    App::new()
        .add_plugins(DefaultPlugins.set(WindowPlugin {
            primary_window: Some(Window {
                title: "isomesh - E-306 aperture gate".into(),
                ..default()
            }),
            ..default()
        }))
        .add_plugins(CommonPlugin)
        .init_gizmo_group::<GateGizmos>()
        .insert_resource(Resolution(samples))
        // Replayed here rather than in `setup`, so every system can take it as a
        // plain `Res`, and after `add_plugins` so the log subscriber is already
        // installed and a failed replay is heard.
        .insert_resource(Replay::run(samples))
        .insert_resource(Pinned(pinned_field()))
        .init_resource::<Sweep>()
        .init_resource::<Stage>()
        .init_resource::<Solved>()
        .add_systems(Startup, setup)
        // **`PreUpdate`, and that is load-bearing rather than a preference.**
        //
        // Two things this example puts on screen are written by systems it does
        // not own, and `Update` gives no ordering against either: the harness's
        // `update_hud` renders `DemoStats`, and its `capture_sequence` both
        // takes the screenshot and advances `Capture::taken`. In `Update` the
        // HUD rendered a frame-old `DemoStats` while the mesh was current, and
        // `advance` read `taken` on either side of the increment — measured, and
        // the committed smoke frames showed a HUD two sweep steps behind the
        // geometry it was describing, which for a demo whose whole claim is
        // "the number on screen is this picture" is the one defect that matters.
        //
        // Running here fixes both by construction: the sweep reads `taken`
        // before it moves, and `DemoStats` is written before the HUD is built,
        // so the mesh, the readout, the grid and the HUD in any one frame are
        // all the same solve. It also puts the rebuild on the frames *between*
        // screenshots rather than sharing one with the readback.
        //
        // After `InputSystems` so a keypress is seen in the frame it happened
        // rather than the next one.
        .add_systems(
            PreUpdate,
            (
                advance,
                rebuild,
                place_body,
                draw_gate,
                paint_matrix,
                report,
            )
                .chain()
                .after(bevy::input::InputSystems),
        )
        .run();
}

/// The field `ISOMESH_FIELD` asks for, if it asks for one.
fn pinned_field() -> Option<usize> {
    let raw = std::env::var("ISOMESH_FIELD").ok()?;
    match raw.trim().parse::<usize>() {
        Ok(index) if index < FIELD_COUNT => Some(index),
        _ => {
            error!("ISOMESH_FIELD={raw} is not one of 0..{FIELD_COUNT}");
            None
        }
    }
}

fn setup(
    mut commands: Commands,
    mut meshes: ResMut<Assets<Mesh>>,
    mut materials: ResMut<Assets<StandardMaterial>>,
    mut config: ResMut<GizmoConfigStore>,
    resolution: Res<Resolution>,
) {
    let (gate, _) = config.config_mut::<GateGizmos>();
    gate.line.width = 2.4;
    gate.depth_bias = -0.6;

    // Translucent and double-sided, or the tube hides the ball travelling
    // inside it -- and a one-sided tube seen from outside shows its own
    // interior, which reads as a rod rather than a hole.
    let channel_material = materials.add(StandardMaterial {
        base_color: Color::srgba(0.32, 0.86, 0.90, 0.30),
        perceptual_roughness: 0.4,
        alpha_mode: AlphaMode::Blend,
        double_sided: true,
        cull_mode: None,
        ..default()
    });
    // Darker than this repo's usual surface grey. A gyroid at 65^3 fills most
    // of the frame, and at the harness's default 0.72 the light HUD text sat on
    // white and was unreadable in the captured frames -- the numbers are the
    // evidence, so the rock loses.
    let rock_material = materials.add(StandardMaterial {
        base_color: Color::srgb(0.26, 0.29, 0.36),
        perceptual_roughness: 0.72,
        metallic: 0.04,
        double_sided: true,
        cull_mode: None,
        ..default()
    });
    let body_material = materials.add(StandardMaterial {
        base_color: Color::srgb(0.92, 0.26, 0.22),
        emissive: LinearRgba::new(0.30, 0.03, 0.02, 1.0),
        perceptual_roughness: 0.35,
        ..default()
    });

    // `Mesh3d::default()` names no asset, so nothing is drawn and nothing is
    // uploaded until the first rebuild. An empty mesh would be worse than
    // nothing: `MeshAllocator` skips a zero-byte vertex buffer and then copies
    // into it anyway, once per frame, in red.
    let surface = commands
        .spawn((
            Mesh3d::default(),
            MeshMaterial3d(channel_material.clone()),
            DemoMesh,
        ))
        .id();

    // A unit sphere scaled to the body radius, so the radius is a transform
    // rather than a mesh rebuild. Not a `DemoMesh`: it is the probe, not the
    // subject, and wireframing it would bury the tube.
    let body = commands
        .spawn((
            Mesh3d(meshes.add(Sphere::new(1.0).mesh().uv(28, 14))),
            MeshMaterial3d(body_material.clone()),
            // Zero scale until the first solve gives it a radius, so a failed
            // rebuild leaves nothing on screen rather than a unit sphere at the
            // origin claiming to be a player.
            Transform::from_scale(Vec3::ZERO),
        ))
        .id();

    // Spawned rather than assumed: `draw_domain` queries for it, and without one
    // the `G` toggle silently does nothing. Filled in by the first rebuild.
    commands.spawn(DemoDomain {
        min: Vec3::splat(-1.0),
        max: Vec3::splat(1.0),
    });

    spawn_matrix(&mut commands);
    spawn_readout(&mut commands);

    commands.insert_resource(Demo {
        engine: Aperture::new(resolution.0),
        mc: MarchingCubes::<f64>::new(),
        buffer: MeshBuffer::<f64>::new(),
        mesh: None,
        surface,
        body,
        channel_material,
        rock_material,
        body_material,
        last_key: None,
        painted_fits: None,
    });
}

/// Frames a capture runs for.
///
/// Read from the environment rather than from [`Capture`], which keeps its
/// length private, because pacing the sweep off the capture is what stops a
/// six-frame smoke test and a ninety-frame clip from both being a still.
fn capture_frames() -> u32 {
    std::env::var("ISOMESH_CAPTURE_FRAMES")
        .ok()
        .and_then(|v| v.parse().ok())
        // The harness's own default.
        .unwrap_or(60)
        .max(1)
}

/// Decide the field and the channel radius for this frame.
///
/// Under capture both come off the captured-frame counter, so a clip of any
/// length shows the gate open and then both gyroids. Interactively the digits
/// pick the field and the channel sweeps on a loop, so the gate opens and closes
/// again without anybody pressing anything.
fn advance(
    time: Res<Time>,
    keys: Res<ButtonInput<KeyCode>>,
    capture: Res<Capture>,
    flags: Res<ViewFlags>,
    pinned: Res<Pinned>,
    mut sweep: ResMut<Sweep>,
    mut stage: ResMut<Stage>,
) {
    if keys.just_pressed(KeyCode::KeyX) {
        sweep.0 = 0.0;
    }

    let (field, t) = if capture.is_active() {
        let total = capture_frames();
        let phase = f64::from(capture.taken) / f64::from(total);
        let phase = phase as f32;
        if phase < SLAB_SHARE {
            (0, (phase / SLAB_SHARE).clamp(0.0, 1.0))
        } else if phase < CAPPED_SHARE {
            (2, 1.0)
        } else {
            (1, 1.0)
        }
    } else {
        if !flags.paused {
            sweep.0 += time.delta_secs();
        }
        let cycle = (sweep.0 / SWEEP_SECONDS).fract();
        (flags.field.min(FIELD_COUNT - 1), cycle)
    };

    stage.field = pinned.0.unwrap_or(field);
    stage.t = t;
    stage.channel_cells = (stage.field == 0).then(|| {
        // Quantised so a rebuild happens about a dozen times a second rather
        // than sixty. Both the step and the cell size are powers of two, so
        // every radius on the sweep is exact and so is its value in cells.
        let raw = CHANNEL_MAX_CELLS * f64::from(t);
        (raw / CHANNEL_STEP_CELLS).round() * CHANNEL_STEP_CELLS
    });
}

/// Sample, solve, mesh, and frame the camera — only when the answer would change.
#[allow(clippy::too_many_arguments)]
fn rebuild(
    stage: Res<Stage>,
    resolution: Res<Resolution>,
    flags: Res<ViewFlags>,
    mut demo: ResMut<Demo>,
    mut solved: ResMut<Solved>,
    mut meshes: ResMut<Assets<Mesh>>,
    mut commands: Commands,
    mut domain: Query<&mut DemoDomain>,
    mut camera: Query<&mut OrbitCamera>,
) {
    // Quarter-cell steps, as an integer, so the key is exact.
    let key = (
        stage.field,
        (stage.channel_cells.unwrap_or(0.0) / CHANNEL_STEP_CELLS).round() as i64,
    );
    if demo.last_key == Some(key) && !flags.remesh_requested {
        return;
    }
    let field_changed = demo.last_key.map(|(f, _)| f) != Some(stage.field);
    demo.last_key = Some(key);

    let samples = resolution.0;
    let demo = demo.as_mut();
    let (grid, measured) = match stage.field {
        0 => {
            let grid = Grid::over(&BoxExact::<f64>::canonical(), samples);
            let cells = stage.channel_cells.unwrap_or(0.0);
            let m = with_slab(&grid, cells, |field| {
                measure(
                    &mut demo.engine,
                    &mut demo.mc,
                    &mut demo.buffer,
                    field,
                    &grid,
                )
            });
            (grid, m)
        }
        1 => {
            // Both gyroids are sampled on the capped field's domain, as P-49
            // measured them, so the two rows are directly comparable.
            let grid = Grid::over(&capped_gyroid::<f64>(), samples);
            let field = Gyroid::<f64>::canonical();
            let m = measure(
                &mut demo.engine,
                &mut demo.mc,
                &mut demo.buffer,
                &field,
                &grid,
            );
            (grid, m)
        }
        _ => {
            let field = capped_gyroid::<f64>();
            let grid = Grid::over(&field, samples);
            let m = measure(
                &mut demo.engine,
                &mut demo.mc,
                &mut demo.buffer,
                &field,
                &grid,
            );
            (grid, m)
        }
    };
    let Some(measured) = measured else {
        return;
    };

    let (lo, hi) = grid.bounds();
    for mut d in &mut domain {
        d.min = lo;
        d.max = hi;
    }

    let air_expected = (stage.field == 0)
        .then(|| stage.channel_cells.and_then(|c| cylinder_air(samples, c)))
        .flatten();

    *solved = Solved {
        field: stage.field,
        grid_cell_size: grid.cell_size,
        domain_min: lo,
        domain_max: hi,
        channel_cells: stage.channel_cells,
        aperture: measured.aperture,
        witness: measured.witness,
        air: measured.air,
        air_expected,
        visited: measured.visited,
        deterministic: measured.deterministic,
        early_exit_sound: measured.early_exit_sound,
        sample_ms: measured.sample_ms,
        solve_ms: measured.solve_ms,
        solve_full_ms: measured.solve_full_ms,
        extract_ms: measured.extract_ms,
        vertices: measured.vertices,
        triangles: measured.triangles,
    };

    if let Some(expected) = air_expected
        && expected != measured.air
    {
        error!(
            "air census disagrees with the closed form on the slab: the union-find found {} air \
             samples where the integer lattice count for radius {:?} cells is {expected}",
            measured.air, stage.channel_cells,
        );
    }

    // The surface. A field with no crossing gets no asset at all rather than an
    // empty one -- see `setup`.
    if measured.triangles == 0 {
        if demo.mesh.take().is_some() {
            commands.entity(demo.surface).insert(Mesh3d::default());
        }
    } else {
        let mesh = to_mesh(&demo.buffer);
        match &demo.mesh {
            // Overwritten in place rather than added and dropped: this runs a
            // dozen times a second and Bevy's mesh slab allocator does not
            // enjoy that many handle churns.
            Some(handle) => {
                if let Some(mut slot) = meshes.get_mut(handle) {
                    *slot = mesh;
                }
            }
            None => {
                let handle = meshes.add(mesh);
                commands.entity(demo.surface).insert(Mesh3d(handle.clone()));
                demo.mesh = Some(handle);
            }
        }
    }

    if field_changed {
        let material = if stage.field == 0 {
            demo.channel_material.clone()
        } else {
            demo.rock_material.clone()
        };
        commands
            .entity(demo.surface)
            .insert(MeshMaterial3d(material));

        // Framed on a field change only, so a mouse drag and a scroll survive
        // and `ISOMESH_SPIN` still adds yaw on top.
        let extent = hi.x - lo.x;
        for mut orbit in &mut camera {
            orbit.focus = (lo + hi) * 0.5;
            if stage.field == 0 {
                // Nearly down z, so the channel runs across the frame and the
                // ball's travel is a horizontal move. Off-axis enough that the
                // tube reads as a tube.
                orbit.radius = extent * 1.42;
                orbit.yaw = 1.24;
                orbit.pitch = 0.30;
            } else {
                // Outside the structure. A gyroid's half-diagonal here is
                // 12.1 units and E-110 found E-109's committed screenshot was a
                // picture of an inner wall for exactly this reason; at 1.18
                // extents the surface filled the frame and buried the HUD.
                orbit.radius = extent * 1.58;
                orbit.yaw = 0.72;
                orbit.pitch = 0.34;
            }
        }
    }

    info!(
        "{} at {samples}^3{}: {} of 15 pairs reachable, air {} samples, {} visited before the \
         early exit; sample {:.2} ms + solve {:.2} ms = {:.2} ms whole against {:.2} ms extract \
         = {:.3}x (marginal {:.3}x, no early exit {:.3}x); deterministic {}, early exit sound {}",
        FIELD_NAMES[stage.field],
        stage
            .channel_cells
            .map_or_else(String::new, |c| format!(", channel r {c:.2} cells")),
        measured.aperture.iter().filter(|a| a.is_some()).count(),
        measured.air,
        measured.visited,
        measured.sample_ms,
        measured.solve_ms,
        measured.sample_ms + measured.solve_ms,
        measured.extract_ms,
        (measured.sample_ms + measured.solve_ms) / measured.extract_ms,
        measured.solve_ms / measured.extract_ms,
        (measured.sample_ms + measured.solve_full_ms) / measured.extract_ms,
        measured.deterministic,
        measured.early_exit_sound,
    );
    info!("    {}", measured.matrix);
}

/// The `f64` extraction as a Bevy mesh.
///
/// Cast rather than re-extracted in `f32`: the numbers on the HUD are M-346's
/// and they are `f64` numbers, so the mesh the picture is drawn from has to be
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

/// Size the ball, colour it, and move it.
///
/// On the slab the channel is a straight line of air along `x`, so travelling
/// the ball down it is literally correct rather than a suggestion. On a gyroid
/// there is no straight path, so the ball parks at the `−x`/`+x` bottleneck
/// witness — the tightest point of the widest route — where "this fits" is a
/// statement about a real place.
fn place_body(
    solved: Res<Solved>,
    stage: Res<Stage>,
    mut demo: ResMut<Demo>,
    mut materials: ResMut<Assets<StandardMaterial>>,
    mut transforms: Query<&mut Transform>,
) {
    let body = solved.body_world() as f32;
    if body <= 0.0 {
        return;
    }
    let fits = solved.fits();

    if demo.painted_fits != Some(fits)
        && let Some(mut material) = materials.get_mut(&demo.body_material)
    {
        material.base_color = if fits {
            Color::srgb(0.22, 0.88, 0.38)
        } else {
            Color::srgb(0.94, 0.26, 0.22)
        };
        material.emissive = if fits {
            LinearRgba::new(0.03, 0.34, 0.08, 1.0)
        } else {
            LinearRgba::new(0.34, 0.03, 0.02, 1.0)
        };
        demo.painted_fits = Some(fits);
    }

    let centre = (solved.domain_min + solved.domain_max) * 0.5;
    let translation = if solved.field == 0 {
        let entry = Vec3::new(solved.domain_min.x - body, centre.y, centre.z);
        if fits {
            let exit = Vec3::new(solved.domain_max.x + body, centre.y, centre.z);
            // Ping-pong, so the ball keeps shuttling through for as long as the
            // gate is open instead of parking on the far side for half the clip.
            // Derived from the sweep's own progress rather than from a timer, so
            // a capture is reproducible frame for frame.
            let opened = (BODY_CELLS / CHANNEL_MAX_CELLS) as f32;
            let after = ((stage.t - opened) / (1.0 - opened).max(1e-3)).clamp(0.0, 1.0);
            let u = (after * BALL_PASSES * 2.0) % 2.0;
            entry.lerp(exit, if u <= 1.0 { u } else { 2.0 - u })
        } else {
            entry
        }
    } else {
        solved.witness[HEADLINE].unwrap_or(centre)
    };

    if let Ok(mut transform) = transforms.get_mut(demo.body) {
        transform.translation = translation;
        transform.scale = Vec3::splat(body);
    }
}

/// The two circles that are the whole comparison, in the plane across the
/// `−x`/`+x` passage.
///
/// Cyan is the aperture, drawn at the sample that set it. The red-or-green one
/// is the body radius. Which circle is inside which is the answer, and it is the
/// same answer the two numbers give.
fn draw_gate(solved: Res<Solved>, mut gizmos: Gizmos<GateGizmos>) {
    let Some(aperture) = solved.aperture[HEADLINE] else {
        return;
    };
    let Some(witness) = solved.witness[HEADLINE] else {
        return;
    };
    // A circle is drawn in the isometry's xy plane, so this rotation puts its
    // normal on world +x -- the axis the headline pair is about.
    let facing = Quat::from_rotation_y(std::f32::consts::FRAC_PI_2);
    let at = Isometry3d::new(witness, facing);
    let body = solved.body_world() as f32;
    let fits = solved.fits();

    gizmos
        .circle(at, aperture as f32, Color::srgb(0.25, 0.88, 0.95))
        .resolution(64);
    gizmos
        .circle(
            at,
            body,
            if fits {
                Color::srgb(0.22, 0.88, 0.38)
            } else {
                Color::srgb(0.94, 0.26, 0.22)
            },
        )
        .resolution(64);
    // The passage itself, so the plane above reads as a cross-section of a
    // route rather than as a free-floating disc.
    gizmos.line(
        Vec3::new(solved.domain_min.x, witness.y, witness.z),
        Vec3::new(solved.domain_max.x, witness.y, witness.z),
        Color::srgb(0.25, 0.88, 0.95).with_alpha(0.5),
    );
}

// ─── the reachability grid ──────────────────────────────────────────────────

/// Side of one matrix cell and the stride between them, in logical pixels.
const CELL_PX: f32 = 15.0;
const CELL_STRIDE: f32 = 17.0;

/// Where the grid's bottom-right corner sits.
const MATRIX_RIGHT: f32 = 16.0;
const MATRIX_BOTTOM: f32 = 22.0;

/// Width reserved for the row labels, left of the grid.
const LABEL_W: f32 = 20.0;

/// One cell of the 6x6.
#[derive(Component)]
struct MatrixCell {
    row: usize,
    col: usize,
}

/// Everything the grid and the readout are made of, so `nohud` hides all of it.
#[derive(Component)]
struct HudPanel;

/// The big number.
#[derive(Component)]
struct Readout;

/// The 6x6 reachability grid, bottom right.
///
/// Root-level absolutely-positioned nodes anchored from the right edge rather
/// than a flex grid, so a cell's position is arithmetic rather than a layout
/// outcome, and [`GlobalZIndex`] rather than spawn order, so the stacking is
/// stated.
fn spawn_matrix(commands: &mut Commands) {
    let grid_w = 6.0 * CELL_STRIDE - (CELL_STRIDE - CELL_PX);
    commands.spawn((
        Node {
            position_type: PositionType::Absolute,
            right: Val::Px(MATRIX_RIGHT - 7.0),
            bottom: Val::Px(MATRIX_BOTTOM - 7.0),
            width: Val::Px(grid_w + LABEL_W + 14.0),
            height: Val::Px(grid_w + 14.0),
            ..default()
        },
        BackgroundColor(Color::srgba(0.04, 0.05, 0.08, 0.86)),
        GlobalZIndex(1),
        HudPanel,
    ));
    commands.spawn((
        Text::new("all 15 face pairs\ndark none  amber tight\ngreen the body fits"),
        TextFont {
            font_size: FontSize::Px(11.0),
            ..default()
        },
        TextColor(Color::srgb(0.80, 0.84, 0.90)),
        Node {
            position_type: PositionType::Absolute,
            right: Val::Px(MATRIX_RIGHT - 7.0),
            bottom: Val::Px(MATRIX_BOTTOM + grid_w + 22.0),
            ..default()
        },
        GlobalZIndex(2),
        HudPanel,
    ));

    for (face, name) in FACE_NAMES.iter().enumerate() {
        // Column labels above the grid, row labels left of it. Same order in
        // both directions, because the matrix is symmetric.
        commands.spawn((
            Text::new(*name),
            TextFont {
                font_size: FontSize::Px(10.0),
                ..default()
            },
            TextColor(Color::srgb(0.72, 0.76, 0.84)),
            Node {
                position_type: PositionType::Absolute,
                right: Val::Px(MATRIX_RIGHT + (5 - face) as f32 * CELL_STRIDE),
                bottom: Val::Px(MATRIX_BOTTOM + grid_w + 4.0),
                width: Val::Px(CELL_PX),
                ..default()
            },
            GlobalZIndex(3),
            HudPanel,
        ));
        commands.spawn((
            Text::new(*name),
            TextFont {
                font_size: FontSize::Px(10.0),
                ..default()
            },
            TextColor(Color::srgb(0.72, 0.76, 0.84)),
            Node {
                position_type: PositionType::Absolute,
                right: Val::Px(MATRIX_RIGHT + grid_w + 2.0),
                bottom: Val::Px(MATRIX_BOTTOM + (5 - face) as f32 * CELL_STRIDE + 2.0),
                ..default()
            },
            GlobalZIndex(3),
            HudPanel,
        ));
    }

    for row in 0..6usize {
        for col in 0..6usize {
            let right = MATRIX_RIGHT + (5 - col) as f32 * CELL_STRIDE;
            let bottom = MATRIX_BOTTOM + (5 - row) as f32 * CELL_STRIDE;
            // A white backing square behind the headline pair, one pixel proud
            // on every side, so `-x/+x` is identifiable without a border
            // component.
            if pair_index(row, col) == Some(HEADLINE) {
                commands.spawn((
                    Node {
                        position_type: PositionType::Absolute,
                        right: Val::Px(right - 2.0),
                        bottom: Val::Px(bottom - 2.0),
                        width: Val::Px(CELL_PX + 4.0),
                        height: Val::Px(CELL_PX + 4.0),
                        ..default()
                    },
                    BackgroundColor(Color::srgb(0.95, 0.96, 1.0)),
                    GlobalZIndex(2),
                    HudPanel,
                ));
            }
            commands.spawn((
                Node {
                    position_type: PositionType::Absolute,
                    right: Val::Px(right),
                    bottom: Val::Px(bottom),
                    width: Val::Px(CELL_PX),
                    height: Val::Px(CELL_PX),
                    ..default()
                },
                BackgroundColor(Color::srgb(0.12, 0.13, 0.17)),
                GlobalZIndex(3),
                MatrixCell { row, col },
                HudPanel,
            ));
        }
    }
}

/// The headline number, centred along the bottom edge where a GIF viewer looks.
fn spawn_readout(commands: &mut Commands) {
    commands
        .spawn((
            Node {
                position_type: PositionType::Absolute,
                left: Val::Px(0.0),
                right: Val::Px(0.0),
                bottom: Val::Px(16.0),
                justify_content: JustifyContent::Center,
                ..default()
            },
            GlobalZIndex(4),
            HudPanel,
        ))
        .with_children(|parent| {
            parent.spawn((
                Text::new(""),
                TextFont {
                    font_size: FontSize::Px(18.0),
                    ..default()
                },
                // `NoWrap`, and it is not cosmetic. In a centring flex row the
                // text measure is handed the container's whole width while the
                // node's own height is resolved before the wrap, so a soft wrap
                // pushed the third line -- the one that says FITS -- off the
                // bottom of the frame. Measured on a 1280x720 capture.
                TextLayout {
                    linebreak: bevy::text::LineBreak::NoWrap,
                    ..default()
                },
                TextColor(Color::srgb(0.94, 0.26, 0.22)),
                BackgroundColor(Color::srgba(0.03, 0.04, 0.07, 0.86)),
                Node {
                    padding: UiRect::axes(Val::Px(10.0), Val::Px(6.0)),
                    ..default()
                },
                Readout,
            ));
        });
}

/// Colour the 36 cells and write the readout.
///
/// **Three colours rather than two, and that is the demo's thesis on screen.**
/// A connectivity query has only "connected" and "not"; amber is the state it
/// cannot name — a pair joined by air that the body still does not fit through.
fn paint_matrix(
    solved: Res<Solved>,
    flags: Res<ViewFlags>,
    mut cells: Query<(&MatrixCell, &mut BackgroundColor)>,
    mut readout: Query<(&mut Text, &mut TextColor), With<Readout>>,
    mut panels: Query<&mut Visibility, With<HudPanel>>,
) {
    let wanted = if flags.hud {
        Visibility::Inherited
    } else {
        Visibility::Hidden
    };
    for mut visibility in &mut panels {
        if *visibility != wanted {
            *visibility = wanted;
        }
    }
    if !flags.hud {
        return;
    }

    const DIAGONAL: Color = Color::srgb(0.09, 0.10, 0.13);
    const NONE: Color = Color::srgb(0.24, 0.13, 0.15);
    const TIGHT: Color = Color::srgb(0.96, 0.68, 0.16);
    const FITS: Color = Color::srgb(0.22, 0.88, 0.38);

    let body = solved.body_world();
    for (cell, mut colour) in &mut cells {
        let wanted = match pair_index(cell.row, cell.col) {
            None => DIAGONAL,
            Some(bit) => match solved.aperture[bit] {
                None => NONE,
                Some(a) if a >= body => FITS,
                Some(_) => TIGHT,
            },
        };
        if colour.0 != wanted {
            colour.0 = wanted;
        }
    }

    let fits = solved.fits();
    // Kept inside 30 characters on purpose. The readout is centred along the
    // bottom and the reachability grid is anchored bottom-right; at 18px a
    // wider string reaches under the grid on a 640-wide frame.
    let text = match (solved.aperture[HEADLINE], solved.cells(HEADLINE)) {
        (Some(world), Some(cells)) => format!(
            "aperture {world:.6}  {cells:.2} cells\n\
             body     {body:.6}  {BODY_CELLS:.2} cells\n{}",
            if fits { "FITS" } else { "TOO TIGHT" }
        ),
        _ => format!(
            "aperture unreachable\n\
             body     {body:.6}  {BODY_CELLS:.2} cells\n\
             NO AIR PATH"
        ),
    };
    let colour = if fits {
        Color::srgb(0.30, 0.94, 0.46)
    } else {
        Color::srgb(0.98, 0.40, 0.34)
    };
    for (mut target, mut text_colour) in &mut readout {
        if target.0 != text {
            target.0.clone_from(&text);
        }
        if text_colour.0 != colour {
            text_colour.0 = colour;
        }
    }
}

/// The HUD. The numbers are the demo.
fn report(
    solved: Res<Solved>,
    replay: Res<Replay>,
    resolution: Res<Resolution>,
    mut stats: ResMut<DemoStats>,
) {
    let samples = resolution.0;
    let body = solved.body_world();
    let whole = solved.sample_ms + solved.solve_ms;
    let full = solved.sample_ms + solved.solve_full_ms;
    let ratio = |ms: f64| {
        if solved.extract_ms > 0.0 {
            ms / solved.extract_ms
        } else {
            0.0
        }
    };
    let reachable: Vec<f64> = {
        let mut v: Vec<f64> = (0..PAIRS).filter_map(|bit| solved.cells(bit)).collect();
        v.sort_by(f64::total_cmp);
        v
    };

    // Every line below is kept inside 76 characters. At the harness's 13px font
    // that is about 600 logical pixels, so the HUD does not wrap at the 640x360
    // a smoke capture uses -- and a wrapped line in a GIF reads as a bug.
    stats.title = format!(
        "E-306  aperture gate - {}  {samples}^3  [1-3] field  [X] restart",
        FIELD_NAMES[solved.field.min(FIELD_COUNT - 1)],
    );
    stats.vertices = solved.vertices;
    stats.triangles = solved.triangles;
    stats.extract_ms = solved.extract_ms;

    stats.extra = vec![
        format!(
            "field     cell {:.6} world   {}",
            solved.grid_cell_size,
            solved.channel_cells.map_or_else(
                || String::from("the air network itself is the channel"),
                |c| format!("Capsule r {c:.2} cells along x"),
            ),
        ),
        String::new(),
        match (solved.aperture[HEADLINE], solved.cells(HEADLINE)) {
            (Some(world), Some(cells)) => format!(
                "-x/+x     aperture {world:.6} = {cells:.4} cells   body {body:.6}   {}",
                if solved.fits() { "FITS" } else { "TOO TIGHT" },
            ),
            _ => String::from("-x/+x     unreachable - no air path joins the x faces"),
        },
        format!(
            "15 pairs  {} reachable   apertures {} cells   body {BODY_CELLS:.2} cells",
            reachable.len(),
            match (reachable.first(), reachable.last()) {
                (Some(lo), Some(hi)) => format!("{lo:.4} to {hi:.4}"),
                _ => String::from("n/a"),
            },
        ),
        match (solved.channel_cells, solved.cells(HEADLINE)) {
            (Some(r), Some(cells)) => format!(
                "exact     channel {r:.4} -> aperture {cells:.4} cells, error {:.6}",
                (cells - r).abs()
            ),
            _ => String::from("exact     n/a - no drilled channel, so no closed form"),
        },
        format!(
            "air       {} samples, {} visited ({:.1}%)   census {}",
            solved.air,
            solved.visited,
            if solved.air == 0 {
                0.0
            } else {
                100.0 * solved.visited as f64 / solved.air as f64
            },
            match solved.air_expected {
                Some(e) if e == solved.air => format!("MATCHES {e}"),
                Some(e) => format!("DISAGREES with {e}"),
                None => String::from("n/a"),
            },
        ),
        format!(
            "checks    determinism {}   early exit {}   replay {}",
            yes_no(solved.deterministic),
            yes_no(solved.early_exit_sound),
            if replay.exact { "exact" } else { "FAILED" },
        ),
        String::new(),
        format!(
            "cost      {:.2} sample + {:.2} solve = {whole:.2} ms / {:.2} ms mesh = {:.3}x",
            solved.sample_ms,
            solved.solve_ms,
            solved.extract_ms,
            ratio(whole),
        ),
        format!(
            "          marginal {:.3}x   no early exit {:.3}x   full sweep {full:.2} ms",
            ratio(solved.solve_ms),
            ratio(full),
        ),
        format!(
            "          M-346 worst case, thin-LTO build: \
             {MEASURED_WHOLE_RATIO:.3}x / {MEASURED_MARGINAL_RATIO:.3}x marginal"
        ),
        replay.line.clone(),
    ];
}

fn yes_no(value: bool) -> &'static str {
    if value { "ok" } else { "FAILED" }
}
