//! E-307 — almost none of a voxel grid is surface, and the bitmap is how you
//! find that out for free.
//!
//! ```bash
//! cd bevy_isomesh && cargo run --example active_cells --release
//! ```
//!
//! **Always `--release`.** This samples the field once for itself and again for
//! the extractor, and it times two predicates five times each.
//!
//! `1`-`5` switch field, `[` and `]` step the resolution ladder, `,` and `.`
//! move the row cursor, `H` hides the extracted surface so the shell can be seen
//! naked. The rest are the shared keys — `G` toggles the grid box, `Space`
//! freezes the cursor, `F12` screenshot.
//!
//! Under `ISOMESH_CAPTURE` it needs no keyboard: the resolution walks the ladder
//! and the cursor scans rows on every captured frame, so the clip is the ratio
//! falling rather than a still. `ISOMESH_FIELD` pins the field, `ISOMESH_SAMPLES`
//! pins the resolution and takes the ladder out of the capture, `ISOMESH_SURFACE=off`
//! starts with the surface hidden.
//!
//! ```bash
//! # 90 frames is 30 per rung, so the clip is exactly one pass up the ladder.
//! # Measured at 1280x720 scaled to 900px: 0.81 MB.
//! ISOMESH_CAPTURE_FRAMES=90 ISOMESH_CAPTURE_EVERY=2 \
//!   scripts/record_gif.sh active_cells docs/gifs/e307.gif
//! ```
//!
//! The camera is still, so the only things moving are the strip, the highlighted
//! row and the numbers — which is exactly what a GIF compresses well, hence the
//! 0.81 MB. Adding parallax costs most of the band: the same recipe with
//! `ISOMESH_SPIN=0.003` measures **3.68 MB**, and **2.83 MB** with `COLORS=64`.
//! Both are inside the 0.7-4.8 MB the committed clips sit within, so the spin is
//! affordable if the cell lattice needs to read as a lattice; it is not needed
//! for the argument, which is a number and a bitmap.
//!
//! **Record it at 1280x720 or wider.** The HUD is twenty-three lines in the
//! upper left and the bit strip is 490 by 120 in the lower left, so the layout
//! wants about **840 by 560** — measured: at 836x671 nothing overlaps and every
//! line is legible. Below that the strip climbs into the HUD. A `640x360`
//! capture still runs and is a fine smoke test, it is just not a clip.
//!
//! Demonstrates **M-337 / P-40**.
//!
//! # What is on screen
//!
//! - **Grey translucent** — the Surface Nets surface. `H` removes it.
//! - **Cyan blocks** — the **active cells**: the cells whose eight corner signs
//!   are not all equal, i.e. the ones the surface actually passes through. This
//!   is the set the mechanism computes.
//! - **Grey wire box** — the whole sampled domain, for scale. The cyan shell is
//!   the fraction of it the extractor has to do any work in.
//! - **White bar** — the run of cells one `u64` word covers, up to sixty-four of
//!   them, all decided in a single pass of about twenty word operations.
//!   **Yellow cages** are the active ones among them.
//! - **The strip, lower left** — that word, drawn. Four rows of `inside` bits for
//!   the four sample rows bounding this cell row, then the three fused words
//!   `any`, `all` and `active` underneath. Sixty-four decisions, three
//!   instructions.
//! - **The big number, upper right** — `active cells / total cells`. It is the
//!   biggest thing on screen because it is the entire justification for the
//!   mechanism: as resolution climbs the absolute count rises and the *fraction*
//!   falls, so the share of the grid that is worth gathering eight corners for
//!   goes to zero.
//!
//! # The ladder is 33 / 64 / 128 samples, and two of the three are M-337's own grids
//!
//! Cells per axis are 32 / 63 / 127, so the strip means something different at
//! each rung and all three are worth seeing:
//!
//! | samples | cells/axis | `u64` words per sample row | the strip is |
//! |---:|---:|---:|---|
//! | 33 | 32 | 1 | **half a word** — bits 32-63 are not cells, and are dimmed |
//! | 64 | 63 | 1 | **one row of cells**, bar bit 63, which `cell_mask` drops |
//! | 128 | 127 | 2 | half a row, and **the word boundary is live** |
//!
//! At 128 samples bit 63 of word 0 has to take its `+x` corner from bit 0 of
//! word 1 — the one real hazard in the mechanism, and the reason M-337 checked
//! its hashes at 256³ rather than at the golden suite's 33³, where every row is
//! one word and the bug cannot show.
//!
//! **`64k + 1` samples is deliberately not on the ladder, and that is a finding
//! rather than a preference.** 65 samples looks like the nicest rung there is —
//! 64 cells per axis, one whole word of cells, the cleanest possible reading of
//! "sixty-four decide at once". But `bit_row` comes from *samples*
//! (`samples.div_ceil(64)`, which is 2) while only `cells.div_ceil(64)` words —
//! one — carry any cells, and `place_vertices` loops `w` over `bit_row`. So a
//! 65-sample grid fuses a whole extra word per cell row and then has
//! [`cell_mask`] throw all 64 of its bits away.
//!
//! Measured through this file's own harness, `ISOMESH_SAMPLES` pinning each grid,
//! `sphere`, median of five within a run and the range over three to five runs:
//!
//! | samples | words/row | words with cells | packed | ratio |
//! |---:|---:|---:|---:|---:|
//! | 64 | 1 | 1 | 0.464-0.499 ns/cell | **5.33-5.76×** |
//! | 65 | 2 | 1 | 0.603-0.655 ns/cell | **4.04-4.38×** |
//! | 128 | 2 | 2 | 0.451-0.470 ns/cell | **5.66-5.76×** |
//! | 129 | 3 | 2 | 0.599-0.637 ns/cell | **4.13-4.56×** |
//!
//! **The packed ranges do not overlap across the `+1`**, so this is the word loop
//! and not run-to-run noise: one wasted word in two costs about 30% of the packed
//! arm and one in three about 32%, on a field whose active fraction moved by 0.08
//! of a percentage point between the two grids. It is **not** a correctness bug —
//! the mask is right and the mesh is unchanged. It is about a fifth of the stage
//! left on the floor at those sizes, and the fix is one expression in `dual.rs`,
//! which is not this file's to make.
//!
//! # The numbers, against M-337
//!
//! M-337 measured `sphere` at 64 / 128 / 256 samples per axis: active
//! **4,730 / 19,010 / 76,778** cells = **1.89% / 0.93% / 0.46%**, and the mesh
//! byte-identical on 12 of 12 grids. Two rungs here are two of those grids, so
//! the counts are comparable digit for digit rather than approximately.
//!
//! Measured live in this window, `f64`, median of five within a run, range over
//! five runs, while Bevy is rendering:
//!
//! | samples | active cells | fraction | M-337 | scalar ns/cell | packed ns/cell | ratio |
//! |---:|---:|---:|---:|---:|---:|---:|
//! | 33 | 1,160 of 32,768 | 3.540% | — | 2.64-2.66 | 0.614-0.627 | **4.23-4.33×** |
//! | 64 | **4,730** of 250,047 | **1.892%** | 4,730 / 1.89% | 2.65-2.72 | 0.464-0.499 | **5.33-5.76×** |
//! | 128 | **19,010** of 2,048,383 | **0.928%** | 19,010 / 0.93% | 2.59-2.66 | 0.451-0.470 | **5.66-5.76×** |
//!
//! **The active counts are M-337's, exactly, to the last digit**, and the
//! fraction halves per doubling of resolution while the absolute count
//! quadruples — which is the whole argument, on screen, in one number.
//!
//! **The stage ratio agrees with M-337's band, and the band is the claim — no
//! single decimal of it is.** The FINDINGS M-337 table gives `sphere` 5.31× at 64
//! samples and 5.49× at 128. `docs/experiments/p-40.csv` has been regenerated
//! since and its figures for those two rows have read 4.91× and 5.15× in one run
//! and 5.23× and 5.09× in another — the bench moves by ±0.3 between
//! regenerations, and this harness's own five runs per rung span 0.43. So the
//! honest reading is a band of roughly **5-5.8×** at 64 and 128 samples, which is
//! what all three sources say and is two and a half times the registered 2× bar.
//! Quoting one decimal from one run as *the* stage ratio, in either direction, is
//! what this file's ranges exist to avoid.
//!
//! **The 33³ rung reads low at 4.23-4.33×, and that is a real effect rather than
//! noise:** 32 cells per axis is half a word, so every word operation is paid for
//! by 32 cells instead of 64 while the bitmap build still pays its one comparison
//! per sample. The mechanism's overhead is per *word*, and a half-empty word is
//! its worst case. It is on the ladder for exactly that reason — the fraction
//! argument is strongest at the top and the cost argument weakest at the bottom,
//! and a demo that only showed the rung where both flatter the mechanism would be
//! an advertisement.
//!
//! # The crossover M-337 left unmeasured is reachable from the digit keys
//!
//! M-337's *"would be shown wrong by"* names one fixture and says it is not
//! measured: *"a field whose active fraction is high enough that the bitmap
//! build's own `n³` comparisons stop paying for themselves … `gyroid` at low
//! resolution is the fixture that would probe it."* Press `2` or `3` at the
//! bottom of the ladder and that is what is on screen:
//!
//! One run per cell, `median of five` within it, so read the ratios to about
//! ±0.2 — the spread measured on `sphere` above:
//!
//! | field | 33³ | 64³ | 128³ |
//! |---|---|---|---|
//! | `sphere` | 3.54% · 4.30× | 1.89% · 5.49× | 0.93% · 5.64× |
//! | `gyroid` | **15.99%** · 4.00× | 8.40% · 5.40× | 4.15% · 5.50× |
//! | `noise_cavity` | **18.85%** · 3.73× | 11.03% · 5.59× | 5.71% · 5.61× |
//! | `fbm_terrain` | 5.98% · 4.16× | 3.27% · 5.69× | 1.65% · 5.66× |
//! | `box_exact` | 4.13% · 4.26× | 2.46% · 5.34× | 1.20% · 5.37× |
//!
//! **The crossover is not reached.** At 18.85% active on the smallest grid on the
//! ladder — twenty times `sphere`'s fraction at 128³, and the worst case for the
//! word loop as well — the packed test is still **3.7× faster**. The fraction
//! costs the mechanism about half a turn of ratio and nothing near its margin,
//! which says the win is not really about skipping the gather at all: it is that
//! one comparison per *sample* replaces eight loads per *cell*, and there are
//! eight times more cell-corners than samples however much surface there is.
//! Whether a fraction exists that closes it is still open; 19% is not it.
//!
//! # The two predicates are timed against each other, here, live
//!
//! [`scalar_active_cells`] is the eight-corner gather `DualMesher::place_vertices`
//! used to run on every cell: load eight samples, count how many are inside,
//! keep the cell when the count is neither 0 nor 8. [`packed_active_cells`] is
//! what the crate does now — [`build_bits`] packs `value < 0` one bit per sample
//! along x, and [`active_word`] fuses the four rows bounding a cell row.
//!
//! Both are given the same buffer, in the same **padded row layout** the crate
//! uses: `samples | 1`, which is M-287's cache-set fix. That is live on this
//! ladder rather than incidental — 64 and 128 samples both pay the pad float per
//! row, and 128 is the exact size (`size[0] = size[1] = 128`, a 512-byte row and
//! a 64 KiB plane) where M-287 measured Surface Nets at **3.37× cost** without
//! it. The packed arm rebuilds the bitmap inside every timed rep, because the
//! `n³` comparisons that build it are part of what it costs.
//!
//! # Both arms produce the same ordered list, and that is checked rather than argued
//!
//! `active &= active - 1` clears the lowest set bit, so the set-bit walk visits
//! a row in ascending `x` — the same order the scalar loop did. That is what
//! keeps vertex creation order and therefore **every index in the mesh**
//! unchanged, which is the whole of M-337's clause three. So the two arms here
//! emit cell indices into two `Vec<u32>`s and the vectors are compared element
//! for element, not by length; the HUD line says which. A count-versus-count
//! check would pass on a permutation, and a permutation is exactly the failure
//! that would move every triangle.
//!
//! # Inside is `value < 0`, and that is deliberately not the IEEE sign bit
//!
//! `-0.0` has its sign bit set and `-0.0 < 0.0` is **false**, so `is_sign_negative`
//! is a *different* predicate, one instruction cheaper, and wrong. M-337 names
//! `box_exact` as where it would bite, because that field is exactly zero across
//! its whole boundary — so `box_exact` is on the digit keys and the HUD counts,
//! per field per rung: how many samples are exactly zero, how many of those are
//! `-0.0`, and how many cells a sign-bit variant of the same packed walk
//! disagrees about.
//!
//! **Two things fall out of measuring it that reading the finding would not tell
//! you, and the second one narrows M-337's own justification.**
//!
//! | field | 33 samples | 64 samples | 128 samples |
//! |---|---:|---:|---:|
//! | `box_exact` | **1,538** exact zeros | **0** | **0** |
//! | `sphere` | 6 | 0 | 0 |
//! | `noise_cavity` | 27 | 0 | 0 |
//! | `gyroid` | 1 | 0 | 0 |
//! | `fbm_terrain` | 0 | 0 | 0 |
//!
//! **First: not one exact zero anywhere is `-0.0`, so the sign-bit variant
//! disagrees on 0 cells on every field at every rung.** That is not a reprieve.
//! `box_sample` is `length(max(q, 0)) + min(max q, 0)`; `|p − c| − h` is `+0.0`
//! when the two are equal, `max` of that against zero is `+0.0`, `sqrt(+0.0)` is
//! `+0.0`, and `+0.0 + +0.0` is `+0.0` — the closed form has no path to a
//! negative zero. The moment a field negates a zero, though — `max(a, −b)`, which
//! is every CSG difference in the crate — `-0.0` is back. So the shortcut is
//! unsafe *and* indistinguishable from safe on the fixtures at hand, which is the
//! worst combination a shortcut can have and exactly why it is named in the
//! pre-registration rather than left to be found.
//!
//! **Second: `box_exact` only has exact zeros when cells per axis is a multiple
//! of four, which none of M-337's own grids are.** The box is `[-1, 1]³` in a
//! `[-2, 2]³` domain, so a sample lands on its face only if `1` is a whole
//! multiple of `4 / cells`. At 32 cells it is, 1,538 times. At 63 and 127 it is
//! not, and the count is **zero** — and 63, 127 and 255 are precisely the grids
//! M-337 timed. So the fixture M-337 cites as the reachable case for a signed
//! zero produces no zero at all on the grids M-337 measured. The reasoning is
//! still right and the decision is still right; the *evidence offered for it* only
//! exists at 33 samples, which is why this example has a 33-sample rung and prints
//! the count instead of asserting the claim.
//!
//! # Surface Nets, and `f64`
//!
//! Surface Nets because that is the extractor M-337's clause two held for:
//! 1.25-1.38× end to end on every grid measured, against Dual Contouring's
//! 1.12-1.18×. The mechanism removes a per-**inactive**-cell cost, so its leverage
//! is `inactive work / total work`, and Dual Contouring's QEF solve per *active*
//! cell makes its denominator bigger. `f64` because M-337's numbers are `f64`
//! numbers; the surface is cast to `f32` on its way into the [`Mesh`] asset and
//! nothing but the picture depends on that.
//!
//! # The cloud is budgeted, and says so
//!
//! One merged cube mesh, not one entity or one gizmo per cell. At 128 samples
//! `noise_cavity` has **116,957** active cells and `gyroid` **85,077**; a gizmo
//! cage each would be a million lines a frame. Past [`MAX_CLOUD_CELLS`] the mesh
//! takes every `k`th cell — `k` is 3 for both of those — and the HUD reports `k`,
//! because a silently decimated picture of a *counting* argument would be the one
//! dishonest thing here. Nothing counted is decimated; only what is drawn.
//!
//! # A resolution change is a hitch, and the noisy fields are the expensive ones
//!
//! A rebuild samples the field once for itself and again for the extractor. On
//! `sphere` that is 6 ms and 26 ms at 128 samples; on `fbm_terrain`, four octaves
//! of value noise per sample, it is **372 ms and 437 ms** — most of a second, once
//! per rung. Under capture that is invisible, because a captured frame is a
//! screenshot and not a deadline. Interactively it is a visible pause on `]`, and
//! it is the field's evaluation cost rather than anything this example does: M-296
//! is the finding that field evaluation dominates, and this is what it looks like.

mod common;

use std::time::Instant;

use bevy::prelude::*;
use bevy_isomesh::MeshBuilder;
use common::{Capture, CommonPlugin, DemoDomain, DemoMesh, DemoStats, OrbitCamera, ViewFlags};
use isomesh::fields::{BoxExact, FbmTerrain, ReferenceField, Sphere, capped_gyroid, noise_cavity};
use isomesh::surface_nets::SurfaceNets;
use isomesh::{MeshBuffer, MeshSink, RuntimeShape3, Sdf};

/// Bits in the word the whole mechanism is built on.
const WORD_BITS: usize = 64;

/// The resolution ladder, in samples per axis.
///
/// Cells per axis are one less: 32, 63, 127. The top two are M-337's own grids,
/// so their active counts are comparable to its table digit for digit. See the
/// module docs for why `64k + 1` is not on here and why 33 is the rung the ratio
/// reads worst at.
const LADDER: [u32; 3] = [33, 64, 128];

/// Where an interactive run starts: 64 samples, one word per row of cells and
/// the first row of M-337's table.
const DEFAULT_RUNG: usize = 1;

/// The fields on the digit keys, in order.
///
/// `sphere` leads because it is the field M-337's table is about. `box_exact`
/// is last and is not filler — it is the field the sign-bit trap would bite on.
const FIELD_COUNT: usize = 5;

/// Timed repetitions per arm, per rebuild. M-337 used five and took the median.
const REPS: usize = 5;

/// Active cells the merged cube cloud will draw before it starts striding.
///
/// 40,000 cubes is 960,000 vertices and 480,000 triangles, which is a rounding
/// error for the renderer and about 40 ms to build on the CPU. `noise_cavity` at
/// 128 samples goes well past it at 116,957 active cells, so the stride exists
/// and is reported rather than being a silent cap.
const MAX_CLOUD_CELLS: usize = 40_000;

/// Cube edge as a fraction of the cell, for the cloud.
///
/// Well under 1 so the lattice has gaps to see the grey surface through. At 1.0
/// the cloud is a solid skin and the picture stops being about cells.
const CLOUD_FILL: f32 = 0.45;

/// Seconds between cursor steps in an interactive run.
const CURSOR_INTERVAL: f32 = 0.09;

// ─── the strip ──────────────────────────────────────────────────────────────

/// Rows in the bit strip: four `inside` words, then `any`, `all`, `active`.
const STRIP_ROWS: usize = 7;

/// One bit's box, and the column stride.
const BIT_W: f32 = 5.0;
const BIT_STRIDE_X: f32 = 6.0;

/// One row's box, and the row stride.
const BIT_H: f32 = 9.0;
const BIT_STRIDE_Y: f32 = 11.0;

/// Where bit 0 of the bottom row sits, in logical pixels from the lower left.
const STRIP_LEFT: f32 = 104.0;
const STRIP_BOTTOM: f32 = 34.0;

/// The panel behind it.
const PANEL_LEFT: f32 = 10.0;
const PANEL_BOTTOM: f32 = 10.0;

/// Unlit, and not-in-the-grid-at-all: a bit past the last cell in the fused rows,
/// or past the last sample in the `inside` rows.
const BIT_OFF: Color = Color::srgb(0.13, 0.15, 0.21);
const BIT_VOID: Color = Color::srgb(0.06, 0.07, 0.10);
/// A set `inside` bit.
const BIT_INSIDE: Color = Color::srgb(0.55, 0.68, 0.95);
/// A set bit of `any`, of `all`, and of `active`.
const BIT_ANY: Color = Color::srgb(0.28, 0.80, 0.45);
const BIT_ALL: Color = Color::srgb(0.96, 0.55, 0.20);
const BIT_ACTIVE: Color = Color::srgb(1.0, 0.92, 0.25);

/// The run of cells the highlighted word decides, up to sixty-four of them.
const WORD_COLOUR: Color = Color::srgb(0.95, 0.96, 1.0);

/// Where the subject sits in frame, as a fraction of the orbit radius.
///
/// Right of centre and a little above it, because the HUD owns the upper left
/// and the bit strip owns the lower left. Applied in the camera's own basis so
/// it survives `ISOMESH_SPIN`, which is the correction E-304 made after E-109
/// photographed its argument with the evidence on top of it.
const SUBJECT_OFFSET: Vec2 = Vec2::new(0.17, -0.07);

/// Orbit radius as a multiple of the domain extent.
///
/// From the field's own `domain()` rather than a constant: the compact four are
/// half-extent 2 and `fbm_terrain` is 8, so a fixed radius would put the camera
/// inside the terrain.
///
/// **1.85 rather than something tighter, because the wire box is half the
/// picture.** The claim is a ratio of the shell to the whole grid, and at 1.35
/// the box ran off every edge of a 1280x720 frame and read as three stray grey
/// lines — the denominator of the argument, invisible. At 1.85 the whole domain
/// is in shot at the cost of the shell being a third of the frame instead of
/// two-thirds, which is the right way round for this demo.
const VIEW_RADIUS_EXTENTS: f32 = 1.85;

// ─── the predicates ─────────────────────────────────────────────────────────

/// Inside, the way every extractor in the crate decides it.
///
/// `cube::is_inside` is `pub(crate)`, so this is the same one line rather than a
/// reference to it. Strictly negative — see the module docs for why this is not
/// `is_sign_negative`.
#[inline]
fn is_inside(value: f64) -> bool {
    value < 0.0
}

/// The grid the field is sampled on, laid out the way `DualMesher` lays it out.
///
/// One definition, so the two arms, the extraction and the picture cannot
/// disagree about which lattice they looked at.
struct Grid {
    /// World position of sample `[0, 0, 0]`.
    origin: [f64; 3],
    cell_size: f64,
    samples: u32,
    /// Cells per axis: `samples - 1`.
    cells: u32,
    /// Samples per row in the buffer, `samples | 1`.
    ///
    /// The crate's own padding (M-287): a row stride that is a multiple of 512
    /// bytes aliases cache sets and cost Surface Nets 3.37× at 128³ on a field
    /// with no surface at all. Forcing it odd removes that, and both arms here
    /// are measured on the layout the crate actually uses — 64 and 128 samples
    /// both take the pad, and 128 is the exact size M-287 measured, so this is
    /// load-bearing on this ladder rather than transcribed for completeness.
    row: usize,
    /// `u64` words per sample row.
    bit_row: usize,
}

impl Grid {
    fn new(origin: [f64; 3], cell_size: f64, samples: u32) -> Self {
        Self {
            origin,
            cell_size,
            samples,
            cells: samples.saturating_sub(1),
            row: samples as usize | 1,
            bit_row: (samples as usize).div_ceil(WORD_BITS),
        }
    }

    /// Where sample `p` lives in the value buffer. `DualMesher::index`.
    #[inline]
    fn index(&self, p: [u32; 3]) -> usize {
        p[0] as usize + self.row * (p[1] as usize + self.samples as usize * p[2] as usize)
    }

    fn point(&self, x: usize, y: usize, z: usize) -> [f64; 3] {
        [
            self.origin[0] + self.cell_size * x as f64,
            self.origin[1] + self.cell_size * y as f64,
            self.origin[2] + self.cell_size * z as f64,
        ]
    }

    /// Cells per axis as a `usize`, which is what every loop below wants.
    #[inline]
    fn cells_usize(&self) -> usize {
        self.cells as usize
    }
}

/// Corner `c` of a cell, as an offset. `cube::corner_offset`, which is
/// `pub(crate)`: bit `i` of the corner index is axis `i`.
#[inline]
const fn corner_offset(corner: usize) -> [u32; 3] {
    [
        (corner & 1) as u32,
        ((corner >> 1) & 1) as u32,
        ((corner >> 2) & 1) as u32,
    ]
}

/// How many samples of a field are exactly zero, and how many of those are
/// `-0.0`.
#[derive(Default, Clone, Copy)]
struct Zeros {
    exact: usize,
    negative: usize,
}

/// Sample the field over the grid, counting exact zeros on the way past.
///
/// Counted here rather than in a second pass over the buffer because the buffer
/// carries padding slots that are written zero and are never read — a filter
/// over it would count them.
fn sample_grid<F: Sdf<Scalar = f64>>(field: &F, grid: &Grid, values: &mut Vec<f64>) -> Zeros {
    let n = grid.samples as usize;
    let pad = grid.row - n;
    values.clear();
    values.reserve(grid.row * n * n);
    let mut zeros = Zeros::default();
    for z in 0..n {
        for y in 0..n {
            for x in 0..n {
                let v = field.sample(grid.point(x, y, z));
                if v == 0.0 {
                    zeros.exact += 1;
                    if v.is_sign_negative() {
                        zeros.negative += 1;
                    }
                }
                values.push(v);
            }
            for _ in 0..pad {
                values.push(0.0);
            }
        }
    }
    zeros
}

/// The predicate the crate used to run: eight loads and eight compares, on every
/// cell, including the ~97% that produce nothing.
///
/// Emits the cell index of every active cell, in the lexicographic order the
/// loop visits them. `DualMesher` did the vertex placement inline where this
/// pushes; both arms push, so the push is not the difference between them.
fn scalar_active_cells(values: &[f64], grid: &Grid, out: &mut Vec<u32>) {
    let c = grid.cells_usize();
    out.clear();
    for z in 0..c {
        for y in 0..c {
            for x in 0..c {
                let base = [x as u32, y as u32, z as u32];
                let mut inside = 0usize;
                for corner in 0..8usize {
                    let o = corner_offset(corner);
                    let s = grid.index([base[0] + o[0], base[1] + o[1], base[2] + o[2]]);
                    if is_inside(values[s]) {
                        inside += 1;
                    }
                }
                if inside != 0 && inside != 8 {
                    out.push((x + y * c + z * c * c) as u32);
                }
            }
        }
    }
}

/// Pack one bit per sample along x, under whichever sign predicate is passed.
///
/// Generic rather than a function pointer so the predicate inlines: the point of
/// the bitmap is that the sign test stops being a branch per corner, and an
/// indirect call in the inner loop would measure the opposite of that. The
/// `is_sign_negative` instantiation is what makes the trap a number on the HUD
/// instead of a paragraph.
fn build_bits<P: Fn(f64) -> bool>(values: &[f64], grid: &Grid, inside: &mut Vec<u64>, bit: P) {
    let sx = grid.samples as usize;
    let rows = sx * sx;
    inside.clear();
    inside.resize(grid.bit_row * rows, 0);

    for row in 0..rows {
        let src = grid.row * row;
        let dst = grid.bit_row * row;
        for w in 0..grid.bit_row {
            let base = w * WORD_BITS;
            let n = (sx - base).min(WORD_BITS);
            let mut word = 0u64;
            for k in 0..n {
                word |= u64::from(bit(values[src + base + k])) << k;
            }
            inside[dst + w] = word;
        }
    }
}

/// Bit `k` is the sign of sample `[64w + k, y, z]`.
#[inline]
fn inside_word(inside: &[u64], grid: &Grid, w: usize, y: usize, z: usize) -> u64 {
    inside[grid.bit_row * (y + grid.samples as usize * z) + w]
}

/// Bit `k` is the sign of sample `[64w + k + 1, y, z]`.
///
/// The high bit comes from the next word, or the cell straddling a word boundary
/// reads its `+x` corner as outside — a hole every 64 cells, invisible at 33 and
/// 64 samples where there is one word per row and fatal at 128 where there are
/// two.
#[inline]
fn inside_word_shifted(inside: &[u64], grid: &Grid, w: usize, y: usize, z: usize) -> u64 {
    let lo = inside_word(inside, grid, w, y, z);
    let hi = if w + 1 < grid.bit_row {
        inside_word(inside, grid, w + 1, y, z)
    } else {
        0
    };
    (lo >> 1) | (hi << 63)
}

/// The three fused words for one cell row, and the four `inside` words they came
/// from.
///
/// This is the mechanism, and it is returned whole rather than reduced to
/// `active` because the strip draws all seven: a viewer who cannot see `any` and
/// `all` has no reason to believe `active`.
#[derive(Clone, Copy, Default)]
struct Fused {
    /// `inside` at `(y,z)`, `(y+1,z)`, `(y,z+1)`, `(y+1,z+1)`.
    rows: [u64; 4],
    any: u64,
    all: u64,
    active: u64,
}

/// Sixty-four active-cell answers in about twenty word operations.
///
/// ```text
/// any    = OR  over the four rows of (a | b)      some corner inside
/// all    = AND over the four rows of (a & b)      every corner inside
/// active = any & !all                             mixed, i.e. crossed
/// ```
#[inline]
fn active_word(inside: &[u64], grid: &Grid, w: usize, y: usize, z: usize) -> Fused {
    let mut out = Fused {
        any: 0,
        all: !0,
        ..Fused::default()
    };
    for dz in 0..2usize {
        for dy in 0..2usize {
            let a = inside_word(inside, grid, w, y + dy, z + dz);
            let b = inside_word_shifted(inside, grid, w, y + dy, z + dz);
            out.rows[dz * 2 + dy] = a;
            out.any |= a | b;
            out.all &= a & b;
        }
    }
    out.active = out.any & !out.all;
    out
}

/// The low `count` bits, with the full-word case named.
///
/// `1u64 << 64` is undefined behaviour's cousin — in Rust it panics in debug and
/// is a masked shift in release, which is worse — so 64 and over is a separate
/// arm rather than an arithmetic accident.
#[inline]
fn low_mask(count: usize) -> u64 {
    if count >= WORD_BITS {
        !0
    } else {
        (1u64 << count) - 1
    }
}

/// Which bits of word `w` are cells that exist.
///
/// The last word of a row runs past the final cell whenever cells per axis is
/// not a multiple of 64 — at 33 samples that is half the word, and at 64 and 128
/// it is exactly the top bit.
#[inline]
fn cell_mask(w: usize, cells_x: usize) -> u64 {
    low_mask(cells_x.saturating_sub(w * WORD_BITS))
}

/// Which bits of word `w` are samples that exist.
///
/// One more than [`cell_mask`]'s, always, and the strip draws the difference:
/// `n` cells need `n + 1` samples, which is the whole reason
/// [`inside_word_shifted`] has to reach into the next word.
#[inline]
fn sample_mask(w: usize, samples: usize) -> u64 {
    low_mask(samples.saturating_sub(w * WORD_BITS))
}

/// The predicate the crate runs now.
///
/// Emits the same cell indices as [`scalar_active_cells`], in the same order:
/// `active &= active - 1` clears the lowest set bit, so a row is walked in
/// ascending `x`.
fn packed_active_cells(inside: &[u64], grid: &Grid, out: &mut Vec<u32>) {
    let c = grid.cells_usize();
    out.clear();
    for z in 0..c {
        for y in 0..c {
            for w in 0..grid.bit_row {
                let mut active = active_word(inside, grid, w, y, z).active & cell_mask(w, c);
                while active != 0 {
                    let x = w * WORD_BITS + active.trailing_zeros() as usize;
                    active &= active - 1;
                    out.push((x + y * c + z * c * c) as u32);
                }
            }
        }
    }
}

/// How many cells two ascending lists disagree about.
///
/// A linear merge, because both walks emit in lexicographic order — which is the
/// property [`packed_active_cells`] exists to preserve, so using it here is not
/// an assumption.
fn disagreements(a: &[u32], b: &[u32]) -> usize {
    let (mut i, mut j, mut count) = (0usize, 0usize, 0usize);
    while i < a.len() && j < b.len() {
        match a[i].cmp(&b[j]) {
            std::cmp::Ordering::Equal => {
                i += 1;
                j += 1;
            }
            std::cmp::Ordering::Less => {
                i += 1;
                count += 1;
            }
            std::cmp::Ordering::Greater => {
                j += 1;
                count += 1;
            }
        }
    }
    count + (a.len() - i) + (b.len() - j)
}

/// Median of a list of timings, which is what M-337 reports.
///
/// A mean over five reps lets one scheduler tick set the answer, and this runs
/// on a machine that is also rendering.
fn median(mut values: Vec<f64>) -> f64 {
    if values.is_empty() {
        return 0.0;
    }
    values.sort_by(f64::total_cmp);
    values[values.len() / 2]
}

// ─── what one rebuild produced ──────────────────────────────────────────────

/// One row's word, kept because the strip draws it.
#[derive(Clone, Copy, Default)]
struct Word {
    w: u32,
    y: u32,
    z: u32,
    fused: Fused,
    /// Which of its bits are cells.
    cells: u64,
    /// Which of its bits are samples. Always one more than [`Self::cells`].
    samples: u64,
}

/// Everything on screen that is not geometry.
///
/// Held as a resource rather than recomputed per frame: a rebuild is tens to
/// hundreds of milliseconds and nothing in it changes until the field or the
/// resolution does.
#[derive(Resource, Default)]
struct Report {
    field_name: &'static str,
    samples: u32,
    cells: u32,
    total_cells: usize,
    bit_row: usize,
    active: usize,
    fraction: f64,
    scalar_ns: f64,
    packed_ns: f64,
    ratio: f64,
    lists_match: bool,
    zeros: Zeros,
    sign_bit_diff: usize,
    sample_ms: f64,
    /// Every word that contains at least one active cell, in `(z, y, w)` order.
    /// The cursor walks this, so the strip is never a blank row.
    words: Vec<Word>,
    origin: Vec3,
    cell_size: f32,
    cloud_stride: usize,
    cloud_drawn: usize,
}

/// Which word the strip and the highlight are showing.
#[derive(Resource, Default)]
struct Cursor {
    index: usize,
    timer: f32,
}

/// Which field is showing.
#[derive(Resource)]
struct Field(usize);

/// Which rung of [`LADDER`] is showing.
#[derive(Resource)]
struct Rung(usize);

/// A resolution pinned by `ISOMESH_SAMPLES`, which takes the ladder out of play.
///
/// Clamped, because the value is arbitrary user input and every cost here is
/// cubic: at 1,000 samples the value buffer alone is 8 GB and Surface Nets'
/// per-cell `[u32; 12]` is another 48. [`MAX_SAMPLES`] is 129 rather than
/// something rounder so that the `64k + 1` finding above stays reproducible from
/// the environment — 65 and 129 are the two grids that show it.
#[derive(Resource)]
struct Pinned(Option<u32>);

/// Below this there is no cell to test.
const MIN_SAMPLES: u32 = 2;

/// Above this a rebuild allocates more than the demo is worth.
const MAX_SAMPLES: u32 = 129;

/// Whether the extracted surface is drawn.
#[derive(Resource)]
struct ShowSurface(bool);

impl Default for ShowSurface {
    /// `ISOMESH_SURFACE=off` starts with the surface hidden, so the naked shell
    /// can be captured without a keyboard — the same contract `ISOMESH_VIEW` and
    /// `ISOMESH_FIELD` provide, for the same reason: a view reachable only by
    /// pressing `H` is a view no committed image can be regenerated from.
    fn default() -> Self {
        Self(!matches!(
            std::env::var("ISOMESH_SURFACE")
                .unwrap_or_default()
                .trim()
                .to_ascii_lowercase()
                .as_str(),
            "off" | "0" | "no" | "hide"
        ))
    }
}

#[derive(Resource)]
struct SurfaceMaterial(Handle<StandardMaterial>);

#[derive(Resource)]
struct CloudMaterial(Handle<StandardMaterial>);

/// The merged cube mesh of active cells.
///
/// Deliberately **not** tagged [`DemoMesh`]: the harness's `W` wireframe reads
/// the mesh back and submits three gizmo lines per triangle, which at 40,000
/// cubes would be 1.4 million lines in a frame.
#[derive(Component)]
struct CloudMesh;

/// The highlighted 64-cell run gets its own group so it draws in front of the
/// translucent surface without dragging the shared wireframe forward with it.
#[derive(Default, Reflect, GizmoConfigGroup)]
struct RowGizmos;

/// One box in the bit strip.
#[derive(Component)]
struct BitCell {
    row: usize,
    col: usize,
}

/// Everything the strip is made of, so `nohud` hides all of it at once.
#[derive(Component)]
struct StripPanel;

/// The line under the strip that says which word this is.
#[derive(Component)]
struct StripFooter;

/// The fraction, in the largest type on screen.
#[derive(Component)]
struct BigFraction;

/// The count under it.
#[derive(Component)]
struct BigCount;

// ─── app ────────────────────────────────────────────────────────────────────

fn main() {
    App::new()
        .add_plugins(DefaultPlugins.set(WindowPlugin {
            primary_window: Some(Window {
                title: "isomesh — E-307 active cells".into(),
                ..default()
            }),
            ..default()
        }))
        .add_plugins(CommonPlugin)
        .init_gizmo_group::<RowGizmos>()
        .insert_resource(Pinned(
            common::samples_override().map(|n| n.clamp(MIN_SAMPLES, MAX_SAMPLES)),
        ))
        .insert_resource(Field(0))
        .insert_resource(Rung(DEFAULT_RUNG))
        .init_resource::<ShowSurface>()
        .init_resource::<Report>()
        .init_resource::<Cursor>()
        .add_systems(Startup, setup)
        .add_systems(
            Update,
            (
                controls,
                rebuild,
                advance_cursor,
                frame_camera,
                update_strip,
                update_big_number,
                report,
                draw_row,
            )
                .chain(),
        )
        .run();
}

/// Captured frames spent on each rung.
///
/// Read from the environment rather than from [`Capture`], which keeps its
/// length private, and the alternative is editing the harness. Three rungs at
/// `frames/3` each means the default 60-frame capture is exactly one pass up
/// the ladder and a six-frame smoke capture still visits all three.
fn capture_frames_per_rung() -> u32 {
    let frames: u32 = std::env::var("ISOMESH_CAPTURE_FRAMES")
        .ok()
        .and_then(|v| v.parse().ok())
        .unwrap_or(60);
    (frames / LADDER.len() as u32).max(1)
}

/// Words the cursor advances per step.
///
/// One per step would move the highlight a couple of cells across a 30-frame rung
/// and read as frozen; the whole list in 30 frames would teleport it and read as
/// noise. Capped at eight, which on a 640x360 capture sweeps the bar visibly
/// along the shell without losing the eye.
fn cursor_step(words: usize) -> usize {
    (words / 240).clamp(1, 8)
}

/// Where the cursor starts after a rebuild: the middle of the list.
///
/// **Not zero, and this was measured rather than guessed.** The word list is
/// built in `(z, y, w)` order, so index 0 is the lowest `z` the surface reaches —
/// the row where it is *tangent* to the grid. There `all` is zero, because no cell
/// has all eight corners inside, so the strip shows `active == any` and the `AND`
/// row is an empty line: a picture of the mechanism with one of its three terms
/// doing nothing. The middle of the list is a row through the thick of the
/// surface, where `any`, `all` and `active` are three different words and
/// `any & !all` visibly subtracts something.
fn cursor_start(words: usize) -> usize {
    words / 2
}

fn setup(
    mut commands: Commands,
    mut materials: ResMut<Assets<StandardMaterial>>,
    mut gizmo_config: ResMut<GizmoConfigStore>,
    mut camera: Query<&mut OrbitCamera>,
) {
    // Off both axes. Down an axis the cell lattice collapses into a grid of
    // squares and the shell stops looking like one cell thick, which is the
    // thing the picture is for.
    for mut orbit in &mut camera {
        orbit.yaw = 0.68;
        orbit.pitch = 0.36;
    }

    let (rows, _) = gizmo_config.config_mut::<RowGizmos>();
    rows.line.width = 2.0;
    rows.depth_bias = -0.6;

    // Translucent and double-sided, or it hides its own subject: the cyan shell
    // straddles the surface, so half of every cube is behind it.
    commands.insert_resource(SurfaceMaterial(materials.add(StandardMaterial {
        base_color: Color::srgba(0.72, 0.76, 0.82, 0.22),
        perceptual_roughness: 0.5,
        alpha_mode: AlphaMode::Blend,
        double_sided: true,
        cull_mode: None,
        ..default()
    })));

    // Opaque, so 40,000 cubes do not need sorting to look like a lattice. The
    // emissive term is what keeps a cell on the dark side of the shell from
    // reading as a hole.
    commands.insert_resource(CloudMaterial(materials.add(StandardMaterial {
        base_color: Color::srgb(0.12, 0.66, 0.86),
        emissive: LinearRgba::rgb(0.02, 0.16, 0.24),
        perceptual_roughness: 0.6,
        ..default()
    })));

    // Spawned rather than assumed: `draw_domain` queries for it, and without one
    // the `G` toggle silently does nothing. Filled in by the first rebuild.
    commands.spawn(DemoDomain {
        min: Vec3::splat(-1.0),
        max: Vec3::splat(1.0),
    });

    spawn_strip(&mut commands);
    spawn_big_number(&mut commands);
}

/// The strip: seven rows of sixty-four boxes, plus labels and a footer.
///
/// Root-level absolutely-positioned nodes rather than a flex grid, so a bit's
/// position is arithmetic rather than a layout outcome — E-305's plot learned
/// that, and 448 boxes is not a place to start guessing about flex.
fn spawn_strip(commands: &mut Commands) {
    let strip_w = WORD_BITS as f32 * BIT_STRIDE_X;
    let top_row = STRIP_BOTTOM + (STRIP_ROWS - 1) as f32 * BIT_STRIDE_Y;

    commands.spawn((
        Node {
            position_type: PositionType::Absolute,
            left: Val::Px(PANEL_LEFT),
            bottom: Val::Px(PANEL_BOTTOM),
            width: Val::Px(STRIP_LEFT - PANEL_LEFT + strip_w + 10.0),
            height: Val::Px(top_row + BIT_H + 20.0 - PANEL_BOTTOM),
            ..default()
        },
        BackgroundColor(Color::srgba(0.05, 0.06, 0.09, 0.66)),
        GlobalZIndex(1),
        StripPanel,
    ));

    commands.spawn((
        Text::new("64 cells in one word:  a=inside(x)  b=a>>1  any=OR(a|b)  all=AND(a&b)  active=any&!all"),
        TextFont {
            font_size: FontSize::Px(11.0),
            ..default()
        },
        TextColor(Color::srgb(0.82, 0.86, 0.92)),
        Node {
            position_type: PositionType::Absolute,
            left: Val::Px(PANEL_LEFT + 6.0),
            bottom: Val::Px(top_row + BIT_H + 4.0),
            ..default()
        },
        GlobalZIndex(2),
        StripPanel,
    ));

    const LABELS: [&str; STRIP_ROWS] = [
        "inside  y, z",
        "inside y+1,z",
        "inside  y,z+1",
        "inside y+1,z+1",
        "any = OR",
        "all = AND",
        "ACTIVE",
    ];
    for (row, label) in LABELS.iter().enumerate() {
        let bottom = STRIP_BOTTOM + (STRIP_ROWS - 1 - row) as f32 * BIT_STRIDE_Y;
        commands.spawn((
            Text::new(*label),
            TextFont {
                font_size: FontSize::Px(10.0),
                ..default()
            },
            TextColor(if row < 4 {
                Color::srgb(0.68, 0.74, 0.86)
            } else {
                Color::srgb(0.92, 0.94, 0.98)
            }),
            Node {
                position_type: PositionType::Absolute,
                left: Val::Px(PANEL_LEFT + 6.0),
                bottom: Val::Px(bottom - 1.0),
                ..default()
            },
            GlobalZIndex(2),
            StripPanel,
        ));
        for col in 0..WORD_BITS {
            commands.spawn((
                Node {
                    position_type: PositionType::Absolute,
                    left: Val::Px(STRIP_LEFT + col as f32 * BIT_STRIDE_X),
                    bottom: Val::Px(bottom),
                    width: Val::Px(BIT_W),
                    height: Val::Px(BIT_H),
                    ..default()
                },
                BackgroundColor(BIT_OFF),
                GlobalZIndex(2),
                BitCell { row, col },
                StripPanel,
            ));
        }
    }

    commands.spawn((
        Text::new(""),
        TextFont {
            font_size: FontSize::Px(10.0),
            ..default()
        },
        TextColor(Color::srgb(0.78, 0.82, 0.90)),
        Node {
            position_type: PositionType::Absolute,
            left: Val::Px(PANEL_LEFT + 6.0),
            bottom: Val::Px(PANEL_BOTTOM + 4.0),
            ..default()
        },
        GlobalZIndex(2),
        StripFooter,
        StripPanel,
    ));
}

/// The fraction, upper right, in 30px type.
///
/// Its own node rather than a HUD line because the shared HUD is one font size
/// and this number is the argument. Right-anchored and absolutely positioned, so
/// it hugs its own content and cannot walk into the HUD as the digits change.
fn spawn_big_number(commands: &mut Commands) {
    commands.spawn((
        Text::new(""),
        TextFont {
            font_size: FontSize::Px(30.0),
            ..default()
        },
        TextColor(Color::srgb(1.0, 0.92, 0.25)),
        Node {
            position_type: PositionType::Absolute,
            top: Val::Px(10.0),
            right: Val::Px(16.0),
            ..default()
        },
        GlobalZIndex(2),
        BigFraction,
        StripPanel,
    ));
    commands.spawn((
        Text::new(""),
        TextFont {
            font_size: FontSize::Px(13.0),
            ..default()
        },
        TextColor(Color::srgb(0.86, 0.90, 0.96)),
        Node {
            position_type: PositionType::Absolute,
            top: Val::Px(46.0),
            right: Val::Px(16.0),
            ..default()
        },
        GlobalZIndex(2),
        BigCount,
        StripPanel,
    ));
}

/// Field, resolution and the surface toggle.
///
/// Under capture the rung follows the frame counter, because an example whose
/// subject only changes on a keypress captures as a still frame. A pinned
/// `ISOMESH_SAMPLES` wins over both, which is the harness's contract.
fn controls(
    keys: Res<ButtonInput<KeyCode>>,
    capture: Res<Capture>,
    flags: Res<ViewFlags>,
    pinned: Res<Pinned>,
    mut field: ResMut<Field>,
    mut rung: ResMut<Rung>,
    mut surface: ResMut<ShowSurface>,
) {
    field.0 = flags.field.min(FIELD_COUNT - 1);

    if pinned.0.is_none() {
        if capture.is_active() {
            rung.0 = (capture.taken / capture_frames_per_rung()) as usize % LADDER.len();
        } else {
            if keys.just_pressed(KeyCode::BracketRight) && rung.0 + 1 < LADDER.len() {
                rung.0 += 1;
            }
            if keys.just_pressed(KeyCode::BracketLeft) && rung.0 > 0 {
                rung.0 -= 1;
            }
        }
    }

    if keys.just_pressed(KeyCode::KeyH) {
        surface.0 = !surface.0;
    }
}

/// Sample the field, time both predicates, extract the surface, build the cloud.
#[allow(clippy::too_many_arguments)]
fn rebuild(
    field: Res<Field>,
    rung: Res<Rung>,
    pinned: Res<Pinned>,
    flags: Res<ViewFlags>,
    mut report: ResMut<Report>,
    mut cursor: ResMut<Cursor>,
    mut stats: ResMut<DemoStats>,
    mut meshes: ResMut<Assets<Mesh>>,
    mut surface_query: Query<&mut Mesh3d, (With<DemoMesh>, Without<CloudMesh>)>,
    mut cloud_query: Query<&mut Mesh3d, (With<CloudMesh>, Without<DemoMesh>)>,
    mut domain: Query<&mut DemoDomain>,
    mut commands: Commands,
    surface_material: Res<SurfaceMaterial>,
    cloud_material: Res<CloudMaterial>,
    mut last: Local<Option<(usize, u32)>>,
) {
    let samples = pinned
        .0
        .unwrap_or_else(|| LADDER[rung.0.min(LADDER.len() - 1)]);
    let key = (field.0, samples);
    if *last == Some(key) && !flags.remesh_requested {
        return;
    }
    *last = Some(key);

    let Some(built) = build(field.0, samples) else {
        return;
    };

    for mut d in &mut domain {
        d.min = built.domain_min;
        d.max = built.domain_max;
    }

    stats.vertices = built.vertices;
    stats.triangles = built.triangles;
    stats.extract_ms = built.extract_ms;
    *report = built.report;
    cursor.index = cursor_start(report.words.len());
    cursor.timer = 0.0;

    let surface = mesh3d(&mut meshes, built.surface);
    if surface_query.is_empty() {
        commands.spawn((
            surface,
            MeshMaterial3d(surface_material.0.clone()),
            DemoMesh,
        ));
    } else {
        for mut mesh in &mut surface_query {
            *mesh = surface.clone();
        }
    }

    let cloud = mesh3d(&mut meshes, built.cloud);
    if cloud_query.is_empty() {
        commands.spawn((cloud, MeshMaterial3d(cloud_material.0.clone()), CloudMesh));
    } else {
        for mut mesh in &mut cloud_query {
            *mesh = cloud.clone();
        }
    }
}

/// A [`Mesh3d`] for a result that may be empty.
///
/// **An empty [`Mesh`] must not become an asset.** `bevy_render`'s
/// `MeshAllocator::allocate_meshes` skips any mesh whose vertex buffer is zero
/// bytes and then copies into it unconditionally, so one empty mesh logs
/// `Use-after-free: attempted to copy element data for an unallocated key` twice
/// a frame, forever. E-305 found that at chunk scale; it is reachable here from
/// `ISOMESH_SAMPLES=2`, where the single cell spans the whole domain and the
/// sphere is nowhere near its corners — measured, four errors before the first
/// screenshot. `Mesh3d::default()` names no asset and draws nothing, which is
/// what an empty result actually wants.
fn mesh3d(meshes: &mut Assets<Mesh>, mesh: Option<Mesh>) -> Mesh3d {
    mesh.map_or_else(Mesh3d::default, |m| Mesh3d(meshes.add(m)))
}

/// Everything one rebuild produced.
///
/// Both meshes are optional for the reason [`mesh3d`] explains.
struct Built {
    surface: Option<Mesh>,
    cloud: Option<Mesh>,
    report: Report,
    domain_min: Vec3,
    domain_max: Vec3,
    vertices: usize,
    triangles: usize,
    extract_ms: f64,
}

/// Dispatch on the field index, then do the work once in [`census`].
///
/// The reference fields are separate types, so a runtime choice has to be a
/// match rather than a loop over a list — the same shape `critical_cells` and
/// `manifold_check` use.
fn build(field: usize, samples: u32) -> Option<Built> {
    match field {
        0 => census(&Sphere::<f64>::canonical(), samples),
        1 => census(&capped_gyroid::<f64>(), samples),
        2 => census(&noise_cavity::<f64>(), samples),
        3 => census(&FbmTerrain::<f64>::canonical(), samples),
        _ => census(&BoxExact::<f64>::canonical(), samples),
    }
}

fn census<F>(field: &F, samples: u32) -> Option<Built>
where
    F: Sdf<Scalar = f64> + ReferenceField,
{
    let (min, max) = field.domain();
    let cells = samples.saturating_sub(1);
    if cells == 0 {
        error!("{samples} samples per axis leaves no cells to test");
        return None;
    }
    let grid = Grid::new(min, (max[0] - min[0]) / f64::from(cells), samples);
    let total_cells = grid.cells_usize().pow(3);

    // ── the grid, sampled once, shared by both arms ──────────────────────────
    let mut values = Vec::new();
    let started = Instant::now();
    let zeros = sample_grid(field, &grid, &mut values);
    let sample_ms = started.elapsed().as_secs_f64() * 1000.0;

    // ── the two predicates, timed against each other ────────────────────────
    //
    // A warm-up outside the loop, so the first rep is not paying for the two
    // `Vec`s' capacity. After it both arms allocate nothing.
    let mut inside = Vec::new();
    let mut scalar_cells = Vec::new();
    let mut packed_cells = Vec::new();
    scalar_active_cells(&values, &grid, &mut scalar_cells);
    build_bits(&values, &grid, &mut inside, is_inside);
    packed_active_cells(&inside, &grid, &mut packed_cells);

    let per_cell =
        |elapsed: std::time::Duration| elapsed.as_secs_f64() * 1.0e9 / total_cells as f64;
    let mut scalar_ns = Vec::with_capacity(REPS);
    let mut packed_ns = Vec::with_capacity(REPS);
    for _ in 0..REPS {
        let started = Instant::now();
        scalar_active_cells(&values, &grid, &mut scalar_cells);
        scalar_ns.push(per_cell(started.elapsed()));

        // The bitmap is rebuilt inside the timed region on purpose: its `n³`
        // comparisons are part of what the packed arm costs, and hoisting them
        // out would be measuring only the half of the mechanism that wins.
        let started = Instant::now();
        build_bits(&values, &grid, &mut inside, is_inside);
        packed_active_cells(&inside, &grid, &mut packed_cells);
        packed_ns.push(per_cell(started.elapsed()));
    }
    let scalar_ns = median(scalar_ns);
    let packed_ns = median(packed_ns);

    // Element for element, not by count: a permutation would pass a count check
    // and would change every index in the mesh.
    let lists_match = scalar_cells == packed_cells;
    if !lists_match {
        let first = scalar_cells
            .iter()
            .zip(packed_cells.iter())
            .position(|(a, b)| a != b);
        error!(
            "the scalar and packed active-cell walks disagree on {} at {samples}^3: \
             {} against {} cells, first difference at position {first:?}. M-337's clause \
             three rests on these being the same ordered list, so every vertex index \
             the extractor would produce is now in question.",
            F::NAME,
            scalar_cells.len(),
            packed_cells.len(),
        );
    }

    // ── the trap, as a number ───────────────────────────────────────────────
    let mut sign_bits = Vec::new();
    let mut sign_cells = Vec::new();
    build_bits(&values, &grid, &mut sign_bits, f64::is_sign_negative);
    packed_active_cells(&sign_bits, &grid, &mut sign_cells);
    let sign_bit_diff = disagreements(&packed_cells, &sign_cells);

    // ── the words worth showing ─────────────────────────────────────────────
    let c = grid.cells_usize();
    let mut words = Vec::new();
    for z in 0..c {
        for y in 0..c {
            for w in 0..grid.bit_row {
                let cells_here = cell_mask(w, c);
                let fused = active_word(&inside, &grid, w, y, z);
                if fused.active & cells_here != 0 {
                    words.push(Word {
                        w: w as u32,
                        y: y as u32,
                        z: z as u32,
                        fused,
                        cells: cells_here,
                        samples: sample_mask(w, samples as usize),
                    });
                }
            }
        }
    }

    // ── the surface ─────────────────────────────────────────────────────────
    let shape = match RuntimeShape3::new([samples; 3]) {
        Ok(shape) => shape,
        Err(error) => {
            error!("grid {samples}^3 rejected: {error}");
            return None;
        }
    };
    let mut buffer = MeshBuffer::<f64>::new();
    let started = Instant::now();
    if let Err(error) =
        SurfaceNets::<f64>::new().extract(field, &shape, grid.origin, grid.cell_size, &mut buffer)
    {
        error!("surface nets failed at {samples}^3 on {}: {error}", F::NAME);
        return None;
    }
    let extract_ms = started.elapsed().as_secs_f64() * 1000.0;

    // ── the cloud ───────────────────────────────────────────────────────────
    let stride = packed_cells.len().div_ceil(MAX_CLOUD_CELLS).max(1);
    let fill = grid.cell_size as f32 * CLOUD_FILL;
    let inset = (grid.cell_size as f32 - fill) * 0.5;
    let origin = Vec3::new(
        grid.origin[0] as f32,
        grid.origin[1] as f32,
        grid.origin[2] as f32,
    );
    let h = grid.cell_size as f32;
    let mut cloud = MeshBuilder::new();
    let mut drawn = 0usize;
    for index in packed_cells.iter().step_by(stride) {
        let i = *index as usize;
        let cell = [i % c, (i / c) % c, i / (c * c)];
        let lo = origin
            + Vec3::new(cell[0] as f32, cell[1] as f32, cell[2] as f32) * h
            + Vec3::splat(inset);
        push_cube(&mut cloud, lo, fill);
        drawn += 1;
    }

    let active = packed_cells.len();
    let fraction = 100.0 * active as f64 / total_cells as f64;
    let ratio = if packed_ns > 0.0 {
        scalar_ns / packed_ns
    } else {
        0.0
    };

    // The HUD is the evidence and a headless capture has no HUD to read. One
    // line per rebuild, so `ISOMESH_CAPTURE` leaves the measurement in the log
    // where a script can hold it against M-337.
    info!(
        "{} at {samples}^3: {active} of {total_cells} cells active ({fraction:.4}%); \
         stage {scalar_ns:.4} -> {packed_ns:.4} ns/cell = {ratio:.3}x; \
         same ordered list {lists_match}; {} exact zeros ({} of them -0.0), \
         sign-bit variant differs on {sign_bit_diff} cells; \
         sample {sample_ms:.1} ms, extract {extract_ms:.1} ms, \
         {} vertices, cloud stride {stride}",
        F::NAME,
        zeros.exact,
        zeros.negative,
        buffer.vertex_count(),
    );

    Some(Built {
        surface: to_mesh(&buffer),
        cloud: finish(cloud),
        report: Report {
            field_name: F::NAME,
            samples,
            cells,
            total_cells,
            bit_row: grid.bit_row,
            active,
            fraction,
            scalar_ns,
            packed_ns,
            ratio,
            lists_match,
            zeros,
            sign_bit_diff,
            sample_ms,
            words,
            origin,
            cell_size: h,
            cloud_stride: stride,
            cloud_drawn: drawn,
        },
        domain_min: Vec3::new(min[0] as f32, min[1] as f32, min[2] as f32),
        domain_max: Vec3::new(max[0] as f32, max[1] as f32, max[2] as f32),
        vertices: buffer.vertex_count(),
        triangles: buffer.triangle_count(),
        extract_ms,
    })
}

/// One axis-aligned cube, six quads, outward normals.
///
/// Merged into one mesh rather than instanced because Bevy 0.19 has no
/// instancing path a `StandardMaterial` can take, and one entity per cell would
/// be 40,000 entities to cull every frame.
fn push_cube(builder: &mut MeshBuilder, min: Vec3, size: f32) {
    for axis in 0..3usize {
        // `(axis, axis+1, axis+2)` is a right-handed cyclic triple, so `u x v`
        // is `+axis` and the corner order below is counter-clockwise seen from
        // `+axis` — which is the winding Bevy takes as front-facing.
        let u = (axis + 1) % 3;
        let v = (axis + 2) % 3;
        for side in 0..2usize {
            let mut normal = [0.0f32; 3];
            normal[axis] = if side == 1 { 1.0 } else { -1.0 };
            let mut quad = [0u32; 4];
            for (corner, slot) in quad.iter_mut().enumerate() {
                let (a, b) = match corner {
                    0 => (0.0, 0.0),
                    1 => (1.0, 0.0),
                    2 => (1.0, 1.0),
                    _ => (0.0, 1.0),
                };
                let mut p = [min.x, min.y, min.z];
                if side == 1 {
                    p[axis] += size;
                }
                p[u] += a * size;
                p[v] += b * size;
                *slot = builder.vertex(p, normal);
            }
            if side == 1 {
                builder.triangle(quad[0], quad[1], quad[2]);
                builder.triangle(quad[0], quad[2], quad[3]);
            } else {
                builder.triangle(quad[0], quad[2], quad[1]);
                builder.triangle(quad[0], quad[3], quad[2]);
            }
        }
    }
}

/// The `f64` extraction as a Bevy mesh, or `None` when it produced nothing.
///
/// Cast rather than re-extracted in `f32`: the numbers on the HUD are `f64`
/// numbers, so the mesh the picture is drawn from has to be the one they were
/// computed alongside.
fn to_mesh(buffer: &MeshBuffer<f64>) -> Option<Mesh> {
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
    finish(builder)
}

/// A builder's mesh, or `None` when nothing was written into it. See [`mesh3d`].
fn finish(builder: MeshBuilder) -> Option<Mesh> {
    if builder.vertex_count() == 0 {
        return None;
    }
    Some(builder.into_mesh())
}

/// Walk the cursor along the words that have something in them.
///
/// Frame-counter driven under capture and clock driven otherwise, so a GIF scans
/// the shell rather than photographing one row of it. `Space` freezes it, which
/// is what makes a paused capture inspectable.
fn advance_cursor(
    time: Res<Time>,
    capture: Res<Capture>,
    keys: Res<ButtonInput<KeyCode>>,
    flags: Res<ViewFlags>,
    report: Res<Report>,
    mut cursor: ResMut<Cursor>,
    mut last_taken: Local<u32>,
) {
    let len = report.words.len();
    if len == 0 {
        cursor.index = 0;
        return;
    }
    if cursor.index >= len {
        cursor.index = cursor_start(len);
    }

    if keys.just_pressed(KeyCode::Period) {
        cursor.index = (cursor.index + 1) % len;
    }
    if keys.just_pressed(KeyCode::Comma) {
        cursor.index = (cursor.index + len - 1) % len;
    }
    if flags.paused {
        return;
    }

    let step = cursor_step(len);
    if capture.is_active() {
        if capture.taken != *last_taken {
            *last_taken = capture.taken;
            cursor.index = (cursor.index + step) % len;
        }
        return;
    }
    cursor.timer += time.delta_secs();
    while cursor.timer >= CURSOR_INTERVAL {
        cursor.timer -= CURSOR_INTERVAL;
        cursor.index = (cursor.index + step) % len;
    }
}

/// Frame the domain, with the subject off centre.
///
/// The offset is applied in the camera's own basis, from the same yaw and pitch
/// the harness builds its transform from, so it is one screen-space nudge
/// however far `ISOMESH_SPIN` has turned.
fn frame_camera(
    report: Res<Report>,
    mut camera: Query<&mut OrbitCamera>,
    domain: Query<&DemoDomain>,
) {
    let Some(d) = domain.iter().next() else {
        return;
    };
    if report.cells == 0 {
        return;
    }
    let centre = (d.min + d.max) * 0.5;
    let radius = (d.max.x - d.min.x) * VIEW_RADIUS_EXTENTS;
    for mut orbit in &mut camera {
        let dir = Vec3::new(
            orbit.yaw.cos() * orbit.pitch.cos(),
            orbit.pitch.sin(),
            orbit.yaw.sin() * orbit.pitch.cos(),
        );
        // `orbit_camera` puts the eye at `focus + dir * radius`, so the view
        // direction is `-dir` and a focus moved along `-right` puts the subject
        // right of centre.
        let forward = -dir;
        let right = forward.cross(Vec3::Y).normalize_or_zero();
        let up = right.cross(forward).normalize_or_zero();
        orbit.focus =
            centre - right * (SUBJECT_OFFSET.x * radius) + up * (SUBJECT_OFFSET.y * radius);
        orbit.radius = radius;
    }
}

/// Paint the strip from the word under the cursor.
///
/// Colours are written only when they differ. An unconditional write marks 448
/// `BackgroundColor`s changed every frame, and Bevy's UI extraction is
/// change-driven, so it would turn a static panel into per-frame work.
fn update_strip(
    report: Res<Report>,
    cursor: Res<Cursor>,
    flags: Res<ViewFlags>,
    mut cells: Query<(&BitCell, &mut BackgroundColor)>,
    mut footer: Query<&mut Text, With<StripFooter>>,
    mut panels: Query<&mut Visibility, With<StripPanel>>,
) {
    let wanted = if flags.hud {
        Visibility::Inherited
    } else {
        Visibility::Hidden
    };
    for mut visible in &mut panels {
        if *visible != wanted {
            *visible = wanted;
        }
    }

    let word = report.words.get(cursor.index).copied().unwrap_or_default();
    for (cell, mut colour) in &mut cells {
        let bit = 1u64 << cell.col;
        let target = match cell.row {
            // Masked by the *sample* extent, not the cell extent: bit 63 of the
            // last word is a real sample that no cell in this word owns, and it
            // is exactly the one `inside_word_shifted` reaches for. Painting it
            // as an unlit cell would hide the reason that function exists.
            0..=3 => lit(
                word.fused.rows[cell.row] & bit != 0,
                word.samples & bit != 0,
                BIT_INSIDE,
            ),
            4 => lit(word.fused.any & bit != 0, word.cells & bit != 0, BIT_ANY),
            5 => lit(word.fused.all & bit != 0, word.cells & bit != 0, BIT_ALL),
            _ => lit(
                word.fused.active & bit != 0,
                word.cells & bit != 0,
                BIT_ACTIVE,
            ),
        };
        if colour.0 != target {
            colour.0 = target;
        }
    }

    let base = word.w as usize * WORD_BITS;
    let cells_here = word.cells.count_ones() as usize;
    let text = if report.words.is_empty() {
        format!(
            "no active cells on {} at {}^3 -- the surface misses this domain entirely",
            report.field_name, report.samples
        )
    } else {
        format!(
            "word {}/{}  y={} z={}  bits 0-{} samples, 0-{} cells = x {}-{}  {} active  row {}/{}",
            word.w + 1,
            report.bit_row,
            word.y,
            word.z,
            word.samples.count_ones().saturating_sub(1),
            cells_here.saturating_sub(1),
            base,
            base + cells_here.saturating_sub(1),
            (word.fused.active & word.cells).count_ones(),
            cursor.index + 1,
            report.words.len(),
        )
    };
    for mut target in &mut footer {
        if target.0 != text {
            target.0.clone_from(&text);
        }
    }
}

/// A bit's colour: lit, unlit, or not in the grid at all.
///
/// The three-way split is the whole reason 33 samples is on the ladder — half of
/// that word is answers about cells that do not exist, and a strip that drew them
/// the same dark grey as an unlit cell would be quietly lying about which bits
/// `cell_mask` keeps.
fn lit(set: bool, exists: bool, colour: Color) -> Color {
    if !exists {
        BIT_VOID
    } else if set {
        colour
    } else {
        BIT_OFF
    }
}

/// The fraction, and the counts under it.
fn update_big_number(
    report: Res<Report>,
    flags: Res<ViewFlags>,
    mut fraction: Query<&mut Text, (With<BigFraction>, Without<BigCount>)>,
    mut count: Query<&mut Text, (With<BigCount>, Without<BigFraction>)>,
) {
    if !flags.hud {
        return;
    }
    let big = format!("{:.3}% active", report.fraction);
    for mut target in &mut fraction {
        if target.0 != big {
            target.0.clone_from(&big);
        }
    }
    let small = format!(
        "{} of {} cells   {}^3 samples",
        thousands(report.active),
        thousands(report.total_cells),
        report.samples
    );
    for mut target in &mut count {
        if target.0 != small {
            target.0.clone_from(&small);
        }
    }
}

/// `1234567` as `1,234,567`.
///
/// The counts here run to eight digits and the argument is that one of them is
/// falling as a fraction while the other climbs; unseparated digit runs make
/// that comparison guesswork.
fn thousands(value: usize) -> String {
    let digits = value.to_string();
    let mut out = String::with_capacity(digits.len() + digits.len() / 3);
    for (i, ch) in digits.chars().enumerate() {
        if i > 0 && (digits.len() - i).is_multiple_of(3) {
            out.push(',');
        }
        out.push(ch);
    }
    out
}

/// The HUD. The numbers are the demo.
fn report(report: Res<Report>, mut stats: ResMut<DemoStats>) {
    // Kept short on purpose. The shared HUD's title is one unconstrained `Text`
    // node, so it wraps at the window width and a wrapped title pushes every
    // line below it down into the geometry -- measured at 836 px wide, which is
    // what this display grants when `ISOMESH_WINDOW` is not honoured.
    stats.title = format!(
        "E-307  active cells - {}   {}^3   [1-5] field  [ ] res  , . row  H surface",
        report.field_name, report.samples,
    );
    stats.extra = vec![
        format!(
            "{:>9} cells/axis   {} cells   {} u64 word{} per sample row   {:.1} ms to sample",
            report.cells,
            thousands(report.total_cells),
            report.bit_row,
            if report.bit_row == 1 { "" } else { "s" },
            report.sample_ms,
        ),
        format!(
            "{:>9} active cells = {:.3}% of the grid",
            thousands(report.active),
            report.fraction,
        ),
        String::from("          M-337 sphere: 1.89% / 0.93% / 0.46% at 64 / 128 / 256 samples"),
        String::new(),
        format!(
            "{:>9.3} ns/cell  eight-corner scalar gather",
            report.scalar_ns
        ),
        format!(
            "{:>9.3} ns/cell  packed u64 word test, bitmap build included",
            report.packed_ns
        ),
        format!(
            "{:>8.2}x  faster   M-337's band for this stage is 5.0-5.5x",
            report.ratio
        ),
        format!(
            "          same ordered list: {}   ({} cells, element for element)",
            if report.lists_match { "YES" } else { "NO" },
            thousands(report.active),
        ),
        String::new(),
        String::from("          inside is `value < 0`, NOT the IEEE sign bit: -0.0 has the"),
        format!(
            "          bit set and -0.0 < 0.0 is false. {} samples here are exactly",
            thousands(report.zeros.exact)
        ),
        format!(
            "          zero, {} of them -0.0; sign-bit variant differs on {} cells",
            thousands(report.zeros.negative),
            thousands(report.sign_bit_diff),
        ),
        format!(
            "          cloud: {} cubes drawn, 1 of every {}",
            thousands(report.cloud_drawn),
            report.cloud_stride,
        ),
    ];
}

/// The highlighted word: one long box round its 64 cells, one cage per active
/// cell inside it.
fn draw_row(
    report: Res<Report>,
    cursor: Res<Cursor>,
    surface: Res<ShowSurface>,
    mut visibility: Query<&mut Visibility, (With<DemoMesh>, Without<StripPanel>)>,
    mut gizmos: Gizmos<RowGizmos>,
) {
    let wanted = if surface.0 {
        Visibility::Inherited
    } else {
        Visibility::Hidden
    };
    for mut visible in &mut visibility {
        if *visible != wanted {
            *visible = wanted;
        }
    }

    let Some(word) = report.words.get(cursor.index) else {
        return;
    };
    let h = report.cell_size;
    let base = word.w as usize * WORD_BITS;
    let span = (report.cells as usize).saturating_sub(base).min(WORD_BITS);
    if span == 0 || h <= 0.0 {
        return;
    }

    let row_lo = report.origin + Vec3::new(base as f32, word.y as f32, word.z as f32) * h;
    box_edges(
        &mut gizmos,
        row_lo,
        row_lo + Vec3::new(span as f32 * h, h, h),
        WORD_COLOUR,
    );

    let mut bits = word.fused.active & word.cells;
    while bits != 0 {
        let k = bits.trailing_zeros() as usize;
        bits &= bits - 1;
        let lo = report.origin + Vec3::new((base + k) as f32, word.y as f32, word.z as f32) * h;
        box_edges(&mut gizmos, lo, lo + Vec3::splat(h), BIT_ACTIVE);
    }
}

/// The twelve edges of a box, at its exact bounds.
///
/// Exact rather than inflated: a cage larger than its cell would make the shell
/// look thicker than it is, which is the one thing a picture of a *fraction*
/// must not do.
fn box_edges(gizmos: &mut Gizmos<RowGizmos>, lo: Vec3, hi: Vec3, colour: Color) {
    // Corner indexing matches the extractor's: bit `i` of the index is axis `i`.
    let corner = |i: usize| {
        Vec3::new(
            if i & 1 == 0 { lo.x } else { hi.x },
            if i & 2 == 0 { lo.y } else { hi.y },
            if i & 4 == 0 { lo.z } else { hi.z },
        )
    };
    for i in 0..8usize {
        for axis in 0..3usize {
            let bit = 1 << axis;
            if i & bit == 0 {
                gizmos.line(corner(i), corner(i | bit), colour);
            }
        }
    }
}
