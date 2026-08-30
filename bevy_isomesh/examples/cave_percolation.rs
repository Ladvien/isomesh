//! E-320 — a generated cave's air is a hundred sealed pockets or one giant
//! component, and the isovalue is what decides which.
//!
//! ```bash
//! cd bevy_isomesh && cargo run --example cave_percolation --release
//! ```
//!
//! Keys: `1` `2` field, `]` digs and `[` fills (one rung of the isovalue ladder
//! each), `S` the 2D slice panel, `W` wireframe, `N` normals, `G` domain box,
//! `H` HUD, `Space` freezes the spin, `F12` screenshot.
//!
//! **Always `--release`.** Switching field or resolution censuses the *whole*
//! ladder — 42 `Air::build` calls plus 42 slice censuses — because the
//! transition isovalue is a property of the sweep and not of one rung. Measured
//! at the default `65^3` on a 5900X: that pre-pass is **123 ms** on
//! `fbm_terrain` and **37 ms** on `noise_cavity`, and one rung after it costs
//! **60 ms** and **22 ms** (an extraction plus both instruments plus the slice).
//! A debug build makes each of those a stall rather than a hitch.
//!
//! Under `ISOMESH_CAPTURE` it needs no keyboard: the isovalue walks the ladder
//! down and back up on every captured frame, so the clip is the giant component
//! appearing and dissolving rather than a still. `ISOMESH_FIELD` pins the field,
//! `ISOMESH_SAMPLES` pins the resolution, `ISOMESH_PANEL=off` starts with the
//! slice panel hidden.
//!
//! ```bash
//! # 82 frames is exactly one down-and-up pass of the 42-rung ladder.
//! # Two clips, because the finding is the contrast between the two fields.
//! # Size not measured here: this host has no display. `record_gif.sh` warns
//! # outside 0.7-4.8 MB.
//! ISOMESH_FIELD=1 ISOMESH_CAPTURE_FRAMES=82 ISOMESH_CAPTURE_EVERY=2 \
//!   scripts/record_gif.sh cave_percolation docs/gifs/e320.gif
//! ISOMESH_FIELD=0 ISOMESH_CAPTURE_FRAMES=82 ISOMESH_CAPTURE_EVERY=2 \
//!   scripts/record_gif.sh cave_percolation docs/gifs/e320-terrain.gif
//! ```
//!
//! Demonstrates **P-176**, and demonstrates the half of it that was falsified as
//! plainly as the half that held.
//!
//! # What it reads, live, at the default `65^3`
//!
//! Recomputed through the public API rather than quoted, and put beside the
//! committed `97^3` row in the panel so the two can be compared:
//!
//! | | `noise_cavity` | `fbm_terrain` |
//! |---|---|---|
//! | 3D onset (live / P-176) | `0.000000` / `0.029340` | `7.981641` / `8.187580` |
//! | 2D onset (live / P-176) | **`none` / `none`** | `7.030127` / `7.272575` |
//! | largest at `iso = 0` (live / P-176) | `0.703597` / `0.814493` | `1.000000` / `1.000000` |
//! | fragmented rungs / giant rungs | 30 / 5 | 1 / 37 |
//! | at `iso = 0`: 3D | 82 components, largest `0.7036` | 1 component |
//! | at `iso = 0`: the slice | 31 components, largest `0.1434` | 1 component |
//!
//! **That bottom pair of rows is the whole finding.** On `noise_cavity` at the
//! isovalue the registration asked about, three dimensions hold seven tenths of
//! the air in one component while the slice through the same field at the same
//! level holds one seventh in its largest of thirty-one. On `fbm_terrain` both
//! are one component, and no comparison is left to make.
//!
//! # What is on screen
//!
//! - **The extracted surface** — Marching Cubes on `f - iso`, so the solid is
//!   `{f < iso}` and the air is `{f >= iso}`. That subtraction is the whole
//!   sweep: `Air::build` takes air to be `value >= 0`, so handing it `f - iso`
//!   asks about the excursion set `{f >= iso}` with no reinterpretation.
//! - **Its colour is the air component behind it.** Every vertex is coloured by
//!   the label of the nearest air sample, which for a Marching Cubes vertex is
//!   the air end of its own grid edge. **Gold is always the largest component**,
//!   whatever its size, so the transition is gold taking over the screen; the
//!   other components cycle a six-colour palette by label. Grey is a vertex with
//!   no air sample beside it at all.
//! - **The panel to the right** — the same field at the same isovalue on one
//!   `z` plane, censused with **4-connectivity**, one quad per air pixel,
//!   coloured the same way. `S` hides it. The cyan rectangle inside the volume
//!   is where that plane was cut from.
//! - **The `sweep` / `giant` / `2D` rows in the HUD** — the whole ladder at
//!   once. `sweep` ramps the largest component's share of the air from `.` to
//!   `@`, `giant` marks with `G` every rung where one component holds more than
//!   half the air, `2D` marks the same test on the slice with `g`, and `^` is the
//!   rung on screen. High isovalue on the left, low on the right.
//!
//! # What the numbers mean
//!
//! P-176 registered its clauses before the harness existed. Duminil-Copin,
//! Rivera, Rodriguez & Vanneuville (`arXiv:2108.08008` / `10.1214/22-aop1594`)
//! prove that for smooth Gaussian fields in `d >= 3` the critical level
//! `l_c(d)` is strictly positive: the excursion set percolates at levels where
//! the two-dimensional analogue, whose critical level is zero, does not.
//!
//! - **C1 — a giant component appears, and the isovalue at which it appears is
//!   reported. HELD.** "Giant" is one component holding **above half** the air
//!   volume. That is also true at the very top of the ladder where the air is a
//!   single voxel, so the reported onset is the **persistent** one: the highest
//!   rung at or below which the giant never disappears again. The HUD carries
//!   the live onset beside the committed one.
//! - **C2 — three dimensions differ qualitatively from a 2D slice of the same
//!   field. SPLIT 42 rows to 42.** The criterion is the theorem's own shape: 3D
//!   has a persistent giant phase **and** the slice has none anywhere in the
//!   swept range. That is what `noise_cavity` does. `fbm_terrain` does not, and
//!   the registration said in advance that it might not: it is hash-based
//!   lattice noise, not a Gaussian field, and the risk was registered as real.
//!   The mechanism is visible here — `fbm_terrain` is `y - h(x, z)`
//!   (`fields/mod.rs:1352-1361`), so its air set is `{y >= h(x, z) + iso}`, the
//!   region **above a graph**. That is one component in three dimensions and one
//!   in two, so no dimension gap can exist to be measured. **This demo shows
//!   both fields on purpose**: a percolation demo that only ever showed the field
//!   where percolation happens would misrepresent how often it happens.
//! - **C3 — the crate's `Air` agrees with an independent union-find. HELD on 84
//!   of 84 rows.** Recomputed live here rather than only cited, because the
//!   colours on screen *are* `Air`'s labels: a second instrument beside it is
//!   what makes the picture evidence. The union-find is deliberately the retired
//!   algorithm — three forward edges per air sample, union by size, path halving
//!   — while `Air` is a flat label array filled by breadth-first flood
//!   (`connectivity.rs:29-46`). The two share the six-neighbour adjacency and
//!   the `value >= 0` test, which is all C3 is entitled to assume.
//!
//! **6 in 3D and 4 in 2D is not this demo's choice.** `Air::neighbours`
//! (`connectivity.rs:610-646`) is the six axis-aligned neighbours, so the second
//! instrument must use six or the two would be comparing different graphs; the
//! slice takes the same rule restricted to a plane, which is four.
//!
//! # The region, and why `noise_cavity` is masked
//!
//! `fbm_terrain` is swept over its own `[-8, 8]^3` domain, unmasked.
//!
//! `noise_cavity` is `NoiseVolume ∩ Sphere{r: 1.5}` over `[-2, 2]^3`, and over
//! the whole box the question stops being about the field: outside the cap every
//! sample has `sphere > 0`, hence `max(noise, sphere) > 0`, hence air at every
//! isovalue at or below zero — and that shell is one connected object wrapping
//! the solid. P-176 measured it: `82.58%` air in 12 components whose largest
//! holds `99.99%`. C1 would then be true by the sampling box rather than by the
//! field. So the region is the field's **own** cap, `cap_sdf <= -0.05`, and the
//! sweep is floored at `iso = -0.05` so that inside the region the cap cannot
//! contribute — `sphere_sdf` is below the floor everywhere in the mask, so
//! `{max(noise, sphere) >= iso}` is exactly `{noise >= iso}`. The HUD reports
//! `cap max` inside the region and complains loudly if it ever reaches the floor.
//!
//! One expression defines that region — [`region_value`] — and it is what the
//! extractor is handed *and* what the census array is filled from, so the
//! picture and the answer cannot disagree.
//!
//! # What this does not show
//!
//! Not a cost measurement. `M-311 / P-23` measured `Air`'s repair cost; this is
//! about what the air region *looks like*. Not the theorem either: these fields
//! are not Gaussian, the grid is finite, and a lattice census of a finite box is
//! not a statement about `l_c(d)`. What it is: the shape the theorem predicts,
//! looked for in this crate's own generated worlds, found in one of them and
//! honestly absent from the other.

mod common;

use std::time::Instant;

use bevy::prelude::*;
use bevy_isomesh::MeshBuilder;
use common::{Capture, CommonPlugin, DemoDomain, DemoMesh, DemoStats, OrbitCamera, ViewFlags};
use isomesh::connectivity::Air;
use isomesh::fields::{FbmTerrain, ReferenceField, Sphere, noise_cavity};
use isomesh::marching_cubes::MarchingCubes;
use isomesh::{MeshBuffer, MeshSink, RuntimeShape3, Sdf};

// ─── the measurement's own constants, from `benches/experiment_p176.rs` ──────

/// Samples per axis when nothing pins it.
///
/// P-176 measured a single `97^3`, and one resolution rather than three because
/// percolation is a property of the field's feature count and not of the
/// sampling rate: refining resolves the same blobs better and does not add
/// blobs. 65 is the interactive rung of the same statement; `ISOMESH_SAMPLES=97`
/// reproduces the committed grid exactly.
const DEFAULT_SAMPLES: u32 = 65;

/// Below this there is no cave to speak of, and the ladder is one value.
const MIN_SAMPLES: u32 = 9;

/// Above this the whole-ladder pre-pass stops being a hitch and becomes a wait.
const MAX_SAMPLES: u32 = 129;

/// Evenly spaced rungs across a field's own sampled value range, as P-176 swept
/// them. `iso = 0` is inserted when the range straddles it, which both fields
/// do, giving 42 rungs each.
const RUNGS: usize = 41;

/// How far inside its own cap `noise_cavity` is measured, and the sweep floor.
const CAP_MARGIN: f64 = 0.05;

/// The cap's radius, from `fields/mod.rs:1242`.
const CAP_RADIUS: f64 = 1.5;

/// The registered giant-component threshold: *above* half the air volume.
const GIANT_SHARE: f64 = 0.5;

/// "Many" small components, for the fragmented-regime control.
const MANY_COMPONENTS: u64 = 8;

/// Air share of the region below which a giant component is a ladder artefact
/// rather than a regime — the top rung admits one voxel, and one voxel is
/// trivially all of its own component.
const REAL_AIR_SHARE: f64 = 0.05;

/// Air pixels a slice needs before its verdict is evidence about dimension.
const SLICE_FLOOR: u64 = 100;

/// How close to zero a rung must be for `iso = 0` to count as already present.
const ZERO_RUNG: f64 = 1e-12;

/// Component sizes written out individually in the HUD.
const TOP_SIZES: usize = 5;

/// Fields on offer. Both of P-176's arms, and the negative one is not optional.
const FIELD_COUNT: usize = 2;

// ─── the committed rows, quoted as citations ────────────────────────────────

/// One field's P-176 row, read from `docs/experiments/p-176.csv`.
///
/// Every field here is a column of that file, and every use of it names P-176 in
/// the HUD line it lands on. The live numbers beside them are recomputed through
/// the public API on whatever grid is on screen; these are the committed `97^3`
/// ones, so a reader can hold one against the other.
struct Cited {
    /// The `field` column, which is also `ReferenceField::NAME`.
    name: &'static str,
    /// `percolation_isovalue` — the refined persistent onset in 3D.
    onset_3d: f64,
    /// `percolation_isovalue_2d`, `None` where the CSV says `none`.
    onset_2d: Option<f64>,
    /// `c2_holds`, which is constant down each field's 42 rows.
    c2: bool,
    /// `largest_component_fraction_at_zero` — the registration's own sentence
    /// was about `iso = 0`, so the ladder contains it and the CSV answers it.
    largest_at_zero: f64,
    /// `percolation_rung_fragmented` and the `isovalue` / `components` /
    /// `largest_component_fraction` of that row: the last rung before the onset.
    fragmented: (f64, u64, f64),
    /// `percolation_rung_giant` and the same three columns: the first rung after
    /// it. The pair is the transition in four numbers.
    giant: (f64, u64, f64),
    /// `single_component_rows` — rungs whose air is one component.
    single_rows: u32,
}

/// The two rows, in field order. `p-176.csv:6` and `:48` onward.
const CITED: [Cited; FIELD_COUNT] = [
    Cited {
        name: "fbm_terrain",
        onset_3d: 8.187_580,
        onset_2d: Some(7.272_575),
        c2: false,
        largest_at_zero: 1.000_000,
        fragmented: (8.468_861, 20, 0.340_507),
        giant: (7.993_365, 18, 0.858_262),
        single_rows: 37,
    },
    Cited {
        name: "noise_cavity",
        onset_3d: 0.029_340,
        onset_2d: None,
        c2: true,
        largest_at_zero: 0.814_493,
        fragmented: (0.036_372, 100, 0.360_519),
        giant: (0.019_098, 81, 0.521_363),
        single_rows: 1,
    },
];

/// Rows in `p-176.csv`. 42 per field, and the denominator of the C2 split.
const CITED_ROWS: u32 = 84;

/// The resolution every cited number was measured at.
const CITED_RESOLUTION: u32 = 97;

/// `outer_shell_air_fraction`, `outer_shell_components` and
/// `outer_shell_largest_fraction`: what `noise_cavity` reads unmasked at
/// `iso = 0`, and the reason the cap mask below is load-bearing.
const CITED_SHELL: (f64, u64, f64) = (0.825_811, 12, 0.999_874);

// ─── colour ─────────────────────────────────────────────────────────────────

/// The largest air component, always, whatever share it holds.
///
/// Reserved rather than cycled, because the transition *is* this colour taking
/// over: a palette that renumbered the giant every rung would draw the same
/// event as a flicker.
const GIANT_COLOUR: [f32; 4] = [1.00, 0.78, 0.22, 1.0];

/// Every other component, by label. Six hues, none of them gold.
const POCKET_COLOURS: [[f32; 4]; 6] = [
    [0.25, 0.70, 0.95, 1.0],
    [0.95, 0.35, 0.45, 1.0],
    [0.45, 0.85, 0.45, 1.0],
    [0.80, 0.45, 0.95, 1.0],
    [0.35, 0.95, 0.85, 1.0],
    [0.95, 0.60, 0.30, 1.0],
];

/// A surface vertex with no air sample beside it: masked out, or a cell the air
/// test rejected on all eight corners.
const ORPHAN_COLOUR: [f32; 4] = [0.42, 0.44, 0.50, 1.0];

/// The ramp the `sweep` row is drawn with: the largest component's share of the
/// air, `.` at nothing and `@` at all of it. ASCII, because the HUD font draws
/// anything else as an empty box.
const RAMP: [char; 8] = ['.', ':', '-', '=', '+', '*', '#', '@'];

// ─── framing ────────────────────────────────────────────────────────────────

/// Where the slice panel sits, in domain widths along `+x`.
const PANEL_GAP: f32 = 1.18;

/// Orbit radius as a multiple of the volume-plus-panel width, which is
/// `1 + PANEL_GAP` domain widths plus the one cell [`panel_bounds`] explains —
/// `2.196 w` at `65^3`.
///
/// Derived rather than tuned by eye, because the machine this was written on has
/// no display. Bevy's `Camera3d` defaults to a 45-degree *vertical* field of
/// view, so at 16:9 the horizontal half-angle has
/// `tan = tan(22.5 deg) * 16/9 = 0.7365`. The subject's half-width is `1.098 w`
/// and [`SUBJECT_OFFSET`] pushes its far edge out by another `0.06 * radius`, so
/// `0.7365 r >= 1.098 w + 0.06 r`, i.e. `r >= 1.623 w`. `1.02 * 2.196 = 2.24 w`
/// clears it with room for the depth the yaw adds; the vertical constraint,
/// `0.4142 r >= 0.5 w`, is slack by a factor of four.
///
/// Not left to that arithmetic either: `the_camera_frames_the_volume_and_the_panel`
/// projects all twelve corners through it, and a `0.70` here fails it.
const VIEW_RADIUS_WIDTHS: f32 = 1.02;

/// Nudge, in the camera's own basis, so the HUD does not sit on the evidence.
///
/// Twenty-odd lines of panel occupy the upper left; centring the subject
/// photographs the argument with its evidence hidden, which is E-112's lesson.
/// Small here rather than `active_cells`' `0.17`, because the subject is two
/// objects side by side and already fills most of the frame's width.
const SUBJECT_OFFSET: Vec2 = Vec2::new(0.06, -0.06);

/// Yaw: nearly face-on down `-z`, off it by about twenty degrees.
///
/// Face-on would be right for the panel and wrong for the volume — a solid read
/// down an axis collapses into a silhouette. Twenty degrees costs the panel six
/// percent of its width and gives the volume back its depth. `dual_contouring_cube`
/// made the same trade in the other direction (`BACKLOG_ARCHIVE.md:178`).
const VIEW_YAW: f32 = std::f32::consts::FRAC_PI_2 + 0.34;

/// Pitch: above the horizon, so a heightfield's air reads as *above* something.
const VIEW_PITCH: f32 = 0.26;

// ─── resources ──────────────────────────────────────────────────────────────

/// Which field, mirrored from [`ViewFlags::field`] so the digit keys work.
#[derive(Resource)]
struct Field(usize);

/// Which rung of the ladder is on screen. Index 0 is the highest isovalue.
#[derive(Resource)]
struct Cursor(usize);

/// A resolution pinned by `ISOMESH_SAMPLES`, which takes the default out of play.
///
/// Clamped: every cost here is cubic, and the whole-ladder pre-pass multiplies
/// it by 42.
#[derive(Resource)]
struct Pinned(Option<u32>);

/// Whether the 2D slice panel is drawn.
#[derive(Resource)]
struct ShowPanel(bool);

impl Default for ShowPanel {
    /// `ISOMESH_PANEL=off` starts with it hidden, so a clip of the volume alone
    /// can be recorded without a keyboard — the same contract `ISOMESH_VIEW` and
    /// `ISOMESH_FIELD` offer, for the same reason: a view reachable only by
    /// pressing a key is a view no committed image can be regenerated from.
    fn default() -> Self {
        Self(!matches!(
            std::env::var("ISOMESH_PANEL")
                .unwrap_or_default()
                .trim()
                .to_ascii_lowercase()
                .as_str(),
            "off" | "0" | "no" | "hide"
        ))
    }
}

/// The surface's material. White, because the vertex colours are the readout.
#[derive(Resource)]
struct SurfaceMaterial(Handle<StandardMaterial>);

/// The slice panel's material. Unlit, so a flat quad facing away from the light
/// still shows the colour it was given.
#[derive(Resource)]
struct PanelMaterial(Handle<StandardMaterial>);

/// The extracted surface.
#[derive(Component)]
struct SurfaceMesh;

/// The slice panel.
///
/// Deliberately **not** [`DemoMesh`]: `W` reads a marked mesh back and submits
/// three gizmo lines per triangle, and this one is two triangles per air pixel.
#[derive(Component)]
struct PanelMesh;

/// The slice plane's outline gets its own group so it draws in front of the
/// surface without dragging the shared wireframe forward with it.
#[derive(Default, Reflect, GizmoConfigGroup)]
struct SliceGizmos;

// ─── the sampled field ──────────────────────────────────────────────────────

/// One rung of the ladder, from the pre-pass.
///
/// Both dimensions, both censused with the crate's rule. The whole ladder is
/// held because the onset is a property of the sweep: "giant at this isovalue"
/// is true at the top of every ladder, where the air is one voxel.
struct Rung {
    /// The isovalue.
    iso: f64,
    /// Air components in 3D, from `Air::components`.
    components: u64,
    /// Air samples in 3D, from `Air::air_samples`.
    air: u64,
    /// The largest component's share of that air.
    largest_fraction: f64,
    /// One component holds above half the air.
    giant: bool,
    /// Air pixels on the slice. P-176's `SLICE_FLOOR` control reads this: a
    /// slice with nothing in it agrees with any verdict about dimension.
    air_2d: u64,
    /// Same test, on the slice.
    giant_2d: bool,
}

/// One field, sampled once, with the region its sweep is confined to and the
/// whole ladder already censused.
///
/// Rebuilt only when the field or the resolution changes; a rung change reuses
/// all of it.
struct Sampled {
    /// Index into [`CITED`], which is also the field key.
    field: usize,
    /// `ReferenceField::NAME`.
    name: &'static str,
    /// Samples per axis, all three equal.
    samples: u32,
    /// The shape `Air::build` and the extractor are handed.
    shape: RuntimeShape3,
    /// World position of sample `[0, 0, 0]`.
    origin: [f64; 3],
    /// Spacing between adjacent samples.
    cell_size: f64,
    /// The sampling box, for `DemoDomain`.
    domain_min: Vec3,
    /// The sampling box, for `DemoDomain`.
    domain_max: Vec3,
    /// The field's own values, `x` fastest. Unshifted.
    base: Vec<f64>,
    /// The cap's values, empty when the field is not capped.
    cap: Vec<f64>,
    /// The cap, for the extractor. `None` when the field is not capped.
    cap_field: Option<Sphere<f64>>,
    /// Samples the sweep may call air at all — `cap_term >= 0`, counted rather
    /// than kept, because [`region_value`] applies the same term itself. This is
    /// the denominator of every air share.
    region: u64,
    /// Human-readable statement of what the region is.
    mask_rule: &'static str,
    /// Largest cap value inside the region, or `-inf` where there is no cap.
    cap_max: f64,
    /// The sweep's floor and ceiling.
    sweep_lo: f64,
    /// The sweep's floor and ceiling.
    sweep_hi: f64,
    /// The isovalues, descending.
    ladder: Vec<f64>,
    /// One entry per rung.
    profile: Vec<Rung>,
    /// The persistent onset in 3D, as a ladder index.
    onset_3d: Option<usize>,
    /// The persistent onset on the slice, as a ladder index.
    onset_2d: Option<usize>,
    /// Rungs with at least [`MANY_COMPONENTS`] components and no giant one.
    fragmented_rows: usize,
    /// Rungs with a giant component over [`REAL_AIR_SHARE`] of the region.
    giant_rows: usize,
    /// Rungs whose slice holds at least [`SLICE_FLOOR`] air pixels.
    ///
    /// Without one, C2 compares a 3D census against an empty plane and its
    /// verdict is a default rather than a measurement (P-176's own control).
    slice_rows: usize,
    /// `z` index of the censused plane.
    slice_z: u32,
    /// Wall clock for the whole pre-pass.
    pre_pass_ms: f64,
}

/// The sampled field, or nothing yet.
#[derive(Resource, Default)]
struct Grid(Option<Sampled>);

// ─── the readout ────────────────────────────────────────────────────────────

/// Everything one rung's census produced, and the numbers the HUD reads.
#[derive(Resource, Default)]
struct Census {
    /// Zero until the first census; [`report`] returns while it is.
    samples: u32,
    /// Index into [`CITED`].
    field: usize,
    /// `ReferenceField::NAME`.
    name: &'static str,
    /// The rung on screen, and how many there are.
    rung: usize,
    /// The rung on screen, and how many there are.
    rungs: usize,
    /// Its isovalue.
    iso: f64,
    /// Air components, from `Air::components`.
    components: u64,
    /// Air samples, from `Air::air_samples`.
    air: u64,
    /// Air samples in the largest component.
    largest: u32,
    /// Its share of the air.
    largest_fraction: f64,
    /// The largest [`TOP_SIZES`] component sizes, descending.
    top_sizes: Vec<u32>,
    /// The independent union-find's component count.
    uf_components: u64,
    /// Its air-sample count.
    uf_air: u64,
    /// Its sorted size multiset equals `Air`'s.
    uf_sizes_match: bool,
    /// `Air::label_count` equals `Air::components`, so no retired label is
    /// hiding a component.
    labels_tight: bool,
    /// Components on the slice.
    components_2d: u64,
    /// Air pixels on the slice.
    air_2d: u64,
    /// The largest slice component's share of them.
    largest_fraction_2d: f64,
    /// Time the census of this rung took, both instruments and both dimensions.
    census_ms: f64,
    /// Time the extraction took.
    extract_ms: f64,
    /// Surface vertices and triangles.
    vertices: usize,
    /// Surface vertices and triangles.
    triangles: usize,
}

impl Census {
    /// C3 for this rung: same count, same air, same sizes, no retired label.
    fn agrees(&self) -> bool {
        self.uf_components == self.components
            && self.uf_air == self.air
            && self.uf_sizes_match
            && self.labels_tight
    }
}

// ─── the region, defined once ───────────────────────────────────────────────

/// The cap's contribution to the region: positive inside it, negative outside.
///
/// The region is `cap_sdf <= -0.05`. P-176's mask is written `cap_sdf < -0.05`,
/// and on the lattice the two are the same set — the CSV's `cap_max_sdf` is
/// `-0.050024`, so no sample lies between them. Written as a value rather than
/// as a predicate because the extractor needs a field, and one definition of
/// "the region" is worth more than an exact reproduction of an open interval.
fn cap_term(cap: f64) -> f64 {
    -(cap + CAP_MARGIN)
}

/// The value whose `>= 0` set is the air the census and the picture both mean.
///
/// `base - iso` is the excursion set `{f >= iso}`; the `min` with [`cap_term`]
/// intersects it with the region. **The single definition**: the extractor is
/// handed this through [`Excursion`] and the census array is filled from it, so
/// the surface on screen bounds exactly the components in the HUD.
fn region_value(base: f64, cap: Option<f64>, iso: f64) -> f64 {
    let air = base - iso;
    match cap {
        Some(c) => air.min(cap_term(c)),
        None => air,
    }
}

/// [`region_value`] as a field, for the extractor.
struct Excursion<'f, F> {
    /// The field being swept.
    field: &'f F,
    /// The level.
    iso: f64,
    /// The cap, when the field is measured inside one.
    cap: Option<Sphere<f64>>,
}

impl<F: Sdf<Scalar = f64>> Sdf for Excursion<'_, F> {
    type Scalar = f64;

    #[inline]
    fn sample(&self, p: [f64; 3]) -> f64 {
        region_value(
            self.field.sample(p),
            self.cap.map(|c| c.sample(p)),
            self.iso,
        )
    }

    #[inline]
    fn gradient(&self, p: [f64; 3]) -> [f64; 3] {
        // Forwarded to whichever branch of the `min` is active, and forwarded at
        // all because M-196 is about exactly this: a private field copy that
        // implements only `sample` loses the analytic gradient and pays seven
        // evaluations per normal. Both fields here have one, and a level shift
        // does not move it.
        if let Some(cap) = self.cap
            && cap_term(cap.sample(p)) < self.field.sample(p) - self.iso
        {
            let g = cap.gradient(p);
            return [-g[0], -g[1], -g[2]];
        }
        self.field.gradient(p)
    }
}

// ─── instrument A: an independent union-find ────────────────────────────────

/// A three-forward-edge union-find over the air set, union by size with path
/// halving.
///
/// P-176's instrument A, and the reason this demo's colours are evidence rather
/// than decoration. Deliberately the *retired* algorithm: `Air` is no longer a
/// union-find at all, so the two share no data structure, no traversal order and
/// no merge rule. They share the six-neighbour adjacency and the `value >= 0`
/// test, which is what C3 is entitled to assume.
///
/// Reused across rungs; `reset` is two `fill`s, so a sweep allocates twice.
struct Uf {
    /// `u32::MAX` where the sample is not air; otherwise the parent index.
    parent: Vec<u32>,
    /// Set size, meaningful only at a root.
    size: Vec<u32>,
}

impl Uf {
    /// One slot per sample, all dead.
    fn with(n: usize) -> Self {
        Self {
            parent: vec![u32::MAX; n],
            size: vec![0; n],
        }
    }

    /// Back to all dead.
    fn reset(&mut self) {
        self.parent.fill(u32::MAX);
        self.size.fill(0);
    }

    /// Whether slot `i` has been made air on this rung.
    fn live(&self, i: usize) -> bool {
        self.parent.get(i).copied().unwrap_or(u32::MAX) != u32::MAX
    }

    /// Start a singleton at `i`.
    fn make(&mut self, i: usize) {
        if let (Some(parent), Some(size)) = (self.parent.get_mut(i), self.size.get_mut(i)) {
            *parent = i as u32;
            *size = 1;
        }
    }

    /// Root of `i`, halving the path on the way. Only ever called on a live slot.
    fn find(&mut self, mut i: u32) -> u32 {
        while self.parent.get(i as usize).copied().unwrap_or(i) != i {
            let parent = self.parent.get(i as usize).copied().unwrap_or(i);
            let grand = self.parent.get(parent as usize).copied().unwrap_or(parent);
            if let Some(slot) = self.parent.get_mut(i as usize) {
                *slot = grand;
            }
            i = grand;
        }
        i
    }

    /// Join the sets of `a` and `b`, the larger set winning.
    fn union(&mut self, a: usize, b: usize) {
        let (mut ra, mut rb) = (self.find(a as u32), self.find(b as u32));
        if ra == rb {
            return;
        }
        let (sa, sb) = (self.root_size(ra), self.root_size(rb));
        if sa < sb {
            core::mem::swap(&mut ra, &mut rb);
        }
        if let Some(slot) = self.parent.get_mut(rb as usize) {
            *slot = ra;
        }
        if let Some(slot) = self.size.get_mut(ra as usize) {
            // Every set is a subset of the lattice, so the sum cannot exceed the
            // sample count and `u32` holds it. Saturating anyway, because an
            // arithmetic panic in a demo is a worse diagnostic than a wrong count
            // that the C3 comparison would then catch.
            *slot = sa.saturating_add(sb);
        }
    }

    /// Size at a root, or zero.
    fn root_size(&self, root: u32) -> u32 {
        self.size.get(root as usize).copied().unwrap_or(0)
    }

    /// Air count and root sizes over the first `n` slots, descending.
    fn harvest(&self, n: usize) -> (u64, Vec<u32>) {
        let mut sizes = Vec::new();
        let mut air = 0u64;
        for (i, &p) in self.parent.iter().take(n).enumerate() {
            if p != u32::MAX {
                air += 1;
            }
            if p == i as u32 {
                sizes.push(self.root_size(i as u32));
            }
        }
        sizes.sort_unstable_by(|a, b| b.cmp(a));
        (air, sizes)
    }

    /// Per-slot component rank over the first `n` slots — 0 for the largest
    /// component, `u32::MAX` where the slot is not air.
    ///
    /// Ranks rather than roots, so the palette is stable in size order rather
    /// than in whatever index the flood happened to reach first.
    fn ranks(&mut self, n: usize) -> Vec<u32> {
        let mut roots: Vec<(u32, u32)> = (0..n)
            .filter(|i| self.parent.get(*i).copied() == Some(*i as u32))
            .map(|i| (self.root_size(i as u32), i as u32))
            .collect();
        // Descending by size, ties by index, so the ranking is deterministic.
        roots.sort_unstable_by(|a, b| b.0.cmp(&a.0).then(a.1.cmp(&b.1)));
        let mut rank_of = vec![u32::MAX; self.parent.len()];
        for (rank, (_, root)) in roots.iter().enumerate() {
            if let Some(slot) = rank_of.get_mut(*root as usize) {
                *slot = rank as u32;
            }
        }
        let mut out = vec![u32::MAX; n];
        for i in 0..n {
            if !self.live(i) {
                continue;
            }
            let root = self.find(i as u32);
            if let (Some(slot), Some(rank)) = (out.get_mut(i), rank_of.get(root as usize)) {
                *slot = *rank;
            }
        }
        out
    }
}

// ─── the censuses ───────────────────────────────────────────────────────────

/// Instrument B: `Air::build` on the region array.
///
/// Returns `(Air, components, air samples, label count, sizes descending)`.
fn census_air(values: &[f64], shape: &RuntimeShape3) -> Option<(Air, u64, u64, u64, Vec<u32>)> {
    let built = match Air::build(values, shape) {
        Ok((built, _repair)) => built,
        Err(error) => {
            error!("E-320: Air::build over the census grid failed: {error}");
            return None;
        }
    };
    let mut sizes: Vec<u32> = (0..built.label_count() as u32)
        .map(|l| built.component_size(l))
        .filter(|s| *s > 0)
        .collect();
    sizes.sort_unstable_by(|a, b| b.cmp(a));
    let components = built.components();
    let air = built.air_samples();
    let labels = built.label_count() as u64;
    Some((built, components, air, labels, sizes))
}

/// Instrument A in 3D: 6-connectivity over `{region_value >= 0}`.
///
/// Scan in index order and, for each air sample, union it with its `-x`, `-y`
/// and `-z` neighbours where those are already air, so each lattice edge is seen
/// exactly once.
fn census_uf_3d(values: &[f64], dims: [u32; 3], uf: &mut Uf) -> (u64, Vec<u32>) {
    let [nx, ny, nz] = dims;
    let sy = nx as usize;
    let sz = (nx as usize) * (ny as usize);
    uf.reset();
    for k in 0..nz as usize {
        for j in 0..ny as usize {
            for i in 0..nx as usize {
                let idx = i + j * sy + k * sz;
                if values.get(idx).copied().unwrap_or(-1.0) < 0.0 {
                    continue;
                }
                uf.make(idx);
                if i > 0 && uf.live(idx - 1) {
                    uf.union(idx, idx - 1);
                }
                if j > 0 && uf.live(idx - sy) {
                    uf.union(idx, idx - sy);
                }
                if k > 0 && uf.live(idx - sz) {
                    uf.union(idx, idx - sz);
                }
            }
        }
    }
    uf.harvest(values.len())
}

/// One slice census.
struct Plane {
    /// Air pixels.
    air: u64,
    /// Components among them.
    components: u64,
    /// The largest one's share of the air.
    largest_fraction: f64,
    /// Per-pixel component rank, `u32::MAX` where solid. Row-major, `x` fastest.
    ranks: Vec<u32>,
}

/// Instrument A on one `z` plane: **4-connectivity**, the same air test.
///
/// Four rather than six because that is the six-neighbour rule restricted to a
/// plane. C2's whole content is the comparison, so the two graphs have to be the
/// same rule in different dimensions and not two separate choices.
fn census_2d(values: &[f64], dims: [u32; 3], z: u32, uf: &mut Uf) -> Plane {
    let [nx, ny, _] = dims;
    let plane = z as usize * (nx as usize) * (ny as usize);
    let pixels = (nx as usize) * (ny as usize);
    uf.reset();
    for j in 0..ny as usize {
        for i in 0..nx as usize {
            let local = i + j * nx as usize;
            if values.get(plane + local).copied().unwrap_or(-1.0) < 0.0 {
                continue;
            }
            uf.make(local);
            if i > 0 && uf.live(local - 1) {
                uf.union(local, local - 1);
            }
            if j > 0 && uf.live(local - nx as usize) {
                uf.union(local, local - nx as usize);
            }
        }
    }
    let (air, sizes) = uf.harvest(pixels);
    let largest = sizes.first().copied().unwrap_or(0);
    Plane {
        air,
        components: sizes.len() as u64,
        largest_fraction: if air == 0 {
            0.0
        } else {
            f64::from(largest) / air as f64
        },
        ranks: uf.ranks(pixels),
    }
}

/// 41 rungs across `[lo, hi]`, descending, with `iso = 0` inserted when the
/// range straddles it and no rung already lands there.
///
/// Spanning the field's own sampled range is what makes both regimes reachable
/// by construction: the top rung admits exactly the argmax sample and the bottom
/// rung admits everything the region allows.
fn ladder(lo: f64, hi: f64) -> Vec<f64> {
    let span = hi - lo;
    let last = (RUNGS - 1) as f64;
    let mut rungs: Vec<f64> = (0..RUNGS).map(|k| hi - span * (k as f64) / last).collect();
    if lo < 0.0 && hi > 0.0 && !rungs.iter().any(|v| v.abs() < ZERO_RUNG) {
        rungs.push(0.0);
    }
    rungs.sort_by(|a, b| b.total_cmp(a));
    rungs
}

/// The persistent onset: the lowest index from which `giant` holds on every
/// remaining rung. `None` when the bottom rung is not giant.
///
/// The persistent form rather than the first hit, because "one component holds
/// more than half the air" is true at the top of every ladder, where the air is
/// a single voxel and that voxel is all of its own component.
fn persistent_onset(profile: &[Rung], giant: impl Fn(&Rung) -> bool) -> Option<usize> {
    let mut onset = None;
    for i in (0..profile.len()).rev() {
        match profile.get(i) {
            Some(rung) if giant(rung) => onset = Some(i),
            _ => break,
        }
    }
    onset
}

// ─── sampling one field ─────────────────────────────────────────────────────

/// Sample `field` on an `n^3` lattice, `x` fastest.
fn sample_grid<S: Sdf<Scalar = f64>>(field: &S, lo: [f64; 3], h: f64, n: u32) -> Vec<f64> {
    let mut out = Vec::with_capacity((n as usize).pow(3));
    for k in 0..n {
        for j in 0..n {
            for i in 0..n {
                out.push(field.sample([
                    lo[0] + h * f64::from(i),
                    lo[1] + h * f64::from(j),
                    lo[2] + h * f64::from(k),
                ]));
            }
        }
    }
    out
}

/// Dispatch on the field index, then do the work once in [`sample_field`].
///
/// The reference fields are separate types, so a runtime choice has to be a
/// match rather than a loop over a list — the same shape `active_cells` and
/// `critical_cells` use.
fn sample(field: usize, samples: u32) -> Option<Sampled> {
    match field {
        0 => sample_field(field, &FbmTerrain::<f64>::canonical(), samples, None),
        _ => sample_field(
            1,
            &noise_cavity::<f64>(),
            samples,
            Some(Sphere::<f64> {
                center: [0.0; 3],
                radius: CAP_RADIUS,
            }),
        ),
    }
}

/// Sample one field, mask it, build its ladder, and census every rung.
///
/// The pre-pass is here rather than per frame because the onset is a property of
/// the whole sweep. It also puts the vacuity control on screen: a sweep that
/// visited only one regime would report an onset that is an artefact of where the
/// ladder stops.
fn sample_field<F>(
    index: usize,
    field: &F,
    samples: u32,
    cap: Option<Sphere<f64>>,
) -> Option<Sampled>
where
    F: Sdf<Scalar = f64> + ReferenceField,
{
    let (min, max) = field.domain();
    if samples < 2 {
        error!("E-320: {samples} samples per axis leaves no cell to march");
        return None;
    }
    let cell_size = (max[0] - min[0]) / f64::from(samples - 1);
    let shape = match RuntimeShape3::new([samples; 3]) {
        Ok(shape) => shape,
        Err(error) => {
            error!("E-320: grid {samples}^3 rejected: {error}");
            return None;
        }
    };

    let started = Instant::now();
    let base = sample_grid(field, min, cell_size, samples);
    if !base.iter().all(|v| v.is_finite()) {
        // `fbm_terrain` declares `FieldBound::Unbounded`, and `v - iso >= 0.0` is
        // false for a NaN — so a non-finite sample would read as cleanly empty
        // air instead of as broken. Said out loud rather than absorbed.
        error!("E-320: {} sampled a non-finite value", F::NAME);
        return None;
    }
    let cap_values = cap.map_or_else(Vec::new, |c| sample_grid(&c, min, cell_size, samples));
    let mask: Vec<bool> = if cap_values.is_empty() {
        vec![true; base.len()]
    } else {
        cap_values.iter().map(|c| cap_term(*c) >= 0.0).collect()
    };
    let region = mask.iter().filter(|m| **m).count() as u64;
    if region == 0 {
        error!(
            "E-320: {}'s measured region is empty at {samples}^3",
            F::NAME
        );
        return None;
    }
    let cap_max = cap_values
        .iter()
        .zip(&mask)
        .filter(|(_, m)| **m)
        .map(|(c, _)| *c)
        .fold(f64::NEG_INFINITY, f64::max);
    let in_region = || base.iter().zip(&mask).filter(|(_, m)| **m).map(|(v, _)| *v);
    let sweep_hi = in_region().fold(f64::NEG_INFINITY, f64::max);
    // The cavity's floor is its cap margin, not its minimum: below `-0.05` the
    // cap itself would be above the isovalue inside the region and the
    // components would be measuring the sphere rather than the cave.
    let sweep_lo = if cap.is_some() {
        -CAP_MARGIN
    } else {
        in_region().fold(f64::INFINITY, f64::min)
    };
    if sweep_hi <= sweep_lo {
        error!(
            "E-320: {}'s sweep range is degenerate ({sweep_lo} to {sweep_hi})",
            F::NAME
        );
        return None;
    }
    if cap.is_some() && cap_max >= sweep_lo {
        error!(
            "E-320: {}'s cap reaches {cap_max} inside the region while the sweep floor is \
             {sweep_lo}, so at the bottom rungs the cap alone makes samples air \
             (P-176's vacuity control)",
            F::NAME
        );
    }
    let sample_ms = started.elapsed().as_secs_f64() * 1000.0;

    let ladder = ladder(sweep_lo, sweep_hi);
    let slice_z = samples / 2;
    let started = Instant::now();
    let mut uf = Uf::with(base.len());
    let mut values = vec![0.0f64; base.len()];
    let mut profile = Vec::with_capacity(ladder.len());
    for iso in &ladder {
        fill_region(&base, &cap_values, *iso, &mut values);
        let (_, components, air, _, sizes) = census_air(&values, &shape)?;
        let largest = sizes.first().copied().unwrap_or(0);
        let largest_fraction = if air == 0 {
            0.0
        } else {
            f64::from(largest) / air as f64
        };
        let plane = census_2d(&values, [samples; 3], slice_z, &mut uf);
        profile.push(Rung {
            iso: *iso,
            components,
            air,
            largest_fraction,
            giant: largest_fraction > GIANT_SHARE,
            air_2d: plane.air,
            giant_2d: plane.largest_fraction > GIANT_SHARE,
        });
    }
    let pre_pass_ms = started.elapsed().as_secs_f64() * 1000.0;

    let onset_3d = persistent_onset(&profile, |r| r.giant);
    let onset_2d = persistent_onset(&profile, |r| r.giant_2d);
    let fragmented_rows = profile
        .iter()
        .filter(|r| !r.giant && r.components >= MANY_COMPONENTS)
        .count();
    let giant_rows = profile
        .iter()
        .filter(|r| r.giant && r.air as f64 / region as f64 >= REAL_AIR_SHARE)
        .count();
    let slice_rows = profile.iter().filter(|r| r.air_2d >= SLICE_FLOOR).count();
    if slice_rows == 0 {
        error!(
            "E-320: no rung of {}'s ladder puts {SLICE_FLOOR} air pixels on the z={slice_z} \
             plane, so the 2D column below is a default rather than a control \
             (P-176's vacuity control)",
            F::NAME
        );
    }

    // The HUD is the evidence and a headless capture has no HUD to read. One
    // line per pre-pass, so `ISOMESH_CAPTURE` leaves the sweep in the log where a
    // script can hold it against p-176.csv.
    info!(
        "E-320 {} at {samples}^3: region {region} samples ({}), sweep {sweep_hi:.6} down to \
         {sweep_lo:.6} in {} rungs; onset 3D {} / 2D {}; {fragmented_rows} fragmented rungs, \
         {giant_rows} giant rungs over {:.0}% air, {slice_rows} populated slices; \
         sample {sample_ms:.1} ms, pre-pass {pre_pass_ms:.1} ms",
        F::NAME,
        if cap.is_some() {
            "cap_sdf <= -0.05"
        } else {
            "full domain"
        },
        ladder.len(),
        onset_label(onset_3d.and_then(|i| ladder.get(i).copied())),
        onset_label(onset_2d.and_then(|i| ladder.get(i).copied())),
        REAL_AIR_SHARE * 100.0,
    );

    Some(Sampled {
        field: index,
        name: F::NAME,
        samples,
        shape,
        origin: min,
        cell_size,
        domain_min: Vec3::new(min[0] as f32, min[1] as f32, min[2] as f32),
        domain_max: Vec3::new(max[0] as f32, max[1] as f32, max[2] as f32),
        base,
        cap: cap_values,
        cap_field: cap,
        region,
        mask_rule: if cap.is_some() {
            "cap_sdf <= -0.05"
        } else {
            "full domain"
        },
        cap_max,
        sweep_lo,
        sweep_hi,
        ladder,
        profile,
        onset_3d,
        onset_2d,
        fragmented_rows,
        giant_rows,
        slice_rows,
        slice_z,
        pre_pass_ms,
    })
}

/// Fill `out` with [`region_value`] at every sample, for one isovalue.
fn fill_region(base: &[f64], cap: &[f64], iso: f64, out: &mut [f64]) {
    for (i, slot) in out.iter_mut().enumerate() {
        let b = base.get(i).copied().unwrap_or(-1.0);
        *slot = region_value(b, cap.get(i).copied(), iso);
    }
}

/// An onset isovalue, or the word the CSV uses when there is none.
fn onset_label(iso: Option<f64>) -> String {
    iso.map_or_else(|| String::from("none"), |v| format!("{v:.6}"))
}

// ─── app ────────────────────────────────────────────────────────────────────

fn main() {
    App::new()
        .add_plugins(DefaultPlugins.set(WindowPlugin {
            primary_window: Some(Window {
                title: "isomesh - E-320 cave percolation".into(),
                ..default()
            }),
            ..default()
        }))
        .add_plugins(CommonPlugin)
        .init_gizmo_group::<SliceGizmos>()
        .insert_resource(Pinned(
            common::samples_override().map(|n| n.clamp(MIN_SAMPLES, MAX_SAMPLES)),
        ))
        .insert_resource(Field(0))
        .insert_resource(Cursor(0))
        .init_resource::<ShowPanel>()
        .init_resource::<Grid>()
        .init_resource::<Census>()
        .add_systems(Startup, setup)
        // `PreUpdate` for E-306's reason: the harness's `update_hud` and its
        // `capture_sequence` both live in `Update` with no ordering against an
        // example's own systems, so in `Update` the HUD would render a frame-old
        // census beside a current isovalue (E-308, E-312).
        .add_systems(PreUpdate, (controls, rebuild, frame_camera, report).chain())
        .add_systems(Update, draw_slice)
        .run();
}

/// Captured frames one down-and-up pass of the ladder takes.
fn capture_period(rungs: usize) -> usize {
    (rungs * 2).saturating_sub(2).max(1)
}

fn setup(
    mut commands: Commands,
    mut materials: ResMut<Assets<StandardMaterial>>,
    mut gizmo_config: ResMut<GizmoConfigStore>,
    mut camera: Query<&mut OrbitCamera>,
) {
    for mut orbit in &mut camera {
        orbit.yaw = VIEW_YAW;
        orbit.pitch = VIEW_PITCH;
    }

    let (slices, _) = gizmo_config.config_mut::<SliceGizmos>();
    slices.line.width = 2.4;
    slices.depth_bias = -0.6;

    // White base colour, or the vertex colours the census wrote would be tinted
    // by it and the palette would stop being readable (E-301). Double-sided
    // because a cave is seen from the inside as often as from the outside.
    commands.insert_resource(SurfaceMaterial(materials.add(StandardMaterial {
        base_color: Color::WHITE,
        perceptual_roughness: 0.55,
        metallic: 0.02,
        double_sided: true,
        cull_mode: None,
        ..default()
    })));

    // Unlit: the panel is a chart, not an object. A lit flat quad facing away
    // from the single directional light draws the whole slice in shadow, and a
    // component's colour is the entire readout.
    commands.insert_resource(PanelMaterial(materials.add(StandardMaterial {
        base_color: Color::WHITE,
        unlit: true,
        double_sided: true,
        cull_mode: None,
        ..default()
    })));

    commands.spawn(DemoDomain {
        min: Vec3::splat(-1.0),
        max: Vec3::splat(1.0),
    });
}

/// Field, isovalue and the panel toggle.
///
/// Under capture the isovalue follows the frame counter, because an example whose
/// subject only changes on a keypress captures as a still frame. The pass is a
/// ping-pong so the clip ends where it started and loops.
fn controls(
    keys: Res<ButtonInput<KeyCode>>,
    capture: Res<Capture>,
    flags: Res<ViewFlags>,
    grid: Res<Grid>,
    mut field: ResMut<Field>,
    mut cursor: ResMut<Cursor>,
    mut panel: ResMut<ShowPanel>,
) {
    field.0 = flags.field.min(FIELD_COUNT - 1);

    let rungs = grid.0.as_ref().map_or(RUNGS, |s| s.ladder.len()).max(1);
    if capture.is_active() {
        let period = capture_period(rungs);
        let phase = capture.taken as usize % period;
        cursor.0 = if phase < rungs { phase } else { period - phase };
    } else {
        if keys.just_pressed(KeyCode::BracketRight) && cursor.0 + 1 < rungs {
            cursor.0 += 1;
        }
        if keys.just_pressed(KeyCode::BracketLeft) && cursor.0 > 0 {
            cursor.0 -= 1;
        }
    }
    cursor.0 = cursor.0.min(rungs - 1);

    if keys.just_pressed(KeyCode::KeyS) {
        panel.0 = !panel.0;
    }
}

/// Sample when the field or resolution changes, census when the rung does,
/// extract, and colour.
#[allow(clippy::too_many_arguments)]
fn rebuild(
    field: Res<Field>,
    cursor: Res<Cursor>,
    pinned: Res<Pinned>,
    flags: Res<ViewFlags>,
    mut grid: ResMut<Grid>,
    mut census: ResMut<Census>,
    mut meshes: ResMut<Assets<Mesh>>,
    mut surfaces: Query<&mut Mesh3d, (With<SurfaceMesh>, Without<PanelMesh>)>,
    mut panels: Query<&mut Mesh3d, (With<PanelMesh>, Without<SurfaceMesh>)>,
    mut domain: Query<&mut DemoDomain>,
    mut commands: Commands,
    surface_material: Res<SurfaceMaterial>,
    panel_material: Res<PanelMaterial>,
    mut uf: Local<Option<Uf>>,
    mut last: Local<Option<(usize, u32, usize)>>,
) {
    let samples = pinned.0.unwrap_or(DEFAULT_SAMPLES);
    let stale = grid
        .0
        .as_ref()
        .is_none_or(|s| s.field != field.0 || s.samples != samples);
    if stale {
        grid.0 = sample(field.0, samples);
        *uf = None;
    }
    let Some(sampled) = grid.0.as_ref() else {
        return;
    };

    // The HUD puts a committed row beside a live one, so the two had better be
    // about the same field. Checked here, once per resample, rather than in
    // `report`, which would say it sixty times a second.
    if stale
        && let Some(cited) = CITED.get(sampled.field)
        && cited.name != sampled.name
    {
        error!(
            "E-320: the cited p-176.csv row is {} but the field on screen is {} -- CITED and \
             the field list have drifted apart, so every quoted number below is the wrong row",
            cited.name, sampled.name
        );
    }

    let rung = cursor.0.min(sampled.ladder.len().saturating_sub(1));
    let key = (sampled.field, sampled.samples, rung);
    if *last == Some(key) && !flags.remesh_requested && !stale {
        return;
    }
    *last = Some(key);

    for mut d in &mut domain {
        d.min = sampled.domain_min;
        d.max = sampled.domain_max;
    }

    let Some(iso) = sampled.ladder.get(rung).copied() else {
        return;
    };
    let solver = uf.get_or_insert_with(|| Uf::with(sampled.base.len()));

    // ── both instruments, on the same array ─────────────────────────────────
    let started = Instant::now();
    let mut values = vec![0.0f64; sampled.base.len()];
    fill_region(&sampled.base, &sampled.cap, iso, &mut values);
    let Some((graph, components, air_samples, labels, sizes)) = census_air(&values, &sampled.shape)
    else {
        return;
    };
    let (uf_air, uf_sizes) = census_uf_3d(&values, [sampled.samples; 3], solver);
    let plane = census_2d(&values, [sampled.samples; 3], sampled.slice_z, solver);
    let census_ms = started.elapsed().as_secs_f64() * 1000.0;

    let largest = sizes.first().copied().unwrap_or(0);
    let largest_fraction = if air_samples == 0 {
        0.0
    } else {
        f64::from(largest) / air_samples as f64
    };

    // ── the surface, and its colour ─────────────────────────────────────────
    let started = Instant::now();
    let buffer = extract(sampled, iso);
    let extract_ms = started.elapsed().as_secs_f64() * 1000.0;
    let ranks = rank_labels(&graph, labels);
    let surface = buffer
        .as_ref()
        .and_then(|b| coloured_surface(b, sampled, &graph, &ranks));
    let panel = slice_panel(sampled, &plane);

    let vertices = buffer.as_ref().map_or(0, |b| b.positions.len());
    let triangles = buffer.as_ref().map_or(0, |b| b.indices.len() / 3);

    *census = Census {
        samples: sampled.samples,
        field: sampled.field,
        name: sampled.name,
        rung,
        rungs: sampled.ladder.len(),
        iso,
        components,
        air: air_samples,
        largest,
        largest_fraction,
        top_sizes: sizes.iter().take(TOP_SIZES).copied().collect(),
        uf_components: uf_sizes.len() as u64,
        uf_air,
        uf_sizes_match: uf_sizes == sizes,
        labels_tight: labels == components,
        components_2d: plane.components,
        air_2d: plane.air,
        largest_fraction_2d: plane.largest_fraction,
        census_ms,
        extract_ms,
        vertices,
        triangles,
    };

    if !census.agrees() {
        // A classification that means "the two instruments disagree about the
        // topology" gets an `error!` and not a colour alone (E-301): every colour
        // on screen is one of `Air`'s labels, so a disagreement makes the picture
        // wrong rather than merely surprising.
        error!(
            "E-320: Air and the union-find disagree at iso {iso:.6} on {} {}^3: \
             {} vs {} components, {} vs {} air samples, sizes match {}, labels {} for {} components",
            sampled.name,
            sampled.samples,
            components,
            uf_sizes.len(),
            air_samples,
            uf_air,
            uf_sizes == sizes,
            labels,
            components,
        );
    }

    info!(
        "E-320 {} {}^3 rung {}/{} iso {iso:+.6}: 3D {components} components, largest \
         {largest_fraction:.4} of {air_samples} air ({}); 2D z={} {} components, largest \
         {:.4} of {} air ({}); Air vs union-find {}; census {census_ms:.1} ms, \
         extract {extract_ms:.1} ms, {vertices} vertices",
        sampled.name,
        sampled.samples,
        rung + 1,
        sampled.ladder.len(),
        verdict(largest_fraction > GIANT_SHARE),
        sampled.slice_z,
        plane.components,
        plane.largest_fraction,
        plane.air,
        verdict(plane.largest_fraction > GIANT_SHARE),
        if census.agrees() { "AGREE" } else { "DISAGREE" },
    );

    // ── publish ─────────────────────────────────────────────────────────────
    let surface = mesh3d(&mut meshes, surface);
    if surfaces.is_empty() {
        commands.spawn((
            surface,
            MeshMaterial3d(surface_material.0.clone()),
            SurfaceMesh,
            DemoMesh,
        ));
    } else {
        for mut mesh in &mut surfaces {
            *mesh = surface.clone();
        }
    }

    let panel = mesh3d(&mut meshes, panel);
    if panels.is_empty() {
        commands.spawn((panel, MeshMaterial3d(panel_material.0.clone()), PanelMesh));
    } else {
        for mut mesh in &mut panels {
            *mesh = panel.clone();
        }
    }
}

/// A [`Mesh3d`] for a result that may be empty.
///
/// **An empty [`Mesh`] must not become an asset.** `bevy_render`'s
/// `MeshAllocator::allocate_meshes` skips any mesh whose vertex buffer is zero
/// bytes and then copies into it unconditionally, so one empty mesh logs
/// `Use-after-free: attempted to copy element data for an unallocated key` twice
/// a frame, forever (E-305, E-307). The top rung of every ladder here reaches it:
/// one air sample makes no surface at all. `Mesh3d::default()` names no asset and
/// draws nothing, which is what an empty result actually wants.
fn mesh3d(meshes: &mut Assets<Mesh>, mesh: Option<Mesh>) -> Mesh3d {
    mesh.map_or_else(Mesh3d::default, |m| Mesh3d(meshes.add(m)))
}

/// Marching Cubes on the region function at one isovalue.
fn extract(sampled: &Sampled, iso: f64) -> Option<MeshBuffer<f64>> {
    let mut buffer = MeshBuffer::<f64>::new();
    let result = match sampled.field {
        0 => MarchingCubes::<f64>::new().extract(
            &Excursion {
                field: &FbmTerrain::<f64>::canonical(),
                iso,
                cap: sampled.cap_field,
            },
            &sampled.shape,
            sampled.origin,
            sampled.cell_size,
            &mut buffer,
        ),
        _ => MarchingCubes::<f64>::new().extract(
            &Excursion {
                field: &noise_cavity::<f64>(),
                iso,
                cap: sampled.cap_field,
            },
            &sampled.shape,
            sampled.origin,
            sampled.cell_size,
            &mut buffer,
        ),
    };
    if let Err(error) = result {
        error!(
            "E-320: marching cubes failed on {} at iso {iso:.6}: {error}",
            sampled.name
        );
        return None;
    }
    if buffer.positions.is_empty() {
        return None;
    }
    Some(buffer)
}

/// Component rank per `Air` label — 0 for the largest, `u32::MAX` for a retired
/// or unissued one.
///
/// By size descending, ties by label, so the palette is stable in size order
/// rather than in whatever order the flood happened to reach the components in.
fn rank_labels(air: &Air, labels: u64) -> Vec<u32> {
    let mut live: Vec<(u32, u32)> = (0..labels as u32)
        .map(|l| (air.component_size(l), l))
        .filter(|(size, _)| *size > 0)
        .collect();
    live.sort_unstable_by(|a, b| b.0.cmp(&a.0).then(a.1.cmp(&b.1)));
    let mut rank = vec![u32::MAX; labels as usize];
    for (r, (_, label)) in live.iter().enumerate() {
        if let Some(slot) = rank.get_mut(*label as usize) {
            *slot = r as u32;
        }
    }
    rank
}

/// The colour a component rank is drawn in.
fn rank_colour(rank: u32) -> [f32; 4] {
    if rank == u32::MAX {
        return ORPHAN_COLOUR;
    }
    if rank == 0 {
        return GIANT_COLOUR;
    }
    POCKET_COLOURS[(rank as usize - 1) % POCKET_COLOURS.len()]
}

/// sRGB as a human picks it into the linear RGBA [`Mesh::ATTRIBUTE_COLOR`]
/// wants. Feeding sRGB in raw renders it washed out (E-208).
fn linear(srgb: [f32; 4]) -> [f32; 4] {
    Color::srgba(srgb[0], srgb[1], srgb[2], srgb[3])
        .to_linear()
        .to_f32_array()
}

/// The `f64` extraction as a Bevy mesh, one colour per vertex.
///
/// Cast rather than re-extracted in `f32`: the numbers on the HUD are `f64`
/// numbers, so the mesh the picture is drawn from has to be the one they were
/// computed alongside (E-307). `f32` findings quoted from an `f64` measurement is
/// a break this repository has already had (E-305, `game_edit_tape_trim.rs:126`).
fn coloured_surface(
    buffer: &MeshBuffer<f64>,
    sampled: &Sampled,
    air: &Air,
    ranks: &[u32],
) -> Option<Mesh> {
    let mut builder = MeshBuilder::new();
    for i in 0..buffer.positions.len() {
        let (Some(p), Some(n)) = (buffer.positions.get(i), buffer.normals.get(i)) else {
            continue;
        };
        builder.vertex(
            [p[0] as f32, p[1] as f32, p[2] as f32],
            [n[0] as f32, n[1] as f32, n[2] as f32],
        );
        let rank = label_at(air, sampled, *p)
            .and_then(|l| ranks.get(l as usize).copied())
            .unwrap_or(u32::MAX);
        builder.colors_mut().push(linear(rank_colour(rank)));
    }
    for t in buffer.indices.as_chunks::<3>().0 {
        let (Some(a), Some(b), Some(c)) = (t.first(), t.get(1), t.get(2)) else {
            continue;
        };
        builder.triangle(*a, *b, *c);
    }
    if builder.vertex_count() == 0 {
        return None;
    }
    Some(builder.into_mesh())
}

/// The `Air` label of the air sample nearest `p`, over the cell holding it.
///
/// A Marching Cubes vertex sits on a grid edge between an inside sample and an
/// air one, so the nearest air corner of its cell **is** the air end of its own
/// edge — exact rather than approximate for every edge vertex, and the closest
/// honest answer for the cell-local centroid vertices A-015 adds. `None` when no
/// corner of the cell is air, which is what [`ORPHAN_COLOUR`] draws.
fn label_at(air: &Air, sampled: &Sampled, p: [f64; 3]) -> Option<u32> {
    let last = sampled.samples.saturating_sub(2);
    let mut base = [0u32; 3];
    for axis in 0..3 {
        let local = (p[axis] - sampled.origin[axis]) / sampled.cell_size;
        base[axis] = if local <= 0.0 {
            0
        } else {
            (local as u32).min(last)
        };
    }
    let mut best: Option<(f64, u32)> = None;
    for corner in 0..8u32 {
        let c = [
            base[0] + (corner & 1),
            base[1] + ((corner >> 1) & 1),
            base[2] + ((corner >> 2) & 1),
        ];
        let Some(label) = air.label_of(c) else {
            continue;
        };
        let mut d = 0.0f64;
        for axis in 0..3 {
            let world = sampled.origin[axis] + f64::from(c[axis]) * sampled.cell_size;
            d += (world - p[axis]) * (world - p[axis]);
        }
        if best.is_none_or(|(bd, _)| d < bd) {
            best = Some((d, label));
        }
    }
    best.map(|(_, label)| label)
}

/// The box the slice panel occupies.
///
/// The plane's own `xy` extent, offset along `+x` by [`PANEL_GAP`] domain
/// widths, at the plane's own `z`. **One definition**, because the panel is
/// drawn from it and the camera is framed from it, and a camera framed from a
/// second copy of this arithmetic is how a demo ends up photographing the edge
/// of its own evidence.
fn panel_bounds(sampled: &Sampled) -> (Vec3, Vec3) {
    let h = sampled.cell_size as f32;
    // `samples` quads of side `h`, not `samples - 1`: the last quad *starts* at
    // the last sample and reaches one cell past it, so the panel is one cell
    // wider than the domain box beside it.
    let extent = sampled.samples as f32 * h;
    let width = sampled.domain_max.x - sampled.domain_min.x;
    let lo = Vec3::new(
        sampled.domain_min.x + width * PANEL_GAP,
        sampled.domain_min.y,
        sampled.origin[2] as f32 + sampled.slice_z as f32 * h,
    );
    (lo, lo + Vec3::new(extent, extent, 0.0))
}

/// The censused plane as one quad per air pixel, beside the volume.
///
/// In the plane's own `xy` orientation and at the plane's own `z`, offset along
/// `+x` so it sits next to the solid rather than inside it. Same colouring rule
/// as the surface, so "gold is the largest component" reads across both.
fn slice_panel(sampled: &Sampled, plane: &Plane) -> Option<Mesh> {
    let n = sampled.samples as usize;
    let h = sampled.cell_size as f32;
    let (lo, _) = panel_bounds(sampled);

    let mut builder = MeshBuilder::new();
    for j in 0..n {
        for i in 0..n {
            let rank = plane.ranks.get(i + j * n).copied().unwrap_or(u32::MAX);
            if rank == u32::MAX {
                continue;
            }
            let colour = linear(rank_colour(rank));
            let (x, y) = (lo.x + i as f32 * h, lo.y + j as f32 * h);
            let base = builder.vertex_count() as u32;
            for (dx, dy) in [(0.0, 0.0), (h, 0.0), (h, h), (0.0, h)] {
                builder.vertex([x + dx, y + dy, lo.z], [0.0, 0.0, 1.0]);
                builder.colors_mut().push(colour);
            }
            builder.triangle(base, base + 1, base + 2);
            builder.triangle(base, base + 2, base + 3);
        }
    }
    if builder.vertex_count() == 0 {
        return None;
    }
    Some(builder.into_mesh())
}

/// Frame the volume and the panel together, from the field's own domain.
///
/// Never a hardcoded radius: a fixed 8.0 put the camera *inside* the gyroid once
/// and the committed screenshot was a picture of an inner wall
/// (`critical_cells.rs:343`). The two fields' domains differ by a factor of four
/// here, so a constant would frame one of them and lose the other.
fn frame_camera(grid: Res<Grid>, mut camera: Query<&mut OrbitCamera>) {
    let Some(sampled) = grid.0.as_ref() else {
        return;
    };
    let (panel_lo, panel_hi) = panel_bounds(sampled);
    let lo = sampled.domain_min.min(panel_lo);
    let hi = sampled.domain_max.max(panel_hi);
    let centre = (lo + hi) * 0.5;
    let radius = (hi.x - lo.x) * VIEW_RADIUS_WIDTHS;
    for mut orbit in &mut camera {
        // `orbit_camera` puts the eye at `focus + dir * radius`, so the view
        // direction is `-dir` and a focus moved along `-right` puts the subject
        // right of centre. Applied in the camera's own basis, so it is one
        // screen-space nudge however far `ISOMESH_SPIN` has turned.
        let dir = Vec3::new(
            orbit.yaw.cos() * orbit.pitch.cos(),
            orbit.pitch.sin(),
            orbit.yaw.sin() * orbit.pitch.cos(),
        );
        let forward = -dir;
        let right = forward.cross(Vec3::Y).normalize_or_zero();
        let up = right.cross(forward).normalize_or_zero();
        orbit.focus =
            centre - right * (SUBJECT_OFFSET.x * radius) + up * (SUBJECT_OFFSET.y * radius);
        orbit.radius = radius;
    }
}

/// Outline the censused plane inside the volume, and show or hide the panel.
///
/// Both here rather than in `rebuild` because both are per-frame view state, and
/// the visibility write is conditional: an unconditional `*visible = ...` marks
/// the component changed every frame and Bevy's visibility propagation is
/// change-driven, so it would turn a toggle nobody pressed into per-frame work.
fn draw_slice(
    grid: Res<Grid>,
    panel: Res<ShowPanel>,
    mut visibility: Query<&mut Visibility, With<PanelMesh>>,
    mut gizmos: Gizmos<SliceGizmos>,
) {
    const CYAN: Color = Color::srgb(0.15, 0.85, 1.0);

    let wanted = if panel.0 {
        Visibility::Inherited
    } else {
        Visibility::Hidden
    };
    for mut visible in &mut visibility {
        if *visible != wanted {
            *visible = wanted;
        }
    }

    let Some(sampled) = grid.0.as_ref() else {
        return;
    };
    if !panel.0 {
        return;
    }
    let z = sampled.origin[2] as f32 + sampled.slice_z as f32 * sampled.cell_size as f32;
    let (min, max) = (sampled.domain_min, sampled.domain_max);
    let corners = [
        Vec3::new(min.x, min.y, z),
        Vec3::new(max.x, min.y, z),
        Vec3::new(max.x, max.y, z),
        Vec3::new(min.x, max.y, z),
    ];
    for i in 0..4 {
        gizmos.line(corners[i], corners[(i + 1) % 4], CYAN);
    }
}

// ─── the HUD ────────────────────────────────────────────────────────────────

/// `GIANT` or `fragmented`, so the regime is a word a reader can check at a
/// glance rather than a fraction they have to compare against a threshold.
fn verdict(giant: bool) -> &'static str {
    if giant { "GIANT" } else { "fragmented" }
}

/// The `sweep` row: the largest component's share of the air, one character per
/// rung, `' '` where there is no air at all.
fn ramp_row(profile: &[Rung]) -> String {
    profile
        .iter()
        .map(|r| {
            if r.air == 0 {
                ' '
            } else {
                let i = (r.largest_fraction * (RAMP.len() - 1) as f64).round();
                RAMP[(i.clamp(0.0, (RAMP.len() - 1) as f64)) as usize]
            }
        })
        .collect()
}

/// A `giant` row: the mark where the test passes, `-` where it does not.
fn giant_row(profile: &[Rung], mark: char, giant: impl Fn(&Rung) -> bool) -> String {
    profile
        .iter()
        .map(|r| if giant(r) { mark } else { '-' })
        .collect()
}

/// The `^` under the rung on screen.
fn cursor_row(rungs: usize, rung: usize) -> String {
    (0..rungs)
        .map(|i| if i == rung { '^' } else { ' ' })
        .collect()
}

/// The HUD. The numbers are the demo.
fn report(census: Res<Census>, grid: Res<Grid>, mut stats: ResMut<DemoStats>) {
    if census.samples == 0 {
        return;
    }
    let Some(sampled) = grid.0.as_ref() else {
        return;
    };
    let Some(cited) = CITED.get(census.field) else {
        return;
    };

    stats.title = format!(
        "E-320  cave percolation - {}   {}^3   [1,2] field  [ ] isovalue  S slice",
        census.name, census.samples,
    );

    // The fixed part of the panel is written here too, so one system owns every
    // number on screen and `rebuild` owns none of them.
    stats.vertices = census.vertices;
    stats.triangles = census.triangles;
    stats.extract_ms = census.extract_ms;

    let total = (census.samples as usize).pow(3);
    let air_share = if sampled.region == 0 {
        0.0
    } else {
        census.air as f64 / sampled.region as f64
    };
    let sizes = census
        .top_sizes
        .iter()
        .map(u32::to_string)
        .collect::<Vec<_>>()
        .join("|");
    let live_3d = onset_of(sampled, sampled.onset_3d);
    let live_2d = onset_of(sampled, sampled.onset_2d);
    let at_zero = sampled
        .profile
        .iter()
        .find(|r| r.iso.abs() < ZERO_RUNG)
        .map_or_else(
            || String::from("not swept"),
            |r| format!("{:.6}", r.largest_fraction),
        );

    stats.extra = vec![
        format!(
            "{:>10.4} isovalue    rung {} of {}    region {} of {} samples ({})",
            census.iso,
            census.rung + 1,
            census.rungs,
            sampled.region,
            total,
            sampled.mask_rule,
        ),
        format!(
            "{:>10} air components    largest holds {} samples, {:.4} of the air    {}",
            census.components,
            census.largest,
            census.largest_fraction,
            verdict(census.largest_fraction > GIANT_SHARE),
        ),
        format!(
            "{:>10} air samples = {air_share:.4} of the region    biggest: {sizes}",
            census.air,
        ),
        String::new(),
        format!(
            "  3D  6-connectivity          {:>6} components   largest {:.4}   {}",
            census.components,
            census.largest_fraction,
            verdict(census.largest_fraction > GIANT_SHARE),
        ),
        format!(
            "  2D  4-connectivity  z={:<4}  {:>6} components   largest {:.4}   {}",
            sampled.slice_z,
            census.components_2d,
            census.largest_fraction_2d,
            verdict(census.largest_fraction_2d > GIANT_SHARE),
        ),
        format!(
            "      6 in 3D and 4 in 2D is Air::neighbours' own rule, not this demo's; \
             {} air pixels in the slice",
            census.air_2d,
        ),
        String::new(),
        format!("  sweep {}", ramp_row(&sampled.profile)),
        format!(
            "  giant {}   high isovalue left, low right",
            giant_row(&sampled.profile, 'G', |r| r.giant)
        ),
        format!(
            "  2D    {}   G / g: one component holds over half the air",
            giant_row(&sampled.profile, 'g', |r| r.giant_2d)
        ),
        format!(
            "        {}   sweep ramps '.' to '@'",
            cursor_row(sampled.profile.len(), census.rung)
        ),
        format!(
            "  onset live 3D {live_3d}   live 2D {live_2d}   (persistent: giant on every rung below)"
        ),
        format!(
            "        {} fragmented rungs, {} giant rungs over {:.0}% air, {} slices over {} \
             pixels -- the sweep visits both regimes and the control is populated",
            sampled.fragmented_rows,
            sampled.giant_rows,
            REAL_AIR_SHARE * 100.0,
            sampled.slice_rows,
            SLICE_FLOOR,
        ),
        String::new(),
        format!(
            "  Air vs an independent union-find: {} == {} components, {} == {} air, \
             sizes {}, labels {} -> {}",
            census.components,
            census.uf_components,
            census.air,
            census.uf_air,
            if census.uf_sizes_match {
                "same"
            } else {
                "DIFFER"
            },
            if census.labels_tight {
                "tight"
            } else {
                "STALE"
            },
            if census.agrees() { "AGREE" } else { "DISAGREE" },
        ),
        String::new(),
        format!(
            "  P-176 at {CITED_RESOLUTION}^3 (docs/experiments/p-176.csv): C1 HELD, \
             C3 HELD on {CITED_ROWS} of {CITED_ROWS} rows,"
        ),
        format!(
            "        C2 split 42/42 -- HELD on noise_cavity, FALSIFIED on fbm_terrain. \
             this field: C2 {}",
            if cited.c2 { "HELD" } else { "FALSIFIED" },
        ),
        format!(
            "        transition 3D {:.6}   2D {}   largest at iso 0: {:.6} (live {at_zero})",
            cited.onset_3d,
            onset_label(cited.onset_2d),
            cited.largest_at_zero,
        ),
        format!(
            "        iso {:.6} -> {} components, largest {:.4}, then iso {:.6} -> {} \
             components, largest {:.4}",
            cited.fragmented.0,
            cited.fragmented.1,
            cited.fragmented.2,
            cited.giant.0,
            cited.giant.1,
            cited.giant.2,
        ),
        format!(
            "        {} of 42 rungs were one component there; unmasked, noise_cavity at iso 0 \
             is {:.4} air in {} components ({:.4}) -- the cap is why",
            cited.single_rows, CITED_SHELL.0, CITED_SHELL.1, CITED_SHELL.2,
        ),
        String::from(
            "        fbm_terrain is hash-based lattice noise, not a Gaussian field, and that",
        ),
        String::from(
            "        was registered in advance as a real risk. Its air is {y >= h(x,z) + iso},",
        ),
        String::from(
            "        the region above a graph: one component in 3D and in 2D alike, so there is",
        ),
        String::from(
            "        no dimension gap there to find. C2's 42/42 is that, honestly reported.",
        ),
        format!(
            "        range {:.6} to {:.6}   cap max in region {:.6} vs floor {:.6}   \
             pre-pass {:.0} ms, census {:.1} ms",
            sampled.sweep_hi,
            sampled.sweep_lo,
            sampled.cap_max,
            sampled.sweep_lo,
            sampled.pre_pass_ms,
            census.census_ms,
        ),
    ];
}

/// A live onset index rendered as its isovalue, or `none`.
fn onset_of(sampled: &Sampled, index: Option<usize>) -> String {
    onset_label(index.and_then(|i| sampled.ladder.get(i).copied()))
}

// ─── tests ──────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use bevy::ecs::system::RunSystemOnce;
    use bevy::time::TimeUpdateStrategy;
    use std::time::Duration;

    /// A fixed frame, so nothing here depends on how fast the machine is.
    const FRAME: Duration = Duration::from_millis(16);

    /// Samples per axis in the harness: **the resolution the demo opens at**.
    ///
    /// Not a smaller one. The verdicts below are about a transition, and a
    /// transition is a property of how well the grid resolves the field's own
    /// features — measured, `noise_cavity` at `17^3` never percolates at any
    /// isovalue, because `h = 0.25` against a `0.29` feature leaves the blobs
    /// disconnected. A test at a resolution the demo never uses would assert a
    /// verdict nobody sees. Both fields, whole ladder, twice over: 0.4 s.
    ///
    /// Inserted as [`Pinned`] rather than through `ISOMESH_SAMPLES`, because
    /// `std::env::set_var` is `unsafe` and this crate's `[lints.rust]` says
    /// `unsafe_code = "forbid"`.
    const TEST_SAMPLES: u32 = DEFAULT_SAMPLES;

    /// The demo's own systems, in an `App` with no window and no renderer.
    ///
    /// This is the closest thing to running the demo that a machine with no
    /// display can do. `frame_camera` and `draw_slice` are left out on purpose:
    /// they want a camera and gizmos, and neither writes a number. `report` is
    /// left out too and run as a one-shot below, which is the same system the
    /// demo runs every frame.
    fn harness(field: usize) -> App {
        let mut app = App::new();
        app.add_plugins(MinimalPlugins)
            .add_plugins(bevy::asset::AssetPlugin::default())
            .init_asset::<Mesh>()
            .insert_resource(TimeUpdateStrategy::ManualDuration(FRAME))
            .init_resource::<ButtonInput<KeyCode>>()
            .init_resource::<ViewFlags>()
            .insert_resource(Capture::default())
            .insert_resource(Pinned(Some(TEST_SAMPLES)))
            .insert_resource(Field(field))
            .insert_resource(Cursor(0))
            .insert_resource(ShowPanel(true))
            .insert_resource(SurfaceMaterial(Handle::default()))
            .insert_resource(PanelMaterial(Handle::default()))
            .init_resource::<Grid>()
            .init_resource::<Census>()
            .init_resource::<DemoStats>()
            .add_systems(Update, (controls, rebuild).chain());
        // `ViewFlags::default` reads `ISOMESH_FIELD`, and `controls` mirrors it
        // over `Field` every frame — so the field has to be set where `controls`
        // reads it, or a developer with that variable exported would test a
        // different field than the one this asked for.
        app.world_mut().resource_mut::<ViewFlags>().field = field;
        app
    }

    /// One frame, with the input clearing `InputPlugin` would have done.
    ///
    /// **`reset_all`, not `clear`.** `clear` drops `just_pressed` and leaves the
    /// key *held*, and `ButtonInput::press` only registers a `just_pressed` for a
    /// key that was not already down — so a second `]` on a still-held key is a
    /// silent no-op and a loop that walks the ladder one rung per press walks
    /// exactly one rung. Measured: the sweep below stopped at rung 1 of 42.
    fn step(app: &mut App) {
        app.update();
        app.world_mut()
            .resource_mut::<ButtonInput<KeyCode>>()
            .reset_all();
    }

    /// Walk down the ladder until the census reports a giant component **over a
    /// region that is really air**, or run out of rungs.
    ///
    /// A stall detector rather than a flat frame count: the rung a field's
    /// transition lands on depends on the field and on the resolution, so a
    /// hardcoded number of `]` presses would be right for one of them and
    /// silently wrong for the other.
    ///
    /// The [`REAL_AIR_SHARE`] floor is P-176's own control and it is what makes
    /// this a stop worth asserting at: the top rung of every ladder admits one
    /// voxel, that voxel is all of its own component, and "largest holds 100% of
    /// the air" is then a fact about the ladder rather than about the field.
    /// Measured without the floor: this stopped at rung 2 of 42 with 125 air
    /// samples out of 274,625.
    ///
    /// Returns the rung it stopped on.
    fn drain_to_giant(app: &mut App) -> usize {
        step(app);
        let rungs = app.world().resource::<Census>().rungs;
        let region = app
            .world()
            .resource::<Grid>()
            .0
            .as_ref()
            .map_or(1, |s| s.region)
            .max(1);
        for _ in 0..rungs {
            let census = app.world().resource::<Census>();
            if census.largest_fraction > GIANT_SHARE
                && census.air as f64 / region as f64 >= REAL_AIR_SHARE
            {
                break;
            }
            app.world_mut()
                .resource_mut::<ButtonInput<KeyCode>>()
                .press(KeyCode::BracketRight);
            step(app);
        }
        app.world().resource::<Census>().rung
    }

    /// The HUD reports the live census, the dimension comparison and the row it
    /// cites — including the clause that was falsified.
    ///
    /// **This is the only way to see this demo's screen on a machine with no
    /// display**, and every line it checks is a line a reader is meant to read:
    /// the two connectivity rules, the 2D control beside the 3D census, the live
    /// agreement between `Air` and a second instrument, and P-176's own verdicts
    /// with the field-specific half of C2 named. A test that only checked
    /// [`Census`]'s fields would pass with any of them missing from the panel.
    #[test]
    fn the_hud_reports_both_dimensions_and_the_split_verdict() {
        for field in 0..FIELD_COUNT {
            let mut app = harness(field);
            let rung = drain_to_giant(&mut app);
            app.world_mut()
                .run_system_once(report)
                .expect("the HUD system");

            let census = app.world().resource::<Census>();
            let name = census.name;
            // Non-vacuity: the population, not its shape. A census over an empty
            // air set agrees with anything and percolates trivially.
            assert!(
                census.air > 1 && census.components > 0,
                "{name}: the census found {} air samples in {} components at rung {rung}, \
                 so every assertion below is over an empty set",
                census.air,
                census.components
            );
            assert!(
                census.largest_fraction > GIANT_SHARE,
                "{name}: no rung of the {}-rung ladder produced a giant component, which is \
                 C1's own claim (stopped at rung {rung}, largest {:.4})",
                census.rungs,
                census.largest_fraction
            );
            // C3, live, on this rung: the colours on screen are `Air`'s labels.
            assert!(
                census.agrees(),
                "{name}: Air and the union-find disagree at rung {rung} -- {} vs {} \
                 components, {} vs {} air samples",
                census.components,
                census.uf_components,
                census.air,
                census.uf_air
            );

            // C2's own shape, recomputed live rather than only quoted: a
            // persistent giant phase in 3D **and none anywhere in the slice**.
            // This is the assertion that keeps the demo honest — it fails if the
            // negative arm ever starts looking like the positive one, which is
            // the one way a demo of a split result can quietly become a lie.
            let cited = CITED.get(field).expect("a cited row per field");
            let (onset_3d, onset_2d) = {
                let grid = app.world().resource::<Grid>();
                let sampled = grid.0.as_ref().expect("the sampled field");
                (sampled.onset_3d, sampled.onset_2d)
            };
            assert!(
                onset_3d.is_some(),
                "{name}: the sweep found no persistent 3D onset, which is C1's claim"
            );
            assert_eq!(
                onset_3d.is_some() && onset_2d.is_none(),
                cited.c2,
                "{name}: the live dimension gap (3D onset {onset_3d:?}, 2D onset \
                 {onset_2d:?}) disagrees with p-176.csv's c2_holds = {}",
                cited.c2
            );

            let lines = app.world().resource::<DemoStats>().extra.clone();
            for line in &lines {
                println!("{line}");
            }

            let three = lines
                .iter()
                .find(|l| l.contains("6-connectivity"))
                .expect("the 3D connectivity line");
            assert!(
                three.contains("components") && three.contains("largest"),
                "{name}: the 3D line stopped carrying the census: {three}"
            );
            let two = lines
                .iter()
                .find(|l| l.contains("4-connectivity"))
                .expect("the 2D connectivity line");
            assert!(
                two.contains(&format!("z={}", TEST_SAMPLES / 2)),
                "{name}: the 2D line does not say which plane it censused: {two}"
            );
            assert!(
                lines
                    .iter()
                    .any(|l| l.contains("6 in 3D and 4 in 2D is Air::neighbours' own rule")),
                "{name}: the HUD stopped stating the connectivity rule"
            );

            let agreement = lines
                .iter()
                .find(|l| l.contains("Air vs an independent union-find"))
                .expect("the agreement line");
            assert!(
                agreement.contains("AGREE") && !agreement.contains("DISAGREE"),
                "{name}: the agreement line does not report agreement: {agreement}"
            );

            let clauses = lines
                .iter()
                .find(|l| l.contains("C1 HELD"))
                .expect("the P-176 verdict line");
            assert!(
                clauses.contains(&format!("C3 HELD on {CITED_ROWS} of {CITED_ROWS} rows")),
                "{name}: the verdict line stopped citing C3's row count: {clauses}"
            );
            let split = lines
                .iter()
                .find(|l| l.contains("C2 split 42/42"))
                .expect("the C2 split line");
            assert!(
                split.contains("HELD on noise_cavity")
                    && split.contains("FALSIFIED on fbm_terrain"),
                "{name}: the C2 line stopped naming which field held: {split}"
            );
            assert!(
                split.contains(&format!(
                    "this field: C2 {}",
                    if field == 1 { "HELD" } else { "FALSIFIED" }
                )),
                "{name}: the C2 line reports the wrong verdict for this field: {split}"
            );

            let transition = lines
                .iter()
                .find(|l| l.contains("transition 3D"))
                .expect("the cited transition line");
            assert!(
                transition.contains(&format!("{:.6}", cited.onset_3d))
                    && transition.contains(&onset_label(cited.onset_2d)),
                "{name}: the transition line stopped quoting p-176.csv's onsets: {transition}"
            );

            assert!(
                lines.iter().any(|l| l.contains("not a Gaussian field")),
                "{name}: the HUD dropped the registered caveat about fbm_terrain"
            );

            let profile = lines
                .iter()
                .find(|l| l.starts_with("  giant "))
                .expect("the giant row");
            assert!(
                profile.contains('G'),
                "{name}: the sweep found no giant rung to mark: {profile}"
            );
        }
    }

    /// The ladder spans the sampled range, descending, and contains `iso = 0`.
    ///
    /// The vacuity control P-176 registered, as a test rather than as a comment:
    /// a sweep that missed one of the two regimes would report an onset that is
    /// an artefact of where the ladder stops, and the registration's own sentence
    /// is about `iso = 0`.
    #[test]
    fn the_ladder_spans_the_range_and_lands_on_zero() {
        let rungs = ladder(-9.6, 9.42);
        assert_eq!(
            rungs.len(),
            RUNGS + 1,
            "a range straddling zero should gain the inserted zero rung: {rungs:?}"
        );
        assert!(
            rungs.windows(2).all(|w| w[0] > w[1]),
            "the ladder is not strictly descending: {rungs:?}"
        );
        assert!(
            rungs.iter().any(|v| v.abs() < ZERO_RUNG),
            "the ladder does not contain iso = 0: {rungs:?}"
        );
        let (first, last) = (
            rungs.first().copied().unwrap_or_default(),
            rungs.last().copied().unwrap_or_default(),
        );
        assert!(
            (first - 9.42).abs() < 1e-12 && (last + 9.6).abs() < 1e-12,
            "the ladder does not span the range it was given: {first} to {last}"
        );

        let no_zero = ladder(0.5, 1.5);
        assert_eq!(
            no_zero.len(),
            RUNGS,
            "a range that does not straddle zero should gain nothing: {no_zero:?}"
        );
    }

    /// Both objects are inside the frustum, on both fields.
    ///
    /// **This exists because nothing else here can see the screen.** Every other
    /// check in this file is about a string; a camera is about geometry, and the
    /// two failures this repository has actually had with one are a hardcoded
    /// radius that put the eye *inside* the subject (`critical_cells.rs:343`,
    /// `BACKLOG_ARCHIVE.md:157`) and a committed GIF that advertised a sweep it
    /// never performed (M-241, `BACKLOG_ARCHIVE.md:836`). Neither is visible in a
    /// HUD line. So the arithmetic in [`frame_camera`] is checked against the
    /// frustum it is derived from, on the two domains that differ by a factor of
    /// four — the twelve corners of the volume box and the four of the panel,
    /// projected into the camera's own basis.
    #[test]
    fn the_camera_frames_the_volume_and_the_panel() {
        // Bevy's `Camera3d` default is a 45-degree *vertical* field of view.
        let tan_v = (std::f32::consts::FRAC_PI_8).tan();
        let tan_h = tan_v * 16.0 / 9.0;

        for field in 0..FIELD_COUNT {
            let mut app = harness(field);
            step(&mut app);
            app.world_mut().spawn(OrbitCamera {
                focus: Vec3::ZERO,
                radius: 1.0,
                yaw: VIEW_YAW,
                pitch: VIEW_PITCH,
            });
            app.world_mut()
                .run_system_once(frame_camera)
                .expect("the framing system");

            let (focus, radius, yaw, pitch) = {
                let mut cameras = app.world_mut().query::<&OrbitCamera>();
                let orbit = cameras.iter(app.world()).next().expect("the orbit camera");
                (orbit.focus, orbit.radius, orbit.yaw, orbit.pitch)
            };
            let grid = app.world().resource::<Grid>();
            let sampled = grid.0.as_ref().expect("the sampled field");
            let name = sampled.name;

            // `orbit_camera`'s own expressions, so this checks the framing rather
            // than a second idea of where the camera is.
            let dir = Vec3::new(
                yaw.cos() * pitch.cos(),
                pitch.sin(),
                yaw.sin() * pitch.cos(),
            );
            let eye = focus + dir * radius;
            let forward = -dir;
            let right = forward.cross(Vec3::Y).normalize_or_zero();
            let up = right.cross(forward).normalize_or_zero();

            let (panel_lo, panel_hi) = panel_bounds(sampled);
            let mut points = Vec::new();
            for corner in 0..8usize {
                points.push(Vec3::new(
                    if corner & 1 == 0 {
                        sampled.domain_min.x
                    } else {
                        sampled.domain_max.x
                    },
                    if corner & 2 == 0 {
                        sampled.domain_min.y
                    } else {
                        sampled.domain_max.y
                    },
                    if corner & 4 == 0 {
                        sampled.domain_min.z
                    } else {
                        sampled.domain_max.z
                    },
                ));
            }
            for corner in 0..4usize {
                points.push(Vec3::new(
                    if corner & 1 == 0 {
                        panel_lo.x
                    } else {
                        panel_hi.x
                    },
                    if corner & 2 == 0 {
                        panel_lo.y
                    } else {
                        panel_hi.y
                    },
                    panel_lo.z,
                ));
            }

            for p in points {
                let v = p - eye;
                let depth = v.dot(forward);
                assert!(
                    depth > 0.0,
                    "{name}: {p:?} is behind the eye at {eye:?} -- the camera is inside its \
                     own subject"
                );
                let (x, y) = (v.dot(right) / depth, v.dot(up) / depth);
                assert!(
                    x.abs() <= tan_h && y.abs() <= tan_v,
                    "{name}: {p:?} projects to ({x:.3}, {y:.3}) against a frustum of \
                     ({tan_h:.3}, {tan_v:.3}) at radius {radius:.2} -- part of the subject is \
                     off screen at 1280x720"
                );
            }
        }
    }
}
