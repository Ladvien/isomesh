//! E-319 — a metric built from the trilinear's own Hessian, beside uniform
//! refinement at matched Hausdorff, and it does not pay.
//!
//! ```bash
//! cd bevy_isomesh && cargo run --example anisotropic_metric --release
//! ```
//!
//! Keys: `1`–`7` pick a field (the eighth, `noise_cavity`, is reachable only as
//! `ISOMESH_FIELD=7`), `[` and `]` move which rung of the ladder is on screen,
//! `W` wireframes both arms, `G` outlines both domains, `H` hides the HUD.
//!
//! Under `ISOMESH_CAPTURE` it needs no keyboard: the demo walks the eight
//! reference fields, ten captured frames each, so `record_gif.sh`'s default 80
//! frames is exactly one pass over the roster and the clip loops. `ISOMESH_FIELD`
//! pins one field and turns the clip into a still of that field;
//! `ISOMESH_SAMPLES` pins the shown rung and must name one of `17`, `25`, `33`.
//!
//! ```bash
//! scripts/record_gif.sh anisotropic_metric docs/gifs/e319.gif
//! ```
//!
//! No clip size is quoted above because none was measured — this file was
//! written on a host with no display, so nothing here has been run with a
//! window. `record_gif.sh:99-105` warns outside 0.7–4.8 MB and is the check.
//!
//! # What is on screen
//!
//! Two meshes of one field, viewed face-on.
//!
//! - **Left, grey — the uniform arm.** `N` samples on every axis, spacing `h`
//!   equal in `x`, `y` and `z`. This is uniform refinement, and it is P-146 C1's
//!   baseline.
//! - **Right, tan — the metric arm.** The *same total sample budget* `N³`, spent
//!   on per-axis counts `(n_x, n_y, n_z)` proportional to the point densities
//!   `√(e_aᵀ M e_a)` that the optimal `L^p` metric prescribes along each world
//!   axis.
//!
//! A small wire box sits at each arm's low domain corner: that arm's grid cell,
//! drawn at its exact size rather than inflated, so the element the metric asked
//! for is visible beside the cube it replaced. On `thin_plate` and `fbm_terrain`
//! it is a pancake; on five of the eight fields the two boxes are the same cube,
//! which is the finding rather than a bug.
//!
//! **Nothing draws a per-cell metric ellipsoid, and that omission is
//! deliberate.** This arm has no per-cell anisotropy to draw, and a picture
//! implying otherwise would misrepresent what was measured — see the last
//! section.
//!
//! # What the number means
//!
//! `M_Lp = D_Lp · det(|H_u|)^(−1/(2p + d)) · |H_u|` with `d = 3`, `p = 2`, from
//! **NASA NTRS 20200003084**, which restates Loseille & Alauzet verbatim (the two
//! SIAM originals are paywalled). `|H_u|` is the Hessian of the field with the
//! absolute value taken on its spectrum, differenced at cell size on the same
//! nineteen-point stencil as `crates/isomesh/benches/common/metric.rs:440`. Every
//! number in the panel is recomputed here, live, through the public API; the
//! three lines that name a `P-` id are quoted from the committed CSVs and are
//! marked as citations.
//!
//! - **Triangles and symmetric Hausdorff, per arm, at the shown rung.** The
//!   Hausdorff comes from `validate::accuracy`, and **both arms are graded on the
//!   uniform arm's seed lattice** — `accuracy` explicitly licenses a seed lattice
//!   that is not the extraction grid (`validate/accuracy.rs:332-337`), and a
//!   comparison in which each arm grades its own homework on its own lattice is
//!   not a comparison.
//! - **The matched-Hausdorff read-off.** A triangle count at a *different* error
//!   is not a comparison either, so both arms are fitted `ln T` against `ln E`
//!   over the ladder and read off at one error `E*` — the uniform arm's finest.
//!   The ratio is `T_metric / T_uniform` there. **P-146 C1 asked for `≤ 0.75`.**
//! - **Aspect ratio, always beside an at-floor cell count.** `metric_lp` floors
//!   each `|eigenvalue|` at `H_FLOOR = 1e-9` before forming the determinant,
//!   because a genuinely flat direction would otherwise send `det^(−1/(2p+3))` to
//!   infinity. Where a direction *is* flat, `aspect_ratio` returns
//!   `|λ|max / H_FLOOR ≈ 1e9`–`1e11`, which is **a restatement of that constant
//!   and not a measured anisotropy** (`benches/common/metric.rs:67-74`). So the
//!   panel never shows a maximum without the count of cells sitting at the floor,
//!   and shows the maximum over the cells that are *not* at the floor on the next
//!   line. The module's own measured counts say why this matters: `box_exact`
//!   1686/1790, `fbm_terrain` 1156/1156 — a heightfield SDF is exactly linear in
//!   `y`, so `∂²f/∂y²` is identically zero — and `gyroid` 1/2945 with a genuine
//!   `5.11e3` (`benches/experiment_p146.rs:216`).
//! - **Metric build over extraction, twice over.** What is timed is `hessian` +
//!   `metric_lp` over the band points and nothing else; band *selection* is not
//!   timed, because the extractor already visits every sample and the band test
//!   is one comparison on a value it already has. **P-146 C3 asked for
//!   `≤ 0.15`.** The clock is a median of five repeats after a warm-up, reported
//!   with the span the repeats spanned and with a word saying whether that span
//!   straddles the bar — `M-280` measured a 1.45× swing between runs of one
//!   binary on this host's governor, so a share averaged into a pass is not a
//!   measurement. Beside it sits the same statement as a **count**: `19` field
//!   samples per band point against the grid's `N³`, which no governor can move.
//!
//!   **Run it `--release`, as the command above says.** Measured on this host,
//!   the release build's shares reproduce p-146's medians across the roster —
//!   `1.924` against `1.952995` on `sphere`, `0.681` against `0.695277` on
//!   `thin_plate` — while the test profile's do not, because
//!   `bevy_isomesh/Cargo.toml:179-180` builds the local crate at `opt-level = 1`
//!   and that lands on the extraction in the share's *denominator* rather than on
//!   the Hessian arithmetic in its numerator. It is also why the `#[cfg(test)]`
//!   gate below is on the count and not on the clock.
//!
//! # What this shows, and what it does not
//!
//! **Both of P-146's headline clauses were falsified, and this demo exists to
//! show that honestly rather than to sell the mechanism.**
//!
//! - **P-146 C1: FALSIFIED.** No measurable field saved 25% of its triangles at
//!   matched Hausdorff — `c1_winners = 0` of a `c1_population = 4`. The measured
//!   `ratio` column runs `sphere 1.000000`, `box_exact 1.000000`,
//!   `torus 1.153211`, `thin_plate 2.298935`: on the one field with real
//!   anisotropy to spend, the metric arm needed **2.3× more** triangles at the
//!   same error, not fewer.
//! - **P-146 C2 is `unmeasurable` on 20 of 40 rows**, and that is not a gap in
//!   the harness. `validate::accuracy` projects along `∇f` and compares against
//!   `|f|`, and `fields/mod.rs:83-84` is explicit that only `FieldBound::Exact`
//!   *"admits a Hausdorff measurement against the field's own values"*. Four of
//!   the eight reference fields are `Exact`; the other four carry the reason
//!   string in place of a distance. C1's population is therefore 4, not 8.
//! - **P-146 C3: FALSIFIED at every rung.** The band is a shell (`∝ N²`) and
//!   extraction is a volume (`∝ N³`), so the share falls as `1/N` — but one band
//!   point costs nineteen field samples plus a Jacobi eigendecomposition against
//!   an extraction cell's roughly one sample, and the measured share runs from
//!   `0.55` to `8.75`. The bar was `0.15`.
//! - **P-147: all three clauses FALSIFIED.** The theory says anisotropy buys the
//!   *constant*, not the exponent. Neither behaved: `exponent_difference` is
//!   `0.108109` on `thin_plate` and `0.352013` on `fbm_terrain` against C1's bar
//!   of `0.1`, and the constant moved the *wrong way* where the exponent held —
//!   `constant_ratio 3.131972` on `thin_plate`. C2's rank correlation between the
//!   constant's improvement and the AM–GM gap the theory names is `-0.095238`.
//! - **P-149 found the mechanism absent, not weak.** 70 of its 112 rows carry
//!   `arms_identical = true`: on `sphere`, `box_exact`, `csg_difference`,
//!   `gyroid` and `noise_cavity` the metric prescribes **the same grid** the
//!   uniform arm already used, so there is no anisotropic arm on those fields at
//!   all. Its pooled correlation is `0.187131` against a bar of `0.7`.
//!
//! **And the one caveat that governs every verdict above.** What is built here —
//! in the bench and in this file — is a **per-axis global** anisotropic sampling
//! grid, **not per-cell anisotropy**. `crates/isomesh/src/` was frozen for Phase
//! 27 and every inherent `extract` in the crate takes a single *scalar*
//! `cell_size` (`marching_cubes/mod.rs:193`), so a per-cell anisotropic mesher is
//! a source change by construction rather than by choice. The consequence is
//! stated in the bench header (`benches/experiment_p146.rs:54-60`) and is
//! repeated here because it is the single most important sentence in either file:
//! **a field whose flat direction *rotates* over the surface has per-cell
//! anisotropy this arm structurally cannot spend, and it reads as a null here
//! while a real anisotropic mesher would win on it.** The `flat axis` line in the
//! panel is that mechanism, measured: it is `y` only on `fbm_terrain`, whose
//! `sample` is exactly linear in `y`. Everywhere else it is `mixed` or `none`,
//! and a `mixed` field is one this construction cannot help. Every verdict on
//! screen is a verdict about *this* construction and about this roster of eight
//! fields — not about metric-based adaptation.

mod common;

use std::f32::consts::FRAC_PI_2;
use std::time::Instant;

use bevy::asset::RenderAssetUsages;
use bevy::mesh::{Indices, PrimitiveTopology};
use bevy::prelude::*;
use common::{Capture, CommonPlugin, DemoDomain, DemoMesh, DemoStats, OrbitCamera, ViewFlags};
use isomesh::fields::{
    BoxExact, FbmTerrain, FieldBound, ReferenceField, Sphere, ThinPlate, Torus, capped_gyroid,
    csg_difference, noise_cavity,
};
use isomesh::marching_cubes::MarchingCubes;
use isomesh::validate::{AccuracyConfig, accuracy};
use isomesh::{MeshBuffer, Real, RuntimeShape3, Sdf};

// ─── the constants, and where each one comes from ────────────────────────────

/// The demo's own ladder, in samples per axis.
///
/// Three rungs rather than P-146's five (`17, 25, 33, 49, 65`), because the two
/// `validate::accuracy` calls per rung are the cost here and a demo must rebuild
/// inside a frame. **So the read-off below need not equal p-146's `ratio` column
/// to the digit** — on the four measurable fields this ladder gives `1.0000`,
/// `1.1959`, `1.0000` and `3.5386` against the committed `1.000000`,
/// `1.153211`, `1.000000` and `2.298935`. The direction and the verdict are the
/// same on every one; the magnitude differs because a shorter ladder extrapolates
/// further, which the `read off by` line on screen says outright.
///
/// Odd throughout: `M-266` proved the canonical grids' odd counts are
/// load-bearing — `thin_plate` is centred on `y = 0` and loses its surface
/// entirely on an even count, which would be a parity result wearing a metric's
/// clothes.
const LADDER: [u32; 3] = [17, 25, 33];

/// The rung the demo opens on: the finest, where the fitted `E*` is read.
const DEFAULT_RUNG: usize = LADDER.len() - 1;

/// The norm `M_Lp` is optimised for. `L²`, p-146's `p_norm` column.
const P_NORM: f64 = 2.0;

/// The global Lagrange constant of the optimal `L^p` metric, folded to 1.
///
/// It multiplies the whole metric field by one scalar and cancels in every ratio
/// reported here, so recording it as a number would be recording a choice of
/// units (`benches/common/metric.rs:46-53`).
const D_LP: f64 = 1.0;

/// The absolute floor applied to each `|Hessian eigenvalue|` before the
/// `det(|H|)^(−1/(2p+d))` factor is formed.
///
/// `1e-9`, in the units of a second derivative. The reference fields' genuine
/// curvatures run from about `1e-2` to about `1e2`, so this sits seven orders
/// below the smallest real one and cannot mask a curvature. **An `aspect_ratio`
/// near `|λ|max / H_FLOOR` is reporting this constant and not a measurement**,
/// which is why the panel never prints a maximum without an at-floor count.
const H_FLOOR: f64 = 1e-9;

/// Maximum cyclic Jacobi sweeps. Three or four is normal for a `3 × 3`; twelve
/// is a cap that makes termination unconditional, not an expectation.
const JACOBI_SWEEPS: usize = 12;

/// Early-exit tolerance for the Jacobi sweep, relative to `‖M‖_F`.
///
/// The Frobenius norm and deliberately not the trace: a Hessian is indefinite
/// and `diag(1, −1, 0)` is an ordinary saddle whose trace is zero, against which
/// the tolerance would be identically zero and the early exit could never fire.
const JACOBI_TOLERANCE: f64 = 1e-14;

/// A component of an eigenvector counts as non-zero for the sign convention
/// above this. Columns are unit length, so the test is already relative.
const SIGN_EPS: f64 = 1e-12;

/// `(i, j)` → index into [`Sym3`]'s six-entry array.
const IDX: [[usize; 3]; 3] = [[0, 1, 2], [1, 3, 4], [2, 4, 5]];

/// A grid sample joins the surface band when `|f| <= BAND_CELLS · h`.
///
/// One cell, not three: the metric is wanted where the mesh is, and a wider
/// shell would charge C3 for points no triangle is near.
const BAND_CELLS: f64 = 1.0;

/// Field samples one [`hessian`] costs: seven on the diagonal stencils sharing
/// one centre, plus twelve corners.
///
/// A **count**, and the only machine-independent half of C3. `19 · band` against
/// the grid's `N³` is the same statement the clock makes, exactly reproducible,
/// and it is a *floor* on the metric's cost because it does not charge the one
/// Jacobi eigendecomposition per band point.
const HESSIAN_SAMPLES: u64 = 19;

/// Timed repeats per stage, after one warm-up.
///
/// Five, and reported as a **median** rather than a mean, because `M-280`
/// measured a 1.45x swing between runs of the same binary on this machine's
/// `amd-pstate-epp` governor and one slow repeat must not move the headline. The
/// min and max are carried beside it, so the panel can say whether the spread
/// straddles C3's bar instead of averaging into a verdict.
///
/// **A single shot is not enough, and that is measured rather than assumed.** The
/// first version of this file timed one call each. On `thin_plate` at `33³` the
/// unwarmed extraction read `1.342 ms` while the metric build read `0.198`
/// against p-146's `0.213`, so the share came out `0.148` and the panel printed
/// `HOLDS` for a clause p-146 falsified on all forty of its rows. The whole error
/// was in the share's *denominator*, and a demo claiming a clause held when the
/// measurement says otherwise is worse than no demo. With the warm-up and the
/// median the same row reads `0.681` in a release build against p-146's
/// `0.695277`.
///
/// Timed at the shown rung only: the other rungs contribute triangles and
/// Hausdorff, which is all the read-off reads, and five repeats of every rung
/// would put a fifth of a second between a keypress and a picture on
/// `fbm_terrain`.
const REPEATS: usize = 5;

/// Fewest samples any axis of the metric grid may carry.
///
/// Odd, and five rather than two: `Error::GridTooSmall` starts at two, and a
/// two-sample axis has one cell and can only place a crossing against the wall.
const MIN_SAMPLES: u32 = 5;

/// P-146 C1's bar: "at least 25% fewer triangles at matched Hausdorff".
const WIN_RATIO: f64 = 0.75;

/// P-146 C3's bar: the metric must cost under 15% of extraction.
const C3_MAX_SHARE: f64 = 0.15;

/// The ladder must move the error instrument by at least this factor on both
/// arms, or a fit through it has no slope to read.
const LADDER_SPAN_FLOOR: f64 = 1.2;

/// `cos(5°)`. A floored eigenvector counts as axis-aligned above this.
const AXIS_ALIGNED_COS: f64 = 0.996_194_698_091_745_5;

/// Captured frames spent on each field.
///
/// Eight fields at ten frames is eighty, which is `record_gif.sh:47`'s default
/// `ISOMESH_CAPTURE_FRAMES`, so the clip is exactly one pass over the roster and
/// ends where it started.
const CAPTURE_FRAMES_PER_FIELD: u32 = 10;

/// The uniform arm's surface, grey. The same pair of hues E-103 uses, so the
/// crate's two side-by-side demos read alike.
const UNIFORM_SURFACE: Color = Color::srgb(0.78, 0.80, 0.86);

/// The metric arm's surface, tan.
const METRIC_SURFACE: Color = Color::srgb(0.86, 0.74, 0.55);

/// The uniform arm's grid-cell outline, a brightened form of its surface so the
/// cell is attributable to the mesh beneath it without a caption.
const UNIFORM_CELL: Color = Color::srgb(0.60, 0.75, 1.00);

/// The metric arm's grid-cell outline.
const METRIC_CELL: Color = Color::srgb(1.00, 0.78, 0.35);

/// A falsified clause's headline colour.
const FALSIFIED: Color = Color::srgb(1.00, 0.45, 0.35);

/// An unmeasurable field's headline colour: neither a pass nor a failure.
const UNMEASURABLE: Color = Color::srgb(1.00, 0.80, 0.30);

/// A clause that held, which no field on this roster reaches.
const HELD: Color = Color::srgb(0.35, 0.95, 0.55);

// ─── the citations ───────────────────────────────────────────────────────────

/// One field's committed figures, copied verbatim from the CSVs.
///
/// Every string here is a **citation**, not a number this file computed: the
/// text is byte-for-byte the cell named in each field's doc comment, and the HUD
/// prints it *beside* the live figure rather than in place of it, which is the
/// house rule (`game_dig.rs:2946-2952`). Where a clause was unmeasurable the
/// CSV's own reason string is carried, because `validate::accuracy` compares
/// against `|f|` and only a `bound()` of `Exact` makes that a distance
/// (`fields/mod.rs:83-84`).
struct Cited {
    /// `ReferenceField::NAME`, and the CSV's `field` column.
    name: &'static str,
    /// `p-146.csv`, `ratio`: triangles at matched symmetric Hausdorff, metric arm
    /// over uniform arm, read off a five-rung ladder. C1 wanted `<= 0.75`.
    p146_ratio: &'static str,
    /// `p-146.csv`, `metric_share` on this field's `33` row: metric build over
    /// extraction, median of five after a warm-up. C3 wanted `<= 0.15`.
    p146_share: &'static str,
    /// `p-147.csv`, `exponent_difference`. C1 wanted `< 0.1` on every smooth
    /// field, the theory's claim being that anisotropy cannot buy an exponent.
    p147_exponent: &'static str,
    /// `p-147.csv`, `constant_ratio`: the metric arm's fitted constant over the
    /// uniform arm's. Below one is an improvement; above one is a penalty.
    p147_constant: &'static str,
}

/// The eight reference fields, in `for_each_reference_field!`'s own order
/// (`fields/mod.rs:212-255`), each beside its committed row.
///
/// The order is load-bearing twice over: it is the order the digit keys index,
/// and it is the order the dispatch in [`build`] matches on. [`measure`] asserts
/// that the field it was handed and the row it was handed agree by name, so a
/// mis-ordered table is a loud failure rather than a demo quoting one field's
/// CSV row against another field's mesh.
const CITED: [Cited; 8] = [
    Cited {
        name: "sphere",
        p146_ratio: "1.000000",
        p146_share: "1.952995",
        p147_exponent: "0.000000",
        p147_constant: "1.000000",
    },
    Cited {
        name: "torus",
        p146_ratio: "1.153211",
        p146_share: "1.755654",
        p147_exponent: "0.018566",
        p147_constant: "1.125249",
    },
    Cited {
        name: "box_exact",
        p146_ratio: "1.000000",
        p146_share: "2.888934",
        p147_exponent: "0.000000",
        p147_constant: "1.000000",
    },
    Cited {
        name: "csg_difference",
        p146_ratio: "unmeasurable:bound=Underestimate",
        p146_share: "4.595306",
        p147_exponent: "0.000000",
        p147_constant: "1.000000",
    },
    Cited {
        name: "thin_plate",
        p146_ratio: "2.298935",
        p146_share: "0.695277",
        p147_exponent: "0.108109",
        p147_constant: "3.131972",
    },
    Cited {
        name: "gyroid",
        p146_ratio: "unmeasurable:bound=Lipschitz",
        p146_share: "3.393888",
        p147_exponent: "0.000000",
        p147_constant: "1.000000",
    },
    Cited {
        name: "fbm_terrain",
        p146_ratio: "unmeasurable:bound=Unbounded",
        p146_share: "1.067020",
        p147_exponent: "0.352013",
        p147_constant: "0.114405",
    },
    Cited {
        name: "noise_cavity",
        p146_ratio: "unmeasurable:bound=Unbounded",
        p146_share: "2.513031",
        p147_exponent: "0.000000",
        p147_constant: "1.000000",
    },
];

/// `p-146.csv`: `c1_winners` reads `0` against a `c1_population` of `4` on all
/// forty rows, and `c3_holds` reads `false` on all forty.
const P146_GLOBAL: &str =
    "P-146 global: c1_winners 0 of population 4, c3_holds false on all 40 rows";

/// `p-147.csv`: `c1_holds`, `c2_holds` and `c3_holds` are `false` on all eight
/// rows, and `c2_rank_correlation` — the constant's improvement against the
/// AM-GM gap the theory names — is `-0.095238`.
const P147_GLOBAL: &str =
    "P-147 global: C1, C2 and C3 all false on 8 of 8 rows, rank corr. -0.095238";

/// `p-149.csv`: 70 of its 112 rows carry `arms_identical=true`, and
/// `pooled_correlation` is `0.187131` against C1's bar of `0.7`.
const P149_GLOBAL: &str =
    "P-149 global: 70 of 112 rows had IDENTICAL arms, pooled corr. 0.187131 (bar 0.7)";

// ─── the metric algebra ──────────────────────────────────────────────────────

/// A symmetric 3x3 matrix, stored as the six upper entries
/// `[xx, xy, xz, yy, yz, zz]`.
///
/// Six `f64`s rather than nine: a Hessian and a metric are symmetric by
/// construction, and storing the lower triangle separately invites the two
/// halves to disagree after a rotation. This is `benches/common/metric.rs:158`
/// reimplemented example-locally — that module is `pub(crate)` inside
/// `crates/isomesh/benches/` and nothing outside it can call it, and Phase 27
/// froze `crates/isomesh/src/`, so a demo of P-146 has to carry the four
/// operations it needs. The eigensolver is reproduced sweep for sweep and
/// tie-break for tie-break, because a divergence would make the numbers on
/// screen something other than p-146's.
///
/// The module's argument assertions are *not* carried across: `metric_lp`'s
/// `p > 0` and `hessian`'s `h > 0` guard against callers, and both arguments
/// here are consts of this file or a domain extent over a positive count.
#[derive(Clone, Copy, Debug)]
struct Sym3([f64; 6]);

impl Sym3 {
    /// Entry `(i, j)`, either triangle. Panics for an index outside `0..3`.
    fn get(&self, i: usize, j: usize) -> f64 {
        self.0[IDX[i][j]]
    }

    /// Jacobi eigendecomposition: eigenvalues ascending, eigenvectors as
    /// columns.
    ///
    /// `vectors[row][col]` is component `row` of eigenvector `col`, so
    /// `values[c]` pairs with column `c`. Eigenvalues come back **ascending by
    /// value, not by magnitude** — a saddle's smallest magnitude is not
    /// `values[0]`, which is why the census below searches for it.
    ///
    /// Ordering is tie-broken lexicographically on the eigenvector column and
    /// then on the original axis index, and each column is signed so its first
    /// component above [`SIGN_EPS`] is positive; without both, a degenerate
    /// spectrum produces whichever ordering the sweep happened to leave and two
    /// calls on numerically identical matrices can disagree.
    ///
    /// Panics on a non-finite entry: a NaN pivot would flow silently into every
    /// eigenvalue and every consumer, and a field that sampled NaN is a defect
    /// upstream of here rather than a case to be handled.
    fn eigen(&self) -> ([f64; 3], [[f64; 3]; 3]) {
        assert!(
            self.0.iter().all(|entry| entry.is_finite()),
            "Sym3::eigen on a non-finite matrix: {self:?}"
        );

        let [xx, xy, xz, yy, yz, zz] = self.0;
        let mut work = [[xx, xy, xz], [xy, yy, yz], [xz, yz, zz]];
        let mut basis = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];

        // ‖M‖_F is invariant under the sweep, so it is fixed once here.
        let frobenius = (xx * xx + yy * yy + zz * zz + 2.0 * (xy * xy + xz * xz + yz * yz)).sqrt();
        let tolerance = JACOBI_TOLERANCE * frobenius;

        for _ in 0..JACOBI_SWEEPS {
            let off = (2.0
                * (work[0][1] * work[0][1] + work[0][2] * work[0][2] + work[1][2] * work[1][2]))
                .sqrt();
            if off <= tolerance {
                break;
            }
            for (row, col) in [(0usize, 1usize), (0, 2), (1, 2)] {
                let pivot = work[row][col];
                if pivot.abs() <= 0.0 {
                    continue;
                }

                // The rotation that annihilates `work[row][col]`. Taking the
                // smaller root keeps |t| <= 1, which is what makes the sweep
                // unconditionally stable; a huge theta sends `t` to zero, i.e.
                // no rotation, which is the correct limit reached without a
                // branch.
                let diff = work[col][col] - work[row][row];
                let theta = diff / (2.0 * pivot);
                let sign = if theta < 0.0 { -1.0 } else { 1.0 };
                let tan = sign / (theta.abs() + (theta * theta + 1.0).sqrt());
                let cos = 1.0 / (tan * tan + 1.0).sqrt();
                let sin = tan * cos;
                let half_tan = sin / (1.0 + cos);

                work[row][row] -= tan * pivot;
                work[col][col] += tan * pivot;
                work[row][col] = 0.0;
                work[col][row] = 0.0;

                // The third index: (0,1)->2, (0,2)->1, (1,2)->0.
                let other = 3 - row - col;
                let from_row = work[other][row];
                let from_col = work[other][col];
                let to_row = from_row - sin * (from_col + half_tan * from_row);
                let to_col = from_col + sin * (from_row - half_tan * from_col);
                work[other][row] = to_row;
                work[row][other] = to_row;
                work[other][col] = to_col;
                work[col][other] = to_col;

                for basis_row in &mut basis {
                    let old_row = basis_row[row];
                    let old_col = basis_row[col];
                    basis_row[row] = cos * old_row - sin * old_col;
                    basis_row[col] = sin * old_row + cos * old_col;
                }
            }
        }

        let raw = [work[0][0], work[1][1], work[2][2]];

        let mut order = [0usize, 1, 2];
        order.sort_by(|&left, &right| {
            raw[left].total_cmp(&raw[right]).then_with(|| {
                for basis_row in &basis {
                    let cmp = basis_row[left].total_cmp(&basis_row[right]);
                    if cmp != std::cmp::Ordering::Equal {
                        return cmp;
                    }
                }
                left.cmp(&right)
            })
        });

        let values = [raw[order[0]], raw[order[1]], raw[order[2]]];
        let mut vectors = [[0.0f64; 3]; 3];
        for (col, &src) in order.iter().enumerate() {
            for (row, vector_row) in vectors.iter_mut().enumerate() {
                vector_row[col] = basis[row][src];
            }
        }

        for col in 0..3 {
            let mut lead = 0usize;
            while lead < 2 && vectors[lead][col].abs() <= SIGN_EPS {
                lead += 1;
            }
            if vectors[lead][col] < 0.0 {
                for row in &mut vectors {
                    row[col] = -row[col];
                }
            }
        }

        (values, vectors)
    }

    /// Rebuild from eigenvalues and eigenvectors: `M = Σ_c λ_c v_c v_cᵀ`.
    ///
    /// `vectors[row][col]` is component `row` of eigenvector `col`, matching
    /// [`Sym3::eigen`]'s output exactly.
    fn from_eigen(values: [f64; 3], vectors: [[f64; 3]; 3]) -> Self {
        let mut out = [0.0f64; 6];
        for (col, &value) in values.iter().enumerate() {
            for i in 0..3 {
                for j in i..3 {
                    out[IDX[i][j]] += value * vectors[i][col] * vectors[j][col];
                }
            }
        }
        Self(out)
    }

    /// max/min `|eigenvalue|`, i.e. the anisotropy of the element the metric
    /// prescribes; `INFINITY` if the minimum is zero.
    ///
    /// The element's edge lengths go as `1/√λ`, so the *length* ratio is the
    /// square root of this. The zero matrix prescribes nothing and reports
    /// `INFINITY` rather than `NaN` — a degenerate metric is infinitely
    /// anisotropic, not undefined.
    fn aspect_ratio(&self) -> f64 {
        let (values, _) = self.eigen();
        let mut lo = f64::INFINITY;
        let mut hi = 0.0f64;
        for value in values {
            let magnitude = value.abs();
            if magnitude < lo {
                lo = magnitude;
            }
            if magnitude > hi {
                hi = magnitude;
            }
        }
        if lo > 0.0 { hi / lo } else { f64::INFINITY }
    }
}

/// Central-difference Hessian of `sdf` at `p`, with step `h` (the cell size).
///
/// Diagonal on the 7-point stencil `(f(p + h eᵢ) − 2 f(p) + f(p − h eᵢ)) / h²`,
/// mixed entries on the 4-point stencil `(f⁺⁺ − f⁺⁻ − f⁻⁺ + f⁻⁻) / 4h²`.
/// Nineteen samples in all, and the same stencil at the same step the crate
/// already differences the trilinear at (`M-65`).
///
/// Both stencils have identically zero truncation error on a quadratic, so `h`
/// should be the **cell size** rather than something much smaller: shrinking it
/// does not reduce a truncation error that is already zero, it only amplifies
/// the cancellation in the subtraction. That is also why every number in this
/// file is `f64` — a second difference divided by `h²` is the one quantity in
/// this crate that cannot afford `f32`, and the narrowing to `f32` happens once,
/// in [`bevy_mesh`], on positions that are going to a GPU.
fn hessian<S: Sdf<Scalar = f64>>(sdf: &S, p: [f64; 3], h: f64) -> Sym3 {
    let at = |offset: [f64; 3]| sdf.sample([p[0] + offset[0], p[1] + offset[1], p[2] + offset[2]]);
    let centre = at([0.0; 3]);
    let inv_h2 = 1.0 / (h * h);

    let mut out = [0.0f64; 6];

    for (axis, slot) in [(0usize, 0usize), (1, 3), (2, 5)] {
        let mut plus = [0.0; 3];
        plus[axis] = h;
        let mut minus = [0.0; 3];
        minus[axis] = -h;
        out[slot] = (at(plus) - 2.0 * centre + at(minus)) * inv_h2;
    }

    for (i, j, slot) in [(0usize, 1usize, 1usize), (0, 2, 2), (1, 2, 4)] {
        let corner = |si: f64, sj: f64| {
            let mut offset = [0.0; 3];
            offset[i] = si * h;
            offset[j] = sj * h;
            at(offset)
        };
        out[slot] = (corner(1.0, 1.0) - corner(1.0, -1.0) - corner(-1.0, 1.0) + corner(-1.0, -1.0))
            * (0.25 * inv_h2);
    }

    Sym3(out)
}

/// The optimal `L^p` metric: `D_Lp · det(|H|)^(−1/(2p + d)) · |H|`, `d = 3`.
///
/// Each `|eigenvalue|` is floored at [`H_FLOOR`] **before** the determinant is
/// formed, so a flat direction produces a very anisotropic metric rather than an
/// infinity. The sign of a curvature does not make an element cheaper, so only
/// the magnitude survives into the metric.
fn metric_lp(hessian: &Sym3) -> Sym3 {
    let (values, vectors) = hessian.eigen();
    let floored = [
        values[0].abs().max(H_FLOOR),
        values[1].abs().max(H_FLOOR),
        values[2].abs().max(H_FLOOR),
    ];
    let determinant = floored[0] * floored[1] * floored[2];
    let exponent = -1.0 / (2.0 * P_NORM + 3.0);
    let factor = D_LP * determinant.powf(exponent);

    Sym3::from_eigen(
        [
            floored[0] * factor,
            floored[1] * factor,
            floored[2] * factor,
        ],
        vectors,
    )
}

// ─── the metric field, and what it says ──────────────────────────────────────

/// The grid samples within one cell of the surface, swept `z`, `y`, `x` with
/// `x` innermost — the crate's order, so the result is reproducible.
fn band_points<F>(field: &F, origin: [f64; 3], h: f64, samples: u32) -> Vec<[f64; 3]>
where
    F: Sdf<Scalar = f64>,
{
    let band = BAND_CELLS * h;
    let mut out = Vec::new();
    for k in 0..samples {
        for j in 0..samples {
            for i in 0..samples {
                let p = [
                    origin[0] + f64::from(i) * h,
                    origin[1] + f64::from(j) * h,
                    origin[2] + f64::from(k) * h,
                ];
                if field.sample(p).abs() <= band {
                    out.push(p);
                }
            }
        }
    }
    out
}

/// What the metric field says about one rung — all of it untimed, because none
/// of it is part of *building* the metric and none of it may be charged to C3.
struct Census {
    /// Band points, and the denominator every other count here needs.
    points: usize,
    /// Points whose smallest `|Hessian eigenvalue|` is at or below [`H_FLOOR`].
    /// Their aspect ratio is the floor talking.
    at_floor: usize,
    /// Points with a floored eigenvalue, at least one un-floored one, and the
    /// floored eigenvector within 5° of a world axis.
    ///
    /// The middle condition is what makes a flat direction *exploitable*: on a
    /// box face all three eigenvalues are floored, the point is flat in every
    /// direction and it prescribes nothing. Counting it would make every box
    /// look like a heightfield.
    exploitable_flat: usize,
    /// Which world axis each exploitable flat direction pointed along.
    axis_hits: [usize; 3],
    /// Maximum and mean aspect ratio over every band point.
    aspect: [f64; 2],
    /// Maximum and mean over the points that are **not** at the floor. Zero when
    /// there are none, which is itself the reading.
    aspect_off_floor: [f64; 2],
    /// Mean `√(e_aᵀ M e_a)` per axis: the metric's own point density along each
    /// world axis, and the only thing the anisotropic split is derived from.
    weights: [f64; 3],
}

/// Census the metric field.
///
/// Recomputes the Hessian because the at-floor test is a question about `H`, not
/// about `M`: [`metric_lp`] scales every eigenvalue by one common factor, so the
/// floor is invisible in the metric alone.
fn census_of<F>(field: &F, points: &[[f64; 3]], metrics: &[Sym3], h: f64) -> Census
where
    F: Sdf<Scalar = f64>,
{
    let mut census = Census {
        points: points.len(),
        at_floor: 0,
        exploitable_flat: 0,
        axis_hits: [0; 3],
        aspect: [0.0; 2],
        aspect_off_floor: [0.0; 2],
        weights: [0.0; 3],
    };
    let mut sum = 0.0f64;
    let mut sum_off_floor = 0.0f64;
    let mut off_floor = 0usize;

    for (&p, m) in points.iter().zip(metrics) {
        let (values, vectors) = hessian(field, p, h).eigen();

        let mut flat = 0usize;
        let mut stiff = 0usize;
        for (index, value) in values.iter().enumerate() {
            if value.abs() < values[flat].abs() {
                flat = index;
            }
            if value.abs() > H_FLOOR {
                stiff += 1;
            }
        }
        let floored = values[flat].abs() <= H_FLOOR;

        let aspect = m.aspect_ratio();
        sum += aspect;
        if aspect > census.aspect[0] {
            census.aspect[0] = aspect;
        }
        if floored {
            census.at_floor += 1;
        } else {
            off_floor += 1;
            sum_off_floor += aspect;
            if aspect > census.aspect_off_floor[0] {
                census.aspect_off_floor[0] = aspect;
            }
        }

        if floored && stiff > 0 {
            let vector = [vectors[0][flat], vectors[1][flat], vectors[2][flat]];
            let mut axis = 0usize;
            for (index, component) in vector.iter().enumerate() {
                if component.abs() > vector[axis].abs() {
                    axis = index;
                }
            }
            if vector[axis].abs() >= AXIS_ALIGNED_COS {
                census.exploitable_flat += 1;
                census.axis_hits[axis] += 1;
            }
        }

        for (axis, weight) in census.weights.iter_mut().enumerate() {
            *weight += m.get(axis, axis).sqrt();
        }
    }

    let n = points.len() as f64;
    census.aspect[1] = sum / n;
    if off_floor > 0 {
        census.aspect_off_floor[1] = sum_off_floor / off_floor as f64;
    }
    for weight in &mut census.weights {
        *weight /= n;
    }
    census
}

/// The dominant axis of the exploitable flat directions, or `mixed` / `none`.
///
/// "Dominant" is a two-thirds majority: a field whose flat directions spread
/// over two or three axes has no axis for the grid to spend a budget on, and
/// calling the largest bucket the answer would launder that into a mechanism.
///
/// **Read it together with the axis-aligned fraction on the same line.** A
/// two-thirds majority of a population of one is still a two-thirds majority, so
/// a field with six exploitable band points out of eighteen thousand can be
/// labelled `y` while having no flat structure whatever. The fraction is the
/// column that says so.
fn flat_axis(census: &Census) -> &'static str {
    if census.exploitable_flat == 0 {
        return "none";
    }
    let names = ["x", "y", "z"];
    for (axis, &hits) in census.axis_hits.iter().enumerate() {
        if hits * 3 >= census.exploitable_flat * 2 {
            return names[axis];
        }
    }
    "mixed"
}

// ─── the metric-driven grid ──────────────────────────────────────────────────

/// Nearest odd integer, at least one. Ties go up, deterministically.
fn round_odd(x: f64) -> u32 {
    let half = ((x - 1.0) * 0.5).round();
    (2.0f64.mul_add(half, 1.0)).max(1.0) as u32
}

/// Per-axis sample counts from the metric's per-axis point densities, at the
/// uniform arm's total budget.
///
/// `n_a ∝ weights[a]` with `∏ n_a = N³`, so the constant is
/// `k = (N³)^{1/3} / geomean(weights)` and the two arms differ in **shape**
/// only. That is what makes a triangle-count comparison between them mean
/// anything, and the `budget` figure on screen is the check.
///
/// **There is exactly one clamp and it is a lower one.** An axis that would fall
/// below [`MIN_SAMPLES`] is pinned there and the remaining budget is *re-solved*
/// over the axes still free, which is what keeps `∏ n_a` on the budget instead
/// of merely near it. At most three rounds, because every round that changes
/// anything pins at least one more axis. No upper clamp: P-146's first version
/// had one and it double-counted the budget by binding on two axes in the same
/// round, measured 5.77× over on `fbm_terrain`
/// (`benches/experiment_p146.rs:604-617`).
///
/// Returns the counts and how many axes were pinned at the floor.
fn metric_grid(weights: [f64; 3], samples: u32) -> ([u32; 3], usize) {
    let budget = f64::from(samples).powi(3);
    let mut pinned = [false; 3];
    let mut n = [samples; 3];

    for _round in 0..3 {
        let free: Vec<usize> = (0..3).filter(|&axis| !pinned[axis]).collect();
        if free.is_empty() {
            break;
        }
        let mut held = 1.0f64;
        for (axis, &fixed) in pinned.iter().enumerate() {
            if fixed {
                held *= f64::from(n[axis]);
            }
        }
        let count = free.len() as f64;
        let target = budget / held;
        let mut logsum = 0.0f64;
        for &axis in &free {
            logsum += weights[axis].ln();
        }
        let geometric_mean = (logsum / count).exp();
        let scale = target.powf(1.0 / count) / geometric_mean;

        let mut newly_pinned = false;
        for &axis in &free {
            let raw = round_odd(scale * weights[axis]);
            if raw < MIN_SAMPLES {
                n[axis] = MIN_SAMPLES;
                pinned[axis] = true;
                newly_pinned = true;
            } else {
                n[axis] = raw;
            }
        }
        if !newly_pinned {
            break;
        }
    }

    (n, pinned.iter().filter(|fixed| **fixed).count())
}

/// The field seen through a per-axis coordinate stretch.
///
/// `sample(q) = f(lo + q ⊙ s)`. Extracting this on a **cubic** grid of
/// `cell_size = h` in `q` is exactly extracting `f` on a rectilinear grid whose
/// physical spacings are `h · s`, which is the only way to reach a rectilinear
/// anisotropic grid through an `extract` that takes a scalar `cell_size`
/// (`marching_cubes/mod.rs:193`).
///
/// `gradient` is left as `Sdf`'s central-difference default, so the emitted
/// normals are in `q` space. [`measure`] maps both the positions and the normals
/// back to world space before anything measures or renders them.
struct Stretched<'a, F> {
    /// The field being stretched.
    field: &'a F,
    /// World position of `q = 0`.
    lo: [f64; 3],
    /// World length of one unit of `q`, per axis.
    s: [f64; 3],
}

impl<F> Sdf for Stretched<'_, F>
where
    F: Sdf<Scalar = f64>,
{
    type Scalar = f64;

    fn sample(&self, q: [f64; 3]) -> f64 {
        self.field.sample([
            q[0].mul_add(self.s[0], self.lo[0]),
            q[1].mul_add(self.s[1], self.lo[1]),
            q[2].mul_add(self.s[2], self.lo[2]),
        ])
    }
}

// ─── the matched-error read-off ──────────────────────────────────────────────

/// Least squares of `ln(triangles)` against `ln(error)`.
struct Fit {
    /// `ln T` at `ln E = 0`.
    intercept: f64,
    /// `d ln T / d ln E`.
    slope: f64,
}

impl Fit {
    /// The fitted triangle count at error `e`.
    fn eval(&self, e: f64) -> f64 {
        self.slope.mul_add(e.ln(), self.intercept).exp()
    }
}

/// Fit one arm's ladder. The caller has already established that the errors span
/// [`LADDER_SPAN_FLOOR`], so the slope's denominator is non-zero.
fn fit_log(points: &[(f64, f64)]) -> Fit {
    let n = points.len() as f64;
    let mut mean_x = 0.0f64;
    let mut mean_y = 0.0f64;
    for &(error, triangles) in points {
        mean_x += error.ln();
        mean_y += triangles.ln();
    }
    mean_x /= n;
    mean_y /= n;

    let mut sxx = 0.0f64;
    let mut sxy = 0.0f64;
    for &(error, triangles) in points {
        let dx = error.ln() - mean_x;
        sxx = dx.mul_add(dx, sxx);
        sxy = dx.mul_add(triangles.ln() - mean_y, sxy);
    }
    let slope = sxy / sxx;

    Fit {
        intercept: slope.mul_add(-mean_x, mean_y),
        slope,
    }
}

/// One field's C1 answer, once both ladders exist.
struct Matched {
    /// `E*`, the error both arms are read off at.
    error: f64,
    /// Fitted triangle counts there, `[uniform, metric]`.
    triangles: [f64; 2],
    /// `triangles[metric] / triangles[uniform]`. C1 wanted [`WIN_RATIO`].
    ratio: f64,
    /// Whether `E*` lies inside the metric arm's own measured error range. False
    /// means the metric arm's count is an **extrapolation**, and a reader is owed
    /// that word.
    interpolated: bool,
}

/// Read both arms off at one matched error, or say why not.
///
/// `E*` is the **uniform arm's finest error** — a concrete measured quantity,
/// and the natural question: "at the error uniform refinement buys at the top of
/// the ladder, what does each arm cost?" Both arms are read from their fits
/// rather than one from a fit and one from a raw point, so neither gets a
/// smoothing the other did not.
fn match_arms(rungs: &[Rung]) -> Result<Matched, &'static str> {
    let mut uniform: Vec<(f64, f64)> = Vec::with_capacity(rungs.len());
    let mut metric: Vec<(f64, f64)> = Vec::with_capacity(rungs.len());
    for rung in rungs {
        let pair = rung.hausdorff?;
        let counts = rung.triangles;
        if counts[0] == 0 || counts[1] == 0 {
            return Err("unmeasurable:empty_mesh");
        }
        if pair[0] <= 0.0 || pair[1] <= 0.0 {
            return Err("unmeasurable:zero_error");
        }
        uniform.push((pair[0], counts[0] as f64));
        metric.push((pair[1], counts[1] as f64));
    }

    let span = |points: &[(f64, f64)]| {
        let mut lo = f64::INFINITY;
        let mut hi = 0.0f64;
        for &(error, _) in points {
            lo = lo.min(error);
            hi = hi.max(error);
        }
        (lo, hi)
    };
    let (uniform_lo, uniform_hi) = span(&uniform);
    let (metric_lo, metric_hi) = span(&metric);
    if uniform_hi < uniform_lo * LADDER_SPAN_FLOOR || metric_hi < metric_lo * LADDER_SPAN_FLOOR {
        return Err("unmeasurable:flat_ladder");
    }

    let error = uniform_lo;
    let triangles = [fit_log(&uniform).eval(error), fit_log(&metric).eval(error)];

    Ok(Matched {
        error,
        triangles,
        ratio: triangles[1] / triangles[0],
        interpolated: error >= metric_lo && error <= metric_hi,
    })
}

// ─── one rung ────────────────────────────────────────────────────────────────

/// Which arm a mesh belongs to.
#[derive(Component, Clone, Copy, Debug, PartialEq, Eq)]
enum Side {
    /// `(N, N, N)` samples, `h` equal on every axis.
    Uniform,
    /// `(n_x, n_y, n_z)` from the metric, at the same total budget.
    Metric,
}

impl Side {
    /// Index into the `[uniform, metric]` pairs every readout here is stored in.
    /// One source of truth for the column order.
    const fn index(self) -> usize {
        self as usize
    }
}

/// Everything one rung of the ladder produced.
struct Rung {
    /// Samples per axis on the uniform arm.
    samples: u32,
    /// The metric arm's per-axis sample counts.
    grid: [u32; 3],
    /// Axes pinned at [`MIN_SAMPLES`] while the budget was re-solved.
    pinned: usize,
    /// `∏ n_a / N³`. The two arms are budget-matched only if this is near one.
    budget_ratio: f64,
    /// `max n_a / min n_a`. One means the metric asked for the uniform grid.
    axis_ratio: f64,
    /// Triangles, `[uniform, metric]`.
    triangles: [usize; 2],
    /// Vertices, `[uniform, metric]`.
    vertices: [usize; 2],
    /// Symmetric Hausdorff `[uniform, metric]`, or why there is none.
    hausdorff: Result<[f64; 2], &'static str>,
    /// What the metric field looked like here.
    census: Census,
}

impl Rung {
    /// Whether the metric prescribed the grid the uniform arm already used.
    fn arms_identical(&self) -> bool {
        self.grid == [self.samples; 3]
    }
}

/// Median, min and max of one stage's repeats, in milliseconds.
///
/// [`REPEATS`] is odd, so the median is an observation rather than an average of
/// two.
#[derive(Clone, Copy)]
struct Timing {
    /// The headline.
    median: f64,
    /// Fastest repeat.
    min: f64,
    /// Slowest repeat.
    max: f64,
}

/// Median, min and max of `ms`, which must be non-empty.
///
/// Sorted by [`f64::total_cmp`]: there is no `partial_cmp` and no `unwrap` in
/// this file, and a NaN repeat would otherwise choose its own position.
fn timing(mut ms: Vec<f64>) -> Timing {
    ms.sort_by(f64::total_cmp);
    let last = ms.len() - 1;
    Timing {
        median: ms[last / 2],
        min: ms[0],
        max: ms[last],
    }
}

/// What C3 charges the metric, at the rung on screen.
///
/// Measured only here and not per rung, because it is the only rung whose cost
/// is reported and because five repeats of every rung would put a fifth of a
/// second between a keypress and a picture on the heaviest field.
struct Cost {
    /// `hessian` + `metric_lp` over the band points, and nothing else. Band
    /// *selection* is outside the clock: the extractor already visits every grid
    /// sample and the band test is one comparison on a value it already has, so
    /// charging the metric for that scan would charge it twice.
    metric: Timing,
    /// Each extraction, `[uniform, metric]`.
    extract: [Timing; 2],
    /// `HESSIAN_SAMPLES · band points`: field samples the metric read.
    metric_samples: u64,
    /// `N³`: samples in the uniform grid the extraction classifies.
    grid_samples: u64,
}

impl Cost {
    /// `metric.median / extract[uniform].median`, which is p-146's
    /// `metric_share`.
    fn share(&self) -> f64 {
        self.metric.median / self.extract[Side::Uniform.index()].median
    }

    /// The widest honest interval the repeats admit: `[min/max, max/min]`, which
    /// is p-146's `metric_share_lo` and `metric_share_hi`.
    fn share_span(&self) -> [f64; 2] {
        let extract = self.extract[Side::Uniform.index()];
        [self.metric.min / extract.max, self.metric.max / extract.min]
    }

    /// Whether that interval falls wholly on one side of C3's bar.
    ///
    /// A share that straddles the bar has not decided the clause on this run,
    /// and saying so is p-146's `c3_row_decisive` column.
    fn decisive(&self) -> bool {
        let [lo, hi] = self.share_span();
        (lo > C3_MAX_SHARE) == (hi > C3_MAX_SHARE)
    }

    /// `metric_samples / grid_samples`: the count form of the same statement,
    /// and a floor, because it does not charge the eigendecompositions.
    fn sample_share(&self) -> f64 {
        self.metric_samples as f64 / self.grid_samples as f64
    }
}

// ─── the whole measurement ───────────────────────────────────────────────────

/// One field's ladder, its verdicts, and the geometry of the shown rung.
struct Measured {
    /// `(field, rung)` this was built for.
    key: (usize, usize),
    /// Which rung of [`LADDER`] is on screen.
    shown: usize,
    /// The committed row this field is cited against.
    cited: &'static Cited,
    /// Every rung, in [`LADDER`] order. Never empty.
    rungs: Vec<Rung>,
    /// The matched-error read-off, or why there is none.
    matched: Result<Matched, &'static str>,
    /// What the shown rung cost, and what C3 is read from.
    cost: Cost,
    /// Wall time of the whole ladder.
    ladder_ms: f64,
    /// The field's domain, unshifted.
    lo: Vec3,
    /// The far corner of the same domain.
    hi: Vec3,
    /// Cell spacing of the shown rung, `[uniform, metric]`.
    spacing: [Vec3; 2],
    /// World `x` each arm is drawn at, `[uniform, metric]`.
    shift: [f32; 2],
}

impl Measured {
    /// The rung whose meshes are on screen.
    ///
    /// [`measure`] builds exactly one rung per [`LADDER`] entry and `shown` is
    /// clamped into that range before it is stored, so the index is proven by
    /// construction; the `min` is what makes that provable at a glance.
    fn shown_rung(&self) -> &Rung {
        &self.rungs[self.shown.min(self.rungs.len() - 1)]
    }
}

/// Why a field's Hausdorff cannot be measured, or `None` when it can.
///
/// `None` rather than a `"measurable"` sentinel, so the caller has one branch
/// and no unreachable arm: `validate::accuracy` compares against `|f|`, and
/// `fields/mod.rs:83-84` says only `Exact` makes that a distance.
fn unmeasurable_reason(bound: FieldBound) -> Option<&'static str> {
    match bound {
        FieldBound::Exact => None,
        FieldBound::Lipschitz { .. } => Some("unmeasurable:bound=Lipschitz"),
        FieldBound::Underestimate { .. } => Some("unmeasurable:bound=Underestimate"),
        FieldBound::Unbounded => Some("unmeasurable:bound=Unbounded"),
    }
}

/// A `MeshBuffer<f64>` as a Bevy mesh, narrowed to `f32` at the boundary.
///
/// The one narrowing in this file, and it happens after every measurement is
/// taken — `Real::as_f32` is lossy for `f64` and grading through it would
/// quantise the mesh before measuring it, which is the trap
/// `precision_f32_vs_f64.rs:351-358` records.
fn bevy_mesh(buffer: &MeshBuffer<f64>) -> Mesh {
    let positions: Vec<[f32; 3]> = buffer
        .positions
        .iter()
        .map(|p| [p[0].as_f32(), p[1].as_f32(), p[2].as_f32()])
        .collect();
    let normals: Vec<[f32; 3]> = buffer
        .normals
        .iter()
        .map(|n| [n[0].as_f32(), n[1].as_f32(), n[2].as_f32()])
        .collect();

    let mut mesh = Mesh::new(
        PrimitiveTopology::TriangleList,
        RenderAssetUsages::default(),
    );
    mesh.insert_attribute(Mesh::ATTRIBUTE_POSITION, positions);
    mesh.insert_attribute(Mesh::ATTRIBUTE_NORMAL, normals);
    mesh.insert_indices(Indices::U32(buffer.indices.clone()));
    mesh
}

/// Measure one reference field across the whole ladder, and build the shown
/// rung's two meshes.
///
/// Both arms are extracted by the same `MarchingCubes::new()` at its shipped
/// defaults and graded by **one instrument**: `validate::accuracy` with the
/// *uniform* rung's `shape`, `origin` and `cell_size` for the seed lattice in
/// both calls.
fn measure<F>(
    field: &F,
    cited: &'static Cited,
    shown: usize,
) -> Result<(Measured, [Mesh; 2]), String>
where
    F: ReferenceField + Sdf<Scalar = f64>,
{
    assert_eq!(
        F::NAME,
        cited.name,
        "E-319's field dispatch and its citation table disagree, so the HUD would quote one \
         field's committed row beside another field's mesh"
    );

    let started = Instant::now();
    let bound = field.bound();
    let unmeasurable = unmeasurable_reason(bound);
    let (lo, hi) = field.domain();
    let extent = hi[0] - lo[0];

    let mut mc = MarchingCubes::<f64>::new();
    let mut uniform = MeshBuffer::<f64>::new();
    let mut metric = MeshBuffer::<f64>::new();
    let mut metrics: Vec<Sym3> = Vec::new();
    let mut rungs: Vec<Rung> = Vec::with_capacity(LADDER.len());
    let mut shown_meshes: Option<[Mesh; 2]> = None;
    let mut shown_cost: Option<Cost> = None;
    let mut spacing = [Vec3::ZERO; 2];

    for (index, samples) in LADDER.into_iter().enumerate() {
        let is_shown = index == shown;
        let h = extent / f64::from(samples - 1);
        let shape = RuntimeShape3::new([samples; 3]).map_err(|error| {
            format!(
                "{}: the {samples}^3 uniform grid was rejected: {error}",
                cited.name
            )
        })?;

        // ── the metric field ────────────────────────────────────────────────
        //
        // Built once for the census, and -- at the rung on screen -- built
        // `REPEATS` more times against a clock, that first build serving as the
        // warm-up. What is timed is `hessian` + `metric_lp` and nothing else;
        // see [`Cost`] for what is deliberately left outside the clock.
        let points = band_points(field, lo, h, samples);
        if points.is_empty() {
            return Err(format!(
                "{} at {samples}^3 put no grid sample within {BAND_CELLS} cell of its surface, \
                 so every aspect ratio and every weight below would be a statistic of an empty \
                 population (M-44)",
                cited.name
            ));
        }
        metrics.reserve(points.len());
        metrics.clear();
        for &p in &points {
            metrics.push(metric_lp(&hessian(field, p, h)));
        }
        let mut metric_timing = None;
        if is_shown {
            let mut repeats = Vec::with_capacity(REPEATS);
            for _ in 0..REPEATS {
                let clock = Instant::now();
                metrics.clear();
                for &p in &points {
                    metrics.push(metric_lp(&hessian(field, p, h)));
                }
                repeats.push(clock.elapsed().as_secs_f64() * 1000.0);
            }
            metric_timing = Some(timing(repeats));
        }
        let census = census_of(field, &points, &metrics, h);

        // ── the uniform arm ─────────────────────────────────────────────────
        uniform.reset();
        mc.extract(field, &shape, lo, h, &mut uniform)
            .map_err(|error| format!("{} at {samples}^3, uniform arm: {error}", cited.name))?;
        let mut uniform_timing = None;
        if is_shown {
            let mut repeats = Vec::with_capacity(REPEATS);
            for _ in 0..REPEATS {
                uniform.reset();
                let clock = Instant::now();
                mc.extract(field, &shape, lo, h, &mut uniform)
                    .map_err(|error| {
                        format!("{} at {samples}^3, uniform arm: {error}", cited.name)
                    })?;
                repeats.push(clock.elapsed().as_secs_f64() * 1000.0);
            }
            uniform_timing = Some(timing(repeats));
        }

        // ── the metric arm, at the same total sample budget ──────────────────
        let (grid, pinned) = metric_grid(census.weights, samples);
        let stretch = [
            extent / f64::from(grid[0] - 1) / h,
            extent / f64::from(grid[1] - 1) / h,
            extent / f64::from(grid[2] - 1) / h,
        ];
        let stretched = Stretched {
            field,
            lo,
            s: stretch,
        };
        let metric_shape = RuntimeShape3::new(grid).map_err(|error| {
            format!(
                "{}: the metric grid {}x{}x{} was rejected: {error}",
                cited.name, grid[0], grid[1], grid[2]
            )
        })?;
        // Timed before the positions are mapped back to world space, because a
        // repeat re-extracts in `q` and the map-back is this file's own work
        // rather than the extractor's.
        metric.reset();
        mc.extract(&stretched, &metric_shape, [0.0; 3], h, &mut metric)
            .map_err(|error| format!("{} at {samples}^3, metric arm: {error}", cited.name))?;
        let mut metric_arm_timing = None;
        if is_shown {
            let mut repeats = Vec::with_capacity(REPEATS);
            for _ in 0..REPEATS {
                metric.reset();
                let clock = Instant::now();
                mc.extract(&stretched, &metric_shape, [0.0; 3], h, &mut metric)
                    .map_err(|error| {
                        format!("{} at {samples}^3, metric arm: {error}", cited.name)
                    })?;
                repeats.push(clock.elapsed().as_secs_f64() * 1000.0);
            }
            metric_arm_timing = Some(timing(repeats));
        }

        for p in &mut metric.positions {
            p[0] = p[0].mul_add(stretch[0], lo[0]);
            p[1] = p[1].mul_add(stretch[1], lo[1]);
            p[2] = p[2].mul_add(stretch[2], lo[2]);
        }
        // A normal is a covector, so the stretch acts on it by the inverse
        // transpose -- component-wise *division* by `s`, not multiplication.
        // p-146 leaves this alone because nothing downstream of a bench reads a
        // normal; here the mesh is shaded, and multiplying would light the
        // metric arm as though it were a different surface.
        for n in &mut metric.normals {
            let raw = [n[0] / stretch[0], n[1] / stretch[1], n[2] / stretch[2]];
            let length = (raw[0] * raw[0] + raw[1] * raw[1] + raw[2] * raw[2]).sqrt();
            if length > 0.0 {
                *n = [raw[0] / length, raw[1] / length, raw[2] / length];
            }
        }

        // ── one instrument, both arms ───────────────────────────────────────
        let hausdorff = match unmeasurable {
            Some(reason) => Err(reason),
            None => {
                let cfg = AccuracyConfig::from_cell_size(h).map_err(|error| {
                    format!(
                        "{} at {samples}^3: cell size {h} rejected: {error}",
                        cited.name
                    )
                })?;
                let uniform_report = accuracy(
                    &uniform.positions,
                    &uniform.indices,
                    field,
                    &shape,
                    lo,
                    &cfg,
                )
                .map_err(|error| {
                    format!("{} at {samples}^3, uniform accuracy: {error}", cited.name)
                })?;
                let metric_report =
                    accuracy(&metric.positions, &metric.indices, field, &shape, lo, &cfg).map_err(
                        |error| format!("{} at {samples}^3, metric accuracy: {error}", cited.name),
                    )?;
                if uniform_report.has_coverage() && metric_report.has_coverage() {
                    Ok([
                        uniform_report.symmetric_hausdorff(),
                        metric_report.symmetric_hausdorff(),
                    ])
                } else {
                    Err("unmeasurable:no_coverage")
                }
            }
        };

        if let (Some(metric_stage), Some(uniform_stage), Some(metric_arm)) =
            (metric_timing, uniform_timing, metric_arm_timing)
        {
            shown_meshes = Some([bevy_mesh(&uniform), bevy_mesh(&metric)]);
            shown_cost = Some(Cost {
                metric: metric_stage,
                extract: [uniform_stage, metric_arm],
                metric_samples: HESSIAN_SAMPLES * points.len() as u64,
                grid_samples: u64::from(samples).pow(3),
            });
            spacing = [
                Vec3::splat(h.as_f32()),
                Vec3::new(
                    (h * stretch[0]).as_f32(),
                    (h * stretch[1]).as_f32(),
                    (h * stretch[2]).as_f32(),
                ),
            ];
        }

        let product = f64::from(grid[0]) * f64::from(grid[1]) * f64::from(grid[2]);
        let axis_hi = f64::from(grid.iter().copied().max().unwrap_or(samples));
        let axis_lo = f64::from(grid.iter().copied().min().unwrap_or(samples));
        rungs.push(Rung {
            samples,
            grid,
            pinned,
            budget_ratio: product / f64::from(samples).powi(3),
            axis_ratio: axis_hi / axis_lo,
            triangles: [uniform.triangle_count(), metric.triangle_count()],
            vertices: [uniform.positions.len(), metric.positions.len()],
            hausdorff,
            census,
        });
    }

    // `shown` is clamped into the ladder by both callers, so exactly one pass of
    // the loop above fills these. A missing pair is a bug in that clamp rather
    // than a state a reader could reach, and it says so instead of unwrapping.
    let (Some(meshes), Some(cost)) = (shown_meshes, shown_cost) else {
        return Err(format!(
            "{}: rung {shown} is not on a ladder of {} rungs",
            cited.name,
            LADDER.len()
        ));
    };
    let matched = match_arms(&rungs);
    let width = extent.as_f32();
    let offset = width * 0.62;

    Ok((
        Measured {
            key: (0, shown),
            shown,
            cited,
            rungs,
            matched,
            cost,
            ladder_ms: started.elapsed().as_secs_f64() * 1000.0,
            lo: Vec3::new(lo[0].as_f32(), lo[1].as_f32(), lo[2].as_f32()),
            hi: Vec3::new(hi[0].as_f32(), hi[1].as_f32(), hi[2].as_f32()),
            spacing,
            shift: [-offset, offset],
        },
        meshes,
    ))
}

/// Dispatch on the field index.
///
/// One arm per entry of [`CITED`], in the same order, which is the order
/// `for_each_reference_field!` yields and the order the digit keys index. The
/// eight fields are eight different types, so this is a `match` rather than a
/// loop.
fn build(key: (usize, usize)) -> Result<(Measured, [Mesh; 2]), String> {
    let (field, shown) = key;
    let cited = &CITED[field];
    let measured = match field {
        1 => measure(&Torus::<f64>::canonical(), cited, shown),
        2 => measure(&BoxExact::<f64>::canonical(), cited, shown),
        3 => measure(&csg_difference::<f64>(), cited, shown),
        4 => measure(&ThinPlate::<f64>::canonical(), cited, shown),
        5 => measure(&capped_gyroid::<f64>(), cited, shown),
        6 => measure(&FbmTerrain::<f64>::canonical(), cited, shown),
        7 => measure(&noise_cavity::<f64>(), cited, shown),
        _ => measure(&Sphere::<f64>::canonical(), cited, shown),
    };
    measured.map(|(mut measured, meshes)| {
        measured.key = key;
        (measured, meshes)
    })
}

// ─── the app ─────────────────────────────────────────────────────────────────

/// What the reader has asked to see.
#[derive(Resource)]
struct Demo {
    /// Index into [`CITED`].
    field: usize,
    /// Index into [`LADDER`].
    shown: usize,
    /// `ISOMESH_FIELD`, re-parsed here with a range check.
    ///
    /// The harness parses that variable into `ViewFlags::field` already, but its
    /// default is `0` and so is indistinguishable from an explicit
    /// `ISOMESH_FIELD=0` — and this demo needs the distinction, because a pin has
    /// to beat the capture walk. Several examples re-parse it for the same reason
    /// (`aperture_gate.rs:1026-1034`).
    pinned_field: Option<usize>,
    /// `ISOMESH_SAMPLES`, as an index into [`LADDER`].
    pinned_rung: Option<usize>,
}

/// One material per arm, so the two sides are nameable in the caption.
#[derive(Resource)]
struct Materials {
    /// The uniform arm's, grey.
    uniform: Handle<StandardMaterial>,
    /// The metric arm's, tan.
    metric: Handle<StandardMaterial>,
}

/// What the last rebuild produced.
///
/// Three named states rather than an `Option` beside an error field: a HUD
/// showing the previous field's numbers under this field's caption is the exact
/// failure `game_mirror_dedup.rs:1350-1353` records, and a rejected rebuild is a
/// thing a reader must be told rather than a reason to keep stale numbers.
///
/// The measured variant is boxed because it is an order larger than a rejection
/// string and `clippy::large_enum_variant` is right about that.
#[derive(Resource, Default)]
enum Readout {
    /// Before the first rebuild.
    #[default]
    Pending,
    /// The last rebuild's numbers.
    Measured(Box<Measured>),
    /// The last rebuild was refused, with the reason.
    Rejected(String),
}

/// The grid cell each arm samples on. Its own group so it is not drawn at the
/// shared wireframe's width or depth bias.
#[derive(Default, Reflect, GizmoConfigGroup)]
struct CellGizmos;

fn main() {
    App::new()
        .add_plugins(DefaultPlugins.set(WindowPlugin {
            primary_window: Some(Window {
                title: "isomesh - E-319 anisotropic metric".into(),
                ..default()
            }),
            ..default()
        }))
        .add_plugins(CommonPlugin)
        .init_gizmo_group::<CellGizmos>()
        .init_resource::<Readout>()
        .insert_resource(Demo {
            field: 0,
            shown: DEFAULT_RUNG,
            pinned_field: pinned_field(),
            pinned_rung: pinned_rung(),
        })
        .add_systems(Startup, setup)
        // `PreUpdate` for E-306's reason, restated at
        // `game_mirror_dedup.rs:1350-1354`: the harness's `update_hud` and
        // `capture_sequence` both run in `Update` with no ordering against an
        // example's own systems, so a `report` there would render a frame-old
        // readout beside a current caption. Chained, because `report` formats
        // exactly what `remesh` just measured.
        .add_systems(PreUpdate, (controls, remesh, report).chain())
        .add_systems(Update, draw_cells)
        .run();
}

/// `ISOMESH_FIELD` as an index into [`CITED`], or `None`.
///
/// Out of range logs once and pins nothing, which is what the other examples
/// that re-parse this variable do.
fn pinned_field() -> Option<usize> {
    let raw = std::env::var("ISOMESH_FIELD").ok()?;
    match raw.parse::<usize>() {
        Ok(index) if index < CITED.len() => Some(index),
        _ => {
            error!(
                "ISOMESH_FIELD={raw} is not one of 0..{}; the field list is {:?}",
                CITED.len() - 1,
                CITED.map(|cited| cited.name)
            );
            None
        }
    }
}

/// `ISOMESH_SAMPLES` as an index into [`LADDER`], or `None`.
///
/// The harness's contract is that a pinned resolution wins over a sweep
/// (`active_cells.rs:1137-1141`), so this is honoured ahead of the capture walk
/// and the keys. It must name a rung: the ladder is what the read-off is fitted
/// through, and a fourth resolution outside it would be a grid with no fit.
fn pinned_rung() -> Option<usize> {
    let samples = common::samples_override()?;
    match LADDER.iter().position(|&rung| rung == samples) {
        Some(index) => Some(index),
        None => {
            error!("ISOMESH_SAMPLES={samples} is not a rung of this demo's ladder {LADDER:?}");
            None
        }
    }
}

fn setup(
    mut commands: Commands,
    mut materials: ResMut<Assets<StandardMaterial>>,
    mut camera: Query<&mut OrbitCamera>,
    mut gizmo_config: ResMut<GizmoConfigStore>,
) {
    for mut orbit in &mut camera {
        // Look straight down -z at the pair, so neither side is nearer than the
        // other and perspective does not make one look bigger than it is. A
        // side-by-side comparison viewed from an angle is not a comparison
        // (E-103). `radius` and `focus` are set per rebuild instead, because
        // both scale with the field's domain and the gyroid's is three and a half
        // times the others'.
        orbit.yaw = FRAC_PI_2;
        orbit.pitch = 0.10;
    }

    // Biased forward and drawn thin: a grid cell is the smallest thing on
    // screen, it sits at a corner *inside* the surface, and an exact cage that
    // hides behind the mesh it belongs to shows nothing.
    let (cells, _) = gizmo_config.config_mut::<CellGizmos>();
    cells.line.width = 2.4;
    cells.depth_bias = -0.8;

    commands.insert_resource(Materials {
        uniform: materials.add(StandardMaterial {
            base_color: UNIFORM_SURFACE,
            perceptual_roughness: 0.45,
            ..default()
        }),
        metric: materials.add(StandardMaterial {
            base_color: METRIC_SURFACE,
            perceptual_roughness: 0.45,
            ..default()
        }),
    });
}

/// Precedence: a pin beats the capture walk, which beats the keys.
fn controls(
    keys: Res<ButtonInput<KeyCode>>,
    flags: Res<ViewFlags>,
    capture: Res<Capture>,
    mut demo: ResMut<Demo>,
) {
    match demo.pinned_rung {
        Some(rung) => demo.shown = rung,
        None => {
            if keys.just_pressed(KeyCode::BracketRight) {
                demo.shown = (demo.shown + 1).min(LADDER.len() - 1);
            }
            if keys.just_pressed(KeyCode::BracketLeft) {
                demo.shown = demo.shown.saturating_sub(1);
            }
        }
    }

    if let Some(field) = demo.pinned_field {
        demo.field = field;
    } else if capture.is_active() {
        demo.field = walked_field(capture.taken);
    } else if flags.field < CITED.len() {
        demo.field = flags.field;
    }
}

/// The field the capture walk is on after `taken` captured frames.
///
/// Eight fields at [`CAPTURE_FRAMES_PER_FIELD`] frames each is eighty, which is
/// `record_gif.sh:47`'s default `ISOMESH_CAPTURE_FRAMES`, so the clip is exactly
/// one pass over the roster and ends where it started.
///
/// Split out of [`controls`] so it can be tested. `Capture::is_active()` is
/// false unless `ISOMESH_CAPTURE` names a directory, and `unsafe_code =
/// "forbid"` means a test cannot set an environment variable -- the same reason
/// `ViewFlags::parse` was split out of `Default` (`common/mod.rs:122-126`). The
/// arithmetic is worth a test on its own account: `M-241` records a GIF that
/// advertised a sweep it never performed, and a single still passed the check
/// that was supposed to catch it.
fn walked_field(taken: u32) -> usize {
    (taken / CAPTURE_FRAMES_PER_FIELD) as usize % CITED.len()
}

#[allow(clippy::too_many_arguments)]
fn remesh(
    demo: Res<Demo>,
    mut readout: ResMut<Readout>,
    mut assets: ResMut<Assets<Mesh>>,
    materials: Res<Materials>,
    mut commands: Commands,
    mut query: Query<(&mut Mesh3d, &mut Transform, &mut DemoDomain, &Side)>,
    mut camera: Query<&mut OrbitCamera>,
    mut last: Local<Option<(usize, usize)>>,
) {
    let key = (
        demo.field.min(CITED.len() - 1),
        demo.shown.min(LADDER.len() - 1),
    );
    if *last == Some(key) {
        return;
    }
    *last = Some(key);

    let (measured, meshes) = match build(key) {
        Ok(pair) => pair,
        Err(reason) => {
            error!("E-319 rebuild refused: {reason}");
            *readout = Readout::Rejected(reason);
            return;
        }
    };
    let [uniform_mesh, metric_mesh] = meshes;

    // The HUD is the evidence and a headless capture has no HUD to read. One
    // line per rebuild, so `ISOMESH_CAPTURE` leaves the read-off in the log
    // where a script can hold it against p-146 (E-203 learned this the hard
    // way: a measurement that only exists on screen cannot be verified from a
    // terminal).
    let rung = measured.shown_rung();
    info!(
        "E-319 {} at {}^3: uniform {}x{}x{} -> {} tri, metric {}x{}x{} -> {} tri; \
         matched ratio {} (p-146 ratio {}); metric/extract {:.4} in [{:.4}, {:.4}] \
         (p-146 {}), and {} metric samples against the grid's {} = {:.4}; \
         aspect max {:.3e} with {} of {} band cells at the H_FLOOR = {:e}; \
         flat axis {}; ladder {:.1} ms",
        measured.cited.name,
        rung.samples,
        rung.samples,
        rung.samples,
        rung.samples,
        rung.triangles[Side::Uniform.index()],
        rung.grid[0],
        rung.grid[1],
        rung.grid[2],
        rung.triangles[Side::Metric.index()],
        measured.matched.as_ref().map_or_else(
            |reason| (*reason).to_string(),
            |m| format!("{:.4}", m.ratio)
        ),
        measured.cited.p146_ratio,
        measured.cost.share(),
        measured.cost.share_span()[0],
        measured.cost.share_span()[1],
        measured.cited.p146_share,
        measured.cost.metric_samples,
        measured.cost.grid_samples,
        measured.cost.sample_share(),
        rung.census.aspect[0],
        rung.census.at_floor,
        rung.census.points,
        H_FLOOR,
        flat_axis(&rung.census),
        measured.ladder_ms,
    );

    let width = measured.hi.x - measured.lo.x;
    for mut orbit in &mut camera {
        orbit.radius = width * 2.75;
        // Raised, so the pair sits low and the panel -- which is two dozen lines
        // here -- is not drawn across the evidence. E-112's lesson and E-109's
        // committed screenshot, applied to a pair rather than to one subject.
        orbit.focus = Vec3::Y * (width * 0.30);
    }

    // `Mesh3d::default()` names no asset and draws nothing, which is what an
    // empty result actually wants: an empty mesh in `Assets` produces
    // `Use-after-free: attempted to copy element data for an unallocated key`
    // twice a frame, forever (`active_cells.rs:1240-1247`).
    let mut handle = |mesh: Mesh, triangles: usize| {
        if triangles == 0 {
            Handle::default()
        } else {
            assets.add(mesh)
        }
    };
    let built = [
        (
            Side::Uniform,
            handle(uniform_mesh, rung.triangles[Side::Uniform.index()]),
            measured.shift[Side::Uniform.index()],
        ),
        (
            Side::Metric,
            handle(metric_mesh, rung.triangles[Side::Metric.index()]),
            measured.shift[Side::Metric.index()],
        ),
    ];

    // Spawn once and swap the handle thereafter, which is what keeps the two
    // `Side` markers -- and so the two materials -- attached to the same halves
    // of the screen across a rebuild (E-103).
    if query.is_empty() {
        for (side, mesh, x) in &built {
            let material = match side {
                Side::Uniform => materials.uniform.clone(),
                Side::Metric => materials.metric.clone(),
            };
            commands.spawn((
                Mesh3d(mesh.clone()),
                MeshMaterial3d(material),
                Transform::from_xyz(*x, 0.0, 0.0),
                DemoDomain {
                    min: measured.lo + Vec3::X * *x,
                    max: measured.hi + Vec3::X * *x,
                },
                DemoMesh,
                *side,
            ));
        }
    } else {
        for (mut mesh, mut transform, mut domain, side) in &mut query {
            if let Some((_, handle, x)) = built.iter().find(|(s, _, _)| s == side) {
                mesh.0 = handle.clone();
                transform.translation.x = *x;
                domain.min = measured.lo + Vec3::X * *x;
                domain.max = measured.hi + Vec3::X * *x;
            }
        }
    }

    *readout = Readout::Measured(Box::new(measured));
}

fn report(readout: Res<Readout>, mut stats: ResMut<DemoStats>) {
    let measured = match &*readout {
        Readout::Pending => return,
        Readout::Rejected(reason) => {
            stats.title = "E-319  anisotropic metric".into();
            stats.banner = Some((format!("REBUILD REFUSED: {reason}"), FALSIFIED));
            stats.extra = vec![format!("no numbers: {reason}")];
            // Cleared, not left: the harness prints these three unconditionally,
            // and the previous field's counts under a REBUILD REFUSED headline
            // would be the exact staleness this enum exists to prevent.
            stats.vertices = 0;
            stats.triangles = 0;
            stats.extract_ms = 0.0;
            return;
        }
        Readout::Measured(measured) => measured,
    };

    let rung = measured.shown_rung();
    let uniform = Side::Uniform.index();
    let metric = Side::Metric.index();
    let ladder: Vec<String> = LADDER.iter().map(u32::to_string).collect();
    let census = &rung.census;

    stats.title = format!(
        "E-319  uniform grid (grey, left)  vs  metric-driven grid (tan, right)\n\
         \x20      field: {}   [1-7] to switch",
        measured.cited.name
    );
    stats.vertices = rung.vertices[uniform] + rung.vertices[metric];
    stats.triangles = rung.triangles[uniform] + rung.triangles[metric];
    stats.extract_ms = measured.cost.extract[uniform].median + measured.cost.extract[metric].median;

    let mut lines = vec![
        format!(
            "ladder {} samples/axis   showing {}^3   [ and ] to move",
            ladder.join(" / "),
            rung.samples
        ),
        String::new(),
        format!("{:<14}{:>14}  {:>14}", "", "uniform", "metric"),
        format!(
            "{:<14}{:>14}  {:>14}",
            "grid",
            format!("{0}x{0}x{0}", rung.samples),
            format!("{}x{}x{}", rung.grid[0], rung.grid[1], rung.grid[2])
        ),
        format!(
            "{:<14}{:>14}  {:>14}",
            "triangles", rung.triangles[uniform], rung.triangles[metric]
        ),
    ];

    lines.push(match rung.hausdorff {
        Ok(pair) => format!(
            "{:<14}{:>14.4e}  {:>14.4e}",
            "sym.hausdorff", pair[uniform], pair[metric]
        ),
        Err(reason) => format!("{:<14}{reason}", "sym.hausdorff"),
    });
    lines.push(format!(
        "{:<14}{:>14.3}  {:>14.3}",
        "extract ms", measured.cost.extract[uniform].median, measured.cost.extract[metric].median
    ));
    lines.push(format!(
        "budget {:.3} x N^3   axis ratio {:.3}   axes pinned at {} samples: {}",
        rung.budget_ratio, rung.axis_ratio, MIN_SAMPLES, rung.pinned
    ));
    lines.push(String::new());

    // ── C3, live: a clock reported with its scatter, and a count beside it ───
    let cost = &measured.cost;
    let share = cost.share();
    let [share_lo, share_hi] = cost.share_span();
    lines.push(format!(
        "metric build {:>8.3} ms = {:>7.3} x extract   bar <= {:.2}   {} (P-146 C3)",
        cost.metric.median,
        share,
        C3_MAX_SHARE,
        if share <= C3_MAX_SHARE {
            "HOLDS"
        } else {
            "FALSIFIED"
        }
    ));
    lines.push(format!(
        "  median of {REPEATS} after a warm-up; the repeats span [{share_lo:.3}, {share_hi:.3}], \
         {}",
        if cost.decisive() {
            "one side of the bar"
        } else {
            "STRADDLING the bar: this run did not decide C3"
        }
    ));
    lines.push(format!(
        "  and as a count: {} x {} band = {} field samples vs the grid's {} = {:.3} x",
        HESSIAN_SAMPLES,
        census.points,
        cost.metric_samples,
        cost.grid_samples,
        cost.sample_share()
    ));
    lines.push(
        "  a floor, since it does not charge the one eigendecomposition per band cell".into(),
    );
    lines.push(String::new());

    // ── C1, live ────────────────────────────────────────────────────────────
    match &measured.matched {
        Ok(matched) => {
            lines.push(format!(
                "matched-Hausdorff read-off: E* = {:.4e}, the uniform arm's finest",
                matched.error
            ));
            lines.push(format!(
                "  T_uniform {:>10.1}   T_metric {:>10.1}   ratio {:.4}   bar <= {:.3}",
                matched.triangles[uniform], matched.triangles[metric], matched.ratio, WIN_RATIO
            ));
            lines.push(format!(
                "  {} (P-146 C1)",
                if matched.ratio <= WIN_RATIO {
                    "HOLDS: the metric arm saved at least a quarter of the triangles"
                } else {
                    "FALSIFIED: the metric arm needs MORE triangles at the same error"
                }
            ));
            lines.push(format!(
                "  read off by {}",
                if matched.interpolated {
                    "interpolation, E* inside the metric arm's own error range"
                } else {
                    "EXTRAPOLATION beyond the metric arm's own error range"
                }
            ));
        }
        Err(reason) => {
            lines.push(format!("matched-Hausdorff read-off: {reason}"));
            lines.push(
                "  validate::accuracy compares against |f|, and only a bound() of Exact".into(),
            );
            lines.push(
                "  makes that a distance (fields/mod.rs:83-84). 4 of 8 fields qualify.".into(),
            );
        }
    }
    lines.push(String::new());

    // ── the aspect ratio, never without its at-floor count ──────────────────
    lines.push(format!(
        "aspect ratio  max {:>10.3e}  mean {:>10.3e}   at floor {} of {} band cells",
        census.aspect[0], census.aspect[1], census.at_floor, census.points
    ));
    lines.push(if census.at_floor == census.points {
        format!(
            "  off the floor: NO CELL IS - all {} of them restate H_FLOOR = {:e}",
            census.points, H_FLOOR
        )
    } else {
        format!(
            "  off the floor max {:>10.3e}  mean {:>10.3e}   over the other {} <- measured",
            census.aspect_off_floor[0],
            census.aspect_off_floor[1],
            census.points - census.at_floor
        )
    });
    lines.push(format!(
        "  H_FLOOR = {:e}, so an at-floor ratio is |lambda|max / H_FLOOR: the floor talking",
        H_FLOOR
    ));
    lines.push(format!(
        "flat axis {}   axis-aligned fraction {:.4}   {}",
        flat_axis(census),
        census.exploitable_flat as f64 / census.points as f64,
        if rung.arms_identical() {
            "IDENTICAL GRIDS: no anisotropy prescribed"
        } else {
            "the metric asked for a different grid"
        }
    ));
    lines.push(String::new());

    // ── the citations ───────────────────────────────────────────────────────
    lines.push(format!(
        "p-146.csv {:<15} ratio {}  metric_share {} at 33^3",
        measured.cited.name, measured.cited.p146_ratio, measured.cited.p146_share
    ));
    lines.push(P146_GLOBAL.into());
    lines.push(format!(
        "p-147.csv {:<15} exponent_difference {}  constant_ratio {}",
        measured.cited.name, measured.cited.p147_exponent, measured.cited.p147_constant
    ));
    lines.push(P147_GLOBAL.into());
    lines.push(P149_GLOBAL.into());

    stats.banner = Some(banner(measured, rung));
    stats.extra = lines;
    stats.keys = Some(
        "[1-7] field (ISOMESH_FIELD=7 for the eighth)   [ and ] rung   [W] wire   [G] boxes\n\
         [Space] pause   [R] re-mesh   [H] HUD   [F12] shot   [Esc] quit"
            .into(),
    );
}

/// The one line above the panel: what this field's arms did, and in what sense.
///
/// Identical grids take precedence over a ratio, because "the metric prescribed
/// the grid we already had" is a stronger and more surprising statement than
/// "the ratio is one", and on five of the eight reference fields it is the whole
/// result (P-149's 70 of 112 rows).
fn banner(measured: &Measured, rung: &Rung) -> (String, Color) {
    if rung.arms_identical() {
        return (
            format!(
                "IDENTICAL GRIDS: the metric prescribed {0}x{0}x{0}, the grid the uniform arm \
                 already used",
                rung.samples
            ),
            UNMEASURABLE,
        );
    }
    match &measured.matched {
        Ok(matched) if matched.ratio <= WIN_RATIO => (
            format!(
                "SAVED {:.1}% of the triangles at matched Hausdorff",
                (1.0 - matched.ratio) * 100.0
            ),
            HELD,
        ),
        Ok(matched) => (
            format!(
                "NO SAVING: {:.3}x the triangles at matched Hausdorff, against C1's bar of {:.2}",
                matched.ratio, WIN_RATIO
            ),
            FALSIFIED,
        ),
        Err(reason) => (format!("NO HAUSDORFF TO MATCH ON: {reason}"), UNMEASURABLE),
    }
}

/// Draw each arm's grid cell at its low domain corner.
///
/// At its **exact** size, not inflated: a cell drawn larger than it is would
/// make every metric look as anisotropic as the reader expected, which is the
/// one thing this picture must not do (E-304's cage earned that rule). The
/// sample *counts* go in the panel rather than on screen as tick marks for the
/// same reason -- ticks at 191 samples over four units smear into a solid line,
/// and a smear is not a density.
fn draw_cells(readout: Res<Readout>, mut gizmos: Gizmos<CellGizmos>) {
    let Readout::Measured(measured) = &*readout else {
        return;
    };
    for (side, colour) in [(Side::Uniform, UNIFORM_CELL), (Side::Metric, METRIC_CELL)] {
        let index = side.index();
        let min = measured.lo + Vec3::X * measured.shift[index];
        cell(&mut gizmos, min, measured.spacing[index], colour);
    }
}

/// The twelve edges of one cell. Bit `i` of the corner index is axis `i`, the
/// same convention the extractor and `common::draw_domain` use.
fn cell(gizmos: &mut Gizmos<CellGizmos>, min: Vec3, size: Vec3, colour: Color) {
    let corner = |i: usize| {
        min + Vec3::new(
            if i & 1 == 0 { 0.0 } else { size.x },
            if i & 2 == 0 { 0.0 } else { size.y },
            if i & 4 == 0 { 0.0 } else { size.z },
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

#[cfg(test)]
mod tests {
    use bevy::ecs::system::RunSystemOnce;

    use super::*;

    /// `thin_plate`, index 4 of [`CITED`].
    ///
    /// The field the whole demo turns on: it is one of the four with a
    /// `bound()` of `Exact`, so it has a real Hausdorff and a real read-off; its
    /// flat direction is real enough to move the grid to `13x191x13`; and it is
    /// the field on which P-146's `ratio` is furthest from a win at
    /// `2.298935`. It is also the cheapest of the four to mesh, which matters
    /// because `cargo test` builds unoptimised and a debug extraction is 37-62x
    /// slower (`M-152`).
    const THIN_PLATE: usize = 4;

    /// The demo's own systems, in an `App` with no window and no renderer.
    ///
    /// This is the closest thing to running the demo that a machine with no
    /// display can do. `setup`, `draw_cells` and the harness's own systems are
    /// left out: the first wants `Assets<StandardMaterial>` and a camera, the
    /// second wants `Gizmos`, and `report` is run below as a one-shot with a
    /// `DemoStats` inserted, which is the same system the demo runs every frame.
    ///
    /// No `TimeUpdateStrategy`: nothing in this demo reads `Time`, so pinning a
    /// clock would pin one nothing consults.
    fn harness(field: usize) -> App {
        let mut app = App::new();
        app.add_plugins(MinimalPlugins)
            .add_plugins(bevy::asset::AssetPlugin::default())
            .init_asset::<Mesh>()
            .init_resource::<ButtonInput<KeyCode>>()
            .init_resource::<ViewFlags>()
            .init_resource::<Readout>()
            .insert_resource(Capture::default())
            .insert_resource(Materials {
                uniform: Handle::default(),
                metric: Handle::default(),
            })
            .insert_resource(Demo {
                field,
                shown: DEFAULT_RUNG,
                // Pinned, so the run is the same whatever `ISOMESH_FIELD` and
                // `ISOMESH_SAMPLES` happen to say in the environment this test
                // was started in. `unsafe_code = "forbid"` means a test cannot
                // set them, so it has to be able to out-rank them.
                pinned_field: Some(field),
                pinned_rung: Some(DEFAULT_RUNG),
            })
            .add_systems(Update, (controls, remesh).chain());
        app
    }

    /// One frame, with the input clearing `InputPlugin` would have done.
    ///
    /// Without it `just_pressed` stays true for ever, and a bracket key would
    /// move the rung on every frame of the test rather than on the one it was
    /// pressed.
    fn step(app: &mut App) {
        app.update();
        app.world_mut()
            .resource_mut::<ButtonInput<KeyCode>>()
            .clear();
    }

    /// Step until the readout is the one the current keys ask for.
    ///
    /// A stall detector rather than a flat frame count: the ladder is
    /// synchronous, so one frame is enough today, and a count that is exactly
    /// right is a count that says nothing when it stops being right. The panic
    /// names what was outstanding.
    fn settle(app: &mut App) {
        const CAP: usize = 8;
        for _ in 0..CAP {
            step(app);
            let key = {
                let demo = app.world().resource::<Demo>();
                (demo.field, demo.shown)
            };
            match app.world().resource::<Readout>() {
                Readout::Measured(measured) if measured.key == key => return,
                Readout::Rejected(reason) => {
                    panic!("E-319 refused the ladder for key {key:?}: {reason}")
                }
                _ => {}
            }
        }
        panic!("the ladder did not settle in {CAP} frames");
    }

    /// The uniform arm is on the left, the metric arm on the right, and a
    /// rebuild swaps handles rather than spawning a second pair.
    ///
    /// The caption says *"uniform grid (grey, left) vs metric-driven grid (tan,
    /// right)"*, and on a machine with no display that sentence is otherwise
    /// unfalsifiable — a sign flip in the offset would leave every number correct
    /// and every word of the caption wrong. The reuse half is the other failure
    /// this wiring has: E-103's spawn-once-then-swap exists because spawning per
    /// rebuild leaves the previous pair in the world, and eight field switches
    /// would end with sixteen overlapping meshes and a HUD describing two of
    /// them.
    #[test]
    fn the_uniform_arm_is_on_the_left_and_a_rebuild_reuses_its_entities() {
        let mut app = harness(THIN_PLATE);
        settle(&mut app);

        let placed = |app: &mut App| -> Vec<(Side, f32, AssetId<Mesh>, Vec3)> {
            let mut out: Vec<(Side, f32, AssetId<Mesh>, Vec3)> = app
                .world_mut()
                .query::<(&Side, &Transform, &Mesh3d, &DemoDomain)>()
                .iter(app.world())
                .map(|(side, transform, mesh, domain)| {
                    (*side, transform.translation.x, mesh.0.id(), domain.min)
                })
                .collect();
            out.sort_by_key(|(side, ..)| side.index());
            out
        };

        let first = placed(&mut app);
        assert_eq!(
            first.len(),
            2,
            "the demo drew {} arms rather than a pair",
            first.len()
        );
        let (uniform_x, metric_x) = (first[0].1, first[1].1);
        assert!(
            uniform_x < 0.0 && metric_x > 0.0,
            "the uniform arm is at x {uniform_x} and the metric arm at x {metric_x}, so the \
             caption's 'uniform (left) vs metric-driven (right)' is the wrong way round"
        );
        assert!(
            (uniform_x + metric_x).abs() < f32::EPSILON,
            "the pair is not symmetric about the origin: {uniform_x} and {metric_x}, so one \
             side sits nearer the camera's centre than the other"
        );
        // Each arm's `DemoDomain` is the field's own box plus that arm's shift,
        // so `[G]` outlines the grid the mesh beside it was sampled on. The
        // unshifted corner is read off the uniform arm and held against both.
        let unshifted = first[0].3.x - first[0].1;
        for (side, x, _, domain_min) in &first {
            assert!(
                (domain_min.x - (unshifted + x)).abs() < 1e-4,
                "{side:?}'s domain box is at {} where its mesh is at {x}, so [G] would outline \
                 the wrong place",
                domain_min.x
            );
        }

        // A different field, driven the way the pin drives it, so the rebuild
        // path under test is the one the demo runs.
        app.world_mut().resource_mut::<Demo>().pinned_field = Some(0);
        settle(&mut app);
        let second = placed(&mut app);
        assert_eq!(
            second.len(),
            2,
            "a rebuild left {} arms in the world instead of reusing the pair",
            second.len()
        );
        for (before, after) in first.iter().zip(&second) {
            assert_eq!(before.0, after.0, "an arm changed sides across a rebuild");
            assert_ne!(
                before.2, after.2,
                "{:?} still points at the mesh of the previous field",
                before.0
            );
        }
    }

    /// The panel shows the at-floor count beside the aspect ratio, and calls the
    /// read-off falsified.
    ///
    /// **This is the only way to see this demo's screen on a machine with no
    /// display**, and the three things it checks are the whole deliverable of
    /// the HUD:
    ///
    /// 1. **`aspect_ratio_max` never appears without an at-floor count.** At a
    ///    flat direction the aspect ratio is `|lambda|max / H_FLOOR`, so on
    ///    `thin_plate` — where nearly every band cell is at the floor — the
    ///    `1e10` on screen is a restatement of `H_FLOOR = 1e-9` and not a
    ///    measured anisotropy. Separate the maximum from its denominator and the
    ///    number lies, which is why this asserts on the *same line* carrying
    ///    both.
    /// 2. **The verdict word is the measured one.** P-146 C1 was falsified and a
    ///    demo implying a win the measurement did not find is worse than no demo,
    ///    so the ratio must be above C1's bar and the line must say so. If a
    ///    future change made the metric arm win here, this fails and someone has
    ///    to reconcile the demo with the committed finding rather than discover
    ///    the disagreement from a GIF.
    /// 3. **The citation is the CSV's text.** The two figures are asserted as
    ///    literals rather than through [`CITED`], because a test that reads the
    ///    same constant the code prints cannot catch a mis-transcribed cell.
    #[test]
    fn the_hud_reports_the_at_floor_count_beside_the_aspect_ratio_and_a_falsified_verdict() {
        let mut app = harness(THIN_PLATE);
        settle(&mut app);
        // `report` is left out of `harness` because it wants a `DemoStats`;
        // this runs it as a one-shot with one inserted, which is the same
        // system the demo runs every frame.
        app.init_resource::<DemoStats>();
        app.world_mut()
            .run_system_once(report)
            .expect("the HUD system");

        let stats = app.world().resource::<DemoStats>();
        let lines = stats.extra.clone();
        let title = stats.title.clone();
        let banner = stats.banner.clone().expect("the headline");
        println!("{title}");
        for line in &lines {
            println!("{line}");
        }
        println!("banner: {}", banner.0);

        // The caption has to name the field the numbers below it belong to, or a
        // reader walking the roster cannot tell which row they are looking at --
        // and a caption disagreeing with the panel above it by one frame is the
        // failure `game_mirror_dedup.rs:1350-1353` records.
        assert!(
            title.contains("E-319") && title.contains(CITED[THIN_PLATE].name),
            "the caption does not name its ticket and its field: {title}"
        );

        let (points, at_floor, ratio) = {
            let readout = app.world().resource::<Readout>();
            let Readout::Measured(measured) = readout else {
                panic!("the ladder settled but the readout is not a measurement");
            };
            let rung = measured.shown_rung();
            let matched = measured
                .matched
                .as_ref()
                .expect("thin_plate is Exact, so it has a matched read-off");
            (rung.census.points, rung.census.at_floor, matched.ratio)
        };

        let aspect = lines
            .iter()
            .find(|line| line.contains("aspect ratio"))
            .expect("the aspect-ratio line");
        assert!(
            aspect.contains(&format!("at floor {at_floor} of {points} band cells")),
            "the aspect-ratio maximum lost its at-floor count, so it is a restatement of \
             H_FLOOR wearing a measurement's clothes: {aspect}"
        );
        assert!(
            lines.iter().any(|line| line.contains("H_FLOOR")),
            "no line explains that an at-floor aspect ratio is the floor talking"
        );

        assert!(
            ratio > WIN_RATIO,
            "the live read-off says the metric arm saved triangles on thin_plate (ratio \
             {ratio}), where p-146's committed `ratio` column says 2.298935. One of the two \
             is stale and a demo must not be the one claiming the win"
        );
        let verdict = lines
            .iter()
            .find(|line| line.contains("P-146 C1"))
            .expect("the C1 verdict line");
        assert!(
            verdict.contains("FALSIFIED"),
            "the C1 verdict line stopped saying what the measurement found: {verdict}"
        );
        assert!(
            banner.0.contains("NO SAVING") && banner.1 == FALSIFIED,
            "the headline does not carry the falsified verdict: {}",
            banner.0
        );

        let cited = lines
            .iter()
            .find(|line| line.contains("p-146.csv"))
            .expect("the p-146 citation line");
        assert!(
            cited.contains("2.298935") && cited.contains("0.695277"),
            "the citation line stopped quoting p-146's thin_plate row: {cited}"
        );
        let global = lines
            .iter()
            .find(|line| line.contains("P-149 global"))
            .expect("the p-149 citation line");
        assert!(
            global.contains("70 of 112"),
            "the citation line stopped quoting p-149's identical-arms count: {global}"
        );
    }

    /// P-146's registered vacuity control on the anisotropic arm: `axis_ratio`
    /// must exceed this on at least one field, or the two arms are the same grid
    /// under two names and every ratio above is measuring noise.
    const AXIS_RATIO_FLOOR: f64 = 1.5;

    /// The registered control's bar on `aspect_ratio_max`, and the amendment
    /// that makes it informative: the maximum must clear this among the cells
    /// that are **not** at the floor.
    const ASPECT_FLOOR: f64 = 3.0;

    /// Every field builds, no field saves a triangle, and both of P-146's
    /// vacuity controls are met live.
    ///
    /// This is the gate that stops the demo from ever implying a win the
    /// measurement did not find. P-146's `c1_winners` column reads `0` against a
    /// `c1_population` of `4` on all forty rows; recomputed here over the whole
    /// roster at the shown rung, the count of measurable fields must still be
    /// four and the count of winners must still be zero. A change that makes one
    /// of them win fails this test, and then either the demo or the committed
    /// finding is stale and a person has to say which — rather than the
    /// disagreement being discovered from a GIF.
    ///
    /// Both controls are asserted because `M-44`'s rule cuts both ways here. The
    /// arms must be genuinely different somewhere (`axis_ratio`), or "no saving"
    /// is a statement about a grid that was never anisotropic; and the aspect
    /// ratio must clear its bar among cells that are **off** the floor, because
    /// a maximum met at the floor is a restatement of `H_FLOOR`
    /// (`benches/common/metric.rs:67-74`). The expected witnesses are
    /// `fbm_terrain` and `thin_plate` for the first and `gyroid` for the second,
    /// the last being Lipschitz rather than a distance field and so having no
    /// exactly-flat direction to floor.
    ///
    /// It also asserts that **both** readings are reachable — a field whose arms
    /// come out identical and a field whose arms do not. Five of the eight are
    /// identical, which is P-149's 70 of 112 rows, and a demo that could only
    /// ever show one of the two would misrepresent how often each arrives.
    #[test]
    fn every_field_reproduces_p146s_verdict_and_its_vacuity_controls() {
        let mut measurable = 0usize;
        let mut winners = 0usize;
        let mut anisotropic_somewhere = false;
        let mut off_floor_anisotropy_somewhere = false;
        let mut identical = 0usize;

        for (index, cited) in CITED.iter().enumerate() {
            let (measured, meshes) = build((index, DEFAULT_RUNG))
                .unwrap_or_else(|reason| panic!("{}: {reason}", cited.name));
            let rung = measured.shown_rung();
            let ratio = measured.matched.as_ref().map_or_else(
                |reason| (*reason).to_string(),
                |m| format!("{:.4}", m.ratio),
            );
            println!(
                "{:>14} {}^3  grid {}x{}x{}  budget {:.3}  axis {:>6.2}  tri {:>6}/{:<6} \
                 aspect {:>9.3e} ({} of {} at floor, off-floor {:>9.3e})  flat {:<5} \
                 share {:>6.3} in [{:.3}, {:.3}]  samples {:>6.3}  ratio {}",
                cited.name,
                rung.samples,
                rung.grid[0],
                rung.grid[1],
                rung.grid[2],
                rung.budget_ratio,
                rung.axis_ratio,
                rung.triangles[Side::Uniform.index()],
                rung.triangles[Side::Metric.index()],
                rung.census.aspect[0],
                rung.census.at_floor,
                rung.census.points,
                rung.census.aspect_off_floor[0],
                flat_axis(&rung.census),
                measured.cost.share(),
                measured.cost.share_span()[0],
                measured.cost.share_span()[1],
                measured.cost.sample_share(),
                ratio,
            );

            // The mesh handed to the renderer is the mesh that was measured. A
            // side-by-side whose picture and numbers come from different
            // extractions is not a comparison.
            for side in [Side::Uniform, Side::Metric] {
                assert_eq!(
                    meshes[side.index()].count_vertices(),
                    rung.vertices[side.index()],
                    "{}: the drawn mesh is not the measured one",
                    cited.name
                );
            }

            assert!(
                (0.25..=4.0).contains(&rung.budget_ratio),
                "{}: the metric arm spent {:.3}x the uniform arm's budget, so the two arms \
                 are two grids rather than two shapes and no triangle count compares",
                cited.name,
                rung.budget_ratio
            );
            assert!(
                rung.census.points > 0,
                "{}: no grid sample landed within {BAND_CELLS} cell of the surface, so every \
                 statistic on this row is a statistic of nothing (M-44)",
                cited.name
            );
            // C3 gated on a **count**, not on the clock. p-146 falsified it on
            // all forty of its rows, and the number that says so without asking
            // what governor this machine is running is `19 x band` against the
            // grid's `N^3`. It is a floor on the metric's cost -- it does not
            // charge the eigendecompositions -- so clearing the bar on the floor
            // is the stronger statement. This is the assertion that catches a
            // single-shot clock printing HOLDS for a falsified clause.
            assert!(
                measured.cost.sample_share() > C3_MAX_SHARE,
                "{}: the metric read {} field samples against the grid's {}, {:.4}x, which is \
                 inside C3's bar of {C3_MAX_SHARE} -- yet p-146's c3_holds column reads false \
                 on all forty of its rows",
                cited.name,
                measured.cost.metric_samples,
                measured.cost.grid_samples,
                measured.cost.sample_share()
            );

            if let Ok(matched) = &measured.matched {
                measurable += 1;
                if matched.ratio <= WIN_RATIO {
                    winners += 1;
                }
            }
            if rung.axis_ratio > AXIS_RATIO_FLOOR {
                anisotropic_somewhere = true;
            }
            if rung.census.at_floor < rung.census.points
                && rung.census.aspect_off_floor[0] > ASPECT_FLOOR
            {
                off_floor_anisotropy_somewhere = true;
            }
            if rung.arms_identical() {
                identical += 1;
            }
        }

        assert_eq!(
            measurable, 4,
            "p-146's c1_population is 4 -- the four reference fields whose bound() is Exact -- \
             and this roster now offers a different number of Hausdorff measurements"
        );
        assert_eq!(
            winners,
            0,
            "a field saved at least {}% of its triangles at matched Hausdorff, where p-146's \
             c1_winners column reads 0 of 4 on all forty rows. A demo implying a win the \
             measurement did not find is worse than no demo, so this is a stop rather than a \
             pleasant surprise",
            (1.0 - WIN_RATIO) * 100.0
        );
        assert!(
            anisotropic_somewhere,
            "no field's metric grid exceeded an axis ratio of {AXIS_RATIO_FLOOR}, so the \
             'metric' arm is the uniform arm everywhere and the verdict above is noise"
        );
        assert!(
            off_floor_anisotropy_somewhere,
            "no field carried an aspect ratio above {ASPECT_FLOOR} among cells that are not at \
             the floor, so every anisotropy on screen is a restatement of H_FLOOR = {H_FLOOR:e}"
        );
        assert!(
            identical > 0 && identical < CITED.len(),
            "{identical} of {} fields had identical arms; the demo must be able to show both \
             readings, because P-149 found 70 of its 112 rows identical and a demo that only \
             ever showed one of the two would misrepresent how often each arrives",
            CITED.len()
        );
    }

    /// The metric arm spends the uniform arm's budget, and equal weights give
    /// back the uniform grid exactly.
    ///
    /// Both halves are load-bearing and neither needs a field.
    ///
    /// The first is what makes a triangle-count comparison mean anything: two
    /// arms on different budgets are two grids, not two shapes. P-146's own
    /// vacuity control puts `∏ n_a / N³` in `[0.25, 4]`, and its harness caught a
    /// version of [`metric_grid`] that was **5.77x** over on a weight vector like
    /// the third case here.
    ///
    /// The second is what makes the demo's `IDENTICAL GRIDS` headline a
    /// measurement rather than a rounding artifact: on five of the eight
    /// reference fields the metric prescribes the uniform grid, and that claim is
    /// only worth making if an isotropic weight vector reproduces `(N, N, N)`
    /// *exactly*.
    #[test]
    fn the_metric_grid_spends_the_uniform_arms_budget() {
        for samples in LADDER {
            let (grid, pinned) = metric_grid([1.0, 1.0, 1.0], samples);
            assert_eq!(
                grid, [samples; 3],
                "equal point densities did not reproduce the uniform grid at {samples}^3, so \
                 an IDENTICAL GRIDS reading would be a rounding artifact"
            );
            assert_eq!(
                pinned, 0,
                "an isotropic split pinned an axis at {samples}^3"
            );
        }

        // A heightfield's `y` weight is ~1e-4 of the other two (`fbm_terrain`);
        // the second case is a slab-like anisotropy; the first is a mild one.
        for weights in [
            [1.0, 0.5, 1.0],
            [1.0, 0.05, 1.0],
            [1.0, 1e-4, 1.0],
            [3.0, 1e-6, 2.0],
        ] {
            for samples in LADDER {
                let (grid, _) = metric_grid(weights, samples);
                let product = f64::from(grid[0]) * f64::from(grid[1]) * f64::from(grid[2]);
                let budget = product / f64::from(samples).powi(3);
                assert!(
                    (0.25..=4.0).contains(&budget),
                    "weights {weights:?} at {samples}^3 gave {grid:?}, {budget}x the budget: \
                     the two arms are then two grids rather than two shapes"
                );
                for axis in grid {
                    assert!(
                        axis >= MIN_SAMPLES,
                        "weights {weights:?} at {samples}^3 gave {grid:?}, below the floor"
                    );
                    assert_eq!(
                        axis % 2,
                        1,
                        "weights {weights:?} at {samples}^3 gave {grid:?}, and an even axis \
                         count loses thin_plate's surface outright (M-266)"
                    );
                }
            }
        }
    }

    /// The capture walk visits every field exactly once and then loops.
    ///
    /// The GIF is the only form in which most readers will meet this demo, and
    /// its whole content is the walk. `M-241` records a clip that *advertised a
    /// sweep it never performed* and a check that passed on a single still, so
    /// the arithmetic gets a gate: `record_gif.sh`'s default eighty captured
    /// frames must cover all eight fields, land on each for the same number of
    /// frames, and return to the first so the loop does not jump.
    #[test]
    fn the_capture_walk_covers_the_roster_once_and_loops() {
        let total = CAPTURE_FRAMES_PER_FIELD * CITED.len() as u32;
        assert_eq!(
            total, 80,
            "the walk no longer fills record_gif.sh's default ISOMESH_CAPTURE_FRAMES of 80, so \
             the clip would stop part-way through the roster or repeat part of it"
        );

        let mut frames = vec![0u32; CITED.len()];
        for taken in 0..total {
            frames[walked_field(taken)] += 1;
        }
        assert!(
            frames
                .iter()
                .all(|&count| count == CAPTURE_FRAMES_PER_FIELD),
            "the walk does not spend equal time on each field: {frames:?}"
        );
        assert_eq!(
            walked_field(0),
            0,
            "the walk does not start on the first field"
        );
        assert_eq!(
            walked_field(total),
            walked_field(0),
            "the walk does not return to its first field, so the GIF jumps at the loop point"
        );
    }
}
