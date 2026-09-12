//! **P-177 — this crate's reach, both halves of it, and whether it is the denominator for resolution.**
//!
//! Ticket: `R-182`. Pre-registered before this harness existed; minted by the
//! discovery loop (`docs/research/2026-09-12-discovery-loop-prompt.md`,
//! iteration 1) from the question ledger's `Q1`.
//!
//! Writes `docs/experiments/p-177.csv`.
//!
//! # The source, and the two halves it licenses
//!
//! Aamari, Kim, Chazal, Michel, Rinaldo & Wasserman, *Estimating the reach of a
//! manifold*, `10.1214/19-ejs1551` — in the corpus, converted, 43 pp., and
//! uncited in this repository before this row. **Theorem 3.4**: for a compact
//! submanifold with reach `τ > 0`, at least one of
//!
//! - **(Global Case)** a *bottleneck* — Definition 3.1's reach-attaining pair
//!   with `‖q₁ − q₂‖ = 2τ`, the two points meeting the same medial ball; and
//! - **(Local Case)** *"there exists `q₀ ∈ M` and an arc-length parametrized
//!   geodesic `γ₀` such that `γ₀(0) = q₀` and `‖γ₀″(0)‖ = 1/τ`"* — a point whose
//!   radius of curvature **is** the reach,
//!
//! holds, and the paper is explicit that *"the global case and the local case of
//! Theorem 3.4 are not mutually exclusive"*. That is what makes the reach
//! computable here **without a medial axis**, which matters because the direct
//! route is closed: `✗27` falsified *"an exactly-zero SDF gradient detects the
//! medial axis"* and `M-325`'s discrete medial score converges only `O(h)`.
//!
//! # What is measured
//!
//! - `tau_local` — `1 / max|κ|` over the extracted vertices, the principal
//!   curvatures being the eigenvalues of the shape operator `(I − nnᵀ)H/|∇f|`
//!   built from central differences of the field's own gradient. This is the
//!   Local Case read off the field rather than off the triangles.
//! - `tau_global` — half the shortest segment joining two surface points whose
//!   normals are **collinear with it at both ends**, in either of the two
//!   orientations (normals facing each other across a gap, or away from each
//!   other across a sheet). Both are reach-attaining configurations: the inner
//!   and outer branches of the medial axis. This is the Global Case, by
//!   Definition 3.1, with no medial axis constructed.
//! - `h_star` — the coarsest spacing whose extracted topology, `(χ, components)`
//!   from `validate_indexed`, equals the finest rung's *and keeps equalling it at
//!   every finer rung*. A field whose topology never settles inside the ladder
//!   records `h_star_reachable = false` and is not scored.
//!
//! # The vacuity control, which runs before any of the eight numbers is read
//!
//! Four closed forms, asserted:
//!
//! | field | `τ` | why |
//! |---|---|---|
//! | `sphere` | **1.0** | `Sphere::canonical` is the unit sphere; a sphere's reach is its radius |
//! | `torus` | **0.3** | `Torus::canonical` is major 1, minor 0.3; `min(minor, major − minor) = 0.3` |
//! | `capsule` | **0.35** | `brush::Capsule`, built here for the control alone: a capsule's reach is its radius |
//! | `box_exact` | **0.0** | a sharp edge has no reach at all, so the estimator must report zero to within the grid |
//!
//! `tau_analytic` and the measured `tau` are both columns, `calibration_gap_h`
//! is their difference in units of `h`, and the run refuses to continue if any of
//! the four misses by more than `2h`. Without that, the other four fields'
//! numbers would be unanchored.
//!
//! # Determinism
//!
//! No RNG. Every field is deterministic, every ladder is arithmetic, and the two
//! reach halves are minima over fixed index ranges compared with `total_cmp`.

#![allow(
    clippy::needless_range_loop,
    clippy::cast_possible_truncation,
    clippy::cast_precision_loss,
    // The AUC's tie detection is an exact-equality question about two reads of
    // one array, which is what `P-173`'s copy of this function also allows.
    clippy::float_cmp,
    clippy::too_many_lines
)]

mod common;

use isomesh::fields::{
    BoxExact, FbmTerrain, ReferenceField, Sphere, ThinPlate, Torus, capped_gyroid, csg_difference,
    noise_cavity,
};
use isomesh::marching_cubes::MarchingCubes;
use isomesh::marching_cubes::ambiguity::joined_mask;
use isomesh::marching_cubes::table::{NO_EDGE, is_inside, segment_links};
use isomesh::validate::{ValidateConfig, validate_indexed};
use isomesh::{MeshBuffer, RuntimeShape3, Sdf};

/// The resolution ladder, in samples per axis. Odd throughout, because `M-266`
/// proved the canonical grids' odd counts are load-bearing.
const LADDER: [u32; 5] = [17, 25, 33, 49, 65];

/// Where `C3`'s per-cell census runs: `P-173` used the same rung for the same
/// two fields, so the populations are comparable.
const C3_SAMPLES: u32 = 33;

/// The step the field's gradient and Hessian are differenced at.
///
/// **Not the grid spacing, and the difference is the instrument's whole
/// correctness.** The reach is a property of the *field's* zero set, and every
/// field here is an analytic function that can be queried at any scale — the
/// grid is where the mesh lives, not where the geometry is. Differencing at `h`
/// caps the measurable curvature at `O(1/h)` and so puts a floor of several
/// cells under any measured reach: the first run of this harness did exactly
/// that and read `tau` **0.25** on `box_exact`, whose reach is **0**, at
/// `h` **0.0625**. `max_abs_curvature_grid_step` records that quantity beside
/// the real one so the difference is visible in the CSV rather than only here.
const CURVATURE_STEP: f64 = 1e-3;

/// How nearly collinear a segment must be with both normals to count as a
/// bottleneck: `cos θ ≥ 1 − COLLINEAR_SLACK`.
///
/// A tolerance is unavoidable — the vertices are a sampling of the surface, not
/// the surface — and it is a recorded number rather than a hidden one. At `0.02`
/// the admitted cone is about 11.5°.
const COLLINEAR_SLACK: f64 = 0.02;

/// How far above the shortest bottleneck a pair may sit and still be kept as a
/// medial-axis sample for `C3`'s score.
const AXIS_BAND: f64 = 4.0;

/// `C3`'s registered bar.
const AUC_BAR: f64 = 0.8;

/// `C2`'s registered bar on `max(h*/τ) / min(h*/τ)`.
const SPREAD_BAR: f64 = 2.0;

/// Vertices above which the `O(n²)` bottleneck search strides the outer loop.
///
/// The minimum is over pairs, so striding the *outer* loop keeps every candidate
/// partner for the points it does visit; `bottleneck_pairs` records how many
/// qualifying pairs were seen so a thinned search cannot be mistaken for a dense
/// one.
const BOTTLENECK_BUDGET: usize = 20_000;

/// A field's sampled grid and the mesh extracted from it.
struct Extraction {
    /// Samples per axis.
    samples: u32,
    /// Grid spacing.
    h: f64,
    /// The sampled values, `x` fastest. Kept so a reader of the CSV can see the
    /// grid the topology came from.
    #[allow(
        dead_code,
        reason = "the grid is carried for readers of the CSV's provenance"
    )]
    values: Vec<f64>,
    /// Mesh vertex positions.
    positions: Vec<[f64; 3]>,
    /// Triangle indices, carried beside the positions they index.
    #[allow(
        dead_code,
        reason = "the mesh is carried whole; only positions are read here"
    )]
    indices: Vec<u32>,
    /// `validate_indexed`'s Euler characteristic.
    euler: i64,
    /// `validate_indexed`'s component count.
    components: u64,
}

/// Sample `field` on its own domain and extract at `samples³`.
fn extract<S>(field: &S, samples: u32, lo: [f64; 3], hi: [f64; 3]) -> Extraction
where
    S: Sdf<Scalar = f64>,
{
    let h = (hi[0] - lo[0]) / f64::from(samples - 1);
    let n = samples as usize;
    let mut values = vec![0.0f64; n * n * n];
    for z in 0..n {
        for y in 0..n {
            let row = n * (y + n * z);
            for x in 0..n {
                values[row + x] = field.sample([
                    h.mul_add(x as f64, lo[0]),
                    h.mul_add(y as f64, lo[1]),
                    h.mul_add(z as f64, lo[2]),
                ]);
            }
        }
    }

    let shape = RuntimeShape3::new([samples; 3]).expect("the ladder fits u32");
    let mut mesh = MeshBuffer::<f64>::new();
    MarchingCubes::new()
        .extract(field, &shape, lo, h, &mut mesh)
        .expect("the ladder's grids are extractable");
    let config = ValidateConfig::from_cell_size(h).expect("the ladder's spacings are positive");
    let report = validate_indexed(&mesh.positions, &mesh.indices, &config);

    Extraction {
        samples,
        h,
        values,
        positions: mesh.positions.clone(),
        indices: mesh.indices.clone(),
        euler: report.euler_characteristic,
        components: report.components,
    }
}

/// The field's gradient at `p`, by central differences at `h`.
fn gradient<S: Sdf<Scalar = f64>>(field: &S, p: [f64; 3], h: f64) -> [f64; 3] {
    let mut out = [0.0; 3];
    for axis in 0..3 {
        let (mut lo, mut hi) = (p, p);
        lo[axis] -= h;
        hi[axis] += h;
        out[axis] = (field.sample(hi) - field.sample(lo)) / (2.0 * h);
    }
    out
}

/// The largest absolute principal curvature of the level set through `p`.
///
/// The shape operator of a level set is `(I − n nᵀ) H / |∇f|` with `H` the
/// Hessian; its two non-trivial eigenvalues are the principal curvatures. The
/// largest is obtained from the trace and the Frobenius norm of the tangential
/// part rather than from an eigen-solver: for a symmetric `2×2` block,
/// `max|κ| = |tr|/2 + sqrt(max(0, ‖·‖²/2 − tr²/4))`.
fn max_abs_curvature<S: Sdf<Scalar = f64>>(field: &S, p: [f64; 3], h: f64) -> Option<f64> {
    let g = gradient(field, p, h);
    let norm = (g[0] * g[0] + g[1] * g[1] + g[2] * g[2]).sqrt();
    if !norm.is_finite() || norm < 1e-9 {
        return None;
    }
    let n = [g[0] / norm, g[1] / norm, g[2] / norm];

    let mut hessian = [[0.0f64; 3]; 3];
    for axis in 0..3 {
        let (mut lo, mut hi) = (p, p);
        lo[axis] -= h;
        hi[axis] += h;
        let glo = gradient(field, lo, h);
        let ghi = gradient(field, hi, h);
        for k in 0..3 {
            hessian[axis][k] = (ghi[k] - glo[k]) / (2.0 * h);
        }
    }
    // Symmetrise: central differences of a gradient are only symmetric up to
    // truncation, and an asymmetric "Hessian" has complex eigenvalues.
    for a in 0..3 {
        for b in 0..a {
            let mean = (hessian[a][b] + hessian[b][a]) / 2.0;
            hessian[a][b] = mean;
            hessian[b][a] = mean;
        }
    }

    // S = (I − n nᵀ) H (I − n nᵀ) / |∇f|, whose rank is at most two and whose
    // non-zero eigenvalues are the principal curvatures.
    let project = |v: [f64; 3]| {
        let dot = v[0] * n[0] + v[1] * n[1] + v[2] * n[2];
        [v[0] - dot * n[0], v[1] - dot * n[1], v[2] - dot * n[2]]
    };
    let mut shape = [[0.0f64; 3]; 3];
    for a in 0..3 {
        let column = project([hessian[0][a], hessian[1][a], hessian[2][a]]);
        for b in 0..3 {
            shape[b][a] = column[b];
        }
    }
    for a in 0..3 {
        let row = project([shape[a][0], shape[a][1], shape[a][2]]);
        shape[a] = row;
    }
    let trace = (shape[0][0] + shape[1][1] + shape[2][2]) / norm;
    let mut frobenius = 0.0;
    for a in 0..3 {
        for b in 0..3 {
            frobenius += (shape[a][b] / norm) * (shape[a][b] / norm);
        }
    }
    let discriminant = (frobenius / 2.0 - trace * trace / 4.0).max(0.0);
    let magnitude = trace.abs() / 2.0 + discriminant.sqrt();
    magnitude.is_finite().then_some(magnitude)
}

/// The Global Case: half the shortest bottleneck, and the medial midpoints found.
///
/// Returns `(tau_global, qualifying_pairs, axis_points)`. `tau_global` is
/// infinite when the surface has no bottleneck at all, which is the correct
/// answer for a sphere and is what makes it curvature-bound.
fn bottleneck<S: Sdf<Scalar = f64>>(
    field: &S,
    positions: &[[f64; 3]],
    step: f64,
    h: f64,
) -> (f64, u64, Vec<[f64; 3]>) {
    let normals: Vec<Option<[f64; 3]>> = positions
        .iter()
        .map(|p| {
            let g = gradient(field, *p, step);
            let norm = (g[0] * g[0] + g[1] * g[1] + g[2] * g[2]).sqrt();
            (norm.is_finite() && norm > 1e-9).then(|| [g[0] / norm, g[1] / norm, g[2] / norm])
        })
        .collect();

    let stride = positions.len().div_ceil(BOTTLENECK_BUDGET).max(1);
    let mut best = f64::INFINITY;
    let mut pairs = 0u64;
    let mut axis_points: Vec<([f64; 3], f64)> = Vec::new();
    let threshold = 1.0 - COLLINEAR_SLACK;
    let mut i = 0usize;
    while i < positions.len() {
        let Some(ni) = normals[i] else {
            i += stride;
            continue;
        };
        let p = positions[i];
        for j in 0..positions.len() {
            if j == i {
                continue;
            }
            let Some(nj) = normals[j] else { continue };
            let q = positions[j];
            let d = [q[0] - p[0], q[1] - p[1], q[2] - p[2]];
            let length = (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt();
            // Two samples of the same triangle are not a bottleneck; anything
            // shorter than one cell is below what the grid can resolve.
            if length <= h {
                continue;
            }
            let u = [d[0] / length, d[1] / length, d[2] / length];
            let ci = u[0] * ni[0] + u[1] * ni[1] + u[2] * ni[2];
            let cj = u[0] * nj[0] + u[1] * nj[1] + u[2] * nj[2];
            // Either orientation: normals facing away from each other across a
            // sheet, or toward each other across a gap. Both are reach-attaining.
            let collinear =
                (ci >= threshold && cj <= -threshold) || (ci <= -threshold && cj >= threshold);
            if !collinear {
                continue;
            }
            pairs += 1;
            let half = length / 2.0;
            if half < best {
                best = half;
            }
            axis_points.push((
                [
                    (p[0] + q[0]) / 2.0,
                    (p[1] + q[1]) / 2.0,
                    (p[2] + q[2]) / 2.0,
                ],
                half,
            ));
        }
        i += stride;
    }
    if best.is_finite() {
        // Keep the midpoints of the SHORTER bottlenecks: those are the
        // medial-axis samples C3's score is a distance to. A midpoint of a pair
        // far above the minimum is still on the medial axis, but it is a branch
        // the reach does not run through.
        let cut = best * AXIS_BAND;
        let keep: Vec<[f64; 3]> = axis_points
            .into_iter()
            .filter(|(_, half)| *half <= cut)
            .map(|(m, _)| m)
            .collect();
        (best, pairs, keep)
    } else {
        (f64::INFINITY, pairs, Vec::new())
    }
}

/// How many vertices `ManifoldDualContouring` would place in this cell.
///
/// `P-173`'s helper, copied: `segment_links(case, joined_mask(corner, 0))` under
/// the default `FaceAmbiguity::Separate`, one vertex per cycle.
fn cycle_count(case: u8, corner: &[f64; 8]) -> u32 {
    let next = segment_links(case, joined_mask(corner, 0));
    let mut visited = 0u16;
    let mut cycles = 0;
    for (start, &link) in next.iter().enumerate() {
        if link == NO_EDGE || visited & (1 << start) != 0 {
            continue;
        }
        let mut current = start as u8;
        while visited & (1 << current) == 0 {
            visited |= 1 << current;
            current = next[current as usize];
        }
        cycles += 1;
    }
    cycles
}

/// Mann–Whitney `AUC`, `P-173`'s statistic, copied so the two rows compare.
fn mann_whitney_auc(scores: &[f64], positive: &[bool]) -> Option<f64> {
    assert_eq!(scores.len(), positive.len(), "one label per score");
    let mut order: Vec<u32> = (0..scores.len() as u32).collect();
    order.sort_by(|a, b| scores[*a as usize].total_cmp(&scores[*b as usize]));
    let mut rank = vec![0.0; scores.len()];
    let mut i = 0usize;
    while i < order.len() {
        let mut j = i + 1;
        while j < order.len() && scores[order[j] as usize] == scores[order[i] as usize] {
            j += 1;
        }
        let average = (i + 1 + j) as f64 / 2.0;
        for slot in &order[i..j] {
            rank[*slot as usize] = average;
        }
        i = j;
    }
    let mut rank_sum = 0.0;
    let mut positives = 0usize;
    for (slot, &is_pos) in positive.iter().enumerate() {
        if is_pos {
            rank_sum += rank[slot];
            positives += 1;
        }
    }
    let negatives = positive.len() - positives;
    if positives == 0 || negatives == 0 {
        return None;
    }
    let u = rank_sum - (positives as f64) * (positives as f64 + 1.0) / 2.0;
    Some(u / (positives as f64 * negatives as f64))
}

/// `C3`: do the second-vertex cells sit nearer a bottleneck axis than the rest?
fn second_vertex_auc<S>(
    field: &S,
    axis_points: &[[f64; 3]],
    lo: [f64; 3],
    hi: [f64; 3],
) -> (Option<f64>, u64)
where
    S: Sdf<Scalar = f64>,
{
    if axis_points.is_empty() {
        return (None, 0);
    }
    let n = C3_SAMPLES as usize;
    let h = (hi[0] - lo[0]) / f64::from(C3_SAMPLES - 1);
    let mut values = vec![0.0f64; n * n * n];
    for z in 0..n {
        for y in 0..n {
            let row = n * (y + n * z);
            for x in 0..n {
                values[row + x] = field.sample([
                    h.mul_add(x as f64, lo[0]),
                    h.mul_add(y as f64, lo[1]),
                    h.mul_add(z as f64, lo[2]),
                ]);
            }
        }
    }

    let mut scores = Vec::new();
    let mut positive = Vec::new();
    let mut second = 0u64;
    for z in 0..n - 1 {
        for y in 0..n - 1 {
            for x in 0..n - 1 {
                let mut corner = [0.0f64; 8];
                for c in 0..8 {
                    let o = [c & 1, (c >> 1) & 1, (c >> 2) & 1];
                    corner[c] = values[(x + o[0]) + n * ((y + o[1]) + n * (z + o[2]))];
                }
                let mut case = 0u8;
                for c in 0..8 {
                    if is_inside(corner[c]) {
                        case |= 1 << c;
                    }
                }
                if case == 0 || case == u8::MAX {
                    continue;
                }
                let centre = [
                    h.mul_add(x as f64 + 0.5, lo[0]),
                    h.mul_add(y as f64 + 0.5, lo[1]),
                    h.mul_add(z as f64 + 0.5, lo[2]),
                ];
                let mut nearest = f64::INFINITY;
                for m in axis_points {
                    let d = [centre[0] - m[0], centre[1] - m[1], centre[2] - m[2]];
                    nearest = nearest.min(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]);
                }
                let cycles = cycle_count(case, &corner);
                if cycles > 1 {
                    second += 1;
                }
                // Nearer must score HIGHER, or the AUC reads the hypothesis
                // backwards: negate the distance.
                scores.push(-nearest.sqrt());
                positive.push(cycles > 1);
            }
        }
    }
    (mann_whitney_auc(&scores, &positive), second)
}

/// One field's whole story.
struct Story {
    /// Name.
    field: &'static str,
    /// Closed-form reach, where one exists.
    analytic: Option<f64>,
    /// Local Case.
    tau_local: f64,
    /// Global Case.
    tau_global: f64,
    /// `min` of the two.
    tau: f64,
    /// Which half binds.
    binding: &'static str,
    /// Qualifying bottleneck pairs seen.
    pairs: u64,
    /// `max|κ|` over the finest rung's vertices.
    curvature: f64,
    /// The same maximum differenced at the **grid** step, for comparison: this
    /// is the quantity that put a floor of several cells under the first run.
    curvature_grid: f64,
    /// Coarsest settled spacing, if the ladder settles.
    h_star: Option<f64>,
    /// The finest rung's topology, for the record.
    euler: i64,
    /// Components at the finest rung.
    components: u64,
    /// `C3`'s AUC, on the two fields it names.
    auc: Option<f64>,
    /// Second-vertex cells at `C3_SAMPLES`.
    second: u64,
    /// The finest rung's spacing.
    h: f64,
}

/// `f64` for the CSV, with infinities spelled rather than printed as `inf`.
fn num(value: f64) -> String {
    if value.is_infinite() {
        String::from("unbounded")
    } else if value.is_nan() {
        String::from("na")
    } else {
        format!("{value:.9}")
    }
}

fn main() {
    if !std::env::args().any(|arg| arg == "--bench") {
        return;
    }
    let prereg = isomesh::experiment!("P-177");

    common::experiment::run(prereg, |run| {
        let mut stories: Vec<Story> = Vec::new();

        macro_rules! study {
            ($name:expr, $field:expr, $domain:expr, $analytic:expr, $c3:expr) => {{
                let field = $field;
                let (lo, hi) = $domain;
                let finest = extract(&field, LADDER[LADDER.len() - 1], lo, hi);
                let mut curvature = 0.0f64;
                let mut curvature_grid = 0.0f64;
                for p in &finest.positions {
                    if let Some(k) = max_abs_curvature(&field, *p, CURVATURE_STEP) {
                        curvature = curvature.max(k);
                    }
                    if let Some(k) = max_abs_curvature(&field, *p, finest.h) {
                        curvature_grid = curvature_grid.max(k);
                    }
                }
                let tau_local = if curvature > 0.0 {
                    1.0 / curvature
                } else {
                    f64::INFINITY
                };
                let (tau_global, pairs, axis_points) =
                    bottleneck(&field, &finest.positions, CURVATURE_STEP, finest.h);
                let tau = tau_local.min(tau_global);
                let binding = if tau_global < tau_local {
                    "global"
                } else if tau_local < tau_global {
                    "local"
                } else {
                    "tied"
                };

                // h*: the coarsest rung that matches the finest and keeps
                // matching at every finer rung.
                let target = (finest.euler, finest.components);
                let mut ladder: Vec<(f64, bool)> = Vec::new();
                for samples in LADDER {
                    let rung = if samples == finest.samples {
                        (finest.h, true)
                    } else {
                        let e = extract(&field, samples, lo, hi);
                        (e.h, (e.euler, e.components) == target)
                    };
                    ladder.push(rung);
                }
                let mut h_star = None;
                for index in 0..ladder.len() {
                    if ladder[index..].iter().all(|(_, matches)| *matches) {
                        h_star = Some(ladder[index].0);
                        break;
                    }
                }

                let (auc, second) = if $c3 {
                    second_vertex_auc(&field, &axis_points, lo, hi)
                } else {
                    (None, 0)
                };

                stories.push(Story {
                    field: $name,
                    analytic: $analytic,
                    tau_local,
                    tau_global,
                    tau,
                    binding,
                    pairs,
                    curvature,
                    curvature_grid,
                    h_star,
                    euler: finest.euler,
                    components: finest.components,
                    auc,
                    second,
                    h: finest.h,
                });
            }};
        }

        let sphere = Sphere::<f64>::canonical();
        let sphere_domain = sphere.domain();
        study!("sphere", sphere, sphere_domain, Some(1.0), false);
        let torus = Torus::<f64>::canonical();
        let torus_domain = torus.domain();
        study!("torus", torus, torus_domain, Some(0.3), false);
        let box_exact = BoxExact::<f64>::canonical();
        let box_domain = box_exact.domain();
        study!("box_exact", box_exact, box_domain, Some(0.0), false);
        // The registration's fourth closed form. `capsule` is NOT one of this
        // crate's eight reference fields -- it is `brush::Capsule`, which
        // `M-487` also used as a bench-local fixture -- so it is built here for
        // the calibration alone and carries `h_star` like any other row.
        study!(
            "capsule",
            isomesh::brush::Capsule::<f64> {
                a: [-0.5, 0.0, 0.0],
                b: [0.5, 0.0, 0.0],
                radius: 0.35,
            },
            ([-2.0f64; 3], [2.0f64; 3]),
            Some(0.35),
            false
        );
        let csg = csg_difference::<f64>();
        let csg_domain = csg.domain();
        study!("csg_difference", csg, csg_domain, None, false);
        let plate = ThinPlate::<f64>::canonical();
        let plate_domain = plate.domain();
        study!("thin_plate", plate, plate_domain, None, false);
        let gyroid = capped_gyroid::<f64>();
        let gyroid_domain = gyroid.domain();
        study!("gyroid", gyroid, gyroid_domain, None, true);
        let terrain = FbmTerrain::<f64>::canonical();
        let terrain_domain = terrain.domain();
        study!("fbm_terrain", terrain, terrain_domain, None, true);
        let cavity = noise_cavity::<f64>();
        let cavity_domain = cavity.domain();
        study!("noise_cavity", cavity, cavity_domain, None, false);

        // ── the vacuity control: scored, printed, and carried as a column ───
        //
        // Not a panic. `✗126 / M-490` established the rule the hard way: an
        // abort here deletes the CSV that is the evidence for *why* the control
        // could not be met, and a control that fails is a finding rather than a
        // reason to have no data.
        let mut calibration_failures: Vec<String> = Vec::new();
        for story in &stories {
            let Some(analytic) = story.analytic else {
                continue;
            };
            let gap = (story.tau.min(1e30) - analytic).abs();
            if gap > 2.0 * story.h {
                calibration_failures.push(format!(
                    "{}: analytic {analytic}, measured {}, gap {gap} > 2h {}",
                    story.field,
                    story.tau,
                    2.0 * story.h
                ));
            }
        }
        let calibration_holds = calibration_failures.is_empty();

        // ── the clauses ────────────────────────────────────────────────────
        let global_bound = stories.iter().filter(|s| s.binding == "global").count();
        let local_bound = stories.iter().filter(|s| s.binding == "local").count();
        let c1 = global_bound >= 2 && local_bound >= 2;

        let ratios: Vec<(&str, f64)> = stories
            .iter()
            .filter(|s| s.tau > 0.0 && s.tau.is_finite())
            .filter_map(|s| s.h_star.map(|h| (s.field, h / s.tau)))
            .collect();
        let spread = if ratios.len() < 2 {
            f64::INFINITY
        } else {
            let hi = ratios.iter().map(|(_, r)| *r).fold(f64::MIN, f64::max);
            let lo = ratios.iter().map(|(_, r)| *r).fold(f64::MAX, f64::min);
            hi / lo
        };
        let c2 = spread <= SPREAD_BAR;

        let c3_rows: Vec<&Story> = stories
            .iter()
            .filter(|s| s.field == "gyroid" || s.field == "fbm_terrain")
            .collect();
        let c3 = c3_rows.iter().all(|s| s.auc.is_some_and(|a| a >= AUC_BAR));

        println!(
            "\nC1 global-bound {global_bound} local-bound {local_bound} -> {c1}\n\
             C2 h*/tau spread {spread:.6} over {} fields -> {c2}\n\
             C3 {:?} -> {c3}",
            ratios.len(),
            c3_rows
                .iter()
                .map(|s| (s.field, s.auc, s.second))
                .collect::<Vec<_>>()
        );
        for (field, ratio) in &ratios {
            println!("   h*/tau  {field:<16} {ratio:.6}");
        }
        println!(
            "VACUITY (calibration): holds {calibration_holds}{}",
            if calibration_holds {
                String::new()
            } else {
                format!("\n   failures: {}", calibration_failures.join("  |  "))
            }
        );
        for story in &stories {
            println!(
                "   tau {:<16} local {:>14} global {:>14} binding {:<7} h* {}",
                story.field,
                num(story.tau_local),
                num(story.tau_global),
                story.binding,
                story.h_star.map_or_else(|| String::from("na"), num)
            );
        }

        for story in &stories {
            run.record(&[
                ("field", String::from(story.field)),
                ("resolution", LADDER[LADDER.len() - 1].to_string()),
                ("tau_local", num(story.tau_local)),
                ("tau_global", num(story.tau_global)),
                ("tau", num(story.tau)),
                ("binding_half", String::from(story.binding)),
                (
                    "tau_analytic",
                    story.analytic.map_or_else(|| String::from("na"), num),
                ),
                (
                    "calibration_gap_h",
                    story.analytic.map_or_else(
                        || String::from("na"),
                        |a| num((story.tau.min(1e30) - a).abs() / story.h),
                    ),
                ),
                ("bottleneck_pairs", story.pairs.to_string()),
                ("max_abs_curvature", num(story.curvature)),
                ("max_abs_curvature_grid_step", num(story.curvature_grid)),
                ("curvature_step", num(CURVATURE_STEP)),
                (
                    "h_star",
                    story.h_star.map_or_else(|| String::from("na"), num),
                ),
                (
                    "h_star_over_tau",
                    match story.h_star {
                        Some(h) if story.tau > 0.0 && story.tau.is_finite() => num(h / story.tau),
                        _ => String::from("na"),
                    },
                ),
                ("h_star_reachable", story.h_star.is_some().to_string()),
                ("topology_components", story.components.to_string()),
                ("euler_characteristic", story.euler.to_string()),
                ("second_vertex_cells", story.second.to_string()),
                (
                    "bottleneck_distance_auc",
                    story.auc.map_or_else(|| String::from("na"), num),
                ),
                ("calibration_holds", calibration_holds.to_string()),
                ("c1_holds", c1.to_string()),
                ("c2_holds", c2.to_string()),
                ("c3_holds", c3.to_string()),
            ]);
        }
    });
}
