//! **P-80 — the LOD residual as a normal map, so coarse chunks keep their
//! detail as shading.**
//!
//! Ticket: R-080. Pre-registered before this harness existed; the registration
//! is in `crates/isomesh/src/experiment.rs` and is not amended here.
//!
//! ```bash
//! cargo bench --bench experiment_p80
//! ```
//!
//! Writes `docs/experiments/p-80.csv`, one row per reference field per LOD level
//! 1–3.
//!
//! # Hypothesis
//!
//! `M-72` measured that a sub-cell feature does not vanish under coarsening, it
//! **aliases**. The registered response: the difference between the fine surface
//! and the coarse one is a vector field the crate can evaluate from the analytic
//! field alone, so bake it into a tangent-space normal map per coarse chunk and
//! the geometry fades while the shading does not.
//!
//! - **C1.** For each coarse-LOD vertex the direction to the nearest
//!   fine-surface point is computable from the field alone — one gradient and
//!   one root refinement, no fine mesh — and agrees with the true nearest point
//!   to within **0.1 cells** on **95%** of vertices at LOD 1 and 2.
//! - **C2.** Mean angular difference between LOD-2 shading-with-map and LOD-0
//!   shading is under **10°** on `fbm_terrain`, against over **25°** without.
//! - **C3.** It must fail on `thin_plate` and `gyroid` at LOD 3 — above **25°**
//!   with the map — because a normal map cannot restore a silhouette.
//!
//! **Vacuity control:** the LOD-0 reference must differ from LOD-2 by more than
//! 25° somewhere, reported as `changed_vertices` and asserted, or C2 is
//! comparing two identical shadings.
//!
//! # The LOD ladder is `M-72`'s, exactly
//!
//! `ChunkLayout::new(64 >> level, h₀, lo).at_lod(level)` over each field's own
//! `ReferenceField::domain`, which is the construction
//! `a_sub_cell_feature_aliases_under_coarsening_rather_than_vanishing` uses. The
//! world extent is therefore fixed and only the spacing changes, so every level
//! meshes the same region — and `thin_plate`'s triangle counts are directly
//! comparable with `M-72`'s **4,088 → 1,016 → 248 → 56**. `lod0_triangles` and
//! `coarse_triangles` carry them so a reader can check the ladder is the same
//! one rather than trusting this comment.
//!
//! # What the instruments are
//!
//! **The predictor (C1's subject).** One field sample, one gradient, one Newton
//! step: `p = v − f(v)·∇f(v)/|∇f(v)|²`. No fine mesh is read.
//!
//! **The reference (C1's ground truth), and it is two references, not one.**
//! `M-289` is the precedent for a reference that was wrong exactly where the
//! measurement was taken, so this harness carries two independent ones and
//! checks both against a closed form before believing either:
//!
//! 1. *Mesh reference* — the exact nearest point on the **LOD-0 triangle mesh**,
//!    closed-form point-vs-triangle over a uniform grid of the fine triangles.
//!    This is literally the registration's "nearest fine-surface point": the
//!    fine surface is what LOD 0 renders. It populates the registered columns
//!    `residual_agree_fraction` and `residual_p95_cells`.
//! 2. *Analytic reference* — the nearest point on the field's **analytic zero
//!    set**, found by a fan of 256 Fibonacci rays plus three rounds of local
//!    angular refinement, each ray marched at `h₀/8` and bisected 50 times. It
//!    touches the gradient nowhere, which is what makes it a valid reference for
//!    a gradient-based predictor. Run on a stride-subsample and reported as
//!    `residual_agree_fraction_analytic` / `residual_p95_analytic_cells`.
//!
//! Both are checked against a **closed form** on the four fields that have one —
//! `sphere`, `torus`, `box_exact`, `thin_plate` — in `closed_form_ref_err_cells`
//! and `closed_form_analytic_err_cells`, and the analytic one is *asserted*
//! against it. `-1` in those columns means the field has no closed form here.
//!
//! **The shading comparison (C2 and C3).** There is no Bevy renderer in this
//! crate (`CLAUDE.md` hard rule 2) and none is faked. What a renderer does at a
//! pixel is: take the triangle covering it, barycentrically interpolate that
//! triangle's vertex normals, and light with the result. So the instrument is
//! that same arithmetic, evaluated on the coarse surface at 16 area-uniform
//! barycentric points per coarse triangle and **area-weighted** — the
//! view-independent form of pixel coverage, stated because no camera exists to
//! supply `|N·V|`. Three normals per sample point:
//!
//! - **LOD-0 shading (the reference):** the LOD-0 mesh's own barycentrically
//!   interpolated vertex normal at the LOD-0 point nearest the sample. Those
//!   vertex normals are the crate's, not this harness's.
//! - **LOD-2 without the map:** the coarse mesh's interpolated vertex normal at
//!   the sample. Again the crate's normals.
//! - **LOD-2 with the map:** a bilinear fetch from the baked map, decoded
//!   through the coarse tangent frame.
//!
//! **Sampling at coarse mesh *vertices* would have measured almost nothing, and
//! that is a fixture defect this harness caught before running.** At a vertex
//! the interpolated normal *is* the vertex normal, which is
//! `normalize(∇f(v))` — the exact analytic normal, all octaves included. On
//! `fbm_terrain` the field is `y − h(x, z)`, exactly linear along a `y` edge, so
//! a coarse vertex on a vertical edge sits exactly on the analytic surface and
//! its normal is exactly right. A vertex-sampled comparison would therefore have
//! reported ≈0° for both arms and "falsified" C2 by measuring the one place
//! where coarsening costs nothing. The loss lives in the **interpolation across
//! a coarse triangle**, which is what an area sample sees and a vertex sample
//! cannot. `changed_coarse_mesh_vertices`, `angle_no_map_at_vertices_deg` and
//! `angle_with_map_at_vertices_deg` report the vertex-only versions beside the
//! registered ones, so the size of that trap is on the row rather than in this
//! comment.
//!
//! # Four attribution arms, because "the map failed" names three different things
//!
//! C2 can miss for reasons that want different answers, so the with-map arm is
//! bracketed:
//!
//! - `angle_direct_deg` — the residual normal evaluated **at the shading point
//!   with no texture at all**: no projection, no quantisation, no tangent frame,
//!   no bilinear filter. The difference from `angle_with_map_deg` is the map's
//!   own cost.
//! - `angle_deep_deg` — the same with `DEEP_STEPS` Newton steps instead of one,
//!   so "one root refinement was not enough" is a measurement rather than a
//!   suspicion.
//! - `lod0_self_angle_deg` — the angle between the **analytic** normal at the
//!   LOD-0 nearest point and the **LOD-0 mesh's own interpolated normal** there.
//!   This is a *floor* under the with-map arm and it is the reachability
//!   question `✗51`'s rule asks: a map baked from the analytic field cannot get
//!   closer to LOD-0 shading than LOD-0 shading is to the analytic field. C1
//!   does not bound it — C1 bounds a *position* to 0.1 cells, and the normal
//!   error that position error implies is `0.1·h·κ`, which on a field whose
//!   curvature radius is under a coarse cell is a large angle regardless.
//! - `angle_with_map_vs_analytic_deg` — the with-map normal against that
//!   analytic normal instead, which separates "the residual is wrong" from "the
//!   reference is smoother than the field it came from".
//!
//! # The map
//!
//! Triplanar, because the crate has no UV parametrisation and adding one would
//! be a source change: three axis-aligned pages of `64 × 64` texels, two `u8`
//! channels each, `3 · 64² · 2 = 24,576` bytes per chunk. 64 texels across the
//! chunk is one texel per LOD-0 cell, which is the Nyquist rate for the detail
//! being restored.
//!
//! A texel is baked in the page whose axis is **dominant** for the coarse normal
//! there, and fetched from the page dominant for the fetch point's coarse
//! normal — one rule at both ends, which is also what keeps the tangent frame
//! well conditioned: with `u` an in-page axis and the page axis dominant,
//! `|u − (u·N)N| ≥ √(2/3) = 0.816` always, so no degenerate-frame branch exists.
//! Texels no dominant-page triangle covered are filled by a one-pass dilation
//! (the standard bake gutter). An untouched texel encodes `(0, 0)`, which decodes
//! to the coarse normal exactly — "no detail here" needs no fallback path.
//!
//! # Controls, each an `assert!` rather than a printed number
//!
//! - **`changed_vertices` > 0** — the registered vacuity control, counted over
//!   the population C2's means are taken over, because that is the population it
//!   exists to guard.
//! - **`map_perturbation_deg` > 0** — the map must actually move the shading
//!   normal. A bug that decoded every texel to the identity would leave both C2
//!   arms equal and score a false HELD wherever the coarse normal happened to be
//!   close; this catches it.
//! - **The analytic reference against the closed form**, asserted under
//!   `0.01 h₀` on the four fields that have one. This is the `M-289` guard: a
//!   wrong ring bound in the nearest-triangle search or a wrong bisection
//!   bracket moves it by orders of magnitude.
//! - **Both meshes non-empty, both sample populations non-empty.**

#![allow(
    clippy::cast_possible_truncation,
    clippy::cast_precision_loss,
    clippy::cast_sign_loss,
    clippy::cast_possible_wrap,
    clippy::needless_range_loop,
    clippy::too_many_arguments,
    clippy::too_many_lines,
    clippy::print_literal
)]

mod common;

use std::cell::Cell;
use std::time::Instant;

use isomesh::chunk::ChunkLayout;
use isomesh::fields::{BoxExact, ReferenceField, Sphere, ThinPlate, Torus};
use isomesh::marching_cubes::MarchingCubes;
use isomesh::{MeshBuffer, Sdf};

// ─── knobs ──────────────────────────────────────────────────────────────────

/// Cells per axis at LOD 0, over each field's whole domain.
///
/// 64 because that is `M-72`'s ladder and `ThinPlate::CANONICAL_CELL_SIZE`'s
/// grid: the plate is built for `2·COMPACT_DOMAIN / 64`, so any other number
/// would measure a differently-aliased plate and stop being comparable with the
/// finding this experiment is built on.
const LOD0_CELLS: u32 = 64;

/// Levels measured. C1 is registered at 1 and 2, C3 at 3.
const LEVELS: [u32; 3] = [1, 2, 3];

/// Texels per axis per triplanar page.
///
/// Equal to [`LOD0_CELLS`], i.e. one texel per LOD-0 cell: the detail being
/// restored is LOD-0 detail, and a map coarser than the detail cannot carry it
/// while a finer one carries noise the fine mesh does not have either.
const MAP_TEXELS: usize = 64;

/// C1's tolerance, in cells **at the level being measured**.
const AGREE_CELLS: f64 = 0.1;

/// C2's ceiling with the map, and C2/C3's floor without, in degrees.
const ANGLE_WITH_MAP_DEG: f64 = 10.0;
const ANGLE_NO_MAP_DEG: f64 = 25.0;

/// Rays in the analytic reference's first-stage fan.
///
/// A locally planar surface at distance `d` in direction `u` is hit at `d/cos θ`
/// by a ray `θ` off it, so `N` Fibonacci directions cost `θ²/2` in relative
/// distance with `θ ≈ ½√(4π/N)`: **0.62%** at 256, which is 0.003 cells on a
/// 0.5-cell residual — thirty times inside C1's 0.1. The refinement below drives
/// it further down; this is the bound that makes the clause reachable at all.
const FAN_DIRS: usize = 256;

/// Local refinement rounds around the best direction, and directions per round.
///
/// **Five, and the count is derived rather than picked.** What the closed-form
/// assertion below compares is a *position*, not a distance, and a residual
/// angular error `θ` displaces the found point tangentially by `d·θ` — so the
/// distance error `d·θ²/2` being negligible is not enough. Each round localises
/// the direction to about `cap/3`, and `cap` shrinks by `FAN_CAP_SHRINK`, so
/// after `n` rounds `θ ≈ 0.15·0.3ⁿ/3`: **4.05e-4 rad at n = 5**. A coarse vertex
/// sits on a straddling coarse edge, so `d ≤ h/2`, which is `4 h₀` at LOD 3 —
/// giving `d·θ ≈ 1.6e-3 h₀` against the `1e-2 h₀` the assertion allows. At three
/// rounds it would have been `1.8e-2 h₀` and the reference would have failed its
/// own check.
const FAN_ROUNDS: usize = 5;
const FAN_ROUND_DIRS: usize = 24;

/// Starting half-angle of the refinement cap, shrunk by `FAN_CAP_SHRINK` each
/// round. 0.15 rad covers the 0.111 rad worst-case Fibonacci gap at 256.
const FAN_CAP: f64 = 0.15;
const FAN_CAP_SHRINK: f64 = 0.3;

/// March step, as a fraction of the LOD-0 cell size.
///
/// A step coarser than the thinnest feature steps over it, and `thin_plate` is
/// 0.4 LOD-0 cells thick — so eight steps per fine cell puts three samples
/// inside the plate.
const MARCH_STEPS_PER_FINE_CELL: f64 = 8.0;

/// How far the fan looks, in cells at the level being measured.
const FAN_REACH_COARSE_CELLS: f64 = 3.0;

/// Bisection halvings once a bracket is found. 50 exhausts `f64` on any bracket
/// this harness can produce.
const BISECT_STEPS: u32 = 50;

/// Coarse vertices given the analytic reference, at most. Stride-subsampled.
const ANALYTIC_CAP: usize = 256;

/// Shading sample points per row, at most. Sub-samples per triangle drop from
/// 16 to 4 to 1 rather than the triangles being strided, because Marching Cubes
/// emits in `z`-major order and a strided triangle set is a set of `z` slabs —
/// a spatially biased sample of a spatially varying quantity.
const SHADING_CAP: usize = 200_000;

/// Bakes timed per row; median reported.
const BAKE_REPS: usize = 3;

/// Newton steps in the attribution arm, against the registered one step.
const DEEP_STEPS: u32 = 3;

// ─── small vector arithmetic ────────────────────────────────────────────────

type V3 = [f64; 3];

fn sub(a: V3, b: V3) -> V3 {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
fn add(a: V3, b: V3) -> V3 {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}
fn mul(a: V3, s: f64) -> V3 {
    [a[0] * s, a[1] * s, a[2] * s]
}
fn dot(a: V3, b: V3) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
fn cross(a: V3, b: V3) -> V3 {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
fn norm2(a: V3) -> f64 {
    dot(a, a)
}
fn length(a: V3) -> f64 {
    norm2(a).sqrt()
}

/// Unit vector, or `None` when there is no direction to return.
///
/// `None` rather than a substituted axis: a zero-length normal is an absence of
/// information and inventing one would put a second execution path under every
/// shading number in this file.
fn unit(a: V3) -> Option<V3> {
    let l = length(a);
    if l > 0.0 { Some(mul(a, 1.0 / l)) } else { None }
}

/// Angle between two unit vectors, in degrees.
fn angle_deg(a: V3, b: V3) -> f64 {
    dot(a, b).clamp(-1.0, 1.0).acos().to_degrees()
}

fn sign_of(x: f64) -> i8 {
    if x > 0.0 {
        1
    } else if x < 0.0 {
        -1
    } else {
        0
    }
}

/// The axis of largest magnitude. First-max, so ties are deterministic.
fn dominant_axis(n: V3) -> usize {
    let a = [n[0].abs(), n[1].abs(), n[2].abs()];
    if a[0] >= a[1] && a[0] >= a[2] {
        0
    } else if a[1] >= a[2] {
        1
    } else {
        2
    }
}

/// The two axes spanning the page whose normal axis is `axis`.
const PAGE_AXES: [[usize; 2]; 3] = [[1, 2], [2, 0], [0, 1]];

/// Tangent frame for a tangent-space map, from the shading normal and the page.
///
/// Well conditioned by construction **provided `axis` is `n`'s dominant axis**:
/// then `|n[u]| ≤ |n[axis]|` and `|n|= 1` force `|n[u]| ≤ 1/√3`, so the
/// projected `u` has length at least `√(2/3) = 0.816`. Bake and fetch both apply
/// that rule, which is why no ill-conditioned-frame branch exists here.
fn tangent_frame(n: V3, axis: usize) -> (V3, V3) {
    let ua = PAGE_AXES[axis][0];
    let mut u = [0.0; 3];
    u[ua] = 1.0;
    let t = unit(sub(u, mul(n, dot(u, n)))).expect("dominant-axis frame is never degenerate");
    (t, cross(n, t))
}

// ─── the predictor: one gradient, one root refinement, no fine mesh ─────────

/// `p = v − f(v)·∇f(v)/|∇f(v)|²`, exactly once.
///
/// Where `|∇f| = 0` the field names no direction to move in and the point is
/// returned unmoved — the same "no information" rule as [`unit`], counted by the
/// caller rather than hidden.
fn project_once<F: Sdf<Scalar = f64>>(field: &F, v: V3) -> V3 {
    let f = field.sample(v);
    let g = field.gradient(v);
    let g2 = norm2(g);
    if g2 > 0.0 {
        add(v, mul(g, -f / g2))
    } else {
        v
    }
}

/// `steps` Newton steps. The attribution arm for whether one step is the limit.
fn project_n<F: Sdf<Scalar = f64>>(field: &F, v: V3, steps: u32) -> V3 {
    let mut p = v;
    for _ in 0..steps {
        p = project_once(field, p);
    }
    p
}

// ─── reference 1: nearest point on the LOD-0 triangle mesh ─────────────────

/// Closest point on a triangle, with its barycentric coordinates.
///
/// Ericson, *Real-Time Collision Detection* (2005) §5.1.5 — the seven-region
/// form, so it is exact on edges and at vertices rather than solving a
/// least-squares system that is singular there.
fn closest_on_tri(p: V3, a: V3, b: V3, c: V3) -> (V3, [f64; 3]) {
    let ab = sub(b, a);
    let ac = sub(c, a);
    let ap = sub(p, a);
    let d1 = dot(ab, ap);
    let d2 = dot(ac, ap);
    if d1 <= 0.0 && d2 <= 0.0 {
        return (a, [1.0, 0.0, 0.0]);
    }
    let bp = sub(p, b);
    let d3 = dot(ab, bp);
    let d4 = dot(ac, bp);
    if d3 >= 0.0 && d4 <= d3 {
        return (b, [0.0, 1.0, 0.0]);
    }
    let vc = d1 * d4 - d3 * d2;
    if vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0 {
        let v = d1 / (d1 - d3);
        return (add(a, mul(ab, v)), [1.0 - v, v, 0.0]);
    }
    let cp = sub(p, c);
    let d5 = dot(ab, cp);
    let d6 = dot(ac, cp);
    if d6 >= 0.0 && d5 <= d6 {
        return (c, [0.0, 0.0, 1.0]);
    }
    let vb = d5 * d2 - d1 * d6;
    if vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0 {
        let w = d2 / (d2 - d6);
        return (add(a, mul(ac, w)), [1.0 - w, 0.0, w]);
    }
    let va = d3 * d6 - d5 * d4;
    if va <= 0.0 && (d4 - d3) >= 0.0 && (d5 - d6) >= 0.0 {
        let w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
        return (add(b, mul(sub(c, b), w)), [0.0, 1.0 - w, w]);
    }
    let denom = 1.0 / (va + vb + vc);
    let v = vb * denom;
    let w = vc * denom;
    (
        add(add(a, mul(ab, v)), mul(ac, w)),
        [1.0 - v - w, v, w],
    )
}

/// A uniform grid over a mesh's triangles, for exact nearest-point queries.
///
/// Degenerate triangles are dropped at build time and counted: `closest_on_tri`
/// divides by `va + vb + vc` in its interior branch, and that sum is zero for a
/// zero-area triangle.
struct MeshGrid {
    lo: V3,
    cell: f64,
    dims: [usize; 3],
    starts: Vec<u32>,
    items: Vec<u32>,
    degenerate: usize,
}

impl MeshGrid {
    fn new(positions: &[V3], indices: &[u32], cell: f64) -> Self {
        let mut lo = [f64::INFINITY; 3];
        let mut hi = [f64::NEG_INFINITY; 3];
        for p in positions {
            for k in 0..3 {
                lo[k] = lo[k].min(p[k]);
                hi[k] = hi[k].max(p[k]);
            }
        }
        for k in 0..3 {
            lo[k] -= cell;
            hi[k] += cell;
        }
        let dims = [0, 1, 2].map(|k| (((hi[k] - lo[k]) / cell).ceil() as usize).max(1));
        let ncells = dims[0] * dims[1] * dims[2];

        let tri_count = indices.len() / 3;
        let mut degenerate = 0;
        let mut boxes: Vec<(u32, [usize; 3], [usize; 3])> = Vec::with_capacity(tri_count);
        for t in 0..tri_count {
            let (a, b, c) = (
                positions[indices[t * 3] as usize],
                positions[indices[t * 3 + 1] as usize],
                positions[indices[t * 3 + 2] as usize],
            );
            if length(cross(sub(b, a), sub(c, a))) <= 0.0 {
                degenerate += 1;
                continue;
            }
            let mut tlo = [0usize; 3];
            let mut thi = [0usize; 3];
            for k in 0..3 {
                let mn = a[k].min(b[k]).min(c[k]);
                let mx = a[k].max(b[k]).max(c[k]);
                tlo[k] = (((mn - lo[k]) / cell).floor() as isize).clamp(0, dims[k] as isize - 1)
                    as usize;
                thi[k] = (((mx - lo[k]) / cell).floor() as isize).clamp(0, dims[k] as isize - 1)
                    as usize;
            }
            boxes.push((t as u32, tlo, thi));
        }

        let mut counts = vec![0u32; ncells + 1];
        let idx = |dims: [usize; 3], x: usize, y: usize, z: usize| {
            (z * dims[1] + y) * dims[0] + x
        };
        for (_, tlo, thi) in &boxes {
            for z in tlo[2]..=thi[2] {
                for y in tlo[1]..=thi[1] {
                    for x in tlo[0]..=thi[0] {
                        counts[idx(dims, x, y, z) + 1] += 1;
                    }
                }
            }
        }
        for i in 0..ncells {
            counts[i + 1] += counts[i];
        }
        let total = counts[ncells] as usize;
        let starts = counts.clone();
        let mut cursor = counts;
        let mut items = vec![0u32; total];
        for (t, tlo, thi) in &boxes {
            for z in tlo[2]..=thi[2] {
                for y in tlo[1]..=thi[1] {
                    for x in tlo[0]..=thi[0] {
                        let c = idx(dims, x, y, z);
                        items[cursor[c] as usize] = *t;
                        cursor[c] += 1;
                    }
                }
            }
        }

        Self {
            lo,
            cell,
            dims,
            starts,
            items,
            degenerate,
        }
    }

    /// Nearest point on the mesh: `(distance, triangle, barycentric, point)`.
    ///
    /// Chebyshev shells outward from the query cell. A shell at index `r ≥ 1`
    /// cannot hold a point closer than `(r − 1)·cell` — the `−1` because the
    /// query point may sit anywhere inside its own cell — so the search stops
    /// when that bound exceeds the best distance found, and cannot stop before
    /// visiting the shell the true nearest triangle is registered in.
    fn nearest(&self, positions: &[V3], indices: &[u32], p: V3) -> (f64, usize, [f64; 3], V3) {
        let base = [0, 1, 2].map(|k| {
            (((p[k] - self.lo[k]) / self.cell).floor() as isize)
                .clamp(0, self.dims[k] as isize - 1)
        });
        let max_r = self.dims[0].max(self.dims[1]).max(self.dims[2]) as isize;
        let mut best = f64::INFINITY;
        let mut best_tri = usize::MAX;
        let mut best_bary = [0.0; 3];
        let mut best_pt = [0.0; 3];

        let visit = |x: isize,
                     y: isize,
                     z: isize,
                     best: &mut f64,
                     bt: &mut usize,
                     bb: &mut [f64; 3],
                     bp: &mut V3| {
            if x < 0
                || y < 0
                || z < 0
                || x >= self.dims[0] as isize
                || y >= self.dims[1] as isize
                || z >= self.dims[2] as isize
            {
                return;
            }
            let c = ((z as usize) * self.dims[1] + y as usize) * self.dims[0] + x as usize;
            let (s, e) = (self.starts[c] as usize, self.starts[c + 1] as usize);
            for &t in &self.items[s..e] {
                let t = t as usize;
                let (a, b, cc) = (
                    positions[indices[t * 3] as usize],
                    positions[indices[t * 3 + 1] as usize],
                    positions[indices[t * 3 + 2] as usize],
                );
                let (q, bary) = closest_on_tri(p, a, b, cc);
                let d = norm2(sub(q, p));
                if d < *best {
                    *best = d;
                    *bt = t;
                    *bb = bary;
                    *bp = q;
                }
            }
        };

        for r in 0..=max_r {
            if r > 0 && best.is_finite() {
                let bound = (r - 1) as f64 * self.cell;
                if bound * bound > best {
                    break;
                }
            }
            for dz in -r..=r {
                for dy in -r..=r {
                    if dz.abs() == r || dy.abs() == r {
                        for dx in -r..=r {
                            visit(
                                base[0] + dx,
                                base[1] + dy,
                                base[2] + dz,
                                &mut best,
                                &mut best_tri,
                                &mut best_bary,
                                &mut best_pt,
                            );
                        }
                    } else if r == 0 {
                        visit(
                            base[0], base[1], base[2], &mut best, &mut best_tri, &mut best_bary,
                            &mut best_pt,
                        );
                    } else {
                        for dx in [-r, r] {
                            visit(
                                base[0] + dx,
                                base[1] + dy,
                                base[2] + dz,
                                &mut best,
                                &mut best_tri,
                                &mut best_bary,
                                &mut best_pt,
                            );
                        }
                    }
                }
            }
        }
        assert!(
            best.is_finite(),
            "the nearest-triangle search found nothing, which means the grid is empty"
        );
        (best.sqrt(), best_tri, best_bary, best_pt)
    }
}

// ─── reference 2: nearest point on the analytic zero set ───────────────────

fn fib_dirs(n: usize) -> Vec<V3> {
    let ga = std::f64::consts::PI * (3.0 - 5.0f64.sqrt());
    (0..n)
        .map(|i| {
            let z = 1.0 - 2.0 * (i as f64 + 0.5) / (n as f64);
            let r = (1.0 - z * z).max(0.0).sqrt();
            let th = ga * i as f64;
            [r * th.cos(), r * th.sin(), z]
        })
        .collect()
}

/// Any orthonormal pair spanning the plane perpendicular to `n`.
fn ortho_basis(n: V3) -> (V3, V3) {
    let a = [n[0].abs(), n[1].abs(), n[2].abs()];
    let k = if a[0] <= a[1] && a[0] <= a[2] {
        0
    } else if a[1] <= a[2] {
        1
    } else {
        2
    };
    let mut u = [0.0; 3];
    u[k] = 1.0;
    let e1 = unit(cross(n, u)).expect("n is unit and u is its least-aligned axis");
    (e1, cross(n, e1))
}

/// First sign change of `f` along `o + t·d`, refined by bisection.
///
/// `sign0` is the sign at `t = 0` and is never zero — the caller handles a point
/// already on the surface, where there is no ray to march.
fn ray_first_root<F: Sdf<Scalar = f64>>(
    field: &F,
    o: V3,
    d: V3,
    t_max: f64,
    step: f64,
    sign0: i8,
) -> Option<f64> {
    let n = (t_max / step).ceil() as usize;
    let mut t_prev = 0.0;
    for i in 1..=n {
        let t = ((i as f64) * step).min(t_max);
        let s = sign_of(field.sample(add(o, mul(d, t))));
        if s != sign0 {
            let (mut lo, mut hi) = (t_prev, t);
            for _ in 0..BISECT_STEPS {
                let m = f64::midpoint(lo, hi);
                if sign_of(field.sample(add(o, mul(d, m)))) == sign0 {
                    lo = m;
                } else {
                    hi = m;
                }
            }
            return Some(f64::midpoint(lo, hi));
        }
        t_prev = t;
    }
    None
}

/// Nearest point on the analytic zero set, by ray fan plus local refinement.
///
/// `None` means no crossing within `t_max` in any probed direction, which is a
/// reportable fact rather than a failure: the surface really is further away
/// than the fan looked.
fn analytic_nearest<F: Sdf<Scalar = f64>>(
    field: &F,
    v: V3,
    dirs: &[V3],
    t_max: f64,
    step: f64,
) -> Option<(V3, f64)> {
    let sign0 = sign_of(field.sample(v));
    if sign0 == 0 {
        return Some((v, 0.0));
    }
    let mut best: Option<(V3, f64)> = None;
    for &d in dirs {
        if let Some(t) = ray_first_root(field, v, d, t_max, step, sign0)
            && best.is_none_or(|(_, bt)| t < bt)
        {
            best = Some((d, t));
        }
    }
    let (mut bd, mut bt) = best?;
    let mut cap = FAN_CAP;
    for _ in 0..FAN_ROUNDS {
        // `bd` is held fixed for the whole round and the round's winner is
        // committed at the end. **The first version updated `bd` inside the
        // loop, and the closed-form assertion caught it**: `e1`/`e2` came from
        // the round's starting direction, so once `bd` moved, `side` stopped
        // being orthogonal to it, `cos·bd + sin·side` stopped being a unit
        // vector, and `t` stopped being a distance — the comparison `t < bt`
        // was then between parameters along vectors of different lengths. The
        // search stalled at about 0.011 rad instead of 4e-4, which on `sphere`
        // at LOD 2 read 0.0214 LOD-0 cells against the closed form.
        let (e1, e2) = ortho_basis(bd);
        let rings = 3;
        let azimuths = FAN_ROUND_DIRS / rings;
        let mut round = (bd, bt);
        for ring in 1..=rings {
            let theta = cap * (ring as f64) / (rings as f64);
            let (st, ct) = (theta.sin(), theta.cos());
            for k in 0..azimuths {
                let phi = std::f64::consts::TAU * (k as f64) / (azimuths as f64);
                let side = add(mul(e1, phi.cos()), mul(e2, phi.sin()));
                let d = unit(add(mul(bd, ct), mul(side, st)))
                    .expect("a unit direction tilted by less than a radian is never zero");
                if let Some(t) = ray_first_root(field, v, d, t_max, step, sign0)
                    && t < round.1
                {
                    round = (d, t);
                }
            }
        }
        (bd, bt) = round;
        cap *= FAN_CAP_SHRINK;
    }
    Some((add(v, mul(bd, bt)), bt))
}

// ─── reference 3: closed forms, to check the other two ─────────────────────

/// Nearest surface point in closed form, where one exists.
#[derive(Clone, Copy)]
enum ClosedForm {
    None,
    Ball { c: V3, r: f64 },
    Ring { c: V3, major: f64, minor: f64 },
    Slab { c: V3, h: V3 },
}

impl ClosedForm {
    fn label(self) -> &'static str {
        match self {
            Self::None => "none",
            Self::Ball { .. } => "sphere",
            Self::Ring { .. } => "torus",
            Self::Slab { .. } => "box",
        }
    }

    /// The closed-form nearest point, or `None` where the formula is undefined
    /// (a query on a sphere's centre, a torus's axis or its core circle).
    fn nearest(self, p: V3) -> Option<V3> {
        match self {
            Self::None => None,
            Self::Ball { c, r } => unit(sub(p, c)).map(|d| add(c, mul(d, r))),
            Self::Ring { c, major, minor } => {
                let d = sub(p, c);
                let s = (d[0] * d[0] + d[2] * d[2]).sqrt();
                if s <= 0.0 {
                    return None;
                }
                let core = [major * d[0] / s, 0.0, major * d[2] / s];
                let q = sub(d, core);
                unit(q).map(|qu| add(add(c, core), mul(qu, minor)))
            }
            Self::Slab { c, h } => {
                let l = sub(p, c);
                let outside = (0..3).any(|k| l[k].abs() > h[k]);
                if outside {
                    let q = [0, 1, 2].map(|k| l[k].clamp(-h[k], h[k]));
                    Some(add(c, q))
                } else {
                    // Interior: leave through the face whose slack is smallest.
                    let mut k = 0;
                    for j in 1..3 {
                        if l[j].abs() - h[j] > l[k].abs() - h[k] {
                            k = j;
                        }
                    }
                    let mut q = l;
                    q[k] = if l[k] >= 0.0 { h[k] } else { -h[k] };
                    Some(add(c, q))
                }
            }
        }
    }
}

/// The closed form for a reference field, taken from the field's own canonical
/// parameters so the two cannot drift apart.
fn closed_form_for(name: &str) -> ClosedForm {
    match name {
        "sphere" => {
            let f = Sphere::<f64>::canonical();
            ClosedForm::Ball {
                c: f.center,
                r: f.radius,
            }
        }
        "torus" => {
            let f = Torus::<f64>::canonical();
            ClosedForm::Ring {
                c: f.center,
                major: f.major,
                minor: f.minor,
            }
        }
        "box_exact" => {
            let f = BoxExact::<f64>::canonical();
            ClosedForm::Slab {
                c: f.center,
                h: f.half_extents,
            }
        }
        "thin_plate" => {
            let f = ThinPlate::<f64>::canonical();
            ClosedForm::Slab {
                c: f.center,
                h: f.half_extents,
            }
        }
        _ => ClosedForm::None,
    }
}

// ─── the map ────────────────────────────────────────────────────────────────

fn encode(x: f64) -> u8 {
    ((x.clamp(-1.0, 1.0) * 127.0).round() as i32 + 128) as u8
}
fn decode(u: u8) -> f64 {
    (f64::from(u) - 128.0) / 127.0
}

#[derive(Default)]
struct BakeStats {
    written: usize,
    dilated: usize,
    conflicts: usize,
    hemisphere_clamped: usize,
    zero_gradient: usize,
}

/// Three axis-aligned pages of tangent-space normals, `MAP_TEXELS²` texels each,
/// two `u8` channels per texel.
struct Triplanar {
    lo: V3,
    extent: f64,
    enc: [Vec<u8>; 3],
    facing: [Vec<f32>; 3],
}

impl Triplanar {
    /// An untouched map. `(0, 0)` decodes to the coarse normal exactly, so an
    /// unbaked texel means "no detail here" without a second code path.
    fn blank(lo: V3, extent: f64) -> Self {
        let identity = encode(0.0);
        Self {
            lo,
            extent,
            enc: [
                vec![identity; MAP_TEXELS * MAP_TEXELS * 2],
                vec![identity; MAP_TEXELS * MAP_TEXELS * 2],
                vec![identity; MAP_TEXELS * MAP_TEXELS * 2],
            ],
            facing: [
                vec![-1.0; MAP_TEXELS * MAP_TEXELS],
                vec![-1.0; MAP_TEXELS * MAP_TEXELS],
                vec![-1.0; MAP_TEXELS * MAP_TEXELS],
            ],
        }
    }

    fn bytes(&self) -> usize {
        self.enc.iter().map(Vec::len).sum()
    }

    /// Continuous texel coordinates of a world point in a page.
    fn texel_of(&self, p: V3, axis: usize) -> (f64, f64) {
        let [s, t] = PAGE_AXES[axis];
        let r = MAP_TEXELS as f64;
        (
            (p[s] - self.lo[s]) / self.extent * r,
            (p[t] - self.lo[t]) / self.extent * r,
        )
    }

    /// Rasterise the coarse mesh into the three pages and write the residual
    /// normal at every covered texel.
    fn bake<F: Sdf<Scalar = f64>>(&mut self, field: &F, coarse: &MeshBuffer<f64>) -> BakeStats {
        let mut st = BakeStats::default();
        for t in 0..coarse.indices.len() / 3 {
            let iv = [
                coarse.indices[t * 3] as usize,
                coarse.indices[t * 3 + 1] as usize,
                coarse.indices[t * 3 + 2] as usize,
            ];
            let pv = iv.map(|i| coarse.positions[i]);
            let nv = iv.map(|i| coarse.normals[i]);
            let face = cross(sub(pv[1], pv[0]), sub(pv[2], pv[0]));
            let Some(face) = unit(face) else { continue };

            for axis in 0..3 {
                let facing = face[axis].abs() as f32;
                let uv: [(f64, f64); 3] = [0, 1, 2].map(|k| self.texel_of(pv[k], axis));
                let area2 = (uv[1].0 - uv[0].0) * (uv[2].1 - uv[0].1)
                    - (uv[2].0 - uv[0].0) * (uv[1].1 - uv[0].1);
                if area2.abs() <= 1e-12 {
                    continue;
                }
                let lo_i = uv.iter().map(|q| q.0).fold(f64::INFINITY, f64::min).floor();
                let hi_i = uv
                    .iter()
                    .map(|q| q.0)
                    .fold(f64::NEG_INFINITY, f64::max)
                    .ceil();
                let lo_j = uv.iter().map(|q| q.1).fold(f64::INFINITY, f64::min).floor();
                let hi_j = uv
                    .iter()
                    .map(|q| q.1)
                    .fold(f64::NEG_INFINITY, f64::max)
                    .ceil();
                let i0 = (lo_i as isize).clamp(0, MAP_TEXELS as isize - 1);
                let i1 = (hi_i as isize).clamp(0, MAP_TEXELS as isize - 1);
                let j0 = (lo_j as isize).clamp(0, MAP_TEXELS as isize - 1);
                let j1 = (hi_j as isize).clamp(0, MAP_TEXELS as isize - 1);

                for j in j0..=j1 {
                    for i in i0..=i1 {
                        let px = i as f64 + 0.5;
                        let py = j as f64 + 0.5;
                        let b1 = ((px - uv[0].0) * (uv[2].1 - uv[0].1)
                            - (uv[2].0 - uv[0].0) * (py - uv[0].1))
                            / area2;
                        let b2 = ((uv[1].0 - uv[0].0) * (py - uv[0].1)
                            - (px - uv[0].0) * (uv[1].1 - uv[0].1))
                            / area2;
                        let b0 = 1.0 - b1 - b2;
                        if b0 < -1e-9 || b1 < -1e-9 || b2 < -1e-9 {
                            continue;
                        }
                        let bary = [b0, b1, b2];
                        let s3 = (0..3).fold([0.0; 3], |acc, k| add(acc, mul(pv[k], bary[k])));
                        let ni = (0..3).fold([0.0; 3], |acc, k| add(acc, mul(nv[k], bary[k])));
                        let Some(nc) = unit(ni) else { continue };
                        // One rule at both ends: a texel belongs to the page its
                        // coarse normal is dominant in.
                        if dominant_axis(nc) != axis {
                            continue;
                        }
                        let idx = (j as usize) * MAP_TEXELS + i as usize;
                        if facing <= self.facing[axis][idx] {
                            continue;
                        }
                        if self.facing[axis][idx] >= 0.0 {
                            st.conflicts += 1;
                        } else {
                            st.written += 1;
                        }

                        let p = project_once(field, s3);
                        let g = field.gradient(p);
                        let (mut x, mut y) = (0.0, 0.0);
                        match unit(g) {
                            None => st.zero_gradient += 1,
                            Some(nd) => {
                                let (tt, bb) = tangent_frame(nc, axis);
                                let (ex, ey, ez) = (dot(nd, tt), dot(nd, bb), dot(nd, nc));
                                if ez > 0.0 {
                                    x = ex;
                                    y = ey;
                                } else {
                                    // A two-channel tangent map cannot hold a
                                    // normal in the lower hemisphere. The closest
                                    // it can hold is the equator.
                                    st.hemisphere_clamped += 1;
                                    let l = (ex * ex + ey * ey).sqrt();
                                    if l > 0.0 {
                                        x = ex / l;
                                        y = ey / l;
                                    }
                                }
                            }
                        }
                        self.enc[axis][idx * 2] = encode(x);
                        self.enc[axis][idx * 2 + 1] = encode(y);
                        self.facing[axis][idx] = facing;
                    }
                }
            }
        }

        // The bake gutter: one dilation pass, so a fetch straddling the edge of
        // a page's coverage lands on real detail rather than on the identity.
        for axis in 0..3 {
            let mut fills: Vec<(usize, u8, u8, f32)> = Vec::new();
            for j in 0..MAP_TEXELS {
                for i in 0..MAP_TEXELS {
                    let idx = j * MAP_TEXELS + i;
                    if self.facing[axis][idx] >= 0.0 {
                        continue;
                    }
                    let mut best = -1.0f32;
                    let mut pick = (0u8, 0u8);
                    for (di, dj) in [(-1isize, 0isize), (1, 0), (0, -1), (0, 1)] {
                        let ni = i as isize + di;
                        let nj = j as isize + dj;
                        if ni < 0 || nj < 0 || ni >= MAP_TEXELS as isize || nj >= MAP_TEXELS as isize
                        {
                            continue;
                        }
                        let n = (nj as usize) * MAP_TEXELS + ni as usize;
                        if self.facing[axis][n] > best {
                            best = self.facing[axis][n];
                            pick = (self.enc[axis][n * 2], self.enc[axis][n * 2 + 1]);
                        }
                    }
                    if best >= 0.0 {
                        fills.push((idx, pick.0, pick.1, best));
                    }
                }
            }
            st.dilated += fills.len();
            for (idx, x, y, f) in fills {
                self.enc[axis][idx * 2] = x;
                self.enc[axis][idx * 2 + 1] = y;
                self.facing[axis][idx] = f;
            }
        }
        st
    }

    /// Bilinear fetch, decoded through the coarse tangent frame.
    ///
    /// Returns the perturbed shading normal and whether any of the four texels
    /// had ever been written.
    fn fetch(&self, p: V3, nc: V3) -> (V3, bool) {
        let axis = dominant_axis(nc);
        let (fs, ft) = self.texel_of(p, axis);
        let (fs, ft) = (fs - 0.5, ft - 0.5);
        let i0 = fs.floor();
        let j0 = ft.floor();
        let (u, v) = (fs - i0, ft - j0);
        let mut x = 0.0;
        let mut y = 0.0;
        let mut hit = false;
        for (di, dj, w) in [
            (0isize, 0isize, (1.0 - u) * (1.0 - v)),
            (1, 0, u * (1.0 - v)),
            (0, 1, (1.0 - u) * v),
            (1, 1, u * v),
        ] {
            let i = (i0 as isize + di).clamp(0, MAP_TEXELS as isize - 1) as usize;
            let j = (j0 as isize + dj).clamp(0, MAP_TEXELS as isize - 1) as usize;
            let idx = j * MAP_TEXELS + i;
            hit |= self.facing[axis][idx] >= 0.0;
            x += w * decode(self.enc[axis][idx * 2]);
            y += w * decode(self.enc[axis][idx * 2 + 1]);
        }
        let z = (1.0 - x * x - y * y).max(0.0).sqrt();
        let (tt, bb) = tangent_frame(nc, axis);
        let n = add(add(mul(tt, x), mul(bb, y)), mul(nc, z));
        (
            unit(n).expect("x² + y² + z² = 1 by construction, so the sum is never zero"),
            hit,
        )
    }
}

// ─── field-evaluation counting, for a machine-independent bake cost ─────────

struct Counted<'a, F> {
    inner: &'a F,
    samples: &'a Cell<u64>,
    gradients: &'a Cell<u64>,
}

impl<F: Sdf<Scalar = f64>> Sdf for Counted<'_, F> {
    type Scalar = f64;
    fn sample(&self, p: V3) -> f64 {
        self.samples.set(self.samples.get() + 1);
        self.inner.sample(p)
    }
    fn gradient(&self, p: V3) -> V3 {
        self.gradients.set(self.gradients.get() + 1);
        self.inner.gradient(p)
    }
}

// ─── the ladder ─────────────────────────────────────────────────────────────

fn mesh_at<F: Sdf<Scalar = f64>>(field: &F, cells: u32, h0: f64, lo: V3, level: u32) -> (MeshBuffer<f64>, f64) {
    let lod = ChunkLayout::<f64>::new(cells, h0, lo)
        .expect("valid layout")
        .at_lod(level)
        .expect("valid level");
    let shape = lod.sample_shape().expect("valid shape");
    let mut out = MeshBuffer::<f64>::new();
    MarchingCubes::<f64>::new()
        .extract(field, &shape, lo, lod.cell_size(), &mut out)
        .expect("extraction");
    (out, lod.cell_size())
}

/// The 16 area-uniform barycentric points of a 4-way triangle subdivision:
/// the centroids of the 10 upright and 6 inverted sub-triangles.
fn bary_pattern(sub: usize) -> Vec<[f64; 3]> {
    match sub {
        1 => vec![[1.0 / 3.0; 3]],
        4 => {
            let mut v = Vec::with_capacity(4);
            for i in 0..2 {
                for j in 0..2 - i {
                    let u = (3.0 * i as f64 + 1.0) / 6.0;
                    let w = (3.0 * j as f64 + 1.0) / 6.0;
                    v.push([u, w, 1.0 - u - w]);
                }
            }
            v.push([2.0 / 6.0, 2.0 / 6.0, 2.0 / 6.0]);
            v
        }
        _ => {
            let mut v = Vec::with_capacity(16);
            for i in 0..4 {
                for j in 0..4 - i {
                    let u = (3.0 * i as f64 + 1.0) / 12.0;
                    let w = (3.0 * j as f64 + 1.0) / 12.0;
                    v.push([u, w, 1.0 - u - w]);
                }
            }
            for i in 0..3 {
                for j in 0..3 - i {
                    let u = (3.0 * i as f64 + 2.0) / 12.0;
                    let w = (3.0 * j as f64 + 2.0) / 12.0;
                    v.push([u, w, 1.0 - u - w]);
                }
            }
            v
        }
    }
}

/// One shading sample on the coarse surface.
struct Shade {
    pos: V3,
    nc: V3,
    weight: f64,
}

fn percentile(sorted: &[f64], q: f64) -> f64 {
    if sorted.is_empty() {
        return f64::NAN;
    }
    let k = ((q * sorted.len() as f64).ceil() as usize).clamp(1, sorted.len()) - 1;
    sorted[k]
}

fn cpu_mhz() -> u64 {
    std::fs::read_to_string("/proc/cpuinfo")
        .ok()
        .and_then(|s| {
            s.lines()
                .find(|l| l.starts_with("cpu MHz"))
                .and_then(|l| l.split(':').nth(1))
                .and_then(|v| v.trim().parse::<f64>().ok())
        })
        .map_or(0, |m| m.round() as u64)
}

type Row = Vec<(&'static str, String)>;

fn run_field<F: ReferenceField<Scalar = f64>>(
    name: &'static str,
    field: &F,
    dirs: &[V3],
    rows: &mut Vec<Row>,
) {
    assert_eq!(name, F::NAME, "the macro's name must be the field's own");
    let (lo, hi) = field.domain();
    let extent = hi[0] - lo[0];
    let h0 = extent / f64::from(LOD0_CELLS);
    let cf = closed_form_for(name);

    let (fine, fine_h) = mesh_at(field, LOD0_CELLS, h0, lo, 0);
    assert!(
        fine.triangle_count() > 0,
        "{name}: the LOD-0 reference mesh is empty, so nothing can be measured against it"
    );
    let grid = MeshGrid::new(&fine.positions, &fine.indices, 2.0 * h0);
    let march = h0 / MARCH_STEPS_PER_FINE_CELL;
    println!(
        "{name}: lod0 h = {fine_h:.6}, {} vertices, {} triangles ({} degenerate dropped)",
        fine.vertex_count(),
        fine.triangle_count(),
        grid.degenerate
    );

    for level in LEVELS {
        let cells = LOD0_CELLS >> level;
        let (coarse, h) = mesh_at(field, cells, h0, lo, level);
        assert!(
            coarse.triangle_count() > 0,
            "{name} lod {level}: the coarse mesh is empty"
        );
        let fan_reach = FAN_REACH_COARSE_CELLS * h;

        // ── C1: the residual, against both references ──────────────────────
        let mut err_mesh: Vec<f64> = Vec::with_capacity(coarse.vertex_count());
        let mut err_an: Vec<f64> = Vec::new();
        let mut mesh_vs_an: Vec<f64> = Vec::new();
        let mut cf_ref = f64::NAN;
        let mut cf_an = f64::NAN;
        let mut cf_pred = f64::NAN;
        let mut zero_grad_vertices = 0usize;
        let mut an_unreached = 0usize;
        let stride = (coarse.vertex_count() / ANALYTIC_CAP).max(1);

        for (vi, &v) in coarse.positions.iter().enumerate() {
            let p = project_once(field, v);
            if norm2(field.gradient(v)) <= 0.0 {
                zero_grad_vertices += 1;
            }
            let (_, _, _, q_mesh) = grid.nearest(&fine.positions, &fine.indices, v);
            err_mesh.push(length(sub(p, q_mesh)) / h);

            if vi % stride == 0 {
                match analytic_nearest(field, v, dirs, fan_reach, march) {
                    None => an_unreached += 1,
                    Some((q_an, _)) => {
                        err_an.push(length(sub(p, q_an)) / h);
                        mesh_vs_an.push(length(sub(q_mesh, q_an)) / h);
                        if let Some(q_cf) = cf.nearest(v) {
                            let e_ref = length(sub(q_mesh, q_cf)) / h0;
                            let e_an = length(sub(q_an, q_cf)) / h0;
                            let e_pred = length(sub(p, q_cf)) / h0;
                            cf_ref = if cf_ref.is_nan() { e_ref } else { cf_ref.max(e_ref) };
                            cf_an = if cf_an.is_nan() { e_an } else { cf_an.max(e_an) };
                            cf_pred = if cf_pred.is_nan() {
                                e_pred
                            } else {
                                cf_pred.max(e_pred)
                            };
                        }
                    }
                }
            }
        }
        assert!(
            !err_mesh.is_empty(),
            "{name} lod {level}: no coarse vertices, so C1 has no population"
        );
        assert!(
            !err_an.is_empty(),
            "{name} lod {level}: the analytic reference reached nothing"
        );
        // The M-289 guard. The analytic reference never touches the gradient and
        // never touches the fine mesh, so on a field with a closed form it must
        // agree with it to the fan's own angular precision. A wrong bisection
        // bracket or a wrong ring bound in the triangle search moves this by
        // orders of magnitude, not by percent.
        if !cf_an.is_nan() {
            assert!(
                cf_an < 0.01,
                "{name} lod {level}: the analytic reference is {cf_an:.6} LOD-0 cells from the \
                 closed form, so it is not a reference"
            );
        }

        err_mesh.sort_unstable_by(f64::total_cmp);
        err_an.sort_unstable_by(f64::total_cmp);
        mesh_vs_an.sort_unstable_by(f64::total_cmp);
        let agree = err_mesh.iter().filter(|e| **e <= AGREE_CELLS).count() as f64
            / err_mesh.len() as f64;
        let agree_an =
            err_an.iter().filter(|e| **e <= AGREE_CELLS).count() as f64 / err_an.len() as f64;

        // ── the bake ───────────────────────────────────────────────────────
        let samples = Cell::new(0u64);
        let gradients = Cell::new(0u64);
        let counted = Counted {
            inner: field,
            samples: &samples,
            gradients: &gradients,
        };
        let mut times: Vec<f64> = Vec::with_capacity(BAKE_REPS);
        let mut map = Triplanar::blank(lo, extent);
        let mut bake = BakeStats::default();
        for rep in 0..BAKE_REPS {
            let mut m = Triplanar::blank(lo, extent);
            samples.set(0);
            gradients.set(0);
            let t0 = Instant::now();
            let st = m.bake(&counted, &coarse);
            times.push(t0.elapsed().as_secs_f64() * 1e3);
            if rep == BAKE_REPS - 1 {
                bake = st;
                map = m;
            }
        }
        times.sort_unstable_by(f64::total_cmp);
        let bake_ms = times[times.len() / 2];
        assert!(
            bake.written > 0,
            "{name} lod {level}: the bake wrote no texels, so the map is blank"
        );

        // ── the shading sample set ─────────────────────────────────────────
        let tris = coarse.indices.len() / 3;
        let subs = if tris * 16 <= SHADING_CAP {
            16
        } else if tris * 4 <= SHADING_CAP {
            4
        } else {
            1
        };
        let pattern = bary_pattern(subs);
        let mut shades: Vec<Shade> = Vec::with_capacity(tris * pattern.len());
        let mut degenerate_samples = 0usize;
        for t in 0..tris {
            let iv = [
                coarse.indices[t * 3] as usize,
                coarse.indices[t * 3 + 1] as usize,
                coarse.indices[t * 3 + 2] as usize,
            ];
            let pv = iv.map(|i| coarse.positions[i]);
            let nv = iv.map(|i| coarse.normals[i]);
            let area = 0.5 * length(cross(sub(pv[1], pv[0]), sub(pv[2], pv[0])));
            if area <= 0.0 {
                degenerate_samples += pattern.len();
                continue;
            }
            let w = area / pattern.len() as f64;
            for b in &pattern {
                let pos = (0..3).fold([0.0; 3], |acc, k| add(acc, mul(pv[k], b[k])));
                let ni = (0..3).fold([0.0; 3], |acc, k| add(acc, mul(nv[k], b[k])));
                match unit(ni) {
                    None => degenerate_samples += 1,
                    Some(nc) => shades.push(Shade { pos, nc, weight: w }),
                }
            }
        }
        assert!(
            !shades.is_empty(),
            "{name} lod {level}: no shading samples, so C2 has no population"
        );

        // ── the shadings, and the floor under the with-map arm ─────────────
        //
        // `a_self` is the angle between the **analytic** normal at the LOD-0
        // nearest point and the LOD-0 *mesh's own* interpolated normal there. It
        // is the reference's distance from the field, and it is the floor under
        // C2's with-map arm: a map baked from the analytic field cannot get
        // closer to LOD-0 shading than LOD-0 shading is to the analytic field.
        // `a_with_an` measures the with-map normal against that analytic normal
        // instead, so the two together say whether the residual is inaccurate or
        // whether the reference is smoother than the field it came from.
        let mut wsum = 0.0;
        let mut a_no = 0.0;
        let mut a_with = 0.0;
        let mut a_direct = 0.0;
        let mut a_deep = 0.0;
        let mut a_perturb = 0.0;
        let mut a_self = 0.0;
        let mut a_with_an = 0.0;
        let mut changed = 0usize;
        let mut fetch_misses = 0usize;
        let mut ref_zero = 0usize;
        for s in &shades {
            let (_, tri, bary, q) = grid.nearest(&fine.positions, &fine.indices, s.pos);
            let rn = (0..3).fold([0.0; 3], |acc, k| {
                add(
                    acc,
                    mul(fine.normals[fine.indices[tri * 3 + k] as usize], bary[k]),
                )
            });
            let Some(reference) = unit(rn) else {
                ref_zero += 1;
                continue;
            };
            let (with_map, hit) = map.fetch(s.pos, s.nc);
            if !hit {
                fetch_misses += 1;
            }
            let direct = unit(field.gradient(project_once(field, s.pos))).unwrap_or(s.nc);
            let deep = unit(field.gradient(project_n(field, s.pos, DEEP_STEPS))).unwrap_or(s.nc);
            let analytic_at_q = unit(field.gradient(q)).unwrap_or(reference);

            let d_no = angle_deg(reference, s.nc);
            wsum += s.weight;
            a_no += s.weight * d_no;
            a_with += s.weight * angle_deg(reference, with_map);
            a_direct += s.weight * angle_deg(reference, direct);
            a_deep += s.weight * angle_deg(reference, deep);
            a_perturb += s.weight * angle_deg(s.nc, with_map);
            a_self += s.weight * angle_deg(reference, analytic_at_q);
            a_with_an += s.weight * angle_deg(analytic_at_q, with_map);
            if d_no > ANGLE_NO_MAP_DEG {
                changed += 1;
            }
        }
        assert!(wsum > 0.0, "{name} lod {level}: zero shading weight");
        let (a_no, a_with, a_direct, a_deep, a_perturb, a_self, a_with_an) = (
            a_no / wsum,
            a_with / wsum,
            a_direct / wsum,
            a_deep / wsum,
            a_perturb / wsum,
            a_self / wsum,
            a_with_an / wsum,
        );

        let c2_reg = name == "fbm_terrain" && level == 2;
        let c3_reg = (name == "thin_plate" || name == "gyroid") && level == 3;

        // The registered vacuity control, counted over exactly the population
        // the two means above are taken over, because that is the population it
        // exists to guard: "or C2 is comparing two identical shadings".
        //
        // **Asserted on the rows a clause is claimed on, and reported on all of
        // them.** The first version asserted it on every row and stopped the run
        // on `sphere` at LOD 1, correctly: a unit sphere at `h = 0.125` is
        // resolved well enough that one coarsening step never moves a shading
        // normal by 25° anywhere. That zero is a *finding* about a smooth field,
        // not a broken fixture, and turning it into a halt would have made the
        // control decide which fields may appear in the dataset. What the
        // control has to stop is a vacuous *verdict*, so it stops one exactly
        // where a verdict is registered — C2 on `fbm_terrain` at LOD 2, C3 on
        // `thin_plate` and `gyroid` at LOD 3 — and `main` additionally asserts
        // the literal registered wording, that the two shadings differ by more
        // than 25° *somewhere*, straight off the emitted rows.
        if c2_reg || c3_reg {
            assert!(
                changed > 0,
                "{name} lod {level}: the LOD-0 reference never differs from LOD-{level} by more \
                 than {ANGLE_NO_MAP_DEG}° anywhere, so this row's registered verdict would be \
                 comparing two identical shadings"
            );
        }
        // And the map must actually move the normal, or a decode bug would score
        // a false HELD wherever the coarse normal was already close.
        assert!(
            a_perturb > 0.0,
            "{name} lod {level}: the map perturbs nothing, so the with-map arm is the no-map arm"
        );

        // Beside the registered control: the same comparison taken at coarse
        // mesh **vertices** only, count and both means. This is the fixture the
        // module header rejects, measured rather than argued: at a vertex the
        // interpolated normal is the vertex normal, which is the exact analytic
        // gradient, so the loss coarsening actually causes is invisible there.
        let mut changed_mesh_vertices = 0usize;
        let mut v_no = 0.0;
        let mut v_with = 0.0;
        let mut v_n = 0usize;
        for (vi, &v) in coarse.positions.iter().enumerate() {
            let (_, tri, bary, _) = grid.nearest(&fine.positions, &fine.indices, v);
            let rn = (0..3).fold([0.0; 3], |acc, k| {
                add(
                    acc,
                    mul(fine.normals[fine.indices[tri * 3 + k] as usize], bary[k]),
                )
            });
            if let Some(reference) = unit(rn) {
                let d = angle_deg(reference, coarse.normals[vi]);
                if d > ANGLE_NO_MAP_DEG {
                    changed_mesh_vertices += 1;
                }
                let (with_map, _) = map.fetch(v, coarse.normals[vi]);
                v_no += d;
                v_with += angle_deg(reference, with_map);
                v_n += 1;
            }
        }
        assert!(v_n > 0, "{name} lod {level}: no coarse vertex had a reference");
        let (v_no, v_with) = (v_no / v_n as f64, v_with / v_n as f64);

        let c1 = agree >= 0.95;
        let c2 = a_with < ANGLE_WITH_MAP_DEG && a_no > ANGLE_NO_MAP_DEG;
        let c3 = a_with > ANGLE_NO_MAP_DEG;
        let c1_reg = level == 1 || level == 2;

        println!(
            "  lod {level} h={h:.5}: {} verts, {} tris | C1 agree {:.4} p95 {:.4} cells \
             (analytic {:.4}/{:.4}) | angles no-map {a_no:.2}° with-map {a_with:.2}° \
             direct {a_direct:.2}° | changed {changed}/{} | bake {bake_ms:.3} ms",
            coarse.vertex_count(),
            coarse.triangle_count(),
            agree,
            percentile(&err_mesh, 0.95),
            agree_an,
            percentile(&err_an, 0.95),
            shades.len(),
        );

        rows.push(vec![
            ("field", name.to_string()),
            ("lod_level", level.to_string()),
            ("coarse_vertices", coarse.vertex_count().to_string()),
            ("residual_agree_fraction", format!("{agree:.6}")),
            (
                "residual_p95_cells",
                format!("{:.6}", percentile(&err_mesh, 0.95)),
            ),
            ("angle_lod0_vs_lod2_no_map", format!("{a_no:.4}")),
            ("angle_lod0_vs_lod2_with_map", format!("{a_with:.4}")),
            ("changed_vertices", changed.to_string()),
            ("bake_ms", format!("{bake_ms:.4}")),
            ("map_bytes_per_chunk", map.bytes().to_string()),
            ("c1_holds", c1.to_string()),
            ("c2_holds", c2.to_string()),
            ("c3_holds", c3.to_string()),
            // ── extras: the same two angles under names that cannot be
            // misread on a row whose level is not 2 (C5's lesson).
            ("angle_no_map_deg", format!("{a_no:.4}")),
            ("angle_with_map_deg", format!("{a_with:.4}")),
            ("angle_direct_deg", format!("{a_direct:.4}")),
            ("angle_deep_deg", format!("{a_deep:.4}")),
            ("map_perturbation_deg", format!("{a_perturb:.4}")),
            // The floor under the with-map arm, and the same arm scored against
            // the analytic normal instead of against the LOD-0 mesh's own.
            ("lod0_self_angle_deg", format!("{a_self:.4}")),
            ("angle_with_map_vs_analytic_deg", format!("{a_with_an:.4}")),
            // The fixture the header rejects, measured: vertices only.
            ("angle_no_map_at_vertices_deg", format!("{v_no:.4}")),
            ("angle_with_map_at_vertices_deg", format!("{v_with:.4}")),
            // ── which clause this row is the registered evidence for
            ("c1_registered_row", c1_reg.to_string()),
            ("c2_registered_row", c2_reg.to_string()),
            ("c3_registered_row", c3_reg.to_string()),
            // ── C1, in full
            ("residual_mean_cells", {
                let m = err_mesh.iter().sum::<f64>() / err_mesh.len() as f64;
                format!("{m:.6}")
            }),
            (
                "residual_max_cells",
                format!("{:.6}", err_mesh.last().copied().unwrap_or(f64::NAN)),
            ),
            (
                "residual_p95_world",
                format!("{:.8}", percentile(&err_mesh, 0.95) * h),
            ),
            (
                "residual_agree_fraction_analytic",
                format!("{agree_an:.6}"),
            ),
            (
                "residual_p95_analytic_cells",
                format!("{:.6}", percentile(&err_an, 0.95)),
            ),
            ("analytic_vertices", err_an.len().to_string()),
            ("analytic_unreached", an_unreached.to_string()),
            (
                "mesh_vs_analytic_p95_cells",
                format!("{:.6}", percentile(&mesh_vs_an, 0.95)),
            ),
            ("zero_gradient_vertices", zero_grad_vertices.to_string()),
            // ── the references' own validation
            ("closed_form", cf.label().to_string()),
            (
                "closed_form_ref_err_cells",
                format!("{:.6}", if cf_ref.is_nan() { -1.0 } else { cf_ref }),
            ),
            (
                "closed_form_analytic_err_cells",
                format!("{:.6}", if cf_an.is_nan() { -1.0 } else { cf_an }),
            ),
            (
                "closed_form_pred_err_cells",
                format!("{:.6}", if cf_pred.is_nan() { -1.0 } else { cf_pred }),
            ),
            // ── the ladder, so a reader can check it is M-72's
            ("cell_size", format!("{h:.8}")),
            ("lod0_cell_size", format!("{h0:.8}")),
            ("coarse_triangles", coarse.triangle_count().to_string()),
            ("lod0_vertices", fine.vertex_count().to_string()),
            ("lod0_triangles", fine.triangle_count().to_string()),
            (
                "lod0_mesh_bytes",
                (fine.vertex_count() * 24 + fine.indices.len() * 4).to_string(),
            ),
            // ── the shading population
            ("shading_samples", shades.len().to_string()),
            ("sub_samples_per_triangle", subs.to_string()),
            ("degenerate_samples", degenerate_samples.to_string()),
            ("reference_normal_missing", ref_zero.to_string()),
            (
                "changed_fraction",
                format!("{:.6}", changed as f64 / shades.len() as f64),
            ),
            (
                "changed_coarse_mesh_vertices",
                changed_mesh_vertices.to_string(),
            ),
            // ── the map's own accounting
            ("map_texels_total", (3 * MAP_TEXELS * MAP_TEXELS).to_string()),
            ("map_texels_written", bake.written.to_string()),
            ("map_texels_dilated", bake.dilated.to_string()),
            ("map_texel_conflicts", bake.conflicts.to_string()),
            (
                "hemisphere_clamped_texels",
                bake.hemisphere_clamped.to_string(),
            ),
            ("zero_gradient_texels", bake.zero_gradient.to_string()),
            ("fetch_misses", fetch_misses.to_string()),
            ("bake_field_samples", samples.get().to_string()),
            ("bake_field_gradients", gradients.get().to_string()),
            (
                "bake_us_per_texel_written",
                format!("{:.4}", bake_ms * 1e3 / bake.written as f64),
            ),
            ("cpu_mhz", cpu_mhz().to_string()),
        ]);
    }
}

fn main() {
    // Cargo passes `--bench` under `cargo bench` and nothing under `cargo test`,
    // which is the discriminator that keeps `--all-targets` from running this in
    // a debug build and overwriting the committed CSV.
    if !std::env::args().any(|a| a == "--bench") {
        return;
    }

    let dirs = fib_dirs(FAN_DIRS);
    common::experiment::run(isomesh::experiment!("P-80"), |run| {
        let mut rows: Vec<Row> = Vec::new();
        isomesh::for_each_reference_field!(f64, |name, field| {
            run_field(name, &field, &dirs, &mut rows);
        });

        // The registered vacuity control in its literal form — "the LOD-0
        // reference must differ from LOD-2 by more than 25 degrees somewhere" —
        // read back off the values about to be written rather than off an
        // internal accumulator, so what is asserted is what the CSV says.
        let column = |row: &Row, key: &str| -> String {
            row.iter()
                .find(|(k, _)| *k == key)
                .map(|(_, v)| v.clone())
                .expect("every row carries every registered column")
        };
        let nonvacuous: Vec<String> = rows
            .iter()
            .filter(|r| {
                column(r, "changed_vertices")
                    .parse::<u64>()
                    .expect("a count")
                    > 0
            })
            .map(|r| format!("{}@{}", column(r, "field"), column(r, "lod_level")))
            .collect();
        println!(
            "\nvacuity control: {} of {} rows have changed_vertices > 0 — {}",
            nonvacuous.len(),
            rows.len(),
            nonvacuous.join(" ")
        );
        assert!(
            !nonvacuous.is_empty(),
            "changed_vertices is zero on every row, so LOD-0 and the coarse levels shade \
             identically everywhere and C2 has nothing to measure"
        );

        for row in rows {
            run.record(&row);
        }
    });
}
