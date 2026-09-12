//! Six closed solids whose boundary genus is fixed by their construction.
//!
//! Ticket: `R-177`, registered as `P-182`. These are `P-140`'s bench-local
//! fixtures (`benches/experiment_p140.rs`) moved into the crate, so that a
//! validity harness has closed fields of genus **1, 2, 3 and 5** whose Euler
//! characteristic is *asserted* rather than recorded — until this module the
//! suite could assert `χ` on `sphere` and `torus` only, and `M-455` recorded
//! that `gyroid`'s topology gate was therefore *"record whatever the mesh
//! said"*.
//!
//! **The prescribed `χ` is arithmetic, not a reading.** Each constructor's doc
//! comment carries the derivation `P-140` used, and
//! [`expected_euler`](super::ReferenceField::expected_euler) returns `2 − 2g`
//! from the type's genus parameter. The measurement that checks it lives in
//! `fields/tests.rs` and shares no code with the derivation.
//!
//! # Construction A — a ball with `g` cylinders drilled straight through it
//!
//! `f(p) = max(|p| − R, −min_i(|(p_x, p_y) − c_i| − r))`: the ball of radius `R`
//! minus the union of `g` infinite `z`-parallel cylinders of radius `r` at
//! lateral centres `c_i`. Both operands are 1-Lipschitz, and `max`, `min` and
//! negation preserve that, so the field is 1-Lipschitz everywhere.
//!
//! Sound when `|c_i − c_j| > 2r` for every pair — the bores are disjoint — and
//! `max_i|c_i| + r < R`, which is what makes every bore exit through the ball's
//! top *and* bottom cap while leaving a shell of positive thickness at the
//! equator. Then the boundary is a sphere with `2g` open disks removed, glued
//! along `2g` circles to `g` open annuli, and Mayer–Vietoris is one line:
//!
//! ```text
//! chi = (2 − 2g) + g·chi(annulus) − 2g·chi(circle) = (2 − 2g) + 0 − 0 = 2 − 2g
//! ```
//!
//! Connected, closed and orientable, so the genus is `(2 − χ)/2 = g`.
//!
//! # Construction B — the closed `t`-neighbourhood of an embedded graph
//!
//! `f(p) = min_e dist(p, e) − t` over the graph's edges, each a closed segment.
//! Exact segment distance, so again 1-Lipschitz.
//!
//! Sound when the neighbourhood is a **regular** neighbourhood: every pair of
//! edges sharing no node is farther apart than `2t`, so their tubes are
//! disjoint, and every pair sharing one node `v` has its tubes meeting only near
//! `v` — two tubes of radius `t` about rays leaving `v` at angle `θ` intersect
//! exactly within `t / sin(θ/2)` of `v` along each ray, so holding that under
//! half the shorter edge keeps the two ends of any one edge from merging with
//! each other. `P-140` computed both numbers for each of the three graphs below
//! (`docs/experiments/p-140.csv`, columns `separation` and `merge_headroom`).
//!
//! A regular neighbourhood of a graph is a handlebody with `χ(N) = χ(G)`, and
//! `χ(∂M) = 2χ(M)` for a compact 3-manifold, so
//!
//! ```text
//! chi = 2·chi(G) = 2·(V − E)      genus = E − V + 1 = b1(G)
//! ```
//!
//! # Why the genus is a type parameter
//!
//! [`ReferenceField::NAME`](super::ReferenceField) is one constant per
//! type and is the key of the golden fixture, so six fields need six types.
//! `DrilledBall<R, 3>` and `ThickenedGraph<R, 5>` are those types, with the
//! genus where the compiler can read it: `expected_euler` is `Some(2 − 2·G)`
//! with no runtime field to drift, and a drilled ball's bore array has exactly
//! `G` entries by construction. As with every other reference field, `NAME`
//! names the **canonical** instance of the type — a `Sphere` of any radius is
//! `"sphere"`, and a `ThickenedGraph<R, 3>` built from any `b₁ = 3` graph is
//! `"graph_k4_g3"`.
//!
//! # What is not claimed
//!
//! A correct `χ` is not a manifoldness certificate. `P-140` measured ten
//! `(field, extractor)` pairs that reach the prescribed `χ` on a non-manifold
//! mesh, and `MeshReport::genus` stays whatever `validate_indexed` computes.

use alloc::vec::Vec;

use super::{BoundedSdf, COMPACT_DOMAIN, FieldBound, ReferenceField, cube_domain};
use crate::vec3::{dot, length, scale, sub};
use crate::{Real, Sdf};

/// Ball radius of every drilled ball: keeps the solid `0.45` clear of the
/// `[-2, 2]³` wall and leaves a `0.52` shell outside the outermost bore.
const BALL_RADIUS: f64 = 1.55;

/// Bore and tube radius, shared by both constructions: `4.5` cells at 33³.
const TUBE: f64 = 0.28;

/// Ring radius of the bores, and the cube graph's half-edge: three bores on a
/// ring of `0.75` are `1.299` apart, so a `0.739` wall survives between them.
const RHO: f64 = 0.75;

// ─── construction A ─────────────────────────────────────────────────────────

/// A ball with `G` `z`-parallel cylinders drilled straight through it.
///
/// Construction A of `P-140`; the module doc carries the Mayer–Vietoris
/// derivation of `χ = 2 − 2G`. The soundness conditions — bores pairwise more
/// than `2·bore` apart, every bore centre within `radius − bore` of the axis —
/// are the caller's to keep, as a [`Torus`](super::Torus) with `minor > major`
/// is the caller's problem; the three canonical constructors keep them.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DrilledBall<R: Real, const G: usize> {
    /// Ball radius.
    pub radius: R,
    /// Bore radius, shared by every bore.
    pub bore: R,
    /// Lateral `(x, y)` centre of each bore. Exactly `G` of them.
    pub bores: [[R; 2]; G],
}

impl<R: Real, const G: usize> DrilledBall<R, G> {
    /// The nearest bore's signed lateral distance, `|(p_x, p_y) − c_i| − bore`,
    /// with the lateral offset it was measured from.
    ///
    /// `INFINITY` with a zero offset when there are no bores, so `max(ball,
    /// −drilled)` is the plain ball.
    #[inline]
    fn nearest_bore(&self, p: [R; 3]) -> (R, [R; 2]) {
        let mut best = (R::INFINITY, [R::ZERO; 2]);
        for c in &self.bores {
            let d = [p[0] - c[0], p[1] - c[1]];
            let rho = (d[0] * d[0] + d[1] * d[1]).sqrt();
            let signed = rho - self.bore;
            if signed < best.0 {
                best = (signed, d);
            }
        }
        best
    }
}

impl<R: Real, const G: usize> Sdf for DrilledBall<R, G> {
    type Scalar = R;

    #[inline]
    fn sample(&self, p: [R; 3]) -> R {
        let ball = length(p) - self.radius;
        let (drilled, _) = self.nearest_bore(p);
        ball.max(-drilled)
    }

    /// Analytic: the active operand's unit gradient. The ball's is `p / |p|`;
    /// a bore's is minus the lateral unit vector from its axis. Both are
    /// undefined on a single line or point (the ball's centre, a bore's axis),
    /// where the zero vector is returned rather than `NaN`.
    #[inline]
    fn gradient(&self, p: [R; 3]) -> [R; 3] {
        let norm = length(p);
        let ball = norm - self.radius;
        let (drilled, d) = self.nearest_bore(p);
        if ball >= -drilled {
            if norm > R::ZERO {
                scale(p, R::ONE / norm)
            } else {
                [R::ZERO; 3]
            }
        } else {
            let rho = drilled + self.bore;
            if rho > R::ZERO {
                [-d[0] / rho, -d[1] / rho, R::ZERO]
            } else {
                [R::ZERO; 3]
            }
        }
    }
}

impl<R: Real, const G: usize> ReferenceField for DrilledBall<R, G> {
    /// The canonical instance's name: a drilled ball of genus 1, 2 or 3 is one
    /// of the three shipped reference fields, and any other bore count is the
    /// family name.
    const NAME: &'static str = match G {
        1 => "ball_drilled_g1",
        2 => "ball_drilled_g2",
        3 => "ball_drilled_g3",
        _ => "ball_drilled",
    };
    fn domain(&self) -> ([R; 3], [R; 3]) {
        cube_domain(COMPACT_DOMAIN)
    }
    fn closed_in_domain(&self) -> bool {
        true
    }
    /// `2 − 2G`, by Mayer–Vietoris (module doc). Asserted, never recorded.
    fn expected_euler(&self) -> Option<i64> {
        // `G` is a bore count, so it fits any `i64` that could hold a mesh.
        #[allow(clippy::cast_possible_wrap)]
        Some(2 - 2 * G as i64)
    }
    /// 1-Lipschitz: `|p| − R` and `|(p_x, p_y) − c| − r` are distances, and
    /// `max`, `min` and negation of 1-Lipschitz functions are 1-Lipschitz. Not
    /// `Exact`: inside the ball near a bore the value is the bore's lateral
    /// distance, which overstates the distance to the nearer of the two caps'
    /// rims, so `|f|` is not the distance to the surface everywhere.
    fn bound(&self) -> FieldBound {
        FieldBound::Lipschitz { l: 1.0 }
    }
}

impl<R: Real, const G: usize> BoundedSdf for DrilledBall<R, G> {
    fn value_bound(&self) -> FieldBound {
        FieldBound::Lipschitz { l: 1.0 }
    }
}

/// `G` points on the circle of radius `rho` in the `z = 0` plane, the first on
/// the `+x` axis.
///
/// **`libm::cos` and `libm::sin` by name, not `a.cos()`.** This crate is
/// `no_std`, but `lib.rs` links `std` under `#[cfg(test)]`, and on a concrete
/// `f64` the inherent `f64::cos` — the *platform's* libm — then shadows
/// [`Real::cos`]. The first version of this function called `a.cos()`; the golden
/// fixture was blessed under test on Linux/glibc, and macOS read 48 of 378 hashes
/// differing in the last ULP of `cos(2π/3)` and `sin(4π/3)` — every combination
/// of `graph_theta_g2` and `ball_drilled_g3`, the two fields with a three-point
/// ring (`✗130`). Naming the function keeps one float backend on every platform,
/// which is the whole reason `libm` is the dependency.
fn ring<R: Real, const G: usize>(rho: f64) -> [[R; 2]; G] {
    core::array::from_fn(|k| {
        // `k < G`, and `G` is 1, 2 or 3 here; the conversion is exact for any
        // bore count a solid could have.
        #[allow(clippy::cast_precision_loss)]
        let a = core::f64::consts::TAU * (k as f64) / (G as f64);
        [
            R::from_f64(rho * libm::cos(a)),
            R::from_f64(rho * libm::sin(a)),
        ]
    })
}

/// Construction A at bore count `G`, with `P-140`'s radii.
fn drilled<R: Real, const G: usize>() -> DrilledBall<R, G> {
    DrilledBall {
        radius: R::from_f64(BALL_RADIUS),
        bore: R::from_f64(TUBE),
        bores: ring(RHO),
    }
}

/// A ball of radius `1.55` with one axial bore of radius `0.28`: genus **1**,
/// `χ = 0`.
///
/// **Mayer–Vietoris.** A sphere with `2` disks removed (`χ = 2 − 2 = 0`) glued
/// along `2` circles (`χ = 0` each) to `1` annulus (`χ = 0`): `0 + 0 − 0 = 0`.
/// **Handlebody.** One handle on a ball is a solid torus, whose boundary is a
/// torus, `χ = 0`. Both give `2 − 2·1 = 0`.
#[must_use]
pub fn ball_drilled_g1<R: Real>() -> DrilledBall<R, 1> {
    DrilledBall {
        radius: R::from_f64(BALL_RADIUS),
        bore: R::from_f64(TUBE),
        bores: [[R::ZERO, R::ZERO]],
    }
}

/// A ball of radius `1.55` with two bores of radius `0.28` on a ring of radius
/// `0.75`: genus **2**, `χ = −2`.
///
/// **Mayer–Vietoris.** A sphere with `4` disks removed (`χ = −2`) glued along
/// `4` circles to `2` annuli: `−2 + 0 − 0 = −2`. **Handlebody.** Two handles on
/// a ball is the genus-2 handlebody, whose boundary is the genus-2 surface,
/// `χ = 2 − 2·2 = −2`. The bores are `1.5` apart, more than `2·0.28`.
#[must_use]
pub fn ball_drilled_g2<R: Real>() -> DrilledBall<R, 2> {
    drilled()
}

/// A ball of radius `1.55` with three bores of radius `0.28` on a ring of
/// radius `0.75`: genus **3**, `χ = −4`.
///
/// **Mayer–Vietoris.** A sphere with `6` disks removed (`χ = −4`) glued along
/// `6` circles to `3` annuli: `−4 + 0 − 0 = −4`. **Handlebody.** Three handles,
/// boundary of genus 3, `χ = 2 − 2·3 = −4`. The bores are `1.299` apart, so a
/// `0.739` wall survives between each pair.
#[must_use]
pub fn ball_drilled_g3<R: Real>() -> DrilledBall<R, 3> {
    drilled()
}

// ─── construction B ─────────────────────────────────────────────────────────

/// Distance from `p` to the closed segment `a → b`, and the point it is
/// measured to.
///
/// Exact, and the clamp is what makes a union of these a union of *capsules*
/// rather than of infinite cylinders.
#[inline]
fn point_segment<R: Real>(p: [R; 3], a: [R; 3], b: [R; 3]) -> (R, [R; 3]) {
    let ab = sub(b, a);
    let ap = sub(p, a);
    let len2 = dot(ab, ab);
    let t = if len2 > R::ZERO {
        (dot(ap, ab) / len2).clamp(R::ZERO, R::ONE)
    } else {
        R::ZERO
    };
    let offset = sub(ap, scale(ab, t));
    (length(offset), offset)
}

/// The closed `tube`-neighbourhood of a graph embedded in `R³`, whose first
/// Betti number is `G`.
///
/// Construction B of `P-140`; the module doc carries the graph-thickening
/// derivation of `χ = 2(V − E) = 2 − 2G`. The genus claim needs
/// `b₁(G) = E − V + 1 = G` for a connected graph and the neighbourhood to be
/// regular (module doc); the three canonical constructors satisfy both, and an
/// edge naming a node that does not exist is skipped rather than dereferenced.
#[derive(Clone, Debug, PartialEq)]
pub struct ThickenedGraph<R: Real, const G: usize> {
    /// Node positions.
    pub nodes: Vec<[R; 3]>,
    /// Edges, as index pairs into `nodes`.
    pub edges: Vec<[usize; 2]>,
    /// Tube radius.
    pub tube: R,
}

impl<R: Real, const G: usize> ThickenedGraph<R, G> {
    /// Distance to the nearest edge and the offset it was measured along.
    ///
    /// `INFINITY` with a zero offset when no edge is valid, so the field is
    /// positive everywhere and the solid is empty.
    #[inline]
    fn nearest_edge(&self, p: [R; 3]) -> (R, [R; 3]) {
        let mut best = (R::INFINITY, [R::ZERO; 3]);
        for e in &self.edges {
            if let (Some(a), Some(b)) = (self.nodes.get(e[0]), self.nodes.get(e[1])) {
                let candidate = point_segment(p, *a, *b);
                if candidate.0 < best.0 {
                    best = candidate;
                }
            }
        }
        best
    }
}

impl<R: Real, const G: usize> Sdf for ThickenedGraph<R, G> {
    type Scalar = R;

    #[inline]
    fn sample(&self, p: [R; 3]) -> R {
        self.nearest_edge(p).0 - self.tube
    }

    /// Analytic: the unit vector from the nearest point of the nearest edge to
    /// `p`. Undefined on the graph itself, where the zero vector is returned.
    #[inline]
    fn gradient(&self, p: [R; 3]) -> [R; 3] {
        let (d, offset) = self.nearest_edge(p);
        if d > R::ZERO && d < R::INFINITY {
            scale(offset, R::ONE / d)
        } else {
            [R::ZERO; 3]
        }
    }
}

impl<R: Real, const G: usize> ReferenceField for ThickenedGraph<R, G> {
    /// The canonical instance's name: the theta, `K4` and cube graphs are the
    /// three shipped reference fields of this family, at `b₁` 2, 3 and 5.
    const NAME: &'static str = match G {
        2 => "graph_theta_g2",
        3 => "graph_k4_g3",
        5 => "graph_cube_g5",
        _ => "thickened_graph",
    };
    fn domain(&self) -> ([R; 3], [R; 3]) {
        cube_domain(COMPACT_DOMAIN)
    }
    fn closed_in_domain(&self) -> bool {
        true
    }
    /// `2 − 2G`, by graph thickening (module doc). Asserted, never recorded.
    fn expected_euler(&self) -> Option<i64> {
        // `G` is a cycle count, so it fits any `i64` that could hold a mesh.
        #[allow(clippy::cast_possible_wrap)]
        Some(2 - 2 * G as i64)
    }
    /// 1-Lipschitz: a point-to-segment distance is a distance and `min`
    /// preserves the constant. Not `Exact`: where two tubes meet at a node the
    /// value is the distance to the nearer *segment*, which overstates the
    /// distance to the merged surface in the junction's re-entrant corner.
    fn bound(&self) -> FieldBound {
        FieldBound::Lipschitz { l: 1.0 }
    }
}

impl<R: Real, const G: usize> BoundedSdf for ThickenedGraph<R, G> {
    fn value_bound(&self) -> FieldBound {
        FieldBound::Lipschitz { l: 1.0 }
    }
}

/// Construction B on one graph, with `P-140`'s tube radius.
///
/// The debug assert is the `b₁ = E − V + 1` identity the genus parameter
/// states; the three graphs below are connected, so it is the whole claim.
fn graphed<R: Real, const G: usize>(
    nodes: Vec<[R; 3]>,
    edges: Vec<[usize; 2]>,
) -> ThickenedGraph<R, G> {
    debug_assert_eq!(
        edges.len() + 1,
        nodes.len() + G,
        "b1 = E - V + 1 must be the genus"
    );
    ThickenedGraph {
        nodes,
        edges,
        tube: R::from_f64(TUBE),
    }
}

/// Two poles joined by three two-segment arcs: `V = 5`, `E = 6`, `b₁ = 2`,
/// genus **2**, `χ = −2`.
///
/// **Graph thickening.** `χ = 2·χ(G) = 2·(V − E) = 2·(5 − 6) = −2`.
/// **Handlebody.** `genus = E − V + 1 = 6 − 5 + 1 = 2`, and `2 − 2·2 = −2`. The
/// two independent cycles are two of the three four-cycles through both poles;
/// the third is their sum. Poles at `z = ±1`, waypoints on a ring of radius `1`
/// in the `z = 0` plane, tube `0.28`.
#[must_use]
pub fn graph_theta_g2<R: Real>() -> ThickenedGraph<R, 2> {
    let pole = R::ONE;
    let mut nodes = alloc::vec![[R::ZERO, R::ZERO, pole], [R::ZERO, R::ZERO, -pole]];
    for w in ring::<R, 3>(1.0) {
        nodes.push([w[0], w[1], R::ZERO]);
    }
    graphed(
        nodes,
        alloc::vec![[0, 2], [0, 3], [0, 4], [1, 2], [1, 3], [1, 4]],
    )
}

/// The tetrahedron graph `K4` at half-diagonal `0.72`: `V = 4`, `E = 6`,
/// `b₁ = 3`, genus **3**, `χ = −4`.
///
/// **Graph thickening.** `χ = 2·(V − E) = 2·(4 − 6) = −4`. **Handlebody.**
/// `genus = 6 − 4 + 1 = 3`, and `2 − 2·3 = −4`. Three of the four triangles are
/// a cycle basis; the fourth is their sum. Tube `0.28`.
#[must_use]
pub fn graph_k4_g3<R: Real>() -> ThickenedGraph<R, 3> {
    let a = R::from_f64(0.72);
    graphed(
        alloc::vec![[a, a, a], [a, -a, -a], [-a, a, -a], [-a, -a, a]],
        alloc::vec![[0, 1], [0, 2], [0, 3], [1, 2], [1, 3], [2, 3]],
    )
}

/// The cube graph `Q3` at half-edge `0.75`: `V = 8`, `E = 12`, `b₁ = 5`, genus
/// **5**, `χ = −8`.
///
/// **Graph thickening.** `χ = 2·(V − E) = 2·(8 − 12) = −8`. **Handlebody.**
/// `genus = 12 − 8 + 1 = 5`, and `2 − 2·5 = −8`. Five of the six faces are a
/// cycle basis; the six four-cycles sum to zero over `GF(2)`. Tube `0.28`.
///
/// **Its adequacy is not monotone in resolution.** `P-140` measured Marching
/// Cubes on this field at χ **−8** at 7³ and 11³ and **nothing meshed** at 9³
/// (`docs/experiments/p-140.csv`): at 9³ over `[-2, 2]³` the cell size is
/// exactly `0.5`, so every tube axis at `±0.75` sits midway between the grid
/// planes at `±0.5` and `±1.0`, the nearest line of samples is `0.25·√2 = 0.354`
/// from every axis, and a tube of radius `0.28` contains no sample at all.
/// `fields/tests.rs` pins that row.
#[must_use]
pub fn graph_cube_g5<R: Real>() -> ThickenedGraph<R, 5> {
    let a = R::from_f64(RHO);
    let coord = |bit: u32, i: u32| if i & bit == 0 { -a } else { a };
    let nodes: Vec<[R; 3]> = (0..8u32)
        .map(|i| [coord(1, i), coord(2, i), coord(4, i)])
        .collect();

    let mut edges: Vec<[usize; 2]> = Vec::with_capacity(12);
    for i in 0..8u32 {
        for bit in [1u32, 2, 4] {
            let j = i ^ bit;
            if j > i {
                edges.push([i as usize, j as usize]);
            }
        }
    }
    graphed(nodes, edges)
}
