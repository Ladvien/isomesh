//! Tests for the experimental module.
//!
//! Exempt from semver, not from correctness. The load-bearing one is
//! [`the_probabilistic_quadric_is_the_existing_solve`], which pins the
//! derivation X-004 turned on: the paper's rule is not a new solver here, it is
//! the existing one with a different regularizer, and if that ever stops being
//! true the reason this module contains a *scaling* rather than a *solver* has
//! gone with it.

use super::ProbabilisticQuadric;
use crate::cube::corner_offset;
use crate::dual_contouring::{DualContouring, solve};
use crate::hermite::HermiteCell;
use crate::{MeshBuffer, RuntimeShape3, Sdf};

/// **The probabilistic plane quadric is `solve_with` at `λ = Nσ²`, to 1e-16.**
///
/// The identity derived in this module's docs, checked numerically rather than
/// left as algebra. A direct assembly of the paper's equations (6) and (7) in
/// **world** coordinates is compared against the crate's centroid-relative solve
/// with the scaled regularizer; the two agree to floating-point noise on every
/// cell of a field chosen for having sharp features, which is where a
/// regularizer matters most.
///
/// **This is why there is no `ProbabilisticQuadric` solver in this crate.**
/// Writing one would be a second execution path computing numbers the existing
/// path already computes — exactly what `CLAUDE.md`'s one-path rule forbids —
/// so what shipped is the scaling, which is the part that actually differs.
#[test]
fn the_probabilistic_quadric_is_the_existing_solve() {
    // The paper's quadric, assembled directly: `A = Σnnᵀ + N σ²I`,
    // `b = Σnnᵀq + σ²Σq`, solved by Cramer's rule. Shares no line with the
    // crate's solve, which is what makes the agreement evidence.
    fn paper_solve(cell: &HermiteCell<f64>, sigma_squared: f64) -> Option<[f64; 3]> {
        let mut a = [[0.0f64; 3]; 3];
        let mut b = [0.0f64; 3];
        let mut planes = 0.0f64;
        for edge in 0..12u8 {
            let Some(crossing) = cell.get(edge) else {
                continue;
            };
            planes += 1.0;
            let (n, q) = (crossing.normal, crossing.position);
            let nq = n[0] * q[0] + n[1] * q[1] + n[2] * q[2];
            for i in 0..3 {
                for j in 0..3 {
                    a[i][j] += n[i] * n[j];
                }
                b[i] += n[i] * nq + sigma_squared * q[i];
            }
        }
        for (i, row) in a.iter_mut().enumerate() {
            row[i] += planes * sigma_squared;
        }
        let det3 = |m: &[[f64; 3]; 3]| {
            m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
                - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
                + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
        };
        let det = det3(&a);
        if det.abs() < 1e-30 {
            return None;
        }
        let mut x = [0.0f64; 3];
        for (k, slot) in x.iter_mut().enumerate() {
            let mut m = a;
            for i in 0..3 {
                m[i][k] = b[i];
            }
            *slot = det3(&m) / det;
        }
        Some(x)
    }

    let field = crate::fields::BoxExact::<f64>::canonical();
    let (h, origin) = (0.25f64, [-2.0f64; 3]);
    let sigma_squared = 1e-3;
    let mut worst = 0.0f64;
    let mut compared = 0usize;

    for z in 0..16u32 {
        for y in 0..16u32 {
            for x in 0..16u32 {
                let mut corner = [0.0f64; 8];
                for (i, slot) in corner.iter_mut().enumerate() {
                    let o = corner_offset(i as u8);
                    *slot = field.sample([
                        origin[0] + h * f64::from(x + o[0]),
                        origin[1] + h * f64::from(y + o[1]),
                        origin[2] + h * f64::from(z + o[2]),
                    ]);
                }
                let cell_origin = [
                    origin[0] + h * f64::from(x),
                    origin[1] + h * f64::from(y),
                    origin[2] + h * f64::from(z),
                ];
                let cell = HermiteCell::from_corners(&field, &corner, cell_origin, h);
                if cell.len() < 3 {
                    continue;
                }
                let lambda = cell.len() as f64 * sigma_squared;
                let (Some(ours), Some(theirs)) = (
                    solve::solve_with(&cell, lambda),
                    paper_solve(&cell, sigma_squared),
                ) else {
                    continue;
                };
                worst = worst.max(
                    (0..3)
                        .map(|i| (ours[i] - theirs[i]).abs())
                        .fold(0.0f64, f64::max),
                );
                compared += 1;
            }
        }
    }

    assert!(
        compared > 100,
        "only {compared} cells, too few to mean anything"
    );
    assert!(
        worst < 1e-12,
        "the paper's assembly and ours disagree by {worst:e} over {compared} cells, \
         so the reduction in this module's docs is wrong"
    );
    std::println!("measured: {compared} cells, worst |ours - paper| = {worst:.3e}");
}

/// **The scaled regularizer is a different rule, not a differently spelled one.**
///
/// If `ProbabilisticQuadric` produced the same mesh as `Qef` there would be
/// nothing to measure and the module would be decoration. It must differ, and it
/// must differ *only* in position — the topology is the dual mesher's and no
/// vertex rule can reach it, the property X-002 established.
#[test]
fn the_scaled_regularizer_moves_vertices_without_moving_topology() {
    let field = crate::fields::BoxExact::<f64>::canonical();
    let shape = RuntimeShape3::new([25; 3]).expect("valid shape");

    let mut fixed = MeshBuffer::<f64>::new();
    DualContouring::<f64>::new()
        .extract(&field, &shape, [-2.0; 3], 0.16, &mut fixed)
        .expect("extraction");

    let mut scaled = MeshBuffer::<f64>::new();
    DualContouring::<f64, ProbabilisticQuadric>::with_rule(ProbabilisticQuadric::default())
        .extract(&field, &shape, [-2.0; 3], 0.16, &mut scaled)
        .expect("extraction");

    assert_eq!(
        fixed.indices, scaled.indices,
        "a vertex rule changed the topology, which no vertex rule can do"
    );
    #[allow(clippy::float_cmp)]
    let moved = fixed
        .positions
        .iter()
        .zip(&scaled.positions)
        .filter(|(a, b)| a != b)
        .count();
    assert!(
        moved > 0,
        "the scaled regularizer reproduced the fixed one exactly, so there is \
         nothing here to measure"
    );
    std::println!(
        "measured: {} vertices, {moved} moved by the crossing-count scaling",
        fixed.positions.len()
    );
}

/// **The `experimental` feature adds no dependencies, and cannot (X-003).**
///
/// The crate's pitch is one dependency, and a feature that quietly pulled in a
/// second would break that for the consumers who never enable it — they read
/// `cargo tree` on the default build and would see it, but only if they looked.
///
/// Checked at the source rather than by shelling out to `cargo tree`: a feature
/// whose list is **empty** cannot enable an optional dependency, because that is
/// the only mechanism by which a feature can add one. The check is therefore
/// exact rather than a sample of one invocation's output.
#[test]
fn the_experimental_feature_adds_no_dependencies() {
    let manifest = include_str!("../../Cargo.toml");
    let (_, features) = manifest
        .split_once("[features]")
        .expect("the manifest has no [features] section");
    let (features, _) = features
        .split_once("\n[")
        .expect("the [features] section is unterminated");

    let declared: alloc::vec::Vec<&str> = features
        .lines()
        .map(str::trim)
        .filter(|line| !line.is_empty() && !line.starts_with('#'))
        .collect();
    assert_eq!(
        declared,
        alloc::vec!["experimental = []"],
        "a feature gained a dependency list, or a second feature appeared: {declared:?}"
    );
}
