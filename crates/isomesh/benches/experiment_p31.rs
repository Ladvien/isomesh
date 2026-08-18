//! **P-31 — is the homotopy certificate ever available in a dug scene?**
//!
//! Ticket: R-032. Pre-registered in the commit before this one.
//!
//! ```bash
//! cargo bench --bench experiment_p31
//! ```
//!
//! Writes `docs/experiments/p-31.csv`.
//!
//! # The design in one line
//!
//! The λ-medial homotopy guarantee needs `λ < wfs`. Measure `wfs` — the
//! minimum air-side `ρ` over discrete critical points — across 20 seeded
//! scenes carved by the crate's own [`BrushStack`] composition, and ask
//! whether it ever clears brush scale.
//!
//! # The instrument, and where its constants come from
//!
//! Critical-point candidates are interior air samples whose voxel-step
//! central-difference `‖∇ρ‖ < 0.5` — the dossier's own `θ > 120°` filter
//! constant (on a two-sheet medial, `‖∇ρ‖ = cos(θ/2)`), reused rather than
//! re-invented — then non-maximum-suppressed over 26-neighbourhoods: a sample
//! survives only when no neighbour reads a strictly smaller magnitude. `wfs`
//! is the minimum `ρ` over survivors, in voxels. The scan is restricted to
//! `|xᵢ| ≤ 1.5`, inside the solid block, so the exterior shell (gradient ≡ 1,
//! boundary partly outside the sampled domain) cannot contribute.
//!
//! # Reachability and inversion in one control
//!
//! A single 20-voxel cavity in generic position has a genuine critical point
//! at its centre with `ρ ≈ 20` voxels. The instrument must report the
//! certificate **available** there (`wfs ≥ 10` voxels, asserted) before it is
//! trusted reporting it absent on the dug scenes. Every dug scene must also
//! report a non-empty critical set and a sane carved fraction (2–60% of the
//! interior), or the minimum is not a measurement.
//!
//! # Counted, not timed
//!
//! Every registered column is a count or a distance in voxels — machine
//! independent. No timing A/B exists, so M-197's interleaving rule does not
//! apply; this note exists so its absence is not read as an oversight.

mod common;

use isomesh::Sdf;
use isomesh::brush::{Brush, BrushStack, Capsule};
use isomesh::normals::central_difference;

/// Samples per axis over `[-2, 2]`; `h = 2⁻⁴`.
const N: usize = 65;
const DOMAIN: f64 = 2.0;
const H: f64 = 2.0 * DOMAIN / 64.0;
/// The registered candidate threshold: the dossier's θ > 120° constant.
const EPSILON: f64 = 0.5;
/// Scan window half-width — inside the solid block, away from its faces.
const ROI: f64 = 1.5;
/// Dug scenes, and carve brushes per scene.
const SCENES: usize = 20;
const BRUSHES: usize = 12;
/// The registered wfs threshold, voxels.
const WFS_BAR: f64 = 2.0;

/// The solid block being dug: an exact box SDF, half-extent 1.9 — negative
/// (solid) across the whole scan window.
struct Block;

impl Sdf for Block {
    type Scalar = f64;

    fn sample(&self, p: [f64; 3]) -> f64 {
        let q = [p[0].abs() - 1.9, p[1].abs() - 1.9, p[2].abs() - 1.9];
        let outside =
            (q[0].max(0.0).powi(2) + q[1].max(0.0).powi(2) + q[2].max(0.0).powi(2)).sqrt();
        let inside = q[0].max(q[1]).max(q[2]).min(0.0);
        outside + inside
    }
}

/// Deterministic LCG, the p26 pattern.
struct Lcg(u64);

impl Lcg {
    fn next_u64(&mut self) -> u64 {
        self.0 = self
            .0
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        self.0
    }

    fn unit(&mut self) -> f64 {
        (((self.next_u64() >> 11) as f64) + 1.0) / (1u64 << 53) as f64
    }

    fn range(&mut self, lo: f64, hi: f64) -> f64 {
        lo + (hi - lo) * self.unit()
    }
}

/// The carve set for one dug scene: capsules (a third of them degenerate to
/// spheres, which `Capsule` handles by construction), all `Subtract`.
fn dug_brushes(scene: usize) -> Vec<Brush<Capsule<f64>>> {
    let mut lcg = Lcg(0x51ED_2701_u64.wrapping_mul(scene as u64 + 1));
    (0..BRUSHES)
        .map(|_| {
            let a = [
                lcg.range(-1.2, 1.2),
                lcg.range(-1.2, 1.2),
                lcg.range(-1.2, 1.2),
            ];
            let sphere = lcg.unit() < 1.0 / 3.0;
            let b = if sphere {
                a
            } else {
                [
                    lcg.range(-1.2, 1.2),
                    lcg.range(-1.2, 1.2),
                    lcg.range(-1.2, 1.2),
                ]
            };
            let radius = lcg.range(0.1, 0.35);
            Brush::subtract(Capsule { a, b, radius })
        })
        .collect()
}

struct SceneResult {
    air: u64,
    critical: u64,
    wfs_voxels: f64,
}

/// Scan one composed field: candidate magnitudes on the ROI lattice, NMS,
/// then the minimum ρ over survivors.
fn scan(field: &impl Sdf<Scalar = f64>) -> SceneResult {
    let lo = ((-ROI + DOMAIN) / H).round() as usize;
    let hi = ((ROI + DOMAIN) / H).round() as usize;
    let side = hi - lo + 1;
    let idx = |i: usize, j: usize, k: usize| (i * side + j) * side + k;
    let mut mag = vec![f64::INFINITY; side * side * side];
    let mut rho = vec![0.0f64; side * side * side];
    let mut air = 0u64;
    for i in 0..side {
        for j in 0..side {
            for k in 0..side {
                let p = [
                    -DOMAIN + ((lo + i) as f64) * H,
                    -DOMAIN + ((lo + j) as f64) * H,
                    -DOMAIN + ((lo + k) as f64) * H,
                ];
                let f = field.sample(p);
                if f <= 0.0 {
                    continue;
                }
                air += 1;
                let g = central_difference(field, p, H);
                mag[idx(i, j, k)] = (g[0] * g[0] + g[1] * g[1] + g[2] * g[2]).sqrt();
                rho[idx(i, j, k)] = f;
            }
        }
    }
    let mut critical = 0u64;
    let mut wfs = f64::INFINITY;
    for i in 1..side - 1 {
        for j in 1..side - 1 {
            for k in 1..side - 1 {
                let m = mag[idx(i, j, k)];
                if m >= EPSILON {
                    continue;
                }
                let mut minimal = true;
                'nms: for di in -1i64..=1 {
                    for dj in -1i64..=1 {
                        for dk in -1i64..=1 {
                            if di == 0 && dj == 0 && dk == 0 {
                                continue;
                            }
                            let n = idx(
                                (i as i64 + di) as usize,
                                (j as i64 + dj) as usize,
                                (k as i64 + dk) as usize,
                            );
                            if mag[n] < m {
                                minimal = false;
                                break 'nms;
                            }
                        }
                    }
                }
                if minimal {
                    critical += 1;
                    wfs = wfs.min(rho[idx(i, j, k)]);
                }
            }
        }
    }
    SceneResult {
        air,
        critical,
        wfs_voxels: wfs / H,
    }
}

fn main() {
    if !std::env::args().any(|arg| arg == "--bench") {
        return;
    }

    let prereg = isomesh::experiment!("P-31");
    common::experiment::run(prereg, |run| {
        let roi_samples = {
            let lo = ((-ROI + DOMAIN) / H).round() as usize;
            let hi = ((ROI + DOMAIN) / H).round() as usize;
            let side = hi - lo + 1;
            (side * side * side) as f64
        };

        let mut rows: Vec<(String, SceneResult)> = Vec::new();

        for scene in 0..SCENES {
            let brushes = dug_brushes(scene);
            let field = BrushStack {
                base: Block,
                brushes: &brushes,
            };
            let r = scan(&field);
            let air_frac = r.air as f64 / roi_samples;
            assert!(
                (0.02..=0.60).contains(&air_frac),
                "scene {scene}: carved fraction {air_frac:.3} outside [2%, 60%] — the fixture is \
                 not a dug scene"
            );
            assert!(
                r.critical > 0,
                "scene {scene}: zero critical points — a minimum over an empty set is not a \
                 measurement"
            );
            rows.push((format!("dug_{scene:02}"), r));
        }

        // The control: one 20-voxel cavity in generic position. The
        // certificate must be reported AVAILABLE here.
        let control_brushes = vec![Brush::subtract(Capsule {
            a: [0.0173, 0.0231, 0.0117],
            b: [0.0173, 0.0231, 0.0117],
            radius: 1.25,
        })];
        let field = BrushStack {
            base: Block,
            brushes: &control_brushes,
        };
        let control = scan(&field);
        assert!(
            control.critical > 0,
            "control cavity: zero critical points — the instrument is blind"
        );
        assert!(
            control.wfs_voxels >= 10.0,
            "control cavity reported wfs {:.1} voxels < 10 — an instrument that cannot say \
             AVAILABLE is not to be trusted saying absent",
            control.wfs_voxels
        );
        rows.push(("control_cavity".to_string(), control));

        println!(
            "{:>16} {:>9} {:>10} {:>10}",
            "scene", "air", "critical", "wfs (vx)"
        );
        let mut below = 0usize;
        for (name, r) in &rows {
            if name.starts_with("dug_") && r.wfs_voxels < WFS_BAR {
                below += 1;
            }
            println!(
                "{:>16} {:>9} {:>10} {:>10.2}",
                name, r.air, r.critical, r.wfs_voxels
            );
            run.record(&[
                ("scene", name.clone()),
                ("air_samples", r.air.to_string()),
                ("critical_points", r.critical.to_string()),
                ("wfs_voxels", format!("{:.3}", r.wfs_voxels)),
                ("epsilon", format!("{EPSILON}")),
                ("samples_per_axis", N.to_string()),
            ]);
        }

        let pct = 100.0 * below as f64 / SCENES as f64;
        println!();
        println!(
            "wfs < {WFS_BAR} voxels on {below}/{SCENES} dug scenes ({pct:.0}%) -- {} (H says > 80%; \
             falsified at ≥ 50% of scenes clearing the bar)",
            if pct > 80.0 {
                "HELD — the homotopy certificate is essentially never available, and the λ-medial \
                 line rests on Hausdorff stability"
            } else if pct <= 50.0 {
                "FALSIFIED — the certificate is comfortably available at brush scale"
            } else {
                "UNDECIDED, loudly — between the registered bounds"
            }
        );
        println!(
            "control: certificate AVAILABLE at {:.1} voxels — the absent verdicts above are \
             trustworthy",
            rows.last().map(|(_, r)| r.wfs_voxels).unwrap_or(f64::NAN)
        );
    });
}
