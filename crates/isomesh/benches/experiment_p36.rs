//! **P-36 — sparse-factor update vs refactor at brush scale: the routing decision.**
//!
//! Ticket: R-035b. Pre-registered in the commit before this one; the
//! substrate is M-333's, verified against a dense reference.
//!
//! ```bash
//! cargo bench --bench experiment_p36
//! ```
//!
//! Writes `docs/experiments/p-36.csv`.
//!
//! # The design in one line
//!
//! Perturb the gyroid chunk's surface inside a radius-4-voxel ball — a
//! *value* edit on a stable vertex-slot pattern, M-318's regime — rebuild the
//! heat operator's changed rows, and race a partial refactorization over the
//! elimination-tree ancestor closure against a full refactorization,
//! interleaved both orders, flops counted beside wall time.
//!
//! # The fixture contract, asserted — and what the first two runs taught
//!
//! The edited mesh must have the **same cell-keyed slot set and the same
//! index buffer** as the base mesh — this experiment is about value updates
//! on a stable pattern, and a pattern change aborts rather than measures.
//! Two field-level bump amplitudes (0.3 and 0.08 voxels) both aborted on
//! that contract before any verdict: near the zero set there are always
//! samples with |f| below any amplitude, so a *field* value edit flips grid
//! edges and changes the slot set — M-318's "appears, vanishes, or moves"
//! is the common case for real brushes, not the edge case. The perturbation
//! therefore acts at the **mesh level**: every vertex inside the radius-4
//! ball is displaced 0.1 voxels along its normal, clamped to stay inside its
//! own cell — changed values, provably stable pattern, which is exactly the
//! component of a brush edit a value-only factor update can serve. The
//! slot-set component needs symbolic repair on top, and the FINDINGS row
//! says so. Changed slots are the displaced vertices; changed operator rows
//! add their one-ring (cotan weights read the ring).
//!
//! # Validity, asserted before any verdict
//!
//! The updated factor must hold the same Frobenius bound as a fresh one; its
//! solve must agree with the refactored solve within 1e-8 relative; and a
//! deliberately skipped closure column must push the residual past the bound
//! — the inversion seen red first.

mod common;

use common::heat;
use isomesh::fields::capped_gyroid;
use isomesh::surface_nets::SurfaceNets;
use isomesh::{MeshBuffer, RuntimeShape3, Sdf};
use std::time::Instant;

const SAMPLES: u32 = 65;
const ORIGIN: [f64; 3] = [-2.0, -2.0, -2.0];
const GRID_H: f64 = 4.0 / 64.0;
/// The registered brush radius: 4 voxels. Displacement amplitude 0.1 voxels
/// along the vertex normal, clamped into the vertex's own cell — the
/// mesh-level edit that satisfies the registered stable-pattern contract
/// after two field-level amplitudes aborted on it (see the module docs).
const BRUSH_R: f64 = 4.0 * GRID_H;
const BRUSH_AMP: f64 = 0.1 * GRID_H;
/// Where the brush lands: the surface vertex nearest this probe.
const PROBE: [f64; 3] = [0.3, 0.2, -0.1];
const REPS: usize = 11;

fn extract(field: &impl Sdf<Scalar = f64>) -> MeshBuffer<f64> {
    let shape = RuntimeShape3::new([SAMPLES; 3]).expect("grid fits");
    let mut sn = SurfaceNets::<f64>::new();
    let mut out = MeshBuffer::<f64>::new();
    sn.extract(field, &shape, ORIGIN, GRID_H, &mut out)
        .expect("extraction");
    out
}

/// Cell identity of a Surface Nets vertex: the cell its position lies in.
fn cell_key(p: [f64; 3]) -> u64 {
    let c = |x: f64| (((x - ORIGIN[0]) / GRID_H).floor() as u64).min(63);
    (c(p[0]) << 32) | (c(p[1]) << 16) | c(p[2])
}

fn clone_factor(f: &heat::LdlFactor) -> heat::LdlFactor {
    heat::LdlFactor {
        n: f.n,
        parent: f.parent.clone(),
        li: f.li.clone(),
        lx: f.lx.clone(),
        d: f.d.clone(),
        flops: 0,
    }
}

fn main() {
    if !std::env::args().any(|arg| arg == "--bench") {
        return;
    }

    let prereg = isomesh::experiment!("P-36");
    common::experiment::run(prereg, |run| {
        let field = capped_gyroid::<f64>();
        let base_mesh = extract(&field);

        // Brush center: the surface vertex nearest the registered probe.
        let center = base_mesh
            .positions
            .iter()
            .min_by(|a, b| {
                let da = (a[0] - PROBE[0]).powi(2)
                    + (a[1] - PROBE[1]).powi(2)
                    + (a[2] - PROBE[2]).powi(2);
                let db = (b[0] - PROBE[0]).powi(2)
                    + (b[1] - PROBE[1]).powi(2)
                    + (b[2] - PROBE[2]).powi(2);
                da.total_cmp(&db)
            })
            .copied()
            .expect("mesh has vertices");
        // Mesh-level edit: displace in-ball vertices along their normals,
        // clamped inside their own cells so the slot keys cannot move.
        let mut edit_mesh = base_mesh.clone();
        for i in 0..edit_mesh.positions.len() {
            let p = edit_mesh.positions[i];
            let d2 = (p[0] - center[0]).powi(2)
                + (p[1] - center[1]).powi(2)
                + (p[2] - center[2]).powi(2);
            if d2 >= BRUSH_R * BRUSH_R {
                continue;
            }
            let n = edit_mesh.normals[i];
            let mut q = [
                p[0] + BRUSH_AMP * n[0],
                p[1] + BRUSH_AMP * n[1],
                p[2] + BRUSH_AMP * n[2],
            ];
            for axis in 0..3 {
                let cell = ((p[axis] - ORIGIN[axis]) / GRID_H).floor();
                let lo = ORIGIN[axis] + cell * GRID_H + 0.01 * GRID_H;
                let hi = ORIGIN[axis] + (cell + 1.0) * GRID_H - 0.01 * GRID_H;
                q[axis] = q[axis].clamp(lo, hi);
            }
            edit_mesh.positions[i] = q;
        }

        // ---- fixture contract: stable slot set, identical topology --------
        assert!(
            base_mesh.positions.len() == edit_mesh.positions.len(),
            "vertex count changed ({} -> {}) — the slot set moved and this \
             fixture is about value updates on a stable pattern",
            base_mesh.positions.len(),
            edit_mesh.positions.len()
        );
        assert!(
            base_mesh.indices == edit_mesh.indices,
            "index buffer changed — the slot set moved; fixture aborts"
        );
        let mut keys_base: Vec<u64> = base_mesh.positions.iter().map(|&p| cell_key(p)).collect();
        let keys_edit: Vec<u64> = edit_mesh.positions.iter().map(|&p| cell_key(p)).collect();
        assert!(
            keys_base == keys_edit,
            "cell keys diverge — vertices crossed cells; fixture aborts"
        );
        keys_base.sort_unstable();
        keys_base.dedup();
        assert!(
            keys_base.len() == base_mesh.positions.len(),
            "cell keys not unique — Surface Nets emits one vertex per cell"
        );

        // ---- changed slots and changed operator rows ----------------------
        #[allow(
            clippy::float_cmp,
            reason = "a slot changed iff its position differs at all; bitwise is the honest reading"
        )]
        let changed: Vec<usize> = (0..base_mesh.positions.len())
            .filter(|&i| base_mesh.positions[i] != edit_mesh.positions[i])
            .collect();
        assert!(
            !changed.is_empty(),
            "reachability: the brush moved no vertex at all"
        );
        let mut in_rows = vec![false; base_mesh.positions.len()];
        for &i in &changed {
            in_rows[i] = true;
        }
        for tri in edit_mesh.indices.chunks_exact(3) {
            let touches = tri
                .iter()
                .any(|&v| changed.binary_search(&(v as usize)).is_ok());
            if touches {
                for &v in tri {
                    in_rows[v as usize] = true;
                }
            }
        }
        let changed_rows: Vec<usize> = (0..in_rows.len()).filter(|&i| in_rows[i]).collect();

        // ---- operators, ordering, base factor ------------------------------
        let (a_base, _) = heat::heat_operator(&base_mesh);
        let (a_edit, _) = heat::heat_operator(&edit_mesh);
        let perm = heat::nested_dissection(&a_base, 8);
        let mut inv = vec![0usize; perm.len()];
        for (new, &old) in perm.iter().enumerate() {
            inv[old] = new;
        }
        let ap_base = heat::permute(&a_base, &perm);
        let ap_edit = heat::permute(&a_edit, &perm);
        let f_base = heat::ldl_factor(&ap_base);
        let seeds: Vec<usize> = changed_rows.iter().map(|&i| inv[i]).collect();
        let closure = heat::ancestor_closure(&f_base.parent, &seeds);
        let closure_rows = closure.iter().filter(|&&m| m).count();

        println!(
            "changed slots {} (H says ≤ 400), changed operator rows {} ({}x), closure {} of {} rows ({:.1}%)",
            changed.len(),
            changed_rows.len(),
            changed_rows.len() as f64 / changed.len().max(1) as f64,
            closure_rows,
            ap_base.n,
            100.0 * closure_rows as f64 / ap_base.n as f64
        );

        // ---- validity oracles, before any timing ---------------------------
        let bound = 100.0 * ap_edit.n as f64 * 2.0f64.powi(-53);
        let mut f_up = clone_factor(&f_base);
        let mut mask = closure.clone();
        heat::refactor_rows(&ap_edit, &mut f_up, &mut mask);
        let res_up = heat::residual_fro(&ap_edit, &f_up);
        assert!(
            res_up <= bound,
            "updated factor residual {res_up:.3e} above the bound {bound:.3e}"
        );
        let f_fresh = heat::ldl_factor(&ap_edit);
        let rhs: Vec<f64> = (0..ap_edit.n).map(|i| ((i % 11) as f64) - 5.0).collect();
        let x_up = heat::solve(&f_up, &rhs);
        let x_new = heat::solve(&f_fresh, &rhs);
        let (mut num, mut den) = (0.0f64, 0.0f64);
        for (a, b) in x_up.iter().zip(&x_new) {
            num += (a - b) * (a - b);
            den += b * b;
        }
        let agree = (num / den).sqrt();
        assert!(
            agree <= 1e-8,
            "update-path and refactor-path solves disagree at {agree:.3e}"
        );
        // Inversion: skip one changed column from the closure — red required.
        let skip = seeds[0];
        let mut bad_mask = closure.clone();
        bad_mask[skip] = false;
        let mut f_bad = clone_factor(&f_base);
        heat::refactor_rows(&ap_edit, &mut f_bad, &mut bad_mask);
        let res_bad = heat::residual_fro(&ap_edit, &f_bad);
        assert!(
            res_bad > bound,
            "inversion failed: skipping changed column {skip} left residual \
             {res_bad:.3e} under the bound — the oracle cannot see a missed row"
        );
        println!(
            "validity: update residual {res_up:.2e} (bound {bound:.2e}), solve agreement {agree:.2e}, \
             skipped-column inversion red at {res_bad:.2e}"
        );

        // ---- the race: interleaved both orders, flops beside wall ----------
        let mut rows_out: Vec<(usize, &str, f64, f64, u64, u64)> = Vec::new();
        for rep in 0..REPS {
            let update = || {
                let mut f = clone_factor(&f_base);
                let mut m = closure.clone();
                let t = Instant::now();
                heat::refactor_rows(&ap_edit, &mut f, &mut m);
                (t.elapsed().as_secs_f64() * 1e3, f.flops)
            };
            let refactor = || {
                let t = Instant::now();
                let f = heat::ldl_factor(&ap_edit);
                (t.elapsed().as_secs_f64() * 1e3, f.flops)
            };
            let (order, (u_ms, u_fl), (r_ms, r_fl)) = if rep % 2 == 0 {
                let u = update();
                let r = refactor();
                ("UR", u, r)
            } else {
                let r = refactor();
                let u = update();
                ("RU", u, r)
            };
            rows_out.push((rep, order, u_ms, r_ms, u_fl, r_fl));
        }

        println!(
            "\n{:>4} {:>6} {:>10} {:>12} {:>13} {:>15}",
            "rep", "order", "update_ms", "refactor_ms", "update_flops", "refactor_flops"
        );
        for &(rep, order, u_ms, r_ms, u_fl, r_fl) in &rows_out {
            println!("{rep:>4} {order:>6} {u_ms:>10.2} {r_ms:>12.2} {u_fl:>13} {r_fl:>15}");
            run.record(&[
                ("rep", rep.to_string()),
                ("order", order.to_string()),
                ("update_ms", format!("{u_ms:.3}")),
                ("refactor_ms", format!("{r_ms:.3}")),
                ("update_flops", u_fl.to_string()),
                ("refactor_flops", r_fl.to_string()),
                ("changed_slots", changed.len().to_string()),
                ("changed_rows", changed_rows.len().to_string()),
                ("closure_rows", closure_rows.to_string()),
            ]);
        }

        let median = |mut v: Vec<f64>| -> f64 {
            v.sort_unstable_by(f64::total_cmp);
            v[v.len() / 2]
        };
        let u_med = median(rows_out.iter().map(|r| r.2).collect());
        let r_med = median(rows_out.iter().map(|r| r.3).collect());
        let wall_ratio = r_med / u_med;
        let flop_ratio = rows_out[0].5 as f64 / rows_out[0].4.max(1) as f64;

        println!();
        let slots_ok = changed.len() <= 400;
        println!(
            "slots: {} -- {} (H says ≤ 400; M-318 extrapolated 346)",
            changed.len(),
            if slots_ok {
                "HELD"
            } else {
                "EXCEEDED — re-scopes M-318's extrapolation"
            }
        );
        let held = wall_ratio >= 20.0 && flop_ratio >= 10.0;
        let dead = wall_ratio < 10.0 || flop_ratio < 10.0;
        println!(
            "race: update {u_med:.2} ms / {} flops vs refactor {r_med:.2} ms / {} flops -- wall {wall_ratio:.1}x, flops {flop_ratio:.1}x -- {}",
            rows_out[0].4,
            rows_out[0].5,
            if held {
                "HELD (≥ 20x wall with the ≥ 10x flop floor)"
            } else if dead {
                "FALSIFIED — the prefactored family is dead for live carving; the intrinsic lane routes to the Closest Point Method"
            } else {
                "H FALSIFIED at the 20x bar, routing UNDECIDED between 10x and 20x — recorded, loudly"
            }
        );
    });
}
