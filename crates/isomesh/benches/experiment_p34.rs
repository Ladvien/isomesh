//! **P-34 — warm-start economics of the feasibility program, counted.**
//!
//! Ticket: R-034b. Pre-registered in the commit before this one; the solver
//! is M-330's, unchanged.
//!
//! ```bash
//! cargo bench --bench experiment_p34
//! ```
//!
//! Writes `docs/experiments/p-34.csv`.
//!
//! # The design in one line
//!
//! A running-bond wall — redundant, so a severing edit can leave a standing
//! structure, which M-330's arch could not — re-solved after twenty edits:
//! ten that leave the load paths alone (weight nudges, small tilts) and ten
//! that sever one (an interior block removed, forces rerouted around the
//! hole). Cost is the alternating-projection iteration count to the 1e-5
//! feasibility decision, cold from zero against warm from the pre-edit
//! solution mapped across by interface identity. The registered claim is the
//! bimodality: non-severing median ratio ≤ 0.15, severing ≥ 0.5, medians
//! more than 3× apart.
//!
//! # Counted, not timed (✗24)
//!
//! Iteration counts are integers and identical on every machine; wall-clock
//! is printed beside them and gates nothing. The decision line is probed
//! every 10 iterations, so counts carry that granularity — recorded, and far
//! finer than the 3× separation under test. No timing A/B is interleaved
//! because no wall-clock comparison is registered (M-197 does not apply).

mod common;

/// Friction coefficient, M-330's.
const MU: f64 = 0.7;
/// Iteration cap and decision threshold, M-330's.
const CAP: usize = 20_000;
const FEASIBLE_BELOW: f64 = 1e-5;
/// Decision-probe granularity, iterations.
const PROBE: usize = 10;
/// Wall: 8 courses of 0.5 × 0.25 blocks, 12 per even course.
const COURSES: usize = 8;
const PER_COURSE: usize = 12;
const WB: f64 = 0.5;
const HB: f64 = 0.25;

#[derive(Clone)]
struct Block {
    x0: f64,
    x1: f64,
    course: usize,
}

/// Interface identity that survives block removal: endpoint blocks by their
/// ORIGINAL ids (`None` = ground), plus which joint kind.
#[derive(Clone, PartialEq)]
struct Iface {
    giver: Option<usize>,
    receiver: usize,
    kind: u8, // 0 = bed/ground (normal +y), 1 = head (normal +x)
    v: [[f64; 2]; 2],
    n: [f64; 2],
    t: [f64; 2],
}

fn build_blocks() -> Vec<Block> {
    let mut blocks = Vec::new();
    for c in 0..COURSES {
        if c % 2 == 0 {
            for i in 0..PER_COURSE {
                blocks.push(Block {
                    x0: i as f64 * WB,
                    x1: (i + 1) as f64 * WB,
                    course: c,
                });
            }
        } else {
            for i in 0..PER_COURSE - 1 {
                blocks.push(Block {
                    x0: WB / 2.0 + i as f64 * WB,
                    x1: WB / 2.0 + (i + 1) as f64 * WB,
                    course: c,
                });
            }
        }
    }
    blocks
}

/// Interfaces for a given alive-set, with identities in original block ids.
fn build_ifaces(blocks: &[Block], alive: &[bool]) -> Vec<Iface> {
    let mut ifaces = Vec::new();
    for (b, blk) in blocks.iter().enumerate() {
        if !alive[b] {
            continue;
        }
        let y0 = blk.course as f64 * HB;
        // Ground under course 0.
        if blk.course == 0 {
            ifaces.push(Iface {
                giver: None,
                receiver: b,
                kind: 0,
                v: [[blk.x0, 0.0], [blk.x1, 0.0]],
                n: [0.0, 1.0],
                t: [1.0, 0.0],
            });
        }
        // Bed joints: this block resting on course-below overlaps.
        if blk.course > 0 {
            for (lo, lob) in blocks.iter().enumerate() {
                if !alive[lo] || lob.course + 1 != blk.course {
                    continue;
                }
                let ox0 = blk.x0.max(lob.x0);
                let ox1 = blk.x1.min(lob.x1);
                if ox1 - ox0 > 1e-9 {
                    ifaces.push(Iface {
                        giver: Some(lo),
                        receiver: b,
                        kind: 0,
                        v: [[ox0, y0], [ox1, y0]],
                        n: [0.0, 1.0],
                        t: [1.0, 0.0],
                    });
                }
            }
        }
        // Head joint with the right-hand neighbour in the same course.
        for (r, rb) in blocks.iter().enumerate() {
            if !alive[r] || rb.course != blk.course {
                continue;
            }
            if (rb.x0 - blk.x1).abs() < 1e-9 {
                ifaces.push(Iface {
                    giver: Some(b),
                    receiver: r,
                    kind: 1,
                    v: [[blk.x1, y0], [blk.x1, y0 + HB]],
                    n: [1.0, 0.0],
                    t: [0.0, 1.0],
                });
            }
        }
    }
    ifaces
}

struct Problem {
    a: Vec<f64>,
    rows: usize,
    cols: usize,
    b: Vec<f64>,
    total_w: f64,
}

/// Assemble the program for the alive-set under gravity tilted by `tilt`,
/// with per-block weight scale factors.
fn assemble(
    blocks: &[Block],
    alive: &[bool],
    ifaces: &[Iface],
    tilt: f64,
    wscale: &[f64],
) -> Problem {
    let mut row_of = vec![None; blocks.len()];
    let mut nb = 0usize;
    for (i, &a) in alive.iter().enumerate() {
        if a {
            row_of[i] = Some(nb);
            nb += 1;
        }
    }
    let rows = 3 * nb;
    let cols = 4 * ifaces.len();
    let mut a = vec![0.0f64; rows * cols];
    let mut b = vec![0.0f64; rows];
    let gv = [tilt.sin(), -tilt.cos()];
    let mut total_w = 0.0;
    for (i, blk) in blocks.iter().enumerate() {
        let Some(r) = row_of[i] else { continue };
        let w = (blk.x1 - blk.x0) * HB * wscale[i];
        total_w += w;
        b[3 * r] = -w * gv[0];
        b[3 * r + 1] = -w * gv[1];
    }
    for (k, f) in ifaces.iter().enumerate() {
        for (v, p) in f.v.iter().enumerate() {
            let col_n = 4 * k + 2 * v;
            let col_t = col_n + 1;
            for (blk, sign) in [(Some(f.receiver), 1.0f64), (f.giver, -1.0f64)] {
                let Some(bi) = blk else { continue };
                let Some(r) = row_of[bi] else { continue };
                let blkref = &blocks[bi];
                let c = [
                    (blkref.x0 + blkref.x1) / 2.0,
                    (blkref.course as f64 + 0.5) * HB,
                ];
                let (fx, fy, tq) = (3 * r, 3 * r + 1, 3 * r + 2);
                a[fx * cols + col_n] += sign * f.n[0];
                a[fx * cols + col_t] += sign * f.t[0];
                a[fy * cols + col_n] += sign * f.n[1];
                a[fy * cols + col_t] += sign * f.t[1];
                let arm = [p[0] - c[0], p[1] - c[1]];
                a[tq * cols + col_n] += sign * (arm[0] * f.n[1] - arm[1] * f.n[0]);
                a[tq * cols + col_t] += sign * (arm[0] * f.t[1] - arm[1] * f.t[0]);
            }
        }
    }
    Problem {
        a,
        rows,
        cols,
        b,
        total_w,
    }
}

#[allow(
    clippy::needless_range_loop,
    reason = "dense matrix kernels index two arrays per loop; iterator forms obscure the algebra"
)]
struct AffineProjector {
    rows: usize,
    cols: usize,
    a: Vec<f64>,
    l: Vec<f64>,
}

#[allow(
    clippy::needless_range_loop,
    reason = "dense matrix kernels index two arrays per loop; iterator forms obscure the algebra"
)]
impl AffineProjector {
    fn new(p: &Problem) -> Self {
        let (rows, cols) = (p.rows, p.cols);
        let mut m = vec![0.0f64; rows * rows];
        for i in 0..rows {
            for j in 0..=i {
                let mut s = 0.0;
                for k in 0..cols {
                    s += p.a[i * cols + k] * p.a[j * cols + k];
                }
                m[i * rows + j] = s;
                m[j * rows + i] = s;
            }
        }
        for j in 0..rows {
            for k in 0..j {
                let ljk = m[j * rows + k];
                if ljk == 0.0 {
                    continue;
                }
                for i in j..rows {
                    m[i * rows + j] -= m[i * rows + k] * ljk;
                }
            }
            let d = m[j * rows + j];
            assert!(d > 1e-12, "A·Aᵀ lost rank at row {j} — formulation failure");
            let root = d.sqrt();
            for i in j..rows {
                m[i * rows + j] /= root;
            }
        }
        AffineProjector {
            rows,
            cols,
            a: p.a.clone(),
            l: m,
        }
    }

    fn project(&self, f: &mut [f64], b: &[f64], scratch: &mut [f64]) {
        let (rows, cols) = (self.rows, self.cols);
        for i in 0..rows {
            let mut s = -b[i];
            for k in 0..cols {
                s += self.a[i * cols + k] * f[k];
            }
            scratch[i] = s;
        }
        for i in 0..rows {
            let mut s = scratch[i];
            for k in 0..i {
                s -= self.l[i * rows + k] * scratch[k];
            }
            scratch[i] = s / self.l[i * rows + i];
        }
        for i in (0..rows).rev() {
            let mut s = scratch[i];
            for k in i + 1..rows {
                s -= self.l[k * rows + i] * scratch[k];
            }
            scratch[i] = s / self.l[i * rows + i];
        }
        for k in 0..cols {
            let mut s = 0.0;
            for i in 0..rows {
                s += self.a[i * cols + k] * scratch[i];
            }
            f[k] -= s;
        }
    }
}

fn project_cone(f: &mut [f64]) {
    for pair in f.chunks_exact_mut(2) {
        let n = pair[0];
        let t = pair[1];
        if t.abs() <= MU * n {
            continue;
        }
        if MU * t.abs() <= -n {
            pair[0] = 0.0;
            pair[1] = 0.0;
        } else {
            let nn = (n + MU * t.abs()) / (1.0 + MU * MU);
            pair[0] = nn;
            pair[1] = t.signum() * MU * nn;
        }
    }
}

#[allow(
    clippy::needless_range_loop,
    reason = "dense matrix kernels index two arrays per loop; iterator forms obscure the algebra"
)]
fn residual_rel(p: &Problem, f: &[f64]) -> f64 {
    let mut e = 0.0f64;
    for i in 0..p.rows {
        let mut s = -p.b[i];
        for k in 0..p.cols {
            s += p.a[i * p.cols + k] * f[k];
        }
        e += s * s;
    }
    e.sqrt() / p.total_w
}

/// Iterations (probed every `PROBE`) until the cone-side residual crosses the
/// feasibility line. Returns (iterations, reached, final f).
fn solve_counted(p: &Problem, start: &[f64]) -> (usize, bool, Vec<f64>) {
    let proj = AffineProjector::new(p);
    let mut f = start.to_vec();
    let mut scratch = vec![0.0f64; p.rows];
    let mut done = None;
    let mut it = 0usize;
    while it < CAP {
        for _ in 0..PROBE {
            proj.project(&mut f, &p.b, &mut scratch);
            project_cone(&mut f);
        }
        it += PROBE;
        if residual_rel(p, &f) < FEASIBLE_BELOW {
            done = Some(it);
            break;
        }
    }
    (done.unwrap_or(CAP), done.is_some(), f)
}

struct Lcg(u64);

impl Lcg {
    fn next_u64(&mut self) -> u64 {
        self.0 = self
            .0
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        self.0
    }

    fn pick(&mut self, n: usize) -> usize {
        (self.next_u64() >> 33) as usize % n
    }
}

fn median(mut v: Vec<f64>) -> f64 {
    assert!(!v.is_empty(), "median of an empty class — reachability");
    v.sort_unstable_by(f64::total_cmp);
    v[v.len() / 2]
}

fn main() {
    if !std::env::args().any(|arg| arg == "--bench") {
        return;
    }

    let prereg = isomesh::experiment!("P-34");
    common::experiment::run(prereg, |run| {
        let blocks = build_blocks();
        let alive = vec![true; blocks.len()];
        let ifaces = build_ifaces(&blocks, &alive);
        let wscale = vec![1.0f64; blocks.len()];
        let base_p = assemble(&blocks, &alive, &ifaces, 0.0, &wscale);
        let zero = vec![0.0f64; base_p.cols];
        let (base_iters, base_ok, base_f) = solve_counted(&base_p, &zero);
        assert!(base_ok, "the base wall did not stand — fixture failure");
        println!(
            "base wall: {} blocks, {} interfaces, feasible in {} iterations",
            blocks.len(),
            ifaces.len(),
            base_iters
        );

        let mut rows: Vec<(String, String, usize, usize, f64, bool)> = Vec::new();
        let mut nonsever = Vec::new();
        let mut sever = Vec::new();
        let mut collapsed = 0usize;

        // --- non-severing: 5 weight nudges + 5 small tilts -----------------
        let mut lcg = Lcg(0xA5A5_1234_5678_9ABC);
        for e in 0..5 {
            let mut ws = wscale.clone();
            ws[lcg.pick(blocks.len())] = 1.10;
            let p = assemble(&blocks, &alive, &ifaces, 0.0, &ws);
            let (ic, okc, _) = solve_counted(&p, &zero);
            let (iw, okw, _) = solve_counted(&p, &base_f);
            assert!(okc && okw, "non-severing edit failed to decide in the cap");
            let ratio = iw as f64 / ic as f64;
            nonsever.push(ratio);
            rows.push((
                format!("weight_{e}"),
                "non_severing".into(),
                ic,
                iw,
                ratio,
                true,
            ));
        }
        for (e, deg) in [0.3f64, -0.3, 0.5, -0.5, 0.7].iter().enumerate() {
            let p = assemble(&blocks, &alive, &ifaces, deg.to_radians(), &wscale);
            let (ic, okc, _) = solve_counted(&p, &zero);
            let (iw, okw, _) = solve_counted(&p, &base_f);
            assert!(okc && okw, "tilt edit failed to decide in the cap");
            let ratio = iw as f64 / ic as f64;
            nonsever.push(ratio);
            rows.push((
                format!("tilt_{e}"),
                "non_severing".into(),
                ic,
                iw,
                ratio,
                true,
            ));
        }

        // --- severing: 10 interior removals, forces rerouted ---------------
        let eligible: Vec<usize> = blocks
            .iter()
            .enumerate()
            .filter(|(_, b)| (2..=5).contains(&b.course) && b.x0 > 0.6 && b.x1 < 5.4)
            .map(|(i, _)| i)
            .collect();
        let mut picked = Vec::new();
        while picked.len() < 10 {
            let c = eligible[lcg.pick(eligible.len())];
            if !picked.contains(&c) {
                picked.push(c);
            }
        }
        for (e, &gone) in picked.iter().enumerate() {
            let mut alive2 = alive.clone();
            alive2[gone] = false;
            let ifaces2 = build_ifaces(&blocks, &alive2);
            let p = assemble(&blocks, &alive2, &ifaces2, 0.0, &wscale);
            // Warm start: carry forces across by interface identity.
            let mut warm = vec![0.0f64; p.cols];
            for (k2, f2) in ifaces2.iter().enumerate() {
                if let Some(k1) = ifaces.iter().position(|f1| {
                    f1.giver == f2.giver && f1.receiver == f2.receiver && f1.kind == f2.kind
                }) {
                    warm[4 * k2..4 * k2 + 4].copy_from_slice(&base_f[4 * k1..4 * k1 + 4]);
                }
            }
            let (ic, okc, _) = solve_counted(&p, &vec![0.0; p.cols]);
            let (iw, okw, _) = solve_counted(&p, &warm);
            if !(okc && okw) {
                collapsed += 1;
                rows.push((
                    format!("remove_{e}"),
                    "collapsed".into(),
                    ic,
                    iw,
                    f64::NAN,
                    false,
                ));
                continue;
            }
            let ratio = iw as f64 / ic as f64;
            sever.push(ratio);
            rows.push((
                format!("remove_{e}"),
                "severing".into(),
                ic,
                iw,
                ratio,
                true,
            ));
        }
        assert!(
            collapsed <= 3,
            "{collapsed} of 10 removals collapsed — the corpus is a fixture failure"
        );
        assert!(
            !sever.is_empty(),
            "no severing edit survived — reachability"
        );

        // --- emit -----------------------------------------------------------
        println!(
            "\n{:>12} {:>14} {:>11} {:>11} {:>8} {:>9}",
            "edit", "class", "iters_cold", "iters_warm", "ratio", "feasible"
        );
        for (edit, class, ic, iw, ratio, ok) in &rows {
            println!(
                "{:>12} {:>14} {:>11} {:>11} {:>8.3} {:>9}",
                edit, class, ic, iw, ratio, ok
            );
            run.record(&[
                ("edit", edit.clone()),
                ("class", class.clone()),
                ("iters_cold", ic.to_string()),
                ("iters_warm", iw.to_string()),
                ("ratio", format!("{ratio:.4}")),
                ("feasible", ok.to_string()),
            ]);
        }

        let m_non = median(nonsever);
        let m_sev = median(sever);
        println!();
        let c1 = m_non <= 0.15;
        let c2 = m_sev >= 0.5;
        let separated = m_sev / m_non.max(1e-9) > 3.0;
        println!(
            "C1 (non-severing cheap): median ratio {m_non:.3} -- {} (H says <= 0.15)",
            if c1 { "HELD" } else { "FALSIFIED" }
        );
        println!(
            "C2 (severing expensive): median ratio {m_sev:.3} -- {} (H says >= 0.5)",
            if c2 { "HELD" } else { "FALSIFIED" }
        );
        println!(
            "bimodality: medians separated {:.1}x -- {} (H says > 3x; unimodal kills the \
             cheap-incremental story)",
            m_sev / m_non.max(1e-9),
            if separated { "HELD" } else { "FALSIFIED" }
        );
        if collapsed > 0 {
            println!("collapsed removals excluded, on the record: {collapsed}");
        }
    });
}
