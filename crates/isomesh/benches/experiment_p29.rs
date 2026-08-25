//! **P-29 — does the wormhole competition reproduce, with the paper's own constants?**
//!
//! Ticket: R-031. Pre-registered in the commit before this one.
//!
//! ```bash
//! cargo bench --bench experiment_p29
//! ```
//!
//! Writes `docs/experiments/p-29.csv`.
//!
//! # The model, verbatim from the source
//!
//! Dreybrodt & Gabrovšek, *Dynamics of wormhole formation in fractured
//! limestones* (`10.5194/hess-23-1995-2019`, READ): a 2-D square net of 1-D
//! fractures (length 200 cm, width 100 cm, aperture a₀ = 0.02 cm), constant
//! heads 15 m → 0 with **periodic** transverse boundaries; cubic-law
//! resistance `R = 12η/(ρg)·Σ Δx/(a³b)`; linear kinetics `F = k·(1 − c/c_eq)`
//! with the composite `k = k₁/(1 + k₁a/(6Dc_eq))`, `k₁ = 4e-11 mol cm⁻² s⁻¹`,
//! `c_eq = 1e-6 mol cm⁻³`, `D = 1e-5 cm² s⁻¹`; closed-form transport per
//! segment (the exponential is the exact solution of the linear law, so there
//! is no Δc stability question); widening `da/dt = 2γF` with γ = calcite's
//! molar volume, 36.93 cm³ mol⁻¹. Seeds are the paper's: `Δa = 1e-9·a₀` on
//! the first ten downstream fractures of a chosen input. The heterogeneous
//! net is the paper's lognormal (0.015–0.025 cm about 0.02, σ = 0.2). The
//! recharge-limited arm is Perne, Covington & Gabrovšek
//! (`10.5194/hess-18-4617-2014`): cap the recharge and *"further expansion of
//! the network is suppressed"*.
//!
//! # Determinism, and why the flow solve is direct
//!
//! Cubic-law contrast reaches (a_max/a_min)³ ~ 1e6 exactly when the dynamics
//! get interesting, which is where a fixed-iteration Krylov solve would turn
//! the pressure field into noise with an unpredictable direction. So: one
//! banded Cholesky per tick (y-major indexing, periodic wrap Δ = 63,
//! bandwidth 65, Dirichlet columns eliminated), exact to float. The time step
//! is a pure function of state — a fixed fraction of the fastest relative
//! widening (10% primary, 5% control arm; the registration's pre-run
//! correction, since the always-fresh dissolution front binds the step and a
//! 1% cap would mean ~1e6 ticks), clamped to [1e-3, 1] yr, with the per-tick
//! assert scaled to the cap. The control arm must reproduce every registered
//! verdict or the run fails loudly. The only randomness is the registered LCG
//! behind the lognormal apertures.
//!
//! # The detector's own red and green, run before any verdict
//!
//! The central-gap statistic (max consecutive gap in sorted ln a, 1% tails
//! dropped, threshold 0.2 = the initial log-sd) must call the t = 0
//! lognormal apertures **unimodal** and a synthetic half-shifted (+4σ) sample
//! **bimodal**, or the run aborts. The recharge arm is the physics inversion:
//! both-arms-bimodal indicts the statistic, both-unimodal the kinetics.
//!
//! # Timed only incidentally
//!
//! The registered clauses are distribution shapes and a concentration share —
//! counts and ratios, machine-independent. Tick cost is recorded because the
//! dossier asked, and gates nothing (✗24). No timing A/B exists, so M-197's
//! interleaving rule does not apply.

mod common;

use std::time::Instant;

/// Nodes per axis; x is the head gradient direction, y is periodic.
const N: usize = 64;
/// Segments per fracture. Δx = 25 cm ≈ 0.8·λ₀ resolves the initial front.
const M: usize = 8;
/// Fracture length and width, cm (the paper's).
const FRAC_LEN: f64 = 200.0;
const FRAC_WIDTH: f64 = 100.0;
const DX: f64 = FRAC_LEN / M as f64;
/// Initial aperture, cm.
const A0: f64 = 0.02;
/// Input head, cm (15 m); output is 0.
const H_IN: f64 = 1500.0;
/// Kinetics (the paper's): k₁ [mol cm⁻² s⁻¹], c_eq [mol cm⁻³], D [cm² s⁻¹].
const K1: f64 = 4e-11;
const CEQ: f64 = 1e-6;
const DIFF: f64 = 1e-5;
/// Water: viscosity [poise], ρ·g [g cm⁻² s⁻²].
const VISC: f64 = 0.01;
const RHO_G: f64 = 981.0;
/// Calcite molar volume, cm³ mol⁻¹ (100.09 g/mol ÷ 2.71 g/cm³) — γ's factor.
const MOLAR_VOL: f64 = 36.93;
const SECS_PER_YEAR: f64 = 3.155_76e7;
/// The registered bimodality threshold: the initial log-sd.
const GAP_THRESHOLD: f64 = 0.2;
/// Lognormal parameters (the paper's heterogeneous net).
const LOGSD: f64 = 0.2;
const A_MIN: f64 = 0.015;
const A_MAX: f64 = 0.025;
/// Seed positions (input y rows) for the homogeneous arm — gaps 17/20/27,
/// inside the paper's interaction range (≤ 30 nodes).
const SEEDS: [usize; 3] = [10, 27, 47];
/// Hard caps: a constant-head arm that has not broken through by here is a
/// fixture failure, loudly.
const MAX_YEARS_HEAD: f64 = 20_000.0;
const MAX_YEARS_RECHARGE: f64 = 20_000.0;
/// Per-tick relative-widening caps: the primary run and the halved control
/// arm whose verdicts must agree (the pre-run correction in the P-29
/// registration; the aperture ODE is self-limiting, so the cap is an accuracy
/// knob, and the control arm is the measurement of that claim).
const DT_CAP: f64 = 0.10;
const DT_CAP_CONTROL: f64 = 0.05;

/// The composite rate constant k(a): surface reaction + diffusion in series.
fn k_of(a: f64) -> f64 {
    K1 / (1.0 + K1 * a / (6.0 * DIFF * CEQ))
}

/// Deterministic LCG, the p26 pattern: top bits only.
struct Lcg(u64);

impl Lcg {
    fn next_u64(&mut self) -> u64 {
        self.0 = self
            .0
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        self.0
    }

    /// Uniform in (0, 1] — the offset makes ln() safe.
    fn unit(&mut self) -> f64 {
        (((self.next_u64() >> 11) as f64) + 1.0) / (1u64 << 53) as f64
    }

    /// Standard normal via Box–Muller (libm through std's f64 methods).
    fn normal(&mut self) -> f64 {
        let u1 = self.unit();
        let u2 = self.unit();
        (-2.0 * u1.ln()).sqrt() * (2.0 * std::f64::consts::PI * u2).cos()
    }
}

/// One fracture: `M` aperture segments between two nodes.
struct Edge {
    from: usize,
    to: usize,
    a: [f64; M],
}

impl Edge {
    fn resistance(&self) -> f64 {
        let mut r = 0.0;
        for a in &self.a {
            r += 12.0 * VISC * DX / (RHO_G * a * a * a * FRAC_WIDTH);
        }
        r
    }

    fn mean_a(&self) -> f64 {
        self.a.iter().sum::<f64>() / M as f64
    }
}

fn node(x: usize, y: usize) -> usize {
    x * N + y
}

/// Build the lattice: horizontal fractures x→x+1 and periodic vertical ones.
fn build_edges(apertures: impl FnMut() -> f64) -> Vec<Edge> {
    let mut ap = apertures;
    let mut edges = Vec::new();
    for x in 0..N - 1 {
        for y in 0..N {
            let a = ap();
            edges.push(Edge {
                from: node(x, y),
                to: node(x + 1, y),
                a: [a; M],
            });
        }
    }
    for x in 0..N {
        for y in 0..N {
            let a = ap();
            edges.push(Edge {
                from: node(x, y),
                to: node(x, (y + 1) % N),
                a: [a; M],
            });
        }
    }
    edges
}

/// Boundary condition of an arm.
enum Boundary {
    /// Dirichlet H_IN on the x = 0 column, 0 on the x = N−1 column.
    ConstantHead,
    /// Fixed influx per x = 0 node; Dirichlet 0 on the x = N−1 column only.
    FixedFlux(f64),
}

/// Banded SPD solve, column-major: `cols[j·(bw+1) + (i−j)]` holds `L[i][j]`
/// for `i ∈ [j, j+bw]` — the left-looking update then runs contiguously over
/// each column, which is what makes an exact solve per tick affordable.
/// Panics loudly if the matrix loses positive definiteness — a model failure,
/// not a condition to paper over.
struct Banded {
    n: usize,
    bw: usize,
    cols: Vec<f64>,
}

impl Banded {
    fn new(n: usize, bw: usize) -> Self {
        Banded {
            n,
            bw,
            cols: vec![0.0; n * (bw + 1)],
        }
    }

    fn at(&mut self, i: usize, j: usize) -> &mut f64 {
        debug_assert!(j <= i && i - j <= self.bw);
        &mut self.cols[j * (self.bw + 1) + (i - j)]
    }

    fn factor(&mut self) {
        let bw = self.bw;
        let w = bw + 1;
        for j in 0..self.n {
            let start = j.saturating_sub(bw);
            for k in start..j {
                let off = j - k;
                let ljk = self.cols[k * w + off];
                if ljk == 0.0 {
                    continue;
                }
                // L[i][j] -= L[i][k]·L[j][k] for i = j .. min(k+bw, j+bw, n−1);
                // both runs are contiguous slices of their columns.
                let len = (w - off).min(self.n - j);
                let (ka, ja) = (k * w + off, j * w);
                for d in 0..len {
                    self.cols[ja + d] -= ljk * self.cols[ka + d];
                }
            }
            let d = self.cols[j * w];
            assert!(d > 0.0, "flow matrix lost positive definiteness at {j}");
            let root = d.sqrt();
            let len = w.min(self.n - j);
            for v in &mut self.cols[j * w..j * w + len] {
                *v /= root;
            }
        }
    }

    fn solve(&self, rhs: &mut [f64]) {
        let bw = self.bw;
        let w = bw + 1;
        // Forward: L y = b, walking each finished y_j down its column.
        for j in 0..self.n {
            let yj = rhs[j] / self.cols[j * w];
            rhs[j] = yj;
            let len = w.min(self.n - j);
            for d in 1..len {
                rhs[j + d] -= self.cols[j * w + d] * yj;
            }
        }
        // Backward: Lᵀ x = y; row i of Lᵀ is column i of L, contiguous.
        for i in (0..self.n).rev() {
            let mut s = rhs[i];
            let len = w.min(self.n - i);
            for d in 1..len {
                s -= self.cols[i * w + d] * rhs[i + d];
            }
            rhs[i] = s / self.cols[i * w];
        }
    }
}

/// Solve the flow, returning per-node heads.
fn solve_heads(edges: &[Edge], boundary: &Boundary) -> Vec<f64> {
    // Unknown nodes: interior columns, plus the input column under FixedFlux.
    let first_unknown_col = match boundary {
        Boundary::ConstantHead => 1,
        Boundary::FixedFlux(_) => 0,
    };
    let cols = (N - 1) - first_unknown_col;
    let n_unknown = cols * N;
    let bw = N;
    let unknown = |nd: usize| -> Option<usize> {
        let x = nd / N;
        if x >= first_unknown_col && x < N - 1 {
            Some((x - first_unknown_col) * N + (nd % N))
        } else {
            None
        }
    };
    let mut mat = Banded::new(n_unknown, bw);
    let mut rhs = vec![0.0; n_unknown];
    if let Boundary::FixedFlux(q) = boundary {
        for y in 0..N {
            rhs[node(0, y)] += q;
        }
    }
    for e in edges {
        let c = 1.0 / e.resistance();
        match (unknown(e.from), unknown(e.to)) {
            (Some(i), Some(j)) => {
                *mat.at(i, i) += c;
                *mat.at(j, j) += c;
                let (hi, lo) = if i > j { (i, j) } else { (j, i) };
                *mat.at(hi, lo) -= c;
            }
            (Some(i), None) | (None, Some(i)) => {
                *mat.at(i, i) += c;
                let fixed = if unknown(e.from).is_none() {
                    e.from
                } else {
                    e.to
                };
                let h_fixed = if fixed / N == 0 { H_IN } else { 0.0 };
                rhs[i] += c * h_fixed;
            }
            (None, None) => {}
        }
    }
    mat.factor();
    mat.solve(&mut rhs);
    let mut heads = vec![0.0; N * N];
    for (nd, h) in heads.iter_mut().enumerate() {
        let x = nd / N;
        if x == N - 1 {
            *h = 0.0;
        } else if x < first_unknown_col {
            *h = H_IN;
        } else {
            *h = rhs[(x - first_unknown_col) * N + (nd % N)];
        }
    }
    heads
}

/// One transport-and-widening tick. Returns (dt_years, max da/a, dissolved
/// mol, per-edge dissolution rate mol/s, per-edge |Q|).
fn tick(
    edges: &mut [Edge],
    heads: &[f64],
    boundary: &Boundary,
    cap: f64,
) -> (f64, f64, f64, Vec<f64>, Vec<f64>) {
    let flows: Vec<f64> = edges
        .iter()
        .map(|e| (heads[e.from] - heads[e.to]) / e.resistance())
        .collect();

    // Node processing order: descending head (ties by index) guarantees every
    // inflow's exit concentration is known before its node is reached.
    let mut order: Vec<usize> = (0..N * N).collect();
    order.sort_unstable_by(|&p, &q| heads[q].total_cmp(&heads[p]).then(p.cmp(&q)));

    // Inflow bookkeeping: (Σ Q·c, Σ Q) per node. Injected recharge is an
    // inflow at c = 0; Dirichlet inputs likewise start fresh water.
    let mut qc = vec![0.0; N * N];
    let mut qsum = vec![0.0; N * N];
    if let Boundary::FixedFlux(q) = boundary {
        for y in 0..N {
            qsum[node(0, y)] += q;
        }
    }
    let mut edges_of: Vec<Vec<usize>> = vec![Vec::new(); N * N];
    for (i, e) in edges.iter().enumerate() {
        edges_of[e.from].push(i);
        edges_of[e.to].push(i);
    }

    let mut rate = vec![[0.0f64; M]; edges.len()];
    let mut dissolution = vec![0.0f64; edges.len()];
    let mut max_rel = 0.0f64;
    for &nd in &order {
        let fresh_input = matches!(boundary, Boundary::ConstantHead) && nd / N == 0;
        let c_node = if fresh_input {
            0.0
        } else if qsum[nd] > 1e-12 {
            (qc[nd] / qsum[nd]).min(CEQ)
        } else {
            CEQ // stagnant water sits at equilibrium and dissolves nothing
        };
        for &ei in &edges_of[nd] {
            let e = &edges[ei];
            let q = flows[ei];
            // Only edges LEAVING this node are transported now.
            let (downstream, qabs) = if q > 0.0 && e.from == nd {
                (e.to, q)
            } else if q < 0.0 && e.to == nd {
                (e.from, -q)
            } else {
                continue;
            };
            if qabs < 1e-12 {
                continue;
            }
            let mut c = c_node;
            let mut diss = 0.0;
            for (s, a) in edges[ei].a.iter().enumerate() {
                let k = k_of(*a);
                let p = 2.0 * (a + FRAC_WIDTH);
                let alpha = k * p * DX / (qabs * CEQ);
                let f_in = k * (1.0 - c / CEQ).max(0.0);
                let mean_f = if alpha > 1e-12 {
                    f_in * (1.0 - (-alpha).exp()) / alpha
                } else {
                    f_in
                };
                rate[ei][s] = 2.0 * mean_f * MOLAR_VOL * SECS_PER_YEAR; // cm/yr
                diss += mean_f * p * DX; // mol/s
                c = CEQ - (CEQ - c) * (-alpha).exp();
            }
            dissolution[ei] = diss;
            qc[downstream] += qabs * c;
            qsum[downstream] += qabs;
        }
    }

    // dt: 1% of the fastest relative widening, a pure function of state.
    for (ei, e) in edges.iter().enumerate() {
        for (s, a) in e.a.iter().enumerate() {
            if rate[ei][s] > 0.0 {
                max_rel = max_rel.max(rate[ei][s] / a);
            }
        }
    }
    let dt = if max_rel > 0.0 {
        (cap / max_rel).clamp(1e-3, 1.0)
    } else {
        1.0
    };

    let mut dissolved = 0.0;
    let mut max_step = 0.0f64;
    for (ei, e) in edges.iter_mut().enumerate() {
        for (s, a) in e.a.iter_mut().enumerate() {
            let da = rate[ei][s] * dt;
            max_step = max_step.max(da / *a);
            *a += da;
        }
        dissolved += dissolution[ei] * dt * SECS_PER_YEAR;
    }
    assert!(
        max_step <= cap * 1.01 + 1e-6,
        "per-tick widening exceeded the {cap} cap: {max_step}"
    );
    let qabs: Vec<f64> = flows.iter().map(|q| q.abs()).collect();
    (dt, max_step, dissolved, dissolution, qabs)
}

/// Is there an input→output path over fractures with mean aperture ≥ 2a₀?
fn broken_through(edges: &[Edge]) -> bool {
    let mut adj: Vec<Vec<usize>> = vec![Vec::new(); N * N];
    for e in edges {
        if e.mean_a() >= 2.0 * A0 {
            adj[e.from].push(e.to);
            adj[e.to].push(e.from);
        }
    }
    let mut seen = vec![false; N * N];
    let mut queue: Vec<usize> = (0..N).map(|y| node(0, y)).collect();
    for &q in &queue {
        seen[q] = true;
    }
    while let Some(u) = queue.pop() {
        if u / N == N - 1 {
            return true;
        }
        for &v in &adj[u] {
            if !seen[v] {
                seen[v] = true;
                queue.push(v);
            }
        }
    }
    false
}

/// Max consecutive gap in sorted ln(mean aperture), 1% tails dropped.
fn central_gap(edges_ln_a: &mut [f64]) -> f64 {
    edges_ln_a.sort_unstable_by(f64::total_cmp);
    let trim = edges_ln_a.len() / 100;
    let core = &edges_ln_a[trim..edges_ln_a.len() - trim];
    core.windows(2).map(|w| w[1] - w[0]).fold(0.0, f64::max)
}

fn gini(values: &mut [f64]) -> f64 {
    values.sort_unstable_by(f64::total_cmp);
    let n = values.len() as f64;
    let total: f64 = values.iter().sum();
    if total <= 0.0 {
        return 0.0;
    }
    let weighted: f64 = values
        .iter()
        .enumerate()
        .map(|(i, v)| (i as f64 + 1.0) * v)
        .sum();
    (2.0 * weighted) / (n * total) - (n + 1.0) / n
}

struct ArmResult {
    ticks: u64,
    years: f64,
    breakthrough_years: Option<f64>,
    dissolved_at_breakthrough: Option<f64>,
    max_gap_ln: f64,
    flux_top10_pct: f64,
    gini_flow: f64,
    max_da_over_a_pct: f64,
    tick_ms_median: f64,
    dissolved: f64,
    t0_inflow_per_node: f64,
}

/// Run one arm to its stopping condition.
fn run_arm(
    edges: &mut [Edge],
    boundary: &Boundary,
    stop_dissolved: Option<f64>,
    max_years: f64,
    cap: f64,
) -> ArmResult {
    let mut years = 0.0;
    let mut ticks = 0u64;
    let mut breakthrough: Option<f64> = None;
    let mut dissolved_at_bt: Option<f64> = None;
    let mut dissolved = 0.0;
    let mut worst_step = 0.0f64;
    let mut tick_ms = Vec::new();
    let mut t0_inflow = 0.0;
    let (last_diss, last_q) = loop {
        let start = Instant::now();
        let heads = solve_heads(edges, boundary);
        if ticks == 0 {
            // t = 0 total inflow across the input column, per node.
            let total: f64 = edges
                .iter()
                .filter(|e| e.from / N == 0 && e.to / N == 1)
                .map(|e| (heads[e.from] - heads[e.to]) / e.resistance())
                .sum();
            t0_inflow = total / N as f64;
        }
        let (dt, step, diss, per_edge_diss, qabs) = tick(edges, &heads, boundary, cap);
        tick_ms.push(start.elapsed().as_secs_f64() * 1e3);
        years += dt;
        ticks += 1;
        dissolved += diss;
        worst_step = worst_step.max(step);
        if breakthrough.is_none() && broken_through(edges) {
            breakthrough = Some(years);
            dissolved_at_bt = Some(dissolved);
        }
        let done = match (breakthrough, stop_dissolved) {
            (_, Some(target)) => dissolved >= target,
            (Some(t), None) => years >= 1.2 * t,
            (None, None) => false,
        };
        if done || years >= max_years {
            break (per_edge_diss, qabs);
        }
    };
    let mut ln_a: Vec<f64> = edges.iter().map(|e| e.mean_a().ln()).collect();
    let gap = central_gap(&mut ln_a);
    let total_diss: f64 = last_diss.iter().sum();
    assert!(
        total_diss > 0.0,
        "reachability: zero dissolution flux at the evaluation tick"
    );
    let mut sorted = last_diss;
    sorted.sort_unstable_by(|x, y| y.total_cmp(x));
    let top = sorted.len() / 10;
    let share: f64 = sorted[..top].iter().sum::<f64>() / total_diss;
    tick_ms.sort_unstable_by(f64::total_cmp);
    let mut last_q = last_q;
    ArmResult {
        ticks,
        years,
        breakthrough_years: breakthrough,
        dissolved_at_breakthrough: dissolved_at_bt,
        max_gap_ln: gap,
        flux_top10_pct: share * 100.0,
        gini_flow: gini(&mut last_q),
        max_da_over_a_pct: worst_step * 100.0,
        tick_ms_median: tick_ms[tick_ms.len() / 2],
        dissolved,
        t0_inflow_per_node: t0_inflow,
    }
}

fn het_apertures(seed: u64) -> impl FnMut() -> f64 {
    let mut lcg = Lcg(seed);
    move || (A0 * (LOGSD * lcg.normal()).exp()).clamp(A_MIN, A_MAX)
}

fn main() {
    if !std::env::args().any(|arg| arg == "--bench") {
        return;
    }

    let prereg = isomesh::experiment!("P-29");
    common::experiment::run(prereg, |run| {
        // ---- detector red/green, before anything is simulated ------------
        let het_seed = 0x9E37_79B9_7F4A_7C15u64;
        let t0: Vec<f64> = {
            let mut sample = het_apertures(het_seed);
            (0..2 * N * N - N).map(|_| sample().ln()).collect()
        };
        let mut t0_sorted = t0.clone();
        let g0 = central_gap(&mut t0_sorted);
        assert!(
            g0 < GAP_THRESHOLD,
            "detector green failed: t=0 lognormal gap {g0} ≥ {GAP_THRESHOLD}"
        );
        let mut shifted = t0.clone();
        let mut lcg = Lcg(7);
        for v in &mut shifted {
            if lcg.next_u64().is_multiple_of(2) {
                *v += 4.0 * LOGSD;
            }
        }
        let gs = central_gap(&mut shifted);
        assert!(
            gs >= GAP_THRESHOLD,
            "detector red failed: synthetic +4σ half-shift gap {gs} < {GAP_THRESHOLD}"
        );
        println!(
            "detector: t0 gap {g0:.5} (< {GAP_THRESHOLD} — green), synthetic gap {gs:.3} (≥ {GAP_THRESHOLD} — red demonstrated)"
        );

        // ---- arm 1: homogeneous with the paper's three seeds --------------
        let mut hom = build_edges(|| A0);
        for &sy in &SEEDS {
            for x in 0..10 {
                let ei = x * N + sy; // horizontal edges are laid out x-major
                hom[ei].a.fill(A0 * (1.0 + 1e-9));
            }
        }
        let hom_result = run_arm(
            &mut hom,
            &Boundary::ConstantHead,
            None,
            MAX_YEARS_HEAD,
            DT_CAP,
        );
        assert!(
            hom_result.breakthrough_years.is_some(),
            "hom_seeded3 never broke through — fixture failure"
        );

        // ---- arm 2: the paper's lognormal heterogeneous net ---------------
        let mut het = build_edges(het_apertures(het_seed));
        let mut t0_gini: Vec<f64> = {
            let heads = solve_heads(&het, &Boundary::ConstantHead);
            het.iter()
                .map(|e| ((heads[e.from] - heads[e.to]) / e.resistance()).abs())
                .collect()
        };
        let g_t0 = gini(&mut t0_gini);
        assert!(g_t0 > 0.0, "t=0 heterogeneous flow Gini must be non-zero");
        let het_result = run_arm(
            &mut het,
            &Boundary::ConstantHead,
            None,
            MAX_YEARS_HEAD,
            DT_CAP,
        );

        // ---- dt control: half the cap, same net, same verdicts required ---
        let mut het_ctl = build_edges(het_apertures(het_seed));
        let ctl_result = run_arm(
            &mut het_ctl,
            &Boundary::ConstantHead,
            None,
            MAX_YEARS_HEAD,
            DT_CAP_CONTROL,
        );
        assert!(
            (ctl_result.max_gap_ln >= GAP_THRESHOLD) == (het_result.max_gap_ln >= GAP_THRESHOLD),
            "dt control changed the bimodal verdict: {} vs {}",
            ctl_result.max_gap_ln,
            het_result.max_gap_ln
        );
        assert!(
            (ctl_result.flux_top10_pct - het_result.flux_top10_pct).abs() <= 2.0,
            "dt control moved C2's share by more than 2 points: {} vs {}",
            ctl_result.flux_top10_pct,
            het_result.flux_top10_pct
        );
        assert!(
            het_result.breakthrough_years.is_some(),
            "heterogeneous never broke through — fixture failure"
        );

        // ---- arm 3: recharge-limited, same net, flux frozen at the head
        // arm's own t=0 inflow, run to matched dissolved volume -------------
        let mut rl = build_edges(het_apertures(het_seed));
        let q = het_result.t0_inflow_per_node;
        let bt_volume = het_result
            .dissolved_at_breakthrough
            .expect("het arm broke through, so its at-breakthrough volume exists");
        let rl_result = run_arm(
            &mut rl,
            &Boundary::FixedFlux(q),
            Some(bt_volume),
            MAX_YEARS_RECHARGE,
            DT_CAP,
        );

        // ---- emit ----------------------------------------------------------
        println!(
            "\n{:>16} {:>6} {:>9} {:>12} {:>10} {:>8} {:>11} {:>9} {:>9} {:>8}",
            "arm",
            "ticks",
            "years",
            "breakthrough",
            "max_gap",
            "bimodal",
            "flux_top10%",
            "gini",
            "max_da/a%",
            "tick_ms"
        );
        for (name, cap, r) in [
            ("hom_seeded3", DT_CAP, &hom_result),
            ("heterogeneous", DT_CAP, &het_result),
            ("het_dt_control", DT_CAP_CONTROL, &ctl_result),
            ("recharge_limited", DT_CAP, &rl_result),
        ] {
            let bimodal = r.max_gap_ln >= GAP_THRESHOLD;
            let bt = r
                .breakthrough_years
                .map_or_else(|| "none".to_string(), |t| format!("{t:.1}"));
            println!(
                "{:>16} {:>6} {:>9.1} {:>12} {:>10.3} {:>8} {:>11.2} {:>9.4} {:>9.3} {:>8.2}",
                name,
                r.ticks,
                r.years,
                bt,
                r.max_gap_ln,
                bimodal,
                r.flux_top10_pct,
                r.gini_flow,
                r.max_da_over_a_pct,
                r.tick_ms_median
            );
            run.record(&[
                ("arm", name.to_string()),
                ("ticks", r.ticks.to_string()),
                ("years", format!("{:.2}", r.years)),
                ("breakthrough_years", bt),
                ("max_gap_ln", format!("{:.4}", r.max_gap_ln)),
                ("bimodal", bimodal.to_string()),
                ("flux_top10_pct", format!("{:.2}", r.flux_top10_pct)),
                ("gini_flow", format!("{:.4}", r.gini_flow)),
                ("max_da_over_a_pct", format!("{:.3}", r.max_da_over_a_pct)),
                ("tick_ms_median", format!("{:.3}", r.tick_ms_median)),
                ("dt_cap", format!("{cap}")),
            ]);
        }

        println!();
        let c1 = hom_result.max_gap_ln >= GAP_THRESHOLD
            && het_result.max_gap_ln >= GAP_THRESHOLD
            && rl_result.max_gap_ln < GAP_THRESHOLD;
        println!(
            "C1 (bimodal split): hom {:.3} / het {:.3} bimodal, recharge {:.3} unimodal -- {}",
            hom_result.max_gap_ln,
            het_result.max_gap_ln,
            rl_result.max_gap_ln,
            if c1 { "HELD" } else { "FALSIFIED" }
        );
        let c2 = hom_result.flux_top10_pct > 90.0 && het_result.flux_top10_pct > 90.0;
        println!(
            "C2 (concentration): top-10% edges carry {:.1}% (hom) / {:.1}% (het) of dissolution flux -- {} (H says > 90%)",
            hom_result.flux_top10_pct,
            het_result.flux_top10_pct,
            if c2 { "HELD" } else { "FALSIFIED" }
        );
        println!(
            "recorded ordering: het breakthrough {:.0} yr vs hom {:.0} yr (paper: 560 vs ~1890) -- {}",
            het_result.breakthrough_years.unwrap_or(f64::NAN),
            hom_result.breakthrough_years.unwrap_or(f64::NAN),
            if het_result.breakthrough_years < hom_result.breakthrough_years {
                "same order"
            } else {
                "ORDER INVERTED"
            }
        );
        println!(
            "recharge arm dissolved {:.3e} mol vs at-breakthrough target {:.3e} mol in {:.0} yr",
            rl_result.dissolved, bt_volume, rl_result.years
        );
    });
}
