//! **P-32 — the modal kill-shot: does one voxel move λ₁ at all?**
//!
//! Ticket: R-033. Pre-registered in the commit before this one.
//!
//! ```bash
//! cargo bench --bench experiment_p32
//! ```
//!
//! Writes `docs/experiments/p-32.csv`.
//!
//! # The design in one line
//!
//! Hexahedral FEM assembled directly on the occupancy grid — the
//! discretisation the sound literature itself uses, no tetrahedralisation —
//! and the smallest eigenpair before and after single-voxel edits. The pitch
//! JND is 0.6% on *frequency*, `f ∝ √λ`, so the registered audibility
//! threshold on λ₁ is **1.2%**.
//!
//! # Why fixed iteration counts are safe here
//!
//! Every reported eigenvalue carries an a-posteriori certificate: with lumped
//! (diagonal) mass, `ε = ‖Kx − λ_R·Mx‖_{M⁻¹}/‖x‖_M` is one pass, and there is
//! an eigenvalue of the pencil within `ε` of `λ_R`. The run asserts
//! `ε ≤ 5e-4·λ_R` — 24× under the decision threshold — on every solve, two
//! deterministic starts must agree within `2ε`, and the certificate is shown
//! able to go red on a deliberately under-converged run before any verdict is
//! read.
//!
//! # Element validation before any assembly is trusted
//!
//! The 24×24 trilinear stiffness (2×2×2 Gauss of `BᵀDB`, `E = 1`, `ν = 0.3`)
//! must be symmetric to 1e-10, annihilate the three translations and three
//! linearised rotations, and — via a Jacobi eigensolve of the free element —
//! hold **exactly six** rigid modes. A full-integration trilinear hex has no
//! spurious zero-energy modes; more than six means the matrix is wrong.
//!
//! # Descoped, deliberately
//!
//! The k ∈ {8..128} mode-count timing sweep the dossier sketched is not here:
//! the registered hypothesis is λ₁-only, deflation is several hundred lines
//! that matter only if the kill-shot survives, and premise-falsifiers-first
//! is this backlog's own ordering. λ comparisons are counts against a
//! threshold, not a timing A/B, so M-197's interleaving rule does not apply.

mod common;

/// Occupancy grid side (cells); nodes are (SIDE+1)³.
const SIDE: usize = 32;
/// Isotropic material: E = 1, ν = 0.3, ρ = 1 — units cancel in Δλ/λ.
const NU: f64 = 0.3;
/// The registered audibility threshold on λ₁ (1.2% = 2 × 0.6% pitch JND).
const AUDIBLE_PCT: f64 = 1.2;
/// The control cavity's registered floor.
const CONTROL_PCT: f64 = 15.0;
/// Certificate bound, relative to λ.
const CERT_REL: f64 = 5e-4;
/// Fixed iteration counts (the certificate, not the counts, is the safety).
const OUTER: usize = 48;
const INNER: usize = 256;

// ---------------------------------------------------------------------------
// Element
// ---------------------------------------------------------------------------

/// Corner offsets in the local node order used throughout: n = dx + 2dy + 4dz.
const CORNERS: [[usize; 3]; 8] = [
    [0, 0, 0],
    [1, 0, 0],
    [0, 1, 0],
    [1, 1, 0],
    [0, 0, 1],
    [1, 0, 1],
    [0, 1, 1],
    [1, 1, 1],
];

/// The trilinear hex stiffness on a unit cube, 2×2×2 Gauss of `BᵀDB`.
#[allow(
    clippy::needless_range_loop,
    reason = "dense matrix kernels index two arrays per loop; iterator forms obscure the algebra"
)]
fn element_stiffness() -> [[f64; 24]; 24] {
    let lambda = NU / ((1.0 + NU) * (1.0 - 2.0 * NU));
    let mu = 1.0 / (2.0 * (1.0 + NU));
    let g = 1.0 / 3.0f64.sqrt();
    let mut ke = [[0.0f64; 24]; 24];
    for gp in 0..8 {
        let xi = [
            g * (2.0 * ((gp) & 1) as f64 - 1.0),
            g * (2.0 * ((gp >> 1) & 1) as f64 - 1.0),
            g * (2.0 * ((gp >> 2) & 1) as f64 - 1.0),
        ];
        // dN/dx for each node: natural derivative × (2/h), h = 1.
        let mut dndx = [[0.0f64; 3]; 8];
        for (n, c) in CORNERS.iter().enumerate() {
            let s = [
                2.0 * c[0] as f64 - 1.0,
                2.0 * c[1] as f64 - 1.0,
                2.0 * c[2] as f64 - 1.0,
            ];
            dndx[n][0] = 2.0 * s[0] * (1.0 + s[1] * xi[1]) * (1.0 + s[2] * xi[2]) / 8.0;
            dndx[n][1] = 2.0 * (1.0 + s[0] * xi[0]) * s[1] * (1.0 + s[2] * xi[2]) / 8.0;
            dndx[n][2] = 2.0 * (1.0 + s[0] * xi[0]) * (1.0 + s[1] * xi[1]) * s[2] / 8.0;
        }
        // B is 6×24; accumulate BᵀDB directly. Strain order:
        // [εxx, εyy, εzz, γxy, γyz, γzx].
        let mut b = [[0.0f64; 24]; 6];
        for n in 0..8 {
            let (dx, dy, dz) = (dndx[n][0], dndx[n][1], dndx[n][2]);
            b[0][3 * n] = dx;
            b[1][3 * n + 1] = dy;
            b[2][3 * n + 2] = dz;
            b[3][3 * n] = dy;
            b[3][3 * n + 1] = dx;
            b[4][3 * n + 1] = dz;
            b[4][3 * n + 2] = dy;
            b[5][3 * n] = dz;
            b[5][3 * n + 2] = dx;
        }
        let mut d = [[0.0f64; 6]; 6];
        for i in 0..3 {
            for j in 0..3 {
                d[i][j] = if i == j { lambda + 2.0 * mu } else { lambda };
            }
            d[i + 3][i + 3] = mu;
        }
        // Ke += Bᵀ D B · (detJ = 1/8, weight 1).
        for i in 0..24 {
            for r in 0..6 {
                if b[r][i] == 0.0 {
                    continue;
                }
                for s in 0..6 {
                    let drs = d[r][s];
                    if drs == 0.0 {
                        continue;
                    }
                    let coeff = b[r][i] * drs / 8.0;
                    for (j, ke_ij) in ke[i].iter_mut().enumerate() {
                        *ke_ij += coeff * b[s][j];
                    }
                }
            }
        }
    }
    ke
}

/// Cyclic Jacobi eigenvalues of a symmetric 24×24 — for the free-element
/// rigid-mode census only.
#[allow(
    clippy::needless_range_loop,
    reason = "dense matrix kernels index two arrays per loop; iterator forms obscure the algebra"
)]
fn jacobi_eigenvalues(mut a: [[f64; 24]; 24]) -> [f64; 24] {
    for _sweep in 0..100 {
        let mut off = 0.0;
        for i in 0..24 {
            for j in i + 1..24 {
                off += a[i][j] * a[i][j];
            }
        }
        if off < 1e-24 {
            break;
        }
        for p in 0..24 {
            for q in p + 1..24 {
                if a[p][q].abs() < 1e-15 {
                    continue;
                }
                let theta = (a[q][q] - a[p][p]) / (2.0 * a[p][q]);
                let t = theta.signum() / (theta.abs() + (theta * theta + 1.0).sqrt());
                let c = 1.0 / (t * t + 1.0).sqrt();
                let s = t * c;
                for k in 0..24 {
                    let akp = a[k][p];
                    let akq = a[k][q];
                    a[k][p] = c * akp - s * akq;
                    a[k][q] = s * akp + c * akq;
                }
                for k in 0..24 {
                    let apk = a[p][k];
                    let aqk = a[q][k];
                    a[p][k] = c * apk - s * aqk;
                    a[q][k] = s * apk + c * aqk;
                }
            }
        }
    }
    let mut eig = [0.0f64; 24];
    for (i, e) in eig.iter_mut().enumerate() {
        *e = a[i][i];
    }
    eig.sort_unstable_by(f64::total_cmp);
    eig
}

#[allow(
    clippy::needless_range_loop,
    reason = "dense matrix kernels index two arrays per loop; iterator forms obscure the algebra"
)]
fn validate_element(ke: &[[f64; 24]; 24]) {
    for i in 0..24 {
        for j in 0..24 {
            assert!(
                (ke[i][j] - ke[j][i]).abs() < 1e-10,
                "element stiffness not symmetric at ({i},{j})"
            );
        }
    }
    // Null vectors: three translations, three linearised rotations u = ω × x.
    let mut nulls: Vec<[f64; 24]> = Vec::new();
    for axis in 0..3 {
        let mut v = [0.0f64; 24];
        for n in 0..8 {
            v[3 * n + axis] = 1.0;
        }
        nulls.push(v);
    }
    for axis in 0..3 {
        let mut v = [0.0f64; 24];
        let w = match axis {
            0 => [1.0, 0.0, 0.0],
            1 => [0.0, 1.0, 0.0],
            _ => [0.0, 0.0, 1.0],
        };
        for (n, c) in CORNERS.iter().enumerate() {
            let x = [c[0] as f64, c[1] as f64, c[2] as f64];
            v[3 * n] = w[1] * x[2] - w[2] * x[1];
            v[3 * n + 1] = w[2] * x[0] - w[0] * x[2];
            v[3 * n + 2] = w[0] * x[1] - w[1] * x[0];
        }
        nulls.push(v);
    }
    for (which, v) in nulls.iter().enumerate() {
        for i in 0..24 {
            let mut s = 0.0;
            for j in 0..24 {
                s += ke[i][j] * v[j];
            }
            assert!(
                s.abs() < 1e-10,
                "element stiffness does not annihilate null vector {which} (row {i}: {s})"
            );
        }
    }
    let eig = jacobi_eigenvalues(*ke);
    let max = eig[23];
    let zeros = eig.iter().filter(|e| e.abs() < 1e-9 * max).count();
    assert!(
        zeros == 6,
        "free element holds {zeros} rigid modes, not 6 — the matrix is wrong"
    );
}

// ---------------------------------------------------------------------------
// Fixture
// ---------------------------------------------------------------------------

/// The carved pillar: a tapering, drifting ellipse column plus a two-cell
/// thick web (the thin feature) sticking out of its +x flank.
fn base_occupancy() -> Vec<bool> {
    let mut occ = vec![false; SIDE * SIDE * SIDE];
    for z in 0..SIDE {
        let zf = z as f64;
        let cx = 12.0 + 0.18 * zf;
        let cy = 14.0 + 0.07 * zf;
        let rx = 9.0 - 0.16 * zf;
        let ry = 7.0 - 0.09 * zf;
        for x in 0..SIDE {
            for y in 0..SIDE {
                let dx = (x as f64 + 0.5 - cx) / rx;
                let dy = (y as f64 + 0.5 - cy) / ry;
                if dx * dx + dy * dy <= 1.0 {
                    occ[cell(x, y, z)] = true;
                }
            }
        }
        // The web: z ∈ [8, 20], two cells thick in y, reaching +x.
        if (8..=20).contains(&z) {
            let x0 = (cx + rx) as usize;
            for x in x0..(x0 + 7).min(SIDE) {
                for wy in 0..2usize {
                    let y = (cy as usize + wy).min(SIDE - 1);
                    occ[cell(x, y, z)] = true;
                }
            }
        }
    }
    occ
}

fn cell(x: usize, y: usize, z: usize) -> usize {
    (z * SIDE + y) * SIDE + x
}

/// Cells of the web (recomputed by the same rule, for edit classification).
fn web_cells(occ: &[bool]) -> Vec<usize> {
    let mut web = Vec::new();
    for z in 8..=20usize {
        let zf = z as f64;
        let cx = 12.0 + 0.18 * zf;
        let cy = 14.0 + 0.07 * zf;
        let rx = 9.0 - 0.16 * zf;
        let x0 = (cx + rx) as usize;
        for x in x0..(x0 + 7).min(SIDE) {
            for wy in 0..2usize {
                let y = (cy as usize + wy).min(SIDE - 1);
                if occ[cell(x, y, z)] {
                    web.push(cell(x, y, z));
                }
            }
        }
    }
    web
}

/// Occupied cells must form one component containing the base layer.
fn assert_connected(occ: &[bool]) {
    let total = occ.iter().filter(|&&o| o).count();
    let mut seen = vec![false; occ.len()];
    let mut queue: Vec<usize> = (0..SIDE * SIDE)
        .map(|i| {
            let (x, y) = (i % SIDE, i / SIDE);
            cell(x, y, 0)
        })
        .filter(|&c| occ[c])
        .collect();
    assert!(!queue.is_empty(), "no occupied cells on the base layer");
    for &q in &queue {
        seen[q] = true;
    }
    let mut reached = queue.len();
    while let Some(c) = queue.pop() {
        let x = c % SIDE;
        let y = (c / SIDE) % SIDE;
        let z = c / (SIDE * SIDE);
        let mut push = |xx: i64, yy: i64, zz: i64| {
            if (0..SIDE as i64).contains(&xx)
                && (0..SIDE as i64).contains(&yy)
                && (0..SIDE as i64).contains(&zz)
            {
                let n = cell(xx as usize, yy as usize, zz as usize);
                if occ[n] && !seen[n] {
                    seen[n] = true;
                    reached += 1;
                    queue.push(n);
                }
            }
        };
        push(x as i64 + 1, y as i64, z as i64);
        push(x as i64 - 1, y as i64, z as i64);
        push(x as i64, y as i64 + 1, z as i64);
        push(x as i64, y as i64 - 1, z as i64);
        push(x as i64, y as i64, z as i64 + 1);
        push(x as i64, y as i64, z as i64 - 1);
    }
    assert!(
        reached == total,
        "occupancy is not one base-anchored component: {reached} of {total}"
    );
}

// ---------------------------------------------------------------------------
// Assembly-free operators
// ---------------------------------------------------------------------------

struct Mesh {
    cells: Vec<usize>,
    /// Per cell: the 24 global DOF indices (usize::MAX = fixed or absent).
    cell_dofs: Vec<[usize; 24]>,
    ndof: usize,
    mass: Vec<f64>,
    kdiag: Vec<f64>,
}

fn build_mesh(occ: &[bool], ke: &[[f64; 24]; 24]) -> Mesh {
    let np = SIDE + 1;
    let node = |x: usize, y: usize, z: usize| (z * np + y) * np + x;
    let mut dof_of = vec![usize::MAX; np * np * np];
    let mut ndof = 0usize;
    let cells: Vec<usize> = (0..occ.len()).filter(|&c| occ[c]).collect();
    for &c in &cells {
        let x = c % SIDE;
        let y = (c / SIDE) % SIDE;
        let z = c / (SIDE * SIDE);
        for corner in &CORNERS {
            let (nx, ny, nz) = (x + corner[0], y + corner[1], z + corner[2]);
            if nz == 0 {
                continue; // base layer: fixed
            }
            let n = node(nx, ny, nz);
            if dof_of[n] == usize::MAX {
                dof_of[n] = ndof;
                ndof += 3;
            }
        }
    }
    let mut mass = vec![0.0f64; ndof];
    let mut kdiag = vec![0.0f64; ndof];
    let mut cell_dofs = Vec::with_capacity(cells.len());
    for &c in &cells {
        let x = c % SIDE;
        let y = (c / SIDE) % SIDE;
        let z = c / (SIDE * SIDE);
        let mut dofs = [usize::MAX; 24];
        for (i, corner) in CORNERS.iter().enumerate() {
            let (nx, ny, nz) = (x + corner[0], y + corner[1], z + corner[2]);
            if nz == 0 {
                continue;
            }
            let d = dof_of[node(nx, ny, nz)];
            for a in 0..3 {
                dofs[3 * i + a] = d + a;
                mass[d + a] += 0.125; // ρ h³ / 8, ρ = h = 1
                kdiag[d + a] += ke[3 * i + a][3 * i + a];
            }
        }
        cell_dofs.push(dofs);
    }
    Mesh {
        cells,
        cell_dofs,
        ndof,
        mass,
        kdiag,
    }
}

fn apply_k(mesh: &Mesh, ke: &[[f64; 24]; 24], v: &[f64], out: &mut [f64]) {
    out.fill(0.0);
    for dofs in &mesh.cell_dofs {
        let mut local = [0.0f64; 24];
        for (l, &d) in dofs.iter().enumerate() {
            if d != usize::MAX {
                local[l] = v[d];
            }
        }
        for (l, &d) in dofs.iter().enumerate() {
            if d == usize::MAX {
                continue;
            }
            let row = &ke[l];
            let mut s = 0.0;
            for (r, lv) in row.iter().zip(local.iter()) {
                s += r * lv;
            }
            out[d] += s;
        }
    }
}

/// Jacobi-preconditioned CG on `K y = b`, fixed iteration count, warm start.
fn pcg(mesh: &Mesh, ke: &[[f64; 24]; 24], b: &[f64], y: &mut [f64], iters: usize) {
    let n = mesh.ndof;
    let mut r = vec![0.0f64; n];
    let mut z = vec![0.0f64; n];
    let mut p = vec![0.0f64; n];
    let mut kp = vec![0.0f64; n];
    apply_k(mesh, ke, y, &mut r);
    for i in 0..n {
        r[i] = b[i] - r[i];
        z[i] = r[i] / mesh.kdiag[i];
        p[i] = z[i];
    }
    let mut rz: f64 = r.iter().zip(&z).map(|(a, b)| a * b).sum();
    for _ in 0..iters {
        if rz.abs() < 1e-300 {
            break;
        }
        apply_k(mesh, ke, &p, &mut kp);
        let pkp: f64 = p.iter().zip(&kp).map(|(a, b)| a * b).sum();
        if pkp <= 0.0 {
            break;
        }
        let alpha = rz / pkp;
        for i in 0..n {
            y[i] += alpha * p[i];
            r[i] -= alpha * kp[i];
        }
        let mut rz_new = 0.0;
        for i in 0..n {
            z[i] = r[i] / mesh.kdiag[i];
            rz_new += r[i] * z[i];
        }
        let beta = rz_new / rz;
        rz = rz_new;
        for i in 0..n {
            p[i] = z[i] + beta * p[i];
        }
    }
}

struct EigResult {
    lambda: f64,
    cert_rel: f64,
}

/// Smallest eigenpair of `K x = λ M x` by inverse power iteration; the
/// certificate is what makes the fixed counts trustworthy.
fn lambda1(mesh: &Mesh, ke: &[[f64; 24]; 24], seed: u64, outer: usize) -> EigResult {
    let n = mesh.ndof;
    let mut x = vec![0.0f64; n];
    if seed == 0 {
        x.fill(1.0);
    } else {
        let mut state = seed;
        for v in &mut x {
            state = state
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            *v = (((state >> 11) as f64) / (1u64 << 53) as f64) - 0.5;
        }
    }
    let mnorm = |v: &[f64]| -> f64 {
        v.iter()
            .zip(&mesh.mass)
            .map(|(a, m)| a * a * m)
            .sum::<f64>()
            .sqrt()
    };
    let nm = mnorm(&x);
    for v in &mut x {
        *v /= nm;
    }
    let mut y = vec![0.0f64; n];
    let mut b = vec![0.0f64; n];
    for _ in 0..outer {
        for i in 0..n {
            b[i] = mesh.mass[i] * x[i];
        }
        pcg(mesh, ke, &b, &mut y, INNER);
        let nm = mnorm(&y);
        assert!(nm > 0.0, "inverse iteration collapsed to zero");
        for i in 0..n {
            x[i] = y[i] / nm;
        }
    }
    let mut kx = vec![0.0f64; n];
    apply_k(mesh, ke, &x, &mut kx);
    let xkx: f64 = x.iter().zip(&kx).map(|(a, b)| a * b).sum();
    let xmx: f64 = x
        .iter()
        .zip(&mesh.mass)
        .map(|(a, m)| a * a * m)
        .sum::<f64>();
    let lambda = xkx / xmx;
    let mut cert = 0.0f64;
    for i in 0..n {
        let r = kx[i] - lambda * mesh.mass[i] * x[i];
        cert += r * r / mesh.mass[i];
    }
    let cert = cert.sqrt() / xmx.sqrt();
    EigResult {
        lambda,
        cert_rel: cert / lambda,
    }
}

/// λ₁ with the full protocol: two deterministic starts, agreement within
/// twice the certificate, both certificates under the registered bound.
fn lambda1_certified(mesh: &Mesh, ke: &[[f64; 24]; 24]) -> EigResult {
    let a = lambda1(mesh, ke, 0, OUTER);
    let b = lambda1(mesh, ke, 0x9E37_79B9_7F4A_7C15, OUTER);
    for (name, r) in [("ones", &a), ("lcg", &b)] {
        assert!(
            r.cert_rel <= CERT_REL,
            "certificate {:.3e} above the registered {CERT_REL:.0e} bound ({name} start) — the \
             fixed counts were too few and the run refuses to report",
            r.cert_rel
        );
    }
    let tol = 2.0 * (a.cert_rel.max(b.cert_rel)) * a.lambda;
    assert!(
        (a.lambda - b.lambda).abs() <= tol.max(1e-12),
        "two starts disagree: {} vs {} (tol {tol})",
        a.lambda,
        b.lambda
    );
    a
}

fn main() {
    if !std::env::args().any(|arg| arg == "--bench") {
        return;
    }

    let prereg = isomesh::experiment!("P-32");
    common::experiment::run(prereg, |run| {
        let ke = element_stiffness();
        validate_element(&ke);
        println!("element: symmetric, nulls annihilated, exactly 6 rigid modes — validated");

        let occ = base_occupancy();
        assert_connected(&occ);
        let web: Vec<usize> = web_cells(&occ);
        assert!(!web.is_empty(), "the web is missing from the fixture");

        let mesh = build_mesh(&occ, &ke);
        let cells0 = mesh.cells.len();

        // Certificate inversion: a deliberately under-converged run must go
        // red, before any verdict is read.
        let bad = lambda1(&mesh, &ke, 0, 2);
        assert!(
            bad.cert_rel > 10.0 * CERT_REL,
            "certificate failed to flag a 2-outer run ({:.3e}) — it cannot be trusted green",
            bad.cert_rel
        );
        println!(
            "certificate inversion: 2-outer run flagged at {:.2e} (bound {CERT_REL:.0e}) — red demonstrated",
            bad.cert_rel
        );

        let base = lambda1_certified(&mesh, &ke);
        println!(
            "baseline: {} cells, {} dof, λ₁ = {:.6e} (cert {:.2e})",
            cells0, mesh.ndof, base.lambda, base.cert_rel
        );

        // Edit sets. Interior: all six face neighbours occupied, ≥ 4 cells
        // from any web cell (Chebyshev), z ≥ 4. Web-adjacent: within 1 cell
        // of a web cell (web cells themselves qualify).
        let is_web_near = |c: usize, dist: i64| -> bool {
            let x = (c % SIDE) as i64;
            let y = ((c / SIDE) % SIDE) as i64;
            let z = (c / (SIDE * SIDE)) as i64;
            web.iter().any(|&w| {
                let wx = (w % SIDE) as i64;
                let wy = ((w / SIDE) % SIDE) as i64;
                let wz = (w / (SIDE * SIDE)) as i64;
                (x - wx).abs().max((y - wy).abs()).max((z - wz).abs()) <= dist
            })
        };
        let interior: Vec<usize> = mesh
            .cells
            .iter()
            .copied()
            .filter(|&c| {
                let x = c % SIDE;
                let y = (c / SIDE) % SIDE;
                let z = c / (SIDE * SIDE);
                z >= 4
                    && x >= 1
                    && y >= 1
                    && x + 1 < SIDE
                    && y + 1 < SIDE
                    && z + 1 < SIDE
                    && occ[cell(x + 1, y, z)]
                    && occ[cell(x - 1, y, z)]
                    && occ[cell(x, y + 1, z)]
                    && occ[cell(x, y - 1, z)]
                    && occ[cell(x, y, z + 1)]
                    && occ[cell(x, y, z - 1)]
                    && !is_web_near(c, 4)
            })
            .collect();
        assert!(
            interior.len() >= 8,
            "only {} strictly-interior candidates",
            interior.len()
        );
        let near_web: Vec<usize> = mesh
            .cells
            .iter()
            .copied()
            .filter(|&c| is_web_near(c, 1))
            .collect();
        assert!(near_web.len() >= 4, "web neighbourhood too small");

        // Deterministic picks: evenly spaced through each candidate list.
        let pick = |list: &[usize], want: usize| -> Vec<usize> {
            (0..want)
                .map(|i| list[(i * list.len()) / want + (list.len() / (2 * want))])
                .collect()
        };
        let interior_edits = pick(&interior, 8);
        let web_edits = pick(&near_web, 4);

        let mut rows: Vec<(String, usize, usize, f64, f64)> = Vec::new();
        rows.push((
            "baseline".to_string(),
            cells0,
            mesh.ndof,
            base.lambda,
            base.cert_rel,
        ));

        let mut run_edit = |name: String, removed: &[usize]| -> f64 {
            let mut occ2 = occ.clone();
            for &c in removed {
                occ2[c] = false;
            }
            assert_connected(&occ2);
            let mesh2 = build_mesh(&occ2, &ke);
            let r = lambda1_certified(&mesh2, &ke);
            rows.push((name, mesh2.cells.len(), mesh2.ndof, r.lambda, r.cert_rel));
            100.0 * (r.lambda - base.lambda).abs() / base.lambda
        };

        let mut interior_deltas = Vec::new();
        for (i, &c) in interior_edits.iter().enumerate() {
            interior_deltas.push(run_edit(format!("interior_{i}"), &[c]));
        }
        let mut web_deltas = Vec::new();
        for (i, &c) in web_edits.iter().enumerate() {
            web_deltas.push(run_edit(format!("near_web_{i}"), &[c]));
        }

        // Control: a cavity of ≥ 20% of the cells, centred mid-pillar,
        // radius grown deterministically until the volume is reached.
        let centre = [14.0, 14.5, 12.0];
        let mut radius = 3.0f64;
        let cavity: Vec<usize> = loop {
            let ball: Vec<usize> = mesh
                .cells
                .iter()
                .copied()
                .filter(|&c| {
                    let x = (c % SIDE) as f64 + 0.5;
                    let y = ((c / SIDE) % SIDE) as f64 + 0.5;
                    let z = (c / (SIDE * SIDE)) as f64 + 0.5;
                    let d =
                        (x - centre[0]).powi(2) + (y - centre[1]).powi(2) + (z - centre[2]).powi(2);
                    d.sqrt() <= radius
                })
                .collect();
            if ball.len() as f64 >= 0.20 * cells0 as f64 {
                break ball;
            }
            radius += 0.5;
            assert!(radius < SIDE as f64, "cavity never reached 20%");
        };
        let control_delta = run_edit("cavity_20pct".to_string(), &cavity);
        assert!(
            control_delta > CONTROL_PCT,
            "control cavity moved λ₁ by only {control_delta:.2}% — the instrument cannot see \
             shifts and the null verdicts are void"
        );

        // ---- emit ----------------------------------------------------------
        println!(
            "\n{:>14} {:>7} {:>8} {:>14} {:>11} {:>9}",
            "edit", "cells", "dof", "lambda1", "delta%", "cert"
        );
        for (name, cells, dof, lambda, cert) in &rows {
            let delta = 100.0 * (lambda - base.lambda).abs() / base.lambda;
            let audible = delta >= AUDIBLE_PCT;
            println!(
                "{:>14} {:>7} {:>8} {:>14.6e} {:>11.4} {:>9.2e}",
                name, cells, dof, lambda, delta, cert
            );
            run.record(&[
                ("edit", name.clone()),
                ("cells", cells.to_string()),
                ("dof", dof.to_string()),
                ("lambda1_base", format!("{:.6e}", base.lambda)),
                ("lambda1_edited", format!("{lambda:.6e}")),
                ("delta_pct", format!("{delta:.4}")),
                ("audible", audible.to_string()),
                ("certificate_rel", format!("{cert:.3e}")),
            ]);
        }

        println!();
        let c1_worst = interior_deltas.iter().copied().fold(0.0f64, f64::max);
        let c1 = c1_worst < AUDIBLE_PCT;
        println!(
            "C1 (interior digs inaudible): worst {:.4}% of λ₁ against the {AUDIBLE_PCT}% threshold \
             -- {}",
            c1_worst,
            if c1 {
                "HELD"
            } else {
                "FALSIFIED — per-edit modal audio earns its ticket"
            }
        );
        let c2_best = web_deltas.iter().copied().fold(0.0f64, f64::max);
        let c2 = c2_best >= AUDIBLE_PCT;
        println!(
            "C2 (near-web digs audible): best {:.4}% -- {}",
            c2_best,
            if c2 {
                "HELD"
            } else {
                "FALSIFIED — even thin-feature edits are inaudible, the direction closes entirely"
            }
        );
        println!(
            "C3 (control): 20% cavity moved λ₁ by {control_delta:.1}% (> {CONTROL_PCT}%) -- the \
             null verdicts above are readings of the edits, not of a numb instrument"
        );
    });
}
