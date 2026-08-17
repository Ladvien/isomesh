//! **R-035a — the prefactored-surface-operator substrate, verified before it
//! is timed.**
//!
//! ```bash
//! cargo bench --bench r035a_substrate
//! ```
//!
//! Writes `docs/measurements/r035a-substrate.csv`. Comparative and
//! deliberately without a `P-` id: correctness against a dense reference is
//! not a hypothesis (the M-322 precedent).
//!
//! Two stages, both loud:
//!
//! 1. **Dense reference.** On a small Surface Nets mesh, the sparse
//!    up-looking LDLᵀ (nested-dissection ordered) must match a dense LDLᵀ of
//!    the permuted heat operator entrywise at 1e-10, and their solves must
//!    agree.
//! 2. **At scale.** On the 64³ chunk mesh, `‖L·D·Lᵀ − A‖_F/‖A‖_F ≤ 100·n·u`,
//!    plus a partial-refactor identity check: refactoring the ancestor
//!    closure of an arbitrary seed set over UNCHANGED values must reproduce
//!    the factor bit-for-bit — the update path shown to be a no-op exactly
//!    when nothing changed.

mod common;

use common::heat;
use isomesh::fields::capped_gyroid;
use isomesh::surface_nets::SurfaceNets;
use isomesh::{MeshBuffer, RuntimeShape3};

fn mesh_at(samples: u32) -> MeshBuffer<f64> {
    let field = capped_gyroid::<f64>();
    let shape = RuntimeShape3::new([samples; 3]).expect("grid fits");
    let h = 4.0 / f64::from(samples - 1);
    let mut sn = SurfaceNets::<f64>::new();
    let mut out = MeshBuffer::<f64>::new();
    sn.extract(&field, &shape, [-2.0, -2.0, -2.0], h, &mut out)
        .expect("extraction");
    out
}

#[allow(
    clippy::needless_range_loop,
    reason = "dense matrix kernels index two arrays per loop; iterator forms obscure the algebra"
)]
fn dense_ldl(a: &heat::SymUpper) -> (Vec<f64>, Vec<f64>) {
    let n = a.n;
    let mut m = vec![0.0f64; n * n];
    for (k, col) in a.cols.iter().enumerate() {
        for &(i, v) in col {
            m[i * n + k] = v;
            m[k * n + i] = v;
        }
    }
    let mut l = vec![0.0f64; n * n];
    let mut d = vec![0.0f64; n];
    for j in 0..n {
        let mut dj = m[j * n + j];
        for k in 0..j {
            dj -= l[j * n + k] * l[j * n + k] * d[k];
        }
        assert!(dj > 0.0, "dense LDL pivot {j} not positive");
        d[j] = dj;
        l[j * n + j] = 1.0;
        for i in j + 1..n {
            let mut s = m[i * n + j];
            for k in 0..j {
                s -= l[i * n + k] * l[j * n + k] * d[k];
            }
            l[i * n + j] = s / dj;
        }
    }
    (l, d)
}

fn main() {
    if !std::env::args().any(|arg| arg == "--bench") {
        return;
    }

    // ---- stage 1: dense reference on a small mesh -------------------------
    let small = mesh_at(13);
    let (a_small, _) = heat::heat_operator(&small);
    let perm = heat::nested_dissection(&a_small, 4);
    let ap = heat::permute(&a_small, &perm);
    let f = heat::ldl_factor(&ap);
    let (dl, dd) = dense_ldl(&ap);
    let n = ap.n;
    let mut worst = 0.0f64;
    for j in 0..n {
        worst = worst.max((f.d[j] - dd[j]).abs() / dd[j].abs());
        for (&r, &v) in f.li[j].iter().zip(f.lx[j].iter()) {
            worst = worst.max((v - dl[r * n + j]).abs());
        }
        // And the reverse direction: dense entries the sparse factor missed.
        for i in j + 1..n {
            let dv = dl[i * n + j];
            if dv.abs() > 1e-12 {
                let sv = match f.li[j].binary_search(&i) {
                    Ok(pos) => f.lx[j][pos],
                    Err(_) => 0.0,
                };
                worst = worst.max((dv - sv).abs());
            }
        }
    }
    assert!(
        worst <= 1e-10,
        "sparse vs dense factor diverges: worst {worst:.3e}"
    );
    let rhs: Vec<f64> = (0..n).map(|i| ((i % 7) as f64) - 3.0).collect();
    let xs = heat::solve(&f, &rhs);
    // Dense solve through the dense factor.
    let mut xd = rhs.clone();
    for i in 0..n {
        for k in 0..i {
            xd[i] -= dl[i * n + k] * xd[k];
        }
    }
    for i in 0..n {
        xd[i] /= dd[i];
    }
    for i in (0..n).rev() {
        for k in i + 1..n {
            xd[i] -= dl[k * n + i] * xd[k];
        }
    }
    let mut solve_diff = 0.0f64;
    for i in 0..n {
        solve_diff = solve_diff.max((xs[i] - xd[i]).abs());
    }
    assert!(solve_diff <= 1e-8, "solve disagreement {solve_diff:.3e}");
    println!(
        "stage 1: {} vertices — sparse matches dense at {worst:.2e}, solves agree at {solve_diff:.2e}",
        small.positions.len()
    );

    // ---- stage 2: the 64³ chunk --------------------------------------------
    let big = mesh_at(65);
    let (a_big, h_bar) = heat::heat_operator(&big);
    let perm = heat::nested_dissection(&a_big, 8);
    let apb = heat::permute(&a_big, &perm);
    let start = std::time::Instant::now();
    let fb = heat::ldl_factor(&apb);
    let factor_ms = start.elapsed().as_secs_f64() * 1e3;
    let nnz: usize = fb.li.iter().map(Vec::len).sum();
    let res = heat::residual_fro(&apb, &fb);
    let bound = 100.0 * apb.n as f64 * 2.0f64.powi(-53);
    assert!(
        res <= bound,
        "factor residual {res:.3e} above the derived bound {bound:.3e}"
    );

    // Partial-refactor identity: unchanged values, arbitrary seeds → the
    // factor must reproduce itself bit-for-bit.
    let seeds: Vec<usize> = (0..apb.n).step_by(97).collect();
    let mut mask = heat::ancestor_closure(&fb.parent, &seeds);
    let masked = mask.iter().filter(|&&m| m).count();
    let mut fb2 = heat::LdlFactor {
        n: fb.n,
        parent: fb.parent.clone(),
        li: fb.li.clone(),
        lx: fb.lx.clone(),
        d: fb.d.clone(),
        flops: 0,
    };
    heat::refactor_rows(&apb, &mut fb2, &mut mask);
    let mut identical = fb2.d == fb.d;
    for j in 0..fb.n {
        identical &= fb2.lx[j] == fb.lx[j] && fb2.li[j] == fb.li[j];
    }
    assert!(
        identical,
        "partial refactor over unchanged values failed to reproduce the factor"
    );

    println!(
        "stage 2: {} vertices, {} triangles, h̄ = {h_bar:.4}, nnz(L) = {nnz}, factor {factor_ms:.1} ms, \
         residual {res:.2e} (bound {bound:.2e}), no-op partial refactor over {masked} closure rows is \
         bit-identical",
        big.positions.len(),
        big.indices.len() / 3
    );

    let csv = format!(
        "# r035a substrate verification\nvertices,triangles,nnz_l,factor_ms,residual_fro,bound,closure_rows\n{},{},{},{:.3},{:.3e},{:.3e},{}\n",
        big.positions.len(),
        big.indices.len() / 3,
        nnz,
        factor_ms,
        res,
        bound,
        masked
    );
    let path = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../docs/measurements/r035a-substrate.csv");
    match std::fs::write(&path, &csv) {
        Ok(()) => println!("wrote {}", path.display()),
        Err(e) => println!("::error:: {}: {e}", path.display()),
    }
}
