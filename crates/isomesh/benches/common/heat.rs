//! The prefactored-surface-operator substrate (R-035a), shared by the
//! verification bench and the update experiment (R-035b).
//!
//! Pieces, in dependency order:
//!
//! - triangle mesh → cotan stiffness `W`, lumped mass `M`, heat operator
//!   `A = M + t·W` with `t = h̄²` (Crane, Weischedel & Wardetzky,
//!   `10.1145/2516971.2516977`, READ — *"the one free parameter is t = h²"*);
//! - **nested dissection by BFS bisection** for the ordering — deliberately
//!   not RCM: a banded ordering's elimination tree is essentially a path, so
//!   any factor repair touches nearly every column and R-035b's verdict would
//!   be an ordering artifact;
//! - elimination tree + up-looking sparse `L·D·Lᵀ` (the compact Davis
//!   `LDL` algorithm), columns of `L` stored sorted by row;
//! - `refactor_partial`: re-run the numeric pass only for rows in the
//!   elimination-tree **ancestor closure** of the changed columns. Sufficient
//!   by the etree's own reachability property: `A(k, j) ≠ 0, j < k` puts `k`
//!   on `j`'s ancestor path, and every `L`-row a recomputed row reads is
//!   itself in the closure. No downdates, nothing to lose definiteness over.
//!
//! Flops are counted with plain `u64`s so R-035b can gate its ratio on work
//! as well as time (a large time ratio over a small flop ratio is an
//! implementation artifact, not a result).

use isomesh::MeshBuffer;

/// Symmetric sparse matrix in upper-column form: column `k` holds entries
/// `(i, v)` with `i ≤ k`, sorted by `i`. This is exactly what the up-looking
/// factorization consumes.
pub(crate) struct SymUpper {
    pub n: usize,
    pub cols: Vec<Vec<(usize, f64)>>,
}

/// The factor `L·D·Lᵀ`: unit-lower `L` by columns (strictly-lower entries,
/// sorted by row), positive `D`.
pub(crate) struct LdlFactor {
    pub n: usize,
    pub parent: Vec<i64>,
    pub li: Vec<Vec<usize>>,
    pub lx: Vec<Vec<f64>>,
    pub d: Vec<f64>,
    pub flops: u64,
}

/// Cotan heat operator on a triangle mesh, plus the mean edge length that
/// sets `t`. Entries are (i, j, w) with `i ≤ j` merged.
pub(crate) fn heat_operator(mesh: &MeshBuffer<f64>) -> (SymUpper, f64) {
    let n = mesh.positions.len();
    let mut diag = vec![0.0f64; n];
    let mut mass = vec![0.0f64; n];
    let mut off: std::collections::BTreeMap<(usize, usize), f64> =
        std::collections::BTreeMap::new();
    let mut edge_len_sum = 0.0f64;
    let mut edge_count = 0u64;
    let pos = |i: usize| mesh.positions[i];
    let sub = |a: [f64; 3], b: [f64; 3]| [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
    let dot = |a: [f64; 3], b: [f64; 3]| a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
    let cross = |a: [f64; 3], b: [f64; 3]| {
        [
            a[1] * b[2] - a[2] * b[1],
            a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0],
        ]
    };
    for tri in mesh.indices.chunks_exact(3) {
        let (a, b, c) = (tri[0] as usize, tri[1] as usize, tri[2] as usize);
        let (pa, pb, pc) = (pos(a), pos(b), pos(c));
        let area2 = {
            let cr = cross(sub(pb, pa), sub(pc, pa));
            dot(cr, cr).sqrt()
        };
        if area2 <= 0.0 {
            continue; // a degenerate sliver contributes nothing it can own
        }
        for i in 0..3 {
            mass[tri[i] as usize] += area2 / 6.0; // (area = area2/2), /3 lumped
        }
        // cot at each corner weights the OPPOSITE edge.
        let corners = [(a, b, c), (b, c, a), (c, a, b)];
        for &(apex, u, v) in &corners {
            let e1 = sub(pos(u), pos(apex));
            let e2 = sub(pos(v), pos(apex));
            let cr = cross(e1, e2);
            let denom = dot(cr, cr).sqrt();
            if denom <= 0.0 {
                continue;
            }
            let cot = dot(e1, e2) / denom;
            let w = cot / 2.0;
            let (lo, hi) = if u < v { (u, v) } else { (v, u) };
            *off.entry((lo, hi)).or_insert(0.0) -= w;
            diag[u] += w;
            diag[v] += w;
        }
        for &(u, v) in &[(a, b), (b, c), (c, a)] {
            let d = sub(pos(u), pos(v));
            edge_len_sum += dot(d, d).sqrt();
            edge_count += 1;
        }
    }
    let h_bar = edge_len_sum / edge_count as f64;
    let t = h_bar * h_bar;
    let mut cols: Vec<Vec<(usize, f64)>> = vec![Vec::new(); n];
    for (i, (&m, &w)) in mass.iter().zip(diag.iter()).enumerate() {
        cols[i].push((i, m + t * w));
    }
    for (&(lo, hi), &w) in &off {
        cols[hi].push((lo, t * w));
    }
    for col in &mut cols {
        col.sort_unstable_by_key(|&(i, _)| i);
    }
    (SymUpper { n, cols }, h_bar)
}

/// Permute a `SymUpper`: `perm[new] = old`.
pub(crate) fn permute(a: &SymUpper, perm: &[usize]) -> SymUpper {
    let n = a.n;
    let mut inv = vec![0usize; n];
    for (new, &old) in perm.iter().enumerate() {
        inv[old] = new;
    }
    let mut cols: Vec<Vec<(usize, f64)>> = vec![Vec::new(); n];
    for (k, col) in a.cols.iter().enumerate() {
        for &(i, v) in col {
            let (ni, nk) = (inv[i], inv[k]);
            let (lo, hi) = if ni <= nk { (ni, nk) } else { (nk, ni) };
            cols[hi].push((lo, v));
        }
    }
    for col in &mut cols {
        col.sort_unstable_by_key(|&(i, _)| i);
    }
    SymUpper { n, cols }
}

/// Nested dissection by recursive BFS bisection over the matrix graph.
/// Returns `perm[new] = old`. Separators are ordered last within each
/// recursion, which is what gives the elimination tree its bushiness.
pub(crate) fn nested_dissection(a: &SymUpper, levels: usize) -> Vec<usize> {
    let n = a.n;
    let mut adj: Vec<Vec<usize>> = vec![Vec::new(); n];
    for (k, col) in a.cols.iter().enumerate() {
        for &(i, _) in col {
            if i != k {
                adj[i].push(k);
                adj[k].push(i);
            }
        }
    }
    let mut order = Vec::with_capacity(n);
    let nodes: Vec<usize> = (0..n).collect();
    dissect(&adj, &nodes, levels, &mut order);
    order
}

fn dissect(adj: &[Vec<usize>], nodes: &[usize], levels: usize, order: &mut Vec<usize>) {
    if levels == 0 || nodes.len() < 32 {
        order.extend_from_slice(nodes);
        return;
    }
    // BFS levels from a pseudo-peripheral start; the median level set is the
    // separator.
    let inset: std::collections::HashSet<usize> = nodes.iter().copied().collect();
    let mut level = vec![usize::MAX; adj.len()];
    let mut start = nodes[0];
    for _ in 0..2 {
        for v in nodes.iter() {
            level[*v] = usize::MAX;
        }
        let mut queue = std::collections::VecDeque::from([start]);
        level[start] = 0;
        let mut last = start;
        while let Some(u) = queue.pop_front() {
            last = u;
            for &w in &adj[u] {
                if inset.contains(&w) && level[w] == usize::MAX {
                    level[w] = level[u] + 1;
                    queue.push_back(w);
                }
            }
        }
        start = last;
    }
    let mut reached: Vec<usize> = nodes
        .iter()
        .copied()
        .filter(|&v| level[v] != usize::MAX)
        .collect();
    let disconnected: Vec<usize> = nodes
        .iter()
        .copied()
        .filter(|&v| level[v] == usize::MAX)
        .collect();
    if reached.len() < nodes.len() / 2 {
        // Heavily disconnected — do not fight it, just emit.
        order.extend_from_slice(nodes);
        return;
    }
    reached.sort_unstable_by_key(|&v| (level[v], v));
    let median = level[reached[reached.len() / 2]];
    let left: Vec<usize> = reached
        .iter()
        .copied()
        .filter(|&v| level[v] < median)
        .collect();
    let right: Vec<usize> = reached
        .iter()
        .copied()
        .filter(|&v| level[v] > median)
        .collect();
    let mut sep: Vec<usize> = reached
        .iter()
        .copied()
        .filter(|&v| level[v] == median)
        .collect();
    sep.extend(disconnected);
    if left.is_empty() || right.is_empty() {
        order.extend_from_slice(nodes);
        return;
    }
    dissect(adj, &left, levels - 1, order);
    dissect(adj, &right, levels - 1, order);
    order.append(&mut sep);
}

/// Elimination tree of a `SymUpper` (Davis, LDL's symbolic half).
pub(crate) fn etree(a: &SymUpper) -> Vec<i64> {
    let n = a.n;
    let mut parent = vec![-1i64; n];
    let mut ancestor = vec![-1i64; n];
    for k in 0..n {
        for &(i, _) in &a.cols[k] {
            let mut i = i as i64;
            while i != -1 && (i as usize) < k {
                let next = ancestor[i as usize];
                ancestor[i as usize] = k as i64;
                if next == -1 {
                    parent[i as usize] = k as i64;
                }
                i = next;
            }
        }
    }
    parent
}

/// Full numeric up-looking LDLᵀ. Panics on a non-positive pivot — the heat
/// operator is SPD or the substrate is wrong.
pub(crate) fn ldl_factor(a: &SymUpper) -> LdlFactor {
    let parent = etree(a);
    let mut f = LdlFactor {
        n: a.n,
        parent,
        li: vec![Vec::new(); a.n],
        lx: vec![Vec::new(); a.n],
        d: vec![0.0; a.n],
        flops: 0,
    };
    let mut mask = vec![true; a.n];
    refactor_rows(a, &mut f, &mut mask);
    f
}

/// Recompute rows `k` with `mask[k]` true, ascending; entries for unmasked
/// rows are left in place. For a full factorization pass `mask = all true`;
/// for R-035b's update, the elimination-tree ancestor closure of the changed
/// columns. `mask` is consumed (cleared) as rows complete.
pub(crate) fn refactor_rows(a: &SymUpper, f: &mut LdlFactor, mask: &mut [bool]) {
    let n = a.n;
    let mut y = vec![0.0f64; n];
    let mut pattern = Vec::with_capacity(64);
    let mut flag = vec![usize::MAX; n];
    let mut flops = 0u64;
    for k in 0..n {
        if !mask[k] {
            continue;
        }
        // Scatter A(0..=k, k) and collect the row pattern by etree walk.
        flag[k] = k;
        let mut top = Vec::new();
        for &(i, v) in &a.cols[k] {
            if i == k {
                y[k] = v;
                continue;
            }
            y[i] += v;
            let mut ii = i;
            let mut chain = Vec::new();
            while flag[ii] != k {
                flag[ii] = k;
                chain.push(ii);
                let p = f.parent[ii];
                if p == -1 || p as usize >= k {
                    break;
                }
                ii = p as usize;
            }
            // chain is root-last; prepend reversed for ascending processing.
            chain.reverse();
            top.push(chain);
        }
        pattern.clear();
        for chain in top {
            for v in chain {
                pattern.push(v);
            }
        }
        pattern.sort_unstable();
        pattern.dedup();

        // Sparse triangular solve across the pattern, ascending.
        let mut dk = y[k];
        y[k] = 0.0;
        for &i in &pattern {
            let yi = y[i];
            y[i] = 0.0;
            // Never skip on yi == 0: a partial refactor may need to OVERWRITE
            // a stale non-zero L(k,i) with an exact zero.
            let lki = yi / f.d[i];
            // y -= L(:,i)·yi over rows in this row's future — only entries
            // with row < k matter for the solve; the row-k entry is L(k,i).
            for (&r, &v) in f.li[i].iter().zip(f.lx[i].iter()) {
                if r < k {
                    y[r] -= v * yi;
                } else if r == k {
                    // stale entry from a previous factorization state; the
                    // subtraction below re-derives it, so skip.
                } else {
                    break;
                }
            }
            flops += 2 * f.li[i].len() as u64;
            dk -= lki * yi;
            // Write L(k,i): replace if the slot exists, append otherwise.
            match f.li[i].binary_search(&k) {
                Ok(pos) => f.lx[i][pos] = lki,
                Err(pos) => {
                    f.li[i].insert(pos, k);
                    f.lx[i].insert(pos, lki);
                }
            }
        }
        assert!(dk > 0.0, "LDL pivot {k} not positive — operator not SPD");
        f.d[k] = dk;
        mask[k] = false;
    }
    f.flops += flops;
}

/// Ancestor closure of `seeds` in the elimination tree, as a mask.
pub(crate) fn ancestor_closure(parent: &[i64], seeds: &[usize]) -> Vec<bool> {
    let mut mask = vec![false; parent.len()];
    for &s in seeds {
        let mut i = s as i64;
        while i != -1 && !mask[i as usize] {
            mask[i as usize] = true;
            i = parent[i as usize];
        }
    }
    mask
}

/// Solve `A x = b` through the factor: `L D Lᵀ x = b`.
#[allow(
    clippy::needless_range_loop,
    reason = "triangular solves index the same vector they update; iterator forms obscure the algebra"
)]
pub(crate) fn solve(f: &LdlFactor, b: &[f64]) -> Vec<f64> {
    let n = f.n;
    let mut x = b.to_vec();
    for i in 0..n {
        let xi = x[i];
        if xi != 0.0 {
            for (&r, &v) in f.li[i].iter().zip(f.lx[i].iter()) {
                x[r] -= v * xi;
            }
        }
    }
    for i in 0..n {
        x[i] /= f.d[i];
    }
    for i in (0..n).rev() {
        let mut s = x[i];
        for (&r, &v) in f.li[i].iter().zip(f.lx[i].iter()) {
            s -= v * x[r];
        }
        x[i] = s;
    }
    x
}

/// ‖L·D·Lᵀ − A‖_F / ‖A‖_F, computed column-by-column against the upper form.
#[allow(
    clippy::needless_range_loop,
    reason = "the reconstruction scatters into a dense column by row index; iterator forms obscure the algebra"
)]
pub(crate) fn residual_fro(a: &SymUpper, f: &LdlFactor) -> f64 {
    // Reconstruct each upper-column of L·D·Lᵀ sparsely: (LDLᵀ)(i,k) for i ≤ k
    // = Σ_j L(i,j)·D(j)·L(k,j), with L unit diagonal.
    let n = a.n;
    let mut num = 0.0f64;
    let mut den = 0.0f64;
    let mut dense = vec![0.0f64; n];
    for k in 0..n {
        // Column k of L·D·Lᵀ restricted to rows ≤ k.
        // Contributions: j = k term: D[k]·L(i,k)? L(i,k)=0 for i<k (unit
        // lower); handle j < k via L(k,j) ≠ 0, plus j = k giving D[k] at i=k.
        dense[k] += f.d[k];
        for j in 0..=k {
            let lkj = if j == k {
                1.0
            } else {
                match f.li[j].binary_search(&k) {
                    Ok(pos) => f.lx[j][pos],
                    Err(_) => 0.0,
                }
            };
            if lkj == 0.0 {
                continue;
            }
            let w = f.d[j] * lkj;
            // rows i ≤ k with L(i,j) ≠ 0: i = j itself (unit), plus stored.
            if j < k {
                dense[j] += w;
            }
            for (&r, &v) in f.li[j].iter().zip(f.lx[j].iter()) {
                if r < k {
                    dense[r] += w * v;
                } else if r == k {
                    if j != k {
                        dense[k] += w * v;
                    }
                } else {
                    break;
                }
            }
        }
        for &(i, v) in &a.cols[k] {
            den += v * v;
            let diff = dense[i] - v;
            num += diff * diff;
            dense[i] = 0.0;
        }
        // Anything left in `dense` on this column is fill vs a structural
        // zero of A — count it as pure error.
        for j in 0..=k {
            if dense[j] != 0.0 {
                num += dense[j] * dense[j];
                dense[j] = 0.0;
            }
        }
    }
    (num / den).sqrt()
}
