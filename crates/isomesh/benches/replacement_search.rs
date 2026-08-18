//! **When a fill severs a cave, how big is the smaller side — and how many of
//! the deleted edges were in the spanning forest at all?**
//!
//! Ticket: R-022b. **Exploratory. Nothing is registered against this run.** It
//! checks the arithmetic of two clauses before either is written into a `P-`,
//! which is the house style this phase arrived at the hard way: it failed twice
//! by registering first (P-23 clause 3, P-24) and got it right five times by
//! measuring first.
//!
//! ```bash
//! cargo bench --bench replacement_search
//! ```
//!
//! Writes `docs/measurements/replacement_search.csv`.
//!
//! # The two clauses
//!
//! M-319 established that R-022b is real: one fill in six changes the air
//! component count, so detect-and-recompute would rebuild at `O(n^3)` every
//! sixth edit. What it did **not** establish is which deletion structure the
//! lattice actually needs, and two numbers from the literature decide that.
//!
//! **Clause 1 — most deletions should be free.** Deleting a *non-tree* edge
//! cannot change connectivity, and HDT distinguishes the two cases in `O(1)`
//! from a hash table; `10.48550/arXiv.2411.11781` names this as the reason HDT
//! stays practical on dense graphs, *"where most of the deletions target
//! non-tree edges"*. A 6-connected lattice has `|E| ~ 3|V|` against a spanning
//! forest of `|V| - c`, so the predicted free fraction is about **two thirds**.
//! That is arithmetic, and arithmetic is checkable.
//!
//! **Clause 2 — and this one decides the whole ticket.** D-Tree
//! (`10.48550/arXiv.2509.14433`, describing [13]) finds a replacement edge by
//! running BFS *on the smaller component of the split*, with **no theoretical
//! guarantee** and good measured behaviour. Its cost is therefore the size of
//! the smaller side. So: when a fill severs a cave passage, does it shed a
//! 40-voxel pocket or half the cave?
//!
//! If the smaller side is reliably tiny, BFS-on-the-smaller-side is cheap and
//! the levelled HDT machinery — `O(log n)` levels, edges pushed down a level to
//! amortise a failed search — is being paid for a case the lattice does not
//! present. That would **shrink** R-022b rather than confirm its shape. If the
//! smaller side is a constant fraction of the domain, it would not.
//!
//! # Method, and why nothing here is maintained incrementally
//!
//! Same discipline as M-319's harness: the spanning forest is rebuilt by BFS
//! before each fill and the component labelling is rebuilt after it. Both are
//! `O(n^3)` and this is a measurement rather than a hot path — and a structure
//! that maintained either could not be used as evidence about whether
//! maintaining it is worthwhile.
//!
//! Removal is the only operation, so a component can **split or vanish but
//! never merge**. That is what lets each post-fill component be attributed to
//! exactly one pre-fill component by looking at any single member.

// Exact equality throughout: the question is whether a value changed at all.
#![allow(
    clippy::float_cmp,
    reason = "the question is whether a value changed at all"
)]

mod common;

use std::fmt::Write as _;

use isomesh::fields::{ReferenceField, noise_cavity};
use isomesh::{RuntimeShape3, Sdf, Shape3};

/// Samples per axis. Matches M-319 so the two tables can be read together.
const RESOLUTIONS: [u32; 3] = [33, 49, 65];

/// Brush radius in samples, fixed across resolutions. Matches M-319.
const BRUSH_RADIUS: f64 = 4.0;

/// Fills applied per resolution. Matches M-319.
const FILLS: usize = 200;

/// Sentinel for "no parent" / "not air".
const NONE: u32 = u32::MAX;

/// The generator, so the fill sequence is the same on every machine — and the
/// same sequence M-319 used, so the split counts must agree.
struct Lcg(u64);

impl Lcg {
    fn next_u64(&mut self) -> u64 {
        self.0 = self
            .0
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        self.0
    }

    fn below(&mut self, n: u32) -> u32 {
        (self.next_u64() >> 33) as u32 % n
    }
}

/// The six axis-aligned neighbours of `i`, as flat indices, skipping any that
/// leave the lattice.
fn neighbours(i: usize, n: usize, out: &mut [usize; 6]) -> usize {
    let x = i % n;
    let y = (i / n) % n;
    let z = i / (n * n);
    let mut count = 0;
    let mut push = |v: usize| {
        if let Some(slot) = out.get_mut(count) {
            *slot = v;
            count += 1;
        }
    };
    if x > 0 {
        push(i - 1);
    }
    if x + 1 < n {
        push(i + 1);
    }
    if y > 0 {
        push(i - n);
    }
    if y + 1 < n {
        push(i + n);
    }
    if z > 0 {
        push(i - n * n);
    }
    if z + 1 < n {
        push(i + n * n);
    }
    count
}

/// A BFS spanning forest of the air region, plus the component labelling it
/// induces. `parent[i]` is `NONE` for a root or a solid sample.
struct Forest {
    parent: Vec<u32>,
    label: Vec<u32>,
    sizes: Vec<u64>,
}

/// Air is the sublevel set, matching `connectivity::Air`.
fn is_air(v: f64) -> bool {
    v >= 0.0
}

fn forest(values: &[f64], n: usize) -> Forest {
    let count = values.len();
    let mut parent = vec![NONE; count];
    let mut label = vec![NONE; count];
    let mut sizes: Vec<u64> = Vec::new();
    let mut queue: Vec<usize> = Vec::new();
    let mut nb = [0usize; 6];

    for start in 0..count {
        let Some(&v) = values.get(start) else {
            continue;
        };
        if !is_air(v) {
            continue;
        }
        if label.get(start).is_some_and(|l| *l != NONE) {
            continue;
        }
        let this = sizes.len() as u32;
        sizes.push(0);
        if let Some(slot) = label.get_mut(start) {
            *slot = this;
        }
        queue.clear();
        queue.push(start);
        let mut seen = 0u64;
        let mut head = 0;
        while let Some(&i) = queue.get(head) {
            head += 1;
            seen += 1;
            let used = neighbours(i, n, &mut nb);
            for &j in nb.iter().take(used) {
                if !values.get(j).copied().is_some_and(is_air) {
                    continue;
                }
                if label.get(j).is_some_and(|l| *l != NONE) {
                    continue;
                }
                if let Some(slot) = label.get_mut(j) {
                    *slot = this;
                }
                if let Some(slot) = parent.get_mut(j) {
                    *slot = i as u32;
                }
                queue.push(j);
            }
        }
        if let Some(slot) = sizes.get_mut(this as usize) {
            *slot = seen;
        }
    }

    Forest {
        parent,
        label,
        sizes,
    }
}

fn main() {
    if !std::env::args().any(|arg| arg == "--bench") {
        return;
    }

    println!(
        "{:>5} {:>9} {:>8} {:>10} {:>9} {:>7} {:>9} {:>9} {:>10}",
        "n", "air", "deleted", "tree-edge", "free %", "splits", "min side", "med side", "max side"
    );
    let mut csv = String::from(
        "samples_per_axis,air_samples,edges_deleted,tree_edges_deleted,free_pct,\
         splits,smaller_side_min,smaller_side_median,smaller_side_max,\
         smaller_side_median_pct_of_air\n",
    );

    for n in RESOLUTIONS {
        let Ok(shape) = RuntimeShape3::new([n; 3]) else {
            continue;
        };
        let nn = n as usize;
        let field = noise_cavity::<f64>();
        let (lo, hi) = field.domain();
        let h = (hi[0] - lo[0]) / f64::from(n - 1);

        let mut values = Vec::with_capacity(shape.element_count());
        for z in 0..n {
            for y in 0..n {
                for x in 0..n {
                    values.push(field.sample([
                        lo[0] + h * f64::from(x),
                        lo[1] + h * f64::from(y),
                        lo[2] + h * f64::from(z),
                    ]));
                }
            }
        }

        let air_samples = values.iter().filter(|v| is_air(**v)).count() as u64;
        let mut rng = Lcg(0x0000_A022_B000_0001);
        let (mut deleted, mut tree_deleted) = (0u64, 0u64);
        let mut smaller_sides: Vec<u64> = Vec::new();
        let mut filled_flag = vec![false; values.len()];
        let mut nb = [0usize; 6];

        for _ in 0..FILLS {
            // The forest as it stands *before* this fill: the edges the fill is
            // about to delete are tree edges or not with respect to this one.
            let before = forest(&values, nn);

            let centre = [rng.below(n), rng.below(n), rng.below(n)];
            let r = BRUSH_RADIUS.ceil() as i64;
            filled_flag.iter_mut().for_each(|f| *f = false);
            let mut touched: Vec<usize> = Vec::new();

            for dz in -r..=r {
                for dy in -r..=r {
                    for dx in -r..=r {
                        if (dx * dx + dy * dy + dz * dz) as f64 > BRUSH_RADIUS * BRUSH_RADIUS {
                            continue;
                        }
                        let p = [
                            i64::from(centre[0]) + dx,
                            i64::from(centre[1]) + dy,
                            i64::from(centre[2]) + dz,
                        ];
                        if p.iter().any(|c| *c < 0 || *c >= i64::from(n)) {
                            continue;
                        }
                        let i = (p[2] as usize * nn + p[1] as usize) * nn + p[0] as usize;
                        if values.get(i).copied().is_some_and(is_air) {
                            if let Some(slot) = filled_flag.get_mut(i) {
                                *slot = true;
                            }
                            touched.push(i);
                        }
                    }
                }
            }

            // Every air-air edge this fill removes, counted once. An edge with
            // both ends filled would otherwise be counted twice, so the higher
            // index yields to the lower.
            for &i in &touched {
                let used = neighbours(i, nn, &mut nb);
                for &j in nb.iter().take(used) {
                    if !values.get(j).copied().is_some_and(is_air) {
                        continue;
                    }
                    if filled_flag.get(j).copied().unwrap_or(false) && j < i {
                        continue;
                    }
                    deleted += 1;
                    let tree = before.parent.get(i).copied() == Some(j as u32)
                        || before.parent.get(j).copied() == Some(i as u32);
                    if tree {
                        tree_deleted += 1;
                    }
                }
            }

            for &i in &touched {
                if let Some(slot) = values.get_mut(i) {
                    *slot = -1.0;
                }
            }

            // Removal never merges, so every post-fill component belongs to
            // exactly one pre-fill component — read off any single member.
            let after = forest(&values, nn);
            let mut pieces: Vec<Vec<u64>> = vec![Vec::new(); before.sizes.len()];
            let mut first_member = vec![NONE; after.sizes.len()];
            for (i, &l) in after.label.iter().enumerate() {
                if l == NONE {
                    continue;
                }
                if let Some(slot) = first_member.get_mut(l as usize)
                    && *slot == NONE
                {
                    *slot = i as u32;
                }
            }
            for (new_label, &member) in first_member.iter().enumerate() {
                if member == NONE {
                    continue;
                }
                let Some(&old) = before.label.get(member as usize) else {
                    continue;
                };
                if old == NONE {
                    continue;
                }
                let Some(&size) = after.sizes.get(new_label) else {
                    continue;
                };
                if let Some(slot) = pieces.get_mut(old as usize) {
                    slot.push(size);
                }
            }
            for group in &pieces {
                if group.len() < 2 {
                    continue;
                }
                // The split shed one or more pieces; D-Tree's cost is the
                // smaller side, so that is what gets recorded.
                if let Some(&min) = group.iter().min() {
                    smaller_sides.push(min);
                }
            }
        }

        smaller_sides.sort_unstable();
        let splits = smaller_sides.len() as u64;
        let pick = |q: usize| smaller_sides.get(q).copied().unwrap_or(0);
        let (min_side, med_side, max_side) = if smaller_sides.is_empty() {
            (0, 0, 0)
        } else {
            (
                pick(0),
                pick(smaller_sides.len() / 2),
                pick(smaller_sides.len() - 1),
            )
        };
        let free_pct = if deleted == 0 {
            0.0
        } else {
            (deleted - tree_deleted) as f64 / deleted as f64 * 100.0
        };
        let med_pct = if air_samples == 0 {
            0.0
        } else {
            med_side as f64 / air_samples as f64 * 100.0
        };

        println!(
            "{n:>5} {air_samples:>9} {deleted:>8} {tree_deleted:>10} {free_pct:>8.1}% \
             {splits:>7} {min_side:>9} {med_side:>9} {max_side:>10}"
        );
        let _ = writeln!(
            csv,
            "{n},{air_samples},{deleted},{tree_deleted},{free_pct:.2},{splits},\
             {min_side},{med_side},{max_side},{med_pct:.4}"
        );
    }

    let path = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../docs/measurements/replacement_search.csv");
    match std::fs::write(&path, &csv) {
        Ok(()) => println!("\nwrote {}", path.display()),
        Err(e) => println!("\n::error:: {}: {e}", path.display()),
    }
}
