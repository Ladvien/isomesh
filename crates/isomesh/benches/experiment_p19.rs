//! **P-19 — is the winding number's on-demand crossover set by row sharing?**
//!
//! Ticket: S-009. Pre-registered in the commit before this one.
//!
//! ```bash
//! cargo bench --bench experiment_p19
//! ```
//!
//! Writes `docs/experiments/p-19.csv`.
//!
//! # Why the obvious prediction is wrong
//!
//! ✗23 removed S-009's original justification — Manifold Dual Contouring does
//! *not* query where it needs to, it pre-samples all `N³` points at
//! `dual.rs:257`. What survives is the question for a genuinely sparse consumer,
//! and M-299 already measured the thing that decides it.
//!
//! `winding_numbers` casts **one ray per grid row**, shared across every sample
//! on it, so an `n³` grid costs `n²` casts. An on-demand field cannot share
//! anything: one ray per query. So the naive crossover — *"on-demand wins below
//! `N³` queries, since that is how many the batch path answers"* — charges the
//! batch path for work it does not do.
//!
//! # Why this measures the two terms rather than simulating on-demand
//!
//! The tempting harness is to call `winding_numbers` on a grid shaped so that
//! every row holds one sample, making rays and queries equal. **That shape does
//! not exist**: the function requires `size ≥ 2` on every axis, so the thinnest
//! grid still shares each ray between two samples. Simulating on-demand would
//! mean building the on-demand field first — and P-19's second falsifier can
//! *close* S-009, so building the thing to find out whether it is worth building
//! is backwards.
//!
//! Instead, decompose the batch cost into the two terms the hypothesis is about.
//! For a grid `[nx, ny, nz]` the work is one ray per `(y, z)` row plus one
//! boundary correction per sample, over a fixed per-call setup:
//!
//! ```text
//! T(nx, ny, nz)  =  s  +  (ny·nz)·r  +  (nx·ny·nz)·p
//! ```
//!
//! Three shapes give three equations in `s`, `r`, `p`:
//!
//! ```text
//! T_min  = T(2, 2, 2)  = s +   4·r +      8·p
//! T_thin = T(2, N, N)  = s +  N²·r +   2N²·p
//! T_full = T(N, N, N)  = s +  N²·r +    N³·p
//! ```
//!
//! `T_full − T_thin` isolates `p` exactly — the ray count is identical in both,
//! so every difference between them is per-sample work. `r` then follows from
//! `T_thin − T_min`, and `s` from either.
//!
//! An on-demand field pays `r + p` per query and the same `s` once, so the
//! crossover is where the marginal costs meet:
//!
//! ```text
//! Q*·(r + p)  =  N²·r + N³·p        ⟹        Q* = (N²·r + N³·p) / (r + p)
//! ```
//!
//! Read off the limits: **`p ≪ r` sends `Q*` to `N²`**, and `p ≫ r` sends it to
//! `N³`. That is the hypothesis, and it is now a ratio of two measured constants
//! rather than a curve fit.
//!
//! # The control is the point
//!
//! The boundary correction runs over **boundary edges**, so `p` scales with how
//! holed the mesh is while `r` does not. A nearly-closed mesh should put `Q*`
//! near `N²`; the same mesh with holes punched in it should push `Q*` up. **If
//! the hole-punched fixture does not move `Q*`, the mechanism named here is
//! wrong even if the nearly-closed number comes out right** — M-279's rule, that
//! a falsifier must separate the hypothesis from its rivals rather than merely be
//! capable of failing.

mod common;

use std::hint::black_box;
use std::time::Instant;

use isomesh::construct::winding::winding_numbers;
use isomesh::fields::Sphere;
use isomesh::marching_cubes::MarchingCubes;
use isomesh::{MeshBuffer, RuntimeShape3, Sdf, Shape3};

/// Grids the crossover is solved at. `17³` is the smallest where `N²` and `N³`
/// are far enough apart to tell the two predictions apart at all.
const RESOLUTIONS: [u32; 3] = [17, 33, 65];

/// Repeats of each timed shape. The full grid at `65³` is not cheap, and the
/// quantity being extracted is a *difference* between shapes, so noise in either
/// lands directly in `p`.
const REPEATS: u32 = 5;

/// A closed sphere mesh, and the same mesh with a fraction of its triangles
/// removed to open it up.
fn fixture(samples: u32, drop_every: Option<usize>) -> (Vec<[f64; 3]>, Vec<u32>) {
    let shape = RuntimeShape3::new([samples; 3]).expect("valid shape");
    let field = Sphere::<f64>::canonical();
    let origin = [-1.5_f64; 3];
    let h = 3.0 / f64::from(samples - 1);

    let mut values = Vec::with_capacity(shape.element_count());
    let size = shape.size();
    for z in 0..size[2] {
        for y in 0..size[1] {
            for x in 0..size[0] {
                values.push(field.sample([
                    origin[0] + h * f64::from(x),
                    origin[1] + h * f64::from(y),
                    origin[2] + h * f64::from(z),
                ]));
            }
        }
    }

    let mut out = MeshBuffer::<f64>::default();
    let sampled =
        isomesh::construct::SampledField::new(&values, &shape, origin, h).expect("wrap grid");
    MarchingCubes::<f64>::new()
        .extract(&sampled, &shape, origin, h, &mut out)
        .expect("extraction");

    let positions = out.positions.clone();
    let indices = match drop_every {
        None => out.indices.clone(),
        // Punch holes by dropping whole triangles. Every dropped triangle opens
        // three edges, so this raises the boundary-edge count sharply without
        // touching the triangle count much -- which is what separates `p` from
        // `r`, since the ray cost scales with triangles and the correction with
        // boundary edges.
        Some(step) => out
            .indices
            .as_chunks::<3>()
            .0
            .iter()
            .enumerate()
            .filter(|(i, _)| !i.is_multiple_of(step))
            .flat_map(|(_, t)| t.iter().copied())
            .collect(),
    };
    (positions, indices)
}

/// Directed edges with no opposite partner — what the correction term runs over.
fn boundary_edge_count(indices: &[u32]) -> usize {
    use std::collections::BTreeMap;
    let mut net: BTreeMap<(u32, u32), i32> = BTreeMap::new();
    for tri in indices.as_chunks::<3>().0 {
        for (u, v) in [(tri[0], tri[1]), (tri[1], tri[2]), (tri[2], tri[0])] {
            if u == v {
                continue;
            }
            let (key, delta) = if u < v { ((u, v), 1) } else { ((v, u), -1) };
            *net.entry(key).or_insert(0) += delta;
        }
    }
    net.values().filter(|&&n| n != 0).count()
}

/// Best-of-`REPEATS` nanoseconds for one grid shape. Minimum rather than mean:
/// the quantity wanted is the cost of the work, and every source of noise on a
/// governed CPU adds to it.
fn time_shape(positions: &[[f64; 3]], indices: &[u32], size: [u32; 3]) -> f64 {
    let shape = RuntimeShape3::new(size).expect("valid shape");
    let origin = [-1.5_f64; 3];
    let h = 3.0 / f64::from(size[0].max(2) - 1);

    let mut best = f64::INFINITY;
    for _ in 0..REPEATS {
        let start = Instant::now();
        let w = winding_numbers(positions, indices, &shape, origin, h).expect("winding");
        black_box(&w);
        let ns = start.elapsed().as_secs_f64() * 1e9;
        if ns < best {
            best = ns;
        }
    }
    best
}

struct Solved {
    setup_ns: f64,
    per_ray_ns: f64,
    per_point_ns: f64,
    crossover: f64,
}

/// Solve `s`, `r`, `p` from the three shapes, then the crossover.
fn solve(n: u32, t_min: f64, t_thin: f64, t_full: f64) -> Solved {
    let n2 = f64::from(n) * f64::from(n);
    let n3 = n2 * f64::from(n);

    // T_full - T_thin = (N^3 - 2N^2) p. Same ray count in both, so the whole
    // difference is per-sample work.
    let per_point_ns = (t_full - t_thin) / (n3 - 2.0 * n2);
    // (T_thin - 2N^2 p) - (T_min - 8 p) = (N^2 - 4) r.
    let per_ray_ns =
        ((t_thin - 2.0 * n2 * per_point_ns) - (t_min - 8.0 * per_point_ns)) / (n2 - 4.0);
    let setup_ns = t_min - 4.0 * per_ray_ns - 8.0 * per_point_ns;

    let marginal = per_ray_ns + per_point_ns;
    let crossover = if marginal > 0.0 {
        (n2 * per_ray_ns + n3 * per_point_ns) / marginal
    } else {
        f64::NAN
    };

    Solved {
        setup_ns,
        per_ray_ns,
        per_point_ns,
        crossover,
    }
}

fn main() {
    if !std::env::args().any(|arg| arg == "--bench") {
        return;
    }

    let prereg = isomesh::experiment!("P-19");
    common::experiment::run(prereg, |run| {
        println!(
            "{:>5} {:>8} {:>9} {:>12} {:>12} {:>12} {:>12} {:>9} {:>9}",
            "n", "tris", "bnd edges", "batch ns", "ray ns", "point ns", "Q*", "Q*/N^2", "Q*/N^3"
        );
        for samples in RESOLUTIONS {
            for drop_every in [None, Some(7usize)] {
                let (positions, indices) = fixture(samples, drop_every);
                let triangles = indices.len() / 3;
                let boundary_edges = boundary_edge_count(&indices);

                let t_min = time_shape(&positions, &indices, [2, 2, 2]);
                let t_thin = time_shape(&positions, &indices, [2, samples, samples]);
                let t_full = time_shape(&positions, &indices, [samples; 3]);

                let s = solve(samples, t_min, t_thin, t_full);
                let n2 = f64::from(samples) * f64::from(samples);
                let n3 = n2 * f64::from(samples);

                println!(
                    "{:>5} {:>8} {:>9} {:>12.0} {:>12.1} {:>12.4} {:>12.0} {:>9.2} {:>9.4}",
                    samples,
                    triangles,
                    boundary_edges,
                    t_full,
                    s.per_ray_ns,
                    s.per_point_ns,
                    s.crossover,
                    s.crossover / n2,
                    s.crossover / n3
                );

                run.record(&[
                    ("samples_per_axis", samples.to_string()),
                    ("triangles", triangles.to_string()),
                    ("boundary_edges", boundary_edges.to_string()),
                    ("batch_total_ns", format!("{:.0}", t_full)),
                    (
                        "on_demand_ns_per_query",
                        format!("{:.4}", s.per_ray_ns + s.per_point_ns),
                    ),
                    ("crossover_queries", format!("{:.0}", s.crossover)),
                    (
                        "crossover_over_n_squared",
                        format!("{:.4}", s.crossover / n2),
                    ),
                    ("crossover_over_n_cubed", format!("{:.6}", s.crossover / n3)),
                ]);

                let _ = s.setup_ns;
            }
        }
    });
}
