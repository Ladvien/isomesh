//! **P-51 — the tangent-sphere constraint has two halves, and both are counted.**
//!
//! Ticket: R-046. Pre-registered in `src/experiment.rs` before this harness
//! existed; the registration is not amended by anything below.
//!
//! ```bash
//! cargo bench --bench experiment_p51
//! ```
//!
//! Writes `docs/experiments/p-51.csv`.
//!
//! # What is being measured
//!
//! Sellan, Batty & Stein (`10.1145/3610548.3618196`) state what a signed
//! distance sample `p` with value `d(p)` actually constrains: the surface must
//! **exclude** the open ball of radius `|d(p)|` around `p`, *and* it must
//! **touch** that ball's boundary sphere at least once. Every extractor in this
//! crate reads `d` as a number to interpolate on an edge and discards both
//! halves. This harness counts how far each half is from being satisfied, as
//! integers over the extractor's own output.
//!
//! **PIERCING — the exclusion half.** For a vertex `v`,
//! `violation(v) = max over p of (|d(p)| − ‖v − p‖)`, normalised by the cell
//! size `h`, and a vertex counts as *pierced* when that exceeds `0.05` cells.
//! Per the registration the samples `p` are those "within one cell of `v`",
//! which here is the 4³ = 64 sample index block `{c−1, c, c+1, c+2}` per axis
//! around `v`'s own cell `c` — the eight corners of that cell plus one
//! neighbouring layer on every side. That is a *superset* of the strict
//! Chebyshev-`h` ball, and a superset can only find piercing, never invent it:
//! a positive `|d(p)| − ‖v − p‖` at *any* range means a mesh vertex sits
//! strictly inside a sample's ball, which is a violation whatever the range.
//!
//! The neighbourhood is a tractability device in the registration, so this file
//! also reports `pierced_exhaustive` — the same quantity with the maximum taken
//! over **all 274,625 samples**. It gates nothing; it exists so that "the
//! neighbourhood hid the piercing" is a checkable claim rather than a worry.
//!
//! **TOUCHING — the tangency half, which no extractor here attempts.** For a
//! sample `p`, `touch(p) = min over mesh vertices v of | ‖v − p‖ − |d(p)| |`,
//! normalised by `h`, and `p` counts as *untouched* above `0.05` cells.
//!
//! # There is no neighbourhood cutoff on the touching half, and that is the point
//!
//! The registration allows restricting the touching search for tractability, and
//! `samples_probed_per_vertex` exists to record the radius. **No restriction was
//! needed and none is applied: the minimum is taken over every vertex of the
//! mesh, for every one of the 274,625 samples.** So the untouched count cannot
//! be an artefact of a cutoff, because there is no cutoff — `touch_search` reads
//! `exhaustive` on every row and `vertices_probed_per_sample` is the full vertex
//! count.
//!
//! That is affordable because the quantity being minimised is monotone in the
//! *squared* distance around the target radius. A vertex can only beat the
//! current best `t` if its distance lands in the open interval
//! `(r − t, r + t)`, i.e. its squared distance lands in `((r−t)², (r+t)²)`, so
//! the inner loop is three subtractions, three multiplies, two adds and two
//! compares, and the square root is taken only on the rare candidate that
//! actually narrows the window. `65³ samples × ~10⁴ vertices` is a few seconds
//! per row at that price, which is cheaper than any correct spatial index and
//! has no cutoff to justify.
//!
//! `samples_probed_per_vertex` is therefore reported as what its name says —
//! **64**, the piercing half's per-vertex sample block — and the touching half's
//! (absent) radius is stated in `touch_search`. One column cannot honestly carry
//! both: the two halves quantify in opposite orders, `max over p for each v`
//! against `min over v for each p`.
//!
//! # The five fields, and one thing the registration got wrong
//!
//! The registration names "the five fields declaring `FieldBound::Exact` --
//! sphere, torus, box_exact, thin_plate, csg_difference". **Four fields declare
//! `FieldBound::Exact`.** `csg_difference` declares
//! `FieldBound::Underestimate { q: 0.5 }`, and `fields/mod.rs` is explicit about
//! why: it is `max(box, −sphere)`, and `max` of two exact distances is not an
//! exact distance. The field list is measured as registered and the actual bound
//! is recorded per row in `field_bound`, so the artefact says which rows rest on
//! an exact distance and which do not.
//!
//! It matters in a stateable direction. `q·d ≤ |f| ≤ d` means `csg_difference`
//! never *overstates* distance, so its balls are too small: piercing is
//! **under**-reported there (a smaller ball is harder to pierce) and tangency is
//! **over**-reported as violated (a sphere strictly inside the true tangency
//! sphere is one the surface cannot reach). Its row is honest about both because
//! `field_bound` names the reason.
//!
//! # Every row prints its own worst case, with the arithmetic
//!
//! A count is not auditable on its own, and the two headline numbers here are
//! both extremes. So each row prints the geometry behind its worst piercing
//! (when there is one) and behind its worst untouched sphere: the sample, the
//! sphere's radius, the vertex, and the distance. Both witnesses are one
//! subtraction away from the reported figure, and the first one this harness
//! produced checks out exactly by hand — see the control's verdict in the
//! ticket's findings row.
//!
//! # Determinism
//!
//! No map iteration, no PRNG, no threads, no clock in any recorded quantity.
//! Every number is a count, a ratio of counts, or a maximum of a geometric
//! quantity, and the two reductions (`max` over samples, `min` over vertices)
//! are order-independent up to ties that carry the same value.
//! `nonfinite_field_samples` counts any `|d(p)|` that is not finite, so a NaN
//! is reported rather than silently losing every comparison it touches.

mod common;

use isomesh::extractor::Extractor;
use isomesh::fields::{
    BoxExact, FieldBound, ReferenceField, Sphere, ThinPlate, Torus, csg_difference,
};
use isomesh::{MeshBuffer, Sdf};

/// Samples per axis. Registered.
const SAMPLES_PER_AXIS: u32 = 65;

/// The gate, in cells. Registered, and derived rather than chosen: M-12
/// measured `h²` convergence, and `sphere` at 65³ has mean error `6.5e-4`
/// against `h = 0.0625`, i.e. 1.0% of a cell. This sits five times above that
/// honest discretisation floor and far below M-27's 0.35 cells.
const THRESHOLD_CELLS: f64 = 0.05;

/// Sample index layers on each side of the cell that produced a vertex.
///
/// One layer out from the cell's own two planes gives four sample planes per
/// axis, `{c−1, c, c+1, c+2}`.
const PROBE_LAYERS: i64 = 1;

/// `4³` — the samples the piercing half probes per vertex.
const SAMPLES_PROBED_PER_VERTEX: usize = 64;

/// The three extractors P-51 covers, named so the registry selection cannot
/// drift into a silent fourth.
const COVERED: [&str; 3] = ["marching_cubes", "surface_nets", "dual_contouring"];

/// One `(field, extractor)` measurement.
struct Row {
    field: &'static str,
    bound: &'static str,
    extractor: &'static str,
    cell_size: f64,
    samples: usize,
    nonfinite: u64,
    vertices: usize,
    triangles: usize,
    /// Piercing, over the registered one-cell neighbourhood.
    pierced: u64,
    worst_piercing: f64,
    /// Piercing, over every sample in the grid. Gates nothing.
    pierced_exhaustive: u64,
    worst_piercing_exhaustive: f64,
    untouched: u64,
    worst_untouched: f64,
}

impl Row {
    /// Pierced vertices per 1,000 vertices.
    fn pierced_per_1k(&self) -> f64 {
        if self.vertices == 0 {
            return 0.0;
        }
        1000.0 * self.pierced as f64 / self.vertices as f64
    }

    /// Pierced vertices per 1,000, over the whole grid rather than one cell.
    fn pierced_exhaustive_per_1k(&self) -> f64 {
        if self.vertices == 0 {
            return 0.0;
        }
        1000.0 * self.pierced_exhaustive as f64 / self.vertices as f64
    }

    /// Untouched samples per 1,000 samples.
    fn untouched_per_1k(&self) -> f64 {
        if self.samples == 0 {
            return 0.0;
        }
        1000.0 * self.untouched as f64 / self.samples as f64
    }
}

/// A ratio of two rates, spelled out when the denominator is zero.
///
/// `0 / 0` is not "one" and `x / 0` is not "large": both are the absence of a
/// ratio, and a CSV that prints `NaN` or `inf` there says so, where a silently
/// substituted number would not.
fn ratio(numerator: f64, denominator: f64) -> String {
    let r = numerator / denominator;
    if r.is_finite() {
        format!("{r:.4}")
    } else if r.is_nan() {
        String::from("undefined")
    } else {
        String::from("inf")
    }
}

/// What a field's declared bound is, as one word for the CSV.
fn bound_name(bound: FieldBound) -> &'static str {
    match bound {
        FieldBound::Exact => "exact",
        FieldBound::Lipschitz { .. } => "lipschitz",
        FieldBound::Underestimate { .. } => "underestimate",
        FieldBound::Unbounded => "unbounded",
    }
}

/// The piercing half's result, with the geometry that produced the worst case.
///
/// The witness is here because the control is the headline: a claim that a
/// vertex sits inside a sample's ball should name the vertex, the sample, the
/// ball's radius and the distance, so a reader can check the subtraction by
/// hand rather than trusting a count.
struct Piercing {
    pierced: u64,
    worst_cells: f64,
    vertex: [f64; 3],
    sample: [f64; 3],
    radius: f64,
    distance: f64,
}

/// The exhaustive scan's result: the touching half, plus the piercing half
/// recomputed with no neighbourhood at all.
///
/// The witness names the sample whose sphere is missed by the widest margin,
/// the sphere's radius, and the vertex that comes closest to it. That is the
/// number this experiment exists to produce, so it should be checkable by
/// subtracting two printed floats rather than by trusting a count.
struct Scan {
    untouched: u64,
    worst_untouched_cells: f64,
    pierced: u64,
    worst_piercing_cells: f64,
    sample: [f64; 3],
    radius: f64,
    nearest_shell: [f64; 3],
    shell_distance: f64,
}

/// The sampled grid: `|d(p)|` at every sample, in `k`-major order.
struct Samples {
    abs_d: Vec<f64>,
    lo: [f64; 3],
    h: f64,
    n: usize,
    nonfinite: u64,
}

impl Samples {
    /// Sample `field` at `SAMPLES_PER_AXIS` per axis over its own domain.
    fn of<F>(field: &F) -> (Self, isomesh::RuntimeShape3)
    where
        F: ReferenceField + Sdf<Scalar = f64>,
    {
        let (shape, lo, h) = common::grid(field, SAMPLES_PER_AXIS);
        let n = SAMPLES_PER_AXIS as usize;
        let mut abs_d = vec![0.0_f64; n * n * n];
        let mut nonfinite = 0;
        for k in 0..n {
            let z = lo[2] + k as f64 * h;
            for j in 0..n {
                let y = lo[1] + j as f64 * h;
                for i in 0..n {
                    let x = lo[0] + i as f64 * h;
                    let a = field.sample([x, y, z]).abs();
                    if !a.is_finite() {
                        nonfinite += 1;
                    }
                    abs_d[(k * n + j) * n + i] = a;
                }
            }
        }
        (
            Self {
                abs_d,
                lo,
                h,
                n,
                nonfinite,
            },
            shape,
        )
    }

    /// Position of the sample at index `(i, j, k)`.
    fn at(&self, i: usize, j: usize, k: usize) -> [f64; 3] {
        [
            self.lo[0] + i as f64 * self.h,
            self.lo[1] + j as f64 * self.h,
            self.lo[2] + k as f64 * self.h,
        ]
    }

    /// PIERCING, over the registered one-cell neighbourhood.
    ///
    /// The worst violation is the raw maximum and may be **negative**: for an
    /// exact field every point of the true surface is at least `|d(p)|` from
    /// `p`, so a mesh that tracks the surface has a negative margin everywhere,
    /// and how close that margin comes to zero is the interesting part.
    fn piercing(&self, verts: &[[f64; 3]]) -> Piercing {
        let inv_h = self.h.recip();
        let last = self.n as i64 - 1;
        let mut found = Piercing {
            pierced: 0,
            worst_cells: if verts.is_empty() {
                0.0
            } else {
                f64::NEG_INFINITY
            },
            vertex: [0.0; 3],
            sample: [0.0; 3],
            radius: 0.0,
            distance: 0.0,
        };
        for v in verts {
            let mut cell = [0_i64; 3];
            for axis in 0..3 {
                let c = ((v[axis] - self.lo[axis]) * inv_h).floor() as i64;
                cell[axis] = c.clamp(0, last - 1);
            }
            let mut violation = f64::NEG_INFINITY;
            let mut at = [0.0_f64; 3];
            let mut radius = 0.0;
            let mut distance = 0.0;
            for k in cell[2] - PROBE_LAYERS..=cell[2] + 1 + PROBE_LAYERS {
                if k < 0 || k > last {
                    continue;
                }
                for j in cell[1] - PROBE_LAYERS..=cell[1] + 1 + PROBE_LAYERS {
                    if j < 0 || j > last {
                        continue;
                    }
                    for i in cell[0] - PROBE_LAYERS..=cell[0] + 1 + PROBE_LAYERS {
                        if i < 0 || i > last {
                            continue;
                        }
                        let (iu, ju, ku) = (i as usize, j as usize, k as usize);
                        let p = self.at(iu, ju, ku);
                        let r = self.abs_d[(ku * self.n + ju) * self.n + iu];
                        let dx = v[0] - p[0];
                        let dy = v[1] - p[1];
                        let dz = v[2] - p[2];
                        let d = (dx * dx + dy * dy + dz * dz).sqrt();
                        if r - d > violation {
                            violation = r - d;
                            at = p;
                            radius = r;
                            distance = d;
                        }
                    }
                }
            }
            let cells = violation * inv_h;
            if cells > THRESHOLD_CELLS {
                found.pierced += 1;
            }
            if cells > found.worst_cells {
                found.worst_cells = cells;
                found.vertex = *v;
                found.sample = at;
                found.radius = radius;
                found.distance = distance;
            }
        }
        found
    }

    /// TOUCHING over every vertex for every sample, and PIERCING over every
    /// sample for every vertex, in one pass.
    ///
    /// Both reductions are **exact over the whole grid** — nothing is cut off —
    /// and both are pruned by the same trick from opposite sides.
    ///
    /// *Touching* keeps a window: a vertex beats the running best `t` only if
    /// `‖v − p‖ ∈ (r − t, r + t)`, so squared distances outside
    /// `((r−t)², (r+t)²)` are rejected without a square root.
    ///
    /// *Piercing* keeps a per-vertex running minimum `m` of `‖v − p‖ − |d(p)|`,
    /// whose negation is the violation. A sample improves it only if
    /// `‖v − p‖ < r + m`, so squared distances at or above `(r + m)²` are
    /// rejected without a square root, and `r + m ≤ 0` rejects the sample
    /// outright. `m` starts at `+∞`, so the first sample always evaluates and
    /// the maximum that comes out is the true one rather than a value clipped at
    /// the gate.
    fn scan(&self, verts: &[[f64; 3]]) -> Scan {
        let inv_h = self.h.recip();
        let mut margin = vec![f64::INFINITY; verts.len()];
        let mut out = Scan {
            untouched: 0,
            worst_untouched_cells: 0.0,
            pierced: 0,
            worst_piercing_cells: 0.0,
            sample: [0.0; 3],
            radius: 0.0,
            nearest_shell: [0.0; 3],
            shell_distance: 0.0,
        };
        for k in 0..self.n {
            for j in 0..self.n {
                for i in 0..self.n {
                    let p = self.at(i, j, k);
                    let r = self.abs_d[(k * self.n + j) * self.n + i];
                    // Touching: the window that can still beat `best`.
                    let mut best = f64::INFINITY;
                    let mut best_at = [0.0_f64; 3];
                    let mut best_distance = 0.0;
                    let mut window_lo = f64::NEG_INFINITY;
                    let mut window_hi = f64::INFINITY;
                    for (vi, v) in verts.iter().enumerate() {
                        let dx = v[0] - p[0];
                        let dy = v[1] - p[1];
                        let dz = v[2] - p[2];
                        let d2 = dx * dx + dy * dy + dz * dz;
                        if d2 > window_lo && d2 < window_hi {
                            let d = d2.sqrt();
                            let t = (d - r).abs();
                            if t < best {
                                best = t;
                                best_at = *v;
                                best_distance = d;
                                let inner = r - t;
                                window_lo = if inner > 0.0 {
                                    inner * inner
                                } else {
                                    f64::NEG_INFINITY
                                };
                                window_hi = (r + t) * (r + t);
                            }
                        }
                        let reach = r + margin[vi];
                        if reach > 0.0 && d2 < reach * reach {
                            let m = d2.sqrt() - r;
                            if m < margin[vi] {
                                margin[vi] = m;
                            }
                        }
                    }
                    let cells = best * inv_h;
                    if cells > THRESHOLD_CELLS {
                        out.untouched += 1;
                    }
                    if cells > out.worst_untouched_cells {
                        out.worst_untouched_cells = cells;
                        out.sample = p;
                        out.radius = r;
                        out.nearest_shell = best_at;
                        out.shell_distance = best_distance;
                    }
                }
            }
        }
        let mut worst_pierce = f64::NEG_INFINITY;
        for m in &margin {
            // `0.0 - m` rather than `-m`: an exactly tangent vertex has
            // `m == 0`, and unary negation would report it as `-0.0`.
            let cells = (0.0 - m) * inv_h;
            if cells > THRESHOLD_CELLS {
                out.pierced += 1;
            }
            if cells > worst_pierce {
                worst_pierce = cells;
            }
        }
        out.worst_piercing_cells = if verts.is_empty() { 0.0 } else { worst_pierce };
        out
    }
}

/// Measure one field with all three covered extractors, appending three rows.
fn measure<F>(field: &F, rows: &mut Vec<Row>)
where
    F: ReferenceField + Sdf<Scalar = f64>,
{
    let (samples, shape) = Samples::of(field);
    let bound = bound_name(field.bound());
    let mut mesh = MeshBuffer::<f64>::new();
    isomesh::for_each_extractor!(f64, |name, extractor| {
        if COVERED.contains(&name) {
            mesh.reset();
            if let Err(e) = extractor.extract_into(field, &shape, samples.lo, samples.h, &mut mesh)
            {
                panic!("P-51: {name} could not extract {}: {e}", F::NAME);
            }
            let pierce = samples.piercing(&mesh.positions);
            let scan = samples.scan(&mesh.positions);
            let row = Row {
                field: F::NAME,
                bound,
                extractor: name,
                cell_size: samples.h,
                samples: samples.abs_d.len(),
                nonfinite: samples.nonfinite,
                vertices: mesh.vertex_count(),
                triangles: mesh.triangle_count(),
                pierced: pierce.pierced,
                worst_piercing: pierce.worst_cells,
                pierced_exhaustive: scan.pierced,
                worst_piercing_exhaustive: scan.worst_piercing_cells,
                untouched: scan.untouched,
                worst_untouched: scan.worst_untouched_cells,
            };
            println!(
                "  {:>16} {:>16}  verts {:>6}  pierced {:>6} ({:>7.3}/1k, worst {:>9.4} cells)  \
                 untouched {:>6} ({:>8.3}/1k, worst {:>8.4} cells)  exhaustive pierced {:>6} \
                 (worst {:>9.4} cells)",
                row.field,
                row.extractor,
                row.vertices,
                row.pierced,
                row.pierced_per_1k(),
                row.worst_piercing,
                row.untouched,
                row.untouched_per_1k(),
                row.worst_untouched,
                row.pierced_exhaustive,
                row.worst_piercing_exhaustive,
            );
            if pierce.pierced > 0 {
                println!(
                    "      pierced worst: vertex [{:.6}, {:.6}, {:.6}] sits {:.6} inside the ball \
                     of sample [{:.6}, {:.6}, {:.6}] — |d| = {:.6}, distance = {:.6}, \
                     violation = {:.6} cells",
                    pierce.vertex[0],
                    pierce.vertex[1],
                    pierce.vertex[2],
                    pierce.radius - pierce.distance,
                    pierce.sample[0],
                    pierce.sample[1],
                    pierce.sample[2],
                    pierce.radius,
                    pierce.distance,
                    pierce.worst_cells,
                );
            }
            println!(
                "      untouched worst: sample [{:.6}, {:.6}, {:.6}] has |d| = {:.6} and its \
                 closest-to-shell vertex [{:.6}, {:.6}, {:.6}] is {:.6} away — miss {:.6}, \
                 {:.6} cells",
                scan.sample[0],
                scan.sample[1],
                scan.sample[2],
                scan.radius,
                scan.nearest_shell[0],
                scan.nearest_shell[1],
                scan.nearest_shell[2],
                scan.shell_distance,
                scan.shell_distance - scan.radius,
                scan.worst_untouched_cells,
            );
            rows.push(row);
        }
    });
}

/// The `pierced_per_1k` of one extractor on one field, or zero if absent.
fn rate_of(rows: &[Row], field: &str, extractor: &str) -> f64 {
    rows.iter()
        .find(|r| r.field == field && r.extractor == extractor)
        .map_or(0.0, Row::pierced_per_1k)
}

/// The `untouched_per_1k` of one extractor on one field, or zero if absent.
fn untouched_of(rows: &[Row], field: &str, extractor: &str) -> f64 {
    rows.iter()
        .find(|r| r.field == field && r.extractor == extractor)
        .map_or(0.0, Row::untouched_per_1k)
}

fn main() {
    let prereg = isomesh::experiment!("P-51");

    common::experiment::run(prereg, |run| {
        let mut rows: Vec<Row> = Vec::new();

        // The control runs and is reported FIRST. An instrument that finds
        // piercing on an axis-aligned planar surface is measuring its own
        // tolerance, and every other number below would be void.
        println!("control — box_exact must report ZERO piercing for all three:");
        measure(&BoxExact::<f64>::canonical(), &mut rows);
        let control_zero = rows
            .iter()
            .filter(|r| r.field == "box_exact")
            .all(|r| r.pierced == 0);
        println!(
            "control_box_exact_zero = {control_zero}  ({})\n",
            if control_zero {
                "instrument sound; the rest of the table stands"
            } else {
                "VOID — the instrument is measuring its own tolerance"
            }
        );

        println!("the remaining four fields:");
        measure(&Sphere::<f64>::canonical(), &mut rows);
        measure(&Torus::<f64>::canonical(), &mut rows);
        measure(&ThinPlate::<f64>::canonical(), &mut rows);
        measure(&csg_difference::<f64>(), &mut rows);

        // C2 is a ratio between rows, so it is computed after every row exists
        // and written onto all three rows of the field it is about.
        let mc_total: u64 = rows
            .iter()
            .filter(|r| r.extractor == "marching_cubes")
            .map(|r| r.pierced)
            .sum();
        let mc_verts: usize = rows
            .iter()
            .filter(|r| r.extractor == "marching_cubes")
            .map(|r| r.vertices)
            .sum();
        let dc_total: u64 = rows
            .iter()
            .filter(|r| r.extractor == "dual_contouring")
            .map(|r| r.pierced)
            .sum();
        let dc_verts: usize = rows
            .iter()
            .filter(|r| r.extractor == "dual_contouring")
            .map(|r| r.vertices)
            .sum();
        let mc_all = 1000.0 * mc_total as f64 / mc_verts as f64;
        let dc_all = 1000.0 * dc_total as f64 / dc_verts as f64;
        let all_fields = ratio(dc_all, mc_all);

        for row in &rows {
            let field_ratio = ratio(
                rate_of(&rows, row.field, "dual_contouring"),
                rate_of(&rows, row.field, "marching_cubes"),
            );
            let untouched_ratio = ratio(
                untouched_of(&rows, row.field, "marching_cubes"),
                untouched_of(&rows, row.field, "dual_contouring"),
            );
            run.record(&[
                ("field", row.field.to_string()),
                ("extractor", row.extractor.to_string()),
                ("samples_per_axis", SAMPLES_PER_AXIS.to_string()),
                ("vertices", row.vertices.to_string()),
                ("samples", row.samples.to_string()),
                ("pierced", row.pierced.to_string()),
                ("pierced_per_1k", format!("{:.4}", row.pierced_per_1k())),
                ("worst_piercing_cells", format!("{:.6}", row.worst_piercing)),
                ("dc_over_mc_ratio", field_ratio),
                ("untouched", row.untouched.to_string()),
                ("untouched_per_1k", format!("{:.4}", row.untouched_per_1k())),
                (
                    "worst_untouched_cells",
                    format!("{:.6}", row.worst_untouched),
                ),
                (
                    "samples_probed_per_vertex",
                    SAMPLES_PROBED_PER_VERTEX.to_string(),
                ),
                ("control_box_exact_zero", control_zero.to_string()),
                ("threshold_cells", format!("{THRESHOLD_CELLS}")),
                // Extras. None of these gate anything.
                ("field_bound", row.bound.to_string()),
                ("triangles", row.triangles.to_string()),
                ("cell_size", format!("{:.6}", row.cell_size)),
                ("touch_search", String::from("exhaustive")),
                ("vertices_probed_per_sample", row.vertices.to_string()),
                ("pierced_exhaustive", row.pierced_exhaustive.to_string()),
                (
                    "pierced_exhaustive_per_1k",
                    format!("{:.4}", row.pierced_exhaustive_per_1k()),
                ),
                (
                    "worst_piercing_exhaustive_cells",
                    format!("{:.6}", row.worst_piercing_exhaustive),
                ),
                ("untouched_mc_over_dc", untouched_ratio),
                ("dc_over_mc_ratio_all_fields", all_fields.clone()),
                ("nonfinite_field_samples", row.nonfinite.to_string()),
            ]);
        }

        println!("\nC1  marching_cubes pierced_per_1k, worst field:");
        let mut c1 = 0.0_f64;
        for row in rows.iter().filter(|r| r.extractor == "marching_cubes") {
            let rate = row.pierced_per_1k();
            println!("      {:>16} {rate:.4}", row.field);
            if rate > c1 {
                c1 = rate;
            }
        }
        println!("    worst = {c1:.4} (gate: < 1)");
        println!(
            "C2  dc_over_mc_ratio, all five fields pooled: {all_fields} \
             (dc {dc_all:.4}/1k over {dc_verts} verts, mc {mc_all:.4}/1k over {mc_verts} verts, \
             gate: >= 20)"
        );
        println!("C3  untouched, every row:");
        for row in &rows {
            println!(
                "      {:>16} {:>16} untouched {:>6} ({:.4}/1k)",
                row.field,
                row.extractor,
                row.untouched,
                row.untouched_per_1k()
            );
        }
    });
}
