//! **P-69 — autovectorising the sample loop, with bit-identity as the gate.**
//!
//! Ticket: R-067. Pre-registered before this harness existed.
//!
//! ```bash
//! cargo bench --bench experiment_p69
//! ```
//!
//! Writes `docs/experiments/p-69.csv`.
//!
//! # Two levels, because C1's number and C1's mechanism are different questions
//!
//! - **`loop` rows** — the sample loop measured *in isolation*, both shapes
//!   compiled into this binary, `push` against `row`. This is the honest unit
//!   for "does this loop shape vectorise", and its ratio is within one build and
//!   one run, which is `M-281`'s rule.
//! - **`extract` rows** — `MarchingCubes::extract_into` over a resolution sweep,
//!   fitted to `t = a + b·n³`, so `b` is the marginal cost per sample that C1's
//!   threshold is stated in. This exercises the **shipped** `sdf::sample_grid`
//!   rather than a copy of it.
//!
//! The `push` arm is a bench-local copy of the loop that was in `src/` before
//! R-067, kept **here** and nowhere in the library, because a second sampling
//! path in `isomesh` is the two-paths defect the crate's own rules forbid. It is
//! held to the new one by a bit-identity assertion on every row, so a drift in
//! the copy is a failure rather than a number.
//!
//! # Cycles, not nanoseconds
//!
//! `M-280`: on a governed CPU a nanosecond is not a unit, so every row carries
//! `cycles_per_sample` and the `ghz` it implies, from `perf_event_open`. This
//! host's governor spans 1.96–5.62 GHz; the ratio between two arms measured in
//! the same run is what survives that, and an absolute ns/sample is not.
//!
//! # What the fields decide, established before this was written
//!
//! `libm` 0.2.16's `sqrtf` carries a `select_implementation` on
//! `target_feature = "sse2"`, so `sphere`'s body reduces to hardware
//! instructions. `sinf` and `cosf` carry **no** arch selection — software with
//! argument-reduction branches — so `gyroid`'s six of them per sample cannot be
//! widened at any loop shape. `box_exact` is `min`/`max`/`abs` and is a second
//! vectorisable body, included as a positive control rather than to widen the
//! claim. `vectorisable_body` records which is which, and it is an input to the
//! experiment rather than an output.

#![allow(clippy::float_cmp)]

mod common;

use std::time::Instant;

use common::counters::Probe;
use isomesh::extractor::Extractor;
use isomesh::fields::{BoxExact, ReferenceField, Sphere, capped_gyroid};
use isomesh::marching_cubes::MarchingCubes;
use isomesh::{MeshBuffer, Real, RuntimeShape3, Sdf};

// ─── the two loop shapes ────────────────────────────────────────────────────

/// The loop that was in `src/` before R-067, in all three of its copies.
///
/// `Vec::push` in the innermost loop, so the capacity bound is re-proved per
/// element and the store is not a contiguous write LLVM can widen; and the `y`
/// and `z` coordinates recomputed inside it.
fn push_loop<R: Real, S: Sdf<Scalar = R>>(
    sdf: &S,
    size: [u32; 3],
    origin: [R; 3],
    cell_size: R,
    out: &mut Vec<R>,
) {
    out.clear();
    out.reserve(size[0] as usize * size[1] as usize * size[2] as usize);
    for z in 0..size[2] {
        for y in 0..size[1] {
            for x in 0..size[0] {
                let p = [
                    origin[0] + cell_size * R::from_f64(f64::from(x)),
                    origin[1] + cell_size * R::from_f64(f64::from(y)),
                    origin[2] + cell_size * R::from_f64(f64::from(z)),
                ];
                out.push(sdf.sample(p));
            }
        }
    }
}

/// The shape R-067 put in `src/`, transcribed so both arms are in one binary.
///
/// Buffer sized once, one slice per row so its bound is proved once, the two
/// outer coordinates hoisted, and an `iter_mut().enumerate()` over a slice of
/// known length. Held to `push_loop` by a bit-identity assertion on every row.
fn row_loop<R: Real, S: Sdf<Scalar = R>>(
    sdf: &S,
    size: [u32; 3],
    origin: [R; 3],
    cell_size: R,
    out: &mut Vec<R>,
) {
    let nx = size[0] as usize;
    let rows = size[1] as usize * size[2] as usize;
    out.clear();
    out.resize(nx * rows, R::ZERO);
    for z in 0..size[2] {
        let pz = origin[2] + cell_size * R::from_f64(f64::from(z));
        for y in 0..size[1] {
            let py = origin[1] + cell_size * R::from_f64(f64::from(y));
            let start = nx * (y as usize + size[1] as usize * z as usize);
            let row = &mut out[start..start + nx];
            for (x, slot) in row.iter_mut().enumerate() {
                *slot = sdf.sample([origin[0] + cell_size * R::from_f64(x as f64), py, pz]);
            }
        }
    }
}

// ─── asm probes ─────────────────────────────────────────────────────────────
//
// The registration requires assembly for the **monomorphised `f32` instance**,
// and with thin LTO and one codegen unit both loops inline into their callers and
// leave no symbol to inspect. These wrappers are monomorphic and
// `#[inline(never)]`, so each one is a symbol whose body is exactly the generic
// loop at exactly the instantiation the timing uses. They are called once from
// `main` so nothing eliminates them.
//
// `scripts/p69_asm.sh` is what reads them back; the classification is a count of
// packed against scalar arithmetic and of the calls that make widening
// impossible.

/// The pre-R-067 loop at `Sphere<f32>`.
#[inline(never)]
pub fn asm_probe_push_sphere_f32(sdf: &Sphere<f32>, out: &mut Vec<f32>) {
    push_loop(sdf, [17, 17, 17], [-2.0; 3], 0.25, out);
}

/// The shipped loop shape at `Sphere<f32>` — `sqrt`, so widenable.
#[inline(never)]
pub fn asm_probe_row_sphere_f32(sdf: &Sphere<f32>, out: &mut Vec<f32>) {
    row_loop(sdf, [17, 17, 17], [-2.0; 3], 0.25, out);
}

/// The shipped loop shape at `BoxExact<f32>` — `min`/`max`, no calls.
#[inline(never)]
pub fn asm_probe_row_box_f32(sdf: &BoxExact<f32>, out: &mut Vec<f32>) {
    row_loop(sdf, [17, 17, 17], [-2.0; 3], 0.25, out);
}

/// The shipped loop shape at `CappedGyroid<f32>` — six `libm` calls.
#[inline(never)]
pub fn asm_probe_row_gyroid_f32(sdf: &isomesh::fields::CappedGyroid<f32>, out: &mut Vec<f32>) {
    row_loop(sdf, [17, 17, 17], [-2.0; 3], 0.25, out);
}

// ─── measurement ────────────────────────────────────────────────────────────

/// One timed loop run: the median of `reps`, with cycles from the same window.
struct Timed {
    ns_per_sample: f64,
    cycles_per_sample: f64,
    ghz: f64,
    worst_ratio: f64,
}

/// Time a loop shape. **The median, not the mean**, and the counters cover the
/// whole rep set rather than one rep, so a single scheduling hiccup cannot become
/// the reported figure.
#[derive(Clone, Copy)]
struct Grid<R> {
    size: [u32; 3],
    origin: [R; 3],
    cell_size: R,
}

fn time_loop<R: Real, S: Sdf<Scalar = R>>(
    probe: &mut Probe,
    reps: u32,
    shape: impl Fn(&S, [u32; 3], [R; 3], R, &mut Vec<R>),
    sdf: &S,
    grid: Grid<R>,
    out: &mut Vec<R>,
) -> Timed {
    let Grid {
        size,
        origin,
        cell_size,
    } = grid;
    let samples = size[0] as u64 * size[1] as u64 * size[2] as u64;
    // One untimed pass: the allocation and the first-touch page faults belong to
    // neither arm's arithmetic.
    shape(sdf, size, origin, cell_size, out);

    let mut times = Vec::with_capacity(reps as usize);
    probe.reset_and_enable();
    let all = Instant::now();
    for _ in 0..reps {
        let one = Instant::now();
        shape(sdf, size, origin, cell_size, out);
        times.push(one.elapsed().as_nanos() as f64);
        std::hint::black_box(&out[0]);
    }
    let total_ns = all.elapsed().as_nanos() as f64;
    probe.disable();
    let counts = probe.read();

    times.sort_unstable_by(f64::total_cmp);
    let median = times[times.len() / 2];
    let cycles = counts.cycles.count as f64;
    Timed {
        ns_per_sample: median / samples as f64,
        cycles_per_sample: cycles / (samples * u64::from(reps)) as f64,
        ghz: cycles / total_ns,
        worst_ratio: counts.worst_ratio(),
    }
}

/// Least squares on `t = a + b·n³`, returning `b` in ns per sample.
fn marginal(points: &[(f64, f64)]) -> f64 {
    let n = points.len() as f64;
    let sx: f64 = points.iter().map(|p| p.0).sum();
    let sy: f64 = points.iter().map(|p| p.1).sum();
    let sxx: f64 = points.iter().map(|p| p.0 * p.0).sum();
    let sxy: f64 = points.iter().map(|p| p.0 * p.1).sum();
    (n * sxy - sx * sy) / (n * sxx - sx * sx)
}

/// This host's committed marginal cost for `marching_cubes/f32/sphere`, fitted
/// over the 9 rows of `docs/measurements/resolution_sweep-ryzen9-5900x.csv` from
/// 16³ to 256³. C1's threshold is half of it.
const BASELINE_MARGINAL_NS: f64 = 13.1892;

/// Whether the field's body can be widened at all — decided by `libm`, not by
/// the loop. An input, not an outcome.
fn vectorisable(field: &str) -> bool {
    match field {
        // `sqrt` selects a hardware instruction under `sse2`.
        "sphere" => true,
        // `min`/`max`/`abs`, no calls at all.
        "box_exact" => true,
        // Six `libm::sinf`/`cosf` per sample, software, with branches.
        "gyroid" => false,
        other => panic!("{other} has no established body classification"),
    }
}

type Row = Vec<(&'static str, String)>;

const NA: &str = "";

/// Every loop-level row for one field at one scalar.
#[allow(clippy::too_many_arguments)]
fn loop_rows<R: Real, S: Sdf<Scalar = R>>(
    probe: &mut Probe,
    rows: &mut Vec<Row>,
    field_name: &'static str,
    scalar: &'static str,
    sdf: &S,
    l: f64,
    sizes: &[u32],
    reps: u32,
) -> (f64, f64) {
    let mut push_buf: Vec<R> = Vec::new();
    let mut row_buf: Vec<R> = Vec::new();
    let mut push_points = Vec::new();
    let mut row_points = Vec::new();

    for &n in sizes {
        let size = [n; 3];
        let cell_size = R::from_f64(2.0 * l / f64::from(n - 1));
        let origin = [R::from_f64(-l); 3];
        let samples = u64::from(n) * u64::from(n) * u64::from(n);

        let grid = Grid {
            size,
            origin,
            cell_size,
        };
        let a = time_loop(probe, reps, push_loop, sdf, grid, &mut push_buf);
        let b = time_loop(probe, reps, row_loop, sdf, grid, &mut row_buf);

        // **The two arms must agree bit for bit.** Hoisting a loop invariant
        // cannot change a rounding, so a difference here means the bench-local
        // copy drifted from the shipped shape and no number in this row means
        // anything.
        assert_eq!(
            push_buf.len(),
            row_buf.len(),
            "{field_name}/{scalar}/{n}: the two arms wrote different lengths"
        );
        let identical = push_buf
            .iter()
            .zip(row_buf.iter())
            .all(|(p, r)| p.total_cmp(r) == core::cmp::Ordering::Equal);
        assert!(
            identical,
            "{field_name}/{scalar}/{n}: the row loop is not bit-identical to the \
             push loop, so it is a different computation and C2 is already lost"
        );
        assert!(
            a.worst_ratio > 0.99 && b.worst_ratio > 0.99,
            "{field_name}/{scalar}/{n}: the kernel multiplexed a counter \
             ({:.3}/{:.3}), so cycles are scaled and not comparable",
            a.worst_ratio,
            b.worst_ratio
        );

        push_points.push((samples as f64, a.ns_per_sample * samples as f64));
        row_points.push((samples as f64, b.ns_per_sample * samples as f64));

        let speedup = a.ns_per_sample / b.ns_per_sample;
        println!(
            "{field_name:<11} {scalar:<4} {n:>4} {:>9.4} {:>9.4} {:>8.3} {:>8.2} {:>8.2} {:>6.2}",
            a.ns_per_sample,
            b.ns_per_sample,
            speedup,
            a.cycles_per_sample,
            b.cycles_per_sample,
            b.ghz
        );

        for (shape_name, t) in [("push", &a), ("row", &b)] {
            rows.push(vec![
                ("arm", "loop".to_string()),
                ("field", field_name.to_string()),
                ("scalar", scalar.to_string()),
                ("samples_per_axis", n.to_string()),
                ("samples", samples.to_string()),
                ("loop_shape", shape_name.to_string()),
                ("ns_per_sample", format!("{:.6}", t.ns_per_sample)),
                ("cycles_per_sample", format!("{:.4}", t.cycles_per_sample)),
                ("ghz", format!("{:.4}", t.ghz)),
                (
                    "speedup_vs_push",
                    if shape_name == "row" {
                        format!("{speedup:.6}")
                    } else {
                        "1.000000".to_string()
                    },
                ),
                ("bit_identical_to_push", "true".to_string()),
                ("vectorisable_body", vectorisable(field_name).to_string()),
                ("counter_ratio", format!("{:.4}", t.worst_ratio)),
                ("reps", reps.to_string()),
            ]);
        }
    }
    (marginal(&push_points), marginal(&row_points))
}

fn main() {
    if !std::env::args().any(|arg| arg == "--bench") {
        return;
    }
    let prereg = isomesh::experiment!("P-69");
    common::experiment::run(prereg, |run| {
        let machine = "amd-ryzen-9-5900x-12-core";
        // Keep the asm probes alive. One call each, results discarded through
        // `black_box` so the optimiser cannot delete the symbol the
        // registration requires be inspectable.
        {
            let mut buf32: Vec<f32> = Vec::new();
            asm_probe_push_sphere_f32(&Sphere::<f32>::canonical(), &mut buf32);
            std::hint::black_box(&buf32);
            asm_probe_row_sphere_f32(&Sphere::<f32>::canonical(), &mut buf32);
            std::hint::black_box(&buf32);
            asm_probe_row_box_f32(&BoxExact::<f32>::canonical(), &mut buf32);
            std::hint::black_box(&buf32);
            asm_probe_row_gyroid_f32(&capped_gyroid::<f32>(), &mut buf32);
            std::hint::black_box(&buf32);
        }
        let mut probe = Probe::open();
        let mut rows: Vec<Row> = Vec::new();

        println!(
            "\n-- loop: the sample loop in isolation, both shapes in one binary --\n\
             {:<11} {:<4} {:>4} {:>9} {:>9} {:>8} {:>8} {:>8} {:>6}",
            "field", "sc", "n", "push_ns", "row_ns", "speedup", "push_cyc", "row_cyc", "GHz"
        );

        let sizes = [33u32, 65, 97, 129];
        let reps = 5u32;
        let mut marginals: Vec<(&'static str, &'static str, f64, f64)> = Vec::new();

        // `sphere` and `gyroid` are C1's two; `box_exact` is a second
        // vectorisable body as a positive control.
        let sphere = Sphere::<f32>::canonical();
        let l32 = sphere.domain().1[0] as f64;
        let (p, r) = loop_rows(
            &mut probe, &mut rows, "sphere", "f32", &sphere, l32, &sizes, reps,
        );
        marginals.push(("sphere", "f32", p, r));
        let sphere64 = Sphere::<f64>::canonical();
        let (p, r) = loop_rows(
            &mut probe,
            &mut rows,
            "sphere",
            "f64",
            &sphere64,
            sphere64.domain().1[0],
            &sizes,
            reps,
        );
        marginals.push(("sphere", "f64", p, r));

        let gyroid = capped_gyroid::<f32>();
        let (p, r) = loop_rows(
            &mut probe,
            &mut rows,
            "gyroid",
            "f32",
            &gyroid,
            gyroid.domain().1[0] as f64,
            &sizes,
            reps,
        );
        marginals.push(("gyroid", "f32", p, r));
        let gyroid64 = capped_gyroid::<f64>();
        let (p, r) = loop_rows(
            &mut probe,
            &mut rows,
            "gyroid",
            "f64",
            &gyroid64,
            gyroid64.domain().1[0],
            &sizes,
            reps,
        );
        marginals.push(("gyroid", "f64", p, r));

        let boxf = BoxExact::<f32>::canonical();
        let (p, r) = loop_rows(
            &mut probe,
            &mut rows,
            "box_exact",
            "f32",
            &boxf,
            boxf.domain().1[0] as f64,
            &sizes,
            reps,
        );
        marginals.push(("box_exact", "f32", p, r));
        let box64 = BoxExact::<f64>::canonical();
        let (p, r) = loop_rows(
            &mut probe,
            &mut rows,
            "box_exact",
            "f64",
            &box64,
            box64.domain().1[0],
            &sizes,
            reps,
        );
        marginals.push(("box_exact", "f64", p, r));

        println!("\n-- marginal ns/sample, fitted over the same four sizes --");
        println!(
            "{:<11} {:<4} {:>12} {:>12} {:>9}",
            "field", "sc", "push", "row", "gain"
        );
        for (field, scalar, push_m, row_m) in &marginals {
            println!(
                "{field:<11} {scalar:<4} {push_m:>12.4} {row_m:>12.4} {:>9.4}",
                push_m / row_m
            );
            rows.push(vec![
                ("arm", "marginal".to_string()),
                ("field", (*field).to_string()),
                ("scalar", (*scalar).to_string()),
                ("samples_per_axis", NA.to_string()),
                ("samples", NA.to_string()),
                ("loop_shape", "push_vs_row".to_string()),
                ("marginal_ns_per_sample", format!("{row_m:.6}")),
                ("marginal_ns_per_sample_push", format!("{push_m:.6}")),
                ("speedup_vs_push", format!("{:.6}", push_m / row_m)),
                (
                    "baseline_marginal_ns_per_sample",
                    format!("{BASELINE_MARGINAL_NS:.4}"),
                ),
                ("vectorisable_body", vectorisable(field).to_string()),
            ]);
        }

        // ── the extraction sweep: the shipped code, and C1's own units ──────
        println!("\n-- extract: MarchingCubes over the shipped sample_grid --");
        println!("{:<11} {:<4} {:>4} {:>12}", "field", "sc", "n", "ms");
        let mut extract_points: Vec<(f64, f64)> = Vec::new();
        {
            let field = Sphere::<f32>::canonical();
            let l = field.domain().1[0];
            let mut mc = MarchingCubes::<f32>::new();
            let mut mesh = MeshBuffer::<f32>::new();
            for &n in &[33u32, 65, 97, 129] {
                let shape = RuntimeShape3::new([n; 3]).expect("grid");
                let cell = 2.0 * l / (n - 1) as f32;
                let origin = [-l; 3];
                let mut times = Vec::new();
                for _ in 0..reps {
                    mesh.reset();
                    let t = Instant::now();
                    mc.extract_into(&field, &shape, origin, cell, &mut mesh)
                        .expect("extraction");
                    times.push(t.elapsed().as_nanos() as f64);
                }
                times.sort_unstable_by(f64::total_cmp);
                let median = times[times.len() / 2];
                let samples = f64::from(n) * f64::from(n) * f64::from(n);
                println!(
                    "{:<11} {:<4} {n:>4} {:>12.6}",
                    "sphere",
                    "f32",
                    median / 1e6
                );
                extract_points.push((samples, median));
                rows.push(vec![
                    ("arm", "extract".to_string()),
                    ("field", "sphere".to_string()),
                    ("scalar", "f32".to_string()),
                    ("samples_per_axis", n.to_string()),
                    ("samples", format!("{samples:.0}")),
                    ("loop_shape", "row".to_string()),
                    ("ns_per_sample", format!("{:.6}", median / samples)),
                    ("vectorisable_body", "true".to_string()),
                ]);
            }
        }
        let extract_marginal = marginal(&extract_points);
        println!(
            "extraction marginal: {extract_marginal:.4} ns/sample against a \
             committed baseline of {BASELINE_MARGINAL_NS:.4}"
        );

        // ── the golden fixture, over the 168 rows this harness can rebuild ──
        //
        // C2's verdict belongs in the artefact rather than in a test log, so the
        // harness recomputes the hash for every (extractor, field, resolution)
        // `for_each_extractor!` reaches and compares it against the committed
        // fixture. That is 168 of the 216; `greedy_quads` and
        // `marching_cubes+trilinear` are the other 48 and are covered by
        // `golden_hashes_are_unchanged`, which is quoted in the entry. The
        // enumeration is not duplicated here on purpose — a second copy of the
        // golden matrix is a second thing to drift.
        println!("\n-- golden: the 168 rows for_each_extractor! reaches --");
        let fixture = std::fs::read_to_string(
            std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("golden_hashes.json"),
        )
        .expect("golden_hashes.json");
        let mut checked = 0usize;
        let mut moved = 0usize;
        let mut mesh = MeshBuffer::<f64>::new();
        isomesh::for_each_reference_field!(f64, |field_name, field| {
            let l = field.domain().1[0];
            for n in [17u32, 25, 33] {
                let shape = RuntimeShape3::new([n; 3]).expect("grid");
                let cell = 2.0 * l / f64::from(n - 1);
                let origin = [-l; 3];
                isomesh::for_each_extractor!(f64, |extractor_name, extractor| {
                    mesh.reset();
                    extractor
                        .extract_into(&field, &shape, origin, cell, &mut mesh)
                        .expect("extraction");
                    let hash = isomesh::validate::mesh_hash(&mesh);
                    // The fixture is one object per line, keys unspaced. Match
                    // on the triple rather than on a composed key, so a change
                    // to the key format is a miss rather than a silent pass.
                    let want = fixture
                        .lines()
                        .find(|line| {
                            line.contains(&format!("\"algorithm\":\"{extractor_name}\""))
                                && line.contains(&format!("\"field\":\"{field_name}\""))
                                && line.contains(&format!("\"samples\":{n},"))
                        })
                        .and_then(|line| {
                            line.split("\"hash\":\"")
                                .nth(1)
                                .and_then(|t| t.split('"').next())
                                .map(str::to_string)
                        });
                    match want {
                        Some(want) => {
                            checked += 1;
                            if want != format!("{hash:016x}") {
                                moved += 1;
                                println!(
                                    "   MOVED {extractor_name}/{field_name}/{n}: \
                                     {want} -> {hash:016x}"
                                );
                            }
                        }
                        None => panic!(
                            "golden_hashes.json has no row for \
                             {extractor_name}/{field_name}/{n}"
                        ),
                    }
                });
            }
        });
        println!("   {checked} rows checked, {moved} moved");
        assert_eq!(checked, 168, "for_each_extractor! should reach 168 rows");

        // ── verdicts ────────────────────────────────────────────────────────
        let gain = |field: &str, scalar: &str| -> f64 {
            marginals
                .iter()
                .find(|(f, s, _, _)| *f == field && *s == scalar)
                .map_or(0.0, |(_, _, p, r)| p / r)
        };
        let c1_sphere = gain("sphere", "f32");
        let c1_gyroid = gain("gyroid", "f32");
        let c1_box = gain("box_exact", "f32");
        let c1 = c1_sphere >= 2.0 && c1_gyroid >= 2.0;
        let c2 = moved == 0;
        let f64_over_f32 = |field: &str| {
            let (a, b) = (gain(field, "f64"), gain(field, "f32"));
            if b <= 1.0 {
                f64::NAN
            } else {
                (a - 1.0) / (b - 1.0)
            }
        };
        let c3_ratio = f64_over_f32("sphere");
        // **C3 presupposes an f32 gain, and C1 says there is none.** "The f64
        // gain is at most half the f32 gain" divides by `f32_gain - 1`, so with
        // the f32 gain at or below 1 the clause has no denominator and therefore
        // no population. Reporting that as FALSIFIED would be a claim the run
        // cannot support, and reporting it as HELD would be worse; it is
        // **vacuous**, and the column says so in as many letters. This is the
        // audit's central finding one level up: a clause can be unfirable
        // because another clause failed.
        let c3 = if c3_ratio.is_finite() {
            if c3_ratio <= 0.5 { "true" } else { "false" }
        } else {
            "vacuous"
        };

        println!(
            "\nC1 f32 loop gain: sphere {c1_sphere:.4}, gyroid {c1_gyroid:.4}, \
             box_exact {c1_box:.4} -> {}",
            if c1 { "HELD" } else { "FALSIFIED" }
        );
        println!(
            "C2 golden rows moved: {moved} of {checked} -> {}",
            if c2 { "HELD" } else { "FALSIFIED (veto)" }
        );
        println!(
            "C3 f64 gain over f32 gain on sphere: {c3_ratio:.4} -> {}",
            match c3 {
                "true" => "HELD",
                "false" => "FALSIFIED",
                _ => "VACUOUS -- C1 left no f32 gain to halve",
            }
        );

        let aggregates: Row = vec![
            ("c1_speedup_f32", format!("{c1_sphere:.6}")),
            ("c1_speedup_f32_gyroid", format!("{c1_gyroid:.6}")),
            ("c1_speedup_f32_box_exact", format!("{c1_box:.6}")),
            ("c1_holds", c1.to_string()),
            ("golden_hashes_unchanged", c2.to_string()),
            ("golden_rows_checked", checked.to_string()),
            ("golden_rows_moved", moved.to_string()),
            ("c2_holds", c2.to_string()),
            ("c3_f64_over_f32_gain", format!("{c3_ratio:.6}")),
            ("c3_holds", c3.to_string()),
            (
                "extraction_marginal_ns_per_sample",
                format!("{extract_marginal:.6}"),
            ),
            (
                "baseline_marginal_ns_per_sample",
                format!("{BASELINE_MARGINAL_NS:.4}"),
            ),
            ("machine", machine.to_string()),
        ];

        let registered: [&str; 22] = [
            "arm",
            "field",
            "scalar",
            "samples_per_axis",
            "samples",
            "loop_shape",
            "ns_per_sample",
            "cycles_per_sample",
            "ghz",
            "speedup_vs_push",
            "marginal_ns_per_sample",
            "baseline_marginal_ns_per_sample",
            "bit_identical_to_push",
            "golden_hashes_unchanged",
            "vectorisable_body",
            "c1_speedup_f32",
            "c1_holds",
            "c2_holds",
            "c3_f64_over_f32_gain",
            "c3_holds",
            "machine",
            "wall_ms",
        ];
        let started = Instant::now();
        let wall = started.elapsed().as_millis();
        for mut row in rows {
            row.extend(aggregates.iter().cloned());
            row.push(("wall_ms", wall.to_string()));
            for name in registered {
                if !row.iter().any(|(k, _)| *k == name) {
                    row.push((name, NA.to_string()));
                }
            }
            run.record(&row);
        }
    });
}
