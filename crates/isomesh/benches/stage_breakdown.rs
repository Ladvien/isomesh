//! M-003 — where the time actually goes, before anyone optimises anything.
//!
//! ```bash
//! cargo bench --bench stage_breakdown
//! ```
//!
//! Writes `docs/measurements/stage_breakdown.csv`.
//!
//! # The question, and the published figure that motivates it
//!
//! A meshing paper reports the *contour*. A consumer pays for a **usable mesh**,
//! which is the contour plus everything that has to happen before a renderer or
//! a physics engine will accept it. The one published comparison this project
//! has of the two puts contouring at 68 ms against halfedge construction at
//! 58 ms — so **the contour was 54% of the job**, and optimising it alone could
//! have bought at most a 1.85x speedup on a workload everyone describes as
//! "meshing".
//!
//! Nothing in this repository knew its own ratio. Every timing here — the
//! shootout, the resolution sweep, `extract` — measures the contour and stops.
//! This measures the rest.
//!
//! # The four stages
//!
//! | stage | what it is | who needs it |
//! |---|---|---|
//! | contour | the extractor | everyone |
//! | normals | [`normals::recompute`] with [`AreaWeightedFaces`] | anyone shading a mesh whose field has no analytic gradient |
//! | weld | [`weld::Welder`] | anyone who needs shared vertices — colliders, simplification, adjacency |
//! | collider | [`collider::readiness`] | anyone handing the mesh to a physics engine |
//!
//! Normals are the one that needs justifying: every extractor here already
//! writes normals from [`Sdf::gradient`], so this stage is **not** on the
//! critical path for the seven reference fields, which all override it. It is
//! the cost a consumer pays whose field does *not* — a sampled volume, an
//! imported scan — and it is measured with the geometric strategy because that
//! is the one such a consumer is forced onto.
//!
//! # What "upload" is, and why there is no bar for it
//!
//! M-003 names five stages and the fifth is upload. There is no bar for it here
//! and that is not an omission:
//!
//! - The CPU half is **a move**. `bevy_isomesh::MeshBuilder` writes extraction
//!   output straight into the arrays a Bevy `Mesh` owns, and `into_mesh`
//!   transfers them rather than copying. There is no pass to time.
//! - The GPU half needs a device, which `crates/isomesh` does not have and must
//!   not acquire — rule 2. It belongs to GPU-001 and later.
//!
//! Stating that is better than inventing a stage: a bar labelled "upload"
//! measuring a `Vec` move would imply a cost that is not there.

use std::fmt::Write as _;
use std::fs;
use std::hint::black_box;
use std::path::PathBuf;
use std::time::{Duration, Instant};

use isomesh::collider;
use isomesh::fields::ReferenceField;
use isomesh::marching_cubes::MarchingCubes;
use isomesh::normals::{self, NormalStrategy};
use isomesh::validate::ValidateConfig;
use isomesh::weld::Welder;
use isomesh::{MeshBuffer, Real, RuntimeShape3, Sdf};

type Scalar = f64;

/// Resolutions to break down. Two, because the *ratio* is the question and a
/// ratio that holds at one grid and not another is the more interesting answer.
const RESOLUTIONS: [u32; 2] = [33, 65];

const WARMUP_RUNS: u32 = 2;
const TIMED_RUNS: u32 = 7;

#[derive(Clone, Copy, PartialEq, Eq)]
enum Stage {
    Contour,
    Normals,
    Weld,
    Collider,
}

impl Stage {
    const ALL: [Self; 4] = [Self::Contour, Self::Normals, Self::Weld, Self::Collider];

    fn name(self) -> &'static str {
        match self {
            Self::Contour => "contour",
            Self::Normals => "normals",
            Self::Weld => "weld",
            Self::Collider => "collider",
        }
    }
}

struct Row {
    field: &'static str,
    samples: u32,
    triangles: usize,
    /// Median milliseconds per stage, in [`Stage::ALL`] order.
    ms: [f64; 4],
}

impl Row {
    fn total(&self) -> f64 {
        self.ms.iter().sum()
    }
    fn share(&self, stage: usize) -> f64 {
        let total = self.total();
        if total > 0.0 {
            100.0 * self.ms[stage] / total
        } else {
            0.0
        }
    }
}

fn median(times: &mut [Duration]) -> f64 {
    times.sort_unstable();
    times[times.len() / 2].as_secs_f64() * 1000.0
}

fn main() {
    // The same guard `resolution_sweep` and `shootout` use: `cargo test
    // --all-targets` re-selects bench targets even with `test = false`, and a
    // debug-build run would overwrite the committed CSV with numbers from
    // whichever machine ran the tests.
    if !std::env::args().any(|arg| arg == "--bench") {
        return;
    }

    println!("stage breakdown — contour / normals / weld / collider, f64, one process\n");
    let mut rows = Vec::new();
    isomesh::for_each_reference_field!(Scalar, |name, field| {
        for samples in RESOLUTIONS {
            if let Some(row) = measure(name, &field, samples) {
                print_row(&row);
                rows.push(row);
            }
        }
    });

    let path = write_csv(&rows);
    println!("\nwrote {}", path.display());
    report(&rows);
}

fn measure<F: Sdf<Scalar = Scalar> + ReferenceField>(
    field_name: &'static str,
    field: &F,
    samples: u32,
) -> Option<Row> {
    let (lo, hi) = field.domain();
    let h = (hi[0] - lo[0]) / Scalar::from_f64(f64::from(samples - 1));
    let shape = RuntimeShape3::new([samples; 3]).ok()?;
    let cfg = ValidateConfig::from_cell_size(h.as_f32().into()).ok()?;

    let mut ms = [0.0f64; 4];

    // Contour. Timed on its own buffer so the later stages start from an
    // identical mesh every run rather than from whatever the previous stage
    // left behind.
    let mut times = Vec::with_capacity((WARMUP_RUNS + TIMED_RUNS) as usize);
    let mut mesh = MeshBuffer::<Scalar>::new();
    for _ in 0..(WARMUP_RUNS + TIMED_RUNS) {
        mesh.reset();
        let start = Instant::now();
        MarchingCubes::<Scalar>::new()
            .extract(field, &shape, lo, h, &mut mesh)
            .ok();
        times.push(start.elapsed());
        black_box(mesh.triangle_count());
    }
    if mesh.triangle_count() == 0 {
        return None;
    }
    times.drain(..WARMUP_RUNS as usize);
    ms[0] = median(&mut times);
    let triangles = mesh.triangle_count();

    // Every later stage is timed from a fresh copy of the contour's output.
    // Welding in particular *mutates* the mesh, so timing it twice in a row on
    // the same buffer would measure a second weld of an already-welded mesh --
    // which is a different and much cheaper operation.
    let pristine = mesh.clone();

    let mut times = Vec::with_capacity((WARMUP_RUNS + TIMED_RUNS) as usize);
    for _ in 0..(WARMUP_RUNS + TIMED_RUNS) {
        let mut work = pristine.clone();
        let start = Instant::now();
        normals::recompute(&mut work, field, NormalStrategy::AreaWeightedFaces).ok();
        times.push(start.elapsed());
        black_box(work.normals.len());
    }
    times.drain(..WARMUP_RUNS as usize);
    ms[1] = median(&mut times);

    let mut times = Vec::with_capacity((WARMUP_RUNS + TIMED_RUNS) as usize);
    for _ in 0..(WARMUP_RUNS + TIMED_RUNS) {
        let mut work = pristine.clone();
        let start = Instant::now();
        Welder::<Scalar>::new()
            .weld(&mut work, isomesh::weld::epsilon_for(h))
            .ok();
        times.push(start.elapsed());
        black_box(work.positions.len());
    }
    times.drain(..WARMUP_RUNS as usize);
    ms[2] = median(&mut times);

    // The collider check runs on the **welded** mesh, because that is the only
    // mesh anyone would hand a physics engine: M-96 records that a validity
    // measurement on unwelded output is meaningless, and G-005's contract is
    // welded, manifold, correctly wound.
    let mut welded = pristine.clone();
    Welder::<Scalar>::new()
        .weld(&mut welded, isomesh::weld::epsilon_for(h))
        .ok();
    let mut times = Vec::with_capacity((WARMUP_RUNS + TIMED_RUNS) as usize);
    for _ in 0..(WARMUP_RUNS + TIMED_RUNS) {
        let start = Instant::now();
        let ready = collider::readiness(&welded, &cfg);
        times.push(start.elapsed());
        black_box(ready.is_usable());
    }
    times.drain(..WARMUP_RUNS as usize);
    ms[3] = median(&mut times);

    Some(Row {
        field: field_name,
        samples,
        triangles,
        ms,
    })
}

fn print_row(row: &Row) {
    print!(
        "{:<16} {:>3}^3 {:>7} tris  ",
        row.field, row.samples, row.triangles
    );
    for (i, stage) in Stage::ALL.iter().enumerate() {
        print!(
            "{} {:>8.3}ms ({:>4.1}%)  ",
            stage.name(),
            row.ms[i],
            row.share(i)
        );
    }
    println!("total {:>8.3}ms", row.total());
}

fn write_csv(rows: &[Row]) -> PathBuf {
    let mut csv = String::from(
        "field,samples,triangles,contour_ms,normals_ms,weld_ms,collider_ms,total_ms,contour_share\n",
    );
    for row in rows {
        let _ = writeln!(
            csv,
            "{},{},{},{:.6},{:.6},{:.6},{:.6},{:.6},{:.4}",
            row.field,
            row.samples,
            row.triangles,
            row.ms[0],
            row.ms[1],
            row.ms[2],
            row.ms[3],
            row.total(),
            row.share(0),
        );
    }
    let path = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../docs/measurements")
        .join("stage_breakdown.csv");
    let _ = fs::create_dir_all(path.parent().unwrap_or(&path));
    let _ = fs::write(&path, csv);
    path
}

/// The headline: what fraction of a usable mesh is the contour?
fn report(rows: &[Row]) {
    if rows.is_empty() {
        return;
    }
    println!("\ncontour as a share of a usable mesh:");
    let mut shares: Vec<f64> = Vec::new();
    for samples in RESOLUTIONS {
        let at: Vec<&Row> = rows.iter().filter(|r| r.samples == samples).collect();
        if at.is_empty() {
            continue;
        }
        let lo = at.iter().map(|r| r.share(0)).fold(f64::MAX, f64::min);
        let hi = at.iter().map(|r| r.share(0)).fold(0.0f64, f64::max);
        let mean = at.iter().map(|r| r.share(0)).sum::<f64>() / at.len() as f64;
        println!("  {samples:>3}^3   {mean:>5.1}%  (range {lo:.1}% .. {hi:.1}%)");
        shares.extend(at.iter().map(|r| r.share(0)));
    }

    // Every stage's share, so the answer is not just about the contour.
    println!("\nmean share per stage, all fields and resolutions:");
    for (i, stage) in Stage::ALL.iter().enumerate() {
        let mean = rows.iter().map(|r| r.share(i)).sum::<f64>() / rows.len() as f64;
        let bar = "#".repeat((mean / 2.0).round() as usize);
        println!("  {:<9} {mean:>5.1}%  {bar}", stage.name());
    }

    let overall = shares.iter().sum::<f64>() / shares.len() as f64;
    println!(
        "\nthe published comparison this ticket cites had the contour at 54% of a\n\
         usable mesh. here it is {overall:.1}%. optimising the contour alone can\n\
         therefore buy at most {:.2}x, and that ceiling is the number to know\n\
         before starting.",
        100.0 / (100.0 - overall)
    );
}
