//! **P-40 — the active-cell test is one bit, so 64 cells decide at once.**
//!
//! Ticket: R-039. Pre-registered in the commit before the mechanism landed.
//!
//! ```bash
//! ISOMESH_BLESS=1 cargo bench --bench experiment_p40   # at the parent commit
//! cargo bench --bench experiment_p40                   # after the mechanism lands
//! ```
//!
//! Writes `docs/experiments/p-40.csv`.
//!
//! # Two arms, measured two different ways, and why
//!
//! **The stage arm is self-contained.** Clause one is about the *predicate*, and
//! a predicate is small enough to hold both versions side by side in one
//! process: this bench samples the field into its own buffer, laid out exactly
//! as `DualMesher` lays out `values` (row stride `size[0] | 1`, A-024), and then
//! runs the scalar eight-corner gather and the packed word test over it. They
//! must agree on the active-cell *list*, not merely the count — an order
//! disagreement is the one failure that would silently change every index.
//!
//! **The extractor arm cannot be.** One binary contains one mesher, so the
//! before-and-after of a whole extraction spans two commits. The baseline is
//! therefore an input, produced once at the parent commit under
//! `ISOMESH_BLESS=1` and committed to
//! `docs/measurements/p40-baseline.csv` — the same idiom `golden_hashes.json`
//! uses, and for the same reason. Without the file this bench **fails** rather
//! than reporting one arm: a speedup with no baseline is not a measurement.
//!
//! The baseline carries a `mesh_hash` per row as well as a time, so clause three
//! is checked on exactly the grids that were timed rather than inferred from the
//! golden suite running on different ones.

mod common;

use std::collections::BTreeMap;
use std::path::PathBuf;
use std::time::Instant;

use isomesh::dual_contouring::DualContouring;
use isomesh::fields::{ReferenceField, Sphere};
use isomesh::surface_nets::SurfaceNets;
use isomesh::validate::mesh_hash;
use isomesh::{MeshBuffer, RuntimeShape3};

const REPS: usize = 5;

/// The grids the clauses name, plus the two that bracket them for context.
const SIZES: [u32; 3] = [64, 128, 256];

/// Where the parent commit's numbers live.
fn baseline_path() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../docs/measurements/p40-baseline.csv")
}

/// `field,samples,extractor` → `(ns_per_sample, mesh_hash)`.
type Baseline = BTreeMap<(String, u32, String), (f64, u64)>;

fn read_baseline() -> Option<Baseline> {
    let text = std::fs::read_to_string(baseline_path()).ok()?;
    let mut map = Baseline::new();
    for line in text
        .lines()
        .filter(|l| !l.starts_with('#') && !l.is_empty())
    {
        let cols: Vec<&str> = line.split(',').collect();
        if cols.len() != 5 || cols[0] == "field" {
            continue;
        }
        let samples = cols[1].parse().ok()?;
        let ns = cols[3].parse().ok()?;
        let hash = cols[4].parse().ok()?;
        map.insert(
            (cols[0].to_string(), samples, cols[2].to_string()),
            (ns, hash),
        );
    }
    Some(map)
}

/// The two fields, by the name they carry in the CSV.
///
/// `sphere` is the canonical unit sphere on `[-2, 2]³`, surface and all.
/// `sphere_surface_free` is the same sphere sampled a long way from itself, so
/// no corner is ever inside and the active path never runs — clause one's field,
/// and the one A-024 used for the same reason.
fn domain_of(field: &str, n: u32) -> ([f64; 3], f64) {
    if field == "sphere_surface_free" {
        ([10.0; 3], 4.0 / f64::from(n - 1))
    } else {
        let (lo, hi) = Sphere::<f64>::canonical().domain();
        (lo, (hi[0] - lo[0]) / f64::from(n - 1))
    }
}

/// `DualMesher`'s own layout, mirrored so the stage arm measures the same
/// addressing the crate does. `row` is `size[0] | 1` for A-024's reason.
struct Grid {
    values: Vec<f64>,
    row: usize,
    size: [u32; 3],
}

impl Grid {
    fn sample(field: &str, n: u32) -> Self {
        let (origin, cell) = domain_of(field, n);
        let sphere = Sphere::<f64>::canonical();
        let size = [n; 3];
        let row = n as usize | 1;
        let mut values = Vec::with_capacity(row * n as usize * n as usize);
        for z in 0..n {
            for y in 0..n {
                for x in 0..n {
                    values.push(isomesh::Sdf::sample(
                        &sphere,
                        [
                            origin[0] + cell * f64::from(x),
                            origin[1] + cell * f64::from(y),
                            origin[2] + cell * f64::from(z),
                        ],
                    ));
                }
                values.resize(values.len() + (row - n as usize), 0.0);
            }
        }
        Self { values, row, size }
    }

    #[inline]
    fn at(&self, x: usize, y: usize, z: usize) -> f64 {
        self.values[x + self.row * (y + self.size[1] as usize * z)]
    }

    fn cells(&self) -> [usize; 3] {
        [
            self.size[0] as usize - 1,
            self.size[1] as usize - 1,
            self.size[2] as usize - 1,
        ]
    }

    /// The scalar predicate: eight loads and eight comparisons, every cell.
    fn active_scalar(&self, out: &mut Vec<u32>) {
        out.clear();
        let c = self.cells();
        for z in 0..c[2] {
            for y in 0..c[1] {
                for x in 0..c[0] {
                    let mut inside = 0u32;
                    for corner in 0..8usize {
                        let v = self.at(
                            x + (corner & 1),
                            y + ((corner >> 1) & 1),
                            z + ((corner >> 2) & 1),
                        );
                        if v < 0.0 {
                            inside += 1;
                        }
                    }
                    if inside != 0 && inside != 8 {
                        out.push((x + c[0] * (y + c[1] * z)) as u32);
                    }
                }
            }
        }
    }

    /// The packed predicate: one bit per sample, 64 cells per fused test.
    fn active_bitmap(&self, bits: &mut Vec<u64>, out: &mut Vec<u32>) {
        let sx = self.size[0] as usize;
        let bit_row = sx.div_ceil(64);
        let rows = self.size[1] as usize * self.size[2] as usize;
        bits.clear();
        bits.resize(bit_row * rows, 0);
        for row in 0..rows {
            let src = self.row * row;
            let dst = bit_row * row;
            for w in 0..bit_row {
                let base = w * 64;
                let n = (sx - base).min(64);
                let mut word = 0u64;
                for k in 0..n {
                    word |= u64::from(self.values[src + base + k] < 0.0) << k;
                }
                bits[dst + w] = word;
            }
        }

        let word =
            |w: usize, y: usize, z: usize| bits[bit_row * (y + self.size[1] as usize * z) + w];
        let shifted = |w: usize, y: usize, z: usize| {
            let lo = word(w, y, z);
            let hi = if w + 1 < bit_row {
                word(w + 1, y, z)
            } else {
                0
            };
            (lo >> 1) | (hi << 63)
        };

        out.clear();
        let c = self.cells();
        for z in 0..c[2] {
            for y in 0..c[1] {
                for w in 0..bit_row {
                    let mut any = 0u64;
                    let mut all = !0u64;
                    for dz in 0..2 {
                        for dy in 0..2 {
                            let a = word(w, y + dy, z + dz);
                            let b = shifted(w, y + dy, z + dz);
                            any |= a | b;
                            all &= a & b;
                        }
                    }
                    let remaining = c[0].saturating_sub(w * 64);
                    let mask = if remaining >= 64 {
                        !0u64
                    } else {
                        (1u64 << remaining) - 1
                    };
                    let mut active = (any & !all) & mask;
                    while active != 0 {
                        let x = w * 64 + active.trailing_zeros() as usize;
                        active &= active - 1;
                        out.push((x + c[0] * (y + c[1] * z)) as u32);
                    }
                }
            }
        }
    }
}

fn median(mut v: Vec<f64>) -> f64 {
    v.sort_by(f64::total_cmp);
    v[v.len() / 2]
}

/// Median nanoseconds per sample for a closure run `REPS` times after a warm-up.
fn timed(samples: f64, mut body: impl FnMut()) -> f64 {
    body();
    let mut runs = Vec::with_capacity(REPS);
    for _ in 0..REPS {
        let t = Instant::now();
        body();
        runs.push(t.elapsed().as_secs_f64() * 1e9 / samples);
    }
    median(runs)
}

/// One whole extraction, timed, plus the hash of what it produced.
fn extract_arm(which: &str, field: &str, n: u32) -> (f64, u64) {
    let (origin, cell) = domain_of(field, n);
    let sphere = Sphere::<f64>::canonical();
    let shape = RuntimeShape3::new([n; 3]).expect("valid shape");
    let samples = f64::from(n).powi(3);
    let mut sn = SurfaceNets::<f64>::new();
    let mut dc = DualContouring::<f64>::new();
    let mut out = MeshBuffer::<f64>::new();

    let mut run = |out: &mut MeshBuffer<f64>| {
        out.reset();
        if which == "surface_nets" {
            sn.extract(&sphere, &shape, origin, cell, out)
                .expect("extraction");
        } else {
            dc.extract(&sphere, &shape, origin, cell, out)
                .expect("extraction");
        }
    };

    run(&mut out);
    let hash = mesh_hash(&out);
    let mut runs = Vec::with_capacity(REPS);
    for _ in 0..REPS {
        let t = Instant::now();
        run(&mut out);
        runs.push(t.elapsed().as_secs_f64() * 1e9 / samples);
    }
    (median(runs), hash)
}

const FIELDS: [&str; 2] = ["sphere_surface_free", "sphere"];
const EXTRACTORS: [&str; 2] = ["surface_nets", "dual_contouring"];

/// Write the parent commit's numbers, for the child commit to read.
fn bless() {
    let mut csv = String::from(
        "# P-40 baseline: the extractor arm as it stood before the bitmap prepass.\n\
         # Produced by `ISOMESH_BLESS=1 cargo bench --bench experiment_p40`, committed,\n\
         # and required by the ordinary run -- a speedup with no baseline is not a\n\
         # measurement. Regenerating it after the mechanism has landed destroys the\n\
         # comparison, which is why it is a committed fixture and not a cache.\n\
         field,samples_per_axis,extractor,ns_per_sample,mesh_hash\n",
    );
    for field in FIELDS {
        for n in SIZES {
            for which in EXTRACTORS {
                let (ns, hash) = extract_arm(which, field, n);
                csv.push_str(&format!("{field},{n},{which},{ns:.4},{hash}\n"));
                println!("blessed {field:>20} {n:>4}³ {which:>16} {ns:7.3} ns/sample");
            }
        }
    }
    std::fs::write(baseline_path(), csv).expect("write baseline");
    println!("\nbaseline → {}", baseline_path().display());
}

fn main() {
    if std::env::var_os("ISOMESH_BLESS").is_some() {
        bless();
        return;
    }

    let prereg = isomesh::experiment!("P-40");
    let baseline = read_baseline().unwrap_or_else(|| {
        panic!(
            "P-40 needs the parent commit's numbers at {}. Produce them with \
             `ISOMESH_BLESS=1 cargo bench --bench experiment_p40` on the commit \
             before the bitmap prepass landed.",
            baseline_path().display()
        )
    });

    common::experiment::run(prereg, |run| {
        for field in FIELDS {
            for n in SIZES {
                let grid = Grid::sample(field, n);
                let cells = grid.cells();
                let cell_count = (cells[0] * cells[1] * cells[2]) as f64;

                let mut scalar_list = Vec::new();
                let mut bitmap_list = Vec::new();
                let mut bits = Vec::new();
                grid.active_scalar(&mut scalar_list);
                grid.active_bitmap(&mut bits, &mut bitmap_list);
                assert_eq!(
                    scalar_list, bitmap_list,
                    "{field} {n}³: the two predicates disagree — the bitmap is not \
                     order-preserving, which would change every vertex index"
                );
                let active_fraction = scalar_list.len() as f64 / cell_count;

                let stage_scalar = timed(cell_count, || grid.active_scalar(&mut scalar_list));
                let stage_bitmap = timed(cell_count, || {
                    grid.active_bitmap(&mut bits, &mut bitmap_list)
                });

                for which in EXTRACTORS {
                    let (after_ns, after_hash) = extract_arm(which, field, n);
                    let key = (field.to_string(), n, which.to_string());
                    let (before_ns, before_hash) = baseline
                        .get(&key)
                        .copied()
                        .unwrap_or_else(|| panic!("baseline has no row for {key:?}"));

                    println!(
                        "{field:>20} {n:>4}³ {which:>16}  stage {stage_scalar:8.4} → \
                         {stage_bitmap:8.4} ns/cell ({:.2}×)   extract {before_ns:7.3} → \
                         {after_ns:7.3} ns/sample ({:.3}×)   active {:.3}%",
                        stage_scalar / stage_bitmap,
                        before_ns / after_ns,
                        active_fraction * 100.0,
                    );

                    run.record(&[
                        ("field", field.to_string()),
                        ("samples_per_axis", n.to_string()),
                        ("active_fraction", format!("{active_fraction:.6}")),
                        ("stage_ns_scalar", format!("{stage_scalar:.4}")),
                        ("stage_ns_bitmap", format!("{stage_bitmap:.4}")),
                        ("stage_ratio", format!("{:.4}", stage_scalar / stage_bitmap)),
                        ("extract_ns_scalar", format!("{before_ns:.4}")),
                        ("extract_ns_bitmap", format!("{after_ns:.4}")),
                        ("extract_ratio", format!("{:.4}", before_ns / after_ns)),
                        ("mesh_identical", (before_hash == after_hash).to_string()),
                        ("extractor", which.to_string()),
                        ("active_cells", scalar_list.len().to_string()),
                        ("cells", (cell_count as u64).to_string()),
                        ("mesh_hash", after_hash.to_string()),
                    ]);
                }
            }
        }
    });
}
