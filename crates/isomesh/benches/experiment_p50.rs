//! **P-50 — the count the ratio was sampling.**
//!
//! Ticket: R-039a. Pre-registered in the commit before this one.
//!
//! ```bash
//! cargo bench --bench experiment_p50
//! ```
//!
//! Writes `docs/experiments/p-50.csv`.
//!
//! # Why this exists
//!
//! P-40's C2 was a wall-clock ratio registered as a threshold. It failed its own
//! committed artefact and then failed three quiet-machine re-runs (M-348), and
//! ✗24 had already earned the rule: *a wall-clock ratio is not a gate; gate the
//! count the ratio samples.*
//!
//! The bitmap prepass does not make an eight-corner gather faster. **It removes
//! gathers.** How many run is an integer, and an integer does not care about a
//! governor, a load average or a machine. So every clause here is an equality
//! over counters, and the two timing columns are reported because they are
//! interesting and gate nothing.
//!
//! # The predicates are re-implemented here, deliberately
//!
//! `DualMesher`'s buffers are private, so both arms are rebuilt in this file over
//! a grid laid out exactly as the crate lays out `values` — row stride
//! `size[0] | 1`, A-024. That makes the two arms comparable to each other, which
//! is what the clauses are about. The check that they also match the *crate* is
//! `p-40.csv`'s `mesh_identical`, which is deterministic and already 12 of 12.

mod common;

use std::time::Instant;

use isomesh::Sdf;
use isomesh::fields::{CappedGyroid, FbmTerrain, ReferenceField, Sphere, capped_gyroid};

const REPS: usize = 5;

/// Which field a row is about, and how to sample it.
#[derive(Clone, Copy)]
enum Field {
    Sphere,
    Gyroid,
    Terrain,
}

impl Field {
    const fn name(self) -> &'static str {
        match self {
            Self::Sphere => "sphere",
            Self::Gyroid => "gyroid",
            Self::Terrain => "fbm_terrain",
        }
    }
}

/// The sampled grid, in `DualMesher`'s own layout.
struct Grid {
    values: Vec<f64>,
    row: usize,
    size: [u32; 3],
    /// Comparisons the bitmap build performs: one per real sample.
    comparisons: u64,
}

impl Grid {
    fn sample(field: Field, n: u32) -> Self {
        let sphere = Sphere::<f64>::canonical();
        let gyroid: CappedGyroid<f64> = capped_gyroid();
        let terrain = FbmTerrain::<f64>::canonical();
        let (lo, hi) = match field {
            Field::Sphere => sphere.domain(),
            Field::Gyroid => gyroid.domain(),
            Field::Terrain => terrain.domain(),
        };
        let cell = (hi[0] - lo[0]) / f64::from(n - 1);

        let row = n as usize | 1;
        let mut values = Vec::with_capacity(row * n as usize * n as usize);
        for z in 0..n {
            for y in 0..n {
                for x in 0..n {
                    let p = [
                        lo[0] + cell * f64::from(x),
                        lo[1] + cell * f64::from(y),
                        lo[2] + cell * f64::from(z),
                    ];
                    values.push(match field {
                        Field::Sphere => sphere.sample(p),
                        Field::Gyroid => gyroid.sample(p),
                        Field::Terrain => terrain.sample(p),
                    });
                }
                values.resize(values.len() + (row - n as usize), 0.0);
            }
        }
        Self {
            values,
            row,
            size: [n; 3],
            comparisons: u64::from(n) * u64::from(n) * u64::from(n),
        }
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

    /// The scalar predicate. Returns the active list and **the number of
    /// eight-corner gathers it performed**, which is the whole point.
    fn active_scalar(&self, out: &mut Vec<u32>) -> u64 {
        out.clear();
        let c = self.cells();
        let mut gathers = 0u64;
        for z in 0..c[2] {
            for y in 0..c[1] {
                for x in 0..c[0] {
                    // One gather per cell, unconditionally: that is the thing
                    // the bitmap removes, and it is counted rather than timed.
                    gathers += 1;
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
        gathers
    }

    /// The packed predicate. Returns the active list, the gathers performed —
    /// one per *active* cell — and the fused word groups evaluated.
    fn active_bitmap(&self, bits: &mut Vec<u64>, out: &mut Vec<u32>) -> (u64, u64) {
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
        // Words that carry a **cell**, not words that carry a sample. The
        // difference is one whole word per cell row at every `64k+1` grid, and
        // it is the defect E-307 found and M-348 fixed in `dual.rs`.
        let cell_words = c[0].div_ceil(64);
        let mut gathers = 0u64;
        let mut word_groups = 0u64;
        for z in 0..c[2] {
            for y in 0..c[1] {
                for w in 0..cell_words {
                    word_groups += 1;
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
                        // A gather runs only for a cell the mask named.
                        gathers += 1;
                        out.push((x + c[0] * (y + c[1] * z)) as u32);
                    }
                }
            }
        }
        (gathers, word_groups)
    }
}

fn median(mut v: Vec<f64>) -> f64 {
    v.sort_by(f64::total_cmp);
    v[v.len() / 2]
}

fn main() {
    let prereg = isomesh::experiment!("P-50");

    // 64k and 64k+1 on both sides of a word boundary, because the word-group
    // prediction is exactly what differs between them.
    const GRIDS: [u32; 5] = [33, 64, 65, 128, 129];

    common::experiment::run(prereg, |run| {
        for field in [Field::Sphere, Field::Gyroid, Field::Terrain] {
            for n in GRIDS {
                let grid = Grid::sample(field, n);
                let c = grid.cells();
                let cells = (c[0] * c[1] * c[2]) as u64;

                let mut scalar_list = Vec::new();
                let mut bitmap_list = Vec::new();
                let mut bits = Vec::new();

                let gathers_scalar = grid.active_scalar(&mut scalar_list);
                let (gathers_bitmap, word_groups) = grid.active_bitmap(&mut bits, &mut bitmap_list);

                let active = scalar_list.len() as u64;
                let same_ordered = scalar_list == bitmap_list;
                let predicted = (c[0].div_ceil(64) * c[1] * c[2]) as u64;

                let ns_scalar = {
                    grid.active_scalar(&mut scalar_list);
                    let mut r = Vec::with_capacity(REPS);
                    for _ in 0..REPS {
                        let t = Instant::now();
                        let g = grid.active_scalar(&mut scalar_list);
                        std::hint::black_box(g);
                        r.push(t.elapsed().as_secs_f64() * 1e9 / cells as f64);
                    }
                    median(r)
                };
                let ns_bitmap = {
                    grid.active_bitmap(&mut bits, &mut bitmap_list);
                    let mut r = Vec::with_capacity(REPS);
                    for _ in 0..REPS {
                        let t = Instant::now();
                        let g = grid.active_bitmap(&mut bits, &mut bitmap_list);
                        std::hint::black_box(g);
                        r.push(t.elapsed().as_secs_f64() * 1e9 / cells as f64);
                    }
                    median(r)
                };

                println!(
                    "{:>12} {n:>4}³  gathers {gathers_scalar:>9} → {gathers_bitmap:>8}  \
                     (cells {cells}, active {active})  words {word_groups} = {predicted}?  \
                     ordered {same_ordered}",
                    field.name()
                );

                run.record(&[
                    ("field", field.name().to_string()),
                    ("samples_per_axis", n.to_string()),
                    ("cells", cells.to_string()),
                    ("active_cells", active.to_string()),
                    ("gathers_scalar", gathers_scalar.to_string()),
                    ("gathers_bitmap", gathers_bitmap.to_string()),
                    ("gathers_equal_cells", (gathers_scalar == cells).to_string()),
                    (
                        "gathers_equal_active",
                        (gathers_bitmap == active).to_string(),
                    ),
                    ("bitmap_comparisons", grid.comparisons.to_string()),
                    ("bitmap_word_groups", word_groups.to_string()),
                    ("word_groups_predicted", predicted.to_string()),
                    ("same_ordered_list", same_ordered.to_string()),
                    ("ns_per_cell_scalar", format!("{ns_scalar:.4}")),
                    ("ns_per_cell_bitmap", format!("{ns_bitmap:.4}")),
                    ("word_groups_match", (word_groups == predicted).to_string()),
                    (
                        "gathers_removed_fraction",
                        format!("{:.6}", 1.0 - active as f64 / cells as f64),
                    ),
                ]);
            }
        }
    });
}
