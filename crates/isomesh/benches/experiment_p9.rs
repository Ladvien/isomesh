//! **P-9 — is the `k`-way weld order-dependent?**
//!
//! Ticket: R-002. Pre-registered at R-000.
//!
//! ```bash
//! cargo bench --bench experiment_p9
//! ```
//!
//! Writes `docs/experiments/p-9.csv`.
//!
//! # Why this threatens something
//!
//! `CLAUDE.md` promises byte-identical output for identical input. Dey, Fan &
//! Wang decompose a `k`-way merge into `k − 1` pairwise merges **in the
//! intermediate complex**, so if the order those happen in is free, the answer
//! may not be a function of the input alone.
//!
//! For `Welder` the order is not free in the obvious sense — first fit walks in
//! **index order** and a vertex joins the *lowest-indexed* representative within
//! `ε`. So "within-bucket merge order" and "the order the bucket's members
//! appear in the buffer" are the same thing, and permuting one is permuting the
//! other. That is what this does.
//!
//! # What is permuted, and what deliberately is not
//!
//! Only the **slots of vertices inside one bucket of coincident vertices**, with
//! the index buffer rewritten to match. Everything else — buffer length, which
//! triangles exist, which bucket each vertex belongs to — is untouched. A
//! wholesale shuffle of the vertex array would change the output trivially and
//! would measure nothing about the weld.
//!
//! # The comparison is byte-identity, per the registration
//!
//! Not "the same up to a tolerance". Coincident vertices differ by about an ulp
//! (M-32), so if the surviving representative changes, the output changes by an
//! ulp — and an ulp is exactly what a byte-identity guarantee is about.

mod common;

use std::collections::BTreeMap;

use isomesh::chunk::{ChunkId, ChunkLayout};
use isomesh::extractor::Extractor;
use isomesh::weld::{Welder, epsilon_for};
use isomesh::{MeshBuffer, Sdf};

/// Cells per chunk. Same fixture as P-8, on purpose: the two experiments are
/// about the same buckets and a different one would make their numbers
/// incomparable. See P-8 for why 18 rather than 8 (M-274).
const CELLS: u32 = 18;
/// See P-8 — `4/35` because a seam is bit-exact only at a power of two (M-32).
const CELL_SIZE: f64 = 4.0 / 35.0;

/// Block origin, placed so the 2×2×2 block is **centred on the field**.
///
/// The eight-chunk corner then sits at the origin and the four-chunk edges cross
/// the surface, which is where a bucket can hold more than two vertices at all.
const ORIGIN: f64 = -(2.0 * CELLS as f64) * CELL_SIZE / 2.0;
/// Permutations per field.
const PERMUTATIONS: u32 = 8;

/// A seeded xorshift64*, written out rather than depended on.
///
/// The crate has no RNG and does not want one for this. What matters is that the
/// sequence is **reproducible from the seed** — an experiment whose shuffle
/// cannot be replayed cannot be re-checked, and that is worse than not shuffling.
struct Rng(u64);

impl Rng {
    const fn new(seed: u64) -> Self {
        // Any non-zero state; xorshift is degenerate at zero.
        Self(seed | 1)
    }

    fn next(&mut self) -> u64 {
        let mut x = self.0;
        x ^= x >> 12;
        x ^= x << 25;
        x ^= x >> 27;
        self.0 = x;
        x.wrapping_mul(0x2545_F491_4F6C_DD1D)
    }

    /// Fisher–Yates, so every ordering is reachable.
    fn shuffle<T>(&mut self, slice: &mut [T]) {
        for i in (1..slice.len()).rev() {
            let j = (self.next() % (i as u64 + 1)) as usize;
            slice.swap(i, j);
        }
    }
}

/// FNV-1a over a mesh's bytes.
///
/// Only distinctness is needed, not cryptographic strength — this counts how
/// many different outputs the permutations produced, and a collision would
/// under-report, which is the safe direction for a hypothesis predicting that
/// the count is greater than one.
fn digest(mesh: &MeshBuffer<f64>) -> u64 {
    let mut h: u64 = 0xcbf2_9ce4_8422_2325;
    let mut eat = |b: u8| {
        h ^= u64::from(b);
        h = h.wrapping_mul(0x1000_0000_01b3);
    };
    for p in &mesh.positions {
        for c in p {
            for b in c.to_le_bytes() {
                eat(b);
            }
        }
    }
    for n in &mesh.normals {
        for c in n {
            for b in c.to_le_bytes() {
                eat(b);
            }
        }
    }
    for i in &mesh.indices {
        for b in i.to_le_bytes() {
            eat(b);
        }
    }
    h
}

/// Eight chunks in a 2×2×2 block, meshed independently and appended.
fn eight_chunks<E: Extractor<f64>>(
    field: &impl Sdf<Scalar = f64>,
    layout: &ChunkLayout<f64>,
    extractor: &mut E,
) -> MeshBuffer<f64> {
    let shape = layout.sample_shape().expect("valid shape");
    let mut joined = MeshBuffer::<f64>::new();
    for z in 0..2 {
        for y in 0..2 {
            for x in 0..2 {
                let id = ChunkId::new([x, y, z]);
                let mut piece = MeshBuffer::<f64>::new();
                extractor
                    .extract_into(
                        field,
                        &shape,
                        layout.sample_origin(id),
                        layout.cell_size(),
                        &mut piece,
                    )
                    .expect("extraction");
                joined.append(&piece).expect("the meshes fit u32");
            }
        }
    }
    joined
}

/// Move each bucket's members among their own slots, rewriting indices to match.
fn permute_within_buckets(
    mesh: &MeshBuffer<f64>,
    buckets: &[Vec<u32>],
    rng: &mut Rng,
) -> MeshBuffer<f64> {
    let n = mesh.positions.len();
    // `to[old] = new`. Identity outside the buckets.
    let mut to: Vec<u32> = (0..n as u32).collect();
    for members in buckets {
        if members.len() < 2 {
            continue;
        }
        let mut shuffled = members.clone();
        rng.shuffle(&mut shuffled);
        for (slot, &member) in members.iter().zip(&shuffled) {
            to[member as usize] = *slot;
        }
    }

    let mut out = MeshBuffer::<f64>::new();
    out.positions = vec![[0.0; 3]; n];
    out.normals = vec![[0.0; 3]; n];
    for (old, &slot) in to.iter().enumerate() {
        let new = slot as usize;
        out.positions[new] = mesh.positions[old];
        out.normals[new] = mesh.normals[old];
    }
    out.indices = mesh.indices.iter().map(|&i| to[i as usize]).collect();
    out
}

fn main() {
    if !std::env::args().any(|arg| arg == "--bench") {
        return;
    }

    let prereg = isomesh::experiment!("P-9");
    common::experiment::run(prereg, |run| {
        let layout = ChunkLayout::<f64>::new(CELLS, CELL_SIZE, [ORIGIN; 3]).expect("valid layout");

        println!(
            "{:<16} {:<28} {:>10} {:>16} {:>9}",
            "field", "extractor", "distinct", "vertex spread", "k≥3"
        );

        isomesh::for_each_reference_field!(f64, |name, field| {
            // Inline blocks, so no `return` in either (M-253).
            isomesh::for_each_extractor!(f64, |ename, extractor| {
                let mesh = eight_chunks(&field, &layout, &mut extractor);
                if !mesh.indices.is_empty() {
                    // Buckets, from one reference weld.
                    let mut probe = mesh.clone();
                    let mut welder = Welder::<f64>::new();
                    welder
                        .weld(&mut probe, epsilon_for(CELL_SIZE))
                        .expect("valid epsilon");
                    let mut groups: BTreeMap<u32, Vec<u32>> = BTreeMap::new();
                    for (v, &to) in welder.remap().iter().enumerate() {
                        groups.entry(to).or_default().push(v as u32);
                    }
                    let buckets: Vec<Vec<u32>> =
                        groups.into_values().filter(|g| g.len() > 1).collect();
                    let big = buckets.iter().filter(|g| g.len() >= 3).count();

                    let mut digests = std::collections::BTreeSet::new();
                    let mut lo = usize::MAX;
                    let mut hi = 0usize;
                    // Seeded from the loop index alone, so the run is
                    // reproducible and every field sees the same permutations.
                    for seed in 0..PERMUTATIONS {
                        let mut rng = Rng::new(u64::from(seed) * 0x9E37_79B9_7F4A_7C15 + 1);
                        let mut shuffled = permute_within_buckets(&mesh, &buckets, &mut rng);
                        let mut w = Welder::<f64>::new();
                        w.weld(&mut shuffled, epsilon_for(CELL_SIZE))
                            .expect("valid epsilon");
                        digests.insert(digest(&shuffled));
                        lo = lo.min(shuffled.positions.len());
                        hi = hi.max(shuffled.positions.len());
                    }

                    println!(
                        "{name:<16} {ename:<28} {:>10} {:>16} {big:>9}",
                        digests.len(),
                        hi - lo
                    );
                    run.record(&[
                        ("field", name.to_string()),
                        ("extractor", ename.to_string()),
                        ("distinct_outputs", digests.len().to_string()),
                        ("vertex_count_spread", (hi - lo).to_string()),
                        ("buckets_of_three_or_more", big.to_string()),
                    ]);
                }
            });
        });
    });
}
