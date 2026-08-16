//! **P-11 — is a seam crack arithmetic, or algorithm?**
//!
//! Ticket: R-004. Pre-registered at R-000.
//!
//! ```bash
//! cargo bench --bench experiment_p11
//! ```
//!
//! Writes `docs/experiments/p-11.csv`.
//!
//! # The three measurements this is built on
//!
//! **M-32** — two adjacent chunks compute their shared plane as `(o + h·cn) + h·n`
//! and as `o + h·(c+1)n`, equal by algebra and not by IEEE, and 22% of random
//! `(origin, h, cells, chunk)` combinations disagree by an ulp. **M-49** — the
//! same identity failing in `cell_of`. **M-73** — a transition cell that offsets
//! from a face origin puts a crossing at `y = -1.11e-16` where the coarse mesh
//! has one at exactly `0`, *"and no weld can close it"*.
//!
//! Three subsystems, one cause: a world position reached by adding a step to a
//! block origin rather than by one multiply from a global integer index. R-004
//! asks what that costs, against what the *algorithm* costs — a coarse block and
//! a fine block genuinely disagree about where the surface is, and no arithmetic
//! fixes that.
//!
//! # The two arms, and why neither is a second code path
//!
//! Both arms call the same functions. They differ only in **what they are given**:
//!
//! | | block origin | sample position |
//! |---|---|---|
//! | `canonical` | the **grid** origin, plus the block's integer base | `o + h·(base + local)` |
//! | `offset` | `world_of_sample(base)`, the block's own corner | `(o + h·base) + h·local` |
//!
//! `TransitionCell::sample` takes exactly this pair — `origin` and an integer
//! `base` — so passing the face's own world origin with a zero base reproduces
//! the pre-M-73 code character for character, using the shipped function. There
//! is no second implementation of anything here, which matters: two copies of a
//! formula agreeing proves only that they were written on the same day.
//!
//! `MarchingCubes::extract` has no integer base, so its canonical arm is
//! obtained by rooting the extraction at the **grid** origin and keeping the
//! triangles of the block's own cells. Marching Cubes is cell-local — a cell's
//! triangles depend on its eight corners and nothing else — so the kept set is
//! the set the block would have produced.
//!
//! That is an argument, and [`clip_agrees_with_the_block`] is the check.
//! At a power-of-two spacing the two arms compute **identical** positions, so
//! any difference there is the clip and nothing else: the clipped extraction and
//! the block's own must have the same triangle count and the same multiset of
//! referenced positions, bit for bit. It runs on every power-of-two row before
//! that row is measured.
//!
//! # Why the fixture is offset from the grid origin
//!
//! A block based at `[0, 0, 0]` computes `o + h·local` under **both** arms, the
//! two coincide, and the experiment measures nothing. That is precisely the trap
//! M-32's first fixture fell into — `h = 4/33` *looks* irregular and lands in the
//! 78% of cases that agree exactly — and Part 5 of `FINDINGS.md` carries the rule
//! it earned. So every block here starts [`OFFSET_CELLS`] cells from the grid
//! origin, and `seam_plane_delta` records whether the two arms actually reached
//! different bits. A row with `seam_plane_delta = 0` is a null fixture and says
//! so in the CSV rather than passing quietly.

mod common;

use isomesh::chunk::ChunkLayout;
use isomesh::fields::{ReferenceField, Sphere, Torus};
use isomesh::marching_cubes::MarchingCubes;
use isomesh::transvoxel::cell::TransitionCell;
use isomesh::validate::{ValidateConfig, validate_features};
use isomesh::weld::{Welder, epsilon_for};
use isomesh::{MeshBuffer, RuntimeShape3, Sdf};

/// Fine cells between the grid origin and the fixture's own minimum corner.
///
/// Non-zero on every axis on purpose — see the module docs. Twelve halves twice
/// (12, 6, 3), which the LOD pairs need, and is not a power of two, so `h·12` is
/// a generic product rather than an exact rescale.
const OFFSET_CELLS: i64 = 12;

/// Level-0 cells across the 4-unit domain. `h = 4/n`.
///
/// Divisible by eight, so both LOD pairs tile: the seam index must be even at
/// the fine level of each pair. Two are powers of two (`32 → 0.125`,
/// `64 → 0.0625`) and three are not (`40 → 0.1`, `48 → 1/12`, `56 → 1/14`),
/// which is the comparison M-32's "recommend power-of-two cell sizes" turns on.
const RESOLUTIONS: [u32; 5] = [32, 40, 48, 56, 64];

/// Which level of the pair the *fine* block sits at; the coarse block is one up.
const FINE_LEVELS: [u32; 2] = [0, 1];

/// Which of the two expressions a block's sample positions are built from.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Arithmetic {
    /// `o + h·(base + local)` — one multiply from a global integer index.
    Canonical,
    /// `(o + h·base) + h·local` — a step added to the block's own corner.
    Offset,
}

impl Arithmetic {
    const ALL: [Self; 2] = [Self::Canonical, Self::Offset];

    const fn name(self) -> &'static str {
        match self {
            Self::Canonical => "canonical",
            Self::Offset => "offset",
        }
    }
}

/// A layout whose only job is [`ChunkLayout::world_of_sample`].
///
/// The crate's own expression rather than a copy of it: this harness exists to
/// measure that expression, so re-typing it here would measure the retyping.
/// `cells` names a chunk size the fixture does not use — the blocks are
/// rectangular and each states its own extent — so it is 1.
fn layout(step: f64, origin: [f64; 3]) -> ChunkLayout<f64> {
    ChunkLayout::<f64>::new(1, step, origin).expect("a positive step is a valid layout")
}

/// Keep the triangles whose cells lie at or above `lo` on every axis.
///
/// Classified on the triangle's **coordinate range**, not on a centroid: the cut
/// is a grid plane, so no cell straddles it and every triangle lies wholly on one
/// side. A cell below spans `[lo − h, lo]` and has `min < lo`; a cell at the cut
/// spans `[lo, lo + h]` and has `max > lo`. Both comparisons are exact on
/// coordinates the extractor produced, where a centroid would be an arithmetic of
/// its own and could round across the boundary.
///
/// `min == max == lo` is the one case the range cannot place: a triangle lying
/// **in** the plane belongs to the cell on either side of it. Those are counted
/// into `ambiguous` and the caller asserts the count is zero rather than picking
/// a side.
fn clip_to_block(mesh: &MeshBuffer<f64>, lo: [f64; 3], ambiguous: &mut usize) -> MeshBuffer<f64> {
    // Exact: this is a question about which cell emitted a triangle, and a
    // tolerance would answer a different one.
    #![allow(clippy::float_cmp)]
    let mut out = MeshBuffer::<f64>::new();
    let mut remap = alloc_remap(mesh.positions.len());

    for tri in mesh.indices.chunks_exact(3) {
        let p = [
            mesh.positions[tri[0] as usize],
            mesh.positions[tri[1] as usize],
            mesh.positions[tri[2] as usize],
        ];
        let mut keep = true;
        let mut in_plane = false;
        for (axis, cut) in lo.iter().enumerate() {
            let least = p[0][axis].min(p[1][axis]).min(p[2][axis]);
            let most = p[0][axis].max(p[1][axis]).max(p[2][axis]);
            if least == *cut && most == *cut {
                in_plane = true;
                keep = false;
            } else if most <= *cut {
                keep = false;
            }
        }
        if in_plane {
            *ambiguous += 1;
        }
        if !keep {
            continue;
        }
        for &v in tri {
            if remap[v as usize] == u32::MAX {
                remap[v as usize] = out.positions.len() as u32;
                out.positions.push(mesh.positions[v as usize]);
                out.normals.push(mesh.normals[v as usize]);
            }
            out.indices.push(remap[v as usize]);
        }
    }
    out
}

/// `u32::MAX`-filled remap of `n` entries.
fn alloc_remap(n: usize) -> Vec<u32> {
    vec![u32::MAX; n]
}

/// Mesh one rectangular block of cells under one arithmetic.
///
/// `base` is the block's minimum corner as a global sample index **at this
/// block's own spacing**; `cells` is what it owns.
fn mesh_block<S: Sdf<Scalar = f64>>(
    field: &S,
    origin: [f64; 3],
    step: f64,
    base: [i64; 3],
    cells: [u32; 3],
    arm: Arithmetic,
    ambiguous: &mut usize,
) -> MeshBuffer<f64> {
    let mut mesh = MeshBuffer::<f64>::new();
    match arm {
        Arithmetic::Offset => {
            let shape = RuntimeShape3::new([cells[0] + 1, cells[1] + 1, cells[2] + 1])
                .expect("the fixture fits u32");
            let lo = layout(step, origin).world_of_sample(base);
            MarchingCubes::<f64>::new()
                .extract(field, &shape, lo, step, &mut mesh)
                .expect("extraction");
        }
        Arithmetic::Canonical => {
            // Rooted at the grid origin, so a local index *is* a global one and
            // `origin + step·local` is the canonical expression by construction.
            let shape = RuntimeShape3::new([
                base[0] as u32 + cells[0] + 1,
                base[1] as u32 + cells[1] + 1,
                base[2] as u32 + cells[2] + 1,
            ])
            .expect("the fixture fits u32");
            let mut full = MeshBuffer::<f64>::new();
            MarchingCubes::<f64>::new()
                .extract(field, &shape, origin, step, &mut full)
                .expect("extraction");
            let lo = layout(step, origin).world_of_sample(base);
            mesh = clip_to_block(&full, lo, ambiguous);
        }
    }
    mesh
}

/// The two-resolution fixture: one grid, two blocks, one seam.
///
/// A struct rather than eight parameters because every part of the measurement
/// needs the same eight, and a fixture whose pieces can be passed in different
/// combinations is a fixture that can be passed in the wrong one.
struct Fixture {
    /// World position of global sample `[0, 0, 0]`, at every level.
    origin: [f64; 3],
    /// Spacing of the fine block; the coarse block is at twice this.
    fine_step: f64,
    /// The fine block's minimum corner, in fine samples, and its cells.
    fine_base: [i64; 3],
    fine_cells: [u32; 3],
    /// The coarse block's minimum corner, in **coarse** samples, and its cells.
    coarse_base: [i64; 3],
    coarse_cells: [u32; 3],
    /// The seam plane's global sample index at the fine spacing.
    seam: i64,
}

impl Fixture {
    /// Spacing of the coarse block. Doubling is exact in IEEE (M-70).
    fn coarse_step(&self) -> f64 {
        self.fine_step + self.fine_step
    }

    /// Every transition cell on the seam face, under one arithmetic.
    ///
    /// One per coarse cell face, spanning two fine cells on each in-plane axis.
    ///
    /// The width is zero, which is Lengyel §4.3's own position: it *"still
    /// produce\[s\] results that seamlessly stitch multiresolution meshes
    /// together"* and costs the shading, not the stitch (M-74). A non-zero width
    /// additionally needs Equation 4.2's inset of the coarse block's boundary
    /// cells, which would put a second variable in a two-arm experiment.
    fn transition_patches<S: Sdf<Scalar = f64>>(
        &self,
        field: &S,
        arm: Arithmetic,
    ) -> (MeshBuffer<f64>, usize) {
        let grid = layout(self.fine_step, self.origin);
        let mut out = MeshBuffer::<f64>::new();
        let mut emitted = 0usize;
        for jw in 0..i64::from(self.coarse_cells[2]) {
            for jv in 0..i64::from(self.coarse_cells[1]) {
                let index = [
                    self.seam,
                    self.fine_base[1] + 2 * jv,
                    self.fine_base[2] + 2 * jw,
                ];
                let cell = match arm {
                    Arithmetic::Canonical => {
                        TransitionCell::sample(field, self.origin, self.fine_step, index, 1, 2, 0.0)
                    }
                    // The pre-M-73 shape: the face's own world origin, and local
                    // offsets added to it.
                    Arithmetic::Offset => TransitionCell::sample(
                        field,
                        grid.world_of_sample(index),
                        self.fine_step,
                        [0, 0, 0],
                        1,
                        2,
                        0.0,
                    ),
                };
                let mut patch = MeshBuffer::<f64>::new();
                cell.emit(field, 0, &mut patch);
                if patch.triangle_count() > 0 {
                    emitted += 1;
                    out.append(&patch).expect("the meshes fit u32");
                }
            }
        }
        (out, emitted)
    }

    /// The seam plane's world `x`, as each side of it computes the value.
    ///
    /// The coarse block and the transition cells name the plane by its global
    /// index; the fine block reaches it by adding `fine_cells[0]` steps to its
    /// own corner. Under `canonical` it does that by index too and the two
    /// collapse to one value — which is the whole of the arm distinction, in one
    /// number.
    fn seam_plane(&self, arm: Arithmetic) -> [f64; 2] {
        let grid = layout(self.fine_step, self.origin);
        let by_index = grid.world_of_sample([self.seam, 0, 0])[0];
        let by_block = match arm {
            Arithmetic::Canonical => by_index,
            Arithmetic::Offset => {
                grid.world_of_sample(self.fine_base)[0]
                    + self.fine_step * f64::from(self.fine_cells[0])
            }
        };
        [by_block, by_index]
    }
}

/// Everything one joined mesh says about its seam.
struct SeamReport {
    /// Boundary edges lying wholly in the seam plane.
    cracks: usize,
    /// The widest hole, in cells — M-106's metric, in a form that does not need
    /// a ray caster.
    ///
    /// A crack is bounded by two **lips**, and its width is how far one lip is
    /// from the other. So: for every endpoint of a seam-plane boundary edge, the
    /// distance to the nearest *other* such endpoint it is not joined to along
    /// the boundary — its own lip's next vertex is a cell away, the opposite lip
    /// is the crack's width away, and the minimum takes whichever is nearer. The
    /// maximum over all of them is the widest place. Zero when there are no
    /// seam-plane boundary edges at all, which is a seam that closed.
    discontinuity: f64,
}

/// How two coincident vertices are made into one.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Merge {
    /// [`epsilon_for`]`(step)` — the crate's one weld policy (T-009).
    Weld,
    /// Bit-identity. What "shared by construction" means, and the only merge a
    /// consumer that never welds — a collider (M-69) — gets for free.
    Exact,
}

/// Collapse bit-identical positions, returning the remap.
///
/// `-0.0` is folded onto `0.0` first. The two are equal and hash differently,
/// and a seam vertex landing on the origin plane is exactly where that would
/// bite — so the fold is a correctness fix, not a tidy-up.
fn merge_exact(mesh: &mut MeshBuffer<f64>) -> Vec<u32> {
    let mut seen: std::collections::BTreeMap<[u64; 3], u32> = std::collections::BTreeMap::new();
    let mut remap = alloc_remap(mesh.positions.len());
    let mut positions = Vec::new();
    let mut normals = Vec::new();
    for (v, p) in mesh.positions.iter().enumerate() {
        let key = [
            (p[0] + 0.0).to_bits(),
            (p[1] + 0.0).to_bits(),
            (p[2] + 0.0).to_bits(),
        ];
        remap[v] = *seen.entry(key).or_insert_with(|| {
            positions.push(*p);
            normals.push(mesh.normals[v]);
            positions.len() as u32 - 1
        });
    }
    mesh.positions = positions;
    mesh.normals = normals;
    let mut indices = Vec::with_capacity(mesh.indices.len());
    for tri in mesh.indices.chunks_exact(3) {
        let t = [
            remap[tri[0] as usize],
            remap[tri[1] as usize],
            remap[tri[2] as usize],
        ];
        // Same rule as `Welder::weld`: a triangle two of whose corners became
        // one vertex has no area left.
        if t[0] != t[1] && t[1] != t[2] && t[2] != t[0] {
            indices.extend_from_slice(&t);
        }
    }
    mesh.indices = indices;
    remap
}

/// Join, merge, and measure the seam.
fn seam_report(
    parts: &[&MeshBuffer<f64>],
    seam_x: [f64; 2],
    fine_step: f64,
    policy: Merge,
) -> (SeamReport, MeshBuffer<f64>) {
    // Exact: a vertex is in the seam plane because an interpolation between two
    // samples at that plane put it there identically, and a tolerance would
    // sweep in vertices merely near it.
    #![allow(clippy::float_cmp)]
    let mut joined = MeshBuffer::<f64>::new();
    for mesh in parts {
        joined.append(mesh).expect("the meshes fit u32");
    }
    match policy {
        Merge::Weld => {
            Welder::<f64>::new()
                .weld(&mut joined, epsilon_for(fine_step))
                .expect("a positive epsilon");
        }
        Merge::Exact => {
            merge_exact(&mut joined);
        }
    }

    let in_seam = |p: [f64; 3]| p[0] == seam_x[0] || p[0] == seam_x[1];
    let cfg = ValidateConfig::from_cell_size(fine_step).expect("a positive cell size");
    let (_report, features) = validate_features(&joined.positions, &joined.indices, &cfg);
    let cracks: Vec<[u32; 2]> = features
        .boundary_edges
        .iter()
        .copied()
        .filter(|[a, b]| {
            in_seam(joined.positions[*a as usize]) && in_seam(joined.positions[*b as usize])
        })
        .collect();

    let mut lips: Vec<u32> = cracks.iter().flat_map(|e| e.iter().copied()).collect();
    lips.sort_unstable();
    lips.dedup();
    let joined_along: std::collections::BTreeSet<[u32; 2]> = cracks.iter().copied().collect();

    let mut discontinuity = 0.0f64;
    for &v in &lips {
        let p = joined.positions[v as usize];
        let mut nearest = f64::INFINITY;
        for &u in &lips {
            let pair = if u < v { [u, v] } else { [v, u] };
            if u == v || joined_along.contains(&pair) {
                continue;
            }
            let q = joined.positions[u as usize];
            let d = (p[0] - q[0]).hypot(p[1] - q[1]).hypot(p[2] - q[2]);
            nearest = nearest.min(d);
        }
        if nearest.is_finite() {
            discontinuity = discontinuity.max(nearest / fine_step);
        }
    }

    (
        SeamReport {
            cracks: cracks.len(),
            discontinuity,
        },
        joined,
    )
}

/// How closely the two sides' seam vertices coincide, **before** any weld.
///
/// This is the arithmetic itself, with the weld's tolerance taken out of the
/// picture: `exact` counts partners that are bit-identical, and `gap` is the
/// worst disagreement among partners the weld would still merge.
struct Coincidence {
    pairs: usize,
    exact: usize,
    gap: f64,
}

fn coincidence(
    fine: &MeshBuffer<f64>,
    others: &[&MeshBuffer<f64>],
    seam_x: [f64; 2],
    fine_step: f64,
) -> Coincidence {
    #![allow(clippy::float_cmp)]
    let in_seam = |p: [f64; 3]| p[0] == seam_x[0] || p[0] == seam_x[1];
    let mut far: Vec<[f64; 3]> = Vec::new();
    for mesh in others {
        far.extend(mesh.positions.iter().copied().filter(|p| in_seam(*p)));
    }

    let epsilon = epsilon_for(fine_step);
    let mut out = Coincidence {
        pairs: 0,
        exact: 0,
        gap: 0.0,
    };
    for p in fine.positions.iter().filter(|p| in_seam(**p)) {
        let mut best = f64::INFINITY;
        for q in &far {
            let d = (p[0] - q[0]).hypot(p[1] - q[1]).hypot(p[2] - q[2]);
            best = best.min(d);
        }
        if best <= epsilon {
            out.pairs += 1;
            if best == 0.0 {
                out.exact += 1;
            } else {
                out.gap = out.gap.max(best);
            }
        }
    }
    out
}

/// The fixture for one field at one resolution and one LOD pair.
///
/// The fine block owns the negative half of the domain on `x` and all of it on
/// `y` and `z`; the coarse block owns the positive half at twice the spacing.
/// The surface crosses the seam, which is what makes the seam measurable at all.
fn fixture<F>(field: &F, cells_per_axis: u32, level: u32) -> Fixture
where
    F: ReferenceField + Sdf<Scalar = f64>,
{
    let (lo, hi) = field.domain();
    let h0 = (hi[0] - lo[0]) / f64::from(cells_per_axis);
    // Level `k` doubles the spacing `k` times, exactly as `ChunkLayout::at_lod`
    // does, and doubling is exact in IEEE (M-70).
    let mut fine_step = h0;
    for _ in 0..level {
        fine_step += fine_step;
    }
    let scale = 1i64 << level;

    // The grid origin sits `OFFSET_CELLS` level-0 cells below the domain, so
    // every block has a non-zero base and the two arms can differ at all.
    let origin = [lo[0] - h0 * OFFSET_CELLS as f64; 3];
    let base = OFFSET_CELLS / scale;
    let across = cells_per_axis / scale as u32;
    let fine_cells = [across / 2, across, across];
    let seam = base + i64::from(fine_cells[0]);

    Fixture {
        origin,
        fine_step,
        fine_base: [base; 3],
        fine_cells,
        coarse_base: [seam / 2, base / 2, base / 2],
        coarse_cells: [fine_cells[0] / 2, across / 2, across / 2],
        seam,
    }
}

/// Every position a triangle names, as bits, sorted — a comparable multiset.
fn referenced_positions(mesh: &MeshBuffer<f64>) -> Vec<[u64; 3]> {
    let mut out: Vec<[u64; 3]> = mesh
        .indices
        .iter()
        .map(|&i| {
            let p = mesh.positions[i as usize];
            [
                (p[0] + 0.0).to_bits(),
                (p[1] + 0.0).to_bits(),
                (p[2] + 0.0).to_bits(),
            ]
        })
        .collect();
    out.sort_unstable();
    out
}

/// Is `step` a power of two? Then `step · k` is exact for every integer `k` in
/// range and the two arms cannot differ (M-32).
fn is_power_of_two(step: f64) -> bool {
    step.is_finite() && step > 0.0 && step.to_bits() & ((1u64 << 52) - 1) == 0
}

/// **The clip's own acceptance.** See the module docs.
///
/// Only meaningful at a power-of-two spacing, where the two arms are the same
/// arithmetic and a difference can only be the clip. Silent elsewhere: at any
/// other spacing the two meshes are *supposed* to differ, which is the
/// experiment.
///
/// # Panics
///
/// If the clipped extraction is not the block's own mesh.
fn clip_agrees_with_the_block<S: Sdf<Scalar = f64>>(
    field: &S,
    origin: [f64; 3],
    step: f64,
    base: [i64; 3],
    cells: [u32; 3],
) {
    if !is_power_of_two(step) {
        return;
    }
    let mut ambiguous = 0usize;
    let clipped = mesh_block(
        field,
        origin,
        step,
        base,
        cells,
        Arithmetic::Canonical,
        &mut ambiguous,
    );
    let direct = mesh_block(
        field,
        origin,
        step,
        base,
        cells,
        Arithmetic::Offset,
        &mut ambiguous,
    );
    assert_eq!(
        clipped.triangle_count(),
        direct.triangle_count(),
        "the clip kept {} triangles where the block itself emits {}",
        clipped.triangle_count(),
        direct.triangle_count()
    );
    assert_eq!(
        referenced_positions(&clipped),
        referenced_positions(&direct),
        "the clip kept the right number of triangles and not the right ones"
    );
}

/// One row: one field, one resolution, one LOD pair, one arithmetic.
fn measure<F>(
    field: &F,
    name: &str,
    cells_per_axis: u32,
    level: u32,
    arm: Arithmetic,
    run: &mut common::experiment::Run,
) where
    F: ReferenceField + Sdf<Scalar = f64>,
{
    let fx = fixture(field, cells_per_axis, level);
    let step = fx.fine_step;
    clip_agrees_with_the_block(field, fx.origin, step, fx.fine_base, fx.fine_cells);
    clip_agrees_with_the_block(
        field,
        fx.origin,
        fx.coarse_step(),
        fx.coarse_base,
        fx.coarse_cells,
    );

    let mut ambiguous = 0usize;
    let fine = mesh_block(
        field,
        fx.origin,
        step,
        fx.fine_base,
        fx.fine_cells,
        arm,
        &mut ambiguous,
    );
    let coarse = mesh_block(
        field,
        fx.origin,
        fx.coarse_step(),
        fx.coarse_base,
        fx.coarse_cells,
        arm,
        &mut ambiguous,
    );
    assert_eq!(
        ambiguous, 0,
        "{name} n={cells_per_axis} level={level}: {ambiguous} triangles lie in a block's own \
         cut plane, so the clip cannot say which cell emitted them"
    );
    assert!(
        fine.triangle_count() > 0 && coarse.triangle_count() > 0,
        "{name} n={cells_per_axis} level={level}: a block missed the surface entirely"
    );

    let (patch, patches) = fx.transition_patches(field, arm);
    let seam_x = fx.seam_plane(arm);
    let seam_plane_delta = (seam_x[0] - seam_x[1]).abs();

    let stitched = [&fine, &coarse, &patch];
    let bare = [&fine, &coarse];
    let (welded, joined) = seam_report(&stitched, seam_x, step, Merge::Weld);
    let (shared, _) = seam_report(&stitched, seam_x, step, Merge::Exact);
    let (welded_bare, _) = seam_report(&bare, seam_x, step, Merge::Weld);
    let (shared_bare, _) = seam_report(&bare, seam_x, step, Merge::Exact);
    let coincide = coincidence(&fine, &[&coarse, &patch], seam_x, step);

    let lod_pair = format!("{level}-{}", level + 1);
    println!(
        "{name:<8} {cells_per_axis:<5} {lod_pair:<4} {:<10} {:>6} {:>6} {:>8} {:>8} {:>6} {:>6} \
         {:>10.3e} {:>10.3e}",
        arm.name(),
        welded.cracks,
        shared.cracks,
        welded_bare.cracks,
        shared_bare.cracks,
        coincide.pairs,
        coincide.exact,
        coincide.gap,
        seam_plane_delta
    );

    run.record(&[
        ("cell_size", format!("{step:.17e}")),
        ("lod_pair", lod_pair),
        ("crack_count", welded.cracks.to_string()),
        ("max_discontinuity", format!("{:.6e}", welded.discontinuity)),
        ("arithmetic", arm.name().to_string()),
        ("field", name.to_string()),
        ("cells_per_axis", cells_per_axis.to_string()),
        ("power_of_two", u8::from(is_power_of_two(step)).to_string()),
        ("crack_count_exact_merge", shared.cracks.to_string()),
        (
            "max_discontinuity_exact_merge",
            format!("{:.6e}", shared.discontinuity),
        ),
        ("crack_count_no_transition", welded_bare.cracks.to_string()),
        (
            "max_discontinuity_no_transition",
            format!("{:.6e}", welded_bare.discontinuity),
        ),
        (
            "crack_count_exact_merge_no_transition",
            shared_bare.cracks.to_string(),
        ),
        ("seam_pairs", coincide.pairs.to_string()),
        ("exact_pairs", coincide.exact.to_string()),
        ("max_pair_gap", format!("{:.3e}", coincide.gap)),
        ("seam_plane_delta", format!("{seam_plane_delta:.3e}")),
        ("patches", patches.to_string()),
        ("triangles", joined.triangle_count().to_string()),
    ]);
}

fn main() {
    if !std::env::args().any(|arg| arg == "--bench") {
        return;
    }

    let prereg = isomesh::experiment!("P-11");
    common::experiment::run(prereg, |run| {
        println!(
            "{:<8} {:<5} {:<4} {:<10} {:>13} {:>17} {:>6} {:>6} {:>10} {:>10}",
            "", "", "", "", "with transition", "no transition", "", "", "", ""
        );
        println!(
            "{:<8} {:<5} {:<4} {:<10} {:>6} {:>6} {:>8} {:>8} {:>6} {:>6} {:>10} {:>10}",
            "field",
            "n",
            "lod",
            "arithmetic",
            "weld",
            "exact",
            "weld",
            "exact",
            "pairs",
            "shared",
            "pair gap",
            "plane Δ"
        );
        for n in RESOLUTIONS {
            for level in FINE_LEVELS {
                for arm in Arithmetic::ALL {
                    measure(&Sphere::<f64>::canonical(), "sphere", n, level, arm, run);
                    measure(&Torus::<f64>::canonical(), "torus", n, level, arm, run);
                }
            }
        }
        println!(
            "\nEach crack column is seam-plane boundary edges, under one merge: `weld` is the \
             crate's\nown epsilon policy, `exact` is bit-identity — what a consumer that never \
             welds gets\n(M-69's collider). `pairs` counts the fine block's seam vertices with a \
             partner within\nthe weld epsilon; `shared` how many of those are bit-identical. \
             `plane Δ` is the\nfixture's own check that the two arms reached different bits at all."
        );
    });
}
