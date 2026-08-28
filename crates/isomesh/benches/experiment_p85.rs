//! **P-85 - the collider's 45%, attributed.**
//!
//! Ticket: R-085. Pre-registered before this harness existed.
//!
//! ```bash
//! cargo bench --bench experiment_p85
//! ```
//!
//! Writes `docs/experiments/p-85.csv`.
//!
//! # SHARE, recomputed before the code was written
//!
//! **This experiment moves nothing, so `✗51`'s reachability arithmetic has no
//! ratio to bite on.** The registration says so itself - *"SHARE: this
//! experiment measures shares and moves nothing"* - and the recomputation
//! confirms there is no `1/(1 − share/factor)` to clear:
//!
//! - **C1** asks whether one of four stages exceeds **50%** of the collider
//!   cost. The vacuity control caps the residual at 5% of the total, so the four
//!   shares sum to at least 0.95. A single share above 0.50 therefore requires
//!   one stage to exceed the other three combined; nothing in the arithmetic
//!   forbids it, so the bar is **reachable**.
//! - **C2** is nominal - *which* stage is largest - and is reachable as soon as
//!   the four stages are separately instrumented, which is the whole harness.
//! - **C3** is a ratio **between two arms** (`gyroid` against `sphere` at fixed
//!   triangle count), not a share of a total, so it is unbounded above and the
//!   **1.5×** bar is reachable.
//!
//! The number this row *informs* is `M-135`'s **45%**: the largest stage of a
//! usable mesh, measured there as `collider::readiness` at 45.0% mean share
//! against contour 29.0%, weld 25.5% and normals 0.4%.
//!
//! # What "the collider stage" is here, and the discrepancy that has to be said
//! out loud
//!
//! `M-135`'s 45% is **not** a collider build. It is `collider::readiness`, the
//! T-001 validity walk, and nothing in this repository builds a `parry3d`
//! `TriMesh` outside tests and benches - `parry3d` is a dev-dependency by
//! design (see `collider.rs`), and `game_dig`'s body resolves against the
//! **field** rather than against triangles. So the four registered stages
//! describe a **prospective** per-chunk handover, and the honest framing is:
//!
//! | stage | what it is here |
//! |---|---|
//! | `handoff` | everything between the `MeshBuffer` and the constructor: the seam **weld** (`M-69`, `✗18`) plus the `collider::readiness` gate, which is `M-135`'s 45% stage |
//! | `copy` | `MeshBuffer` positions to `Vec<parry3d::math::Vector>` plus `collider::triangle_indices` |
//! | `construct` | `TriMesh::new` **minus** the BVH build it performs internally |
//! | `bvh` | `TriMesh::rebuild_bvh`'s own body, `Bvh::from_iter(BvhBuildStrategy::Binned, …)` |
//!
//! The weld is inside `handoff` because the registration's own hypothesis names
//! it as a candidate answer (*"or the weld that `M-69` and `✗18` argue
//! about"*), and a stage named in the hypothesis that is outside the total
//! cannot be the answer to it. `weld_ms` and `readiness_ms` are recorded
//! separately so `handoff` is not a black box.
//!
//! # Why `construct` is a difference and `bvh` is a probe
//!
//! `parry3d` 0.30.2 builds the BVH **eagerly inside the constructor and offers
//! no way not to**. `TriMesh::new` is `with_flags(v, i, TriMeshFlags::empty())`,
//! which checks for an empty index buffer, moves the two buffers into the
//! struct, calls `set_flags(empty())` - which builds nothing - and then, because
//! `result.bvh.is_empty()`, calls `rebuild_bvh()`
//! (`parry3d-0.30.2/src/shape/trimesh.rs:724-751`). There is no `TriMesh`
//! without a BVH.
//!
//! So `bvh_ms` is measured by a probe that reproduces `rebuild_bvh`'s body
//! exactly - the same `enumerate` over the same `[u32; 3]` triples, the same
//! `Triangle::local_aabb`, the same `BvhBuildStrategy::Binned`
//! (`trimesh.rs:1163-1175`) - over the buffers the constructor was just handed,
//! read back through `trimesh.vertices()` and `trimesh.indices()`. The probe
//! runs **outside** the total window, and `construct_ms` is the constructor's
//! wall time minus it. `✗52` did the same thing for `submit = indirect − execute`
//! and said so; this says so too, and adds two things that entry did not have:
//! `trimesh_new_ms` is recorded, so the subtraction is reversible from the CSV,
//! and `construct_ms` is **not clamped** - a negative value would mean the
//! constructor's non-BVH part is below the noise floor, which is information
//! rather than an error.
//!
//! # The vacuity control, and exactly how far its teeth reach
//!
//! Registered: *"the four stages must sum to the measured total within 5%,
//! reported as a residual column, or the decomposition is missing a stage."*
//!
//! `residual_ms = total_ms − (copy + construct + bvh + handoff)` is recorded and
//! `assert!`ed at 5% of `total_ms`. Because `construct + bvh` is identically
//! `trimesh_new_ms`, the residual reduces to
//! `total − weld − readiness − copy − trimesh_new`, which is precisely "is
//! there a segment of the real build nobody timed".
//!
//! **A zero residual has to prove it could have been non-zero (`M-44`).** The
//! proof is a named, plausible omission rather than an arithmetic identity: the
//! weld is the stage a four-item list called *copy / construct / bvh / gate*
//! does not contain, and it is the stage the registration's hypothesis puts in
//! contention. `residual_share_without_weld` is recorded and asserted to
//! **exceed** the 5% bar, so the instrument is shown firing on the omission it
//! exists to catch. `min_stage_share` and `stages_above_bar` state the limit
//! from the other side: a stage below 5% of the total could be dropped without
//! the residual noticing, and the CSV names which ones those are rather than
//! implying the control is total.
//!
//! # The timing protocol, and the two defects that shaped it
//!
//! **Nine timed repeats after three warmups, with the repeat as the OUTER loop
//! and the arm as the inner one.** All six arms are built inside each repeat, so
//! a machine excursion lands across every arm and a ratio between two of them
//! survives it. That is not a precaution, it is a repair: on the first
//! sequential version, with five sibling agents building on this host, the
//! sphere arm's wall time moved **58%** between two runs of the same binary
//! while its cycle count moved **1.5%** - `M-280` arriving uninvited. The price
//! of interleaving is that each build starts on a cache the other five arms just
//! evicted, which raises every absolute by about **1.4x** uniformly and leaves
//! the shares where they were; and it is what a chunked pipeline building chunk
//! after chunk actually does.
//!
//! **The row is the median repeat, reported whole, not five independent
//! medians.** The second version took a median per stage and computed
//! `residual = median(total) − Σ median(stage)`; that read **+4.33%** on
//! `gyroid` at 65³ - inside the registered bar and pointing at a missing stage -
//! while `residual_share_worst_rep`, the largest residual any *single* repeat
//! produced, read **0.000%**. Every repeat's stages summed to its own total
//! exactly; the 4.33% was five medians disagreeing about which repeat they came
//! from, and it would have made the registered columns fail to add up in the
//! CSV. So the reported repeat is the one whose **total** is the median, and
//! `copy + construct + bvh + handoff + residual == total` on every row to the
//! bit. Per-stage medians over all nine repeats are kept as `*_median_ms`, and
//! `median_mixing_share` states how far apart the two views are.
//!
//! # C3's fixture: fixed triangle count
//!
//! `gyroid` and `sphere` cannot be compared at fixed triangle count by holding
//! the grid, because the two fields have wildly different surface area per unit
//! volume - `gyroid` is a surface everywhere inside its cap, `sphere` is one
//! shell. So the **sphere's resolution is searched**: a monotone bisection over
//! samples per axis for the resolution whose triangle count is closest to the
//! `gyroid` row it is paired with, then a ±3 scan around it because a triangle
//! count is only nearly monotone in resolution.
//!
//! **The matched quantity is the POST-WELD count**, and the harness's own
//! `assert_eq!` is why: matching on the contour's count put 10,632 in the search
//! against 10,628 on the row, because welding collapses duplicate vertices and a
//! triangle whose corners collapse together stops being a triangle. Parry is
//! handed the welded mesh, so that is what "fixed triangle count" has to mean.
//! `triangle_mismatch` and `triangle_mismatch_share` are recorded and the share
//! is asserted under 3%; `cost_per_triangle` is the same control done
//! dimensionlessly and is what C3's verdict is read from.
//!
//! The residual confound is named rather than hidden: the matched sphere grid is
//! **larger** than the gyroid grid it matches, so its contour touched more cells
//! and its cache state entering the collider build is not identical. The
//! collider stage's input is the mesh and not the grid, so the confound is
//! limited to cache residency; `chunk_cells` on the row says which grid it was.
//!
//! # Controls, each an assertion rather than a printed number
//!
//! - **The residual, at 5%** - the registered control, above.
//! - **The residual fires** - `residual_share_without_weld > 0.05`.
//! - **The mesh is real** - `ready.is_usable()` and a non-zero triangle count on
//!   every arm, so no row reports a fast time for building a collider out of
//!   nothing.
//! - **The seam fixture reaches a seam** - the two-chunk probe must report a
//!   non-zero duplicate-vertex count before the weld and strictly fewer
//!   boundary edges after it, which is `M-69`'s measurement (36 duplicates, 180
//!   boundary edges, 72 of them the seam's) reproduced per field.
//! - **The stage ordering is stable** - the largest stage is recomputed per
//!   repeat and `reps_agreeing_on_largest` records how many of the nine agree
//!   with the row. A row whose ordering is decided by one repeat is not an
//!   attribution.
//! - **The gyroid arm is the mesh the sphere was matched to** - `assert_eq!` on
//!   the timed arm's welded triangle count against the fixture phase's. This is
//!   the control that found the pre-weld/post-weld mismatch above.
//! - **The matched triangle counts really match** - `triangle_mismatch_share`
//!   under 3%; measured **0.677%** at the 33³ pair and **0.223%** at 65³.
//! - **The machine held still** - `counted_over_median`, the counted build's
//!   wall time over the reported repeat's. A row where it is far from 1 was
//!   measured through someone else's build.
//! - **Nothing was multiplexed** - the counted pass asserts
//!   `Counts::worst_ratio() >= MIN_TIME_RATIO`, so `ghz` is a reading and not an
//!   extrapolation.
//!
//! # Units, and why there is a clock on the row
//!
//! `M-280`: on a governed CPU a nanosecond is not a unit. Every row carries
//! `ghz`, computed as cycles divided by nanoseconds over the same span, plus
//! `cycles_per_triangle`, which is clock-independent. `cost_per_triangle` is
//! **microseconds per triangle** and moves with the governor; the shares and the
//! ratios do not, which is why the verdicts are read from those, and
//! `c3_ratio_cycles` carries C3's ratio in cycles beside the registered
//! wall-clock one. Every arm is measured in one binary, one process and one
//! interleaved set of repeats (`M-281`).

#![allow(
    clippy::cast_possible_truncation,
    clippy::cast_precision_loss,
    clippy::cast_sign_loss,
    clippy::too_many_lines
)]

mod common;

#[cfg(target_os = "linux")]
mod experiment {
    use std::hint::black_box;
    use std::time::{Duration, Instant};

    use parry3d::math::Vector;
    use parry3d::partitioning::{Bvh, BvhBuildStrategy};
    use parry3d::shape::{TriMesh, Triangle};

    use isomesh::chunk::{ChunkId, ChunkLayout};
    use isomesh::collider::{self, ColliderReadiness};
    use isomesh::fields::{FbmTerrain, ReferenceField, Sphere, capped_gyroid};
    use isomesh::marching_cubes::MarchingCubes;
    use isomesh::validate::ValidateConfig;
    use isomesh::weld::{self, Welder};
    use isomesh::{MeshBuffer, RuntimeShape3, Sdf};

    use crate::common::counters::{MIN_TIME_RATIO, Probe};
    use crate::common::experiment::Run;
    use crate::common::grid;

    /// `f32`, because that is the only width `parry3d` has.
    ///
    /// `M-135` measured its 45% at `f64` and this is the one deliberate
    /// deviation from that fixture. It is the right one: `parry3d::math::Vector`
    /// is `glam::Vec3`, so an `f64` mesh would pay a narrowing in `copy` that no
    /// shipped path pays, and `bevy_isomesh` is `f32` throughout. The stages that
    /// dominate - the weld's hash grid and the validator's edge walk - are driven
    /// by indices rather than by scalar width, so the comparison with `M-135`
    /// survives the change.
    type Scalar = f32;

    /// The registered arms: 33³ and 65³ chunks.
    ///
    /// Samples per axis, which is what `M-135`'s own fixture and
    /// `common::grid` mean by 33³ - so `cells_per_axis` on the row is 32 and 64.
    const RESOLUTIONS: [u32; 2] = [33, 65];
    /// Builds discarded before timing, so the constructor's first-touch page
    /// faults are not in the medians.
    const WARMUP: usize = 3;
    /// Timed builds per arm, median taken.
    ///
    /// **Nine, and the number is justified by a column rather than by taste.**
    /// C1 and C2 are decided by an *ordering*, not by a magnitude, so what has
    /// to be stable is which stage is largest - and `M-337`'s re-audit is the
    /// standing warning about a single measurement on a governed CPU. The
    /// harness recomputes the largest stage for each of the nine and records
    /// `reps_agreeing_on_largest`, so the claim "the ordering is stable" is a
    /// number in the file.
    const REPS: usize = 9;
    /// The registered residual bar.
    const RESIDUAL_BAR: f64 = 0.05;
    /// C1's bar.
    const DOMINANT_BAR: f64 = 0.50;
    /// C3's bar.
    const FIELD_PENALTY_BAR: f64 = 1.5;
    /// How far the matched sphere may miss the gyroid's triangle count.
    const MATCH_TOLERANCE: f64 = 0.03;

    /// The four stages, in the order the registration names them.
    const STAGES: [&str; 4] = ["copy", "construct", "bvh", "handoff"];

    /// One build's stage times, and what the mesh was.
    struct Stages {
        weld: Duration,
        readiness: Duration,
        copy: Duration,
        /// `TriMesh::new`'s whole wall time, BVH included.
        construct_and_bvh: Duration,
        /// The window the four stages must sum to.
        total: Duration,
        ready: ColliderReadiness,
        vertices: usize,
        triangles: usize,
        /// Bytes written by `copy`. The read side is the same size.
        copy_bytes: usize,
    }

    /// The real per-chunk collider build, timed segment by segment.
    ///
    /// `work` arrives already cloned, because the clone is the harness's cost and
    /// not the build's. The `TriMesh` is returned rather than dropped so that
    /// freeing it lands outside `total`.
    fn build_once(
        mut work: MeshBuffer<Scalar>,
        cell_size: Scalar,
        cfg: &ValidateConfig,
    ) -> (Stages, TriMesh) {
        let t_total = Instant::now();

        // `handoff`, first half: the seam weld. A chunked collider must be
        // welded or parry reads the seam as a hole (M-69: 72 boundary edges per
        // seam, closed by the weld and nothing else).
        let t = Instant::now();
        Welder::<Scalar>::new()
            .weld(&mut work, weld::epsilon_for(cell_size))
            .expect("the welder accepts a mesh Marching Cubes just produced");
        let weld_t = t.elapsed();

        // `handoff`, second half: the readiness gate. This is M-135's 45% stage.
        let t = Instant::now();
        let ready = collider::readiness(&work, cfg);
        let readiness_t = t.elapsed();

        // `copy`: the conversion `collider.rs`'s module docs write out, verbatim.
        let t = Instant::now();
        let vertices: Vec<Vector> = work
            .positions
            .iter()
            .map(|p| Vector::new(p[0], p[1], p[2]))
            .collect();
        let indices = collider::triangle_indices(&work);
        let copy_t = t.elapsed();

        // `construct` + `bvh`: parry's constructor, which builds the BVH inside
        // itself and cannot be asked not to.
        let copy_bytes = vertices.len() * size_of::<Vector>() + indices.len() * size_of::<[u32; 3]>();
        let t = Instant::now();
        let trimesh = TriMesh::new(vertices, indices)
            .expect("parry accepts a non-empty index buffer, which readiness already asserted");
        let construct_and_bvh = t.elapsed();

        let total = t_total.elapsed();
        let stages = Stages {
            weld: weld_t,
            readiness: readiness_t,
            copy: copy_t,
            construct_and_bvh,
            total,
            ready,
            vertices: work.positions.len(),
            triangles: work.triangle_count(),
            copy_bytes,
        };
        (stages, trimesh)
    }

    /// `TriMesh::rebuild_bvh`'s body, over the buffers the constructor holds.
    ///
    /// Transcribed from `parry3d-0.30.2/src/shape/trimesh.rs:1163-1175`: the same
    /// `enumerate`, the same `Triangle::local_aabb`, the same
    /// `BvhBuildStrategy::Binned`. Run outside the total window, so subtracting
    /// it from the constructor's wall time leaves the residual measuring untimed
    /// work rather than this probe.
    fn bvh_once(trimesh: &TriMesh) -> Duration {
        let vertices = trimesh.vertices();
        let indices = trimesh.indices();
        let t = Instant::now();
        let bvh = Bvh::from_iter(
            BvhBuildStrategy::Binned,
            indices.iter().enumerate().map(|(i, idx)| {
                let aabb = Triangle::new(
                    vertices[idx[0] as usize],
                    vertices[idx[1] as usize],
                    vertices[idx[2] as usize],
                )
                .local_aabb();
                (i, aabb)
            }),
        );
        let elapsed = t.elapsed();
        black_box(&bvh);
        elapsed
    }

    /// Contour one field at `samples` per axis.
    fn contour<F: Sdf<Scalar = Scalar> + ReferenceField>(
        field: &F,
        samples: u32,
    ) -> (MeshBuffer<Scalar>, Scalar) {
        let (shape, lo, cell_size): (RuntimeShape3, [Scalar; 3], Scalar) = grid(field, samples);
        let mut mesh = MeshBuffer::<Scalar>::new();
        MarchingCubes::<Scalar>::new()
            .extract(field, &shape, lo, cell_size, &mut mesh)
            .expect("extraction on a reference field's own grid");
        (mesh, cell_size)
    }

    /// The triangle count a collider build actually receives.
    ///
    /// **Post-weld, and the difference is not cosmetic.** C3's fixture matched
    /// the sphere's resolution to the gyroid's *contour* count on its first run
    /// and the `assert_eq!` in `run` caught it: the search said 10,632 and the
    /// row said 10,628, because welding collapses duplicate vertices and a
    /// triangle whose corners collapse onto one another becomes a repeated-index
    /// triangle the count no longer includes (`M-185`'s mechanism). Four
    /// triangles at 33³ and four at 65³ on `gyroid`, two at each on
    /// `fbm_terrain` - small, and the wrong quantity is still the wrong
    /// quantity: what parry is handed is the welded mesh, so that is what
    /// "fixed triangle count" has to mean.
    fn welded_triangles<F: Sdf<Scalar = Scalar> + ReferenceField>(
        field: &F,
        samples: u32,
    ) -> usize {
        let (mut mesh, cell_size) = contour(field, samples);
        Welder::<Scalar>::new()
            .weld(&mut mesh, weld::epsilon_for(cell_size))
            .expect("the welder accepts a mesh Marching Cubes just produced");
        mesh.triangle_count()
    }

    /// The same contour, timed, so the Godot claim can be checked in one build.
    ///
    /// Godot Voxel Tools' performance documentation is the registration's cited
    /// independent corroboration - *"creating a collider from a mesh is actually
    /// much more expensive than meshing itself (about 3 to 5 times)"*, with no
    /// absolute timings and no hardware. Checking it against `M-135`'s committed
    /// `contour_ms` would be a cross-build ratio and `M-281` forbids that, so the
    /// contour is re-timed here, in this binary, beside the collider build it is
    /// the denominator of.
    fn contour_timed<F: Sdf<Scalar = Scalar> + ReferenceField>(
        field: &F,
        samples: u32,
    ) -> (MeshBuffer<Scalar>, Scalar, f64) {
        let (shape, lo, cell_size): (RuntimeShape3, [Scalar; 3], Scalar) = grid(field, samples);
        let mut mesh = MeshBuffer::<Scalar>::new();
        let mut times = Vec::with_capacity(REPS);
        for rep in 0..(WARMUP + REPS) {
            mesh.reset();
            let t = Instant::now();
            MarchingCubes::<Scalar>::new()
                .extract(field, &shape, lo, cell_size, &mut mesh)
                .expect("extraction on a reference field's own grid");
            let elapsed = t.elapsed();
            if rep >= WARMUP {
                times.push(elapsed.as_secs_f64() * 1000.0);
            }
            black_box(mesh.triangle_count());
        }
        (mesh, cell_size, median(times))
    }

    /// `M-69`'s seam, per field, at this cell size.
    ///
    /// Two chunks of the same `samples − 1` cells, adjacent in X, placed so the
    /// seam plane is the domain centre and the slab is centred transversely -
    /// otherwise a compact field's surface misses the seam entirely and the
    /// column is `M-44`'s zero. Returns `(seam boundary edges, duplicate
    /// vertices before the weld)`, where the first is
    /// `boundary_edges(unwelded) − boundary_edges(welded)`: exactly the quantity
    /// `M-69` isolated when it found 180 boundary edges falling to 108, the 72
    /// removed being the seam and the 108 remaining being the slab's own open
    /// border.
    fn seam<F: Sdf<Scalar = Scalar> + ReferenceField>(field: &F, samples: u32) -> (u64, u64) {
        let cells = samples - 1;
        let (lo, hi) = field.domain();
        let cell_size = (hi[0] - lo[0]) / cells as Scalar;
        let extent = cell_size * cells as Scalar;
        let layout = ChunkLayout::<Scalar>::new(
            cells,
            cell_size,
            [-extent, -extent * 0.5, -extent * 0.5],
        )
        .expect("a chunk layout at the field's own cell size");
        let shape = layout.sample_shape().expect("the chunk's sample shape");

        let mut joined = MeshBuffer::<Scalar>::new();
        for id in [ChunkId::new([0, 0, 0]), ChunkId::new([1, 0, 0])] {
            let mut chunk = MeshBuffer::<Scalar>::new();
            MarchingCubes::<Scalar>::new()
                .extract(
                    field,
                    &shape,
                    layout.sample_origin(id),
                    layout.cell_size(),
                    &mut chunk,
                )
                .expect("extraction of one seam chunk");
            joined
                .append(&chunk)
                .expect("two chunks fit the u32 index space");
        }

        let cfg = ValidateConfig::from_cell_size(f64::from(cell_size))
            .expect("a positive finite cell size");
        let before = collider::readiness(&joined, &cfg);
        assert!(
            before.duplicate_vertices > 0,
            "the two chunks met with no duplicate vertices at all, so the seam fixture is not \
             exercising a seam and seam_boundary_edges would be M-44's zero"
        );
        Welder::<Scalar>::new()
            .weld(&mut joined, weld::epsilon_for(cell_size))
            .expect("the welder accepts the joined slab");
        let after = collider::readiness(&joined, &cfg);
        assert!(
            after.boundary_edges < before.boundary_edges,
            "welding closed no boundary edge, so the seam carried none: {before:?} -> {after:?}"
        );
        (
            before.boundary_edges - after.boundary_edges,
            before.duplicate_vertices,
        )
    }

    /// Everything one arm produced.
    struct Arm {
        field: &'static str,
        arm: &'static str,
        samples: u32,
        cell_size: f64,
        vertices: usize,
        triangles: usize,
        copy_ms: f64,
        construct_ms: f64,
        bvh_ms: f64,
        handoff_ms: f64,
        total_ms: f64,
        weld_ms: f64,
        readiness_ms: f64,
        trimesh_new_ms: f64,
        residual_ms: f64,
        residual_share: f64,
        residual_share_worst_rep: f64,
        residual_share_without_weld: f64,
        largest_stage: &'static str,
        largest_share: f64,
        /// The largest of the **five** stages the harness can separate: the four
        /// registered ones with `handoff` split back into weld and gate.
        ///
        /// Recorded because C1's verdict turns on the registration's decision to
        /// call "anything between" one stage, and a reader has to be able to see
        /// that from the file rather than from prose.
        finest_largest_stage: &'static str,
        finest_largest_share: f64,
        min_stage_share: f64,
        stages_above_bar: usize,
        reps_agreeing: usize,
        cost_per_triangle: f64,
        contour_ms: f64,
        /// The wall time of the **counted** build, the one `cycles` came from.
        ///
        /// Recorded so a reader can see whether the machine held still: if this
        /// disagrees with the median of the nine timed reps, something outside
        /// this process moved during the row and the wall-clock columns carry it.
        /// `M-280` says report the clock; this says report whether the clock was
        /// the only thing moving.
        counted_total_ms: f64,
        cycles: u64,
        ghz: f64,
        cycles_per_triangle: f64,
        copy_bytes: usize,
        copy_gb_per_s: f64,
        seam_boundary_edges: u64,
        seam_duplicate_vertices: u64,
        unwelded_duplicate_vertices: u64,
        unwelded_boundary_edges: u64,
        degenerate_triangles: u64,
        boundary_edges: u64,
        duplicate_vertices: u64,
        /// Per-stage medians over all nine repeats, in `STAGES` order plus weld
        /// and readiness: the robust view, kept beside the median repeat's own
        /// numbers so both are in the file.
        medians: [f64; 6],
    }

    fn median(mut xs: Vec<f64>) -> f64 {
        xs.sort_unstable_by(f64::total_cmp);
        xs[xs.len() / 2]
    }

    /// Which of the four is largest, and its share.
    fn largest(stages: [f64; 4], total: f64) -> (&'static str, f64) {
        let mut best = 0;
        for i in 1..4 {
            if stages[i] > stages[best] {
                best = i;
            }
        }
        (STAGES[best], stages[best] / total)
    }

    /// The finest split the harness can make: `handoff` back into its two halves.
    ///
    /// Not a fifth registered stage - the registration names four and this
    /// harness reports four. This is the answer to "and if `handoff` were not one
    /// stage?", which is the one question C1's verdict is sensitive to, recorded
    /// so the sensitivity is in the artefact.
    fn finest_largest(stages: [f64; 5], total: f64) -> (&'static str, f64) {
        const FINEST: [&str; 5] = ["copy", "construct", "bvh", "weld", "readiness"];
        let mut best = 0;
        for i in 1..5 {
            if stages[i] > stages[best] {
                best = i;
            }
        }
        (FINEST[best], stages[best] / total)
    }

    /// One arm's fixture: everything that can be built before the clock matters.
    ///
    /// Field-typed work ends here. Six of these exist at once, so the timed loop
    /// can visit all six per repeat without needing the fields again.
    struct Prepared {
        field: &'static str,
        arm: &'static str,
        samples: u32,
        cell_size: Scalar,
        cfg: ValidateConfig,
        pristine: MeshBuffer<Scalar>,
        contour_ms: f64,
        seam_boundary_edges: u64,
        seam_duplicate_vertices: u64,
        /// What a **single** chunk's mesh looks like *before* the weld.
        ///
        /// This is the number that decides how to read `weld_share`. The seam
        /// probe proves a two-chunk join needs the weld (`M-69`); it says nothing
        /// about one chunk, and `duplicate_vertices` on the row is measured after
        /// the weld and is therefore zero by construction - `M-44`'s zero exactly.
        /// So the unwelded mesh is read once, here, untimed.
        unwelded_duplicate_vertices: u64,
        unwelded_boundary_edges: u64,
    }

    /// One repeat's five numbers, in milliseconds.
    struct Rep {
        weld: f64,
        readiness: f64,
        copy: f64,
        trimesh_new: f64,
        bvh: f64,
        total: f64,
    }

    fn prepare<F: Sdf<Scalar = Scalar> + ReferenceField>(
        field: &F,
        name: &'static str,
        arm: &'static str,
        samples: u32,
    ) -> Prepared {
        let (pristine, cell_size, contour_ms) = contour_timed(field, samples);
        assert!(
            pristine.triangle_count() > 0,
            "{name} at {samples}³ contoured nothing, so every stage below would be timing an \
             empty buffer"
        );
        let cfg = ValidateConfig::from_cell_size(f64::from(cell_size))
            .expect("a positive finite cell size");
        let (seam_boundary_edges, seam_duplicate_vertices) = seam(field, samples);
        // What one chunk looks like before the weld, so `weld_share` can be read.
        let unwelded = collider::readiness(&pristine, &cfg);
        assert!(
            unwelded.is_usable(),
            "{name} at {samples}³ contoured a mesh no engine would take: {unwelded:?}"
        );
        Prepared {
            field: name,
            arm,
            samples,
            cell_size,
            cfg,
            pristine,
            contour_ms,
            seam_boundary_edges,
            seam_duplicate_vertices,
            unwelded_duplicate_vertices: unwelded.duplicate_vertices,
            unwelded_boundary_edges: unwelded.boundary_edges,
        }
    }

    /// Build one arm's collider once, and read its clocks.
    fn one_rep(p: &Prepared) -> (Rep, Stages) {
        let (stages, trimesh) = build_once(p.pristine.clone(), p.cell_size, &p.cfg);
        let bvh = bvh_once(&trimesh);
        assert!(
            stages.ready.is_usable(),
            "{} at {}³ produced a mesh no engine would take: {:?}",
            p.field,
            p.samples,
            stages.ready
        );
        let ms = |d: Duration| d.as_secs_f64() * 1000.0;
        let rep = Rep {
            weld: ms(stages.weld),
            readiness: ms(stages.readiness),
            copy: ms(stages.copy),
            trimesh_new: ms(stages.construct_and_bvh),
            bvh: ms(bvh),
            total: ms(stages.total),
        };
        (rep, stages)
    }

    /// One arm's counted build: cycles, and the nanoseconds they were spent in.
    struct Counted {
        cycles: u64,
        ns: f64,
    }

    /// Build one arm once with the hardware counters running.
    ///
    /// Called from its **own interleaved wave**, not from `finish`. The first
    /// version counted inside `finish`, which runs after every timed wave has
    /// ended - so the counted build saw a machine in a different state than the
    /// repeats it is compared with, and `counted_over_median` read **1.36** on a
    /// run whose residuals were 0.00002%. A control that fires on its own
    /// scheduling is not a control.
    fn count_once(p: &Prepared, probe: &mut Probe) -> Counted {
        probe.reset_and_enable();
        let (stages, trimesh) = build_once(p.pristine.clone(), p.cell_size, &p.cfg);
        probe.disable();
        let counts = probe.read();
        drop(trimesh);
        assert!(
            counts.worst_ratio() >= MIN_TIME_RATIO,
            "{} at {}³: a counter was multiplexed at ratio {:.4}, so `ghz` would be an \
             extrapolation",
            p.field,
            p.samples,
            counts.worst_ratio()
        );
        Counted {
            cycles: counts.cycles.count,
            ns: stages.total.as_secs_f64() * 1e9,
        }
    }

    /// Turn one arm's repeats into a row, with every control on the way.
    fn finish(p: &Prepared, reps: &[Rep], last: &Stages, counted: &Counted) -> Arm {
        let (name, samples) = (p.field, p.samples);
        let counted_ns = counted.ns;

        // ── the median repeat, reported whole ───────────────────────────────
        //
        // **Five independent medians are not one observation, and the residual
        // control caught that too.** The first version took a median per stage
        // and computed `residual = median(total) − Σ median(stage)`. Under
        // sibling load that read **+4.33%** on `gyroid` at 65³ - a hair inside
        // the registered 5% bar and pointing at a missing stage - while
        // `residual_share_worst_rep`, the largest residual any *single* repeat
        // produced, read **0.000%**. So every repeat's four stages summed to its
        // own total exactly and the 4.33% was the five medians disagreeing about
        // which repeat they came from. That is a real defect: the registered
        // columns would not add up, and a reader recomputing the residual from
        // the CSV would get a different number than the CSV states.
        //
        // The fix is to pick the repeat whose **total** is the median and report
        // that repeat's five stages. It is still a median - of the quantity the
        // shares are denominated in - and the row is now internally exact:
        // `copy + construct + bvh + handoff + residual == total`, to the bit.
        // Per-stage medians over all nine repeats are kept as `*_median_ms`
        // columns, so the robust view is in the file and nothing is lost.
        let mut order: Vec<usize> = (0..reps.len()).collect();
        order.sort_unstable_by(|a, b| reps[*a].total.total_cmp(&reps[*b].total));
        let m = &reps[order[reps.len() / 2]];
        let (weld, readiness, copy, trimesh_new, bvh, total) =
            (m.weld, m.readiness, m.copy, m.trimesh_new, m.bvh, m.total);
        let construct = trimesh_new - bvh;
        let handoff = weld + readiness;

        let residual = total - (copy + construct + bvh + handoff);
        let residual_share = residual / total;
        let (largest_stage, largest_share) = largest([copy, construct, bvh, handoff], total);
        let (finest_largest_stage, finest_largest_share) =
            finest_largest([copy, construct, bvh, weld, readiness], total);
        let shares = [copy / total, construct / total, bvh / total, handoff / total];
        let min_stage_share = shares.iter().copied().fold(f64::INFINITY, f64::min);
        let stages_above_bar = shares.iter().filter(|s| **s > RESIDUAL_BAR).count();
        let reps_agreeing = reps
            .iter()
            .filter(|r| {
                largest(
                    [r.copy, r.trimesh_new - r.bvh, r.bvh, r.weld + r.readiness],
                    r.total,
                )
                .0 == largest_stage
            })
            .count();
        let worst_residual_share = reps
            .iter()
            .map(|r| {
                ((r.total - (r.weld + r.readiness + r.copy + r.trimesh_new)) / r.total).abs()
            })
            .fold(0.0f64, f64::max);

        // ── the registered vacuity control ──────────────────────────────────
        //
        // Asserted twice: on the row, and on the worst of the nine repeats. The
        // second is the stronger statement - it says *no* repeat had untimed work
        // in it, not merely that the reported one did not.
        assert!(
            residual_share.abs() <= RESIDUAL_BAR,
            "{name} at {samples}³: the four stages miss the total by {:.2}% \
             (residual {residual:.4} ms of {total:.4} ms), so the decomposition is missing a \
             stage — copy {copy:.4}, construct {construct:.4}, bvh {bvh:.4}, handoff {handoff:.4}",
            residual_share * 100.0
        );
        assert!(
            worst_residual_share <= RESIDUAL_BAR,
            "{name} at {samples}³: the worst of {} repeats misses its own total by {:.2}%, so at \
             least one build had a segment nobody timed",
            reps.len(),
            worst_residual_share * 100.0
        );
        // ── and the proof it could have fired (`M-44`) ───────────────────────
        //
        // The weld is the stage a four-item list called copy/construct/bvh/gate
        // does not have, and the one the registration's hypothesis puts in
        // contention. If omitting it does not breach the bar, the residual is a
        // control that cannot catch the omission it exists for.
        let residual_share_without_weld = (residual + weld) / total;
        assert!(
            residual_share_without_weld > RESIDUAL_BAR,
            "{name} at {samples}³: dropping the weld from the decomposition leaves the residual \
             at {:.2}%, inside the 5% bar — so the residual control could not have caught a \
             missing weld and is vacuous on this row",
            residual_share_without_weld * 100.0
        );

        let triangles = last.triangles;
        Arm {
            field: name,
            arm: p.arm,
            samples,
            cell_size: f64::from(p.cell_size),
            vertices: last.vertices,
            triangles,
            copy_ms: copy,
            construct_ms: construct,
            bvh_ms: bvh,
            handoff_ms: handoff,
            total_ms: total,
            weld_ms: weld,
            readiness_ms: readiness,
            trimesh_new_ms: trimesh_new,
            residual_ms: residual,
            residual_share,
            residual_share_worst_rep: worst_residual_share,
            residual_share_without_weld,
            largest_stage,
            largest_share,
            finest_largest_stage,
            finest_largest_share,
            min_stage_share,
            stages_above_bar,
            reps_agreeing,
            cost_per_triangle: total * 1000.0 / triangles as f64,
            contour_ms: p.contour_ms,
            counted_total_ms: counted_ns / 1e6,
            cycles: counted.cycles,
            ghz: counted.cycles as f64 / counted_ns,
            cycles_per_triangle: counted.cycles as f64 / triangles as f64,
            copy_bytes: last.copy_bytes,
            copy_gb_per_s: last.copy_bytes as f64 / (copy / 1000.0) / 1e9,
            seam_boundary_edges: p.seam_boundary_edges,
            seam_duplicate_vertices: p.seam_duplicate_vertices,
            unwelded_duplicate_vertices: p.unwelded_duplicate_vertices,
            unwelded_boundary_edges: p.unwelded_boundary_edges,
            degenerate_triangles: last.ready.degenerate_triangles,
            boundary_edges: last.ready.boundary_edges,
            duplicate_vertices: last.ready.duplicate_vertices,
            medians: [
                median(reps.iter().map(|r| r.copy).collect()),
                median(reps.iter().map(|r| r.trimesh_new - r.bvh).collect()),
                median(reps.iter().map(|r| r.bvh).collect()),
                median(reps.iter().map(|r| r.weld + r.readiness).collect()),
                median(reps.iter().map(|r| r.total).collect()),
                median(reps.iter().map(|r| r.weld).collect()),
            ],
        }
    }

    /// The sphere resolution whose triangle count is closest to `target`.
    ///
    /// Monotone bisection, then a ±3 scan: a Marching Cubes triangle count grows
    /// as the square of the resolution but is only *nearly* monotone in it, so
    /// the bisection lands near the optimum and the scan finds it.
    fn sphere_samples_for(target: usize) -> u32 {
        let sphere = Sphere::<Scalar>::canonical();
        let count = |samples: u32| welded_triangles(&sphere, samples);

        let (mut lo, mut hi) = (9u32, 201u32);
        assert!(
            count(hi) >= target,
            "a unit sphere at {hi}³ emits fewer than the {target} triangles the gyroid arm needs, \
             so C3's fixed-triangle-count fixture is out of reach on this grid range"
        );
        while lo + 1 < hi {
            let mid = (lo + hi) / 2;
            if count(mid) < target {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        let mut best = hi;
        let mut best_gap = count(hi).abs_diff(target);
        for samples in hi.saturating_sub(3).max(9)..=(hi + 3) {
            let gap = count(samples).abs_diff(target);
            if gap < best_gap {
                best = samples;
                best_gap = gap;
            }
        }
        best
    }

    type Row = Vec<(&'static str, String)>;

    /// C3's comparison, computed once and written to **both** rows of a pair.
    ///
    /// **Every quantity here is gyroid ÷ sphere, in that order, on every row.**
    /// The first version computed the ratio as `self ÷ other` inside `row`, so
    /// `c3_ratio` read 1.026 on the gyroid row and 0.975 on the sphere row - the
    /// same column meaning two opposite things depending on which row you were
    /// looking at, which is `P-64`'s corrupt-but-plausible CSV in a new costume.
    /// C3's claim has a direction ("gyroid's collider costs at least 1.5x
    /// sphere's"), so the column has one too.
    struct Pairing {
        /// The *other* arm, as seen from the row being written.
        pair_field: &'static str,
        pair_chunk_cells: u32,
        /// Gyroid triangles minus sphere triangles.
        mismatch: i64,
        mismatch_share: f64,
        /// Gyroid microseconds per triangle divided by sphere's.
        ratio: f64,
        /// The same ratio in cycles per triangle.
        ratio_cycles: f64,
        held: bool,
    }

    fn row(a: &Arm, c1: bool, c2: bool, c3_all: bool, pair: Option<&Pairing>) -> Row {
        let mut r: Row = vec![
            ("field", a.field.to_string()),
            ("chunk_cells", a.samples.to_string()),
            ("triangles", a.triangles.to_string()),
            ("copy_ms", format!("{:.6}", a.copy_ms)),
            ("construct_ms", format!("{:.6}", a.construct_ms)),
            ("bvh_ms", format!("{:.6}", a.bvh_ms)),
            ("handoff_ms", format!("{:.6}", a.handoff_ms)),
            ("total_ms", format!("{:.6}", a.total_ms)),
            ("residual_ms", format!("{:.6}", a.residual_ms)),
            ("largest_stage", a.largest_stage.to_string()),
            ("largest_share", format!("{:.6}", a.largest_share)),
            ("cost_per_triangle", format!("{:.6}", a.cost_per_triangle)),
            ("seam_boundary_edges", a.seam_boundary_edges.to_string()),
            (
                "degenerate_triangles",
                a.degenerate_triangles.to_string(),
            ),
            ("c1_holds", c1.to_string()),
            ("c2_holds", c2.to_string()),
            // ── extras: the attribution, reversible from the file ───────────
            //
            // `c3_holds` is written **per pair** rather than as one aggregate,
            // because C3 is a statement about a pair and there are two of them.
            // The conjunction is here beside it.
            ("c3_holds_all_pairs", c3_all.to_string()),
            ("arm", a.arm.to_string()),
            ("scalar", "f32".to_string()),
            ("cells_per_axis", (a.samples - 1).to_string()),
            ("cell_size", format!("{:.6}", a.cell_size)),
            ("vertices", a.vertices.to_string()),
            ("weld_ms", format!("{:.6}", a.weld_ms)),
            ("readiness_ms", format!("{:.6}", a.readiness_ms)),
            ("trimesh_new_ms", format!("{:.6}", a.trimesh_new_ms)),
            ("copy_share", format!("{:.6}", a.copy_ms / a.total_ms)),
            (
                "construct_share",
                format!("{:.6}", a.construct_ms / a.total_ms),
            ),
            ("bvh_share", format!("{:.6}", a.bvh_ms / a.total_ms)),
            ("handoff_share", format!("{:.6}", a.handoff_ms / a.total_ms)),
            ("weld_share", format!("{:.6}", a.weld_ms / a.total_ms)),
            (
                "readiness_share",
                format!("{:.6}", a.readiness_ms / a.total_ms),
            ),
            (
                "bvh_share_of_construct",
                format!("{:.6}", a.bvh_ms / a.trimesh_new_ms),
            ),
            ("copy_median_ms", format!("{:.6}", a.medians[0])),
            ("construct_median_ms", format!("{:.6}", a.medians[1])),
            ("bvh_median_ms", format!("{:.6}", a.medians[2])),
            ("handoff_median_ms", format!("{:.6}", a.medians[3])),
            ("total_median_ms", format!("{:.6}", a.medians[4])),
            ("weld_median_ms", format!("{:.6}", a.medians[5])),
            (
                "median_mixing_share",
                format!(
                    "{:.6}",
                    (a.medians[4]
                        - (a.medians[0] + a.medians[1] + a.medians[2] + a.medians[3]))
                        / a.medians[4]
                ),
            ),
            (
                "finest_largest_stage",
                a.finest_largest_stage.to_string(),
            ),
            (
                "finest_largest_share",
                format!("{:.6}", a.finest_largest_share),
            ),
            ("contour_ms", format!("{:.6}", a.contour_ms)),
            (
                "collider_over_contour",
                format!("{:.6}", a.total_ms / a.contour_ms),
            ),
            (
                "counted_total_ms",
                format!("{:.6}", a.counted_total_ms),
            ),
            (
                "counted_over_median",
                format!("{:.6}", a.counted_total_ms / a.total_ms),
            ),
            ("bvh_dominates", (a.largest_stage == "bvh").to_string()),
            ("residual_share", format!("{:.6}", a.residual_share)),
            (
                "residual_share_worst_rep",
                format!("{:.6}", a.residual_share_worst_rep),
            ),
            (
                "residual_share_without_weld",
                format!("{:.6}", a.residual_share_without_weld),
            ),
            ("min_stage_share", format!("{:.6}", a.min_stage_share)),
            ("stages_above_bar", a.stages_above_bar.to_string()),
            ("reps", REPS.to_string()),
            ("reps_agreeing_on_largest", a.reps_agreeing.to_string()),
            ("copy_bytes", a.copy_bytes.to_string()),
            ("copy_gb_per_s", format!("{:.4}", a.copy_gb_per_s)),
            ("cycles", a.cycles.to_string()),
            ("ghz", format!("{:.4}", a.ghz)),
            (
                "cycles_per_triangle",
                format!("{:.2}", a.cycles_per_triangle),
            ),
            (
                "seam_duplicate_vertices",
                a.seam_duplicate_vertices.to_string(),
            ),
            (
                "unwelded_duplicate_vertices",
                a.unwelded_duplicate_vertices.to_string(),
            ),
            (
                "unwelded_boundary_edges",
                a.unwelded_boundary_edges.to_string(),
            ),
            ("boundary_edges", a.boundary_edges.to_string()),
            ("duplicate_vertices", a.duplicate_vertices.to_string()),
        ];
        match pair {
            Some(p) => {
                r.push(("c3_holds", p.held.to_string()));
                r.push(("c3_pair", p.pair_field.to_string()));
                r.push(("c3_pair_chunk_cells", p.pair_chunk_cells.to_string()));
                r.push(("triangle_mismatch", p.mismatch.to_string()));
                r.push((
                    "triangle_mismatch_share",
                    format!("{:.6}", p.mismatch_share),
                ));
                r.push(("c3_ratio", format!("{:.6}", p.ratio)));
                // The same ratio in cycles. `M-280`: on a governed and shared
                // machine the wall clock is the weaker instrument, and this row's
                // own `counted_over_median` says when it moved. C3's verdict is
                // read from `c3_ratio` because `cost_per_triangle` is the
                // registered column; `c3_ratio_cycles` is what makes the verdict
                // trustworthy, and both are in the file.
                r.push(("c3_ratio_cycles", format!("{:.6}", p.ratio_cycles)));
            }
            None => {
                r.push(("c3_holds", "NA".to_string()));
                r.push(("c3_pair", "NA".to_string()));
                r.push(("c3_pair_chunk_cells", "NA".to_string()));
                r.push(("triangle_mismatch", "NA".to_string()));
                r.push(("triangle_mismatch_share", "NA".to_string()));
                r.push(("c3_ratio", "NA".to_string()));
                r.push(("c3_ratio_cycles", "NA".to_string()));
            }
        }
        r
    }

    fn report(a: &Arm) {
        println!(
            "{:<12} {:>4}³  tri {:>7}  total {:>8.4} ms  |  copy {:>7.4} ({:>5.1}%)  \
             construct {:>7.4} ({:>5.1}%)  bvh {:>7.4} ({:>5.1}%)  handoff {:>7.4} ({:>5.1}%)",
            a.field,
            a.samples,
            a.triangles,
            a.total_ms,
            a.copy_ms,
            100.0 * a.copy_ms / a.total_ms,
            a.construct_ms,
            100.0 * a.construct_ms / a.total_ms,
            a.bvh_ms,
            100.0 * a.bvh_ms / a.total_ms,
            a.handoff_ms,
            100.0 * a.handoff_ms / a.total_ms,
        );
        println!(
            "             largest {} at {:.1}%, {}/{} reps agree; residual {:+.4} ms \
             ({:+.3}%), worst rep {:.3}%, without the weld {:.2}%; weld {:.4} / readiness {:.4}; \
             bvh is {:.1}% of the constructor; {:.4} GHz, {:.1} cycles/triangle; \
             seam edges {}, degenerate {}",
            a.largest_stage,
            100.0 * a.largest_share,
            a.reps_agreeing,
            REPS,
            a.residual_ms,
            100.0 * a.residual_share,
            100.0 * a.residual_share_worst_rep,
            100.0 * a.residual_share_without_weld,
            a.weld_ms,
            a.readiness_ms,
            100.0 * a.bvh_ms / a.trimesh_new_ms,
            a.ghz,
            a.cycles_per_triangle,
            a.seam_boundary_edges,
            a.degenerate_triangles,
        );
        println!(
            "             finest split: largest {} at {:.1}%; contour {:.4} ms, so the collider \
             build is {:.3}x the contour (Godot's docs say 3-5x)",
            a.finest_largest_stage,
            100.0 * a.finest_largest_share,
            a.contour_ms,
            a.total_ms / a.contour_ms,
        );
    }

    pub(crate) fn run(out: &mut Run) {
        let mut probe = Probe::open();

        let fbm = FbmTerrain::<Scalar>::canonical();
        let gyroid = capped_gyroid::<Scalar>();
        let sphere = Sphere::<Scalar>::canonical();

        // ── the untimed fixture phase ───────────────────────────────────────
        //
        // Everything that is not a clock happens here: contours, weld epsilons,
        // validator configs, and the search for C3's matched sphere resolution.
        // The first version of this harness searched for that resolution
        // *between* the gyroid arm and the sphere arm it is compared with, and
        // the consequence was measured rather than imagined: with five other
        // agents building on this host, the sphere arm's wall time moved **58%**
        // between two runs of the same binary while its cycle count moved 1.5%.
        println!("── the untimed fixture phase ─────────────────────────────────\n");
        let mut fixtures: Vec<Prepared> = Vec::new();
        let mut targets: Vec<(u32, usize, u32)> = Vec::new();
        for samples in RESOLUTIONS {
            let target = welded_triangles(&gyroid, samples);
            let at = sphere_samples_for(target);
            println!("gyroid {samples}³ welds to {target} triangles → sphere at {at}³");
            targets.push((samples, target, at));
        }
        for (samples, _, at) in &targets {
            fixtures.push(prepare(&fbm, "fbm_terrain", "c1c2", *samples));
            fixtures.push(prepare(&gyroid, "gyroid", "c1c2", *samples));
            fixtures.push(prepare(&sphere, "sphere", "c3_match", *at));
        }

        // ── the timed phase: every arm inside every repeat ───────────────────
        //
        // **The repeat is the outer loop and the arm is the inner one, and that
        // ordering is the whole point.** `M-281` says compare within one build
        // and one run; on a shared machine that is not enough, because "one run"
        // can still contain a two-second episode that lands entirely on one arm.
        // Measuring all six arms inside each of the nine repeats puts every
        // machine excursion across every arm, so a ratio between two arms
        // survives it even when neither absolute number does. The price is that
        // each build starts on a cache the other five arms have just evicted -
        // which is what a chunked pipeline building chunk after chunk actually
        // does, so it is the more representative state as well as the fairer one.
        println!("\n── the timed phase: {} arms interleaved inside every repeat ──\n", fixtures.len());
        let mut reps: Vec<Vec<Rep>> = fixtures.iter().map(|_| Vec::with_capacity(REPS)).collect();
        let mut last: Vec<Option<Stages>> = fixtures.iter().map(|_| None).collect();
        for rep in 0..(WARMUP + REPS) {
            for (i, p) in fixtures.iter().enumerate() {
                let (r, stages) = one_rep(p);
                if rep >= WARMUP {
                    reps[i].push(r);
                }
                last[i] = Some(stages);
            }
        }

        // One more wave, counted, inside the same interleaving so the counted
        // build sees the machine the repeats saw. `M-280`'s clock, and
        // `counted_over_median` is what says the two agree.
        let counted: Vec<Counted> = fixtures.iter().map(|p| count_once(p, &mut probe)).collect();

        let mut c1c2 = Vec::new();
        let mut pairs = Vec::new();
        for (i, p) in fixtures.iter().enumerate() {
            let stages = last[i].as_ref().expect("WARMUP + REPS is non-zero");
            let a = finish(p, &reps[i], stages, &counted[i]);
            report(&a);
            if a.arm == "c3_match" {
                let (samples, target, _) = targets
                    .iter()
                    .copied()
                    .find(|(_, _, at)| *at == a.samples)
                    .expect("every c3_match arm was prepared from a target");
                let share = (a.triangles as i64 - target as i64).abs() as f64
                    / target.max(a.triangles) as f64;
                assert!(
                    share <= MATCH_TOLERANCE,
                    "the closest sphere resolution ({}³, {} triangles) misses gyroid {samples}³'s \
                     {target} by {:.2}%, above the {:.0}% this fixture allows — C3 would be \
                     comparing two different triangle counts",
                    a.samples,
                    a.triangles,
                    share * 100.0,
                    MATCH_TOLERANCE * 100.0
                );
                pairs.push((samples, a));
            } else {
                if a.field == "gyroid" {
                    let (_, target, _) = targets
                        .iter()
                        .copied()
                        .find(|(s, _, _)| *s == a.samples)
                        .expect("every gyroid arm has a target");
                    assert_eq!(
                        a.triangles, target,
                        "the timed gyroid arm at {}³ welded to a different triangle count than the \
                         fixture phase did, so the sphere resolution was matched to a mesh that is \
                         not the one C3 compares",
                        a.samples
                    );
                }
                c1c2.push(a);
            }
        }

        // ── verdicts ────────────────────────────────────────────────────────
        //
        // C1 is over the registered arms only: the sphere rows are C3's control
        // and were never part of "33³ and 65³ chunks on fbm_terrain and gyroid".
        let c1 = c1c2.iter().all(|a| a.largest_share > DOMINANT_BAR);
        // C2 as registered: the prediction is the triangle copy. `bvh_dominates`
        // carries the registered falsifier separately, because a third stage
        // leading falsifies the prediction without vindicating the folklore.
        let c2 = c1c2.iter().all(|a| a.largest_stage == "copy");

        println!("\n── verdicts ─────────────────────────────────────────────────\n");
        for a in &c1c2 {
            println!(
                "C1 {} {}³: largest {} at {:.4} → {}",
                a.field,
                a.samples,
                a.largest_stage,
                a.largest_share,
                if a.largest_share > DOMINANT_BAR {
                    "over 50%"
                } else {
                    "under 50%"
                }
            );
        }
        println!(
            "C1 one stage over 50% on all four registered arms: {}",
            if c1 { "HELD" } else { "FALSIFIED" }
        );
        // C1's verdict is sensitive to exactly one fixture decision - that
        // "anything between" is one stage - and the sensitivity belongs in the
        // output rather than in the reader's head.
        //
        // **The aggregation is the MINIMUM over arms, and the first version of
        // this line took the maximum.** C1 reads "one stage is over 50% at 33³
        // *and* 65³ on fbm_terrain *and* gyroid", so it is a conjunction and its
        // weakest arm decides it - a max reports the best arm and calls the
        // clause held, which is the same one-sidedness `F-002` is about. The two
        // answers differ here: the finest split's shares are 0.481 to 0.501, so
        // the max says HELD and the min says FALSIFIED, and only the second is
        // what C1 asks.
        let finest_min = c1c2
            .iter()
            .map(|a| a.finest_largest_share)
            .fold(f64::INFINITY, f64::min);
        let finest_max = c1c2
            .iter()
            .map(|a| a.finest_largest_share)
            .fold(0.0f64, f64::max);
        println!(
            "C1 under the finest split this harness can make (handoff back into weld + gate), the \
             largest stage is {} on every arm and its share spans {:.4}-{:.4} → C1 would be {} on \
             that decomposition, decided by the weakest arm",
            c1c2[0].finest_largest_stage,
            finest_min,
            finest_max,
            if finest_min > DOMINANT_BAR {
                "HELD"
            } else {
                "FALSIFIED"
            }
        );
        println!(
            "C2 the dominant stage is the triangle copy: {} (bvh dominates on {} of {} arms)",
            if c2 { "HELD" } else { "FALSIFIED" },
            c1c2.iter().filter(|a| a.largest_stage == "bvh").count(),
            c1c2.len()
        );
        assert!(
            c1c2.iter()
                .all(|a| a.finest_largest_stage == c1c2[0].finest_largest_stage),
            "the finest split's largest stage is not the same on every arm, so the line above \
             names one arm's answer as if it were all of them"
        );

        // The registration's cited independent corroboration, checked inside one
        // build. Godot Voxel Tools: "creating a collider from a mesh is actually
        // much more expensive than meshing itself (about 3 to 5 times)", with no
        // absolute timings and no hardware.
        for a in c1c2.iter().chain(pairs.iter().map(|(_, p)| p)) {
            println!(
                "Godot's 3-5x, {} {}³: collider build {:.4} ms against contour {:.4} ms = {:.3}x{}",
                a.field,
                a.samples,
                a.total_ms,
                a.contour_ms,
                a.total_ms / a.contour_ms,
                if a.total_ms < a.contour_ms {
                    " — INVERTED, the contour is the larger cost"
                } else {
                    ""
                }
            );
        }

        // ── C3, computed once per pair in one direction ─────────────────────
        let mut pairings: Vec<(u32, Pairing, Pairing)> = Vec::new();
        for (samples, sphere_arm) in &pairs {
            let gyroid_arm = c1c2
                .iter()
                .find(|a| a.field == "gyroid" && a.samples == *samples)
                .expect("the gyroid arm at this resolution was just measured");
            let mismatch = gyroid_arm.triangles as i64 - sphere_arm.triangles as i64;
            let mismatch_share = mismatch.abs() as f64
                / gyroid_arm.triangles.max(sphere_arm.triangles) as f64;
            let ratio = gyroid_arm.cost_per_triangle / sphere_arm.cost_per_triangle;
            let ratio_cycles = gyroid_arm.cycles_per_triangle / sphere_arm.cycles_per_triangle;
            let held = ratio >= FIELD_PENALTY_BAR;
            println!(
                "C3 gyroid {samples}³ ({} tri, {:.6} µs/tri, {} seam edges, {} degenerate) against \
                 sphere {}³ ({} tri, {:.6} µs/tri, {} seam edges, {} degenerate): {ratio:.4}x → {}",
                gyroid_arm.triangles,
                gyroid_arm.cost_per_triangle,
                gyroid_arm.seam_boundary_edges,
                gyroid_arm.degenerate_triangles,
                sphere_arm.samples,
                sphere_arm.triangles,
                sphere_arm.cost_per_triangle,
                sphere_arm.seam_boundary_edges,
                sphere_arm.degenerate_triangles,
                if held { "HELD" } else { "FALSIFIED" }
            );
            println!(
                "   the same ratio in cycles, which no governor and no neighbour moves: \
                 {:.1} against {:.1} cycles/triangle = {ratio_cycles:.4}x; the two arms' counted \
                 builds tracked their medians at {:.3} and {:.3}; triangle counts differ by \
                 {mismatch} ({:.3}%)",
                gyroid_arm.cycles_per_triangle,
                sphere_arm.cycles_per_triangle,
                gyroid_arm.counted_total_ms / gyroid_arm.total_ms,
                sphere_arm.counted_total_ms / sphere_arm.total_ms,
                mismatch_share * 100.0,
            );
            // Two views of one comparison: each row names the *other* arm, and
            // every number is gyroid ÷ sphere on both.
            pairings.push((
                *samples,
                Pairing {
                    pair_field: sphere_arm.field,
                    pair_chunk_cells: sphere_arm.samples,
                    mismatch,
                    mismatch_share,
                    ratio,
                    ratio_cycles,
                    held,
                },
                Pairing {
                    pair_field: gyroid_arm.field,
                    pair_chunk_cells: gyroid_arm.samples,
                    mismatch,
                    mismatch_share,
                    ratio,
                    ratio_cycles,
                    held,
                },
            ));
        }
        let c3 = pairings.iter().all(|(_, g, _)| g.held);
        println!(
            "C3 gyroid costs at least 1.5x sphere per triangle: {}",
            if c3 { "HELD" } else { "FALSIFIED" }
        );

        // Rows in fixture order - fbm, gyroid, sphere at 33³ then the same at
        // 65³ - so the file reads the way the run does.
        let mut rows: Vec<Row> = Vec::new();
        for a in c1c2.iter().chain(pairs.iter().map(|(_, p)| p)) {
            let pairing = match a.field {
                // The gyroid row's view: its pair is the sphere.
                "gyroid" => pairings
                    .iter()
                    .find(|(s, _, _)| *s == a.samples)
                    .map(|(_, g, _)| g),
                // The sphere row's view: its pair is the gyroid, located by the
                // sphere resolution the gyroid-side view names.
                "sphere" => pairings
                    .iter()
                    .find(|(_, g, _)| g.pair_chunk_cells == a.samples)
                    .map(|(_, _, s)| s),
                _ => None,
            };
            rows.push(row(
                a,
                if a.arm == "c1c2" {
                    c1
                } else {
                    a.largest_share > DOMINANT_BAR
                },
                if a.arm == "c1c2" {
                    c2
                } else {
                    a.largest_stage == "copy"
                },
                c3,
                pairing,
            ));
        }

        for r in &rows {
            out.record(r);
        }
    }
}

fn main() {
    if !std::env::args().any(|arg| arg == "--bench") {
        return;
    }
    let prereg = isomesh::experiment!("P-85");

    #[cfg(target_os = "linux")]
    common::experiment::run(prereg, experiment::run);

    #[cfg(not(target_os = "linux"))]
    {
        eprintln!(
            "{} puts a clock on every row — `ghz` and `cycles_per_triangle`, which come\n\
             from hardware performance counters — because M-280 established that a\n\
             nanosecond on a governed CPU is not a unit. This platform has no\n\
             `perf_event_open`. Refusing rather than reporting stage shares with no\n\
             clock behind them. Run it on Linux; `perf_event_paranoid = 2` is\n\
             permissive enough and no root is needed.",
            prereg.id
        );
        std::process::exit(1);
    }
}
