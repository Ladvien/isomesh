//! **P-121 — what fraction of extraction is bit work — runs first, gates Group A.**
//!
//! Ticket: R-121. Pre-registered before this harness existed.
//!
//! ```bash
//! cargo bench --bench experiment_p121
//! ```
//!
//! Writes `docs/experiments/p-121.csv`. **Linux only**, for `experiment_p12`'s
//! reason: the counters are `perf_event_open`, and a verdict about integer work
//! cannot be carried by a clock (`✗24`, `M-281`).
//!
//! # What was missing
//!
//! `✗51` is the incident this row exists to prevent. `P-69`'s C1 asked for 2× on
//! the marginal extraction cost from restructuring the **sample** loop, and the
//! loop is **11.6%** of that marginal — so the clause's ceiling was
//! `1/(1 − 0.116/2) = 1.061×` and it was unreachable however well the loop
//! vectorised. Nobody had measured the share first. Group A (`P-103`, `P-104`,
//! `P-106`) is about to propose three bit-parallel mechanisms for the
//! *classification* stage, and there was no number anywhere in the repository
//! saying how much of extraction that stage is.
//!
//! There was no *instrument* either, and that is the harder gap. `M-135`'s
//! `docs/measurements/stage_breakdown.csv` decomposes the **pipeline**
//! (`contour / normals / weld / collider`), not the extractor. `experiment_p15`
//! decomposes `DualMesher` into `sample / place / emit_quads` by *regression*
//! over grid shapes, which cannot see a classify stage at all, because
//! classification and placement share one loop. `M-279` reports whole-extraction
//! instructions per sample and no split whatever. So this harness had to build
//! the decomposition, and `crates/isomesh/src/**` is read-only for the phase.
//!
//! # The mirror, and why its first job is to agree
//!
//! The five stages do not exist as five loops in the crate: `marching_cubes`
//! computes the case index, reads the table, interpolates the crossing and
//! pushes the triangle **inside one cell body** (`mod.rs:254-380`), and `dual`
//! fuses classification into `place_vertices` (`dual.rs:487-547`). Counter
//! windows around five things that are one thing is not possible.
//!
//! So the pipeline is rebuilt **bench-local**, one pass per stage, the way
//! `experiment_p40` rebuilds `DualMesher`'s layout in its own `Grid`. Every
//! formula comes from the crate — `table::{CASES, EDGE_AXIS, EDGE_CORNERS,
//! edge_offset, is_inside}`, `HermiteCell::from_corners` and `solve::solve_with`
//! are all public and are *called*, not re-derived. Only `corner_offset`,
//! `corner_position`, `place`, `vec3::length`, `vec3::scale`, `unit_gradient`
//! and the `ToCell` clamp are private, and those are copied verbatim rather than
//! made `pub`.
//!
//! A mirror is worth nothing on its own, so **every row asserts that the mirror
//! emits the shipped extractor's mesh bit for bit** — positions, normals and
//! indices, compared as bit patterns, in creation order. That is `M-279`'s rule:
//! a new instrument's first job is to agree with the old one where they overlap,
//! and it is what licenses reading a share off the mirror at all.
//! `mesh_identical_to_shipped` is a column and it is asserted;
//! `agreement_ratio` reports the mirror's cycles against the shipped
//! extractor's, on the same grid, in the same run.
//!
//! # Where the stage boundaries fall, and why
//!
//! The registration fixes the six names. It does not, by itself, say which side
//! of the boundary the **cell walk** falls on, and that decision moves both C2
//! and C3, so it is made here and stated:
//!
//! - **`classify` owns the walk over every cell.** The eight-sign case index is
//!   what the walk is *for*, it is computed for all `(n−1)³` cells whether or not
//!   the surface is there, and the compaction that turns "this cell is mixed"
//!   into a list is the same integer pass. For the dual path `classify` also owns
//!   the bit-packing (`dual.rs:359-381`), the fused `any & !all` word test
//!   (`:424`) and the crossed-grid-edge scan that decides *which quads exist*
//!   (`:697-713`) — all three are sign tests, and none of them pushes a triangle.
//! - **`emit` owns only what happens on a cell that produces something**: the
//!   table lookup, the edge-cache probe and index assignment, the triangle and
//!   quad push. It iterates the compacted list, so on a field with no surface it
//!   iterates nothing.
//!
//! Putting the walk in `emit` instead would have made C3 false by construction —
//! a per-cell `continue` costs the same whether or not there is a surface — and
//! it would have credited classification with none of the work it does.
//!
//! `emit` is measured as **two prefix steps** whose sum is the registered
//! `cycles_emit`: `cycles_emit_prepare`, the surface-independent buffer
//! initialisation (`marching_cubes/mod.rs:250-251` resizes `edge_vertices` to
//! `3·sample_count` slots of `u32::MAX` on *every* extract — 3.3 MB at 65³, a
//! cost this instrument names for the first time), and `cycles_emit_walk`, the
//! productive remainder. C3 is scored on `emit_walk`, because a memset that runs
//! either way cannot testify about separability.
//!
//! # How a stage is measured: prefixes, not isolation
//!
//! **A stage is never timed on its own.** Isolation changes what a stage costs,
//! and three instruments were built before this one to establish that rather
//! than assume it:
//!
//! - *Stage batches* — `sample × inner`, then `classify × inner` — run with the
//!   stage's own arrays already resident. At 33³ the mirror's working set is
//!   about 600 KB against a 512 KB L2, so the batched decomposition came out
//!   **6% cheaper** than the extraction it was meant to decompose.
//! - *One window per stage per pass* puts a kernel round trip inside every
//!   stage's counted region. That over-attributed **+4.5% to +6.4%** on the
//!   cheapest fields (`sphere` and `sphere_surface_free` under Marching Cubes,
//!   where a whole extraction is 26 cycles per cell) and lifted the tiny
//!   surface-free `emit_walk` to 5.5% of `sphere`'s, which is instrument noise
//!   wearing C3's clothes.
//! - Measured on its own, `sample` at `sphere 33³ f32` reads 5.19 cycles per
//!   sample in one shape and 6.74 in the other. Both cannot be its share.
//!
//! So the instrument is [`Mc::prefix`] and [`Dc::prefix`]: `cut[k]` is a counter
//! window over the first `k` stages of a pipeline pass, batched over
//! `inner_reps` passes, and a stage is `cut[k] − cut[k−1]`. Every prefix runs the
//! stages in pipeline order with the pipeline's own cache and predictor state,
//! every prefix window has exactly one counter boundary — amortised over
//! `inner_reps` passes and cancelling in the difference — and `cut[stages]` is
//! the whole pipeline. Because the stages telescope, the decomposition partitions
//! the extraction **by construction**.
//!
//! Medians are taken **per quantity** rather than per repetition: the cuts are
//! monotone in `k`, and one repetition disturbed by another process on the
//! machine moves one cut and therefore two stages. Medianing each cut over
//! [`REPS`] repetitions keeps the telescoping exact and makes no stage hostage
//! to one reading.
//!
//! # Residual: what it is, and why it is an absolute value
//!
//! `cycles_total` is a **separate** pair of counter windows over whole pipeline
//! passes — one before the prefix sweep and one after, averaged so a
//! monotonically drifting clock cancels — and `cycles_residual` is
//! `cycles_total − Σ five stages`, which is `cycles_total − cut[stages]`. The two
//! measure the *same function body* in different windows, so the residual is the
//! instrument's own reproducibility, and C1 is the clause that says the
//! accounting means something. That is weaker than it would be if the stages
//! could be measured independently without changing them — they cannot, and the
//! three attempts above are why — and it is stated here rather than implied.
//! The other half of C1, *"any stage reading zero on a fixture where it must
//! run"*, is asserted per row and not merely recorded.
//!
//! Being a difference of two measurements, the residual can come out
//! **negative**. A signed `residual_share < 0.05` would then pass on a 3%
//! *overshoot*, which is exactly as much a failure to account for the total as a
//! 3% shortfall. So `residual_share` is `|residual| / total` and the signed form
//! is beside it as `residual_signed_share`.
//!
//! `c1_holds` is scored on **the row's** `residual_share`, which is the
//! registration's wording — *"under 5% ... on EVERY row"*. The worst single
//! repetition of any row is reported beside it as `residual_share_rep_worst`,
//! because a 30 ms window on a machine running other work swings much further
//! than the clause's bar and hiding that would be worse than reporting it.
//!
//! Every window is batched over `inner_reps` repetitions, so the ~28
//! `perf_event` system calls a window costs land outside it and cannot inflate
//! anything. `inner_reps` is chosen per row to make one window about
//! [`TARGET_BATCH_NS`] nanoseconds, and it is a column.
//!
//! # SHARE
//!
//! **This row is the share instrument.** It moves nothing and claims nothing;
//! it produces the denominator the rest of Group A is denominated in. Its own
//! clauses therefore have no `1/(1 − share/factor)` ceiling, and each one's
//! reachable share is a column rather than an argument:
//!
//! - **C1's share is `residual_share`, and it is the whole clause.** The bar is
//!   5% of `cycles_total`, per row. Nothing else in the file can be believed if
//!   this fails, which is why it is clause one.
//! - **C2's share is `integer_share` = `(cycles_classify + cycles_emit) /
//!   cycles_total`**, over the eight *reference* fields at 65³.
//!   `sphere_surface_free` is a control, not a reference field, and is excluded
//!   from C2's maximum. The bar is 0.15 on at least one of them.
//!   `max_reference_integer_share_at_65` is a column, so the closure decision
//!   for `P-103`, `P-104` and `P-106` can be re-checked from one cell of the CSV
//!   instead of re-derived from seventy-two rows.
//! - **C3's share is the surface-free arm against its own `sphere` arm** at the
//!   same resolution, scalar and extractor. C3 holds for such a group when all
//!   of these hold:
//!   1. `active_cells == 0`, `vertices == 0` and `triangles == 0` on the
//!      surface-free arm — the emit stage genuinely produced nothing;
//!   2. `cycles_classify > 0` there, and `c3_classify_hold` — its classify cycles
//!      over `sphere`'s — is at least [`CLASSIFY_HOLD_FLOOR`], because
//!      classification is per-cell, surface-independent, and must **not**
//!      collapse;
//!   3. `c3_emit_walk_collapse` — surface-free `cycles_emit_walk` over
//!      `sphere`'s — is below [`EMIT_COLLAPSE_CEILING`], **and**
//!      `c3_emit_walk_collapse_instructions` is too. The second is the
//!      load-bearing one: on the surface-free arm `emit_walk` is the difference
//!      of two nearly equal windows, so its *cycle* ratio has an unbounded
//!      relative noise floor, while an instruction count is deterministic and
//!      cache-independent (`M-279`, `experiment_p15`). Requiring both can only
//!      make the clause harder;
//!   4. `cycles_sample` and `cycles_classify` are non-zero on **every** row,
//!      which is the registration's own vacuity control and is asserted rather
//!      than merely recorded.
//!
//! `float_share` is `(sample + interpolate + solve) / total` and
//! `integer_share` is `(classify + emit) / total`, exactly as registered.
//!
//! # `ns_per_cell` and `ghz` are provenance, not verdicts
//!
//! `M-280` and `M-281`: on a governed CPU a nanosecond is not a unit. Every
//! clause here reads cycles or instructions. `ns_per_cell` and `ghz` are on the
//! row so a later reader can see what clock the cycles were taken at, and no
//! clause consults either.

#![allow(
    clippy::cast_possible_truncation,
    clippy::cast_precision_loss,
    clippy::float_cmp,
    clippy::too_many_arguments,
    clippy::too_many_lines
)]

mod common;

#[cfg(target_os = "linux")]
mod experiment {
    use std::hint::black_box;
    use std::time::Instant;

    use isomesh::dual_contouring::{CLAMP_EPSILON, DualContouring, solve};
    use isomesh::fields::Sphere;
    use isomesh::hermite::HermiteCell;
    use isomesh::marching_cubes::MarchingCubes;
    use isomesh::marching_cubes::table::{
        CASES, CENTROID_BASE, EDGE_AXIS, EDGE_CORNERS, MAX_CENTROIDS, edge_offset, is_centroid,
        is_inside,
    };
    use isomesh::{MeshBuffer, Real, RuntimeShape3, Sdf};

    use crate::common::counters::{MIN_TIME_RATIO, Probe};

    // ─── the fixture ───────────────────────────────────────────────────────

    /// The registered resolutions, in samples per axis.
    const RESOLUTIONS: [u32; 2] = [33, 65];
    /// Measured repetitions per row. The median by `cycles_total` is reported;
    /// the **worst** residual over all of them is what C1 is scored on.
    const REPS: usize = 9;
    /// Untimed passes before anything is counted, so a batch is steady state.
    const WARMUP: usize = 2;
    /// How long one counter window should last, in nanoseconds.
    ///
    /// A window costs about 28 `perf_event` system calls, all of them *outside*
    /// the counted region, so this does not exist to amortise overhead — it
    /// exists so that two independently taken windows over the same work differ
    /// by noise rather than by clock drift, which is what the residual measures.
    const TARGET_BATCH_NS: f64 = 30_000_000.0;
    /// Ceiling on the batch, so a cheap row cannot take minutes.
    const MAX_INNER: usize = 4096;
    /// The registered surface-free control: `experiment_p40.rs:84-93`'s field —
    /// the canonical sphere sampled a long way from itself, so no cell is
    /// active. `M-337` and `✗28` both scored on it.
    const SURFACE_FREE: &str = "sphere_surface_free";
    /// The eight reference fields C2's maximum is taken over.
    const REFERENCE_FIELDS: usize = 8;

    /// The most classify may fall by on the surface-free arm. C3, fixed here.
    const CLASSIFY_HOLD_FLOOR: f64 = 0.5;
    /// The most `emit_walk` may survive on the surface-free arm. C3, fixed here.
    const EMIT_COLLAPSE_CEILING: f64 = 0.05;
    /// C1's bar, from the registration.
    const RESIDUAL_CEILING: f64 = 0.05;
    /// C2's bar, from the registration.
    const INTEGER_SHARE_BAR: f64 = 0.15;

    // ─── private crate mechanisms, copied rather than made `pub` ────────────

    /// `cube::corner_offset`. Private, and `src/**` is read-only this phase.
    #[inline]
    const fn corner_offset(corner: u8) -> [u32; 3] {
        [
            (corner & 1) as u32,
            ((corner >> 1) & 1) as u32,
            ((corner >> 2) & 1) as u32,
        ]
    }

    /// `cube::place`: the centred frame, spelled exactly once in the crate and
    /// exactly once here.
    #[inline]
    fn place<R: Real>(lo: R, hi: R, d: R) -> R {
        (lo + hi) * R::HALF + (hi - lo) * d
    }

    /// `vec3::length`. Left-to-right summation, which is what makes the mirror's
    /// normals bit-identical to the crate's rather than merely close.
    #[inline]
    fn length<R: Real>(a: [R; 3]) -> R {
        (a[0] * a[0] + a[1] * a[1] + a[2] * a[2]).sqrt()
    }

    /// `marching_cubes::unit_gradient`, over `vec3::scale`.
    #[inline]
    fn unit_gradient<R: Real, S: Sdf<Scalar = R>>(sdf: &S, position: [R; 3]) -> [R; 3] {
        let g = sdf.gradient(position);
        let inv = length(g).recip();
        [g[0] * inv, g[1] * inv, g[2] * inv]
    }

    /// `marching_cubes::corner_position`.
    #[inline]
    fn corner_position<R: Real>(
        base: [u32; 3],
        corner: u8,
        origin: [R; 3],
        cell_size: R,
    ) -> [R; 3] {
        let o = corner_offset(corner);
        [
            origin[0] + cell_size * R::from_f64(f64::from(base[0] + o[0])),
            origin[1] + cell_size * R::from_f64(f64::from(base[1] + o[1])),
            origin[2] + cell_size * R::from_f64(f64::from(base[2] + o[2])),
        ]
    }

    /// `dual_contouring::apply_clamp` under `Clamp::ToCell`, which is
    /// `DualContouring::new`'s default and therefore the shipped path.
    #[inline]
    fn clamp_to_cell<R: Real>(x: [R; 3], cell_origin: [R; 3], cell_size: R) -> [R; 3] {
        let half = cell_size * R::HALF;
        let inset = half * R::from_f64(1.0 - CLAMP_EPSILON);
        let mut out = x;
        for (axis, slot) in out.iter_mut().enumerate() {
            let centre = cell_origin[axis] + half;
            *slot = slot.clamp(centre - inset, centre + inset);
        }
        out
    }

    // ─── the grid ──────────────────────────────────────────────────────────

    /// One cubic grid: `n` samples per axis, spanning `n − 1` cells.
    #[derive(Clone, Copy)]
    struct Grid<R: Real> {
        n: u32,
        origin: [R; 3],
        cell_size: R,
    }

    impl<R: Real> Grid<R> {
        #[inline]
        fn samples(self) -> usize {
            let n = self.n as usize;
            n * n * n
        }

        #[inline]
        fn cells(self) -> u32 {
            self.n - 1
        }

        #[inline]
        fn cell_count(self) -> usize {
            let c = self.cells() as usize;
            c * c * c
        }

        /// `RuntimeShape3::linearize`, which is the layout Marching Cubes
        /// samples into: no row padding, so the stride is the row itself.
        #[inline]
        fn sample_index(self, p: [u32; 3]) -> usize {
            let n = self.n as usize;
            p[0] as usize + n * (p[1] as usize + n * p[2] as usize)
        }

        /// `dual::row_stride`: `size[0] | 1`, A-024's odd row.
        #[inline]
        fn dual_row(self) -> usize {
            self.n as usize | 1
        }

        /// `DualMesher::index`, which is **not** `linearize` — one bit of
        /// difference and 3.37× at 128³ (`M-287`).
        #[inline]
        fn dual_index(self, p: [u32; 3]) -> usize {
            p[0] as usize + self.dual_row() * (p[1] as usize + self.n as usize * p[2] as usize)
        }

        #[inline]
        fn cell_index(self, p: [u32; 3]) -> usize {
            let c = self.cells() as usize;
            p[0] as usize + c * (p[1] as usize + c * p[2] as usize)
        }

        #[inline]
        fn cell_origin(self, base: [u32; 3]) -> [R; 3] {
            [
                self.origin[0] + self.cell_size * R::from_f64(f64::from(base[0])),
                self.origin[1] + self.cell_size * R::from_f64(f64::from(base[1])),
                self.origin[2] + self.cell_size * R::from_f64(f64::from(base[2])),
            ]
        }
    }

    /// One cell the surface passes through, as its grid base.
    ///
    /// `case` is Marching Cubes' eight-sign index; the dual path leaves it zero,
    /// because its classification produces the *word* answer and never needs the
    /// per-cell byte. Neither path stores the linear cell index, because neither
    /// shipped path re-reads one — the dual's `cell_first[index] = slot` write
    /// consumes it where it is computed.
    #[derive(Clone, Copy)]
    struct Cell {
        base: [u32; 3],
        case: u8,
    }

    // ─── the Marching Cubes mirror ─────────────────────────────────────────

    /// One vertex the emit walk allocated, and what `interpolate` must compute
    /// for it: an edge crossing when `code < CENTROID_BASE`, otherwise cycle
    /// centroid `code − CENTROID_BASE` of `case`.
    #[derive(Clone, Copy)]
    struct Job {
        base: [u32; 3],
        code: u8,
        case: u8,
    }

    /// `MarchingCubes`, one pass per stage.
    #[derive(Default)]
    struct Mc<R: Real> {
        values: Vec<R>,
        active: Vec<Cell>,
        /// `MarchingCubes::edge_vertices`: one `u32` slot per (sample, axis).
        edge_vertices: Vec<u32>,
        jobs: Vec<Job>,
        positions: Vec<[R; 3]>,
        normals: Vec<[R; 3]>,
        indices: Vec<u32>,
    }

    impl<R: Real> Mc<R> {
        /// `sdf::sample_grid` with `row_stride == size[0]`, which is what
        /// `MarchingCubes::extract` passes.
        fn sample<S: Sdf<Scalar = R>>(&mut self, sdf: &S, g: Grid<R>) {
            self.values.clear();
            self.values.reserve(g.samples());
            for z in 0..g.n {
                for y in 0..g.n {
                    for x in 0..g.n {
                        self.values.push(sdf.sample([
                            g.origin[0] + g.cell_size * R::from_f64(f64::from(x)),
                            g.origin[1] + g.cell_size * R::from_f64(f64::from(y)),
                            g.origin[2] + g.cell_size * R::from_f64(f64::from(z)),
                        ]));
                    }
                }
            }
            black_box(&self.values);
        }

        /// The eight-sign case index for every cell, and the compaction.
        ///
        /// `case == 0` and `case == 255` are exactly the cells the table gives
        /// no triangles for — `validate_table`'s `cut_edges_without_triangles`
        /// is zero over all 256 entries — so dropping them here is the same
        /// decision `extract`'s `entry.count == 0 { continue }` makes, taken
        /// without reading the table.
        fn classify(&mut self, g: Grid<R>) {
            self.active.clear();
            let c = g.cells();
            for z in 0..c {
                for y in 0..c {
                    for x in 0..c {
                        let base = [x, y, z];
                        let mut case = 0u8;
                        for corner in 0..8u8 {
                            let o = corner_offset(corner);
                            let s =
                                g.sample_index([base[0] + o[0], base[1] + o[1], base[2] + o[2]]);
                            if is_inside(self.values[s]) {
                                case |= 1 << corner;
                            }
                        }
                        if case != 0 && case != u8::MAX {
                            self.active.push(Cell { base, case });
                        }
                    }
                }
            }
            black_box(&self.active);
        }

        /// `mod.rs:250-251`: the edge cache, cleared and resized every extract.
        fn emit_prepare(&mut self, g: Grid<R>) {
            self.edge_vertices.clear();
            self.edge_vertices.resize(g.samples() * 3, u32::MAX);
            self.jobs.clear();
            self.indices.clear();
            black_box(&self.edge_vertices);
        }

        /// The table lookup, the edge-cache indexing and the triangle push.
        ///
        /// Vertex indices are allocated in exactly `extract`'s order — a cell's
        /// cycle centroids first, then each triangle's three codes in order,
        /// creating on a cache miss — so `jobs[i]` is vertex `i` and the mesh
        /// comes out bit-identical.
        fn emit_walk(&mut self, g: Grid<R>) {
            let mut next = 0u32;
            for cell in &self.active {
                let entry = &CASES[cell.case as usize];
                let mut centroid = [0u32; MAX_CENTROIDS];
                for (k, slot) in centroid
                    .iter_mut()
                    .enumerate()
                    .take(entry.centroids as usize)
                {
                    *slot = next;
                    next += 1;
                    self.jobs.push(Job {
                        base: cell.base,
                        code: CENTROID_BASE + k as u8,
                        case: cell.case,
                    });
                }
                for tri in &entry.triangles[..entry.count as usize] {
                    for &code in tri {
                        let index = if is_centroid(code) {
                            centroid[(code - CENTROID_BASE) as usize]
                        } else {
                            let axis = EDGE_AXIS[code as usize] as usize;
                            let lo = EDGE_CORNERS[code as usize][0];
                            let o = corner_offset(lo);
                            let lo_sample = g.sample_index([
                                cell.base[0] + o[0],
                                cell.base[1] + o[1],
                                cell.base[2] + o[2],
                            ]);
                            let key = lo_sample * 3 + axis;
                            let cached = self.edge_vertices[key];
                            if cached == u32::MAX {
                                let fresh = next;
                                next += 1;
                                self.edge_vertices[key] = fresh;
                                self.jobs.push(Job {
                                    base: cell.base,
                                    code,
                                    case: cell.case,
                                });
                                fresh
                            } else {
                                cached
                            }
                        };
                        self.indices.push(index);
                    }
                }
            }
            black_box(&self.indices);
        }

        /// Crossing positions and normals, in vertex-index order.
        fn interpolate<S: Sdf<Scalar = R>>(&mut self, sdf: &S, g: Grid<R>) {
            self.positions.clear();
            self.normals.clear();
            self.positions.reserve(self.jobs.len());
            self.normals.reserve(self.jobs.len());
            for job in &self.jobs {
                let position = if is_centroid(job.code) {
                    // `mod.rs:320-355`: a cycle's edges are the non-centroid
                    // corners of the triangles naming it, summed in table order.
                    let entry = &CASES[job.case as usize];
                    let mut sum = [R::ZERO; 3];
                    let mut n = 0u32;
                    for tri in &entry.triangles[..entry.count as usize] {
                        if tri[0] != job.code {
                            continue;
                        }
                        let p = crossing(&self.values, g, job.base, tri[1]);
                        sum = [sum[0] + p[0], sum[1] + p[1], sum[2] + p[2]];
                        n += 1;
                    }
                    let scale = R::from_f64(f64::from(n)).recip();
                    [sum[0] * scale, sum[1] * scale, sum[2] * scale]
                } else {
                    crossing(&self.values, g, job.base, job.code)
                };
                self.positions.push(position);
                self.normals.push(unit_gradient(sdf, position));
            }
            black_box(&self.positions);
        }

        /// Stages, in pipeline order.
        const STAGES: usize = 5;

        /// The first `upto` stages of one pipeline pass.
        ///
        /// **The whole instrument is prefixes of this function.** A stage is
        /// measured as the difference between two prefix windows, never in
        /// isolation, because isolation changes what the stage costs: measured
        /// on its own, `sample` at `sphere 33³ f32` reads 5.19 cycles per sample
        /// against 6.74 in a batch of its own and the two cannot both be its
        /// share of an extraction. A prefix window runs the stages in pipeline
        /// order with the pipeline's cache and predictor state, has exactly one
        /// counter boundary — amortised over `inner` passes — and the boundary
        /// cancels in the difference. `upto == STAGES` is the whole pipeline, and
        /// it is the same function body `cycles_total` is measured over, which is
        /// what makes the residual a measure of the instrument's own
        /// reproducibility rather than of an artefact it introduced.
        fn prefix<S: Sdf<Scalar = R>>(&mut self, sdf: &S, g: Grid<R>, upto: usize) {
            if upto >= 1 {
                self.sample(sdf, g);
            }
            if upto >= 2 {
                self.classify(g);
            }
            if upto >= 3 {
                self.emit_prepare(g);
            }
            if upto >= 4 {
                self.emit_walk(g);
            }
            if upto >= 5 {
                self.interpolate(sdf, g);
            }
        }

        fn pipeline<S: Sdf<Scalar = R>>(&mut self, sdf: &S, g: Grid<R>) {
            self.prefix(sdf, g, Self::STAGES);
        }
    }

    /// `marching_cubes::edge_position` with `crossing_refinement == 0`, which is
    /// `MarchingCubes::new`'s default: `refine_crossing` returns `d0` unchanged
    /// at zero steps, so there is no field evaluation on this path.
    #[inline]
    fn crossing<R: Real>(values: &[R], g: Grid<R>, base: [u32; 3], edge: u8) -> [R; 3] {
        let [lo_corner, hi_corner] = EDGE_CORNERS[edge as usize];
        let corner = |c: u8| {
            let o = corner_offset(c);
            values[g.sample_index([base[0] + o[0], base[1] + o[1], base[2] + o[2]])]
        };
        let d = edge_offset(corner(lo_corner), corner(hi_corner));
        let lo_pos = corner_position(base, lo_corner, g.origin, g.cell_size);
        let hi_pos = corner_position(base, hi_corner, g.origin, g.cell_size);
        [
            place(lo_pos[0], hi_pos[0], d),
            place(lo_pos[1], hi_pos[1], d),
            place(lo_pos[2], hi_pos[2], d),
        ]
    }

    // ─── the Dual Contouring mirror ────────────────────────────────────────

    /// One crossed grid edge: the quad the emit stage will push, found by the
    /// sign scan in `classify`.
    #[derive(Clone, Copy)]
    struct Quad {
        p: [u32; 3],
        axis: u8,
        inside0: bool,
    }

    /// `DualContouring` with the default `Qef` rule, one pass per stage.
    #[derive(Default)]
    struct Dc<R: Real> {
        /// `DualMesher::values`, on the **odd** row stride.
        values: Vec<R>,
        /// `DualMesher::inside`, one bit per sample along `x`.
        inside: Vec<u64>,
        bit_row: usize,
        active: Vec<Cell>,
        /// `DualMesher::cell_first`, doubling as the active flag.
        cell_first: Vec<u32>,
        quads: Vec<Quad>,
        hermite: Vec<HermiteCell<R>>,
        positions: Vec<[R; 3]>,
        normals: Vec<[R; 3]>,
        indices: Vec<u32>,
    }

    impl<R: Real> Dc<R> {
        fn sample<S: Sdf<Scalar = R>>(&mut self, sdf: &S, g: Grid<R>) {
            let row = g.dual_row();
            let pad = row - g.n as usize;
            self.values.clear();
            self.values.reserve(row * g.n as usize * g.n as usize);
            for z in 0..g.n {
                for y in 0..g.n {
                    for x in 0..g.n {
                        self.values.push(sdf.sample([
                            g.origin[0] + g.cell_size * R::from_f64(f64::from(x)),
                            g.origin[1] + g.cell_size * R::from_f64(f64::from(y)),
                            g.origin[2] + g.cell_size * R::from_f64(f64::from(z)),
                        ]));
                    }
                    for _ in 0..pad {
                        self.values.push(R::ZERO);
                    }
                }
            }
            black_box(&self.values);
        }

        #[inline]
        fn inside_word(&self, w: usize, y: usize, z: usize, n: usize) -> u64 {
            self.inside[self.bit_row * (y + n * z) + w]
        }

        #[inline]
        fn inside_word_shifted(&self, w: usize, y: usize, z: usize, n: usize) -> u64 {
            let lo = self.inside_word(w, y, z, n);
            let hi = if w + 1 < self.bit_row {
                self.inside_word(w + 1, y, z, n)
            } else {
                0
            };
            (lo >> 1) | (hi << 63)
        }

        #[inline]
        fn active_word(&self, w: usize, y: usize, z: usize, n: usize) -> u64 {
            let mut any = 0u64;
            let mut all = !0u64;
            for dz in 0..2usize {
                for dy in 0..2usize {
                    let a = self.inside_word(w, y + dy, z + dz, n);
                    let b = self.inside_word_shifted(w, y + dy, z + dz, n);
                    any |= a | b;
                    all &= a & b;
                }
            }
            any & !all
        }

        #[inline]
        fn cell_mask(w: usize, cells_x: usize) -> u64 {
            let remaining = cells_x.saturating_sub(w * 64);
            if remaining >= 64 {
                !0
            } else {
                (1u64 << remaining) - 1
            }
        }

        /// Every sign test the dual path makes: the bit-packing, the fused
        /// `any & !all` active-cell word, the set-bit walk that fixes vertex
        /// order, and the crossed-grid-edge scan that decides which quads exist.
        fn classify(&mut self, g: Grid<R>) {
            let n = g.n as usize;
            let row = g.dual_row();

            // `dual.rs:359-381`.
            self.bit_row = n.div_ceil(64);
            self.inside.clear();
            self.inside.resize(self.bit_row * n * n, 0);
            for r in 0..n * n {
                let src = row * r;
                let dst = self.bit_row * r;
                for w in 0..self.bit_row {
                    let base = w * 64;
                    let take = (n - base).min(64);
                    let mut word = 0u64;
                    for k in 0..take {
                        word |= u64::from(is_inside(self.values[src + base + k])) << k;
                    }
                    self.inside[dst + w] = word;
                }
            }

            // `dual.rs:487-497`, in ascending `x`, so vertex order is the
            // scalar loop's.
            self.cell_first.clear();
            self.cell_first.resize(g.cell_count(), u32::MAX);
            self.active.clear();
            let cells = g.cells() as usize;
            let cell_words = cells.div_ceil(64);
            let mut slot = 0u32;
            for z in 0..cells {
                for y in 0..cells {
                    for w in 0..cell_words {
                        let mut bits = self.active_word(w, y, z, n) & Self::cell_mask(w, cells);
                        while bits != 0 {
                            let x = w * 64 + bits.trailing_zeros() as usize;
                            bits &= bits - 1;
                            let base = [x as u32, y as u32, z as u32];
                            let index = g.cell_index(base);
                            self.cell_first[index] = slot;
                            slot += 1;
                            self.active.push(Cell { base, case: 0 });
                        }
                    }
                }
            }

            // `dual.rs:697-713`, three axes in order, so quad order is the
            // shipped order.
            self.quads.clear();
            for axis in 0..3usize {
                let u = (axis + 1) % 3;
                let v = (axis + 2) % 3;
                let mut p = [0u32; 3];
                for a in 0..g.n - 1 {
                    for b in 1..g.cells() {
                        for c in 1..g.cells() {
                            p[axis] = a;
                            p[u] = b;
                            p[v] = c;
                            let s0 = g.dual_index(p);
                            let mut q = p;
                            q[axis] += 1;
                            let s1 = g.dual_index(q);
                            let inside0 = is_inside(self.values[s0]);
                            if inside0 == is_inside(self.values[s1]) {
                                continue;
                            }
                            self.quads.push(Quad {
                                p,
                                axis: axis as u8,
                                inside0,
                            });
                        }
                    }
                }
            }
            black_box(&self.quads);
        }

        /// `HermiteCell::from_corners`: one crossing position and one gradient
        /// per cut edge. The crate's own function, called rather than copied.
        fn interpolate<S: Sdf<Scalar = R>>(&mut self, sdf: &S, g: Grid<R>) {
            self.hermite.clear();
            self.hermite.reserve(self.active.len());
            for cell in &self.active {
                let mut corner = [R::ZERO; 8];
                for (c, value) in corner.iter_mut().enumerate() {
                    let o = corner_offset(c as u8);
                    *value = self.values[g.dual_index([
                        cell.base[0] + o[0],
                        cell.base[1] + o[1],
                        cell.base[2] + o[2],
                    ])];
                }
                self.hermite.push(HermiteCell::from_corners(
                    sdf,
                    &corner,
                    g.cell_origin(cell.base),
                    g.cell_size,
                ));
            }
            black_box(&self.hermite);
        }

        /// `solve::solve_with` plus `Clamp::ToCell`, which is `Qef::place`.
        fn solve(&mut self, g: Grid<R>) {
            self.positions.clear();
            self.positions.reserve(self.hermite.len());
            let lambda = R::from_f64(solve::LAMBDA);
            for (cell, hermite) in self.active.iter().zip(&self.hermite) {
                let cell_origin = g.cell_origin(cell.base);
                let x = solve::solve_with(hermite, lambda)
                    .expect("an active cell has at least one crossing");
                self.positions
                    .push(clamp_to_cell(x, cell_origin, g.cell_size));
            }
            black_box(&self.positions);
        }

        /// `DualMesher::emit_vertices`' normal: the field's gradient at the
        /// solved position. Float work, and the second half of `interpolate`.
        fn vertex_normals<S: Sdf<Scalar = R>>(&mut self, sdf: &S) {
            self.normals.clear();
            self.normals.reserve(self.positions.len());
            for &position in &self.positions {
                let g = sdf.gradient(position);
                let inv = length(g).recip();
                self.normals.push([g[0] * inv, g[1] * inv, g[2] * inv]);
            }
            black_box(&self.normals);
        }

        fn emit_prepare(&mut self) {
            self.indices.clear();
            self.indices.reserve(self.quads.len() * 6);
            black_box(&self.indices);
        }

        /// `dual.rs:717-758`: four cells around the edge, then two triangles.
        ///
        /// `cell_edge_slot` is not mirrored, because under `Qef` every active
        /// cell's `push_whole_cell` sets all twelve slots to 0
        /// (`dual.rs:94-99`) and `slot_vertex[k] == k` for a reset sink — so the
        /// corner index *is* `cell_first`. That is a claim about the shipped
        /// code, and it is checked rather than argued: the mesh comparison is
        /// bit-exact and fails if it is wrong anywhere.
        fn emit_walk(&mut self, g: Grid<R>) {
            for quad in &self.quads {
                let axis = quad.axis as usize;
                let u = (axis + 1) % 3;
                let v = (axis + 2) % 3;
                let mut corners = [0u32; 4];
                for (slot, (du, dv)) in
                    corners
                        .iter_mut()
                        .zip([(0u32, 0u32), (1, 0), (1, 1), (0, 1)])
                {
                    let mut cell = quad.p;
                    cell[u] -= du;
                    cell[v] -= dv;
                    *slot = self.cell_first[g.cell_index(cell)];
                }
                if quad.inside0 {
                    self.indices
                        .extend_from_slice(&[corners[0], corners[1], corners[2]]);
                    self.indices
                        .extend_from_slice(&[corners[0], corners[2], corners[3]]);
                } else {
                    self.indices
                        .extend_from_slice(&[corners[0], corners[2], corners[1]]);
                    self.indices
                        .extend_from_slice(&[corners[0], corners[3], corners[2]]);
                }
            }
            black_box(&self.indices);
        }

        /// Stages, in pipeline order.
        const STAGES: usize = 7;

        /// [`Mc::prefix`]'s counterpart. Seven stages, because "crossing
        /// positions and normals" needs two windows on this path.
        fn prefix<S: Sdf<Scalar = R>>(&mut self, sdf: &S, g: Grid<R>, upto: usize) {
            if upto >= 1 {
                self.sample(sdf, g);
            }
            if upto >= 2 {
                self.classify(g);
            }
            if upto >= 3 {
                self.interpolate(sdf, g);
            }
            if upto >= 4 {
                self.solve(g);
            }
            if upto >= 5 {
                self.vertex_normals(sdf);
            }
            if upto >= 6 {
                self.emit_prepare();
            }
            if upto >= 7 {
                self.emit_walk(g);
            }
        }

        fn pipeline<S: Sdf<Scalar = R>>(&mut self, sdf: &S, g: Grid<R>) {
            self.prefix(sdf, g, Self::STAGES);
        }
    }

    // ─── counting ──────────────────────────────────────────────────────────

    /// Cycles and instructions from one or more windows.
    #[derive(Clone, Copy, Default)]
    struct Counted {
        cycles: f64,
        instructions: f64,
    }

    impl Counted {
        fn plus(self, other: Self) -> Self {
            Self {
                cycles: self.cycles + other.cycles,
                instructions: self.instructions + other.instructions,
            }
        }

        fn scaled(self, by: f64) -> Self {
            Self {
                cycles: self.cycles * by,
                instructions: self.instructions * by,
            }
        }

        /// The mean of two windows over the same work.
        ///
        /// `cycles_total` is taken twice — once before the stage windows and once
        /// after — and averaged, so that a clock drifting monotonically over the
        /// time a repetition takes cancels to first order instead of landing
        /// entirely in the residual.
        fn mean(self, other: Self) -> Self {
            self.plus(other).scaled(0.5)
        }
    }

    /// One counter window, undivided.
    ///
    /// The `perf_event` system calls are all outside the counted region, so they
    /// cannot inflate a stage or the residual. What *is* inside it is the cache
    /// and branch-predictor state the kernel round trip left behind, and that is
    /// the whole content of the residual — see the module docs.
    fn raw_window(probe: &mut Probe, body: impl FnOnce()) -> (Counted, f64) {
        probe.reset_and_enable();
        let started = Instant::now();
        body();
        let nanos = started.elapsed().as_nanos() as f64;
        probe.disable();
        let counted = probe.read();
        assert!(
            counted.worst_ratio() >= MIN_TIME_RATIO,
            "a counter ran only {:.1}% of the time it was enabled, so its value is an \
             extrapolation rather than a measurement",
            counted.worst_ratio() * 100.0
        );
        (
            Counted {
                cycles: counted.cycles.count as f64,
                instructions: counted.instructions.count as f64,
            },
            nanos,
        )
    }

    /// [`raw_window`] over `inner` repetitions of `body`, divided by `inner`.
    fn window(probe: &mut Probe, inner: usize, mut body: impl FnMut()) -> (Counted, f64) {
        let scale = 1.0 / inner as f64;
        let (counted, nanos) = raw_window(probe, || {
            for _ in 0..inner {
                body();
            }
        });
        (counted.scaled(scale), nanos * scale)
    }

    /// [`window`], for the callers that do not want the clock.
    fn cycles_of(probe: &mut Probe, inner: usize, body: impl FnMut()) -> Counted {
        window(probe, inner, body).0
    }

    /// The most stages any mirror has, so one array type covers both.
    const MAX_STAGES: usize = 7;

    /// One repetition, as measured: the prefix cuts, and two views of the whole.
    ///
    /// `cut[k]` is the cost of the first `k` stages of a pipeline pass, so
    /// `cut[0]` is zero, `cut[stages]` is the whole pipeline, and every stage is
    /// a difference of neighbouring cuts.
    #[derive(Clone, Copy, Default)]
    struct RepRaw {
        total: Counted,
        total_ns: f64,
        cut: [Counted; MAX_STAGES + 1],
        shipped: Counted,
    }

    /// The median of a set of readings.
    ///
    /// Taken **per quantity** rather than per repetition, which is a deliberate
    /// departure from `experiment_p15`'s median-run. The quantities here are
    /// *prefix cuts*, monotone in `k` by construction, and a single repetition
    /// disturbed by another process on the machine moves one cut and therefore
    /// two stages. Medianing each cut independently keeps the accounting exact —
    /// the stages still telescope to `cut[stages]`, because they are differences
    /// of the medianed cuts — while making no stage hostage to one repetition.
    fn median(pick: &dyn Fn(&RepRaw) -> f64, reps: &[RepRaw]) -> f64 {
        let mut values: Vec<f64> = reps.iter().map(pick).collect();
        values.sort_by(|a, b| a.total_cmp(b));
        values[values.len() / 2]
    }

    /// [`median`] over both counters of one quantity.
    fn median_counted(pick: &dyn Fn(&RepRaw) -> Counted, reps: &[RepRaw]) -> Counted {
        Counted {
            cycles: median(&|r| pick(r).cycles, reps),
            instructions: median(&|r| pick(r).instructions, reps),
        }
    }

    /// Everything one repetition measured.
    #[derive(Clone, Copy, Default)]
    struct Rep {
        /// One whole mirrored extraction in pipeline order, uninstrumented
        /// inside: the registered `cycles_total`, and the number comparable to
        /// the shipped extractor. Mean of the windows taken before and after the
        /// instrumented passes.
        total: Counted,
        total_ns: f64,
        sample: Counted,
        classify: Counted,
        interpolate: Counted,
        solve: Counted,
        emit: Counted,
        emit_prepare: Counted,
        shipped: Counted,
    }

    impl Rep {
        fn stages(self) -> f64 {
            self.sample.cycles
                + self.classify.cycles
                + self.interpolate.cycles
                + self.solve.cycles
                + self.emit.cycles
        }

        fn residual(self) -> f64 {
            self.total.cycles - self.stages()
        }

        fn residual_share(self) -> f64 {
            self.residual().abs() / self.total.cycles
        }
    }

    // ─── one row ───────────────────────────────────────────────────────────

    /// Which extractor a row decomposes.
    #[derive(Clone, Copy, PartialEq, Eq)]
    enum Which {
        MarchingCubes,
        DualContouring,
    }

    impl Which {
        fn name(self) -> &'static str {
            match self {
                Self::MarchingCubes => "marching_cubes",
                Self::DualContouring => "dual_contouring",
            }
        }
    }

    /// One measured `(field, resolution, scalar, extractor)`.
    struct Row {
        field: &'static str,
        resolution: u32,
        scalar: &'static str,
        extractor: Which,
        rep: Rep,
        residual_worst: f64,
        inner: usize,
        cells: usize,
        active_cells: usize,
        vertices: usize,
        triangles: usize,
        mesh_identical: bool,
        /// C3 is a comparison between two rows, so these are filled in once
        /// every row exists.
        classify_hold: f64,
        emit_collapse: f64,
        emit_collapse_instructions: f64,
        c3_holds: bool,
    }

    impl Row {
        fn emit_walk(&self) -> f64 {
            self.rep.emit.cycles - self.rep.emit_prepare.cycles
        }

        /// `emit_walk` in **instructions**, which is deterministic and
        /// cache-independent (`M-279`, and `experiment_p15`'s "load-bearing
        /// half"). C3's cycle ratio is a difference of two nearly equal windows
        /// on the surface-free arm, so its relative noise is unbounded; the
        /// instruction ratio is not, and C3 requires **both**.
        fn emit_walk_instructions(&self) -> f64 {
            self.rep.emit.instructions - self.rep.emit_prepare.instructions
        }

        fn float_share(&self) -> f64 {
            (self.rep.sample.cycles + self.rep.interpolate.cycles + self.rep.solve.cycles)
                / self.rep.total.cycles
        }

        fn integer_share(&self) -> f64 {
            (self.rep.classify.cycles + self.rep.emit.cycles) / self.rep.total.cycles
        }

        /// The mirror's cost against the shipped extractor's, on the same grid
        /// in the same run. `M-279`'s agreement check in cycles; the bit-exact
        /// mesh comparison is the other half.
        fn agreement_ratio(&self) -> f64 {
            self.rep.total.cycles / self.rep.shipped.cycles
        }

        fn group(&self) -> (u32, &'static str, Which) {
            (self.resolution, self.scalar, self.extractor)
        }
    }

    /// Bit for bit, as bit patterns rather than as values — `M-279`'s agreement
    /// check, and the licence for reading a share off a mirror.
    fn same<R: Real>(a: &[[R; 3]], b: &[[R; 3]]) -> bool {
        a.len() == b.len()
            && a.iter()
                .zip(b)
                .all(|(p, q)| (0..3).all(|k| p[k].as_f64().to_bits() == q[k].as_f64().to_bits()))
    }

    /// Measure one row.
    fn measure<R, S>(
        field: &'static str,
        scalar: &'static str,
        which: Which,
        n: u32,
        sdf: &S,
        origin: [R; 3],
        cell_size: R,
    ) -> Row
    where
        R: Real,
        S: Sdf<Scalar = R>,
    {
        let g = Grid {
            n,
            origin,
            cell_size,
        };
        let shape = RuntimeShape3::new([n; 3]).expect("the fixture fits u32");

        let mut mc = Mc::<R>::default();
        let mut dc = Dc::<R>::default();
        let mut shipped_mc = MarchingCubes::<R>::new();
        let mut shipped_dc = DualContouring::<R>::new();
        let mut out = MeshBuffer::<R>::new();

        macro_rules! shipped {
            () => {{
                out.reset();
                match which {
                    Which::MarchingCubes => shipped_mc
                        .extract(sdf, &shape, origin, cell_size, &mut out)
                        .expect("extraction"),
                    Which::DualContouring => shipped_dc
                        .extract(sdf, &shape, origin, cell_size, &mut out)
                        .expect("extraction"),
                }
                black_box(&out);
            }};
        }
        macro_rules! mirror {
            () => {
                match which {
                    Which::MarchingCubes => mc.pipeline(sdf, g),
                    Which::DualContouring => dc.pipeline(sdf, g),
                }
            };
        }

        for _ in 0..WARMUP {
            shipped!();
            mirror!();
        }

        // ── the batch, chosen from a timed pass ──────────────────────────────
        let started = Instant::now();
        mirror!();
        let pass_ns = started.elapsed().as_nanos() as f64;
        let inner = ((TARGET_BATCH_NS / pass_ns.max(1.0)).ceil() as usize).clamp(1, MAX_INNER);

        // ── the agreement the whole file rests on ───────────────────────────
        shipped!();
        let (vertices, triangles, active_cells, mesh_identical) = {
            let (positions, normals, indices) = match which {
                Which::MarchingCubes => (&mc.positions, &mc.normals, &mc.indices),
                Which::DualContouring => (&dc.positions, &dc.normals, &dc.indices),
            };
            let identical = same(positions, &out.positions)
                && same(normals, &out.normals)
                && indices.as_slice() == out.indices.as_slice();
            assert!(
                identical,
                "{field} {n}^3 {scalar} {}: the mirror's mesh differs from the shipped \
                 extractor's ({} vs {} vertices, {} vs {} indices) — every share below would \
                 be a share of a different algorithm",
                which.name(),
                positions.len(),
                out.positions.len(),
                indices.len(),
                out.indices.len()
            );
            let active = match which {
                Which::MarchingCubes => mc.active.len(),
                Which::DualContouring => dc.active.len(),
            };
            (positions.len(), indices.len() / 3, active, identical)
        };

        // ── REPS repetitions ────────────────────────────────────────────────
        let mut probe = Probe::open();
        let mut reps: Vec<RepRaw> = Vec::with_capacity(REPS);
        let stages = match which {
            Which::MarchingCubes => Mc::<R>::STAGES,
            Which::DualContouring => Dc::<R>::STAGES,
        };
        for _ in 0..REPS {
            let shipped_counts = cycles_of(&mut probe, inner, || shipped!());

            // `cycles_total` is a window over whole pipeline passes with nothing
            // instrumented inside. Taken twice, either side of the prefix sweep,
            // and averaged, so a monotonically drifting clock cancels.
            let (before, before_ns) = window(&mut probe, inner, || mirror!());

            let mut cut = [Counted::default(); MAX_STAGES + 1];
            for (upto, slot) in cut.iter_mut().enumerate().take(stages + 1).skip(1) {
                *slot = cycles_of(&mut probe, inner, || match which {
                    Which::MarchingCubes => mc.prefix(sdf, g, upto),
                    Which::DualContouring => dc.prefix(sdf, g, upto),
                });
            }

            let after = window(&mut probe, inner, || mirror!()).0;

            reps.push(RepRaw {
                total: before.mean(after),
                total_ns: before_ns,
                cut,
                shipped: shipped_counts,
            });
        }

        // Worst repetition, reported rather than scored on: it is the honest
        // measure of how much the machine moved under the instrument.
        let residual_worst = reps
            .iter()
            .map(|r| (r.total.cycles - r.cut[stages].cycles).abs() / r.total.cycles)
            .fold(0.0f64, f64::max);

        let mut cut = [Counted::default(); MAX_STAGES + 1];
        for (k, slot) in cut.iter_mut().enumerate().take(stages + 1).skip(1) {
            *slot = median_counted(&|r| r.cut[k], &reps);
        }
        let step = |k: usize| Counted {
            cycles: cut[k].cycles - cut[k - 1].cycles,
            instructions: cut[k].instructions - cut[k - 1].instructions,
        };
        let (sample, classify, interpolate, solve, emit, emit_prepare) = match which {
            // sample, classify, emit_prepare, emit_walk, interpolate.
            Which::MarchingCubes => (
                step(1),
                step(2),
                step(5),
                Counted::default(),
                step(3).plus(step(4)),
                step(3),
            ),
            // sample, classify, crossings, solve, normals, emit_prepare,
            // emit_walk. "Crossing positions and normals" is one registered
            // stage; on this path the vertex normal can only be taken after the
            // solve, so it is two prefix steps added together.
            Which::DualContouring => (
                step(1),
                step(2),
                step(3).plus(step(5)),
                step(4),
                step(6).plus(step(7)),
                step(6),
            ),
        };
        let rep = Rep {
            total: median_counted(&|r| r.total, &reps),
            total_ns: median(&|r| r.total_ns, &reps),
            sample,
            classify,
            interpolate,
            solve,
            emit,
            emit_prepare,
            shipped: median_counted(&|r| r.shipped, &reps),
        };

        // The registration's own vacuity control, asserted rather than merely
        // recorded: a stage that reads zero where it must run means the
        // instrument cannot see what it claims to.
        assert!(
            rep.sample.cycles > 0.0 && rep.classify.cycles > 0.0,
            "{field} {n}^3 {scalar} {}: sample or classify read zero cycles",
            which.name()
        );

        Row {
            field,
            resolution: n,
            scalar,
            extractor: which,
            rep,
            residual_worst,
            inner,
            cells: g.cell_count(),
            active_cells,
            vertices,
            triangles,
            mesh_identical,
            classify_hold: f64::NAN,
            emit_collapse: f64::NAN,
            emit_collapse_instructions: f64::NAN,
            c3_holds: false,
        }
    }

    /// Every `(field, resolution, scalar, extractor)` this row measures: eight
    /// reference fields plus the surface-free control, twice over in scalar and
    /// twice over in extractor, at both registered resolutions.
    fn sweep() -> Vec<Row> {
        let mut rows = Vec::new();
        for n in RESOLUTIONS {
            for which in [Which::MarchingCubes, Which::DualContouring] {
                isomesh::for_each_reference_field!(f32, |name, field| {
                    let (_, origin, cell_size) = crate::common::grid(&field, n);
                    rows.push(measure(name, "f32", which, n, &field, origin, cell_size));
                });
                {
                    let sphere = Sphere::<f32>::canonical();
                    let cell = 4.0f32 / (n - 1) as f32;
                    rows.push(measure(
                        SURFACE_FREE,
                        "f32",
                        which,
                        n,
                        &sphere,
                        [10.0; 3],
                        cell,
                    ));
                }
                isomesh::for_each_reference_field!(f64, |name, field| {
                    let (_, origin, cell_size) = crate::common::grid(&field, n);
                    rows.push(measure(name, "f64", which, n, &field, origin, cell_size));
                });
                {
                    let sphere = Sphere::<f64>::canonical();
                    let cell = 4.0f64 / f64::from(n - 1);
                    rows.push(measure(
                        SURFACE_FREE,
                        "f64",
                        which,
                        n,
                        &sphere,
                        [10.0; 3],
                        cell,
                    ));
                }
            }
        }
        rows
    }

    pub(crate) fn run(run: &mut crate::common::experiment::Run) {
        let mut rows = sweep();

        // ── C3: the surface-free control against its own `sphere` arm ────────
        let c3: Vec<(f64, f64, f64, bool)> = rows
            .iter()
            .map(|row| {
                let group = row.group();
                let sphere = rows
                    .iter()
                    .find(|r| r.group() == group && r.field == "sphere")
                    .expect("every group has a sphere arm");
                let free = rows
                    .iter()
                    .find(|r| r.group() == group && r.field == SURFACE_FREE)
                    .expect("every group has a surface-free arm");
                let hold = free.rep.classify.cycles / sphere.rep.classify.cycles;
                let collapse = free.emit_walk() / sphere.emit_walk();
                let collapse_i = free.emit_walk_instructions() / sphere.emit_walk_instructions();
                let holds = free.active_cells == 0
                    && free.vertices == 0
                    && free.triangles == 0
                    && free.rep.classify.cycles > 0.0
                    && hold >= CLASSIFY_HOLD_FLOOR
                    && collapse < EMIT_COLLAPSE_CEILING
                    && collapse_i < EMIT_COLLAPSE_CEILING;
                (hold, collapse, collapse_i, holds)
            })
            .collect();
        for (row, &(hold, collapse, collapse_i, holds)) in rows.iter_mut().zip(&c3) {
            row.classify_hold = hold;
            row.emit_collapse = collapse;
            row.emit_collapse_instructions = collapse_i;
            row.c3_holds = holds;
        }

        // ── C2: the maximum integer share over the eight reference fields ────
        let at_65: Vec<&Row> = rows
            .iter()
            .filter(|r| r.resolution == 65 && r.field != SURFACE_FREE)
            .collect();
        assert_eq!(
            at_65.len(),
            REFERENCE_FIELDS * 2 * 2,
            "eight reference fields x two scalars x two extractors at 65^3"
        );
        let max_integer_share = at_65
            .iter()
            .map(|r| r.integer_share())
            .fold(0.0f64, f64::max);
        let c2_holds = max_integer_share >= INTEGER_SHARE_BAR;
        // C1 is scored on **the row's** `residual_share`, which is what the
        // registration says: "residual_share under 5% of measured cycles on
        // EVERY row". `residual_worst` — the worst single repetition of any row —
        // is reported beside it, because on a machine running other work a
        // single repetition of a 30 ms window can swing far further than the
        // clause's bar, and hiding that would be worse than reporting it.
        let residual_row_worst = rows
            .iter()
            .map(|r| r.rep.residual_share())
            .fold(0.0f64, f64::max);
        let residual_rep_worst = rows.iter().map(|r| r.residual_worst).fold(0.0f64, f64::max);
        let c1_holds = residual_row_worst < RESIDUAL_CEILING;
        let c3_all = rows.iter().all(|r| r.c3_holds);

        // ── the tables ──────────────────────────────────────────────────────
        println!(
            "{:<20} {:>3} {:>4} {:>16} {:>6} {:>9} {:>7} {:>7} {:>7} {:>7} {:>7} {:>8} {:>7} \
             {:>7} {:>6}",
            "field",
            "n",
            "R",
            "extractor",
            "inner",
            "cyc/cell",
            "smpl%",
            "clsf%",
            "intp%",
            "solv%",
            "emit%",
            "resid%",
            "float%",
            "int%",
            "agree"
        );
        for r in &rows {
            let total = r.rep.total.cycles;
            let pct = |x: f64| x / total * 100.0;
            println!(
                "{:<20} {:>3} {:>4} {:>16} {:>6} {:>9.2} {:>7.2} {:>7.2} {:>7.2} {:>7.2} \
                 {:>7.2} {:>+8.2} {:>7.2} {:>7.2} {:>6.3}",
                r.field,
                r.resolution,
                r.scalar,
                r.extractor.name(),
                r.inner,
                total / r.cells as f64,
                pct(r.rep.sample.cycles),
                pct(r.rep.classify.cycles),
                pct(r.rep.interpolate.cycles),
                pct(r.rep.solve.cycles),
                pct(r.rep.emit.cycles),
                pct(r.rep.residual()),
                r.float_share() * 100.0,
                r.integer_share() * 100.0,
                r.agreement_ratio()
            );
        }

        println!(
            "\nC1 residual: worst |residual| / total over the 72 rows {residual_row_worst:.4} \
             (bar {RESIDUAL_CEILING}) -> {}",
            if c1_holds { "HELD" } else { "FALSIFIED" }
        );
        println!(
            "    worst single repetition of any row: {residual_rep_worst:.4} — reported, not \
             scored on; a 30 ms window on a machine running other work moves this far more than \
             the clause's bar"
        );
        println!(
            "C2 integer share: max over the eight reference fields at 65^3 = \
             {max_integer_share:.4} (bar {INTEGER_SHARE_BAR}) -> {}",
            if c2_holds { "HELD" } else { "FALSIFIED" }
        );
        println!(
            "    -> {}",
            if c2_holds {
                "Group A proceeds: P-103, P-104 and P-106 have a reachable denominator"
            } else {
                "Group A closes: P-103, P-104 and P-106 are Amdahl-dead at the measured share"
            }
        );

        println!(
            "\n{:<3} {:>4} {:>16} {:>14} {:>14} {:>14} {:>9} {:>9} {:>6}",
            "n",
            "R",
            "extractor",
            "classify_hold",
            "emit_cyc",
            "emit_instr",
            "free_act",
            "free_tri",
            "C3"
        );
        for r in rows.iter().filter(|r| r.field == SURFACE_FREE) {
            println!(
                "{:<3} {:>4} {:>16} {:>14.4} {:>14.6} {:>14.6} {:>9} {:>9} {:>6}",
                r.resolution,
                r.scalar,
                r.extractor.name(),
                r.classify_hold,
                r.emit_collapse,
                r.emit_collapse_instructions,
                r.active_cells,
                r.triangles,
                r.c3_holds
            );
        }
        println!(
            "C3 separability: classify_hold >= {CLASSIFY_HOLD_FLOOR}, emit_walk collapse < \
             {EMIT_COLLAPSE_CEILING} in BOTH cycles and instructions, nothing emitted -> {}",
            if c3_all { "HELD" } else { "FALSIFIED" }
        );

        for r in &rows {
            let total = r.rep.total.cycles;
            run.record(&[
                ("field", r.field.to_string()),
                ("resolution", r.resolution.to_string()),
                ("scalar", r.scalar.to_string()),
                ("extractor", r.extractor.name().to_string()),
                ("cycles_total", format!("{total:.1}")),
                ("cycles_sample", format!("{:.1}", r.rep.sample.cycles)),
                ("cycles_classify", format!("{:.1}", r.rep.classify.cycles)),
                (
                    "cycles_interpolate",
                    format!("{:.1}", r.rep.interpolate.cycles),
                ),
                ("cycles_solve", format!("{:.1}", r.rep.solve.cycles)),
                ("cycles_emit", format!("{:.1}", r.rep.emit.cycles)),
                ("cycles_residual", format!("{:.1}", r.rep.residual())),
                ("residual_share", format!("{:.6}", r.rep.residual_share())),
                ("float_share", format!("{:.6}", r.float_share())),
                ("integer_share", format!("{:.6}", r.integer_share())),
                (
                    "instructions_total",
                    format!("{:.1}", r.rep.total.instructions),
                ),
                (
                    "instructions_classify",
                    format!("{:.1}", r.rep.classify.instructions),
                ),
                (
                    "instructions_emit",
                    format!("{:.1}", r.rep.emit.instructions),
                ),
                ("cells", r.cells.to_string()),
                ("active_cells", r.active_cells.to_string()),
                ("c1_holds", c1_holds.to_string()),
                ("c2_holds", c2_holds.to_string()),
                ("c3_holds", r.c3_holds.to_string()),
                // ── extra columns (M-273) ──
                (
                    "residual_signed_share",
                    format!("{:.6}", r.rep.residual() / total),
                ),
                ("residual_share_worst", format!("{:.6}", r.residual_worst)),
                (
                    "cycles_emit_prepare",
                    format!("{:.1}", r.rep.emit_prepare.cycles),
                ),
                ("cycles_emit_walk", format!("{:.1}", r.emit_walk())),
                (
                    "instructions_emit_prepare",
                    format!("{:.1}", r.rep.emit_prepare.instructions),
                ),
                (
                    "instructions_emit_walk",
                    format!("{:.1}", r.emit_walk_instructions()),
                ),
                (
                    "instructions_sample",
                    format!("{:.1}", r.rep.sample.instructions),
                ),
                (
                    "instructions_interpolate",
                    format!("{:.1}", r.rep.interpolate.instructions),
                ),
                (
                    "instructions_solve",
                    format!("{:.1}", r.rep.solve.instructions),
                ),
                (
                    "cycles_per_cell_mirror",
                    format!("{:.4}", total / r.cells as f64),
                ),
                (
                    "cycles_per_cell_shipped",
                    format!("{:.4}", r.rep.shipped.cycles / r.cells as f64),
                ),
                ("agreement_ratio", format!("{:.4}", r.agreement_ratio())),
                ("mesh_identical_to_shipped", r.mesh_identical.to_string()),
                (
                    "instructions_per_cell_total",
                    format!("{:.4}", r.rep.total.instructions / r.cells as f64),
                ),
                (
                    "ipc_total",
                    format!("{:.4}", r.rep.total.instructions / total),
                ),
                ("vertices", r.vertices.to_string()),
                ("triangles", r.triangles.to_string()),
                (
                    "active_fraction",
                    format!("{:.6}", r.active_cells as f64 / r.cells as f64),
                ),
                ("inner_reps", r.inner.to_string()),
                (
                    "ns_per_cell",
                    format!("{:.4}", r.rep.total_ns / r.cells as f64),
                ),
                ("ghz", format!("{:.4}", total / r.rep.total_ns)),
                ("c3_classify_hold", format!("{:.4}", r.classify_hold)),
                ("c3_emit_walk_collapse", format!("{:.6}", r.emit_collapse)),
                (
                    "c3_emit_walk_collapse_instructions",
                    format!("{:.6}", r.emit_collapse_instructions),
                ),
                (
                    "residual_share_row_worst",
                    format!("{residual_row_worst:.6}"),
                ),
                (
                    "residual_share_rep_worst",
                    format!("{residual_rep_worst:.6}"),
                ),
                (
                    "max_reference_integer_share_at_65",
                    format!("{max_integer_share:.6}"),
                ),
            ]);
        }
    }
}

fn main() {
    if !std::env::args().any(|arg| arg == "--bench") {
        return;
    }

    let prereg = isomesh::experiment!("P-121");

    #[cfg(target_os = "linux")]
    common::experiment::run(prereg, experiment::run);

    // `experiment_p12`'s precedent, and the registration names it: this row's
    // instrument *is* `perf_event_open`, so off Linux there is nothing to
    // degrade to. A recorded zero would be a fabricated share, and a fabricated
    // share is what `✗51` cost.
    #[cfg(not(target_os = "linux"))]
    {
        eprintln!(
            "{} decomposes extraction with hardware performance counters, and this platform has \
             no `perf_event_open`. There is no clock substitute: M-281 forbids a nanosecond \
             carrying this verdict.",
            prereg.id
        );
        std::process::exit(1);
    }
}
