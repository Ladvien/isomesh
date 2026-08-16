//! **P-12 — is the dual's superlinear cost the crossed-edge gather?**
//!
//! Ticket: R-005. Pre-registered at R-000.
//!
//! ```bash
//! cargo bench --bench experiment_p12
//! ```
//!
//! Writes `docs/experiments/p-12.csv`. **Linux only** — see below.
//!
//! # The thing that has never been explained
//!
//! **M-21:** Surface Nets is not `O(n³)` over `16³…256³` and Marching Cubes is;
//! per-sample cost rises where the other's falls and flattens. **M-45:** it
//! reproduces on a second machine and gets *worse* there — Zen 3's Surface Nets
//! runs `37.38 → 49.08 ns/sample` against Marching Cubes' `15.18 → 13.19` — so
//! it is not one cache hierarchy. Two machines, same shape, mechanism unknown.
//! That is O-11, half-answered since T-006.
//!
//! P-12 names a mechanism: the dual gathers the **four cells around a crossed
//! edge**, which are at offsets `1`, `n` and `n²` in the sample grid, so the
//! working set of one gather grows with the grid while Marching Cubes' eight
//! corners stay within two rows. If that is it, cache misses per sample rise
//! with `n` for the dual and stay flat for Marching Cubes.
//!
//! # The falsifier names two alternatives, and both are measured here
//!
//! *"Flat miss rates, pointing at branch misprediction or allocation instead."*
//! So this reads `BRANCH_MISSES` in the same session, and `PAGE_FAULTS` as the
//! allocation proxy — fresh memory has to be faulted in on first touch, so an
//! allocation-driven cost shows there. A hand-written counting allocator is not
//! available: `unsafe_code = "forbid"` is a workspace lint and `GlobalAlloc`
//! cannot be implemented without `unsafe`, so page faults are the measurement
//! that exists rather than the one that would be ideal.
//!
//! # Why this bench is Linux-only, and what it does elsewhere
//!
//! Hardware counters come from `perf_event_open`, which is a Linux system call.
//! It needs no privileges here — `perf_event_paranoid = 2` still permits
//! user-space counting of one's own process — but it has no macOS equivalent
//! that a bench can call, and `perf` itself is not installed on this machine
//! either.
//!
//! On any other platform this **refuses and exits non-zero** rather than
//! reporting a fabricated column. That is the shape `CLAUDE.md` settles on for
//! the mesh-shader probe: a capability check that fails loudly is one path
//! choosing itself by measurement, where a silent substitute would be two.
//!
//! # The fixture is the resolution sweep's, deliberately
//!
//! Same field, same scalar, same [`common::grid`], same warmup and median rule
//! as `benches/resolution_sweep.rs`, so these rows compose with
//! `docs/measurements/resolution_sweep-ryzen9-5900x.csv` rather than describing
//! a different experiment. Four resolutions the ticket names, plus the points
//! around 128 that M-21's per-sample spike sits on.

mod common;

#[cfg(target_os = "linux")]
mod experiment {
    use std::hint::black_box;
    use std::time::Instant;

    use perf_event::events::{Cache, CacheOp, CacheResult, Hardware, Software};
    use perf_event::{Builder, Counter};

    use isomesh::dual_contouring::DualContouring;
    use isomesh::extractor::Extractor;
    use isomesh::fields::{ReferenceField, Sphere};
    use isomesh::marching_cubes::MarchingCubes;
    use isomesh::surface_nets::SurfaceNets;
    use isomesh::{MeshBuffer, RuntimeShape3};

    use crate::common;

    /// `f32`, because that is what a game passes and what the resolution sweep
    /// measured. Comparing against `f64` numbers would compare two experiments.
    type Scalar = f32;

    /// Samples per axis.
    ///
    /// The ticket asks for 96/128/192/256. The rest are here because M-21 and
    /// M-45 both show a per-sample **spike at 128³** on two unrelated machines
    /// and nobody has followed it; four points cannot tell a spike from a step.
    /// **127 and 129 flank it deliberately**: their working sets differ from
    /// 128's by 2%, so if only 128 misbehaves the cause is the power-of-two
    /// stride and not the size of anything.
    const RESOLUTIONS: [u32; 12] = [48, 64, 96, 112, 127, 128, 129, 144, 160, 192, 224, 256];

    /// The control that separates P-12's two candidates, and the whole reason
    /// this sweep has a second field.
    ///
    /// A sphere of radius 10 over a domain of half-extent 2 is **entirely
    /// inside**: every sample is negative, no grid edge is crossed, no vertex is
    /// placed and no quad is walked. The dense per-cell arrays are still
    /// allocated, filled and scanned in full.
    ///
    /// So the gather P-12 blames happens `O(n²)` times on `sphere` and **zero**
    /// times on `empty`, while the `O(n³)` dense state is identical. If the
    /// miss rate survives the control, the gather is not what is being paid for.
    const EMPTY_RADIUS: Scalar = 10.0;

    /// **The orientation control**, which turns a reading of the source into a
    /// measurement.
    ///
    /// `DualMesher::emit_quads` sweeps each of the three edge axes with the
    /// **innermost** loop over `v = (axis + 2) % 3`, so the axis-0 pass walks
    /// `values` at stride `nx·ny`, the axis-2 pass at stride `nx`, and only the
    /// axis-1 pass runs along `x`. That is a reading of the code; these three
    /// shapes are the check.
    ///
    /// Same sample count, same field, no surface — only the *order* memory is
    /// visited in differs. If the traversal is what costs, the shape whose
    /// strided passes reach furthest is the most expensive by a wide margin; if
    /// the cost is isotropic, all three are the same and the reading is wrong.
    ///
    /// Dimensions are deliberately **not** powers of two, so the 128³ aliasing
    /// this same sweep finds cannot contaminate the comparison.
    ///
    /// **Two sizes, and the small one is a control on the control.** The first
    /// three hold 4.3 M samples — a 17 MB `values` array, inside this machine's
    /// 32 MB L3 — where no traversal order can miss and the comparison
    /// therefore cannot discriminate. The second three hold 16.7 M, the same as
    /// the 256³ row, where it is out of cache and the order is free to matter.
    /// Running only the small one would have reported "orientation does not
    /// matter" from a fixture that could not have said anything else.
    const ORIENTATIONS: [[u32; 3]; 6] = [
        [68, 252, 252],
        [252, 68, 252],
        [252, 252, 68],
        [68, 496, 496],
        [496, 68, 496],
        [496, 496, 68],
    ];

    /// Untimed runs first, so the numbers are steady-state re-meshing rather
    /// than first-touch page faults on freshly grown scratch. **The page-fault
    /// column depends on this**: without it every row would report the cost of
    /// growing the buffers once.
    const WARMUP_RUNS: u32 = 2;

    /// Counted runs per configuration. The median by wall time is reported, and
    /// its own counters with it — mixing one run's time with another's misses
    /// would describe a run that never happened.
    const TIMED_RUNS: usize = 5;

    /// Below this, a counter was multiplexed and its value is an extrapolation.
    ///
    /// Zen 3 has six general-purpose counters and this opens five plus one
    /// software event, so nothing should be scheduled out. Asserted rather than
    /// hoped: a multiplexed count is a scaled estimate, and reporting one as a
    /// measurement is the kind of thing this file exists to prevent.
    const MIN_TIME_RATIO: f64 = 0.99;

    /// One counter, its label, and what it read.
    struct Reading {
        count: u64,
        /// `time_running / time_enabled`. Below 1 means the kernel had to share
        /// the counter and scaled the result.
        ratio: f64,
    }

    /// The six hardware events and one software event, opened together.
    ///
    /// Zen 3 has six general-purpose counters, so this is exactly full and
    /// nothing should be multiplexed; [`MIN_TIME_RATIO`] is what says so rather
    /// than hoping. `STALLED_CYCLES_BACKEND` would be the seventh and the one
    /// that would settle where the cycles go — it is **not available on this
    /// machine**, `perf_event_open` answering ENOENT, because AMD does not map
    /// the generic event.
    struct Probe {
        cycles: Counter,
        instructions: Counter,
        cache_misses: Counter,
        l1d_read_misses: Counter,
        /// Transparent huge pages are `always` here, so a 67 MB array is ~34
        /// pages and this should stay near zero. Measured rather than assumed,
        /// because "the working set grew" and "the page walk grew" are different
        /// mechanisms with the same symptom.
        dtlb_read_misses: Counter,
        branch_misses: Counter,
        page_faults: Counter,
    }

    /// What one counted run produced.
    struct Counts {
        cycles: Reading,
        instructions: Reading,
        cache_misses: Reading,
        l1d_read_misses: Reading,
        dtlb_read_misses: Reading,
        branch_misses: Reading,
        page_faults: Reading,
    }

    impl Probe {
        /// Open every counter, or say which one the kernel refused.
        ///
        /// # Panics
        ///
        /// If any counter cannot be opened. An experiment that silently drops an
        /// event reports a column it did not measure.
        fn open() -> Self {
            let hardware = |kind: Hardware| {
                Builder::new()
                    .kind(kind)
                    .build()
                    .unwrap_or_else(|e| panic!("perf_event_open for {kind:?}: {e}"))
            };
            Self {
                cycles: hardware(Hardware::CPU_CYCLES),
                instructions: hardware(Hardware::INSTRUCTIONS),
                cache_misses: hardware(Hardware::CACHE_MISSES),
                l1d_read_misses: Builder::new()
                    .kind(Cache {
                        which: perf_event::events::WhichCache::L1D,
                        operation: CacheOp::READ,
                        result: CacheResult::MISS,
                    })
                    .build()
                    .unwrap_or_else(|e| panic!("perf_event_open for L1D read miss: {e}")),
                dtlb_read_misses: Builder::new()
                    .kind(Cache {
                        which: perf_event::events::WhichCache::DTLB,
                        operation: CacheOp::READ,
                        result: CacheResult::MISS,
                    })
                    .build()
                    .unwrap_or_else(|e| panic!("perf_event_open for dTLB read miss: {e}")),
                branch_misses: hardware(Hardware::BRANCH_MISSES),
                page_faults: Builder::new()
                    .kind(Software::PAGE_FAULTS)
                    .build()
                    .unwrap_or_else(|e| panic!("perf_event_open for page faults: {e}")),
            }
        }

        fn each(&mut self) -> [&mut Counter; 7] {
            [
                &mut self.cycles,
                &mut self.instructions,
                &mut self.cache_misses,
                &mut self.l1d_read_misses,
                &mut self.dtlb_read_misses,
                &mut self.branch_misses,
                &mut self.page_faults,
            ]
        }

        fn reset_and_enable(&mut self) {
            for counter in self.each() {
                counter.reset().expect("reset");
                counter.enable().expect("enable");
            }
        }

        fn disable(&mut self) {
            for counter in self.each() {
                counter.disable().expect("disable");
            }
        }

        fn read(&mut self) -> Counts {
            fn one(counter: &mut Counter) -> Reading {
                let read = counter.read_count_and_time().expect("read");
                Reading {
                    count: read.count,
                    ratio: if read.time_enabled == 0 {
                        1.0
                    } else {
                        read.time_running as f64 / read.time_enabled as f64
                    },
                }
            }
            Counts {
                cycles: one(&mut self.cycles),
                instructions: one(&mut self.instructions),
                cache_misses: one(&mut self.cache_misses),
                l1d_read_misses: one(&mut self.l1d_read_misses),
                dtlb_read_misses: one(&mut self.dtlb_read_misses),
                branch_misses: one(&mut self.branch_misses),
                page_faults: one(&mut self.page_faults),
            }
        }
    }

    impl Counts {
        /// The worst scheduling ratio of the six.
        fn worst_ratio(&self) -> f64 {
            [
                self.cycles.ratio,
                self.instructions.ratio,
                self.cache_misses.ratio,
                self.l1d_read_misses.ratio,
                self.dtlb_read_misses.ratio,
                self.branch_misses.ratio,
                self.page_faults.ratio,
            ]
            .into_iter()
            .fold(1.0f64, f64::min)
        }
    }

    /// One timed, counted run.
    struct Run {
        nanos: u128,
        counts: Counts,
        triangles: usize,
    }

    /// Measure one extractor at one resolution on one field.
    fn measure<E: Extractor<Scalar>>(
        extractor: &mut E,
        field: &Sphere<Scalar>,
        samples: u32,
    ) -> Run {
        let (shape, origin, cell_size) = common::grid(field, samples);
        measure_on(extractor, field, &shape, origin, cell_size)
    }

    /// Measure one extractor on an explicit grid.
    fn measure_on<E: Extractor<Scalar>>(
        extractor: &mut E,
        field: &Sphere<Scalar>,
        shape: &RuntimeShape3,
        origin: [Scalar; 3],
        cell_size: Scalar,
    ) -> Run {
        let mut mesh = MeshBuffer::<Scalar>::new();

        for _ in 0..WARMUP_RUNS {
            extractor
                .extract_into(field, shape, origin, cell_size, &mut mesh)
                .expect("extraction");
            black_box(&mesh);
        }

        let mut runs: Vec<Run> = Vec::with_capacity(TIMED_RUNS);
        let mut probe = Probe::open();
        for _ in 0..TIMED_RUNS {
            probe.reset_and_enable();
            let started = Instant::now();
            extractor
                .extract_into(field, shape, origin, cell_size, &mut mesh)
                .expect("extraction");
            let nanos = started.elapsed().as_nanos();
            probe.disable();
            black_box(&mesh);
            runs.push(Run {
                nanos,
                counts: probe.read(),
                triangles: mesh.triangle_count(),
            });
        }
        runs.sort_by_key(|r| r.nanos);
        runs.swap_remove(TIMED_RUNS / 2)
    }

    /// Run the sweep and record it.
    pub(crate) fn run(run: &mut crate::common::experiment::Run) {
        let fields = [
            ("sphere", Sphere::<Scalar>::canonical()),
            (
                "empty",
                Sphere::<Scalar> {
                    center: [0.0; 3],
                    radius: EMPTY_RADIUS,
                },
            ),
        ];
        println!(
            "{:<8} {:<16} {:>5} {:>10} {:>10} {:>10} {:>10} {:>10} {:>9}",
            "field",
            "extractor",
            "n",
            "ns/sample",
            "LLC miss/s",
            "L1D miss/s",
            "br miss/s",
            "pf/sample",
            "triangles"
        );
        for (field_name, field) in &fields {
            for samples in RESOLUTIONS {
                for (name, measured) in [
                    ("marching_cubes", {
                        let mut e = MarchingCubes::<Scalar>::new();
                        measure(&mut e, field, samples)
                    }),
                    ("surface_nets", {
                        let mut e = SurfaceNets::<Scalar>::new();
                        measure(&mut e, field, samples)
                    }),
                    ("dual_contouring", {
                        let mut e = DualContouring::<Scalar>::new();
                        measure(&mut e, field, samples)
                    }),
                ] {
                    let n = f64::from(samples);
                    let total = n * n * n;
                    let ratio = measured.counts.worst_ratio();
                    assert!(
                        ratio >= MIN_TIME_RATIO,
                        "{name} at {samples}³: a counter ran only {:.1}% of the time it was \
                         enabled, so its value is a scaled estimate rather than a measurement",
                        ratio * 100.0
                    );

                    let per = |c: &Reading| c.count as f64 / total;
                    let ns_per_sample = measured.nanos as f64 / total;
                    println!(
                        "{field_name:<8} {name:<16} {samples:>5} {ns_per_sample:>10.3} \
                         {:>10.4} {:>10.4} {:>10.4} {:>10.6} {:>9}",
                        per(&measured.counts.cache_misses),
                        per(&measured.counts.l1d_read_misses),
                        per(&measured.counts.branch_misses),
                        per(&measured.counts.page_faults),
                        measured.triangles,
                    );

                    run.record(&[
                        ("samples", format!("{total:.0}")),
                        ("extractor", name.to_string()),
                        (
                            "cache_misses_per_sample",
                            format!("{:.6}", per(&measured.counts.cache_misses)),
                        ),
                        ("ns_per_sample", format!("{ns_per_sample:.4}")),
                        ("field", (*field_name).to_string()),
                        ("samples_per_axis", samples.to_string()),
                        ("shape", format!("{samples}x{samples}x{samples}")),
                        (
                            "l1d_read_misses_per_sample",
                            format!("{:.6}", per(&measured.counts.l1d_read_misses)),
                        ),
                        (
                            "dtlb_read_misses_per_sample",
                            format!("{:.6}", per(&measured.counts.dtlb_read_misses)),
                        ),
                        (
                            "branch_misses_per_sample",
                            format!("{:.6}", per(&measured.counts.branch_misses)),
                        ),
                        (
                            "page_faults_per_sample",
                            format!("{:.8}", per(&measured.counts.page_faults)),
                        ),
                        (
                            "instructions_per_sample",
                            format!("{:.3}", per(&measured.counts.instructions)),
                        ),
                        (
                            "cycles_per_sample",
                            format!("{:.3}", per(&measured.counts.cycles)),
                        ),
                        (
                            "ipc",
                            format!(
                                "{:.3}",
                                measured.counts.instructions.count as f64
                                    / measured.counts.cycles.count as f64
                            ),
                        ),
                        ("triangles", measured.triangles.to_string()),
                        ("counter_time_ratio", format!("{ratio:.4}")),
                    ]);
                }
            }
        }
        orientations(run);
        println!(
            "\n`f32`, the resolution sweep's own grid, median of {TIMED_RUNS} runs after \
             {WARMUP_RUNS} warmups.\n`empty` is the same code on a field with no surface at all: \
             the dense O(n³) state is identical and\nthe crossed-edge gather happens zero times. \
             `LLC miss/s` is PERF_COUNT_HW_CACHE_MISSES per sample;\n`pf/sample` is page faults, \
             the allocation proxy."
        );
    }

    /// Same sample count, three axis orders, no surface. See [`ORIENTATIONS`].
    fn orientations(run: &mut crate::common::experiment::Run) {
        let field = Sphere::<Scalar> {
            center: [0.0; 3],
            radius: EMPTY_RADIUS,
        };
        let (lo, hi) = ReferenceField::domain(&field);
        println!(
            "\n{:<8} {:<16} {:>14} {:>12} {:>10} {:>10} {:>10}",
            "field", "extractor", "shape", "samples", "ns/sample", "LLC miss/s", "L1D miss/s"
        );
        for size in ORIENTATIONS {
            let shape = RuntimeShape3::new(size).expect("the fixture fits u32");
            // One spacing per group, from its longest axis, so the three grids
            // of a group cover comparable world volume and only their strides
            // differ.
            let longest = size[0].max(size[1]).max(size[2]);
            let cell_size = (hi[0] - lo[0]) / longest as Scalar;
            let label = format!("{}x{}x{}", size[0], size[1], size[2]);
            for (name, measured) in [
                ("marching_cubes", {
                    let mut e = MarchingCubes::<Scalar>::new();
                    measure_on(&mut e, &field, &shape, lo, cell_size)
                }),
                ("surface_nets", {
                    let mut e = SurfaceNets::<Scalar>::new();
                    measure_on(&mut e, &field, &shape, lo, cell_size)
                }),
                ("dual_contouring", {
                    let mut e = DualContouring::<Scalar>::new();
                    measure_on(&mut e, &field, &shape, lo, cell_size)
                }),
            ] {
                let total = f64::from(size[0]) * f64::from(size[1]) * f64::from(size[2]);
                let ratio = measured.counts.worst_ratio();
                assert!(
                    ratio >= MIN_TIME_RATIO,
                    "{name} on {label}: a counter ran only {:.1}% of the time it was enabled",
                    ratio * 100.0
                );
                let per = |c: &Reading| c.count as f64 / total;
                let ns_per_sample = measured.nanos as f64 / total;
                println!(
                    "{:<8} {name:<16} {label:>14} {total:>12.0} {ns_per_sample:>10.3} {:>10.4} \
                     {:>10.4}",
                    "empty",
                    per(&measured.counts.cache_misses),
                    per(&measured.counts.l1d_read_misses),
                );
                run.record(&[
                    ("samples", format!("{total:.0}")),
                    ("extractor", name.to_string()),
                    (
                        "cache_misses_per_sample",
                        format!("{:.6}", per(&measured.counts.cache_misses)),
                    ),
                    ("ns_per_sample", format!("{ns_per_sample:.4}")),
                    ("field", "empty".to_string()),
                    ("samples_per_axis", "0".to_string()),
                    ("shape", label.clone()),
                    (
                        "l1d_read_misses_per_sample",
                        format!("{:.6}", per(&measured.counts.l1d_read_misses)),
                    ),
                    (
                        "dtlb_read_misses_per_sample",
                        format!("{:.6}", per(&measured.counts.dtlb_read_misses)),
                    ),
                    (
                        "branch_misses_per_sample",
                        format!("{:.6}", per(&measured.counts.branch_misses)),
                    ),
                    (
                        "page_faults_per_sample",
                        format!("{:.8}", per(&measured.counts.page_faults)),
                    ),
                    (
                        "instructions_per_sample",
                        format!("{:.3}", per(&measured.counts.instructions)),
                    ),
                    (
                        "cycles_per_sample",
                        format!("{:.3}", per(&measured.counts.cycles)),
                    ),
                    (
                        "ipc",
                        format!(
                            "{:.3}",
                            measured.counts.instructions.count as f64
                                / measured.counts.cycles.count as f64
                        ),
                    ),
                    ("triangles", measured.triangles.to_string()),
                    ("counter_time_ratio", format!("{ratio:.4}")),
                ]);
            }
        }
        println!(
            "\nSame sample count, no surface, three axis orders. `emit_quads` walks `values` at \
             stride\n`nx·ny` on its axis-0 pass and `nx` on its axis-2 pass, so only the order \
             differs here."
        );
    }
}

fn main() {
    if !std::env::args().any(|arg| arg == "--bench") {
        return;
    }

    let prereg = isomesh::experiment!("P-12");

    #[cfg(target_os = "linux")]
    common::experiment::run(prereg, experiment::run);

    #[cfg(not(target_os = "linux"))]
    {
        eprintln!(
            "{} needs hardware performance counters, and this platform has no `perf_event_open`.\n\
             Refusing rather than writing a `cache_misses_per_sample` column that was not measured.\n\
             Run it on Linux; `perf_event_paranoid = 2` is permissive enough and no root is needed.",
            prereg.id
        );
        std::process::exit(1);
    }
}
