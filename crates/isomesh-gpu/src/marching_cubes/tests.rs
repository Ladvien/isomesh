//! The GPU kernel against the CPU extractor, on the same grid and the same
//! table.
//!
//! Exact comparison where the two really must agree — triangle counts, and the
//! table itself. Geometry is compared as a *set of points* rather than
//! index-for-index, because the CPU shares vertices through a grid-edge cache
//! and this emits a soup; that divergence is documented on the module rather
//! than smoothed over here.
#![allow(clippy::float_cmp)]

use isomesh::fields::{BoxExact, Sphere, Torus};
use isomesh::marching_cubes::{MarchingCubes, table};
use isomesh::{MeshBuffer, RuntimeShape3, Sdf};

use super::{MarchingCubesGpu, case_table_bytes};
use crate::headless::Gpu;
use crate::{FieldBuffer, GridParams, read_buffer, read_buffer_u32};

fn gpu() -> &'static Gpu {
    crate::headless::shared()
}

/// The CPU's answer on the same grid, for comparison.
fn cpu_mesh<F: Sdf<Scalar = f32>>(field: &F, grid: GridParams) -> MeshBuffer<f32> {
    let mut out = MeshBuffer::new();
    let shape = RuntimeShape3::new(grid.samples()).expect("valid shape");
    MarchingCubes::<f32>::new()
        .extract(field, &shape, grid.origin(), grid.cell_size(), &mut out)
        .expect("cpu extraction");
    out
}

/// Rule 5's guarantee, made mechanical: the bytes handed to the shader are
/// exactly `isomesh`'s table, unpacked again and compared entry by entry.
///
/// If this passes there is no second case table to drift, which is the whole
/// reason the table is uploaded rather than written in WGSL.
#[test]
fn the_uploaded_table_is_isomeshs_own() {
    let bytes = case_table_bytes();
    assert_eq!(bytes.len(), 256 * (1 + table::MAX_TRIANGLES) * 4);

    let word = |i: usize| {
        u32::from_le_bytes([
            bytes[i * 4],
            bytes[i * 4 + 1],
            bytes[i * 4 + 2],
            bytes[i * 4 + 3],
        ])
    };
    let stride = 1 + table::MAX_TRIANGLES;
    for (case, entry) in table::CASES.iter().enumerate() {
        let header = word(case * stride);
        assert_eq!(header & 0xff, u32::from(entry.count), "case {case} count");
        assert_eq!(
            (header >> 8) & 0xff,
            u32::from(entry.centroids),
            "case {case} centroids"
        );
        for (t, tri) in entry.triangles.iter().enumerate() {
            let packed = word(case * stride + 1 + t);
            assert_eq!(packed & 0xff, u32::from(tri[0]), "case {case} tri {t}.0");
            assert_eq!(
                (packed >> 8) & 0xff,
                u32::from(tri[1]),
                "case {case} tri {t}.1"
            );
            assert_eq!(
                (packed >> 16) & 0xff,
                u32::from(tri[2]),
                "case {case} tri {t}.2"
            );
        }
    }
}

/// The headline: the GPU emits the same number of triangles as the CPU, on
/// every reference field it is run against.
///
/// Triangle count is the strongest thing that *must* be equal. Both sides read
/// the same samples, classify with the same table and iterate the same cells,
/// so a disagreement means the kernel misread the table or the grid — and it
/// would be invisible in a picture.
#[test]
fn the_triangle_count_matches_the_cpu_on_every_field() {
    let gpu = gpu();
    let mc = MarchingCubesGpu::new(gpu.device(), gpu.queue()).expect("pipeline");

    let grid = GridParams::new([33; 3], [-2.0; 3], 0.125).expect("valid grid");
    let fields: [(&str, &dyn Sdf<Scalar = f32>); 3] = [
        ("sphere", &Sphere::<f32>::canonical()),
        ("torus", &Torus::<f32>::canonical()),
        (
            "box_exact",
            &BoxExact::<f32> {
                center: [0.0; 3],
                half_extents: [1.1, 0.9, 0.7],
            },
        ),
    ];

    for (name, field) in fields {
        let buffer = FieldBuffer::sampled(gpu.device(), gpu.queue(), grid, &field).expect("upload");
        let mesh = mc
            .extract(gpu.device(), gpu.queue(), &buffer)
            .expect("extract");
        let cpu = cpu_mesh(&field, grid);

        assert_eq!(
            mesh.triangle_count(),
            cpu.indices.len() / 3,
            "{name}: gpu and cpu disagree on triangle count"
        );
        assert!(mesh.triangle_count() > 0, "{name}: nothing was extracted");
        assert_eq!(mesh.positions.len(), mesh.normals.len());
    }
}

/// Every GPU vertex is a CPU vertex to within **one ULP per axis**, and most
/// are bit-identical.
///
/// Both sides evaluate the same expressions — `t = a / (a - b)` and
/// `lo + (hi - lo) * t`, over `origin + h * index` — on identical `f32`
/// samples, so bit-equality was the expectation. It does not hold, and the size
/// of the miss is the finding: **1 ULP, on 6% of vertices** (M-142). WGSL
/// permits a multiply-add to be contracted into a fused one and this adapter's
/// compiler takes that permission, which rounds once where the CPU rounds
/// twice.
///
/// So the assertion is the measured bound rather than equality, and the bound
/// is tight: two ULPs would fail this.
///
/// The neighbour search is over the 27 combinations of ±1 ULP per axis rather
/// than a nearest-point scan, which makes it exact instead of a tolerance and
/// linear instead of quadratic.
#[test]
fn every_gpu_vertex_is_a_cpu_vertex_to_within_one_ulp() {
    let gpu = gpu();
    let mc = MarchingCubesGpu::new(gpu.device(), gpu.queue()).expect("pipeline");
    let grid = GridParams::new([33; 3], [-2.0; 3], 0.125).expect("valid grid");
    let field = Sphere::<f32>::canonical();

    let buffer = FieldBuffer::sampled(gpu.device(), gpu.queue(), grid, &field).expect("upload");
    let mesh = mc
        .extract(gpu.device(), gpu.queue(), &buffer)
        .expect("extract");
    let cpu = cpu_mesh(&field, grid);

    let key = |p: &[f32; 3]| [p[0].to_bits(), p[1].to_bits(), p[2].to_bits()];
    let cpu_points: std::collections::HashSet<[u32; 3]> = cpu.positions.iter().map(key).collect();

    // One ULP up or down is the neighbouring bit pattern, for any finite float
    // of either sign as long as it is not a zero crossing -- and these are
    // surface positions on a sphere of radius 1, so none are near zero in a way
    // that matters. Signed-magnitude means "one step toward zero" is a
    // decrement in magnitude, which `wrapping_add` on the raw bits gives for
    // positives and mirrors for negatives; both directions are tried, so the
    // sign handling cannot bias the result.
    let step = |bits: u32, delta: i32| -> u32 { bits.wrapping_add(delta as u32) };

    let mut exact = 0usize;
    let mut within = 0usize;
    let mut strangers = 0usize;
    for p in &mesh.positions {
        let k = key(p);
        if cpu_points.contains(&k) {
            exact += 1;
            continue;
        }
        let mut found = false;
        for dx in [-1i32, 0, 1] {
            for dy in [-1i32, 0, 1] {
                for dz in [-1i32, 0, 1] {
                    let probe = [step(k[0], dx), step(k[1], dy), step(k[2], dz)];
                    if cpu_points.contains(&probe) {
                        found = true;
                    }
                }
            }
        }
        if found {
            within += 1;
        } else {
            strangers += 1;
        }
    }

    println!(
        "gpu vertices: {exact} bit-exact, {within} within 1 ULP, {strangers} further ({} total)",
        mesh.positions.len()
    );
    assert_eq!(
        strangers, 0,
        "{strangers} gpu vertices are more than one ULP from any cpu vertex"
    );
    // Non-zero on both sides, or the bound above is being asserted over an
    // empty set on one of them.
    assert!(
        exact > 0,
        "no vertex agreed exactly -- suspect a formula, not rounding"
    );
    assert!(
        within > 0,
        "every vertex agreed exactly -- this adapter does not contract, and the \
         1-ULP finding should be re-read before being trusted elsewhere"
    );
}

/// The same input twice gives the same bytes.
///
/// The two-pass design exists for this: an atomic bump allocator would order
/// the output by whichever workgroup arrived first.
#[test]
fn extraction_is_deterministic() {
    let gpu = gpu();
    let mc = MarchingCubesGpu::new(gpu.device(), gpu.queue()).expect("pipeline");
    let grid = GridParams::new([25; 3], [-2.0; 3], 0.16).expect("valid grid");
    let field = Torus::<f32>::canonical();
    let buffer = FieldBuffer::sampled(gpu.device(), gpu.queue(), grid, &field).expect("upload");

    let first = mc
        .extract(gpu.device(), gpu.queue(), &buffer)
        .expect("extract");
    let second = mc
        .extract(gpu.device(), gpu.queue(), &buffer)
        .expect("extract");
    // Geometry, explicitly. `GpuMesh` is not `PartialEq` precisely because it
    // carries wall-clock timings, which never repeat and are not part of what
    // "the same mesh" means.
    assert_eq!(first.positions, second.positions, "positions differ");
    assert_eq!(first.normals, second.normals, "normals differ");
}

/// Normals point away from the solid, which is the convention every consumer
/// reads them under.
///
/// Checked against the field's analytic gradient rather than against the CPU's
/// normals: the divergence between central differences and the analytic
/// gradient is the documented one (M-65), so the useful assertion is that they
/// agree in *direction*, well inside that measured 0.46° bound.
#[test]
fn normals_face_away_from_the_solid() {
    let gpu = gpu();
    let mc = MarchingCubesGpu::new(gpu.device(), gpu.queue()).expect("pipeline");
    let grid = GridParams::new([33; 3], [-2.0; 3], 0.125).expect("valid grid");
    let field = Sphere::<f32>::canonical();
    let buffer = FieldBuffer::sampled(gpu.device(), gpu.queue(), grid, &field).expect("upload");
    let mesh = mc
        .extract(gpu.device(), gpu.queue(), &buffer)
        .expect("extract");

    let mut worst = 1.0f32;
    for (p, n) in mesh.positions.iter().zip(&mesh.normals) {
        let g = field.gradient(*p);
        let len = (g[0] * g[0] + g[1] * g[1] + g[2] * g[2]).sqrt();
        assert!(len > 0.0, "the test field has no gradient at {p:?}");
        let dot = (n[0] * g[0] + n[1] * g[1] + n[2] * g[2]) / len;
        worst = worst.min(dot);
    }
    assert!(
        worst > 0.99,
        "worst normal agreement with the analytic gradient was {worst}, i.e. {:.2} deg",
        worst.acos().to_degrees()
    );
}

/// An empty field produces an empty mesh rather than an error or a panic.
#[test]
fn a_field_with_no_surface_extracts_nothing() {
    let gpu = gpu();
    let mc = MarchingCubesGpu::new(gpu.device(), gpu.queue()).expect("pipeline");
    let grid = GridParams::new([9; 3], [10.0; 3], 0.1).expect("valid grid");
    let buffer = FieldBuffer::sampled(gpu.device(), gpu.queue(), grid, &Sphere::<f32>::canonical())
        .expect("upload");

    let mesh = mc
        .extract(gpu.device(), gpu.queue(), &buffer)
        .expect("extract");
    assert_eq!(mesh.triangle_count(), 0);
    assert!(mesh.positions.is_empty());
}

/// Bit-identity is a property of the **cell size**, not of the port, and the
/// bound that survives either way is a distance.
///
/// At `h = 0.125` almost every vertex is bit-identical; at `h = 0.1` almost
/// none are. Both are correct — every vertex is still within a rounding
/// distance of a CPU vertex — but a test written only at a power of two would
/// report near-perfect agreement and hide a real defect, which is exactly what
/// happened: `GridParams::sample_position` used `mul_add` where `isomesh` uses
/// `origin + h * i`, and `h = 0.125` made the two forms bit-identical (M-143).
///
/// So the durable assertion is the *offset*, in cells, and it is checked at a
/// spacing where the two forms can disagree.
#[test]
fn the_agreement_is_a_distance_and_holds_at_a_non_power_of_two_spacing() {
    let gpu = gpu();
    let mc = MarchingCubesGpu::new(gpu.device(), gpu.queue()).expect("pipeline");
    let field = Sphere::<f32>::canonical();

    for (samples, cell) in [(33u32, 0.125f32), (41, 0.1)] {
        let grid = GridParams::new([samples; 3], [-2.0; 3], cell).expect("valid grid");
        let buffer = FieldBuffer::sampled(gpu.device(), gpu.queue(), grid, &field).expect("upload");
        let mesh = mc
            .extract(gpu.device(), gpu.queue(), &buffer)
            .expect("extract");
        let cpu = cpu_mesh(&field, grid);

        assert_eq!(
            mesh.triangle_count(),
            cpu.indices.len() / 3,
            "h = {cell}: triangle counts differ"
        );

        // Nearest CPU vertex through a spatial hash: a ULP probe only answers
        // "within k" for the k you thought to ask for.
        let bucket = cell * 0.01;
        let cell_of = |p: &[f32; 3]| {
            [
                (p[0] / bucket).floor() as i64,
                (p[1] / bucket).floor() as i64,
                (p[2] / bucket).floor() as i64,
            ]
        };
        let mut index: std::collections::HashMap<[i64; 3], Vec<[f32; 3]>> =
            std::collections::HashMap::new();
        for p in &cpu.positions {
            index.entry(cell_of(p)).or_default().push(*p);
        }

        let mut exact = 0usize;
        let mut worst = 0.0f32;
        for p in &mesh.positions {
            let mut best = f32::INFINITY;
            let home = cell_of(p);
            for dx in -1..=1 {
                for dy in -1..=1 {
                    for dz in -1..=1 {
                        let at = [home[0] + dx, home[1] + dy, home[2] + dz];
                        for q in index.get(&at).into_iter().flatten() {
                            let d = ((p[0] - q[0]).powi(2)
                                + (p[1] - q[1]).powi(2)
                                + (p[2] - q[2]).powi(2))
                            .sqrt();
                            best = best.min(d);
                        }
                    }
                }
            }
            if best == 0.0 {
                exact += 1;
            }
            worst = worst.max(best / cell);
        }

        // Rounding is many orders of magnitude below a cell; a disagreement
        // about geometry is not.
        assert!(
            worst < 1e-4,
            "h = {cell}: worst offset {worst:e} cells is too large to be rounding"
        );
        println!(
            "h = {cell}: {exact} of {} bit-identical, worst offset {worst:e} cells",
            mesh.positions.len()
        );
    }
}

/// `extract_buffers` is the same extraction with the last wait removed.
///
/// Same triangle count as the read-back path, and `geometry_readback_ms` is
/// zero by construction rather than by luck — if it were not, the "drawing
/// needs no read-back" claim in `mesh_render` would be false.
#[test]
fn extract_buffers_matches_extract_and_skips_the_geometry_readback() {
    let gpu = gpu();
    let mc = MarchingCubesGpu::new(gpu.device(), gpu.queue()).expect("pipeline");
    let grid = GridParams::new([33; 3], [-2.0; 3], 0.125).expect("valid grid");
    let field = Sphere::<f32>::canonical();
    let buffer = FieldBuffer::sampled(gpu.device(), gpu.queue(), grid, &field).expect("upload");

    let read_back = mc
        .extract(gpu.device(), gpu.queue(), &buffer)
        .expect("extract");
    let on_gpu = mc
        .extract_buffers(gpu.device(), gpu.queue(), &buffer)
        .expect("extract_buffers");

    assert_eq!(
        on_gpu.triangles as usize,
        read_back.triangle_count(),
        "the two entry points disagree about how many triangles there are"
    );
    assert_eq!(
        on_gpu.timings.geometry_readback_ms, 0.0,
        "extract_buffers read the geometry back after all"
    );
    // And the read-back path really does pay that cost, or the assertion above
    // is comparing against a stage nobody runs.
    assert!(
        read_back.timings.geometry_readback_ms > 0.0,
        "extract did not time its geometry read-back"
    );

    // The buffers are bindable as storage, which is what a mesh shader needs.
    assert!(
        on_gpu
            .positions
            .usage()
            .contains(wgpu::BufferUsages::STORAGE),
        "positions cannot be bound to a shader"
    );
    assert!(
        on_gpu.normals.usage().contains(wgpu::BufferUsages::STORAGE),
        "normals cannot be bound to a shader"
    );
}

/// An empty surface still hands back bindable buffers.
///
/// wgpu rejects a zero-sized binding, so a caller that always binds -- which a
/// renderer does -- would fail on an empty chunk if these were zero-length.
#[test]
fn empty_extraction_still_yields_bindable_buffers() {
    let gpu = gpu();
    let mc = MarchingCubesGpu::new(gpu.device(), gpu.queue()).expect("pipeline");
    let grid = GridParams::new([9; 3], [10.0; 3], 0.1).expect("valid grid");
    let buffer = FieldBuffer::sampled(gpu.device(), gpu.queue(), grid, &Sphere::<f32>::canonical())
        .expect("upload");

    let on_gpu = mc
        .extract_buffers(gpu.device(), gpu.queue(), &buffer)
        .expect("extract_buffers");
    assert_eq!(on_gpu.triangles, 0);
    assert!(
        on_gpu.positions.size() > 0,
        "a zero-sized binding is invalid"
    );
    assert!(on_gpu.normals.size() > 0);
}

/// The zero-read-back path agrees with the synchronising one, triangle for
/// triangle and byte for byte.
///
/// `extract_indirect` never learns the total, so the only way to check it is to
/// read what it produced afterwards — which is exactly the explicit,
/// caller-chosen check its contract describes.
#[test]
fn the_indirect_path_produces_the_same_geometry() {
    let gpu = gpu();
    let mc = MarchingCubesGpu::new(gpu.device(), gpu.queue()).expect("pipeline");
    let grid = GridParams::new([33; 3], [-2.0; 3], 0.125).expect("grid");
    let field = Sphere::<f32>::canonical();
    let buffer = FieldBuffer::sampled(gpu.device(), gpu.queue(), grid, &field).expect("upload");

    let sync = mc
        .extract_buffers(gpu.device(), gpu.queue(), &buffer)
        .expect("extract_buffers");
    // Generous budget: this test is about agreement, not truncation.
    let indirect = mc
        .extract_indirect(gpu.device(), gpu.queue(), &buffer, sync.triangles * 2)
        .expect("extract_indirect");

    // The total the GPU wrote, read back only because a test has to look.
    let total = read_buffer_u32(gpu.device(), gpu.queue(), &indirect.total, 4).expect("total")[0];
    assert_eq!(
        total, sync.triangles,
        "the indirect path counted differently"
    );

    // And the draw arguments it derived from it.
    let args = read_buffer_u32(gpu.device(), gpu.queue(), &indirect.indirect, 12).expect("args");
    assert_eq!(args[0], total.div_ceil(32), "workgroup count");
    assert_eq!([args[1], args[2]], [1, 1], "y and z must be 1");

    let params =
        read_buffer_u32(gpu.device(), gpu.queue(), &indirect.draw_params, 4).expect("params");
    assert_eq!(params[0], total, "the draw uniform must carry the total");

    // The geometry itself, bit for bit over the triangles that exist.
    let floats = u64::from(total) * 9 * 4;
    let a = read_buffer(gpu.device(), gpu.queue(), &sync.positions, floats).expect("sync");
    let b = read_buffer(gpu.device(), gpu.queue(), &indirect.positions, floats).expect("indirect");
    assert_eq!(
        a.iter().map(|v| v.to_bits()).collect::<Vec<_>>(),
        b.iter().map(|v| v.to_bits()).collect::<Vec<_>>(),
        "the two paths wrote different positions"
    );
}

/// A budget smaller than the surface truncates, and says so where the caller
/// can find it.
///
/// The important half is the second: truncation is not silent, because
/// `total` holds the real count and comparing it against the budget is the
/// documented check. The first half is that the shader stops at the buffer's
/// own length rather than writing past it — which without the `arrayLength`
/// bound would be memory corruption rather than a short mesh.
#[test]
fn a_budget_smaller_than_the_surface_truncates_detectably() {
    let gpu = gpu();
    let mc = MarchingCubesGpu::new(gpu.device(), gpu.queue()).expect("pipeline");
    let grid = GridParams::new([33; 3], [-2.0; 3], 0.125).expect("grid");
    let field = Sphere::<f32>::canonical();
    let buffer = FieldBuffer::sampled(gpu.device(), gpu.queue(), grid, &field).expect("upload");

    let full = mc
        .extract_buffers(gpu.device(), gpu.queue(), &buffer)
        .expect("extract");
    let budget = full.triangles / 4;
    let tight = mc
        .extract_indirect(gpu.device(), gpu.queue(), &buffer, budget)
        .expect("extract_indirect");

    let total = read_buffer_u32(gpu.device(), gpu.queue(), &tight.total, 4).expect("total")[0];
    assert_eq!(
        total, full.triangles,
        "the total must report the real surface, not the budget"
    );
    assert!(
        total > tight.budget,
        "the fixture did not actually overflow: {total} against a budget of {}",
        tight.budget
    );

    // The draw the GPU wrote for itself must stop where the emit pass stopped
    // writing: at the budget, not at the total. Before this clamp the indirect
    // args dispatched ceil(total/32) workgroups over budget-sized buffers --
    // 4x past the end of everything the emit pass wrote.
    let args = read_buffer_u32(gpu.device(), gpu.queue(), &tight.indirect, 12).expect("args");
    assert_eq!(
        args[0],
        tight.budget.div_ceil(32),
        "workgroups must cover the budget only"
    );
    let dp = read_buffer_u32(gpu.device(), gpu.queue(), &tight.draw_params, 4).expect("params");
    assert_eq!(
        dp[0], tight.budget,
        "the draw uniform must stop at the budget"
    );
}

/// The two WGSL constants that ARE hand-transcribed -- `EDGE_CORNERS` and the
/// `corner_offset` bit rule -- pinned to the same source of truth as the
/// uploaded table. The 12 pairs are parsed back out of the shader text and
/// compared to `cube::EDGE_CORNERS` (public through `marching_cubes::table`),
/// in order, because the position is the edge index.
///
/// Device-free on purpose, like `the_uploaded_table_is_isomeshs_own`: CI has
/// no GPU, and a drifted transcription would otherwise ship invisibly.
#[test]
fn the_wgsl_edge_corners_are_isomeshs_own() {
    let src = crate::MARCHING_CUBES_WGSL;
    let start = src
        .find("const EDGE_CORNERS")
        .expect("marching_cubes.wgsl lost `const EDGE_CORNERS`");
    let body = &src[start..];
    let end = body.find(");").expect("`const EDGE_CORNERS` never closes");
    let mut rest = &body[..end];

    let mut pairs: Vec<[u8; 2]> = Vec::new();
    while let Some(open) = rest.find("vec2<u32>(") {
        let args = &rest[open + "vec2<u32>(".len()..];
        let close = args.find(')').expect("unclosed vec2 in EDGE_CORNERS");
        let mut nums = args[..close].split(',').map(|n| {
            n.trim()
                .strip_suffix('u')
                .and_then(|digits| digits.parse::<u8>().ok())
                .expect("EDGE_CORNERS entries are `<digits>u`")
        });
        pairs.push([
            nums.next().expect("vec2 lost its low corner"),
            nums.next().expect("vec2 lost its high corner"),
        ]);
        rest = &args[close..];
    }

    assert_eq!(
        pairs,
        table::EDGE_CORNERS,
        "the WGSL EDGE_CORNERS transcription drifted from cube::EDGE_CORNERS"
    );
}

/// The corner bit layout, pinned as text: corner `c` at
/// `(c & 1, (c >> 1) & 1, (c >> 2) & 1)` is the crate's public corner
/// numbering contract (documented on `cube::CORNER_COUNT`), and this is the
/// exact expression the shader's `corner_offset` must compute. Text rather
/// than values because `cube::corner_offset` itself is `pub(crate)`.
#[test]
fn the_wgsl_corner_offset_is_the_documented_bit_layout() {
    assert!(
        crate::MARCHING_CUBES_WGSL.contains("vec3<u32>(c & 1u, (c >> 1u) & 1u, (c >> 2u) & 1u)"),
        "marching_cubes.wgsl's corner_offset no longer matches the corner bit layout"
    );
}

/// R-069's instrument, shown measuring rather than assumed to.
///
/// The whole point of `with_timestamps` is that `execute` stops being a
/// wall-clock guess, so this asserts the four things that would make it a
/// column named and not measured: the feature is reachable on this adapter, the
/// period is positive, both compute passes report a span, and each span is
/// strictly positive. A device that advertises `TIMESTAMP_QUERY` and resolves
/// zeros would pass a `>= 0` check and mean nothing.
#[test]
fn the_timestamp_instrument_reports_a_positive_span_for_each_compute_pass() {
    let Ok(gpu) = Gpu::with_timestamps() else {
        // Not a silent skip: the adapter's own answer is the finding, and this
        // host's is on record in P-71's registration (RTX 3090, all three
        // timestamp features true). A machine without it should not fail a
        // suite it cannot run, but it must say so.
        std::println!("TIMESTAMP_QUERY unavailable on this adapter; span check skipped");
        return;
    };
    let mc =
        MarchingCubesGpu::with_timestamps(gpu.device(), gpu.queue()).expect("timestamped pipeline");

    // No passes recorded yet, and that is a different answer from "no query set".
    let before = mc
        .take_timestamps(gpu.device(), gpu.queue())
        .expect("resolve")
        .expect("this extractor has a query set");
    assert!(before.spans.is_empty(), "nothing has dispatched yet");
    assert!(before.period_ns > 0.0, "a zero period cannot time anything");

    let grid = GridParams::new([17; 3], [-2.0; 3], 0.25).expect("grid");
    let field = FieldBuffer::sampled(gpu.device(), gpu.queue(), grid, &Sphere::<f32>::canonical())
        .expect("field");
    let mesh = mc
        .extract(gpu.device(), gpu.queue(), &field)
        .expect("extraction");
    assert!(
        !mesh.positions.is_empty(),
        "the fixture must produce geometry"
    );

    let spans = mc
        .take_timestamps(gpu.device(), gpu.queue())
        .expect("resolve")
        .expect("this extractor has a query set");
    assert!(spans.complete, "the query set overflowed");
    assert_eq!(
        spans.spans.len(),
        2,
        "one span per compute pass: count and emit — got {:?}",
        spans.spans
    );
    for span in &spans.spans {
        assert!(
            span.ms > 0.0,
            "{} measured {} ms: a zero span is a driver that advertises the \
             feature and does not implement it",
            span.label,
            span.ms
        );
    }
    // Resolving clears, so a second read is empty rather than a repeat.
    let after = mc
        .take_timestamps(gpu.device(), gpu.queue())
        .expect("resolve")
        .expect("query set")
        .spans;
    assert!(after.is_empty(), "resolve must clear: got {after:?}");
    std::println!(
        "timestamp period {:.4} ns; spans {:?}",
        spans.period_ns,
        spans.spans
    );
}
