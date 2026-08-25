//! The wasm module behind the front page's interactive demo: eight reference
//! fields, five extractors, and the validity report, in about 100 KB.
//!
//! # Why this crate exists at all
//!
//! The nine playable demos on the site are Bevy builds, and each one is ~36 MB
//! because it carries a renderer, an asset system and a window backend. **None of
//! that is `isomesh`.** The core crate is `no_std` with `libm` as its only normal
//! dependency, and this module is what that costs on its own: field in, triangles
//! out, with a hand-written WebGL2 renderer in `web/lite.js` on the other side of
//! the ABI.
//!
//! # The ABI is `extern "C"` over linear memory, and there is no `wasm-bindgen`
//!
//! [`MeshBuffer`]'s three fields are `pub` precisely so a consumer can read them
//! without a copy. Typed-array views over `instance.exports.memory` are both
//! smaller than bindgen's glue and fewer moving parts, so the exports below hand
//! out pointers and counts and the JS builds the views.
//!
//! **Every pointer returned here is valid only until the next [`iso_mesh`]
//! call**, which may reallocate the buffers or grow the module's memory. Growth
//! detaches `memory.buffer`, so a `Float32Array` captured before a re-mesh
//! silently reads a dead buffer; `lite.js` re-creates every view after every
//! call for that reason.
//!
//! # Every export is panic-free, and that is a hard rule rather than a habit
//!
//! `[profile.release]` sets `panic = "abort"`, so a panic is not an exception the
//! page can catch -- it is a dead module and a blank canvas with nothing on the
//! console to explain it. So there is no `unwrap`, no unchecked indexing, no
//! arithmetic that can overflow a cast, every input is clamped or rejected, and a
//! poisoned lock returns the same sentinel a refusal does.
//!
//! A refusal is one path: [`iso_mesh`] returns `0`, and the mesh and the report
//! are both cleared. It never leaves the previous mesh's numbers on screen beside
//! a triangle count of zero.

use std::sync::Mutex;

use isomesh::dual_contouring::DualContouring;
use isomesh::fields::{
    BoxExact, CappedGyroid, CsgDifference, FbmTerrain, NoiseCavity, ReferenceField, Sphere,
    ThinPlate, Torus, capped_gyroid, csg_difference, noise_cavity,
};
use isomesh::manifold_dual_contouring::ManifoldDualContouring;
use isomesh::marching_cubes::MarchingCubes;
use isomesh::marching_tetrahedra::MarchingTetrahedra;
use isomesh::surface_nets::SurfaceNets;
use isomesh::validate::{ValidateConfig, validate};
use isomesh::{MeshBuffer, RuntimeShape3, Sdf};

/// A boxed field, which is what lets one extractor instantiation serve all eight.
///
/// `crates/isomesh/src/sdf.rs`'s blanket `impl<S: Sdf + ?Sized> Sdf for Box<S>`
/// is what makes this work: `extract` is generic over the field, so eight
/// concrete types would monomorphise every extractor eight times -- five copies
/// of Marching Cubes becomes forty, in a module whose whole point is its size.
/// One virtual call per sample is the price, and this demo is not the benchmark;
/// `crates/isomesh/benches` is, and it uses concrete types for exactly that
/// reason.
type BoxedField = Box<dyn Sdf<Scalar = f32>>;

/// The dropdown labels for the fields, in the order [`resolve_field`] resolves.
///
/// Taken from the trait constants rather than written out, so the site cannot
/// disagree with the crate about what a field is called.
/// `the_field_names_are_the_registry_in_order` holds this against
/// `for_each_reference_field!`.
const FIELD_NAMES: [&str; 8] = [
    <Sphere<f32> as ReferenceField>::NAME,
    <Torus<f32> as ReferenceField>::NAME,
    <BoxExact<f32> as ReferenceField>::NAME,
    <CsgDifference<f32> as ReferenceField>::NAME,
    <ThinPlate<f32> as ReferenceField>::NAME,
    <CappedGyroid<f32> as ReferenceField>::NAME,
    <FbmTerrain<f32> as ReferenceField>::NAME,
    <NoiseCavity<f32> as ReferenceField>::NAME,
];

/// The dropdown labels for the extractors, in the order [`extract_into`] matches.
///
/// String literals rather than constants because `isomesh` publishes the registry
/// as one array of names rather than a name per type, and
/// `the_extractor_names_partition_the_registry` proves this list plus
/// [`NOT_OFFERED`] is exactly that array -- so a seventh extractor in the core
/// crate fails this crate's tests rather than quietly missing from the page.
const EXTRACTOR_NAMES: [&str; 5] = [
    "marching_cubes",
    "marching_tetrahedra",
    "surface_nets",
    "dual_contouring",
    "manifold_dual_contouring",
];

/// The registry entries this demo deliberately does not offer, and why.
///
/// - **`marching_cubes+decider`** is not a sixth algorithm. It is
///   [`MarchingCubes`] with `FaceAmbiguity::AsymptoticDecider` set, so offering it
///   as a peer of the others would present a configuration as an implementation.
/// - **`subgrid_marching_tetrahedra`** is ~196× Marching Cubes, and this runs on
///   the browser's main thread. Seconds of block on a front page is a broken
///   page, not a slow one.
///
/// `GreedyQuads` is absent from both lists because it is absent from the registry
/// itself: it is `isomesh::extractor::UNREGISTERED`, an axis-aligned blocky
/// surface rather than an isosurface, and
/// `greedy_quads_is_excluded_by_the_crate_not_by_this_demo` checks that this
/// remains the crate's decision rather than becoming ours.
///
/// `cfg(test)` because the exclusions are a decision recorded in prose and
/// checked by a test; nothing in the shipped module reads them, and a name list
/// compiled into a wasm binary for no reader is bytes on someone's front page.
#[cfg(test)]
const NOT_OFFERED: [&str; 2] = ["marching_cubes+decider", "subgrid_marching_tetrahedra"];

/// The narrowest grid worth meshing: eight cells per axis.
const MIN_SAMPLES: u32 = 9;

/// The widest grid this will mesh on the main thread.
///
/// 49³ is `subgrid`'s registered sampling resolution and the shootout's middle
/// rung, so a number read off this page is comparable with the committed CSVs.
/// Marching Tetrahedra at 49³ is ~3× the triangles of Marching Cubes, which is
/// the real ceiling here: past this the *renderer* is what stalls, not the
/// mesher.
const MAX_SAMPLES: u32 = 49;

/// Field names. The `kind` argument to [`iso_name`] and [`iso_name_len`].
const KIND_FIELD: u32 = 0;

/// Extractor names. The `kind` argument to [`iso_name`] and [`iso_name_len`].
const KIND_EXTRACTOR: u32 = 1;

/// Everything the page can read, behind one lock.
///
/// A `Mutex` rather than a `static mut` because `Mutex::new` is const, so the
/// static needs no initialiser call and no `unsafe` block -- which is what keeps
/// `#[unsafe(no_mangle)]` the only `unsafe` token in this crate. The browser's
/// main thread is the only caller, so the lock is never contended; it is here to
/// make the aliasing sound rather than to coordinate anything.
static DEMO: Mutex<Demo> = Mutex::new(Demo::new());

/// The last mesh, the camera framing it needs, and its validity report.
#[derive(Debug)]
struct Demo {
    /// Positions, normals and indices. Reused across calls, never reallocated
    /// once it has grown -- which is the whole reason `MeshBuffer` is
    /// caller-provided.
    mesh: MeshBuffer<f32>,
    /// The domain's centre, so the camera has something to orbit.
    centre: [f32; 3],
    /// Half the domain's longest axis, so the camera has a distance.
    extent: f32,
    /// `referenced_vertices − edges + faces`.
    euler: i32,
    /// Edges incident to three or more faces. Zero is the claim; a non-zero
    /// count on `gyroid` under a dual method is the measured, documented
    /// one-vertex-per-cell defect rather than a bug in this page.
    non_manifold_edges: u32,
    /// Edges incident to exactly one face. Non-zero for `fbm_terrain`, which
    /// exits through the sides of its own domain by construction.
    boundary_edges: u32,
    /// Triangles with numerically zero area. Recorded, never gated: Marching
    /// Cubes emits slivers whenever a corner value sits near zero.
    degenerate_triangles: u32,
}

impl Demo {
    /// An empty demo that has allocated nothing.
    const fn new() -> Self {
        Self {
            mesh: MeshBuffer::new(),
            centre: [0.0; 3],
            extent: 1.0,
            euler: 0,
            non_manifold_edges: 0,
            boundary_edges: 0,
            degenerate_triangles: 0,
        }
    }

    /// Drop the mesh and the report together.
    ///
    /// Called on every refusal, so a zero triangle count is never displayed
    /// beside the previous mesh's χ.
    fn clear(&mut self) {
        self.mesh.reset();
        self.centre = [0.0; 3];
        self.extent = 1.0;
        self.euler = 0;
        self.non_manifold_edges = 0;
        self.boundary_edges = 0;
        self.degenerate_triangles = 0;
    }
}

/// Read one value out of the demo, or `sentinel` if the lock is poisoned.
///
/// A poisoned lock means an earlier call panicked, and under `panic = "abort"`
/// that cannot have happened -- the module would be gone. It is handled anyway
/// because the alternative is `unwrap`, and an `unwrap` in an export is a panic
/// waiting for an input nobody thought of.
fn read<T>(sentinel: T, f: impl FnOnce(&Demo) -> T) -> T {
    match DEMO.lock() {
        Ok(demo) => f(&demo),
        Err(_) => sentinel,
    }
}

/// `usize` to `u32`, saturating.
///
/// Saturation rather than a wrap: a truncated count would make the page draw
/// part of a mesh and report the wrong size for it, which is worse than drawing
/// a clamped one. Unreachable in practice -- `MAX_SAMPLES` bounds the mesh far
/// below `u32::MAX` -- and here so that the bound is enforced rather than
/// assumed.
fn count(n: usize) -> u32 {
    u32::try_from(n).unwrap_or(u32::MAX)
}

/// `u64` to `u32`, saturating, for the report's counters.
fn counter(n: u64) -> u32 {
    u32::try_from(n).unwrap_or(u32::MAX)
}

/// One of the eight reference fields, boxed, with the domain it is meant to be
/// sampled over.
///
/// The domain is read from the concrete type *before* boxing, because
/// [`ReferenceField`] is not part of the `dyn Sdf` the extractors take -- and it
/// is the field that knows where its surface is, not this module.
fn resolve_field(index: u32) -> Option<(BoxedField, [f32; 3], [f32; 3])> {
    fn boxed<F>(field: F) -> (BoxedField, [f32; 3], [f32; 3])
    where
        F: ReferenceField + Sdf<Scalar = f32> + 'static,
    {
        let (lo, hi) = field.domain();
        (Box::new(field), lo, hi)
    }

    Some(match index {
        0 => boxed(Sphere::<f32>::canonical()),
        1 => boxed(Torus::<f32>::canonical()),
        2 => boxed(BoxExact::<f32>::canonical()),
        3 => boxed(csg_difference::<f32>()),
        4 => boxed(ThinPlate::<f32>::canonical()),
        5 => boxed(capped_gyroid::<f32>()),
        6 => boxed(FbmTerrain::<f32>::canonical()),
        7 => boxed(noise_cavity::<f32>()),
        _ => return None,
    })
}

/// Run one extractor over one boxed field, or `None` if the index names no
/// extractor or the extraction is refused.
///
/// A refusal is returned rather than absorbed. `isomesh` never substitutes a
/// degraded mesh, and neither does this.
fn extract_into(
    extractor: u32,
    field: &BoxedField,
    shape: &RuntimeShape3,
    origin: [f32; 3],
    cell_size: f32,
    mesh: &mut MeshBuffer<f32>,
) -> Option<()> {
    let extracted = match extractor {
        0 => MarchingCubes::<f32>::new().extract(field, shape, origin, cell_size, mesh),
        1 => MarchingTetrahedra::<f32>::new().extract(field, shape, origin, cell_size, mesh),
        2 => SurfaceNets::<f32>::new().extract(field, shape, origin, cell_size, mesh),
        3 => DualContouring::<f32>::new().extract(field, shape, origin, cell_size, mesh),
        4 => ManifoldDualContouring::<f32>::new().extract(field, shape, origin, cell_size, mesh),
        _ => return None,
    };
    extracted.ok()
}

/// Mesh `field` with `extractor` at `samples` per axis, and return the triangle
/// count.
///
/// `0` means refused, and a refusal clears the mesh and the report both. There
/// are four ways to be refused and no fifth: an index that names no field, an
/// index that names no extractor, a grid whose index space `RuntimeShape3`
/// rejects, and an extraction the algorithm itself declines.
///
/// `samples` is clamped into `9..=49` rather than rejected, because it arrives
/// from a range slider and clamping is what a slider means. The grid is derived
/// exactly as every bench in the repository derives it -- `cell_size =
/// (hi[0] − lo[0]) / (samples − 1)`, origin `lo` -- so the counts this returns
/// are comparable with the committed CSVs rather than merely plausible.
#[unsafe(no_mangle)]
pub extern "C" fn iso_mesh(field: u32, extractor: u32, samples: u32) -> u32 {
    let Ok(mut demo) = DEMO.lock() else {
        return 0;
    };
    demo.clear();

    let Some((sdf, lo, hi)) = resolve_field(field) else {
        return 0;
    };
    let samples = samples.clamp(MIN_SAMPLES, MAX_SAMPLES);
    let Ok(shape) = RuntimeShape3::new([samples; 3]) else {
        return 0;
    };
    // `samples >= MIN_SAMPLES` holds by the clamp above, so the divisor is at
    // least 8 and this is neither a division by zero nor a subtraction that can
    // underflow.
    let cell_size = (hi[0] - lo[0]) / (samples - 1) as f32;
    let Ok(config) = ValidateConfig::from_cell_size(f64::from(cell_size)) else {
        return 0;
    };

    let demo = &mut *demo;
    if extract_into(extractor, &sdf, &shape, lo, cell_size, &mut demo.mesh).is_none() {
        demo.clear();
        return 0;
    }

    let report = validate(&demo.mesh, &config);
    demo.centre = [
        (lo[0] + hi[0]) * 0.5,
        (lo[1] + hi[1]) * 0.5,
        (lo[2] + hi[2]) * 0.5,
    ];
    demo.extent = ((hi[0] - lo[0]).max(hi[1] - lo[1]).max(hi[2] - lo[2])) * 0.5;
    demo.euler = i32::try_from(report.euler_characteristic).unwrap_or(i32::MAX);
    demo.non_manifold_edges = counter(report.non_manifold_edges);
    demo.boundary_edges = counter(report.boundary_edges);
    demo.degenerate_triangles = counter(report.degenerate_triangles);
    count(demo.mesh.indices.len() / 3)
}

/// The vertex positions: `3 × iso_vertex_count()` floats, `xyz` interleaved.
///
/// Null if the lock is poisoned. Dangling-but-aligned when the mesh is empty,
/// which is what `Vec::as_ptr` returns for an unallocated buffer -- so the count
/// is what says whether there is anything to read, never the pointer.
#[unsafe(no_mangle)]
pub extern "C" fn iso_positions() -> *const f32 {
    read(std::ptr::null(), |demo| {
        demo.mesh.positions.as_ptr().cast::<f32>()
    })
}

/// The vertex normals: `3 × iso_vertex_count()` floats, unit length, parallel to
/// [`iso_positions`].
#[unsafe(no_mangle)]
pub extern "C" fn iso_normals() -> *const f32 {
    read(std::ptr::null(), |demo| {
        demo.mesh.normals.as_ptr().cast::<f32>()
    })
}

/// The triangle indices: `iso_index_count()` `u32`s, flat triples.
#[unsafe(no_mangle)]
pub extern "C" fn iso_indices() -> *const u32 {
    read(std::ptr::null(), |demo| demo.mesh.indices.as_ptr())
}

/// Vertices in the last mesh.
#[unsafe(no_mangle)]
pub extern "C" fn iso_vertex_count() -> u32 {
    read(0, |demo| count(demo.mesh.positions.len()))
}

/// Indices in the last mesh, i.e. three per triangle.
#[unsafe(no_mangle)]
pub extern "C" fn iso_index_count() -> u32 {
    read(0, |demo| count(demo.mesh.indices.len()))
}

/// The domain's centre: three floats, for the camera to orbit.
#[unsafe(no_mangle)]
pub extern "C" fn iso_centre() -> *const f32 {
    read(std::ptr::null(), |demo| demo.centre.as_ptr())
}

/// Half the domain's longest axis, for the camera's distance.
///
/// `1.0` when nothing has been meshed, so a camera initialised before the first
/// [`iso_mesh`] is merely wrong about the scale rather than at the origin looking
/// at itself.
#[unsafe(no_mangle)]
pub extern "C" fn iso_extent() -> f32 {
    read(1.0, |demo| demo.extent)
}

/// `referenced_vertices − edges + faces` for the last mesh.
///
/// `2` for a sphere, `0` for a torus, and the number to watch when the extractor
/// changes under a fixed field: a dual method that shares one vertex between two
/// surface sheets moves this.
#[unsafe(no_mangle)]
pub extern "C" fn iso_euler() -> i32 {
    read(0, |demo| demo.euler)
}

/// Edges incident to three or more faces in the last mesh. Zero is the claim.
#[unsafe(no_mangle)]
pub extern "C" fn iso_non_manifold_edges() -> u32 {
    read(0, |demo| demo.non_manifold_edges)
}

/// Edges incident to exactly one face in the last mesh.
#[unsafe(no_mangle)]
pub extern "C" fn iso_boundary_edges() -> u32 {
    read(0, |demo| demo.boundary_edges)
}

/// Triangles of numerically zero area in the last mesh.
#[unsafe(no_mangle)]
pub extern "C" fn iso_degenerate_triangles() -> u32 {
    read(0, |demo| demo.degenerate_triangles)
}

/// How many fields [`iso_mesh`] accepts.
#[unsafe(no_mangle)]
pub extern "C" fn iso_field_count() -> u32 {
    count(FIELD_NAMES.len())
}

/// How many extractors [`iso_mesh`] accepts.
#[unsafe(no_mangle)]
pub extern "C" fn iso_extractor_count() -> u32 {
    count(EXTRACTOR_NAMES.len())
}

/// The name table `kind` selects, or `None` for any other `kind`.
fn names(kind: u32) -> Option<&'static [&'static str]> {
    match kind {
        KIND_FIELD => Some(&FIELD_NAMES),
        KIND_EXTRACTOR => Some(&EXTRACTOR_NAMES),
        _ => None,
    }
}

/// The UTF-8 bytes of one name, for `kind` `0` (fields) or `1` (extractors).
///
/// Null for an unknown `kind` or an out-of-range `index`. The bytes are
/// `'static` and outlive every call, unlike the mesh pointers.
///
/// This exists so the two `<select>` elements have **one** source of truth. A
/// duplicated label list in JavaScript is a thing that rots silently -- the
/// dropdown says `gyroid` and the module meshes `fbm_terrain` -- and a count
/// check would not catch it because the counts would still agree.
#[unsafe(no_mangle)]
pub extern "C" fn iso_name(kind: u32, index: u32) -> *const u8 {
    let Some(table) = names(kind) else {
        return std::ptr::null();
    };
    match table.get(index as usize) {
        Some(name) => name.as_ptr(),
        None => std::ptr::null(),
    }
}

/// The byte length of the name [`iso_name`] returns, or `0` if there is none.
#[unsafe(no_mangle)]
pub extern "C" fn iso_name_len(kind: u32, index: u32) -> u32 {
    let Some(table) = names(kind) else {
        return 0;
    };
    match table.get(index as usize) {
        Some(name) => count(name.len()),
        None => 0,
    }
}

#[cfg(test)]
mod tests;
