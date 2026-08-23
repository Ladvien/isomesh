//! T-007 — golden hashes over every (algorithm, field, resolution).
//!
//! A hash of the exact bits of every position, normal and index, committed as
//! `golden_hashes.json` beside the manifest. It answers one question the rest of
//! the suite cannot: **did the output change at all?**
//!
//! Everything else here checks a *property*. The validity harness checks
//! topology and would shrug at a vertex that moved. The accuracy harness checks
//! distance to the surface and would shrug at a vertex that moved by an ulp. The
//! property suite generates fresh meshes every run and has no memory. So a change
//! that is topologically identical, geometrically indistinguishable and
//! statistically invisible — a reordered summation, a different rounding on one
//! axis, a case table entry rotated onto an equivalent one — passes everything
//! and reaches a consumer's golden data as a silent diff.
//!
//! This is the file that notices.
//!
//! # Regenerating
//!
//! ```bash
//! ISOMESH_BLESS=1 cargo test -p isomesh --lib golden
//! ```
//!
//! Then **read the diff** before committing it. A golden fixture regenerated
//! without looking is worse than no fixture, because it converts an alarm into a
//! rubber stamp.
//!
//! # Why the hash is hand-rolled
//!
//! FNV-1a, ten lines, no dependency. `std`'s `DefaultHasher` is explicitly
//! documented as unstable across Rust releases, so a fixture built on it would
//! break on a toolchain bump and say "the mesher changed" — the one thing a
//! regression fixture must never lie about.
//!
//! The hash consumes raw IEEE bit patterns rather than formatted decimals, so it
//! distinguishes `+0.0` from `-0.0` and never depends on float formatting. That
//! matters here: a sign flip on a zero coordinate is exactly what a reordered
//! summation produces, and it is the class of change T-004 exists to catch.
//!
//! # Why `f64` only
//!
//! [`Real`](crate::Real) exposes no bit access — `total_cmp` compares, it does
//! not extract — so a generic hash would need `to_bits` added to a sealed core
//! trait, which this ticket does not justify. `f64` is what every existing test
//! meshes in, and it is the width where a change in the *algorithm* shows up
//! undisguised by rounding. E-112 is the ticket that measures `f32`'s behaviour
//! specifically.
//!
//! # What portability this asserts
//!
//! CI runs Linux and macOS, and these hashes are committed once. The crate routes
//! every transcendental through `libm` rather than `std` precisely so that they
//! agree — that decision is recorded in `CLAUDE.md` with this file as its
//! intended proof. If a platform ever disagrees, that is the finding.

use alloc::format;
use alloc::string::String;
use alloc::vec::Vec;

use crate::dual_contouring::DualContouring;
use crate::fields::ReferenceField;
use crate::greedy_quads::GreedyQuads;
use crate::manifold_dual_contouring::ManifoldDualContouring;
use crate::marching_cubes::{FaceAmbiguity, MarchingCubes};
use crate::marching_tetrahedra::MarchingTetrahedra;
use crate::surface_nets::SurfaceNets;
use crate::{MeshBuffer, RuntimeShape3, Sdf};

/// Resolutions every field is hashed at, in samples per axis.
///
/// Three, spanning coarse to ordinary. A single resolution would miss a change
/// that only fires on a case the coarse grid never reaches.
const RESOLUTIONS: [u32; 3] = [17, 25, 33];

/// The hash lives in [`crate::validate::mesh_hash`], not here.
///
/// It was written here first, and moved out when `benches/experiment_p38.rs`
/// needed to answer the same question — *did a change to a private buffer's
/// layout move any output?* — because two copies of a hash are two answers to
/// one question. The rationale for FNV-1a, for hashing bit patterns rather than
/// values, and for `f64` only, moved with it.
use crate::validate::mesh_hash;

/// The extractors, by name.
///
/// `marching_cubes+decider` is `MarchingCubes` with
/// [`FaceAmbiguity::AsymptoticDecider`] rather than a separate type, so it is a
/// fourth *row* here and not a fourth algorithm. Its hashes match plain
/// `marching_cubes`'s on every field where no ambiguous face occurs, which on the
/// reference set is five of the seven — that agreement is itself pinned by the
/// fixture.
#[derive(Clone, Copy, PartialEq, Eq)]
enum Algorithm {
    GreedyQuads,
    MarchingCubes,
    MarchingCubes33,
    /// Marching Cubes 33 with the **interior** rule as well as the face one —
    /// A-002b's, and the only row that meshes a tunnel as a tunnel.
    MarchingCubes33Trilinear,
    MarchingTetrahedra,
    SurfaceNets,
    DualContouring,
    ManifoldDualContouring,
    /// Subgrid Marching Tetrahedra, at a **fixed** 1D sampling resolution.
    ///
    /// M-95 is why the resolution is nailed down rather than derived from the
    /// grid: raising `samples` leaves the topology identical and moves the
    /// positions by around `1e-12`, because bisection converges to *an* ulp of a
    /// root and which one depends on the bracket it started from. A hash over
    /// this extractor is therefore only reproducible if the sampling is part of
    /// the fixture's definition, and [`SUBGRID_SAMPLES`] is that definition.
    SubgridMarchingTetrahedra,
}

/// The 1D sampling resolution the subgrid golden hashes are taken at.
///
/// Changing this changes every subgrid row in the fixture, and should be treated
/// like changing a resolution: deliberate, and re-baselined in the same commit.
const SUBGRID_SAMPLES: u32 = 16;

impl Algorithm {
    const ALL: [Self; 9] = [
        Self::GreedyQuads,
        Self::MarchingCubes,
        Self::MarchingCubes33,
        Self::MarchingCubes33Trilinear,
        Self::MarchingTetrahedra,
        Self::SurfaceNets,
        Self::DualContouring,
        Self::ManifoldDualContouring,
        Self::SubgridMarchingTetrahedra,
    ];

    fn name(self) -> &'static str {
        match self {
            Self::GreedyQuads => "greedy_quads",
            Self::MarchingCubes => "marching_cubes",
            Self::MarchingCubes33 => "marching_cubes+decider",
            Self::MarchingCubes33Trilinear => "marching_cubes+trilinear",
            Self::MarchingTetrahedra => "marching_tetrahedra",
            Self::SurfaceNets => "surface_nets",
            Self::DualContouring => "dual_contouring",
            Self::ManifoldDualContouring => "manifold_dual_contouring",
            Self::SubgridMarchingTetrahedra => "subgrid_marching_tetrahedra",
        }
    }
}

/// One row of the fixture.
struct Entry {
    algorithm: &'static str,
    field: &'static str,
    samples: u32,
    vertices: usize,
    triangles: usize,
    hash: u64,
}

fn extract<F>(algorithm: Algorithm, field: &F, samples: u32) -> MeshBuffer<f64>
where
    F: Sdf<Scalar = f64> + ReferenceField,
{
    let (lo, hi) = field.domain();
    let cell_size = (hi[0] - lo[0]) / f64::from(samples - 1);
    let shape = RuntimeShape3::new([samples; 3]).expect("valid shape");
    let mut out = MeshBuffer::<f64>::new();
    match algorithm {
        Algorithm::GreedyQuads => GreedyQuads::<f64>::new()
            .extract(field, &shape, lo, cell_size, &mut out)
            .expect("extraction"),
        Algorithm::MarchingCubes => MarchingCubes::<f64>::new()
            .extract(field, &shape, lo, cell_size, &mut out)
            .expect("extraction"),
        Algorithm::MarchingCubes33 => {
            let mut mc = MarchingCubes::<f64>::new();
            mc.set_face_ambiguity(FaceAmbiguity::AsymptoticDecider);
            mc.extract(field, &shape, lo, cell_size, &mut out)
                .expect("extraction");
        }
        Algorithm::MarchingCubes33Trilinear => {
            let mut mc = MarchingCubes::<f64>::new();
            mc.set_face_ambiguity(FaceAmbiguity::AsymptoticDecider);
            mc.set_interior_ambiguity(crate::marching_cubes::InteriorAmbiguity::Trilinear);
            mc.extract(field, &shape, lo, cell_size, &mut out)
                .expect("extraction");
        }
        Algorithm::MarchingTetrahedra => MarchingTetrahedra::<f64>::new()
            .extract(field, &shape, lo, cell_size, &mut out)
            .expect("extraction"),
        Algorithm::SurfaceNets => SurfaceNets::<f64>::new()
            .extract(field, &shape, lo, cell_size, &mut out)
            .expect("extraction"),
        Algorithm::DualContouring => DualContouring::<f64>::new()
            .extract(field, &shape, lo, cell_size, &mut out)
            .expect("extraction"),
        Algorithm::ManifoldDualContouring => ManifoldDualContouring::<f64>::new()
            .extract(field, &shape, lo, cell_size, &mut out)
            .expect("extraction"),
        Algorithm::SubgridMarchingTetrahedra => {
            crate::subgrid::extract::SubgridMarchingTetrahedra::<f64>::new(SUBGRID_SAMPLES)
                .expect("a positive sampling resolution")
                .extract(field, &shape, lo, cell_size, &mut out)
                .expect("extraction");
        }
    }
    out
}

/// Every combination, in a fixed order: field (as the reference sweep declares
/// them), then algorithm, then resolution.
fn compute_all() -> Vec<Entry> {
    let mut entries = Vec::new();
    crate::for_each_reference_field!(f64, |name, field| {
        for algorithm in Algorithm::ALL {
            for samples in RESOLUTIONS {
                let mesh = extract(algorithm, &field, samples);
                entries.push(Entry {
                    algorithm: algorithm.name(),
                    field: name,
                    samples,
                    vertices: mesh.vertex_count(),
                    triangles: mesh.triangle_count(),
                    hash: mesh_hash(&mesh),
                });
            }
        }
    });
    entries
}

fn fixture_path() -> std::path::PathBuf {
    std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("golden_hashes.json")
}

fn render(entries: &[Entry]) -> String {
    let mut out = String::from("[\n");
    for (i, e) in entries.iter().enumerate() {
        let comma = if i + 1 == entries.len() { "" } else { "," };
        out.push_str(&format!(
            "  {{\"algorithm\":\"{}\",\"field\":\"{}\",\"samples\":{},\"vertices\":{},\"triangles\":{},\"hash\":\"{:016x}\"}}{comma}\n",
            e.algorithm, e.field, e.samples, e.vertices, e.triangles, e.hash
        ));
    }
    out.push_str("]\n");
    out
}

/// Pull one value out of a fixture line.
///
/// A hand-rolled scanner rather than a JSON parser, because a JSON parser is a
/// dependency and this file is written by [`render`] two functions above — the
/// grammar is one line, fixed key order, no nesting and no escapes.
fn field_of<'a>(line: &'a str, key: &str) -> Option<&'a str> {
    let needle = format!("\"{key}\":");
    let start = line.find(&needle)? + needle.len();
    let rest = &line[start..];
    let rest = rest.strip_prefix('"').unwrap_or(rest);
    let end = rest.find(['"', ',', '}'])?;
    Some(&rest[..end])
}

#[cfg(test)]
mod tests;
