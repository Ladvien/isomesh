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

use crate::dc::DualContouring;
use crate::fields::ReferenceField;
use crate::mc::MarchingCubes;
use crate::sn::SurfaceNets;
use crate::{MeshBuffer, RuntimeShape3, Sdf};

/// Resolutions every field is hashed at, in samples per axis.
///
/// Three, spanning coarse to ordinary. A single resolution would miss a change
/// that only fires on a case the coarse grid never reaches.
const RESOLUTIONS: [u32; 3] = [17, 25, 33];

/// FNV-1a, 64-bit. Specified rather than chosen: the constants are the published
/// ones and the algorithm is four lines, so a reader can confirm it rather than
/// trust it.
struct Fnv(u64);

impl Fnv {
    const OFFSET_BASIS: u64 = 0xcbf2_9ce4_8422_2325;
    const PRIME: u64 = 0x0000_0100_0000_01b3;

    fn new() -> Self {
        Self(Self::OFFSET_BASIS)
    }

    fn write(&mut self, bytes: &[u8]) {
        for &b in bytes {
            self.0 ^= u64::from(b);
            self.0 = self.0.wrapping_mul(Self::PRIME);
        }
    }

    fn write_u64(&mut self, v: u64) {
        self.write(&v.to_le_bytes());
    }

    fn write_f64(&mut self, v: f64) {
        // Bits, not the value: `+0.0 == -0.0` compares equal but hashes
        // differently, which is the distinction T-004 cares about.
        self.write_u64(v.to_bits());
    }

    fn finish(self) -> u64 {
        self.0
    }
}

/// Hash a mesh's exact contents.
///
/// Counts are hashed first, so two meshes cannot collide by one being a prefix
/// of the other — a truncated index buffer is a defect this must catch.
fn hash_mesh(mesh: &MeshBuffer<f64>) -> u64 {
    let mut h = Fnv::new();
    h.write_u64(mesh.positions.len() as u64);
    h.write_u64(mesh.normals.len() as u64);
    h.write_u64(mesh.indices.len() as u64);
    for p in &mesh.positions {
        for v in p {
            h.write_f64(*v);
        }
    }
    for n in &mesh.normals {
        for v in n {
            h.write_f64(*v);
        }
    }
    for i in &mesh.indices {
        h.write_u64(u64::from(*i));
    }
    h.finish()
}

/// The three extractors, by name.
#[derive(Clone, Copy, PartialEq, Eq)]
enum Algorithm {
    MarchingCubes,
    SurfaceNets,
    DualContouring,
}

impl Algorithm {
    const ALL: [Self; 3] = [Self::MarchingCubes, Self::SurfaceNets, Self::DualContouring];

    fn name(self) -> &'static str {
        match self {
            Self::MarchingCubes => "mc",
            Self::SurfaceNets => "sn",
            Self::DualContouring => "dc",
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
        Algorithm::MarchingCubes => MarchingCubes::<f64>::new()
            .extract(field, &shape, lo, cell_size, &mut out)
            .expect("extraction"),
        Algorithm::SurfaceNets => SurfaceNets::<f64>::new()
            .extract(field, &shape, lo, cell_size, &mut out)
            .expect("extraction"),
        Algorithm::DualContouring => DualContouring::<f64>::new()
            .extract(field, &shape, lo, cell_size, &mut out)
            .expect("extraction"),
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
                    hash: hash_mesh(&mesh),
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
