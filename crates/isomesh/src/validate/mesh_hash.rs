//! Did the output change at all?
//!
//! Every other report here checks a *property*. The validity harness checks
//! topology and would shrug at a vertex that moved; the accuracy harness checks
//! distance to the surface and would shrug at a vertex that moved by an ulp; the
//! property suite generates fresh meshes every run and has no memory. So a change
//! that is topologically identical, geometrically indistinguishable and
//! statistically invisible — a reordered summation, a different rounding on one
//! axis, a case table entry rotated onto an equivalent one — passes everything
//! and reaches a consumer's golden data as a silent diff.
//!
//! This is the function that notices, and it is public because a consumer with
//! committed golden data needs exactly the same answer this crate's own T-007
//! fixture needs. `golden.rs` calls it; so does `benches/experiment_p38.rs`,
//! which has to show that a change to a *private buffer's layout* moved no
//! output. Two copies of a hash are two answers to "did it change", so there is
//! one.
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
//! undisguised by rounding.

use crate::MeshBuffer;

/// FNV-1a, 64-bit. Specified rather than chosen: the constants are the published
/// ones and the algorithm is four lines, so a reader can confirm it rather than
/// trust it.
struct Fnv(u64);

impl Fnv {
    const OFFSET_BASIS: u64 = 0xcbf2_9ce4_8422_2325;
    const PRIME: u64 = 0x0000_0100_0000_01b3;

    const fn new() -> Self {
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

    const fn finish(self) -> u64 {
        self.0
    }
}

/// Hash a mesh's exact contents.
///
/// Counts are hashed first, so two meshes cannot collide by one being a prefix
/// of the other — a truncated index buffer is a defect this must catch.
///
/// ```
/// use isomesh::fields::Sphere;
/// use isomesh::marching_cubes::MarchingCubes;
/// use isomesh::validate::mesh_hash;
/// use isomesh::{MeshBuffer, RuntimeShape3};
///
/// let shape = RuntimeShape3::new([17; 3])?;
/// let mut out = MeshBuffer::<f64>::new();
/// MarchingCubes::<f64>::new().extract(&Sphere::canonical(), &shape, [-2.0; 3], 0.25, &mut out)?;
///
/// // The same mesh hashes the same way, and that is the whole contract.
/// assert_eq!(mesh_hash(&out), mesh_hash(&out.clone()));
/// # Ok::<(), isomesh::Error>(())
/// ```
#[must_use]
pub fn mesh_hash(mesh: &MeshBuffer<f64>) -> u64 {
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
