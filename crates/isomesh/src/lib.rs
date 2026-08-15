//! Engine-agnostic isosurface extraction. Field in, triangles out.
//!
//! ```
//! use isomesh::marching_cubes::MarchingCubes;
//! use isomesh::fields::Sphere;
//! use isomesh::{MeshBuffer, RuntimeShape3};
//!
//! let field = Sphere::<f32>::canonical();          // any `Sdf` -- yours, or one of eight here
//! let shape = RuntimeShape3::new([33; 3])?;        // 33 samples per axis, so 32 cells
//! let mut mesh = MeshBuffer::<f32>::new();         // reused across calls; never reallocated
//!
//! MarchingCubes::<f32>::new()
//!     .extract(&field, &shape, [-2.0; 3], 0.125, &mut mesh)?;
//!
//! assert!(mesh.triangle_count() > 0);
//! # Ok::<(), isomesh::Error>(())
//! ```
//!
//! `mesh.positions` is `Vec<[f32; 3]>` and `mesh.indices` is `Vec<u32>` — hand
//! them to a renderer, a physics engine, or an exporter. There is no mesh type to
//! learn and no math library in the way.
//!
//! # Why it looks like this
//!
//! `isomesh` has to serve both a real-time voxel game and a CAD tool, and that
//! single constraint decides most of its design: no math library appears in a
//! public signature, output buffers are caller-provided and reusable, and the
//! scalar type is generic over `f32` and `f64`.
//!
//! # Conventions
//!
//! These hold for every type and every algorithm in this crate. They are stated
//! once here and repeated on the items they constrain, because each one is a
//! convention rather than a fact, and a mismatch across a module boundary
//! produces output that looks plausible and is wrong.
//!
//! - **Sign.** *Negative is inside.* A field is negative strictly inside the
//!   solid, positive strictly outside, and the surface is the zero level set.
//! - **Handedness.** Right-handed. `x × y = +z`.
//! - **Winding.** Counter-clockwise viewed from outside the solid.
//! - **Normals.** Point away from the solid, i.e. along increasing field value.
//! - **Index order.** `x` varies fastest: `i = x + y·sx + z·sx·sy`.
//!
//! # Float math
//!
//! Every transcendental goes through [`libm`], unconditionally — there is no
//! `std` fast path. Two reasons. One execution path: a `#[cfg(feature = "std")]`
//! fork would mean two float backends and two sets of results for the same
//! input. And determinism: `std`'s `sin`/`cos` are the platform's, are not
//! correctly rounded, and differ between macOS and Linux, which would make
//! committed golden hashes platform-specific. It costs nothing at run time —
//! `libm::sqrtf` compiles to `fsqrt` on aarch64+neon and `sqrtss` on
//! x86-64+sse2.
//!
//! # `no_std`
//!
//! This crate is `no_std` unconditionally, not `cfg_attr`-conditionally. A
//! conditional attribute would put `std`'s prelude in scope for the default
//! build, so a stray `use std::collections::HashMap` would compile here and
//! break only in some later `no_std` build. Unconditional makes that a compile
//! error at the moment it is typed.

#![no_std]

extern crate alloc;

// The test harness itself needs `std`; the library never does.
#[cfg(test)]
extern crate std;

pub mod brush;
pub mod chunk;
pub mod collider;
pub mod dual_contouring;
/// Speculative algorithms, behind the `experimental` feature.
///
/// Off by default, exempt from semver, and exempt from nothing else — see the
/// module's own docs.
#[cfg(feature = "experimental")]
pub mod experimental;
pub mod extractor;
pub mod fields;
pub mod greedy_quads;
pub mod hermite;
pub mod manifold_dual_contouring;
pub mod marching_cubes;
pub mod marching_tetrahedra;
pub mod normals;
pub mod orient;
pub mod paint;
pub mod subgrid;
pub mod surface_nets;
pub mod transvoxel;
pub mod validate;
pub mod weld;

#[cfg(test)]
mod golden;
#[cfg(test)]
mod property;

mod cube;
pub mod dual;
mod equivariant;
mod error;
mod mesh;
mod real;
mod sdf;
mod shape;
mod vec3;

pub use error::{Error, Result};
pub use mesh::{MeshBuffer, MeshSink};
pub use real::Real;
pub use sdf::Sdf;
pub use shape::{ConstShape3, RuntimeShape3, Shape3};

/// Compiles the README's example as a doctest, without putting the README into
/// these docs.
///
/// A code sample on a crate's landing page is the first thing a reader tries and
/// the first thing to rot, because nothing builds it. `cfg(doctest)` means this
/// item exists only when rustdoc is collecting doctests, so the example is
/// checked on every `cargo test` and appears in no rendered documentation.
#[cfg(doctest)]
#[doc = include_str!("../README.md")]
pub struct ReadmeExample;
