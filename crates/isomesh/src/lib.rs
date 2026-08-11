//! Engine-agnostic isosurface extraction. Field in, triangles out.
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

pub mod fields;

mod mesh;
mod real;
mod sdf;
mod shape;

pub use mesh::{MeshBuffer, MeshSink};
pub use real::Real;
pub use sdf::Sdf;
pub use shape::{ConstShape3, RuntimeShape3, Shape3};
