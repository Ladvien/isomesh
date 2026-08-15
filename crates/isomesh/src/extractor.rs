//! One trait over every extraction algorithm, and one list of them.
//!
//! # Why this exists
//!
//! Every algorithm in this crate already had the *same* signature — field,
//! shape, origin, cell size, caller-provided sink, `Result<()>` — and nothing
//! said so. The consequence was measured at X-001: the algorithms are named by
//! hand in `benches/shootout.rs` (26 references), `src/property/extraction.rs`
//! (12), `benches/resolution_sweep.rs` (9), `benches/extract.rs` (10) and
//! `benches/stage_breakdown.rs` (2). **Adding a ninth algorithm was an `O(N)`
//! edit across benches, property tests and examples**, and the failure mode is
//! not a compile error — it is a bench that silently measures six algorithms
//! while its header says seven.
//!
//! [`Extractor`] is that shared shape, and [`for_each_extractor!`](crate::for_each_extractor) is the single
//! list.
//!
//! # This is a shape, not a hierarchy
//!
//! The trait does one thing: name the call every algorithm already answers. It
//! deliberately has **no** `new`, no configuration and no name.
//!
//! **Configuration lives in the registry entry rather than in an associated
//! type**, which is a deviation from X-001 as written and is the simpler answer
//! to the case that motivated it. `marching_cubes+decider` is not a seventh
//! algorithm; it is [`MarchingCubes`] with
//! [`FaceAmbiguity::AsymptoticDecider`]
//! set, so a `const NAME` on the type cannot name it and an associated `Config`
//! would have to be threaded through every caller to reach it. A registry entry
//! is an expression, so it carries its own configuration and its own name:
//! Dual Contouring's Hermite settings, subgrid's sampling resolution and a
//! future Transvoxel LOD are each expressible without the trait knowing they
//! exist.
//!
//! That also leaves the seam X-002 wants: a variant is another entry, not
//! another branch inside one.
//!
//! # Why a macro and not a `Vec<Box<dyn Extractor>>`
//!
//! [`extract_into`](Extractor::extract_into) is generic over the field and the
//! sink, so the trait is not object-safe and could only be made so by fixing
//! both — which would force one scalar type and one sink on every caller, and
//! this crate is generic over `f32` and `f64` on purpose. The macro
//! monomorphises each entry instead, so a sweep pays no dispatch cost and keeps
//! the reused-buffer path rule 6 requires.
//!
//! It is also the pattern this crate already uses for the other list that has to
//! stay in one place —
//! [`for_each_reference_field!`](crate::for_each_reference_field) — so there is
//! one idiom for "enumerate the things the suite must cover", not two.

use crate::marching_cubes::{FaceAmbiguity, MarchingCubes};
use crate::mesh::MeshSink;
use crate::real::Real;
use crate::sdf::Sdf;
use crate::shape::Shape3;

/// The call every extraction algorithm in this crate answers.
///
/// `field` is sampled on a grid of `shape` samples whose lowest corner sits at
/// `origin`, spaced `cell_size` apart, and the triangles are written into `out`.
///
/// # The signature is not negotiable
///
/// `out` is caller-provided and reusable, per `CLAUDE.md` rule 6 — the real
/// workload re-meshes thousands of chunks per edit, so an algorithm that
/// returned a freshly allocated mesh would allocate on every one of them. And
/// nothing here names a math library, per rule 1: `[R; 3]` rather than any
/// `Vec3`, so a consumer using two crates with incompatible `glam` versions can
/// still call this.
pub trait Extractor<R: Real> {
    /// Extract the zero level set of `field` into `out`.
    ///
    /// # Errors
    ///
    /// Whatever the underlying algorithm reports — a grid too small to hold a
    /// cell, an index space exhausted, or a configuration the published
    /// construction does not define. Errors are returned rather than absorbed;
    /// no implementation substitutes a degraded mesh.
    fn extract_into<S, M>(
        &mut self,
        field: &S,
        shape: &impl Shape3,
        origin: [R; 3],
        cell_size: R,
        out: &mut M,
    ) -> crate::Result<()>
    where
        S: Sdf<Scalar = R>,
        M: MeshSink<Scalar = R>;
}

/// Implement [`Extractor`] by forwarding to the inherent `extract`.
///
/// Every algorithm's inherent method already has this exact signature, so each
/// impl is a forward and nothing else. Written as a macro so that seven
/// identical bodies cannot drift apart from each other.
macro_rules! forward_extractor {
    ($($ty:path),+ $(,)?) => {
        $(
            impl<R: Real> Extractor<R> for $ty {
                fn extract_into<S, M>(
                    &mut self,
                    field: &S,
                    shape: &impl Shape3,
                    origin: [R; 3],
                    cell_size: R,
                    out: &mut M,
                ) -> crate::Result<()>
                where
                    S: Sdf<Scalar = R>,
                    M: MeshSink<Scalar = R>,
                {
                    self.extract(field, shape, origin, cell_size, out)
                }
            }
        )+
    };
}

forward_extractor!(
    crate::marching_cubes::MarchingCubes<R>,
    crate::marching_tetrahedra::MarchingTetrahedra<R>,
    crate::surface_nets::SurfaceNets<R>,
    crate::dual_contouring::DualContouring<R>,
    crate::manifold_dual_contouring::ManifoldDualContouring<R>,
    crate::greedy_quads::GreedyQuads<R>,
    crate::subgrid::extract::SubgridMarchingTetrahedra<R>,
);

/// The sampling resolution the subgrid entry runs at.
///
/// Subgrid Marching Tetrahedra is the one entry whose constructor can fail, and
/// this is the value that makes it succeed. Named rather than inlined because
/// M-98's 70× figure is *about* this number — it is the `16` in
/// `6 tets × 6 edges × 16 samples = 576` field evaluations per cell — so a
/// change here moves a published measurement.
pub const SUBGRID_SAMPLES: u32 = 16;

/// Every extractor the suite covers, run once each.
///
/// Binds `$name` to the entry's label and `$extractor` to a configured, mutable
/// instance, then evaluates the body once per entry. **This is the single list**
/// that X-001 exists to create: a new algorithm is one line here and nothing
/// else.
///
/// ```
/// # use isomesh::{for_each_extractor, MeshBuffer, RuntimeShape3};
/// # use isomesh::extractor::Extractor;
/// # use isomesh::fields::Sphere;
/// let field = Sphere::<f64>::canonical();
/// let shape = RuntimeShape3::new([17; 3]).expect("valid shape");
/// let mut out = MeshBuffer::<f64>::new();
/// for_each_extractor!(f64, |name, extractor| {
///     out.reset();
///     extractor
///         .extract_into(&field, &shape, [-2.0; 3], 0.25, &mut out)
///         .expect("extraction");
///     assert!(out.triangle_count() > 0, "{name} meshed a sphere to nothing");
/// });
/// ```
///
/// # `marching_cubes+decider` is an entry, not an algorithm
///
/// It is the same type with a different face rule, which is why entries are
/// expressions rather than types. See this module's header.
#[macro_export]
macro_rules! for_each_extractor {
    ($scalar:ty, |$name:ident, $extractor:ident| $body:block) => {{
        {
            let $name = "marching_cubes";
            #[allow(unused_mut)]
            let mut $extractor = $crate::marching_cubes::MarchingCubes::<$scalar>::new();
            $body
        }
        {
            let $name = "marching_cubes+decider";
            #[allow(unused_mut)]
            let mut $extractor = {
                let mut mesher = $crate::marching_cubes::MarchingCubes::<$scalar>::new();
                mesher.set_face_ambiguity($crate::marching_cubes::FaceAmbiguity::AsymptoticDecider);
                mesher
            };
            $body
        }
        {
            let $name = "marching_tetrahedra";
            #[allow(unused_mut)]
            let mut $extractor = $crate::marching_tetrahedra::MarchingTetrahedra::<$scalar>::new();
            $body
        }
        {
            let $name = "surface_nets";
            #[allow(unused_mut)]
            let mut $extractor = $crate::surface_nets::SurfaceNets::<$scalar>::new();
            $body
        }
        {
            let $name = "dual_contouring";
            #[allow(unused_mut)]
            let mut $extractor = $crate::dual_contouring::DualContouring::<$scalar>::new();
            $body
        }
        {
            let $name = "manifold_dual_contouring";
            #[allow(unused_mut)]
            let mut $extractor =
                $crate::manifold_dual_contouring::ManifoldDualContouring::<$scalar>::new();
            $body
        }
        {
            let $name = "subgrid_marching_tetrahedra";
            #[allow(unused_mut)]
            let mut $extractor =
                $crate::subgrid::extract::SubgridMarchingTetrahedra::<$scalar>::new(
                    $crate::extractor::SUBGRID_SAMPLES,
                )
                .expect("SUBGRID_SAMPLES is a positive sampling resolution");
            $body
        }
    }};
}

/// The entries [`for_each_extractor!`](crate::for_each_extractor) visits, in order.
///
/// Kept beside the macro so a reader can see the list without expanding it, and
/// so `the_registry_and_the_macro_agree` can check the two against each other —
/// a name list that drifts from the thing it names is worse than no list.
pub const ALL_EXTRACTORS: [&str; 7] = [
    "marching_cubes",
    "marching_cubes+decider",
    "marching_tetrahedra",
    "surface_nets",
    "dual_contouring",
    "manifold_dual_contouring",
    "subgrid_marching_tetrahedra",
];

/// [`Extractor`] impls that are deliberately **not** in [`for_each_extractor!`](crate::for_each_extractor),
/// and why.
///
/// An impl missing from the registry is normally a mistake, so
/// `every_extractor_impl_is_registered_or_excused` fails on one. This is the
/// escape hatch, and it costs a name and a reason rather than silence.
///
/// - **`GreedyQuads`** — it fits the trait and is not an isosurface extractor.
///   It classifies whole cells as solid or empty and emits the axis-aligned
///   faces between them, so the output is *"a Minecraft surface rather than an
///   isosurface"*. Sweeping it beside the others would compare a blocky mesh's
///   Hausdorff error and manifoldness against surfaces that interpolate, which
///   is a category error rather than a measurement. It keeps the impl because a
///   caller may legitimately want to drive it through the same shape.
pub const UNREGISTERED: [&str; 1] = ["GreedyQuads"];

/// Silences the unused-import warning `FaceAmbiguity` and `MarchingCubes` would
/// otherwise raise: they are named only inside the exported macro, which expands
/// in the *caller's* crate and so does not count as a use here.
#[allow(dead_code)]
fn macro_paths_are_used(_: Option<FaceAmbiguity>, _: Option<MarchingCubes<f32>>) {}

#[cfg(test)]
mod tests;
