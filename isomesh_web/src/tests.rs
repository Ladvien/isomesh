//! What these tests are for, and it is not coverage.
//!
//! This crate is a *label list plus a dispatch table* over someone else's
//! registry. The failure mode that matters is not a wrong triangle -- the golden
//! hashes in `crates/isomesh` own that -- it is the page offering `gyroid` and
//! meshing `fbm_terrain`, or a seventh extractor landing in the core crate and
//! silently never appearing in the dropdown. Both are invisible in a browser and
//! invisible in a build that succeeds. So the load-bearing tests here hold this
//! crate's two `const` arrays against `isomesh`'s own registries, and the rest
//! pin the exports' contract at their edges.

use std::sync::{Mutex, PoisonError};

use isomesh::extractor::{ALL_EXTRACTORS, UNREGISTERED};
use isomesh::for_each_reference_field;

use super::*;

/// Serialises the tests that read [`DEMO`] after writing it.
///
/// The exports are global by construction -- that is what a wasm module is -- so
/// two `#[test]` threads calling [`iso_mesh`] interleave and each sees the
/// other's mesh. Poisoning is stepped over rather than unwrapped: a failing test
/// holding this lock should fail once, not turn every sibling into a panic that
/// names the wrong problem.
static SERIAL: Mutex<()> = Mutex::new(());

#[test]
fn the_field_names_are_the_registry_in_order() {
    let mut registry = Vec::new();
    for_each_reference_field!(f32, |name, field| {
        // `field` is bound by the macro and every arm must use it, or the
        // expansion warns on eight unused bindings.
        let _ = field.domain();
        registry.push(name);
    });
    assert_eq!(
        registry, FIELD_NAMES,
        "the dropdown's labels and their order must be the crate's, not a copy \
         of it that has drifted"
    );
}

#[test]
fn the_extractor_names_partition_the_registry() {
    let mut offered: Vec<&str> = EXTRACTOR_NAMES.into_iter().chain(NOT_OFFERED).collect();
    offered.sort_unstable();
    let mut registry: Vec<&str> = ALL_EXTRACTORS.into_iter().collect();
    registry.sort_unstable();
    assert_eq!(
        offered, registry,
        "every registry entry is either offered on the page or named in \
         NOT_OFFERED with a reason -- a new extractor must be a decision, not an \
         omission"
    );
    assert_eq!(
        EXTRACTOR_NAMES.len() + NOT_OFFERED.len(),
        ALL_EXTRACTORS.len()
    );
}

#[test]
fn greedy_quads_is_excluded_by_the_crate_not_by_this_demo() {
    assert!(UNREGISTERED.contains(&"GreedyQuads"));
    assert!(
        !ALL_EXTRACTORS.contains(&"greedy_quads"),
        "if greedy quads ever joins the registry, NOT_OFFERED's comment about it \
         stops being true and this crate has to say why it is still absent"
    );
}

#[test]
fn resolve_field_covers_exactly_the_named_fields() {
    for index in 0..count(FIELD_NAMES.len()) {
        let resolved = resolve_field(index);
        assert!(
            resolved.is_some(),
            "field {index} is named but not resolvable"
        );
        let (_, lo, hi) = resolved.expect("checked on the line above");
        for axis in 0..3 {
            assert!(
                hi[axis] > lo[axis],
                "field {index} has an empty domain on axis {axis}"
            );
        }
    }
    assert!(resolve_field(count(FIELD_NAMES.len())).is_none());
}

#[test]
fn the_name_tables_end_where_the_counts_say() {
    for (kind, len) in [
        (KIND_FIELD, iso_field_count()),
        (KIND_EXTRACTOR, iso_extractor_count()),
    ] {
        for index in 0..len {
            assert!(!iso_name(kind, index).is_null());
            assert!(iso_name_len(kind, index) > 0);
        }
        assert!(iso_name(kind, len).is_null());
        assert_eq!(iso_name_len(kind, len), 0);
    }
    // Not a field table and not an extractor table, so there is no table.
    assert!(iso_name(2, 0).is_null());
    assert_eq!(iso_name_len(2, 0), 0);
}

#[test]
fn every_field_meshes_with_every_extractor() {
    let _serial = SERIAL.lock().unwrap_or_else(PoisonError::into_inner);
    for field in 0..iso_field_count() {
        for extractor in 0..iso_extractor_count() {
            let triangles = iso_mesh(field, extractor, 17);
            assert!(
                triangles > 0,
                "{} × {} meshed to nothing",
                FIELD_NAMES[field as usize],
                EXTRACTOR_NAMES[extractor as usize]
            );
            assert_eq!(iso_index_count(), triangles * 3);
            assert!(iso_vertex_count() > 0);
            assert!(iso_extent() > 0.0);
            assert!(!iso_positions().is_null());
            assert!(!iso_normals().is_null());
            assert!(!iso_indices().is_null());
        }
    }
}

#[test]
fn a_refusal_clears_the_mesh_and_the_report() {
    let _serial = SERIAL.lock().unwrap_or_else(PoisonError::into_inner);
    assert!(iso_mesh(0, 0, 17) > 0, "the sphere is the baseline");
    assert!(iso_vertex_count() > 0);

    for (field, extractor) in [(iso_field_count(), 0), (0, iso_extractor_count())] {
        assert_eq!(iso_mesh(field, extractor, 17), 0);
        assert_eq!(iso_vertex_count(), 0, "a refusal leaves no mesh behind");
        assert_eq!(iso_index_count(), 0);
        assert_eq!(iso_euler(), 0, "nor a previous run's report");
        assert_eq!(iso_non_manifold_edges(), 0);
        assert_eq!(iso_boundary_edges(), 0);
        assert_eq!(iso_degenerate_triangles(), 0);
    }
}

#[test]
fn samples_are_clamped_rather_than_refused() {
    let _serial = SERIAL.lock().unwrap_or_else(PoisonError::into_inner);
    let floor = iso_mesh(0, 0, 0);
    assert!(floor > 0, "a slider at zero means the narrowest grid");
    let at_min = iso_mesh(0, 0, MIN_SAMPLES);
    assert_eq!(floor, at_min);

    let ceiling = iso_mesh(0, 0, u32::MAX);
    assert!(ceiling > 0, "a slider past the end means the widest grid");
    let at_max = iso_mesh(0, 0, MAX_SAMPLES);
    assert_eq!(ceiling, at_max);
    assert!(
        at_max > at_min,
        "the widest grid must produce more triangles than the narrowest, or the \
         clamp is collapsing the range"
    );
}

#[test]
fn the_sphere_is_a_sphere_and_the_torus_is_a_torus() {
    let _serial = SERIAL.lock().unwrap_or_else(PoisonError::into_inner);
    // Marching Cubes on a closed field at a resolution that resolves it: χ is
    // the field's analytic value, and the HUD on the page is reading exactly
    // this. A number that is merely plausible would not catch a field/extractor
    // index swap; these two would.
    assert!(iso_mesh(0, 0, 33) > 0);
    assert_eq!(iso_euler(), 2, "a sphere");
    assert_eq!(iso_non_manifold_edges(), 0);
    assert_eq!(iso_boundary_edges(), 0);

    assert!(iso_mesh(1, 0, 33) > 0);
    assert_eq!(iso_euler(), 0, "a torus");
    assert_eq!(iso_non_manifold_edges(), 0);
    assert_eq!(iso_boundary_edges(), 0);
}
