//! Most of this needs no GPU. The last test does, because "the WGSL is valid"
//! is not a claim a string comparison can make.

use super::{Composer, FEATURES};
use crate::Error;

fn composer(modules: &[(&str, &'static str)]) -> Composer {
    let mut composer = Composer::new();
    for (name, source) in modules {
        composer.insert(name, source);
    }
    composer
}

#[test]
fn text_without_directives_passes_through() {
    let c = composer(&[("a", "fn f() {}\n  indented\n")]);
    assert_eq!(
        c.compose("a", &[]).as_deref(),
        Ok("fn f() {}\n  indented\n")
    );
}

#[test]
fn an_include_pastes_the_module_in_place() {
    let c = composer(&[
        ("dep", "DEP\n"),
        ("root", "before\n#include <dep>\nafter\n"),
    ]);
    assert_eq!(
        c.compose("root", &[]).as_deref(),
        Ok("before\nDEP\nafter\n")
    );
}

/// The semantics WGSL forces: two modules can both depend on a shared header
/// and the header appears once. A duplicated function is a hard error in WGSL,
/// so the alternative is not a style preference.
#[test]
fn a_module_included_twice_appears_once() {
    let c = composer(&[
        ("shared", "SHARED\n"),
        ("left", "#include <shared>\nLEFT\n"),
        ("right", "#include <shared>\nRIGHT\n"),
        ("root", "#include <left>\n#include <right>\n"),
    ]);
    let out = c.compose("root", &[]).expect("composes");
    assert_eq!(out.matches("SHARED").count(), 1, "got:\n{out}");
    assert!(out.contains("LEFT") && out.contains("RIGHT"));
}

/// Include-once would quietly terminate a cycle. It is still an error, because
/// a cycle is a question about which module owns what and absorbing it hides
/// the question.
#[test]
fn a_cycle_is_an_error_rather_than_being_absorbed() {
    let c = composer(&[("a", "#include <b>\n"), ("b", "#include <a>\n")]);
    assert_eq!(
        c.compose("a", &[]).err(),
        Some(Error::ShaderCircularInclude {
            name: String::from("a")
        })
    );
}

#[test]
fn a_missing_module_is_named() {
    let c = composer(&[("root", "#include <absent>\n")]);
    assert_eq!(
        c.compose("root", &[]).err(),
        Some(Error::ShaderModuleMissing {
            name: String::from("absent")
        })
    );
    assert_eq!(
        c.compose("also_absent", &[]).err(),
        Some(Error::ShaderModuleMissing {
            name: String::from("also_absent")
        })
    );
}

#[test]
fn ifdef_and_ifndef_select_opposite_branches() {
    let c = composer(&[("root", "#ifdef X\nYES\n#else\nNO\n#endif\n")]);
    assert_eq!(c.compose("root", &["X"]).as_deref(), Ok("YES\n"));
    assert_eq!(c.compose("root", &[]).as_deref(), Ok("NO\n"));

    let c = composer(&[("root", "#ifndef X\nABSENT\n#else\nPRESENT\n#endif\n")]);
    assert_eq!(c.compose("root", &["X"]).as_deref(), Ok("PRESENT\n"));
    assert_eq!(c.compose("root", &[]).as_deref(), Ok("ABSENT\n"));
}

#[test]
fn regions_nest() {
    let source = "#ifdef A\nA1\n#ifdef B\nAB\n#else\nAnotB\n#endif\nA2\n#else\nnotA\n#endif\n";
    let c = composer(&[("root", source)]);
    assert_eq!(
        c.compose("root", &["A", "B"]).as_deref(),
        Ok("A1\nAB\nA2\n")
    );
    assert_eq!(c.compose("root", &["A"]).as_deref(), Ok("A1\nAnotB\nA2\n"));
    assert_eq!(c.compose("root", &["B"]).as_deref(), Ok("notA\n"));
    assert_eq!(c.compose("root", &[]).as_deref(), Ok("notA\n"));
}

/// A directive inside a region that is switched off must still be *parsed* for
/// nesting, or its `#endif` closes the wrong region. This is the bug a naive
/// "skip every line while disabled" implementation has.
#[test]
fn a_disabled_region_still_tracks_nesting() {
    let source = "#ifdef ON\nkept\n#ifdef OFF\ndropped\n#endif\nstill kept\n#endif\ntail\n";
    let c = composer(&[("root", source)]);
    assert_eq!(
        c.compose("root", &["ON"]).as_deref(),
        Ok("kept\nstill kept\ntail\n")
    );
    assert_eq!(c.compose("root", &[]).as_deref(), Ok("tail\n"));
}

/// An include inside a switched-off region must not be pasted — and must not be
/// marked as included either, or a later live include of the same module
/// silently produces nothing.
#[test]
fn an_include_in_a_dead_branch_is_neither_pasted_nor_consumed() {
    let c = composer(&[
        ("dep", "DEP\n"),
        (
            "root",
            "#ifdef NEVER\n#include <dep>\n#endif\n#include <dep>\n",
        ),
    ]);
    assert_eq!(c.compose("root", &[]).as_deref(), Ok("DEP\n"));
}

#[test]
fn unbalanced_directives_are_refused_with_a_line_number() {
    for (source, line) in [
        ("a\n#endif\n", 2),
        ("a\n#else\n", 2),
        ("#ifdef A\nbody\n", 3),
        ("#ifdef A\n#else\n#else\n#endif\n", 3),
        ("#ifdef\n#endif\n", 1),
        ("#ifdef A B\n#endif\n", 1),
    ] {
        assert_eq!(
            composer(&[("root", source)]).compose("root", &[]).err(),
            Some(Error::ShaderDirective {
                module: String::from("root"),
                line
            }),
            "accepted {source:?}"
        );
    }
}

#[test]
fn a_malformed_include_is_refused() {
    for source in ["#include dep\n", "#include <>\n", "#include <dep\n"] {
        assert!(
            matches!(
                composer(&[("dep", "D\n"), ("root", source)])
                    .compose("root", &[])
                    .err(),
                Some(Error::ShaderDirective { .. })
            ),
            "accepted {source:?}"
        );
    }
}

#[test]
fn composition_is_deterministic() {
    let c = Composer::with_builtins();
    assert_eq!(c.compose("grid", &[]), c.compose("grid", &[]));
    assert_eq!(c.module_names(), ["grid"]);
}

/// The builtin module carries the layout the CPU side packs, so the names the
/// rest of this crate will call are pinned here rather than discovered missing
/// at pipeline-creation time.
#[test]
fn the_grid_module_declares_what_the_cpu_side_packs() {
    let out = Composer::with_builtins()
        .compose("grid", &[])
        .expect("composes");
    for expected in [
        "struct GridParams",
        "samples: vec4<u32>",
        "placement: vec4<f32>",
        "fn grid_index",
        "fn grid_position",
        "fn grid_cells",
        "fn grid_contains",
    ] {
        assert!(out.contains(expected), "grid.wgsl is missing `{expected}`");
    }
}

/// "It is valid WGSL" is not a claim `contains` can make, so this one asks a
/// driver.
///
/// A no-GPU `naga` check over every permutation is GPU-003 and is the version
/// that belongs in CI. This is the version available today, and it is stronger
/// per-shader: it is the same path a real pipeline takes.
#[test]
fn the_grid_module_is_valid_wgsl_on_a_real_device() {
    let gpu =
        crate::headless::Gpu::new().expect("a GPU adapter -- no software fallback, by design");
    let source = Composer::with_builtins()
        .compose("grid", &[])
        .expect("composes");

    // A shader module has to be *used* before some backends validate it, and
    // grid.wgsl declares only functions and a struct. Give it an entry point
    // that calls every one of them, so nothing is dead-stripped before it is
    // checked.
    let probe = format!(
        "{source}
@group(0) @binding(0) var<uniform> params: GridParams;
@group(0) @binding(1) var<storage, read_write> out: array<f32>;

@compute @workgroup_size(1)
fn main(@builtin(global_invocation_id) id: vec3<u32>) {{
    if (!grid_contains(params, id)) {{ return; }}
    let p = grid_position(params, id);
    let cells = grid_cells(params);
    out[grid_index(params, id)] = p.x + p.y + p.z + f32(cells.x + grid_sample_count(params));
}}
"
    );

    // wgpu 29 hands back a guard rather than pairing push with a device method,
    // so the scope cannot be left open by an early return.
    let scope = gpu.device().push_error_scope(wgpu::ErrorFilter::Validation);
    let module = gpu
        .device()
        .create_shader_module(wgpu::ShaderModuleDescriptor {
            label: Some("grid.wgsl validation probe"),
            source: wgpu::ShaderSource::Wgsl(probe.into()),
        });
    let error = crate::block_on::block_on(scope.pop());
    assert!(error.is_none(), "grid.wgsl failed validation: {error:?}");
    drop(module);
}

// -- GPU-003: the validation sweep -------------------------------------------
//
// No GPU, no device, no display. This is the half that belongs in CI, and it
// catches the entire class of "compiles on my Vulkan driver, explodes on DX12"
// before an adapter is ever opened.

/// Parse and validate one composed source with the same `naga` `wgpu` carries.
///
/// Returns the error as a string rather than propagating it, because the
/// interesting output of a sweep is *which variant* failed and what it said,
/// and a `?` would stop at the first one.
fn validate(source: &str) -> Result<(), String> {
    let module = naga::front::wgsl::parse_str(source)
        .map_err(|e| format!("parse: {}", e.emit_to_string(source)))?;
    naga::valid::Validator::new(
        naga::valid::ValidationFlags::all(),
        naga::valid::Capabilities::all(),
    )
    .validate(&module)
    .map(|_| ())
    .map_err(|e| format!("validate: {e:?}"))
}

/// The cross product is `modules × 2^features`, in a stable order.
///
/// Asserted with synthetic features because the real list is empty today, and a
/// combinatorial harness that has only ever been run at n = 0 has not been run.
#[test]
fn the_variant_sweep_is_the_full_cross_product() {
    let c = composer(&[("a", "A\n"), ("b", "B\n")]);

    assert_eq!(c.variants(&[]).len(), 2, "two modules, no features");

    let variants = c.variants(&["X", "Y"]);
    assert_eq!(variants.len(), 2 * 4, "two modules x 2^2 subsets");

    // Stable order, and every subset present exactly once per module.
    let for_a: Vec<Vec<&str>> = variants
        .iter()
        .filter(|(m, _)| *m == "a")
        .map(|(_, d)| d.clone())
        .collect();
    assert_eq!(
        for_a,
        vec![vec![], vec!["X"], vec!["Y"], vec!["X", "Y"]],
        "subsets must be ascending bitmask order"
    );
    assert_eq!(c.variants(&["X", "Y"]), c.variants(&["X", "Y"]));
}

/// The sweep itself, over what this crate actually ships.
#[test]
fn every_shader_permutation_validates() {
    let composer = Composer::with_builtins();
    let variants = composer.variants(FEATURES);

    // M-44: a gate that has only ever passed is indistinguishable from one that
    // cannot fail, and a sweep over an empty set passes beautifully. Pin the
    // size against what it is derived from.
    assert_eq!(
        variants.len(),
        composer.module_names().len() * (1usize << FEATURES.len()),
        "the sweep is not the full cross product"
    );
    assert!(!variants.is_empty(), "the sweep covered nothing");

    let mut failures = Vec::new();
    for (module, defines) in &variants {
        match composer.compose(module, defines) {
            Ok(source) => {
                if let Err(why) = validate(&source) {
                    failures.push(format!("{module} {defines:?}: {why}"));
                }
            }
            Err(why) => failures.push(format!("{module} {defines:?}: compose: {why}")),
        }
    }
    assert!(
        failures.is_empty(),
        "{} of {} variants failed:\n{}",
        failures.len(),
        variants.len(),
        failures.join("\n")
    );
}

/// And the sweep can fail, which is the only way to know it is a gate.
///
/// M-44 again, applied to this test rather than to the code it guards: a
/// validator that accepts everything would make the test above green forever.
#[test]
fn the_validator_rejects_invalid_wgsl() {
    assert!(validate("fn f() -> f32 { return 1.0; }\n").is_ok());
    assert!(
        validate("fn f() -> f32 { return no_such_thing(); }\n").is_err(),
        "the validator accepted a call to an undeclared function"
    );
    assert!(
        validate("this is not wgsl\n").is_err(),
        "the validator accepted prose"
    );
}
