//! The probe, and the second question it cannot answer by itself.

use super::probe_mesh_shaders as probe;

/// Print what this machine says, and assert only what is true everywhere.
///
/// Deliberately **not** asserting that mesh shaders are present: the answer is
/// a property of the machine, and a test that demanded `true` would fail on
/// exactly the hardware the ticket wants reported on. What is asserted is that
/// the probe *ran* and that its two claims are consistent — `usable` without
/// `advertised` would mean the probe is reporting a device it never opened.
#[test]
fn the_probe_reports_every_adapter() {
    let reports = probe();
    for report in &reports {
        println!("{}", report.summary());
    }
    assert!(
        !reports.is_empty(),
        "no adapter at all -- this machine cannot answer the question"
    );
    for report in &reports {
        assert!(!report.name.is_empty(), "an adapter came back with no name");
        // The sub-features are refinements of the base one, so advertising a
        // refinement without it would mean the bits were read wrong.
        assert!(
            report.advertised || !(report.multiview || report.points),
            "{}: advertises a mesh-shader sub-feature without mesh shaders",
            report.name
        );
    }
}

/// What naga 29 will and will not accept, checked against the toolchain rather
/// than assumed from a doc.
///
/// The architecture doc records that WGSL mesh-shader *"frontend parsing is
/// done"*. That is right, and the first version of this test still failed —
/// because it wrote `@mesh` where naga requires **`@mesh(<global>)`**, naming a
/// mesh output variable that must live in the `workgroup` address space
/// (`front/wgsl/parse/mod.rs:1922` and `valid/interface.rs:1531`).
///
/// So this asserts the two things that can be stated exactly, and stops short
/// of a whole mesh shader: writing one means deriving the output struct naga
/// infers, and inventing that is GPU-008's job to do properly rather than this
/// probe's to guess.
#[test]
fn naga_implements_the_wgsl_mesh_shader_extension() {
    // The enable directive alone. If the extension were merely *named* and not
    // implemented, naga reports it as unimplemented and this fails -- which is
    // the distinction the doc's "parsing is done" claim rests on.
    let enabled = "enable wgpu_mesh_shader;\n";
    assert!(
        naga::front::wgsl::parse_str(enabled).is_ok(),
        "naga 29 does not implement the wgpu_mesh_shader enable extension"
    );

    // A task entry point, which takes no attribute arguments. Its dispatch is
    // a **builtin output** -- `@builtin(mesh_task_size)`,
    // `crate::BuiltIn::MeshTaskSize` at `parse/conv.rs:128` -- and not a
    // function call. The first version of this test called a
    // `dispatchMeshWorkgroups` that does not exist, which is what asking the
    // toolchain rather than recalling an API is for.
    let task = r"
enable wgpu_mesh_shader;

@task @workgroup_size(1)
fn task_main() -> @builtin(mesh_task_size) vec3<u32> {
    return vec3<u32>(1u, 1u, 1u);
}
";
    match naga::front::wgsl::parse_str(task) {
        Ok(_) => println!("naga 29 parses a @task entry point"),
        Err(why) => panic!("naga 29 rejects @task: {}", why.emit_to_string(task)),
    }

    // And the shape of the requirement, recorded by asserting the failure: a
    // bare `@mesh` is not enough.
    let bare_mesh = r"
enable wgpu_mesh_shader;

@mesh @workgroup_size(1)
fn mesh_main() {}
";
    assert!(
        naga::front::wgsl::parse_str(bare_mesh).is_err(),
        "`@mesh` without an output global now parses -- the note above is stale"
    );
}
