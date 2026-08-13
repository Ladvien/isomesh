//! GPU-007 — print what this machine says about mesh shaders, and stop.
//!
//! ```bash
//! cargo run -p isomesh-gpu --example mesh_shader_probe
//! ```
//!
//! No Bevy, no window, no shader. `CLAUDE.md` records two sources that
//! disagree about Metal and instructs that a capability probe report the truth
//! before anything is built on mesh shaders. This is that probe, and it is a
//! plain `wgpu` program on purpose: the question is about the adapter, and
//! answering it through a game engine would make the engine part of the answer.
//!
//! **It reports; it does not decide.** Every line below is either a bit read
//! off an adapter or a quotation with its source.

fn main() {
    let reports = isomesh_gpu::probe_mesh_shaders();

    println!("adapters on this machine:");
    if reports.is_empty() {
        println!("  (none -- this machine cannot answer the question)");
    }
    for report in &reports {
        println!("  {}", report.summary());
    }

    let any = reports.iter().any(|r| r.advertised);
    println!();
    println!("EXPERIMENTAL_MESH_SHADER advertised by at least one adapter: {any}");

    println!();
    println!("What the bits do not tell you, with sources:");
    println!();
    println!("  1. Enabling it needs `unsafe`. `ExperimentalFeatures::enabled()` is a");
    println!("     `const unsafe fn` (wgpu-types 29.0.4, src/tokens.rs), and its `disabled()`");
    println!("     counterpart is documented as \"uses of Features prefixed with EXPERIMENTAL");
    println!("     are disallowed\". This workspace sets `unsafe_code = \"forbid\"`, so no");
    println!("     device here can request it.");
    println!();
    println!("  2. WGSL mesh shaders are Vulkan-only. wgpu-types 29.0.4 on the feature:");
    println!("     \"Naga is only supported on vulkan. On other platforms you will have to");
    println!("     use passthrough shaders.\" So on Metal and DX12 the WGSL this crate");
    println!("     composes cannot be a mesh shader at all.");
    println!();
    println!("  3. naga 29 does implement the WGSL extension. `enable wgpu_mesh_shader;`");
    println!("     parses, `@task` parses, and `@mesh` requires `@mesh(<global>)` naming an");
    println!("     output variable in the `workgroup` address space (front/wgsl/parse/mod.rs");
    println!("     :1922, valid/interface.rs:1531).");
    println!();
    println!("  4. Metal is NOT answered by this run unless a Metal adapter appears above.");
    println!("     The macOS case is the one CLAUDE.md flags as unverified; run this there.");
}
