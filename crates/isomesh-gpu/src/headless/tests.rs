//! What this machine actually is, reported rather than assumed.

use super::Gpu;

/// Opening a device is the precondition for every other GPU test here, so it
/// gets its own test — otherwise a machine with no adapter fails five tests
/// with five confusing messages instead of one clear one.
#[test]
fn a_device_opens_without_a_window() {
    let gpu = Gpu::new().expect("a GPU adapter -- there is no software fallback, by design");
    let report = gpu.report();

    println!(
        "adapter: {} ({:?}, {:?}) driver {}",
        report.name, report.backend, report.device_type, report.driver
    );
    println!(
        "max storage binding {} bytes, max workgroups/dim {}",
        report.max_storage_buffer_binding_size, report.max_compute_workgroups_per_dimension
    );

    assert!(!report.name.is_empty());
    // Compute needs both of these to be non-zero to dispatch anything at all,
    // and GPU-004's grid sizing reads them.
    assert!(report.max_storage_buffer_binding_size > 0);
    assert!(report.max_compute_workgroups_per_dimension > 0);
}

/// The no-fallback promise in `Gpu::new`'s docs, asserted rather than trusted.
///
/// A CPU adapter here would mean every timing this crate ever reports is off by
/// orders of magnitude while looking merely slow.
#[test]
fn the_adapter_is_not_a_software_rasteriser() {
    let gpu = Gpu::new().expect("a GPU adapter");
    assert_ne!(
        gpu.report().device_type,
        wgpu::DeviceType::Cpu,
        "opened a software adapter: {}",
        gpu.report().name
    );
}

/// The coverage hole GPU-009 accepted, pinned so the documentation cannot rot
/// into a belief.
///
/// If a future wgpu stops gating experimental features behind an `unsafe`
/// token, or this crate's lint policy changes, this fails — and the table on
/// [`Gpu`] saying mesh shaders cannot run headlessly needs re-reading rather
/// than quietly staying wrong.
#[test]
fn the_headless_device_has_no_experimental_features() {
    let gpu = Gpu::new().expect("a GPU adapter");
    let features = gpu.device().features();

    assert!(
        !features.contains(wgpu::Features::EXPERIMENTAL_MESH_SHADER),
        "the headless device now has mesh shaders -- GPU-009 accepted a coverage \
         hole on the grounds that it cannot, and that reasoning needs re-reading"
    );
    // The adapter advertising it while the device lacks it is the whole shape
    // of the hole: the hardware can, and this crate will not ask.
    let advertised = crate::probe_mesh_shaders()
        .iter()
        .any(|report| report.advertised);
    if advertised {
        println!("adapter advertises mesh shaders; this device does not request them, by design");
    }
}
