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
