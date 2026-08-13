//! What this machine actually says about mesh shaders.
//!
//! # The disagreement this exists to settle
//!
//! `CLAUDE.md` records two sources that contradict each other: *"wgpu's spec
//! table lists MSL as planned while the tracking issue says the Metal HAL
//! backend merged"*, and instructs that nothing be built on mesh shaders
//! before a capability probe reports what is actually there.
//!
//! **Both sources are right, and they are about different layers.** From
//! `wgpu-types` 29.0.4's own source, on `EXPERIMENTAL_MESH_SHADER`:
//!
//! > Supported platforms:
//! > - Vulkan (with `VK_EXT_mesh_shader`)
//! > - DX12
//! > - Metal
//! >
//! > **Naga is only supported on vulkan. On other platforms you will have to
//! > use passthrough shaders.**
//!
//! So the *feature* reaches Metal and the *WGSL compiler* does not. On Metal a
//! caller must hand wgpu pre-compiled MSL rather than the WGSL this crate
//! composes — which makes mesh shaders a fork in the shader pipeline, not a
//! flag on it. That is the fact GPU-008 has to be designed around, and it is a
//! documentation claim rather than a measurement, so it carries its source.
//!
//! # Enabling it requires `unsafe`, and this workspace forbids `unsafe`
//!
//! The probe **cannot** report whether a device opens with mesh shaders, and
//! that is a finding rather than a gap. `ExperimentalFeatures::enabled()` is a
//! `const unsafe fn` in `wgpu-types` 29.0.4 — its `disabled()` counterpart is
//! documented as *"uses of `Features` prefixed with EXPERIMENTAL are
//! disallowed"* — so requesting `EXPERIMENTAL_MESH_SHADER` at device creation
//! needs an `unsafe` block acknowledging possible UB.
//!
//! This workspace sets `unsafe_code = "forbid"`. So GPU-008 cannot be written
//! here at all without that policy changing, which is a decision rather than a
//! task. See M-146.
//!
//! An earlier version of this probe requested the feature with
//! `ExperimentalFeatures::disabled()` and reported `usable: false` for an
//! adapter that advertises it — measuring its own configuration and calling it
//! a hardware limit. The field is gone rather than fixed, because without
//! `unsafe` there is nothing honest for it to say.
//!
//! # What is reported
//!
//! The feature bits, for **every adapter on every backend**, not just the one
//! [`headless::Gpu`](crate::headless::Gpu) would pick — the question is what
//! this *machine* can do, and one high-performance adapter would miss a second
//! GPU or a second backend.

use crate::block_on::block_on;

/// One adapter's answer.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct MeshShaderReport {
    /// Driver's name for the device.
    pub name: String,
    /// Backend this adapter was reached through.
    pub backend: wgpu::Backend,
    /// Discrete, integrated, virtual, CPU, or unknown.
    pub device_type: wgpu::DeviceType,
    /// The adapter lists `EXPERIMENTAL_MESH_SHADER`.
    pub advertised: bool,
    /// `EXPERIMENTAL_MESH_SHADER_MULTIVIEW`.
    pub multiview: bool,
    /// `EXPERIMENTAL_MESH_SHADER_POINTS`.
    pub points: bool,
}

impl MeshShaderReport {
    /// One line, for a table.
    #[must_use]
    pub fn summary(&self) -> String {
        format!(
            "{:<40} {:<8} {:<14} advertised {:<5} multiview {:<5} points {}",
            self.name,
            format!("{:?}", self.backend),
            format!("{:?}", self.device_type),
            self.advertised,
            self.multiview,
            self.points,
        )
    }
}

/// Probe every adapter on every backend this build supports.
///
/// Returns one report per adapter, in enumeration order. An empty result means
/// no adapter at all, which is itself an answer and not an error — a machine
/// with no GPU has honestly reported that it has no mesh shaders.
#[must_use]
pub fn probe_mesh_shaders() -> Vec<MeshShaderReport> {
    let instance = wgpu::Instance::new(wgpu::InstanceDescriptor::new_without_display_handle());
    // `enumerate_adapters` is a future in wgpu 29, like the rest of adapter
    // acquisition. Same thread park as everywhere else in this crate.
    block_on(instance.enumerate_adapters(wgpu::Backends::all()))
        .into_iter()
        .map(|adapter| {
            let info = adapter.get_info();
            let features = adapter.features();
            MeshShaderReport {
                name: info.name,
                backend: info.backend,
                device_type: info.device_type,
                advertised: features.contains(wgpu::Features::EXPERIMENTAL_MESH_SHADER),
                multiview: features.contains(wgpu::Features::EXPERIMENTAL_MESH_SHADER_MULTIVIEW),
                points: features.contains(wgpu::Features::EXPERIMENTAL_MESH_SHADER_POINTS),
            }
        })
        .collect()
}

#[cfg(test)]
mod tests;
