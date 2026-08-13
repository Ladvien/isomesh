//! The refusal path, which is the one this crate can actually run.
//!
//! `headless::Gpu` opens a device with `Features::empty()` and cannot ask for
//! `EXPERIMENTAL_MESH_SHADER` without `unsafe` (M-146, GPU-009). That makes it
//! the perfect fixture for the *unsupported* case — the half of GPU-008's
//! acceptance that says "never panics on an unsupported adapter" — and it means
//! the success path is exercised by the Bevy example rather than here.
#![allow(clippy::float_cmp)]

use super::MeshShaderRenderer;
use crate::headless::Gpu;
use crate::{Composer, Error};

/// The acceptance criterion, on the device this crate can build.
///
/// An error value, not an abort, and not a silently different pipeline.
#[test]
fn an_unsupported_device_is_refused_rather_than_crashed() {
    let gpu = Gpu::new().expect("a GPU adapter -- no software fallback, by design");

    assert!(
        !MeshShaderRenderer::is_supported(gpu.device()),
        "the headless device reports mesh shaders -- it requests Features::empty(), so either \
         wgpu changed or this test is no longer measuring the unsupported case"
    );

    let built = MeshShaderRenderer::new(gpu.device(), wgpu::TextureFormat::Rgba8UnormSrgb, None);
    assert_eq!(
        built.err(),
        Some(Error::MeshShadersUnavailable),
        "an unsupported device must be refused with an error"
    );
}

/// The shader is registered, so GPU-003's no-GPU sweep validates it.
///
/// Without this the mesh shader would be checked only where a capable device
/// exists, which is exactly the coverage GPU-003 was built to avoid depending
/// on.
#[test]
fn the_mesh_shader_is_a_builtin_and_composes() {
    let composer = Composer::with_builtins();
    assert!(
        composer.module_names().contains(&"mesh_render"),
        "mesh_render is not registered, so the validation sweep does not cover it"
    );

    let source = composer.compose("mesh_render", &[]).expect("composes");
    for expected in [
        "enable wgpu_mesh_shader;",
        "@mesh(mesh_out)",
        "@builtin(vertex_count)",
        "@builtin(primitives)",
        "fn draw_mesh",
        "fn shade",
    ] {
        assert!(
            source.contains(expected),
            "mesh_render.wgsl lost `{expected}`"
        );
    }
}

/// The uniform packings, which have no device to check them.
#[test]
fn the_camera_uniform_is_column_major_and_64_bytes() {
    // A matrix whose entries name their own position, so a transposition is
    // visible rather than plausible.
    let mut m = [[0.0f32; 4]; 4];
    for (c, column) in m.iter_mut().enumerate() {
        for (r, value) in column.iter_mut().enumerate() {
            *value = (c * 4 + r) as f32;
        }
    }
    let bytes = MeshShaderRenderer::camera_bytes(m);
    assert_eq!(bytes.len(), 64);

    let word = |i: usize| {
        f32::from_le_bytes([
            bytes[i * 4],
            bytes[i * 4 + 1],
            bytes[i * 4 + 2],
            bytes[i * 4 + 3],
        ])
    };
    // Column-major: element [c][r] lands at word c*4 + r, which is what WGSL
    // reads a `mat4x4<f32>` as.
    for c in 0..4 {
        for r in 0..4 {
            assert_eq!(
                word(c * 4 + r),
                (c * 4 + r) as f32,
                "element [{c}][{r}] moved"
            );
        }
    }
}

#[test]
fn the_draw_uniform_is_padded_to_sixteen_bytes() {
    let bytes = MeshShaderRenderer::draw_bytes(3704);
    assert_eq!(bytes.len(), 16, "a uniform binding is padded to 16 bytes");
    assert_eq!(
        u32::from_le_bytes([bytes[0], bytes[1], bytes[2], bytes[3]]),
        3704
    );
    assert!(
        bytes[4..].iter().all(|b| *b == 0),
        "padding must be zero, not whatever was on the stack"
    );
}
