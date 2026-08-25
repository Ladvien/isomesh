// Interpolated triplanar mapping for isosurface terrain.
//
// An isosurface has no natural parameterisation, so there is no correct UV to
// put on a vertex. `bevy_isomesh`'s `MeshBuilder` emits a dominant-axis planar
// projection as a stand-in and its own doc comment names the cost: it seams
// visibly wherever the dominant axis flips. This samples all three planes and
// interpolates between them by the normal instead, so there is no seam to see.
//
// The projection is a function of world position alone, never of the mesh, which
// is what makes it continuous across a chunk boundary -- two chunks meshed
// independently agree on the texture as exactly as they agree on the surface.

#import bevy_pbr::{
    pbr_fragment::pbr_input_from_standard_material,
    pbr_functions::alpha_discard,
}

#ifdef PREPASS_PIPELINE
#import bevy_pbr::{
    prepass_io::{VertexOutput, FragmentOutput},
    pbr_deferred_functions::deferred_output,
}
#else
#import bevy_pbr::{
    forward_io::{VertexOutput, FragmentOutput},
    pbr_functions::{apply_pbr_lighting, main_pass_post_lighting_processing},
}
#endif

struct Triplanar {
    // x: world units per texture tile. y: blend sharpness. z and w are unused.
    //
    // One `vec4` rather than two `f32`s, and that is not laziness: WebGL2
    // requires a uniform struct to be 16-byte aligned, and a `vec4` is 16 bytes
    // by construction. Bevy's own `extended_material` example pads with three
    // `#[cfg(feature = "webgl2")]` fields, which cannot work outside the `bevy`
    // package because `webgl2` is not a feature of this one -- the cfg would be
    // permanently false and the uniform would be 4 bytes where this expects 16.
    settings: vec4<f32>,
}

// `#{MATERIAL_BIND_GROUP}` is a shader-def, not a literal: `bevy_pbr` pushes it
// from `MATERIAL_BIND_GROUP_INDEX`, which is 3 in 0.19. Bindings start at 100
// because the base `StandardMaterial` owns 0..=12 with Bevy's default features
// (0..=20 with every texture feature on, and 0..31 reserved by the bindless
// index table), and `ExtendedMaterial`'s `AsBindGroup` silently *drops* an
// extension entry whose number collides -- the base wins the filter, and the
// symptom is a pipeline that compiles and samples nothing.
@group(#{MATERIAL_BIND_GROUP}) @binding(100) var<uniform> triplanar: Triplanar;
@group(#{MATERIAL_BIND_GROUP}) @binding(101) var albedo_roughness_texture: texture_2d<f32>;
@group(#{MATERIAL_BIND_GROUP}) @binding(102) var albedo_roughness_sampler: sampler;
@group(#{MATERIAL_BIND_GROUP}) @binding(103) var normal_ao_texture: texture_2d<f32>;
@group(#{MATERIAL_BIND_GROUP}) @binding(104) var normal_ao_sampler: sampler;

// The interpolation. `abs(n)` because a plane does not care which way through it
// the surface faces; the power sharpens the transition so a 45-degree face is a
// blend rather than a smear, and normalising keeps the total energy at one so a
// blended fragment is not darker than an axis-aligned one.
fn plane_weights(n: vec3<f32>, sharpness: f32) -> vec3<f32> {
    let w = pow(abs(n), vec3(sharpness));
    return w / max(w.x + w.y + w.z, 1e-5);
}

@fragment
fn fragment(
    in: VertexOutput,
    @builtin(front_facing) is_front: bool,
) -> FragmentOutput {
    var pbr_input = pbr_input_from_standard_material(in, is_front);

    let p = in.world_position.xyz / triplanar.settings.x;
    let n = normalize(pbr_input.world_normal);
    let w = plane_weights(n, triplanar.settings.y);

    // One coordinate pair per plane, each the two world axes that plane spans.
    let uv_x = p.zy;
    let uv_y = p.xz;
    let uv_z = p.xy;

    let ar_x = textureSample(albedo_roughness_texture, albedo_roughness_sampler, uv_x);
    let ar_y = textureSample(albedo_roughness_texture, albedo_roughness_sampler, uv_y);
    let ar_z = textureSample(albedo_roughness_texture, albedo_roughness_sampler, uv_z);
    let ar = ar_x * w.x + ar_y * w.y + ar_z * w.z;

    let na_x = textureSample(normal_ao_texture, normal_ao_sampler, uv_x);
    let na_y = textureSample(normal_ao_texture, normal_ao_sampler, uv_y);
    let na_z = textureSample(normal_ao_texture, normal_ao_sampler, uv_z);

    pbr_input.material.base_color = vec4(
        pbr_input.material.base_color.rgb * ar.rgb,
        pbr_input.material.base_color.a,
    );
    // Clamped off zero: a perfectly smooth fragment is a mirror, and a mirror of
    // a one-light scene is a black hole in the rock.
    pbr_input.material.perceptual_roughness = clamp(ar.a, 0.05, 1.0);
    pbr_input.diffuse_occlusion = vec3(na_x.a * w.x + na_y.a * w.y + na_z.a * w.z);

    // Whiteout blend, and Bevy's own normal mapping is unreachable rather than
    // merely unused: `apply_normal_mapping` needs a TBN from
    // `calculate_tbn_mikktspace`, which needs `ATTRIBUTE_TANGENT`, which the
    // mesher does not emit -- so `pbr_fragment.wgsl`'s entire
    // `#ifdef VERTEX_TANGENTS` block is compiled out and `pbr_input.N` has to be
    // written here. Each plane's tangent normal is added to the geometric
    // normal's other two components and swizzled back into world space; `abs` on
    // the third component keeps the sign coming from the geometry rather than
    // from the texture.
    var t_x = na_x.xyz * 2.0 - 1.0;
    var t_y = na_y.xyz * 2.0 - 1.0;
    var t_z = na_z.xyz * 2.0 - 1.0;
    t_x = vec3(t_x.xy + n.zy, abs(t_x.z) * n.x);
    t_y = vec3(t_y.xy + n.xz, abs(t_y.z) * n.y);
    t_z = vec3(t_z.xy + n.xy, abs(t_z.z) * n.z);
    pbr_input.N = normalize(t_x.zyx * w.x + t_y.xzy * w.y + t_z.xyz * w.z);
    pbr_input.clearcoat_N = pbr_input.N;

    pbr_input.material.base_color = alpha_discard(pbr_input.material, pbr_input.material.base_color);

#ifdef PREPASS_PIPELINE
    let out = deferred_output(in, pbr_input);
#else
    var out: FragmentOutput;
    out.color = apply_pbr_lighting(pbr_input);
    out.color = main_pass_post_lighting_processing(pbr_input, out.color);
#endif

    return out;
}
