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
    // x: world units per texture tile. y: blend sharpness. z: the forced array
    // layer -- negative blends the terrain layers by slope and depth, and `>= 0`
    // samples that one layer and nothing else, which is what the sandbox walls
    // are. w is unused.
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
@group(#{MATERIAL_BIND_GROUP}) @binding(101) var albedo_roughness_texture: texture_2d_array<f32>;
@group(#{MATERIAL_BIND_GROUP}) @binding(102) var albedo_roughness_sampler: sampler;
@group(#{MATERIAL_BIND_GROUP}) @binding(103) var normal_ao_texture: texture_2d_array<f32>;
@group(#{MATERIAL_BIND_GROUP}) @binding(104) var normal_ao_sampler: sampler;

// The interpolation. `abs(n)` because a plane does not care which way through it
// the surface faces; the power sharpens the transition so a 45-degree face is a
// blend rather than a smear, and normalising keeps the total energy at one so a
// blended fragment is not darker than an axis-aligned one.
fn plane_weights(n: vec3<f32>, sharpness: f32) -> vec3<f32> {
    let w = pow(abs(n), vec3(sharpness));
    return w / max(w.x + w.y + w.z, 1e-5);
}

// Layer order, matching the stack in `textures/PROVENANCE.md`. These *are* the
// array slice indices, so this list and that document's table are one fact
// written twice: stack the image bottom-up instead and the walls come out
// grassy with nothing to report it.
const LAYER_GRASS: i32 = 0;
const LAYER_DIRT_SURFACE: i32 = 1;
const LAYER_DIRT_DEEP: i32 = 2;
// Grass only where the surface faces up: `n.y` is the cosine off vertical, so
// this is grass on anything shallower than ~35 degrees, none steeper than ~57.
const GRASS_SLOPE_LO: f32 = 0.55;
const GRASS_SLOPE_HI: f32 = 0.82;
// `Ground`'s top lives in y = [-0.5, 0.5], so an undisturbed surface reads ~0.98
// shallow and a tunnel two units down reads pure deep dirt. Absolute world
// height rather than a re-transcription of `Ground`'s height formula into WGSL:
// two statements of one field drift apart, and the symptom of forgetting to move
// these when the amplitude changes is grass appearing inside a tunnel.
const SHALLOW_Y_LO: f32 = -1.6;
const SHALLOW_Y_HI: f32 = -0.4;

// One layer's triplanar albedo/roughness. Called from uniform control flow only.
fn layer_ar(l: i32, uv_x: vec2<f32>, uv_y: vec2<f32>, uv_z: vec2<f32>, w: vec3<f32>) -> vec4<f32> {
    return textureSample(albedo_roughness_texture, albedo_roughness_sampler, uv_x, l) * w.x
         + textureSample(albedo_roughness_texture, albedo_roughness_sampler, uv_y, l) * w.y
         + textureSample(albedo_roughness_texture, albedo_roughness_sampler, uv_z, l) * w.z;
}

// One plane's normal/AO, blended across the three terrain layers. Per plane
// rather than per fragment, because the whiteout blend below swizzles each plane
// individually: a normal blended across planes has already lost which plane it
// came from and cannot be swizzled back.
fn layers_na(uv: vec2<f32>, lw: vec3<f32>) -> vec4<f32> {
    return textureSample(normal_ao_texture, normal_ao_sampler, uv, LAYER_GRASS) * lw.x
         + textureSample(normal_ao_texture, normal_ao_sampler, uv, LAYER_DIRT_SURFACE) * lw.y
         + textureSample(normal_ao_texture, normal_ao_sampler, uv, LAYER_DIRT_DEEP) * lw.z;
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

    // Layer weights. `lw` needs no renormalising: `up*s + (1-up)*s + (1-s) == 1`
    // identically, for every `up` and every `s`.
    var lw = vec3(0.0, 0.0, 1.0);
    if (triplanar.settings.z < 0.0) {
        let up = smoothstep(GRASS_SLOPE_LO, GRASS_SLOPE_HI, n.y);
        let shallow = smoothstep(SHALLOW_Y_LO, SHALLOW_Y_HI, in.world_position.y);
        lw = vec3(up * shallow, (1.0 - up) * shallow, 1.0 - shallow);
    }
    var ar: vec4<f32>;
    var na_x: vec4<f32>;
    var na_y: vec4<f32>;
    var na_z: vec4<f32>;
    // The branch is on a **uniform**, so both arms are uniform control flow and
    // `textureSample`'s implicit derivatives stay legal. A per-fragment branch
    // around a sample is a WGSL validation error, which is why the blend always
    // pays for all three layers rather than skipping one whose weight is small.
    if (triplanar.settings.z >= 0.0) {
        let l = i32(triplanar.settings.z);
        ar = layer_ar(l, uv_x, uv_y, uv_z, w);
        na_x = textureSample(normal_ao_texture, normal_ao_sampler, uv_x, l);
        na_y = textureSample(normal_ao_texture, normal_ao_sampler, uv_y, l);
        na_z = textureSample(normal_ao_texture, normal_ao_sampler, uv_z, l);
    } else {
        ar = layer_ar(LAYER_GRASS, uv_x, uv_y, uv_z, w) * lw.x
           + layer_ar(LAYER_DIRT_SURFACE, uv_x, uv_y, uv_z, w) * lw.y
           + layer_ar(LAYER_DIRT_DEEP, uv_x, uv_y, uv_z, w) * lw.z;
        na_x = layers_na(uv_x, lw);
        na_y = layers_na(uv_y, lw);
        na_z = layers_na(uv_z, lw);
    }

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
